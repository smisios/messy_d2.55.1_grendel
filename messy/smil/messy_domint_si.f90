!> mz_sg_20200724 | f90

!! DOMINT: a sub-module for calculating domain integrals
!! SI level implementation for ECHAM5
!!
!! Author: S.Gromov (MPIC), 2018-2020
!!
!! CHANGELOG:
!! v1.0 (20200724): initial implementation, atmospheric vertical/horizontal domains
!! v1.1 (20201102): added surface domains, added wildcards support for channel names

#if defined(__GFORTRAN__) || defined(__G95__)
#define GF_VERSION (__GNUC__)
#if (GF_VERSION > 6)
#define string character(:),allocatable
#else
#define string character(LEN=250)
#endif
#else
#define string character(:),allocatable
#endif

#define _CCHSS_ call channel_halt(substr,status)
!! #define DEBUG /* uncomment for debug */

module messy_domint_si

#ifdef ECHAM5

  use messy_main_constants_mem, only: STRLEN_XLONG
  use messy_main_timer_event,   only: time_event, io_time_event
  use messy_main_blather_bi,    only: start_message_bi, end_message_bi, &
                                      info_bi, warning_bi, error_bi
  use messy_main_channel,       only: STRLEN_CHANNEL, STRLEN_OBJECT
  use messy_domint

  implicit none
  private
  save

! ----- types and data structures -----

! --- namelist entries

! list of tasks to integrate
  type t_domint_taskinfo
    integer                     :: comp = 0     ! domain composition no. (e.g. 1 = atmosphere zonal, 2 = ...)
    character(len=STRLEN_CHANNEL+1+STRLEN_OBJECT) :: &
      var = '', &  ! 'ch:ob' reference to the variable/weighting fields
      wgt = '', &
      info = ''    ! optional info string containing name:caption:unit
    real(dp)                    :: sf = 1._dp   ! optional scaling factor
  end type t_domint_taskinfo
  public :: t_domint_taskinfo
  integer, parameter            :: domint_ntask_max = 500
  type(t_domint_taskinfo), &
    dimension(domint_ntask_max) :: domint       ! task namelist entries

! input field references set (add more custom req. input here)
  type t_domint_inputset
    type(t_domint_input) ::  &
    sinlat, sinlon, press,   &  ! for vert./horiz. location
    theta, tph_idx, blh_idx, &  ! for LS, TP & BL detection
    landfrac,                &  ! for land/ocean detection
    ozone                       ! for O3-based tropopause
  ! .....
  end type t_domint_inputset

! domint domain/tasks integration struct
  type t_domint
  ! calculation trigger frequency, default: every step
    character(len=STRLEN_CHANNEL) &      ! customisable name
                            :: name = ''
    logical                 :: ldisable = .false.
    type(io_time_event)     :: trigcalc = io_time_event(1,'steps','first',0)
    character(len=STRLEN_XLONG) &
                            :: input     ! input fields ref. str ( input1=channel:object; input2=channel:object; ... )
  end type t_domint
  integer, parameter        :: domint_ndom = 2  ! no of implemented domain kinds
  type(t_domint)            :: domain(1:domint_ndom)

! work arrays for calculation @CL
  type t_domint_work
    logical                 :: ltrigcalc = .false., lupdated = .false.
    type(time_event)        :: event
    integer                 :: reprid
    type(t_domint_comp)     :: comp      ! composition
    type(t_domint_inputset) :: input
    type(t_domint_task), pointer &
                            :: tasks(:)  ! integration tasks
  end type t_domint_work
  type(t_domint_work)       :: work(1:domint_ndom)


! --- etc.
  real(dp)                  :: thres_land  = 0.90_dp, &  ! land/coast definition threshold - yield ~23%(land)+~6%(coast)=
                               thres_coast = 0.45_dp     !   ~29% [at T42] (current land/ocean estimate is 29.2%/70.8%)
  real(dp), parameter       :: FLAG_UNDEF = -999._dp     ! undefined value flag

! ----- domain compositions -----

! #1: atmospheric zonal-vertical domains composition
  integer, parameter     :: domint_atm_VZ = 1                ! composition #
  integer                :: nvdom, nzdom                     ! working arrays for domain sizes
  real(dp), allocatable  :: vdom_mids(:), vdom_bnds(:), &    ! domain coords/bounds
                            zdom_mids(:), zdom_bnds(:)
  character(len=*), parameter :: &                           ! domain information captions
    vdom_info = '(k=1-11) 1:ATM 2:T100 3:TROP 4:TO3 5:MS STRAT=(6:US 7:LS) 8:TP TROP=(9:FT 10:BL) 11:SRF;' &
              //' T100 | TO3: troposphere at P>100hPa | O3<150 ppb', &
    zdom_info = '(j=1-11) 1:GLOB 2:SH 3:ANT (4:ETSH 5:IT 6:ETNH) 7:ARC 8:NH (9:LAND 10:OCEAN 11:COAST);' &
              //' GLOB=(NH+SH)=(ETNH+IT+ETSH)=(LAND+OCEAN+COAST)'

! #2: surface domins composition
  integer, parameter     :: domint_surf = 2                  ! composition #
  integer                :: nsdom                            ! working arrays for domain sizes
  real(dp), allocatable  :: sdom_mids(:), sdom_bnds(:)       ! domain coords/bounds
  character(len=*), parameter :: &                           ! domain information captions
    sdom_info = '(j=1-11) 1:GLOB 2:SH 3:ANT (4:ETSH 5:IT 6:ETNH) 7:ARC 8:NH (9:LAND 10:OCEAN 11:COAST);' &
              //' GLOB=(NH+SH)=(ETNH+IT+ETSH)=(LAND+OCEAN+COAST)'

! ----- interface -----
  public :: domint_initialize
  public :: domint_init_memory
  public :: domint_init_coupling
  public :: domint_global_start
  public :: domint_global_end
  public :: domint_free_memory

contains

! ----------------------------------------------------------------------

  subroutine domint_initialize

    use messy_main_mpi_bi,     only: p_parallel_io, p_io, p_bcast
    use messy_main_timer_bi,   only: p_bcast_event
    use messy_main_tools,      only: find_next_free_unit

    implicit none

    character(len=*), parameter :: substr = modstr//'_initialize'
    integer                     :: iou    ! I/O unit
    integer                     :: status ! error status
    integer                     :: jt, jd

    call start_message_bi(modstr, 'initialize', substr)

  ! read namelist
    if (p_parallel_io) then
      iou = find_next_free_unit(100,200)
      call domint_read_nml_cpl(status, iou)
      if (status /= 0) call error_bi(' ',substr)
    endif

  ! broadcast tasks
    do jt = 1, domint_ntask_max
      call p_bcast(domint(jt)%comp, p_io)
      call p_bcast(domint(jt)%var,  p_io)
      call p_bcast(domint(jt)%wgt,  p_io)
      call p_bcast(domint(jt)%info, p_io)
      call p_bcast(domint(jt)%sf,   p_io)
    enddo

  ! process domain info/input
    do jd = 1, domint_ndom
    ! parse input string for parameters
      if (p_parallel_io) then
#ifdef DEBUG
        print *,'domint_initialize(',jd,') input = ',trim(domain(jd)%input)
#endif
        work(jd)%input%press%ref    = get_par(trim(domain(jd)%input),'press')
        work(jd)%input%sinlat%ref   = get_par(trim(domain(jd)%input),'sinlat')
        work(jd)%input%sinlon%ref   = get_par(trim(domain(jd)%input),'sinlon')
        work(jd)%input%theta%ref    = get_par(trim(domain(jd)%input),'theta')
        work(jd)%input%tph_idx%ref  = get_par(trim(domain(jd)%input),'tph_idx')
        work(jd)%input%blh_idx%ref  = get_par(trim(domain(jd)%input),'blh_idx')
        work(jd)%input%landfrac%ref = get_par(trim(domain(jd)%input),'landfrac')
        work(jd)%input%ozone%ref    = get_par(trim(domain(jd)%input),'ozone')
      endif
    ! broadcast
      call p_bcast(domain(jd)%ldisable, p_io)
      call p_bcast(domain(jd)%name, p_io)
      call p_bcast_event(domain(jd)%trigcalc, p_io)
      call p_bcast(work(jd)%input%press%ref, p_io)
      call p_bcast(work(jd)%input%sinlat%ref, p_io)
      call p_bcast(work(jd)%input%sinlon%ref, p_io)
      call p_bcast(work(jd)%input%theta%ref, p_io)
      call p_bcast(work(jd)%input%tph_idx%ref, p_io)
      call p_bcast(work(jd)%input%blh_idx%ref, p_io)
      call p_bcast(work(jd)%input%landfrac%ref, p_io)
      call p_bcast(work(jd)%input%ozone%ref, p_io)
    enddo

  ! done
    call end_message_bi(modstr, 'initialize', substr)

  contains

    function get_par(inp, par)
      string                       :: get_par
      character(len=*), intent(in) :: inp, par    ! parameter to search
      string  :: p, q
      integer :: jp
      get_par = ''
      jp = 1; p = '*'
      do while ( p.ne.'' )
        p = str_field(inp,';',jp); q = str_field(p,'=',1)
        if ( q.eq.par ) then
          get_par = str_field(p,'=',2)
          exit
        endif
        jp = jp + 1
      enddo
    end function get_par

  end subroutine domint_initialize

! ----------------------------------------------------------------------

  subroutine domint_init_memory

    use messy_main_grid_def_mem_bi,  only: nproma, nlev, ngpblks
    use messy_main_timer_bi,         only: timer_event_init
    use messy_main_channel_error_bi, only: channel_halt
    use messy_main_channel,          only: new_channel, new_channel_object, new_attribute

    implicit none

    intrinsic :: trim

    character(len=*), parameter   :: substr = modstr//'_init_memory'
    integer                       :: status, jd
    character(len=STRLEN_CHANNEL) :: chn
    character(len=STRLEN_OBJECT)  :: obn
    real(dp), pointer             :: mem(:,:,:,:)

  ! ----- set up domains

    call start_message_bi(modstr, 'initialize domains', substr)

  ! process domain info/input
    do jd = 1, domint_ndom

      if ( domain(jd)%ldisable ) cycle

      select case(jd)

      case(domint_atm_VZ) ! atmosphere, vertical-zonal
      ! assign name (if not given already)
        if ( trim(domain(1)%name).eq.'' ) domain(jd)%name = 'atm-VZ'
      ! define VZ domains and new output representation
        call domint_create_domain_repr_atm_VZ(status, work(jd)%reprid)
      ! allocate memory at core level
        status = domint_allocate_comp(work(jd)%comp, (/nvdom,nzdom/), (/nproma,nlev,ngpblks/))
        if ( status .ne. 0 ) &
          call error_bi('problem during allocation of atm_VZ domain comp', substr)
      ! create the output channel
        chn = trim(modstr)//'_'//trim(domain(jd)%name)
        call new_channel(status, trim(chn), reprid=work(jd)%reprid); _CCHSS_
      ! domains info (as global attribute)
        call new_attribute(status, trim(chn), 'domains_vertical', c=vdom_info); _CCHSS_
        call new_attribute(status, trim(chn), 'domains_horizontal', c=zdom_info); _CCHSS_

      case(domint_surf) ! surface
      ! assign name (if not given already)
        if ( trim(domain(jd)%name).eq.'' ) domain(jd)%name = 'surf'
      ! define surfacedomains and new output representation
        call domint_create_domain_repr_surf(status, work(jd)%reprid)
      ! allocate memory at core level
        status = domint_allocate_comp(work(jd)%comp, (/nsdom/), (/nproma,ngpblks/))
        if ( status .ne. 0 ) &
          call error_bi('problem during allocation of surf domain comp', substr)
      ! create the output channel
        chn = trim(modstr)//'_'//trim(domain(jd)%name)
        call new_channel(status, trim(chn), reprid=work(jd)%reprid); _CCHSS_
      ! domains info (as global attribute)
        call new_attribute(status, trim(chn), 'domains_surface', c=sdom_info); _CCHSS_

      case default
        call error_bi('unknown domain kind encountered',substr)

      end select

    ! no. of cells
      obn = 'cells'
      mem(1:work(jd)%comp%ddims(1), &
          1:work(jd)%comp%ddims(2), &
          1:work(jd)%comp%ddims(3), 1:1) => work(jd)%comp%cells(1:work(jd)%comp%ndom)
      call new_channel_object(status, trim(chn), trim(obn), mem=mem); _CCHSS_
      call new_attribute(status, trim(chn), trim(obn), 'long_name',     c='domain gridcell count' ); _CCHSS_
      call new_attribute(status, trim(chn), trim(obn), 'units',         c='count'); _CCHSS_
      call new_attribute(status, trim(chn), trim(obn), 'missing_value', r=FLAG_UNDEF); _CCHSS_
      call new_attribute(status, trim(chn), trim(obn), '_FillValue',    r=FLAG_UNDEF); _CCHSS_

    ! set calculation trigger
      call timer_event_init(work(jd)%event, domain(jd)%trigcalc, &
                             modstr//'_trigcalc_'//trim(domain(jd)%name), 'next')
    enddo

  ! finally
    call end_message_bi(modstr, 'initialize domains', substr)

  end subroutine domint_init_memory

! ----------------------------------------------------------------------

  subroutine domint_init_coupling

    use messy_main_mpi_bi,           only: p_parallel_io
    USE messy_main_tools,            only: match_wild
    use messy_main_channel,          only: NCHANNEL, get_channel_name, get_channel_info, &
                                           new_channel_object, get_attribute, new_attribute
    use messy_main_channel_error_bi, only: channel_halt

    implicit none

    character(len=*), parameter :: substr = modstr//'_init_cpl'
    integer                     :: jd, jt, jc, jo, jl, status, di
    string                      :: chn, obn, s, u
    character(len=STRLEN_OBJECT), &
                        pointer :: obns(:) => null()
    character(len=STRLEN_CHANNEL) &
                                :: chns
    real(dp)                    :: r

    call start_message_bi(modstr, 'initialize coupling', substr)

    domain_loop: do jd = 1, domint_ndom

    ! skip diabled domains
      if ( domain(jd)%ldisable ) cycle

    ! get pointers to the [required] fields

      select case(jd)

      case(domint_atm_VZ) ! atmosphere vertica/horizontal
      ! required
        call assign_input('pressure',       work(jd)%input%press,  .true., p3=.true.)
        call assign_input('sin(latitude)',  work(jd)%input%sinlat, .true., p2=.true.)
      !!call assign_input('sin(longitude)', work(jd)%input%sinlon, .true., p2=.true.)
      ! optional
        call assign_input('pot. temperature',       work(jd)%input%theta,    .false., p3=.true.)
        call assign_input('tropop. layer idx.',     work(jd)%input%tph_idx,  .false., p2=.true.)
        call assign_input('bound. hgt. layer idx.', work(jd)%input%blh_idx,  .false., p2=.true.)
        call assign_input('land fraction',          work(jd)%input%landfrac, .false., p2=.true.)
        call assign_input('ozone mixing ratio',     work(jd)%input%ozone,    .false., p3=.true.)

      case(domint_surf) ! surface
      ! required
        call assign_input('sin(latitude)',  work(jd)%input%sinlat, .true., p2=.true.)
      !!call assign_input('sin(longitude)', work(jd)%input%sinlon, .true., p2=.true.)
      ! optional
        call assign_input('land fraction',  work(jd)%input%landfrac, .false., p2=.true.)

      case default
        call error_bi('unknown domain kind encountered',substr)

      end select

    enddo domain_loop

  ! scan task request list for entries, add using wildcard accoriding to composition
    tasklist_loop: do jl = 1, domint_ntask_max

    ! skip empty/not implemented domain entries
      if ( ( domint(jl)%comp.lt.1 ).or.( domint(jl)%comp.gt.domint_ndom ) ) cycle
      jd = domint(jl)%comp

    ! split channel:object names
      chn = str_field(trim(domint(jl)%var),':',1)
      obn = str_field(trim(domint(jl)%var),':',2)

    ! skip empty records
      if ( (chn.eq.'') .or. (obn.eq.'') ) cycle

#ifdef DEBUG
      print *, substr//': task #',jd,' chn =',chn,' obn =',obn,' size(obns) =',size(obns)
#endif

    ! scan for matching channels/objects
      channel_loop: do jc = 1, NCHANNEL

      ! reset channel name (in case it hase wildcard)
        chn = str_field(trim(domint(jl)%var),':',1)

        chns = ''
        call get_channel_name(status,jc,chns); if ( status /= 0 ) cycle
        if (.not.match_wild(chn,trim(adjustl(chns)))) cycle

        chn = trim(adjustl(chns))

      ! scan for channel objects
        call get_channel_info(status, chn, ONAMES=obns)
        if (status /= 3003) then ! channel (name) does not exist
          call channel_halt(substr, status)
        else
          call warning_bi(' ... channel '//chn//' does not exist, although matches wildcard. (?)', substr)
          cycle
        endif

        object_loop: do jo = 1, size(obns)
          if ( match_wild(obn, trim(obns(jo))) ) then

          ! it's a match, add a task
            jt = domint_allocate_task(work(jd)%comp,work(jd)%tasks)

          ! data refs.
            work(jd)%tasks(jt)%var%ref = trim(chn//':'//trim(obns(jo)))
            work(jd)%tasks(jt)%wgt%ref = trim(domint(jl)%wgt)

          ! info
            work(jd)%tasks(jt)%name = chn//'_'//trim(obns(jo))
            work(jd)%tasks(jt)%capt = ''
            work(jd)%tasks(jt)%unit = ''

          ! assign input (domain-dependent)
            select case(jd)
            case(domint_atm_VZ)
              call assign_input('int:'//work(jd)%tasks(jt)%name, work(jd)%tasks(jt)%var, .true., p3=.true.)
            case(domint_surf)
              call assign_input('int:'//work(jd)%tasks(jt)%name, work(jd)%tasks(jt)%var, .true., p2=.true.)
            case default
              call error_bi('unknown domain kind encountered',substr)
            end select
            work(jd)%tasks(jt)%info = trim(work(jd)%tasks(jt)%var%ref)

          ! if variable info is present, use it
            if ( work(jd)%tasks(jt)%var%capt/='' ) then
              work(jd)%tasks(jt)%capt = work(jd)%tasks(jt)%var%capt
            else
              work(jd)%tasks(jt)%capt = trim(obns(jo))
            endif
            work(jd)%tasks(jt)%unit = work(jd)%tasks(jt)%var%unit

          ! weighting
            if ( trim(work(jd)%tasks(jt)%wgt%ref)/='' ) then
            ! assign input (domain-dependent)
              select case(jd)
              case(domint_atm_VZ)
                call assign_input('\--- wgt:'//work(jd)%tasks(jt)%name, work(jd)%tasks(jt)%wgt, .true., p3=.true.)
              case(domint_surf)
                call assign_input('\--- wgt:'//work(jd)%tasks(jt)%name, work(jd)%tasks(jt)%wgt, .true., p2=.true.)
              case default
                call error_bi('unknown domain kind encountered',substr)
              end select
              work(jd)%tasks(jt)%lwgt = .true.
              work(jd)%tasks(jt)%info = work(jd)%tasks(jt)%info//' weighted by '//trim(work(jd)%tasks(jt)%wgt%ref)
            ! update units
              if ( work(jd)%tasks(jt)%wgt%unit/='' ) &
                work(jd)%tasks(jt)%unit = work(jd)%tasks(jt)%unit//'*'//work(jd)%tasks(jt)%wgt%unit
            endif

          ! if custom info is present, use it
            s = str_field(domint(jl)%info,':',1); if ( s/='' ) work(jd)%tasks(jt)%name = s
            s = str_field(domint(jl)%info,':',2); if ( s/='' ) work(jd)%tasks(jt)%capt = s
            s = str_field(domint(jl)%info,':',3); if ( s/='' ) work(jd)%tasks(jt)%unit = s

          ! scaling factor
            work(jd)%tasks(jt)%sf = domint(jl)%sf
            if ( work(jd)%tasks(jt)%sf .ne. 1._dp ) &
              work(jd)%tasks(jt)%info = work(jd)%tasks(jt)%info//' scaled by '//trim(num2str(r=work(jd)%tasks(jt)%sf,fmt='(E15.7)'))

          ! create output object
            s = trim(modstr)//'_'//trim(domain(jd)%name)
            u = work(jd)%tasks(jt)%name
            call new_channel_object(status, s, u,           mem=work(jd)%tasks(jt)%int_mem); _CCHSS_
            call new_attribute(status, s, u, 'long_name',     c=work(jd)%tasks(jt)%capt); _CCHSS_
            call new_attribute(status, s, u, 'units',         c=work(jd)%tasks(jt)%unit); _CCHSS_
            call new_attribute(status, s, u, 'info',          c=work(jd)%tasks(jt)%info); _CCHSS_
            call new_attribute(status, s, u, 'missing_value', r=FLAG_UNDEF); _CCHSS_
            call new_attribute(status, s, u, '_FillValue',    r=FLAG_UNDEF); _CCHSS_
          ! add .molarmass attribute, if present in var
            call get_attribute(status, chn, trim(obns(jo)), 'molarmass', r=r)
            if ( status.eq.0 ) then
              call new_attribute(status, s, u, 'molarmass', r=r); _CCHSS_
            endif
          endif ! if match_wild()

        enddo object_loop

      ! cleanup
        if (associated(obns)) then
          deallocate(obns)
          nullify(obns)
        endif

      enddo channel_loop

    enddo tasklist_loop

  ! done
    call end_message_bi(modstr, 'initialize coupling', substr)

  contains

    subroutine assign_input(what, inp, musthave, p2, p3)

      use messy_main_channel,       only: get_channel_object, get_attribute

      character(len=*), intent(in)        :: what         ! parameter name (for info)
      type(t_domint_input), intent(inout) :: inp          ! 'channel:object' reference, flag for found
      logical, intent(in)                 :: musthave
      logical, intent(in), optional       :: p2, p3

      integer                             :: status
      string                              :: msg, chn, obn
      character(len=STRLEN_XLONG)         :: s

      msg = ' ... '//modstr//'::assign_input - '//trim(what)//' from <'//trim(inp%ref)//'> ... '

    ! split channel:object names
      chn = str_field(trim(inp%ref),':',1)
      obn = str_field(trim(inp%ref),':',2)

      if ( (chn/='').and.(obn/='') ) then
        if ( present(p2).eqv.present(p3) ) then
          msg = msg//'error: need either p2 or p3 rank field pointer as the parameter.'
          inp%is = .false.
        else
          if (present(p2)) call get_channel_object(status, chn, obn, p2=inp%p2)
          if (present(p3)) call get_channel_object(status, chn, obn, p3=inp%p3)
          if ( status/=0 ) then
            msg = msg//'not found - variables/domains identified using this variable will not be calculated'
            call info_bi(msg, substr)
            inp%is = .false.
          else
            inp%is = .true.
          ! get caption & unit, if present
            inp%unit = ''; inp%capt = ''
            call get_attribute(status, chn, obn, 'units', c=s)
            if ( status.eq.0 ) inp%unit = trim(s)
            call get_attribute(status, chn, obn, 'long_name', c=s)
            if ( status.eq.0 ) inp%capt = trim(s)
            msg = msg//'('//inp%capt//') ['//inp%unit//'] - ok'
          endif
        endif
      else   ! problem with input ref.
        inp%is = .false.
        msg = msg//'channel and/or object names are not specified/recognised (should be "channel:object").'
      endif

#ifndef DEBUG
        if ( p_parallel_io ) &
#endif
          write(*,*) msg

      if ( (.not.inp%is).and.musthave ) &
        call error_bi('could not initialise required input ('//trim(what)//'), stop', substr)

    end subroutine assign_input

  end subroutine domint_init_coupling

! ----------------------------------------------------------------------

  subroutine domint_create_domain_repr_atm_VZ(status, reprid_domain)

    use messy_main_channel_bi,         only: DC_BC!, SCALAR
    use messy_main_channel_error_bi,   only: channel_halt
    use messy_main_channel_dimensions, only: get_dimension_info, new_dimension, &
                                             add_dimension_variable, add_dimension_variable_att
    use messy_main_channel_repr,       only: get_representation_info,new_representation, &
                                             set_representation_decomp, &
                                             IRANK, PIOTYPE_COL!, AUTO

    implicit none

    integer, intent(out)        :: status, reprid_domain

    ! LOCAL
    character(len=*), parameter :: substr = modstr//'_create_domain_repr_atm_VZ'
    character(len=*), parameter :: repr_name = 'GP_2D_DOMAIN_ATM_VZ', &
                                   dim_vdom_name = 'domain_vertical', &
                                   dim_zdom_name = 'domain_horizontal'
    integer                     :: dimid_vdom, dimid_zdom!, dimid_bnds, &
                                 ! reprid_vdom_bnd, reprid_zdom_bnd
   !integer                     :: dim_vdom_len, dim_zdom_len
    integer                     :: status_dv, status_dz

    ! PARALLEL DECOMPOSITION
    integer                          :: nseg = 0
    integer, dimension(:,:), pointer :: start => NULL()
    integer, dimension(:,:), pointer :: cnt   => NULL()
    integer, dimension(:,:), pointer :: meml  => NULL()
    integer, dimension(:,:), pointer :: memu  => NULL()

  ! adopted an example of messy_main_rnd_bi

  ! first check whether dimensions/representation exist
    call get_representation_info(status, repr_name, reprid_domain)
    call get_dimension_info(status_dv, dim_vdom_name, dimid_vdom, nvdom)
    call get_dimension_info(status_dz, dim_zdom_name, dimid_zdom, nzdom)
    if (.not.((status==2003).and.(status_dv==905).and.(status_dz==905))) then
      print *,trim(substr)//': error checking for dimension/representation for '//repr_name
      print *,'  status = ',status,'  status_dv = ',status_dv,'  status_dz = ',status_dz
      status = status+status_dv+status_dz
      return
    endif

  ! vertical domains (applicable to the global 2D+ model only) - aligned to log(P) levels
    nvdom = 11 ! ATM, T100, T, TO3, MS, US, LS, TP, FT, BL, SRF
    allocate(vdom_mids(1:nvdom))
    allocate(vdom_bnds(1:nvdom+1))
  !                      ATM  T100 T    TO3 MS   US     LS     TP     FT     BL    SRF
    vdom_bnds(:) = (/ 0.1, 0.2, 0.3, 0.5, 1.0, 1e2, 100e2, 200e2, 240e2, 750e2, 1e5, 1.25e5 /)
    vdom_bnds(:) = log10(vdom_bnds(:))
    vdom_mids(:) = (vdom_bnds(1:nvdom)+vdom_bnds(2:nvdom+1))/2.
    vdom_mids(1) = -1.        ! place ATM at 10^-1
    vdom_mids(nvdom) = 5      !       SRF at 10^5
  ! vdom_ATM=1 _T100=2 _T=3 _TO3=4 _MS=5 _US=6 _LS=7 _TP=8 _FT=9 _BL=10 _SRF=11;

  ! zonal domains (applicable to the global 2D+ model only) - aligned to latitudes
    nzdom = 11 ! GL, NH, ARC, ETNH, IT, ETSH, ANT, SH, LAND, OCEAN, COAST
    allocate(zdom_mids(1:nzdom))
    allocate(zdom_bnds(1:nzdom+1))
  !                        GLOB   SH    ANT     ETSH    IT      ETNH    ARC   NH     LAND  OCEAN COAST
   !zdom_mids(:) = (/      -105., -95., -78.,   -45.,   0.,     +45.,   +78., +95.,  +102., +106., +110.   /)
   !zdom_bnds(:) = (/ -110., -110., -90., -66.56, -23.44, +23.44, +66.56, +90., +100., +104., +108., +112. /)
    zdom_mids(:) = (/      -120.,-100., -78.,   -45.,   0.,     +45.,   +78., +100., +120., +140., +160.   /)
    zdom_bnds(:) = (/ -130., -110., -90., -66.56, -23.44, +23.44, +66.56, +90., +110., +130., +150., +170. /)
  ! zdom_GLOB=1 _NH=2 _ARC=3 _ETNH=4 _IT=5 _ETSH=6 _ANT=7 _SH=8 _LAND=9 _OCEAN=10 _COAST=11;

  ! create new dimensions
  !
!!  call new_dimension(status, dimid_bnds, 'bnds', 2); _CCHSS_
    call new_dimension(status_dv, dimid_vdom, dim_vdom_NAME, nvdom)
    call channel_halt(substr, status_dv)
    call new_dimension(status_dz, dimid_zdom, dim_zdom_NAME, nzdom)
    call channel_halt(substr, status_dz)
  !
  ! vertical
    call add_dimension_variable(status, dim_vdom_NAME, dim_vdom_NAME, vdom_mids); _CCHSS_
    call add_dimension_variable_att(status, dim_vdom_NAME, dim_vdom_NAME, 'units', c='log10(press)'); _CCHSS_
    call add_dimension_variable_att(status, dim_vdom_NAME, dim_vdom_NAME, 'long_name', c='vertical domain'); _CCHSS_
    call add_dimension_variable_att(status, dim_vdom_NAME, dim_vdom_NAME, 'domains', c=vdom_info); _CCHSS_
    call add_dimension_variable_att(status, dim_vdom_NAME, dim_vdom_NAME, 'axis', c='Z'); _CCHSS_
    call add_dimension_variable_att(status, dim_vdom_NAME, dim_vdom_NAME, 'positive', c='down'); _CCHSS_
!!  call add_dimension_variable_att(status, dim_vdom_NAME, dim_vdom_NAME, 'bounds', c=dim_vdom_NAME//"_bnds")
!!  _CCHSS_
  !
  ! horizontal
    call add_dimension_variable(status, dim_zdom_NAME, dim_zdom_NAME, zdom_mids); _CCHSS_
    call add_dimension_variable_att(status, dim_zdom_NAME, dim_zdom_NAME, 'units', c='degrees_north'); _CCHSS_
    call add_dimension_variable_att(status, dim_zdom_NAME, dim_zdom_NAME, 'long_name', c='horizontal domain'); _CCHSS_
    call add_dimension_variable_att(status, dim_zdom_NAME, dim_zdom_NAME, 'domains', c=zdom_info); _CCHSS_
    call add_dimension_variable_att(status, dim_zdom_NAME, dim_zdom_NAME, 'axis', c='Y'); _CCHSS_
!!  call add_dimension_variable_att(status, dim_zdom_NAME, dim_zdom_NAME, 'bounds', c=dim_zdom_NAME//"_bnds")
!!  _CCHSS_

  ! create new representation(s)
  !
    call new_representation(status, reprid_domain &
           , trim(repr_name)                           &
           , rank = 2, link = 'xx--', dctype = DC_BC        &
           , dimension_ids = (/ dimid_zdom, dimid_vdom /)  &
           , ldimlen       = (/ nzdom,      nvdom           /)  &
!!         , nbounds       = (/ AUTO,       AUTO            /)  &
           , axis = 'ZY--'                                  &
           ); _CCHSS_
  !
!!  call new_representation(status, reprid_vdom_bnd &
!!         , trim(repr_name//'_VDOM_BND')                   &
!!         , rank = 2, link = '-x-x', dctype = DC_AG        &
!!         , dimension_ids = (/ dimid_vdom, dimid_bnds /)  &
!!         , ldimlen       = (/ AUTO,            AUTO            /)  &
!!         , nbounds       = (/ AUTO,            AUTO            /)  &
!!         , axis = '-Y-N'                                  &
!!         ); _CCHSS_
!!!
!!  call new_representation(status, reprid_zdom_bnd &
!!         , trim(repr_name//'_ZDOM_BND')                   &
!!         , rank = 2, link = '--xx', dctype = DC_AG        &
!!         , dimension_ids = (/ dimid_zdom, dimid_bnds /)  &
!!         , ldimlen       = (/ AUTO,            AUTO            /)  &
!!         , nbounds       = (/ AUTO,            AUTO            /)  &
!!         , axis = '--ZN'                                  &
!!         ); _CCHSS_

  ! parallel I/O ...
    nseg = 1
    allocate(start(nseg,IRANK))
    allocate(cnt(nseg,IRANK))
    allocate(meml(nseg,IRANK))
    allocate(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:)  = 1
    meml(:,:) = 1
    memu(:,:) = 1

    cnt(:,1)   = nvdom
    memu(:,1)  = nvdom
    cnt(:,2)   = nzdom
    memu(:,2)  = nzdom

    call set_representation_decomp(status, reprid_domain &
           , start, cnt, memu, meml, .false., piotype=PIOTYPE_COL)
    _CCHSS_

    deallocate(start) ; nullify(start)
    deallocate(cnt)   ; nullify(cnt)
    deallocate(meml)  ; nullify(meml)
    deallocate(memu)  ; nullify(memu)
  ! ... parallel I/O

  ! cleaup
    deallocate(vdom_mids)
    deallocate(vdom_bnds)
    deallocate(zdom_mids)
    deallocate(zdom_bnds)

    status = 0

  end subroutine domint_create_domain_repr_atm_VZ

! ----------------------------------------------------------------------

  subroutine domint_create_domain_repr_surf(status, reprid_domain)

    use messy_main_channel_bi,         only: DC_BC!, SCALAR
    use messy_main_channel_error_bi,   only: channel_halt
    use messy_main_channel_dimensions, only: get_dimension_info, new_dimension, &
                                             add_dimension_variable, add_dimension_variable_att
    use messy_main_channel_repr,       only: get_representation_info,new_representation, &
                                             set_representation_decomp, &
                                             IRANK, PIOTYPE_COL!, AUTO

    implicit none

    integer, intent(out)        :: status, reprid_domain

    ! LOCAL
    character(len=*), parameter :: substr = modstr//'_create_domain_repr_surf'
    character(len=*), parameter :: repr_name = 'GP_2D_DOMAIN_SURF', &
                                   dim_sdom_name = 'domain_surface'
    integer                     :: dimid_sdom!, reprid_sdom_bnd, dim_sdom_len
    integer                     :: status_ds

    ! PARALLEL DECOMPOSITION
    integer                          :: nseg = 0
    integer, dimension(:,:), pointer :: start => NULL()
    integer, dimension(:,:), pointer :: cnt   => NULL()
    integer, dimension(:,:), pointer :: meml  => NULL()
    integer, dimension(:,:), pointer :: memu  => NULL()

  ! adopted an example of messy_main_rnd_bi

  ! first check whether dimensions/representation exist
    call get_representation_info(status, repr_name, reprid_domain)
    call get_dimension_info(status_ds, dim_sdom_name, dimid_sdom, nsdom)
    if (.not.((status==2003).and.(status_ds==905))) then
      print *,trim(substr)//': error checking for dimension/representation for '//repr_name
      print *,'  status = ',status,'  status_ds = ',status_ds
      status = status+status_ds
      return
    endif

  ! horizontal domains (applicable to the global 2D+ model only) - aligned to latitudes
    nsdom = 11 ! GL, NH, ARC, ETNH, IT, ETSH, ANT, SH, LAND, OCEAN, COAST
    allocate(sdom_mids(1:nsdom))
    allocate(sdom_bnds(1:nsdom+1))
  !                        GLOB   SH    ANT     ETSH    IT      ETNH    ARC   NH     LAND  OCEAN COAST
   !sdom_mids(:) = (/      -105., -95., -78.,   -45.,   0.,     +45.,   +78., +95.,  +102., +106., +110.   /)
   !sdom_bnds(:) = (/ -110., -110., -90., -66.56, -23.44, +23.44, +66.56, +90., +100., +104., +108., +112. /)
    sdom_mids(:) = (/      -120.,-100., -78.,   -45.,   0.,     +45.,   +78., +100., +120., +140., +160.   /)
    sdom_bnds(:) = (/ -130., -110., -90., -66.56, -23.44, +23.44, +66.56, +90., +110., +130., +150., +170. /)
  ! sdom_GLOB=1 _NH=2 _ARC=3 _ETNH=4 _IT=5 _ETSH=6 _ANT=7 _SH=8 _LAND=9 _OCEAN=10 _COAST=11;

  ! create new dimensions
    call new_dimension(status_ds, dimid_sdom, dim_sdom_NAME, nsdom)
    call channel_halt(substr, status_ds)
  ! surface
    call add_dimension_variable(status, dim_sdom_NAME, dim_sdom_NAME, sdom_mids); _CCHSS_
    call add_dimension_variable_att(status, dim_sdom_NAME, dim_sdom_NAME, 'units', c='degrees_north'); _CCHSS_
    call add_dimension_variable_att(status, dim_sdom_NAME, dim_sdom_NAME, 'long_name', c='horizontal domain'); _CCHSS_
    call add_dimension_variable_att(status, dim_sdom_NAME, dim_sdom_NAME, 'domains', c=sdom_info); _CCHSS_
    call add_dimension_variable_att(status, dim_sdom_NAME, dim_sdom_NAME, 'axis', c='Y'); _CCHSS_
!!  call add_dimension_variable_att(status, dim_sdom_NAME, dim_sdom_NAME, 'bounds', c=dim_sdom_NAME//"_bnds")
!!  _CCHSS_

  ! create new representation(s)
  !
    call new_representation(status, reprid_domain &
           , trim(repr_name)                           &
           , rank = 1, link = 'x---', dctype = DC_BC        &
           , dimension_ids = (/ dimid_sdom /)  &
           , ldimlen       = (/ nsdom      /)  &
!!         , nbounds       = (/ AUTO       /)  &
           , axis = 'Y---'                                  &
           ); _CCHSS_
  !
!!  call new_representation(status, reprid_zdom_bnd &
!!         , trim(repr_name//'_SDOM_BND')                   &
!!         , rank = 2, link = 'x---', dctype = DC_AG        &
!!         , dimension_ids = (/ dimid_sdom, dimid_bnds /)  &
!!         , ldimlen       = (/ AUTO,       AUTO       /)  &
!!         , nbounds       = (/ AUTO,       AUTO       /)  &
!!         , axis = '--ZN'                                  &
!!         ); _CCHSS_

  ! parallel I/O ...
    nseg = 1
    allocate(start(nseg,IRANK))
    allocate(cnt(nseg,IRANK))
    allocate(meml(nseg,IRANK))
    allocate(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:)  = 1
    meml(:,:) = 1
    memu(:,:) = 1

    cnt(:,1)   = nsdom
    memu(:,1)  = nsdom

    call set_representation_decomp(status, reprid_domain &
           , start, cnt, memu, meml, .false., piotype=PIOTYPE_COL)
    _CCHSS_

    deallocate(start) ; nullify(start)
    deallocate(cnt)   ; nullify(cnt)
    deallocate(meml)  ; nullify(meml)
    deallocate(memu)  ; nullify(memu)
  ! ... parallel I/O

  ! cleaup
    deallocate(sdom_mids)
    deallocate(sdom_bnds)

    status = 0

  end subroutine domint_create_domain_repr_surf

! ----------------------------------------------------------------------

  subroutine domint_update_domain_masks_atm_VZ(work)

    use messy_main_grid_def_mem_bi,       only: nlev!, jrow, kproma

    implicit none

    type(t_domint_work), intent(inout) :: work ! domint work struct

    integer, parameter  :: &  ! domain custom indices
      vdom_ATM=1, &
      vdom_T100=2, vdom_T=3, vdom_TO3=4, &
      vdom_MS=5, vdom_US=6, vdom_LS=7, &
      vdom_TP=8, vdom_FT=9, vdom_BL=10, vdom_SRF=11, &
      zdom_GLOB=1, &
      zdom_SH=2, zdom_ANT=3, zdom_ETSH=4, zdom_IT=5, zdom_ETNH=6, zdom_ARC=7, zdom_NH=8, &
      zdom_LAND=9, zdom_OCEAN=10, zdom_COAST=11
    real(dp), parameter :: deg2rad = 4._dp*atan(1._dp)/180._dp ! degrees 2 radians conversion f-r
    real(dp), parameter :: sin_tc = sin(23.44_dp*deg2rad), &   ! tropic of cancer/capricorn
                           sin_pc = sin(66.56_dp*deg2rad)      ! polar circle
    integer             :: jd, jl, jr, jp, tp_il, bl_il

  ! already updated?
    if (work%lupdated) return

  ! invariable domains initialisation
    if (.not.work%comp%linvinit) then
    !
    ! initialising domain flags
      work%comp%mask(:,:,:,1) = .true.                 ! entire atmosphere (ATM-GLOB, seq. no. 1)
      work%comp%mask(:,:,:,2:work%comp%ndom) = .false. ! all other
    !
    ! calculating zonal flags in ATM, lev=1
    ! NH/SH
      work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_NH,1))   = (work%input%sinlat%p2(:,:).gt.0)
      work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_SH,1))   = (work%input%sinlat%p2(:,:).lt.0)
    ! ETNH/IT/SH
      work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_ETNH,1)) = (work%input%sinlat%p2(:,:).gt.sin_tc)
      work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_IT,1))   = (work%input%sinlat%p2(:,:).le.sin_tc) .and. &
                                                                    (work%input%sinlat%p2(:,:).ge.(-sin_tc))
      work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_ETSH,1)) = (work%input%sinlat%p2(:,:).lt.(-sin_tc))
    ! ARC/ANT
      work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_ARC,1))  = (work%input%sinlat%p2(:,:).gt.sin_pc)
      work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_ANT,1))  = (work%input%sinlat%p2(:,:).lt.(-sin_pc))
    ! LAND/OCEAN/COAST (warning: move this clause out of this if-endif block if land-sea mask is not invariant!)
      if ( work%input%landfrac%is ) then
        where( work%input%landfrac%p2(:,:).lt.thres_coast )
          work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_OCEAN,1)) = .true.
        endwhere
        where( work%input%landfrac%p2(:,:).ge.thres_coast )
          work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_COAST,1)) = .true.
        endwhere
        where( work%input%landfrac%p2(:,:).ge.thres_land )
          work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_COAST,1)) = .false.
          work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,zdom_LAND,1))  = .true.
        endwhere
      endif
    !
    ! + propagating to all levels
      do jl=2, nlev
        work%comp%mask(:,jl,:,work%comp%seqno(vdom_ATM,:,1)) = work%comp%mask(:,1,:,work%comp%seqno(vdom_ATM,:,1))
      enddo
    !
    ! flag invariant masks initialised
      work%comp%linvinit = .true.
    endif

  ! ----- first populating zonal flags from ATM to remaining vertical domains -----
    do jd=2, nvdom
      work%comp%mask(:,1:nlev,:,work%comp%seqno(jd,:,1)) = work%comp%mask(:,1:nlev,:,work%comp%seqno(vdom_ATM,:,1))
    enddo
  ! defining surface
    work%comp%mask(:,1:nlev-1,:,work%comp%seqno(vdom_SRF,:,1)) = .false.

  ! ----- then modulating vertical domain flags according to available data -----
    do jd=1, nzdom

    ! pressure-based
    ! mesosphere, MS, above 1hPa
      where( work%input%press%p3(:,:,:).lt.1e2_dp )
      elsewhere
        work%comp%mask(:,:,:,work%comp%seqno(vdom_MS,jd,1)) = .false.
      endwhere
    ! upper strat, US, 100hPa-1hPa
      where( ( work%input%press%p3(:,:,:) .ge. 1e2_dp ).and.( work%input%press%p3(:,:,:) .lt. 100e2_dp ) )
      elsewhere
        work%comp%mask(:,:,:,work%comp%seqno(vdom_US,jd,1)) = .false.
      endwhere
    ! T100, troposphere below 100hPa
      where( work%input%press%p3(:,:,:) .ge. 100e2_dp )
      elsewhere
        work%comp%mask(:,:,:,work%comp%seqno(vdom_T100,jd,1)) = .false.
      endwhere
    ! O3-based tropoposphere ( stratosphere at O3>=150 ppb, p<500 hPa )
      if ( work%input%ozone%is ) then  ! O3 mixing ratio-based troposphere
        where( (( work%input%ozone%p3(:,:,:) .ge. 150e-9_dp ).and.( work%input%press%p3(:,:,:) .lt. 500e2_dp )) &
              .or.( work%input%press%p3(:,:,:) .lt. 10e2_dp ) )
          work%comp%mask(:,:,:,work%comp%seqno(vdom_TO3,jd,1)) = .false.
        endwhere
      else
        work%comp%mask(:,:,:,work%comp%seqno(vdom_TO3,jd,1)) = .false.       ! set work. off
      endif

    ! domains using tropopause/bundary layer heights
      if ( work%input%tph_idx%is .or. work%input%blh_idx%is ) then
        do jr=1, size(work%comp%mask,3)   !   row, sx in pdef SMCL
          do jp=1, size(work%comp%mask,1) ! proma, s1 in pdef SMCL

          ! tropopause/boundary-layer index-based
            if ( work%input%tph_idx%is ) then
              tp_il = work%input%tph_idx%p2(jp,jr)  ! local tropopause layer index
            ! tropopause
              work%comp%mask(jp,1:tp_il-1,   jr,work%comp%seqno(vdom_TP,jd,1)) = .false.
              work%comp%mask(jp,tp_il+1:nlev,jr,work%comp%seqno(vdom_TP,jd,1)) = .false.
            ! troposphere
              work%comp%mask(jp,1:tp_il-1,   jr,work%comp%seqno(vdom_T,jd,1)) = .false.

              if ( work%input%theta%is ) then  ! lower stratosphere (380K downto tropopause)
                where( work%input%theta%p3(jp,1:tp_il,jr).gt.380_dp )
                  work%comp%mask(jp,1:tp_il,jr,work%comp%seqno(vdom_LS,jd,1)) = .false.    ! theta above 380K
                endwhere
                work%comp%mask(jp,tp_il+1:nlev,jr,work%comp%seqno(vdom_LS,jd,1)) = .false. ! TP to surface
              endif
            endif

            if ( work%input%blh_idx%is ) then
              bl_il = work%input%blh_idx%p2(jp,jr)  ! local boundary layer index
            ! boundary layer
              work%comp%mask(jp,1:bl_il-1,jr,work%comp%seqno(vdom_BL,jd,1)) = .false.
            ! free troposphere
              if ( work%input%tph_idx%is ) then
                work%comp%mask(jp,1:tp_il-1, jr,work%comp%seqno(vdom_FT,jd,1)) = .false.
                work%comp%mask(jp,bl_il:nlev,jr,work%comp%seqno(vdom_FT,jd,1)) = .false.
              endif
            endif

          enddo ! jp
        enddo  ! jr
      endif  ! if (l_tph.or.l_blh) then

    ! nullifying domains for which no data is available
      if ( .not.work%input%tph_idx%is ) then
        work%comp%mask(:,:,:,work%comp%seqno(vdom_LS,jd,1)) = .false.       ! set LS work. off
        work%comp%mask(:,:,:,work%comp%seqno(vdom_TP,jd,1)) = .false.       ! set TP work. off
        work%comp%mask(:,:,:,work%comp%seqno(vdom_T,jd,1))  = .false.       ! set T work. off
        work%comp%mask(:,:,:,work%comp%seqno(vdom_FT,jd,1)) = .false.       ! set FT work. off
      endif
      if ( .not.work%input%theta%is ) &
        work%comp%mask(:,:,:,work%comp%seqno(vdom_LS,jd,1)) = .false.       ! set LS work. off
      if ( .not.work%input%blh_idx%is ) then
        work%comp%mask(:,:,:,work%comp%seqno(vdom_FT,jd,1)) = .false.       ! set FT work. off
        work%comp%mask(:,:,:,work%comp%seqno(vdom_BL,jd,1)) = .false.       ! set BL work. off
      endif

    enddo

  ! flag ad updated
    work%lupdated = .true.

  end subroutine domint_update_domain_masks_atm_VZ

! ----------------------------------------------------------------------

  subroutine domint_update_domain_masks_surf(work)

    use messy_main_grid_def_mem_bi,       only: nlev!, jrow, kproma

    implicit none

    type(t_domint_work), intent(inout) :: work ! domint work struct

    integer, parameter  :: &  ! domain custom indices
      sdom_GLOB=1, &
      sdom_SH=2, sdom_ANT=3, sdom_ETSH=4, sdom_IT=5, sdom_ETNH=6, sdom_ARC=7, sdom_NH=8, &
      sdom_LAND=9, sdom_OCEAN=10, sdom_COAST=11
    real(dp), parameter :: deg2rad = 4._dp*atan(1._dp)/180._dp ! degrees 2 radians conversion f-r
    real(dp), parameter :: sin_tc = sin(23.44_dp*deg2rad), &   ! tropic of cancer/capricorn
                           sin_pc = sin(66.56_dp*deg2rad)      ! polar circle
    integer             :: jd, jl, jr, jp, tp_il, bl_il

  ! already updated?
    if (work%lupdated) return

  ! invariable domains initialisation
    if (.not.work%comp%linvinit) then
    !
    ! initialising domain flags
      work%comp%mask(:,:,1,1) = .true.                 ! whole model domain (GLOB, seq. no. 1)
      work%comp%mask(:,:,1,2:work%comp%ndom) = .false. ! all other
    !
    ! calculating zonal flags in ATM, lev=1
    ! NH/SH
      work%comp%mask(:,:,1,work%comp%seqno(sdom_NH,1,1))   = (work%input%sinlat%p2(:,:).gt.0)
      work%comp%mask(:,:,1,work%comp%seqno(sdom_SH,1,1))   = (work%input%sinlat%p2(:,:).lt.0)
    ! ETNH/IT/SH
      work%comp%mask(:,:,1,work%comp%seqno(sdom_ETNH,1,1)) = (work%input%sinlat%p2(:,:).gt.sin_tc)
      work%comp%mask(:,:,1,work%comp%seqno(sdom_IT,1,1))   = (work%input%sinlat%p2(:,:).le.sin_tc) .and. &
                                                             (work%input%sinlat%p2(:,:).ge.(-sin_tc))
      work%comp%mask(:,:,1,work%comp%seqno(sdom_ETSH,1,1)) = (work%input%sinlat%p2(:,:).lt.(-sin_tc))
    ! ARC/ANT
      work%comp%mask(:,:,1,work%comp%seqno(sdom_ARC,1,1))  = (work%input%sinlat%p2(:,:).gt.sin_pc)
      work%comp%mask(:,:,1,work%comp%seqno(sdom_ANT,1,1))  = (work%input%sinlat%p2(:,:).lt.(-sin_pc))

    ! LAND/OCEAN/COAST (warning: move this clause out of this if-endif block if land-sea mask is not invariant!)
      if ( work%input%landfrac%is ) then
        work%comp%mask(:,:,1,work%comp%seqno(sdom_LAND,1,1))  = (work%input%landfrac%p2(:,:).ge.thres_land)
        work%comp%mask(:,:,1,work%comp%seqno(sdom_OCEAN,1,1)) = (work%input%landfrac%p2(:,:).lt.thres_coast)
        work%comp%mask(:,:,1,work%comp%seqno(sdom_COAST,1,1)) = &
          .not.( work%comp%mask(:,:,1,work%comp%seqno(sdom_LAND,1,1)) .or. &
                 work%comp%mask(:,:,1,work%comp%seqno(sdom_OCEAN,1,1)) )
      endif
    !
    ! flag invariant masks initialised
      work%comp%linvinit = .true.
    endif

  ! flag ad updated
    work%lupdated = .true.

  end subroutine domint_update_domain_masks_surf

! ----------------------------------------------------------------------

  subroutine domint_global_start

    use messy_main_timer_bi, only: event_state
    use messy_main_timer,    only: next_date

    implicit none

    character(len=*), parameter :: substr = modstr//'_global_start'
    integer :: jd

  ! update calculation triggers from timer
    do jd = 1, domint_ndom
      if ( .not.domain(jd)%ldisable ) work(jd)%ltrigcalc = event_state(work(jd)%event, next_date)
    enddo



  end subroutine domint_global_start

! ----------------------------------------------------------------------

  subroutine domint_global_end

    use messy_main_mpi_bi, only: &
#ifdef DEBUG
                                 p_parallel_io, &
#endif
                                 p_sum

    implicit none

    character(len=*), parameter :: substr = modstr//'_global_end'
    integer :: status
    integer :: jd, jt
    real(dp), dimension(:,:), allocatable :: itd

  ! integrate domain-representation-wise
    domain_loop: do jd = 1, domint_ndom

    ! time to integrate?
      if (.not.work(jd)%ltrigcalc) cycle

      work(jd)%lupdated = .false.

      select case(jd)

      case(domint_atm_VZ) ! atmosphere, vertical-zonal
        call domint_update_domain_masks_atm_VZ(work(jd))

      case(domint_surf) ! surface
        call domint_update_domain_masks_surf(work(jd))

      case default
        call error_bi('unknown domain kind encountered',substr)

      end select

#ifdef DEBUG
      if ( p_parallel_io ) &
        print *,'domint_global_end(',jd,'): integration'
#endif

    ! integrate
      call domint_integrate(status, work(jd)%comp, work(jd)%tasks)

    ! sum/broadcast result over all PEs
      allocate(itd(size(work(jd)%tasks),work(jd)%comp%ndom))
    ! optimisation: gather all tasks results -> p_sum -> distribute
      do jt = 1, size(work(jd)%tasks)
        itd(jt,:) = (work(jd)%tasks(jt)%int_dom(:))
      enddo
      itd(:,:) = p_sum(itd(:,:))
      do jt = 1, size(work(jd)%tasks)
        work(jd)%tasks(jt)%int_dom(:) = itd(jt,:)
      enddo
      deallocate(itd)

    ! cells count
      work(jd)%comp%cells(:) = p_sum(work(jd)%comp%cells(:))

    ! put missing flags for domains with zero cells
      do jt = 1, size(work(jd)%tasks)
        where ( work(jd)%comp%cells(:) .lt. 1 )
          work(jd)%tasks(jt)%int_dom(:) = FLAG_UNDEF
        endwhere
      enddo

    enddo domain_loop

  end subroutine domint_global_end

! ----------------------------------------------------------------------

  subroutine domint_free_memory
    implicit none
    integer :: jd

    do jd = 1, domint_ndom
      call domint_deallocate(work(jd)%comp, work(jd)%tasks)
    enddo

  end subroutine domint_free_memory

! ----------------------------------------------------------------------

  subroutine domint_read_nml_cpl(status, iou)
    use messy_main_tools, only: read_nml_open, read_nml_check, read_nml_close
    implicit none
    integer, intent(out) :: status   ! error status
    integer, intent(in)  :: iou      ! I/O unit
    character(len=*), parameter :: substr = modstr//'_read_nml_cpl'
    namelist /cpl/ domint, domain, thres_land, thres_coast
    logical              :: lex      ! file exists ?
    integer              :: fstat    ! file status
    status = 1
    call read_nml_open(lex, substr, iou, 'cpl', modstr)
    if (.not.lex) return             ! <modstr>.nml does not exist
    read(iou, NML=cpl, iostat=fstat)
    call read_nml_check(fstat, substr, iou, 'cpl', modstr)
    if (fstat /= 0) return           ! error while reading namelist
    call read_nml_close(substr, iou, modstr)
    status = 0                       ! no ERROR
  end subroutine domint_read_nml_cpl

! ----------------------------------------------------------------------

#endif

end module messy_domint_si
