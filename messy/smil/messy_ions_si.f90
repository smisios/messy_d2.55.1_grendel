#include "messy_main_ppd_bi.inc" 
! defines ranks and indices for different base models
! Submodel interface

!!#> ... !#< today, to continue construction after interruption

module messy_ions_si

  use messy_main_constants_mem, only: wp=>dp, dp, strlen_medium
  use messy_ions
  use messy_main_channel, only: t_chaobj_cpl
  use messy_main_tools,   only: PTR_1D_ARRAY

  implicit none
  private
  save

  public :: ions_initialize
!  public :: ions_new_tracer
  public :: ions_init_memory
  public :: ions_init_coupling
  public :: ions_physc
  public :: ions_free_memory

  integer :: idx_si(2) ! id for tracer
  integer :: idx_Rndecay(n_Rndec) ! id for Radon decay tracers
  logical :: lions

  ! pointer
  real(kind=wp), dimension(:,:,:), pointer :: total_ipr => null() ! total ion pair producation rate
  real(kind=wp), dimension(:,:,:), pointer :: small_ions_neg => null() ! small positive ion concentration (cm-3)
  real(kind=wp), dimension(:,:,:), pointer :: small_ions_pos => null() ! small negative ion concentration (cm-3)
  real(kind=wp), dimension(:,:,:), pointer :: krec => null()  ! ion ion recombination rate constant
  real(kind=wp), dimension(:,:,:), pointer :: radon_ipr => null() ! ion pair production rate due to Radon decay  
  real(kind=wp), dimension(:,:,:), pointer :: gcr_ipr => null()  ! ion pair production rate due to GCR
  real(kind=wp), dimension(:,:,:), pointer :: aero_cs => null()  ! aerosol ion interaction, condensation rate of ions to aerosol (cm3 s-1)
  real(kind=wp), dimension(:,:,:,:), pointer :: aero_r => null() ! aerosol radius (m)
  real(kind=wp), dimension(:,:,:,:), pointer :: aerocon => null() ! aerosol number concentration (cm-3)
  real(kind=wp), dimension(:,:), pointer :: pc_cha => null()  ! geomagnetic cut off rigidity, measure of shielding potential against GCR
!  real(kind=wp), dimension(:,:,:,:), pointer :: Rn_decay => null() ! dimension must match the dimension in dradon_si

  ! strings for coupling
  character(len=strlen_medium), public :: Rn_cpl(2) = '' ! Radon decay chain
  character(len=5), dimension(n_Rndec), parameter:: decay_chain_name = &
           & (/'Rn222', 'Po218', 'Pb214', 'Bi214', 'Pb210'/)

  ! coupling
  type(t_chaobj_cpl) :: igrf   ! coupling for the geomagnetic field time series
  type(t_chaobj_cpl) :: phi    ! coupling for the GCR modulation time series

  real(dp), dimension(:), pointer, public :: phi_now => null() ! GCR modulation at current time, linear interpolation
  real(dp), dimension(:), pointer, public :: igrf_now => null() ! parameters of geomagnetic field current time, linera interpolation

  ! coupling to GEC or other models
  type(t_chaobj_cpl) :: cpl_ipr_Rn
  type(t_chaobj_cpl) :: cpl_ipr_gcr
  type(t_chaobj_cpl) :: cpl_ipr_total
  type(t_chaobj_cpl) :: cpl_ionc      ! take or provide ion concentration from/to other submodel
  type(t_chaobj_cpl) :: cpl_aero_cs   ! aerosol-ion interaction uptake of small ions to aerosol particles
  type(t_chaobj_cpl) :: cpl_aero_r    ! aerosol radius (cm)     
  type(t_chaobj_cpl) :: cpl_aerocon   ! aerosol concentration (cm-3)
  logical, public    :: lssion = .true.  ! if true calculate steady state ion concentrations
  logical, public    :: ltotalipr  = .true.  ! if true calculate the total ion pair production rate

  integer :: si_idx(2) ! small ions tracer index

  ! Look up tables
  type(ptr_1d_array),           dimension(:), pointer :: dimaxe  => null()
  integer,                      dimension(:), pointer :: dimlen  => null()
  character(len=strlen_medium), dimension(:), pointer :: dimname => null()
  character(len=strlen_medium), dimension(:), pointer :: dimunit => null()

contains

subroutine ions_initialize
  ! read in namelist files and initialise all global variables

  use messy_main_mpi_bi,     only: p_parallel_io, p_io, p_bcast
  use messy_main_blather_bi, only: error_bi, start_message_bi, end_message_bi
  use messy_main_tools,      only: find_next_free_unit
  use messy_main_import_lt,  only: get_lookup_table

  implicit none

  ! LOCAL
  character(len=*), parameter :: substr = 'ions_initialize'
  integer :: iou    ! I/O unit
  integer :: status ! error status
  integer :: ltrank ! for look up table
  integer :: jr     ! rank index for look up table

  ! Read name list file
  if(p_parallel_io) then
    iou = find_next_free_unit(100,200)  ! what is this doing?
    call read_ionisation_nml(status, iou)
    if(status /= 0) call error_bi(' ',substr)
    if(status == 0) lions = .true.
  end if

  call p_bcast(lions, p_io)
  call p_bcast(gcr_method, p_io)
  call p_bcast(lqradon, p_io)
  call p_bcast(lgcr, p_io)
  call p_bcast(laero, p_io)


  if(p_parallel_io) then
    iou = find_next_free_unit(100,200)  
    call ions_read_nml_cpl(status, iou)
    if(status /= 0) call error_bi(' ',substr)
  end if

!  call p_bcast(Rd_cpl, p_io) ! Radon currently fixed,
  CALL p_bcast(igrf%cha, p_io) 
  CALL p_bcast(igrf%obj, p_io) 
  CALL p_bcast(phi%cha, p_io) 
  CALL p_bcast(phi%obj, p_io) 

  CALL p_bcast(lssion,  p_io) 
  CALL p_bcast(ltotalipr,   p_io)

  CALL p_bcast(cpl_ipr_Rn%cha,  p_io)
  CALL p_bcast(cpl_ipr_gcr%cha, p_io) 
  CALL p_bcast(cpl_ipr_total%cha, p_io)
  CALL p_bcast(cpl_ionc%cha, p_io)
  CALL p_bcast(cpl_aero_cs%cha, p_io)
  CALL p_bcast(cpl_aero_r%cha, p_io)
  CALL p_bcast(cpl_aerocon%cha, p_io)

  CALL p_bcast(cpl_ipr_Rn%obj,  p_io)
  CALL p_bcast(cpl_ipr_gcr%obj, p_io) 
  CALL p_bcast(cpl_ipr_total%obj, p_io)
  CALL p_bcast(cpl_ionc%obj, p_io)
  CALL p_bcast(cpl_aero_cs%obj, p_io)
  CALL p_bcast(cpl_aero_r%obj, p_io)
  CALL p_bcast(cpl_aerocon%obj, p_io)


  ! read the look up table if 
  CALL start_message_bi(modstr, 'LOOKUP TABLE INITIALISATION', substr)

  CALL get_lookup_table(status             &
       , 'CRII', p3=CRII_TABLES, rank=ltrank        &
       , dimlen=dimlen, dimaxe=dimaxe, dimname=dimname, dimunit=dimunit)

  if(status /= 0) CALL error_bi('CRII lookup table not available',substr)
  if(ltrank /= 3) CALL error_bi('CRII lookup table has rank /= 3',substr)

  do jr=1, ltrank
    select case(trim(dimname(jr)))
      case("H")
        if (jr /= 1) CALL error_bi('CRII LT; atmospheric depth H not rank 1',substr)
        H => dimaxe(jr)%ptr(:) 
      case("PHI")
        if (jr /= 2) CALL error_bi('CRII LT; PHI not rank 2',substr)
!        if 
        PHI_AV => dimaxe(jr)%ptr(:)
      case("PC")
        if(jr /= 3) CALL error_bi('CRII LT; PC not rank 3 is ', substr)
        PC => dimaxe(jr)%ptr(:)
      case default
        CALL error_bi('CRII lookup table with unknonw axis: '// trim(dimname(jr)), substr)
    end select
  end do

  CALL end_message_bi(modstr, 'LOOKUP TABLE INITIALISATION', substr)


end subroutine ions_initialize



! subroutine ions_new_tracer
!   ! IONS ARE CURRENTLY NOT CONSIDERED TRACERS
!   ! initialise the two small ion tracers
!   ! small_ions_neg : negative small ions
!   ! small_ions_pos : positive small ions
!   use messy_main_mpi_bi,        only: p_parallel_io
!   use messy_main_tracer_mem_bi, only: GPTRSTR
!   use messy_main_tracer,        only: new_tracer, set_tracer
! 
!   implicit none
!   ! Local
! 
! end subroutine ions_new_tracer


subroutine ions_init_memory
  ! Here we will start creating memory leaks ;)
  ! actually we define all channels needed to make the small ion concentration 
  ! and ion pair production rate known to all other models
  ! Base model interface
  use messy_main_blather_bi,       only: start_message_bi, end_message_bi
  use messy_main_mpi_bi,           only: p_parallel_io
  use messy_main_channel_error_bi, only: channel_halt
  use messy_main_channel_bi,       only: GP_3D_MID, GP_2D_HORIZONTAL 
  ! MESSy
  use messy_main_channel,          only: new_channel, new_channel_object &
                                       , new_attribute

  implicit none

  ! local
  integer :: status
  character(len=*), parameter :: substr="ions_init_memory"

  if(.not. lions) return 

  call start_message_bi(modstr, "CHANNEL DEFINITION", substr)

  ! define new channel
  call new_channel(status, modstr, reprid=GP_3D_MID)
  call channel_halt(substr, status)

  if(trim(cpl_ipr_total%cha)=='') then
    if(p_parallel_io) write(*,*) ' ... total ion pair production rate'
    call new_channel_object(status, modstr, 'total_ipr', p3 = total_ipr, lrestreq=.TRUE.)
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'total_ipr', 'long_name', c='Total Ion Pair Production Rate' )
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'total_ipr', 'units', c='cm-3 s-1' )
    call channel_halt(substr, status)
  end if 

  if(trim(cpl_ionc%cha)=='') then 
    if(p_parallel_io) write(*,*) ' ... small_ions_neg'
    call new_channel_object(status, modstr, 'small_ions_neg', p3 = small_ions_neg, lrestreq=.TRUE.)
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'small_ions_neg', 'long_name', c='Concentration of small ions' )
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'small_ions_neg', 'units', c='cm-3' )
    call channel_halt(substr, status)
  end if

!  if(trim(cpl_ionc%cha)=='') then ! other models don't distinguish between positive and negative ions
    if(p_parallel_io) write(*,*) ' ... small_ions_pos'
    call new_channel_object(status, modstr, 'small_ions_pos', p3 = small_ions_pos, lrestreq=.TRUE.)
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'small_ions_pos', 'long_name', c='Concentration of small ions' )
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'small_ions_pos', 'units', c='cm-3' )
    call channel_halt(substr, status)
!  end if

  if(p_parallel_io) write(*,*) ' ... recombination rate constant'
  call new_channel_object(status, modstr, 'krec', p3 = krec, lrestreq=.TRUE.)
  call channel_halt(substr, status)
  call new_attribute(status, modstr, 'krec', 'long_name', c='Recombination rate constant' )
  call channel_halt(substr, status)
  call new_attribute(status, modstr, 'krec', 'units', c='cm3 s-1' )
  call channel_halt(substr, status)

  if(lqradon .and. trim(cpl_ipr_Rn%cha)=='' ) then
    if(p_parallel_io) write(*,*) '... ion pair production from Radon decay' 
    call new_channel_object(status, modstr, 'radon_ipr', p3 = radon_ipr, lrestreq=.TRUE.)
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'radon_ipr', 'long_name', c='Ion Pair Production Rate from Radon decay' )
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'radon_ipr', 'units', c='cm-3 s-1' )
    call channel_halt(substr, status)
  end if

  if(lgcr .and. trim(cpl_ipr_gcr%cha)=='' ) then
    if(p_parallel_io) write(*,*) '... ion pair production from Radon decay'
    call new_channel_object(status, modstr, 'gcr_ipr', p3 = gcr_ipr, lrestreq=.TRUE.)
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'gcr_ipr', 'long_name', c='Ion Pair Production Rate from GCR' )
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'gcr_ipr', 'long_name', c='cm-3 s-1')
  end if

  if(laero .and. trim(cpl_aero_cs%cha)=='') then
    if(p_parallel_io) write(*,*) '... ion uptake by aerosol particles'
    call new_channel_object(status, modstr, 'aero_cs', p3 = aero_cs, lrestreq=.TRUE.)
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'aero_cs', 'long_name', c='Ion uptake by aerosol particles')
    call channel_halt(substr, status)
    call new_attribute(status, modstr, 'aero_cs', 'long_name', c='cm3 s-1')
  end if


  if(p_parallel_io) write(*,*) ' ... geomagnetic cut off rigidity'
  call new_channel(status, modstr//'_2d', reprid=GP_2D_HORIZONTAL)
  call channel_halt(substr, status)

  call new_channel_object(status, modstr//'_2d', 'PC', p2 = pc_cha, lrestreq=.TRUE.)
  call channel_halt(substr, status)
  call new_attribute(status, modstr//'_2d', 'PC', 'long_name', c='Geomagnetic cut off rigidity' )
  call channel_halt(substr, status)
  call new_attribute(status, modstr//'_2d', 'PC', 'units', c='GV' )
  call channel_halt(substr, status)


  call end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

end subroutine ions_init_memory


subroutine ions_init_coupling
  ! set the pointers to variables from other models.
  ! for this model the radon concentration is most important.
  use messy_main_blather_bi,       only: start_message_bi, end_message_bi
  use messy_main_channel,          only: get_channel_object, get_channel_info  &
                                       , get_channel_object_info               &
                                       , get_channel_object_dimvar
  use messy_main_channel_error_bi, only: channel_halt 
  use messy_main_tracer,           only: get_tracer
  use messy_main_mpi_bi,           only: p_parallel_io
  use messy_main_blather_bi,       only: error_bi, info_bi
  use messy_main_tracer_mem_bi,    only: ti_gp, GPTRSTR
  use messy_main_tools,            only: PTR_1D_ARRAY
  use messy_main_constants_mem,    only: STRLEN_ULONG

  implicit none

  ! local variables here only
  character(len=*), parameter :: substr="ions_init_coupling"
  integer :: status, ierr

  type(PTR_1D_ARRAY), dimension(:), pointer          :: dvs => null()  
  character(len=STRLEN_ULONG), dimension(:), pointer :: units => null()

  integer :: i


  call start_message_bi(modstr, 'COUPLING', substr)

  do i=1, size(decay_chain_name)
    call get_tracer(ierr, GPTRSTR, trim(decay_chain_name(i)), idx=idx_Rndecay(i))
  end do

! coupling to other ion models
  if(trim(cpl_ionc%cha)/='') then
    call info_bi('Coupling ion concentration: '//cpl_ionc%cha//'  '//cpl_ionc%obj)
    call get_channel_object(status, &
        & trim(cpl_ionc%cha), trim(cpl_ionc%obj), p3=small_ions_neg)
    call channel_halt(substr, status)
  end if 

  if(trim(cpl_ipr_gcr%cha)/='') then
    call info_bi('Coupling ion concentration: '//cpl_ipr_gcr%cha//'  '//cpl_ipr_gcr%obj)
    call get_channel_object(status, &
        & trim(cpl_ipr_gcr%cha), trim(cpl_ipr_gcr%obj), p3=gcr_ipr)
    call channel_halt(substr, status)
  end if

  if(trim(cpl_ipr_Rn%cha)/='') then
    call info_bi('Coupling ion concentration: '//cpl_ipr_Rn%cha//'  '//cpl_ipr_Rn%obj)
    call get_channel_object(status, &
        & trim(cpl_ipr_Rn%cha), trim(cpl_ipr_Rn%obj), p3=radon_ipr)
    call channel_halt(substr, status)
  end if

  if(trim(cpl_ipr_total%cha)/='') then
    call info_bi('Coupling ion concentration: '//cpl_ipr_total%cha//'  '//cpl_ipr_total%obj)
    call get_channel_object(status, &
        & trim(cpl_ipr_total%cha), trim(cpl_ipr_total%obj), p3=total_ipr)
    call channel_halt(substr, status)
  end if

! coupling to solar cylce and GCR 
  if (trim(phi%cha) /= '') then
    call info_bi('Looking for GCR DATA ... ')
    call info_bi('       channel: '//phi%cha)
    call info_bi('       object : '//phi%obj)
   
    call get_channel_object(status &
        & , trim(phi%cha), trim(phi%obj), p1=phi_now)
    call channel_halt(substr, status)

    call get_channel_object_dimvar(status &
        & , trim(phi%cha), trim(phi%obj) &
        & , dvs, units)
    call channel_halt(substr, status)

    do i=1, size(dvs)
       if (p_parallel_io) &
         & write(*,*) 'DIMVAR ',i,' [',trim(units(i)),']: ',dvs(i)%ptr
    end do

  endif


! coupling to geomagnetic field 
  if (trim(igrf%cha) /= '') then
    call info_bi('Looking for Geomagnetic DATA ... ')
    call info_bi('       channel: '//igrf%cha)
    call info_bi('       object : '//igrf%obj)
   
    call get_channel_object(status &
        & , trim(igrf%cha), trim(igrf%obj), p1=igrf_now)
    call channel_halt(substr, status)

    call get_channel_object_dimvar(status &
        & , trim(igrf%cha), trim(igrf%obj) &
        & , dvs, units)
    call channel_halt(substr, status)

    do i=1, size(dvs)
       if (p_parallel_io) &
         & write(*,*) 'DIMVAR ',i,' [',trim(units(i)),']: ',dvs(i)%ptr
    end do

   endif


! coupling to aerosol model
   if(trim(cpl_aero_cs%cha) /= '') then
     call info_bi('Take ION UPTAKE by aerosols from ...')
     call info_bi('       channel: '//cpl_aero_cs%cha)
     call info_bi('       object : '//cpl_aero_cs%obj)

     call get_channel_object(status, trim(cpl_aero_cs%cha), trim(cpl_aero_cs%obj), p3=aero_cs) 
     call channel_halt(substr, status)
   end if

   if(trim(cpl_aero_r%cha) /= '') then
     call info_bi('Take Aerosol radius from ...')
     call info_bi('       channel: '//cpl_aero_r%cha)
     call info_bi('       object : '//cpl_aero_r%obj)

     call get_channel_object(status, trim(cpl_aero_r%cha), trim(cpl_aero_r%obj), p4=aero_r)
     call channel_halt(substr, status)
   end if

   if(trim(cpl_aerocon%cha) /= '') then
     call info_bi('Take aerosols number concentration from ...')
     call info_bi('       channel: '//cpl_aerocon%cha)
     call info_bi('       object : '//cpl_aerocon%obj)

     call get_channel_object(status, trim(cpl_aerocon%cha), trim(cpl_aerocon%obj), p4=aerocon)
     call channel_halt(substr, status)
   end if

   call end_message_bi(modstr, 'COUPLING', substr)

end subroutine ions_init_coupling


subroutine ions_physc
  ! calculate ion physics here
  ! determine: 
  ! 1.) ion pair production rate from Radon decay
  ! 2.) ion pair production rate from GCR
  ! 3.) total ion pair production rate
  ! 4.) recombination rate constant
  ! 5.) steady state small ion concentration

  use messy_ions, only : usoskin_crii, recom_brasseur, steady_state_ions

  use messy_main_constants_mem, only : k_B, pi, Rdry => rd
  use messy_main_data_bi,       only : tm1, tte_3d, press_3d ! mid-level pressures [Pa]
  USE messy_main_grid_def_mem_bi, ONLY:jrow, kproma, nlev
  USE messy_main_grid_def_bi,     ONLY: philat_2d, philon_2d
  USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1
  use messy_main_timer,           only : tmst => time_step_len, year, month, day

  implicit none

  real(kind=wp) :: radlon(kproma), radlat(kproma)
  real(kind=wp) :: m(kproma,nlev), hgt(kproma, nlev)
  real(kind=wp) :: zpress(kproma,nlev), ztemp(kproma,nlev)
  real(kind=wp) :: mass2vol(kproma,nlev) ! conversion of g-1 to cm-3
  real(kind=wp) :: Rn_decay(n_Rndec)
  ! only local and act on pointers defined above
  integer :: jk, jp, jm, i
  ! needed for Usoskin parameterisation
  real(kind=wp) :: B_0, D, L(3), EMF, X(3) ! PHI_NOW, GH_NOW(8),
  ! local for the first start up

  ! initialise values with 0.0_dp if calculated by the ions submodel 
  if(associated(radon_ipr)) radon_ipr(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0.0_wp 
  if(associated(gcr_ipr))   gcr_ipr(_RI_XYZ__(1:kproma,jrow,1:nlev))   = 0.0_wp
  if(associated(total_ipr)) total_ipr(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0.0_wp
  if(associated(aero_cs))   aero_cs(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0.0_wp 
  if(associated(krec)) krec(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0.0_wp
  pc_cha(1:kproma, jrow) = 0.0_wp
!  write(*,*) "radon_ipr and gcr_ipr set to zero"

  ! 1) --- some calculations before -----------------------------------------------

!mz_se<
 ! aero_r(nproma, nlev, nmod, ngpblks)
 ! dim     1,     2,    3,    4,
  if(laero) then
    do jm =1,ubound(aero_r,3) ! nmod
      do jk=1, nlev ! ubound(aero_r,2) ?
        do jp=1, kproma ! ubound(aeror_r,1) ?
          aero_cs(_RI_XYZ__(jp,jrow,jk)) = aero_cs(_RI_XYZ__(jp,jrow,jk)) + &
                              & ion_aerosol_TZ06(aero_r(_RI_XYZN_(jp,jrow,jk,jm))*1.0e6 ) * aerocon(_RI_XYZN_(jp,jrow,jk,jm))
                                                     !EMAC: jp, jk,  jm, jrow
        end do
      end do
    end do
  end if 
!mz_se>

  if(lgcr .or. lqradon .or. ltotalipr) then
    ztemp(1:kproma,1:nlev) = tm1(_RI_XYZ__(1:kproma,jrow,1:nlev)) + &
                         & tte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) * tmst
    zpress(1:kproma,1:nlev) = press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))

    ! calculate total number density and convert from #/m3 to #/cm3
    m = zpress/(k_B*ztemp)*1e-6_wp ! number density in # cm-3
    radlat = philat_2d(1:kproma,jrow)/180.0_wp*pi ! latitude in rad
    radlon = philon_2d(1:kproma,jrow)/180.0_wp*pi ! longitude in rad

    if(lgcr .and. gcr_method==2) then 
      hgt = 0.010195_wp*zpress
      mass2vol = m*4.810582e-23_wp ! zpress/ztemp / Rdry / 1e3_dp * 1e-6_dp ! alternatively use gridcell mass / grid volume
      call geomag(igrf_now, B_0, D, L, EMF, X)
    end if
  end if

  ! 2.)  Individual contributions
  if( lqradon .or. lgcr ) then ! don't even enter the loop if neither calculated
    do jk=1, nlev
      do jp=1, kproma
        if(lqradon) then 
           do i=1, size(idx_Rndecay)
             Rn_decay(i) = (pxtte(jp,jk,idx_Rndecay(i))*tmst +  &
                         & pxtm1(jp,jk,idx_Rndecay(i)))*m(jp,jk)*lam(i) !* conversion factor cm-3 to Bq k_dec = -tau_1/2 * ln(0.5)
           end do
           where(Rn_decay < 0.0_wp) Rn_decay = 0.0_wp
           radon_ipr(_RI_XYZ__(jp,jrow,jk)) = decay_ipr(Rn_decay, Rn_chain_ions) 
        end if
        if(lgcr) then
          pc_cha(jp, jrow) = CALC_PC(B_0, D, igrf_now, EMF, X,(pi/2.0_wp)-radlat(jp),radlon(jp))
          gcr_ipr(_RI_XYZ__(jp,jrow,jk)) = mass2vol(jp,jk) * usoskin_CRII(hgt(jp,jk), phi_now(1), pc_cha(jp, jrow))
        end if
      end do
    end do
  end if

  ! 3.) now the total ion pair production rate
  if(ltotalipr) &
     & total_ipr(_RI_XYZ__(1:kproma,jrow,1:nlev)) = radon_ipr(_RI_XYZ__(1:kproma,jrow,1:nlev)) + gcr_ipr(_RI_XYZ__(1:kproma,jrow,1:nlev))

  ! 4.) recombination rate constant (Elemental function)
!  if(lssion .or. laero) &
    krec(_RI_XYZ__(1:kproma,jrow,1:nlev)) = recom_brasseur(ztemp,m)

  ! 5.) steady state ion concentration (Elemental function), neglecting losses from nucleation and condensation onto aerosols
  ! these are calculated in the nucleation module
  if(lssion) then !.and. .not. laero) then
    small_ions_neg(_RI_XYZ__(1:kproma,jrow,1:nlev)) = steady_state_ions(total_ipr(_RI_XYZ__(1:kproma,jrow,1:nlev)), krec(_RI_XYZ__(1:kproma,jrow,1:nlev)))
    small_ions_pos(_RI_XYZ__(1:kproma,jrow,1:nlev)) = small_ions_neg(_RI_XYZ__(1:kproma,jrow,1:nlev)) ! charge balance, potentially wrong
  end if

end subroutine ions_physc


subroutine ions_free_memory
! clean up
! not needed here
  implicit none

  nullify(dimaxe)
  nullify(dimlen)
  nullify(dimname)
  nullify(dimunit)

end subroutine ions_free_memory


!
!----------- LOCAL SUBROUTINES ---------------------------------------------
!

subroutine ions_read_nml_cpl(status, iou)

  ! read the coupling information from the namelist files
  use messy_main_tools, only: read_nml_open, read_nml_check, read_nml_close

  implicit none
  integer, intent(in) :: iou
  integer, intent(out) :: status

  ! LOCAL
  character(len=*), parameter :: substr='ions_read_nml_cpl'
  logical                     :: lex      ! file exists ?
  integer                     :: fstat    ! file status

  namelist /CPL/ igrf, phi, lssion, ltotalipr, & 
                & cpl_ipr_Rn,  cpl_ipr_gcr, cpl_ipr_total, &
                & cpl_ionc, cpl_aero_cs, cpl_aero_r, cpl_aerocon

  status = 1
         
  call read_nml_open(lex, substr, iou, 'CPL', modstr)
  if (.not.lex) return    ! <modstr>.nml does not exist

  read(iou, nml=CPL, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'CPL', modstr)
  if (fstat /= 0) return  ! error while reading namelist

  call read_nml_close(substr, iou, modstr)
  status = 0

  ! sanity
  call ions_sanity_nml_cpl(status, igrf, 'GEOMAGNETIC FIELD')
  call ions_sanity_nml_cpl(status, phi, 'GCR MODULATION')

  status = 0 ! NO ERROR


end subroutine ions_read_nml_cpl


subroutine ions_sanity_nml_cpl(status, ch, danam)

  use messy_main_blather_bi, only: error_bi, warning_bi, info_bi
  implicit none

  integer, intent(out) :: status ! 0= all clear, 1= bad 
  type(t_chaobj_cpl), intent(in) :: ch
  character(len=*), intent(in) :: danam ! data name
  ! local
  character(len=*), parameter :: substr='sanity_read_nml_cpl'

  status = 1
  ! sanity checks
  if(trim(ch%cha) == '') then
     call warning_bi('empty channel name for'//danam//';', substr)
     status = 1
     return
  else
     call info_bi(danam//' channel :'//ch%cha)
     if (trim(ch%obj) == '') THEN
       call warning_bi('ERROR: empty channel object name for'//danam//';', substr)
       status = 1
       return
     else
       call info_bi(danam//' object  :'//ch%obj)
     end if
  end if

end subroutine ions_sanity_nml_cpl


end module messy_ions_si

