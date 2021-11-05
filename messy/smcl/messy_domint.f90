!> mz_sg_20200730 | f90

!! DOMINT: a sub-module for calculating domain integrals
!! CL level implementation
!!
!! Author: S.Gromov (MPIC), 2018-2020

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

module messy_domint

  use messy_main_constants_mem, only: dp, STRLEN_LONG, STRLEN_MEDIUM

  implicit none
  private

  intrinsic :: null

  public :: dp

  public :: domint_allocate_comp
  public :: domint_allocate_task
  public :: domint_integrate
  public :: domint_deallocate
  public :: str_field
  public :: num2str

! module global parameters
  character(len=*), parameter, public :: modstr = 'domint'
  character(len=*), parameter, public :: modver = '1.1'

! string lengths a-la CHANNEL
  integer, parameter, public :: STRLEN_A_LA_CHANNEL = 2*STRLEN_LONG + 1
  integer, parameter, public :: STRLEN_A_LA_OBJECT  = 2*STRLEN_MEDIUM + 1 + 4


! domain composition kind
  type t_domint_comp
  ! info, etc.
    integer           :: ddims(1:3) = 0              ! up to 3 domain dimensions (e.g. (/nvdom,nzdom,1/) for atmosphere)
    integer           :: gdims(1:3) = 0              ! up to 3D grid dimensions (e.g. (/nproma,nlev,nrow/))
    integer           :: ndom = 0                    ! total no. of domains
    logical           :: linvinit = .false.          ! flag for invariable domains initialised
  ! fields
  ! - (proma,lev,row,dom_seqno)
    logical, pointer  :: mask(:,:,:,:) => null()     ! domain mask (.true. when gridpoint is in the given domain)
  ! - (x,x,x)
    integer, pointer  :: seqno(:,:,:) => null()      ! sequential decomposition of domain indices, up to 3 dimensions of domains
  ! - (dom_seqno)
    real(dp), pointer :: cells(:) => null()          ! domain cells count
  end type t_domint_comp
  public :: t_domint_comp

! domint input info type
  type t_domint_input
    character(len=STRLEN_A_LA_CHANNEL+1+STRLEN_A_LA_OBJECT) &
                      :: ref = ''                    ! 'channel:object' reference
    logical           :: is = .false.                ! flag for presence
    string            :: capt, unit                  ! caption & unit, if available
    real(dp), pointer :: p1(:)     => null()         ! rank=1 field
    real(dp), pointer :: p2(:,:)   => null()         ! rank=2 field
    real(dp), pointer :: p3(:,:,:) => null()         ! rank=3 field
  end type t_domint_input
  public t_domint_input

! task operation data array
  type t_domint_task
  ! info, flags, etc.
    logical           :: lwgt = .false.              ! .true. if we do weighting
    string            :: name, capt, unit, info      ! info fields
    real(dp)          :: sf = 1._dp                  ! scaling factor
  ! fields
  ! - (e.g. proma,lev,row)
    type(t_domint_input) :: var, wgt                    ! input variable & weighting
  ! - (dom_seqno)
    real(dp), pointer :: int_dom(:) => null()        ! domain integrals
    real(dp), pointer :: int_mem(:,:,:,:) => null()  ! for channel output
  end type t_domint_task
  public :: t_domint_task

contains

! ----------------------------------------------------------------------
! add domain composition, allocate memory for it
  integer function domint_allocate_comp(comp, ddims, gdims)

    implicit none
    intrinsic :: size

    type(t_domint_comp), &
          intent(inout) :: comp
    integer, intent(in) :: ddims(:)  ! dimension sizes for domains, up to 3 dims (e.g. 2 (/nvdom,nzdom/) for atmosphere)
    integer, intent(in) :: gdims(:)  ! grid decomposition (e.g. nproma, nlev, ngpblks=nrow )

    integer :: n, s1, s2, s3
    integer :: addims(1:3), agdims(1:3)

  ! check input dimensions
    if ((size(ddims).lt.1).or.(size(ddims).gt.3)) then
      print *,'/// domain dimensions input is wrong (should be 1 to 3 dimensions)'
      domint_allocate_comp = -2
      return
    endif
    if ((size(gdims).lt.1).or.(size(gdims).gt.3)) then
      print *,'/// grid dimensions input is wrong (should be 1 to 3 dimensions)'
      domint_allocate_comp = -3
      return
    endif

  ! transfer dimensions and size
    addims(:) = 1
    do n=1, size(ddims)
      addims(n) = max(ddims(n),0)
    enddo
    agdims(:) = 1
    do n=1, size(gdims)
      agdims(n) = max(gdims(n),0)
    enddo

  ! check total dimensions product
    if (product(addims).lt.1) then
      print *,'/// some or all domain dimensions are not positive'
      domint_allocate_comp = -4
      return
    endif
    comp%ddims(:) = addims(:)
    comp%ndom = product(addims)
    if (product(agdims).lt.1) then
      print *,'/// some or all grid dimensions are not positive'
      domint_allocate_comp = -5
      return
    endif
    comp%gdims(:) = agdims(:)

  ! linear decomposition of domain indices
    allocate(comp%seqno(addims(1),addims(2),addims(3)))
    do s1=1, addims(1)
      do s2=1, addims(2)
        do s3=1, addims(3)
          comp%seqno(s1,s2,s3) = s3+addims(3)*((s2-1)+(s1-1)*addims(2))  ! s1 varies fastest
        enddo
      enddo
    enddo

  ! domain mask & cells count
    allocate(comp%mask(agdims(1),agdims(2),agdims(3),comp%ndom))  ! TODO: so far 3D fields only
    allocate(comp%cells(comp%ndom))

    domint_allocate_comp = 0 ! all OK

  end function domint_allocate_comp

! ----------------------------------------------------------------------

! add/allocate memory for int. task
  integer function domint_allocate_task(comp,tasks)

    implicit none
    intrinsic :: size, associated

    type(t_domint_task), intent(inout), pointer :: tasks(:)
    type(t_domint_comp), intent(in)             :: comp
    type(t_domint_task), pointer                :: xtmp(:) => null()
    integer                                     :: xn

  ! adding new entry to task
    xn = 0
    if (associated(tasks)) then
      xn = size(tasks)
      allocate(xtmp(1:xn))
      xtmp(:) = tasks(:)
      deallocate(tasks)
    endif
    allocate(tasks(xn+1))
    if (associated(xtmp)) then
      tasks(1:xn) = xtmp(1:xn)
      deallocate(xtmp)
    endif

  ! allocate memory for the integrals
    allocate(tasks(xn+1)%int_dom(comp%ndom))
  ! create mem pointer for channel output
    tasks(xn+1)%int_mem(1:comp%ddims(1),1:comp%ddims(2),1:comp%ddims(3),1:1) => tasks(xn+1)%int_dom(1:comp%ndom)

  ! all OK, return task no.
    domint_allocate_task = size(tasks)

  end function domint_allocate_task

! ----------------------------------------------------------------------

  subroutine domint_integrate(status, comp, tasks, ntask)

    implicit none

    integer, intent(out)            :: status
    type(t_domint_comp), intent(in) :: comp      ! composition
    type(t_domint_task), intent(in) :: tasks(:)  ! tasks to do
    integer, intent(in), optional   :: ntask     ! (optional) no. of task

    character(len=*), parameter     :: substr = 'domint_integrate'
    integer                         :: jd, jt

    intrinsic :: merge, sum, size

    domain_loop: do jd = 1, comp%ndom

    ! calculate # of cells
      comp%cells(jd) = sum( merge(1., 0., comp%mask(:,:,:,jd)) )

    ! do tasks
      task_loop: do jt = 1, size(tasks)

      ! if only a certain task is requested
        if ( present(ntask) ) then
          if ( jt.ne.ntask ) cycle
        endif

      ! integrate
      ! 2D
        if ( associated(tasks(jt)%var%p2) ) then
          if ( tasks(jt)%lwgt ) then
          ! weighted integral
            tasks(jt)%int_dom(jd) = &
              sum( merge(tasks(jt)%var%p2(:,:), 0.0_dp, comp%mask(:,:,1,jd)) * tasks(jt)%wgt%p2(:,:) )
          else
          ! plain integral
            tasks(jt)%int_dom(jd) = &
              sum( merge(tasks(jt)%var%p2(:,:), 0.0_dp, comp%mask(:,:,1,jd)) )
          endif
        endif
      ! 3D
        if ( associated(tasks(jt)%var%p3) ) then
          if ( tasks(jt)%lwgt ) then
          ! weighted integral
            tasks(jt)%int_dom(jd) = &
              sum( merge(tasks(jt)%var%p3(:,:,:), 0.0_dp, comp%mask(:,:,:,jd)) * tasks(jt)%wgt%p3(:,:,:) )
          else
          ! plain integral
            tasks(jt)%int_dom(jd) = &
              sum( merge(tasks(jt)%var%p3(:,:,:), 0.0_dp, comp%mask(:,:,:,jd)) )
          endif
        endif

      ! scaling factor
        tasks(jt)%int_dom(jd) = tasks(jt)%int_dom(jd) * tasks(jt)%sf

      enddo task_loop

    enddo domain_loop

  ! done
    status = 0

  end subroutine domint_integrate

! ----------------------------------------------------------------------

  subroutine domint_deallocate(comp,tasks)

    implicit none

    type(t_domint_comp), intent(inout) :: comp
    type(t_domint_task), intent(inout), pointer :: tasks(:)

    integer :: n

  ! tasks
    do n = 1, size(tasks)
      deallocate(tasks(n)%int_dom)
    enddo
    deallocate(tasks)

  ! composition
    deallocate(comp%cells)
    deallocate(comp%mask)
    deallocate(comp%seqno)

  end subroutine domint_deallocate

! ----------------------------------------------------------------------

! extract nth field from a string using given separator, trimmed by default
  function str_field(str, fsep, fnum, notrim)
    implicit none
    intrinsic trim
    string                        :: str_field
    character(len=*),  intent(in) :: str     !> input string
    character,         intent(in) :: fsep    !> field separator
    integer,           intent(in) :: fnum    !> requested field no.
    logical, optional, intent(in) :: notrim  !> set to .true. to avoid trimming by default
    integer                       :: cn, fn
    string                        :: istr    !> working string
    logical                       :: inotrim
    istr = str//fsep  ! add separator in case a string without separators is given
    inotrim = .false.; if (present(notrim)) inotrim = notrim
    str_field = ''; fn = 0; cn = 0
    do while (cn.lt.len(istr))
      cn = cn + 1
      if ( istr(cn:cn).eq.fsep ) then
        fn = fn + 1
        if ( fn.eq.fnum ) then
          if ( .not.inotrim ) str_field = trim(adjustl(str_field))
          return
        endif
        str_field = ''
      else
        str_field = str_field//istr(cn:cn)
      endif
    enddo
  end function str_field

! ----------------------------------------------------------------------

! number-to-string conversion
  function num2str(i,r,l,fmt)
    integer,      intent(in), optional :: i
    real(dp),     intent(in), optional :: r
    logical,      intent(in), optional :: l
    character(*), intent(in), optional :: fmt
    character(len=256)                 :: afmt, astr
    logical                            :: fnp
    string                             :: num2str         ! 'dynamic' strings version
!!  character(len=32)                  :: num2str_        ! 'static' version (req. macros with trim())
    fnp = .not.(present(fmt)); if (.not.fnp) afmt = trim(fmt)
    if ( present(i)  ) then ; if (fnp) afmt = '(I0)';   write(astr,trim(afmt)) i;  endif
    if ( present(r)  ) then ; if (fnp) afmt = '(F0.1)'; write(astr,trim(afmt)) r;  endif
    if ( present(l)  ) then ; if (fnp) afmt = '(L1)';   write(astr,trim(afmt)) l;  endif
    num2str = trim(astr)
  end function num2str

! ----------------------------------------------------------------------

end module messy_domint
