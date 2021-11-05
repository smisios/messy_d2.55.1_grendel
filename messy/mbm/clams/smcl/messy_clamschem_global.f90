module messy_clamschem_global

  USE messy_clams_global,        ONLY: DP, prec, mdi, specnamelen, filenamelen
  USE messy_cmn_photol_mem,      ONLY: IP_MAX

  !---------------------------------------------------------------------------
  ! Constants
  !---------------------------------------------------------------------------

  real(PREC), parameter, public :: pi = 3.14159265358979324d0
  real(PREC), parameter, public :: todeg = 180.0d0/pi
  real(PREC), parameter, public :: torad = pi/180.0d0 
  real(PREC), parameter, public :: r = 287.0 
  real(PREC), parameter, public :: re = 6.371E3
  real(PREC), parameter, public :: g = 9.81
  real(PREC), parameter, public :: cp = 1004.0  
  real(PREC), parameter, public :: rovercp = r/cp
  real(PREC), parameter, public :: roverg = r/g
  real(PREC), parameter, public :: p0 = 1013.15
  real(PREC), parameter, public :: boltz = 1.38066E-23
  real(PREC), parameter, public :: rmmair = 28.9644
  real(PREC), parameter, public :: u = 1.66056E-27
  real(PREC), parameter, public :: airm = rmmair*u
  real(PREC), parameter, public :: avogad = 6.022E23

  ! mdi aus messy_clams_global nutzen:
  !real(prec), parameter         :: missing_value = -1e30
  real(prec), parameter         :: missing_value = mdi

  integer, parameter            :: maxpacks=1000
  integer, parameter            :: maxtraj=100000

  integer, parameter :: jpdim = 48   ! number of implemented photolysis rates in dissoc

  integer, parameter :: nhetspec = 9  ! heterogenous variables on SPECARR 
  integer, parameter :: nhetpar  = 11  ! heterogenous variables on channel HETPAR
!  integer, parameter :: nhetvar  = 22  ! all heterogenous variables

  integer, parameter :: twod_nlats=18
  integer, parameter :: twod_npress=34
  integer, parameter :: twod_nmonths=12
  integer, parameter :: nspecs=4

  character*8,dimension(nspecs),parameter::twodspecs=(/'OH      ','O(1D)   ','Cl      ','HO2     '/)

  !---------------------------------------------------------------------------
  !  Type definitions
  !---------------------------------------------------------------------------

  type vartype
     character(20) :: name
     real(prec), dimension(:), allocatable :: values
     real(prec), dimension(:), allocatable :: startvalues
     real(prec), dimension(:), allocatable :: endvalues
     real(prec), dimension(:), allocatable :: val0
     real(prec), dimension(:), allocatable :: val1
  end type vartype
  
  TYPE rc_type
     CHARACTER(len=30)    :: name  
     CHARACTER(len=50)    :: units
     CHARACTER(len=80)    :: longname
     CHARACTER(len=120)   :: descr
     REAL(PREC), DIMENSION(:), POINTER :: values
  END TYPE rc_type

  TYPE rate_type
     CHARACTER(len=9)                :: jname_messy
     CHARACTER(len=10)               :: reactant
     CHARACTER(len=10)               :: product1
     CHARACTER(len=10)               :: product2
     CHARACTER(len=50)               :: units
     CHARACTER(len=80)               :: longname
     REAL(DP), DIMENSION(:), POINTER :: values
  END TYPE rate_type

  TYPE hetvartype
     CHARACTER(len=SPECNAMELEN)  :: name  
     CHARACTER(len=80)           :: longname
     CHARACTER(len=50)           :: units
  END TYPE hetvartype

  !---------------------------------------------------------------------------
  ! Variables
  !---------------------------------------------------------------------------

  integer, dimension(maxpacks)  :: jpnl_count, jpnl_offset
  integer                       :: ipart=1
  integer                       :: npacks=1
  integer                       :: ntraj
  integer                       :: ntraj_rank

  INTEGER                       :: rank    ! number of local process
  INTEGER                       :: ntasks  ! number of desired processes
  INTEGER                       :: ierr    ! MPI error


  ! CTRL namelist:
  INTEGER                :: timestep_chem       ! timestep for calling chem (seconds)
  INTEGER                :: ncdt                ! internal chemistry timestep (seconds)
  logical                :: IODUMP = .FALSE.    !
  logical                :: IODUMPo = .FALSE.   !
  CHARACTER(filenamelen) :: dsn_twodavg         !
  logical                :: rates = .FALSE.     ! write out  rates
  logical                :: const = .FALSE.     ! write out rate constants
  logical                :: hetparam = .FALSE.  ! write out parameter for heterogenous chemistry
  logical                :: emrates = .FALSE.   !
  ! messy_clamschem_messy_clamschem_asad_mod.f90: nfphot
  ! messy_clamschem_messy_clamschem_asad_mod_clams: lhet, lphotol, chemdata_type


  !LOGICAL           :: chem_first_time
!!$  LOGICAL           ::  asad_gfirst ! op_pj_20170110 moved to messy_clams_global
  logical           :: dates30=.false.
  logical           :: emit=.false.

  integer, dimension(jpdim) :: ip_messy ! corresponding MESSy photolysis indices from messy_cmn_photol_mem

  TYPE(rate_type), DIMENSION(:), POINTER :: DISSOC_RATE => NULL()

  real(PREC),dimension(twod_nlats,twod_npress,twod_nmonths,nspecs),save::twod_avg_all
  real(PREC),dimension(twod_nlats,twod_npress,twod_nmonths),       save::twod_avg
  real(PREC),dimension(twod_nlats),                                save::twod_lats
  real(PREC),dimension(twod_npress),                               save::twod_press,twod_lnp
  logical,                                                         save::twod_data_read=.false.

  character(filenamelen) :: twodavgfile


  real(kind=DP),dimension(:,:),allocatable :: ftr

  real(kind=DP),dimension(:,:),allocatable :: ftr_ini

  logical,      dimension(:),  allocatable :: missing_index

  integer,       dimension(:), allocatable :: ftrindex
  integer                                  :: nfamily
  character(specnamelen), dimension(:), allocatable :: fnames

 

  integer,      dimension(:),  allocatable :: therm_flag


  !  cindex is an index array of the length of the 
  !  number of photo reactions in ASAD (namely jppj) 
  !  the indices point to the number of the photol reaction in 
  !  the provided photlysis scheme
  integer,      dimension(:),  allocatable :: cindex
  integer                                  :: iphot


  ! used to store the zenith angle (used in `photol')
  real*8,       dimension(:),  allocatable :: zangle


  integer,      dimension(:),  allocatable :: slenb, slent, slenj, slenh


  type(hetvartype), dimension(nhetspec+nhetpar) :: hetvar


  !     parth is the fraction of species in gas phase (HCl,HNO3,H2O),
  !     wt weight fraction of soluable gases (S,N,Cl)
  !     ar area of 4 different surfaces(ICE,NAT,SAT,liquid)
  !     aerosol number density
  real(prec), dimension(:),     allocatable :: densaero, aerh2so4
  real(prec), dimension(:,:),   allocatable :: con, wt, ar
  real(prec), dimension(:,:,:), allocatable :: parth
  integer,    dimension(:),     allocatable :: shindex  
  integer,    dimension(:,:),   allocatable :: chindex  
  logical,    dimension(:),     allocatable :: sedinucl


  ! ----------------------------------
  ! used to store the trajectory position data time, latitude, longitude, 
  ! pressure, temperature and zenith angle etc for communication between 
  ! these routines within clamschem 
  ! dynamic, emissn, chem, get_const_2d, mixinit
  ! ---------------------------------
  real(prec),dimension(:),allocatable::  lats,lons,temps,press,sza,theta,slt
  
  real(kind(0.d0))             :: js
  real(prec)                   :: js_monthmean
  integer,dimension(7)         :: date_time



end module messy_clamschem_global
