!********************************************************************************!
!
! MODULE messy_clams_global
! -------------------------
!
! Module contains global declaration for clams packages
!
!********************************************************************************!

MODULE messy_clams_global

  !---------------------------------------------------------------------------
  ! Constants
  !---------------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: pi, radius_earth

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: prec = dp

  REAL(PREC), PARAMETER   ::  mdi=-1.00D+30  ! missing data indicator
  REAL(PREC), PARAMETER   ::  eps=1E-6   !used by comparisons of real values


  ! real, parameter :: universal_gas_constant = 8.3143 ! [J/(K mol)]
  real, parameter :: k_boltzmann = 1.380662e-23 ! [J/K]
  real, parameter :: avogadro = 6.022e23 ! [molec/mol]
  real, parameter :: viscosity = 1.72e-5      ! [kg/(m s)]
  
  ! Molekülmassen
  real, parameter :: masse_h2o = 18.e-3  ! [kg/mol]
  real, parameter :: masse_luft = 28.97e-3 ! [kg/mol] 

  real, parameter :: rho_ice = 0.917*1.e3  ! [kg/m^3]
  
  INTEGER, PARAMETER :: SPECNAMELEN = 15

  INTEGER, PARAMETER :: filenamelen = 250

  INTEGER, PARAMETER :: max_params = 100  ! max. number of output parameters

  INTEGER, PARAMETER :: nmaxcltr = 500  ! max. number of clams messy-tracers

  INTEGER, PARAMETER :: maxspec = 500 ! max. number of species (in chem + mix / tracer)

  REAL(DP), PARAMETER :: sample_interval = 0.

  ! op_sb_20191022
  REAL(dp), PARAMETER :: cmcell = 1._dp
  !---------------------------------------------------------------------------
  !  Type definitions
  !---------------------------------------------------------------------------

  TYPE species_type
     CHARACTER(len=SPECNAMELEN)        :: name  
     CHARACTER(len=2)                  :: ctype
     CHARACTER(len=50)                 :: units
     CHARACTER(len=80)                 :: longname
     REAL(PREC), DIMENSION(:), POINTER :: values
  END TYPE species_type

  TYPE species_type_3d
     CHARACTER(len=SPECNAMELEN)        :: name  
     CHARACTER(len=2)                  :: ctype
     CHARACTER(len=50)                 :: units
     CHARACTER(len=80)                 :: longname
     REAL(PREC), DIMENSION(:,:,:), POINTER :: values
  END TYPE species_type_3d

  TYPE :: timetype
     INTEGER :: year, month, day, sec
  END TYPE timetype

  TYPE datatype
     CHARACTER(len=80)                     :: name  
     CHARACTER(len=80)                     :: units
     CHARACTER(len=80)                     :: longname
     REAL(PREC), DIMENSION(:,:,:), POINTER :: values
  END TYPE datatype

  TYPE paramtype
     CHARACTER(len=80)                 :: name  
     CHARACTER(len=80)                 :: units
     CHARACTER(len=80)                 :: longname
     REAL(PREC), DIMENSION(:), POINTER :: values
  END TYPE paramtype

  ! ju_ec_20181108+
  ! One could add in other things like the tracer indexes for the messy tracer
  ! set or for the clams specarr (the clams parcel data storage).
  TYPE cl_grid_tracer
      CHARACTER(len=SPECNAMELEN) :: name     ! tracer names
      LOGICAL :: couple             ! Whether tracer is coupled. Specifically,
                                    ! that means whether the tracer values
                                    ! for clams PARCELS should be set to the
                                    ! values for their EMAC cells.
      INTEGER :: lower              ! Lower boundary for coupling.
      INTEGER :: upper              ! Upper boundary for coupling.
      CHARACTER(len=20) :: fill_ch  ! channel for fill values
      CHARACTER(len=20) :: fill_ob  ! object for fill values
      ! For the fill values, if the fill channel is "#NAN" then the fill_ob
      ! will be ignored and the fill value will be -1.
  END TYPE cl_grid_tracer
  ! ju_ec_20181108-

  !---------------------------------------------------------------------------
  ! Variables
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  !     CTRL namelist:
  !---------------------------------------------------------------------------

  ! first_initfile:    filename of first initial positions file
  ! initfile:          filename of initial positions file
  ! met_dir:           directory with meteorological files
  ! met_prefix:        prefix of met. file (isen_xxx, ecmwf, ...)
  ! theta_dir:         directory with meteorological files on theta levels
  ! theta_prefix:      prefix of met. file (isen_xxx, ...) on theta levels
  ! username:          user name
  ! levdotname:        short name of variable containing THETA-DOT information

  CHARACTER(filenamelen), PUBLIC :: first_initfile='' 
  CHARACTER(filenamelen), PUBLIC :: initfile=''      
  CHARACTER(filenamelen), PUBLIC :: met_dir=''       
  CHARACTER(30),          PUBLIC :: met_prefix=''  
  CHARACTER(filenamelen), PUBLIC :: theta_dir=''   
  CHARACTER(30),          PUBLIC :: theta_prefix=''  
  CHARACTER(30),          PUBLIC :: username
  CHARACTER(20),          PUBLIC :: levdotname   

  ! rres:              
  ! rres_shuffle:      
  ! buffersize:        default buffersize = 1 MB
  ! met_freq:          frequency for met. files (in hours)
  ! timestep_clamsout: timestep for clams output
  ! lperpetuum:        switch for perpetuum runs (reset of BA (mean age) 
  !                               at the beginning of new year)
  ! ldiagout:          switch diagnostic output on/off
  ! use_3d:            use THETA-DOT Information of FILE or set theta-dot to Zero 
  ! resume_run:        resume run with new init file (use dnparts from init-file)
  ! nparams:           number of output parameters 

  REAL(DP),       PUBLIC :: rres = 0.0_dp
  REAL(DP),       PUBLIC :: rres_shuffle = 0.0_dp

  INTEGER,        PUBLIC :: buffersize = 1048576
  INTEGER,        PUBLIC :: met_freq     
  INTEGER,        PUBLIC :: timestep_clamsout = 24    

  LOGICAL,        PUBLIC :: lperpetuum = .false.     
  LOGICAL,        PUBLIC :: ldiagout   = .false.     
  LOGICAL,        PUBLIC :: use_3d     = .TRUE. 
  LOGICAL,        PUBLIC :: resume_run = .FALSE. 

  INTEGER,        PUBLIC :: nparams      

  ! init_h2o_emac:     switch, if CLaMS H2O is initialized from EMAC (T) or initfile (F)
  LOGICAL,        PUBLIC :: init_h2o_emac = .false.  

  ! CTRL namelist gridding variables
  ! ju_ec_20181108+
  LOGICAL,        PUBLIC :: clams_gridding=.FALSE.  ! Grid clams data to the echam grid
  LOGICAL,        PUBLIC :: clams_grid_verbose=.FALSE.  ! Verbosity of clams gridding
  INTEGER,        PUBLIC :: n_cltr = 0 ! number of clams messy-tracers 
  TYPE(cl_grid_tracer), DIMENSION(nmaxcltr), PUBLIC :: cl_grid_tracers
  ! ju_ec_20181108-

  ! paramnames can not be allocatable, because it is read from a namelist !
  CHARACTER(len=SPECNAMELEN), PUBLIC :: paramnames(max_params)

 
  !---------------------------------------------------------------------------

  INTEGER,        PUBLIC :: nparams_old

  CHARACTER(len=SPECNAMELEN), PUBLIC :: paramnames_old(max_params)
  CHARACTER(len=SPECNAMELEN), DIMENSION(:), POINTER :: specnames 


  ! Channel objects fuer veraenderliche Namelist-Variablen:
  REAL(PREC), DIMENSION(:), POINTER :: dnparts_co
  REAL(PREC), DIMENSION(:), POINTER :: dnparts_max_co
  REAL(PREC),               POINTER :: grid_switch_co

  ! Channel objects fuer naechste Datenzeit (fuer Restarts speichern)
  REAL(PREC),               POINTER  :: pre_year_co, pre_month_co, &
                                        pre_day_co, pre_sec_co


  ! pre_metfile, fut_metfile :  names of the two meteorological files currently used 
  !                              (data time before and data time after current time)
  ! fut_year, fut_month, fut_day, fut_sec: time of next met. file (fut_metfile)
  ! pre_year, pre_month, pre_day, pre_sec: time of previous met. file (pre_metfile) 
  CHARACTER(filenamelen),PUBLIC  :: pre_metfile, fut_metfile
  CHARACTER(filenamelen),PUBLIC  :: pre_thetafile, fut_thetafile
  INTEGER,               PUBLIC  :: fut_year, fut_month, fut_day, fut_sec 
  INTEGER                        :: pre_year, pre_month, pre_day, pre_sec


  ! UDT,VDT,WDT:    winds on current time (used in TRAJ)
  ! UFUT,VFUT,WFUT: winds on future data time (used in TRAJ and SEDI)
  ! leveldt:        THETA/ZETA/PRESS values on current time 
  !                 (for model data, used in TRAJ)
  ! levelfut:       THETA/ZETA/PRESS values on future data time 
  !                 (for model data, use in TRAJ and SEDI)
  ! dlevdzdt:       delta_level/delta_z on current time (dummy variable in TRAJ)
  ! dlevdzfut:      delta_level/delta_z on future data time (used in SEDI)

  ! Used in TRAJ:
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: UDT       => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: VDT       => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: WDT       => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: leveldt   => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: dlevdzdt  => NULL()

  ! Used in TRAJ and SEDI:
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: UFUT      => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: VFUT      => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: WFUT      => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: levelfut  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: dlevdzfut => NULL()

 
  ! for trajectory calculation (used in TRAJ and SEDI)
  INTEGER, POINTER :: levelno(:)

  ! corr_thetadot: true, if correction file for thetadot is used
  ! corr_3d: true, if correction file has three dimensions (time,lat and level)
  !          otherwise correction file has two dimension (time and level)
  ! corrfile : filename of correction file (for dzeta/dt)
  LOGICAL :: corr_thetadot = .false.     !->CTRL-Namelist
  LOGICAL :: corr_3d = .false.
  CHARACTER(filenamelen) :: corrfile=''

  ! time_corr: containing data of variable "time" in correction file
  ! theta_corr: containing data of variable theta|zeta|press in correction file
  ! lat_corr: containing data of variable "lat" in correction file
  ! thetadot_corr: containing data of variable "CORR" in correction file (if it is two-dimensional)
  ! thetadot_corr3d: containing data of variable "CORR" in correction file (if it is three-dimensional)
  real(dp),        dimension(:),     pointer :: time_corr        => NULL()
  real(prec),      dimension(:),     pointer :: theta_corr       => NULL()
  real(prec),      dimension(:),     pointer :: lat_corr         => NULL()
  real(prec),      dimension(:,:),   pointer :: thetadot_corr    => NULL()
  real(prec),      dimension(:,:,:), pointer :: thetadot_corr3d  => NULL()
  


  ! Workspace:
  INTEGER, PUBLIC  :: nparts = 0
  INTEGER, PUBLIC  :: nparts_max  = 0
  INTEGER, PUBLIC  :: dnparts_max = 0
  INTEGER, PUBLIC  :: dnparts     = 0
  INTEGER, PUBLIC  :: dnparts_max_shuffle = 0
  INTEGER, PUBLIC  :: nparts_max_shuffle = 0

  !jug DISSOC dimenstions
  INTEGER, PUBLIC :: jpslev = 0
  INTEGER, PUBLIC :: jpschi = 0
  INTEGER, PUBLIC :: jpwave = 0
  INTEGER, PUBLIC :: jplats = 0


  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: idx => NULL()
   
 
  TYPE(datatype), DIMENSION(:), POINTER, PUBLIC :: PREDATA => NULL()
  TYPE(datatype), DIMENSION(:), POINTER, PUBLIC :: FUTDATA => NULL()


!!!!! Deklarationen nach messy_clams_si verschoben:
!!$  REAL(PREC), DIMENSION(:), POINTER :: LAT        => NULL()
!!$  REAL(PREC), DIMENSION(:), POINTER :: LON        => NULL()
!!$  REAL(PREC), DIMENSION(:), POINTER :: LEV        => NULL()
!!$  REAL(DP),   DIMENSION(:), POINTER :: JULSEC     => NULL()
!!$
!!$  REAL(PREC), DIMENSION(:), POINTER :: LAT_OLD        => NULL()
!!$  REAL(PREC), DIMENSION(:), POINTER :: LON_OLD        => NULL()
!!$  REAL(PREC), DIMENSION(:), POINTER :: LEV_OLD        => NULL()
!!$
!!$  REAL(PREC), DIMENSION(:), POINTER :: LAT_OLD_MIX    => NULL()
!!$  REAL(PREC), DIMENSION(:), POINTER :: LON_OLD_MIX    => NULL()
!!$
!!$  REAL(DP),   POINTER               :: JULTIME


  
  REAL(DP), DIMENSION(:), POINTER, PUBLIC :: E5_LAT    => NULL()
  REAL(DP), DIMENSION(:), POINTER, PUBLIC :: E5_LON    => NULL()
  REAL(DP), DIMENSION(:), POINTER, PUBLIC :: E5_LEVEL  => NULL()

  

  TYPE(paramtype),    DIMENSION(:), POINTER :: PARAM      => NULL()

  TYPE(paramtype),    DIMENSION(:), POINTER :: PARAM_OLD  => NULL()

  TYPE(SPECIES_TYPE), DIMENSION(:), POINTER :: SPECARR    => NULL()

  ! nx, ny, nz: grid size in meteorological dataset (lon, lat, lev)
  ! longrid:    longitude grid in met. file (nx)
  ! latgrid:    latitude grid in met. file (ny)
  ! levelgrid:  levels (in met. file) (nz)
  ! longrid_rad: longitude grid in met. file in radians (nx+1)
  ! latgrid_rad: latitude grid in met. file in radians (0:ny+1)
  INTEGER                             :: nx, ny, nz
  INTEGER                             :: ntheta
  REAL(PREC),dimension(:),    pointer :: longrid      => NULL()
  REAL(PREC),dimension(:),    pointer :: latgrid      => NULL()
  REAL(PREC),dimension(:),    pointer :: longrid_rad  => NULL()
  REAL(PREC),dimension(:),    pointer :: latgrid_rad  => NULL()
  REAL(PREC),dimension(:),    pointer :: levelgrid    => NULL()
  REAL(PREC),dimension(:),    pointer :: thetagrid    => NULL()

!!$ use_modellev: true, if model levels are used  
  ! level_is_vertcoor: true, if level of trajectories (vertical coordinate in pos-files)
  !                    is also vertical coordinate in windfiles
  ! asc_level:    true, if coordinate variable in pos-files is in ascending order 
  ! asc_lat:      true, if latitudes in windfile are in ascending order
  ! loglev:       interpolate level log. linear (true) or linear (false)
  ! logpress:     interpolate PRESS linear + log. linear (true) or only linear (false)
  !logical :: use_modellev = .false.
  logical :: level_is_vertcoor = .true.
  logical :: asc_level    = .true.
  logical :: asc_lat      = .false.
  logical :: loglev       = .false.
  logical :: logpress     = .true.
              
  
  INTEGER :: NSPEC=0, NCHEMSPEC=0

  ! ntasks : number of desired processes
  ! rank   : number of local process
  INTEGER :: ntasks=1
  INTEGER :: rank=0

  ! vertical coordinate in initfile 
  CHARACTER(30) :: init_vertcoorname=''

  ! vertical coordinate in wind files
  CHARACTER(30) :: met_vertcoorname=''

  LOGICAL :: dates30 = .false.

  ! irdday  : interval between windfile data in days
  ! irdsec  : interval between windfile data in seconds
  INTEGER :: irdday, irdsec

  ! Switches for active/inactive event_state of submodels
  LOGICAL :: lmixevent, lbmixevent, lcirrusevent, lchemevent, lclamsoutevent, &
             lsedievent, ltracerevent, ldeepconvevent

  INTEGER :: H2O_index, H2O_100_index, IWC_index, IWC_100_index, CLWC_index
  
  ! coupled to ECHAM5 ?
  LOGICAL :: lcoupled = .false.

  ! op_pj_20170110+
  LOGICAL           ::  asad_gfirst ! used in CLAMSCHEM and CLAMSMIX
  ! op_pj_20170110-

  ! TIME SETTING
  INTEGER :: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
  INTEGER :: YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT
  INTEGER :: MINUTE_NEXT, SECOND_NEXT
  INTEGER :: YEAR_START, MONTH_START, DAY_START, HOUR_START
  INTEGER :: MINUTE_START, SECOND_START
  REAL(PREC) :: delta_time
  LOGICAL    :: lstart, lfirst_cycle, lresume


END MODULE messy_clams_global

