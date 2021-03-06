module messy_clamscirrus_global

  USE messy_clams_global, ONLY: rank, ntasks, DP, filenamelen

  character :: name*80

  character(filenamelen) :: trajfile, multiple_iwcfile
  character (len=10) :: input_filename ='cirrus.inp'

  real(DP) :: timestep_init

  real(DP),allocatable :: t_rad(:),densice_fit(:)

  real(DP),allocatable :: t_massendichte(:),t_masse(:),t_vol(:)
  real(DP),allocatable :: mean_free_path(:),slip_correction(:),terminal_velocity(:)
  real(DP),allocatable :: sed_laenge(:)

  !  logical, allocatable :: mdi(:)

  logical :: with_chem, use_traj

  logical :: troposphere = .false.   ! if true, Cirrus is restricted to the troposphere

  logical, parameter :: use_sh = .FALSE. ! in Messy mixing ratio of water vapor
  !  is used, not specific humidity

  integer :: ntimes_traj, ntimes_chem

  real(DP) :: pi

  !  real(DP), parameter :: delta_t=21600. ! [s]
  !  real(DP), parameter :: delta_t=86400. ! [s]
  !  real(DP), parameter :: characteristic_height = 650. ! [m]
  real(DP) :: characteristic_height, characteristic_height_100, sat_crit
  integer  :: type_ice_fit
  integer  :: freeze_out_type = 0 

  REAL(DP), PARAMETER   ::  mdi=-1.0000000E+30  ! missing data indicator
  REAL(DP), PARAMETER   ::  eps=1.E-6   !used by comparisons of real(DP) values

  ! Molekuelmassen

  real(DP), parameter :: masse_luft = 28.9644e-3 ! [kg/mol]
  real(DP), parameter :: masse_h2o = 18.0153e-3 ! [kg/mol]

  real(DP), parameter :: rhoice = 917  ! [kg/m^3]
  real(DP), parameter :: densice= 1.e6 ! [1/m^3]

  ! Umrechnungsfaktoren + Konstanten

  real(DP), parameter :: ctoa = 7.336E+21 
  real(DP), parameter :: universal_gas_constant=8.3143 ![J/(K mol)]
  real(DP), parameter :: k_boltzmann = 1.380662e-23 ! [J/K]
  real(DP), parameter :: avogadro=6.022e23 ! [molec/mol]
  real(DP), parameter :: viscosity = 1.72e-5      ! [kg/(m s)]

  ! Temperatur-Offset

  real(DP) :: delta_temp

  !   Variablen fuer die Laufzeitmessung

  character :: time*5
  integer ::finish,start,clock(3)
  real(DP)    ::seconds
  integer :: elements(8)   ! Nimmt die Zeitwerte auf

  real(kind(0d0)) :: current_time(1)

  real(kind(0d0)) :: starttime_traj, endtime_traj, starttime_iwc, endtime_iwc

  integer :: timestep_cirrus ! Timestep for regular call of cirrus in minutes 
                             !(not after MIX)


end module messy_clamscirrus_global
