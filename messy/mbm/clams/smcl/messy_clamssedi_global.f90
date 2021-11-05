!********************************************************************************!
!
! MODULE messy_clamssedi_global
! -------------
! Global definitions for program sedi
!
!********************************************************************************!

MODULE messy_clamssedi_global

  ! USE messy_main_constants_mem, ONLY: DP
  USE messy_clams_global, only: dp, prec
  
  implicit none

  !---------------------------------------------------------------------------
  !  Type definitions
  !---------------------------------------------------------------------------

  type particle_type     ! tbox (old name)                                             
     real(prec), dimension(:), pointer :: lat 
     real(prec), dimension(:), pointer :: lon
     real(prec), dimension(:), pointer :: lev        ! theta
     real(prec), dimension(:), pointer :: pressure
     real(prec), dimension(:), pointer :: temperature
     real(prec), dimension(:), pointer :: radius
     real(prec), dimension(:), pointer :: density
     real(prec), dimension(:), pointer :: sedimentation
     real(prec), dimension(:), pointer :: tsv
     real(prec), dimension(:), pointer :: hno3
     real(prec), dimension(:), pointer :: h2o
     real(prec), dimension(:), pointer :: sice
     real(prec), dimension(:), pointer :: snat
     real(prec), dimension(:), pointer :: icebin
     real(prec), dimension(:), pointer :: natbin
     real(prec), dimension(:), pointer :: airparcel_density_change
     ! real(prec), dimension(:), pointer :: lev_dot
     real(prec), dimension(:), pointer :: class
     real(prec), dimension(:), pointer :: particle_id

     real(prec), dimension(:), pointer :: natbin_diff
     real(prec), dimension(:), pointer :: icebin_diff
     integer, dimension(:), pointer :: tr_ind_up
     integer, dimension(:), pointer :: tr_ind_down
     integer, dimension(:), pointer :: lev_up
     integer, dimension(:), pointer :: lev_down
     real(prec), dimension(:), pointer :: snatmax
     real(prec), dimension(:), pointer :: sicemax
     real(prec), dimension(:), pointer :: tmin
  end type particle_type

  type airparcel_type   ! cbox (old name)
     real(prec), dimension(:), pointer :: lat
     real(prec), dimension(:), pointer :: lon
     real(prec), dimension(:), pointer :: lev
     real(prec), dimension(:), pointer :: hno3
     real(prec), dimension(:), pointer :: h2o
     real(prec), dimension(:), pointer :: natbin
     real(prec), dimension(:), pointer :: icebin
     real(prec), dimension(:), pointer :: natbin_diff
     real(prec), dimension(:), pointer :: icebin_diff
     real(prec), dimension(:), pointer :: diff_h2o
     real(prec), dimension(:), pointer :: diff_hno3
     real(prec), dimension(:), pointer :: temp
     real(prec), dimension(:), pointer :: press
     real(prec), dimension(:), pointer :: aer_h2so4

     real(prec), dimension(:,:), pointer :: coor        ! cartesian coordinats
     integer, dimension(:), pointer :: ntriang          ! counter for neighboring triangles
     integer, dimension(:,:), pointer :: triang_ids     ! indices of neighboring trinagles
     integer, dimension(:), pointer :: ilev             ! level index
  end type airparcel_type

  type triangle_type    ! triangle_cell
     integer, dimension(:,:,:), pointer :: airparcel_indices
     integer, dimension(:,:), pointer   :: nparticles_nat    
     integer, dimension(:,:), pointer   :: nparticles_ice    
  end type triangle_type

  !---------------------------------------------------------------------------
  ! Constants
  !---------------------------------------------------------------------------
  real(prec), parameter :: hno3_background = 20.e-9
  real(prec), parameter :: h2o_background = 5.e-6
  real(prec), parameter :: radius_default = .1e-6

  ! Max. Anzahl angrenzender triangle (Nachbarn)
  integer, parameter :: max_nb = 60

  ! Molekülmassen
  real, parameter :: masse_nat = 117.e-3 ! [kg/mol]
  real, parameter :: masse_hno3 = 63.e-3 ! [kg/mol]

  ! Werte der Standardatmosphäre
  real, parameter :: p_0 = 1013.25 ![hPa]
  real, parameter :: T_0 = 273.15  ![K]

  ! Minimal anzunehmende Werte für Wasser und HNO3 in Mischungsverhältnis
  real, parameter :: hno3_minval = 1.e-11
  real, parameter :: h2o_minval = 1.e-8

  !---------------------------------------------------------------------------
  ! Variables
  !---------------------------------------------------------------------------

  ! CTRL namelist:
  integer            :: nfine
  real(prec)         :: densnat_val, factor
  real (prec)        :: lat_min, lat_max
  real (prec)        :: lev_start, lev_end
  integer            :: timestep_sedi 
  integer            :: nparticle_max ! Max. Anzahl particles; Dimension fuer Channel
  character (len=80) :: ice_nuc_table, nat_nuc_table
  character (len=80) :: part_init_file = ''

  !
  ! Variables used for trajectory calculation in sedi:
  !
  ! UDT_sedi,VDT_sedi,WDT_sedi: winds at time of last sedi call 
  ! leveldt_sedi:               THETA/ZETA/PRESS values on time of last sedi call (for model data)
  ! dlevdzdt_sedi:              DTHETA/DZ at time of last sedi call 
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: UDT_sedi  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: VDT_sedi  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: WDT_sedi  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: leveldt_sedi => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: dlevdzdt_sedi   => NULL()

!!!!! verschoben nach messy_clamssedi_traj:
!!$  real(DP), dimension(:,:,:),pointer :: u0, v0, w0, dlevdz0, level0           ! winds at time t
!!$  real(DP), dimension(:,:,:),pointer :: umid, vmid, wmid,dlevdzmid, levelmid  ! winds at time t+delta_t/2

  ! character(150)      :: out_filename


  integer :: nlevs       ! number of levels (in lev_grid)
  integer :: pos_nlevs   ! number of levels (in pos_lev_grid)
  real(prec) :: lev_min, lev_max, lat_down, lat_up, r_coarse, r_high

  real(prec), dimension(:),   pointer :: lev_grid   => NULL() ! (nlevs+2), vertical grid (level values)
  real(prec), dimension(:,:), pointer :: lev_window => NULL() ! (2,nlevs)

  ! Arrays used in pos_sedi, including only levels between lev_start and lev_end 
  real(prec), dimension(:),   pointer :: pos_lev_grid  => NULL() ! (pos_nlevs), vertical grid (level values)
  real(prec), dimension(:),   pointer :: pos_lev_delta => NULL() ! (pos_nlevs), array of the layer-thickness
  real(prec), dimension(:),   pointer :: pos_r_grid    => NULL() ! (pos_nlevs), horizontal resolution 


  ! Werden NAT-Rocks am Anfang oder im Verlauf der Simulation initialisiert?
  logical :: nat_init_ex=.false. !!!!,nat_add_ex=.true.
  
  ! (local) start and end index in the global particles arrays for each processor
  integer :: part_ind_start, part_ind_end

  integer :: dnparticles = 0         ! local number of particles 
  integer :: nparticles = 0          ! anzahl_luftmassen=0
  integer :: nparticles_old = 0      ! anzahl_luftmassen_old=0
  integer :: nparticles_added = 0

  integer :: part_id_max = 0  

  real(prec), pointer :: part_id_max_co  ! Channelobject 
  real(prec), pointer :: nparticles_co  ! Channelobject 

  integer, dimension(:), allocatable :: nairparcels50_lev_init ! nairparcels_lev  !(nlevs)

  integer, dimension(:), allocatable :: nairparcels_lev        !(nlevs)
  integer, dimension(:), allocatable :: nairparcels50_lev      !(nlevs)

  integer, dimension(:), allocatable :: ntriang_lev            !(nlevs)
  integer                            :: max_ntriang_lev        ! max number of triangles per level


  logical :: nhemi                ! northern or southern hemisphere

  ! temperature fluctuations [0=no; 1=deterministic; 2=random]
  integer :: flaggary

  ! Timestep calculation nat and ice in hours [1 or 24]
  integer :: nat_tstep, ice_tstep

  ! used for the calculation of ice and NAT nucleation
  integer, parameter :: ntbins = 31
  integer, parameter :: nbins = 1000
  real, allocatable  :: snat_table(:,:), xnnat_table(:)
  real, allocatable  :: sice_table(:,:), xnice_table(:)
  
  ! cnt neg hno3/h2o in airparcel after sedi
  integer :: cnt_neg_hno3 = 0
  integer :: cnt_neg_h2o = 0

  type(particle_type)  :: particles
  type(airparcel_type) :: airparcels
  type(triangle_type)  :: triangles

END MODULE messy_clamssedi_global
