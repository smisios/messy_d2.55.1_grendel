Module messy_clamsmix_global

  USE messy_clams_global, only: sp, dp, prec, specnamelen, maxspec

  integer, parameter :: max_nb=15
  !  integer, parameter :: strlen = 10

  type,public :: ap                        ! air parcel (ap)
     real(prec)               :: lon       ! its longitude (UKMO-notation)
     real(prec)               :: lat       ! its latitude      (-''-)
     real(prec)               :: lev       ! its pot. temp.
     real(dp)                 :: time_init ! its julian sec
     real(prec), dimension(3) :: coor      ! kartesian coordinates of ap on a unit sphere
     real(prec), dimension(3) :: coor_old  ! old kart. coordinates for anisentropic mixing
     real                     :: theta
     real                     :: theta_old
     real                     :: bvf
     real                     :: bvf_old
     !real(prec)               :: temp      ! its temperature
     !real(prec)               :: press     ! its pressure
     real(prec)               :: tracer    ! its tracer mixing ratio
     integer                  :: state     ! mixing state: 1-pure elimination, 
                                           !               2-pure adaptation, 
                                           !               3-adaptation induced elimination
     integer                  :: state_vert! vertical mixing state: 1-resolved mixing
                                           !                        2-unreloved mixing      
     logical                  :: subset    ! 
     real(prec), dimension(:),pointer   :: c  ! mixing ratios of chemical species
     integer                  :: nb        ! # of Voronoi neighbours
     integer,dimension(:),pointer :: ind   ! their indicies
  end type ap


  ! short version of ap
  type,public :: ap_short                  ! air parcel (ap)
     real(prec)       :: lon       ! its longitude (UKMO-notation)
     real(prec)       :: lat       ! its latitude      (-''-)
     real(prec)       :: lev       ! its pot. temp.
     real(dp)         :: time_init ! its julian sec
     real(prec)       :: temp      ! its temperature
     real(prec)       :: press     ! its pressure
  end type ap_short

  type :: triangle                    ! triangle of vorticies resulting from 
                                      ! Delaunay-triangulation
    integer, dimension(4)    :: ind   ! indicies of its vertices
  end type 

  type :: adapt_set
    integer                  :: no_steps      ! number of iteration for implicite mixing
!    integer                  :: dim           ! 2d/3d mixing
    integer                  :: nlevs       ! number of levels
    real(prec)               :: lat_up        ! lower value of the considered pv-window    
    real(prec)               :: lat_down      ! upper value of the considered pv-window    
    real(prec)               :: lat_min       ! lower value of the considered lat-window    
    real(prec)               :: lat_max       ! upper value of the considered lat-window    
    real(prec)               :: lev_min     ! lower value of the considered lev_window
    real(prec)               :: lev_max     ! upper value of the considered lev_window
    real(prec)               :: lexp          ! critical Ljapunov exponent (1/day)
    real(prec)               :: timestep      
    real(prec), dimension(:), allocatable:: fac_min       ! mixing partition:     fac_min < r/r_mean < fac_max
    real(prec), dimension(:), allocatable :: fac_max       ! 
    real(prec)               :: fac_limit_outside  ! miximal distance between aps 
                                                   ! (as a factor of r_mean) which can be mixed
    real(prec)               :: fac_limit_inside
    real(prec)               :: fac_limit_lev_down
    real(prec)               :: fac_limit_lev_up
    real(prec)               :: fac_eliminate ! factor to change r_min in the elimination-loops 
    real(prec)               :: fac_bvf_min   ! minimal bvf value
    real(prec)               :: r_mean_c      ! mean distance between the APs in the coarse grid
    real(prec), dimension(:), allocatable :: r_min_c     ! lower limit (0 < r_min_c < r_mean) in the coarse grid
    real(prec), dimension(:), allocatable :: r_max_c     ! upper limit (r_mean < r_max_c) in the coarse grid
    real(prec), dimension(:), allocatable :: r_lim_c_inside ! total limit (r_max_c < r_lim_c) in the coarse grid
    real(prec), dimension(:), allocatable :: r_lim_c_outside
    real(prec)               :: r_mean_h      ! mean distance between the APs in the coarse grid
    real(prec), dimension(:), allocatable :: r_min_h       ! lower limit (0 < r_min_c < r_mean) 
                                              ! in the high resolution grid
    real(prec), dimension(:), allocatable :: r_max_h       ! upper limit (0 < r_min_c < r_mean)
                                              ! in the high resolution grid
    real(prec), dimension(:), allocatable :: r_lim_h_inside ! total limit (r_max_h < r_lim_h) in the high resol grid
    real(prec), dimension(:), allocatable :: r_lim_h_outside
    real(prec)               :: delta_lev   ! effective lev thickness
    real(prec)               :: interp_lev  ! adjustement factor for vertical interpolation
    real(prec)               :: r_dev         ! max. rel. deviation between n and n_old 
                                              ! (used for implicit mixing)
    real(prec), dimension(:),pointer  :: lev_grid  => NULL() ! vertical grid (lev values)
    real(prec), dimension(:),pointer  :: lev_delta => NULL() ! array of the layer-thickness
    real(prec), dimension(:),pointer  :: r_grid    => NULL() ! horizontal resolution 
 end type adapt_set

  type :: gridtype
     real(prec) :: lat0, lon0, dlat, dlon
     integer    :: nlat, nlon
  end type gridtype


  real(prec), dimension(:),  pointer :: lev_min_act, lev_max_act, lev_delta_act
  real(prec), dimension(:),  pointer :: l_min_act, l_max_act, l_delta_act

  integer,    dimension(:,:),pointer :: levelrange    ! (2,0:ntasks-1)

 
  real(prec)        :: lev_min, lev_max

  INTEGER           :: nlevs, l_nlevs

  CHARACTER(3)      :: vert_mix_param = 'WET'
  
  LOGICAL           :: ctrl_out = .false.  ! write additional control output (on standard output device)
  LOGICAL           :: dates30  = .false.


! op_pj_20160606+
!!$  type(adapt_set)   :: adapt_par
  type(adapt_set),SAVE   :: adapt_par
! op_pj_20160606-


  INTEGER                :: nmixspec
  CHARACTER(SPECNAMELEN), DIMENSION(maxspec) :: mixspec 

  
  !----------------------------------------------------------
  !     CTRL namelist:
  !----------------------------------------------------------
  
  REAL(PREC)        :: lexp              
  REAL(PREC)        :: fac_limit_outside
  REAL(PREC)        :: fac_limit_inside
  REAL(PREC)        :: fac_limit_lev_down
  REAL(PREC)        :: fac_limit_lev_up
  REAL(PREC)        :: fac_eliminate       = 1.0
  REAL(PREC)        :: fac_bvf_min         = 0.0001
  REAL(PREC)        :: delta_lev
  REAL(PREC)        :: r_dev

  INTEGER           :: no_steps     = 1   
  INTEGER           :: grid_switch  = 0        
  INTEGER           :: nintervals   = 1    
  INTEGER           :: timestep_mix = 24

  INTEGER           :: switch_mixing = 1  ! 0 = nothing
                                          ! 1 = implicit hor_mixing
                                          ! 2 = hor_mixing+vert_mixing

  

  

end Module messy_clamsmix_global
