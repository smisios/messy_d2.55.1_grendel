MODULE grid

USE domain
!use perf_mod

implicit none
SAVE
!INTEGER, PARAMETER :: crmvars=7

INTEGER, PARAMETER :: SHR_KIND_R8 = selected_real_kind(12) 
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
        
!INTEGER, PARAMETER :: nx = nx_gl/nsubdomains_x
!INTEGER, PARAMETER :: ny = ny_gl/nsubdomains_y 
!INTEGER, PARAMETER :: nz = nz_gl+1
!INTEGER, PARAMETER :: nzm = nz-1
        
!INTEGER, PARAMETER :: nsubdomains = nsubdomains_x * nsubdomains_y

!LOGICAL, PARAMETER :: RUN3D = ny_gl.gt.1
!LOGICAL, PARAMETER :: RUN2D = .not.RUN3D

!INTEGER, PARAMETER :: nxp1 = nx + 1
!INTEGER, PARAMETER :: nyp1 = ny + 1 * YES3D
!INTEGER, PARAMETER :: nxp2 = nx + 2
!INTEGER, PARAMETER :: nyp2 = ny + 2 * YES3D
!INTEGER, PARAMETER :: nxp3 = nx + 3
!INTEGER, PARAMETER :: nyp3 = ny + 3 * YES3D
!INTEGER, PARAMETER :: nxp4 = nx + 4
!INTEGER, PARAMETER :: nyp4 = ny + 4 * YES3D

!INTEGER, PARAMETER :: dimx1_u = -1
!INTEGER, PARAMETER :: dimx2_u = nxp3
!INTEGER, PARAMETER :: dimy1_u = 1-2*YES3D
!INTEGER, PARAMETER :: dimy2_u = nyp2
!INTEGER, PARAMETER :: dimx1_v = -1
!INTEGER, PARAMETER :: dimx2_v = nxp2
!INTEGER, PARAMETER :: dimy1_v = 1-2*YES3D
!INTEGER, PARAMETER :: dimy2_v = nyp3
!INTEGER, PARAMETER :: dimx1_w = -1
!INTEGER, PARAMETER :: dimx2_w = nxp2
!INTEGER, PARAMETER :: dimy1_w = 1-2*YES3D
!INTEGER, PARAMETER :: dimy2_w = nyp2
!INTEGER, PARAMETER :: dimx1_s = -2
!INTEGER, PARAMETER :: dimx2_s = nxp3
!INTEGER, PARAMETER :: dimy1_s = 1-3*YES3D
!INTEGER, PARAMETER :: dimy2_s = nyp3

!INTEGER, PARAMETER :: ncols = nx*ny

INTEGER :: crmvars

INTEGER :: nx 
INTEGER :: ny
INTEGER :: nz
INTEGER :: nzm
INTEGER :: nsubdomains
        
LOGICAL :: RUN3D
LOGICAL :: RUN2D

INTEGER :: nxp1
INTEGER :: nyp1
INTEGER :: nxp2
INTEGER :: nyp2
INTEGER :: nxp3
INTEGER :: nyp3
INTEGER :: nxp4
INTEGER :: nyp4

INTEGER :: dimx1_u
INTEGER :: dimx2_u
INTEGER :: dimy1_u
INTEGER :: dimy2_u
INTEGER :: dimx1_v
INTEGER :: dimx2_v
INTEGER :: dimy1_v
INTEGER :: dimy2_v
INTEGER :: dimx1_w
INTEGER :: dimx2_w
INTEGER :: dimy1_w
INTEGER :: dimy2_w
INTEGER :: dimx1_s
INTEGER :: dimx2_s
INTEGER :: dimy1_s
INTEGER :: dimy2_s

INTEGER :: ncols

REAL(dp) :: dx ! grid spacing in x direction
REAL(dp) :: dy ! grid spacing in y direction
REAL(dp) :: dz ! grid spacing in z direction

REAL(dp), POINTER, DIMENSION(:) :: z      => NULL() !  height of the pressure levels above surface (m)
REAL(dp), POINTER, DIMENSION(:) :: pres   => NULL() !  pressure (mb) at scalar levels
REAL(dp), POINTER, DIMENSION(:) :: zi     => NULL() !  height of the interface levels (m)
REAL(dp), POINTER, DIMENSION(:) :: presi  => NULL() !  pressure (mb) at interface levels

REAL(dp), POINTER, DIMENSION(:) :: adz    => NULL() !  ratio of the grid spacing to dz for pressure levels
REAL(dp), POINTER, DIMENSION(:) :: adzw   => NULL() !  ratio of the grid spacing to dz for w levels
REAL(dp), POINTER, DIMENSION(:) :: grdf_x => NULL() !  grid factor for eddy diffusion in x
REAL(dp), POINTER, DIMENSION(:) :: grdf_y => NULL() !  grid factor for eddy diffusion in y
REAL(dp), POINTER, DIMENSION(:) :: grdf_z => NULL() !  grid factor for eddy diffusion in z

REAL(dp) ::  at,bt,ct ! coefficients for the Adams-Bashforth scheme 
REAL(dp) ::  dt       ! dynamical timestep
REAL(dp) ::  dtn      ! current dynamical timestep (can be smaller than dt)
REAL(dp) ::  dt3(3)   ! dynamical timesteps for three most recent time steps
REAL(dp) ::  time     ! current time in sec.
REAL(dp) ::  day0     ! starting day (including fraction)
REAL(dp) ::  day      ! current day (including fraction)
REAL(dp) ::  dtfactor ! dtn/dt

INTEGER nstep   ! current number of performed time steps 
INTEGER nstop   ! time step number to stop the integration
INTEGER nelapse ! time step number to elapse before stoping
INTEGER na, nb, nc ! indeces for swapping the rhs arrays for AB scheme
INTEGER ncycle  ! number of subcycles over the dynamical timestep
INTEGER icycle  ! current subcycle 
INTEGER nadams  ! the order of the AB scheme (should be kept at 3)        
INTEGER nstat   ! the interval in time steps to compute statistics
INTEGER nstatis ! the interval between substeps to compute statistics
INTEGER nstatfrq! frequency of computing statistics 
INTEGER nprint  ! frequency of printing a listing (steps)
INTEGER nrestart! switch to control starting/restarting of the model
LOGICAL restart_sep ! write separate restart files for sub-domains
INTEGER nrestart_skip ! number of skips of writing restart (default 0)
LOGICAL output_sep ! write separate 3D and 2D files for sub-domains
INTEGER nrad    ! frequency of calling the radiation routines
LOGICAL save3Dbin ! save 3D data in binary format(no 2-byte compression)
LOGICAL save3Dsep ! use separate file for each time point for2-model
LOGICAL save2Dsep !write a separate file for each time point for 2D horizontal fields
LOGICAL save2Davg !flag to time-average 2D output fields (default .false.)
INTEGER nsave3D ! frequency of writting 3D fields (steps)
INTEGER nsave3Dstart ! timestep to start writting 3D fields
INTEGER nsave3Dend   ! timestep to end writting 3D fields
real    qnsave3D !threshold manimum cloud water(kg/kg) to save 3D fields
LOGICAL dogzip3D ! gzip compress a 3D output file   
LOGICAL dogzip2D ! gzip compress a 2D output file if save2Dsep=.true.   
LOGICAL save2Dbin !save 2D data in binary format, rather than compressed
INTEGER nsave2D ! frequency of writting 2D fields (steps)
INTEGER nsave2Dstart ! timestep to start writting 2D fields
INTEGER nsave2Dend   ! timestep to end writting 2D fields
character *40 caseid! 8-symbol id-string to identify a run      
character *40 case  ! 8-symbol id-string to identify a case-name        
LOGICAL dostatis! flag to permit the gathering of statistics
LOGICAL dostatisrad! flag to permit the gathering of radiation statistics
INTEGER nensemble ! the number of subensemble set of perturbations
LOGICAL notopened2D ! flag to see if the 2D output datafile is opened   
LOGICAL notopened3D ! flag to see if the 3D output datafile is opened   
character *256 rundatadir ! path to data directory containing data files needed to run
INTEGER perturb_type  ! type of initial noise in setperturb()

!   Flags:

LOGICAL CEM     ! flag for Cloud Ensemble Model
LOGICAL LES     ! flag for Large-Eddy Simulation
LOGICAL OCEAN   ! flag indicating that surface is water
LOGICAL LAND    ! flag indicating that surface is land
LOGICAL SFC_FLX_FXD  ! surface sensible flux is fixed
LOGICAL SFC_TAU_FXD! surface drag is fixed

!       Multitasking staff:     
          
INTEGER rank   ! rank of the current subdomain task (default 0) 
INTEGER ranknn ! rank of the "northern" subdomain task
INTEGER rankss ! rank of the "southern" subdomain task
INTEGER rankee ! rank of the "eastern"  subdomain task
INTEGER rankww ! rank of the "western"  subdomain task
INTEGER rankne ! rank of the "north-eastern" subdomain task
INTEGER ranknw ! rank of the "north-western" subdomain task
INTEGER rankse ! rank of the "south-eastern" subdomain task
INTEGER ranksw ! rank of the "south-western" subdomain task
LOGICAL dompi  ! LOGICAL switch to do multitasking
LOGICAL masterproc ! .true. if rank.eq.0 
        
!   LOGICAL switches and flags:

LOGICAL   dodamping, doupperbound, docloud, doprecip, &
          dolongwave, doshortwave, dosgs, dosubsidence, &
          docoriolis, dosurface, dolargescale, doradforcing, &
          dosfcforcing, doradsimple, donudging_uv, donudging_tq, & 
          dosmagor, doscalar, doensemble, doxy, dowallx, dowally, docup, &
          docolumn, doperpetual, doseasons, doradhomo, dosfchomo, &
          doisccp, dodynamicocean, dosolarconstant, dotracers, dosmoke 

! For dosolarconstant simulations, allow solar constant and zenith
! angle to be specified individually
REAL(dp) :: solar_constant  ! solar constant (in W/m2)
REAL(dp) :: zenith_angle    ! zenith angle (in degrees)

LOGICAL doSAMconditionals, dosatupdnconditionals
!bloss: option for reading data from a SCAM netcdf input file
LOGICAL doscamiopdata
LOGICAL :: isInitialized_scamiopdata = .false.
character(len=120) iopfile
LOGICAL dozero_out_day0 ! set day0 to zero  for ideal cases

! SCAM uses omega instead of w for large-scale vertical motion.  
!   Converted to w in forcing()
LOGICAL :: wgls_holds_omega = .false.

INTEGER nstatmom ! frequency of writting statistical moment fields (steps)
INTEGER nstatmomstart ! timestep to start writting statistical moment fields
INTEGER nstatmomend   ! timestep to end writting statistical moment fields
LOGICAL notopenedmom ! flag to see if the statistical moment file is opened
LOGICAL savemomsep ! use one file with stat moments  for each time point for 2-model runs
LOGICAL savemombin ! save statistical moment data in binary format(no 2-byte compression)

INTEGER nmovie ! frequency of writting movie fields (steps)
INTEGER nmoviestart ! timestep to start writting statistical moment fields
INTEGER nmovieend   ! timestep to end writting statistical moment fields

!integer, parameter :: npressureslabs = nsubdomains
integer :: npressureslabs

integer :: nzslab
integer :: nx2, ny2
integer :: n3i, n3j

END MODULE grid

! ====================================================================

  SUBROUTINE set_crm_dims(crm_nx_nml,crm_ny_nml,crm_nlev_nml,crm_size_nml,crm_time_nml,crm_Y3D_nml, crm_nvars_nml)
  
    use grid

    IMPLICIT NONE

    INTEGER  :: crm_nx_nml, crm_ny_nml, crm_nlev_nml, crm_nvars_nml
    REAL(dp) :: crm_size_nml, crm_time_nml
    INTEGER  :: crm_Y3D_nml

    crmvars = crm_nvars_nml

    crm_nx = crm_nx_nml
    crm_ny = crm_ny_nml
    crm_nz = crm_nlev_nml
    
    crm_dx = crm_size_nml
    crm_dy = crm_size_nml
    crm_dt = crm_time_nml

    dx     = crm_dx
    dy     = crm_dy
    dt     = crm_dt

    YES3DVAL = crm_Y3D_nml

    YES3D = YES3DVAL   ! Domain dimensionality: 1 - 3D, 0 - 2D
    nx_gl = crm_nx     ! Number of grid points in X
    ny_gl = crm_ny     ! Number of grid points in Y
    nz_gl = crm_nz 

    nsubdomains_x  = 1
    nsubdomains_y  = 1
    nsubdomains    = nsubdomains_x * nsubdomains_y
    npressureslabs = nsubdomains
    
    nx  = nx_gl/nsubdomains_x
    ny  = ny_gl/nsubdomains_y 
    nz  = nz_gl+1
    nzm = nz-1

    RUN3D = ny_gl.gt.1
    RUN2D = .not.RUN3D

    nxp1 = nx + 1
    nyp1 = ny + 1 * YES3D
    nxp2 = nx + 2
    nyp2 = ny + 2 * YES3D
    nxp3 = nx + 3
    nyp3 = ny + 3 * YES3D
    nxp4 = nx + 4
    nyp4 = ny + 4 * YES3D

    dimx1_u = -1
    dimx2_u = nxp3
    dimy1_u = 1-2*YES3D
    dimy2_u = nyp2
    dimx1_v = -1
    dimx2_v = nxp2
    dimy1_v = 1-2*YES3D
    dimy2_v = nyp3
    dimx1_w = -1
    dimx2_w = nxp2
    dimy1_w = 1-2*YES3D
    dimy2_w = nyp2
    dimx1_s = -2
    dimx2_s = nxp3
    dimy1_s = 1-3*YES3D
    dimy2_s = nyp3

    ncols = nx*ny
    
    nzslab = max(1,nzm / npressureslabs)
    nx2=nx_gl+2
    ny2=ny_gl+2*YES3D
    n3i=3*nx_gl/2+1
    n3j=3*ny_gl/2+1

  END SUBROUTINE set_crm_dims

! ====================================================================
