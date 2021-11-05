! -*- f90 -*-
  USE mo_control,               ONLY : e5_ngl=>ngl, e5_nhgl=>nhgl, e5_nlon=>nlon
  
  IMPLICIT NONE
  SAVE

  REAL(dp), PARAMETER :: apzero = 100000.0_dp
  REAL(dp), PARAMETER :: h0     = 7._dp   ! km, scale height

  ! namelist variables (CPL_MBM)
  LOGICAL :: ladd_tte = .FALSE.    ! TRUE: add temperature tendency
  LOGICAL :: lv_echam = .FALSE.    ! TRUE: ECHAM vertical grid
  LOGICAL :: lzonal_mean = .FALSE. ! TRUE: use zonal mean representation 
  INTEGER :: nn   = 42             ! max meridional wave number for m=0.
  INTEGER :: nlev = 90             ! number of vertical levels
  !
  INTEGER :: nlon   = 0
  INTEGER :: nlat   = 0
  INTEGER :: nlevp1 = 0  ! nlev + 1
  INTEGER :: nvclev = 0  ! number of levels with vertical coefficients.
  INTEGER :: nglon  = 0 ! nlon   ! GP: NUMBER OF LON (1 if lcolumn=T)
  INTEGER :: nglat  = 0 ! nlat   ! GP: NUMBER OF LAT (1 if lcolumn=T)
  INTEGER :: nmp1   = 0 ! nn + 1 !   max zonal wave number + 1.
  INTEGER :: nsp    = 0 ! number of spectral coefficients.
  INTEGER :: nproma = 0
  INTEGER :: ngpblks= 0
  INTEGER :: npromz = 0
  INTEGER :: kproma = 0
  INTEGER :: jrow   = 0
  INTEGER, ALLOCATABLE, DIMENSION(:), TARGET :: glat  ! global latitude index N->S
  INTEGER, ALLOCATABLE, DIMENSION(:), TARGET :: glon  ! offset to global longitude
  INTEGER,  DIMENSION(1)     :: glats ! start values of latitudes
  INTEGER,  DIMENSION(1)     :: glate ! end values of latitudes
  INTEGER,  DIMENSION(1)     :: glons ! start values of longitudes
  INTEGER,  DIMENSION(1)     :: glone ! end values of longitudes

  REAL(dp), ALLOCATABLE, DIMENSION(:) :: vct ! vertical coefficients 

  LOGICAL, PARAMETER :: lmidatm = .TRUE.
