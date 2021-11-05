! -*- f90 -*-

  INTEGER :: nlev = 2         ! number of vertical levels: BL (level 2) and FT (level 1)
  INTEGER :: nlon = 1
  INTEGER :: nlat = 1

  INTEGER :: nproma           ! vector length
  INTEGER :: ngpblks          ! number of vector rows
  INTEGER :: npromz           ! vector length of last row
  INTEGER :: nlevp1           ! nlev +1 (0 number interface levels)
  INTEGER :: nvclev           ! number of levels with vertical coefficients

  INTEGER :: kproma           ! vector length of current row
  INTEGER :: jrow             ! row loop index

  REAL(dp), ALLOCATABLE :: vct(:) ! vertical coefficients table

  ! imitate ECHAM5
  LOGICAL :: lmidatm = .FALSE. ! switch for photolysis in mesosphere: not relevant for VERTICO
