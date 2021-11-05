! -*- f90 -*-
  ! LOCALIZED PARAMETERS
  ! - GRID CONTROL HERE FIXED VALUES normally determined by grid definition
  INTEGER :: nlev = 1          ! number of vertical levels
  INTEGER :: nlon = 1
  INTEGER :: nlat = 1

  INTEGER :: nproma           ! vector length
  INTEGER :: ngpblks          ! number of vector rows
  INTEGER :: npromz           ! vector length of last row

  INTEGER :: kproma           ! vector length of current row
  INTEGER :: jrow             ! row loop index
  INTEGER :: je               !

