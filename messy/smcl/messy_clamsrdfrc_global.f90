MODULE messy_clamsrdfrc_global

  USE messy_clams_global, ONLY: PREC

  !    ap_g_all:   lat, lon, lev, time fuer alle x (inneres Gitter),
  !                wird aus file_g eingelesen
  !    ap_s_all:   lat, lon, lev, time fuer alle o (aeuﬂeres Gitter),
  !                wird aus file_s eingelesen
  !    data_s_all: Parameterwerte fuer alle o (z.B. CH4, O3,...),
  !                wird aus file_s eingelesen

  TYPE :: ap                            
     REAL(PREC)                 :: lon        ! its longitude (UKMO-notation)
     REAL(PREC)                 :: lat        ! its latitude  (UKMO-notation)
     REAL(PREC)                 :: lev      ! its pot. temp.
     REAL(KIND(0d0))      :: time       ! its julian sec  
  END TYPE ap
  
  TYPE :: datetype
     INTEGER :: year, month, day, sec
  END TYPE datetype

  TYPE(ap), DIMENSION(:), POINTER :: ap_g_all, ap_s_all

  REAL(PREC), DIMENSION(:,:), POINTER :: data_s_all
  REAL(PREC), DIMENSION(:,:), POINTER :: data_g_tot

  REAL(PREC), PARAMETER :: earth_radius=6371.
  
  REAL(PREC), PARAMETER :: mdi=-1.E30  ! missing data indicator

  INTEGER, PARAMETER :: max_nb=150 ! max. Anz. angrenzender Triangeln

!!!!!
!  integer, parameter :: strlen = 20
  INTEGER, PARAMETER :: strlen = 10

  REAL(PREC)    :: lat_min, lat_max  ! aus i3d.inp eingelesen
  REAL(PREC)    :: lat_down, lat_up, r_coarse, r_high ! aus init-file (file_s)

  REAL(PREC)    :: max_dist

  !    npart_g:    Dimension von ap_g_all
  !    npart_s:    Dimension von ap_s_all
  !    dim_tr:     Number of triangles
  INTEGER :: npart_s, npart_g, nwts, ntags, dim_tr

  INTEGER :: weighted

  LOGICAL :: north_hemisphere, one_file, ctrl_out=.TRUE.

  ! file_s: Initial positions of APs and their dynamical properties 
  ! file_g: experiment data set (aircraft or satellite data)
  CHARACTER(3)   :: prefix_s
  CHARACTER(160) :: file_s, file_g, file_out_dir, file_out, dir_s
  CHARACTER(160) :: file_g2 !EMAC

  REAL(PREC),       DIMENSION(:,:), POINTER :: lev_window

END MODULE messy_clamsrdfrc_global
