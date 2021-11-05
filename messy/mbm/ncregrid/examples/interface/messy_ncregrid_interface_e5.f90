! -*- f90 -*-
! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_INTERFACE
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: NCREGRID_SET_MESSAGEMODE
  PUBLIC :: INTERFACE_GEOHYBGRID

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE NCREGRID_SET_MESSAGEMODE

  USE MESSY_NCREGRID_BASE, ONLY: MSGMODE &
                               , MSGMODE_S, MSGMODE_E, MSGMODE_VL  &
                               , MSGMODE_W, MSGMODE_VM, MSGMODE_I

  IMPLICIT NONE

  ! RESET MESSAGE MODUS
  MSGMODE = MSGMODE_S + MSGMODE_E + MSGMODE_VL ! &
!          + MSGMODE_W + MSGMODE_VM + MSGMODE_I

END SUBROUTINE NCREGRID_SET_MESSAGEMODE
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INTERFACE_GEOHYBGRID(g, ok)

  ! ECHAM5/MESSy
!!$  USE messy_main_bmluse_bi,    ONLY: get_date_components  &
!!$                                   , next_date, start_date, TC_get
  USE messy_main_data_bi,      ONLY: nvclev, vct, nlon, ngl, nhgl, nlev &
                                   , gl_gmu, aps 
  USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START &
                                   , HOUR_START, MINUTE_START, SECOND_START &
                                   , YEAR, MONTH, DAY &
                                   , HOUR, MINUTE, SECOND

  ! MESSy
  USE messy_main_constants_mem, ONLY: api=>pi
  USE messy_main_timer,         ONLY: time_span_d

  ! NCREGRID MODULES
  USE MESSY_NCREGRID_BASE
  USE MESSY_NCREGRID_NETCDF
  USE MESSY_NCREGRID_GEOHYB

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(OUT) :: g
  LOGICAL          , INTENT(OUT) :: ok

  !LOCAL
  REAL(DP) :: dts
  INTEGER  :: i,j
  INTEGER  :: status
  CHARACTER(LEN=100)       :: tunit
  REAL(DP), DIMENSION(4,2) :: ranges

  ! INIT
  CALL INIT_GEOHYBGRID(g)

  g%file = 'ECHAM_5_GEO-HYBRID-GRID'       ! Filename
  g%t    = 0                               ! time step

  ! LONGITUDE (MID) ...
  g%lonm%name  = 'lon'
  g%lonm%id    = NULL_VARID
  g%lonm%xtype = NF90_FLOAT
  ! ... dimensions
  g%lonm%ndims = 1
  ALLOCATE(g%lonm%dim(g%lonm%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
  g%lonm%dim(1)%name  = 'lon'
  g%lonm%dim(1)%id    = NULL_DIMID
  g%lonm%dim(1)%len   = nlon
  g%lonm%dim(1)%fuid  = .false.
  g%lonm%dim(1)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(g%lonm%dat, g%lonm%ndims, (/g%lonm%dim(1)%len/) &
                   ,VTYPE_REAL)
  DO i=1, nlon
     g%lonm%dat%vr(i) = 360. * REAL(i-1) / REAL(nlon)
  END DO
  ! ... attributes
  CALL ADD_NCATT(g%lonm, 'long_name', vs='longitude')
  CALL ADD_NCATT(g%lonm, 'units', vs='degrees_east')
  ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
  ranges(1,2) = RGEMPTY

  ! LATITUDE (MID) ...
  g%latm%name  = 'lat'
  g%latm%id    = NULL_VARID
  g%latm%xtype = NF90_FLOAT
  ! ... dimensions
  g%latm%ndims = 1
  ALLOCATE(g%latm%dim(g%latm%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
  g%latm%dim(1)%name  = 'lat'
  g%latm%dim(1)%id    = NULL_DIMID
  g%latm%dim(1)%len   = ngl
  g%latm%dim(1)%fuid  = .false.
  g%latm%dim(1)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(g%latm%dat, g%latm%ndims, (/g%latm%dim(1)%len/) &
                   ,VTYPE_REAL)
  DO j=1,nhgl
     g%latm%dat%vr(j)=ASIN(gl_gmu(j))*180./api
     g%latm%dat%vr(ngl+1-j)=-g%latm%dat%vr(j)
  END DO
  ! ... attributes
  CALL ADD_NCATT(g%latm, 'long_name', vs='latitude')
  CALL ADD_NCATT(g%latm, 'units', vs='degrees_north')
  ranges(2,1) = -90.0_dp
  ranges(2,2) =  90.0_dp

  ! HYBRID-A-COEFFICIENTS (INTERFACES) ...
  g%hyai%name  = 'hyai'
  g%hyai%id    = NULL_VARID
  g%hyai%xtype = NF90_FLOAT
  ! ... dimensions
  g%hyai%ndims = 1
  ALLOCATE(g%hyai%dim(g%hyai%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,3)
  g%hyai%dim(1)%name  = 'ilev'
  g%hyai%dim(1)%id    = NULL_DIMID
  g%hyai%dim(1)%len   = nlev+1
  g%hyai%dim(1)%fuid  = .false.
  g%hyai%dim(1)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(g%hyai%dat, g%hyai%ndims, (/g%hyai%dim(1)%len/) &
                   ,VTYPE_REAL)
  g%hyai%dat%vr(:) = vct(1:nvclev)
  ! ... attributes
  CALL ADD_NCATT(g%hyai, 'long_name'                              &
               ,vs='hybrid-A-coefficients at layer interfaces')
  CALL ADD_NCATT(g%hyai, 'units', vs='1')
  ranges(3,1) = 0.0_dp
  ranges(3,2) = 0.0_dp

  ! HYBRID-B-COEFFICIENTS (INTERFACES) ...
  g%hybi%name  = 'hybi'
  g%hybi%id    = NULL_VARID
  g%hybi%xtype = NF90_FLOAT
  ! ... dimensions
  g%hybi%ndims = 1
  ALLOCATE(g%hybi%dim(g%hybi%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,4)
  g%hybi%dim(1) = g%hyai%dim(1)
  ! ... data
  CALL INIT_NARRAY(g%hybi%dat, g%hybi%ndims, (/g%hybi%dim(1)%len/) &
                   ,VTYPE_REAL)
  g%hybi%dat%vr(:) = vct(nvclev+1:2*nvclev)
  ! ... attributes
  CALL ADD_NCATT(g%hybi, 'long_name'                              &
               ,vs='hybrid-B-coefficients at layer interfaces')
  CALL ADD_NCATT(g%hybi, 'units', vs='1')
  ranges(4,1) = 0.0_dp
  ranges(4,2) = 1.0_dp

  ! SURFACE PRESSURE
  g%ps%name  = 'ps'
  g%ps%id    = NULL_VARID
  g%ps%xtype = NF90_FLOAT
  ! ... dimensions
  g%ps%ndims = 2
  ALLOCATE(g%ps%dim(g%ps%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)
  g%ps%dim(1) = g%lonm%dim(1)
  g%ps%dim(2) = g%latm%dim(1)
  ! ... data
  CALL INIT_NARRAY(g%ps%dat, g%ps%ndims                  &
                   ,(/g%ps%dim(1)%len, g%ps%dim(2)%len/) &
                   ,VTYPE_REAL)
  IF (ASSOCIATED(aps)) THEN
     DO i=1, g%ps%dim(1)%len
        DO j=1, g%ps%dim(2)%len
           g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                                 ,(/i,j/)))                           &
!qqq
!           = aps(i,j)
           = 101325.0
        END DO
     END DO
  ELSE
     DO i=1, g%ps%dim(1)%len
        DO j=1, g%ps%dim(2)%len
           g%ps%dat%vr(POSITION((/g%ps%dim(1)%len, g%ps%dim(2)%len/)  &
                                 ,(/i,j/)))                           &
           = 101325.0
        END DO
     END DO
  END IF
  ! ... attributes
  CALL ADD_NCATT(g%ps, 'long_name', vs='surface pressure')
  CALL ADD_NCATT(g%ps, 'units', vs='Pa')

  ! REFERENCE PRESSURE ...
  g%p0%name  = 'p0'
  g%p0%id    = NULL_VARID
  g%p0%xtype = NF90_FLOAT
  ! ... dimensions
  g%p0%ndims = 0
  ! ... data
  CALL INIT_NARRAY(g%p0%dat, 1, (/ 1 /), VTYPE_REAL)
  g%p0%dat%vr(1) = 1.0
  ! ... attributes
  CALL ADD_NCATT(g%p0, 'long_name', vs='reference pressure')
  CALL ADD_NCATT(g%p0, 'units', vs='Pa')

  ! TIME (MID) ...
  g%timem%name  = 'time'
  g%timem%id    = NULL_VARID
  g%timem%xtype = NF90_FLOAT
  ! ... dimensions
  g%timem%ndims = 1
  ALLOCATE(g%timem%dim(g%timem%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
  g%timem%dim(1)%name  = 'time'
  g%timem%dim(1)%id    = NULL_DIMID
  g%timem%dim(1)%len   = 1
  g%timem%dim(1)%fuid  = .true.
  g%timem%dim(1)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(g%timem%dat, g%timem%ndims, (/g%timem%dim(1)%len/) &
                   ,VTYPE_DOUBLE)
  ! TIME: SECONDS SINCE MODEL START
  ! mz_pj_20080327+
!!$  CALL TC_get(start_date, start_day, start_sec)
!!$  CALL TC_get(next_date,  next_day,  next_sec)
!!$  g%timem%dat%vd(1) = 86400.*(next_day-start_day)+next_sec-start_sec
! mz_pj_20090519+
!!$  CALL time_span_s(dts   &
!!$       , YEAR_START, MONTH_START, DAY_START &
!!$       , HOUR_START, MINUTE_START, SECOND_START  &
!!$       , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
!!$  g%timem%dat%vd(1) = REAL(dts, dp)  
  CALL time_span_d(dts   &
       , YEAR_START, MONTH_START, DAY_START &
       , HOUR_START, MINUTE_START, SECOND_START  &
       , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
  g%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds
! mz_pj_20090519-
  ! ... attributes
  CALL ADD_NCATT(g%timem, 'long_name'                              &
               ,vs='time in seconds since model start')
  ! mz_pj_20080327+
!!$  CALL get_date_components(start_date,pyr,pmo,pdy,phr,pmn,pse) ! mz_pj_20080327
!!$  WRITE(tunit, &
!!$       '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
!!$       pyr,pmo,pdy,phr,pmn,pse
  WRITE(tunit, &
       '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
       YEAR_START, MONTH_START, DAY_START &
       , HOUR_START, MINUTE_START, SECOND_START
  ! mz_pj_20080327-
  CALL ADD_NCATT(g%timem, 'units', vs=TRIM(tunit))

  ! CALCULATE INTs from MIDs
  CALL COMPLETE_GEOHYBGRID(g, ranges)

  ! INTERFACE OK
  ok = .true.

END SUBROUTINE INTERFACE_GEOHYBGRID
! ------------------------------------------------------------------

! ******************************************************************
! ------------------------------------------------------------------
END MODULE MESSY_NCREGRID_INTERFACE
! ------------------------------------------------------------------
! ******************************************************************
