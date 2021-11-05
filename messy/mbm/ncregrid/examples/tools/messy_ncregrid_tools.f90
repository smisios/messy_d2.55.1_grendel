#ifndef _SX

! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_TOOLS
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, October 2002
! ******************************************************************

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  INTEGER,           PARAMETER   :: NCCNTMAXNLEN = 40                  !   ! 
  CHARACTER(LEN=24), PARAMETER   :: fstr = '(a80,1x,a40,1x,4(i8,1x))'  ! <-V
  INTEGER,           PARAMETER   :: NCCRSTRL = 80                      !   |

  ! TYPE DECLARATION TO HOLD COUNTER INFORMATION
  TYPE NCRGCNT
     CHARACTER(LEN=NCCNTMAXNLEN) :: name = ''
     INTEGER                     :: start = 1
     INTEGER                     :: step  = 1
     INTEGER                     :: reset = 1
     INTEGER                     :: current = 1
  END TYPE NCRGCNT
  
  ! LIST TO STORE COUNTER INFORMATION FOR RESTART
  TYPE T_NCRGCNT_LIST
     CHARACTER(NCCRSTRL)           :: mname  = ''
     TYPE(NCRGCNT)                 :: this
     TYPE(T_NCRGCNT_LIST), POINTER :: next => NULL()
  END TYPE T_NCRGCNT_LIST

  LOGICAL, PUBLIC, SAVE :: LMODE2DH = .TRUE.

  ! ====================================================================
  ! CONCAT. LIST OF RGT-EVENTS
  TYPE(T_NCRGCNT_LIST), POINTER, SAVE :: GRGTLIST => NULL()
  ! ====================================================================
  
  PUBLIC  :: RGTOOL_CONVERT     ! CONVERTS NCVAR INTO 4D-ARRAY
  PUBLIC  :: RGTOOL_CONVERT_DAT2VAR ! CONVERTS 4D-ARRAY INTO NCVAR
  PUBLIC  :: RGTOOL_READ_NCVAR  ! NCREGRID ONE STEP OF ONE FIELD
  PUBLIC  :: RGTOOL_READ_NCFILE ! NCREGRID ONE STEP OF ONE FILE
  PUBLIC  :: RGTOOL_G2C         ! CONVERTS GRID-INFORMATION INTO ARRAYS

  ! COUNTER INFORMATION HANDLING
  PUBLIC  :: NCRGCNT            ! TYPE STRUCT FOR COUNTER INFORMATION I/O
  PUBLIC  :: RGTOOL_NCRGCNT_RST ! MANAGES COUNTER INFORMATION RESTART I/O
  !
  PUBLIC  :: CLEAN_NCRGCNT_LIST
  PUBLIC  :: WRITE_NCRGCNT_LIST
  PUBLIC  :: READ_NCRGCNT_LIST
  PUBLIC  :: GET_NEXT_NCRGCNT
  !
  !PRIVATE :: loc_ncrgcnt
  !PRIVATE :: new_ncrgcnt
  !PRIVATE :: ncrgcnt_error_str

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_CONVERT(var, dat, grid, order)

    ! CONVERTS NCREGRID OUTPUT OF TYPE N-ARRAY (var) ON GRID
    ! grid TO 4D-ARRAY (dat)
    ! THE OPTIONAL STRING order DEFINES THE ORDER OF DIMENSIONS
    ! (DEFAULT: 'xyzn')
    !
    ! Author: Patrick Joeckel, MPIC, October 2002

    ! REGRID MODULES
    USE MESSY_NCREGRID_BASE
    USE MESSY_NCREGRID_NETCDF
    USE MESSY_NCREGRID_GEOHYB

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, INDEX, PRESENT, PRODUCT, REAL, SIZE, SUM, TRIM

    ! I/O
    TYPE (ncvar),      INTENT(IN)           :: var   ! nc-variable
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: dat   ! DATA ON DESTINATION GRID
    TYPE (geohybgrid), INTENT(IN)           :: grid  ! grid information
    CHARACTER(LEN=4),  INTENT(IN), OPTIONAL :: order ! DEFAULT: 'xyzn'

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_CONVERT'
    INTEGER :: ovec(4)          ! order of (x,y,z,n) for output field
    INTEGER :: ldvar(4)         ! ordered dim-lenghts in var
    INTEGER :: dpos(3)          ! position of lon, lat, lev in variable
    INTEGER :: nrvdim           ! number of recognized dimensions in variable
    INTEGER :: npdlv            ! product of dimension-lengths in variable
    CHARACTER(LEN=6) :: ostr(4) ! order-string for output
    INTEGER :: i                ! counter
    INTEGER :: nvec(4)          ! counter parameter for dimensions
    INTEGER :: vtype            ! type of variable
    INTEGER, DIMENSION(:), POINTER :: vec   ! element vector
    INTEGER :: status           ! status flag

    ! INIT
    NULLIFY(vec)
    ovec = (/ 1,2,3,4 /)                               ! DEFAULT
    IF (LMODE2DH) THEN
       ostr = (/ 'x/LON ', 'y/LAT ', 'z/LEV ', 'n/PAR '/) ! DEFAULT
    ELSE
       ostr = (/ 'x/COL ', 'y/DUM ', 'z/LEV ', 'n/PAR '/) ! DEFAULT
    END IF
    IF (PRESENT(order)) THEN
       IF (LMODE2DH) THEN
          ovec(1) = INDEX(order, 'x')
          ostr(ovec(1)) = 'x/LON '
          ovec(2) = INDEX(order, 'y')
          ostr(ovec(2)) = 'y/LAT '
       ELSE
          ovec(1) = INDEX(order, 'x')
          ostr(ovec(1)) = 'x/COL '
          ovec(2) = INDEX(order, 'y')
          ostr(ovec(2)) = 'y/DUM '
       END IF
       ovec(3) = INDEX(order, 'z')
       ostr(ovec(3)) = 'z/LEV '
       ovec(4) = INDEX(order, 'n')
       ostr(ovec(4)) = 'n/PAR '
       IF (SUM(ovec) /= 10) THEN
          CALL RGMSG(substr, RGMLE, &
               'ERROR IN '''//TRIM(order)//'''-STRING AT SUBROUTINE CALL !')
       END IF
    END IF
    !
    IF (ASSOCIATED(dat)) THEN
       DEALLOCATE(dat, STAT=status)
       CALL ERRMSG(substr,status,1)
    END IF
    NULLIFY(dat)
    !
    dpos(:) = 0
    nrvdim = 0
    ldvar(:) = 1
    IF (ASSOCIATED(var%dat%dim)) THEN
       npdlv = PRODUCT(var%dat%dim)
    ELSE
       RETURN
    END IF

    ! LOOP OVER VARIABLE DIMENSIONS AND CHECK NAMES, LENGTHS, AND POSITIONS
    DO i=1, var%ndims
       IF (LMODE2DH) THEN
          IF (ASSOCIATED(grid%lonm%dim)) THEN
             IF (QCMP_NCDIM(var%dim(i), grid%lonm%dim(1)) > 1) THEN
                nrvdim = nrvdim+1
                ldvar(ovec(1)) = var%dim(i)%len
                dpos(1) = i
             END IF
          END IF
          IF (ASSOCIATED(grid%latm%dim)) THEN
             IF (QCMP_NCDIM(var%dim(i), grid%latm%dim(1)) > 1) THEN
                nrvdim = nrvdim+1
                ldvar(ovec(2)) = var%dim(i)%len
                dpos(2) = i
             END IF
          END IF
       ELSE
          IF (ASSOCIATED(grid%col%dim)) THEN
             IF (QCMP_NCDIM(var%dim(i), grid%col%dim(1)) > 1) THEN
                nrvdim = nrvdim+1
                ldvar(ovec(1)) = var%dim(i)%len
                dpos(1) = i
             END IF
          END IF
       END IF

       IF (ASSOCIATED(grid%hyam%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%hyam%dim(1)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(3)) = var%dim(i)%len
             dpos(3) = i
          END IF
       ELSE
          IF (ASSOCIATED(grid%hybm%dim)) THEN
             IF (QCMP_NCDIM(var%dim(i), grid%hybm%dim(1)) > 1) THEN
                nrvdim = nrvdim+1
                ldvar(ovec(3)) = var%dim(i)%len
                dpos(3) = i
             END IF
          END IF
       END IF
    END DO

    ldvar(ovec(4)) = npdlv/(ldvar(ovec(1))*ldvar(ovec(2))*ldvar(ovec(3)))

    CALL RGMSG(substr, RGMLI, &
         'VARIABLE '''//TRIM(var%name)//''' IS ',var%ndims, &
         '-DIMENSIONAL')
    CALL RGMSG(substr, RGMLIC, '(',var%dat%dim,')')
    CALL RGMSG(substr, RGMLIC, &
         'DATA DIMENSION-LENGHTS ARE: ')
    DO i=1, SIZE(ostr)
       CALL RGMSG(substr, RGMLIC, '  '//ostr(i)//': ',ldvar(i),' ')
    END DO
    IF (LMODE2DH) THEN
       CALL RGMSG(substr, RGMLIC, &
            'SOURCE POSITIONS (x/LON, y/LAT, z/LEV) ARE: ',dpos,' ')
    ELSE
       CALL RGMSG(substr, RGMLIC, &
            'SOURCE POSITIONS (x/COL, y/DUM, z/LEV) ARE: ',dpos,' ')
    END IF

    ! ALLOCATE DATA SPACE
    ALLOCATE(dat(ldvar(1),ldvar(2),ldvar(3),ldvar(4)), STAT=status)
    CALL ERRMSG(substr,status,1)

    ! COPY DATA (LOOP OVER ELEMENTS)
    vtype = var%dat%type
    DO i=1, npdlv
       CALL ELEMENT(var%dat%dim, i, vec)
       nvec((ovec(4))) = PRODUCT(vec)
       IF (dpos(1) /= 0) THEN
          nvec((ovec(1))) = vec(dpos(1))
          nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(1))
       ELSE
          nvec(ovec(1)) = 1
       END IF
       IF (dpos(2) /= 0) THEN
          nvec((ovec(2))) = vec(dpos(2))
          nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(2))
       ELSE
          nvec(ovec(2)) = 1
       END IF
       IF (dpos(3) /= 0) THEN
          nvec((ovec(3))) = vec(dpos(3))
          nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(3))
       ELSE
          nvec(ovec(3)) = 1
       END IF
       SELECT CASE(vtype)
       CASE(VTYPE_DOUBLE)
          dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
               = var%dat%vd(i)
       CASE(VTYPE_REAL)
          dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
               = var%dat%vr(i)
       CASE(VTYPE_INT)
          dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
               = REAL(var%dat%vi(i),dp)
       CASE(VTYPE_BYTE)
          dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
               = REAL(var%dat%vb(i),dp)
       CASE(VTYPE_CHAR)
          CALL RGMSG(substr, RGMLE, &
               'DATA OF TYPE CHAR NOT SUPPORTED !')
       CASE(VTYPE_UNDEF)
          CALL RGMSG(substr, RGMLE, &
               'UNDEFINED TYPE OF DATA !')
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, &
               'UNRECOGNIZED TYPE OF DATA !')
       END SELECT
       !
       DEALLOCATE(vec, STAT=status)
       NULLIFY(vec)
       CALL ERRMSG(substr,status,2)
    END DO

  END SUBROUTINE RGTOOL_CONVERT
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_CONVERT_DAT2VAR(var, dat, vname, grid, order)

    ! CONVERTS 4D-ARRAY (dat) TO TYPE N-ARRAY (var) ON GRID
    ! THE OPTIONAL STRING order DEFINES THE ORDER OF DIMENSIONS
    ! (DEFAULT: 'xyzn')
    !
    ! Author: Andreas Baumgaertner, MPIC, July 2010

    ! REGRID MODULES
    USE MESSY_NCREGRID_BASE
    USE MESSY_NCREGRID_NETCDF
    USE MESSY_NCREGRID_GEOHYB

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, INDEX, PRESENT, PRODUCT, SIZE, SUM, TRIM

    ! I/O
    TYPE (ncvar),      INTENT(OUT)          :: var   ! nc-variable
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: dat   ! DATA ON DESTINATION GRID
    CHARACTER(LEN=*), INTENT(IN)            :: vname     ! name of variable
    TYPE (geohybgrid), INTENT(IN)           :: grid  ! grid information
    CHARACTER(LEN=4),  INTENT(IN), OPTIONAL :: order ! DEFAULT: 'xyzn'

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_CONVERT_DAT2VAR'
    INTEGER :: ovec(4)          ! order of (x,y,z,n) for input field
!!$    INTEGER :: ldvar(4)         ! ordered dim-lenghts in var
    INTEGER :: dpos(3)          ! position of lon, lat, lev in variable
    INTEGER :: nrvdim           ! number of recognized dimensions in variable
    INTEGER :: npdlv            ! product of dimension-lengths in variable
!!$    CHARACTER(LEN=6) :: ostr(4) ! order-string for output
    INTEGER :: i,j                ! counter
    INTEGER :: nvec(4)          ! counter parameter for dimensions
    INTEGER, DIMENSION(:), POINTER :: vec   ! element vector
    INTEGER, DIMENSION(:), ALLOCATABLE :: zdimlen  ! local dim. lengths
    INTEGER :: status           ! status flag

    ! INIT
    CALL INIT_NCVAR(var)
    CALL RENAME_NCVAR(var, vname)
    var%xtype = NF90_DOUBLE ! xtype ?
    ! CREATE VARIABLE DIMENSIONS FROM GRID 
    nrvdim = 0
    IF (LMODE2DH) THEN
       IF (QDEF_NCVAR(grid%lonm)) nrvdim = nrvdim + 1
       IF (QDEF_NCVAR(grid%latm)) nrvdim = nrvdim + 1
    ELSE
       IF (QDEF_NCVAR(grid%col)) nrvdim = nrvdim + 1
    END IF
    IF (QDEF_NCVAR(grid%hyam) .OR. QDEF_NCVAR(grid%hybm)) &
         nrvdim = nrvdim + 1
    IF (QDEF_NCVAR(grid%timem)) nrvdim = nrvdim + 1

    var%ndims = nrvdim
    ALLOCATE(var%dim(var%ndims))
    
    NULLIFY(vec)
    ovec = (/ 1,2,3,4 /)                               ! DEFAULT
!!$    IF (LMODE2DH) THEN
!!$       ostr = (/ 'x/LON ', 'y/LAT ', 'z/LEV ', 'n/PAR '/) ! DEFAULT
!!$    ELSE
!!$       ostr = (/ 'x/COL ', 'y/DUM ', 'z/LEV ', 'n/PAR '/) ! DEFAULT
!!$    END IF
    IF (PRESENT(order)) THEN
       IF (LMODE2DH) THEN
          ovec(1) = INDEX(order, 'x')
!!$          ostr(ovec(1)) = 'x/LON '
          ovec(2) = INDEX(order, 'y')
!!$          ostr(ovec(2)) = 'y/LAT '
       ELSE
          ovec(1) = INDEX(order, 'x')
!!$          ostr(ovec(1)) = 'x/COL '
          ovec(2) = INDEX(order, 'y')
!!$          ostr(ovec(2)) = 'y/DUM '
       END IF
       ovec(3) = INDEX(order, 'z')
!!$       ostr(ovec(3)) = 'z/LEV '
       ovec(4) = INDEX(order, 'n')
!!$       ostr(ovec(4)) = 'n/PAR '
       IF (SUM(ovec) /= 10) THEN
          CALL RGMSG(substr, RGMLE, &
               'ERROR IN '''//TRIM(order)//'''-STRING AT SUBROUTINE CALL !')
       END IF
    END IF
    !
    !
    dpos(:) = 0
    nrvdim = 0
 
    ! LOOP OVER VARIABLE DIMENSIONS AND CHECK NAMES, LENGTHS, AND POSITIONS
    IF (LMODE2DH) THEN
       IF (ASSOCIATED(grid%lonm%dim)) THEN
          nrvdim = nrvdim + 1
          CALL COPY_NCDIM(var%dim(nrvdim), grid%lonm%dim(1))
!!$          ldvar(ovec(1)) = var%dim(nrvdim)%len
          dpos(1) = 1 ! BUT THIS SHOULD BE DETERMINED BY ORDER
       END IF
       IF (ASSOCIATED(grid%latm%dim)) THEN
          nrvdim = nrvdim + 1
          CALL COPY_NCDIM(var%dim(nrvdim), grid%latm%dim(1))
!!$          ldvar(ovec(2)) = var%dim(nrvdim)%len
          dpos(2) = 2
       END IF
    ELSE
       IF (ASSOCIATED(grid%col%dim)) THEN
          nrvdim = nrvdim + 1
          CALL COPY_NCDIM(var%dim(nrvdim), grid%col%dim(1))
!!$          ldvar(ovec(1)) = var%dim(nrvdim)%len
          dpos(1) = 1 ! BUT THIS SHOULD BE DETERMINED BY ORDER
       END IF
    END IF

    IF (ASSOCIATED(grid%hyam%dim)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(var%dim(nrvdim), grid%hyam%dim(1))
!!$       ldvar(ovec(3)) = var%dim(nrvdim)%len
       dpos(3) = 3
    ELSE
       IF (ASSOCIATED(grid%hybm%dim)) THEN
          nrvdim = nrvdim + 1
          CALL COPY_NCDIM(var%dim(nrvdim), grid%hyam%dim(1))
!!$          ldvar(ovec(3)) = var%dim(nrvdim)%len
          dpos(3) = 3
       END IF
    END IF
    IF (ASSOCIATED(grid%timem%dim)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(var%dim(nrvdim), grid%timem%dim(1))
    END IF

    var%ustep = grid%t
       

    ! PASS INFO AND DATA TO NARRAY
    IF (var%ndims > 0) THEN
       ALLOCATE(zdimlen(var%ndims), STAT=status)
       CALL ERRMSG(substr,status,13)
    ELSE
       ALLOCATE(zdimlen(1), STAT=status)
       CALL ERRMSG(substr,status,14)
       zdimlen(1) = 1
    END IF
    DO i=1, var%ndims
       zdimlen(i) = var%dim(i)%len  ! = 1 for UNLIM-DIM-ID
    END DO
    j = SIZE(zdimlen)
    CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_DOUBLE)

    IF (ASSOCIATED(var%dat%dim)) THEN
       npdlv = PRODUCT(var%dat%dim)
    ELSE
       RETURN
    END IF

    ! COPY DATA (LOOP OVER ELEMENTS)
    DO i=1, npdlv
       CALL ELEMENT(var%dat%dim, i, vec)
       nvec((ovec(4))) = PRODUCT(vec)
       IF (dpos(1) /= 0) THEN
          nvec((ovec(1))) = vec(dpos(1))
          nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(1))
       ELSE
          nvec(ovec(1)) = 1
       END IF
       IF (dpos(2) /= 0) THEN
          nvec((ovec(2))) = vec(dpos(2))
          nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(2))
       ELSE
          nvec(ovec(2)) = 1
       END IF
       IF (dpos(3) /= 0) THEN
          nvec((ovec(3))) = vec(dpos(3))
          nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(3))
       ELSE
          nvec(ovec(3)) = 1
       END IF
       var%dat%vd(i) = &
            dat(nvec(1), nvec(2), nvec(3), nvec(4)) 
       !
       DEALLOCATE(vec, STAT=status)
       NULLIFY(vec)
       CALL ERRMSG(substr,status,2)
    END DO

  END SUBROUTINE RGTOOL_CONVERT_DAT2VAR
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_READ_NCVAR(iou, nmlfile, vname, t, var   &
                                ,grid, lrg, lrgx, lrgy, lrgz, lok)

    ! PERFORMS ONE RE-GRIDDING (TIME-)STEP (t) FOR ONE VARIABLE
    ! (vname) AND RETURNS THE RESULT OF TYPE NCVAR (var)
    ! AND OPTIONALLY THE GRID (grid)
    ! LOGICAL I/O UNIT  : iou          (for namelist-file)
    ! NAMELIST-FILE     : nmlfile
    ! DIMENSION-SWITCHES: lrgx, lrgy, lrgz
    ! REGRIDDING-SWITCH : lrg
    ! SUCCESS-FLAG      : lok
    !
    ! Author: Patrick Joeckel, MPICH, October 2002

    ! NCREGRID INTERFACE
    USE messy_ncregrid_control,   ONLY: RG_CTRL, RG_NML, RG_STATUS  &
                                      , RG_SCAN, RG_PROC, RG_STOP   &
                                      , NML_NEXT,NML_STAY           &
                                      , RGSTAT_STOP                 &
                                      , REGRID_CONTROL
    USE messy_ncregrid_netcdf,    ONLY: ncvar, copy_ncvar, init_ncvar
    USE messy_ncregrid_geohyb,    ONLY: geohybgrid, copy_geohybgrid
    ! USE messy_ncregrid_diag

    IMPLICIT NONE

    INTRINSIC :: PRESENT, SIZE, TRIM

    ! I/O
    INTEGER, INTENT (IN)               :: iou          ! logical I/O unit
    CHARACTER(LEN=*), INTENT(IN)       :: nmlfile      ! namelist-file
    CHARACTER(LEN=*), INTENT(IN)       :: vname        ! variable name
    INTEGER, INTENT (IN)               :: t            ! netCDF time step
    TYPE(ncvar)                        :: var          ! nc-variable
    TYPE (geohybgrid), OPTIONAL        :: grid         ! output grid info
    LOGICAL, INTENT(IN), OPTIONAL      :: lrg          ! regrid really ?
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgx         ! regrid in x
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgy         ! regrid in y
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgz         ! regrid in z
    LOGICAL, INTENT (OUT), OPTIONAL    :: lok          ! OK?

    ! LOCAL
    ! OUTPUT OF REGRIDDING PROCEDURE
    TYPE (ncvar), DIMENSION(:), POINTER :: rvar    ! list of variables
    TYPE (geohybgrid)                   :: rgrid   ! output grid info
    INTEGER                             :: tt      ! time step
    LOGICAL                             :: llrg    ! regrid really ?
    INTEGER                             :: pcnt    ! procedure counter
    INTEGER                             :: i       ! counter

    ! INITIALIZE
    IF (PRESENT(lok)) THEN
       lok = .FALSE.
    END IF
    IF (PRESENT(lrg)) THEN
       llrg = lrg
    ELSE
       llrg = .true.  ! DEFAULT
    END IF
    !
    NULLIFY(rvar)
    !
    ! USE 1st TIME STEP FOR SCANNING, SINCE A NAMELIST FILE
    ! WITH MORE THAN ONE NAMELIST MUST NOT NECESSARILY
    ! CONTAIN TIME AXES OF EQUAL LENGTH ...
    tt = 1
    !
    ! COUNT PROCEDURE LEVEL: 0 : START, VARIABLE NOT YET FOUND IN NAMELIST FILE
    !                        1 : VARIABLE FOUND IN NAMELSIT FILE, BUT NOT YET
    !                            PROCESSED
    !                        2 : VARIABLE FOUND IN NAMELIST FILE, PROCESSING
    !                            ALREADY PERFORMED
    pcnt = 0

    ! START REGRIDDING
    ! REGRIDDING CONTROLLED BY INPUT VARIABELS
    RG_CTRL = RG_SCAN   ! SCAN-MODE, NO REGRIDDING
    RG_NML  = NML_NEXT  ! READ NEXT NAMELIST FROM FILE
    !
    DO ! endless DO loop (must be terminated with EXIT)

       CALL REGRID_CONTROL(RG_CTRL, RG_NML, RG_STATUS           &
                           , TRIM(nmlfile)                      &
                           , iounit = iou                       &
                           , lrgx=lrgx, lrgy=lrgy, lrgz=lrgz    &
                           , tmin=tt, tmax=tt, tstep=1, tret=tt &
                           , var = rvar                         &
                           , grid = rgrid                       &
                          )
       IF (RG_STATUS == RGSTAT_STOP) EXIT

       DO i=1, SIZE(rvar)
          IF (TRIM(rvar(i)%name) == TRIM(vname)) THEN
             ! VARIABLE FOUND IN NAMELIST:
             ! SET PROCEDURE LEVEL
             pcnt = pcnt + 1
             EXIT
          END IF
       END DO

       ! VARIABLE NOT YET FOUND IN NAMELIST FILE:
       ! SCAN NEXT NAMELIST IN NAMELIST FILE
       IF (pcnt == 0) THEN
          RG_NML  = NML_NEXT   ! ... TRY NEXT NAMELIST
          RG_CTRL = RG_SCAN    ! ... SCAN FOR VARIABLE
          tt = 1               ! ... USING 1st TIME STEP
       END IF

       ! VARIABLE FOUND IN NAMELIST FILE:
       ! SET SWITCHES FOR PROCESSING
       IF (pcnt == 1) THEN    ! IMMEDIATELY AFTER FILE HAS BEEN FOUND ...
          RG_NML  = NML_STAY  ! ... DO NOT READ NEXT NAMELIST FROM FILE
          tt = t              ! ... PROCESS REQUESTET TIME STEP
          IF (llrg) THEN
             RG_CTRL = RG_PROC  ! ... PERFORM REGRIDDING
          ELSE
             RG_CTRL = RG_SCAN  ! ... SCAN AGAIN (= IMPORT RAW DATA)
          END IF
       END IF

       ! FILENAME FOUND IN NAMELIST FILE AND PROCESSING FINISHED:
       ! SET SWITCHES TO END LOOP
       IF (pcnt == 2) THEN
          RG_CTRL = RG_STOP              ! ... TERMINATE ...
                                         ! ... REGRIDDING PROCEDURE CORRECTLY
          CALL COPY_NCVAR(var, rvar(i))  ! ... RETURN VARIABLE
          IF (PRESENT(grid)) THEN
             CALL COPY_GEOHYBGRID(grid, rgrid)   ! ... RETURN GRID
          END IF
          IF (PRESENT(lok)) THEN
             lok = .TRUE.                        ! ... REPORT SUCCESS
          END IF
       END IF

       DO i=1, SIZE(rvar)
          CALL INIT_NCVAR(rvar(i))
       END DO
       DEALLOCATE(rvar)
       NULLIFY(rvar)

    END DO  ! ENDLESS DO-LOOP ...
    ! END REGRIDDING

  END SUBROUTINE RGTOOL_READ_NCVAR
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_READ_NCFILE(iou, nmlfile, fname, t, var        &
                                ,grid, lrg, lrgx, lrgy, lrgz, lok)

    ! PERFORMS ONE RE-GRIDDING (TIME-)STEP (t) FOR ONE FILE
    ! (fname) AND RETURNS LIST OF TYPE NCVAR (var)
    ! AND OPTIONALLY THE GRID (grid)
    ! LOGICAL I/O UNIT  : iou          (for namelist-file)
    ! NAMELIST-FILE     : nmlfile
    ! DIMENSION-SWITCHES: lrgx, lrgy, lrgz
    ! REGRIDDING-SWITCH : lrg
    ! SUCCESS-FLAG      : lok
    !
    ! IF fname = '' THEN PROCESS 1st NAMELIST ONLY !!!
    !
    ! Author: Patrick Joeckel, MPICH, October 2002

    ! NCREGRID INTERFACE
    USE messy_ncregrid_control,   ONLY: RG_CTRL, RG_NML, RG_STATUS  &
                                      , RG_SCAN, RG_PROC, RG_STOP   &
                                      , NML_NEXT,NML_STAY           &
                                      , RGSTAT_STOP                 &
                                      , REGRID_CONTROL
    USE messy_ncregrid_netcdf,    ONLY: ncvar, copy_ncvar, init_ncvar
    USE messy_ncregrid_geohyb,    ONLY: geohybgrid, copy_geohybgrid
    USE messy_ncregrid_base,      ONLY: errmsg
    ! USE messy_ncregrid_diag

    IMPLICIT NONE

    INTRINSIC :: PRESENT, SIZE, TRIM

    ! I/O
    INTEGER, INTENT (IN)               :: iou          ! logical I/O unit
    CHARACTER(LEN=*), INTENT(IN)       :: nmlfile      ! namelist-file
    CHARACTER(LEN=*), INTENT(IN)       :: fname        ! filename
    INTEGER, INTENT (IN)               :: t            ! netCDF time step
    TYPE(ncvar), DIMENSION(:), POINTER :: var          ! nc-variables
    TYPE (geohybgrid), OPTIONAL        :: grid         ! output grid info
    LOGICAL, INTENT(IN), OPTIONAL      :: lrg          ! regrid really ?
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgx         ! regrid in x
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgy         ! regrid in y
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgz         ! regrid in z
    LOGICAL, INTENT (OUT), OPTIONAL    :: lok          ! OK?

    ! LOCAL
    ! OUTPUT OF REGRIDDING PROCEDURE
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_READ_NCFILE'
    TYPE (ncvar), DIMENSION(:), POINTER :: rvar    ! list of variables
    TYPE (geohybgrid)                   :: rgrid   ! output grid info
    INTEGER                             :: status  ! memory status
    INTEGER                             :: tt      ! time step
    LOGICAL                             :: llrg    ! regrid really ?
    INTEGER                             :: pcnt    ! procedure counter
    INTEGER                             :: i       ! counter

    ! INITIALIZE
    IF (PRESENT(lok)) THEN
       lok = .FALSE.
    END IF
    IF (PRESENT(lrg)) THEN
       llrg = lrg
    ELSE
       llrg = .true.  ! DEFAULT
    END IF
    !
    NULLIFY(rvar)
    !
    ! USE 1st TIME STEP FOR SCANNING, SINCE A NAMELIST FILE
    ! WITH MORE THAN ONE NAMELIST MUST NOT NECESSARILY
    ! CONTAIN TIME AXES OF EQUAL LENGTH ...
    tt = 1
    !
    ! COUNT PROCEDURE LEVEL: 0 : START, FILENAME NOT YET FOUND IN NAMELIST FILE
    !                        1 : FILENAME FOUND IN NAMELSIT FILE, BUT NOT YET
    !                            PROCESSED
    !                        2 : FILENAME FOUND IN NAMELIST FILE, PROCESSING
    !                            ALREADY PERFORMED
    pcnt = 0

    ! START REGRIDDING
    ! REGRIDDING CONTROLLED BY INPUT VARIABELS
    RG_CTRL = RG_SCAN   ! SCAN-MODE, NO REGRIDDING
    RG_NML  = NML_NEXT  ! READ NEXT NAMELIST FROM FILE
    !
    DO ! endless DO loop (must be terminated with EXIT)

       CALL REGRID_CONTROL(RG_CTRL, RG_NML, RG_STATUS           &
                           , TRIM(nmlfile)                      &
                           , iounit = iou                       &
                           , lrgx=lrgx, lrgy=lrgy, lrgz=lrgz    &
                           , tmin=tt, tmax=tt, tstep=1, tret=tt &
                           , var = rvar                         &
                           , grid = rgrid                       &
                          )
       IF (RG_STATUS == RGSTAT_STOP) EXIT

       ! FILENAME FOUND IN NAMELIST OR ALREADY PROCESSED
       IF ((TRIM(rgrid%file) == TRIM(fname)).OR.(pcnt == 1) &
            .OR. (TRIM(fname) == '') ) THEN
          ! SET PROCEDURE LEVEL
          pcnt = pcnt + 1
       END IF

       ! FILENAME NOT YET FOUND IN NAMELIST FILE:
       ! SCAN NEXT NAMELIST IN NAMELIST FILE
       IF (pcnt == 0) THEN
          RG_NML  = NML_NEXT   ! ... TRY NEXT NAMELIST
          RG_CTRL = RG_SCAN    ! ... SCAN FOR VARIABLE
          tt = 1               ! ... USING 1st TIME STEP
       END IF

       ! FILENAME FOUND IN NAMELIST FILE:
       ! SET SWITCHES FOR PROCESSING
       IF (pcnt == 1) THEN    ! IMMEDIATELY AFTER FILE HAS BEEN FOUND ...
          RG_NML  = NML_STAY  ! ... DO NOT READ NEXT NAMELIST FROM FILE
          tt = t              ! ... PROCESS REQUESTET TIME STEP
          IF (llrg) THEN
             RG_CTRL = RG_PROC  ! ... PERFORM REGRIDDING
          ELSE
             RG_CTRL = RG_SCAN  ! ... SCAN AGAIN (= IMPORT RAW DATA)
          END IF
       END IF

       ! FILENAME FOUND IN NAMELIST FILE AND PROCESSING FINISHED:
       ! SET SWITCHES TO END LOOP
       IF (pcnt == 2) THEN
          RG_CTRL = RG_STOP              ! ... TERMINATE ...
                                         ! ... REGRIDDING PROCEDURE CORRECTLY
          ALLOCATE(var(SIZE(rvar)), STAT=status)
          CALL ERRMSG(substr,status,1)
          DO i=1, SIZE(rvar)
             CALL COPY_NCVAR(var(i), rvar(i))    ! ... RETURN VARIABLES
          END DO
          IF (PRESENT(grid)) THEN
             CALL COPY_GEOHYBGRID(grid, rgrid)   ! ... RETURN GRID
          END IF
          IF (PRESENT(lok)) THEN
             lok = .TRUE.                        ! ... REPORT SUCCESS
          END IF
       END IF

       DO i=1, SIZE(rvar)
          CALL INIT_NCVAR(rvar(i))
       END DO
       DEALLOCATE(rvar)
       NULLIFY(rvar)

    END DO  ! ENDLESS DO-LOOP ...
    ! END REGRIDDING

  END SUBROUTINE RGTOOL_READ_NCFILE
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_G2C(grid, hyam, hybm, p0, ps         &
                        ,hyai, hybi                      &
                        ,latm, lonm                      &
                        ,lati, loni, nlat, nlon, nlev)

    ! CONVERTS GRID INFORMATION (grid) TO ARRAYS
    ! hyam, hybm, hyai, hybi, latm, lati, lonm, loni, p0 -> 1D
    ! ps                                                 -> 2D
    !
    ! nlat, nlon, nlev are OPTIONAL parameters in order to specify
    ! the array-length (box-mid), in case the grid does not contain
    ! the respective dimension
    !
    ! Author: Patrick Joeckel, MPICH, October 2002

    USE messy_ncregrid_geohyb,    ONLY: geohybgrid
    USE messy_ncregrid_base,      ONLY: errmsg

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT, SIZE

    ! I/O
    TYPE (geohybgrid), INTENT(IN)           :: grid
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: hyam   ! hybrid-A-coeff.
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: hybm   ! hybrid-B-coeff.
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: p0     ! reference pressure
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL :: ps     ! surface pressure
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: hyai   ! hybrid-A-coeff.
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: hybi   ! hybrid-B-coeff.
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: latm   ! latitude
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: lonm   ! longitude
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: lati   ! latitude
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: loni   ! longitude
    INTEGER,           INTENT(IN), OPTIONAL :: nlat
    INTEGER,           INTENT(IN), OPTIONAL :: nlon
    INTEGER,           INTENT(IN), OPTIONAL :: nlev

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_G2C'
    CHARACTER(LEN=4), PARAMETER :: order='xzny'
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: hha, hp0, hhb, hps, hlat, hlon
    INTEGER                     :: status

    NULLIFY(hha)
    NULLIFY(hp0)
    NULLIFY(hhb)
    NULLIFY(hps)
    NULLIFY(hlat)
    NULLIFY(hlon)

    IF (PRESENT(hyam)) THEN
       CALL RGTOOL_CONVERT(grid%hyam,hha,grid, order=order)
       IF (ASSOCIATED(hha)) THEN
          ALLOCATE(hyam(SIZE(hha,2)), STAT=status)
          CALL ERRMSG(substr,status,1)
          hyam(:) = hha(1,:,1,1)
          DEALLOCATE(hha, STAT=status)
          NULLIFY(hha)
          CALL ERRMSG(substr,status,2)
       ELSE
          IF (PRESENT(nlev)) THEN
             ALLOCATE(hyam(nlev), STAT=status)
             CALL ERRMSG(substr,status,3)
             hyam(:) = 0.0
          ELSE
             NULLIFY(hyam)
          END IF
       END IF
    END IF

    IF (PRESENT(hyai)) THEN
       CALL RGTOOL_CONVERT(grid%hyai,hha,grid, order=order)
       IF (ASSOCIATED(hha)) THEN
          ALLOCATE(hyai(SIZE(hha,2)), STAT=status)
          CALL ERRMSG(substr,status,1)
          hyai(:) = hha(1,:,1,1)
          DEALLOCATE(hha, STAT=status)
          NULLIFY(hha)
          CALL ERRMSG(substr,status,2)
       ELSE
          IF (PRESENT(nlev)) THEN
             ALLOCATE(hyai(nlev+1), STAT=status)
             CALL ERRMSG(substr,status,3)
             hyai(:) = 0.0
          ELSE
             NULLIFY(hyai)
          END IF
       END IF
    END IF

    IF (PRESENT(p0)) THEN
       CALL RGTOOL_CONVERT(grid%p0, hp0, grid, order=order)
       IF (ASSOCIATED(hp0)) THEN
          ALLOCATE(p0(1), STAT=status)
          CALL ERRMSG(substr,status,4)
          p0 = hp0(1,1,1,1)
          DEALLOCATE(hp0, STAT=status)
          NULLIFY(hp0)
          CALL ERRMSG(substr,status,5)
       ELSE
          ALLOCATE(p0(1), STAT=status)
          CALL ERRMSG(substr,status,6)
          p0(1) = 0.0
       ENDIF
    END IF

    IF (PRESENT(hybm)) THEN
       CALL RGTOOL_CONVERT(grid%hybm, hhb, grid, order=order)
       IF (ASSOCIATED(hhb)) THEN
          ALLOCATE(hybm(SIZE(hhb,2)), STAT=status)
          CALL ERRMSG(substr,status,7)
          hybm(:) = hhb(1,:,1,1)
          DEALLOCATE(hhb, STAT=status)
          NULLIFY(hhb)
          CALL ERRMSG(substr,status,8)
       ELSE
          IF (PRESENT(nlev)) THEN
             ALLOCATE(hybm(nlev), STAT=status)
             CALL ERRMSG(substr,status,9)
             hybm(:) = 0.0
          ELSE
             NULLIFY(hybm)
          END IF
       END IF
    END IF

    IF (PRESENT(hybi)) THEN
       CALL RGTOOL_CONVERT(grid%hybi, hhb, grid, order=order)
       IF (ASSOCIATED(hhb)) THEN
          ALLOCATE(hybi(SIZE(hhb,2)), STAT=status)
          CALL ERRMSG(substr,status,7)
          hybi(:) = hhb(1,:,1,1)
          DEALLOCATE(hhb, STAT=status)
          NULLIFY(hhb)
          CALL ERRMSG(substr,status,8)
       ELSE
          IF (PRESENT(nlev)) THEN
             ALLOCATE(hybi(nlev+1), STAT=status)
             CALL ERRMSG(substr,status,9)
             hybi(:) = 0.0
          ELSE
             NULLIFY(hybi)
          END IF
       END IF
    END IF

    IF (PRESENT(ps)) THEN
       CALL RGTOOL_CONVERT(grid%ps, hps, grid, order=order)
       IF (ASSOCIATED(hps)) THEN
          ALLOCATE(ps(SIZE(hps,1),SIZE(hps,4)), STAT=status)
          CALL ERRMSG(substr,status,10)
          ps(:,:) = hps(:,1,1,:)
          DEALLOCATE(hps, STAT=status)
          NULLIFY(hps)
          CALL ERRMSG(substr,status,11)
       ELSE
          IF (PRESENT(nlat).AND.PRESENT(nlon)) THEN
             ALLOCATE(ps(nlon,nlat), STAT=status)
             CALL ERRMSG(substr,status,12)
             ps(:,:) = 0.0
          ELSE
             NULLIFY(ps)
          END IF
       END IF
    END IF

    IF (LMODE2DH) THEN
       IF (PRESENT(latm)) THEN
          CALL RGTOOL_CONVERT(grid%latm, hlat, grid, order=order)
          IF (ASSOCIATED(hlat)) THEN
             ALLOCATE(latm(SIZE(hlat,4)), STAT=status)
             CALL ERRMSG(substr,status,13)
             latm(:) = hlat(1,1,1,:)
             DEALLOCATE(hlat, STAT=status)
             NULLIFY(hlat)
             CALL ERRMSG(substr,status,14)
          ELSE
             IF (PRESENT(nlat)) THEN
                ALLOCATE(latm(nlat), STAT=status)
                CALL ERRMSG(substr,status,15)
                latm(:) = 0.0
             ELSE
                NULLIFY(latm)
             END IF
          END IF
       END IF

       IF (PRESENT(lati)) THEN
          CALL RGTOOL_CONVERT(grid%lati, hlat, grid, order=order)
          IF (ASSOCIATED(hlat)) THEN
             ALLOCATE(lati(SIZE(hlat,4)), STAT=status)
             CALL ERRMSG(substr,status,13)
             lati(:) = hlat(1,1,1,:)
             DEALLOCATE(hlat, STAT=status)
             NULLIFY(hlat)
             CALL ERRMSG(substr,status,14)
          ELSE
             IF (PRESENT(nlat)) THEN
                ALLOCATE(lati(nlat+1), STAT=status)
                CALL ERRMSG(substr,status,15)
                lati(:) = 0.0
             ELSE
                NULLIFY(lati)
             END IF
          END IF
       END IF

       IF (PRESENT(lonm)) THEN
          CALL RGTOOL_CONVERT(grid%lonm, hlon, grid, order=order)
          IF (ASSOCIATED(hlon)) THEN
             ALLOCATE(lonm(SIZE(hlon,1)), STAT=status)
             CALL ERRMSG(substr,status,16)
             lonm(:) = hlon(:,1,1,1)
             DEALLOCATE(hlon, STAT=status)
             NULLIFY(hlon)
             CALL ERRMSG(substr,status,17)
          ELSE
             IF (PRESENT(nlon)) THEN
                ALLOCATE(lonm(nlon), STAT=status)
                CALL ERRMSG(substr,status,18)
                lonm(:) = 0.0
             ELSE
                NULLIFY(lonm)
             END IF
          END IF
       END IF

       IF (PRESENT(loni)) THEN
          CALL RGTOOL_CONVERT(grid%loni, hlon, grid, order=order)
          IF (ASSOCIATED(hlon)) THEN
             ALLOCATE(loni(SIZE(hlon,1)), STAT=status)
             CALL ERRMSG(substr,status,16)
             loni(:) = hlon(:,1,1,1)
             DEALLOCATE(hlon, STAT=status)
             NULLIFY(hlon)
             CALL ERRMSG(substr,status,17)
          ELSE
             IF (PRESENT(nlon)) THEN
                ALLOCATE(loni(nlon+1), STAT=status)
                CALL ERRMSG(substr,status,18)
                loni(:) = 0.0
             ELSE
                NULLIFY(loni)
             END IF
          END IF
       END IF
    END IF

  END SUBROUTINE RGTOOL_G2C
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE CLEAN_NCRGCNT_LIST(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_ncrgcnt_list), POINTER :: ai => NULL()
    TYPE(t_ncrgcnt_list), POINTER :: ae => NULL()

    ! INIT
    status = 1
    
    ai => GRGTLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       !
       ae => ai
       ai => ai%next
       !
       DEALLOCATE(ae)
       NULLIFY(ae)
       !
    END DO    

    NULLIFY(GRGTLIST)

    status = 0

  END SUBROUTINE CLEAN_NCRGCNT_LIST
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE WRITE_NCRGCNT_LIST(status, fname, iou)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: fname
    INTEGER,          INTENT(IN)  :: iou

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr = 'WRITE_NCRGCNT_LIST'
    TYPE(t_ncrgcnt_list), POINTER :: ai => NULL()
    TYPE(t_ncrgcnt_list), POINTER :: ae => NULL()

    ! EMPTY LIST
    IF (.NOT. ASSOCIATED (GRGTLIST)) THEN
       status = 0
       RETURN
    END IF

    WRITE(*,*) substr//&
         &': WRITING COUNTER LIST TO FILE '''//TRIM(fname)//''' ...'

    OPEN(iou, file=TRIM(fname), status='UNKNOWN')
    
    ai => GRGTLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       !
       ae => ai
       ai => ai%next
       !
       WRITE(*,*) '--> ',TRIM(ae%mname)//' @ '//TRIM(ae%this%name)
       WRITE(*,*) '    ', ae%this%start, &
            ae%this%step, ae%this%reset, ae%this%current
       !
       WRITE(iou, fstr) ae%mname, ae%this%name, ae%this%start, &
            ae%this%step, ae%this%reset, ae%this%current
       !
    END DO    

    CLOSE(iou)

    WRITE(*,*) substr//&
         &': ... END OF COUNTER LIST REACHED !'

    status = 0

  END SUBROUTINE WRITE_NCRGCNT_LIST
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE READ_NCRGCNT_LIST(status, fname, iou)

    USE messy_ncregrid_base, ONLY: RGMSG, RGMLE, RGMLEC, RGMLI, RGMLIC

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: fname
    INTEGER,          INTENT(IN)  :: iou

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'READ_NCRGCNT_LIST'
    LOGICAL                  :: lex
    INTEGER                  :: fstat
    TYPE(ncrgcnt), POINTER   :: cntptr => NULL()
    CHARACTER(LEN=NCCRSTRL)  :: mname
    TYPE(ncrgcnt)            :: ch
    
    ! EMPTY LIST
    IF (.NOT. ASSOCIATED (GRGTLIST)) THEN
       status = 0
       RETURN
    END IF

    CALL RGMSG(substr, RGMLI, &
         'OPENING COUNTER FILE '''//TRIM(fname)//''' ...')

    INQUIRE(file=TRIM(fname), exist=lex)
    IF (.NOT.lex) THEN  ! ERROR
       CALL RGMSG(substr, RGMLE, &
            'COUNTER FILE NOT FOUND: '//TRIM(fname) )
       status = 5020 ! COUNTER FILE NOT FOUND
       RETURN
    END IF

    OPEN(iou, file=TRIM(fname), status='OLD')

    DO
       READ(iou,fstr, IOSTAT=fstat) &
            mname, ch%name,ch%start,ch%step,ch%reset,ch%current
       IF (fstat /= 0) EXIT
       !
       CALL RGMSG(substr, RGMLIC, &
            '... LOOKING FOR COUNTER '''//TRIM(mname)//' @ '&
            &//TRIM(ch%name)//''' ...')
       !
       CALL loc_ncrgcnt(status, TRIM(mname), TRIM(ch%name), cntptr)
       IF ((status /= 0) .AND. (status /= 5005)) THEN
          CALL ncrgcnt_halt(substr, status)
       END IF
       !
       IF (status == 0) THEN
          CALL RGMSG(substr, RGMLIC, '... UPDATING COUNTER:')
          ! UPDATE
          cntptr%start   = ch%start
          cntptr%step    = ch%step
          cntptr%reset   = ch%reset
          cntptr%current = ch%current
          CALL RGMSG(substr, RGMLIC, &
               '         START  : ',cntptr%start,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         STEP   : ',cntptr%step,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         RESET  : ',cntptr%reset,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         CURRENT: ',cntptr%current,' ')
       ELSE
          ! COUNTER DOES NOT EXIST
          CALL RGMSG(substr, RGMLIC, '... COUNTER NOT FOUND')
       END IF
       !
    END DO

    IF (fstat < 0) THEN
       ! OK: END OF FILE REACHED
       CALL RGMSG(substr, RGMLIC, &
            'END OF FILE '''//TRIM(fname)//''' REACHED !')
    ELSE
       ! ERROR
       CLOSE(iou)
       CALL RGMSG(substr, RGMLE, &
            'READ ERROR READING FROM FILE '''//&
            &TRIM(fname)//''' !' , .false.)
       CALL RGMSG(substr, RGMLEC, &
            'STATUS: ',fstat,' ')
    END IF

    CLOSE(iou)

    status = 0

  END SUBROUTINE READ_NCRGCNT_LIST
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_NCRGCNT_RST(mname, start, restart, event, c, lout, linit)

    ! MANAGES I/O OF COUNTER INFORMATION AT START AND RESTART
    !
    ! Author: Patrick Joeckel, MPICH, September 2005

    USE messy_ncregrid_base, ONLY: RGMSG, RGMLI, RGMLIC

    IMPLICIT NONE

    INTRINSIC :: PRESENT, TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)    :: mname
    LOGICAL,          INTENT(IN)    :: start   ! .true. at first time step
    LOGICAL,          INTENT(IN)    :: restart ! .true. at first time step of
                                               ! rerun
    LOGICAL,          INTENT(IN)    :: event   ! .true. on event
    TYPE(NCRGCNT),    INTENT(INOUT) :: c       ! counter - struct
    LOGICAL,          INTENT(IN)           :: lout
    LOGICAL,          INTENT(IN), OPTIONAL :: linit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_NCRGCNT_RST'
    TYPE(NCRGCNT),    POINTER   :: cntptr
    LOGICAL                     :: lf
    INTEGER                     :: status
    LOGICAL                     :: zlinit

    ! DO NOTHING IN CASE ...
    IF (.NOT.(start.OR.restart.OR.event)) THEN
       IF (lout) CALL RGMSG(substr, RGMLI, &
            'RGTEVENT '''//TRIM(mname)//' @ '//TRIM(c%name)//''' NOT ACTIVE !')
       RETURN
    ELSE
       IF (lout) CALL RGMSG(substr, RGMLI, &
            'RGTEVENT '''//TRIM(mname)//' @ '//TRIM(c%name)//''' ACTIVE !')
    END IF

    ! INIT
    IF (PRESENT(linit)) THEN
       zlinit = linit
    ELSE
       zlinit = .FALSE. ! DEFAULT
    END IF

    ! LOCATE COUNTER IN LIST
    IF (lout) CALL RGMSG(substr, RGMLIC, &
         '... LOOKING FOR COUNTER '''//TRIM(mname)//' @ '//&
         &TRIM(c%name)//''' ...')
    CALL loc_ncrgcnt(status, TRIM(mname), TRIM(c%name), cntptr)
    IF ((status /= 0) .AND. (status /= 5005)) THEN
       IF (lout) CALL ncrgcnt_halt(substr, status)
    END IF
    lf = (status == 0) ! FOUND ?
    IF (lf) THEN
       ! UPDATE FROM LIST
       IF (lout) CALL RGMSG(substr, RGMLIC, '... COPYING COUNTER FROM LIST:')
       c%start   = cntptr%start
       c%step    = cntptr%step
       c%reset   = cntptr%reset
       c%current = cntptr%current
       IF (lout) THEN
          CALL RGMSG(substr, RGMLIC, &
               '         START  : ',c%start,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         STEP   : ',c%step,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         RESET  : ',c%reset,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         CURRENT: ',c%current,' ')
       END IF
    END IF

    ! START
    IF (zlinit) THEN
       IF (lf) THEN
          IF (lout) CALL ncrgcnt_halt(substr, 5004) ! COUNTER EXISTS
       ELSE
          IF (lout) CALL RGMSG(substr, RGMLIC, &
               '... SAVING COUNTER '''//TRIM(mname)//' @ '//TRIM(c%name)&
               &//''' IN LIST (INIT)')
          ! SAVE IN LIST
          CALL new_ncrgcnt(status, TRIM(mname), c)
          IF (lout) CALL ncrgcnt_halt(substr, status)
          ! RE-SET POINTER
          CALL loc_ncrgcnt(status, TRIM(mname), TRIM(c%name), cntptr)
          IF (lout) CALL ncrgcnt_halt(substr, status)
          lf = .TRUE. ! NOW IT EXISTS
          !
       END IF
    END IF

    ! EVENT
    IF (event.AND.(.NOT.(start))) THEN
       IF (.NOT. lf) THEN
          IF (lout) CALL ncrgcnt_halt(substr, 5005) ! COUNTER DOES NOT EXIST
       ELSE
          IF (lout) CALL RGMSG(substr, RGMLIC, &
               '... UPDATING COUNTER '''//TRIM(mname)//' @ '//TRIM(c%name)&
               &//''' IN LIST (EVENT)')
          !
          ! UPDATE
          c%current = c%current + c%step                   ! INCREMENT
          IF (c%current > c%reset) c%current = c%start     ! RESET
          !
          ! SAVE IN LIST
          cntptr%start   = c%start
          cntptr%step    = c%step
          cntptr%reset   = c%reset
          cntptr%current = c%current
          !
          IF (lout)  THEN
             CALL RGMSG(substr, RGMLIC, &
                  '         START  : ',cntptr%start,' ')
             CALL RGMSG(substr, RGMLIC, &
                  '         STEP   : ',cntptr%step,' ')
             CALL RGMSG(substr, RGMLIC, &
                  '         RESET  : ',cntptr%reset,' ')
             CALL RGMSG(substr, RGMLIC, &
                  '         CURRENT: ',cntptr%current,' ')
          END IF
       END IF
    END IF

  END SUBROUTINE RGTOOL_NCRGCNT_RST
! ----------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE GET_NEXT_NCRGCNT(last, cntptr)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    LOGICAL,              INTENT(OUT) :: last
    TYPE(ncrgcnt),        POINTER     :: cntptr

    ! LOCAL
    TYPE(t_ncrgcnt_list), POINTER, SAVE :: ai  => NULL()
    TYPE(t_ncrgcnt_list), POINTER, SAVE :: ae  => NULL()
    INTEGER, PARAMETER                  :: MODE_INIT = 1
    INTEGER, PARAMETER                  :: MODE_CONT = 2
    INTEGER,                       SAVE :: MODE = MODE_INIT

    ! INIT
    last = .FALSE.

    DO
       SELECT CASE (MODE)
       CASE(MODE_INIT)
          ai => GRGTLIST
          MODE = MODE_CONT
       CASE(MODE_CONT)
          IF (.NOT. ASSOCIATED(ai)) EXIT
          ae => ai
          ai => ai%next
          cntptr => ae%this
          RETURN
       END SELECT
    END DO

    MODE = MODE_INIT
    last = .TRUE.

  END SUBROUTINE GET_NEXT_NCRGCNT
  ! -------------------------------------------------------------------

! ----------------------------------------------------------------------
! PRIVATE SUBROUTINES
! ----------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_ncrgcnt(status, mname, cname, cntptr)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, LEN_TRIM, TRIM, ASSOCIATED

    ! I/O
    INTEGER,              INTENT(OUT) :: status
    CHARACTER(LEN=*),     INTENT(IN)  :: mname
    CHARACTER(LEN=*),     INTENT(IN)  :: cname
    TYPE(ncrgcnt),        POINTER     :: cntptr

    ! LOCAL
    TYPE(t_ncrgcnt_list), POINTER :: ai  => NULL()
    TYPE(t_ncrgcnt_list), POINTER :: ae  => NULL()
    LOGICAL                       :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(cntptr)

    ! CHECKS
    IF (LEN_TRIM(ADJUSTL(mname)) > NCCRSTRL) THEN
       status = 5003   ! NAME TOO LONG
       RETURN
    END IF
    !
    IF (LEN_TRIM(ADJUSTL(cname)) > NCCRSTRL) THEN
       status = 5006   ! NAME TOO LONG
       RETURN
    END IF

    ! CHECK, IF IT EXISTS
    ai => GRGTLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF ( (TRIM(ADJUSTL(mname)) == TRIM(ai%mname)) .AND. &
            (TRIM(ADJUSTL(cname)) == TRIM(ai%this%name)) ) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO    

    IF (lex) THEN
       cntptr => ai%this
    ELSE
       status = 5005   ! DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_ncrgcnt
  ! -------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE new_ncrgcnt(status, mname, cnt)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, LEN_TRIM, ASSOCIATED, TRIM, NULL

    ! I/O
    INTEGER,          INTENT(OUT)  :: status
    CHARACTER(LEN=*), INTENT(IN)   :: mname
    TYPE(NCRGCNT),    INTENT(IN)   :: cnt

    ! LOCAL
    TYPE(t_ncrgcnt_list), POINTER :: ai     => NULL()
    TYPE(t_ncrgcnt_list), POINTER :: ae     => NULL()
    TYPE(ncrgcnt),        POINTER :: cntptr => NULL()
    INTEGER                       :: zstat

    IF (LEN_TRIM(ADJUSTL(mname)) > NCCRSTRL) THEN
       status = 5003   ! SUBMODEL NAME TOO LONG
       RETURN
    END IF

    CALL loc_ncrgcnt(zstat, TRIM(mname), TRIM(cnt%name), cntptr)
    IF (zstat /= 5005) THEN  ! COUNTER DOES NOT EXIST (IS OK HERE !)
       IF (zstat == 0) THEN  ! COUNTER EXISTS ALREADY
          status = 5004      ! COUNTER EXISTS ALREADY
       ELSE
          status = zstat     ! ERROR
       END IF
       RETURN
    END IF

    ! GOTO END OF LIST
    ai => GRGTLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO

    ! ADD NEW
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(GRGTLIST)) THEN
       GRGTLIST => ai                 ! SET POINTER TO FIRST ELEMENT
    ELSE
       ae%next => ai                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF
       
    ! SET VALUES
    ai%mname        = TRIM(ADJUSTL(mname))
    !
    ai%this%name    = TRIM(ADJUSTL(cnt%name))
    ai%this%start   = cnt%start
    ai%this%step    = cnt%step
    ai%this%reset   = cnt%reset
    ai%this%current = cnt%current

    status = 0

  END SUBROUTINE new_ncrgcnt
  ! ----------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE ncrgcnt_halt(substr, status)

    ! MESSy
    USE messy_main_constants_mem,  ONLY: STRLEN_VLONG
    USE messy_ncregrid_base,       ONLY: RGMSG, RGMLE
    
    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status
    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr
    
    IF (status == 0) RETURN
    
    errstr = ncrgcnt_error_str(status)
    
    CALL RGMSG(substr, RGMLE, TRIM(errstr))
    
  END SUBROUTINE ncrgcnt_halt
  ! -------------------------------------------------------------------
  
  ! -------------------------------------------------------------------
  FUNCTION ncrgcnt_error_str(status)
    
    USE messy_main_tools,         ONLY: int2str
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG
    
    IMPLICIT NONE
    
    ! I/O
    CHARACTER(LEN=STRLEN_VLONG)           :: ncrgcnt_error_str
    INTEGER,                   INTENT(IN) :: status
    
    ! LOCAL
    CHARACTER(LEN=4) :: echar = '    '
    
    ncrgcnt_error_str = ''
    
    CALL int2str(echar, status, cpad='0', cerr='*')
    
    SELECT CASE(echar)
       ! NO ERROR
    CASE('0000')
       ncrgcnt_error_str = 'E'//echar//': NO ERROR'

    CASE('5003')
       ncrgcnt_error_str = 'E'//echar//': SUBMODEL NAME TOO LONG'
    CASE('5004')
       ncrgcnt_error_str = 'E'//echar//': COUNTER EXISTS ALREADY'
    CASE('5005')
       ncrgcnt_error_str = 'E'//echar//': COUNTER DOES NOT EXIST'
    CASE('5006')
       ncrgcnt_error_str = 'E'//echar//': COUNTER NAME TOO LONG'

    CASE('5020')
       ncrgcnt_error_str = 'E'//echar//': COUNTER FILE NOT FOUND'

    CASE DEFAULT
       ncrgcnt_error_str = 'E'//echar//': UNKONW ERROR STATUS'

    END SELECT

  END FUNCTION ncrgcnt_error_str
  ! -------------------------------------------------------------------

! ******************************************************************
! ------------------------------------------------------------------
END MODULE MESSY_NCREGRID_TOOLS
! ------------------------------------------------------------------
! ******************************************************************

#else

! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_TOOLS
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, October 2002
! ******************************************************************

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  INTEGER,           PARAMETER   :: NCCNTMAXNLEN = 40                  !   ! 
  CHARACTER(LEN=24), PARAMETER   :: fstr = '(a80,1x,a40,1x,4(i8,1x))'  ! <-V
  INTEGER,           PARAMETER   :: NCCRSTRL = 80                      !   |

  ! TYPE DECLARATION TO HOLD COUNTER INFORMATION
  TYPE NCRGCNT
     CHARACTER(LEN=NCCNTMAXNLEN) :: name = ''
     INTEGER                     :: start = 1
     INTEGER                     :: step  = 1
     INTEGER                     :: reset = 1
     INTEGER                     :: current = 1
  END TYPE NCRGCNT
  
  ! LIST TO STORE COUNTER INFORMATION FOR RESTART
  TYPE T_NCRGCNT_LIST
     CHARACTER(NCCRSTRL)           :: mname  = ''
     TYPE(NCRGCNT)                 :: this
     TYPE(T_NCRGCNT_LIST), POINTER :: next => NULL()
  END TYPE T_NCRGCNT_LIST

  ! ====================================================================
  ! CONCAT. LIST OF RGT-EVENTS
  TYPE(T_NCRGCNT_LIST), POINTER, SAVE :: GRGTLIST => NULL()
  ! ====================================================================
  
  PUBLIC  :: RGTOOL_CONVERT     ! CONVERTS NCVAR INTO 4D-ARRAY
  PUBLIC  :: RGTOOL_READ_NCVAR  ! NCREGRID ONE STEP OF ONE FIELD
  PUBLIC  :: RGTOOL_READ_NCFILE ! NCREGRID ONE STEP OF ONE FILE
  PUBLIC  :: RGTOOL_G2C         ! CONVERTS GRID-INFORMATION INTO ARRAYS

  ! COUNTER INFORMATION HANDLING
  PUBLIC  :: NCRGCNT            ! TYPE STRUCT FOR COUNTER INFORMATION I/O
  PUBLIC  :: RGTOOL_NCRGCNT_RST ! MANAGES COUNTER INFORMATION RESTART I/O
  !
  PUBLIC  :: CLEAN_NCRGCNT_LIST
  PUBLIC  :: WRITE_NCRGCNT_LIST
  PUBLIC  :: READ_NCRGCNT_LIST
  PUBLIC  :: GET_NEXT_NCRGCNT
  !
  !PRIVATE :: loc_ncrgcnt
  !PRIVATE :: new_ncrgcnt
  !PRIVATE :: ncrgcnt_error_str

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_CONVERT(var, dat, grid, order)

    ! CONVERTS NCREGRID OUTPUT OF TYPE N-ARRAY (var) ON GRID
    ! grid TO 4D-ARRAY (dat)
    ! THE OPTIONAL STRING order DEFINES THE ORDER OF DIMENSIONS
    ! (DEFAULT: 'xyzn')
    !
    ! Author: Patrick Joeckel, MPIC, October 2002

    ! REGRID MODULES
    USE MESSY_NCREGRID_BASE
    USE MESSY_NCREGRID_NETCDF
    USE MESSY_NCREGRID_GEOHYB

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, INDEX, PRESENT, PRODUCT, REAL, SIZE, SUM, TRIM

    ! I/O
    TYPE (ncvar),      INTENT(IN)           :: var   ! nc-variable
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: dat   ! DATA ON DESTINATION GRID
    TYPE (geohybgrid), INTENT(IN)           :: grid  ! grid information
    CHARACTER(LEN=4),  INTENT(IN), OPTIONAL :: order ! DEFAULT: 'xyzn'

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_CONVERT'
    INTEGER :: ovec(4)          ! order of (x,y,z,n) for output field
    INTEGER :: ldvar(4)         ! ordered dim-lenghts in var
    INTEGER :: dpos(3)          ! position of lon, lat, lev in variable
    INTEGER :: nrvdim           ! number of recognized dimensions in variable
    INTEGER :: npdlv            ! product of dimension-lengths in variable
    CHARACTER(LEN=6) :: ostr(4) ! order-string for output
    INTEGER :: i                ! counter
    INTEGER :: nvec(4)          ! counter parameter for dimensions
    INTEGER :: vtype            ! type of variable
    INTEGER, DIMENSION(:), POINTER :: vec   ! element vector
    INTEGER :: status           ! status flag

    INTEGER, DIMENSION(:,:), POINTER      :: vechelp  ! element vector of PS
    INTEGER                               :: mhelp
    INTEGER , DIMENSION(:,:), ALLOCATABLE :: daccelementhelp
    ! counter parameter for dimensions
    INTEGER , DIMENSION(:,:), ALLOCATABLE :: nvec_help
    NULLIFY(vechelp)

    ! INIT
    NULLIFY(vec)
    ovec = (/ 1,2,3,4 /)                               ! DEFAULT
    ostr = (/ 'x/LON ', 'y/LAT ', 'z/LEV ', 'n/PAR '/) ! DEFAULT
    IF (PRESENT(order)) THEN
       ovec(1) = INDEX(order, 'x')
       ostr(ovec(1)) = 'x/LON '
       ovec(2) = INDEX(order, 'y')
       ostr(ovec(2)) = 'y/LAT '
       ovec(3) = INDEX(order, 'z')
       ostr(ovec(3)) = 'z/LEV '
       ovec(4) = INDEX(order, 'n')
       ostr(ovec(4)) = 'n/PAR '
       IF (SUM(ovec) /= 10) THEN
          CALL RGMSG(substr, RGMLE, &
               'ERROR IN '''//TRIM(order)//'''-STRING AT SUBROUTINE CALL !')
       END IF
    END IF
    !
    IF (ASSOCIATED(dat)) THEN
       DEALLOCATE(dat, STAT=status)
       CALL ERRMSG(substr,status,1)
    END IF
    NULLIFY(dat)
    !
    dpos(:) = 0
    nrvdim = 0
    ldvar(:) = 1
    IF (ASSOCIATED(var%dat%dim)) THEN
       npdlv = PRODUCT(var%dat%dim)
    ELSE
       RETURN
    END IF

    ! LOOP OVER VARIABLE DIMENSIONS AND CHECK NAMES, LENGTHS, AND POSITIONS
    DO i=1, var%ndims
       IF (ASSOCIATED(grid%lonm%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%lonm%dim(1)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(1)) = var%dim(i)%len
             dpos(1) = i
          END IF
       END IF
       IF (ASSOCIATED(grid%latm%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%latm%dim(1)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(2)) = var%dim(i)%len
             dpos(2) = i
          END IF
       END IF
       IF (ASSOCIATED(grid%hyam%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%hyam%dim(1)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(3)) = var%dim(i)%len
             dpos(3) = i
          END IF
       ELSE
          IF (ASSOCIATED(grid%hybm%dim)) THEN
             IF (QCMP_NCDIM(var%dim(i), grid%hybm%dim(1)) > 1) THEN
                nrvdim = nrvdim+1
                ldvar(ovec(3)) = var%dim(i)%len
                dpos(3) = i
             END IF
          END IF
       END IF
    END DO

    ldvar(ovec(4)) = npdlv/(ldvar(ovec(1))*ldvar(ovec(2))*ldvar(ovec(3)))

    CALL RGMSG(substr, RGMLI, &
         'VARIABLE '''//TRIM(var%name)//''' IS ',var%ndims, &
         '-DIMENSIONAL')
    CALL RGMSG(substr, RGMLIC, '(',var%dat%dim,')')
    CALL RGMSG(substr, RGMLIC, &
         'DATA DIMENSION-LENGHTS ARE: ')
    DO i=1, SIZE(ostr)
       CALL RGMSG(substr, RGMLIC, '  '//ostr(i)//': ',ldvar(i),' ')
    END DO
    CALL RGMSG(substr, RGMLIC, &
         'SOURCE POSITIONS (x/LON, y/LAT, z/LEV) ARE: ',dpos,' ')

    ! ALLOCATE DATA SPACE
    ALLOCATE(dat(ldvar(1),ldvar(2),ldvar(3),ldvar(4)), STAT=status)
    CALL ERRMSG(substr,status,1)

    ! COPY DATA (LOOP OVER ELEMENTS)
    vtype = var%dat%type
    ALLOCATE(nvec_help(npdlv,4), STAT=status)

    SELECT CASE(SIZE(var%dat%dim))

    CASE(1) !SIZE(dim)=1

      ALLOCATE(vechelp(npdlv,1), STAT=status)

      SELECT CASE(vtype)
      CASE(VTYPE_DOUBLE)

        DO i=1, npdlv
           vechelp(i,1) = i
  
           nvec_help(i,(ovec(4))) = vechelp(i,1)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
              dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                   , nvec_help(i,4)) &
                   = var%dat%vd(i)
        END DO

      CASE(VTYPE_REAL)

        DO i=1, npdlv
           vechelp(i,1) = i
  
           nvec_help(i,(ovec(4))) = vechelp(i,1)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
              dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                   , nvec_help(i,4)) &
                   = var%dat%vr(i)
        END DO

      CASE(VTYPE_INT)

        DO i=1, npdlv
           vechelp(i,1) = i
  
           nvec_help(i,(ovec(4))) = vechelp(i,1)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
              dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                   , nvec_help(i,4)) &
                   = REAL(var%dat%vi(i),dp)
        END DO

      CASE(VTYPE_BYTE)

        DO i=1, npdlv
           vechelp(i,1) = i
  
           nvec_help(i,(ovec(4))) = vechelp(i,1)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
              dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                   , nvec_help(i,4)) &
                   = REAL(var%dat%vb(i),dp)
        END DO

      END SELECT

      DEALLOCATE(vechelp, STAT=status)
      NULLIFY(vechelp)

    CASE(2) !SIZE(dim)=2

      ALLOCATE(daccelementhelp(npdlv,2))
      ALLOCATE(vechelp(npdlv,2), STAT=status)

      SELECT CASE(vtype)
      CASE(VTYPE_DOUBLE)

        DO i=1, npdlv
  
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp
  
           nvec_help(i,(ovec(4))) = vechelp(i,1)*vechelp(i,2)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
            dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                 , nvec_help(i,4)) &
                 = var%dat%vd(i)
        END DO

      CASE(VTYPE_REAL)

        DO i=1, npdlv
  
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp
  
           nvec_help(i,(ovec(4))) = vechelp(i,1)*vechelp(i,2)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
            dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                 , nvec_help(i,4)) &
                 = var%dat%vr(i)
        END DO

      CASE(VTYPE_INT)

        DO i=1, npdlv
  
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp
  
           nvec_help(i,(ovec(4))) = vechelp(i,1)*vechelp(i,2)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                , nvec_help(i,4)) &
                = REAL(var%dat%vi(i),dp)
        END DO

      CASE(VTYPE_BYTE)

        DO i=1, npdlv
  
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp
  
           nvec_help(i,(ovec(4))) = vechelp(i,1)*vechelp(i,2)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                , nvec_help(i,4)) &
                = REAL(var%dat%vb(i),dp)
        END DO

      END SELECT

      DEALLOCATE(daccelementhelp)
      DEALLOCATE(vechelp, STAT=status)
      NULLIFY(vechelp)

    CASE(3) !SIZE(dim)=3

      ALLOCATE(daccelementhelp(npdlv,3))
      ALLOCATE(vechelp(npdlv,3), STAT=status)

      SELECT CASE(vtype)
      CASE(VTYPE_DOUBLE)

        DO i=1, npdlv

           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*var%dat%dim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           nvec_help(i,(ovec(4))) = vechelp(i,1)*vechelp(i,2)*vechelp(i,3)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                , nvec_help(i,4)) &
                = var%dat%vd(i)
        END DO

      CASE(VTYPE_REAL)

        DO i=1, npdlv

           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*var%dat%dim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           nvec_help(i,(ovec(4))) = vechelp(i,1)*vechelp(i,2)*vechelp(i,3)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                , nvec_help(i,4)) &
                = var%dat%vr(i)
        END DO

      CASE(VTYPE_INT)

        DO i=1, npdlv

           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*var%dat%dim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           nvec_help(i,(ovec(4))) = vechelp(i,1)*vechelp(i,2)*vechelp(i,3)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                , nvec_help(i,4)) &
                = REAL(var%dat%vi(i),dp)
        END DO

      CASE(VTYPE_BYTE)

        DO i=1, npdlv

           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*var%dat%dim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           nvec_help(i,(ovec(4))) = vechelp(i,1)*vechelp(i,2)*vechelp(i,3)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                , nvec_help(i,4)) &
                = REAL(var%dat%vb(i),dp)
        END DO

      END SELECT

      DEALLOCATE(daccelementhelp)
      DEALLOCATE(vechelp, STAT=status)
      NULLIFY(vechelp)

    CASE(4) !SIZE(dim)=4

      ALLOCATE(daccelementhelp(npdlv,4))
      ALLOCATE(vechelp(npdlv,4), STAT=status)

      SELECT CASE(vtype)
      CASE(VTYPE_DOUBLE)

        DO i=1, npdlv
  
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           vechelp(i,4) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*var%dat%dim(2)
           daccelementhelp(i,4) = daccelementhelp(i,3)*var%dat%dim(3)
           vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
           mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           nvec_help(i,(ovec(4))) = &
                vechelp(i,1)*vechelp(i,2)*vechelp(i,3)*vechelp(i,4)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3) &
                , nvec_help(i,4)) &
                = var%dat%vd(i)
        END DO

      CASE(VTYPE_REAL)

        DO i=1, npdlv
  
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           vechelp(i,4) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*var%dat%dim(2)
           daccelementhelp(i,4) = daccelementhelp(i,3)*var%dat%dim(3)
           vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
           mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           nvec_help(i,(ovec(4))) = &
                vechelp(i,1)*vechelp(i,2)*vechelp(i,3)*vechelp(i,4)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3), nvec_help(i,4)) &
                = var%dat%vr(i)
        END DO

      CASE(VTYPE_INT)

        DO i=1, npdlv
  
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           vechelp(i,4) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*var%dat%dim(2)
           daccelementhelp(i,4) = daccelementhelp(i,3)*var%dat%dim(3)
           vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
           mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           nvec_help(i,(ovec(4))) = &
                vechelp(i,1)*vechelp(i,2)*vechelp(i,3)*vechelp(i,4)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3), nvec_help(i,4)) &
                = REAL(var%dat%vi(i),dp)
        END DO

      CASE(VTYPE_BYTE)

        DO i=1, npdlv
  
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           vechelp(i,4) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*var%dat%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*var%dat%dim(2)
           daccelementhelp(i,4) = daccelementhelp(i,3)*var%dat%dim(3)
           vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
           mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           nvec_help(i,(ovec(4))) = &
                vechelp(i,1)*vechelp(i,2)*vechelp(i,3)*vechelp(i,4)
           IF (dpos(1) /= 0) THEN
              nvec_help(i,(ovec(1))) = vechelp(i,dpos(1))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(1))
           ELSE
              nvec_help(i,ovec(1)) = 1
           END IF
           IF (dpos(2) /= 0) THEN
              nvec_help(i,(ovec(2))) = vechelp(i,dpos(2))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(2))
           ELSE
              nvec_help(i,ovec(2)) = 1
           END IF
           IF (dpos(3) /= 0) THEN
              nvec_help(i,(ovec(3))) = vechelp(i,dpos(3))
              nvec_help(i,ovec(4)) = nvec_help(i,ovec(4))/nvec_help(i,ovec(3))
           ELSE
              nvec_help(i,ovec(3)) = 1
           END IF
           dat(nvec_help(i,1), nvec_help(i,2), nvec_help(i,3), nvec_help(i,4)) &
                = REAL(var%dat%vb(i),dp)
        END DO

      END SELECT

      DEALLOCATE(daccelementhelp)
      DEALLOCATE(vechelp, STAT=status)
      NULLIFY(vechelp)


    CASE DEFAULT
      DO i=1, npdlv
         CALL ELEMENT(var%dat%dim, i, vec)
         nvec((ovec(4))) = PRODUCT(vec)
         IF (dpos(1) /= 0) THEN
            nvec((ovec(1))) = vec(dpos(1))
            nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(1))
         ELSE
            nvec(ovec(1)) = 1
         END IF
         IF (dpos(2) /= 0) THEN
            nvec((ovec(2))) = vec(dpos(2))
            nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(2))
         ELSE
            nvec(ovec(2)) = 1
         END IF
         IF (dpos(3) /= 0) THEN
            nvec((ovec(3))) = vec(dpos(3))
            nvec(ovec(4)) = nvec(ovec(4))/nvec(ovec(3))
         ELSE
            nvec(ovec(3)) = 1
         END IF
         SELECT CASE(vtype)
         CASE(VTYPE_DOUBLE)
            dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
                 = var%dat%vd(i)
         CASE(VTYPE_REAL)
            dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
                 = var%dat%vr(i)
         CASE(VTYPE_INT)
            dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
                 = REAL(var%dat%vi(i),dp)
         CASE(VTYPE_BYTE)
            dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
                 = REAL(var%dat%vb(i),dp)
         CASE(VTYPE_CHAR)
            CALL RGMSG(substr, RGMLE, &
                 'DATA OF TYPE CHAR NOT SUPPORTED !')
         CASE(VTYPE_UNDEF)
            CALL RGMSG(substr, RGMLE, &
                 'UNDEFINED TYPE OF DATA !')
         CASE DEFAULT
            CALL RGMSG(substr, RGMLE, &
                 'UNRECOGNIZED TYPE OF DATA !')
       END SELECT
       !
       DEALLOCATE(vec, STAT=status)
       NULLIFY(vec)
       CALL ERRMSG(substr,status,2)
      END DO
    END SELECT
    DEALLOCATE(nvec_help)
    NULLIFY(nvec_help)

  END SUBROUTINE RGTOOL_CONVERT
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_READ_NCVAR(iou, nmlfile, vname, t, var   &
                                ,grid, lrg, lrgx, lrgy, lrgz, lok)

    ! PERFORMS ONE RE-GRIDDING (TIME-)STEP (t) FOR ONE VARIABLE
    ! (vname) AND RETURNS THE RESULT OF TYPE NCVAR (var)
    ! AND OPTIONALLY THE GRID (grid)
    ! LOGICAL I/O UNIT  : iou          (for namelist-file)
    ! NAMELIST-FILE     : nmlfile
    ! DIMENSION-SWITCHES: lrgx, lrgy, lrgz
    ! REGRIDDING-SWITCH : lrg
    ! SUCCESS-FLAG      : lok
    !
    ! Author: Patrick Joeckel, MPICH, October 2002

    ! NCREGRID INTERFACE
    USE messy_ncregrid_control,   ONLY: RG_CTRL, RG_NML, RG_STATUS  &
                                      , RG_SCAN, RG_PROC, RG_STOP   &
                                      , NML_NEXT,NML_STAY           &
                                      , RGSTAT_STOP                 &
                                      , REGRID_CONTROL
    USE messy_ncregrid_netcdf,    ONLY: ncvar, copy_ncvar
    USE messy_ncregrid_geohyb,    ONLY: geohybgrid, copy_geohybgrid
    ! USE messy_ncregrid_diag

    IMPLICIT NONE

    INTRINSIC :: PRESENT, SIZE, TRIM

    ! I/O
    INTEGER, INTENT (IN)               :: iou          ! logical I/O unit
    CHARACTER(LEN=*), INTENT(IN)       :: nmlfile      ! namelist-file
    CHARACTER(LEN=*), INTENT(IN)       :: vname        ! variable name
    INTEGER, INTENT (IN)               :: t            ! netCDF time step
    TYPE(ncvar)                        :: var          ! nc-variable
    TYPE (geohybgrid), OPTIONAL        :: grid         ! output grid info
    LOGICAL, INTENT(IN), OPTIONAL      :: lrg          ! regrid really ?
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgx         ! regrid in x
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgy         ! regrid in y
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgz         ! regrid in z
    LOGICAL, INTENT (OUT), OPTIONAL    :: lok          ! OK?

    ! LOCAL
    ! OUTPUT OF REGRIDDING PROCEDURE
    TYPE (ncvar), DIMENSION(:), POINTER :: rvar    ! list of variables
    TYPE (geohybgrid)                   :: rgrid   ! output grid info
    INTEGER                             :: tt      ! time step
    LOGICAL                             :: llrg    ! regrid really ?
    INTEGER                             :: pcnt    ! procedure counter
    INTEGER                             :: i       ! counter

    ! INITIALIZE
    IF (PRESENT(lok)) THEN
       lok = .FALSE.
    END IF
    IF (PRESENT(lrg)) THEN
       llrg = lrg
    ELSE
       llrg = .true.  ! DEFAULT
    END IF
    !
    NULLIFY(rvar)
    !
    ! USE 1st TIME STEP FOR SCANNING, SINCE A NAMELIST FILE
    ! WITH MORE THAN ONE NAMELIST MUST NOT NECESSARILY
    ! CONTAIN TIME AXES OF EQUAL LENGTH ...
    tt = 1
    !
    ! COUNT PROCEDURE LEVEL: 0 : START, VARIABLE NOT YET FOUND IN NAMELIST FILE
    !                        1 : VARIABLE FOUND IN NAMELSIT FILE, BUT NOT YET
    !                            PROCESSED
    !                        2 : VARIABLE FOUND IN NAMELIST FILE, PROCESSING
    !                            ALREADY PERFORMED
    pcnt = 0

    ! START REGRIDDING
    ! REGRIDDING CONTROLLED BY INPUT VARIABELS
    RG_CTRL = RG_SCAN   ! SCAN-MODE, NO REGRIDDING
    RG_NML  = NML_NEXT  ! READ NEXT NAMELIST FROM FILE
    !
    DO ! endless DO loop (must be terminated with EXIT)

       CALL REGRID_CONTROL(RG_CTRL, RG_NML, RG_STATUS           &
                           , TRIM(nmlfile)                      &
                           , iounit = iou                       &
                           , lrgx=lrgx, lrgy=lrgy, lrgz=lrgz    &
                           , tmin=tt, tmax=tt, tstep=1, tret=tt &
                           , var = rvar                         &
                           , grid = rgrid                       &
                          )
       IF (RG_STATUS == RGSTAT_STOP) EXIT

       DO i=1, SIZE(rvar)
          IF (TRIM(rvar(i)%name) == TRIM(vname)) THEN
             ! VARIABLE FOUND IN NAMELIST:
             ! SET PROCEDURE LEVEL
             pcnt = pcnt + 1
             EXIT
          END IF
       END DO

       ! VARIABLE NOT YET FOUND IN NAMELIST FILE:
       ! SCAN NEXT NAMELIST IN NAMELIST FILE
       IF (pcnt == 0) THEN
          RG_NML  = NML_NEXT   ! ... TRY NEXT NAMELIST
          RG_CTRL = RG_SCAN    ! ... SCAN FOR VARIABLE
          tt = 1               ! ... USING 1st TIME STEP
       END IF

       ! VARIABLE FOUND IN NAMELIST FILE:
       ! SET SWITCHES FOR PROCESSING
       IF (pcnt == 1) THEN    ! IMMEDIATELY AFTER FILE HAS BEEN FOUND ...
          RG_NML  = NML_STAY  ! ... DO NOT READ NEXT NAMELIST FROM FILE
          tt = t              ! ... PROCESS REQUESTET TIME STEP
          IF (llrg) THEN
             RG_CTRL = RG_PROC  ! ... PERFORM REGRIDDING
          ELSE
             RG_CTRL = RG_SCAN  ! ... SCAN AGAIN (= IMPORT RAW DATA)
          END IF
       END IF

       ! FILENAME FOUND IN NAMELIST FILE AND PROCESSING FINISHED:
       ! SET SWITCHES TO END LOOP
       IF (pcnt == 2) THEN
          RG_CTRL = RG_STOP              ! ... TERMINATE ...
                                         ! ... REGRIDDING PROCEDURE CORRECTLY
          CALL COPY_NCVAR(var, rvar(i))  ! ... RETURN VARIABLE
          IF (PRESENT(grid)) THEN
             CALL COPY_GEOHYBGRID(grid, rgrid)   ! ... RETURN GRID
          END IF
          IF (PRESENT(lok)) THEN
             lok = .TRUE.                        ! ... REPORT SUCCESS
          END IF
       END IF

    END DO  ! ENDLESS DO-LOOP ...
    ! END REGRIDDING

  END SUBROUTINE RGTOOL_READ_NCVAR
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_READ_NCFILE(iou, nmlfile, fname, t, var        &
                                ,grid, lrg, lrgx, lrgy, lrgz, lok)

    ! PERFORMS ONE RE-GRIDDING (TIME-)STEP (t) FOR ONE FILE
    ! (fname) AND RETURNS LIST OF TYPE NCVAR (var)
    ! AND OPTIONALLY THE GRID (grid)
    ! LOGICAL I/O UNIT  : iou          (for namelist-file)
    ! NAMELIST-FILE     : nmlfile
    ! DIMENSION-SWITCHES: lrgx, lrgy, lrgz
    ! REGRIDDING-SWITCH : lrg
    ! SUCCESS-FLAG      : lok
    !
    ! IF fname = '' THEN PROCESS 1st NAMELIST ONLY !!!
    !
    ! Author: Patrick Joeckel, MPICH, October 2002

    ! NCREGRID INTERFACE
    USE messy_ncregrid_control,   ONLY: RG_CTRL, RG_NML, RG_STATUS  &
                                      , RG_SCAN, RG_PROC, RG_STOP   &
                                      , NML_NEXT,NML_STAY           &
                                      , RGSTAT_STOP                 &
                                      , REGRID_CONTROL
    USE messy_ncregrid_netcdf,    ONLY: ncvar, copy_ncvar
    USE messy_ncregrid_geohyb,    ONLY: geohybgrid, copy_geohybgrid
    USE messy_ncregrid_base,      ONLY: errmsg
    ! USE messy_ncregrid_diag

    IMPLICIT NONE

    INTRINSIC :: PRESENT, SIZE, TRIM

    ! I/O
    INTEGER, INTENT (IN)               :: iou          ! logical I/O unit
    CHARACTER(LEN=*), INTENT(IN)       :: nmlfile      ! namelist-file
    CHARACTER(LEN=*), INTENT(IN)       :: fname        ! filename
    INTEGER, INTENT (IN)               :: t            ! netCDF time step
    TYPE(ncvar), DIMENSION(:), POINTER :: var          ! nc-variables
    TYPE (geohybgrid), OPTIONAL        :: grid         ! output grid info
    LOGICAL, INTENT(IN), OPTIONAL      :: lrg          ! regrid really ?
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgx         ! regrid in x
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgy         ! regrid in y
    LOGICAL, INTENT(IN), OPTIONAL      :: lrgz         ! regrid in z
    LOGICAL, INTENT (OUT), OPTIONAL    :: lok          ! OK?

    ! LOCAL
    ! OUTPUT OF REGRIDDING PROCEDURE
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_READ_NCFILE'
    TYPE (ncvar), DIMENSION(:), POINTER :: rvar    ! list of variables
    TYPE (geohybgrid)                   :: rgrid   ! output grid info
    INTEGER                             :: status  ! memory status
    INTEGER                             :: tt      ! time step
    LOGICAL                             :: llrg    ! regrid really ?
    INTEGER                             :: pcnt    ! procedure counter
    INTEGER                             :: i       ! counter

    ! INITIALIZE
    IF (PRESENT(lok)) THEN
       lok = .FALSE.
    END IF
    IF (PRESENT(lrg)) THEN
       llrg = lrg
    ELSE
       llrg = .true.  ! DEFAULT
    END IF
    !
    NULLIFY(rvar)
    !
    ! USE 1st TIME STEP FOR SCANNING, SINCE A NAMELIST FILE
    ! WITH MORE THAN ONE NAMELIST MUST NOT NECESSARILY
    ! CONTAIN TIME AXES OF EQUAL LENGTH ...
    tt = 1
    !
    ! COUNT PROCEDURE LEVEL: 0 : START, FILENAME NOT YET FOUND IN NAMELIST FILE
    !                        1 : FILENAME FOUND IN NAMELSIT FILE, BUT NOT YET
    !                            PROCESSED
    !                        2 : FILENAME FOUND IN NAMELIST FILE, PROCESSING
    !                            ALREADY PERFORMED
    pcnt = 0

    ! START REGRIDDING
    ! REGRIDDING CONTROLLED BY INPUT VARIABELS
    RG_CTRL = RG_SCAN   ! SCAN-MODE, NO REGRIDDING
    RG_NML  = NML_NEXT  ! READ NEXT NAMELIST FROM FILE
    !
    DO ! endless DO loop (must be terminated with EXIT)

       CALL REGRID_CONTROL(RG_CTRL, RG_NML, RG_STATUS           &
                           , TRIM(nmlfile)                      &
                           , iounit = iou                       &
                           , lrgx=lrgx, lrgy=lrgy, lrgz=lrgz    &
                           , tmin=tt, tmax=tt, tstep=1, tret=tt &
                           , var = rvar                         &
                           , grid = rgrid                       &
                          )
       IF (RG_STATUS == RGSTAT_STOP) EXIT

       ! FILENAME FOUND IN NAMELIST OR ALREADY PROCESSED
       IF ((TRIM(rgrid%file) == TRIM(fname)).OR.(pcnt == 1) &
            .OR. (TRIM(fname) == '') ) THEN
          ! SET PROCEDURE LEVEL
          pcnt = pcnt + 1
       END IF

       ! FILENAME NOT YET FOUND IN NAMELIST FILE:
       ! SCAN NEXT NAMELIST IN NAMELIST FILE
       IF (pcnt == 0) THEN
          RG_NML  = NML_NEXT   ! ... TRY NEXT NAMELIST
          RG_CTRL = RG_SCAN    ! ... SCAN FOR VARIABLE
          tt = 1               ! ... USING 1st TIME STEP
       END IF

       ! FILENAME FOUND IN NAMELIST FILE:
       ! SET SWITCHES FOR PROCESSING
       IF (pcnt == 1) THEN    ! IMMEDIATELY AFTER FILE HAS BEEN FOUND ...
          RG_NML  = NML_STAY  ! ... DO NOT READ NEXT NAMELIST FROM FILE
          tt = t              ! ... PROCESS REQUESTET TIME STEP
          IF (llrg) THEN
             RG_CTRL = RG_PROC  ! ... PERFORM REGRIDDING
          ELSE
             RG_CTRL = RG_SCAN  ! ... SCAN AGAIN (= IMPORT RAW DATA)
          END IF
       END IF

       ! FILENAME FOUND IN NAMELIST FILE AND PROCESSING FINISHED:
       ! SET SWITCHES TO END LOOP
       IF (pcnt == 2) THEN
          RG_CTRL = RG_STOP              ! ... TERMINATE ...
                                         ! ... REGRIDDING PROCEDURE CORRECTLY
          ALLOCATE(var(SIZE(rvar)), STAT=status)
          CALL ERRMSG(substr,status,1)
          DO i=1, SIZE(rvar)
             CALL COPY_NCVAR(var(i), rvar(i))    ! ... RETURN VARIABLES
          END DO
          IF (PRESENT(grid)) THEN
             CALL COPY_GEOHYBGRID(grid, rgrid)   ! ... RETURN GRID
          END IF
          IF (PRESENT(lok)) THEN
             lok = .TRUE.                        ! ... REPORT SUCCESS
          END IF
       END IF

    END DO  ! ENDLESS DO-LOOP ...
    ! END REGRIDDING

  END SUBROUTINE RGTOOL_READ_NCFILE
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_G2C(grid, hyam, hybm, p0, ps         &
                        ,hyai, hybi                      &
                        ,latm, lonm                      &
                        ,lati, loni, nlat, nlon, nlev)

    ! CONVERTS GRID INFORMATION (grid) TO ARRAYS
    ! hyam, hybm, hyai, hybi, latm, lati, lonm, loni, p0 -> 1D
    ! ps                                                 -> 2D
    !
    ! nlat, nlon, nlev are OPTIONAL parameters in order to specify
    ! the array-length (box-mid), in case the grid does not contain
    ! the respective dimension
    !
    ! Author: Patrick Joeckel, MPICH, October 2002

    USE messy_ncregrid_geohyb,    ONLY: geohybgrid
    USE messy_ncregrid_base,      ONLY: errmsg

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT, SIZE

    ! I/O
    TYPE (geohybgrid), INTENT(IN)           :: grid
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: hyam   ! hybrid-A-coeff.
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: hybm   ! hybrid-B-coeff.
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: p0     ! reference pressure
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL :: ps     ! surface pressure
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: hyai   ! hybrid-A-coeff.
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: hybi   ! hybrid-B-coeff.
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: latm   ! latitude
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: lonm   ! longitude
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: lati   ! latitude
    REAL(dp), DIMENSION(:),   POINTER, OPTIONAL :: loni   ! longitude
    INTEGER,           INTENT(IN), OPTIONAL :: nlat
    INTEGER,           INTENT(IN), OPTIONAL :: nlon
    INTEGER,           INTENT(IN), OPTIONAL :: nlev

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_G2C'
    CHARACTER(LEN=4), PARAMETER :: order='xzny'
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: hha, hp0, hhb, hps, hlat, hlon
    INTEGER                     :: status

    NULLIFY(hha)
    NULLIFY(hp0)
    NULLIFY(hhb)
    NULLIFY(hps)
    NULLIFY(hlat)
    NULLIFY(hlon)

    IF (PRESENT(hyam)) THEN
       CALL RGTOOL_CONVERT(grid%hyam,hha,grid, order=order)
       IF (ASSOCIATED(hha)) THEN
          ALLOCATE(hyam(SIZE(hha,2)), STAT=status)
          CALL ERRMSG(substr,status,1)
          hyam(:) = hha(1,:,1,1)
          DEALLOCATE(hha, STAT=status)
          NULLIFY(hha)
          CALL ERRMSG(substr,status,2)
       ELSE
          IF (PRESENT(nlev)) THEN
             ALLOCATE(hyam(nlev), STAT=status)
             CALL ERRMSG(substr,status,3)
             hyam(:) = 0.0
          ELSE
             NULLIFY(hyam)
          END IF
       END IF
    END IF

    IF (PRESENT(hyai)) THEN
       CALL RGTOOL_CONVERT(grid%hyai,hha,grid, order=order)
       IF (ASSOCIATED(hha)) THEN
          ALLOCATE(hyai(SIZE(hha,2)), STAT=status)
          CALL ERRMSG(substr,status,1)
          hyai(:) = hha(1,:,1,1)
          DEALLOCATE(hha, STAT=status)
          NULLIFY(hha)
          CALL ERRMSG(substr,status,2)
       ELSE
          IF (PRESENT(nlev)) THEN
             ALLOCATE(hyai(nlev+1), STAT=status)
             CALL ERRMSG(substr,status,3)
             hyai(:) = 0.0
          ELSE
             NULLIFY(hyai)
          END IF
       END IF
    END IF

    IF (PRESENT(p0)) THEN
       CALL RGTOOL_CONVERT(grid%p0, hp0, grid, order=order)
       IF (ASSOCIATED(hp0)) THEN
          ALLOCATE(p0(1), STAT=status)
          CALL ERRMSG(substr,status,4)
          p0 = hp0(1,1,1,1)
          DEALLOCATE(hp0, STAT=status)
          NULLIFY(hp0)
          CALL ERRMSG(substr,status,5)
       ELSE
          ALLOCATE(p0(1), STAT=status)
          CALL ERRMSG(substr,status,6)
          p0(1) = 0.0
       ENDIF
    END IF

    IF (PRESENT(hybm)) THEN
       CALL RGTOOL_CONVERT(grid%hybm, hhb, grid, order=order)
       IF (ASSOCIATED(hhb)) THEN
          ALLOCATE(hybm(SIZE(hhb,2)), STAT=status)
          CALL ERRMSG(substr,status,7)
          hybm(:) = hhb(1,:,1,1)
          DEALLOCATE(hhb, STAT=status)
          NULLIFY(hhb)
          CALL ERRMSG(substr,status,8)
       ELSE
          IF (PRESENT(nlev)) THEN
             ALLOCATE(hybm(nlev), STAT=status)
             CALL ERRMSG(substr,status,9)
             hybm(:) = 0.0
          ELSE
             NULLIFY(hybm)
          END IF
       END IF
    END IF

    IF (PRESENT(hybi)) THEN
       CALL RGTOOL_CONVERT(grid%hybi, hhb, grid, order=order)
       IF (ASSOCIATED(hhb)) THEN
          ALLOCATE(hybi(SIZE(hhb,2)), STAT=status)
          CALL ERRMSG(substr,status,7)
          hybi(:) = hhb(1,:,1,1)
          DEALLOCATE(hhb, STAT=status)
          NULLIFY(hhb)
          CALL ERRMSG(substr,status,8)
       ELSE
          IF (PRESENT(nlev)) THEN
             ALLOCATE(hybi(nlev+1), STAT=status)
             CALL ERRMSG(substr,status,9)
             hybi(:) = 0.0
          ELSE
             NULLIFY(hybi)
          END IF
       END IF
    END IF

    IF (PRESENT(ps)) THEN
       CALL RGTOOL_CONVERT(grid%ps, hps, grid, order=order)
       IF (ASSOCIATED(hps)) THEN
          ALLOCATE(ps(SIZE(hps,1),SIZE(hps,4)), STAT=status)
          CALL ERRMSG(substr,status,10)
          ps(:,:) = hps(:,1,1,:)
          DEALLOCATE(hps, STAT=status)
          NULLIFY(hps)
          CALL ERRMSG(substr,status,11)
       ELSE
          IF (PRESENT(nlat).AND.PRESENT(nlon)) THEN
             ALLOCATE(ps(nlon,nlat), STAT=status)
             CALL ERRMSG(substr,status,12)
             ps(:,:) = 0.0
          ELSE
             NULLIFY(ps)
          END IF
       END IF
    END IF

    IF (PRESENT(latm)) THEN
       CALL RGTOOL_CONVERT(grid%latm, hlat, grid, order=order)
       IF (ASSOCIATED(hlat)) THEN
          ALLOCATE(latm(SIZE(hlat,4)), STAT=status)
          CALL ERRMSG(substr,status,13)
          latm(:) = hlat(1,1,1,:)
          DEALLOCATE(hlat, STAT=status)
          NULLIFY(hlat)
          CALL ERRMSG(substr,status,14)
       ELSE
          IF (PRESENT(nlat)) THEN
             ALLOCATE(latm(nlat), STAT=status)
             CALL ERRMSG(substr,status,15)
             latm(:) = 0.0
          ELSE
             NULLIFY(latm)
          END IF
       END IF
    END IF

    IF (PRESENT(lati)) THEN
       CALL RGTOOL_CONVERT(grid%lati, hlat, grid, order=order)
       IF (ASSOCIATED(hlat)) THEN
          ALLOCATE(lati(SIZE(hlat,4)), STAT=status)
          CALL ERRMSG(substr,status,13)
          lati(:) = hlat(1,1,1,:)
          DEALLOCATE(hlat, STAT=status)
          NULLIFY(hlat)
          CALL ERRMSG(substr,status,14)
       ELSE
          IF (PRESENT(nlat)) THEN
             ALLOCATE(lati(nlat+1), STAT=status)
             CALL ERRMSG(substr,status,15)
             lati(:) = 0.0
          ELSE
             NULLIFY(lati)
          END IF
       END IF
    END IF

    IF (PRESENT(lonm)) THEN
       CALL RGTOOL_CONVERT(grid%lonm, hlon, grid, order=order)
       IF (ASSOCIATED(hlon)) THEN
          ALLOCATE(lonm(SIZE(hlon,1)), STAT=status)
          CALL ERRMSG(substr,status,16)
          lonm(:) = hlon(:,1,1,1)
          DEALLOCATE(hlon, STAT=status)
          NULLIFY(hlon)
          CALL ERRMSG(substr,status,17)
       ELSE
          IF (PRESENT(nlon)) THEN
             ALLOCATE(lonm(nlon), STAT=status)
             CALL ERRMSG(substr,status,18)
             lonm(:) = 0.0
          ELSE
             NULLIFY(lonm)
          END IF
       END IF
    END IF

    IF (PRESENT(loni)) THEN
       CALL RGTOOL_CONVERT(grid%loni, hlon, grid, order=order)
       IF (ASSOCIATED(hlon)) THEN
          ALLOCATE(loni(SIZE(hlon,1)), STAT=status)
          CALL ERRMSG(substr,status,16)
          loni(:) = hlon(:,1,1,1)
          DEALLOCATE(hlon, STAT=status)
          NULLIFY(hlon)
          CALL ERRMSG(substr,status,17)
       ELSE
          IF (PRESENT(nlon)) THEN
             ALLOCATE(loni(nlon+1), STAT=status)
             CALL ERRMSG(substr,status,18)
             loni(:) = 0.0
          ELSE
             NULLIFY(loni)
          END IF
       END IF
    END IF

  END SUBROUTINE RGTOOL_G2C
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE CLEAN_NCRGCNT_LIST(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_ncrgcnt_list), POINTER :: ai => NULL()
    TYPE(t_ncrgcnt_list), POINTER :: ae => NULL()

    ! INIT
    status = 1
    
    ai => GRGTLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       !
       ae => ai
       ai => ai%next
       !
       DEALLOCATE(ae)
       NULLIFY(ae)
       !
    END DO    

    NULLIFY(GRGTLIST)

    status = 0

  END SUBROUTINE CLEAN_NCRGCNT_LIST
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE WRITE_NCRGCNT_LIST(status, fname, iou)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: fname
    INTEGER,          INTENT(IN)  :: iou

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr = 'WRITE_NCRGCNT_LIST'
    TYPE(t_ncrgcnt_list), POINTER :: ai => NULL()
    TYPE(t_ncrgcnt_list), POINTER :: ae => NULL()

    ! EMPTY LIST
    IF (.NOT. ASSOCIATED (GRGTLIST)) THEN
       status = 0
       RETURN
    END IF

    WRITE(*,*) substr//&
         &': WRITING COUNTER LIST TO FILE '''//TRIM(fname)//''' ...'

    OPEN(iou, file=TRIM(fname), status='UNKNOWN')
    
    ai => GRGTLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       !
       ae => ai
       ai => ai%next
       !
       WRITE(*,*) '--> ',TRIM(ae%mname)//' @ '//TRIM(ae%this%name)
       WRITE(*,*) '    ', ae%this%start, &
            ae%this%step, ae%this%reset, ae%this%current
       !
       WRITE(iou, fstr) ae%mname, ae%this%name, ae%this%start, &
            ae%this%step, ae%this%reset, ae%this%current
       !
    END DO    

    CLOSE(iou)

    WRITE(*,*) substr//&
         &': ... END OF COUNTER LIST REACHED !'

    status = 0

  END SUBROUTINE WRITE_NCRGCNT_LIST
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE READ_NCRGCNT_LIST(status, fname, iou)

    USE messy_ncregrid_base, ONLY: RGMSG, RGMLE, RGMLEC, RGMLI, RGMLIC

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: fname
    INTEGER,          INTENT(IN)  :: iou

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'READ_NCRGCNT_LIST'
    LOGICAL                  :: lex
    INTEGER                  :: fstat
    TYPE(ncrgcnt), POINTER   :: cntptr => NULL()
    CHARACTER(LEN=NCCRSTRL)  :: mname
    TYPE(ncrgcnt)            :: ch
    
    ! EMPTY LIST
    IF (.NOT. ASSOCIATED (GRGTLIST)) THEN
       status = 0
       RETURN
    END IF

    CALL RGMSG(substr, RGMLI, &
         'OPENING COUNTER FILE '''//TRIM(fname)//''' ...')

    INQUIRE(file=TRIM(fname), exist=lex)
    IF (.NOT.lex) THEN  ! ERROR
       CALL RGMSG(substr, RGMLE, &
            'COUNTER FILE NOT FOUND: '//TRIM(fname) )
       status = 5020 ! COUNTER FILE NOT FOUND
       RETURN
    END IF

    OPEN(iou, file=TRIM(fname), status='OLD')

    DO
       READ(iou,fstr, IOSTAT=fstat) &
            mname, ch%name,ch%start,ch%step,ch%reset,ch%current
       IF (fstat /= 0) EXIT
       !
       CALL RGMSG(substr, RGMLIC, &
            '... LOOKING FOR COUNTER '''//TRIM(mname)//' @ '&
            &//TRIM(ch%name)//''' ...')
       !
       CALL loc_ncrgcnt(status, TRIM(mname), TRIM(ch%name), cntptr)
       IF ((status /= 0) .AND. (status /= 5005)) THEN
          CALL ncrgcnt_halt(substr, status)
       END IF
       !
       IF (status == 0) THEN
          CALL RGMSG(substr, RGMLIC, '... UPDATING COUNTER:')
          ! UPDATE
          cntptr%start   = ch%start
          cntptr%step    = ch%step
          cntptr%reset   = ch%reset
          cntptr%current = ch%current
          CALL RGMSG(substr, RGMLIC, &
               '         START  : ',cntptr%start,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         STEP   : ',cntptr%step,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         RESET  : ',cntptr%reset,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         CURRENT: ',cntptr%current,' ')
       ELSE
          ! COUNTER DOES NOT EXIST
          CALL RGMSG(substr, RGMLIC, '... COUNTER NOT FOUND')
       END IF
       !
    END DO

    IF (fstat < 0) THEN
       ! OK: END OF FILE REACHED
       CALL RGMSG(substr, RGMLIC, &
            'END OF FILE '''//TRIM(fname)//''' REACHED !')
    ELSE
       ! ERROR
       CLOSE(iou)
       CALL RGMSG(substr, RGMLE, &
            'READ ERROR READING FROM FILE '''//&
            &TRIM(fname)//''' !' , .false.)
       CALL RGMSG(substr, RGMLEC, &
            'STATUS: ',fstat,' ')
    END IF

    CLOSE(iou)

    status = 0

  END SUBROUTINE READ_NCRGCNT_LIST
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_NCRGCNT_RST(mname, start, restart, event, c, lout, linit)

    ! MANAGES I/O OF COUNTER INFORMATION AT START AND RESTART
    !
    ! Author: Patrick Joeckel, MPICH, September 2005

    USE messy_ncregrid_base, ONLY: RGMSG, RGMLI, RGMLIC

    IMPLICIT NONE

    INTRINSIC :: PRESENT, TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)    :: mname
    LOGICAL,          INTENT(IN)    :: start   ! .true. at first time step
    LOGICAL,          INTENT(IN)    :: restart ! .true. at first time step of
                                               ! rerun
    LOGICAL,          INTENT(IN)    :: event   ! .true. on event
    TYPE(NCRGCNT),    INTENT(INOUT) :: c       ! counter - struct
    LOGICAL,          INTENT(IN)           :: lout
    LOGICAL,          INTENT(IN), OPTIONAL :: linit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_NCRGCNT_RST'
    TYPE(NCRGCNT),    POINTER   :: cntptr
    LOGICAL                     :: lf
    INTEGER                     :: status
    LOGICAL                     :: zlinit

    ! DO NOTHING IN CASE ...
    IF (.NOT.(start.OR.restart.OR.event)) THEN
       IF (lout) CALL RGMSG(substr, RGMLI, &
            'RGTEVENT '''//TRIM(mname)//' @ '//TRIM(c%name)//''' NOT ACTIVE !')
       RETURN
    ELSE
       IF (lout) CALL RGMSG(substr, RGMLI, &
            'RGTEVENT '''//TRIM(mname)//' @ '//TRIM(c%name)//''' ACTIVE !')
    END IF

    ! INIT
    IF (PRESENT(linit)) THEN
       zlinit = linit
    ELSE
       zlinit = .FALSE. ! DEFAULT
    END IF

    ! LOCATE COUNTER IN LIST
    IF (lout) CALL RGMSG(substr, RGMLIC, &
         '... LOOKING FOR COUNTER '''//TRIM(mname)//' @ '//&
         &TRIM(c%name)//''' ...')
    CALL loc_ncrgcnt(status, TRIM(mname), TRIM(c%name), cntptr)
    IF ((status /= 0) .AND. (status /= 5005)) THEN
       IF (lout) CALL ncrgcnt_halt(substr, status)
    END IF
    lf = (status == 0) ! FOUND ?
    IF (lf) THEN
       ! UPDATE FROM LIST
       IF (lout) CALL RGMSG(substr, RGMLIC, '... COPYING COUNTER FROM LIST:')
       c%start   = cntptr%start
       c%step    = cntptr%step
       c%reset   = cntptr%reset
       c%current = cntptr%current
       IF (lout) THEN
          CALL RGMSG(substr, RGMLIC, &
               '         START  : ',c%start,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         STEP   : ',c%step,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         RESET  : ',c%reset,' ')
          CALL RGMSG(substr, RGMLIC, &
               '         CURRENT: ',c%current,' ')
       END IF
    END IF

    ! START
    IF (zlinit) THEN
       IF (lf) THEN
          IF (lout) CALL ncrgcnt_halt(substr, 5004) ! COUNTER EXISTS
       ELSE
          IF (lout) CALL RGMSG(substr, RGMLIC, &
               '... SAVING COUNTER '''//TRIM(mname)//' @ '//TRIM(c%name)&
               &//''' IN LIST (INIT)')
          ! SAVE IN LIST
          CALL new_ncrgcnt(status, TRIM(mname), c)
          IF (lout) CALL ncrgcnt_halt(substr, status)
          ! RE-SET POINTER
          CALL loc_ncrgcnt(status, TRIM(mname), TRIM(c%name), cntptr)
          IF (lout) CALL ncrgcnt_halt(substr, status)
          lf = .TRUE. ! NOW IT EXISTS
          !
       END IF
    END IF

    ! EVENT
    IF (event.AND.(.NOT.(start))) THEN
       IF (.NOT. lf) THEN
          IF (lout) CALL ncrgcnt_halt(substr, 5005) ! COUNTER DOES NOT EXIST
       ELSE
          IF (lout) CALL RGMSG(substr, RGMLIC, &
               '... UPDATING COUNTER '''//TRIM(mname)//' @ '//TRIM(c%name)&
               &//''' IN LIST (EVENT)')
          !
          ! UPDATE
          c%current = c%current + c%step                   ! INCREMENT
          IF (c%current > c%reset) c%current = c%start     ! RESET
          !
          ! SAVE IN LIST
          cntptr%start   = c%start
          cntptr%step    = c%step
          cntptr%reset   = c%reset
          cntptr%current = c%current
          !
          IF (lout)  THEN
             CALL RGMSG(substr, RGMLIC, &
                  '         START  : ',cntptr%start,' ')
             CALL RGMSG(substr, RGMLIC, &
                  '         STEP   : ',cntptr%step,' ')
             CALL RGMSG(substr, RGMLIC, &
                  '         RESET  : ',cntptr%reset,' ')
             CALL RGMSG(substr, RGMLIC, &
                  '         CURRENT: ',cntptr%current,' ')
          END IF
       END IF
    END IF

  END SUBROUTINE RGTOOL_NCRGCNT_RST
! ----------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE GET_NEXT_NCRGCNT(last, cntptr)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    LOGICAL,              INTENT(OUT) :: last
    TYPE(ncrgcnt),        POINTER     :: cntptr

    ! LOCAL
    TYPE(t_ncrgcnt_list), POINTER, SAVE :: ai  => NULL()
    TYPE(t_ncrgcnt_list), POINTER, SAVE :: ae  => NULL()
    INTEGER, PARAMETER                  :: MODE_INIT = 1
    INTEGER, PARAMETER                  :: MODE_CONT = 2
    INTEGER,                       SAVE :: MODE = MODE_INIT

    ! INIT
    last = .FALSE.

    DO
       SELECT CASE (MODE)
       CASE(MODE_INIT)
          ai => GRGTLIST
          MODE = MODE_CONT
       CASE(MODE_CONT)
          IF (.NOT. ASSOCIATED(ai)) EXIT
          ae => ai
          ai => ai%next
          cntptr => ae%this
          RETURN
       END SELECT
    END DO

    MODE = MODE_INIT
    last = .TRUE.

  END SUBROUTINE GET_NEXT_NCRGCNT
  ! -------------------------------------------------------------------

! ----------------------------------------------------------------------
! PRIVATE SUBROUTINES
! ----------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_ncrgcnt(status, mname, cname, cntptr)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, LEN_TRIM, TRIM, ASSOCIATED

    ! I/O
    INTEGER,              INTENT(OUT) :: status
    CHARACTER(LEN=*),     INTENT(IN)  :: mname
    CHARACTER(LEN=*),     INTENT(IN)  :: cname
    TYPE(ncrgcnt),        POINTER     :: cntptr

    ! LOCAL
    TYPE(t_ncrgcnt_list), POINTER :: ai  => NULL()
    TYPE(t_ncrgcnt_list), POINTER :: ae  => NULL()
    LOGICAL                       :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(cntptr)

    ! CHECKS
    IF (LEN_TRIM(ADJUSTL(mname)) > NCCRSTRL) THEN
       status = 5003   ! NAME TOO LONG
       RETURN
    END IF
    !
    IF (LEN_TRIM(ADJUSTL(cname)) > NCCRSTRL) THEN
       status = 5006   ! NAME TOO LONG
       RETURN
    END IF

    ! CHECK, IF IT EXISTS
    ai => GRGTLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF ( (TRIM(ADJUSTL(mname)) == TRIM(ai%mname)) .AND. &
            (TRIM(ADJUSTL(cname)) == TRIM(ai%this%name)) ) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO    

    IF (lex) THEN
       cntptr => ai%this
    ELSE
       status = 5005   ! DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_ncrgcnt
  ! -------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE new_ncrgcnt(status, mname, cnt)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, LEN_TRIM, ASSOCIATED, TRIM, NULL

    ! I/O
    INTEGER,          INTENT(OUT)  :: status
    CHARACTER(LEN=*), INTENT(IN)   :: mname
    TYPE(NCRGCNT),    INTENT(IN)   :: cnt

    ! LOCAL
    TYPE(t_ncrgcnt_list), POINTER :: ai     => NULL()
    TYPE(t_ncrgcnt_list), POINTER :: ae     => NULL()
    TYPE(ncrgcnt),        POINTER :: cntptr => NULL()
    INTEGER                       :: zstat

    IF (LEN_TRIM(ADJUSTL(mname)) > NCCRSTRL) THEN
       status = 5003   ! SUBMODEL NAME TOO LONG
       RETURN
    END IF

    CALL loc_ncrgcnt(zstat, TRIM(mname), TRIM(cnt%name), cntptr)
    IF (zstat /= 5005) THEN  ! COUNTER DOES NOT EXIST (IS OK HERE !)
       IF (zstat == 0) THEN  ! COUNTER EXISTS ALREADY
          status = 5004      ! COUNTER EXISTS ALREADY
       ELSE
          status = zstat     ! ERROR
       END IF
       RETURN
    END IF

    ! GOTO END OF LIST
    ai => GRGTLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO

    ! ADD NEW
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(GRGTLIST)) THEN
       GRGTLIST => ai                 ! SET POINTER TO FIRST ELEMENT
    ELSE
       ae%next => ai                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF
       
    ! SET VALUES
    ai%mname            = TRIM(ADJUSTL(mname))
    !
    ai%this%name    = TRIM(ADJUSTL(cnt%name))
    ai%this%start   = cnt%start
    ai%this%step    = cnt%step
    ai%this%reset   = cnt%reset
    ai%this%current = cnt%current

    status = 0

  END SUBROUTINE new_ncrgcnt
  ! ----------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE ncrgcnt_halt(substr, status)

    ! MESSy
    USE messy_ncregrid_base,       ONLY: RGMSG, RGMLE
    USE messy_main_constants_mem,  ONLY: STRLEN_VLONG

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status
    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr
    
    IF (status == 0) RETURN
    
    errstr = ncrgcnt_error_str(status)
    
    CALL RGMSG(substr, RGMLE, TRIM(errstr))
    
  END SUBROUTINE ncrgcnt_halt
  ! -------------------------------------------------------------------
  
  ! -------------------------------------------------------------------
  FUNCTION ncrgcnt_error_str(status)
    
    USE messy_main_tools,         ONLY: int2str
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG

    IMPLICIT NONE
    
    ! I/O
    CHARACTER(LEN=STRLEN_VLONG)           :: ncrgcnt_error_str
    INTEGER,                   INTENT(IN) :: status
    
    ! LOCAL
    CHARACTER(LEN=4) :: echar = '    '
    
    ncrgcnt_error_str = ''
    
    CALL int2str(echar, status, cpad='0', cerr='*')
    
    SELECT CASE(echar)
       ! NO ERROR
    CASE('0000')
       ncrgcnt_error_str = 'E'//echar//': NO ERROR'

    CASE('5003')
       ncrgcnt_error_str = 'E'//echar//': SUBMODEL NAME TOO LONG'
    CASE('5004')
       ncrgcnt_error_str = 'E'//echar//': COUNTER EXISTS ALREADY'
    CASE('5005')
       ncrgcnt_error_str = 'E'//echar//': COUNTER DOES NOT EXIST'
    CASE('5006')
       ncrgcnt_error_str = 'E'//echar//': COUNTER NAME TOO LONG'

    CASE('5020')
       ncrgcnt_error_str = 'E'//echar//': COUNTER FILE NOT FOUND'

    CASE DEFAULT
       ncrgcnt_error_str = 'E'//echar//': UNKONW ERROR STATUS'

    END SELECT

  END FUNCTION ncrgcnt_error_str
  ! -------------------------------------------------------------------

! ******************************************************************
! ------------------------------------------------------------------
END MODULE MESSY_NCREGRID_TOOLS
! ------------------------------------------------------------------
! ******************************************************************

#endif
