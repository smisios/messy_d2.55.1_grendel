MODULE MESSY_MAIN_GRID_TOOLS

! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
!         Astrid  Kerkweg, UniMz, Mainz, 2012-2013
!         -> re-structured for more general GRID application
! ******************************************************************

  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: dp, r_e=> radius_earth

  IMPLICIT NONE

  PUBLIC  :: RGTOOL_CONVERT         ! CONVERTS NCVAR INTO 4D-ARRAY
  PUBLIC  :: RGTOOL_CONVERT_DAT2VAR ! CONVERTS 4D-ARRAY INTO NCVAR
  PUBLIC  :: RGTOOL_CONVERT_DAT2PREDEFVAR ! CONVERTS 4D-ARRAY INTO NCVAR
  PUBLIC  :: RGTOOL_G2C             ! CONVERTS GRID-INFORMATION INTO ARRAYS
  PUBLIC  :: ALINE_ARRAY
  PUBLIC  :: DEALINE_ARRAY
  PUBLIC  :: SET_SURFACE_PRESSURE   ! SET CURRENT SURFACE PRESSURE FOR VERTICAL
                                    ! GRID um_ak_20141007
  PUBLIC  :: UVROT2UV

CONTAINS

! ----------------------------------------------------------------------
! ******************************************************************
! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_CONVERT(var, dat, grid, order)

    ! CONVERTS NCREGRID OUTPUT OF TYPE N-ARRAY (var) ON GRID
    ! grid TO 4D-ARRAY (dat)
    ! THE OPTIONAL STRING order DEFINES THE ORDER OF DIMENSIONS
    ! (DEFAULT: 'xyzn')
    !
    ! Author: Patrick Joeckel, MPIC, October 2002

    USE MESSY_MAIN_GRID_NETCDF,  ONLY: t_ncvar, RGMSG , RGMLE, RGMLI, RGMLIC &
                                     , ERRMSG, QCMP_NCDIM, VTYPE_DOUBLE    &
                                     , VTYPE_INT, VTYPE_BYTE, VTYPE_CHAR   &
                                     , VTYPE_UNDEF, VTYPE_REAL, ELEMENT
    USE MESSY_MAIN_GRID,         ONLY: t_geohybgrid
  
    ! REGRID MODULES

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, INDEX, PRESENT, PRODUCT, REAL, SIZE, SUM, TRIM

    ! I/O
    TYPE (t_ncvar),    INTENT(IN)           :: var   ! nc-variable
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: dat   ! DATA ON DESTINATION GRID
    TYPE(t_geohybgrid),INTENT(IN)           :: grid  ! grid information
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
    INTEGER :: ix, iy           ! ub_ak_20170523

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
       IF (ASSOCIATED(grid%clonm%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%clonm%dim(1)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(1)) = var%dim(i)%len
             dpos(1) = i
          END IF
       END IF
       IF (ASSOCIATED(grid%clatm%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%clatm%dim(2)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(2)) = var%dim(i)%len
             dpos(2) = i
          END IF
       END IF
       IF (ASSOCIATED(grid%rlonm%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%rlonm%dim(1)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(1)) = var%dim(i)%len
             dpos(1) = i
          END IF
       END IF
       IF (ASSOCIATED(grid%rlatm%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%rlatm%dim(1)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(2)) = var%dim(i)%len
             dpos(2) = i
          END IF
       END IF
       IF (ASSOCIATED(grid%ulonm%dim)) THEN
          IF (QCMP_NCDIM(var%dim(i), grid%ulonm%dim(1)) > 1) THEN
             nrvdim = nrvdim+1
             ldvar(ovec(1)) = var%dim(i)%len
             dpos(1) = i
          END IF
          IF (SIZE(grid%ulonm%dim)>1) THEN
             IF (QCMP_NCDIM(var%dim(i), grid%ulonm%dim(2)) > 1) THEN
                nrvdim = nrvdim+1
                ldvar(ovec(2)) = var%dim(i)%len
                dpos(2) = i
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
    ! ADDITIONAL CHECK FOR HEIGHT GRIDS
    IF (dpos(3)  == 0 ) THEN
       DO i=1, var%ndims    
          ! EXCLUDE ALREADY ASSIGN HORIZONTAL DIMENSIONS
          ifpos: IF (i /= dpos(1) .AND. i /= dpos(2)) THEN
             IF (ASSOCIATED(grid%pressm%dim)) THEN
                DO ix = 1,SIZE(grid%pressm%dim)
                   IF (QCMP_NCDIM(var%dim(i), grid%pressm%dim(ix)) > 1) THEN
                      nrvdim = nrvdim+1
                      ldvar(ovec(3)) = var%dim(i)%len
                      dpos(3) = i
                      EXIT
                   END IF
                END DO
                IF (dpos(3) == 0 ) THEN
                   ! no dimension of pressm fits third var%dim. CHECK if
                   ! pressi exists, if not, check pressm%dim(ix) +1
                   IF (.NOT. ASSOCIATED(grid%pressi%dim)) THEN
                      DO ix = 1,SIZE(grid%pressm%dim)
                         DO iy = 1, var%ndims  
                            ! EXCLUDE HORIZONTAL DIM VARS
                            IF (QCMP_NCDIM(var%dim(iy), grid%pressm%dim(ix)) > 1) CYCLE
                            IF (grid%pressm%dim(ix)%len +1 == var%dim(i)%len) THEN
                               nrvdim = nrvdim+1
                               ldvar(ovec(3)) = var%dim(i)%len 
                               dpos(3) = i
                               EXIT
                            END IF
                         END DO
                         IF (dpos(3) > 0) EXIT
                      END DO
                   END IF
                END IF
                 !IF (.NOT. ASSOCIATED(grid%pressi%dim) THEN
             END IF
             IF (ASSOCIATED(grid%pressi%dim) .AND. dpos(3) == 0) THEN
                DO ix = 1,SIZE(grid%pressi%dim)
                   IF (QCMP_NCDIM(var%dim(i), grid%pressi%dim(ix)) > 1) THEN
                      nrvdim = nrvdim+1
                      ldvar(ovec(3)) = var%dim(i)%len 
                      dpos(3) = i
                      EXIT
                   END IF
                END DO
                IF (dpos(3) == 0 ) THEN
                   ! no dimension of pressi fits third var%dim. CHECK if
                   ! pressm exists, if not, check pressi%dim(ix) -1
                   IF (.NOT. ASSOCIATED(grid%pressm%dim)) THEN
                      DO ix = 1,SIZE(grid%pressi%dim)
                         DO iy = 1, var%ndims  
                            ! EXCLUDE HORIZONTAL DIM VARS
                            IF (QCMP_NCDIM(var%dim(iy), grid%pressi%dim(ix)) > 1) CYCLE
                            IF (grid%pressi%dim(ix)%len -1 == var%dim(i)%len) THEN
                               nrvdim = nrvdim+1
                               ldvar(ovec(3)) = var%dim(i)%len 
                               dpos(3) = i
                               EXIT
                            END IF
                         END DO
                         IF (dpos(3) > 0) EXIT
                      END DO
                   END IF
                END IF
             END IF
          END IF ifpos
       END DO
    END IF

    ldvar(ovec(4)) = npdlv/(ldvar(ovec(1))*ldvar(ovec(2))*ldvar(ovec(3)))

    CALL RGMSG(substr, RGMLI, 'VARIABLE '''//TRIM(var%name)//''' IS ' &
         ,var%ndims,'-DIMENSIONAL')
    CALL RGMSG(substr, RGMLIC, '(',var%dat%dim,')')
    CALL RGMSG(substr, RGMLIC, 'DATA DIMENSION-LENGHTS ARE: ')
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
               = REAL(var%dat%vr(i),dp)
       CASE(VTYPE_INT)
          dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
               = REAL(var%dat%vi(i),dp)
       CASE(VTYPE_BYTE)
          dat(nvec(1), nvec(2), nvec(3), nvec(4)) &
               = REAL(var%dat%vb(i),dp)
       CASE(VTYPE_CHAR)
          CALL RGMSG(substr, RGMLE, 'DATA OF TYPE CHAR NOT SUPPORTED !')
       CASE(VTYPE_UNDEF)
          CALL RGMSG(substr, RGMLE, 'UNDEFINED TYPE OF DATA !')
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF DATA !')
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

    USE MESSY_MAIN_GRID_NETCDF,  ONLY: t_ncvar, t_ncdim, RGMSG , RGMLE     &
                                     , ERRMSG, INIT_NCVAR, RENAME_NCVAR    &
                                     , QDEF_NCVAR, COPY_NCDIM, INIT_NARRAY &
                                     , ELEMENT, NF90_DOUBLE, VTYPE_DOUBLE  &
                                     , QCMP_NCDIM, INIT_NCDIM
    USE MESSY_MAIN_GRID,         ONLY: t_geohybgrid

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, INDEX, PRESENT, PRODUCT, SIZE, SUM, TRIM

    ! I/O
    TYPE (t_ncvar),      INTENT(INOUT)       :: var   ! nc-variable 
    REAL(DP), DIMENSION(:,:,:,:), POINTER    :: dat   ! DATA ON DESTINATION GRID
    CHARACTER(LEN=*),    INTENT(IN)          :: vname ! name of variable
    TYPE (t_geohybgrid), INTENT(IN)          :: grid  ! grid information
    CHARACTER(LEN=4),    INTENT(IN),OPTIONAL :: order ! DEFAULT: 'xyzn'

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_CONVERT_DAT2VAR'
    INTEGER :: ovec(4)          ! order of (x,y,z,n) for input field
    INTEGER :: dpos(3)          ! position of lon, lat, lev in variable
    INTEGER :: nrvdim           ! number of recognized dimensions in variable
    INTEGER :: npdlv            ! product of dimension-lengths in variable
    INTEGER :: i,j              ! counter
    INTEGER :: nvec(4)          ! counter parameter for dimensions
    INTEGER, DIMENSION(:), POINTER     :: vec     ! element vector
    INTEGER, DIMENSION(:), ALLOCATABLE :: zdimlen ! local dim. lengths
    INTEGER :: status           ! status flag
    INTEGER :: ix
    INTEGER :: length(5)
    TYPE (t_ncdim), DIMENSION(5) :: ncdim

    ! INIT
    CALL INIT_NCVAR(var)
    CALL RENAME_NCVAR(var, vname)
    var%xtype = NF90_DOUBLE ! xtype ?
    length(:)=1
    DO ix = 1, 4
       CALL INIT_NCDIM(ncdim(ix))
    END DO
    ! CREATE VARIABLE DIMENSIONS FROM GRID 
    nrvdim = 0
    ! 1st horizontal dimension
    IF (QDEF_NCVAR(grid%lonm)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(1),grid%lonm%dim(1))
       length(1) = grid%lonm%dim(1)%len
    ELSE IF (QDEF_NCVAR(grid%clonm)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(1),grid%clonm%dim(1))
       length(1) = grid%clonm%dim(1)%len
    ELSE IF (QDEF_NCVAR(grid%ulonm)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(1),grid%lonm%dim(1))
       length(1) = grid%ulonm%dim(1)%len
    ENDIF
    ! 2nd horizontal dimensiom
    IF (QDEF_NCVAR(grid%latm)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(2),grid%latm%dim(1))
       length(2) = grid%latm%dim(1)%len
    ELSE IF (QDEF_NCVAR(grid%clatm)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(2),grid%clatm%dim(2))
       length(2) = grid%clatm%dim(2)%len
    ELSE IF (QDEF_NCVAR(grid%ulatm)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(2),grid%ulatm%dim(1))
       length(2) = grid%ulatm%dim(1)%len
    ENDIF
    ! VERTICAL dimension
    IF (QDEF_NCVAR(grid%hyam)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(3),grid%hyam%dim(1))
       length(3) = grid%hyam%dim(1)%len
    ELSE IF (QDEF_NCVAR(grid%hybm)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(3),grid%hybm%dim(1))
       length(3) = grid%hyam%dim(1)%len
    ELSE IF (QDEF_NCVAR(grid%hyai)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(3),grid%hyai%dim(1))
       length(3) = grid%hyai%dim(1)%len - 1
    ELSE IF (QDEF_NCVAR(grid%hybi)) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(ncdim(3),grid%hybi%dim(1))
       length(3) = grid%hybi%dim(1)%len - 1
    ELSE IF (QDEF_NCVAR(grid%pressm)) THEN
       nrvdim = nrvdim + 1
       DO ix = 1,SIZE(grid%pressm%dim)
          IF (QCMP_NCDIM(ncdim(1), grid%pressm%dim(ix)) > 1) THEN
             CYCLE
          ELSE IF (QCMP_NCDIM(ncdim(2), grid%pressm%dim(ix)) > 1) THEN
             CYCLE
          ELSE
             EXIT
          END IF
       END DO
       CALL  COPY_NCDIM(ncdim(3),grid%pressm%dim(ix))
       length(3) = grid%pressm%dim(ix)%len
    ELSE IF  (QDEF_NCVAR(grid%pressi)) THEN
       nrvdim = nrvdim + 1
       DO ix = 1,SIZE(grid%pressi%dim)
          IF (QCMP_NCDIM(ncdim(1), grid%pressi%dim(ix)) > 1) THEN
             CYCLE
          ELSE IF (QCMP_NCDIM(ncdim(2), grid%pressi%dim(ix)) > 1) THEN
             CYCLE
          ELSE
             EXIT
          END IF
       END DO
       CALL  COPY_NCDIM(ncdim(3),grid%pressi%dim(1))
       length(3) = grid%pressi%dim(ix)%len - 1
    ENDIF
    IF (SIZE(dat)/(length(1)*length(2)*length(3)) > 1) THEN
       nrvdim = nrvdim + 1
       length(4) = SIZE(dat)/(length(1)*length(2)*length(3))
       ncdim(4)%len = length(4)
       ncdim(4)%name = 'PARAM '
       ncdim(4)%id = 99
    END IF
    IF (QDEF_NCVAR(grid%timem)) THEN
       nrvdim = nrvdim + 1
       IF (ncdim(4)%len > 0) THEN
          length(5) = 1
          !CALL COPY_NCDIM(ncdim(5),grid%timem%dim(1))
       ELSE
          length(4) = 1
          !CALL COPY_NCDIM(ncdim(4),grid%timem%dim(1))
       END IF
    END IF

    var%ndims = nrvdim
    ALLOCATE(var%dim(var%ndims))
    
    NULLIFY(vec)
    ovec = (/ 1,2,3,4 /)                                ! DEFAULT
    !ostr = (/ 'x/LON ', 'y/LAT ', 'z/LEV ', 'n/PAR '/) ! DEFAULT
    IF (PRESENT(order)) THEN
       ovec(1) = INDEX(order, 'x')
       !ostr(ovec(1)) = 'x/LON '
       ovec(2) = INDEX(order, 'y')
       !ostr(ovec(2)) = 'y/LAT '
       ovec(3) = INDEX(order, 'z')
       !ostr(ovec(3)) = 'z/LEV '
       ovec(4) = INDEX(order, 'n')
       !ostr(ovec(4)) = 'n/PAR '
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

    IF (TRIM(ncdim(1)%name) /= '' ) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(var%dim(nrvdim),ncdim(1))
       dpos(1) = 1
    END IF
    IF (TRIM(ncdim(2)%name) /= '' ) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(var%dim(nrvdim),ncdim(2))
       dpos(2) = 2
    END IF
    IF (TRIM(ncdim(3)%name) /= '' ) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(var%dim(nrvdim),ncdim(3))
       dpos(3) = 3
    END IF
    IF (TRIM(ncdim(4)%name) /= '' ) THEN
       nrvdim = nrvdim + 1
       CALL COPY_NCDIM(var%dim(nrvdim),ncdim(4))
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
       zdimlen(i) = length(i)   ! = 1 for UNLIM-DIM-ID
    END DO
    j = SIZE(zdimlen)
    CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_DOUBLE)
    DEALLOCATE(zdimlen)

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
    DO ix = 1, 4
       CALL INIT_NCDIM(ncdim(ix))
    END DO

  END SUBROUTINE RGTOOL_CONVERT_DAT2VAR
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_CONVERT_DAT2PREDEFVAR(stat, var, dat &
                                                , order_var, order_dat)

    ! CONVERTS 4D-ARRAY (dat) TO TYPE N-ARRAY (var%dat) ON 
    ! for a variable var with predefined dimensions
    ! THE OPTIONAL STRING order DEFINES THE ORDER OF DIMENSIONS
    ! (DEFAULT: 'xyzn')
    !
    ! Author: Astrid Kerkweg, Uni Bonn, Dezember 2016

    USE MESSY_MAIN_GRID_NETCDF,  ONLY: t_ncvar, RGMSG , RGMLE, ERRMSG    &
                                     , INIT_NARRAY, ELEMENT, NF90_DOUBLE &
                                     , VTYPE_DOUBLE, NULL_XTYPE

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, INDEX, PRESENT, PRODUCT, SIZE, SUM, TRIM

    ! I/O
    INTEGER, INTENT(OUT)                     :: stat
    TYPE (t_ncvar),      INTENT(INOUT)       :: var   ! nc-variable 
    ! um_ak_20140312-
    REAL(DP), DIMENSION(:,:,:,:), POINTER    :: dat   ! DATA ON DESTINATION GRID
    CHARACTER(LEN=4),    INTENT(IN)          :: order_var ! DEFAULT: 'xyzn'
    CHARACTER(LEN=4),    INTENT(IN),OPTIONAL :: order_dat ! DEFAULT: 'xyzn'

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_CONVERT_DAT2PREDEFVAR'
    INTEGER :: ovec(4)          ! order of (x,y,z,n) for input field
    INTEGER :: dpos(4)          ! position of lon, lat, lev in variable
    INTEGER :: npdlv            ! product of dimension-lengths in variable
    !CHARACTER(LEN=6) :: ostr(4)! order-string for output
    INTEGER :: i,j              ! counter
    INTEGER :: nvec(4)          ! counter parameter for dimensions
    INTEGER, DIMENSION(:), POINTER     :: vec     ! element vector
    INTEGER, DIMENSION(:), ALLOCATABLE :: zdimlen ! local dim. lengths
    INTEGER                            :: status
    ! INIT
    stat = 0
    IF (.NOT. ASSOCIATED(var%dim) .AND. var%xtype == NULL_XTYPE) THEN
       stat = 4000
       RETURN
    END IF

    NULLIFY(vec)
    ovec = (/ 1,2,3,4 /)                                ! DEFAULT
    !ostr = (/ 'x/LON ', 'y/LAT ', 'z/LEV ', 'n/PAR '/) ! DEFAULT
    IF (PRESENT(order_dat)) THEN
       ovec(1) = INDEX(order_dat, 'x')
       !ostr(ovec(1)) = 'x/LON '
       ovec(2) = INDEX(order_dat, 'y')
       !ostr(ovec(2)) = 'y/LAT '
       ovec(3) = INDEX(order_dat, 'z')
       !ostr(ovec(3)) = 'z/LEV '
       ovec(4) = INDEX(order_dat, 'n')
       !ostr(ovec(4)) = 'n/PAR '
       IF (SUM(ovec) /= 10) THEN
          CALL RGMSG(substr, RGMLE, &
               'ERROR IN '''//TRIM(order_dat)//'''-STRING AT SUBROUTINE CALL !')
       END IF
    END IF
    !
    dpos(:) = 0
    dpos(1) = INDEX(order_var, 'x')
    dpos(2) = INDEX(order_var, 'y')
    dpos(3) = INDEX(order_var, 'z')
    dpos(4) = INDEX(order_var, 'n')
     IF (SUM(ovec) /= 10) THEN
       CALL RGMSG(substr, RGMLE, &
            'ERROR IN '''//TRIM(order_var)//'''-STRING AT SUBROUTINE CALL !')
    END IF
    !
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
    var%xtype = NF90_DOUBLE
    CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_DOUBLE)
    DEALLOCATE(zdimlen) ! um_ak_20140312

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

  END SUBROUTINE RGTOOL_CONVERT_DAT2PREDEFVAR
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE RGTOOL_G2C(grid, hyam, hybm, p0, ps         &
                        ,hyai, hybi                      &
                        ,latm, lonm                      &
                        ,lati, loni, nlat, nlon, nlev    &
                        ,order)

    ! CONVERTS GRID INFORMATION (grid) TO ARRAYS
    ! hyam, hybm, hyai, hybi, latm, lati, lonm, loni, p0 -> 1D
    ! ps                                                 -> 2D
    !
    ! nlat, nlon, nlev are OPTIONAL parameters in order to specify
    ! the array-length (box-mid), in case the grid does not contain
    ! the respective dimension
    !
    ! Author: Patrick Joeckel, MPICH, October 2002

    USE MESSY_MAIN_GRID,         ONLY: t_geohybgrid
    USE MESSY_MAIN_GRID_NETCDF,  ONLY: ERRMSG

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT, SIZE

    ! I/O
    TYPE (t_geohybgrid),      INTENT(IN)        :: grid
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
    INTEGER,               INTENT(IN), OPTIONAL :: nlat
    INTEGER,               INTENT(IN), OPTIONAL :: nlon
    INTEGER,               INTENT(IN), OPTIONAL :: nlev
    CHARACTER(LEN=4),      INTENT(IN), OPTIONAL :: order  ! order of dimensions

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_G2C'
    CHARACTER(LEN=4)            :: loc_order
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: hha, hp0, hhb, hps, hlat, hlon
    INTEGER                     :: status

    NULLIFY(hha)
    NULLIFY(hp0)
    NULLIFY(hhb)
    NULLIFY(hps)
    NULLIFY(hlat)
    NULLIFY(hlon)
    
    IF (PRESENT(order)) THEN
       loc_order = order
    ELSE
       loc_order = 'xyzn'
    END IF

    IF (PRESENT(hyam)) THEN
       CALL RGTOOL_CONVERT(grid%hyam,hha,grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%hyai,hha,grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%p0, hp0, grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%hybm, hhb, grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%hybi, hhb, grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%ps, hps, grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%latm, hlat, grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%lati, hlat, grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%lonm, hlon, grid, order=loc_order)
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
       CALL RGTOOL_CONVERT(grid%loni, hlon, grid, order=loc_order)
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

  ! ===========================================================================
  SUBROUTINE ALINE_ARRAY(status, var2d, var1d)

    IMPLICIT NONE

    ! I/O
    INTEGER,                  INTENT(OUT) :: status
    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: var2d
    REAL(dp), DIMENSION(:),   POINTER     :: var1d
    ! LOCAL
    INTEGER :: i,j,n
    INTEGER :: idim
    INTEGER :: isize

    status = 3999

    ! ALLOCATE(var1d(isize))
    isize= SIZE(var2d)
    IF (ASSOCIATED(var1d)) THEN
       IF (SIZE(var1d) /= isize) THEN
          DEALLOCATE(var1d)
          ALLOCATE(var1d(isize))
       ENDIF
    ELSE
       ALLOCATE(var1d(isize))
    END IF

    idim = SIZE(var2d,1)
    DO i=1, idim
       DO j=1, SIZE(var2d,2)
          n = (j-1) * idim + i
          var1d(n) = var2d(i,j)
       END DO
    END DO

    status = 0

  END SUBROUTINE ALINE_ARRAY
  ! ========================================================================

  ! ===========================================================================
  SUBROUTINE DEALINE_ARRAY(status, dim1, dim2, var1d, var2d)

    IMPLICIT NONE

    ! I/O
    INTEGER,                 INTENT(OUT) :: status
    INTEGER,                 INTENT(IN)  :: dim1, dim2
    REAL(dp), DIMENSION(:),  INTENT(IN)  :: var1d
    REAL(dp), DIMENSION(:,:), POINTER    :: var2d
    ! LOCAL
    INTEGER :: i,j,n

    status = 3999

    IF (dim1 *dim2 /= SIZE(var1d)) THEN
       status = 3022 ! sizes incompatible
       RETURN
    ENDIF

    IF (ASSOCIATED(var2d)) DEALLOCATE(var2d)
    ALLOCATE(var2d(dim1,dim2))
    DO i=1, dim1
       DO j=1, dim2
          n = (j-1) * dim1 + i
          var2d(i,j) = var1d(n)
       END DO
    END DO

    status = 0
    
  END SUBROUTINE DEALINE_ARRAY
  ! ========================================================================
  ! ----------------------------------------------------------------------
  ! moved here from mmd2way_child

  SUBROUTINE SET_SURFACE_PRESSURE(status, grid, press, lcalc)

    USE MESSY_MAIN_GRID,         ONLY: t_geohybgrid
    USE MESSY_MAIN_GRID_NETCDF,  ONLY: NULL_VARID, NF90_FLOAT, VTYPE_DOUBLE &
                                     , POSITION, ERRMSG, INIT_NARRAY        &
                                     , INIT_NCVAR, COPY_NCDIM

    IMPLICIT NONE

    ! I/O
    INTEGER,                  INTENT(OUT)          :: status
    TYPE(t_geohybgrid),       INTENT(INOUT)        :: grid
    REAL(dp), DIMENSION(:,:), INTENT(IN)           :: press
    ! limit variable pressure to subparts
    LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lcalc

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'set_surface_pressure'
    INTEGER :: i, j

    CALL INIT_NCVAR(grid%ps)
    grid%ps%name  = 'ps'
    grid%ps%id    = NULL_VARID
    grid%ps%xtype = NF90_FLOAT
    ! ... dimensions
    IF (.NOT. ASSOCIATED(grid%ps%dim)) THEN
       grid%ps%ndims = 2
       ALLOCATE(grid%ps%dim(grid%ps%ndims), STAT=status)
       CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)
    ELSE
       IF (SIZE(grid%ps%dim) /= 2) THEN
          status = 2100 
          RETURN
          ! 'ps grid dimensions do not match'
       END IF
    END IF

    ! horizontal non-rotated grid has been defined 
    CALL COPY_NCDIM(grid%ps%dim(1), grid%lonm%dim(1))
    CALL COPY_NCDIM(grid%ps%dim(2), grid%latm%dim(1))

    ! ... data
    CALL INIT_NARRAY(grid%ps%dat, grid%ps%ndims       &
         , (/grid%ps%dim(1)%len, grid%ps%dim(2)%len/) &
         , VTYPE_DOUBLE)
    DO i=1,grid%ps%dim(1)%len
       DO j=1, grid%ps%dim(2)%len
          IF (PRESENT(lcalc)) THEN
             IF (lcalc(i,j)) THEN
                grid%ps%dat%vd(&
                     POSITION((/grid%ps%dim(1)%len, grid%ps%dim(2)%len/)&
                     ,(/i,j/)))  = press(i,j)
             ELSE
                grid%ps%dat%vd(&
                     POSITION((/grid%ps%dim(1)%len, grid%ps%dim(2)%len/) &
                     ,(/i,j/)))   =  101325.0_dp
             END IF
          ELSE
             grid%ps%dat%vd(&
                  POSITION((/grid%ps%dim(1)%len, grid%ps%dim(2)%len/)&
                  ,(/i,j/)))  = press(i,j)
          END IF
       END DO
    END DO

    status = 0

  END SUBROUTINE SET_SURFACE_PRESSURE
! ---------------------------------------------------------------------------
! subroutine is copied from cosmo/src/utilities.f90
SUBROUTINE uvrot2uv (urot, vrot, rlat, rlon, pollat, pollon, u, v)
!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the rotated system
!   to the real geographical system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

! Parameter list:
REAL (DP), INTENT (IN)          ::    &
  urot, vrot,     & ! wind components in the rotated grid
  rlat, rlon,     & ! latitude and longitude in the true geographical system
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

REAL (DP), INTENT (OUT)         ::    &
  u, v              ! wind components in the true geographical system

! Local variables

REAL (DP)                       ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm

! added PARAMETER attribute to enable inlining and vectorization 
! (otherwise zrpi18 / zpir18 get an implicit "SAVE" attribute which hinders inlining!)
REAL (DP), PARAMETER   ::    &
  zrpi18 = 57.2957795_dp,         & !
  zpir18 = 0.0174532925_dp
!------------------------------------------------------------------------------
! Begin subroutine uvrot2uv
!------------------------------------------------------------------------------

! Converting from degree to radians
  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)
  zlonp   = (pollon-rlon) * zpir18
  zlat    =         rlat  * zpir18

  zarg1   = zcospol*SIN(zlonp)
  zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
  znorm   = 1.0/SQRT(zarg1**2 + zarg2**2)

! Convert the u- and v-components
  u       =   urot*zarg2*znorm + vrot*zarg1*znorm
  v       = - urot*zarg1*znorm + vrot*zarg2*znorm

END SUBROUTINE uvrot2uv

!------------------------------------------------------------------------------

SUBROUTINE uvrot2uv_vec(urot, vrot, rlat, rlon, pollat, pollon, idim, jdim,u,v)

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the rotated
!   system to the real geographical system. This is the vectorized form
!   of the routine above, i.e. the computation is for a whole 2D field.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Parameter list:
INTEGER       , INTENT(IN)     ::    &
  idim, jdim        ! dimensions of the field

REAL (KIND=dp), INTENT (IN) ::    &
  urot  (idim,jdim), & ! wind components in the true geographical system
  vrot  (idim,jdim)    !

REAL (KIND=dp), INTENT (IN)    ::    &
  rlat(idim,jdim),& ! coordinates in the true geographical system
  rlon(idim,jdim),& !
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

REAL (KIND=dp), INTENT (OUT) ::    &
  u  (idim,jdim), & ! wind components in the true geographical system
  v  (idim,jdim)    !

! Local variables
REAL (KIND=dp)                 ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zugeo, zvgeo

INTEGER                        ::    i, j

! added PARAMETER attribute to enable inlining and vectorization
! (otherwise zrpi18 / zpir18 get an implicit "SAVE" attribute which hinders inlining!)
REAL (KIND=dp),     PARAMETER            ::    &
  zrpi18 = 57.2957795_dp,                      & !
  zpir18 = 0.0174532925_dp

!------------------------------------------------------------------------------
! Begin Subroutine uvrot2uv_vec
!------------------------------------------------------------------------------

!$acc data present(u, v, rlat, rlon)

! Converting from degree to radians
  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)

  !$acc parallel
  !$acc loop gang
  DO j = 1, jdim
    !$acc loop vector private (zlonp, zlat, zarg1, zarg2, znorm, zugeo, zvgeo)
    DO i = 1, idim

      zlonp   = (pollon-rlon(i,j)) * zpir18
      zlat    =         rlat(i,j)  * zpir18

      zarg1   = zcospol*SIN(zlonp)
      zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
      znorm   = 1.0_dp/SQRT(zarg1**2 + zarg2**2)

      ! Convert the u- and v-components
      zugeo   =  urot(i,j)*zarg2*znorm + vrot(i,j)*zarg1*znorm
      zvgeo   = -urot(i,j)*zarg1*znorm + vrot(i,j)*zarg2*znorm
      u(i,j) = zugeo
      v(i,j) = zvgeo

    ENDDO
  ENDDO
  !$acc end parallel

!$acc end data

END SUBROUTINE uvrot2uv_vec

!==============================================================================
END MODULE MESSY_MAIN_GRID_TOOLS
