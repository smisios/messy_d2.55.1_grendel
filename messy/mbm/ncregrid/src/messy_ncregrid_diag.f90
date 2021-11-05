! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_DIAG
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  USE netcdf,                    ONLY: NF90_FLOAT, NF90_DOUBLE  &
                                     , NF90_SHORT, NF90_BYTE    &
                                     , NF90_INT, NF90_CHAR
  USE MESSY_NCREGRID_BASE,       ONLY: narray, axis             &
                                     , VTYPE_UNDEF, VTYPE_INT   &
                                     , VTYPE_REAL, VTYPE_DOUBLE &
                                     , VTYPE_BYTE, VTYPE_CHAR   &
                                     , RGMSG                    &
                                     , RGMLVL, RGMLVLC          &
                                     , RGMLVM, RGMLVMC          &
                                     , RGMLI, RGMLIC            &
                                     , RGMLW, RGMLWC            &
                                     , RGMLE, RGMLEC
  USE MESSY_NCREGRID_NETCDF,     ONLY: ncdim, ncvar, ncatt      &
                                     , NULL_VARID, NULL_XTYPE   &
                                     , NULL_DIMID
  USE MESSY_NCREGRID_GEOHYB,     ONLY: geohybgrid

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, MAXVAL, MINVAL, PRESENT, SIZE, TRIM
  PRIVATE   :: ASSOCIATED, MAXVAL, MINVAL, PRESENT, SIZE, TRIM

  INTERFACE STDOUT
     MODULE PROCEDURE WRITE_NARRAY
     MODULE PROCEDURE WRITE_AXIS
     MODULE PROCEDURE WRITE_NCDIM
     MODULE PROCEDURE WRITE_NCVAR
     MODULE PROCEDURE WRITE_NCATT
     MODULE PROCEDURE WRITE_GEOHYBGRID
  END INTERFACE

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE WRITE_NARRAY(a, ldout)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(IN)     :: a
  LOGICAL, OPTIONAL, INTENT(IN) :: ldout

  ! LOCAL
  INTEGER :: vtype
  LOGICAL :: l

  IF (PRESENT(ldout)) THEN
     l = ldout
  ELSE
     l = .true.
  END IF

  WRITE(*,*) '--- NARRAY -----------'
  WRITE(*,*) 'Number of dimensions :',a%n
  IF (ASSOCIATED(a%dim)) THEN
     WRITE(*,*) 'Length of dimensions :',a%dim
  END IF

  vtype = a%type

  SELECT CASE(vtype)
  CASE(VTYPE_INT)
     WRITE(*,*) 'Type                 : INTEGER'
     WRITE(*,*) 'Size                 : ',SIZE(a%vi)
     WRITE(*,*) 'Range                : ',MINVAL(a%vi),MAXVAL(a%vi)
     IF (l) THEN
        WRITE(*,*) 'Values               : ',a%vi
     END IF
  CASE(VTYPE_REAL)
     WRITE(*,*) 'Type                 : REAL'
     WRITE(*,*) 'Size                 : ',SIZE(a%vr)
     WRITE(*,*) 'Range                : ',MINVAL(a%vr),MAXVAL(a%vr)
     IF (l) THEN
        WRITE(*,*) 'Values               : ',a%vr
     END IF
  CASE(VTYPE_DOUBLE)
     WRITE(*,*) 'Type                 : DOUBLEPRECISION'
     WRITE(*,*) 'Size                 : ',SIZE(a%vd)
     WRITE(*,*) 'Range                : ',MINVAL(a%vd),MAXVAL(a%vd)
     IF (l) THEN
        WRITE(*,*) 'Values               : ',a%vd
     END IF
  CASE(VTYPE_CHAR)
     WRITE(*,*) 'Type                 : CHARACTER'
     WRITE(*,*) 'Size                 : ',SIZE(a%vc)
     IF (l) THEN
        WRITE(*,*) 'Values               : ',a%vc
     END IF
  CASE(VTYPE_BYTE)
     WRITE(*,*) 'Type                 : BYTE'
     WRITE(*,*) 'Size                 : ',SIZE(a%vb)
     WRITE(*,*) 'Range                : ',MINVAL(a%vb),MAXVAL(a%vb)
     IF (l) THEN
        WRITE(*,*) 'Values               : ',a%vb
     END IF
  CASE(VTYPE_UNDEF)
     WRITE(*,*) 'Type                 : UNDEFINED'
     WRITE(*,*) 'Values               : <EMPTY>'
  CASE DEFAULT
     WRITE(*,*) 'Unrecognized type of data !'
  END SELECT

  WRITE(*,*) '----------------------'

END SUBROUTINE WRITE_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE WRITE_AXIS(a, ldout)

  IMPLICIT NONE

  ! I/O
  TYPE (axis), INTENT(IN) :: a
  LOGICAL, OPTIONAL, INTENT(IN) :: ldout

  ! LOCAL
  LOGICAL :: l

  IF (PRESENT(ldout)) THEN
     l = ldout
  ELSE
     l = .true.
  END IF

  WRITE(*,*) '=== AXIS =================='
  IF (a%lm) THEN
     WRITE(*,*) 'Axis is modulo'
  ELSE
     WRITE(*,*) 'Axis is non-modulo'
  END IF
  WRITE(*,*) 'Number of dependencies    :',a%ndp
  IF (ASSOCIATED(a%dep)) THEN
     WRITE(*,*) 'Dependencies              :',a%dep
  ELSE
     WRITE(*,*) 'Dependencies              : <EMPTY>'
  END IF
  CALL WRITE_NARRAY(a%dat, l)
  WRITE(*,*) '==========================='

END SUBROUTINE WRITE_AXIS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE WRITE_NCDIM(dim)

  IMPLICIT NONE

  ! I/O
  TYPE (ncdim), INTENT(IN) :: dim

  WRITE(*,*) 'NCDIM>----------------'
  WRITE(*,*) 'DIMENSION: ',TRIM(dim%name)
  WRITE(*,*) '       ID: ',dim%id
  WRITE(*,*) '      LEN: ',dim%len
  IF (dim%fuid) THEN
     WRITE(*,*)  '         : -> UNLIMITED DIMENSION'
  END IF
  IF (dim%varid /= NULL_VARID) THEN
     WRITE(*,*)  '   VAR-ID: ',dim%varid
  END IF
  WRITE(*,*) '<NCDIM----------------'

END SUBROUTINE WRITE_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE WRITE_NCATT(att, ldout)

  IMPLICIT NONE

  ! I/O
  TYPE (ncatt), INTENT(IN)      :: att
  LOGICAL, OPTIONAL, INTENT(IN) :: ldout ! OUTPUT data

  ! LOCAL
  LOGICAL :: l

  IF (PRESENT(ldout)) THEN
     l = ldout
  ELSE
     l = .true.
  END IF

  WRITE(*,*) 'ATTRIBUTE: ',TRIM(att%name)
  WRITE(*,*) '      NO.: ',att%num
  WRITE(*,*) '      LEN: ',att%len
  WRITE(*,*) '   VAR-ID: ',att%varid
  SELECT CASE(att%xtype)
     CASE(NF90_BYTE)
        WRITE(*,*) '     TYPE: NF90_BYTE'
     CASE(NF90_SHORT)
        WRITE(*,*) '     TYPE: NF90_SHORT'
     CASE(NF90_CHAR)
        WRITE(*,*) '     TYPE: NF90_CHAR'
     CASE(NF90_INT)
        WRITE(*,*) '     TYPE: NF90_INT'
     CASE(NF90_FLOAT)
        WRITE(*,*) '     TYPE: NF90_FLOAT'
     CASE(NF90_DOUBLE)
        WRITE(*,*) '     TYPE: NF90_DOUBLE'
     CASE DEFAULT
        WRITE(*,*) '     TYPE: *** UNKNOWN OR UNDEFINED ***'
  END SELECT

  WRITE(*,*)  'ATT-VALUE(S): '
  CALL WRITE_NARRAY(att%dat, l)

END SUBROUTINE WRITE_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE WRITE_NCVAR(var, ldout)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar), INTENT(IN)      :: var
  LOGICAL, OPTIONAL, INTENT(IN) :: ldout ! OUTPUT DATA ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'WRITE_NCVAR'
  INTEGER :: i
  LOGICAL :: l

  IF (PRESENT(ldout)) THEN
     l = ldout
  ELSE
     l = .true.
  END IF

  WRITE(*,*) 'NCVAR>----------------------------------------------'
  WRITE(*,*) 'VARIABLE : ',TRIM(var%name)
  WRITE(*,*) '       ID: ',var%id
  SELECT CASE(var%xtype)
     CASE(NF90_BYTE)
        WRITE(*,*) '     TYPE: NF90_BYTE'
     CASE(NF90_SHORT)
        WRITE(*,*) '     TYPE: NF90_SHORT'
     CASE(NF90_CHAR)
        WRITE(*,*) '     TYPE: NF90_CHAR'
     CASE(NF90_INT)
        WRITE(*,*) '     TYPE: NF90_INT'
     CASE(NF90_FLOAT)
        WRITE(*,*) '     TYPE: NF90_FLOAT'
     CASE(NF90_DOUBLE)
        WRITE(*,*) '     TYPE: NF90_DOUBLE'
     CASE DEFAULT
        WRITE(*,*) '     TYPE: *** UNKNOWN OR UNDEFINED ***'
  END SELECT

  ! DIMENSIONS
  WRITE(*,*) '    NDIMS: ',var%ndims
  ! UNLIMITED DIM ?
  IF (var%uid /= -1) THEN
     WRITE(*,*) '    U-ID : ',var%uid
     WRITE(*,*) '    step : ',var%ustep
  END IF
  ! ATTRIBUTES
  WRITE(*,*) '    NATTS: ',var%natts

  IF (.NOT.ASSOCIATED(var%dim).AND.(var%ndims > 0)) THEN
     CALL RGMSG(substr, RGMLE, 'NUMBER OF DIMENSIONS IS ', &
          var%ndims,' ', .false.)
     CALL RGMSG(substr, RGMLEC, 'BUT DIM-ARRAY IS NOT ASSOCIATED !')
  END IF

  IF ((SIZE(var%dim) /= var%ndims).AND.(SIZE(var%dim) /= 1)) THEN
     CALL RGMSG(substr, RGMLE, &
          'NUMBER OF DIMENSIONS IS NOT CONFORM WITH', .false.)
     CALL RGMSG(substr, RGMLEC, 'ARRAY SIZE (',SIZE(var%dim),')')
  END IF

  DO i=1, var%ndims
     CALL WRITE_NCDIM(var%dim(i))
  END DO

  IF (.NOT.ASSOCIATED(var%att).AND.(var%natts > 0)) THEN
     CALL RGMSG(substr, RGMLE, 'NUMBER OF ATTRIBUTES IS ', &
          var%natts,' ', .false.)
     CALL RGMSG(substr, RGMLEC, 'BUT ATT-ARRAY IS NOT ASSOCIATED !')
  END IF

  IF ((SIZE(var%att) /= var%natts).AND.(SIZE(var%att) /= 1)) THEN
     CALL RGMSG(substr, RGMLE, &
          'NUMBER OF ATTRIBUTES IS NOT CONFORM WITH', .false.)
     CALL RGMSG(substr, RGMLEC, 'ARRAY SIZE (',SIZE(var%att),')')
  END IF

  DO i=1, var%natts
     CALL WRITE_NCATT(var%att(i))
  END DO

  WRITE(*,*) 'VAR-VALUE(S): '
  CALL WRITE_NARRAY(var%dat, l)

  WRITE(*,*) '<NCVAR----------------------------------------------'

END SUBROUTINE WRITE_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE WRITE_GEOHYBGRID(grid, ldout)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(IN) :: grid
  LOGICAL, OPTIONAL, INTENT(IN) :: ldout

  ! LOCAL
  LOGICAL :: l

  IF (PRESENT(ldout)) THEN
     l = ldout
  ELSE
     l = .true.
  END IF

  WRITE(*,*) 'GEOHYBGRID>*****************************************'
  WRITE(*,*) 'FILE     : ',TRIM(grid%file)
  WRITE(*,*) 'time-step: ',grid%t

  CALL WRITE_NCVAR(grid%lonm, l)
  CALL WRITE_NCVAR(grid%loni, l)
  CALL WRITE_NCVAR(grid%latm, l)
  CALL WRITE_NCVAR(grid%lati, l)
  CALL WRITE_NCVAR(grid%col, l)
  CALL WRITE_NCVAR(grid%hyam, l)
  CALL WRITE_NCVAR(grid%hyai, l)
  CALL WRITE_NCVAR(grid%p0, l)
  CALL WRITE_NCVAR(grid%hybm, l)
  CALL WRITE_NCVAR(grid%hybi, l)
  CALL WRITE_NCVAR(grid%ps, l)
  CALL WRITE_NCVAR(grid%timem, l)
  CALL WRITE_NCVAR(grid%timei, l)

  WRITE(*,*) '<GEOHYBGRID*****************************************'

END SUBROUTINE WRITE_GEOHYBGRID
! ------------------------------------------------------------------

! ******************************************************************
END MODULE MESSY_NCREGRID_DIAG
! ******************************************************************
