!*************************************************************************
MODULE messy_offemis
!*************************************************************************

  ! MODULE FOR EMISSION OF TRACERS (IMPORTED 'OFFLINE' FROM
  ! EXTERNAL netCDF FILES)
  ! MESSy-SMCL
  !
  ! Authors:
  !    Patrick Joeckel, MPICH, March 2004
  !    Astrid Kerkweg, UNI-MZ, Oct 2009, modified from OFFLEM to OFFEMIS
  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'offemis'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'

  ! SUBROUTINES
  PUBLIC :: parse_str
  PUBLIC :: parse_zstr
  PUBLIC :: parse_datatstr ! op_mm_20140116

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE parse_str(status, strlen, str, icl, ilg, igp) ! op_sb_20191024

    USE messy_main_tools, ONLY: strcrack

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    INTEGER,          INTENT(INOUT) :: icl      ! LG method (CLaMS)
    INTEGER,          INTENT(INOUT) :: ilg      ! LG method (ATTILA)
    INTEGER,          INTENT(INOUT) :: igp      ! GP method

    ! LOCAL
    CHARACTER(LEN=strlen), POINTER  :: sl1(:)
    CHARACTER(LEN=strlen), POINTER  :: sl2(:)
    INTEGER :: n, m
    INTEGER :: i
    INTEGER :: iostat

    status = 1 ! ERROR

    NULLIFY(sl1)
    NULLIFY(sl2)

    CALL strcrack(str, ';', sl1, n)
    DO i=1, n

       CALL strcrack(sl1(i), '=', sl2, m)
       IF (SIZE(sl2) == 2) THEN
          IF (TRIM(ADJUSTL(sl2(2))) == '') THEN
             status = 2    ! EMPTY SPECIFICATION
             RETURN
          END IF
       END IF

       SELECT CASE(TRIM(ADJUSTL(sl2(1))))
          CASE('LG')
             READ(sl2(2),*,IOSTAT=iostat) ilg
             IF (iostat /= 0) THEN
                status = 4  ! ERROR IN READING INTEGER
                RETURN
             END IF

             IF ((ilg < 0) .OR. (ilg > 4)) THEN
                status = 6
                RETURN
             ENDIF

           CASE('CL')
             READ(sl2(2),*,IOSTAT=iostat) icl
             IF (iostat /= 0) THEN
                status = 44  ! ERROR IN READING INTEGER
                RETURN
             END IF
             IF ((icl < 0) .OR. (icl > 4)) THEN
                status = 66
                RETURN
             ENDIF

          CASE('GP')
             READ(sl2(2),*,IOSTAT=iostat) igp
             IF (iostat /= 0) THEN
                status = 5  ! ERROR IN READING INTEGER
                RETURN
             END IF
             IF ((igp < 0) .OR. (igp > 2)) THEN
                status = 7
                RETURN
             ENDIF

          CASE DEFAULT
             status = 3 ! UNKNOWN SPECIFIER
             RETURN
       END SELECT

    END DO

    ! CLEAN UP
    IF (ASSOCIATED(sl1)) DEALLOCATE(sl1)
    IF (ASSOCIATED(sl2)) DEALLOCATE(sl2)

    status = 0 ! NO ERROR

  END SUBROUTINE parse_str
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE parse_zstr(status, strlen, str, z)

    USE messy_main_tools, ONLY: strcrack

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    REAL(DP), DIMENSION(:), POINTER :: z        ! level information

    ! LOCAL
    CHARACTER(LEN=strlen), POINTER  :: sl1(:)
    INTEGER :: j, l
    INTEGER :: iostat

    status = 1 ! ERROR

    NULLIFY(sl1)

    IF (ASSOCIATED(z)) THEN
       DEALLOCATE(z)
       NULLIFY(z)
    END IF

    CALL strcrack(str, ',', sl1, l)
    ALLOCATE(z(l))
    DO j=1, l
       READ(sl1(j),*,IOSTAT=iostat) z(j)
       IF (iostat /= 0) THEN
          status = 6  ! ERROR IN READING REAL
          RETURN
       END IF
    END DO

    IF (ASSOCIATED(sl1)) DEALLOCATE(sl1)

    status = 0 ! NO ERROR

  END SUBROUTINE parse_zstr
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
 SUBROUTINE parse_datatstr(status, strlen, str, name, subname, scale, n)

    !AUTHOR:
    !       Mariano Mertens, DLR, 2013 
    !
    !modified version from parse_f2tstr (onemis) written by
    ! Patrick Joeckel, MPICH, 2004
    ! Astrid  Kerkweg, MPICH, 2004

    USE messy_main_tools,         ONLY: strcrack, str2num
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: name    ! (OUT)
    CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: subname ! (OUT)
    REAL(DP),                         DIMENSION(:), POINTER :: scale   ! (OUT)
    ! LOCAL
    CHARACTER(LEN=strlen),               POINTER     :: sl1(:)
    CHARACTER(LEN=strlen),               POINTER     :: sl2(:)
    CHARACTER(LEN=2*STRLEN_MEDIUM+1)                 :: tmpstr = ''
    INTEGER :: n, m
    INTEGER :: i
    INTEGER :: ix

    status = 0 ! NO ERROR
    n = 0

    NULLIFY(sl1)
    NULLIFY(sl2)

    !the different entries are seperated by a ;
    !Each entry then consists of tracername_subname,scalingfactor
    CALL strcrack(str, ';', sl1, n) 
    ! -> n strings with tracer[_subname][,scaling];
    
    IF (n == 0) THEN
       status = -1
       RETURN
    END IF

    ALLOCATE(name(n))
    name(:) = ''
    ALLOCATE(subname(n))
    subname(:) = ''
    ALLOCATE(scale(n))
    scale = 1.0_DP
    
    ! LOOP OVER ALL TRACERS
    DO i=1, n

       ! CHECK FOR SCALING FACTOR (AFTER ',')
       CALL strcrack(sl1(i), ',', sl2, m) 

       tmpstr = ''
       SELECT CASE(m)
          CASE(0)
             ! empty string
             status = 1 
             return
          CASE(1)
             ! no ',' -> no scaling factor
             tmpstr = sl2(1)
             scale(i) = 1.0_dp
             status = 0
          CASE(2)
             ! one ',' -> scaling factor
             tmpstr = sl2(1)
             CALL str2num(sl2(2),scale(i), status)
             IF (status /= 0) RETURN
          CASE DEFAULT
             ! more than one ','
             status = 3
             return
       END SELECT

       ix = INDEX(tmpstr, '_')
       IF (ix > 0) THEN
          name(i)     = TRIM(tmpstr(:ix-1))
          subname(i)  = TRIM(tmpstr(ix+1:))
       ELSE
          name(i)     = TRIM(tmpstr)
       END IF
       
    END DO

    ! CLEAN UP
    IF (ASSOCIATED(sl1)) DEALLOCATE(sl1)
    IF (ASSOCIATED(sl2)) DEALLOCATE(sl2)

    status = 0 ! NO ERROR
    
  END SUBROUTINE parse_datatstr
!-------------------------------------------------------------------------

!*************************************************************************
END MODULE messy_offemis
!*************************************************************************
