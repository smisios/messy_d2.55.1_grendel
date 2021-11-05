! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_MAIN_GRID_NETCDF
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
!         Astrid  Kerkweg, UniMz, Mainz, 2012-2013
!         -> re-structured for more general GRID application
! ******************************************************************

  USE NETCDF
  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: dp, sp, i4, i8

  IMPLICIT NONE
  PUBLIC

  INTRINSIC :: ASSOCIATED, SIZE, TRIM, PRESENT, MIN, PRODUCT    &
             , LEN_TRIM, INT, MAXVAL, MINVAL, REAL, NULL, IAND

  PRIVATE   :: ASSOCIATED, SIZE, TRIM, PRESENT, MIN, PRODUCT    &
             , LEN_TRIM, INT, MAXVAL, MINVAL, REAL, NULL, IAND

  INTEGER, PARAMETER :: GRD_MAXSTRLEN = 200
  INTEGER, PARAMETER :: NULL_DIMID = -1
  INTEGER, PARAMETER :: NULL_VARID = -1
  INTEGER, PARAMETER :: NULL_XTYPE = -1

  ! DATA TYPES
  INTEGER, PARAMETER :: VTYPE_UNDEF  = 0
  INTEGER, PARAMETER :: VTYPE_INT    = 1
  INTEGER, PARAMETER :: VTYPE_REAL   = 2
  INTEGER, PARAMETER :: VTYPE_DOUBLE = 3
  INTEGER, PARAMETER :: VTYPE_BYTE   = 4
  INTEGER, PARAMETER :: VTYPE_CHAR   = 5

  ! MESSAGE LEVEL TYPES
  INTEGER, PARAMETER :: RGMLE   = 0   ! ERROR
  INTEGER, PARAMETER :: RGMLEC  = 1   ! ERROR CONTINUED
  INTEGER, PARAMETER :: RGMLVL  = 2   ! LITTLE VERBOSE
  INTEGER, PARAMETER :: RGMLVLC = 3   ! LITTLE VERBOSE CONTINUED
  INTEGER, PARAMETER :: RGMLW   = 4   ! WARNING
  INTEGER, PARAMETER :: RGMLWC  = 5   ! WARNING CONTINUED
  INTEGER, PARAMETER :: RGMLVM  = 6   ! MEDIUM VERBOSE
  INTEGER, PARAMETER :: RGMLVMC = 7   ! MEDIUM VERBOSE CONTINUED
  INTEGER, PARAMETER :: RGMLI   = 8   ! INFO
  INTEGER, PARAMETER :: RGMLIC  = 9   ! INFO CONTINUED

  ! MESSAGE OUTPUT LEVEL
  INTEGER, PARAMETER :: MSGMODE_S  =  0  ! SILENT
  INTEGER, PARAMETER :: MSGMODE_E  =  1  ! ERROR MESSAGES
  INTEGER, PARAMETER :: MSGMODE_VL =  2  ! LITTLE VERBOSE
  INTEGER, PARAMETER :: MSGMODE_W  =  4  ! WARNING MESSAGES
  INTEGER, PARAMETER :: MSGMODE_VM =  8  ! MEDIUM VERBOSE
  INTEGER, PARAMETER :: MSGMODE_I  = 16  ! INFO MESSAGES
  INTEGER, SAVE      :: MSGMODE    = MSGMODE_E
                                         !MSGMODE_S + MSGMODE_E + MSGMODE_VL &
                                         !+ MSGMODE_W + MSGMODE_VM + MSGMODE_I
  INTEGER, SAVE      :: MY_RANK    = -1

  TYPE t_narray
     ! n-dimenional array as 1D (LINEAR) array (REAL)
     INTEGER                              :: type = VTYPE_UNDEF
     INTEGER                              :: n = 0        ! number of dimensions
     INTEGER     , DIMENSION(:), POINTER  :: dim => NULL()! dim. vector
     REAL (SP)   , DIMENSION(:), POINTER  :: vr => NULL() ! real values
     REAL (DP)   , DIMENSION(:), POINTER  :: vd => NULL() ! double values
     INTEGER (I8), DIMENSION(:), POINTER  :: vi => NULL() ! integer values
     INTEGER (I4), DIMENSION(:), POINTER  :: vb => NULL() ! byte values
     CHARACTER,    DIMENSION(:), POINTER  :: vc => NULL() ! char. values
  END TYPE t_narray

  TYPE t_ncdim
     CHARACTER(LEN=GRD_MAXSTRLEN) :: name = ''          ! name of dimension
     INTEGER                      :: id   = NULL_DIMID  ! dimension ID
     INTEGER                      :: len  = 0           ! length of dimension
     LOGICAL                      :: fuid = .false.     ! flag for "UNLIMITED"
                                                        ! ... dimension
     INTEGER                      :: varid = NULL_VARID ! variable ID,
                                                        ! ... if coordinate var
  END TYPE t_ncdim

  TYPE t_ncatt
     CHARACTER(LEN=GRD_MAXSTRLEN) :: name  = ''         ! attribute name
     INTEGER                      :: ID    = 0          ! number of attribute
     INTEGER                      :: xtype = NULL_XTYPE ! type of attribute
     INTEGER                      :: varid = NULL_VARID ! variable ID or
                                                        ! ... NF90_GLOBAL
     INTEGER                      :: len   = 0          ! length of attribute
     TYPE(t_narray)               :: dat                ! content of attribute
  END TYPE t_ncatt

  TYPE t_ncvar
     CHARACTER(LEN=GRD_MAXSTRLEN)         :: name  = ''         ! variable name
     INTEGER                              :: id    = NULL_VARID ! variable ID
     INTEGER                              :: xtype = NULL_XTYPE ! type of
                                                               ! ... variable
     INTEGER                              :: ndims = 0          ! number of
                                                               ! ... dimensions
     TYPE(t_ncdim), DIMENSION(:), POINTER :: dim => NULL()      ! netCDF
                                                               ! ... dimensions
     INTEGER                              :: uid   = NULL_DIMID ! unlimited
                                                               ! ... dim ID
     INTEGER                              :: ustep = 0  ! step along unlim. ID
     INTEGER                              :: natts = 0  ! number of attributes
     TYPE(t_ncatt), DIMENSION(:), POINTER :: att => NULL() ! list of attributes
     TYPE(t_narray)                       :: dat           ! content of variable
  END TYPE t_ncvar

  ! for multi-netcdf descriptor files
  TYPE t_multinc
     ! already scanned
     LOGICAL                                             :: l = .FALSE.
     ! number of files
     INTEGER                                             :: nf = 0
     ! list of filenames
     CHARACTER(LEN=GRD_MAXSTRLEN), DIMENSION(:), POINTER :: flist => NULL()
     ! number of time steps per file
     INTEGER,                      DIMENSION(:), POINTER :: nt => NULL()
  END type t_multinc

  TYPE t_ncatt_array
     INTEGER                                  :: natts = 0
     TYPE(t_ncatt), DIMENSION(:),   POINTER   :: atts  => NULL()
  END type t_ncatt_array

!! NOTE: DOES NOT WORK PROPERLY FOR SOME COMPILERS ...
!  INTERFACE ASSIGNMENT (=)
!     MODULE PROCEDURE COPY_NCDIM
!     MODULE PROCEDURE COPY_NCATT
!     MODULE PROCEDURE COPY_NCVAR
!  END INTERFACE

INTERFACE QSORT
   ! THE GOOD OLD (RECURSIVE) QUICKSORT ALGORITHM FOR LINEAR ARRAYS
   MODULE PROCEDURE QSORT_I    ! INTEGER
#if !(defined(__SX__))
   MODULE PROCEDURE QSORT_B    ! BYTE
   MODULE PROCEDURE QSORT_R    ! REAL
#endif
   MODULE PROCEDURE QSORT_D    ! DOUBLE PRECISION
END INTERFACE

INTERFACE RGMSG                    ! MESSAGE OUTPUT
   MODULE PROCEDURE RGMSG_C
   MODULE PROCEDURE RGMSG_I
   MODULE PROCEDURE RGMSG_IA
   MODULE PROCEDURE RGMSG_R
   MODULE PROCEDURE RGMSG_D
END INTERFACE

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE INIT_NCDIM(dim)

  ! initialise netcdf dimension struct t_ncdim

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncdim), INTENT(INOUT) :: dim

  dim%name  = ''
  dim%id    = NULL_DIMID
  dim%len   = 0
  dim%fuid  = .false.
  dim%varid = NULL_VARID

END SUBROUTINE INIT_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_NCATT(att)

  ! initialise netcdf attribute struct t_ncatt

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncatt), INTENT(INOUT) :: att

  att%name  = ''
  att%ID    = 0
  att%xtype = NULL_XTYPE
  att%len   = 0
  att%varid = NULL_VARID

  CALL INIT_NARRAY(att%dat)

END SUBROUTINE INIT_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_NCVAR(var)

  ! initialise netcdf variable struct t_ncvar

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar), INTENT(INOUT) :: var

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER   :: substr = 'INIT_NCVAR'
  INTEGER                       :: i
  INTEGER                       :: status

  var%name  = ''
  var%id    = NULL_VARID
  var%xtype = NULL_XTYPE
  var%uid   = NULL_DIMID
  var%ustep = 0

  IF (ASSOCIATED(var%dim)) THEN
     DO i=1, SIZE(var%dim)
        CALL INIT_NCDIM(var%dim(i))
     END DO
     DEALLOCATE(var%dim, STAT=status)
     CALL ERRMSG(substr,status,1)
  END IF
  NULLIFY(var%dim)
  var%ndims = 0

  IF (ASSOCIATED(var%att)) THEN
     DO i=1, SIZE(var%att)
        CALL INIT_NCATT(var%att(i))
     END DO
     DEALLOCATE(var%att, STAT=status)
     CALL ERRMSG(substr,status,2)
  END IF
  NULLIFY(var%att)
  var%natts = 0

  CALL INIT_NARRAY(var%dat)

END SUBROUTINE INIT_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_NCDIM(dd, ds)

  ! copy netcdf dimension struct t_ncdim

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncdim), INTENT(INOUT) :: dd   ! destination
  TYPE (t_ncdim), INTENT(IN)    :: ds   ! source

  CALL INIT_NCDIM(dd)

  dd%name  = TRIM(ds%name)
  dd%id    = ds%id
  dd%len   = ds%len
  dd%fuid  = ds%fuid
  dd%varid = ds%varid

END SUBROUTINE COPY_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_NCATT(ad, as)

  ! copy netcdf attribute struct t_ncatt

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncatt), INTENT(INOUT) :: ad   ! destination
  TYPE (t_ncatt), INTENT(IN)    :: as   ! source

  CALL INIT_NCATT(ad)

  ad%name  = TRIM(as%name)
  ad%ID   = as%ID
  ad%xtype = as%xtype
  ad%len   = as%len
  ad%varid = as%varid

  CALL COPY_NARRAY(ad%dat, as%dat)

END SUBROUTINE COPY_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_NCVAR(vd, vs)

  ! copa netcdf variable struct t_ncvar

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar), INTENT(INOUT) :: vd   ! destination
  TYPE (t_ncvar), INTENT(IN)    :: vs   ! source

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'COPY_NCVAR'
  INTEGER :: i
  INTEGER :: status

  CALL INIT_NCVAR(vd)

  vd%name  = TRIM(vs%name)
  vd%id    = vs%id
  vd%xtype = vs%xtype
  vd%uid   = vs%uid
  vd%ustep = vs%ustep

  vd%ndims = vs%ndims
  IF (ASSOCIATED(vs%dim)) THEN
     ALLOCATE(vd%dim(SIZE(vs%dim)), STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, vs%ndims
        CALL COPY_NCDIM(vd%dim(i), vs%dim(i))
     END DO
  END IF

  vd%natts = vs%natts
  IF (ASSOCIATED(vs%att)) THEN
     ALLOCATE(vd%att(SIZE(vs%att)), STAT=status)
     CALL ERRMSG(substr,status,2)
     DO i=1, vs%natts
        CALL COPY_NCATT(vd%att(i), vs%att(i))
     END DO
  END IF

  CALL COPY_NARRAY(vd%dat, vs%dat)

END SUBROUTINE COPY_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
LOGICAL FUNCTION QDEF_NCVAR(var)

  ! test if netcdf variable struct (t_ncvar) is defined

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar), INTENT(IN)  :: var

  QDEF_NCVAR = (var%dat%type /= VTYPE_UNDEF)

END FUNCTION QDEF_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------

LOGICAL FUNCTION QCMP_NARRAY(a1, a2)

  ! Compare two structures of type t_narray

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(IN) :: a1, a2  ! arrays to compare
  ! LOCAL
  INTEGER                     :: i

  QCMP_NARRAY = .FALSE.

  IF (a1%n    /= a2%n )   RETURN
     IF (SIZE(a1%dim) /= SIZE(a2%dim)) RETURN
     DO i = 1, SIZE(a1%dim)
        IF (a1%dim(i) /= a2%dim(i)) RETURN
     END DO

  IF (a1%type /= a2%type) RETURN
  SELECT CASE(a1%type)
  CASE(VTYPE_UNDEF)
     ! if both narrays have undefined variable they are equal ...
  CASE(VTYPE_INT)
     IF (SIZE(a1%vi) /= SIZE(a2%vi))  RETURN
     DO i=1,SIZE(a1%vi)
        IF (a1%vi(i) /= a2%vi(i)) RETURN
     END DO
  CASE(VTYPE_REAL)
     IF (SIZE(a1%vr) /= SIZE(a2%vr)) RETURN
     DO i=1,SIZE(a1%vr)
        IF (a1%vr(i) /= a2%vr(i)) RETURN
     END DO
  CASE(VTYPE_DOUBLE)
     IF (SIZE(a1%vd) /= SIZE(a2%vd)) RETURN
     DO i=1,SIZE(a1%vd)
        IF (a1%vd(i) /= a2%vd(i)) RETURN
     END DO
  CASE(VTYPE_BYTE)
     IF (SIZE(a1%vb) /= SIZE(a2%vb)) RETURN
     DO i=1,SIZE(a1%vb)
        IF (a1%vb(i) /= a2%vb(i)) RETURN
     END DO
  CASE(VTYPE_CHAR)
    IF (SIZE(a1%vc) /= SIZE(a2%vc)) RETURN
     DO i=1,SIZE(a1%vc)
        IF (a1%vc(i) /= a2%vc(i)) RETURN
     END DO
   END SELECT

  QCMP_NARRAY = .TRUE.

END FUNCTION QCMP_NARRAY

 !*****************************************************

LOGICAL FUNCTION QCMP_NCATT(att1, att2)

  ! Compare two structures of type t_ncatt

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncatt), INTENT(IN) :: att1, att2  ! arrays to compare

  QCMP_NCATT = .FALSE.

  IF (att1%xtype       /= att2%xtype )      RETURN
  IF (att1%varid       /= att2%varid)       RETURN
  IF (att1%len         /= att2%len)         RETURN

  IF (.NOT. QCMP_NARRAY(att1%dat,att2%dat)) RETURN

  QCMP_NCATT = .TRUE.

END FUNCTION QCMP_NCATT

 !*****************************************************

LOGICAL FUNCTION QCMP_NCVAR(v1, v2)

  ! Compare two structures of type t_ncvar

  IMPLICIT NONE

  ! I/O
  TYPE (t_NCVAR), INTENT(IN) :: v1, v2  ! variables to compare

  ! LOCAL
  INTEGER :: i,j

  QCMP_NCVAR = .FALSE.

  IF (v1%xtype /= v2%xtype) RETURN
  IF (v1%ndims /= v2%ndims) RETURN
  DO i = 1, v1%ndims
     IF (QCMP_NCDIM(v1%dim(i),v2%dim(i)) < 2) RETURN
  END DO

  IF (v1%uid   /= v2%uid)   RETURN
  IF (v1%ustep /= v2%ustep) RETURN
  IF (v1%natts /= v2%natts) RETURN
  ! FIND EVERY ATTRIBUT OF V1 in V2
  DO i=1, v1%natts
     DO j=1, v2%natts
        IF (QCMP_NCATT(v1%att(i),v2%att(j))) EXIT
     END DO
     IF (j > v2%natts) RETURN
  END DO

  IF (.NOT. QCMP_NARRAY(v1%dat, v2%dat) ) RETURN

  QCMP_NCVAR = .TRUE.

END FUNCTION QCMP_NCVAR

! ------------------------------------------------------------------

! ------------------------------------------------------------------
INTEGER  FUNCTION QCMP_NCDIM(d1, d2)

  ! Compare two structures of type t_ncdim

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncdim), INTENT(IN) :: d1, d2  ! dimensions to compare

  QCMP_NCDIM = 0

  IF ((d1%len == 0).OR.(d2%len == 0)) THEN
     RETURN
  END IF

  IF (d1%len == d2%len) QCMP_NCDIM = QCMP_NCDIM + 1
  IF ( ((TRIM(d1%name) == TRIM(d2%name)).AND.(TRIM(d1%name) /= '')).OR. &
       ((d1%id == d2%id).AND.(d1%id /= NULL_DIMID)) )                   &
       QCMP_NCDIM = QCMP_NCDIM + 1

END FUNCTION QCMP_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IMPORT_NCDIM(dim, dimname, dimid, file, ncid)

  ! read structure of type t_ncdim from netcdf file

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncdim),   INTENT(INOUT)         :: dim      ! dimension
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: dimname  ! name of dimension
  INTEGER,          INTENT(IN),  OPTIONAL :: dimid    ! dimension ID
  INTEGER,          INTENT(IN),  OPTIONAL :: ncid     ! netCDF file ID
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: file     ! filename

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMPORT_NCDIM'
  INTEGER                     :: lncid    ! local netCDF file-ID
  INTEGER                     :: uid      ! unlimited ID
  INTEGER                     :: status   ! netCDF status

  ! CHECK CALL OF SUBROUTINE
  IF ((.NOT.PRESENT(dimname)).AND.(.NOT.PRESENT(dimid))) THEN
     CALL RGMSG(substr, RGMLE, 'NAME OR ID OF DIMENSION MUST BE SPECIFIED !')
  END IF
  !
  IF ((.NOT.PRESENT(ncid)).AND.(.NOT.PRESENT(file))) THEN
     CALL RGMSG(substr, RGMLE, &
          'netCDF-FILE-ID OR FILENAME MUST BE SPECIFIED !')
  END IF

  ! INIT
  CALL INIT_NCDIM(dim)

  ! OPEN FILE
  IF (PRESENT(ncid)) THEN
     lncid = ncid
  ELSE
     ! OPEN FILE
     CALL NFERR(substr                                 &
          ,nf90_open(TRIM(file),NF90_NOWRITE,lncid)    &
          ,1)
  END IF

  IF (PRESENT(dimname)) THEN
     ! GET DIM-ID from DIMNAME
     dim%name = TRIM(dimname)
     IF (TRIM(dim%name) /= '') THEN
        CALL NFERR(substr                                 &
             ,nf90_inq_dimid(lncid,TRIM(dim%name),dim%id) &
             ,2)
     ELSE
        CALL RGMSG(substr, RGMLE, 'INVALID NAME OF DIMENSION !')
     END IF
  ELSE
     dim%id = dimid
     ! GET DIM-NAME from DIM-ID
     CALL NFERR(substr                                           &
          ,nf90_Inquire_Dimension(lncid, dim%id, name=dim%name)  &
          ,3)
  END IF

  ! GET DIM-LEN
  CALL NFERR(substr                                              &
             ,nf90_Inquire_Dimension(lncid, dim%id, len=dim%len) &
             ,4)

  ! CHECK IF UNLIMITED
  CALL NFERR(substr                                   &
             ,nf90_Inquire(lncid, unlimitedDimID=uid) &
             ,5)
  IF (uid == -1) uid = NULL_DIMID
  IF ((uid == dim%id).AND.(uid /= NULL_DIMID)) THEN
     dim%fuid = .true.
  ELSE
     dim%fuid = .false.
  END IF

  ! CHECK FOR VARIABLE WITH SAME NAME (= COORDINATE VARIABLE)
  status = nf90_inq_varid(lncid, TRIM(dim%name), dim%varid)
  IF (status /= NF90_NOERR) THEN
     dim%varid = NULL_VARID
  END IF

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL NFERR(substr       &
          ,nf90_close(lncid) &
          ,6)
  END IF

END SUBROUTINE IMPORT_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_NCDIM(dim, file, ncid)

  ! write structure of type t_ncdim to netcdf file

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncdim),   INTENT(INOUT)           :: dim
  CHARACTER(LEN=*), INTENT(IN),   OPTIONAL  :: file ! filename
  INTEGER,          INTENT(IN),   OPTIONAL  :: ncid ! netCDF file ID

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'EXPORT_NCDIM'
  LOGICAL :: lex      ! file exists ?
  INTEGER :: lncid    ! local netCDF file-ID
  INTEGER :: uid      ! unlimited ID
  INTEGER :: status   ! netCDF status
  INTEGER :: zdimlen  ! local dimension length

  ! CHECK CALL OF SUBROUTINE
  IF ((.NOT.PRESENT(ncid)).AND.(.NOT.PRESENT(file))) THEN
     CALL RGMSG(substr, RGMLE, &
          'netCDF-FILE-ID OR FILENAME MUST BE SPECIFIED !')
  END IF

  ! CHECK DIMNAME
  IF (TRIM(dim%name) == '') THEN
     CALL RGMSG(substr, RGMLE, 'DIMENSION IS NAMELESS !')
  END IF

  IF (PRESENT(ncid)) THEN
     lncid = ncid
     ! DEFINE MODE
     CALL RGMSG(substr, RGMLI, &
          'netCDF-FILE WITH ID ',lncid,': SWITCHING TO DEFINE MODE !')
     status = nf90_redef(lncid)

  ELSE
     ! CHECK IF FILE EXISTS
     INQUIRE(file=TRIM(file), exist=lex)
     IF (lex) THEN   ! FILE EXISTS
        CALL RGMSG(substr, RGMLI, &
             'netCDF-FILE '''//TRIM(file)//''' EXISTS ! OPENING ...')
        CALL NFERR(substr                            &
             ,nf90_open(TRIM(file),NF90_WRITE,lncid) &
             ,2)
        CALL NFERR(substr        &
             ,nf90_redef(lncid)  &
             ,3)
     ELSE            ! NEW FILE
        CALL RGMSG(substr, RGMLI, &
             'netCDF-FILE '''//TRIM(file)//''' DOES NOT EXIST ! NEW FILE ...')
        CALL NFERR(substr                                &
             ,nf90_create(TRIM(file),NF90_CLOBBER,lncid) &
             ,4)
     END IF
  END IF

  CALL RGMSG(substr, RGMLI, &
       'EXPORTING DIMENSION '''//TRIM(dim%name)//''' ...' )

  ! CHECK IF DIMENSION EXISTS ALREADY
  status = nf90_inq_dimid(lncid,TRIM(dim%name), dim%id)
  IF (status == NF90_NOERR) THEN  ! DIM exists
     CALL RGMSG(substr, RGMLIC, ' ... DIMENSION EXISTS')
     ! CHECK IF UNLIMITED ID
     CALL NFERR(substr                              &
          ,nf90_Inquire(lncid, unlimitedDimID=uid)  &
          ,5)
     IF (uid == -1) uid = NULL_DIMID
     ! GET (CURRENT) LENGTH
     CALL NFERR(substr                                        &
          ,nf90_Inquire_Dimension(lncid, dim%id, len=zdimlen) &
          ,6)
     IF ((uid == dim%id).AND.(uid /= NULL_DIMID)) THEN   ! UNLIMITED
        CALL RGMSG(substr, RGMLIC, ' ... IS UNLIMITED')
        IF (.NOT.dim%fuid) THEN
           CALL RGMSG(substr, RGMLW, &
                'EXPORT DIMENSION HAS NO U-ID FLAG ! SETTING  ...')
           dim%fuid = .true.
        END IF
     ELSE                                        ! LIMITED
        IF (dim%fuid) THEN
           CALL RGMSG(substr, RGMLW, &
                'EXPORT DIMENSION HAS U-ID FLAG ! UNSETTING  ...')
           dim%fuid = .false.
           IF (zdimlen == dim%len) THEN
              CALL RGMSG(substr, RGMLIC, ' ... HAS CORRECT LENGTH')
           ELSE
              CALL RGMSG(substr, RGMLIC, ' ... HAS DIFFERENT LENGTH')
              CALL RGMSG(substr, RGMLE, &
                   'DIMENSION '''//TRIM(dim%name)//'''', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'IN FILE '''//TRIM(file)//'''', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'EXISTS ALREADY WITH LENGTH ',zdimlen,' !')
           END IF
        END IF
     END IF                                      ! (UN)LIMITED
  ELSE                            ! DIM does not exist
     IF (dim%fuid) THEN
        CALL RGMSG(substr, RGMLIC, ' ... NEW (UNLIMITED)')
        CALL NFERR(substr                       &
             ,nf90_def_dim(lncid,TRIM(dim%name) &
             ,NF90_UNLIMITED, dim%id)           &
             ,7)
     ELSE
        CALL RGMSG(substr, RGMLIC, ' ... NEW (CONSTANT LENGTH)')
        CALL NFERR(substr                       &
             ,nf90_def_dim(lncid,TRIM(dim%name) &
             ,dim%len, dim%id)                  &
             ,8)
     END IF
  END IF                          ! DIM exists ?

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL NFERR(substr        &
          ,nf90_close(lncid)  &
          ,9)
  END IF

END SUBROUTINE EXPORT_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IMPORT_NCATT(att, varname, varid, attname, attid   &
                         ,file, ncid                          &
                         ,lnostop)

  ! read structure of type t_ncatt from netcdf file

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncatt),   INTENT(INOUT)          :: att      ! attribute
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL  :: varname  ! variable name
  INTEGER,          INTENT(IN),  OPTIONAL  :: varid    ! netCDF variable ID
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL  :: attname  ! attribute name
  INTEGER,          INTENT(IN),  OPTIONAL  :: attid    ! attribute ID
  INTEGER,          INTENT(IN),  OPTIONAL  :: ncid     ! netCDF file ID
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL  :: file     ! filename
  LOGICAL,          INTENT(IN),  OPTIONAL  :: lnostop  ! do not stop if
                                                       ! att. does not exist

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER      :: substr = 'IMPORT_NCATT'
  INTEGER                          :: lncid    ! local netCDF file-ID
  INTEGER                          :: znatts   ! local number of attributes
  CHARACTER(LEN=GRD_MAXSTRLEN)     :: lvarname ! local variable name
  CHARACTER(LEN=GRD_MAXSTRLEN)     :: lattname ! local attribute name
  CHARACTER(LEN=100*GRD_MAXSTRLEN) :: astr = ''
  INTEGER                          :: RGMLH, RGMLHC    ! RG-message-level switch
  INTEGER                          :: i
  LOGICAL                          :: lstp             ! STOP ON ERROR IF TRUE
  INTEGER                          :: status

  ! CHECK CALL OF SUBROUTINE
  IF ((.NOT.PRESENT(varname)).AND.(.NOT.PRESENT(varid))) THEN
     CALL RGMSG(substr, RGMLE, 'NAME OR ID OF VARIABLE MUST BE SPECIFIED !')
  END IF
  !
  IF ((.NOT.PRESENT(attname)).AND.(.NOT.PRESENT(attid ))) THEN
     CALL RGMSG(substr, RGMLE, &
          'NAME OR NUMBER OF ATTRIBUTE MUST BE SPECIFIED !')
  END IF
  !
  IF ((.NOT.PRESENT(ncid)).AND.(.NOT.PRESENT(file))) THEN
     CALL RGMSG(substr, RGMLE, &
          'netCDF-FILE-ID OR FILENAME MUST BE SPECIFIED !')
  END IF

  ! INIT
  CALL INIT_NCATT(att)
  IF (PRESENT(lnostop)) THEN
     lstp = .NOT.lnostop
  ELSE
     lstp = .true.             ! DEFAULT: STOP IF ATTRIBUTE IS NOT PRESENT
  END IF
  ! SET RG-MESSAGE LEVEL
  IF (lstp) THEN
     RGMLH  = RGMLE   ! ERROR AND STOP
     RGMLHC = RGMLEC  ! ERROR AND STOP
  ELSE
     RGMLH  = RGMLW   ! WARNING AND GO ON
     RGMLHC = RGMLWC  ! WARNING AND GO ON
  END IF

  ! OPEN FILE
  IF (PRESENT(ncid)) THEN
     lncid = ncid
  ELSE
     ! OPEN FILE
     CALL NFERR(substr                              &
          ,nf90_open(TRIM(file),NF90_NOWRITE,lncid) &
          ,1)
  END IF

  ! GET VAR-ID
  IF (PRESENT(varname)) THEN
     lvarname = TRIM(varname)
     ! GET VAR-ID from VARNAME
     IF (TRIM(lvarname) /= '') THEN
        CALL NFERR(substr                                    &
             ,nf90_inq_varid(lncid,TRIM(lvarname),att%varid) &
             ,2)
     ELSE
        CALL RGMSG(substr, RGMLE, 'INVALID VARIABLE NAME !')
     END IF
  ELSE
     ! GET VAR-NAME from VAR-ID
     att%varid = varid
     IF (varid /= NF90_GLOBAL) THEN
        ! GET VAR-NAME from VAR-ID
        CALL NFERR(substr                                            &
             ,nf90_Inquire_Variable(lncid, att%varid, name=lvarname) &
             ,3)
     ELSE
        lvarname = ''
     END IF
  END IF

  ! GET ATTID
  IF (PRESENT(attname)) THEN         ! ATTRIBUTE NAME
     lattname = TRIM(attname)
     att%name = TRIM(attname)
     ! GET NUM, TYPE, LEN  from ATTNAME
     IF (TRIM(lattname) /= '') THEN
        status = nf90_Inquire_Attribute(lncid, att%varid   &
                 ,lattname                                 &
                 ,att%xtype, att%len, att%ID)
        IF (status /= NF90_NOERR) THEN
           CALL RGMSG(substr, RGMLH, &
                'ATTRIBUTE '''//TRIM(att%name)//'''', .false.)
           CALL RGMSG(substr, RGMLHC, &
                'OF VARIABLE WITH ID ',att%varid,' ', .false.)
           CALL RGMSG(substr, RGMLHC, 'DOES NOT EXIST !', lstp)
           !
           ! IF lstp = .false. : JUMP BACK
           ! CLOSE FILE
           IF (.NOT.PRESENT(ncid)) THEN
              CALL NFERR(substr       &
                   ,nf90_close(lncid) &
                   ,4)
           END IF
           RETURN
        END IF
     ELSE
        CALL RGMSG(substr, RGMLE, 'INVALID ATTRIBUTE NAME !')
     END IF
  ELSE                              ! ATTRIBUTE NUMBER
     ! GET NUMBER OF ATTRIBUTES FOR VARIABLE
     CALL NFERR(substr                                           &
          ,nf90_Inquire_Variable(lncid, att%varid, natts=znatts) &
          ,5)
     ! CHECK RANGE
     IF (attid  > znatts) THEN
        CALL RGMSG(substr, RGMLH, &
             'ATTRIBUTE WITH NUMBER',attid ,' ',.false.)
        CALL RGMSG(substr, RGMLHC, &
             'OF VARIABLE WITH ID ',att%varid,' ', .false.)
        CALL RGMSG(substr, RGMLHC, 'DOES NOT EXIST !', lstp)
        !
        ! IF lstp = .false. : JUMP BACK
        ! CLOSE FILE
        IF (.NOT.PRESENT(ncid)) THEN
           CALL NFERR(substr       &
                ,nf90_close(lncid) &
                ,6)
        END IF
        RETURN
     END IF

     att%ID = attid
     ! GET ATT-NAME FROM VAR-ID AND NUMBER
     CALL NFERR(substr                       &
          ,nf90_inq_attname(lncid,att%varid  &
          ,att%ID, att%name)                &
          ,7)
     ! GET ATT-LENGTH and TYPE
     CALL NFERR(substr                                        &
          ,nf90_Inquire_Attribute(lncid, att%varid, att%name  &
          ,xtype=att%xtype, len=att%len)                      &
          ,8)
  END IF                           ! ATTRIBUTE NAME OR NUMBER

  ! READ VALUE
  SELECT CASE(att%xtype)
  CASE(NF90_BYTE)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_BYTE )
     CALL NFERR(substr                     &
          ,nf90_get_att(lncid, att%varid   &
          ,TRIM(att%name)                  &
          ,values=att%dat%vb)              &
          ,9)
  CASE(NF90_SHORT)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_INT )
     CALL NFERR(substr                    &
          ,nf90_get_att(lncid, att%varid  &
          ,TRIM(att%name)                 &
          ,values=att%dat%vi)             &
          ,10)
  CASE(NF90_CHAR)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_CHAR )
     CALL NFERR(substr                    &
          ,nf90_get_att(lncid, att%varid  &
          ,TRIM(att%name)                 &
!          ,values=att%dat%vc)             &
          ,values = astr)                 &
          ,11)
     DO i=1, MIN(att%len, 100*GRD_MAXSTRLEN)
        att%dat%vc(i) = astr(i:i)
     END DO
  CASE(NF90_INT)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_INT )
     CALL NFERR(substr                    &
          ,nf90_get_att(lncid, att%varid  &
          ,TRIM(att%name)                 &
          ,values=att%dat%vi)             &
          ,12)
  CASE(NF90_FLOAT)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_REAL )
     CALL NFERR(substr                    &
          ,nf90_get_att(lncid, att%varid  &
          ,TRIM(att%name)                 &
          ,values=att%dat%vr)             &
          ,13)
  CASE(NF90_DOUBLE)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_DOUBLE )
     CALL NFERR(substr                    &
          ,nf90_get_att(lncid, att%varid  &
          ,TRIM(att%name)                 &
          ,values=att%dat%vd)             &
          ,14)
  CASE DEFAULT
     CALL RGMSG(substr, RGMLW,  'UNKNOWN ATTRIBUTE TYPE ',att%xtype,' ')
     CALL RGMSG(substr, RGMLWC, 'OF VARIABLE WITH ID ',att%varid,' ')
     CALL RGMSG(substr, RGMLWC, 'IN FILE '''//TRIM(file)//''' !')
  END SELECT

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL NFERR(substr       &
          ,nf90_close(lncid) &
          ,15)
  END IF

END SUBROUTINE IMPORT_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_NCATT(att, file, ncid, clobber)

  ! write structure of type t_ncatt to netcdf file

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncatt),   INTENT(INOUT)           :: att     ! attribute
  CHARACTER(LEN=*), INTENT(IN),   OPTIONAL  :: file    ! filename
  INTEGER,          INTENT(IN),   OPTIONAL  :: ncid    ! netCDF file ID
  LOGICAL,          INTENT(IN),   OPTIONAL  :: clobber ! overwrite ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'EXPORT_NCATT'
  INTEGER :: lncid    ! local netCDF file ID
  INTEGER :: zxtype   ! local attribute type
  INTEGER :: zattid   ! local attribute number
  INTEGER :: status   ! netcdf status
  LOGICAL :: lex      ! file exists ?
  LOGICAL :: lclobber ! overwrite ?
  INTEGER :: alen     ! attribute length

  ! CHECK CALL OF SUBROUTINE
  IF ((.NOT.PRESENT(ncid)).AND.(.NOT.PRESENT(file))) THEN
     CALL RGMSG(substr, RGMLE, &
          'netCDF-FILE-ID OR FILENAME MUST BE SPECIFIED !')
  END IF

  ! CHECK ATTNAME AND VARID
  IF ((TRIM(att%name) == '').OR.(att%varid == NULL_VARID)) THEN
     CALL RGMSG(substr, RGMLE, &
          'NAME AND VARIABLE-ID OF ATTRIBUTE MUST BE SPECIFIED !')
  END IF

  IF (PRESENT(clobber)) THEN
     lclobber = clobber
  ELSE
     lclobber = .false.   ! DEFAULT: DO NOT OVERWRITE
  END IF

  IF (PRESENT(ncid)) THEN
     lncid = ncid
     ! DEFINE MODE
     CALL RGMSG(substr, RGMLI, &
          'netCDF-FILE WITH ID ',lncid,': SWITCHING TO DEFINE MODE !')
     status = nf90_redef(lncid)
  ELSE
     ! CHECK IF FILE EXISTS
     INQUIRE(file=TRIM(file), exist=lex)
     IF (lex) THEN   ! FILE EXISTS
        ! OPEN FILE
        CALL RGMSG(substr, RGMLI, &
             'FILE '''//TRIM(file)//''' EXISTS ! OPENING ...')
        CALL NFERR(substr                            &
             ,nf90_open(TRIM(file),NF90_WRITE,lncid) &
             ,2)
        CALL NFERR(substr       &
             ,nf90_redef(lncid) &
             ,3)
     ELSE
        CALL RGMSG(substr, RGMLE, &
             'FILE '''//TRIM(file)//''' NOT FOUND !')
        ! STOP
     END IF
  END IF

  ! CHECK IF ATTRIBUTE EXISTS
  status = nf90_Inquire_Attribute(lncid,att%varid         &
            ,att%name, xtype=zxtype, attnum =zattid )
  IF (status == NF90_NOERR) THEN ! attribute exists
     CALL RGMSG(substr, RGMLW, &
          'ATTRIBUTE '''//TRIM(att%name)//'''')
     CALL RGMSG(substr, RGMLWC, &
          'OF VARIABLE WITH ID ',att%varid,' ')
     IF (PRESENT(file)) THEN
        CALL RGMSG(substr, RGMLWC, &
             'IN FILE '''//TRIM(file)//'''')
     ELSE
        CALL RGMSG(substr, RGMLWC, &
             'IN FILE WITH ID ',lncid,' ')
     END IF
     CALL RGMSG(substr, RGMLWC, &
          'EXISTS ALREADY WITH NUMBER ',zattid , ' ')
     CALL RGMSG(substr, RGMLWC, &
          'AND TYPE ', zxtype,' !')
     lex = .true.
  ELSE                           ! attribute does not exist
     CALL RGMSG(substr, RGMLI, &
          'WRITING NEW ATTRIBUTE '''//TRIM(att%name)//'''')
     CALL RGMSG(substr, RGMLIC, &
          'OF VARIABLE with ID ',att%varid,' ')
     IF (PRESENT(file)) THEN
        CALL RGMSG(substr, RGMLIC, &
             'TO FILE '''//TRIM(file)//''' !')
     ELSE
        CALL RGMSG(substr, RGMLIC, &
             'TO FILE WITH ID',lncid,' !')
     END IF
     lex = .false.
  END IF

  IF ((.NOT.lex).OR.(lclobber.AND.(zxtype == att%xtype))) THEN
     IF (lclobber.AND.(zxtype == att%xtype)) THEN
        CALL RGMSG(substr, RGMLIC, 'OVERWRITE !')
     END IF
     !
     alen = PRODUCT(att%dat%dim)
     !

     IF (TRIM(att%name) == '_FillValue') THEN
        CALL RGMSG(substr, RGMLW, &
             '------------ === *** !!! *** === -----------------')
        CALL RGMSG(substr, RGMLWC, &
             'ATTRIBUTE ''_FillValue'' CURRENTLY NOT SUPPORTED !')
        CALL RGMSG(substr, RGMLWC, &
             '------------ === *** !!! *** === -----------------')
     ELSE

        SELECT CASE(att%xtype)
        CASE(NF90_BYTE)
           CALL NFERR(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vb(1:alen))                     &
                ,4)
        CASE(NF90_SHORT)
           CALL NFERR(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vi(1:alen))                     &
                ,5)
        CASE(NF90_CHAR)
           IF (alen >= 1) THEN
              CALL NFERR(substr                             &
                   ,nf90_put_att(lncid, att%varid, att%name &
                   ,string(att%dat%vc))                     &
                   ,6)
           ELSE
              CALL NFERR(substr                             &
                   ,nf90_put_att(lncid, att%varid, att%name &
                   ,' ')                                    &
                   ,6)
           END IF
        CASE(NF90_INT)
           CALL NFERR(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vi(1:alen))                     &
                ,7)
        CASE(NF90_FLOAT)
           CALL NFERR(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vr(1:alen))                     &
                ,8)
        CASE(NF90_DOUBLE)
           CALL NFERR(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vd(1:alen))                     &
                ,9)
        CASE DEFAULT
           CALL RGMSG(substr, RGMLW, 'ATTRIBUTE '''//TRIM(att%name)//'''' )
           CALL RGMSG(substr, RGMLWC, 'OF VARIABLE WITH ID ',att%varid,' ')
           CALL RGMSG(substr, RGMLWC, 'IS OF UNKNOWN TYPE ',att%xtype,' ')
           CALL RGMSG(substr, RGMLWC, 'AND CANNOT BE WRITTEN !')
        END SELECT
  END IF ! ATTRIBUTE _FillValue

  END IF                         ! attribute exists ?

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL NFERR(substr       &
          ,nf90_close(lncid) &
          ,16)
  END IF

END SUBROUTINE EXPORT_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IMPORT_NCVAR(var, ustep, varname, varid, file, ncid, setuid &
                           , pstart, pcount, hdimids)

  ! read variable structure of type t_ncvar from netcdf file

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar),   INTENT(INOUT)           :: var      ! variable
  INTEGER,          INTENT(IN),  OPTIONAL   :: ustep    ! step along unlim. DIM
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL   :: varname  ! variable name
  INTEGER,          INTENT(IN),  OPTIONAL   :: varid    ! variable ID
  INTEGER,          INTENT(IN),  OPTIONAL   :: ncid     ! netCDF file ID
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL   :: file     ! filename
  INTEGER,          INTENT(IN),  OPTIONAL   :: setuid   ! set this as unlim. ID
  INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: pstart
  INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: pcount
  INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: hdimids 
  ! if (pstart(i) + pcount(i)) > var%dim(i)%len a "modulo" axis is assumed and
  ! a two-stage method for reading in from the netcdf is applied

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER        :: substr = 'IMPORT_NCVAR'
  INTEGER                            :: lncid    ! netCDF file-ID
  INTEGER, DIMENSION(:), ALLOCATABLE :: zdimvec  ! local dim. ID vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: zdimlen  ! local dim. lengths
  INTEGER                            :: noudims  ! number of non-ulim. dim.s
  INTEGER                            :: i,j      ! counter
  INTEGER, DIMENSION(:), ALLOCATABLE :: start, zcount, stride, map
  INTEGER                            :: status
  INTEGER                            :: ulen     ! length of unlimited ID
  INTEGER                            :: uid      ! ID of unlimited dimension
  INTEGER, DIMENSION(:), ALLOCATABLE :: fcount    ! count number written in first pass
  LOGICAL, DIMENSION(:), ALLOCATABLE :: laxmodulo ! axis is modulo, apply two-stage read-in
  TYPE (t_narray)                    :: varr1      ! variable
  TYPE (t_narray)                    :: varr2      ! variable
  INTEGER                            :: map2, map3, map4
  INTEGER                            :: map2a, map3a, map4a
  INTEGER                            :: map2b, map3b, map4b
  INTEGER, DIMENSION(:), ALLOCATABLE :: zdimlenA  ! local dim. lengths
  INTEGER, DIMENSION(:), ALLOCATABLE :: zdimlen2  ! local dim. lengths
  INTEGER, DIMENSION(5)              :: idx, jdim
  INTEGER                            :: i1, i2, i3,i4, n, n2, ix
  INTEGER, DIMENSION(5)              :: dimall, dim1, dim2
  INTEGER                            :: hdim1, hdim2
  INTEGER, DIMENSION(2)              :: hdids = -99

  ! CHECK CALL OF SUBROUTINE
  IF ((.NOT.PRESENT(varname)).AND.(.NOT.PRESENT(varid))) THEN
     CALL RGMSG(substr, RGMLE, &
          'NAME OR ID OF VARIABLE MUST BE SPECIFIED !')
  END IF
  !
  IF ((.NOT.PRESENT(ncid)).AND.(.NOT.PRESENT(file))) THEN
     CALL RGMSG(substr, RGMLE, &
          'netCDF-FILE-ID OR FILENAME MUST BE SPECIFIED !')
  END IF
  IF (PRESENT(ustep)) THEN
     IF (ustep < 1) THEN
        CALL RGMSG(substr, RGMLE, 'U-ID STEP MUST BE > 0!')
     END IF
     var%ustep = ustep
  END IF

  IF (PRESENT(hdimids)) THEN
     hdids = hdimids
  ELSE
     IF (PRESENT(pstart) .OR. PRESENT(pcount)) THEN
        CALL RGMSG(substr, RGMLE, &
             'pstart / pcount present also requires hdimids!')
     END IF
     hdids = -99
  END IF

  ! INIT
  CALL INIT_NCVAR(var)

  ! OPEN FILE
  IF (PRESENT(ncid)) THEN
     lncid = ncid
  ELSE
     ! OPEN FILE
     CALL NFERR(substr                              &
          ,nf90_open(TRIM(file),NF90_NOWRITE,lncid) & 
          ,1)
  END IF

  IF (PRESENT(varname)) THEN
     ! GET VAR-ID from VARNAME
     var%name = TRIM(varname)
     IF (TRIM(var%name) /= '') THEN
        CALL NFERR(substr                                 &
             ,nf90_inq_varid(lncid,TRIM(var%name),var%id) &
             ,2)
     ELSE
        write (0,*) 'INVALID VARIABLE NAME !', TRIM(varname)
        CALL RGMSG(substr, RGMLE, 'INVALID VARIABLE NAME ! '//TRIM(varname))
     END IF
  ELSE
     var%id = varid
     ! GET VAR-NAME from VAR-ID
     CALL NFERR(substr                                          &
          ,nf90_Inquire_Variable(lncid, var%id, name=var%name)  &
          ,3)
  END IF
!  write (0,*) '**********************************************'
!  write (0,*) 'IMPORT_NCVAR VARIABLE NAME: ', TRIM(var%name),PRESENT(hdimids)
!  write (0,*) '**********************************************'
  ! GET INFOs
  CALL NFERR(substr                                         &
       ,nf90_Inquire_Variable(lncid, var%id                 & 
       ,xtype=var%xtype, ndims=var%ndims, natts=var%natts ) &
       ,4)

  ! CHECK FOR UNLIMITED DIM ID
  CALL NFERR(substr                                   &
             ,nf90_Inquire(lncid, unlimitedDimID=uid) &
             ,5)
  IF (uid == -1) uid = NULL_DIMID
  ! SET optionally unlimitedID, if not in file
  IF ((uid == NULL_DIMID).AND.PRESENT(setuid)) THEN
     uid = setuid
  END IF

  ! ALLOCATE MEMORY AND GET DIMENSION INFO
  IF (var%ndims > 0) THEN
     ALLOCATE(zdimvec(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,1)
     CALL NFERR(substr                          &
          ,nf90_Inquire_Variable(lncid, var%id  & 
          ,dimids=zdimvec)                      &
          ,6)
     !
     noudims = var%ndims
     DO i=1, var%ndims
        IF ((zdimvec(i) == uid).AND.(uid /= NULL_DIMID)) THEN
           var%uid = uid
           noudims = noudims - 1
        END IF
     END DO
     !
     hdim1 = -99
     hdim2 = -99
     ALLOCATE(var%dim(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,2)
     DO i=1, var%ndims
        CALL IMPORT_NCDIM(var%dim(i), dimid=zdimvec(i), ncid=lncid)
        IF (zdimvec(i) == var%uid) THEN
           ulen = var%dim(i)%len
           var%dim(i)%fuid = .true.
           IF (PRESENT(ustep)) THEN
              var%dim(i)%len = 1     ! ONLY 1 SLICE ALONG UNLIM-DIM
           END IF
        END IF
        IF (var%dim(i)%id==hdids(1)) hdim1 = i
        IF (var%dim(i)%id==hdids(2)) hdim2 = i
     END DO
     IF (hdim1 <0) hdim1 = 1
     IF (hdim2 <0) hdim2 = 2
     DEALLOCATE(zdimvec, STAT=status)
     CALL ERRMSG(substr,status,3)
  END IF

  ! ALLOCATE MEMORY AND GET ATTRIBUTE INFO
  IF (var%natts > 0) THEN
     ALLOCATE(var%att(var%natts), STAT=status)
     CALL ERRMSG(substr,status,4)
     DO i=1, var%natts   ! LOOP OVER ATTRIBUTES
        CALL IMPORT_NCATT(var%att(i),varname=TRIM(var%name) &
                          ,attid =i,ncid=lncid)
!    ! ALTERNATIVE:
!        CALL IMPORT_NCATT(var%att(i),varid=var%id) &
!                          ,attid =i,ncid=lncid)
     END DO  ! LOOP OVER ATTRIBUTES
  END IF

  ! ALLOCATE MEMORY AND GET DATA
  ! (IF USTEP IS PRESENT, then ONLY ONE SLICE ALONG THE UNLIMITED DIM)
  ! VECTORS
!!$  dim = 1
  ! initialise modulo marker as .FALSE., cover all number of dimensions (ELSE branch!)
  IF (var%ndims > 0) THEN        ! DIM > 0
     ALLOCATE(start(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,5)
     ALLOCATE(zcount(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,6)
     ALLOCATE(stride(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,7)
     ALLOCATE(map(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,8)

     ALLOCATE(laxmodulo(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,9)
     laxmodulo = .FALSE.
     ALLOCATE(fcount(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,10)
     fcount = 0
!CDIR NOVECTOR
     DO i=1, var%ndims
        ! set indices for pstart / pcount
        IF (i==hdim1) ix=1
        IF (i==hdim2) ix=2
        
        IF ((var%dim(i)%id == var%uid).AND.(var%uid /= NULL_DIMID).AND.   &
             PRESENT(ustep)) THEN
           ! UNLIMITED DIMENSION
           IF (ustep <= ulen) THEN
              start(i)  = ustep
              zcount(i) = 1
              stride(i) = 1
           ELSE
              CALL RGMSG(substr, RGMLE, &
                   'REQUESTET U-ID STEP (',ustep,')', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'IS LARGER THAN AVAILABLE (',ulen,') !')
           END IF
        ELSE
           ! LIMITED DIMENSION
           IF (PRESENT(pstart) .AND. (i == hdim1 .OR. i == hdim2) ) THEN
              start(i) = pstart(ix)
           ELSE
              start(i)  = 1 
           ENDIF
           IF (PRESENT(pcount) .AND. (i == hdim1 .OR. i == hdim2) ) THEN
              IF (pcount(ix) > 0) THEN
                 zcount(i) = pcount(ix)
              ELSE
                 zcount(i) = var%dim(i)%len
              ENDIF
           ELSE
              zcount(i) = var%dim(i)%len
           ENDIF
           stride(i) = 1
        END IF
!!$        dim = dim * zcount(i)
        IF (i == 1) THEN
           map(i) = 1
        ELSE
           map(i) = map(i-1)*zcount(i-1)
        END IF
        ! test for two-stage read-in of modulo axis
        ! restrict zcount to reasonable length for first step
        IF ((i == hdim1 .OR. i == hdim2) &
             .AND. .NOT.(var%dim(i)%id == var%uid)) THEN
           IF ((start(i) + zcount(i) - 1) > var%dim(i)%len) THEN
              laxmodulo(i) = .TRUE.
              fcount(i) = zcount(i)
              zcount(i) = var%dim(i)%len - start(i) + 1
           END IF
        END IF
        IF (var%ndims >=2) THEN
           IF (laxmodulo(hdim1) .AND. laxmodulo(hdim2)) THEN
              CALL RGMSG(substr, RGMLE, 'GRID CAN NOT HAVE 2 MODULO AXES.')
           END IF
        END IF
     END DO
  ELSE
     ALLOCATE(start(1), STAT=status)
     CALL ERRMSG(substr,status,9)
     ALLOCATE(zcount(1), STAT=status)
     CALL ERRMSG(substr,status,10)
     ALLOCATE(stride(1), STAT=status)
     CALL ERRMSG(substr,status,11)
     ALLOCATE(map(1), STAT=status)
     CALL ERRMSG(substr,status,12)
     start(1) = 1
     zcount(1) = 1
     stride(1) = 1
     map(1) = 1
  END IF

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
     IF (PRESENT(pcount) .AND. (i == hdim1 .OR. i == hdim2)) THEN
        ! set indices for pstart / pcount
        IF (i==hdim1) ix=1
        IF (i==hdim2) ix=2
        IF (pcount(ix) > 0) THEN
           zdimlen(i) = zcount(i)
        ELSE
           zdimlen(i) = var%dim(i)%len  ! = 1 for UNLIM-DIM-ID
        END IF
     ELSE
        zdimlen(i) = var%dim(i)%len     ! = 1 for UNLIM-DIM-ID
     END IF
  END DO
  j = SIZE(zdimlen)

  ! ALLOCATE DATA MEMORY AND GET DATA
  SELECT CASE(var%xtype)
  CASE(NF90_BYTE)
     CALL INIT_NARRAY(varr1, j, zdimlen, VTYPE_BYTE)
     CALL NFERR(substr                                             &
                ,nf90_get_var(lncid, var%id, varr1%vb           &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,7)
  CASE(NF90_SHORT)
     CALL INIT_NARRAY(varr1, j, zdimlen, VTYPE_INT)
     CALL NFERR(substr                                             &
                ,nf90_get_var(lncid, var%id, varr1%vi            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,8)     
  CASE(NF90_CHAR)
     CALL INIT_NARRAY(varr1, j, zdimlen, VTYPE_CHAR)
     CALL NFERR(substr                                             &
                ,nf90_get_var(lncid, var%id, varr1%vc            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,9)           
  CASE(NF90_INT)
     CALL INIT_NARRAY(varr1, j, zdimlen, VTYPE_INT)
     CALL NFERR(substr                                             &
                ,nf90_get_var(lncid, var%id, varr1%vi            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,10)
  CASE(NF90_FLOAT)
     CALL INIT_NARRAY(varr1, j, zdimlen, VTYPE_REAL)
     CALL NFERR(substr                                             &
                ,nf90_get_var(lncid, var%id, varr1%vr            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,11)
  CASE(NF90_DOUBLE)
     CALL INIT_NARRAY(varr1, j, zdimlen, VTYPE_DOUBLE)
     CALL NFERR(substr                                             &
                ,nf90_get_var(lncid, var%id, varr1%vd            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,12)
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, &
          'UNKNOWN VARIABLE TYPE ',var%xtype,' ', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'OF VARIABLE '''//TRIM(var%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(file)//''' !')
  END SELECT

  ! Second step of two-stage read-in for wrapped-around modulo axes
  IF (ANY(laxmodulo)) THEN
     DO ix = 1,MIN(var%ndims,2)
        IF (ix == 1) i = hdim1
        IF (ix == 2) i = hdim2
        IF (laxmodulo(i)) THEN
           start(i)  = 1
           ! adjust the map according to the already read data
           zcount(i) = fcount(i) - zcount(i)
           
           DO i1 = 1, var%ndims
              IF (i1 == 1) THEN
                 map(i1) = 1
              ELSE
                 map(i1) = map(i1-1)*zcount(i1-1)
              END IF
           END DO
        END IF
     END DO
     ALLOCATE(zdimlen2(var%ndims))
     DO i=1, var%ndims
        IF (PRESENT(pcount) .AND. (i==hdim1 .OR. i==hdim2)) THEN
           ! set indices for pstart / pcount
           IF (i==hdim1) ix=1
           IF (i==hdim2) ix=2
           IF (pcount(ix) > 0) THEN
              zdimlen2(i) = zcount(i)
           ELSE
              zdimlen2(i) = var%dim(i)%len  ! = 1 for UNLIM-DIM-ID
           END IF
        ELSE
           zdimlen2(i) = var%dim(i)%len     ! = 1 for UNLIM-DIM-ID
        END IF
     END DO
     j = SIZE(zdimlen2)

     SELECT CASE(var%xtype)
     CASE(NF90_BYTE)
        CALL INIT_NARRAY(varr2, j, zdimlen2, VTYPE_BYTE)
        CALL NFERR(substr                                                                     &
           &      ,nf90_get_var(lncid, var%id, varr2%vb                                     &
           &                   ,start=start, count=zcount, stride=stride,map=map)             &
           &      ,13)
     CASE(NF90_SHORT)
        CALL INIT_NARRAY(varr2, j, zdimlen2, VTYPE_INT)

        CALL NFERR(substr                                                                    & 
           &      ,nf90_get_var(lncid, var%id, var%dat%vi                                     &
           &                   ,start=start, count=zcount, stride=stride,map=map)             &
           &      ,14)
     CASE(NF90_CHAR)
        CALL INIT_NARRAY(varr2, j, zdimlen2, VTYPE_CHAR)
        CALL NFERR(substr                                                                     &
           &      ,nf90_get_var(lncid, var%id, varr2%vc                                     &
           &                   ,start=start, count=zcount, stride=stride,map=map)             &
           &      ,15)
     CASE(NF90_INT)
        CALL INIT_NARRAY(varr2, j, zdimlen2, VTYPE_INT)

        CALL NFERR(substr                                                                     &
           &      ,nf90_get_var(lncid, var%id, varr2%vi                                     &
           &                   ,start=start, count=zcount, stride=stride, map=map)            &
           &      ,16)

     CASE(NF90_FLOAT)
        CALL INIT_NARRAY(varr2, j, zdimlen2, VTYPE_REAL)
        CALL NFERR(substr                                                                     &
           &      ,nf90_get_var(lncid, var%id, varr2%vr                                     &
           &                   ,start=start, count=zcount, stride=stride, map=map)            &
           &      ,17)
     CASE(NF90_DOUBLE)
        CALL INIT_NARRAY(varr2, j, zdimlen2, VTYPE_DOUBLE)
        CALL NFERR(substr                                                                     &
           &      ,nf90_get_var(lncid, var%id, varr2%vd                                     &
           &                   ,start=start, count=zcount, stride=stride, map=map)            &
           &      ,18)
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE,  'UNKNOWN VARIABLE TYPE ',var%xtype,' ', .false.)
        CALL RGMSG(substr, RGMLEC, 'OF VARIABLE '''//TRIM(var%name)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, 'IN FILE '''//TRIM(file)//''' !')
     END SELECT

     dimall(:) = 1
     dim1(:)  = 1
     dim2(:)  = 1
     ALLOCATE(zdimlenA(var%ndims))
     DO i=1, var%ndims
        IF (PRESENT(pcount) .AND. (i==hdim1 .OR. i == hdim2)) THEN
           ! set indices for pstart / pcount
           IF (i==hdim1) ix=1
           IF (i==hdim2) ix=2
           IF (pcount(ix) > 0) THEN
              zdimlenA(i) = pcount(ix)
           ELSE
              zdimlenA(i) = var%dim(i)%len  ! = 1 for UNLIM-DIM-ID
           END IF
        ELSE
           zdimlenA(i) = var%dim(i)%len     ! = 1 for UNLIM-DIM-ID
        END IF
     END DO
     j = SIZE(zdimlenA)
     
     SELECT CASE(var%xtype)
     CASE(NF90_BYTE)
        CALL INIT_NARRAY(var%dat, j, zdimlenA, VTYPE_BYTE)
     CASE(NF90_SHORT)
        CALL INIT_NARRAY(var%dat, j, zdimlenA, VTYPE_INT)
     CASE(NF90_CHAR)
        CALL INIT_NARRAY(var%dat, j, zdimlenA, VTYPE_CHAR)
     CASE(NF90_INT)
        CALL INIT_NARRAY(var%dat, j, zdimlenA, VTYPE_INT)
     CASE(NF90_FLOAT)
        CALL INIT_NARRAY(var%dat, j, zdimlenA, VTYPE_REAL)
     CASE(NF90_DOUBLE)
        CALL INIT_NARRAY(var%dat, j, zdimlenA, VTYPE_DOUBLE)
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE,                              &
             'N-ARRAY OF UNRECOGNIZED TYPE CANNOT BE COPIED !' )          
     END SELECT

     dimall(1:SIZE(zdimlenA)) = zdimlenA(:)
     dim1(1:SIZE(zdimlen))    = zdimlen(:)
     dim2(1:SIZE(zdimlen2))   = zdimlen2(:)
 
     ! IMPLICIT ASSUMPTION TIME IS ALWaYS LAST (5th dim for 4D data and can 
     ! be ignored)
     DO i4 = 1, dimall(4)
        idx(4) = i4
        map4= dimall(3)*dimall(2)*dimall(1)
        map4a = dimall(3)*dim1(2)*dim1(1)
        map4b = dimall(3)*dim2(2)*dim2(1)
        DO i3 = 1, dimall(3)
           idx(3) = i3
           map3= dimall(2)*dimall(1)
           map3a = dim1(2)*dim1(1)
           map3b = dim2(2)*dim2(1)
           DO i2 = 1, dimall(2)
              idx(2) = i2
              map2 = dimall(1)
              map2a = dim1(1)
              map2b = dim2(1)
              DO i1 = 1, dimall(1)
                 idx(1) = i1
                 n = (idx(4)-1) * map4 + (idx(3)-1) * map3 + &
                      (idx(2)-1) *map2 + idx(1)
                 IF (idx(hdim1) > dim1(hdim1)) THEN
                    ! second chunck
                    jdim(:) = 0
                    jdim(hdim1) = dim1(hdim1)
                    n2 =(idx(4)-jdim(4)-1)*map4b +(idx(3)-jdim(3)-1)*map3b &
                         +(idx(2)-jdim(2)-1)*map2b + (idx(1)-jdim(1))
                    SELECT CASE(var%dat%type)
                    CASE(VTYPE_INT)
                       var%dat%vi(n) = varr2%vi(n2)
                    CASE(VTYPE_REAL)
                       var%dat%vr(n) = varr2%vr(n2)
                    CASE(VTYPE_DOUBLE)
                       var%dat%vd(n) = varr2%vd(n2)
                    CASE(VTYPE_CHAR)
                       var%dat%vc(n) = varr2%vc(n2)
                    CASE(VTYPE_BYTE)
                       var%dat%vb(n) = varr2%vb(n2)
                    END SELECT
                 ELSE
                    IF (idx(hdim2) > dim1(hdim2)) THEN
                       jdim(:) = 0
                       jdim(hdim2) = dim1(hdim2)
                       n2 =(idx(4)-jdim(4)-1)*map4b +(idx(3)-jdim(3)-1)*map3b &
                            +(idx(2)-jdim(2)-dim1(2)-1)*map2b + idx(1)-jdim(1)
                       var%dat%vd(n) = varr2%vd(n2)
                       SELECT CASE(var%dat%type)
                       CASE(VTYPE_INT)
                          var%dat%vi(n) = varr2%vi(n2)
                       CASE(VTYPE_REAL)
                          var%dat%vr(n) = varr2%vr(n2)
                       CASE(VTYPE_DOUBLE)
                          var%dat%vd(n) = varr2%vd(n2)
                       CASE(VTYPE_CHAR)
                          var%dat%vc(n) = varr2%vc(n2)
                       CASE(VTYPE_BYTE)
                          var%dat%vb(n) = varr2%vb(n2)
                       END SELECT
                    ELSE
                       n2 =(idx(4)-1)*map4a +(idx(3)-1)*map3a &
                            +(idx(2)-1)*map2a + idx(1)
                       SELECT CASE(var%dat%type)
                       CASE(VTYPE_INT)
                          var%dat%vi(n) = varr1%vi(n2)
                       CASE(VTYPE_REAL)
                          var%dat%vr(n) = varr1%vr(n2)
                       CASE(VTYPE_DOUBLE)
                          var%dat%vd(n) = varr1%vd(n2)
                       CASE(VTYPE_CHAR)
                          var%dat%vc(n) = varr1%vc(n2)
                       CASE(VTYPE_BYTE)
                          var%dat%vb(n) = varr1%vb(n2)
                       END SELECT
                    END IF
                 END IF
              END DO
           END DO
        END DO
     END DO

     ! SET CORRECT DIMENSION LENGTH
     DO i=1,var%ndims
        var%dim(i)%len = zdimlenA(i)
     END DO

     CALL INIT_NARRAY(varr2)
     DEALLOCATE(zdimlen2, zdimlenA)
  ELSE
     CALL COPY_NARRAY(var%dat, varr1)
     ! SET CORRECT DIMENSION LENGTH
     DO i=1,var%ndims
        var%dim(i)%len = zdimlen(i)
     END DO
  END IF ! ANY(laxmodulo)
  CALL INIT_NARRAY(varr1)

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL NFERR(substr       &
          ,nf90_close(lncid) &
          ,19)
  END IF

  ! CLEAN UP
  DEALLOCATE(start, zcount, stride, map, zdimlen, STAT=status)
  IF (ALLOCATED(fcount))    DEALLOCATE(fcount)
  IF (ALLOCATED(laxmodulo)) DEALLOCATE(laxmodulo)
  CALL ERRMSG(substr,status,15)

END SUBROUTINE IMPORT_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_NCVAR(var, file, ncid)

  ! write variable  structure of type t_ncvar to netcdf file

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar),   INTENT(INOUT)           :: var     ! variable
  CHARACTER(LEN=*), INTENT(IN),   OPTIONAL  :: file    ! netCDF filename
  INTEGER,          INTENT(IN),   OPTIONAL  :: ncid    ! netCDF file ID

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER        :: substr = 'EXPORT_NCVAR'
  LOGICAL                            :: lex      ! file exists ?
  INTEGER                            :: lncid    ! local netCDF file-ID
  INTEGER                            :: i        ! counter
  INTEGER, DIMENSION(:), ALLOCATABLE :: start, zcount, stride, map
  INTEGER                            :: dim      ! 1D dimension of variable
  INTEGER                            :: status   ! netcdf status
  INTEGER                            :: zndims   ! local number of dimensions
  INTEGER                            :: zxtype   ! local var./att. type
  INTEGER                            :: ulen     ! length of unlimited ID
  INTEGER                            :: ustep    ! step along unlimite ID
  INTEGER, DIMENSION(:), ALLOCATABLE :: zdids  ! local dimension IDs
  INTEGER, DIMENSION(:), ALLOCATABLE :: zdlen  ! local dimension lengths
  CHARACTER(LEN=GRD_MAXSTRLEN), DIMENSION(:), ALLOCATABLE :: zdnames

  ! CHECK CALL OF SUBROUTINE
  IF ((.NOT.PRESENT(ncid)).AND.(.NOT.PRESENT(file))) THEN
     CALL RGMSG(substr, RGMLE, &
          'netCDF-FILE-ID OR FILENAME MUST BE SPECIFIED !')
  END IF

  ! CHECK VARNAME, VARID
  IF (TRIM(var%name) == '') THEN
     CALL RGMSG(substr, RGMLE, 'VARIABLE IS NAMELESS !')
  END IF

  ! INIT
  ustep = 0
  ulen = 0

  IF (PRESENT(ncid)) THEN
     lncid = ncid
     ! DEFINE MODE
     CALL RGMSG(substr, RGMLI, &
          'netCDF-FILE WITH ID ',lncid,': SWITCHING TO DEFINE MODE !' )
     CALL NFERR(substr       &
          ,nf90_redef(lncid) &
          ,1)
  ELSE
     ! CHECK IF FILE EXISTS
     INQUIRE(file=TRIM(file), exist=lex)
     IF (lex) THEN   ! FILE EXISTS
        CALL RGMSG(substr, RGMLI, &
             'FILE '''//TRIM(file)//''' EXISTS ! OPENING ...')
        CALL NFERR(substr                            &
             ,nf90_open(TRIM(file),NF90_WRITE,lncid) &
             ,2)
        CALL NFERR(substr       &
             ,nf90_redef(lncid) &
             ,3)
     ELSE            ! NEW FILE
        CALL RGMSG(substr, RGMLI, &
             'FILE '''//TRIM(file)//''' DOES NOT EXIST ! NEW FILE ...')
        CALL NFERR(substr                                &
             ,nf90_create(TRIM(file),NF90_CLOBBER,lncid) &
             ,4)
     END IF
  END IF

  CALL RGMSG(substr, RGMLI, &
       'EXPORTING DIMENSIONS OF VARIABLE '''//TRIM(var%name)//''' ...')

  ! DEFINE REQUIRED DIMENSIONS FOR VARIABLE AND SET UID
  var%uid = NULL_DIMID
  IF (var%ndims > 0) THEN
     DO i=1, var%ndims   ! LOOP OVER DIMENSIONS
        CALL EXPORT_NCDIM(var%dim(i), ncid=lncid)
        IF (var%dim(i)%fuid) THEN
           var%uid = var%dim(i)%id
           ustep = var%ustep
        END IF
     END DO
  END IF

  CALL RGMSG(substr, RGMLI, &
       'EXPORTING VARIABLE '''//TRIM(var%name)//''' ...')
  !
  ! DEFINE VARIABLE
  status = nf90_inq_varid(lncid,TRIM(var%name),var%id)
  IF (status == NF90_NOERR) THEN  ! VAR exists
     CALL RGMSG(substr, RGMLIC, ' ... VARIABLE EXISTS')
     ! CHECK FOR COMPATIBILITY
     ! type and number of dimensions
     CALL NFERR(substr                &
          ,nf90_Inquire_Variable(lncid,var%id  &
          , xtype = zxtype ,ndims = zndims)    &
          ,5)
     ! type
     IF (zxtype == var%xtype) THEN
     CALL RGMSG(substr, RGMLIC, ' ... IS OF CORRECT TYPE')
     ELSE
        CALL RGMSG(substr, RGMLIC, ' ... IS OF DIFFERENT TYPE')
        CALL RGMSG(substr, RGMLE, 'VARIABLE '''//TRIM(var%name)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, 'IN FILE '''//TRIM(file)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, 'EXISTS ALREADY WITH TYPE ',zxtype,' !')
     END IF
     ! number of dimensions
     IF (zndims == var%ndims) THEN
        CALL RGMSG(substr, RGMLIC, ' ... HAS CORRECT NUMBER OF DIMENSIONS')
     ELSE
        CALL RGMSG(substr, RGMLIC, ' ... HAS DIFFERENT NUMBER OF DIMENSIONS')
        CALL RGMSG(substr, RGMLE, 'VARIABLE '''//TRIM(var%name)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, 'IN FILE '''//TRIM(file)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, 'EXISTS ALREADY ',zndims,' DIMENSIONS !')
     END IF
     !
     ALLOCATE(zdids(zndims), zdlen(zndims), zdnames(zndims), STAT=status)
     CALL ERRMSG(substr,status,1)
     CALL NFERR(substr                  &
          ,nf90_Inquire_Variable(Lncid,var%id    &
          ,dimids = zdids)                       &
          ,6)
     DO i=1, zndims ! LOOP OVER DIMENSIONS
        IF (zdids(i) /= var%dim(i)%id) THEN
           CALL RGMSG(substr, RGMLIC, ' ... HAS DIFFERENT DIMENSION ID')
           CALL RGMSG(substr, RGMLE,  &
                'VARIABLE '''//TRIM(var%name)//'''', .false.)
           CALL RGMSG(substr, RGMLEC, &
                'IN FILE '''//TRIM(file)//'''', .false.)
           CALL RGMSG(substr, RGMLEC, &
                'EXISTS ALREADY WITH DIM-ID',zdids(i),' ', .false.)
           CALL RGMSG(substr, RGMLEC, 'AT POSITION ',i,' !')
        END IF
        CALL RGMSG(substr, RGMLIC, &
             ' ... HAS CORRECT DIMENSION ID AT POSITION ',i,' ')
        !
        CALL NFERR(substr                           &
             ,nf90_Inquire_Dimension(lncid,zdids(i) &
             ,name = zdnames(i), len = zdlen(i))    &
             ,6)
        ! dimension length (if not unlim-dim)
        IF (var%dim(i)%fuid) THEN ! UNLIMITED DIM.
           ulen = zdlen(i)
           ! CHECK IF ustep < ulen
           IF ((ustep < ulen).AND.(ustep > 1)) THEN
              CALL RGMSG(substr, RGMLW,  'REQUESTET U-ID STEP (',ustep,')')
              CALL RGMSG(substr, RGMLWC, 'IS LOWER THAN OR EQUAL TO ONE')
              CALL RGMSG(substr, RGMLWC, &
                   'WHICH IS ALREADY PRESENT (',ulen,') !')
              CALL RGMSG(substr, RGMLWC, 'DATA WILL BE OVERWRITTEN !')
           END IF
        ELSE
           IF (zdlen(i) == var%dim(i)%len) THEN
              CALL RGMSG(substr, RGMLIC, &
                   ' ... HAS CORRECT DIMENSION LENGHT AT POSITION ',i,' ')
           ELSE
              CALL RGMSG(substr, RGMLIC, &
                   ' ... HAS DIFFERENT DIMENSION LENGTH')
              CALL RGMSG(substr, RGMLE, &
                   'VARIABLE '''//TRIM(var%name)//'''', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'IN FILE '''//TRIM(file)//'''', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'EXISTS ALREADY WITH DIM-LEN ',zdlen(i),' ', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'AT POSITION ',i,' !')
           END IF
           ! dimension name
           IF (TRIM(zdnames(i)) == TRIM(var%dim(i)%name)) THEN
              CALL RGMSG(substr, RGMLIC, &
                   ' ... HAS CORRECT DIMENSION NAME AT POSITION ',i,' ')
           ELSE
              CALL RGMSG(substr, RGMLIC, &
                   ' ... HAS DIFFERENT DIMENSION NAME')
              CALL RGMSG(substr, RGMLE, &
                   'VARIABLE '''//TRIM(var%name)//'''', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'IN FILE '''//TRIM(file)//'''', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'EXISTS ALREADY WITH DIM-NAME '''//TRIM(zdnames(i))//'''',&
                   .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'AT POSITION ',i,' !')
           END IF
        END IF ! (non-unlim-dim)
     END DO    ! LOOP OVER DIMENSIONS
     DEALLOCATE(zdids, zdlen, zdnames, STAT=status)
     CALL ERRMSG(substr,status,2)
  ELSE                            ! VAR does not exist
     CALL RGMSG(substr, RGMLIC, ' ... NEW')
     IF (var%ndims > 0) THEN
        ALLOCATE(zdids(var%ndims), STAT=status)
        CALL ERRMSG(substr,status,3)
        DO i=1, var%ndims
           zdids(i) = var%dim(i)%id
        END DO
        CALL NFERR(substr                       &
             ,nf90_def_var(lncid,TRIM(var%name) &
             ,var%xtype, zdids, var%id)         &
             ,7)
        DEALLOCATE(zdids, STAT=status)
        CALL ERRMSG(substr,status,4)
     ELSE
        CALL NFERR(substr                       &
             ,nf90_def_var(lncid,TRIM(var%name) &
             ,xtype=var%xtype, varid=var%id)    &
             ,8)
     END IF
  END IF    ! VAR exists ?

  CALL RGMSG(substr, RGMLI, &
       'EXPORTING ATTRIBUTES OF VARIABLE '''//TRIM(var%name)//''' ...')

  ! DEFINE ATTRIBUTES
  IF (var%natts > 0) THEN
     DO i=1, var%natts  ! LOOP OVER ATTRIBUTES
        var%att(i)%varid = var%id
        CALL EXPORT_NCATT(var%att(i), ncid=lncid)
     END DO
  END IF

  ! SWITCH FROM DEFINITION TO DATA MODE
  CALL NFERR(substr,nf90_enddef(lncid),9)

  CALL RGMSG(substr, RGMLI, &
       'EXPORTING DATA OF VARIABLE '''//TRIM(var%name)//''' ...')

  ! ALLOCATE VECTORS
  dim = 1
  IF (var%ndims > 0) THEN        ! DIM > 0
     ALLOCATE(start(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,5)
     ALLOCATE(zcount(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,6)
     ALLOCATE(stride(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,7)
     ALLOCATE(map(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,8)
     DO i=1, var%ndims
        IF ((var%dim(i)%id == var%uid).AND.(var%uid /= NULL_DIMID).AND. &
             (ustep > 0)) THEN  ! UNLIMITED DIMENSION
           start(i)  = ustep
           zcount(i) = 1
           stride(i) = 1
        ELSE
           ! LIMITED DIMENSION OR (USTEP <= 0) (-> COMPLETE VARIABLE)
           start(i)  = 1
           zcount(i) = var%dim(i)%len
           stride(i) = 1
        END IF
        dim = dim * zcount(i)
        IF (i == 1) THEN
           map(i) = 1
        ELSE
           map(i) = map(i-1)*zcount(i-1)
        END IF
     END DO
  ELSE
     ALLOCATE(start(1), STAT=status)
     CALL ERRMSG(substr,status,9)
     ALLOCATE(zcount(1), STAT=status)
     CALL ERRMSG(substr,status,10)
     ALLOCATE(stride(1), STAT=status)
     CALL ERRMSG(substr,status,11)
     ALLOCATE(map(1), STAT=status)
     CALL ERRMSG(substr,status,12)
     start(1)  = 1
     zcount(1) = 1
     stride(1) = 1
     map(1)    = 1
  END IF

  ! WRITE DATA
  SELECT CASE (var%xtype)
  CASE(NF90_BYTE)
     CALL NFERR(substr                                             &
                ,nf90_put_var(lncid, var%id, var%dat%vb            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,10)
  CASE(NF90_SHORT)
     CALL NFERR(substr                                             &
                ,nf90_put_var(lncid, var%id, var%dat%vi            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,11)
  CASE(NF90_CHAR)
     CALL NFERR(substr                                             &
                ,nf90_put_var(lncid, var%id, var%dat%vc            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,12)
  CASE(NF90_INT)
     CALL NFERR(substr                                             &
                ,nf90_put_var(lncid, var%id, var%dat%vi            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,13)
  CASE(NF90_FLOAT)
     CALL NFERR(substr                                             &
                ,nf90_put_var(lncid, var%id, var%dat%vr            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,14)
  CASE(NF90_DOUBLE)
     CALL NFERR(substr                                             &
                ,nf90_put_var(lncid, var%id, var%dat%vd            &
                         ,start=start, count=zcount, stride=stride &
                         ,map = map)                               &
                         ,15)
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, &
          'UNKNOWN VARIABLE TYPE ',var%xtype,' ', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'OF VARIABLE '''//TRIM(var%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'CANNOT BE EXPORTED TO FILE '''//TRIM(file)//''' !')
  END SELECT

  ! CLOSE NETCDF FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL NFERR(substr,nf90_close(lncid),16)
  END IF

  ! CLEAN UP
  DEALLOCATE(start, zcount, stride, map, STAT=status)
  CALL ERRMSG(substr,status,13)

  CALL RGMSG(substr, RGMLIC, &
       '... END EXPORTING VRIABLE '''//TRIM(var%name)//''' !')

END SUBROUTINE EXPORT_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE ADD_NCATT(var, name, replace, vs, vr, vd, vi, vb)

  ! add/replace attribute of name 'name' to variable of type t_ncvar

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar),            INTENT(INOUT)          :: var
  CHARACTER(LEN=*),          INTENT(IN)             :: name
  LOGICAL,                   INTENT(IN),   OPTIONAL :: replace
  CHARACTER(LEN=*),          INTENT(IN),   OPTIONAL :: vs
  REAL   (SP), DIMENSION(:), INTENT(IN),   OPTIONAL :: vr
  REAL   (DP), DIMENSION(:), INTENT(IN),   OPTIONAL :: vd
  INTEGER(I8), DIMENSION(:), INTENT(IN),   OPTIONAL :: vi
  INTEGER(I4), DIMENSION(:), INTENT(IN),   OPTIONAL :: vb

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER               :: substr = 'ADD_NCATT'
  INTEGER                                   :: i
  INTEGER                                   :: natts  ! number of attributes
  INTEGER                                   :: num    ! number of new attribute
  TYPE (t_ncatt), DIMENSION(:), ALLOCATABLE :: att  ! for saving
  INTEGER                                   :: status
  LOGICAL                                   :: lrpl   ! replace ?
  INTEGER                                   :: narg
  INTEGER                                   :: vtype
  INTEGER                                   :: nf90type

  ! CHECK SUBROUTINE CALL
  IF (PRESENT(replace)) THEN
     lrpl = replace
  ELSE
     lrpl = .false.
  END IF

  ! NUMBER OF ARGUMENTS
  narg = 0
  IF (PRESENT(vs)) THEN
     vtype = VTYPE_CHAR
     nf90type = NF90_CHAR
     narg = narg + 1
  END IF
  IF (PRESENT(vi)) THEN
     vtype = VTYPE_INT
     nf90type = NF90_INT
     narg = narg + 1
  END IF
  IF (PRESENT(vr)) THEN
     vtype = VTYPE_REAL
     nf90type = NF90_FLOAT
     narg = narg + 1
  END IF
  IF (PRESENT(vd)) THEN
     vtype = VTYPE_DOUBLE
     nf90type = NF90_DOUBLE
     narg = narg + 1
  END IF
  IF (PRESENT(vb)) THEN
     vtype = VTYPE_BYTE
     nf90type = NF90_BYTE
     narg = narg + 1
  END IF

  IF (narg /= 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'TYPE OF ARGUMENT AT SUBROUTINE CALL MUST BE ONE OF', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'INTEGER, REAL, DOUBLE PRECISION, OR STRING !')
  END IF

  ! LOOP OVER ATTRIBUTES AND CHECK, IF EXISTENT
  natts = 0
  num = 0
  DO i=1, var%natts
     natts = natts + 1
     IF (TRIM(var%att(i)%name) == TRIM(name)) THEN  ! FOUND
        num = natts
        EXIT
     END IF
  END DO

  IF (num == 0) THEN       ! NOT EXISTENT
     natts = natts + 1     ! ADD ONE
     var%natts = natts     ! UPDATE
     num = natts           ! NEW ONE AT END
     IF (natts > 1) THEN  ! NEW IS NOT THE FIRST ATTRIBUTE
        ! SAVE ATTRIBUTES
        ALLOCATE(att(natts-1), STAT=status)
        CALL ERRMSG(substr,status,1)
        DO i=1, natts-1
           CALL COPY_NCATT(att(i), var%att(i))
           CALL INIT_NCATT(var%att(i))
        END DO
        DEALLOCATE(var%att, STAT=status)
        CALL ERRMSG(substr,status,2)
        NULLIFY(var%att)
     END IF
     ALLOCATE(var%att(natts), STAT=status)  ! NEW SPACE
     CALL ERRMSG(substr,status,3)
     IF (natts > 1) THEN ! NEW IS NOT THE FIRST ATTRIBUTE
        ! RESTORE OLD ATTRIBUTES
        DO i=1, natts-1
           CALL COPY_NCATT(var%att(i), att(i))
           CALL INIT_NCATT(att(i))
        END DO
        DEALLOCATE(att, STAT=status)
        CALL ERRMSG(substr,status,4)
     END IF
  END IF  ! NOT EXISTENT

  IF ((TRIM(var%att(num)%name) == TRIM(name)).AND.(.NOT.lrpl)) THEN
     CALL RGMSG(substr, RGMLW, 'ATTRIBUTE '''//TRIM(name)//'''')
     CALL RGMSG(substr, RGMLWC,'OF VARIABLE '''//TRIM(var%name)//'''')
     CALL RGMSG(substr, RGMLWC,&
          'EXISTS ALREADY WITH NUMBER ',var%att(num)%ID,' ')
     CALL RGMSG(substr, RGMLWC,&
          'AND TYPE ',var%att(num)%xtype,' ')
     CALL RGMSG(substr, RGMLWC, 'AND SHOULD BE REPLACED !')
     RETURN
  END IF

  IF ((TRIM(var%att(num)%name) == TRIM(name)).AND.(lrpl)) THEN
     CALL RGMSG(substr, RGMLI, 'REPLACING ATTRIBUTE '''//TRIM(name)//'''')
     CALL RGMSG(substr, RGMLIC,'OF VARIABLE '''//TRIM(var%name)//''' !')

     ! DEALLOCATE MEMORY OF REPLACED ATTRIBUTE
     CALL INIT_NARRAY(var%att(num)%dat)
  END IF

  ! SET ATTRIBUTE
  var%att(num)%name  = TRIM(name)
  var%att(num)%ID    = num
  var%att(num)%xtype = nf90type
  var%att(num)%varid = var%id

  SELECT CASE(vtype)
  CASE(VTYPE_INT)
     var%att(num)%len = SIZE(vi)
     ALLOCATE(var%att(num)%dat%vi(var%att(num)%len), STAT=status)
     CALL ERRMSG(substr,status,5)
     var%att(num)%dat%vi(:) = vi(:)
     var%att(num)%dat%type = VTYPE_INT
  CASE(VTYPE_REAL)
     var%att(num)%len = SIZE(vr)
     ALLOCATE(var%att(num)%dat%vr(var%att(num)%len), STAT=status)
     CALL ERRMSG(substr,status,6)
     var%att(num)%dat%vr(:) = vr(:)
     var%att(num)%dat%type = VTYPE_REAL
  CASE(VTYPE_DOUBLE)
     var%att(num)%len = SIZE(vd)
     ALLOCATE(var%att(num)%dat%vd(var%att(num)%len), STAT=status)
     CALL ERRMSG(substr,status,7)
     var%att(num)%dat%vd(:) = vd(:)
     var%att(num)%dat%type = VTYPE_DOUBLE
  CASE(VTYPE_BYTE)
     var%att(num)%len = SIZE(vb)
     ALLOCATE(var%att(num)%dat%vb(var%att(num)%len), STAT=status)
     CALL ERRMSG(substr,status,8)
     var%att(num)%dat%vb(:) = vb(:)
     var%att(num)%dat%type = VTYPE_BYTE
  CASE(VTYPE_CHAR)
     var%att(num)%len = LEN_TRIM(vs)
     ALLOCATE(var%att(num)%dat%vc(var%att(num)%len), STAT=status)
     CALL ERRMSG(substr,status,9)
     DO i=1, LEN_TRIM(vs)
        var%att(num)%dat%vc(i) = vs(i:i)
     END DO
     var%att(num)%dat%type = VTYPE_CHAR
  !CASE(VTYPE_UNDEF)
  !CASE DEFAULT
  END SELECT

  var%att(num)%dat%n = 1
  ALLOCATE(var%att(num)%dat%dim(1), STAT=status)
  CALL ERRMSG(substr,status,10)
  var%att(num)%dat%dim(1) = var%att(num)%len

END SUBROUTINE ADD_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SCAN_NCVAR(var, file, ncid)

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar), DIMENSION(:), POINTER              :: var   ! variables
  CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: file  ! filename
  INTEGER,                      INTENT(IN), OPTIONAL :: ncid  ! netCDF ID

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SCAN_NCVAR'
  INTEGER :: lncid  ! local netCDF file ID
  INTEGER :: nvar   ! number of variables
  INTEGER :: i      ! counter
  INTEGER :: status

  ! CHECK CALL OF SUBROUTINE
  IF ((.NOT.PRESENT(ncid)).AND.(.NOT.PRESENT(file))) THEN
     CALL RGMSG(substr, RGMLE, &
          'netCDF-FILE-ID OR FILENAME MUST BE SPECIFIED !')
  END IF

  ! OPEN FILE
  IF (PRESENT(ncid)) THEN
     lncid = ncid
  ELSE
     ! OPEN FILE
     CALL NFERR(substr                              &
          ,nf90_open(TRIM(file),NF90_NOWRITE,lncid) &
          ,1)
  END IF

  ! GET VARIABLES
  CALL NFERR(substr                            &
       ,nf90_Inquire(lncid, nVariables = nvar) &
       ,2)

  ! ALLOCATE SPACE
  ALLOCATE(var(nvar), STAT=status)
  CALL ERRMSG(substr,status,1)

  ! IMPORT VARIABLES
  DO i=1, nvar
     CALL IMPORT_NCVAR(var(i), 1, varid=i, ncid=lncid)
     ! SAVE SPACE
     CALL INIT_NARRAY(var(i)%dat)
  END DO

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL NFERR(substr       &
          ,nf90_close(lncid) &
          ,3)
  END IF

END SUBROUTINE SCAN_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RENAME_NCVAR(var, newname)

  ! rename variable of type t_ncvar

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar),   INTENT(INOUT) :: var
  CHARACTER(LEN=*), INTENT(IN)    :: newname

  IF (TRIM(var%name) == TRIM(newname)) RETURN
  var%name = TRIM(newname)

END SUBROUTINE RENAME_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IDX2FRAC_NCVAR(vi, vf, vtype, ixf)

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar), INTENT(IN)           :: vi     ! variable with index field
  TYPE (t_ncvar), INTENT(INOUT)        :: vf     ! variable with index fraction
  INTEGER,        INTENT(IN), OPTIONAL :: vtype
  INTEGER,        INTENT(IN), OPTIONAL :: ixf(2) ! range of index classes

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER        :: substr = 'IDX2FRAC_NCVAR'
  INTEGER                            :: i
  INTEGER                            :: qtype
  INTEGER                            :: lvtype
  INTEGER                            :: status
  INTEGER (I8)                       :: xmin, xmax  ! min. and max. index
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec
  INTEGER, DIMENSION(:), ALLOCATABLE :: svec
  INTEGER                            :: len         ! length of data field
  INTEGER, DIMENSION(:), POINTER     :: vivec => NULL() ! position vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: vfvec       ! position vector
  INTEGER                            :: dnew        ! new dim
  LOGICAL                            :: missval_set ! is missing value set in nc file
  INTEGER                            :: missval     ! value of missing data
  INTEGER                            :: tval        ! temoprary value

  ! INIT
  missval_set = .FALSE.
  missval = -999
  NULLIFY(vivec)

  ! INTERFACE
  IF (PRESENT(vtype)) THEN
     IF ((vtype /= VTYPE_REAL).AND.(vtype /= VTYPE_DOUBLE)) THEN
        CALL RGMSG(substr, RGMLE, &
             'INDEX FRACTION MUST BE OF TYPE REAL OR DOUBLE PRECISION !')
     ELSE
        lvtype = vtype
     END IF
  ELSE
     lvtype = VTYPE_REAL  ! DEFAULT
  END IF

  ! GET INDEX RANGE
  ! priority is as follows:
  ! First determine values from the data, which is the least reliable option
  ! as NO syncing between PEs is implemented
  ! This can be overwritten by the variable attribute RG_INDEX_RANGE
  ! This is overwritten via regrid namelists (by passing the optional ixf
  ! parameter to this subroutine)
  qtype = vi%dat%type
!CDIR BEGIN NOVECTOR
  SELECT CASE(qtype)
  CASE(VTYPE_INT)
     xmin = MINVAL(vi%dat%vi)
     xmax = MAXVAL(vi%dat%vi)
  CASE(VTYPE_REAL)
     xmin = INT(MINVAL(vi%dat%vr))
     xmax = INT(MAXVAL(vi%dat%vr))
  CASE(VTYPE_DOUBLE)
     xmin = INT(MINVAL(vi%dat%vd))
     xmax = INT(MAXVAL(vi%dat%vd))
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, 'TYPE OF VARIABLE IS CHAR !')
  CASE(VTYPE_BYTE)
     xmin = INT(MINVAL(vi%dat%vb))
     xmax = INT(MAXVAL(vi%dat%vb))
  CASE(VTYPE_UNDEF)
     CALL RGMSG(substr, RGMLE, 'TYPE OF VARIABLE IS UNDEFINED !')
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF VARIABLE !')
  END SELECT
!CDIR END

  ! GET INDEX RANGE AND MISSING VALUE FROM ATTRIBUTE
  DO i=1, vi%natts
     IF (TRIM(vi%att(i)%name) == 'RG_INDEX_RANGE') THEN
        xmin = vi%att(i)%dat%vi(1)
        xmax = vi%att(i)%dat%vi(2)
     END IF
     IF ((TRIM(vi%att(i)%name) == 'missing_value') .OR. (TRIM(vi%att(i)%name) == '_FillValue')) THEN
        missval_set = .TRUE.
        IF (vi%att(i)%dat%type == VTYPE_REAL) &
           & missval = INT(vi%att(i)%dat%vr(1))
        IF (vi%att(i)%dat%type == VTYPE_DOUBLE) &
           & missval = INT(vi%att(i)%dat%vd(1))
        IF (vi%att(i)%dat%type == VTYPE_INT) &
           & missval = INT(vi%att(i)%dat%vi(1))
        IF (vi%att(i)%dat%type == VTYPE_BYTE) &
           & missval = INT(vi%att(i)%dat%vb(1))
     END IF
  END DO

  ! min and max from optional parameter
  IF (PRESENT(ixf)) THEN
     xmin = ixf(1)
     xmax = ixf(2)
  END IF

  ! INIT
  CALL INIT_NCVAR(vf)
  vf%name  = TRIM(vi%name)
  vf%uid   = vi%uid
  vf%ustep = vi%ustep
  IF (lvtype == VTYPE_REAL)   vf%xtype = NF90_FLOAT
  IF (lvtype == VTYPE_DOUBLE) vf%xtype = NF90_DOUBLE

  ! DIMENSIONS
  vf%ndims = vi%ndims + 1
  ALLOCATE(vf%dim(vf%ndims), STAT=status)
  CALL ERRMSG(substr,status,1)
  ALLOCATE(dimvec(vf%ndims), svec(vf%ndims), STAT=status)
  CALL ERRMSG(substr,status,2)
  !
  ALLOCATE(vfvec(vf%ndims), STAT=status)
  CALL ERRMSG(substr,status,3)
  vfvec(:) = 0

  ! DISTRIBUTE DIMENSIONS
  dnew = vf%ndims   ! NEW DIMENSION AT END OF LIST
  DO i=1, vi%ndims
     IF (vi%dim(i)%fuid) THEN   ! UNLIMITED ID
        ! ULIM-DIM MUST BE LAST IN LIST ...
        vf%dim(dnew) = vi%dim(i)
        dimvec(dnew) = vi%dim(i)%len
        ! ... AND NEW DIMENSION TAKES POS. OF UNLIM
        dnew = i
     ELSE
        vf%dim(i) = vi%dim(i)
        dimvec(i) = vf%dim(i)%len
     END IF
  END DO
  ! NEW DIMENSION
  CALL INIT_NCDIM(vf%dim(dnew))
  vf%dim(dnew)%len = INT(xmax-xmin+1)
  dimvec(dnew) = vf%dim(dnew)%len
  vf%dim(dnew)%name = TRIM(vi%name)//'_idx'

  ! ATTRIBUTES
  vf%natts = vi%natts
  ALLOCATE(vf%att(vf%natts), STAT=status)
  CALL ERRMSG(substr,status,4)
  DO i=1, vi%natts
     CALL COPY_NCATT(vf%att(i), vi%att(i))
  END DO
  CALL ADD_NCATT(vf, 'RG_INDEX_RANGE', vi = (/xmin, xmax/) )

  ! DATA
  CALL INIT_NARRAY(vf%dat, vf%ndims, dimvec, lvtype)
  len = PRODUCT(vi%dat%dim)
  !
  DO i=1, len
     svec(:) = dimvec(:)
     svec(dnew) = svec(vf%ndims)
     CALL ELEMENT(svec(1:vi%ndims),i,vivec)
     vfvec(1:vi%ndims) = vivec(:)
     vfvec(vf%ndims) = vfvec(dnew)
     tval = GET_NARRAY_ELEMENT_I(vi%dat, vivec)
     vfvec(dnew) = INT(tval - xmin +1)
     IF (.NOT. (missval_set .AND. (tval == missval))) THEN
        IF ((tval < xmin) .OR. (tval > xmax)) &
           & CALL RGMSG(substr, RGMLE,  'INPUT DATA OUT OF SPECIFIED RANGE FOR IXF!')
        SELECT CASE(lvtype)
        CASE(VTYPE_REAL)
           vf%dat%vr(POSITION(dimvec,vfvec)) = 1.0
        CASE(VTYPE_DOUBLE)
           vf%dat%vd(POSITION(dimvec,vfvec)) = REAL(1.0, DP)
        CASE(VTYPE_INT)
           vf%dat%vi(POSITION(dimvec,vfvec)) = 1
        CASE(VTYPE_CHAR)
           vf%dat%vc(POSITION(dimvec,vfvec)) = 'X'
        CASE(VTYPE_BYTE)
           vf%dat%vb(POSITION(dimvec,vfvec)) = 1
        CASE(VTYPE_UNDEF)
           CALL RGMSG(substr, RGMLE, 'TYPE OF VARIABLE IS UNDEFINED !')
        CASE DEFAULT
           CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF VARIABLE !')
        END SELECT
     END IF
     DEALLOCATE(vivec, STAT=status)
     CALL ERRMSG(substr,status,5)
     NULLIFY(vivec)
  END DO

  ! CLEAN UP
  DEALLOCATE(dimvec, svec, vfvec, STAT=status)
  CALL ERRMSG(substr,status,6)

CONTAINS

  INTEGER FUNCTION GET_NARRAY_ELEMENT_I(na, vec)

    IMPLICIT NONE

    ! I/O
    TYPE (t_narray),             INTENT(IN) :: na
    INTEGER,       DIMENSION(:), INTENT(IN) :: vec

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'GET_NARRAY_ELEMENT_I'
    INTEGER :: i
    INTEGER :: vtype
    INTEGER :: pos

    IF (SIZE(vec) /= na%n) THEN
       CALL RGMSG(substr, RGMLE, 'NUMBER OF DIMENSIONS MISMATCH !')
    END IF

    DO i=1, na%n
       IF ((vec(i) < 1).OR.(vec(i) > na%dim(i))) THEN
       CALL RGMSG(substr, RGMLE, &
            'POSITION VECTOR ELEMENT ',i,' OUT OF RANGE:')
             CALL RGMSG(substr, RGMLEC, 'DIM: ',na%dim,' ', .false.)
             CALL RGMSG(substr, RGMLEC, 'VEC: ',vec,' ')
       END IF
    END DO

    pos   = POSITION(na%dim, vec)
    vtype = na%type

    SELECT CASE(vtype)
    CASE(VTYPE_INT)
       GET_NARRAY_ELEMENT_I = INT(na%vi(pos))
    CASE(VTYPE_REAL)
       GET_NARRAY_ELEMENT_I = INT(na%vr(pos))
    CASE(VTYPE_DOUBLE)
       GET_NARRAY_ELEMENT_I = INT(na%vd(pos))
    CASE(VTYPE_CHAR)
       CALL RGMSG(substr, RGMLE, 'TYPE OF N-ARRAY IS CHAR !')
    CASE(VTYPE_BYTE)
       GET_NARRAY_ELEMENT_I = INT(na%vb(pos))
    CASE(VTYPE_UNDEF)
       CALL RGMSG(substr, RGMLE, 'TYPE OF N-ARRAY IS UNDEFINED !')
    CASE DEFAULT
       CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF N-ARRAY !')
    END SELECT

  END FUNCTION GET_NARRAY_ELEMENT_I

END SUBROUTINE IDX2FRAC_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE MAXFRAC2IDX_NCVAR(vf, vi, vtype)

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar), INTENT(IN)         :: vf    ! index fractions
  TYPE (t_ncvar), INTENT(INOUT)      :: vi    ! index
  INTEGER,      INTENT(IN), OPTIONAL :: vtype

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER        :: substr = 'MAXFRAC2IDX_NCVAR'
  INTEGER                            :: lvtype
  INTEGER (I8)                       :: xr(2)   ! xmin, xmax (index range)
  INTEGER                            :: i
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec
  INTEGER, DIMENSION(:), ALLOCATABLE :: svec
  INTEGER, DIMENSION(:), POINTER     :: vivec => NULL() ! position vector
  TYPE (t_narray)                    :: n1      ! LAST dim. extracted
  TYPE (t_narray)                    :: n1i     ! SORT INDEX
  INTEGER                            :: len     ! length of index data
  INTEGER                            :: xlen    ! length of last dim
  INTEGER                            :: status
  INTEGER                            :: dnew    ! pos. of new dimension

  ! INTERFACE
  IF (PRESENT(vtype)) THEN
     lvtype = vtype
  ELSE
     lvtype = VTYPE_REAL  ! DEFAULT
  END IF

  ! GET INDEX RANGE FROM ATTRIBUTE
  xr(:) = 0
  DO i=1, vf%natts
     IF (TRIM(vf%att(i)%name) == 'RG_INDEX_RANGE') THEN
        xr(:) = vf%att(i)%dat%vi(1:2)
     END IF
  END DO

  ! INIT
  CALL INIT_NCVAR(vi)
  vi%name  = TRIM(vf%name)
  vi%uid   = vf%uid
  vi%ustep = vf%ustep
  SELECT CASE(lvtype)
  CASE(VTYPE_REAL)
     vi%xtype = NF90_FLOAT
  CASE(VTYPE_DOUBLE)
     vi%xtype = NF90_DOUBLE
  CASE(VTYPE_INT)
     vi%xtype = NF90_INT
  CASE(VTYPE_BYTE)
     vi%xtype = NF90_BYTE
  CASE(VTYPE_CHAR)
     vi%xtype = NF90_CHAR
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, 'INDEX TYPE IS UNDEFINED OR UNRECOGNIZED !')
  END SELECT

  ! DIMENSIONS
  dnew = vf%ndims
  vi%ndims = vf%ndims - 1
  IF (vi%ndims > 0) THEN
     ALLOCATE(vi%dim(vi%ndims), STAT=status)
     CALL ERRMSG(substr,status,1)
     ALLOCATE(dimvec(vi%ndims), svec(vf%ndims), STAT=status)
     CALL ERRMSG(substr,status,2)
     DO i=1, vi%ndims
        vi%dim(i) = vf%dim(i)
        dimvec(i) = vf%dim(i)%len
     END DO
     IF (vf%dim(vf%ndims)%fuid) THEN   ! UNLIMITED ID IS LAST IN LIST
        dnew = vf%ndims - 1            ! THE NEW IS THEN LAST BUT ONE
        vi%dim(vi%ndims) = vf%dim(vf%ndims)
        dimvec(vi%ndims) = vf%dim(vf%ndims)%len
     END IF
  END IF

  ! ATTRIBUTES
  vi%natts = vf%natts
  ALLOCATE(vi%att(vi%natts), STAT=status)
  CALL ERRMSG(substr,status,3)
  DO i=1, vi%natts
     CALL COPY_NCATT(vi%att(i), vf%att(i))
  END DO

  ! DATA
  CALL INIT_NARRAY(vi%dat, vi%ndims, dimvec, lvtype)
  len = PRODUCT(dimvec)
  xlen = vf%dim(dnew)%len
  !
  DO i=1, len
     CALL ELEMENT(dimvec,i,vivec)
     svec(:) = 0
     svec(1:vi%ndims) = vivec(:)
     svec(vf%ndims) = svec(dnew)
     svec(dnew) = 0
     CALL EXTRACT_NARRAY_1DIM(vf%dat,svec,n1)
     CALL SORT_NARRAY(n1,n1i)
     SELECT CASE(lvtype)
     CASE(VTYPE_INT)
        vi%dat%vi(i) = n1i%vi(xlen) - 1 + xr(1)
     CASE(VTYPE_REAL)
        vi%dat%vr(i) = REAL(n1i%vi(xlen) - 1 + xr(1), SP)
     CASE(VTYPE_DOUBLE)
        vi%dat%vd(i) = REAL(n1i%vi(xlen) - 1 + xr(1), DP)
     CASE(VTYPE_CHAR)
        CALL RGMSG(substr, RGMLE, 'N-ARRAY IS OF TYPE CHAR !')
     CASE(VTYPE_BYTE)
        vi%dat%vb(i) = INT(n1i%vi(xlen), I4) - INT(1, I4) + INT(xr(1), I4)
     CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, 'TYPE OF N-ARRAY IS UNDEFINED !')
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF N-ARRAY!')
     END SELECT
     DEALLOCATE(vivec, STAT=status)
     CALL ERRMSG(substr,status,4)
     NULLIFY(vivec)
     CALL INIT_NARRAY(n1)
     CALL INIT_NARRAY(n1i)
  END DO

  ! CLEAN UP
  DEALLOCATE(dimvec, svec, STAT=status)
  CALL ERRMSG(substr,status,5)

CONTAINS

  SUBROUTINE EXTRACT_NARRAY_1DIM(na, vec, n1)

    IMPLICIT NONE

    ! I/O
    TYPE (t_narray),            INTENT(IN)  :: na
    INTEGER,      DIMENSION(:), INTENT(IN)  :: vec ! vec. with 0 at ...
                                                   ! ... dim. to extract
    TYPE (t_narray),            INTENT(INOUT) :: n1

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'EXTRACT_NARRAY_1DIM'
    INTEGER                     :: vtype
    INTEGER                     :: i
    INTEGER                     :: nvec(na%n)
    INTEGER                     :: pos
    INTEGER                     :: zcount  ! how many zero's
    INTEGER                     :: zpos    ! where is the zero

    IF (SIZE(vec) /= na%n ) THEN
       CALL RGMSG(substr, RGMLE, 'NUMBER OF DIMENSIONS MISMATCH !')
    END IF

    zcount = 0
    zpos = 0
    DO i=1, na%n
       IF (vec(i) /= 0) THEN  ! THIS IS DIM TO EXTRACT
          IF ((vec(i) < 1).OR.(vec(i) > na%dim(i))) THEN
             CALL RGMSG(substr, RGMLE, &
                  'POSITION VECTOR ELEMENT ',i,' OUT OF RANGE:', .false.)
             CALL RGMSG(substr, RGMLEC, 'DIM: ',na%dim,' ', .false.)
             CALL RGMSG(substr, RGMLEC, 'VEC: ',vec,' ')
          END IF
       ELSE
          zcount = zcount + 1
          zpos = i
       END IF
    END DO

    IF (zcount /= 1) THEN
       CALL RGMSG(substr, RGMLE, 'ONE ENTRY IN VECTOR MUST BE ZERO !')
    END IF

    vtype = na%type
    CALL INIT_NARRAY(n1, 1, (/ na%dim(zpos) /), vtype)

    nvec(:) = vec(:)
    DO i=1, na%dim(zpos)
       nvec(zpos) = i
       pos = POSITION(na%dim, nvec)
       !
       SELECT CASE(vtype)
       CASE(VTYPE_INT)
          n1%vi(i) = na%vi(pos)
       CASE(VTYPE_REAL)
          n1%vr(i) = na%vr(pos)
       CASE(VTYPE_DOUBLE)
          n1%vd(i) = na%vd(pos)
       CASE(VTYPE_CHAR)
          n1%vc(i) = na%vc(pos)
       CASE(VTYPE_BYTE)
          n1%vb(i) = na%vb(pos)
       CASE(VTYPE_UNDEF)
          CALL RGMSG(substr, RGMLE, 'TYPE OF N-ARRAY IS UNDEFINED !')
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF N-ARRAY !')
       END SELECT
    END DO

  END SUBROUTINE EXTRACT_NARRAY_1DIM

END SUBROUTINE MAXFRAC2IDX_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PRINT_NCVAR(var, str, long, verbose_level)

  ! write out variable of type t_ncvar
  ! str can be used to determine location of print statement

  IMPLICIT NONE

   TYPE (t_ncvar)               :: var ! variable
   CHARACTER(LEN=*), INTENT(IN) :: str
   LOGICAL,          INTENT(IN), OPTIONAL :: long
   INTEGER,          INTENT(IN), OPTIONAL :: verbose_level

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER  :: substr = 'PRINT_NCVAR'
   INTEGER                      :: i
   INTEGER                      :: iverb

  ! verbose_level
  ! 100: print all
  !  80: all but no attributes
  iverb = 100
  IF (PRESENT(verbose_level) ) iverb = verbose_level

   write (0,*) substr, ' ', TRIM(str),'+++++++++++++++++++++++++++++++++++++'
   write (0,*) substr, ' ', TRIM(str),'+++++++++++++++++++++++++++++++++++++'
   write (0,*) substr, ' ', TRIM(str), ' VAR NAME: ',      TRIM(var%name)
   write (0,*) substr, ' ', TRIM(str),'+++++++++++++++++++++++++++++++++++++'
   write (0,*) substr, ' ', TRIM(str), ' VAR ID / TYPE: ',  var%id, ' ',var%xtype
   write (0,*) substr, ' ', TRIM(str), ' VAR DIMENSIONS: ', var%ndims
   DO i=1,var%ndims
      write (0,*) substr, ' ', TRIM(str), ' VAR DIM: ',i, TRIM(var%dim(i)%name)&
           , var%dim(i)%id , var%dim(i)%len, var%dim(i)%fuid, var%dim(i)%varid
   ENDDO
   write (0,*) substr, ' ', TRIM(str), ' VAR uid: ',        var%uid
   IF (iverb > 80) THEN
   write (0,*) substr, ' ', TRIM(str),'+++++++++++++++++++++++++++++'
   write (0,*) substr, ' ', TRIM(str), ' VAR ATTRIBUTES: ', var%natts
   DO i=1,var%natts
      write (0,*) substr,' ',TRIM(str), ' VAR ATTS: ',i, TRIM(var%att(i)%name) &
           , var%att(i)%ID, var%att(i)%xtype, var%att(i)%varid &
           , var%att(i)%len
      CALL PRINT_NARRAY(var%att(i)%dat,  TRIM(str)//TRIM(var%name))
      write (0,*) substr, ' ', TRIM(str),'+++++++++++++++++++++++++++++'
   ENDDO
   ENDIF

   write (0,*) substr, ' ', TRIM(str),'+++++++++++++++++++++++++++++++++++++'
   write (0,*) substr,  ' ',TRIM(str), ' VAR DATA: '
   CALL PRINT_NARRAY(var%dat,   TRIM(str)//TRIM(var%name), long)
   write (0,*) substr, ' ', TRIM(str),'+++++++++++++++++++++++++++++++++++++'
   write (0,*) substr, ' ', TRIM(str),'+++++++++++++++++++++++++++++++++++++'


END SUBROUTINE PRINT_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PRINT_NARRAY(array, str, long)

  ! write out array of type t_narray
  ! str can be used to determine location of print statement

  IMPLICIT NONE

   TYPE (t_narray),  INTENT(IN)  :: array ! array
   CHARACTER(LEN=*), INTENT(IN)  :: str
   LOGICAL,          INTENT(IN), OPTIONAL :: long

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER  :: substr = 'PRINT_NARRAY'
   INTEGER                      :: i
   LOGICAL                      :: llong

   IF (PRESENT(long)) THEN
      llong = long
   ELSE
      llong= .FALSE.
   ENDIF

   write (0,*) substr, ' ',TRIM(str),' # DIMENSION', array%type,' ', array%n
   DO i = 1, array%n
      write (0,*) substr,' ', i, 'DIMENSION ',TRIM(str), array%dim(i)
   ENDDO
   SELECT CASE (array%type)
      CASE(VTYPE_UNDEF)
         ! nothing to print
      CASE(VTYPE_INT)
         write (0,*) substr,' ', TRIM(str),i, ' MINVAL/MAXVAL INT'&
              , MINVAL(array%vi), MAXVAL(array%vi)
      CASE(VTYPE_REAL)
         IF (llong) THEN
            write (0,*) substr,' ', TRIM(str),i, SIZE(array%vr), 'DATA ',array%vr
         ELSE
            write (0,*) substr,' ', TRIM(str),i, ' MINVAL/MAXVAL REAL'&
                 , SIZE(array%vr), MINVAL(array%vr), MAXVAL(array%vr)
         ENDIF
      CASE(VTYPE_DOUBLE)
         IF (llong) THEN
            write (0,*) substr,' ', TRIM(str),i, SIZE(array%vd), 'DATA ',array%vd
         ELSE
            write (0,*) substr,' ', TRIM(str),i, ' MINVAL/MAXVAL DOUBLE'&
                 , SIZE(array%vd), MINVAL(array%vd), MAXVAL(array%vd)
         END IF
      CASE(VTYPE_BYTE)
         write (0,*) substr,' ', TRIM(str),i, ' ',array%vb
      CASE(VTYPE_CHAR)
         write (0,*) substr,' ', TRIM(str),i, ' ',array%vc

      END SELECT

 END SUBROUTINE PRINT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PRINT_NCDIM(dim, str)

  ! write out dimension of type t_ncdim

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncdim),   INTENT(IN) :: dim
  CHARACTER(LEN=*), INTENT(IN) :: str

  write (0,*)  str,' NCDIM NAME',  TRIM(dim%name)
  write (0,*)  str,' NCDIM ID',    dim%id
  write (0,*)  str,' NCDIM LEN',   dim%len
  write (0,*)  str,' NCDIM FIUD',  dim%fuid
  write (0,*)  str,' NCDIM VARID', dim%varid

END SUBROUTINE PRINT_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PRINT_NCATT(att, str)

  ! write out dimension of type t_ncdim

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncatt),   INTENT(IN) :: att
  CHARACTER(LEN=*), INTENT(IN) :: str

  write (0,*) '*************************************************'
  write (0,*) '***** PRINT NCATT  ******************************'
  write (0,*) '*************************************************'
  write (0,*) 'NCATT NAME',  str, TRIM(att%name)
  write (0,*) 'NCATT ID',    str, att%ID
  write (0,*) 'NCATT XTYPE', str, att%xtype
  write (0,*) 'NCATT LEN',   str, att%len
  write (0,*) 'NCATT VARID', str, att%varid

  IF (att%xtype /= NULL_XTYPE) THEN
     CALL PRINT_NARRAY(att%dat, 'ATT '//str)
  END IF
  write (0,*) '*************************************************'

END SUBROUTINE PRINT_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXTRACT_NCATT(att, type, name, i, c, r)

  IMPLICIT NONE
  INTRINSIC :: INT, MIN, REAL, LEN

  ! I/O
  TYPE (t_ncatt),   INTENT(IN)    :: att
  INTEGER,          INTENT(OUT)   :: type ! 1: integer, 2: string, 3: real(dp)
  CHARACTER(LEN=*), INTENT(OUT)   :: name
  INTEGER,          INTENT(INOUT) :: i
  CHARACTER(LEN=*), INTENT(INOUT) :: c
  REAL(DP),         INTENT(INOUT) :: r

  ! LOCAL
  INTEGER :: n

  name = att%name

  SELECT CASE(att%xtype)
  CASE(NF90_BYTE)
     type = 1
     i = INT(att%dat%vb(1))
  CASE(NF90_SHORT)
     type = 1
     i = INT(att%dat%vi(1))
  CASE(NF90_CHAR)
     type = 2
     c = ''
     DO n=1, MIN(att%len, LEN(c))
        c(n:n) = att%dat%vc(n)
     END DO
  CASE(NF90_INT)
     type = 1
     i = INT(att%dat%vi(1))
  CASE(NF90_FLOAT)
     type = 3
     r = REAL(att%dat%vr(1), dp)
  CASE(NF90_DOUBLE)
     type = 3
     r =att%dat%vd(1)
  CASE DEFAULT
     type = 0
     !
  END SELECT

END SUBROUTINE EXTRACT_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE NFERR(routine, status, no)

  ! improved netcdf error output routines

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: status
  INTEGER,          INTENT(IN), OPTIONAL  :: no

  IF (status /= NF90_NOERR) THEN
     IF (PRESENT(no)) THEN
        CALL RGMSG(routine, RGMLE, &
             'netCDF ERROR AT POSITION ',no,' !', .false.)
     ELSE
        CALL RGMSG(routine, RGMLE, &
             'netCDF ERROR !', .false.)
     END IF
     CALL RGMSG(routine, RGMLEC, &
          'netCDF ERROR-STATUS: ',status,' ', .false.)
     CALL RGMSG(routine, RGMLEC, nf90_strerror(status))
  END IF

END SUBROUTINE NFERR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
CHARACTER(LEN=100*GRD_MAXSTRLEN) FUNCTION string(c)

  IMPLICIT NONE

  ! I/O
  CHARACTER, DIMENSION(:), POINTER :: c

  ! LOCAL
  INTEGER :: i

  string = ''

  IF (.NOT. ASSOCIATED(c)) RETURN

  DO i=1, MIN(SIZE(c), 100*GRD_MAXSTRLEN)
     string(i:i) = c(i)
  END DO

END FUNCTION string
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_NARRAY(na, n, dim, qtype)

  ! initialise array of type narray
  ! INPUT:
  !   - n :   number of dimensions
  !   - dim:  length of dimensions
  !   - type: type of array (REAL, REAL(dp), INTEGER, CHAR, BYTE)

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray),       INTENT(INOUT)          :: na
  INTEGER,               INTENT(IN),   OPTIONAL :: n
  INTEGER, DIMENSION(:), INTENT(IN),   OPTIONAL :: dim
  INTEGER,                             OPTIONAL :: qtype

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'INIT_NARRAY'
  INTEGER                     :: status
  INTEGER                     :: len

  na%type = VTYPE_UNDEF

  na%n = 0
  IF (ASSOCIATED(na%dim)) THEN
     IF (SIZE(na%dim) > 0) THEN
        DEALLOCATE(na%dim, STAT=status)
        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(na%dim)
  END IF

  IF (ASSOCIATED(na%vi)) THEN
     IF (SIZE(na%vi) > 0) THEN
        DEALLOCATE(na%vi, STAT=status)
        CALL ERRMSG(substr,status,2)
     END IF
     NULLIFY(na%vi)
  END IF

  IF (ASSOCIATED(na%vr)) THEN
     IF (SIZE(na%vr) > 0) THEN
        DEALLOCATE(na%vr, STAT=status)
        CALL ERRMSG(substr,status,3)
     END IF
     NULLIFY(na%vr)
  END IF

  IF (ASSOCIATED(na%vd)) THEN
     IF (SIZE(na%vd) > 0) THEN
        DEALLOCATE(na%vd, STAT=status)
        CALL ERRMSG(substr,status,4)
     END IF
     NULLIFY(na%vd)
  END IF

  IF (ASSOCIATED(na%vc)) THEN
     IF (SIZE(na%vc) > 0) THEN
        DEALLOCATE(na%vc, STAT=status)
        CALL ERRMSG(substr,status,5)
     END IF
     NULLIFY(na%vc)
  END IF

  IF (ASSOCIATED(na%vb)) THEN
     IF (SIZE(na%vb) > 0) THEN
        DEALLOCATE(na%vb, STAT=status)
        CALL ERRMSG(substr,status,6)
     END IF
     NULLIFY(na%vb)
  END IF

  IF (PRESENT(dim) .AND. PRESENT(n)) THEN
     IF (SIZE(dim) /= n) THEN
        CALL RGMSG(substr, RGMLE,                &
             'NUMBER OF DIMENSIONS MISMATCH !' )
     END IF
  END IF

  len = 0
  IF (PRESENT(n)) THEN
     na%n = n
     ALLOCATE(na%dim(n), STAT=status)
     CALL ERRMSG(substr,status,7)
     na%dim(:) = 0
  END IF

  IF (PRESENT(dim)) THEN
     IF (.NOT. ASSOCIATED(na%dim)) THEN
        na%n = SIZE(dim)
        ALLOCATE(na%dim(na%n), STAT=status)
        CALL ERRMSG(substr,status,8)
     END IF
     na%dim(:) = dim(:)
     len = PRODUCT(na%dim(:))
  END IF

  IF (PRESENT(qtype).AND.(len > 0)) THEN
     SELECT CASE(qtype)
     CASE(VTYPE_INT)
        na%type = VTYPE_INT
        ALLOCATE(na%vi(len), STAT=status)
        CALL ERRMSG(substr,status,9)
        na%vi(:) = 0
     CASE(VTYPE_REAL)
        na%type = VTYPE_REAL
        ALLOCATE(na%vr(len), STAT=status)
        CALL ERRMSG(substr,status,10)
        na%vr(:) = 0.0
     CASE(VTYPE_DOUBLE)
        na%type = VTYPE_DOUBLE
        ALLOCATE(na%vd(len), STAT=status)
        CALL ERRMSG(substr,status,11)
        na%vd(:) = REAL(0.0, DP)
     CASE(VTYPE_CHAR)
        na%type = VTYPE_CHAR
        ALLOCATE(na%vc(len), STAT=status)
        CALL ERRMSG(substr,status,12)
        na%vc(:) = ' '
     CASE(VTYPE_BYTE)
        na%type = VTYPE_BYTE
        ALLOCATE(na%vb(len), STAT=status)
        CALL ERRMSG(substr,status,13)
        na%vb(:) = 0
     CASE(VTYPE_UNDEF)
        na%type = VTYPE_UNDEF
        ! DO NOTHING, KEEP UNDEFINED
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE,                           &
             'REQUESTED TYPE FOR N-ARRAY IS UNRECOGNIZED !' )
     END SELECT
  END IF

END SUBROUTINE INIT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_NARRAY(d, s)

  ! copy array of type t_narray from source (s) to destination (d)

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(INOUT) :: d
  TYPE (t_narray), INTENT(IN)    :: s

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER  :: substr = 'COPY_NARRAY'
  INTEGER                      :: n
  INTEGER                      :: vtype
  INTEGER                      :: status

  ! INIT
  n = 0

  CALL INIT_NARRAY(d)

  d%n = s%n
  IF (ASSOCIATED(s%dim)) THEN
     ALLOCATE(d%dim(d%n),STAT=status)
     CALL ERRMSG(substr,status,1)
     d%dim(:) = s%dim(:)
  END IF

  vtype = s%type
  SELECT CASE(vtype)
  CASE(VTYPE_INT)
     d%type = VTYPE_INT
     n = SIZE(s%vi)
     IF (n > 0) THEN
        ALLOCATE(d%vi(n),STAT=status)
        CALL ERRMSG(substr,status,2)
        d%vi(:) = s%vi(:)
     END IF
  CASE(VTYPE_REAL)
     d%type = VTYPE_REAL
     n = SIZE(s%vr)
     IF (n > 0) THEN
        ALLOCATE(d%vr(n),STAT=status)
        CALL ERRMSG(substr,status,3)
        d%vr(:) = s%vr(:)
     END IF
  CASE(VTYPE_DOUBLE)
     d%type = VTYPE_DOUBLE
     n = SIZE(s%vd)
     IF (n > 0) THEN
        ALLOCATE(d%vd(n),STAT=status)
        CALL ERRMSG(substr,status,4)
        d%vd(:) = s%vd(:)
     END IF
  CASE(VTYPE_CHAR)
     d%type = VTYPE_CHAR
     n = SIZE(s%vc)
     IF (n > 0) THEN
        ALLOCATE(d%vc(n),STAT=status)
        CALL ERRMSG(substr,status,5)
        d%vc(:) = s%vc(:)
     END IF
  CASE(VTYPE_BYTE)
     d%type = VTYPE_BYTE
     n = SIZE(s%vb)
     IF (n > 0) THEN
        ALLOCATE(d%vb(n),STAT=status)
        CALL ERRMSG(substr,status,6)
        d%vb(:) = s%vb(:)
     END IF
  CASE(VTYPE_UNDEF)
     d%type = VTYPE_UNDEF
     ! DO NOTHING, EMPTY N-ARRAY IS COPIED
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE,                              &
          'N-ARRAY OF UNRECOGNIZED TYPE CANNOT BE COPIED !' )
  END SELECT

END SUBROUTINE COPY_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SORT_NARRAY(na, nx, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(INOUT)          :: na
  TYPE (t_narray), INTENT(INOUT)          :: nx
  LOGICAL        , INTENT(IN)   ,OPTIONAL :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SORT_NARRAY'
  INTEGER         :: vtype
  LOGICAL         :: lrev
  INTEGER         :: n, i
  TYPE (t_narray) :: nh

  IF (PRESENT(reverse)) THEN
     lrev = reverse
  ELSE
     lrev = .false.  ! DEFAULT
  END IF

  IF (na%n == 0) THEN
     CALL RGMSG(substr, RGMLW, 'EMPTY ARRAY ! NOTHING TO DO !')
     RETURN
  END IF

  IF (na%n > 1) THEN
     CALL RGMSG(substr, RGMLW, &
          'SORTING A ',na%n,'-DIMENSIONAL N-ARRAY AS LINEAR ARRAY !')
  END IF

  IF (lrev) THEN

     vtype = nx%type
     IF (vtype /= VTYPE_INT) THEN
        CALL RGMSG(substr, RGMLE, 'INDEX N-ARRAY MUST BE OF TYPE INTEGER !')
     END IF
     n = SIZE(nx%vi)
     CALL COPY_NARRAY(nh, na)   ! COPY TO BE RE-ORDERED
     vtype = na%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        DO i=1,n
           na%vr(nx%vi(i)) = nh%vr(i)
        END DO
     CASE(VTYPE_DOUBLE)
        DO i=1,n
           na%vd(nx%vi(i)) = nh%vd(i)
        END DO
     CASE(VTYPE_INT)
        DO i=1,n
           na%vi(nx%vi(i)) = nh%vi(i)
        END DO
     CASE(VTYPE_BYTE)
        DO i=1,n
           na%vb(nx%vi(i)) = nh%vb(i)
        END DO
     CASE(VTYPE_CHAR)
     CALL RGMSG(substr,RGMLE,'UN-SORTING OF TYPE CHAR IS NOT IMPLEMENTED !')
     CASE(VTYPE_UNDEF)
     CALL RGMSG(substr,RGMLE,'ARRAY OF UNDEFINED TYPE CANNOT BE UN-SORTED !')
     CASE DEFAULT
   CALL RGMSG(substr,RGMLE,'ARRAY OF UNRECOGNIZED TYPE CANNOT BE UN-SORTED !')
     END SELECT
     ! CLEAN UP
     CALL INIT_NARRAY(nh)

  ELSE

     CALL INIT_NARRAY(nx, na%n, na%dim, VTYPE_INT)
     DEALLOCATE(nx%vi)
     NULLIFY(nx%vi)
     !
     vtype = na%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        CALL QSORT_R(na%vr, nx%vi)
     CASE(VTYPE_DOUBLE)
        CALL QSORT_D(na%vd, nx%vi)
     CASE(VTYPE_INT)
        CALL QSORT_I(na%vi, nx%vi)
     CASE(VTYPE_BYTE)
        CALL QSORT_B(na%vb, nx%vi)
     CASE(VTYPE_CHAR)
     CALL RGMSG(substr,RGMLE,'SORTING OF TYPE CHAR IS NOT IMPLEMENTED !')
     CASE(VTYPE_UNDEF)
     CALL RGMSG(substr,RGMLE,'ARRAY OF UNDEFINED TYPE CANNOT BE SORTED !')
     CASE DEFAULT
     CALL RGMSG(substr,RGMLE,'ARRAY OF UNRECOGNIZED TYPE CANNOT BE SORTED !')
     END SELECT

  END IF ! (reverse ?)

END SUBROUTINE SORT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE REORDER_NARRAY(na, nx)

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(INOUT) :: na   ! n-array to reorder
  TYPE (t_narray), INTENT(IN)    :: nx   ! index n-array

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'REORDER_NARRAY'
  TYPE (t_narray)             :: nh   ! copy of na
  INTEGER                     :: vtype
  INTEGER                     :: i, n

  IF (na%n == 0) THEN
     CALL RGMSG(substr, RGMLW, 'EMPTY ARRAY ! NOTHING TO DO !')
     RETURN
  END IF

  vtype = nx%type
  IF (vtype /= VTYPE_INT) THEN
     CALL RGMSG(substr, RGMLE, 'INDEX N-ARRAY MUST BE OF TYPE INT !')
  END IF

  IF (na%n > 1) THEN
     CALL RGMSG(substr, RGMLW, 'REORDERING A ',na%n, &
          '-DIMENSIONAL N-ARRAY AS LINEAR ARRAY !')
  END IF

  IF (na%n /= nx%n) THEN
     CALL RGMSG(substr, RGMLE, 'DIMENSION MISMATCH BETWEEN N-ARRAY (',&
          na%n,')' , .false.)
     CALL RGMSG(substr, RGMLEC, 'AND INDEX N-ARRAY (',nx%n,') !')
  END IF

  IF (PRODUCT(na%dim) /= PRODUCT(nx%dim)) THEN
     CALL RGMSG(substr, RGMLE, 'LENGTH OF N-ARRAY MISMATCH BETWEEN', .false.)
     CALL RGMSG(substr, RGMLEC, 'N-ARRAY (',PRODUCT(na%dim),') AND', .false.)
     CALL RGMSG(substr, RGMLEC, 'INDEX N-ARRAY (',PRODUCT(nx%dim),') !')
  END IF

  vtype = na%type
  CALL COPY_NARRAY(nh, na)
  n = PRODUCT(na%dim)

  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
     DO i=1,n
        na%vr(i) = nh%vr(nx%vi(i))
     END DO
  CASE(VTYPE_DOUBLE)
     DO i=1,n
        na%vd(i) = nh%vd(nx%vi(i))
     END DO
  CASE(VTYPE_INT)
     DO i=1,n
        na%vi(i) = nh%vi(nx%vi(i))
     END DO
  CASE(VTYPE_BYTE)
     DO i=1,n
        na%vb(i) = nh%vb(nx%vi(i))
     END DO
  CASE(VTYPE_CHAR)
     DO i=1,n
        na%vc(i) = nh%vc(nx%vi(i))
     END DO
  CASE(VTYPE_UNDEF)
  CALL RGMSG(substr, RGMLE, 'ARRAY OF UNDEFINED TYPE CANNOT BE UN-SORTED !')
  CASE DEFAULT
  CALL RGMSG(substr, RGMLE, 'ARRAY OF UNRECOGNIZED TYPE CANNOT BE UN-SORTED !')
  END SELECT

  ! CLEAN UP
  CALL INIT_NARRAY(nh)

END SUBROUTINE REORDER_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE DOUBLE_NARRAY(na)

  ! convert data of the variable of type t_narray to double precision

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(INOUT) :: na

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'DOUBLE_NARRAY'
  INTEGER                     :: vtype
  INTEGER                     :: status

  vtype = na%type
  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
     ALLOCATE(na%vd(SIZE(na%vr)), STAT=status)
     CALL ERRMSG(substr,status,1)
     na%vd(:) = REAL(na%vr(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vr, STAT=status)
     CALL ERRMSG(substr,status,2)
     NULLIFY(na%vr)
  CASE(VTYPE_DOUBLE)
     ! NOTHING TO DO
  CASE(VTYPE_INT)
     ALLOCATE(na%vd(SIZE(na%vi)), STAT=status)
     CALL ERRMSG(substr,status,3)
     na%vd(:) = REAL(na%vi(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vi, STAT=status)
     CALL ERRMSG(substr,status,4)
     NULLIFY(na%vi)
  CASE(VTYPE_BYTE)
     ALLOCATE(na%vd(SIZE(na%vi)), STAT=status)
     CALL ERRMSG(substr,status,5)
     na%vd(:) = REAL(na%vb(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vb, STAT=status)
     CALL ERRMSG(substr,status,6)
     NULLIFY(na%vb)
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, 'CHAR CANNOT BE CONVERTED TO DOUBLEPRECISION !')
  CASE(VTYPE_UNDEF)
     CALL RGMSG(substr, RGMLE, 'UNDEFINED N-ARRAY TYPE !')
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED N-ARRAY TYPE !')
  END SELECT

END SUBROUTINE DOUBLE_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE REAL_SP_NARRAY(na)

  ! convert data of the variable of type t_narray to single precision

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(INOUT) :: na

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER    :: substr = 'DOUBLE_NARRAY'
  INTEGER                        :: vtype
  INTEGER                        :: status

  vtype = na%type
  SELECT CASE(vtype)
  CASE(VTYPE_DOUBLE)
     ALLOCATE(na%vr(SIZE(na%vd)), STAT=status)
     CALL ERRMSG(substr,status,1)
     na%vr(:) = REAL(na%vd(:), SP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vr, STAT=status)
     CALL ERRMSG(substr,status,2)
     NULLIFY(na%vd)
  CASE(VTYPE_REAL)
     ! NOTHING TO DO
  CASE(VTYPE_INT)
     ALLOCATE(na%vr(SIZE(na%vi)), STAT=status)
     CALL ERRMSG(substr,status,3)
     na%vr(:) = REAL(na%vi(:), SP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vi, STAT=status)
     CALL ERRMSG(substr,status,4)
     NULLIFY(na%vi)
  CASE(VTYPE_BYTE)
     ALLOCATE(na%vr(SIZE(na%vi)), STAT=status)
     CALL ERRMSG(substr,status,5)
     na%vr(:) = REAL(na%vb(:), SP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vb, STAT=status)
     CALL ERRMSG(substr,status,6)
     NULLIFY(na%vb)
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, 'CHAR CANNOT BE CONVERTED TO DOUBLEPRECISION !')
  CASE(VTYPE_UNDEF)
     CALL RGMSG(substr, RGMLE, 'UNDEFINED N-ARRAY TYPE !')
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED N-ARRAY TYPE !')
  END SELECT

END SUBROUTINE REAL_SP_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SCALE_NARRAY(na, sc)

  ! scale data of the variable of type t_narray by sc

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(INOUT) :: na  ! N-array
  REAL(dp)       , INTENT(IN)    :: sc  ! scaling factor

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SCALE_NARRAY'
  INTEGER                     :: vtype
  INTEGER                     :: status
  INTEGER                     :: i

  vtype = na%type

  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
     na%vr(:) = REAL(REAL(na%vr(:), DP) * sc, SP)
  CASE(VTYPE_DOUBLE)
     na%vd(:) = na%vd(:) * sc
  CASE(VTYPE_INT)
     CALL RGMSG(substr, RGMLI, 'N-ARRAY OF TYPE INTEGER CONVERTED TO REAL !')
     ALLOCATE(na%vr(SIZE(na%vi)), STAT=status)
     na%type = VTYPE_REAL
     CALL ERRMSG(substr,status,1)
     DO i=1, SIZE(na%vi)
        na%vr(i) = REAL(na%vi(i), SP) * REAL(sc, SP)
     END DO
     DEALLOCATE(na%vi, STAT=status)
     CALL ERRMSG(substr,status,2)
     NULLIFY(na%vi)
  CASE(VTYPE_BYTE)
     CALL RGMSG(substr, RGMLI, 'N-ARRAY OF TYPE BYTE CONVERTED TO REAL !')
     ALLOCATE(na%vr(SIZE(na%vb)), STAT=status)
     na%type = VTYPE_REAL
     CALL ERRMSG(substr,status,3)
     DO i=1, SIZE(na%vb)
        na%vr(i) = REAL(INT(na%vb(i), I8), SP) * REAL(sc, SP)
     END DO
     DEALLOCATE(na%vb, STAT=status)
     CALL ERRMSG(substr,status,4)
     NULLIFY(na%vb)
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, 'N-ARRAY OF TYPE CHAR CANNOT BE SCALED !')
  CASE(VTYPE_UNDEF)
     CALL RGMSG(substr, RGMLE, 'UNDEFINED N-ARRAY TYPE !')
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED N-ARRAY TYPE !')
  END SELECT

END SUBROUTINE SCALE_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CAT_NARRAY(na, nb)

  IMPLICIT NONE

  ! I/O
  TYPE(t_narray), INTENT(INOUT) :: na
  TYPE(t_narray), INTENT(in)    :: nb

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'CAT_NARRAY'
  TYPE(t_narray)              :: nh
  INTEGER                     :: vtype1, vtype2
  INTEGER                     :: n,m,i

  vtype1 = na%type
  vtype2 = nb%type

  IF (vtype2 == VTYPE_UNDEF) THEN
     CALL RGMSG(substr, RGMLE, 'N-ARRAY TO BE APPENDED MUST BE DEFINED !')
  ELSE
     IF (vtype1 == VTYPE_UNDEF) THEN
        CALL COPY_NARRAY(na, nb)
        RETURN
     ELSE
        IF (vtype1 /= vtype2) THEN
           CALL RGMSG(substr, RGMLE, 'N-ARRAYS MUST BE OF SAME TYPE !')
        END IF
     END IF
     IF ((na%n /= 1).OR.(nb%n /= 1)) THEN
        CALL RGMSG(substr, RGMLE, 'N-ARRAYS MUST BE 1-DIMENSIONAL !')
     END IF
  END IF

  n = na%dim(1)
  m = nb%dim(1)
  CALL INIT_NARRAY(nh, 1, (/ n+m /), vtype2)

  SELECT CASE(vtype2)
  CASE(VTYPE_REAL)
     DO i=1, n
        nh%vr(i) = na%vr(i)
     END DO
     DO i=1, m
        nh%vr(n+i) = nb%vr(i)
     END DO
  CASE(VTYPE_DOUBLE)
     DO i=1, n
        nh%vd(i) = na%vd(i)
     END DO
     DO i=1, m
        nh%vd(n+i) = nb%vd(i)
     END DO
  CASE(VTYPE_INT)
     DO i=1, n
        nh%vi(i) = na%vi(i)
     END DO
     DO i=1, m
        nh%vi(n+i) = nb%vi(i)
     END DO
  CASE(VTYPE_BYTE)
     DO i=1, n
        nh%vb(i) = na%vb(i)
     END DO
     DO i=1, m
        nh%vb(n+i) = nb%vb(i)
     END DO
  CASE(VTYPE_CHAR)
     DO i=1, n
        nh%vc(i) = na%vc(i)
     END DO
     DO i=1, m
        nh%vc(n+i) = nb%vc(i)
     END DO
  CASE(VTYPE_UNDEF)
     ! NOTHING TO DO FOR UNDEFINED ARRAYS
  CASE DEFAULT
     ! NOTHING TO DO FOR UNRECOGNIZED ARRAYS
  END SELECT

  CALL COPY_NARRAY(na, nh)
  ! CLEAN UP
  CALL INIT_NARRAY(nh)

END SUBROUTINE CAT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
INTEGER FUNCTION POSITION(vdim, vec)

! THIS FUNCTION CALCULATES THE POSITION NUMBER IN A
! 1-D (LINEAR) ARRAY, GIVEN THAT THE ARRAY SHOULD BE INTERPRETED
! AS N-D ARRAY WITH DIMENSIONS
! dim = (d1, d2, d3, ..., dN)
! OF THE ELEMENT
! vec = (v1, v2, v3, ..., vN)
!

  IMPLICIT NONE

  ! I/O
  INTEGER, DIMENSION(:), INTENT(IN) :: vdim
  INTEGER, DIMENSION(:), INTENT(IN) :: vec

  ! LOCAL
  INTEGER :: i
  INTEGER :: n
  INTEGER :: dacc

  IF (SIZE(vdim) /= SIZE(vec)) THEN
     POSITION = 0
     RETURN
  END IF

  DO i=1, SIZE(vdim)
     IF (vec(i) > vdim(i)) THEN
        POSITION = 0
        RETURN
     END IF
  END DO

  n = vec(1)
  dacc = 1
  DO i=2,SIZE(vdim)
     dacc = dacc*vdim(i-1)
     n = n + dacc*(vec(i)-1)
  END DO

  POSITION = n

END FUNCTION POSITION
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE ELEMENT(dim, n, vec)

! THIS SUBROUTINE CALCULATES THE ELEMENT VECTOR
! vec = (v1, v2, v3, ..., vN)
! OF THE ELEMENT WITH POSITION n IN A
! 1-D (LINEAR) ARRAY, GIVEN THAT THE ARRAY SHOULD BE INTERPRETED
! AS N-D ARRAY WITH DIMENSIONS
! dim = (d1, d2, d3, ..., dN)
!

  IMPLICIT NONE

  ! I/O
  INTEGER, DIMENSION(:), INTENT(IN)  :: dim   ! dimension vector
  INTEGER,               INTENT(IN)  :: n     ! element in linear array
  INTEGER, DIMENSION(:), POINTER     :: vec   ! element vector

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'ELEMENT'
  INTEGER                             :: m    ! COPY of n
  INTEGER                             :: i    ! counter
  INTEGER , DIMENSION(:), ALLOCATABLE :: dacc
  INTEGER                             :: l    ! length of dim
  INTEGER                             :: status

  m = n
  l = SIZE(dim)
  IF (ASSOCIATED(vec)) THEN
     IF (SIZE(vec) > 0) THEN
        DEALLOCATE(vec, STAT=status)
        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(vec)
  END IF
  ALLOCATE(vec(l), STAT=status)
  CALL ERRMSG(substr,status,2)
  vec(:) = 0

  ALLOCATE(dacc(l))
  dacc(1) = 1
  DO i=2, l
     dacc(i) = dacc(i-1)*dim(i-1)
  END DO

  IF (m > dacc(l)*dim(l)) RETURN

  DO i=l, 2, -1
     vec(i) = (m-1)/dacc(i)+1
     m = m - (vec(i)-1)*dacc(i)
  END DO
  vec(1) = m

  DEALLOCATE(dacc, stat=STATUS)
  CALL ERRMSG(substr,status,3)

END SUBROUTINE ELEMENT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_I(data,idx,ileft,iright)

  IMPLICIT NONE

  ! I/O
  INTEGER (I8), DIMENSION(:), INTENT(INOUT)        :: data   ! data to sort
  INTEGER (I8), DIMENSION(:), POINTER              :: idx  ! index list
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ileft, iright

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'QSORT_I'
  INTEGER (I8)  :: temp            ! temporal data
  INTEGER       :: n               ! LENGTH OF LIST
  INTEGER (I8)  :: left, right
  INTEGER (I8)  :: i, last, apu
  INTEGER (I8)  :: ti              ! temporal index
  INTEGER       :: status

  ! INIT
  n = SIZE(data)
  IF (.NOT.ASSOCIATED(idx)) THEN
     ALLOCATE(idx(n), STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, n
        idx(i) = i
     END DO
  END IF

  IF (PRESENT(ileft)) THEN
     left = ileft
  ELSE
     left = 1
  END IF
  IF (PRESENT(iright)) THEN
     right = iright
  ELSE
     right = n
  END IF

  IF(left >= right) RETURN

  apu = (left+right)/2

  temp = data(left)
  data(left) = data(apu)
  data(apu) = temp

  ti = idx(left)
  idx(left) = idx(apu)
  idx(apu) = ti

  last = left

  DO i=left+1,right
     IF(data(i) < data(left)) THEN
        last = last+1
        temp = data(last)
        data(last) = data(i)
        data(i) = temp

        ti = idx(last)
        idx(last) = idx(i)
        idx(i) = ti
     ENDIF
  END DO

  temp = data(left)
  data(left) = data(last)
  data(last) = temp

  ti = idx(left)
  idx(left) = idx(last)
  idx(last) = ti

  CALL QSORT_I(data,idx,left,last-1)
  CALL QSORT_I(data,idx,last+1,right)

  RETURN

END SUBROUTINE QSORT_I
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_B(data,idx,ileft,iright)

  IMPLICIT NONE

  ! I/O
  INTEGER (I4), DIMENSION(:), INTENT(INOUT)        :: data   ! data to sort
  INTEGER (I8), DIMENSION(:), POINTER              :: idx  ! index list
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ileft, iright

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'QSORT_B'
  INTEGER (I4)  :: temp            ! temporal data
  INTEGER       :: n               ! LENGTH OF LIST
  INTEGER (I8)  :: left, right
  INTEGER (I8)  :: i, last, apu
  INTEGER (I8)  :: ti              ! temporal index
  INTEGER       :: status

  ! INIT
  n = SIZE(data)
  IF (.NOT.ASSOCIATED(idx)) THEN
     ALLOCATE(idx(n), STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, n
        idx(i) = i
     END DO
  END IF

  IF (PRESENT(ileft)) THEN
     left = ileft
  ELSE
     left = 1
  END IF
  IF (PRESENT(iright)) THEN
     right = iright
  ELSE
     right = n
  END IF

  IF(left >= right) RETURN

  apu = (left+right)/2

  temp = data(left)
  data(left) = data(apu)
  data(apu) = temp

  ti = idx(left)
  idx(left) = idx(apu)
  idx(apu)  = ti

  last = left

  DO i=left+1,right
     IF(data(i) < data(left)) THEN
        last = last+1
        temp = data(last)
        data(last) = data(i)
        data(i) = temp
        ti = idx(last)
        idx(last) = idx(i)
        idx(i) = ti
     ENDIF
  END DO

  temp = data(left)
  data(left) = data(last)
  data(last) = temp

  ti = idx(left)
  idx(left) = idx(last)
  idx(last) = ti

  CALL QSORT_B(data,idx,left,last-1)
  CALL QSORT_B(data,idx,last+1,right)

  RETURN

END SUBROUTINE QSORT_B
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_R(data,idx,ileft,iright)

  IMPLICIT NONE

  ! I/O
  REAL (SP),    DIMENSION(:), INTENT(INOUT)        :: data   ! data to sort
  INTEGER (I8), DIMENSION(:), POINTER              :: idx  ! index list
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ileft, iright

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'QSORT_R'
  REAL (SP)    :: temp            ! temporal data
  INTEGER      :: n               ! LENGTH OF LIST
  INTEGER (I8) :: left, right
  INTEGER (I8) :: i, last, apu
  INTEGER (I8) :: ti              ! temporal index
  INTEGER      :: status

  ! INIT
  n = SIZE(data)
  IF (.NOT.ASSOCIATED(idx)) THEN
     ALLOCATE(idx(n),STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, n
        idx(i) = i
     END DO
  END IF

  IF (PRESENT(ileft)) THEN
     left = ileft
  ELSE
     left = 1
  END IF
  IF (PRESENT(iright)) THEN
     right = iright
  ELSE
     right = n
  END IF

  IF(left >= right) RETURN

  apu = (left+right)/2

  temp = data(left)
  data(left) = data(apu)
  data(apu) = temp

  ti = idx(left)
  idx(left) = idx(apu)
  idx(apu) = ti

  last = left

  DO i=left+1,right
     IF(data(i) < data(left)) THEN
        last = last+1
        temp = data(last)
        data(last) = data(i)
        data(i) = temp
        ti = idx(last)
        idx(last) = idx(i)
        idx(i) = ti
      ENDIF
  END DO

  temp = data(left)
  data(left) = data(last)
  data(last) = temp

  ti = idx(left)
  idx(left) = idx(last)
  idx(last) = ti

  CALL QSORT_R(data,idx,left,last-1)
  CALL QSORT_R(data,idx,last+1,right)

  RETURN

END SUBROUTINE QSORT_R
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_D(data,idx,ileft,iright)

  IMPLICIT NONE

  ! I/O
  REAL (DP),    DIMENSION(:), INTENT(INOUT)        :: data   ! data to sort
  INTEGER (I8), DIMENSION(:), POINTER              :: idx  ! index list
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ileft, iright

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'QSORT_D'
  REAL (DP)    :: temp            ! temporal data
  INTEGER      :: n               ! LENGTH OF LIST
  INTEGER (I8) :: left, right
  INTEGER (I8) :: i, last, apu
  INTEGER (I8) :: ti              ! temporal index
  INTEGER      :: status

  ! INIT
  n = SIZE(data)
  IF (.NOT.ASSOCIATED(idx)) THEN
     ALLOCATE(idx(n), STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, n
        idx(i) = i
     END DO
  END IF

  IF (PRESENT(ileft)) THEN
     left = ileft
  ELSE
     left = 1
  END IF
  IF (PRESENT(iright)) THEN
     right = iright
  ELSE
     right = n
  END IF

  IF(left >= right) RETURN

  apu = (left+right)/2

  temp = data(left)
  data(left) = data(apu)
  data(apu) = temp

  ti = idx(left)
  idx(left) = idx(apu)
  idx(apu) = ti

  last = left

  DO i=left+1,right
     IF(data(i) < data(left)) THEN
        last = last+1
        temp = data(last)
        data(last) = data(i)
        data(i) = temp
        ti = idx(last)
        idx(last) = idx(i)
        idx(i) = ti
     ENDIF
  END DO

  temp = data(left)
  data(left) = data(last)
  data(last) = temp

  ti = idx(left)
  idx(left) = idx(last)
  idx(last) = ti

  CALL QSORT_D(data,idx,left,last-1)
  CALL QSORT_D(data,idx,last+1,right)

  RETURN

END SUBROUTINE QSORT_D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE ERRMSG(routine, status, pos)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)  :: routine
  INTEGER,          INTENT(IN)  :: status
  INTEGER,          INTENT(IN)  :: pos

  IF (status == 0) THEN
     RETURN
  ELSE
     CALL RGMSG(routine, RGMLE,  'ERROR STATUS ',status,' ', .false.)
     CALL RGMSG(routine, RGMLEC, 'AT POSITION ',pos,' !')
  END IF

END SUBROUTINE ERRMSG
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_C(routine, level, c, lstop)

#if ! defined(NOMPI)
  USE MESSY_MAIN_GRID_MPI,  ONLY: GRID_ABORT
#endif

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  LOGICAL           :: llstop
  INTEGER           :: iou
  LOGICAL           :: opened
  CHARACTER(LEN=22) :: errfname = ''
  CHARACTER(LEN=22) :: endfname = ''

  ! INIT
  IF (PRESENT(lstop)) THEN
     llstop = lstop
  ELSE
     IF ((level == RGMLE).OR.(level == RGMLEC)) THEN
        llstop = .true.      ! STOP ON ERROR
     ELSE
        llstop = .false.
     END IF
  END IF

  IF ((level == RGMLE).OR.(level == RGMLEC)) THEN
     DO iou=100,300
        INQUIRE(unit=iou,opened=opened)
        IF (.NOT.opened) EXIT
     END DO
     IF (MY_RANK < 0) THEN
        errfname = 'ERROR.grid_netcdf'
        endfname = 'END'
     ELSE
        WRITE(errfname,'(a17,a1,i4.4)') 'ERROR.grid_netcdf','.',MY_RANK
        WRITE(endfname,'(a3,i4.4)') 'END',MY_RANK
     ENDIF
  END IF

  SELECT CASE(level)
  CASE(RGMLE)   ! ERROR MESSAGE
     IF (IAND(MSGMODE, MSGMODE_E) == MSGMODE_E) THEN
        WRITE(*,*) '*** ',TRIM(routine),' ERROR: '
        WRITE(*,*) '    ',TRIM(c)
        OPEN(iou, FILE=TRIM(errfname), STATUS='UNKNOWN')
        WRITE(iou,*) '*** ',TRIM(routine),' ERROR: '
        WRITE(iou,*) '    ',TRIM(c)
        CLOSE(iou)
        OPEN(iou, FILE=TRIM(endfname), STATUS='UNKNOWN')
        WRITE(iou,*) 'ERROR! See ',TRIM(errfname)
        CLOSE(iou)
     END IF
  CASE(RGMLEC) ! ERROR MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_E) == MSGMODE_E) THEN
        WRITE(*,*) '    ',TRIM(c)
        OPEN(iou, FILE=TRIM(errfname), STATUS='OLD',POSITION='APPEND')
        WRITE(iou,*) '    ',TRIM(c)
        CLOSE(iou)
     END IF
  CASE(RGMLVL)  ! LITTLE VERBOSE
     IF (IAND(MSGMODE, MSGMODE_VL) == MSGMODE_VL) THEN
        WRITE(*,*) TRIM(c)
     END IF
  CASE(RGMLVLC) ! LITTLE VERBOSE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_VL) == MSGMODE_VL) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLW) ! WARNING MESSAGE
     IF (IAND(MSGMODE, MSGMODE_W) == MSGMODE_W) THEN
        WRITE(*,*) '+++ ',TRIM(routine),' WARNING: '
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLWC) ! WARNING MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_W) == MSGMODE_W) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLVM)  ! MEDIUM VERBOSE
     IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
        WRITE(*,*) TRIM(c)
     END IF
  CASE(RGMLVMC) ! MEDIUM VERBOSE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLI)  ! INFO MESSAGE
     IF (IAND(MSGMODE, MSGMODE_I) == MSGMODE_I) THEN
        WRITE(*,*) '=== ',TRIM(routine),' INFO: '
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLIC) ! INFO MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_I) == MSGMODE_I) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE DEFAULT
  END SELECT

#if ! defined(NOMPI)
  IF (llstop) CALL GRID_ABORT('grid')
#else
  IF (llstop) STOP
#endif

END SUBROUTINE RGMSG_C
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_I(routine, level, c1, i, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  INTEGER,          INTENT(IN)            :: i
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: istr = ''

  WRITE(istr,*) TRIM(c1),i,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(istr), lstop)

END SUBROUTINE RGMSG_I
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_IA(routine, level, c1, i, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*),      INTENT(IN)            :: routine
  INTEGER,               INTENT(IN)            :: level
  CHARACTER(LEN=*),      INTENT(IN)            :: c1
  INTEGER, DIMENSION(:), INTENT(IN)            :: i
  CHARACTER(LEN=*),      INTENT(IN)            :: c2
  LOGICAL,               INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: istr = ''

  WRITE(istr,*) TRIM(c1),i,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(istr), lstop)

END SUBROUTINE RGMSG_IA
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_R(routine, level, c1, r, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  REAL,             INTENT(IN)            :: r
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: rstr = ''

  WRITE(rstr,*) TRIM(c1),r,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(rstr), lstop)

END SUBROUTINE RGMSG_R
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_D(routine, level, c1, d, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  REAL(dp),         INTENT(IN)            :: d
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: rstr = ''

  WRITE(rstr,*) TRIM(c1),d,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(rstr), lstop)

END SUBROUTINE RGMSG_D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE MAIN_GRID_SET_MESSAGEMODE(SETVALUE)

  ! set MESSAGE mode

  IMPLICIT NONE

  INTEGER, INTENT(IN), OPTIONAL :: SETVALUE

  IF (PRESENT(SETVALUE)) THEN
     MSGMODE = SETVALUE
  ELSE
     ! RESET MESSAGE MODUS
     MSGMODE = MSGMODE_S + MSGMODE_E + MSGMODE_VL  &
          + MSGMODE_W + MSGMODE_VM + MSGMODE_I
  ENDIF

END SUBROUTINE MAIN_GRID_SET_MESSAGEMODE
! ------------------------------------------------------------------

! ------------------------------------------------------------------
 REAL(dp) FUNCTION NARRAYMAXVAL (array)

   IMPLICIT NONE

   TYPE(t_narray), INTENT(IN) :: array

   SELECT CASE (array%type)
   CASE(VTYPE_UNDEF)
      ! nothing to print
   CASE(VTYPE_INT)
      NARRAYMAXVAL = REAL(MAXVAL(array%vi), dp)
   CASE(VTYPE_REAL)
      NARRAYMAXVAL = REAL(MAXVAL(array%vr), dp)
   CASE(VTYPE_DOUBLE)
      NARRAYMAXVAL = MAXVAL(array%vd)
   CASE DEFAULT
      NARRAYMAXVAL = -9999.99_dp
   END SELECT

   RETURN

 END FUNCTION NARRAYMAXVAL
! ------------------------------------------------------------------

! ------------------------------------------------------------------
 REAL(dp) FUNCTION NARRAYMINVAL (array)

   IMPLICIT NONE

   TYPE(t_narray), INTENT(IN) :: array

   SELECT CASE (array%type)
   CASE(VTYPE_UNDEF)
      ! nothing to print
   CASE(VTYPE_INT)
      NARRAYMINVAL = REAL(MINVAL(array%vi), dp)
   CASE(VTYPE_REAL)
      NARRAYMINVAL = REAL(MINVAL(array%vr), dp)
   CASE(VTYPE_DOUBLE)
      NARRAYMINVAL = MINVAL(array%vd)
   CASE DEFAULT
      NARRAYMINVAL = 9999.99_dp
   END SELECT

   RETURN

 END FUNCTION NARRAYMINVAL
! ------------------------------------------------------------------

! ------------------------------------------------------------------
 SUBROUTINE multinetcdf(status, fin, tname, fout, tl, tf, mc)

   USE messy_main_tools, ONLY: find_next_free_unit, dir_and_base

   IMPLICIT NONE

   INTRINSIC SUM

   ! I/O
   INTEGER,                      INTENT(OUT)   :: status
   CHARACTER(LEN=GRD_MAXSTRLEN), INTENT(IN)    :: fin   ! filename
   CHARACTER(LEN=GRD_MAXSTRLEN), INTENT(IN)    :: tname ! name of time variable
   CHARACTER(LEN=GRD_MAXSTRLEN), INTENT(OUT)   :: fout  ! filename
   INTEGER,                      INTENT(IN)    :: tl   ! time step index in list
   INTEGER,                      INTENT(OUT)   :: tf   ! time step index in file
   TYPE(t_multinc),              INTENT(INOUT) :: mc   ! file list

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER :: substr = 'multinetcdf'
   LOGICAL :: lex
   INTEGER :: ncstat
   INTEGER :: ncid
   INTEGER :: unlimitedDimId
   INTEGER :: iou
   INTEGER :: fstat
   INTEGER :: i
   INTEGER :: mt
   CHARACTER(LEN=GRD_MAXSTRLEN) :: path
   CHARACTER(LEN=GRD_MAXSTRLEN) :: base
   TYPE(t_ncvar)                :: var

   status = 0

   not_yet_scanned: IF (.NOT. mc%l) THEN
      ! check, if fin is present
      INQUIRE(file=TRIM(fin), exist=lex)
      IF (.NOT. lex) THEN
         CALL RGMSG(substr, RGMLE, &
              'FILE '//TRIM(fin)//' NOT FOUND!')
         status = 1
         RETURN
      END IF
      !
      ! check, if fin is netcdf file or multi-netcdf descriptor file
      ncstat = NF90_OPEN(TRIM(fin), NF90_NOWRITE, ncid)
      is_nc: IF (ncstat == NF90_NOERR) THEN
         ! file is netcdf file
         fout  = TRIM(fin)
         tf    = tl
         mc%nf = 0
         CALL NFERR(substr, NF90_CLOSE(ncid), 1)
         status = 0
         RETURN
      ELSE
         ! separate path and basenam
         CALL dir_and_base(fin, path, base)
         !
         ! file is NOT a netcdf file -> assuming multi-descriptor file
         iou = find_next_free_unit(100, 200)
         OPEN(unit=iou,file=TRIM(fin))
         ! COUNT LINES
         mc%nf = 0
         lines1: DO
            READ(iou, *, IOSTAT=fstat)
            IF (fstat < 0) EXIT
            mc%nf = mc%nf + 1
         END DO lines1
         ! rewind to beginning of file
         REWIND(unit=iou)
         !
         ! allocate memory
         ALLOCATE(mc%flist(mc%nf))
         ALLOCATE(mc%nt(mc%nf))
         !
         ! read data
         lines2: DO i=1, mc%nf
            READ(iou, *, IOSTAT=fstat) mc%flist(i)
            IF (fstat /= 0) THEN
               CLOSE(iou)
               CALL RGMSG(substr, RGMLE, &
                    'READ ERROR (FILE '//TRIM(fin)//')!')
               status = 2
               RETURN
            END IF
         END DO lines2
         !
         ! close file
         CLOSE(iou)
         !
         ! set full path and
         ! read number of time steps from netcdf files
         file_loop: DO i=1, mc%nf
            mc%flist(i) = TRIM(path)//'/'//TRIM(mc%flist(i))
            timevar: IF (TRIM(ADJUSTL(tname)) == '') THEN
               ! TRY UNLIMITED DIMENSION
               CALL RGMSG(substr, RGMLI, &
                    'empty time variable, searching for UNLIMITED'//&
                    &' dimension in file '//TRIM(mc%flist(i)) )
               CALL NFERR(substr &
                    , NF90_OPEN(TRIM(mc%flist(i)), NF90_NOWRITE, ncid) &
                    , 2)
               CALL NFERR(substr &
                    , NF90_INQUIRE(ncid, unlimitedDimId=unlimitedDimId) &
                    , 3)
               IF (unlimitedDimId < 0) THEN
                  CALL NFERR(substr, NF90_CLOSE(ncid), 4)
                  CALL RGMSG(substr, RGMLE, &
                    'no UNLIMITED dimension in file'//TRIM(mc%flist(i)) )
                  status = 3
                  RETURN
               ELSE
                  CALL NFERR(substr &
                       , NF90_INQUIRE_DIMENSION(&
                           ncid, unlimitedDimId, len=mc%nt(i)) &
                       , 5)
                  CALL NFERR(substr, NF90_CLOSE(ncid), 6)
               END IF
            ELSE
               ! READ timem VARIABLE
               CALL RGMSG(substr, RGMLI, &
                    'reading variable '//TRIM(ADJUSTL(tname))//' from file '//&
                    &TRIM(mc%flist(i)) )
               CALL INIT_NCVAR(var)
               CALL IMPORT_NCVAR(var, varname = TRIM(ADJUSTL(tname)) &
                    , file = TRIM(mc%flist(i)))
               IF (var%ndims /= 1) THEN
                  CALL RGMSG(substr, RGMLE, &
                    TRIM(ADJUSTL(tname))//' in file'//&
                    &TRIM(mc%flist(i))//' is multi-dimensional')
                  status = 4
                  RETURN
               ELSE
                  mc%nt(i) = var%dim(1)%len
               END IF
               CALL INIT_NCVAR(var)
            END IF timevar
            !
            CALL RGMSG(substr, RGMLI, &
                 ' ... found ',mc%nt(i),' steps in file '//TRIM(mc%flist(i)) )
            !
         END DO file_loop
         !
      END IF is_nc

      mc%l = .TRUE.

   END IF not_yet_scanned

   ! ALWAYS AFTER SCAN IS COMPLETE:

   IF (mc%nf == 0) THEN
      ! file is netcdf file, see above
      fout  = TRIM(fin)
      tf    = tl
      status = 0
      RETURN
   END IF

   mt = 0 ! counter of cumulative time steps
   DO i=1, mc%nf
      mt = mt + mc%nt(i)
      IF (tl <= mt) EXIT
   END DO
   IF (i > mc%nf) THEN ! requested time step index out of range
      CALL RGMSG(substr, RGMLE, &
           'requested time step', tl, 'is out of range'//&
           &' in file '//TRIM(fin))
      status = 5
      RETURN
   END IF
   fout = TRIM(mc%flist(i))
   tf   = tl - SUM(mc%nt(1:i-1))

 END SUBROUTINE multinetcdf
! ------------------------------------------------------------------

! ------------------------------------------------------------------
 SUBROUTINE init_multinc(mc)

   IMPLICIT NONE

   ! I/O
   TYPE(t_multinc), INTENT(INOUT) :: mc

   mc%l  = .FALSE.
   mc%nf = 0
   IF (ASSOCIATED(mc%flist)) THEN
      DEALLOCATE(mc%flist); NULLIFY(mc%flist)
   END IF
   IF (ASSOCIATED(mc%nt)) THEN
      DEALLOCATE(mc%nt); NULLIFY(mc%nt)
   END IF

 END SUBROUTINE init_multinc
! ------------------------------------------------------------------

! ******************************************************************
END MODULE MESSY_MAIN_GRID_NETCDF
! ******************************************************************
