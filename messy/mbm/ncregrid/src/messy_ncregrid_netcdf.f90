! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_NETCDF
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  USE netcdf
  USE MESSY_NCREGRID_BASE

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, SIZE, TRIM, PRESENT, MIN, PRODUCT &
             , LEN_TRIM, INT, MAXVAL, MINVAL, REAL, NULL

  PRIVATE   :: ASSOCIATED, SIZE, TRIM, PRESENT, MIN, PRODUCT &
             , LEN_TRIM, INT, MAXVAL, MINVAL, REAL, NULL

  INTEGER, PARAMETER :: GRD_MAXSTRLEN = 200
  INTEGER, PARAMETER :: NULL_DIMID = -1
  INTEGER, PARAMETER :: NULL_VARID = -1
  INTEGER, PARAMETER :: NULL_XTYPE = -1

  TYPE ncdim
     CHARACTER(LEN=GRD_MAXSTRLEN) :: name = ''          ! name of dimension
     INTEGER                      :: id   = NULL_DIMID  ! dimension ID
     INTEGER                      :: len  = 0           ! length of dimension
     LOGICAL                      :: fuid = .false.     ! flag for "UNLIMITED"
                                                        ! ... dimension
     INTEGER                      :: varid = NULL_VARID ! variable ID,
                                                        ! ... if coordinate var
  END TYPE ncdim
  
  TYPE ncatt
     CHARACTER(LEN=GRD_MAXSTRLEN) :: name  = ''         ! attribute name
     INTEGER                      :: num   = 0          ! number of attribute
     INTEGER                      :: xtype = NULL_XTYPE ! type of attribute
     INTEGER                      :: varid = NULL_VARID ! variable ID or
                                                        ! ... NF90_GLOBAL
     INTEGER                      :: len   = 0          ! length of attribute
     TYPE(narray)                 :: dat                ! content of attribute
  END TYPE ncatt

  TYPE ncvar
     CHARACTER(LEN=GRD_MAXSTRLEN)        :: name  = ''         ! variable name
     INTEGER                             :: id    = NULL_VARID ! variable ID
     INTEGER                             :: xtype = NULL_XTYPE ! type of
                                                               ! ... variable
     INTEGER                             :: ndims = 0          ! number of
                                                               ! ... dimensions
     TYPE (ncdim), DIMENSION(:), POINTER :: dim => NULL()      ! netCDF
                                                               ! ... dimensions
     INTEGER                             :: uid   = NULL_DIMID ! unlimited
                                                               ! ... dim ID
     INTEGER                             :: ustep = 0  ! step along unlim. ID
     INTEGER                             :: natts = 0  ! number of attributes
     TYPE (ncatt), DIMENSION(:), POINTER :: att => NULL() ! list of attributes
     TYPE(narray)                        :: dat           ! content of variable
  END TYPE ncvar

!! NOTE: DOES NOT WORK PROPERLY FOR SOME COMPILERS ...
!  INTERFACE ASSIGNMENT (=)
!     MODULE PROCEDURE COPY_NCDIM
!     MODULE PROCEDURE COPY_NCATT
!     MODULE PROCEDURE COPY_NCVAR
!  END INTERFACE

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE INIT_NCDIM(dim)
  
  IMPLICIT NONE

  ! I/O
  TYPE (ncdim), INTENT(INOUT) :: dim

  dim%name  = ''
  dim%id    = NULL_DIMID
  dim%len   = 0
  dim%fuid  = .false.
  dim%varid = NULL_VARID

END SUBROUTINE INIT_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_NCATT(att)

  IMPLICIT NONE
  
  ! I/O
  TYPE (ncatt), INTENT(INOUT) :: att

  att%name  = ''
  att%num   = 0
  att%xtype = NULL_XTYPE
  att%len   = 0
  att%varid = NULL_VARID

  CALL INIT_NARRAY(att%dat)

END SUBROUTINE INIT_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_NCVAR(var)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar), INTENT(INOUT) :: var

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'INIT_NCVAR'
  INTEGER :: i 
  INTEGER :: status

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

  IMPLICIT NONE

  ! I/O
  TYPE (ncdim), INTENT(INOUT) :: dd   ! destination
  TYPE (ncdim), INTENT(IN)    :: ds   ! source

  dd%name  = TRIM(ds%name)
  dd%id    = ds%id
  dd%len   = ds%len
  dd%fuid  = ds%fuid
  dd%varid = ds%varid  

END SUBROUTINE COPY_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_NCATT(ad, as)

  ! I/O
  TYPE (ncatt), INTENT(INOUT) :: ad   ! destination
  TYPE (ncatt), INTENT(IN)    :: as   ! source

  ad%name  = TRIM(as%name)
  ad%num   = as%num
  ad%xtype = as%xtype
  ad%len   = as%len
  ad%varid = as%varid

  CALL COPY_NARRAY(ad%dat, as%dat)

END SUBROUTINE COPY_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_NCVAR(vd, vs)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar), INTENT(INOUT) :: vd   ! destination
  TYPE (ncvar), INTENT(IN)    :: vs   ! source

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

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar), INTENT(IN)  :: var

  QDEF_NCVAR = (var%dat%type /= VTYPE_UNDEF)
  ! QDEF_NCVAR = (var%id /= NULL_VARID)

END FUNCTION QDEF_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
INTEGER  FUNCTION QCMP_NCDIM(d1, d2)

  IMPLICIT NONE

  ! I/O
  TYPE (ncdim), INTENT(IN) :: d1, d2  ! dimensions to compare

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

  IMPLICIT NONE

  ! I/O
  TYPE (ncdim),     INTENT(INOUT)         :: dim      ! dimension
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: dimname  ! name of dimension
  INTEGER,          INTENT(IN),  OPTIONAL :: dimid    ! dimension ID
  INTEGER,          INTENT(IN),  OPTIONAL :: ncid     ! netCDF file ID
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: file     ! filename

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMPORT_NCDIM'
  INTEGER :: lncid    ! local netCDF file-ID
  INTEGER :: uid      ! unlimited ID
  INTEGER :: status   ! netCDF status

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
     CALL nferr(substr                                 &
          ,nf90_open(TRIM(file),NF90_NOWRITE,lncid)    & 
          ,1)
  END IF

  IF (PRESENT(dimname)) THEN
     ! GET DIM-ID from DIMNAME
     dim%name = TRIM(dimname)
     IF (TRIM(dim%name) /= '') THEN
        CALL nferr(substr                                 &
             ,nf90_inq_dimid(lncid,TRIM(dim%name),dim%id) &
             ,2)
     ELSE
        CALL RGMSG(substr, RGMLE, 'INVALID NAME OF DIMENSION !')
     END IF
  ELSE
     dim%id = dimid
     ! GET DIM-NAME from DIM-ID
     CALL nferr(substr                                           &
          ,nf90_Inquire_Dimension(lncid, dim%id, name=dim%name)  &
          ,3)
  END IF

  ! GET DIM-LEN
  CALL nferr(substr                                              &
             ,nf90_Inquire_Dimension(lncid, dim%id, len=dim%len) &
             ,4)
  
  ! CHECK IF UNLIMITED
  CALL nferr(substr                                   &
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
     CALL nferr(substr       &
          ,nf90_close(lncid) &
          ,6)
  END IF

END SUBROUTINE IMPORT_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_NCDIM(dim, file, ncid)
  
  IMPLICIT NONE

  ! I/O
  TYPE (ncdim),     INTENT(INOUT)           :: dim
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
!     CALL nferr(substr        &
!          ,nf90_redef(lncid)  & 
!          ,1)
     status = nf90_redef(lncid)

  ELSE
     ! CHECK IF FILE EXISTS
     INQUIRE(file=TRIM(file), exist=lex)
     IF (lex) THEN   ! FILE EXISTS
        CALL RGMSG(substr, RGMLI, &
             'netCDF-FILE '''//TRIM(file)//''' EXISTS ! OPENING ...')
        CALL nferr(substr                            &
             ,nf90_open(TRIM(file),NF90_WRITE,lncid) & 
             ,2)
        CALL nferr(substr        &
             ,nf90_redef(lncid)  & 
             ,3)
     ELSE            ! NEW FILE
        CALL RGMSG(substr, RGMLI, &
             'netCDF-FILE '''//TRIM(file)//''' DOES NOT EXIST ! NEW FILE ...')
        CALL nferr(substr                                &
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
     CALL nferr(substr                              &
          ,nf90_Inquire(lncid, unlimitedDimID=uid)  &
          ,5)
     IF (uid == -1) uid = NULL_DIMID
     ! GET (CURRENT) LENGTH
     CALL nferr(substr                                        &
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
        CALL nferr(substr                       &
             ,nf90_def_dim(lncid,TRIM(dim%name) &
             ,NF90_UNLIMITED, dim%id)           &
             ,7)
     ELSE
        CALL RGMSG(substr, RGMLIC, ' ... NEW (CONSTANT LENGTH)')
        CALL nferr(substr                       &
             ,nf90_def_dim(lncid,TRIM(dim%name) &
             ,dim%len, dim%id)                  &
             ,8)
     END IF
  END IF                          ! DIM exists ?

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL nferr(substr        &
          ,nf90_close(lncid)  &
          ,9)
  END IF

END SUBROUTINE EXPORT_NCDIM
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IMPORT_NCATT(att, varname, varid, attname, attnum  &
                         ,file, ncid                          &
                         ,lnostop)

  IMPLICIT NONE

  ! I/O
  TYPE (ncatt),     INTENT(INOUT)          :: att      ! attribute
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL  :: varname  ! variable name
  INTEGER,          INTENT(IN),  OPTIONAL  :: varid    ! netCDF variable ID
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL  :: attname  ! attribute name
  INTEGER,          INTENT(IN),  OPTIONAL  :: attnum   ! attribute number
  INTEGER,          INTENT(IN),  OPTIONAL  :: ncid     ! netCDF file ID
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL  :: file     ! filename
  LOGICAL,          INTENT(IN),  OPTIONAL  :: lnostop  ! do not stop if
                                                       ! att. does not exist

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMPORT_NCATT'
  INTEGER :: lncid    ! local netCDF file-ID
  INTEGER :: znatts   ! local number of attributes
  CHARACTER(LEN=GRD_MAXSTRLEN)     :: lvarname ! local variable name
  CHARACTER(LEN=GRD_MAXSTRLEN)     :: lattname ! local attribute name
  CHARACTER(LEN=100*GRD_MAXSTRLEN) :: astr = ''
  INTEGER :: RGMLH, RGMLHC    ! RG-message-level switch
  INTEGER :: i
  LOGICAL :: lstp             ! STOP ON ERROR IF TRUE
  INTEGER :: status

  ! CHECK CALL OF SUBROUTINE
  IF ((.NOT.PRESENT(varname)).AND.(.NOT.PRESENT(varid))) THEN
     CALL RGMSG(substr, RGMLE, 'NAME OR ID OF VARIABLE MUST BE SPECIFIED !')
  END IF
  !
  IF ((.NOT.PRESENT(attname)).AND.(.NOT.PRESENT(attnum))) THEN
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
     CALL nferr(substr                              &
          ,nf90_open(TRIM(file),NF90_NOWRITE,lncid) & 
          ,1)
  END IF

  ! GET VAR-ID
  IF (PRESENT(varname)) THEN
     lvarname = TRIM(varname)
     ! GET VAR-ID from VARNAME
     IF (TRIM(lvarname) /= '') THEN
        CALL nferr(substr                                    &
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
        CALL nferr(substr                                            &
             ,nf90_Inquire_Variable(lncid, att%varid, name=lvarname) &
             ,3)
     ELSE
        lvarname = ''
     END IF
  END IF

  ! GET ATTNUM
  IF (PRESENT(attname)) THEN         ! ATTRIBUTE NAME
     lattname = TRIM(attname)
     att%name = TRIM(attname)
     ! GET NUM, TYPE, LEN  from ATTNAME
     IF (TRIM(lattname) /= '') THEN
        status = nf90_Inquire_Attribute(lncid, att%varid   &
                 ,lattname                                 &
                 ,att%xtype, att%len, att%num)
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
              CALL nferr(substr       &
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
     CALL nferr(substr                                           &
          ,nf90_Inquire_Variable(lncid, att%varid, natts=znatts) &
          ,5)
     ! CHECK RANGE
     IF (attnum > znatts) THEN
        CALL RGMSG(substr, RGMLH, &
             'ATTRIBUTE WITH NUMBER',attnum,' ',.false.)
        CALL RGMSG(substr, RGMLHC, &
             'OF VARIABLE WITH ID ',att%varid,' ', .false.)
        CALL RGMSG(substr, RGMLHC, 'DOES NOT EXIST !', lstp)
        !
        ! IF lstp = .false. : JUMP BACK
        ! CLOSE FILE
        IF (.NOT.PRESENT(ncid)) THEN
           CALL nferr(substr       &
                ,nf90_close(lncid) &
                ,6)
        END IF
        RETURN
     END IF

     att%num = attnum
     ! GET ATT-NAME FROM VAR-ID AND NUMBER
     CALL nferr(substr                       &
          ,nf90_inq_attname(lncid,att%varid  &
          ,att%num, att%name)                &
          ,7)
     ! GET ATT-LENGTH and TYPE
     CALL nferr(substr                                        &
          ,nf90_Inquire_Attribute(lncid, att%varid, att%name  &
          ,xtype=att%xtype, len=att%len)                      &
          ,8)
  END IF                           ! ATTRIBUTE NAME OR NUMBER

  ! READ VALUE
  SELECT CASE(att%xtype)
  CASE(NF90_BYTE)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_BYTE )
     CALL nferr(substr                     &
          ,nf90_get_att(lncid, att%varid   &
          ,TRIM(att%name)                  &
          ,values=att%dat%vb)              &
          ,9)
  CASE(NF90_SHORT)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_INT )
     CALL nferr(substr                    &
          ,nf90_get_att(lncid, att%varid  &
          ,TRIM(att%name)                 &
          ,values=att%dat%vi)             &
          ,10)
  CASE(NF90_CHAR)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_CHAR )
     CALL nferr(substr                    &
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
     CALL nferr(substr                    &
          ,nf90_get_att(lncid, att%varid  &
          ,TRIM(att%name)                 &
          ,values=att%dat%vi)             &
          ,12)
  CASE(NF90_FLOAT)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_REAL )
     CALL nferr(substr                    &
          ,nf90_get_att(lncid, att%varid  &
          ,TRIM(att%name)                 &
          ,values=att%dat%vr)             &
          ,13)
  CASE(NF90_DOUBLE)
     CALL INIT_NARRAY(att%dat, 1, (/att%len/), VTYPE_DOUBLE )
     CALL nferr(substr                    &
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
     CALL nferr(substr       &
          ,nf90_close(lncid) &
          ,15)
  END IF

END SUBROUTINE IMPORT_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_NCATT(att, file, ncid, clobber)

  IMPLICIT NONE

  ! I/O
  TYPE (ncatt),     INTENT(INOUT)           :: att     ! attribure
  CHARACTER(LEN=*), INTENT(IN),   OPTIONAL  :: file    ! filename
  INTEGER,          INTENT(IN),   OPTIONAL  :: ncid    ! netCDF file ID
  LOGICAL,          INTENT(IN),   OPTIONAL  :: clobber ! overwrite ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'EXPORT_NCATT'
  INTEGER :: lncid    ! local netCDF file ID
  INTEGER :: zxtype   ! local attribute type
  INTEGER :: zattnum  ! local attribute number
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
!     CALL nferr(substr        &
!          ,nf90_redef(lncid)  & 
!          ,1)
     status = nf90_redef(lncid)
  ELSE
     ! CHECK IF FILE EXISTS
     INQUIRE(file=TRIM(file), exist=lex)
     IF (lex) THEN   ! FILE EXISTS
        ! OPEN FILE
        CALL RGMSG(substr, RGMLI, &
             'FILE '''//TRIM(file)//''' EXISTS ! OPENING ...')
        CALL nferr(substr                            &
             ,nf90_open(TRIM(file),NF90_WRITE,lncid) & 
             ,2)
        CALL nferr(substr       &
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
            ,att%name, xtype=zxtype, attnum=zattnum)
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
          'EXISTS ALREADY WITH NUMBER ',zattnum, ' ')
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

!qqq
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
           CALL nferr(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vb(1:alen))                     &
                ,4)
        CASE(NF90_SHORT)
           CALL nferr(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vi(1:alen))                     &
                ,5)
        CASE(NF90_CHAR)
           IF (alen >= 1) THEN
              CALL nferr(substr                             &
                   ,nf90_put_att(lncid, att%varid, att%name &
                   ,string(att%dat%vc))                     &
                   ,6)
           ELSE
              CALL nferr(substr                             &
                   ,nf90_put_att(lncid, att%varid, att%name &
                   ,' ')                                    &
                   ,6)
           END IF
        CASE(NF90_INT)
           CALL nferr(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vi(1:alen))                     &
                ,7)              
        CASE(NF90_FLOAT)
           CALL nferr(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vr(1:alen))                     &
                ,8)
        CASE(NF90_DOUBLE)
           CALL nferr(substr                             &
                ,nf90_put_att(lncid, att%varid, att%name &
                ,att%dat%vd(1:alen))                     &
                ,9)
        CASE DEFAULT
           CALL RGMSG(substr, RGMLW, 'ATTRIBUTE '''//TRIM(att%name)//'''' )
           CALL RGMSG(substr, RGMLWC, 'OF VARIABLE WITH ID ',att%varid,' ')
           CALL RGMSG(substr, RGMLWC, 'IS OF UNKNOWN TYPE ',att%xtype,' ')
           CALL RGMSG(substr, RGMLWC, 'AND CANNOT BE WRITTEN !')
        END SELECT

!qqq
  END IF ! ATTRIBUTE _FillValue

  END IF                         ! attribute exists ?

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL nferr(substr       &
          ,nf90_close(lncid) &
          ,16)
  END IF

END SUBROUTINE EXPORT_NCATT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IMPORT_NCVAR(var, ustep, varname, varid, file, ncid, setuid)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar),     INTENT(INOUT)           :: var      ! variable
  INTEGER,          INTENT(IN),  OPTIONAL   :: ustep    ! step along unlim. DIM
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL   :: varname  ! variable name
  INTEGER,          INTENT(IN),  OPTIONAL   :: varid    ! variable ID
  INTEGER,          INTENT(IN),  OPTIONAL   :: ncid     ! netCDF file ID
  CHARACTER(LEN=*), INTENT(IN),  OPTIONAL   :: file     ! filename
  INTEGER,          INTENT(IN),  OPTIONAL   :: setuid   ! set this as unlim. ID

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMPORT_NCVAR'
  INTEGER                            :: lncid    ! netCDF file-ID
  INTEGER, DIMENSION(:), ALLOCATABLE :: zdimvec  ! local dim. ID vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: zdimlen  ! local dim. lengths
  INTEGER                            :: noudims  ! number of non-ulim. dim.s
  INTEGER                            :: i,j      ! counter
  INTEGER, DIMENSION(:), ALLOCATABLE :: start, count, stride, map
  INTEGER                            :: dim      ! 1D dimension of variable
  INTEGER                            :: status
  INTEGER                            :: ulen     ! length of unlimited ID
  INTEGER                            :: uid      ! ID of unlimited dimension

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

  ! INIT
  CALL INIT_NCVAR(var)

  ! OPEN FILE
  IF (PRESENT(ncid)) THEN
     lncid = ncid
  ELSE
     ! OPEN FILE
     CALL nferr(substr                              &
          ,nf90_open(TRIM(file),NF90_NOWRITE,lncid) & 
          ,1)
  END IF

  IF (PRESENT(varname)) THEN
     ! GET VAR-ID from VARNAME
     var%name = TRIM(varname)
     IF (TRIM(var%name) /= '') THEN
        CALL nferr(substr                                 &
             ,nf90_inq_varid(lncid,TRIM(var%name),var%id) &
             ,2)
     ELSE
        CALL RGMSG(substr, RGMLE, 'INVALID VARIABLE NAME !')
     END IF
  ELSE
     var%id = varid
     ! GET VAR-NAME from VAR-ID
     CALL nferr(substr                                          &
          ,nf90_Inquire_Variable(lncid, var%id, name=var%name)  &
          ,3)
  END IF

  ! GET INFOs
  CALL nferr(substr                                         &
       ,nf90_Inquire_Variable(lncid, var%id                 & 
       ,xtype=var%xtype, ndims=var%ndims, natts=var%natts ) &
       ,4)

  ! CHECK FOR UNLIMITED DIM ID
  CALL nferr(substr                                   &
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
     CALL nferr(substr                          &
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
     END DO
     DEALLOCATE(zdimvec, STAT=status)
     CALL ERRMSG(substr,status,3)
  END IF

  ! ALLOCATE MEMORY AND GET ATTRIBUTE INFO
  IF (var%natts > 0) THEN
     ALLOCATE(var%att(var%natts), STAT=status)
     CALL ERRMSG(substr,status,4)
     DO i=1, var%natts   ! LOOP OVER ATTRIBUTES
        CALL IMPORT_NCATT(var%att(i),varname=TRIM(var%name) &
                          ,attnum=i,ncid=lncid)
!    ! ALTERNATIVE:
!        CALL IMPORT_NCATT(var%att(i),varid=var%id) &
!                          ,attnum=i,ncid=lncid)
     END DO  ! LOOP OVER ATTRIBUTES
  END IF

  ! ALLOCATE MEMORY AND GET DATA
  ! (IF USTEP IS PRESENT, then ONLY ONE SLICE ALONG THE UNLIMITED DIM)
  ! VECTORS
  dim = 1
  IF (var%ndims > 0) THEN        ! DIM > 0
     ALLOCATE(start(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,5)
     ALLOCATE(count(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,6)
     ALLOCATE(stride(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,7)
     ALLOCATE(map(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,8)
!CDIR NOVECTOR
     DO i=1, var%ndims
        IF ((var%dim(i)%id == var%uid).AND.(var%uid /= NULL_DIMID).AND.   &
             PRESENT(ustep)) THEN
           ! UNLIMITED DIMENSION
           IF (ustep <= ulen) THEN
              start(i) = ustep
              count(i) = 1
              stride(i) = 1
           ELSE
              CALL RGMSG(substr, RGMLE, &
                   'REQUESTET U-ID STEP (',ustep,')', .false.)
              CALL RGMSG(substr, RGMLEC, &
                   'IS LARGER THAN AVAILABLE (',ulen,') !')
           END IF
        ELSE
           ! LIMITED DIMENSION
           start(i) = 1
           count(i) = var%dim(i)%len
           stride(i) = 1
        END IF
        dim = dim * count(i)
        IF (i == 1) THEN
           map(i) = 1
        ELSE
           map(i) = map(i-1)*count(i-1)
        END IF
     END DO
  ELSE
     ALLOCATE(start(1), STAT=status)
     CALL ERRMSG(substr,status,9)
     ALLOCATE(count(1), STAT=status)
     CALL ERRMSG(substr,status,10)
     ALLOCATE(stride(1), STAT=status)
     CALL ERRMSG(substr,status,11)
     ALLOCATE(map(1), STAT=status)
     CALL ERRMSG(substr,status,12)
     start(1) = 1
     count(1) = 1
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
     zdimlen(i) = var%dim(i)%len  ! = 1 for UNLIM-DIM-ID
  END DO
  j = SIZE(zdimlen)

  ! ALLOCATE DATA MEMORY AND GET DATA
  SELECT CASE(var%xtype)
  CASE(NF90_BYTE)
     CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_BYTE)
     CALL nferr(substr                                            &
                ,nf90_get_var(lncid, var%id, var%dat%vb           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,7)
  CASE(NF90_SHORT)
     CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_INT)
     CALL nferr(substr                                            &
                ,nf90_get_var(lncid, var%id, var%dat%vi           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,8)     
  CASE(NF90_CHAR)
     CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_CHAR)
     CALL nferr(substr                                            &
                ,nf90_get_var(lncid, var%id, var%dat%vc           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,9)           
  CASE(NF90_INT)
     CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_INT)
     CALL nferr(substr                                            &
                ,nf90_get_var(lncid, var%id, var%dat%vi           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,10)
  CASE(NF90_FLOAT)
     CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_REAL)
     CALL nferr(substr                                            &
                ,nf90_get_var(lncid, var%id, var%dat%vr           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,11)
  CASE(NF90_DOUBLE)
     CALL INIT_NARRAY(var%dat, j, zdimlen, VTYPE_DOUBLE)
     CALL nferr(substr                                            &
                ,nf90_get_var(lncid, var%id, var%dat%vd           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,12)
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, &
          'UNKNOWN VARIABLE TYPE ',var%xtype,' ', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'OF VARIABLE '''//TRIM(var%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(file)//''' !')
  END SELECT

  ! CLOSE FILE
  IF (.NOT.PRESENT(ncid)) THEN
     CALL nferr(substr       &
          ,nf90_close(lncid) &
          ,13)
  END IF

  ! CLEAN UP
  DEALLOCATE(start, count, stride, map, zdimlen, STAT=status)
  CALL ERRMSG(substr,status,15)

END SUBROUTINE IMPORT_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_NCVAR(var, file, ncid)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar),     INTENT(INOUT)           :: var     ! variable
  CHARACTER(LEN=*), INTENT(IN),   OPTIONAL  :: file    ! netCDF filename
  INTEGER,          INTENT(IN),   OPTIONAL  :: ncid    ! netCDF file ID

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'EXPORT_NCVAR'
  LOGICAL                            :: lex      ! file exists ?
  INTEGER                            :: lncid    ! local netCDF file-ID
  INTEGER                            :: i        ! counter
  INTEGER, DIMENSION(:), ALLOCATABLE :: start, count, stride, map
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
     CALL nferr(substr       &
          ,nf90_redef(lncid) &
          ,1)
  ELSE
     ! CHECK IF FILE EXISTS
     INQUIRE(file=TRIM(file), exist=lex)
     IF (lex) THEN   ! FILE EXISTS
        CALL RGMSG(substr, RGMLI, &
             'FILE '''//TRIM(file)//''' EXISTS ! OPENING ...')
        CALL nferr(substr                            &
             ,nf90_open(TRIM(file),NF90_WRITE,lncid) & 
             ,2)
        CALL nferr(substr       &
             ,nf90_redef(lncid) &
             ,3)
     ELSE            ! NEW FILE
        CALL RGMSG(substr, RGMLI, &
             'FILE '''//TRIM(file)//''' DOES NOT EXIST ! NEW FILE ...')
        CALL nferr(substr                                &
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
     CALL nferr(substr                &
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
     CALL nferr(substr                  &
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
        CALL nferr(substr                           &
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
                   ' ... HAS DIFFERENT DIMENSION LENGHT')
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
        CALL nferr(substr                       &
             ,nf90_def_var(lncid,TRIM(var%name) &
             ,var%xtype, zdids, var%id)         &
             ,7)
        DEALLOCATE(zdids, STAT=status)
        CALL ERRMSG(substr,status,4)
     ELSE
        CALL nferr(substr                       &
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
  CALL nferr(substr,nf90_enddef(lncid),9)

  CALL RGMSG(substr, RGMLI, &
       'EXPORTING DATA OF VARIABLE '''//TRIM(var%name)//''' ...')

  ! ALLOCATE VECTORS
  dim = 1
  IF (var%ndims > 0) THEN        ! DIM > 0
     ALLOCATE(start(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,5)
     ALLOCATE(count(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,6)
     ALLOCATE(stride(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,7)
     ALLOCATE(map(var%ndims), STAT=status)
     CALL ERRMSG(substr,status,8)
     DO i=1, var%ndims
        IF ((var%dim(i)%id == var%uid).AND.(var%uid /= NULL_DIMID).AND. &
             (ustep > 0)) THEN  ! UNLIMITED DIMENSION
           start(i) = ustep
           count(i) = 1
           stride(i) = 1
        ELSE
           ! LIMITED DIMENSION OR (USTEP <= 0) (-> COMPLETE VARIABLE)
           start(i) = 1
           count(i) = var%dim(i)%len
           stride(i) = 1
        END IF
        dim = dim * count(i)
        IF (i == 1) THEN
           map(i) = 1
        ELSE
           map(i) = map(i-1)*count(i-1)
        END IF
     END DO
  ELSE
     ALLOCATE(start(1), STAT=status)
     CALL ERRMSG(substr,status,9)
     ALLOCATE(count(1), STAT=status)
     CALL ERRMSG(substr,status,10)
     ALLOCATE(stride(1), STAT=status)
     CALL ERRMSG(substr,status,11)
     ALLOCATE(map(1), STAT=status)
     CALL ERRMSG(substr,status,12)
     start(1) = 1
     count(1) = 1
     stride(1) = 1
     map(1) = 1
  END IF

  ! WRITE DATA
  SELECT CASE (var%xtype)
  CASE(NF90_BYTE)
     CALL nferr(substr                                            &
                ,nf90_put_var(lncid, var%id, var%dat%vb           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,10)
  CASE(NF90_SHORT)
     CALL nferr(substr                                            &
                ,nf90_put_var(lncid, var%id, var%dat%vi           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,11)     
  CASE(NF90_CHAR)
     CALL nferr(substr                                            &
                ,nf90_put_var(lncid, var%id, var%dat%vc           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,12)           
  CASE(NF90_INT)
     CALL nferr(substr                                            &
                ,nf90_put_var(lncid, var%id, var%dat%vi           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,13)
  CASE(NF90_FLOAT)
     CALL nferr(substr                                            &
                ,nf90_put_var(lncid, var%id, var%dat%vr           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
                         ,14)
  CASE(NF90_DOUBLE)
     CALL nferr(substr                                            &
                ,nf90_put_var(lncid, var%id, var%dat%vd           &
                         ,start=start, count=count, stride=stride &
                         ,map = map)                              &
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
     CALL nferr(substr,nf90_close(lncid),16)
  END IF

  ! CLEAN UP
  DEALLOCATE(start, count, stride, map, STAT=status)
  CALL ERRMSG(substr,status,13)

  CALL RGMSG(substr, RGMLIC, &
       '... END EXPORTING VRIABLE '''//TRIM(var%name)//''' !')

END SUBROUTINE EXPORT_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE ADD_NCATT(var, name, replace, vs, vr, vd, vi, vb)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar),               INTENT(INOUT)          :: var
  CHARACTER(LEN=*),           INTENT(IN)             :: name
  LOGICAL,                    INTENT(IN),   OPTIONAL :: replace
  CHARACTER(LEN=*),           INTENT(IN),   OPTIONAL :: vs
  REAL (SP), DIMENSION(:),    INTENT(IN),   OPTIONAL :: vr
  REAL (DP), DIMENSION(:),    INTENT(IN),   OPTIONAL :: vd
  INTEGER (I8), DIMENSION(:), INTENT(IN),   OPTIONAL :: vi
  INTEGER (I4), DIMENSION(:), INTENT(IN),   OPTIONAL :: vb

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'ADD_NCATT'
  INTEGER                                 :: i
  INTEGER                                 :: natts  ! number of attributes
  INTEGER                                 :: num    ! number of new attribute
  TYPE (ncatt), DIMENSION(:), ALLOCATABLE :: att  ! for saving
  INTEGER                                 :: status
  LOGICAL                                 :: lrpl   ! replace ?
  INTEGER                                 :: narg
  INTEGER                                 :: vtype
  INTEGER                                 :: nf90type

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
          'EXISTS ALREADY WITH NUMBER ',var%att(num)%num,' ')
     CALL RGMSG(substr, RGMLWC,&
          'AND TYPE ',var%att(num)%xtype,' ')
     CALL RGMSG(substr, RGMLWC, 'AND SHOULD BE REPLACED !')
     RETURN
  END IF

  IF ((TRIM(var%att(num)%name) == TRIM(name)).AND.(lrpl)) THEN
     CALL RGMSG(substr, RGMLI, 'REPLACING ATTRIBUTE '''//TRIM(name)//'''')
     CALL RGMSG(substr, RGMLIC,'OF VARIABLE '''//TRIM(var%name)//''' !')
  END IF

  ! SET ATTRIBUTE
  var%att(num)%name  = TRIM(name)
  var%att(num)%num   = num
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
  TYPE (ncvar), DIMENSION(:)  , POINTER              :: var   ! variables
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
     CALL nferr(substr                              &
          ,nf90_open(TRIM(file),NF90_NOWRITE,lncid) & 
          ,1)
  END IF

  ! GET VARIABLES
  CALL nferr(substr                            &
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
     CALL nferr(substr       &
          ,nf90_close(lncid) & 
          ,3)
  END IF

END SUBROUTINE SCAN_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RENAME_NCVAR(var, newname)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar),     INTENT(INOUT) :: var
  CHARACTER(LEN=*), INTENT(IN)    :: newname

  IF (TRIM(var%name) == TRIM(newname)) RETURN
  var%name = TRIM(newname)

END SUBROUTINE RENAME_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IDX2FRAC_NCVAR(vi, vf, vtype)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar), INTENT(IN)           :: vi   ! variable with index field
  TYPE (ncvar), INTENT(INOUT)        :: vf   ! variable with index fraction
  INTEGER,      INTENT(IN), OPTIONAL :: vtype

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IDX2FRAC_NCVAR'
  INTEGER                            :: i
  INTEGER                            :: qtype
  INTEGER                            :: lvtype
  INTEGER                            :: status
  INTEGER (I8)                       :: xmin, xmax  ! min. and max. index
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec
  INTEGER, DIMENSION(:), ALLOCATABLE :: svec
  INTEGER                            :: len         ! length of data field
  INTEGER, DIMENSION(:), POINTER     :: vivec       ! position vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: vfvec       ! position vector
  INTEGER                            :: dnew        ! new dim

  ! INIT
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

  ! INIT
  CALL INIT_NCVAR(vf)
  vf%name  = TRIM(vi%name)!//'_f'
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
  vf%dim(dnew)%len = xmax-xmin+1
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
!  CALL ADD_NCATT(vf, 'RG_TYPE', replace=.true., vs = 'INT' )

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
     vfvec(dnew) = GET_NARRAY_ELEMENT_I(vi%dat, vivec) - xmin + 1
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
    TYPE (narray),               INTENT(IN) :: na
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
       GET_NARRAY_ELEMENT_I = na%vi(pos)
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
  TYPE (ncvar), INTENT(IN)           :: vf    ! index fractions
  TYPE (ncvar), INTENT(INOUT)        :: vi    ! index
  INTEGER,      INTENT(IN), OPTIONAL :: vtype

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'MAXFRAC2IDX_NCVAR'
  INTEGER                            :: lvtype
  INTEGER (I8)                       :: xr(2)   ! xmin, xmax (index range)
  INTEGER                            :: i
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec
  INTEGER, DIMENSION(:), ALLOCATABLE :: svec
  INTEGER, DIMENSION(:), POINTER     :: vivec   ! position vector
  TYPE (narray)                      :: n1      ! LAST dim. extracted
  TYPE (narray)                      :: n1i     ! SORT INDEX
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
  vi%name  = TRIM(vf%name)!//'_i'
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
!  CALL ADD_NCATT(vi, 'RG_TYPE', replace=.true., vs = 'IDX' )

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
  END DO

  ! CLEAN UP
  DEALLOCATE(dimvec, svec, STAT=status)
  CALL ERRMSG(substr,status,5)

CONTAINS

  SUBROUTINE EXTRACT_NARRAY_1DIM(na, vec, n1)

    IMPLICIT NONE

    ! I/O
    TYPE (narray),              INTENT(IN)  :: na
    INTEGER,      DIMENSION(:), INTENT(IN)  :: vec ! vec. with 0 at ...
                                                   ! ... dim. to extract
    TYPE (narray),              INTENT(INOUT) :: n1

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'EXTRACT_NARRAY_1DIM'
    INTEGER :: vtype
    INTEGER :: i
    INTEGER :: nvec(na%n)
    INTEGER :: pos
    INTEGER :: zcount  ! how many zero's
    INTEGER :: zpos    ! where is the zero

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
SUBROUTINE NFERR(routine, status, no)

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

! ******************************************************************
END MODULE MESSY_NCREGRID_NETCDF
! ******************************************************************
