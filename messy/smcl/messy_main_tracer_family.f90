! ************************************************************************
MODULE messy_main_tracer_family
! ************************************************************************

  ! MODULE FOR TRACER FAMILIES  (MESSy-SMCL)
  ! 
  ! Authors:
  !    Patrick Joeckel,   MPICH, January 2004, June 2007
  !    Astrid Kerkweg,    MPICH, May 2004,     June 2007
  !    Joachim Buchholz,  MPICH, November 2004

  USE messy_main_constants_mem, ONLY: DP, TRACNAMELEN => STRLEN_MEDIUM
  USE messy_main_tracer,        ONLY: STRLEN_TRSET, NSETID, TRSET, t_trinfo &
       , I_advect, I_convect, I_vdiff, I_sedi, I_scav, I_drydep &
       , I_mix

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  PUBLIC :: DP

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr = 'tracer_family'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodver = '2.3'

  ! MAX. NUMBER OF TRACER FAMILIES
  INTEGER, PARAMETER, PUBLIC :: NMAXTFAM       = 50
  ! MAX. NUMBER OF TRACERS PER FAMILY
  INTEGER, PARAMETER, PUBLIC :: NMAXTRAC       = 100
  ! MAX LENGTH OF REAL-WEIGHT AS STRING
  INTEGER, PARAMETER         :: MAXRSTR        = 12

  ! MAX LENGTH OF INPUT TRACER-NAME STRING (name_subname:weight)
  INTEGER, PARAMETER :: MAXIOSTR = 2*TRACNAMELEN+1 + 1 + MAXRSTR

  INTEGER, PARAMETER, PUBLIC :: FTYPE_UNDEF = 0
  INTEGER, PARAMETER, PUBLIC :: FTYPE_TRAN  = 1
  INTEGER, PARAMETER, PUBLIC :: FTYPE_SCAL  = 2

  ! FOR NAMELIST ...
  TYPE IO_TFAM
     CHARACTER(LEN=STRLEN_TRSET) :: set = ''
     INTEGER                     :: type = FTYPE_UNDEF ! type of tracer family
     ! ONLY FOR FTYPE_SCAL; if FALSE, no re-scaling to sum
     LOGICAL                     :: l_rescale = .TRUE.
     CHARACTER(len=TRACNAMELEN)  :: name=''     ! name of tracer family
     CHARACTER(len=TRACNAMELEN)  :: subname=''  ! subname of tracer family
     ! tracers: 'name' OR 'name_subname' OR
     !          'name:weight' OR 'name_subname:weight'
     CHARACTER(LEN=MAXIOSTR), DIMENSION(NMAXTRAC) :: tracer=''
  END TYPE IO_TFAM
  PUBLIC :: IO_TFAM

  ! FOR INTERNAL WORKSPACE ...
  TYPE TFAM
     TYPE(IO_TFAM)                 :: IO       ! from namelist input
     INTEGER                       :: nt       ! number of tracers in family
     INTEGER                       :: fidt     ! index of family in tracer list
     INTEGER, DIMENSION(NMAXTRAC)  :: idt      ! indices of tracers in list
     TYPE(t_trinfo), DIMENSION(NMAXTRAC) :: ti ! tracer attributes
     REAL(DP), DIMENSION(NMAXTRAC) :: weight   ! weight (atoms or mass)
  END TYPE TFAM
  PUBLIC :: TFAM
  !
  TYPE L_PTR_1D_ARRAY
     LOGICAL, DIMENSION(:), POINTER :: ptr => NULL()
  END TYPE L_PTR_1D_ARRAY
  PUBLIC :: L_PTR_1D_ARRAY

  ! CTRL-NAMELIST
  LOGICAL, PUBLIC :: l_verbose = .FALSE.
  INTEGER, PUBLIC :: i_diag_jrow = 1
  INTEGER, PUBLIC :: i_diag_pe = 0
  TYPE(IO_TFAM), DIMENSION(NMAXTFAM), SAVE          :: TF

  ! WORKSPACE
  TYPE(TFAM),    DIMENSION(NMAXTFAM), PUBLIC, SAVE        :: XTF
  INTEGER,                            PUBLIC              :: NTF
  ! .TRUE. -> TRACERS ACTIVE ; .FALSE. -> FAMILIES ACTIVE
  TYPE(L_PTR_1D_ARRAY), DIMENSION(:), POINTER, SAVE, PUBLIC :: &
       IS_TRACER => NULL()

!!$  ! A GLOBAL SWITCH
!!$  LOGICAL, SAVE :: LFAMILY = .FALSE.

  PUBLIC :: tracfamily_init
  PUBLIC :: tracfamily_newtrac
  PUBLIC :: tracfamily_initmode
  PUBLIC :: tracfamily_freemem
  PUBLIC :: tracfamily_meta
  PUBLIC :: tracfamily_meta_1_t2f
  PUBLIC :: tracer_family_read_nml_ctrl  
  !PRIVATE :: parse_str

CONTAINS

  ! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_init(status)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracfamily_init'
    INTEGER :: i, j

    status = 0

    ! GET NUMBER OF ENTRIES
    NTF = 1
    DO i=1, NMAXTFAM        ! LOOP OVER FAMILIES
       ! INIT
       XTF(i)%fidt = 0
       !
       ! CHECK IF EMPTY
       IF (TRIM(TF(i)%set) == '') CYCLE
       IF (TRIM(TF(i)%name) == '') CYCLE
       IF (TF(i)%type == FTYPE_UNDEF) CYCLE
       !
       XTF(NTF)%IO%set = TF(i)%set
       XTF(NTF)%IO%type = TF(i)%type
       XTF(NTF)%IO%l_rescale = TF(i)%l_rescale
       XTF(NTF)%IO%name = TRIM(TF(i)%name)
       XTF(NTF)%IO%subname = TRIM(TF(i)%subname)
       !
       IF (XTF(NTF)%IO%subname == '') THEN
          WRITE(*,*) 'FAMILY '''//TRIM(XTF(NTF)%IO%name)//&
               &''' (SET '''//TRIM(TF(i)%set)//''', '//&
               &'TYPE ',XTF(NTF)%IO%type,') REQUESTS MEMBERS:'
       ELSE
          WRITE(*,*) 'FAMILY ''',TRIM(XTF(NTF)%IO%name)//&
               &'_'//TRIM(XTF(NTF)%IO%subname)//&
               &''' (SET '''//TRIM(TF(i)%set)//''', '//&
               &'TYPE ',XTF(NTF)%IO%type,') REQUESTS MEMBERS:'
       END IF
       
       IF (XTF(NTF)%IO%type == FTYPE_SCAL) THEN
          IF (XTF(NTF)%IO%l_rescale) THEN
             WRITE(*,*) ' ... re-scaling: ON'
          ELSE
             WRITE(*,*) ' ... re-scaling: OFF'
          END IF
       END IF
       
       XTF(NTF)%nt = 1
       DO j=1, NMAXTRAC     ! LOOP OVER TRACERS
          ! INIT
          XTF(i)%idt(j) = 0
          !
          IF (TRIM(TF(i)%tracer(j)) == '') CYCLE
          ! PARSE STRING FOR WEIGHT
          CALL parse_str(status, MAXIOSTR, TRIM(TF(i)%tracer(j)) &
               , XTF(NTF)%IO%tracer(XTF(NTF)%nt)                 &
               , XTF(NTF)%weight(XTF(NTF)%nt))
          WRITE(*,*) ' ... ',TRIM(XTF(NTF)%IO%tracer(XTF(NTF)%nt)) &
               ,XTF(NTF)%weight(XTF(NTF)%nt)
          !
          SELECT CASE(status)
          CASE(0)
             ! OK
          CASE(1)
             WRITE(*,*) 'ERROR IN READING REAL (WEIGHT)'
          CASE(2)
             WRITE(*,*) 'MORE THAN ONE '':'' IN STRING'
          CASE DEFAULT
             WRITE(*,*) 'UNKNOWN ERROR STATUS'
          END SELECT
          IF (status /= 0) RETURN
          !
          XTF(NTF)%nt =  XTF(NTF)%nt + 1
       END DO               ! LOOP OVER TRACERS
       XTF(NTF)%nt = XTF(NTF)%nt - 1
       WRITE(*,*) ' ---> ', XTF(NTF)%nt,' tracers requested'
       !
       NTF = NTF + 1
    END DO                  ! LOOP OVER FAMILIES
    NTF = NTF - 1

  END SUBROUTINE tracfamily_init
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_newtrac(status, ldiagout)

    USE messy_main_tracer,        ONLY: new_tracer, get_tracer &
                                      , get_tracer_set, t_trinfo_list &
                                      , full2base_sub &
                                      , OFF, FAMILY, STRLEN_TRSET
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG
    USE messy_main_tools,         ONLY: split_name_domain

    IMPLICIT NONE

    INTRINSIC :: NULL, ASSOCIATED, TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status
    LOGICAL, INTENT(IN)  :: ldiagout

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER        :: substr = 'tracfam_newtrac'
    LOGICAL, DIMENSION(:), ALLOCATABLE :: is_member ! member of family ?
    INTEGER                            :: i,j, n
    INTEGER                            :: ierr, idt
    INTEGER                            :: nvalid ! number of valid tracers
    INTEGER                            :: idt_att
    TYPE(t_trinfo)                     :: save_tratt, save_fam
    CHARACTER(LEN=2*TRACNAMELEN+1)     :: famfull
    CHARACTER(LEN=TRACNAMELEN)         :: basename, subname
    INTEGER                            :: ntrac
    TYPE(t_trinfo_list),       POINTER :: ztrlist => NULL()
    TYPE(t_trinfo_list),       POINTER :: til => NULL()
    CHARACTER(LEN=STRLEN_VLONG)        :: fmembers = ''
    CHARACTER(LEN=12)                  :: weightstr = ''
    CHARACTER(LEN=STRLEN_TRSET)        :: zname
    INTEGER                            :: domain
    
    ! INIT
    status = 0

!!$    LFAMILY = .TRUE.

    set_loop: DO n=1, NSETID

       ! INIT
       CALL get_tracer_set(status, TRSET(n)%name, trlist=ztrlist &
            , ntrac=ntrac)
       IF (status /= 0) RETURN

       ! NO EMPTY SETS ...
       IF (ntrac == 0) CYCLE

       ALLOCATE(is_member(ntrac))
       is_member(:) = .FALSE.

       family_loop: DO i = 1, NTF  ! LOOP OVER FAMILIES

          fmembers = ''

          IF (TRIM(XTF(i)%IO%set)  == '') CYCLE
          IF (TRIM(XTF(i)%IO%name) == '') CYCLE

          ! ONLY THIS SET ...
          CALL split_name_domain(status, TRIM(TRSET(n)%name) &
               , zname, domain)
          status = 0
          IF ( (TRIM(XTF(i)%IO%set)  /= TRIM(TRSET(n)%name)) .AND. &
               (TRIM(XTF(i)%IO%set)  /= TRIM(zname)) ) CYCLE
          
          ! CHECK NAME OF TRACER FAMILY
          CALL get_tracer(ierr, TRIM(TRSET(n)%name), TRIM(XTF(i)%IO%name) &
               , idx=idt, subname=TRIM(XTF(i)%IO%subname) )

          IF (XTF(i)%IO%subname == '') THEN
             famfull = TRIM(XTF(i)%IO%name)
          ELSE
             famfull = TRIM(XTF(i)%IO%name)//'_'//TRIM(XTF(i)%IO%subname)
          ENDIF
       
          IF (ierr == 0) THEN ! TRACER EXISTS ALREADY
             IF (ldiagout) &
                  WRITE(*,*)' '''//TRIM(famfull)//''''//&
                  &' - tracer exists already in set '''&
                  &//TRIM(TRSET(n)%name)//'''... skipping'
             CYCLE
          ELSE
             IF (ldiagout) &
                  WRITE(*,*) ' tracer family '''//TRIM(famfull)&
                  &//''' in set '''//TRIM(TRSET(n)%name)//''' contains ...'
          END IF

          ! INITIALIZE LOOP OVER TRACERS
          nvalid  = 0
          idt     = 0
          idt_att = 0

          ! LOOP OVER TRACERS
          tracer_loop: DO j=1, XTF(i)%nt

             CALL full2base_sub(status, TRIM(XTF(i)%IO%tracer(j)) &
                  , basename, subname)
             IF (status /= 0) RETURN
             CALL get_tracer(ierr, TRIM(TRSET(n)%name), basename &
                  , subname=subname, idx=idt)
             IF (ierr /= 0) THEN
                ! TRACER DOES NOT EXIST
                IF (ldiagout) &
                     WRITE(*,*) '  ... '//TRIM(XTF(i)%IO%tracer(j))//&
                     &' - tracer not found in set '''&
                     &//TRIM(TRSET(n)%name)//'''... skipping'
                CYCLE
             ELSE
                ! TRACER EXISTS
                IF (.NOT. ( (XTF(i)%IO%type == FTYPE_SCAL) &
                     .AND. (.NOT. XTF(i)%IO%l_rescale)) ) THEN
                   ! IF FTYPE_SCAL AND re-scaling is switched OFF,
                   ! a tracer might be member of more than one family
                   IF (is_member(idt)) THEN
                      ! TRACER IS ALREADY MEMBER OF A TRACER FAMILY
                      IF (ldiagout) &
                           WRITE(*,*) &
                           '  ... '//TRIM(XTF(i)%IO%tracer(j))//&
                           &' - tracer already in (other) family ... skipping'
                      CYCLE
                   END IF
                END IF
             END IF ! TRACER EXISTS

             ! TRACER IS VALID
             IF (ldiagout) &
                  WRITE(*,*) '  ... '//TRIM(XTF(i)%IO%tracer(j)) &
                  , ' ( weight = ',XTF(i)%weight(j),')'
             ! ... NOW IT BELONGS TO A FAMILY
             is_member(idt) = .TRUE.
             ! ... SAVE TRACER NUMBER
             XTF(i)%idt(j) = idt

             ! ... COUNT VALID TRACERS
             nvalid = nvalid + 1

             write(weightstr,'(f5.2)') XTF(i)%weight(j)
             IF (TRIM(subname) == '') THEN
                IF (TRIM(fmembers) == '') THEN
                   fmembers = TRIM(weightstr)//'*'//TRIM(basename)
                ELSE
                   fmembers = TRIM(fmembers)//' + '//TRIM(weightstr)&
                        &//'*'//TRIM(basename)
                ENDIF
             ELSE
                IF (TRIM(fmembers) == '') THEN
                   fmembers = TRIM(weightstr)//'*'//TRIM(basename)&
                        &//'_'//TRIM(subname)
                ELSE
                   fmembers = TRIM(fmembers)//' + '//TRIM(weightstr)&
                        &//'*'//TRIM(basename)//'_'//TRIM(subname)
                ENDIF
             ENDIF

             ! ... SET POINTER TO TRACER INFO STRUCTURE
             til => ztrlist
             DO
                IF (.NOT. ASSOCIATED(til)) EXIT
                IF (til%info%ident%idx == idt) EXIT
                til => til%next
             END DO
             ! ... SAVE ATTRIBUTES OF FIRST VALID TRACER
             IF (idt_att == 0) THEN
                idt_att = idt
                save_tratt = til%info
             END IF
             !
             ! SAVE TRACER INFO
             XTF(i)%ti(j) = til%info

          END DO tracer_loop ! LOOP OVER TRACERS

          IF (nvalid == 0) THEN  ! NO VALID TRACER
             XTF(i)%fidt = 0
             IF (ldiagout) &
                  WRITE(*,*) &
                  '    empty tracer family '''//TRIM(famfull)//&
                  &'''... skipping'
             CYCLE
          END IF

          ! CREATE NEW TRACER (FAMILY)
          CALL new_tracer(status, TRIM(TRSET(n)%name) &
               , TRIM(XTF(i)%IO%name), submodstr      &
               , subname  = TRIM(XTF(i)%IO%subname)   &
               , type = FAMILY                        &
               , longname = TRIM(fmembers)            &
               , idx=idt)
          IF (status /= 0) RETURN
       
          ! SET POINTER TO TRACER INFO STRUCTURE (FAMILY)
          til => ztrlist
          DO
             IF (.NOT. ASSOCIATED(til)) EXIT
             IF (til%info%ident%idx == idt) EXIT
             til => til%next
          END DO
       
          XTF(i)%fidt = idt
       
          ! SAVE FAMILY ATTRIBUTES
          save_fam = til%info
       
          ! COPY ATTRIBUTES FROM FIRST VALID TRACER TO FAMILY TRACER
          til%info = save_tratt
       
          ! RE-OVERWRITE SOME FAMILY ATTRIBUTES
          til%info%ident%basename    = save_fam%ident%basename
          til%info%ident%subname     = save_fam%ident%subname
          til%info%ident%fullname    = save_fam%ident%fullname
          til%info%ident%longname    = save_fam%ident%longname
          til%info%ident%submodel    = save_fam%ident%submodel
          til%info%ident%idx         = save_fam%ident%idx
          til%info%ident%type        = save_fam%ident%type
       
          IF (XTF(i)%IO%type == FTYPE_TRAN) THEN
             til%info%meta%cask_i(I_advect)  = save_fam%meta%cask_i(I_advect)
             til%info%meta%cask_i(I_convect) = save_fam%meta%cask_i(I_convect)
             til%info%meta%cask_i(I_vdiff)   = save_fam%meta%cask_i(I_vdiff)
          END IF

          IF ((XTF(i)%IO%type == FTYPE_SCAL) .AND. &
               (.NOT. XTF(i)%IO%l_rescale)) THEN
             til%info%meta%cask_i(I_advect)     = OFF 
             til%info%meta%cask_i(I_convect)    = OFF 
             til%info%meta%cask_i(I_vdiff)      = OFF 
             til%info%meta%cask_i(I_sedi)       = OFF 
             til%info%meta%cask_i(I_scav)       = OFF 
             til%info%meta%cask_i(I_drydep)     = OFF 
          END IF

       END DO family_loop ! LOOP OVER FAMILIES

       DEALLOCATE(is_member)

    END DO set_loop

  END SUBROUTINE tracfamily_newtrac
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_initmode(ll_dim)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! LOCAL
    INTEGER, INTENT(IN) :: ll_dim ! 'local loop' dimension
    INTEGER :: n, s

    IF (.NOT. ASSOCIATED(IS_TRACER)) &  ! only once
         ALLOCATE(IS_TRACER(NSETID))

    DO n=1, NSETID

       ! NO EMPTY SETS ...
       IF (TRSET(n)%ntrac == 0) CYCLE

       s = SIZE(TRSET(n)%xt, ll_dim)             ! 'local loop' dimension
       IF (.NOT. ASSOCIATED(IS_TRACER(n)%ptr)) & ! only once
            ALLOCATE(IS_TRACER(n)%ptr(s))
       IS_TRACER(n)%ptr(:) = .TRUE.

    END DO

  END SUBROUTINE tracfamily_initmode
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_meta(status, flag, callstr, setname, ldiagout)

    USE messy_main_blather, ONLY: start_message, end_message
    USE messy_main_tracer,  ONLY: t_trinfo_tp, get_tracer_set, STRLEN_TRSET
    USE messy_main_tools,   ONLY: split_name_domain

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: flag
    CHARACTER(LEN=*), INTENT(IN)  :: callstr
    CHARACTER(LEN=*), INTENT(IN)  :: setname
    LOGICAL,          INTENT(IN)  :: ldiagout

    ! LOCAL
    CHARACTER(LEN=*),  PARAMETER             :: substr = 'tracfamily_meta'
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti
    INTEGER :: ntrac
    INTEGER :: i,j
    CHARACTER(LEN=STRLEN_TRSET) :: zname
    INTEGER :: domain

    CALL get_tracer_set(status, setname, ti=ti, ntrac=ntrac)
    IF (status /= 0) RETURN

    ! NO EMPTY SETS ...
    IF (ntrac == 0) RETURN

    SELECT CASE(flag)
    CASE(1)
       !
       IF (ldiagout) &
            CALL start_message(submodstr, &
            'resetting meta information of family-members',callstr)
       
       family_loop1: DO i = 1, NTF  ! LOOP OVER FAMILIES

          ! ONLY THIS SET ...
          CALL split_name_domain(status, TRIM(setname) &
               , zname, domain)
          status = 0
          IF ( (TRIM(XTF(i)%IO%set) /= TRIM(setname)) .AND. &
               (TRIM(XTF(i)%IO%set) /= TRIM(zname)) ) CYCLE

          IF (XTF(i)%fidt == 0) CYCLE

          ! CONVERSION ONLY FOR TRANSPORT-FAMILIES
          IF (XTF(i)%IO%type /= FTYPE_TRAN) CYCLE

          tracer_loop1: DO j=1, XTF(i)%nt  ! LOOP OVER TRACERS
             IF (XTF(i)%idt(j) == 0) CYCLE
          
             ti(XTF(i)%idt(j))%tp = XTF(i)%ti(j)

          END DO tracer_loop1

       END DO family_loop1

       IF (ldiagout) &
            CALL end_message(submodstr, &
            'resetting meta information of family-members',callstr)
       !
    CASE(2)
       !
       IF (ldiagout) &
            CALL start_message(submodstr, &
            'setting meta information of family-members to fraction',callstr)
       
       family_loop2: DO i = 1, NTF  ! LOOP OVER FAMILIES

          ! ONLY THIS SET ...
          CALL split_name_domain(status, TRIM(setname) &
               , zname, domain)
          status = 0
          IF ( (TRIM(XTF(i)%IO%set) /= TRIM(setname)) .AND. &
               (TRIM(XTF(i)%IO%set) /= TRIM(zname)) ) CYCLE

          IF (XTF(i)%fidt == 0) CYCLE
          
          ! CONVERSION ONLY FOR TRANSPORT-FAMILIES
          IF (XTF(i)%IO%type /= FTYPE_TRAN) CYCLE

          tracer_loop2: DO j=1, XTF(i)%nt  ! LOOP OVER TRACERS
             IF (XTF(i)%idt(j) == 0) CYCLE
             
             CALL tracfamily_meta_1_t2f(ti(XTF(i)%idt(j))%tp, &
                  XTF(i)%IO%name, XTF(i)%IO%subname )
             
          END DO tracer_loop2

       END DO family_loop2

       IF (ldiagout) &
            CALL end_message(submodstr, &
            'setting meta information of family-members to fraction',callstr)
       !
    CASE DEFAULT
       !
       status = 2010 ! UNKNOWN CONVERSION FLAG
       !
    END SELECT

  END SUBROUTINE tracfamily_meta
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_freemem

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER :: n
    
    DO n=1, NSETID
       IF (ASSOCIATED(IS_TRACER(n)%ptr)) THEN
          DEALLOCATE(IS_TRACER(n)%ptr)
          NULLIFY(IS_TRACER(n)%ptr)
       END IF
    END DO
    DEALLOCATE(IS_TRACER) ; NULLIFY(IS_TRACER)

  END SUBROUTINE tracfamily_freemem
  ! ----------------------------------------------------------------------

  ! **********************************************************************
  ! PRIVATE ROUTINES
  ! **********************************************************************

  ! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_meta_1_t2f(ti, name, subname)

    USE messy_main_tracer, ONLY: OFF
    
    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    TYPE(t_trinfo),   INTENT(INOUT) :: ti
    CHARACTER(LEN=*), INTENT(IN)    :: name, subname  ! tracer family
    
    ! IDENTIFICATION
    ! ti%ident%basename
    ! ti%ident%subname
    ! ti%ident%fullname
    ! ti%ident%longname
    IF (TRIM(subname) == '') THEN
       ti%ident%unit      = 'frac. of '//TRIM(name)
    ELSE
       ti%ident%unit      = 'frac. of '//TRIM(name)//'_'//TRIM(subname)
    END IF
    ! ti%ident%submodel
    ! ti%ident%idx
    ! ti%ident%medium
    ! ti%ident%quantity
    ! ti%ident%type

    ! PROCESSES
    ti%meta%cask_i(I_advect)     = OFF
    ti%meta%cask_i(I_convect)    = OFF
    ti%meta%cask_i(I_vdiff)      = OFF
    ti%meta%cask_i(I_drydep)     = OFF
    ti%meta%cask_i(I_sedi)       = OFF
    ti%meta%cask_i(I_scav)       = OFF
    ti%meta%cask_i(I_mix)        = OFF
    ! ti%meta%cask_i(I_integrate)  = OFF
    ! ti%meta%cask_i(I_timefilter) = OFF
    ! ti%meta%cask_i(I_force_col)

  END SUBROUTINE tracfamily_meta_1_t2f
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE tracer_family_read_nml_ctrl(status, iou)

    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer, ONLY: modstr

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracer_family_read_nml_ctrl'

    NAMELIST /CTRL_FAMILY/ l_verbose, i_diag_jrow, i_diag_pe, TF

    ! LOCAL
    LOGICAL :: lex      ! file exists ?
    INTEGER :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CTRL_FAMILY', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_FAMILY, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_FAMILY', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE tracer_family_read_nml_ctrl
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE parse_str(status, strlen, str, name, weight)

    USE messy_main_tools, ONLY: strcrack

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    CHARACTER(LEN=*), INTENT(OUT)   :: name     ! tracer name
    REAL(DP),         INTENT(OUT)   :: weight   ! weight

    ! LOCAL
    CHARACTER(LEN=strlen), POINTER  :: sl1(:)
    INTEGER :: n
    INTEGER :: iostat

    NULLIFY(sl1)

    CALL strcrack(str, ':', sl1, n)
#if defined (__G95__)
    IF (n>1) READ(sl1(2),*,IOSTAT=iostat) weight
#endif

    SELECT CASE(n)
       CASE(1)
          ! NO ':' -> DEFAULT WEIGHT
          name = TRIM(sl1(1))
          weight = 1.0_DP ! DEFAULT
          status = 0
       CASE(2)
          ! name:weight
          name = TRIM(sl1(1))
          READ(sl1(2),*,IOSTAT=iostat) weight
          IF (iostat /= 0) THEN
             status = 1  ! ERROR IN READING REAL
          ELSE
             status = 0   
          END IF
       CASE DEFAULT
          ! ERROR: MORE THEN ONE ':' IN STRING
          status = 2
    END SELECT

    ! CLEAN UP
    IF (ASSOCIATED(sl1)) THEN
       DEALLOCATE(sl1) ; NULLIFY(sl1)
    END IF

  END SUBROUTINE parse_str
  ! ----------------------------------------------------------------------

! ************************************************************************
END MODULE messy_main_tracer_family
! ************************************************************************
