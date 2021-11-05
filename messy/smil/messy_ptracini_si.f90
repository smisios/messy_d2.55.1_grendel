#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_ptracini_si
! **********************************************************************
  ! Prognostic TRACer INItialisation
  !
  ! INTERFACE FOR MESSy/SMIL
  !
  ! Author: Christiane Hofmann, UniMz, November 2011 - June 2013
  !         Astrid Kerkweg, UniMZ, November 2011

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM,dp
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY, PTR_2D_ARRAY
  USE messy_main_data_bi,       ONLY: l2tls
  USE messy_ptracini

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: NULL

  ! NAMELIST LIMITS
  INTEGER, PARAMETER :: NMAXTRAC  =  5  ! MAX. NUMBER OF TRACERS PER SET
  INTEGER, PARAMETER :: NMAXCRIT  = 10  ! MAX. NUMBER OF CRITs PER SET
  INTEGER, PARAMETER :: NMAXSETS  = 20  ! MAX. NUMBER OF SETS (TRINI)

  ! USER INTERFACE
  ! ... define tracers, criterions,logicals ...
  TYPE COMP
     REAL(dp)           :: const = -9999._dp ! const to which field is comp to 
     TYPE(t_chaobj_cpl) :: compfld           ! fld to which field is comp to
  END TYPE COMP

  TYPE CRITDEF
     TYPE(t_chaobj_cpl) :: field             ! field which has to meet crit
     CHARACTER(LEN=2)   :: operator = ''     ! operator
     TYPE(COMP)         :: comp              ! const. value or field 
  END type CRITDEF

  TYPE STR2
     CHARACTER(LEN=STRLEN_MEDIUM) :: name       = '' ! NAME TRACER
     CHARACTER(LEN=STRLEN_MEDIUM) :: subname    = '' ! SUBNAME TRACER
  END TYPE STR2
  
  TYPE NML_INPUT 
     !time-crits
     INTEGER                            :: EVENT_START(6) = 0 !date of ini 
     TYPE(io_time_event)                :: EMIS_IOEVENT       !intervall of ini 
     LOGICAL                            :: INI1STEP           !if T, date of 
                                                              !ini is used
     LOGICAL                            :: RELAX              !relax to coarser
                                                              !instance 
     LOGICAL                            :: IFNOT              !use invers crits
     LOGICAL                            :: DOMAIN             !ini only in dom
     TYPE(STR2), DIMENSION(NMAXTRAC)    :: TRACER  
     TYPE(CRITDEF), DIMENSION(NMAXCRIT) :: CRIT
  END type NML_INPUT
  
  TYPE(NML_INPUT), DIMENSION(NMAXSETS)  :: TRINI

  ! NUMBER OF INITIAL SETS
  INTEGER :: NUMSETS = 0


  ! WORKSPACE
  ! ... representation independent
  TYPE SETDEF
     INTEGER            :: EVENT_START(6) = 0 
     TYPE(io_time_event):: EMIS_IOEVENT = io_time_event(1, 'steps','first',0)
     TYPE(time_event)   :: EMIS_EVENT
     LOGICAL            :: lemis    = .FALSE. !is set T, if time-crit fullfilled
     LOGICAL            :: lfirst   = .TRUE.  !forces ini 2nd-step if leapfrog
     LOGICAL            :: linistart = .FALSE. !special treatment: if 1st step 
     LOGICAL            :: lini2nd  = .TRUE.   !and leapfrog (dt(1)=0.5*dt(n))
     LOGICAL            :: INI1STEP = .FALSE. 
     LOGICAL            :: RELAX    = .FALSE.
     LOGICAL            :: lrelax   = .FALSE.
     TYPE(time_event)   :: RELAX_EVENT
     LOGICAL            :: IFNOT    = .FALSE.
     LOGICAL            :: DOMAIN   = .FALSE.
     INTEGER            :: NUMTRAC  = 0
     TYPE(STR2)                  , DIMENSION(:), ALLOCATABLE :: TRACER 
     INTEGER                     , DIMENSION(:), ALLOCATABLE :: idt 
     INTEGER            :: NUMCRIT  = 0
     TYPE(CRITDEF)               , DIMENSION(:), ALLOCATABLE :: CRIT    
     TYPE(PTR_3D_ARRAY)          , DIMENSION(:), POINTER     :: CRT
     TYPE(PTR_2D_ARRAY)          , DIMENSION(:), POINTER     :: CRT_2D
     TYPE(PTR_3D_ARRAY)          , DIMENSION(:), POINTER     :: CRT_comp
     TYPE(PTR_2D_ARRAY)          , DIMENSION(:), POINTER     :: CRT_comp2D
  END TYPE SETDEF
     
  TYPE(SETDEF), DIMENSION(:), ALLOCATABLE :: SET
  !calculate additional fields 
  LOGICAL               :: lwind   = .FALSE.
  LOGICAL               :: lthetae = .FALSE.

  LOGICAL                          :: ldomain = .FALSE.!is set true if DOMAINFLD
                                                       !is defined
  CHARACTER(LEN=STRLEN_MEDIUM)     :: DOMAINFLD(2) = '' 
  REAL(dp), DIMENSION(:,:),POINTER :: input_DOMAINFLD 
  
  ! NEW CHANNEL OBJECTS: (DIAGNOSED FIELDS)
  REAL(DP), DIMENSION(:,:,:), POINTER :: wind        ! windspeed
  REAL(DP), DIMENSION(:,:,:), POINTER :: thetae      ! equiv. potential temp.

  PUBLIC :: ptracini_initialize
  PUBLIC :: ptracini_init_memory
  PUBLIC :: ptracini_init_coupling
  PUBLIC :: ptracini_global_start
  PUBLIC :: ptracini_vdiff
  PUBLIC :: ptracini_physc

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE ptracini_initialize

    ! ptracini MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    ! BML/MESSy
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_timer_bi,      ONLY: p_bcast_event
    ! SMCL
    USE messy_main_tools,         ONLY: find_next_free_unit
   
    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ptracini_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i1, i2, i3, idx
    LOGICAL                     :: lfield(NMAXSETS)
    

    CALL start_message_bi(modstr,'INITIALISATION ',substr)

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL ptracini_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    !DETERMINE AND ALLOCATE SETS ON IO PE  
    IF (p_parallel_io) THEN
       IF ((TRIM(DOMAINFLD(1)) /= '') .OR. & 
            (TRIM(DOMAINFLD(2)) /= ''))    &
            ldomain=.TRUE.
       ! DETERMINE NUMBER OF SETS
       NUMSETS = 0 
       sets1: DO i1=1,NMAXSETS
          IF ((TRINI(i1)%INI1STEP) .AND. (SUM(TRINI(i1)%EVENT_START) == 0 )) &
               CALL error_bi(substr, 'EVENT START DATE NOT GIVEN IN NAMELIST')
          IF ((TRINI(i1)%DOMAIN) .AND. (.NOT. ldomain)) &
               CALL error_bi('NO DOMAINFLD DEFINED', substr)
          lfield(i1) = .FALSE.
          trac1: DO i2= 1, NMAXTRAC
             IF (TRIM(TRINI(i1)%TRACER(i2)%name) =='') CYCLE
             IF (TRIM(TRINI(i1)%TRACER(i2)%subname) == '') THEN
                WRITE(*,*) '####','TRACER =' &
                     , TRIM(TRINI(i1)%TRACER(i2)%name),'####'
             ELSE
                WRITE(*,*) '####','TRACER =' &
                     , TRIM(TRINI(i1)%TRACER(i2)%name)//'_'// &
                     (TRINI(i1)%TRACER(i2)%subname),'####'
             END IF
             ! TRACER exists
             field1: DO i3 = 1,NMAXCRIT
                IF (TRIM(TRINI(i1)%CRIT(i3)%field%CHA) == '' .OR. &
                     TRIM(TRINI(i1)%CRIT(i3)%field%OBJ) == '') THEN
                   CYCLE
                ELSE
                   ! FIELD exists
                   lfield(i1) = .TRUE.
                   IF (TRIM(TRINI(i1)%CRIT(i3)%operator) /= '>'   .AND. &
                        TRIM(TRINI(i1)%CRIT(i3)%operator) /= '>=' .AND. &
                        TRIM(TRINI(i1)%CRIT(i3)%operator) /= '<'  .AND. &
                        TRIM(TRINI(i1)%CRIT(i3)%operator) /= '<=' .AND. &
                        TRIM(TRINI(i1)%CRIT(i3)%operator) /= 'IF' ) & 
                        CALL error_bi('INCORRECT OPERATOR', substr)
                   !no comp switched
                   IF (TRINI(i1)%CRIT(i3)%comp%const < -9990._dp .AND. &  
                        (TRIM(TRINI(i1)%CRIT(i3)%comp%compfld%CHA) == '' .OR. &
                        TRIM(TRINI(i1)%CRIT(i3)%comp%compfld%OBJ) == ''))     & 
                        CALL error_bi('NO CONST OR COMPFLD DEFINED', substr)
                   !const and fld switched 
                   IF (TRINI(i1)%CRIT(i3)%comp%const > -9990._dp .AND. &  
                        (TRIM(TRINI(i1)%CRIT(i3)%comp%compfld%CHA) /= '' .OR. &
                        TRIM(TRINI(i1)%CRIT(i3)%comp%compfld%OBJ) /= ''))     & 
                        CALL error_bi('CONST AND COMPFLD DEFINED', substr)
                   CYCLE
                ENDIF
             END DO field1
             IF (lfield(i1)) NUMSETS = NUMSETS + 1
             IF (lfield(i1)) EXIT
          END DO trac1
       END DO sets1
       
       ! ALLOCATE SET
       ALLOCATE(SET(NUMSETS))

       ! DETERMINE TRACER AND FIELDS FOR EACH SET 
       NUMSETS = 0
       sets2: DO i1=1,NMAXSETS
          IF (.NOT. lfield(i1)) CYCLE
          NUMSETS = NUMSETS + 1
          SET(NUMSETS)%EVENT_START  = TRINI(i1)%EVENT_START
          SET(NUMSETS)%EMIS_IOEVENT = TRINI(i1)%EMIS_IOEVENT
          SET(NUMSETS)%INI1STEP     = TRINI(i1)%INI1STEP 
          SET(NUMSETS)%RELAX        = TRINI(i1)%RELAX
          SET(NUMSETS)%IFNOT        = TRINI(i1)%IFNOT
          SET(NUMSETS)%DOMAIN       = TRINI(i1)%DOMAIN
          SET(NUMSETS)%NUMTRAC      = 0
          ! DETERMINE NUMBER OF TRACER
          trac2a: DO i2= 1, NMAXTRAC
             IF (TRIM(TRINI(i1)%TRACER(i2)%name) =='') THEN 
                CYCLE
             ELSE
             ! TRACER exists
                SET(NUMSETS)%NUMTRAC = SET(NUMSETS)%NUMTRAC + 1
             END IF
          END DO trac2a
          
          ALLOCATE(SET(NUMSETS)%TRACER(SET(NUMSETS)%NUMTRAC))
          ALLOCATE(SET(NUMSETS)%idt(SET(NUMSETS)%NUMTRAC))
          ! SORT TRACER BY DEDICATING INDICES
          idx = 0
          trac2b: DO i2= 1, NMAXTRAC
             IF (TRIM(TRINI(i1)%TRACER(i2)%name) =='') THEN
                CYCLE
             ELSE
                ! TRACER exists
                idx = idx + 1
                SET(NUMSETS)%TRACER(idx)%name = TRINI(i1)%TRACER(i2)%name
                SET(NUMSETS)%TRACER(idx)%subname=TRINI(i1)%TRACER(i2)%subname
             END IF
          END DO trac2b
          ! DETERMINE NUMBER OF FIELDS
          SET(NUMSETS)%NUMCRIT = 0
          field2a: DO i3 = 1,NMAXCRIT
             IF (TRIM(TRINI(i1)%CRIT(i3)%field%CHA) == '' .OR. &
                  TRIM(TRINI(i1)%CRIT(i3)%field%OBJ) == '') THEN
                CYCLE
             ELSE
                ! FIELD exists
                SET(NUMSETS)%NUMCRIT = SET(NUMSETS)%NUMCRIT + 1
             ENDIF
          END DO field2a

          ALLOCATE(SET(NUMSETS)%CRIT(SET(NUMSETS)%NUMCRIT))
          ALLOCATE(SET(NUMSETS)%CRT(SET(NUMSETS)%NUMCRIT))
          ALLOCATE(SET(NUMSETS)%CRT_2D(SET(NUMSETS)%NUMCRIT))
          ALLOCATE(SET(NUMSETS)%CRT_comp(SET(NUMSETS)%NUMCRIT))
          ALLOCATE(SET(NUMSETS)%CRT_comp2D(SET(NUMSETS)%NUMCRIT))
                  
          ! SORT FIELDS BY DEDICATING INDICES
          idx = 0
          field2b: DO i3 = 1,NMAXCRIT
             IF (TRIM(TRINI(i1)%CRIT(i3)%field%CHA) == '' .OR. &
                  TRIM(TRINI(i1)%CRIT(i3)%field%OBJ) == '') THEN
                CYCLE
             ELSE
                idx = idx + 1 
                ! FIELD exists
                SET(NUMSETS)%CRIT(idx)%field%CHA = TRINI(i1)%CRIT(i3)%field%CHA
                SET(NUMSETS)%CRIT(idx)%field%OBJ = TRINI(i1)%CRIT(i3)%field%OBJ
                SET(NUMSETS)%CRIT(idx)%operator  = TRINI(i1)%CRIT(i3)%operator
                SET(NUMSETS)%CRIT(idx)%comp%const = TRINI(i1)%CRIT(i3)%comp%const
                SET(NUMSETS)%CRIT(idx)%comp%compfld%CHA = &
                     TRINI(i1)%CRIT(i3)%comp%compfld%CHA
                SET(NUMSETS)%CRIT(idx)%comp%compfld%OBJ = &
                     TRINI(i1)%CRIT(i3)%comp%compfld%OBJ 
                
                NULLIFY(SET(NUMSETS)%CRT(idx)%ptr)
                NULLIFY(SET(NUMSETS)%CRT_2D(idx)%ptr)
                NULLIFY(SET(NUMSETS)%CRT_comp(idx)%ptr)
                NULLIFY(SET(NUMSETS)%CRT_comp2D(idx)%ptr)
             ENDIF
          END DO field2b
       END DO sets2
    END IF 
    
    ! BROADCAST SETS, TRACERS AND FIELDS AND ALLOCATE on NON IO PEs
    CALL p_bcast(NUMSETS,  p_io)
    IF (.NOT. p_parallel_io)  ALLOCATE(SET(NUMSETS))
    
    DO i1=1,NUMSETS
       IF (p_parallel_io) THEN
          WRITE(*,*) '#################################'
          WRITE(*,*) 'SET=',      i1
       ENDIF
       CALL p_bcast(SET(i1)%EVENT_START, p_io)
       CALL p_bcast_event(SET(i1)%EMIS_IOEVENT, p_io)
       CALL p_bcast(SET(i1)%INI1STEP,    p_io)
       CALL p_bcast(SET(i1)%RELAX,       p_io)
       CALL p_bcast(SET(i1)%IFNOT,       p_io)
       CALL p_bcast(SET(i1)%DOMAIN,      p_io)
       CALL p_bcast(SET(i1)%NUMTRAC,     p_io)
       CALL p_bcast(SET(i1)%NUMCRIT,     p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(SET(i1)%TRACER(SET(i1)%NUMTRAC))
          ALLOCATE(SET(i1)%idt(SET(i1)%NUMTRAC))
          ALLOCATE(SET(i1)%CRIT(SET(i1)%NUMCRIT))
          ALLOCATE(SET(i1)%CRT(SET(i1)%NUMCRIT))
          ALLOCATE(SET(i1)%CRT_2D(SET(i1)%NUMCRIT)) 
          ALLOCATE(SET(i1)%CRT_comp(SET(i1)%NUMCRIT))
          ALLOCATE(SET(i1)%CRT_comp2D(SET(i1)%NUMCRIT))
       ENDIF

       DO i2=1,SET(i1)%NUMTRAC
          CALL p_bcast(SET(i1)%TRACER(i2)%name,    p_io)
          CALL p_bcast(SET(i1)%TRACER(i2)%subname, p_io)
       END DO
       DO i3=1,SET(i1)%NUMCRIT
          
          CALL p_bcast(SET(i1)%CRIT(i3)%field%CHA,    p_io)
          CALL p_bcast(SET(i1)%CRIT(i3)%field%OBJ,    p_io)
          CALL p_bcast(SET(i1)%CRIT(i3)%operator,     p_io)
          CALL p_bcast(SET(i1)%CRIT(i3)%comp%const,   p_io)
          CALL p_bcast(SET(i1)%CRIT(i3)%comp%compfld%CHA, p_io)
          CALL p_bcast(SET(i1)%CRIT(i3)%comp%compfld%OBJ, p_io)
          IF (p_parallel_io) THEN
             WRITE(*,*) '-----------------------------------'
             WRITE(*,*) 'CHA=',   SET(i1)%CRIT(i3)%field%CHA
             WRITE(*,*) 'OBJ=',   SET(i1)%CRIT(i3)%field%OBJ
             WRITE(*,*) 'ifnot=', SET(i1)%IFNOT
             WRITE(*,*) 'domain=',SET(i1)%DOMAIN
             WRITE(*,*) 'operator=',  SET(i1)%CRIT(i3)%operator
             IF ((SET(i1)%CRIT(i3)%comp%const) < -9990._dp) THEN
                WRITE(*,*) 'compfld_CHA=', SET(i1)%CRIT(i3)%comp%compfld%CHA
                WRITE(*,*) 'compfld_OBJ=', SET(i1)%CRIT(i3)%comp%compfld%OBJ
             ELSE
                WRITE(*,*) 'const=', SET(i1)%CRIT(i3)%comp%const
             ENDIF
          END IF
       ENDDO
       
    ENDDO

    CALL p_bcast(lwind,       p_io)
    CALL p_bcast(lthetae,     p_io)
    CALL p_bcast(ldomain,     p_io)
    CALL p_bcast(DOMAINFLD(1),p_io)
    CALL p_bcast(DOMAINFLD(2),p_io)
  
    CALL end_message_bi(modstr,'INITIALISATION ',substr)

  END SUBROUTINE ptracini_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ptracini_init_memory
    ! PTRACINI MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! define ptracini specific channel(s) and allocate memory
    !
    ! Author: Patrick Joeckel, MPICH, Sep 2003
    
    ! ECHAM5/MESSy
    USE messy_main_timer_bi,         ONLY: timer_event_init
    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'ptracini_init_memory'
    INTEGER                     :: status ! error status
    INTEGER                     :: i1

     ! INITIALIZE EVENT FOR READING OF EXTERNAL DATA 
    CALL start_message_bi(modstr, 'SET EVENTS', substr) 
    
    setloop:DO i1= 1, NUMSETS
       IF (.not.SET(i1)%INI1STEP)THEN
          CALL timer_event_init(SET(i1)%EMIS_EVENT,SET(i1)%EMIS_IOEVENT,&
               'EMIS_EVENT', 'present')
#ifdef COSMO          
          !set relax for cont. ini. in cosmo not before 2nd timestep because 
          !tracer fld in coarser inst. would be zero at this point at 1st step 
          CALL timer_event_init(SET(i1)%RELAX_EVENT,SET(i1)%EMIS_IOEVENT,&
               'RELAX_EVENT', 'previous')
#endif          
       ELSE
          WRITE(*,*) 'NO EVENTS TO SET FOR SET=', i1
       END IF
    END DO setloop
    CALL end_message_bi(modstr, 'SET EVENTS', substr)
   
     
    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)
    ! define new channel
    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    ! PTRACINI DIAGNOSTIC OUTPUT
    IF (lwind) THEN
       CALL new_channel_object(status, modstr, &
            'wind', p3=wind)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'wind', 'long_name', c='horizontal windspeed' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'wind', 'units', c='m/s' )
       CALL channel_halt(substr, status)
    END IF
    IF (lthetae) THEN
       CALL new_channel_object(status, modstr, &
            'thetae', p3=thetae)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'thetae', 'long_name', c='equivalent potential temperature' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , 'thetae', 'units', c='K' )
       CALL channel_halt(substr, status)
    END IF
    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)
   
  END SUBROUTINE ptracini_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ptracini_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR 
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_tracer,           ONLY: get_tracer
    USE messy_main_channel,          ONLY: get_channel_object 
    
  
    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    INTEGER                          :: status
    CHARACTER(LEN=*), PARAMETER      :: substr = 'ptracini_init_coupling'
    INTEGER                          :: i1, i2, i3
    
    CALL start_message_bi(modstr, 'COUPLING INITIALIZATION', substr)
    status = 0
   
    IF (ldomain) THEN  !CPL DOMAINFLD
       CALL get_channel_object(status, cname=domainfld(1), oname=domainfld(2) &
            , p2=input_DOMAINFLD)
       IF (status /= 0) THEN
          WRITE(*,*) 'error in get_channel_object'
       END IF
       CALL channel_halt(substr, status)
    END IF
    
    DO i1=1,NUMSETS
       ! COUPLE TRACER
       DO i2=1,SET(i1)%NUMTRAC
          ! get tracer
          CALL get_tracer(status, GPTRSTR, TRIM(SET(i1)%TRACER(i2)%name) &
               , subname=TRIM(SET(i1)%TRACER(i2)%subname)                &
               , idx = SET(i1)%idt(i2))                                  
          IF (status /= 0) THEN
             WRITE(*,*) 'TRACER ', TRIM(SET(i1)%TRACER(i2)%name),'_' &
                  ,TRIM(SET(i1)%TRACER(i2)%subname), 'NOT FOUND'
             CALL tracer_halt(substr,status)
          ENDIF
       END DO
       
       ! COUPLE FIELDS
       DO i3=1,SET(i1)%NUMCRIT
          CALL get_channel_object(status, SET(i1)%CRIT(i3)%field%CHA &
               , SET(i1)%CRIT(i3)%field%OBJ, p2=SET(i1)%CRT_2D(i3)%ptr)
           IF(status /= 0) THEN
             IF (status == 2019) THEN
                CALL get_channel_object(status, SET(i1)%CRIT(i3)%field%CHA &
                     , SET(i1)%CRIT(i3)%field%OBJ, p3=SET(i1)%CRT(i3)%ptr)
              ENDIF
              CALL channel_halt(substr, status)
           ENDIF
       ! IF COMP = FLD, COUPLE COMPFLD
          IF ((SET(i1)%CRIT(i3)%comp%const) < -9990._dp) THEN
             CALL get_channel_object(status, SET(i1)%CRIT(i3)%comp%compfld%CHA &
                  , SET(i1)%CRIT(i3)%comp%compfld%OBJ &
                  , p2=SET(i1)%CRT_comp2D(i3)%ptr)
             IF(status /= 0) THEN
                IF (status == 2019) THEN
                   CALL get_channel_object(status &
                        , SET(i1)%CRIT(i3)%comp%compfld%CHA &
                        , SET(i1)%CRIT(i3)%comp%compfld%OBJ &
                        , p3=SET(i1)%CRT_comp(i3)%ptr)
                END IF
                CALL channel_halt(substr, status)
             END IF
          END IF
          IF (ASSOCIATED(SET(i1)%CRT_2D(i3)%ptr) .AND. &
               ASSOCIATED(SET(i1)%CRT_comp(i3)%ptr))   & 
               !case (field=2D and fldcomp=3D) should be impossible   
               CALL error_bi('field3d comp to compfld2d', substr)
       END DO
    END DO
  
    CALL end_message_bi(modstr, 'COUPLING INITIALIZATION', substr)

   END SUBROUTINE ptracini_init_coupling
   ! -------------------------------------------------------------------

   ! -------------------------------------------------------------------
   SUBROUTINE ptracini_global_start

     USE messy_main_timer_bi,      ONLY: event_state
     USE messy_main_timer,         ONLY: current_date,gregor2julian,      &
                                         YEAR, MONTH, DAY, HOUR, MINUTE,  &
                                         SECOND,lstart
#ifdef COSMO
     USE messy_main_timer,           ONLY: previous_date, time_step_len
     USE messy_main_tracer,          ONLY: set_tracer, ON, I_RELAX
     USE messy_main_tracer_tools_bi, ONLY: tracer_halt 
     USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
#endif     
     IMPLICIT NONE

     ! LOCAL
     CHARACTER(LEN=*), PARAMETER :: substr = 'ptracini_global_start'
     INTEGER                     :: i1
     REAL(dp)                    :: NOW, NOW_START
     REAL(DP), PARAMETER         :: tiny =  1e-8_dp
     INTEGER                     :: YEAR_START, MONTH_START, DAY_START
     INTEGER                     :: HOUR_START,MINUTE_START, SECOND_START
#ifdef COSMO
     INTEGER                     :: status
     INTEGER                     :: i2
     INTEGER                     :: SECOND_RELAX
     REAL(dp)                    :: NOW_RELAX
#endif
     
     CALL start_message_bi(modstr, 'ptracini global_start', substr)
     setloop:DO i1= 1, NUMSETS
        !set lemis 
        IF (SET(i1)%INI1STEP) THEN
           YEAR_START   = SET(i1)%EVENT_START(1)
           MONTH_START  = SET(i1)%EVENT_START(2)
           DAY_START    = SET(i1)%EVENT_START(3)
           HOUR_START   = SET(i1)%EVENT_START(4)
           MINUTE_START = SET(i1)%EVENT_START(5)
           SECOND_START = SET(i1)%EVENT_START(6)
           
           NOW = gregor2julian(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)  
           NOW_START = gregor2julian(YEAR_START, MONTH_START, DAY_START, &
                HOUR_START,MINUTE_START, SECOND_START)
           IF ( ABS(NOW-NOW_START) .le. tiny ) THEN
              SET(i1)%lemis = .TRUE.
           ELSE
              SET(i1)%lemis = .FALSE.
           END IF
           
           IF (.not. l2tls) then !leapfrog => ini also 2nd step
              IF (.not.SET(i1)%lemis) THEN 
                 IF (.not.SET(i1)%lfirst) THEN !lfirst_default=T 
                    SET(i1)%lemis  = .TRUE.    !ini also 2nd step 
                    SET(i1)%lfirst = .TRUE.
                    IF (SET(i1)%linistart) THEN
                       SET(i1)%lini2nd   = .TRUE.
                       SET(i1)%linistart = .FALSE.
                    END IF
                 END IF
              ELSE
                 IF (SET(i1)%lfirst) THEN !ini first step
                    SET(i1)%lfirst  = .FALSE.
                    SET(i1)%lini2nd = .FALSE.
                    IF (lstart) THEN
                       SET(i1)%linistart = .TRUE.
                    END IF
                 END IF
              END IF
           END IF
        ELSE  !ini1step=F 
           SET(i1)%lemis = event_state(SET(i1)%EMIS_EVENT, current_date)
        END IF
           

#ifdef COSMO
        ! set relaxation event
        IF (SET(i1)%RELAX) THEN 
           IF (SET(i1)%INI1STEP) THEN
              SECOND_RELAX = SECOND_START+INT(time_step_len)
              NOW_RELAX = gregor2julian(YEAR_START, MONTH_START, DAY_START, &
                   HOUR_START,MINUTE_START, SECOND_RELAX)
              !ub_ch_20170501+
              IF ( ABS(NOW-NOW_RELAX) .le. tiny) then 
              !IF (NOW .eq. NOW_RELAX) THEN
              !ub_ch_20170501-
                 SET(i1)%lrelax = .TRUE.
              ELSE
                 SET(i1)%lrelax = .FALSE.  
              END IF
           ELSE !ini1step=F 
              IF (SET(i1)%lemis) THEN
                 SET(i1)%lrelax = &
                      event_state(SET(i1)%RELAX_EVENT, previous_date)
              END IF
           END IF
           !set relaxation property for tracer        
           IF (SET(i1)%lrelax) THEN
              DO i2=1,SET(i1)%NUMTRAC
                 CALL set_tracer(status, GPTRSTR, SET(i1)%idt(i2), I_RELAX, ON)
                 CALL tracer_halt(substr, status)
              END DO
              SET(i1)%RELAX=.FALSE.
           END IF
        END IF
#endif

     END DO setloop
     CALL end_message_bi(modstr, 'ptracini global_start', substr)

   END SUBROUTINE ptracini_global_start
! *********************************************************

   SUBROUTINE ptracini_vdiff
     ! ECHAM5/MESSy
     USE messy_main_grid_def_mem_bi, ONLY: kproma, jrow
     USE messy_main_data_bi,     ONLY: tm1, press_3d, qm1_3d, um1, vm1
  
     ! MESSy
     USE messy_main_constants_mem, ONLY: cp_air,rd,alv
     
     IMPLICIT NONE
     INTEGER :: jp
     !calculate additional fields
     IF (lwind) THEN 
        DO jp=1, kproma
           CALL calc_wind(wind(_RI_XYZ__(1:kproma,jrow,:))           &
                ,um1(_RI_XYZ__(1:kproma,jrow,:)),vm1(_RI_XYZ__(1:kproma,jrow,:)))
        END DO
     END IF
     IF (lthetae) THEN 
        DO jp=1, kproma
           CALL calc_thetae(thetae(_RI_XYZ__(1:kproma,jrow,:))  &
                ,tm1(_RI_XYZ__(1:kproma,jrow,:))                &
                ,press_3d(_RI_XYZ__(1:kproma,jrow,:))           &
                ,qm1_3d(_RI_XYZ__(1:kproma,jrow,:)),rd,cp_air,alv)
        END DO
     END IF

   END SUBROUTINE ptracini_vdiff

   ! -------------------------------------------------------------------
   SUBROUTINE ptracini_physc
     
     ! MESSY BMIL
     USE messy_main_grid_def_mem_bi, ONLY: kproma,jrow,nlev,ngpblks,nproma
     USE messy_main_data_bi,       ONLY: eps
     USE messy_main_tracer_mem_bi, ONLY: xtte=>qxtte,xtm1=>qxtm1
     USE messy_main_timer,         ONLY: time_step_len
     
     IMPLICIT NONE

     INTEGER :: i1, i2, i3
     INTEGER :: jp, jk
     INTEGER :: idt
     LOGICAL :: lexit,lcompfld, lcompfld2D, loutfield
     
     REAL(dp), DIMENSION(:,:,:), POINTER     :: comppointer => NULL()  
     REAL(dp), DIMENSION(:,:,:), POINTER     :: comparray => NULL()
     REAL(DP)                                :: dt, conc
    
     
     ALLOCATE(comparray(_RI_XYZ__(nproma,ngpblks,nlev)))
     comparray(:,:,:) = 0._dp
     setloop:DO i1= 1, NUMSETS
        IF (.not.SET(i1)%lemis) CYCLE setloop
        
        jploop: DO jp = 1, kproma 
           conc = 1.e-7
           loutfield=.FALSE.
           IF (SET(i1)%DOMAIN) THEN !does gridpoint belong to DOMAINFLD?
              IF (.not. SET(i1)%IFNOT) THEN
                 IF (input_DOMAINFLD(jp,jrow) .eq. 0._dp) CYCLE jploop
              ELSE
                 IF (input_DOMAINFLD(jp,jrow) .eq. 0._dp) THEN
                    loutfield=.TRUE.
                 END IF
              END IF
           END IF
           jkloop: DO jk = 1, nlev   
              lexit      = .FALSE.
              conc       = 1.e-7 ! ub_ak_20190705
              critloop: DO i3 = 1, SET(i1)%NUMCRIT
                 lcompfld   = .FALSE.
                 lcompfld2D = .FALSE.
                 IF (loutfield) THEN
                    lexit = .TRUE.
                    EXIT !ini anyway outside of compfield (ifnot=T)
                 ELSEIF (ASSOCIATED(SET(i1)%CRT(i3)%ptr)) THEN
                    IF (ASSOCIATED(SET(i1)%CRT_comp2D(i3)%ptr)) THEN
                       !comp = fld2D => SPREAD!
                       lcompfld2D = .TRUE.
#ifdef ECHAM5
                       comparray(:,:,:)= &              
                            SPREAD(SET(i1)%CRT_comp2D(i3)%ptr,2,nlev)
#else
                       comparray(:,:,:)= &             
                            SPREAD(SET(i1)%CRT_comp2D(i3)%ptr,3,nlev)
#endif
                       comppointer => comparray
                    ELSEIF (ASSOCIATED(SET(i1)%CRT_comp(i3)%ptr)) THEN
                       !comp = fld3D
                       lcompfld = .TRUE.
                       comppointer => SET(i1)%CRT_comp(i3)%ptr
                    END IF
                    
                    SELECT CASE (TRIM(SET(i1)%CRIT(i3)%operator))
                    CASE('>')
                       ! exit loop "if not >"    
                       IF ((lcompfld) .OR. (lcompfld2D)) THEN
                          IF (SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk)) <= &
                               comppointer(_RI_XYZ__(jp,jrow,jk))) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       ELSE !comp=const
                          IF (SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk)) <= &
                               SET(i1)%CRIT(i3)%comp%const) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       END IF
                    CASE('>=')
                       IF ((lcompfld) .OR. (lcompfld2D)) THEN
                          IF (SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk)) < &     
                               comppointer(_RI_XYZ__(jp,jrow,jk))) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       ELSE
                          IF (SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk)) < &
                               SET(i1)%CRIT(i3)%comp%const) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       END IF
                    CASE('<')
                       IF ((lcompfld) .OR. (lcompfld2D)) THEN
                          IF (SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk)) >= &   
                               comppointer(_RI_XYZ__(jp,jrow,jk))) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       ELSE
                          IF (SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk)) >= &
                               SET(i1)%CRIT(i3)%comp%const) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       END IF
                    CASE('<=')
                       IF ((lcompfld) .OR. (lcompfld2D)) THEN
                          IF (SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk)) > &     
                               comppointer(_RI_XYZ__(jp,jrow,jk))) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       ELSE
                          IF (SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk)) > &
                               SET(i1)%CRIT(i3)%comp%const) THEN
                             lexit = .TRUE.
                             EXIT
                          ENDIF
                       END IF
                    CASE('IF')
                       conc = SET(i1)%CRT(i3)%ptr(_RI_XYZ__(jp,jrow,jk))
                    CASE DEFAULT
                       ! SHOULD NEVER BE REACHED
                       ! NOTHING TO DO
                    END SELECT
                 ELSE
                    !field = 2d:
                    ! => no spread necessary 
                    SELECT CASE (TRIM(SET(i1)%CRIT(i3)%operator))
                    CASE('>')
                       ! exit loop "if not >"  
                       IF ((SET(i1)%CRIT(i3)%comp%const) < -9990._dp) THEN
                          !comp=fld_2d
                          IF (SET(i1)%CRT_2D(i3)%ptr(jp,jrow) <= &     
                               SET(i1)%CRT_comp2D(i3)%ptr(jp,jrow)) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       ELSE
                          !comp=const
                          IF (SET(i1)%CRT_2D(i3)%ptr(jp,jrow) <= &     
                               SET(i1)%CRIT(i3)%comp%const) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       END IF
                    CASE('>=')
                       IF ((SET(i1)%CRIT(i3)%comp%const) < -9990._dp) THEN
                          !comp=fld_2d
                          IF (SET(i1)%CRT_2D(i3)%ptr(jp,jrow) < &     
                                  SET(i1)%CRT_comp2D(i3)%ptr(jp,jrow)) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       ELSE
                          !comp=const
                          IF (SET(i1)%CRT_2D(i3)%ptr(jp,jrow) < &     
                               SET(i1)%CRIT(i3)%comp%const) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       END IF
                    CASE('<')
                       IF ((SET(i1)%CRIT(i3)%comp%const) < -9990._dp) THEN
                          !comp=fld_2d
                          IF (SET(i1)%CRT_2D(i3)%ptr(jp,jrow) >= &     
                                  SET(i1)%CRT_comp2D(i3)%ptr(jp,jrow)) THEN
                             lexit = .TRUE.
                             EXIT
                          END IF
                       ELSE
                          !comp=const
                          IF (SET(i1)%CRT_2D(i3)%ptr(jp,jrow) >= &     
                               SET(i1)%CRIT(i3)%comp%const) THEN
                             lexit = .TRUE.
                                EXIT
                             END IF
                          END IF
                       CASE('<=')
                          IF ((SET(i1)%CRIT(i3)%comp%const) < -9990._dp) THEN
                             !comp=fld_2d
                             IF (SET(i1)%CRT_2D(i3)%ptr(jp,jrow) > &     
                                  SET(i1)%CRT_comp2D(i3)%ptr(jp,jrow)) THEN
                                lexit = .TRUE.
                                EXIT
                             END IF
                          ELSE
                             !comp=const
                             IF (SET(i1)%CRT_2D(i3)%ptr(jp,jrow) > &     
                                  SET(i1)%CRIT(i3)%comp%const) THEN
                                lexit = .TRUE.
                                EXIT
                             END IF
                          END IF
                       CASE('IF') ! ub_ak_20190705
                          conc = SET(i1)%CRT_2D(i3)%ptr(jp,jrow) 
                       CASE DEFAULT
                          ! SHOULD NEVER BE REACHED
                          ! NOTHING TO DO
                       END SELECT
                    END IF
                 END DO critloop
                 ! no initialisation in this box
                 IF (lexit) THEN   !crit not fullfilled
                    IF (.not. SET(i1)%IFNOT) CYCLE
                 ELSE !crit fullfilled
                    IF (SET(i1)%IFNOT) CYCLE
                 END IF !crit fullfilled

                 ! NOW initialise TRACER
                 tracloop: DO i2=1, SET(i1)%NUMTRAC
                    idt  = SET(i1)%idt(i2)
                    dt   = time_step_len
                    !conc = 1.e-7_dp
                    IF (.not. l2tls) then !leapfrog => ini also 2nd step  
                       IF (SET(i1)%INI1STEP) THEN
                          IF (SET(i1)%linistart) THEN !lstart=T
                             xtte(_RI_X_ZN_(jp,jk,idt)) = conc/dt
                          ELSE IF (SET(i1)%lini2nd) THEN
                             xtte(_RI_X_ZN_(jp,jk,idt)) = (conc*111._dp/100._dp)/dt
                          ELSE !not lstart
                             IF (.not.SET(i1)%lfirst) THEN !=> ini. 1st step
                                xtte(_RI_X_ZN_(jp,jk,idt)) = (conc*20._dp/19._dp)/dt
                             ELSE !=> initialisation 2nd step 
                                xtte(_RI_X_ZN_(jp,jk,idt)) = (conc*18._dp/19._dp)/dt
                             END IF
                          END IF
                       ELSE !ini each step => concentr. is set back to conc
                          xtm1(_RI_X_ZN_(jp,jk,idt)) = (1._DP -eps)*conc                 
                          xtte(_RI_X_ZN_(jp,jk,idt)) = (-1._DP * xtm1(_RI_X_ZN_(jp, jk,idt)) + &
                               conc)/dt
                       END IF
                    ELSE !no leapfrog
                       xtte(_RI_X_ZN_(jp,jk,idt)) = xtte(_RI_X_ZN_(jp,jk,idt))+ &
                            (conc - (xtm1(_RI_X_ZN_(jp, jk,idt)) +    &
                            xtte(_RI_X_ZN_(jp,jk,idt))*dt))/dt  
                    END IF
                 END DO tracloop
              END DO jkloop
           END DO jploop
        END DO setloop
        
        
        ! CLEAN UP
        DEALLOCATE(comparray)
     
   END SUBROUTINE ptracini_physc
   ! -------------------------------------------------------------------
   SUBROUTINE ptracini_read_nml_cpl(status, iou)

     ! ptracini MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
     !
     ! read namelist for 'coupling' to ECHAM5
     !
     ! Author: Patrick Joeckel, MPICH, Dec 2004

     ! MESSy
     USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
     
    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status    ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'ptracini_read_nml_cpl'
    
    NAMELIST /CPL/ DOMAINFLD,lwind,lthetae,TRINI

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    
    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    ! NOTE: already at definition
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    CALL read_nml_close(substr, iou, modstr)
    

    status = 0  ! no ERROR

  END SUBROUTINE ptracini_read_nml_cpl
  ! -------------------------------------------------------------------
  SUBROUTINE calc_wind(wind,u,v)
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: wind  ! horizontal windspeed [m/s]
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: u     ! U-component of wind [m/s]
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: v     ! V-component of wind [m/s]

    INTRINSIC :: sqrt
    
    wind=sqrt(u**2+ v**2)
    
  END SUBROUTINE calc_wind
  ! -------------------------------------------------------------------
  SUBROUTINE calc_thetae(thetae,t,p,q,r_d,cp_d,lh_v)
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: thetae! equi. pot. temp. [K]
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: t     ! temperature [K]
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: p     ! pressure [Pa]
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: q     !specific humidity [kg/kg]
    REAL(DP), INTENT(IN), OPTIONAL :: r_d  ! gas const. of dry air [J/K/kg]
    REAL(DP), INTENT(IN), OPTIONAL :: cp_d ! spec.heat of dry air at
                                           ! constant pressure [J/K/kg]
    REAL(dp), INTENT(IN), OPTIONAL :: lh_v ! latent heat for vaporisation[J/kg]
  
    thetae = (100000._dp/p)**(r_d/cp_d)*t*exp(lh_v/cp_d*q/t) 
    
  END SUBROUTINE calc_thetae
  !--------------------------------------------------------------------
  
! **********************************************************************
END MODULE messy_ptracini_si
! **********************************************************************
