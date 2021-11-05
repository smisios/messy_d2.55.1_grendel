! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL OASIS3-MCT
!
! Author : Astrid Kerkweg, Uni Bonn, 2018
!
!
! 0) configure with COSMO
!./configure --enable-COSMO --enable-MESSYMMD --enable-OASIS3MCT --disable-ECHAM5
! 1) configure with COSMO and CLM
!./configure --enable-COSMO --enable-MESSYMMD --enable-OASIS3MCT --disable-ECHAM5 --enable-CLM
! 2) configure MECO(n) system with COSMO and CLM using COSMO master
!./configure --enable-COSMO --enable-I2CINC --enable-MESSYMMD --enable-OASIS3MCT --disable-ECHAM5 --enable-CLM
! 3) configure MECO(n) system including CLM
!./configure --enable-COSMO --enable-I2CINC --enable-MESSYMMD --enable-OASIS3MCT --enable-CLM
!
! **********************************************************************
! TODO:
! --- RESTART TEST

! FESOM:
!   - masks  (calculate,  make channel objects?)
!       cpl_oce_mask
!       cpl_ice_mask, cpl_wat_mask
!   - 2 grids!
!   - input : a2o field ?
!   - out hardcoded block  ( calc cpl_ice_mask, cpl_wat_mask)
!  WARNING: the extra communication exchanging net fluxes is omitted (but required for conservative
!           remapping. This needs to be made obsolete
! **********************************************************************

#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_oasis3mct_si
! **********************************************************************
!!$#ifdef COSMO
!!$#define WRITEGRIDPARALLEL 1
!!$#endif

 ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: error_bi
#ifdef OASIS3MCT
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      warning_bi
  USE messy_main_data_bi,       ONLY: lcpl_land
#endif

  USE messy_main_blather_bi,    ONLY: error_bi
  ! SMCL
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, t_chaobj_cpl
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM &
                                    , STRLEN_VLONG, iouerr
  USE messy_oasis3mct

#ifdef OASIS3MCT
  USE messy_main_channel,       ONLY: STRLEN_OBJECT

  ! EXTERNAL LIBRARY
  USE MOD_OASIS
#endif

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! COMPONENT
#ifdef OASIS3MCT
  INTEGER                :: compid    = 0       ! component ID
  INTEGER                :: cnlen      = 0      ! compname length

  ! component participates in coupling
  INTEGER                :: local_comm = -9999  ! communicator for component
  INTEGER, ALLOCATABLE, DIMENSION(:) :: paral

  ! partition ID
  INTEGER                     :: part_id = 0

  ! start and end index of coupled field
  INTEGER               :: iis, iie, jjs, jje, xdim_tot, ydim_tot
#ifdef ECHAM5
  INTEGER               :: nsegs    ! number of segments
  INTEGER               :: nsegslen ! length of segments
#endif

  ! DEFINE ' EXTERNAL' land-sea fraction
  ! for COSMO with I2CINC we need to read this in OASIS3MCT submodel, as it is
  ! known too late in COSMO due to dependencies with MMD
  CHARACTER(LEN=STRLEN_VLONG)  :: slf_file = ' '
  CHARACTER(LEN=STRLEN_OBJECT) :: slf_name = ' '
  CHARACTER(LEN=STRLEN_OBJECT) :: lon_name = ' '
  CHARACTER(LEN=STRLEN_OBJECT) :: lat_name = ' '
#endif
  ! ------------------------------------------------------------------
  ! COUPLE MODEL
  ! ------------------------------------------------------------------
  TYPE t_cpl_model_io
     CHARACTER(LEN=STRLEN_CHANNEL) :: modname= ''
     CHARACTER(LEN=STRLEN_MEDIUM)  :: locsnd = ''
     CHARACTER(LEN=STRLEN_MEDIUM)  :: locrcv = ''
     INTEGER                       :: imask  = -99
  END type t_cpl_model_io

  TYPE t_cpl_model
     CHARACTER(LEN=STRLEN_CHANNEL) :: modname   = ''
     INTEGER                       :: locsndSR  = -99 ! subroutine to send data
     INTEGER                       :: locsndBE  = -99 ! call at begin / end of subroutine
     INTEGER                       :: locsndord = -99 ! order of send / receive
     INTEGER                       :: locrcvSR  = -99 ! subroutine to recv data
     INTEGER                       :: locrcvBE  = -99 ! call at begin / end of subroutine
     INTEGER                       :: locrcvord = -99 ! order of send / receive
     LOGICAL                       :: lactive   = .FALSE.
     INTEGER                       :: masktype  = -99
     LOGICAL,DIMENSION(:,:), POINTER :: mask => NULL()
  END type t_cpl_model

  INTEGER,              PARAMETER                 :: MAXNUM_MODELS = 20

#ifdef OASIS3MCT
  ! LOCATION FOR COMPONENT COUPLING
  TYPE(t_cpl_model_io), DIMENSION(MAXNUM_MODELS)  :: LOCCOUPLE
  TYPE(t_cpl_model),    ALLOCATABLE, DIMENSION(:) :: LOCCPL

  ! NUMBER OF COMPONENT MODELS (dimension of LOCCPL)
  INTEGER  :: NUMMOD
#endif
  ! ------------------------------------------------------------------
  ! COUPLE FIELD
  ! ------------------------------------------------------------------
 TYPE t_cpl_rfield_io
#ifdef OASIS3MCT
     CHARACTER(LEN=ic_lvar)       :: oasisname = ' '
#endif
     TYPE(t_chaobj_cpl)           :: messyname
     LOGICAL                      :: lcreatine = .FALSE. ! create if not exists
     ! Mask field 0: all points are used, 1: land points 2: sea points
     INTEGER                      :: imask = -99         ! mask field
     ! ACTION STRING DEFINES THE POSTPROCESSING OF THE FIELD
     ! non-empty variable needs specific treatment
     !    - ABS      use absolute value
     !    - x:fac    multiply by "fac", e.g. "-1"
     !    - +:const  add constant
     CHARACTER(LEN=STRLEN_VLONG)  :: postproc    = ' ' ! variable needs to
     ! DEFINES HOW TO INTERPRET THE OASIS MULTI-LEVEL INFORMATION
     ! FORMAT:  'Z;N:nnum;T:tlev' => respective order is important !
     !           'Z' can be without specifier => nlev
     CHARACTER(LEN=STRLEN_MEDIUM) :: oasaxis   = ' '
     !
     ! SPECIFIES LEVEL of variable , enables rank reduction
     ! timelev:  nnow: -1 ; nnew: -2 ; nold: -3
     CHARACTER(LEN=STRLEN_MEDIUM) :: levelstr  = ' '     ! if variable with time
                                                         ! be calculated
     CHARACTER(LEN=STRLEN_MEDIUM) :: repr      =' '
  END type t_cpl_rfield_io

 TYPE t_cpl_sfield_io
#ifdef OASIS3MCT
     CHARACTER(LEN=ic_lvar)       :: oasisname = ' '
#endif
     TYPE(t_chaobj_cpl)           :: messyname
     ! DEFINES HOW TO INTERPRET THE OASIS MULTI-LEVEL INFORMATION
     ! FORMAT:  'Z;N:nnum;T:tlev' => respective order is important !
     !           'Z' can be without specifier => nlev
     CHARACTER(LEN=STRLEN_MEDIUM) :: oasaxis   = ' '
     !
     ! SPECIFIES LEVEL of variable , enables rank reduction
     ! timelev:  nnow: -1 ; nnew: -2 ; nold: -3
     CHARACTER(LEN=STRLEN_MEDIUM) :: levelstr  = ' '     ! if variable with time
                                                         ! be calculated
   END type t_cpl_sfield_io


   TYPE t_postproc_info
      CHARACTER(LEN=STRLEN_MEDIUM) :: action = ' '
      REAL(dp)                     :: number = -9999.e-36
      REAL(dp)                     :: limit  = -9999.e-36
   END type t_postproc_info

  TYPE t_cpl_rfield
#ifdef OASIS3MCT
     CHARACTER(LEN=ic_lvar)       :: oasisname = ' '
#endif
     TYPE(t_chaobj_cpl)           :: messyname
     LOGICAL                      :: lcreatine = .FALSE.  ! create if not exists
     CHARACTER(LEN=STRLEN_MEDIUM) :: repr      = ' '
     INTEGER                      :: vertlev   = -99
     INTEGER                      :: numblev   = -99
     INTEGER                      :: timelev   = -99
     LOGICAL                      :: lcalc     = .FALSE.
     INTEGER                      :: cplcompID = -99 ! id of coupled component
                                               ! corresponds to LOCCOUPLE IDX
     INTEGER                      :: reprID    = -99
     INTEGER                      :: oasvarID  = -99
     INTEGER                      :: rank = -99
     REAL(dp), POINTER, DIMENSION(:,:)     :: store => NULL()
     REAL(dp), POINTER, DIMENSION(:,:)     :: ptr2d => NULL()
     REAL(dp), POINTER, DIMENSION(:,:,:)   :: ptr3d => NULL()
     INTEGER                      :: idt = -99 ! TEST
     ! POST-PROCESS FIELD
     ! mask field?
     ! 0: all points
     ! 1: land points only
     ! 2: sea points only
     INTEGER                                          :: imask
     TYPE(t_postproc_info), DIMENSION(:), ALLOCATABLE :: pp
 END type t_cpl_rfield

  TYPE t_cpl_sfield
#ifdef OASIS3MCT
     CHARACTER(LEN=ic_lvar)       :: oasisname = ' '
#endif
     TYPE(t_chaobj_cpl)           :: messyname
     INTEGER                      :: cplcompID = -99 ! id of coupled component
                                               ! corresponds to LOCCOUPLE IDX
     INTEGER                      :: oasvarID  = -99
     INTEGER                      :: vertlev   = -99
     INTEGER                      :: numblev   = -99
     INTEGER                      :: timelev   = -99
     INTEGER                      :: rank = -99
     REAL(dp), POINTER, DIMENSION(:,:)     :: ptr2d => NULL()
     REAL(dp), POINTER, DIMENSION(:,:,:)   :: ptr3d => NULL()
     INTEGER                      :: idt = -99 ! TEST
 END type t_cpl_sfield

#ifdef OASIS3MCT
  INTEGER              , PARAMETER                 :: MAXNUM_FIELDS = 100
  TYPE(t_cpl_rfield_io), DIMENSION(MAXNUM_FIELDS)  :: RFIELD
  TYPE(t_cpl_rfield),    ALLOCATABLE, DIMENSION(:) :: RFLD
  TYPE(t_cpl_sfield_io), DIMENSION(MAXNUM_FIELDS)  :: SFIELD
  TYPE(t_cpl_sfield),    ALLOCATABLE, DIMENSION(:) :: SFLD

  ! NUMBER OF FIELDS COMPONENT MODELS (dimension of FLD)
  INTEGER  :: NUMRFLD
  INTEGER  :: NUMSFLD

  ! LOGICAL LANDMASK
  LOGICAL,  DIMENSION(:,:), POINTER   :: llandmask  => NULL()
  ! ------------------------------------------------------------------

  ! DATA REQUIRED FOR SPECIFIC POSTPROCESSINGS
  INTEGER :: isoil_external_ppSR = -99
  INTEGER :: isoil_external_ppBE = -99
  INTEGER :: isoil_external_cplfreq = -99
  ! ------------------------------------------------------------------
#endif


  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_e5.f90
  PUBLIC :: oasis3mct_setup         ! setup OASIS interface
  PUBLIC :: oasis3mct_initialize    ! initialize submodel
  PUBLIC :: oasis3mct_init_memory   ! request memory
  PUBLIC :: oasis3mct_init_coupler  ! wrapper init_memory / init_coupling
  PUBLIC :: oasis3mct_init_coupling ! set pointers for coupling to BM and other SMs
  PUBLIC :: oasis3mct_time            ! currently used in ECHAM5
  PUBLIC :: oasis3mct_init_loop     ! currently used in COSMO
  PUBLIC :: oasis3mct_global_start  ! currently used in ECHAM5
  PUBLIC :: oasis3mct_convec        ! currently used in COSMO
  PUBLIC :: oasis3mct_write_output   ! currently used in ECHAM (for FESOM2)


  PUBLIC :: oasis3mct_finalize  ! TODO: implement CALL to subroutine!!

  ! PRIVATE SUBROUTINES
  ! PRIVATE :: oasis3mct_read_nml_cpl
  ! PRIVATE :: oasis3mct_read_nml_cpl_extdata
  ! PRIVATE :: oasis3mct_sndrcv    ! Subroutine for data exchange
                                   ! CALL from each possible entry point

  ! PRIVATE :: parse_location

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################
  ! ====================================================================
  SUBROUTINE oasis3mct_setup(lcoupled)

#ifdef MESSYMMD
    ! MMD
    USE mmd_handle_communicator, ONLY: MMD_instance_name  &
                                     , MMD_numberofNonMMDcoupledmodels
#endif
    IMPLICIT NONE

    ! I/O
    LOGICAL,          INTENT(in)  :: lcoupled
#ifdef OASIS3MCT
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr='oasis3mct_setup'
    INTEGER                       :: status
    LOGICAL                       :: lcpl = .TRUE.

    CALL start_message_bi(modstr,'SETUP',substr)  ! log-output

    IF (MMD_numberofNonMMDcoupledmodels() == 0) THEN
       lcpl = .FALSE.
    ELSE

       compname=MMD_instance_name()
       write (*,*) 'OASIS couple-name of model: ', TRIM(compname),lcoupled

       CALL oasis_init_comp(compid, compname, status, lcoupled)
    END IF

    IF (.NOT. lcoupled .OR. .NOT. lcpl) THEN
        CALL end_message_bi(modstr,'SETUP - NOT COUPLED',substr)
       RETURN
    END IF
    IF (status /= OASIS_OK) THEN
       CALL OASIS_abort(cd_routine='oasis3mct_setup' &
            , cd_message='ERROR in OASIS_init_comp', rcode=status)
    END IF

    CALL oasis_get_localcomm(local_comm,status)
    IF (status /= OASIS_OK) CALL OASIS_abort(cd_routine='oasis3mct_setup' &
            , cd_message='ERROR in OASIS_get_localcomm',rcode=status)

    cnlen = LEN_TRIM(compname)

    CALL end_message_bi(modstr,'SETUP',substr)  ! log-output

#endif
  END SUBROUTINE oasis3mct_setup
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE oasis3mct_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

#ifdef OASIS3MCT
    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,          ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_grid_def_mem_bi, ONLY: nlev
    USE messy_main_data_bi,         ONLY: L_IS_CHILD
    USE messy_main_tools,           ONLY: find_next_free_unit, strcrack, str2num
    USE messy_main_control,         ONLY: entry_name, MEPS_BEGIN
#endif

    IMPLICIT NONE

    ! LOCAL
#ifdef OASIS3MCT
    CHARACTER(LEN=*), PARAMETER :: substr = 'oasis3mct_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: ix, iy
    INTEGER                     :: ioas
    INTEGER,                  DIMENSION(ig_final_nfield) :: cplcompid
    INTEGER,                  DIMENSION(ig_final_nfield) :: cplsubid
    INTEGER,                  DIMENSION(ig_final_nfield) :: cpllevid
    CHARACTER(LEN=ic_lvar),   DIMENSION(ig_final_nfield) :: oasisname
    INTEGER                     :: jx, jxlen
    LOGICAL                                       :: lfound = .FALSE.
    CHARACTER(LEN=ic_lvar), DIMENSION(:), POINTER :: tmp1  => NULL()
    CHARACTER(LEN=ic_lvar), DIMENSION(:), POINTER :: tmp2  => NULL()
    INTEGER                                        :: nl, ni, is , il
    INTEGER                                        :: naxlen, taxlen, zaxlen
    CHARACTER(LEN=3)                               :: axstr = '   '
    INTEGER                                        :: ppsize

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CPL namelist
    IF (p_parallel_io) THEN                      ! read only on I/O-PE
       iou = find_next_free_unit(100,200)        ! find next free I/O unit
       CALL oasis3mct_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF

    ! ANALYSE LOCCOUPLE (setup of coupling models)
    CALL start_message_bi(modstr,'ANALYSE LOCCOUPLE',substr)  ! log-output
    nummod = 0
    IF (p_parallel_io) THEN
       DO ix = 1, MAXNUM_MODELS
          IF ( TRIM(ADJUSTL(LOCCOUPLE(ix)%modname)) /= '') THEN
             nummod = nummod +1
          END IF
       END DO

       ALLOCATE(LOCCPL(nummod))
       nummod   = 0
       DO ix = 1, MAXNUM_MODELS
          IF ( TRIM(ADJUSTL(LOCCOUPLE(ix)%modname)) == '') CYCLE
          nummod = nummod + 1
          LOCCPL(nummod)%modname = TRIM(ADJUSTL(LOCCOUPLE(ix)%modname))
          CALL parse_location(status, LOCCOUPLE(ix)%locsnd &
               , LOCCPL(nummod)%locsndSR,LOCCPL(nummod)%locsndBE &
               , LOCCPL(nummod)%locsndord)
          CALL parse_location(status, LOCCOUPLE(ix)%locrcv &
               , LOCCPL(nummod)%locrcvSR,LOCCPL(nummod)%locrcvBE &
               , LOCCPL(nummod)%locrcvord)

          write (*,*) 'COUPLE MODEL ', TRIM(LOCCPL(nummod)%modname)
          IF (LOCCPL(nummod)%locsndBE == MEPS_BEGIN) THEN
             write (*,*) '.. send at begin in ', entry_name(LOCCPL(nummod)%locsndSR)
          ELSE
             write (*,*) '.. send at end of ', entry_name(LOCCPL(nummod)%locsndSR)
          END IF
          IF (LOCCPL(nummod)%locrcvBE == MEPS_BEGIN) THEN
             write (*,*) '.. receive at begin in ', entry_name(LOCCPL(nummod)%locrcvSR)
          ELSE
             write (*,*) '.. receive at end of ', entry_name(LOCCPL(nummod)%locrcvSR)
          END IF
          write (*,*) '.. order send / receive ' &
               , LOCCPL(nummod)%locsndord, LOCCPL(nummod)%locrcvord
          !
          LOCCPL(nummod)%masktype = LOCCOUPLE(ix)%imask
       END DO
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs

    CALL p_bcast(nummod, p_io)
    IF (.NOT. p_parallel_io) ALLOCATE(LOCCPL(nummod))
    DO ix = 1, nummod
       CALL p_bcast(LOCCPL(ix)%modname,   p_io)
       CALL p_bcast(LOCCPL(ix)%locsndSR,  p_io)
       CALL p_bcast(LOCCPL(ix)%locsndBE,  p_io)
       CALL p_bcast(LOCCPL(ix)%locsndord, p_io)
       CALL p_bcast(LOCCPL(ix)%locrcvSR,  p_io)
       CALL p_bcast(LOCCPL(ix)%locrcvBE,  p_io)
       CALL p_bcast(LOCCPL(ix)%locrcvord, p_io)
       CALL p_bcast(LOCCPL(ix)%masktype,  p_io)
       LOCCPL(ix)%lactive = .FALSE.
       LOCCPL(ix)%mask => NULL()
    END DO
    CALL end_message_bi(modstr,'ANALYSE LOCCOUPLE',substr)  ! log-output

    !------------------------------------------------------------------
    ! ANALYSE RFIELDs
    CALL start_message_bi(modstr,'ANALYSE RFIELDS',substr)  ! log-output
    IF (p_parallel_io) THEN
       cplsubid  = -99
       cpllevid  = -99
       cplcompID = -99
       oasisname = ' '
       numrfld = 0
       DO ioas = 1, ig_final_nfield
          lfound = .FALSE.

          ! split OASIS name in level and basename
          CALL strcrack(namdstfld(ioas),'_',tmp1,nl)

          ! CHECK if OASISname exists in namcouple namelist.
          ! IF not, field does not have to be taken into account
          DO ix = 1, MAXNUM_FIELDS
             IF ( TRIM(ADJUSTL(RFIELD(ix)%OASISname)) /= '') THEN

                IF  ((TRIM(tmp1(1)) ==                   &
                     TRIM(ADJUSTL(RFIELD(ix)%oasisname)) )       &
                     .AND. &
                     (namdstgrd(ioas)(1:cnlen) == TRIM(compname))      &
                     ) THEN

                   IF ( TRIM(ADJUSTL(RFIELD(ix)%messyname%cha)) == 'oasis3mct' &
                        .OR. RFIELD(ix)%lcreatine) THEN
                      ! representation info is required for channel object
                      ! creation
                      IF (TRIM(ADJUSTL(RFIELD(ix)%repr)) == '') THEN
                         CYCLE
                      END IF
                   END IF

                   DO jx = 1,nummod
                      jxlen = LEN_TRIM(LOCCPL(jx)%modname)
                      IF (TRIM(LOCCPL(jx)%modname) == namsrcgrd(ioas)(1:jxlen) &
                           ) THEN
                         IF (nl > 1) THEN
                            CALL str2num(tmp1(2),cpllevid(ioas), status)
                            IF (status /= 0) CALL error_bi( &
                                 'ERROR in string conversion', substr)
                         END IF
                         LOCCPL(jx)%lactive = .TRUE.
                         oasisname(ioas) = TRIM(namdstfld(ioas))
                         cplcompid(ioas) = jx
                         cplsubid(ioas)   = ix
                         ! establish postprocessing flags
                         ! for Community Land Model coupling
                         IF (LOCCPL(jx)%modname(1:3) =='clm') THEN
                            lcpl_land = .TRUE.
                            isoil_external_ppSR =  LOCCPL(jx)%locrcvSR
                            isoil_external_ppBE =  LOCCPL(jx)%locrcvBE
                         END IF
                         EXIT
                      END IF
                   END DO
                   IF (cplcompid(ioas) > 0) THEN
                      numrfld = numrfld +1
                      lfound = .TRUE.
                      EXIT
                   END IF
                END IF
             END IF
          END DO

          IF (.NOT. lfound .AND. namdstgrd(ioas)(1:cnlen) == TRIM(compname)) THEN
             write(iouerr,*) 'OASIS namelist field: ',TRIM(namdstfld(ioas))
             CALL error_bi('no matching entry for receive field found!', substr)
          END IF
          IF (ASSOCIATED(tmp1)) DEALLOCATE(tmp1)
          NULLIFY(tmp1)

       END DO

       !----------------------------------------------------------

       write (*,*)  '-> NUMBER RECEIVE FIELDS FOUND ',numrfld
       ALLOCATE(RFLD(numrfld))
       numrfld   = 0

       DO ioas = 1, ig_final_nfield

          IF ( cplcompid(ioas) < 0) CYCLE
          numrfld = numrfld + 1
          RFLD(numrfld)%OASISname     = TRIM(ADJUSTL(oasisname(ioas)))
          RFLD(numrfld)%MESSyname%CHA = &
               TRIM(ADJUSTL(RFIELD(cplsubid(ioas))%MESSyname%CHA))
          RFLD(numrfld)%MESSyname%OBJ = &
               TRIM(ADJUSTL(RFIELD(cplsubid(ioas))%MESSyname%OBJ))
          RFLD(numrfld)%lcreatine     = RFIELD(cplsubid(ioas))%lcreatine
          RFLD(numrfld)%imask         = RFIELD(cplsubid(ioas))%imask
          RFLD(numrfld)%repr          = &
               TRIM(ADJUSTL(RFIELD(cplsubid(ioas))%repr))


          ! LOGFILE OUTPUT INFO
          write (*,*) '--> ASSOCIATE OASIS RECEIVE FIELD ' &
               , TRIM(RFLD(numrfld)%OASISname)
          write (*,*) '    to channel object '&
               , TRIM(RFLD(numrfld)%messyname%cha), ' ' &
               , TRIM(RFLD(numrfld)%messyname%obj)

          write (*,*) '    create if channel object does not exist? '&
                          , RFLD(numrfld)%lcreatine
          write (*,*) '    create in representation: ', TRIM(RFLD(numrfld)%repr)

          IF (RFLD(numrfld)%imask == 0) THEN
             write (*,*) ' Field will not be masked: '
          ELSE IF (RFLD(numrfld)%imask == 1) THEN
             write (*,*) ' Field will be restricted to land points: '
          ELSE IF (RFLD(numrfld)%imask == 2) THEN
             write (*,*) ' Field will be restricted to sea points: '
          ELSE
             CALL error_bi('UNKOWN MASK FOR RECEIVE FIELD', substr)
          END IF

          ! analyse postprocessing string
          CALL parse_postproc(RFIELD(cplsubid(ioas))%postproc &
               , RFLD(numrfld))

          ! set cplcompID
          RFLD(numrfld)%cplcompID    = cplcompid(ioas)

           write (*,*) ' FIELD is received from component ' &
               , TRIM(LOCCPL(RFLD(numrfld)%cplcompID)%modname)

           ! INVESTIGATE level string
           ! 1. is OASIS field multi-level source ?
           cpllev: IF (cpllevid(ioas) > 0) THEN
              ! analyse how to interpret cpllevid (T,Z,N -AXIS?)

               IF (TRIM(ADJUSTL(RFIELD(cplsubid(ioas))%oasaxis)) /= '') THEN
                 CALL strcrack(RFIELD(cplsubid(ioas))%oasaxis,';',tmp1,nl)

                  DO il = 1, nl
                    axstr = '   '
                    CALL strcrack(tmp1(il),':', tmp2, ni)

                    IF (ni > 0 .OR. ni < 3) THEN
                       SELECT CASE (tmp2(1))
                       CASE('Z')
                          IF (SIZE(tmp2) == 1) THEN
                             zaxlen = nlev
                          ELSE
                             CALL str2num(tmp2(2), zaxlen, status)
                             IF (status /= 0)  CALL error_bi( &
                                 'error in string conversion! '//tmp2(2), substr)
                          ENDIF
                          axstr(il:il) = 'Z'
                       CASE('N')
                          IF (SIZE(tmp2) <= 1) CALL error_bi(&
                               'n-axis requires length specification', substr)
                          CALL str2num(tmp2(2), naxlen, status)
                          IF (status /= 0)  CALL error_bi( &
                               'error in string conversion! '//tmp2(2), substr)
                          axstr(il:il) = 'N'
                       CASE('T')
                          IF (SIZE(tmp2) <= 1) CALL error_bi(&
                               't-axis requires length specification', substr)
                          CALL str2num(tmp2(2), taxlen, status)
                          IF (status /= 0)  CALL error_bi( &
                               'error in string conversion! '//tmp2(2), substr)
                          axstr(il:il) = 'T'
                       CASE DEFAULT
                          CALL error_bi(' UNKNOWN AXIS IDENTIFIER! '//tmp2(1) &
                               , substr)
                       END SELECT
                    END IF
                    IF (ASSOCIATED(tmp2)) DEALLOCATE(tmp2)
                    NULLIFY(tmp2)
                 END DO
                 IF (ASSOCIATED(tmp1)) DEALLOCATE(tmp1)
                 NULLIFY(tmp1)

                 SELECT CASE(axstr)
                 CASE ('Z  ')
                    IF (cpllevid(ioas) > zaxlen) THEN
                       write(iouerr,*) 'ERROR: level dimensions: ' &
                            ,cpllevid(ioas), zaxlen
                       CALL error_bi('level dimensions do not match', substr)
                    END IF
                    RFLD(numrfld)%vertlev=cpllevid(ioas)
                 CASE ('T  ')
                    IF (cpllevid(ioas) > taxlen) THEN
                       write(iouerr,*) 'ERROR: t-axis length : ' &
                            ,cpllevid(ioas), taxlen
                       CALL error_bi('t-axis length do not match', substr)
                    END IF
                    RFLD(numrfld)%timelev=cpllevid(ioas)
                 CASE ('N  ')
                    IF (cpllevid(ioas) > naxlen) THEN
                       write(iouerr,*) 'ERROR: t-axis length : '&
                            ,cpllevid(ioas), naxlen
                       CALL error_bi('n-axis length do not match', substr)
                    END IF
                    RFLD(numrfld)%numblev=cpllevid(ioas)
                 CASE ('ZT ')
                    RFLD(numrfld)%vertlev= INT(cpllevid(ioas) / zaxlen)
                    RFLD(numrfld)%timelev= MOD(cpllevid(ioas),zaxlen) + 1
                    IF (RFLD(numrfld)%timelev > taxlen) THEN
                       write(iouerr,*) 'ZT-axis dimensions: ' &
                            ,cpllevid(ioas), zaxlen, taxlen
                    END IF
                 CASE ('ZN ')
                    RFLD(numrfld)%vertlev= INT(cpllevid(ioas) / zaxlen)
                    RFLD(numrfld)%numblev= MOD(cpllevid(ioas),zaxlen) + 1
                    IF (RFLD(numrfld)%numblev > naxlen) THEN
                       write(iouerr,*) 'ZN-axis dimensions: ' &
                            ,cpllevid(ioas), zaxlen, naxlen
                    END IF
                 CASE ('TZ ')
                    RFLD(numrfld)%timelev= INT(cpllevid(ioas) / taxlen)
                    RFLD(numrfld)%vertlev= MOD(cpllevid(ioas),taxlen) + 1
                    IF (RFLD(numrfld)%vertlev > zaxlen) THEN
                       write(iouerr,*) 'TZ-axis dimensions: ' &
                            ,cpllevid(ioas), taxlen, zaxlen
                    END IF
                 CASE ('NZ ')
                    RFLD(numrfld)%numblev= INT(cpllevid(ioas) /naxlen)
                    RFLD(numrfld)%vertlev= MOD(cpllevid(ioas),naxlen) + 1
                    IF (RFLD(numrfld)%vertlev > zaxlen) THEN
                       write(iouerr,*) 'NZ-axis dimensions: ' &
                            ,cpllevid(ioas), naxlen, zaxlen
                    END IF
                 CASE ('TN ')
                    RFLD(numrfld)%timelev= INT(cpllevid(ioas) / taxlen)
                    RFLD(numrfld)%numblev= MOD(cpllevid(ioas),taxlen) + 1
                    IF (RFLD(numrfld)%numblev > naxlen) THEN
                       write(iouerr,*) 'TN-axis dimensions: ' &
                            ,cpllevid(ioas), taxlen, naxlen
                    END IF
                 CASE ('NT ')
                    RFLD(numrfld)%numblev= INT(cpllevid(ioas) / naxlen)
                    RFLD(numrfld)%timelev= MOD(cpllevid(ioas),naxlen) + 1
                    IF (RFLD(numrfld)%timelev > taxlen) THEN
                       write(iouerr,*) 'NT-axis dimensions: ' &
                            ,cpllevid(ioas), naxlen, taxlen
                    END IF
              CASE DEFAULT
                 CALL error_bi('unknown OASIS-AXIS string! '//axstr, substr)
              END SELECT

           END IF

        END IF cpllev

        CALL strcrack(TRIM(RFIELD(cplsubid(ioas))%levelstr), ';', tmp1, nl)

        nl1: IF (nl > 0) THEN
           IF (TRIM(tmp1(1)) /= '') THEN
              DO is = 1, nl
                 CALL strcrack(tmp1(is),':', tmp2, ni)
                 IF (ni /=2) CALL error_bi('levelstr is not correct', substr)
                 SELECT CASE (TRIM(tmp2(1)))
                 CASE ('Z')
                    CALL str2num(tmp2(2), il, status)
                    IF ( RFLD(numrfld)%vertlev < 0) THEN
                       IF (il == 0) THEN
                          SFLD(numsfld)%vertlev = nlev
                       ELSE
                          SFLD(numsfld)%vertlev = il
                       END IF
                    ELSE
                       IF (RFLD(numrfld)%vertlev /= il) THEN
                          write(iouerr,*) ' namelist settings inconsistent: ' &
                               , RFLD(numrfld)%vertlev, il
                          CALL error_bi(&
                           ' namelist settings inconsistent for vertical level'&
                           , substr)
                       END IF
                    END IF
                 CASE ('T')
                    CALL str2num(tmp2(2), il, status)
                    IF ( RFLD(numrfld)%timelev < 0) THEN
                       RFLD(numrfld)%timelev = il
                    ELSE
                       IF (RFLD(numrfld)%timelev /= il) THEN
                          write(iouerr,*) ' namelist settings inconsistent: ' &
                               , RFLD(numrfld)%timelev, il
                          CALL error_bi(&
                               ' namelist settings inconsistent for time level' &
                               , substr)
                       END IF
                    END IF
                 CASE ('N')
                    CALL str2num(tmp2(2), il, status)
                    IF ( RFLD(numrfld)%numblev < 0) THEN
                       RFLD(numrfld)%numblev = il
                    ELSE
                       IF (RFLD(numrfld)%numblev /= il) THEN
                          write(iouerr,*) ' namelist settings inconsistent: ' &
                               , RFLD(numrfld)%numblev, il
                          CALL error_bi(&
                            ' namelist settings inconsistent for number level' &
                               , substr)
                       END IF
                    END IF

                 CASE DEFAULT
                    CALL error_bi('unknown level specifier!' //TRIM(tmp2(1)) &
                         , substr)
                 END SELECT

                 IF (ASSOCIATED(tmp2)) DEALLOCATE(tmp2)
                 NULLIFY(tmp2)
              END DO
              IF (ASSOCIATED(tmp1)) DEALLOCATE(tmp1)
              NULLIFY(tmp1)
           END IF
        END IF nl1

        IF ( RFLD(numrfld)%vertlev > 0) &
             write (*,*) '    ...  TO vertical level: ' &
             , RFLD(numrfld)%vertlev

        IF ( RFLD(numrfld)%timelev > -10) &
             write (*,*) '   ...  TO time level: ' &
             , RFLD(numrfld)%timelev

        IF ( RFLD(numrfld)%numblev > 0) &
             write (*,*) '    ... TO numb level: ' &
             , RFLD(numrfld)%numblev

        IF ( RFLD(numrfld)%vertlev > 0 .AND. RFLD(numrfld)%timelev > -10 &
             .AND. RFLD(numrfld)%numblev > 0) THEN
           write(iouerr,*) 'RLFD: TOO many ranks: ', numrfld, RFLD(numrfld)%vertlev &
               , RFLD(numrfld)%timelev , RFLD(numrfld)%numblev
           CALL error_bi( &
             ' RFLD: too many rank reductions requested ', substr)
        END IF
     END DO
  END IF
  ! BROADCAST RFIELD from I/O-PE to ALL OTHER PEs
  CALL p_bcast(lcpl_land,      p_io)
  CALL p_bcast(isoil_external_ppSR, p_io)
  CALL p_bcast(isoil_external_ppBE, p_io)

  CALL p_bcast(numrfld, p_io)
  IF (.NOT. p_parallel_io) ALLOCATE(RFLD(numrfld))
  DO ix = 1, numrfld
     CALL p_bcast(RFLD(ix)%OASISname,      p_io)
     CALL p_bcast(RFLD(ix)%MESSyname%cha,  p_io)
     CALL p_bcast(RFLD(ix)%MESSyname%obj,  p_io)
     CALL p_bcast(RFLD(ix)%lcreatine,      p_io)
     CALL p_bcast(RFLD(ix)%imask,          p_io)
     CALL p_bcast(RFLD(ix)%repr,           p_io)
     CALL p_bcast(RFLD(ix)%cplcompID,      p_io)
     CALL p_bcast(RFLD(ix)%vertlev,        p_io)
     CALL p_bcast(RFLD(ix)%timelev,        p_io)
     CALL p_bcast(RFLD(ix)%numblev,        p_io)
     ppsize=0
     IF (p_parallel_io) THEN
        IF (ALLOCATED(RFLD(ix)%pp)) ppsize = SIZE(RFLD(ix)%pp)
     END IF
     CALL p_bcast(ppsize, p_io)
     IF (ppsize > 0) THEN
        IF (.NOT. p_parallel_io)  ALLOCATE(RFLD(ix)%pp(ppsize))
        DO iy = 1, ppsize
           CALL p_bcast(RFLD(ix)%pp(iy)%number,   p_io)
           CALL p_bcast(RFLD(ix)%pp(iy)%limit,    p_io)
           CALL p_bcast(RFLD(ix)%pp(iy)%action,   p_io)
        END DO
     END IF
  END DO

  CALL end_message_bi(modstr,'ANALYSE RFIELDS',substr)  ! log-output

    !------------------------------------------------------------------
  ! ANALYSE SFIELDs
  CALL start_message_bi(modstr,'ANALYSE SFIELDS',substr)  ! log-output
  IF (p_parallel_io) THEN
     cplsubid  = -99
     cpllevid  = -99
     cplcompID = -99
     oasisname = ' '
     numsfld = 0
     ! CHECK if OASISname exists in namcouple namelist.
     ! IF not, field does not have to be taken into account
     DO ioas = 1, ig_final_nfield
        lfound = .FALSE.
        CALL strcrack(namsrcfld(ioas),'_',tmp1,nl)

        DO ix = 1, MAXNUM_FIELDS

           IF ( TRIM(ADJUSTL(SFIELD(ix)%OASISname)) /= '') THEN
              IF (TRIM(ADJUSTL(SFIELD(ix)%oasisname)) == &
                   TRIM(tmp1(1)) .AND.           &
                   namsrcgrd(ioas)(1:cnlen) == TRIM(compname)  &
                   ) THEN

                 DO jx = 1,nummod
                    jxlen = LEN_TRIM(LOCCPL(jx)%modname)
                    IF (TRIM(LOCCPL(jx)%modname) == namdstgrd(ioas)(1:jxlen) &
                         ) THEN
                       IF (nl > 1) THEN
                          CALL str2num(tmp1(2),cpllevid(ioas), status)
                          IF (status /= 0) CALL error_bi( &
                               'ERROR in string conversion', substr)
                       END IF
                       LOCCPL(jx)%lactive = .TRUE.
                       oasisname(ioas) = namsrcfld(ioas)
                       cplcompid(ioas) = jx
                       cplsubid(ioas)  = ix
                       EXIT
                    END IF
                 END DO
                 IF (cplcompid(ioas) > 0) THEN
                    numsfld = numsfld +1
                    lfound = .TRUE.
                    EXIT
                 END IF
              END IF
            END IF
        END DO
        IF (ASSOCIATED(tmp1)) DEALLOCATE(tmp1)
        NULLIFY(tmp1)
        IF (.NOT. lfound .AND. (namsrcgrd(ioas)(1:cnlen) == TRIM(compname)))THEN
           write(iouerr,*) 'OASIS namelist field: ',namsrcfld(ioas)
           CALL error_bi('no matching entriy for source field found!', substr)
        END IF
     END DO

     ! -------------------------------------------------------------

     write (*,*)  ' -> NUMBER SEND FIELDS FOUND ',numsfld
     ALLOCATE(SFLD(numsfld))
     numsfld   = 0

     DO ioas = 1, ig_final_nfield


        IF ( cplcompid(ioas) < 0) CYCLE
        numsfld = numsfld + 1

        SFLD(numsfld)%OASISname     = TRIM(ADJUSTL(oasisname(ioas)))
        SFLD(numsfld)%MESSyname%CHA = &
             TRIM(ADJUSTL(SFIELD(cplsubid(ioas))%MESSyname%CHA))
        SFLD(numsfld)%MESSyname%OBJ = &
             TRIM(ADJUSTL(SFIELD(cplsubid(ioas))%MESSyname%OBJ))

        write (*,*) '--> ASSOCIATE OASIS SOURCE FIELD '  &
             , TRIM(SFLD(numsfld)%OASISname)
        write (*,*) '    to channel object '         &
             , TRIM(SFLD(numsfld)%messyname%cha), ' '&
             , TRIM(SFLD(numsfld)%messyname%obj)

        ! set cplcompID
        SFLD(numsfld)%cplcompID    = cplcompid(ioas)
        write (*,*) ' FIELD is sent to component ' &
             , TRIM(LOCCPL(SFLD(numsfld)%cplcompID)%modname)

        ! INVESTIGATE level string
        ! 1. is OASIS field multi-level source ?

        cpllev2: IF (cpllevid(ioas) > 0) THEN
           ! analyse how to interpret cpllevid (T,Z,N -AXIS?)

           IF (TRIM(ADJUSTL(SFIELD(cplsubid(ioas))%oasaxis)) /= '') THEN
              CALL strcrack(SFIELD(cplsubid(ioas))%oasaxis,';',tmp1,nl)

               DO il = 1, nl
                 axstr = '   '
                 CALL strcrack(tmp1(il),':', tmp2, ni)

                 IF (ni > 0 .OR. ni < 3) THEN
                    SELECT CASE (tmp2(1))
                    CASE('Z')
                       IF (SIZE(tmp2) == 1) zaxlen = nlev
                       axstr(il:il) = 'Z'
                    CASE('N')
                       IF (SIZE(tmp2) <= 1) CALL error_bi(&
                            'n-axis requires length specification', substr)
                       CALL str2num(tmp2(2), naxlen, status)
                       IF (status /= 0)  CALL error_bi( &
                            'error in string conversion! '//tmp2(2), substr)
                       axstr(il:il) = 'N'
                    CASE('T')
                       IF (SIZE(tmp2) <= 1) CALL error_bi(&
                            't-axis requires length specification', substr)
                       CALL str2num(tmp2(2), taxlen, status)
                       IF (status /= 0)  CALL error_bi( &
                            'error in string conversion! '//tmp2(2), substr)
                       axstr(il:il) = 'T'
                    CASE DEFAULT
                       CALL error_bi(' UNKNOWN AXIS IDENTIFIER! '//tmp2(1) &
                            , substr)
                    END SELECT
                 END IF
                IF (ASSOCIATED(tmp2)) DEALLOCATE(tmp2)
                NULLIFY(tmp2)
              END DO
              IF (ASSOCIATED(tmp1)) DEALLOCATE(tmp1)
              NULLIFY(tmp1)

              SELECT CASE(axstr)
              CASE ('Z  ')
                 IF (cpllevid(ioas) > zaxlen) THEN
                    write(iouerr,*) 'ERROR: level dimensions: ' &
                         ,cpllevid(ioas), zaxlen
                    CALL error_bi('level dimensions do not match', substr)
                 END IF
                 SFLD(numsfld)%vertlev=cpllevid(ioas)
              CASE ('T  ')
                 IF (cpllevid(ioas) > taxlen) THEN
                    write(iouerr,*) 'ERROR: t-axis length : ' &
                         ,cpllevid(ioas), taxlen
                    CALL error_bi('t-axis length do not match', substr)
                 END IF
                 SFLD(numsfld)%timelev=cpllevid(ioas)
              CASE ('N  ')
                 IF (cpllevid(ioas) > naxlen) THEN
                    write(iouerr,*) 'ERROR: t-axis length : ',cpllevid(ioas), naxlen
                    CALL error_bi('n-axis length do not match', substr)
                 END IF
                 SFLD(numsfld)%numblev=cpllevid(ioas)
              CASE ('ZT ')
                 SFLD(numsfld)%vertlev= INT(cpllevid(ioas) / zaxlen)
                 SFLD(numsfld)%timelev= MOD(cpllevid(ioas),zaxlen) + 1
                 IF (SFLD(numsfld)%timelev > taxlen) THEN
                    write(iouerr,*) 'ZT-axis dimensions: ' &
                         ,cpllevid(ioas), zaxlen, taxlen
                 END IF
              CASE ('ZN ')
                 SFLD(numsfld)%vertlev= INT(cpllevid(ioas) / zaxlen)
                 SFLD(numsfld)%numblev= MOD(cpllevid(ioas),zaxlen) + 1
                 IF (SFLD(numsfld)%numblev > naxlen) THEN
                    write(iouerr,*) 'ZN-axis dimensions: ' &
                         ,cpllevid(ioas), zaxlen, naxlen
                 END IF
              CASE ('TZ ')
                 SFLD(numsfld)%timelev= INT(cpllevid(ioas) / taxlen)
                 SFLD(numsfld)%vertlev= MOD(cpllevid(ioas),taxlen) + 1
                 IF (SFLD(numsfld)%vertlev > zaxlen) THEN
                    write(iouerr,*) 'TZ-axis dimensions: ' &
                         ,cpllevid(ioas), taxlen, zaxlen
                 END IF
              CASE ('NZ ')
                 SFLD(numsfld)%numblev= INT(cpllevid(ioas) /naxlen)
                 SFLD(numsfld)%vertlev= MOD(cpllevid(ioas),naxlen) + 1
                 IF (SFLD(numsfld)%vertlev > zaxlen) THEN
                    write(iouerr,*) 'NZ-axis dimensions: ' &
                         ,cpllevid(ioas), naxlen, zaxlen
                 END IF
              CASE ('TN ')
                 SFLD(numsfld)%timelev= INT(cpllevid(ioas) / taxlen)
                 SFLD(numsfld)%numblev= MOD(cpllevid(ioas),taxlen) + 1
                 IF (SFLD(numsfld)%numblev > naxlen) THEN
                    write(iouerr,*) 'TN-axis dimensions: ' &
                         ,cpllevid(ioas), taxlen, naxlen
                 END IF
              CASE ('NT ')
                 SFLD(numsfld)%numblev= INT(cpllevid(ioas) / naxlen)
                 SFLD(numsfld)%timelev= MOD(cpllevid(ioas),naxlen) + 1
                 IF (SFLD(numsfld)%timelev > taxlen) THEN
                    write(iouerr,*) 'NT-axis dimensions: ' &
                         ,cpllevid(ioas), naxlen, taxlen
                 END IF
              CASE DEFAULT
                 CALL error_bi('unknown OASIS-AXIS string! '//axstr, substr)
              END SELECT

           END IF

        END IF cpllev2


        CALL strcrack(TRIM(SFIELD(cplsubid(ioas))%levelstr), ';', tmp1, nl)
        nl2: IF (nl > 0) THEN
           IF (TRIM(tmp1(1)) /= '') THEN
              DO is = 1, nl
                 CALL strcrack(tmp1(is),':', tmp2, ni)
                 IF (ni /=2) CALL error_bi('levelstr is not correct', substr)
                 SELECT CASE (TRIM(tmp2(1)))
                 CASE ('Z')
                    CALL str2num(tmp2(2), il, status)
                    IF ( SFLD(numsfld)%vertlev < 0) THEN
                       IF (il == 0) THEN
                          SFLD(numsfld)%vertlev = nlev
                       ELSE
                          SFLD(numsfld)%vertlev = il
                       END IF
                    ELSE
                       IF (SFLD(numsfld)%vertlev /= il) THEN
                          write(iouerr,*) ' namelist settings inconsistent: ' &
                               , SFLD(numsfld)%vertlev, il
                          CALL error_bi(&
                          ' namelist settings inconsistent for vertical level' &
                             , substr)
                       END IF
                    END IF
                 CASE ('T')
                    CALL str2num(tmp2(2), il, status)
                    IF ( SFLD(numsfld)%timelev < 0) THEN
                       SFLD(numsfld)%timelev = il
                    ELSE
                       IF (SFLD(numsfld)%timelev /= il) THEN
                          write(iouerr,*) ' namelist settings inconsistent: ' &
                               , SFLD(numsfld)%timelev, il
                          CALL error_bi(&
                               ' namelist settings inconsistent for time level'&
                               , substr)
                       END IF
                    END IF
                 CASE ('N')
                    CALL str2num(tmp2(2), il, status)
                    IF ( SFLD(numsfld)%numblev < 0) THEN
                       SFLD(numsfld)%numblev = il
                    ELSE
                       IF (SFLD(numsfld)%numblev /= il) THEN
                          write(iouerr,*) ' namelist settings inconsistent: ' &
                               , SFLD(numsfld)%numblev, il
                          CALL error_bi(&
                          ' namelist settings inconsistent for number level' &
                          , substr)
                       END IF
                    END IF

                 CASE DEFAULT
                    CALL error_bi('unknown level specifier!' //TRIM(tmp2(1)) &
                         , substr)
                 END SELECT

                 IF (ASSOCIATED(tmp2)) DEALLOCATE(tmp2)
                 NULLIFY(tmp2)
              END DO
              IF (ASSOCIATED(tmp1)) DEALLOCATE(tmp1)
              NULLIFY(tmp1)
           END IF
        END IF nl2
        IF ( SFLD(numsfld)%vertlev > 0) &
             write (*,*) '         ... TO vertical level: ' &
             , SFLD(numsfld)%vertlev

        IF ( SFLD(numsfld)%timelev > -10) &
             write (*,*) '         ...  TO time level: '  &
             , SFLD(numsfld)%timelev

        IF ( SFLD(numsfld)%numblev > 0) &
             write (*,*) '        ...  TO number level: ' &
             , SFLD(numsfld)%numblev

        IF ( SFLD(numsfld)%vertlev > 0 .AND. SFLD(numsfld)%timelev > -10 &
             .AND. SFLD(numsfld)%numblev > 0) CALL error_bi( &
             ' SFLD: too many rank reductions requested ', substr)

     END DO
     write (*,*) 'NUMBER OF SOURCE FIELDS FOUND ', numsfld
  END IF

  ! BROADCAST SFIELD from I/O-PE to ALL OTHER PEs
  CALL p_bcast(numsfld, p_io)
  IF (.NOT. p_parallel_io) ALLOCATE(SFLD(numsfld))
  DO ix = 1, numsfld
     CALL p_bcast(SFLD(ix)%OASISname,      p_io)
     CALL p_bcast(SFLD(ix)%MESSyname%cha,  p_io)
     CALL p_bcast(SFLD(ix)%MESSyname%obj,  p_io)
     CALL p_bcast(SFLD(ix)%cplcompID,      p_io)
     CALL p_bcast(SFLD(ix)%vertlev,        p_io)
     CALL p_bcast(SFLD(ix)%timelev,        p_io)
     CALL p_bcast(SFLD(ix)%numblev,        p_io)
  END DO
  CALL end_message_bi(modstr,'ANALYSE SFIELDS',substr)  ! log-output

   ! ----------------------------------------------------------
  ! BROADCAST lactive (set during RFIELD AND SFIELD analyses)

    DO ix = 1, nummod
       CALL p_bcast(LOCCPL(ix)%lactive,  p_io)
    END DO

   ! ----------------------------------------------------------

    IF (L_IS_CHILD) THEN
       ! READ CPL_EXTDATA namelist
       IF (p_parallel_io) THEN                  ! read only on I/O-PE
          iou = find_next_free_unit(100,200)    ! find next free I/O unit
          CALL oasis3mct_read_nml_cpl_extdata(status, iou)  ! read namelist
          ! terminate if error
          IF (status /= 0) &
               CALL error_bi('Error in reading CPL_COSMO namelist',substr)
       END IF
       CALL p_bcast(slf_file,  p_io)
       CALL p_bcast(slf_name,  p_io)
       CALL p_bcast(lon_name,  p_io)
       CALL p_bcast(lat_name,  p_io)
    END IF

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output
#endif
  END SUBROUTINE oasis3mct_initialize
  ! ====================================================================

 SUBROUTINE oasis3mct_init_memory

#if !defined(COSMO) && !defined(ECHAM5)
   USE messy_main_blather_bi, ONLY: error_bi
#endif

   IMPLICIT NONE

   CHARACTER(LEN=*), PARAMETER :: substr = 'oasis3mct_init_memory'

#ifdef COSMO
   CALL oasis3mct_write_grids_cosmo
#elif ECHAM5
   CALL oasis3mct_write_grids_echam
#else
    CALL error_bi('coupling indices not implemented for basemodel', substr)
#endif

 END SUBROUTINE oasis3mct_init_memory

  ! ====================================================================
  SUBROUTINE oasis3mct_init_coupler(flag)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: flag

    SELECT CASE(flag)
    CASE (1)
       CALL oasis3mct_init_memory
    CASE (2)
       CALL oasis3mct_init_coupling
    END SELECT

  END SUBROUTINE oasis3mct_init_coupler
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE oasis3mct_init_coupling

#ifdef OASIS3MCT
    ! MESSy/BMIL
    USE  MESSY_MAIN_CHANNEL_ERROR_BI, ONLY: channel_halt
    USE  MESSY_MAIN_CHANNEL_BI,       ONLY: GP_2D_HORIZONTAL
    USE  MESSY_MAIN_MPI_BI,           ONLY: p_parallel_io
#ifdef ECHAM5
    USE messy_main_mpi_bi,            ONLY: dcl
#endif
    ! MESSy/SMCL
    USE messy_main_channel_repr,  ONLY: get_representation_info
    USE MESSY_MAIN_CHANNEL,       ONLY: get_channel_object       &
                                      , new_channel_object       &
                                      , get_channel_object_info  &
                                      , get_channel_object_slice &
                                      , new_channel
#endif

    IMPLICIT NONE

#ifdef OASIS3MCT
    CHARACTER(LEN=*), PARAMETER :: substr = 'oasis3mct_init_coupling'
    INTEGER                     :: ix
    INTEGER                     :: status
    INTEGER                     :: reprid
    INTEGER                     :: var_nodims(2)
    INTEGER                     :: var_shape(4)
    INTEGER                     :: cplfreq(1)

   CALL start_message_bi(modstr,'  ',substr)

   CALL new_channel(status, modstr)
   CALL channel_halt(substr, status)

    ! ==============================================================
    ! GET POINTERS TO SEND FIELDS
    ! ==============================================================
    DO ix = 1, NUMSFLD
       write(iouerr,*) substr, ' SFIELD ',TRIM(SFLD(ix)%messyname%CHA) &
            ,' ',TRIM(SFLD(ix)%messyname%OBJ)

       ! GET ADDITIONAL REPRESENTATION INFO AND DATA ACCESS FOR AUTOMATIC
       ! FIELD CONTROL IN SEND / RECEIVE
       ! DEFINE AXIS STRING
       ! 1. get representation id
       CALL get_channel_object_info(status                              &
            , TRIM(SFLD(ix)%messyname%CHA),TRIM(SFLD(ix)%messyname%OBJ) &
            , reprid=reprid)
       CALL channel_halt(substr,status)
       ! 2. get axis string
       CALL get_representation_info(status, ' ', reprid, rank=SFLD(ix)%rank)
       CALL channel_halt(substr, status)

       IF (SFLD(ix)%vertlev < 0 .AND. SFLD(ix)%timelev < -10 .AND. &
            SFLD(ix)%numblev < 0) THEN
          ! NO REDUCTION OF LEVEL
          IF (SFLD(ix)%rank == 2) THEN
             CALL get_channel_object(status &
                  , TRIM(SFLD(ix)%messyname%cha), TRIM(SFLD(ix)%messyname%obj) &
                  , p2=SFLD(ix)%ptr2d)
             IF (status /= 0) &
                  write(iouerr,*) substr, TRIM(SFLD(ix)%messyname%cha),' ' &
                  , TRIM(SFLD(ix)%messyname%obj)
             CALL channel_halt(substr,status)
             IF (p_parallel_io) write (*,*)   TRIM(SFLD(ix)%OASISNAME),':'&
             ,TRIM(SFLD(ix)%messyname%cha)//' '//TRIM(SFLD(ix)%messyname%obj) &
             ,': got 2d pointer for 2d field'
          ELSE
             write(iouerr,*) ' FIELD: ', TRIM(SFLD(ix)%OASISNAME) &
                  , SFLD(ix)%rank
             write(iouerr,*) ' Level indicators:', SFLD(ix)%vertlev &
                  , SFLD(ix)%timelev , SFLD(ix)%numblev
             CALL error_bi('FIELD has too many ranks 1! ',substr)
          END IF
       ELSE IF (SFLD(ix)%rank == 2) THEN
          CALL error_bi('rank reduction not possible for 2D Field! ', substr)
       ELSE IF (SFLD(ix)%timelev < -10) THEN
          ! NO TIME LEVEL REDUCTION
          CALL get_channel_object_slice(status &
                  , TRIM(SFLD(ix)%messyname%cha), TRIM(SFLD(ix)%messyname%obj) &
                  , p2=SFLD(ix)%ptr2d &
                  , zslice=SFLD(ix)%vertlev, nslice=SFLD(ix)%numblev)
             IF (status /= 0)  write(iouerr,*) substr                  &
                  , TRIM(SFLD(ix)%messyname%cha),' '               &
                  , TRIM(SFLD(ix)%messyname%obj), SFLD(ix)%vertlev &
                  , SFLD(ix)%numblev
          CALL channel_halt(substr,status)
             IF (p_parallel_io) write (*,*)   TRIM(SFLD(ix)%OASISNAME),':'&
             ,TRIM(SFLD(ix)%messyname%cha)//' '//TRIM(SFLD(ix)%messyname%obj) &
             ,': got 2d pointer for 3d/4d field: z=',SFLD(ix)%vertlev,'n=' &
             , SFLD(ix)%numblev
       ELSE
          ! => timelevel reduction
          ! 3d POINTER REQUIRED as time level reduction in COSMO only possible
          ! during runtime
           CALL get_channel_object_slice(status &
                  , TRIM(SFLD(ix)%messyname%cha), TRIM(SFLD(ix)%messyname%obj) &
                  , p3=SFLD(ix)%ptr3d &
                  , zslice=SFLD(ix)%vertlev, nslice=SFLD(ix)%numblev)
             IF (status /= 0)  write(iouerr,*) substr                   &
                  , TRIM(SFLD(ix)%messyname%cha),' '                &
                  , TRIM(SFLD(ix)%messyname%obj), SFLD(ix)%vertlev  &
                  , SFLD(ix)%numblev
          CALL channel_halt(substr,status)
             IF (p_parallel_io) write (*,*)   TRIM(SFLD(ix)%OASISNAME),':'   &
            ,TRIM(SFLD(ix)%messyname%cha)//' '//TRIM(SFLD(ix)%messyname%obj) &
            ,': got 3d pointer for time level reduction of field during runtime'&
            , SFLD(ix)%timelev
       END IF

    END DO
    ! ==============================================================
    ! GET POINTERS TO RECEIVE FIELDS
    ! MAKE MEMORY FOR OASIS TRACERS
    ! ==============================================================

    DO ix = 1, NUMRFLD
       ! CREATE CHANNEL OBJECT FOR CONTINUOUS APPLICATION OF EXCHANGE FIELD
       ! As OASIS fields are always 2D use here GP_2D_HORIZONTAL
       CALL get_channel_object(status, modstr, TRIM(RFLD(ix)%messyname%obj)    &
                  , p2=RFLD(ix)%store)
       IF (status /= 0) THEN
          CALL new_channel_object(status, modstr, TRIM(RFLD(ix)%messyname%obj) &
            , p2 = RFLD(ix)%store, reprid = GP_2D_HORIZONTAL, lrestreq=.TRUE.)
          CALL channel_halt(substr,status)
       END IF

       CALL get_channel_object(status                                  &
           , TRIM(RFLD(ix)%messyname%CHA),TRIM(RFLD(ix)%messyname%OBJ))
       if_stat: IF (status /= 0) THEN
          if_create: IF (.NOT. RFLD(ix)%lcreatine) THEN
            write(iouerr,*) 'ERROR: RFLD CHANNEL OBJECT '&
                  , TRIM(RFLD(ix)%messyname%cha)&
                  , TRIM(RFLD(ix)%messyname%obj),' does not exist'
             CALL channel_halt(substr,status)
          END IF if_create
       END IF if_stat

       IF (p_parallel_io) write (*,*)  &
            ' get representation info for channel object: ' &
            , TRIM(RFLD(ix)%messyname%CHA),' ',TRIM(RFLD(ix)%messyname%OBJ)
       ! GET ADDITIONAL REPRESENTATION INFO AND DATA ACCESS FOR AUTOMATIC
       ! FIELD CONTROL IN SEND / RECEIVE
       ! DEFINE AXIS STRING
       ! 1. get representation id
       CALL get_channel_object_info(status                              &
            , TRIM(RFLD(ix)%messyname%CHA),TRIM(RFLD(ix)%messyname%OBJ) &
            , reprid=reprid)
       CALL channel_halt(substr,status)
       ! 2. get axis string
       CALL get_representation_info(status, ' ', reprid, rank=RFLD(ix)%rank)
       CALL channel_halt(substr, status)


       IF (RFLD(ix)%vertlev < 0 .AND. RFLD(ix)%timelev < -10 .AND. &
            RFLD(ix)%numblev < 0) THEN
          ! NO REDUCTION OF LEVEL
          IF (RFLD(ix)%rank == 2) THEN
             CALL get_channel_object(status &
                  , TRIM(RFLD(ix)%messyname%cha), TRIM(RFLD(ix)%messyname%obj) &
                  , p2=RFLD(ix)%ptr2d)
             CALL channel_halt(substr,status)
          ELSE
             write(iouerr,*) 'FIELD: ',  TRIM(RFLD(ix)%messyname%cha) &
                  , TRIM(RFLD(ix)%messyname%obj)
             CALL error_bi('FIELD has too many ranks 2 ! ',substr)
          END IF
       ELSE IF (RFLD(ix)%rank == 2) THEN
          CALL error_bi('rank reduction not possible for 2D Field! ', substr)
       ELSE IF (RFLD(ix)%timelev < -10) THEN
          ! NO TIME LEVEL REDUCTION
          CALL get_channel_object_slice(status &
                  , TRIM(RFLD(ix)%messyname%cha), TRIM(RFLD(ix)%messyname%obj) &
                  , p2=RFLD(ix)%ptr2d &
                  , zslice=RFLD(ix)%vertlev, nslice=RFLD(ix)%numblev)
          CALL channel_halt(substr,status)
       ELSE
          ! => timelevel reduction
          ! 3d POINTER REQUIRED as time level reduction in COSMO only possible
          ! during runtime
          IF (RFLD(ix)%vertlev > 0 .OR. RFLD(ix)%numblev > 0 ) THEN
             CALL get_channel_object_slice(status &
                  , TRIM(RFLD(ix)%messyname%cha), TRIM(RFLD(ix)%messyname%obj) &
                  , p3=RFLD(ix)%ptr3d &
                  , zslice=RFLD(ix)%vertlev, nslice=RFLD(ix)%numblev)
             CALL channel_halt(substr,status)
          ELSE
             CALL get_channel_object(status &
                  , TRIM(RFLD(ix)%messyname%cha), TRIM(RFLD(ix)%messyname%obj) &
                  , p3=RFLD(ix)%ptr3d)
             CALL channel_halt(substr,status)
          END IF
       END IF


    END DO



    ! ==============================================================
    !
    ! DEFINE OASIS VARIABLES
    !
    ! ==============================================================

    ! Define oasis dimension variables
    var_nodims(1)   = 2 ! number of ranks, 1 or 2
    var_nodims(2)   = 1 ! number of bundles, in OASIS3-MCT always 1
    var_shape(1)    = 1 ! minimum index , in OASIS3-MCT always 1
    var_shape(3)    = 1 ! minimum index , in OASIS3-MCT always 1
#ifdef COSMO
    var_shape(2)    = paral(3) ! maximum index (dimension length)
    var_shape(4)    = paral(4) ! maximum index (dimension length)
#elif ECHAM5
    var_shape(2)    = dcl%nglon ! maximum index (dimension length)
    var_shape(4)    = dcl%nglat ! maximum index (dimension length)
#endif
    write(iouerr,*) ' DEFINE SOURCE VARIABLES!'
    ! A) SEND FIELDS
    DO ix = 1, numsfld
       write(iouerr,*) ' DEFINE SOURCE VARIABLE: ', ix, TRIM(SFLD(ix)%oasisname)
       CALL oasis_def_var( SFLD(ix)%oasvarid, TRIM(SFLD(ix)%oasisname) &
            , part_id, var_nodims, OASIS_Out, var_shape      &
            , OASIS_DOUBLE, status)
       IF (status /= Oasis_ok) &
            CALL error_bi('ERROR in OASIS_DEF_VAR for SFLD',substr)

    END DO

    ! B) RECEIVE FIELDS
    DO ix = 1, numrfld
       write(iouerr,*) ' DEFINE RECEIVE VARIABLE: ', ix, TRIM(RFLD(ix)%oasisname)
       CALL oasis_def_var( RFLD(ix)%oasvarid, TRIM(RFLD(ix)%oasisname) &
            , part_id, var_nodims, OASIS_In, var_shape      &
            , OASIS_DOUBLE, status)
       IF (status /= Oasis_ok) &
            CALL error_bi('ERROR in OASIS_DEF_VAR for RFLD',substr)

    END DO

    !CALL oasis_set_debug(20) !KEEP useful for new coupling implementations
    CALL oasis_enddef (status)
    IF (status /= Oasis_ok) CALL error_bi('ERROR in OASIS_ENDDEF',substr)
! --------------------------------------------

      DO ix = 1, numrfld

       ! GET COUPLING FREQUENCY FOR COmmunity Land Model
       ! (ALBEDO coupled in both coupling schemes)
       IF ( (TRIM(RFLD(ix)%oasisname) == 'COSALBED') .AND. &
           (LOCCPL(RFLD(ix)%cplcompID)%modname(1:3) == 'clm')) THEN

          CALL oasis_get_freqs(RFLD(ix)%oasvarID, OASIS_in,1 &
               ,cplfreq, status)

          isoil_external_cplfreq = cplfreq(1)
          write(iouerr,*) 'SOIL_EXT CPLFREQ ', isoil_external_cplfreq
       END IF

    END DO


    CALL end_message_bi(modstr,'  ',substr)
#endif

  END SUBROUTINE oasis3mct_init_coupling
 ! ====================================================================

  SUBROUTINE oasis3mct_init_loop(callID,flag)

    ! currently called by COSMO

    INTEGER, INTENT(IN) :: callID
    INTEGER, INTENT(IN) :: flag

    CALL oasis3mct_sndrcv(callID,flag)

  END SUBROUTINE oasis3mct_init_loop

! ====================================================================

  SUBROUTINE oasis3mct_time(callID,flag)

    ! currently called by ECHAM5 for FESOM coupling

    INTEGER, INTENT(IN) :: callID
    INTEGER, INTENT(IN) :: flag

    CALL oasis3mct_sndrcv(callID,flag)

  END SUBROUTINE oasis3mct_time

 ! ====================================================================

  SUBROUTINE oasis3mct_global_start(callID,flag)

    ! currently called by ECHAM5

    INTEGER, INTENT(IN) :: callID
    INTEGER, INTENT(IN) :: flag

    CALL oasis3mct_sndrcv(callID,flag)

  END SUBROUTINE oasis3mct_global_start

 ! ====================================================================

  SUBROUTINE oasis3mct_convec(callID,flag)

    ! currently called by COSMO

    INTEGER, INTENT(IN) :: callID
    INTEGER, INTENT(IN) :: flag

    CALL oasis3mct_sndrcv(callID,flag)

  END SUBROUTINE oasis3mct_convec

 ! ====================================================================
 ! ====================================================================

  SUBROUTINE oasis3mct_write_output(callID,flag)

    ! currently called from ECHAM

    INTEGER, INTENT(IN) :: callID
    INTEGER, INTENT(IN) :: flag

    CALL oasis3mct_sndrcv(callID,flag)

  END SUBROUTINE oasis3mct_write_output

 ! ====================================================================

  SUBROUTINE oasis3mct_finalize(lcoupled)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lcoupled


#ifdef OASIS3MCT
    CHARACTER(LEN=*), PARAMETER :: substr='oasis3mct_finalize'
    INTEGER  :: ix

    CALL start_message_bi(modstr,'TERMINATE OASIS',substr)  ! log-output
    CALL oasis_terminate !(status)
    CALL end_message_bi(modstr,'TERMINATE OASIS',substr)  ! log-output

    IF (.NOT. lcoupled) RETURN

    CALL start_message_bi(modstr,'DEALLOC SUBMODEL MEMORY',substr)  ! log-output

    IF (ASSOCIATED(llandmask)) DEALLOCATE(llandmask)
    DO ix = 1, SIZE(LOCCPL)
       IF (ASSOCIATED(LOCCPL(ix)%mask)) DEALLOCATE(LOCCPL(ix)%mask)
    END DO
    IF (ALLOCATED(LOCCPL)) DEALLOCATE(LOCCPL)
    DO ix = 1, SIZE(RFLD)
       IF (ALLOCATED(RFLD(ix)%pp)) DEALLOCATE (RFLD(ix)%pp)
    END DO
    IF (ALLOCATED(RFLD)) DEALLOCATE(RFLD)
    IF (ALLOCATED(SFLD)) DEALLOCATE(SFLD)

    DEALLOCATE(paral)

    CALL end_message_bi(modstr,'DEALLOC SUBMODEL MEMORY',substr)  ! log-output
#endif

  END SUBROUTINE oasis3mct_finalize

 ! ====================================================================



  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################


  SUBROUTINE oasis3mct_sndrcv(callID,flag)

#ifdef OASIS3MCT
    USE messy_main_timer, ONLY: time_span_s &
                              , YEAR_START, MONTH_START, DAY_START, HOUR_START &
                              , MINUTE_START, SECOND_START, MILLISECOND_START  &
                              , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND         &
                              , MILLISECOND
    USE messy_main_control,  ONLY: MEPS_ONE

#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: callID
    INTEGER, INTENT(IN) :: flag

#ifdef OASIS3MCT
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'oasis3mct_sndrcv'
    INTEGER                     :: ix
    INTEGER                     :: date
    REAL(dp)                    :: timespan

    ! calculate date in seconds
    ! TODO: make unit chosable seconds, minutes, hours, days
    CALL time_span_s(timespan &
         , YEAR_START, MONTH_START, DAY_START, HOUR_START  &
         , MINUTE_START, SECOND_START, MILLISECOND_START   &
         , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, MILLISECOND )
    date = NINT(timespan)

    ! 1    FIRST CYCLE
    ! 1.1  CHECK SEND FIELDS (order = 1)
       ! 1.1.1 CHECK IF FIELD SHOULD BE COUPLED BY CALLING SUBROUTINE
       !       at correct location (Begin /end)
    DO ix = 1, numsfld
       ! write(iouerr,*) 'S1 ',TRIM(SFLD(ix)%messyname%obj) &
       !,' CALL ', CALLID, LOCCPL(SFLD(ix)%cplcompID)%locsndSR &
       !,  'FLAG ',flag,LOCCPL(SFLD(ix)%cplcompID)%locsndBE    &
       !, LOCCPL(SFLD(ix)%cplcompID)%locsndord, SFLD(ix)%cplcompID
       IF (callID /= LOCCPL(SFLD(ix)%cplcompID)%locsndSR) CYCLE
       IF (flag /= LOCCPL(SFLD(ix)%cplcompID)%locsndBE .AND. flag /= MEPS_ONE)&
            CYCLE
       ! 1.1.2 CHECK IF order = 1 (called before receive, if not => second call)
       IF (LOCCPL(SFLD(ix)%cplcompID)%locsndord /= 1)     CYCLE
       ! 1.1.3 CALL OASIS send routine
       CALL oasis_send(SFLD(ix))
    END DO

    ! 1.2  CHECK RECEIVE FIELDS (order = 1 or 2)
    DO ix = 1, numrfld
       ! 1.2.1 CHECK IF FIELD SHOULD BE COUPLED BY CALLING SUBROUTINE
       !       at correct location (Begin /end)
       IF (callID /= LOCCPL(RFLD(ix)%cplcompID)%locrcvSR) CYCLE
       IF (flag /= LOCCPL(RFLD(ix)%cplcompID)%locrcvBE .AND. flag /= MEPS_ONE) CYCLE
       ! 1.2.2 CALL OASIS receive routine

       ! RECEIVE data from OASIS and postprocess field according to namelist
       CALL oasis_rcv(RFLD(ix))

    END DO

    ! 1.3  CHECK SEND FIELDS (order = 2)
    DO ix = 1, numsfld
       ! 1.3.1 CHECK IF FIELD SHOULD BE COUPLED BY CALLING SUBROUTINE
       !       at correct location (Begin /end)
       IF (callID /= LOCCPL(SFLD(ix)%cplcompID)%locsndSR) CYCLE
       IF (flag /= LOCCPL(SFLD(ix)%cplcompID)%locsndBE .AND. flag /= MEPS_ONE) CYCLE
       ! 1.3.2 CHECK IF order = 2 (called after receive)
       IF (LOCCPL(SFLD(ix)%cplcompID)%locsndord /= 2)     CYCLE

       ! 1.3.3 CALL OASIS send routine
       CALL oasis_send(SFLD(ix))

    END DO

    ! MODEL and coupling SPECIFIC CALCULATIONS after coupling
#ifdef COSMO
    ! For COSMO-CLM coupling:

   IF (lcpl_land) THEN
      IF (isoil_external_ppSR == callID .AND. &
           ( isoil_external_ppBE == flag .OR. flag == MEPS_ONE)) THEN
            CALL postproc_cosmo_clm
      END IF
   END IF
#endif
#ifdef ECHAM
      CALL prostproc_echam_fesom
#endif
    RETURN

  CONTAINS

    SUBROUTINE oasis_send(SFLD)

#ifdef COSMO
      USE messy_main_data_bi, ONLY: nnow, nnew
#endif

      IMPLICIT NONE

      ! IN/OUT
      TYPE(t_cpl_sfield) :: SFLD
      ! LOCAL
      CHARACTER(LEN=*), PARAMETER       :: substr = 'oasis_send'
      INTEGER                           :: info
      REAL(dp), DIMENSION(:,:), POINTER :: ptr2d => NULL()
      REAL(dp), DIMENSION(:,:), POINTER :: xptr2d => NULL()
#ifdef ECHAM5
      INTEGER                           :: iseg, jseg, idx, iylen, iy,jy
#endif

!      write(iouerr,*) 'OASIS-SEND:01 ', TRIM(SFLD%oasisname)
!           , SFLD%timelev, date, SFLD%oasvarID
      IF (ASSOCIATED(SFLD%ptr2d)) THEN
         ptr2d => SFLD%ptr2d(iis:iie,jjs:jje)
      ELSE
         IF (SFLD%timelev >0) THEN
            ptr2d => SFLD%ptr3d(iis:iie,jjs:jje,SFLD%timelev)
#ifdef COSMO
         ELSE IF (SFLD%timelev == -1) THEN
            ptr2d => SFLD%ptr3d(iis:iie,jjs:jje,nnow)
         ELSE IF (SFLD%timelev == -2) THEN
            ptr2d => SFLD%ptr3d(iis:iie,jjs:jje,nnew)
#endif
         ELSE
            CALL error_bi(' INVALID TIMELEVEL', substr)
         END IF
      END IF

#ifndef ECHAM5
      xptr2d => ptr2d
#else
      ALLOCATE(xptr2d(nsegslen,nsegs))
      iylen = iie - iis + 1
      DO jy = jjs, jje
         DO iy = iis , iie
            idx = (jy -1) * iylen + iy
            jseg = INT( (idx-1) / nsegslen) + 1
            iseg = idx - nsegslen * (jseg-1)
            xptr2d(iseg,jseg) = ptr2d(iy,jy)
         END DO
      END DO
#endif
!      write (*,*) 'KKK OASIS3MCT SEND: date [s]: ' &
!           , date,'FLD: ',TRIM(SFLD%oasisname)
      CALL oasis_put (SFLD%oasvarID, date, xptr2d, info)

      ! CLEAN UP
      NULLIFY(ptr2d)
#ifdef ECHAM5
      DEALLOCATE(xptr2d)
#endif
      NULLIFY(xptr2d)

    END SUBROUTINE oasis_send

    ! ----------------------------------------------------

    SUBROUTINE oasis_rcv(RFLD)

#ifdef COSMO
      USE messy_main_data_bi, ONLY: nnow, nnew
#endif

      IMPLICIT NONE

      ! IN/OUT
      TYPE(t_cpl_rfield) :: RFLD

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER       :: substr = 'oasis_rcv'
      INTEGER                           :: info
      REAL(dp), DIMENSION(:,:), POINTER :: rcvptr  => NULL()
      INTEGER                           :: iy, jy

#ifdef ECHAM5
      REAL(dp), DIMENSION(:,:), POINTER :: rcvptr2 => NULL()
      INTEGER                           :: iseg, jseg, idx, iylen
#endif
!      write(iouerr,*) 'KKK OASIS-RECV:00 ', TRIM(RFLD%oasisname),' date [s] ',date

      ! currently only 2D fields (always xy are exchanged via oasis)

      ALLOCATE(rcvptr(iis:iie,jjs:jje))

#ifndef ECHAM5
      CALL oasis_get (RFLD%oasvarID, date, rcvptr, info)
#else
      ALLOCATE(rcvptr2(nsegslen,nsegs))
      CALL oasis_get (RFLD%oasvarID, date, rcvptr2, info)
      iylen = iie-iis+1 ! should be nproma for echam
      DO jy = jjs,jje
         DO iy = iis,iie
            idx  = (jy-1)*iylen + iy
            jseg = INT( (idx -1) / nsegslen) + 1
            iseg = idx - (jseg-1) * nsegslen
            rcvptr(iy,jy) = rcvptr2(iseg, jseg)
         END DO
      END DO
      DEALLOCATE(rcvptr2)
#endif

      IF (info == OASIS_Recvd) THEN

         CALL postprocess_field(RFLD, rcvptr)

         DO iy = iis,iie
            DO jy = jjs,jje
               RFLD%store(iy,jy) = rcvptr(iy,jy)
            END DO
         END DO

      END IF
!         write(iouerr,*) 'OASIS-RECEIVED BEFORE PP: ',  TRIM(RFLD%oasisname) &
!              , MINVAL(rcvptr), MAXVAL(rcvptr)


      IF (ASSOCIATED(RFLD%ptr2d)) THEN
         DO iy = iis,iie
            DO jy = jjs,jje
               IF (  LOCCPL(RFLD%CplCompID)%mask(iy,jy) ) THEN
                  RFLD%ptr2d(iy,jy) = RFLD%store(iy,jy)
               END IF
            END DO
         END DO
      ELSE
         IF (RFLD%timelev > 0) THEN
            DO iy = iis,iie
               DO jy = jjs,jje
                  IF (  LOCCPL(RFLD%CplCompID)%mask(iy,jy) ) THEN
                     RFLD%ptr3d(iy,jy,RFLD%timelev) = RFLD%store(iy,jy)
                  END IF
               END DO
            END DO
#ifdef COSMO
         ELSE IF (RFLD%timelev == -1) THEN
            DO iy = iis,iie
               DO jy = jjs,jje
                  IF (  LOCCPL(RFLD%CplCompID)%mask(iy,jy) ) THEN
                     RFLD%ptr3d(iy,jy,nnow) = RFLD%store(iy,jy)
                  END IF
               END DO
            END DO
         ELSE IF (RFLD%timelev == -2) THEN
            DO iy = iis,iie
               DO jy = jjs,jje
                  IF (  LOCCPL(RFLD%CplCompID)%mask(iy,jy) ) THEN
                     RFLD%ptr3d(iy,jy,nnew) = RFLD%store(iy,jy)
                  END IF
               END DO
            END DO
#endif
         ELSE
            CALL error_bi(' UNVALID TIMELEVEL', substr)
         END IF
      END IF

      DEALLOCATE(rcvptr)
      NULLIFY(rcvptr)

    END SUBROUTINE oasis_rcv
    ! ----------------------------------------------------


    SUBROUTINE postproc_cosmo_clm

#if defined(COSMO) && !defined(COSMOv509)
      ! BMIL
      USE messy_main_data_bi,         ONLY: nnow, nnew
      USE messy_main_grid_def_mem_bi, ONLY: ke, ie, je
      USE messy_main_grid_def_bi,     ONLY: hhl
      USE messy_main_data_bi,         ONLY:   &
#ifndef COSMOv5s5
                                             t_g, t_s, qv_s                 &
#else
                                             t_g, t_s, qv_s                 &
#endif
                                           , p0, pp, ps                &
                                           , u, v, t, tch, tcm              &
                                           , shfl_s, lhfl_s, umfl_s, vmfl_s &
                                           , idt_qv, vel_min
      USE messy_main_tracer_mem_bi,     ONLY: xtm1

      !  SMCL
      USE messy_main_timer,         ONLY: dt => time_step_len
      USE messy_main_constants_mem, ONLY: g, r_d => rd, cp_d=> cp_air &
                                        , lh_v => alv, rvd_m_o=> vtmpc1


      IMPLICIT NONE


      INTEGER :: nx
      INTEGER :: i, j, im1 ,jm1
      REAL(dp) :: zpp, dzke, zvbke, ztvb, coef_lim
      REAL(dp), DIMENSION(ie,je) :: zpianf
      REAL(dp), DIMENSION(ie,je) :: zpia
      REAL(dp)                   :: tcm_epsi
      REAL(dp)                   :: ztch
      REAL(dp), PARAMETER        :: gr    = 1._dp / g
      REAL(dp), PARAMETER        :: rdocp = r_d / cp_d

      ! Statement function zfpi for calculation of the exner function
      ! where the dummy argument zppa is pressure

      nx = nnow

      DO j = jjs, jje
         jm1 = MAX(j-1,1)
         DO i = iis, iie
            im1 = MAX(i-1,1)

            ifland: IF (llandmask(i,j)) THEN

            t_g(i,j,nx)   = t_g(i,j,nnew)           !Surface temperture
            t_s(i,j,nx)   = t_g(i,j,nnew)           !Surface temperature
            t_s(i,j,nnew) = t_g(i,j,nnew)

            !*******************************************************************
            !***  Compute COSMO transfer coefficients              ************
            !    itype_turb = 3 , imode_turb = 1
            !   (Neumann boundary conditions for heat and moisture transport  *
            !      at the lower boundary (speciffied fluxes)       ************
            !*******************************************************************


            zpp         = p0(i,j,ke) + pp(i,j,ke,nx)
            dzke        = hhl(i,j,ke) - hhl(i,j,ke+1)
            zpia(i,j)   = zfpi( zpp, rdocp)
            zpianf(i,j) = zfpi( ps(i,j,nx), rdocp)
            zvbke       = MAX( 0.5_dp *SQRT((u(i,j,ke,nx) + u(im1,j,ke,nx))**2 &
                                +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2), vel_min )
            ztvb        = t_s (i,j,nx)*(1.0_dp + rvd_m_o*qv_s(i,j,nx))

            ! Sensible heat flux inversion: derive new tch from clm flux
            tch(i,j) = -shfl_s (i,j) / (                          &
                 (zvbke*ps(i,j,nx)/(r_d*ztvb))*cp_d*              &
                    ( t_g(i,j,nx) - zpianf(i,j)*t(i,j,ke,nx)/zpia(i,j) ) )


            !-------------------------------------------------------------------
            ! calculate effective surface humidity based on CLM latent heat flux
            !-------------------------------------------------------------------
            ztch  = tch(i,j)*zvbke*ps(i,j,nx)/(r_d*ztvb)
            IF (ABS(ztch) >= 1.e-20_dp) THEN
               qv_s(i,j,nx) = &
                    xtm1(_RI_XYZN_(i,j,ke,idt_qv)) - lhfl_s(i,j) / (lh_v * ztch)
            ELSE
               qv_s(i,j,nx) = xtm1(_RI_XYZN_(i,j,ke,idt_qv))
            END IF
            ! limit qv_s to positive values (no negative surface humidity ;-)
            qv_s(i,j,nx) = min(max(1.E-5_dp,qv_s(i,j,nx)),1.E-1_dp)

            ! set qv_s(nnew)=qv_s(nnow) similarly to what is done in TERRA this
            ! will ensure that the flux rediagnosed in COSMO is the same as here
            qv_s(i,j,nnew) = qv_s(i,j,nx)


            !----------------------------------------------------------------
            ! momentum flux inversion: derive new tcm from CLM momentum flux
            !----------------------------------------------------------------
            !
            ! correct momentum fluxes for wind direction:
            ! u/v winds are always positive in CLM3.5??
            ! momentum flx are always negative in CLM4?
             if (u(i,j,ke,nx)<0.) umfl_s(i,j) = -umfl_s(i,j)
             if (v(i,j,ke,nx)<0.) vmfl_s(i,j) = -vmfl_s(i,j)

             ! umfl , vmfl from CLM are in mass points and
             ! tcm, tch is also in mass point
             if (u(i,j,ke,nx) == 0) then
               tcm(i,j) = 0._dp !stable condition
            ELSE
               tcm(i,j) = umfl_s(i,j) / &
                     (zvbke*g*ps(i,j,nx)/(r_d*ztvb) * gr* u(i,j,ke,nx))
            ENDIF


            ! implement limiter so exchange coeff. is not
            !     too large (Oliver Fuhrer, 2007)
             ! NOTE: might add a check if tcm<0 here!
             tcm_epsi=0.5_dp
             coef_lim=tcm_epsi*dzke/(zvbke*dt)
             if (abs(tcm(i,j))>coef_lim) then
                tcm(i,j)=sign(coef_lim,tcm(i,j))
             end if

             ! lw radiation is coupled through updated t_g
             ! (no change in src_radiation.f90)
             ! sw radiation is coupled via albedo alb_rad
             ! (see change in src_radiation.f90)
          ELSE
             ! Update SSTs here because tgcom is not called anymore
             t_g(i,j,nnew) = t_s(i,j,nnew)
          ENDIF ifland                      ! MASK for land-points only
      ENDDO                        ! i loop
   ENDDO                          ! j loop



! COSMO
#endif
    END SUBROUTINE postproc_cosmo_clm

    REAL(DP) FUNCTION zfpi(zppa, rdocp)

      ! function zfpi for calculation of the exner function
      ! where the  argument zppa is pressure

      REAL(dp), INTENT(IN) :: zppa
      REAL(dp)             :: rdocp

      zfpi = (1.E-5_dp*zppa)**rdocp

    END FUNCTION zfpi

#ifdef ECHAM5

    SUBROUTINE postproc_echam_fesom

      USE messy_main_data_bi,       ONLY: sni, siced, seaice
      USE messy_main_constants_mem, ONLY: rho_h2o, rho_snow

      IMPLICIT NONE

      siced = siced/MAX(seaice,0.01_dp)
      sni = rho_snow/rho_h2o*sni/MAX(seaice,0.01_dp)

    END SUBROUTINE postproc_echam_fesom
#endif
! ECHAM5

! OASIS3MCT
#endif
  END SUBROUTINE oasis3mct_sndrcv


  ! ====================================================================
  SUBROUTINE oasis3mct_read_nml_cpl(status, iou)

    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

#ifdef OASIS3MCT
    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close
#endif
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

#ifdef OASIS3MCT
    NAMELIST /CPL/ LOCCOUPLE, RFIELD, SFIELD

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='oasis3mct_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

#else
    status = 0
#endif

  END SUBROUTINE oasis3mct_read_nml_cpl
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE oasis3mct_read_nml_cpl_extdata(status, iou)

    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------
#ifdef OASIS3MCT
    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

#ifdef OASIS3MCT
    NAMELIST /CPL_EXTDATA/ slf_file, slf_name, lon_name, lat_name

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='oasis3mct_read_nml_cpl_extdata'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    slf_file = ' '
    slf_name = ' '
    lon_name = ' '
    lat_name = ' '

    CALL read_nml_open(lex, substr, iou, 'CPL_EXTDATA', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_EXTDATA, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_EXTDATA', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

#else
    status = 0
#endif

  END SUBROUTINE oasis3mct_read_nml_cpl_extdata
  ! ====================================================================

  ! =========================================================================
  SUBROUTINE parse_postproc(string, FLD)

    USE messy_main_tools,         ONLY: strcrack, str2num
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),   INTENT(IN)    :: string
    TYPE(t_cpl_rfield), INTENT(INOUT) :: FLD

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER            :: substr ="parse_postproc"
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER  :: str1(:)  => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER  :: str2(:)  => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER  :: str3(:)  => NULL()
    INTEGER                                :: num1, num2, num3
    INTEGER                                :: ix
    INTEGER                                :: status

    status = -99
    num1 = -1
    num2 = -1
    num2 = -1

    IF (TRIM(ADJUSTL(string)) == ' ') RETURN

    write (*,*) 'FIELD requires postprocessing '

    ! CRACK STRING
    CALL strcrack(string,';',str1,num1)
    IF (num1 < 0 .OR. num1 > 33) THEN
       CALL error_bi('postproc string must contain colon seperated entries'&
            , substr)
    END IF

    ALLOCATE(FLD%pp(num1))
    DO ix = 1, num1
       CALL strcrack(str1(ix),':',str2,num2)
       SELECT CASE(TRIM(ADJUSTL(str2(1))))
       CASE('ABS')
          FLD%pp(ix)%action='ABS'
          FLD%pp(ix)%number=-9999.e-36
          FLD%pp(ix)%limit =-9999.e-36
          write (*,*) ' take absolute value'
       CASE('+')
          FLD%pp(ix)%action='+'
          CALL str2num(str2(2),FLD%pp(ix)%number, status)
          FLD%pp(ix)%limit =-9999.e-36
          IF (status /= 0) CALL error_bi( &
               'ERROR in string conversion +', substr)
          write (*,*) ' add constant: ', FLD%pp(ix)%number
       CASE('x')
          FLD%pp(ix)%action='x'
          CALL str2num(str2(2),FLD%pp(ix)%number, status)
          FLD%pp(ix)%limit =-9999.e-36
          IF (status /= 0) CALL error_bi( &
               'ERROR in string conversion x ', substr)
          write (*,*) ' multiply with factor: ', FLD%pp(ix)%number
       CASE('GT')
          FLD%pp(ix)%action='GT'
          CALL strcrack(str2(2),'=',str3,num3)
          CALL str2num(str3(1),FLD%pp(ix)%limit,  status)
          CALL str2num(str3(2),FLD%pp(ix)%number, status)
          IF (status /= 0) CALL error_bi( &
               'ERROR in string conversion GT ', substr)
          write (*,*) 'Greater than limit: ', FLD%pp(ix)%limit &
               ,' set to ', FLD%pp(ix)%number
          DEALLOCATE(str3); NULLIFY(str3)
          num3 = -1
       CASE('LT')
          FLD%pp(ix)%action='LT'
          CALL strcrack(str2(2),'=',str3,num3)
          CALL str2num(str3(1),FLD%pp(ix)%limit,  status)
          CALL str2num(str3(2),FLD%pp(ix)%number, status)
          IF (status /= 0) CALL error_bi( &
               'ERROR in string conversion LT ', substr)
          write (*,*) 'LESS than limit: ', FLD%pp(ix)%limit &
               ,' set to ', FLD%pp(ix)%number
          DEALLOCATE(str3); NULLIFY(str3)
          num3 = -1
       CASE('GE')
          FLD%pp(ix)%action='GT'
          CALL strcrack(str2(2),'=',str3,num3)
          CALL str2num(str3(1),FLD%pp(ix)%limit,  status)
          CALL str2num(str3(2),FLD%pp(ix)%number, status)
          IF (status /= 0) CALL error_bi( &
               'ERROR in string conversion GT ', substr)
          write (*,*) 'Greater/Equal than limit: ', FLD%pp(ix)%limit &
               ,' set to ', FLD%pp(ix)%number
          DEALLOCATE(str3); NULLIFY(str3)
          num3 = -1
       CASE('LE')
          FLD%pp(ix)%action='LT'
          CALL strcrack(str2(2),'=',str3,num3)
          CALL str2num(str3(1),FLD%pp(ix)%limit,  status)
          CALL str2num(str3(2),FLD%pp(ix)%number, status)
          IF (status /= 0) CALL error_bi( &
               'ERROR in string conversion LT ', substr)
          write (*,*) 'LESS/EQUAL than limit: ', FLD%pp(ix)%limit &
               ,' set to ', FLD%pp(ix)%number
          DEALLOCATE(str3); NULLIFY(str3)
          num3 = -1
       CASE('EXP')
          FLD%pp(ix)%action='EXP'
          CALL str2num(str2(2),FLD%pp(ix)%number, status)
          FLD%pp(ix)%limit =-9999.e-36
          IF (status /= 0) CALL error_bi( &
               'ERROR in string conversion EXP ', substr)
          write (*,*) ' exponentiate with: ', FLD%pp(ix)%number
       CASE DEFAULT
          CALL error_bi('Postproc string contains unknown options', substr)
       END SELECT

       DEALLOCATE(str2); NULLIFY(str2)
       num2 = -1
    END DO


    DEALLOCATE(str1); NULLIFY(str1)
    num1 = -1

  END SUBROUTINE parse_postproc
! **********************************************************************
  SUBROUTINE postprocess_field(FLD, cptr)

    IMPLICIT NONE

    TYPE(t_cpl_rfield), INTENT(INOUT) :: FLD
    REAL(dp), DIMENSION(:,:), POINTER :: cptr

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER            :: substr ="postprocess_field"
    INTEGER :: ix, jx, ii

    IF (.NOT. ALLOCATED(FLD%pp)) RETURN

    DO ii = 1, SIZE(FLD%pp)

       SELECT CASE(TRIM(FLD%pp(ii)%action))
       CASE('ABS')
          cptr = ABS(cptr)
       CASE('+')
          cptr = cptr + FLD%pp(ii)%number
       CASE('x')
          cptr = cptr * FLD%pp(ii)%number
       CASE('EXP')
          cptr = cptr**FLD%pp(ii)%number
       CASE('GT')
           DO ix = LBOUND(cptr,1), UBOUND(cptr,1)
             DO jx = LBOUND(cptr,2), UBOUND(cptr,2)
                IF (cptr(ix,jx) > FLD%pp(ii)%limit) &
                     cptr(ix,jx) = FLD%pp(ii)%number
             END DO
          END DO
       CASE('LT')
          DO ix = LBOUND(cptr,1), UBOUND(cptr,1)
             DO jx = LBOUND(cptr,2), UBOUND(cptr,2)
                IF (cptr(ix,jx) < FLD%pp(ii)%limit) &
                     cptr(ix,jx) = FLD%pp(ii)%number
             END DO
          END DO
       CASE('GE')
          DO ix = LBOUND(cptr,1), UBOUND(cptr,1)
             DO jx = LBOUND(cptr,2), UBOUND(cptr,2)

                IF (cptr(ix,jx) >= FLD%pp(ii)%limit) &
                     cptr(ix,jx) = FLD%pp(ii)%number
             END DO
          END DO
       CASE('LE')
          DO ix = LBOUND(cptr,1), UBOUND(cptr,1)
             DO jx = LBOUND(cptr,2), UBOUND(cptr,2)
                IF (cptr(ix,jx) <= FLD%pp(ii)%limit) &
                     cptr(ix,jx) = FLD%pp(ii)%number
             END DO
          END DO
       CASE DEFAULT
          CALL error_bi('Postproc string contains unknown options', substr)
       END SELECT
    END DO

  END SUBROUTINE postprocess_field
! **********************************************************************
#ifdef COSMO
  ! ====================================================================
  SUBROUTINE oasis3mct_write_grids_cosmo

#ifdef OASIS3MCT
    ! BMIL
    USE messy_main_mpi_bi,          ONLY: isubpos,  nboundlines, my_cart_id
    USE messy_main_grid_def_mem_bi, ONLY: ie_tot, je_tot, ie, je         &
                                        , istart, jstart                 &
                                        , istartpar, jstartpar           &
                                        , iendpar, jendpar               &
                                        , startlon_tot, startlat_tot     &
                                        , dlon, dlat                     &
                                        , pollon, pollat, polgam         &
                                        , startlon, startlat
    USE messy_main_grid_def_bi,     ONLY: gboxarea_2d
    USE messy_main_data_bi,         ONLY: slf, slm, alake
    USE MESSY_MAIN_GRID_TRAFO,      ONLY: phirot2phi,  rlarot2rla
    USE netcdf
#endif
    IMPLICIT NONE

#ifdef OASIS3MCT
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'oasis3mct_write_grids_cosmo'
    INTEGER                     :: status

    !vector of integers describing the local partition in the global index space
    INTEGER                             :: dummy_flag
    REAL(dp), DIMENSION(:,:),   POINTER :: lon      => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: lat      => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: lon_ptr  => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: lat_ptr  => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: loni     => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: lati     => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: loni_ptr => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: lati_ptr => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: corner_lon  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: corner_lat  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: corner_slon => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: corner_slat => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: tmp_mask    => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: tmp_imask   => NULL()
    REAL(dp)                            :: rlon
    REAL(dp)                            :: rlat
    REAL(dp)                            :: intfac = 100000._dp
    INTEGER                             :: istartlon, istartlat
    INTEGER                             :: idlon, idlat, itmp
    INTEGER                             :: illsub, jllsub
    INTEGER                             :: ix, jx, iy, jy, iix, jjx
    INTEGER,  DIMENSION(:,:), POINTER   :: imask  => NULL()
    REAL(dp), DIMENSION(:,:), POINTER   :: area   => NULL()
    CHARACTER(LEN=1)                    :: addstr= ' '

    INTEGER :: fileid, varid
    INTEGER :: ndims
    INTEGER :: dimids(4) = 0
    INTEGER :: iedim, jedim
    INTEGER :: iss, jss

    REAL(dp), DIMENSION(:),   POINTER :: vlon => NULL()
    REAL(dp), DIMENSION(:),   POINTER :: vlat => NULL()

    CALL start_message_bi(modstr,'  ',substr)

    ! READ sea-land fraction from external data file
    IF ((TRIM(slf_file) /= '') .AND. (TRIM(slf_name) /= '') &
          .AND. (TRIM(lon_name) /= '') .AND. (TRIM(lat_name) /= '')) THEN

          ! OPEN FILE
          status = nf90_open(TRIM(slf_file), nf90_nowrite, fileid)
          IF (status /= NF90_noerr) &
               CALL error_bi('error opening slf-file', substr)


          ! INQUIRE LON / LAT
          status = nf90_inq_varid (fileid, TRIM(lon_name), varid)
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing lon-var varid', substr)

          status =  NF90_inquire_variable(fileid, varid, dimids = dimids(1:1))
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing lon-var dimids', substr)

          status = NF90_inquire_dimension (fileid, dimids(1), len=iedim)
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing lon-var dim1', substr)

          ALLOCATE(vLON(iedim))

          status=NF90_get_var (fileid, varid, vLON, start=(/1/), count=(/iedim/))

          iss = -999
          DO ix = 1,iedim
             IF (vLON(ix) >= startlon) THEN
                iss = ix
                EXIT
             END IF
          END DO
          IF (iss == -999) CALL error_bi('wrong input data : LON', substr)
          DEALLOCATE(vLON)

          status = nf90_inq_varid (fileid, TRIM(lat_name), varid)
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing lon-var varid', substr)

          status =  NF90_inquire_variable(fileid, varid, dimids = dimids(1:1))
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing lon-var dimids', substr)

          status = NF90_inquire_dimension (fileid, dimids(1), len=jedim)
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing lon-var dim1', substr)

          ALLOCATE(vLAT(jedim))

          status=NF90_get_var (fileid, varid, vLAT, start=(/1/), count=(/jedim/))

          jss = -999
          DO ix = 1,jedim
             IF (vLAT(ix) >= startlat) THEN
                jss = ix
                EXIT
             END IF
          END DO
          IF (jss == -999) CALL error_bi('wrong input data : LAT', substr)
          DEALLOCATE(vLAT)

          ! INQUIRE VARIABLE
          status = nf90_inq_varid (fileid, TRIM(slf_name), varid)
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing slf-var varid', substr)

          status =  NF90_inquire_variable(fileid, varid, ndims = ndims)
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing slf-var ndims', substr)


          status =  NF90_inquire_variable(fileid, varid, dimids = dimids(1:ndims))
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing slf-var dimids', substr)

          status = NF90_inquire_dimension (fileid, dimids(1), len=iedim)
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing slf-var dim1', substr)
          status = NF90_inquire_dimension (fileid, dimids(2), len=jedim)
          IF (status /= NF90_noerr) &
               CALL error_bi('error inquireing slf-var dim2', substr)

         ! READ VARIABLE
          status = NF90_get_var (fileid, varid, slf,     &
               start=(/iss,jss/), count=(/ie, je/))
          IF (status /= NF90_noerr) &
               CALL error_bi('error reading slf-var ', substr)

          ! CLOSE FILE
          status = nf90_close(fileid)
          IF (status /= NF90_noerr)  &
               CALL error_bi('error closing slf-file', substr)
   END IF

   ! ---------------------------------------------------------------
    iis = istartpar
    iie = iendpar
    jjs = jstartpar
    jje = jendpar
    xdim_tot = ie_tot
    ydim_tot = je_tot

   CALL start_message_bi(modstr,'DEFINE PARTITION',substr)
    !  CALL oasis_set_debug(20)

    ALLOCATE(paral(5))
    paral(CLIM_Strategy) = CLIM_box ! box partition
    ! the upper left corner global offset
    paral(CLIM_offset) = xdim_tot * ( isubpos(my_cart_id,2) - (jstart-jjs) -1) &
                + ( isubpos(my_cart_id,1) - (istart-iis) - 1)
    paral(CLIM_Length) = iie - istartpar + 1 ! the local extent in x
    paral(4) = jje - jjs + 1 ! the local extent in y
    paral(5) = xdim_tot                  ! the global extent in x

    CALL oasis_def_partition(part_id, paral, status)
    CALL   end_message_bi(modstr,'DEFINE PARTITION',substr)

    !=========================================================================

    CALL start_message_bi(modstr,'DEFINE GRIDS',substr)

    ! access basemodelgrid

    ! get box mid point lon / lat
    ! including a boundary of 1 at each side for interfaces of staggered grids
    ALLOCATE(lon(0:ie+1,0:je+1))
    ALLOCATE(lat(0:ie+1,0:je+1))

    illsub = isubpos(my_cart_id,1) - nboundlines - 1
    jllsub = isubpos(my_cart_id,2) - nboundlines - 1

    istartlon = NINT(startlon_tot * intfac)
    istartlat = NINT(startlat_tot * intfac)
    idlon     = NINT(dlon * intfac)
    idlat     = NINT(dlat * intfac)

    DO ix = 0, ie+1
       itmp = istartlon + idlon * (illsub + ix -1)
       rlon = REAL(itmp,dp) / intfac

       DO jx = 0, je+1
          ! j index 0
          itmp = istartlat + idlat * (jllsub + jx - 1)
          rlat =  REAL(itmp,dp) / intfac

          lon(ix,jx) = rlarot2rla ( rlat , rlon , pollat, pollon, polgam)
          lat(ix,jx) = phirot2phi ( rlat , rlon , pollat, polgam)

       END DO
    END DO

    lon_ptr => lon(iis:iie,jjs:jje)
    lat_ptr => lat(iis:iie,jjs:jje)

    ! get box interface lon / lat
    ALLOCATE(loni(ie+1,je+1))
    ALLOCATE(lati(ie+1,je+1))

    DO ix = 1, ie+1
       itmp = istartlon + idlon * (illsub + ix -1-0.5)
       rlon = REAL(itmp,dp) / intfac

       DO jx = 1, je+1
          ! j index 0
          itmp = istartlat + idlat * (jllsub + jx - 1-0.5)
          rlat =  REAL(itmp,dp) / intfac

          loni(ix,jx) = rlarot2rla ( rlat , rlon , pollat, pollon, polgam)
          lati(ix,jx) = phirot2phi ( rlat , rlon , pollat, polgam)
       END DO
    END DO


    loni_ptr => loni(iis:iie,jjs:jje)
    lati_ptr => lati(iis:iie,jjs:jje)

    ! calculate corners
    ALLOCATE(corner_lon(paral(3),paral(4),4))
    ALLOCATE(corner_lat(paral(3),paral(4),4))
    ALLOCATE(corner_slon(paral(3),paral(4),4))
    ALLOCATE(corner_slat(paral(3),paral(4),4))

    ! corner lower left
    ! -- mid point grid
    corner_lon(:,:,1) = loni(iis:iie,jjs:jje)
    corner_lat(:,:,1) = lati(iis:iie,jjs:jje)
    ! -- staggered grid
    corner_slon(:,:,1) = lon(iis:iie,jjs:jje)
    corner_slat(:,:,1) = lat(iis:iie,jjs:jje)

    ! corner lower right
    corner_lon(:,:,2) = loni(iis+1:iie+1,jjs:jje)
    corner_lat(:,:,2) = lati(iis+1:iie+1,jjs:jje)
    ! -- staggered grid
    corner_slon(:,:,2) = lon(iis+1:iie+1,jjs:jje)
    corner_slat(:,:,2) = lat(iis+1:iie+1,jjs:jje)

    ! corner upper right
    corner_lon(:,:,3) = loni(iis+1:iie+1,jjs+1:jje+1)
    corner_lat(:,:,3) = lati(iis+1:iie+1,jjs+1:jje+1)
    ! -- staggered grid
    corner_slon(:,:,3) = lon(iis+1:iie+1,jjs+1:jje+1)
    corner_slat(:,:,3) = lat(iis+1:iie+1,jjs+1:jje+1)

    ! corner upper left
    corner_lon(:,:,4) = loni(iis:iie,jjs+1:jje+1)
    corner_lat(:,:,4) = lati(iis:iie,jjs+1:jje+1)
    ! -- staggered grid
    corner_slon(:,:,4) = lon(iis:iie,jjs+1:jje+1)
    corner_slat(:,:,4) = lat(iis:iie,jjs+1:jje+1)

   CALL   end_message_bi(modstr,'DEFINE GRIDS',substr)
! ==========================================================================

!   CALL start_message_bi(modstr,'DEFINE MASKS',substr)

   CALL oasis_start_grids_writing (dummy_flag)

    ! DEFINE LOGICAL LANDMASK
    ALLOCATE(llandmask(SIZE(slf,1),SIZE(slf,2)))
    llandmask = .FALSE.
    DO ix = 1,SIZE(llandmask,1)
       DO jx = 1,SIZE(llandmask,2)
          IF (slf(ix,jx) >= 0.5_dp) llandmask(ix,jx)=.TRUE.
       END DO
    END DO

    ALLOCATE(imask(iis:iie,jjs:jje))
   ! unmasked by default
   imask = 0

   ALLOCATE(area(iis:iie,jjs:jje))
   ! CLM^2 set area = 1
   ! true area required, if conservative interpolation used
   area = 1._dp

   DO ix = 1, NUMMOD
      IF (LOCCPL(ix)%lactive) THEN
!         write (*,*) 'WRITE MASK TYPE ',LOCCPL(ix)%masktype &
!              ,' FOR ',TRIM(LOCCPL(ix)%modname)

         SELECT CASE(LOCCPL(ix)%masktype)
         CASE(0)
            ! INVOLVE ALL POINTS
            imask = 0
            area(iis:iie,jjs:jje)= gboxarea_2d(iis:iie,jjs:jje)
            addstr=' '
         CASE(1)
            ! ACTIVATE (0) ONLY POINTS OVER LAND / DEACTIVATE (1) SEA POINTS
            imask = 1
            DO iix = iis, iie
               DO jjx = jjs, jje
                  IF (llandmask(iix,jjx)) imask(iix,jjx)  = 0
                  area(iix,jjx)= gboxarea_2d(iix,jjx)
               END DO
            END DO
            addstr='l'
         CASE(2)
            ! ACTIVATE (0) ONLY POINTS OVER SEA / DEACTIVATE (1) LAND POINTS
            imask = 0
            DO iix = iis, iie
               DO jjx = jjs, jje
                  IF (llandmask(iix,jjx)) imask(iix,jjx)  = 1
                  area(iix,jjx)= gboxarea_2d(iix,jjx)
               END DO
            END DO
            addstr='s'
         CASE(3)
            ! ACTIV: ALL WATER POINTS (SEA AND LAKE) / DEACTIV. (1) LAND POINTS
            imask = 0
            ALLOCATE(tmp_mask(iis:iie,jjs:jje))
            tmp_mask(iis:iie,jjs:jje) =  &
                 slf(iis:iie,jjs:jje) + alake(iis:iie,jjs:jje)
            DO iix = iis, iie
               DO jjx = jjs, jje
                  IF (tmp_mask(iix,jjx) > 0.9999999_dp) imask(iix,jjx) = 1
                  area(iix,jjx)= gboxarea_2d(iix,jjx)
                  IF (tmp_mask(iix,jjx) > 0.0_dp .AND. imask(iix, jjx) == 0) &
                       area(iix,jjx) = (1.0_dp-tmp_mask(iix,jjx)) * area(iix,jjx)
                  IF (imask(iix,jjx) == 1) area(iix,jjx) = 0.0_dp
               END DO
            END DO

            DEALLOCATE(tmp_mask)
            addstr='w'
         CASE(4)
            ! ACTIV: ALL WATER POINTS (SEA AND LAKE) / DEACTIV. (1) LAND POINTS
            imask = 0
            ALLOCATE(tmp_imask(iis:iie,jjs:jje))
            tmp_imask(iis:iie,jjs:jje) = INT(slm(iis:iie,jjs:jje))
            DO iix = iis, iie
               DO jjx = jjs, jje
                  IF (alake(iix,jjx) > 0.0001_dp) THEN
                     imask(iix,jjx) = 1
                     area(iix,jjx)= 0._dp
                  ELSE
                     imask(iix,jjx) = 0
                     area(iix,jjx)= gboxarea_2d(iix,jjx)
                  END IF
               END DO
            END DO

            DEALLOCATE(tmp_imask)
            addstr='k'
         CASE DEFAULT
            ! ACTIVATE ALL POINTS
            imask = 0
            CALL warning_bi('mask not yet implemented => activate all points'&
                 , substr)
         END SELECT

         CALL oasis_write_grid (TRIM(compname)//TRIM(addstr)  &
              ,  xdim_tot, ydim_tot,  lon_ptr,  lat_ptr,  part_id)
         CALL oasis_write_corner(TRIM(compname)//TRIM(addstr) &
              , xdim_tot, ydim_tot,4, corner_lon,  corner_lat,  part_id)
         CALL oasis_write_mask (TRIM(compname)//TRIM(addstr)  &
              ,  xdim_tot, ydim_tot, imask, part_id)
         CALL oasis_write_area (TRIM(compname)//TRIM(addstr)  &
              , xdim_tot, ydim_tot, area, part_id)

         ! SAVE MASK for POSTPROCESSING
         ALLOCATE(LOCCPL(ix)%mask(iis:iie,jjs:jje))
         LOCCPL(ix)%mask = .TRUE.
         DO iy = iis, iie
            DO jy = jjs, jje
               IF (imask(iy,jy) == 1) LOCCPL(ix)%mask(iy,jy) = .FALSE.
            END DO
         END DO
      END IF
   END DO

   DEALLOCATE(imask); NULLIFY(imask)
   DEALLOCATE(area);  NULLIFY(area)
   DEALLOCATE(lon);   NULLIFY(lon)
   DEALLOCATE(lat);   NULLIFY(lat)
   IF(ASSOCIATED(loni)) THEN
      DEALLOCATE(loni); NULLIFY(loni)
   END IF
   IF(ASSOCIATED(lati)) THEN
      DEALLOCATE(lati); NULLIFY(lati)
   END IF
   NULLIFY(lon_ptr)
   NULLIFY(loni_ptr)
   NULLIFY(lat_ptr)
   NULLIFY(lati_ptr)
   DEALLOCATE(corner_lon);
   NULLIFY(corner_lon)
   DEALLOCATE(corner_lat);
   NULLIFY(corner_lat)
   IF (ASSOCIATED(corner_slon)) DEALLOCATE(corner_slon)
   NULLIFY(corner_slon)
   IF (ASSOCIATED(corner_slat)) DEALLOCATE(corner_slat)
   NULLIFY(corner_slat)

   CALL   end_message_bi(modstr,'DEFINE MASKS, AREA',substr)

! ==========================================================================

   CALL oasis_terminate_grids_writing()

! ==========================================================================
   CALL   end_message_bi(modstr,'  ',substr)
#endif
! OASIS
! ==========================================================================
 END SUBROUTINE oasis3mct_write_grids_cosmo
#endif
!COSMO

#ifdef ECHAM5
  ! ====================================================================
  SUBROUTINE oasis3mct_write_grids_echam

#ifdef OASIS3MCT
    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, gather_field, dcl

    USE messy_main_data_bi,          ONLY: slf, slm, alake
    USE messy_main_grid_def_mem_bi,  ONLY: nlon, ngl, nproma, npromz, ngpblks
    USE messy_main_grid_def_bi,      ONLY: philon, philat, gridarea
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL
    USE messy_main_channel_repr,     ONLY: t_representation, get_representation
    USE messy_main_channel,          ONLY: get_channel_object
#endif

    IMPLICIT NONE

#ifdef OASIS3MCT
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'oasis3mct_write_grids_echam'
    INTEGER                     :: status

    !vector of integers describing the local partition in the global index space
    INTEGER                             :: dummy_flag
    REAL(dp), DIMENSION(:,:),   POINTER :: lon        => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: lat        => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: corner_lon => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: corner_lat => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER :: tmp_mask   => NULL()
    INTEGER                             :: ix, jx, iy, jy, kproma
    INTEGER,  DIMENSION(:,:), POINTER   :: imask  => NULL()
    REAL(dp), DIMENSION(:,:), POINTER   :: area   => NULL()
    CHARACTER(LEN=1)                    :: addstr= ' '
    !INTEGER                            :: nsegs ! number of segments
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: p4_slf   => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: gl_slf   => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: p4_slm   => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: gl_slm   => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: p4_alake => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: gl_alake => NULL()
    TYPE(t_representation),       POINTER   :: repr     => NULL()

    iis = 1
    jjs = 1
    iie = nproma  !dcl%nglon
    jje = ngpblks !dcl%nglat

    xdim_tot = nlon
    ydim_tot = ngl

    CALL start_message_bi(modstr,'DEFINE PARTITION',substr)
    !  CALL oasis_set_debug(20)

    nsegs    = dcl%nglat
    nsegslen = dcl%nglon

    ALLOCATE(paral(2 + 2 * nsegs))
    paral(:) = -99
    paral(CLIM_Strategy) = CLIM_Orange
    paral(CLIM_Segments) = nsegs
    DO ix = 1, nsegs
       paral(CLIM_Segments + ix*2-1) =  (dcl%glat(ix)-1) * nlon +  dcl%glon(ix)
       paral(CLIM_Segments + ix*2  ) =   nsegslen
    END DO

    CALL oasis_def_partition(part_id, paral, status)
    CALL   end_message_bi(modstr,'DEFINE PARTITION',substr)

    !=========================================================================
    CALL start_message_bi(modstr,'DEFINE GRIDS',substr)
    ! get global masks
    CALL get_representation(status, GP_2D_HORIZONTAL, repr)
    CALL channel_halt(substr, status)

    ALLOCATE (gl_slf  (repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ))
    ALLOCATE (gl_slm  (repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ))
    ALLOCATE (gl_alake(repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ))

    CALL get_channel_object(status, 'g3b', 'slf',   p4=p4_slf)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'g3b', 'slm',   p4=p4_slm)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'g3b', 'alake', p4=p4_alake)
    CALL channel_halt(substr, status)

    CALL gather_field(gl_slf,   repr%gdimlen, p4_slf)
    CALL gather_field(gl_slm,   repr%gdimlen, p4_slm)
    CALL gather_field(gl_alake, repr%gdimlen, p4_alake)

    NULLIFY(p4_slf)
    NULLIFY(p4_slm)
    NULLIFY(p4_alake)

    ! define grid in global (!) decomposition
    IF (p_parallel_io) THEN

       ALLOCATE(lon(nlon,ngl))
       ALLOCATE(lat(nlon,ngl))

       DO ix = 1, nlon
          lat(ix,:) = philat(:)
       ENDDO
       DO jx = 1, ngl
          lon(:,jx) = philon(:)
       ENDDO
       WHERE (lon(:,:) < 0._dp)
          lon = lon + 360._dp
       END WHERE
       WHERE (lon(:,:) >= 360._dp)
          lon = lon - 360._dp
       END WHERE

       ALLOCATE(corner_lon(nlon,ngl,4))
       ALLOCATE(corner_lat(nlon,ngl,4))

       DO ix = 1, nlon-1
          corner_lon(ix,:,1) = 0.5_dp * (lon(ix+1,1) + lon(ix,1))
          corner_lon(ix,:,4) = 0.5_dp * (lon(ix+1,1) + lon(ix,1))
       ENDDO

       DO ix = 2, nlon
          corner_lon(ix,:,2) = 0.5_dp * (lon(ix-1,1) + lon(ix,1))
          corner_lon(ix,:,3) = 0.5_dp * (lon(ix-1,1) + lon(ix,1))
       ENDDO

       corner_lon(nlon,:,1) = 0.5_dp * (lon(1,1) + lon(nlon,1) + 360._dp)
       corner_lon(nlon,:,4) = 0.5_dp * (lon(1,1) + lon(nlon,1) + 360._dp)
       corner_lon(   1,:,2) = 0.5_dp * (lon(nlon,1) + lon(1,1) - 360._dp)
       corner_lon(   1,:,3) = 0.5_dp * (lon(nlon,1) + lon(1,1) - 360._dp)

       IF (corner_lon(nlon,1,1) >= 360._dp) &
            corner_lon(nlon,:,1) = corner_lon(nlon,:,1) - 360._dp
       IF (corner_lon(nlon,1,4) >= 360._dp) &
            corner_lon(nlon,:,4) = corner_lon(nlon,:,4) - 360._dp
       IF (corner_lon(nlon,1,2) < 0._dp)    &
            corner_lon(nlon,:,2) = corner_lon(nlon,:,2) + 360._dp
       IF (corner_lon(nlon,1,3) < 0._dp)    &
            corner_lon(nlon,:,3) = corner_lon(nlon,:,3) + 360._dp

       DO jx = 1, ngl-1
          corner_lat(:,jx,3) = 0.5_dp * (lat(1,jx+1) + lat(1,jx))
          corner_lat(:,jx,4) = 0.5_dp * (lat(1,jx+1) + lat(1,jx))
       ENDDO

       DO jx = 2, ngl
          corner_lat(:,jx,1) = 0.5_dp * (lat(1,jx-1) + lat(1,jx))
          corner_lat(:,jx,2) = 0.5_dp * (lat(1,jx-1) + lat(1,jx))
       ENDDO

       corner_lat(:,  1,1) =  90._dp
       corner_lat(:,  1,2) =  90._dp
       corner_lat(:,ngl,3) = -90._dp
       corner_lat(:,ngl,4) = -90._dp

!  ==========================================================================

       CALL oasis_start_grids_writing (dummy_flag)

       ALLOCATE(imask(nlon, ngl))
       ALLOCATE(area(nlon,ngl))
       ! unmasked by default
       imask = 0
       area = 1._dp

       DO ix = 1, NUMMOD
          IF (LOCCPL(ix)%lactive) THEN
             !         write (*,*) 'WRITE MASK TYPE ',LOCCPL(ix)%masktype &
             !              ,' FOR ',TRIM(LOCCPL(ix)%modname)
             SELECT CASE(LOCCPL(ix)%masktype)
             CASE(0)
                ! INVOLVE ALL POINTS
                imask = 0
                ! set area here
                DO jx = 1,ngl
                   area(:,jx) = gridarea(jx)
                END DO
                addstr=' '
             CASE(1)
                ! ACTIVATE (0) ONLY POINTS OVER LAND / DEACTIVATE (1) SEA POINTS
                imask = 0
                WHERE (gl_slf(:,:,1,1) < 0.5_dp) imask(:,:) = 1
                ! set area here
                DO jx = 1,ngl
                   area(:,jx) = gridarea(jx)
                END DO
                addstr='l'
             CASE(2)
                ! ACTIVATE (0) ONLY POINTS OVER SEA / DEACTIVATE (1) LAND POINTS
                imask = 0
                WHERE (gl_slf(:,:,1,1) >= 0.5_dp) imask(:,:) = 1
                ! set area here
                DO jx = 1,ngl
                   area(:,jx) = gridarea(jx)
                END DO
                addstr='s'
             CASE(3)
                ! ACTIV: ALL WATER POINTS (SEA AND LAKE) / DEACTIV. (1) LAND POINTS
                imask = 0
                ALLOCATE(tmp_mask(nlon,ngl))
                tmp_mask(:,:) = gl_slf(:,:,1,1) + gl_alake(:,:,1,1)
                WHERE (tmp_mask(:,:) > 0.9999999_dp) imask(:,:) = 1
                ! calc area
                DO jx = 1,ngl
                   area(:,jx) = gridarea(jx)
                END DO
                WHERE (tmp_mask(:,:) > 0.0_dp .AND. imask(:,:) == 0)
                   area(:,:) = (1.0_dp-tmp_mask(:,:)) * area(:,:)
                END WHERE
                WHERE (imask(:,:) == 1)
                   area(:,:) = 0.0_dp
                END WHERE
                DEALLOCATE(tmp_mask)
                addstr='w'
             CASE(4)
                ! ACTIV: ALL WATER POINTS  SLM/ DEACTIV. (1) LAND POINTS
                imask(:,:) = INT(gl_slm(:,:,1,1))
                WHERE (gl_alake(:,:,1,1) > 0.0001_dp) imask(:,:) = 1
                ! calc area
                DO jx = 1,ngl
                   area(:,jx) = gridarea(jx)
                END DO
                WHERE (imask(:,:) == 1)
                   area(:,:) = 0.0_dp
                END WHERE
                addstr='k'
             CASE DEFAULT
                ! ACTIVATE ALL POINTS
                imask = 0
                CALL warning_bi(&
                     'mask not yet implemented => activate all points'&
                     , substr)
             END SELECT

             CALL oasis_write_grid (TRIM(compname)//TRIM(addstr)  &
                  ,  xdim_tot, ydim_tot,  lon,  lat)
             CALL oasis_write_corner(TRIM(compname)//TRIM(addstr) &
                  , xdim_tot, ydim_tot,4, corner_lon,  corner_lat)
             CALL start_message_bi(modstr,'DEFINE MASKS',substr)
             CALL oasis_write_mask (TRIM(compname)//TRIM(addstr)  &
                  ,  xdim_tot, ydim_tot, imask)
             CALL oasis_write_area (TRIM(compname)//TRIM(addstr)  &
                  , xdim_tot, ydim_tot, area)
          END IF
       END DO

       CALL oasis_terminate_grids_writing()

       DEALLOCATE(imask);      NULLIFY(imask)
       DEALLOCATE(area);       NULLIFY(area)
       DEALLOCATE(lon);        NULLIFY(lon)
       DEALLOCATE(lat);        NULLIFY(lat)
       DEALLOCATE(corner_lon); NULLIFY(corner_lon)
       DEALLOCATE(corner_lat); NULLIFY(corner_lat)
    END IF ! p_parallel_io

    DEALLOCATE (gl_slf);    NULLIFY(gl_slf)
    DEALLOCATE (gl_slm);    NULLIFY(gl_slm)
    DEALLOCATE (gl_alake);  NULLIFY(gl_alake)
    NULLIFY(repr)

    CALL   end_message_bi(modstr,'DEFINE MASKS, AREA',substr)

! ==========================================================================

    DO ix = 1, NUMMOD
      IF (LOCCPL(ix)%lactive) THEN
         ALLOCATE(imask(nproma,ngpblks))
         SELECT CASE(LOCCPL(ix)%masktype)
         CASE(0)
            ! INVOLVE ALL POINTS
            imask = 0
         CASE(1)
            ! ACTIVATE (0) ONLY POINTS OVER LAND / DEACTIVATE (1) SEA POINTS
            imask = 0
            WHERE (slf(:,:) < 0.5_dp) imask(:,:) = 1
         CASE(2)
            ! ACTIVATE (0) ONLY POINTS OVER SEA / DEACTIVATE (1) LAND POINTS
            imask = 0
            WHERE (slf(:,:) >= 0.5_dp) imask(:,:) = 1
         CASE(3)
            ! ACTIV: ALL WATER POINTS (SEA AND LAKE) / DEACTIV. (1) LAND POINTS
            imask = 0
            WHERE (slf(:,:) + alake(:,:) > 0.9999999_dp) imask(:,:) = 1
         CASE(4)
            ! ACTIV: ALL WATER POINTS  SLM/ DEACTIV. (1) LAND POINTS
            imask(:,:) = INT(slm(:,:))
            WHERE (alake(:,:) > 0.0001_dp) imask(:,:) = 1
         CASE DEFAULT
            ! ACTIVATE ALL POINTS
            imask = 0
            CALL warning_bi('mask not yet implemented => activate all points'&
                 , substr)
         END SELECT


         ALLOCATE(LOCCPL(ix)%mask(nproma,ngpblks))
         LOCCPL(ix)%mask = .FALSE.
         IF (LOCCPL(ix)%modname(1:5) == 'fesom') THEN
            ! special case for different definitions of ocean mask
            DO jy = 1, ngpblks
               kproma = nproma
               IF (jy == ngpblks) kproma = npromz
               DO iy = 1, kproma
                  IF (slf(iy,jy) < 1._dp .AND. alake(iy,jy) == 0._dp) &
                       LOCCPL(ix)%mask(iy,jy) = .TRUE.
               END DO
            END DO
         ELSE
            DO jy = 1, ngpblks
               kproma = nproma
               IF (jy == ngpblks) kproma = npromz
               DO iy = 1, kproma
                  IF (imask(iy,jy) == 0) LOCCPL(ix)%mask(iy,jy) = .TRUE.
               END DO
            END DO
         END IF
         DEALLOCATE(imask); NULLIFY(imask)
      END IF
   END DO
! ==========================================================================
   CALL   end_message_bi(modstr,'  ',substr)

! ==========================================================================
! ECHAM
#endif
  END SUBROUTINE oasis3mct_write_grids_echam

  ! ====================================================================
! OASIS3MCT
#endif
  SUBROUTINE parse_location(status, string, isubroutine, ibegend, iorder)

    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
    USE messy_main_tools,         ONLY: strcrack, str2num
    USE messy_main_control,       ONLY: MEP_ENTRY_MAX, entry_name &
                                      , MEPS_BEGIN, MEPS_END


    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: string
    INTEGER,          INTENT(OUT) :: isubroutine
    INTEGER,          INTENT(OUT) :: ibegend
    INTEGER,          INTENT(OUT) :: iorder

    ! LOCAL
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER  :: strchunk(:)  => NULL()
    INTEGER                                :: chunknum
    INTEGER                                :: ix

    ! INITIALIZE
    status = -1
    ibegend = MEPS_END ! set default
    iorder  = -99

    ! CRACK STRING
    CALL strcrack(string,':',strchunk,chunknum)
    IF (chunknum < 1 .OR. chunknum > 3) THEN
       status = 100  ! location string must contain 1-3 colon seperated entries
       RETURN
    END IF

    ! Find subroutine name and index
    DO ix=1,MEP_ENTRY_MAX
       IF (TRIM(ADJUSTL(strchunk(1))) == TRIM(entry_name(ix))) THEN
          isubroutine = ix
          EXIT
       END IF
    END DO
    IF (ix > MEP_ENTRY_MAX) THEN
       status = 101 ! parse_location: subroutine names do not much
       RETURN
    END IF

    IF (chunknum < 2) RETURN

    ! couple at subroutine start or end
    IF (TRIM(ADJUSTL(strchunk(2))) == 'B') THEN
       ibegend = MEPS_BEGIN
    ELSE IF (TRIM(ADJUSTL(strchunk(2))) == 'E') THEN
       ibegend = MEPS_END
    ELSE
       status = 102 ! parse_location: second chunk must equal B or E
       RETURN
    END IF

    IF (chunknum < 3) RETURN
    ! fix sequence
    CALL str2num(strchunk(3),iorder)

    DEALLOCATE(strchunk); NULLIFY(strchunk)
    chunknum = 0

  END SUBROUTINE parse_location

 ! ====================================================================
END MODULE messy_oasis3mct_si
! **********************************************************************
