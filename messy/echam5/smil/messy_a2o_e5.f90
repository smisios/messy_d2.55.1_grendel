! ************************************************************************
MODULE messy_a2o_e5
! ************************************************************************

  ! AUTHOR:
  !  Pozzer Andrea, MPICH, August 2007

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_main_timer_bi,      ONLY: p_bcast_event, timer_event_init
  USE messy_main_timer_event,   ONLY: io_time_event, TRIG_FIRST   &
                                    , time_event
  USE messy_mpiom_mem_e5,       ONLY: ie,je,ie_g, je_g
  USE messy_mpiom,              ONLY: L_COUPLING
  USE messy_main_grid_def_mem_bi, ONLY: nlat=>ngl, nlon, nlev
  USE messy_main_data_bi,       ONLY:lcouple
  USE messy_a2o
  ! um_ak_20130502+
  ! USE messy_main_gridtrafo_scrip
  USE messy_main_grid_trafo_scrp_base
  ! um_ak_20130502-

  IMPLICIT NONE
  PRIVATE

  REAL(dp), PARAMETER :: zepsec = 1.E-12_dp

  LOGICAL :: LRUNSM = .TRUE.  ! RUN THIS SUBMODEL
  LOGICAL :: LRUNSM_A2O = .FALSE.  ! RUN A2O TRANSFORMATION
  LOGICAL :: LRUNSM_O2A = .FALSE.  ! RUN O2A TRANSFORMATION

  ! column related
  INTEGER, PARAMETER, PUBLIC :: A2O_CSUM = 1
  INTEGER, PARAMETER, PUBLIC :: A2O_CAVE = 2
  INTEGER, PARAMETER, PUBLIC :: A2O_CMAX = 3
  INTEGER, PARAMETER, PUBLIC :: A2O_CLLEV = 4
  INTEGER, PARAMETER, PUBLIC :: A2O_CHLEV = 5
  
  ! time related
  INTEGER, PARAMETER, PUBLIC :: A2O_TSUM = 1
  INTEGER, PARAMETER, PUBLIC :: A2O_TAVE = 2


  !TIME MANAGMENT!
  LOGICAL, PUBLIC, SAVE ::  l_trig_a2o = .TRUE.  
  LOGICAL, PUBLIC, SAVE ::  l_trig_o2a = .TRUE.  
  TYPE(io_time_event), PUBLIC, SAVE :: trig_a2o = & 
       io_time_event (1,'days',TRIG_FIRST,0)
  TYPE(io_time_event), PUBLIC, SAVE :: trig_o2a = & 
       io_time_event (1,'days',TRIG_FIRST,0)
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_a2o
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_o2a
  
  ! #####################
  ! ATMOSPHERE -> OCEAN
  ! #####################

  TYPE IO_A2O
     CHARACTER(LEN=STRLEN_OBJECT)  :: name = ''
     CHARACTER(LEN=STRLEN_CHANNEL) :: a_channel = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: a_object = ''
     INTEGER                       :: cmethod = A2O_CLLEV
     INTEGER                       :: tmethod = A2O_TAVE
     LOGICAL                       :: lmask = .FALSE.
     REAL(DP)                      :: min_value = -1.E34_dp
     REAL(DP)                      :: max_value = 1.E34_dp
  END TYPE IO_A2O

  ! interpolation specific
  TYPE INT_A2O
     CHARACTER(char_len) :: source_grid = ''
     CHARACTER(char_len) :: dest_grid   = ''
     CHARACTER(char_len) :: map_method  = ''
     CHARACTER(char_len) :: normalize_opt = ''
  END TYPE INT_A2O

  LOGICAL, DIMENSION(4,4,3) :: INTCALCA2O = .FALSE.

  TYPE INT_XA2O
     INTEGER :: source_grid 
     INTEGER :: dest_grid   
     INTEGER :: map_method  
     INTEGER :: norm_opt 
  END TYPE INT_XA2O

  ! TRANSFORMATION INFORMATIONS
  TYPE A2O_transformation
    INTEGER                          :: nn                      ! number of neighbors
    INTEGER, DIMENSION(:),  POINTER  :: e5_ilon_g   => NULL()   ! e5 lon index global
    INTEGER, DIMENSION(:),  POINTER  :: e5_ilat_g   => NULL()   ! e5 lat index global
    REAL(DP), DIMENSION(:), POINTER  :: weight      => NULL()   ! weight per point
    REAL(DP), POINTER                :: frac        => NULL()   ! fractional area per point
  END TYPE A2O_transformation

  TYPE T_A2O
     TYPE(IO_A2O)   :: io
     TYPE(INT_XA2O) :: interp
     TYPE(A2O_transformation),DIMENSION(:,:), POINTER   :: a2o_tr
     REAL(DP), DIMENSION(:,:),   POINTER :: o_2d  => NULL() ! OCEAN
     REAL(DP), DIMENSION(:,:),   POINTER :: a_2d  => NULL() ! ATMOSPHERE
     REAL(DP), DIMENSION(:,:,:), POINTER :: a_3d  => NULL() ! ATMOSPHERE
     REAL(DP), DIMENSION(:,:),   POINTER :: a_tmp => NULL() ! ATMOSPHERE
     LOGICAL                             :: ok = .FALSE.
     INTEGER                             :: rep
  END TYPE T_A2O

  INTEGER, PARAMETER :: NMAXA2O = 100
  TYPE(IO_A2O), DIMENSION(NMAXA2O), SAVE :: A2O
  INTEGER,                          SAVE :: NA2O = 0
  TYPE(INT_A2O),DIMENSION(NMAXA2O), SAVE :: INTA2O
  TYPE(T_A2O),  DIMENSION(NMAXA2O), SAVE :: XA2O
  REAL(DP), POINTER                      :: r_reset_XA2O ! 1 reset, 0 nothing
  REAL(DP),POINTER                       :: dt_a2o ! time

  ! FOR VECTORS
  TYPE IOV_A2O
     CHARACTER(LEN=STRLEN_OBJECT)  :: x_name = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: y_name = ''
  END TYPE IOV_A2O

  TYPE V_A2O
     TYPE(IOV_A2O) :: io
     REAL(DP), DIMENSION(:,:),   POINTER :: x_2d  => NULL() ! x-direction field
     REAL(DP), DIMENSION(:,:),   POINTER :: y_2d  => NULL() ! y-direction field
     LOGICAL                       :: lexist = .FALSE. 
  END TYPE V_A2O

  TYPE(IOV_A2O), DIMENSION(NMAXA2O/2), SAVE :: VEC_A2O
  INTEGER,                             SAVE :: NVEC_A2O = 0
  TYPE(V_A2O), DIMENSION(NMAXA2O/2),   SAVE :: XVEC_A2O

  REAL(DP), DIMENSION(:,:),   POINTER :: vec_tmp_x => NULL() ! TEMPORARY vector
  REAL(DP), DIMENSION(:,:),   POINTER :: vec_tmp_y => NULL() ! TEMPORARY vector


  ! #####################
  ! OCEAN -> ATMOSPHERE
  ! #####################

  TYPE IO_O2A
     CHARACTER(LEN=STRLEN_OBJECT)  :: name = ''
     CHARACTER(LEN=STRLEN_CHANNEL) :: o_channel = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: o_object = ''
     INTEGER                       :: cmethod = A2O_CHLEV
     INTEGER                       :: tmethod = A2O_TAVE
     LOGICAL                       :: lmask = .FALSE.
     REAL(DP)                      :: min_value = -1.E34_dp
     REAL(DP)                      :: max_value = 1.E34_dp
  END TYPE IO_O2A

  ! interpolation specific
  TYPE INT_O2A
     CHARACTER(char_len) :: source_grid = ''
     CHARACTER(char_len) :: dest_grid   = ''
     CHARACTER(char_len) :: map_method  = ''
     CHARACTER(char_len) :: normalize_opt = ''
  END TYPE INT_O2A

  LOGICAL, DIMENSION(4,4,3) :: INTCALCO2A = .FALSE.

  TYPE INT_XO2A
     INTEGER :: source_grid 
     INTEGER :: dest_grid   
     INTEGER :: map_method  
     INTEGER :: norm_opt 
  END TYPE INT_XO2A

  ! TRANSFORMATION INFORMATIONS
  TYPE O2A_transformation
    INTEGER                          :: nn                      ! number of neighbors
    INTEGER, DIMENSION(:),  POINTER  :: om_ilon_g   => NULL()   ! mpiom lon index global
    INTEGER, DIMENSION(:),  POINTER  :: om_ilat_g   => NULL()   ! mpiom lat index global
    REAL(DP), DIMENSION(:), POINTER  :: weight      => NULL()   ! weight per point
    REAL(DP), POINTER                :: frac        => NULL()   ! fractional area per point
  END TYPE O2A_transformation

  TYPE T_O2A
     TYPE(IO_O2A)   :: io
     TYPE(INT_XO2A) :: interp
     TYPE(O2A_transformation),DIMENSION(:,:), POINTER   :: o2a_tr
     REAL(DP), DIMENSION(:,:),   POINTER :: a_2d  => NULL() ! ATMOSPHERE
     REAL(DP), DIMENSION(:,:),   POINTER :: o_2d  => NULL() ! OCEAN
     REAL(DP), DIMENSION(:,:,:), POINTER :: o_3d  => NULL() ! OCEAN
     REAL(DP), DIMENSION(:,:),   POINTER :: o_tmp => NULL() ! OCEAN
     LOGICAL                             :: ok = .FALSE.
     INTEGER                             :: rep 
  END TYPE T_O2A
  
  INTEGER, PARAMETER :: NMAXO2A = 100
  TYPE(IO_O2A), DIMENSION(NMAXO2A), SAVE :: O2A
  INTEGER,                          SAVE :: NO2A = 0
  TYPE(INT_O2A),DIMENSION(NMAXA2O), SAVE :: INTO2A
  TYPE(T_O2A),  DIMENSION(NMAXO2A), SAVE :: XO2A
  REAL(DP), POINTER                      :: r_reset_XO2A ! 1 reset, 0 nothing
  REAL(DP), POINTER                 :: dt_o2a ! time

  ! FOR VECTORS
  TYPE IOV_O2A
     CHARACTER(LEN=STRLEN_OBJECT)  :: x_name = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: y_name = ''
  END TYPE IOV_O2A

  TYPE V_O2A
     TYPE(IOV_O2A) :: io
     REAL(DP), DIMENSION(:,:),   POINTER :: x_2d  => NULL() ! x-direction field
     REAL(DP), DIMENSION(:,:),   POINTER :: y_2d  => NULL() ! y-direction field
     INTEGER                       :: i_x ! index x
     INTEGER                       :: i_y ! index y
     LOGICAL                       :: lexist = .FALSE. 
  END TYPE V_O2A

  TYPE(IOV_O2A), DIMENSION(NMAXO2A/2), SAVE :: VEC_O2A
  INTEGER,                             SAVE :: NVEC_O2A = 0
  TYPE(V_O2A), DIMENSION(NMAXO2A/2),   SAVE :: XVEC_O2A


  ! needed to create a2o channel for the first time
  LOGICAL                     :: lcfirst = .TRUE.

  !----------------------------
  !ATMOSPHERE COUPLING 
  !---------------------------

  REAL(dp), DIMENSION(:,:),   POINTER :: tsw_a2o  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: seaice_a2o  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: siced_a2o  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: sni_a2o  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ocu_a2o  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ocv_a2o  => NULL()

  !----------------------------
  !OCEAN COUPLING 
  !---------------------------
! POINTERS FOR EXTERNAL CHANNEL OBJECTS (replacing USE from messy_main_data_bi)
  ! NEW OBJECT
  REAL(dp), DIMENSION(:,:),   POINTER :: awhea     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: awfre     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: aifre     => NULL()

  !E5VDIFF
  REAL(dp), DIMENSION(:,:),   POINTER :: ahflw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ahfsw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evapw     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: evapi     => NULL()
  !ECHAM5
  REAL(dp), DIMENSION(:,:),   POINTER :: radflxw   => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ssfl       => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ssfc       => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: rsfl       => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: rsfc       => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: seaice     => NULL()


  PUBLIC :: a2o_initialize
  PUBLIC :: a2o_init_memory
  PUBLIC :: a2o_init_coupling
  PUBLIC :: a2o_global_start
  PUBLIC :: a2o_global_end
  PUBLIC :: a2o_free_memory


CONTAINS

  ! ######################################################################
  ! PUBLIC SUBROUTINES
  ! ######################################################################
  
  ! ---------------------------------------------------------------------

  SUBROUTINE a2o_initialize

    ! a2o MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Pozzer Andrea, MPICH, Aug 2007
    
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_channel_bi, ONLY: GP_3D_MPIOM, REPR_UNDEF
!    USE messy_mpiom_e5,        ONLY: DT
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'a2o_initialize'
    INTEGER                 :: iou    ! I/O unit
    INTEGER                 :: status ! error status
    INTEGER                 :: i

    CALL start_message_bi(modstr, 'INITIALISATION', substr)

    LRUNSM = (GP_3D_MPIOM /= REPR_UNDEF)
    IF (.NOT. LRUNSM) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'MPIOM REPRESENTATION IS UNDEFINED.'
          WRITE(*,*) 'SUBMODEL CANNOT BE RUN.'
       END IF
       CALL end_message_bi(modstr, 'INITIALISATION', substr)
       CALL finish(substr)
    END IF
    
    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL a2o_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    CALL p_bcast(l_verbose,        p_io)

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL a2o_read_nml_cpl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

!----------- EVEN INITIALISATION --------------------------

    CALL p_bcast_event (trig_a2o, p_io)
    CALL p_bcast_event (trig_o2a, p_io)

!  NOT needs to exchange data with shorter time step than the
! longest timestep of the two models. HENCE:
!   IF (DT < delta_time) THEN 
!      exchange_time = DT
!   ELSE
!      exchange_time = delta_time
!   ENDIF
!   IF (trig_a2o < exchange_time) then trig_a2o=exchange_time 
!   IF (trig_o2a < exchange_time) then trig_o2a=exchange_time 
!    trig_a2o = io_time_event ( INT(DT_steps), 'steps',TRIG_FIRST, 0)
!    trig_o2a = io_time_event ( INT(DT_steps), 'steps',TRIG_FIRST, 0)

    CALL timer_event_init(ev_trig_a2o, trig_a2o, &
           'a2o computation', 'present')
    CALL timer_event_init(ev_trig_o2a, trig_o2a, &
           'o2a computation', 'present')

!---------DEFINE FIELDS FOR COUPLING (A-->O)---------------------

    IF (p_parallel_io) THEN

       ! ### (1) A2O ###
       ! GET NUMBER OF ENTRIES
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'TRANSFORMATION ATMOSPHERE -> OCEAN:'
       WRITE(*,*) '------------------------------------------------------'
       ! COPY DATA
       NA2O = 1
       DO i=1, NMAXA2O

          IF (TRIM(A2O(i)%name) == '') CYCLE
          XA2O(NA2O)%io%name       = TRIM(A2O(i)%name)

          IF (TRIM(A2O(i)%a_channel) == '') CYCLE
          XA2O(NA2O)%io%a_channel = TRIM(A2O(i)%a_channel)

          IF (TRIM(A2O(i)%a_object) == '') CYCLE
          XA2O(NA2O)%io%a_object  = TRIM(A2O(i)%a_object)

          IF (TRIM(INTA2O(i)%source_grid) == '') CYCLE

          IF (TRIM(INTA2O(i)%dest_grid) == '') CYCLE

          IF (TRIM(INTA2O(i)%map_method) == '') CYCLE

          IF (TRIM(INTA2O(i)%normalize_opt) == '') CYCLE

          WRITE(*,*) '  ', TRIM(XA2O(NA2O)%io%name),' <- ', &
               TRIM(XA2O(NA2O)%io%a_channel), &
               '(',TRIM(XA2O(NA2O)%io%a_object),') ...'

          SELECT CASE (A2O(i)%cmethod)
          CASE(A2O_CSUM)
             WRITE(*,*) '   ... SUM OF THE COLUMN'
          CASE(A2O_CAVE)
             WRITE(*,*) '   ... AVERAGE IN THE COLUMN'
          CASE(A2O_CMAX)
             WRITE(*,*) '   ... MAXIMUM IN THE COLUMN'
          CASE(A2O_CLLEV)
             WRITE(*,*) '   ... LOWEST LEVEL IN THE COLUMN'
          CASE(A2O_CHLEV)
             WRITE(*,*) '   ... HIGHEST LEVEL IN THE COLUMN'
          CASE DEFAULT
             WRITE(*,*) '   ... UNKNOWN TRANSFORMATION METHOD '//&
                  &'... skipping'
             CYCLE
          END SELECT
          XA2O(NA2O)%io%cmethod = A2O(i)%cmethod

          SELECT CASE (A2O(i)%tmethod)
          CASE(A2O_TSUM)
             WRITE(*,*) '   ... SUM IN TIME '
          CASE(A2O_TAVE)
             WRITE(*,*) '   ... AVERAGE IN TIME'
          CASE DEFAULT
             WRITE(*,*) '   ... UNKNOWN TRANSFORMATION METHOD '//&
                  &'... skipping'
             CYCLE
          END SELECT
          XA2O(NA2O)%io%tmethod = A2O(i)%tmethod

          XA2O(NA2O)%io%lmask = A2O(i)%lmask
          IF (XA2O(NA2O)%io%lmask) THEN
             WRITE(*,*) '   ... mask correction : YES'
          ELSE
             WRITE(*,*) '   ... mask correction : NO'
          END IF

          WRITE(*,*)'VALUE FORCE TO BE BETWEEN ', &
                      A2O(i)%min_value,' AND ', A2O(i)%max_value

          XA2O(NA2O)%io%min_value = A2O(i)%min_value
          XA2O(NA2O)%io%max_value = A2O(i)%max_value

          select case(INTA2O(i)%source_grid) 
          case('atm')
            XA2O(NA2O)%interp%source_grid = atm 
            WRITE(*,'(A)', ADVANCE='NO') '  ATM -->'
          case('oces')
            XA2O(NA2O)%interp%source_grid = oces 
            WRITE(*,*) '  wrong source grid (Atmosphere to Ocean?)'
            CALL finish(substr)
          case('oceu')
            XA2O(NA2O)%interp%source_grid = oceu 
            WRITE(*,*) '  wrong source grid (Atmosphere to Ocean?)'
            CALL finish(substr)
          case('ocev')
            XA2O(NA2O)%interp%source_grid = ocev 
            WRITE(*,*) '  wrong source grid (Atmosphere to Ocean?)'
            CALL finish(substr)
          case default
          WRITE(*,*) 'unknown source grid'
          CALL finish(substr)
          end select

          select case(INTA2O(i)%dest_grid) 
          case('atm')
            XA2O(NA2O)%interp%dest_grid = atm 
            WRITE(*,*) ' ATM'
            WRITE(*,*) '  wrong destination grid (Atmosphere to Ocean?)'
            CALL finish(substr)
          case('oces')
            XA2O(NA2O)%interp%dest_grid = oces 
            WRITE(*,*) ' OCES'
          case('oceu')
            XA2O(NA2O)%interp%dest_grid = oceu 
            WRITE(*,*) ' OCEU'
          case('ocev')
            XA2O(NA2O)%interp%dest_grid = ocev 
            WRITE(*,*) ' OCEV'
          case default
          WRITE(*,*) 'unknown destination grid'
          CALL finish(substr)
          end select

          WRITE(*,*)'INTERPOLATION  METHOD: '
          select case(INTA2O(i)%map_method)
          case ('conservative')
            XA2O(NA2O)%interp%map_method = map_type_conserv
            WRITE(*,*) ' conservative'
          case ('bilinear')
            XA2O(NA2O)%interp%map_method = map_type_bilinear
            WRITE(*,*) ' bilinear'
          case ('bicubic')
            XA2O(NA2O)%interp%map_method = map_type_bicubic
            WRITE(*,*) ' bicubic'
          case ('distwgt')
            XA2O(NA2O)%interp%map_method = map_type_distwgt
            WRITE(*,*) ' distance weighted'
          case default
            WRITE(*,*) 'unknown mapping method'
            CALL finish(substr)
          end select

          WRITE(*,*)'NORMALIZATION OPTION: '
          select case(INTA2O(i)%normalize_opt(1:4))
          case ('none')
            XA2O(NA2O)%interp%norm_opt = norm_opt_none
            WRITE(*,*) ' none  '
          case ('frac')
            XA2O(NA2O)%interp%norm_opt = norm_opt_frcarea
            WRITE(*,*) ' fractional area  '
          case ('dest')
            XA2O(NA2O)%interp%norm_opt = norm_opt_dstarea
            WRITE(*,*) ' destination area  '
          case default
            WRITE(*,*) 'unknown normalization option'
            CALL finish(substr)
          end select

          ! NEXT ENTRY
          NA2O = NA2O + 1
          WRITE(*,*) '------------------------------------------------------'

       END DO
       NA2O = NA2O - 1

       ! ### (2) VEC_A2O ###
       ! GET NUMBER OF ENTRIES
       NVEC_A2O = 1
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'VECTOR IN TRANSFORMATION ATMOSPHERE -> OCEAN:'
       WRITE(*,*) '------------------------------------------------------'
       DO i=1, NMAXA2O/2
          IF (TRIM(VEC_A2O(i)%x_name) == '') CYCLE
          XVEC_A2O(NVEC_A2O)%io%x_name = TRIM(VEC_A2O(i)%x_name)
          IF (TRIM(VEC_A2O(i)%y_name) == '') CYCLE
          XVEC_A2O(NVEC_A2O)%io%y_name = TRIM(VEC_A2O(i)%y_name)
          WRITE(*,*) 'X-direction:', TRIM(XVEC_A2O(NVEC_A2O)%io%x_name)
          WRITE(*,*) 'Y-direction:', TRIM(XVEC_A2O(NVEC_A2O)%io%y_name)
          NVEC_A2O = NVEC_A2O + 1
          WRITE(*,*) '------------------------------------------------------'
       END DO
       NVEC_A2O = NVEC_A2O - 1

       ! ### (3) O2A ###
       ! GET NUMBER OF ENTRIES
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'TRANSFORMATION OCEAN -> ATMOSPHERE:'
       WRITE(*,*) '------------------------------------------------------'
       ! COPY DATA
       NO2A = 1
       DO i=1, NMAXO2A

          IF (TRIM(O2A(i)%name) == '') CYCLE
          XO2A(NO2A)%io%name       = TRIM(O2A(i)%name)

          IF (TRIM(O2A(i)%o_channel) == '') CYCLE
          XO2A(NO2A)%io%o_channel = TRIM(O2A(i)%o_channel)

          IF (TRIM(O2A(i)%o_object) == '') CYCLE
          XO2A(NO2A)%io%o_object  = TRIM(O2A(i)%o_object)

          IF (TRIM(INTO2A(i)%source_grid) == '') CYCLE

          IF (TRIM(INTO2A(i)%dest_grid) == '') CYCLE

          IF (TRIM(INTO2A(i)%map_method) == '') CYCLE

          IF (TRIM(INTO2A(i)%normalize_opt) == '') CYCLE

          WRITE(*,*) '  ', TRIM(XO2A(NO2A)%io%name),' <- ', &
               TRIM(XO2A(NO2A)%io%o_channel), &
               '(',TRIM(XO2A(NO2A)%io%o_object),') ...'

          SELECT CASE (O2A(i)%cmethod)
          CASE(A2O_CSUM)
             WRITE(*,*) '   ... SUM OF THE COLUMN'
          CASE(A2O_CAVE)
             WRITE(*,*) '   ... AVERAGE IN THE COLUMN'
          CASE(A2O_CMAX)
             WRITE(*,*) '   ... MAXIMUM IN THE COLUMN'
          CASE(A2O_CLLEV)
             WRITE(*,*) '   ... LOWEST LEVEL IN THE COLUMN'
          CASE(A2O_CHLEV)
             WRITE(*,*) '   ... HIGHEST LEVEL IN THE COLUMN'
          CASE DEFAULT
             WRITE(*,*) '   ... UNKNOWN TRANSFORMATION METHOD '//&
                  &'... skipping'
             CYCLE
          END SELECT
          XO2A(NO2A)%io%cmethod = O2A(i)%cmethod

          SELECT CASE (O2A(i)%tmethod)
          CASE(A2O_TSUM)
             WRITE(*,*) '   ... SUM IN TIME '
          CASE(A2O_TAVE)
             WRITE(*,*) '   ... AVERAGE IN TIME'
          CASE DEFAULT
             WRITE(*,*) '   ... UNKNOWN TRANSFORMATION METHOD '//&
                  &'... skipping'
             CYCLE
          END SELECT
          XO2A(NO2A)%io%tmethod = O2A(i)%tmethod
          
          XO2A(NO2A)%io%lmask = O2A(i)%lmask
          IF (XO2A(NO2A)%io%lmask) THEN
             WRITE(*,*) '   ... mask correction : YES'
          ELSE
             WRITE(*,*) '   ... mask correction : NO'
          END IF

          WRITE(*,*) ' VALUE FORCE TO BE BETWEEN ', &
                      O2A(i)%min_value,' AND ', O2A(i)%max_value

          XO2A(NO2A)%io%min_value = O2A(i)%min_value
          XO2A(NO2A)%io%max_value = O2A(i)%max_value

          select case(INTO2A(i)%source_grid) 
          case('atm')
            XO2A(NO2A)%interp%source_grid = atm 
            WRITE(*,*) '  wrong source grid (Ocean to Atmosphere?)'
            CALL finish(substr)
          case('oces')
            XO2A(NO2A)%interp%source_grid = oces 
            WRITE(*,'(A)', ADVANCE='NO') '  OCES -->'
          case('oceu')
            XO2A(NO2A)%interp%source_grid = oceu 
            WRITE(*,'(A)', ADVANCE='NO') '  OCEU -->'
          case('ocev')
            XO2A(NO2A)%interp%source_grid = ocev 
            WRITE(*,'(A)', ADVANCE='NO') '  OCEV -->'
          case default
          WRITE(*,*) 'unknown source grid'
          CALL finish(substr)
          end select

          select case(INTO2A(i)%dest_grid) 
          case('atm')
            XO2A(NO2A)%interp%dest_grid = atm 
            WRITE(*,*) ' ATM'
          case('oces')
            XO2A(NO2A)%interp%dest_grid = oces 
            WRITE(*,*) '  wrong destination grid (Ocean to Atmosphere?)'
            CALL finish(substr)
          case('oceu')
            XO2A(NO2A)%interp%dest_grid = oceu 
            WRITE(*,*) '  wrong destination grid (Ocean to Atmosphere?)'
            CALL finish(substr)
          case('ocev')
            XO2A(NO2A)%interp%dest_grid = ocev 
            WRITE(*,*) '  wrong destination grid (Ocean to Atmosphere?)'
            CALL finish(substr)
          case default
          WRITE(*,*) 'unknown destination grid'
          CALL finish(substr)
          end select

          WRITE(*,*)'INTERPOLATION  METHOD: '
          select case(INTO2A(i)%map_method)
          case ('conservative')
            XO2A(NO2A)%interp%map_method = map_type_conserv
            WRITE(*,*) ' conservative'
          case ('bilinear')
            XO2A(NO2A)%interp%map_method = map_type_bilinear
            WRITE(*,*) ' bilinear'
          case ('bicubic')
            XO2A(NO2A)%interp%map_method = map_type_bicubic
            WRITE(*,*) ' bicubic'
          case ('distwgt')
            XO2A(NO2A)%interp%map_method = map_type_distwgt
            WRITE(*,*) ' distance weighted'
          case default
            WRITE(*,*) 'unknown mapping method'
            CALL finish(substr)
          end select

          WRITE(*,*)'NORMALIZATION OPTION: '
          select case(INTO2A(i)%normalize_opt(1:4))
          case ('none')
            XO2A(NO2A)%interp%norm_opt = norm_opt_none
            WRITE(*,*) ' none  '
          case ('frac')
            XO2A(NO2A)%interp%norm_opt = norm_opt_frcarea
            WRITE(*,*) ' fractional area  '
          case ('dest')
            XO2A(NO2A)%interp%norm_opt = norm_opt_dstarea
            WRITE(*,*) ' destination area  '
          case default
            WRITE(*,*) 'unknown normalization option'
            CALL finish(substr)
          end select

          ! NEXT ENTRY
          NO2A = NO2A + 1
          WRITE(*,*) '------------------------------------------------------'

       END DO
       NO2A = NO2A - 1

       ! ### (4) VEC_O2A ###
       ! GET NUMBER OF ENTRIES
       NVEC_O2A = 1
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'VECTOR IN TRANSFORMATION OCEAN --> ATMOSPHERE:'
       WRITE(*,*) '------------------------------------------------------'
       DO i=1, NMAXO2A/2
          IF (TRIM(VEC_O2A(i)%x_name) == '') CYCLE
          XVEC_O2A(NVEC_O2A)%io%x_name = TRIM(VEC_O2A(i)%x_name)
          IF (TRIM(VEC_O2A(i)%y_name) == '') CYCLE
          XVEC_O2A(NVEC_O2A)%io%y_name = TRIM(VEC_O2A(i)%y_name)
          WRITE(*,*) 'X-direction:', TRIM(XVEC_O2A(NVEC_O2A)%io%x_name)
          WRITE(*,*) 'Y-direction:', TRIM(XVEC_O2A(NVEC_O2A)%io%y_name)
          NVEC_O2A = NVEC_O2A + 1
          WRITE(*,*) '------------------------------------------------------'
       END DO
       NVEC_O2A = NVEC_O2A - 1

    ENDIF

    ! BROADCAST RESULTS

    CALL p_bcast(NA2O, p_io)
    CALL p_bcast(NVEC_A2O, p_io)
    CALL p_bcast(NO2A, p_io)
    CALL p_bcast(NVEC_O2A, p_io)

    CALL p_bcast(num_neighbors, p_io)
    CALL p_bcast(north_thresh,p_io)
    CALL p_bcast(south_thresh,p_io)
    CALL p_bcast(max_subseg,p_io) 
    CALL p_bcast(max_iter,p_io) 
    CALL p_bcast(converge,p_io)    

    DO i=1, NA2O
       CALL p_bcast(XA2O(i)%io%name, p_io)
       CALL p_bcast(XA2O(i)%io%a_channel, p_io)
       CALL p_bcast(XA2O(i)%io%a_object, p_io)
       CALL p_bcast(XA2O(i)%io%cmethod, p_io)
       CALL p_bcast(XA2O(i)%io%tmethod, p_io)
       CALL p_bcast(XA2O(i)%io%lmask, p_io)
       CALL p_bcast(XA2O(i)%io%max_value, p_io)
       CALL p_bcast(XA2O(i)%io%min_value, p_io)
    END DO

    DO i=1, NA2O
        CALL p_bcast(XA2O(i)%interp%source_grid,p_io)
        CALL p_bcast(XA2O(i)%interp%dest_grid,p_io)
        CALL p_bcast(XA2O(i)%interp%map_method,p_io)
        CALL p_bcast(XA2O(i)%interp%norm_opt,p_io)
    ENDDO

    DO i=1, NVEC_A2O
       CALL p_bcast(XVEC_A2O(i)%io%x_name, p_io)
       CALL p_bcast(XVEC_A2O(i)%io%y_name, p_io)
    END DO

    DO i=1, NO2A 
       CALL p_bcast(XO2A(i)%io%name, p_io)
       CALL p_bcast(XO2A(i)%io%o_channel, p_io)
       CALL p_bcast(XO2A(i)%io%o_object, p_io)
       CALL p_bcast(XO2A(i)%io%cmethod, p_io)
       CALL p_bcast(XO2A(i)%io%tmethod, p_io)
       CALL p_bcast(XO2A(i)%io%lmask, p_io)
       CALL p_bcast(XO2A(i)%io%max_value, p_io)
       CALL p_bcast(XO2A(i)%io%min_value, p_io)
    END DO

    DO i=1, NO2A
        CALL p_bcast(XO2A(i)%interp%source_grid,p_io)
        CALL p_bcast(XO2A(i)%interp%dest_grid,p_io)
        CALL p_bcast(XO2A(i)%interp%map_method,p_io)
        CALL p_bcast(XO2A(i)%interp%norm_opt,p_io)
    ENDDO
    DO i=1, NVEC_O2A
       CALL p_bcast(XVEC_O2A(i)%io%x_name, p_io)
       CALL p_bcast(XVEC_O2A(i)%io%y_name, p_io)
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NA2O+NO2A,' TRANSFORMATIONS(S) INITIALIZED !'
       WRITE(*,*) ' ---> ',NVEC_A2O+NVEC_O2A,' VECTOR TRANSFORMATIONS(S) INITIALIZED !'
    END IF

    LRUNSM = (NA2O + NO2A > 0)

    ! FINALLY CALCULATE WHICH INTERPOLATION ROUTINES MUST BE CALCULATED
    DO i=1, NA2O
       INTCALCA2O(XA2O(i)%interp%dest_grid,       &
                  XA2O(i)%interp%map_method,XA2O(i)%interp%norm_opt) =.TRUE.
    ENDDO
    DO i=1, NO2A
       INTCALCO2A(XO2A(i)%interp%source_grid,       &
                  XO2A(i)%interp%map_method,XO2A(i)%interp%norm_opt) =.TRUE.
    ENDDO

    CALL end_message_bi(modstr, 'INITIALISATION', substr)

  END SUBROUTINE a2o_initialize

  ! ---------------------------------------------------------------------
  SUBROUTINE a2o_init_memory

! op_pj_20160613+
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL, GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'a2o_init_memory'
    INTEGER                     :: status

    IF (lcfirst) THEN
       lcfirst = .FALSE.
       CALL new_channel(status, modstr, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
    ENDIF
    ! ------------------------------------------------------------------------
    ! ... required for A2O/MPIOM (a2o.nml), HD
    CALL new_channel_object(status, modstr,  'awhea', &
         p2=awhea, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awhea', &
         'long_name', c='net heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awhea', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ... required for A2O/MPIOM (a2o.nml), HD
    CALL new_channel_object(status, modstr,  'awfre', &
         p2=awfre, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awfre', &
         'long_name', c='water flux into ocean')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awfre', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ... required for A2O/MPIOM (a2o.nml), HD
    CALL new_channel_object(status, modstr,  'aifre', &
         p2=aifre, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aifre', &
         'long_name', c='downward snow flux where sea-ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aifre', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------

  END SUBROUTINE a2o_init_memory
  ! ---------------------------------------------------------------------

  SUBROUTINE a2o_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, message,finish
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_3D_1LEV,     & 
                                           GP_2D_HORIZONTAL,GP_3D_INT,&
                                           GP_3D_MPIOM, GP_2D_MPIOM,  &
                                           SCALAR
    ! MESSy
    USE messy_main_channel,         ONLY: new_channel, new_channel_object &
                                        , new_attribute, get_attribute    &
                                        , get_channel_object              &
                                        , get_channel_object_info         &
                                        , get_channel_info
    USE messy_main_constants_mem,   ONLY: STRLEN_ULONG
    USE messy_main_timer,           ONLY: delta_time, lstart

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    LOGICAL                     :: lfirst_A2O = .TRUE.
    LOGICAL                     :: lfirst_O2A = .TRUE.
    CHARACTER(LEN=*), PARAMETER :: substr = 'a2o_init_coupling'
    INTEGER                     :: status
    INTEGER                     :: i,j,n,ii,l,k,nn
    INTEGER                     :: reprid
    CHARACTER(LEN=STRLEN_ULONG) :: unit

    IF (.NOT.LRUNSM) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! (1) Atmosphere -> Ocean
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'TRANSFORMATION ATMOSPHERE -> OCEAN:'
       WRITE(*,*) '------------------------------------------------------'
    END IF

    loop_a2o: DO i=1, NA2O

       IF (XA2O(i)%ok) CYCLE

       CALL get_channel_object_info(status &
            , TRIM(XA2O(i)%io%a_channel), TRIM(XA2O(i)%io%a_object) &
            , reprid=reprid)
       XA2O(i)%rep = reprid    

       IF (XA2O(i)%rep == GP_3D_MID .or. XA2O(i)%rep == GP_3D_1LEV) THEN
          ! CHECK AVAILABILITY OF OCEAN CHANNEL/OBJECT
          CALL get_channel_object(status &
               , TRIM(XA2O(i)%io%a_channel), TRIM(XA2O(i)%io%a_object) &
               , p3=XA2O(i)%a_3d &
               )
          IF (status /= 0) THEN
             CALL message('  ',TRIM(XA2O(i)%io%a_channel)//' - '//&
                  &TRIM(XA2O(i)%io%a_object)//' not found ... skipping' )
             CYCLE
          ELSE
             CALL message('  ',TRIM(XA2O(i)%io%a_channel)//' - '//&
                  &TRIM(XA2O(i)%io%a_object)//' found ...' )
             CALL message('  ',' ... GP_3D_MID/GP_3D_1LEV representation ... ')
          ENDIF
       ELSE IF(XA2O(i)%rep == GP_3D_INT) THEN
          ! CHECK AVAILABILITY OF OCEAN CHANNEL/OBJECT
          CALL get_channel_object(status &
               , TRIM(XA2O(i)%io%a_channel), TRIM(XA2O(i)%io%a_object) &
               , p3=XA2O(i)%a_3d &
               )
          IF (status /= 0) THEN
             CALL message('  ',TRIM(XA2O(i)%io%a_channel)//' - '//&
                  &TRIM(XA2O(i)%io%a_object)//' not found ... skipping' )
             CYCLE
          ELSE
             CALL message('  ',TRIM(XA2O(i)%io%a_channel)//' - '//&
                  &TRIM(XA2O(i)%io%a_object)//' found ...' )
             CALL message('  ',' ... GP_3D_INT representation ... ')
          ENDIF
       ELSE IF(XA2O(i)%rep == GP_2D_HORIZONTAL) THEN
          ! CHECK AVAILABILITY OF OCEAN CHANNEL/OBJECT
          CALL get_channel_object(status &
               , TRIM(XA2O(i)%io%a_channel), TRIM(XA2O(i)%io%a_object) &
               , p2=XA2O(i)%a_2d &
               )
          IF (status /= 0) THEN
             CALL message('  ',TRIM(XA2O(i)%io%a_channel)//' - '//&
                  &TRIM(XA2O(i)%io%a_object)//' not found ... skipping' )
             CYCLE
          ELSE
             CALL message('  ',TRIM(XA2O(i)%io%a_channel)//' - '//&
                  &TRIM(XA2O(i)%io%a_object)//' found ...' )
             CALL message('  ',' ... GP_2D_HORIZONTAL representation ... ')
          ENDIF
       ELSE
          CALL message('  ',' ... wrong representation ... skipping')
          CYCLE
       ENDIF

       ! EVERYTHIG IS OK NOW
       XA2O(i)%ok = .TRUE.

       ! WE HAVE AT LEAST 1 FIELD PRESENT : 
       ! WE MUST RUN THE A2O TRANSFORMATION
       LRUNSM_A2O =.TRUE.

       ! CREATE CHANNEL AND OBJECTS FOR THE FIRST OBJECT
       IF (lcfirst) THEN
          lcfirst = .FALSE.
          CALL new_channel(status, modstr)
          CALL channel_halt(substr, status)
       ENDIF

       IF (lfirst_A2O) THEN
          lfirst_A2O = .FALSE.
          CALL new_channel_object(status, modstr, 'r_reset_XA2O' &
               ,p0=r_reset_XA2O, REPRID=SCALAR, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          r_reset_XA2O = 0.0_dp
          CALL new_channel_object(status, modstr, 'dt_a2o' &
               ,p0=dt_a2o, REPRID=SCALAR, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
       ENDIF

       ! CREATE NEW CHANNEL OBJECT
       CALL new_channel_object(status, modstr    &
            , TRIM(XA2O(i)%io%name), p2=XA2O(i)%o_2d, & 
            reprid=GP_2D_MPIOM, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)

       ! ALLOCATE FILEDS FOR TEMPORARY STORING
       CALL new_channel_object(status, modstr    &
            , TRIM(XA2O(i)%io%name)//'_tmp', p2=XA2O(i)%a_tmp, & 
            reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)

       ! ADD ATTRIBUTES
       ! ... UNITS
       CALL get_attribute(status &
            , TRIM(XA2O(i)%io%a_channel), TRIM(XA2O(i)%io%a_object) &
            , 'units', c=unit)
       IF (status == 0) THEN
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'units',c = TRIM(unit) )
       END IF
       ! ... ORIGIN
       CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
            , 'long_name', c='A2O transformation of '//&
            &TRIM(XA2O(i)%io%a_channel)//' - '&
            &//TRIM(XA2O(i)%io%a_object) )
       CALL channel_halt(substr, status)
       ! ... COLUMN  METHOD
       SELECT CASE(XA2O(i)%io%cmethod)
       CASE(A2O_CSUM)
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'column_method', c='SUM')
       CASE(A2O_CAVE)
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'column_method', c='AVERAGE')
       CASE(A2O_CMAX)
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'column_method', c='MAX')
       CASE(A2O_CHLEV)
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'column_method', c='HIGHEST LEV')
       CASE(A2O_CLLEV)
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'column_method', c='LOWEST LEV')
       END SELECT
       CALL channel_halt(substr, status)
       ! ... TIME  METHOD
       SELECT CASE(XA2O(i)%io%tmethod)
       CASE(A2O_TSUM)
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'time_method', c='SUM')
       CASE(A2O_TAVE)
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'time_method', c='AVERAGE')
       END SELECT
       CALL channel_halt(substr, status)

       ! ... MASS CONSERVATION
       IF (XA2O(i)%io%lmask) THEN
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'mask_correction', c='YES')
       ELSE
          CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
               , 'mask_correction', c='NO')
       END IF
       CALL channel_halt(substr, status)

       ! ... MAX-MIN VALUE
       CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
            , 'max_value', r=XA2O(i)%io%max_value)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(XA2O(i)%io%name) &
            , 'min_value', r=XA2O(i)%io%min_value)
       CALL channel_halt(substr, status)
       
    END DO loop_a2o


    ! (2) Atmosphere vector -> Ocean vector
    loop_vec_a2o: DO i=1, NVEC_A2O
          CALL get_channel_object(status &
               , modstr, TRIM(XVEC_A2O(i)%io%x_name), p2=XVEC_A2O(i)%x_2d )
          IF (status /= 0) THEN
             CALL message('  ','vector component x ('//TRIM(XVEC_A2O(i)%io%x_name)//') &
                &not found ... skipping' )
             CYCLE
          ELSE 
             CALL message('  ','vector component x ('//TRIM(XVEC_A2O(i)%io%x_name)//') &
                &found !' )
          ENDIF
          CALL get_channel_object(status &
               , modstr, TRIM(XVEC_A2O(i)%io%y_name), p2=XVEC_A2O(i)%y_2d )
          IF (status /= 0) THEN
             CALL message('  ','vector component y ('//TRIM(XVEC_A2O(i)%io%y_name)//') &
                &not found ... skipping' )
             CYCLE
          ELSE
             CALL message('  ','vector component y ('//TRIM(XVEC_A2O(i)%io%y_name)//') &
                &found !' )
          END IF
          XVEC_A2O(i)%lexist=.TRUE.
    ENDDO loop_vec_a2o

    ! (3) Ocean -> Atmosphere
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'TRANSFORMATION OCEAN -> ATMOSPHERE:'
       WRITE(*,*) '------------------------------------------------------'
    END IF

    loop_o2a: DO i=1, NO2A

       IF (XO2A(i)%ok) CYCLE

       CALL get_channel_object_info(status &
            , TRIM(XO2A(i)%io%o_channel), TRIM(XO2A(i)%io%o_object) &
            , reprid=reprid)
       XO2A(i)%rep = reprid    

       IF (XO2A(i)%rep == GP_3D_MPIOM) THEN
          ! CHECK AVAILABILITY OF OCEAN CHANNEL/OBJECT
          CALL get_channel_object(status &
               , TRIM(XO2A(i)%io%o_channel), TRIM(XO2A(i)%io%o_object) &
               , p3=XO2A(i)%o_3d &
               )
          IF (status /= 0) THEN
             CALL message('  ',TRIM(XO2A(i)%io%o_channel)//' - '//&
                  &TRIM(XO2A(i)%io%o_object)//' not found ... skipping' )
             CYCLE
          ELSE
             CALL message('  ',TRIM(XO2A(i)%io%o_channel)//' - '//&
                  &TRIM(XO2A(i)%io%o_object)//' found ...' )
             CALL message('  ',' ... GP_3D_MPIOM representation ... ')
          ENDIF
       ELSE IF (XO2A(i)%rep == GP_2D_MPIOM) THEN
          ! CHECK AVAILABILITY OF OCEAN CHANNEL/OBJECT
          CALL get_channel_object(status &
               , TRIM(XO2A(i)%io%o_channel), TRIM(XO2A(i)%io%o_object) &
               , p2=XO2A(i)%o_2d &
               )
          IF (status /= 0) THEN
             CALL message('  ',TRIM(XO2A(i)%io%o_channel)//' - '//&
                  &TRIM(XO2A(i)%io%o_object)//' not found ... skipping' )
             CYCLE
          ELSE
             CALL message('  ',TRIM(XO2A(i)%io%o_channel)//' - '//&
                  &TRIM(XO2A(i)%io%o_object)//' found ...' )
             CALL message('  ',' ... GP_2D_MPIOM representation ... ')
          ENDIF
       ELSE
          CALL message('  ',' ... wrong representation ... skipping')
          CYCLE
       END IF

       ! EVERYTHIG IS OK NOW
       XO2A(i)%ok = .TRUE.

       ! WE HAVE AT LEAST 1 FIELD PRESENT : 
       ! WE MUST RUN THE O2A TRANSFORMATION
       LRUNSM_O2A =.TRUE.

       ! CREATE CHANNEL AND OBJECTS FOR THE FIRST OBJECT
       IF (lcfirst) THEN
          lcfirst = .FALSE.
          CALL new_channel(status, modstr)
          CALL channel_halt(substr, status)
       ENDIF

       IF (lfirst_O2A) THEN
          lfirst_O2A = .FALSE.
          CALL new_channel_object(status, modstr, 'r_reset_XO2A'    &
               ,p0=r_reset_XO2A, REPRID=SCALAR, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          r_reset_XO2A = 0.0_dp
          CALL new_channel_object(status, modstr, 'dt_o2a' &
               ,p0=dt_o2a, REPRID=SCALAR, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
       ENDIF

       ! CREATE NEW CHANNEL OBJECT
       CALL new_channel_object(status, modstr    &
            , TRIM(XO2A(i)%io%name), p2=XO2A(i)%a_2d, & 
            reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)

       ! ALLOCATE FILEDS FOR TEMPORARY STORING
       CALL new_channel_object(status, modstr    &
            , TRIM(XO2A(i)%io%name)//'_tmp', p2=XO2A(i)%o_tmp, & 
            reprid=GP_2D_MPIOM, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)

       ! ADD ATTRIBUTES
       ! ... UNITS
       CALL get_attribute(status &
            , TRIM(XO2A(i)%io%o_channel), TRIM(XO2A(i)%io%o_object) &
            , 'units', c=unit)
       IF (status == 0) THEN
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'units',c = TRIM(unit) )
       END IF
       ! ... ORIGIN
       CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
            , 'long_name', c='O2A transformation of '//&
            &TRIM(XO2A(i)%io%o_channel)//' - '&
            &//TRIM(XO2A(i)%io%o_object) )
       CALL channel_halt(substr, status)
       ! ... COLUMN METHOD
       SELECT CASE(XO2A(i)%io%cmethod)
       CASE(A2O_CSUM)
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'column_method', c='SUM')
       CASE(A2O_CAVE)
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'column_method', c='AVERAGE')
       CASE(A2O_CMAX)
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'column_method', c='MAX')
       CASE(A2O_CHLEV)
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'column_method', c='HIGHEST LEV')
       CASE(A2O_CLLEV)
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'column_method', c='LOWEST LEV')
       END SELECT
       CALL channel_halt(substr, status)
       ! ... TIME METHOD
       SELECT CASE(XO2A(i)%io%tmethod)
       CASE(A2O_TSUM)
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'time_method', c='SUM')
       CASE(A2O_TAVE)
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'time_method', c='AVERAGE')
       END SELECT
       CALL channel_halt(substr, status)
       ! ... mask CORRECTION
       IF (XO2A(i)%io%lmask) THEN
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'mask_correction', c='YES')
       ELSE
          CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
               , 'mask_correction', c='NO')
       END IF
       CALL channel_halt(substr, status)

       ! ... MAX-MIN VALUE
       CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
            , 'max_value', r=XO2A(i)%io%max_value)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(XO2A(i)%io%name) &
            , 'min_value', r=XO2A(i)%io%min_value)
       CALL channel_halt(substr, status)
       
    END DO loop_o2a

    ! (4) Ocean vector -> Atmosphere vector
    loop_vec_o2a: DO i=1, NVEC_O2A
          CALL get_channel_object(status &
               , modstr, TRIM(XVEC_O2A(i)%io%x_name), p2=XVEC_O2A(i)%x_2d )
          IF (status /= 0) THEN
             CALL message('  ','vector component x ('//TRIM(XVEC_O2A(i)%io%x_name)//') &
                &not found ... skipping' )
             CYCLE
          ELSE 
             CALL message('  ','vector component x ('//TRIM(XVEC_O2A(i)%io%x_name)//') &
                &found !' )
          ENDIF
          CALL get_channel_object(status &
               , modstr, TRIM(XVEC_O2A(i)%io%y_name), p2=XVEC_O2A(i)%y_2d )
          IF (status /= 0) THEN
             CALL message('  ','vector component y ('//TRIM(XVEC_O2A(i)%io%y_name)//') &
                &not found ... skipping' )
             CYCLE
          ELSE 
             CALL message('  ','vector component y ('//TRIM(XVEC_O2A(i)%io%y_name)//') &
                &found !' )
          END IF
          XVEC_O2A(i)%lexist=.TRUE.
          DO ii=1, NO2A
            IF (XO2A(ii)%io%name == TRIM(XVEC_O2A(i)%io%x_name)) XVEC_O2A(i)%i_x = ii
            IF (XO2A(ii)%io%name == TRIM(XVEC_O2A(i)%io%y_name)) XVEC_O2A(i)%i_y = ii
          ENDDO
    ENDDO loop_vec_o2a

! here we create channel needed for the coupling atmosphere to ocean
    IF(L_COUPLING) THEN
         CALL get_channel_object(status          &
              , oname='ahflw', cname='e5vdiff'  & 
              ,p2=ahflw)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING ahflw from e5vdiff... trying with vertex"
            CALL get_channel_object(status          &
                  , oname='ahflw', cname='vertex'  & 
                  ,p2=ahflw)
            IF (status /= 0) THEN 
               WRITE(*,*) "NO channel ahflw from e5vdiff or vertex"
               CALL finish(substr)
            ENDIF
         ENDIF
         CALL get_channel_object(status          &
              , oname='ahfsw', cname='e5vdiff'  & 
              ,p2=ahfsw)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING ahfsw from e5vdiff... trying with vertex"
            CALL get_channel_object(status          &
                  , oname='ahfsw', cname='vertex'  & 
                  ,p2=ahfsw)
            IF (status /= 0) THEN 
               WRITE(*,*) "NO channel ahflw from e5vdiff or vertex"
               CALL finish(substr)
            ENDIF
         ENDIF
         CALL get_channel_object(status          &
              , oname='evapw', cname='e5vdiff'  & 
              ,p2=evapw)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING evapw from e5vdiff... trying with vertex"
            CALL get_channel_object(status          &
                  , oname='evapw', cname='vertex'  & 
                  ,p2=evapw)
            IF (status /= 0) THEN 
               WRITE(*,*) "NO channel ahflw from e5vdiff or vertex"
               CALL finish(substr)
            ENDIF
         ENDIF
         CALL get_channel_object(status          &
              , oname='evapi', cname='e5vdiff'  & 
              ,p2=evapi)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING evapi from e5vdiff... trying with vertex"
            CALL get_channel_object(status          &
                  , oname='evapi', cname='vertex'  & 
                  ,p2=evapi)
            IF (status /= 0) THEN 
               WRITE(*,*) "NO channel ahflw from e5vdiff or vertex"
               CALL finish(substr)
            ENDIF
         ENDIF
         CALL get_channel_object(status          &
              , oname='radflxw', cname='ECHAM5'  & 
              ,p2=radflxw)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING radflxw from ECHAM5"
            CALL finish(substr)
         ENDIF
         CALL get_channel_object(status          &
              , oname='rsfl_2d', cname='ECHAM5'  & 
              ,p2=rsfl)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING rsfl_2d from ECHAM5"
            CALL finish(substr)
         ENDIF
         CALL get_channel_object(status          &
              , oname='rsfc_2d', cname='ECHAM5'  & 
              ,p2=rsfc)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING ssfl_2d from ECHAM5"
            CALL finish(substr)
         ENDIF
         CALL get_channel_object(status          &
              , oname='ssfl_2d', cname='ECHAM5'  & 
              ,p2=ssfl)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING ssfl_2d from ECHAM5"
            CALL finish(substr)
         ENDIF
         CALL get_channel_object(status          &
              , oname='ssfc_2d', cname='ECHAM5'  & 
              ,p2=ssfc)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING ssfc_2d from ECHAM5"
            CALL finish(substr)
         ENDIF
         CALL get_channel_object(status          &
              , oname='seaice', cname='ECHAM5'  & 
              ,p2=seaice)
         IF (status /= 0) THEN
            WRITE(*,*) "MISSING seaice from ECHAM5"
            CALL finish(substr)
         ENDIF
    ENDIF
! here we get channel needed for the coupling ocean to atmosphere
    IF(lcouple) THEN
         ! sst
          CALL get_channel_object(status          &
               , oname='tho', cname='a2o'  & 
               ,p2=tsw_a2o)
         IF (status /= 0) THEN
         WRITE(*,*) "MISSING tho from a2o"
         CALL finish(substr)
         ENDIF
         ! sea ice thikness
          CALL get_channel_object(status          &
               , oname='sictho', cname='a2o'  & 
               ,p2=siced_a2o)
         IF (status /= 0) THEN
         WRITE(*,*) "MISSING sictho from a2o"
         CALL finish(substr)
         ENDIF
         ! ice compactness
          CALL get_channel_object(status          &
               , oname='sicomo', cname='a2o'  & 
               ,p2=seaice_a2o)
         IF (status /= 0) THEN
         WRITE(*,*) "MISSING sicomo from a2o"
         CALL finish(substr)
         ENDIF
         ! snow thikness
          CALL get_channel_object(status          &
               , oname='sicsno', cname='a2o'  & 
               ,p2=sni_a2o)
         IF (status /= 0) THEN
         WRITE(*,*) "MISSING sicsno from a2o"
         CALL finish(substr)
         ENDIF
         ! socu ocean westward velocity
          CALL get_channel_object(status          &
               , oname='socu', cname='a2o'  & 
               ,p2=ocu_a2o)
         IF (status /= 0) THEN
         WRITE(*,*) "MISSING socu from a2o"
         CALL finish(substr)
         ENDIF
         ! socv ocean eastward velocity
          CALL get_channel_object(status          &
               , oname='socv', cname='a2o'  & 
               ,p2=ocv_a2o)
         IF (status /= 0) THEN
         WRITE(*,*) "MISSING socv from a2o"
         CALL finish(substr)
         ENDIF

    ENDIF

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE a2o_init_coupling

  ! ---------------------------------------------------------------------

  SUBROUTINE a2o_global_start

    USE messy_main_timer,            ONLY: lstart, lresume, delta_time
    USE messy_main_grid_def_mem_bi,  ONLY: ngpblks, nproma, npromz
    USE messy_main_grid_def_bi,      ONLY: gboxarea_2d, philat_2d
    USE messy_main_data_bi,          ONLY: slf,slm           &
                                         , tsw, siced, sni   &
                                         , ocu, ocv, alake
    USE messy_main_timer_bi,         ONLY: event_state
    USE messy_main_timer,            ONLY: current_date, previous_date
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, & 
                                           GP_2D_HORIZONTAL,GP_3D_INT, &
                                           GP_3D_MPIOM, GP_2D_MPIOM
    USE messy_mpiom_mem_e5,          ONLY: dlxp, dlyp &
                                         , bounds_exch, gather_arr, DT &
                                         , allgather_arr ! mz_bk_20120728
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_pe,p_bcast
    USE messy_main_transform_bi,     ONLY: trp_gpdc_gpgl
#if defined(MPIOM_13B)
    USE messy_mpiom_mem_e5,          ONLY: rhosno,rhowat, weto
#elif defined(MPIOM_2000)
    USE messy_mpiom_mem_e5,          ONLY: rhoref_snow,rhoref_water, weto
#endif

    ! LOCAL
    INTEGER                     :: l,i,j,x,y,n, kproma
    INTEGER                     :: maxlev
    REAL(dp) :: sum_a          
    REAL(dp) :: sum_o          

    REAL(dp) :: factor, tsw_old, conti ! mz_ss_20141203: tsw_old, conti 
    
    REAL(dp), DIMENSION(:,:), POINTER :: agf=>NULL()         ! result
    REAL(dp), DIMENSION(:,:), POINTER :: ogf=>NULL()         ! result

    REAL(dp), PARAMETER :: radius_earth = 6371000.0_dp ! radius of the Earth [m]

!-----------------------------------------------------------------
!             CALCULATION OF THE INTERPOLATION
!           BETWEEN ATMOSPHERIC AND OCEANIC GRID
!
!       must be done here: ECHAM5 lsf (land sea fraction)
!                    not read in before 
!-----------------------------------------------------------------

    IF (lstart.or.lresume) THEN
            CALL GRID_INTERPOLATION
    END IF

    IF (r_reset_XA2O.eq.1_dp.or.lstart) THEN
      DO i=1, NA2O
         IF ( XA2O(i)%ok ) XA2O(i)%a_tmp = 0.0_dp
      ENDDO
      r_reset_XA2O = 0.0_dp
      dt_a2o = 0.0_dp
    ENDIF

    IF (r_reset_XO2A.eq.1_dp.or.lstart) THEN
      DO i=1, NO2A
         IF (XO2A(i)%ok) XO2A(i)%o_tmp = 0.0_dp
      ENDDO
      r_reset_XO2A =  0.0_dp 
      dt_o2a = 0.0_dp
    ENDIF

    l_trig_a2o  = event_state(ev_trig_a2o, current_date) 
    l_trig_o2a  = event_state(ev_trig_o2a, previous_date) 

!-----------------------------------------------------------------
!               OCEAN ---> ATMOSPHERE  
!-----------------------------------------------------------------
      IF (LRUNSM_O2A) THEN ! LRUNSM_O2A

        !------------------FIELD COLLECTION ---------------------------!

        dt_o2a = dt_o2a+delta_time

        DO i=1, NO2A
           IF (.NOT. XO2A(i)%ok) CYCLE
            IF (XO2A(i)%rep == GP_3D_MPIOM) THEN
               maxlev=SIZE(XO2A(i)%o_3d,DIM=3)
               SELECT CASE (XO2A(i)%io%cmethod)
                 CASE(A2O_CSUM)
                    DO l=1,maxlev
                       XO2A(i)%o_tmp(1:IE,1:JE) =      &
                          XO2A(i)%o_tmp(1:IE,1:JE) +   &
                          XO2A(i)%o_3d(1:IE,1:JE,l)*delta_time
                    ENDDO
                 CASE(A2O_CAVE)
                    DO l=1,maxlev
                       XO2A(i)%o_tmp(1:IE,1:JE) =      &
                          XO2A(i)%o_tmp(1:IE,1:JE) +   &
                          XO2A(i)%o_3d(1:IE,1:JE,l)*delta_time
                    ENDDO
                     XO2A(i)%o_tmp(1:IE,1:JE) =      &
                     XO2A(i)%o_tmp(1:IE,1:JE)/maxlev
                 CASE(A2O_CMAX)
                     XO2A(i)%o_tmp(1:IE,1:JE) =     &
                        XO2A(i)%o_tmp(1:IE,1:JE) +   &
                        MAXVAL(XO2A(i)%o_3d(1:IE,1:JE,1:maxlev),DIM=3)* &
                                delta_time
                 CASE(A2O_CLLEV)
                     XO2A(i)%o_tmp(1:IE,1:JE) =    &
                        XO2A(i)%o_tmp(1:IE,1:JE) +   &
                        XO2A(i)%o_3d(1:IE,1:JE,maxlev)*delta_time
                 CASE(A2O_CHLEV)
                     XO2A(i)%o_tmp(1:IE,1:JE) =    &
                        XO2A(i)%o_tmp(1:IE,1:JE) +   &
                        XO2A(i)%o_3d(1:IE,1:JE,1)*delta_time
                 END SELECT
            ELSE IF (XO2A(i)%rep == GP_2D_MPIOM) THEN
                XO2A(i)%o_tmp(1:IE,1:JE) =      &
                XO2A(i)%o_tmp(1:IE,1:JE) +      &
                    XO2A(i)%o_2d(1:IE,1:JE)*delta_time
            ENDIF

        END DO

        IF(l_trig_o2a.or.lstart) THEN !l_trig_o2a

        !------------------VECTOR TRANSFORMATION ---------------------------!
    
           DO i=1, NVEC_O2A
              IF (XVEC_O2A(i)%lexist) THEN
                 CALL rotate_2_ne(XO2A(XVEC_O2A(i)%i_x)%o_tmp,XO2A(XVEC_O2A(i)%i_y)%o_tmp,ie,je)
              ENDIF
           ENDDO
    
        !------------------FIELD TRANSFORMATION ---------------------------!
    
           DO i=1, NO2A
                
              IF (.NOT. XO2A(i)%ok) CYCLE

                 ! average the field if necessary
                 IF (XO2A(i)%io%tmethod == A2O_TAVE) THEN
                        XO2A(i)%o_tmp(1:IE,1:JE) =              &
                           XO2A(i)%o_tmp(1:IE,1:JE)/dt_o2a    
                 ENDIF
                ALLOCATE(ogf(IE_G,JE_G))
                ! mz_bk_20120728+
!!$                CALL gather_arr(XO2A(i)%o_tmp,ogf,0)
!!$                CALL p_bcast(ogf,0)
                CALL allgather_arr(XO2A(i)%o_tmp,ogf)
                ! mz_bk_20120728-

                DO y=1,ngpblks
                kproma = nproma
                IF (y==ngpblks) kproma = npromz
                DO x=1,kproma
                 XO2A(i)%a_2d(x,y)=0.0_dp
                 DO n=1,XO2A(i)%o2a_tr(x,y)%nn
                     XO2A(i)%a_2d(x,y)= XO2A(i)%a_2d(x,y)+    & 
                       ogf(XO2A(i)%o2a_tr(x,y)%om_ilon_g(n),  &
                           XO2A(i)%o2a_tr(x,y)%om_ilat_g(n))* &
                           XO2A(i)%o2a_tr(x,y)%weight(n)
                 ENDDO
                 factor=1.0_dp
                 IF (XO2A(i)%interp%map_method.eq.map_type_conserv.and. &
                     XO2A(i)%interp%norm_opt.eq.norm_opt_dstarea) THEN   
                     IF (XO2A(i)%o2a_tr(x,y)%frac.gt.tiny) THEN
                         factor = 1._dp/(XO2A(i)%o2a_tr(x,y)%frac)
                     ELSE
                         factor = 0.0_dp
                     ENDIF
                 ENDIF
                 IF (XO2A(i)%interp%map_method.eq.map_type_conserv.and. &
                     XO2A(i)%interp%norm_opt.eq.norm_opt_none)     THEN
                     IF (XO2A(i)%o2a_tr(x,y)%frac.gt.tiny) THEN
                         factor = 1._dp/(XO2A(i)%o2a_tr(x,y)%frac*           &
                                 (gboxarea_2d(x,y)/(radius_earth**2))) 
                     ELSE
                         factor = 0.0_dp
                     ENDIF
                 ENDIF
                 XO2A(i)%a_2d(x,y)= XO2A(i)%a_2d(x,y)*factor
                ENDDO !1,kproma
                ENDDO !1,ngpblks

                DEALLOCATE(ogf)

               CALL check_max_min(XO2A(i)%a_2d(1:nproma,1:ngpblks),    &
                       XO2A(i)%io%max_value,XO2A(i)%io%min_value)
           
           ENDDO
           
           r_reset_XO2A = 1.0_dp
    
           IF (lcouple) THEN ! here ECHAM5 the MPIOM values.
!mz_ap_20080918 :  Only the open water value are modified.
!              WHERE(slf(1:nproma,1:ngpblks).LT.1.0_dp.and.alake(1:nproma,1:ngpblks).eq.0._dp) &
              DO j=1,ngpblks
                kproma = nproma
                IF (j==ngpblks) kproma = npromz
                DO i=1,kproma
                  !IF ((slf(i,j)+alake(i,j)).LT.1.0_dp-zepsec) THEN 
! mz_ss_20141203+
! Note: there might be lakes next to the coast!
!!$                  IF (slf(i,j).lt.1.0_dp.and.alake(i,j).eq.0.0_dp) THEN 
!!$                    tsw   (i,j) = tsw_a2o   (i,j)
!!$                    ocu   (i,j) = ocu_a2o   (i,j)   
!!$                    ocv   (i,j) = ocv_a2o   (i,j)   
!!$                    seaice(i,j) = seaice_a2o(i,j)
                  IF (slf(i,j).lt.1.0_dp) THEN
                    tsw_old = tsw(i,j)
                    conti = alake(i,j)+slf(i,j)
                    conti = min(conti,1.0_dp)

                    tsw   (i,j) = conti*tsw(i,j)+(1._dp-conti)*tsw_a2o(i,j)
                    ocu   (i,j) = conti*ocu(i,j)+(1._dp-conti)*ocu_a2o(i,j)   
                    ocv   (i,j) = conti*ocv(i,j)+(1._dp-conti)*ocv_a2o(i,j)   
                    seaice(i,j) = conti*seaice(i,j)+(1._dp-conti)*seaice_a2o(i,j) 
!mz_ss_20141203-
                    IF (seaice_a2o(i,j) > 0.0_dp ) THEN
                      !sicomo --> seaice_a2o  : ice compactness (fraction of ice) 
                      !sicsno --> sni_a2o    : snow thikness
                      !sictho --> siced_a2o  : ice thikness 
                      ! IN ECHAM
                      !siced  = ice_depth = ice thickness <-- siced_a2o <-- sictho
                      !sni    = water equivalent of snow on ice <-- sni_a2o <-- sicsno
                      !seaice = fraction of ice <-- seaice <-- sicomo
                      siced(i,j) = MIN(1000.0_dp,(siced_a2o(i,j)/seaice_a2o(i,j)))
                      sni(i,j)   = MIN(1000.0_dp,(sni_a2o(i,j)/seaice_a2o(i,j)))   
#if defined(MPIOM_13B)
                      sni(i,j)   = sni(i,j)*rhosno/rhowat 
#elif defined(MPIOM_2000)
                      sni(i,j)   = sni(i,j)*rhoref_snow/rhoref_water 
#endif
                    ELSE
                      siced(i,j)  = 0.0_dp
                      sni(i,j)    = 0.0_dp
                      seaice(i,j) = 0.0_dp
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
           ENDIF
    
        ENDIF !l_trig_a2o

      ENDIF ! if (LRUNSM_O2A)
                 
  END SUBROUTINE a2o_global_start

  ! ---------------------------------------------------------------------

  SUBROUTINE a2o_global_end

    USE messy_main_grid_def_bi,      ONLY: gboxarea_2d
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, npromz &
                                         , ngpblks, nlev    
    USE messy_main_timer_bi,         ONLY: event_state
    USE messy_main_timer,            ONLY: current_date, time_days &
                                         , delta_time, lstart 
    USE messy_mpiom_e5,              ONLY: DT_steps
    USE messy_main_constants_mem,    ONLY: STRLEN_ULONG
    USE messy_main_constants_mem,    ONLY: rhoh2o=>rho_H2O
    USE messy_main_transform_bi,     ONLY: trp_gpdc_gpgl
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, & 
                                           GP_3D_1LEV,                 &
                                           GP_2D_HORIZONTAL,GP_3D_INT, &
                                           GP_3D_MPIOM, GP_2D_MPIOM
    USE messy_mpiom_mem_e5,          ONLY: dlxp, dlyp, weto,  &
                                           bounds_exch, gather_arr, DT

#if defined(MPIOM_2000)
! mz_bk_20110216+
    USE messy_mpiom_mem_e5,       ONLY: dlxp_g, dlyp_g, dlxu_g, dlyu_g, &
                                        dlxv_g, dlyv_g
! mz_bk_20110216-
#endif
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_pe,p_bcast

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    REAL(dp), PARAMETER :: radius_earth = 6371000.0_dp ! radius of the Earth [m]
    CHARACTER(LEN=*), PARAMETER :: substr = 'a2o_init_coupling'
    INTEGER                     :: status
    INTEGER                     :: l,i,x,y,n,j,kproma
    INTEGER                     :: maxlev
    CHARACTER(LEN=STRLEN_ULONG) :: unit

    REAL(dp) :: sum_a          
    REAL(dp) :: sum_o          

    REAL(dp) :: factor          

    REAL(dp), DIMENSION(:,:), POINTER :: agf=>NULL()         ! result
    REAL(dp), DIMENSION(:,:), POINTER :: ogf=>NULL()         ! result

!-----------------------------------------------------------------
!    GATHER INFORMATION FROM ECHAM5 (previously done in physc.f90)
!-----------------------------------------------------------------

   DO j=1,ngpblks
     kproma = nproma
     IF (j==ngpblks) kproma = npromz
     DO i=1,kproma
        awhea(i,j) = ahflw(i,j)+ahfsw(i,j)+radflxw(i,j) 
        awfre(i,j) = ((rsfl(i,j)+rsfc(i,j)) + &
              (ssfl(i,j)+ssfc(i,j)+evapw(i,j))*         &
              (1._dp-seaice(i,j))) / &
              rhoh2o
        aifre(i,j) = ((ssfl(i,j)+ssfc(i,j) + &
              evapi(i,j)) * (seaice(i,j))) / rhoh2o
     ENDDO
   ENDDO
     
    IF (.NOT.LRUNSM) RETURN

!-----------------------------------------------------------------
!               ATMOSPHERE --->  OCEAN
!-----------------------------------------------------------------
      IF (LRUNSM_A2O) THEN ! LRUNSM_A2O

        !------------------FIELD COLLECTION ---------------------------!

        dt_a2o = dt_a2o+delta_time

        DO i=1, NA2O
           IF (.NOT. XA2O(i)%ok) CYCLE
            IF (XA2O(i)%rep == GP_3D_MID .or. XA2O(i)%rep == GP_3D_INT .or. XA2O(i)%rep == GP_3D_1LEV) THEN
               maxlev=SIZE(XA2O(i)%a_3d,DIM=2)
               SELECT CASE (XA2O(i)%io%cmethod)
                 CASE(A2O_CSUM)
                    DO y=1,ngpblks
                      kproma = nproma
                      IF (y==ngpblks) kproma = npromz
                      DO x=1,kproma
                        DO l=1,maxlev
                           XA2O(i)%a_tmp(x,y) =      &
                              XA2O(i)%a_tmp(x,y) +   &
                              XA2O(i)%a_3d(x,l,y)*delta_time
                        ENDDO
                      ENDDO
                    ENDDO
                 CASE(A2O_CAVE)
                    DO y=1,ngpblks
                      kproma = nproma
                      IF (y==ngpblks) kproma = npromz
                      DO x=1,kproma
                        DO l=1,maxlev
                           XA2O(i)%a_tmp(x,y) =      &
                              XA2O(i)%a_tmp(x,y) +   &
                              XA2O(i)%a_3d(x,l,y)*delta_time
                        ENDDO
                        XA2O(i)%a_tmp(x,y) =      &
                        XA2O(i)%a_tmp(x,y)/maxlev
                      ENDDO  
                    ENDDO  
                 CASE(A2O_CMAX)
                    DO y=1,ngpblks
                      kproma = nproma
                      IF (y==ngpblks) kproma = npromz
                      DO x=1,kproma
                        XA2O(i)%a_tmp(x,y) =     &
                           XA2O(i)%a_tmp(x,y) +   &
                           MAXVAL(XA2O(i)%a_3d(x,1:maxlev,y)) &
                           *delta_time
                      ENDDO
                    ENDDO
                 CASE(A2O_CLLEV)
                    DO y=1,ngpblks
                      kproma = nproma
                      IF (y==ngpblks) kproma = npromz
                      DO x=1,kproma
                        XA2O(i)%a_tmp(x,y) =    &
                           XA2O(i)%a_tmp(x,y) +   &
                           XA2O(i)%a_3d(x,maxlev,y)*delta_time
                      ENDDO
                    ENDDO
                 CASE(A2O_CHLEV)
                    DO y=1,ngpblks
                      kproma = nproma
                      IF (y==ngpblks) kproma = npromz
                      DO x=1,kproma
                        XA2O(i)%a_tmp(x,y) =    &
                           XA2O(i)%a_tmp(x,y) +   &
                           XA2O(i)%a_3d(x,1,y)*delta_time
                      ENDDO
                    ENDDO
                 END SELECT
            ELSE IF (XA2O(i)%rep == GP_2D_HORIZONTAL) THEN
              DO y=1,ngpblks
                kproma = nproma
                IF (y==ngpblks) kproma = npromz
                DO x=1,kproma
                  XA2O(i)%a_tmp(x,y) =      &
                  XA2O(i)%a_tmp(x,y) +      &
                      XA2O(i)%a_2d(x,y)*delta_time
                ENDDO
              ENDDO
            ENDIF
        END DO
    
        IF(l_trig_a2o.or.lstart) THEN !l_trig_a2o
    
        !------------------FIELD TRANSFORMATION ---------------------------!
    
           DO i=1, NA2O
              IF (.NOT. XA2O(i)%ok) CYCLE
                 ! average the field if necessary
                 IF (XA2O(i)%io%tmethod == A2O_TAVE) THEN
                       DO y=1,ngpblks
                         kproma = nproma
                         IF (y==ngpblks) kproma = npromz
                         DO x=1,kproma
                            XA2O(i)%a_tmp(x,y) =              &
                               XA2O(i)%a_tmp(x,y)/dt_a2o    
                         ENDDO
                      ENDDO
                 ENDIF
                CALL trp_gpdc_gpgl(1, XA2O(i)%a_tmp, agf)
                DO x=1,IE
                DO y=1,JE
                 XA2O(i)%o_2d(x,y)=0.0_dp
                 DO n=1,XA2O(i)%a2o_tr(x,y)%nn
                     XA2O(i)%o_2d(x,y)= XA2O(i)%o_2d(x,y)+    & 
                       agf(XA2O(i)%a2o_tr(x,y)%e5_ilon_g(n),  &
                           XA2O(i)%a2o_tr(x,y)%e5_ilat_g(n))* &
                       XA2O(i)%a2o_tr(x,y)%weight(n)
                 ENDDO
                 factor=1.0_dp
                 IF (XA2O(i)%interp%map_method.eq.map_type_conserv.and. &
                     XA2O(i)%interp%norm_opt.eq.norm_opt_dstarea)  THEN      
                     IF (XA2O(i)%a2o_tr(x,y)%frac.gt.tiny) THEN
                         factor = 1/(XA2O(i)%a2o_tr(x,y)%frac)
                     ELSE
                         factor=0.0_dp
                     ENDIF
                 ENDIF
                 IF (XA2O(i)%interp%map_method.eq.map_type_conserv.and. &
                     XA2O(i)%interp%norm_opt.eq.norm_opt_none)     THEN     
                     IF (XA2O(i)%a2o_tr(x,y)%frac.gt.tiny) THEN
                         factor = 1/(XA2O(i)%a2o_tr(x,y)%frac*       &
                                 (dlxp(x,y)*dlyp(x,y)/(radius_earth**2))) 
                     ELSE
                         factor=0.0_dp
                     ENDIF
                 ENDIF
                 XA2O(i)%o_2d(x,y)= XA2O(i)%o_2d(x,y)*factor
                ENDDO
                ENDDO
               !---BOUNDS EXCHANGE! FOR OCEAN! ------------------
#if defined(MPIOM_13B)
               CALL bounds_exch('p',XA2O(i)%o_2d)
#elif defined(MPIOM_2000)
               CALL bounds_exch(1,'p',XA2O(i)%o_2d)
#endif  
               CALL check_max_min(XA2O(i)%o_2d(1:ie,1:je),    &
                       XA2O(i)%io%max_value,XA2O(i)%io%min_value)

               r_reset_XA2O = 1.0_dp

           ENDDO
           
        !------------------VECTOR TRANSFORMATION ---------------------------!
    
           DO i=1, NVEC_A2O
             IF (XVEC_A2O(i)%lexist) THEN
                   ALLOCATE (vec_tmp_x(ie,je))  
                   ALLOCATE (vec_tmp_y(ie,je))  
                      vec_tmp_x = XVEC_A2O(i)%x_2d
                      vec_tmp_y = XVEC_A2O(i)%y_2d  
                      CALL rotate_u(XVEC_A2O(i)%x_2d,vec_tmp_y,ie,je)
                      CALL rotate_v(vec_tmp_x,XVEC_A2O(i)%y_2d ,ie,je)
#if defined(MPIOM_13B)
                   CALL bounds_exch('u',XVEC_A2O(i)%x_2d)
                   CALL bounds_exch('v',XVEC_A2O(i)%y_2d)
#elif defined(MPIOM_2000)
                   CALL bounds_exch(1,'u',XVEC_A2O(i)%x_2d)
                   CALL bounds_exch(1,'v',XVEC_A2O(i)%y_2d)
#endif
                   DEALLOCATE (vec_tmp_x)  
                   DEALLOCATE (vec_tmp_y)  
             ENDIF
           ENDDO
    
        ENDIF !l_trig_a2o.or.lstart
    
      ENDIF ! if (LRUNSM_A2O)


  END SUBROUTINE a2o_global_end

  ! ---------------------------------------------------------------------
  SUBROUTINE a2o_free_memory

    USE messy_main_grid_def_mem_bi,   ONLY : nproma, ngpblks

    IMPLICIT NONE
    INTRINSIC ASSOCIATED

    INTEGER :: n,i,j

    DO n=1,NA2O
      DO i=1,IE
        DO j=1,JE
          IF (ASSOCIATED(XA2O(n)%a2o_tr(i,j)%e5_ilon_g)) THEN
            DEALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilon_g)
            NULLIFY(XA2O(n)%a2o_tr(i,j)%e5_ilon_g)
          ENDIF
          IF (ASSOCIATED(XA2O(n)%a2o_tr(i,j)%e5_ilat_g)) THEN
            DEALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilat_g)
            NULLIFY(XA2O(n)%a2o_tr(i,j)%e5_ilat_g)
          ENDIF
          IF (ASSOCIATED(XA2O(n)%a2o_tr(i,j)%weight)) THEN
            DEALLOCATE(XA2O(n)%a2o_tr(i,j)%weight)
            NULLIFY(XA2O(n)%a2o_tr(i,j)%weight)
          ENDIF
        ENDDO
      ENDDO
      IF (ASSOCIATED(XA2O(n)%a2o_tr)) THEN
        DEALLOCATE(XA2O(n)%a2o_tr)
        NULLIFY(XA2O(n)%a2o_tr)
      ENDIF
    ENDDO

    DO n=1,NO2A
      DO i=1,nproma
        DO j=1,ngpblks
          IF (ASSOCIATED(XO2A(n)%o2a_tr(i,j)%om_ilon_g)) THEN
            DEALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilon_g)
            NULLIFY(XO2A(n)%o2a_tr(i,j)%om_ilon_g)
          ENDIF
          IF (ASSOCIATED(XO2A(n)%o2a_tr(i,j)%om_ilat_g)) THEN
            DEALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilat_g)
            NULLIFY(XO2A(n)%o2a_tr(i,j)%om_ilat_g)
          ENDIF
          IF (ASSOCIATED(XO2A(n)%o2a_tr(i,j)%weight)) THEN
            DEALLOCATE(XO2A(n)%o2a_tr(i,j)%weight)
            NULLIFY(XO2A(n)%o2a_tr(i,j)%weight)
          ENDIF
          IF (ASSOCIATED(XO2A(n)%o2a_tr(i,j)%frac)) THEN
            DEALLOCATE(XO2A(n)%o2a_tr(i,j)%frac)
            NULLIFY(XO2A(n)%o2a_tr(i,j)%frac)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    IF (ASSOCIATED(XO2A(n)%o2a_tr)) THEN
      DEALLOCATE(XO2A(n)%o2a_tr)
      NULLIFY(XO2A(n)%o2a_tr)
    ENDIF


  END SUBROUTINE a2o_free_memory

  ! ---------------------------------------------------------------------

  ! ######################################################################
  ! PRIVATE SUBROUTINES
  ! ######################################################################

SUBROUTINE rotate_2_ne(u_i,v_j,ix,iy)
!----------------------------------------------------------------------
!    BASED ON rotate_2_north_east (mpiom/src)
!
!     rotation of vectors: in ocean models with rotated grids velocity
!     vectors are given in the direction of grid lines and rows. they 
!     have to be rotated in latitudinal and longitudinal direction.
!
!     note: this routine assumes positive meridional flow for a flow
!           from grid point(i,j) to grid point(i,j+1) and positive 
!           zonal flow for a flow from grid point(i,j) to point(i+1,j).
!           this is not the case for mpi-om!
!
!           if this routine is used to rotate data of mpi-om, the 
!           logical change_sign_v needs to be true.
! note here for the coupling fields u-i,v_j are on the non-verlapping
! (ie-2) grid, furthermore, the velocity fields were previously
! interpolated onto the scalar points !
!
!h.haak: 07.10.2005 vectorisation and omp directives      
!----------------------------------------------------------------------
  !
  USE mo_param1
  USE mo_commo1,   ONLY: gila,giph,weto
!  USE mo_units
  !MPIOM-core
  USE messy_mpiom_mem_e5, ONLY : bounds_exch
  !-----------------------------------------------------------------------
  !     local variables
  !-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: ix,iy
  REAL(dp), INTENT(INOUT) ::  u_i(ix,iy),v_j(ix,iy)             ! vector component in i-direction
  REAL(dp) lat(ix,iy),lon(ix,iy)     ! latitudes and longitudes
  REAL(dp) u_lon(ix,iy),v_lat(ix,iy) ! vector component in logitudinal direction 

  REAL(dp) dlat_i, dlat_j,dlon_i,dlon_j,dist_i,dist_j
  REAL(dp) lat_factor,pi
!  real(dp) absold,absnew

  INTEGER i,j,ip1,im1,jp1,jm1

  LOGICAL change_sign_u,change_sign_v
  !-----------------------------------------------------------------------!
  !     specification whether change in sign is needed for the input arrays
  change_sign_u=.FALSE.
  change_sign_v=.TRUE.

  !     transformation to radians
  !     -------------------------
!  IF (ix .NE. ie_g-2) WRITE(*,*) 'alarm boundary in rotate!'
  pi = 3.14159265359

#ifdef _MESSY_OMP
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,ip1,jp1,im1,jm1,dlat_i,dlat_j,dlon_i,dlon_j,dist_i,dist_j)

!$OMP DO
#endif
  DO j=1,iy
     DO i=1,ix
        lat(i,j)=giph(2*i,2*j)
        lon(i,j)=gila(2*i,2*j)
     END DO
  END DO
#ifdef _MESSY_OMP
!$OMP END DO
#endif

  !     initialization
  !     --------------
#ifdef _MESSY_OMP
!$OMP DO
#endif
  DO i = 1, ix
     DO j = 1, iy
        v_lat(i,j) = 0.0_dp 
        u_lon(i,j) = 0.0_dp 
     END DO
  END DO
#ifdef _MESSY_OMP
!$OMP END DO
#endif

  IF (change_sign_u) THEN
#ifdef _MESSY_OMP
!$OMP DO
#endif
     DO i = 1, ix
        DO j = 1, iy
           u_i(i,j)=u_i(i,j)*(-1.0_dp)
        END DO
     END DO
#ifdef _MESSY_OMP
!$OMP END DO
#endif
  ENDIF


  IF (change_sign_v) THEN
#ifdef _MESSY_OMP
!$OMP DO
#endif
     DO i = 1, ix
        DO j = 1, iy
           v_j(i,j)=v_j(i,j)*(-1.0_dp)
        END DO
     END DO
#ifdef _MESSY_OMP
!$OMP END DO
#endif
  ENDIF


  !     rotation
  !     --------
#ifdef _MESSY_OMP
!$OMP DO
#endif
  DO i = 2, ix-1
     DO j = 2, iy-1

        ip1 = i + 1
        im1 = i - 1
        jp1 = j + 1
        jm1 = j - 1

! only for global... here not needed
!        IF (ip1 > ix) ip1 = ip1 - ix ! the 0-meridian
!        IF (im1 < 1 ) im1 = ix
!        IF (jp1 > iy) THEN           ! treatment of the last..
!           jp1 = j
!        ENDIF
!        IF (jm1 < 1 ) THEN ! .. and the fist grid-row
!           jm1 = j
!        ENDIF

        !                 difference in latitudes
        dlat_i = lat(ip1,j) - lat(im1,j)
        dlat_j = lat(i,jp1) - lat(i,jm1)

        !                 difference in longitudes                  
        dlon_i = lon(ip1,j) - lon(im1,j)
        IF (dlon_i >   pi)  dlon_i = dlon_i - (2._dp*pi)
        IF (dlon_i < (-pi)) dlon_i = dlon_i + (2._dp*pi)
        dlon_j = lon(i,jp1) - lon(i,jm1)
        IF (dlon_j >   pi)  dlon_j = dlon_j - (2._dp*pi)
        IF (dlon_j < (-pi)) dlon_j = dlon_j + (2._dp*pi)

        lat_factor = COS(lat(i,j))
        dlon_i = dlon_i * lat_factor
        dlon_j = dlon_j * lat_factor

        !                 projection by scalar product
        !                 ----------------------------
        u_lon(i,j) = u_i(i,j)*dlon_i + v_j(i,j)*dlat_i
        v_lat(i,j) = u_i(i,j)*dlon_j + v_j(i,j)*dlat_j

        dist_i = SQRT(dlon_i**2+dlat_i**2)
        dist_j = SQRT(dlon_j**2+dlat_j**2)
        IF (dist_i /= 0._dp .AND. dist_j /= 0._dp) THEN
           u_lon(i,j) = u_lon(i,j)/dist_i
           v_lat(i,j) = v_lat(i,j)/dist_j
        ELSE
           u_lon(i,j) = 0.0_dp 
           v_lat(i,j) = 0.0_dp  
        ENDIF

        !                  absold = sqrt(u_i(i,j)**2 + v_j(i,j)**2)
        !                  absnew = sqrt(u_lon(i,j)**2 + v_lat(i,j)**2)
        !                  print*, absold, absnew

        !                 test orthogonality
        !                 ------------------
                          if ((dlon_i*dlon_j+dlat_j*dlat_i) > 0.1) then             
                             write(*,*) 'A2O (rotate_2_ne): orthogonal? ', i, j,               &
                                 (dlon_i*dlon_j+dlat_j*dlat_i)
                          endif

     END DO
  END DO
#ifdef _MESSY_OMP
!$OMP END DO
#endif

  ! write back to input field!
#ifdef _MESSY_OMP
!$OMP DO
#endif
  DO i=2,ix-1
     DO j=2,iy-1
        u_i(i,j)=u_lon(i,j)*weto(i,j,1)
        v_j(i,j)=v_lat(i,j)*weto(i,j,1)
     END DO
  END DO
#ifdef _MESSY_OMP
!$OMP END DO

!$OMP END PARALLEL
#endif

#if defined(MPIOM_13B)
  CALL bounds_exch('p',u_i,'rotate_2_ne 1')
  CALL bounds_exch('p',v_j,'rotate_2_ne 2')
#elif defined(MPIOM_2000)
  CALL bounds_exch(1,'p',u_i,'rotate_2_ne 1')
  CALL bounds_exch(1,'p',v_j,'rotate_2_ne 2')
#endif

END SUBROUTINE rotate_2_ne

!-------------------------------------------------------------------------

  SUBROUTINE a2o_read_nml_cpl(status, iou)

    ! a2o MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Pozzer Andrea, MPICH, Dec 2007

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'a2o_read_nml_cpl'

    NAMELIST /CPL/ trig_a2o, trig_o2a,               & ! exchange time scale
                   num_neighbors,                    & ! interpolation general
                   north_thresh, south_thresh,       &
                   max_subseg, max_iter, converge,   &
                   A2O, INTA2O, VEC_A2O, O2A, INTO2A, VEC_O2A

    ! LOCAL
    LOGICAL        :: lex      ! file exists ?
    INTEGER        :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE a2o_read_nml_cpl
!-----------------------------------------------------------------------

  SUBROUTINE restrict_radiants(ARRAYLON,ARRAYLAT)

    ! a2o MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)

    IMPLICIT NONE

    ! I/O
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: ARRAYLON     
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: ARRAYLAT     

!-----------------------------------------------------------------------
!
!     convert longitudes to 0,2pi interval
!
!-----------------------------------------------------------------------

     where (ARRAYLON .gt. pi2)  ARRAYLON = ARRAYLON - pi2
     where (ARRAYLON .lt. zero) ARRAYLON = ARRAYLON + pi2

!-----------------------------------------------------------------------
!
!     make sure input latitude range is within the machine values
!     for +/- pi/2 
!
!-----------------------------------------------------------------------

      where (ARRAYLAT >  pih) ARRAYLAT =  pih
      where (ARRAYLAT < -pih) ARRAYLAT = -pih

  END SUBROUTINE restrict_radiants

!-----------------------------------------------------------------------

 SUBROUTINE GRID_INTERPOLATION

    ! ECHAM5/MESSy
    USE messy_main_blather_bi,   ONLY: info_bi

    USE messy_main_data_bi,         ONLY: alake, slf
    USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks, nlev   &
                                        , npromz
    USE messy_main_grid_def_bi,    ONLY: philon_2d, philat_2d      & 
                                         ! obsolete, to be removed
                                       , philon, philat, gridarea  &
                                       , gboxarea_2d

    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl
    !MPIOM-core
    ! mz_bk_20110220+
#if defined (MPIOM_13B)
! op_pj_20110407+
!    USE messy_mpiom_mem_e5
    USE messy_main_mpi_bi,        ONLY: p_io!, p_bcast ! op_pj_20110407
    USE messy_mpiom_mem_e5,       p_io_wrong => p_io
! op_pj_20110407-
#elif defined (MPIOM_2000)
    USE messy_main_mpi_bi,        ONLY: p_io, p_bcast ! op_bk_20110802
    USE messy_mpiom_mem_e5, ONLY: dlxp_g, dlyp_g, dlxu_g, dlyu_g, dlxv_g      &
                                , dlyv_g, amsuo, amsue, gila_g, giph_g        &
                                , weto_g                                      &
                                , gather_arr, scatter_arr, bounds_exch
#endif
    ! mz_bk_20110220-

    IMPLICIT NONE
    ! --------------------------!
    ! grid calculation array    !
    ! --------------------------!
    REAL(dp), PARAMETER :: radius_earth = 6371000.0_dp ! radius of the Earth [m]
    INTEGER :: links
    INTEGER,PARAMETER :: nc = 4 ! number of corners
    ! ECHAM5
    INTEGER   :: atm_add
    REAL(dp), DIMENSION(:,:), POINTER  :: local  => NULL() 
    REAL(dp), DIMENSION(:,:), POINTER  :: global => NULL() 
    ! global 2D grid
    REAL(dp), TARGET  :: glatmlon(nlon,nlat)     
    REAL(dp), TARGET  :: glatmlat(nlon,nlat)     
    REAL(dp), TARGET  :: glatmclon(nc,nlon,nlat) 
    REAL(dp), TARGET  :: glatmclat(nc,nlon,nlat) 
    REAL(dp), TARGET  :: glatmarea(nlon,nlat)    
    REAL(dp), TARGET  :: glatmmask(nlon,nlat)    
    REAL(dp), TARGET  :: glatmslf(nlon,nlat) 
    ! LOCAL 2D grid
    REAL(dp), TARGET  :: lcatmlon(nproma,ngpblks)     
    REAL(dp), TARGET  :: lcatmlat(nproma,ngpblks)     
    REAL(dp), TARGET  :: lcatmclon(nc,nproma,ngpblks) 
    REAL(dp), TARGET  :: lcatmclat(nc,nproma,ngpblks) 
    REAL(dp), TARGET  :: lcatmarea(nproma,ngpblks)  
    REAL(dp), TARGET  :: lcatmmask(nproma,ngpblks)  
    REAL(dp), TARGET  :: lcatmslf(nproma,ngpblks) 
    ! GLOBAL 1D grid
    REAL(dp),  TARGET:: g1atmlon(nlon*nlat)  
    REAL(dp),  TARGET:: g1atmlat(nlon*nlat)  
    REAL(dp),  TARGET:: g1atmclon(nc,nlon*nlat) 
    REAL(dp),  TARGET:: g1atmclat(nc,nlon*nlat) 
    REAL(dp),  TARGET:: g1atmarea(nlon*nlat)   
    LOGICAL ,  TARGET:: g1atmmask(nlon*nlat)   
    REAL(dp),  TARGET:: g1atmslf(nlon*nlat) 
    ! LOCAL 1D grid
    REAL(dp),  TARGET:: l1atmlon(nproma*ngpblks)     
    REAL(dp),  TARGET:: l1atmlat(nproma*ngpblks)     
    REAL(dp),  TARGET:: l1atmclon(nc,nproma*ngpblks) 
    REAL(dp),  TARGET:: l1atmclat(nc,nproma*ngpblks) 
    REAL(dp),  TARGET:: l1atmarea(nproma*ngpblks)   
    LOGICAL ,  TARGET:: l1atmmask(nproma*ngpblks)   
    REAL(dp),  TARGET:: l1atmslf(nproma*ngpblks) 

    REAL(dp), DIMENSION(:,:), POINTER  :: gl_slf => NULL() ! sea mask related
    REAL(dp), DIMENSION(:,:), POINTER  :: gl_alake => NULL() ! sea mask related

    ! MPIOM
    REAL(dp)  :: lon(2*ie_g,2*je_g)   ! longitudes of doubled array
    REAL(dp)  :: lat(2*ie_g,2*je_g)   ! latitudes of doubled array
    REAL(dp)  :: lonc(2*ie_g,2*je_g)   ! longitudes of doubled array
    REAL(dp)  :: latc(2*ie_g,2*je_g)   ! latitudes of doubled array
    REAL(dp)  :: mask(ie_g,je_g)
    INTEGER   :: oce_add
    INTEGER   :: ip1,im1          ! i+1, i-1
    REAL(dp)  :: pi180
#if defined(MPIOM_13B)
    REAL(dp), DIMENSION(:,:), POINTER  :: dlxp_g => NULL()
    REAL(dp), DIMENSION(:,:), POINTER  :: dlyp_g => NULL()
    REAL(dp), DIMENSION(:,:), POINTER  :: dlxu_g => NULL()
    REAL(dp), DIMENSION(:,:), POINTER  :: dlyu_g => NULL()
    REAL(dp), DIMENSION(:,:), POINTER  :: dlxv_g => NULL() 
    REAL(dp), DIMENSION(:,:), POINTER  :: dlyv_g => NULL()
#endif
    REAL(dp), DIMENSION(:,:), POINTER  :: AMSUE_G_L1(:,:) => NULL()
    REAL(dp), DIMENSION(:,:), POINTER  :: AMSUO_G_L1(:,:) => NULL()
    ! global 2D grid, scalar grid, vector v grid, vector u grid
    REAL(dp)  :: glocelons(ie_g,je_g)    , glocelonu(ie_g,je_g)    ,  glocelonv(ie_g,je_g)      
    REAL(dp)  :: glocelats(ie_g,je_g)    , glocelatu(ie_g,je_g)    ,  glocelatv(ie_g,je_g)      
    REAL(dp)  :: gloceclons(nc,ie_g,je_g), gloceclonu(nc,ie_g,je_g),  gloceclonv(nc,ie_g,je_g)  
    REAL(dp)  :: gloceclats(nc,ie_g,je_g), gloceclatu(nc,ie_g,je_g),  gloceclatv(nc,ie_g,je_g)  
    REAL(dp)  :: gloceareas(ie_g,je_g)   , gloceareau(ie_g,je_g)   ,  gloceareav(ie_g,je_g)     
    REAL(dp)  :: glocemasks(ie_g,je_g)   , glocemasku(ie_g,je_g)   ,  glocemaskv(ie_g,je_g)     
    ! LOCAL 2D grid
    REAL(dp)  :: lcocelons(ie,je)    ,  lcocelonu(ie,je)    , lcocelonv(ie,je)    
    REAL(dp)  :: lcocelats(ie,je)    ,  lcocelatu(ie,je)    , lcocelatv(ie,je)    
    REAL(dp)  :: lcoceclons(nc,ie,je),  lcoceclonu(nc,ie,je), lcoceclonv(nc,ie,je)
    REAL(dp)  :: lcoceclats(nc,ie,je),  lcoceclatu(nc,ie,je), lcoceclatv(nc,ie,je)
    REAL(dp)  :: lcoceareas(ie,je)   ,  lcoceareau(ie,je)   , lcoceareav(ie,je)      
    REAL(dp)  :: lcocemasks(ie,je)   ,  lcocemasku(ie,je)   , lcocemaskv(ie,je)     
    ! GLOBAL 1D grid
    REAL(dp),  TARGET:: g1ocelons((ie_g-2)*je_g)    ,  g1ocelonu((ie_g-2)*je_g)    , g1ocelonv((ie_g-2)*je_g)     
    REAL(dp),  TARGET:: g1ocelats((ie_g-2)*je_g)    ,  g1ocelatu((ie_g-2)*je_g)    , g1ocelatv((ie_g-2)*je_g)      
    REAL(dp),  TARGET:: g1oceclons(nc,(ie_g-2)*je_g),  g1oceclonu(nc,(ie_g-2)*je_g), g1oceclonv(nc,(ie_g-2)*je_g)  
    REAL(dp),  TARGET:: g1oceclats(nc,(ie_g-2)*je_g),  g1oceclatu(nc,(ie_g-2)*je_g), g1oceclatv(nc,(ie_g-2)*je_g) 
    REAL(dp),  TARGET:: g1oceareas((ie_g-2)*je_g)   ,  g1oceareau((ie_g-2)*je_g)   , g1oceareav((ie_g-2)*je_g)    
    LOGICAL ,  TARGET:: g1ocemasks((ie_g-2)*je_g)   ,  g1ocemasku((ie_g-2)*je_g)   , g1ocemaskv((ie_g-2)*je_g)      
    ! LOCAL 1D grid  
    REAL(dp),  TARGET:: l1ocelons(ie*je)    ,  l1ocelonu(ie*je)    , l1ocelonv(ie*je)    
    REAL(dp),  TARGET:: l1ocelats(ie*je)    ,  l1ocelatu(ie*je)    , l1ocelatv(ie*je)    
    REAL(dp),  TARGET:: l1oceclons(nc,ie*je),  l1oceclonu(nc,ie*je), l1oceclonv(nc,ie*je)
    REAL(dp),  TARGET:: l1oceclats(nc,ie*je),  l1oceclatu(nc,ie*je), l1oceclatv(nc,ie*je) 
    REAL(dp),  TARGET:: l1oceareas(ie*je)   ,  l1oceareau(ie*je)   , l1oceareav(ie*je)   
    LOGICAL ,  TARGET:: l1ocemasks(ie*je)   ,  l1ocemasku(ie*je)   , l1ocemaskv(ie*je)   

    INTEGER                     :: i,j,n,ii,l,k,nn
    INTEGER                     :: row, column

    !----------------------------------------------------------
    !                    GRID CALCULATION
    ! here we precalculate the grid transformation based on
    ! grid properties
    !----------------------------------------------------------

!!$    babystep = 0.0001_dp ! op_pj_20130703
    babystep = 0.00001_dp ! op_pj_20140218 according to modification
                          ! by A. Pozzer in 2.42y-ocean 
                          ! (messy/smcl/messy_main_gridtrafo_scrip.f90)

    !----------------------------------------------------------
    !  GRID CALCULATION !
    !----------------------------------------------------------

      CALL trp_gpdc_gpgl(1, slf,   gl_slf)
      CALL trp_gpdc_gpgl(1, alake, gl_alake)

      ALLOCATE(AMSUE_G_L1(ie_g,je_g),AMSUO_G_L1(ie_g,je_g))
#if defined(MPIOM_13B)
      ALLOCATE(dlxp_g(ie_g,je_g), dlyp_g(ie_g,je_g)       &
              ,dlxu_g(ie_g,je_g), dlyu_g(ie_g,je_g)       &
              ,dlxv_g(ie_g,je_g), dlyv_g(ie_g,je_g))
      CALL gather_arr(dlxp,dlxp_g,p_io)
      CALL gather_arr(dlyp,dlyp_g,p_io)
      CALL gather_arr(dlxu,dlxu_g,p_io)
      CALL gather_arr(dlyu,dlyu_g,p_io)
      CALL gather_arr(dlxv,dlxv_g,p_io)
      CALL gather_arr(dlyv,dlyv_g,p_io)
#endif

      CALL p_bcast(dlxp_g,p_io)
      CALL p_bcast(dlyp_g,p_io)
      CALL p_bcast(dlxu_g,p_io)
      CALL p_bcast(dlyu_g,p_io)
      CALL p_bcast(dlxv_g,p_io)
      CALL p_bcast(dlyv_g,p_io)

      CALL gather_arr(AMSUO(:,:,1),AMSUO_G_L1,p_io)
      CALL gather_arr(AMSUE(:,:,1),AMSUE_G_L1,p_io)
      CALL p_bcast(AMSUO_G_L1,p_io)
      CALL p_bcast(AMSUE_G_L1,p_io)

      !-----------------------------------------
      ! ATMOSPHERIC GRID
      !-----------------------------------------

        !----------------------------------
        ! GLOBAL 2D
        !----------------------------------
        DO i = 1, nlon
           glatmlat(i,:) = philat(:)
        ENDDO
        DO j = 1, nlat
           glatmlon(:,j) = philon(:)
        ENDDO
        WHERE (glatmlon(:,:) < 0._dp)
           glatmlon = glatmlon + 360._dp
        END WHERE
        WHERE (glatmlon(:,:) >= 360._dp)
           glatmlon = glatmlon - 360._dp
        END WHERE
   !
   !-- create 3d arrays of grid cell corner longitudes and latitudes
   !   (grid cell corners must be written in counterclockwise sense)
   !      The writing of corners is optional. If they are missing in the grids
   !      file they will be calculated by scrip. In this case it is better to
   !      calculate the corners by the model. (Scrip-corner will not reach the
   !      poles.) 

   ! order grid centers: 
   ! N-->S, E-->W
   ! corners:
   ! 3 2
   ! 4 1
   !
   
        DO i = 1, nlon-1
           glatmclon(1,i,:) = 0.5_dp * (glatmlon(i+1,1) + glatmlon(i,1))
           glatmclon(2,i,:) = 0.5_dp * (glatmlon(i+1,1) + glatmlon(i,1))
        ENDDO
        DO i = 2, nlon
           glatmclon(4,i,:) = 0.5_dp * (glatmlon(i-1,1) + glatmlon(i,1))
           glatmclon(3,i,:) = 0.5_dp * (glatmlon(i-1,1) + glatmlon(i,1))
        ENDDO
        glatmclon(1,nlon,:) = 0.5_dp * (glatmlon(1,1) + glatmlon(nlon,1) + 360._dp)
        glatmclon(2,nlon,:) = 0.5_dp * (glatmlon(1,1) + glatmlon(nlon,1) + 360._dp)
        glatmclon(4,   1,:) = 0.5_dp * (glatmlon(nlon,1) + glatmlon(1,1) - 360._dp)
        glatmclon(3,   1,:) = 0.5_dp * (glatmlon(nlon,1) + glatmlon(1,1) - 360._dp)
        IF (glatmclon(1,nlon,1) >= 360._dp)  &
             glatmclon(1,nlon,:) = glatmclon(1,nlon,:) - 360._dp
        IF (glatmclon(2,nlon,1) >= 360._dp)  &
             glatmclon(2,nlon,:) = glatmclon(2,nlon,:) - 360._dp
        IF (glatmclon(4,nlon,1) < 0._dp)     &
             glatmclon(4,nlon,:) = glatmclon(4,nlon,:) + 360._dp
        IF (glatmclon(3,nlon,1) < 0._dp)     &
             glatmclon(3,nlon,:) = glatmclon(3,nlon,:) + 360._dp
   
        DO j = 1, nlat-1
           glatmclat(1,:,j) = 0.5_dp * (glatmlat(1,j+1) + glatmlat(1,j))
           glatmclat(4,:,j) = 0.5_dp * (glatmlat(1,j+1) + glatmlat(1,j))
        ENDDO
        DO j = 2, nlat
           glatmclat(3,:,j) = 0.5_dp * (glatmlat(1,j-1) + glatmlat(1,j))
           glatmclat(2,:,j) = 0.5_dp * (glatmlat(1,j-1) + glatmlat(1,j))
        ENDDO
        glatmclat(1,:,nlat) = -90._dp
        glatmclat(4,:,nlat) = -90._dp
        glatmclat(3,:,1) = 90._dp
        glatmclat(2,:,1) = 90._dp
   
        glatmslf = 0._dp
        DO j = 1, nlat
           glatmslf(:,j) = gl_slf(:,j) + gl_alake(:,j)
           ! FROM m^2 to Steradian!!!!
           glatmarea(:,j) = gridarea(j) / (radius_earth**2)
        ENDDO

        ! TODO:
        !special tratment for North Pole:
        !IF(p_pe==p_io) THEN
        !  DO i=1,nlat
        !  write(*,*) i,glatmarea(1,i)
        !  ENDDO
        !ENDIF

   
   !
   !-- create 2d integer sea land mask
   !
        glatmmask(:,:) = 0._dp
!        WHERE (glatmslf(:,:) > 0.999999_dp)
        WHERE (glatmslf(:,:).GT.1.0_dp-zepsec) 
           glatmmask(:,:) = 1._dp
        END WHERE

   !vg     WHERE (slf_sn(:,:) > 0.)
!        WHERE (glatmslf(:,:) > 0._dp .AND. glatmmask(:,:) == 0)
!           glatmarea(:,:) = (1-glatmslf(:,:)) * glatmarea(:,:)
!        END WHERE

        !----------------------------------
        ! LOCAL 2D
        !----------------------------------

        local => lcatmlon
        global=> glatmlon
        call trp_gpdc_gpgl(-1,local,global)
        local => lcatmlat
        global=> glatmlat
        call trp_gpdc_gpgl(-1,local,global)
        DO n=1,4
          local  => lcatmclon(n,:,:)
          global => glatmclon(n,:,:)
          call trp_gpdc_gpgl(-1,local,global)
          local => lcatmclat(n,:,:) 
          global=> glatmclat(n,:,:) 
          call trp_gpdc_gpgl(-1,local,global)
        ENDDO
        local => lcatmarea
        global=> glatmarea
        call trp_gpdc_gpgl(-1,local,global)
        local => lcatmmask
        global=> glatmmask
        call trp_gpdc_gpgl(-1,local,global)

        !----------------------------------
        ! GLOBAL 1D
        !----------------------------------

        DO j = 1, nlat
          DO i = 1, nlon
            atm_add = (j-1)*nlon + i
            g1atmlat(atm_add)  = glatmlat(i,j)  * deg2rad
            g1atmlon(atm_add)  = glatmlon(i,j)  * deg2rad
            g1atmarea(atm_add) = glatmarea(i,j) 
            g1atmmask(atm_add) = .TRUE.

            IF (glatmmask(i,j) == 1._dp) g1atmmask(atm_add) = .FALSE. 
            DO n = 1,4
               g1atmclat(n,atm_add) = glatmclat(n,i,j) * deg2rad
               g1atmclon(n,atm_add) = glatmclon(n,i,j) * deg2rad
            ENDDO
          ENDDO
        ENDDO

        call restrict_radiants(g1atmlon(:),   g1atmlat(:))
        DO n=1,4
          call restrict_radiants(g1atmclon(n,:),g1atmclat(n,:))
        ENDDO

        !----------------------------------
        ! LOCAL 1D
        !----------------------------------

        DO j = 1, ngpblks
          DO i = 1, nproma
            atm_add = (j-1)*nproma + i

            l1atmlat(atm_add)  = lcatmlat(i,j)  * deg2rad
            l1atmlon(atm_add)  = lcatmlon(i,j)  * deg2rad
            l1atmarea(atm_add) =  lcatmarea(i,j) 

            l1atmmask(atm_add) = .TRUE.
            IF (lcatmmask(i,j) == 1._dp) l1atmmask(atm_add) = .FALSE. 
            DO n = 1,4
               l1atmclat(n,atm_add) = lcatmclat(n,i,j) * deg2rad
               l1atmclon(n,atm_add) = lcatmclon(n,i,j) * deg2rad
            ENDDO
          ENDDO
        ENDDO

        call restrict_radiants(l1atmlon(:),   l1atmlat(:))
        DO n=1,4
          call restrict_radiants(l1atmclon(n,:),l1atmclat(n,:))
        ENDDO

! skip array rest
        DO j = 1, ngpblks
          DO i = 1, nproma
            atm_add = (j-1)*nproma + i
            IF (j==ngpblks.and.i.gt.npromz) THEN
              l1atmlat(atm_add)  = -90_dp * deg2rad
              l1atmlon(atm_add)  =  0._dp  * deg2rad
              l1atmmask(atm_add) = .FALSE. 
              l1atmarea(atm_add) = 0._dp 
              DO n = 1,4
                 l1atmclat(n,atm_add) = -90_dp * deg2rad
                 l1atmclon(n,atm_add) =  0._dp * deg2rad
              ENDDO
            ENDIF
          ENDDO
        ENDDO

      !-----------------------------------------
      ! OCEAN GRID
      !-----------------------------------------

     !-- create arrays of longitudes and latitudes
          DO j = 1, 2*je_g
             lon(:,j)=gila_g(:,j)
             lat(:,j)=giph_g(:,j)
          ENDDO
     !
     !--  convert from radiant to degree
     !
          pi180=180._dp/pi
          lon(:,:)=lon(:,:) * pi180
          lat(:,:)=lat(:,:) * pi180
     
          WHERE (lon(:,:) < 0._dp)
             lon(:,:) = lon(:,:) + 360._dp
          END WHERE
     !
     !--  extract scalar/vector grid points
     !
     !     2*ij                            
     !      :                              s: scalar
     !      :                              u: vector-u
     !      6   u  s  u  s  u  s           v: vector-v
     !      5   c  v  c  v  c  v           c: grid cell corners of scalars
     !      4   u  s  u  s  u  s  
     !      3   c  v  c  v  c  v  
     !      2   u  s  u  s  u  s           Line 0 and 1 are identical
     !      1   c  v  c  v  c  v  
     !         
     !          1  2  3  4  5  6 ... 2*ie
     !
     !
     !     corner:
     !     2 3      => It is a rotate grid  => 1 4
     !     1 4      => (East->West)         => 2 3
     !

          DO i = 1, ie_g
             DO j = 1, je_g
     !--     scalar
                glocelats(i,j) = lat(i*2,j*2)
                glocelons(i,j) = lon(i*2,j*2)
     !--     vector - u                  
                glocelatu(i,j) = lat(i*2-1,j*2)
                glocelonu(i,j) = lon(i*2-1,j*2)
     !--     vector - v                  
                glocelatv(i,j) = lat(i*2,j*2-1)
                glocelonv(i,j) = lon(i*2,j*2-1)
     !--     corners of scalar grid cells                  
                latc(i,j) = lat(i*2-1,j*2-1)
                lonc(i,j) = lon(i*2-1,j*2-1)
             ENDDO
          ENDDO
     !
     !--  create corner arrays for SCRIP interpolation
     !
          DO i = 1, ie_g
     !
             im1 = i - 1
             ip1 = i + 1
             IF (im1 == 0) im1 = ie_g-2
             IF (ip1 == ie_g+1) ip1 = 3

!-------------------------------------------------
! REVERSED  to  ORIGINAL from mo_coupling.f90
!------------------------------------------------
     !
     !--     scalar
     !
             DO j = 1, je_g-1
                gloceclons(1,i,j) = lonc(i  ,j)
                gloceclons(2,i,j) = lonc(i  ,j+1)
                gloceclons(3,i,j) = lonc(ip1,j+1)
                gloceclons(4,i,j) = lonc(ip1,j)
                gloceclats(1,i,j) = latc(i  ,j)
                gloceclats(2,i,j) = latc(i  ,j+1)
                gloceclats(3,i,j) = latc(ip1,j+1)
                gloceclats(4,i,j) = latc(ip1,j)
             ENDDO
             gloceclons(1,i,je_g) = lonc     (i  ,je_g)
             gloceclons(2,i,je_g) = glocelonu(i  ,je_g)
             gloceclons(3,i,je_g) = glocelonu(ip1,je_g)
             gloceclons(4,i,je_g) = lonc     (ip1,je_g)
             gloceclats(1,i,je_g) = latc     (i  ,je_g)
             gloceclats(2,i,je_g) = glocelatu(i  ,je_g)
             gloceclats(3,i,je_g) = glocelatu(ip1,je_g)
             gloceclats(4,i,je_g) = latc     (ip1,je_g)

     !
     !--     vector - u
     !
             DO j = 1, je_g-1
                gloceclonu(1,i,j) = glocelonv(im1,j  )
                gloceclonu(2,i,j) = glocelonv(im1,j+1)
                gloceclonu(3,i,j) = glocelonv(i  ,j+1)
                gloceclonu(4,i,j) = glocelonv(i  ,j  )
                gloceclatu(1,i,j) = glocelatv(im1,j  )
                gloceclatu(2,i,j) = glocelatv(im1,j+1)
                gloceclatu(3,i,j) = glocelatv(i  ,j+1)
                gloceclatu(4,i,j) = glocelatv(i  ,j  )
             ENDDO
             gloceclonu(1,i,je_g) = glocelonv(im1,je_g)
             gloceclonu(2,i,je_g) = glocelons(im1,je_g)
             gloceclonu(3,i,je_g) = glocelons(i  ,je_g)
             gloceclonu(4,i,je_g) = glocelonv(i  ,je_g)
             gloceclatu(1,i,je_g) = glocelatv(im1,je_g)
             gloceclatu(2,i,je_g) = glocelats(im1,je_g)
             gloceclatu(3,i,je_g) = glocelats(i  ,je_g)
             gloceclatu(4,i,je_g) = glocelatv(i  ,je_g)
     !
     !--     vector - v
     !
             DO j = 2, je_g
                gloceclonv(1,i,j) = glocelonu(i  ,j-1)
                gloceclonv(2,i,j) = glocelonu(i  ,j  )
                gloceclonv(3,i,j) = glocelonu(ip1,j  )
                gloceclonv(4,i,j) = glocelonu(ip1,j-1)
                gloceclatv(1,i,j) = glocelatu(i  ,j-1)
                gloceclatv(2,i,j) = glocelatu(i  ,j  )
                gloceclatv(3,i,j) = glocelatu(ip1,j  )
                gloceclatv(4,i,j) = glocelatu(ip1,j-1)
             ENDDO
             gloceclonv(1,i,1) = lonc     (i  ,1)
             gloceclonv(2,i,1) = glocelonu(i  ,1)
             gloceclonv(3,i,1) = glocelonu(ip1,1)
             gloceclonv(4,i,1) = lonc     (ip1,1)
             gloceclatv(1,i,1) = latc     (i  ,1)
             gloceclatv(2,i,1) = glocelatu(i  ,1)
             gloceclatv(3,i,1) = glocelatu(ip1,1)
             gloceclatv(4,i,1) = latc     (ip1,1)
          ENDDO

          glocemasks(:,:) = 0
          mask(1:ie_g,1:je_g)=weto_g(1:ie_g, 1:je_g, 1)
          WHERE (mask(:,:) == 0)
             glocemasks = 1._dp
          END WHERE
     
          glocemasku(:,:) = 0._dp
          mask(1:ie_g,1:je_g)=AMSUO_G_L1(1:ie_g, 1:je_g)
          WHERE (mask(:,:) == 0)
             glocemasku = 1._dp
          END WHERE
     
          glocemaskv(:,:) = 0._dp
          mask(1:ie_g,1:je_g)=AMSUE_G_L1(1:ie_g, 1:je_g)
          WHERE (mask(:,:) == 0)
             glocemaskv = 1._dp
          END WHERE
       
     !
     !-- create area array
     !
     ! from m^2 to steradian
           gloceareas(:,:)=dlxp_g(1:ie_g,1:je_g)*dlyp_g(1:ie_g,1:je_g)/(radius_earth**2)
           gloceareau(:,:)=dlxu_g(1:ie_g,1:je_g)*dlyu_g(1:ie_g,1:je_g)/(radius_earth**2)
           gloceareav(:,:)=dlxv_g(1:ie_g,1:je_g)*dlyv_g(1:ie_g,1:je_g)/(radius_earth**2)

        !----------------------------------
        ! LOCAL 2D
        !----------------------------------

          call scatter_arr(glocelons,lcocelons, p_io)
          call scatter_arr(glocelats,lcocelats, p_io)
          call scatter_arr(glocelonu,lcocelonu, p_io)
          call scatter_arr(glocelatu,lcocelatu, p_io)
          call scatter_arr(glocelonv,lcocelonv, p_io)
          call scatter_arr(glocelatv,lcocelatv, p_io)
          do n=1,4
            call scatter_arr(gloceclons(n,:,:),lcoceclons(n,:,:),p_io)
            call scatter_arr(gloceclats(n,:,:),lcoceclats(n,:,:),p_io)
            call scatter_arr(gloceclonu(n,:,:),lcoceclonu(n,:,:),p_io)
            call scatter_arr(gloceclatu(n,:,:),lcoceclatu(n,:,:),p_io)
            call scatter_arr(gloceclonv(n,:,:),lcoceclonv(n,:,:),p_io)
            call scatter_arr(gloceclatv(n,:,:),lcoceclatv(n,:,:),p_io)
          enddo
          call scatter_arr(gloceareas,lcoceareas, p_io)
          call scatter_arr(gloceareau,lcoceareau, p_io)
          call scatter_arr(gloceareav,lcoceareav, p_io)
          call scatter_arr(glocemasks,lcocemasks, p_io)
          call scatter_arr(glocemasku,lcocemasku, p_io)
          call scatter_arr(glocemaskv,lcocemaskv, p_io)

#if defined(MPIOM_13B)
          CALL bounds_exch('p',lcocelons)
          CALL bounds_exch('p',lcocelats)
          CALL bounds_exch('u+',lcocelonu)
          CALL bounds_exch('u+',lcocelatu)
          CALL bounds_exch('v+',lcocelonu)
          CALL bounds_exch('v+',lcocelatu)
          do n=1,4
            CALL bounds_exch('p', lcoceclons(n,:,:))
            CALL bounds_exch('p', lcoceclats(n,:,:))
            CALL bounds_exch('u+', lcoceclonu(n,:,:))
            CALL bounds_exch('u+', lcoceclatu(n,:,:))
            CALL bounds_exch('v+', lcoceclonv(n,:,:))
            CALL bounds_exch('v+', lcoceclatv(n,:,:))
          enddo
          call bounds_exch('p', lcoceareas)
          call bounds_exch('u+',lcoceareau)
          call bounds_exch('v+',lcoceareav)
          call bounds_exch('p', lcocemasks)
          call bounds_exch('u+',lcocemasku)
          call bounds_exch('v+',lcocemaskv)
#elif defined(MPIOM_2000)
          CALL bounds_exch(1,'p',lcocelons)
          CALL bounds_exch(1,'p',lcocelats)
          CALL bounds_exch(1,'u+',lcocelonu)
          CALL bounds_exch(1,'u+',lcocelatu)
          CALL bounds_exch(1,'v+',lcocelonu)
          CALL bounds_exch(1,'v+',lcocelatu)
          do n=1,4
            CALL bounds_exch(1,'p', lcoceclons(n,:,:))
            CALL bounds_exch(1,'p', lcoceclats(n,:,:))
            CALL bounds_exch(1,'u+', lcoceclonu(n,:,:))
            CALL bounds_exch(1,'u+', lcoceclatu(n,:,:))
            CALL bounds_exch(1,'v+', lcoceclonv(n,:,:))
            CALL bounds_exch(1,'v+', lcoceclatv(n,:,:))
          enddo
          call bounds_exch(1,'p', lcoceareas)
          call bounds_exch(1,'u+',lcoceareau)
          call bounds_exch(1,'v+',lcoceareav)
          call bounds_exch(1,'p', lcocemasks)
          call bounds_exch(1,'u+',lcocemasku)
          call bounds_exch(1,'v+',lcocemaskv)
#endif

        !----------------------------------
        ! GLOBAL 1D
        !----------------------------------

        DO i = 2, ie_g-1
          DO j = 1, je_g
            oce_add = (j-1)*(ie_g-2) + (i-1)
            g1ocelats(oce_add)  = glocelats(i,j)  * deg2rad
            g1ocelons(oce_add)  = glocelons(i,j)  * deg2rad
            g1ocelatu(oce_add)  = glocelatu(i,j)  * deg2rad
            g1ocelonu(oce_add)  = glocelonu(i,j)  * deg2rad
            g1ocelatv(oce_add)  = glocelatv(i,j)  * deg2rad
            g1ocelonv(oce_add)  = glocelonv(i,j)  * deg2rad
            g1oceareas(oce_add) = gloceareas(i,j)  
            g1oceareau(oce_add) = gloceareau(i,j)  
            g1oceareav(oce_add) = gloceareav(i,j)  

            g1ocemasks(oce_add) = .TRUE.
            g1ocemasku(oce_add) = .TRUE.
            g1ocemaskv(oce_add) = .TRUE.
            IF (glocemasks(i,j) == 1._dp) g1ocemasks(oce_add) = .FALSE. 
            IF (glocemasku(i,j) == 1._dp) g1ocemasku(oce_add) = .FALSE. 
            IF (glocemaskv(i,j) == 1._dp) g1ocemaskv(oce_add) = .FALSE. 
            DO n = 1,4
               g1oceclats(n,oce_add) = gloceclats(n,i,j) * deg2rad
               g1oceclons(n,oce_add) = gloceclons(n,i,j) * deg2rad
               g1oceclatu(n,oce_add) = gloceclatu(n,i,j) * deg2rad
               g1oceclonu(n,oce_add) = gloceclonu(n,i,j) * deg2rad
               g1oceclatv(n,oce_add) = gloceclatv(n,i,j) * deg2rad
               g1oceclonv(n,oce_add) = gloceclonv(n,i,j) * deg2rad
            ENDDO
          ENDDO
        ENDDO

        call restrict_radiants(g1ocelons (:),  g1ocelats (:))
        call restrict_radiants(g1ocelonv (:),  g1ocelatv (:))
        call restrict_radiants(g1ocelonu (:),  g1ocelatu (:))
        DO n = 1,4
          call restrict_radiants(g1oceclonu(n,:),g1oceclatu(n,:))
          call restrict_radiants(g1oceclonv(n,:),g1oceclatv(n,:))
          call restrict_radiants(g1oceclons(n,:),g1oceclats(n,:))
        ENDDO

        !----------------------------------
        ! LOCAL 1D
        !----------------------------------

        DO i = 1, ie
          DO j = 1, je
            oce_add = (j-1)*ie + i
            l1ocelats(oce_add)  = lcocelats(i,j)  * deg2rad
            l1ocelons(oce_add)  = lcocelons(i,j)  * deg2rad
            l1ocelatu(oce_add)  = lcocelatu(i,j)  * deg2rad
            l1ocelonu(oce_add)  = lcocelonu(i,j)  * deg2rad
            l1ocelatv(oce_add)  = lcocelatv(i,j)  * deg2rad
            l1ocelonv(oce_add)  = lcocelonv(i,j)  * deg2rad
            l1oceareas(oce_add) = lcoceareas(i,j) 
            l1oceareau(oce_add) = lcoceareau(i,j) 
            l1oceareav(oce_add) = lcoceareav(i,j) 
            l1ocemasks(oce_add) = .TRUE.
            l1ocemasku(oce_add) = .TRUE.
            l1ocemaskv(oce_add) = .TRUE.
            IF (lcocemasks(i,j) == 1._dp) l1ocemasks(oce_add) = .FALSE. 
            IF (lcocemasku(i,j) == 1._dp) l1ocemasku(oce_add) = .FALSE. 
            IF (lcocemaskv(i,j) == 1._dp) l1ocemaskv(oce_add) = .FALSE. 
            DO n = 1,4
               l1oceclats(n,oce_add) = lcoceclats(n,i,j) * deg2rad
               l1oceclons(n,oce_add) = lcoceclons(n,i,j) * deg2rad
               l1oceclatu(n,oce_add) = lcoceclatu(n,i,j) * deg2rad
               l1oceclonu(n,oce_add) = lcoceclonu(n,i,j) * deg2rad
               l1oceclatv(n,oce_add) = lcoceclatv(n,i,j) * deg2rad
               l1oceclonv(n,oce_add) = lcoceclonv(n,i,j) * deg2rad
            ENDDO
          ENDDO
        ENDDO

        call restrict_radiants(l1ocelons (:),  l1ocelats (:))
        call restrict_radiants(l1ocelonv (:),  l1ocelatv (:))
        call restrict_radiants(l1ocelonu (:),  l1ocelatu (:))
        DO n = 1,4
          call restrict_radiants(l1oceclonu(n,:),l1oceclatu(n,:))
          call restrict_radiants(l1oceclonv(n,:),l1oceclatv(n,:))
          call restrict_radiants(l1oceclons(n,:),l1oceclats(n,:))
        ENDDO

    !----------------------------------------------------------
    !    SCRIP interpolation
    !----------------------------------------------------------
    ! has to be defined:
    ! grid#_size
    ! (must be in radiant!) : 
    ! grid#_center_lon                       in radiants!
    ! grid#_center_lat                             " 
    ! grid#_corner_lat(number_of_corner,:)         "
    ! grid#_corner_lon(number_of_corner,:)         "
    ! grid#_mask
    ! grid#_area_in (we always used luse_grid#_area)
    ! grid#_frac  set to zero!

    !----------------------------------------------------------
    !  TRANSFORMATION CALCULATION: ATMOSPHERE --> OCEAN
    !----------------------------------------------------------

    IF (LRUNSM_A2O) THEN
       !-----------------------------------------------!
       ! SOURCE GRID: SAME FOR ALL TRANSFORMATION (atm)
       !-----------------------------------------------!
       ! source and destination dimensions do not change
       ! source grid informations
       grid1_mask => g1atmmask
       grid1_center_lat => g1atmlat
       grid1_center_lon => g1atmlon
       grid1_area_in => g1atmarea
       grid1_corner_lat => g1atmclat
       grid1_corner_lon => g1atmclon
       !-----------------------------------------------!

       IF (ANY(MASK=INTCALCA2O(:,map_type_conserv,norm_opt_none))) THEN
               CALL info_bi( ' calculate conservative')
               CALL info_bi( ' normalization : none')
               luse_grid_centers = .false.
               map_type = map_type_conserv
               norm_opt = norm_opt_none
               DO ii=1,4
                 IF (INTCALCA2O(ii,map_type_conserv,norm_opt_none)) THEN
                   IF (ii.eq.oces) THEN
                      CALL info_bi(' destination grid OCES')
                      grid2_mask => l1ocemasks
                      grid2_center_lat => l1ocelats
                      grid2_center_lon => l1ocelons
                      grid2_area_in => l1oceareas
                      grid2_corner_lat => l1oceclats
                      grid2_corner_lon => l1oceclons
                   ELSE IF (ii.eq.oceu) THEN
                      CALL info_bi(' destination grid OCEU')
                      grid2_mask => l1ocemasku
                      grid2_center_lat => l1ocelatu
                      grid2_center_lon => l1ocelonu
                      grid2_area_in => l1oceareau
                      grid2_corner_lat => l1oceclatu
                      grid2_corner_lon => l1oceclonu
                   ELSE IF (ii.eq.ocev) THEN
                      CALL info_bi(' destination grid OCEV')
                      grid2_mask => l1ocemaskv
                      grid2_center_lat => l1ocelatv
                      grid2_center_lon => l1ocelonv
                      grid2_area_in => l1oceareav
                      grid2_corner_lat => l1oceclatv
                      grid2_corner_lon => l1oceclonv
                   ENDIF
                     call remap_init(4,nlon,nlat,4,ie,je)
                     call bounds_calc
                     call remap_vars(1)
                     call remap_conserv
                     call remap_vars(2)
                     call check_mask
                     !distribute:
                     DO n=1,NA2O
                       IF (XA2O(n)%interp%map_method.eq.map_type_conserv.and. &
                           XA2O(n)%interp%norm_opt.eq.norm_opt_none) THEN
                           IF (XA2O(n)%interp%dest_grid.eq.ii) THEN
                            ALLOCATE(XA2O(n)%a2o_tr(ie,je))
                            DO i = 1, ie
                               DO j = 1, je
                                oce_add = (j-1)*(ie) + i
                                links=0
                                do l=1,num_links_map1  
                                 IF (grid2_add_map1(l).eq.oce_add) links=links+1
                                end do
                                IF (XA2O(n)%io%lmask.and..not.grid2_mask(oce_add)) links=0
                                XA2O(n)%a2o_tr(i,j)%nn = links
                                IF (links==0) CYCLE
                                ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links))
                                ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links))
                                ALLOCATE(XA2O(n)%a2o_tr(i,j)%weight(links))
                                ALLOCATE(XA2O(n)%a2o_tr(i,j)%frac)
                                links=0
                                do l=1,num_links_map1  
                                 IF (grid2_add_map1(l).eq.oce_add) THEN
                                   links=links+1
                                   row =  MOD(grid1_add_map1(l),nlon)
                                   IF (row == 0 ) THEN
                                     row = nlon
                                   ELSE
                                     row = row
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/nlon)+1
                                   IF (row == nlon) column=column-1
                                   XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links)=row  
                                   XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links)=column
                                   XA2O(n)%a2o_tr(i,j)%weight(links)=wts_map1(1,l)
                                   XA2O(n)%a2o_tr(i,j)%frac=grid2_frac(oce_add)
                                 ENDIF
                                end do
                               ENDDO
                             ENDDO
                           ENDIF
                       ENDIF
                     ENDDO
                     call remap_dealloc
                 ENDIF
               ENDDO
       ENDIF
       IF (ANY(MASK=INTCALCA2O(:,map_type_conserv,norm_opt_dstarea))) THEN
               CALL info_bi( ' calculate conservative')
               CALL info_bi( ' normalization : destarea')
               luse_grid_centers = .false.
               map_type = map_type_conserv
               norm_opt = norm_opt_dstarea
               ! grid init(repeat for each dest grid)
               DO ii=1,4
                 IF (INTCALCA2O(ii,map_type_conserv,norm_opt_dstarea)) THEN
                   IF (ii.eq.oces) THEN
                      CALL info_bi(' destination grid OCES')
                      grid2_mask => l1ocemasks
                      grid2_center_lat => l1ocelats
                      grid2_center_lon => l1ocelons
                      grid2_area_in => l1oceareas
                      grid2_corner_lat => l1oceclats
                      grid2_corner_lon => l1oceclons
                   ELSE IF (ii.eq.oceu) THEN
                      CALL info_bi(' destination grid OCEU')
                      grid2_mask => l1ocemasku
                      grid2_center_lat => l1ocelatu
                      grid2_center_lon => l1ocelonu
                      grid2_area_in => l1oceareau
                      grid2_corner_lat => l1oceclatu
                      grid2_corner_lon => l1oceclonu
                   ELSE IF (ii.eq.ocev) THEN
                      CALL info_bi(' destination grid OCEV')
                      grid2_mask => l1ocemaskv
                      grid2_center_lat => l1ocelatv
                      grid2_center_lon => l1ocelonv
                      grid2_area_in => l1oceareav
                      grid2_corner_lat => l1oceclatv
                      grid2_corner_lon => l1oceclonv
                   ENDIF
                     call remap_init(4,nlon,nlat,4,ie,je)
                     call bounds_calc
                     call remap_vars(1)
                     call remap_conserv
                     call remap_vars(2)
                     call check_mask
                      !distribute:
                      DO n=1,NA2O
                        IF (XA2O(n)%interp%map_method.eq.map_type_conserv.and. &
                            XA2O(n)%interp%norm_opt.eq.norm_opt_dstarea) THEN
                            IF (XA2O(n)%interp%dest_grid.eq.ii) THEN
                             ALLOCATE(XA2O(n)%a2o_tr(ie,je))
                             DO i = 1, ie
                                DO j = 1, je
                                 oce_add = (j-1)*ie + i
                                 links=0
                                 do l=1,num_links_map1  
                                  IF (grid2_add_map1(l).eq.oce_add) links=links+1
                                 end do
                                 IF (XA2O(n)%io%lmask.and..not.grid2_mask(oce_add)) links=0
                                 XA2O(n)%a2o_tr(i,j)%nn = links
                                 IF (links==0) CYCLE
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%weight(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%frac)
                                 links=0
                                 do l=1,num_links_map1  
                                 IF (grid2_add_map1(l).eq.oce_add) THEN
                                   links=links+1
                                   row =  MOD(grid1_add_map1(l),nlon)
                                   IF (row == 0 ) THEN
                                     row = nlon
                                   ELSE
                                     row = row
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/nlon)+1
                                   IF (row == nlon) column=column-1
                                   XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links)=row  
                                   XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links)=column
                                   XA2O(n)%a2o_tr(i,j)%weight(links)=wts_map1(1,l)
                                   XA2O(n)%a2o_tr(i,j)%frac=grid2_frac(oce_add)
                                 ENDIF
                                 end do
                                ENDDO
                              ENDDO
                            ENDIF
                        ENDIF
                      ENDDO
                     call remap_dealloc
                 ENDIF
               ENDDO
       ENDIF
       IF (ANY(MASK=INTCALCA2O(:,map_type_conserv,norm_opt_frcarea))) THEN
               CALL info_bi( ' calculate conservative')
               CALL info_bi( ' normalization : fracarea')
               luse_grid_centers = .false.
               map_type = map_type_conserv
               norm_opt = norm_opt_frcarea
               ! grid init(repeat for each dest grid)
               DO ii=1,4
                 IF (INTCALCA2O(ii,map_type_conserv,norm_opt_frcarea)) THEN
                   IF (ii.eq.oces) THEN
                      CALL info_bi(' destination grid OCES')
                      grid2_mask => l1ocemasks
                      grid2_center_lat => l1ocelats
                      grid2_center_lon => l1ocelons
                      grid2_area_in => l1oceareas
                      grid2_corner_lat => l1oceclats
                      grid2_corner_lon => l1oceclons
                   ELSE IF (ii.eq.oceu) THEN
                      CALL info_bi(' destination grid OCEU')
                      grid2_mask => l1ocemasku
                      grid2_center_lat => l1ocelatu
                      grid2_center_lon => l1ocelonu
                      grid2_area_in => l1oceareau
                      grid2_corner_lat => l1oceclatu
                      grid2_corner_lon => l1oceclonu
                   ELSE IF (ii.eq.ocev) THEN
                      CALL info_bi(' destination grid OCEV')
                      grid2_mask => l1ocemaskv
                      grid2_center_lat => l1ocelatv
                      grid2_center_lon => l1ocelonv
                      grid2_area_in => l1oceareav
                      grid2_corner_lat => l1oceclatv
                      grid2_corner_lon => l1oceclonv
                   ENDIF
                     call remap_init(4,nlon,nlat,4,ie,je)
                     call bounds_calc
                     call remap_vars(1)
                     call remap_conserv
                     call remap_vars(2)
                     call check_mask
                      !distribute:
                      DO n=1,NA2O
                        IF (XA2O(n)%interp%map_method.eq.map_type_conserv.and. &
                            XA2O(n)%interp%norm_opt.eq.norm_opt_frcarea) THEN
                            IF (XA2O(n)%interp%dest_grid.eq.ii) THEN
                             ALLOCATE(XA2O(n)%a2o_tr(ie,je))
                             DO i = 1, ie
                                DO j = 1, je
                                 oce_add = (j-1)*ie + i
                                 links=0
                                 do l=1,num_links_map1  
                                  IF (grid2_add_map1(l).eq.oce_add) links=links+1
                                 end do
                                 IF (XA2O(n)%io%lmask.and..not.grid2_mask(oce_add)) links=0
                                 XA2O(n)%a2o_tr(i,j)%nn = links
                                 IF (links==0) CYCLE
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%weight(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%frac)
                                 links=0
                                 do l=1,num_links_map1  
                                 IF (grid2_add_map1(l).eq.oce_add) THEN
                                   links=links+1
                                   row =  MOD(grid1_add_map1(l),nlon)
                                   IF (row == 0 ) THEN
                                     row = nlon
                                   ELSE
                                     row = row
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/nlon)+1
                                   IF (row == nlon) column=column-1
                                   XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links)=row  
                                   XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links)=column
                                   XA2O(n)%a2o_tr(i,j)%weight(links)=wts_map1(1,l)
                                   XA2O(n)%a2o_tr(i,j)%frac=grid2_frac(oce_add)
                                 ENDIF
                                 end do
                                ENDDO
                              ENDDO
                            ENDIF
                        ENDIF
                      ENDDO
                     call remap_dealloc
                 ENDIF
               ENDDO
       ENDIF
       IF (ANY(MASK=INTCALCA2O(:,map_type_bilinear,:))) THEN
               CALL info_bi ( ' calculate bilinear')
               luse_grid_centers = .true.
               map_type = map_type_bilinear
               norm_opt = norm_opt_none
               ! grid init(repeat for each dest grid)
               DO ii=1,4
                 IF (INTCALCA2O(ii,map_type_bilinear,norm_opt_none)) THEN
                   IF (ii.eq.oces) THEN
                      CALL info_bi(' destination grid OCES')
                      grid2_mask => l1ocemasks
                      grid2_center_lat => l1ocelats
                      grid2_center_lon => l1ocelons
                      grid2_area_in => l1oceareas
                      grid2_corner_lat => l1oceclats
                      grid2_corner_lon => l1oceclons
                   ELSE IF (ii.eq.oceu) THEN
                      CALL info_bi(' destination grid OCEU')
                      grid2_mask => l1ocemasku
                      grid2_center_lat => l1ocelatu
                      grid2_center_lon => l1ocelonu
                      grid2_area_in => l1oceareau
                      grid2_corner_lat => l1oceclatu
                      grid2_corner_lon => l1oceclonu
                   ELSE IF (ii.eq.ocev) THEN
                      CALL info_bi(' destination grid OCEV')
                      grid2_mask => l1ocemaskv
                      grid2_center_lat => l1ocelatv
                      grid2_center_lon => l1ocelonv
                      grid2_area_in => l1oceareav
                      grid2_corner_lat => l1oceclatv
                      grid2_corner_lon => l1oceclonv
                   ENDIF
                      call remap_init(4,nlon,nlat,4,ie,je)
                      call bounds_calc
                      call remap_vars(1)
                      call remap_bilin
                      call remap_vars(2)
                      call check_mask
                      !distribute:
                      DO n=1,NA2O
                        IF (XA2O(n)%interp%map_method.eq.map_type_bilinear.and. &
                            XA2O(n)%interp%norm_opt.eq.norm_opt_none) THEN
                            IF (XA2O(n)%interp%dest_grid.eq.ii) THEN
                             ALLOCATE(XA2O(n)%a2o_tr(ie,je))
                             DO i = 1, ie
                                DO j = 1, je
                                 oce_add = (j-1)*ie + i
                                 links=0
                                 do l=1,num_links_map1  
                                  IF (grid2_add_map1(l).eq.oce_add) links=links+1
                                 end do
                                 IF (XA2O(n)%io%lmask.and..not.grid2_mask(oce_add)) links=0
                                 XA2O(n)%a2o_tr(i,j)%nn = links
                                 IF (links==0) CYCLE
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%weight(links))
                                 links=0
                                 do l=1,num_links_map1  
                                 IF (grid2_add_map1(l).eq.oce_add) THEN
                                   links=links+1
                                   row =  MOD(grid1_add_map1(l),nlon)
                                   IF (row == 0 ) THEN
                                     row = nlon
                                   ELSE
                                     row = row
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/nlon)+1
                                   IF (row == nlon) column=column-1
                                   XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links)=row  
                                   XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links)=column
                                   XA2O(n)%a2o_tr(i,j)%weight(links)=wts_map1(1,l)
                                 ENDIF
                                 end do
                                ENDDO
                              ENDDO
                            ENDIF
                        ENDIF
                      ENDDO
                      call remap_dealloc
                 ENDIF
               ENDDO
       ENDIF
       IF (ANY(INTCALCA2O(:,map_type_bicubic,:))) THEN
               CALL info_bi ( ' calculate bicubic')
               luse_grid_centers = .true.
               map_type = map_type_bicubic
               norm_opt = norm_opt_none
               ! grid init(repeat for each dest grid)
               DO ii=1,4
                 IF (INTCALCA2O(ii,map_type_bicubic,norm_opt_none)) THEN
                   IF (ii.eq.oces) THEN
                      CALL info_bi ( ' destination grid OCES')
                      grid2_mask => l1ocemasks
                      grid2_center_lat => l1ocelats
                      grid2_center_lon => l1ocelons
                      grid2_area_in => l1oceareas
                      grid2_corner_lat => l1oceclats
                      grid2_corner_lon => l1oceclons
                   ELSE IF (ii.eq.oceu) THEN
                      CALL info_bi ( ' destination grid OCEU')
                      grid2_mask => l1ocemasku
                      grid2_center_lat => l1ocelatu
                      grid2_center_lon => l1ocelonu
                      grid2_area_in => l1oceareau
                      grid2_corner_lat => l1oceclatu
                      grid2_corner_lon => l1oceclonu
                   ELSE IF (ii.eq.ocev) THEN
                      CALL info_bi ( ' destination grid OCEV')
                      grid2_mask => l1ocemaskv
                      grid2_center_lat => l1ocelatv
                      grid2_center_lon => l1ocelonv
                      grid2_area_in => l1oceareav
                      grid2_corner_lat => l1oceclatv
                      grid2_corner_lon => l1oceclonv
                   ENDIF
                      call remap_init(4,nlon,nlat,4,ie,je)
                      call bounds_calc
                      call remap_vars(1)
                      call remap_bicub
                      call remap_vars(2)
                      call check_mask
                      !distribute:
                      DO n=1,NA2O
                        IF (XA2O(n)%interp%map_method.eq.map_type_bicubic.and. &
                            XA2O(n)%interp%norm_opt.eq.norm_opt_none) THEN
                            IF (XA2O(n)%interp%dest_grid.eq.ii) THEN
                             ALLOCATE(XA2O(n)%a2o_tr(ie,je))
                             DO i = 1, ie
                                DO j = 1, je
                                 oce_add = (j-1)*ie + i
                                 links=0
                                 do l=1,num_links_map1  
                                  IF (grid2_add_map1(l).eq.oce_add) links=links+1
                                 end do
                                 IF (XA2O(n)%io%lmask.and..not.grid2_mask(oce_add)) links=0
                                 XA2O(n)%a2o_tr(i,j)%nn = links
                                 IF (links==0) CYCLE
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%weight(links))
                                 links=0
                                 do l=1,num_links_map1  
                                 IF (grid2_add_map1(l).eq.oce_add) THEN
                                   links=links+1
                                   row =  MOD(grid1_add_map1(l),nlon)
                                   IF (row == 0 ) THEN
                                     row = nlon
                                   ELSE
                                     row = row
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/nlon)+1
                                   IF (row == nlon) column=column-1
                                   XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links)=row  
                                   XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links)=column
                                   XA2O(n)%a2o_tr(i,j)%weight(links)=wts_map1(1,l)
                                 ENDIF
                                 end do
                                ENDDO
                              ENDDO
                            ENDIF
                        ENDIF
                      ENDDO
                      call remap_dealloc
                 ENDIF
               ENDDO
       ENDIF
       IF (ANY(INTCALCA2O(:,map_type_distwgt,:))) THEN
               CALL info_bi ( ' calculate distwgt')
               luse_grid_centers = .true.
               map_type = map_type_distwgt
               norm_opt = norm_opt_none
               DO ii=1,4
                 IF (INTCALCA2O(ii,map_type_distwgt,norm_opt_none)) THEN
                   IF (ii.eq.oces) THEN
                      CALL info_bi ( ' destination grid OCES')
                      grid2_mask => l1ocemasks
                      grid2_center_lat => l1ocelats
                      grid2_center_lon => l1ocelons
                      grid2_area_in => l1oceareas
                      grid2_corner_lat => l1oceclats
                      grid2_corner_lon => l1oceclons
                   ELSE IF (ii.eq.oceu) THEN
                      CALL info_bi ( ' destination grid OCEU')
                      grid2_mask => l1ocemasku
                      grid2_center_lat => l1ocelatu
                      grid2_center_lon => l1ocelonu
                      grid2_area_in => l1oceareau
                      grid2_corner_lat => l1oceclatu
                      grid2_corner_lon => l1oceclonu
                   ELSE IF (ii.eq.ocev) THEN
                      CALL info_bi ( ' destination grid OCEV')
                      grid2_mask => l1ocemaskv
                      grid2_center_lat => l1ocelatv
                      grid2_center_lon => l1ocelonv
                      grid2_area_in => l1oceareav
                      grid2_corner_lat => l1oceclatv
                      grid2_corner_lon => l1oceclonv
                   ENDIF
                      call remap_init(4,nlon,nlat,4,ie,je)
                      call bounds_calc
                      call remap_vars(1)
                      call remap_distwgt
                      call remap_vars(2)
                      ! done on-line....
                      !call check_mask
                      !distribute:
                      DO n=1,NA2O
                        IF (XA2O(n)%interp%map_method.eq.map_type_distwgt.and. &
                            XA2O(n)%interp%norm_opt.eq.norm_opt_none) THEN
                            IF (XA2O(n)%interp%dest_grid.eq.ii) THEN
                             ALLOCATE(XA2O(n)%a2o_tr(ie,je))
                             DO i = 1, ie
                                DO j = 1, je
                                 oce_add = (j-1)*ie + i
                                 links=0
                                 do l=1,num_links_map1  
                                  IF (grid2_add_map1(l).eq.oce_add) links=links+1
                                 end do
                                 IF (XA2O(n)%io%lmask.and..not.grid2_mask(oce_add)) links=0
                                 XA2O(n)%a2o_tr(i,j)%nn = links
                                 IF (links==0) CYCLE
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links))
                                 ALLOCATE(XA2O(n)%a2o_tr(i,j)%weight(links))
                                 links=0
                                 do l=1,num_links_map1  
                                 IF (grid2_add_map1(l).eq.oce_add) THEN
                                   links=links+1
                                   row =  MOD(grid1_add_map1(l),nlon)
                                   IF (row == 0 ) THEN
                                     row = nlon
                                   ELSE
                                     row = row
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/nlon)+1
                                   IF (row == nlon) column=column-1
                                   XA2O(n)%a2o_tr(i,j)%e5_ilon_g(links)=row  
                                   XA2O(n)%a2o_tr(i,j)%e5_ilat_g(links)=column
                                   XA2O(n)%a2o_tr(i,j)%weight(links)=wts_map1(1,l)
                                 ENDIF
                                 end do
                                ENDDO
                              ENDDO
                            ENDIF
                        ENDIF
                      ENDDO
                      call remap_dealloc
                 ENDIF
               ENDDO
       ENDIF
      
      CALL info_bi ('transformation ATMOSPHERE --> OCEAN LOCKED !' )

    ENDIF

    !----------------------------------------------------------
    ! GRID CALCULATION: OCEAN --> ATMOSPHERE
    !----------------------------------------------------------

    IF (LRUNSM_O2A) THEN !lrunsum
       !-----------------------------------------------!
       ! DESTINATION GRID: SAME FOR ALL TRANSFORMATION (atm)
       !-----------------------------------------------!
       grid2_mask => l1atmmask
       grid2_center_lat => l1atmlat
       grid2_center_lon => l1atmlon
       grid2_area_in => l1atmarea
       grid2_corner_lat => l1atmclat
       grid2_corner_lon => l1atmclon
       !-----------------------------------------------!

      IF (ANY(MASK=INTCALCO2A(:,map_type_conserv,norm_opt_none))) THEN
              CALL info_bi ( ' calculate conservative')
              CALL info_bi ( ' normalization : none')
              luse_grid_centers = .false.
              map_type = map_type_conserv
              norm_opt = norm_opt_none
              DO ii=1,4
                IF (INTCALCO2A(ii,map_type_conserv,norm_opt_none)) THEN
                  IF (ii.eq.oces) THEN
                     CALL info_bi ( ' source grid OCES')
                     grid1_corners = 4
                     grid1_mask => g1ocemasks
                     grid1_center_lat => g1ocelats
                     grid1_center_lon => g1ocelons
                     grid1_area_in => g1oceareas
                     grid1_corner_lat => g1oceclats
                     grid1_corner_lon => g1oceclons
                  ELSE IF (ii.eq.oceu) THEN
                     CALL info_bi ( ' source grid OCEU')
                     grid1_corners = 4
                     grid1_mask => g1ocemasku
                     grid1_center_lat => g1ocelatu
                     grid1_center_lon => g1ocelonu
                     grid1_area_in => g1oceareau
                     grid1_corner_lat => g1oceclatu
                     grid1_corner_lon => g1oceclonu
                  ELSE IF (ii.eq.ocev) THEN
                     CALL info_bi ( ' source grid OCEV')
                     grid1_corners = 4
                     grid1_mask => g1ocemaskv
                     grid1_center_lat => g1ocelatv
                     grid1_center_lon => g1ocelonv
                     grid1_area_in => g1oceareav
                     grid1_corner_lat => g1oceclatv
                     grid1_corner_lon => g1oceclonv
                  ENDIF
                    call remap_init(4,ie_g-2,je_g,4,nproma,ngpblks)
                    call bounds_calc
                    call remap_vars(1)
                    call remap_conserv
                    call remap_vars(2)
                    call check_mask
                    !distribute:
                    DO n=1,NO2A
                      IF (XO2A(n)%interp%map_method.eq.map_type_conserv.and. &
                          XO2A(n)%interp%norm_opt.eq.norm_opt_none) THEN
                          IF (XO2A(n)%interp%source_grid.eq.ii) THEN
                           ALLOCATE(XO2A(n)%o2a_tr(nproma,ngpblks))
                           DO j = 1, ngpblks
                             DO i = 1, nproma
                               atm_add = (j-1)*nproma + i
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) links=links+1
                               end do
                               IF (XO2A(n)%io%lmask.and..not.grid2_mask(atm_add)) links=0
                               XO2A(n)%o2a_tr(i,j)%nn = links
                               IF (links==0) CYCLE
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilon_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilat_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%weight(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%frac)
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) THEN
                                   links=links+1
                                   row =MOD(grid1_add_map1(l),(ie_g-2))
                                   IF (row == 0) THEN
                                     row = (ie_g-2)+1 !!!!!
                                   ELSE
                                     row = row+1
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/(ie_g-2))+1
                                   IF (row == (ie_g-2)+1) column=column-1
                                   XO2A(n)%o2a_tr(i,j)%om_ilon_g(links)=row
                                   XO2A(n)%o2a_tr(i,j)%om_ilat_g(links)=column
                                   XO2A(n)%o2a_tr(i,j)%weight(links)=wts_map1(1,l)
                                   XO2A(n)%o2a_tr(i,j)%frac=grid2_frac(atm_add)
                                   IF (j==ngpblks.and.i.gt.npromz) XO2A(n)%o2a_tr(i,j)%weight(links)=0._dp 
                                ENDIF
                               end do
                              ENDDO
                            ENDDO
                          ENDIF
                      ENDIF
                    ENDDO
                    call remap_dealloc
                ENDIF
              ENDDO
      ENDIF
      IF (ANY(MASK=INTCALCO2A(:,map_type_conserv,norm_opt_dstarea))) THEN
              CALL info_bi ( ' calculate conservative')
              CALL info_bi ( ' normalization : destarea')
              luse_grid_centers = .false.
              map_type = map_type_conserv
              norm_opt = norm_opt_dstarea
              DO ii=1,4
                IF (INTCALCO2A(ii,map_type_conserv,norm_opt_dstarea)) THEN
                  IF (ii.eq.oces) THEN
                     CALL info_bi ( ' source grid OCES')
                     grid1_corners = 4
                     grid1_mask => g1ocemasks
                     grid1_center_lat => g1ocelats
                     grid1_center_lon => g1ocelons
                     grid1_area_in => g1oceareas
                     grid1_corner_lat => g1oceclats
                     grid1_corner_lon => g1oceclons
                  ELSE IF (ii.eq.oceu) THEN
                     CALL info_bi ( ' source grid OCEU')
                     grid1_corners = 4
                     grid1_mask => g1ocemasku
                     grid1_center_lat => g1ocelatu
                     grid1_center_lon => g1ocelonu
                     grid1_area_in => g1oceareau
                     grid1_corner_lat => g1oceclatu
                     grid1_corner_lon => g1oceclonu
                  ELSE IF (ii.eq.ocev) THEN
                     CALL info_bi ( ' source grid OCEV')
                     grid1_corners = 4
                     grid1_mask => g1ocemaskv
                     grid1_center_lat => g1ocelatv
                     grid1_center_lon => g1ocelonv
                     grid1_area_in => g1oceareav
                     grid1_corner_lat => g1oceclatv
                     grid1_corner_lon => g1oceclonv
                  ENDIF
                    call remap_init(4,ie_g-2,je_g,4,nproma,ngpblks)
                    call bounds_calc
                    call remap_vars(1)
                    call remap_conserv
                    call remap_vars(2)
                    call check_mask
                    !distribute:
                    DO n=1,NO2A
                      IF (XO2A(n)%interp%map_method.eq.map_type_conserv.and. &
                          XO2A(n)%interp%norm_opt.eq.norm_opt_dstarea) THEN
                          IF (XO2A(n)%interp%source_grid.eq.ii) THEN
                           ALLOCATE(XO2A(n)%o2a_tr(nproma,ngpblks))
                           DO j = 1, ngpblks
                             DO i = 1, nproma
                               atm_add = (j-1)*nproma + i
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) links=links+1
                               end do
                               IF (XO2A(n)%io%lmask.and..not.grid2_mask(atm_add)) links=0
                               XO2A(n)%o2a_tr(i,j)%nn = links
                               IF (links==0) CYCLE
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilon_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilat_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%weight(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%frac)
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) THEN
                                   links=links+1
                                   row =MOD(grid1_add_map1(l),(ie_g-2))
                                   IF (row == 0) THEN
                                     row = (ie_g-2)+1
                                   ELSE
                                     row = row+1
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/(ie_g-2))+1
                                   IF (row == (ie_g-2)+1) column=column-1
                                   XO2A(n)%o2a_tr(i,j)%om_ilon_g(links)=row
                                   XO2A(n)%o2a_tr(i,j)%om_ilat_g(links)=column
                                   XO2A(n)%o2a_tr(i,j)%weight(links)=wts_map1(1,l)
                                   XO2A(n)%o2a_tr(i,j)%frac=grid2_frac(atm_add)
                                   IF (j==ngpblks.and.i.gt.npromz) XO2A(n)%o2a_tr(i,j)%weight(links)=0._dp 
                                ENDIF
                               end do
                              ENDDO
                            ENDDO
                          ENDIF
                      ENDIF
                    ENDDO
                    call remap_dealloc
                ENDIF
              ENDDO
      ENDIF
      IF (ANY(MASK=INTCALCO2A(:,map_type_conserv,norm_opt_frcarea))) THEN
              CALL info_bi ( ' calculate conservative')
              CALL info_bi ( ' normalization : fracarea')
              luse_grid_centers = .false.
              map_type = map_type_conserv
              norm_opt = norm_opt_frcarea
              DO ii=1,4
                IF (INTCALCO2A(ii,map_type_conserv,norm_opt_frcarea)) THEN
                  IF (ii.eq.oces) THEN
                     CALL info_bi ( ' source grid OCES')
                     grid1_corners = 4
                     grid1_mask => g1ocemasks
                     grid1_center_lat => g1ocelats
                     grid1_center_lon => g1ocelons
                     grid1_area_in => g1oceareas
                     grid1_corner_lat => g1oceclats
                     grid1_corner_lon => g1oceclons
                  ELSE IF (ii.eq.oceu) THEN
                     CALL info_bi ( ' source grid OCEU')
                     grid1_corners = 4
                     grid1_mask => g1ocemasku
                     grid1_center_lat => g1ocelatu
                     grid1_center_lon => g1ocelonu
                     grid1_area_in => g1oceareau
                     grid1_corner_lat => g1oceclatu
                     grid1_corner_lon => g1oceclonu
                  ELSE IF (ii.eq.ocev) THEN
                     CALL info_bi ( ' source grid OCEV')
                     grid1_corners = 4
                     grid1_mask => g1ocemaskv
                     grid1_center_lat => g1ocelatv
                     grid1_center_lon => g1ocelonv
                     grid1_area_in => g1oceareav
                     grid1_corner_lat => g1oceclatv
                     grid1_corner_lon => g1oceclonv
                  ENDIF
                    call remap_init(4,ie_g-2,je_g,4,nproma,ngpblks)
                    call bounds_calc
                    call remap_vars(1)
                    call remap_conserv
                    call remap_vars(2)
                    call check_mask
                    !distribute:
                    DO n=1,NO2A
                      IF (XO2A(n)%interp%map_method.eq.map_type_conserv.and. &
                          XO2A(n)%interp%norm_opt.eq.norm_opt_frcarea) THEN
                          IF (XO2A(n)%interp%source_grid.eq.ii) THEN
                           ALLOCATE(XO2A(n)%o2a_tr(nproma,ngpblks))
                           DO j = 1, ngpblks
                             DO i = 1, nproma
                               atm_add = (j-1)*nproma + i
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) links=links+1
                               end do
                               IF (XO2A(n)%io%lmask.and..not.grid2_mask(atm_add)) links=0
                               XO2A(n)%o2a_tr(i,j)%nn = links
                               IF (links==0) CYCLE
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilon_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilat_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%weight(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%frac)
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) THEN
                                   links=links+1
                                   row =MOD(grid1_add_map1(l),(ie_g-2))
                                   IF (row == 0) THEN
                                     row = (ie_g-2)+1
                                   ELSE
                                     row = row+1
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/(ie_g-2))+1
                                   IF (row == (ie_g-2)+1) column=column-1
                                   XO2A(n)%o2a_tr(i,j)%om_ilon_g(links)=row
                                   XO2A(n)%o2a_tr(i,j)%om_ilat_g(links)=column
                                   XO2A(n)%o2a_tr(i,j)%weight(links)=wts_map1(1,l)
                                   XO2A(n)%o2a_tr(i,j)%frac=grid2_frac(atm_add)
                                   IF (j==ngpblks.and.i.gt.npromz) XO2A(n)%o2a_tr(i,j)%weight(links)=0._dp 
                                ENDIF
                               end do
                              ENDDO
                            ENDDO
                          ENDIF
                      ENDIF
                    ENDDO
                    call remap_dealloc
                ENDIF
              ENDDO
      ENDIF
      IF (ANY(MASK=INTCALCO2A(:,map_type_bilinear,:))) THEN
              CALL info_bi ( ' calculate bilinear')
              luse_grid_centers = .true.
              map_type = map_type_bilinear
              norm_opt = norm_opt_none
              DO ii=1,4
                IF (INTCALCO2A(ii,map_type_bilinear,norm_opt_none)) THEN
                  IF (ii.eq.oces) THEN
                     CALL info_bi ( ' source grid OCES')
                     grid1_corners = 4
                     grid1_mask => g1ocemasks
                     grid1_center_lat => g1ocelats
                     grid1_center_lon => g1ocelons
                     grid1_area_in => g1oceareas
                     grid1_corner_lat => g1oceclats
                     grid1_corner_lon => g1oceclons
                  ELSE IF (ii.eq.oceu) THEN
                     CALL info_bi ( ' source grid OCEU')
                     grid1_corners = 4
                     grid1_mask => g1ocemasku
                     grid1_center_lat => g1ocelatu
                     grid1_center_lon => g1ocelonu
                     grid1_area_in => g1oceareau
                     grid1_corner_lat => g1oceclatu
                     grid1_corner_lon => g1oceclonu
                  ELSE IF (ii.eq.ocev) THEN
                     CALL info_bi ( ' source grid OCEV')
                     grid1_corners = 4
                     grid1_mask => g1ocemaskv
                     grid1_center_lat => g1ocelatv
                     grid1_center_lon => g1ocelonv
                     grid1_area_in => g1oceareav
                     grid1_corner_lat => g1oceclatv
                     grid1_corner_lon => g1oceclonv
                  ENDIF
                    call remap_init(4,ie_g-2,je_g,4,nproma,ngpblks)
                    call bounds_calc
                    call remap_vars(1)
                    call remap_bilin
                    call remap_vars(2)
                    call check_mask
                    !distribute:
                    DO n=1,NO2A
                      IF (XO2A(n)%interp%map_method.eq.map_type_bilinear.and. &
                          XO2A(n)%interp%norm_opt.eq.norm_opt_none) THEN
                          IF (XO2A(n)%interp%source_grid.eq.ii) THEN
                           ALLOCATE(XO2A(n)%o2a_tr(nproma,ngpblks))
                           DO j = 1, ngpblks
                             DO i = 1, nproma
                               atm_add = (j-1)*nproma + i
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) links=links+1
                               end do
                               IF (XO2A(n)%io%lmask.and..not.grid2_mask(atm_add)) links=0
                               XO2A(n)%o2a_tr(i,j)%nn = links
                               IF (links==0) CYCLE
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilon_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilat_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%weight(links))
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) THEN
                                   links=links+1
                                   row =MOD(grid1_add_map1(l),(ie_g-2))
                                   IF (row == 0) THEN
                                     row = (ie_g-2)+1
                                   ELSE
                                     row = row+1
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/(ie_g-2))+1
                                   IF (row == (ie_g-2)+1) column=column-1
                                   XO2A(n)%o2a_tr(i,j)%om_ilon_g(links)=row
                                   XO2A(n)%o2a_tr(i,j)%om_ilat_g(links)=column
                                   XO2A(n)%o2a_tr(i,j)%weight(links)=wts_map1(1,l)
                                   IF (j==ngpblks.and.i.gt.npromz) XO2A(n)%o2a_tr(i,j)%weight(links)=0._dp 
                                ENDIF
                               end do
                              ENDDO
                            ENDDO
                          ENDIF
                      ENDIF
                    ENDDO
                    call remap_dealloc
                ENDIF
              ENDDO
      ENDIF
      IF (ANY(INTCALCO2A(:,map_type_bicubic,:))) THEN
              CALL info_bi ( ' calculate bicubic')
              luse_grid_centers = .true.
              map_type = map_type_bicubic
              norm_opt = norm_opt_none
              ! grid init(repeat for each dest grid)
              DO ii=1,4
                IF (INTCALCO2A(ii,map_type_bicubic,norm_opt_none)) THEN
                  IF (ii.eq.oces) THEN
                     CALL info_bi ( ' source grid OCES')
                     grid1_corners = 4
                     grid1_mask => g1ocemasks
                     grid1_center_lat => g1ocelats
                     grid1_center_lon => g1ocelons
                     grid1_area_in => g1oceareas
                     grid1_corner_lat => g1oceclats
                     grid1_corner_lon => g1oceclons
                  ELSE IF (ii.eq.oceu) THEN
                     CALL info_bi ( ' source grid OCEU')
                     grid1_corners = 4
                     grid1_mask => g1ocemasku
                     grid1_center_lat => g1ocelatu
                     grid1_center_lon => g1ocelonu
                     grid1_area_in => g1oceareau
                     grid1_corner_lat => g1oceclatu
                     grid1_corner_lon => g1oceclonu
                  ELSE IF (ii.eq.ocev) THEN
                     CALL info_bi ( ' source grid OCEV')
                     grid1_corners = 4
                     grid1_mask => g1ocemaskv
                     grid1_center_lat => g1ocelatv
                     grid1_center_lon => g1ocelonv
                     grid1_area_in => g1oceareav
                     grid1_corner_lat => g1oceclatv
                     grid1_corner_lon => g1oceclonv
                  ENDIF
                    call remap_init(4,ie_g-2,je_g,4,nproma,ngpblks)
                    call bounds_calc
                    call remap_vars(1)
                    call remap_bicub
                    call remap_vars(2)
                    call check_mask
                    !distribute:
                    DO n=1,NO2A
                      IF (XO2A(n)%interp%map_method.eq.map_type_bicubic.and. &
                          XO2A(n)%interp%norm_opt.eq.norm_opt_none) THEN
                          IF (XO2A(n)%interp%source_grid.eq.ii) THEN
                           ALLOCATE(XO2A(n)%o2a_tr(nproma,ngpblks))
                           DO j = 1, ngpblks
                             DO i = 1, nproma
                               atm_add = (j-1)*nproma + i
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) links=links+1
                               end do
                               IF (XO2A(n)%io%lmask.and..not.grid2_mask(atm_add)) links=0
                               XO2A(n)%o2a_tr(i,j)%nn = links
                               IF (links==0) CYCLE
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilon_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilat_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%weight(links))
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) THEN
                                   links=links+1
                                   row =MOD(grid1_add_map1(l),(ie_g-2))
                                   IF (row == 0) THEN
                                     row = (ie_g-2)+1
                                   ELSE
                                     row = row+1
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/(ie_g-2))+1
                                   IF (row == (ie_g-2)+1) column=column-1
                                   XO2A(n)%o2a_tr(i,j)%om_ilon_g(links)=row
                                   XO2A(n)%o2a_tr(i,j)%om_ilat_g(links)=column
                                   XO2A(n)%o2a_tr(i,j)%weight(links)=wts_map1(1,l)
                                   IF (j==ngpblks.and.i.gt.npromz) XO2A(n)%o2a_tr(i,j)%weight(links)=0._dp 
                                ENDIF
                               end do
                              ENDDO
                            ENDDO
                          ENDIF
                      ENDIF
                    ENDDO
                    call remap_dealloc
                ENDIF
              ENDDO
              ! and distribute to all the fields
      ENDIF
      IF (ANY(INTCALCO2A(:,map_type_distwgt,:))) THEN
              CALL info_bi ( ' calculate distwgt')
              luse_grid_centers = .true.
              map_type = map_type_distwgt
              norm_opt = norm_opt_none
              DO ii=1,4
                IF (INTCALCO2A(ii,map_type_distwgt,norm_opt_none)) THEN
                  IF (ii.eq.oces) THEN
                     CALL info_bi ( ' source grid OCES')
                     grid1_corners = 4
                     grid1_mask => g1ocemasks
                     grid1_center_lat => g1ocelats
                     grid1_center_lon => g1ocelons
                     grid1_area_in => g1oceareas
                     grid1_corner_lat => g1oceclats
                     grid1_corner_lon => g1oceclons
                  ELSE IF (ii.eq.oceu) THEN
                     CALL info_bi ( ' source grid OCEU')
                     grid1_corners = 4
                     grid1_mask => g1ocemasku
                     grid1_center_lat => g1ocelatu
                     grid1_center_lon => g1ocelonu
                     grid1_area_in => g1oceareau
                     grid1_corner_lat => g1oceclatu
                     grid1_corner_lon => g1oceclonu
                  ELSE IF (ii.eq.ocev) THEN
                     CALL info_bi ( ' source grid OCEV')
                     grid1_corners = 4
                     grid1_mask => g1ocemaskv
                     grid1_center_lat => g1ocelatv
                     grid1_center_lon => g1ocelonv
                     grid1_area_in => g1oceareav
                     grid1_corner_lat => g1oceclatv
                     grid1_corner_lon => g1oceclonv
                  ENDIF
                    call remap_init(4,ie_g-2,je_g,4,nproma,ngpblks)
                    call bounds_calc
                    call remap_vars(1)
                    call remap_distwgt
                    call remap_vars(2)
                    ! not needed... done on-line
                    !  call check_mask
                    !distribute:
                    DO n=1,NO2A
                      IF (XO2A(n)%interp%map_method.eq.map_type_distwgt.and. &
                          XO2A(n)%interp%norm_opt.eq.norm_opt_none) THEN
                          IF (XO2A(n)%interp%source_grid.eq.ii) THEN
                           ALLOCATE(XO2A(n)%o2a_tr(nproma,ngpblks))
                           DO j = 1, ngpblks
                             DO i = 1, nproma
                               atm_add = (j-1)*nproma + i
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) links=links+1
                               end do
                               IF (XO2A(n)%io%lmask.and..not.grid2_mask(atm_add)) links=0
                               XO2A(n)%o2a_tr(i,j)%nn = links
                               IF (links==0) CYCLE
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilon_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%om_ilat_g(links))
                               ALLOCATE(XO2A(n)%o2a_tr(i,j)%weight(links))
                               links=0
                               do l=1,num_links_map1  
                                IF (grid2_add_map1(l).eq.atm_add) THEN
                                   links=links+1
                                   row =MOD(grid1_add_map1(l),(ie_g-2))
                                   IF (row == 0) THEN
                                     row = (ie_g-2)+1
                                   ELSE
                                     row = row+1
                                   ENDIF
                                   column = INT(grid1_add_map1(l)/(ie_g-2))+1
                                   IF (row == (ie_g-2)+1) column=column-1
                                   XO2A(n)%o2a_tr(i,j)%om_ilon_g(links)=row
                                   XO2A(n)%o2a_tr(i,j)%om_ilat_g(links)=column
                                   XO2A(n)%o2a_tr(i,j)%weight(links)=wts_map1(1,l)
                                   IF (j==ngpblks.and.i.gt.npromz) XO2A(n)%o2a_tr(i,j)%weight(links)=0._dp 
                                ENDIF
                               end do
                              ENDDO
                            ENDDO
                          ENDIF
                      ENDIF
                    ENDDO
                    call remap_dealloc
                ENDIF
              ENDDO
              ! and distribute to all the fields
      ENDIF

      CALL info_bi ('transformation OCEAN --> ATMOSPHERE LOCKED !' )

    ENDIF !lrunsum

#ifdef MPIOM_13B
      DEALLOCATE(dlxp_g, dlyp_g       &
              ,dlxu_g, dlyu_g       &
              ,dlxv_g, dlyv_g)
#endif
      DEALLOCATE(AMSUE_G_L1,AMSUO_G_L1)


   END SUBROUTINE GRID_INTERPOLATION

!-----------------------------------------------------------------------

  SUBROUTINE check_mask

    USE messy_main_blather_bi, ONLY: info_bi
    USE messy_main_mpi_bi,     ONLY: p_parallel_io
    ! a2o MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)

    IMPLICIT NONE
! ---------------------------- Local declarations ----------------------
!
      INTEGER                 ::  ila_nneiadd  ! Nearest-neighbor address

      INTEGER                 :: &
          ib_dst,                &! INDEX loop for the distance grid
          ib_src,                &! INDEX loop for the source grid
          ib_neigh,              &! INDEX loop on the corresponding src pts
          ib_links                ! INDEX loop for the links      
 
      INTEGER                 :: &
          nb_Vmm,                &! number of Mtt points
          ntotland,              &! number of land points
          ntotoce,               &! number of oceanic points
          ntotngh                 ! number of corresponding source points
 
      INTEGER                 :: &
          beg_links,             &! begining of the serie of links
          nb_links                ! number of links
 
      REAL(DP) ::                &
          coslat,                &! cosinus of the latitude
          sinlat,                &! sinus of the latitude
          coslon,                &! cosinus of the longitude
          sinlon,                &! sinus of the longitude
          distance,              &
          dist_min

      LOGICAL :: lfirst_Vmm
 
      INTEGER   :: n


!        Treating Vmm points   V
!        -------------------  m m
!     The target point is a non-masked Valid point while the source points 
!         are all masked points. Use of 4 non-masked nearest neighbours.


      nb_Vmm = 0
      beg_links = 1

!  -- Loop all other the target points
      DO ib_dst = 1, grid2_size

      lfirst_Vmm = .TRUE.

!  -- If the point is a sea point
        IF (grid2_mask(ib_dst)) THEN

            ntotland = 0
            ntotoce = 0
            ntotngh = 0

            DO ib_links = 1, num_links_map1
              IF (grid2_add_map1(ib_links) .eq. ib_dst) THEN
                  ntotngh = ntotngh+1
                  IF (.not. grid1_mask(grid1_add_map1(ib_links))) THEN
                      ntotland = ntotland + 1
                  ELSE IF  (grid1_mask(grid1_add_map1(ib_links))) THEN
                      ntotoce = ntotoce + 1
                  ELSE
                      WRITE (*, *) 'Pb with ocean mask with Mtt 1'
                  END IF
              ENDIF
            END DO

!  -- If all the src points are land, treat it !
            IF (ntotland .EQ. ntotngh) THEN
              nb_Vmm = nb_Vmm + 1

              coslat = cos(grid2_center_lat(ib_dst))
              sinlat = sin(grid2_center_lat(ib_dst))
              coslon = cos(grid2_center_lon(ib_dst))
              sinlon = sin(grid2_center_lon(ib_dst))

              dist_min = bignum
              ila_nneiadd = 0
              DO ib_src = 1, grid1_size
                IF (grid1_mask(ib_src)) THEN
                    distance = acos(coslat*cos(grid1_center_lat(ib_src))* &
                        (coslon*cos(grid1_center_lon(ib_src)) +           &
                        sinlon*sin(grid1_center_lon(ib_src)))+            &
                        sinlat*sin(grid1_center_lat(ib_src)))
                    IF (distance < dist_min) THEN
                        ila_nneiadd = ib_src
                        dist_min = distance
                    ENDIF
                ENDIF
              END DO

              IF (ntotngh.eq.0) THEN
                num_links_map1  = num_links_map1 + 1
                call resize_remap_vars(1, 1)
                grid1_add_map1(num_links_map1) = ila_nneiadd
                grid2_add_map1(num_links_map1) = ib_dst
                wts_map1(1,num_links_map1) = 1.0_dp
              ELSE
                DO ib_links = 1, num_links_map1
                  IF (grid2_add_map1(ib_links) .eq. ib_dst) THEN
                     IF (lfirst_Vmm) THEN
                       lfirst_Vmm =.FALSE.
                       grid1_add_map1(ib_links) = ila_nneiadd
                       wts_map1(1,beg_links) = 1.0_dp
                     ELSE
                       grid1_add_map1(ib_links) = 0
                       wts_map1(1,beg_links) = 0.0_dp
                     ENDIF
                  ENDIF
                END DO
              ENDIF
            ENDIF
        ENDIF
      END DO

  END SUBROUTINE check_mask

!-----------------------------------------------------------------------
END MODULE messy_a2o_e5
