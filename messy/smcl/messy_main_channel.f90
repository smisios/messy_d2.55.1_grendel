! **********************************************************************
MODULE messy_main_channel
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

!!$#if (defined MPI) && (defined PNETCDF)
!!$#define ZPNETCDF
!!$#endif

  USE messy_main_constants_mem,      ONLY: DP, SP, I8       &
                                         , STRLEN_MEDIUM, STRLEN_ULONG &
                                         , STRLEN_LONG
  USE messy_main_channel_attributes, ONLY: t_attribute_list, add_attribute &
                                         , return_attribute, AF_NONE
  USE messy_main_channel_dimensions, ONLY: DIMID_UNDEF
  USE messy_main_channel_repr,       ONLY: IRANK, REPR_UNDEF &
                                         , t_representation &
                                         , get_representation
  USE messy_main_channel_mem,        ONLY: l_dom, n_dom, dom_unbound &
                                         , dom_current
  USE messy_main_tools,              ONLY: PTR_4D_ARRAY, split_name_domain
#ifdef HAVE_FORPY
  USE forpy_mod, ONLY: dict
#endif

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: HUGE, NULL, RESHAPE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'channel'
  CHARACTER(len=*), PARAMETER, PUBLIC :: modver = '2.4.3'

  PUBLIC :: DP, SP, I8, REPR_UNDEF, DIMID_UNDEF

  ! check with values in CDI library
  INTEGER, PARAMETER :: max_z_axes = 34
  INTEGER, PARAMETER :: CDI_UNDEFID= -1

  ! FILE-TYPES
  INTEGER, PARAMETER, PUBLIC :: FTYPE_UNDEFINED = -1
  INTEGER, PARAMETER, PUBLIC :: FTYPE_ASCII     = 1
!!$#ifdef ZPNETCDF
  INTEGER, PARAMETER, PUBLIC :: FTYPE_NETCDF    = 2
  INTEGER, PARAMETER, PUBLIC :: FTYPE_PNETCDF   = 3
!!$#else
!!$  INTEGER, PARAMETER, PUBLIC :: FTYPE_NETCDF    = 2
!!$!!$  INTEGER, PARAMETER, PUBLIC :: FTYPE_PNETCDF   = 3
!!$  INTEGER, PARAMETER, PUBLIC :: FTYPE_PNETCDF   = 2
!!$#endif
  INTEGER, PARAMETER, PUBLIC :: FTYPE_GRIB      = 4
  INTEGER, PARAMETER, PUBLIC :: FTYPE_HDF4      = 5
  INTEGER, PARAMETER, PUBLIC :: FTYPE_HDF5      = 6
  INTEGER, PARAMETER, PUBLIC :: FTYPE_CDI_NC    = 7
  INTEGER, PARAMETER, PUBLIC :: FTYPE_FORPY     = 8
  INTEGER, PARAMETER, PUBLIC :: FTYPE_MAXIMUM   = 8
  ! NOTE: DO NOT SET FTYPE_UNDEFINED AS DEFAULT
!!$#ifdef ZPNETCDF
!!$  INTEGER, PARAMETER, PUBLIC :: FTYPE_DEFAULT   = FTYPE_PNETCDF
!!$#else
  INTEGER, PARAMETER, PUBLIC :: FTYPE_DEFAULT   = FTYPE_NETCDF
!!$#endif
  CHARACTER(LEN=4), DIMENSION(FTYPE_MAXIMUM), PARAMETER, PUBLIC :: &
       FTYPE_EXT_TEXT = (/'.asc','.nc ','.nc ','.grb','.hdf','.hdf','.nc ' &
                         ,'    ' /)

  ! INPUT/OUTPUT MODES
  INTEGER, PARAMETER, PUBLIC :: IOMODE_OUT = 1  ! OUTPUT
  INTEGER, PARAMETER, PUBLIC :: IOMODE_RST = 2  ! RESTART
  INTEGER, PARAMETER, PUBLIC :: IOMODE_MAX = 2  ! NUMBER OF O-MODES
  CHARACTER(LEN=*), DIMENSION(IOMODE_MAX), PARAMETER, PUBLIC :: &
       IOMODE_TEXT = (/'output ', 'restart'/)
  ! ACCESS MODES
  INTEGER, PARAMETER, PUBLIC :: AMODE_WRITE = 1
  INTEGER, PARAMETER, PUBLIC :: AMODE_READ  = 2
  INTEGER, PARAMETER, PUBLIC :: AMODE_INIT  = 3

  ! PRIMARY AND SECONDARY DATA TYPE
  INTEGER, PARAMETER, PUBLIC :: SND_UNDEF = 0
  INTEGER, PARAMETER, PUBLIC :: SND_INS = 1  !!! DO NOT CHANGE POSITION OF _INS
  INTEGER, PARAMETER, PUBLIC :: SND_AVE = 2
  INTEGER, PARAMETER, PUBLIC :: SND_STD = 3
  INTEGER, PARAMETER, PUBLIC :: SND_MIN = 4
  INTEGER, PARAMETER, PUBLIC :: SND_MAX = 5
  INTEGER, PARAMETER, PUBLIC :: SND_CNT = 6
  INTEGER, PARAMETER, PUBLIC :: SND_CAV = 7
  INTEGER, PARAMETER, PUBLIC :: SND_MAXLEN = 7
  CHARACTER(LEN=4), DIMENSION(SND_MAXLEN, IOMODE_MAX), PARAMETER, PUBLIC :: &
       SND_TEXT  = RESHAPE ( &
       (/ '    ', '_ave', '_std', '_min', '_max', '_cnt', '_cav',     &
          '    ', '_x1 ', '_x2 ', '_min', '_max', '_cnt', '_csm' /),  &
       (/ SND_MAXLEN, IOMODE_MAX /) )

  ! ATTRIBUTE FLAGS FOR RESTART HANDLING
  INTEGER, PARAMETER, PUBLIC :: AF_RST_NONE = AF_NONE ! NOT REQUIRED
  INTEGER, PARAMETER, PUBLIC :: AF_RST_CMP  = 1       ! COMPARE
  INTEGER, PARAMETER, PUBLIC :: AF_RST_INP  = 2       ! INPUT

  ! netCDF IDs
  INTEGER, PARAMETER, PUBLIC :: NC_ID_UNDEF = -1

  ! STRING LENGTHS
  INTEGER, PARAMETER, PUBLIC :: STRLEN_CHANNEL = 2*STRLEN_LONG + 1
  INTEGER, PARAMETER, PUBLIC :: STRLEN_OBJECT  = 2*STRLEN_MEDIUM + 1 + 4

  ! ====================================================================
  ! CHANNEL OBJECTS
  ! ====================================================================
  ! FOR &CPL namelists
  TYPE t_chaobj_cpl
     CHARACTER(LEN=STRLEN_CHANNEL) :: cha  = ''
     CHARACTER(LEN=STRLEN_OBJECT ) :: obj  = ''
  END TYPE t_chaobj_cpl

  ! MEMORY MANAGEMENT
  TYPE t_channel_object_mem
     !
     ! MEMORY USAGE
     INTEGER(I8) :: usage       = 0_I8  ! primary data section
     INTEGER(I8) :: usage_2nd   = 0_I8  ! secondary data section
     !
     ! FLAGS FOR INTERNAL USE
     LOGICAL     :: lalloc  = .FALSE.   ! AUTOMATIC MEMORY ALLOCATION
  END TYPE t_channel_object_mem

  ! I/O MANAGEMENT
  TYPE t_channel_object_io
     !
     ! RESTART HANDLING
     LOGICAL :: lrestart      = .FALSE. ! OUTPUT TO RESTART FILE ?
     LOGICAL :: lignore       = .FALSE. ! IGNORE lrestreq ?
     !
     ! OUTPUT FLAGS
     LOGICAL, DIMENSION(SND_MAXLEN) :: lout = .FALSE.
     ! SPECIAL FOR CONDITIONAL COUNTER (CNT) / AVAERAGE (CAV)
     REAL(DP), DIMENSION(2)         :: range = &
          (/ -HUGE(0.0_DP), HUGE(0.0_DP) /)
     !
  END TYPE t_channel_object_io

  ! netCDF I/O (internal use only !)
  TYPE t_channel_object_netcdf
     ! variable ID
     INTEGER                        :: varid        = NC_ID_UNDEF
     ! dimension IDs
     INTEGER                        :: dimid(IRANK) = NC_ID_UNDEF
     ! IDs OF SECONDARY VARIABLES
     INTEGER, DIMENSION(:), POINTER :: svarid => NULL()
  END TYPE t_channel_object_netcdf

  TYPE t_channel_object_cdi
     ! CDI output?
     LOGICAL                        :: lout                = .FALSE.
     ! variable ID
     INTEGER                        :: varid               = CDI_UNDEFID
     ! grid ID
     INTEGER                        :: gridid              = CDI_UNDEFID
     ! zaxis ID
     INTEGER                        :: zaxisid             = CDI_UNDEFID
     !
     ! IDs OF SECONDARY VARIABLES
     INTEGER, DIMENSION(:), POINTER :: svarid => NULL()
  END TYPE t_channel_object_cdi

#ifdef HAVE_FORPY
  TYPE t_channel_object_forpy
     TYPE(dict), DIMENSION(SND_MAXLEN) :: d
  END type t_channel_object_forpy
#endif
  !
  ! +++ ADD OTHER OUTPUT FORMATS HERE

  TYPE t_channel_object_int
     !
     ! OUTPUT AND RESTART
     LOGICAL :: lout = .FALSE.     ! ANY OUTPUT ?
     LOGICAL :: lrst = .FALSE.     ! ANY RESTART ?
     LOGICAL :: lign = .FALSE.     ! IGNORE lrestreq ?
     LOGICAL :: lref = .FALSE.     ! IS reference ?
     ! EXPORT DATA ?
     LOGICAL, DIMENSION(SND_MAXLEN, IOMODE_MAX) :: lexp = .FALSE.
     ! ... MEMORY MANAGEMENT FOR PRIMARY AND SECONDARY DATA
     INTEGER, DIMENSION(SND_MAXLEN) :: i2nd = 0  ! INDEX IN 2ndary DATA
     INTEGER :: n2nd = 0  ! SECONDARY DATA DIMENSION
     !
     ! MISC
     ! - FIELD HAS BEEN SET FROM RESTART FILE
     LOGICAL, DIMENSION(SND_MAXLEN) :: lrestart_read = .FALSE.
     !
     ! netCDF
     TYPE(t_channel_object_netcdf), DIMENSION(IOMODE_MAX) :: netcdf
     !
     ! CDI
     TYPE(t_channel_object_cdi), DIMENSION(IOMODE_MAX) :: cdi
     !
#ifdef HAVE_FORPY
     TYPE(t_channel_object_forpy) :: forpy
#endif
     !
     ! +++ ADD OTHER OUTPUT FORMATS HERE
     !
  END TYPE t_channel_object_int

  TYPE t_channel_object
     CHARACTER(LEN=STRLEN_OBJECT)     :: name  = ''     ! NAME
     TYPE(t_attribute_list), POINTER  :: att => NULL()  ! OBJECT ATTRIBUTES
     TYPE(t_representation), POINTER  :: repr => NULL() ! REPRESENTATION
     TYPE(t_channel_object_mem)       :: memory         ! MEMORY MANAGEMENT
     TYPE(t_channel_object_io)        :: io             ! I/O MANAGEMENT
     ! ABSOLUTELY REQUIRED IN RESTART FILE ?
     LOGICAL                          :: lrestreq = .FALSE.
     ! object contents is not time dependent: lstatic = .FALSE.
     LOGICAL                          :: lstatic = .FALSE.
     ! object contents needs always to be ouput in double precision
     LOGICAL                          :: ldp     = .FALSE.
     ! ONLY instantaneous output makes sense (i.e. for additional attribute
     ! variales, which can no be defined as dimension variables
     LOGICAL                          :: linst   = .FALSE.
     !
     TYPE(t_channel_object_int)       :: int            ! FOR INTERNAL USE
     REAL(DP), DIMENSION(:,:,:,:),     POINTER :: DATA => NULL() ! DATA
     TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER :: sdat => NULL() ! 2ndary DATA
     ! POINTER TO REGION WITHOUT BOUNDARIES FOR I/O
     REAL(DP), DIMENSION(:,:,:,:),     POINTER :: ioptr => NULL()
  END TYPE t_channel_object

  TYPE t_channel_object_list
!     PRIVATE
     TYPE(t_channel_object)               :: this
     TYPE(t_channel_object_list), POINTER :: next => NULL()
  END TYPE t_channel_object_list
  ! ====================================================================

  ! ====================================================================
  ! CHANNELS
  ! ====================================================================
  TYPE t_channel_def
     !
     ! DEFAULT REPRESENTATION
     INTEGER                    :: reprid   = REPR_UNDEF
     ! DEFAULT: NOT REQUIRED IN RESTART
     LOGICAL                    :: lrestreq = .FALSE.
     ! DEFAULT OUTPUT FLAGS
     TYPE (t_channel_object_io) :: io
     !
  END TYPE t_channel_def

  ! netCDF I/O (internal use only !)
  TYPE t_channel_netcdf
     INTEGER  :: fileID     = NC_ID_UNDEF
     INTEGER  :: dimid_time = NC_ID_UNDEF
     INTEGER  :: io_pe = 0                 ! IO PE FOR THIS CHANNEL
  END TYPE t_channel_netcdf

#ifdef HAVE_FORPY
  TYPE t_channel_forpy
     INTEGER    :: io_pe     = 0        ! pe responsible for this channel
     TYPE(dict) :: ga                   ! storage for global attributes
     TYPE(dict) :: time                 ! time information
  END type t_channel_forpy
#endif

  TYPE t_channel_cdi
     INTEGER  :: fileID     = CDI_UNDEFID
     INTEGER  :: vlistID    = CDI_UNDEFID
     INTEGER  :: taxisID    = CDI_UNDEFID

     INTEGER  :: CellGridID = CDI_UNDEFID
     INTEGER  :: EdgeGridID = CDI_UNDEFID
     INTEGER  :: VertGridID = CDI_UNDEFID

     INTEGER  :: ZaxisID(max_z_axes) = CDI_UNDEFID
     !
     INTEGER  :: io_pe = 0 ! IO PE FOR THIS CHANNEL
  END TYPE t_channel_cdi
  !
  ! +++ ADD OTHER OUTPUT FORMATS HERE

  TYPE t_channel_io
     ! OUTPUT FILE TYPE
     INTEGER, DIMENSION(IOMODE_MAX)  :: ftype  = FTYPE_DEFAULT
     ! NO. OF TIME STEPS PER FILE (OUT)
     INTEGER                         :: ntpf   = 0
     !
  END TYPE t_channel_io

  TYPE t_channel_int
     !
     ! OUTPUT AND RESTART
     LOGICAL  :: lout      = .FALSE.         ! OUTPUT ? ANY OBJECT ?
     LOGICAL  :: lrst      = .FALSE.         ! RESTART ? ANY OBJECT ?
     LOGICAL  :: lrestreq  = .FALSE.         ! ANY OBJECT REQUIRED IN RESTART
     LOGICAL  :: lign      = .FALSE.         ! IGNORE ALL (!) lrestreq ?
     ! - OUTPUT FILENAME : <EXPERIMENT (15)>_YYYYMMDD_HHMM_<CHANNEL>.<ext>
     ! - RESTART FILENAME: restart_<CHANNEL>.<ext>
     ! <EXPERIMENT (14)>_YYYYMMDD_HHMMSSsss_<CHANNEL>_D<DOM>.<ext>
     CHARACTER(LEN=34+STRLEN_CHANNEL+8), DIMENSION(IOMODE_MAX) :: fname = ''
     !
     ! TIMER
     LOGICAL  :: lout_now   = .FALSE.        ! TIME MANAGER: OUTPUT NOW?
     LOGICAL  :: lrst_now   = .FALSE.        !               RESTART NOW?
     INTEGER  :: ntpfcnt    = 0              ! COUNTER OF TIME STEPS PER FILE
     LOGICAL  :: lnew_file  = .FALSE.        ! OPEN NEW FILE ?
     REAL(DP) :: tslo       = 0.0_DP         ! time [s] since last output
     ! FORCE OUTPUT/RESTART (to be set by set_channel_output)
     LOGICAL  :: lforce_out  = .FALSE.
     LOGICAL  :: lforce_rst  = .FALSE.
     ! SUPPRESS OUTPUT (to be set by set_channel_output)
     LOGICAL  :: lsuppr_out  = .FALSE.
     LOGICAL  :: lsuppr_rst  = .FALSE.
     ! FORCE NEW FILE
     LOGICAL  :: lforce_newfile = .FALSE.
     !
     !
     ! netCDF
     TYPE(t_channel_netcdf), DIMENSION(IOMODE_MAX) :: netcdf
     !
#ifdef HAVE_FORPY
     ! forpy: only for output, not for restart
     TYPE(t_channel_forpy)                         :: forpy
#endif
     ! CDI
     TYPE(t_channel_cdi), DIMENSION(IOMODE_MAX)    :: cdi
     !
     ! +++ ADD OTHER OUTPUT FORMATS HERE
     !
  END TYPE t_channel_int

  TYPE t_channel
     ! IDENTIFICATION
     CHARACTER(LEN=STRLEN_CHANNEL)    :: name  = ''       ! NAME
     INTEGER                          :: id    = 0        ! ID
     INTEGER                          :: dom_id = -1      ! DOMAIN
     TYPE(t_attribute_list), POINTER  :: att   => NULL()  ! CHANNEL ATTRIBUTES
     TYPE(t_channel_object_mem)       :: memory           ! MEMORY MANAGEMENT
     TYPE(t_channel_io)               :: io               ! I/O
     TYPE(t_channel_def)              :: default          ! OBJECT DEFAULTS
     ! INTERNAL
     TYPE(t_channel_int)              :: int              ! FOR INTERNAL USE
     !
     ! CHANNEL OBJECTS
     TYPE(t_channel_object_list), POINTER :: list => NULL()
  END TYPE t_channel

  TYPE t_channel_list
     ! PRIVATE
     TYPE(t_channel)               :: this
     TYPE(t_channel_list), POINTER :: next => NULL()
  END TYPE t_channel_list
  ! ====================================================================

  ! ====================================================================
  ! GLOBAL ATTRIBUTES COMMON TO ALL CHANNELS
  TYPE(t_attribute_list), POINTER, SAVE, PUBLIC :: GATT => NULL()
  ! CONCT. LIST OF CHANNELS / OBJECTS
  TYPE(t_channel_list),   POINTER, SAVE, PUBLIC :: GCHANNELLIST  => NULL()
  INTEGER,                         SAVE, PUBLIC :: NCHANNEL = 0
  ! CHANNELS ALREADY FIXATED
  LOGICAL,                         SAVE, PUBLIC :: LFIXATE = .FALSE.
  ! ====================================================================

  ! ====================================================================
  ! FOR CTRL-NAMELIST
  TYPE t_channel_object_ctrl
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = ''  ! CHANNEL NAME
     CHARACTER(LEN=STRLEN_OBJECT)  :: oname = ''  ! OBJECT NAME (TRACER !)
     TYPE(t_channel_object_io)     :: io
  END TYPE t_channel_object_ctrl
  !
  TYPE t_channel_ctrl
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = ''   ! CHANNEL NAME
     TYPE(t_channel_io)            :: cio
     TYPE(t_channel_object_io)     :: oio
  END TYPE t_channel_ctrl
  !
  TYPE t_channel_new
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = ''   ! CHANNEL NAME
  END TYPE t_channel_new
  !
  TYPE t_channel_ref_new
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname1 = ''  ! CHANNEL NAME  SOURCE
     CHARACTER(LEN=STRLEN_OBJECT)  :: oname1 = ''  ! OBJECT NAME SOURCE
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname2 = ''  ! CHANNEL NAME  NEW
     CHARACTER(LEN=STRLEN_OBJECT)  :: oname2 = ''  ! OBJECT NAME NEW
  END TYPE t_channel_ref_new
  !
  TYPE t_chaobj_att_new
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = ''   ! CHANNEL NAME
     CHARACTER(LEN=STRLEN_OBJECT)  :: oname = ''   ! OBJECT NAME
     CHARACTER(LEN=STRLEN_ULONG)   :: attname = '' ! ATRRIBUTE NAME
     INTEGER                       :: atttype = 0  ! STRING, INTEGER, REAL
     CHARACTER(LEN=STRLEN_ULONG)   :: attval  = '' ! ATRRIBUTE VALUE
     LOGICAL                       :: lforce  = .FALSE. ! OVERWRITE?
     LOGICAL                       :: lout    = .TRUE.  ! OUTPUT?
  END type t_chaobj_att_new
  !
  CHARACTER(LEN=14), PUBLIC, SAVE :: EXP_NAME = ''
  CHARACTER(LEN=STRLEN_ULONG), PUBLIC, SAVE :: EXEC_CHECKSUM = 'unknown'
  !
  LOGICAL, PUBLIC, SAVE :: L_FLUSH_IOBUFFER = .TRUE.
  ! 0: only error messages,
  ! 1: initialize and finalize
  ! 2: 1 and time loop
  ! 3: 2 with additional WARNINGS for DEBUGGING
  ! 4: 3 with additional expensive (CPU-time !) output for debugging
  INTEGER, PUBLIC, SAVE :: I_VERBOSE_LEVEL  = 2
  !
  TYPE(t_channel_ctrl),                                SAVE  :: OUT_DEFAULT
  INTEGER,                    DIMENSION(FTYPE_MAXIMUM),SAVE  :: OUT_PREC = &
       (/1, 1, 1, 1, 1, 1, 1, 1/)
  INTEGER, PARAMETER  :: NMAXCHANNELS =  700
  TYPE(t_channel_ctrl),        DIMENSION(NMAXCHANNELS),  SAVE  :: OUT_CHANNEL
  INTEGER, PARAMETER  :: NMAXOBJECTS  = 1000
  TYPE(t_channel_object_ctrl), DIMENSION(NMAXOBJECTS),   SAVE  :: OUT_OBJECT
  INTEGER, PARAMETER  :: NMAXADDCHANNEL = 50
  TYPE(t_channel_new),         DIMENSION(NMAXADDCHANNEL),SAVE  :: ADD_CHANNEL
  INTEGER, PARAMETER  :: NMAXADDREF    = 1100
  TYPE(t_channel_ref_new),     DIMENSION(NMAXADDREF),    SAVE  :: ADD_REF
  INTEGER, PARAMETER  :: NMAXADDATT    = 1000
  TYPE(t_chaobj_att_new),      DIMENSION(NMAXADDATT),    SAVE  :: ADD_ATT
  ! since modify_attributes is called twice, already applied attributes
  ! need to be skipped at 2nd call, otherwise lforce=.FALSE. will always
  ! cause an error ...
  LOGICAL, DIMENSION(NMAXADDATT),    SAVE  :: LADD_ATT = .FALSE.

  LOGICAL, PUBLIC, SAVE :: l_cheat = .FALSE.

  PUBLIC :: NMAXCHANNELS, NMAXOBJECTS, NMAXADDCHANNEL, NMAXADDREF
  PUBLIC :: NMAXADDATT
  PUBLIC :: OUT_DEFAULT, OUT_CHANNEL, OUT_OBJECT, ADD_CHANNEL, ADD_REF
  PUBLIC :: ADD_ATT
  PUBLIC :: OUT_PREC
  PUBLIC :: t_channel_list, t_channel                                 &
       , t_channel_int, t_channel_def, t_channel_io, t_channel_netcdf &
       , t_channel_cdi                                                &
       , t_channel_new, t_channel_ref_new, t_channel_ctrl
  PUBLIC :: t_channel_object_list, t_channel_object                   &
       , t_channel_object_int, t_channel_object_io, t_channel_object_netcdf &
       , t_channel_object_cdi                                         &
       , t_channel_object_mem, t_channel_object_ctrl
  PUBLIC :: t_chaobj_cpl
  PUBLIC :: t_chaobj_att_new
  ! ====================================================================

  ! PUBLIC INTERFACES ==================================================
  !
  ! - ATTRIBUTES
  INTERFACE new_attribute
     MODULE PROCEDURE add_attribute
     MODULE PROCEDURE new_global_attribute
     MODULE PROCEDURE new_channel_attribute
     MODULE PROCEDURE new_channel_object_attribute
  END INTERFACE
  PUBLIC :: new_attribute    ! add   GLOBAL, CHANNEL, or OBJECT ATTRIBUTES
  !
  INTERFACE get_attribute
     MODULE PROCEDURE return_attribute
     MODULE PROCEDURE get_global_attribute
     MODULE PROCEDURE get_channel_attribute
     MODULE PROCEDURE get_channel_object_attribute
  END INTERFACE
  PUBLIC :: get_attribute    ! add   GLOBAL, CHANNEL, or OBJECT ATTRIBUTES
  !
  INTERFACE write_attribute
     MODULE PROCEDURE write_global_attributes
     MODULE PROCEDURE write_channel_attributes
     MODULE PROCEDURE write_channel_object_attributes
  END INTERFACE
  PUBLIC :: write_attribute ! write GLOBAL, CHANNEL, or OBJECT ATTRIBUTES
  !
  ! - CHANNELS
  PUBLIC :: new_channel           ! define new CHANNEL
  INTERFACE write_channel
     MODULE PROCEDURE write_channel_by_name
     MODULE PROCEDURE write_channel_all
  END INTERFACE
  PUBLIC :: write_channel
  PUBLIC :: get_channel_info
  PUBLIC :: get_channel_name
  PUBLIC :: set_channel_output           ! force/suppress output
  PUBLIC :: get_channel_output           ! info on output
  PUBLIC :: set_channel_newfile          ! force new output file
  !
  ! - CHANNEL OBJECTS
  PUBLIC :: new_channel_object           ! define new channel OBJECT
  PUBLIC :: get_channel_object           ! get POINTER to OBJECT-memory
  PUBLIC :: get_channel_object_info      ! selected information
  PUBLIC :: new_channel_object_reference ! channel memory reference
  PUBLIC :: set_channel_object_restreq   ! set lrestreq = .TRUE.
  PUBLIC :: set_channel_object_inst      ! set int/inst = .TRUE.
  PUBLIC :: get_channel_object_dimvar    ! get dimension variables and units
  PUBLIC :: get_channel_object_dimvalue  ! get data of one dimension variable
  !                                      !  selected by axis
  PUBLIC :: get_channel_object_slice     ! get 2D or 3D (horizontal) Slices
  PUBLIC :: copy_channel_object_atts     ! copy all object attributes
  !
  ! - OVERALL SETUP
  PUBLIC :: main_channel_read_nml_ctrl
  !
  PUBLIC :: modify_attributes            ! CTRL ADD_ATT
  !
  PUBLIC :: fixate_channels              ! consistency check and 2ndary memory
  !                                      ! (once after initialization)
  PUBLIC :: trigger_channel_output       ! trigger next output
  !                                      ! (once every time step)
  PUBLIC :: update_channels              ! update 2ndary variables
  !                                      ! (-> see flag)
  PUBLIC :: clean_channels               ! cleanup memory
  !
  PUBLIC :: update_cell_methods_attribute  ! update CF specific attribute
  PUBLIC :: delete_cell_methods_attribute  ! delete CF specific attribute
  ! ====================================================================

  ! PRIVATE INTERFACES =================================================
  !
  ! - ATTRIBUTES
  ! NONE
  !
  ! - CHANNELS
  INTERFACE loc_channel                    ! locate pointer to channel
     MODULE PROCEDURE loc_channel_by_name  ! locate pointer to channel by NAME
     MODULE PROCEDURE loc_channel_by_id    ! locate pointer to channel by ID
  END INTERFACE
  PUBLIC :: loc_channel
  !PRIVATE :: write_channel_by_ptr
  !
  INTERFACE write_channel_object
     MODULE PROCEDURE write_channel_object_by_name
     MODULE PROCEDURE write_channel_object_by_ptr
  END INTERFACE
  !
  INTERFACE trigger_channel_output
     MODULE PROCEDURE trigger_channel_output_full
     MODULE PROCEDURE trigger_channel_output_part
  END INTERFACE trigger_channel_output
  ! - CHANNEL OBJECTS
  !PRIVATE :: loc_channel_object ! locate pinter to channel-object (by NAME)
  !PRIVATE :: write_channel_object
  ! ====================================================================

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! **********************************************************************
  ! ATTRIBUTES
  ! **********************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE new_global_attribute(status &
       , ganame, i, r, c, loverwrite, iflag, lout)

    USE messy_main_channel_attributes, ONLY: add_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: ganame
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag
    LOGICAL,                   INTENT(IN), OPTIONAL :: lout

    CALL add_attribute(status, GATT, TRIM(ganame), i, c, r &
         , loverwrite, iflag, lout)

  END SUBROUTINE new_global_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel_attribute(status, cname &
       , caname, i, r, c, loverwrite, iflag, dom_id, lout)

    USE messy_main_channel_attributes, ONLY: add_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: cname
    CHARACTER(LEN=*),          INTENT(IN)           :: caname
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag
    INTEGER,                   INTENT(IN), OPTIONAL :: dom_id
    LOGICAL,                   INTENT(IN), OPTIONAL :: lout

    ! LOCAL
    TYPE(t_channel), POINTER :: channel => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel, dom_id)
    IF (status /= 0) RETURN

    CALL add_attribute(status, channel%att, TRIM(caname), i, c, r &
         , loverwrite, iflag, lout)

  END SUBROUTINE new_channel_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel_object_attribute(status, cname, oname &
       , oaname, i, r, c, loverwrite, iflag, dom_id, lout)

    USE messy_main_channel_attributes, ONLY: add_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: cname
    CHARACTER(LEN=*),          INTENT(IN)           :: oname
    CHARACTER(LEN=*),          INTENT(IN)           :: oaname
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag
    INTEGER,                   INTENT(IN), OPTIONAL :: dom_id
    LOGICAL,                   INTENT(IN), OPTIONAL :: lout

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object => NULL()

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object, dom_id)
    IF (status /= 0) RETURN

    CALL add_attribute(status, object%att, TRIM(oaname), i, c, r &
         , loverwrite, iflag, lout)

  END SUBROUTINE new_channel_object_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_global_attributes(status)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status

    WRITE(*,*) '=== GLOBAL ATTRIBUTES: ============================='
    CALL write_attribute(status, GATT)
    IF (status /= 0) RETURN
    WRITE(*,*) '===================================================='

  END SUBROUTINE write_global_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_attributes(status, cname, dom_id)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel), POINTER :: channel => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel, dom_id)
    IF (status /= 0) RETURN

    CALL write_attribute(status, channel%att)

  END SUBROUTINE write_channel_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_object_attributes(status, cname, oname, dom_id)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    CHARACTER(LEN=*), INTENT(IN)  :: oname
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object => NULL()

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object, dom_id)
    IF (status /= 0) RETURN

    CALL write_attribute(status, object%att)

  END SUBROUTINE write_channel_object_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_global_attribute(status &
       , ganame, i, r, c, iflag, lout)

    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)           :: status
    CHARACTER(LEN=*),          INTENT(IN)            :: ganame
    INTEGER,                   INTENT(OUT), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(OUT), OPTIONAL :: c
    REAL(DP),                  INTENT(OUT), OPTIONAL :: r
    INTEGER,                   INTENT(OUT), OPTIONAL :: iflag
    LOGICAL,                   INTENT(OUT), OPTIONAL :: lout

    CALL return_attribute(status, GATT, TRIM(ganame), i, c, r, iflag, lout)

  END SUBROUTINE get_global_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_attribute(status, cname &
       , caname, i, c, r, iflag, dom_id, lout)

    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)           :: status
    CHARACTER(LEN=*),          INTENT(IN)            :: cname
    CHARACTER(LEN=*),          INTENT(IN)            :: caname
    INTEGER,                   INTENT(OUT), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(OUT), OPTIONAL :: c
    REAL(DP),                  INTENT(OUT), OPTIONAL :: r
    INTEGER,                   INTENT(OUT), OPTIONAL :: iflag
    INTEGER,                   INTENT(IN),  OPTIONAL :: dom_id
    LOGICAL,                   INTENT(OUT), OPTIONAL :: lout

    ! LOCAL
    TYPE(t_channel), POINTER :: channel => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel, dom_id)
    IF (status /= 0) RETURN

    CALL return_attribute(status, channel%att, TRIM(caname) &
         , i, c, r, iflag, lout)

  END SUBROUTINE get_channel_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object_attribute(status, cname, oname &
       , oaname, i, r, c, iflag, dom_id, lout)

    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: cname
    CHARACTER(LEN=*),          INTENT(IN)           :: oname
    CHARACTER(LEN=*),          INTENT(IN)           :: oaname
    INTEGER,                   INTENT(OUT), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(OUT), OPTIONAL :: c
    REAL(DP),                  INTENT(OUT), OPTIONAL :: r
    INTEGER,                   INTENT(OUT), OPTIONAL :: iflag
    INTEGER,                   INTENT(IN),  OPTIONAL :: dom_id
    LOGICAL,                   INTENT(OUT), OPTIONAL :: lout

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object => NULL()

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object, dom_id)
    IF (status /= 0) RETURN

    CALL return_attribute(status, object%att, TRIM(oaname) &
         , i, c, r, iflag, lout)

  END SUBROUTINE get_channel_object_attribute
  ! -------------------------------------------------------------------

  ! **********************************************************************
  ! CHANNELS
  ! **********************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel(status, cname, reprid, lrestreq, dom_id)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)          :: status
    CHARACTER(LEN=*), INTENT(IN)           :: cname    ! CHANNEL NAME
    INTEGER,          INTENT(IN), OPTIONAL :: reprid   ! REPRESENTATION ID
    LOGICAL,          INTENT(IN), OPTIONAL :: lrestreq ! REQUIRED IN RESTART
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel_list), POINTER :: ai     => NULL()
    TYPE(t_channel_list), POINTER :: ae     => NULL()
    TYPE(t_channel),      POINTER :: channel => NULL()
    INTEGER                       :: zstat
    CHARACTER(LEN=STRLEN_CHANNEL) :: zname

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    CALL loc_channel(zstat, GCHANNELLIST, cname, channel, dom_id)
    IF (zstat /= 3003) THEN  ! CHANNEL (NAME) DOES NOT EXIST (IS OK HERE !)
       IF (zstat == 0) THEN  ! CHANNEL EXISTS ALREADY
          status = 3002      ! CHANNEL EXISTS ALREADY
       ELSE
          status = zstat    ! ERROR
       END IF
       RETURN
    END IF

    ! GOTO END OF LIST
    ai => GCHANNELLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO

    ! ADD NEW
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(GCHANNELLIST)) THEN
       GCHANNELLIST => ai              ! SET POINTER TO FIRST OBJECT
    ELSE
       ae%next => ai                   ! SET NEXT POINTER OF LAST OBJECT
       !                               ! TO NEW OBJECT
    END IF

    ! SET VALUES
    ai%this%name = TRIM(ADJUSTL(cname))

    ! SET DEFAULT REPRESENTATION
    IF (PRESENT(reprid))   ai%this%default%reprid   = reprid
    ! SET DEFAULT REQUIREMENT FOR RESTART
    IF (PRESENT(lrestreq)) ai%this%default%lrestreq = lrestreq

    ! COUNT AND SET ID
    NCHANNEL   = NCHANNEL + 1
    ai%this%id = NCHANNEL

    ! SET DOMAIN
    IF (l_dom) THEN
       IF (PRESENT(dom_id)) THEN
          ai%this%dom_id = dom_id
       ELSE
          CALL split_name_domain(status, TRIM(cname), zname, ai%this%dom_id)
          IF (status /= 0) THEN
             ai%this%dom_id = dom_current
          ELSE
             ai%this%name = TRIM(ADJUSTL(zname))
          END IF
       END IF
    ELSE
       ai%this%dom_id = dom_current
    ENDIF

    status = 0

  END SUBROUTINE new_channel
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_by_ptr(status, channel)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,         INTENT(OUT) :: status
    TYPE(t_channel), POINTER     :: channel

    ! LOCAL
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel_object),      POINTER :: object

    IF (.NOT. ASSOCIATED(channel)) THEN
       status = 3006 ! CHANNEL POINTER NOT ASSOCIATED
       RETURN
    END IF

    WRITE(*,*) '========================================================================'
    WRITE(*,*) ' NAME          : ',channel%name
    WRITE(*,*) ' ID            : ',channel%id
    WRITE(*,*) ' DOMAIN        : ',channel%dom_id
    WRITE(*,*) ' OUT-FILE-TYPE : ',channel%io%ftype(IOMODE_OUT)
    IF (channel%io%ntpf > 0) THEN
       WRITE(*,*) ' STEPS PER FILE: ',channel%io%ntpf
    ELSE
       WRITE(*,*) ' STEPS PER FILE: ','(event triggered)'
    ENDIF
    WRITE(*,*) ' ANY OUTPUT  ? : ',channel%int%lout
    WRITE(*,*) ' RST-FILE-TYPE : ',channel%io%ftype(IOMODE_RST)
    WRITE(*,*) ' ANY RESTART ? : ',channel%int%lrst
    WRITE(*,*) ' IGNORE ?      : ',channel%int%lign
    WRITE(*,*) ' MEMORY        : ',channel%memory%usage,' + ' &
         , channel%memory%usage_2nd
    WRITE(*,*) ' ATTRIBUTES    : '
    CALL write_attribute(status, channel%att)
    IF (status /= 0) RETURN

    WRITE(*,*) '------------------------------------------------------------------------'
    WRITE(*,'(1x,a24,1x,a7,1x,a13,1x,a4,1x,a20)') &
         '                        ','--RST--','---OUTPUT----','REPR','-------MEMORY-------'
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','R','R','I','U','I','A','S','M','M','C','C'
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','E','E','G','S','N','V','T','I','A','N','A'
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','S','S','N','E','S','E','D','N','X','T','V'
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','T','T','O','R',' ',' ',' ',' ',' ',' ',' '
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','A','R','R',' ',' ',' ',' ',' ',' ',' ',' '
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','R','E','E',' ',' ',' ',' ',' ',' ',' ',' '
    WRITE(*,'(1x,a24,1x,11(a1,1x),a4,1x,a1,1x,2(a8,1x))') &
         'NAME                    ','T','Q',' ',' ',' ',' ',' ',' ',' ',' ',' ','REPR','M','   MEM_1','  MEM_02'

    WRITE(*,*) '------------------------------------------------------------------------'

    le => channel%list
    object_loop: DO
       IF (.NOT. ASSOCIATED(le)) EXIT

       object => le%this

       CALL write_channel_object_by_ptr(status, object)
       IF (status /= 0) RETURN

       le => le%next
    END DO object_loop

    WRITE(*,*) '========================================================================'

    status = 0

  END SUBROUTINE write_channel_by_ptr
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_by_name(status, cname, dom_id)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel), POINTER :: channel

    CALL loc_channel(status, GCHANNELLIST, cname, channel, dom_id)
    IF (status /= 0) RETURN

    CALL write_channel_by_ptr(status, channel)

  END SUBROUTINE write_channel_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_all(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER, INTENT(OUT)  :: status

    ! LOCAL
    TYPE(t_channel_list),         POINTER :: ls
    TYPE(t_channel),              POINTER :: channel

    WRITE(*,*) '+++ CHANNELS: ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    ! EMPTY LIST ?
    IF (.NOT. ASSOCIATED(GCHANNELLIST)) THEN

       WRITE(*,*) '*** CHANNEL LIST IS EMPTY ***'

    ELSE

       ls => GCHANNELLIST
       channel_loop: DO
          IF (.NOT. ASSOCIATED(ls)) EXIT

          channel => ls%this

          CALL write_channel_by_ptr(status, channel)
          IF (status /= 0) RETURN

          ls => ls%next
       END DO channel_loop

    END IF

    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    status = 0

  END SUBROUTINE write_channel_all
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_info(status, cname, LDIMS, LREPRS, ONAMES &
       , pick, dom_id)

    USE messy_main_channel_dimensions,      ONLY: NDIM
    USE messy_main_channel_repr,            ONLY: NREP, IRANK

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, PRESENT, TRIM

    ! I/O
    INTEGER,               INTENT(OUT)          :: status
    CHARACTER(LEN=*),      INTENT(IN)           :: cname   ! CHANNEL NAME
    LOGICAL, DIMENSION(:), POINTER,    OPTIONAL :: LDIMS   ! DIMENSION FLAGS
    LOGICAL, DIMENSION(:), POINTER,    OPTIONAL :: LREPRS  ! REPR. FLAGS
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER, OPTIONAL :: ONAMES
    CHARACTER(LEN=*),      INTENT(IN), OPTIONAL :: pick
    INTEGER,               INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel),             POINTER :: channel
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: i
    INTEGER                              :: nvarcnt
    INTEGER                              :: ipick
    LOGICAL                              :: zl

    ! INIT
    IF (PRESENT(LDIMS)) THEN
       IF (ASSOCIATED(LDIMS)) DEALLOCATE(LDIMS)
       ALLOCATE(LDIMS(NDIM))
       LDIMS(:) = .FALSE.
    END IF
    !
    IF (PRESENT(LREPRS)) THEN
       IF (ASSOCIATED(LREPRS)) DEALLOCATE(LREPRS)
       ALLOCATE(LREPRS(NREP))
       LREPRS(:) = .FALSE.
    END IF
    !
    IF (PRESENT(ONAMES)) THEN
       IF (ASSOCIATED(ONAMES)) DEALLOCATE(ONAMES)
       NULLIFY(ONAMES)
    END IF
    !
    IF (PRESENT(pick)) THEN
       SELECT CASE(TRIM(ADJUSTL(pick)))
          CASE('restart')
             ipick = 0
          CASE('output')
             ipick = 1
          CASE('all')
             ipick = 2
          CASE DEFAULT
       END SELECT
    ELSE
       ipick = 2    ! DEFAULT: 'all'
    END IF

    CALL loc_channel(status, GCHANNELLIST, cname, channel, dom_id)
    IF (status /= 0) RETURN

    nvarcnt = 0
    le => channel%list
    object_loop: DO
       IF (.NOT. ASSOCIATED(le)) EXIT

       object => le%this

       SELECT CASE(ipick)
          CASE(0) ! RESTART
             zl = object%int%lrst
          CASE(1) ! OUTPUT
             zl = object%int%lout
          CASE(2) ! ALL
             zl = .TRUE.
       END SELECT

       IF (PRESENT(LREPRS)) LREPRS(object%repr%id) = zl

       IF (PRESENT(LDIMS)) THEN
          DO i=1, IRANK
             IF (ASSOCIATED(object%repr%dim(i)%ptr)) &
                  LDIMS(object%repr%dim(i)%ptr%id) = zl
          END DO
       END IF

       IF (zl) nvarcnt = nvarcnt + 1

       le => le%next
    END DO object_loop

    IF (PRESENT(ONAMES)) THEN
       ALLOCATE(ONAMES(nvarcnt))
       ONAMES(:) = ''
       nvarcnt = 0
       le => channel%list
       object_loop2: DO
          IF (.NOT. ASSOCIATED(le)) EXIT

          object => le%this

          SELECT CASE(ipick)
          CASE(0) ! RESTART
             zl = object%int%lrst
          CASE(1) ! OUTPUT
             zl = object%int%lout
          CASE(2) ! ALL
             zl = .TRUE.
          END SELECT

          IF (zl) THEN
             nvarcnt = nvarcnt + 1
             ONAMES(nvarcnt) = TRIM(object%name)
          END IF

          le => le%next
       END DO object_loop2
    END IF

    status = 0

  END SUBROUTINE get_channel_info
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_name(status, id, cname, dom_id)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                       INTENT(OUT) :: status
    INTEGER,                       INTENT(IN)  :: id
    CHARACTER(LEN=STRLEN_CHANNEL), INTENT(OUT) :: cname
    ! we can also optionally return the domain of the channel
    INTEGER,                       INTENT(OUT), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel),                   POINTER :: channel  => NULL()

    ! INIT
    cname = ''

    CALL loc_channel_by_id(status, GCHANNELLIST, id, channel)
    IF (status /= 0) RETURN

    cname = TRIM(channel%name)
    IF (PRESENT(dom_id)) &
         & dom_id = channel%dom_id

  END SUBROUTINE get_channel_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_channel_output(status, cname, lout, lrst, dom_id)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    LOGICAL,          INTENT(IN)  :: lout
    LOGICAL,          INTENT(IN), OPTIONAL  :: lrst
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel),                   POINTER :: channel  => NULL()
    LOGICAL :: zlrst

    CALL loc_channel(status, GCHANNELLIST, cname, channel, dom_id)
    IF (status /= 0) RETURN

    IF (PRESENT(lrst)) THEN
       zlrst = lrst    ! .TRUE.:  force / suppress restart output
    ELSE
       zlrst = .FALSE. ! DEFAULT: force / suppress output
    END IF

    IF (zlrst) THEN
       ! force / suppress restart output
       IF (lout) THEN
          channel%int%lforce_rst = .TRUE.
       ELSE
          channel%int%lsuppr_rst = .TRUE.
       END IF
    ELSE
       ! force / suppress output
       IF (lout) THEN
          channel%int%lforce_out = .TRUE.
          channel%int%lsuppr_out = .FALSE.
       ELSE
          channel%int%lsuppr_out = .TRUE.
          channel%int%lforce_out = .FALSE.
       END IF
    END IF

  END SUBROUTINE set_channel_output
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_output(status, cname, lout, lstat, dom_id)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    LOGICAL,          INTENT(OUT) :: lout
    LOGICAL,          INTENT(OUT), OPTIONAL :: lstat
    INTEGER,          INTENT(IN),  OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel),                   POINTER :: channel  => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel, dom_id)
    IF (status /= 0) RETURN

    lout = channel%int%lout .AND. &
            (channel%int%lout_now .OR. channel%int%lforce_out) .AND. &
            (.NOT. channel%int%lsuppr_out)

    IF (PRESENT(lstat)) THEN
       ! .TRUE. IF AVE, STD, ... etc. are requested
       lstat = (channel%memory%usage_2nd > 0)
    END IF

  END SUBROUTINE get_channel_output
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_channel_newfile(status, cname, lnew, dom_id)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    LOGICAL,          INTENT(IN)  :: lnew
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel),                   POINTER :: channel  => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel, dom_id)
    IF (status /= 0) RETURN

    channel%int%lforce_newfile = lnew

  END SUBROUTINE set_channel_newfile
  ! -------------------------------------------------------------------

  ! **********************************************************************
  ! CHANNEL OBJECTS
  ! **********************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel_object(status, cname, oname &
       , p0, p1, p2, p3, p4, mem                     &
       , reprid                                      &
       , lrestreq                                    &
       , lstatic                                     &
       , ldp                                         &
       , linst                                       &
       , dom_id)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, PRESENT, TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)          :: status
    CHARACTER(LEN=*),             INTENT(IN)           :: cname ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)           :: oname ! OBJECT NAME
    REAL(DP),                     POINTER,    OPTIONAL :: p0    ! POINTER ...
    REAL(DP), DIMENSION(:),       POINTER,    OPTIONAL :: p1    ! ...
    REAL(DP), DIMENSION(:,:),     POINTER,    OPTIONAL :: p2    ! ...
    REAL(DP), DIMENSION(:,:,:),   POINTER,    OPTIONAL :: p3    ! ...
    REAL(DP), DIMENSION(:,:,:,:), POINTER,    OPTIONAL :: p4    ! ... TO MEMORY
    REAL(DP), DIMENSION(:,:,:,:), POINTER,    OPTIONAL :: mem   ! EXT. MEMORY
    !
    ! REPRESENTATION ID
    INTEGER,                      INTENT(IN), OPTIONAL :: reprid
    ! REQ. IN RESTART FILE ?
    LOGICAL,                      INTENT(IN), OPTIONAL :: lrestreq
    LOGICAL,                      INTENT(IN), OPTIONAL :: lstatic
    LOGICAL,                      INTENT(IN), OPTIONAL :: linst
    LOGICAL,                      INTENT(IN), OPTIONAL :: ldp
    INTEGER,                      INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel),             POINTER :: channel  => NULL()
    TYPE(t_channel_object),      POINTER :: object => NULL()
    TYPE(t_channel_object_list), POINTER :: ai      => NULL()
    TYPE(t_channel_object_list), POINTER :: ae      => NULL()
    INTEGER                              :: zstat   ! local status
    INTEGER                              :: i1,i2,i3,i4
    LOGICAL                              :: lshape

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    CALL loc_channel(status, GCHANNELLIST, TRIM(cname), channel, dom_id)
    IF (status /= 0) RETURN

    CALL get_channel_object(zstat, TRIM(cname), TRIM(oname), dom_id=dom_id)
    IF (zstat == 0) THEN
       IF (I_VERBOSE_LEVEL >= 3) THEN
          WRITE(*,*) '**** WARNING **** ', &
               'get_channel_object: CHANNEL - OBJECT ', &
               TRIM(cname),' - ',TRIM(oname),' EXISTS ALREADY!'
       END IF
       status = 3102    ! CHANNEL OBJECT EXISTS ALREADY
       RETURN
    END IF

    ! GOTO END OF LIST
    ai => channel%list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO

    ! ADD NEW OBJECT TO LIST
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(channel%list)) THEN
       channel%list => ai               ! SET POINTER TO FIRST OBJECT
    ELSE
       ae%next     => ai                ! SET NEXT POINTER OF LAST OBJECT
       !                                ! TO NEW OBJECT
    END IF

    object => ai%this

    ! SET NAME
    object%name = TRIM(oname)
    ! SET REPRESENTATION
    IF (PRESENT(reprid)) THEN
       CALL get_representation(status, reprid, object%repr)
       IF (status /= 0) RETURN
    ELSE
       ! CHECK DEFAULT REPRESENTATION OF CHANNEL
       IF (channel%default%reprid == REPR_UNDEF) THEN
          status = 3104 ! CHANNEL DEFAULT REPRESENTATION IS UNDEFINED
          RETURN
       ELSE
          CALL get_representation(status, channel%default%reprid, object%repr)
          IF (status /= 0) RETURN
       END IF
    END IF

    ! MEMORY MANAGEMENT
    object%memory%lalloc = (.NOT. PRESENT(mem))
    IF (object%memory%lalloc) THEN
       i1 = object%repr%ldimlen(1)
       i2 = object%repr%ldimlen(2)
       i3 = object%repr%ldimlen(3)
       i4 = object%repr%ldimlen(4)
       ALLOCATE(object%data(i1,i2,i3,i4), STAT=zstat)
       IF (zstat /= 0) THEN
          status = 1000 ! MEMORY ALLOCATION FAILED
          RETURN
       END IF
       object%data(:,:,:,:) = 0.0_DP
       ! OBJECT
       object%memory%usage = INT(SIZE(object%data), I8)
       ! CHANNEL
       channel%memory%usage = channel%memory%usage + object%memory%usage
    ELSE
       ! CHECK SIZE
       i1 = SIZE(mem, 1)
       i2 = SIZE(mem, 2)
       i3 = SIZE(mem, 3)
       i4 = SIZE(mem, 4)
       lshape = (i1 == object%repr%ldimlen(1)) .AND. &
                (i2 == object%repr%ldimlen(2)) .AND. &
                (i3 == object%repr%ldimlen(3)) .AND. &
                (i4 == object%repr%ldimlen(4))
       IF (.NOT. lshape) THEN
          status = 3105  ! SHAPE OF MEMORY NOT CONFORM WITH REPRESENTATION
          RETURN
       ELSE
          object%data => mem(:,:,:,:)
       END IF
    END IF

    IF (PRESENT(lstatic)) THEN
       object%lstatic = lstatic
       object%linst   = object%linst .OR. object%lstatic
    END IF
    IF (PRESENT(linst)) THEN
       object%linst   = linst .OR. object%linst
    END IF
    IF (PRESENT(ldp)) THEN
       object%ldp = ldp
    ELSE
       object%ldp = .FALSE. ! default
    ENDIF

    ! SET DEFAULTS
    ! AVAILABILITY IN RESTART FILE IS MANDATORY ?
    IF (PRESENT(lrestreq)) THEN
       ! USER DEFINED PER OBJECT
       object%lrestreq = lrestreq
    ELSE
       ! USER DEFINED DEFAULT PER CHANNEL (DEFAULT: .FALSE.)
       object%lrestreq = channel%default%lrestreq
    END IF
    !
    object%io = channel%default%io

    ! POINTER TO MEMORY FOR I/O (WITHIN POTENTIAL BOUNDARIES)
    IF (object%repr%bounds%lbounds) THEN
       ! IN CASE OF BOUNDARIES
       object%ioptr => object%data( &
            object%repr%bounds%nbounds(1)+1: &
            object%repr%ldimlen(1)-object%repr%bounds%nbounds(1), &
            object%repr%bounds%nbounds(2)+1: &
            object%repr%ldimlen(2)-object%repr%bounds%nbounds(2), &
            object%repr%bounds%nbounds(3)+1: &
            object%repr%ldimlen(3)-object%repr%bounds%nbounds(3), &
            object%repr%bounds%nbounds(4)+1: &
            object%repr%ldimlen(4)-object%repr%bounds%nbounds(4)  &
            )
    ELSE
       object%ioptr => object%data
    ENDIF

    ! SET POINTER
    CALL get_channel_object(status, TRIM(cname), TRIM(oname) &
         , p0, p1, p2, p3, p4, dom_id=dom_id)

  END SUBROUTINE new_channel_object
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object(status, cname, oname, p0, p1, p2, p3, p4 &
       , linner &
       , dom_id)

    USE messy_main_channel_repr,  ONLY: repr_getptr

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)          :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)           :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)           :: oname
    REAL(DP),                     POINTER,    OPTIONAL :: p0
    REAL(DP), DIMENSION(:),       POINTER,    OPTIONAL :: p1
    REAL(DP), DIMENSION(:,:),     POINTER,    OPTIONAL :: p2
    REAL(DP), DIMENSION(:,:,:),   POINTER,    OPTIONAL :: p3
    REAL(DP), DIMENSION(:,:,:,:), POINTER,    OPTIONAL :: p4
    LOGICAL,                      INTENT(IN), OPTIONAL :: linner
    INTEGER,                      INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel_object),       POINTER :: object => NULL()

    LOGICAL :: zlinner

    IF (PRESENT(linner)) THEN
       zlinner = linner
    ELSE
       zlinner = .FALSE.   ! default
    ENDIF

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object, dom_id)
    IF (status /= 0) THEN
       IF (I_VERBOSE_LEVEL >= 3) THEN
          WRITE(*,*) '**** WARNING **** ', &
               'get_channel_object: CHANNEL - OBJECT ', &
               TRIM(cname),' - ',TRIM(oname),' NOT FOUND!'
       END IF
       RETURN
    END IF

    ! return only inner part of object with bounds
    ! repr_getptr always returns the pointer p0...p4 to the entire input
    ! pointer (3rd argument). This means zlinner has no meaning if
    ! the reprsentation (2nd argument) has no boundaries. In case
    ! boundaries are present
    ! - and linner is requested: We use the already present ioptr
    !                            to the inner section.
    ! - and linner is not requested (i.e. we want the field incl bounds):
    !                            We use the full data pointer and force
    !                            repr_getptr to remap the indices to
    !                            e.g. (-1,...)

    IF (.NOT.zlinner) THEN
       CALL repr_getptr(status, object%repr, object%data,  p0, p1, p2, p3, p4 &
                      , .FALSE.)
    ELSE
       CALL repr_getptr(status, object%repr, object%ioptr, p0, p1, p2, p3, p4 &
                      , .TRUE.)
    END IF

    IF (status /= 0) RETURN

  END SUBROUTINE get_channel_object
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object_info(status, cname, oname &
       , lrestart_read , reprid, axis, nbounds, lstatic   &
       , dom_id)

    IMPLICIT NONE

    INTRINSIC :: PRESENT, TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)           :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: oname
    LOGICAL,                      INTENT(OUT), OPTIONAL :: lrestart_read
    INTEGER,                      INTENT(OUT), OPTIONAL :: reprid
    CHARACTER(LEN=IRANK),         INTENT(OUT), OPTIONAL :: axis
    INTEGER, DIMENSION(IRANK),    INTENT(OUT), OPTIONAL :: nbounds
    LOGICAL,                      INTENT(OUT), OPTIONAL :: lstatic
    INTEGER,                      INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel_object),       POINTER :: object => NULL()

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object, dom_id)
    IF (status /= 0) RETURN

    IF (PRESENT(lrestart_read)) &
         lrestart_read = object%int%lrestart_read(SND_INS)

    IF (PRESENT(reprid)) &
         reprid = object%repr%id

    IF (PRESENT(axis)) &
         axis = object%repr%axis

    IF (PRESENT(nbounds)) &
         nbounds(:) = object%repr%bounds%nbounds(:)

    IF (PRESENT(lstatic)) &
         lstatic = object%lstatic

  END SUBROUTINE get_channel_object_info
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel_object_reference(status, &
       cname1, oname1, cname2, oname2, lcopyatt, dom_id1, dom_id2)

    USE messy_main_channel_attributes, ONLY: copy_attribute_list

    IMPLICIT NONE

    INTRINSIC :: PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)          :: status
    CHARACTER(LEN=*), INTENT(IN)           :: cname1   ! ORIGINAL CHANNEL NAME
    CHARACTER(LEN=*), INTENT(IN)           :: oname1   ! ORIGINAL OBJECT NAME
    CHARACTER(LEN=*), INTENT(IN)           :: cname2   ! NEW CHANNEL NAME
    CHARACTER(LEN=*), INTENT(IN)           :: oname2   ! NEW OBJECT NAME
    LOGICAL,          INTENT(IN), OPTIONAL :: lcopyatt ! COPY ALL ATTRIBUTES?
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id1, dom_id2

    ! LOCAL
    TYPE(t_channel_object),       POINTER  :: object1 => NULL()
    TYPE(t_channel_object),       POINTER  :: object2 => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: mem
    LOGICAL                                :: zlcopyatt

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    ! INIT
    IF (PRESENT(lcopyatt)) THEN
       zlcopyatt = lcopyatt
    ELSE
       zlcopyatt = .FALSE.   ! DEFAULT
    END IF

    CALL loc_channel_object(status, TRIM(cname1), TRIM(oname1), object1 &
         , dom_id1)
    IF (status /= 0) RETURN

    ! SHARE PRIMARY MEMORY
    ! SAME REPRESENTATAION
    mem => object1%data(:,:,:,:)
    CALL new_channel_object(status, TRIM(cname2), TRIM(oname2), mem=mem &
         ,reprid = object1%repr%id, dom_id = dom_id2)
    IF (status /= 0) RETURN

    ! PRE-SET I/O SETTINGS
    CALL loc_channel_object(status, TRIM(cname2), TRIM(oname2), object2 &
         , dom_id2)
    IF (status /= 0) RETURN

    object2%lrestreq    = .FALSE.       ! SHARED MEMORY !!!
    !object2%io          = object1%io
    object2%io%lrestart = .FALSE.       ! SHARED MEMORY !!!

    object2%int%lref    = .TRUE.        ! REFERENCE

    object2%linst       =  object1%linst

    object2%lstatic     = object1%lstatic

    object2%ldp      = object1%ldp

    ! COPY ALL ATTRIBUTES
    IF (zlcopyatt) THEN
       CALL copy_attribute_list(status, object1%att, object2%att)
       IF (status /= 0) RETURN
    END IF

    ! ADD SPECIAL ATTRIBUTE
    CALL new_channel_object_attribute(status, TRIM(cname2), TRIM(oname2) &
         , 'REFERENCE_TO', c=TRIM(cname1)//': '//TRIM(oname1) &
         , loverwrite=.TRUE., dom_id = dom_id2, lout=.TRUE.  )

  END SUBROUTINE new_channel_object_reference
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_channel_object_restreq(status, cname, oname, dom_id)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)           :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: oname
    INTEGER,                      INTENT(IN),  OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel_object),       POINTER :: object => NULL()

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object, dom_id)
    IF (status /= 0) RETURN

    object%lrestreq = .TRUE.

  END SUBROUTINE set_channel_object_restreq
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_channel_object_inst(status, cname, oname, dom_id)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)           :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: oname
    INTEGER,                      INTENT(IN),  OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel_object),       POINTER :: object => NULL()

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object, dom_id)
    IF (status /= 0) RETURN

    object%linst = .TRUE.

  END SUBROUTINE set_channel_object_inst
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object_dimvar(status, cname, oname, dva, units &
       , levtype, dom_id)

    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG
    USE messy_main_tools,              ONLY: PTR_1D_ARRAY
    USE messy_main_channel_dimvar,     ONLY: t_dimvar, get_dimvar
    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM, ADJUSTL

    ! I/O
    INTEGER,                      INTENT(OUT)           :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: oname
    ! POINTER TO DIMENSION VARIABLE DATA ARRAYS ! INTENT(OUT)
    TYPE (PTR_1D_ARRAY), DIMENSION(:), POINTER          :: dva
    ! UNITS  ! INTENT(OUT)
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER, OPTIONAL :: units
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER, OPTIONAL :: levtype
    INTEGER,                      INTENT(IN), OPTIONAL  :: dom_id

    ! LOCAL
    INTEGER :: reprid
    TYPE(t_representation), POINTER  :: repr
    INTEGER :: i
    TYPE(t_dimvar),      POINTER :: dimvar => NULL()

    ! CLEAN FIRST
    IF (ASSOCIATED(dva)) DEALLOCATE(dva)
    NULLIFY(dva)

    IF (PRESENT(units)) THEN
       IF (ASSOCIATED(units)) DEALLOCATE(units)
       NULLIFY(units)
    END IF

    CALL get_channel_object_info(status, cname, oname, reprid=reprid &
         , dom_id=dom_id)
    IF (status /= 0) RETURN

    CALL get_representation(status, reprid, repr)
    IF (status /= 0) RETURN

    ! ALLOCATE
    ALLOCATE(dva(repr%rank))
    DO i=1, repr%rank
       NULLIFY(dva(i)%ptr)
    END DO

    IF (PRESENT(units)) THEN
       ALLOCATE(units(repr%rank))
       DO i=1, repr%rank
          units(i) = ''
       END DO
    END IF

    IF (PRESENT(levtype)) THEN
       IF (ASSOCIATED(levtype)) THEN
          DEALLOCATE(levtype); NULLIFY(levtype)
       END IF
       ALLOCATE(levtype(repr%rank))
       DO i=1, repr%rank
          levtype(i) = ''
       END DO
    END IF

    dimension_loop: DO i=1, repr%rank

       CALL get_dimvar(status, repr%dim(i)%ptr%var &
            , TRIM(ADJUSTL(repr%dim(i)%ptr%var%this%name)), dimvar)
       ! 953: DIMENSION VARIABLE DOES NOT EXIST (OK HERE!)
       IF ((status /= 953) .AND. (status /= 0)) RETURN

       dva(i)%ptr => dimvar%val(:)

       IF (PRESENT(units)) THEN
          CALL return_attribute(status, dimvar%att, 'units', c=units(i))
          ! 805: ATTRIBUTE DOES NOT EXIST (OK HERE!)
          IF ((status /= 805) .AND. (status /= 0)) RETURN
       ENDIF
       IF (PRESENT(levtype)) THEN
          CALL return_attribute(status, dimvar%att, 'levtype', c=levtype(i))
          ! 805: ATTRIBUTE DOES NOT EXIST (OK HERE!)
          IF ((status /= 805) .AND. (status /= 0)) RETURN
       ENDIF

    END DO dimension_loop

    status = 0

  END SUBROUTINE get_channel_object_dimvar
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object_dimvalue(status, cname, oname, data, axis &
       , unit, levtype, dom_id)

    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG
    USE messy_main_tools,              ONLY: PTR_1D_ARRAY
    USE messy_main_channel_dimvar,     ONLY: t_dimvar, get_dimvar
    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM, ADJUSTL

    ! I/O
    INTEGER,                     INTENT(OUT)          :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),            INTENT(IN)           :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),            INTENT(IN)           :: oname
    ! POINTER TO DIMENSION VARIABLE DATA ARRAYS ! INTENT(OUT)
    REAL(dp),    DIMENSION(:),   POINTER              :: data
    ! Character defining  axis string for dimvar
    CHARACTER(LEN=1),            INTENT(IN)           :: axis
    ! UNITS  ! INTENT(OUT)
    CHARACTER(LEN=STRLEN_ULONG),             OPTIONAL :: unit
    CHARACTER(LEN=STRLEN_ULONG),             OPTIONAL :: levtype
    INTEGER,                     INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    INTEGER                                            :: reprid
    INTEGER                                            :: i
    LOGICAL                                            :: lfound
    TYPE(t_representation),                   POINTER  :: repr   => NULL()
    TYPE (PTR_1D_ARRAY),         DIMENSION(:), POINTER :: dva    => NULL()
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: units  => NULL()
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: levtypes  => NULL()

    ! CLEAN FIRST
    IF (ASSOCIATED(data)) DEALLOCATE(data)
    NULLIFY(data)

    IF (PRESENT(unit)) unit = ''

    CALL get_channel_object_info(status, cname, oname, reprid=reprid &
         , dom_id=dom_id)
    IF (status /= 0) RETURN

    CALL get_representation(status, reprid, repr)
    IF (status /= 0) RETURN

    CALL get_channel_object_dimvar(status, cname, oname, dva, units &
       , levtypes, dom_id)
    IF (status /= 0) RETURN

    lfound = .FALSE.
    DO i = 1, repr%rank
       IF (repr%axis(i:i) == axis) THEN
          IF (lfound) THEN
             ! DIMVALUE: REQUESTED AXIS TWO TIMES AVAILABLE FOR CHANNEL OBJECT
             status = 959
             RETURN
          END IF
          lfound = .TRUE.
          ALLOCATE(data(SIZE(dva(i)%ptr)))
          data(:) = dva(i)%ptr(:)
          IF (PRESENT(unit)) unit = units(i)
          IF (PRESENT(levtype)) levtype = levtypes(i)
       END IF
    END DO

    IF (ASSOCIATED(dva)) THEN
       DO i=1,  repr%rank
          NULLIFY(dva(i)%ptr)
       END DO
       DEALLOCATE(dva)
       NULLIFY(dva)
    END IF
    IF (ASSOCIATED(units)) THEN
       DEALLOCATE(units)
       NULLIFY(units)
    END IF
    IF (ASSOCIATED(levtypes)) THEN
       DEALLOCATE(levtypes)
       NULLIFY(levtypes)
    END IF

   IF (.NOT. lfound) THEN
       status = 958 ! DIMVALUE: REQUESTED AXIS NOT AVAILABLE
    ELSE
       status = 0
    END IF

  END SUBROUTINE get_channel_object_dimvalue
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE copy_channel_object_atts(status, chsrc, obsrc, chdst, obdst)

    USE messy_main_channel_attributes, ONLY: copy_attribute_list

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: chsrc   ! source channel name
    CHARACTER(LEN=*), INTENT(IN)  :: obsrc   ! source object name
    CHARACTER(LEN=*), INTENT(IN)  :: chdst   ! destination channel name
    CHARACTER(LEN=*), INTENT(IN)  :: obdst   ! destination object name

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object_src => NULL()
    TYPE(t_channel_object), POINTER :: object_dst => NULL()

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    CALL loc_channel_object(status, TRIM(chsrc), TRIM(obsrc), object_src)
    IF (status /= 0) RETURN

    CALL loc_channel_object(status, TRIM(chdst), TRIM(obdst), object_dst)
    IF (status /= 0) RETURN

    ! COPY ALL ATTRIBUTES
    CALL copy_attribute_list(status, object_src%att, object_dst%att)
    IF (status /= 0) RETURN

    ! ADD SPECIAL ATTRIBUTE
    CALL new_channel_object_attribute(status, TRIM(chdst), TRIM(obdst) &
         , 'ATT_COPY_FROM', c=TRIM(chsrc)//'::'//TRIM(obsrc) &
         , loverwrite=.TRUE., lout=.TRUE.  )
    IF (status /= 0) RETURN

    status = 0

  END SUBROUTINE copy_channel_object_atts
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_object_by_ptr(status, object)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_channel_object), POINTER     :: object

    IF (.NOT. ASSOCIATED(object)) THEN
       status = 3107 ! CHANNEL OBJECT POINTER NOT ASSOCIATED
       RETURN
    END IF

    WRITE(*,'(1x,a24,1x,11(L1,1x),i4,1x,L1,1x,2(i8,1x))') &
         object%name    &
         , object%int%lrst                        & ! ANY RESTART ?
         , object%lrestreq                        & ! REQUIRED IN RESTART
         , object%int%lign                        & ! IGNORE lrestreq ?
         , object%io%lrestart                     & ! USER DEFINED
         , object%io%lout(:)                                 & ! WHICH OUTPUT
         , object%repr%id                                    &
         , object%memory%lalloc                              &
         , object%memory%usage, object%memory%usage_2nd

    CALL write_attribute(status, object%att)

    status = 0

  END SUBROUTINE write_channel_object_by_ptr
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_object_by_name(status, cname, oname, dom_id)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    CHARACTER(LEN=*), INTENT(IN)  :: oname
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object

    CALL loc_channel_object(status, cname, oname, object, dom_id)
    IF (status /= 0) RETURN

    CALL write_channel_object_by_ptr(status, object)

  END SUBROUTINE write_channel_object_by_name
  ! -------------------------------------------------------------------

  ! **********************************************************************
  ! ALL CHANNELS
  ! **********************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE fixate_channels(status, flag)

    USE messy_main_tools,         ONLY: match_wild &
                                      , domains_from_string
    USE messy_main_channel_error, ONLY: channel_error_str
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, ANY, SIZE, TRIM, INDEX

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: flag

    ! LOCAL
    TYPE(t_channel_list),        POINTER :: ls
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel),             POINTER :: channel
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: zstat
    INTEGER                              :: i,i1,i2,i3,i4
    INTEGER                              :: js, je
    INTEGER                              :: jsnd
    INTEGER                              :: m
    CHARACTER(LEN=STRLEN_CHANNEL)        :: s1 = ''
    CHARACTER(LEN=STRLEN_CHANNEL)        :: s2 = ''
    CHARACTER(LEN=STRLEN_OBJECT)         :: e1 = ''
    CHARACTER(LEN=STRLEN_OBJECT)         :: e2 = ''
    CHARACTER(LEN=STRLEN_VLONG)          :: errstr
    ! search pointer for wildcards in source objects
    LOGICAL                              :: l_wild_exist
    TYPE(t_channel_list),        POINTER :: s_ls
    TYPE(t_channel_object_list), POINTER :: s_le
    TYPE(t_channel),             POINTER :: s_channel
    TYPE(t_channel_object),      POINTER :: s_object
    CHARACTER(LEN=STRLEN_CHANNEL)        :: cname = ''
    INTEGER                              :: do, do1, do2
    LOGICAL                              :: lok, lallch
    INTEGER,       DIMENSION(:), POINTER :: domnums  => NULL() &
                                          , domnums2 => NULL()
    INTEGER                              :: nd, num, num2

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    SELECT CASE (flag)
    CASE (1)
    ! CREATE NEW CHANNELS
    DO js=1, NMAXADDCHANNEL
       IF (TRIM(ADJUSTL(ADD_CHANNEL(js)%cname)) == '') CYCLE
       IF (l_dom) THEN
          ! if user input is "D00" an unbound domain is created,
          ! the unbound domain has to be given explicitly (extend=.FALSE.)
          CALL domains_from_string(status &
               , TRIM(ADJUSTL(ADD_CHANNEL(js)%cname)), n_dom, num, cname &
               , dnums=domnums, extend=.FALSE.)
          DO nd=1,SIZE(domnums)
             CALL new_channel(status, TRIM(cname), dom_id=domnums(nd))
             IF (status /= 0) RETURN
          END DO
          DEALLOCATE(domnums)
          NULLIFY(domnums)
       ELSE
          CALL new_channel(status, TRIM(ADJUSTL(ADD_CHANNEL(js)%cname)))
          IF (status /= 0) RETURN
       END IF
    END DO

    ls => GCHANNELLIST
    channel_loop_1: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       ! -------------------------------------------
       ! CREATE NEW CHANNEL OBJECT REFERENCES
       ! -------------------------------------------
       DO je=1, NMAXADDREF
          IF (TRIM(ADJUSTL(ADD_REF(je)%cname1)) == '') CYCLE
          IF (TRIM(ADJUSTL(ADD_REF(je)%cname2)) == '') CYCLE

          IF (l_dom) THEN
             CALL domains_from_string(status &
                  , TRIM(ADJUSTL(ADD_REF(je)%cname1)), n_dom, num &
                  , s1, dnums=domnums, extend=.TRUE.)
             IF (ANY(domnums == channel%dom_id)) THEN
                do1 = channel%dom_id
             ELSE
                CYCLE
             END IF
          ELSE
             s1 = TRIM(ADJUSTL(ADD_REF(je)%cname1))
             do1 = dom_unbound
          ENDIF
          e1 = TRIM(ADJUSTL(ADD_REF(je)%oname1))
          IF (TRIM(e1) == '') CYCLE

          IF (l_dom) THEN
             CALL domains_from_string(status &
                  , TRIM(ADJUSTL(ADD_REF(je)%cname2)), n_dom, num2 &
                  , s2, dnums=domnums2, extend=.TRUE.)
             IF (ANY(domnums2 == channel%dom_id)) THEN
                do2 = channel%dom_id
             ELSE
                CYCLE
             END IF
          ELSE
             s2 = TRIM(ADJUSTL(ADD_REF(je)%cname2))
             do2 = dom_unbound
          ENDIF
          IF ( TRIM(ADJUSTL(ADD_REF(je)%oname2)) == '' ) THEN
             e2 = TRIM(e1)                          ! USE SAME OBJECT NAME
          ELSE
             e2 = TRIM(ADJUSTL(ADD_REF(je)%oname2)) ! NEW NAME
          END IF
          !
          ! WILDCARD-MATCH, BUT NO SELF REFERENCE WITH IDENTICAL OBJECT NAME
          lok = ( match_wild( TRIM(s2), TRIM(channel%name) ) .AND. &
                  ( .NOT. &
                    ( (TRIM(channel%name) == TRIM(s1)) .AND. &
                      (TRIM(e2) == TRIM(e1)) &
                    ) &
                  ) &
                )

          IF (lok) THEN
             ! enable wildcard matches in the source object name
             l_wild_exist = (INDEX(e1, '?') > 0) .OR. (INDEX(e1, '*') > 0)
             if_wc_ext: IF (l_wild_exist) THEN

                s_ls => GCHANNELLIST
                channel_loop2: DO
                   IF (.NOT. ASSOCIATED(s_ls)) EXIT

                   s_channel => s_ls%this
                   IF (TRIM(s_channel%name) == TRIM(s1)) THEN

                      IF (l_dom) THEN
                         ! for domain models, the domain must also match
                         IF (s_channel%dom_id /= do1) THEN
                            s_ls => s_ls%next
                            CYCLE
                         ENDIF
                      END IF

                      s_le => s_channel%list
                      object_loop2: DO
                         IF (.NOT. ASSOCIATED(s_le)) EXIT

                         s_object => s_le%this

                         IF (match_wild(TRIM(e1),TRIM(s_object%name))) THEN
                            !Note: when the source object name includes a
                            !      wildcard, the object names in the
                            !      destination channel will be the same, as no
                            !      specific name can be attributed for wildcard
                            !      strings. => e2 == s_object%name
                            CALL new_channel_object_reference(status       &
                                 , TRIM(s1), TRIM(s_object%name)           &
                                 , TRIM(channel%name), TRIM(s_object%name) &
                                 , lcopyatt = .TRUE.             &
                                 , dom_id1=do1             &
                                 , dom_id2=channel%dom_id)
                            SELECT CASE(status)
                            CASE(0)
                               ! OK: CONTINUE
                            CASE(3003, 3103)
                               ! 3003: SOURCE CHANNEL DOES NOT EXIST
                               ! 3103: SOURCE CHANNEL OBJECT DOES NOT EXIST
                               WRITE(*,*) '**** WARNING CALLING ', &
                                    'new_channel_object_reference: ' &
                                    , TRIM(s1),'(', do1, ') / ' &
                                    , TRIM(s_object%name),' / '&
                                    , TRIM(channel%name), '('&
                                    , channel%dom_id, ') / ' &
                                    , TRIM(s_object%name)
                               errstr = channel_error_str(status)
                               WRITE(*,*) TRIM(errstr)
                            CASE DEFAULT
                               ! SEVERE ERROR
                               WRITE(*,*) '**** ERROR CALLING ' &
                                    , 'new_channel_object_reference: ' &
                                    , TRIM(s1),'(', do1, ') / '&
                                    , TRIM(s_object%name),' / '&
                                    , TRIM(channel%name), '('&
                                    , channel%dom_id, ') / ' &
                                    , TRIM(s_object%name)
                               errstr = channel_error_str(status)
                               WRITE(*,*) TRIM(errstr)
                               RETURN
                            END SELECT
                         END IF ! match wild source object
                         s_le => s_le%next
                     END DO object_loop2
                   END IF
                   s_ls => s_ls%next
                END DO channel_loop2
             ELSE
                CALL new_channel_object_reference(status, &
                     TRIM(s1), TRIM(e1), TRIM(channel%name), TRIM(e2), &
                     lcopyatt = .TRUE., &
                     dom_id1=do1, dom_id2=channel%dom_id)
                SELECT CASE(status)
                CASE(0)
                   ! OK: CONTINUE
                CASE(3003, 3103)
                   ! 3003: SOURCE CHANNEL DOES NOT EXIST
                   ! 3103: SOURCE CHANNEL OBJECT DOES NOT EXIST
                   WRITE(*,*) &
                        '*** WARNING CALLING new_channel_object_reference: ' &
                        , TRIM(s1), '(', do1, ') / ', TRIM(e1), ' / '   &
                        , TRIM(channel%name), '(', channel%dom_id, ') / ' &
                        , TRIM(e2)
                   errstr = channel_error_str(status)
                   WRITE(*,*) TRIM(errstr)
                CASE DEFAULT
                   ! SEVERE ERROR
                   WRITE(*,*) &
                        '*** ERROR CALLING new_channel_object_reference: ' &
                        , TRIM(s1), '(', do1, ') / ', TRIM(e1), ' / '   &
                        , TRIM(channel%name), '(', channel%dom_id, ') / ' &
                        , TRIM(e2)
                   errstr = channel_error_str(status)
                   WRITE(*,*) TRIM(errstr)
                   RETURN
                END SELECT
             ENDIF if_wc_ext ! wildcard exists
          ENDIF

          IF (l_dom) THEN
             DEALLOCATE(domnums, domnums2)
             NULLIFY(domnums)
             NULLIFY(domnums2)
          END IF
       END DO

       ! -------------------------------------------
       ! (1) SET EVERYTHING TO DEFAULT VALUE
       ! -------------------------------------------
       channel%io          = OUT_DEFAULT%cio
       channel%default%io  = OUT_DEFAULT%oio

       ! -------------------------------------------
       ! (2) SET TO SPECIFIC, IF AVAIALABLE
       ! -------------------------------------------
       DO js=1, NMAXCHANNELS
          IF (TRIM(OUT_CHANNEL(js)%cname) == '') CYCLE
          IF (l_dom) THEN
             ! get list of domain extensions from namelist input and match
             ! when any domain number in the given list corresponds to domain
             ! id of actual channel
             CALL domains_from_string(status &
                  , TRIM(OUT_CHANNEL(js)%cname), n_dom, num, cname &
                  , dnums=domnums, extend=.TRUE.)
             IF (match_wild(TRIM(cname), TRIM(channel%name)) .AND. &
                  ANY(domnums == channel%dom_id)) THEN
                channel%io         = OUT_CHANNEL(js)%cio
                channel%default%io = OUT_CHANNEL(js)%oio
             END IF
             DEALLOCATE(domnums)
             NULLIFY(domnums)
          ELSE
             ! if no domain is given, match all domains
             IF (match_wild(TRIM(OUT_CHANNEL(js)%cname), TRIM(channel%name))) &
                  THEN
                channel%io         = OUT_CHANNEL(js)%cio
                channel%default%io = OUT_CHANNEL(js)%oio
             END IF
          END IF
       END DO

       ! -------------------------------------------
       ! (3) INITIALIZE INTERNAL SETUP
       ! -------------------------------------------
       channel%int%lign = .TRUE.  ! NEEDED FOR INITIALIZATION LOGICS BELOW

       le => channel%list
       object_loop_1: DO
          IF (.NOT. ASSOCIATED(le)) EXIT

          object => le%this

          ! -------------------------------------------
          ! (1) SET EVERYTHING TO DEFAULT VALUE
          ! -------------------------------------------
          object%io = channel%default%io

          ! -------------------------------------------
          ! (2) SET TO SPECIFIC, IF AVAIALABLE
          ! -------------------------------------------
          DO je=1, NMAXOBJECTS
             IF ( (TRIM(OUT_OBJECT(je)%cname) == '') .OR.      &
                  (TRIM(OUT_OBJECT(je)%oname) == '')) CYCLE
             IF ( match_wild(TRIM(OUT_OBJECT(je)%cname), TRIM(channel%name) ) &
                  .AND. &
                  match_wild(TRIM(OUT_OBJECT(je)%oname), TRIM(object%name)) ) &
                  THEN
                  object%io = OUT_OBJECT(je)%io
             END IF
          END DO
          le => le%next
       END DO object_loop_1

        ls => ls%next
    END DO channel_loop_1

    CASE(2)

       ls => GCHANNELLIST
       channel_loop: DO
          IF (.NOT. ASSOCIATED(ls)) EXIT

          channel => ls%this

          le => channel%list
          object_loop: DO
             IF (.NOT. ASSOCIATED(le)) EXIT

             object => le%this
          ! -------------------------------------------
          ! (3) INTERNAL SETUP
          ! -------------------------------------------
          ! CHANNEL OBJECT EXPORT
          !
          ! CORRECT FOR OBJECTS WITH ONLY INSTANTANEOUS OUTPUT
          IF (object%linst) THEN
             object%io%lout(1) = ANY(object%io%lout(:))
             object%io%lout(2:SND_MAXLEN) = .FALSE.
          END IF
          ! - OUTPUT
          !   (USER DEFINED)
          object%int%lexp(:,IOMODE_OUT) = object%io%lout(:)
          ! - ANY OUTPUT ?
          object%int%lout = ANY(object%int%lexp(:, IOMODE_OUT))
          !
          ! - RESTART
          !   (REQUESTED STATISTICS ALWAYS REQUIRES RESTART OUTPUT ...
          object%int%lexp(:,IOMODE_RST) = object%io%lout(:)
          !   ... INSTANTANEOUS VALUE ONLY ON REQUEST
          object%int%lexp(SND_INS,IOMODE_RST) = &
               object%io%lrestart .OR. & ! FORCED BY USER OR ...
               object%lrestreq           ! ... REQUIRED (CODE !)
          !   ... x2 (STD) REQUIRES x1 (AVE)
          object%int%lexp(SND_AVE,IOMODE_RST) = &
               object%int%lexp(SND_AVE,IOMODE_RST) .OR. &
               object%int%lexp(SND_STD,IOMODE_RST)
          !   ... csm (CAV) REQUIRES cnt (CNT)
          object%int%lexp(SND_CNT,IOMODE_RST) = &
               object%int%lexp(SND_CAV,IOMODE_RST) .OR. &
               object%int%lexp(SND_CNT,IOMODE_RST)
          ! - ANY RESTART ?
          object%int%lrst = ANY(object%int%lexp(:, IOMODE_RST))
          ! - IGNORE lrestreq
          object%int%lign = object%io%lignore !.OR. &
          !     ( (.NOT. object%lrestreq) .AND. (.NOT. object%io%lrestart) )

          ! - 2ndary MEMORY MANAGEMENT (1)
          DO jsnd=2, SND_MAXLEN
             IF (object%int%lexp(jsnd, IOMODE_RST)) THEN
                object%int%n2nd = object%int%n2nd + 1
                object%int%i2nd(jsnd) = object%int%n2nd
             END IF
          END DO
          !
          ! - 2ndary MEMORY MANAGEMENT (2)
          i1 = SIZE(object%data, 1)
          i2 = SIZE(object%data, 2)
          i3 = SIZE(object%data, 3)
          i4 = SIZE(object%data, 4)
          ALLOCATE(object%sdat(object%int%n2nd), STAT=zstat)
          IF (zstat /= 0) THEN
             status = 1000
             RETURN
          END IF
          DO i=1, object%int%n2nd
             ALLOCATE(object%sdat(i)%ptr(i1,i2,i3,i4), STAT=zstat)
             IF (zstat /= 0) THEN
                status = 1000
                RETURN
             END IF
             object%sdat(i)%ptr(:,:,:,:) = 0.0_DP
             object%memory%usage_2nd = &
                  object%memory%usage_2nd + &
                  INT( SIZE(object%sdat(i)%ptr), I8 )
          END DO

          ! CHANNEL
          channel%memory%usage_2nd = &
               channel%memory%usage_2nd + object%memory%usage_2nd
          ! - ANY OUTPUT ?
          channel%int%lout = channel%int%lout .OR. object%int%lout
          ! - ANY OBJECT REQUIRED IN RESTART ?
          channel%int%lrestreq = channel%int%lrestreq .OR. object%lrestreq
          ! - ANY RESTART ?
          channel%int%lrst = channel%int%lrst .OR. object%int%lrst
          ! - IGNORE ALL lrestreq ? (ONLY IF ALL ARE TRUE !!!)
          channel%int%lign = channel%int%lign .AND. object%int%lign
          ! -------------------------------------------

          ! -------------------------------------------
          ! OUTPUT I/O
          ! -------------------------------------------
          DO m = 1, IOMODE_MAX
             SELECT CASE(channel%io%ftype(m))
             CASE (FTYPE_NETCDF, FTYPE_PNETCDF)
                !
                ! netCDF
                !
                ALLOCATE(object%int%netcdf(m)%svarid(object%int%n2nd))
                object%int%netcdf(m)%svarid(:) = NC_ID_UNDEF
                !
                !
             CASE (FTYPE_CDI_NC)
                !
                ! CDI
                !
                ALLOCATE(object%int%cdi(m)%svarid(object%int%n2nd))
                object%int%cdi(m)%svarid(:) = CDI_UNDEFID
                !
                !
             CASE (FTYPE_FORPY)
                !
                ! FORPY
                ! force to trigger "new file" in every time step ...
                channel%io%ntpf = 1
                !
                !CASE(...)
                ! +++ ADD OTHER OUTPUT FORMATS HERE
             CASE DEFAULT
                !
             END SELECT
          END DO

          le => le%next
       END DO object_loop

       ! -------------------------------------------
       ! OUTPUT TIMER
       ! -------------------------------------------
       ! TRIGGER NEW FILE AT BEGINNING
       channel%int%lnew_file = .TRUE.

       ls => ls%next
    END DO channel_loop

    CALL modify_attributes(status)
    IF (status /= 0) RETURN

    ! SET FIXATED FLAG
    LFIXATE = .TRUE.

    END SELECT

    status = 0

  END SUBROUTINE fixate_channels
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE modify_attributes(status)

    USE messy_main_tools,              ONLY: domains_from_string &
                                           , str2num, match_wild
    USE messy_main_channel_attributes, ONLY: TYPE_INTEGER, TYPE_STRING &
                                           , TYPE_REAL_DP

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status

    ! LOCAL
    INTEGER :: i
    CHARACTER(LEN=STRLEN_CHANNEL)        :: cname  = ''
    CHARACTER(LEN=STRLEN_CHANNEL)        :: cxname = ''
    INTEGER                              :: do
    INTEGER                              :: ai
    REAL(DP)                             :: ar
    LOGICAL                              :: lc ! channel (T) or object (F)
    TYPE(t_channel_list), POINTER        :: ls  => NULL()
    INTEGER, DIMENSION(:), POINTER       :: domnums => NULL()
    INTEGER                              :: num
    LOGICAL                              :: lwild

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    attribute_loop: DO i=1, NMAXADDATT

       IF (TRIM(ADD_ATT(i)%cname) == '')   CYCLE
       IF (TRIM(ADD_ATT(i)%attname) == '') CYCLE

       ! attribute has already been applied at first call ...
       IF (LADD_ATT(i)) CYCLE

       lc = TRIM(ADD_ATT(i)%oname) == ''

       ! set domain to dom_unbound (for single domain models)
       ! in case of multi-domain model, it will be (re)set later in the
       ! channel loop
       do = dom_unbound
       IF (l_dom) THEN
          CALL domains_from_string(status &
               , TRIM(ADD_ATT(i)%cname), n_dom, num, cxname &
               , dnums=domnums, extend=.TRUE.)
       ELSE
          cxname = TRIM(ADD_ATT(i)%cname)
       END IF
       !
       ! To allow for wild_cards in the channel name, we have to loop
       ! here over all channels to find matches
       ls => GCHANNELLIST
       channel_loop: DO
          IF (.NOT. ASSOCIATED(ls)) EXIT
          IF (match_wild(TRIM(cxname),TRIM(ls%this%name)) ) THEN
             IF(l_dom) THEN
                ! cycle and check next channel if we cannot
                ! match domain in multi-domain model case
                IF (ANY(domnums==ls%this%dom_id)) THEN
                   do = ls%this%dom_id
                ELSE
                   ls => ls%next
                   CYCLE
                END IF
             END IF
             cname = TRIM(ls%this%name)
             SELECT CASE(ADD_ATT(i)%atttype)
                !
             CASE(TYPE_INTEGER)
                CALL str2num(TRIM(ADD_ATT(i)%attval), ai, status)
                IF (status /= 0) THEN
                   status = 0100
                   RETURN
                ENDIF
                IF (lc) THEN
                   CALL new_attribute(status, TRIM(cname) &
                        , TRIM(ADD_ATT(i)%attname), i=ai  &
                        , loverwrite = ADD_ATT(i)%lforce, dom_id=do &
                        , lout=ADD_ATT(i)%lout)
                ELSE
                   CALL new_attribute(status, TRIM(cname)&
                        , TRIM(ADD_ATT(i)%oname)         &
                        , TRIM(ADD_ATT(i)%attname), i=ai &
                        , loverwrite = ADD_ATT(i)%lforce, dom_id=do &
                        , lout=ADD_ATT(i)%lout)
                ENDIF
                !
             CASE(TYPE_STRING)
                !
                IF (lc) THEN
                   CALL new_attribute(status, TRIM(cname) &
                        , TRIM(ADD_ATT(i)%attname), c=TRIM(ADD_ATT(i)%attval) &
                        , loverwrite = ADD_ATT(i)%lforce, dom_id=do &
                        , lout=ADD_ATT(i)%lout)
                ELSE
                   CALL new_attribute(status, TRIM(cname) &
                        , TRIM(ADD_ATT(i)%oname) &
                        , TRIM(ADD_ATT(i)%attname), c=TRIM(ADD_ATT(i)%attval) &
                        , loverwrite = ADD_ATT(i)%lforce, dom_id=do &
                        , lout=ADD_ATT(i)%lout)
                ENDIF
                !
             CASE(TYPE_REAL_DP)
                !
                CALL str2num(TRIM(ADD_ATT(i)%attval), ar, status)
                IF (status /= 0) THEN
                   status = 0100
                   RETURN
                ENDIF
                IF (lc) THEN
                   CALL new_attribute(status, TRIM(cname) &
                        , TRIM(ADD_ATT(i)%attname), r=ar &
                        , loverwrite = ADD_ATT(i)%lforce, dom_id=do &
                        , lout=ADD_ATT(i)%lout)
                ELSE
                   CALL new_attribute(status, TRIM(cname)&
                        , TRIM(ADD_ATT(i)%oname)         &
                        , TRIM(ADD_ATT(i)%attname), r=ar &
                        , loverwrite = ADD_ATT(i)%lforce, dom_id=do &
                        , lout=ADD_ATT(i)%lout)
                ENDIF
                !
             CASE DEFAULT
                status = 806
                RETURN
             END SELECT
             !
             SELECT CASE(status)
             CASE(0, 3003, 3103)
                ! THIS IS OK HERE:
                ! 3003: CHANNEL DOES NOT EXIST
                ! 3103: CHANNEL OBJECT DOES NOT EXIST
             CASE DEFAULT
                RETURN
             END SELECT

             ! attribute has been successfully applied
             LADD_ATT(i) = (status == 0)
          END IF
          ls => ls%next
       END DO channel_loop
       IF (l_dom) THEN
          DEALLOCATE(domnums)
          NULLIFY(domnums)
       END IF
    END DO attribute_loop

    status = 0

  END SUBROUTINE modify_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trigger_channel_output_full(status, lnow, ltnf, lforce_new)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    LOGICAL, DIMENSION(:), INTENT(IN)  :: lnow       ! output now
    LOGICAL, DIMENSION(:), INTENT(IN)  :: ltnf       ! new file now
    LOGICAL,               INTENT(IN)  :: lforce_new ! new files for all now

    ! LOCAL
    TYPE(t_channel_list),         POINTER :: ls
    TYPE(t_channel),              POINTER :: channel
    LOGICAL                               :: lfull

    ! CHECKS
    IF (SIZE(lnow) /= NCHANNEL) THEN
       status = 3005
       RETURN
    END IF

    IF (SIZE(ltnf) /= NCHANNEL) THEN
       status = 3005
       RETURN
    END IF

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       ! SET TRIGGER FOR NEW FILE part C
       ! in case of nested domain this needs to be set in the farthest
       ! nested domain (smallest time step) for all domains, otherwise
       ! the force information is lost. Thus this part needed to be
       ! shifted before the CYCLE if dom_id /= dom_current
       channel%int%lnew_file = channel%int%lnew_file .OR. lforce_new

       IF (channel%dom_id /= dom_current) THEN
          ls => ls%next
          CYCLE
       END IF

       ! -------------------------------------------
       channel%int%lout_now = channel%int%lout .AND. &
            (lnow(channel%id) .OR. channel%int%lforce_out) .AND. &
            (.NOT. channel%int%lsuppr_out)

       ! CHECK IF NEW FILE IS REQUIRED BECAUSE OLD IS FULL
       lfull = .FALSE.
       ! INCREMENT COUNTER
       IF (channel%int%lout_now) THEN
          channel%int%ntpfcnt   = channel%int%ntpfcnt + 1
       END IF
       IF (channel%io%ntpf >0) THEN
          ! COUNTER
          lfull                 = (channel%int%ntpfcnt > channel%io%ntpf)
       ELSE
          ! EVENT TRIGGERED
          lfull = ltnf(channel%id)
       END IF
          ! SET TRIGGER FOR NEW FILE
          !  a) KEEP STATUS (RESET IN channel_finish_io (ONLY AFTER OUTPUT))
          !  b) OLD IS FULL
          !  c) FORCED (e.g., after restart) AND ANY OUTPUT REQUESTED
          !  d) FORCED BY USER (e.g. EVENT) AND ANY OUTPUT REQUESTED
       channel%int%lnew_file = channel%int%lnew_file .OR. &  ! a)
            lfull .OR.                                    &  ! b)
            (lforce_new .AND. channel%int%lout) .OR.      &  ! c)
            (channel%int%lforce_newfile .AND. channel%int%lout) ! d)
       ! RESET COUNTER IF NEW FILE IS TRIGGERED
       IF (channel%int%lnew_file) channel%int%ntpfcnt = 1

       ! RESET (only for one time step)
       channel%int%lforce_out = .FALSE.
       channel%int%lsuppr_out = .FALSE.
       channel%int%lforce_newfile = .FALSE.
       ! -------------------------------------------

       ! FORCE or SUPPRRESS restart output?
       channel%int%lrst_now = (channel%int%lrst .OR. channel%int%lforce_rst) &
            .AND. (.NOT. channel%int%lsuppr_rst)
       !  RESET (only for one time step)
       channel%int%lforce_rst = .FALSE.
       channel%int%lsuppr_rst = .FALSE.

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE trigger_channel_output_full
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trigger_channel_output_part(status, lnow)

     IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    LOGICAL, DIMENSION(:), INTENT(IN)  :: lnow       ! output now

    ! LOCAL
    TYPE(t_channel_list),         POINTER :: ls
    TYPE(t_channel),              POINTER :: channel

    ! CHECKS
    IF (SIZE(lnow) /= NCHANNEL) THEN
       status = 3005
       RETURN
    END IF

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       IF (channel%dom_id /= dom_current) THEN
          ls => ls%next
          CYCLE
       END IF

       channel%int%lout_now = channel%int%lout .AND. lnow(channel%id)

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE trigger_channel_output_part
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE update_channels(status, flag, dtime)

    USE messy_main_constants_mem, ONLY: FLAGGED_BAD

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, SQRT, MIN, MAX

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    INTEGER,  INTENT(IN)  :: flag
    REAL(DP), INTENT(IN)  :: dtime  ! time step length

    ! LOCAL
    TYPE(t_channel_list),        POINTER :: ls
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel),             POINTER :: channel
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: n, m
    REAL(DP)                             :: quot
    LOGICAL                              :: zlout_now
    INTEGER                              :: jsnd
    ! MUST BE TRUE IN FIRST STEP OF MODEL START AND IN FIRST STEP
    ! AFTER RESTART (lstart . OR. lresume)
    LOGICAL, SAVE                        :: lfirst = .TRUE.

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       IF (channel%dom_id /= dom_current) THEN
          ls => ls%next
          CYCLE
       END IF

       ! -------------------------------------------
       SELECT CASE(flag)
       CASE(1)
          ! ACCUMULATE FIELDS --------------------------------
          !
          ! TIME INTERVAL
          channel%int%tslo = channel%int%tslo + dtime
          !
       CASE(2)
          ! PREPARE FOR OUTPUT --------------------------------
          !
          IF (channel%int%lout_now) THEN
             IF (channel%int%tslo > 0.0_DP) THEN
                quot = 1.0_DP / channel%int%tslo
             ELSE
                quot = 1.0_DP
             END IF
          END IF
          !
       CASE(3)
          ! RESET AFTER OUTPUT --------------------------------
          !
          zlout_now = channel%int%lout_now  ! -> SAVE FOR OBJECTS BELOW
          IF (channel%int%lout_now) THEN
             ! TIME INTERVAL
             channel%int%tslo = 0.0_DP
             ! OUTPUT FLAG
             channel%int%lout_now = .FALSE.
          END IF
          !
       CASE DEFAULT
          ! ERROR
          !
          STATUS = 3106
          RETURN
          !
       END SELECT
       ! -------------------------------------------

       le => channel%list
       object_loop: DO
          IF (.NOT. ASSOCIATED(le)) EXIT

          object => le%this

          ! -------------------------------------------
          SELECT CASE(flag)
          CASE(1)
             ! ACCUMULATE FIELDS --------------------------------
             !
             ! AVERAGE: SUM X
             n = object%int%i2nd(SND_AVE)
             IF (n /= 0) THEN
                object%sdat(n)%ptr(:,:,:,:) = object%sdat(n)%ptr(:,:,:,:) &
                     + object%data(:,:,:,:) * dtime
             END IF
             !
             ! STANDARD DEVIATION: SUM X^2
             n = object%int%i2nd(SND_STD)
             IF (n /= 0) THEN
                object%sdat(n)%ptr(:,:,:,:) = object%sdat(n)%ptr(:,:,:,:) &
                  + object%data(:,:,:,:)**2 * dtime
             END IF
             ! MINIMUM
             n = object%int%i2nd(SND_MIN)
             IF (n /= 0) THEN
                IF (channel%int%tslo > dtime) THEN
                   ! in first step after start or restart ...
                   ! ... data from restart file must have been read ...
                   ! (if not, reset with object)
                   IF ( lfirst .AND. &
                        (.NOT. object%int%lrestart_read(SND_MIN)) ) THEN
                      object%sdat(n)%ptr(:,:,:,:) = object%data(:,:,:,:)
                   ELSE
                      object%sdat(n)%ptr(:,:,:,:) = &
                           MIN(object%sdat(n)%ptr(:,:,:,:), &
                               object%data(:,:,:,:))
                   END IF
                ELSE
                   ! LAST STEP WAS OUTPUT STEP => INITIALIZE
                   object%sdat(n)%ptr(:,:,:,:) = object%data(:,:,:,:)
                ENDIF
             END IF
             ! MAXIMUM
             n = object%int%i2nd(SND_MAX)
             IF (n /= 0) THEN
                IF (channel%int%tslo > dtime) THEN
                   ! in first step after start or restart ...
                   ! ... data from restart file must have been read ...
                   ! (if not, reset with object)
                   IF ( lfirst .AND. &
                        (.NOT. object%int%lrestart_read(SND_MAX)) ) THEN
                      object%sdat(n)%ptr(:,:,:,:) = object%data(:,:,:,:)
                   ELSE
                      object%sdat(n)%ptr(:,:,:,:) = &
                           MAX(object%sdat(n)%ptr(:,:,:,:), &
                               object%data(:,:,:,:))
                   END IF
                ELSE
                   ! LAST STEP WAS OUTPUT STEP => INITIALIZE
                   object%sdat(n)%ptr(:,:,:,:) = object%data(:,:,:,:)
                ENDIF
             END IF
             !
             ! CONDITION COUNTER
             n = object%int%i2nd(SND_CNT)
             IF (n /= 0) THEN
                WHERE( (object%data(:,:,:,:) >= object%io%range(1)) .AND. &
                     (object%data(:,:,:,:) <= object%io%range(2)) )
                       object%sdat(n)%ptr = object%sdat(n)%ptr + 1.0_DP
                END WHERE
             END IF
             !
             ! CONDITIONAL AVERAGE: SUM
             n = object%int%i2nd(SND_CAV)
             IF (n /= 0) THEN
                 WHERE( (object%data(:,:,:,:) >= object%io%range(1)) .AND. &
                        (object%data(:,:,:,:) <= object%io%range(2)) )
                        object%sdat(n)%ptr = object%sdat(n)%ptr + object%data
                 END WHERE
             END IF
             !
          CASE (2)
             ! PREPARE FOR OUTPUT --------------------------------
             IF (channel%int%lout_now) THEN
                ! AVERAGE
                n = object%int%i2nd(SND_AVE)
                m = n
                IF (n > 0) THEN
                   object%sdat(n)%ptr(:,:,:,:) = &
                        object%sdat(n)%ptr(:,:,:,:) * quot
                END IF
                ! STANDARD DEVIATION
                n = object%int%i2nd(SND_STD)
                IF (n > 0) THEN
                   object%sdat(n)%ptr(:,:,:,:) = &
                        SQRT( ABS (object%sdat(n)%ptr(:,:,:,:)*quot   &
                                 - object%sdat(m)%ptr(:,:,:,:)**2 ) )
                END IF
                ! CONDITIONAL AVERAGE
                n = object%int%i2nd(SND_CAV)
                m = object%int%i2nd(SND_CNT)
                IF (n > 0) THEN
                   WHERE (object%sdat(m)%ptr(:,:,:,:) > 0.0_DP)
                      object%sdat(n)%ptr = &
                           object%sdat(n)%ptr / &
                           object%sdat(m)%ptr
                   ELSEWHERE
                      !object%sdat(n)%ptr = REAL(-HUGE(0.0_SP),DP)
                      object%sdat(n)%ptr = FLAGGED_BAD
                   END WHERE
                END IF
             END IF
             !
          CASE (3)
             ! RESET AFTER OUTPUT --------------------------------
             !
             lfirst = .FALSE.  ! must be only .TRUE. in first time step after
             !                 ! start or restart
             !
             IF (zlout_now) THEN
                DO jsnd=2, SND_MAXLEN
                   n = object%int%i2nd(jsnd)
                   IF (n /= 0) THEN
                      IF ((jsnd == SND_MIN) .OR. (jsnd == SND_MAX)) THEN
                         ! SND_MIN and SND_MAX must not be reset to zero;
                         ! set to actual value instead. Two cases can happen
                         ! (note that restart-files are written later than
                         !  the output files):
                         ! 1) no restart output in this time step:
                         !    With channel%int%tslo > dtime
                         !    the value will be re-set (again) to object%data
                         !    in CASE(1) above.
                         ! 2) restart output in this time step:
                         !    The actual object value will be saved in the
                         !    restart file. After restart, CASE(1) above
                         !    will evaluate a=MIN(a,a) and b=MAX(b,b),
                         !    respectively.
                         object%sdat(n)%ptr(:,:,:,:) = object%data(:,:,:,:)
                      ELSE
                         ! SND_AVE, SND_STD, SND_CAV, SND_CNT: reset to zero
                         object%sdat(n)%ptr(:,:,:,:) = 0.0_DP
                      END IF
                      ! content does not longer come from restart file
                      ! -> (re)start calculation with actual value
                      !    (see CASE(1) above)
                      object%int%lrestart_read(jsnd) = .FALSE.
                   END IF
                END DO
                !
             END IF
             !
          CASE DEFAULT
             ! ERROR
             !
             STATUS = 3106
             RETURN
             !
          END SELECT
          ! -------------------------------------------

          le => le%next
       END DO object_loop

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE update_channels
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE clean_channels(status)

    USE messy_main_channel_attributes, ONLY: clean_attribute_list

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,  INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_channel_list),         POINTER :: lsai, lsae
    TYPE(t_channel_object_list),  POINTER :: leai, leae
    TYPE(t_channel),              POINTER :: channel
    TYPE(t_channel_object),       POINTER :: object
    INTEGER                               :: zstat
    INTEGER                               :: n, m

    status = 0
    zstat = 0

    lsai => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(lsai)) EXIT

       channel => lsai%this

       lsae => lsai
       lsai => lsai%next

       leai => channel%list

       object_loop: DO
          IF (.NOT. ASSOCIATED(leai)) EXIT

          object => leai%this

          leae => leai
          leai => leai%next

          ! -------------------------------------------
          ! ATTRIBUTES
          CALL clean_attribute_list(status, object%att)
          IF (status /= 0) RETURN
          !
          ! MEMORY: DATA
          !
          IF (object%memory%lalloc) THEN
             object%memory%usage = object%memory%usage &
                  - INT( SIZE(object%data), I8 )
             channel%memory%usage = channel%memory%usage &
                  - INT( SIZE(object%data), I8 )
             DEALLOCATE(object%data, STAT=zstat)
             IF (zstat /= 0) THEN
                status = 1001 ! MEMORY DEALLOCATION FAILED
                RETURN
             END IF
             IF (object%memory%usage /= 0_I8) THEN
                status = 3108 ! CHANNEL OBJECT PRIMARY MEMORY ERROR
                RETURN
             END IF
          END IF
          !
          ! MEMORY: SDAT
          DO n=1, object%int%n2nd
             object%memory%usage_2nd = object%memory%usage_2nd &
                  - INT(SIZE(object%sdat(n)%ptr), I8)
             channel%memory%usage_2nd = channel%memory%usage_2nd &
                  - INT(SIZE(object%sdat(n)%ptr), I8)
             DEALLOCATE(object%sdat(n)%ptr, STAT=zstat)
             IF (zstat /= 0) THEN
                status = 1001 ! MEMORY DEALLOCATION FAILED
                RETURN
             END IF
          END DO

          IF (ASSOCIATED(object%sdat)) THEN
             DEALLOCATE(object%sdat, STAT=zstat)
             NULLIFY(object%sdat)
          END IF
          IF (zstat /= 0) THEN
             status = 1001 ! MEMORY DEALLOCATION FAILED
             RETURN
          END IF

          IF (object%memory%usage_2nd /= 0_I8) THEN
             status = 3109 ! CHANNEL OBJECT SECONDARY MEMORY ERROR
             RETURN
          END IF
          ! -------------------------------------------

          ! -------------------------------------------
          ! IO
          ! -------------------------------------------
          DO m=1, IOMODE_MAX
             ! netCDF
             IF (ASSOCIATED(object%int%netcdf(m)%svarid)) &
                  DEALLOCATE(object%int%netcdf(m)%svarid)
             !
             ! CDI
             IF (ASSOCIATED(object%int%cdi(m)%svarid)) &
                  DEALLOCATE(object%int%cdi(m)%svarid)
             !
             ! +++ ADD OTHER OUTPUT FORMATS HERE
             !
          END DO

          DEALLOCATE(leae)
          NULLIFY(leae)

       END DO object_loop

       ! -------------------------------------------
       ! MEMORY: DATA
       IF (channel%memory%usage /= 0_I8) THEN
          status = 3108 ! CHANNEL SECONDARY MEMORY ERROR
          RETURN
       END IF
       ! MEMORY: SDAT
       IF (channel%memory%usage_2nd /= 0_I8) THEN
          status = 3109 ! CHANNEL SECONDARY MEMORY ERROR
          RETURN
       END IF
       !
       ! ATTRIBUTES
       CALL clean_attribute_list(status, channel%att)
       IF (status /= 0) RETURN
       !
       ! LIST OBJECT
       DEALLOCATE(lsae)
       NULLIFY(lsae)
       !
       ! COUNT
       NCHANNEL = NCHANNEL - 1
       ! -------------------------------------------

    END DO channel_loop

    NULLIFY(GCHANNELLIST)

    IF (NCHANNEL /= 0) THEN
       status = 1004
       RETURN
    END IF

    status = 0

  END SUBROUTINE clean_channels
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_read_nml_ctrl(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ EXP_NAME, EXEC_CHECKSUM, L_FLUSH_IOBUFFER, I_VERBOSE_LEVEL &
         , ADD_CHANNEL, ADD_REF, OUT_DEFAULT, OUT_PREC, OUT_CHANNEL &
         , OUT_OBJECT, ADD_ATT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status
    INTEGER                     :: js

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    ! FOR SOME STRANGE REASONS THIS IS REQUIRED FOR LF8.1
    DO js=1, NMAXADDCHANNEL
       ADD_CHANNEL(js)%cname = ''
    END DO
    L_FLUSH_IOBUFFER = .TRUE.
    I_VERBOSE_LEVEL = 1

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_channel_read_nml_ctrl
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_channel_by_name(status, list, name, channel, dom_id)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, LEN_TRIM, TRIM

    ! I/O
    INTEGER,              INTENT(OUT) :: status
    TYPE(t_channel_list), POINTER     :: list
    CHARACTER(LEN=*),     INTENT(IN)  :: name
    TYPE(t_channel),      POINTER     :: channel
    INTEGER,              INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel_list), POINTER :: ai  => NULL()
    TYPE(t_channel_list), POINTER :: ae  => NULL()
    LOGICAL                       :: lex
    CHARACTER(LEN=STRLEN_CHANNEL) :: zname
    INTEGER                       :: jg

    zname = ADJUSTL(TRIM(name))

    IF (l_dom) THEN
       IF (PRESENT(dom_id)) THEN
          jg = dom_id
       ELSE
          CALL split_name_domain(status, TRIM(name), zname, jg)
          IF (status /= 0) jg = dom_current
       END IF
    ELSE
       jg = dom_current
    END IF

    ! INIT
    lex = .FALSE.
    NULLIFY(channel)

    ! CHECKS
    IF (TRIM(zname) == '') THEN
       status = 3010 ! CHANNEL NAME IS EMPTY
       RETURN
    END IF
    IF (LEN_TRIM(ADJUSTL(name)) > STRLEN_CHANNEL) THEN
       status = 3001  ! CHANNEL NAME TOO LONG
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 3003  ! CHANNEL (NAME) DOES NOT EXIST
       RETURN
    END IF

    ! CHECK, IF IT EXISTS
    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF ((TRIM(ADJUSTL(zname)) == TRIM(ai%this%name)) .AND. &
            & (ai%this%dom_id == jg)) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO

    IF (lex) THEN
       channel => ai%this
    ELSE
       status = 3003  ! CHANNEL DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_channel_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_channel_by_id(status, list, id, channel)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,              INTENT(OUT)          :: status
    TYPE(t_channel_list), POINTER              :: list
    INTEGER,              INTENT(IN)           :: id
    TYPE(t_channel),      POINTER              :: channel

    ! LOCAL
    TYPE(t_channel_list), POINTER :: ai  => NULL()
    TYPE(t_channel_list), POINTER :: ae  => NULL()
    LOGICAL                       :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(channel)

    ! CHECKS
    IF (id <= 0) THEN
       status = 3007  ! INVALID CHANNEL ID
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 3004  ! CHANNEL (ID) DOES NOT EXIST
       RETURN
    END IF

    ! CHECK, IF IT EXISTS
    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF (id == ai%this%id) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO

    IF (lex) THEN
       channel => ai%this
    ELSE
       status = 3003  ! CHANNEL DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_channel_by_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_channel_object(status, cname, oname, object, dom_id)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, LEN_TRIM, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)   :: status
    CHARACTER(LEN=*), INTENT(IN)    :: cname
    CHARACTER(LEN=*), INTENT(IN)    :: oname
    TYPE(t_channel_object), POINTER :: object
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id

    ! LOCAL
    TYPE(t_channel),             POINTER :: channel => NULL()
    TYPE(t_channel_object_list), POINTER :: ai     => NULL()
    TYPE(t_channel_object_list), POINTER :: ae     => NULL()
    LOGICAL                              :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(object)

    CALL loc_channel(status, GCHANNELLIST, TRIM(cname), channel, dom_id)
    IF (status /= 0) RETURN

    ! CHECKS
    IF (TRIM(oname) == '') THEN
       status = 3110 ! CHANNEL OBJECT NAME IS EMPTY
    END IF
    IF (LEN_TRIM(ADJUSTL(oname)) > STRLEN_OBJECT) THEN
       status = 3101   ! CHANNEL OBJECT NAME TOO LONG
       RETURN
    END IF

    ! CHECK, IF NAME EXISTS
    ai => channel%list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF (TRIM(ADJUSTL(oname)) == TRIM(ai%this%name)) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO

    IF (lex) THEN
       object => ai%this
    ELSE
       status = 3103     ! CHANNEL OBJECT DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_channel_object
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object_slice(status, cname, oname, p2, p3 &
       , linner, dom_id, zslice, nslice, tslice)

    USE messy_main_channel_repr,  ONLY: repr_getptr

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)          :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)           :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)           :: oname
    REAL(DP), DIMENSION(:,:),     POINTER,    OPTIONAL :: p2
    REAL(DP), DIMENSION(:,:,:),   POINTER,    OPTIONAL :: p3
    LOGICAL,                      INTENT(IN), OPTIONAL :: linner
    INTEGER,                      INTENT(IN), OPTIONAL :: dom_id
    INTEGER,                      INTENT(IN), OPTIONAL :: zslice
    INTEGER,                      INTENT(IN), OPTIONAL :: nslice
    INTEGER,                      INTENT(IN), OPTIONAL :: tslice

    ! LOCAL
    TYPE(t_channel_object),       POINTER :: object  => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: lp4 => NULL()

    LOGICAL :: zlinner
    INTEGER :: ix, zid, nid, tid, sid
    LOGICAL :: lz, ln, lt

    IF (PRESENT(linner)) THEN
       zlinner = linner
    ELSE
       zlinner = .FALSE.   ! default
    ENDIF

    IF (PRESENT(zslice)) THEN
       IF (zslice < 0) THEN
          lz = .FALSE.
       ELSE
          lz = .TRUE.
       END IF
    ELSE
       lz = .FALSE.
    END IF
    IF (PRESENT(nslice)) THEN
       IF (nslice < 0) THEN
          ln = .FALSE.
       ELSE
          ln = .TRUE.
       END IF
    ELSE
       ln = .FALSE.
    END IF
    IF (PRESENT(tslice)) THEN
       IF (tslice < 0) THEN
          lt = .FALSE.
       ELSE
          lt = .TRUE.
       END IF
    ELSE
       lt = .FALSE.
    END IF

    IF (lz .AND. ln .AND. lt) THEN
       status = 2054 ! three diminisher are too much!
       RETURN
    END IF

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object, dom_id)
    IF (status /= 0) THEN
       IF (I_VERBOSE_LEVEL >= 3) THEN
          WRITE(*,*) '**** WARNING **** ', &
               'reduce_channel_object: CHANNEL - OBJECT ', &
               TRIM(cname),' - ',TRIM(oname),' NOT FOUND!'
       END IF
       RETURN
    END IF

    IF (object%repr%rank < 3) THEN
       status = 2050 ! RANK OF CHANNEL OBJECT TOO SMALL!
       RETURN
    END IF

    ! return only inner part of object with bounds
    ! repr_getptr always returns the pointer p0...p4 to the entire input
    ! pointer (3rd argument). This means zlinner has no meaning if
    ! the reprsentation (2nd argument) has no boundaries. In case
    ! boundaries are present
    ! - and linner is requested: We use the already present ioptr
    !                            to the inner section.
    ! - and linner is not requested (i.e. we want the field incl bounds):
    !                            We use the full data pointer and force
    !                            repr_getptr to remap the indices to
    !                            e.g. (-1,...)


    IF (.NOT.zlinner) THEN
       CALL repr_getptr(status, object%repr, object%data, p4=lp4 &
            , linner=.FALSE.)
    ELSE
       CALL repr_getptr(status, object%repr, object%ioptr, p4=lp4 &
            , linner=.TRUE.)
    END IF

    zid = -99
    nid = -99
    tid = -99
    sid = -99
    DO ix = 1, IRANK
       IF (object%repr%axis(ix:ix) == 'Z') THEN
          zid = ix
       ELSE IF (object%repr%axis(ix:ix) == 'N') THEN
          nid = ix
       ELSE IF (object%repr%axis(ix:ix) == 'T') THEN
          tid = ix
       ELSE IF (object%repr%axis(ix:ix) == '-') THEN
          sid = ix
       END IF
    END DO

    IF (lz .AND. zid < 0 ) THEN
       status = 2051 ! OBJECT DOES NOT HAVE A Z-AXIS!'
       RETURN
    END IF
    IF (ln .AND. nid < 0 ) THEN
       status = 2052 !OBJECT DOES NOT HAVE A N-AXIS!'
       RETURN
    END IF
    IF (lt .AND. tid < 0 ) THEN
       status = 2053 !OBJECT DOES NOT HAVE A T-AXIS!'
       RETURN
    END IF

    IF ( (lz .AND. ln) .or. (lt .AND. lz) .or. (lt .AND. ln)) THEN
       IF (PRESENT(p3)) THEN
          status = 2045 ! REDUCTION BY TWO RANKS implies 2D PTR
          RETURN
       END IF
       IF (PRESENT(p2)) THEN
          IF (.NOT. lt) CALL set_2D_ptr(status,p2, lp4, zid, nid, zslice &
               , nslice)
          IF (.NOT. lz) CALL set_2D_ptr(status,p2, lp4, tid, nid, tslice &
               , nslice)
          IF (.NOT. ln) CALL set_2D_ptr(status,p2, lp4, zid, tid, zslice &
               , tslice)
          IF (status /=0) RETURN
       END IF
    ELSE IF (lz) THEN
       IF (PRESENT(p2)) THEN
          CALL set_2D_ptr(status,p2, lp4, zid, sid, zslice,1)
          IF (status /=0) RETURN
       END IF
       IF (PRESENT(p3)) THEN
          CALL set_3D_ptr(status,p3, lp4, zid, zslice)
          IF (status /=0) RETURN
       END IF
    ELSE IF (ln) THEN
       IF (PRESENT(p2)) THEN
          CALL set_2D_ptr(status,p2, lp4, sid, nid, 1, nslice)
          IF (status /=0) RETURN
       END IF
       IF (PRESENT(p3)) THEN
          CALL set_3D_ptr(status,p3, lp4, nid, nslice)
          IF (status /=0) RETURN
       END IF
    ELSE IF (lt) THEN
       IF (PRESENT(p2)) THEN
          CALL set_2D_ptr(status,p2, lp4, sid, tid, 1, tslice)
          IF (status /=0) RETURN
       END IF
       IF (PRESENT(p3)) THEN
          CALL set_3D_ptr(status,p3, lp4, tid, tslice)
          IF (status /=0) RETURN
       END IF
    END IF
    IF (status /= 0) RETURN

    ! -----------------------------------------------------
    CONTAINS

      SUBROUTINE set_2D_ptr(status, p2, p4, id1, id2, sl1, sl2)

        IMPLICIT NONE

        INTEGER, INTENT(OUT)                  :: status
        REAL(dp), DIMENSION(:,:),     POINTER :: p2
        REAL(dp), DIMENSION(:,:,:,:), POINTER :: p4
        INTEGER, INTENT(IN)                   :: id1
        INTEGER, INTENT(IN)                   :: id2
        INTEGER, INTENT(IN)                   :: sl1
        INTEGER, INTENT(IN)                   :: sl2

        IF (id1 <= 0 .OR. id2 < 0) THEN
           status = 2060
           RETURN
        END IF

        IF (id1 == 1) THEN
           IF (id2 == 2) THEN
              p2 => p4(sl1,sl2,:,:)
           ELSE IF (id2 == 3) THEN
              p2 => p4(sl1,:,sl2,:)
           ELSE IF (id2 == 4) THEN
              p2 => p4(sl1,:,:,sl2)
           END IF
        ELSE IF  (id1 == 2) THEN
           IF (id2 == 1) THEN
              p2 => p4(sl2,sl1,:,:)
           ELSE IF (id2 == 3) THEN
              p2 => p4(:,sl1,sl2,:)
           ELSE IF (id2 == 4) THEN
              p2 => p4(:,sl1,:,sl2)
           END IF
        ELSE IF  (id1 == 3) THEN
           IF (id2 == 1) THEN
              p2 => p4(sl2,:,sl1,:)
           ELSE IF (id2 == 2) THEN
              p2 => p4(:,sl2,sl1,:)
           ELSE IF (id2 == 4) THEN
              p2 => p4(:,:,sl1,sl2)
           END IF
        ELSE IF  (id1 == 4) THEN
           IF (id2 == 1) THEN
              p2 => p4(sl2,:,:,sl1)
           ELSE IF (id2 == 2) THEN
              p2 => p4(:,sl2,:,sl1)
           ELSE IF (id2 == 3) THEN
              p2 => p4(:,:,sl2,sl1)
           END IF
        ENDIF

        status = 0

      END SUBROUTINE set_2D_ptr

      ! ------------------------------------------------

      SUBROUTINE set_3D_ptr(status, p3, p4, id, slice)

        IMPLICIT NONE

        INTEGER, INTENT(OUT)                  :: status
        REAL(dp), DIMENSION(:,:,:),   POINTER :: p3
        REAL(dp), DIMENSION(:,:,:,:), POINTER :: p4
        INTEGER, INTENT(IN)                   :: id
        INTEGER, INTENT(IN)                   :: slice

        IF (id <= 0) THEN
           status = 2060
           RETURN
        END IF

        IF (id == 1) THEN
           p3 => p4(slice,:,:,:)
        ELSE IF (id == 2) THEN
           p3 => p4(:,slice,:,:)
        ELSE IF (id == 3) THEN
           p3 => p4(:,:,slice,:)
        ELSE
           p3 => p4(:,:,:,slice)
        END IF

      END SUBROUTINE set_3D_ptr


  END SUBROUTINE get_channel_object_slice
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE update_cell_methods_attribute(status, list, jsnd, lstatic)

    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG
    USE messy_main_channel_attributes, ONLY: return_attribute &
                                           , add_attribute    &
                                           , t_attribute_list

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)              :: status
    TYPE(t_attribute_list),  POINTER  :: list
    INTEGER, INTENT(IN)               :: jsnd
    LOGICAL, INTENT(IN)               :: lstatic

    ! LOCAL
    CHARACTER(LEN=STRLEN_ULONG)       :: str
    CHARACTER(LEN=STRLEN_ULONG)       :: str_ori

    str     = ''
    str_ori = ''

    CALL return_attribute(status, list, 'cell_methods', c=str_ori)
    IF (status == 805 ) THEN
       str_ori = ''
    ELSE IF (status /=0 ) THEN
       RETURN
    END IF

    str = TRIM(ADJUSTL(str_ori))

    SELECT CASE(jsnd)
    CASE(SND_INS)
       IF (str(1:4) /= 'time') THEN
          ! OK
          IF (.NOT. lstatic) THEN
             str = 'time: point '//TRIM(ADJUSTL(str_ori))
          END IF
       END IF
    CASE(SND_AVE)
       str = 'time: mean '//TRIM(ADJUSTL(str_ori))
    CASE(SND_STD)
       str = 'time: standard_deviation '//TRIM(ADJUSTL(str_ori))
    CASE(SND_MIN)
       str = 'time: minimum '//TRIM(ADJUSTL(str_ori))
    CASE(SND_MAX)
       str = 'time: maximum '//TRIM(ADJUSTL(str_ori))
    CASE(SND_CNT)
       str = 'time: event_count '//TRIM(ADJUSTL(str_ori))
    CASE(SND_CAV)
       str = 'time: event_mean '//TRIM(ADJUSTL(str_ori))
    END SELECT

    IF (TRIM(str) /= '') THEN
       CALL add_attribute (status, list, 'cell_methods' &
            , c=str, loverwrite=.TRUE.)
       IF (status /=0) THEN
          RETURN
       END IF
    END IF
    IF (TRIM(str_ori) /= '') THEN
       CALL add_attribute (status, list, 'cell_methods_grid' &
            , c=TRIM(str_ori), loverwrite=.TRUE., lout=.FALSE.)
       IF (status /=0) RETURN
    END IF

    status = 0

  END SUBROUTINE update_cell_methods_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE delete_cell_methods_attribute(status, list)

    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG

    USE messy_main_channel_attributes, ONLY: return_attribute &
                                           , add_attribute    &
                                           , delete_attribute &
                                           , t_attribute_list

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)              :: status
    TYPE(t_attribute_list),  POINTER  :: list

    ! LOCAL
    CHARACTER(LEN=STRLEN_ULONG) :: str

    str = ''

    CALL delete_attribute(status, list, 'cell_methods')
    IF (status == 805) THEN ! attribute does not exist
       status = 0
       RETURN
    ELSE IF (status /= 0) THEN
       RETURN
    END IF

    CALL return_attribute(status, list, 'cell_methods_grid', c=str)
    IF (status == 805 ) THEN
       status = 0
       RETURN
    ELSE IF (status /=0 ) THEN
       RETURN
    END IF

    IF (TRIM(str) /= '') THEN
       CALL add_attribute (status, list, 'cell_methods' &
            , c=TRIM(str), loverwrite=.TRUE.)
       IF (status /=0) RETURN

       CALL delete_attribute(status, list, 'cell_methods_grid')
    END IF

  END SUBROUTINE delete_cell_methods_attribute
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel
! **********************************************************************
