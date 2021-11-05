! **********************************************************************
MODULE messy_main_channel_io
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

  USE messy_main_channel

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: NULL

  PUBLIC :: initialize_parallel_io ! FOR PARALLEL I/O
  PUBLIC :: channel_init_restart ! INITIALIZE DATE FROM RESTART
  PUBLIC :: channel_init_io      ! OPEN FILE FOR READ/WRITE
  PUBLIC :: channel_write_header ! INITIALIZE OUTPUT FILE (WRITE HEADER)
  PUBLIC :: channel_write_time   ! WRITE TIME INFORMATION
  PUBLIC :: channel_write_data   ! WRITE DATA TO OUTPUT FILE
  PUBLIC :: channel_finish_io    ! FLUSH / CLOSE FILE
  !
  PUBLIC :: channel_read_data    ! READ DATA FROM RESTART FILE
  !
  PUBLIC :: channel_dist_io      ! DISTRIBUTE I/O IN PARALLEL
  !
  !PRIVATE :: setio_channels       ! get/set I/O info of channels

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE initialize_parallel_io(status, p_pe, p_io, p_all_comm, p_nprocs)

    USE messy_main_channel_pnetcdf,     ONLY: ch_pnetcdf_init_pio
#ifdef HAVE_FORPY
    USE messy_main_channel_forpy,       ONLY: ch_forpy_init_pio
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: p_pe
    INTEGER, INTENT(IN)  :: p_io
    INTEGER, INTENT(IN)  :: p_all_comm
    INTEGER, INTENT(IN)  :: p_nprocs

    CALL ch_pnetcdf_init_pio(status, p_pe, p_io, p_all_comm)

#ifdef HAVE_FORPY
       CALL ch_forpy_init_pio(p_nprocs, p_pe)
#endif

  END SUBROUTINE initialize_parallel_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_init_restart(status, lp, lp_io, fname_base, rstatt)

    USE messy_main_channel_attributes, ONLY: t_attribute_list
    USE messy_main_channel_netcdf,     ONLY: ch_netcdf_init_rst
#ifdef HAVE_PNETCDF
    USE messy_main_channel_pnetcdf,    ONLY: ch_pnetcdf_init_rst
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    LOGICAL,                INTENT(OUT) :: lp
    LOGICAL,                INTENT(IN)  :: lp_io
    CHARACTER(LEN=*),       INTENT(IN)  :: fname_base
    TYPE(t_attribute_list), POINTER     :: rstatt     ! INTENT(IN)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'channel_init_restart'
    LOGICAL                     :: lex=.false.
    INTEGER                     :: i

    ! CHECK IF RESTART FILE IS PRESENT
    DO i=1,  FTYPE_MAXIMUM
       IF (lp_io .AND. (I_VERBOSE_LEVEL >= 1)) &
            WRITE(*,*) substr,': checking for file '//&
            &TRIM(fname_base)//TRIM(FTYPE_EXT_TEXT(i))//' ...'
       INQUIRE(file=TRIM(fname_base)//TRIM(FTYPE_EXT_TEXT(i)), exist=lex)
       IF (lex) EXIT
    END DO

    ! READ ATTRIBUTES FROM FILE
    ! ------------------------------------------------
    SELECT CASE(i)
    CASE(FTYPE_UNDEFINED)
       status = 3200 ! OUTPUT FILE TYPE UNDEFINED
    CASE(FTYPE_ASCII)
       status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
    CASE(FTYPE_NETCDF)
#else
    CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
       !
       lp = .FALSE.
       !
       IF (lp_io) THEN
          CALL ch_netcdf_init_rst(status, &
               TRIM(fname_base)//TRIM(FTYPE_EXT_TEXT(i)), rstatt)
       ELSE
          status = 0
       END IF
       !
#ifdef HAVE_PNETCDF
    CASE(FTYPE_PNETCDF)
       !
       lp = .TRUE.
       !
       CALL ch_pnetcdf_init_rst(status, &
            TRIM(fname_base)//TRIM(FTYPE_EXT_TEXT(i)), rstatt)
       !
#endif
    CASE(FTYPE_GRIB)
       status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
    CASE(FTYPE_HDF4)
       status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
    CASE(FTYPE_HDF5)
       status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_CDI
    CASE(FTYPE_CDI_NC)
       status = 0
#endif
#ifdef HAVE_FORPY
    CASE(FTYPE_FORPY)
       !
       ! FORPY cannot be used for restart files
       status = 5001
       !
#endif
    CASE DEFAULT
       status = 3206 ! RESTART FILE REQUIRED BUT NOT PRESENT
    END SELECT
    !
    IF (status /= 0) RETURN
    ! ------------------------------------------------

  END SUBROUTINE channel_init_restart
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_init_io(status, IOMODE, fname, AMODE, pe, att, chname &
       , dom_id, lclose, llp_io)

    USE messy_main_channel,         ONLY: new_attribute, AF_RST_CMP &
                                        , AF_RST_INP, get_attribute
    USE messy_main_channel_attributes, ONLY: t_attribute_list
    USE messy_main_channel_netcdf,  ONLY: ch_netcdf_init_io
#ifdef HAVE_PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_init_io
#endif

#ifdef HAVE_CDI
    USE messy_main_channel_cdi,     ONLY: ch_cdi_init_io
#endif
#ifdef HAVE_FORPY
    USE messy_main_channel_forpy,   ONLY: ch_forpy_init_io
#endif
    USE messy_main_tools,           ONLY: int2str
    USE messy_main_channel_mem,     ONLY: dom_current, dom_unbound

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,                INTENT(OUT)       :: status
    INTEGER,                INTENT(IN)        :: IOMODE
    CHARACTER(LEN=*),       INTENT(IN)        :: fname
    INTEGER,                INTENT(IN)        :: AMODE
    INTEGER,                INTENT(IN)       :: pe
    TYPE(t_attribute_list), POINTER, OPTIONAL :: att ! INTENT(IN)
    CHARACTER(LEN=*),       INTENT(IN), OPTIONAL :: chname
    INTEGER,                INTENT(IN), OPTIONAL :: dom_id
    LOGICAL,                INTENT(IN), OPTIONAL :: lclose
    LOGICAL,                INTENT(IN), OPTIONAL :: llp_io

    ! LOCAL
    !CHARACTER(LEN=*),      PARAMETER :: substr = 'channel_init_io'
    TYPE(t_channel_list),  POINTER   :: ls
    TYPE(t_channel),       POINTER   :: channel
    LOGICAL                          :: lskip
    LOGICAL                          :: lp
    LOGICAL                          :: lp_io
    INTEGER                          :: jg
    CHARACTER(LEN=2)                 :: domstr
    LOGICAL                          :: llclose

    IF (PRESENT(lclose)) THEN
       llclose = lclose
    ELSE
       llclose = .FALSE.
    END IF

    IF (PRESENT(dom_id)) THEN
       jg = dom_id
    ELSE
       jg = dom_current
    END IF

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       ! ------------------------------------------------
       SELECT CASE(AMODE)
       CASE(AMODE_READ)
          ! READ ALL AVAILABLE DATA
          lskip = .FALSE.
          IF (PRESENT(chname)) THEN
             IF ((TRIM(chname) /= TRIM(channel%name)) .OR. &
                  & (channel%dom_id /= jg)) lskip = .TRUE.
          ENDIF
       CASE(AMODE_WRITE)
          ! NEW OUTPUT FILE ? RESTART FILE ?
          SELECT CASE(IOMODE)
          CASE(IOMODE_OUT)
             lskip = .NOT. (channel%int%lnew_file .AND. channel%int%lout_now) &
                  .OR. (channel%dom_id /= jg)
          CASE(IOMODE_RST)
             lskip = .NOT. channel%int%lrst_now
          END SELECT
       END SELECT
       IF (lskip) THEN
          ls => ls%next
          CYCLE
       END IF
       ! ------------------------------------------------

       ! ------------------------------------------------
       IF (channel%io%ftype(IOMODE) == FTYPE_UNDEFINED) THEN
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
          RETURN
       END IF

       IF (channel%io%ftype(IOMODE) > FTYPE_MAXIMUM) THEN
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
          RETURN
       END IF

       ! => domain is added to output file name, if not on unbound domain
       IF (channel%dom_id > dom_unbound) THEN
          CALL int2str(domstr, channel%dom_id)
          channel%int%fname(IOMODE)= &
               TRIM(fname)//TRIM(channel%name)//"_D"//TRIM(domstr)//&
               &FTYPE_EXT_TEXT(channel%io%ftype(IOMODE))
       ELSE
          channel%int%fname(IOMODE)= &
               TRIM(fname)//TRIM(channel%name)//&
               &FTYPE_EXT_TEXT(channel%io%ftype(IOMODE))
       END IF
       ! ------------------------------------------------

       ! ------------------------------------------------
       CALL new_attribute(status, channel%att &
            , 'channel_name', c=TRIM(channel%name)&
            , loverwrite=.TRUE., iflag = AF_RST_CMP, lout=.TRUE.)
       IF (status /= 0) RETURN

       CALL new_attribute(status, channel%att  &
            , 'channel_file_type', c=TRIM(IOMODE_TEXT(IOMODE)) &
            , loverwrite=.TRUE., iflag = AF_RST_CMP, lout=.TRUE.)
       IF (status /= 0) RETURN

       CALL new_attribute(status, channel%att &
            , 'channel_file_name', c=TRIM(channel%int%fname(IOMODE)) &
            , loverwrite=.TRUE., lout=.TRUE.)
       IF (status /= 0) RETURN

       IF (IOMODE == IOMODE_RST) THEN
          !
          CALL new_attribute(status, channel%att &
               , 'channel_time_slo', r=channel%int%tslo &
               , loverwrite=.TRUE., iflag = AF_RST_INP, lout=.TRUE.)
          IF (status /= 0) RETURN
       END IF

       ! ------------------------------------------------

       ! ------------------------------------------------
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          IF (.NOT. PRESENT(llp_io)) THEN
            ! IF IO_PE DEPENDING ON TASK
            lp_io = (pe == channel%int%netcdf(IOMODE)%io_pe)
          ELSE
            lp_io = llp_io
          ENDIF
          !
          lp = .FALSE.
          !
          IF (lp_io) THEN
             CALL ch_netcdf_init_io(status, IOMODE, channel, AMODE, att &
                  , lclose)
          ELSE
             status = 0
          END IF
          !
#ifdef HAVE_PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          lp = .TRUE.
          !
          CALL ch_pnetcdf_init_io(status, IOMODE, channel, AMODE, att &
               , lclose)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_CDI
       CASE(FTYPE_CDI_NC)
          lp = .FALSE.
          IF (.NOT. PRESENT(llp_io)) THEN
            ! IF IO_PE DEPENDING ON TASK
            lp_io = (pe == channel%int%cdi(IOMODE)%io_pe)
          ELSE
            lp_io = llp_io
          ENDIF
          IF (lp_io) THEN
             CALL ch_cdi_init_io(status, IOMODE, channel, AMODE, att &
                  lclose)
          ELSE
             status = 0
          END IF
#endif
#ifdef HAVE_FORPY
       CASE(FTYPE_FORPY)
          lp = .FALSE.
          CALL ch_forpy_init_io(status, IOMODE, channel, AMODE)
#endif
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       ! ------------------------------------------------

       ! ------------------------------------------------
       IF ( (lp_io .OR. lp) .AND. .NOT. llclose) THEN
          IF (AMODE == AMODE_READ) THEN
             CALL get_attribute(status, channel%att, 'channel_time_slo' &
                  , r=channel%int%tslo)
             IF (status /= 0) RETURN
          END IF
       END IF
       ! ------------------------------------------------

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE channel_init_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_write_header(status, IOMODE, DIMID_TIME, pe, att)

    USE messy_main_channel_dimensions, ONLY: get_dimension, t_dimension
    USE messy_main_channel_attributes, ONLY: t_attribute_list
    USE messy_main_channel_netcdf,     ONLY: ch_netcdf_write_header
#ifdef HAVE_PNETCDF
    USE messy_main_channel_pnetcdf,    ONLY: ch_pnetcdf_write_header
#endif
#ifdef HAVE_CDI
    USE messy_main_channel_cdi,        ONLY: ch_cdi_write_header
#endif
#ifdef HAVE_FORPY
    USE messy_main_channel_forpy,      ONLY: ch_forpy_write_header
#endif
    USE messy_main_channel_mem,        ONLY: dom_current

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT)       :: status
    INTEGER,                INTENT(IN)        :: IOMODE
    INTEGER,                INTENT(IN)        :: DIMID_TIME
    INTEGER,                INTENT(IN)        :: pe
    TYPE(t_attribute_list), POINTER, OPTIONAL :: att ! INTENT(IN)

    ! LOCAL
    !CHARACTER(LEN=*),    PARAMETER :: substr = 'channel_write_header'
    TYPE(t_dimension),     POINTER :: dim_time => NULL()
    TYPE(t_channel_list),  POINTER :: ls       => NULL()
    TYPE(t_channel),       POINTER :: channel  => NULL()
    LOGICAL                        :: lskip
    LOGICAL                        :: lp_io

    ! SET TIME DIMENSION
    CALL get_dimension(status, DIMID_TIME, dim_time)
    IF (status /= 0) RETURN

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       ! ------------------------------------------------
       ! NEW OUTPUT FILE ? RESTART FILE ?
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          lskip = .NOT. (channel%int%lnew_file .AND. channel%int%lout_now)
          lskip = lskip .OR. (channel%dom_id /= dom_current)
       CASE(IOMODE_RST)
          lskip = .NOT. channel%int%lrst_now
       END SELECT
       IF (lskip) THEN
          ls => ls%next
          CYCLE
       END IF
       ! ------------------------------------------------

       ! ------------------------------------------------
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          ! IF IO_PE DEPENDING ON TASK
          lp_io = (pe == channel%int%netcdf(IOMODE)%io_pe)
          !
          IF (lp_io) THEN
             CALL ch_netcdf_write_header(status, IOMODE, channel &
                  , dim_time, att)
          ELSE
             status = 0
          END IF
          !
#ifdef HAVE_PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          CALL ch_pnetcdf_write_header(status, IOMODE, channel &
               , dim_time, att)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_CDI
       CASE(FTYPE_CDI_NC)
          lp_io = (pe == channel%int%cdi(IOMODE)%io_pe)
          IF (lp_io) THEN
             CALL ch_cdi_write_header(status, IOMODE, channel, att)
          ELSE
             status = 0
          END IF
#endif
#ifdef HAVE_FORPY
       CASE(FTYPE_FORPY)
          CALL ch_forpy_write_header(status, IOMODE, channel)
#endif
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       ! ------------------------------------------------

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE channel_write_header
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_write_time(status, IOMODE, DIMID_TIME, pe)

    USE messy_main_channel_dimensions, ONLY: get_dimension, t_dimension
    USE messy_main_channel_netcdf,     ONLY: ch_netcdf_write_time
#ifdef HAVE_PNETCDF
    USE messy_main_channel_pnetcdf,    ONLY: ch_pnetcdf_write_time
#endif
#ifdef HAVE_CDI
    USE messy_main_channel_cdi,        ONLY: ch_cdi_write_time
#endif
#ifdef HAVE_FORPY
    USE messy_main_channel_forpy,      ONLY: ch_forpy_write_time
#endif
    USE messy_main_channel_mem,        ONLY: dom_current

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: IOMODE
    INTEGER, INTENT(IN)  :: DIMID_TIME
    INTEGER, INTENT(IN)  :: pe

    ! LOCAL
    !CHARACTER(LEN=*),   PARAMETER :: substr = 'channel_write_time'
    TYPE(t_dimension),    POINTER :: dim_time => NULL()
    TYPE(t_channel_list), POINTER :: ls       => NULL()
    TYPE(t_channel),      POINTER :: channel  => NULL()
    LOGICAL                       :: lp_io

    IF (IOMODE /= IOMODE_OUT) RETURN

    ! SET TIME DIMENSION
    CALL get_dimension(status, DIMID_TIME, dim_time)
    IF (status /= 0) RETURN

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       ! ------------------------------------------------
       ! OUTPUT NOW ?
       IF ( (.NOT. channel%int%lout_now) .OR. &
            (channel%dom_id /= dom_current) ) THEN
          ls => ls%next
          CYCLE
       END IF
       ! ------------------------------------------------

       ! ------------------------------------------------
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          ! IF IO_PE DEPENDING ON TASK
          lp_io = (pe == channel%int%netcdf(IOMODE)%io_pe)
          !
          IF (lp_io) THEN
             CALL ch_netcdf_write_time(status, channel, dim_time)
          ELSE
             status = 0
          END IF
          !
#ifdef HAVE_PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          CALL ch_pnetcdf_write_time(status, channel, dim_time)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_CDI
       CASE(FTYPE_CDI_NC)
          ! IF IO_PE DEPENDING ON TASK
          lp_io = (pe == channel%int%cdi(IOMODE)%io_pe)
          IF (lp_io) THEN
             CALL ch_cdi_write_time(status, channel)
          ELSE
             status = 0
          END IF
#endif
#ifdef HAVE_FORPY
       CASE(FTYPE_FORPY)
          CALL ch_forpy_write_time(status, channel, dim_time)
#endif
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       ! ------------------------------------------------

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE channel_write_time
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_write_data(status, lp, IOMODE, lexit, ptr, reprid &
               , p_io_c, pe)

    USE messy_main_channel_repr,   ONLY: REPR_UNDEF, repr_reorder
    USE messy_main_channel_netcdf, ONLY: ch_netcdf_write_data
#ifdef HAVE_PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_write_data
#endif
#ifdef HAVE_CDI
    USE messy_main_channel_cdi,     ONLY: ch_cdi_write_data
#endif
#ifdef HAVE_FORPY
    USE messy_main_channel_forpy,   ONLY: ch_forpy_write_data
#endif
    USE messy_main_channel_mem,     ONLY: dom_current

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                      INTENT(OUT)   :: status
    LOGICAL,                      INTENT(OUT)   :: lp
    INTEGER,                      INTENT(IN)    :: IOMODE
    LOGICAL,                      INTENT(OUT)   :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER       :: ptr
    INTEGER,                      INTENT(OUT)   :: reprid
    ! check for IO PE of each CHANNEL and return to bi_decompose
    INTEGER,                      INTENT(OUT)   :: p_io_c
    INTEGER,                      INTENT(IN)    :: pe

    ! LOCAL
    INTEGER, PARAMETER :: MODE_INITIALIZE   = 0
    INTEGER, PARAMETER :: MODE_NEXT_CHANNEL = 1
    INTEGER, PARAMETER :: MODE_NEXT_OBJECT  = 2
    INTEGER, PARAMETER :: MODE_NEXT_DATA    = 3
    INTEGER, PARAMETER :: MODE_OUTPUT       = 4

    !CHARACTER(LEN=*),     PARAMETER     :: substr = 'channel_write_data'
    INTEGER,                       SAVE :: MODE = MODE_INITIALIZE
    TYPE(t_channel_list), POINTER, SAVE :: ls
    TYPE(t_channel),      POINTER, SAVE :: channel
    TYPE(t_channel_object_list), POINTER, SAVE :: le
    TYPE(t_channel_object),      POINTER, SAVE :: object
    ! OUTPUT DATA TYPE
    INTEGER,                              SAVE :: jsnd = 0
    ! INDEX IN SECONDARY DATA POINTER
    INTEGER,                              SAVE :: i2nd = 0
    REAL(DP), DIMENSION(:,:,:,:), POINTER      :: zptr => NULL()
    LOGICAL                                    :: loutstatic
    LOGICAL                                    :: lp_io

    ! INIT
    lexit = .FALSE.

    DO

    SELECT CASE(mode)

    CASE(MODE_INITIALIZE)
       !
       ! INIT
       NULLIFY(ptr)
       !
       ls => GCHANNELLIST
       !
       reprid = REPR_UNDEF
       NULLIFY(le)
       NULLIFY(channel)
       NULLIFY(object)
       jsnd = 0
       i2nd = 0
       !
       MODE = MODE_NEXT_CHANNEL
       !
    CASE(MODE_NEXT_CHANNEL)
       !
       ! LOOK FOR NEXT CHANNEL WITH OUTPUT
       DO
          IF (.NOT. ASSOCIATED(ls)) THEN
             lexit = .TRUE.                  ! NO MORE CHANNEL
             lp    = .FALSE.                 !!$! ... always broadcast lexit ...
             MODE = MODE_INITIALIZE          ! NEXT OUTPUT
             p_io_c = -1
             status = 0
             RETURN
          ELSE
             channel => ls%this
             ! OUTPUT OR RESTART
             IF ( (IOMODE == IOMODE_OUT .AND. channel%int%lout_now &
                  .AND. (channel%dom_id == dom_current)) .OR. &
                  (IOMODE == IOMODE_RST .AND. channel%int%lrst_now) ) THEN
                ! O.K. CHANNEL FOR OUTPUT
                le => channel%list
                MODE = MODE_NEXT_OBJECT
                EXIT
             END IF
          END IF
          ls => ls%next
       END DO
       !
    CASE(MODE_NEXT_OBJECT)
       !
       ! LOOK FOR NEXT OBJECT WITH OUTPUT
       DO
          IF (.NOT. ASSOCIATED(le)) THEN  ! NO MORE OBJECT
             ls => ls%next
             MODE = MODE_NEXT_CHANNEL
             EXIT
          ELSE
             object => le%this
             loutstatic = (.NOT. object%lstatic) .OR. &
                          (       object%lstatic .AND. channel%int%lnew_file)
             IF ( (((IOMODE == IOMODE_OUT) .AND. object%int%lout) &
                  .AND. loutstatic ) &
                  .OR. &
                  ((IOMODE == IOMODE_RST) .AND. object%int%lrst) ) THEN
                ! O.K. OBJECT FOR OUTPUT
                jsnd = 1
                MODE = MODE_NEXT_DATA
                EXIT
             END IF
          END IF
          le => le%next
       END DO
       !
    CASE(MODE_NEXT_DATA)
       !
       ! LOOK FOR NEXT DATA WITH OUTPUT
       DO
          IF (jsnd > SND_MAXLEN) THEN ! NO MORE DATA AVAILABLE
             le => le%next
             MODE = MODE_NEXT_OBJECT
             EXIT
          ELSE
             !
             SELECT CASE(jsnd)
             CASE(SND_INS)
                ptr => object%ioptr(:,:,:,:)
             CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX, ...
                i2nd = object%int%i2nd(jsnd)
                IF (i2nd > 0) THEN
                   ptr => object%sdat(i2nd)%ptr(:,:,:,:)
                ELSE
                   NULLIFY(ptr)
                END IF
             END SELECT
             !
             IF (object%int%lexp(jsnd, IOMODE)) THEN
                ! OUTPUT OF PRIMARY/SECONDARY DATA REQUESTED
                reprid = object%repr%id
                MODE = MODE_OUTPUT
                IF (I_VERBOSE_LEVEL >= 5) write (*,*) 'OUT OBJECT ',object%name
                !
                SELECT CASE(channel%io%ftype(IOMODE))
                CASE(FTYPE_UNDEFINED)
                   status = 3200 ! OUTPUT FILE TYPE UNDEFINED
                CASE(FTYPE_ASCII)
                   status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
                CASE(FTYPE_NETCDF)
#else
                CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
                   !
                   p_io_c = channel%int%netcdf(IOMODE)%io_pe
                   lp = .FALSE.
                   status = 0
                   !
#ifdef HAVE_PNETCDF
                CASE(FTYPE_PNETCDF)
                   !
                   lp = .TRUE.
                   status = 0
                   !
#endif
                CASE(FTYPE_GRIB)
                   status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
                CASE(FTYPE_HDF4)
                   status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
                CASE(FTYPE_HDF5)
                   status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_CDI
                CASE(FTYPE_CDI_NC)
                   p_io_c = channel%int%cdi(IOMODE)%io_pe
                   lp=.FALSE.
                   status=0
#endif
#ifdef HAVE_FORPY
                CASE(FTYPE_FORPY)
                   p_io_c = channel%int%forpy%io_pe
                   lp=.FALSE.
                   status=0
#endif
                CASE DEFAULT
                   status = 3202 ! OUTPUT FILE TYPE UNKNOWN
                END SELECT
                !
                RETURN ! JUMP BACK TO GATHER ON ONE PE
             END IF
          END IF
          jsnd = jsnd + 1
       END DO
       !
    CASE(MODE_OUTPUT)
       !
!!$       IF (ASSOCIATED(ptr)) THEN  ! I/O - PE
!!$       END IF
       !
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          ! IF IO_PE DEPENDING ON TASK
          lp_io = (pe == channel%int%netcdf(IOMODE)%io_pe)
          p_io_c = channel%int%netcdf(IOMODE)%io_pe
          lp = .FALSE.
          !
          IF (lp_io) THEN
             ! CHECK, IF OBJECT HAS A "gather-able" SHAPE
             IF (ASSOCIATED(ptr)) THEN
                CALL repr_reorder(status, 1, lp, object%repr, ptr, zptr)
                IF (status /= 0) RETURN
                !
                CALL ch_netcdf_write_data(status, IOMODE, channel &
                     , object, zptr, jsnd, i2nd)
             ELSE
                status=0
             ENDIF
          ELSE
             status = 0
          END IF
          !
#ifdef HAVE_PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          lp = .TRUE.
          !
          ! OBJECT OUTPUT IS NOT IMPLEMENTED YET
          IF (ASSOCIATED(ptr)) THEN
             CALL repr_reorder(status, 1, lp, object%repr, ptr, zptr)
             IF (status /= 0) RETURN
             !
             CALL ch_pnetcdf_write_data(status, IOMODE, channel &
                  , object, zptr, jsnd, i2nd)
          ELSE
             status=0
          ENDIF
             !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_CDI
       CASE(FTYPE_CDI_NC)
          !
          ! IF IO_PE DEPENDING ON TASK
          lp_io = (pe == channel%int%cdi(IOMODE)%io_pe)
          p_io_c = channel%int%cdi(IOMODE)%io_pe
          lp = .FALSE.
          !
          IF (lp_io) THEN
             ! CHECK, IF OBJECT HAS A "gather-able" SHAPE
             IF (object%int%cdi(IOMODE)%lout) THEN
                IF (ASSOCIATED(ptr)) THEN
                   CALL ch_cdi_write_data(status, IOMODE, channel, object &
                        &               , ptr, jsnd, i2nd)
                END IF
             END IF
             status=0
          ELSE
             status = 0
          END IF
#endif
#ifdef HAVE_FORPY
       CASE(FTYPE_FORPY)
          !
          p_io_c = channel%int%forpy%io_pe
          lp = .FALSE.
          !
          ! OBJECT OUTPUT IS NOT IMPLEMENTED YET
          IF (ASSOCIATED(ptr)) THEN
             CALL repr_reorder(status, 1, lp, object%repr, ptr, zptr)
             IF (status /= 0) RETURN
             !
             CALL ch_forpy_write_data(status, IOMODE, channel &
                  , object, zptr, jsnd)
          ELSE
             status=0
          ENDIF
          !
#endif
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       ! CLEAN UP
       IF (ASSOCIATED(zptr)) THEN
          DEALLOCATE(zptr)
          NULLIFY(zptr)
       END IF
       !
       jsnd = jsnd + 1
       MODE = MODE_NEXT_DATA
       RETURN
       !
    END SELECT

    END DO

  END SUBROUTINE channel_write_data
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_finish_io(status,  IOMODE, lclose, pe, chname, dom_id &
           , llp_io)

    USE messy_main_channel_netcdf, ONLY: ch_netcdf_finish_io
#ifdef HAVE_PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_finish_io
#endif
!!$#ifdef HAVE_CDI
    USE messy_main_channel_cdi, ONLY: ch_cdi_finish_io
!!$#endif
#ifdef HAVE_FORPY
    USE messy_main_channel_forpy, ONLY: ch_forpy_finish_io
#endif
    USE messy_main_channel_mem, ONLY: dom_current

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    INTEGER,  INTENT(IN)  :: IOMODE
    LOGICAL,  INTENT(IN)  :: lclose
    INTEGER,  INTENT(IN)  :: pe
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: chname
    INTEGER,          INTENT(IN), OPTIONAL :: dom_id
    LOGICAL,          INTENT(IN), OPTIONAL :: llp_io

    ! LOCAL
    !CHARACTER(LEN=*),     PARAMETER :: substr = 'channel_finish_io'
    TYPE(t_channel_list), POINTER   :: ls
    TYPE(t_channel),      POINTER   :: channel
    LOGICAL                         :: lskip
    LOGICAL                         :: lcycle
    INTEGER                         :: jg
    LOGICAL                         :: lp_io

    IF (PRESENT(dom_id)) THEN
       jg = dom_id
    ELSE
       jg = dom_current
    END IF

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       lcycle = .FALSE.
       ! ------------------------------------------------
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          lskip = .NOT. channel%int%lout_now
          lskip = lskip .OR. (channel%dom_id /= jg)
       CASE(IOMODE_RST)
          lskip = .NOT. channel%int%lrst_now
          IF (PRESENT(chname)) THEN
             IF ((TRIM(chname) /= TRIM(channel%name)) .OR. &
                  & (channel%dom_id /= jg)) lcycle = .TRUE.
          ENDIF
       END SELECT
       IF ((lskip .AND. (.NOT.lclose)) .OR. lcycle) THEN
          ls => ls%next
          CYCLE
       END IF
       ! ------------------------------------------------

       ! RESET TRIGGER FOR NEW FILE
       IF (IOMODE == IOMODE_OUT) channel%int%lnew_file = .FALSE.

       ! ------------------------------------------------
       ! OPEN FILE ? (= NAME ASSIGNED ?)
       !
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          IF (.NOT. PRESENT(llp_io)) THEN
            ! IF IO_PE DEPENDING ON TASK
            lp_io = (pe == channel%int%netcdf(IOMODE)%io_pe)
          ELSE
            lp_io = llp_io
          ENDIF
          !
          IF (lp_io) THEN
             CALL ch_netcdf_finish_io(status, IOMODE, channel, lclose)
          ELSE
             status = 0
          END IF
          !
#ifdef HAVE_PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          CALL ch_pnetcdf_finish_io(status, IOMODE, channel, lclose)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
!!$#ifdef HAVE_CDI
       CASE(FTYPE_CDI_NC)
          IF (.NOT. PRESENT(llp_io)) THEN
            ! IF IO_PE DEPENDING ON TASK
            lp_io = (pe == channel%int%cdi(IOMODE)%io_pe)
          ELSE
            lp_io = llp_io
          ENDIF
          IF (lp_io) THEN
             CALL ch_cdi_finish_io(status, IOMODE, channel, lclose)
          ELSE
             status = 0
          END IF
!!$#endif
#ifdef HAVE_FORPY
       CASE(FTYPE_FORPY)
          CALL ch_forpy_finish_io(status, IOMODE, channel, lclose)
#endif
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       !
       ! ------------------------------------------------

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE channel_finish_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_read_data(status, lp_io, IOMODE, lexit, ptr, reprid, lp &
       , chname, dom_id &
       , lskipinput, filename, objname, chaname, lignore, lrestreq)

    USE messy_main_channel_repr,   ONLY: REPR_UNDEF, repr_reorder
    USE messy_main_channel_netcdf, ONLY: ch_netcdf_read_data
#ifdef HAVE_PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_read_data
#endif
    USE messy_main_channel_mem,     ONLY: dom_current

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                      INTENT(OUT)   :: status
    LOGICAL,                      INTENT(IN)    :: lp_io
    INTEGER,                      INTENT(IN)    :: IOMODE
    LOGICAL,                      INTENT(OUT)   :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER       :: ptr
    INTEGER,                      INTENT(OUT)   :: reprid
    LOGICAL,                      INTENT(OUT)   :: lp
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL      :: chname
    INTEGER,          INTENT(IN), OPTIONAL      :: dom_id
    LOGICAL,          INTENT(IN), OPTIONAL      :: lskipinput
    CHARACTER(LEN=34+STRLEN_CHANNEL+8), INTENT(OUT), OPTIONAL :: filename
    CHARACTER(LEN=STRLEN_OBJECT+4),     INTENT(OUT), OPTIONAL :: objname
    CHARACTER(LEN=STRLEN_CHANNEL),      INTENT(OUT), OPTIONAL :: chaname
    LOGICAL,          INTENT(OUT), OPTIONAL      :: lrestreq
    LOGICAL,          INTENT(OUT), OPTIONAL      :: lignore

    ! LOCAL
    INTEGER, PARAMETER :: MODE_INITIALIZE   = 0
    INTEGER, PARAMETER :: MODE_NEXT_CHANNEL = 1
    INTEGER, PARAMETER :: MODE_NEXT_OBJECT  = 2
    INTEGER, PARAMETER :: MODE_NEXT_DATA    = 3
    INTEGER, PARAMETER :: MODE_INPUT        = 4
    INTEGER, PARAMETER :: MODE_DISTRIBUTE   = 5

    !CHARACTER(LEN=*),     PARAMETER     :: substr = 'channel_read_data'
    INTEGER,                       SAVE :: MODE = MODE_INITIALIZE
    TYPE(t_channel_list), POINTER, SAVE :: ls
    TYPE(t_channel),      POINTER, SAVE :: channel
    TYPE(t_channel_object_list), POINTER, SAVE :: le
    TYPE(t_channel_object),      POINTER, SAVE :: object
    ! OUTPT DATA TYPE
    INTEGER,                              SAVE :: jsnd = 0
    ! INDEX IN SECONDARY DATA POINTER
    INTEGER,                              SAVE :: i2nd = 0
    REAL(DP), DIMENSION(:,:,:,:), POINTER      :: zptr => NULL()
    INTEGER                                    :: jg
    LOGICAL                                    :: lskpinp = .FALSE.

    ! ONLY FOR RESTART FILES
    IF (IOMODE /= IOMODE_RST) THEN
       status = 3205 ! NO INPUT OF OUTPUT FILES
       RETURN
    END IF

    ! INIT
    lexit = .FALSE.

    IF (PRESENT(dom_id)) THEN
       jg = dom_id
    ELSE
       jg = dom_current
    END IF

    IF (PRESENT(lskipinput)) THEN
       lskpinp = lskipinput
    ELSE
       lskpinp = .FALSE.
    END IF

    DO

    SELECT CASE(mode)

    CASE(MODE_INITIALIZE)
       !
       ! INIT
       IF (ASSOCIATED(ptr)) DEALLOCATE(ptr)
       NULLIFY(ptr)
       !
       ls => GCHANNELLIST
       !
       reprid = REPR_UNDEF
       NULLIFY(le)
       NULLIFY(channel)
       NULLIFY(object)
       jsnd = 0
       i2nd = 0
       !
       MODE = MODE_NEXT_CHANNEL
       !
    CASE(MODE_NEXT_CHANNEL)
       !
       ! LOOK FOR NEXT CHANNEL
       IF (.NOT. ASSOCIATED(ls)) THEN
          lexit = .TRUE.                  ! NO MORE CHANNEL
          lp    = .FALSE.                 !!$! ... always broadcast lexit ...
          MODE = MODE_INITIALIZE          ! NEXT INPUT
          status = 0
          RETURN
       ELSE
          channel => ls%this
          le => channel%list
          MODE = MODE_NEXT_OBJECT
          IF (PRESENT(chname)) THEN
             IF ((TRIM(chname) /= TRIM(channel%name)) .OR. &
                  & (channel%dom_id /= jg))THEN
                ls => ls%next
                MODE = MODE_NEXT_CHANNEL
             ENDIF
          ENDIF
       END IF
       !
    CASE(MODE_NEXT_OBJECT)
       !
       ! LOOK FOR NEXT OBJECT
       IF (.NOT. ASSOCIATED(le)) THEN  ! NO MORE OBJECT
          ls => ls%next
          MODE = MODE_NEXT_CHANNEL
       ELSE
          object => le%this
          jsnd = 1
          MODE = MODE_NEXT_DATA
       END IF
       !
    CASE(MODE_NEXT_DATA)
       !
       ! LOOK FOR NEXT DATA
       DO
          IF (jsnd > SND_MAXLEN) THEN ! NO MORE DATA AVAILABLE
             le => le%next
             MODE = MODE_NEXT_OBJECT
             EXIT
          ELSE
             IF (object%int%lrestart_read(jsnd)) THEN
                ! DATA was alread read
                jsnd = jsnd + 1
                MODE = MODE_NEXT_DATA
                EXIT
             END IF

             reprid = object%repr%id
             MODE = MODE_INPUT
             !
             SELECT CASE(jsnd)
             CASE(SND_INS)
                ! PRIMARY DATA
                EXIT
             CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX
                ! SEONDARY DATA, IF PRESENT
                i2nd = object%int%i2nd(jsnd)
                IF (i2nd > 0) THEN
                   EXIT
                END IF
             END SELECT
             !
          END IF
          jsnd = jsnd + 1
       END DO
       !
    CASE(MODE_INPUT)
       !
       IF (lskpinp) THEN
          status = 0
          ! RETURN TO READ and SCATTER ON ALL PEs AND DISTRIBUTE AFTERWARDS
          MODE      = MODE_DISTRIBUTE
          !write (*,*) 'SET FILENAME ', PRESENT(filename) &
          ! , TRIM(channel%int%fname(IOMODE)), object%lrestreq, object%int%lign
          IF (PRESENT(filename)) THEN
             filename  = TRIM(channel%int%fname(IOMODE))
          END IF
          IF (PRESENT(objname)) THEN
             objname   = TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))
          END IF
          IF (PRESENT(chaname)) THEN
             chaname   = TRIM(channel%name)
          END IF
          IF (PRESENT(lrestreq)) lrestreq = object%lrestreq
          ! LFIXATE is needed here, as messy_read_restart needs
          ! to be called twice in ICON. The first time BEFORE (!)
          ! messy_init_coupling and therefore before fixate_channels, i.e.
          ! in the first call it is still unknown, if lignore for an object
          ! in the restart files is .TRUE. or .FALSE.
          ! => LOGIC: assume lignore in the first call, interrupt in the second
          !   CALL if lignore=F
          IF (PRESENT(lignore)) lignore = object%int%lign .OR. .NOT. LFIXATE
          RETURN
       END IF

       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          lp = .FALSE.
          !
          IF (lp_io) THEN  ! I/O - PE
             CALL ch_netcdf_read_data(status, IOMODE, channel &
                  , object, zptr, jsnd, i2nd)
          ELSE
             ! NOTHING TO DO FOR NON-I/O PE
             IF (ASSOCIATED(ptr)) DEALLOCATE(ptr)
             NULLIFY(ptr)
             status = 0
          END IF
          !
#ifdef HAVE_PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          lp = .TRUE.
          !
          CALL ch_pnetcdf_read_data(status, IOMODE, channel &
               , object, zptr, jsnd, i2nd)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_CDI
       CASE(FTYPE_CDI_NC)
          status = 0
#endif
#ifdef HAVE_FORPY
       CASE(FTYPE_FORPY)
          status = 5002
#endif
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       !
       ! REORDER DATA
       IF (ASSOCIATED(zptr)) THEN
          CALL repr_reorder(status, -1, lp, object%repr, ptr, zptr)
          IF (status /= 0) RETURN
          ! CLEAN UP
          DEALLOCATE(zptr)
          NULLIFY(zptr)
       END IF
       !
       ! RETURN TO SCATTER ON ALL PEs AND DISTRIBUTE AFTERWARDS
       MODE = MODE_DISTRIBUTE
       RETURN
       !
    CASE(MODE_DISTRIBUTE)
       !
       IF (ASSOCIATED(ptr)) THEN
          ! DISTRIBUTE ON ALL PEs AFTER SCATTER
          SELECT CASE(jsnd)
          CASE(SND_INS)
             object%ioptr(:,:,:,:) = ptr(:,:,:,:)
          CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX, ...
             object%sdat(i2nd)%ptr(:,:,:,:) = ptr(:,:,:,:)
          END SELECT
          !
          ! FLAG DATA AS READ
          object%int%lrestart_read(jsnd) = .TRUE.
          !
       ELSE
          object%int%lrestart_read(jsnd) = .FALSE.
       END IF
       !
       ! NEXT DATA
       jsnd = jsnd + 1
       MODE = MODE_NEXT_DATA
       RETURN
       !
    END SELECT

    END DO

  END SUBROUTINE channel_read_data
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_dist_io(status, nodeid, p_io, ppiotype)

    IMPLICIT NONE

    ! I/O
    INTEGER,                           INTENT(OUT) :: status
    INTEGER, DIMENSION(:),             INTENT(IN)  :: nodeid
    INTEGER,                           INTENT(IN)  :: p_io
    INTEGER, DIMENSION(FTYPE_MAXIMUM), INTENT(IN)  :: ppiotype
    ! LOCAL
    LOGICAL, DIMENSION(:), POINTER :: lout
#ifndef HAVE_CDI .AND. HAVE_FORPY
    INTEGER, PARAMETER, DIMENSION(1) :: ft = (/FTYPE_NETCDF/)
#else
    INTEGER, PARAMETER, DIMENSION(1) :: ft = (/FTYPE_NETCDF &
#ifdef HAVE_FORPY
                ,FTYPE_FORPY &
#endif
#ifdef HAVE_CDI
                ,FTYPE_CDI_NC &
#endif
    /)
#endif
    INTEGER :: jft
    INTEGER, DIMENSION(:), ALLOCATABLE :: pe

    status = 0

    NULLIFY(lout)

    DO jft=1, SIZE(ft)
       ! retrieve list of output/restart requests for specific filetype
       CALL setio_channels(status, ft(jft), get=lout)
       IF (status /=0) RETURN
       !
       ! distribute according to selected algorithm
       CALL distribute(status, pe, lout, nodeid, ppiotype(ft(jft)), p_io)
       IF (status /= 0) RETURN
       DEALLOCATE(lout); NULLIFY(lout)
       !
       CALL setio_channels(status, ft(jft), set=pe)
       IF (status /=0) RETURN
       DEALLOCATE(pe)
       !
    END DO

  CONTAINS

    SUBROUTINE distribute(stat, zpe, zlout, znodeid, piotype, zp_io)

      IMPLICIT NONE
      INTRINSIC :: MAXVAL

      ! I/O
      INTEGER,                            INTENT(OUT) :: stat
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: zpe
      LOGICAL, DIMENSION(:),              INTENT(IN)  :: zlout
      INTEGER, DIMENSION(:),              INTENT(IN)  :: znodeid
      INTEGER,                            INTENT(IN)  :: piotype
      INTEGER,                            INTENT(IN)  :: zp_io

      ! LOCAL
      INTEGER :: nt   ! number of tasks
      INTEGER :: nn   ! number of nodes
      INTEGER :: ntpn ! number of tasks per node
      INTEGER :: i    ! loop counter for channel (output)
      INTEGER :: pe
      INTEGER :: oc   ! offset counter

      stat = 0
      ALLOCATE(zpe(SIZE(zlout)))

      SELECT CASE(piotype)
      CASE(0)
         ! OLD DEFAULT: all on the standard I/O PE
         zpe(:) = zp_io
         !
      CASE(1)
         !
         ! PRESET WITH DEFAULT
         zpe(:) = zp_io
         !
         ! round robin of first task on each node
         !
         ! number of tasks and nodes
         nt = SIZE(znodeid)
         nn = MAXVAL(znodeid) + 1 ! node counting starts at 0 !
         ! tasks per node
         ntpn = nt / nn
         !
         ! first task ID (PE)
         pe = 0
         !
         oc = 0 ! offset counter
         !
         DO i=1, SIZE(zlout)
            IF (zlout(i)) THEN
               zpe(i) = pe
               IF (zpe(i) > (nt-1)) THEN
                  status = 3214 ! out of range, must never be reached!
                  RETURN
               ENDIF
               pe = pe + ntpn
               ! reset, if out of range
               IF (pe > nt-1) THEN       ! task counting starts at 0 !
                  oc = oc + 1
                  IF (oc > (ntpn-1)) THEN  ! task counting starts at 0 !
                     oc = 0
                  END IF
                  pe = oc
               END IF
            END IF
         END DO
         !
      CASE DEFAULT
         !
         stat = 3213
         !
      END SELECT

    END SUBROUTINE distribute

  END SUBROUTINE channel_dist_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE setio_channels(status, FTYPE, get, set)

    ! THIS ROUTINE IS TO
    ! - get: retrieve a LOGICAL list of channels with output via a
    !        selected (FTYPE) I/O backend (FTYPE)
    ! - set: assign a list of PEs for semi-parallel I/O
    !        (for those backends for which this is possible)

    USE messy_main_channel, ONLY: new_attribute

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: FTYPE  ! I/O backend (file type)
    LOGICAL, DIMENSION(:), POINTER,    OPTIONAL :: get !  OUT
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: set

    ! LOCAL
    TYPE(t_channel_list), POINTER :: ls
    TYPE(t_channel),      POINTER :: channel

    status = 0

    IF (PRESENT(get)) THEN
       IF (ASSOCIATED(get)) THEN
          DEALLOCATE(get); NULLIFY(get)
       END IF
       ALLOCATE(get(NCHANNEL))
       get(:) = .FALSE.
    END IF

    IF (PRESENT(set)) THEN
       IF (SIZE(set) /= NCHANNEL) THEN
          status = 3005 ! NUMBER OF CHANNELS MISMATCH
          RETURN
       END IF
    END IF

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT
       channel => ls%this

       IF (PRESENT(get)) THEN
          get(channel%id) = &
               (channel%io%ftype(IOMODE_OUT) == FTYPE .AND. channel%int%lout) &
               .OR. &
               (channel%io%ftype(IOMODE_RST) == FTYPE .AND. channel%int%lrst)
       END IF

       IF (PRESENT(set)) THEN
          SELECT CASE(FTYPE)
          CASE(FTYPE_UNDEFINED)
             status = 3200 ! OUTPUT FILE TYPE UNDEFINED
          CASE(FTYPE_ASCII)
             status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_PNETCDF
          CASE(FTYPE_NETCDF)
#else
          CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
             channel%int%netcdf(IOMODE_OUT)%io_pe = set(channel%id)
             channel%int%netcdf(IOMODE_RST)%io_pe = set(channel%id)

             CALL new_attribute(status, channel%att     &
                  , 'channel_io_pe', i=set(channel%id)  &
                  , loverwrite=.TRUE., lout=.TRUE.)
             IF (status /= 0) RETURN

#ifdef HAVE_PNETCDF
          CASE(FTYPE_PNETCDF)
#endif
             ! NOT APPLICABLE
             status = 0
          CASE(FTYPE_GRIB)
             status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
          CASE(FTYPE_HDF4)
             status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
          CASE(FTYPE_HDF5)
             status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef HAVE_CDI
          CASE(FTYPE_CDI_NC)
             channel%int%cdi(IOMODE_OUT)%io_pe = set(channel%id)
             channel%int%cdi(IOMODE_RST)%io_pe = set(channel%id)
             ! NOT APPLICABLE
             status = 0
#endif
#ifdef HAVE_FORPY
          CASE(FTYPE_FORPY)
             channel%int%forpy%io_pe = set(channel%id)
#endif
          CASE DEFAULT
             status = 3212 ! CHANNEL OUTPUT FILE TYPE NOT AVAILABLE
          END SELECT
       END IF

       ls => ls%next
    END DO channel_loop

  END SUBROUTINE setio_channels
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_io
! **********************************************************************
