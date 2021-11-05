! **********************************************************************
MODULE messy_main_channel_forpy
! **********************************************************************

#ifdef HAVE_FORPY

  USE messy_main_constants_mem, ONLY: dp
  USE forpy_mod

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER         :: my_pe
  LOGICAL         :: is_worker = .FALSE.
  LOGICAL         :: lactive   = .FALSE.  ! python kernel launched?
  type(module_py) :: pymod
  type(list)      :: paths

  ! channel.py: - contains all required functions
  !

  ! MAIN ENTRY POINTS TO BE CALLED FROM messy_main_channel_io
  PUBLIC :: ch_forpy_init_pio ! initialistion of distributed workspace
  !                           ! channels are distributed among pe-s
  PUBLIC :: ch_forpy_init_io  ! launch python kernel
  !
  PUBLIC :: ch_forpy_write_header ! add global attributes
  PUBLIC :: ch_forpy_write_time   ! add time information
  PUBLIC :: ch_forpy_write_data   ! convert channel objects (as dictionary)
  PUBLIC :: ch_forpy_finish_io  !
  !
  !PRIVATE :: fpy_master      ! throw entire channel (as dictionary to python)
  !                           ! and call corresponding python function
  !PRIVATE :: fpy_att2dict    ! convert list of attributes to python dictionary
  !PRIVATE :: fpy_dim2dict    ! convert dimension to python dictionary
  !PRIVATE :: fpy_dimvar2dict ! convert dimension variable to python dictionary
  !PRIVATE :: fpy_repr2dict   ! convert representation to python dictionary
  !PRIVATE :: fpy_channel_object2dict
  !PRIVATE :: fpy_channel_object_data2dict
  !PRIVATE :: fpy_channel2dict

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE ch_forpy_init_pio(p_nprocs, p_pe)

    ! I/O
    INTEGER, INTENT(IN) :: p_nprocs
    INTEGER, INTENT(IN) :: p_pe

    my_pe = p_pe

  END SUBROUTINE ch_forpy_init_pio
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_forpy_init_io(status, IOMODE, channel, AMODE)

    USE messy_main_channel, ONLY: t_channel, IOMODE_OUT, AMODE_WRITE

    IMPLICIT NONE

    ! I/O
    INTEGER,           INTENT(OUT) :: status
    INTEGER,           INTENT(IN)  :: IOMODE, AMODE
    TYPE(t_channel),   POINTER     :: channel    ! INTENT(INOUT)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_forpy_init_io'
    INTEGER                     :: ierror

    status = 0
    IF (IOMODE /= IOMODE_OUT)  RETURN
    IF (AMODE  /= AMODE_WRITE) RETURN

    is_worker = is_worker .OR. (channel%int%forpy%io_pe == my_pe)

    IF (is_worker .AND. (.NOT. lactive)) THEN

       ierror = forpy_initialize()

       ! forpy_initialize returns NO_NUMPY_ERROR if numpy could not be imported
       IF (ierror == NO_NUMPY_ERROR) then
          status = 5003
          RETURN
       END IF
       IF (ierror /= 0) THEN
          status = 5004
          return
       END if

       ! APPEND . TO ENABLE LOADING OWN MODULES FROM WORKDIR
       ierror = get_sys_path(paths)
       if (ierror /= 0) THEN
          status = 5005
          RETURN
       endif
       ierror = paths%append(".")
       if (ierror /= 0) THEN
          status = 5006
          RETURN
       endif
       ! NOTE: Only one python kernel can be launched on each task.
       !       This implies that only one python moduel can be loaded.
       !       Yet, the python module can conatin several functions.
       !       The latter are selected by the name of the channel, see
       !       fpy_master ...
       ierror = import_py(pymod, "channel") ! channel.py
       if (ierror /= 0) THEN
          status = 5007
          RETURN
       endif

       lactive = .TRUE.
    END IF

  END SUBROUTINE ch_forpy_init_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_forpy_write_header(status, IOMODE, channel)


    USE messy_main_channel,            ONLY: t_channel, IOMODE_OUT, GATT
    USE messy_main_channel_attributes, ONLY: t_attribute_list

    IMPLICIT NONE

    ! I/O
    INTEGER,                      INTENT(OUT)   :: status
    INTEGER,                      INTENT(IN)    :: IOMODE
    TYPE(t_channel),              INTENT(INOUT) :: channel

    ! LOCAL
    INTEGER :: ierror

    status = 0
    IF (IOMODE /= IOMODE_OUT) RETURN
    IF (channel%int%forpy%io_pe /= my_pe) RETURN

    CALL fpy_att2dict(ierror, GATT, channel%int%forpy%ga)
    IF (ierror /= 0) THEN
       status = 5012
       RETURN
    END IF

  END SUBROUTINE ch_forpy_write_header
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_forpy_write_time(status, channel, dim_time)

    USE messy_main_channel,            ONLY: t_channel
    USE messy_main_channel_dimensions, ONLY: t_dimension

    IMPLICIT NONE

    ! I/O
    INTEGER,           INTENT(OUT) :: status
    TYPE(t_channel),   POINTER     :: channel ! INTENT(INOUT)
    TYPE(t_dimension), POINTER     :: dim_time

    status = 0
    !IF (IOMODE /= IOMODE_OUT) RETURN
    IF (channel%int%forpy%io_pe /= my_pe) RETURN

    CALL fpy_dim2dict(status, dim_time, channel%int%forpy%time)

  END SUBROUTINE ch_forpy_write_time
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_forpy_write_data(status, IOMODE, channel, object, ptr, jsnd)

    USE messy_main_channel, ONLY: t_channel, t_channel_object, IOMODE_OUT &
                                , SND_TEXT

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)   :: status
    INTEGER,                      INTENT(IN)    :: IOMODE
    TYPE(t_channel),              INTENT(IN)    :: channel ! INTENT(IN)
    TYPE(t_channel_object),       INTENT(INOUT) :: object
    REAL(DP), DIMENSION(:,:,:,:), POINTER       :: ptr
    INTEGER,                      INTENT(IN)    :: jsnd  ! OUTPUT DATA TYPE

    status = 0
    IF (IOMODE /= IOMODE_OUT) RETURN
    IF (channel%int%forpy%io_pe /= my_pe) RETURN

    CALL fpy_channel_object2dict(status, object, object%int%forpy%d(jsnd) &
         , TRIM(SND_TEXT(jsnd, IOMODE)))
    IF (status /= 0) RETURN

    CALL fpy_channel_object_data2dict(status, ptr &
         , object%repr%rank                       &
         , object%int%forpy%d(jsnd))

  END SUBROUTINE ch_forpy_write_data
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_forpy_finish_io(status, IOMODE, channel, lfinal)

    USE messy_main_channel, ONLY: t_channel, t_channel_object_list &
                                , t_channel_object, IOMODE_OUT, SND_MAXLEN

    IMPLICIT NONE

    ! I/O
    INTEGER,           INTENT(OUT) :: status
    INTEGER,           INTENT(IN)  :: IOMODE
    TYPE(t_channel),   POINTER     :: channel    ! INTENT(INOUT)
    LOGICAL,           INTENT(IN)  :: lfinal

    ! LOCAL
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: jsnd

    status = 0
    IF (IOMODE /= IOMODE_OUT) RETURN
    IF (channel%int%forpy%io_pe /= my_pe) RETURN

    IF (.NOT. lfinal) THEN
       CALL fpy_master(status, channel)
       IF (status /= 0) RETURN
    END IF

    le => channel%list
    object_loop: DO
       IF (.NOT. ASSOCIATED(le)) EXIT
       object => le%this
       DO jsnd = 1, SND_MAXLEN
          IF ( .NOT. is_null(object%int%forpy%d(jsnd)) ) &
               CALL object%int%forpy%d(jsnd)%destroy()
       END DO
       le => le%next
    END DO object_loop
!    IF (.NOT. is_null(channel%int%forpy%ga)) &
!         CALL channel%int%forpy%ga%destroy()
    IF (.NOT. is_null(channel%int%forpy%time)) &
         CALL channel%int%forpy%time%destroy()

    IF (lfinal) THEN

       IF (is_worker .AND. lactive ) THEN
          CALL paths%destroy
          CALL pymod%destroy
          CALL forpy_finalize
          lactive = .FALSE.
       END IF

    END IF

  END SUBROUTINE ch_forpy_finish_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! HELPER ROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE fpy_master(status, channel)

    USE messy_main_channel, ONLY: t_channel, STRLEN_CHANNEL
    USE messy_main_tools,   ONLY: strcrack

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,         INTENT(OUT)   :: status
    TYPE(t_channel), INTENT(IN)    :: channel ! INTENT(IN)

    ! LOCAL
    INTEGER      :: ierror
    type(dict)   :: d
    type(tuple)  :: args
    type(object) :: status_obj
    INTEGER      :: n, j

    status = 0
    IF (channel%int%forpy%io_pe /= my_pe) RETURN

    CALL fpy_channel2dict(ierror, channel, d)
    if (ierror /= 0) then
       status = 5011
       return
    endif

    ierror = tuple_create(args, 1)
    if (ierror /= 0) then
       status = 5008
       return
    endif
    ierror = args%setitem(0, d)
    if (ierror /= 0) then
       status = 5008
       return
    endif

    ierror = call_py(status_obj, pymod, 'fpy_main', args)
    if (ierror /= 0) then
       status = 5009
       return
    endif

    ierror = cast(status, status_obj)
    if (ierror /= 0) then
       status = 5010
       return
    endif

    CALL status_obj%destroy
    CALL d%destroy
    CALL args%destroy

    status = 0

  END SUBROUTINE fpy_master
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE fpy_att2dict(status, a, d)

    USE messy_main_channel_attributes, ONLY: t_attribute_list, t_attribute &
                                           , TYPE_STRING, TYPE_INTEGER &
                                           , TYPE_REAL_DP

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_attribute_list), POINTER     :: a ! INTENT(IN)
    TYPE(dict),             INTENT(OUT) :: d

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'fpy_att2dict'
    INTEGER                         :: ierror
    TYPE(t_attribute_list), POINTER :: ai
    TYPE(t_attribute_list), POINTER :: ae
    TYPE(t_attribute),      POINTER :: att

    status = 0
    NULLIFY(ai, ae, att)

    ierror = dict_create(d) ! Python: d = {}
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ai => a
    NULLIFY(ae)
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       !
       att => ai%this
       !
       SELECT CASE(att%type)
       CASE(TYPE_STRING)
          ierror = d%setitem(TRIM(att%name), TRIM(att%c))
       CASE(TYPE_INTEGER)
          ierror = d%setitem(TRIM(att%name), att%i)
       CASE(TYPE_REAL_DP)
          ierror = d%setitem(TRIM(att%name), att%r)
       END SELECT
       !
       if (ierror /= 0) then
          CALL nerr(ierror, status, substr)
          RETURN
       end if
       !
       ae => ai
       ai => ai%next
    END DO

  END SUBROUTINE fpy_att2dict
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE fpy_dim2dict(status, di, d)

    USE messy_main_channel_dimensions, ONLY: t_dimension
    USE messy_main_channel_dimvar,     ONLY: t_dimvar, t_dimvar_list

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_dimension),      INTENT(IN)  :: di
    TYPE(dict),             INTENT(OUT) :: d

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'fpy_dim2dict'
    INTEGER    :: ierror
    TYPE(list) :: l
    TYPE(dict) :: dv
    TYPE(t_dimvar_list), POINTER :: ai
    TYPE(t_dimvar_list), POINTER :: ae
    TYPE(t_dimvar),      POINTER :: dimvar

    status = 0
    NULLIFY(ai, ae, dimvar)

    ierror = dict_create(d) ! Python: d = {}
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("name", TRIM(di%name))
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("domain_id", di%dom_id)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("length", di%len)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("is_time", di%ltime)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = list_create(l)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ai => di%var
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       !
       dimvar => ai%this
       CALL fpy_dimvar2dict(ierror, dimvar, dv)
       if (ierror /= 0) then
          CALL nerr(ierror, status, substr)
          RETURN
       end if

       ierror = l%append(dv)
       if (ierror /= 0) then
          CALL nerr(ierror, status, substr)
          RETURN
       end if

       CALL dv%destroy()
       !
       ae => ai
       ai => ai%next
    END DO

    ierror = d%setitem("dimvars", l)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    CALL l%destroy()

  END SUBROUTINE fpy_dim2dict
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE fpy_dimvar2dict(status, v, d)

    USE messy_main_channel_dimvar, ONLY: t_dimvar

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,        INTENT(OUT) :: status
    TYPE(t_dimvar), INTENT(IN)  :: v
    TYPE(dict),     INTENT(OUT) :: d

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'fpy_dimvar2dict'
    INTEGER       :: ierror
    TYPE(dict)    :: da
    type(ndarray) :: dv

    status = 0

    ierror = dict_create(d) ! Python: d = {}
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("name", TRIM(v%name))
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    IF (ASSOCIATED(v%val)) THEN
       ierror = ndarray_create(dv, v%val)
       if (ierror /= 0) then
          CALL nerr(ierror, status, substr)
          RETURN
       end if
       ierror = d%setitem("data", dv)
       if (ierror /= 0) then
          CALL nerr(ierror, status, substr)
          RETURN
       end if
    END IF

    CALL fpy_att2dict(ierror, v%att, da)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if
    ierror = d%setitem("attributes", da)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    CALL dv%destroy()
    CALL da%destroy()

  END SUBROUTINE fpy_dimvar2dict
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE fpy_repr2dict(status, r, d)

    USE messy_main_channel_repr, ONLY: t_representation

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_representation), INTENT(IN)  :: r
    TYPE(dict),             INTENT(OUT) :: d

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'fpy_repr2dict'
    INTEGER       :: ierror
    TYPE(dict)    :: dd
    type(tuple)   :: gdim
    type(tuple)   :: di
    INTEGER       :: i, j

    status = 0

    ierror = dict_create(d) ! Python: d = {}
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("name", TRIM(r%name))
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("rank", r%rank)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("link", TRIM(r%link))
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("axis", TRIM(r%axis))
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = tuple_create(gdim, r%rank)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    DO i=1, r%rank
       ierror = gdim%setitem(i-1,r%gdimlen(i))
       if (ierror /= 0) then
          CALL nerr(ierror, status, substr)
          RETURN
       end if
    END DO
    ierror = d%setitem("gdimlen", gdim)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = tuple_create(di, r%rank)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if
    j = 0
    DO i=1, r%rank
       IF (ASSOCIATED(r%dim(i)%ptr)) THEN
          CALL fpy_dim2dict(status, r%dim(i)%ptr, dd)
          IF (status /= 0) RETURN
          ierror = di%setitem(j, dd)
          if (ierror /= 0) then
             CALL nerr(ierror, status, substr)
             RETURN
          end if
          CALL dd%destroy()
          j = j + 1
       END IF
    END DO
    ierror = d%setitem("dimensions", di)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    CALL gdim%destroy()
    CALL di%destroy()

  END SUBROUTINE fpy_repr2dict
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE fpy_channel_object2dict(status, co, d, suffix)

    USE messy_main_channel,            ONLY: t_channel_object

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_channel_object), INTENT(IN)  :: co
    TYPE(dict),             INTENT(OUT) :: d
    CHARACTER(LEN=*),       INTENT(IN)  :: suffix

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'fpy_channel_object2dict'
    INTEGER    :: ierror
    TYPE(dict) :: da, dr

    status = 0

    ierror = dict_create(d) ! Python: d = {}
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("name", TRIM(co%name)//TRIM(suffix))
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    CALL fpy_att2dict(ierror, co%att, da)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("attributes", da)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    CALL fpy_repr2dict(status, co%repr, dr)
    IF (status /= 0) RETURN

    ierror = d%setitem("representation", dr)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ! co%data needs to be added after gather on BMIL ...

    CALL da%destroy()
    CALL dr%destroy()

  END SUBROUTINE fpy_channel_object2dict
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE fpy_channel_object_data2dict(status, dat, irank, d)

    IMPLICIT NONE

    ! I/O
    INTEGER,                      INTENT(OUT)   :: status
    REAL(DP), DIMENSION(:,:,:,:), POINTER       :: dat ! INTENT(IN)
    INTEGER,                      INTENT(IN)    :: irank
    TYPE(dict),                   INTENT(INOUT) :: d

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'fpy_channel_object_data2dict'
    INTEGER       :: ierror
    type(ndarray) :: nd
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: p4
    REAL(DP), DIMENSION(:,:,:),   POINTER :: p3
    REAL(DP), DIMENSION(:,:),     POINTER :: p2
    REAL(DP), DIMENSION(:),       POINTER :: p1
    REAL(DP),                     POINTER :: p0

    status = 0

    SELECT CASE(irank)
    CASE(0)
       p0 => dat(1,1,1,1)
       ierror = ndarray_create(nd, (/p0/))
    CASE(1)
       p1 => dat(:,1,1,1)
       ierror = ndarray_create(nd, p1)
    CASE(2)
       p2 => dat(:,:,1,1)
       ierror = ndarray_create(nd, p2)
    CASE(3)
       p3 => dat(:,:,:,1)
       ierror = ndarray_create(nd, p3)
    CASE(4)
       p4 => dat(:,:,:,:)
       ierror = ndarray_create(nd, p4)
    END SELECT
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("data", nd)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

  END SUBROUTINE fpy_channel_object_data2dict
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE fpy_channel2dict(status, c, d)

    USE messy_main_channel, ONLY: t_channel, t_channel_object &
                                , t_channel_object_list, SND_MAXLEN

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER,         INTENT(OUT) :: status
    TYPE(t_channel), INTENT(IN)  :: c
    TYPE(dict),      INTENT(OUT) :: d

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'fpy_channel2dict'
    INTEGER    :: ierror
    TYPE(dict) :: da
    type(list) :: ob
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: jsnd

    status = 0

    ierror = dict_create(d) ! Python: d = {}
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("name", TRIM(c%name))
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("domain_id", c%dom_id)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("global_attributes", c%int%forpy%ga)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    CALL fpy_att2dict(status, c%att, da)
    IF (status /= 0) RETURN

    ierror = d%setitem("channel_attributes", da)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = d%setitem("time",c%int%forpy%time)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    ierror = list_create(ob)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if
    !
    le => c%list
    object_loop: DO
       IF (.NOT. ASSOCIATED(le)) EXIT
       object => le%this
       DO jsnd = 1, SND_MAXLEN
          IF ( .NOT. is_null(object%int%forpy%d(jsnd)) ) THEN
             ierror = ob%append(object%int%forpy%d(jsnd))
             if (ierror /= 0) then
                CALL nerr(ierror, status, substr)
                RETURN
             end if
          END IF
       END DO
       le => le%next
    END DO object_loop

    ierror = d%setitem("objects", ob)
    if (ierror /= 0) then
       CALL nerr(ierror, status, substr)
       RETURN
    end if

    CALL da%destroy()
    CALL ob%destroy()

  END SUBROUTINE fpy_channel2dict
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE nerr(ierror, sout, substr)

    USE messy_main_constants_mem, ONLY: iouerr
    USE messy_main_tools,         ONLY: int2str

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(IN)  :: ierror
    INTEGER,          INTENT(OUT) :: sout
    CHARACTER(LEN=*), INTENT(IN)  :: substr

    ! LOCAL
    INTEGER :: iou
    CHARACTER(LEN=4) :: pestr
    LOGICAL :: opened

    IF (ierror == 0) RETURN

    sout = 5013

    WRITE(iouerr,*) TRIM(substr),': *** forpy ERROR on PE ' &
         , my_pe,' : ', ierror
    CALL err_print

    DO iou=100,300
       INQUIRE(unit=iou,opened=opened)
       IF (.NOT.opened) EXIT
    END DO
    CALL int2str(pestr, my_pe, '0', '_')
    OPEN(iou, FILE='ERROR.'//pestr//'.forpy', STATUS='UNKNOWN')
    WRITE(iou,*) TRIM(substr),': *** forpy ERROR: '
    CLOSE(iou)

  END SUBROUTINE nerr
  ! -------------------------------------------------------------------

#endif

! **********************************************************************
END MODULE messy_main_channel_forpy
! **********************************************************************
