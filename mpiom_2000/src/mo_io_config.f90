!> module mo_io_config
!> provides an abstraction layer so that unique unit
!> numbers can be retrieved without filling a global table
!> cross-referenced by humans
!>
!> use one of the next_free_unit* functions to find a number suitable for
!> use as logical unit number (lun) and not yet opened or reserved
!> this is further referred to as available
!>
!> @author Thomas Jahns <jahns@dkrz.de>
!> @date 2009-04-09
MODULE mo_io_config
  IMPLICIT NONE
  PRIVATE
  INTEGER :: default_min_lun=1, default_max_lun=100
  INTEGER, PARAMETER :: default_min_reserved_lun=1, &
       default_max_reserved_lun=10
  INTEGER :: min_reserved_lun=default_min_reserved_lun, &
       max_reserved_lun=default_max_reserved_lun
  LOGICAL, ALLOCATABLE, PRIVATE :: reserved_luns(:)
  PUBLIC add_lun_reservation, remove_lun_reservation, add_lun_reservations, &
       remove_lun_reservations, set_reserved_luns, next_free_unit_in_range, &
       next_free_unit, default_min_lun, default_max_lun, setup_io_config
CONTAINS

  !> allocate array of reserved luns, reserves unit from 1 to 10 by default
  SUBROUTINE setup_io_config
    ALLOCATE(reserved_luns(default_min_reserved_lun:default_max_reserved_lun))
    reserved_luns=.TRUE.
  END SUBROUTINE setup_io_config

  !> excludes lun from list of automatically probed luns
  !> @param lun number of lun to be marked as unavailable
  SUBROUTINE add_lun_reservation(lun)
    INTEGER, INTENT(in) :: lun
    CALL adjust_reserved_map(lun)
    reserved_luns(lun) = .TRUE.
  END SUBROUTINE add_lun_reservation

  !> includes lun in list of automatically probed luns
  !> (i.e. cancels previous reservation)
  !> @param lun number of lun to be marked as available
  SUBROUTINE remove_lun_reservation(lun)
    INTEGER, INTENT(in) :: lun
    IF (lun .GE. min_reserved_lun .AND. lun .LE. max_reserved_lun) THEN
      reserved_luns(lun) = .TRUE.
    END IF
  END SUBROUTINE remove_lun_reservation

  !> excludes luns from list of automatically probed luns
  !> @param luns array holding numbers of luns to be marked as unavailable
  !> @param nluns size of array luns
  SUBROUTINE add_lun_reservations(luns, nluns)
    INTEGER, INTENT(in) :: nluns
    INTEGER, DIMENSION(nluns), INTENT(in) :: luns
    INTEGER :: i
    DO i=1,nluns
      CALL add_lun_reservation(luns(i))
    END DO
  END SUBROUTINE add_lun_reservations

  !> includes luns in list of automatically probed luns
  !> @param luns array holding numbers of luns to be marked as available
  !> @param nluns size of array luns
  SUBROUTINE remove_lun_reservations(luns, nluns)
    INTEGER, INTENT(in) :: nluns
    INTEGER, DIMENSION(nluns), INTENT(in) :: luns
    INTEGER :: i
    DO i=1,nluns
      CALL remove_lun_reservation(luns(i))
    END DO
  END SUBROUTINE remove_lun_reservations

  !> set list of luns excluded from automatically probed luns
  !> @param luns array holding numbers of luns to be marked as unavailable
  !> @param nluns size of array luns
  SUBROUTINE set_reserved_luns(luns, nluns)
    INTEGER, INTENT(in) :: nluns
    INTEGER, DIMENSION(nluns), INTENT(in) :: luns
    INTEGER :: i
    max_reserved_lun=MAXVAL(luns)
    min_reserved_lun=MINVAL(luns)
    DEALLOCATE(reserved_luns)
    ALLOCATE(reserved_luns(min_reserved_lun:max_reserved_lun))
    DO i=1,nluns
      CALL add_lun_reservation(luns(i))
    END DO
  END SUBROUTINE set_reserved_luns

  !> returns next available unit in range [min_lun,max_lun] in param unit
  !> @param min_lun do not search for luns less than this
  !> @param max_lun do not search for luns greater than this
  !> @param found set to .false. if no available unit in
  !>   given range could be found
  !> @param unit set to number of available unit iff found is .true.
  SUBROUTINE next_free_unit_in_range(min_lun, max_lun, found, unit)
    INTEGER, INTENT(in) :: min_lun, max_lun
    INTEGER, INTENT(out) :: unit
    LOGICAL, INTENT(out) :: found
    LOGICAL :: opened
    INTEGER :: i
    found = .FALSE.
    DO i=min_lun,max_lun
      IF ((i .LE. max_reserved_lun .AND. i .GE. min_reserved_lun)) THEN
        IF (reserved_luns(i)) CYCLE
      END IF
      INQUIRE(unit=i,opened=opened)
      IF (.NOT.opened) THEN
        unit = i
        found = .TRUE.
        EXIT
      END IF
    END DO
  END SUBROUTINE next_free_unit_in_range

  !> find currently available unit number, the range is determined
  !> by default_min_lun and default_max_lun
  !> @return next available unit
  FUNCTION next_free_unit()
    USE mo_parallel, ONLY: stop_all
    INTEGER :: next_free_unit
    LOGICAL :: found
    CALL next_free_unit_in_range(default_min_lun, &
         default_max_lun, found, next_free_unit)
    IF (.NOT. found) THEN
      CALL stop_all('No free logical unit available for open.')
    END IF
  END FUNCTION next_free_unit

  !> find currently available unit number,  the range is determined
  !> by default_min_lun and default_max_lun and the returned unit is also
  !> immediately reserved
  !> @return next available unit
  FUNCTION reserve_and_get_next_free_unit()
    INTEGER :: reserve_and_get_next_free_unit
    reserve_and_get_next_free_unit = next_free_unit()
    CALL add_lun_reservation(reserve_and_get_next_free_unit)
  END FUNCTION reserve_and_get_next_free_unit

  SUBROUTINE adjust_reserved_map(lun)
    INTEGER, INTENT(in) :: lun
    LOGICAL, ALLOCATABLE :: reserved_luns_copy(:)
    INTEGER :: newmin, newmax
    IF (ALLOCATED(reserved_luns)) THEN
      IF (lun < min_reserved_lun .OR. lun > max_reserved_lun) THEN
        ALLOCATE(reserved_luns_copy(min_reserved_lun:max_reserved_lun))
        reserved_luns_copy = reserved_luns
        newmin=MIN(min_reserved_lun, lun)
        newmax=MAX(max_reserved_lun, lun)
        DEALLOCATE(reserved_luns)
        ALLOCATE(reserved_luns(newmin:newmax))
        reserved_luns=.false.
        reserved_luns(min_reserved_lun:max_reserved_lun)=reserved_luns_copy
        min_reserved_lun=newmin
        max_reserved_lun=newmax
      END IF
    ELSE
      ALLOCATE(reserved_luns(lun:lun))
    END IF
  END SUBROUTINE adjust_reserved_map

! the next function is a port from ECHAM but needs the backtrace
! generating function 'finish' not currently implemented in MPIOM
#if 0
  !> find currently available unit number in given range
  !> @param istart return unit number not less than istart
  !> @param istop return unit number not greater than istop
  !> @return next available unit
  FUNCTION find_next_free_unit(istart,istop) RESULT(unit)
    INTEGER :: istart, istop, unit
    LOGICAL :: found
    INTEGER :: i
    CHARACTER(256) :: info

    CALL next_free_unit_in_range(istart, istop, found, unit)

    IF (.NOT. found) THEN
       WRITE(info,'(a,i2.2,a,i2.2,a)') &
         'No unit in range <',istart,':',istop,'> free.'
       CALL finish('find_next_free_unit',info)
    END IF

  END FUNCTION find_next_free_unit
#endif

END MODULE mo_io_config
