!>
!! Linked list containing input/output settings per data file.
!!
MODULE mo_file_list

#ifndef NO_NEW_IO

  use mo_average_list
  use mo_varlist
  use mo_model_time, only: is_partial_day_supported
  use mo_parallel, only: stop_all

  ! Backend module
#ifndef NOCDI
  USE mo_io_backend_cdi, ONLY: vlist_new, vlist_get_taxis
#endif
  implicit none
#include "pointer_intent_macro.inc"

  private

  !>
  !! Constants for type of averaging.
  !!
  integer, parameter :: offset_span = 100
  ! Factors must be binary disjunct i.e. powers of 2
  INTEGER, PARAMETER :: snapshot_factor = 1
  INTEGER, PARAMETER :: unmasked_factor = 2
  INTEGER, PARAMETER :: hires_factor = 4

  ! Caution: new file types must also be checked in file_list_add!
  INTEGER, PARAMETER, public :: file_restart = 99
  INTEGER, PARAMETER, public :: file_final = 98
  INTEGER, PARAMETER, public :: file_initial = 90
  INTEGER, PARAMETER, public :: file_daily = 1
  INTEGER, PARAMETER, public :: file_monthly = 2
  INTEGER, PARAMETER, public :: file_annual = 3
  INTEGER, PARAMETER, public :: file_timestep = 4
  INTEGER, PARAMETER, public :: file_12h = 5
  INTEGER, PARAMETER, public :: file_6h = 6
  INTEGER, PARAMETER, public :: file_3h = 7
  INTEGER, PARAMETER, public :: file_2h = 8
  INTEGER, PARAMETER, public :: file_1h = 9

  !>
  !! Linked list containing input/output settings per data file.
  !!
  TYPE, PUBLIC :: file_list
     logical :: on_io_node = .false.
     type(file_item), pointer :: head => null() ! , private
     !type(file_item), pointer, private :: last => null()
     type(file_item), pointer :: last => null()
  end type file_list

  !>
  !! Linked list element containing input/output settings per data file.
  !!
  TYPE, PUBLIC :: file_item
     INTEGER              :: TYPE = 0
     LOGICAL              :: snapshot = .false.
     LOGICAL              :: unmasked = .false.
     LOGICAL              :: hires = .false.
     CHARACTER(len=128)   :: name = ''
     CHARACTER(len=3)     :: FORMAT = ''
     INTEGER              :: streamid = -1
     INTEGER              :: vlistid = -1
     INTEGER              :: taxisid = -1
     INTEGER              :: nstep = 0
     TYPE (average_list)  :: avglist
     !TYPE (file_item), POINTER, private :: next => NULL()
     TYPE (file_item), POINTER :: next => NULL()
  END TYPE file_item

  !>
  !! Structure describing name list format for I/O configuration.
  !!
  !! This structure is only needed for reading I/O configuration from name list.
  !! As soon as the name list is read, data is re-formatted into lists.
  !!
  TYPE, PUBLIC :: file_data
     INTEGER              :: TYPE=0
     CHARACTER(len=128)   :: name=''
     CHARACTER(len=3)     :: FORMAT=''
     INTEGER              :: codes(255)=0
  END TYPE file_data

  ! Declare public elements

  PUBLIC :: file_list_set_on_io_node
  PUBLIC :: file_list_tail
  PUBLIC :: file_list_is_not_empty
  PUBLIC :: file_list_add
  PUBLIC :: file_item_to_string

CONTAINS

  !>
  !! Determine whether list instance actually performs I/O.
  !!
  !! Datatype default is false i.e. no I/O is performed.
  !! Currently exactly one list instance must have this flag set.
  !!
  SUBROUTINE file_list_set_on_io_node(this, on_io_node)

    TYPE(file_list), INTENT(inout) :: this
    logical, intent(in) :: on_io_node

    this%on_io_node = on_io_node

  end subroutine

  !>
  !! Get all but the first element of an I/O list.
  !!
  FUNCTION file_list_tail(this)

    TYPE(file_list), INTENT(in) :: this
    TYPE(file_list) :: file_list_tail

    file_list_tail%head => this%head%next

  END FUNCTION file_list_tail

  !>
  !! Check if I/O list contains an element.
  !!
  FUNCTION file_list_is_not_empty(this)

    TYPE(file_list), INTENT(in) :: this
    LOGICAL :: file_list_is_not_empty

    file_list_is_not_empty = ASSOCIATED(this%head)

  END FUNCTION file_list_is_not_empty

  !>
  !! Add an element to I/O configuration using data from name list.
  !!
  SUBROUTINE file_list_add(this, data, the_varlist)

    TYPE(file_list), intent(inout) :: this !list to which data is added
    TYPE(file_data), INTENT(in) :: data ! new list element
    TYPE(varlist), POINTERINTENT(in) :: the_varlist

    TYPE(file_item), POINTER :: current ! temporary list containing new list element
    type(varlist_element), pointer :: var

    integer :: i, file_type, file_options
    character(1000) :: message
    logical :: snapshot, unmasked, hires, supported

    ! Check pre-conditions

    ! Extract base type and options from type number
    file_type = modulo(data%type, offset_span)
    file_options = data%type/offset_span

    ! Set snapshot flag if type number is within proper offset range.
    ! Assume implicit snapshot request for restart and timestep
    snapshot = &
         iand(file_options, snapshot_factor) > 0 .or. &
         file_type == file_restart .or. &
         file_type == file_final .or. &
         file_type == file_initial .or. &
         file_type == file_timestep

    ! Set unmasked flag if type number is within proper offset range.
    ! Implicitly assume masking for restart files.
    unmasked = &
         iand(file_options, unmasked_factor) > 0 .or. &
         file_type == file_restart

    ! Set high resolution flag if type number is within proper offset range.
    ! Implicitly assume hires for restart files.
    hires = &
         iand(file_options, hires_factor) > 0 .or. &
         file_type == file_restart

    ! Supported file types
    select case(file_type)
    case(file_annual, file_monthly, file_daily, file_timestep, &
         file_restart, file_initial, file_final)
                    supported = .true.
    case(file_12h); supported = is_partial_day_supported(2)
    case(file_6h);  supported = is_partial_day_supported(4)
    case(file_3h);  supported = is_partial_day_supported(8)
    case(file_2h);  supported = is_partial_day_supported(12)
    case(file_1h);  supported = is_partial_day_supported(24)
    case default;   supported = .false.
    end select

    ! Check for fatal error
    ! File type is not supported or partial day is incompatible.
    if(.not. supported) then
        WRITE(message, '(A, I0, A, A, A)') &
             "Oops: file type ", data%type, " for file '", trim(data%name), &
             "' is not supported or incompatible with model time step."
        CALL stop_all(trim(message))
    end if

    ! Create new list element structure

    ALLOCATE(current) ! create new list element
    if(file_list_is_not_empty(this)) then
        this%last%next => current ! append new list element to old list
    else
        this%head => current ! make new element first of list
    end if
    this%last => current ! make new list element last of list

    ! Add data to new list element

    current%type = file_type
    current%name = data%name
    current%format = data%format
    current%snapshot = snapshot
    current%unmasked = unmasked
    current%hires = hires

#ifndef NOCDI
! FIXME: this needs to be reenabled once a default implementation of
! vlist_var_new is available
    ! Call backend routines
    IF(this%on_io_node) THEN
      current%vlistid = vlist_new()
      current%taxisid = vlist_get_taxis(current%vlistid)
    END IF
#endif
    ! Create dependent lists

    ! Make sure we only allocate memory for relevant codes.
    DO i = 1, COUNT(data%codes /= 0)
      var => varlist_get_element_ref(the_varlist, data%codes(i))
      if(associated(var)) then
         ! Check use of scaled variables in restart files.
         ! Restart files should use unaltered data, so factoring is not applied.
         ! This causes inconsistent unit information,
         ! so these variables may not be used in restart files.
         if(current%type == file_restart .and. var%factor /= 1.0_dp) then
            WRITE(message, '(A, I0, A, A, A)') &
                 "Oops: variable definition for code ", var%icode, &
                 " has a scaling factor and may not be used in restart file '", &
                 trim(current%name), "'."
            CALL stop_all(trim(message))
         end if
         CALL average_list_add(current%avglist, var, this%on_io_node, &
              current%vlistid, current%hires, current%snapshot)
      end if
    END DO

  END SUBROUTINE file_list_add

  !>
  !! Format list for output.
  !!
  function file_item_to_string(this, verbose) result(the_string)

    type(file_item), intent(in) :: this
    logical, intent(in), optional :: verbose
    character(2048) :: the_string

    logical :: do_verbose
    integer :: pos, length
    character(2048) :: buffer

    ! Default to non-verbose output.
    do_verbose = merge(verbose, .false., present(verbose))

    the_string = ''
    pos = 1

    ! prefix
    if(do_verbose) then
        the_string(pos:pos+8) = 'file_item'
        pos = pos + 9
    end if
    the_string(pos:pos) = '('
    pos = pos + 1

    ! %name
    if(do_verbose) then
        the_string(pos:pos+5) = 'name: '
        pos = pos + 6
    end if
    length = len_trim(this%name)
    the_string(pos:pos) = "'"
    the_string(pos+1:pos+length) = this%name(1:length)
    the_string(pos+length+1:pos+length+3) = "'"
    pos = pos + length+2

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %format
    if(do_verbose) then
        the_string(pos:pos+7) = 'format: '
        pos = pos + 8
    end if
    length = len_trim(this%format)
    the_string(pos:pos) = "'"
    the_string(pos+1:pos+length) = this%format(1:length)
    the_string(pos+length+1:pos+length+1) = "'"
    pos = pos + length+2

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %snapshot
    if(this%snapshot) then
       the_string(pos:pos+7) = 'snapshot'
       pos = pos + 8
    else
       the_string(pos:pos+7) = 'averaged'
       pos = pos + 8
    end if

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %unmasked
    if(this%unmasked) then
       the_string(pos:pos+7) = 'unmasked'
       pos = pos + 8
    else
       the_string(pos:pos+5) = 'masked'
       pos = pos + 6
    end if

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %hires
    if(this%hires) then
       the_string(pos:pos+1) = 'hi'
       pos = pos + 2
    else
       the_string(pos:pos+1) = 'lo'
       pos = pos + 2
    end if

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %avglist
    if(do_verbose) then
        the_string(pos:pos+8) = 'avglist: '
        pos = pos + 9
    end if
    buffer = average_list_to_string(this%avglist)
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! suffix
    the_string(pos:pos) = ')'

  end function file_item_to_string

#endif/*ndef NO_NEW_IO */

END MODULE mo_file_list
