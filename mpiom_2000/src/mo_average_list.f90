 !>
!! Linked list containing averaging information per variable code.
!!
MODULE mo_average_list

#ifndef NO_NEW_IO

  USE mo_kind
  USE mo_varlist

  ! Backend module
#ifndef NOCDI
  USE mo_io_backend_cdi, ONLY: vlist_var_new
#endif
  IMPLICIT NONE
#include "pointer_intent_macro.inc"

  !>
  !! Linked list containing averaging information per variable code.
  !!
  TYPE average_list
     type(average_item), pointer :: head => null() ! , private
     !type(average_item), pointer, private :: last => null()
     type(average_item), pointer :: last => null()
  END TYPE average_list

  !>
  !! Data element containing averaging information per variable code.
  !!
  TYPE average_item
     INTEGER :: code = -1
     INTEGER :: varid = -1
     REAL(dp) :: missval = -9e33_dp
     REAL(dp), ALLOCATABLE :: field(:,:,:)
     !TYPE(average_item), POINTER, private :: next => NULL()
     TYPE(average_item), POINTER :: next => NULL()
   END TYPE average_item

CONTAINS

  !>
  !! Get all but the first element of the list.
  !!
  FUNCTION average_list_tail(this)

    TYPE(average_list), INTENT(in) :: this
    TYPE(average_list) :: average_list_tail

    average_list_tail%head => this%head%next

  END FUNCTION average_list_tail

  !>
  !! Check if list contains an element.
  !!
  FUNCTION average_list_is_not_empty(this)

    TYPE(average_list), INTENT(in) :: this
    LOGICAL :: average_list_is_not_empty

    average_list_is_not_empty = ASSOCIATED(this%head)

  END FUNCTION average_list_is_not_empty

  !>
  !! Add a new element to list.
  !!
  SUBROUTINE average_list_add(this, var, on_io_node, vlist, high_precision, &
       snapshot)

    TYPE (average_list), INTENT(inout) :: this
    TYPE(varlist_element), POINTERINTENT(inout) :: var
    LOGICAL, INTENT(in) :: on_io_node
    INTEGER, INTENT(in) :: vlist
    LOGICAL, INTENT(in) :: high_precision
    LOGICAL, INTENT(in) :: snapshot

    TYPE (average_item), POINTER :: current ! temporary list containing new list element

    ! Bail out if code is already in list
    if(average_list_contains_code(this, var%icode)) then
        return
    end if
    
    ! Create new list element structure

    ALLOCATE(current) ! create new list element
    IF(average_list_is_not_empty(this)) THEN
      this%last%next => current ! append new list element to old list
    ELSE
      this%head => current ! make new element first of list
    END IF
    this%last => current ! make new list element last of list

    ! Add data to new list element

    current%code = var%icode
    if(.not. snapshot .or. var%column_integrate .or. var%surface_integrate) then
       ALLOCATE(current%field(var%i_shape(1), var%i_shape(2), var%i_shape(3)))
       ! @todo This should not be necessary, but is needed for column_integration?
       !       (otherwise GRIB encode overflow)
       ! current%field = 0._dp
    end if

    ! Call backend routines.
    IF(on_io_node) THEN
#ifndef NOCDI
! FIXME: this needs to be reenabled once a default implementation of
! vlist_var_new is available
        current%varid = vlist_var_new(vlist, var%gridid, var%zaxisid, &
            trim(var%name), trim(var%std_name), trim(var%unit), &
            current%code, current%missval, high_precision)
#endif
    END IF

    ! Register variable
    var%registered = .true.

  END SUBROUTINE average_list_add

  !>
  !! Format list for output.
  !!
  !! @todo should check for overflow of output string
  !!
  function average_list_to_string(this) result(the_string)

    TYPE(average_list), intent(in) :: this
    character(2048) :: the_string

    logical :: not_first
    integer :: pos, length
    character(10) :: buffer
    TYPE(average_list) :: current

    the_string = ''

    not_first = .false.
    pos = 1

    the_string(pos:pos) = '['
    pos = pos + 1

    current = this
    do while(average_list_is_not_empty(current))

        if(not_first) then
            the_string(pos:pos+1) = ', '
            pos = pos + 2
        else
            not_first = .true.
        end if

        write(buffer, '(I0)') current%head%code
        length = len_trim(buffer)
        the_string(pos:pos+length-1) = buffer(1:length)
        pos = pos + length

        current = average_list_tail(current)
    end do

    the_string(pos:pos) = ']'

  end function
  
  !>
  !! Search for element with given code number
  !!
  function average_list_contains_code(this, code)
    type(average_list), intent(in) :: this
    integer, intent(in) :: code
    logical :: average_list_contains_code

    type(average_list) :: current

    current = this
    do while(average_list_is_not_empty(current))

        if(current%head%code == code) then
            average_list_contains_code = .true.
            return
        end if 

        current = average_list_tail(current)
    end do

    average_list_contains_code = .false.

  end function average_list_contains_code

#endif/*ndef NO_NEW_IO */

END MODULE mo_average_list
