  !>
  !! Simple linked list implementation.
  !!
  
  !>
  !! Linked list element
  !!
  TYPE :: list_cell
     PRIVATE
     TYPE(list_item) :: item
     TYPE(list_cell), POINTER :: next => NULL()
  END TYPE list_cell
  
  !>
  !! Linked list as separate data type to allow for painless initialization.
  !!
  TYPE :: list_type
     PRIVATE
     TYPE(list_cell), POINTER :: head => NULL()
     TYPE(list_cell), POINTER :: last => NULL()
  END TYPE list_type
  
  !>
  !! Private methods
  !!
  PRIVATE :: &
     list_head_internal, &
     list_tail_internal, &
     list_is_not_empty_internal, &
     list_get_item_by_id_internal, &
     list_add_internal, &
     list_to_string_internal
     
  !>
  !! Public interfaces
  !!
  INTERFACE list_head
     MODULE PROCEDURE list_head_internal
  END INTERFACE list_head

  INTERFACE list_tail
     MODULE PROCEDURE list_tail_internal
  END INTERFACE list_tail
  
  INTERFACE list_is_not_empty
     MODULE PROCEDURE list_is_not_empty_internal
  END INTERFACE list_is_not_empty

  INTERFACE list_get_item_by_id
     MODULE PROCEDURE list_get_item_by_id_internal
  END INTERFACE list_get_item_by_id
  
  INTERFACE list_add
     MODULE PROCEDURE list_add_internal
  END INTERFACE list_add
  
  INTERFACE list_to_string
     MODULE PROCEDURE list_to_string_internal
  END INTERFACE list_to_string
  
CONTAINS

  !>
  !! Get contents of first list element.
  !!
  FUNCTION list_head_internal(this)

    TYPE(list_type), INTENT(in) :: this
    TYPE(list_item) :: list_head_internal

    list_head_internal = this%head%item

  END FUNCTION list_head_internal
  
  !>
  !! Get all but the first element of the list.
  !!
  FUNCTION list_tail_internal(this)

    TYPE(list_type), INTENT(in) :: this
    TYPE(list_type) :: list_tail_internal

    list_tail_internal%head => this%head%next
    list_tail_internal%last => this%last

  END FUNCTION list_tail_internal

  !>
  !! Check if list contains an element.
  !!
  FUNCTION list_is_not_empty_internal(this)

    TYPE(list_type), INTENT(in) :: this
    LOGICAL :: list_is_not_empty_internal

    list_is_not_empty_internal = ASSOCIATED(this%head)

  END FUNCTION list_is_not_empty_internal

  !>
  !! Get a list element by item ID.
  !!
  FUNCTION list_get_item_by_id_internal(this, id, item)

    TYPE(list_type), INTENT(in) :: this
    INTEGER, INTENT(in) :: id
    TYPE(list_item), INTENT(out) :: item
    LOGICAL :: list_get_item_by_id_internal

    TYPE(list_type) :: current

    current = this
    DO WHILE(list_is_not_empty_internal(current))
       item = list_head_internal(current)
       IF(list_item_get_id(item) == id) THEN
          list_get_item_by_id_internal = .true.
          RETURN
       END IF
       current = list_tail_internal(current)
    END DO

    list_get_item_by_id_internal = .false.

  END FUNCTION list_get_item_by_id_internal

  !>
  !! Add a new element to list.
  !!
  SUBROUTINE list_add_internal(this, item)

    TYPE (list_type), INTENT(inout) :: this
    TYPE(list_item), INTENT(in) :: item

    TYPE (list_cell), POINTER :: current ! temporary list containing new list element

    ! Create new list element structure

    ALLOCATE(current) ! create new list element
    IF(list_is_not_empty_internal(this)) THEN
       this%last%next => current ! append new list element to old list
    ELSE
       this%head => current ! make new element first of list
    END IF
    this%last => current ! make new list element last of list

    ! Add data to new list element

    current%item = item

  END SUBROUTINE list_add_internal

  !>
  !! Format list for output.
  !!
  !! @todo should check for overflow of output string
  !!
  FUNCTION list_to_string_internal(this) RESULT(the_string)

    TYPE(list_type), INTENT(in) :: this
    CHARACTER(1024) :: the_string

    LOGICAL :: not_first
    INTEGER :: pos, length
    CHARACTER(1024) :: buffer
    TYPE(list_type) :: current

    the_string = ''

    not_first = .FALSE.
    pos = 1

    the_string(pos:pos) = '['
    pos = pos + 1

    current = this
    DO WHILE(list_is_not_empty_internal(current))

       IF(not_first) THEN
          the_string(pos:pos+1) = ', '
          pos = pos + 2
       ELSE
          not_first = .TRUE.
       END IF

       buffer = list_item_to_string(list_head_internal(current))
       length = LEN_TRIM(buffer)
       the_string(pos:pos+length-1) = buffer(1:length)
       pos = pos + length

       current = list_tail_internal(current)
    END DO

    the_string(pos:pos) = ']'

  END FUNCTION list_to_string_internal
