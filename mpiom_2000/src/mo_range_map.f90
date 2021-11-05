!>
!! mpiom-latest-cmip5 - mo_range_map 
!!
!! @author Karl-Hermann Wieners, Max Planck Institute for Meteorology
!!
!! $Revision$
!! $Date$
!!
module mo_range_map

    implicit none
    
    type list_item
        integer :: key_from, key_to, value
    end type
    
include 'mo_linked_list.inc'

    !>
    !! Create combined key from integer range.
    !!
    function key(key_from, key_to) 
        integer, intent(in) :: key_from, key_to
        integer :: key
        key = key_from*100000 + key_to
    end function key
    
    !>
    !! Checks range of key generation input.
    !!
    !! Allows range from 0-21473 to 0-99999
    !! to stay with 4 byte integer range for combined key.
    !!
    function is_range_invalid(key_from, key_to) 
        integer, intent(in) :: key_from, key_to
        logical :: is_range_invalid
        is_range_invalid = key_from < 0 .or. 21473 < key_from .or. &
                           key_to < 0 .or. 99999 < key_to
    end function is_range_invalid

    !>
    !! Define unique ID for list elements.
    !!
    function list_item_get_id(this)
        type(list_item), intent(in) :: this
        integer :: list_item_get_id
        list_item_get_id = key(this%key_from, this%key_to)
    end function list_item_get_id
    
  !>
  !! Format list item for output.
  !!
  !! @todo should check for overflow of output string
  !!
  FUNCTION list_item_to_string(this) RESULT(the_string)

    TYPE(list_item), INTENT(in) :: this
    CHARACTER(50) :: the_string
    
    write(the_string, "('((', I0, ',', I0, '), ', I0, ')')") &
         this%key_from, this%key_to, this%value

  END FUNCTION list_item_to_string
  
  subroutine range_map_add(this, key_from, key_to, value)
      type(list_type), intent(inout) :: this
      integer, intent(in) :: key_from, key_to, value
      if(is_range_invalid(key_from, key_to)) &
           stop 'mo_range_map::range_map_add: key size exceeded'
      call list_add(this, list_item(key_from, key_to, value))
  end subroutine range_map_add
  
  function range_map_get(this, key_from, key_to, value)
      type(list_type), intent(inout) :: this
      integer, intent(in) :: key_from, key_to
      integer, intent(out) :: value
      logical :: range_map_get
      type(list_item) :: item
      range_map_get = list_get_item_by_id(this, key(key_from, key_to), item)
      if(range_map_get) then
          value = item%value
      end if
  end function range_map_get
  
end module mo_range_map
