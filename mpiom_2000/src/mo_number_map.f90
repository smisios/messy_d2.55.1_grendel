!>
!! mpiom-latest-cmip5 - mo_number_map 
!!
!! @author Karl-Hermann Wieners, Max Planck Institute for Meteorology
!!
!! $Revision$
!! $Date$
!!
module mo_number_map

    implicit none
    
    type list_item
        integer :: key, value
    end type
    
include 'mo_linked_list.inc'

    function list_item_get_id(this)
        type(list_item), intent(in) :: this
        integer :: list_item_get_id
        list_item_get_id = this%key
    end function list_item_get_id
    
  !>
  !! Format list item for output.
  !!
  !! @todo should check for overflow of output string
  !!
  FUNCTION list_item_to_string(this) RESULT(the_string)

    TYPE(list_item), INTENT(in) :: this
    CHARACTER(30) :: the_string
    
    write(the_string, "('(', I0, ', ', I0, ')')") this%key, this%value

  END FUNCTION list_item_to_string
  
  subroutine number_map_add(this, key, value)
      type(list_type), intent(inout) :: this
      integer, intent(in) :: key, value
      call list_add(this, list_item(key, value))
  end subroutine number_map_add
  
  function number_map_get(this, key, value)
      type(list_type), intent(inout) :: this
      integer, intent(in) :: key
      integer, intent(out) :: value
      logical :: number_map_get
      type(list_item) :: item
      number_map_get = list_get_item_by_id(this, key, item)
      if(number_map_get) then
          value = item%value
      end if
  end function number_map_get
  
end module mo_number_map
