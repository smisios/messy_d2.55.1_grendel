!>
!! Short description goes here, don't forget period at end.
!!
!! @todo add module description for mo_postprocess.f90
!!
module mo_postprocess

  ! Numerics
  use mo_kind, only: dp
  use mo_parallel, only: have_g_js
  use mo_param1, only: ie, je
  ! Physical system properties
  use mo_commo1, only: weto, dlxp, dlyp
  use mo_grid, only: thkcello
  ! Simulation options
  use mo_commo1, only: lbounds_exch_tp
 
  implicit none
  private
  
  public :: column_integral_field, level_integral_vector, surface_integral

contains

  pure function surface_sum(source, level)
  
     real(dp), intent(in) :: source(:,:)
     integer, intent(in), optional :: level
     real(dp) :: surface_sum
  
     integer :: i_east, i_west, j_north, j_south, halo, halo_north
     integer :: my_level
     
     if(present(level)) then
        my_level = level
     else
        my_level = 1
     end if
    
     halo = 1
     halo_north = merge(2, 1, lbounds_exch_tp .and. have_g_js)
     i_west = 1 + halo
     i_east = ie - halo
     j_north = 1 + halo_north
     j_south = je - halo
     
     surface_sum = sum(source(i_west:i_east,j_north:j_south), &
                       weto(i_west:i_east,j_north:j_south,my_level) > .5_dp)
     
  end function surface_sum
    
  function surface_integral(source)

    real(dp), intent(in) :: source(:,:)
    real(dp) :: surface_integral

    surface_integral = surface_sum(dlxp * dlyp * source)

  end function surface_integral
  
  function level_integral_vector(source)
    
     real(dp), intent(in) :: source(:,:,:)
     real(dp) :: level_integral_vector(size(source,3))
     
     integer :: k
     
     forall(k=1:size(source,3))
         level_integral_vector(k) = surface_sum(dlxp * dlyp * source(:,:,k), k)
     end forall
     
  end function level_integral_vector

  !>
  !! Compute column integral field for 3D variable given per volume.
  !!
  function column_integral_field(source)

    real(dp), intent(in) :: source(:,:,:)
    real(dp) :: column_integral_field(size(source,1),size(source,2))

    column_integral_field = &
         sum( thkcello(:,:,1:size(source,3)) * source, 3, &
              weto(:,:,1:size(source,3)) > .5_dp )

  end function column_integral_field

end module mo_postprocess
