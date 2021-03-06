module messy_vertdiff_camtrbmtnstress

  USE messy_main_constants_mem, ONLY: r8=>dp, &
                                      karman => c_vKar, &      ! von Karman constant
                                      gravit => g , &          ! Acceleration due to gravity
                                      rair   => rd             ! Gas constant for dry air

  implicit none
  private      
  save

  public compute_tms                          ! Full routine

  ! ------------ !
  ! Private data !
  ! ------------ !

  real(r8), parameter :: horomin= 1._r8       ! Minimum value of subgrid orographic height for mountain stress [ m ]
  real(r8), parameter :: z0max  = 100._r8     ! Maximum value of z_0 for orography [ m ]
  real(r8), parameter :: dv2min = 0.01_r8     ! Minimum shear squared [ m2/s2 ]
  real(r8), public    :: orocnst =  1.0D0     ! Converts from standard deviation to height [ no unit ]
  real(r8), public    :: z0fac   = 0.075D0    ! Factor determining z_0 from orographic standard deviation [ no unit ] 

contains

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine compute_tms( pcols    , pver    , ncol    ,                     &
                          u        , v       , t       , pmid    , exner   , &
                          zm       , sgh     , ksrf    , taux    , tauy    , & 
                          landfrac )

    !------------------------------------------------------------------------------ !
    ! Turbulent mountain stress parameterization                                    !  
    !                                                                               !
    ! Returns surface drag coefficient and stress associated with subgrid mountains !
    ! For points where the orographic variance is small ( including ocean ),        !
    ! the returned surface drag coefficient and stress is zero.                     !
    !                                                                               !
    ! Lastly arranged : Sungsu Park. Jan. 2010.                                     !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments ! 
    ! ---------------------- !

    integer,  intent(in)  :: pcols                 ! Number of columns dimensioned
    integer,  intent(in)  :: pver                  ! Number of model layers
    integer,  intent(in)  :: ncol                  ! Number of columns actually used

    real(r8), intent(in)  :: u(:,:)         ! Layer mid-point zonal wind [ m/s ]
    real(r8), intent(in)  :: v(:,:)         ! Layer mid-point meridional wind [ m/s ]
    real(r8), intent(in)  :: t(:,:)         ! Layer mid-point temperature [ K ]
    real(r8), intent(in)  :: pmid(:,:)      ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)  :: exner(:,:)     ! Layer mid-point exner function [ no unit ]
    real(r8), intent(in)  :: zm(:,:)        ! Layer mid-point height [ m ]
    real(r8), intent(in)  :: sgh(:)            ! Standard deviation of orography [ m ]
    real(r8), intent(in)  :: landfrac(:)       ! Land fraction [ fraction ]
    
    real(r8), intent(out) :: ksrf(:)           ! Surface drag coefficient [ kg/s/m2 ]
    real(r8), intent(out) :: taux(:)           ! Surface zonal      wind stress [ N/m2 ]
    real(r8), intent(out) :: tauy(:)           ! Surface meridional wind stress [ N/m2 ]


    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer  :: i                                  ! Loop index
    integer  :: kb, kt                             ! Bottom and top of source region
    
    real(r8) :: horo                               ! Orographic height [ m ]
    real(r8) :: z0oro                              ! Orographic z0 for momentum [ m ]
    real(r8) :: dv2                                ! (delta v)**2 [ m2/s2 ]
    real(r8) :: ri                                 ! Richardson number [ no unit ]
    real(r8) :: stabfri                            ! Instability function of Richardson number [ no unit ]
    real(r8) :: rho                                ! Density [ kg/m3 ]
    real(r8) :: cd                                 ! Drag coefficient [ no unit ]
    real(r8) :: vmag                               ! Velocity magnitude [ m /s ]

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
       
    do i = 1, ncol

     ! determine subgrid orgraphic height ( mean to peak )

       horo = orocnst * sgh(i)

     ! No mountain stress if horo is too small

       if( horo < horomin ) then

           ksrf(i) = 0._r8
           taux(i) = 0._r8
           tauy(i) = 0._r8

       else

         ! Determine z0m for orography

           z0oro = min( z0fac * horo, z0max )

         ! Calculate neutral drag coefficient

           cd = ( karman / log( ( zm(i,pver) + z0oro ) / z0oro) )**2

         ! Calculate the Richardson number over the lowest 2 layers

           kt  = pver - 1
           kb  = pver
           dv2 = max( ( u(i,kt) - u(i,kb) )**2 + ( v(i,kt) - v(i,kb) )**2, dv2min )

         ! Modification : Below computation of Ri is wrong. Note that 'Exner' function here is
         !                inverse exner function. Here, exner function is not multiplied in
         !                the denominator. Also, we should use moist Ri not dry Ri.
         !                Also, this approach using the two lowest model layers can be potentially
         !                sensitive to the vertical resolution.  
         ! OK. I only modified the part associated with exner function.

           ri  = 2._r8 * gravit * ( t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb) ) * ( zm(i,kt) - zm(i,kb) ) &
                                / ( ( t(i,kt) * exner(i,kt) + t(i,kb) * exner(i,kb) ) * dv2 )

         ! ri  = 2._r8 * gravit * ( t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb) ) * ( zm(i,kt) - zm(i,kb) ) &
         !                      / ( ( t(i,kt) + t(i,kb) ) * dv2 )

         ! Calculate the instability function and modify the neutral drag cofficient.
         ! We should probably follow more elegant approach like Louis et al (1982) or Bretherton and Park (2009) 
         ! but for now we use very crude approach : just 1 for ri < 0, 0 for ri > 1, and linear ramping.

           stabfri = max( 0._r8, min( 1._r8, 1._r8 - ri ) )
           cd      = cd * stabfri

         ! Compute density, velocity magnitude and stress using bottom level properties

           rho     = pmid(i,pver) / ( rair * t(i,pver) ) 
           vmag    = sqrt( u(i,pver)**2 + v(i,pver)**2 )
           ksrf(i) = rho * cd * vmag * landfrac(i)
           taux(i) = -ksrf(i) * u(i,pver)
           tauy(i) = -ksrf(i) * v(i,pver)

       end if

    end do
    
    return
  end subroutine compute_tms

end module messy_vertdiff_camtrbmtnstress
