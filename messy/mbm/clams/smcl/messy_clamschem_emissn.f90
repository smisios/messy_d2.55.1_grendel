Module messy_clamschem_emissn


contains


!     emissn   - Emissions scheme
!    use parameterization by Heaps (Planet. Space Sci., 26, 513-517, 1978)
!    to get ion pair production rate by cosmic rays
!    and efficiency 1.5 NO molecules produced per ionization 
!    (after Crutzen et al., Science 189, 457-459, 1975)

!
!     written by V. Cals 
!     last modification: 3-Jul-2001 (jug)


subroutine emissn
  USE messy_clams_global,         ONLY: prec
  !use messy_clamschem_phys_const, only: torad
  use messy_clamschem_asad_mod,   only: speci, nlemit, nemit, emr, tnd
  use messy_clamschem_global,     only: ntraj, torad, &
                                        lats, lons, press
  implicit none

  real(PREC) :: coslambda, lambda, n
  integer       :: i, j
  
  ! Koordinaten des magnetischen Nordpols (2000)
  real(PREC), parameter::  lat0 = 81.65 * torad,  lon0 = 277.36 * torad

  ! value for cosmic ray maximum (solar cycle minimum) 
  real(PREC), parameter::   A = 1.74E-18 ,   B = 2.84E-17
  ! alternative for cosmic ray minimum (solar cycle maximum) 
  !real(PREC), parameter::  A = 1.74E-18 ,   B = 1.93E-17

  ! efficiency NO molecules produced per ionization 
  !( after Crutzen et al., Science 189, 457-459, 1975)
  real(PREC), parameter::   eff=1.5

  do i=1, ntraj
     coslambda = cos(lat0)*cos(lats(i)*torad)*cos(lons(i)*torad-lon0) + &
          sin(lat0)*sin(lats(i)*torad)
     ! coslamda =  cos(Winkel) im Bogenmass zw magn NP und Bezugspunkt
     
     !write(*,*)'coslambda=',coslambda
     lambda = asin(coslambda)
     !write(*,*)'lat,lon,lambda=',lats(i),lons(i),lambda/torad
     do j=1, nemit
        emr(i,nlemit(j)) = 0
        select case (trim(speci(nlemit(j))))
        case('NO')
           if( press(i) .lt. 1000. ) then
              ! zweite Formel  > 30km (p <= 1000 Pa)
              n = 0.6+0.8*cos(lambda)
              emr(i,nlemit(j)) = eff*(A+(B*(sin(lambda))**4))*&
                   (3.03E+17**(1-n))*(tnd(i)**n)
           else if( press(i) .lt. 7000. ) then
              ! erste Formel  18-30km  (p <= 7000 Pa)
              emr(i,nlemit(j)) = eff*(A+(B*(sin(lambda))**4))*tnd(i)
           end if
        end select
     end do
  end do

  return
end subroutine emissn



end Module messy_clamschem_emissn
