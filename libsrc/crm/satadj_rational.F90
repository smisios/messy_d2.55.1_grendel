real function satadj_water(t_in,qv_in,pref)

  ! this function finds the (possibly new) temperature which results from
  !   saturation adjustment with respect to water saturation.  The
  !   conservation of moist static energy during the saturation
  !   adjustment process is the basis for the scheme.  The moist
  !   static energy is evaluated without the potential energy terms as:
  !
  !     h(T)/Cp = T + (L/Cp)*qv = T + (L/Cp)*min(qsat(T,pref),qv_in)
  !
  ! The adjustment is formulated as a root finding problem for the
  ! function
  !
  !     f(T) = h(T) - h(T_in)
  !
  ! which is solved using a combined secant/bisection method.
  !
  ! written by Peter Blossey, January 2007.

!!$  use params, only: fac_cond

  implicit none

  real, parameter :: fac_cond = 2.5104e+06/1004.

  real, intent(in) :: t_in, qv_in, pref

  real :: t0, t1, t2, t3, tpos, tneg ! guesses for temperature at saturation
  real :: r0, r1, r2, r3, rpos, rneg ! residuals associated with those guesses

  real :: satratio, alpha, beta, qsw_in
  integer :: count

  real, external :: qsatw_crm

  qsw_in = qsatw_crm(t_in,pref)
  satratio = qv_in/qsw_in 
!!$  write(*,*) 'satratio = ', satratio

  if(satratio.le.1.) then
     ! if input water vapor qv_in is less than saturation value, return
     !   with temperature at saturation as output.
     satadj_water = t_in
     return
  end if

  ! if we are supersaturated, we need to do saturation adjustment.

  ! generate two initial guesses for the temperature
  t0 = t_in
  t1 = min(343.15,t_in + fac_cond*(satratio-1.)*qv_in)

!!$  write(*,*) 'initial t guesses = ', t0, t1

  ! find the value/residual of the function f(T) that we want a root of.
  r0 = t0-t_in + fac_cond*min(0.,qsw_in-qv_in)
  r1 = t1-t_in + fac_cond*min(0.,qsatw_crm(t1,pref)-qv_in)

!!$  write(*,*) 'initial residuals = ', t0, t1

  tneg = t0 ! tneg should have rneg < 0
  tpos = t1 ! tpos should have rpos < 0

  rneg = r0
  rpos = r1

  ! new guess from secant
  t2 = t1 - r1*(t0-t1)/(r0-r1)

  ! if new guess is outside of [tneg,tpos], use bisection
  if((t2.gt.tpos).or.(t2.lt.tneg)) t2 = 0.5*(tneg+tpos)

  ! compute residual
  r2 = t2-t_in + fac_cond*min(0.,qsatw_crm(t2,pref)-qv_in)

!!$  write(*,*) 'temperatures = ', t2, t1, t0
!!$  write(*,*) 'residuals = ', r2, r1, r0 

  if(r2.ge.0.) then
     tpos = t2 ! update upper bound for bisection
     rpos = r2
  else
     tneg = t2 ! update lower bound for bisection
     rneg = r2
  end if

  do count = 4,100 ! count is number of function evaluations

     ! use combination bisection/secant/rational interpolation
     !   method to find root of f(T).  The rational interpolation
     !   is taken from Bus & Dekker (1975).

     if(rneg*rpos>0.) then
        write(*,*) 'Error in saturation adjustment in micro_init.'
        write(*,*) 'Need two temperature guesses that bracket saturation.'
        write(*,*) 'T guesses = ', tneg, tpos
        write(*,*) 'f(T guess)= ', rneg, rpos
        EXIT
     end if

     ! new guess from rational function interpolation 
     !   (Bus & Dekker ACM Trans on Math Software, vol. 1, p. 334, 1975.
     !    See equation 3.1.5.)
     alpha = r2*(r1-r0)/(t1-t0)
     beta = r1*(r2-r0)/(t2-t0)
     t3 = t1 - beta*(t1-t2)/(beta-alpha)

     ! if new guess is outside of [tneg,tpos], try secant and then bisection
     if((t3.gt.tpos).or.(t3.lt.tneg)) then
        t3 = t2 - r2*(t1-t2)/(r1-r2) ! secant
        
        if((t3.gt.tpos).or.(t3.lt.tneg)) then ! if secant is out of bounds...
           t3 = 0.5*(tneg+tpos) !bisection
        end if
     end if

     ! compute residual
     r3 = t3-t_in + fac_cond*min(0.,qsatw_crm(t3,pref)-qv_in)

!!$     write(*,*) 'count = ', count
!!$     write(*,*) 'temperatures = ', t3, t2, t1
!!$     write(*,*) 'residuals = ', r3, r2, r1 
!!$     write(*,*)

     if((abs(r3).lt.4.e-7*t_in).or.(abs(t3-t2).lt.4.e-7*t_in)) then
        t2 = t3
!        write(*,*) 'Required # of function evals for sat adj = ', count
        EXIT  ! drop out of loop and exit
     end if

     if(r3.ge.0.) then
        tpos = t3 ! update upper bound for bisection
        rpos = r3
     else
        tneg = t3 ! update lower bound for bisection
        rneg = r3
     end if

     t0 = t1 ! update oldest guess
     r0 = r1

     t1 = t2 ! update old guess 
     r1 = r2

     t2 = t3 ! update newer guess 
     r2 = r3

     if(count.gt.100) then
        write(*,*) &
             'Warning: Exceeded 100 iterations in satadj_water w/o converging'
        EXIT
     end if

  end do

  satadj_water = t2

end function satadj_water
