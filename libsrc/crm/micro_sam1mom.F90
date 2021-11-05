module micro_sam1mom

! module for original SAM bulk microphysics
! Marat Khairoutdinov, 2006

use micro_params
use microphysics
use grid, only: dp

implicit none

PRIVATE :: dp
!----------------------------------------------------------------------
!!! required definitions:

!integer, parameter :: nmicro_fields = 2   ! total number of prognostic water vars

!!! microphysics prognostic variables are storred in this array:

!REAL(dp), POINTER, DIMENSION(:,:,:,:) :: micro_field     => NULL()

!integer, parameter :: flag_wmass(nmicro_fields) = (/1,1/)
!integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
!integer, parameter :: index_cloud_ice = 1   ! index for cloud ice (sedimentation)
!integer, parameter :: flag_precip(nmicro_fields) = (/0,1/)

! both variables correspond to mass, not number
!integer, parameter :: flag_number(nmicro_fields) = (/0,0/)

! SAM1MOM 3D microphysical fields are output by default.
!integer, parameter :: flag_micro3Dout(nmicro_fields) = (/0,0/)

!REAL(dp), POINTER, DIMENSION(:,:,:) :: fluxbmk  => NULL() ! surface flux of tracers
!REAL(dp), POINTER, DIMENSION(:,:,:) :: fluxtmk  => NULL() ! top boundary flux of tracers

!!! these arrays are needed for output statistics:

!REAL(dp), POINTER, DIMENSION(:,:) :: mkwle   => NULL() ! resolved vertical flux
!REAL(dp), POINTER, DIMENSION(:,:) :: mkwsb   => NULL() ! SGS vertical flux
!REAL(dp), POINTER, DIMENSION(:,:) :: mkadv   => NULL() ! tendency due to vertical advection
!REAL(dp), POINTER, DIMENSION(:,:) :: mklsadv => NULL() ! tendency due to large-scale vertical advection
!REAL(dp), POINTER, DIMENSION(:,:) :: mkdiff  => NULL() ! tendency due to vertical diffusion

!======================================================================
! UW ADDITIONS

!bloss: arrays with names/units for microphysical outputs in statistics.
!CHARACTER(LEN=3) , POINTER, DIMENSION(:) :: mkname
!CHARACTER(LEN=80), POINTER, DIMENSION(:) :: mklongname
!CHARACTER(LEN=10), POINTER, DIMENSION(:) :: mkunits
!REAL(dp), POINTER, DIMENSION(:) :: mkoutputscale

! END UW ADDITIONS
!======================================================================

!------------------------------------------------------------------
! Optional (internal) definitions)

! make aliases for prognostic variables:
! note that the aliases should be local to microphysics

!REAL(dp), POINTER, DIMENSION(:,:,:) :: q   ! total nonprecipitating water
!REAL(dp), POINTER, DIMENSION(:,:,:) :: qp  ! total precipitating water

!!!! necessary??? !!!!!
!equivalence (q(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
!equivalence (qp(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,2))
!!!!!!!!!!!!!!!!!!!!!!!

!REAL(dp), POINTER, DIMENSION(:,:,:) :: qn     ! cloud condensate (liquid + ice)
!REAL(dp), POINTER, DIMENSION(:)     :: qpsrc  ! source of precipitation microphysical processes
!REAL(dp), POINTER, DIMENSION(:)     :: qpevp  ! sink of precipitating water due to evaporation

!REAL(dp) :: vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysics options from prm file

! um_hr_20190315 -> moved micro_setparm to microphysics.F90

!subroutine micro_setparm()
!  ! no user-definable options in SAM1MOM microphysics.
!end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:

subroutine micro_init_sam1mom()

  use vars, only: q0, docloud, doprecip, nrestart, dosmoke

  implicit none

  integer k

  index_water_vapor = 1      ! index for variable that has water vapor
  index_cloud_ice   = 1      ! index for cloud ice (sedimentation)
  flag_wmass        = (/1,1/)
  flag_precip       = (/0,1/)
  flag_number       = (/0,0/)
  flag_micro3Dout   = (/0,0/)

  a_bg = 1./(tbgmax-tbgmin)
  a_pr = 1./(tprmax-tprmin)
  a_gr = 1./(tgrmax-tgrmin)

  if(nrestart.eq.0) then

     fluxbmk = 0.
     fluxtmk = 0.

     if(docloud) then
       call micro_diagnose_sam1mom()
     end if
     if(dosmoke) then
       call micro_diagnose_sam1mom()
     end if

  end if

  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.

  qpsrc = 0.
  qpevp = 0.

  mkname(1) = 'QT'
  mklongname(1) = 'TOTAL WATER (VAPOR + CONDENSATE)'
  mkunits(1) = 'g/kg'
  mkoutputscale(1) = 1.e3

  mkname(2) = 'QP'
  mklongname(2) = 'PRECIPITATING WATER'
  mkunits(2) = 'g/kg'
  mkoutputscale(2) = 1.e3

end subroutine micro_init_sam1mom

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
subroutine micro_flux_sam1mom()

  use vars, only: fluxbq, fluxtq

  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)

end subroutine micro_flux_sam1mom

!----------------------------------------------------------------------
!!! compute local microphysics processes (beyond advection and SGS diffusion):
!
subroutine micro_proc_sam1mom()

  use vars, only: nstep,dt,icycle,docloud,doprecip,dosmoke

  ! Update bulk coefficient
  if(doprecip.and.icycle.eq.1) call precip_init() 
  
  if(docloud) then
     call cloud_crm(micro_field(:,:,:,1),micro_field(:,:,:,2))
     if(doprecip) call precip_proc(micro_field(:,:,:,1),micro_field(:,:,:,2))
     call micro_diagnose_sam1mom()
  end if
  if(dosmoke) then
     call micro_diagnose_sam1mom()
  end if
  
end subroutine micro_proc_sam1mom

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine micro_diagnose_sam1mom()
 
  use vars

  real omn, omp
  integer i,j,k

  do k=1,nzm
     do j=1,ny
        do i=1,nx
           qv(i,j,k)  = micro_field(i,j,k,1) - qn(i,j,k)
           omn        = max(0._dp,min(1._dp,(tabs(i,j,k)-tbgmin)*a_bg))
           qcl(i,j,k) = qn(i,j,k)*omn
           qci(i,j,k) = qn(i,j,k)*(1.-omn)
           omp        = max(0._dp,min(1._dp,(tabs(i,j,k)-tprmin)*a_pr))
           qpl(i,j,k) = micro_field(i,j,k,2)*omp
           qpi(i,j,k) = micro_field(i,j,k,2)*(1.-omp)
        end do
     end do
  end do

end subroutine micro_diagnose_sam1mom

!----------------------------------------------------------------------
!!! function to compute terminal velocity for precipitating variables:
! In this particular case there is only one precipitating variable.

real function term_vel_qp(i,j,k,ind)
  
  use vars
  integer, intent(in) :: i,j,k,ind
  real wmax, omp, omg, qrr, qss, qgg

  term_vel_qp = 0.
  if(micro_field(i,j,k,2).gt.qp_threshold) then
    omp = max(0._dp,min(1._dp,(tabs(i,j,k)-tprmin)*a_pr))
    if(omp.eq.1.) then
       term_vel_qp = vrain*(rho(k)*micro_field(i,j,k,2))**crain
    elseif(omp.eq.0.) then
       omg = max(0._dp,min(1._dp,(tabs(i,j,k)-tgrmin)*a_gr))
       qgg=omg*micro_field(i,j,k,2)
       qss=micro_field(i,j,k,2)-qgg
       term_vel_qp = (omg*vgrau*(rho(k)*qgg)**cgrau &
                                 +(1.-omg)*vsnow*(rho(k)*qss)**csnow)
    else
       omg = max(0._dp,min(1._dp,(tabs(i,j,k)-tgrmin)*a_gr))
       qrr=omp*micro_field(i,j,k,2)
       qss=micro_field(i,j,k,2)-qrr
       qgg=omg*qss
       qss=qss-qgg
       term_vel_qp = (omp*vrain*(rho(k)*qrr)**crain &
                     +(1.-omp)*(omg*vgrau*(rho(k)*qgg)**cgrau &
                          +(1.-omg)*vsnow*(rho(k)*qss)**csnow))
    endif
  end if  
end function term_vel_qp

!----------------------------------------------------------------------
!!! compute sedimentation 
!
subroutine micro_precip_fall_sam1mom()
  
  use vars
  use params, only : pi

  real omega(nx,ny,nzm)
  integer ind
  integer i,j,k

  crain = b_rain / 4.
  csnow = b_snow / 4.
  cgrau = b_grau / 4.
  vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain
  vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
  vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

 do k=1,nzm
  do j=1,ny
   do i=1,nx
       omega(i,j,k) = max(0._dp,min(1._dp,(tabs(i,j,k)-tprmin)*a_pr))
   end do
  end do
 end do

 call precip_fall(micro_field(:,:,:,2), term_vel_qp, 2, omega, ind)

 do j=1,ny
   do i=1,nx
     if(micro_field(i,j,1,2).gt.1.e-6) s_ar=s_ar+dtfactor
   end do
 end do


end subroutine micro_precip_fall_sam1mom

!----------------------------------------------------------------------
end module micro_sam1mom
!----------------------------------------------------------------------
