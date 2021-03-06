module messy_vertdiff_tools

!--------------------------------------------------------------------!
! Module Overview:                                                   !
!                                                                    !
! This module provides an interface to wv_sat_methods, providing     !
! saturation vapor pressure and related calculations to CAM.         !
!                                                                    !
! The original wv_saturation codes were introduced by J. J. Hack,    !
! February 1990. The code has been extensively rewritten since then, !
! including a total refactoring in Summer 2012.                      !
!                                                                    !
!--------------------------------------------------------------------!
! Methods:                                                           !
!                                                                    !
! Pure water/ice saturation vapor pressures are calculated on the    !
! fly, with the specific method determined by a runtime option.      !
! Mixed phase SVP is interpolated from the internal table, estbl,    !
! which is created during initialization.                            !
!                                                                    !
! The default method for calculating SVP is determined by a namelist !
! option, and used whenever svp_water/ice or qsat are called.        !
!                                                                    !
!--------------------------------------------------------------------!

  USE messy_main_constants_mem, ONLY: r8=>dp, &
                               cpair  => cp_air , &     ! Specific heat of dry air
                               gravit => g , &     ! Acceleration due to gravity
                               g , &
                               rair   => rd , &     ! Gas constant for dry air
                               zvir   => vtmpc1  , &     ! rh2o/rair - 1
                               latvap => alv , &     ! Latent heat of vaporization
                               latice => alf, &     ! Latent heat of fusion
                               vk => c_vKar, &     ! von Karman constant
                               mwdry  => M_air , &     ! Molecular weight of dry air
                               avogad => N_A_kmol, &     ! Avogadro's number
!!$                               boltz  => k_B, &      ! Boltzman's constant
                               epsilo => MM_eps, &
                               rh2o => rv
                               !tms_orocnst,&  ! turbulent mountain stress parameter
                               !tms_z0fac      ! Factor determining z_0 from orographic standard deviation [no unit]
!!$use shr_kind_mod, only: r8 => shr_kind_r8
!!$use physconst,    only: epsilo, &
!!$                        latvap, &
!!$                        latice, &
!!$                        rh2o,   &
!!$                        cpair,  &
!!$                        tmelt,  &
!!$                        h2otrip

!!$use wv_sat_methods, only: &
!!$     svp_to_qsat => wv_sat_svp_to_qsat, &
!!$     svp_to_qmmr => wv_sat_svp_to_qmmr

implicit none
private
save

! Public interfaces
! Namelist, initialization, finalization
public wv_sat_readnl
public wv_sat_init
public wv_sat_final

! Saturation vapor pressure calculations
public svp_water
public svp_ice
  
! Mixed phase (water + ice) saturation vapor pressure table lookup
public estblf

public wv_sat_svp_to_qsat
public wv_sat_svp_to_qmmr

! Subroutines that return both SVP and humidity
! Optional arguments do temperature derivatives
public qsat           ! Mixed phase
public qsat_water     ! SVP over water only
public qsat_ice       ! SVP over ice only

! Wet bulb temperature solver
public findsp_vc

! Unusual/non-default methods.
! Zhang-McFarlane scheme uses mass mixing ratio rather than
! specific humidity.
public qmmr


! from pbl_utils:
public virtem
public calc_ustar
public calc_obklen

! Data

! This value is slightly high, but it seems to be the value for the
! steam point of water originally (and most frequently) used in the
! Goff & Gratch scheme.
real(r8), parameter :: tboil = 373.16_r8

! Table of saturation vapor pressure values (estbl) from tmin to
! tmax+1 Kelvin, in one degree increments.  ttrice defines the
! transition region, estbl contains a combination of ice & water
! values.
! Make these public parameters in case another module wants to see the
! extent of the table.
  real(r8), public, parameter :: tmin = 127.16_r8
  real(r8), public, parameter :: tmax = 375.16_r8

  real(r8), parameter :: ttrice = 20.00_r8  ! transition range from es over H2O to es over ice

  integer :: plenest                             ! length of estbl
! op_pj_20140401+
!!$  real(r8), save, allocatable :: estbl(:)              ! table values of saturation vapor pressure
  real(r8), allocatable :: estbl(:)              ! table values of saturation vapor pressure
! op_pj_20140401-

  real(r8),PARAMETER :: omeps=1.0_r8-epsilo      ! 1.0_r8 - epsilo

  real(r8) :: c3         ! parameter used by findsp

  ! Set coefficients for polynomial approximation of difference
  ! between saturation vapor press over water and saturation pressure
  ! over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial is
  ! valid in the range -40 < t < 0 (degrees C).
  real(r8) :: pcf(5) = (/ &
       5.04469588506e-01_r8, &
       -5.47288442819e+00_r8, &
       -3.67471858735e-01_r8, &
       -8.95963532403e-03_r8, &
       -7.78053686625e-05_r8 /)

!   --- Degree 6 approximation ---
!  real(r8) :: pcf(6) = (/ &
!       7.63285250063e-02, &
!       5.86048427932e+00, &
!       4.38660831780e-01, &
!       1.37898276415e-02, &
!       2.14444472424e-04, &
!       1.36639103771e-06 /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! original CAM code: wv_sat_methods.F90

real(r8) :: tmelt   ! Melting point of water at 1 atm (K)
real(r8) :: h2otrip ! Triple point temperature of water (K)
!!$real(r8) :: tboil   ! Boiling point of water at 1 atm (K)

!!$real(r8) :: ttrice  ! Ice-water transition range

!!$real(r8) :: epsilo  ! Ice-water transition range
!!$real(r8) :: omeps   ! 1._r8 - epsilo

! Indices representing individual schemes
integer, parameter :: Invalid_idx = -1
integer, parameter :: OldGoffGratch_idx = 0
integer, parameter :: GoffGratch_idx = 1
integer, parameter :: MurphyKoop_idx = 2
integer, parameter :: Bolton_idx = 3

! Index representing the current default scheme.
integer, parameter :: initial_default_idx = GoffGratch_idx
integer :: default_idx = initial_default_idx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! original CAM code: pbl_utils

real(r8), parameter :: ustar_min = 0.01_r8

!!$real(r8) :: g         ! acceleration of gravity
!!$real(r8) :: vk        ! Von Karman's constant
!!$real(r8) :: cpair     ! specific heat of dry air
!!$real(r8) :: rair      ! gas constant for dry air
!!$real(r8) :: zvir      ! rh2o/rair - 1


contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

subroutine wv_sat_readnl(nlfile)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Get runtime options for wv_saturation.                         !
  !------------------------------------------------------------------!

!!$  use wv_sat_methods, only: wv_sat_get_scheme_idx, &
!!$                            wv_sat_valid_idx, &
!!$                            wv_sat_set_default

!!$  use spmd_utils,      only: masterproc
!!$  use namelist_utils,  only: find_group_name
!!$  use units,           only: getunit, freeunit
!!$  use mpishorthand
!!$  use abortutils,      only: endrun

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
   
  ! Local variables
  integer :: unitn, ierr

  character(len=32) :: wv_sat_scheme = "GoffGratch"

  character(len=*), parameter :: subname = 'wv_sat_readnl'

!!$  namelist /wv_sat_nl/ wv_sat_scheme
  !-----------------------------------------------------------------------------

!!$  if (masterproc) then
!!$     unitn = getunit()
!!$     open( unitn, file=trim(nlfile), status='old' )
!!$     call find_group_name(unitn, 'wv_sat_nl', status=ierr)
!!$     if (ierr == 0) then
!!$        read(unitn, wv_sat_nl, iostat=ierr)
!!$        if (ierr /= 0) then
!!$           call endrun(subname // ':: ERROR reading namelist')
!!$           return
!!$        end if
!!$     end if
!!$     close(unitn)
!!$     call freeunit(unitn)
!!$
!!$  end if
!!$
!!$#ifdef SPMD
!!$  call mpibcast(wv_sat_scheme, len(wv_sat_scheme) , mpichar, 0, mpicom)
!!$#endif
!!$
!!$  if (.not. wv_sat_set_default(wv_sat_scheme)) then
!!$     call endrun('wv_sat_readnl :: Invalid wv_sat_scheme.')
!!$     return
!!$  end if

end subroutine wv_sat_readnl

subroutine wv_sat_init(status)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Initialize module (e.g. setting parameters, initializing the   !
  !   SVP lookup table).                                             !
  !------------------------------------------------------------------!

!!$  use wv_sat_methods, only: wv_sat_methods_init, &
!!$                            wv_sat_get_scheme_idx, &
!!$                            wv_sat_valid_idx
!!$  use spmd_utils,     only: masterproc
!!$  use cam_logfile,    only: iulog
!!$  use abortutils,     only: endrun
!!$  use shr_assert_mod, only: shr_assert_in_domain
!!$  use error_messages, only: handle_errmsg

  integer, intent(out) :: status

  ! For wv_sat_methods error reporting.
  character(len=256) :: errstring

  ! For generating internal SVP table.
  real(r8) :: t         ! Temperature
  integer  :: i         ! Increment counter

  status = 0

!!$  ! Precalculated because so frequently used.
!!$  omeps  = 1.0_r8 - epsilo

  ! Transition range method is only valid for transition temperatures at:
  ! -40 deg C < T < 0 deg C
!!$  call shr_assert_in_domain(ttrice, ge=0._r8, le=40._r8, varname="ttrice",&
!!$       msg="wv_sat_init: Invalid transition temperature range.")

! This parameter uses a hardcoded 287.04_r8?
  c3 = 287.04_r8*(7.5_r8*log(10._r8))/cpair

! Init "methods" module containing actual SVP formulae.

!!$  call wv_sat_methods_init(r8, tmelt, h2otrip, tboil, ttrice, &
!!$       epsilo, errstring)
!!$
!!$  call handle_errmsg(errstring, subname="wv_sat_methods_init")

  ! Add two to make the table slightly too big, just in case.
  plenest = ceiling(tmax-tmin) + 2

  ! Allocate SVP table.
  allocate(estbl(plenest), stat=status)
  if (status /= 0) then 
!!$     call endrun('wv_sat_init :: ERROR allocating saturation vapor pressure table')
     return
  end if

  do i = 1, plenest
    estbl(i) = svp_trans(tmin + real(i-1,r8))
  end do

!!$  if (masterproc) then
     write(*,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
!!$  end if

end subroutine wv_sat_init

subroutine wv_sat_final
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Deallocate global variables in module.                         !
  !------------------------------------------------------------------!
!!$  use abortutils,   only: endrun

  integer :: status

  if (allocated(estbl)) then

     deallocate(estbl, stat=status)

     if (status /= 0) then
!!$        call endrun('wv_sat_final :: ERROR deallocating table')
        return
     end if

  end if

end subroutine wv_sat_final

!---------------------------------------------------------------------
! DEFAULT SVP FUNCTIONS
!---------------------------------------------------------------------

! Compute saturation vapor pressure over water
elemental function svp_water(t) result(es)

!!$  use wv_sat_methods, only: &
!!$       wv_sat_svp_water

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = wv_sat_svp_water(T)

end function svp_water

! Compute saturation vapor pressure over ice
elemental function svp_ice(t) result(es)

!!$  use wv_sat_methods, only: &
!!$       wv_sat_svp_ice

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = wv_sat_svp_ice(T)

end function svp_ice

! Compute saturation vapor pressure with an ice-water transition
elemental function svp_trans(t) result(es)

!!$  use wv_sat_methods, only: &
!!$       wv_sat_svp_trans

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = wv_sat_svp_trans(T)

end function svp_trans

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Does linear interpolation from nearest values found
! in the table (estbl).
elemental function estblf(t) result(es)

  real(r8), intent(in) :: t ! Temperature 
  real(r8) :: es            ! SVP (Pa)

  integer  :: i         ! Index for t in the table
  real(r8) :: t_tmp     ! intermediate temperature for es look-up

  real(r8) :: weight ! Weight for interpolation

  t_tmp = max(min(t,tmax)-tmin, 0._r8)   ! Number of table entries above tmin
  i = int(t_tmp) + 1                     ! Corresponding index.
  weight = t_tmp - aint(t_tmp, r8)       ! Fractional part of t_tmp (for interpolation).
  es = (1._r8 - weight)*estbl(i) + weight*estbl(i+1)

end function estblf

! Get enthalpy based only on temperature
! and specific humidity.
elemental function tq_enthalpy(t, q, hltalt) result(enthalpy)

  real(r8), intent(in) :: t      ! Temperature
  real(r8), intent(in) :: q      ! Specific humidity
  real(r8), intent(in) :: hltalt ! Modified hlat for T derivatives

  real(r8) :: enthalpy

  enthalpy = cpair * t + hltalt * q
  
end function tq_enthalpy

!---------------------------------------------------------------------
! LATENT HEAT OF VAPORIZATION CORRECTIONS
!---------------------------------------------------------------------

elemental subroutine no_ip_hltalt(t, hltalt)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate latent heat of vaporization of pure liquid water at  !
  !   a given temperature.                                           !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t        ! Temperature
  ! Outputs
  real(r8), intent(out) :: hltalt  ! Appropriately modified hlat

  hltalt = latvap

  ! Account for change of latvap with t above freezing where
  ! constant slope is given by -2369 j/(kg c) = cpv - cw
  if (t >= tmelt) then
     hltalt = hltalt - 2369.0_r8*(t-tmelt)
  end if

end subroutine no_ip_hltalt

elemental subroutine calc_hltalt(t, hltalt, tterm)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate latent heat of vaporization of water at a given      !
  !   temperature, taking into account the ice phase if temperature  !
  !   is below freezing.                                             !
  !   Optional argument also calculates a term used to calculate     !
  !   d(es)/dT within the water-ice transition range.                !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t        ! Temperature
  ! Outputs
  real(r8), intent(out) :: hltalt  ! Appropriately modified hlat
  ! Term to account for d(es)/dT in transition region.
  real(r8), intent(out), optional :: tterm

  ! Local variables
  real(r8) :: tc      ! Temperature in degrees C
  real(r8) :: weight  ! Weight for es transition from water to ice
  ! Loop iterator
  integer :: i

  if (present(tterm)) tterm = 0.0_r8

  call no_ip_hltalt(t,hltalt)
  if (t < tmelt) then
     ! Weighting of hlat accounts for transition from water to ice.
     tc = t - tmelt

     if (tc >= -ttrice) then
        weight = -tc/ttrice

        ! polynomial expression approximates difference between es
        ! over water and es over ice from 0 to -ttrice (C) (max of
        ! ttrice is 40): required for accurate estimate of es
        ! derivative in transition range from ice to water
        if (present(tterm)) then
           do i = size(pcf), 1, -1
              tterm = pcf(i) + tc*tterm
           end do
           tterm = tterm/ttrice
        end if

     else
        weight = 1.0_r8
     end if

     hltalt = hltalt + weight*latice

  end if

end subroutine calc_hltalt

!---------------------------------------------------------------------
! OPTIONAL OUTPUTS
!---------------------------------------------------------------------

! Temperature derivative outputs, for qsat_*
elemental subroutine deriv_outputs(t, p, es, qs, hltalt, tterm, &
     gam, dqsdt)

  ! Inputs
  real(r8), intent(in) :: t      ! Temperature
  real(r8), intent(in) :: p      ! Pressure
  real(r8), intent(in) :: es     ! Saturation vapor pressure
  real(r8), intent(in) :: qs     ! Saturation specific humidity
  real(r8), intent(in) :: hltalt ! Modified latent heat
  real(r8), intent(in) :: tterm  ! Extra term for d(es)/dT in
                                 ! transition region.

  ! Outputs
  real(r8), intent(out), optional :: gam      ! (hltalt/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt    ! (d(qs)/dt)

  ! Local variables
  real(r8) :: desdt        ! d(es)/dt
  real(r8) :: dqsdt_loc    ! local copy of dqsdt

  if (qs == 1.0_r8) then
     dqsdt_loc = 0._r8
  else
     desdt = hltalt*es/(rh2o*t*t) + tterm
     dqsdt_loc = qs*p*desdt/(es*(p-omeps*es))
  end if

  if (present(dqsdt)) dqsdt = dqsdt_loc
  if (present(gam))   gam   = dqsdt_loc * (hltalt/cpair)

end subroutine deriv_outputs

!---------------------------------------------------------------------
! QSAT (SPECIFIC HUMIDITY) PROCEDURES
!---------------------------------------------------------------------

elemental subroutine qmmr(t, p, es, qm)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !     Provide saturation mass mixing ratio.                        !
  !                                                                  !
  ! Note that qmmr is a ratio over dry air, whereas qsat is a ratio  !
  ! over total mass. These differ mainly at low pressures, where     !
  ! SVP may exceed the actual pressure, and therefore qmmr can blow  !
  ! up.                                                              !
  !------------------------------------------------------------------!
!!$  use wv_sat_methods, only: &
!!$       wv_sat_svp_water


  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass)

  es = wv_sat_svp_water(t)

  qm = wv_sat_svp_to_qmmr(es, p)

  ! Ensures returned es is consistent with limiters on qmmr.
  es = min(es, p)

end subroutine qmmr

elemental subroutine qsat(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Look up and return saturation vapor pressure from precomputed  !
  !   table, then calculate and return saturation specific humidity. !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives
  real(r8) :: tterm        ! Account for d(es)/dT in transition region

  es = estblf(t)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

  ! Calculate optional arguments.
  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! "generalized" analytic expression for t derivative of es
     ! accurate to within 1 percent for 173.16 < t < 373.16
     call calc_hltalt(t, hltalt, tterm)

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     call deriv_outputs(t, p, es, qs, hltalt, tterm, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat

elemental subroutine qsat_water(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

!!$  use wv_sat_methods, only: wv_sat_qsat_water

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives

  call wv_sat_qsat_water(t, p, es, qs)

  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! "generalized" analytic expression for t derivative of es
     ! accurate to within 1 percent for 173.16 < t < 373.16
     call no_ip_hltalt(t, hltalt)

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     ! For pure water/ice transition term is 0.
     call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat_water

elemental subroutine qsat_ice(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

!!$  use wv_sat_methods, only: wv_sat_qsat_ice

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives

  call wv_sat_qsat_ice(t, p, es, qs)

  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! For pure ice, just add latent heats.
     hltalt = latvap + latice

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     ! For pure water/ice transition term is 0.
     call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat_ice

!---------------------------------------------------------------------
! FINDSP (WET BULB TEMPERATURE) PROCEDURES
!---------------------------------------------------------------------

subroutine findsp_vc(q, t, p, use_ice, tsp, qsp)

!!$  use cam_logfile,  only: iulog
!!$  use abortutils,   only: endrun
  ! Wrapper for findsp which is 1D and handles the output status.
  ! Changing findsp to elemental restricted debugging output.
  ! If that output is needed again, it's preferable *not* to copy findsp,
  ! but to change the existing version.

  ! input arguments
  real(r8), intent(in) :: q(:)        ! water vapor (kg/kg)
  real(r8), intent(in) :: t(:)        ! temperature (K)
  real(r8), intent(in) :: p(:)        ! pressure    (Pa)
  logical,  intent(in) :: use_ice     ! flag to include ice phase in calculations

  ! output arguments
  real(r8), intent(out) :: tsp(:)     ! saturation temp (K)
  real(r8), intent(out) :: qsp(:)     ! saturation mixing ratio (kg/kg)

  integer :: status(size(q))   ! flag representing state of output
                               ! 0 => Successful convergence
                               ! 1 => No calculation done: pressure or specific
                               !      humidity not within usable range
                               ! 2 => Run failed to converge
                               ! 4 => Temperature fell below minimum
                               ! 8 => Enthalpy not conserved

  integer :: n, i

  n = size(q)

  call findsp(q, t, p, use_ice, tsp, qsp, status)

  ! Currently, only 2 and 8 seem to be treated as fatal errors.
  do i = 1,n
     if (status(i) == 2) then
        write(*,*) ' findsp not converging at i = ', i
        write(*,*) ' t, q, p ', t(i), q(i), p(i)
        write(*,*) ' tsp, qsp ', tsp(i), qsp(i)
!!$        call endrun ('wv_saturation::FINDSP -- not converging')
     else if (status(i) == 8) then
        write(*,*) ' the enthalpy is not conserved at i = ', i
        write(*,*) ' t, q, p ', t(i), q(i), p(i)
        write(*,*) ' tsp, qsp ', tsp(i), qsp(i)
!!$        call endrun ('wv_saturation::FINDSP -- enthalpy is not conserved')
     endif
  end do

end subroutine findsp_vc

elemental subroutine findsp (q, t, p, use_ice, tsp, qsp, status)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!     find the wet bulb temperature for a given t and q
!     in a longitude height section
!     wet bulb temp is the temperature and spec humidity that is 
!     just saturated and has the same enthalpy
!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method: 
! a Newton method is used
! first guess uses an algorithm provided by John Petch from the UKMO
! we exclude points where the physical situation is unrealistic
! e.g. where the temperature is outside the range of validity for the
!      saturation vapor pressure, or where the water vapor pressure
!      exceeds the ambient pressure, or the saturation specific humidity is 
!      unrealistic
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
!
!     input arguments
!

  real(r8), intent(in) :: q        ! water vapor (kg/kg)
  real(r8), intent(in) :: t        ! temperature (K)
  real(r8), intent(in) :: p        ! pressure    (Pa)
  logical,  intent(in) :: use_ice  ! flag to include ice phase in calculations
!
! output arguments
!
  real(r8), intent(out) :: tsp      ! saturation temp (K)
  real(r8), intent(out) :: qsp      ! saturation mixing ratio (kg/kg)
  integer,  intent(out) :: status   ! flag representing state of output
                                    ! 0 => Successful convergence
                                    ! 1 => No calculation done: pressure or specific
                                    !      humidity not within usable range
                                    ! 2 => Run failed to converge
                                    ! 4 => Temperature fell below minimum
                                    ! 8 => Enthalpy not conserved
!
! local variables
!
  integer, parameter :: iter = 8    ! max number of times to iterate the calculation
  integer :: l                      ! iterator

  real(r8) es                   ! sat. vapor pressure
  real(r8) gam                  ! change in sat spec. hum. wrt temperature (times hltalt/cpair)
  real(r8) dgdt                 ! work variable
  real(r8) g                    ! work variable
  real(r8) hltalt               ! lat. heat. of vap.
  real(r8) qs                   ! spec. hum. of water vapor

! work variables
  real(r8) t1, q1, dt, dq
  real(r8) qvd
  real(r8) r1b, c1, c2
  real(r8), parameter :: dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
  real(r8), parameter :: dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
  real(r8) enin, enout

  ! Saturation specific humidity at this temperature
  if (use_ice) then
     call qsat(t, p, es, qs)
  else
     call qsat_water(t, p, es, qs)
  end if

  ! make sure a meaningful calculation is possible
  if (p <= 5._r8*es .or. qs <= 0._r8 .or. qs >= 0.5_r8 &
       .or. t < tmin .or. t > tmax) then
     status = 1
     ! Keep initial parameters when conditions aren't suitable
     tsp = t
     qsp = q
     enin = 1._r8
     enout = 1._r8

     return
  end if

  ! Prepare to iterate
  status = 2

  ! Get initial enthalpy
  if (use_ice) then
     call calc_hltalt(t,hltalt)
  else
     call no_ip_hltalt(t,hltalt)
  end if
  enin = tq_enthalpy(t, q, hltalt)

  ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
  c1 = hltalt*c3
  c2 = (t + 36._r8)**2
  r1b = c2/(c2 + c1*qs)
  qvd = r1b * (q - qs)
  tsp = t + ((hltalt/cpair)*qvd)

  ! Generate qsp, gam, and enout from tsp.
  if (use_ice) then
     call qsat(tsp, p, es, qsp, gam=gam, enthalpy=enout)
  else
     call qsat_water(tsp, p, es, qsp, gam=gam, enthalpy=enout)
  end if

  ! iterate on first guess
  do l = 1, iter

     g = enin - enout
     dgdt = -cpair * (1 + gam)

     ! New tsp
     t1 = tsp - g/dgdt
     dt = abs(t1 - tsp)/t1
     tsp = t1

     ! bail out if past end of temperature range
     if ( tsp < tmin ) then
        tsp = tmin
        ! Get latent heat and set qsp to a value
        ! that preserves enthalpy.
        if (use_ice) then
           call calc_hltalt(tsp,hltalt)
        else
           call no_ip_hltalt(tsp,hltalt)
        end if
        qsp = (enin - cpair*tsp)/hltalt
        enout = tq_enthalpy(tsp, qsp, hltalt)
        status = 4
        exit
     end if

     ! Re-generate qsp, gam, and enout from new tsp.
     if (use_ice) then
        call qsat(tsp, p, es, q1, gam=gam, enthalpy=enout)
     else
        call qsat_water(tsp, p, es, q1, gam=gam, enthalpy=enout)
     end if
     dq = abs(q1 - qsp)/max(q1,1.e-12_r8)
     qsp = q1

     ! if converged at this point, exclude it from more iterations
     if (dt < dttol .and. dq < dqtol) then
        status = 0
        exit
     endif
  end do

  ! Test for enthalpy conservation
  if (abs((enin-enout)/(enin+enout)) > 1.e-4_r8) status = 8

end subroutine findsp

! This portable module contains all CAM methods for estimating
! the saturation vapor pressure of water.
!
! wv_saturation provides CAM-specific interfaces and utilities
! based on these formulae.
!
! Typical usage of this module:
!
! Init:
! call wv_sat_methods_init(r8, tmelt, h2otrip, tboil, errstring)
!
! Get scheme index from a name string:
! scheme_idx = wv_sat_get_scheme_idx(scheme_name)
! if (.not. wv_sat_valid_idx(scheme_idx)) <throw some error>
!
! Get pressures:
! es = wv_sat_svp_water(t, scheme_idx)
! es = wv_sat_svp_ice(t, scheme_idx)
!
! Use ice/water transition range:
! es = wv_sat_svp_trice(t, ttrice, scheme_idx)
!
! Note that elemental functions cannot be pointed to, nor passed
! as arguments. If you need to do either, it is recommended to
! wrap the function so that it can be given an explicit (non-
! elemental) interface.





!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

! Get physical constants
!!$subroutine wv_sat_methods_init(kind, tmelt_in, h2otrip_in, tboil_in, &
!!$     ttrice_in, epsilo_in, errstring)
!!$  integer, intent(in) :: kind
!!$  real(r8), intent(in) :: tmelt_in
!!$  real(r8), intent(in) :: h2otrip_in
!!$  real(r8), intent(in) :: tboil_in
!!$  real(r8), intent(in) :: ttrice_in
!!$  real(r8), intent(in) :: epsilo_in
!!$  character(len=*), intent(out)  :: errstring
!!$
!!$  errstring = ' '
!!$
!!$  if (kind /= r8) then
!!$     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
!!$          kind,' was input kind but ',r8,' is internal kind.'
!!$     return
!!$  end if
!!$
!!$  if (ttrice_in < 0._r8) then
!!$     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
!!$          ttrice_in,' was input for ttrice, but negative range is invalid.'
!!$     return
!!$  end if
!!$
!!$  tmelt = tmelt_in
!!$  h2otrip = h2otrip_in
!!$  tboil = tboil_in
!!$  ttrice = ttrice_in
!!$  epsilo = epsilo_in
!!$
!!$  omeps = 1._r8 - epsilo
!!$
!!$end subroutine wv_sat_methods_init

! Look up index by name.
pure function wv_sat_get_scheme_idx(name) result(idx)
  character(len=*), intent(in) :: name
  integer :: idx
  
  select case (name)
  case("GoffGratch")
     idx = GoffGratch_idx
  case("MurphyKoop")
     idx = MurphyKoop_idx
  case("OldGoffGratch")
     idx = OldGoffGratch_idx
  case("Bolton")
     idx = Bolton_idx
  case default
     idx = Invalid_idx
  end select

end function wv_sat_get_scheme_idx

! Check validity of an index from the above routine.
pure function wv_sat_valid_idx(idx) result(status)
  integer, intent(in) :: idx
  logical :: status

  status = (idx /= Invalid_idx)

end function wv_sat_valid_idx

! Set default scheme (otherwise, Goff & Gratch is default)
! Returns a logical representing success (.true.) or
! failure (.false.).
function wv_sat_set_default(name) result(status)
  character(len=*), intent(in) :: name
  logical :: status

  ! Don't want to overwrite valid default with invalid,
  ! so assign to temporary and check it first.
  integer :: tmp_idx

  tmp_idx = wv_sat_get_scheme_idx(name)

  status = wv_sat_valid_idx(tmp_idx)

  if (status) default_idx = tmp_idx

end function wv_sat_set_default

! Reset default scheme to initial value.
! The same thing can be accomplished with wv_sat_set_default;
! the real reason to provide this routine is to reset the
! module for testing purposes.
subroutine wv_sat_reset_default()

  default_idx = initial_default_idx

end subroutine wv_sat_reset_default

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Get saturation specific humidity given pressure and SVP.
! Specific humidity is limited to range 0-1.
elemental function wv_sat_svp_to_qsat(es, p) result(qs)

  real(r8), intent(in) :: es  ! SVP
  real(r8), intent(in) :: p   ! Current pressure.
  real(r8) :: qs

  ! If pressure is less than SVP, set qs to maximum of 1.
  if ( (p - es) <= 0._r8 ) then
     qs = 1.0_r8
  else
     qs = epsilo*es / (p - omeps*es)
  end if

end function wv_sat_svp_to_qsat

! Get saturation "mass mixing ratio".
! It is almost always preferable to use saturation
! specific humidity rather than this mixing ratio,
! which blows up at low pressure.
elemental function wv_sat_svp_to_qmmr(es, p) result(qmmr)

  real(r8), intent(in) :: es  ! SVP
  real(r8), intent(in) :: p   ! Current pressure.
  real(r8) :: qmmr

  ! When this function is used in regions of very low
  ! pressure, we set it to "huge" rather than using
  ! the small denominator.
  if ( (p - es) < epsilon(1.0_r8)**2 ) then
     qmmr = huge(1.0_r8)
  else
     qmmr = epsilo*es / (p - es)
  end if

end function wv_sat_svp_to_qmmr

elemental subroutine wv_sat_qsat_water(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_water(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_water

elemental subroutine wv_sat_qsat_ice(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_ice(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_ice

elemental subroutine wv_sat_qsat_trans(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_trans(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_trans

!---------------------------------------------------------------------
! SVP INTERFACE FUNCTIONS
!---------------------------------------------------------------------

elemental function wv_sat_svp_water(t, idx) result(es)
  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(r8) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(GoffGratch_idx)
     es = GoffGratch_svp_water(t)
  case(MurphyKoop_idx)
     es = MurphyKoop_svp_water(t)
  case(OldGoffGratch_idx)
     es = OldGoffGratch_svp_water(t)
  case(Bolton_idx)
     es = Bolton_svp_water(t)
  end select

end function wv_sat_svp_water

elemental function wv_sat_svp_ice(t, idx) result(es)
  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(r8) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  select case (use_idx)
  case(GoffGratch_idx)
     es = GoffGratch_svp_ice(t)
  case(MurphyKoop_idx)
     es = MurphyKoop_svp_ice(t)
  case(OldGoffGratch_idx)
     es = OldGoffGratch_svp_ice(t)
  case(Bolton_idx)
     es = Bolton_svp_water(t)
  end select

end function wv_sat_svp_ice

elemental function wv_sat_svp_trans(t, idx) result (es)

  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx

  real(r8) :: es

  real(r8) :: esice      ! Saturation vapor pressure over ice
  real(r8) :: weight     ! Intermediate scratch variable for es transition

!
! Water
!
  if (t >= (tmelt - ttrice)) then
     es = wv_sat_svp_water(t,idx)
  else
     es = 0.0_r8
  end if

!
! Ice
!
  if (t < tmelt) then

     esice = wv_sat_svp_ice(t,idx)

     if ( (tmelt - t) > ttrice ) then
        weight = 1.0_r8
     else
        weight = (tmelt - t)/ttrice
     end if

     es = weight*esice + (1.0_r8 - weight)*es
  end if

end function wv_sat_svp_trans

!---------------------------------------------------------------------
! SVP METHODS
!---------------------------------------------------------------------

! Goff & Gratch (1946)

elemental function GoffGratch_svp_water(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! uncertain below -70 C
  es = 10._r8**(-7.90298_r8*(tboil/t-1._r8)+ &
       5.02808_r8*log10(tboil/t)- &
       1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/tboil))-1._r8)+ &
       8.1328e-3_r8*(10._r8**(-3.49149_r8*(tboil/t-1._r8))-1._r8)+ &
       log10(1013.246_r8))*100._r8

end function GoffGratch_svp_water

elemental function GoffGratch_svp_ice(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! good down to -100 C
  es = 10._r8**(-9.09718_r8*(h2otrip/t-1._r8)-3.56654_r8* &
       log10(h2otrip/t)+0.876793_r8*(1._r8-t/h2otrip)+ &
       log10(6.1071_r8))*100._r8

end function GoffGratch_svp_ice

! Murphy & Koop (2005)

elemental function MurphyKoop_svp_water(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! (good for 123 < T < 332 K)
  es = exp(54.842763_r8 - (6763.22_r8 / t) - (4.210_r8 * log(t)) + &
       (0.000367_r8 * t) + (tanh(0.0415_r8 * (t - 218.8_r8)) * &
       (53.878_r8 - (1331.22_r8 / t) - (9.44523_r8 * log(t)) + &
       0.014025_r8 * t)))

end function MurphyKoop_svp_water

elemental function MurphyKoop_svp_ice(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  ! (good down to 110 K)
  es = exp(9.550426_r8 - (5723.265_r8 / t) + (3.53068_r8 * log(t)) &
       - (0.00728332_r8 * t))

end function MurphyKoop_svp_ice

! Old CAM implementation, also labelled Goff & Gratch (1946)

! The water formula differs only due to compiler-dependent order of
! operations, so differences are roundoff level, usually 0.

! The ice formula gives fairly close answers to the current
! implementation, but has been rearranged, and uses the
! 1 atm melting point of water as the triple point.
! Differences are thus small but above roundoff.

! A curious fact: although using the melting point of water was
! probably a mistake, it mildly improves accuracy for ice svp,
! since it compensates for a systematic error in Goff & Gratch.

elemental function OldGoffGratch_svp_water(t) result(es)
  real(r8), intent(in) :: t
  real(r8) :: es
  real(r8) :: ps, e1, e2, f1, f2, f3, f4, f5, f

  ps = 1013.246_r8
  e1 = 11.344_r8*(1.0_r8 - t/tboil)
  e2 = -3.49149_r8*(tboil/t - 1.0_r8)
  f1 = -7.90298_r8*(tboil/t - 1.0_r8)
  f2 = 5.02808_r8*log10(tboil/t)
  f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
  f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
  f5 = log10(ps)
  f  = f1 + f2 + f3 + f4 + f5

  es = (10.0_r8**f)*100.0_r8
  
end function OldGoffGratch_svp_water

elemental function OldGoffGratch_svp_ice(t) result(es)
  real(r8), intent(in) :: t
  real(r8) :: es
  real(r8) :: term1, term2, term3

  term1 = 2.01889049_r8/(tmelt/t)
  term2 = 3.56654_r8*log(tmelt/t)
  term3 = 20.947031_r8*(tmelt/t)

  es = 575.185606e10_r8*exp(-(term1 + term2 + term3))
  
end function OldGoffGratch_svp_ice

! Bolton (1980)
! zm_conv deep convection scheme contained this SVP calculation.
! It appears to be from D. Bolton, 1980, Monthly Weather Review.
! Unlike the other schemes, no distinct ice formula is associated
! with it. (However, a Bolton ice formula exists in CLUBB.)

! The original formula used degrees C, but this function
! takes Kelvin and internally converts.

elemental function Bolton_svp_water(t) result(es)
  real(r8),parameter :: c1 = 611.2_r8
  real(r8),parameter :: c2 = 17.67_r8
  real(r8),parameter :: c3 = 243.5_r8

  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  es = c1*exp( (c2*(t - tmelt))/((t - tmelt)+c3) )

end function Bolton_svp_water

!-----------------------------------------------------------------------!
! Module to hold PBL-related subprograms that may be used with multiple !
! different vertical diffusion schemes.                                 !
!                                                                       !
! Public subroutines:                                                   !
!     
!     calc_obklen                                                       !
!                                                                       !
!------------------ History --------------------------------------------!
! Created: Apr. 2012, by S. Santos                                      !
!-----------------------------------------------------------------------!

elemental subroutine calc_ustar( t,    pmid, taux, tauy, &
                                 rrho, ustar)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate ustar and bottom level density (necessary for      !
  !  Obukhov length calculation).                                         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t         ! surface temperature
  real(r8), intent(in) :: pmid      ! midpoint pressure (bottom level)
  real(r8), intent(in) :: taux      ! surface u stress [N/m2]
  real(r8), intent(in) :: tauy      ! surface v stress [N/m2]

  real(r8), intent(out) :: rrho     ! 1./bottom level density
  real(r8), intent(out) :: ustar    ! surface friction velocity [m/s]

  rrho = rair * t / pmid
  ustar = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )
  
end subroutine calc_ustar

elemental subroutine calc_obklen( ths,  thvs, qflx, shflx, rrho, ustar, &
                                  khfs, kqfs, kbfs, obklen)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate Obukhov length and kinematic fluxes.               !
  !-----------------------------------------------------------------------!

  real(r8), intent(in)  :: ths           ! potential temperature at surface [K]
  real(r8), intent(in)  :: thvs          ! virtual potential temperature at surface
  real(r8), intent(in)  :: qflx          ! water vapor flux (kg/m2/s)
  real(r8), intent(in)  :: shflx         ! surface heat flux (W/m2)

  real(r8), intent(in)  :: rrho          ! 1./bottom level density [ m3/kg ]
  real(r8), intent(in)  :: ustar         ! Surface friction velocity [ m/s ]
  
  real(r8), intent(out) :: khfs          ! sfc kinematic heat flux [mK/s]
  real(r8), intent(out) :: kqfs          ! sfc kinematic water vapor flux [m/s]
  real(r8), intent(out) :: kbfs          ! sfc kinematic buoyancy flux [m^2/s^3]
  real(r8), intent(out) :: obklen        ! Obukhov length
  
  ! Need kinematic fluxes for Obukhov:
  khfs = shflx*rrho/cpair
  kqfs = qflx*rrho
  kbfs = khfs + zvir*ths*kqfs
  
  ! Compute Obukhov length:
  obklen = -thvs * ustar**3 / (g*vk*(kbfs + sign(1.e-10_r8,kbfs)))

end subroutine calc_obklen

elemental real(r8) function virtem(t,q)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate virtual temperature from temperature and specific  !
  !  humidity.                                                            !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t, q

  virtem = t * (1.0_r8 + zvir*q)

end function virtem

end module messy_vertdiff_tools
