! *************************************************************************
! *************************************************************************
! *************************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL SVOC
!
! Author:
! Mega Octaviani, MPIC, 03.2017 (last modified)
!
! Modifications according to MESSy standard (Patrick Joeckel, DLR, Mar 2019):
! - lfirst_day renamed to lstart
! - MESSYTENDENCY commented, because implementation is incomplete and
!   does not even compile
! *************************************************************************
! *************************************************************************
! *************************************************************************

MODULE messy_svoc

USE messy_main_constants_mem, ONLY: DP, pi
USE messy_main_blather,       ONLY: start_message, end_message

IMPLICIT NONE

INTRINSIC :: TRIM, MAX, MIN, SUM, LOG, LOG10, EXP

PRIVATE
SAVE

PUBLIC :: DP

CHARACTER(LEN = *), PARAMETER, PUBLIC :: modstr = 'svoc'
CHARACTER(LEN = *), PARAMETER, PUBLIC :: modver = '1.1'

LOGICAL, PUBLIC :: L_GP = .TRUE.  ! gridpoint calculation
LOGICAL, PUBLIC :: L_LG = .FALSE. ! lagrangian calculation
LOGICAL, PUBLIC :: l_svocpart = .TRUE. ! switch for gas-particle partitioning
LOGICAL, PUBLIC :: l_svocvola = .TRUE. ! switch for surface volatilization
LOGICAL, PUBLIC :: l_glacier = .FALSE. ! switch for glacier compartment
LOGICAL, PUBLIC :: l_landsnow = .FALSE. ! switch for snow compartment
LOGICAL, PUBLIC :: l_pahderiv = .FALSE. ! switch for nitro-pah derivative

PUBLIC :: svoc_read_nml_ctrl
PUBLIC :: svoc_volatilizations
PUBLIC :: svoc_partition
PUBLIC :: svoc_degradation
PUBLIC :: svoc_drydep
PUBLIC :: svoc_wetdep
PUBLIC :: svoc_hetchem_phe
PUBLIC :: svoc_hetchem_bap
!! mz_jw_20160918+
PUBLIC :: svoc_hetchem_npah
!! mz_jw_20160918-

CONTAINS

! *************************************************************************

!! mz_jw_20160918+
! *************************************************************************

SUBROUTINE svoc_hetchem_npah(kproma, &
                             jno2, &
                             pxtp1p, &
                             zdcdt)

!!! photocatalysed loss of particulate bound npah

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN) :: kproma
REAL(DP), DIMENSION(kproma), INTENT(IN) :: jno2 ! j_no2
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pxtp1p ! mixing ratio of particle
REAL(DP), DIMENSION(kproma), INTENT(OUT) :: zdcdt ! decay rate

! local
INTEGER :: jl
REAL(DP), PARAMETER :: scal = 0.05

zdcdt(:) = 0._dp

DO jl = 1, kproma
   zdcdt(jl) = -scal * jno2(jl) * pxtp1p(jl)
END DO

END SUBROUTINE svoc_hetchem_npah

! *************************************************************************
!! mz_jw_20160918-

! *************************************************************************

SUBROUTINE svoc_hetchem_bap(kproma, &
                            param_bapo3, &
                            press, &
                            relhum, &
                            temp, &
                            pco3, &
                            pxtp1p, &
                            zdcdt)

!!!

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN) :: kproma
INTEGER, INTENT(IN) :: param_bapo3 ! 1: poeschl, 2: kwamena, 3: ROI(T), 4: kahan
REAL(DP), DIMENSION(kproma), INTENT(IN) :: press ! air pressure
REAL(DP), DIMENSION(kproma), INTENT(IN) :: relhum ! relative humidity
REAL(DP), DIMENSION(kproma), INTENT(IN) :: temp ! air temperature
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pco3 ! o3 concentration
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pxtp1p ! mixing ratio of particle
REAL(DP), DIMENSION(kproma), INTENT(OUT) :: zdcdt ! decay rate

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_hetchem_bap'
 ! rkmax: maximum pseudo-first-order BaP decay rate coefficient [s-1]
REAL(DP), PARAMETER :: rkmax_1 = 0.015_dp ! poeschl scheme
REAL(DP), PARAMETER :: rkmax_2 = 0.060_dp ! kwamena scheme
 ! bko3: langmuir adsoprtion equilibrium constant for o3
REAL(DP), PARAMETER :: bko3_1 = 2.8e-13_dp ! poeschl scheme
REAL(DP), PARAMETER :: bko3_2 = 2.8e-15_dp ! kwamena scheme
 ! fitting parameters for kahan scheme
REAL(DP), PARAMETER :: ra = 5.5e-3_dp  ! s-1
REAL(DP), PARAMETER :: rb = 2.8e+15_dp  ! molec cm-3
REAL(DP) :: rkact ! rate constant [s-1]
INTEGER :: jl

!CALL start_message(modstr, 'BaP HETEROGENEOUS OXIDATION BY O3', substr)

zdcdt(:) = 0._dp   ! decay [mol(svoc) mol(air)-1 s-1]

DO jl = 1, kproma
   rkact = 0._dp

   IF (param_bapo3 == 1) THEN !!! poeschl scheme
      rkact = rkmax_1 * bko3_1 * pco3(jl) / (1 + bko3_1 * pco3(jl))

   ELSE IF (param_bapo3 == 2) THEN !!! kwamena scheme
      rkact = rkmax_2 * bko3_2 * pco3(jl) / (1 + bko3_2 * pco3(jl))

   ELSE IF (param_bapo3 == 3) THEN !!! ROI(T) scheme
      CALL svoc_bapo3_roi(press(jl), &
                          relhum(jl), &
                          temp(jl), &
                          pco3(jl), &
                          rkact)

   ELSE IF (param_bapo3 == 4) THEN !!! kahan scheme
      rkact = ra * pco3(jl) / (rb + pco3(jl))
   
   END IF

   zdcdt(jl) = -rkact * pxtp1p(jl)
END DO

!CALL end_message(modstr, 'BaP HETEROGENEOUS OXIDATION BY O3', substr)

END SUBROUTINE svoc_hetchem_bap
 
! *************************************************************************

SUBROUTINE svoc_hetchem_phe(kproma, &
                            rtotarea, &
                            pco3, &
                            pxtp1p, &
                            zdcdt)
                       
!!! heterogeneous oxidation of phenanthrene with ozone

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN) :: kproma
REAL(DP), DIMENSION(kproma), INTENT(IN) :: rtotarea ! [cm2 cm-3]
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pco3 ! o3 concentration
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pxtp1p ! mixing ratio of particle
REAL(DP), DIMENSION(kproma), INTENT(OUT) :: zdcdt ! decay rate

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_hetchem_phe'
REAL(DP), PARAMETER :: rko3p = 6.2e-17_dp ! [cm4 molec-1 s-1]
INTEGER :: jl

!CALL start_message(modstr, 'PHE HETEROGENEOUS OXIDATION BY O3', substr)

zdcdt(:) = 0._dp

DO jl = 1, kproma
   ! decay = [cm4 molec-1 s-1] * [cm2 cm-3] * [molec cm-3] * [mol mol-1]
   ! decay = [mol mol-1 s-1]
   IF (rtotarea(jl) .GE. 1e-20_dp) &
      zdcdt(jl) = -rko3p * rtotarea(jl) * pco3(jl) * pxtp1p(jl)
END DO

!CALL end_message(modstr, 'PHE HETEROGENEOOUS OXIDATION BY O3', substr)

END SUBROUTINE svoc_hetchem_phe

! *************************************************************************

SUBROUTINE svoc_wetdep(kproma, klev, dtime, &
                       molmass, &
                       loland, &
                       loglac, &
                       pcvs, &
                       pvgrat, &
                       pwetflxg_ls, &
                       pwetflxg_cv, &
                       pwetflxp, &
                       zburd, &
                       zdepo)

!!!

USE messy_main_constants_mem, ONLY: N_A, M_air

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN) :: kproma, klev
REAL(DP), INTENT(IN) :: dtime
REAL(DP), INTENT(IN) :: molmass ! [g mol-1]
LOGICAL, DIMENSION(kproma), INTENT(IN) :: loland ! logical mask land
LOGICAL, DIMENSION(kproma), INTENT(IN) :: loglac ! logical mask glacier
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pcvs ! snow fraction
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pvgrat ! vegetation fraction
 ! pwetflxg_ls: wet depo flux of gas from stratiform rain [molec m-2 s-1]
 ! pwetflxg_cv: wet depo flux of gas from convective rain [molec m-2 s-1]
 ! pwetflxp: total wet flux of particles [mol(trac) mol(air)-1 kg(air) m-2 s-1]
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pwetflxg_ls
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pwetflxg_cv
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pwetflxp
REAL(DP), DIMENSION(kproma, klev), INTENT(INOUT) :: zburd ! burden
REAL(DP), DIMENSION(kproma), INTENT(OUT) :: zdepo ! deposition

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_wetdep'
REAL(DP), PARAMETER :: g2kg = 1.e-3_dp
REAL(DP), DIMENSION(kproma) :: wetflux_g ! wetdep flux for gas
REAL(DP), DIMENSION(kproma) :: wetflux_p ! wetdep flux for particle
REAL(DP), DIMENSION(kproma) :: wetflux_tot ! total wetdep flux
INTEGER :: jl, jk

!CALL start_message(modstr, 'PAH WET DEPOSITION', substr)

! convert wetflux to [kg(trac) m-2 s-1]
wetflux_g = (pwetflxg_ls(1: kproma) + pwetflxg_cv(1: kproma)) * &
            molmass * g2kg / N_A

wetflux_p = pwetflxp(1: kproma) * molmass / M_air

wetflux_tot = wetflux_g + wetflux_p

DO jl = 1, kproma
   IF (loland(jl)) THEN ! land
      IF (loglac(jl)) THEN
         IF (l_glacier) THEN
            jk = klev - 4
            zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * dtime

         ELSE
            jk = klev
            zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * &
                           (1._dp - pvgrat(jl)) * dtime

            jk = klev - 1
            zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * &
                           pvgrat(jl) * dtime
         END IF

      ELSE ! not loglac
         IF (l_landsnow) THEN
            jk = klev
            zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * &
                           (1._dp - pcvs(jl)) * (1._dp - pvgrat(jl)) * dtime

            jk = klev - 1
            zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * &
                           (1._dp - pcvs(jl)) * pvgrat(jl) * dtime

            jk = klev - 3
            zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * &
                            pcvs(jl) * dtime

         ELSE
            jk = klev
            zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * &
                           (1._dp - pvgrat(jl)) * dtime
            
            jk = klev - 1
            zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * &
                           pvgrat(jl) * dtime
         END IF
      END IF

   ELSE ! ocean
      jk = klev - 2
      zburd(jl, jk) = zburd(jl, jk) + wetflux_tot(jl) * dtime
   END IF

   zdepo(jl) = wetflux_tot(jl)
END DO

!CALL end_message(modstr, 'PAH WET DEPOSITION', substr)

END SUBROUTINE svoc_wetdep

! *************************************************************************

SUBROUTINE svoc_drydep(kproma, klev, nmod, &
                       dtime, tstep, &
                       param_soilv, &
                       molmass, &
                       ha, hb, &
                       loland, &
                       loglac, &
                       pcvs, &
                       pvgrat, &
                       ptslm1, &
                       vd_g, &
                       vd_p, &
                       zdz, &
                       zdp, &
                       pxtp1g, &
                       pxtp1p, &
                       zxtems_g, &
                       zxtems_p, &
                       zburd, &
                       zdepo)

!!! 

USE messy_main_constants_mem, ONLY: R_gas, g, M_air

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN) :: kproma, klev, nmod
REAL(DP), INTENT(IN) :: dtime, tstep
INTEGER, INTENT(IN) :: param_soilv ! 1: smit, 2: jury
REAL(DP), INTENT(IN) :: molmass ! [g mol-1]
REAL(DP), INTENT(IN) :: ha  ! henry parameter A [M atm-1]
REAL(DP), INTENT(IN) :: hb  ! henry parameter B [K]
LOGICAL, DIMENSION(kproma), INTENT(IN) :: loland ! logical mask land
LOGICAL, DIMENSION(kproma), INTENT(IN) :: loglac ! logical mask glacier
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pcvs ! snow fraction
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pvgrat ! vegetation fraction
REAL(DP), DIMENSION(kproma), INTENT(IN) :: ptslm1 ! sfc temp at t-dt
REAL(DP), DIMENSION(kproma), INTENT(IN) :: vd_g ! ddep vel of gas [cm s-1]
REAL(DP), DIMENSION(kproma, nmod), INTENT(IN) :: vd_p ! ddep vel of particle
REAL(DP), DIMENSION(kproma), INTENT(IN) :: zdz ! layer thickness [m]
REAL(DP), DIMENSION(kproma), INTENT(IN) :: zdp ! delta pressure
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pxtp1g ! vmr of gas
REAL(DP), DIMENSION(kproma, nmod), INTENT(IN) :: pxtp1p ! vmr of particle
REAL(DP), DIMENSION(kproma), INTENT(INOUT) :: zxtems_g ! gas emission
REAL(DP), DIMENSION(kproma, nmod), INTENT(INOUT) :: zxtems_p ! particle emission
REAL(DP), DIMENSION(kproma, klev), INTENT(INOUT) :: zburd ! burden
REAL(DP), DIMENSION(kproma), INTENT(OUT) :: zdepo ! deposition flux

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_drydep'
REAL(DP), PARAMETER :: tstd = 298.15_dp  ! standard temperature [K]
REAL(DP), PARAMETER :: cnv = 101.325   ! convert M atm-1 to mol m-3 Pa-1
REAL(DP), PARAMETER :: cm2m = 1.e-2_dp ! convert meter to cm
REAL(DP), DIMENSION(kproma) :: ddepflux_g ! drydep flux for gas
REAL(DP), DIMENSION(kproma, nmod) :: ddepflux_p ! drydep flux for particle
REAL(DP), DIMENSION(kproma) :: ddepflux_tot ! total drydep flux
REAL(DP) :: zqvd, vd_int, vde_g, vde_p
REAL(DP) :: tmp1, tmp2
REAL(DP) :: henry, kaw, zlfrac
INTEGER :: jl, jk, jm

!CALL start_message(modstr, 'PAH DRY DEPOSITION', substr)

ddepflux_g(:) = 0._dp
ddepflux_p(:, :) = 0._dp

DO jl = 1, kproma
   ! calculate effective dry deposition velocity for gas
   ! following eq 28 in kerkweg et al. (2006) 
   zqvd = 1._dp / MAX(1.e-20_dp, vd_g(jl) * cm2m)
   vd_int = -1._dp / (zqvd * zdz(jl) / tstep)
   vde_g = zdz(jl) / tstep * (1._dp - EXP(vd_int))

   ! drydep flux in [mol(trac) mol(air)-1 kg(air) m-2 s-1]
   ddepflux_g(jl) = pxtp1g(jl) * zdp(jl) * vde_g / (g * zdz(jl))

   ! follow function drydep_posfinit
   ! ensure positive tracer concentrations by adjusting ddepflux
   ! -> gas
   tmp1 = pxtp1g(jl) - ddepflux_g(jl) * tstep / (zdp(jl) / g)

   IF (tmp1 < 0._dp) THEN
      tmp2 = pxtp1g(jl) * (zdp(jl) / g) / tstep
   ELSE
      tmp2 = ddepflux_g(jl)
   END IF

   ! -> ddep updates emission with ddepflux, the following step
   ! -> is to correct the emission value when ddepflux is negative
   zxtems_g(jl) = zxtems_g(jl) + ddepflux_g(jl) - tmp2

   ! -> correct dry depo flux of gas
   ddepflux_g(jl) = tmp2

   DO jm = 1, nmod
      ! calculate effective dry deposition velocity for particle
      zqvd = 1._dp / MAX(1.e-20_dp, vd_p(jl, jm) * cm2m)
      vd_int = -1._dp / (zqvd * zdz(jl) / tstep)
      vde_p = zdz(jl) / tstep * (1._dp - EXP(vd_int))

      ! drydep flux in [mol(trac) mol(air)-1 kg(air) m-2 s-1]
      ddepflux_p(jl, jm) = pxtp1p(jl, jm) * zdp(jl) * vde_p / (g * zdz(jl))

      ! follow function drydep_posfinit
      ! ensure positive tracer concentrations by adjusting ddepflux
      ! -> particle
      tmp1 = pxtp1p(jl, jm) - ddepflux_p(jl, jm) * tstep / (zdp(jl) / g)

      IF (tmp1 < 0._dp) THEN
         tmp2 = pxtp1p(jl, jm) * (zdp(jl) / g) / tstep
      ELSE
         tmp2 = ddepflux_p(jl, jm)
      END IF

      ! -> correct drydep flux of particles
      ddepflux_p(jl, jm) = tmp2

      ! -> update emision of particles
      zxtems_p(jl, jm) = zxtems_p(jl, jm) - ddepflux_p(jl, jm)
   END DO

   ! total drydep flux in kg(trac) m-2 s-1
   IF (param_soilv == 2 .AND. loland(jl)) THEN  ! param_soilv = Jury
      ! add liquid phase dry deposition over land to compensate
      ! liquid phase diffusion in the jury soil model
      henry = ha * EXP(hb * &
              ((1._dp / ptslm1(jl)) - (1._dp / tstd))) / cnv

      ! air-liquid partitioning coefficient [-]
      kaw = 1._dp / (R_gas * ptslm1(jl) * henry)

      zlfrac = 1._dp / kaw
      IF (zlfrac > 1._dp) zlfrac = 1._dp
      IF (zlfrac < 0._dp) zlfrac = 0._dp

      ddepflux_tot(jl) = (ddepflux_g(jl) + &
                        SUM(ddepflux_p(jl, 1: nmod)) * &
                        (1._dp + zlfrac)) * molmass / M_air

   ELSE
      ddepflux_tot(jl) = (ddepflux_g(jl) + SUM(ddepflux_p(jl, 1: nmod))) * &
                         molmass / M_air
   END IF

   IF (loland(jl)) THEN
      IF (loglac(jl)) THEN
         IF (l_glacier) THEN
            jk = klev - 4
            zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * dtime

         ELSE
            jk = klev
            zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * &
                           (1._dp - pvgrat(jl)) * dtime

            jk = klev - 1
            zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * &
                           pvgrat(jl) * dtime
         END IF
      ELSE ! not loglac
         IF (l_landsnow) THEN
            jk = klev
            zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * &
                           (1._dp - pcvs(jl)) * (1._dp - pvgrat(jl)) * dtime

            jk = klev - 1
            zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * &
                           (1._dp - pcvs(jl)) * pvgrat(jl) * dtime

            jk = klev - 3
            zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * &
                            pcvs(jl) * dtime

         ELSE
            jk = klev
            zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * &
                           (1._dp - pvgrat(jl)) * dtime

            jk = klev - 1
            zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * &
                           pvgrat(jl) * dtime
         END IF
      END IF

   ELSE ! ocean
      ! if l_svocvola = true and airsea is switched on:
      ! airsea computes gaseous  net fluxes (vola - depo), so in fact
      ! only ddepflux_p is used to update ocean burden below
      ! (ddepflux_tot = ddepflux_p because vd_g = 0 => ddepflux_g = 0)
      jk = klev - 2
      zburd(jl, jk) = zburd(jl, jk) + ddepflux_tot(jl) * dtime
   END IF

   zdepo(jl) = ddepflux_tot(jl)
END DO

!CALL end_message(modstr, 'PAH DRY DEPOSITION', substr)

END SUBROUTINE svoc_drydep

! *************************************************************************

SUBROUTINE svoc_degradation(kproma, klev, dtime, &
                            lstart, &
                            rksoil, rkocean, &
                            loland, &
                            ptslm1, &
                            pmld, &
                            zburd, &
                            zdegr)

!!! adjust svoc reservoirs in soil, vegetation, and ocean surfaces
!!! and compute svoc loss due to decay at the surfaces
!!!
!!! follows first-order decay rate with the rate coefficient
!!! is assumed to double per 10K temperature increase

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN) :: kproma, klev
REAL(DP), INTENT(IN) :: dtime
LOGICAL, INTENT(IN) :: lstart
REAL(DP), INTENT(IN) :: rksoil, rkocean ! degr rate coeff
LOGICAL, DIMENSION(kproma), INTENT(IN) :: loland ! logical mask land
REAL(DP), DIMENSION(kproma), INTENT(IN) :: ptslm1 ! sfc temp at t-dt
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pmld ! oce mixed layer depth [m]
REAL(DP), DIMENSION(kproma, klev), INTENT(INOUT) :: zburd ! burden
REAL(DP), DIMENSION(kproma, klev), INTENT(INOUT) :: zdegr ! degradation

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_degradation'
REAL(DP), PARAMETER :: tstds = 298.15_dp ! standard temp of soil surf
REAL(DP), PARAMETER :: tstdo = 273.15_dp ! standard temp of oce surf
REAL(DP) :: zrate ! temp-dependent degr rate coeff [s-1]
INTEGER :: jl, jk

!CALL start_message(modstr, 'PAH BIOTIC DEGRADATION', substr)

DO jl = 1, kproma
   IF (loland(jl)) THEN
      zrate = rksoil * 2._dp ** ((ptslm1(jl) - tstds) / 10._dp)

      ! degradation and burden in soil
      jk = klev
      zdegr(jl, jk) = zburd(jl, jk) * zrate
      zburd(jl, jk) = zburd(jl, jk) * EXP(-zrate * dtime)

      ! in vegetation
      jk = klev - 1
      zdegr(jl, jk) = zburd(jl, jk) * zrate
      zburd(jl, jk) = zburd(jl, jk) * EXP(-zrate * dtime)

   ELSE
      ! in ocean
      jk = klev - 2

      IF (.NOT. lstart) THEN
         IF (pmld(jl) > 0._dp) THEN
            zrate = rkocean * 2._dp ** ((ptslm1(jl) - tstdo) / 10._dp)

            zdegr(jl, jk) = zburd(jl, jk) * zrate
            zburd(jl, jk) = zburd(jl, jk) * EXP(-zrate * dtime)
         END IF
      END IF ! lstart

   END IF ! loland
END DO
         
!CALL end_message(modstr, 'PAH BIOTIC DEGRADATION', substr)

END SUBROUTINE svoc_degradation

! *************************************************************************

SUBROUTINE svoc_partition(kproma, klev, &
                          param_part, &
                          ha, hb, &
                          rvapp, rhvap, &
                          rkoab, rkoam, &
                          rhabc, & !! mz_jw_20170220
                          rpplfer, &
                          ptair, &
                          prhum, &
                          pctsp, &
                          pcwsom, &
                          pcwiom, &
                          pcbc, &
                          pcso4, &
                          pcss, &
                          rarea, &
                          zkp, &
                          ztheta)

!!! calculate fraction of particle svocs using different parameterization

USE messy_main_constants_mem, ONLY: R_gas

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN) :: kproma, klev ! grid dimension
INTEGER, INTENT(IN) :: param_part ! 1: L+L, 2: Fin, 3: J+P
REAL(DP), INTENT(IN) :: ha, hb ! henry parameters
REAL(DP), INTENT(IN) :: rvapp ! saturated vapor pressure [Pa] at 298 K
REAL(DP), INTENT(IN) :: rhvap ! heat of vaporization [J mol-1] at 298 K
REAL(DP), INTENT(IN) :: rkoab  ! intercept of temp-dependent rlogkoa
REAL(DP), INTENT(IN) :: rkoam  ! slope of temp-dependent rlogkoa [K]
!! mz_jw_20170220+
REAL(DP), INTENT(IN) :: rhabc !! enthalpy of adsorption on bc [J mol-1]
!! mz_jw_20170220-
REAL(DP), DIMENSION(6), INTENT(IN) :: rpplfer ! pp-lfers constants
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: ptair ! air temp.
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: prhum ! rel. humidity
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: pctsp ! tsp conc [ug m-3]
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: pcwsom ! ws om conc [ug m-3]
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: pcwiom ! wi om conc [ug m-3]
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: pcbc ! bc conc [ug m-3]
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: pcso4 ! sulphate conc [ug m-3]
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: pcss ! sea-salt conc [ug m-3]
REAL(DP), DIMENSION(kproma, klev), INTENT(IN) :: rarea ! [cm2 cm-3]
REAL(DP), DIMENSION(kproma, klev), INTENT(OUT) :: zkp ! partition coeff. [m3 ug-1]
REAL(DP), DIMENSION(kproma, klev), INTENT(OUT) :: ztheta ! particle fract.

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_partition'
REAL(DP) :: zfinizio, zjunge, zlohlam, zpplfers
REAL(DP) :: totocconc
INTEGER :: jl, jk

!CALL start_message(modstr, 'GAS-PARTICLE PARTITIONING', substr)

zkp = 0._dp 
ztheta = 0._dp  ! without partitioning to aerosols

IF (l_svocpart) THEN
   DO jk = 1, klev
      DO jl = 1, kproma
         totocconc = pcwsom(jl, jk) + pcwiom(jl, jk)

         IF (param_part == 1) THEN ! lohmann-lammel
            CALL svoc_part_lohlam(rkoab, rkoam, &
                          !! if ksa is from ksw and kaw (obsolete)
                                 !ha, hb, &
                                 !rlogksw, &
                          !! if ksa is from vap press
                                  rvapp, rhvap, &
                                  ptair(jl, jk), &
                                  totocconc, &
                                  pcbc(jl, jk), &
                                  zkp(jl, jk), &
                                  zlohlam)

            IF (pctsp(jl, jk) > 0._dp) THEN
               zkp(jl, jk) = zkp(jl, jk) / pctsp(jl, jk)
            ELSE
               zkp(jl, jk) = 0._dp
            END IF

            ztheta(jl, jk) = zlohlam

         ELSE IF (param_part == 2) THEN ! finizio
            CALL svoc_part_finizio(rkoab, rkoam, &
                                   ptair(jl, jk), &
                                   pctsp(jl, jk), &
                                   zkp(jl, jk), &
                                   zfinizio)

            ztheta(jl, jk) = zfinizio

         ELSE IF (param_part == 3) THEN ! junge-pankow
            CALL svoc_part_junge(rvapp, rhvap, &
                                 ptair(jl, jk), &
                                 rarea(jl, jk), &
                                 zjunge)

            ztheta(jl, jk) = zjunge

         ELSE IF (param_part == 4) THEN ! pp-lfers
            CALL svoc_part_pplfers(rhvap, &
                                   rhabc, & !! mz_jw_20170220
                                   rpplfer, &
                                   ptair(jl, jk), &
                                   prhum(jl, jk), &
                                   pcbc(jl, jk), &
                                   pcwsom(jl, jk), &
                                   pcwiom(jl, jk), &
                                   pcso4(jl, jk), &
                                   pcss(jl, jk), &
                                   zkp(jl, jk), &
                                   zpplfers)

            IF (pctsp(jl, jk) > 0._dp) THEN
               zkp(jl, jk) = zkp(jl, jk) / pctsp(jl, jk)
            ELSE
               zkp(jl, jk) = 0._dp
            END IF

            ztheta(jl, jk) = zpplfers
         END IF

      END DO ! kproma
   END DO ! klev
END IF  ! l_svocpart

!CALL end_message(modstr, 'GAS-PARTICLE PARTITIONING', substr)

END SUBROUTINE svoc_partition

! *************************************************************************

SUBROUTINE svoc_volatilizations(kproma, klev, dtime, &
                                param_soilv, &
                                molmass, &
                                rwsol, rhsol, &
                                rvapp, rhvap, &
                                rlogkow, rhsub, &
                                ha, hb, &
                                loland, &
                                loglac, &
                                pforest, &
                                pvgrat, &
                                psn, &
                                pgld, &
                                pcvs, &
                                pws, &
                                pwsmx, &
                                ptslm1, &
                                pfoms, &
                                prhos, &
                                zxtems, &
                                zcvs_old, &
                                zburd, &
                                zvola)

!!! calculate volatilization from soil, vegetation, glacier and snow

USE messy_main_constants_mem, ONLY: R_gas, M_air

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN) :: kproma, klev
REAL(DP), INTENT(IN) :: dtime
INTEGER, INTENT(IN) :: param_soilv ! 1: smit, 2: jury
REAL(DP), INTENT(IN) :: molmass ! [g mol-1]
REAL(DP), INTENT(IN) :: rwsol ! water solubility [mg L-1] at 298 K
REAL(DP), INTENT(IN) :: rhsol ! heat of solution [J mol-1]
REAL(DP), INTENT(IN) :: rvapp ! saturated vapor pressure [Pa] at 298 K
REAL(DP), INTENT(IN) :: rhvap ! heat of vaporization [J mol-1]
REAL(DP), INTENT(IN) :: rlogkow ! octanol-water partition coefficient
REAL(DP), INTENT(IN) :: rhsub ! heat of sublimation [J mol-1]
REAL(DP), INTENT(IN) :: ha, hb ! henry parameters
LOGICAL, DIMENSION(kproma), INTENT(IN) :: loland ! logical mask land
LOGICAL, DIMENSION(kproma), INTENT(IN) :: loglac ! logical mask glacier
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pforest ! forest fraction
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pvgrat ! vegetation fraction
REAL(DP), DIMENSION(kproma), INTENT(IN) :: psn ! snow depth [m water eq]
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pgld ! glacier depth (incl. snow)
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pcvs ! snow fraction
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pws ! surface soil wetness
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pwsmx ! soil fld capacity
REAL(DP), DIMENSION(kproma), INTENT(IN) :: ptslm1 ! sfc temp at t-dt
REAL(DP), DIMENSION(kproma), INTENT(IN) :: pfoms ! fraction of om in soil
REAL(DP), DIMENSION(kproma), INTENT(IN) :: prhos ! dry bulk density of soil
REAL(DP), DIMENSION(kproma), INTENT(INOUT) :: zxtems ! tracer emission
REAL(DP), DIMENSION(kproma), INTENT(INOUT) :: zcvs_old ! snow frac at t-dt
REAL(DP), DIMENSION(kproma, klev), INTENT(INOUT) :: zburd ! burden [kg(trac) m-2]
REAL(DP), DIMENSION(kproma, klev), INTENT(OUT) :: zvola ! vola [kg(trac) m-2 s-1]

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_volatilizations'
REAL(DP), PARAMETER :: tstd = 298.15_dp ! standard temperature [K]
REAL(DP), PARAMETER :: cnv = 101.325_dp ! convert M atm-1 to mol m-3 Pa-1
REAL(DP), PARAMETER :: zrho_sn = 330._dp ! density of snow in [kg m-3]
REAL(DP), PARAMETER :: porefr = 5.e-1_dp ! volume fraction of air in soil
REAL(DP) :: henry, kaw, pksl, plkoc, volg, volw
REAL(DP) :: zvolglac, zvolsnow, zvolveg, zvolsoil
REAL(DP) :: enclf
INTEGER :: jl, jk

!CALL start_message(modstr, 'LAND SURFACE VOLATILIZATIONS', substr)

DO jl = 1, kproma
   IF (loland(jl)) THEN
      zvolsnow = 0._dp ! always 0 except for snow calculation
      zvolglac = 0._dp ! always 0 except for glacier calculation

      !!! volatilization from vegetation (smit 1997,1998)
      CALL svoc_vola_veg(rhvap, rvapp, &
                         ptslm1(jl), &
                         zvolveg)

      ! henry coefficient [mol m-3 Pa-1]
      henry = ha * EXP(hb * &
              (1._dp / ptslm1(jl) - 1._dp / tstd)) / cnv

      ! air-liquid partitioning coefficient [-]
      kaw = 1._dp / (R_gas * ptslm1(jl) * henry)

      ! solid organic carbon-water partitioning coefficient [-]
      ! method 1: if kom [m3/kg] is available
      !  pksl = rkom * pfoms(jl)
      ! method 2: if kom is not available, first derive koc [mL/g]
      !  from kow based on their regression relation, then use 
      !  kom = 0.56 * koc * 1.e-3 [m3/kg] following mackay et al. (2006)
      ! a. mackay and boethling (2010) => for pahs, as used here
      !    koc = 10 ** (0.823 * rlogkow + 0.727) [mL/g]
      ! b. rao and davidson (1980) => for pesticides, as used in
      !    jury et al. (1983)
      !    koc = 10 ** (1.029 * rlogkow - 0.18) [mL/g]
      plkoc = 0.823_dp * rlogkow + 0.727_dp !!! log koc
      pksl = 0.56 * 10 ** (plkoc) * 1.e-3_dp * pfoms(jl)

      ! volume fraction of gaseous organics [-]
      volg = (1._dp - pws(jl) / pwsmx(jl)) * porefr

      ! volume fraction of water 
      volw = porefr - volg

      !!! volatilization from soil
      IF (param_soilv == 1) THEN
         !!! follows smit (1997, 1998)
         CALL svoc_vola_soil_smit(dtime, &
                                  kaw, pksl, &
                                  volg, volw, &
                                  pws(jl), pwsmx(jl), &
                                  prhos(jl), &
                                  zvolsoil)

      ELSE IF (param_soilv == 2) THEN
         !!! follows jury (1983)
         CALL svoc_vola_soil_jury(dtime, &
                                  pforest(jl), prhos(jl), &
                                  porefr, kaw, pksl, &
                                  volg, volw, &
                                  zvolsoil)
      END IF ! param_soilv

      IF (loglac(jl)) THEN
         IF (l_glacier) THEN
            jk = klev - 4

            !!! volatilization from glacier (wania, 1997)
            CALL svoc_vola_glac_snow(loglac(jl), dtime, &
                                     molmass, &
                                     zrho_sn, &
                                     ha, hb, &
                                     rhsub, rhsol, &
                                     rwsol, &
                                     pgld(jl), &
                                     zburd(jl, jk), &
                                     ptslm1(jl), &
                                     zvolglac)
         END IF
 
      ELSE ! not loglac
         !!! volatilization from snow over land (wania, 1997)
         !!! requirements: no glacier, snow cover frac > 0, and 
         !!! landsnow calc is desired
         IF (pcvs(jl) .GT. 0._dp .AND. l_landsnow) THEN
            jk = klev - 3
            CALL svoc_vola_glac_snow(loglac(jl), dtime, &
                                     molmass, &
                                     zrho_sn, &
                                     ha, hb, &
                                     rhsub, rhsol, &
                                     rwsol, &
                                     psn(jl), &
                                     zburd(jl, jk), &
                                     ptslm1(jl), &
                                     zvolsnow)
         END IF
      END IF

      !!! update emission, compartment burden and volatilization
      ! convert flux unit to emission unit by multiply with mw_air/mw_svoc
      ! flux -> kg(trac) m-2 s-1;
      ! zxtems -> mol(trac) mol(air)-1 kg(air) m-2 s-1
      IF (loglac(jl)) THEN
         IF (l_glacier) THEN
            jk = klev - 4
            zxtems(jl) = zxtems(jl) + &
                         zvolglac * zburd(jl, jk) * M_air / molmass

            zvola(jl, jk) = zburd(jl, jk) * zvolglac
            zburd(jl, jk) = zburd(jl, jk) - zburd(jl, jk) * zvolglac * dtime

         ELSE ! not l_glacier
            jk = klev - 1
            zxtems(jl) = zxtems(jl) + (zvolveg * zburd(jl, jk) + &
                         zvolsoil * zburd(jl, klev)) * M_air / molmass

            zvola(jl, jk) = zburd(jl, jk) * zvolveg
            zburd(jl, jk) = zburd(jl, jk) - zburd(jl, jk) * zvolveg * dtime

            zvola(jl, klev) = zburd(jl, klev) * zvolsoil
            zburd(jl, klev) = zburd(jl, klev) - &
                              zburd(jl, klev) * zvolsoil * dtime
         END IF

      ELSE ! not loglac
         IF (l_landsnow) THEN
            jk = klev - 3
            zxtems(jl) = zxtems(jl) + ((zvolveg * zburd(jl, klev - 1) + &
                         zvolsoil * zburd(jl, klev)) * (1._dp - pcvs(jl)) + &
                         zvolsnow * zburd(jl, jk)) * M_air / molmass

            zvola(jl, jk) = zburd(jl, jk) * zvolsnow
            zburd(jl, jk) = zburd(jl, jk) - zburd(jl, jk) * zvolsnow * dtime

            jk = klev - 1
            zvola(jl, jk) = zburd(jl, jk) * zvolveg * (1._dp - pcvs(jl))
            zburd(jl, jk) = zburd(jl, jk) - zburd(jl, jk) * zvolveg * &
                           (1._dp - pcvs(jl)) * dtime

            jk = klev
            zvola(jl, jk) = zburd(jl, jk) * zvolsoil * (1._dp - pcvs(jl))
            zburd(jl, jk) = zburd(jl, jk) - zburd(jl, jk) * zvolsoil * &
                           (1._dp - pcvs(jl)) * dtime

            ! distribution if snow is melting
            IF (pcvs(jl) .LT. zcvs_old(jl)) THEN
               jk = klev - 3

               zburd(jl, klev) = zburd(jl, klev) + zburd(jl, jk) * &
                              (1._dp - pcvs(jl) / zcvs_old(jl)) * &
                              (1._dp - pvgrat(jl))

               zburd(jl, klev - 1) = zburd(jl, klev - 1) + zburd(jl, jk) * &
                               (1._dp - pcvs(jl) / zcvs_old(jl)) * pvgrat(jl)

               zburd(jl, jk) = zburd(jl, jk) * pcvs(jl) / zcvs_old(jl)
            END IF

         ELSE ! not l_landsnow
            jk = klev - 1
            zxtems(jl) = zxtems(jl) + (zvolveg * zburd(jl, jk) + &
                         zvolsoil * zburd(jl, klev)) * M_air / molmass

            zvola(jl, jk) = zburd(jl, jk) * zvolveg
            zburd(jl, jk) = zburd(jl, jk) - zburd(jl, jk) * zvolveg * dtime

            jk = klev
            zvola(jl, jk) = zburd(jl, jk) * zvolsoil
            zburd(jl, jk) = zburd(jl, jk) - zburd(jl, jk) * zvolsoil * dtime
         END IF
      END IF

   END IF ! loland
END DO
      
!CALL end_message(modstr, 'LAND SURFACE VOLATILIZATIONS', substr)

END SUBROUTINE svoc_volatilizations

! *************************************************************************

SUBROUTINE svoc_read_nml_ctrl(status, iou)

!!! read CTRL namelist

USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

IMPLICIT NONE

! I/O
INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
INTEGER, INTENT(OUT) :: status ! error status

NAMELIST /CTRL/ L_GP, &
                L_LG, &
                l_svocpart, &
                l_svocvola, &
                l_glacier, &
                l_landsnow, &
                l_pahderiv

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_read_nml_ctrl'
LOGICAL :: lex   ! check if file exists
INTEGER :: fstat ! file status

!CALL start_message(modstr, 'READ NML CTRL', substr)

status = 1 ! error

CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
IF (.not. lex) RETURN  ! <modstr>.nml does not exist

READ(iou, NML = CTRL, IOSTAT = fstat)
CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
IF (fstat /= 0) RETURN  ! error while reading namelist

CALL read_nml_close(substr, iou, modstr)

status = 0  ! no error

!CALL end_message(modstr, 'READ NML CTRL', substr)

END SUBROUTINE svoc_read_nml_ctrl

! *************************************************************************
! PRIVATE SUBROUTINES
! *************************************************************************

SUBROUTINE svoc_bapo3_roi(press, &
                          relhum, &
                          temp, &
                          o3conc, &
                          rconst)

!!! on-particle bap oxidation with ozone based on ROI(T) scheme

USE messy_main_constants_mem, ONLY: R_gas, N_A

IMPLICIT NONE

! I/O
REAL(DP), INTENT(IN) :: press, relhum, temp
REAL(DP), INTENT(IN) :: o3conc
REAL(DP), INTENT(OUT) :: rconst

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_bapo3_roi'
INTEGER, PARAMETER :: nrharr = 3  !! no. relhum array
INTEGER, PARAMETER :: ntmarr = 13  !! no. temp array
REAL(DP), PARAMETER :: toppb = 1e9

REAL(DP), PARAMETER :: &
   midrh_arr(nrharr) = (/ 10._dp, 50._dp, 70._dp /)

REAL(DP), PARAMETER :: &
   midtm_arr(ntmarr) = (/ 313._dp, 308._dp, 303._dp, &
                          298._dp, 296._dp, 288._dp, &
                          283._dp, 278._dp, 273._dp, &
                          268._dp, 263._dp, 258._dp, 253._dp /)

REAL(DP), DIMENSION(nrharr, ntmarr) :: cbase_arr

DATA cbase_arr(1, 1: ntmarr) / 0.00016051_dp, 0.00011024_dp,  6.75e-05_dp, &
                                 2.30e-05_dp,  -1.90e-05_dp, -2.65e-05_dp, &
                                -4.28e-05_dp,  -1.12e-05_dp,  3.39e-06_dp, &
                                 0._dp, 0._dp, 0._dp, 0._dp /

DATA cbase_arr(2, 1: ntmarr) /  0.00016526_dp, 0.00011635_dp, 0.000082449_dp, &
                               0.000054404_dp, 0.000026441_dp,  2.202e-06_dp, &
                                 7.836e-07_dp,  -4.790e-07_dp, -7.019e-06_dp, &
                                 0._dp, 0._dp, 0._dp, 0._dp /
 
DATA cbase_arr(3, 1: ntmarr) / 0.00016692_dp, 0.00011854_dp, 0.000085186_dp, &
                               0.000059366_dp, 0.000048442_dp, 0.000025472_dp, &
                               0.000015771_dp, 9.493e-06_dp, 6.849e-06, &
                                5.311e-06_dp, 3.831e-06_dp, 8.928e-07_dp, &
                                1.157e-06_dp /

REAL(DP), DIMENSION(nrharr, ntmarr) :: cmax_arr

DATA cmax_arr(1, 1: ntmarr) / 0.017433_dp, 0.011941_dp, 0.0073069_dp, &
                              0.0030959_dp, 0.00092774_dp, 0.0003086_dp, &
                              0.00010305_dp, 3.49e-05_dp, 1.03e-05_dp, &
                              0._dp, 0._dp, 0._dp, 0._dp /

DATA cmax_arr(2, 1: ntmarr) / 0.017991_dp, 0.012676_dp, 0.0086796_dp, &
                              0.0057148_dp, 0.0030781_dp, 0.001096_dp, &
                              0.00025979_dp, 0.00007189_dp, 0.000021764_dp, &
                              0._dp, 0._dp, 0._dp, 0._dp /

DATA cmax_arr(3, 1: ntmarr) / 0.018187_dp, 0.012945_dp, 0.0090248_dp, &
                              0.0063313_dp, 0.0051402_dp, 0.0027209_dp, &
                              0.001676_dp, 0.0010024_dp, 0.00056362_dp, &
                              0.00031025_dp, 0.00014791_dp, 0.00010517_dp, &
                              0.000025839_dp /

REAL(DP), DIMENSION(nrharr, ntmarr) :: crate_arr

DATA crate_arr(1, 1: ntmarr) /  0.6666_dp, 0.66699_dp, 0.65303_dp, &
                               0.58645_dp, 0.44643_dp, 0.34062_dp, &
                               0.19055_dp, 0.14182_dp, 0.38165_dp, &
                               0._dp, 0._dp, 0._dp, 0._dp /

DATA crate_arr(2, 1: ntmarr) / 0.67213_dp, 0.67801_dp, 0.68829_dp, &
                               0.68755_dp, 0.61771_dp, 0.5436_dp, &
                               0.55929_dp, 0.56406_dp, 0.46381_dp, &
                               0._dp, 0._dp, 0._dp, 0._dp /

DATA crate_arr(3, 1: ntmarr)  / 0.67404_dp, 0.68188_dp, 0.69595_dp, &
                                0.70671_dp,  0.7002_dp, 0.71101_dp, &
                                0.70568_dp, 0.69809_dp, 0.70441_dp, &
                                0.71294_dp, 0.73128_dp, 0.58043_dp, &
                                0.67322_dp /
 
REAL(DP), DIMENSION(nrharr, ntmarr) :: cxhlf_arr

DATA cxhlf_arr(1, 1: ntmarr) / 1185.5_dp, 1206.1_dp, 1114.5_dp, &
                               702.11_dp, 159.16_dp, 83.811_dp, &
                                4.411_dp, 4.3287_dp, 884.64_dp, &
                               0._dp, 0._dp, 0._dp, 0._dp /

DATA cxhlf_arr(2, 1: ntmarr) / 1206.3_dp, 1238.6_dp, 1217.6_dp, &
                               1211.9_dp, 796.87_dp, 583.32_dp, &
                               74.982_dp, 13.122_dp, 1.052_dp, &
                               0._dp, 0._dp, 0._dp, 0._dp /

DATA cxhlf_arr(3, 1: ntmarr) / 1213.5_dp, 1250.8_dp, 1234.1_dp, &
                               1273.2_dp, 1251.8_dp, 1289.4_dp, &
                               1271.7_dp, 1243.1_dp, 1107.8_dp, &
                               1008.8_dp, 733.11_dp, 1882.3_dp, 334.89_dp /

REAL(DP) :: conv, o3ppb, rhact
REAL(DP), DIMENSION(nrharr) :: rhlow, rhupp
REAL(DP), DIMENSION(ntmarr) :: tmlow, tmupp
INTEGER :: tmpos, rhpos
INTEGER :: i, j

DO i = 1, nrharr
   IF (i > 1) THEN
      rhlow(i) = midrh_arr(i) - (midrh_arr(i) - midrh_arr(i - 1)) / 2
   ELSE
      rhlow(i) = 0._dp
   END IF

   IF (i < nrharr) THEN
      rhupp(i) = midrh_arr(i) + (midrh_arr(i + 1) - midrh_arr(i)) / 2
   ELSE
      rhupp(i) = 100._dp
   END IF
END DO

rhact = MAX(0._dp, MIN(99._dp, relhum))

DO j = ntmarr, 1, -1
   IF (j < ntmarr) THEN
      tmlow(j) = midtm_arr(j) - (midtm_arr(j) - midtm_arr(j + 1)) / 2
   ELSE
      tmlow(j) = 100._dp
   END IF

   IF (j > 1) THEN
      tmupp(j) = midtm_arr(j) + (midtm_arr(j - 1) - midtm_arr(j)) / 2
   ELSE
      tmupp(j) = 350._dp
   END IF
END DO

rhpos = 0
tmpos = 0

DO i = 1, nrharr
   IF (rhact >= rhlow(i) .AND. rhact < rhupp(i)) rhpos = i

   DO j = 1, ntmarr
      IF (temp >= tmlow(j) .AND. temp < tmupp(j)) tmpos = j
   END DO
END DO

! convert o3 conc from molec cm-3 to ppb
o3ppb = o3conc * toppb / (press / (R_gas * temp) * 1e-6_dp * N_A)

! decay rate [s-1]
IF (o3ppb > 0.4_dp) THEN
   rconst = cbase_arr(rhpos, tmpos) + &
         (cmax_arr(rhpos, tmpos) - cbase_arr(rhpos, tmpos)) / &
         (1 + (cxhlf_arr(rhpos, tmpos) / o3ppb) ** crate_arr(rhpos, tmpos))
ELSE
   rconst = 0._dp
END IF

END SUBROUTINE svoc_bapo3_roi

! *************************************************************************

SUBROUTINE svoc_part_lohlam(rkoab, rkoam, &
                        !! if ksa is from ksw and kaw (obsolete)
                            !ha, hb, &
                            !rlogksw, &
                        !! if ksa is from vapor pressure
                            rvapp, rhvap, &
                            tair, &
                            concom, &
                            concbc, &
                            kp, &
                            theta)

!!! gas-particle partitioning considering absorption into organic matter
!!! and adsoprtion to black carbon (lohmann-lammel, 2004)

USE messy_main_constants_mem, ONLY: R_gas

IMPLICIT NONE

! I/O
REAL(DP), INTENT(IN) :: rkoab, rkoam
!REAL(DP), INTENT(IN) :: ha, hb
!REAL(DP), INTENT(IN) :: rlogksw
REAL(DP), INTENT(IN) :: rvapp, rhvap
REAL(DP), INTENT(IN) :: tair
REAL(DP), INTENT(IN) :: concom, concbc
REAL(DP), INTENT(OUT) :: kp, theta

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_part_lohlam'
REAL(DP), PARAMETER :: ssaec = 18.21_dp  ! m2 g-1 (same as in pplfer)
REAL(DP), PARAMETER :: tstd = 298.15_dp
REAL(DP), PARAMETER :: kg2ug = 1e+9_dp
REAL(DP), PARAMETER :: cnv = 101.325_dp  ! to convert M atm-1 to mol m-3 Pa-1
REAL(DP), PARAMETER :: tvap = 298.15_dp
!REAL(DP) :: henry, kaw
REAL(DP) :: zvapp
REAL(DP) :: koa, ksa

! temperature dependency for koa [-]
koa = 10._dp ** (rkoab + rkoam / tair)

! ksa [L kg -1]
! method1: as ratio of ksw and kaw
!henry = ha * EXP(hb * ((1._dp / tair) - (1._dp / tstd))) / cnv
!kaw = 1._dp / (R_gas * tair * henry)  ! air-liquid partitioning coef. [-]
!ksa = 10._dp ** (rlogksw) / kaw

! method2: as a function of vapor pressure (van Noort, 2003)
! liquid vapor pressure at environment temperature is interpolated
! from vapor pressure at 298 kelvin
zvapp = rvapp * EXP(-rhvap / R_gas * ((1._dp / tair) - (1._dp / tvap)))
ksa = 10._dp ** (-0.85_dp * LOG10(zvapp) + 8.94_dp - LOG10(998._dp / ssaec))

! partitioning coefficient (sehili and lammel, 2007; eq. 5)
kp = koa * 1.22e-12_dp * concom + ksa * 1.e-12_dp * concbc

! particle mass fraction
theta = kp / (1._dp + kp)

END SUBROUTINE svoc_part_lohlam

! *************************************************************************

SUBROUTINE svoc_part_finizio(rkoab, rkoam, &
                             tair, &
                             conctsp, &
                             kp, &
                             theta)

!!! partitioning through absorption of gas-phase chemical into
!!! an organic film; partitioning coefficient [mg ug-1] is derived
!!! from plots of log kp vs log koa (finizio, 1997)

IMPLICIT NONE

! I/O
REAL(DP), INTENT(IN) :: rkoab, rkoam
REAL(DP), INTENT(IN) :: tair, conctsp
REAL(DP), INTENT(OUT) :: kp, theta

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_part_finizio'
REAL(DP) :: koa

! temperature dependency for koa [-]
koa = 10._dp ** (rkoab + rkoam / tair)

! partitioning coefficient [m3 ug-1]
kp = 10._dp ** (.79_dp * LOG10(koa) - 10.01)

! particle mass fraction
theta = conctsp * kp / (kp * conctsp + 1._dp)

END SUBROUTINE svoc_part_finizio

! *************************************************************************

SUBROUTINE svoc_part_junge(rvapp, rhvap, &
                           tair, &
                           rarea, &
                           theta)

!!! partitioning process considers only adsorption onto aerosol surfaces
!!! (junge 1977; pankow 1987). fraction of particles depends on vapor
!!! pressure and total aerosol surface area per unit volume

USE messy_main_constants_mem, ONLY: R_gas

IMPLICIT NONE

! I/O
REAL(DP), INTENT(IN) :: rvapp, rhvap
REAL(DP), INTENT(IN) :: tair
REAL(DP), INTENT(IN) :: rarea
REAL(DP), INTENT(OUT) :: theta

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_part_junge'
REAL(DP), PARAMETER :: tvap = 298.15_dp
REAL(DP) :: zvapp, zvaps

! sub-cooled liquid vapor pressure at environment temperature is
! interpolated from liquid vapor pressure at 298 kelvin
zvapp = rvapp * EXP(-rhvap / R_gas * ((1._dp / tair) - (1._dp / tvap)))

! particle mass fraction
IF (rarea .GE. 1.e-20_dp) THEN
   ! a constant adopted from pankow (1987) for pahs (1.3 torr cm or
   ! 173.32 Pa cm3 cm-2) overestimates theta for semi-volatile compounds
   ! (tested for anthracene), hence a lower constant value originally
   ! applied in junge (17.2 Pa cm3 cm-2) is used
   theta = 17.2_dp * rarea / (17.2_dp * rarea + zvapp)
ELSE
   theta = 0._dp
END IF

END SUBROUTINE svoc_part_junge

! *************************************************************************

SUBROUTINE svoc_part_pplfers(rhvap, &
                             rhabc, & !! mz_jw_20170220
                             rpplfer, &
                             tair, &
                             rhum, &
                             concbc, &
                             concwsom, &
                             concwiom, &
                             concso4, &
                             concss, &
                             kp, &
                             theta)

!!! partitioning from absorption and adsorption contributions
!!! based on multiple linear regression models
!!! given: substance solute descriptors and system parameters

USE messy_main_constants_mem, ONLY: R_gas

IMPLICIT NONE

! I/O
REAL(DP), INTENT(IN) :: rhvap, rhabc
REAL(DP), DIMENSION(6), INTENT(IN) :: rpplfer
REAL(DP), INTENT(IN) :: tair, rhum
REAL(DP), INTENT(IN) :: concbc, concwsom, concwiom
REAL(DP), INTENT(IN) :: concso4, concss
REAL(DP), INTENT(OUT) :: kp, theta

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_part_pplfers'
REAL(DP), PARAMETER :: t15 = 288.15_dp
REAL(DP), PARAMETER :: t25 = 298.15_dp
REAL(DP), PARAMETER :: midrh_arr(3) = (/20._dp, 40._dp, 60._dp/)
! specific surface area [m2 g-1]
REAL(DP), PARAMETER :: ssa_ec = 18.21_dp    !!! elemental carbon
REAL(DP), PARAMETER :: ssa_nh42so4 = 0.1_dp !!! ammonium sulfate
REAL(DP), PARAMETER :: ssa_nacl = 0.1_dp  !!! sodium chloride
! density [g m-3]
REAL(DP), PARAMETER :: rho_dmso = 1.1e6_dp !!! dmso
REAL(DP), PARAMETER :: rho_hexd = 0.77e6_dp !!! hexadecane
! molecular weight [g mol-1]
REAL(DP), PARAMETER :: mw_nh42so4 = 132.14_dp
REAL(DP), PARAMETER :: mw_so4mm = 96.07_dp
! system parameters for soot-air partition (adsorption) coeff.
REAL(DP), PARAMETER :: tsoot = t15
REAL(DP), PARAMETER :: se1 = 0._dp
REAL(DP), PARAMETER :: ss1 = 0._dp
REAL(DP), PARAMETER :: sa1 = 2.7_dp
REAL(DP), PARAMETER :: sb1 = 2.45_dp
REAL(DP), PARAMETER :: sv1 = 0._dp
REAL(DP), PARAMETER :: sl1 = 1.09_dp
REAL(DP), PARAMETER :: sc1 = -8.47_dp
! system parameters for ammonium sulfate-air partition (adsorption) coeff.
REAL(DP), PARAMETER :: tso4 = t15
REAL(DP), PARAMETER :: se2(3) = (/0._dp, 0._dp, 0._dp/)
REAL(DP), PARAMETER :: ss2(3) = (/0._dp, 0._dp, 0._dp /)
REAL(DP), PARAMETER :: sa2(3) = (/2.46_dp, 2.46_dp, 2.13_dp/)
REAL(DP), PARAMETER :: sb2(3) = (/5.23_dp, 5.23_dp, 5.34_dp/)
REAL(DP), PARAMETER :: sv2(3) = (/0._dp, 0._dp, 0._dp/)
REAL(DP), PARAMETER :: sl2(3) = (/0.90_dp, 0.89_dp, 0.88_dp/)
REAL(DP), PARAMETER :: sc2(3) = (/-8.47_dp, -8.47_dp, -8.47_dp/)
! system parameters for sodium chloride-air partition (adsorption) coeff.
REAL(DP), PARAMETER :: tss = t15
REAL(DP), PARAMETER :: se3(3) = (/0._dp, 0._dp, 0._dp/)
REAL(DP), PARAMETER :: ss3(3) = (/0._dp, 0._dp, 0._dp/)
REAL(DP), PARAMETER :: sa3(3) = (/3.12_dp, 2.94_dp, 2.86_dp/)
REAL(DP), PARAMETER :: sb3(3) = (/4.77_dp, 4.82_dp, 4.82_dp/)
REAL(DP), PARAMETER :: sv3(3) = (/0._dp, 0._dp, 0._dp/)
REAL(DP), PARAMETER :: sl3(3) = (/0.87_dp, 0.86_dp, 0.84_dp/)
REAL(DP), PARAMETER :: sc3(3) = (/-8.47_dp, -8.47_dp, -8.47_dp/)
! system parameters for dry dmso-air partition (absorption) coeff.
REAL(DP), PARAMETER :: tdmso = t25
REAL(DP), PARAMETER :: se4 = -0.223_dp
REAL(DP), PARAMETER :: ss4 = 2.903_dp
REAL(DP), PARAMETER :: sa4 = 5.036_dp
REAL(DP), PARAMETER :: sb4 = 0._dp
REAL(DP), PARAMETER :: sv4 = 0._dp
REAL(DP), PARAMETER :: sl4 = 0.719_dp
REAL(DP), PARAMETER :: sc4 = -0.556_dp
! system parameters for dry dmso-air phase transfer enthalpy
REAL(DP), PARAMETER :: see4 = 0._dp
REAL(DP), PARAMETER :: sss4 = -19.041_dp
REAL(DP), PARAMETER :: saa4 = -47.799_dp
REAL(DP), PARAMETER :: sbb4 = -5.521_dp
REAL(DP), PARAMETER :: svv4 = -0.746_dp
REAL(DP), PARAMETER :: sll4 = -6.189_dp
REAL(DP), PARAMETER :: scc4 = -2.390_dp
! system parameters for polyurethane ether-air partition (absorption) coeff.
REAL(DP), PARAMETER :: tpu = t15
REAL(DP), PARAMETER :: se5 = 0._dp 
REAL(DP), PARAMETER :: ss5 = 1.69_dp 
REAL(DP), PARAMETER :: sa5 = 3.66_dp
REAL(DP), PARAMETER :: sb5 = 0._dp
REAL(DP), PARAMETER :: sv5 = 0.36_dp
REAL(DP), PARAMETER :: sl5 = 0.71_dp
REAL(DP), PARAMETER :: sc5 = -0.15_dp
! system parameters for polyurethane ether-air phase transfer enthalpy
REAL(DP), PARAMETER :: see5 = 0._dp
REAL(DP), PARAMETER :: sss5 = -17.6_dp
REAL(DP), PARAMETER :: saa5 = -46.6_dp
REAL(DP), PARAMETER :: sbb5 = 0._dp
REAL(DP), PARAMETER :: svv5 = -12.8_dp
REAL(DP), PARAMETER :: sll5 = -4.3_dp
REAL(DP), PARAMETER :: scc5 = 2.7_dp
! system parameters for hexadecane-air partition (absorption) coeff.
REAL(DP), PARAMETER :: thd = t25
REAL(DP), PARAMETER :: se6 = 0._dp
REAL(DP), PARAMETER :: ss6 = 0._dp
REAL(DP), PARAMETER :: sa6 = 0._dp
REAL(DP), PARAMETER :: sb6 = 0._dp
REAL(DP), PARAMETER :: sv6 = 0._dp
REAL(DP), PARAMETER :: sl6 = 1._dp
REAL(DP), PARAMETER :: sc6 = 0._dp
! system parameters for hexadecane-air phase transfer enthalphy
REAL(DP), PARAMETER :: see6 = 0._dp
REAL(DP), PARAMETER :: sss6 = -3.20_dp
REAL(DP), PARAMETER :: saa6 = -1.51_dp
REAL(DP), PARAMETER :: sbb6 = 2.58_dp
REAL(DP), PARAMETER :: svv6 = -10.86_dp
REAL(DP), PARAMETER :: sll6 = -6.79_dp
REAL(DP), PARAMETER :: scc6 = -2.97_dp

! local
INTEGER :: i, bin
REAL(DP), DIMENSION(3) :: rhlow, rhupp
REAL(DP) :: rbige, rbigs, rbiga, rbigb, rbigv, rbigl
REAL(DP) :: kpec, kpso4, kpss, kpdmso, kppu, kphd
REAL(DP) :: conc_nh42so4
REAL(DP) :: rhec, rhso4, rhss, rhdmso, rhpu, rhhd
REAL(DP) :: kpec_t, kpso4_t, kpss_t, kpdmso_t, kppu_t, kphd_t
REAL(DP) :: rhact

! solute descriptors
rbige = rpplfer(1)
rbigs = rpplfer(2)
rbiga = rpplfer(3)
rbigb = rpplfer(4)
rbigv = rpplfer(5)
rbigl = rpplfer(6)

! find rh bin
DO i = 1, 3
   IF (i > 1) THEN
      rhlow(i) = midrh_arr(i) - (midrh_arr(i) - midrh_arr(i - 1)) / 2
   ELSE
      rhlow(i) = 0._dp
   END IF

   IF (i < 3) THEN
      rhupp(i) = midrh_arr(i) + (midrh_arr(i + 1) - midrh_arr(i)) / 2
   ELSE
      rhupp(i) = 100._dp
   END IF
END DO

rhact = MAX(0._dp, MIN(99._dp, rhum))

bin = 0

DO i = 1, 3
   IF (rhact >= rhlow(i) .AND. rhact < rhupp(i)) THEN
      bin = i
      EXIT
   END IF
END DO

! partitioning due to adsorption [m3(air) m-2(surface)] in:
! (1) soot surface-air partition system
kpec = 10._dp ** (sl1 * rbigl + sa1 * rbiga + sb1 * rbigb + sc1)

! (2) ammonium sulfate surface-air
kpso4 = 10._dp ** (sl2(bin) * rbigl + sa2(bin) * rbiga + &
   sb2(bin) * rbigb + sc2(bin))

! (3) sodium chloride surface-air
kpss = 10._dp ** (sl3(bin) * rbigl + sa3(bin) * rbiga + &
   sb3(bin) * rbigb + sc3(bin))

! partitioning due to absorption to organic fraction in
! (1) dmso-air system [L(air) L-1(solvent)]
kpdmso = 10._dp ** (sl4 * rbigl + ss4 * rbigs + sa4 * rbiga + &
   sb4 * rbigb + se4 * rbige + sc4)

! (2) polyurethane-air system [L(air) kg-1(pu)]
kppu = 10._dp ** (sl5 * rbigl + ss5 * rbigs + sa5 * rbiga + &
   sb5 * rbigb + sv5 * rbigv + sc5)

! (3) hexadecane-air system [L(air) L-1(solvent)]
kphd = 10._dp ** (sl6 * rbigl + ss6 * rbigs + sa6 * rbiga + &
   sb6 * rbigb + se6 * rbige + sc6)

! calculate enthalpy of phase transfer [J mol-1]
!! mz_jw_20170220
if (rhabc > 0._dp) then
   rhec = -rhabc
else
   rhec = -rhvap
endif
!! mz_jw_20170220
rhso4 = (-10.2_dp * LOG10(kpso4) - 89.6_dp) * 1.e3_dp
rhss = (-10.2_dp * LOG10(kpss) - 89.6_dp) * 1.e3_dp
rhdmso = (sll4 * rbigl + sss4 * rbigs + saa4 * rbiga + &
      sbb4 * rbigb + svv4 * rbigv + scc4) * 1.e3_dp
rhpu = (sll5 * rbigl + sss5 * rbigs + saa5 * rbiga + &
      svv5 * rbigv + scc5) * 1.e3_dp
rhhd = (sll6 * rbigl + sss6 * rbigs + saa6 * rbiga + &
   sbb6 * rbigb + svv6 * rbigv + scc6) * 1.e3_dp

! temperature correction of indiv. partitioning coefficients
kpec_t = kpec * EXP(-rhec / R_gas * ((1._dp / tair) - (1._dp / tsoot)))
kpso4_t = kpso4 * EXP(-rhso4 / R_gas * ((1._dp / tair) - (1._dp / tso4)))
kpss_t =  kpss * EXP(-rhss / R_gas * ((1._dp / tair) - (1._dp / tss)))
kpdmso_t = kpdmso * EXP(-rhdmso / R_gas * ((1._dp / tair) - (1._dp / tdmso)))
kppu_t = kppu * EXP(-rhpu / R_gas * ((1._dp / tair) - (1._dp / tpu)))
kphd_t = kphd * EXP(-rhhd / R_gas * ((1._dp / tair) - (1._dp / thd)))

! total partitoning coefficient [-]
conc_nh42so4 = concso4 * mw_nh42so4 / mw_so4mm

kp = (kpec_t * ssa_ec * concbc * 1.e-6_dp + &
     kpso4_t * ssa_nh42so4 * conc_nh42so4 * 1.e-6_dp + &
     kpss_t * ssa_nacl * concss * 1.e-6_dp) + &
     (kpdmso_t * concwsom / rho_dmso * 1.e-6_dp + &
     kppu_t * 0.20_dp * concwiom * 1.e-12_dp + &
     kphd_t * 0.80_dp * concwiom / rho_hexd * 1.e-6_dp)

! particle mass fraction
theta = kp / (1._dp + kp)

END SUBROUTINE svoc_part_pplfers

! *************************************************************************

SUBROUTINE svoc_vola_glac_snow(loglac, dtime, &
                               molmass, &
                               rho_sn, &
                               ha, hb, &
                               rhsub, rhsol, &
                               rwsol, &
                               depth, &
                               burd, &
                               tsurf, &
                               zvolf)

!!! volatilization from glacier and snow (wania, 1997)

USE messy_main_constants_mem, ONLY: R_gas, rho_H2O

IMPLICIT NONE

! I/O
LOGICAL, INTENT(IN) :: loglac
REAL(DP), INTENT(IN) :: dtime
REAL(DP), INTENT(IN) :: molmass, rho_sn
REAL(DP), INTENT(IN) :: ha, hb
REAL(DP), INTENT(IN) :: rhsub, rhsol, rwsol
REAL(DP), INTENT(IN) :: depth, burd, tsurf
REAL(DP), INTENT(OUT) :: zvolf

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_vola_glac_snow'
REAL(DP), PARAMETER :: tstd = 298.15_dp
REAL(DP), PARAMETER :: tvap = 298.15_dp 
REAL(DP), PARAMETER :: t20 = 293.15_dp
REAL(DP), PARAMETER :: cnv = 101.325_dp ! to convert M atm-1 to mol m-3 Pa-1
REAL(DP) :: zsnow_h, moldifa, moldifw
REAL(DP) :: aden, bden, zu5, zu6, zu7
REAL(DP) :: henry, henry_t20, zwsol_t20, lkia_t20
REAL(DP) :: rhads, lkia, pkia, pzi
REAL(DP) :: asa, vsa, vsl
REAL(DP) :: pdvl, pdval, fug, pMTC

vsa = 0.3_dp ! volume fraction of air in snowpack
vsl = 0.1_dp ! volume fraction of liquid in snowpack

zsnow_h = depth * rho_H2O / rho_sn  ! snow cover depth

IF (loglac) THEN
   IF (zsnow_h .LT. 0.1_dp) THEN
      zsnow_h = 0.1_dp ! min snow cover for calculation
   ELSE IF (zsnow_h .GT. 0.5_dp) THEN
      zsnow_h = 0.5_dp ! max before transfer to deep snow/ice
   END IF
ELSE ! landsnow
   IF (zsnow_h .LT. 0.05_dp) &
      zsnow_h = 0.05_dp ! min snow cover
END IF

! molecular diffusivity in air is  calculated from the molar mass
! (schwarzenbach et al., 2nd ed., 2003, fig 18.9) [cm2 s-1],
! while diffusivity in water is assumed to be 1e4 less than that in air
moldifa = (1.55_dp / molmass ** 0.65_dp) / 1.e4_dp ! [m2 s-1]
moldifw = 1.e-4_dp * moldifa ! [m2 s-1]

aden = (vsa + vsl) ** 2
bden = LOG(2._dp) * zsnow_h
zu5 = moldifw * (vsl ** (10._dp / 3._dp) / aden) / bden
zu6 = moldifa * (vsa ** (10._dp / 3._dp) / aden) / bden
zu7 = 5._dp / 3600._dp ! [m s-1] (wania, 1997)

! ice surface-air partition coefficients (kia)
! kia at 20-deg derived from kwa (RTH) and water solubility [mol m-3]
! temperature dependence of the solubility is described by van't Hoff eq.
henry_t20 = ha * EXP(hb * ((1._dp / t20) - (1._dp / tstd))) / cnv

zwsol_t20 = EXP(-rhsol / R_gas * &
                ((1._dp / t20) - (1._dp / tstd)) + LOG(rwsol))

lkia_t20 = -0.769_dp * LOG10(zwsol_t20) - 5.97_dp + &
      LOG10(R_gas * t20 * henry_t20) ! log(kia) at 20-deg celcius

! extrapolate kia to other temperature
! rhads = enthalphy of adsorption, estimated from enthaly of sublimation
rhads = 0.878_dp * rhsub
lkia = lkia_t20 + rhads / (2.303_dp * R_gas) * (1._dp / tsurf - 1._dp / t20)
pkia = 10._dp**(lkia)

! capacity term [mol m-2 Pa-1] for ice-air interface
pzi = pkia / (R_gas * tsurf)

! diffusion value for chemical loss due to volatilization to
! the atmosphere [mol Pa-1 s-1] based on eq 1 in wania (1997)
! area covered by snowpack is assumed to be 1 m2
henry = ha * EXP(hb * ((1._dp / tsurf) - (1._dp / tstd))) / cnv

pdvl = 1._dp / (zu7 / (R_gas * tsurf))  + &
       1._dp / (zu5 * henry + zu6 / (R_gas * tsurf))

pdval = 1._dp / pdvl

! fugacity follows Galy and Wania (2004) initially based on Wania (1997)
!
! area of snow-air interface per volume snowpack (asa, in m-1)
! asa = specific surface area * density of snowmelt water
! specific surface area = 0.1 m2 g-1 during snow accumulation period,
!   and 0.1 to 0.01 m2 g-1 with linear decrease during snowmelt period;
!   here, assumed to be 0.05 m2 g-1
! density of snowmelt water = 0.433 g cm-3 = 4.33 x 1e5 g m-3
! hence, asa = 4.33e5_dp * 0.05_dp
! but here follows mpi-mctm:
! specific surface area = 0.025 m2 g-1; density of snowpack = 7.e5 g m-3
asa = 7.e5_dp * 0.025_dp
 
fug = burd / (vsa * zsnow_h / (R_gas * tsurf) + &
      henry * vsl * zsnow_h + pzi * asa * zsnow_h) ! [Pa kg mol-1 m-2]
pMTC = pdval * fug ! [kg m-2 s-1]

! if all burden are readily volatilized
IF (burd .LT. pMTC * dtime) THEN
   pMTC = burd / dtime
END IF

IF (burd .GT. 0._dp) THEN
   zvolf = pMTC / burd ! s-1
ELSE
   zvolf = 0._dp
END IF

END SUBROUTINE svoc_vola_glac_snow

! *************************************************************************

SUBROUTINE svoc_vola_veg(rhvap, rvapp, &
                         tsurf, &
                         zvolf)

!!! volatilization from vegetation (smit 1998)

USE messy_main_constants_mem, ONLY: R_gas

IMPLICIT NONE

! I/O
REAL(DP), INTENT(IN) :: rhvap, rvapp
REAL(DP), INTENT(IN) :: tsurf
REAL(DP), INTENT(OUT) :: zvolf

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_vola_veg'
REAL(DP), PARAMETER :: tvap = 298.15_dp
REAL(DP) :: ztdecv, zvapp
REAL(DP) :: zcvv

! conversion factors for decay in in vegetation (7 days) [s-1]
ztdecv = 1._dp/7._dp/24._dp/3600._dp

! vapor pressure [Pa]
zvapp = rvapp * EXP(-rhvap / R_gas * &
        ((1._dp / tsurf) - (1._dp / tvap)))

IF (zvapp >= 1.03e-2_dp) zvapp = 1.03e-2_dp
IF (zvapp <= 5.2655e-7_dp) zvapp = 5.2655e-7_dp

! accumulated volatilization during 7 days after exposure [% of stored]
zcvv = 10._dp ** (1.528_dp + 46.6e-2_dp * LOG10(1.e3_dp * zvapp))

! volatilization rate [s-1]
zvolf = (zcvv / 100._dp) * ztdecv

END SUBROUTINE svoc_vola_veg

! *************************************************************************

SUBROUTINE svoc_vola_soil_jury(dtime, &
                               pforest, prhos, &
                               zporefr, zkaw, zpksl, &
                               zvolg, zvolw, &
                               zvolf)

!!! volatilization from soil (jury et al., 1983 and 1990)
!!!
!!! assumptions (as applied by cousins et al., 1997):
!!! - homogeneously contaminated soil layers
!!! - thickness of the stagnant boundary layer above soil surface is zero
!!! - no water evaporation

IMPLICIT NONE 

! I/O
REAL(DP), INTENT(IN) :: dtime
REAL(DP), INTENT(IN) :: pforest, prhos
REAL(DP), INTENT(IN) :: zporefr, zkaw, zpksl
REAL(DP), INTENT(IN) :: zvolg, zvolw
REAL(DP), INTENT(OUT) :: zvolf

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_vola_soil_jury'
REAL(DP) :: hsoil
REAL(DP) :: B_air, B_wat
REAL(DP) :: D_gas, D_wat, D_eff

! hsoil: average soil depth (m)
! hsoil = 0.1 m, forest (CoZmo-POP)
! hsoil = 0.2 m, agriculture (CoZmo-POP)
! make use of pforest to weight hsoil within gcell
hsoil = pforest * 0.1 + (1 - pforest) * 0.2
IF (hsoil .EQ. 0._dp) hsoil = 0.2_dp

! diffusion coefficients
! constants below are representative for most pesticides but it
! varies only slightly for other similar molecular weight organic
! compounds (jury, 1983)
B_air = 0.018_dp / 3600._dp ! [m2 s-1]; Boynton & Brattain, 1929
B_wat = 0.0000018_dp / 3600._dp ! [m2 s-1]; Bruins, 1929

! B_air and B_wat are corrected for tortuosity according to
! the millington-quirk model
D_gas = (B_air * zvolg ** (10._dp / 3._dp) / zporefr ** 2._dp)
D_wat = (B_wat * zvolw ** (10._dp / 3._dp) / zporefr ** 2._dp)

! effective diffusion coefficient (jury et al., 1990)
D_eff = (D_gas * zkaw + D_wat) / &
        (prhos * zpksl + zvolw + zvolg * zkaw) 

! volatilization rate [s-1] (jury et al., 1990)
zvolf  = SQRT(D_eff / (pi * dtime)) * &
         (1 - EXP(-(hsoil ** 2._dp) / (4._dp * D_eff * dtime))) / hsoil

END SUBROUTINE svoc_vola_soil_jury

! *************************************************************************

SUBROUTINE svoc_vola_soil_smit(dtime, &
                               zkaw, zpksl, &
                               zvolg, zvolw, &
                               pws, pwsmx, &
                               prhos, &
                               zvolf)

!!! volatilization from soil (smit 1997)
!!! 
!!! empiric relationships for normal-moist and dry soil as
!!! functions of CV (cumulative volatilization) and fraction
!!! of the substance in gas phase
!!!
!!! note: formulation here follows semeena et al. (2006)

USE messy_main_constants_mem, ONLY: R_gas

IMPLICIT NONE

! I/O
REAL(DP), INTENT(IN) :: dtime
REAL(DP), INTENT(IN) :: prhos
REAL(DP), INTENT(IN) :: zkaw, zpksl
REAL(DP), INTENT(IN) :: zvolg, zvolw
REAL(DP), INTENT(IN) :: pws, pwsmx
REAL(DP), INTENT(OUT) :: zvolf

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_vola_soil_smit'
REAL(DP) :: ztdecs, zkwa, zfsvocg
REAL(DP) :: zcvs

! conversion factors for decay in soil (3 weeks) [s-1]
ztdecs = 1._dp/21._dp/24._dp/3600._dp

! fraction of svoc in gas phase
zkwa = 1._dp / zkaw
zfsvocg = zvolg / (zvolg + zkwa * (zvolw + prhos * zpksl))
IF (zfsvocg > 1._dp) zfsvocg = 1._dp

! accumulative volatilization during 3 weeks after exposure [% of stored]
! threshold of normal-moist vs dry (@10% soil water of field capacity)
IF (pws / pwsmx > 1.e-1_dp) THEN  ! normal to moist soil
   IF (zfsvocg <= 6.3347e-9_dp) zfsvocg = 6.3347e-9_dp
   zcvs = 71.9_dp + 11.6_dp * LOG10(100._dp * zfsvocg)
ELSE ! dry
   IF (zfsvocg <= 1.9953e-7_dp) zfsvocg = 1.9953e-7_dp
   zcvs = 42.3_dp + 9.0_dp * LOG10(100._dp * zfsvocg)
END IF

! volatilization rate [s-1]
zvolf = (zcvs / 100._dp) * ztdecs

END SUBROUTINE svoc_vola_soil_smit

! *************************************************************************

END MODULE messy_svoc
