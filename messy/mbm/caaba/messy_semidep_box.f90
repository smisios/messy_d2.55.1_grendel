! Time-stamp: <2019-10-16 15:14:23 sander>

! SEMIDEP = Simplified EMIssion and DEPosition

! Authors:
! Rolf Sander,    MPICH, Mainz, 2003-2007
! Hella Riede,    MPICH, Mainz, 2007

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************

MODULE messy_semidep_box

  USE messy_mecca_kpp            ! dp, ind_*
  USE messy_main_constants_mem,  ONLY: OneDay, N_A
  USE caaba_mem,                 ONLY: model_time, zmbl, timesteplen, &
#ifdef __GFORTRAN__
                                      C_ => &
#endif
                                                                      C
#ifdef __GFORTRAN__    /* newer gfortran versions have a run-time bug of    */
#define C C_           /* module-level import of C(:), so it has to be      */
#define c C_           /* renamed into something else, e.g. C_(:)           */
#endif

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: semidep_physc

CONTAINS

  !***************************************************************************

  SUBROUTINE semidep_physc

    CALL emission ! emission
    CALL drydep   ! dry deposition

  END SUBROUTINE semidep_physc

  !***************************************************************************

  SUBROUTINE emission

    ! Emission and deposition are currently calculated with Euler forward.
    ! This should be changed if numerical problems occur.

    USE messy_main_constants_mem, ONLY: MH, MC ! for OOMPH only
    USE caaba_mem,                ONLY: emission_scenario, cair,         &
                                        l_injectNOx, t_NOxon, t_NOxoff,  &
                                        model_start_day
    IMPLICIT NONE
    REAL(dp) :: fct, fct2, fct3
    LOGICAL, SAVE :: O3_injected = .FALSE.

    !-------------------------------------------------------------------------

    ! activate next line to inject a chemical species into the system:
    ! CALL injection

    ! emission rates [molec/(cm2*s)] conv. to [molec/cm3]
    fct = timesteplen / (100. * zmbl)
    ! conversion 1 nmol/m2/day = 7E5 mcl/cm2/s:
    fct2 = fct * 1E-9 * N_A * 1E-4 / OneDay
    ! conversion from ng/(m2*s²)to g/(cm2*s):
    fct3 = 1.E-4 * 1.E-9 ! cm2/m2 * g/ng

    SELECT CASE (TRIM(emission_scenario))
    CASE ('CUMULUS')
      CALL emission_cumulus
    CASE ('FF_ARCTIC','FF_ANTARCTIC')
      CALL emission_ff
    CASE ('OOMPH')
      CALL emission_oomph
    CASE ('LAB')
      CALL emission_lab
    CASE ('LAB_C15')
      CALL emission_lab_c15
    CASE ('MBL')
      CALL emission_mbl
    CASE ('MOM')
      CALL emission_mom
    CASE ('ISO')
      CALL emission_iso
    CASE ('TAG')
      CALL emission_tag
    CASE ('')
      CALL emission_default
    CASE DEFAULT
      PRINT *, 'ERROR, emission_scenario '//TRIM(emission_scenario)// &
        ' is not defined'
      STOP 1
    END SELECT

    ! additional NO(x) EMISSIONS, l_injectNOx in namelist
    IF (l_injectNOx .AND. &
      (model_time >= (t_NOxon+model_start_day)*OneDay)   .AND. &
      (model_time <= (t_NOxoff+model_start_day)*OneDay)) THEN
      PRINT *, 'NO emissions'
      IF (ind_NO      /= 0) c(ind_NO)     = c(ind_NO)  + 1.0E11 * fct
    ENDIF

    !-------------------------------------------------------------------------

  CONTAINS

    !-------------------------------------------------------------------------

    SUBROUTINE emission_default
      IF (ind_O3     /= 0) c(ind_O3)     = c(ind_O3)     + 5.E10 * fct ! ref0203, p. 6699
    END SUBROUTINE emission_default

    !-------------------------------------------------------------------------

    SUBROUTINE emission_mbl
      IF (ind_O3     /= 0) c(ind_O3)     = c(ind_O3)     + 5.E10 * fct ! ref0203, p. 6699
      IF (ind_NO     /= 0) c(ind_NO)     = c(ind_NO)     + 1.5E9 * fct
      IF (ind_NH3    /= 0) c(ind_NH3)    = c(ind_NH3)    + 1.0E9 * fct
      IF (ind_DMS    /= 0) c(ind_DMS)    = c(ind_DMS)    + 2.0E9 * fct
      IF (ind_CHBr3  /= 0) c(ind_CHBr3)  = c(ind_CHBr3)  + 221.4 * fct2 ! Lucy
      IF (ind_CH3I   /= 0) c(ind_CH3I)   = c(ind_CH3I)   + 6.0E6 * fct ! ref0897
      IF (ind_C3H7I  /= 0) c(ind_C3H7I)  = c(ind_C3H7I)  + 1.0E7 * fct ! ref0897
      ! emissions of short-lived iodocarbons are scaled by a factor of 10 to
      ! indicate that they are not evenly distributed throughout the whole mbl
      IF (ind_CH2ClI /= 0) c(ind_CH2ClI) = c(ind_CH2ClI) + 185. *10.*fct2 !Lucy
      IF (ind_CH2I2  /= 0) c(ind_CH2I2)  = c(ind_CH2I2)  + 54.7 *10.*fct2 !Lucy
      !IF (ind_CH2ClI /= 0) c(ind_CH2ClI) = c(ind_CH2ClI) + 2.0E7 * fct ! ref0897
      !IF (ind_CH2I2  /= 0) c(ind_CH2I2)  = c(ind_CH2I2)  + 3.0E7 * fct ! ref0897
    END SUBROUTINE emission_mbl

    !-------------------------------------------------------------------------

    SUBROUTINE emission_lab
      ! O3 injection at 12:26:54 - 11:23:39 = 1 h, 3 min, 15 s after model_start (CIMS)
      IF ((.NOT.O3_injected).AND. &
        !(model_time >= model_start_day*OneDay+1.*3600.+3.*60.+15.)) THEN  
        (model_time >= model_start_day*OneDay+10020.)) THEN  !02:47 Start at 9619 20171106 Fe
        IF (ind_O3 /= 0) c(ind_O3) = 620E-09 * cair
        O3_injected = .TRUE.
      ENDIF
    END SUBROUTINE emission_lab

    !-------------------------------------------------------------------------

    SUBROUTINE emission_lab_c15

      USE caaba_mem, ONLY: temp, press
      USE messy_main_constants_mem, ONLY: R_gas
      REAL :: flow_1, flow_2, flow_3, flow_4, flow_t
      REAL :: flow_O2, flow_N2, flow_O3, flow_BCARY
      REAL :: V_chamber, conv1, conv2
      INTEGER :: i

      ! flows in ml/min, converted to m3/s:
      conv1 = 1E-6 / 60.
      flow_1 = 1200. * conv1 ! synthetic air (80 % N2 and 20 % O2)
      flow_2 = 1200. * conv1 ! synthetic air with 1 Âµmol/mol ozone
      flow_3 =  192. * conv1 ! N2 with 30. nmol/mol sesquiterpene
      flow_4 =  192. * conv1 ! N2 with 24.312 Âµmol/mol OH producer or scavenger
      flow_t = flow_1 + flow_2 + flow_3 + flow_4 ! total flow
      ! chamber volume [L], converted to [m3]:
      V_chamber = 100. * 1E-3

      ! all species are removed from the chamber (only H2O is constant):
      DO i = 1,NSPEC
        IF (TRIM(SPC_NAMES(i))/='H2O') THEN
          c(i) = c(i) * (1. - timesteplen * flow_t / V_chamber)
        ENDIF
      ENDDO

      flow_O2    = 0.2      * (flow_1 + flow_2)
      flow_O3    = 2.32E-6  * flow_2
      flow_BCARY = 739.5E-9 * flow_3
      ! the rest is nitrogen:
      flow_N2    = flow_t - (flow_O2+flow_O3+flow_BCARY)

      ! conv2 converts a flux in [m3/s] to a concentration increase in
      ! the chamber per timestep in [molec/cm3]:
      conv2 = timesteplen * N_A * press / ( R_gas * temp * (1E6*V_chamber) )

      IF (ind_O2    /= 0) c(ind_O2)    = c(ind_O2)    + flow_O2    * conv2
      IF (ind_N2    /= 0) c(ind_N2)    = c(ind_N2)    + flow_N2    * conv2
      IF (ind_O3    /= 0) c(ind_O3)    = c(ind_O3)    + flow_O3    * conv2
      IF (ind_BCARY /= 0) c(ind_BCARY) = c(ind_BCARY) + flow_BCARY * conv2

    END SUBROUTINE emission_lab_c15

    !-------------------------------------------------------------------------

    SUBROUTINE emission_mom
      IF (ind_NO     /= 0) c(ind_NO)     = c(ind_NO)     + 3.33E9 * fct
    END SUBROUTINE emission_mom

    !-------------------------------------------------------------------------

    SUBROUTINE emission_iso
#ifdef MECCA_TAG

#include "messy_mecca_tag_parameters.inc"
#include "./mecca/tag/caaba_exp.inc"

      use caaba_mem,           only: cair, temp
      use messy_mecca_tag_box, only: mecca_tag_depos

#ifdef ISO_STRAT
      real(dp) :: iwish, cfct      ! desired MR, corr. factor

      iwish = 1e-9_dp * cair  ! 1 ppb of H2O
      if ( C(ind_H2O) .gt. 0. ) then
        cfct = ( iwish / C(ind_H2O) )          ! corr. factor
        C(ind_H2O) = C(ind_H2O) * cfct         ! correcting regular
        call mecca_tag_depos(ind_H2O, cfct)    !    and all tagged
      endif

      iwish = 551e-9_dp * cair  ! 551 ppb of H2
      if ( C(ind_H2) .gt. 0. ) then
        cfct = ( iwish / C(ind_H2) )           ! corr. factor
        C(ind_H2) = C(ind_H2) * cfct           ! correcting regular
        call mecca_tag_depos(ind_H2, cfct)     !    and all tagged
      endif
#endif

#ifdef xxxEXP_STOZ
      C(ind_FSO3) = c(ind_FSO3) + 1e-3 * C(ind_O3)
      C(ind_FTO3) = c(ind_FTO3) + 1e-2 * C(ind_O3)
#endif


#ifdef CIO_YEUNG_EXP
    ! we simulate Yeung et al. [JGR 2014] exp. lab conditions
    ! decrease T every 5 mins by 20 degrees

      if ( (mod(model_time-model_start_day*OneDay,5.*60.) <= 1.) &
                .and. (model_time-model_start_day*OneDay > 60. ) ) temp = temp - 20.
#endif

#ifdef tag_FN
    ! C(ind_O3)     = C(ind_O3)     + 5.E10 * fct ! ref0203, p. 6699
    ! C(ind_NH3)    = C(ind_NH3)    + 1.0E9 * fct

    ! emit NO around noon as the road source (and track this emitted NO in NO2 and other species)
!     if ( abs(model_time-(model_start_day+1./24.)*OneDay) <= 5. ) then
      if ( abs(model_time-(130.+15./24.)*OneDay) <= 5. ) then
        C(ind_NO)     = C(ind_NO)     + 1.0E-12 * cair * timesteplen
        C(ind_FNNO)   = C(ind_FNNO)   + 1.0E-12 * cair * timesteplen
        print *,"1 ppt/s plum of NO over the road"
      endif
#endif // tag_FN

#ifdef tag_COSK
!      print *,C(ind_CO),C(ind_E13CO),C(ind_E17CO),C(ind_E18CO)
      C(ind_CO)     = C(ind_CO)     + 1e-11 * cair
      C(ind_E17CO)  = C(ind_E17CO)  + 1e-11 * cair
      C(ind_E13CO)  = C(ind_E13CO)  + 1e-11 * cair
      C(ind_E18CO)  = C(ind_E18CO)  + 1e-11 * cair
      C(ind_E17pCO) = C(ind_E17pCO) + 1e-11 * cair
      C(ind_E13pCO) = C(ind_E13pCO) + 1e-11 * cair
      C(ind_E18pCO) = C(ind_E18pCO) + 1e-11 * cair
#endif

#endif
    END SUBROUTINE emission_iso

    !-------------------------------------------------------------------------

    SUBROUTINE emission_tag
#ifdef MECCA_TAG
      ! emissions for tagging are to be treated here
#endif
    END SUBROUTINE emission_tag

    !-------------------------------------------------------------------------

    SUBROUTINE emission_ff
      IF (ind_NO     /= 0) c(ind_NO)     = c(ind_NO)     + 1.0E8 * fct
    END SUBROUTINE emission_ff

    !-------------------------------------------------------------------------

    SUBROUTINE emission_oomph
      ! OOMPH MARINE BACKGROUND
      IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3) = c(ind_CH3COCH3) + 1.4E7 * fct ! 4Igen
      !IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3) = c(ind_CH3COCH3) + 1.1E7 * fct ! 4GEN
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)   = c(ind_CH3CHO)   + 1.55E9 * fct ! 4Igen
      !IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)   = c(ind_CH3CHO)   + 1.3E9 * fct ! 4GEN
      IF (ind_DMS    /= 0) c(ind_DMS)    = c(ind_DMS)    + 2.9E9 * fct ! 4Igen
      !IF (ind_DMS    /= 0) c(ind_DMS)    = c(ind_DMS)    + 2.55E9 * fct ! 4GEN
      !IF (ind_DMS    /= 0) c(ind_DMS)    = c(ind_DMS)    + 3.5E9 * fct ! 3rd gen
      !IF (ind_DMS    /= 0) c(ind_DMS)    = c(ind_DMS)    + 2.0E9 * fct
      IF (ind_C5H8   /= 0) c(ind_C5H8)   = c(ind_C5H8)   + 5._dp * N_A/(5._dp*MC+8._dp*MH) * fct3 * fct ! 4Igen
      !IF (ind_C5H8   /= 0) c(ind_C5H8)   = c(ind_C5H8)   + 4.25_dp * N_A/(5_dp*MC+8_dp*MH) * fct3 * fct ! 4gen
      IF (ind_NH3    /= 0) c(ind_NH3)    = c(ind_NH3)    + 1.0E9 * fct
      IF (ind_NO     /= 0) c(ind_NO)     = c(ind_NO)     + 0.1E9 * fct
      IF (ind_O3     /= 0) c(ind_O3)     = c(ind_O3)     + 5.E10 * fct
      IF (ind_CHBr3  /= 0) c(ind_CHBr3)  = c(ind_CHBr3)  + 221.4 * fct2 ! Lucy
      IF (ind_CH3I   /= 0) c(ind_CH3I)   = c(ind_CH3I)   + 6.0E6 * fct ! ref0897
      IF (ind_C3H7I  /= 0) c(ind_C3H7I)  = c(ind_C3H7I)  + 1.0E7 * fct ! ref0897
      ! emissions of short-lived iodocarbons are scaled by a factor of 10 to
      ! indicate that they are not evenly distributed throughout the whole mbl
      IF (ind_CH2ClI /= 0) c(ind_CH2ClI) = c(ind_CH2ClI) + 185. *10.*fct2 !Lucy
      IF (ind_CH2I2  /= 0) c(ind_CH2I2)  = c(ind_CH2I2)  + 54.7 *10.*fct2 !Lucy
      !IF (ind_CH2ClI /= 0) c(ind_CH2ClI) = c(ind_CH2ClI) + 2.0E7 * fct ! ref0897
      !IF (ind_CH2I2  /= 0) c(ind_CH2I2)  = c(ind_CH2I2)  + 3.0E7 * fct ! ref0897
    END SUBROUTINE emission_oomph

    !-------------------------------------------------------------------------

    !mz_hr_20130425+
    SUBROUTINE emission_cumulus
    !mz_hr_20130514+ fixed O3, corrected mixing ratios and time
    !mz_hr_20140127+ rev: for paper revision
    ! emis only HCHO or only CH3OOH or NO + HCHO or NO + CH3OOH
    ! additional NO(x) EMISSIONS via injectNOx
      !IF (l_injectNOx .AND. &
      !  (model_time >= (t_NOxon+model_start_day)*OneDay)   .AND. &
      !  (model_time <= (t_NOxoff+model_start_day)*OneDay)) THEN
      !PRINT *, 'emissions'
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  82.6E9 * fct ! bg4 + only HCHO, 400 ppt
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  56.2E9 * fct ! bg4 + only HCHO, 300 ppt
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  108.7E9 * fct ! bg4 + only HCHO, 500 ppt
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 70.7E9 * fct ! only CH3OOH -> 400 ppt HCHO
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 49.1E9 * fct ! only CH3OOH -> 300 ppt HCHO
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 91E9 * fct ! only CH3OOH -> 500 ppt HCHO
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  40.5E9 * fct ! isop only -> 400 ppt HCHO
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  29.8E9 * fct ! isop only -> 300 ppt HCHO
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  50.9E9 * fct ! isop only -> 500 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 3.0E9 * fct ! NO only -> 1.3 ppb NO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 3.33E9 * fct ! NO only, init 46 ppt NO2 -> 0.8 ppb NO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 2.52E9 * fct ! NO only, init 60 ppt NO2 -> 0.8 ppb NO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 1.7E9 * fct ! NO only, init 100 ppt NO2 -> 0.8 ppb NO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 4.33E9 * fct ! NO only, init 170 ppt NO2 -> 1.7 ppb NO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 2.0E9 * fct ! NO only, eq of init2.0ppbNO init70pptNO2 -> 1.7 ppb NO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 1.3E9 * fct ! NO only, init 70 ppt NO2 -> 1.7 ppb NO

      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  84.6E9 * fct ! B4 + HCHO only, 0>400 ppt or 400>400ppt
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  83.6E9 * fct ! B4 + HCHO only, 400 ppt*0.99
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  85.6E9 * fct ! B4 + HCHO only, 400 ppt*1.01
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 70.1E9 * fct ! only CH3OOH -> 400 ppt HCHO
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 69.3E9 * fct ! only CH3OOH -> 400 ppt*0.99 HCHO
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 70.9E9 * fct ! only CH3OOH -> 400 ppt*1.01 HCHO
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  40.3E9 * fct ! isop only -> 400 ppt HCHO
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  39.9E9 * fct ! isop only -> 400ppt*0.99 HCHO
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  40.8E9 * fct ! isop only -> 400ppt*1.01 HCHO

      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 4.6E9 * fct ! B4 + NO only !B4> 1.3 ppb NO, SS with 5 specs (no HO2)
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 3.7E9 * fct ! B4 + NO only !B4> 1.3*0.99 ppb NO, SS 5 specs noHO2
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 5.3E9 * fct ! B4 + NO only !B4> 1.3*1.01 ppb NO, SS 5 specs noHO2
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 2.1E9 * fct ! B4 + NO only !B4> 0.8 ppb NO, SS 5 specs noHO2
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 1.6E9 * fct ! B4 + NO only !B4> 0.8*0.99 ppb NO, SS with 5 specs
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 2.8E9 * fct ! B4 + NO only !B4> 0.8*1.01 ppb NO, SS with 5 specs
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 4.8E9 * fct ! B4 + NO only !B4> 1.7 ppb NO, SS with 5 specs (no HO2)
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 4.6E9 * fct ! B4 + NO only !B4> 1.7*0.99 ppb NO, SS with 5 specs (no HO2), no *1.01 stable
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 3.7E9 * fct ! B4 + NO only !B4> 1.17 ppb NO, SS with 5 specs (no HO2)

      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  120.9E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 10.0E9 * fct   ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  118.4E9 * fct ! 1.3*0.99 ppb NO + 400*0.99 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +    8.9E9 * fct   !  !1.3*0.99 ppb NO + 400*0.99 ppt HCHO
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  123.6E9 * fct ! 1.3*1.01 ppb NO + 400*1.01 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +   11.3E9 * fct   !  !1.3*1.01 ppb NO + 400*1.01 ppt HCHO
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) +  91.9E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_NO /= 0)     c(ind_NO)    = c(ind_NO)      +  8.5E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) +  90.3E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_NO /= 0)     c(ind_NO)    = c(ind_NO)      +  8.1E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) +  94.2E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_NO /= 0)     c(ind_NO)    = c(ind_NO)      +  8.8E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  35.5E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_NO /= 0)     c(ind_NO) = c(ind_NO)   +   7.4E9 * fct ! 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  34.9E9 * fct ! 0.99 1.3 ppb NO + 400 ppt HCHO
      !IF (ind_NO /= 0)     c(ind_NO) = c(ind_NO)   +   7.3E9 * fct ! 0.99 1.3 ppb NO + 400 ppt HCHO

      !!! combinations of NO and VOC
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  111.0E9 * fct ! 0.8 ppb NO + 400 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 6.0E9 * fct
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  66.6E9 * fct ! 0.8 ppb NO + 300 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 3.7E9 * fct
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  160.5E9 * fct ! 0.8 ppb NO + 500 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 7.5E9 * fct
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  177.6E9 * fct ! 1.3 ppb NO + 500 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 11.6E9 * fct
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  70.0E9 * fct ! 1.3 ppb NO + 300 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 8.0E9 * fct
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  123.4E9 * fct ! 1.65 ppb NO + 400 ppt HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 9.5E9 * fct
      IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  71.0E9 * fct ! 1.7 ppb NO + 500 ppt HCHO
      IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 8.0E9 * fct
      !ENDIF

      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) + 117.8E9 * fct ! HCHO + NO
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 91E9 * fct ! CH3OOH + NO
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) + 102E9 * fct ! only HCHO !old J_NO2
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) + 108E9 * fct ! only HCHO into coudmed or bg2CAR
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  87E9 * fct ! only HCHO into bg2
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  83E9 * fct ! only HCHO into bg2CARlowN
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  72E9 * fct ! only HCHO into bg2CARlowN -10%
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  93E9 * fct ! only HCHO into bg2CARlowN +10%
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  80E9 * fct ! only HCHO into bg3CARlowN
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  70E9 * fct ! only HCHO into bg3CARlowN -10%
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  90E9 * fct ! only HCHO into bg3CARlowN +10%
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  121E9 * fct ! NO+HCHO, init bg3CARhiN_hiVOC
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  6.5E9 * fct ! NO+HCHO, bg3CARhiN_hiVOC, 265K
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  0.0E9 * fct ! NO+HCHO, bg3CARhiN_hiVOC, 272K
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  120E9 * fct ! NO+HCHO, init bg3CARhiN_hiVOC ave
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) +  50E9 * fct ! rev NO+HCHO, init bg3CARlowN
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 80E9 * fct ! only !CH3OOH
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 121E9 * fct ! same !as HCHO (HCHO+NO)
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 70E9 * fct ! only CH3OOH, bg2CARlowN init
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 62E9 * fct ! only CH3OOH, bg2CARlowN init -10%
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 79E9 * fct ! only CH3OOH, bg2CARlowN init +10%
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 68E9 * fct ! only CH3OOH, bg3CARlowN init
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 60E9 * fct ! only CH3OOH, bg3CARlowN init -10%
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) + 77E9 * fct ! only CH3OOH, bg3CARlowN init +10%
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) +  92E9 * fct ! NO+CH3OOH, init bg3CARhiN_hiVOC
      !IF (ind_CH3OOH /= 0) c(ind_CH3OOH) = c(ind_CH3OOH) +  103E9 * fct ! NO+CH3OOH, init bg3CARhiN_hiVOC ave
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  40.1E9 * fct ! isop only, init bg3CARhiN_hiVOC
      !IF (ind_C5H8 /= 0) c(ind_C5H8) = c(ind_C5H8) +  36.1E9 * fct ! NO+isop, init bg3CARhiN_hiVOC
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 2.6E9 * fct ! bg3CARhiNloVOC, only NO emis
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 3.2E9 * fct ! bg3CARhiNhiVOC, only NO emis
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 2.6E9 * fct ! bg3CARhiNhiVOC, only NO emis, ave = 1.316ppb
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 1.35E9 * fct ! bg3CARhiNhiVOC, only NO emis -10%
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 5.7E9 * fct ! bg3CARhiNhiVOC, only NO emis +10%
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 54.0E9 * fct ! bg3CARhiNhiVOC, only NO emis maxNO, min so that 1st step NO NOT down
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 1.0E9 * fct ! bg3CARhiNloVOC, only NO emis -10%
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 6.56E9 * fct ! bg3CARhiNloVOC, only NO emis +10%
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 12.7E9 * fct ! bg3CARhiNloVOC, only NO emis, maxNO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +  8.5E9 * fct ! bg3CARhiN_hiVOC, NO+HCHO
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +  8.05E9 * fct ! bg3CARhiN_hiVOC, NO+isop
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +  11.1E9 * fct ! bg3CARhiN_hiVOC, NO+HCHO, 265K
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +  11.9E9 * fct ! bg3CARhiN_hiVOC, NO+HCHO, 272K
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +  7.55E9 * fct ! bg3CARhiN_hiVOC, NO+HCHO ave
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +  8.2E9 * fct ! bg3CARhiN_hiVOC, NO+CH3OOH
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   +  7.55E9 * fct ! bg3CARhiN_hiVOC, NO+CH3OOH ave
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 12.9E9 * fct ! HCHO fixed
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 12.9E9 * fct ! CH3OOH fixed
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 14.4E9 * fct ! (HCHO or CH3OOH) + NO thresh 1e-6 rises later
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 14.4E9 * fct * (1 - 7.6E-2 - 6.7E-3) ! Bhetanabhotla lightning ratios
      !IF (ind_NO2 /= 0) c(ind_NO2)   = c(ind_NO2)  + 14.4E9 * fct * 7.6E-2 ! Bhetanabhotla ratio for flashes
      !IF (ind_HONO /= 0) c(ind_HONO) = c(ind_HONO) + 14.4E9 * fct * 6.7E-3 ! Bhetanabhotla ratio for flashes
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 6.7E9 * fct ! only NO thresh 1e-6, rises later
      !IF (ind_NO /= 0) c(ind_NO)     = c(ind_NO)   + 6.3E9 * fct ! new J_NO2 = 0.025 only NO thresh 1e-6, rises later
      !IF (ind_NO2 /= 0) c(ind_NO2)   = c(ind_NO2)  + 6.7E9 * fct * 7.6E-2 ! Bhetanabhotla ratio for flashes
      !IF (ind_HONO /= 0) c(ind_HONO) = c(ind_HONO) + 6.7E9 * fct * 6.7E-3 ! Bhetanabhotla ratio for flashes
      !IF (ind_NO /= 0) c(ind_NO) = c(ind_NO) + 0.6E9 * fct ! HCHO + NO thresh 1e-7 too low rises later
    ! emissions with unfixed O3
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) + 9.E10 * fct ! only HCHO or with NO
      !IF (ind_NO /= 0) c(ind_NO) = c(ind_NO) + 10.E9 * fct ! if with 16.E10 HCHO
    ! use these emissions if O3 in set of fixed species:
      !IF (ind_HCHO /= 0) c(ind_HCHO) = c(ind_HCHO) + 16.E10 * fct ! only HCHO or with NO
      !IF (ind_NO /= 0) c(ind_NO) = c(ind_NO) + 15.E9 * fct ! if with 16.E10 HCHO
      !!IF (ind_NO /= 0) c(ind_NO) = c(ind_NO) + 4.5E9 * fct ! if only NO emis
    END SUBROUTINE emission_cumulus
    !mz_hr_20130425-

    !-------------------------------------------------------------------------

    SUBROUTINE injection

      USE caaba_mem, ONLY: model_start_day, cair
      LOGICAL, SAVE :: l_not_yet_injected = .TRUE.
      REAL(dp) :: c_injection, t_injection
      INTEGER :: ind_species

      ! define species, injection time and mixing ratio:
      ind_species = ind_NO2
      t_injection = model_start_day*OneDay + 43200. ! 12 h after model start
      c_injection = 1E-9 * cair ! 1 nmol/mol

      ! injection happens at the specified time (or slightly later,
      ! depending on timesteplen)
      IF (model_time >= t_injection .AND. l_not_yet_injected) THEN
        c(ind_species) = c_injection
        l_not_yet_injected = .FALSE.
      ENDIF

    END SUBROUTINE injection

    !-------------------------------------------------------------------------

  END SUBROUTINE emission

  !***************************************************************************

  SUBROUTINE drydep

    ! deposition velocities [cm/s]

    USE caaba_mem, ONLY: drydep_scenario
    IMPLICIT NONE
    REAL(dp) :: fct

    fct = timesteplen / (zmbl * 100.)

    SELECT CASE (TRIM(drydep_scenario))
    CASE ('FF_ARCTIC','FF_ANTARCTIC')
      CALL drydep_ff
    CASE ('OOMPH','MBL')
      CALL drydep_mbl
    CASE ('LAB','LAB_C15')
      CALL drydep_lab
    CASE ('MOM')
      CALL drydep_mom
    CASE ('VOLCANO')
      CALL dilute
      ! CALL dilute_once
    CASE ('')
      CALL drydep_default
    CASE DEFAULT
      PRINT *, 'ERROR, drydep_scenario '//TRIM(drydep_scenario)// &
        ' is not defined'
      STOP 1
    END SELECT

    !-------------------------------------------------------------------------

  CONTAINS

    !-------------------------------------------------------------------------

    SUBROUTINE drydep_ff
      ! no deposition for frostflower model setup
    END SUBROUTINE drydep_ff

    !-------------------------------------------------------------------------

    SUBROUTINE drydep_mbl
      ! default values for mbl:
      IF (ind_O3       /= 0) c(ind_O3)       = (1.-fct*0.04) * c(ind_O3)
      IF (ind_H2O2     /= 0) c(ind_H2O2)     = (1.-fct*0.5)  * c(ind_H2O2)
      IF (ind_NH3      /= 0) c(ind_NH3)      = (1.-fct*0.1)  * c(ind_NH3)
      IF (ind_NO2      /= 0) c(ind_NO2)      = (1.-fct*0.1)  * c(ind_NO2)
      IF (ind_N2O5     /= 0) c(ind_N2O5)     = (1.-fct*1.0)  * c(ind_N2O5)
      IF (ind_HNO3     /= 0) c(ind_HNO3)     = (1.-fct*2.0)  * c(ind_HNO3)
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)    = (1.-fct*0.08) * c(ind_CH3OH) ! 4GEN, Jacob2005
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)   = (1.-fct*0.1)  * c(ind_CH3OOH)
      IF (ind_HCHO     /= 0) c(ind_HCHO)     = (1.-fct*0.5)  * c(ind_HCHO)
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)    = (1.-fct*1.0)  * c(ind_HCOOH)
      !IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3) = (1.-fct*0.1)  * c(ind_CH3COCH3) ! 4GEN, Jacob2005
      IF (ind_HCl      /= 0) c(ind_HCl)      = (1.-fct*2.0)  * c(ind_HCl)
      IF (ind_HOCl     /= 0) c(ind_HOCl)     = (1.-fct*1.0)  * c(ind_HOCl)
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)    = (1.-fct*1.0)  * c(ind_ClNO3)
      IF (ind_HBr      /= 0) c(ind_HBr)      = (1.-fct*2.0)  * c(ind_HBr)
      IF (ind_HOBr     /= 0) c(ind_HOBr)     = (1.-fct*1.0)  * c(ind_HOBr)
      IF (ind_BrNO3    /= 0) c(ind_BrNO3)    = (1.-fct*1.0)  * c(ind_BrNO3)
      IF (ind_I2O2     /= 0) c(ind_I2O2)     = (1.-fct*1.0)  * c(ind_I2O2)
      IF (ind_HI       /= 0) c(ind_HI)       = (1.-fct*1.0)  * c(ind_HI)
      IF (ind_HOI      /= 0) c(ind_HOI)      = (1.-fct*1.0)  * c(ind_HOI)
      IF (ind_INO2     /= 0) c(ind_INO2)     = (1.-fct*1.0)  * c(ind_INO2)
      IF (ind_INO3     /= 0) c(ind_INO3)     = (1.-fct*1.0)  * c(ind_INO3)
      IF (ind_SO2      /= 0) c(ind_SO2)      = (1.-fct*0.5)  * c(ind_SO2)
      IF (ind_H2SO4    /= 0) c(ind_H2SO4)    = (1.-fct*2.0)  * c(ind_H2SO4)
      IF (ind_CH3SO3H  /= 0) c(ind_CH3SO3H)  = (1.-fct*1.0)  * c(ind_CH3SO3H)
      IF (ind_DMSO     /= 0) c(ind_DMSO)     = (1.-fct*1.0)  * c(ind_DMSO)
    END SUBROUTINE drydep_mbl

    !-------------------------------------------------------------------------

    SUBROUTINE drydep_lab
      IF (ind_O3       /= 0) c(ind_O3)       = (1.-fct*0.04) * c(ind_O3)
    END SUBROUTINE drydep_lab

    !-------------------------------------------------------------------------

    SUBROUTINE drydep_mom
      IF (ind_O3       /= 0) c(ind_O3)       = (1.-fct*0.5)  * c(ind_O3)
      IF (ind_H2O2     /= 0) c(ind_H2O2)     = (1.-fct*0.5)  * c(ind_H2O2)
      IF (ind_NH3      /= 0) c(ind_NH3)      = (1.-fct*0.1)  * c(ind_NH3)
      IF (ind_NO2      /= 0) c(ind_NO2)      = (1.-fct*0.1)  * c(ind_NO2)
      IF (ind_N2O5     /= 0) c(ind_N2O5)     = (1.-fct*1.0)  * c(ind_N2O5)
      IF (ind_HNO3     /= 0) c(ind_HNO3)     = (1.-fct*2.0)  * c(ind_HNO3)
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)    = (1.-fct*0.08) * c(ind_CH3OH)
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)   = (1.-fct*0.1)  * c(ind_CH3OOH)
      IF (ind_HCHO     /= 0) c(ind_HCHO)     = (1.-fct*0.5)  * c(ind_HCHO)
    END SUBROUTINE drydep_mom

    !-------------------------------------------------------------------------

    SUBROUTINE drydep_default
      ! no deposition for simple default setup
    END SUBROUTINE drydep_default

    !-------------------------------------------------------------------------

    SUBROUTINE dilute

      ! dilute air in a plume by adding fresh air to it

      ! set l_ignore_relhum=T when dilution should also affect c(ind_H2O)
      ! TODO to use dilute with traject, note that at this point trajectory
      ! data and cair has not been updated yet.

      USE caaba_mem,                ONLY: cair, timesteplen
      USE messy_main_constants_mem, ONLY: OneDay
      USE messy_mecca_kpp           ! ind_*, NSPEC
      IMPLICIT NONE
      REAL(DP) :: factor
      REAL(DP), PARAMETER :: &
        exchng = 1. / (3.*OneDay) ! exchange with fresh air [1/s]
      INTEGER :: i

      ! add fresh air:
      factor = timesteplen * exchng
      ! adjustment for operator splitting:
      factor = factor / (1.-factor)
      !PRINT *, 'dilution factor = ', factor
      IF (ind_O3  /= 0) c(ind_O3)  = c(ind_O3)  + factor *  80.E-09 * cair
      IF (ind_O2  /= 0) c(ind_O2)  = c(ind_O2)  + factor * 210.E-03 * cair
      IF (ind_N2  /= 0) c(ind_N2)  = c(ind_N2)  + factor * 780.E-03 * cair
      IF (ind_NO2 /= 0) c(ind_NO2) = c(ind_NO2) + factor *   5.E-09 * cair
      IF (ind_CH4 /= 0) c(ind_CH4) = c(ind_CH4) + factor *  1.8E-06 * cair
      IF (ind_CO  /= 0) c(ind_CO)  = c(ind_CO)  + factor *  70.E-09 * cair
      IF (ind_CO2 /= 0) c(ind_CO2) = c(ind_CO2) + factor * 350.E-06 * cair

      ! remove plume air:
      factor = 1. - timesteplen * exchng
      DO i = 1,NSPEC
        ! all species are lost (except for reaction-rate pseudo species RR*)
        IF (INDEX(SPC_NAMES(i),'RR') /= 1) THEN
          !PRINT *,'dilution loss of ', TRIM(SPC_NAMES(i))
          c(i) = c(i) * factor
        ELSE
          !PRINT *,'no dilution loss of ', TRIM(SPC_NAMES(i))
        ENDIF
      ENDDO

    END SUBROUTINE dilute

    !-------------------------------------------------------------------------

    SUBROUTINE dilute_once

      ! dilute air in a plume by adding fresh air to it
      USE caaba_mem,                ONLY: cair, model_time, model_start_day
      USE messy_main_constants_mem, ONLY: OneDay
      USE messy_mecca_kpp           ! ind_*, NSPEC
      IMPLICIT NONE
      LOGICAL, SAVE :: l_already_diluted = .FALSE.
      REAL(DP), PARAMETER :: t_dilution = 30. ! [s] enter value here
      REAL(DP), PARAMETER :: factor = 1000. ! enter value here
      REAL(DP), PARAMETER :: factor2 = (factor-1.)/factor
      INTEGER :: i

      IF (.NOT.l_already_diluted) THEN
        IF (model_time >= t_dilution+model_start_day*OneDay) THEN
          DO i = 1,NSPEC
            ! all species are diluted
            ! (except for reaction-rate pseudo species RR*)
            IF (INDEX(SPC_NAMES(i),'RR') /= 1) THEN
              PRINT *,'dilution loss of ', TRIM(SPC_NAMES(i))
              c(i) = c(i) / factor
            ELSE
              PRINT *,'no dilution loss of ', TRIM(SPC_NAMES(i))
            ENDIF
          ENDDO
          ! add fresh air:
          IF (ind_O2 /= 0) c(ind_O2) = c(ind_O2) + factor2 * 210.E-03 * cair
          IF (ind_N2 /= 0) c(ind_N2) = c(ind_N2) + factor2 * 780.E-03 * cair
          l_already_diluted = .TRUE.
        ENDIF
      ENDIF

    END SUBROUTINE dilute_once

    !-------------------------------------------------------------------------

  END SUBROUTINE drydep

  !***************************************************************************

END MODULE messy_semidep_box

!*****************************************************************************
