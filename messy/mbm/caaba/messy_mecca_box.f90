!*****************************************************************************
! Time-stamp: <2020-07-10 16:12:12 b302010>
!*****************************************************************************

! MECCA: Module Efficiently Calculating the Chemistry of the Atmosphere

! Authors:
! Rolf Sander,    MPICH, Mainz, 2003-...: original code
! Astrid Kerkweg, MPICH, Mainz, 2003-2007: halogen/aerosol chemistry
! Hella Riede,    MPICH, Mainz, 2007:

!*****************************************************************************

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

MODULE messy_mecca_box

  USE caaba_mem,                ONLY: model_start_day, model_time,           &
                                      model_end, timesteplen, init_spec,     &
#ifdef __GFORTRAN__
                                      C_ => &
#endif
                                      C,                                     &
                                      istatus, ierrf,                        &
                                      cair, temp, press, relhum, spechum,    &
                                      cair_old, l_ignore_relhum, l_hum_emac, &
                                      l_relhum_wmo, l_psat_liquid,           &
                                      l_ff, Ca_precip, l_steady_state_stop,  &
                                      l_RRconc, l_skeleton, photrat_channel
#ifdef __GFORTRAN__    /* newer gfortran versions have a run-time bug of    */
#define C C_           /* module-level import of C(:), so it has to be      */
#define c C_           /* renamed into something else, e.g. C_(:)           */
#endif
  USE messy_cmn_photol_mem      ! IP_MAX, ip_*, jname
  USE messy_mecca_kpp           ! DP, ind_*, kpp_integrate, rconst,
                                ! APN, initialize_indexarrays,
                                ! SPC_NAMES, EQN_TAGS, NSPEC, NREACT, batchfile
#ifdef MECCA_TAG
  USE messy_mecca_tag_box,      ONLY: mecca_tag_init, mecca_tag_preprocess, &
                                      mecca_tag_postprocess, &
                                      mecca_tag_result, mecca_tag_finish
#endif
  USE messy_main_constants_mem, ONLY: N_A, rho_H2O, M_H2O, R_gas,            &
                                      OneDay, STRLEN_VLONG,                  &
                                      pi, HLINE2
  USE messy_mecca,              ONLY: l_aero, l_tag, modstr

  USE caaba_io,                 ONLY: nf, open_input_file, nf90_inquire,     &
                                      nf90_inquire_variable, nf90_inq_varid, &
                                      nf90_get_var, open_output_file,        &
                                      write_output_file, close_file

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: NBL = 1 ! N_block_length
  ! declarations for all variables that are transported to the SMCL via fill
  ! subroutines (or via kpp_integrate for "C")
  REAL(DP) :: jx(IP_MAX) = 0.
  REAL(DP), DIMENSION(:), ALLOCATABLE :: xaer
  REAL(DP) :: cvfac(APN)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: lwc
  REAL(DP) :: k_exf(APN,NSPEC)
  REAL(DP) :: k_exb(APN,NSPEC)
  REAL(DP) :: k_exf_N2O5(APN)
  REAL(DP) :: k_exf_ClNO3(APN)
  REAL(DP) :: k_exf_BrNO3(APN)
  REAL(DP) :: dummy_khet_Tr(IHT_MAX) ! dummy, not needed in box model
  REAL(DP) :: dummy_khet_St(IHS_MAX) ! dummy, not needed in box model
  ! (for C, cair, and temp, see caaba_mem.f90)
  ! (see also important notes about temp, press, and cair in gas.eqn!)

  REAL(DP), DIMENSION(NSPEC) :: output
  CHARACTER(LEN=20) :: c_unit(NSPEC) ! c unit
  CHARACTER(LEN=20), PARAMETER :: ratesa_unit(NREACT) = 'cm-3s-1'
  INTEGER :: ncid_aero, ncid_tracer, ncid_spec
  INTEGER :: ncid_ratesa

  ! sea water composition, relative to total salt:
  REAL(DP) :: HCO3m_rel, NO3m_rel, Clm_rel, Brm_rel, Im_rel, IO3m_rel, &
    SO4mm_rel

  REAL(DP), DIMENSION(:), ALLOCATABLE :: csalt, exchng, radius
  REAL(DP), DIMENSION(:), ALLOCATABLE :: c0_NH4p, c0_Nap, c0_FeOHp, c0_FeCl2p
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    c0_HCO3m, c0_NO3m, c0_Clm, c0_Brm, c0_Im, c0_IO3m, c0_SO4mm, c0_HSO4m

  REAL(DP) :: mcexp(MAX_MCEXP) ! Monte-Carlo factor

  PRIVATE
  PUBLIC :: mecca_init   ! initialize aero chemistry
  PUBLIC :: mecca_physc  ! calculate chemistry
  PUBLIC :: mecca_result ! print results
  PUBLIC :: mecca_finish ! close files

CONTAINS

  !***************************************************************************

  SUBROUTINE mecca_init

    !  mz_hr_20160509 convert from spechum => all issues with psat
    !  parameterization and relhum definition dealt with in rel2spechum (and
    !  spec2relhum), relhum <-> spechum handled in caaba_read_nml (caaba.f90),
    !  and for TRAJECT mode in get_model_start_end and get_physc_data
    !  (messy_traject_box.f90)
    USE messy_main_tools,         ONLY: str
    USE messy_mecca_kpp,          ONLY: REQ_MCFCT, initialize
    USE messy_mecca,              ONLY: mecca_read_nml_ctrl
    USE messy_main_constants_mem, ONLY: MH, MC, MO, MS, MCl, MBr

    IMPLICIT NONE

    INTRINSIC :: TRIM

    CHARACTER(LEN=20)   :: lwc_unit(APN)
    CHARACTER(LEN=7)    :: lwc_name(APN)
    INTEGER, PARAMETER  :: iou = 999   ! I/O unit
    INTEGER             :: status ! error status
    INTEGER             :: i, jb

    REAL(DP), PARAMETER :: rho_sw = 1025._dp    ! density of sea water [kg/m3]
    REAL(DP) :: csw_HCO3m, csw_NO3m, csw_Clm, csw_Brm, csw_Im, csw_IO3m, &
                csw_SO4mm, csw_total

    ALLOCATE(xaer(APN))
    ALLOCATE(radius(APN))
    ALLOCATE(lwc(APN))
    ALLOCATE(csalt(APN))
    ALLOCATE(c0_NH4p(APN))
    ALLOCATE(c0_Nap(APN))
    ALLOCATE(c0_FeOHp(APN))
    ALLOCATE(c0_FeCl2p(APN))
    ALLOCATE(c0_HCO3m(APN))
    ALLOCATE(c0_NO3m(APN))
    ALLOCATE(c0_Clm(APN))
    ALLOCATE(c0_Brm(APN))
    ALLOCATE(c0_Im(APN))
    ALLOCATE(c0_IO3m(APN))
    ALLOCATE(c0_SO4mm(APN))
    ALLOCATE(c0_HSO4m(APN))
    ALLOCATE(exchng(APN))

    ! read mecca ctrl namelist:
    CALL mecca_read_nml_ctrl(status, iou)
    IF (status /= 0) STOP 1

    ! read kpp ctrl namelist:
    CALL initialize_kpp_ctrl(status, iou, modstr)
    IF (status /= 0) STOP 1

    ! csw_* = concentration in sea water [mol/L]
    ! mass fraction w_i [kg/kg] (Tab. 4 of ref2307 = Millero et al. 2008)
    ! HCO3- + 2 CO3-- = 0.10481 + 2*0.01434 = 0.13349 g/kg
    csw_HCO3m = 0.00013349 * rho_sw / (MH+MC+MO*3.) * (1.-Ca_precip)
    csw_NO3m  = 1E-6 ! see Fig. 1D in ref2428 = Duce et al. 2008
    ! activate next line for Southern Ocean:
    !csw_NO3m  = 25E-6 ! ref2428 Southern Ocean
    csw_Clm   = 0.01935271 * rho_sw / MCl
    csw_Brm   = 0.00006728 * rho_sw / MBr
    csw_Im    = 7.4E-8  ! Tsunogai, Tab. 1a, ref0845
    csw_IO3m  = 2.64E-7 ! Tsunogai, Tab. 1a, ref0845
    csw_SO4mm = 0.00271235 * rho_sw / (MS+MO*4.)
    csw_total = csw_HCO3m + csw_NO3m + csw_Clm + csw_Brm + csw_Im &
      + csw_IO3m + csw_SO4mm

    ! sea water composition, relative to total salt
    HCO3m_rel = csw_HCO3m / csw_total
    NO3m_rel  = csw_NO3m  / csw_total
    Clm_rel   = csw_Clm   / csw_total
    Brm_rel   = csw_Brm   / csw_total
    Im_rel    = csw_Im    / csw_total
    IO3m_rel  = csw_IO3m  / csw_total
    SO4mm_rel = csw_SO4mm / csw_total

    IF (l_aero) THEN
      PRINT *, 'Sea water composition, relative to total salt:'
      PRINT *, 'HCO3m_rel = ', HCO3m_rel
      PRINT *, 'NO3m_rel  = ', NO3m_rel
      PRINT *, 'Clm_rel   = ', Clm_rel
      PRINT *, 'Brm_rel   = ', Brm_rel
      PRINT *, 'Im_rel    = ', Im_rel
      PRINT *, 'IO3m_rel  = ', IO3m_rel
      PRINT *, 'SO4mm_rel = ', SO4mm_rel
      PRINT *, 'SUM       = ', &
        HCO3m_rel+NO3m_rel+Clm_rel+Brm_rel+Im_rel+IO3m_rel+SO4mm_rel
    ENDIF

    C(:) = 0. ! default value unless explicitly initialized
    CALL initialize

    ! Monte-Carlo:
    IF (REQ_MCFCT) THEN
      CALL define_mcexp(status, mcexp)
      IF (status /= 0) STOP 1
      DO i=1, MAX_MCEXP
        WRITE(*,'(A,I4,A,F12.8)') ' mcexp(', i, ') = ', mcexp(i)
        ! mz_rs_20100727+
        ! In this example, only mcexp(40) is used, all other mcexp(:)
        ! are set to zero. Look at mecca.eqn to find which reaction is
        ! affect by mcexp(40). This can be used to make Monte-Carlo
        ! simulations where only one (or a few) rate coefficients are
        ! modified.
        ! IF (i/=40) mcexp(i) = 0.                 ! modify 1 rate coefficient
        ! IF ((i/=40).OR.(i/=50)) mcexp(i) = 0.    ! modify 2 rate coefficients
        ! mz_rs_20100727-
      ENDDO
      PRINT *, 'mcexp(avg)  = ', SUM(mcexp)/MAX_MCEXP
    ENDIF

    CALL x0 ! initial gas phase mixing ratios

    IF (l_aero) THEN
      ! define radius, lwc, etc.:
      CALL define_aerosol
      ! define gas-aq physicochemical constants:
      CALL mecca_aero_init_gasaq(l_print=.TRUE.)
    ENDIF

#ifdef MECCA_TAG
    IF (l_tag)  CALL mecca_tag_init
#endif
    ! ------------------------------------------------------------------------

    ! open output files and write headers:

    ! define units for reaction rates RR* and species:
    DO i = 1,NSPEC
      IF (INDEX(SPC_NAMES(i),'RR') == 1) THEN
        c_unit(i) = 'mol/mol/s'
      ELSE
        c_unit(i) = 'mol/mol'
      ENDIF
    ENDDO
    CALL open_output_file(ncid_tracer, 'caaba_mecca', SPC_NAMES, c_unit)

    CALL open_output_file(ncid_ratesa, 'caaba_ratesa', EQN_TAGS, ratesa_unit)

    IF (l_aero) THEN
      DO jb = 1, APN
        lwc_name(jb) = 'lwc_a'//str(jb,'(I2.2)')
      ENDDO
      lwc_unit(:) = 'm3/m3'
      CALL open_output_file(ncid_aero, 'caaba_mecca_aero', lwc_name, lwc_unit)
    ENDIF

    ! ------------------------------------------------------------------------

    PRINT *, HLINE2
    PRINT *, '*** Info about MECCA and KPP:'
    PRINT *, gas_spc_file,     ' (', gas_spc_file_sum,     ')'
    PRINT *, aqueous_spc_file, ' (', aqueous_spc_file_sum, ')'
    PRINT *, mecca_spc_file,   ' (', mecca_spc_file_sum,   ')'
    PRINT *, gas_eqn_file,     ' (', gas_eqn_file_sum,     ')'
    PRINT *, aqueous_eqn_file, ' (', aqueous_eqn_file_sum, ')'
    PRINT *, mecca_eqn_file,   ' (', mecca_eqn_file_sum,   ')'
    PRINT *, 'timestamp    = ', timestamp
    PRINT *, 'batchfile    = ', batchfile
    PRINT *, 'rplfile      = ', rplfile
    PRINT *, 'wanted       = ', wanted
    PRINT *, 'diagtracfile = ', diagtracfile
    PRINT *, 'tag          = ', tag
    PRINT *, 'kppoption    = ', kppoption
    PRINT *, 'KPP_HOME     = ', KPP_HOME
    PRINT *, 'KPP_version  = ', KPP_version
    PRINT *, 'integr       = ', integr

    ! ------------------------------------------------------------------------

  END SUBROUTINE mecca_init

  !***************************************************************************

  !***************************************************************************

    SUBROUTINE mecca_aero_init_gasaq(l_print)

    USE messy_main_tracer
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, HLINE2
    USE messy_mecca_aero,         ONLY: Henry_T0, Henry_Tdep, molar_mass &
                                      , alpha_T0, alpha_Tdep

    IMPLICIT NONE
    INTEGER :: icp ! Index ChemProp
    INTEGER :: jn
    LOGICAL, INTENT(IN), OPTIONAL :: l_print

    IF (PRESENT(l_print)) THEN
      IF (l_print) THEN
      PRINT *, HLINE2
      PRINT *, "         Henry's law coefficients "// &
        "and accommodation coefficients"
      PRINT *, HLINE2
      PRINT *, 'species           Henry_T0 Henry_Tdep'// &
        '   alpha_T0 alpha_Tdep         M'
      PRINT *, '                   [M/atm]        [K]'// &
        '        [1]        [K]   [g/mol]'
      PRINT *, HLINE2
      ENDIF
    ENDIF
    DO jn = 1,NSPEC
      icp = get_chemprop_index(TRIM(SPC_NAMES(jn)))
      IF (icp /=0) THEN
        Henry_T0(jn)   = chemprop(icp)%cask_r(R_Henry_T0)
        Henry_Tdep(jn) = chemprop(icp)%cask_r(R_Henry_Tdep)
        alpha_T0(jn)   = chemprop(icp)%cask_r(R_alpha_T0)
        alpha_Tdep(jn) = chemprop(icp)%cask_r(R_alpha_Tdep)
        molar_mass(jn) = chemprop(icp)%cask_r(R_MOLARMASS) / 1000. ! [kg/mol]
        IF (PRESENT(l_print)) THEN
          IF (l_print) THEN
            WRITE(*,'(1X,A15)', ADVANCE='NO') SPC_NAMES(jn)
            IF (Henry_T0(jn)>=0.) THEN
              WRITE(*,'(ES11.2)', ADVANCE='NO') Henry_T0(jn)
            ELSE
              WRITE(*,'(A)', ADVANCE='NO') '   --------'
            ENDIF
            IF (Henry_Tdep(jn)>=0.) THEN
              WRITE(*,'(F11.0)', ADVANCE='NO') Henry_Tdep(jn)
            ELSE
              WRITE(*,'(A)', ADVANCE='NO') '     ------'
            ENDIF
            IF (alpha_T0(jn)>=0.) THEN
              WRITE(*,'(ES11.2)', ADVANCE='NO') alpha_T0(jn)
            ELSE
              WRITE(*,'(A)', ADVANCE='NO') '   --------'
            ENDIF
            IF (alpha_Tdep(jn)>=0.) THEN
              WRITE(*,'(F11.0)', ADVANCE='NO') alpha_Tdep(jn)
            ELSE
              WRITE(*,'(A)', ADVANCE='NO') '     ------'
            ENDIF
            WRITE(*,'(F10.2)') 1000.*molar_mass(jn) ! [g/mol]
          ENDIF
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE mecca_aero_init_gasaq


  !***************************************************************************

  REAL(DP) FUNCTION r_wet_koehler(nions, T, rh)

    ! For an aerosol particle with a given amount of ions [mol] in it at
    ! a given temperature and relative humidity, calculate the wet
    ! radius according to the Koehler equation.
    USE messy_main_constants_mem, ONLY: MH, MO

    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: nions ! amount of substance of ions [mol]
    REAL(DP), INTENT(IN) :: T     ! temperature [K]
    REAL(DP), INTENT(IN) :: rh    ! relative humidity [1]
    REAL(DP) :: Akoehl, Bkoehl ! Koehler coefficients
    REAL(DP) :: a, b, c, d, p, q, DET, u3, v3, y ! for solving cubic equation
    REAL(DP) :: rhtest, nionstest, error ! error check
    REAL(DP), PARAMETER :: MH2O = 1E-3*(MH*2.+MO) ! [kg/mol]
    REAL(DP), PARAMETER :: sigma = 7.49e-2 ! surface tension [J/m2]

    ! Koehler equation
    Akoehl = 2. * sigma * MH2O / (rho_H2O * R_gas * T)
    Bkoehl = nions * MH2O / (4./3.*pi*rho_H2O)
    ! analytical solution of a cubic equation
    ! WARNING: does not work for rh approx 1
    ! (rh-1) * R^3 - Akoehl * R^2         + Bkoehl = 0
    ! a      * R^3 + b      * R^2 + c * R + d      = 0 ! general cubic eqn
    a     = (rh-1)
    b     = -Akoehl
    c     = 0.
    d     = Bkoehl
    p     = c/a - b**2/(3.*a**2)
    q     = 2./27. * (b/a)**3 - b*c/(3.*a**2) + d/a
    DET   = SQRT(q**2/4.+p**3/27.)
    u3    = -q/2. + DET
    v3    = -q/2. - DET
    y     = u3**(1./3.) + v3**(1./3.)
    r_wet_koehler = y - b/(3.*a) ! [m]
    ! check if the result is correct:
    rhtest = 1. + Akoehl/r_wet_koehler - Bkoehl/r_wet_koehler**3
    error = ABS(rh-rhtest)/rh
    IF (error>1E-2) PRINT *, rh, rhtest, error
    ! check if n_ions is correct:
    nionstest = nions_koehler(r_wet_koehler, T, rh)
    error = ABS(nions-nionstest)/nions
    IF (error>1E-2) PRINT *, nions, nionstest, error

  END FUNCTION r_wet_koehler

  !***************************************************************************

  REAL(DP) FUNCTION nions_koehler(r_wet, T, rh)

    ! For an aerosol particle with given wet radius at a given
    ! temperature and relative humidity, calculate the amount of ions
    ! [mol] in it according to the Koehler equation.
    USE messy_main_constants_mem, ONLY: MH, MO

    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: r_wet ! wet radius [m]
    REAL(DP), INTENT(IN) :: T     ! temperature [K]
    REAL(DP), INTENT(IN) :: rh    ! relative humidity [1]
    REAL(DP) :: Akoehl, Bkoehl ! Koehler coefficients
    REAL(DP), PARAMETER :: MH2O = 1E-3*(MH*2.+MO) ! [kg/mol]
    REAL(DP), PARAMETER :: sigma = 7.49e-2 ! surface tension [J/m2]

    ! Koehler equation: (rh-1) = Akoehl / R - Bkoehl / R^3
    ! <=> Bkoehl = Akoehl * R^2 - (rh-1) * R^3
    Akoehl = 2. * sigma * MH2O / (rho_H2O * R_gas * T)
    Bkoehl = Akoehl * r_wet**2 - (rh-1) * r_wet**3
    nions_koehler = Bkoehl * (4./3.*pi*rho_H2O) / MH2O

  END FUNCTION nions_koehler

  !***************************************************************************

  SUBROUTINE define_aerosol

    USE caaba_mem, ONLY: degree_lat, zmix, zmbl, init_scenario
    IMPLICIT NONE

    !LOCAL
    INTEGER :: jb
    REAL    :: scalefactor, numberconc, nsalt

    CALL initialize_indexarrays

    ! ------------------------------------------------------------------------
    SELECT CASE (APN)
    CASE (1)
      xaer(1)     = 1.
      ! CASE(1) can be used for either chamber or volcano aerosol:
      SELECT CASE (TRIM(init_scenario))
      CASE ('LAB') ! chamber aerosol:
      radius(1)   = 2.0E-7
      lwc(1)      = 1.2E-9
      csalt(1)    = 6.1 * lwc(1) * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
      !c0_FeOHp(1) = 0.073 * csalt(1)
      c0_FeCl2p(1) = 0.073 * csalt(1)
      c0_NH4p(1)  = 0.
      c0_HCO3m(1) = 0.
      c0_NO3m(1)  = 0.
      c0_Clm(1)   = 0.877 * csalt(1)
      c0_Brm(1)   = 0.05 * csalt(1)
      c0_Im(1)    = 0.
      c0_IO3m(1)  = 0.
      c0_SO4mm(1) = 0.
      c0_HSO4m(1) = 0.
      exchng(1)   = 0.
      CASE ('VOLCANO') ! volcano aerosol (values from N. Bobrowski, pers. comm.):
        radius(1)   = 1.4E-6
        numberconc  = 2E6 ! [m-3]
        lwc(1)      = numberconc * 4.*PI/3. * radius(1)**3
        ! divide nions by 2 because NaHSO4 produces 2 ions: Na+ and HSO4-:
        nsalt       = nions_koehler(radius(1),temp,relhum) / 2.
        csalt(1)    = nsalt * numberconc * N_A / 1.E6 ! [mcl/cm3(air)]
        c0_NH4p(1)  = 0.
        c0_HCO3m(1) = 0.
        c0_NO3m(1)  = 0.
        c0_Clm(1)   = 0.
        c0_Brm(1)   = 0.
        c0_Im(1)    = 0.
        c0_IO3m(1)  = 0.
        c0_SO4mm(1) = 0.
        c0_HSO4m(1) = csalt(1)
        exchng(1)   = 0.
      CASE DEFAULT
        PRINT *, 'ERROR, init_scenario '//TRIM(init_scenario)//' is not defined'
        STOP 1
      END SELECT
    CASE (2) ! all values for APN=2 are from MOCCA
      ! index 1 = sulfate aerosol
      xaer(1)     = 1.
      radius(1)   = 0.0882E-6
      lwc(1)      = 1.08E-12
      csalt(1)    = 3.74074074074 * lwc(1) * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
      c0_NH4p(1)  = csalt(1)
      c0_HCO3m(1) = 0.
      c0_NO3m(1)  = 0.
      c0_Clm(1)   = 0.
      c0_Brm(1)   = 0.
      c0_Im(1)    = 0.
      c0_IO3m(1)  = 0.
      c0_SO4mm(1) = 0.
      c0_HSO4m(1) = csalt(1)
      exchng(1)   = 1. / (7.*OneDay) ! exchange with fresh aerosol [1/s]
      ! index 2 = sea-salt aerosol
      IF (l_ff) THEN
        xaer(2)   = 0. ! frost flower aerosol will be switched on later
        radius(2) = 1.0E-6   ! assumed, see email 23 Dec 04
        IF (degree_lat>0.) THEN
          ! Arctic:
          lwc(2)    = 5.0E-10  ! assumed for fresh frost flower "plume"
        ELSE
          ! Antarctic:
          lwc(2)    = 3.0E-10  ! assumed for fresh frost flower "plume"
        ENDIF
        csalt(2)  = 5. * lwc(2) * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
      ELSE
        xaer(2)   = 1.
        radius(2) = 1.67E-6  ! radius [m]
        lwc(2)    = 3.04E-11 ! liquid water content [m3/m3]
        csalt(2)  = 5.3 * lwc(2) * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
      ENDIF
      c0_NH4p(2)  = 0.
      c0_HCO3m(2) = HCO3m_rel * csalt(2)
      c0_NO3m(2)  = 0.
      c0_Clm(2)   = Clm_rel   * csalt(2)
      c0_Brm(2)   = Brm_rel   * csalt(2)
      c0_Im(2)    = Im_rel    * csalt(2)
      c0_IO3m(2)  = IO3m_rel  * csalt(2)
      c0_SO4mm(2) = SO4mm_rel * csalt(2)
      c0_HSO4m(2) = 0.
      exchng(2)   = 1. / (2.*OneDay) ! exchange with fresh aerosol [1/s]
      ! ----------------------------------------------------------------------
    CASE (3) ! 1 and 2 = same as for CASE(2); 3 = ocean surface
      ! index 1 = sulfate aerosol
      xaer(1)     = 1.
      radius(1)   = 0.0882E-6
      lwc(1)      = 1.08E-12
      csalt(1)    = 3.74074074074 * lwc(1) * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
      c0_NH4p(1)  = csalt(1)
      c0_HCO3m(1) = 0.
      c0_NO3m(1)  = 0.
      c0_Clm(1)   = 0.
      c0_Brm(1)   = 0.
      c0_Im(1)    = 0.
      c0_IO3m(1)  = 0.
      c0_SO4mm(1) = 0.
      c0_HSO4m(1) = csalt(1)
      exchng(1)   = 1. / (7.*OneDay) ! exchange with fresh aerosol [1/s]
      ! index 2 = sea-salt aerosol
      xaer(2)   = 1.
      radius(2) = 1.67E-6  ! radius [m]
      lwc(2)    = 3.04E-11 ! liquid water content [m3/m3]
      csalt(2)  = 5.3 * lwc(2) * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
      c0_NH4p(2)  = 0.
      c0_HCO3m(2) = HCO3m_rel * csalt(2)
      c0_NO3m(2)  = 0.
      c0_Clm(2)   = Clm_rel   * csalt(2)
      c0_Brm(2)   = Brm_rel   * csalt(2)
      c0_Im(2)    = Im_rel    * csalt(2)
      c0_IO3m(2)  = IO3m_rel  * csalt(2)
      c0_SO4mm(2) = SO4mm_rel * csalt(2)
      c0_HSO4m(2) = 0.
      exchng(2)   = 1. / (2.*OneDay) ! exchange with fresh aerosol [1/s]
      ! index 3 = ocean surface
      xaer(3)   = 1.
      radius(3) = -999. ! radius [m] negative dummy value for ocean surface
      lwc(3)    = zmix / zmbl ! liquid water content [m3/m3]
      ! the concentration of sea water is about 0.53 mol/L
      csalt(3)  = 0.53 * lwc(3) * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
      c0_NH4p(3)  = 0.
      c0_HCO3m(3) = HCO3m_rel * csalt(3)
      c0_NO3m(3)  = NO3m_rel  * csalt(3)
      c0_Clm(3)   = Clm_rel   * csalt(3)
      c0_Brm(3)   = Brm_rel   * csalt(3)
      c0_Im(3)    = Im_rel    * csalt(3)
      c0_IO3m(3)  = IO3m_rel  * csalt(3)
      c0_SO4mm(3) = SO4mm_rel * csalt(3)
      c0_HSO4m(3) = 0.
      exchng(3)   = 1. / (2.*OneDay) ! exchange with deeper layers [1/s]

      ! ----------------------------------------------------------------------
    CASE (5)
      ! IMPORTANT: for testing purposes, the whole aerosol distribution
      ! has been enlarged by a scalefactor!!!
      scalefactor = 5.
      !
      ! index 1-2 = sulfate aerosol
      xaer(1:2)     = 1.
      ! derived from Ntot=2.80E+08 m-3, R=5.50E-08 m, and lg(sigma)=0.111
      ! radius [m]
      radius(1)  = 5.300E-08
      radius(2)  = 7.745E-08
      ! liquid water content [m3/m3]
      lwc(1)     = 9.412E-14 * scalefactor
      lwc(2)     = 1.397E-13 * scalefactor
      ! cations and anions
      csalt(1:2)    = 3.7 * lwc(1:2)*N_A/1.E3 ! [mol/L] conv to [mcl/cm3(air)]
      c0_NH4p(1:2)  = csalt(1:2)
      c0_HCO3m(1:2) = 0.
      c0_NO3m(1:2)  = 0.
      c0_Clm(1:2)   = 0.
      c0_Brm(1:2)   = 0.
      c0_Im(1:2)    = 0.
      c0_IO3m(1:2)  = 0.
      c0_SO4mm(1:2) = 0.
      c0_HSO4m(1:2) = csalt(1:2)
      exchng(1:2)   = 1. / (7.*OneDay) ! exchange with fresh aerosol [1/s]
      ! index 3-5 = sea-salt aerosol
      xaer(3:5)     = 1.
      ! derived from Ntot=6.60E+05 m-3, R=8.70E-07 m, and lg(sigma)=0.191
      ! radius [m]
      radius(3)  = 6.467E-07
      radius(4)  = 1.414E-06
      radius(5)  = 3.092E-06
      ! liquid water content [m3/m3]
      lwc(3)     = 3.938E-13 * scalefactor
      lwc(4)     = 3.066E-12 * scalefactor
      lwc(5)     = 1.865E-12 * scalefactor
      ! cations and anions
      csalt(3:5)    = 5.3 * lwc(3:5)*N_A/1.E3 ! [mol/L] conv to [mcl/cm3(air)]
      c0_NH4p(3:5)  = 0.
      c0_HCO3m(3:5) = HCO3m_rel * csalt(3:5)
      c0_NO3m(3:5)  = 0.
      c0_Clm(3:5)   = Clm_rel   * csalt(3:5)
      c0_Brm(3:5)   = Brm_rel   * csalt(3:5)
      c0_Im(3:5)    = Im_rel    * csalt(3:5)
      c0_IO3m(3:5)  = IO3m_rel  * csalt(3:5)
      c0_SO4mm(3:5) = SO4mm_rel * csalt(3:5)
      c0_HSO4m(3:5) = 0.
      exchng(3:5)   = 1. / (2.*OneDay) ! exchange with fresh aerosol [1/s]
      ! ----------------------------------------------------------------------
    CASE (8)
      ! IMPORTANT: for testing purposes, the whole aerosol distribution
      ! has been enlarged by a scalefactor!!!
      scalefactor = 1.
      !
      !---------------------------------
      ! define sea salt aerosol bin 1-7
      !---------------------------------
      xaer(1:7)     = 1.
      ! radius [m]
      radius(1)  = 0.230E-06
      radius(2)  = 0.445E-06
      radius(3)  = 0.900E-06
      radius(4)  = 1.650E-06
      radius(5)  = 3.550E-06
      radius(6)  = 7.500E-06
      radius(7)  = 14.50E-06
      ! liquid water content [m3/m3]
      lwc(1)     = 0.3E-12 * scalefactor
      lwc(2)     = 0.9E-12 * scalefactor
      lwc(3)     = 3.4E-12 * scalefactor
      lwc(4)     = 5.4E-12 * scalefactor
      lwc(5)     = 13.9E-12 * scalefactor
      lwc(6)     = 16.5E-12 * scalefactor
      lwc(7)     = 7.9E-12 * scalefactor
      ! cations and anions
      csalt(1:7)    = 5.3 * lwc(1:7)*N_A/1.E3 ! [mol/L] conv to [mcl/cm3(air)]
      c0_NH4p(1:7)  = 0.
      c0_HCO3m(1:7) = HCO3m_rel * csalt(1:7)
      c0_NO3m(1:7)  = 0.
      c0_Clm(1:7)   = Clm_rel   * csalt(1:7)
      c0_Brm(1:7)   = Brm_rel   * csalt(1:7)
      c0_Im(1:7)    = Im_rel    * csalt(1:7)
      c0_IO3m(1:7)  = IO3m_rel  * csalt(1:7)
      c0_SO4mm(1:7) = SO4mm_rel * csalt(1:7)
      c0_HSO4m(1:7) = 0.
      exchng(1:7)   = 1. / (2.*OneDay) ! exchange with fresh aerosol [1/s]
      ! TODO: calculate explicitly the dry deposition velocity and a typical
      ! mixing height for each aerosol size and obtain the aerosol lifetime
      ! from those
      !---------------------------------
      ! define sulfate aerosol bin 8
      !---------------------------------
      xaer(8)     = 1.
      ! radius [m]
      radius(8)   = 0.0882E-6
      ! liquid water content [m3/m3]
      lwc(8)      = 1.8E-12
      ! cations and anions
      csalt(8)    = 3.7 * lwc(8) * N_A / 1.E3 ! [mol/L], conv to [mcl/cm3(air)]
      c0_NH4p(8)  = csalt(8)
      c0_HCO3m(8) = 0.
      c0_NO3m(8)  = 0.
      c0_Clm(8)   = 0.
      c0_Brm(8)   = 0.
      c0_Im(8)    = 0.
      c0_IO3m(8)  = 0.
      c0_SO4mm(8) = 0.
      c0_HSO4m(8) = csalt(8)
      exchng(8)   = 1. / (7.*OneDay) ! exchange with fresh aerosol [1/s]
      ! ----------------------------------------------------------------------
    CASE (10)
      ! index 1-5 = sulfate aerosol
      xaer(1:5)     = 1.
      ! derived from Ntot=2.80E+08 m-3, R=5.50E-08 m, and lg(sigma)=0.111
      ! radius [m]
      radius(1)  = 3.626E-08
      radius(2)  = 5.300E-08
      radius(3)  = 7.745E-08
      radius(4)  = 1.132E-07
      radius(5)  = 1.654E-07
      ! liquid water content [m3/m3]
      lwc(1)     = 9.949E-15
      lwc(2)     = 9.412E-14
      lwc(3)     = 1.397E-13
      lwc(4)     = 3.093E-14
      lwc(5)     = 9.147E-16
      ! cations and anions
      csalt(1:5)    = 3.7 * lwc(1:5)*N_A/1.E3 ! [mol/L] conv to [mcl/cm3(air)]
      c0_NH4p(1:5)  = csalt(1:5)
      c0_HCO3m(1:5) = 0.
      c0_NO3m(1:5)  = 0.
      c0_Clm(1:5)   = 0.
      c0_Brm(1:5)   = 0.
      c0_Im(1:5)    = 0.
      c0_IO3m(1:5)  = 0.
      c0_SO4mm(1:5) = 0.
      c0_HSO4m(1:5) = csalt(1:5)
      exchng(1:5)   = 1. / (7.*OneDay) ! exchange with fresh aerosol [1/s]
      ! index 6-10 = sea-salt aerosol
      xaer(6:10)     = 1.
      ! derived from Ntot=6.60E+05 m-3, R=8.70E-07 m, and lg(sigma)=0.191
      ! radius [m]
      radius(6)  = 2.957E-07
      radius(7)  = 6.467E-07
      radius(8)  = 1.414E-06
      radius(9)  = 3.092E-06
      radius(10) = 6.762E-06
      ! liquid water content [m3/m3]
      lwc(6)     = 4.151E-15
      lwc(7)     = 3.938E-13
      lwc(8)     = 3.066E-12
      lwc(9)     = 1.865E-12
      lwc(10)    = 6.702E-14
      ! cations and anions
      csalt(6:10)    = 5.3 * lwc(6:10)*N_A/1.E3 ! [mol/L] conv to [mcl/cm3(air)]
      c0_NH4p(6:10)  = 0.
      c0_HCO3m(6:10) = HCO3m_rel * csalt(6:10)
      c0_NO3m(6:10)  = 0.
      c0_Clm(6:10)   = Clm_rel   * csalt(6:10)
      c0_Brm(6:10)   = Brm_rel   * csalt(6:10)
      c0_Im(6:10)    = Im_rel    * csalt(6:10)
      c0_IO3m(6:10)  = IO3m_rel  * csalt(6:10)
      c0_SO4mm(6:10) = SO4mm_rel * csalt(6:10)
      c0_HSO4m(6:10) = 0.
      exchng(6:10)   = 1. / (2.*OneDay) ! exchange with fresh aerosol [1/s]
      ! ----------------------------------------------------------------------
    ! CASE (12)
      ! ...
      ! ----------------------------------------------------------------------
    CASE DEFAULT
      PRINT *, 'ERROR: No aerosol definition available for ', APN, ' aerosol sizes.'
      STOP 1
    END SELECT
    ! ------------------------------------------------------------------------

    DO jb=1, APN
      IF (ind_H2O_a(jb)   /= 0) c(ind_H2O_a(jb))   = &
        rho_H2O * lwc(jb) * 1000./M_H2O * N_A/1.E6
      IF (ind_NH4p_a(jb)  /= 0) c(ind_NH4p_a(jb))  = c0_NH4p(jb)
      IF (ind_FeOHp_a(jb) /= 0) c(ind_FeOHp_a(jb)) = c0_FeOHp(jb)
      IF (ind_FeCl2p_a(jb)/= 0) c(ind_FeCl2p_a(jb))= c0_FeCl2p(jb)
      IF (ind_HCO3m_a(jb) /= 0) c(ind_HCO3m_a(jb)) = c0_HCO3m(jb)
      IF (ind_NO3m_a(jb)  /= 0) c(ind_NO3m_a(jb))  = c0_NO3m(jb)
      IF (ind_Clm_a(jb)   /= 0) c(ind_Clm_a(jb))   = c0_Clm(jb)
      IF (ind_Brm_a(jb)   /= 0) c(ind_Brm_a(jb))   = c0_Brm(jb)
      IF (ind_Im_a(jb)    /= 0) c(ind_Im_a(jb))    = c0_Im(jb)
      IF (ind_IO3m_a(jb)  /= 0) c(ind_IO3m_a(jb))  = c0_IO3m(jb)
      IF (ind_SO4mm_a(jb) /= 0) c(ind_SO4mm_a(jb)) = c0_SO4mm(jb)
      IF (ind_HSO4m_a(jb) /= 0) c(ind_HSO4m_a(jb)) = c0_HSO4m(jb)
      ! Na+ is a generic cation to keep the ion balance
      IF (ind_Nap_a(jb) /= 0) THEN
        c0_Nap(jb) = 0.
        IF (ind_NH4p_a(jb)  /= 0) c0_Nap(jb) = c0_Nap(jb) - c(ind_NH4p_a(jb))
        IF (ind_FeOHp_a(jb) /= 0) c0_Nap(jb) = c0_Nap(jb) - c(ind_FeOHp_a(jb))
        IF (ind_FeCl2p_a(jb)/= 0) c0_Nap(jb) = c0_Nap(jb) - c(ind_FeCl2p_a(jb))
        IF (ind_HCO3m_a(jb) /= 0) c0_Nap(jb) = c0_Nap(jb) + c(ind_HCO3m_a(jb))
        IF (ind_NO3m_a(jb)  /= 0) c0_Nap(jb) = c0_Nap(jb) + c(ind_NO3m_a(jb))
        IF (ind_Clm_a(jb)   /= 0) c0_Nap(jb) = c0_Nap(jb) + c(ind_Clm_a(jb))
        IF (ind_Brm_a(jb)   /= 0) c0_Nap(jb) = c0_Nap(jb) + c(ind_Brm_a(jb))
        IF (ind_Im_a(jb)    /= 0) c0_Nap(jb) = c0_Nap(jb) + c(ind_Im_a(jb))
        IF (ind_IO3m_a(jb)  /= 0) c0_Nap(jb) = c0_Nap(jb) + c(ind_IO3m_a(jb))
        IF (ind_SO4mm_a(jb) /= 0) c0_Nap(jb) = c0_Nap(jb) + 2.*c(ind_SO4mm_a(jb))
        IF (ind_HSO4m_a(jb) /= 0) c0_Nap(jb) = c0_Nap(jb) + c(ind_HSO4m_a(jb))
        c(ind_Nap_a(jb)) = c0_Nap(jb)
      ENDIF
    ENDDO

    ! cvfac: conversion factor dm^3(aq)/mol => cm^3(air)/molecule
    cvfac(:) = 1.E3 / ( N_A * lwc(:) )

    ! ------------------------------------------------------------------------

    ! print aerosol properties:
    PRINT *, HLINE2
    PRINT *, '              aerosol properties'
    PRINT *, HLINE2
    PRINT *, 'APN         r       LWC         N         A   1 mol/L =   1 mol/L ='
    PRINT *, '          [m]   [m3/m3]    [1/m3]   [m2/m3]   [mcl/cm3]   [mol/mol]'
    PRINT *, HLINE2
    DO jb=1, APN
      WRITE(*,'(A2,I2.2,4(ES10.2),2(ES12.4))') ' A', &
        jb, radius(jb), lwc(jb), &
        lwc(jb)*3./(4.*pi*radius(jb)**3), &  ! number
        3.*lwc(jb)/radius(jb), &             ! surface
        1. / cvfac(jb), &                    ! [mol/L] --> [mcl/cm3]
        lwc(jb) * 1E3 * R_gas * temp / press ! [mol/L] --> [mol/mol]
    ENDDO
    PRINT *, HLINE2
    PRINT *, '              aerosol anion concentrations [mol/L]'
    PRINT *, HLINE2
    PRINT *, 'APN    HCO3-     NO3-      Cl-      Br-       I-'// &
             '     IO3-    SO4--    HSO4-'
    PRINT *, HLINE2
    DO jb=1, APN
      WRITE(*,'(A2,I2.2,7(ES9.2))') ' A', jb, &
        c0_HCO3m(jb) *1E3/(lwc(jb)*N_A), &
        c0_NO3m(jb)  *1E3/(lwc(jb)*N_A), &
        c0_Clm(jb)   *1E3/(lwc(jb)*N_A), &
        c0_Brm(jb)   *1E3/(lwc(jb)*N_A), &
        c0_Im(jb)    *1E3/(lwc(jb)*N_A), &
        c0_IO3m(jb)  *1E3/(lwc(jb)*N_A), &
        c0_SO4mm(jb) *1E3/(lwc(jb)*N_A), &
        c0_HSO4m(jb) *1E3/(lwc(jb)*N_A)
    ENDDO
    PRINT *, HLINE2
    PRINT *, '              aerosol cation concentrations [mol/L]'
    PRINT *, HLINE2
    PRINT *, 'APN     NH4+      Na+'
    PRINT *, HLINE2
    DO jb=1, APN
      WRITE(*,'(A2,I2.2,2(ES9.2))') ' A', jb, &
        c0_NH4p(jb)  *1E3/(lwc(jb)*N_A), &
        c0_Nap(jb)   *1E3/(lwc(jb)*N_A)
    ENDDO

  END SUBROUTINE define_aerosol

  !***************************************************************************

  SUBROUTINE x0

    USE caaba_io,                 ONLY: nf90_inquire_dimension,              &
                                        nf90_inq_dimid, nf90_get_att
    USE caaba_mem,                ONLY: init_scenario, degree_lat, time0_jul
    USE messy_main_constants_mem, ONLY: OneDay
    USE messy_main_tools,         ONLY: ucase, mr2spechum, spec2relhum,      &
                                        spechum2mr, cair_c
    USE messy_main_timer,         ONLY: eval_time_str, gregor2julian

    IMPLICIT NONE

    INTRINSIC :: TRIM

    LOGICAL :: l_init_tpoint    ! init_spec only 1 time point?
    !INTEGER                     :: i, n_var, n_dim, varid_x, ct_spc
    INTEGER :: i, n_var, n_dim, varid_x, ct_spc, tpoint, len_time, dimid_time
    !INTEGER   :: status
    INTEGER :: status, varid_time, st0year, st0month, st0day, st0hour,       &
               st0min, st0sec
    !REAL(DP) :: mr_x
    REAL(DP) :: mr_x, mr_x1, ratio1, ratio2, time1, time2, stuf, stime0_jul
    CHARACTER(LEN=33) :: time_string     ! 'seconds since 2000-01-01 00:00:00'
    CHARACTER(LEN=STRLEN_VLONG) :: name_x, name_spc ! chem species' names
    CHARACTER(LEN=STRLEN_VLONG) :: timename ! time name
    
    ! ------------------------------------------------------------------------

    !> initialize some mixing ratios
    !> values in mol/mol, cair converts to particles/cm3

    SELECT CASE (TRIM(init_scenario))
    CASE ('') ! default if no scenario selected in namelist:
      CALL x0_simple
    CASE ('ZEROAIR')
      CALL x0_zeroair
    CASE ('CUMULUS')
      CALL x0_cumulus
    CASE ('FF_ANTARCTIC')
      CALL x0_ff_antarctic
    CASE ('FF_ARCTIC')
      CALL x0_ff_arctic
    CASE ('FREE_TROP')
      CALL x0_free_trop
    CASE ('HOOVER')
      CALL x0_hoover
    CASE ('LAB','LAB_C15')
      CALL x0_lab
    CASE ('MBL')
      CALL x0_mbl
    CASE ('MOM')
      CALL x0_mom
    CASE ('MTCHEM')
      CALL x0_mtchem
    CASE ('OOMPH')
      CALL x0_oomph
    CASE ('ISO')
      CALL x0_iso
    CASE ('TAG')
      CALL x0_tag
    CASE ('STRATO')
      ! choose one:
      CALL x0_strato10
      !CALL x0_strato20
    ! op_ff_20160309+
    CASE ('TROPOPAUSE')
       CALL x0_tropopause
    CASE ('LOW_STRATO')
       CALL x0_low_strato
    CASE ('MID_STRATO')
       CALL x0_mid_strato
    CASE ('HIGH_STRATO')
       CALL x0_high_strato
    ! op_ff_20160309-
    CASE ('VOLCANO')
      CALL x0_volcano
    CASE DEFAULT
      PRINT *, 'ERROR, init_scenario '//TRIM(init_scenario)//' is not defined'
      STOP 1
    END SELECT

    ! ------------------------------------------------------------------------

    ! read init at model_time (not just first value)
    !> initialization of chemical species from external file
    IF (TRIM(init_spec)/="") THEN
      PRINT *, HLINE2
      PRINT *, 'External chemical initialization:'
      PRINT *, HLINE2

      ! default: init_spec contains just 1 point in time
      l_init_tpoint = .TRUE.

      ! open init_spec file, get file ID
      CALL open_input_file(ncid_spec, init_spec)
      ! get ID of unlimited dimension:
      CALL nf(nf90_inquire(ncid_spec, unlimitedDimId = dimid_time))
      IF (dimid_time == -1) THEN
        ! there is no unlimited dimension, thus only 1 dataset:
        len_time = 1
      ELSE
        ! get name and length of unlimited dimension:
        CALL nf(nf90_inquire_dimension(ncid_spec, dimid_time, len = len_time, &
          name = timename), "inquire_dimension for dimid_time")
      ENDIF
      ! if more than one time point, consider time axis of init_spec
      IF (len_time > 1) THEN
        PRINT *, "Time axis points considered: len_time = ", len_time
        WRITE(*,*) "Interpolating chemical initialization to model start"
        ! init_spec contains > 1 time point
        l_init_tpoint = .FALSE.

        ! get time variable ID and time units
        CALL nf(nf90_inq_varid(ncid_spec, TRIM(timename), varid_time))
        CALL nf(nf90_get_att(ncid_spec, varid_time, "units", time_string))
        ! determine time unit conversion factor to seconds
        CALL eval_time_str(status, time_string, stuf, &
               st0year, st0month, st0day, st0hour, st0min, st0sec)
        IF (status/=0) THEN
          WRITE(*,*) "ERROR in x0: eval_time_str"
          STOP 1
        ENDIF

        ! convert init_spec start into Julian day
        stime0_jul = gregor2julian(st0year, st0month, st0day, st0hour,&
          st0min, st0sec)

        ! check if init_spec includes model start time
        CALL nf(nf90_get_var(ncid_spec, varid_time, time1, &
          start = (/ 1 /)))
        time1 = time1 * stuf + stime0_jul*OneDay - time0_jul*OneDay
        CALL nf(nf90_get_var(ncid_spec, varid_time, time2, &
          start = (/ len_time /)))
        time2 = time2 * stuf + stime0_jul*OneDay - time0_jul*OneDay
        IF (time1 > model_time .OR. time2 < model_time) THEN
          WRITE(*,*) "ERROR x0: init_spec does not cover model start time"
          WRITE(*,*) "model start:     ", model_time
          WRITE(*,*) "init_spec start: ", time1
          WRITE(*,*) "init_spec end:   ", time2
          STOP 1
        ENDIF

        ! find out between which time points starting point is
        DO tpoint = 2, len_time, 1
          CALL nf(nf90_get_var(ncid_spec, varid_time, time2, &
            start = (/ tpoint /)))
          ! transform time to seconds with general time origin
          time2 = time2 * stuf + stime0_jul*OneDay - time0_jul*OneDay

          IF (time2 >= model_time) THEN
            ! tpoint-1 <= model_time < tpoint
            CALL nf(nf90_get_var(ncid_spec, varid_time, time1, &
              start = (/ tpoint-1 /)))
            time1 = time1 * stuf + stime0_jul*OneDay - time0_jul*OneDay

            ! calculate interpolation factors
            ratio2 = (model_time - time1)/(time2 - time1)
            ratio1 = 1._dp - ratio2
            print *, "mo_time = ", model_time
            print *, "time1   = ", time1
            print *, "time2   = ", time2
            print *, "ratio1  = ", ratio1
            print *, "ratio2  = ", ratio2
            EXIT ! exit loop when interpolation weights found
          ENDIF ! time2 >= model_time
        ENDDO ! loop over time points of init_spec file
      ENDIF ! > 1 time point in init_spec

      !> initialize chemical species from external file
      ! get number of dimensions and variables (=specs)
      CALL nf(nf90_inquire(ncid_spec, n_dim, n_var))

      varid_x = 1
      DO WHILE (varid_x <= n_var) ! loop over ext init species
        ! get name of external variable
        CALL nf(nf90_inquire_variable(ncid_spec, varid_x, name_x))
        ! convert name to uppercase for comparison
        CALL ucase(name_x)
        ct_spc = 1
        DO WHILE (ct_spc <= NSPEC) ! loop over chemical species (kpp)
          name_spc = SPC_NAMES(ct_spc)
          ! convert to uppercase for comparison:
          CALL ucase(name_spc)
          ! if external spec found in meccanism
          IF (TRIM(name_spc) == TRIM(name_x)) THEN
            ! get external mixing ratio at time point after model_time
            IF (l_init_tpoint) THEN
              CALL nf(nf90_get_var(ncid_spec, varid_x, mr_x))
            ELSE
              ! interpolate on init_spec time axis
              ! get external mixing ratio at time point after model_time
              CALL nf(nf90_get_var(ncid_spec, varid_x, mr_x, &
                start = (/ tpoint /)))
              ! get external mixing ratio at time point before model_time
              CALL nf(nf90_get_var(ncid_spec, varid_x, mr_x1, &
                start = (/ tpoint-1 /)))
              ! use ratios to interpolate mixing ratio between points
              mr_x = ratio1 * mr_x1 + ratio2 * mr_x
            ENDIF

            ! put overwrite warning in same line
            IF (c(ct_spc)>0.) THEN
              WRITE(*,'(A,A16,ES10.2,A,ES10.2,A,ES10.2,A)') &
                '     initialized: ', &
                name_x, mr_x, ' mol/mol   = ', mr_x*cair, &
                ' mcl/cm3   (overwriting', c(ct_spc)/cair, ' mol/mol)'
            ELSE
              WRITE(*,'(A,A16,ES10.2,A,ES10.2,A)') &
                '     initialized: ', &
                name_x, mr_x, ' mol/mol   = ', mr_x*cair, ' mcl/cm3'
            ENDIF
            c(ct_spc) = mr_x * cair
            EXIT
          ELSE
            ct_spc = ct_spc + 1
            IF (ct_spc .GT. NSPEC) THEN
              WRITE(*,'(A,A16)') ' NOT_initialized: ', name_x
            ENDIF
          ENDIF
        ENDDO ! loop over chem spec
        varid_x = varid_x + 1
      ENDDO ! extloop
      CALL close_file(ncid_spec)
    ENDIF

    ! after x0: c(H2O) -> update cair, spechum, relhum
    IF (l_ignore_relhum) THEN
      ! update cair with current c(H2O) after x0
      cair = cair_c(c(ind_H2O), temp, press, l_hum_emac)

      ! update concentrations with new cair
      C(:) = C(:) * cair/cair_old
      cair_old = cair

      ! update spechum and relhum with new c(H2O) and cair
      spechum = mr2spechum(status, c(ind_H2O)/cair)
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR in x0: mr2spechum'
        STOP 1
      ENDIF

      relhum  = spec2relhum(status, spechum, temp, press, &
                            l_hum_emac, l_psat_liquid, l_relhum_wmo)
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR in x0: spec2relhum'
        STOP 1
      ENDIF
    ELSE
      ! convert from/to spechum
      !> convert specific humidity to water vapor concentration,
      !>   overwrite initialized c(H2O)
      !c(ind_H2O) = cair * &
      !  rh2mr(status,relhum,temp,press,l_hum_emac,l_relhum_wmo)
      c(ind_H2O) = cair * spechum2mr(status,spechum)
      IF (status > 1) THEN
        WRITE(*,*) "ERROR in x0: spechum2mr"
        STOP 1
      ENDIF
      !IF (status /= 0) STOP 1
    ENDIF

    ! ------------------------------------------------------------------------

    PRINT *, HLINE2
    PRINT *, 'Initial gas-phase mixing ratios and concentrations:'
    PRINT *, HLINE2
    DO i = 1,NSPEC
      IF (c(i)>0) THEN
        WRITE(*,'(2A,ES10.2,A,ES10.2,A)') ' ', SPC_NAMES(i), &
          c(i)/cair, ' mol/mol   = ', c(i), ' mcl/cm3'
      ENDIF
    ENDDO

    ! ------------------------------------------------------------------------

  CONTAINS

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_simple
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_CO       /= 0) c(ind_CO)      =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
    END SUBROUTINE x0_simple

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_zeroair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
    END SUBROUTINE x0_zeroair

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_cumulus
    ! cloud med: 10s CARIBIC data, Flight353, 16 Aug 2011 18:20-18:35 UTC
    ! bg = "fake" bg from over continental U.S.
    ! bg2 = proposed bg for NO, HCHO, CH3OH
    ! bg3 = NO equil lower at 153 ppt
    ! bg4 = paper revision
      IF (ind_NO       /= 0) c(ind_NO)       =    0.440E-09 * cair ! bg4 loj start -> reach 130 ppt at equil
      !IF (ind_NO       /= 0) c(ind_NO)       =    1.316E-09 * cair ! cloud med
      !IF (ind_NO       /= 0) c(ind_NO)       =    0.350E-09 * cair ! bg2 -> reach 193 ppt at equil
      !IF (ind_NO       /= 0) c(ind_NO)       =    0.280E-09 * cair ! bg3 -> reach 153 ppt at equil
      !IF (ind_H2O      /= 0) c(ind_H2O)      =  187.185E-06 * cair ! cloud med
      IF (ind_H2O      /= 0) c(ind_H2O)      =  108.000E-06 * cair ! bg4 loj fix
      IF (ind_O3       /= 0) c(ind_O3)       =   38.000E-09 * cair ! bg4 loj fix + cloud med revised
      !IF (ind_O3       /= 0) c(ind_O3)       =   35.265E-09 * cair ! cloud med
      !IF (ind_O3       /= 0) c(ind_O3)       =   66.000E-09 * cair ! bg (U.S. influence?)
      !IF (ind_CO       /= 0) c(ind_CO)      =   73.050E-09 * cair ! cloud med orig
      IF (ind_CO       /= 0) c(ind_CO)       =   73.000E-09 * cair ! bg4 loj fix cloud med
      IF (ind_CH4      /= 0) c(ind_CH4)      = 1807.050E-09 * cair ! bg4 loj fix, cloud med
      !IF (ind_CH4      /= 0) c(ind_CH4)      = 1814.000E-09 * cair ! bg (U.S. influence?)
      !IF (ind_C5H8     /= 0) c(ind_C5H8)     =   50.000E-12 * cair
      !IF (ind_CO2      /= 0) c(ind_CO2)     =  389.615E-06 * cair ! cloud med orig
      IF (ind_CO2      /= 0) c(ind_CO2)      =  390.000E-06 * cair ! bg4 loj fix + cloud med
    ! from Neumaier KIT [Crawford1999b]: methanol = 1.5*acetone, here: factor 3.3
      IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3) =  300.000E-12 * cair ! bg4 loj fix + cloud med (noisy data)
      !IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3) =  665.000E-12 * cair ! bg (U.S. influence?)
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)    =  929.500E-12 * cair ! bg4 loj fix cloud med and before cloud
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)    =  400.000E-12 * cair ! (U.S. influence?)
    ! from WAS
      IF (ind_N2O      /= 0) c(ind_N2O)      =  325.830E-09 * cair ! fix, cloud med
      IF (ind_C2H2     /= 0) c(ind_C2H2)     =   25.140E-12 * cair ! bg4 loj fix, cloud med
      IF (ind_C2H6     /= 0) c(ind_C2H6)     =  400.010E-12 * cair ! bg4 loj fix, cloud med
      !IF (ind_C2H6     /= 0) c(ind_C2H6)     =  790.010E-12 * cair ! bg (U.S.  !influence?)
      IF (ind_C3H8     /= 0) c(ind_C3H8)     =   10.760E-12 * cair ! bg4 loj fix, cloud med
      !IF (ind_C3H8     /= 0) c(ind_C3H8)     =   88.000E-12 * cair ! bg (U.S.  !influence?)
      IF (ind_NC4H10   /= 0) c(ind_NC4H10)   =    1.850E-12 * cair ! bg4 loj fix !cloud med
      !not in meccanism IF (ind_SF6      /= 0) c(ind_SF6)      =    7.420E-12 * cair ! cloud med
      !not in meccanism IF (ind_CH3Cl    /= 0) c(ind_CH3Cl)    =  563.880E-12 * cair ! cloud med
      !not in meccanism IF (ind_CH3Cl    /= 0) c(ind_CH3Cl)    =  600.000E-12 * cair ! bg
    ! from DOAS
      IF (ind_HONO     /= 0) c(ind_HONO)     =    1.000E-12 * cair ! bg4 loj start
      !IF (ind_HONO     /= 0) c(ind_HONO)     =   37.000E-12 * cair ! cloud med
      !IF (ind_HONO     /= 0) c(ind_HONO)     =   46.000E-12 * cair ! cloud med orig
      IF (ind_HCHO     /= 0) c(ind_HCHO)     =   70.000E-12 * cair ! bg4 loj start
      !IF (ind_HCHO     /= 0) c(ind_HCHO)     =  500.000E-12 * cair ! cloud med orig
      !IF (ind_HCHO     /= 0) c(ind_HCHO)     =  400.000E-12 * cair ! cloud med
      !IF (ind_HCHO     /= 0) c(ind_HCHO)     =  130.000E-12 * cair ! new cloud med DOAS + profile EMAC
      !IF (ind_HCHO     /= 0) c(ind_HCHO)     =   50.000E-12 * cair ! bg2 [Singh2009]
      !IF (ind_NO2      /= 0) c(ind_NO2)      =  220.000E-12 * cair ! cloud med
      !IF (ind_NO2      /= 0) c(ind_NO2)      =   20.000E-12 * cair ! bg2 -> reach slightly below 200 ppt NO at equil
      !IF (ind_NO2      /= 0) c(ind_NO2)      =   20.000E-12 * cair ! bg equil
      !IF (ind_NO2      /= 0) c(ind_NO2)      =   14.000E-12 * cair ! bg DOAS/Scia
    ! additional for mecca
      !not in meccanism IF (ind_SO2      /= 0) c(ind_SO2)      =    0.000 * cair
      IF (ind_O2       /= 0) c(ind_O2)       =  210.000E-03 * cair ! fix
      IF (ind_N2       /= 0) c(ind_N2)       =  780.000E-03 * cair ! fix
    ! bg2 values from [Mari2000], [Fried2008b], [Snow2007] SONEX: all 30 ppt
    ! CH3OOH bg, Snow 80 ppt H2O2 bg
      !IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)   =   600.000E-12 * cair ! fixed to get 400 ppt HCHO
      !IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)   =   490.000E-12 * cair ! fixed to get 400 ppt HCHO + NO emis
    ! from [Batenburg2012]: lowest value at 16 degN: 510 ppb
      IF (ind_H2 /= 0)       c(ind_H2)       =  510.000E-09 * cair ! fix
    ! [Keim2008] tropical UT over Brazil
      !let equil instead IF (ind_HNO3 /= 0)     c(ind_HNO3)     =  400.000E-12 * cair
    ! from equilibration in this study
      IF (ind_PAN /= 0)      c(ind_PAN)      =   80.000E-12 * cair ! bg4 loj start
      !IF (ind_PAN /= 0)      c(ind_PAN)      =   54.000E-12 * cair ! NO emis
      !IF (ind_PAN /= 0)      c(ind_PAN)      =   85.000E-12 * cair ! bg equil
      !IF (ind_HCOOH /= 0)    c(ind_HCOOH)    =    9.000E-12 * cair ! bg equil
    END SUBROUTINE x0_cumulus

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_ff_antarctic
      IF (ind_O3       /= 0) c(ind_O3)      =  30.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_NO       /= 0) c(ind_NO)      =  2.0E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     =  2.0E-12 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    =  50.E-12 * cair
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)  =  10.E-12 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 170.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_DMS      /= 0) c(ind_DMS)     =  10.E-12 * cair
      IF (ind_SO2      /= 0) c(ind_SO2)     =  30.E-12 * cair
      IF (ind_CH3I     /= 0) c(ind_CH3I)    =  2.0E-12 * cair
      IF (ind_CH3Br    /= 0) c(ind_CH3Br)   =  5.0E-12 * cair
      IF (ind_C2H4     /= 0) c(ind_C2H4)    =  10.E-12 * cair
      IF (ind_C2H2     /= 0) c(ind_C2H2)    =  10.E-12 * cair
      IF (ind_C2H6     /= 0) c(ind_C2H6)    = 300.E-12 * cair
      IF (ind_Hg       /= 0) c(ind_Hg)      = 1.68E-13 * cair
    END SUBROUTINE x0_ff_antarctic

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_ff_arctic
      IF (ind_O3       /= 0) c(ind_O3)      =  40.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_NO       /= 0) c(ind_NO)      =  10.E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     =  10.E-12 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 200.E-12 * cair
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)  = 100.E-12 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 170.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_DMS      /= 0) c(ind_DMS)     =  10.E-12 * cair
      IF (ind_SO2      /= 0) c(ind_SO2)     = 100.E-12 * cair
      IF (ind_CH3I     /= 0) c(ind_CH3I)    =  2.0E-12 * cair
      IF (ind_CH3Br    /= 0) c(ind_CH3Br)   =  5.0E-12 * cair
      IF (ind_C2H4     /= 0) c(ind_C2H4)    =   26E-12 * cair ! ref1737
      ! see also: C2H4 = 100E-12 ! ref0351, Tab.1, 3 Apr
      IF (ind_C2H2     /= 0) c(ind_C2H2)    =  329E-12 * cair ! ref1737
      ! see also: C2H2 = 840E-12 ! ref0351, Tab.1, 3 Apr
      IF (ind_C2H6     /= 0) c(ind_C2H6)    =  2.0E-09 * cair
      IF (ind_Hg       /= 0) c(ind_Hg)      = 1.68E-13 * cair
    END SUBROUTINE x0_ff_arctic

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_hoover
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.8E-06 * cair
      IF (ind_H2       /= 0) c(ind_H2)      = 0.6E-06 * cair
    END SUBROUTINE x0_hoover

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_free_trop
      ! average init from CARIBIC2 trajectories 5-day back init
      ! flight: 20060706_CAN_FRA_157_EMAC
      IF (ind_CH3OH       /= 0) c(ind_CH3OH)       =  7.23043E-10 * cair
      IF (ind_HCHO        /= 0) c(ind_HCHO)        =  2.53104E-10 * cair
      IF (ind_CO          /= 0) c(ind_CO)          =  7.29323E-08 * cair
      IF (ind_HCOOH       /= 0) c(ind_HCOOH)       =  4.13365E-11 * cair
      IF (ind_CH3CHO      /= 0) c(ind_CH3CHO)      =  2.00554E-11 * cair
      IF (ind_CH3CO2H     /= 0) c(ind_CH3CO2H)     =  6.93619E-11 * cair
      IF (ind_CH3COCH3    /= 0) c(ind_CH3COCH3)    =  3.34549E-10 * cair
      IF (ind_ACETOL      /= 0) c(ind_ACETOL)      =  5.64005E-11 * cair
      IF (ind_MGLYOX      /= 0) c(ind_MGLYOX)      =  2.03489E-11 * cair
      IF (ind_BIACET      /= 0) c(ind_BIACET)      =  4.19846E-13 * cair
      IF (ind_Br          /= 0) c(ind_Br)          =  3.38058E-14 * cair
      IF (ind_Br2         /= 0) c(ind_Br2)         =  3.87407E-16 * cair
      IF (ind_BrO         /= 0) c(ind_BrO)         =  2.39308E-13 * cair
      IF (ind_HBr         /= 0) c(ind_HBr)         =  5.70085E-13 * cair
      IF (ind_HOBr        /= 0) c(ind_HOBr)        =  2.86479E-13 * cair
      IF (ind_BrNO2       /= 0) c(ind_BrNO2)       =  0.          * cair
      IF (ind_BrNO3       /= 0) c(ind_BrNO3)       =  6.53232E-13 * cair
      IF (ind_BrCl        /= 0) c(ind_BrCl)        =  9.05211E-17 * cair
      IF (ind_Cl          /= 0) c(ind_Cl)          =  7.12196E-17 * cair
      IF (ind_Cl2         /= 0) c(ind_Cl2)         =  3.46324E-17 * cair
      IF (ind_ClO         /= 0) c(ind_ClO)         =  8.18961E-14 * cair
      IF (ind_HCl         /= 0) c(ind_HCl)         =  5.60562E-11 * cair
      IF (ind_HOCl        /= 0) c(ind_HOCl)        =  2.5891E-13  * cair
      IF (ind_Cl2O2       /= 0) c(ind_Cl2O2)       =  3.47767E-17 * cair
      IF (ind_OClO        /= 0) c(ind_OClO)        =  1.33148E-16 * cair
      IF (ind_ClNO3       /= 0) c(ind_ClNO3)       =  3.81274E-12 * cair
      IF (ind_CCl4        /= 0) c(ind_CCl4)        =  8.9102E-11  * cair
      IF (ind_CH3Cl       /= 0) c(ind_CH3Cl)       =  4.55034E-10 * cair
      IF (ind_CH3CCl3     /= 0) c(ind_CH3CCl3)     =  1.50035E-11 * cair
      IF (ind_CF2Cl2      /= 0) c(ind_CF2Cl2)      =  7.87179E-10 * cair
      IF (ind_CFCl3       /= 0) c(ind_CFCl3)       =  2.42555E-10 * cair
      IF (ind_CH2ClBr     /= 0) c(ind_CH2ClBr)     =  9.40689E-14 * cair
      IF (ind_CHCl2Br     /= 0) c(ind_CHCl2Br)     =  5.51093E-14 * cair
      IF (ind_CHClBr2     /= 0) c(ind_CHClBr2)     =  4.54419E-14 * cair
      IF (ind_CH2Br2      /= 0) c(ind_CH2Br2)      =  1.02541E-12 * cair
      IF (ind_CH3Br       /= 0) c(ind_CH3Br)       =  8.78744E-15 * cair
      IF (ind_CHBr3       /= 0) c(ind_CHBr3)       =  4.8533E-13  * cair
      IF (ind_CF3Br       /= 0) c(ind_CF3Br)       =  2.43175E-12 * cair
      IF (ind_CF2ClBr     /= 0) c(ind_CF2ClBr)     =  4.20382E-12 * cair
      IF (ind_CH4         /= 0) c(ind_CH4)         =  1.74813E-06 * cair
      IF (ind_C2H6        /= 0) c(ind_C2H6)        =  5.04266E-10 * cair
      IF (ind_C2H4        /= 0) c(ind_C2H4)        =  1.51334E-11 * cair
      IF (ind_C3H8        /= 0) c(ind_C3H8)        =  4.6215E-11  * cair
      IF (ind_C3H6        /= 0) c(ind_C3H6)        =  1.68475E-12 * cair
      IF (ind_NC4H10      /= 0) c(ind_NC4H10)      =  7.76293E-11 * cair
      IF (ind_MVK         /= 0) c(ind_MVK)         =  3.90326E-11 * cair
      IF (ind_MEK         /= 0) c(ind_MEK)         =  3.98487E-11 * cair
      IF (ind_C5H8        /= 0) c(ind_C5H8)        =  8.45986E-12 * cair
      IF (ind_HgO         /= 0) c(ind_HgO)         =  3.57349E-16 * cair
      IF (ind_HgBr2       /= 0) c(ind_HgBr2)       =  9.89811E-16 * cair
      IF (ind_ClHgBr      /= 0) c(ind_ClHgBr)      =  1.13395E-16 * cair
      IF (ind_BrHgOBr     /= 0) c(ind_BrHgOBr)     =  7.30999E-15 * cair
      IF (ind_ClHgOBr     /= 0) c(ind_ClHgOBr)     =  1.1631E-15  * cair
      IF (ind_HgCl        /= 0) c(ind_HgCl)        =  6.52374E-17 * cair
      IF (ind_HgBr        /= 0) c(ind_HgBr)        =  7.25411E-16 * cair
      IF (ind_Hg          /= 0) c(ind_Hg)          =  1.75258E-13 * cair
      IF (ind_NACA        /= 0) c(ind_NACA)        =  2.85441E-12 * cair
      IF (ind_MPAN        /= 0) c(ind_MPAN)        =  5.31725E-12 * cair
      IF (ind_IC3H7NO3    /= 0) c(ind_IC3H7NO3)    =  8.94082E-13 * cair
      IF (ind_LC4H9NO3    /= 0) c(ind_LC4H9NO3)    =  8.32727E-12 * cair
      IF (ind_ISON        /= 0) c(ind_ISON)        =  9.0462E-12  * cair
      IF (ind_N2          /= 0) c(ind_N2)          =  0.78        * cair
      IF (ind_NH3         /= 0) c(ind_NH3)         =  1.2068E-10  * cair
      IF (ind_N2O         /= 0) c(ind_N2O)         =  3.16661E-07 * cair
      IF (ind_NO          /= 0) c(ind_NO)          =  4.65026E-11 * cair
      IF (ind_NO2         /= 0) c(ind_NO2)         =  1.25056E-10 * cair
      IF (ind_NO3         /= 0) c(ind_NO3)         =  3.1458E-12  * cair
      IF (ind_N2O5        /= 0) c(ind_N2O5)        =  4.37511E-12 * cair
      IF (ind_HONO        /= 0) c(ind_HONO)        =  6.12167E-13 * cair
      IF (ind_HNO3        /= 0) c(ind_HNO3)        =  1.33769E-10 * cair
      IF (ind_HNO4        /= 0) c(ind_HNO4)        =  2.77194E-11 * cair
      IF (ind_PAN         /= 0) c(ind_PAN)         =  3.1475E-10  * cair
      IF (ind_NH2OH       /= 0) c(ind_NH2OH)       =  1.31453E-12 * cair
      IF (ind_NHOH        /= 0) c(ind_NHOH)        =  1.00251E-11 * cair
      IF (ind_HNO         /= 0) c(ind_HNO)         =  9.47677E-14 * cair
      IF (ind_NH2         /= 0) c(ind_NH2)         =  8.9126E-18  * cair
      IF (ind_O3P         /= 0) c(ind_O3P)         =  2.31275E-15 * cair
      IF (ind_O2          /= 0) c(ind_O2)          =  0.21        * cair
      IF (ind_O3          /= 0) c(ind_O3)          =  1.60052E-07 * cair ! also lower strat in average
      IF (ind_H2          /= 0) c(ind_H2)          =  5.35521E-07 * cair
      IF (ind_OH          /= 0) c(ind_OH)          =  7.07567E-14 * cair
      IF (ind_HO2         /= 0) c(ind_HO2)         =  9.74215E-13 * cair
      IF (ind_H2O2        /= 0) c(ind_H2O2)        =  6.10601E-10 * cair
      IF (ind_H2O         /= 0) c(ind_H2O)         =  0.0038177   * cair
      IF (ind_CH3O2       /= 0) c(ind_CH3O2)       =  1.1287E-12  * cair
      IF (ind_CH3OOH      /= 0) c(ind_CH3OOH)      =  2.37819E-10 * cair
      IF (ind_C2H5O2      /= 0) c(ind_C2H5O2)      =  1.21222E-14 * cair
      IF (ind_C2H5OOH     /= 0) c(ind_C2H5OOH)     =  9.14669E-12 * cair
      IF (ind_CH3CO3      /= 0) c(ind_CH3CO3)      =  9.65657E-14 * cair
      IF (ind_CH3CO3H     /= 0) c(ind_CH3CO3H)     =  5.72797E-11 * cair
      IF (ind_IC3H7O2     /= 0) c(ind_IC3H7O2)     =  3.568E-15   * cair
      IF (ind_IC3H7OOH    /= 0) c(ind_IC3H7OOH)    =  1.98726E-12 * cair
      IF (ind_LHOC3H6O2   /= 0) c(ind_LHOC3H6O2)   =  2.6452E-14  * cair
      IF (ind_LHOC3H6OOH  /= 0) c(ind_LHOC3H6OOH)  =  4.19966E-12 * cair
      IF (ind_CH3COCH2O2  /= 0) c(ind_CH3COCH2O2)  =  4.96396E-15 * cair
      IF (ind_HYPERACET   /= 0) c(ind_HYPERACET)   =  1.95516E-12 * cair
      IF (ind_LC4H9O2     /= 0) c(ind_LC4H9O2)     =  1.79502E-14 * cair
      IF (ind_LC4H9OOH    /= 0) c(ind_LC4H9OOH)    =  9.50258E-12 * cair
      IF (ind_MVKO2       /= 0) c(ind_MVKO2)       =  7.3364E-14  * cair
      IF (ind_MVKOOH      /= 0) c(ind_MVKOOH)      =  2.49326E-11 * cair
      IF (ind_LMEKO2      /= 0) c(ind_LMEKO2)      =  6.33844E-15 * cair
      IF (ind_LMEKOOH     /= 0) c(ind_LMEKOOH)     =  3.52573E-12 * cair
      IF (ind_ISO2        /= 0) c(ind_ISO2)        =  2.85094E-14 * cair
      IF (ind_ISOOH       /= 0) c(ind_ISOOH)       =  7.02518E-12 * cair
      IF (ind_SO2         /= 0) c(ind_SO2)         =  6.63355E-11 * cair
      IF (ind_H2SO4       /= 0) c(ind_H2SO4)       =  5.51045E-14 * cair
      IF (ind_CH3SO3H     /= 0) c(ind_CH3SO3H)     =  6.24971E-11 * cair
      IF (ind_DMS         /= 0) c(ind_DMS)         =  5.13875E-12 * cair
      IF (ind_DMSO        /= 0) c(ind_DMSO)        =  7.64715E-14 * cair
      IF (ind_CH3SO2      /= 0) c(ind_CH3SO2)      =  4.44594E-17 * cair
      IF (ind_CH3SO3      /= 0) c(ind_CH3SO3)      =  1.33623E-14 * cair
      IF (ind_CO2         /= 0) c(ind_CO2)         =  0.000382388 * cair
      IF (ind_CH3I        /= 0) c(ind_CH3I)        =  4.39702E-14 * cair
    END SUBROUTINE x0_free_trop

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_lab
      IF (ind_O2       /= 0) c(ind_O2)      = 210E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780E-03 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     = 3.5E-09 * cair
      IF (ind_DMP      /= 0) c(ind_DMP)     = 22.E-09 * cair
      IF (ind_DMB      /= 0) c(ind_DMB)     = 16.E-09 * cair
      IF (ind_TM5      /= 0) c(ind_TM5)     = 12.E-09 * cair
      IF (ind_TOLUENE  /= 0) c(ind_TOLUENE) = 19.E-09 * cair
      !IF (ind_O3      /= 0) c(ind_O3)      = 670E-09 * cair
      !IF (ind_CH4     /= 0) c(ind_CH4)     = 500E-09 * cair
    END SUBROUTINE x0_lab

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_mbl
      IF (ind_H2       /= 0) c(ind_H2)      =   1.E-06 * cair
      IF (ind_O3       /= 0) c(ind_O3)      =  25.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      !IF (ind_NO      /= 0) c(ind_NO)      =  10.E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     =  20.E-12 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 300.E-12 * cair
      !IF (ind_CH3CHO  /= 0) c(ind_CH3CHO)  =  1.0E-10 * cair
      IF (ind_CO       /= 0) c(ind_CO)      =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_DMS      /= 0) c(ind_DMS)     =  60.E-12 * cair
      IF (ind_SO2      /= 0) c(ind_SO2)     =  90.E-12 * cair
      IF (ind_CH3I     /= 0) c(ind_CH3I)    =  2.0E-12 * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)    = 600.E-12 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)     = 200.E-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)    =  5.0E-12 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)     =  40.E-12 * cair
      IF (ind_C3H7I    /= 0) c(ind_C3H7I)   =  1.0E-12 * cair
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 300.E-12 * cair
      !IF (ind_C5H8    /= 0) c(ind_C5H8)    =  1.0E-09 * cair
      IF (ind_Hg       /= 0) c(ind_Hg)      = 1.68E-13 * cair

      ! examples for initializing aqueous-phase species:
      ! (the index must not be greater than APN)
      ! IF (ind_NH4p_a(2)  /=0) c(ind_NH4p_a(2))  = 300.E-12 * cair
      ! 1 nmol/L DMS:
      ! IF (ind_DMS_a(3) /=0) c(ind_DMS_a(3)) = 1E-9 * lwc(3) * N_A / 1.E3 * cair

    END SUBROUTINE x0_mbl

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_mom

      ! data from Tab. 2 of Taraborrelli et al., ACP, 9, 2751-2777 (2009):
      IF (ind_H2O2     /= 0) c(ind_H2O2)     =   7.E-09 * cair
      IF (ind_O3       /= 0) c(ind_O3)       =  30.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)       = 210.E-03 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)      = 100.E-12 * cair
      IF (ind_NO       /= 0) c(ind_NO)       =  10.E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)      = 100.E-12 * cair
      IF (ind_HONO     /= 0) c(ind_HONO)     =  40.E-14 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)     =  5.0E-12 * cair
      IF (ind_N2       /= 0) c(ind_N2)       = 780.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)      =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)     =  5.0E-09 * cair
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)    = 500.E-12 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)   =  4.0E-09 * cair
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)    = 350.E-12 * cair
      IF (ind_CO       /= 0) c(ind_CO)       = 100.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)      = 350.E-06 * cair
      IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H)  =  2.0E-09 * cair
      IF (ind_CH3CO3H  /= 0) c(ind_CH3CO3H)  =  1.5E-09 * cair
      IF (ind_ACETOL   /= 0) c(ind_ACETOL)   =  4.0E-09 * cair
      IF (ind_MGLYOX   /= 0) c(ind_MGLYOX)   = 500.E-12 * cair
      IF (ind_C5H8     /= 0) c(ind_C5H8)     =  2.0E-09 * cair
      IF (ind_PAN      /= 0) c(ind_PAN)      = 100.E-12 * cair
      ! for testing the terpene mechanism:
      IF (ind_APINENE  /= 0) c(ind_APINENE)  = 100.E-12 * cair
      IF (ind_BPINENE  /= 0) c(ind_BPINENE)  = 100.E-12 * cair
      IF (ind_CARENE   /= 0) c(ind_CARENE)   = 100.E-12 * cair
      IF (ind_SABINENE /= 0) c(ind_SABINENE) = 100.E-12 * cair
      IF (ind_CAMPHENE /= 0) c(ind_CAMPHENE) = 100.E-12 * cair
      ! replace missing species in MOZART by similar compounds:
      IF (batchfile=="mozart") THEN
        PRINT *, "adjusting initial values for mozart..."
        c(ind_HNO3)     =  40.E-14 * cair + c(ind_HNO3)    ! add HONO
        c(ind_CH3CO2H)  = 350.E-12 * cair + c(ind_CH3CO2H) ! add HCOOH
        c(ind_LTERP)    = 500.E-12 * cair ! APINENE, BPINENE, CARENE, SABINENE, CAMPHENE
      ENDIF
      ! replace missing species in CB05BASCOE by similar compounds:
      IF (batchfile=="cb05bascoe") THEN
        PRINT *, "adjusting initial values for cb05bascoe..."
        c(ind_HNO3)     =  40.E-14 * cair + c(ind_HNO3) ! add HONO
        c(ind_MACO2H)   =  3.5E-09 * cair ! CH3CO2H, CH3CO3H
        c(ind_CH3COCH3) =  4.0E-09 * cair ! ACETOL
        c(ind_LTERP)    = 500.E-12 * cair ! APINENE, BPINENE, CARENE, SABINENE, CAMPHENE
      ENDIF
      ! replace missing species in MCM and JAM by similar compounds:
      IF ((batchfile=="mcm").OR.(batchfile=="jam")) THEN
        PRINT *, "adjusting initial values for mcm and jam..."
        c(ind_APINENE)  = 100.E-12 * cair + c(ind_APINENE) ! add CARENE
        c(ind_BPINENE)  = 200.E-12 * cair + c(ind_BPINENE) ! add SABINENE, CAMPHENE
      ENDIF

    END SUBROUTINE x0_mom

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_oomph
      ! fixed species:
      IF (ind_CO2      /= 0) c(ind_CO2)     = 382.E-06 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.75E-06 * cair ! JW

      ! def:   default as in usual mbl setup
      ! Heard: North Atlantic Campaign NAMBLEX
      ! JW:    J. Williams' guess
      ! RS:    Rolf Sander's guess
      ! Sa:    Sander's lit compilation
      ! Singh: Hanwant Singh's missions 15, 16, 18 (remote clean air)
      !        in tropical Pacific = ref0314
      ! Wa:    Warneck: Nat Atm + info (page) = ref0067
      ! mod<#> different modifications to reach steady state earlier

      ! REMOTE MARINE BACKGROUND
      PRINT *, 'OOMPH init: marine background'
      IF (ind_ACETOL   /= 0) c(ind_ACETOL)  = 205.E-12 * cair ! 3rdgen
      !IF (ind_ACETOL   /= 0) c(ind_ACETOL)  = 250.E-12 * cair ! new
      IF (ind_HYPERACET    /= 0) c(ind_HYPERACET)   =  0.74E-12 * cair ! 4Igen
      !IF (ind_HYPERACET    /= 0) c(ind_HYPERACET)   =  0.7E-12 * cair ! 4gen
      !IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)  = 100.E-12 * cair ! def(off)
      IF (ind_C3H7I    /= 0) c(ind_C3H7I)   =  0.6E-12 * cair ! 4Igen
      !IF (ind_C3H7I    /= 0) c(ind_C3H7I)   =  0.7E-12 * cair ! mod5
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)  = 100.E-12 * cair ! 4gen
      IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3)= 100.E-12 * cair ! 4gen
      IF (ind_MGLYOX /= 0) c(ind_MGLYOX)= 230.E-12 * cair ! 4Igen
      !IF (ind_MGLYOX /= 0) c(ind_MGLYOX)= 280.E-12 * cair ! 3rdgen
      !IF (ind_MGLYOX /= 0) c(ind_MGLYOX)= 340.E-12 * cair ! new
      IF (ind_CH3I     /= 0) c(ind_CH3I)    =  1.24E-12 * cair ! 4Igen
      !IF (ind_CH3I     /= 0) c(ind_CH3I)    =  1.3E-12 * cair ! mod8
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 500.E-12 * cair ! 4gen
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 1.47E-09 * cair ! 3rdgen
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 1.44E-09 * cair ! mod10
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 800.E-12 * cair ! mod9+Sinha
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 300.E-12 * cair ! def+JW
      !IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 0.E-12 * cair ! 4Igentest
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 590.E-12 * cair ! 3rdgen
      !IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 530.E-12 * cair ! mod8
      !IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 0.E-09 * cair ! 4Igentest
      IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 1.34E-09 * cair ! 4Igen
      !IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 1.54E-09 * cair ! 4gen
      !IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 1.48E-09 * cair ! 3rdgen
      !IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 1.50E-09 * cair ! new
      IF (ind_CO       /= 0) c(ind_CO)      =  35.E-09 * cair ! 4Igen
      !IF (ind_CO       /= 0) c(ind_CO)      = 115.E-09 * cair ! 3rdgen, 4gen
      !IF (ind_CO       /= 0) c(ind_CO)      = 113.E-09 * cair ! mod6
      IF (ind_DMS      /= 0) c(ind_DMS)     = 200.E-12 * cair ! 4gen
      !IF (ind_DMS      /= 0) c(ind_DMS)     = 130.E-12 * cair ! 3rdgen
      !IF (ind_DMS      /= 0) c(ind_DMS)     =  30.E-12 * cair ! mod8
      IF (ind_DMSO     /= 0) c(ind_DMSO)    =  18.E-12 * cair ! 4Igen
      IF (ind_H2       /= 0) c(ind_H2)      =   1.E-06 * cair ! def
      IF (ind_H2O2     /= 0) c(ind_H2O2)    = 220.E-12 * cair ! 4Igen
      !IF (ind_H2O2     /= 0) c(ind_H2O2)    = 260.E-12 * cair ! mod6
      !IF (ind_H2O2     /= 0) c(ind_H2O2)    = 552.E-12 * cair ! Singh
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 220.E-12 * cair ! 4Igen
      !IF (ind_HCHO     /= 0) c(ind_HCHO)    = 200.E-12 * cair ! 3rdgen
      !IF (ind_HCHO     /= 0) c(ind_HCHO)    = 180.E-12 * cair ! mod3
      !IF (ind_HCHO     /= 0) c(ind_HCHO)    = 200.E-12 * cair ! Sa
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)   =   7.E-12 * cair ! 4Igen
      !IF (ind_HCOOH    /= 0) c(ind_HCOOH)   =   8.E-12 * cair ! 3rdgen
      !IF (ind_HCOOH    /= 0) c(ind_HCOOH)   =  11.E-12 * cair ! mod5
      IF (ind_HCl      /= 0) c(ind_HCl)     =  30.E-12 * cair ! 4gen
      !IF (ind_HCl      /= 0) c(ind_HCl)     =  40.E-12 * cair ! def,3gen
      !IF (ind_HCl      /= 0) c(ind_HCl)     = 100.E-12 * cair ! Sa
      IF (ind_HNO3     /= 0) c(ind_HNO3)    = 0.15E-12 * cair ! mod3
      IF (ind_ISON     /= 0) c(ind_ISON)    =  2.4E-12 * cair ! 4Igen
      !IF (ind_ISON     /= 0) c(ind_ISON)    =  3.4E-12 * cair ! 4gen
      !IF (ind_ISON     /= 0) c(ind_ISON)    =  2.2E-12 * cair ! 3rdgen
      !IF (ind_ISON     /= 0) c(ind_ISON)    =  2.8E-12 * cair ! new
      IF (ind_ISOOH    /= 0) c(ind_ISOOH)   =  30.E-12 * cair ! 4Igen
      !IF (ind_ISOOH    /= 0) c(ind_ISOOH)   =  40.E-12 * cair ! new
      !IF (ind_ISOOH    /= 0) c(ind_ISOOH)   =  40.E-12 * cair ! new
      IF (ind_C5H8     /= 0) c(ind_C5H8)    =  50.E-12 * cair ! new
      IF (ind_MVK      /= 0) c(ind_MVK)     = 130.E-12 * cair ! 4Igen
      !IF (ind_MVK      /= 0) c(ind_MVK)     = 160.E-12 * cair ! 3rdgen
      !IF (ind_MVK      /= 0) c(ind_MVK)     = 200.E-12 * cair ! new
      !IF (ind_MVKOOH   /= 0) c(ind_MVKOOH)  = 0.E-12 * cair ! 4Igentest
      IF (ind_MVKOOH   /= 0) c(ind_MVKOOH)  = 180.E-12 * cair ! 4Igen
      !IF (ind_MVKOOH   /= 0) c(ind_MVKOOH)  = 260.E-12 * cair ! 3rdgen
      !IF (ind_MVKOOH   /= 0) c(ind_MVKOOH)  = 295.E-12 * cair ! new
      IF (ind_NACA     /= 0) c(ind_NACA)    =  1.6E-12 * cair ! 4Igen
      !IF (ind_NACA     /= 0) c(ind_NACA)    =  2.1E-12 * cair ! 4gen
      !IF (ind_NACA     /= 0) c(ind_NACA)    = 1.35E-12 * cair ! 3rdgen
      !IF (ind_NACA     /= 0) c(ind_NACA)    = 1.65E-12 * cair ! new
      IF (ind_NH3      /= 0) c(ind_NH3)     = 170.E-12 * cair ! 4Igen
      !IF (ind_NH3      /= 0) c(ind_NH3)     = 200.E-12 * cair ! 4gen
      !IF (ind_NH3      /= 0) c(ind_NH3)     = 150.E-12 * cair ! 3rdgen
      !IF (ind_NH3      /= 0) c(ind_NH3)     = 250.E-12 * cair ! mod8
      !IF (ind_NO       /= 0) c(ind_NO)      =  10.E-12 * cair ! def(off)
      !IF (ind_NO       /= 0) c(ind_NO)      =  0.3E-12 * cair ! JW (<<Singh)
      !IF (ind_NO2      /= 0) c(ind_NO2)     = 1.0E-12 * cair ! JW (<<Singh)
      IF (ind_NO2      /= 0) c(ind_NO2)     =  2.E-12 * cair ! 4gen
      !IF (ind_NO2      /= 0) c(ind_NO2)     =  1.6E-12 * cair ! mod3
      IF (ind_O3       /= 0) c(ind_O3)      = 10.4E-09 * cair ! 4Igen
      !IF (ind_O3       /= 0) c(ind_O3)      = 10.0E-09 * cair ! 3rdgen
      !IF (ind_O3       /= 0) c(ind_O3)      = 10.6E-09 * cair ! mod6
      IF (ind_CH3CO3H      /= 0) c(ind_CH3CO3H)     = 300.E-12 * cair ! 4Igen
      !IF (ind_CH3CO3H      /= 0) c(ind_CH3CO3H)     = 460.E-12 * cair ! 4gen
      !IF (ind_CH3CO3H      /= 0) c(ind_CH3CO3H)     = 390.E-12 * cair ! 3rdgen
      !IF (ind_CH3CO3H      /= 0) c(ind_CH3CO3H)     = 420.E-12 * cair ! new
      IF (ind_PAN      /= 0) c(ind_PAN)     =   2.E-12 * cair ! 4gen
      !IF (ind_PAN      /= 0) c(ind_PAN)     =   1.E-12 * cair ! new
      IF (ind_SO2      /= 0) c(ind_SO2)     = 135.E-12 * cair ! 4Igen
      !IF (ind_SO2      /= 0) c(ind_SO2)     = 130.E-12 * cair ! 3rdgen
      !IF (ind_SO2      /= 0) c(ind_SO2)     =  35.E-12 * cair ! mod8

      ! NOT TAKEN INTO ACCOUNT EVEN THOUGH PRESENT:
      !IF (ind_C2H2     /= 0) c(ind_C2H2)     = 28.7E-12 * cair ! Singh
      !IF (ind_C2H4     /= 0) c(ind_C2H4)     =  21.E-12 * cair ! Singh
      !IF (ind_C2H6     /= 0) c(ind_C2H6)     = 28.7E-12 * cair ! Singh
      !IF (ind_C3H6     /= 0) c(ind_C3H6)     = 11.4E-12 * cair ! Singh
      !IF (ind_C3H8     /= 0) c(ind_C3H8)     =  16.E-12 * cair ! Singh
      !IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3) = 300.E-12 * cair ! JW
      !IF (ind_PAN      /= 0) c(ind_PAN)      =   2.E-12 * cair ! Singh
    END SUBROUTINE x0_oomph

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_iso
#ifdef MECCA_TAG
      use messy_mecca_tag_common

      real(dp) :: R18O, R17O, d18O, d17O

#include "messy_mecca_tag_parameters.inc"
#include "./mecca/tag/caaba_exp.inc"

      IF (ind_CO2      /= 0) c(ind_CO2)     = 382.E-06 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.80E-06 * cair

      ! Photochem.smog parcel chemically degrading in the low troposphere
      ! + oxygen isotopes

#ifdef EXP_PM_TANS
      PRINT *, 'ISO init: PM:PI+Tans experiment test'
     !IF (ind_CH4      /= 0) c(ind_CH4)     = 1.80E-06 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     = 0.0

      IF (ind_OH       /= 0) c(ind_OH)      = 1.e6
      IF (ind_O1D      /= 0) c(ind_O1D)     = 1.e0
      IF (ind_Cl       /= 0) c(ind_Cl)      = 1.e2
      ! "Tans" tracers
      IF (ind_CH4t     /= 0) c(ind_CH4t    )     = c(ind_CH4)
      IF (ind_CH4c     /= 0) c(ind_CH4c    )     = c(ind_CH4)
      IF (ind_CH4c0    /= 0) c(ind_CH4c0   )     = c(ind_CH4)
      IF (ind_CH4m     /= 0) c(ind_CH4m    )     = c(ind_CH4)
      ! "Spivakovsky" OH
      IF (ind_OHc      /= 0) c(ind_OHc     )     = c(ind_OH)
      IF (ind_OHc0     /= 0) c(ind_OHc0    )     = c(ind_OH)/2_dp
      ! modulated OH
      IF (ind_OHm      /= 0) c(ind_OHm     )     = c(ind_OH)*2_dp
      ! eff. enrichment
      IF (ind_CH4ee    /= 0) c(ind_CH4     )     = c(ind_CH4)
#endif

#ifdef ISO_TROP
      PRINT *, 'ISO init: tropospheric composition'
      ! Photochem. smog parcel chemically degrading in the low troposphere
      IF (ind_ACETOL   /= 0) c(ind_ACETOL)    = 250.E-12 * cair
      IF (ind_HYPERACET/= 0) c(ind_HYPERACET) = 0.74E-12 * cair ! 4Igen
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)    = 100.E-12 * cair ! 4gen
      IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3)  = 100.E-12 * cair ! 4gen
      IF (ind_MGLYOX   /= 0) c(ind_MGLYOX)    = 230.E-12 * cair ! 4Igen
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)     = 500.E-12 * cair ! 4gen
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)    = 590.E-12 * cair ! 3rdgen
      IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H)   = 1.34E-09 * cair ! 4Igen
      IF (ind_CO       /= 0) c(ind_CO)        =  80.E-09 * cair ! 4Igen  ! 200
      IF (ind_DMS      /= 0) c(ind_DMS)       = 200.E-12 * cair ! 4gen
      IF (ind_DMSO     /= 0) c(ind_DMSO)      =  18.E-12 * cair ! 4Igen
      IF (ind_H2       /= 0) c(ind_H2)        =   1.E-06 * cair ! def
      IF (ind_H2O2     /= 0) c(ind_H2O2)      = 220.E-12 * cair ! 4Igen
      IF (ind_HCHO     /= 0) c(ind_HCHO)      = 220.E-12 * cair ! 4Igen
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)     =   7.E-12 * cair ! 4Igen
      IF (ind_HNO3     /= 0) c(ind_HNO3)      = 0.15E-12 * cair ! mod3
      IF (ind_NH3      /= 0) c(ind_NH3)       = 170.E-12 * cair ! 4Igen
      IF (ind_NO2      /= 0) c(ind_NO2)       =  10.E-12 * cair ! 4gen 2
      IF (ind_O3       /= 0) c(ind_O3)        =  25.E-09 * cair ! 4Igen 10.2
      IF (ind_CH3CO3H  /= 0) c(ind_CH3CO3H)   = 300.E-12 * cair ! 4Igen
      IF (ind_PAN      /= 0) c(ind_PAN)       =   2.E-12 * cair ! 4gen
      IF (ind_SO2      /= 0) c(ind_SO2)       = 135.E-12 * cair ! 4Igen
#endif

#if defined(ISO_STRAT)||defined(EXP_STOZ)
      PRINT *, 'ISO init: stratospheric parcel'
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)    =   1.E-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   =   1.E-10 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair

      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.60E-06 * cair
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)  = 500.E-12 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  30.E-09 * cair
      IF (ind_H2       /= 0) c(ind_H2)     = 550.E-09 * cair
    ! IF (ind_H2O2     /= 0) c(ind_H2O2)   = 220.E-12 * cair
    ! IF (ind_HCHO     /= 0) c(ind_HCHO)   = 220.E-12 * cair
    ! IF (ind_HCOOH    /= 0) c(ind_HCOOH)  =   7.E-12 * cair
    ! IF (ind_NH3      /= 0) c(ind_NH3)    = 170.E-12 * cair
    ! IF (ind_NO2      /= 0) c(ind_NO2)    =   2.E-12 * cair
      IF (ind_O3       /= 0) c(ind_O3)     = 500.E-09 * cair
      IF (ind_Cl       /= 0) c(ind_Cl)     =   1.E-21 * cair
    ! IF (ind_HCl      /= 0) c(ind_HCl)    = 400.E-12 * cair
#ifdef EXP_STOZ
      IF (ind_O3       /= 0) c(ind_O3)     = 0._dp !100.E-09 * cair
#endif

#ifdef HTRANSTEST
      ! H transfer test
      IF (ind_A4       /= 0) c(ind_A4)     = 1.6E-09 * cair
      IF (ind_Z3       /= 0) c(ind_Z3)     = 1.6E-09 * cair
#endif
#endif

#ifdef CIO_YEUNG_EXP
      C(:) = 0._dp
      C(ind_CIOIEXBOOST) = 1._dp

    ! we simulate Yeung et al. [JGR 2014] exp. lab conditions
      cair = 1.e17  ! cm^-3

    ! init Yeung's mech isotopologues
      IF (ind_iOO /= 0) c(ind_iOO)    = cair

      d18O = -35._dp !-36e-3_dp  ! -36 permil exp.
      d17O = 0._dp !
      R18O = (1.+d18O/1e3)*VSMOW_18O
      R17O = (1.+d17O/1e3)*VSMOW_17O

      IF (ind_iOP /= 0) C(ind_iOP)    = C(ind_iOO) * isofrac3r( 0._dp, R17O, 0._dp, R18O, 2) - R17O**2
      IF (ind_iOQ /= 0) C(ind_iOQ)    = C(ind_iOO) * isofrac3r( 0._dp, R18O, 0._dp, R17O, 2) - R18O**2
      IF (ind_iOO /= 0) C(ind_iOO)    = C(ind_iOO) - ( C(ind_iOQ)+C(ind_iOP) ) / 2.

      IF (ind_iPP /= 0) C(ind_iPP)    = C(ind_iOO) * R17O**2
      IF (ind_iQQ /= 0) C(ind_iQQ)    = C(ind_iOO) * R18O**2
      IF (ind_iPQ /= 0) C(ind_iPQ)    = C(ind_iOO) * 2.*R17O*R18O

    ! init MECCA counterparts
      IF (ind_O2  /= 0) C(ind_O2)     = cair
#endif

#ifdef CIO_ATOMTR_TEST
      if (ind_O3P       /= 0) C(ind_O3P)    = 1e-12_dp*cair
      if (ind_CIOIEXBOOST/=0) C(ind_CIOIEXBOOST) = 1._dp
#endif

#endif
    END SUBROUTINE x0_iso

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_tag
#ifdef MECCA_TAG
      IF (ind_CO2      /= 0) c(ind_CO2)     = 382.E-06 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.60E-06 * cair

    ! Photochem. parcel chemically degrading in the stratosphere
      PRINT *, 'TAG init: scenarios for kinetic tagging'
#endif
    END SUBROUTINE x0_tag

    ! ------------------------------------------------------------------------

    ! mz_ab_20091111+
    SUBROUTINE x0_strato20

      ! stratosphere, 20 hPa
      ! from scout02 (ProSECCO simulation)

      IF (ind_H        /= 0) c(ind_H)      =   1.E-12 * cair
      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)    =   1.E-15 * cair
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)    =   1.E-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   =   1.E-10 * cair
      IF (ind_CL       /= 0) c(ind_CL)     =   1.E-21 * cair
      IF (ind_CLO      /= 0) c(ind_CLO)    =   1.E-15 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   =  40.E-12 * cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  =   1.E-13 * cair
      IF (ind_CL2      /= 0) c(ind_CL2)    =   1.E-13 * cair
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  =   1.E-12 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)    =  1.3E-07 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  1.4E-08 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) =  12.E-12 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  =   6.E-10 * cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  =  1.4E-11 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) =   1.E-12 * cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  =   1.E-12 * cair
      IF (ind_CCL4     /= 0) c(ind_CCL4)   =   1.E-12 * cair
      IF (ind_CH3CCL3  /= 0) c(ind_CH3CCL3)=   1.E-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)    =   1.E-12 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
      IF (ind_H2       /= 0) c(ind_H2)     =   5.E-07 * cair
      IF (ind_O3       /= 0) c(ind_O3)     =   4.E-06 * cair
      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   =   7.E-11 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 350.E-06 * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 450.E-12 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 400.E-12 * cair
      ! additional for mecca
      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair

    END SUBROUTINE x0_strato20
    ! mz_ab_20091111-

    ! ------------------------------------------------------------------------

    ! mz_ab_20091111+
    SUBROUTINE x0_strato10

      ! stratosphere, 10 hPa
      IF (ind_H        /= 0) c(ind_H)      =   1.E-16 * cair
      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)    =   1.E-15 * cair
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)    =   1.E-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   =   1.5E-10 * cair
      IF (ind_CL       /= 0) c(ind_CL)     =   1.E-21 * cair
      IF (ind_CLO      /= 0) c(ind_CLO)    =   1.E-15 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   =  40.E-12 * cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  =   1.E-13 * cair
      IF (ind_CL2      /= 0) c(ind_CL2)    =   1.E-13 * cair
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  =   1.E-12 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)    =  1.3E-07 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  1.4E-08 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) =  12.E-12 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  =   6.E-10 * cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  =  1.4E-11 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) =   1.E-12 * cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  =   1.E-12 * cair
      IF (ind_CCL4     /= 0) c(ind_CCL4)   =   1.E-12 * cair
      IF (ind_CH3CCL3  /= 0) c(ind_CH3CCL3)=   1.E-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)    = 4.251E-06 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
      IF (ind_H2       /= 0) c(ind_H2)     =   5.E-07 * cair
      IF (ind_O3       /= 0) c(ind_O3)     =   8.E-06 * cair
      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   =   7.E-11 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 350.E-06 * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 180.E-12 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 400.E-12 * cair
      ! additional for mecca
      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair

    END SUBROUTINE x0_strato10
    ! mz_ab_20091111-

    ! ------------------------------------------------------------------------

    ! mz_ab_20091111+
    SUBROUTINE x0_mtchem

!!$      ! mesosphere, 0.01 hPa
!!$      IF (ind_H        /= 0) c(ind_H)      = 1.876E-07* cair
!!$      IF (ind_OH       /= 0) c(ind_OH)     = 1.240E-08* cair
!!$      IF (ind_HO2      /= 0) c(ind_HO2)    = 5.012E-09* cair
!!$      IF (ind_N        /= 0) c(ind_N)      = 2.114E-10* cair
!!$      IF (ind_NO3      /= 0) c(ind_NO3)    = 3.083E-21* cair
!!$      IF (ind_N2O5     /= 0) c(ind_N2O5)   = 9.072E-27* cair
!!$      IF (ind_HNO4     /= 0) c(ind_HNO4)   = 4.754E-18* cair
!!$      IF (ind_CL       /= 0) c(ind_CL)     = 8.001E-11* cair
!!$      IF (ind_CLO      /= 0) c(ind_CLO)    = 1.564E-13* cair
!!$      IF (ind_HOCl     /= 0) c(ind_HOCl)   = 7.015E-15* cair
!!$      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  = 1.558E-25* cair
!!$      IF (ind_CL2      /= 0) c(ind_CL2)    = 1.E-13* cair ! ????
!!$      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  = 1.E-12* cair! ????
!!$      IF (ind_N2O      /= 0) c(ind_N2O)    = 7.077E-11* cair
!!$      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 1.685E-12* cair
!!$      IF (ind_CO       /= 0) c(ind_CO)     = 9.862E-07* cair
!!$      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) = 3.558E-13* cair
!!$      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  = 1.460E-22* cair
!!$      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  = 2.854E-38* cair
!!$      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) = 2.776E-17* cair
!!$      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  = 2.212E-14* cair
!!$      IF (ind_HNO3     /= 0) c(ind_HNO3)   = 2.084E-13* cair
!!$      IF (ind_H2O      /= 0) c(ind_H2O)    = 4.567E-06* cair
!!$      IF (ind_O3P      /= 0) c(ind_O3P)    = 5.350E-06* cair
!!$      IF (ind_O1D      /= 0) c(ind_O1D)    = 3.143E-14* cair
!!$      IF (ind_H2       /= 0) c(ind_H2)     = 1.310E-06* cair
!!$      IF (ind_O3       /= 0) c(ind_O3)     = 7.823E-08* cair
!!$      IF (ind_NO       /= 0) c(ind_NO)     = 2.454E-09* cair
!!$      IF (ind_NO2      /= 0) c(ind_NO2)    = 1.685E-12* cair
!!$      IF (ind_CH4      /= 0) c(ind_CH4)    = 1.113E-07* cair
!!$      IF (ind_HCHO     /= 0) c(ind_HCHO)   = 1.417E-12* cair
!!$      IF (ind_CO       /= 0) c(ind_CO)     = 9.862E-07* cair
!!$      IF (ind_CO2      /= 0) c(ind_CO2)    = 3.641E-04* cair
!!$      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 1.169E-10* cair
!!$      IF (ind_HCl      /= 0) c(ind_HCl)    = 3.342E-09* cair
!!$      ! additional for mecca
!!$      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
!!$      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
!!$      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair
!!$      IF (ind_em       /= 0) c(ind_em)     = 1.2E-14 * cair
!!$      IF (ind_NOp      /= 0) c(ind_NOp)     = 1.2E-10 * cair
!!$      IF (ind_O2p      /= 0) c(ind_O2p)     = 1.2E-10 * cair
!!$      IF (ind_Op      /= 0) c(ind_Op)     = 1.2E-10 * cair
!!$      IF (ind_Np      /= 0) c(ind_Np)     = 1.2E-10 * cair
!!$      IF (ind_H2O      /= 0) c(ind_H2O)     = 1.2E-7 * cair

      ! from Miriam Sinnhuber, 0.51 Pa
      IF (ind_H        /= 0) c(ind_H)      = 7.69140E-07* cair
      IF (ind_OH       /= 0) c(ind_OH)     = 6.53740E-10* cair
      IF (ind_HO2      /= 0) c(ind_HO2)    = 2.09240E-10* cair
      IF (ind_N        /= 0) c(ind_N)      = 2.114E-10* cair
      IF (ind_NO3      /= 0) c(ind_NO3)    = 2.07260E-24* cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   = 8.59140E-30* cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   = 4.754E-18* cair
      IF (ind_CL       /= 0) c(ind_CL)     = 8.09210E-10* cair
      IF (ind_CLO      /= 0) c(ind_CLO)    = 1.28810E-13* cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   = 8.07360E-18* cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  = 4.67980E-27* cair
      IF (ind_CL2      /= 0) c(ind_CL2)    = 1.E-13* cair ! ????
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  = 1.E-12* cair! ????
      IF (ind_N2O      /= 0) c(ind_N2O)    = 7.84100E-10* cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 4.64940E-13* cair
      IF (ind_CO       /= 0) c(ind_CO)     = 3.11580E-05* cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) = 3.558E-13* cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  = 3.84100E-26* cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  = 2.854E-38* cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) = 2.776E-17* cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  = 2.212E-14* cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   = 1.14160E-18* cair
      IF (ind_O3P      /= 0) c(ind_O3P)    = 1.18100E-04* cair
      IF (ind_O1D      /= 0) c(ind_O1D)    = 9.87100E-14* cair
      IF (ind_H2       /= 0) c(ind_H2)     = 5.29290E-06* cair
      IF (ind_O3       /= 0) c(ind_O3)     = 8.45990E-08* cair
      IF (ind_NO       /= 0) c(ind_NO)     = 2.05670E-08* cair
      IF (ind_NO2      /= 0) c(ind_NO2)    = 4.20810E-14* cair
      IF (ind_CH4      /= 0) c(ind_CH4)    = 2.63900E-09* cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   = 1.417E-12* cair
      IF (ind_CO       /= 0) c(ind_CO)     = 3.11580E-05 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 3.54360E-04* cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 1.169E-10* cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 2.39140E-09* cair
      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair
!!$      IF (ind_em       /= 0) c(ind_em)     = 1.2E-14 * cair
!!$      IF (ind_NOp      /= 0) c(ind_NOp)     = 1.2E-10 * cair
!!$      IF (ind_O2p      /= 0) c(ind_O2p)     = 1.2E-10 * cair
!!$      IF (ind_Op      /= 0) c(ind_Op)     = 1.2E-10 * cair
!!$      IF (ind_Np      /= 0) c(ind_Np)     = 1.2E-10 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)     = 3.36210E-09 * cair
      !mz_hr_20160515+ '* cair' was missing
      !IF (ind_N       /= 0) c(ind_N)     =9.27150E-10
      IF (ind_N       /= 0) c(ind_N)     =9.27150E-10 * cair
      !mz_hr_20160515-

    END SUBROUTINE x0_mtchem
    ! mz_ab_20091111-

    ! ------------------------------------------------------------------------

    ! op_ff_20160309+
    SUBROUTINE x0_tropopause
      ! Reference concentrations of air parcel at:
      ! lat      = 0
      ! lon      = 0
      ! month    = 3
      ! pressure = 226.32hPa

      ! values were taken from a climatology
      ! of the ESCiMo simulation RC1SD-base-07

      ! tr_transp
      IF (ind_CO2      /= 0) c(ind_CO2)     = 0.00036392  * cair

      ! tr_NOx_NOy
      IF (ind_N        /= 0) c(ind_N)       = 0 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)     = 4.18354e-14 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)    = 2.76229e-12 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)    = 3.3254e-11  * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)    = 6.57861e-11 * cair
      IF (ind_NO       /= 0) c(ind_NO)      = 3.56444e-11 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)     = 3.12893e-07 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 0.78        * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     = 5.75031e-11 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)     = 1.82963e-11 * cair

      ! tr_Ox_HOx
      IF (ind_H        /= 0) c(ind_H)       = 4.9586e-20  * cair
      IF (ind_OH       /= 0) c(ind_OH)      = 2.10936e-13 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)     = 6.11271e-12 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)     = 0.00017312  * cair
      IF (ind_O3P      /= 0) c(ind_O3P)     = 1.28702e-15 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)     = 5.63444e-21 * cair
      IF (ind_H2       /= 0) c(ind_H2)      = 5.54021e-07 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 0.21        * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)    = 4.57365e-10 * cair
      IF (ind_O3       /= 0) c(ind_O3)      = 5.90356e-08 * cair

      ! tr_hycarbs
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.71916e-06 * cair

      ! tr_alks
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 9.80751e-10 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 1.00918e-07 * cair
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)   = 4.96023e-11 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 9.80851e-11 * cair

      ! tr_chlorine
      IF (ind_Cl       /= 0) c(ind_Cl)      = 3.21085e-17 * cair
      IF (ind_ClO      /= 0) c(ind_ClO)     = 9.72203e-15 * cair
      IF (ind_OClO     /= 0) c(ind_OClO)    = 7.33498e-18 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)    = 5.01393e-14 * cair
      IF (ind_Cl2O2    /= 0) c(ind_Cl2O2)   = 7.65851e-20 * cair
      IF (ind_Cl2      /= 0) c(ind_Cl2)     = 2.25151e-19 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)   = 1.0305e-13  * cair
      IF (ind_HCl      /= 0) c(ind_HCl)     = 6.78836e-12 * cair

      ! tr_perox
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)   = 9.99738e-13 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 1.37036e-10 * cair

      ! tr_halocarbs
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)   = 2.65308e-10 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2)  = 5.32336e-10 * cair
      IF (ind_CH3Cl    /= 0) c(ind_CH3Cl)   = 5.01803e-10 * cair
      IF (ind_CCl4     /= 0) c(ind_CCl4)    = 1.01461e-10 * cair
      IF (ind_CH3CCl3  /= 0) c(ind_CH3CCl3) = 7.31265e-11 * cair

      ! tr_sulphur
      IF (ind_SO2      /= 0) c(ind_SO2)     = 7.54215e-12 * cair
    END SUBROUTINE x0_tropopause

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_low_strato
      ! Reference concentrations of air parcel at:
      ! lat      = 0
      ! lon      = 0
      ! month    = 3
      ! pressure = 120.45hPa

      ! values were taken from a climatology
      ! of the ESCiMo simulation RC1SD-base-07

      ! tr_transp
      IF (ind_CO2      /= 0) c(ind_CO2)     = 0.000363703 * cair

      ! tr_NOx_NOy
      IF (ind_N        /= 0) c(ind_N)       = 4.3354e-21  * cair
      IF (ind_NO3      /= 0) c(ind_NO3)     = 8.80458e-15 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)    = 2.32101e-12 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)    = 7.64665e-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)    = 9.31376e-11 * cair
      IF (ind_NO       /= 0) c(ind_NO)      = 2.0396e-10  * cair
      IF (ind_N2O      /= 0) c(ind_N2O)     = 3.12665e-07 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 0.78        * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     = 1.97631e-10 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)     = 1.71703e-11 * cair

      ! tr_Ox_HOx
      IF (ind_H        /= 0) c(ind_H)       = 1.13608e-19 * cair
      IF (ind_OH       /= 0) c(ind_OH)      = 4.71945e-13 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)     = 2.11008e-12 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)     = 4.89619e-06 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)     = 3.47866e-15 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)     = 1.16954e-20 * cair
      IF (ind_H2       /= 0) c(ind_H2)      = 5.53817e-07 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 0.21        * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)    = 1.65937e-10 * cair
      IF (ind_O3       /= 0) c(ind_O3)      = 7.73848e-08 * cair

      ! tr_hycarbs
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.71615e-06 * cair

      ! tr_alks
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 8.33597e-10 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 8.19542e-08 * cair
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)   = 4.17248e-11 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 5.36401e-11 * cair

      ! tr_chlorine
      IF (ind_Cl       /= 0) c(ind_Cl)      = 1.63985e-16 * cair
      IF (ind_ClO      /= 0) c(ind_ClO)     = 1.16066e-14 * cair
      IF (ind_OClO     /= 0) c(ind_OClO)    = 4.01014e-17 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)    = 1.65232e-14 * cair
      IF (ind_Cl2O2    /= 0) c(ind_Cl2O2)   = 3.25006e-19 * cair
      IF (ind_Cl2      /= 0) c(ind_Cl2)     = 7.79006e-19 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)   = 8.56242e-14 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)     = 1.21202e-11 * cair

      ! tr_perox
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)   = 1.86593e-13 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 3.77783e-11 * cair

      ! tr_halocarbs
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)   = 2.64856e-10 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2)  = 5.31734e-10 * cair
      IF (ind_CH3Cl    /= 0) c(ind_CH3Cl)   = 4.99562e-10 * cair
      IF (ind_CCl4     /= 0) c(ind_CCl4)    = 1.01214e-10 * cair
      IF (ind_CH3CCl3  /= 0) c(ind_CH3CCl3) = 7.30846e-11 * cair

      ! tr_sulphur
      IF (ind_SO2      /= 0) c(ind_SO2)     = 6.30675e-12 * cair
    END SUBROUTINE x0_low_strato

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_mid_strato
      ! Reference concentrations of air parcel at:
      ! lat      = 0
      ! lon      = 0
      ! month    = 3
      ! pressure = 25.11hPa

      ! values were taken from a climatology
      ! of the ESCiMo simulation RC1SD-base-07

      ! tr_transp
      IF (ind_CO2      /= 0) c(ind_CO2)     = 0.000359093 * cair

      ! tr_NOx_NOy
      IF (ind_N        /= 0) c(ind_N)       = 3.88566e-18 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)     = 3.95911e-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)    = 2.85476e-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)    = 1.3086e-10  * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)    = 2.7466e-09  * cair
      IF (ind_NO       /= 0) c(ind_NO)      = 3.76835e-10 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)     = 2.53761e-07 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 0.78        * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     = 1.10275e-09 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)     = 4.13072e-12 * cair

      ! tr_Ox_HOx
      IF (ind_H        /= 0) c(ind_H)       = 7.28588e-19 * cair
      IF (ind_OH       /= 0) c(ind_OH)      = 1.51964e-12 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)     = 1.11571e-11 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)     = 3.83427e-06 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)     = 9.46588e-12 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)     = 6.04968e-18 * cair
      IF (ind_H2       /= 0) c(ind_H2)      = 5.71968e-07 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 0.21        * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)    = 6.02973e-11 * cair
      IF (ind_O3       /= 0) c(ind_O3)      = 5.70224e-06 * cair

      ! tr_hycarbs
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.47326e-06 * cair

      ! tr_alks
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 1.90836e-13 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 1.55145e-08 * cair
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)   = 2.52686e-16 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 3.84811e-11 * cair

      ! tr_chlorine
      IF (ind_Cl       /= 0) c(ind_Cl)      = 1.89024e-14 * cair
      IF (ind_ClO      /= 0) c(ind_ClO)     = 4.78724e-11 * cair
      IF (ind_OClO     /= 0) c(ind_OClO)    = 8.43409e-13 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)    = 3.85529e-11 * cair
      IF (ind_Cl2O2    /= 0) c(ind_Cl2O2)   = 8.24974e-14 * cair
      IF (ind_Cl2      /= 0) c(ind_Cl2)     = 4.56747e-14 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)   = 4.25805e-10 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)     = 1.02317e-09 * cair

      ! tr_perox
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)   = 4.62244e-13 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 5.98058e-12 * cair

      ! tr_halocarbs
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)   = 8.22788e-11 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2)  = 4.0098e-10  * cair
      IF (ind_CH3Cl    /= 0) c(ind_CH3Cl)   = 3.17313e-10 * cair
      IF (ind_CCl4     /= 0) c(ind_CCl4)    = 1.78299e-11 * cair
      IF (ind_CH3CCl3  /= 0) c(ind_CH3CCl3) = 1.66882e-11 * cair

      ! tr_sulphur
      IF (ind_SO2      /= 0) c(ind_SO2)     = 1.1427e-18  * cair
    END SUBROUTINE x0_mid_strato

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_high_strato
      ! Reference concentrations of air parcel at:
      ! lat      = 0
      ! lon      = 0
      ! month    = 3
      ! pressure = 2.77hPa

      ! values were taken from a climatology
      ! of the ESCiMo simulation RC1SD-base-07

      ! tr_transp
      IF (ind_CO2      /= 0) c(ind_CO2)     = 0.000356097 * cair

      ! tr_NOx_NOy
      IF (ind_N        /= 0) c(ind_N)       = 1.02777e-14 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)     = 1.15451e-10 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)    = 5.13655e-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)    = 7.00639e-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)    = 1.63597e-10 * cair
      IF (ind_NO       /= 0) c(ind_NO)      = 7.71614e-09 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)     = 4.647e-08   * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 0.78        * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     = 8.19046e-09 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)     = 8.06315e-15 * cair

      ! tr_Ox_HOx
      IF (ind_H        /= 0) c(ind_H)       = 7.76142e-14 * cair
      IF (ind_OH       /= 0) c(ind_OH)      = 1.22653e-10 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)     = 8.22616e-11 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)     = 5.45918e-06 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)     = 6.9964e-09  * cair
      IF (ind_O1D      /= 0) c(ind_O1D)     = 1.90289e-15 * cair
      IF (ind_H2       /= 0) c(ind_H2)      = 5.38355e-07 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 0.21        * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)    = 4.50636e-11 * cair
      IF (ind_O3       /= 0) c(ind_O3)      = 5.6866e-06  * cair

      ! tr_hycarbs
      IF (ind_CH4      /= 0) c(ind_CH4)     = 7.71146e-07 * cair

      ! tr_alks
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 1.20392e-13 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 5.72401e-08 * cair
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)   = 1.84384e-26 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 2.0852e-10  * cair

      ! tr_chlorine
      IF (ind_Cl       /= 0) c(ind_Cl)      = 2.44578e-12 * cair
      IF (ind_ClO      /= 0) c(ind_ClO)     = 2.62561e-10 * cair
      IF (ind_OClO     /= 0) c(ind_OClO)    = 8.50329e-12 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)    = 8.86478e-11 * cair
      IF (ind_Cl2O2    /= 0) c(ind_Cl2O2)   = 2.83879e-15 * cair
      IF (ind_Cl2      /= 0) c(ind_Cl2)     = 9.38255e-14 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)   = 6.58272e-11 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)     = 2.50195e-09 * cair

      ! tr_perox
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)   = 4.20119e-12 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 3.5432e-12  * cair

      ! tr_halocarbs
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)   = 4.94957e-18 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2)  = 2.3584e-11  * cair
      IF (ind_CH3Cl    /= 0) c(ind_CH3Cl)   = 3.38793e-11 * cair
      IF (ind_CCl4     /= 0) c(ind_CCl4)    = 7.8721e-24  * cair
      IF (ind_CH3CCl3  /= 0) c(ind_CH3CCl3) = 2.73592e-21 * cair

      ! tr_sulphur
      IF (ind_SO2      /= 0) c(ind_SO2)     = 3.57895e-23 * cair
    END SUBROUTINE x0_high_strato
    ! op_ff_20160309-

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_volcano
      IF (ind_CH4      /= 0) c(ind_CH4)     = 0.
      IF (ind_O3       /= 0) c(ind_O3)      = 0.
      IF (ind_CO       /= 0) c(ind_CO)      = 0.
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 700.E-06 * cair
      IF (ind_SO2      /= 0) c(ind_SO2)     =  40.E-06 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)     =  20.E-06 * cair
      IF (ind_HBr      /= 0) c(ind_HBr)     =   2.E-06 * cair
      IF (ind_HI       /= 0) c(ind_HI)      =  0.1E-06 * cair
    END SUBROUTINE x0_volcano

    ! ------------------------------------------------------------------------

  END SUBROUTINE x0

  !***************************************************************************

  SUBROUTINE mecca_physc

    USE messy_mecca_kpp,         ONLY: REQ_MCFCT
    USE messy_mecca,             ONLY: steady_state_reached
    USE messy_mecca_aero,        ONLY: mecca_aero_calc_k_ex, &
                                       mecca_aero_calc_k_ex_ocean
    USE messy_main_tools,        ONLY: spechum2mr, mr2spechum, spec2relhum,&
                                       cair_c
!!#D radjimt +
    USE messy_radjimt,           ONLY: jx_radjimt => jx
!!#D radjimt -
    USE messy_readj,             ONLY: jx_readj   => jx
    USE messy_sappho,            ONLY: jx_sappho  => jx
    USE messy_dissoc,            ONLY: dissoc_2d
    USE messy_jval,              ONLY: jval_2d
    USE caaba_mem,               ONLY: l_skipkpp, photol_clev, zmix, &
                                       l_ignore_relhum

    LOGICAL :: l_het(APN) =.TRUE.
    REAL(DP), DIMENSION(NBL,NSPEC) :: cbl
    INTEGER :: ip, status

    IF (l_ff) THEN
      ! at day 4 switch frost flowers on:
      IF (model_time >= (model_start_day+4.)*OneDay) THEN
        xaer(2) = 1.
      ENDIF
    ENDIF

    CALL fill_temp(status, SPREAD(temp,1,NBL))
    CALL fill_cair(status, SPREAD(cair,1,NBL))
    CALL fill_press(status, SPREAD(press,1,NBL)) ! mz_pj_20080716
    IF (REQ_MCFCT) THEN
      CALL fill_mcexp(status, SPREAD(mcexp,1,NBL))
    ENDIF

    ! The following values are only needed for the mesosphere. Since
    ! CAABA has only one type of temperature, it is also used for
    ! temp_elec and temp_ion:
    CALL fill_temp_elec(status, SPREAD(temp,1,NBL))
    CALL fill_temp_ion(status, SPREAD(temp,1,NBL))

    ! dummy values:
    dummy_khet_St(:) = 0.
    dummy_khet_Tr(:) = 0.
    CALL fill_khet_Tr(status, SPREAD(dummy_khet_Tr,1,NBL))
    CALL fill_khet_St(status, SPREAD(dummy_khet_St,1,NBL))

    SELECT CASE (TRIM(photrat_channel))
    CASE('dissoc')
      DO ip=1, IP_MAX
        jx(ip) = dissoc_2d(ip)%ptr(1,photol_clev)
      ENDDO
    CASE('jval')
      DO ip=1, IP_MAX
        jx(ip) = jval_2d(ip)%ptr(1,photol_clev)
      ENDDO
!!#D radjimt +
    CASE('radjimt')
      DO ip=1, IP_MAX
        jx(ip) = jx_radjimt(photol_clev,1,ip,1)
      ENDDO
!!#D radjimt -
    CASE('readj')
      DO ip=1, IP_MAX
        jx(ip) = jx_readj(ip)
      ENDDO
    CASE('sappho')
      DO ip=1, IP_MAX
        jx(ip) = jx_sappho(ip)
      ENDDO
    ENDSELECT

    ! transfer of jx in mecca to jx in kpp
    CALL fill_jx(status, SPREAD(jx,1,NBL))

    IF (l_aero) THEN
      CALL aerosol_exchng ! exchange with fresh aerosol
      ! first, calculate exchange coefficients for aerosols:
      CALL mecca_aero_calc_k_ex( &
        radius, temp, press, l_het, xaer, lwc, c, &         ! in
        k_exf, k_exb, k_exf_N2O5, k_exf_ClNO3, k_exf_BrNO3) ! out
      ! next, calculate exchange coefficients for ocean surface:
      CALL mecca_aero_calc_k_ex_ocean(xaer, radius, temp, zmix, k_exf, k_exb)
      CALL fill_lwc(status, SPREAD(lwc,1,NBL))
      CALL fill_cvfac(status, SPREAD(cvfac,1,NBL))
      CALL fill_xaer(status, SPREAD(xaer,1,NBL))
      CALL fill_k_exf(status, SPREAD(k_exf,1,NBL))
      CALL fill_k_exb(status, SPREAD(k_exb,1,NBL))
      CALL fill_k_exf_N2O5(status, SPREAD(k_exf_N2O5,1,NBL))
      CALL fill_k_exf_ClNO3(status, SPREAD(k_exf_ClNO3,1,NBL))
      CALL fill_k_exf_BrNO3(status, SPREAD(k_exf_BrNO3,1,NBL))
    ENDIF

#ifdef MECCA_TAG
    IF (l_tag) CALL mecca_tag_preprocess
#endif

    CALL check_range('before kpp:',c(:))
    IF (.NOT.l_skipkpp) THEN
      c(:) = MAX(c(:),0._DP) ! set negative values to zero
      cbl = SPREAD(c,1,NBL) ! add one dummy dimension
      CALL kpp_integrate(timesteplen,cbl,ierrf=ierrf,istatus=istatus) ! main kpp call
      c = cbl(1,:)          ! remove the dummy dimension
      CALL check_range('after kpp: ',c(:))
    ENDIF

#ifdef MECCA_TAG
    IF (l_tag) CALL mecca_tag_postprocess
#endif

    ! update cair, spechum, relhum, c(H2O)
    IF (l_ignore_relhum) THEN
      ! update cair with current c(H2O) after kpp_integrate
      cair = cair_c(c(ind_H2O), temp, press, l_hum_emac)

      ! update concentrations with new cair
      C(:) = C(:) * cair/cair_old
      cair_old = cair

      ! update spechum and relhum with new c(H2O) and cair
      spechum = mr2spechum(status, c(ind_H2O)/cair)
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR: mecca_physc: mr2spechum'
        STOP 1
      ENDIF

      relhum  = spec2relhum(status, spechum, temp, press, &
                            l_hum_emac, l_psat_liquid, l_relhum_wmo)
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR: mecca_physc: spec2relhum'
        STOP 1
      ENDIF
    ELSE
      ! update c(H2O) with current spechum, overwriting KPP value
      !c(ind_H2O) = cair * &
      !  rh2mr(status,relhum,temp,press,l_hum_emac,l_relhum_wmo)
      c(ind_H2O) = cair * spechum2mr(status,spechum)
      IF (status > 1) THEN
        WRITE(*,*) "ERROR in mecca_physc: spechum2mr"
        STOP 1
      ENDIF
      !IF (status /= 0) STOP 1
    ENDIF

    IF (l_steady_state_stop) THEN
      IF (steady_state_reached(c(:),timesteplen)) THEN
        ! print time at end of integration (same as written to output file)
        !PRINT *, 'steady-state reached after ', &
        !  model_time/OneDay - model_start_day, ' days.'
        WRITE(*,'(A,F6.3,A)') ' steady-state reached after ', &
          (model_time+timesteplen)/OneDay - model_start_day, ' days.'
        ! adapt model_end to stop at end of this time step
        model_end = model_time + timesteplen
      ENDIF
    ENDIF

  CONTAINS

    !-------------------------------------------------------------------------

    SUBROUTINE aerosol_exchng

      ! Aerosol exchange is currently calculated with Euler forward.
      ! This should be changed if numerical problems occur.

      USE messy_main_tools, ONLY: str

      IMPLICIT NONE

      REAL(DP) :: factor
      INTEGER :: i, jb
      REAL :: delta_lwc, csalt ! chamber
      REAL :: tau_aerosol = 7200. ! lifetime  2700 (45) --> 7200 (120)

      DO jb = 1,APN

        IF (l_ff) THEN
          ! (no additional cations and anions from fresh particles
          ! for frostflower model setup)
        ELSE
          ! additional cations and anions from fresh particles:
          factor = timesteplen * exchng(jb)
          ! adjustment for operator splitting:
          factor = factor / (1.-factor)
          IF (ind_NH4p_a(jb)  /=0) &
            c(ind_NH4p_a(jb))  = c(ind_NH4p_a(jb))  + factor * c0_NH4p(jb)
          IF (ind_Nap_a(jb)   /=0) &
            c(IND_Nap_a(jb))   = c(IND_Nap_a(jb))   + factor * c0_Nap(jb)
          IF (ind_FeOHp_a(jb) /=0) &
            c(IND_FeOHp_a(jb)) = c(IND_FeOHp_a(jb)) + factor * c0_FeOHp(jb)
          IF (ind_FeCl2p_a(jb)/=0) &
            c(IND_FeCl2p_a(jb)) = c(IND_FeCl2p_a(jb)) + factor * c0_FeCl2p(jb)
          IF (ind_HCO3m_a(jb) /=0) &
            c(ind_HCO3m_a(jb)) = c(ind_HCO3m_a(jb)) + factor * c0_HCO3m(jb)
          IF (ind_NO3m_a(jb)  /=0) &
            c(ind_NO3m_a(jb))  = c(ind_NO3m_a(jb))  + factor * c0_NO3m(jb)
          IF (ind_Clm_a(jb)   /=0) &
            c(ind_Clm_a(jb))   = c(ind_Clm_a(jb))   + factor * c0_Clm(jb)
          IF (ind_Brm_a(jb)   /=0) &
            c(ind_Brm_a(jb))   = c(ind_Brm_a(jb))   + factor * c0_Brm(jb)
          IF (ind_Im_a(jb)    /=0) &
            c(ind_Im_a(jb))    = c(ind_Im_a(jb))    + factor * c0_Im(jb)
          IF (ind_IO3m_a(jb)  /=0) &
            c(ind_IO3m_a(jb))  = c(ind_IO3m_a(jb))  + factor * c0_IO3m(jb)
          IF (ind_SO4mm_a(jb) /=0) &
            c(ind_SO4mm_a(jb)) = c(ind_SO4mm_a(jb)) + factor * c0_SO4mm(jb)
          IF (ind_HSO4m_a(jb) /=0) &
            c(ind_HSO4m_a(jb)) = c(ind_HSO4m_a(jb)) + factor * c0_HSO4m(jb)
        ENDIF

        ! loss of cations and anions via particle sedimentation
        ! only if xaer=1, otherwise factor=1
        factor = 1. - xaer(jb) * timesteplen * exchng(jb)
        DO i = 1,NSPEC
          ! all species whose names contain '_a##' are lost via
          ! particle sedimentation (except for RR*)
          IF (INDEX(SPC_NAMES(i),'_a'//TRIM(str(jb,'(I2.2)'))) /= 0) THEN
            IF (INDEX(SPC_NAMES(i),'RR') /= 1) THEN
              !PRINT *,'sedi of ', TRIM(SPC_NAMES(i)), ' from bin ', jb
              c(i) = c(i) * factor
            ELSE
              !PRINT *,'no sedi of ', TRIM(SPC_NAMES(i)), ' from bin ', jb
            ENDIF
          ENDIF
        ENDDO

!!$        !----------
!!$
!!$        ! aerosol injection 1
!!$        !start_aerosol = 8.*60.+35. !11:32:14-11:23:39=8:35
!!$        !end_aerosol = 56.*60.+05.  !12:19:44-11:23:39=56:05 --> delta 00:47:30
!!$        start_aerosol = 6360. !01:46:00
!!$        end_aerosol = 9420.  !02:37:00 --> delta 00:51:00
!!$        IF ((model_time >= model_start_day*OneDay+start_aerosol) .AND. &
!!$          (model_time <= model_start_day*OneDay+end_aerosol)) THEN
!!$          ! add 2E-8 mol/mol HONO together with aerosol:
!!$          IF (ind_HONO /= 0) c(ind_HONO) = &
!!$            c(ind_HONO) + 2E-8 * cair * timesteplen/(end_aerosol-start_aerosol)
!!$          ! calc new aerosol addition:
!!$          delta_lwc = timesteplen * &
!!$            ((lwc(jb)/tau_aerosol) + (2.9E-10 / (end_aerosol-start_aerosol)))
!!$          csalt  = 6.1 * delta_lwc * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
!!$          lwc(jb) = lwc(jb) + delta_lwc
!!$          c(IND_Nap_a(jb)) = c(IND_Nap_a(jb)) + 0.99 * csalt
!!$          c(IND_FeOHp_a(jb)) = c(IND_FeOHp_a(jb)) + 0.01 * csalt
!!$          c(ind_Clm_a(jb)) = c(ind_Clm_a(jb)) + 0.95 * csalt
!!$          c(ind_Brm_a(jb)) = c(ind_Brm_a(jb)) + 0.05 * csalt
!!$          ! update cvfac with new LWC:
!!$          cvfac(jb) = 1.E3 / ( N_A * lwc(jb) )
!!$        ENDIF
!!$        
!!$        !----------
!!$
!!$        ! aerosol injection 2
!!$        !start_aerosol = 3.*3600.+51.*60.+55. !15:15:34-11:23:39=3:51:55
!!$        !end_aerosol = 4.*3600.+33.*60.+29.  !15:57:08-11:23:39=4:33:29 --> delta 00:41:34
!!$        start_aerosol = 21131. 
!!$        end_aerosol = 21132.
!!$        IF ((model_time >= model_start_day*OneDay+start_aerosol) .AND. &
!!$          (model_time <= model_start_day*OneDay+end_aerosol)) THEN
!!$          ! add 2E-8 mol/mol HONO together with aerosol:
!!$          IF (ind_HONO /= 0) c(ind_HONO) = &
!!$            c(ind_HONO) + 2E-8 * cair * timesteplen/(end_aerosol-start_aerosol)
!!$          ! calc new aerosol addition:
!!$          delta_lwc = timesteplen * &
!!$            ((lwc(jb)/tau_aerosol) + (2.9E-10 / (end_aerosol-start_aerosol)))
!!$          csalt  = 6.1 * delta_lwc * N_A / 1.E3 ! mol/L -> mcl/cm3(air)
!!$          lwc(jb) = lwc(jb) + delta_lwc
!!$          c(IND_Nap_a(jb)) = c(IND_Nap_a(jb)) + 0.99 * csalt
!!$          c(IND_FeOHp_a(jb)) = c(IND_FeOHp_a(jb)) + 0.01 * csalt
!!$          c(ind_Clm_a(jb)) = c(ind_Clm_a(jb)) + 0.95 * csalt
!!$          c(ind_Brm_a(jb)) = c(ind_Brm_a(jb)) + 0.05 * csalt
!!$          ! update cvfac with new LWC:
!!$          cvfac(jb) = 1.E3 / ( N_A * lwc(jb) )
!!$        ENDIF
!!$
!!$        !----------
!!$        ! aerosol loss:
!!$        ! loss of cations and anions via particle sedimentation
!!$        factor = 1. - timesteplen / tau_aerosol
!!$        DO i = 1,NSPEC
!!$          ! all species whose names contain '_a##' are lost via
!!$          ! particle sedimentation (except for RR*)
!!$          IF (INDEX(SPC_NAMES(i),'_a'//TRIM(str(jb,'(I2.2)'))) /= 0) THEN
!!$            IF (INDEX(SPC_NAMES(i),'RR') /= 1) THEN
!!$              c(i) = c(i) * factor
!!$            ENDIF
!!$          ENDIF
!!$        ENDDO
!!$        ! lwc decreases due to particle loss:
!!$        lwc(jb) = lwc(jb) * factor
!!$        ! update cvfac with new LWC:
!!$        cvfac(jb) = 1.E3 / ( N_A * lwc(jb) )
!!$
!!$        !----------

        IF (l_ff) THEN
          ! lwc decreases due to particle loss:
          lwc(jb) = lwc(jb) * factor
          ! update cvfac with new LWC:
          cvfac(jb) = 1.E3 / ( N_A * lwc(jb) )
        ENDIF

        ! The aqueous-phase concentration of water may have changed in
        ! the sedimentation code above. To get the correct value, it is
        ! determined from the current LWC:
        IF (ind_H2O_a(jb) /= 0) c(ind_H2O_a(jb)) = &
          rho_H2O * lwc(jb) * 1000./M_H2O * N_A/1.E6

      ENDDO

    END SUBROUTINE aerosol_exchng

    !-------------------------------------------------------------------------

    SUBROUTINE check_range(infostring,conc)

      ! print a warning if a concentration is not in the correct range

      CHARACTER(LEN=*), INTENT(IN) :: infostring
      REAL(DP),         INTENT(IN) :: conc(:) ! tracer concentration
      INTEGER :: jt

      INTRINSIC :: SIZE

      tracer_loop: DO jt=1,SIZE(conc)
        IF (SPC_NAMES(jt)(1:2) == "RR") CYCLE ! no checks for reaction rates
        IF (INDEX(SPC_NAMES(jt),"_a03")/=0) CYCLE ! no checks for ocean conc
        wrong_conc: IF ((conc(jt)<0._DP).OR.(conc(jt)>cair)) THEN
          WRITE(*,'(2A,F10.0,A,1PG12.3E3,2A)') infostring, &
            ' time =', model_time, &
            ', c =', conc(jt), ' mcl/cm3 for ', TRIM(SPC_NAMES(jt))
        ENDIF wrong_conc
      ENDDO tracer_loop

    END SUBROUTINE check_range

    !-------------------------------------------------------------------------

  END SUBROUTINE mecca_physc

  !***************************************************************************

  SUBROUTINE mecca_result

    IMPLICIT NONE
    REAL(DP), SAVE, DIMENSION(NSPEC) :: c_old = 0.
    INTEGER :: i

    IF (l_aero) CALL write_output_file(ncid_aero, model_time, lwc)

    output = c/cair
    DO i = 1,NSPEC
      IF (INDEX(SPC_NAMES(i),'RR') == 1) THEN
        ! for the reaction rates RR*, divide the difference to previous
        ! value by the time step length:
        IF (l_RRconc) THEN ! [concentration/time]
          output(i) = (c(i)-c_old(i)) / (timesteplen)
        ELSE ! [mixing ratio/time]
          output(i) = (c(i)-c_old(i)) / (cair*timesteplen)
        ENDIF
      ENDIF
    ENDDO
    c_old = c
    CALL write_output_file(ncid_tracer, model_time, output)

    ! the reaction rates A are from SUBROUTINE Fun in messy_mecca_kpp_function
    CALL write_output_file(ncid_ratesa, model_time, A/cair)

#ifdef MECCA_TAG
    IF (l_tag) CALL mecca_tag_result(model_time)
#endif

  END SUBROUTINE mecca_result

  !***************************************************************************

  SUBROUTINE mecca_finish

    USE messy_main_tools, ONLY: find_next_free_unit
    IMPLICIT NONE
    INTEGER :: ncid_c_end, ncid_k_end, unit, i
    CHARACTER(LEN=20), PARAMETER :: k_unit(NREACT) = '(cm^3)^(1-N)/s'

    ! write final values of reaction rates to file:
    unit = find_next_free_unit(100,200)
    OPEN(unit, FILE='caaba_mecca_a_end.dat', status='UNKNOWN')
    DO i = 1,NREACT
      WRITE (unit,*) A(i) ! from messy_mecca_kpp_function
    ENDDO
    CLOSE(unit)

    ! write final values of c (c_end) to nc file:
    CALL open_output_file(ncid_c_end, 'caaba_mecca_c_end', &
      SPC_NAMES, c_unit)
    CALL write_output_file(ncid_c_end, model_time, output)
    CALL close_file(ncid_c_end)

    ! write final values of rconst (k_end) to nc file:
    CALL open_output_file(ncid_k_end, 'caaba_mecca_k_end', &
      EQN_TAGS, k_unit)
! op_pj_20160218+: The next line prohibits compilation of CAABA
!                  when a vectorized solver has been selected. In that
!                  case, rconst has rank 2 but the first dimension length
!                  is 1 (only one box in CAABA). A possible solution is
!                  therefore:
!!$ CALL write_output_file(ncid_k_end, model_time, rconst)
    CALL write_output_file(ncid_k_end, model_time, PACK(rconst,.TRUE.))
! op_pj_20160218-
    CALL close_file(ncid_k_end)

    IF (REQ_MCFCT) THEN
      ! create ferret *.jnl files:
      unit = find_next_free_unit(100,200)
      ! plot histogram of all rate coefficients:
      OPEN(unit, FILE='_histogram_k.jnl', status='UNKNOWN')
      DO i = 1,NREACT
        WRITE (unit,*) 'GO frequency_histogram_rs ', TRIM(EQN_TAGS(i))
      ENDDO
      WRITE (unit,*) 'GO newpage'
      CLOSE(unit)
      ! make plots for all species:
      OPEN(unit, FILE='_scatterplot1.jnl', status='UNKNOWN')
      DO i = 1,NSPEC
        WRITE (unit,*) 'GO _scatterplot2 ', TRIM(SPC_NAMES(i))
      ENDDO
      CLOSE(unit)
      ! plot against all rate coefficients:
      OPEN(unit, FILE='_scatterplot2.jnl', status='UNKNOWN')
      WRITE (unit,*) 'GO frequency_histogram_rs ($1)'
      DO i = 1,NREACT
        WRITE (unit,*) 'GO scatterplot_mc ($1) ', TRIM(EQN_TAGS(i))
      ENDDO
      WRITE (unit,*) 'GO newpage'
      CLOSE(unit)
    ENDIF

    IF (l_aero) CALL close_file(ncid_aero)

    CALL close_file(ncid_tracer)

    CALL close_file(ncid_ratesa)

#ifdef MECCA_TAG
    IF (l_tag) CALL mecca_tag_finish
#endif

    DEALLOCATE(c0_HCO3m, c0_NO3m, c0_Clm, c0_Brm, c0_Im, c0_IO3m, &
      c0_SO4mm, c0_HSO4m)
    DEALLOCATE(xaer, radius, lwc, csalt, c0_NH4p, c0_Nap, c0_FeOHp, c0_FeCl2p, exchng)

  END SUBROUTINE mecca_finish

  !***************************************************************************
  ! PRIVATE ROUTINES
  !***************************************************************************

  SUBROUTINE define_mcexp(status, mcexp)

    USE messy_main_rnd,   ONLY: RND_MTW_GAUSS, rnd_init, rnd_number, rnd_finish
    USE messy_mecca,      ONLY: mcexp_seed

    IMPLICIT NONE
    INTEGER,                INTENT(OUT) :: status
    REAL(DP), DIMENSION(:), INTENT(OUT) :: mcexp
    INTEGER :: id_rnd

    ! assign a set of normally distributed random numbers to mcexp:
    CALL rnd_init(status, id_rnd, RND_MTW_GAUSS, mcexp_seed)
    IF (status/=0) RETURN
    CALL rnd_number(id_rnd, mcexp(:))
    CALL rnd_finish(id_rnd)

  END SUBROUTINE define_mcexp

END MODULE messy_mecca_box

!*****************************************************************************
