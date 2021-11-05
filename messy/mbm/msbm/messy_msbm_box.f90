MODULE messy_msbm_box
!-----
! DESCRIPTION
! Interface layer for msbm box model, analog messy_msbm_e5.f90
!
! AUTHOR
! Joachim Buchholz, Max Planck Institute for Chemistry, Mainz, Germany
!
! LAST CHANGES
! 13. September 2004, J. Buchholz
!-----

  USE messy_msbm, ONLY: dp

  IMPLICIT NONE
  PRIVATE


  !-----
  ! module variables which have to be initialised
  !-----
  REAL(dp) :: TEMP, PRESS, H2SO4_t0, &
              H2O_t0, HNO3_t0, &
              HBr_t0, HOBr_t0, BrNO3_t0, &
              HCl_t0, HOCl_t0, ClNO3_t0
  LOGICAL :: LCalcChem
  INTEGER :: phase

  NAMELIST /input/ TEMP, PRESS, H2SO4_t0, &
                   H2O_t0, HNO3_t0, &
                   HBr_t0, HOBr_t0, BrNO3_t0, &
                   HCl_t0, HOCl_t0, ClNO3_t0, &
                   LCalcChem, &
                   phase

  !-----
  ! module variables for output
  !-----
  REAL(dp) :: H2O_l, H2O_i, H2O_g, &
              HNO3_g, HNO3_n, HNO3_l, &
              HCl_g, HCl_l, &
              HOCl_g, HOCl_l, &
              HBr_g, HBr_l, &
              HOBr_g, HOBr_l, &
              N_solid, r_s, vel, &
              A_l, rmean_l
  REAL(dp) :: khet1,  khet2,  khet3,  khet4,  khet5,  &
              khet6,  khet7,  khet8,  khet9,  khet10, &
              khet11, khet12, khet13, khet14, khet15, &
              khet16, khet17, khet18, khet19, khet20, &
              khet21, khet22, khet23, khet24, khet25, &
              khet26, khet27, khet28, khet29, khet30
  INTEGER :: i_val


  PUBLIC :: msbm_initialize    ! checks the phase of strat. particles
  PUBLIC :: msbm_physc         ! program to calculate composition and
                              ! het. reaction rates of PSCs
  PUBLIC :: msbm_output

CONTAINS

!=============================================================================

SUBROUTINE  msbm_initialize

  USE messy_msbm, ONLY: msbm_read_nml_ctrl
  
  IMPLICIT NONE
  
  INTEGER :: status, iou
  
  iou=find_next_free_unit(100,200)
  CALL msbm_read_nml_ctrl(status,iou)   ! read /CTRL/
  IF (status/=0) THEN
    write (*,*) 'msbm_read_nml_ctrl failed with status', status
    STOP
  END IF
  
  iou=find_next_free_unit(100,200)
  open(iou,file='input.nml')
  read(unit=iou,NML=input)

  write (*,*) 'Input namelist:'
  write (*,*) '---------------'
  write (*,NML=input)
  write (*,*)



  CONTAINS

    FUNCTION find_next_free_unit(istart,istop) RESULT(unit)
      !-----
      ! based on a function in mo_filename.f90 from L. Kornblueh, I. Kirchner, 
      ! and A. Rhodin
      !-----
      INTEGER :: istart, istop, unit
      LOGICAL :: found, opened
      INTEGER :: i

      found = .FALSE.
      DO i=istart,istop
         INQUIRE(unit=i,opened=opened)
         IF (.NOT.opened) THEN
            unit = i
            found = .TRUE.
            EXIT
         END IF
      END DO

      IF (.NOT. found) THEN
        WRITE(*,'(a,i2.2,a,i2.2,a)') &
           'No unit in range <',istart,':',istop,'> free.'
        STOP
      END IF
    END FUNCTION find_next_free_unit

END SUBROUTINE msbm_initialize

!=============================================================================

SUBROUTINE msbm_physc

  USE messy_msbm, ONLY: &
    mz_psc_phase,                                                     &
    mz_psc_ice_H2O, mz_psc_nat_HNO3,                                  &
    mz_psc_liq_bH2SO4b, mz_psc_liq_bHNO3b, mz_psc_liq_bH2SO4,         &
    mz_psc_liq_bHNO3, mz_psc_liq_partHNO3, mz_psc_liq_wnen,           &
    mz_psc_liq_wHNO3, mz_psc_liq_wH2SO4, mz_psc_liq_hHCl,             &
    mz_psc_liq_hHBr, mz_psc_liq_hHOCl, mz_psc_liq_hHOBr,              &
    mz_psc_liq_partHBr, mz_psc_liq_wHCl,                              &
    mz_psc_liq_wHOCl, mz_psc_liq_wHOBr,                               &
    mz_psc_liq_H2O, mz_psc_liq_HCl,                                   &
    mz_psc_liq_HOCl, mz_psc_liq_HOBr,                                 &
    mz_psc_N_solid, mz_psc_r_solid,                                   &
    mz_psc_dim_liquid, mz_psc_surface_liquid,                         &
    mz_psc_het_nat1, mz_psc_het_nat2, mz_psc_het_nat3,                &
    mz_psc_het_nat4, mz_psc_het_nat5, mz_psc_het_nat6,                &
    mz_psc_het_nat7, mz_psc_het_nat8, mz_psc_het_nat9,                &
    mz_psc_het_nat10, mz_psc_het_nat11,                               &
    mz_psc_het_ice12, mz_psc_het_ice13, mz_psc_het_ice14,             &
    mz_psc_het_ice15, mz_psc_het_ice16, mz_psc_het_ice17,             &
    mz_psc_het_ice18, mz_psc_het_ice19, mz_psc_het_ice20,             &
    mz_psc_het_ice21, mz_psc_het_ice22,                               &
    mz_psc_het_liq23, mz_psc_het_liq24, mz_psc_het_liq25,             &
    mz_psc_het_liq26, mz_psc_het_liq27, mz_psc_het_liq28,             &
    mz_psc_het_liq29, mz_psc_het_liq30,                               &
    mz_psc_liq_check,                                                 &
    mz_psc_diff, mz_psc_density,                                      &
    mz_psc_vel

  IMPLICIT NONE
  
  REAL(dp) :: H2O_gl, HNO3_gl, H2SO4_l
  REAL(dp) :: partHNO3, partHBr
  REAL(dp) :: bHNO3, bH2SO4, bHNO3b, bH2SO4b
  REAL(dp) :: wHNO3, wH2SO4, wHCl, wHOCl, wHOBr, wnen
  REAL(dp) :: hHCl, hHBr, hHOCl, hHOBr
  REAL(dp), PARAMETER :: sigmaaero=1.8_dp
  REAL(dp) :: diff, dens

  !----------------------------
  !Treatment of solid particles
  !----------------------------
  !-----
  ! determine for thermodynamically stable phases
  !-----
  phase = mz_psc_phase(TEMP,           &
                       PRESS, H2O_t0,  &
                       HNO3_t0, phase, &
                       0.0_dp, 0.0_dp)

  !-----
  ! calculate partitioning between gas phase, ice and NAT
  !-----
  H2O_i   = mz_psc_ice_H2O(phase,   &
                       TEMP, PRESS, &
                       H2O_t0)
  H2O_gl  = H2O_t0-H2O_i
  H2O_l   = 0.0_dp
  H2O_g   = H2O_gl
  HNO3_n  = mz_psc_nat_HNO3(phase,     &
                       TEMP, PRESS,    &
                       H2O_t0, HNO3_t0)
  HNO3_gl = HNO3_t0 - HNO3_n
  HNO3_l  = 0.0_dp
  HNO3_g  = HNO3_gl

  IF (phase==3 .OR. phase==2) THEN
    !-----
    ! calculate size of ice and NAT particles
    !-----
    N_solid = mz_psc_N_solid(phase,       &
                             TEMP, PRESS, &
                             H2O_i, HNO3_n)
    r_s = mz_psc_r_solid(phase,          &
                         TEMP, PRESS,    &
                         N_solid, H2O_i, &
                         HNO3_n)
  ELSE   ! IF (phase==3 .OR. phase==2) THEN ...
    r_s     = 0.0_dp
  END IF   ! IF (phase==3 .OR. phase==2) THEN ...

  !------------------------------
  ! Treatment of liquid particles
  !------------------------------
  !-----
  ! check conditions wether STS parameterisation holds
  !-----
  i_val=mz_psc_liq_check(TEMP,       &
                    PRESS, H2O_g,    &
                    HNO3_g, H2SO4_t0)

  !-----
  ! calculate partitioning (fraction of HX in the gas phase),
  ! henrys law coefficients and aerosolcomposition
  !-----
  bH2SO4b = mz_psc_liq_bH2SO4b(TEMP,       &
                        PRESS, H2O_gl)
  bHNO3b  = mz_psc_liq_bHNO3b(TEMP,        &
                        PRESS, H2O_gl)
  bH2SO4  = mz_psc_liq_bH2SO4(TEMP,        &
                        PRESS, H2SO4_t0,   &
                        H2O_gl, HNO3_gl,   &
                        bH2SO4b, bHNO3b)
  bHNO3   = mz_psc_liq_bHNO3(TEMP,         &
                        PRESS, H2O_gl,     &
                        HNO3_gl, bHNO3b,   &
                        bH2SO4b, bH2SO4)
  partHNO3= mz_psc_liq_partHNO3(TEMP,      &
                        PRESS, H2O_gl,     &
                        HNO3_gl, bHNO3,    &
                        bH2SO4, bH2SO4b)
  wnen    = mz_psc_liq_wnen(TEMP,          &
                        bHNO3, bH2SO4,     &
                        bH2SO4b, partHNO3)
  wHNO3   = mz_psc_liq_wHNO3(TEMP,         &
                        bHNO3, bH2SO4,     &
                        bH2SO4b, partHNO3, &
                        wnen)
  wH2SO4  = mz_psc_liq_wH2SO4(TEMP,        &
                        bH2SO4, bH2SO4b,   &
                        partHNO3, wnen)
  hHCl    = mz_psc_liq_hHCl(TEMP,          &
                        wHNO3, wH2SO4,     &
                        wnen)
  hHBr    = mz_psc_liq_hHBr(TEMP,          &
                        wHNO3, wH2SO4,     &
                        wnen)
  hHOCl   = mz_psc_liq_hHOCl(TEMP,         &
                        bH2SO4, bHNO3)
  hHOBr   = mz_psc_liq_hHOBr(hHOCl)
  partHBr = mz_psc_liq_partHBr(PRESS,      &
                        HBr_t0, H2SO4_t0,  &
                        hHBr, bH2SO4)
  wHCl    = mz_psc_liq_wHCl(PRESS,         &
                        HCl_t0, H2SO4_t0,  &
                        hHCl, bH2SO4,      &
                        wnen)
  wHOCl   = mz_psc_liq_wHOCl(PRESS,        &
                        HOCl_t0, H2SO4_t0, &
                        hHOCl, bH2SO4,     &
                        wnen)
  wHOBr   = mz_psc_liq_wHOBr(PRESS,        &
                        HOBr_t0, H2SO4_t0, &
                        hHOBr, bH2SO4,     &
                        wnen)

  !-----
  ! calculate new gas and liquid phase mixing ratios
  !-----
  H2SO4_l = H2SO4_t0
  H2O_l   = mz_psc_liq_H2O(wH2SO4, wHNO3,  &
                        H2SO4_l, H2O_t0)
  H2O_g   = H2O_gl-H2O_l
  HNO3_l  = (1.0_dp - partHNO3)            &
                       *HNO3_gl
  ! gas phase depletion by HNO3 uptake:
  HNO3_g  = HNO3_gl - HNO3_l
  HCl_l   = mz_psc_liq_HCl(wHCl,           &
                        H2SO4_l, wH2SO4,   &
                        HCl_t0)
  HCl_g   = HCl_t0 - HCl_l
  HOCl_l  = mz_psc_liq_HOCl(wHOCl,         &
                        H2SO4_l, wH2SO4,   &
                        HOCl_t0)
  HOCl_g  = HOCl_t0 - HOCl_l
  HBr_l   = (1.0_dp - partHBr)             &
                       *HBr_t0
  HBr_g   = HBr_t0 - HBr_l
  HOBr_l  = mz_psc_liq_HOBr(wHOBr,         &
                        H2SO4_l, wH2SO4,   &
                        HOBr_t0)
  HOBr_g  = HOBr_t0 - HOBr_l

  IfLCalcChem1: IF (LCalcChem) THEN
    !-----
    ! calculate composition dependent size distribution
    !-----
    A_l     = mz_psc_surface_liquid(TEMP, &
                     PRESS, H2O_l,        &
                     HNO3_l, H2SO4_l)
    rmean_l = mz_psc_dim_liquid(A_l,      &
                     sigmaaero)

    !-----
    ! calculate heterogeneous reaction rates
    !-----
    dens  = mz_psc_density(wH2SO4,         &
                        wHNO3, TEMP)
    diff  = mz_psc_diff(TEMP,              &
                        dens, bHNO3,       &
                        bH2SO4)
    khet23= mz_psc_het_liq23(TEMP,         &
                        PRESS, HCl_t0,     &
                        HOCl_t0, bH2SO4,   &
                        bHNO3, wHCl,       &
                        dens, hHOCl,       &
                        A_l, rmean_l,      &
                        diff)
    khet24= mz_psc_het_liq24(TEMP,         &
                        PRESS, H2O_g,      &
                        HCl_t0, ClNO3_t0,  &
                        wHCl, dens,        &
                        A_l, rmean_l)
    khet25= mz_psc_het_liq25(TEMP,         &
                        PRESS, H2O_g,      &
                        wHCl, dens,        &
                        A_l, rmean_l)
    khet26= mz_psc_het_liq26(TEMP,         &
                        PRESS, H2O_g,      &
                        A_l)
    khet27= mz_psc_het_liq27(TEMP,         &
                        PRESS, HCl_t0,     &
                        HOBr_t0, bH2SO4,   &
                        bHNO3, wHCl,       &
                        dens, hHOBr,       &
                        A_l, rmean_l,      &
                        diff)
    khet28= mz_psc_het_liq28(TEMP,         &
                        PRESS, HBr_t0,     &
                        HOBr_t0, partHBr,  &
                        bH2SO4, bHNO3,     &
                        wHOBr, dens,       &
                        hHBr, hHOBr,       &
                        A_l, rmean_l,      &
                        diff)
    khet29= mz_psc_het_liq29(TEMP,         &
                        PRESS, HBr_t0,     &
                        HOCl_t0, partHBr,  &
                        bH2SO4, bHNO3,     &
                        wHOCl, dens,       &
                        hHOCl, hHBr,       &
                        A_l, rmean_l,      &
                        diff)
    khet30= mz_psc_het_liq30(TEMP,         &
                        PRESS, H2O_g,      &
                        A_l)
  END IF IfLCalcChem1

  IfLCalcChemPhase2: IF (LCalcChem .AND. phase==2) THEN
    !-----
    ! calculate het. reaction rates on nat particles
    !-----
    khet1= mz_psc_het_nat1(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet2= mz_psc_het_nat2(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, H2O_g)
    khet3= mz_psc_het_nat3(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet4= mz_psc_het_nat4(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet5= mz_psc_het_nat5(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, H2O_g)
    khet6= mz_psc_het_nat6(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, HBr_g)
    khet7= mz_psc_het_nat7(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet8= mz_psc_het_nat8(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, HBr_g)
    khet9= mz_psc_het_nat9(TEMP,        &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet10= mz_psc_het_nat10(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HBr_g)
    khet11= mz_psc_het_nat11(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, H2O_g)
  ELSEIF (LCalcChem .AND. phase/=2) THEN IfLCalcChemPhase2
    !-----
    ! Set reaction rates for reactions on nat to zero if there is no nat,
    ! otherwise values from previous time steps would persist.
    !-----
    khet1   = 0.0_dp
    khet2   = 0.0_dp
    khet3   = 0.0_dp
    khet4   = 0.0_dp
    khet5   = 0.0_dp
    khet6   = 0.0_dp
    khet7   = 0.0_dp
    khet8   = 0.0_dp
    khet9   = 0.0_dp
    khet10  = 0.0_dp
    khet11  = 0.0_dp
  END IF IfLCalcChemPhase2

  IfLCalcChemPhase3: IF (LCalcChem .AND. phase==3) THEN
    !-----
    ! assume, that even though solid ptcls are an ice-NAT mixture
    ! the het. react. are similar to those on a pure ice surface
    !-----
    khet12= mz_psc_het_ice12(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet13= mz_psc_het_ice13(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, H2O_g)
    khet14= mz_psc_het_ice14(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet15= mz_psc_het_ice15(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet16= mz_psc_het_ice16(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, H2O_g)
    khet17= mz_psc_het_ice17(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HBr_g)
    khet18= mz_psc_het_ice18(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet19= mz_psc_het_ice19(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HBr_g)
    khet20= mz_psc_het_ice20(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HCl_g)
    khet21= mz_psc_het_ice21(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, HBr_g)
    khet22= mz_psc_het_ice22(TEMP,      &
                          PRESS, r_s,   &
                          N_solid, H2O_g)
  ELSEIF (LCalcChem .AND. phase/=3) THEN IfLCalcChemPhase3
    !-----
    ! Set reaction rates for heterogeneous chemical reactions on ice to
    ! zero if there is no ice, otherwise values from previous time
    ! steps would persist.
    ! Note that nullifying reaction rates is unnecessary if
    ! LCalcChem==.false., however, it does not cause problems either.
    !-----
    khet12  = 0.0_dp
    khet13  = 0.0_dp
    khet14  = 0.0_dp
    khet15  = 0.0_dp
    khet16  = 0.0_dp
    khet17  = 0.0_dp
    khet18  = 0.0_dp
    khet19  = 0.0_dp
    khet20  = 0.0_dp
    khet21  = 0.0_dp
    khet22  = 0.0_dp
  END IF IfLCalcChemPhase3
  
  !-----
  ! calculation of sedimentation velocity
  !-----
  vel = mz_psc_vel(PRESS*100.0_dp, TEMP, r_s)

END SUBROUTINE msbm_physc

!=============================================================================

SUBROUTINE msbm_output

  write (*,*) 'Microphysical data:'
  write (*,*) '-------------------'
  write (*,*) 'phase   = ', phase
  write (*,*) 'N_solid = ', N_solid
  write (*,*) 'r_s     = ', r_s
  write (*,*) 'vel     = ', vel
  write (*,*) 'i_val   = ', i_val
  write (*,*) 'A_l     = ', A_l
  write (*,*) 'rmean_l = ', rmean_l
  write (*,*)

  write (*,*) 'Tracer data:'
  write (*,*) '------------'
  write (*,*) 'H2O_l  = ', H2O_l
  write (*,*) 'H2O_i  = ', H2O_i
  write (*,*) 'H2O_g  = ', H2O_g
  write (*,*) 'HNO3_l = ', HNO3_l
  write (*,*) 'HNO3_n = ', HNO3_n
  write (*,*) 'HNO3_g = ', HNO3_g
  write (*,*) 'HCl_l  = ', HCl_l
  write (*,*) 'HCl_g  = ', HCl_g
  write (*,*) 'HOCl_l = ', HOCl_l
  write (*,*) 'HOCl_g = ', HOCl_g
  write (*,*) 'HBr_l  = ', HBr_l
  write (*,*) 'HBr_g  = ', HBr_g
  write (*,*) 'HOBr_l = ', HOBr_l
  write (*,*) 'HOBr_g = ', HOBr_g
  write (*,*)

  write (*,*) 'Heterogeneous reaction rates:'
  write (*,*) '-----------------------------'
  write (*,*) 'khet1  = ', khet1
  write (*,*) 'khet2  = ', khet2
  write (*,*) 'khet3  = ', khet3
  write (*,*) 'khet4  = ', khet4
  write (*,*) 'khet5  = ', khet5
  write (*,*) 'khet6  = ', khet6
  write (*,*) 'khet7  = ', khet7
  write (*,*) 'khet8  = ', khet8
  write (*,*) 'khet9  = ', khet9
  write (*,*) 'khet10 = ', khet10
  write (*,*) 'khet11 = ', khet11
  write (*,*) 'khet12 = ', khet12
  write (*,*) 'khet13 = ', khet13
  write (*,*) 'khet14 = ', khet14
  write (*,*) 'khet15 = ', khet15
  write (*,*) 'khet16 = ', khet16
  write (*,*) 'khet17 = ', khet17
  write (*,*) 'khet18 = ', khet18
  write (*,*) 'khet19 = ', khet19
  write (*,*) 'khet20 = ', khet20
  write (*,*) 'khet21 = ', khet21
  write (*,*) 'khet22 = ', khet22
  write (*,*) 'khet23 = ', khet23
  write (*,*) 'khet24 = ', khet24
  write (*,*) 'khet25 = ', khet25
  write (*,*) 'khet26 = ', khet26
  write (*,*) 'khet27 = ', khet27
  write (*,*) 'khet28 = ', khet28
  write (*,*) 'khet29 = ', khet29
  write (*,*) 'khet30 = ', khet30
  write (*,*)

END SUBROUTINE msbm_output

!=============================================================================

END MODULE messy_msbm_box
