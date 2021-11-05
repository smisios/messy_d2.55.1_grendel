! ==============================================================================
! {%TAG}_box
! generated: {%TIMEDATE}
!
! this module is created automatically by imtag tool - do not edit
!
! maintenance routines for budgeting configurations (isotopes)
! level: smil boxmodel
!
! {$TAG_INFO} ! this is a template file for imtag utility
!
! [Gromov, MPIC, 2007-2020]
! ==============================================================================

! - general tagging parameters (as conditional defines) ------------------------

#include "{%CMODEL}_tag_parameters.inc"

! ------------------------------------------------------------------------------

! {$CONF_PARAM}

module {%CMODEL}_{%TAG}_box

  use messy_mecca_kpp     ! dp, ... nreact, nspec, ind_*, SPC_NAMES, EQN_TAGS
  use caaba_io,           only: open_output_file, write_output_file, close_file
  use caaba_mem,          only: temp, cair, press, &
#ifdef __GFORTRAN__
                                C_ => &
#endif
                                      C
#ifdef __GFORTRAN__    /* newer gfortran versions have a run-time bug of    */
#define C C_           /* module-level import of C(:), so it has to be      */
#define c C_           /* renamed into something else, e.g. C_(:)           */
#endif

  use {%CMODEL}_tag_common ! common routines
  use {%CMODEL}_{%TAG}     ! SMCL routines

  implicit none

! netcdf handle for deltas, conc., etc. output
  integer :: ncid_{%TAG}

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
! output array: isotope hydrogen: 2 species and d2H (+TH)     = (NSPEC+1)*3
!                                TH(regular),                 +1
!                                d0TH(reg), d0TH(iso), d0d2HTH, +3
!                                NREJCT                       +1
  real(dp)            :: D{%ATOM}out(({%NSPEC}+1)*(2+1)+1+3+1)
  real(dp)            :: d2H({%NSPEC}), d2HTH

! total budget verification
  real(dp)            :: d2HTH0
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
! output array: tagged carbons: 2 species and d13C (+TC)     = (NSPEC+1)*3
!                                TC(regular),                 +1
!                                d0TC(reg), d0TC(iso), d0d13CTC, +3
!                                NREJCT                       +1
  real(dp)            :: D{%ATOM}out(({%NSPEC}+1)*(2+1)+1+3+1)
  real(dp)            :: d13C({%NSPEC}), d13CTC

! total budget verification
  real(dp)            :: d13CTC0
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
! output array: tagged oxygen: 3 species, d17O,d18O,DC17O (+to) = (NSPEC+1)*6
!                               to(regular),
!                               d0TO(reg),d0TO(iso),d0D17TO,d0d18TO,d0DC17OTO,
!                               NREJCT
  real(dp)            :: D{%ATOM}out(({%NSPEC}+1)*(3+3)+1+5+1)
  real(dp)            :: d17O({%NSPEC}), d18O({%NSPEC}), DC17O({%NSPEC}), &
                         d17OTO, d18OTO, DC17OTO

! total budget verification
  real(dp)            :: d17OTO0, d18OTO0, DC17OTO0
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
! totals
  real(dp)            :: D{%A}out(({%NSPEC}*2)+5)
! class fractions
  real(dp)            :: F({%NSPEC},{%NCLASS})
  real(dp)            :: PMC({%NSPEC}), PMCT     ! radiocarbon content in pMC


! total stable/radio C budget
  real(dp)            :: PMCT0
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
! totals
  real(dp)            :: D{%A}out(({%NSPEC})*{%NCLASS}+5)   ! +1
! class fractions
  real(dp)            :: F({%NSPEC},{%NCLASS})
  real(dp)            :: FT({%NCLASS})
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}

! total budget verification
  real(dp)            :: T{%A}0(0:{%NCLASS})

! -----------------------------------------------------------------------------

  public {%TAG}_x0
  public {%TAG}_f0
  public {%TAG}_emis
  public {%TAG}_depos
  public {%TAG}_pmix
  public {%TAG}_set
  public {%TAG}_postprocess
  public {%TAG}_calcdeltas
  public {%TAG}_init
  public {%TAG}_result
  public {%TAG}_finish

! ==============================================================================

contains

! ==============================================================================

  subroutine {%TAG}_x0()

    implicit none

  ! tracers mixing ratios initialization (x0)

    integer :: i

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>case:REM}
! {$x0} [%n%]  (%    d13C({%TAG}_@) = $%)
! n - # of the class, i.e. 2 for d13 / 2 for d17, 3 for d18
! @ - species name; $ - init. value; # - class no
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<case:REM}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
! {$x0} [%2%]  (%    d2H({%TAG}_@) = $%)

#ifdef ZERO_TEST
    d2H(:) = 0.0_dp
#endif

#ifdef UNIT_DELTAPERMIL
  ! 1H, 2H through delta and regular species:

    d2H(:) = d2H(:) / {%TAG}_ufac      ! de-permilizing (if pm units are used)

    C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * &
      isofrac2a(d2H(:), Rstd_2H, {%NQATOM}(:))
    C({%RSIND}(:,2)) = C({%RSIND}(:,0)) * &
      isofrac2r(d2H(:), Rstd_2H, {%NQATOM}(:))
#endif
#ifdef UNIT_FRACMIN
  ! 1H, 2H through minor fraction and regular species:

     C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * (1.0_dp - d2H(:))
     C({%RSIND}(:,2)) = C({%RSIND}(:,0)) * d2H(:)
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
! {$x0} [%2%]  (%    d13C({%TAG}_@) = $%)

#ifdef ZERO_TEST
    d13C(:) = 0.0_dp
#endif

#ifdef UNIT_DELTAPERMIL
  ! 12C, 13C through delta and regular species:

    d13C(:) = d13C(:) / {%TAG}_ufac    ! de-permilizing (if pm units are used)

    C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * &
      isofrac2a(d13C(:), Rstd_13C, {%NQATOM}(:))
    C({%RSIND}(:,2)) = C({%RSIND}(:,0)) * &
      isofrac2r(d13C(:), Rstd_13C, {%NQATOM}(:))
#endif
#ifdef UNIT_FRACMIN
  ! 12C, 13C through minor fraction and regular species:

     C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * (1.0_dp - d13C(:))
     C({%RSIND}(:,2)) = C({%RSIND}(:,0)) * d13C(:)
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}

#ifndef ONLY_MINOR
! {$x0} [%2%]  (%    d17O({%TAG}_@) = $%)
! {$x0} [%3%]  (%    d18O({%TAG}_@) = $%)
#else
! {$x0} [%1%]  (%    d17O({%TAG}_@) = $%)
! {$x0} [%2%]  (%    d18O({%TAG}_@) = $%)
#endif

#ifdef ZERO_TEST
    d17O(:) = 0.0_dp
    d18O(:) = 0.0_dp
#endif

#ifdef UNIT_DELTAPERMIL
  ! 16O, 17O, 18O through delta and regular species:

    d17O(:) = d17O(:) / {%TAG}_ufac    ! de-permilizing (if pm units are used)
    d18O(:) = d18O(:) / {%TAG}_ufac

#ifndef ONLY_MINOR
    C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * &
      isofrac3a(d17O(:), Rstd_17O, d18O(:), Rstd_18O, {%NQATOM}(:))
    C({%RSIND}(:,2)) = C({%RSIND}(:,0)) * &
      isofrac3r(d17O(:), Rstd_17O, d18O(:), Rstd_18O, {%NQATOM}(:))
    C({%RSIND}(:,3)) = C({%RSIND}(:,0)) * &
      isofrac3r(d18O(:), Rstd_18O, d17O(:), Rstd_17O, {%NQATOM}(:))
#else
    C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * &
      isofrac3r(d17O(:), Rstd_17O, d18O(:), Rstd_18O, {%NQATOM}(:))
    C({%RSIND}(:,2)) = C({%RSIND}(:,0)) * &
      isofrac3r(d18O(:), Rstd_18O, d17O(:), Rstd_17O, {%NQATOM}(:))
#endif

#endif
#ifdef UNIT_FRACMIN
  ! 16O, 17O, 18O through minor fractions and regular species:

    C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * (1.0_dp - (d17O(:) + d18O(:)))
    C({%RSIND}(:,2)) = C({%RSIND}(:,0)) * d17O(:)
    C({%RSIND}(:,3)) = C({%RSIND}(:,0)) * d18O(:)
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
#ifdef INIUNIT_PMC
  ! initialising using specified pMC values given in cfg

! {$x0} [%#%]  (%    PMC({%TAG}_@) = $%)

    C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * isofrac2r_pMC(PMC(:), {%NQATOM}(:))
#else
  ! initialising using fractions given in cfg

! {$x0} [%#%]  (%    C({%RSIND}({%TAG}_@,#)) = C({%RSIND}({%TAG}_@,0)) * $%)
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
  ! initialising using fractions given in cfg

! {$x0} [%#%]  (%    C({%RSIND}({%TAG}_@,#)) = C({%RSIND}({%TAG}_@,0)) * $%)

-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}

#ifdef NULL_TEST
  ! minors are initialized emptied

#ifndef CLASSES_1
    C({%RSIND}(:,1)) = C({%RSIND}(:,0))
    C({%RSIND}(:,2:{%NISO})) = 0.0_dp
#else
    C({%RSIND}(:,1)) = 0.0_dp
#endif

#endif

  ! updating total {%ATOM} in the system
    call {%TAG}_calctotals(C)
    call {%TAG}_calcdeltas()

    T{%A}0(:) = T{%A}(:)

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
    d2HTH0 = d2HTH
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
    d13CTC0 = d13CTC
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
    d17OTO0 = d17OTO
    d18OTO0 = d18OTO
    DC17OTO0 = DC17OTO
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
    PMCT0 = PMCT
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}

  end subroutine {%TAG}_x0



! -----------------------------------------------------------------------------

  subroutine {%TAG}_f0()

    implicit none

  ! tracers mixing ratios initialization (f0)

    integer :: i

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>case:REM}
! {$f0} [%n%]  (%    d13C({%TAG}_@) = $%)
! n - # of the class, i.e. 2 for d13 / 2 for d17, 3 for d18
! @ - species name; $ - init. value; # - class no
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<case:REM}

  ! exit if there is no fixed species
    if ({%NFIX} .lt. 1) return

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
! {$f0} [%2%]  (%    d2H({%TAG}_@) = $%)

#ifdef UNIT_DELTAPERMIL
  ! 1H, 2H through delta and regular species:
    d2H({%FSIND}(:,0)) = d2H({%FSIND}(:,0)) / {%TAG}_ufac  ! de-permilizing (if pm units are used)

    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * isofrac2a(d2H({%FSIND}(:,0)), Rstd_2H, {%NQATOM}({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * isofrac2r(d2H({%FSIND}(:,0)), Rstd_2H, {%NQATOM}({%FSIND}(:,0)))
#endif
#ifdef UNIT_FRACMIN
  ! 1H, 2H through minor fraction and regular species:

    C({%FSIND}(:),1)) = C({%FSIND}(:),0) * (1.0_dp - d2H({%FSIND}(:,0)))
    C({%FSIND}(:),2)) = C({%FSIND}(:),0) *           d2H({%FSIND}(:,0))
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
! {$f0} [%2%]  (%    d13C({%TAG}_@) = $%)

#ifdef UNIT_DELTAPERMIL
  ! 12C, 13C through delta and regular species:
    d13C({%FSIND}(:,0)) = d13C({%FSIND}(:,0)) / {%TAG}_ufac  ! de-permilizing (if pm units are used)

    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * isofrac2a(d13C({%FSIND}(:,0)), Rstd_13C, {%NQATOM}({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * isofrac2r(d13C({%FSIND}(:,0)), Rstd_13C, {%NQATOM}({%FSIND}(:,0)))
#endif
#ifdef UNIT_FRACMIN
  ! 12C, 13C through minor fraction and regular species:
    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * (1.0_dp - d13C({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) *           d13C({%FSIND}(:,0))
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}

#ifndef ONLY_MINOR
! {$f0} [%2%]  (%    d17O({%TAG}_@) = $%)
! {$f0} [%3%]  (%    d18O({%TAG}_@) = $%)
#else
! {$f0} [%1%]  (%    d17O({%TAG}_@) = $%)
! {$f0} [%2%]  (%    d18O({%TAG}_@) = $%)
#endif

#ifdef UNIT_DELTAPERMIL
  ! 16O, 17O, 18O through delta and regular species:
    d17O({%FSIND}(:,0)) = d17O({%FSIND}(:,0)) / {%TAG}_ufac    ! de-permilizing (if pm units are used)
    d18O({%FSIND}(:,0)) = d18O({%FSIND}(:,0)) / {%TAG}_ufac

#ifndef ONLY_MINOR
    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * isofrac3a(d17O({%FSIND}(:,0)), Rstd_17O, d18O({%FSIND}(:,0)), Rstd_18O, {%NQATOM}({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * isofrac3r(d17O({%FSIND}(:,0)), Rstd_17O, d18O({%FSIND}(:,0)), Rstd_18O, {%NQATOM}({%FSIND}(:,0)))
    C({%FSIND}(:,3)) = C({%FSIND}(:,0)) * isofrac3r(d18O({%FSIND}(:,0)), Rstd_18O, d17O({%FSIND}(:,0)), Rstd_17O, {%NQATOM}({%FSIND}(:,0)))
#else
    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * isofrac3r(d17O({%FSIND}(:,0)), Rstd_17O, d18O({%FSIND}(:,0)), Rstd_18O, {%NQATOM}({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * isofrac3r(d18O({%FSIND}(:,0)), Rstd_18O, d17O({%FSIND}(:,0)), Rstd_17O, {%NQATOM}({%FSIND}(:,0)))
#endif
#endif
#ifdef UNIT_FRACMIN
  ! 16O, 17O, 18O through minor fractions and regular species:

    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * (1.0_dp - (d17O({%FSIND}(:,0)) + d18O({%FSIND}(:,0))))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * d17O({%FSIND}(:,0))
    C({%FSIND}(:,3)) = C({%FSIND}(:,0)) * d18O({%FSIND}(:,0))
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
#ifdef INIUNIT_PMC
  ! initialising using specified pMC values given in cfg

! {$f0} [%#%]  (%    PMC({%TAG}_@) = $%)

    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * isofrac2r_pMC(PMC({%FSIND}(:,0)), {%NQATOM}({%FSIND}(:,0)))
#else
  ! initialising using fractions given in cfg

! {$f0} [%#%]  (%    C({%RSIND}({%TAG}_@,#)) = C({%RSIND}({%TAG}_@,0)) * $%)
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
  ! initialising using fractions given in cfg

! {$f0} [%#%]  (%    C({%RSIND}({%TAG}_@,#)) = C({%RSIND}({%TAG}_@,0)) * $%)

-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}
  ! updating total {%ATOM} in the system

    call {%TAG}_calctotals(C)
    call {%TAG}_calcdeltas()

    T{%A}0(:) = T{%A}(:)

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
    d2HTH0 = d2HTH
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
    d13CTC0 = d13CTC
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
    d17OTO0 = d17OTO
    d18OTO0 = d18OTO
    DC17OTO0 = DC17OTO
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
    PMCT0 = PMCT
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}

  end subroutine {%TAG}_f0



! -----------------------------------------------------------------------------
  subroutine {%TAG}_emis(ind_r, amount, deltas)

    implicit none

    integer,  intent(in)    :: ind_r
    real(dp), intent(in)    :: amount
    real(dp), intent(in)    :: deltas(:)
    integer                 :: ind_t, ci

  ! getting tagging index
    call {%TAG}_ind_t(ind_r, ind_t)
    if (ind_t .lt. 1) return

! uncomment to manage emission of regular species via {%TAG} as well
!    C({%RSIND}(ind_t,0)) = C({%RSIND(ind_t,0)) + amount

->>- + isotopic part ++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
    ! 1H
    C({%RSIND}(ind_t,1)) = C({%RSIND}(ind_t,1)) + amount * &
      isofrac2a(deltas(1)/{%TAG}_ufac, Rstd_2H, {%NQATOM}(ind_t))
    ! 2H
    C({%RSIND}(ind_t,2)) = C({%RSIND}(ind_t,2)) + amount * &
      isofrac2r(deltas(1)/{%TAG}_ufac, Rstd_2H, {%NQATOM}(ind_t))
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
    ! 12C
    C({%RSIND}(ind_t,1)) = C({%RSIND}(ind_t,1)) + amount * &
      isofrac2a(deltas(1)/{%TAG}_ufac, Rstd_13C, {%NQATOM}(ind_t))
    ! 13C
    C({%RSIND}(ind_t,2)) = C({%RSIND}(ind_t,2)) + amount * &
      isofrac2r(deltas(1)/{%TAG}_ufac, Rstd_13C, {%NQATOM}(ind_t))
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
#ifndef ONLY_MINOR
    C({%RSIND}(ind_t,1)) = C({%RSIND}(ind_t,1)) + amount *  &
      isofrac3a(deltas(1)/{%TAG}_ufac, Rstd_17O, &
                deltas(2)/{%TAG}_ufac, Rstd_18O, {%NQATOM}(ind_t))
    C({%RSIND}(ind_t,2)) = C({%RSIND}(ind_t,2)) + amount *  &
      isofrac3r(deltas(1)/{%TAG}_ufac, Rstd_17O, &
                deltas(2)/{%TAG}_ufac, Rstd_18O, {%NQATOM}(ind_t))
    C({%RSIND}(ind_t,3)) = C({%RSIND}(ind_t,3)) + amount *  &
      isofrac3r(deltas(2)/{%TAG}_ufac, Rstd_18O, &
                deltas(1)/{%TAG}_ufac, Rstd_17O, {%NQATOM}(ind_t))
#else
    C({%RSIND}(ind_t,1)) = C({%RSIND}(ind_t,1)) + amount *  &
      isofrac3r(deltas(1)/{%TAG}_ufac, Rstd_17O, &
                deltas(2)/{%TAG}_ufac, Rstd_18O, {%NQATOM}(ind_t))
    C({%RSIND}(ind_t,2)) = C({%RSIND}(ind_t,2)) + amount *  &
      isofrac3r(deltas(2)/{%TAG}_ufac, Rstd_18O, &
                deltas(1)/{%TAG}_ufac, Rstd_17O, {%NQATOM}(ind_t))
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- + radiocarbon ++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
  ! emission with pMC (in deltas(1)) w.r.t. standard d13C = -19 per mil comp.
    C({%RSIND}(ind_t,1)) = C({%RSIND}(ind_t,1)) + amount * &
      isofrac2r_pMC(deltas(1), {%NQATOM}(ind_t))
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}
->>- + fractional tagging +++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
  ! deltas represent fractions here
    do ci = 1, {%NCLASS}
      C({%RSIND}(ind_t,ci)) = C({%RSIND}(ind_t,ci)) + amount * deltas(ci)
    end do
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}

  end subroutine {%TAG}_emis



! -----------------------------------------------------------------------------

  subroutine {%TAG}_depos(ind_r, factor)

    implicit none

    integer,  intent(in)    :: ind_r
    real(dp), intent(in)    :: factor
    integer                 :: ind_t

  ! getting tagging index
    call {%TAG}_ind_t(ind_r, ind_t)
    if (ind_t .lt. 1) return

  ! simple deposition routine, introduces no KIE during the deposition
    C({%RSIND}(ind_t,1:{%NISO})) = C({%RSIND}(ind_t,1:{%NISO})) * factor

  end subroutine {%TAG}_depos



! -----------------------------------------------------------------------------

  subroutine {%TAG}_pmix(TSL, dilF, ind_r, mix_amount, mix_deltas)

    implicit none

  ! pseudo-mixing of species ind_r with background concentration mix_amount of
  ! mix_deltas composition within TSL timestep with dilF dilution factor [1/s]

    real(dp), intent(in)    :: TSL, dilF      ! timestep length, dilution factor
    integer,  intent(in)    :: ind_r          ! reg. spec. index
    real(dp), intent(in)    :: mix_amount     ! backgr. concentration
    real(dp), intent(in)    :: mix_deltas(:)  ! backgr. deltas
    real(dp)                :: corr, tot
    integer                 :: ind_t

  ! getting tagging index
    call {%TAG}_ind_t(ind_r, ind_t)
    if (ind_t .lt. 1) return

  ! buget to correct to
    tot = sum(C({%RSIND}(ind_t,1:{%NCLASS})))
    corr = tot + ( mix_amount - tot ) * min( TSL * dilF, 1.0_dp )

  ! emission of background iso-composition
    call {%TAG}_emis(ind_r, mix_amount * TSL * dilF, mix_deltas)
    tot = sum(C({%RSIND}(ind_t,1:{%NCLASS})))

  ! removal preserving current composition
    if (tot .gt. 0.0_dp) then
      call {%TAG}_depos(ind_r, corr / tot)
    else
      C({%RSIND}(ind_t,1:{%NCLASS})) = 0.0_dp
      print *,'{%TAG}_pmix(',TSL,' ,',dilF,' ,',trim(SPC_NAMES(ind_r)),' ,', &
                 mix_amount,' ,',mix_deltas,'): mixing to nothing/negative'
    endif

  end subroutine {%TAG}_pmix



! -----------------------------------------------------------------------------

  subroutine {%TAG}_set(ind_r, init_amount, init_deltas)

    implicit none

  ! sets isotopic counterparts of the species ind_r
  ! with given amount and isotope composition
  ! beware: regular is not affected

    integer,  intent(in)    :: ind_r          ! reg. spec. index
    real(dp), intent(in)    :: init_amount     ! concentration
    real(dp), intent(in)    :: init_deltas(:)  ! deltas
    integer                 :: ind_t

  ! getting tagging index
    call {%TAG}_ind_t(ind_r, ind_t)
    if (ind_t .lt. 1) return

  ! zeroing
    C({%RSIND}(ind_t,1:{%NCLASS})) = 0.0_dp

  ! emitting required iso-composition
    call {%TAG}_emis(ind_r, init_amount, init_deltas)

  end subroutine {%TAG}_set



! -----------------------------------------------------------------------------

  subroutine {%TAG}_postprocess(skip_extra_calc)

    implicit none

  ! skip extra calculations (totals+delta), use only for optimisation purposes
    logical, optional, intent(in) :: skip_extra_calc

    integer  :: i
    real(dp) :: chkamnt

  ! calculating the number of specs falling below THRES
    {%TAG}_NREJCT = 0
    do i = 1, {%NSPEC}
      chkamnt = sum( C({%RSIND}(i,1:{%NISO})) )

      if (chkamnt .lt. THRES) then
        {%TAG}_NREJCT = {%TAG}_NREJCT + 1
!       C(RDCIND(i,0:{%NISO})) = 0.0_dp ! UNDEF
      endif
    enddo           ! ndspec cycle

    if ( present(skip_extra_calc) ) then
      if ( skip_extra_calc ) return
    endif

  ! every-step deltas/totals update
    call {%TAG}_calctotals(C)
    call {%TAG}_calcdeltas()

  end subroutine {%TAG}_postprocess



! -----------------------------------------------------------------------------

  subroutine {%TAG}_calcdeltas()

    implicit none

    integer  :: i
    real(dp) :: tot
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
    real(dp) :: f1H
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
    real(dp) :: f12C
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
    real(dp) :: f16O
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}

->>- + isotopic part ++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
  ! calculating new delta-2H values
    do i = 1, {%NSPEC}
      if (C({%RSIND}(i,1)) .gt. 0.0_dp) then
!      if ( ( C({%RSIND}(i,1))+C({%RSIND}(i,2)) ) .gt. THRES) then
        d2H(i) = delta2( C({%RSIND}(i,1)), C({%RSIND}(i,2)), &
                          Rstd_2H, {%NQATOM}(i) ) * {%TAG}_ufac
      else
        d2H(i) = UNDEF
      endif
    enddo        ! NISPEC-cycle

    ! total hydrogen
    if (T{%A}(1) /= 0.0_dp) then
      d2HTH = delta2( T{%A}(1), T{%A}(2), Rstd_2H, 1 ) * {%TAG}_ufac
    else
      d2HTH = UNDEF
    endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
  ! calculating new delta-13C values
    do i = 1, {%NSPEC}
      if (C({%RSIND}(i,1)) .gt. 0.0_dp) then
!      if ( ( C({%RSIND}(i,1))+C({%RSIND}(i,2)) ) .gt. THRES) then
        d13C(i) = delta2( C({%RSIND}(i,1)), C({%RSIND}(i,2)), &
                          Rstd_13C, {%NQATOM}(i) ) * {%TAG}_ufac
      else
        d13C(i) = UNDEF
      endif
    enddo        ! NISPEC-cycle

    ! total carbon
    if (T{%A}(1) /= 0.0_dp) then
      d13CTC = delta2( T{%A}(1), T{%A}(2), Rstd_13C, 1 ) * {%TAG}_ufac
    else
      d13CTC = UNDEF
    endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
  ! calculating new delta-17O, delta-18O, cap.delta-17O values
    do i = 1, {%NSPEC}
#ifndef ONLY_MINOR
      if (C({%RSIND}(i,1)) .gt. 0.0_dp) then
!      if ( ( C({%RSIND}(i,1))+C({%RSIND}(i,2))+C({%RSIND}(i,3)) ) .gt. THRES) then
        d17O(i)  = delta3( C({%RSIND}(i,1)), C({%RSIND}(i,2)), C({%RSIND}(i,3)), &
                           Rstd_17O, {%NQATOM}(i) )
        d18O(i)  = delta3( C({%RSIND}(i,1)), C({%RSIND}(i,3)), C({%RSIND}(i,2)), &
                           Rstd_18O, {%NQATOM}(i) )
#else
      tot = C({%RSIND}(i,0)) - ( C({%RSIND}(i,1)) + C({%RSIND}(i,2)) )
      if (tot .gt. 0.0_dp) then
        d17O(i)  = delta3( tot, C({%RSIND}(i,1)), C({%RSIND}(i,2)), &
                           Rstd_17O, {%NQATOM}(i) )
        d18O(i)  = delta3( tot, C({%RSIND}(i,2)), C({%RSIND}(i,1)), &
                           Rstd_18O, {%NQATOM}(i) )
#endif
        if ( d18O(i) .gt. -1.0_dp ) then
          DC17O(i) = (d17O(i)+1.0_dp)/(d18O(i)+1.0_dp)**MDFSL_O - 1.0_dp
          DC17O(i) = DC17O(i) * {%TAG}_ufac
        else
          DC17O(i) = UNDEF
        endif
        d17O(i) = d17O(i) * {%TAG}_ufac
        d18O(i) = d18O(i) * {%TAG}_ufac
      else
        d17O(i)  = UNDEF
        d18O(i)  = UNDEF
        DC17O(i) = UNDEF
      endif
    enddo        ! NISPEC-cycle

    ! total oxygen
#ifndef ONLY_MINOR
    if (T{%A}(1) /= 0.0_dp) then
      d17OTO  = delta3( T{%A}(1), T{%A}(2), T{%A}(3), Rstd_17O, 1 )
      d18OTO  = delta3( T{%A}(1), T{%A}(3), T{%A}(2), Rstd_18O, 1 )
#else
    tot = T{%A}(0) - ( T{%A}(2) + T{%A}(3) )
    if (tot /= 0.0_dp) then
      d17OTO  = delta3( tot, T{%A}(2), T{%A}(3), Rstd_17O, 1 )
      d18OTO  = delta3( tot, T{%A}(3), T{%A}(2), Rstd_18O, 1 )
#endif
      if ( d18OTO .gt. -1.0_dp ) then
        DC17OTO = (d17OTO+1.0_dp)/(d18OTO+1.0_dp)**MDFSL_O - 1.0_dp
        DC17OTO = DC17OTO * {%TAG}_ufac
      else
        DC17OTO = UNDEF
      endif
      d17OTO  = d17OTO * {%TAG}_ufac
      d18OTO  = d18OTO * {%TAG}_ufac
    else
      d17OTO  = UNDEF
      d18OTO  = UNDEF
      DC17OTO = UNDEF
    endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- + radiocarbon ++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
  ! calculating new pMC values
    PMC(:) = ratio_pMC( C({%RSIND}(:,0)), C({%RSIND}(:,1)), {%NQATOM}(:) )
  ! total pMC (scrambled C)
    PMCT = ratio_pMC( T{%A}(0), T{%A}(1), 1 )
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}
->>- + fractional tagging +++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
  ! calculating fractions here (!) w.r.t. original mech
    F(:,:) = UNDEF
#ifndef OPT_FTOT_WRTTAG
  ! total: careful, accounts for atom number in molecule
    tot = sum( C({%RSIND}(:,0)) * {%NQATOM}(:) )
#else
    tot = sum( sum(C({%RSIND}(:,1:{%NCLASS})),dim=2) * {%NQATOM}(:) )
#endif
    do i = 1, {%NCLASS}
      where (C({%RSIND}(:,0)) .NE. 0.0_dp)
        F(:,i) = C({%RSIND}(:,i)) / C({%RSIND}(:,0))
      endwhere
    ! fraction of total
      FT(i) = sum( C({%RSIND}(:,i)) * {%NQATOM}(:) ) / tot
    enddo
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}

  end subroutine {%TAG}_calcdeltas



! -----------------------------------------------------------------------------

! output file for tagged species info
  subroutine {%TAG}_init()

    implicit none

  ! print configuration info / parameters
    call {%TAG}_info()

  ! set relative tolerances by default equial to those of tagged species
    call {%TAG}_set_atol_rtol()

! TODO: put additional tracers/variables+units after INIT_TRAC, INIT_UNIT

    call open_output_file(ncid_{%TAG}, 'caaba_mecca_{%TAG}', &
      (/   &
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
! {$TAG_SPECS} [%I1@%]
       , &
! {$TAG_SPECS} [%I2@%]
       , &
! {$TAG_SPECS} [%d2@%]
       , &
{$ELSA}       'TH_R', 'I1TH', 'I2TH', 'd2TH' &
       , &
{$ELSA}       'd0TH_R', 'd0TH', 'd0D2TH' &
       , &
{$ELSA}       'NREJCT' &
       /), (/   &
! {$TAG_SPECS} [%mol/mol%]
       , &
! {$TAG_SPECS} [%mol/mol%]
       , &
! {$TAG_SPECS} [%o/oo%]
       , &
{$ELSA}       'atoms', 'atoms', 'atoms', 'o/oo' &
       , &
{$ELSA}       'atoms', 'atoms', 'o/oo' &
       , &
{$ELSA}       'specs' &
       /), (/   &
! {$TAG_SPECS} [%\@SR^1@%]
       , &
! {$TAG_SPECS} [%\@SR^2@%]
       , &
! {$TAG_SPECS} [%\@SGd\@SRD(@)%]
       , &
{$ELSA}       '@SRTH_R (regular mech)', '@SRT^1^2C', '@SRT^1^3C', '@SGdD@SR(TH)' &
       , &
{$ELSA}       '@SGD@SR_t_0(TH_R)', '@SGD@SR_t_0(TH_D)', '@SGD@SR_t_0(@SGdD@SR(TH_D)) ' &
       , &
{$ELSA}       '@SRnumber of rejected species' &
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
! {$TAG_SPECS} [%I12@%]
       , &
! {$TAG_SPECS} [%I13@%]
       , &
! {$TAG_SPECS} [%d13@%]
       , &
{$ELSA}       'TC_R', 'I12TC', 'I13TC', 'd13TC' &
       , &
{$ELSA}       'd0TC_R', 'd0TC', 'd0d13TC' &
       , &
{$ELSA}       'NREJCT' &
       /), (/   &
! {$TAG_SPECS} [%mol/mol%]
       , &
! {$TAG_SPECS} [%mol/mol%]
       , &
! {$TAG_SPECS} [%o/oo%]
       , &
{$ELSA}       'atoms', 'atoms', 'atoms', 'o/oo' &
       , &
{$ELSA}       'atoms', 'atoms', 'o/oo' &
       , &
{$ELSA}       'specs' &
       /), (/   &
! {$TAG_SPECS} [%\@SR^1^2@%]
       , &
! {$TAG_SPECS} [%\@SR^1^3@%]
       , &
! {$TAG_SPECS} [%\@SGd\@SR^1^3C(@)%]
       , &
{$ELSA}       '@SRTC_R (regular mech)', '@SRT^1^2C', '@SRT^1^3C', '@SGd@SR^1^3C(TC)' &
       , &
{$ELSA}       '@SGD@SR_t_0(TC_R)', '@SGD@SR_t_0(TC_D)', '@SGD@SR_t_0(@SGd@SR^1^3C(TC_D)) ' &
       , &
{$ELSA}       '@SRnumber of rejected species' &
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
! {$TAG_SPECS} [%I16@%]
       , &
! {$TAG_SPECS} [%I17@%]
       , &
! {$TAG_SPECS} [%I18@%]
       , &
! {$TAG_SPECS} [%d18@%]
       , &
! {$TAG_SPECS} [%d17@%]
       , &
! {$TAG_SPECS} [%DC17@%]
       , &
{$ELSA}       'TO_R', 'I16TO', 'I17TO', 'I18TO', &
{$ELSA}       'd18TO', 'd17TO', 'DC17TO' &
       , &
{$ELSA}       'd0TO_R', 'd0TO', &
{$ELSA}       'd0d18TO', 'd0D17TO', 'd0DC17TO' &
       , &
{$ELSA}       'NREJCT' &
       /), (/   &
! {$TAG_SPECS} [%mol/mol%]
       , &
! {$TAG_SPECS} [%mol/mol%]
       , &
! {$TAG_SPECS} [%mol/mol%]
       , &
! {$TAG_SPECS} [%o/oo%]
       , &
! {$TAG_SPECS} [%o/oo%]
       , &
! {$TAG_SPECS} [%o/oo%]
       , &
{$ELSA}       'atoms', 'atoms', 'atoms', 'atoms', &
{$ELSA}       'o/oo', 'o/oo', 'o/oo' &
       , &
{$ELSA}       'atoms', 'atoms', &
{$ELSA}       'o/oo', 'o/oo', 'o/oo' &
       , &
{$ELSA}       'specs'  &
       /), (/   &
! {$TAG_SPECS} [%\@SR^1^6@%]
       , &
! {$TAG_SPECS} [%\@SR^1^7@%]
       , &
! {$TAG_SPECS} [%\@SR^1^8@%]
       , &
! {$TAG_SPECS} [%\@SGd\@SR^1^8O(@)%]
       , &
! {$TAG_SPECS} [%\@SGd\@SR^1^7O(@)%]
       , &
! {$TAG_SPECS} [%\@SGD\@SR^1^7O(@)%]
       , &
{$ELSA}       '@SRTO_R', '@SRT^1^6O', '@SRT^1^7O', '@SRT^1^8O', &
{$ELSA}       '@SGd@SR^1^8O(to)', '@SGd@SR^1^7O(to)', '@SGD@SR^1^7O(to)' &
       , &
{$ELSA}       '@SGD@SR_t_0(TO_R)', '@SGD@SR_t_0(TO_D)', &
{$ELSA}       '@SGD@SR_t_0(@SGd@SR^1^8O(TO_D))', '@SGD@SR_t_0(@SGd@SR^1^7O(TO_D))', '@SGD@SR_t_0(@SGD@SR^1^7O(TO_D))' &
       , &
{$ELSA}       '@SRnumber of rejected species'  &
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- + radiocarbon ++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
! {$TAG_SPECS} [%I14@%]
       , &
! {$TAG_SPECS} [%pMC_@%]
       , &
{$ELSA}       'TC_R', 'I14TC', 'pMCT', 'pMCT0' &
       , &
{$ELSA}       'NREJCT' &
       /), (/   &
! {$TAG_SPECS} [%mol/mol%]
       , &
! {$TAG_SPECS} [%pMC%]
       , &
{$ELSA}       'atoms', 'atoms', 'pMC', 'pMC' &
       , &
{$ELSA}       'specs' &
       /), (/   &
! {$TAG_SPECS} [%\@SR^1^4@%]
       , &
! {$TAG_SPECS} [%\@SR^1^4^C^/^1^2^CR(@)%]
       , &
{$ELSA}       '@SRTC_R (regular mech)', '@SRT^1^4C', '@SR^1^4^C^/^1^2^CR(TC)', '@SR^1^4^C^/^1^2^CR_t_0(TC)' &
       , &
{$ELSA}       '@SRnumber of rejected species' &
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}
->>- + fractional tagging +++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
{$ELSA}       'T_R', 'T_1' &
       , &
{$ELSA}       'd0T_R', 'd0T_1' &
       , &
{$ELSA}       'NREJCT' &
       , &
! {$TAG_SPECS} [%f$_@%]
       /), (/   &
{$ELSA}       'mol/mol', 'mol/mol' &
       , &
{$ELSA}       'mol/mol', 'mol/mol' &
       , &
{$ELSA}       'steps' &
       , &
! {$TAG_SPECS} [%$ share%]
       /), (/   &
{$ELSA}       '@SRT_R (regular mech)', '@SRT_1' &
       , &
{$ELSA}       '@SGD@SR_t_0(T_R)', '@SGD@SR_t_0(T_1)' &
       , &
{$ELSA}       '@SRnumber of rejected species' &
       , &
! {$TAG_SPECS} [%\@SRf($) in @%]
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}
       /) )

  end subroutine {%TAG}_init



! -----------------------------------------------------------------------------

  subroutine {%TAG}_result(model_time)

    implicit none

    real(dp), intent(in) :: model_time
    integer              :: i, j

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
    do i = 1, {%NSPEC}
      D{%A}out(i)            = C({%RSIND}(i,1))/cair
      D{%A}out({%NSPEC}+i)   = C({%RSIND}(i,2))/cair
      D{%A}out({%NSPEC}*2+i) = d2H(i)
    enddo

  ! totals
    D{%A}out({%NSPEC}*3+1) = T{%A}(0)
    D{%A}out({%NSPEC}*3+2) = T{%A}(1)
    D{%A}out({%NSPEC}*3+3) = T{%A}(2)
    D{%A}out({%NSPEC}*3+4) = d2HTH

  ! totals verification
    D{%A}out({%NSPEC}*3+5) = T{%A}(0)-T{%A}(0)                ! d0TH_R  = TH_R - TH_R(t=0)
    D{%A}out({%NSPEC}*3+6) = sum(T{%A}(1:2))-TH0(0)           ! d0TH    = (T1H+T2H) - TH_R(t=0)
    D{%A}out({%NSPEC}*3+7) = (d2HTH-d2HTH0)                   ! d0d2HTH = d2HTH - d2HTH(t=0)
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
    do i = 1, {%NSPEC}
      D{%A}out(i)           = C({%RSIND}(i,1))/cair
      D{%A}out({%NSPEC}+i)   = C({%RSIND}(i,2))/cair
      D{%A}out({%NSPEC}*2+i) = d13C(i)
    enddo

  ! totals
    D{%A}out({%NSPEC}*3+1) = T{%A}(0)
    D{%A}out({%NSPEC}*3+2) = T{%A}(1)
    D{%A}out({%NSPEC}*3+3) = T{%A}(2)
    D{%A}out({%NSPEC}*3+4) = d13CTC

  ! totals verification
    D{%A}out({%NSPEC}*3+5) = T{%A}(0)-TC0(0)                  ! d0TC_R   = TC_R - TC_R(t=0)
    D{%A}out({%NSPEC}*3+6) = sum(T{%A}(1:2))-TC0(0)           ! d0TC     = (T12C+T13C) - TC_R(t=0)
    D{%A}out({%NSPEC}*3+7) = (d13CTC-d13CTC0)                 ! d0d13CTC = d13CTC - d13CTC(t=0)
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
    do i = 1, {%NSPEC}
#ifndef ONLY_MINOR
      D{%A}out(i)            = C({%RSIND}(i,1))/cair
      D{%A}out({%NSPEC}+i)   = C({%RSIND}(i,2))/cair
      D{%A}out({%NSPEC}*2+i) = C({%RSIND}(i,3))/cair
#else
      D{%A}out(i)           = (C({%RSIND}(i,0))-( C({%RSIND}(i,1)) + C({%RSIND}(i,2)) ))/cair
      D{%A}out({%NSPEC}+i)   = C({%RSIND}(i,1))/cair
      D{%A}out({%NSPEC}*2+i) = C({%RSIND}(i,2))/cair
#endif
      D{%A}out({%NSPEC}*3+i) = d18O(i)
      D{%A}out({%NSPEC}*4+i) = d17O(i)
      D{%A}out({%NSPEC}*5+i) = DC17O(i)
    enddo

  ! totals
    D{%A}out({%NSPEC}*6+1) = T{%A}(0)
#ifndef ONLY_MINOR
    D{%A}out({%NSPEC}*6+2) = T{%A}(1)
#else
    D{%A}out({%NSPEC}*6+2) = T{%A}(0) - sum(T{%A}(2:3))
#endif
    D{%A}out({%NSPEC}*6+3) = T{%A}(2)
    D{%A}out({%NSPEC}*6+4) = T{%A}(3)
    D{%A}out({%NSPEC}*6+5) = d18OTO
    D{%A}out({%NSPEC}*6+6) = d17OTO
    D{%A}out({%NSPEC}*6+7) = DC17OTO

  ! totals verification
    D{%A}out({%NSPEC}*6+8) = T{%A}(0)-TO0(0)        ! d0TO_R    = TO_R - TO_R(t=0)
#ifndef ONLY_MINOR
    D{%A}out({%NSPEC}*6+9) = sum(T{%A}(1:3))-TO0(0) ! d0TO = (T16O+T17O+T18O) - TO_R(t=0)
#else
    D{%A}out({%NSPEC}*6+9) = D{%A}out({%NSPEC}*6+2)-TO0  ! d0TO = (T16O+T17O+T18O) - TO_R(t=0)
#endif
    D{%A}out({%NSPEC}*6+10) = (d18OTO-d18OTO0)     ! d0d18OTO  = d18OTO - d18OTO(t=0)
    D{%A}out({%NSPEC}*6+11) = (d17OTO-d17OTO0)     ! d0d18OTO  = d17OTO - d17OTO(t=0)
    D{%A}out({%NSPEC}*6+12) = (DC17OTO-DC17OTO0)   ! d0DC17OTO = DC17OTO - DC18OTO(t=0)

  ! last value is NREJCT
    D{%A}out(UBOUND(D{%A}out)) = real({%TAG}_NREJCT)
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC}
    do i = 1, {%NSPEC}
      D{%A}out({%NSPEC}*1+i) = C({%RSIND}(i,1))/cair
      D{%A}out({%NSPEC}*1+i) = PMC(i)
    enddo

  ! totals
    D{%A}out({%NSPEC}*2+1) = T{%A}(0)
    D{%A}out({%NSPEC}*2+2) = T{%A}(1)
    D{%A}out({%NSPEC}*2+3) = PMCT
  ! totals verification
    D{%A}out({%NSPEC}*2+4) = (T{%A}(0)-T{%A}0(0))
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
  ! totals
    D{%A}out(1) = T{%A}(0)
    D{%A}out(2) = T{%A}(1)
  ! totals verification
    D{%A}out(3) = T{%A}(0)-T{%A}0(0)
    D{%A}out(4) = T{%A}(1)-T{%A}0(0)
  ! NREJCT
    D{%A}out(5) = real({%TAG}_NREJCT)
  ! fractions
    do j = 1, {%NCLASS}
      do i = 1, {%NSPEC}
        D{%A}out(5+{%NSPEC}*(j-1)+i) = F(i,j)
      enddo
    enddo
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}

    call write_output_file(ncid_{%TAG}, model_time, D{%A}out)

  end subroutine {%TAG}_result



! -----------------------------------------------------------------------------

  subroutine {%TAG}_finish()

    call close_file(ncid_{%TAG})

  end subroutine {%TAG}_finish



! - some cfg cheks ------------------------------------------------------------

#ifndef UNIT_DELTAPERMIL
#ifndef UNIT_FRACMIN
#ifndef ZERO_TEST
#ifndef NULL_TEST
 FATAL: (init)units are not defined, check the parameters
#endif
#endif
#endif
#endif

! -----------------------------------------------------------------------------

end module messy_mecca_{%TAG}_box

! *****************************************************************************

