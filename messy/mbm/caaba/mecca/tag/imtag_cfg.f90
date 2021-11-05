! ==============================================================================
! {%TAG} core
! generated: {%TIMEDATE}
!
! this module is created automatically by imtag tool - do not edit
!
! maintenance routines for budgeting configurations (isotopes)
! level: smcl
!
! {$TAG_INFO} ! this is a template file for imtag utility
!
! [Gromov, MPIC, 2007-2020]
! ==============================================================================

! - general tagging parameters (as conditional defines) -----------------------

#include "{%CMODEL}_tag_parameters.inc"

! ------------------------------------------------------------------------------

! {$CONF_PARAM}

module {%CMODEL}_{%TAG}

  use messy_mecca_kpp     ! dp, ... nreact, nspec, ind_*, SPC_NAMES, EQN_TAGS
  use {%CMODEL}_tag_common

  implicit none

! treshold value: below it, species might stop to sink to the others
! (but can receive still)
  real(dp), parameter :: THRES = 1.0E-40_dp * 2.5047E+19_dp
!                                ?          * mean cair

! isotope standards, scales, units
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
! reference standard for stable H isotopes
! VSMOW scale
  real(dp), parameter :: Rstd_2H = VSMOW_2H
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
! reference standard for stable C isotopes
! VPDB scale
  real(dp), parameter :: Rstd_13C = VPDB_13C
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
! reference standard for stable O isotopes
#ifdef Rstd_VPDBCO2
! VPDB-CO2 scale
  real(dp), parameter :: Rstd_17O = VPDB_17O
  real(dp), parameter :: Rstd_18O = VPDB_18O
#else
! VSMOW scale
  real(dp), parameter :: Rstd_17O = VSMOW_18O
  real(dp), parameter :: Rstd_18O = VSMOW_18O
#endif
! (mass-dependent) fractionation slope
  real(dp), parameter :: MDFSL_O = MDFSL_MWL
 !real(dp), parameter :: MDFSL_O = MDFSL_LVE
 !real(dp), parameter :: MDFSL_O = MDFSL_CO2
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
! unit factor
  real(dp), parameter :: {%TAG}_ufac = &
#ifdef UNIT_DELTAPERMIL
    1e3_dp
#elif UNIT_DELTAPERMEG
    1e6_dp
#else
    1.0_dp
#endif

! ------------------------------------------------------------------------------

! here constants and tagged species indices are to be defined
! {$TRAC_DECL} [%{%TAG}_@%]

! ------------------------------------------------------------------------------

! total (scrambled) burdens of regular and isotope-distinguished atoms
  real(dp)            :: T{%A}(0:{%NCLASS})

! no. of "rejected" species (below given threshold)
  integer             :: {%TAG}_NREJCT

! ------------------------------------------------------------------------------

  private NKRSPEC, KRSIND

  public {%TAG}_ind_t
  public {%TAG}_set_atol_rtol
  public {%TAG}_calctotals
  public {%TAG}_correct
  public {%TAG}_correct2reg
  public {%TAG}_correct2tag
  public {%TAG}_resetPTs
  public {%TAG}_info

! ==============================================================================

contains

! ==============================================================================

!>\brief
!! returns corresponding tagging index ind_t for a regular ind_r

  subroutine {%TAG}_ind_t(ind_r, ind_t)


    implicit none

    integer, intent(in)  :: ind_r
    integer, intent(out) :: ind_t
    integer              :: i

    ind_t = -1
    do i = 1, {%NSPEC}
      if ({%RSIND}(i,0) .eq. ind_r) then
        ind_t = i
        return
      endif
    enddo

  end subroutine {%TAG}_ind_t



! -----------------------------------------------------------------------------
!> \brief
!! set absolute and reltive tolerance for tagged species
  subroutine {%TAG}_set_atol_rtol(ind_r, atol_in, rtol_in)

    use messy_mecca_kpp, only: atol, rtol, SPC_NAMES

    implicit none

  ! indices
    integer, intent(in), optional :: ind_r
    integer                       :: ind_t, i

  ! make sure you use proper array sizes (i.e. 1:{%NCLASS})
    real(dp), intent(in), optional :: atol_in(:), rtol_in(:)

  ! if no index specified, adjust all species according
  ! to tolerances of regular species
    if (.not.present(ind_r)) then
      do i = 1, {%NCLASS}
        where (atol({%RSIND}(:,0)).ne.0._dp)
          atol({%RSIND}(:,i)) = atol({%RSIND}(:,0))
        endwhere
        where (rtol({%RSIND}(:,0)).ne.0._dp)
          rtol({%RSIND}(:,i)) = rtol({%RSIND}(:,0))
        endwhere
      enddo

#ifdef DEBUG
      print *,"{%TAG}_set_atol_rtol(): tolerances"
      do i = 1,{%NSPEC}
        print '(a,*(e7.1E1))', &
                trim(SPC_NAMES({%RSIND}(i,0)))//' :', &
                atol({%RSIND}(i,:)), rtol({%RSIND}(i,:))
      enddo
#endif
      return
    endif

  ! get the tagging index
    call {%TAG}_ind_t(ind_r, ind_t)
    if (ind_t.le.0) return

    if (present(atol_in)) then
    ! use supplied tolerance value
        atol({%RSIND}(ind_t,1:{%NCLASS})) = atol_in
    else
    ! by default, set tolerance equal to that of regular species
      if (atol({%RSIND}(ind_t,0) ).ne.0._dp) &
          atol({%RSIND}(ind_t,1:{%NCLASS})) = atol({%RSIND}(ind_t,0))
    endif

    if (present(rtol_in)) then
    ! use supplied tolerance value
      rtol({%RSIND}(ind_t,1:{%NCLASS})) = rtol_in
    else
    ! by default, set tolerance equal to that of regular species
      if (rtol({%RSIND}(ind_t,0)).ne.0._dp) &
          rtol({%RSIND}(ind_t,1:{%NCLASS})) = rtol({%RSIND}(ind_t,0))
    endif

#ifdef DEBUG
    print '(a,*(e7.1E1))', &
      "{%TAG}_set_atol_rtol(): adjusted tolerances for "//trim(SPC_NAMES(ind_r))//' :',atol({%RSIND}(i,:)),rtol({%RSIND}(i,:))
#endif

  end subroutine {%TAG}_set_atol_rtol



! -----------------------------------------------------------------------------
!> \brief
!! calculation of the total number of {%ATOM} atoms, according to each species' composition

  subroutine {%TAG}_calctotals(C)

    implicit none

    real(dp), intent(inout) :: C(:)     !< operational vector of concentrations
    integer  :: i


->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
  ! isotopic tagging

  ! TC is calculated for regular mech
    T{%A}(0) = sum( C({%RSIND}(:,0))*real({%NQATOM}(:),dp) )

  ! major isotopologues of abundant isotope
    T{%A}(1) = sum( C({%RSIND}(:,1))*real({%NQATOM}(:),dp) )

    do i = 2, {%NISO}

    ! rare isotope of minor isotopologues
      T{%A}(i) = sum( C({%RSIND}(:,i)) )

    ! adding rare isotope of minor to total of major
      T{%A}(1) = T{%A}(1) + &
        sum( C({%RSIND}(:,i))*real({%NQATOM}(:)-1,dp) )

    enddo
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:RC.+}
  ! radiocarbon
  ! total 14C
    T{%A}(1) = sum( C({%RSIND}(:,1)) )
  ! total stable C
    T{%A}(0) = sum( C({%RSIND}(:,0))*real({%NQATOM}(:),dp) ) + &
               sum( C({%RSIND}(:,1))*real({%NQATOM}(:)-1,dp) )
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:RC.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
  ! fractional tagging

  ! total for regular
    T{%A}(0) = sum( C({%RSIND}(:,0))*real({%NQATOM}(:),dp) )

  ! totals for tagged
    do i = 1, {%NCLASS}
      T{%A}(i) = sum( C({%RSIND}(:,i))*real({%NQATOM}(:),dp) )
    enddo
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}

#ifdef DEBUG
    print *,'{%TAG}_calctotals: passed'
#endif
#ifdef DEEPDEBUG
    print *,'{%TAG}_calctotals: T{%ATOM} (R+TAG): ',T{%A}(:))
#endif

  end subroutine {%TAG}_calctotals



! -----------------------------------------------------------------------------
!>\brief
!! correction of total isotopologues budget to "regular" species budget

  subroutine {%TAG}_correct(C)

    implicit none

    real(dp), intent(inout) :: C(:)     ! operational vector of concentrations
    integer  :: i
    real(dp) :: total


#ifdef OPT_NEG_FILTER
  ! in case anything shoots to negative, correcting
    do i = 1, {%NCLASS}
      where ( C({%RSIND}(1:,i)) .lt. 0.0_dp )
        C({%RSIND}(:,i)) = 0.0_dp
      endwhere
    enddo
#endif

#ifdef OPT_EXC_FILTER
  ! in case fractional classes are defined, checking overshooting with
  ! relation to the original species only
    do i = 1, {%NCLASS}
      where ( C({%RSIND}(:,i)) .gt. C({%RSIND}(:,0)) )
        C({%RSIND}(:,i)) = C({%RSIND}(:,0))
      endwhere
    enddo
#endif

#ifdef OPT_NO_CORR
  ! no correction option
    return
#endif

#ifdef OPT_CORR_2TAG
  ! correction with tagged mech as reference
    call {%TAG}_correct2tag(C)
#else
  ! normal correction to the regular mech
    call {%TAG}_correct2reg(C)
#endif

  end subroutine {%TAG}_correct

! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
!>\brief
!! correction of total isotopologues budget to "regular" species budget

  subroutine {%TAG}_correct2reg(C)

    implicit none

    real(dp), intent(inout) :: C(:)     ! operational vector of concentrations
    integer  :: i
    real(dp) :: total

#ifdef CLASSES_1

#ifdef OPT_OSHOOT_FILTER
  ! in case one class is defined, checking overshooting only
    where ( C({%RSIND}(:,1)) .gt. C({%RSIND}(:,0)) )
      C({%RSIND}(:,1)) = C({%RSIND}(:,0))
    endwhere
#endif

#ifdef DEBUG
    print *,'{%TAG}_correct2reg: overshoot correction performed (one class)'
#endif
    return
#endif

#ifdef ONLY_MINOR
#ifdef DEBUG
    print *,'{%TAG}_correct2reg: skipped (only minor isotopologues)'
#endif
    return
#endif

#ifdef OPT_USE_KRSIND
  ! here is the ver. with corr. of only KIE-rel species to regular
  ! correcting only species related to KIE in this meccanism
    do i = 1, NKRSPEC
      total = sum ( C({%RSIND}(KRSIND(i),1:{%NISO})) )
      if (total .gt. 0.0_dp) then
        C({%RSIND}(KRSIND(i),1:{%NISO})) = &
          ( C({%RSIND}(KRSIND(i),1:{%NISO})) * &
            C({%RSIND}(KRSIND(i),0)) ) / total
      else
        C({%RSIND}(KRSIND(i),1:{%NISO})) = 0.0_dp
      endif
    enddo
#else
  ! here is the ver. with corr. of ALL species to regular
    do i = 1, {%NSPEC}
      total = sum( C({%RSIND}(i,1:{%NISO})) )
      if (total .gt. 0.0_dp) then
        C({%RSIND}(i,1:{%NISO})) = ( C({%RSIND}(i,1:{%NISO})) * C({%RSIND}(i,0)) ) / total
      else
        C({%RSIND}(i,1:{%NISO})) = 0.0_dp
      endif
    enddo
#endif

#ifdef DEEPDEBUG
    print *,'{%TAG}_correct2reg: passed'
#endif

  end subroutine {%TAG}_correct2reg



! -----------------------------------------------------------------------------
!>\brief
!! correction of "regular" species budget to the total isotopologues budget

  subroutine {%TAG}_correct2tag(C)

    implicit none

    real(dp), intent(inout) :: C(:)     ! operational vector of concentrations
    integer                 :: i

#ifdef CLASSES_1
  ! in case one class is defined, quitting
#ifdef DEBUG
    print *,'{%TAG}_correct2tag: no correction performed (one class)'
#endif
    return
#endif

#ifdef OPT_USE_KRSIND
  ! here is the ver. with corr. of only KIE-rel species to regular
    do i = 1, NKRSPEC
      C({%RSIND}(KRSIND(i),0)) = sum( C({%RSIND}(KRSIND(i),1:{%NISO})) )
    enddo
#else
  ! here is the ver. with corr. of ALL species to regular
    do i = 1, {%NSPEC}
      C({%RSIND}(i,0)) = sum( C({%RSIND}(i,1:{%NISO})) )
    enddo
#endif

#ifdef DEEPDEBUG
    print *,'{%TAG}_correct2tag: passed'
#endif

  end subroutine {%TAG}_correct2tag



! -----------------------------------------------------------------------------
!>\brief
!! passive tracers initialization (reset) routine

  subroutine {%TAG}_resetPTs(C)


    implicit none

    real(dp), intent(inout) :: C(:)  !< operational vector of concentrations

! {x$RESET_PTs}
! - currently disabled with use of DRPT{%ATOM}IND()

#ifdef USE_PT
    C(DRPT{%ATOM}IND(:)) = 0.0_dp    ! <-- boxmodel syntax
#endif

#ifdef DEEPDEBUG
    print *,'{%TAG}_resetPTs: passed'
#endif

  end subroutine {%TAG}_resetPTs



! -----------------------------------------------------------------------------
!>\brief
!! issue configuration info / warnings

  subroutine {%TAG}_info()


    implicit none

    character(len=:), allocatable :: si
    character, parameter          :: cr = char(10)

    si = ''

#ifdef OPT_EXC_FILTER
    si = si//"  OPT_EXC_FILTER  : check & correct overshooting for fractional tagging"//cr
#endif
#ifdef OPT_NEG_FILTER
    si = si//"  OPT_NEG_FILTER  : check & correct negative overshooting"//cr
#endif
#ifdef OPT_LOW_FILTER
    si = si//"  OPT_LOW_FILTER  : tagging - treshold cutoff optimisation"//cr
#endif
#ifdef OPT_NO_CORR
    si = si//"  OPT_NO_CORR     : switch off regular <-> tagged mechs correction"//cr
#endif
#ifdef OPT_C2R_FILTER
    si = si//"  OPT_C2R_FILTER  : filter only largely deviated species in correct2reg"//cr
#endif
#ifdef OPT_CORR_2TAG
    si = si//"  OPT_CORR_2TAG   : (if) correction is done (then) with tagged mech as a reference"//cr
#endif
#ifdef OPT_FTOT_WRTTAG
    si = si//"  OPT_FTOT_WRTTAG : calculate fractions of totals w.r.t. to the tagged mech (default is regular)"//cr
#endif
#ifdef OPT_USE_KRSIND
    si = si//"  OPT_USE_KRSIND  : use kie-relates species indices for correction"//cr
#endif
#ifdef TAG_OPT_SKIP_XTRABOXPPROC
    si = si//"  TAG_OPT_SKIP_XTRABOXPPROC : by def. skip extra calc. (delta+totals) in the postprocessing (box-model opt.)"//cr
#endif

#ifdef DEBUG
    si = si//"  DEBUG           : debug enabled"
#endif
#ifdef DEEPDEBUG
    si = si//"  DEEPDEBUG       : deep debug with lots of repetitive messages enabled"//cr
#endif
#ifdef ZERO_TEST
    si = si//"  ZERO_TEST       : rare isotopologues are initialized with 0 permil rel. to standard"//cr
#endif
#ifdef NULL_TEST
    si = si//"  NULL_TEST       : rare isotopologues / fractions are initialized emptied"//cr
#endif
#ifdef INITFAST_E5
    si = si//"  INITFAST_E5     : fast first-step E5M1/EMAC initialization of configurations"//cr
#endif

    if ( si.ne.'' ) print '(a)', "{%CMODEL}_{%TAG}: warning, following options are enabled:"//cr//si

  end subroutine {%TAG}_info

! -----------------------------------------------------------------------------


end module messy_mecca_{%TAG}

! *****************************************************************************

