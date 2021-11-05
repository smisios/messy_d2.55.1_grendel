! ==============================================================================
! {%TAG}_box
! generated: {%TIMEDATE}
!
! this module is created automatically by imtag tool - do not edit
!
! maintenance routines for budgeting configurations (fractional tagging)
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

  use messy_mecca_kpp     ! dp, nreact, nspec, ind_*, SPC_NAMES, EQN_TAGS
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

! treshold value: below it, species might stop to sink to the others
! (but can receive still)
  real(dp), parameter :: THRES = 1.0E-40_dp * 2.5047E+19_dp
!                                ?          * mean cair

! -----------------------------------------------------------------------------

! here constants and tagged species indices are to be defined
! {$TRAC_DECL} [%ind_@%] <-- boxmodel syntax  (%{%TAG}_@%) <-- isotracers syntax

! -----------------------------------------------------------------------------

! no. of "rejected" species (under threshold)
  integer            :: {%TAG}_NREJCT

! output array: minor fractions (+total's fractions),
!               total concentrations (+regular) + NREJCT
  real(dp)           :: TOUT(({%NSPEC}+1)*({%NCLASS}-1)+ &
                             {%NCLASS}+1+1)

! classes fractions
  real(dp)           :: CF({%QSPEC},{%NCLASS})

! totals: concentration & budget fractions
  real(dp)           :: TCC({%NCLASS}), TCF({%NCLASS})

! -----------------------------------------------------------------------------

  public {%TAG}_x0
  public {%TAG}_emis
  public {%TAG}_depos
  public {%TAG}_pmix
  public {%TAG}_process
  public {%TAG}_calctotals
  public {%TAG}_calcfractions
  public {%TAG}_correct2reg
  public {%TAG}_correct2tag
! public {%TAG}_fudge
  public {%TAG}_resetPTs
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

#ifndef UNIT_FRACMIN
 FATAL: initialization unit is not fracmin, please check configuration and former
#endif

! {$x0} [%#%] (%    CF({%TAG}_@,#) = $%)

#ifdef ZERO_TEST
    CF(:,1) = 1.0_dp
    do i = 2, {%NCLASS}
      CF(:,i) = 0.0_dp
    enddo
#endif

  ! initializing isotopologues concentration according to "regular", then

  ! setting all classes fractions
    do i = 1, {%NDSPEC}
      CF(i,1) = 1.0_dp - sum(CF(i,2:{%NCLASS}))
    enddo

  ! initializing tagged tracers according to the fractions
    do i = 1, {%NCLASS}
      C({%RDIND}(:,i)) = C({%RDIND}(:,0)) * CF(:,i)
    enddo

  ! updating totals in the system
    call {%TAG}_calctotals
    call {%TAG}_calcfractions

  end subroutine {%TAG}_x0



! -----------------------------------------------------------------------------

  subroutine {%TAG}_emis(ind_d, amount, fracs)

    implicit none

    integer,  intent(in)    :: ind_d
    real(dp), intent(in)    :: amount
    real(dp), intent(in)    :: fracs(:)

  ! filtering possible dummies
    if ((ind_d .lt. 1) .or. (ind_d .gt. {%NDSPEC})) return

! uncomment to manage emission only through {%TAG}
!    C({%RDIND}(ind_d,0)) = C({%RDIND(ind_d,0)) + amount

! emission of corresponding amount fractions into the box
    C({%RDIND}(ind_d,1:{%NCLASS})) = C({%RDIND}(ind_d,1:{%NCLASS})) + amount * fracs(:)

  end subroutine {%TAG}_emis



! -----------------------------------------------------------------------------

  subroutine {%TAG}_depos(ind_d, factor)

    implicit none

    integer,  intent(in)    :: ind_d
    real(dp), intent(in)    :: factor

  ! filtering possible dummies
    if ((ind_d .lt. 1) .or. (ind_d .gt. {%NDSPEC})) return

  ! simple deposition routine, introduces no selective deposition

    C({%RDIND}(ind_d,1:{%NCLASS})) = C({%RDIND}(ind_d,1:{%NCLASS})) * factor

  end subroutine {%TAG}_depos



! -----------------------------------------------------------------------------

  subroutine {%TAG}_pmix(TSL, dilF, ind_d, mix_amount, mix_fracs)

    implicit none

  ! pseudo-mixing of species ind_d with background concentration mix_amount of
  ! mix_deltas composition within TSL timestep with dilF dilution factor [1/s]

    real(dp), intent(in)    :: TSL, dilF      ! timestep length, dilution factor
    integer,  intent(in)    :: ind_d          ! spec. index
    real(dp), intent(in)    :: mix_amount     ! backgr. concentration
    real(dp), intent(in)    :: mix_fracs(:)   ! backgr. deltas
    real(dp)                :: corr, tot

  ! filtering possible dummies
    if ((ind_d .lt. 1) .or. (ind_d .gt. {%NDSPEC})) return

  ! buget to correct to
    tot = sum(C(RD{%A}IND(ind_d,1:{%NCLASS})))
    corr = tot + ( mix_amount - tot ) * TSL * dilF

  ! emission of background iso-composition
    call {%TAG}_emis(ind_d, mix_amount * TSL * dilF, mix_fracs)
    tot = sum(C(RD{%A}IND(ind_d,1:{%NCLASS})))

  ! removal preserving current composition
    C(RD{%A}IND(ind_d,1:{%NCLASS})) = C(RD{%A}IND(ind_d,1:{%NCLASS})) / tot * corr

  end subroutine {%TAG}_pmix



! -----------------------------------------------------------------------------

  subroutine {%TAG}_process()

    implicit none

    integer  :: i, s
    real(dp) :: chkamnt

  ! calculating the number of specs falling below THRES
    {%TAG}_NREJCT = 0
    do i = 1, {%NDSPEC}
      chkamnt = sum(C({%RDIND}(i,1:{%NCLASS})))
      if (chkamnt .lt. THRES) then
        {%TAG}_NREJCT = {%TAG}_NREJCT + 1
      endif
    enddo           ! ndspec cycle

  ! every-step fractions/totals update
    call {%TAG}_calctotals
    call {%TAG}_calcfractions

  end subroutine {%TAG}_process



! -----------------------------------------------------------------------------

  subroutine {%TAG}_calctotals()

    implicit none

    integer  :: i

  ! here the number of total molecules is calculated from each species composition

  ! careful, {%ABBT}_T is calculated from regular!
    C(ind_{%CONF}T{%A}) = sum( C({%RDIND}(:,0)) )

  ! classes concentrations
    do i = 1, NDCLASS
      TCC(i) = sum( C(RDIND(:,i)) )
    enddo

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:O3F}
    C(ind_O3F_N_T) = TCC(1)
    C(ind_O3F_Z_T) = TCC(2)
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:O3F}

  end subroutine {%TAG}_calctotals



! -----------------------------------------------------------------------------

  subroutine {%TAG}_calcfractions()

    implicit none

    integer  :: i
    real(dp) :: tot

  ! calculating new delta-13C values
    do i = 1, {%NDSPEC}
      tot = sum(C({%RDIND}(i,1:{%NCLASS})))
      if (tot .gt. 0.0_dp) then
        CF(i,:) = C({%RDIND}(i,1:{%NCLASS})) / tot
      else
        CF(i,:) = UNDEF
      endif
    enddo        ! NISPEC-cycle

    ! totals
    tot = sum(TCC(:))
    if (tot .gt. 0.0_dp) then
      TCF(:) = TCC(:) / tot
    else
      TCF(:) = UNDEF
    endif

  end subroutine {%TAG}_calcfractions



! -----------------------------------------------------------------------------

! correction of total isotopomers budget to "regular" species budget

  subroutine {%TAG}_correct2reg()

    implicit none

    integer  :: i
    real(dp) :: tot

#ifdef CLASSES_1
  ! in case one class is defined, quitting
#ifdef DEBUG
    print *,'{%TAG}_correct2reg: no correction performed (one class)'
#endif
    return
#endif

  ! here is the ver. with corr. of ALL species to regular

    do i = 1, {%NDSPEC}
      tot = sum(C({%RDIND}(i,1:{%NCLASS})))
      if (tot .le. 0.0_dp) then
        C({%RDIND}(i,1:{%NCLASS})) = 0.0_dp
      else
        C({%RDIND}(i,1:{%NCLASS})) = ( C({%RDIND}(i,1:{%NCLASS})) * C({%RDIND}(i,0)) ) / tot
      endif
    enddo

  end subroutine {%TAG}_correct2reg



! -----------------------------------------------------------------------------

! correction of "regular" species budget to the total isotopologues budget

  subroutine {%TAG}_correct2tag()

    implicit none

    integer  :: i

#ifdef CLASSES_1
  ! in case one class is defined, quitting
#ifdef DEBUG
    print *,'{%TAG}_correct2tag: no correction performed (one class)'
#endif
    return
#endif

  ! here is the ver. with corr. of ALL species to regular
    do i = 1, NDSPEC
      C({%RDIND}(i,0)) = sum(C({%RDIND}(i,1:{%NCLASS})))
    enddo

!    C({%RDIND}(:,0)) = sum(C({%RDIND}(:,1:{%NCLASS})),dim=2)

  end subroutine {%TAG}_correct2tag



! -----------------------------------------------------------------------------

  subroutine {%TAG}_resetPTs()

  ! production tracers initialization (reset) routine

    implicit none

! {x$RESET_PTs}
! - currently disabled with use of DRPT{%ATOM}IND()

#ifdef USE_PT
    C(DRPT{%ATOM}IND(:)) = 0.0_dp    ! <-- boxmodel syntax
#endif

  end subroutine {%TAG}_resetPTs



! -----------------------------------------------------------------------------

! output file for tagged species info
  subroutine {%TAG}_init()

    implicit none

! TODO: put additional tracers/variables+units after INIT_TRAC, INIT_UNIT

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:O3F}
    call open_output_file(ncid_{%TAG}, 'caaba_mecca_{%TAG}', &
      (/   &
! {$TAG_SPECS} [%fO3_@%]
       , &
{$ELSA}       'fON_T', 'fO3_T', &
{$ELSA}       'TON', 'TO3', 'TOR' &
       , &
{$ELSA}       'NREJCT' &
       /), (/   &
! {$TAG_SPECS} [%frac $%]
       , &
{$ELSA}       'frac', 'frac', &
{$ELSA}       'mol/mol', 'mol/mol', 'mol/mol' &
       , &
{$ELSA}       'specs' &
       /), (/   &
! {$TAG_SPECS} [%@SRf_O_3(@)%]
       , &
{$ELSA}       '@SRf_N_O_N_-_O_3(to)', '@SRf_O_3(to)', &
{$ELSA}       '@SRTO_N_O_N_-_O_3', '@SRT_O_3', '@SRTO (regular)' &
       , &
{$ELSA}       '@SRnumber of rejected species' &
       /) )
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:O3F}

  end subroutine {%TAG}_init



! -----------------------------------------------------------------------------

  subroutine {%TAG}_result(model_time)

    implicit none

    real(dp), intent(in) :: model_time
    integer              :: i

  ! last value is a common parameter
    TOUT(UBOUND(TOUT)) = real({%TAG}_NREJCT)

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:O3F}
! output array: minor fractions (+total's fractions),
!               total concentrations (+regular) + NREJCT
    do i = 2, {%NCLASS}
      TOUT((i-2)*{%NSPEC}+1:(i-1)*{%NSPEC}) = CF(1:{%NSPEC},i)
    enddo

    TOUT(({%NSPEC})*({%NCLASS}-1)+1: &
         ({%NSPEC})*({%NCLASS}-1)+{%NCLASS}) = TCF(:)

    TOUT(({%NSPEC})*({%NCLASS}-1)+{%NCLASS}+1: &
         ({%NSPEC})*({%NCLASS}-1)+{%NCLASS}+{%NCLASS}) = TCC(:)

    TOUT(({%NSPEC})*({%NCLASS}-1)+{%NCLASS}+{%NCLASS}+1) = C(ind_{%CONF}T{%A})

    call write_output_file(ncid_{%TAG}, model_time, D{%A}out)
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:O3F}

  end subroutine {%TAG}_result



! -----------------------------------------------------------------------------

  subroutine {%TAG}_finish()

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:O3F}
    call close_file(ncid_{%TAG})
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:O3F}

  end subroutine {%TAG}_finish



! -----------------------------------------------------------------------------

end module {%CMODEL}_{%TAG}_box

! *****************************************************************************

