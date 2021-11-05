! ==============================================================================
! {%CMODEL}_tag_box
! generated: {%TIMEDATE}
!
! this module is created automatically by imtag tool - do not edit
!
! inter-configuration driver module
! level: smil boxmodel
!
! {$TAG_INFO} ! this is a template for cfg boxmodel module
!
! [Gromov, MPI-C, 2007-2020]
! ==============================================================================

! - general tagging parameters (as conditional defines) ------------------------

#ifdef MECCA_TAG
#include "{%CMODEL}_tag_parameters.inc"
#endif

! ------------------------------------------------------------------------------

module {%CMODEL}_tag_box

  use messy_mecca_kpp     ! dp, ...
#ifdef MECCA_TAG
  use caaba_mem,          only: timesteplen, &  ! interfaced box (caaba)
#ifdef __GFORTRAN__
                                C_ => &
#endif
                                      C
#ifdef __GFORTRAN__    /* newer gfortran versions have a run-time bug of    */
#define C C_           /* module-level import of C(:), so it has to be      */
#define c C_           /* renamed into something else, e.g. C_(:)           */
#endif
  use {%CMODEL}_tag_common  ! inter-configuration module: common utils, consts

! configurations linked
! {$CONF_LIST} [%  use {%CMODEL}_@%]
! {$CONF_LIST} [%  use {%CMODEL}_@_box%]

  implicit none
  save

  public mecca_tag_x0
  public mecca_tag_f0
  public mecca_tag_preprocess       ! routines to be called before mecca operates
  public mecca_tag_postprocess      !                       after
  public mecca_tag_emis             ! proc. for emission in all configurations
  public mecca_tag_depos            !           removal
  public mecca_tag_init
  public mecca_tag_result
  public mecca_tag_finish
#endif

! ==============================================================================

contains

#ifdef MECCA_TAG
! ------------------------------------------------------------------------------

  subroutine mecca_tag_x0()

    implicit none

! {$CONF_LIST} [%    call @_x0()%]

#ifdef DEBUG
    print *,'tag_x0: passed'
#endif

  end subroutine mecca_tag_x0


! ------------------------------------------------------------------------------

  subroutine mecca_tag_f0

    implicit none

! {$CONF_LIST} [%    call @_f0%]

#ifdef DEBUG
    print *,'tag_f0: passed'
#endif

  end subroutine mecca_tag_f0


! ------------------------------------------------------------------------------

  subroutine mecca_tag_preprocess()

    use {%CMODEL}_tag_common ! common routines
!   use {%CMODEL}_{%TAG}     ! SMCL routines
! {$CONF_LIST} [%    use {%CMODEL}_@%]

    implicit none

  ! adjusting "fixed" species before integration starts
    call mecca_tag_f0()

  ! resetting PTs before integration starts
    call mecca_tag_resetPTs(C)

#ifdef tag_FO17
  ! non-anomalously enriched main O reservoirs
  ! O2
    C(ind_FO17O2) = 0.0_dp ! C(ind_O2)
  ! water?
    C(ind_FO17H2O) = 0.0_dp ! C(ind_H2O)
#endif

#ifdef tag_FO17v2
  ! non-anomalously enriched main O reservoirs
  ! O2
  !# C(ind_FO17O2) = C(ind_O2)
  ! water?
  !# C(ind_FO17H2O) = C(ind_H2O)
#endif

#ifdef tag_FO3
  ! ozone fraction tagging - assigning tagged O3 to the regular one
  C(tag_FO3_RSIND(tag_FO3_O3,{%NCLASS})) = C(tag_FO3_RSIND(tag_FO3_O3,0))
#endif


#ifdef DEBUG
    print *,'tag_preprocess: passed'
#endif

  end subroutine mecca_tag_preprocess



! ------------------------------------------------------------------------------

  subroutine mecca_tag_postprocess(skip_extra_calc)

    implicit none

    logical, optional, intent(in) :: skip_extra_calc
    logical                       :: skip_extra_calc_act

  ! by default skip extra calculations (totals+delta)
  ! use only for optimisation purposes, when delta values are not required every timestep
#ifdef TAG_OPT_SKIP_XTRABOXPPROC
    skip_extra_calc_act = .true.
#else
    skip_extra_calc_act = .false.
#endif
    if ( present(skip_extra_calc) ) skip_extra_calc_act = skip_extra_calc

! {$CONF_LIST} [%    call @_postprocess(skip_extra_calc_act)%]

  ! by default, to remove deviation due to KIEs, correcting
! {$CONF_LIST} [%    call @_correct(C)%]

  ! converting PTs values (integral) into average reaction rates
    call mecca_tag_intPTs2arr(C, timesteplen)

#ifdef DEBUG
    print *,'tag_postprocess: passed'
#endif

  end subroutine mecca_tag_postprocess



! ------------------------------------------------------------------------------

  subroutine mecca_tag_emis(ind_r, amount, deltas)

    implicit none

    integer,  intent(in)    :: ind_r
    real(dp), intent(in)    :: amount
    real(dp), intent(in)    :: deltas(:)

! {$CONF_LIST} [%    call @_emis(ind_r, amount, deltas)%]

#ifdef DEBUG
    print *,'tag_emis(',trim(SPC_NAMES(ind_r)),', ',amount,', ',deltas,'): passed'
#endif

  end subroutine mecca_tag_emis



! ------------------------------------------------------------------------------

  subroutine mecca_tag_depos(ind_r, factor)

    implicit none

    integer,  intent(in)    :: ind_r
    real(dp), intent(in)    :: factor

! {$CONF_LIST} [%    call @_depos(ind_r, factor)%]

#ifdef DEBUG
    print *,'tag_depos(',trim(SPC_NAMES(ind_r)),', ',factor,'): passed'
#endif

  end subroutine mecca_tag_depos



! ------------------------------------------------------------------------------

  subroutine mecca_tag_init()

    implicit none

    character(len=4096) :: info

  ! scan PTs if there are any
    call mecca_tag_scanPTs(info)
    print *, trim(info)

  ! initialise tagged tracers
    call mecca_tag_x0()

! {$CONF_LIST} [%    call @_init()%]

#ifdef DEBUG
    print *,'tag_init: passed'
#endif

  end subroutine mecca_tag_init



! ------------------------------------------------------------------------------

  subroutine mecca_tag_result(model_time)

    implicit none

    real(dp), intent(in) :: model_time

! {$CONF_LIST} [%    call @_result(model_time)%]

#ifdef DEBUG
    print *,'tag_result: passed'
#endif

  end subroutine mecca_tag_result



! ------------------------------------------------------------------------------

  subroutine mecca_tag_finish()

    implicit none

! {$CONF_LIST} [%    call @_finish%]

#ifdef DEBUG
    print *,'tag_finish: passed'
#endif

  end subroutine mecca_tag_finish



! ------------------------------------------------------------------------------
#endif

end module {%CMODEL}_tag_box

! ******************************************************************************

