! ==============================================================================
! {%CMODEL}_tag_si
! generated: {%TIMEDATE}
!
! this module is created automatically by imtag tool - do not edit
!
! inter-configuration driver module
! level: smil si
!
! {$TAG_INFO} ! this is a template for tagging configuration
!
! [Gromov, MPI-C, 2007-2020]
! ==============================================================================

! - general tagging parameters (as conditional defines) -----------------------

#include "{%CMODEL}_tag_parameters.inc"

! -----------------------------------------------------------------------------


module {%CMODEL}_tag_si

! mecca/kpp
  use messy_mecca_kpp, only: dp

! ECHAM5/MESSy
  use messy_main_mpi_bi,        only: p_parallel_io
  use messy_main_blather_bi,    only: start_message_bi, end_message_bi
  use messy_main_data_bi,       only: nproma, ngpblks

! inter-configuration module: common utils, consts
  use {%CMODEL}_tag_common
! configurations linked
! {$CONF_LIST} [%  use {%CMODEL}_@\n  use {%CMODEL}_@_si\n%]

  implicit none
  save

  integer, parameter      :: is = 1
  private is

  public mecca_tag_init_coupling
  public mecca_tag_preprocess       ! routines to be called before mecca operates
  public mecca_tag_postprocess      !                       after
  public mecca_tag_calc_xtte4scav   ! scavenging tendencies adjustment
  public mecca_tag_sub_regtracname  ! tag->reg tracer name substitution for drydep
  public mecca_tag_sub_regtracno    ! tag->reg tracer no substitution for drydep
  public mecca_tag_init
  public mecca_tag_result

#if defined(tag_FO17) || defined(tag_IO) || defined(ATOM_O)
  integer :: tag_IO_top_MDF_H2O_level = 0
#endif

! =============================================================================

contains

! -----------------------------------------------------------------------------

  subroutine mecca_tag_init_coupling()

  ! ECHAM5/MESSy
    use messy_main_data_bi,       only: vct, nlev, nlevp1, nvclev

    implicit none

#if defined(tag_FO17) || defined(tag_IO)
  ! for the additional tag_FO17 parameterisation
    real(dp) :: alt, hypi(nlevp1), h_a(nvclev), h_b(nvclev), sfpress
    integer  :: jk
#endif

  ! LOCAL
    character(len=*), parameter :: substr = 'mecca_tag_init_coupling'
    character(len=4096) :: info

    integer :: status

  ! ----------------------------------------------------------------------------

  ! assigning regular to tagged tracers

    call start_message_bi(submodstr, 'SCANNING TRACERS', substr)

  ! scanning tracers for each confiuration
! {$CONF_LIST} [%    call @_scan_tracs()%]

  ! scan PTs if there are any
    call mecca_tag_scanPTs(info)
    if (p_parallel_io) print *, trim(info)

    call end_message_bi(submodstr, 'SCANNING TRACERS', substr)

  ! ----------------------------------------------------------------------------

#if defined(tag_FO17) || defined(tag_IO)
  ! additional tag_FO17/IO parameterisation

  ! borrowed from scav: determining the level of 200hPa (tropopause) for
  ! the D17O tropospheric water approximation

    call start_message_bi(submodstr, 'tag_FO17/IO parameterisation', substr)

    alt = 2.E4_dp       ! desired level in Pa, (= 200 hPa)

    do jk = 1, nvclev
      h_a(jk) = vct(jk)
      h_b(jk) = vct(jk+nvclev)
    enddo

    sfpress = 1.E5_dp   ! reference pressure (= 1000 hPa)

    do jk = 1, nlev+1
      hypi(jk) = h_a(jk) + h_b(jk) * sfpress
    enddo

    tag_IO_top_MDF_H2O_level = nlev + 1

    do jk = 1, nlev
      if (hypi(jk) < alt .and. hypi(jk+1) >= alt) then
        tag_IO_top_MDF_H2O_level = jk
        exit
      endif
    enddo

    if (p_parallel_io) &
      print*, '  ', substr, ': H2O MIF is reset in the levels located below ', &
        alt/100, ' hPa  (nos. ', tag_IO_top_MDF_H2O_level, ' to ', nlev, ')'

    call end_message_bi(submodstr, 'tag_FO17 parameterisation', substr)
#endif

  end subroutine mecca_tag_init_coupling


! -----------------------------------------------------------------------------

! routines to be called before mecca operates
  subroutine mecca_tag_preprocess(Conc)

    use messy_main_data_bi,       only: nlev, kproma
    use messy_main_timer,         only: lstart, current_time_step

    implicit none

    real(dp), intent(inout) :: Conc(:,:)
    integer                 :: jb, jk, jp    ! counters
#ifdef tag_CIO
#ifdef CIO_IEXBOOST
    integer, parameter      :: CIO_boost_steps = 420   ! isotope O2 kinetics boost steps
    integer, parameter      :: CIO_boost_exp = 12      !   and boost exponent
#endif
#ifdef CIO_IEXINIT
    integer, parameter      :: CIO_init_step = 42      ! init ratios to equilibrium values at this step
#endif
#endif

    jb = 0
    level_loop: do jk = 1, nlev
      kproma_loop: do jp = 1, kproma
        jb = jb + 1

      ! fast initializing configurations here at the first timestep (INITFAST_E5)
        if (lstart) then
! {$CONF_LIST} [%        call @_x1(Conc(jb,:))%]
        endif

      ! performing mass correction w.r.t. the regular mechanism
! {$CONF_LIST} [%        call @_correct(Conc(jb,:)) %]

      ! adjusting "fixed" species w.r.t. the regular mechanism
! {$CONF_LIST} [%        call @_f1(Conc(jb,:)) %]

      ! resetting PTs before integration starts
        call mecca_tag_resetPTs(Conc(jb,:))

#ifdef tag_CIO
#ifdef CIO_IEXBOOST
      ! boost clumped isotope O2 kinetics by towards expected equil. value
        if (current_time_step.lt.CIO_boost_steps) then
          Conc(jb,ind_CIOIEXBOOST) = 10._dp**CIO_boost_exp
        else
#endif
          Conc(jb,ind_CIOIEXBOOST) = 1._dp
#ifdef CIO_IEXBOOST
        endif
#endif
#ifdef CIO_IEXINIT
      ! boost clumped isotope O2 kinetics by towards expected equil. value
        if (current_time_step.eq.CIO_init_step) then
          Conc(jb,ind_CIOIEXINIT) = 1._dp
        else
#endif
          Conc(jb,ind_CIOIEXINIT) = 0_dp
#ifdef CIO_IEXINIT
        endif
#endif
#endif

#ifdef tag_FO17
      ! tag_FO17 : non-anomalously enriched main O reservoirs
      !
      ! O2
        Conc(jb,ind_FO17O2) = 0.0_dp
      !
      ! water in the troposphere
        if (jk .gt. tag_IO_top_MDF_H2O_level) then
          Conc(jb,ind_FO17H2O) = 0.0_dp
        endif
      !
      ! CO2 at the surface
        if (jk .eq. nlev) then
          Conc(jb,ind_FO17CO2) = 0.0_dp
        endif
#endif

#ifdef tag_IO
      ! tag_IO : main O reservoirs
      !
      ! O2       2 =         12.08   23.88      ; 2005.RCMS19.Barkan,Luz
        Conc(jb,ind_I16O2) = Conc(jb,ind_O2) * &
                             isofrac3a(12.08e-3_dp, VSMOW_17O, 23.88e-3_dp, VSMOW_18O, 2)
        Conc(jb,ind_I17O2) = Conc(jb,ind_O2) * &
                             isofrac3r(12.08e-3_dp, VSMOW_17O, 23.88e-3_dp, VSMOW_18O, 2)
        Conc(jb,ind_I18O2) = Conc(jb,ind_O2) * &
                             isofrac3r(23.88e-3_dp, VSMOW_18O, 12.08e-3_dp, VSMOW_17O, 2)
      !
      ! so far, zeroed-MDF water in the troposphere
        if (jk .gt. tag_IO_top_MDF_H2O_level) then
          Conc(jb,ind_I16H2O) = Conc(jb,ind_H2O) * &
                                isofrac3a(0.0_dp, VSMOW_17O, 0.0_dp, VSMOW_18O, 1)
          Conc(jb,ind_I17H2O) = Conc(jb,ind_H2O) * &
                                isofrac3r(0.0_dp, VSMOW_17O, 0.0_dp, VSMOW_18O, 1)
          Conc(jb,ind_I18H2O) = Conc(jb,ind_H2O) * &
                                isofrac3r(0.0_dp, VSMOW_18O, 0.0_dp, VSMOW_17O, 1)
        endif
      !
      ! CO2 at the surface
      !#if (jk .eq. nlev) then
      !#  Conc(jb,ind_FO17CO2) = 0.0_dp
      !#endif
#endif

      enddo kproma_loop
    enddo level_loop

#ifdef tag_CIO
#ifdef CIO_IEXBOOST
    if ( p_parallel_io .and. (current_time_step.lt.CIO_boost_steps) ) &
      print *,'mecca_tag_preprocess(si): tag_CIO - boosting kinetics by 10^',CIO_boost_exp
#else
    if ( p_parallel_io .and. lstart ) &
      print *,'mecca_tag_preprocess(si): tag_CIO - not boosting kinetics'
#endif
#ifdef CIO_IEXINIT
    if ( p_parallel_io .and. (current_time_step.eq.CIO_init_step) ) &
      print *,'mecca_tag_preprocess(si): tag_CIO - initialising ratios to equilibrium values, step: ',CIO_init_step
#endif
#endif

#ifdef DEBUG
    if (p_parallel_io) print *,'mecca_tag_preprocess(si): passed'
#endif

  end subroutine mecca_tag_preprocess



! -----------------------------------------------------------------------------

! routines to be called after mecca operates
  subroutine mecca_tag_postprocess(Conc)

    use messy_main_data_bi,       only: nlev, kproma, nproma
    use messy_main_timer,         only: time_step_len

    implicit none

    real(dp), intent(inout) :: Conc(:,:)
    integer                 :: jb, jk, jp    ! counters

    jb = 0
    level_loop: do jk = 1, nlev
      kproma_loop: do jp = 1, kproma
        jb = jb + 1

      ! skipped: {$zzzCONF_LIST} [%    call @_process()%] this is to call when iso_si will be ready
      ! skipped: calculating the number of specs falling below THRES

      ! by default, to remove deviation due to KIEs, correcting
      ! performing mass correction w.r.t. the regular mechanism
! {$CONF_LIST} [%        call @_correct(Conc(jb,:)) %]

      ! every-configuration totals update
! {$CONF_LIST} [%        call @_calctotals(Conc(jb,:)) %]

      ! converting PTs values (integral) into average reaction rates
        call mecca_tag_intPTs2arr(Conc(jb,:), time_step_len)

      enddo kproma_loop
    enddo level_loop

#ifdef DEBUG
    if (p_parallel_io) print *,'mecca_tag_postprocess(si): passed'
#endif

  end subroutine mecca_tag_postprocess



! -----------------------------------------------------------------------------

  subroutine mecca_tag_calc_xtte4scav(xtte_scav, max_lev_scav, pxtp1, kproma)

    use messy_main_data_bi,       only: nlev, nproma
    use messy_main_tracer_mem_bi, only: ntrac => ntrac_gp

    implicit none

  ! input: calculated tendencies for regulars in scav
    real(dp), intent(inout) :: xtte_scav(nproma,nlev,ntrac)   ! scav_e5: allocate(xtte_scav(nproma,nlev,ntrac)), ntrac=>ntrac_gp
  ! scav. calculation levels constraint, kproma
    integer, intent(in)     :: max_lev_scav, kproma
  ! tracer field provided by scav
    real(dp), intent(in)    :: pxtp1(nproma,nlev,ntrac)

  ! calling each configuration for tendency update
! {$CONF_LIST} [%    call @_calc_xtte4scav(xtte_scav, max_lev_scav, pxtp1, kproma) %]

  end subroutine mecca_tag_calc_xtte4scav



! -----------------------------------------------------------------------------

  subroutine mecca_tag_sub_regtracname(trindex, reg_trname)

    implicit none

  ! tracer referring index
    integer, intent(in)             :: trindex
  ! tracer name to substitute
    character(len=*), intent(inout) :: reg_trname

  ! calling each configuration to try to find the reg tracer name
! {$CONF_LIST} [%    if (@_sub_regtracname(trindex, reg_trname)) return%]

  end subroutine mecca_tag_sub_regtracname


! -----------------------------------------------------------------------------

  subroutine mecca_tag_sub_regtracno(trindex, reg_trindex)

    implicit none

  ! tracer referring index
    integer, intent(in)    :: trindex
  ! tracer no to substitute
    integer, intent(inout) :: reg_trindex

  ! calling each configuration to try to find the reg tracer name
! {$CONF_LIST} [%    if (@_sub_regtracno(trindex, reg_trindex)) return%]

  end subroutine mecca_tag_sub_regtracno


! -----------------------------------------------------------------------------

  subroutine mecca_tag_init() !!(Conc)

#ifdef DEBUG
    use messy_mecca_kpp,   only: atol, rtol, SPC_NAMES
#endif

    implicit none

  !!real(dp), intent(inout) :: Conc(:,:)

  ! some routines for configuration initialization can be called here

  ! set relative tolerances by default like in the original mech
! {$CONF_LIST} [%    call @_set_atol_rtol()%]

  ! print configuration info / parameters
    if (p_parallel_io) then
! {$CONF_LIST} [%    call @_info()%]
    endif

#ifdef DEBUG
    if (p_parallel_io) then
      print *,"mecca_tag_init -> after *_set_atol_rtol()"

! {$CONF_LIST} [%        print *,"@:"\n      do i = 1,@_NTSPEC\n        print *,SPC_NAMES(@_RSIND(i,0)), ':', &\n          atol( @_RSIND(i,:) ), rtol(@_RSIND(i,:))\n      enddo\n%]
    endif
#endif

  end subroutine mecca_tag_init



! -----------------------------------------------------------------------------

  subroutine mecca_tag_result(model_time)

    implicit none

    real(dp), intent(in) :: model_time

  ! routines for additional data output (si) for configurations to be here

#ifdef DEBUG
    if (p_parallel_io) print *,'mecca_tag_result(si): passed'
#endif

  end subroutine mecca_tag_result



! ---------------------------------------------------------------------------

end module {%CMODEL}_tag_si

! ***************************************************************************

