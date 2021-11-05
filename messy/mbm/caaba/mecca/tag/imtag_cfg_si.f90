! ==============================================================================
! {%TAG}_si
! generated: {%TIMEDATE}
!
! this module is created automatically by imtag tool - do not edit
!
! some maintenance routines for budgeting configurations (isotopes)
! level: smil si
!
! {$TAG_INFO} ! this is a template file for imtag utility
!
! [Gromov, MPIC, 2007-2020]
! ==============================================================================

! - general tagging parameters (as conditional defines) -----------------------

#include "{%CMODEL}_tag_parameters.inc"
!#define DEBUG_XTTE4SCAV

! ------------------------------------------------------------------------------

! {$CONF_PARAM}

module {%CMODEL}_{%TAG}_si

! MECCA
  use messy_mecca_kpp, only: dp
  use {%CMODEL}_tag_common
  use {%CMODEL}_{%TAG}

! BML/MESSy
  use messy_main_tracer_mem_bi, only: ntrac_gp, ti_gp ! ntrac_lg, ti_lg
  use messy_main_mpi_bi,        only: p_parallel_io

  implicit none
  save

! list of EMAC tracer indices corresponding to {%TAG}_<spec>
  integer             :: {%TAG}_DTI({%NSPEC},0:{%NCLASS})
! flags
  logical             :: {%TAG}_scan_tracs_done = .false., & ! tracers are scanned
                         {%TAG}_x1_info_output = .false.     ! info about INITFAST_E5 output

  public {%TAG}_scan_tracs
  public {%TAG}_x1               ! quick init species MR (eq. of x0 in CAABA)
  public {%TAG}_f1               ! process "fixed" specs (eq. of f0 in CAABA)
  public {%TAG}_calc_xtte4scav
  public {%TAG}_sub_regtracname
  public {%TAG}_sub_regtracno

! ==============================================================================

contains

! -----------------------------------------------------------------------------

! scanning tracers list and index original <-> tagged

  subroutine {%TAG}_scan_tracs

    ! ECHAM5/MESSy
    use messy_main_tracer_mem_bi,   only: GPTRSTR, ntrac_gp
    use messy_main_tracer_tools_bi, only: tracer_halt
    use messy_main_mpi_bi,          only: p_parallel_io
    ! MESSy
    use messy_main_tracer,        only: get_tracer_list, &
                                        t_trinfo_tp, I_TAG_REG_IDT
    use messy_main_constants_mem, only: STRLEN_MEDIUM

    implicit none
    character(len=*), parameter :: substr = '{%TAG}_scan_tracs'

    integer            :: status, js, jc
    integer, pointer   :: idt_temp(:) => NULL()
    TYPE(t_trinfo_tp), pointer, &
      dimension (:)    :: ti => NULL()
    character(len=127) :: info, ts

  ! checking if tracer indices were scanned
    if ({%TAG}_scan_tracs_done) then
      if (p_parallel_io) print *, '  ',substr,': tracers were scanned already'
      return
    endif

    if (p_parallel_io) print *, '  ',substr,': searching for corresponding regular <-> tagging tracer indices'

  ! getting the pointer to the gridpoint tracer info
    ti => ti_gp

    do js = 1, {%NSPEC}

      ! looking for original tracer
      call get_tracer_list(status, GPTRSTR, trim(SPC_NAMES({%RSIND}(js,0))), idt_temp)
      call tracer_halt(substr, status)

      ! storing tracer index, getting 1st entry (basename)
      {%TAG}_DTI(js,0) = idt_temp(1)
      ! checking if a reasonable index has arrived
      if (({%TAG}_DTI(js,0) .gt. ntrac_gp) .or. ({%TAG}_DTI(js,0) .lt. 1)) {%TAG}_DTI(js,0) = -1

      ! some info
      write (ts,'(I0)') {%TAG}_DTI(js,0)
      info = trim(SPC_NAMES({%RSIND}(js,0)))//'('//trim(ts)//') <->'

      do jc = 1, {%NCLASS}
        ! getting tagged tracers
        call get_tracer_list(status, GPTRSTR, trim(SPC_NAMES({%RSIND}(js,jc))), idt_temp)
        call tracer_halt(substr, status)

        ! storing jcth tagged tracer index, getting 1st entry (basename)
        {%TAG}_DTI(js,jc) = idt_temp(1)
        ! checking if a reasonable index has arrived
        if (({%TAG}_DTI(js,jc) .gt. ntrac_gp) .or. ({%TAG}_DTI(js,jc) .lt. 1)) then
        ! if not, abandon all
          {%TAG}_DTI(js,jc) = -1
        else
        ! if yes, update the I_TAG_REG_IDT container with the regular tracer idx
          ti({%TAG}_DTI(js,jc))%tp%meta%cask_i(I_TAG_REG_IDT) = {%TAG}_DTI(js,0)
        endif
        ! output what is found
        write (ts,'(I0)') {%TAG}_DTI(js,jc)
        info = trim(info)//' '//trim(SPC_NAMES({%RSIND}(js,jc)))//'('//trim(ts)//')'
      enddo

      if (p_parallel_io) write(*,*) '  '//trim(info)

    enddo

  ! setting the flag
    {%TAG}_scan_tracs_done = .true.

  end subroutine {%TAG}_scan_tracs



! -----------------------------------------------------------------------------

! tracers mixing ratios initialization (former x1, x0 in CAABA)

  subroutine {%TAG}_x1(C)

    implicit none

    real(dp), intent(inout) :: C(:)

    integer  :: i
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
    real(dp) :: d2H(1:{%NSPEC})
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
    real(dp) :: d13C(1:{%NSPEC})
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
    real(dp) :: d17O(1:{%NSPEC}), d18O(1:{%NSPEC})
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}
#ifndef INITFAST_E5
  ! if not set in cfg. properties, leave
    return
#endif
    if (p_parallel_io .and. .not.({%TAG}_x1_info_output) ) then
      print *,' {%TAG}_x1(): INITFAST_E5 triggered initialisation'
      {%TAG}_x1_info_output = .true.
    endif
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>case:REM}
! {$f0} [%n%]  (%    d13C({%TAG}_@) = $%)
! n - # of the class, i.e. 2 for d13 / 2 for d17, 3 for d18
! @ - species name; $ - init. value
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<case:REM}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
! {$x0} [%2%]  (%    d2H({%TAG}_@) = $%)

#ifdef ZERO_TEST
    d2H(:) = 0.0_dp
#endif

#ifdef UNIT_DELTAPERMIL
  ! 1H, 2H through delta and regular species:
    d2H(:) = d2H(:) / {%TAG}_ufac         ! de-permilizing
    C({%RSIND}(:,1)) = C({%RSIND}(:,0)) * &
      isofrac2a(d2H(:), Rstd_2H, {%NQATOM}(:))
    C({%RSIND}(:,2)) = C({%RSIND}(:,0)) * &
      isofrac2r(d2H(:), Rstd_2H, {%NQATOM}(:))
#endif
#ifdef UNIT_FRACMIN
  ! 1H, 2H through rare fraction and regular species:
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
    d13C(:) = d13C(:) / {%TAG}_ufac         ! de-permilizing
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
    d17O(:) = d17O(:) / {%TAG}_ufac         ! de-permilizing
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
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
  ! initializing using fractions given in cfg

! {$x0} [%#%]  (%    C({%RSIND}({%TAG}_@,#)) = C({%RSIND}({%TAG}_@,0)) * $%)

-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CPAR:FRAC}

  end subroutine {%TAG}_x1



! -----------------------------------------------------------------------------

! adjust tracers that are indicated as fixed (eq. of f0 in CAABA)

  subroutine {%TAG}_f1(C)

    implicit none

    real(dp), intent(inout) :: C(:)
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
    real(dp) :: d2H(1:{%NSPEC})
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
    real(dp) :: d13C(1:{%NSPEC})
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:C}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:O}
    real(dp) :: d17O(1:{%NSPEC}), d18O(1:{%NSPEC})
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:O}
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<CONF:I.+}

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>case:REM}
! {$f0} [%n%]  (%    d13C({%TAG}_@) = $%)
! n - # of the class, i.e. 2 for d13 / 2 for d17, 3 for d18
! @ - species name; $ - init. value
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<case:REM}

  ! exit if there is no fixed species
    if ({%NFIX} .lt. 1) return

->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CONF:I.+}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:H}
! {$f0} [%2%]  (%    d2H({%TAG}_@) = $%)

#ifdef UNIT_DELTAPERMIL
  ! 1H, 2H through delta and regular species:
    d2H({%FSIND}(:,0)) = d2H({%FSIND}(:,0)) / {%TAG}_ufac         ! de-permilizing

    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * isofrac2a(d2H({%FSIND}(:,0)), Rstd_2H, {%NQATOM}({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * isofrac2r(d2H({%FSIND}(:,0)), Rstd_2H, {%NQATOM}({%FSIND}(:,0)))
#endif
#ifdef UNIT_FRACMIN
  ! 12C, 13C through minor fraction and regular species:

    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * (1.0_dp - d2H({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * d2H({%FSIND}(:,0))
#endif
-<<- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {<ATOM:H}
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>ATOM:C}
! {$f0} [%2%]  (%    d13C({%TAG}_@) = $%)

#ifdef UNIT_DELTAPERMIL
  ! 12C, 13C through delta and regular species:
    d13C({%FSIND}(:,0)) = d13C({%FSIND}(:,0)) / {%TAG}_ufac         ! de-permilizing

    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * isofrac2a(d13C({%FSIND}(:,0)), Rstd_13C, {%NQATOM}({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * isofrac2r(d13C({%FSIND}(:,0)), Rstd_13C, {%NQATOM}({%FSIND}(:,0)))
#endif
#ifdef UNIT_FRACMIN
  ! 12C, 13C through minor fraction and regular species:

    C({%FSIND}(:,1)) = C({%FSIND}(:,0)) * (1.0_dp - d13C({%FSIND}(:,0)))
    C({%FSIND}(:,2)) = C({%FSIND}(:,0)) * d13C({%FSIND}(:,0))
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
    d17O({%FSIND}(:,0)) = d17O({%FSIND}(:,0)) / {%TAG}_ufac         ! de-permilizing
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
->>- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ {>CPAR:FRAC}
  ! initializing using fractions given in cfg

! {$f0} [%#%]  (%    C({%RSIND}({%TAG}_@,#)) = C({%RSIND}({%TAG}_@,0)) * $%)

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

  end subroutine {%TAG}_f1



! -----------------------------------------------------------------------------

  subroutine {%TAG}_calc_xtte4scav(xtte_scav, max_lev_scav, pxtp1, kproma)

    use messy_main_data_bi,       only: nlev, nproma
    use messy_main_tracer_mem_bi, only: ntrac => ntrac_gp

    implicit none

  ! input: calculated tendencies for regulars in scav
    real(dp), intent(inout) :: xtte_scav(nproma,nlev,ntrac)   ! scav_e5: allocate(xtte_scav(nproma,nlev,ntrac)), ntrac=>ntrac_gp
  ! scav. calculation levels constraint, kproma
    integer, intent(in)     :: max_lev_scav, kproma
  ! tracer field provided by scav
    real(dp), intent(in)    :: pxtp1(nproma,nlev,ntrac)

    integer  :: jk, jl, jt, js, jc, jtc
    real(dp) :: xtte_val, pxtp1_val
#ifdef DEBUG_XTTE4SCAV
    real(dp) :: tmp
#endif

  ! checking if tracer indices were scanned
    if (.not.{%TAG}_scan_tracs_done) call {%TAG}_scan_tracs

    loop_nlev: do jk = max_lev_scav, nlev

    ! cycling only over the set of tracers of interest
      loop_spec: do js = 1, {%NSPEC}

      ! regular tracer index
        jt = {%TAG}_DTI(js,0)
      ! false index protection
        if (jt .le. 0) cycle

#ifdef DEBUG
#ifdef tag_FO17
        if ((jk .eq. nlev) .and. (js .eq. tag_FO17_SO2)) then
           write(*,*) '- tag_FO17_calc_xtte4scav debug - FO17_SO2 ----------------------------------------------'
           write(*,*) 'proma, xtte_val, pxtp1_val, pxtp1 (class) -> ratio'
        endif
#endif
#endif

        loop_proma: do jl = 1, kproma
!!$            if (zxtp1(jl,jk,jt)*cm(jl,jk) < -1.e-15_dp) &
!!$               print*, "WARNING, ls scav negative",js,jl,jk,jt,jrow,&
!!$               zxtp1(jl,jk,jt)*cm(jl,jk), pxtp1(jl,jk,jt)
          xtte_val = xtte_scav(jl,jk,jt)
          pxtp1_val = pxtp1(jl,jk,jt)

          if (pxtp1_val .le. 0.0_dp) then  ! checking whether projected value is negative (overshot)
            xtte_val = 0._dp               !   then xtte_scav for tagged species will be zeroed
            pxtp1_val = 1._dp
#ifdef DEBUG_XTTE4SCAV
            write(*,*) "{%TAG}_calc_xtte4scav: #WARNING# pxtp1(jl,jk,jt) <= 0.0: ",jl,jk,jt,pxtp1(jl,jk,jt)
#endif
          endif

          loop_class: do jc = 1, {%NCLASS}

          ! jc class tracer
            jtc = {%TAG}_DTI(js,jc)
          ! false index protection
            if (jtc .le. 0) cycle

          ! adjusting the tendency
#ifdef DEBUG_XTTE4SCAV
            tmp = pxtp1(jl,jk,jtc) / pxtp1_val
            if (tmp .gt. 1.0_dp) then
              write(*,*) "{%TAG}_calc_xtte4scav: #WARNING# pxtp1(jl,jk,jtc)/pxtp1_val > 1.0: ",tmp,jl,jk,jc,jtc,jt
            endif
#endif
            xtte_scav(jl,jk,jtc) = xtte_val * &
               (pxtp1(jl,jk,jtc) / pxtp1_val)  ! weighting by the class fraction
#ifdef DEBUG
#ifdef tag_FO17
          ! some debug (only for FO17 configuration)
            if ((jk .eq. nlev) .and. (js .eq. tag_FO17_SO2)) &
              write(*,*) jl, ', ', xtte_val, ', ', pxtp1_val, ', ', pxtp1(jl,jk,jtc), &
                             ' ( ',jc,') -> ', pxtp1(jl,jk,jtc)/pxtp1_val
#endif
#endif
          enddo loop_class
        enddo loop_proma
      enddo loop_spec
    enddo loop_nlev

  end subroutine {%TAG}_calc_xtte4scav


! -----------------------------------------------------------------------------

  logical function {%TAG}_sub_regtracname(trindex, reg_trname)

#ifdef DEBUG
    use messy_main_mpi_bi,        only: p_parallel_io
#endif

    implicit none

  ! tag tracer index
    integer, intent(in)             :: trindex
  ! reg tracer name to substitute
    character(len=*), intent(inout) :: reg_trname

    integer  :: js, jc

  ! checking if tracer indices were scanned
    if (.not.{%TAG}_scan_tracs_done) call {%TAG}_scan_tracs

#ifdef DEBUG
    if (p_parallel_io) print *,'    {%TAG}_sub_regtracname( ',trindex,', ',reg_trname,'): '
#endif

    do js = 1, {%NSPEC}
    ! exiting immediately if the regular tracer met
      if ( trindex .eq. {%TAG}_DTI(js,0) ) then
        {%TAG}_sub_regtracname = .true.
#ifdef DEBUG
        if (p_parallel_io) print *,'REGULAR IDENTIFIED'
#endif
        return
      endif

    ! checking tagged tracers' indices
      do jc = 1, {%NCLASS}
        if ( trindex .eq. {%TAG}_DTI(js,jc) ) then
        ! substituting with the regular name
          reg_trname = trim(SPC_NAMES({%RSIND}(js,0)))   ! this relies on equal naming in tracer and mecca !
          {%TAG}_sub_regtracname = .true.
#ifdef DEBUG
          if (p_parallel_io) print *,'TAGGED IDENTIFIED, SUB: ', reg_trname
#endif
          return
        endif
      enddo
    enddo

  ! no substitution was found
    {%TAG}_sub_regtracname = .false.

#ifdef DEBUG
    if (p_parallel_io) print *,'NOT IDENTIFIED!'
#endif

  end function {%TAG}_sub_regtracname


! -----------------------------------------------------------------------------

  logical function {%TAG}_sub_regtracno(trindex, reg_trindex)

#ifdef DEBUG
    use messy_main_mpi_bi,        only: p_parallel_io
#endif

    implicit none

  ! tracer referring index
    integer, intent(in)    :: trindex
  ! tracer name to substitute
    integer, intent(inout) :: reg_trindex

    integer  :: js, jc

  ! checking if tracer indices were scanned
    if (.not.{%TAG}_scan_tracs_done) call {%TAG}_scan_tracs

#ifdef DEBUG
    if (p_parallel_io) print *,'    {%TAG}_sub_regtracno( ',trindex,', ?): '
#endif

    do js = 1, {%NSPEC}
    ! exiting immediately if the regular tracer met
      if ( trindex .eq. {%TAG}_DTI(js,0) ) then
        {%TAG}_sub_regtracno = .true.
        reg_trindex = trindex
#ifdef DEBUG
        if (p_parallel_io) print *,'REGULAR IDENTIFIED'
#endif
        return
      endif

    ! checking tagged tracers' indices
      do jc = 1, {%NCLASS}
        if ( trindex .eq. {%TAG}_DTI(js,jc) ) then
        ! substituting with the regular name
          reg_trindex = {%TAG}_DTI(js,0)
          {%TAG}_sub_regtracno = .true.
#ifdef DEBUG
          if (p_parallel_io) print *,'TAGGED IDENTIFIED, SUB: ', reg_trindex
#endif
          return
        endif
      enddo
    enddo

  ! no substitution was found
    {%TAG}_sub_regtracno = .false.

#ifdef DEBUG
    if (p_parallel_io) print *,'NOT IDENTIFIED'
#endif

  end function {%TAG}_sub_regtracno


end module {%CMODEL}_{%TAG}_si

! *****************************************************************************
