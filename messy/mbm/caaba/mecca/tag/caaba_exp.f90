!=============================================================================
! Experiment utils unit for boxmodel simulations with CAABA/MECCA:
!
! - mixing of the species concentrations with give background values
!
! - 0D chem-nudging (using the EHCAM/MECCA EVAL run output)
!
! - emission of a certain species from offline emission data from E5/M1
!   (prepared using icogesb scripts)
!
! - Monte-Carlo mod. for caaba simulating isoCO2 in Indian BL->CARIIBIC-2
!
!
!
! [Gromov, MPIC, 2007-2016]
!===============================================================================

#include "caaba_exp.inc"

module caaba_exp

  use messy_mecca_kpp            ! dp, ...
  use messy_main_constants_mem, only: R_gas, N_A
  use netcdf

  implicit none

  save

#ifdef OFFLEMB
! ----- offlem for boxmodel ---------------------------------------
  integer :: ncid_nudge, ncid_offlemb
#endif

#ifdef iCMCb
! ----- isoCO2 MC box runs ---------------------------------------

#ifndef NOMPI
  include 'mpif.h'
#endif

! ----- C2 data & parameters
  integer, parameter   :: C2_vars = 31, C2_samples = 28, &
                          C2_match_from = 11, &
                          C2_match_upto = 18  ! 10-19  ! range of samples no. to match
  integer, dimension(C2_samples) :: C2_no
  character(len=19), dimension(C2_samples) :: C2_date
  real(dp), dimension(C2_samples) :: &
    C2_lat, C2_lon, C2_alt, C2_strato, &
    C2_CO, C2_COe, C2_O3, C2_O3e, C2_H2O, C2_H2O_cloud, C2_Hg, &
    C2_N2O, C2_N2Oe, C2_SF6, C2_SF6e, C2_CH4, C2_CH4e, &
    C2_CO2, C2_CO2e, C2_iCO2, C2_iCO2e, C2_d13C, C2_d13Ce, C2_d18O, C2_d18Oe, &
    C2_C2H6, C2_C2H2, C2_C3H8, C2_C6H6, C2_CH3Cl

! ----- Monte-Carlo apparatus
  integer, parameter :: iCMCb_hit_sets = 2                 !> no. of independent sets of cases to hit (e.g. isotopes & trace gases)
  integer, parameter :: iCMCb_hit_cases = C2_match_upto    !> no. of cases in each set (e.g. no of obs. samples)
  integer, parameter :: iCMCb_req_hits = 1000!0            !> required no. of successful realisations per case (stop trigger + max. saved in common container)
  integer, parameter :: iCMCb_max_hits = 2*iCMCb_req_hits  !> no. of successful realisations per case (max. to save in separate container)

  integer            :: iCMCb_seed, iCMCb_harvest,      &  !> seed & current harvest no.
                        iCMCb_tries, iCMCb_real            !> #s of tries & realisations (tries ended correctly, not necessarily hit)
  integer, dimension(0:iCMCb_hit_cases,0:iCMCb_hit_sets) :: & !> zero indices are reserved for totals / intra-set hits (must be rather rare!)
                        iCMCb_hits, hits_delta, hits_tmp   !> # of hits in each set + temp. buf. for hits got/received

  integer            :: id_rnd_unif, id_rnd_norm           !> IDs of the random and uniform generators
  integer, parameter :: iCMCb_npar_norm = 14               !> # of normally-distr. paramerets ( 2xiCO, discr.d13C, precip-eq.d18O, CO+GHG )
  integer, parameter :: iCMCb_npar_unif = 11               !> # of uniformly-distr. parameters ( rate coeffs., CO+GHG @BL, +2 test )
  real(dp)           :: iCMCb_mcfc_norm(iCMCb_npar_norm)   !> normally-distr. paramerets
  real(dp)           :: iCMCb_mcfc_unif(iCMCb_npar_unif)   !> uniformly-distr. parameters
  integer, parameter :: iCMCb_mpitag_done  = 101, &        !> mpi communication tags
                        iCMCb_mpitag_hits  = 201, &
                        iCMCb_mpitag_out   = 301, &
                        iCMCb_mpitag_stats = 401


! simulation-specific variables
  logical            :: iCMCb_is_SS                        !> signals steady-state, i.e. proper realisation
  integer            :: iCMCb_SS_seqmatches                !> subsequent matches counter used to detect the SS
  integer, parameter :: iCMCb_set_FM = 0, &                !> case sets consts for convenience (0 for intra-case, i.e. "full-match")
                        iCMCb_set_IC = 1, &                !> - isotope CO2
                        iCMCb_set_TG = 2                   !> - trace gases

! parameters/variables not existing in mecca/caaba
  real(dp)           :: d18O_PC, d13C_RC                   !> isotope signatures - parameters in iCMCb

! ----- Output
! parameters & result strings for CO2 isotopes and other trace gases
  character(len=255) :: iCMCb_oprstr, iCMCb_rstrIC, iCMCb_rstrTG         !> used for text debug output
  integer            :: iCMCb_iounit(0:iCMCb_hit_cases,0:iCMCb_hit_sets) !> output files units

! output parameters & resp. units
  integer, parameter :: iCMCb_outNo = 78
  character(len=9), dimension(iCMCb_outNo), parameter :: iCMCb_outN = &
    (/ 'time_hit ', &
       'pe       ','harvest  ','tries    ','real     ', &
       'hits     ','smp_hit  ','discr_13C','d18O_eq  ', &
       'rk_diff  ','rk_fix   ','rk_resp  ','rk_mixL  ','rk_mixH  ', &
       'RLL      ','d13RLL   ','d18RLL   ','D17RLL   ', &
       'RHL      ','d13RHL   ','d18RHL   ','D17RHL   ', &
       'PC       ','d13PC    ','d18PC    ','D17PC    ', &
       'FC       ','d13FC    ','d18FC    ','D17FC    ', &
       'RC       ','d13RC    ','d18RC    ','D17RC    ', &
       'BL       ','d13BL    ','d18BL    ','D17BL    ', &
       'CM       ','d13CM    ','d18CM    ','D17CM    ', &
       'rf_indi  ','d13f_indi','d18f_indi','D17f_indi', &
       'rf_redi  ','d13f_redi','d18f_redi','D17f_redi', &
       'rf_fix   ','d13f_fix ','d18f_fix ','D17f_fix ', &
       'rf_resp  ','d13f_resp','d18f_resp','D17f_resp', &
       'rf_entL  ','rf_detL  ','rf_entH  ','rf_detH  ', &
       'CO_RHL   ','CO_RLL   ','CO_BL    ','CO_CM    ', &
       'CH4_RHL  ','CH4_RLL  ','CH4_BL   ','CH4_CM   ', &
       'N2O_RHL  ','N2O_RLL  ','N2O_BL   ','N2O_CM   ', &
       'SF6_RHL  ','SF6_RLL  ','SF6_BL   ','SF6_CM   ' /)
  character(len=9), dimension(iCMCb_outNo), parameter :: iCMCb_outU = &
    (/ 'sec 2015+', &
       'no.      ','no.      ','no.      ','no.      ', &
       'no.      ','no.C2.smp','o/oo     ','o/oo VPDB', &
       's/s      ','s/s      ','s/s      ','s/s      ','s/s      ', &
       'ppm      ','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ppm      ','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ppm      ','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ppm      ','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ppm      ','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ppm      ','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ppm      ','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ml/s/ml/s','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ml/s/ml/s','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ml/s/ml/s','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ml/s/ml/s','o/oo VPDB','o/oo VPDB','o/oo     ', &
       'ml/s/ml/s','ml/s/ml/s','ml/s/ml/s','ml/s/ml/s', &
       'ppb      ','ppb      ','ppb      ','ppb      ', &
       'ppb      ','ppb      ','ppb      ','ppb      ', &
       'ppb      ','ppb      ','ppb      ','ppb      ', &
       'ppt      ','ppt      ','ppt      ','ppt      ' /)
  real(dp)           :: iCMCb_out(iCMCb_outNo)
#endif

contains



#ifdef iCMCb

! ----- reads C2 data for iCMCb runs ---------------------------------------

  subroutine iCMCb_read_C2_data

    use messy_main_tools, only: find_next_free_unit
    use caaba_mem,        only: cair

    use, intrinsic  :: iso_fortran_env

    implicit none

    character(len=*), parameter :: ifname = 'C2_isoCO2_2015-08.dat' !> C2 data filename
    integer                     :: i, iou, rc, rn

    iou = find_next_free_unit(100,200)
    open(unit=iou, file=trim(ifname), status='old')
  ! skip header lines (2)
    do i = 1,2
      read(iou, *)
    enddo

    rn = 0
    readloop: do i=1, C2_samples
      read(iou, *, iostat=rc) &
        C2_no(i), C2_date(i), C2_lat(i), C2_lon(i), C2_alt(i), C2_strato(i), &
        C2_CO(i), C2_COe(i), C2_O3(i), C2_H2O(i), C2_H2O_cloud(i), C2_Hg(i), &
        C2_N2O(i), C2_N2Oe(i), C2_SF6(i), C2_SF6e(i), C2_CH4(i), C2_CH4e(i), &
        C2_CO2(i), C2_CO2e(i), C2_iCO2(i), C2_iCO2e(i), C2_d13C(i), C2_d13Ce(i), C2_d18O(i), C2_d18Oe(i), &
        C2_C2H6(i), C2_C2H2(i), C2_C3H8(i), C2_C6H6(i), C2_CH3Cl(i)

#ifdef DEBUG
      write(*,'(I0,10(F10.3))'), C2_no(i), C2_CO2(i), C2_CO2e(i), C2_d13C(i), C2_d13Ce(i), C2_d18O(i), C2_d18Oe(i), C2_CH4(i), C2_CH4e(i), C2_CO(i), C2_COe(i)
#endif
      if ( rc /= 0 ) then
        if ( rc == iostat_end ) then
          exit readloop
        else
          print *, 'iCMCb_read_C2_data: error reading "'//trim(ifname)//'", code: ', rc
          stop
        endif
      else
        rn = rn + 1
      endif
    enddo readloop
    print *, 'iCMCb_read_C2_data: read ',rn,' sample record(s)'

  ! converting MRs to conc.
    C2_CO2(:)  = C2_CO2(:)  * 1e-06 * cair
    C2_CO2e(:) = C2_CO2e(:) * 1e-06 * cair
    C2_CO(:)   = C2_CO(:)   * 1e-09 * cair
    C2_COe(:)  = C2_COe(:)  * 1e-09 * cair
    C2_O3(:)   = C2_O3(:)   * 1e-09 * cair
    C2_O3e(:)  = max(0.02*C2_O3(:),1.0) * 1e-09 * cair  ! 1 ppbv or 2% (the higher), http://wiki.caribic-atmospheric.com/index.php?title=Ozone_(O3)_analyser
    C2_CH4(:)  = C2_CH4(:)  * 1e-09 * cair
    C2_CH4e(:) = C2_CH4e(:) * 1e-09 * cair
    C2_N2O(:)  = C2_N2O(:)  * 1e-09 * cair
    C2_N2Oe(:) = C2_N2Oe(:) * 1e-09 * cair
    C2_SF6(:)  = C2_SF6(:)  * 1e-12 * cair
    C2_SF6e(:) = C2_SF6e(:) * 1e-12 * cair
!   C2_xxx(:)  = C2_xxx(:)  * 1e-06 * cair
!   C2_xxxe(:) = C2_xxxe(:) * 1e-06 * cair

    close(iou)
#ifdef DEBUG
    print *,'iCMCb_read_C2_data: finish'
#endif

  end subroutine iCMCb_read_C2_data


! ----- returns current time ---------------------------------------------------
!> current time s of julian days
  real(dp) function iCMCb_time()
    use messy_main_timer, only: julian_day
    implicit none
    intrinsic real, date_and_time
    integer    :: atime(8)
    real(dp)   :: result
  ! getting current time (borrowed from qtimer)
    call date_and_time(values=atime)
    result = julian_day( &
            (real(atime(3), dp) +  &              ! day
            real(atime(5), dp)/24.0_dp + &        ! hours -> days
            real(atime(6), dp)/1440.0_dp + &      ! minutes -> days
            real(atime(7), dp)/86400.0_dp + &     ! seconds -> days
            real(atime(8), dp)/86400000.0_dp ) &  ! milliseconds -> days
            , atime(2), atime(1) )
    iCMCb_time = (result - 2457023.5_dp)*(60_dp*60_dp*24_dp)      ! -> seconds from 01-01-2015 00:00:00
  end function iCMCb_time


! ----- MCb initialisation ------------------------------------------------------
  subroutine iCMCb_init(seed_in)

    use messy_main_tools, only: find_next_free_unit
    use mo_mpi
    use caaba_io,         only: time_string

    use messy_main_rnd,   only: RND_MTW, RND_MTW_GAUSS, &
                                RND_LUX, RND_LUX_GAUSS, &
                                rnd_init

    implicit none

    character(len=*), parameter   :: modstr = 'caaba-amc'

    integer, intent(in), optional :: seed_in
    integer                       :: seed_bcast
    integer                       :: i, rc, rn, re1, re2
    character(len=1024)           :: str, sic, stg
    real                          :: rseed


  ! ----- initialising parallel environment ------------------------------------
    call p_start(modstr)
    if (p_parallel_io) then
      print *,'iCMCb_init@master: hello from pe = ',p_pe, '/', p_nprocs
    else
      print *,'iCMCb_init@slave: greeting from pe = ',p_pe, '/', p_nprocs
    ! call sleep(1) ! this holds slaves back, so that the master initialises the output
    endif


  ! ----- initialising random apparatus ----------------------------------------
    if (p_parallel_io) then
       if ( present(seed_in) ) then
       ! using predefined seed
         iCMCb_seed = seed_in
         print *,'iCMCb_init@master: predefined seed (e.g. mcexp_seed) = ', iCMCb_seed
       else
       ! getting randomised seed
       ! simple clock
        !call system_clock(iCMCb_seed)
       ! xor:PID+clock
         call init_random_seed()
         call RANDOM_NUMBER(rseed)
         iCMCb_seed = int(HUGE(1)*rseed)
         print *,'iCMCb_init@master: no predefined seed, using random (xor:PID+clock) = ', iCMCb_seed
       endif
    endif

    call p_bcast(iCMCb_seed, p_io)  ! sending/receiving @slaves

#ifdef DEBUG
    if ( .not.p_parallel_io) print *,'iCMCb_init@slave(',p_pe,'): received seed = ', iCMCb_seed
#endif

  ! init RNGs
  ! normally dist. randoms
    call rnd_init(re1, id_rnd_norm, RND_MTW_GAUSS, iCMCb_seed)
  ! uniformly dist. randoms
    call rnd_init(re2, id_rnd_unif, RND_MTW, iCMCb_seed)
    if ( (re1+re2)/=0 ) then
      write (*,*) 'iCMCb_init: rnd_init error(s): ', re1, re2
      stop
    endif

  ! used to flag the start
    iCMCb_SS_seqmatches = -1

  ! number of harvests, tries, realisations and hits
    iCMCb_harvest = 0
    iCMCb_tries = 0
    iCMCb_real = 0
    iCMCb_hits(:,:) = 0


  ! ----- custom init, e.g. read obs. samples ----------------------------------
    call iCMCb_read_C2_data


  ! ----- done for slave
    if ( .not.p_parallel_io ) return


  ! ----- preparing output -----------------------------------------------------
  ! tweaking time origin for netcdf files to monitor real time w/suff. precision
    time_string = 'seconds since 2015-01-01 00:00:00'

  ! files
  ! correct realisations (steady-state)
    call iCMCb_init_output( 1, iCMCb_set_FM, 'iCMCb-ss')

  ! intra-case (full-match) realisations
    call iCMCb_init_output( 0, iCMCb_set_FM, 'iCMCb-fm')

  ! set hit containers (all cases)
  !//  do i = 1, iCMCb_hit_sets, ...
    call iCMCb_init_output( 0, iCMCb_set_IC, 'iCMCb-ic')
    call iCMCb_init_output( 0, iCMCb_set_TG, 'iCMCb-tg')

  ! case hit containers
    do i = C2_match_from, C2_match_upto
      write(str,'(I2.2)'), i
      call iCMCb_init_output(i, iCMCb_set_IC, 'iCMCb-ic'//trim(str))
      call iCMCb_init_output(i, iCMCb_set_TG, 'iCMCb-tg'//trim(str))
    enddo
#ifdef DEBUG
    print *,'iCMCb_init: output initialised'
#endif

    if (p_parallel_io) print *,'iCMCb_init: done'

  contains

    subroutine iCMCb_init_output(cn,sn,fname)

      use caaba_io

      implicit none

      integer, intent(in) :: cn, &       !> case no (sample no.)
                             sn          !> set no (IC, TG, etc. => FM )
      character(*), intent(in) :: fname  !> output filename
      integer :: i
      character(len=2048) :: s

    ! getting free unit
      iCMCb_iounit(cn,sn) = find_next_free_unit(100,200)

#ifdef iCMCb_netcdf
      call open_output_file(iCMCb_iounit(cn,sn), fname, iCMCb_outN, iCMCb_outU, master=p_parallel_io)
#else
    ! master initialises the output
      if (p_parallel_io) then
        open(unit=iCMCb_iounit(cn,sn), file=fname//'.dat', status='replace') !, share='denynone')
      ! writing variable names
        s = 'descriptor \'//NEW_LINE('A')
        do i=lbound(iCMCb_outN,1), ubound(iCMCb_outN,1)
          s = s//trim(iCMCb_outN(i))//'	'
        enddo
        write(iCMCb_iounit(cn,sn), *), trim(str)
        close(iCMCb_iounit(cn,sn))
    ! slave 
      else
      ! not every compiler allows opening shared :[
      ! open(unit=iCMCb_iounit(cn,sn), file=fname//'.dat', status='old') !, share='denynone')
      endif
#endif

    end subroutine iCMCb_init_output

  end subroutine iCMCb_init


! ----- deinit, close output ---------------------------------------------------
  subroutine iCMCb_finish
    use caaba_io, only: close_file
    use mo_mpi, only: p_pe, p_nprocs, p_stop, p_parallel_io
    use messy_main_rnd,   only: rnd_finish
    implicit none
    integer :: i

  ! ----- print final stats
    call iCMCb_print_stats
    print *,''

  ! ----- de-init PE
    if (p_parallel_io) then
      print *,'iCMCb_finish: good-bye from the master, pe = ',p_pe, '/', p_nprocs
    ! call sleep(1) ! this holds slaves back, so that the master finishes the last
    else
      print *,'iCMCb_finish: ciao from a slave, pe = ',p_pe, '/', p_nprocs
    endif

#ifdef DEBUG
    print *,'iCMCb_finish: done, calling p_stop'
#endif
    call p_stop

  ! ----- de-init RNGs
    call rnd_finish(id_rnd_norm)
    call rnd_finish(id_rnd_unif)

  ! ----- exit for a slave
    if (.not.p_parallel_io) return

  ! ----- closing output -------------------------------------------------------
    do i = 0, 0
#ifdef iCMCb_netcdf
      call close_file(iCMCb_iounit(i,iCMCb_set_IC), master=p_parallel_io)
      call close_file(iCMCb_iounit(i,iCMCb_set_TG), master=p_parallel_io)
#else
      close(iCMCb_iounit(i,iCMCb_set_IC))
      close(iCMCb_iounit(i,iCMCb_set_TG))
#endif
    enddo
    call close_file(iCMCb_iounit(1,iCMCb_set_FM), master=p_parallel_io)
    call close_file(iCMCb_iounit(0,iCMCb_set_FM), master=p_parallel_io)

    do i = C2_match_from, C2_match_upto
#ifdef iCMCb_netcdf
      call close_file(iCMCb_iounit(i,iCMCb_set_IC), master=p_parallel_io)
      call close_file(iCMCb_iounit(i,iCMCb_set_TG), master=p_parallel_io)
#else
      close(iCMCb_iounit(i,iCMCb_set_IC), master=p_parallel_io)
      close(iCMCb_iounit(i,iCMCb_set_TG), master=p_parallel_io)
#endif
    enddo

  end subroutine iCMCb_finish


! ----- output stats -----------------------------------------------------------
  subroutine iCMCb_print_stats
    use mo_mpi
    implicit none
#ifndef NOMPI
    integer      :: p_status(MPI_STATUS_SIZE) !> standard information of MPI_RECV
    integer      :: p_source                  !> slave no.
#endif
    logical      :: incoming
    integer :: i, r_tries
    character(len=1024) :: sic, stg, sfm, s

    if ( p_parallel_io ) then

      sic=''
      stg=''
      sfm=''

      do i = C2_match_from, C2_match_upto
        write(s,'(I0)'), iCMCB_hits(i,iCMCb_set_IC)
        sic = trim(sic)//'-'//trim(s)
        write(s,'(I0)'), iCMCB_hits(i,iCMCb_set_TG)
        stg = trim(stg)//'-'//trim(s)
        write(s,'(I0)'), iCMCB_hits(i,iCMCb_set_FM)
        sfm = trim(sfm)//'-'//trim(s)
      enddo
      sic = trim(sic)//'-'
      stg = trim(stg)//'-'
      sfm = trim(sfm)//'-'

#ifdef iCMCb_show_match_progress
      print *,''
#endif
      write(*,'(3(A,I0),3(2(A,I0),A))', advance="no"), &
        ' #', p_pe,' > try:',iCMCb_tries/1000,'k real:',iCMCb_real, &
        ' ic:',iCMCb_hits(0,iCMCb_set_IC),'/',sum(iCMCb_hits(1:iCMCb_hit_cases,iCMCb_set_IC)),'('//trim(sic)//')', &
        ' tg:',iCMCb_hits(0,iCMCb_set_TG),'/',sum(iCMCb_hits(1:iCMCb_hit_cases,iCMCb_set_TG)),'('//trim(stg)//')', &
        ' fm:',iCMCb_hits(0,iCMCb_set_FM),'/',sum(iCMCb_hits(1:iCMCb_hit_cases,iCMCb_set_FM)),'('//trim(sfm)//')'

#ifndef NOMPI
    ! listening to incoming stats from slaves
      do
        incoming = .false.
        call MPI_iprobe(MPI_ANY_SOURCE, iCMCb_mpitag_stats, p_all_comm, &
                        incoming, p_status, p_error)
#ifdef DEBUG
        if (p_error/=MPI_SUCCESS) then
          write (nerr,'(a,i4,a,i4,a)') ' MPI_iprobe on pe #', p_pe, &
                      ' for tag ', iCMCb_mpitag_stats, ' failed.'
          write (nerr,'(a,i4)') ' Error = ', p_error
          call p_abort
        endif
#endif
      ! there is an update from a slave
        if (incoming) then
          p_source = p_status(MPI_SOURCE)
        ! receiving stats from a slave
          call p_recv(r_tries, p_source, iCMCb_mpitag_stats)
        ! and print right away
          write(*,'(A,I0,A,I0,A)',advance="no"), ' #', p_source,'(',r_tries/1000,'k)'
        else
          print *
          exit  ! repeating until all slaves are updated (no incoming requests)
        endif
      enddo
    else
    ! send stats to master
      call p_isend(iCMCb_tries, p_io, iCMCb_mpitag_stats)  ! sending tries
#endif
    endif

  end subroutine iCMCb_print_stats


! ----- generate new initial state of parameters -------------------------------
  subroutine iCMCb_new_x0

    use caaba_mem,          only: C, temp, cair, press
    use messy_main_rnd,     only: rnd_number
    use mo_mpi,             only: p_pe, p_nprocs
    use messy_mecca_tag_box

    implicit none

  ! harvest new vector of randoms
  !   cycle-harvesting until this PE's turn is reached (mod of harvest by #PEs = PE+1)
    do
      call rnd_number(id_rnd_norm, iCMCb_mcfc_norm(:))
      call rnd_number(id_rnd_unif, iCMCb_mcfc_unif(:))
      iCMCb_harvest = iCMCb_harvest + 1
      if ( mod(iCMCb_harvest,p_nprocs) .eq. p_pe ) exit
    enddo
#ifdef DEBUG
    print *,'iCMCb_new_x0: pe =',p_pe,'harvest =',iCMCb_harvest
#endif

  ! translate vector into parameters
  !
  ! ----- isotope CO2 -----
  !
  ! normally distributed:
  !
  ! #1-3/4-6: Remote CO2 [R(H/L)L]
  !
  ! Uncertainties are from NOAA ESRL CCCGASN
  ! ftp://aftp.cmdl.noaa.gov/data/trace_gases/co2/flask/surface/README_surface_flask_co2.html
  ! ftp://aftp.cmdl.noaa.gov/data/trace_gases/co2c13/flask/surface/README_surface_flask_co2c13.html
  ! averages & uncertainties are derived from NOAA ERSL GMD event data from 01-14/08/2008 (see GMD.xlsx)
    !
    C(ind_RLL) = ( 383.832 + 0.245 * ( iCMCb_mcfc_norm(01) ) ) * 1e-06_dp * cair
    d13C(tag_IC_RLL) = ( -8.205 - 0.05 + 0.051 * ( iCMCb_mcfc_norm(02) ) )
    d18O(tag_IO_RLL) = (  0.853 + 0.084 * ( iCMCb_mcfc_norm(03) ) )
    d17O(tag_IO_RLL) = delta17Opm( d18O(tag_IO_RLL), -0.1_dp )   ! D17O is SET to -0.1 per mil
    !
    call tag_IC_set( ind_RLL, C(ind_RLL), (/ d13C(tag_IC_RLL) /) )
    call tag_IO_set( ind_RLL, C(ind_RLL), (/ d17O(tag_IO_RLL), d18O(tag_IO_RLL) /) )
    !
    C(ind_RHL) = ( 384.270 + 0.283 * ( iCMCb_mcfc_norm(04) ) ) * 1e-06_dp * cair
    d13C(tag_IC_RHL) = ( -8.202 - 0.05 + 0.058 * ( iCMCb_mcfc_norm(05) ) )
    d18O(tag_IO_RHL) = (  0.318 + 0.109 * ( iCMCb_mcfc_norm(06) ) )
    d17O(tag_IO_RHL) = delta17Opm( d18O(tag_IO_RHL), -0.1_dp )   ! D17O is SET to -0.1 per mil
    !
    call tag_IC_set( ind_RHL, C(ind_RHL), (/ d13C(tag_IC_RHL) /) )
    call tag_IO_set( ind_RHL, C(ind_RHL), (/ d17O(tag_IO_RHL), d18O(tag_IO_RHL) /) )

  ! #9-12: Biosphere-exchange CO2 parameters
  !
  ! Photosynthetic [PC] / respired [RC] CO2
  ! + d18O of both [PC] and [RC] is equal to that of CO2 equilibrated with veg./soil water (-6 o/oo, Assonov, PC)
!   d18O_PC   = (  -6.0 +  1.0 * ( iCMCb_mcfc_norm(xx) ) )    ! d18O: SD of -+1.0 per mil is chosen arbitrarily so far
    d18O_PC   = ( -20.0 + 20.0 * ( iCMCb_mcfc_unif(06) ) )
  !
  ! Average plant 13C discrimination value (-19 o/oo, Assonov, PC)
  ! - currently simulated explicitely
!   discr_13C = (  23.0 +  3.00 * ( iCMCb_mcfc_norm(xx) ) )
!   discr_13C = (  10.0 + 20.00 * ( iCMCb_mcfc_unif(xx) ) )
  !
  ! d13C of respired CO2 [RC]
!   d13C_RC   = d13C(tag_IC_RLL) - discr_13C                  ! d13C: inherits uncertainty via RLL d13C and discrimination
    d13C_RC   = ( -35.0 + 25.0 * ( iCMCb_mcfc_unif(07) ) )
  !
  ! MR of chloroplast CO2 [PC]
  ! - using the RLL MR and 13C discr. to estimate chloroplast MR [PC]
  !   (@@ Farquhar's model,, see ref. in [2015.TG5.7.Affek&Yakir])
  !   a = 4.4 o/oo, b = 29 o/oo, Ca = C(ind_RLL)
  !   [PC] is Cc
!   C(ind_PC)       = C(ind_RLL) * ( discr_13C - discr_13C_a )/( discr_13C_b - discr_13C_a )
    C(ind_PC)       = 0.0
  !
    d18O(tag_IO_PC) = d18O_PC
    d17O(tag_IO_PC) = delta17Opm( d18O_PC, 0.0_dp )            ! D17O is SET to 0 (precipitation => MWL)
    call tag_IO_set( ind_PC, C(ind_PC), (/ d17O(tag_IO_PC), d18O(tag_IO_PC) /) )
  !
  ! MR of respired CO2 [RC]
    C(ind_RC)       = C(ind_RLL)     ! arbitrarily set to that of [PC] or [RLL] init. value
    d13C(tag_IC_RC) = d13C_RC
    d18O(tag_IO_RC) = d18O_PC        ! same signature as for [PC]
    d17O(tag_IO_RC) = delta17Opm( d18O(tag_IO_RC), 0.0_dp )    ! D17O is SET to 0 (precipitation => MWL)
    call tag_IC_set( ind_RC, C(ind_RC), (/ d13C(tag_IC_RC) /) )
    call tag_IO_set( ind_RC, C(ind_RC), (/ d17O(tag_IO_RC), d18O(tag_IO_RC) /) )

  ! uniformly distributed:
  ! -> actual unknowns, MC factors spread 0-1
  !
   !iCMCb_mcfc_unif(:) = ( iCMCb_mcfc_unif(:) - 0.5_dp )  ! this centers factors around 0
  !
  ! advection const (reference)
    k_adv = 1./(24.*60.*60.)  ! -- gives ppm/day throughput, reference rate
  !
  ! Farquhar model-related set
   !k_diff = k_adv * 3.0 * iCMCb_mcfc_unif(01)              ! in/retro diffusion to chloroplast rates const. (k_resp)
   !k_fix = k_diff * ( C(ind_RLL)/C(ind_PC) - 1. )          ! fixation rate (k_fix = k_diff * (Ca/Cc-1) (see @@ above))
   !k_resp = k_adv * ( 0.50 + 1.00 * iCMCb_mcfc_unif(02) )  ! respiration rate const. (k_resp)
  !
  ! A "rather random" set
    k_fix  = k_adv * ( 0.10 +  3.0 * iCMCb_mcfc_unif(01) )   ! fixation rate (k_fix)
    k_diff = k_fix * ( 0.01 +  5.0 * iCMCb_mcfc_unif(02) )   ! in/retro diffusion to chloroplast rates const. (k_resp) - made w.r.t. k_fix to yield uniform discr_13C
    k_resp = k_adv * ( 0.10 +  3.0 * iCMCb_mcfc_unif(03) )   ! respiration rate const. (k_resp)

  ! mixing processes
   !k_mixL = k_adv * ( 0.10 + 14.9 * iCMCb_mcfc_unif(04) )   ! HL air mixing rate const. (k_mixL) , 1:10 to 15:1  ! &&&
   !k_mixH = k_adv * ( 0.10 + 14.9 * iCMCb_mcfc_unif(05) )   ! LL air Mixing rate const. (k_mixH) , 1:10 to 15:1
    k_mixL = k_adv * 15._dp**( 2.0_dp * iCMCb_mcfc_unif(04) - 1.0_dp )   ! HL air mixing rate const. (k_mixL) , 1:15 to 15:1 ratio-uniform
    k_mixH = k_adv * 15._dp**( 2.0_dp * iCMCb_mcfc_unif(05) - 1.0_dp )   ! LL air Mixing rate const. (k_mixH) , 1:15 to 15:1 ratio-uniform

  ! ----- trace gases -----
  !
  ! uniformly distributed:
  ! Unknown non-CO2 concentrations in BL
    C(ind_BL_CO)  = (   80.0 + 500.0 * ( iCMCb_mcfc_unif(08) ) ) * 1e-09_dp * cair
    C(ind_BL_CH4) = ( 1800.0 + 400.0 * ( iCMCb_mcfc_unif(09) ) ) * 1e-09_dp * cair
    C(ind_BL_N2O) = (  316.0 +  10.0 * ( iCMCb_mcfc_unif(10) ) ) * 1e-09_dp * cair
    C(ind_BL_SF6) = (    6.4 +   0.6 * ( iCMCb_mcfc_unif(11) ) ) * 1e-12_dp * cair
   !C(ind_BL_O3)  = (   30.0 +   5.0 * ( iCMCb_mcfc_unif(xx) ) ) * 1e-09_dp * cair
  !
  ! normally distributed:
  ! Observed remote CO+GHGs / Averages & uncertainties are from NOAA ESRL GMD (MLO)
  !
  ! RHL (MLO)
    C(ind_RHL_CO)  = (   65.60  + 2.4   * ( iCMCb_mcfc_norm(07) ) ) * 1e-09_dp * cair
    C(ind_RHL_CH4) = ( 1779.38  + 1.40  * ( iCMCb_mcfc_norm(08) ) ) * 1e-09_dp * cair
    C(ind_RHL_N2O) = (  322.29  + 0.42  * ( iCMCb_mcfc_norm(09) ) ) * 1e-09_dp * cair
    C(ind_RHL_SF6) = (    6.518 + 0.027 * ( iCMCb_mcfc_norm(10) ) ) * 1e-12_dp * cair
   !C(ind_RHL_O3)  = (   OZ.0   + 0.OZ  * ( iCMCb_mcfc_norm(xx) ) ) * 1e-09_dp * cair
  !
  ! RLL (PON/PBL)
    C(ind_RLL_CO)  = (   95.0   + 20.0  * ( iCMCb_mcfc_norm(11) ) ) * 1e-09_dp * cair
    C(ind_RLL_CH4) = ( 1785.0   + 15.0  * ( iCMCb_mcfc_norm(12) ) ) * 1e-09_dp * cair
    C(ind_RLL_N2O) = (  326.5   +  2.0  * ( iCMCb_mcfc_norm(13) ) ) * 1e-09_dp * cair
    C(ind_RLL_SF6) = (    6.3   +  0.1  * ( iCMCb_mcfc_norm(14) ) ) * 1e-12_dp * cair
   !C(ind_RLL_O3)  = (   OZ.0   + 0.OZ  * ( iCMCb_mcfc_norm(xx) ) ) * 1e-09_dp * cair

#ifdef DEBUG
    iCMCb_oprstr = ""
    write(iCMCb_oprstr,"(2(I0,'	'),16(F0.3,'	'))"), iCMCb_tries, iCMCb_seed, &
               C(ind_RLL)/cair*1e6, d13C(tag_IC_RLL), d18O(tag_IO_RLL), &
               C(ind_RHL)/cair*1e6, d13C(tag_IC_RHL), d18O(tag_IO_RHL), &
               discr_13C, d18O_PC, d13C_RC, &
               k_diff/k_adv, k_fix/k_adv, k_resp/k_adv, k_mixL/k_adv, k_mixH/k_adv
#ifdef DEEPDEBUG
    print *, 'iCMCb_new_x0:	',trim(iCMCb_oprstr)
#else
    print *, trim(iCMCb_oprstr)
#endif
#endif

  end subroutine iCMCb_new_x0


! ----- checks whether multiple pairs of values match within tolerance ---------
  logical function iCMCb_check_pairs(cp_no, cp)

    implicit none
    intrinsic dabs, real

    integer, intent(in)                      :: cp_no        ! # of check-pairs
    real(dp), dimension(3,cp_no), intent(in) :: cp           ! check-pairs(+tolerances) themselves
    integer, dimension(cp_no)                :: unmatch      ! no. of non-matches

    where ( dabs(cp(1,:)-cp(2,:)) .gt. cp(3,:) )
      unmatch(:) = 1
    elsewhere
      unmatch(:) = 0
    endwhere

    if ( sum(unmatch(:)) .eq. 0 ) then
      iCMCb_check_pairs = .true.
    else
      iCMCb_check_pairs = .false.
    endif

  end function iCMCb_check_pairs


! ----- prepare (fill up) the output array -------------------------------------
  subroutine iCMCb_prep_output

    use mo_mpi, only: p_pe
    use caaba_mem, only: C, temp, cair, press
    use messy_mecca_tag_box

    implicit none

    real(dp)             :: ppm, ppb, ppt, f_advin

  ! concentrations -> MRs
    ppm = 1e6 / cair
    ppb = ppm * 1e3
    ppt = ppb * 1e3

  ! calculating discr_13C
  ! using Farquhar model
  !   discr_13C = discr_13C_a + ( discr_13C_b - discr_13C_a ) * C(ind_PC)/C(ind_BL)
  ! actual value = diff.b/w fixed CO2 and BL
      discr_13C = d13C(tag_IC_BL)-d13C(tag_IC_FC)

  ! advective influx (used to calculate reduced fluxes)
    f_advin = C(ind_RLL)*k_adv

  ! output array
    iCMCb_out(:) = (/ iCMCb_time(), &
      real(p_pe,dp), real(iCMCb_harvest,dp), real(iCMCb_tries,dp), real(iCMCb_real,dp), &
      -1.0_dp, -1.0_dp, discr_13C, d18O_PC, &
      k_diff/k_adv, k_fix/k_diff, k_resp/k_adv, k_mixL/k_adv, k_mixH/k_adv, &
      ppm*C(ind_RLL), d13C(tag_IC_RLL), d18O(tag_IO_RLL), DC17O(tag_IO_RLL), &
      ppm*C(ind_RHL), d13C(tag_IC_RHL), d18O(tag_IO_RHL), DC17O(tag_IO_RHL), &
      ppm*C(ind_PC),  d13C(tag_IC_PC),  d18O(tag_IO_PC),  DC17O(tag_IO_PC),  &
      ppm*C(ind_FC),  d13C(tag_IC_FC),  d18O(tag_IO_FC),  DC17O(tag_IO_FC),  &
      ppm*C(ind_RC),  d13C(tag_IC_RC),  d18O(tag_IO_RC),  DC17O(tag_IO_RC),  &
      ppm*C(ind_BL),  d13C(tag_IC_BL),  d18O(tag_IO_BL),  DC17O(tag_IO_BL),  &
      ppm*C(ind_CM),  d13C(tag_IC_CM),  d18O(tag_IO_CM),  DC17O(tag_IO_CM),  &
      C(ind_FXindi)/f_advin, d13C(tag_IC_FXindi), d18O(tag_IO_FXindi), DC17O(tag_IO_FXindi), &
      C(ind_FXredi)/f_advin, d13C(tag_IC_FXredi), d18O(tag_IO_FXredi), DC17O(tag_IO_FXredi), &
      C(ind_FXfix)/f_advin,  d13C(tag_IC_FXfix),  d18O(tag_IO_FXfix),  DC17O(tag_IO_FXfix),  &
      C(ind_FXresp)/f_advin, d13C(tag_IC_FXresp), d18O(tag_IO_FXresp), DC17O(tag_IO_FXresp), &
      C(ind_FXentL)/f_advin, C(ind_FXdetL)/f_advin, C(ind_FXentH)/f_advin, C(ind_FXdetH)/f_advin, &
      ppb*C(ind_RHL_CO),  ppb*C(ind_RLL_CO),  ppb*C(ind_BL_CO),  ppb*C(ind_CM_CO),  &
      ppb*C(ind_RHL_CH4), ppb*C(ind_RLL_CH4), ppb*C(ind_BL_CH4), ppb*C(ind_CM_CH4), &
      ppb*C(ind_RHL_N2O), ppb*C(ind_RLL_N2O), ppb*C(ind_BL_N2O), ppb*C(ind_CM_N2O), &
      ppt*C(ind_RHL_SF6), ppt*C(ind_RLL_SF6), ppt*C(ind_BL_SF6), ppt*C(ind_CM_SF6) &
      /)

  end subroutine iCMCb_prep_output


! ----- find matching cases ----------------------------------------------------
  subroutine iCMCb_find_matching_cases

    use caaba_mem, only: C!, temp, cair, press
    use messy_mecca_tag_box

    implicit none

    integer              :: sn, cn            ! set & case nos.

    if ( iCMCb_is_SS ) then

    ! hit counters
      hits_delta(:,:) = 0

      do cn = C2_match_from, C2_match_upto    ! TODO: -> 1, iCMCb_hit_cases

      ! CO2 + isotopes
#ifdef DEEPDEBUG
        print *,cn,cair,(/ (/ C(ind_CM), C2_CO2(cn), C2_CO2e(cn) /), &
                      (/ d13C(tag_IC_CM), C2_d13C(cn), C2_d13Ce(cn) /), &
                      (/ d18O(tag_IO_CM), C2_d18O(cn), C2_d18Oe(cn) /) /)
#endif
#ifdef iCMCb_no_d13C_match
        if ( iCMCb_check_pairs(2, (/ (/ C(ind_CM),       C2_CO2(cn),  2.*C2_CO2e(cn)  /) &
                                   , (/ d18O(tag_IO_CM), C2_d18O(cn), 2.*C2_d18Oe(cn) /) &
#else
        if ( iCMCb_check_pairs(3, (/ (/ C(ind_CM),       C2_CO2(cn),  2.*C2_CO2e(cn)  /) &
                                   , (/ d18O(tag_IO_CM), C2_d18O(cn), 2.*C2_d18Oe(cn) /) &
                                   , (/ d13C(tag_IC_CM), C2_d13C(cn), 2.*C2_d13Ce(cn) /) &
#endif
                                   /) ) ) then
          hits_delta(cn,iCMCb_set_IC) = 1
        endif

      ! CO/O3 + GHGs
#ifdef iCMCb_no_SF6_match
        if ( iCMCb_check_pairs(3, (/ (/ C(ind_CM_CH4), C2_CH4(cn), 1.*C2_CH4e(cn) /) &
                                   , (/ C(ind_CM_N2O), C2_N2O(cn), 1.*C2_N2Oe(cn) /) &
                                   , (/ C(ind_CM_CO),  C2_CO(cn),  1.*C2_COe(cn)  /) &
#else
        if ( iCMCb_check_pairs(4, (/ (/ C(ind_CM_CH4), C2_CH4(cn), 1.*C2_CH4e(cn) /) &
                                   , (/ C(ind_CM_N2O), C2_N2O(cn), 1.*C2_N2Oe(cn) /) &
                                   , (/ C(ind_CM_CO),  C2_CO(cn),  1.*C2_COe(cn)  /) &
                                   , (/ C(ind_CM_SF6), C2_SF6(cn), 1.*C2_SF6e(cn) /) &
#endif
                                   /) ) ) then
          hits_delta(cn,iCMCb_set_TG) = 1

        ! full match?
          if ( hits_delta(cn,iCMCb_set_IC) .gt. 0 ) then
            hits_delta(cn,iCMCb_set_FM) = 1
          endif

        endif

      enddo

#ifdef iCMCb_show_match_progress
    ! updating totals
      do sn = 0, iCMCb_hit_sets
        hits_delta( 0, sn) = sum(hits_delta(1:iCMCb_hit_cases, sn))
      enddo

    ! progress indicator
      if     ( hits_delta(0,iCMCb_set_FM) .gt. 0 ) then
        write(*,'(A)',advance="no"),'*'
      elseif ( hits_delta(0,iCMCb_set_IC) .gt. 0 ) then
        write(*,'(A)',advance="no"),'i'
      elseif ( hits_delta(0,iCMCb_set_TG) .gt. 0 ) then
        write(*,'(A)',advance="no"),'t'
      else
        write(*,'(A)',advance="no"),'.'
      endif
    else

    ! SS was not reached
      write(*,'(A)',advance="no"),'-'
#endif

    endif

  ! all tries
  ! if ( iCMCb_real .le. iCMCb_req_hits ) call iCMCb_output(-1,2,iCMCb_real,s)

  end subroutine iCMCb_find_matching_cases


! ----- sync hit statistic and output results ----------------------------------
  subroutine iCMCb_sync_n_output

    use caaba_io
    use mo_mpi

    implicit none

    logical      :: incoming
    integer      :: sn, cn                    !> set and case nos.
#ifndef NOMPI
    integer      :: p_status(MPI_STATUS_SIZE) !> standard information of MPI_RECV
    integer      :: p_source                  !> slave no.

  ! ----- sync in case of parallel run -----
  ! slaves first needs to sync hit nos. with others
    if ( .not.p_parallel_io ) then

      if ( sum(hits_delta(:,:)) .eq. 0 ) return          ! nothing found? not sending data

      call iCMCb_prep_output                             ! fill-up the output buffer
      call p_isend(hits_delta, p_io, iCMCb_mpitag_hits)  ! sending delta-hits
      call p_isend(iCMCb_out, p_io, iCMCb_mpitag_out)    ! sending the output

      return                                             ! leaving

    else

    ! master first processes slaves
      hits_tmp = hits_delta   ! saving master hits (received hits/data are output in a similar fashion)

    ! listening to incoming results from slaves
      do
        incoming = .false.
        call MPI_iprobe(MPI_ANY_SOURCE, iCMCb_mpitag_hits, p_all_comm, &
                        incoming, p_status, p_error)
#ifdef DEBUG
        if (p_error/=MPI_SUCCESS) then
          write (nerr,'(a,i4,a,i4,a)') ' MPI_iprobe on pe #', p_pe, &
                      ' for tag ', iCMCb_mpitag_hits, ' failed.'
          write (nerr,'(a,i4)') ' Error = ', p_error
          call p_abort
        endif
#endif
      ! there is an update from a slave
        if (incoming) then
          p_source = p_status(MPI_SOURCE)
          hits_tmp(:,:) = 0
        ! receiving delta-hits from slave
          call p_recv(hits_delta, p_source, iCMCb_mpitag_hits)
        ! and output data
          call p_recv(iCMCb_out, p_source, iCMCb_mpitag_out)
        ! output & update hits
          call iCMCb_output_local

        else
          exit  ! repeating until all slaves are updated (no incoming requests)
        endif
      enddo

    ! restoring master's hits
      hits_delta = hits_tmp

    endif
#endif

  ! this is for master's results
    if ( sum(hits_delta(:,:)) .eq. 0 ) return   ! if no output

  ! fill-up output buffer
    call iCMCb_prep_output
  ! output & update hits
    call iCMCb_output_local


  contains

  ! ----- local output + hits update -----
    subroutine iCMCb_output_local

      implicit none

    ! augment the hits
      do sn = 1, iCMCb_hit_sets
        do cn = C2_match_from, C2_match_upto   ! TODO: -> 1, iCMCb_hit_cases
          if ( hits_delta(cn,sn) .gt. 0 ) then
          ! updating hit counters for the set/case
            iCMCb_hits(cn,sn) = iCMCb_hits(cn,sn) + hits_delta(cn,sn)
          ! updating output data
            iCMCb_out(6) = iCMCb_hits(cn,sn)  ! hit no.
            iCMCb_out(7) = cn                 ! hit sample no.
          ! output to case-container
            if ( iCMCb_hits(cn,sn) .le. iCMCb_max_hits ) &
              call iCMCb_output_buf(cn,sn)
          ! and to set-common container
            if ( iCMCb_hits(cn,sn) .le. iCMCb_req_hits ) then
              iCMCb_hits(0,sn) = iCMCb_hits(0,sn) + hits_delta(cn,sn)   ! totals
              call iCMCb_output_buf(0,sn)
            endif
          endif
        enddo
      enddo

    ! full-match samples: output only to common cont.
      do cn = C2_match_from, C2_match_upto   ! TODO: -> 1, iCMCb_hit_cases
        if ( hits_delta(cn,iCMCb_set_FM) .gt. 0 ) then
          iCMCb_hits(cn,iCMCb_set_FM) = iCMCb_hits(cn,iCMCb_set_FM) + hits_delta(cn,iCMCb_set_FM)
          iCMCb_hits( 0,iCMCb_set_FM) = iCMCb_hits( 0,iCMCb_set_FM) + hits_delta(cn,iCMCb_set_FM)
          iCMCb_out(6) = iCMCb_hits(cn,iCMCb_set_FM) ! hit no.
          iCMCb_out(7) = cn                          ! hit sample no.
          call iCMCb_output_buf(0,iCMCb_set_FM)
        endif
      enddo

    end subroutine iCMCb_output_local

    subroutine iCMCb_output_buf(cn, sn)
      implicit none
      integer, intent(in) :: cn, sn
      real(dp) :: out_time

      out_time = iCMCb_time()
#ifdef iCMCb_netcdf
#ifdef DEBUG
      print *,'iCMCb_output: out_time = ',out_time,' p> sn: ',sn,'s> cn:',cn,'c> rec:',iCMCb_hits(cn,sn)
      iCMCb_out(6) = iCMCb_hits(cn,sn)  ! recno
#endif
!if (cn.eq.0) print *,'&&& out_time = ',out_time,' p> sn: ',sn,'s> cn:',cn,'c> rec:',iCMCb_hits(cn,sn)
      call write_output_file(iCMCb_iounit(cn,sn), out_time, iCMCb_out, master=p_parallel_io) !, recno=iCMCb_hits(cn,sn))
#else
     !open(unit=iCMCb_iounit(caseno,setno), file=fname//'.dat', status='replace') !, share='denynone')
     !write...
     !close...
     !...write(str,"(3(I0,'	'),6(F0.5,'	'),18(F0.5,'	'),12(F0.5,'	'))"), &
     !   iCMCb_tries, iCMCb_real, iCMCb_harvest, & ...
#endif

    end subroutine iCMCb_output_buf

  end subroutine iCMCb_sync_n_output


! ----- main loop --------------------------------------------------------------
  subroutine iCMCb_physc

    use caaba_mem,              only: model_time, model_start, model_end, timesteplen
    use mo_mpi
    use messy_mecca_tag_IC_box
    use messy_mecca_tag_IO_box

    implicit none
    intrinsic mod

    integer            :: seqmatch, cn
    integer, parameter :: reqseqmatch = 5
    integer, save      :: stat_quantum
    logical            :: req_hits_found4all = .false.
#ifndef NOMPI
    integer            :: p_status(MPI_STATUS_SIZE) !> standard information of MPI_RECV
#endif

  ! first try
    if ( iCMCb_tries .eq. 0 ) then
      call iCMCb_new_x0
      iCMCb_tries = 1
      stat_quantum = 100
#ifdef DEEPDEBUG
      print *, 'iCMCb_physc: 1st init'
      call tag_IC_set( ind_RC, C(ind_RC), (/ d13C(tag_IC_RC) /) )
      print *, 'd13C(tag_IC_RC) = ',d13C(tag_IC_RC)
#endif
    endif


  ! nudging plant/respired CO2 d18O/d13C
    call tag_IO_set(ind_PC, C(ind_PC), (/ d18O_PC, delta17Opm( d18O_PC, 0._dp ) /) )  ! D17O is SET to 0 (MWL?)
    call tag_IC_set(ind_RC, C(ind_RC), (/ d13C_RC /) )
#ifdef DEEPDEBUG
    write(*, '(2(A,F10.3,A,F10.3))'), 'iCMCb_physc: nudging curr. -> to / d13C(tag_IC_RC) = ', &
      d13C(tag_IC_RC), '->', d13C_RC, ' / d18O(tag_IO_PC) = ',d18O(tag_IO_PC), '->', d18O_PC
#endif

  ! checking for steady-state: UA = CM ?
    if ( iCMCb_check_pairs(3, (/ (/ C(ind_UA),       C(ind_CM),       1.0e-09_dp*cair /), &
                                 (/ d13C(tag_IC_UA), d13C(tag_IC_CM), 0.001_dp /), &
                                 (/ d18O(tag_IO_UA), d18O(tag_IO_CM), 0.001_dp /) /) ) ) then
      iCMCb_ss_seqmatches = iCMCb_ss_seqmatches + 1
    else
      iCMCb_ss_seqmatches = 0
    endif

  ! SS reached?
    iCMCb_is_SS = ( iCMCb_ss_seqmatches.ge.reqseqmatch )

  ! SS/end situation
    if ( iCMCb_is_SS .or. (model_time.ge.model_end) ) then

      if ( iCMCb_is_SS ) then             ! SS/end reached, checking match, starting new realisation
#ifdef DEEPDEBUG
        print *,'iCMCb_physc: SS reached for #',iCMCb_tries
#endif
        iCMCb_real = iCMCb_real + 1  ! next SS realisation
      else
#ifdef DEEPDEBUG
        print *,'iCMCb_physc: end reached for #',iCMCb_tries
#else

#endif
      endif

    ! check for matches
      call iCMCb_find_matching_cases

    ! calling sync and output via master
      call iCMCb_sync_n_output

    ! print stats occasionally
!     if ( mod(iCMCb_tries, 120 + 0*(iCMCb_req_hits/100) ) .eq. 0 ) then
      if ( mod(iCMCb_tries, stat_quantum) .eq. 0 ) then
        call iCMCb_print_stats
        if ( iCMCb_tries .gt. 2000 ) stat_quantum = 1000
      endif

      if ( p_parallel_io ) then
      ! check for req. statistic
        req_hits_found4all = .true.

        !if ( iCMCb_hit(0,iCMCb_set_IC) .gt. ( C2_match_upto - C2_match_from + 1 ) * iCMCb_req_hits ) &
        do cn = C2_match_from, C2_match_upto
          if ( iCMCb_hits(cn,iCMCb_set_IC) .lt. iCMCb_req_hits ) &
             req_hits_found4all = .false.
        enddo
        if ( sum(iCMCb_hits(1:iCMCb_hit_cases,iCMCb_set_TG)) .lt. iCMCb_req_hits * ( C2_match_upto - C2_match_from + 1 ) ) &
           req_hits_found4all = .false.

#ifndef NOMPI
        if ( req_hits_found4all ) then  ! sending "done" signal to slaves
          do cn = 1, p_nprocs-1
            call p_send(req_hits_found4all, cn, iCMCb_mpitag_done)
          enddo
        endif
      else
      ! check for "done" signal from master
        req_hits_found4all = .false.
        call MPI_iprobe(MPI_ANY_SOURCE, iCMCb_mpitag_done, p_all_comm, &
                        req_hits_found4all, p_status, p_error)
#endif
      endif

      if ( req_hits_found4all ) then 
      ! done, triggering exit from the outer (caaba's) loop
        model_time = model_end
      else
      ! prepare new try
        iCMCb_tries = iCMCb_tries + 1     ! new try no.
        iCMCb_ss_seqmatches = 0           ! reset match counter
        model_time = model_start          ! reset time for CAABA
        call iCMCb_new_x0                 ! get new initial cond.
      endif

    endif ! integrating further

  end subroutine iCMCb_physc

#endif



#ifndef DONOTINCLUDETHIS
! checks if file exists
  logical function file_exists(fname)
    character(*), intent(in) :: fname
    logical res
    inquire(file=fname, exist=res)
    file_exists = res
  end function file_exists

! better randomises seed
  subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t

            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed

! ########################################################################################################33
! EOF iCMCb
! ########################################################################################################33
#endif









#ifdef OFFLEMB
!*******************************************************************************

! copied from old mo_netcdf module from MECCA releases for compatibility reasons

  subroutine nf(status) ! turns nf90_* function into subroutine + checks status
    integer :: status
    if (status /= nf90_noerr) then
      write (*,*) 'netcdf error: ', nf90_strerror(status)
      stop
    endif
  end subroutine nf


! = pseudo-mixing ==============================================================

  subroutine pmix(TSL, dilF, ind, mix2_amount)

  ! pseudo-mixing of species ind_d with background concentration mix2_amount
  ! within TSL timestep with dilF dilution factor [1/s]

    use caaba_mem, only: C

    implicit none

    real(dp), intent(in) :: TSL, dilF
    integer, intent(in)  :: ind
    real(dp), intent(in) :: mix2_amount

  ! pseudo-mixing of ind
    C(ind) = C(ind) + ( mix2_amount - C(ind) ) * TSL * dilF

  end subroutine pmix

! = nudging ====================================================================

  subroutine nudge_init(filedata_nc)  ! , clat, clon)

    implicit none

    character(*), intent(in)      :: filedata_nc
!    real(dp), intent(out)         :: clat, clon
    integer                       :: varid

    ! open nc with eval19 data for reading
    print *, 'nudge_init: ouvre de fichier: '//filedata_nc//'.nc'
    call nf(nf90_open(filedata_nc//'.nc', nf90_nowrite, ncid_nudge))

    ! fixer les coordonnees geographiques,

    ! latitude
!    call nf(nf90_inq_varid(ncid_nudge, "CLAT", varid))
!    call nf(nf90_get_var(ncid_nudge, varid, clat, start = (/1/)))

    ! longitude
!    call nf(nf90_inq_varid(ncid_nudge, "CLON", varid))
!    call nf(nf90_get_var(ncid_nudge, varid, clon, start = (/1/)))

!    print *, ' nouvelles cordonnees geographiques: ',clat,'N ',clon,' E'

  end subroutine nudge_init

  !-----------------------------------------------------------------------------

  subroutine nudge_close

    implicit none

    call nf(nf90_close(ncid_nudge))

  end subroutine nudge_close

  !-----------------------------------------------------------------------------

  subroutine nudge_spec(model_time, rel2nc_start, model2nc_scale, spec_list_ind, nudc_list)

    use messy_mecca_kpp_monitor, only: SPC_NAMES
    use messy_main_constants_mem, only: R_gas, N_A, STRLEN_SHORT, STRLEN_MEDIUM
    use netcdf

    use caaba_mem, only: C, temp, cair, press

    implicit none

    integer, intent(in)           :: spec_list_ind(:)      ! liste des substances
    real(dp), intent(in)          :: nudc_list(:)          ! avec les coefficients de nudge

    real(dp), intent(in)          :: model_time, rel2nc_start, model2nc_scale
!                                                              ^ how many model steps fits into 1 nc step
!                                                                for monthly data should be = 365.24 / 12, z.B.
    real(dp)                      :: frac, dum
    integer                       :: timestep_nc, s


  ! calculating corresponding time frame in nc and fraction for linear interpolation
    timestep_nc = 1 + int( ( model_time - rel2nc_start ) / model2nc_scale )
    frac = model_time / model2nc_scale - int( model_time / model2nc_scale )

!    print *, 'nudge_spec (',model_time,',',rel2nc_start,' (',model2nc_scale,'): nc pas # ',timestep_nc,', frac = ',frac

    ! prendre la temperature
    dum = get_val_i("TM1")
    if ( dum .lt. 400 ) temp = get_val_i("TM1")      ! sorting UNDEF through dum...
    cair = (N_A/1.E6) * press / (R_gas*temp) ! cair = c(air) in [mcl/cc]

#ifdef DEBUG
    print *,'nudging the temperature: ',temp,' K'
#endif

    ! et le reste des substances
    do s = 1, ubound(spec_list_ind,1)
#ifdef DEBUG
      print *,'nudging ',trim(SPC_NAMES(spec_list_ind(s))),' (',nudc_list(s),'): ', &
        C(spec_list_ind(s))/cair*1E9,' -> ', &
        get_val_i(SPC_NAMES(spec_list_ind(s)))*1E9,' => ',&
        ( 1.0_dp - nudc_list(s) ) * C(spec_list_ind(s))/cair*1e9 + &
          nudc_list(s) * ( get_val_i(SPC_NAMES(spec_list_ind(s)))*1e9 )
#endif
      C(spec_list_ind(s)) = ( 1.0_dp - nudc_list(s) ) * C(spec_list_ind(s)) + &
                            nudc_list(s) * ( get_val_i(SPC_NAMES(spec_list_ind(s))) * cair )
    enddo

! spec_list_ind:
!
! (/ ind_CH4, ind_CO, ind_CO2, ind_HCHO, ind_OH, ind_HO2, &
!    ind_NO, ind_NO2, ind_NO3, ind_N2O, ind_N2O5, ind_O3, &
!    ind_CL2, ind_SO2, ind_CH3O2, ind_CH3CHO, ind_CH3OH, ind_CH3OOH, &
!    ind_C2H4, ind_C2H6, ind_C3H6, ind_C3H8, ind_NC4H10, &
!    ind_ISO2, ind_C5H8, ind_PA, ind_PAN, ind_MVKO2, ind_LMEKO2 /)

  contains

    real(dp) function get_val_i(varname)

      implicit none

      character(*), intent(in)    :: varname
      integer                     :: varid
      real(dp)                    :: val_ava, val_sui

      if ( nf90_inq_varid(ncid_nudge, varname, varid) == nf90_noerr ) then
        call nf(nf90_get_var(ncid_nudge, varid, val_ava, start = (/timestep_nc/)))
        call nf(nf90_get_var(ncid_nudge, varid, val_sui, start = (/timestep_nc+1/)))
        ! l'iterpolation simple (linéaire)
        get_val_i = val_ava + frac * (val_sui - val_ava)
      else
        ! variable not found
!        print *, "get_val_i(",varname,"): not found"
        get_val_i = 0.0_dp    ! just in case
      endif

    end function get_val_i

  end subroutine nudge_spec

! = offline emission from E5/M1 data ===========================================

  subroutine offlemb_init(filedata_nc)

    implicit none

    character(*), intent(in)      :: filedata_nc
    integer                       :: varid

  ! open nc with eval19 offlem data for reading
    print *, 'offlemb_init: ouvre de fichier: '//filedata_nc//'.nc'
    call nf(nf90_open(filedata_nc//'.nc', nf90_nowrite, ncid_offlemb))

  end subroutine offlemb_init

  !-----------------------------------------------------------------------------

  subroutine offlemb_close

    implicit none

    call nf(nf90_close(ncid_offlemb))

  end subroutine offlemb_close

  !-----------------------------------------------------------------------------

! performs offline emission based on E5/M1 eval19 offlem setup
! using data prepared by icogesb  scripts

! tagging configurations switches

  subroutine offlemb_perform(model_time, rel2nc_start, model2nc_scale, fct)

    use messy_mecca_kpp_monitor, only: SPC_NAMES
    use caaba_mem,               only: C

    use netcdf

#ifdef tag_IC
    use messy_mecca_tag_IC_box
#endif
#ifdef tag_FCF
    use messy_mecca_tag_FCF_box
#endif
#ifdef tag_FCB
    use messy_mecca_tag_FCB_box
#endif
#ifdef tag_FO17v2
    use messy_mecca_tag_FO17v2_box
#endif


#ifdef tag_IC
    use messy_mecca_tag_IC_box
#endif
#ifdef tag_IO
    use messy_mecca_tag_IO_box
#endif
#ifdef tag_O3F
    use messy_mecca_tag_O3F_box
#endif


    implicit none

    real(dp), intent(in)   :: model_time, rel2nc_start, model2nc_scale, fct
!                                                             ^ how many model steps fits into 1 nc step
!                                                               for monthly data should be = 365.24 / 12, z.B.

    integer                :: varid, timestep_nc
    real(dp)               :: val_ava, val_sui, val_cur, val_tot, frac
    character(len=12)      :: spcname, req

  ! exp "undefined" flag
    real(dp), parameter    :: eUNDEF = -1E33_dp

  ! - data section -------------------------------------------------------------

  ! specs list   (to search for in emission file)
    integer            :: iSIL
    integer, parameter :: nSIL = 17
    integer            :: SIL(nSIL) = &
      (/ ind_CH4, &
         ind_CO, &
         ind_HCHO, &
         ind_CH3OH, &
         ind_HCOOH, &
         ind_C2H4, &
         ind_C2H6, &
         ind_C3H6, &
         ind_C3H8, &
         ind_NC4H10, &
         ind_CH3CHO, &
         ind_CH3COCH3, &
         ind_CH3CO2H, &
         ind_MEK, &
         ind_NO, &
         ind_SO2, &
         ind_C5H8 /)

  ! Dummy is instead of CO

  ! emission classes
    integer            :: iOEC
    integer, parameter :: nOEC = 7
    character(len=10)  :: OEC(nOEC) = (/ "BB   ", &
                                         "BF   ", &
                                         "FF   ", &
                                         "L43  ", &
                                         "LAND ", &
                                         "OCE  ", &
                                         "SHIPS" /)

  ! warn: xBB turns BB emissions off

  !        BB          BF        FF         L43        LAND       OCE      SHIPS

  ! corresponding 13C signatures for classes emissions
  ! area: highNH
    real(dp), parameter :: S13C(nOEC,nSIL) = RESHAPE( &
       (/ &
          eUNDEF, -35.00_dp,    eUNDEF, -53.00_dp,    eUNDEF,    eUNDEF, -27.50_dp,  &
       -24.50_dp, -27.50_dp, -27.50_dp, -27.50_dp, -27.00_dp, -13.50_dp, -27.50_dp,  &
       -24.50_dp, -27.50_dp, -27.50_dp, -27.50_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &
       -24.50_dp, -27.50_dp, -27.50_dp, -27.50_dp, -27.00_dp,    eUNDEF,    eUNDEF,  &
       -24.50_dp, -27.50_dp, -27.50_dp, -27.50_dp, -27.00_dp,    eUNDEF,    eUNDEF,  &
       -24.50_dp, -27.50_dp, -22.20_dp, -27.50_dp, -27.00_dp, -20.00_dp, -27.50_dp,  &
       -24.50_dp, -27.50_dp, -27.40_dp, -27.50_dp,    eUNDEF, -20.00_dp, -27.50_dp,  &
       -24.50_dp, -27.50_dp, -25.20_dp, -27.50_dp, -27.00_dp,    eUNDEF, -27.50_dp,  &
       -24.50_dp, -27.50_dp, -27.70_dp, -27.50_dp,    eUNDEF, -20.00_dp, -27.50_dp,  &
       -24.50_dp, -27.50_dp, -30.60_dp, -27.50_dp,    eUNDEF, -20.00_dp, -27.50_dp,  &
       -24.50_dp, -27.50_dp, -27.50_dp, -27.50_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &
       -24.50_dp, -27.50_dp, -27.50_dp, -27.50_dp, -27.00_dp,    eUNDEF,    eUNDEF,  &
       -24.50_dp, -27.50_dp, -27.50_dp, -27.50_dp, -27.00_dp,    eUNDEF,    eUNDEF,  &
       -24.50_dp, -27.50_dp, -27.50_dp, -27.50_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &
          eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,  &
          eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,  &
          eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF, -27.50_dp,    eUNDEF,    eUNDEF   &
        /), (/ nOEC, nSIL /) )

  ! corresponding 18O signatures for classes emissions
  ! area: highNH
    real(dp), parameter :: S18O(nOEC,nSIL) = RESHAPE( &
       (/ &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,  &         ! ind_CH4, &
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,  -8.00_dp,   0.00_dp,  23.50_dp,   &        ! ind_CO, &
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_HCHO, &
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,  -8.00_dp,    eUNDEF,    eUNDEF,   &        ! ind_CH3OH, &
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,  -8.00_dp,    eUNDEF,    eUNDEF,   &        ! ind_HCOOH, &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_C2H4, &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_C2H6, &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_C3H6, &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_C3H8, &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_NC4H10, &
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &         ! ind_CH3CHO,
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,  -8.00_dp,    eUNDEF,    eUNDEF,  &         ! ind_CH3COCH3
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,  -8.00_dp,    eUNDEF,    eUNDEF,  &         ! ind_CH3CO2H,
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &         ! ind_MEK, &
       17.15_dp,  17.20_dp,  23.50_dp,  17.20_dp,   0.00_dp,    eUNDEF,  23.50_dp,  &         ! ind_NO, &
         eUNDEF,  17.20_dp,  23.50_dp,  17.20_dp,    eUNDEF,    eUNDEF,  23.50_dp,  &         ! ind_SO2, &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF   &         ! ind_C5H8 /)
        /), (/ nOEC, nSIL /) )

!#  ! corresponding 18O signatures for classes emissions
!#  ! area: highNH
!#    real(dp), parameter :: S18O(nOEC,nSIL) = RESHAPE( &
!#       (/ &
!#         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,  &         ! ind_CH4, &
!#       17.15_dp,   0.00_dp,  23.50_dp,  17.20_dp,  -8.00_dp,   0.00_dp,  23.50_dp,   &        ! ind_CO, &
!#       17.15_dp,   0.00_dp,  23.50_dp,  17.20_dp,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_HCHO, &
!#       17.15_dp,   0.00_dp,  23.50_dp,  17.20_dp,  -8.00_dp,    eUNDEF,    eUNDEF,   &        ! ind_CH3OH, &
!#       17.15_dp,   0.00_dp,  23.50_dp,  17.20_dp,  -8.00_dp,    eUNDEF,    eUNDEF,   &        ! ind_HCOOH, &
!#         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_C2H4, &
!#         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_C2H6, &
!#         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_C3H6, &
!#         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_C3H8, &
!#         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &        ! ind_NC4H10, &
!#       17.15_dp,   0.00_dp,  23.50_dp,  17.20_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &         ! ind_CH3CHO,
!#       17.15_dp,   0.00_dp,  23.50_dp,  17.20_dp,  -8.00_dp,    eUNDEF,    eUNDEF,  &         ! ind_CH3COCH3
!#       17.15_dp,   0.00_dp,  23.50_dp,  17.20_dp,  -8.00_dp,    eUNDEF,    eUNDEF,  &         ! ind_CH3CO2H,
!#       17.15_dp,   0.00_dp,  23.50_dp,  17.20_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &         ! ind_MEK, &
!#       17.15_dp,  23.50_dp,  23.50_dp,  17.20_dp,   0.00_dp,    eUNDEF,  23.50_dp,  &         ! ind_NO, &
!#         eUNDEF,  23.50_dp,  23.50_dp,  17.20_dp,    eUNDEF,    eUNDEF,  23.50_dp,  &         ! ind_SO2, &
!#         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF   &         ! ind_C5H8 /)
!#        /), (/ nOEC, nSIL /) )

    real(dp), parameter :: S17Ocap(nOEC,nSIL) = RESHAPE( &
       (/ &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &
         0.0_dp,   0.00_dp,    0.0_dp,    0.00_dp,  0.00_dp,   0.00_dp,    0.0_dp,   &
         0.0_dp,   0.00_dp,    0.0_dp,    0.00_dp,   eUNDEF,    eUNDEF,    eUNDEF,   &
         0.0_dp,   0.00_dp,    0.0_dp,    0.00_dp,  0.00_dp,    eUNDEF,    eUNDEF,   &
         0.0_dp,   0.00_dp,    0.0_dp,    0.00_dp,  0.00_dp,    eUNDEF,    eUNDEF,   &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,   &
         0.0_dp,   0.00_dp,    0.0_dp,   0.00_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &
         0.0_dp,   0.00_dp,    0.0_dp,   0.00_dp,   0.00_dp,    eUNDEF,    eUNDEF,  &
         0.0_dp,   0.00_dp,    0.0_dp,   0.00_dp,   0.00_dp,    eUNDEF,    eUNDEF,  &
         0.0_dp,   0.00_dp,    0.0_dp,   0.00_dp,    eUNDEF,    eUNDEF,    eUNDEF,  &
         0.0_dp,   0.00_dp,    0.0_dp,    0.0_dp,   0.00_dp,    eUNDEF,    0.0_dp,  &
         eUNDEF,   0.00_dp,    0.0_dp,    0.0_dp,    eUNDEF,    eUNDEF,    0.0_dp,  &
         eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF,    eUNDEF   &
        /), (/ nOEC, nSIL /) )


  ! - code section -------------------------------------------------------------

  ! calculating corresponding time frame in nc and fraction for linear interpolation
    timestep_nc = 1 + int( ( model_time - rel2nc_start ) / model2nc_scale )
    frac = model_time / model2nc_scale - int( model_time / model2nc_scale )

#ifdef DEBUG
  ! print *, 'offlemb (',model_time,',',rel2nc_start,' (',model2nc_scale,'): nc pas # ',timestep_nc,', frac = ',frac
#endif

  ! species cycle
    do iSIL = 1, nSIL

      if ( SIL(iSIL) .eq. 0 ) cycle

      val_tot = 0.0_dp       ! total emission of the species (inter-class sum)

    ! classes cycle
      do iOEC = 1, nOEC

      ! acquiring variable with name SPEC_CLASS, z.B.  CH4_BB,
      ! if it exists in the datafile, processing

        spcname = trim(SPC_NAMES(SIL(iSIL)))

        if (spcname == "NO") spcname = "NOX"        ! for NO we're searching for NOX emission record

        if ( nf90_inq_varid(ncid_offlemb, trim(spcname)//"_"//OEC(iOEC), varid) == nf90_noerr ) then

          call nf(nf90_get_var(ncid_offlemb, varid, val_ava, start = (/ cyc(timestep_nc,   1, 12, 12) /)))
          call nf(nf90_get_var(ncid_offlemb, varid, val_sui, start = (/ cyc(timestep_nc+1, 1, 12, 12) /)))

        ! l'iterpolation simple (linéaire)
          val_cur = val_ava + frac * (val_sui - val_ava)

 !         if (spcname == "C5H8") val_cur = val_cur * 20        ! isop emission scale
!          if (spcname == "NOX") val_cur = val_cur * 0        ! isop emission scale
        ! val_cur = val_cur / 9.448_dp        ! scaling emissions to a factor of 9.448 (NH versus remNH in CO)
          val_cur = val_cur / 12_dp       ! scaling emissions for 70N
!          val_cur = val_cur / 7_dp        ! scaling emissions for 30N

          val_tot = val_tot + val_cur

#ifdef DEBUG
          print *, trim(spcname)//"_"//OEC(iOEC),': (',varid,') ',val_sui,' <-> ',val_ava,' >> ',val_cur
#endif

        ! TODO: call here emission according to the class

#ifdef tag_IC
        ! checking if the value for the class is assigned in input data
          if (S13C(iOEC,iSIL) /=eUNDEF) then

          ! emitting 12C/13C mixture converted from m2 to cm2
            call tag_IC_emis(SIL(iSIL), val_cur * fct / 1.0E4_dp, (/ S13C(iOEC,iSIL) /) )
!            print *,'#tag_IC#: ',trim(SPC_NAMES(SIL(iSIL))), ' > ',trim(OEC(iOEC)), ':', S13C(iOEC,iSIL)

          else
!            print *,'offlemb_perform: #tag_IC# >eUNDEF signature for', &
!                                                        ' spec: ',trim(SPC_NAMES(SIL(iSIL))), &
!                                                       ' class: ',trim(OEC(iOEC))
          endif
#endif
#ifdef tag_FCF
        ! emission class FF - n=3
          if (iOEC == 3) then
          ! checking if the value for the class is assigned in input data
            if (S13C(iOEC,iSIL) /= eUNDEF) then

            ! emitting FCF species converted from m2 to cm2
            call tag_FCF_emis(SIL(iSIL), val_cur * fct / 1.0E4_dp, (/ 1.0_dp /) )
!              print *,'#tag_FCF#: ',trim(SPC_NAMES(SIL(iSIL))), ' > ',trim(OEC(iOEC)), ':', val_cur * fct / 1.0E4_dp

            else
!              print *,'offlemb_perform: #tag_FCF# >eUNDEF signature for', &
!                                                          ' spec: ',trim(SPC_NAMES(SIL(iSIL))), &
!                                                         ' class: ',trim(OEC(iOEC))
            endif
          endif
#endif
#ifdef tag_IC
        ! checking if the value for the class is assigned in input data
          if (S13C(iOEC,iSIL) /=eUNDEF) then

          ! emitting 12C/13C mixture converted from m2 to cm2
            call tag_IC_emis(SIL(iSIL), val_cur * fct / 1.0E4_dp, (/ S13C(iOEC,iSIL) /) )
!            print *,'#tag_IC#: ',trim(SPC_NAMES(SIL(iSIL))), ' > ',trim(OEC(iOEC)), ':', S13C(iOEC,iSIL)

          else
!            print *,'offlemb_perform: #tag_IC# >eUNDEF signature for', &
!                                                        ' spec: ',trim(SPC_NAMES(SIL(iSIL))), &
!                                                       ' class: ',trim(OEC(iOEC))
          endif
#endif

#ifdef tag_IO
        ! checking if the value for the class is assigned in input data
          if (S18O(iOEC,iSIL) /=eUNDEF) then

          ! emitting 18O/17O/16O mixture converted from m2 to cm2
            call tag_IO_emis(SIL(iSIL), val_cur * fct / 1.0E4_dp, &
                             (/ S17Ocap(iOEC,iSIL) + MDFSL_O * S18O(iOEC,iSIL), S18O(iOEC,iSIL) /) )
!            print *,'#tag_IO#: ',trim(SPC_NAMES(SIL(iSIL))), ' > ',trim(OEC(iOEC)), ':', S18O(iOEC,iSIL),' / ',S17Ocap(iOEC,iSIL)

          else
!            print *,'offlemb_perform: #tag_IO# >eUNDEF signature for', &
!                                                        ' spec: ',trim(SPC_NAMES(SIL(iSIL))), &
!                                                       ' class: ',trim(OEC(iOEC))
          endif
#endif
#ifdef tag_IO
        ! checking if the value for the class is assigned in input data
          if (S18O(iOEC,iSIL) /=eUNDEF) then

          ! emitting 18O/17O/16O mixture converted from m2 to cm2
            call tag_IO_emis(SIL(iSIL), val_cur * fct / 1.0E4_dp, &
                             (/ S17Ocap(iOEC,iSIL) + MDFSL_O * S18O(iOEC,iSIL), S18O(iOEC,iSIL) /) )
!            print *,'#tag_IO#: ',trim(SPC_NAMES(SIL(iSIL))), ' > ',trim(OEC(iOEC)), ':', S18O(iOEC,iSIL),' / ',S17Ocap(iOEC,iSIL)

          else
!            print *,'offlemb_perform: #tag_IO# >eUNDEF signature for', &
!                                                        ' spec: ',trim(SPC_NAMES(SIL(iSIL))), &
!                                                       ' class: ',trim(OEC(iOEC))
          endif
#endif
#ifdef tag_FO17v2
        ! checking if the value for the class is assigned in input data
          if (S18O(iOEC,iSIL) /= eUNDEF) then

          ! emitting 18O/17O/16O mixture converted from m2 to cm2
            call tag_FO17v2_emis( SIL(iSIL), val_cur * fct / 1.0E4_dp , (/ 1.0_dp /) )
#ifdef DEBUG
            print *,'#tag_FO17v2#: ',trim(SPC_NAMES(SIL(iSIL))), ' > ',trim(OEC(iOEC)), ':', val_cur
#endif
          else
#ifdef DEBUG
            print *,'offlemb_perform: #tag_IO# >eUNDEF signature for', &
                                                        ' spec: ',trim(SPC_NAMES(SIL(iSIL))), &
                                                       ' class: ',trim(OEC(iOEC))
#endif
          endif
#endif

        endif

      enddo

  ! emission of the regular species
    if (SIL(iSIL) /= 0) C(SIL(iSIL)) = C(SIL(iSIL)) + val_tot * fct / 1.0E4_dp   ! converting from m2 to cm2
#ifdef DEBUG
    if (SIL(iSIL) /= 0) print *, spcname, ':', C(SIL(iSIL)), ' + ', val_tot   ! converting from m2 to cm2
#endif
    enddo

!    ! et le reste des substances
!    do s = 1, ubound(spec_list_ind,1)
!      C(spec_list_ind(s)) = ( 1.0_dp - nudc ) * C(spec_list_ind(s)) + &
!                            nudc * ( get_val_i(SPC_NAMES(spec_list_ind(s))) * cair )
!    enddo

  contains

  ! cycling function between 1 and given parameter
    integer function cyc(cinp, cmin, cmax, cper)

       implicit none

       integer, intent(in) :: cinp, cmin, cmax, cper
       integer             :: cn

       cn = cinp; do while (cn .lt. cmin); cn = cn + cper; enddo
       do while (cn .gt. cmax); cn = cn - cper; enddo; cyc = cn;

     end function cyc

!
!    real(dp) function get_val_i(varname)
!
!      implicit none
!
!      character(*), intent(in)    :: varname
!      integer                     :: varid
!      real(dp)                    :: val_ava, val_sui
!
!      call nf(nf90_inq_varid(ncid_offlemb, varname, varid))
!      call nf(nf90_get_var(ncid_offlemb, varid, val_ava, start = (/timestep_nc/)))
!      call nf(nf90_get_var(ncid_offlemb, varid, val_sui, start = (/timestep_nc+1/)))
!
!
!    end function get_val_i
!

  end subroutine offlemb_perform

  !-----------------------------------------------------------------------------
#endif


end module caaba_exp

