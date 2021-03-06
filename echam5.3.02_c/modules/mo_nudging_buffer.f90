MODULE mo_nudging_buffer
!BOP

  ! !MODULE: mo_nudging_buffer (layer 4)

  ! !DESCRIPTION: 
  !    define work space arrays for nudging

  ! !REVISION HISTORY: 
  ! Ingo Kirchner, MPI Hamburg, April-2001
  ! Ingo Kirchner, MPI Hamburg, Aug-2002, revision

  ! !LAST CHANGE: 
  ! Joachim Buchholz, MPI Mainz, 29. July 2003

  ! !USES:
  USE mo_kind,        ONLY: dp, i4
  USE mo_linked_list, ONLY: t_stream
  USE mo_decomposition, ONLY: ldc => local_decomposition
!BOX
  IMPLICIT NONE
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: NudgingStreamInit
  PUBLIC :: NudgingAllocMem
  PUBLIC :: NudgingDeallocMem
  PUBLIC :: NudgingAllocRefMem
  PUBLIC :: NudgingDeallocRefMem

  PUBLIC :: NdgInitCounter
  PUBLIC :: NdgSetCounter
  PUBLIC :: NdgRemoveCounter
  PUBLIC :: NdgCorrBuffer
  PUBLIC :: NdgCleanBuffer
!BOX
  ! INPUT: reference data

  REAL(kind=dp), POINTER, DIMENSION(:,:,:), SAVE :: &
       stobs0, sdobs0, svobs0,              &! time level 1
       stobs1, sdobs1, svobs1,              &! time level 2
       stobs2, sdobs2, svobs2,              &! time level 3
       stobs3, sdobs3, svobs3,              &! time level 4
       stobs,  sdobs,  svobs                 ! time interpolated reference data
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, SAVE :: &
       buf_ref0, buf_ref1, buf_ref2, buf_ref3, buf_ref

  ! INPUT: fields with constant part of pattern correlation (LEVEL,2)

  REAL(kind=dp), ALLOCATABLE, DIMENSION(:), SAVE :: &
       a_pat_lnp,   b_pat_lnp,   a_nrm_lnp
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:), SAVE :: &
       a_pat_tem, b_pat_tem, a_nrm_tem, &
       a_pat_div, b_pat_div, a_nrm_div, &
       a_pat_vor, b_pat_vor, a_nrm_vor

!EOX
  ! STREAM definition
  TYPE(t_stream), POINTER :: nudg
#ifdef MESSY
  ! mz_ht_20050421+
  PUBLIC :: NUDG_SP2GP

  TYPE(t_stream), POINTER :: nudg_gp
  TYPE nudg_set
    REAL(dp), DIMENSION(:,:,:), POINTER :: nudge_gpfields => NULL()
  END TYPE nudg_set
  TYPE(nudg_set), DIMENSION(:), POINTER :: gp_nudge     => NULL()
  TYPE nudg_set_ori
    REAL(dp), DIMENSION(:,:,:), POINTER :: nudge_spfields => NULL()
  END TYPE nudg_set_ori
  TYPE(nudg_set_ori), DIMENSION(:), POINTER :: sp_nudge     => NULL()
  
  INTEGER, PARAMETER :: pot_outp_num = 28
  INTEGER, PARAMETER, DIMENSION(pot_outp_num) :: pot_outp_idx =&
                  (/1,2,3,4,5,6,7,8,9,10,11,12,13,&
                    14,15,16,17,18,19,20,27,28,29,&
                    30,43,44,45,46/)
  INTEGER :: output_num
  INTEGER, ALLOCATABLE :: outp_idx(:)
  ! mz_ht_20050421-
#endif
!BOX
  INTEGER                 :: nio_index    ! nudging stream time event index

  ! OUTPUT: nudging tendencies instantaneous (??ten) and mean (??tac)

  REAL(kind=dp), POINTER, DIMENSION(:,:,:), SAVE :: &
       ssten, stten, sdten, svten, sstac, sttac, sdtac, svtac
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET :: buf_sfc

  ! SSTEN, SSTAC means surface fields
  INTEGER :: nssfc = 1  ! number of surface spectral fields
                        ! 1  log surface pressure

  ! OUTPUT (optional): bias outside the nudged region
  ! 

  REAL(kind=dp), POINTER, DIMENSION(:,:,:), SAVE :: &
       isano, itano, idano, ivano, asano, atano, adano, avano
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET :: buf_ano

  INTEGER(i4), SAVE :: indg_accu = 0   !ik temporary fix

  ! OUTPUT (optional): fast mode fraction of initial tendency difference

  REAL(kind=dp), POINTER, DIMENSION(:,:,:) :: &
       stfast_a   , sdfast_a   , svfast_a   ,   &! time level a
       stfast_b   , sdfast_b   , svfast_b   ,   &! time level b
       stfast_accu, sdfast_accu, svfast_accu     ! accumulated
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET :: buf_fm

  INTEGER(i4), SAVE :: ifast_accu = 0 ! control accumulation of fast mode part
  LOGICAL, SAVE :: &
       lfill_a = .FALSE.,    & ! data for level a available?
       lfill_b = .FALSE.       ! data for level b available?

  ! OUTPUT (optional): systematic initial tendency error (SITE) calculation

  ! accumulated Systematic Initial Tendency Errors
  REAL(kind=dp), POINTER, DIMENSION(:,:,:), SAVE :: stsite, sdsite, svsite
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET :: buf_site

  ! INTERNAL: SITE buffer

  REAL(kind=dp), POINTER, DIMENSION(:,:,:), SAVE :: &
       st_o_n0, sd_o_n0, sv_o_n0,  &! time level 1 (reference)
       st_o_n1, sd_o_n1, sv_o_n1,  &! time level 2 (reference)
       st_m_n0, sd_m_n0, sv_m_n0,  &! time level 1 (model)
       st_m_n1, sd_m_n1, sv_m_n1    ! time level 2 (model)

  ! INTERNAL: nudging weights for observed data [1/sec], each for one level

  REAL(kind=dp)               :: nudgpo = 0.0
  REAL(kind=dp), ALLOCATABLE, SAVE  :: nudgto(:), nudgdo(:), nudgvo(:)

  ! INTERNAL: local fields with nudging coefficients A and B
  !           size = (NoSpectralCoeff * NoLevel)

  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: &
       nudgsa, nudgsb,     &! log pressure and other surface fields
       nudgta, nudgtb,     &! temperature
       nudgda, nudgdb,     &! divergence
       nudgva, nudgvb       ! vorticity

  ! INTERNAL: flag field,  set to 1 inside the nudging region

  REAL(kind=dp), ALLOCATABLE, TARGET, SAVE :: flagn(:,:,:)

  ! local definitions
  INTEGER(i4), PUBLIC, SAVE :: isiteaccu = 0 ! control accumulation of SITEs
  LOGICAL, PUBLIC, SAVE :: lsite_n0 = .FALSE., lsite_n1 = .FALSE.

  CHARACTER(len=256) :: mess

!EOX
  TYPE codes
        ! local definition of code table elements
    INTEGER            :: gc   ! GRIB code number
    CHARACTER(len=10)  :: sn   ! short name
    CHARACTER(len=100) :: ln   ! long name
    CHARACTER(len=30)  :: ud   ! unit definition
  END TYPE codes

! mz_hr_20150505: correction of units: in part1 all multiplied with 's'

  ! nudging code table, part 1: 
  TYPE(codes), PARAMETER, DIMENSION(20) :: nct_part1 = (/ &
     !  1 instantaneous residue: 
     codes(111,'NIPSFC','NDG log Psfc (inst)',         'ln(Pa)'), &
     codes(112,'NITEMP','NDG temperature (inst)',      'K'     ), &
     codes(113,'NIDIV', 'NDG divergence (inst)',       '1/s'   ), &
     codes(114,'NIVOR', 'NDG vorticity (inst)',        '1/s'   ), &
     !  5 mean residue: 
     codes(115,'NAPSFC','NDG log Psfc (mean)',         'ln(Pa)'), &
     codes(116,'NATEMP','NDG temperature (mean)',      'K'     ), &
     codes(117,'NADIV', 'NDG divergence (mean)',       '1/s'   ), &
     codes(118,'NAVOR', 'NDG vorticity (mean)',        '1/s'   ), &
     !  9 reference data: 
     codes( 31,'OPSFC', 'OBS log Psfc (inst)',         'ln(Pa)'), &
     codes( 32,'OTEMP', 'OBS temperature (inst)',      'K'     ), &
     codes( 33,'ODIV',  'OBS divergence (inst)',       '1/s'   ), &
     codes( 34,'OVOR',  'OBS vorticity (inst)',        '1/s'   ), &
     ! 13 instantaneous  bias: 
     codes( 21,'DIPSFC','MOD-OBS log Psfc (inst)',     'ln(Pa)'), &
     codes( 22,'DITEMP','MOD-OBS temperature (inst)',  'K'     ), &
     codes( 23,'DIDIV', 'MOD-OBS divergence (inst)',   '1/s'   ), &
     codes( 24,'DIVOR', 'MOD-OBS vorticity (inst)',    '1/s'   ), &
     ! 17  mean bias: 
     codes( 25,'DAPSFC','MOD-OBS log Psfc (mean)',     'ln(Pa)'), &
     codes( 26,'DATEMP','MOD-OBS temperature (mean)',  'K'     ), &
     codes( 27,'DADIV', 'MOD-OBS divergence (mean)',   '1/s'   ), &
     codes( 28,'DAVOR', 'MOD-OBS vorticity (mean)',    '1/s'   ) /)

  ! nudging code table, part 2: 
  TYPE(codes), PARAMETER, DIMENSION(26) :: nct_part2 = (/ &
     ! 21  FM workspace: 
     codes( -1,'FMTPA', 'FM STP buffer A',             'X/s'     ), &
     codes( -1,'FMDA',  'FM DIV buffer A',             '1/s^2'   ), &
     codes( -1,'FMVA',  'FM VOR buffer A',             '1/s^2'   ), &
     codes( -1,'FMTPB', 'FM STP buffer B',             'X/s'     ), &
     codes( -1,'FMDB',  'FM DIV buffer B',             '1/s^2'   ), &
     codes( -1,'FMVB',  'FM VOR buffer B',             '1/s^2'   ), &
     ! 27 FM diagnostics: 
     codes( 35,'FMPSFC','FM-Filter log Psfc (mean)',   'ln(Pa)/s'), &
     codes( 36,'FMTEMP','FM-Filter temperature (mean)','K/s'     ), &
     codes( 37,'FMDIV', 'FM-Filter divergence (mean)', '1/s^2'   ), &
     codes( 38,'FMVOR', 'FM-Filter vorticity (mean)',  '1/s^2'   ), &
     ! 31 SITE worspace: 
     codes( -1,'MTP0',  'STP M-buffer t0',             'X/s'     ), &
     codes( -1,'MD0',   'DIV M-buffer t0',             '1/s^2'   ), &
     codes( -1,'MV0',   'VOR M-buffer t0',             '1/s^2'   ), &
     codes( -1,'MTP1',  'STP M-buffer t1',             'X/s'     ), &
     codes( -1,'MD1',   'DIV M-buffer t1',             '1/s^2'   ), &
     codes( -1,'MV1',   'VOR M-buffer t1',             '1/s^2'   ), &
     codes( -1,'OTP0',  'STP R-buffer t0',             'X/s'     ), &
     codes( -1,'OD0',   'DIV R-buffer t0',             '1/s^2'   ), &
     codes( -1,'OV0',   'VOR R-buffer t0',             '1/s^2'   ), &
     codes( -1,'OTP1',  'STP R-buffer t1',             'X/s'     ), &
     codes( -1,'OD1',   'DIV R-buffer t1',             '1/s^2'   ), &
     codes( -1,'OV1',   'VOR R-buffer t1',             '1/s^2'   ), &
     ! 43 SITE diagnostics: 
     codes(121,'SIPSFC','SITE log Psfc (mean)',        'ln(Pa)/s'), &
     codes(122,'SITEMP','SITE temperature (mean)',     'K/s'     ), &
     codes(123,'SIDIV', 'SITE divergence (mean)',      '1/s^2'   ), &
     codes(124,'SIVOR', 'SITE vorticity (mean)',       '1/s^2'   ) /)

  ! complete nudging code table: 
  TYPE(codes), PARAMETER, DIMENSION(46) :: nct = (/nct_part1, nct_part2/)

!BOX
!EOX
!EOP
CONTAINS
  !=============================================================================
!BOP
  !
  ! !IROUTINE:  NudgingStreamInit
  ! !INTERFACE:

  SUBROUTINE NudgingStreamInit
    ! !DESCRIPTION:
    ! initialize nudging output stream

    ! !USES:
    USE mo_exception,         ONLY: message, finish
    USE mo_nudging_constants, ONLY: lnudgwobs, lsite
#ifndef MESSY
    USE mo_memory_base,       ONLY: new_stream, GRIB, &
                                    SPECTRAL, default_stream_setting, &
                                    add_stream_element, SURFACE
#else
    ! mz_ht_20050421+
    USE mo_memory_base,       ONLY: new_stream, GRIB, &
                                    SPECTRAL, default_stream_setting, &
                                    add_stream_element, SURFACE,      &
                                    GAUSSIAN, NETCDF      
    ! mz_ht_20050421-
#endif
    USE mo_control,           ONLY: nsp, nlev, lnmi, nlevp1
    USE mo_grib,              ONLY: nudging_table
!EOP
!BOC
!BOX

    INTEGER :: snsp, ierr

    REAL(kind=dp), POINTER :: p4(:,:,:,:), p3(:,:,:)
#ifdef MESSY
    ! mz_ht_20050421+   
    INTEGER :: level,i,j  
    ALLOCATE(sp_nudge(46))
    ! mz_ht_20050421-
#endif
    snsp = ldc%snsp

!!!    nct(1) = 111 , 'NIPSFC', 'NDG log Psfc (inst)','ln(Pa)/s'

    CALL new_stream( nudg, 'nudg', &         ! define nudging stream defaults
         lpost=.TRUE.,  lpout=.FALSE., lrerun=.TRUE., &
         filetype=GRIB )

    nio_index = nudg%post_idx

    CALL default_stream_setting ( nudg, &
         bits = 24, &
         lpost = .TRUE., lrerun = .FALSE., &
         laccu = .FALSE., contnorest = .TRUE., &
         ldims=(/nlev, 2, snsp/), &
         gdims=(/nlev, 2,  nsp/), &
         repr = SPECTRAL, table=nudging_table  )

    !--- define output stream

    ALLOCATE (buf_sfc(nssfc,2,snsp,2))
    ssten => buf_sfc(:,:,:,1)
    sstac => buf_sfc(:,:,:,2)

#ifdef MESSY
    ! mz_ht_20050421+    

    CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .FALSE.)
    p4 => buf_sfc(1:1,:,:,1:1)
    CALL add_stream_element(nudg, nct(1)%sn, p3, p4=p4, code=nct(1)%gc, &
         longname=nct(1)%ln, units=nct(1)%ud,                           &
         ldims=(/1, 2, snsp/), gdims=(/1, 2, nsp/),                     &
         klev=1, leveltype=SURFACE)
    sp_nudge(1)%nudge_spfields =>  buf_sfc(1:1,:,:,1)! mz_ht_20050421
    CALL add_stream_element(nudg, nct(2)%sn, stten,  code=nct(2)%gc,    &
         longname=nct(2)%ln, units=nct(2)%ud)
    sp_nudge(2)%nudge_spfields =>  stten! mz_ht_20050421
    CALL add_stream_element(nudg, nct(3)%sn, sdten,  code=nct(3)%gc,    &
         longname=nct(3)%ln, units=nct(3)%ud)
    sp_nudge(3)%nudge_spfields =>  sdten! mz_ht_20050421
    CALL add_stream_element(nudg, nct(4)%sn, svten,  code=nct(4)%gc,    &
         longname=nct(4)%ln, units=nct(4)%ud)
    sp_nudge(4)%nudge_spfields =>  svten! mz_ht_20050421

    CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .TRUE.)
    p4 => buf_sfc(1:1,:,:,2:2)
    CALL add_stream_element(nudg, nct(5)%sn, p3, p4=p4, code=nct(5)%gc, &
         longname=nct(5)%ln, units=nct(5)%ud,                           &
!         ldims=(/1, 2, snsp/) ,gdims=(/1, 2, nsp/),                     &
         klev=1, leveltype=SURFACE)
    sp_nudge(5)%nudge_spfields =>  buf_sfc(1:1,:,:,2)! mz_ht_20050421
    CALL add_stream_element(nudg, nct(6)%sn, sttac,  code=nct(6)%gc,    &
         longname=nct(6)%ln, units=nct(6)%ud)
    sp_nudge(6)%nudge_spfields =>  sttac! mz_ht_20050421
    CALL add_stream_element(nudg, nct(7)%sn, sdtac,  code=nct(7)%gc,    &
         longname=nct(7)%ln, units=nct(7)%ud)
    sp_nudge(7)%nudge_spfields =>  sdtac! mz_ht_20050421
    CALL add_stream_element(nudg, nct(8)%sn, svtac,  code=nct(8)%gc,    &
         longname=nct(8)%ln, units=nct(8)%ud)
    sp_nudge(8)%nudge_spfields =>  svtac! mz_ht_20050421

    sstac(:,:,:) = 0.0_dp
    sttac(:,:,:) = 0.0_dp
    sdtac(:,:,:) = 0.0_dp
    svtac(:,:,:) = 0.0_dp

    CALL message('',' Memory for TEN, TAC and SITE allocated.')

    CALL default_stream_setting ( nudg, &
         contnorest = .TRUE., &
         lpost = .TRUE., lrerun = .FALSE.)

    ALLOCATE(buf_ref(nlevp1,2,snsp,3), stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingInit', &
       'Allocation error observed fields failure.')
    stobs => buf_ref(:    ,:,:,1)
    sdobs => buf_ref(:nlev,:,:,2)
    svobs => buf_ref(:nlev,:,:,3)
    buf_ref (:,:,:,:) = 0.0_dp

    p4 => buf_ref(nlevp1:nlevp1,:,:,1:1)
    CALL add_stream_element(nudg, nct( 9)%sn, p3, p4=p4, code=nct( 9)%gc, &
         longname=nct( 9)%ln, units=nct( 9)%ud,                           &
         klev=1, leveltype=SURFACE)
    sp_nudge(9)%nudge_spfields =>  buf_ref(nlevp1:nlevp1,:,:,1)! mz_ht_20050421
    p4 => buf_ref(      :nlev,  :,:,1:1)
    CALL add_stream_element(nudg, nct(10)%sn, p3, p4=p4, code=nct(10)%gc, &
         longname=nct(10)%ln, units=nct(10)%ud)
    sp_nudge(10)%nudge_spfields =>  buf_ref(1:nlev,:,:,1)! mz_ht_20050421
    p4 => buf_ref(      :nlev,  :,:,2:2)
    CALL add_stream_element(nudg, nct(11)%sn, p3, p4=p4, code=nct(11)%gc, &
         longname=nct(11)%ln, units=nct(11)%ud)
    sp_nudge(11)%nudge_spfields =>  buf_ref(1:nlev,:,:,2)! mz_ht_20050421
    p4 => buf_ref(      :nlev,  :,:,3:3)
    CALL add_stream_element(nudg, nct(12)%sn, p3, p4=p4, code=nct(12)%gc, &
         longname=nct(12)%ln, units=nct(12)%ud)
    sp_nudge(12)%nudge_spfields =>  buf_ref(1:nlev,:,:,3)! mz_ht_20050421

    CALL message('',' Memory for reference data allocated.')

    ! diagnostics output
    IF (lnudgwobs) THEN
      ALLOCATE(buf_ano(nssfc,2,snsp,2))
      isano => buf_ano(:,:,:,1)
      asano => buf_ano(:,:,:,2)

      CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .FALSE.)
      p4 => buf_ano(1:1,:,:,1:1)
      CALL add_stream_element(nudg, nct(13)%sn, p3, p4=p4, code=nct(13)%gc, &
           longname=nct(13)%ln, units=nct(13)%ud,                           &
           klev=1,leveltype=SURFACE)
      sp_nudge(13)%nudge_spfields =>  buf_ano(1:1,:,:,1)! mz_ht_20050421
      CALL add_stream_element(nudg, nct(14)%sn, itano, code=nct(14)%gc,     &
           longname=nct(14)%ln, units=nct(14)%ud)
      sp_nudge(14)%nudge_spfields =>  itano! mz_ht_20050421
      CALL add_stream_element(nudg, nct(15)%sn, idano, code=nct(15)%gc,     &
           longname=nct(15)%ln, units=nct(15)%ud)
      sp_nudge(15)%nudge_spfields =>  idano! mz_ht_20050421
      CALL add_stream_element(nudg, nct(16)%sn, ivano, code=nct(16)%gc,     &
           longname=nct(16)%ln, units=nct(16)%ud)
      sp_nudge(16)%nudge_spfields =>  ivano! mz_ht_20050421

      CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .TRUE.)
      ! mz_ht_20050421+
!      p4 => buf_ano(1:1,:,:,1:1)       ! should be wrong average value wanted ?!?
      p4 => buf_ano(1:1,:,:,2:2)      
      ! mz_ht_20050421-
      CALL add_stream_element(nudg, nct(17)%sn, p3, p4=p4, code=nct(17)%gc, &
           longname=nct(17)%ln,units=nct(17)%ud, &
           klev=1,leveltype=SURFACE)
      sp_nudge(17)%nudge_spfields =>  buf_ano(1:1,:,:,2)! mz_ht_20050421
      CALL add_stream_element(nudg, nct(18)%sn, atano, code=nct(18)%gc,     &
           longname=nct(18)%ln, units=nct(18)%ud)
      sp_nudge(18)%nudge_spfields =>  atano! mz_ht_20050421
      CALL add_stream_element(nudg, nct(19)%sn, adano, code=nct(19)%gc,     &
           longname=nct(19)%ln, units=nct(19)%ud)
      sp_nudge(19)%nudge_spfields =>  adano! mz_ht_20050421
      CALL add_stream_element(nudg, nct(20)%sn, avano, code=nct(20)%gc,     &
           longname=nct(20)%ln, units=nct(20)%ud)
      sp_nudge(20)%nudge_spfields =>  avano! mz_ht_20050421

      asano(:,:,:) = 0.0_dp
      atano(:,:,:) = 0.0_dp
      adano(:,:,:) = 0.0_dp
      avano(:,:,:) = 0.0_dp

      CALL message('',' Memory for MOD-REF diagnostics allocated.')
    END IF

    !--- allocate memory for fast modes
    !    needed for tendency diagnostics using SNMI
    IF (lnmi) THEN
      CALL default_stream_setting ( nudg, lpost = .FALSE., lrerun = .TRUE.)
      CALL add_stream_element(nudg, nct(21)%sn, stfast_a,        &
           longname=nct(21)%ln,units=nct(21)%ud,                 &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2,  nsp/), &
           klev=nlevp1)
      sp_nudge(21)%nudge_spfields =>  stfast_a! mz_ht_20050421
      CALL add_stream_element(nudg, nct(22)%sn, sdfast_a,        &
           longname=nct(22)%ln,units=nct(22)%ud)
      sp_nudge(22)%nudge_spfields =>  sdfast_a! mz_ht_20050421  
      CALL add_stream_element(nudg, nct(23)%sn, svfast_a,        &
           longname=nct(23)%ln,units=nct(23)%ud)
      sp_nudge(23)%nudge_spfields =>  svfast_a! mz_ht_20050421  
      CALL add_stream_element(nudg, nct(24)%sn, stfast_b,        &
           longname=nct(24)%ln,units=nct(24)%ud,                 &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2,  nsp/), &
           klev=nlevp1)
      sp_nudge(24)%nudge_spfields =>  stfast_b! mz_ht_20050421 
      CALL add_stream_element(nudg, nct(25)%sn, sdfast_b,        &
           longname=nct(25)%ln,units=nct(25)%ud)
      sp_nudge(25)%nudge_spfields =>  sdfast_b! mz_ht_20050421
      CALL add_stream_element(nudg, nct(26)%sn, svfast_b,        &
           longname=nct(26)%ln,units=nct(26)%ud)
      sp_nudge(26)%nudge_spfields =>  svfast_b! mz_ht_20050421

      ALLOCATE(buf_fm(nlevp1,2,snsp,3))
      stfast_accu => buf_fm(:,    :,:,1)
      sdfast_accu => buf_fm(:nlev,:,:,2)
      svfast_accu => buf_fm(:nlev,:,:,3)
      buf_fm(:,:,:,:) = 0.0_dp

      CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .TRUE.)
      p4 => buf_fm(nlevp1:nlevp1,:,:,1:1)
      CALL add_stream_element(nudg, nct(27)%sn, p3, p4=p4, code=nct(27)%gc, &
           longname=nct(27)%ln, units=nct(27)%ud,                           &
           klev=1, leveltype=SURFACE)
      sp_nudge(27)%nudge_spfields => buf_fm(nlevp1:nlevp1,:,:,1)! mz_ht_20050421
      p4 => buf_fm(      :nlev,  :,:,1:1)
      CALL add_stream_element(nudg, nct(28)%sn, p3, p4=p4, code=nct(28)%gc, &
           longname=nct(28)%ln, units=nct(28)%ud)
      sp_nudge(28)%nudge_spfields =>  buf_fm(1:nlev,  :,:,1)! mz_ht_20050421
      p4 => buf_fm(      :nlev,  :,:,2:2)
      CALL add_stream_element(nudg, nct(29)%sn, p3, p4=p4, code=nct(29)%gc, &
           longname=nct(29)%ln, units=nct(29)%ud)
      sp_nudge(29)%nudge_spfields =>  buf_fm(1:nlev,  :,:,2)! mz_ht_20050421
      p4 => buf_fm(      :nlev,  :,:,3:3)
      CALL add_stream_element(nudg, nct(30)%sn, p3, p4=p4, code=nct(30)%gc, &
           longname=nct(30)%ln, units=nct(30)%ud)
      sp_nudge(30)%nudge_spfields =>  buf_fm(1:nlev,  :,:,3)! mz_ht_20050421

      CALL message('',' Memory for NMI diagnostics allocated.')
    END IF

    !--- allocate the memory for additional diagnostics
    IF (lsite) THEN
      CALL default_stream_setting ( nudg, lpost = .FALSE., lrerun = .TRUE.)
      CALL add_stream_element(nudg, nct(31)%sn, st_m_n0,         &
           longname=nct(31)%ln, units=nct(31)%ud,                &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2,  nsp/), &
           klev=nlevp1)
      sp_nudge(31)%nudge_spfields =>  st_m_n0! mz_ht_20050421
      CALL add_stream_element(nudg, nct(32)%sn, sd_m_n0,         &
           longname=nct(32)%ln, units=nct(32)%ud)
      sp_nudge(32)%nudge_spfields =>  sd_m_n0! mz_ht_20050421
      CALL add_stream_element(nudg, nct(33)%sn, sv_m_n0,         &
           longname=nct(33)%ln, units=nct(33)%ud)
      sp_nudge(33)%nudge_spfields =>  sv_m_n0! mz_ht_20050421
      CALL add_stream_element(nudg, nct(34)%sn, st_m_n1,         &
           longname=nct(34)%ln, units=nct(34)%ud,                &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2, nsp/),  &
           klev=nlevp1)
      sp_nudge(34)%nudge_spfields =>  st_m_n1! mz_ht_20050421
      CALL add_stream_element(nudg, nct(35)%sn, sd_m_n1,         &
           longname=nct(35)%ln, units=nct(35)%ud)
      sp_nudge(35)%nudge_spfields =>  sd_m_n1! mz_ht_20050421
      CALL add_stream_element(nudg, nct(36)%sn, sv_m_n1,         &
           longname=nct(36)%ln, units=nct(36)%ud)
      sp_nudge(36)%nudge_spfields =>  sv_m_n1! mz_ht_20050421
      CALL add_stream_element(nudg, nct(37)%sn, st_o_n0,         &
           longname=nct(37)%ln, units=nct(37)%ud, &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2,  nsp/), &
           klev=nlevp1)
      sp_nudge(37)%nudge_spfields =>  st_o_n0! mz_ht_20050421
      CALL add_stream_element(nudg, nct(38)%sn, sd_o_n0,         &
           longname=nct(38)%ln, units=nct(38)%ud)
      sp_nudge(38)%nudge_spfields =>  sd_o_n0! mz_ht_20050421
      CALL add_stream_element(nudg, nct(39)%sn, sv_o_n0,         &
           longname=nct(39)%ln, units=nct(39)%ud)
      sp_nudge(39)%nudge_spfields =>  sv_o_n0! mz_ht_20050421
      CALL add_stream_element(nudg, nct(40)%sn, st_o_n1,         &
           longname=nct(40)%ln, units=nct(40)%ud, &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2, nsp/),  &
           klev=nlevp1)
      sp_nudge(40)%nudge_spfields =>  st_o_n1! mz_ht_20050421
      CALL add_stream_element(nudg, nct(41)%sn, sd_o_n1,         &
           longname=nct(41)%ln, units=nct(41)%ud)
      sp_nudge(41)%nudge_spfields =>  sd_o_n1! mz_ht_20050421
      CALL add_stream_element(nudg, nct(42)%sn, sv_o_n1,         &
           longname=nct(42)%ln, units=nct(42)%ud)
      sp_nudge(42)%nudge_spfields =>  sv_o_n1! mz_ht_20050421

      ALLOCATE(buf_site(nlevp1,2,snsp,3))
      stsite => buf_site(:,    :,:,1)
      sdsite => buf_site(:nlev,:,:,2)
      svsite => buf_site(:nlev,:,:,3)
      buf_site(:,:,:,:) = 0.0_dp

      CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .TRUE.)
      p4 => buf_site(nlevp1:nlevp1,:,:,1:1)
      CALL add_stream_element(nudg, nct(43)%sn, p3, p4=p4, code=nct(43)%gc, &
           longname=nct(43)%ln, units=nct(43)%ud,                           &
           klev=1, leveltype=SURFACE)
      sp_nudge(43)%nudge_spfields =>  buf_site(nlevp1:nlevp1,  :,:,1)! mz_ht_20050421
      p4 => buf_site(      :nlev,  :,:,1:1)
      CALL add_stream_element(nudg, nct(44)%sn, p3, p4=p4, code=nct(44)%gc, &
           longname=nct(44)%ln, units=nct(44)%ud)
      sp_nudge(44)%nudge_spfields =>  buf_site(1:nlev,  :,:,1)! mz_ht_20050421
      p4 => buf_site(      :nlev,  :,:,2:2)
      CALL add_stream_element(nudg, nct(45)%sn, p3, p4=p4, code=nct(45)%gc, &
           longname=nct(45)%ln, units=nct(45)%ud)
      sp_nudge(45)%nudge_spfields =>  buf_site(1:nlev,  :,:,2)! mz_ht_20050421
      p4 => buf_site(      :nlev,  :,:,3:3)
      CALL add_stream_element(nudg, nct(46)%sn, p3, p4=p4, code=nct(46)%gc, &
           longname=nct(46)%ln, units=nct(46)%ud)
      sp_nudge(46)%nudge_spfields =>  buf_site(1:nlev,  :,:,3)! mz_ht_20050421

      CALL message('',' Memory for SITE diagnostics allocated.')
    END IF

    CALL MESSAGE('mo_nudging_init', 'Building Gridpoint Nudging diagnostic stream (start)')

    CALL new_stream( nudg_gp, 'nudg_gp',             &    ! define nudging_gp stream defaults
         lpost=.TRUE., lrerun=.FALSE.,               &
         filetype=NETCDF )
    CALL default_stream_setting ( nudg_gp, lpost=.TRUE., lrerun=.FALSE., &
      LACCU = .FALSE., CONTNOREST = .TRUE., repr= GAUSSIAN)

    
    j=0
    DO i=1, pot_outp_num
      if (ASSOCIATED(sp_nudge(pot_outp_idx(i))%nudge_spfields)) j=j+1
    END DO

    output_num=j
    ALLOCATE(outp_idx(output_num))
    ALLOCATE(gp_nudge(output_num))
    
    j=0
    DO i=1, pot_outp_num
      if (ASSOCIATED(sp_nudge(pot_outp_idx(i))%nudge_spfields)) then
        j=j+1
        outp_idx(j) = pot_outp_idx(i)
        SELECT CASE(outp_idx(j))        
        CASE(1,5,9,13,17,27,43)
          level = 1
        CASE DEFAULT
          level = nlev
        END SELECT
        CALL add_stream_element(nudg_gp, nct(pot_outp_idx(i))%sn,           &
             gp_nudge(j)%nudge_gpfields, longname=nct(pot_outp_idx(i))%ln,  &
             units=nct(pot_outp_idx(i))%ud, klev=level) 
      ENDIF
    END DO
      
    CALL MESSAGE('mo_nudging_init', 'Building Gridpoint Nudging diagnostic stream (finished)')
    ! mz_ht_20050421-
#else
    
    CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .FALSE.)
    p4 => buf_sfc(1:1,:,:,1:1)
    CALL add_stream_element(nudg, nct(1)%sn, p3, p4=p4, code=nct(1)%gc, &
         longname=nct(1)%ln, units=nct(1)%ud,                           &
         ldims=(/1, 2, snsp/), gdims=(/1, 2, nsp/),                     &
         klev=1, leveltype=SURFACE)
    CALL add_stream_element(nudg, nct(2)%sn, stten,  code=nct(2)%gc,    &
         longname=nct(2)%ln, units=nct(2)%ud)
    CALL add_stream_element(nudg, nct(3)%sn, sdten,  code=nct(3)%gc,    &
         longname=nct(3)%ln, units=nct(3)%ud)
    CALL add_stream_element(nudg, nct(4)%sn, svten,  code=nct(4)%gc,    &
         longname=nct(4)%ln, units=nct(4)%ud)

    CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .TRUE.)
    p4 => buf_sfc(1:1,:,:,2:2)
    CALL add_stream_element(nudg, nct(5)%sn, p3, p4=p4, code=nct(5)%gc, &
         longname=nct(5)%ln, units=nct(5)%ud,                           &
         ldims=(/1, 2, snsp/) ,gdims=(/1, 2, nsp/),                     &
         klev=1, leveltype=SURFACE)
    CALL add_stream_element(nudg, nct(6)%sn, sttac,  code=nct(6)%gc,    &
         longname=nct(6)%ln, units=nct(6)%ud)
    CALL add_stream_element(nudg, nct(7)%sn, sdtac,  code=nct(7)%gc,    &
         longname=nct(7)%ln, units=nct(7)%ud)
    CALL add_stream_element(nudg, nct(8)%sn, svtac,  code=nct(8)%gc,    &
         longname=nct(8)%ln, units=nct(8)%ud)

    sstac(:,:,:) = 0.0_dp
    sttac(:,:,:) = 0.0_dp
    sdtac(:,:,:) = 0.0_dp
    svtac(:,:,:) = 0.0_dp

    CALL message('',' Memory for TEN, TAC and SITE allocated.')

    CALL default_stream_setting ( nudg, &
         contnorest = .TRUE., &
         lpost = .TRUE., lrerun = .FALSE.)

    ALLOCATE(buf_ref(nlevp1,2,snsp,3), stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingInit', &
       'Allocation error observed fields failure.')
    stobs => buf_ref(:    ,:,:,1)
    sdobs => buf_ref(:nlev,:,:,2)
    svobs => buf_ref(:nlev,:,:,3)
    buf_ref (:,:,:,:) = 0.0_dp

    p4 => buf_ref(nlevp1:nlevp1,:,:,1:1)
    CALL add_stream_element(nudg, nct( 9)%sn, p3, p4=p4, code=nct( 9)%gc, &
         longname=nct( 9)%ln, units=nct( 9)%ud,                           &
         klev=1, leveltype=SURFACE)
    p4 => buf_ref(      :nlev,  :,:,1:1)
    CALL add_stream_element(nudg, nct(10)%sn, p3, p4=p4, code=nct(10)%gc, &
         longname=nct(10)%ln, units=nct(10)%ud)
    p4 => buf_ref(      :nlev,  :,:,2:2)
    CALL add_stream_element(nudg, nct(11)%sn, p3, p4=p4, code=nct(11)%gc, &
         longname=nct(11)%ln, units=nct(11)%ud)
    p4 => buf_ref(      :nlev,  :,:,3:3)
    CALL add_stream_element(nudg, nct(12)%sn, p3, p4=p4, code=nct(12)%gc, &
         longname=nct(12)%ln, units=nct(12)%ud)

    CALL message('',' Memory for reference data allocated.')

    ! diagnostics output
    IF (lnudgwobs) THEN
      ALLOCATE(buf_ano(nssfc,2,snsp,2))
      isano => buf_ano(:,:,:,1)
      asano => buf_ano(:,:,:,2)

      CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .FALSE.)
      p4 => buf_ano(1:1,:,:,1:1)
      CALL add_stream_element(nudg, nct(13)%sn, p3, p4=p4, code=nct(13)%gc, &
           longname=nct(13)%ln, units=nct(13)%ud,                           &
           klev=1,leveltype=SURFACE)
      CALL add_stream_element(nudg, nct(14)%sn, itano, code=nct(14)%gc,     &
           longname=nct(14)%ln, units=nct(14)%ud)
      CALL add_stream_element(nudg, nct(15)%sn, idano, code=nct(15)%gc,     &
           longname=nct(15)%ln, units=nct(15)%ud)
      CALL add_stream_element(nudg, nct(16)%sn, ivano, code=nct(16)%gc,     &
           longname=nct(16)%ln, units=nct(16)%ud)

      CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .TRUE.)
      p4 => buf_ano(1:1,:,:,1:1)
      CALL add_stream_element(nudg, nct(17)%sn, p3, p4=p4, code=nct(17)%gc, &
           longname=nct(17)%ln,units=nct(17)%ud, &
           klev=1,leveltype=SURFACE)
      CALL add_stream_element(nudg, nct(18)%sn, atano, code=nct(18)%gc,     &
           longname=nct(18)%ln, units=nct(18)%ud)
      CALL add_stream_element(nudg, nct(19)%sn, adano, code=nct(19)%gc,     &
           longname=nct(19)%ln, units=nct(19)%ud)
      CALL add_stream_element(nudg, nct(20)%sn, avano, code=nct(20)%gc,     &
           longname=nct(20)%ln, units=nct(20)%ud)

      asano(:,:,:) = 0.0_dp
      atano(:,:,:) = 0.0_dp
      adano(:,:,:) = 0.0_dp
      avano(:,:,:) = 0.0_dp

      CALL message('',' Memory for MOD-REF diagnostics allocated.')
    END IF

    !--- allocate memory for fast modes
    !    needed for tendency diagnostics using SNMI
    IF (lnmi) THEN
      CALL default_stream_setting ( nudg, lpost = .FALSE., lrerun = .TRUE.)
      CALL add_stream_element(nudg, nct(21)%sn, stfast_a,        &
           longname=nct(21)%ln,units=nct(21)%ud,                 &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2,  nsp/), &
           klev=nlevp1)
      CALL add_stream_element(nudg, nct(22)%sn, sdfast_a,        &
           longname=nct(22)%ln,units=nct(22)%ud)
      CALL add_stream_element(nudg, nct(23)%sn, svfast_a,        &
           longname=nct(23)%ln,units=nct(23)%ud)
      CALL add_stream_element(nudg, nct(24)%sn, stfast_b,        &
           longname=nct(24)%ln,units=nct(24)%ud,                 &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2,  nsp/), &
           klev=nlevp1)
      CALL add_stream_element(nudg, nct(25)%sn, sdfast_b,        &
           longname=nct(25)%ln,units=nct(25)%ud)
      CALL add_stream_element(nudg, nct(26)%sn, svfast_b,        &
           longname=nct(26)%ln,units=nct(26)%ud)

      ALLOCATE(buf_fm(nlevp1,2,snsp,3))
      stfast_accu => buf_fm(:,    :,:,1)
      sdfast_accu => buf_fm(:nlev,:,:,2)
      svfast_accu => buf_fm(:nlev,:,:,3)
      buf_fm(:,:,:,:) = 0.0_dp

      CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .TRUE.)
      p4 => buf_fm(nlevp1:nlevp1,:,:,1:1)
      CALL add_stream_element(nudg, nct(27)%sn, p3, p4=p4, code=nct(27)%gc, &
           longname=nct(27)%ln, units=nct(27)%ud,                           &
           klev=1, leveltype=SURFACE)
      p4 => buf_fm(      :nlev,  :,:,1:1)
      CALL add_stream_element(nudg, nct(28)%sn, p3, p4=p4, code=nct(28)%gc, &
           longname=nct(28)%ln, units=nct(28)%ud)
      p4 => buf_fm(      :nlev,  :,:,2:2)
      CALL add_stream_element(nudg, nct(29)%sn, p3, p4=p4, code=nct(29)%gc, &
           longname=nct(29)%ln, units=nct(29)%ud)
      p4 => buf_fm(      :nlev,  :,:,3:3)
      CALL add_stream_element(nudg, nct(30)%sn, p3, p4=p4, code=nct(30)%gc, &
           longname=nct(30)%ln, units=nct(30)%ud)

      CALL message('',' Memory for NMI diagnostics allocated.')
    END IF

    !--- allocate the memory for additional diagnostics
    IF (lsite) THEN
      CALL default_stream_setting ( nudg, lpost = .FALSE., lrerun = .TRUE.)
      CALL add_stream_element(nudg, nct(31)%sn, st_m_n0,         &
           longname=nct(31)%ln, units=nct(31)%ud,                &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2,  nsp/), &
           klev=nlevp1)
      CALL add_stream_element(nudg, nct(32)%sn, sd_m_n0,         &
           longname=nct(32)%ln, units=nct(32)%ud)
      CALL add_stream_element(nudg, nct(33)%sn, sv_m_n0,         &
           longname=nct(33)%ln, units=nct(33)%ud)
      CALL add_stream_element(nudg, nct(34)%sn, st_m_n1,         &
           longname=nct(34)%ln, units=nct(34)%ud,                &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2, nsp/),  &
           klev=nlevp1)
      CALL add_stream_element(nudg, nct(35)%sn, sd_m_n1,         &
           longname=nct(35)%ln, units=nct(35)%ud)
      CALL add_stream_element(nudg, nct(36)%sn, sv_m_n1,         &
           longname=nct(36)%ln, units=nct(36)%ud)
      CALL add_stream_element(nudg, nct(37)%sn, st_o_n0,         &
           longname=nct(37)%ln, units=nct(37)%ud, &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2,  nsp/), &
           klev=nlevp1)
      CALL add_stream_element(nudg, nct(38)%sn, sd_o_n0,         &
           longname=nct(38)%ln, units=nct(38)%ud)
      CALL add_stream_element(nudg, nct(39)%sn, sv_o_n0,         &
           longname=nct(39)%ln, units=nct(39)%ud)
      CALL add_stream_element(nudg, nct(40)%sn, st_o_n1,         &
           longname=nct(40)%ln, units=nct(40)%ud, &
           ldims=(/nlevp1, 2, snsp/), gdims=(/nlevp1, 2, nsp/),  &
           klev=nlevp1)
      CALL add_stream_element(nudg, nct(41)%sn, sd_o_n1,         &
           longname=nct(41)%ln, units=nct(41)%ud)
      CALL add_stream_element(nudg, nct(42)%sn, sv_o_n1,         &
           longname=nct(42)%ln, units=nct(42)%ud)

      ALLOCATE(buf_site(nlevp1,2,snsp,3))
      stsite => buf_site(:,    :,:,1)
      sdsite => buf_site(:nlev,:,:,2)
      svsite => buf_site(:nlev,:,:,3)
      buf_site(:,:,:,:) = 0.0_dp

      CALL default_stream_setting ( nudg, lpost = .TRUE., lrerun = .TRUE.)
      p4 => buf_site(nlevp1:nlevp1,:,:,1:1)
      CALL add_stream_element(nudg, nct(43)%sn, p3, p4=p4, code=nct(43)%gc, &
           longname=nct(43)%ln, units=nct(43)%ud,                           &
           klev=1, leveltype=SURFACE)
      p4 => buf_site(      :nlev,  :,:,1:1)
      CALL add_stream_element(nudg, nct(44)%sn, p3, p4=p4, code=nct(44)%gc, &
           longname=nct(44)%ln, units=nct(44)%ud)
      p4 => buf_site(      :nlev,  :,:,2:2)
      CALL add_stream_element(nudg, nct(45)%sn, p3, p4=p4, code=nct(45)%gc, &
           longname=nct(45)%ln, units=nct(45)%ud)
      p4 => buf_site(      :nlev,  :,:,3:3)
      CALL add_stream_element(nudg, nct(46)%sn, p3, p4=p4, code=nct(46)%gc, &
           longname=nct(46)%ln, units=nct(46)%ud)

      CALL message('',' Memory for SITE diagnostics allocated.')
    END IF

#endif
  END SUBROUTINE NudgingStreamInit
!EOX
!EOC
  !=============================================================================
!BOP
! !IROUTINE:  NudgingAllocMem
! !INTERFACE:

  SUBROUTINE NudgingAllocMem

    ! !DESCRIPTION: 
    ! allocate memory for nudging, define nudging stream

    ! !USES:
    USE mo_exception,         ONLY: message, finish
    USE mo_nudging_constants, ONLY: lnudgpat
    USE mo_control,           ONLY: nlev, nlevp1
!EOP
!BOC
!BOX

    INTEGER :: snsp, ierr

    INTRINSIC allocated

    snsp = ldc%snsp

!EOX
    !--- set work space nudging coefficients
!BOX
    IF (.NOT. ALLOCATED(nudgda)) ALLOCATE(           &
         nudgta(nlev,2,snsp), nudgda(nlev,2,snsp),   &
         nudgva(nlev,2,snsp), nudgsa(nssfc,2,snsp) )
    IF (.NOT. ALLOCATED(nudgdb)) ALLOCATE(           &
         nudgtb(nlev,2,snsp), nudgdb(nlev,2,snsp),   &
         nudgvb(nlev,2,snsp), nudgsb(nssfc,2,snsp) )
    CALL message('',' Memory for nudging weights allocated.')
!EOX
    !--- additional fields for pattern assimilation
    IF (lnudgpat) THEN
       !-- global fields from external data files
!BOX
       ALLOCATE(&
            a_pat_lnp       (4), &
            a_pat_tem(nlevp1,4), &
            a_pat_div(nlev,  4), &
            a_pat_vor(nlev,  4), stat=ierr)
       IF (ierr /= 0) CALL finish('NudgingInit', &
            'Allocation pattern coefficients A_PAT failure.')
!EOX
       !--- norm of pattern
!BOX
       ALLOCATE(&
            a_nrm_lnp       (4), &
            a_nrm_tem(nlevp1,4), &
            a_nrm_div(nlev,  4), &
            a_nrm_vor(nlev,  4), stat=ierr)
       IF (ierr /= 0) CALL finish('NudgingInit', &
            'Allocation pattern coefficients A_NRM failure.')

       !--- local field, corrected each time step
       ALLOCATE(&
            b_pat_lnp     (2), &
            b_pat_tem(nlev,2), &
            b_pat_div(nlev,2), &
            b_pat_vor(nlev,2), stat=ierr)
       IF (ierr /= 0) CALL finish('NudgingInit', &
            'Allocation pattern coefficients B_PAT failure.')

       CALL message('',' Memory for pattern nudging allocated.')
    END IF

    CALL message('','Nudging memory allocated.')
    CALL message('','')

  END SUBROUTINE NudgingAllocMem
!EOX
!EOC
  !=============================================================================
!BOP
  ! !IROUTINE:  NudgingDeallocMem
  ! !INTERFACE:

  SUBROUTINE NudgingDeallocMem

    ! !DESCRIPTION: 
    ! deallocate nudging memory

    ! !USES:
    USE mo_exception,         ONLY: message, finish
    USE mo_nudging_constants, ONLY: lsite, lnudgpat, lnudgwobs
    USE mo_control,           ONLY: lnmi
!EOP
!BOC
!BOX
    INTEGER :: ierr

    DEALLOCATE(flagn, nudgto, nudgdo, nudgvo, stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
         'Error while deallocating flagn, nudgto, nudgdo, nudgvo.')
    DEALLOCATE(nudgta, nudgtb, stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
         'Error while deallocating nudgta, nudgtb.')
    DEALLOCATE(nudgda, nudgdb, stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
         'Error while deallocating nudgda, nudgdb.')
    DEALLOCATE(nudgva, nudgvb, stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
         'Error while deallocating nudgva, nudgvb.')
    DEALLOCATE(nudgsa, nudgsb, stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
         'Error while deallocating nudgsa, nudgsb.')
    CALL message('','  Remove nudging coefficient work space.')

    DEALLOCATE(buf_sfc, buf_ref, stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
         'Error while deallocating buf_sfc, buf_ref.')
    
    IF (lnudgwobs) DEALLOCATE(buf_ano, stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
         'Error while deallocating buf_ano.')
    CALL message('','  Remove TEN, TAC and diagnostic fields.')

    IF (lsite) THEN
      DEALLOCATE(buf_site, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
           'Error while deallocating buf_site.')
      CALL message('','  Remove SITE fields.')
    END IF

    IF (lnmi) THEN
      DEALLOCATE(buf_fm, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
           'Error while deallocating buf_fm.')
      CALL message('','  Remove NMI diagnostic fields.')
    END IF

    IF (lnudgpat) THEN
      DEALLOCATE(a_pat_lnp, a_pat_tem, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
           'Error while deallocating a_pat_lnp, a_pat_tem.')
      DEALLOCATE(a_pat_div, a_pat_vor, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
           'Error while deallocating a_pat_div, a_pat_vor.')
      DEALLOCATE(a_nrm_lnp, a_nrm_tem, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
           'Error while deallocating a_nrm_lnp, a_nrm_tem.')
      DEALLOCATE(a_nrm_div, a_nrm_vor, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
           'Error while deallocating a_nrm_div, a_nrm_vor.')
      DEALLOCATE(b_pat_lnp, b_pat_tem, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
           'Error while deallocating b_pat_lnp, b_pat_tem.')
      DEALLOCATE(b_pat_div, b_pat_vor, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingDeallocMem', &
           'Error while deallocating b_pat_div, b_pat_vor.')
      CALL message('','  Remove pattern nudging fields.')
    END IF

#ifdef MESSY
    ! mz_ht_20050421+
    DEALLOCATE (gp_nudge)
    DEALLOCATE (sp_nudge)   
    DEALLOCATE (outp_idx)
    ! mz_ht_20050421-
#endif
    CALL message('','Nudging memory deallocated.')

  END SUBROUTINE NudgingDeallocMem
!EOX
!EOC

  !=============================================================================
!BOP
  !
  ! !IROUTINE: NudgingAllocRefMem
  ! !INTERFACE:

  SUBROUTINE NudgingAllocRefMem

    ! !DESCRIPTION:
    ! allocate memory for reference data

    ! !USES:
    USE mo_exception,     ONLY: message, finish
    USE mo_control,       ONLY: nlev, nlevp1
!EOP
!BOC
!BOX
    INTEGER :: snsp, ierr
    
    snsp = ldc%snsp

    !--- initialize reference data fields
    ALLOCATE(buf_ref0(nlevp1,2,snsp,3), stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingInit', &
       'Allocation error observed fields 0 failure.')
    stobs0 => buf_ref0(:    ,:,:,1)
    sdobs0 => buf_ref0(:nlev,:,:,2)
    svobs0 => buf_ref0(:nlev,:,:,3)

    ALLOCATE(buf_ref1(nlevp1,2,snsp,3), stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingInit', &
       'Allocation error observed fields 1 failure.')
    stobs1 => buf_ref1(:    ,:,:,1)
    sdobs1 => buf_ref1(:nlev,:,:,2)
    svobs1 => buf_ref1(:nlev,:,:,3)

    ALLOCATE(buf_ref2(nlevp1,2,snsp,3), stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingInit', &
       'Allocation error observed fields 2 failure.')
    stobs2 => buf_ref2(:    ,:,:,1)
    sdobs2 => buf_ref2(:nlev,:,:,2)
    svobs2 => buf_ref2(:nlev,:,:,3)

    ALLOCATE(buf_ref3(nlevp1,2,snsp,3), stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingInit', &
       'Allocation error observed fields 3 failure.')
    stobs3 => buf_ref3(:    ,:,:,1)
    sdobs3 => buf_ref3(:nlev,:,:,2)
    svobs3 => buf_ref3(:nlev,:,:,3)

    buf_ref0(:,:,:,:) = 0.0_dp
    buf_ref1(:,:,:,:) = 0.0_dp
    buf_ref2(:,:,:,:) = 0.0_dp
    buf_ref3(:,:,:,:) = 0.0_dp

    CALL message('',' Memory for observed fields allocated.')

  END SUBROUTINE NudgingAllocRefMem
!EOX
!EOC
  !=============================================================================
!BOP
  !
  ! !IROUTINE:  NudgingDeallocRefMem
  ! !INTERFACE:

  SUBROUTINE NudgingDeallocRefMem

    ! !DESCRIPTION:
    ! clean up reference data buffer
    !
    ! !USES:
    USE mo_exception,         ONLY: message, finish
!EOP
!BOC
!BOX
    INTEGER :: ierr

    NULLIFY(sdobs0, svobs0, stobs0)
    NULLIFY(sdobs1, svobs1, stobs1)
    NULLIFY(sdobs2, svobs2, stobs2)
    NULLIFY(sdobs3, svobs3, stobs3)

    DEALLOCATE(buf_ref0, buf_ref1, buf_ref2, buf_ref3, stat=ierr)
    IF (ierr /= 0) CALL finish('NudgingDeallocRefMem', &
         'Error while deallocating buf_ref0, buf_ref1, buf_ref2, buf_ref3.')
    CALL message('','Remove nudging reference field memory.')

  END SUBROUTINE NudgingDeallocRefMem
!EOX
!EOC
  !=============================================================================
!BOP
  !
  ! !IROUTINE:  NdgInitCounter
  ! !INTERFACE:

  SUBROUTINE NdgInitCounter

    ! !DESCRIPTION:
    ! The procedure initialises internal accumulation counters. The counters 
    ! are coded as the global average of divergence in the rerun files. The
    ! counters are used for additional diagnostics (NMI filter active, 
    ! and/or SITE analysis).
    !
    ! !USES:
    USE mo_exception,         ONLY: message, finish
    USE mo_mpi,               ONLY: p_parallel, p_parallel_io, p_io, p_bcast
    USE mo_nudging_constants, ONLY: lnudgdbx, lsite
    USE mo_control,           ONLY: lnmi
!EOP
!BOC
!BOX

   INTRINSIC int

    indg_accu = 0         !ik temporary

    ! initialize flags and counter
    ifast_accu = 0
    lfill_a = .FALSE.
    lfill_b = .FALSE.

    isiteaccu = 0
    lsite_n0 = .FALSE.
    lsite_n1 = .FALSE.

    ! decode flags from rerun fields
    IF (p_parallel_io .AND. ldc%nsnm0 > 0) THEN
      IF (ldc%snn0(1) == 0) THEN
        IF (lnmi) THEN
          lfill_a = (sdfast_a(1,1,1) < 0.0_dp)
          IF (lfill_a) ifast_accu = - (INT(sdfast_a(1,1,1)) + 1)
          lfill_b = (sdfast_b(1,1,1) < 0.0_dp)
        END IF
        IF (lsite) THEN
          lsite_n0 = (sd_o_n0(1,1,1) < 0.0_dp)
          IF (lsite_n0) isiteaccu  = - (INT(sd_o_n0(1,1,1)) + 1)
          lsite_n1 = (sd_o_n1(1,1,1) < 0.0_dp)
        END IF
      ELSE
        CALL finish('NdgInitCounter','first wave number not available ')
      END IF
    END IF

    IF (p_parallel) THEN
      CALL p_bcast(indg_accu, p_io)

      CALL p_bcast(ifast_accu, p_io)
      CALL p_bcast(lfill_a, p_io)
      CALL p_bcast(lfill_b, p_io)
      
      CALL p_bcast(isiteaccu, p_io)
      CALL p_bcast(lsite_n0, p_io)
      CALL p_bcast(lsite_n1, p_io)
    END IF

    IF (lnudgdbx) CALL message('NdgInitCounter','counter initialised')

  END SUBROUTINE NdgInitCounter
!EOX
!EOC
  !=============================================================================
!BOP
  !
  ! !IROUTINE:  NdgSetCounter
  ! !INTERFACE:

  SUBROUTINE NdgSetCounter

    ! !DESCRIPTION:
    ! code the counter values
    !
    ! !USES:
    USE mo_exception,         ONLY: message, finish
    USE mo_mpi,               ONLY: p_parallel_io
    USE mo_nudging_constants, ONLY: lnudgdbx, lsite
    USE mo_control,           ONLY: lnmi
!EOP
!BOC
!BOX

   INTRINSIC real

    IF (p_parallel_io .AND. ldc%nsnm0 > 0) THEN
      IF (ldc%snn0(1) == 0) THEN
        IF (lnmi) THEN
          IF (lfill_a) sdfast_a(1,1,1) = -( REAL(ifast_accu,dp) + 1.0_dp )
          IF (lfill_b) sdfast_b(1,1,1) = -1.0_dp
        END IF
        IF (lsite) THEN
          IF (lsite_n0) sd_o_n0(1,1,1) = -( REAL(isiteaccu,dp) + 1.0_dp )
          IF (lsite_n1) sd_o_n1(1,1,1) = -1.0_dp
        END IF
      ELSE
        CALL finish('NdgSetCounter','first wave number not available ')
      END IF
    END IF

    IF (lnudgdbx) CALL message('NdgSetCounter','accumulation counter set')
    
  END SUBROUTINE NdgSetCounter
!EOX
!EOC
  !=============================================================================
!BOP
  !
  ! !IROUTINE:  NdgRemoveCounter
  ! !INTERFACE:

  SUBROUTINE NdgRemoveCounter

    ! !DESCRIPTION:
    ! decode the counter values
    !
    ! !USES:
    USE mo_exception,         ONLY: message, finish
    USE mo_mpi,               ONLY: p_parallel_io
    USE mo_nudging_constants, ONLY: lnudgdbx, lsite
    USE mo_control,           ONLY: lnmi
!EOP
!BOC
!BOX
! mz_pj_20050720+
    !IF (p_parallel_io) THEN
    !  IF (ldc%sm(1) == 0 .AND. ldc%snn0(1) == 0) THEN
    IF (p_parallel_io .AND. ldc%nsnm0 > 0) THEN
      IF (ldc%snn0(1) == 0) THEN
! mz_pj_20050720-
        IF (lnmi) THEN
          sdfast_a(1,1,1) = 0.0_dp
          sdfast_b(1,1,1) = 0.0_dp
        END IF

        IF (lsite) THEN
          sd_o_n0(1,1,1) = 0.0_dp
          sd_o_n1(1,1,1) = 0.0_dp
        END IF
      ELSE
        CALL finish('NdgRemoveCounter','first wave number not available ')
      END IF

      IF (lnudgdbx) CALL message('NdgRemoveCounter', &
         'accumulation counter removed')

    END IF

  END SUBROUTINE NdgRemoveCounter
!EOX
!EOC
  !=============================================================================
!BOP
  ! !IROUTINE:  NdgCorrBuffer
  ! !INTERFACE:

  SUBROUTINE NdgCorrBuffer

    ! !DESCRIPTION: 
    ! write nudging diagnostics data

    ! !USES:
    USE mo_nudging_constants, ONLY: lnudgdbx, lnudg_run, lsite, lnudgwobs
    USE mo_exception,         ONLY: message
    !USE mo_time_control,      ONLY: get_interval_seconds, &
    !                                 ev_putdata, &
    !                                 delta_time
    USE mo_time_control,      ONLY: delta_time
    USE mo_control,           ONLY: lnmi, nlev
!EOP
!BOC
!BOX
    REAL(kind=dp)  :: tfact, tfactnmi, tfactsite
    LOGICAL        :: lfaststore, lsitestore
    INTEGER        :: isec

    INTRINSIC int, real

    IF (.NOT. lnudg_run) RETURN

    CALL message('','Store nudging diagnostics ...')

    ! **** proportionality factor, output given in VALUE/sec *******************
    tfact = 1.0_dp

    !    isec  = get_interval_seconds(ev_putdata(nio_index))
    isec = INT(REAL(indg_accu,dp)*delta_time)   !ik temporary

    IF (isec > 0) tfact = 1.0_dp/isec

    lfaststore = ( lnmi .AND. ifast_accu>0 )
    ! rescale the factor
    IF (lfaststore) tfactnmi = 1/(ifast_accu*2.*delta_time)

    lsitestore = lsite .AND. isiteaccu>0
    IF (lsitestore) tfactsite = 1.0_dp/REAL(isiteaccu)

    IF (lnudgdbx) THEN
      WRITE(mess,*) 'correct accumulated DATA with = ',tfact
      CALL message('NudgingOut',mess)
      IF (lfaststore) THEN
        WRITE(mess,*) 'correct fast tendency with = ',tfactnmi
        CALL message('NudgingOut',mess)
      END IF
      IF (lsitestore) THEN
        WRITE(mess,*) 'correct site diagnostics with = ',tfactsite
        CALL message('NudgingOut',mess)
      END IF
    END IF

    sstac(:,:,:) = sstac(:,:,:)*tfact
    sttac(:,:,:) = sttac(:,:,:)*tfact
    sdtac(:,:,:) = sdtac(:,:,:)*tfact
    svtac(:,:,:) = svtac(:,:,:)*tfact

    IF (lnudgwobs) THEN
      asano(:,:,:) = asano(:,:,:)*tfact
      atano(:,:,:) = atano(:,:,:)*tfact
      adano(:,:,:) = adano(:,:,:)*tfact
      avano(:,:,:) = avano(:,:,:)*tfact
    END IF

    IF (lfaststore) THEN
      stfast_accu(:,:,:) = flagn(:,    :,:)*stfast_accu(:,:,:)*tfactnmi
      sdfast_accu(:,:,:) = flagn(:nlev,:,:)*sdfast_accu(:,:,:)*tfactnmi
      svfast_accu(:,:,:) = flagn(:nlev,:,:)*svfast_accu(:,:,:)*tfactnmi
    END IF

    IF (lsitestore) THEN
      stsite(:,:,:) = stsite(:,:,:)*tfactsite
      sdsite(:,:,:) = sdsite(:,:,:)*tfactsite
      svsite(:,:,:) = svsite(:,:,:)*tfactsite
    END IF

#ifdef MESSY
    CALL NUDG_SP2GP
#endif

  END SUBROUTINE NdgCorrBuffer
!EOX
!EOC
  !=============================================================================
!BOP
  ! !IROUTINE:  NdgCleanBuffer
  ! !INTERFACE:

  SUBROUTINE NdgCleanBuffer

    ! !DESCRIPTION: 
    ! clean up accumulated arrays

    ! !USES:
    USE mo_nudging_constants, ONLY: lnudg_run, lsite, lnudgwobs
    USE mo_control,           ONLY: lnmi
!EOP
!BOC
!BOX

    IF (.NOT. lnudg_run) RETURN

    ! Remark: ??TAC and ??ANO can handled with the laccu flag of the
    ! stream interface, but it will work until the revision of the
    ! time control not well for the first interval in the case, where
    ! nuging will start between two stream output dates
    ! temporary fix: local counter in nudging, will fail with changes
    ! of delta_time

    buf_sfc(:,:,:,2) = 0.0_dp
    sttac(:,:,:) = 0.0_dp
    sdtac(:,:,:) = 0.0_dp
    svtac(:,:,:) = 0.0_dp

    IF (lnudgwobs) THEN
      buf_ano(:,:,:,2) = 0.0_dp
      atano(:,:,:) = 0.0_dp
      adano(:,:,:) = 0.0_dp
      avano(:,:,:) = 0.0_dp
    END IF
    indg_accu = 0

    IF (lnmi) THEN
      buf_fm(:,:,:,:) = 0.0_dp
      ifast_accu  = 0
    END IF
    IF (lsite) THEN
      buf_site(:,:,:,:) = 0.0_dp
      isiteaccu  = 0
    END IF

  END SUBROUTINE NdgCleanBuffer
!EOX
!EOC
  !================================================================ 
#ifdef MESSY
  ! mz_ht_20050421+
  SUBROUTINE Nudg_SP2GP

    USE messy_main_transform_bi,    ONLY: sp2gp

    INTEGER :: i

    DO i=1,output_num
      CALL SP2GP(sp_nudge(outp_idx(i))%nudge_spfields  &
           , gp_nudge(i)%nudge_gpfields, .false.)
    ENDDO

  END SUBROUTINE Nudg_SP2GP

!============================================================================
  ! mz_ht_20050421-
#endif
END MODULE mo_nudging_buffer



