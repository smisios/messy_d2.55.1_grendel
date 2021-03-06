#include "messy_main_ppd_bi.inc"

!NOTE: test performance: tendency for all tracers in one step or individal
#define ALLATONCE
!*****************************************************************************
!*****************************************************************************

! Authors:
! Rolf Sander,     MPICH, 2001-2015
! Astrid Kerkweg,  MPICH, 2002-2004
! Patrick Joeckel, MPICH Mar 2004:  strict separation of AERO from MECCA
! Patrick Joeckel, MPICH May 2005:  khet*_3d and j_*_3d converted to
!                                   PTR_3D_ARRAY and 1D-ARRAY, respectively
! Andrea Pozzer,   MPICH Sep 2005:  Lagrangian extension
! Astrid Kerkweg,  UNI_MZ, 2009:    e5 => si

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

MODULE messy_mecca_si

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
    info_bi, warning_bi, error_bi
  USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp, ntrac_lg, ti_lg
#ifdef MESSYTENDENCY
  !tendency budget
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
    mtend_get_start_l,      &
    mtend_add_l,            &
    mtend_register,         &
    mtend_id_tracer,        &
    mtend_id_q                 ! op_pj_20160510
#endif
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_VLONG, M_air, M_H2O
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY, PTR_1D_ARRAY, str

  ! SUBMODEL MECCA
  USE messy_mecca_mem_si     ! jk, jp, idt_*
  USE messy_mecca            ! l_aero, l_force_khet, l_kpp_debug, modstr,
                             ! kpromanlev_to_nbl, nbl_to_nblx, nblx_to_nbl
  USE messy_mecca_khet,      ONLY: l_troposphere, l_stratosphere
  USE messy_mecca_kpp, ONLY: &
    ! Several variables could be imported from messy_mecca_kpp or from
    ! messy_mecca002_kpp etc. In some cases, they are imported from
    ! messy_mecca_kpp because they are the same in all messy_mecca*_kpp
    ! files. In other cases, diagnostic output is only produced for the
    ! first mechanism and messy_mecca002_kpp etc. are ignored.    
    ! same because from gen.c in kpp:
    DP,                  & !
    ! same because from gas.eqn:
    IHT_MAX,             & ! update_physc_gp and _lg
    IHS_MAX,             & ! update_physc_gp and _lg
    ! same because from messy_cmn_photol_mem:
    IP_MAX,              & ! for declaration of photrat_gp and _lg
    jname,               & ! mecca_init_coupling_gp and _lg
    ! same because it must be "4":
    kppoption,           & ! mecca_initialize
    ! different but diagnostic output only shows first mechanism:
    IERR_NAMES,          & ! mecca_physc
    rplfile,             & ! mecca_init_memory
    wanted,              & ! mecca_init_memory
    diagtracfile,        & ! mecca_init_memory
    gas_eqn_file,        & ! mecca_init_memory
    batchfile              ! mecca_init_memory
  USE messy_mecca_poly_si, ONLY: &
    NMAXMECCA,           &
    ! value depends on all mechanisms:
    REQ_HET,             & ! mecca_initialize
    REQ_PHOTRAT            ! mecca_init_coupling_gp and _lg
  ! SUBSUBMODELS
  USE messy_mecca_aero_si    ! mecca_aero_* interface level subroutines
  USE messy_mecca_khet_si    ! aerochem_cpl, khet_[Tr/St]_3d, mecca_khet_*, ...
#ifdef MECCA_TAG
  USE messy_mecca_tag_si     ! tagging process routines
#endif

  IMPLICIT NONE
  INTRINSIC :: NULL

  CHARACTER(LEN=STRLEN_MEDIUM) :: photrat_channel_gp = ''
  CHARACTER(LEN=STRLEN_MEDIUM) :: photrat_channel_lg = ''
  LOGICAL :: l_khet
  LOGICAL :: l_gp = .FALSE.
  LOGICAL :: l_lg = .FALSE.
  LOGICAL :: lcheck_range = .FALSE.
  LOGICAL, SAVE :: l_skipkpp_gp = .FALSE.   ! switch for skip chemistry for gp
  LOGICAL, SAVE :: l_skipkpp_lg = .FALSE.   ! switch for skip chemistry for lg
  LOGICAL, SAVE :: l_polymecca ! use only one MECCA chemistry mechanism?
  ! Which aerosol submodel delivers the properties (radius, sigma) for
  ! the Pseudo Aerosol tracers (see gas.tex) ...
  CHARACTER(LEN=STRLEN_MEDIUM), SAVE :: c_pa_asm = ''
  ! ... and which mode ?
  INTEGER, SAVE :: i_pa_amode = 0
  REAL(DP), PARAMETER :: vtmpc1 = M_air / M_H2O - 1._dp
  ! op_pj_20160510+
  INTEGER :: i_H2O_feedback = 0 ! no feedback to hydrological cycle
  !                             ! (1: feedback from GP, 2: feedback from LG)
  LOGICAL :: lfq_gp = .FALSE.   ! feedback calculation for spec. humidity (GP)
  LOGICAL :: lfq_lg = .FALSE.   ! feedback calculation for spec. humidity (LG)
  ! op_pj_20160510-
  ! op_pj_20161220+
  LOGICAL :: l_sync_H2O_gp = .FALSE. ! synchronize H2O tracer with specific 
  !                                  ! humidity (q)
  LOGICAL :: l_sync_H2O_lg = .FALSE. ! ...
  ! op_pj_20161220-

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
  INTEGER :: my_handle_corr
#endif

  NAMELIST /CPL/ lcheck_range, photrat_channel_gp, photrat_channel_lg, &
    l_lg, l_gp, l_skipkpp_gp, l_skipkpp_lg, &
    c_pa_asm, i_pa_amode, &
    i_H2O_feedback ! op_pj_20160510

  PRIVATE
  PUBLIC :: mecca_initialize    ! initialize aero chemistry
  PUBLIC :: mecca_new_tracer    ! define tracers
  PUBLIC :: mecca_init_memory   ! new channels and initialize kpp variables
  PUBLIC :: mecca_init_tracer   ! initialize tracers
  PUBLIC :: mecca_init_coupling ! coupling to echam5
  PUBLIC :: mecca_vdiff
  PUBLIC :: mecca_physc         ! calculate chemistry
  PUBLIC :: mecca_global_end    ! calculate chemistry
  PUBLIC :: mecca_free_memory

  REAL(DP), DIMENSION(:,:,:), POINTER, SAVE :: out_nsteps_gp
  REAL(DP), DIMENSION(:),     POINTER, SAVE :: out_nsteps_lg

  INTEGER :: zidt_H2O_gp
  INTEGER :: zidt_H2O_lg

  ! LOWER ERROR BOUNDARY FOR CHECK_RANGE
  REAL(DP), PARAMETER :: negerr = -1.e-15_dp

  ! PRESS_PARCEL(NCELLS) = pressure coordinate of the parcel
  REAL(DP), POINTER, DIMENSION(:) :: PRESS_PARCEL
  ! 2d field for density  transformation
  REAL(DP), DIMENSION(:,:,:), POINTER :: tmp_temp_gp => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: tmp_sphum_gp => NULL()
  REAL(DP), DIMENSION(:), POINTER :: tmp_temp_lg => NULL()
  REAL(DP), DIMENSION(:), POINTER :: tmp_sphum_lg => NULL()
  ! op_pj_20160511+
  REAL(DP), DIMENSION(:),     POINTER :: qtend_lg => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: qtend_gp => NULL()
  ! op_pj_20160511-

  ! pointer to channel objects:
  TYPE(PTR_3D_ARRAY), PUBLIC, DIMENSION(IP_MAX), SAVE :: photrat_gp
  TYPE(PTR_1D_ARRAY), PUBLIC, DIMENSION(IP_MAX), SAVE :: photrat_lg
  TYPE(PTR_3D_ARRAY), PUBLIC,                    SAVE :: meccanum_gp

  !mz_bs_20150702+
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER, SAVE :: mtskiptend
  REAL(DP), DIMENSION(:,:,:), POINTER, SAVE :: mtskiptend_q => NULL() ! op_pj_20160511
  LOGICAL, DIMENSION(:), POINTER, SAVE    :: l_trac_mtskip => NULL()
  LOGICAL, SAVE                           :: lmtskip = .FALSE.
  REAL(DP), POINTER, SAVE                 :: rmode_mtskip    => NULL()
  REAL(DP), POINTER, SAVE                 :: tsl_phys_mtskip => NULL()
  !mz_bs_20150702-

CONTAINS

  !***************************************************************************

  SUBROUTINE mecca_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi, ONLY: p_parallel_io, p_io, p_bcast
    ! MESSY
    USE messy_main_tools,  ONLY: find_next_free_unit
    USE messy_mecca_poly_si, ONLY: initialize_poly

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    CALL start_message_bi(modstr,'INITIALISATION',modstr)

    CALL initialize_poly ! this also reads CTRL_KPP 

    IF (kppoption /= '4') THEN
      ! ECHAM5/MESSy - MECCA only usable with kp4
      ! Reason: Native KPP variables (temp, press, ...) are scalars.
      !         KP4 converts them into arrays of rank 1,
      !         fills these with all grid-boxes of one PE,
      !         and puts the loop over grid-boxes as innermost loop of the
      !         integrator.
      !         Thus, native KPP would require a SMIL with fill-routines
      !         etc. inside an outermost loop over grid-boxes.
      !
      CALL error_bi('Sorry. Please recode with kp4 and recompile.', substr)
    ENDIF

#ifdef MECCA_TAG
    ! modify default; note that running a TAGged meccanism with l_tag=.FALSE.
    ! is a rare exception ...
    l_tag = .TRUE.
#endif

    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL mecca_read_nml_ctrl(status, iou)
      IF (status /= 0) CALL error_bi('error in mecca_read_nml_ctrl',substr)
    ENDIF
    ! BROADCAST RESULTS
    CALL p_bcast(l_aero,       p_io)
    CALL p_bcast(l_force_khet, p_io)
    CALL p_bcast(l_kpp_debug,  p_io)
    CALL p_bcast(l_tag,        p_io)

    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL mecca_read_nml_cpl(status, iou)
      IF (status /= 0) CALL error_bi('error in mecca_read_nml_cpl',substr)
    ENDIF
    CALL p_bcast(lcheck_range,       p_io)
    CALL p_bcast(photrat_channel_gp, p_io)
    CALL p_bcast(photrat_channel_lg, p_io)
    CALL p_bcast(l_gp,               p_io)
    CALL p_bcast(l_lg,               p_io)
    CALL p_bcast(l_skipkpp_gp,       p_io)
    CALL p_bcast(l_skipkpp_lg,       p_io)
    CALL p_bcast(c_pa_asm,           p_io)
    CALL p_bcast(i_pa_amode,         p_io)
    CALL p_bcast(i_H2O_feedback,     p_io) ! op_pj_20160510

    l_khet = l_force_khet .OR. REQ_HET
    CALL info_bi('l_force_khet = '//str(l_force_khet))
    CALL info_bi('REQ_HET      = '//str(REQ_HET))
    CALL info_bi('l_khet       = '//str(l_khet))
    IF (l_khet) CALL mecca_khet_initialize

    IF (l_aero) CALL mecca_aero_initialize

    CALL end_message_bi(modstr,'INITIALISATION',modstr)

  END SUBROUTINE mecca_initialize

  !***************************************************************************

  SUBROUTINE mecca_new_tracer

  USE messy_mecca_poly_si, ONLY: mecca_new_tracer_gp

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_new_tracer'

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#ifdef MECCA_TAG
    my_handle_corr = mtend_get_handle(modstr//'_corr')
#endif
#endif

    CALL start_message_bi(modstr, 'TRACER INITIALIZATION', substr)
    IF (l_gp) CALL mecca_new_tracer_gp(c_pa_asm, i_pa_amode)
#if defined(ECHAM5)
    IF (l_lg) CALL mecca_new_tracer_lg
#endif

    CALL end_message_bi(modstr, 'TRACER INITIALIZATION', substr)

  END SUBROUTINE mecca_new_tracer

  !***************************************************************************

#if defined(ECHAM5)
  SUBROUTINE mecca_new_tracer_lg

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR, LGTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, AIR, ON, OFF, &
      set_tracer, &
      AEROSOL, AMOUNTFRACTION, &
      MODAL, BIN, &
      I_ADVECT, I_CONVECT, &
      I_VDIFF, &
      I_DRYDEP, I_SEDI, &
      I_SCAV, I_MIX, &
      I_FORCE_COL, I_INTEGRATE,&
      I_TIMEFILTER, I_FORCE_INIT, &
      I_AEROSOL_METHOD, I_AEROSOL_MODE,&
      I_AEROSOL_SOL, S_AEROSOL_MODEL, &
      R_MOLARMASS, &
      R_DRYREAC_SF, R_VINI, &
      R_AEROSOL_DENSITY
    USE messy_main_constants_mem, ONLY: MH, MC, MN, MNa, MO, MS, MCl, MBr, &
      MI, MF, MHg

    ! When using a small (e.g. gas-phase only) chemistry mechanism,
    ! forcheck will say
    ! that AEROSOL, AMOUNTFRACTION, LGTRSTR, ON, and OFF are not used
    ! ("named constant not used"). Ignore this info and leave these
    ! variables in the ONLY list.

    IMPLICIT NONE

    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_new_tracer_lg'
    CHARACTER(LEN=*), PARAMETER :: setname = LGTRSTR

    CALL start_message_bi(modstr, 'LAGRANGIAN TRACER REQUEST', substr)
    ! The following INCLUDE statement contains all "CALL new_tracer" commands
    ! for all species that are used in the kpp chemistry scheme.
    ! Note that MECCA will create a tracer for H2O only, if MECCA_TAG is
    ! active and the mechanism contains H2O.
    INCLUDE 'messy_mecca_trac_si.inc'
    CALL end_message_bi(modstr, 'LAGRANGIAN TRACER REQUEST', substr)

  END SUBROUTINE mecca_new_tracer_lg
#endif

  !***************************************************************************

  SUBROUTINE mecca_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: new_attribute
    USE messy_main_tools,            ONLY: strcrack
    USE messy_mecca_poly_si,         ONLY: py_initialize
#ifdef MECCA_TAG
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer,           ONLY: get_tracer
#endif

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_init_memory'
    INTEGER                     :: status
    INTEGER                     :: n
    CHARACTER(LEN=STRLEN_VLONG), DIMENSION(:), POINTER :: el => NULL()
#ifdef MESSYTENDENCY
    INTEGER                     :: idt

    CALL mtend_register(my_handle, mtend_id_tracer)
! op_pj_20160511+
    IF (i_H2O_feedback > 0) THEN
       CALL mtend_register(my_handle, mtend_id_q)
    ENDIF
! op_pj_20160511-
#ifdef MECCA_TAG
    CALL get_tracer(status, GPTRSTR, 'H2O', idx=idt)
    IF (status == 0) then
!       CALL mtend_register(my_handle_corr, mtend_id_tracer, idt=idt)
       CALL mtend_register(my_handle_corr, idt)
    ENDIF
#endif
#endif

    CALL start_message_bi(modstr, 'GLOBAL ATTRIBUTES', substr)
    CALL strcrack(gas_eqn_file, ' ', el, n)
    IF (n > 0) THEN
      CALL new_attribute(status, 'MESSy_'//modstr//'_eqn', c= TRIM(el(n)))
      CALL channel_halt(substr, status)
      IF (ASSOCIATED(el)) DEALLOCATE(el); NULLIFY(el)
    ENDIF
    CALL new_attribute(status, 'MESSy_'//modstr//'_bat', c= batchfile)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'MESSy_'//modstr//'_rpl', c= rplfile)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'MESSy_'//modstr//'_nism', c= wanted)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'MESSy_'//modstr//'_diag', c= diagtracfile)
    CALL channel_halt(substr, status)
    CALL end_message_bi(modstr, 'GLOBAL ATTRIBUTES', substr)

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)
    IF (l_khet) CALL mecca_khet_init_memory(L_LG, L_GP)
    IF (l_gp) CALL mecca_init_memory_gp
#if defined(ECHAM5)
    IF (l_lg) CALL mecca_init_memory_lg
#endif
    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL start_message_bi(modstr, 'INITIALIZE KPP VARIABLES', substr)
    CALL py_initialize
#ifdef MECCA_TAG
    IF (l_tag) CALL mecca_tag_init  ! call mecca_tag init on SI level
#endif
    CALL end_message_bi(modstr, 'INITIALIZE KPP VARIABLES', substr)

  END SUBROUTINE mecca_init_memory

  !***************************************************************************

  SUBROUTINE mecca_init_memory_gp

    USE messy_main_tools,            ONLY: strcrack, match_wild
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
      new_attribute

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_init_memory_gp'
    INTEGER :: status

    CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr//'_gp', 'nsteps', &
      p3=out_nsteps_gp)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_gp', 'nsteps', &
      'long_name', c='number of KPP substeps')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_gp', 'nsteps', &
      'units', c=' ')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object '//modstr//'_gp/'// &
      'nsteps'//' was created')

    !-------------------------------------------------------------------------

    ! define and/or modify output channels for aero
    IF (l_aero) CALL mecca_aero_init_memory

  END SUBROUTINE mecca_init_memory_gp

  !***************************************************************************
#if defined(ECHAM5)
  SUBROUTINE mecca_init_memory_lg

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: LG_ATTILA
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, ngpblks, nlev
    USE messy_main_tracer_mem_bi,    ONLY: NCELL
    ! MESSy
    USE messy_main_tools,        ONLY: strcrack, match_wild
    USE messy_main_channel,      ONLY: new_channel, new_channel_object, &
      new_attribute


    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_init_memory_lg'
    INTEGER :: status

    ALLOCATE(tmp_temp_gp(nproma,nlev,ngpblks))
    ALLOCATE(tmp_sphum_gp(nproma,nlev,ngpblks))
    ALLOCATE(tmp_temp_lg(NCELL))
    ALLOCATE(tmp_sphum_lg(NCELL))
    ! op_pj_20160511+
!!$    IF (lfq_lg) THEN  ! lfq_lg is set later in init_tracer
       ALLOCATE(qtend_lg(NCELL))
       ALLOCATE(qtend_gp(nproma,nlev,ngpblks))
!!$    END IF
    ! op_pj_20160511-

    CALL new_channel(status, modstr//'_lg', reprid=LG_ATTILA)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr//'_lg', 'nsteps', &
      p1=out_nsteps_lg)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_lg', 'nsteps', &
      'long_name', c='number of KPP substeps')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_lg', 'nsteps', &
      'units', c=' ')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object '//modstr//'_lg/'// &
      'nsteps'//' was created')

    !-------------------------------------------------------------------------

    ! define and/or modify output channels for aero
    !   IF (l_aero) CALL mecca_aero_init_memory

  END SUBROUTINE mecca_init_memory_lg
#endif
  !***************************************************************************

  SUBROUTINE mecca_init_coupling(iflag) ! mz_bs_20150702: iflag

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)         :: iflag   ! mz_bs_20150702

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_init_coupling'

    SELECT CASE(iflag) ! mz_bs_20150702

       CASE(1)         ! mz_bs_20150702 
          CALL start_message_bi(modstr, 'COUPLING', substr)
          IF (l_khet .AND. .NOT. l_skipkpp_gp)  &
               CALL mecca_khet_init_coupling
          IF (l_gp) CALL mecca_init_coupling_gp(iflag) ! mz_bs_20150702
#if defined(ECHAM5)
          IF (l_lg) CALL mecca_init_coupling_lg
#endif
#ifdef MECCA_TAG
          IF (l_tag) CALL mecca_tag_init_coupling
#endif
          CALL end_message_bi(modstr, 'COUPLING', substr)

       CASE(2) ! mz_bs_20150702+
          CALL start_message_bi(modstr, 'COUPLING PART 2', substr)
          IF (l_gp) CALL mecca_init_coupling_gp(iflag)
          CALL end_message_bi(modstr, 'COUPLING PART 2', substr)

       END SELECT ! mz_bs_20150702-

  END SUBROUTINE mecca_init_coupling

  !***************************************************************************
  
  SUBROUTINE mecca_init_coupling_gp(iflag) !mz_bs_20150702: iflag

    ! MESSy
    USE messy_main_channel, ONLY: get_channel_object, get_channel_info

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)         :: iflag  !mz_bs_20150702: iflag

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_init_coupling_gp'
    INTEGER           :: status
    INTEGER           :: ip
    REAL(DP), POINTER :: rmecca => NULL()     !mz_bs_20150702

    CALL start_message_bi(modstr, 'TEST COUPLING',substr)

    SELECT CASE (iflag) !mz_bs_20150702
    CASE(1)             !mz_bs_20150702

      ! ---------------------------------------------------------------------

      ! coupling to CHEMGLUE:
      CALL get_channel_info(status, 'chemglue_gp')
      ! if chemglue_gp channel does exist, use several MECCA mechanisms:
      l_polymecca = (status == 0)
      IF (l_polymecca) THEN
        CALL get_channel_object(status, 'chemglue_gp', 'meccanum', &
          p3=meccanum_gp%ptr)
        IF (status /= 0) CALL error_bi &
          ('channel object meccanum does not exist in chemglue_gp', substr)
        IF (l_aero) CALL error_bi( &
          'polymecca does not work with l_aero (yet?).', substr)
        IF (l_tag) CALL error_bi( &
          'polymecca does not work with l_tag (yet?).', substr)
      ENDIF

      ! ---------------------------------------------------------------------

      ! coupling to aerosols:
      IF (l_aero)  CALL mecca_aero_init_coupling

      ! ---------------------------------------------------------------------

      ! coupling to photolysis:
      IF (REQ_PHOTRAT .AND. .NOT. l_skipkpp_gp) THEN
        CALL get_channel_info(status, TRIM(photrat_channel_gp))
        IF (status /= 0) THEN
          CALL info_bi('current MECCA scheme requires the channel '&
            &//TRIM(photrat_channel_gp))
          CALL error_bi('Switch on the corresponding submodel',substr)
        ENDIF
        DO ip=1, IP_MAX
          CALL get_channel_object(status, TRIM(photrat_channel_gp), &
            'J_'//TRIM(jname(ip)), p3=photrat_gp(ip)%ptr)
          IF (status /= 0) &
            CALL warning_bi( 'J_'//TRIM(jname(ip))//&
            ' not in channel '//TRIM(photrat_channel_gp),substr)
        ENDDO
      ENDIF

      ! ---------------------------------------------------------------------

      !mz_bs_20150702+
      !Logic for MTSKIP is set here
      lmtskip = .FALSE.
      CALL get_channel_object(status,'mtskip', 'rmecca',p0=rmecca)
      lmtskip = (status == 0)
      ! Note: not with AND in same line, because rmecca might not be associated
      IF (lmtskip) THEN
        lmtskip = (NINT(rmecca)==1)
      ENDIF
      !mz_bs_20150702-

      ! ---------------------------------------------------------------------

    CASE(2) !mz_bs_20150702+
      ! NOTHING TO DO, IF MTSKIP IS OFF
    END SELECT !mz_bs_20150702-

    IF (lmtskip) CALL mecca_mtskip_init_coupling(iflag)

    CALL end_message_bi(modstr, 'TEST COUPLING',substr)

  END SUBROUTINE mecca_init_coupling_gp

  !***************************************************************************

#if defined(ECHAM5)
  SUBROUTINE mecca_init_coupling_lg

    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info
    USE messy_main_channel_error_bi, ONLY: channel_halt

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_init_coupling_lg'
    INTEGER          :: status
    INTEGER          :: ip

    !       IF (l_aero) CALL mecca_aero_init_coupling

    CALL start_message_bi(modstr, 'TEST COUPLING',substr)

    ! REQ_HET generated by xmecca!
    ! NO HET CHEMISTRY FOR THE MOMENT!!!!!

    ! mz_rs_20150223+
    ! IF ((REQ_AEROSOL).AND.(.NOT.l_aero) .AND. .NOT. l_skipkpp_lg) &
    !   CALL error_bi('current MECCA scheme requires l_aero=T',substr)
    ! this cannot happen anymore because mecca_aero="OFF" is not allowed
    ! anymore, see mecca_read_nml_ctrl in messy_mecca.f90
    ! mz_rs_20150223-

    IF (REQ_PHOTRAT .AND. .NOT. l_skipkpp_lg) THEN
      CALL get_channel_info(status, TRIM(photrat_channel_lg))
      IF (status /= 0) THEN
        CALL info_bi('current MECCA scheme requires the channel '// &
          TRIM(photrat_channel_lg))
        CALL error_bi('Switch on the corresponding submodel',substr)
      ENDIF
      DO ip=1, IP_MAX
        CALL get_channel_object(status, TRIM(photrat_channel_lg), &
          'J_'//TRIM(jname(ip)), p1=photrat_lg(ip)%ptr)
        IF (status /= 0) &
          CALL warning_bi( 'J_'//TRIM(jname(ip))//&
          ' not in channel '//TRIM(photrat_channel_lg),substr)
      ENDDO
    ENDIF

    !CALL info_bi('import information from attila submodel')
    CALL get_channel_object(status, cname='attila', oname='PPRESS' &
      , p1=PRESS_PARCEL )
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr, 'TEST COUPLING',substr)

  END SUBROUTINE mecca_init_coupling_lg
#endif
  !***************************************************************************

  SUBROUTINE mecca_init_tracer

    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, xt, LGTRSTR, xt_a
    USE messy_main_timer,         ONLY: lstart
    USE messy_main_tracer,        ONLY: get_tracer

    IMPLICIT NONE

    INTRINSIC MAX, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_init_tracer'
    INTEGER                     :: status   ! status flag
    INTEGER                     :: jt       ! tracer index

    CALL start_message_bi(modstr, 'MECCA tracer initialization', substr)

    IF (L_GP) THEN
      ! check if H2O tracer exists (if yes, it is probably from H2O submodel)
      CALL get_tracer(status, GPTRSTR, 'H2O',idx=zidt_H2O_gp)
      IF (status == 0) THEN
        CALL info_bi('Using H2O tracer for MECCA chemistry.', substr)
        ! op_pj_20160510+
        IF (i_H2O_feedback == 1) THEN
#ifdef MECCA_TAG
           CALL warning_bi('MECCA_TAG: Forcing feedback to specific humidity (q), although H2O tracer is present.', substr)
           lfq_gp = .TRUE.
#else
           CALL warning_bi('Feedback to specific humidity (q) is disabled, if H2O tracer is present.', substr)
           lfq_gp = .FALSE.
#endif
        ENDIF
        ! op_pj_20160510-
      ELSE
        CALL info_bi('No H2O tracer. Using specific humidity (q) instead.', substr)
        zidt_H2O_gp = 0
        ! op_pj_20160510+
        IF (i_H2O_feedback == 1) THEN
           CALL info_bi('Feedback to specific humidity (q) is enabled.', substr)
           lfq_gp = .TRUE.
        ELSE
           lfq_gp = .FALSE.
        END IF
        ! op_pj_20160510-
      ENDIF
    ENDIF

    ! op_pj_20161220+
    ! ---------------------------------------------------------------------
#ifdef MECCA_TAG
    ! Note that in case of MECCA_TAG the H2O tracer (if required) is
    ! defined by MECCA.
    l_sync_H2O_gp = (zidt_H2O_gp > 0)
#endif
    ! ---------------------------------------------------------------------
    ! op_pj_20161220-

    IF (L_LG) THEN
      ! check if H2O tracer exists (if yes, it is probably from H2O submodel)
      CALL get_tracer(status, LGTRSTR, 'H2O',idx=zidt_H2O_lg)
      IF (status == 0) THEN
        CALL info_bi('Using H2O tracer for MECCA chemistry.', substr)
        ! op_pj_20160510+
        IF (i_H2O_feedback == 2) THEN
#ifdef MECCA_TAG
           CALL warning_bi('MECCA_TAG: Forcing feedback to specific humidity (q), although H2O tracer is present.', substr)
        lfq_lg = .TRUE.
#else
           CALL warning_bi('Feedback to specific humidity is disabled, if H2O tracer is present.', substr)
#endif
        ENDIF
        lfq_lg = .FALSE.
        ! op_pj_20160510-
      ELSE
        CALL info_bi('No H2O tracer. Using specific humidity (q) instead.', substr)
        zidt_H2O_lg = 0
        ! op_pj_20160510+
        IF (i_H2O_feedback == 2) THEN
           CALL info_bi('Feedback to specific humidity (q) is enabled.', substr)
           lfq_lg = .TRUE.
        ELSE
           lfq_lg = .FALSE.
        END IF
        ! op_pj_20160510-
      ENDIF
    ENDIF

    ! op_pj_20161220+
#ifdef MECCA_TAG
    ! Note that in case of MECCA_TAG the H2O tracer (if required) is
    ! defined by MECCA.
    l_sync_H2O_lg = (zidt_H2O_lg > 0)
#endif
    ! op_pj_20161220-

    ! op_pj_20160511+
    ! NOTE: according to the logic above, lfq_gp and lfq_lg cannot be
    !       .TRUE. at the same time, since i_H2O_feedback can either be 1 or 2
    ! op_pj_20160511-

    ! set negative values of mecca tracers to zero
    IF (lstart) THEN ! don't do this after a restart !

      IF (L_GP) THEN
        gp_tracer_loop: DO jt = 1, ntrac_gp
          IF (TRIM(ti_gp(jt)%tp%ident%submodel) == 'mecca') THEN
            xt(_RI_XYZN_(:,:,:,jt)) = MAX(xt(_RI_XYZN_(:,:,:,jt)),0._dp)
          ENDIF
        ENDDO gp_tracer_loop
      ENDIF

      IF (L_LG) THEN
        lg_tracer_loop: DO jt = 1, ntrac_lg
          IF (TRIM(ti_lg(jt)%tp%ident%submodel) == 'mecca') THEN
            xt_a(:,1,jt,1) = MAX(xt_a(:,1,jt,1),0._dp)
          ENDIF
        ENDDO lg_tracer_loop
      ENDIF

    ENDIF

    CALL end_message_bi(modstr, 'MECCA tracer initialization', substr)

  END SUBROUTINE mecca_init_tracer

  !***************************************************************************

  SUBROUTINE mecca_vdiff(flag)

    INTEGER, INTENT(IN) :: flag

    IF (l_aero) CALL mecca_aero_vdiff(flag)

  END SUBROUTINE mecca_vdiff

  !***************************************************************************

  SUBROUTINE mecca_physc

    ! This subroutine is the interface between ECHAM5 and MECCA inside the
    ! time loop. Before "CALL kpp_integrate", several variables (which
    ! are declared in messy_mecca_kpp_g_mem.f90) must be updated:
    !
    ! conc(NSPEC)   = concentrations of NSPEC chemical species [mcl/cm3]
    ! temp       = temperature [K]
    ! press      = pressure [Pa]
    ! cair       = c(air) [mcl/cm^3]
    ! jx(ip_*)   = several (1st order) photolysis rate coefficients [1/s]
    ! khet_Tr(iht_*) = troposheric het. rate coefficients [1/s]
    ! khet_St(ihs_*) = stratospheric het. rate coefficients [cm3/(mcl*s)]
    !
    ! for the subsubmodel mecca_aero, some additional variables are required
    ! for the accumulation soluble (AS) and the coarse soluble (CS) mode:
    !
    ! xaer(AS:CS)        = aerosol in mode AS/CS exists [yes=1/no=0]
    ! cvfac(AS:CS)       = conversion factor from [L(aq)/mol] to [cm^3(air)/mcl]
    ! lwc(AS:CS)         = liquid water content [m3/m3]
    ! k_exf(AS:CS,NSPEC) = exchange rate coefficients, forward [1/s]
    ! k_exb(AS:CS,NSPEC) = exchange rate coefficients, backwards [1/s]
    ! k_exf_*(AS:CS)     = several special exchange rate coefficients [1/s]


    USE messy_main_timer,            ONLY: time_step_len
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, jrow, kproma
#if defined(ECHAM5) || defined(CESM1)
    USE messy_main_grid_def_bi,     ONLY:ilon, ilat
#endif
    USE messy_main_grid_def_bi,     ONLY: philon_2d, philat_2d
    USE messy_main_data_bi,         ONLY: press_3d, tm1_3d, tte_3d,  &
                                          qm1_3d, qte_3d
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_tracer_mem_bi,    ONLY: pxtte => qxtte, pxtm1 => qxtm1
    USE messy_main_tracer_mem_bi,    ONLY: ntrac_gp, ti_gp
    USE messy_main_constants_mem,    ONLY: N_A, R_gas, M_air, M_H2O, OneDay
    USE messy_main_data_bi,          ONLY: qm1
    ! MECCA
    USE messy_mecca_poly_si, ONLY: &
      mr2c, c2mr, kpp_integrate, NSPEC, l_fixed_step, &
      fill_cair, fill_jx, fill_khet_St, fill_khet_Tr, fill_press, fill_temp

    IMPLICIT NONE

    INTRINSIC MAX, MOD, INT, ASSOCIATED

    !mz_bs_20150702+
    INTEGER              :: isw_mtskip
    REAL(DP)             :: ztmst
    !mz_bs_20150702-

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_physc'
    REAL(DP), PARAMETER :: scvmr = M_air/M_H2O ! op_pj_20160510
    INTEGER  :: jt ! tracer loop index
    INTEGER  :: jpm ! polymecca loop index
    INTEGER  :: ntrac
    ! tracer mixing ratio:
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zmrbc ! bc = before chemistry
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zmrac ! ac = after chemistry
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zmrbc_nblx ! bc = before chemistry
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zmrac_nblx ! ac = after chemistry
#ifdef MESSYTENDENCY
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zmr_tmp_3d
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: zmr_tmp_2d
#endif
    ! op_pj_20160510+
    ! FOR WATER VAPOUR FEEDBACK, IF NO H2O TRACER IS PRESENT
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: zwvac
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: zwvac_nblx
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: zwvmrac_nblx
    ! op_pj_20160510-

    ! DEFINE 1D temp, press, cair here, as they are needed in
    !   mecca_aero_update_physc too.
    !   (see also important notes about temp, press, and cair in gas.eqn!)
    REAL(DP), DIMENSION(:), ALLOCATABLE :: temp_nbl
    REAL(DP), DIMENSION(:), ALLOCATABLE :: press_nbl
    REAL(DP), DIMENSION(:), ALLOCATABLE :: cair_nbl
    REAL(DP), DIMENSION(:), ALLOCATABLE :: cair_nblx
    REAL(DP), DIMENSION(:), ALLOCATABLE :: sphum_nbl  ! m(H2O)/m(air) [kg/kg]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: sphum_nblx
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: conc
    ! local, 1-dimensional, integer and logical version of meccanum:
    INTEGER, DIMENSION(:), ALLOCATABLE :: meccanum_nbl
    LOGICAL, DIMENSION(:), ALLOCATABLE :: l_meccanum_nbl

    ! info from KPP integration
    INTEGER, DIMENSION(:), ALLOCATABLE :: xNacc_nbl
    INTEGER, DIMENSION(:), ALLOCATABLE :: xNacc_nblx
    INTEGER, DIMENSION(:), ALLOCATABLE :: xNrej_nbl
    INTEGER, DIMENSION(:), ALLOCATABLE :: xNrej_nblx
    INTEGER, DIMENSION(:), ALLOCATABLE :: ierr_nbl
    INTEGER, DIMENSION(:), ALLOCATABLE :: ierr_nblx

    INTEGER :: jb
    INTEGER :: nbl  ! block length
    INTEGER :: nblx ! block length (current mechanism eXtracted)
    INTEGER :: ind_H2O

    ! declarations from previous update_physc_gp:
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: khet_Tr
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: khet_St
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: jx
    INTEGER :: ip, jj, status

    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ztend
    INTEGER :: idt

#ifdef MESSYTENDENCY
    REAL(DP), DIMENSION(kproma,nlev) :: diff_te
#endif

    !-------------------------------------------------------!
    !                   GRID  POINT                         !
    !-------------------------------------------------------!

    !mz_bs_20150702+
    IF (lmtskip) THEN 
       isw_mtskip = NINT(rmode_mtskip)
       ztmst = tsl_phys_mtskip
    ELSE 
       isw_mtskip = 0
       ztmst = time_step_len
    ENDIF
    IF (isw_mtskip == 2) RETURN
    !mz_bs_20150702-

    ! op_pj_20161220+
    IF (l_sync_H2O_gp) THEN
       idt = zidt_H2O_gp
          pxtm1(_RI_X_ZN_(1:kproma,:,idt)) = scvmr * (qm1(_RI_XYZ__(1:kproma,jrow,:)) / &
               (1.0_DP - qm1(_RI_XYZ__(1:kproma,jrow,:))))
          !
          ! CONSISTENT TENDENCIES (BEFORE CHEMISTRY)
          !
          !         d(H2O)   d             q                   dq/dt
          ! H2Ote = ------ = --  (scvmr * ----- ) = scvmr * -----------
          !           dt     dt           1 - q               (1-q)^2
          !
#ifndef MESSYTENDENCY
!          pxtte(_RI_X_ZN_(1:kproma,:,idt)) = scvmr * qte_3d(_RI_XYZ__(1:kproma,jrow,:)) / &
!               ( 1.0_DP - ( qm1(_RI_XYZ__(1:kproma,jrow,:)) + &
!               qte_3d(_RI_XYZ__(1:kproma,jrow,:)) *time_step_len ) )**2
          pxtte(_RI_X_ZN_(1:kproma,:,idt)) = scvmr * qte_3d(_RI_XYZ__(1:kproma,jrow,:)) / &
               ( ( 1.0_DP - ( qm1(_RI_XYZ__(1:kproma,jrow,:)) + &
               qte_3d(_RI_XYZ__(1:kproma,jrow,:)) *time_step_len ) ) * &
               ( 1.0_DP - qm1(_RI_XYZ__(1:kproma,jrow,:)) ) )
#else
          ! difference qte (converted) and pxtte
!          diff_te(:,:) =  scvmr * qte_3d(_RI_XYZ__(1:kproma,jrow,:)) / &
!               ( 1.0_DP - ( qm1(_RI_XYZ__(1:kproma,jrow,:)) + &
!               qte_3d(_RI_XYZ__(1:kproma,jrow,:)) *time_step_len ) )**2 - &
!               pxtte(_RI_X_ZN_(1:kproma,:,idt))
          diff_te(:,:) =  scvmr * qte_3d(_RI_XYZ__(1:kproma,jrow,:)) / &
               ( ( 1.0_DP - ( qm1(_RI_XYZ__(1:kproma,jrow,:)) + &
               qte_3d(_RI_XYZ__(1:kproma,jrow,:)) *time_step_len ) ) * &
               ( 1.0_DP - qm1(_RI_XYZ__(1:kproma,jrow,:)) ) ) - &
               pxtte(_RI_X_ZN_(1:kproma,:,idt))
          ! add tendency
!!$          CALL mtend_add_l(my_handle_corr, mtend_id_tracer, &
!!$               px=diff_te, idt=zidt_H2O_gp)   
          CALL mtend_add_l(my_handle_corr,zidt_H2O_gp, px=diff_te)   
#endif
!!$       ENDIF
    ENDIF
    ! op_pj_20161220-

    nbl = kproma*nlev

    IF (l_gp) THEN
      test_if_skip_gp: IF (.NOT.l_skipkpp_gp) THEN
         IF (l_khet) CALL mecca_khet_physc_gp
         ntrac = ntrac_gp
        ALLOCATE(temp_nbl(nbl))
        ALLOCATE(press_nbl(nbl))
        ALLOCATE(sphum_nbl(nbl))
        ALLOCATE(cair_nbl(nbl))
        ALLOCATE(xNacc_nbl(nbl))
        ALLOCATE(xNrej_nbl(nbl))
        ALLOCATE(ierr_nbl(nbl))
        ALLOCATE(zmrbc(nbl,ntrac))
        ALLOCATE(zmrac(nbl,ntrac))
        ALLOCATE(meccanum_nbl(nbl))
        ALLOCATE(l_meccanum_nbl(nbl))
        ALLOCATE(ztend(kproma,nlev))     ! op_pj_20150430
        IF (lfq_gp) ALLOCATE(zwvac(nbl)) ! op_pj_20160511

        ! for ECHAM5 (see bmil/messy_main_ppd_bi.inc):
        temp_nbl(:)  = kpromanlev_to_nbl(tm1_3d(_RI_XYZ__(1:kproma,jrow,:)) &
          + tte_3d(_RI_XYZ__(1:kproma,jrow,:)) * time_step_len)
        press_nbl(:) = kpromanlev_to_nbl(press_3d(_RI_XYZ__(1:kproma,jrow,:)))
        sphum_nbl(:) = kpromanlev_to_nbl(MAX(qm1_3d(_RI_XYZ__(1:kproma,jrow,:)) + &
          qte_3d(_RI_XYZ__(1:kproma,jrow,:)) * time_step_len, 0._dp)) ! specific humidity
        cair_nbl(:)  = (N_A/1.E6_dp) * press_nbl(:) &
          / (R_gas*temp_nbl(:)*(1.0_dp+vtmpc1*sphum_nbl(:)))

        IF (l_polymecca) THEN
          meccanum_nbl(:) = NINT(kpromanlev_to_nbl(&
               meccanum_gp%ptr(_RI_XYZ__(1:kproma,jrow,:))))
        ELSE
          meccanum_nbl(:) = 1 ! default is mechanism number 1
        ENDIF

        ! ------------------------------------------------------------------
        ! estimate tracer mixing ratio (mr) before chemistry (bc)
#ifndef MESSYTENDENCY
        DO jt=1,ntrac
#ifdef ALLATONCE
          zmrbc(:,jt) = kpromanlev_to_nbl(pxtm1(_RI_X_ZN_(1:kproma,:,jt)) + &
            pxtte(_RI_X_ZN_(1:kproma,:,jt)) * time_step_len)
#else
          jb = 0
          DO jk=1,nlev
            DO jp=1,kproma
              jb = jb+1
              zmrbc(jb,jt) = pxtm1(_RI_X_ZN_(jp,jk,jt)) + pxtte(_RI_X_ZN_(jp,jk,jt))*time_step_len
            ENDDO
          ENDDO
#endif
        ENDDO
#else
#ifdef ALLATONCE
        ALLOCATE(zmr_tmp_3d(_RI_X_ZN_(kproma,nlev,ntrac_gp)))
        zmr_tmp_3d(:,:,:) = 0.0_dp
        CALL mtend_get_start_l(mtend_id_tracer, v0t=zmr_tmp_3d)
        DO jt=1,ntrac_gp
! ub_ak_20190613 zmrbc(:,jt) = kpromanlev_to_nbl(zmr_tmp_3d(:,:,jt))
          zmrbc(:,jt) = kpromanlev_to_nbl(zmr_tmp_3d(_RI_X_ZN_(:,:,jt)))
        ENDDO
        DEALLOCATE(zmr_tmp_3d)
#else
        ALLOCATE(zmr_tmp_2d(kproma,nlev))
        DO jt=1,ntrac_gp
!          CALL mtend_get_start_l(mtend_id_tracer, v0=zmr_tmp_2d, idt=jt)
          CALL mtend_get_start_l(jt, v0=zmr_tmp_2d)
          zmrbc(:,jt) = kpromanlev_to_nbl(zmr_tmp_2d(:,:))
        ENDDO
        DEALLOCATE(zmr_tmp_2d)
#endif
#endif
        IF (lcheck_range) THEN
          CALL check_range_gp('before kpp:',zmrbc)
        ENDIF
        ! set negative values to zero:
        zmrbc(:,:) = MAX(zmrbc(:,:),0._dp)
        zmrac(:,:) = zmrbc(:,:)
        IF (lfq_gp) zwvac(:) = sphum_nbl(:) ! op_pj_20160511

        ! --------------------------------------------

        polymecca_loop: DO jpm=1,NMAXMECCA
          l_meccanum_nbl(:) = (meccanum_nbl(:)==jpm)
          nblx = COUNT(l_meccanum_nbl) ! nbl eXtracted for current mechanism
          IF (nblx == 0) CYCLE
          ALLOCATE(sphum_nblx(nblx))
          ALLOCATE(cair_nblx(nblx))
          ALLOCATE(xNacc_nblx(nblx))
          ALLOCATE(xNrej_nblx(nblx))
          ALLOCATE(ierr_nblx(nblx))
          ALLOCATE(conc(nblx, NSPEC(jpm)))
          ALLOCATE(zmrbc_nblx(nblx,ntrac))
          ALLOCATE(zmrac_nblx(nblx,ntrac))
          ! op_pj_20160511+
          IF (lfq_gp) THEN
             ALLOCATE(zwvac_nblx(nblx))
             ALLOCATE(zwvmrac_nblx(nblx))
          END IF
          ! op_pj_20160511-

          conc(:,:) = 0.0_dp

          DO jt=1,ntrac
            zmrbc_nblx(:,jt) = nbl_to_nblx(zmrbc(:,jt),l_meccanum_nbl,nblx)
          ENDDO
          zmrac_nblx = zmrbc_nblx
!!$       IF (lfq_gp) zwvac_nblx(:) = sphum_nblx(:) ! op_pj_20160511

          CALL fill_temp (jpm, status, nbl_to_nblx(temp_nbl,l_meccanum_nbl,nblx))
          IF (status /= 0) CALL error_bi('fill_temp array size', substr)

          CALL fill_press (jpm, status, nbl_to_nblx(press_nbl,l_meccanum_nbl,nblx))
          IF (status /= 0) CALL error_bi('fill_press array size', substr)

          cair_nblx = nbl_to_nblx(cair_nbl,l_meccanum_nbl,nblx)
          CALL fill_cair (jpm, status, cair_nblx)
          IF (status /= 0) CALL error_bi('fill_cair array size', substr)

          sphum_nblx = nbl_to_nblx(sphum_nbl,l_meccanum_nbl,nblx)
          IF (lfq_gp) zwvac_nblx(:) = sphum_nblx(:) ! op_pj_20160511

          if_khet: IF (l_khet) THEN ! define heterogeneous rate constants
            ALLOCATE(khet_Tr(nblx,IHT_MAX))
            IF ((l_force_khet).OR.(l_troposphere)) THEN ! troposphere
              DO jj=1, IHT_MAX
                khet_Tr(:,jj) = nbl_to_nblx( &
                  kpromanlev_to_nbl( &
                  khet_Tr_3d(xsm_cpl,jj)%ptr(_RI_XYZ__(1:kproma,jrow,:))), &
                  l_meccanum_nbl,nblx)
              ENDDO
            ELSE
              khet_Tr(:,:) = 0._dp
            ENDIF
            CALL fill_khet_Tr (jpm, status, khet_Tr)
            IF (status /= 0) CALL error_bi('fill_khet_Tr array size', substr)
            DEALLOCATE(khet_Tr)
            ALLOCATE(khet_St(nblx,IHS_MAX))
            IF ((l_force_khet).OR.(l_stratosphere)) THEN ! stratosphere
              DO jj=1, IHS_MAX
                khet_St(:,jj) = nbl_to_nblx( &
                  kpromanlev_to_nbl( &
                  khet_St_3d(jj)%ptr(_RI_XYZ__(1:kproma,jrow,:))), &
                  l_meccanum_nbl,nblx)
              ENDDO
            ELSE
              khet_St(:,:) = 0._dp
            ENDIF
            CALL fill_khet_St (jpm, status, khet_St)
            IF (status /= 0) CALL error_bi('fill_khet_St array size', substr)
            DEALLOCATE(khet_St)
          ENDIF if_khet

          ! define photolysis rate constants for kpp
          ! If photrat_gp(ip)%ptr is not associated, because
          !  - the tracer does not exist, or
          !  - the photolysis submodel does not provide it
          ! then jx = 0.
          ALLOCATE(jx(nblx,IP_MAX))
          jx(:,:) = 0.0_dp
          DO ip=1, IP_MAX
            IF (ASSOCIATED(photrat_gp(ip)%ptr)) THEN
              jx(:,ip) = nbl_to_nblx( &
                kpromanlev_to_nbl(photrat_gp(ip)%ptr(_RI_XYZ__(1:kproma,jrow,:))), &
                l_meccanum_nbl,nblx)
            ENDIF
          ENDDO
          CALL fill_jx (jpm, status, jx)
          IF (status /= 0) CALL error_bi('fill_jx array size', substr)
          DEALLOCATE(jx)

          CALL mr2c(jpm, zmrbc_nblx, cair_nblx, conc, ind_H2O) ! mr to conc
          IF (ind_H2O /= 0) THEN ! mz_ht_20150112
            IF (zidt_H2O_gp/=0) THEN ! H2O from tracer:
              conc(:,ind_H2O) = cair_nblx(:) * zmrbc_nblx(:,zidt_H2O_gp)
            ELSE ! H2O from humidity:
              conc(:,ind_H2O) = cair_nblx(:) * &
                scvmr * sphum_nblx(:) / (1._dp - sphum_nblx(:))
            ENDIF
          ENDIF ! mz_ht_20150112

          IF (l_aero) THEN ! doesn't work with polymecca (yet?)
            ! qqq change to: temp_nblx, press_nblx, cair_nblx, conc, zmrbc_nblx)
            CALL mecca_aero_update_physc(nblx, &
              temp_nbl, press_nbl, cair_nbl, conc, zmrbc_nblx)
          ENDIF
#ifdef MECCA_TAG
          IF (l_tag) & ! doesn't work with polymecca (yet?)
            CALL mecca_tag_preprocess(conc)
#endif

          ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
          CALL kpp_integrate(jpm, ztmst, conc, ierrf=ierr_nblx, &
            xNacc=xNacc_nblx, xNrej=xNrej_nblx, l_debug=l_kpp_debug, PE=p_pe)
          ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP

#ifdef MECCA_TAG
          IF (l_tag) & ! doesn't work with polymecca (yet?)
            CALL mecca_tag_postprocess(conc)
#endif
          IF ((.NOT. lcheck_range) .AND. (.NOT. l_fixed_step(jpm))) THEN
            conc(:,:) = MAX(conc(:,:), 0.0_dp)
          ENDIF

          IF (l_aero) CALL mecca_aero_save_pH(conc) ! doesn't work with polymecca (yet?)
          CALL c2mr(jpm, zmrac_nblx, conc, cair_nblx) ! conc to mr
          IF (ind_H2O /= 0) THEN ! op_pj_20160510
             ! Note: H2O tracer can be present for two reasons:
             !       (1) H2O tracer is defined by H2O submodel
             !       (2) H2O tracer is defined by MECCA, if MECCA_TAG is used
             !       In case (1), the conversion must be explicit here, since
             !       it is not contained in the c2mr include file.
             !       In case (2), the calculation is performed twice ...
             IF (zidt_H2O_gp/=0) & ! if H2O tracer exists, define zmrac:
                  zmrac_nblx(:,zidt_H2O_gp) = conc(:,ind_H2O) / cair_nblx(:)
             ! op_pj_20160511+
             ! Feedback should work without H2O tracer
             ! (no H2O submodel, no MECCA_TAG) and with H2O tracer
             ! (H2O submodel or MECCA_TAG).
             IF (lfq_gp) THEN
                ! H2O mixing ratio
                zwvmrac_nblx(:) = conc(:,ind_H2O) / cair_nblx(:)
                ! specific humidity: q = H2O/(scvmr + H2O)
                zwvac_nblx(:) = zwvmrac_nblx(:)/(scvmr + zwvmrac_nblx(:))
             END IF
             ! op_pj_20160511-
          ENDIF ! op_pj_20160510

          DO jt=1,ntrac
            CALL nblx_to_nbl(zmrac_nblx(:,jt), zmrac(:,jt), l_meccanum_nbl)
          ENDDO
          CALL nblx_to_nbl(xNacc_nblx, xNacc_nbl, l_meccanum_nbl)
          CALL nblx_to_nbl(xNrej_nblx, xNrej_nbl, l_meccanum_nbl)
          CALL nblx_to_nbl(ierr_nblx,  ierr_nbl,  l_meccanum_nbl)
          ! op_pj_20160511+
          IF (lfq_gp) CALL nblx_to_nbl(zwvac_nblx(:), zwvac(:), l_meccanum_nbl)
          ! op_pj_20160511-

          DEALLOCATE(sphum_nblx)
          DEALLOCATE(cair_nblx)
          DEALLOCATE(conc)
          DEALLOCATE(zmrbc_nblx)
          DEALLOCATE(zmrac_nblx)
          DEALLOCATE(xNacc_nblx)
          DEALLOCATE(xNrej_nblx)
          DEALLOCATE(ierr_nblx)
          ! op_pj_20160511+
          IF (lfq_gp) THEN
             DEALLOCATE(zwvac_nblx)
             DEALLOCATE(zwvmrac_nblx)
          END IF
          ! op_pj_20160511-
        ENDDO polymecca_loop

        ! --------------------------------------------

        jb=0
        DO jk=1,nlev
          DO jp=1,kproma
            jb = jb +1
            out_nsteps_gp(_RI_XYZ__(jp,jrow,jk)) = REAL(xNacc_nbl(jb)+xNrej_nbl(jb), dp)
            IF (ierr_nbl(jb) < -5) THEN
              IF (l_kpp_debug) &
                CALL write_error(p_pe,jp,jk,jrow,philon_2d(jp,jrow) &
                , philat_2d(jp,jrow))
#if defined(ECHAM5) || defined(CESM1)
              WRITE(*,*) substr,': ',IERR_NAMES(ierr_nbl(jb)), &
                ' PE=',p_pe,' JP=',jp,' JK=',jk,' JROW=',jrow , &
                ' LON=',philon_2d(jp,jrow),' (',ilon(jp,jrow),')', &
                ' LAT=',philat_2d(jp,jrow),' (',ilat(jp,jrow),')'
#endif
#ifdef COSMO
              WRITE(*,*) substr,': ',IERR_NAMES(ierr_nbl(jb)), &
                ' PE=',p_pe,' JP=',jp,' JK=',jk,' JROW=',jrow , &
                ' LON=',philon_2d(jp,jrow), &
                ' LAT=',philat_2d(jp,jrow)
#endif
            ELSE IF (ierr_nbl(jb) < 0) THEN
              CALL error_bi(IERR_NAMES(ierr_nbl(jb)),substr)
            ENDIF
          ENDDO
        ENDDO

        IF (lcheck_range) THEN
          CALL check_range_gp('after kpp:',zmrac)
        ENDIF
        IF (l_aero) CALL mecca_aero_diag_si(zmrbc,zmrac, nbl)
#ifndef MESSYTENDENCY
        DO jt=1,ntrac
          ztend(:,:) = 0.0_dp ! op_pj_20150430
          jb = 0
          DO jk=1,nlev
            DO jp=1,kproma
              jb = jb+1
              IF (ierr_nbl(jb) < 0) CYCLE ! no tendency for problematic boxes
              ! op_pj_20150707+
              ztend(jp,jk) = (zmrac(jb,jt) - zmrbc(jb,jt)) / ztmst
              pxtte(_RI_X_ZN_(jp,jk,jt)) = pxtte(_RI_X_ZN_(jp,jk,jt)) + ztend(jp,jk)
              ! op_pj_20150707-
            ENDDO
          ENDDO
          !mz_bs_20150702+
          IF ((lmtskip) .AND. (isw_mtskip == 1)) THEN
            IF ( l_trac_mtskip(jt) ) THEN
              mtskiptend(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
                mtskiptend(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
                ztend(1:kproma,:)
            ENDIF
          ENDIF
          !mz_bs_20150702-
        ENDDO
        ! op_pj_20160511+
        IF (lfq_gp) THEN
           ztend(:,:) = 0.0_dp
           jb = 0
           DO jk=1,nlev
              DO jp=1,kproma
                 jb = jb+1
                 IF (ierr_nbl(jb) < 0) CYCLE ! no tendency for problematic boxes
                 ztend(jp,jk) = ( zwvac(jb) - sphum_nbl(jb) ) / ztmst
                 qte_3d(_RI_XYZ__(jp,jrow,jk)) = qte_3d(_RI_XYZ__(jp,jrow,jk)) + ztend(jp,jk)
              END DO
           END DO
           IF ((lmtskip) .AND. (isw_mtskip == 1)) THEN
              mtskiptend_q(_RI_XYZ__(1:kproma,jrow,:)) = &
                   mtskiptend_q(_RI_XYZ__(1:kproma,jrow,:)) + ztend(1:kproma,:)
           END IF
        ENDIF
        ! op_pj_20160511-
#else
#ifdef ALLATONCE
        ALLOCATE(zmr_tmp_3d(_RI_X_ZN_(kproma,nlev,ntrac_gp)))
        zmr_tmp_3d(:,:,:) = 0.0_dp
        DO jt=1, ntrac_gp
          jb = 0
          DO jk=1,nlev
            DO jp=1,kproma
              jb = jb+1
              IF (ierr_nbl(jb) < 0) CYCLE ! no tendency for problematic boxes
              zmr_tmp_3d(_RI_X_ZN_(jp,jk,jt)) = &
                (zmrac(jb,jt) - zmrbc(jb,jt)) / ztmst
            ENDDO
          ENDDO
        ENDDO
        CALL mtend_add_l(my_handle, mtend_id_tracer, pxt=zmr_tmp_3d)
        !mz_bs_20150702+
        IF ((lmtskip) .AND. (isw_mtskip == 1)) THEN
          DO jt=1, ntrac
            IF ( l_trac_mtskip(jt) ) THEN
              mtskiptend(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
                mtskiptend(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
                zmr_tmp_3d(_RI_X_ZN_(1:kproma,1:nlev,jt))
            ENDIF
          ENDDO
        ENDIF
        !mz_bs_20150702-
        DEALLOCATE(zmr_tmp_3d)
#else
        ALLOCATE(zmr_tmp_2d(kproma,nlev))
        DO jt=1, ntrac_gp
          zmr_tmp_2d(:,:) = 0.0_dp
          jb = 0
          DO jk=1,nlev
            DO jp=1,kproma
              jb = jb+1
              IF (ierr_nbl(jb) < 0) CYCLE ! no tendency for problematic boxes
              zmr_tmp_2d(jp,jk) = (zmrac(jb,jt) - zmrbc(jb,jt)) / ztmst
            ENDDO
          ENDDO
!          CALL mtend_add_l(my_handle, mtend_id_tracer, px=zmr_tmp_2d, idt=jt)
          CALL mtend_add_l(my_handle, jt, px=zmr_tmp_2d)
          !mz_bs_20150702+
          IF ((lmtskip) .AND. (isw_mtskip == 1)) THEN
            IF ( l_trac_mtskip(jt) ) THEN
              mtskiptend(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
                mtskiptend(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + zmr_tmp_2d(1:kproma,:)
            ENDIF
          ENDIF
          !mz_bs_20150702-
        ENDDO
        DEALLOCATE(zmr_tmp_2d)
#endif
! ALLATONCE
        ! op_pj_20160511+
        ! tendency q by using difference in specific humidity
        IF (lfq_gp) THEN
           ALLOCATE(zmr_tmp_2d(kproma,nlev))
           zmr_tmp_2d(:,:) = 0.0_dp
           jb = 0
           DO jk=1,nlev
              DO jp=1,kproma
                 jb = jb+1
                 IF (ierr_nbl(jb) < 0) CYCLE ! no tendency for problematic boxes
                 zmr_tmp_2d(jp,jk) = ( zwvac(jb) - sphum_nbl(jb) ) / ztmst
              ENDDO
           ENDDO
           CALL mtend_add_l(my_handle, mtend_id_q, px=zmr_tmp_2d)
           IF ((lmtskip) .AND. (isw_mtskip == 1)) THEN
              mtskiptend_q(_RI_XYZ__(1:kproma,jrow,:)) = &
                      mtskiptend_q(_RI_XYZ__(1:kproma,jrow,:)) + zmr_tmp_2d(1:kproma,:)
           ENDIF
           DEALLOCATE(zmr_tmp_2d)
        END IF
        ! op_pj_20160511-
#endif
! MESSYTENDENCY
        DEALLOCATE(cair_nbl)
        DEALLOCATE(ierr_nbl)
        DEALLOCATE(l_meccanum_nbl)
        DEALLOCATE(meccanum_nbl)
        DEALLOCATE(press_nbl)
        DEALLOCATE(sphum_nbl)
        DEALLOCATE(temp_nbl)
        DEALLOCATE(xNacc_nbl)
        DEALLOCATE(xNrej_nbl)
        DEALLOCATE(zmrac)
        DEALLOCATE(zmrbc)
        DEALLOCATE(ztend)             ! op_pj_20150430
        IF (lfq_gp) DEALLOCATE(zwvac) ! op_pj_20160511-

        IF (l_aero) CALL mecca_aero_dealloc

      ENDIF test_if_skip_gp
    ENDIF

  CONTAINS

    !-------------------------------------------------------------------------

    SUBROUTINE check_range_gp(infostring,mixrat)

      CHARACTER(LEN=*), INTENT(IN) :: infostring
      REAL(DP),         INTENT(IN) :: mixrat(:,:) ! tracer mixing ratio
      INTEGER :: jt, jb, jk, jp

      INTRINSIC :: SIZE, TRIM

      jb = 0

      level_loop: DO jk=1,nlev
        kproma_loop: DO jp=1,kproma
          jb = jb+1

          tracer_loop: DO jt = 1, ntrac_gp
            wrong_mr: IF (((mixrat(jb,jt) < negerr) .OR. &
              (mixrat(jb, jt) > 1.0_dp)) &
              .AND.(TRIM(ti_gp(jt)%tp%ident%submodel) == 'mecca')) THEN
              WRITE(*,'(2A,I3,A,I3,A,I3,A,1PG12.3E3,2A)') infostring, &
#if defined(ECHAM5)
                ' ilon =',ilon(jp,jrow), &
                ', ilat =',ilat(jp,jrow), &
#else
                ' ilon =',jp, ', ilat =',jrow, &
#endif
                ', lev =', jk, &
                ', y =', mixrat(jb, jt),' mol/mol for ', &
                TRIM(ti_gp(jt)%tp%ident%fullname)
            ENDIF wrong_mr
          ENDDO tracer_loop

        ENDDO kproma_loop
      ENDDO level_loop

    END SUBROUTINE check_range_gp

    !-------------------------------------------------------------------------

  END SUBROUTINE mecca_physc

  !***************************************************************************

  SUBROUTINE mecca_global_end
!!#D attila +
#if defined(ECHAM5)
    USE messy_main_timer,            ONLY: time_step_len
    USE messy_main_data_bi,          ONLY: tm1_3d, tte_3d,  &
                                           qm1_3d, qte_3d
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_tracer_mem_bi,    ONLY: ntrac_lg, ti_lg, NCELL, &
      pxtte_a => qxtte_a, &
      pxtm1_a => qxtm1_a
    ! MESSy
    USE messy_attila_tools_e5,     ONLY: gp2lg_e5, lg2gp_e5, LG2GP_SUM 
    USE messy_main_constants_mem,  ONLY: N_A, R_gas, M_air, M_H2O, OneDay
    ! MECCA, MECCA_AERO
    USE messy_mecca_kpp,           ONLY: NSPEC, kpp_integrate
! op_pj_20160511+
#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,   ONLY: mtend_add_g, mtend_id_q
#endif
! op_pj_20160511-

    IMPLICIT NONE

    INTRINSIC MAX, MOD, INT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_global_end'
    REAL(DP), PARAMETER :: scvmr = M_air/M_H2O ! op_pj_20160510
    INTEGER  :: jt ! tracer loop index
    INTEGER  :: ntrac
    ! tracer mixing ratio before chemistry
    REAL(DP), DIMENSION(NCELL,ntrac_lg) :: zmrbc
    ! tracer mixing ratio after chemistry
    REAL(DP), DIMENSION(NCELL,ntrac_lg) :: zmrac
    REAL(DP), DIMENSION(NCELL)          :: zwvac, zwvbc ! op_pj_20160511

    ! DEFINE 1D temp, press, cair here, as they are needed in
    !        mecca_aero_update_physc too.
    !REAL(DP), DIMENSION(NCELL) :: press
    REAL(DP), DIMENSION(NCELL) :: cair
    REAL(DP), DIMENSION(NCELL,NSPEC) :: conc

    ! info from KPP integration
    INTEGER, DIMENSION(NCELL) :: xNacc
    INTEGER, DIMENSION(NCELL) :: xNrej
    INTEGER, DIMENSION(NCELL) :: ierr

    INTEGER          :: jb

    !-------------------------------------------------------!
    !                        LAGRANGE                       !
    !-------------------------------------------------------!
    IF (l_lg) THEN

      ! op_pj_20161220+
      IF (l_sync_H2O_lg) THEN
         CALL error_bi(&
              'synchronisation of q and H2O not yet implemented for LG' &
              , substr)
      ENDIF
      ! op_pj_20161220-

      ! we transform here so we have less calculation pro time step
      ! NB it is used as well in mecca_khet_physc_lg
      !temperature
      tmp_temp_gp(:,:,:)=tm1_3d(:,:,:)+     &
        tte_3d(:,:,:)* time_step_len
      !specific humidity
      tmp_sphum_gp(:,:,:)=MAX(qm1_3d(:,:,:)+     &
        qte_3d(:,:,:)*time_step_len, 0._dp)
      CALL gp2lg_e5(tmp_temp_gp, tmp_temp_lg)
      CALL gp2lg_e5(tmp_sphum_gp, tmp_sphum_lg)

      IF (l_khet) CALL mecca_khet_physc_lg(tmp_temp_lg,PRESS_PARCEL)
      ntrac = ntrac_lg

      test_if_skip_lg: IF (.NOT.l_skipkpp_lg) THEN

        CALL update_physc_lg

        ! estimate tracer mixing ratio (mr) before chemistry (bc)
        zmrbc(:,:) = pxtm1_a(:,:) + pxtte_a(:,:) * time_step_len

        IF (lcheck_range) CALL check_range_lg('before kpp:',zmrbc)

        ! set negative values to zero
        zmrbc(:,:) = MAX(zmrbc(:,:),0._dp)
        zmrac(:,:) = zmrbc(:,:)
        conc(:,:) = 0._dp ! default
        CALL mr2c(cair) ! convert mr to conc
        !IF (l_aero) &
        ! CALL mecca_aero_update_physc( NCELL, temp, press, cair &
        !                             , conc, zmrbc)

        CALL kpp_integrate(time_step_len, conc, ierrf=ierr &
          , xNacc=xNacc, xNrej=xNrej)

        out_nsteps_lg(:) = REAL(xNacc(:)+xNrej(:), dp)

        ! ERROR CONTROL
        DO jb=1, NCELL
          IF (ierr(jb) < -5) THEN
            WRITE(*,*) substr,': ',IERR_NAMES(ierr(jb)), &
              ' PE=',p_pe,' CELL=',jb
          ELSE IF (ierr(jb) < 0) THEN
            CALL error_bi(IERR_NAMES(ierr(jb)),substr)
          ENDIF
        ENDDO

        !IF (l_aero) CALL mecca_aero_save_pH(conc)

        CALL c2mr(1._dp/cair) ! convert conc to mr

        IF (lcheck_range) CALL check_range_lg('after kpp: ',zmrac)

        ! add chemical mixing ratio (mr) tendency (te) to
        ! total tendency
        DO jb=1, NCELL
          IF (ierr(jb) < 0) CYCLE ! no tendency for problematic boxes
          DO jt=1,ntrac_lg
            pxtte_a(jb,jt) = pxtte_a(jb,jt) &
              + (zmrac(jb,jt) - zmrbc(jb,jt)) / time_step_len
          ENDDO
        ENDDO

        ! op_pj_20160511+
        IF (lfq_lg) THEN
           qtend_lg = ( zwvac(:)/(scvmr+zwvac(:)) - &
                zwvbc(:)/(scvmr+zwvbc(:)) ) &
                / time_step_len
           ! qqq: it needs to be checked, if this transformation is
           !      appropriate
           CALL lg2gp_e5(qtend_lg, qtend_gp, LG2GP_SUM, fill_value = 0.0_dp)
#ifndef MESSYTENDENCY
           qte_3d(:,:,:) = qte_3d(:,:,:) + qtend_gp(:,:,:)
#else
           CALL mtend_add_g(my_handle, mtend_id_q, px=qtend_gp)
#endif
        ENDIF
        ! op_pj_20160511-

      ENDIF test_if_skip_lg
    ENDIF

  CONTAINS
    !-------------------------------------------------------------------------

    SUBROUTINE check_range_lg(infostring,mixrat)

      CHARACTER(LEN=*), INTENT(IN) :: infostring
      REAL(DP),         INTENT(IN) :: mixrat(:,:) ! tracer mixing ratio
      INTEGER :: jt, jp

      INTRINSIC :: SIZE, TRIM

      cell_loop: DO jp=1, NCELL

        tracer_loop: DO jt = 1, ntrac_lg
          wrong_mr: IF ( ((mixrat(jp,jt) < negerr) &
            .OR. (mixrat(jp, jt) > 1.0_dp) ) &
            .AND.(TRIM(ti_lg(jt)%tp%ident%submodel) == 'mecca')) THEN
            WRITE(*,'(2A,I8,A,1PG12.3E3,2A)') infostring, &
              ' CELL =',jp, &
              ', y =', mixrat(jp, jt),' mol/mol for ', &
              TRIM(ti_lg(jt)%tp%ident%fullname)
          ENDIF wrong_mr
        ENDDO tracer_loop

      ENDDO cell_loop

    END SUBROUTINE check_range_lg
    !-------------------------------------------------------------------------

    SUBROUTINE update_physc_lg

      USE messy_mecca_kpp, ONLY: &
        fill_temp, fill_cair, fill_press, fill_jx, fill_khet_Tr, fill_khet_St
      INTRINSIC MAX, ASSOCIATED

      CHARACTER(LEN=*), PARAMETER :: substr = 'update_physc_lg'
      REAL(DP), DIMENSION(NCELL,IHT_MAX) :: khet_Tr
      REAL(DP), DIMENSION(NCELL,IHS_MAX) :: khet_St
      REAL(DP), DIMENSION(NCELL,IP_MAX)  :: jx
      INTEGER :: ip, jj, status

      CALL fill_temp(status, tmp_temp_lg)
      IF (status /= 0) &
        CALL error_bi('array size of temp (LG) is incompatible to fill' &
        , substr)

      cair(:)    = (N_A/1.E6_dp) * press_parcel(:) / &
        (R_gas*tmp_temp_lg(:)*(1.0_dp+vtmpc1*tmp_sphum_lg(:)))
      ! previous definition was without considering humidity:
      ! (N_A/1.E6_dp) * press / (R_gas*temp)
      CALL fill_cair(status, cair)
      IF (status /= 0) &
        CALL error_bi('array size of cair (LG) is incompatible to fill' &
        , substr)

      CALL fill_press(status, press_parcel)
      IF (status /= 0) &
        CALL error_bi('array size of press (LG) is incompatible to fill' &
        , substr)

      ! NB we repeat the calculation made in GP....
      ! so we have independent module!
      IF (l_khet) THEN ! define heterogeneous rate constants for kpp
        IF ((l_force_khet).OR.(l_troposphere)) THEN ! troposphere
          DO jj=1, IHT_MAX
            khet_Tr(:,jj) = khet_Tr_1d(xsm_cpl,jj)%ptr(:)
          ENDDO
        ELSE
          khet_Tr(:,:) = 0._dp
        ENDIF
        IF ((l_force_khet).OR.(l_stratosphere)) THEN ! stratosphere
          DO jj=1, IHS_MAX
            khet_St(:,jj) = khet_St_1d(jj)%ptr(:)
          ENDDO
        ELSE
          khet_St(:,:) = 0._dp
        ENDIF
        CALL fill_khet_Tr(status, khet_Tr)
        IF (status /= 0) &
          CALL error_bi('array size of khet_Tr (LG) is incompatible to fill' &
          , substr)
        CALL fill_khet_St(status, khet_St)
        IF (status /= 0) &
          CALL error_bi('array size of khet_St (LG) is incompatible to fill' &
          , substr)
      ENDIF

      ! define photolysis rate constants for kpp
      ! j_XYZ_3d pointers are not associated when the corresponding tracer
      ! does not exist. In this case, j_XYZ will not be used. It keeps
      ! its initial value of 0.
      DO ip=1, IP_MAX
        IF (ASSOCIATED(photrat_lg(ip)%ptr)) jx(:,ip) = photrat_lg(ip)%ptr(:)
      ENDDO
      CALL fill_jx(status, jx)
      IF (status /= 0) &
        CALL error_bi('array size of jx (LG) is incompatible to fill' &
        , substr)

    END SUBROUTINE update_physc_lg
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    SUBROUTINE mr2c(c_air)
      ! convert mixing ratio [mol/mol] to concentration [mcl/cc]

      USE messy_mecca_kpp  ! without ONLY to get all ind_*
                           ! if somebody enlarges the mechanism

      REAL(DP), INTENT(IN) :: c_air(:)

      ! special case: H2O
      IF (ind_H2O /= 0) THEN ! mz_ht_20150112
        IF (zidt_H2O_lg/=0) THEN
          ! The tracer H2O is defined. It is assumed here that it contains
          ! the correct molar mixing ratio for H2O.
          conc(:,ind_H2O) = c_air(:) * zmrbc(:,zidt_H2O_lg)
          ! ... this can actually at the time being never be reached ...
          IF (lfq_lg) zwvbc(:) = zmrbc(:,zidt_H2O_lg) ! op_pj_20160511
        ELSE
          ! The tracer H2O is not defined. Obtain H2O from humidity.
          !qqq C(ind_H2O) = cair * sphum * scvmr
          conc(:,ind_H2O) = c_air(:) * scvmr * tmp_sphum_lg(:) &
            / (1._dp - tmp_sphum_lg(:))
          ! op_pj_20160511+
          IF (lfq_lg) zwvbc(:) = scvmr * tmp_sphum_lg(:) &
               / (1._dp - tmp_sphum_lg(:))
          ! op_pj_20160511-
        ENDIF
      ENDIF ! mz_ht_20150112

      ! the following INCLUDE statement contains all conversions
      ! for all species that are used in the kpp chemistry scheme
      INCLUDE 'messy_mecca_mr2c_si.inc'

    END SUBROUTINE mr2c

    !-------------------------------------------------------------------------

    SUBROUTINE c2mr(c_air)
      ! convert concentration [mcl/cc] to mixing ratio [mol/mol]

      USE messy_mecca_kpp  ! without ONLY to get all ind_*
                           ! if somebody enlarges the mechanism

      REAL(DP), INTENT(IN) :: c_air(:)
      REAL(DP), DIMENSION(NCELL) :: riac ! 1/c(air)

      riac(:) = 1._dp/c_air(:)
      ! the following INCLUDE statement contains all conversions
      ! for all species that are used in the kpp chemistry scheme
      INCLUDE 'messy_mecca_c2mr_si.inc'

      ! If the tracer H2O is defined, its chemical tendency must
      ! eventually be put into xtte (via zmrac)
      ! op_pj_20150112+
      ! IF (zidt_H2O_lg/=0) zmrac(:,zidt_H2O_lg) = riac(:) * conc(:,ind_H2O)
      IF ((zidt_H2O_lg/=0) .AND. (ind_H2O /=0)) &
        zmrac(:,zidt_H2O_lg) = riac(:) * conc(:,ind_H2O)
      ! op_pj_20150112-
      ! op_pj_20160511+
      IF ( (lfq_lg) .AND. (ind_H2O /=0) ) THEN
         zwvac(:) = riac(:) * conc(:,ind_H2O)
      END IF
      ! op_pj_20160511-

    END SUBROUTINE c2mr

    !-------------------------------------------------------------------------
#endif
!!#D attila -
  END SUBROUTINE mecca_global_end

  !***************************************************************************

  SUBROUTINE mecca_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    IF (l_lg) THEN
      IF (ASSOCIATED(tmp_temp_gp))  DEALLOCATE(tmp_temp_gp)
      IF (ASSOCIATED(tmp_sphum_gp)) DEALLOCATE(tmp_sphum_gp)
      IF (ASSOCIATED(tmp_temp_lg))  DEALLOCATE(tmp_temp_lg)
      IF (ASSOCIATED(tmp_sphum_lg)) DEALLOCATE(tmp_sphum_lg)
      ! op_pj_20160511+
      IF (ASSOCIATED(qtend_lg)) DEALLOCATE(qtend_lg)
      IF (ASSOCIATED(qtend_gp)) DEALLOCATE(qtend_gp)
      ! op_pj_20160511-
    ENDIF

    IF (l_khet) CALL mecca_khet_free_memory

    !mz_bs_20150702+
    IF (ASSOCIATED(l_trac_mtskip))  DEALLOCATE(l_trac_mtskip)
    IF (ASSOCIATED(mtskiptend))     DEALLOCATE(mtskiptend)
    !mz_bs_20150702-

  END SUBROUTINE mecca_free_memory

  !***************************************************************************

  SUBROUTINE mecca_read_nml_cpl(status, iou)

    ! MECCA MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    ! read namelist for 'coupling' to channel containing reaction rates
    ! Author: Astrid Kerkweg, MPICH, Sep 2004

    ! MESSy
    USE messy_main_tools,          ONLY: read_nml_open, read_nml_check &
      , read_nml_close
    USE messy_main_tracer_mem_bi,  ONLY: NGCELL

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    CHARACTER(LEN=3)            :: istr

    INTRINSIC TRIM

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    IF ((l_lg) .AND. (NGCELL > 0)) THEN
      l_lg = .TRUE.
    ELSE
      IF (l_lg) THEN
        CALL info_bi('l_lg = T in namelist')
        CALL info_bi('However no Lagrangian scheme activated (attila).')
        CALL info_bi(' ... setting l_lg = F')
      ENDIF
      l_lg = .FALSE.
    ENDIF

    IF (l_skipkpp_gp) THEN
      CALL info_bi('******************************************************')
      CALL info_bi(' KPP CHEMISTRY FOR GRID POINT INTEGRATION SWITCHED OFF')
      CALL info_bi('******************************************************')
    ENDIF

    IF (l_skipkpp_lg) THEN
      CALL info_bi('******************************************************')
      CALL info_bi(' KPP CHEMISTRY FOR LAGRANGE INTEGRATION SWITCHED OFF  ')
      CALL info_bi('******************************************************')
    ENDIF

    WRITE(istr,'(i3)') i_pa_amode
    CALL info_bi('AEROSOL SUBMODEL FOR PSEUDO AEROSOL TRACERS (if present): '&
      &//TRIM(c_pa_asm)//'; MODE='//istr)

  END SUBROUTINE mecca_read_nml_cpl

  !*****************************************************************************
!!$  SUBROUTINE error_output(istatus, rstatus, ierr, jp, jk, radius, zmrbc, sphum)
!!$
!!$    USE messy_main_data_bi,   ONLY: philon_2d, philat_2d &
!!$                                  , jrow, current_time_step
!!$    USE messy_main_tools_bi,  ONLY: find_next_free_unit
!!$
!!$    USE messy_main_constants_mem, ONLY:M_air,M_H2O
!!$
!!$    USE messy_mecca_kpp_monitor,  ONLY: SPC_NAMES
!!$    USE messy_mecca_kpp_parameters
!!$    IMPLICIT NONE
!!$
!!$    INTEGER,  INTENT(IN) :: istatus(20), ierr ! info from KPP integration
!!$    REAL(DP), INTENT(IN) :: rstatus(20)       ! info from KPP integration
!!$
!!$    INTEGER,  INTENT(IN) :: jp, jk
!!$    REAL(DP), INTENT(IN) :: radius(APN)
!!$    REAL(DP), INTENT(IN) :: zmrbc(:)
!!$    REAL(DP), INTENT(IN) :: sphum
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=250)  :: filename
!!$    INTEGER              :: iou               ! output unit
!!$    INTEGER              :: status            ! status
!!$!    INTEGER              :: iostat            ! status
!!$    INTEGER              :: i,j
!!$    CHARACTER(len=10)    :: strlon, strlat, strlev, strtime
!!$
!!$    i=philon_2d(jp,jrow)
!!$    j=90+philat_2d(jp,jrow)
!!$
!!$    strlon  = str(i)
!!$    strlat  = str(j)
!!$    strlev  = str(jk)
!!$    strtime = str(current_time_step)
!!$
!!$
!!$    write(filename,*) 'mecca_'//TRIM(strtime)//'_'//TRIM(strlon)//'_'//TRIM(strlat)//'_L'//TRIM(strlev)//'.txt'
!!$
  !    write(*,*) 'A'//TRIM(filename)
  !    write(*,*) 'A'//ADJUSTL(filename)
  !    write(*,*) 'A'//ADJUSTL(TRIM(filename))
!!$
!!$    iou = find_next_free_unit(100,200)
!!$
!!$    OPEN(IOU, iostat=status, FILE =TRIM(ADJUSTL(filename))&
!!$         , STATUS='NEW', ACTION= 'WRITE')
!!$
!!$    write(IOU,*) 'KPP ERRORSTATUS ierr:'
!!$    write(IOU,*) ierr
!!$    write(IOU,*)
!!$    write(IOU,*) 'KPP ERRORSTATUS ISTATUS:'
!!$    DO i=1,20
!!$       write(IOU,*) istatus(i)
!!$    ENDDO
!!$    write(IOU,*)
!!$    write(IOU,*) 'KPP ERRORSTATUS RSTATUS:'
!!$    DO i=1,20
!!$       write(IOU,*) rstatus(i)
!!$    ENDDO
!!$
!!$    write(IOU,*) 'BOX: CURRENT-TIME-STEP, LON, LAT, LEV'
!!$    write(IOU,*) CURRENT_TIME_STEP
!!$    write(IOU,*) philon_2d(jp,jrow)
!!$    write(IOU,*) philat_2d(jp,jrow)
!!$    write(IOU,*) jk
!!$
!!$    write(IOU,*) ' PHYSICS: temperature, pressure, cair'
!!$    write(IOU,*) temp
!!$    write(IOU,*) press
!!$    write(IOU,*) cair
!!$    write(IOU,*)
!!$
!!$    write(IOU,*) ' PHOTOLYSIS RATES (index ip, see kpp_global)'
!!$    DO i=1,59
!!$       write(IOU,*) jx(i)
!!$    ENDDO
!!$
!!$
!!$    write(IOU,*) 'SELECTED MECHANISM'
!!$    write(IOU,*) WANTED
!!$    write(IOU,*)
!!$
!!$    write(IOU,*) 'NUMBER OF SPECIES:'
!!$    write(IOU,*) NSPEC
!!$    write(IOU,*)
!!$    write(IOU,*)  'NAMES OF SPECIES:'
!!$    DO i=1,NSPEC
!!$       write(IOU,*)  SPC_NAMES(i)
!!$    ENDDO
!!$    write(IOU,*)
!!$
!!$      !*********************************************************************
!!$      ! CALCULATE concentrations
!!$      !*********************************************************************
!!$
!!$      IF (zidt_H2O/=0) THEN
!!$        ! The tracer H2O is defined. It is assumed here that it contains
!!$        ! the correct molar mixing ratio for H2O.
!!$        conc(ind_H2O) = cair * zmrbc(zidt_H2O)
!!$      ELSE
!!$        ! The tracer H2O is not defined. Obtain H2O from humidity.
!!$        !qqq conc(ind_H2O) = cair * sphum * scvmr
!!$        conc(ind_H2O) = cair * scvmr * sphum / (1._dp - sphum)
!!$      ENDIF
!!$
!!$      ! the following INCLUDE statement contains all conversions
!!$      ! for all species that are used in the kpp chemistry scheme
!!$      INCLUDE 'messy_mecca_mr2c_si.inc'
!!$      !*********************************************************************
!!$      !*********************************************************************
!!$
!!$    write(IOU,*) 'concentrations (molec/cm^3(air) before KPP'
!!$    DO i=1,NSPEC
!!$       write(IOU,*) conc(i)
!!$    ENDDO
!!$    write(IOU,*)
!!$
!!$    write(IOU,*) 'AEROSOL DATA:'
!!$    write(IOU,*) 'RADII'
!!$    DO i=1,APN
!!$       write(IOU,*) radius(i)
!!$    ENDDO
!!$    write(IOU,*) 'lwc'
!!$    DO i=1,APN
!!$       write(IOU,*) lwc(i)
!!$    ENDDO
!!$    write(IOU,*) 'xaer'
!!$    DO i=1,APN
!!$       write(IOU,*) xaer(i)
!!$    ENDDO
!!$    write(IOU,*) 'lwc'
!!$    DO i=1,APN
!!$       write(IOU,*) cvfac(i)
!!$    ENDDO
!!$    write(IOU,*)
!!$
!!$    CLOSE(IOU)
!!$
!!$  END SUBROUTINE error_output
!!$!****************************************************************************
  SUBROUTINE write_error(PE,jp,jk,jrow,lon,lat)

    USE messy_main_tools,      ONLY: find_next_free_unit

    INTEGER,  INTENT(IN) :: PE, jp,jk,jrow
    REAL(DP), INTENT(IN) :: lon,lat

    INTEGER, SAVE :: NUM =0
    INTEGER       :: iou

    CHARACTER(LEN=250)  :: filename
    CHARACTER(LEN=1000) :: strnum
    CHARACTER(LEN=1000) :: strPE

    NUM = NUM + 1

    strnum=str(NUM)
    strPE=str(PE)

    WRITE(filename,*) 'error_PE'//TRIM(STRPE)//'_'//TRIM(STRNUM)//'.txt'

    iou = find_next_free_unit(100,200)
    OPEN(IOU,FILE =TRIM(ADJUSTL(filename))&
      ,STATUS='NEW',ACTION= 'WRITE')

    WRITE(IOU,*) ' PE=',pe,' JP=',jp,' JK=',jk,' JROW=',jrow , &
      ' LON=',lon, ' LAT=',lat
    WRITE(IOU,*)
    CLOSE(IOU)

  END SUBROUTINE write_error

  !***************************************************************************

  !mz_bs_20150702+
  SUBROUTINE mecca_mtskip_init_coupling(iflag)

    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer,           ONLY: set_tracer, get_tracer, ON, I_MTSKIP
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iflag
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_mtskip_init_coupling'
    INTEGER          :: status
    INTEGER          :: jt, counter

    SELECT CASE(iflag)
    CASE(1)

      CALL start_message_bi(modstr, 'COUPLING TO MTSKIP (1)',substr) 

      ! set pointers to mtskip channel objects
      ! actual mode of mtskip (0,1,2)
      CALL get_channel_object(status,'mtskip', 'rmode', p0=rmode_mtskip)
      CALL channel_halt(substr, status)
      ! time step length of mtskip
      CALL get_channel_object(status,'mtskip', 'tslphys', p0=tsl_phys_mtskip)
      CALL channel_halt(substr, status)

      ! allocate additional memory
      ALLOCATE(l_trac_mtskip(ntrac_gp))
      l_trac_mtskip(1:ntrac_gp) = .FALSE.

      ! Note: %ptr(:,:,:) => NULL() does hardly consume any memory
      ALLOCATE(mtskiptend(ntrac_gp))
      DO jt=1,ntrac_gp
        IF ( TRIM(ti_gp(jt)%tp%ident%submodel) .EQ. 'mecca' ) THEN
          CALL set_tracer(status, GPTRSTR, jt, I_MTSKIP, ON)
          l_trac_mtskip(jt) = .TRUE.
        ENDIF
      ENDDO

      ! SPECIAL CASE H2O
      CALL get_tracer(status, GPTRSTR, 'H2O',idx=jt)
      IF ( status == 0 ) THEN
        CALL set_tracer(status, GPTRSTR, jt, I_MTSKIP, ON)
        l_trac_mtskip(jt) = .TRUE.
      ENDIF

      CALL end_message_bi(modstr, 'COUPLING TO MTSKIP (1)',substr) 

    CASE(2)

      CALL start_message_bi(modstr, 'COUPLING TO MTSKIP (2)',substr) 

      counter = 0
      DO jt = 1, ntrac_gp          
        IF ( l_trac_mtskip(jt) ) THEN
          CALL get_channel_object(status,'mtskip', &
            TRIM(TRIM('mskt')//TRIM(ti_gp(jt)%tp%ident%fullname)), &
            p3=mtskiptend(jt)%ptr)
          CALL channel_halt(substr, status)
        ENDIF
      ENDDO

      ! op_pj_20160511+
      IF ((lfq_gp) .OR. (lfq_lg)) THEN
         CALL get_channel_object(status, 'mtskip', 'msktq', p3=mtskiptend_q)
         CALL channel_halt(substr, status)
      ENDIF
      ! op_pj_20160511-

      CALL start_message_bi(modstr, 'COUPLING TO MTSKIP (2)',substr) 

    END SELECT

  END SUBROUTINE mecca_mtskip_init_coupling
  !***************************************************************************
  !mz_bs_20150702-

END MODULE messy_mecca_si

!*****************************************************************************
