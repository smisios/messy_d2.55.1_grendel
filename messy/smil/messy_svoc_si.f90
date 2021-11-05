#include "messy_main_ppd_bi.inc"

! *************************************************************************
! *************************************************************************
! *************************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL SVOC
! 
! Author:
! Mega Octaviani, MPIC, 03.2017 (last modified)
!
! *************************************************************************
! *************************************************************************
! *************************************************************************

MODULE messy_svoc_si

! op_pj_20190326:
#if defined(ECHAM5) 

USE messy_main_constants_mem,    ONLY: dp
USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi, &
                                       error_bi, warning_bi, info_bi
USE messy_main_channel,          ONLY: t_chaobj_cpl
USE messy_main_tools,            ONLY: PTR_2D_ARRAY, PTR_3D_ARRAY, &
                                       PTR_4D_ARRAY
USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
USE messy_svoc ! smcl

IMPLICIT NONE

INTRINSIC :: NULL, TRIM, MAX, MIN, ABS, LOG, EXP
PRIVATE
SAVE

INTEGER, PUBLIC, PARAMETER :: nmaxsvoc = 100  ! max number of svocs allowed
INTEGER, PUBLIC, PARAMETER :: ncoeff = 6 ! number of pp-lfers coefficients
INTEGER, PUBLIC :: nsvoc = 0  ! number of svocs

! declarations for channel svoc_gp
  ! kp_m: partition coeff. by mode
  ! part_m: particle fraction by mode
  ! kp: partition coeff. (bulk)
  ! part: particle fraction (bulk)
  ! burd: compartment burden (1: soil, 2: veg., 3: oce., 4: snow,
  !       5: glacier, 6: atm)
  ! degr: degradation (1: soil, 2: veg., 3: oce.)
  ! depo: deposition (1: dry, 2: wet, 3: total)
  ! vola: volatilization (1: soil, 2: veg., 3: oce., 4: snow, 5: glacier)
  ! cwat: concentration of tracer in water
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: kp_m(:, :)
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: part_m(:, :)
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: kp(:)
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: part(:)
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: burd(:)
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: degr(:)
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: depo(:)
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: vola(:)
TYPE(PTR_2D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: cwat(:)

REAL(dp), DIMENSION(:, :, :), PUBLIC, POINTER :: depo_p => NULL()
REAL(dp), DIMENSION(:, :, :), PUBLIC, POINTER :: degr_p => NULL()
REAL(dp), DIMENSION(:, :, :), PUBLIC, POINTER :: burd_p => NULL()
REAL(dp), DIMENSION(:, :, :), PUBLIC, POINTER :: vola_p => NULL()
REAL(dp), DIMENSION(:, :, :), PUBLIC, POINTER :: kp_p => NULL()
REAL(dp), DIMENSION(:, :, :), PUBLIC, POINTER :: part_p => NULL()
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: cwat_p => NULL()

! declarations for channel svoc_inp
  ! appl: application-emission (prescribed emission)
  ! ocn_mld: ocean mixed layer depth
  ! cvs_old: land snow cover fraction at t-dt
  ! foms: fraction of organic matter in soil
  ! rhos: dry bulk density of soil
TYPE(PTR_2D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: appl(:)

REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: appl_p => NULL()
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: ocn_mld => NULL()
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: cvs_old => NULL()
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: foms => NULL()
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: rhos => NULL()

! read ss flux data (as mode)
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: in_ss_as_flux => NULL()

! read input data
TYPE(PTR_2D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: in_emiss(:)
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: in_omld => NULL()
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: in_foms => NULL()
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: in_rhos => NULL()

! coupling to other sm
  ! sigma: standard deviation of aerosols
  ! rwet_aer: wet radius of aerosols (m)
  ! aernl: aerosol number (cm-3; only for gmxe)
  ! vddepgas: dry deposition velocity of gaseous svoc
  ! vddepaer: dry deposition velocity of aerosols
  ! flux_airsea: air-sea transfer velocity
  ! wetflx_ls: total scavenging from stratiform precipitation (gas)
  ! wetflx_cv: total scavenging from convective precipitation (gas)
  ! wetflx_aer_tot: total scavenging of aerosols
  ! ttescav_bc: total scavenging tendency of bc
  ! ttescav_bc_cv: convective scavenging tendency of bc
  ! ttescav_bc_ls: large-scale scavenging tendency of bc
  ! ttescav_bc_ev: cloud evaporation tendency of bc
  ! tot_ttescav_cv: convective scavenging tendency of svoc particle (bulk)
  ! tot_ttescav_ls: large-scale scavenging tendency of svoc particle (bulk)
  ! ttescav_m_cv: convective scavenging tendency of svoc particle (by mode)
  ! ttescav_m_ls: large-scale scavenging tendency of svoc particle (by mode)
  ! trac_field: tracer field from cvtrans
  ! col: column from cvtrans (0-1)
REAL(dp), DIMENSION(:), PUBLIC, POINTER :: sigma => NULL()
REAL(dp), DIMENSION(:, :, :, :), PUBLIC, POINTER :: rwet_aer => NULL()
REAL(dp), DIMENSION(:, :, :, :), PUBLIC, POINTER :: aernl_cpl => NULL()
TYPE(PTR_2D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: vddepgas(:)
TYPE(PTR_2D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: vddepaer(:)
TYPE(PTR_2D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: flux_airsea(:)
REAL(dp), DIMENSION(:, :), PUBLIC, POINTER :: flux_airsea_p => NULL()
TYPE(PTR_2D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: wetflx_ls(:)
TYPE(PTR_2D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: wetflx_cv(:)
REAL(dp), DIMENSION(:, :), PUBLIC, ALLOCATABLE :: wetflx_aer_tot
TYPE(PTR_3D_ARRAY), PUBLIC, ALLOCATABLE, TARGET :: ttescav_bc(:)
REAL(dp), DIMENSION(:, :, :), PUBLIC, ALLOCATABLE :: ttescav_bc_cv
REAL(dp), DIMENSION(:, :, :), PUBLIC, ALLOCATABLE :: ttescav_bc_ls
REAL(dp), DIMENSION(:, :, :), PUBLIC, ALLOCATABLE :: ttescav_bc_ev
REAL(dp), DIMENSION(:, :, :), PUBLIC, ALLOCATABLE :: tot_ttescav_cv
REAL(dp), DIMENSION(:, :, :), PUBLIC, ALLOCATABLE :: tot_ttescav_ls
REAL(dp), DIMENSION(:, :, :, :), PUBLIC, ALLOCATABLE :: ttescav_m_cv
REAL(dp), DIMENSION(:, :, :, :), PUBLIC, ALLOCATABLE :: ttescav_m_ls
REAL(dp), DIMENSION(:, :, :, :), PUBLIC, POINTER :: trac_field => NULL()
REAL(dp), PUBLIC, POINTER :: col => NULL()
!! mz_jw_20160920+
REAL(dp), DIMENSION(:,:,:), PUBLIC, POINTER :: J_in_NO2 => NULL()
REAL(dp), DIMENSION(:,:,:), PUBLIC, POINTER :: J_out_NO2 => NULL()
!! mz_jw_20160920-
REAL(dp), DIMENSION(:,:), POINTER :: wind10_2d => NULL()

LOGICAL :: lcolumn = .FALSE., cvtrans_use = .FALSE.

! ozone index
INTEGER, PUBLIC :: idt_o3

! aerosol properties
INTEGER, PUBLIC, PARAMETER :: nmod = 7  !!! as in m7/gmxe
CHARACTER(LEN = 2), DIMENSION(nmod), PUBLIC :: cmode

DATA cmode / 'ns', & ! nucleation soluble
             'ks', & ! aitken soluble
             'as', & ! accumulation soluble
             'cs', & ! coarse soluble
             'ki', & ! aitken insoluble
             'ai', & ! accumulation insoluble
             'ci' /  ! coarse insoluble

! aerosol number (quantity in number density, 1/mol)
CHARACTER(LEN = 1), DIMENSION(nmod), PUBLIC :: bname_num = 'N'
INTEGER, DIMENSION(nmod), PUBLIC :: idt_nnum 

! see messy_gmxe_parameters_ori.inc to see list of permitted modes
INTEGER, PUBLIC, PARAMETER :: naer = 44
CHARACTER(LEN = STRLEN_MEDIUM), DIMENSION(naer), PUBLIC :: bname_aer
CHARACTER(LEN = 2), DIMENSION(naer), PUBLIC :: sname_aer
INTEGER, DIMENSION(naer), PUBLIC :: idt_maer

DATA bname_aer / 'SO4mm',  'SO4mm',  'SO4mm', 'SO4mm', &
                 'SO4mm',  'SO4mm',  'SO4mm',          & ! SO4--
                 'HSO4m',  'HSO4m',  'HSO4m', 'HSO4m', &
                 'HSO4m',  'HSO4m',  'HSO4m',          & ! HSO4-
                 'BC',     'BC',     'BC',    'BC',    & ! black carbon
                 'OC',     'OC',     'OC',    'OC',    & ! organic carbon
                 'SS',     'SS',     'SS',             & ! sea salt
                 'DU',     'DU',     'DU',    'DU',    & ! mineral dust
!l_oc_aging in gmxe.nml must set to .TRUE. (and num_wsoc = 5)
                 'WSOC01', 'WSOC01', 'WSOC01',         &
                 'WSOC02', 'WSOC02', 'WSOC02',         &
                 'WSOC03', 'WSOC03', 'WSOC03',         &
                 'WSOC04', 'WSOC04', 'WSOC04',         &
                 'WSOC05', 'WSOC05', 'WSOC05' /

DATA sname_aer / 'ns',     'ks',     'as',    'cs',     &
                 'ki',     'ai',     'ci',              &
                 'ns',     'ks',     'as',    'cs',     &
                 'ki',     'ai',     'ci',              &
                 'ki',     'ks',     'as',    'cs',     &
                 'ki',     'ks',     'as',    'cs',     &
                 'ks',     'as',     'cs',              &
                 'ai',     'ci',     'as',    'cs',     &
!l_oc_aging in gmxe.nml must set to .TRUE. (and num_wsoc = 5)
                 'ks',     'as',     'cs',              &
                 'ks',     'as',     'cs',              &
                 'ks',     'as',     'cs',              &
                 'ks',     'as',     'cs',              &
                 'ks',     'as',     'cs' /

! variables on cpl-namelist
LOGICAL, PUBLIC :: l_mode_partition
INTEGER, PUBLIC :: param_part
INTEGER, PUBLIC :: param_soilv
INTEGER, PUBLIC :: param_bapo3
CHARACTER(LEN = STRLEN_MEDIUM), PUBLIC :: aermod_str

  ! information on import variables
TYPE(t_chaobj_cpl) :: imp_ss_as_flux
TYPE(t_chaobj_cpl) :: imp_om_soil
TYPE(t_chaobj_cpl) :: imp_rho_soil
TYPE(t_chaobj_cpl) :: imp_mld_oce

  ! information on svoc physicochemical properties
CHARACTER(LEN = STRLEN_MEDIUM), DIMENSION(nmaxsvoc), PUBLIC :: SVOC_NAME
TYPE(t_chaobj_cpl), DIMENSION(nmaxsvoc), PUBLIC :: EMISS_IN
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: MOLMASS
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: MOLVOL
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RKSOIL
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RKOCEAN
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RWSOL
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RVAPP
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RHSOL
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RHVAP
!! mz_jw_20170220+
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RHABC
!! mz_jw_20170220-
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RHSUB
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RLOGKOW
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RKOAM
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RKOAB
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: KAHENRY
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: KBHENRY
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RLOSS
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RSPRAY
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: RDENS
REAL(dp), DIMENSION(nmaxsvoc), PUBLIC :: REMISP
REAL(dp), DIMENSION(nmaxsvoc, nmod), PUBLIC :: REMIS_MOD
REAL(dp), DIMENSION(nmaxsvoc, ncoeff), PUBLIC :: RPPLFER

LOGICAL, PUBLIC :: l_nsvoc = .FALSE.
LOGICAL, PUBLIC :: l_nsvoc_lg_tot = .FALSE.
LOGICAL, PUBLIC :: l_nsvoc_gp_tot = .FALSE.
LOGICAL, DIMENSION(nmaxsvoc), PUBLIC :: l_svoc_lg = .FALSE.
LOGICAL, DIMENSION(nmaxsvoc), PUBLIC :: l_svoc_gp = .TRUE.
INTEGER, DIMENSION(nmaxsvoc), PUBLIC :: idt_svocg = 0 ! svoc gas index
INTEGER, DIMENSION(nmaxsvoc), PUBLIC :: idt_svocp = 0 ! svoc particle (as bulk) index
INTEGER, DIMENSION(nmaxsvoc, nmod), PUBLIC :: idt_svocm = 0 ! svoc particle (by mode) index
LOGICAL, DIMENSION(nmaxsvoc), PUBLIC :: l_emiss = .FALSE. ! primary emission data

! public routines
PUBLIC :: svoc_initialize    ! initialize submodel, read namelists
PUBLIC :: svoc_new_tracer    ! define svoc particles as tracers
PUBLIC :: svoc_init_memory   ! channel definition and local memory alloc
PUBLIC :: svoc_init_coupling ! set pointers for coupling to bm and sm
PUBLIC :: svoc_vdiff         ! distribute emissions, partitioning, exchange
PUBLIC :: svoc_convec        ! convective scavenging of aerosols
PUBLIC :: svoc_physc         ! simple heterogeneous chemistry, wet deposition
PUBLIC :: svoc_free_memory   ! local memory deallocation

CONTAINS

! *************************************************************************
! PUBLIC SUBROUTINES
! *************************************************************************

SUBROUTINE svoc_initialize

!!! to read CTRL- and CPL- namelists and broadcast variables

USE messy_main_mpi_bi,          ONLY: p_parallel, p_parallel_io, p_io, &
                                      p_bcast !!$, message ! op_pj_20190326
USE messy_main_tools,           ONLY: find_next_free_unit
!!$USE mo_util_string,             ONLY: separator ! op_pj_20190326
USE messy_main_tracer_mem_bi,   ONLY: NGCELL


IMPLICIT NONE

CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_initialize'
INTEGER :: status ! error status
INTEGER :: iou    ! I/O unit
INTEGER :: jt, jm ! counter

CALL start_message_bi(modstr, 'SVOC INITIALIZATION', substr)

! op_pj_20190327+
#ifdef MESSYTENDENCY
CALL error_bi('SVOC is not enabled for MESSYTENDENCY',substr)
#endif
! op_pj_20190327-

! initialize CTRL
IF (p_parallel_io) THEN
   iou = find_next_free_unit(100, 200)
   CALL svoc_read_nml_ctrl(status, iou)
   IF (status /= 0) CALL error_bi('Error in reading CTRL namelist', substr)
END IF

IF (L_LG .AND. (NGCELL > 0)) THEN
   L_LG = .TRUE.
ELSE
   IF (L_LG) THEN
      WRITE(*, *) 'L_LG = TRUE in namelist'
      WRITE(*, *) 'But no Lagrangian scheme is activated'
      WRITE(*, *) '... setting L_LG = F'
   END IF

   L_LG = .FALSE.
END IF

! broadcast over processors
IF (p_parallel) THEN
   CALL p_bcast(L_GP, p_io)
   CALL p_bcast(L_LG, p_io)
   CALL p_bcast(l_svocpart, p_io)
   CALL p_bcast(l_svocvola, p_io)
   CALL p_bcast(l_glacier, p_io)
   CALL p_bcast(l_landsnow, p_io)
   CALL p_bcast(l_pahderiv, p_io)
END IF

! setup default values for CPL
l_mode_partition = .FALSE.
param_part = 1
param_soilv = 1
param_bapo3 = 1
aermod_str = ''

imp_ss_as_flux%cha = ''
imp_ss_as_flux%obj = ''
imp_om_soil%cha = ''
imp_om_soil%obj = ''
imp_rho_soil%cha = ''
imp_rho_soil%obj = ''
imp_mld_oce%cha = ''
imp_mld_oce%obj = ''

DO jt = 1, nmaxsvoc
   SVOC_NAME(jt) = ''
   EMISS_IN(jt)%cha = ''
   EMISS_IN(jt)%obj = ''
   MOLMASS(jt) = 0._dp
   MOLVOL(jt) = 0._dp
   RKSOIL(jt) = 0._dp
   RKOCEAN(jt) = 0._dp
   RWSOL(jt) = 0._dp
   RVAPP(jt) = 0._dp
   RHSOL(jt) = 0._dp
   RHVAP(jt) = 0._dp
   RHSUB(jt) = 0._dp
   RLOGKOW(jt) = 1._dp
   RKOAM(jt) = 0._dp
   RKOAB(jt) = 0._dp
   KAHENRY(jt) = 1._dp
   KBHENRY(jt) = 0._dp
   RLOSS(jt) = 1._dp
   RSPRAY(jt) = 0._dp
   RDENS(jt) = 0._dp
   REMISP(jt) = 0._dp
!! mz_jw_20170220+
   RHABC(jt) = 0._dp
!! mz_jw_20170220-
   REMIS_MOD(jt, :) = 0._dp
   RPPLFER(jt, :) = 0._dp
END DO

! initialize CPL
IF (p_parallel_io) THEN
   iou = find_next_free_unit(100, 200)
   CALL svoc_read_nml_cpl(status, iou)
   IF (status /= 0) CALL error_bi('Error in reading CPL namelist', substr)

   DO jt = 1, nmaxsvoc
      IF (TRIM(SVOC_NAME(jt)) == '') CYCLE
      nsvoc = nsvoc + 1

      SVOC_NAME(nsvoc) = TRIM(SVOC_NAME(jt))
      EMISS_IN(nsvoc)%cha = EMISS_IN(jt)%cha
      EMISS_IN(nsvoc)%obj = EMISS_IN(jt)%obj
      MOLMASS(nsvoc) = MOLMASS(jt)
      MOLVOL(nsvoc) = MOLVOL(jt)
      RKSOIL(nsvoc) = RKSOIL(jt)
      RKOCEAN(nsvoc) = RKOCEAN(jt)
      RWSOL(nsvoc) = RWSOL(jt)
      RVAPP(nsvoc) = RVAPP(jt)
      RHSOL(nsvoc) = RHSOL(jt)
      RHVAP(nsvoc) = RHVAP(jt)
      RHSUB(nsvoc) = RHSUB(jt)
      RLOGKOW(nsvoc) = RLOGKOW(jt)
      RKOAM(nsvoc) = RKOAM(jt)
      RKOAB(nsvoc) = RKOAB(jt)
      KAHENRY(nsvoc) = KAHENRY(jt)
      KBHENRY(nsvoc) = KBHENRY(jt)
      RLOSS(nsvoc) = RLOSS(jt)
      RSPRAY(nsvoc) = RSPRAY(jt)
      RDENS(nsvoc) = RDENS(jt)
      REMISP(nsvoc) = REMISP(jt)
!! mz_jw_20170220+
      RHABC(nsvoc) = RHABC(jt)
!! mz_jw_20170220-
      REMIS_MOD(nsvoc, :) = REMIS_MOD(jt, 1: nmod)
      RPPLFER(nsvoc, :) = RPPLFER(jt, 1: ncoeff)
   END DO
END IF

! broadcast over processors
CALL p_bcast(l_mode_partition, p_io)
CALL p_bcast(param_part, p_io)
CALL p_bcast(param_soilv, p_io)
CALL p_bcast(param_bapo3, p_io)
CALL p_bcast(aermod_str, p_io)
CALL p_bcast(imp_ss_as_flux%cha, p_io)
CALL p_bcast(imp_ss_as_flux%obj, p_io)
CALL p_bcast(imp_om_soil%cha, p_io)
CALL p_bcast(imp_om_soil%obj, p_io)
CALL p_bcast(imp_rho_soil%cha, p_io)
CALL p_bcast(imp_rho_soil%obj, p_io)
CALL p_bcast(imp_mld_oce%cha, p_io)
CALL p_bcast(imp_mld_oce%obj, p_io)

CALL p_bcast(nsvoc, p_io)

DO jt = 1, nsvoc
   CALL p_bcast(SVOC_NAME(jt), p_io)
   CALL p_bcast(EMISS_IN(jt)%cha, p_io)
   CALL p_bcast(EMISS_IN(jt)%obj, p_io)
   CALL p_bcast(MOLMASS(jt), p_io)
   CALL p_bcast(MOLVOL(jt), p_io)
   CALL p_bcast(RKSOIL(jt), p_io)
   CALL p_bcast(RKOCEAN(jt), p_io)
   CALL p_bcast(RWSOL(jt), p_io)
   CALL p_bcast(RVAPP(jt), p_io)
   CALL p_bcast(RHSOL(jt), p_io)
   CALL p_bcast(RHVAP(jt), p_io)
   CALL p_bcast(RHSUB(jt), p_io)
   CALL p_bcast(RLOGKOW(jt), p_io)
   CALL p_bcast(RKOAM(jt), p_io)
   CALL p_bcast(RKOAB(jt), p_io)
   CALL p_bcast(KAHENRY(jt), p_io)
   CALL p_bcast(KBHENRY(jt), p_io)
   CALL p_bcast(RLOSS(jt), p_io)
   CALL p_bcast(RSPRAY(jt), p_io)
   CALL p_bcast(RDENS(jt), p_io)
   CALL p_bcast(REMISP(jt), p_io)
!! mz_jw_20170220+
   CALL p_bcast(RHABC(jt), p_io)
!! mz_jw_20170220-
   CALL p_bcast(REMIS_MOD(jt, :), p_io)
   CALL p_bcast(RPPLFER(jt, :), p_io)
END DO

IF (p_parallel_io) THEN
   WRITE(*, *) '-----------------------------------------------------'
   WRITE(*, *) '-----------------------------------------------------'
   WRITE(*, *) 'SVOC submodel setup: '

   IF (L_GP) THEN
      WRITE(*, *) 'SVOC IN GRIDPOINT SPACE : ON'
   ELSE
      WRITE(*, *) 'SVOC IN GRIDPOINT SPACE : OFF'
   END IF

   IF (L_LG) THEN
      WRITE(*, *) 'SVOC IN LAGRANGIAN SPACE : ON'
   ELSE
      WRITE(*, *) 'SVOC IN LAGRANGIAN SPACE : OFF'
   END IF

   IF (.NOT. l_svocpart) THEN
      WRITE(*, *) '   Gas-particle partitioning: OFF'
   ELSE
      WRITE(*, *) '   Gas-particle partitioning: ON'
   END IF

   IF (.NOT. l_svocvola) THEN
      WRITE(*, *) '   Surface volatilizations: OFF'
   ELSE
      WRITE(*, *) '   Surface volatilizations: ON'
   END IF

   IF (.NOT. l_glacier) THEN
      WRITE(*, *) '   No glacier compartment'
   ELSE
      WRITE(*, *) '   Include glacier compartment'
   END IF

   IF (.NOT. l_landsnow) THEN
      WRITE(*, *) '   No landsnow compartment'
   ELSE
      WRITE(*, *) '   Include snow compartment'
   END IF

   IF (.NOT. l_pahderiv) THEN
      WRITE(*, *) '   No Nitro-PAHs'
   ELSE
      WRITE(*, *) '   Include Nitro-PAHs'
   ENDIF

   IF (.NOT. l_mode_partition) THEN
      WRITE(*, *) '   Use Bulk Partitioning'
   ELSE
      WRITE(*, *) '   Use Partitioning by Mode'
   END IF

   WRITE(*, *) '   Partitioning scheme: ', param_part              
   WRITE(*, *) '   Soil model: ', param_soilv
   WRITE(*, *) '   On-particle BaP Oxidation by O3: ', param_bapo3
   WRITE(*, *) '   Aerosol model:', TRIM(aermod_str)

   WRITE(*, *) '   Number of SVOCs: ', nsvoc

   DO jt = 1, nsvoc
      WRITE(*, *) '      NO. NAME = ', jt, TRIM(SVOC_NAME(jt))
   END DO

! op_pj_20190326+
!!$   CALL message('', separator)
!!$   CALL message('', separator)
   WRITE(*, *) '-------------------------------------------------------------'
   WRITE(*, *) '-------------------------------------------------------------'
! op_pj_20190326-
END IF

CALL end_message_bi(modstr, 'SVOC INITIALIZATION', substr)

END SUBROUTINE svoc_initialize

! *************************************************************************

SUBROUTINE svoc_new_tracer

!!! request svoc particle tracers and prescribe their properties

USE messy_main_mpi_bi,          ONLY: p_parallel_io
USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
USE messy_main_tracer_tools_bi, ONLY: tracer_halt ! error handling
USE messy_main_tracer,          ONLY: new_tracer, set_tracer, &
                                      AEROSOL, AMOUNTFRACTION, &
                                      R_MOLARMASS, R_AEROSOL_DENSITY, &
                                      I_AEROSOL_METHOD, I_AEROSOL_MODE, &
                                      S_AEROSOL_MODEL, I_DRYDEP, I_SEDI, &
                                      I_SCAV, MODAL, ON, OFF

IMPLICIT NONE

CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_new_tracer'
INTEGER :: status
INTEGER :: jt, jm

CALL start_message_bi(modstr, 'SVOC PARTICLE TRACER DEFINED', substr)

IF (.NOT. l_mode_partition) THEN
   DO jt = 1, nsvoc
      status = 0

      CALL new_tracer(status, GPTRSTR, TRIM(SVOC_NAME(jt)) // 'part', &
                      modstr, idx = idt_svocp(jt), &
                      type = 0, unit = 'mol/mol', &
                      medium = AEROSOL, quantity = AMOUNTFRACTION)
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_svocp(jt), R_MOLARMASS, &
                      r = MOLMASS(jt))
      CALL tracer_halt(substr, status)

      ! register the tracer in a new aerosol model
      CALL set_tracer(status, GPTRSTR, idt_svocp(jt), S_AEROSOL_MODEL, &
                      TRIM(modstr))
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_svocp(jt), I_AEROSOL_METHOD, MODAL)
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_svocp(jt), I_AEROSOL_MODE, 1)
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_svocp(jt), I_DRYDEP, OFF)
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_svocp(jt), I_SEDI, OFF)
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_svocp(jt), I_SCAV, OFF)
      CALL tracer_halt(substr, status)

      IF (p_parallel_io) &
         WRITE(*, *) ' ... New defined tracer:', &
            TRIM(SVOC_NAME(jt)) // 'part'

!!$#ifdef MESSYTENDENCY
!!$      CALL mtend_register(my_handle, mtend_id_tracer, idt = idt_svocp(jt))
!!$#endif
   END DO

ELSE !!! define new tracer per aerosol mode 
   DO jt = 1, nsvoc
      DO jm = 1, nmod
         status = 0

         CALL new_tracer(status, GPTRSTR, TRIM(SVOC_NAME(jt)) // &
                         TRIM(cmode(jm)), modstr, &
                         idx = idt_svocm(jt, jm), &
                         type = 0, unit = 'mol/mol', &
                         medium = AEROSOL, quantity = AMOUNTFRACTION)
         CALL tracer_halt(substr, status)

         CALL set_tracer(status, GPTRSTR, idt_svocm(jt, jm), R_MOLARMASS, &
                         r = MOLMASS(jt))
         CALL tracer_halt(substr, status)

         CALL set_tracer(status, GPTRSTR, idt_svocm(jt, jm), &
                         R_AEROSOL_DENSITY, r = RDENS(jt))
         CALL tracer_halt(substr, status)

         ! register the tracer in a new aerosol model
         CALL set_tracer(status, GPTRSTR, idt_svocm(jt, jm), &
                         S_AEROSOL_MODEL, TRIM(modstr))
         CALL tracer_halt(substr, status)

         CALL set_tracer(status, GPTRSTR, idt_svocm(jt, jm), &
                         I_AEROSOL_METHOD, MODAL)
         CALL tracer_halt(substr, status)

         CALL set_tracer(status, GPTRSTR, idt_svocm(jt, jm), &
                         I_AEROSOL_MODE, jm)
         CALL tracer_halt(substr, status)

         CALL set_tracer(status, GPTRSTR, idt_svocm(jt, jm), I_DRYDEP, OFF)
         CALL tracer_halt(substr, status)

         CALL set_tracer(status, GPTRSTR, idt_svocm(jt, jm), I_SEDI, OFF)
         CALL tracer_halt(substr, status)

         CALL set_tracer(status, GPTRSTR, idt_svocm(jt, jm), I_SCAV, OFF)
         CALL tracer_halt(substr, status)

         IF (p_parallel_io) &
            WRITE(*, *) ' ... New defined tracer:', &
               TRIM(SVOC_NAME(jt)) // TRIM(cmode(jm))

!!$#ifdef MESSYTENDENCY
!!$         CALL mtend_register(my_handle, mtend_id_tracer, idt = idt_svocm(jt, jm))
!!$#endif

      END DO
   END DO
END IF

CALL end_message_bi(modstr, 'SVOC PARTICLE TRACER DEFINED', substr)

END SUBROUTINE svoc_new_tracer

! *************************************************************************

SUBROUTINE svoc_init_memory

!!! initialize memory

USE messy_main_mpi_bi,           ONLY: p_parallel_io
USE messy_main_grid_def_mem_bi,  ONLY: nproma, nlev, ngpblks
USE messy_main_tracer_mem_bi,    ONLY: ntrac_gp, ti_gp, GPTRSTR
USE messy_main_blather_bi,       ONLY: info_bi
USE messy_main_channel_error_bi, ONLY: channel_halt
USE messy_main_channel_bi,       ONLY: GP_3D_MID, &
                                       GP_2D_HORIZONTAL, DC_BC, DC_GP, &
                                       DIMID_LON, DIMID_LAT, DIMID_LEV, &
                                       gp_nseg, gp_start, gp_cnt, &
                                       gp_meml, gp_memu
USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                       new_attribute
USE messy_main_channel_repr,     ONLY: new_representation, AUTO, &
                                       set_representation_decomp, &
                                       IRANK, PIOTYPE_COL
USE messy_main_channel_dimensions, ONLY: new_dimension

IMPLICIT NONE

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_init_memory'
CHARACTER(LEN = *), PARAMETER :: modstr_gp = TRIM(modstr) // '_' // GPTRSTR
INTEGER :: status ! error status
INTEGER :: jt, js, jm ! counter
CHARACTER(len = 1) :: cjm

! channel management
INTEGER :: DIMID_NMODE
INTEGER :: REPR_SVOC_1D, REPR_SVOC_4D

! parallel decomposition
INTEGER :: nseg = 0
INTEGER :: nmodsvoc
INTEGER, DIMENSION(:, :), POINTER :: start => NULL()
INTEGER, DIMENSION(:, :), POINTER :: cnt   => NULL()
INTEGER, DIMENSION(:, :), POINTER :: meml  => NULL()
INTEGER, DIMENSION(:, :), POINTER :: memu  => NULL()

! channel objects for svoc aerosols
REAL(dp), DIMENSION(:, :, :, :), POINTER :: wetrad_p => NULL()
REAL(dp), DIMENSION(:, :, :, :), POINTER :: aerdens_p => NULL()
REAL(dp), DIMENSION(:), POINTER :: sigma_p => NULL()
REAL(dp) :: aero_sigma = 2._dp !!! dummy parameter
REAL(dp) :: aero_rad = 1e-6_dp   !!! dummy parameter
REAL(dp) :: aero_dens = 0._dp  !!! dummy parameter

CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

IF (.NOT. l_mode_partition) THEN
   nmodsvoc = 1
ELSE
   nmodsvoc = nmod
END IF

! new representations
CALL new_dimension(status, DIMID_NMODE, 'SVOC_NMODE', nmodsvoc)
CALL channel_halt(substr, status)

CALL new_representation(status, REPR_SVOC_1D, 'REPR_SVOC_1D', &
                        rank = 1, link = 'x---', dctype = DC_BC, &
                        dimension_ids = (/ DIMID_NMODE /), &
                        ldimlen = (/ AUTO /), &
                        axis = 'N---')
CALL channel_halt(substr, status)

nseg = 1
ALLOCATE(start(nseg, IRANK))
ALLOCATE(cnt(nseg, IRANK))
ALLOCATE(meml(nseg, IRANK))
ALLOCATE(memu(nseg, IRANK))

start(:, :) = 1
cnt(:, :) = 1
meml(:, :) = 1
memu(:, :) = 1

cnt(:, 1) = nmodsvoc
memu(:, 1) = nmodsvoc

CALL set_representation_decomp(status, REPR_SVOC_1D, &
                               start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
CALL channel_halt(substr, status)

DEALLOCATE(start) ; NULLIFY(start)
DEALLOCATE(cnt)   ; NULLIFY(cnt)
DEALLOCATE(meml)  ; NULLIFY(meml)
DEALLOCATE(memu)  ; NULLIFY(memu)

CALL new_representation(status, REPR_SVOC_4D, 'REPR_SVOC_4D', &
                        rank = 4, link = 'xxxx', dctype = DC_GP, &
                        dimension_ids = (/ DIMID_LON, DIMID_LEV, &
                        DIMID_NMODE, DIMID_LAT /), &
                        ldimlen = (/ nproma, AUTO, AUTO, ngpblks /), &
                        output_order = (/ 3, 1, 4, 2 /), axis = 'XZNY')

nseg = gp_nseg

ALLOCATE(start(nseg, IRANK))
ALLOCATE(cnt(nseg, IRANK))
ALLOCATE(meml(nseg, IRANK))
ALLOCATE(memu(nseg, IRANK))

start(:, :) = gp_start(:, :)
cnt(:, :) = gp_cnt(:, :)
meml(:, :) = gp_meml(:, :)
memu(:, :) = gp_memu(:, :)

cnt(:,3) = nmodsvoc
memu(:,3) = nmodsvoc

CALL set_representation_decomp(status, REPR_SVOC_4D, &
                              start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
CALL channel_halt(substr, status)

DEALLOCATE(start) ; NULLIFY(start)
DEALLOCATE(cnt)   ; NULLIFY(cnt)
DEALLOCATE(meml)  ; NULLIFY(meml)
DEALLOCATE(memu)  ; NULLIFY(memu)

! get tracer information
IF (L_GP) THEN
   DO jt = 1, nsvoc
      l_svoc_gp(jt) = .FALSE.
   END DO

   DO jt = 1, nsvoc
      DO js = 1, ntrac_gp
         IF(ti_gp(js)%tp%ident%basename /= SVOC_NAME(jt)) CYCLE

         WRITE(*, *) ti_gp(js)%tp%ident%basename, SVOC_NAME(jt)
         idt_svocg(jt) = js
         l_svoc_gp(jt) = .TRUE.
         l_nsvoc_gp_tot = .TRUE.
      END DO
   END DO
END IF

! check if you have tracers for the calc.
IF (l_nsvoc_gp_tot) THEN
   l_nsvoc = .TRUE.
ELSE
   L_GP = .FALSE.
END IF

IF (l_nsvoc_lg_tot) THEN
   l_nsvoc = .TRUE.
ELSE
   L_LG = .FALSE.
END IF

IF (.NOT. l_nsvoc) CALL info_bi('No tracer for SVOC application!', substr)

! allocate scav arrays
ALLOCATE(wetflx_aer_tot(nsvoc, nproma))
ALLOCATE(ttescav_m_cv(nsvoc, nproma, nlev, nmod))
ALLOCATE(ttescav_m_ls(nsvoc, nproma, nlev, nmod))
ALLOCATE(tot_ttescav_cv(nsvoc, nproma, nlev))
ALLOCATE(tot_ttescav_ls(nsvoc, nproma, nlev))
ALLOCATE(ttescav_bc_cv(4, nproma, nlev))
ALLOCATE(ttescav_bc_ls(4, nproma, nlev))
ALLOCATE(ttescav_bc_ev(4, nproma, nlev))

! construct new channel: svoc_gp
ALLOCATE(kp_m(nsvoc, nmod))
ALLOCATE(part_m(nsvoc, nmod))
ALLOCATE(kp(nsvoc))
ALLOCATE(part(nsvoc))
ALLOCATE(burd(nsvoc))
ALLOCATE(vola(nsvoc))
ALLOCATE(depo(nsvoc))
ALLOCATE(degr(nsvoc))
ALLOCATE(cwat(nsvoc))

CALL new_channel(status, modstr_gp, reprid = GP_3D_MID, lrestreq = .TRUE.)
CALL channel_halt(substr, status)

DO jt = 1, nsvoc
   ! add channel object: partitioning coefficient [m3 ug-1]
   CALL new_channel_object(status, modstr_gp, &
                           'Kp_' // TRIM(SVOC_NAME(jt)), &
                           p3 = kp(jt)%ptr)
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'Kp_' // TRIM(SVOC_NAME(jt)), &
                      'long_name', c = 'partition coeff. (bulk)')
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'Kp_' // TRIM(SVOC_NAME(jt)), &
                      'units', c = 'm3 ug-1')
   CALL channel_halt(substr, status)
 
   IF (l_mode_partition) THEN
      DO jm = 1, nmod
         WRITE (cjm, '(i1)') jm

         CALL new_channel_object(status, modstr_gp, &
                                 'Kp_' // TRIM(SVOC_NAME(jt)) // &
                                 '_' // TRIM(cmode(jm)), &
                                 p3 = kp_m(jt, jm)%ptr)
         CALL channel_halt(substr, status)

         CALL new_attribute(status, modstr_gp, &
                            'Kp_' // TRIM(SVOC_NAME(jt)) // &
                            '_' // TRIM(cmode(jm)), 'long_name', &
                            c = 'partition coeff. of mode ' // TRIM(cjm))
         CALL channel_halt(substr, status)

         CALL new_attribute(status, modstr_gp, &
                            'Kp_' // TRIM(SVOC_NAME(jt)) // &
                            '_' // TRIM(cmode(jm)), &
                            'units', c = 'm3 ug-1')
         CALL channel_halt(substr, status)
      END DO
   END IF

   IF (p_parallel_io) WRITE(*, *) '... Kp added to channel ', modstr_gp

   ! add channel object: particle fraction

   CALL new_channel_object(status, modstr_gp, &
                           'PART_' // TRIM(SVOC_NAME(jt)), &
                           p3 = part(jt)%ptr)
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'PART_' // TRIM(SVOC_NAME(jt)), &
                      'long_name', c = 'particle (bulk) fraction')
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'PART_' // TRIM(SVOC_NAME(jt)), &
                      'units', c = '-')
   CALL channel_halt(substr, status)
    
   IF (l_mode_partition) THEN
      DO jm = 1, nmod
         WRITE (cjm, '(i1)') jm

         CALL new_channel_object(status, modstr_gp, &
                                 'PART_' // TRIM(SVOC_NAME(jt)) // &
                                 '_' // TRIM(cmode(jm)), &
                                 p3 = part_m(jt, jm)%ptr)
         CALL channel_halt(substr, status)

         CALL new_attribute(status, modstr_gp, &
                            'PART_' // TRIM(SVOC_NAME(jt)) // &
                            '_' // TRIM(cmode(jm)), 'long_name', &
                            c = 'particle fraction at mode ' // TRIM(cjm))
         CALL channel_halt(substr, status)

         CALL new_attribute(status, modstr_gp, &
                            'PART_' // TRIM(SVOC_NAME(jt)) // &
                            '_' // TRIM(cmode(jm)), &
                            'units', c = '-')
         CALL channel_halt(substr, status)
      END DO
   END IF

   IF (p_parallel_io) WRITE(*, *) '... PART added to channel ', modstr_gp

   ! add channel object: compartmental burden
   CALL new_channel_object(status, modstr_gp, &
                           'BURD_' // TRIM(SVOC_NAME(jt)), &
                           p3 = burd(jt)%ptr)
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'BURD_' // TRIM(SVOC_NAME(jt)), &
                      'long_name', c = 'compartmental burden')
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'BURD_' // TRIM(SVOC_NAME(jt)), &
                      'units', c = 'kg(trac) m-2')
   CALL channel_halt(substr, status)

   IF (p_parallel_io) WRITE(*, *) '... BURD added to channel ', modstr_gp

   ! add channel object: volatilization
   CALL new_channel_object(status, modstr_gp, &
                           'VOLA_' // TRIM(SVOC_NAME(jt)), &
                           p3 = vola(jt)%ptr)
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'VOLA_' // TRIM(SVOC_NAME(jt)), &
                      'long_name', c = 'volatilization at surfaces')
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'VOLA_' // TRIM(SVOC_NAME(jt)), &
                      'units', c = 'kg(trac) m-2 s-1')
   CALL channel_halt(substr, status)

   IF (p_parallel_io) WRITE(*, *) '... VOLA added to channel ', modstr_gp

   ! add channel object: deposition
   CALL new_channel_object(status, modstr_gp, &
                           'DEPO_' // TRIM(SVOC_NAME(jt)), &
                           p3 = depo(jt)%ptr)
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'DEPO_' // TRIM(SVOC_NAME(jt)), &
                      'long_name', c = 'deposition')
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'DEPO_' // TRIM(SVOC_NAME(jt)), &
                      'units', c = 'kg(trac) m-2 s-1')
   CALL channel_halt(substr, status)

   IF (p_parallel_io) WRITE(*, *) '... DEPO added to channel ', modstr_gp

   ! add channel object: degradation
   CALL new_channel_object(status, modstr_gp, &
                           'DEGR_' // TRIM(SVOC_NAME(jt)), &
                           p3 = degr(jt)%ptr)
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'DEGR_' // TRIM(SVOC_NAME(jt)), &
                      'long_name', c = 'degradation')
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'DEGR_' // TRIM(SVOC_NAME(jt)), &
                      'units', c = 'kg(trac) m-2 s-1')
   CALL channel_halt(substr, status)

   IF (p_parallel_io) WRITE(*, *) '... DEGR added to channel ', modstr_gp

   ! add channel object: svoc concentration in water
   CALL new_channel_object(status, modstr_gp, &
                           'CWAT_' // TRIM(SVOC_NAME(jt)), &
                           reprid = GP_2D_HORIZONTAL, &
                           p2 = cwat(jt)%ptr)
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'CWAT_' // TRIM(SVOC_NAME(jt)), &
                      'long_name', c = 'water concentration')
   CALL channel_halt(substr, status)

   CALL new_attribute(status, modstr_gp, &
                      'CWAT_' // TRIM(SVOC_NAME(jt)), &
                      'units', c = 'kmol m-3')
   CALL channel_halt(substr, status)

   IF (p_parallel_io) WRITE(*, *) '... CWAT added to channel ', modstr_gp
END DO

! add channel object: standard deviation of aerosol mode
CALL new_channel_object(status, modstr_gp, &
                        'sigma', reprid = REPR_SVOC_1D, &
                        p1 = sigma_p)
CALL channel_halt(substr, status)

CALL new_attribute(status, modstr_gp, 'sigma', 'long_name', &
                   c = 'standard deviation of aerosol mode')
CALL channel_halt(substr, status)

CALL new_attribute(status, modstr_gp, 'sigma', 'units', c = '-')
CALL channel_halt(substr, status)

IF (p_parallel_io) WRITE(*, *) '... sigma added to channel ', modstr_gp

! add channel object: const. ambient aerosol radius
CALL new_channel_object(status, modstr_gp, 'wetradius', &
                        reprid = REPR_SVOC_4D, p4 = wetrad_p)
CALL channel_halt(substr, status)

CALL new_attribute(status, modstr_gp, 'wetradius', 'long_name', &
                   c = 'const. ambient aerosol radius')
CALL channel_halt(substr, status)

CALL new_attribute(status, modstr_gp, 'wetradius', 'units', c = 'm')
CALL channel_halt(substr, status)

IF (p_parallel_io) WRITE(*, *) '... wetradius added to channel ', modstr_gp

! add channel object: aerosol density
CALL new_channel_object(status, modstr_gp, 'densaer', &
                        reprid = REPR_SVOC_4D, p4 = aerdens_p)
CALL channel_halt(substr, status)

CALL new_attribute(status, modstr_gp, 'densaer', &
                   'long_name', c = 'aerosol density')
CALL channel_halt(substr, status)

CALL new_attribute(status, modstr_gp, 'densaer', 'units', c = 'kg m-3')
CALL channel_halt(substr, status)

IF (p_parallel_io) WRITE(*, *) '... densaer added to channel ', modstr_gp

! initialize radius, sigma, and density
DO jm = 1, nmodsvoc
   wetrad_p(:, :, jm, :) = aero_rad
   aerdens_p(:, :, jm, :) = aero_dens
   sigma_p(jm) = aero_sigma
END DO

! construct new channel: svoc_inp
ALLOCATE(appl(nsvoc))

CALL new_channel(status, TRIM(modstr) // '_inp', reprid = GP_2D_HORIZONTAL, &
                 lrestreq = .TRUE.)
CALL channel_halt(substr, status)

DO jt = 1, nsvoc
   ! add channel object: application/emission
   CALL new_channel_object(status, TRIM(modstr) // '_inp', &
                           'appl_' // TRIM(SVOC_NAME(jt)), &
                           p2 = appl(jt)%ptr)
   CALL channel_halt(substr, status)

   CALL new_attribute(status, TRIM(modstr) // '_inp', &
                      'appl_' // TRIM(SVOC_NAME(jt)), 'long_name', &
                      c = 'input application of ' // TRIM(SVOC_NAME(jt)))
   CALL channel_halt(substr, status)

   CALL new_attribute(status, TRIM(modstr) // '_inp', &
                      'appl_' // TRIM(SVOC_NAME(jt)), 'units', &
                      c = 'mol(trac) mol(air)-1 kg(air) m-2 s-1')
   CALL channel_halt(substr, status)

   IF (p_parallel_io) &
      WRITE(*, *) '... appl added to channel ', TRIM(modstr) // '_inp'
END DO

! add channel object: ocean mixed layer depth
CALL new_channel_object(status, TRIM(modstr) // '_inp', 'ocn_mld', &
                        p2 = ocn_mld)
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'ocn_mld', & 
                  'long_name',  c = 'ocean mixed layer depth')
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'ocn_mld', &
                   'units', c = 'm')
CALL channel_halt(substr, status)

IF (p_parallel_io) &
   WRITE(*, *) '... ocn_mld added to channel ', TRIM(modstr) // '_inp'

! add channel object: fraction of organic matter in soil
CALL new_channel_object(status, TRIM(modstr) // '_inp', 'foms', p2 = foms) 
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'foms', &
                  'long_name',  c = 'fraction of organic matter in soil')
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'foms', &
                   'units', c = 'kg(om) kg(solid)-1')
CALL channel_halt(substr, status)

IF (p_parallel_io) &
   WRITE(*, *) '... foms added to channel ', TRIM(modstr) // '_inp'

! add channel object: dry bulk density of soil
CALL new_channel_object(status, TRIM(modstr) // '_inp', 'rhos', p2 = rhos)
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'rhos', &
                  'long_name',  c = 'dry bulk density of soil')
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'rhos', &
                   'units', c = 'kg(solid) m-3')
CALL channel_halt(substr, status)

IF (p_parallel_io) &
   WRITE(*, *) '... rhos added to channel ', TRIM(modstr) // '_inp'

! add channel object: land snow fraction at t-dt
CALL new_channel_object(status, TRIM(modstr) // '_inp', 'cvs_old', &
                        p2 = cvs_old)
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'cvs_old', &
                  'long_name',  c = 'old snow fraction')
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'cvs_old', &
                   'units', c = 'fraction')
CALL channel_halt(substr, status)

IF (p_parallel_io) &
   WRITE(*, *) '... cvs_old added to channel ', TRIM(modstr) // '_inp'

!! mz_jw_20161019+
! add channel object: NO2 photolysis rate
CALL new_channel_object(status, TRIM(modstr) // '_inp', 'J_NO2', &
                        reprid = GP_3D_MID, p3 = J_out_NO2)
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'J_NO2', &
                  'long_name',  c = 'J_NO2 values')
CALL channel_halt(substr, status)

CALL new_attribute(status, TRIM(modstr) // '_inp', 'J_NO2', &
                   'units', c = '1/s')
CALL channel_halt(substr, status)

IF (p_parallel_io) &
    WRITE(*, *) '... J_NO2 added to channel ', TRIM(modstr) // '_inp'
!! mz_jw_20161019-

CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

END SUBROUTINE svoc_init_memory

! *************************************************************************

SUBROUTINE svoc_init_coupling

!!!  allocate coupling variables

USE messy_main_mpi_bi,           ONLY: p_parallel_io
USE messy_main_tracer_mem_bi,    ONLY: ti_gp, GPTRSTR
USE messy_main_channel_error_bi, ONLY: channel_halt
USE messy_main_tracer,           ONLY: R_MOLARMASS, get_tracer
USE messy_main_channel,          ONLY: get_channel_object, get_channel_info, &
                                       get_channel_object_info, STRLEN_CHANNEL, &
                                       STRLEN_OBJECT
USE messy_main_data_bi,          ONLY: basemodstr=>modstr

IMPLICIT NONE

CHARACTER(LEN = *), PARAMETER   :: substr = 'svoc_init_coupling'
CHARACTER(LEN = STRLEN_CHANNEL) :: airsea_channel = 'airsea_gp'
CHARACTER(LEN = STRLEN_CHANNEL) :: aermod_channel = ''
CHARACTER(LEN = STRLEN_CHANNEL) :: drydep_channel = 'ddep_gp'
CHARACTER(LEN = STRLEN_CHANNEL) :: scav_channel = 'scav_gp'
!! mz_jw_20160920+
CHARACTER(LEN = STRLEN_CHANNEL) :: jval_channel = 'jval_gp'
!! mz_jw_20160920-
CHARACTER(LEN = STRLEN_OBJECT)  :: name = ''
CHARACTER(LEN = 2) :: txt = ''
INTEGER :: status, reprid
INTEGER :: idt_dummy = 0
INTEGER :: idt, jt, jm, js

INTRINSIC :: SIZE

CALL start_message_bi(modstr, 'COUPLING INITIALIZATION', substr)

CALL info_bi('use molar mass in tracer meta information', substr)

DO jt = 1, nsvoc
   IF (l_svoc_gp(jt)) THEN

      IF (ti_gp(idt_svocg(jt))%tp%meta%cask_r(R_MOLARMASS) .LT. 1._dp) THEN
         CALL info_bi('problem with molar mass in tracer definition!', &
                       substr)
         CALL info_bi('use molar mass defined in namelist', substr)
      ELSE
         MOLMASS(jt) = ti_gp(idt_svocg(jt))%tp%meta%cask_r(R_MOLARMASS)
      END IF
   END IF

   IF (MOLMASS(jt) .LT. 1._dp) THEN
      CALL error_bi('molar mass below limit (<1.0)', substr)
   END IF
END DO

! import sea spray aerosol fluxes (kg m-2 s-1)
CALL info_bi('IMPORT SS FLUXES', substr)

CALL get_channel_info(status, TRIM(imp_ss_as_flux%cha))
IF (status /= 0) THEN
   IF (p_parallel_io) WRITE(*, *) 'channel soil not present!'
   CALL error_bi('no channel soil defined in namelist!',substr)

ELSE
   CALL get_channel_object_info(status, &
                                cname = TRIM(imp_ss_as_flux%cha), &
                                oname = TRIM(imp_ss_as_flux%obj), &
                                reprid = reprid)

   IF (status /= 0) THEN
      IF (p_parallel_io) WRITE(*, *) 'object ss_as_flux not present!'
      CALL error_bi('no object ss_as_flux defined in namelist!', substr)
   ELSE
      CALL get_channel_object(status, &
                              TRIM(imp_ss_as_flux%cha), &
                              TRIM(imp_ss_as_flux%obj), &
                              p2 = in_ss_as_flux)
      CALL channel_halt(substr, status)
   ENDIF
END IF

! import organic matter content in soil [kg (OM) m-3]
CALL info_bi('IMPORT OM SOIL CONTENT', substr)

CALL get_channel_info(status, TRIM(imp_om_soil%cha))
IF (status /= 0) THEN
   IF (p_parallel_io) WRITE(*, *) 'channel soil not present!'
   CALL error_bi('no channel soil defined in namelist!',substr)

ELSE
   CALL get_channel_object_info(status, &
                                cname = TRIM(imp_om_soil%cha), &
                                oname = TRIM(imp_om_soil%obj), &
                                reprid = reprid)

   IF (status /= 0) THEN
      IF (p_parallel_io) WRITE(*, *) 'object foms not present!'
      CALL error_bi('no object foms defined in namelist!', substr)
   ELSE
      CALL get_channel_object(status, &
                              TRIM(imp_om_soil%cha), &
                              TRIM(imp_om_soil%obj), &
                              p2 = in_foms)
      CALL channel_halt(substr, status)
   ENDIF
END IF

! import dry bulk density of soil [kg solid m-3]
CALL info_bi('IMPORT DRY BULK DENS', substr)

CALL get_channel_info(status, TRIM(imp_rho_soil%cha))
IF (status /= 0) THEN
   IF (p_parallel_io) WRITE(*, *) 'channel soil not present!'
   CALL error_bi('no channel soil defined in namelist!',substr)
ELSE

   CALL get_channel_object_info(status, &
                                cname = TRIM(imp_rho_soil%cha), &
                                oname = TRIM(imp_rho_soil%obj), &
                                reprid = reprid)

   IF (status /= 0) THEN
      IF (p_parallel_io) WRITE(*, *) 'object rhos not present!'
      CALL error_bi('no object rhos defined in namelist!', substr)
   ELSE
      CALL get_channel_object(status, &
                              TRIM(imp_rho_soil%cha), &
                              TRIM(imp_rho_soil%obj), &
                              p2 = in_rhos)
      CALL channel_halt(substr, status)
   ENDIF
END IF

! import mixed layer depth for fake ocean [m]
CALL info_bi('IMPORT MLD OCEAN', substr)

CALL get_channel_info(status, TRIM(imp_mld_oce%cha))
IF (status /= 0) THEN
   IF (p_parallel_io) WRITE(*, *) 'channel oce not present!'
   CALL error_bi('no channel oce defined in namelist!',substr)
ELSE

   CALL get_channel_object_info(status, &
                                cname = TRIM(imp_mld_oce%cha), &
                                oname = TRIM(imp_mld_oce%obj), &
                                reprid = reprid)

   IF (status /= 0) THEN
      IF (p_parallel_io) WRITE(*, *) 'object ocn_mld not present!'
      CALL error_bi('no object ocn_mld defined in namelist!', substr)
   ELSE
      CALL get_channel_object(status, &
                              TRIM(imp_mld_oce%cha), &
                              TRIM(imp_mld_oce%obj), &
                              p2 = in_omld)
      CALL channel_halt(substr, status)
   ENDIF
END IF
! import 10m wind
   CALL get_channel_object(status, TRIM(basemodstr), 'wind10', p2=wind10_2d)
   CALL channel_halt(substr, status)

! import svoc emissions [kg(trac) m-2 s-1]
CALL info_bi("IMPORT SVOC EMISSIONS", substr)

ALLOCATE(in_emiss(nsvoc))

DO jt = 1, nsvoc
   IF (p_parallel_io) WRITE(*, *) 'SVOC NAME: ', TRIM(SVOC_NAME(jt))

   IF (TRIM(EMISS_IN(jt)%cha) == '') THEN
      IF (p_parallel_io) &
         WRITE(*, *) 'no channel present for ', TRIM(SVOC_NAME(jt))

      CALL info_bi('use zero application-emission', substr)
   ELSE
      CALL get_channel_info(status, TRIM(EMISS_IN(jt)%cha))

      IF (status /= 0) THEN
         IF (p_parallel_io) &
            WRITE(*, *) 'no channel emiss defined in namelist for ', &
                        TRIM(SVOC_NAME(jt))

         CALL info_bi('use zero application-emission', substr)

      ELSE
         CALL get_channel_object_info(status, &
                                      cname = TRIM(EMISS_IN(jt)%cha), &
                                      oname = TRIM(EMISS_IN(jt)%obj), &
                                      reprid = reprid)
  
         IF (status /= 0) THEN
            IF (p_parallel_io) &
               WRITE(*, *) 'object ', TRIM(EMISS_IN(jt)%obj), &
                ' not present in channel for ', TRIM(SVOC_NAME(jt))

            CALL info_bi('use zero application-emission', substr)
 
            in_emiss(jt)%ptr = 0._dp
         ELSE
            CALL get_channel_object(status, &
                                    TRIM(EMISS_IN(jt)%cha), &
                                    TRIM(EMISS_IN(jt)%obj), &
                                    p2 = in_emiss(jt)%ptr)
            CALL channel_halt(substr, status)

            l_emiss(jt) = .TRUE.
         END IF
      END IF
   END IF
END DO
      
! import standard deviation of aerosol distribution
! by mode from aerosol sm (m7 or gmxe)
CALL info_bi('IMPORT SIGMA FROM AEROSOL SM', substr)

aermod_channel = TRIM(aermod_str) // '_' // GPTRSTR

CALL get_channel_info(status, TRIM(aermod_channel))
IF (status /= 0) THEN
   IF (p_parallel_io) &
      WRITE(*, *) 'Aerosol module is not running! (no coupling)'
ELSE
   CALL get_channel_object_info(status, &
                                cname = TRIM(aermod_channel), &
                                oname = 'sigma', &
                                reprid = reprid)

   IF (status /= 0) THEN
      IF (p_parallel_io) &
         WRITE(*, *) 'object sigma not present in channel', aermod_channel
   ELSE
      CALL get_channel_object(status, &
                              TRIM(aermod_channel), &
                              'sigma', &
                              p1 = sigma)
      CALL channel_halt(substr, status)
   END IF
END IF

IF (nmod /= SIZE(sigma, 1)) THEN
   CALL error_bi('nmod /= nmodes defined in aerosol module', substr)
END IF

! import aerosol wet radius from aerosol sm (m7 or gmxe) [m]
CALL info_bi('IMPORT AEROSOL WET RADIUS', substr)

CALL get_channel_info(status, TRIM(aermod_channel))
IF (status /= 0) THEN
   IF (p_parallel_io) &
      WRITE(*, *) 'Aerosol module is not running! (no coupling)'
ELSE
   CALL get_channel_object_info(status, &
                                cname = TRIM(aermod_channel), &
                                oname = 'wetradius', &
                                reprid = reprid)

   IF (status /= 0) THEN
      IF (p_parallel_io) &
         WRITE(*, *) 'object wetradius not present in channel', &
                     aermod_channel
   ELSE
      CALL get_channel_object(status, &
                              TRIM(aermod_channel), &
                              'wetradius', &
                              p4 = rwet_aer)
      CALL channel_halt(substr, status)
   END IF
END IF

! if gmxe is used, import aerosol number [cm-3]
IF (TRIM(aermod_str) == 'gmxe') THEN
   CALL info_bi('IMPORT AEROSOL NUMBER', substr)

   CALL get_channel_info(status, TRIM(aermod_channel))
   IF (status /= 0) THEN
      IF (p_parallel_io) &
         WRITE(*, *) 'GMXe is not running! (no coupling)'
   ELSE
      CALL get_channel_object_info(status, &
                                   cname = TRIM(aermod_channel), &
                                   oname = 'anumber', &
                                   reprid = reprid)

      IF (status /= 0) THEN
         IF (p_parallel_io) &
            WRITE(*, *) 'object aernum not present in gmxe_gp channel'
      ELSE
         CALL get_channel_object(status, &
                                 TRIM(aermod_channel), &
                                 'anumber', &
                                 p4 = aernl_cpl)
         CALL channel_halt(substr, status)
      END IF
   END IF
END IF

! request ozone tracer
CALL get_tracer(status, GPTRSTR, 'O3', idx = idt_o3)
IF (status /= 0) CALL error_bi('Tracer O3 not found', substr)

! request aerosol tracers
DO jm = 1, nmod
   CALL get_tracer(status, GPTRSTR, TRIM(bname_num(jm)), &
                   subname = TRIM(cmode(jm)), idx = idt)
   
   IF (status /= 0) THEN
      idt_nnum(jm) = idt_dummy
   ELSE 
      idt_nnum(jm) = idt
   END IF
END DO

idt_maer = idt_dummy

DO jt = 1, naer
   CALL get_tracer(status, GPTRSTR, TRIM(bname_aer(jt)), &
                   subname = TRIM(sname_aer(jt)), idx = idt)

   IF (status == 0) idt_maer(jt) = idt
END DO

! import dry deposition velocity from ddep sm
CALL info_bi('IMPORT DRY DEPOSITION VELOCITY', substr)

ALLOCATE(vddepaer(nmod)) ! aerosols [m s-1]
ALLOCATE(vddepgas(nsvoc)) ! gaseous svocs [cm s-1]

CALL get_channel_info(status, TRIM(drydep_channel))
IF (status /= 0) THEN
   IF (p_parallel_io) &
      WRITE(*, *) 'DDEP is not running! (no coupling)'
ELSE
   DO jm = 1, nmod
      WRITE(txt, '(i2.2)') jm

      name = TRIM(aermod_str) // '_v_mass_m' // TRIM(txt)

      CALL get_channel_object_info(status, &
                                   cname = TRIM(drydep_channel), &
                                   oname = name, &
                                   reprid = reprid)

      IF (status /= 0) THEN
         IF (p_parallel_io) &
            WRITE(*, *) 'object ' // TRIM(aermod_str) // '_v_mass_m' // &
                    TRIM(txt) // 'not present in channel!'
      ELSE
         CALL get_channel_object(status, &
                                 TRIM(drydep_channel), &
                                 name, &
                                 p2 = vddepaer(jm)%ptr)
         CALL channel_halt(substr, status)
      END IF
   END DO

   DO jt = 1, nsvoc
      CALL get_channel_object_info(status, &
                                   cname = TRIM(drydep_channel), &
                                   oname = 'Vd_' // TRIM(SVOC_NAME(jt)), &
                                   reprid = reprid)

      IF (status /= 0) THEN
         IF (p_parallel_io) & 
            WRITE(*, *) 'object Vd_' // TRIM(SVOC_NAME(jt)) // &
                        ' not present in channel'
      ELSE
         CALL get_channel_object(status, &
                                 TRIM(drydep_channel), &
                                 'Vd_' // TRIM(SVOC_NAME(jt)), &
                                  p2 = vddepgas(jt)%ptr)
         CALL channel_halt(substr, status)
      END IF
   END DO
END IF

ALLOCATE(flux_airsea(nsvoc))

IF (l_svocvola) THEN
   ! import air-sea flux from airsea sm [mol m-2 s-1]
   CALL info_bi('IMPORT AIRSEA FLUX', substr)

   CALL get_channel_info(status, TRIM(airsea_channel))
   IF (status /= 0) THEN
      IF (p_parallel_io) &
         WRITE(*, *) 'AIRSEA is not running! (no coupling)'
   ELSE
      DO jt = 1, nsvoc
         CALL get_channel_object_info(status, &
                                      cname = TRIM(airsea_channel), &
                                      oname = 'flux_' // TRIM(SVOC_NAME(jt)), &
                                      reprid = reprid)

         IF (status /= 0) THEN
            IF (p_parallel_io) &
               WRITE(*, *) 'object flux_' // TRIM(SVOC_NAME(jt)) // &
                           ' not present in channel', TRIM(airsea_channel) 
         ELSE
            CALL get_channel_object(status, &
                                    TRIM(airsea_channel), &
                                    'flux_' // TRIM(SVOC_NAME(jt)), &
                                    p2 = flux_airsea(jt)%ptr)
            CALL channel_halt(substr, status)
         END IF
      END DO
   END IF

END IF

! import wet deposition fluxes from scav for gaseous svoc [molecules m-2 s-1]
CALL info_bi('IMPORT DEPOSITION FLUX FROM STRATIFORM RAIN', substr)

ALLOCATE(wetflx_ls(nsvoc*2))  !!! gas and liquid phases

CALL get_channel_info(status, TRIM(scav_channel))
IF (status /= 0) THEN
   IF (p_parallel_io) & 
      WRITE(*, *) 'SCAV is not running! (no coupling)'
ELSE
   DO jt = 1, nsvoc
      CALL get_channel_object_info(status, &
                          cname = TRIM(scav_channel), &
                          oname = 'wetflx_ls_' // TRIM(SVOC_NAME(jt)), &
                          reprid = reprid)

      IF (status /= 0) THEN
         IF (p_parallel_io) &
            WRITE(*, *) 'object wetflx_ls_' // TRIM(SVOC_NAME(jt)) // &
                        ' not present in channel ', TRIM(scav_channel)
      ELSE
         CALL get_channel_object(status, &
                          TRIM(scav_channel), &
                          'wetflx_ls_' // TRIM(SVOC_NAME(jt)), &
                          p2 = wetflx_ls(jt)%ptr)
         CALL channel_halt(substr, status)
      END IF
   END DO

   DO jt = 1, nsvoc
      CALL get_channel_object_info(status, &
                          cname = TRIM(scav_channel), &
                          oname = 'wetflx_ls_' // TRIM(SVOC_NAME(jt)) // '_l', &
                          reprid = reprid)
      
      IF (status /= 0) THEN
         IF (p_parallel_io) &
            WRITE(*, *) 'object wetflx_ls_' // TRIM(SVOC_NAME(jt)) // '_l' // &
                        ' not present in channel ', TRIM(scav_channel)
      ELSE
         CALL get_channel_object(status, &
                          TRIM(scav_channel), &
                          'wetflx_ls_' // TRIM(SVOC_NAME(jt)) // '_l', &
                          p2 = wetflx_ls(jt+nsvoc)%ptr)
         CALL channel_halt(substr, status)
      END IF
   END DO
END IF

CALL info_bi('IMPORT DEPOSITION FLUX FROM CONVECTIVE RAIN', substr)

ALLOCATE(wetflx_cv(nsvoc*2))

CALL get_channel_info(status, TRIM(scav_channel))
IF (status /= 0) THEN
   IF (p_parallel_io) &
      WRITE(*, *) 'SCAV is not running! (no coupling)'
ELSE
   DO jt = 1, nsvoc
      CALL get_channel_object_info(status, &
                          cname = TRIM(scav_channel), &
                          oname = 'wetflx_cv_' // TRIM(SVOC_NAME(jt)), &
                          reprid = reprid)

      IF (status /= 0) THEN
         IF (p_parallel_io) &
            WRITE(*, *) 'object wetflx_cv_' // TRIM(SVOC_NAME(jt)) // &
                        ' not present in channel ', TRIM(scav_channel)
      ELSE
         CALL get_channel_object(status, &
                          TRIM(scav_channel), &
                          'wetflx_cv_' // TRIM(SVOC_NAME(jt)), &
                          p2 = wetflx_cv(jt)%ptr)
         CALL channel_halt(substr, status)
      END IF
   END DO

   DO jt = 1, nsvoc
      CALL get_channel_object_info(status, &
                          cname = TRIM(scav_channel), &
                          oname = 'wetflx_cv_' // TRIM(SVOC_NAME(jt)) // '_l', &
                          reprid = reprid)

      IF (status /= 0) THEN
         IF (p_parallel_io) &
            WRITE(*, *) 'object wetflx_cv_' // TRIM(SVOC_NAME(jt)) // '_l' // &
                        ' not present in channel ', TRIM(scav_channel)
      ELSE
         CALL get_channel_object(status, &
                          TRIM(scav_channel), &
                          'wetflx_cv_' // TRIM(SVOC_NAME(jt)) // '_l', &
                          p2 = wetflx_cv(jt+nsvoc)%ptr)
         CALL channel_halt(substr, status)
      END IF
   END DO
END IF

! import scavenging tendencies of black carbon from scav [mol mol s-1]
CALL info_bi('IMPORT SCAVENGING TENDENCY OF BLACK CARBON', substr)

ALLOCATE(ttescav_bc(4))  !!! only import for BC_ki, BC_ks, BC_as, and BC_cs

js = 0

CALL get_channel_info(status, TRIM(scav_channel))
IF (status /= 0) THEN
   IF (p_parallel_io) &
      WRITE(*, *) 'SCAV is not running! (no coupling)'

ELSE
   DO jt = 1, naer
      IF (TRIM(bname_aer(jt)) == 'BC') THEN
         js = js + 1

         CALL get_channel_object_info(status, &
                          cname = TRIM(scav_channel), &
                          oname = 'BC_' // TRIM(sname_aer(jt)) // '_scte', &
                          reprid = reprid)

         IF (status /= 0) THEN
            IF (p_parallel_io) &
               WRITE(*, *) 'object BC_' // TRIM(sname_aer(jt)) // '_scte' // &
                           ' not present in channel ', TRIM(scav_channel)
         ELSE
            CALL get_channel_object(status, &
                             TRIM(scav_channel), &
                             'BC_' // TRIM(sname_aer(jt)) // '_scte', &
                             p3 = ttescav_bc(js)%ptr)
            CALL channel_halt(substr, status)
         END IF
      END IF
   END DO
END IF

! import tracer field from cvtrans [mol(trac) mol(air)]
CALL info_bi('IMPORT TRACER FIELD FROM CVTRANS', substr)

CALL get_channel_info(status, 'cvtrans')
cvtrans_use = (status == 0)

IF (status /= 0) THEN
   IF (p_parallel_io) &
      WRITE(*, *) 'CVTRANS is not running! (no coupling'

ELSE
   CALL get_channel_object_info(status, &
                                cname = 'cvtrans', oname = 'trac_field', &
                                reprid = reprid)

   IF (status /= 0) THEN
      IF (p_parallel_io) &
         WRITE(*, *) 'object trac_field not found in channel cvtrans'
   ELSE
      CALL get_channel_object(status, &
                              'cvtrans', 'trac_field', &
                              p4 = trac_field)
      CALL channel_halt(substr, status)
   END IF

   CALL get_channel_object_info(status, &
                                cname = 'cvtrans', oname = 'column', &
                                reprid = reprid)
 
   IF (status /= 0) THEN
      IF (p_parallel_io) &
         WRITE(*, *) 'object column not found in channel cvtrans'
   ELSE
      CALL get_channel_object(status, &
                              'cvtrans', 'column', &
                              p0 = col)
      CALL channel_halt(substr, status)

      IF (col > 0._dp) lcolumn = .TRUE.
   END IF
END IF

!! mz_jw_20160920+
! import NO2 photolysis rates [s-1]
CALL info_bi('IMPORT PHOTOLYSIS RATE OF NO2 (J_NO2)', substr)

CALL get_channel_object_info(status, &
                             cname = TRIM(jval_channel), &
                             oname = 'J_NO2', &
                             reprid = reprid)

IF (status /= 0) THEN
   CALL error_bi('J_NO2 not found', substr)
ELSE
   CALL get_channel_object(status, &
                           TRIM(jval_channel), &
                           'J_NO2', p3 = J_in_NO2)

   CALL channel_halt(substr, status)
END IF
!! mz_jw_20160920-

CALL end_message_bi(modstr, 'COUPLING INITIALIZATION', substr)

END SUBROUTINE svoc_init_coupling

! *************************************************************************

SUBROUTINE svoc_vdiff(iflag)

!!! ...

USE messy_main_mpi_bi,          ONLY: p_parallel_io
USE messy_main_grid_def_mem_bi, ONLY: nproma, kproma, nlev, jrow
USE messy_main_grid_def_bi,     ONLY: deltaz
USE messy_main_data_bi,       ONLY: tslm1, & ! surface temp o. land (t-dt)
                                    !wind10_2d, & ! wind speed at 10 m
                                    forest, & ! forest coverage
                                    vgrat, & ! vegetation fraction
                                    sn, & ! snow depth on soil
                                    gld, & ! glacier depth (incl. snow)
                                    cvs, & ! fractional snow cover
                                    cvw, & ! wet skin fraction
                                    ws, & ! surface soil wetness 
                                    wsmx, & ! field capacity of soil
                                    tsw, & ! surface temp o. water
                                    press_3d, & ! full-level pressure
                                    pressi_3d, & ! interface pressure
                                    tm1_3d, & ! temperature at t-dt
                                    tte_3d, & ! temp tendencies [K s-1]
                                    qm1_3d, & ! specific humidity at t-dt
                                    qte_3d, & ! spc hum tendencies
                                    loland_2d, & ! logical mask land 
                                    loglac_2d !!$, & ! logical mask glacier
#if defined(ECHAM5)
USE messy_main_data_bi,       ONLY: pxtems ! surface emission of tracer
                                    ! [mol(trac) mol(air)-1 kg(air) m-2 s-1]
#endif

!!$#ifndef MESSYTENDENCY
USE messy_main_tracer_mem_bi, ONLY: qxtm1, & ! vmr at t-dt
                                    qxtte, & ! vmr tendencies
                                    ti_gp
!!$#else
!!$USE messy_main_tracer_mem_bi, ONLY: ti_gp
!!$#endif

USE messy_main_timer,         ONLY: delta_time, lresume, lstart, &
                                    time_step_len

USE messy_main_tracer,           ONLY: R_MOLARMASS
!!$USE mo_time_control,             ONLY: lfirst_day ! op_pj_20190326
USE messy_main_constants_mem,    ONLY: R_gas, pi, M_air, g, &
                                       rd, & ! 1000. * R_gas/M_air (J/K/kg)
                                       vtmpc1 ! M_air / M_H2O - 1.

IMPLICIT NONE

CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_vdiff'
INTEGER, INTENT(IN) :: iflag
REAL(dp), PARAMETER :: kg2ug = 1.e+9_dp
REAL(dp), PARAMETER :: m2cm = 1.e2_dp
REAL(dp), PARAMETER :: g2kg = 1.e-3_dp
!!$#ifdef MESSYTENDENCY
REAL(dp) :: zxtte(kproma, nlev)
!!$#endif
REAL(dp), DIMENSION(:, :), POINTER :: zxtems => NULL()
REAL(dp), DIMENSION(:, :), POINTER :: vddep_p => NULL()
REAL(dp), DIMENSION(nproma, nlev, nmod) :: rarea, rarea_bc
REAL(dp), DIMENSION(nproma, nlev) :: rtot_area, rtot_area_bc
REAL(dp), DIMENSION(nproma) :: rho_air, zdz, zdp
REAL(dp), DIMENSION(nproma) :: oceflux
REAL(dp), DIMENSION(nproma) :: xtp1_g, emiss_g
REAL(dp), DIMENSION(nproma, nmod) :: zvdrydep, emiss_m, xtp1_m
REAL(dp) :: spc_hum, temp, mol2mass
REAL(dp) :: cmr2ras, zaernl, zaerml
REAL(dp) :: ztotmass, zbcmass
REAL(dp) :: om_ss_frac
CHARACTER(LEN = 2) :: txt = ''
INTEGER  :: jl, jk, jt, jr, jm, jn
INTEGER  :: idt, idt2

CALL start_message_bi(modstr, 'SVOC VDIFF ROUTINES', substr)

#ifdef ECHAM5
   zxtems => pxtems(:, 1, :,jrow)
#endif

SELECT CASE (iflag)
CASE (1) ! before calling airsea, ddep, and m7/gmxe
   !!!
   !!! initialize pointers
   !!!

   DO jt = 1, nsvoc
      burd_p => burd(jt)%ptr
      vola_p => vola(jt)%ptr
      degr_p => degr(jt)%ptr
      depo_p => depo(jt)%ptr
      kp_p => kp(jt)%ptr
      part_p => part(jt)%ptr
      cwat_p => cwat(jt)%ptr

      IF (lstart) THEN
         burd_p(1: kproma, :,jrow) = 0._dp
         vola_p(1: kproma, :,jrow) = 0._dp
         degr_p(1: kproma, :,jrow) = 0._dp
         depo_p(1: kproma, :,jrow) = 0._dp
         kp_p(1: kproma, :,jrow) = 0._dp
         part_p(1: kproma, :,jrow) = 0._dp
         cwat_p(1: kproma ,jrow) = 0._dp
      END IF

      IF (l_mode_partition) THEN
         DO jm = 1, nmod
            kp_p => kp_m(jt, jm)%ptr
            part_p => part_m(jt, jm)%ptr
 
            IF (lstart) THEN
               kp_p(1: kproma, :,jrow) = 0._dp
               part_p(1: kproma, :,jrow) = 0._dp
            END IF
         END DO
      END IF
   END DO

   !!!
   !!! application-emission (primary emissions)
   !!!

   DO jt = 1, nsvoc
      IF (l_emiss(jt)) THEN
         appl_p => appl(jt)%ptr
         burd_p => burd(jt)%ptr

         idt = idt_svocg(jt)  ! gas

         ! create factor MW/M_air to convert volume mixing ratio
         ! [mol mol-1] into mass mixing ratio [kg kg-1]
         mol2mass = ti_gp(idt)%tp%meta%cask_r(R_MOLARMASS) / M_air

         ! convert emission unit from [kg(trac) m-2 s-1] to
         ! [mol(trac) mol(air) kg(air) m-2 s-1] as in pxtems
         appl_p = in_emiss(jt)%ptr / mol2mass

         DO jl = 1, kproma
            ! update burden in soil and vegetation (kg m-2)
            IF (loland_2d(jl ,jrow)) THEN
               jk = nlev ! soil
               burd_p(jl, jk,jrow) = burd_p(jl, jk,jrow) + &
                        in_emiss(jt)%ptr(jl ,jrow) * (1._dp - RSPRAY(jt)) * &
                        (1._dp - RLOSS(jt)) * delta_time

               jk = nlev - 1 ! vegetation
               burd_p(jl, jk,jrow) = burd_p(jl, jk,jrow) + &
                       in_emiss(jt)%ptr(jl ,jrow) *  RSPRAY(jt) * &
                       (1._dp - RLOSS(jt)) * delta_time
            END IF

            ! update emissions
            zxtems(jl, idt) = appl_p(jl ,jrow) * RLOSS(jt) * &
                              (1._dp - REMISP(jt))

            IF (.NOT. l_mode_partition) THEN
               idt2 = idt_svocp(jt) ! particle (bulk)
               zxtems(jl, idt2) = appl_p(jl ,jrow) * RLOSS(jt) * REMISP(jt)
            ELSE
               DO jm = 1, nmod
                  idt2 = idt_svocm(jt, jm)

                  zxtems(jl, idt2) = appl_p(jl ,jrow) * RLOSS(jt) * &
                                     (REMISP(jt) * REMIS_MOD(jt, jm))
               END DO
            END IF
         END DO
      END IF
   END DO

   !!!
   !!! distribute sea-spray (as mode) fluxes into sea-spray aerosols and OM
   !!!
   !!! organic matter fraction in sea spray is based on a wind-dependent empirical
   !!! relationship (Gantt et al., ACP, 2011)
   !!!
   !!! avoid calculation at first time step due to the absence of winds
   !!!
  
   IF (.NOT. lstart) THEN
      DO jl = 1, kproma
         om_ss_frac = 0._dp

         IF (wind10_2d(jl ,jrow) > 0._dp) THEN
            om_ss_frac = 0.5_dp * &
                  (0.78_dp / (1._dp + 0.03_dp * &
                  EXP(0.48_dp * wind10_2d(jl ,jrow))) + &
                  0.24_dp / (1._dp + 0.05_dp * &
                  EXP(0.38_dp * wind10_2d(jl ,jrow))))
         END IF

         om_ss_frac = MIN(1._dp, om_ss_frac)

         ! OM in sea-spray fluxes; added into OC_as tracer
         DO jr = 1, naer
            IF (TRIM(bname_aer(jr)) == 'OC' .AND. &
               TRIM(sname_aer(jr)) == 'as') &
               idt = idt_maer(jr)
         END DO

         zxtems(jl, idt) = zxtems(jl, idt) + &
               om_ss_frac * in_ss_as_flux(jl ,jrow) * &
               M_air / ti_gp(idt)%tp%meta%cask_r(R_MOLARMASS)

         ! update sea-spray aerosol flux
         in_ss_as_flux(jl ,jrow) = (1 - om_ss_frac) * in_ss_as_flux(jl ,jrow)
      END DO
   END IF

   !!!
   !!! volatilization from soil, veg., snow, glacier (secondary emissions)
   !!!

   IF (l_svocvola) THEN
      foms = in_foms
      rhos = in_rhos

      DO jt = 1, nsvoc
         burd_p => burd(jt)%ptr
         vola_p => vola(jt)%ptr

         idt = idt_svocg(jt)  ! gas

         CALL svoc_volatilizations(kproma, nlev, delta_time, &
                                   param_soilv, &
                                   MOLMASS(jt), &
                                   RWSOL(jt), RHSOL(jt), &
                                   RVAPP(jt), RHVAP(jt), &
                                   RLOGKOW(jt), RHSUB(jt), &
                                   KAHENRY(jt), KBHENRY(jt), &
                                   loland_2d(1: kproma ,jrow), &
                                   loglac_2d(1: kproma ,jrow), &
                                   forest(1: kproma ,jrow), &
                                   vgrat(1: kproma ,jrow), &
                                   sn(1: kproma ,jrow), &
                                   gld(1: kproma ,jrow), &
                                   cvs(1: kproma ,jrow), &
                                   ws(1: kproma ,jrow), &
                                   wsmx(1: kproma ,jrow), &
                                   tslm1(1: kproma ,jrow), &
                                   foms(1: kproma ,jrow), &
                                   rhos(1: kproma ,jrow), &
                                   zxtems(1: kproma, idt), &
                                   cvs_old(1: kproma ,jrow), &
                                   burd_p(1: kproma, 1:nlev,jrow), &
                                   vola_p(1: kproma, 1:nlev,jrow))
      END DO
   END IF

   !!!
   !!! biotic degradation in soil, vegetation and ocean
   !!!

   ocn_mld = in_omld

   DO jt = 1, nsvoc
      burd_p => burd(jt)%ptr
      degr_p => degr(jt)%ptr

      CALL svoc_degradation(kproma, nlev, delta_time, &
!!$                         lfirst_day, & ! op_pj_20190326
                            lstart, & ! op_pj_20190326
                            RKSOIL(jt), RKOCEAN(jt), &
                            loland_2d(1: kproma ,jrow), &
                            tslm1(1: kproma ,jrow), &
                            ocn_mld(1: kproma ,jrow), &
                            burd_p(1: kproma, 1:nlev,jrow), &
                            degr_p(1: kproma, 1:nlev,jrow))
   END DO

   !!!
   !!! svoc concentration in ocean [kmol m-3]
   !!! required by airsea
   !!!

   DO jt = 1, nsvoc
      cwat_p => cwat(jt)%ptr
      burd_p => burd(jt)%ptr

      jk = nlev - 2  ! ocean

      DO jl = 1, kproma
         IF (.NOT. loland_2d(jl ,jrow)) THEN
            IF (ocn_mld(jl ,jrow) > 0._dp) THEN
               !!! kmol mol-1 * kg(trac) m-2 / &
               !!! (g(trac) mol(trac)-1 * kg g-1 * m) = kmol m-3
               cwat_p(jl ,jrow) = 1.e-3 * burd_p(jl, jk,jrow) / &
                                    (MOLMASS(jt) * g2kg) / ocn_mld(jl ,jrow)
            ELSE
               cwat_p(jl ,jrow) = 0._dp
            END IF
         END IF

      END DO
   END DO

CASE (2) ! after calling airsea, m7/gmxe and ddep
   !!!
   !!! compute intermediate variables
   !!!

   ! delta pressure at lowest level
   zdp(1: kproma) = pressi_3d(1: kproma, nlev+1,jrow) - &
                    pressi_3d(1: kproma, nlev,jrow)

   ! layer thickness [m]
   zdz(1: kproma) = deltaz(1: kproma, nlev,jrow)

   !!!
   !!! update storage, vola, and depo over ocean
   !!! need coupling with airsea sm
   !!!

   DO jt = 1, nsvoc
      jk = nlev - 2  ! ocean
      idt = idt_svocg(jt)

      ! create factor MW/M_air to convert volume mixing ratio
      ! [mol mol-1] into mass mixing ratio [kg kg-1]
      mol2mass = ti_gp(idt)%tp%meta%cask_r(R_MOLARMASS) / M_air

      vola_p => vola(jt)%ptr
      depo_p => depo(jt)%ptr
      flux_airsea_p => flux_airsea(jt)%ptr ! [mol m-2 s-1]

      oceflux(:) = 0._dp

      DO jl = 1, kproma
         IF (.NOT. loland_2d(jl ,jrow)) THEN
            ! convert flux from [mol(trac) m-2 s-1] to 
            ! [mol(trac) mol(air)-1 kg(air) m-2 s-1]

            IF (l_svocvola) &
               oceflux(jl) = flux_airsea_p(jl ,jrow) * M_air * g2kg

            IF (oceflux(jl) > 0._dp) THEN ! volatilization
               ! if all burden are readily volatilized
               IF (burd_p(jl, jk,jrow) / mol2mass .LT. &
                   oceflux(jl) * delta_time) THEN
                  ! correction on zxtems
                  zxtems(jl, idt) = zxtems(jl, idt) - oceflux(jl) + &
                             burd_p(jl, jk,jrow) / (mol2mass * delta_time)

                  ! change oceflux value
                  oceflux(jl) = burd_p(jl, jk,jrow) / (mol2mass * delta_time)
               END IF

               vola_p(jl, jk,jrow) = oceflux(jl) * mol2mass

            ELSE ! dry deposition
               depo_p(jl, nlev,jrow) = -oceflux(jl) * mol2mass
            END IF

            burd_p(jl, jk,jrow) = burd_p(jl, jk,jrow) - &
                               oceflux(jl) * mol2mass * delta_time
         END IF
      END DO  ! kproma
   END DO  ! nsvoc

   !!!
   !!! dry deposition
   !!! need coupling with ddep sm
   !!!

   ! calculate total area of aerosols and black carbon at the lowest level
   rarea = 0._dp
   rarea_bc = 0._dp
   rtot_area_bc = 0._dp

   IF (.NOT. l_mode_partition) THEN
      DO jl = 1, kproma
         jk = nlev

         spc_hum = qm1_3d(jl, jk,jrow) + qte_3d(jl, jk,jrow) * time_step_len
         temp = tm1_3d(jl, jk,jrow) + tte_3d(jl, jk,jrow) * time_step_len

         ! density of moist air [kg m-3]
         ! the factor 1+vtmpc1*spc_hum accounts for moist correction
         rho_air(jl) = press_3d(jl, jk,jrow) / &
                       (temp * rd * (1._dp + vtmpc1 * spc_hum))

         DO jm = 1, nmod
            ztotmass = 0._dp
            zbcmass = 0._dp

            ! convert count median radius to radius of average surface based
            ! on the hatch-choate conversion (see hinds, equation 4.52)
            ! radius of average surface is used because svoc is only adsorbed
            ! on (bc) aerosol surface
            cmr2ras = EXP(1.0_dp * (LOG(sigma(jm))) ** 2)

            ! calculate total number of aerosols
            IF (aermod_str == 'gmxe') THEN
               zaernl =  aernl_cpl(jl, jk,jm,jrow)
            ELSE IF (aermod_str == 'm7') THEN
               idt = idt_nnum(jm)

               ! convert the number from mol-1 to cm-3 using M_air
               ! 1e-6:  factor to transform m-3 to cm-3
               zaernl = rho_air(jl) * 1.e-6_dp * (qxtm1(jl, jk,idt) + &
                        qxtte(jl, jk,idt) * time_step_len) / (M_air * g2kg)
            END IF

            ! calculate total surface area of aerosol per cm3 of air [cm2 cm-3]
            ! rwet should be converted to radius of avg surface using cmr2ras
            IF (rwet_aer(jl, jk,jm,jrow) > 0._dp) THEN
               rarea(jl, jk, jm) = 4._dp * pi * &
                                  ((rwet_aer(jl, jk,jm,jrow) * m2cm * &
                                  cmr2ras) ** 2) * zaernl
            END IF

            ! calculate total aerosol mass at selected modes (ks, as, ki, cs)
            IF (TRIM(cmode(jm)) == 'ki' .OR. TRIM(cmode(jm)) == 'ks' .OR. &
               TRIM(cmode(jm)) == 'as' .OR. TRIM(cmode(jm)) == 'cs') THEN

               txt = TRIM(cmode(jm))

               DO jr = 1, naer
                  idt = idt_maer(jr)
                  zaerml = 0._dp

                  ! create factor MW/M_air to convert volume mixing ratio
                  ! [mol mol-1] into mass mixing ratio [kg kg-1]
                  mol2mass = ti_gp(idt)%tp%meta%cask_r(R_MOLARMASS) / M_air

                  IF (TRIM(sname_aer(jr)) == txt) THEN
                     ! unit of zaerml is ug m-3
                     zaerml = rho_air(jl) * kg2ug * mol2mass * &
                              (qxtm1(jl, jk,idt) + qxtte(jl, jk,idt) * &
                              time_step_len)

                     ! remove negative values
                     zaerml = MAX(zaerml, 0._dp)

                     IF (TRIM(bname_aer(jr)) == 'BC') zbcmass = zaerml
                  END IF

                  ztotmass = ztotmass + zaerml
               END DO
   
               ! calculate area of bc at lowest layer [cm2 cm-3]
               IF (ztotmass > 0._dp) &
                  rarea_bc(jl, jk, jm) = rarea(jl, jk, jm) * zbcmass / ztotmass
            END IF

            rtot_area_bc(jl, jk) = rtot_area_bc(jl, jk) + rarea_bc(jl, jk, jm)
         END DO
      END DO
   END IF  ! .not. l_mode_partition
   
   ! calculate area weighted dry deposition velocity
   zvdrydep = 0._dp ! [cm s-1]

   DO jm = 1, nmod
      vddep_p => vddepaer(jm)%ptr
      
      IF (.NOT. l_mode_partition) THEN
         ! zvd = (vd(bc_ks) * area(bc_ks) + vd(bc_as) * area(bc_as) +
         !       vd(bc_ki) * area(bc_ki) + vd(bc_cs) + area(bc_cs)) / 
         !       totarea(bc)
         DO jl = 1, kproma
            IF (rtot_area_bc(jl, nlev) .GE. 1.e-20_dp) THEN
               zvdrydep(jl, 1) = zvdrydep(jl, 1) + vddep_p(jl ,jrow) * &
                                 m2cm * rarea_bc(jl, nlev, jm) / &
                                 rtot_area_bc(jl, nlev)
            END IF
         END DO

      ELSE
         zvdrydep(1: kproma, jm) = vddep_p(1: kproma ,jrow)
      END IF
   END DO

   DO jt = 1, nsvoc 
      idt = idt_svocg(jt)  ! gas
 
      burd_p => burd(jt)%ptr
      depo_p => depo(jt)%ptr
      vddep_p => vddepgas(jt)%ptr

      jk = nlev

      xtp1_g(1: kproma) = qxtm1(1: kproma, jk,idt) + &
                          qxtte(1: kproma, jk,idt) * time_step_len + &
                          zxtems(1: kproma, idt) / (zdp(1: kproma) / g) * &
                          time_step_len
                          
      ! remove negative values of concentration
      POSTV1: FORALL (jl = 1: kproma, xtp1_g(jl) < 0._dp)
                xtp1_g(jl) = 0._dp
             END FORALL POSTV1

      emiss_g(1: kproma) = zxtems(1: kproma, idt)

      xtp1_m = 0._dp
      emiss_m = 0._dp

      IF (.NOT. l_mode_partition) THEN
         idt2 = idt_svocp(jt) ! as bulk

         xtp1_m(1: kproma, 1) = qxtm1(1: kproma, jk,idt2) + &
                           qxtte(1: kproma, jk,idt2) * time_step_len + &
                           zxtems(1: kproma, idt2) / (zdp(1: kproma) / g) * &
                           time_step_len

         ! remove negative values of concentration
         POSTV2: FORALL (jl = 1: kproma, xtp1_m(jl, 1) < 0._dp)
                   xtp1_m(jl, 1) = 0._dp
                END FORALL POSTV2

         emiss_m(1: kproma, 1) = zxtems(1: kproma, idt2)

         CALL svoc_drydep(kproma, nlev, 1, &
                          delta_time, time_step_len, &
                          param_soilv, &
                          MOLMASS(jt), &
                          KAHENRY(jt), KBHENRY(jt), &
                          loland_2d(1: kproma ,jrow), &
                          loglac_2d(1: kproma ,jrow), &
                          cvs(1: kproma ,jrow), &
                          vgrat(1: kproma ,jrow), &
                          tslm1(1: kproma ,jrow), &
                          vddep_p(1: kproma ,jrow), &
                          zvdrydep(1: kproma, 1), &
                          zdz(1: kproma), &
                          zdp(1: kproma), &
                          xtp1_g(1: kproma), &
                          xtp1_m(1: kproma, 1), &
                          emiss_g(1: kproma), &
                          emiss_m(1: kproma, 1), &
                          burd_p(1: kproma, 1:nlev,jrow), &
                          depo_p(1: kproma, nlev,jrow))

         zxtems(1: kproma, idt) = emiss_g(1: kproma)
         zxtems(1: kproma, idt2) = emiss_m(1: kproma, 1)

      ELSE
         DO jm = 1 , nmod
            idt2 = idt_svocm(jt, jm)
 
            xtp1_m(1: kproma, jm) = qxtm1(1: kproma, jk,idt2) + &
                           qxtte(1: kproma, jk,idt2) * time_step_len + &
                           zxtems(1: kproma, idt2) / (zdp(1: kproma) / g) * &
                           time_step_len

            ! remove negative values of concentration
            POSTV3: FORALL (jl = 1: kproma, xtp1_m(jl, jm) < 0._dp)
                      xtp1_m(jl, jm) = 0._dp
                   END FORALL POSTV3

            emiss_m(1: kproma, jm) = zxtems(1: kproma, idt2)
         END DO
           
         CALL svoc_drydep(kproma, nlev, nmod, &
                          delta_time, time_step_len, &
                          param_soilv, &
                          MOLMASS(jt), &
                          KAHENRY(jt), KBHENRY(jt), &
                          loland_2d(1: kproma ,jrow), &
                          loglac_2d(1: kproma ,jrow), &
                          cvs(1: kproma ,jrow), &
                          vgrat(1: kproma ,jrow), &
                          tslm1(1: kproma ,jrow), &
                          vddep_p(1: kproma ,jrow), &
                          zvdrydep(1: kproma, 1: nmod), &
                          zdz(1: kproma), &
                          zdp(1: kproma), &
                          xtp1_g(1: kproma), &
                          xtp1_m(1: kproma, 1: nmod), &
                          emiss_g(1: kproma), &
                          emiss_m(1: kproma, 1: nmod), &
                          burd_p(1: kproma, 1:nlev,jrow), &
                          depo_p(1: kproma, nlev,jrow))

         zxtems(1: kproma, idt) = emiss_g(1: kproma)

         DO jm = 1, nmod
            idt2 = idt_svocm(jt, jm)
            zxtems(1: kproma, idt2) = emiss_m(1: kproma, jm)
         END DO

      END IF ! l_mode_partition
   END DO ! nsvoc

   DO jl = 1, kproma
      IF (.NOT. loglac_2d(jl ,jrow)) THEN
         cvs_old(jl ,jrow) = cvs(jl ,jrow)
      END IF
   END DO

CASE DEFAULT
   ! should never reach here
END SELECT

CALL end_message_bi(modstr, 'SVOC VDIFF ROUTINES', substr)

END SUBROUTINE svoc_vdiff

! *************************************************************************

SUBROUTINE svoc_convec

!!!

USE messy_main_grid_def_mem_bi, ONLY: nproma, kproma, nlev, jrow
USE messy_main_grid_def_bi,     ONLY: deltaz
USE messy_main_data_bi,       ONLY: tm1_3d, tte_3d, &
                                    qm1_3d, qte_3d, &
                                    press_3d, pressi_3d
USE messy_main_timer,         ONLY: time_step_len
USE messy_main_tracer_mem_bi, ONLY: qxtm1, & ! mole frac at t-dt
                                    qxtte    ! tracer tendencies
USE messy_main_constants_mem, ONLY: g, rd, & ! rd = 1000. * R_gas/M_air (J/K/kg)
                                    vtmpc1 ! M_air / M_H2O - 1.


IMPLICIT NONE

CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_convec'
REAL(dp), DIMENSION(nsvoc, nproma, nlev, nmod) :: ttescav_m
REAL(dp), DIMENSION(nsvoc, nproma, nlev) :: tot_ttescav
REAL(dp), DIMENSION(nproma, nlev) :: wetflx_aer, zdz
REAL(dp), DIMENSION(nproma) :: xtp1_p, xbc1_p, rho_air
REAL(dp) :: spc_hum, temp
REAL(dp) :: tot_flxpos, tot_flxneg, fac
INTEGER :: jl, jk, jt, js, jr, jm, jn
INTEGER :: idt, idt2 

CALL start_message_bi(modstr, 'SVOC CONVEC ROUTINES', substr)

!!!
!!! update tracer tendency due to convective scavenging
!!! need coupling with scav sm
!!!

wetflx_aer_tot = 0._dp
tot_ttescav_cv = 0._dp
tot_ttescav = 0._dp
ttescav_m_cv = 0._dp
ttescav_m = 0._dp

DO jt = 1, nsvoc
   DO js = 1, 4
      ttescav_bc_cv(js, :, :) = 0._dp
   END DO

   wetflx_aer(:, :) = 0._dp
   zdz(:, :) = 0._dp

   ! layer thickness [m]
   zdz(:, 1: nlev) = deltaz(1: kproma, 1:nlev,jrow)

   DO jk = 1, nlev
      xtp1_p(:) = 0._dp
      xbc1_p(:) = 0._dp

      DO jl = 1, kproma
         spc_hum = qm1_3d(jl, jk,jrow) + &
                   qte_3d(jl, jk,jrow) * time_step_len

         temp = tm1_3d(jl, jk,jrow) + &
                tte_3d(jl, jk,jrow) * time_step_len

         ! density of moist air [kg m-3]
         ! the factor 1+vtmpc1*spc_hum accounts for moist correction
         rho_air(jl) = press_3d(jl, jk,jrow) / &
                      (temp * rd * (1._dp + vtmpc1 * spc_hum))

         js = 0
  
         DO jr = 1, naer
            IF (TRIM(bname_aer(jr)) == 'BC') THEN
               js = js + 1
               idt = idt_maer(jr)

               DO jm = 1, nmod
                  IF (TRIM(sname_aer(jr)) == TRIM(cmode(jm))) THEN
                     jn = jm
                     EXIT
                  END IF
               END DO

               ttescav_bc_cv(js, jl, jk) = ttescav_bc(js)%ptr(jl, jk,jrow)

               IF (l_mode_partition) THEN
                  ! bc conc before scav_cv
                  IF (lcolumn) THEN
                     xbc1_p(jl) = trac_field(jl, jk,idt,jrow) - &
                           ttescav_bc_cv(js, jl, jk) * time_step_len
                  ELSE
                     xbc1_p(jl) = qxtm1(jl, jk,idt) + (qxtte(jl, jk,idt) - &
                           ttescav_bc_cv(js, jl, jk)) * time_step_len
                  END IF

                  idt = idt_svocm(jt, jn)

                  IF (lcolumn) THEN
                     xtp1_p(jl) = trac_field(jl, jk,idt,jrow)
                  ELSE
                     xtp1_p(jl) = qxtm1(jl, jk,idt) + &
                                  qxtte(jl, jk,idt) * time_step_len
                  END IF

                  IF (xbc1_p(jl) > 1.e-34_dp .AND. xtp1_p(jl) > 1.e-34_dp) &
                     ttescav_m(jt, jl, jk, jn) = &
                        xtp1_p(jl) * ttescav_bc_cv(js, jl, jk) / xbc1_p(jl)

                  IF (xtp1_p(jl) + ttescav_m(jt, jl, jk, jn) * &
                      time_step_len < 0._dp) ttescav_m(jt, jl, jk, jn) = &
                      -xtp1_p(jl) / time_step_len

               ELSE ! bulk
                  ! bc conc before scav_cv (sum all modes)
                  IF (lcolumn) THEN
                     xbc1_p(jl) = MAX(xbc1_p(jl), 0._dp) + &
                           trac_field(jl, jk,idt,jrow) - &
                           ttescav_bc_cv(js, jl, jk) * time_step_len
                  
                  ELSE
                     xbc1_p(jl) = MAX(xbc1_p(jl), 0._dp) + qxtm1(jl, jk,idt) + &
                           (qxtte(jl, jk,idt) - ttescav_bc_cv(js, jl, jk)) * &
                           time_step_len
                  END IF
               END IF

            END IF
         END DO ! naer

         IF (.NOT. l_mode_partition) THEN
            idt = idt_svocp(jt)

            IF (lcolumn) THEN
               xtp1_p(jl) = trac_field(jl, jk,idt,jrow)
            ELSE
               xtp1_p(jl) = qxtm1(jl, jk,idt) + &
                            qxtte(jl, jk,idt) * time_step_len
            END IF

            IF (xbc1_p(jl) > 1.e-34_dp .AND. xtp1_p(jl) > 1.e-34_dp) &
               tot_ttescav(jt, jl, jk) = &
                  xtp1_p(jl) * SUM(ttescav_bc_cv(:, jl, jk)) / xbc1_p(jl)

            IF (xtp1_p(jl) + tot_ttescav(jt, jl, jk) * &
                time_step_len < 0._dp) tot_ttescav(jt, jl, jk) = &
                -xtp1_p(jl) / time_step_len
         END IF

      END DO ! kproma
   END DO ! nlev

   DO jl = 1, kproma
      tot_flxpos = 0._dp
      tot_flxneg = 0._dp
      fac = 1.0_dp

      DO jk = 1, nlev
         IF (l_mode_partition) THEN
            DO jm = 1, nmod
               IF (ttescav_m(jt, jl, jk, jm) > 0._dp) THEN
                  tot_flxpos = tot_flxpos + ttescav_m(jt, jl, jk, jm) * &
                               rho_air(jl) * zdz(jl, jk)
               ELSE
                  tot_flxneg = tot_flxneg + ttescav_m(jt, jl, jk, jm) * &
                               rho_air(jl) * zdz(jl, jk)
               END IF
            END DO
      
         ELSE
            IF (tot_ttescav(jt, jl, jk) > 0._dp) THEN
               tot_flxpos = tot_flxpos + tot_ttescav(jt, jl, jk) * &
                            rho_air(jl) * zdz(jl, jk)
            ELSE
               tot_flxneg = tot_flxneg + tot_ttescav(jt, jl, jk) * &
                            rho_air(jl) * zdz(jl, jk)
            END IF
         END IF
      END DO

      IF (tot_flxpos > ABS(tot_flxneg) .AND. ABS(tot_flxneg) > 0._dp) &
         fac = ABS(tot_flxneg) / tot_flxpos

      DO jk = 1, nlev
         IF (l_mode_partition) THEN
            DO jm = 1, nmod
               IF (ttescav_m(jt, jl, jk, jm) > 0._dp) &
                 ttescav_m(jt, jl, jk, jm) = ttescav_m(jt, jl, jk, jm) * fac
            END DO

         ELSE
            IF (tot_ttescav(jt, jl, jk) > 0._dp) &
               tot_ttescav(jt, jl, jk) = tot_ttescav(jt, jl, jk) * fac
         END IF
      END DO
   END DO

   DO jk = 1, nlev
      IF (.NOT. l_mode_partition) THEN
         idt = idt_svocp(jt)

         ! update tracer tendency
         IF (lcolumn) THEN
            trac_field(1: kproma, jk,idt,jrow) = trac_field(1: kproma, jk,idt,jrow) + &
                                tot_ttescav(jt, 1: kproma, jk) * time_step_len
          
         ELSE
            qxtte(1: kproma, jk,idt) = qxtte(1: kproma, jk,idt) + &
                                tot_ttescav(jt, 1: kproma, jk)
         END IF

         ! wet deposition flux [mol(trac) mol(air)-1 kg(air) m-2 s-1]
         wetflx_aer(1: kproma, jk) = tot_ttescav(jt, 1: kproma, jk) * &
                                rho_air(1: kproma) * zdz(1: kproma, jk)

         tot_ttescav_cv(jt, :, jk) = tot_ttescav(jt, 1: kproma, jk)
      ELSE
         DO jm = 1, nmod
            idt = idt_svocm(jt, jm)

            ! update tracer tendency
            IF (lcolumn) THEN
               trac_field(1: kproma, jk,idt,jrow) = trac_field(1: kproma, jk,idt,jrow) + &
                                ttescav_m(jt, 1: kproma, jk, jm) * time_step_len
            ELSE
               qxtte(1: kproma, jk,idt) = qxtte(1: kproma, jk,idt) + &
                                ttescav_m(jt, 1: kproma, jk, jm)
            END IF

            ! wet deposition flux [mol(trac) mol(air)-1 kg(air) m-2 s-1]
            wetflx_aer(1: kproma, jk) = wetflx_aer(1: kproma, jk) + &
                    ttescav_m(jt, 1: kproma, jk, jm) * rho_air(1: kproma) * &
                    zdz(1: kproma, jk)

            ttescav_m_cv(jt, :, jk, jm) = ttescav_m(jt, 1: kproma, jk, jm)
         END DO
      END IF

      wetflx_aer_tot(jt, 1: kproma) = wetflx_aer_tot(jt, 1: kproma) + &
                           wetflx_aer(1: kproma, jk)
   END DO !nlev
END DO ! nsvoc

CALL end_message_bi(modstr, 'SVOC CONVEC ROUTINES', substr)

END SUBROUTINE svoc_convec

! *************************************************************************

SUBROUTINE svoc_physc(iflag)

!!!

USE messy_main_grid_def_mem_bi, ONLY: nproma, kproma, nlev, jrow
USE messy_main_grid_def_bi,     ONLY: deltaz
USE messy_main_data_bi,         ONLY: loland_2d, loglac_2d, &
                                      cvs, vgrat, &
                                      tm1_3d, tte_3d, &
                                      qm1_3d, qte_3d, &
                                      press_3d, pressi_3d, &
                                      rhum_3d
USE messy_main_tracer,        ONLY: R_MOLARMASS
USE messy_main_timer,         ONLY: delta_time, time_step_len, lstart
USE messy_main_tracer_mem_bi, ONLY: qxtm1, & ! mole frac at t-dt
                                    qxtte    ! tracer tendencies
USE messy_main_constants_mem, ONLY: pi, g, M_air, R_gas, N_A, &
                                    rd, & ! 1000. * R_gas/M_air (J/K/kg)
                                    vtmpc1 ! M_air / M_H2O - 1.
USE messy_main_tracer_mem_bi, ONLY: ti_gp


IMPLICIT NONE

CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_physc'
INTEGER, INTENT(IN) :: iflag
REAL(dp), PARAMETER :: kg2ug = 1.e+9_dp
REAL(dp), PARAMETER :: oc2om = 1._dp
!REAL(dp), PARAMETER :: oc2om = 1.724
REAL(dp), PARAMETER :: g2kg = 1.e-3_dp
REAL(dp), PARAMETER :: m2cm = 1.e+2_dp
REAL(dp), DIMENSION(nsvoc, nproma, nlev, nmod) :: ttescav_m
REAL(dp), DIMENSION(nsvoc, nproma, nlev) :: tot_ttescav
REAL(dp), DIMENSION(nproma, nlev, nmod) :: conctsp, concwsom, concwiom
REAL(dp), DIMENSION(nproma, nlev, nmod) :: concbc, concbc_wgt
REAL(dp), DIMENSION(nproma, nlev, nmod) :: concso4mm, concss
REAL(dp), DIMENSION(nproma, nlev, nsvoc) :: zxttetot, zxtm1tot
REAL(dp), DIMENSION(nproma, nlev, nmod) :: rarea
REAL(dp), DIMENSION(nproma, nlev) :: tot_concwsom, tot_concwiom
REAL(dp), DIMENSION(nproma, nlev) :: tot_conctsp, tot_concbc
REAL(dp), DIMENSION(nproma, nlev) :: tot_concso4mm, tot_concss
REAL(dp), DIMENSION(nproma, nlev) :: rtot_area
REAL(dp), DIMENSION(nproma, nlev) :: pzo3, zdz
REAL(dp), DIMENSION(nproma, nlev) :: wetflx_aer
REAL(dp), DIMENSION(nproma, nlev) :: kptot, tsptot
REAL(dp), DIMENSION(nproma, nmod) :: concbc_m
REAL(dp), DIMENSION(nproma) :: xtp1_g, xtp1_p, xbc1_p
REAL(dp), DIMENSION(nproma) :: dcdt, abur, rho_air
REAL(dp), DIMENSION(nproma) :: tot_bc, tot_svoc
REAL(dp), DIMENSION(nproma) :: wetflx_g_ls, wetflx_g_cv
REAL(dp), DIMENSION(nmod) :: bc_xtm1, bc_xtte, bc_xtm1_wgt, bc_xtte_wgt
REAL(dp), DIMENSION(nmod) :: oc_ws_xtm1, oc_wi_xtm1, oc_ws_xtte, oc_wi_xtte
REAL(dp), DIMENSION(nmod) :: so4mm_xtm1, so4mm_xtte, ss_xtm1, ss_xtte
REAL(dp), DIMENSION(nmod) :: tsp_xtm1, tsp_xtte
CHARACTER(LEN = 2) :: txt = ''
CHARACTER(LEN = 24) :: trname
REAL(dp) :: spc_hum, temp, conv, mol2mass
REAL(dp) :: cmr2ras, zaernl, zaerml
REAL(dp) :: ztotmass, zbcmass
REAL(dp) :: tot_flxpos, tot_flxneg, fac
INTEGER :: jl, jk, jt, js, jr, jm, jn
INTEGER :: idt, idt2, idx

CALL start_message_bi(modstr, 'SVOC PHYSC ROUTINES', substr)

ttescav_m = 0._dp
tot_ttescav = 0._dp

SELECT CASE (iflag)
CASE (1) ! after scav_physc(1) and before m7/gmxe_physc
   !!!
   !!! update tracer tendency due to large-scale scavenging
   !!! need coupling with scav sm
   !!!

   IF (lstart) THEN
      wetflx_aer_tot = 0._dp
      ttescav_bc_cv = 0._dp
      tot_ttescav_cv = 0._dp
      ttescav_m_cv = 0._dp
   END IF

   tot_ttescav_ls = 0._dp
   ttescav_m_ls = 0._dp

   DO jt = 1, nsvoc
      DO js = 1, 4
         ttescav_bc_ls(js, :, :) = 0._dp
      END DO

      wetflx_aer(:, :) = 0._dp
      zdz(:, :) = 0._dp

      ! layer thickness [m]
      zdz(:, 1: nlev) = deltaz(1: kproma, 1:nlev,jrow)

      DO jk = 1, nlev
         xtp1_p(:) = 0._dp
         xbc1_p(:) = 0._dp

         DO jl = 1, kproma
            spc_hum = qm1_3d(jl, jk,jrow) + &
                      qte_3d(jl, jk,jrow) * time_step_len

            temp = tm1_3d(jl, jk,jrow) + &
                   tte_3d(jl, jk,jrow) * time_step_len

            ! density of moist air [kg m-3]
            ! the factor 1+vtmpc1*spc_hum accounts for moist correction
            rho_air(jl) = press_3d(jl, jk,jrow) / &
                         (temp * rd * (1._dp + vtmpc1 * spc_hum))

            js = 0

            DO jr = 1, naer
               IF (TRIM(bname_aer(jr)) == 'BC') THEN
                  js = js + 1
                  idt = idt_maer(jr)

                  DO jm = 1, nmod
                     IF (TRIM(sname_aer(jr)) == TRIM(cmode(jm))) THEN
                        jn = jm
                        EXIT
                     END IF
                  END DO

                  ttescav_bc_ls(js, jl, jk) = ttescav_bc(js)%ptr(jl, jk,jrow) - &
                                              ttescav_bc_cv(js, jl, jk)

                  IF (l_mode_partition) THEN
                     ! bc conc after scav_ls
                     xbc1_p(jl) = qxtm1(jl, jk,idt) + (qxtte(jl, jk,idt) - &
                           ttescav_bc_ls(js, jl, jk)) * time_step_len

                     idt2 = idt_svocm(jt, jn)

                     xtp1_p(jl) = qxtm1(jl, jk,idt2) + &
                                  qxtte(jl, jk,idt2) * time_step_len
 
                     IF (xbc1_p(jl) > 1.e-34_dp .AND. xtp1_p(jl) > 1.e-34_dp) &
                        ttescav_m(jt, jl, jk, jn) = &
                           xtp1_p(jl) * ttescav_bc_ls(js, jl, jk) / xbc1_p(jl)

                     IF (xtp1_p(jl) + ttescav_m(jt, jl, jk, jn) * &
                        time_step_len < 0._dp) ttescav_m(jt, jl, jk, jn) = &
                        -xtp1_p(jl) / time_step_len

                  ELSE ! bulk
                     ! bc conc (sum all modes)
                     xbc1_p(jl) = MAX(xbc1_p(jl), 0._dp) + &
                        qxtm1(jl, jk,idt) + (qxtte(jl, jk,idt) - &
                        ttescav_bc_ls(js, jl, jk)) * time_step_len
                  END IF
                
               END IF
            END DO ! naer

            IF (.NOT. l_mode_partition) THEN
               idt2 = idt_svocp(jt)

               xtp1_p(jl) = qxtm1(jl, jk,idt2) + &
                            qxtte(jl, jk,idt2) * time_step_len

               IF (xbc1_p(jl) > 1.e-34_dp .AND. xtp1_p(jl) > 1.e-34_dp) &
                  tot_ttescav(jt, jl, jk) = &
                      xtp1_p(jl) * SUM(ttescav_bc_ls(:, jl, jk)) / xbc1_p(jl)

               IF (xtp1_p(jl) + tot_ttescav(jt, jl, jk) * &
                  time_step_len < 0._dp) tot_ttescav(jt, jl, jk) = &
                     -xtp1_p(jl) / time_step_len
            END IF
         END DO ! kproma
      END DO ! nlev

      DO jl = 1, kproma
         tot_flxpos = 0._dp
         tot_flxneg = 0._dp
         fac = 1.0_dp

         DO jk = 1, nlev
            IF (l_mode_partition) THEN
               DO jm = 1, nmod
                  IF (ttescav_m(jt, jl, jk, jm) > 0._dp) THEN
                     tot_flxpos = tot_flxpos + ttescav_m(jt, jl, jk, jm) * &
                                  rho_air(jl) * zdz(jl, jk)
                  ELSE
                     tot_flxneg = tot_flxneg + ttescav_m(jt, jl, jk, jm) * &
                                  rho_air(jl) * zdz(jl, jk)
                  END IF
               END DO

            ELSE
               IF (tot_ttescav(jt, jl, jk) > 0._dp) THEN
                  tot_flxpos = tot_flxpos + tot_ttescav(jt, jl, jk) * &
                               rho_air(jl) * zdz(jl, jk)
               ELSE
                  tot_flxneg = tot_flxneg + tot_ttescav(jt, jl, jk) * &
                               rho_air(jl) * zdz(jl, jk)
               END IF
            END IF
         END DO

         IF (tot_flxpos > ABS(tot_flxneg) .AND. ABS(tot_flxneg) > 0._dp) &
            fac = ABS(tot_flxneg) / tot_flxpos

         DO jk = 1, nlev
            IF (l_mode_partition) THEN
               DO jm = 1, nmod
                  IF (ttescav_m(jt, jl, jk, jm) > 0._dp) &
                     ttescav_m(jt, jl, jk, jm) = ttescav_m(jt, jl, jk, jm) * fac
               END DO

            ELSE
               IF (tot_ttescav(jt, jl, jk) > 0._dp) &
                  tot_ttescav(jt, jl, jk) = tot_ttescav(jt, jl, jk) * fac
            END IF
         END DO
      END DO

      DO jk = 1, nlev
         IF (.NOT. l_mode_partition) THEN
            idt = idt_svocp(jt)

            ! update tracer tendency
            qxtte(1: kproma, jk,idt) = qxtte(1: kproma, jk,idt) + &
                       tot_ttescav(jt, 1: kproma, jk)

            ! wet deposition flux [mol(trac) mol(air)-1 kg(air) m-2 s-1]
            wetflx_aer(1: kproma, jk) = tot_ttescav(jt, 1: kproma, jk) * &
                       rho_air(1: kproma) * zdz(1: kproma, jk)

            tot_ttescav_ls(jt, :, jk) = tot_ttescav(jt, 1: kproma, jk)

         ELSE
            DO jm = 1, nmod
               idt2 = idt_svocm(jt, jm)

               ! update tracer tendency
               qxtte(1: kproma, jk,idt2) = qxtte(1: kproma, jk,idt2) + &
                       ttescav_m(jt, 1: kproma, jk, jm)

               ! wet deposition flux [mol(trac) mol(air)-1 kg(air) m-2 s-1]
               wetflx_aer(1: kproma, jk) = wetflx_aer(1: kproma, jk) + &
                       ttescav_m(jt, 1: kproma, jk, jm) * rho_air(1: kproma) * &
                       zdz(1: kproma, jk)

               ttescav_m_ls(jt, :, jk, jm) = ttescav_m(jt, 1: kproma, jk, jm)
            END DO
         END IF

         wetflx_aer_tot(jt, 1: kproma) = wetflx_aer_tot(jt, 1: kproma) + &
                              wetflx_aer(1: kproma, jk)
      END DO ! nlev
   END DO ! nsvoc

CASE (2) ! before calling mecca_physc but after m7/gmxe_physc
   !!!
   !!! chemical degradation of svoc particle
   !!!

   rarea = 0._dp
   rtot_area = 0._dp

   DO jk = 1, nlev
      DO jl = 1, kproma
         spc_hum = qm1_3d(jl, jk,jrow) + qte_3d(jl, jk,jrow) * time_step_len
         temp = tm1_3d(jl, jk,jrow) + tte_3d(jl, jk,jrow) * time_step_len

         rho_air(jl) = press_3d(jl, jk,jrow) / &
                      (temp * rd * (1._dp + vtmpc1 * spc_hum))

         conv = (press_3d(jl, jk,jrow) / (R_gas * temp)) * 1.e-6_dp * N_A

         ! convert offline o3 concentration from vmr [mol mol-1] into
         ! number density [molec cm-3]
         idt = idt_o3
         pzo3(jl, jk) = conv * (qxtm1(jl, jk,idt) + &
                        qxtte(jl, jk,idt) * time_step_len)

         pzo3(jl, jk) = MAX(0._dp, pzo3(jl, jk))

         DO jm = 1, nmod
            ! convert count median radius to radius of average surface
            cmr2ras = EXP(1.0_dp * (LOG(sigma(jm))) ** 2)

            ! calculate total number of aerosol in each mode
            IF (aermod_str == 'gmxe') THEN
               zaernl =  aernl_cpl(jl, jk,jm,jrow)
            ELSE IF (aermod_str == 'm7') THEN
               idt = idt_nnum(jm)

               ! convert the number from mol-1 to cm-3 using M_air
               ! 1.e-6:  factor to transform m-3 to cm-3
               zaernl = rho_air(jl) * 1.e-6_dp * (qxtm1(jl, jk,idt) + &
                     qxtte(jl, jk,idt) * time_step_len) / (M_air * g2kg)
            END IF

            ! remove negative values
            zaernl = MAX(zaernl, 0._dp)

            ! calculate total surface area of aerosols per cm3 of air [cm2 cm-3]
            ! rwet_aer: count median radius, should be corrected to radius of 
            ! average surface using cmr2ras and rwet should be in cm (*100)
            IF (rwet_aer(jl, jk,jm,jrow) > 0._dp) THEN
               rarea(jl, jk, jm) = 4._dp * pi * &
                                 ((rwet_aer(jl, jk,jm,jrow) * m2cm * &
                                  cmr2ras) ** 2) * zaernl

               rtot_area(jl, jk) = rtot_area(jl, jk) + rarea(jl, jk, jm)
            END IF

         END DO ! nmod
      END DO ! kproma
   END DO ! nlev

   J_out_NO2 = J_in_NO2

   IF (.NOT. l_mode_partition) THEN
      DO jt = 1, nsvoc
         idt = idt_svocp(jt)
         trname = TRIM(SVOC_NAME(jt))

         DO jk = 1, nlev
            ! calculate particle concentration
            xtp1_p(1: kproma) = qxtm1(1: kproma, jk,idt) + &
                                qxtte(1: kproma, jk,idt) * time_step_len

            ! discard negative values
            POSTV4: FORALL (jl = 1: kproma, xtp1_p(jl) < 0._dp)
                      xtp1_p(jl) = 0._dp
                   END FORALL POSTV4

            ! oxidation with o3, and the oxidant is treated as constant
            dcdt(:) = 0._dp

            IF (trname(1:3) == 'BaP') THEN
               CALL svoc_hetchem_bap(kproma, &
                                     param_bapo3, &
                                     press_3d(1: kproma, jk,jrow), &
                                     rhum_3d(1: kproma, jk,jrow), &
                                     tm1_3d(1: kproma, jk,jrow) + &
                                     tte_3d(1: kproma, jk,jrow) * &
                                     time_step_len, &
                                     pzo3(1: kproma, jk), &
                                     xtp1_p(1: kproma), &
                                     dcdt(1: kproma))

            ELSE IF (trname(1:3) == 'PHE') THEN
               CALL svoc_hetchem_phe(kproma, &
                                     rtot_area(1: kproma, jk), &
                                     pzo3(1: kproma, jk), &
                                     xtp1_p(1: kproma), &
                                     dcdt(1: kproma))
            !! mz_jw_20160918+
            ELSE IF (l_pahderiv .AND. &
                   .NOT. l_emiss(jt) ) THEN ! npah formation from parent pah
               CALL svoc_hetchem_npah(kproma, &
                                      J_in_NO2(1: kproma, jk,jrow), &
                                      xtp1_p(1: kproma), &
                                      dcdt(1: kproma))
            !! mz_jw_20160918-
            END IF

            ! update tracer tendency
            qxtte(1: kproma, jk,idt) = qxtte(1: kproma, jk,idt) + &
                                        dcdt(1: kproma)
         END DO
      END DO
   
   ELSE
      DO jt = 1, nsvoc
         trname = TRIM(SVOC_NAME(jt))

         DO jm = 1, nmod
            idt = idt_svocm(jt, jm)

            DO jk = 1, nlev
               ! calculate particle concentration
               xtp1_p(1: kproma) = qxtm1(1: kproma, jk,idt) + &
                                   qxtte(1: kproma, jk,idt) * time_step_len

               ! discard negative values
               POSTV5: FORALL (jl = 1: kproma, xtp1_p(jl) < 0._dp)
                         xtp1_p(jl) = 0._dp
                      END FORALL POSTV5

               ! oxidation with o3, and the oxidant is treated as constant
               dcdt(:) = 0._dp

               IF (trname(1:3) == 'BaP') THEN
                  CALL svoc_hetchem_bap(kproma, &
                                        param_bapo3, &
                                        press_3d(1: kproma, jk,jrow), &
                                        rhum_3d(1: kproma, jk,jrow), &
                                        tm1_3d(1: kproma, jk,jrow) + &
                                        tte_3d(1: kproma, jk,jrow) * &
                                        time_step_len, &
                                        pzo3(1: kproma, jk), &
                                        xtp1_p(1: kproma), &
                                        dcdt(1: kproma))

               ELSE IF (trname(1:3) == 'PHE') THEN
                  CALL svoc_hetchem_phe(kproma, &
                                        rarea(1: kproma, jk, jm), &
                                        pzo3(1: kproma, jk), &
                                        xtp1_p(1: kproma), &
                                        dcdt(1: kproma))
               !! mz_jw_20160918+
               ELSE IF (l_pahderiv .AND. &
                       .NOT. l_emiss(jt) ) THEN ! npah formation from parent pah
                  CALL svoc_hetchem_npah(kproma, &
                                         J_in_NO2(1: kproma, jk,jrow), &
                                         xtp1_p(1: kproma), &
                                         dcdt(1: kproma))
               !! mz_jw_20160918-
               END IF

               ! update tracer tendency
               qxtte(1: kproma, jk,idt) = qxtte(1: kproma, jk,idt) + &
                                           dcdt(1: kproma)
            END DO ! nlev

         END DO ! nmod

      END DO ! nsvoc
   END IF ! l_mode_partition

CASE (3)  ! after calling scav_physc(2)
   !!!
   !!! update tracer tendency due to cloud evaporation and 
   !!! calculate total wet deposition flux
   !!! need coupling with scav sm
   !!!

   DO jt = 1, nsvoc
      DO js = 1, 4
         ttescav_bc_ev(js, :, :) = 0._dp
      END DO

      wetflx_aer(:, :) = 0._dp
      wetflx_g_ls(:) = 0._dp
      wetflx_g_cv(:) = 0._dp
      zdz(:, :) = 0._dp

      ! layer thickness [m]
      zdz(:, 1: nlev) = deltaz(1: kproma, 1:nlev,jrow)

      !!!  calculate total scavenging tendency of particles
      DO jk = 1, nlev
         xtp1_p(:) = 0._dp
         xbc1_p(:) = 0._dp

         DO jl = 1, kproma
            spc_hum = qm1_3d(jl, jk,jrow) + &
                      qte_3d(jl, jk,jrow) * time_step_len

            temp = tm1_3d(jl, jk,jrow) + &
                   tte_3d(jl, jk,jrow) * time_step_len

            ! density of moist air [kg m-3]
            ! the factor 1+vtmpc1*spc_hum accounts for moist correction
            rho_air(jl) = press_3d(jl, jk,jrow) / &
                          (temp * rd * (1._dp + vtmpc1 * spc_hum))

            js = 0

            DO jr = 1, naer
               IF (TRIM(bname_aer(jr)) == 'BC') THEN
                  js = js + 1
                  idt = idt_maer(jr)

                  DO jm = 1, nmod
                     IF (TRIM(sname_aer(jr)) == TRIM(cmode(jm))) THEN
                        jn = jm
                        EXIT
                     END IF
                  END DO
                        
                  ttescav_bc_ev(js, jl, jk) = ttescav_bc(js)%ptr(jl, jk,jrow) - &
                           ttescav_bc_cv(js, jl, jk) - ttescav_bc_ls(js, jl, jk)

                  IF (l_mode_partition) THEN
                     ! bc total fluxes (here I used tendency although rho_air
                     ! in svoc_convec and svoc_physc is different)
                     xbc1_p(jl) = SUM(ttescav_bc_cv(1: js, jl, jk)) + &
                               SUM(ttescav_bc_ls(1: js, jl, jk)) - &
                               SUM(ttescav_bc_ev(1: js-1, jl, jk))

                     xbc1_p(jl) = MIN(0._dp, xbc1_p(jl))
       
                     idt2 = idt_svocm(jt, jn)

                     IF (jn /= 5) THEN !!! no evaporation at soluble mode
                        xtp1_p(jl) = SUM(ttescav_m_cv(jt, jl, jk, 2: jn)) + &
                                  SUM(ttescav_m_ls(jt, jl, jk, 2: jn)) + &
                                  ttescav_m_cv(jt, jl, jk, 5) + &
                                  ttescav_m_ls(jt, jl, jk, 5) - & 
                                  SUM(ttescav_m(jt, jl, jk, 2: jn-1))

                        xtp1_p(jl) = MIN(0._dp, xtp1_p(jl))
                     END IF

                     IF (ABS(xbc1_p(jl)) > 1.e-34_dp .AND. &
                        ABS(xtp1_p(jl)) > 1.e-34_dp) ttescav_m(jt, jl, jk, jn) = &
                           xtp1_p(jl) * ttescav_bc_ev(js, jl, jk) / xbc1_p(jl) 
                  END IF
               END IF
            END DO ! naer

            IF (.NOT. l_mode_partition) THEN
               xbc1_p(jl) = SUM(ttescav_bc_cv(:, jl, jk)) + &
                            SUM(ttescav_bc_ls(:, jl, jk)) 

               xtp1_p(jl) = tot_ttescav_cv(jt, jl, jk) + &
                            tot_ttescav_ls(jt, jl, jk) 

               IF (ABS(xbc1_p(jl)) > 1.e-34_dp .AND. &
                   ABS(xtp1_p(jl)) > 1.e-34_dp) tot_ttescav(jt, jl, jk) = &
                       xtp1_p(jl) * SUM(ttescav_bc_ev(:, jl, jk)) / xbc1_p(jl)
            END IF
         END DO ! kproma
      END DO ! nlev

      DO jl = 1, kproma
         tot_flxpos = 0._dp
         tot_flxneg = 0._dp
         fac = 1.0_dp

         DO jk = 1, nlev
            IF (l_mode_partition) THEN
               DO jm = 1, nmod
                  IF (ttescav_m_cv(jt, jl, jk, jm) + &
                     ttescav_m_ls(jt, jl, jk, jm) + &
                     ttescav_m(jt, jl, jk, jm) > 0._dp) THEN
                     tot_flxpos = tot_flxpos + (ttescav_m_cv(jt, jl, jk, jm) + &
                          ttescav_m_ls(jt, jl, jk, jm) + &
                          ttescav_m(jt, jl, jk, jm)) * rho_air(jl) * zdz(jl, jk)
                  ELSE
                     tot_flxneg = tot_flxneg + (ttescav_m_cv(jt, jl, jk, jm) + &
                          ttescav_m_ls(jt, jl, jk, jm) + &
                          ttescav_m(jt, jl, jk, jm)) * rho_air(jl) * zdz(jl, jk)
                  END IF
               END DO
   
            ELSE
               IF (tot_ttescav_cv(jt, jl, jk) + &
                   tot_ttescav_ls(jt, jl, jk) + &
                   tot_ttescav(jt, jl, jk) > 0._dp) THEN
                  tot_flxpos = tot_flxpos + (tot_ttescav_cv(jt, jl, jk) + &
                      tot_ttescav_ls(jt, jl, jk) + tot_ttescav(jt, jl, jk)) * &
                      rho_air(jl) * zdz(jl, jk)
               ELSE
                  tot_flxneg = tot_flxneg + (tot_ttescav_cv(jt, jl, jk) + &
                      tot_ttescav_ls(jt, jl, jk) + tot_ttescav(jt, jl, jk)) * &
                      rho_air(jl) * zdz(jl, jk)
               END IF
            END IF
         END DO

         IF (tot_flxpos > ABS(tot_flxneg) .AND. ABS(tot_flxneg) > 0._dp) &
            fac = ABS(tot_flxneg) / tot_flxpos

         DO jk = 1, nlev
            IF (l_mode_partition) THEN
               DO jm = 1, nmod
                  IF (ttescav_m_cv(jt, jl, jk, jm) + &
                     ttescav_m_ls(jt, jl, jk, jm) + &
                     ttescav_m(jt, jl, jk, jm) > 0._dp) &
                     ttescav_m(jt, jl, jk, jm) = ttescav_m(jt, jl, jk, jm) * fac
               END DO

            ELSE
               IF (tot_ttescav_cv(jt, jl, jk) + &
                  tot_ttescav_ls(jt, jl, jk) + &
                  tot_ttescav(jt, jl, jk) > 0._dp) &
                  tot_ttescav(jt, jl, jk) = tot_ttescav(jt, jl, jk) * fac
           END IF
         END DO
      END DO

      DO jk = 1, nlev
         IF (.NOT. l_mode_partition) THEN
            idt = idt_svocp(jt)

            ! update tracer tendency
            qxtte(1: kproma, jk,idt) = qxtte(1: kproma, jk,idt) + &
                                        tot_ttescav(jt, 1: kproma, jk)

            ! wet deposition flux [mol(trac) mol(air)-1 kg(air) m-2 s-1]
            wetflx_aer(1: kproma, jk) = tot_ttescav(jt, 1: kproma, jk) * &
                             rho_air(1: kproma) * zdz(1: kproma, jk)
         ELSE
            DO jm = 1, nmod
               idt = idt_svocm(jt, jm)
 
               ! update tracer tendency
               qxtte(1: kproma, jk,idt) = qxtte(1: kproma, jk,idt) + &
                                           ttescav_m(jt, 1: kproma, jk, jm)

               ! wet deposition flux [mol(trac) mol(air)-1 kg(air) m-2 s-1]
               wetflx_aer(1: kproma, jk) = wetflx_aer(1: kproma, jk) + &
                       ttescav_m(jt, 1: kproma, jk, jm) * rho_air(1: kproma) * &
                       zdz(1: kproma, jk)
            END DO
         END IF

         wetflx_aer_tot(jt, 1: kproma) = wetflx_aer_tot(jt, 1: kproma) + &
                               wetflx_aer(1: kproma, jk)
      END DO ! nlev

      wetflx_aer_tot(jt, 1: kproma) = -1._dp * wetflx_aer_tot(jt, 1: kproma)

      POSTV6: FORALL (jl = 1: kproma, wetflx_aer_tot(jt, jl) < 0._dp)
               wetflx_aer_tot(jt, jl) = 0._dp
            END FORALL POSTV6

      !!! update deposition and surface burden 
      depo_p => depo(jt)%ptr
      burd_p => burd(jt)%ptr

      jk = nlev - 1

      wetflx_g_ls(1: kproma) = wetflx_ls(jt)%ptr(1: kproma ,jrow) + &
                               wetflx_ls(jt+nsvoc)%ptr(1: kproma ,jrow)

      wetflx_g_cv(1: kproma) = wetflx_cv(jt)%ptr(1: kproma ,jrow) + &
                               wetflx_cv(jt+nsvoc)%ptr(1: kproma ,jrow)

      CALL svoc_wetdep(kproma, nlev, delta_time, &
                       MOLMASS(jt), &
                       loland_2d(1: kproma ,jrow), &
                       loglac_2d(1: kproma ,jrow), &
                       cvs(1: kproma ,jrow), &
                       vgrat(1: kproma ,jrow), &
                       wetflx_g_ls(1: kproma), &
                       wetflx_g_cv(1: kproma), &
                       wetflx_aer_tot(jt, 1: kproma), &
                       burd_p(1: kproma, 1:nlev,jrow), &
                       depo_p(1: kproma, jk,jrow))

      ! update total deposition (dry + wet)
      ! deposition lev3 (dry+wet) = lev2 (wet) + lev1 (dry)
      depo_p(1: kproma, jk-1,jrow) = depo_p(1: kproma, jk,jrow) + &
                                     depo_p(1: kproma, jk+1,jrow)
   END DO ! nsvoc

   !!!
   !!! gas-particle partitioning
   !!!

   rarea = 0._dp
   rtot_area = 0._dp

   DO jk = 1, nlev
      DO jl = 1, kproma
         spc_hum = qm1_3d(jl, jk,jrow) + &
                   qte_3d(jl, jk,jrow) * time_step_len

         temp = tm1_3d(jl, jk,jrow) + &
                tte_3d(jl, jk,jrow) * time_step_len

         ! density of moist air [kg m-3]
         rho_air(jl) = press_3d(jl, jk,jrow) / &
                       (temp * rd * (1._dp + vtmpc1 * spc_hum))

         DO jm = 1, nmod
            ! convert count median radius to radius of average surface based
            ! on the hatch-choate conversion (see hinds, equation 4.52)
            cmr2ras = EXP(1.0_dp * (LOG(sigma(jm))) ** 2)

            ! calculate total number of aerosol in each mode
            IF (aermod_str == 'gmxe') THEN
               zaernl =  aernl_cpl(jl, jk,jm,jrow)
            ELSE IF (aermod_str == 'm7') THEN
               idt = idt_nnum(jm)

               ! convert the number from mol-1 to cm-3 using M_air
               ! 1.e-6:  factor to transform m-3 to cm-3
               zaernl = rho_air(jl) * 1.e-6_dp * (qxtm1(jl, jk,idt) + &
                        qxtte(jl, jk,idt) * time_step_len) / (M_air * g2kg)
            END IF
         
            ! remove negative values
            zaernl = MAX(zaernl, 0._dp)

            ! calculate total surface area of aerosols per cm3 of air [cm2 cm-3]
            ! rwet_aer: count median radius, should be corrected to radius of 
            ! average surface using cmr2ras and rwet should be in cm (*100)
            IF (rwet_aer(jl, jk,jm,jrow) > 0._dp) THEN
               rarea(jl, jk, jm) = 4._dp * pi * &
                                 ((rwet_aer(jl, jk,jm,jrow) * m2cm * &
                                  cmr2ras) ** 2) * zaernl

               rtot_area(jl, jk) = rtot_area(jl, jk) + rarea(jl, jk, jm)
            END IF
         END DO ! nmod

         tsp_xtm1 = 0._dp
         tsp_xtte = 0._dp
         oc_ws_xtm1 = 0._dp
         oc_ws_xtte = 0._dp
         oc_wi_xtm1 = 0._dp
         oc_wi_xtte = 0._dp
         bc_xtm1 = 0._dp
         bc_xtte = 0._dp
         bc_xtm1_wgt = 0._dp
         bc_xtte_wgt = 0._dp
         so4mm_xtm1 = 0._dp
         so4mm_xtte = 0._dp
         ss_xtm1 = 0._dp
         ss_xtte = 0._dp

         DO jr = 1, naer
            idt = idt_maer(jr)

            ! create factor MW/M_air to convert volume mixing ratio
            ! [mol mol-1] into mass mixing ratio [kg kg-1]
            mol2mass = ti_gp(idt)%tp%meta%cask_r(R_MOLARMASS) / M_air

            DO jm = 1, nmod
               IF (TRIM(sname_aer(jr)) == TRIM(cmode(jm))) THEN
                  jn = jm
                  EXIT
               END IF
            END DO

            tsp_xtm1(jn) = tsp_xtm1(jn) + qxtm1(jl, jk,idt) * mol2mass
            tsp_xtte(jn) = tsp_xtte(jn) + qxtte(jl, jk,idt) * mol2mass

            IF (TRIM(bname_aer(jr)) == 'OC'     .OR. &
                TRIM(bname_aer(jr)) == 'WSOC01' .OR. & 
                TRIM(bname_aer(jr)) == 'WSOC02' .OR. & 
                TRIM(bname_aer(jr)) == 'WSOC03' .OR. & 
                TRIM(bname_aer(jr)) == 'WSOC04' .OR. & 
                TRIM(bname_aer(jr)) == 'WSOC05') THEN

               IF (TRIM(sname_aer(jr)) == 'ks' .OR. &
                  TRIM(sname_aer(jr)) == 'as' .OR. &
                  TRIM(sname_aer(jr)) == 'cs') THEN    !!! soluble   
                  oc_ws_xtm1(jn) = oc_ws_xtm1(jn) + qxtm1(jl, jk,idt) * mol2mass
                  oc_ws_xtte(jn) = oc_ws_xtte(jn) + qxtte(jl, jk,idt) * mol2mass

               ELSE  !!! insoluble
                  oc_wi_xtm1(jn) = oc_wi_xtm1(jn) + qxtm1(jl, jk,idt) * mol2mass
                  oc_wi_xtte(jn) = oc_wi_xtte(jn) + qxtte(jl, jk,idt) * mol2mass
               END IF

            ELSE IF (TRIM(bname_aer(jr)) == 'BC') THEN
               bc_xtm1(jn) = qxtm1(jl, jk,idt) * mol2mass
               bc_xtte(jn) = qxtte(jl, jk,idt) * mol2mass

               ! weighted average of bc; the weight is defined as the ratio of
               ! aerosol surface area of each mode to the total surface area
               IF (rtot_area(jl, jk) .GE. 1.e-20_dp) THEN
                  bc_xtm1_wgt(jn) = qxtm1(jl, jk,idt) * &
                        rarea(jl, jk, jn) / rtot_area(jl, jk) * mol2mass

                  bc_xtte_wgt(jn) = qxtte(jl, jk,idt) * &
                        rarea(jl, jk, jn) / rtot_area(jl, jk) * mol2mass
               END IF
            
            ELSE IF (TRIM(bname_aer(jr)) == 'SO4mm') THEN
               so4mm_xtm1(jn) = qxtm1(jl, jk,idt) * mol2mass
               so4mm_xtte(jn) = qxtte(jl, jk,idt) * mol2mass

            ELSE IF (TRIM(bname_aer(jr)) == 'SS') THEN
               ss_xtm1(jn) = qxtm1(jl, jk,idt) * mol2mass
               ss_xtte(jn) = qxtte(jl, jk,idt) * mol2mass
            END IF
         END DO ! naer

         DO jm = 1, nmod
            ! concentration of tsp, om, and bc per mode [ug m-3]
            conctsp(jl, jk, jm) = rho_air(jl) * kg2ug * &
                       (tsp_xtm1(jm) + tsp_xtte(jm) * time_step_len)

            concwsom(jl, jk, jm) = rho_air(jl) * kg2ug * oc2om * &
                         (oc_ws_xtm1(jm) + oc_ws_xtte(jm) * time_step_len)

            concwiom(jl, jk, jm) = rho_air(jl) * kg2ug * oc2om * &
                         (oc_wi_xtm1(jm) + oc_wi_xtte(jm) * time_step_len)

            concbc(jl, jk, jm) = rho_air(jl) * kg2ug * &
                       (bc_xtm1(jm) + bc_xtte(jm) * time_step_len)

            concbc_wgt(jl, jk, jm) = rho_air(jl) * kg2ug * &
                       (bc_xtm1_wgt(jm) + bc_xtte_wgt(jm) * time_step_len)

            concso4mm(jl, jk, jm) = rho_air(jl) * kg2ug * &
                       (so4mm_xtm1(jm) + so4mm_xtte(jm) * time_step_len)

            concss(jl, jk, jm) = rho_air(jl) * kg2ug * &
                       (ss_xtm1(jm) + ss_xtte(jm) * time_step_len)

            ! remove negative values
            conctsp(jl, jk, jm) = MAX(0._dp, conctsp(jl, jk, jm))
            concwsom(jl, jk, jm) = MAX(0._dp, concwsom(jl, jk, jm))
            concwiom(jl, jk, jm) = MAX(0._dp, concwiom(jl, jk, jm))
            concbc(jl, jk, jm) = MAX(0._dp, concbc(jl, jk, jm))
            concbc_wgt(jl, jk, jm) = MAX(0._dp, concbc_wgt(jl, jk, jm))
         END DO ! nmod

         ! total concentrations [ug m-3]
         tot_conctsp(jl, jk) = SUM(conctsp(jl, jk, 1: nmod))
         tot_concwsom(jl, jk) = SUM(concwsom(jl, jk, 1: nmod))
         tot_concwiom(jl, jk) = SUM(concwiom(jl, jk, 1: nmod))
         tot_concso4mm(jl, jk) = SUM(concso4mm(jl, jk, 1: nmod))
         tot_concss(jl, jk) = SUM(concss(jl, jk, 1: nmod))

         ! total bc conc is weighted average of bc multiplied by 
         ! nmodes of bc (4)
         tot_concbc(jl, jk) = 4._dp * SUM(concbc_wgt(jl, jk, 1: nmod))

      END DO ! kproma
   END DO ! nlev

   DO jt = 1, nsvoc
      IF (.NOT. l_mode_partition) THEN
         kp_p => kp(jt)%ptr
         part_p => part(jt)%ptr

         CALL svoc_partition(kproma, nlev, &
                             param_part, &
                             KAHENRY(jt), KBHENRY(jt), &
                             RVAPP(jt), RHVAP(jt), &
                             RKOAB(jt), RKOAM(jt), &
                             RHABC(jt), & !! mz_jw_20170220
                             RPPLFER(jt, 1: ncoeff), &
                             tm1_3d(1: kproma, 1:nlev,jrow) + &
                             tte_3d(1: kproma, 1:nlev,jrow) * &
                             time_step_len, &
                             rhum_3d(1: kproma, 1:nlev,jrow), &
                             tot_conctsp(1: kproma, 1: nlev), &
                             tot_concwsom(1: kproma, 1: nlev), &
                             tot_concwiom(1: kproma, 1: nlev), &
                             tot_concbc(1: kproma, 1: nlev), &
                             tot_concso4mm(1: kproma, 1: nlev), &
                             tot_concss(1: kproma, 1: nlev), &
                             rtot_area(1: kproma, 1: nlev), &
                             kp_p(1: kproma, 1:nlev,jrow), &
                             part_p(1: kproma, 1:nlev,jrow))

         !!! redistribute gas and particle
         idt = idt_svocg(jt)  ! gas
         idt2 = idt_svocp(jt) ! particle

         ! total tracer at t-dt
         zxtm1tot(:, :, jt) = qxtm1(:, :,idt) + qxtm1(:, :,idt2)

         ! total tracer tendency
         zxttetot(:, :, jt) = qxtte(:, :,idt) + qxtte(:, :,idt2)

         DO jk = 1, nlev
            ! gas
            qxtm1(1: kproma, jk,idt) = zxtm1tot(1: kproma, jk, jt) * &
                                 (1._dp - part_p(1: kproma, jk,jrow))

            qxtte(1: kproma, jk,idt) = zxttetot(1: kproma, jk, jt) * &
                                 (1._dp - part_p(1: kproma, jk,jrow))

            ! particle
            qxtm1(1: kproma, jk,idt2) = zxtm1tot(1: kproma, jk, jt) * &
                                 part_p(1: kproma, jk,jrow)

            qxtte(1: kproma, jk,idt2) = zxttetot(1: kproma, jk, jt) * &
                                 part_p(1: kproma, jk,jrow)
         END DO

      ELSE ! partition by mode
         idt = idt_svocg(jt)

         zxtm1tot(:, :, jt) = qxtm1(:, :,idt)
         zxttetot(:, :, jt) = qxtte(:, :,idt)

         kptot(:, :) = 0._dp
         tsptot(:, :) = 0._dp

         DO jm = 1, nmod
            kp_p => kp_m(jt, jm)%ptr
            part_p => part_m(jt, jm)%ptr

            CALL svoc_partition(kproma, nlev, &
                                param_part, &
                                KAHENRY(jt), KBHENRY(jt), &
                                RVAPP(jt), RHVAP(jt), &
                                RKOAB(jt), RKOAM(jt), &
                                RHABC(jt), & !! mz_jw_20170220
                                RPPLFER(jt, 1: ncoeff), &
                                tm1_3d(1: kproma, 1:nlev,jrow) + &
                                tte_3d(1: kproma, 1:nlev,jrow) * &
                                time_step_len, &
                                rhum_3d(1: kproma, 1:nlev,jrow), &
                                conctsp(1: kproma, 1: nlev, jm), &
                                concwsom(1: kproma, 1: nlev, jm), &
                                concwiom(1: kproma, 1: nlev, jm), &
                                concbc(1: kproma, 1: nlev, jm), &
                                concso4mm(1: kproma, 1: nlev, jm), &
                                concss(1: kproma, 1: nlev, jm), &
                                rarea(1: kproma, 1: nlev, jm), &
                                kp_p(1: kproma, 1:nlev,jrow), &
                                part_p(1: kproma, 1:nlev,jrow))

            !!! redistribute particle by mode
            idt2 = idt_svocm(jt, jm)

            DO jk = 1, nlev
               DO jl = 1, kproma

                  IF ((qxtm1(jl, jk,idt) + &
                      qxtte(jl, jk,idt) * time_step_len + &
                      qxtm1(jl, jk,idt2) + &
                      qxtte(jl, jk,idt2) * time_step_len) > 1.e-25_dp .OR. &
                      lstart) THEN

                     zxtm1tot(jl, jk, jt) = zxtm1tot(jl, jk, jt) + qxtm1(jl, jk,idt2)
                     zxttetot(jl, jk, jt) = zxttetot(jl, jk, jt) + qxtte(jl, jk,idt2)

                     qxtm1(jl, jk,idt2) = &
                           (qxtm1(jl, jk,idt) + qxtm1(jl, jk,idt2)) * &
                           part_p(jl, jk,jrow)

                     qxtte(jl, jk,idt2) = &
                           (qxtte(jl, jk,idt) + qxtte(jl, jk,idt2)) * &
                           part_p(jl, jk,jrow)

                     kptot(jl, jk) = kptot(jl, jk) + &
                         part_p(jl, jk,jrow) / (1 - part_p(jl, jk,jrow))

                     tsptot(jl, jk) = tsptot(jl, jk) + conctsp(jl, jk, jm)
                  END IF

               END DO ! kproma
            END DO ! nlev
         END DO  ! nmod

         !!! redistribute gas using total particle fraction
         kp_p => kp(jt)%ptr
         part_p => part(jt)%ptr

         DO jk = 1, nlev
            DO jl = 1, kproma

               IF (zxtm1tot(jl, jk, jt) + &
                  zxttetot(jl, jk, jt) * time_step_len > 1.e-25_dp .OR. &
                  lstart) THEN

                  part_p(jl, jk,jrow) = kptot(jl, jk) / (1 + kptot(jl, jk))

                  qxtm1(jl, jk,idt) = zxtm1tot(jl, jk, jt) * &
                                       (1._dp - part_p(jl, jk,jrow))

                  qxtte(jl, jk,idt) = zxttetot(jl, jk, jt) * &
                                       (1._dp - part_p(jl, jk,jrow))

                  IF (tsptot(jl, jk) > 0._dp) &
                     kp_p(jl, jk,jrow) = kptot(jl, jk) / tsptot(jl, jk)
               END IF

            END DO
         END DO

      END IF ! l_mode_partition
   END DO  ! nsvoc

   !!!
   !!! update atmospheric burden [kg(trac) m-2]
   !!!

   ! layer thickness [m]
   zdz(:, 1: nlev) = deltaz(1: kproma, 1:nlev,jrow)

   DO jt = 1, nsvoc
      idt = idt_svocg(jt)  ! gas

      burd_p => burd(jt)%ptr
      abur(:) = 0._dp

      DO jl = 1, kproma
         DO jk = 1, nlev
            spc_hum = qm1_3d(jl, jk,jrow) + &
                     qte_3d(jl, jk,jrow) * time_step_len

            temp = tm1_3d(jl, jk,jrow) + &
                   tte_3d(jl, jk,jrow) * time_step_len

            ! moist air density [kg m-3]
            rho_air(jl) = press_3d(jl, jk,jrow) / &
                         (temp * rd * (1._dp + vtmpc1 * spc_hum))

            IF (qxtm1(jl, jk,idt) + &
                qxtte(jl, jk,idt) * time_step_len < 0._dp) &
               qxtte(jl, jk,idt) = -qxtm1(jl, jk,idt) / time_step_len

            xtp1_g(jl) = qxtm1(jl, jk,idt) + &
                         qxtte(jl, jk,idt) * time_step_len

            IF (.NOT. l_mode_partition) THEN
               idt2 = idt_svocp(jt) ! particle (bulk)

               IF (qxtm1(jl, jk,idt2) + &
                   qxtte(jl, jk,idt2) * time_step_len < 0._dp) &
                  qxtte(jl, jk,idt2) = -qxtm1(jl, jk,idt2) / time_step_len

               xtp1_p(jl) = qxtm1(jl, jk,idt2) + &
                            qxtte(jl, jk,idt2) * time_step_len
            ELSE
               xtp1_p(jl) = 0._dp

               DO jm = 1, nmod
                  idt2 = idt_svocm(jt, jm) ! mode-specific particle

                  IF (qxtm1(jl, jk,idt2) + &
                      qxtte(jl, jk,idt2) * time_step_len < 0._dp) &
                     qxtte(jl, jk,idt2) = -qxtm1(jl, jk,idt2) / time_step_len

                  xtp1_p(jl) = xtp1_p(jl) + qxtm1(jl, jk,idt2) + &
                               qxtte(jl, jk,idt2) * time_step_len
               END DO
            END IF

            abur(jl) = abur(jl) + (xtp1_g(jl) + xtp1_p(jl)) * &
                       rho_air(jl) * zdz(jl, jk) * MOLMASS(jt) / M_air
         END DO
      END DO

      jk = nlev - 5  ! atmosphere
      burd_p(1: kproma, jk,jrow) = abur(1: kproma) 
   END DO

CASE DEFAULT
   ! should never reach here
END SELECT

CALL end_message_bi(modstr, 'SVOC PHYSC ROUTINES', substr)

END SUBROUTINE svoc_physc

! *************************************************************************

SUBROUTINE svoc_free_memory

IMPLICIT NONE

DEALLOCATE(part_m)
DEALLOCATE(part)
DEALLOCATE(burd)
DEALLOCATE(vola)
DEALLOCATE(depo)
DEALLOCATE(degr)
DEALLOCATE(cwat)
DEALLOCATE(appl)
DEALLOCATE(in_emiss)
DEALLOCATE(vddepaer)
DEALLOCATE(vddepgas)
DEALLOCATE(flux_airsea)
DEALLOCATE(wetflx_ls)
DEALLOCATE(wetflx_cv)
DEALLOCATE(wetflx_aer_tot)
DEALLOCATE(ttescav_bc)
DEALLOCATE(ttescav_bc_cv)
DEALLOCATE(ttescav_bc_ls)
DEALLOCATE(ttescav_bc_ev)
DEALLOCATE(ttescav_m_ls)

END SUBROUTINE svoc_free_memory

! *************************************************************************
! PRIVATE SUBROUTINES
! *************************************************************************

SUBROUTINE svoc_read_nml_cpl(status, iou)

!!! to read namelist for coupling to ECHAM5 and other submodels

USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

IMPLICIT NONE

! I/O
INTEGER, INTENT(OUT) :: status  ! error status
INTEGER, INTENT(IN)  :: iou     ! I/O unit

NAMELIST /CPL/ l_mode_partition, &
               param_part, &
               param_soilv, &
               param_bapo3, &
               aermod_str, &
               imp_ss_as_flux, &
               imp_om_soil, &
               imp_rho_soil, &
               imp_mld_oce, &
               SVOC_NAME, &
               EMISS_IN, &
               MOLMASS, &
               MOLVOL, &
               RKSOIL, &
               RKOCEAN, &
               RWSOL, &
               RVAPP, &
               RHSOL, &
               RHVAP, &
               RHSUB, &
               RLOGKOW, &
               RKOAM, &
               RKOAB, &
               KAHENRY, &
               KBHENRY, &
               RLOSS, &
               RSPRAY, &
               RDENS, &
               REMISP, &
               REMIS_MOD, &
               RPPLFER, &
               RHABC !! mz_jw_20170220-

! local
CHARACTER(LEN = *), PARAMETER :: substr = 'svoc_read_nml_cpl'
LOGICAL :: lex   ! check if file exists
INTEGER :: fstat ! file status
INTEGER :: jt    ! counter

status = 1 ! error

CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
IF (.NOT. lex) RETURN ! <modstr>.nml does not exist

READ(iou, NML = CPL, IOSTAT = fstat)
CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
IF (fstat /= 0) RETURN ! error while reading namelist

CALL read_nml_close(substr, iou, modstr)

status = 0 ! no error

END SUBROUTINE svoc_read_nml_cpl

! *************************************************************************
#endif
! :op_pj_20190326

END MODULE messy_svoc_si
