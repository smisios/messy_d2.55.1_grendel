#include "messy_main_ppd_bi.inc"

!---------------------------------------------------------------------------
! INFOS:
! - aerosol "phase":
!   - flt_phase = cumulative real value from stream:
!     3.0 = for ice (and also nat, liq)
!     2.0 = for nat (and also liq)
!     1.0 = for liquid
!     0.0 = outside of psc region
!   - phase = same, but integer and only for jrow (2D-field)
! - the PSC region:
!   - flt_pscreg = real value (0.0 or 1.0)
!   - val_psc = same as flt_pscreg but a logical (F or T)
! - further infos are labeled with "RS:"
!---------------------------------------------------------------------------

!============================================================================
! DESCRIPTION:
! This module is the interface (called by messy_main_control) between the 
! ECHAM5 base model and the msbm module.
!
! AUTHORS:
! S. Meilinger, Max Planck Institute for Chemistry, Mainz, Germany
!   (2002: original code)
! J. Buchholz, Max Planck Institute for Chemistry, Mainz, Germany
!   (2003, 2004, 2005: psc sedimentation, vectorisation, major revisions)
! O. Kirner, Karlsruhe Institute of Technology (KIT), Germany
!   (2010, kinetic nat parameterisation)
!
! LAST CHANGES:
! 25 February 2005
!    by J. Buchholz, Max-Planck-Institute for Chemistry, Mainz, Germany
! 27 April 2005
!    by P. Joeckel, Max-Planck-Institute for Chemistry, Mainz, Germany
!    - adapted to channel objects
!    - khet1-khet30 as PTR_3D_ARRAY
! 25 November 2005
!    by B. Steil, Max-Planck-Institute for Chemistry, Mainz, Germany
! 02 May 2007
!    by P. Joeckel, Max-Planck-Institute for Chemistry, Mainz, Germany
!    - double code (FORALL/WHERE vs. DO/IF) removed
!      NOTE: according to R. Hatzky, RZG, there was a F90-standard inconform
!            coding at the position WHERE(lCalcChemArray) ...; the standard
!            requires that all fields within the WHERE-construct need to
!            have the same shape as the logical array in the WHERE statement;
!            this was not the case, since the pointers were not associated
!            in case lCalcChem = .FALSE.
! 10 November 2009
!    by Ch. Bruehl, MPI-C, Mainz
!    - bug fixes for calculation of lower boundary of psc 
!      region in messy_psc_e5.f90 (wrong units in call of
!      T_ice and T_nat corrected)
! 10 November 2009
!    by R. Deckert, DLR-Oberpfaffenhofen, Wessling, Germany
!    PSC offline mode (namelist nml/psc.nml)
!    - l_feedback = .true. for simulations with dynamical-chemical 
!      feedback (online mode, default) 
!    - l_feedback = .false. for no dynamical-chemical feedback (offline mode)
!    - Offline mode accomplished using pre-defined climatological 
!      predef_HNO3_tot (= HNO3_gas + HNO3_liq + HNO3_nat) field
!      from a separate online simulation 
!      (namelist nml/import_grid/psc_offline.nml). 
!      Note that HNO3 (= HNO3_gas + HNO3_liq) and HNO3_nat are usually
!      standard output in tracer_gp channel.
!
!         - pre-defined fields are used in the calculation of psc region
!           (flt_pscreg) 
!         - " " " used during water repartitioning/sedimentation. 
!           implementation mainly with additional loop over criterion 
!             tloop==1 (online mode) 
!             tloop==2 (offline mode) 
!         - additional offline NAT as clim_HNO3_nat channel object
!         - special treatment of r_solid and flt_phase (see below)
!    - PSC-related chemical/physical changes in HNO3 and HNO3_nat are 
!      left un-touched as 
!      these do not feedback on dynamics (if settings in msbm.nml, h2o.nml, 
!      rad.nml etc. represent QCTM mode) 
! 12 November 2009
!    by K.-D. Gottschaldt, DLR-Oberpfaffenhofen, Wessling, Germany
!    additional diagnostic output
! 15 February 2017
!    MESSYTENDENCY & MECCA_TAG by P. Joeckel, DLR-IPA
!---------------------------------------------------------------------------

MODULE messy_msbm_si

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,   ONLY: start_message_bi, end_message_bi, info_bi
  ! MESSy
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY        !   (! ka_ok_20100118)
  USE messy_msbm,               ONLY: modstr, dp, IHS_MAX, NSB ! ka_ok_20100118)
  ! mz_pj_20070209+
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL &
                                    , t_chaobj_cpl

  ! mz_pj_20070209-
  ! op_pj_20170215+
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
    mtend_get_start_l,      & ! op_pj_20171127
    mtend_add_l,            &
    mtend_register,         &
    mtend_id_tracer,        &
    mtend_id_q, mtend_id_xl, mtend_id_xi
#endif
  ! op_pj_20170215-
  ! op_pj_20170215+
#ifdef MECCA_TAG
  ! this code is included only when tagging is used
  ! mecca_tag routine to process tendencies for related tagging tracers
  USE messy_mecca_tag_si,       ONLY: mecca_tag_calc_xtte4scav
#endif
  ! op_pj_20170215-

  IMPLICIT NONE

  !-----------------------------------------------------------------
  ! Everything is PRIVATE, except when explicitely stated otherwise:
  !-----------------------------------------------------------------
  PRIVATE

  INTRINSIC NULL

  !-----
  ! 3D-fields (from channel "msbm") of heterogeneous reaction rates
  !-----
  INTEGER, PARAMETER                              :: NKHET = 36
  TYPE(PTR_3D_ARRAY), DIMENSION(NKHET), SAVE      :: khet
  CHARACTER(LEN=34),  DIMENSION(NKHET), PARAMETER :: LONGNAME = (/ &
       ! NAT:
       'khet ClNO3 + HCl on nat           ', &  ! khet(1)
       'khet ClNO3 + H2O on nat           ', &  ! khet(2)
       'khet HOCl + HCl on nat            ', &  ! khet(3)
       'khet N2O5 + HCl on nat            ', &  ! khet(4)
       'khet N2O5 + H2O on nat            ', &  ! khet(5)
       'khet ClNO3 + HBr on nat           ', &  ! khet(6)
       'khet BrNO3 + HCl on nat           ', &  ! khet(7)
       'khet HOCl + HBr on nat            ', &  ! khet(8)
       'khet HOBr + HCl on nat            ', &  ! khet(9)
       'khet HOBr + HBr on nat            ', &  ! khet(10)
       'khet BrNO3 + H2O on nat           ', &  ! khet(11)
       ! ICE:
       'khet ClNO3 + HCl on ice           ', &  ! khet(12)
       'khet ClNO3 + H2O on ice           ', &  ! khet(13)
       'khet HOCl + HCl on ice            ', &  ! khet(14)
       'khet N2O5 + HCl on ice            ', &  ! khet(15)
       'khet N2O5 + H2O on ice            ', &  ! khet(16)
       'khet ClNO3 + HBr on ice           ', &  ! khet(17)
       'khet BrNO3 + HCl on ice           ', &  ! khet(18)
       'khet HOCl + HBr on ice            ', &  ! khet(19)
       'khet HOBr + HCl on ice            ', &  ! khet(20)
       'khet HOBr + HBr on ice            ', &  ! khet(21)
       'khet BrNO3 + H2O on ice           ', &  ! khet(22)
       ! LIQ:
       'khet HOCl + HCl on liquid aerosol ', &  ! khet(23)
       'khet ClNO3 + HCl on liquid aerosol', &  ! khet(24)
       'khet ClNO3 + H2O on liquid aerosol', &  ! khet(25)
       'khet N2O5 + H2O on liquid aerosol ', &  ! khet(26)
       'khet HOBr + HCl on liquid aerosol ', &  ! khet(27)
       'khet HOBr + HBr on liquid aerosol ', &  ! khet(28)
       'khet HOCl + HBr on liquid aerosol ', &  ! khet(29)
       'khet BrNO3 + H2O on liquid aerosol', &  ! khet(30)
       ! Hg (NAT, ICE, LIQ):
       'khet Hg    + H2O on nat           ', &  ! khet(31)
       'khet Hg    + H2O on ice           ', &  ! khet(32)
       'khet Hg    + H2O on liquid aerosol', &  ! khet(33)
       'khet RGM   + H2O on nat           ', &  ! khet(34)
       'khet RGM   + H2O on ice           ', &  ! khet(35)
       'khet RGM   + H2O on liquid aerosol'  /) ! khet(36)

  !-----
  ! 3D-fields (from channel "msbm") for amount-of-substance ratios of PSC 
  ! relevant species HX in air / (mol/mol)
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: &
    H2SO4    => NULL(), &
    HNO3_liq => NULL(), &
    HNO3_nat => NULL(), &
    HNO3_gas => NULL(), &
    HCl_gas  => NULL(), &
    HCl_liq  => NULL(), &
    HOCl_gas => NULL(), &
    HOCl_liq => NULL(), &
    HBr_gas  => NULL(), &
    HBr_liq  => NULL(), &
    HOBr_gas => NULL(), &
    HOBr_liq => NULL()

  ! ka_ok_20100118+
  TYPE (PTR_3D_ARRAY), DIMENSION(NSB), SAVE :: HNO3_natsize
  ! ka_ok_20100118-

  !-----
  ! 2d-field for tropopause index (details via namelist)
  !-----
  REAL(dp), DIMENSION(:,:), POINTER :: tp_i => NULL()
  !-----
  ! 3D-field (from channel "msbm") of stratosphere region indicators
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: flt_stratreg => NULL()
  !-----
  ! 3D-field (from channel "msbm") of PSC region indicators
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: flt_pscreg => NULL()
  !-----
  ! 3D-fields (from channel "msbm") of PSC phase indicator
  ! flt_phase = 3.0 for ice, 2.0 for nat, 1.0 for liquid,
  !             and 0.0 outside of psc relevant region
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: flt_phase => NULL()
  !-----
  ! 3D-fields (from channel "msbm") for number density 
  ! of solid aerosol particles in air / (1/m**3)
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: N_solid => NULL()
  !-----
  ! 3D-fields (from channel "msbm") for solid particle radii
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: r_solid => NULL()
  !-----
  ! 3D-fields (from channel "msbm") for sedimentation velocity
  !-----
! ka_ok_20100118+
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: v_sed => NULL()
  ! if (KINPAR) only for ice
  REAL(dp), DIMENSION(:,:,:), POINTER :: v_sed_ice => NULL()
  ! number density of ice particles in air / (1/m**3)
  REAL(dp), DIMENSION(:,:,:), POINTER :: N_ice => NULL()
  ! ice particle radii
  REAL(dp), DIMENSION(:,:,:), POINTER :: r_ice => NULL()
  ! number density of NAT particles in air / (1/m**3)
  REAL(dp), DIMENSION(:,:,:), POINTER :: N_NAT => NULL()
  ! NAT particle radii
  REAL(dp), DIMENSION(:,:,:), POINTER :: r_NAT => NULL()
  ! mean sedimentation velocity of NAT particles (only if KINPAR=.true.)
  REAL(dp), DIMENSION(:,:,:), POINTER :: v_sed_NAT => NULL()

  REAL(dp), DIMENSION(:,:,:), POINTER :: HNO3_Change_Sed => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: SedStep_Ice_noarr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: Ice_Change_Sed_noarr => NULL()
! ka_ok_20100118-
  !-----
  ! 3D-fields (from channel "msbm") for liquid parameterisation validity flag
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: flt_val => NULL()
  !-----
  ! 3D-fields (from channel "msbm") for liquid aerosol surface area density
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: A_liq => NULL()
  !-----
  ! 3D-fields (from channel "msbm") for liquid aerosol surface median radius
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: r_SurfMed => NULL()

  ! op_rd_20100108+
  ! switch for offline mode
  LOGICAL :: l_feedback = .TRUE.
  ! additional channel objects for offline- mode (l_feedback=.FALSE.):
  ! HNO3 climatology imported via import_grid
  REAL(dp), DIMENSION(:,:,:), POINTER :: predef_HNO3_tot => NULL() ! [mol/mol]
  ! NAT calculated from offline HNO3 climatology
  REAL(dp), DIMENSION(:,:,:), POINTER :: clim_HNO3_nat   => NULL() ! [mol/mol]
  ! decide in loops, which HNO3_tot to be used
  REAL(dp), DIMENSION(:,:), POINTER :: tmp_HNO3_tot => NULL()
  ! op_rd_20100108-

  !op_kg_20091112+
  !-----
  ! diagnostic 3D-fields for ice loss due to sedimentation
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: loss_nat => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: prod_nat => NULL()
  !op_kg_20091112- 

  !-----
  ! IDs for gas phase tracers
  !-----
  INTEGER :: idt_HNO3=0,                          &
             idt_HCl=0,  idt_HOCl=0, idt_ClNO3=0, &
             idt_HBr=0,  idt_HOBr=0, idt_BrNO3=0
  !-----
  !ID for particle phase HNO3
  !-----
  INTEGER :: idt_HNO3_nat=0
  INTEGER, DIMENSION(NSB) :: idt_HNO3_natsize ! ka_ok_20100118
  !-----
  ! switch for tracer initialisation
  !-----
!!$  LOGICAL :: tracer_init_required = .false.

  ! mz_pj_20070209+
  CHARACTER(LEN=STRLEN_CHANNEL) :: Tropop_Channel = ''
  CHARACTER(LEN=STRLEN_CHANNEL) :: Tropop_Index = ''
  ! mz_pj_20070209-

  ! CPL namelist
  ! - name of imported H2SO4 climatology
  TYPE(t_chaobj_cpl), SAVE :: c_H2SO4clim
  ! - name of predefined total HNO3 (required, if l_feedback = .FLASE.)
  TYPE(t_chaobj_cpl), SAVE :: c_predef_HNO3_tot

  !-----
  ! switch for calculation of reaction rates
  !-----
  LOGICAL, SAVE :: LCalcChem = .TRUE. ! mz_bs_20050713: SAVE added

  !-----
  ! temperature shift for sensitivity studies
  !-----
  REAL(dp) :: TempShift = 0.0_dp
  ! op_pj_20091013+
  ! lower boundary of PSC region [deg N] (SH, NH)
  REAL(dp), DIMENSION(2) :: r_lat = (/ -55.0_dp, 45.0_dp /)
  ! lower boundary of PSC region [Pa] (SH, NH)
  REAL(dp), DIMENSION(2) :: r_lb = (/ 18000.0_dp, 18000.0_dp /)
  ! middle boundary of PSC region [Pa] (SH, NH)
  REAL(dp), DIMENSION(2) :: r_mb = (/ 14000.0_dp, 10000.0_dp /)
  ! upper boundary of PSC region [PA] (SH, NH)
  REAL(dp), DIMENSION(2) :: r_ub = (/ 2000.0_dp, 2000.0_dp /)
  ! op_pj_20091013-

  ! mz_rs_20060123+
  TYPE(PTR_3D_ARRAY), PUBLIC, DIMENSION(IHS_MAX), SAVE :: khet_St_3d
  ! mz_rs_20060123-

  ! op_pj_20170215+
#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif
  ! op_pj_20170215-

  PUBLIC :: LCalcChem, TempShift
  PUBLIC :: msbm_initialize    ! checks the phase of strat. particles
  PUBLIC :: msbm_new_tracer    ! define PSC-specific tracers
  PUBLIC :: msbm_init_memory   ! allocate memory
  PUBLIC :: msbm_init_coupling ! mz_pj_20070131
!!$  PUBLIC :: msbm_init_tracer   ! initialize PSC-specific tracers
  PUBLIC :: msbm_local_start   ! determine psc relevant region (flt_pscreg)
  PUBLIC :: msbm_physc         ! program to calculate composition and
                              ! het. reaction rates of PSCs
  PUBLIC :: msbm_local_end

CONTAINS

!=============================================================================

  SUBROUTINE  msbm_initialize

    ! BMIL/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,   ONLY: error_bi
    ! MESSy
    USE messy_main_tools,        ONLY: find_next_free_unit
    USE messy_msbm,              ONLY: msbm_read_nml_ctrl, &
                                       LAdvectIceNat, LHomNucNAT, &
                                       NatFormThreshold, &
                                       minKhet, maxKhet, &
                                       SupSatIce, &
                                       r_min, N_solid_max, SedScheme, &
                                       KinPar ! ka_ok_20100118


    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: substr = 'msbm_initialize'
    INTEGER                     :: status, iou

    !-----
    ! read namelist CTRL
    !-----
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL msbm_read_nml_ctrl(status, iou)
      IF (status /= 0) CALL error_bi('call msbm_read_nml_ctrl failed',substr)
    END IF
    CALL p_bcast(LAdvectIceNat, p_io)
    CALL p_bcast(LHomNucNAT, p_io)
    CALL p_bcast(NatFormThreshold, p_io)
    CALL p_bcast(minKhet, p_io)
    CALL p_bcast(maxKhet, p_io)
    CALL p_bcast(SupSatIce, p_io)
    CALL p_bcast(r_min, p_io)
    CALL p_bcast(N_solid_max, p_io)
    CALL p_bcast(SedScheme, p_io)
    CALL p_bcast(KinPar, p_io) ! ka_ok_20100118
    
    !-----
    ! read namelist CPL
    !-----
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL msbm_read_nml_cpl(status, iou)   ! read /CPL/
      IF (status /= 0) CALL error_bi('call msbm_read_nml_cpl failed',substr)
    END IF
    CALL p_bcast(LCalcChem,p_io)
    CALL p_bcast(TempShift,p_io)
    ! mz_pj_20070209+
    CALL p_bcast(Tropop_Channel,p_io)
    CALL p_bcast(Tropop_Index,p_io)
    ! mz_pj_20070209-
    CALL p_bcast(c_H2SO4clim%CHA,p_io)
    CALL p_bcast(c_H2SO4clim%OBJ,p_io)
    CALL p_bcast(c_predef_HNO3_tot%CHA,p_io)
    CALL p_bcast(c_predef_HNO3_tot%OBJ,p_io)
    CALL p_bcast(r_lat ,p_io)      ! op_pj_20091013
    CALL p_bcast(r_lb, p_io)       ! op_pj_20091013
    CALL p_bcast(r_mb, p_io)       ! op_pj_20091013
    CALL p_bcast(r_ub, p_io)       ! ka_ok_20100129
    CALL p_bcast(l_feedback,p_io)  ! op_rd_20100108

  END SUBROUTINE msbm_initialize

!=============================================================================

  ! mz_pj_20070131+
  SUBROUTINE msbm_init_coupling

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi &
                                         , error_bi
    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_info, get_channel_object

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'msbm_init_coupling'
    INTEGER :: status

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

    CALL get_channel_info(status, 'psc')
    IF (status == 0) CALL error_bi( &
         'PSC and MSBM cannot be applied simultaneously', substr)

    ! mz_pj_20070209+
    ! GET POINTER TO TROPOPAUSE INDEX
    IF (p_parallel_io) THEN
       WRITE(*,*) '  INITIALIZING TROPOPAUSE LEVEL INDEX ...'
    END IF
    !
    CALL get_channel_object(status &
         , TRIM(Tropop_Channel), TRIM(Tropop_Index), p2=tp_i)
    CALL channel_halt(substr, status)
    !
    IF (p_parallel_io) THEN
       WRITE(*,*) '  ... OK: ',&
            TRIM(Tropop_Channel)//' - '// TRIM(Tropop_Index)
    END IF
    ! mz_pj_20070209-

    IF (.NOT. l_feedback) THEN
       CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output
       CALL get_channel_object(status, TRIM(c_predef_HNO3_tot%CHA) &
            , TRIM(c_predef_HNO3_tot%OBJ), p3=predef_HNO3_tot)
       IF (status /= 0) &
            CALL error_bi( &
            'channel-object for HNO3(total) climatology not available' &
            , substr)
    END IF

    CALL get_channel_object(status, TRIM(c_H2SO4clim%CHA) &
         , TRIM(c_H2SO4clim%OBJ) &
         , p3=H2SO4)
    IF (status /= 0) &
         CALL error_bi('channel-object for H2SO4 climatology not available' &
         , substr)

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE msbm_init_coupling
  ! mz_pj_20070131-

!=============================================================================

  SUBROUTINE msbm_new_tracer

  !-----------------------------------------------------------------------------
  !This subroutine defines psc-specific tracers
  !-----------------------------------------------------------------------------
    ! ECHAM5/MESSy
    USE messy_main_constants_mem,   ONLY: MBr, MCl, MH, MN, MO
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, get_tracer,       &
                                        AIR,AEROSOL,AMOUNTFRACTION,   &
                                        ON, OFF, I_SCAV,              &
                                        I_DRYDEP, R_molarmass, R_pss, &
                                        set_tracer
    USE messy_msbm,               ONLY: KinPar, NSB ! ka_ok_20100118

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='msbm_new_tracer'
    INTEGER :: i_err
    INTEGER :: status
    ! ka_ok_20100118+
    INTEGER          :: jt
    CHARACTER(LEN=1) :: istr = ' '
    ! ka_ok_20100118-

    CALL start_message_bi(modstr, 'defining psc-relevant tracers', substr)

    CALL get_tracer(i_err, GPTRSTR, 'HNO3', idx=idt_HNO3)
    IF (i_err.ne.0) THEN

       CALL new_tracer(status, GPTRSTR, &
            'HNO3','msbm',idx=idt_HNO3, &
            unit='mol/mol', &
            longname='HNO3 molar mixing ratio', &
            medium=AIR, &
            quantity=AMOUNTFRACTION)
       CALL tracer_halt(substr,status)

      CALL set_tracer(status, GPTRSTR, idt_HNO3, I_DRYDEP,  ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HNO3, I_SCAV,    ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HNO3, R_molarmass, MH+MN+MO*3.0_dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HNO3, R_pss,     1.0e4_dp)
      CALL tracer_halt(substr, status)

      CALL tracer_halt(substr,status)
    END IF

    CALL get_tracer(i_err, GPTRSTR, 'HNO3', subname='nat', idx=idt_HNO3_nat)
    IF (i_err.ne.0) THEN

      CALL new_tracer(status, GPTRSTR, &
        'HNO3','msbm',subname='nat',idx=idt_HNO3_nat, &
         unit='mol/mol', &
         longname='NAT phase HNO3 molar mixing ratio', &
         medium=AEROSOL, &
         quantity=AMOUNTFRACTION)
      CALL tracer_halt(substr,status)
       
      CALL set_tracer(status, GPTRSTR, idt_HNO3_nat, I_DRYDEP,  OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HNO3_nat, I_SCAV,    OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HNO3_nat, R_molarmass,MH+MN+MO*3.0_dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HNO3_nat, R_pss,     1.0_dp)
      CALL tracer_halt(substr, status)

    END IF

    CALL get_tracer(i_err, GPTRSTR, 'HCl', idx=idt_HCl)
    IF (i_err.ne.0) THEN

      CALL new_tracer(status, GPTRSTR, &
        'HCl', 'msbm', idx=idt_HCl, &
         unit='mol/mol', &
         longname='HCl molar mixing ratio', &
         medium=AIR, &
         quantity=AMOUNTFRACTION)
      CALL tracer_halt(substr,status)

      CALL set_tracer(status, GPTRSTR, idt_HCL, I_DRYDEP,    ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HCL, I_SCAV,      ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HCL, R_molarmass, MH+MCl)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HCL, R_pss,      1.e14_dp)
      CALL tracer_halt(substr, status)
    END IF

    CALL get_tracer(i_err, GPTRSTR, 'HOCl', idx=idt_HOCl)
    IF (i_err.ne.0) THEN

       CALL new_tracer(status, GPTRSTR, &
        'HOCl','msbm', idx=idt_HOCl, &
         unit='mol/mol', &
         longname='HOCl molar mixing ratio', &
         medium=AIR, &
         quantity=AMOUNTFRACTION)
      CALL tracer_halt(substr,status)

      CALL set_tracer(status, GPTRSTR, idt_HOCL, I_DRYDEP,    ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HOCL, I_SCAV,      ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HOCL, R_molarmass, MH+MO+MCl)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HOCL, R_pss,       6.7e2_dp)
      CALL tracer_halt(substr, status)
  
    END IF

    CALL get_tracer(i_err, GPTRSTR, 'ClNO3', idx=idt_ClNO3)
    IF (i_err.ne.0) THEN

      CALL new_tracer(status, GPTRSTR, &
        'ClNO3','msbm', idx=idt_ClNO3, &
         unit='mol/mol', &
         longname='ClNO3 molar mixing ratio', &
         medium=AIR, &
         quantity=AMOUNTFRACTION)
      CALL tracer_halt(substr,status)
      
      CALL set_tracer(status, GPTRSTR, idt_ClNO3, I_DRYDEP,    ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_ClNO3, I_SCAV,      ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_ClNO3, R_molarmass, MCl+MN+MO*3.0_dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_ClNO3, R_pss,       1.e30_dp)
      CALL tracer_halt(substr, status)

    END IF

    CALL get_tracer(i_err, GPTRSTR, 'HBr', idx=idt_HBr)
    IF (i_err.ne.0) THEN

     CALL new_tracer(status, GPTRSTR, &
        'HBr', 'msbm', idx=idt_HBr, &
         unit='mol/mol', &
         longname='HBr molar mixing ratio', &
         medium=AIR, &
         quantity=AMOUNTFRACTION)
      CALL tracer_halt(substr,status)
      
      CALL set_tracer(status, GPTRSTR, idt_HBr, I_DRYDEP,    ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HBr, I_SCAV,      ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HBr, R_molarmass, MH+MBr)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HBr, R_pss,       1.3e17_dp)
      CALL tracer_halt(substr, status)
    END IF

    CALL get_tracer(i_err, GPTRSTR, 'HOBr', idx=idt_HOBr)
    IF (i_err.ne.0) THEN

      CALL new_tracer(status, GPTRSTR, &
        'HOBr','msbm', idx=idt_HOBr, &
         unit='mol/mol', &
         longname='HOBr molar mixing ratio', &
         medium=AIR, &
         quantity=AMOUNTFRACTION)
      CALL tracer_halt(substr,status)

      CALL set_tracer(status, GPTRSTR, idt_HOBr, I_DRYDEP,    ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HOBr, I_SCAV,      ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HOBr, R_molarmass, MH+MO+MBr)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_HOBr, R_pss,       9.1e1_dp)
      CALL tracer_halt(substr, status)
    END IF

    CALL get_tracer(i_err, GPTRSTR, 'BrNO3', idx=idt_BrNO3)
    IF (i_err.ne.0) THEN

       CALL new_tracer(status, GPTRSTR, &
        'BrNO3','msbm', idx=idt_BrNO3, &
         unit='mol/mol', &
         longname='BrNO3 molar mixing ratio', &
         medium=AIR, &
         quantity=AMOUNTFRACTION)
      CALL tracer_halt(substr,status)

      CALL set_tracer(status, GPTRSTR, idt_BrNO3, I_DRYDEP,    ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_BrNO3, I_SCAV,      ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_BrNO3, R_molarmass, MBr+MN+MO*3.0_dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_BrNO3, R_pss,       1.e30_dp)
      CALL tracer_halt(substr, status)
    END IF

    ! ka_ok_20100118+
    IF (KinPar) THEN
       DO jt=1, NSB
          WRITE(istr,'(i1)') jt
          CALL get_tracer(i_err, GPTRSTR, 'HNO3', subname='natsize'//istr, &
               idx=idt_HNO3_natsize(jt))
          IF (i_err.ne.0) THEN

             CALL new_tracer(status, GPTRSTR, &
                  'HNO3','msbm',subname='natsize'//istr,idx=idt_HNO3_natsize(jt), &
                  unit='mol/mol', &
                  longname='NAT in sizebin '//istr//', HNO3 molar mixing ratio', &
                  medium=AEROSOL, &
                  quantity=AMOUNTFRACTION)
             CALL tracer_halt(substr,status)

             CALL set_tracer(status, GPTRSTR, idt_HNO3_natsize(jt), I_DRYDEP,OFF)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, GPTRSTR, idt_HNO3_natsize(jt), I_SCAV,  OFF)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, GPTRSTR, idt_HNO3_natsize(jt), R_molarmass &
                  , MH+MN+MO*3.0_dp)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, GPTRSTR, idt_HNO3_natsize(jt), R_pss,1.0_dp)
             CALL tracer_halt(substr, status)
          END IF
       END DO
    ENDIF
    ! ka_ok_20100118-

    CALL end_message_bi(modstr, 'defining psc-relevant tracers', substr)

  END SUBROUTINE msbm_new_tracer

!=============================================================================

  SUBROUTINE msbm_init_memory
  !-----------------------------------------------------------------------------
  !This subroutine defines and/or modifies psc-specific output channels:
  !1) 3D-fields of all PSC relevant species (HX=H2O, HNO3, HCl, HOCl, HBr, HOBr)
  !   distinguishing between the different phases (gas, liquid and ice/nat)
  !   UNITS: mixing ratio (mol HX/mol air)
  !2) 3D-field of phase indicator:
  !   flt_phase=1. if only liquid particles present
  !   flt_phase=2. if only liquid and NAT particles present
  !   flt_phase=3. if liquid, NAT and ice particles present
  !   flt_phase=0. outside of psc relevant region
  !3) 3D-field of solid particle radii: r_solid
  !4) 3D-field of flags indicating whether a grid box is in the region where
  !   polar stratospheric cloud are calculated: flt_pscreg
  !5) 3D-field of flags indicating whether the parameterisation for liquid
  !   ternary solutions is valid: flt_val=1.
  !-----------------------------------------------------------------------------
    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    ! MESSY
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_msbm,                  ONLY: khet_St_name, KinPar, NSB

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'msbm_init_memory'
    INTEGER          :: status
    INTEGER          :: i
    CHARACTER(LEN=2) :: istr2
    CHARACTER(LEN=STRLEN_MEDIUM) :: name
    ! ka_ok_20100118+
!!$    REAL(DP), POINTER :: r_kinpar => NULL()
    INTEGER           :: jt
    CHARACTER(LEN=1)  :: istr
    ! ka_ok_20100118-

    CALL start_message_bi(modstr,'defining channel for msbm', substr)
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    ! -----
    ! Define a new output channel
    ! -----
    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    
    ! -----
    ! PSC specific molecule information
    ! -----
    CALL new_channel_object(status, modstr, 'HNO3_liq', p3=HNO3_liq )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HNO3_liq' &
         , 'long_name', c='liquid phase HNO3')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HNO3_liq', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HNO3_nat', p3=HNO3_nat )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HNO3_nat' &
         , 'long_name', c='NAT phase HNO3')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HNO3_nat', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HNO3_gas', p3=HNO3_gas )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HNO3_gas' &
         , 'long_name', c='gas phase HNO3')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HNO3_gas', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HCl_liq', p3=HCl_liq )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HCl_liq' &
         , 'long_name', c='liquid phase HCl')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HCl_liq', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HCl_gas', p3=HCl_gas )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HCl_gas' &
         , 'long_name', c='gas phase HCl')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HCl_gas', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HBr_liq', p3=HBr_liq )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HBr_liq' &
         , 'long_name', c='liquid phase HBr')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HBr_liq', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HBr_gas', p3=HBr_gas )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HBr_gas' &
         , 'long_name', c='gas phase HBr')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HBr_gas', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HOCl_liq', p3=HOCl_liq )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HOCl_liq' &
         , 'long_name', c='liquid phase HOCl')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HOCl_liq', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HOCl_gas', p3=HOCl_gas )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HOCl_gas' &
         , 'long_name', c='gas phase HOCl')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HOCl_gas', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HOBr_liq', p3=HOBr_liq )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HOBr_liq' &
         , 'long_name', c='liquid phase HOBr')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HOBr_liq', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'HOBr_gas', p3=HOBr_gas )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HOBr_gas' &
         , 'long_name', c='gas phase HOBr')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'HOBr_gas', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    !-----
    ! stratosphere region indicator
    !-----
    CALL new_channel_object(status, modstr, 'STRAT_region', p3=flt_stratreg )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'STRAT_region' &
         , 'long_name', c='flag indicating stratosphere region')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'STRAT_region', 'units', c=' ')
    CALL channel_halt(substr, status)

    !-----
    ! psc region indicator
    !-----
    CALL new_channel_object(status, modstr, 'PSC_region', p3=flt_pscreg )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PSC_region' &
         , 'long_name', c='flag indicating PSC relevant region')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PSC_region', 'units', c=' ')
    CALL channel_halt(substr, status)

    !-----
    ! Phase indicators:
    !   flt_phase=1. if only liquid particles present
    !   flt_phase=2. if only liquid and NAT particles present
    !   flt_phase=3. if liquid, NAT and ice particles present
    !   flt_phase=0. outside of psc relevant region
    !-----
    CALL new_channel_object(status, modstr, 'phase', p3=flt_phase &
         , lrestreq=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'phase' &
         , 'long_name', c='liq, liq+nat, liq+nat+ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'phase', 'units', c=' ')
    CALL channel_halt(substr, status)

    !-----
    ! solid particle information
    !-----
    CALL new_channel_object(status, modstr, 'N_solid', p3=N_solid )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'N_solid' &
         , 'long_name', c='solid particle number density')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'N_solid', 'units', c='1/m**3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'r_solid', p3=r_solid )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'r_solid' &
         , 'long_name', c='Radius of ice/NAT particles')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'r_solid', 'units', c='m')
    CALL channel_halt(substr, status)

! ka_ok_20100118+
    CALL new_channel_object(status, modstr, 'v_sed_ice', p3=v_sed_ice )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v_sed_ice' &
         , 'long_name', c='sedimentation velocity of ice particles')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v_sed_ice', 'units', c='m/s')
    CALL channel_halt(substr, status)
! ka_ok_20100118-

    !-----
    ! liquid particle information
    !-----
    CALL new_channel_object(status, modstr, 'val', p3=flt_val )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'val' &
         , 'long_name', c='validity of ternary liq. param.')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'val', 'units', c=' ')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'A_liq', p3=A_liq )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'A_liq' &
         , 'long_name', c='liquid aerosol surface area density')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'A_liq', 'units', c='cm**2/cm**3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'r_SurfMed', p3=r_SurfMed )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'r_SurfMed' &
         , 'long_name', c='liquid aerosol surface median radius')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'r_SurfMed', 'units', c='cm')
    CALL channel_halt(substr, status)

    IF (LCalcChem) THEN
      !-----
      ! heterogeneous reaction rates
      !-----
      DO i=1, NKHET
         istr2 = '  '
         IF (i < 10) THEN
            WRITE(istr2,'(i1)') i
         ELSE
            WRITE(istr2,'(i2)') i
         END IF
         CALL new_channel_object(status, modstr, 'khet'//TRIM(istr2) &
              , p3=khet(i)%ptr )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr,  'khet'//TRIM(istr2) &
              , 'long_name', c=TRIM(LONGNAME(i)) )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr,  'khet'//TRIM(istr2) &
              , 'units', c='cm**3/s' )
      END DO
    END IF

    ! mz_rs_20060123+
    ! stratospheric rate coefficients
    DO i=1, IHS_MAX
      name = 'khet_'//TRIM(khet_St_name(i))
      CALL new_channel_object(status, modstr, name, &
        p3=khet_St_3d(i)%PTR)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, name, &
        'long_name', c=name)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, name, 'units', c='cm**3/s')
      CALL channel_halt(substr, status)
      CALL info_bi('channel/object '//modstr//'/'//TRIM(name)//' was created')
    ENDDO
    ! mz_rs_20060123-

    ! op_rd_20100108+
    IF (.NOT. l_feedback) THEN
       CALL new_channel_object(status, modstr, 'clim_HNO3_nat' &
            , p3 = clim_HNO3_nat)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'clim_HNO3_nat' &
            , 'long_name', c='climatological NAT phase HNO3')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'clim_HNO3_nat', 'units', c='mol/mol')
       CALL channel_halt(substr, status)
    END IF
    ! op_rd_20100108-

    !op_kg_20091112+
    CALL new_channel_object(status, modstr, 'PROD_NAT', p3 = prod_nat)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PROD_NAT' &
         , 'long_name', c='Prod of NAT due to sedimentation ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PROD_NAT', 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LOSS_NAT', p3 = loss_nat)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LOSS_NAT' &
     , 'long_name', c='Loss of NAT due to sedimentation ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LOSS_NAT', 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)
    !op_kg_20091112-

    ! ka_ok_20100118+
    CALL new_channel_object(status, modstr, 'r_ice', p3 = r_ice)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'r_ice' &
         , 'long_name', c='Radius of ice particles')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'r_ice', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'v_sed_NAT', p3 = v_sed_NAT)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v_sed_NAT' &
         , 'long_name', c='medium sedimentation velocity of NAT particles')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v_sed_NAT', 'units', c='m/s')
    CALL channel_halt(substr, status)

    IF (KinPar) THEN

       DO jt=1, NSB
          WRITE(istr, '(i1)') jt

          CALL new_channel_object(status, modstr, 'HNO3_natsize'//istr &
               , p3 = HNO3_natsize(jt)%ptr)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'HNO3_natsize'//istr &
               , 'long_name', c='NAT in sizebin '//istr//' HNO3')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'HNO3_natsize'//istr &
               , 'units', c='mol/mol')
          CALL channel_halt(substr, status)
       END DO

       CALL new_channel_object(status, modstr, 'N_ice', p3 = N_ice)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'N_ice' &
            , 'long_name', c='ice particle number density')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'N_ice', 'units', c='1/m**3')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'N_NAT', p3 = N_NAT)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'N_NAT' &
            , 'long_name', c='number density of NAT particle')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'N_NAT' , 'units', c='1/m**3')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'r_NAT', p3 = r_NAT)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'r_NAT' &
            , 'long_name', c='Radius of NAT particles')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'r_NAT', 'units', c='m')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'HNO3_Change_Sed' &
            , p3 = HNO3_Change_Sed)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'HNO3_Change_Sed' &
            , 'long_name', c='sedimentation change of HNO3 ')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'HNO3_Change_Sed' &
            , 'units', c='mol/mol')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'SedStep_ICE' &
            , p3 = SedStep_ICE_noarr)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'SedStep_ICE' &
            , 'long_name', c='sedimentation step of ICE particles')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'SedStep_ICE' &
            , 'units', c='Pa')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'Ice_Change_Sed' &
            , p3 = Ice_Change_Sed_noarr)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'Ice_Change_Sed' &
            , 'long_name', c='sedimentation change of ICE')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'Ice_Change_Sed' &
            , 'units', c='mol/mol')
       CALL channel_halt(substr, status)
    END IF
    ! ka_ok_20100118-

    ! op_pj_20170215+
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_q)
    CALL mtend_register(my_handle, mtend_id_xl)
    CALL mtend_register(my_handle, mtend_id_xi)
    !
#ifdef MECCA_TAG
    !  tagged tracer IDs a priori not known; thus allow all tracers ...
    CALL mtend_register(my_handle, mtend_id_tracer)
    !
#else
    !

    CALL mtend_register(my_handle, idt_HNO3)
    CALL mtend_register(my_handle, idt_HCl)
    CALL mtend_register(my_handle, idt_HOCl)
    CALL mtend_register(my_handle, idt_ClNO3)
    CALL mtend_register(my_handle, idt_HBr)
    CALL mtend_register(my_handle, idt_HOBr)
    CALL mtend_register(my_handle, idt_BrNO3)
    !
    CALL mtend_register(my_handle, idt_HNO3_nat)
    !
    IF (KinPar) THEN
       CALL mtend_register(my_handle, mtend_id_vec=idt_HNO3_natsize(:))
    END IF
#endif
#endif
    ! op_pj_20170215-

    CALL end_message_bi(modstr,'defining channel for msbm', substr)

  END SUBROUTINE msbm_init_memory

!=============================================================================

!!$SUBROUTINE msbm_init_tracer
!!$
!!$  ! ECHAM5/MESSy
!!$  USE messy_main_tracer_bi, ONLY: tracer_init
!!$
!!$  IMPLICIT NONE
!!$
!!$  IF (tracer_init_required) CALL tracer_init(modstr)
!!$
!!$END SUBROUTINE msbm_init_tracer
!!$
!!$!=============================================================================

!=============================================================================

SUBROUTINE msbm_local_start
!-------------------------------------------------------------------------
! This subroutine checks which grid boxes are relevant for msbm calculation
! and sets the msbm channel object flt_pscreg accordingly.
!-------------------------------------------------------------------------

  ! ECHAM5/MESSy
  USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1, ntrac_gp
  USE messy_main_grid_def_mem_bi,ONLY: nlev, jrow, kproma, nproma
  USE messy_main_grid_def_bi,    ONLY: philat_2d
  USE messy_main_data_bi,       ONLY: qm1_3d => qm1    &
                                    , ledith           & ! op_pj_20180725
#ifdef ECHAM5
                                    , xlm1_3d => xlm1  &
                                    , xim1_3d => xim1  &
#endif
#if defined(COSMO) || defined(MESSYDWARF)
                                    , xlm1_3d  &
                                    , xim1_3d  &
#endif
#ifdef CESM1
                                    , xlm1_3d  &
                                    , xim1_3d  &
#endif
                                    , qte_3d    &
                                    , xlte_3d   &
                                    , xite_3d   &
                                    , press_3d
  USE messy_main_timer,         ONLY: time_step_len  
  USE messy_main_constants_mem, ONLY: M_air, M_H2O ! mz_hr_20120524
  USE messy_msbm,               ONLY: T_nat, T_ice


  IMPLICIT none

  ! mz_hr_20120524+
  REAL(dp), PARAMETER                      :: mmratio = M_air/M_H2O
  ! mz_hr_20120524-
  REAL(dp), DIMENSION(nproma,nlev)         :: H2O_tot, Tice, Tnat
  REAL(dp), DIMENSION(nproma,nlev), TARGET :: HNO3_tot  ! op_pj_20100108
  INTEGER,  DIMENSION(nproma)              :: index_arr ! mz_bs_20060112
  INTEGER :: jk, jp
  INTEGER :: idt, idt2
  LOGICAL :: zl ! op_pj_20180725

  INTRINSIC MAX, TINY, MINLOC, SUM

  flt_pscreg(_RI_XYZ__(1:kproma,jrow,:))= 0.0_dp

  H2O_tot(1:kproma,:) =                            &
    max(( xim1_3d(_RI_XYZ__(1:kproma,jrow,:))                &
         +xite_3d(_RI_XYZ__(1:kproma,jrow,:))*time_step_len  &
         +xlm1_3d(_RI_XYZ__(1:kproma,jrow,:))                &
         +xlte_3d(_RI_XYZ__(1:kproma,jrow,:))*time_step_len  &
         +qm1_3d(_RI_XYZ__(1:kproma,jrow,:))                 &
         +qte_3d(_RI_XYZ__(1:kproma,jrow,:))*time_step_len)  &
        * mmratio, tiny(1.0_dp))

  idt  = idt_HNO3_nat
  idt2 = idt_HNO3
  HNO3_tot(1:kproma,:) =                           &
    max(( pxtm1(_RI_X_ZN_(1:kproma,:,idt))                  &
         +pxtte(_RI_X_ZN_(1:kproma,:,idt)) *time_step_len   &
         +pxtm1(_RI_X_ZN_(1:kproma,:,idt2))                 &
         +pxtte(_RI_X_ZN_(1:kproma,:,idt2))*time_step_len), &
        tiny(1.0_dp))

  !-----
  ! The rest of the routine defines flt_pscreg depending on the latitude, 
  ! the height, and T_nat compared to T_ice:
  !-----

  !-----
  ! Set default for Tice, Tnat so that Tice>Tnat if not explicitely set
  ! otherwise.
  !-----
  Tice = 1.0_dp
  Tnat = 0.0_dp

  !-----
  ! set Tice<Tnat if we are north of r_lat(2) degree northern latitude
  ! and above r_mb(2) hPa
  !-----
  DO jk=1,nlev
     DO jp=1,kproma
        IF( (philat_2d(jp,jrow) > r_lat(2)).AND.&
             (press_3d(_RI_XYZ__(jp,jrow,jk)) < r_mb(2)) ) THEN
           Tice(jp,jk)=0.0_dp
           Tnat(jp,jk)=1.0_dp
        ENDIF
     ENDDO
  ENDDO

  !-----
  ! set Tice<Tnat if we are south of r_lat(1) degree southern latitude
  ! and above r_mb(1) hPa
  !-----
  DO jk=1,nlev
     DO jp=1,kproma
        IF( (philat_2d(jp,jrow) < r_lat(1)) .AND. &
             (press_3d(_RI_XYZ__(jp,jrow,jk)) < r_mb(1)) ) THEN
           Tice(jp,jk)=0.0_dp
           Tnat(jp,jk)=1.0_dp
        ENDIF
     ENDDO
  ENDDO

  !-----
  ! Calculate Tice and Tnat if we are north of r_lat(2) 
  ! northern latitude, below r_mb(2) hPa, and above r_lb(2) hPa
  !-----
  ! op_rd_20100108+
  IF (l_feedback) THEN
     tmp_HNO3_tot => HNO3_tot(:,:)
  ELSE
     tmp_HNO3_tot => predef_HNO3_tot(_RI_XYZ__(:,jrow,:))
  END IF
  ! op_rd_20100108-

  DO jk=1,nlev
     DO jp=1,kproma
        IF( (philat_2d(jp,jrow) >  r_lat(2))  .AND. &
             (press_3d(_RI_XYZ__(jp,jrow,jk)) >= r_mb(2))   .AND. &
             (press_3d(_RI_XYZ__(jp,jrow,jk)) <= r_lb(2)) ) THEN
           Tice(jp,jk) = T_ice(press_3d(_RI_XYZ__(jp,jrow,jk)) * &
! mz_cb_20100108+
! bug fix
!!$                H2O_tot(jp,jk)/1013.0_dp)
                H2O_tot(jp,jk)/100.0_dp)
! mz_cb_20100108-
           Tnat(jp,jk) = T_nat(press_3d(_RI_XYZ__(jp,jrow,jk)) * &
! mz_cb_20100108+
!!$                H2O_tot(jp,jk)/1013.0_dp, &
                H2O_tot(jp,jk)/101325.0_dp, &
! mz_cb_20100108-
! op_rd_20100108+
! mz_cb_20100108+
                press_3d(_RI_XYZ__(jp,jrow,jk)) * tmp_HNO3_tot(jp,jk)/101325.0_dp)
! mz_cb_20100108-
! op_rd_20100108-
        ENDIF
     ENDDO
  ENDDO

  !-----
  ! Calculate Tice and Tnat if we are south of r_lat(1) 
  ! southern latitude, below r_mb(1) and above r_lb(1) hPa
  !-----
  DO jk=1,nlev
     DO jp=1,kproma
        IF( (philat_2d(jp,jrow) < r_lat(1))  .AND. &
             (press_3d(_RI_XYZ__(jp,jrow,jk)) >= r_mb(1))  .AND. &
             (press_3d(_RI_XYZ__(jp,jrow,jk)) <= r_lb(1)) ) THEN
           Tice(jp,jk) = T_ice(press_3d(_RI_XYZ__(jp,jrow,jk)) * &
! mz_cb_20100108+
!!$                H2O_tot(jp,jk)/1013.0_dp)
                H2O_tot(jp,jk)/100.0_dp)
! mz_cb_20100108-
           Tnat(jp,jk) = T_nat(press_3d(_RI_XYZ__(jp,jrow,jk)) * &
! op_cb_20100108+
!!$                H2O_tot(jp,jk)/1013.0_dp, &
                H2O_tot(jp,jk)/101325.0_dp, &
! op_cb_20100108-
! op_rd_20100108+
! op_cb_20100108+
                press_3d(_RI_XYZ__(jp,jrow,jk)) * tmp_HNO3_tot(jp,jk)/101325.0_dp)
! op_cb_20100108-
! op_rd_20100108-
        ENDIF
     ENDDO
  ENDDO

  !-----
  ! First guess for PSC region: we are in the the PSC region where Tice<Tnat. 
  !-----
  DO jk=1,nlev
     DO jp=1,kproma
        IF( Tice(jp,jk)<Tnat(jp,jk) ) THEN
           flt_pscreg(_RI_XYZ__(jp,jrow,jk))=1.0_dp
        ELSE
           flt_pscreg(_RI_XYZ__(jp,jrow,jk))=0.0_dp
        ENDIF
     ENDDO
  ENDDO

  !-----
  ! Avoid not-psc-boxes above psc-boxes
  !-----
  index_arr(:) = 0
  DO jk=1,nlev
     DO jp=1,kproma
        IF( (flt_pscreg(_RI_XYZ__(jp,jrow,jk)) < 0.001_dp) .AND. &
             (index_arr(jp) .EQ. 0) ) index_arr(jp) = jk
     ENDDO
  ENDDO
  DO jk=1,nlev
     DO jp=1,kproma
        IF( jk > index_arr(jp) ) flt_pscreg(_RI_XYZ__(jp,jrow,jk))=0.0_dp
     ENDDO
  ENDDO

  !-----
  ! Now remove highest levels from psc region
  !-----
  DO jk=1,nlev
     DO jp=1,kproma
! ka_ok_20100129+
!!$     ! mz_pj_20060224 1500->2000
!!$     IF( press_3d(jp,jk,jrow)<=2000.0_dp ) flt_pscreg(jp,jk,jrow) = 0.0_dp
        ! upper boundary above r_ub(1) and south of r_lat(1) degree 
        IF( (philat_2d(jp,jrow) < r_lat(1)) .AND. &
            (press_3d(_RI_XYZ__(jp,jrow,jk)) <= r_ub(1)) ) THEN
           flt_pscreg(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
        ENDIF
        ! upper boundary above r_ub(2) and north of r_lat(2) degree 
        IF( (philat_2d(jp,jrow) > r_lat(2)) .AND. &
            (press_3d(_RI_XYZ__(jp,jrow,jk)) <= r_ub(2)) ) THEN
           flt_pscreg(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
        ENDIF
! ka_ok_20100129-
     ENDDO
  ENDDO

  ! mz_pj_20070209+
  ! define flt_stratreg:
  DO jk=1,nlev
    DO jp=1,kproma
! op_pj_20180725+
      IF (.NOT. ledith) THEN
         zl = (flt_pscreg(_RI_XYZ__(jp,jrow,jk))>0.5) .OR. (jk<NINT(tp_i(jp,jrow)))
      ELSE
         !ka_sv_20170509+
         zl = ( (flt_pscreg(_RI_XYZ__(jp,jrow,jk))>0.5) .OR. (jk<NINT(tp_i(jp,jrow))) ) &
              .AND. (press_3d(_RI_XYZ__(jp,jrow,jk)) >= 1.0_dp)
         !ka_sv_20170509-
      ENDIF
      IF (zl) THEN
! op_pj_20180725-
        flt_stratreg(_RI_XYZ__(jp,jrow,jk)) = 1.0_dp
      ELSE
        flt_stratreg(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
      ENDIF
    ENDDO
  ENDDO
  ! mz_pj_20070209-

END SUBROUTINE msbm_local_start

!=============================================================================

SUBROUTINE msbm_physc
  !--------------------------------------------------------------------------
  !This subroutine calculates wether PSCs exist and
  ! which surface area and composition they have,
  ! as a necessary information to calculate heterogeneous reaction
  ! rates (msbm_chem).
  ! 
  !USE associated variables: 
  ! tm1_3d: temperature, [tm1_3d] = K
  ! press_3d: pressure (in the middle of grid boxes), [press_3d] = Pa
  ! pressi_3d: pressure (at level interfaces), [pressi_3d] = Pa
  !
!!$! q: mass of water vapour divided by mass of dry air, [q] = kg/kg
  ! q is actually in kg/kg moist air, hence should be divided by mass of 
  ! moist air ! op_ff_20170808
  !
  ! xl: mass of liquid water divided by mass of dry air, [xl] = kg/kg
  ! xi: mass of water ice divided by mass of dry air, [xi] = kg/kg
  !--------------------------------------------------------------------------
  ! ECHAM5/MESSy
  USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1, ntrac_gp
  USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma, nproma
  USE messy_main_data_bi,       ONLY: press_3d        &
                                    , pressi_3d       &
                                    , tm1_3d, tte_3d  &
                                    , qm1_3d => qm1   &
                                    , aclc            & ! op_pj_20171127
#ifdef ECHAM5
                                    , xlm1_3d => xlm1  &
                                    , xim1_3d => xim1  &
#endif
#if defined(COSMO) || defined(MESSYDWARF)
                                    , xlm1_3d  &
                                    , xim1_3d  &
#endif
#ifdef CESM1
                                    , xlm1_3d  &
                                    , xim1_3d  &
#endif
                                    , qte_3d   &
                                    , xlte_3d  &
                                    , xite_3d 
  USE messy_main_timer,         ONLY: time_step_len
  USE messy_main_constants_mem, ONLY: MHg, M_air, M_H2O & ! mz_hr_20120524
                                    , ccwmin ! op_pj_20171127
  USE messy_msbm, ONLY:                                                &
    mz_psc_phase,                                                     &
    mz_psc_ice_H2O, mz_psc_nat_HNO3,                                  &
    mz_psc_liq_bH2SO4b, mz_psc_liq_bHNO3b, mz_psc_liq_bH2SO4,         &
    mz_psc_liq_bHNO3, mz_psc_liq_partHNO3, mz_psc_liq_wnen,           &
    mz_psc_liq_wHNO3, mz_psc_liq_wH2SO4, mz_psc_liq_hHCl,             &
    mz_psc_liq_hHBr, mz_psc_liq_hHOCl, mz_psc_liq_hHOBr,              &
    mz_psc_liq_partHBr, mz_psc_liq_wHCl,                              &
    mz_psc_liq_wHOCl, mz_psc_liq_wHOBr,                               &
    mz_psc_liq_H2O, mz_psc_liq_HCl,                                   &
    mz_psc_liq_HOCl, mz_psc_liq_HOBr,                                 &
    mz_psc_N_solid, mz_psc_r_solid,                                   &
    mz_psc_dim_liquid, mz_psc_surface_liquid,                         &
    mz_psc_het_nat1, mz_psc_het_nat2, mz_psc_het_nat3,                &
    mz_psc_het_nat4, mz_psc_het_nat5, mz_psc_het_nat6,                &
    mz_psc_het_nat7, mz_psc_het_nat8, mz_psc_het_nat9,                &
    mz_psc_het_nat10, mz_psc_het_nat11,                               &
    mz_psc_het_ice12, mz_psc_het_ice13, mz_psc_het_ice14,             &
    mz_psc_het_ice15, mz_psc_het_ice16, mz_psc_het_ice17,             &
    mz_psc_het_ice18, mz_psc_het_ice19, mz_psc_het_ice20,             &
    mz_psc_het_ice21, mz_psc_het_ice22,                               &
    mz_psc_het_liq23, mz_psc_het_liq24, mz_psc_het_liq25,             &
    mz_psc_het_liq26, mz_psc_het_liq27, mz_psc_het_liq28,             &
    mz_psc_het_liq29, mz_psc_het_liq30,                               &
    mz_psc_het_liq_gen, mz_psc_het_sol_gen,                           &
    mz_psc_liq_check,                                                 &
    mz_psc_diff, mz_psc_density,                                      &
    mz_psc_vel, mz_psc_SedStep,                                       &
    mz_psc_sed,                                                       &
    ihs_N2O5_H2O,   ihs_HOCl_HCl,   ihs_ClNO3_HCl,   &
    ihs_ClNO3_H2O,  ihs_N2O5_HCl,   ihs_ClNO3_HBr,   &
    ihs_BrNO3_HCl,  ihs_HOCl_HBr,   ihs_HOBr_HCl,    &
    ihs_HOBr_HBr,   ihs_BrNO3_H2O,  ihs_Hg, ihs_RGM, &
    ! ka_ok_20100118+
    k_para_nat, mz_psc_N_ice, mz_psc_r_ice, mz_NAT_vel,               &
    diff_H2O, diff_HNO3, speed_HNO3, press_HNO3, press_H2O,           & 
    press_HNO3_over_NAT, press_H2OS, press_HNO3S, svapn_fkt, NSB, KinPar
    ! ka_ok_20100118-

  IMPLICIT none

  ! mz_hr_20120524+
  REAL(dp), PARAMETER :: mmratio  = M_air/M_H2O
  REAL(dp), PARAMETER :: rmmratio = M_H2O/M_air
  ! mz_hr_20120524-
  INTEGER :: jp, jk, i, tloop, tloop_index ! op_rd_20100108: tloop* added
  INTEGER :: jt                            ! ka_ok_20100118

  REAL(dp), DIMENSION(nproma,nlev) :: TEMP !temperature [K]
  REAL(dp), DIMENSION(nproma,nlev) :: PRESS !pressure [mbar=hPa]

  ! op_pj_20171127+
  REAL(DP), DIMENSION(nproma, nlev) :: xlp1, xip1
  LOGICAL :: lo, lo1
  ! op_pj_20171127-

  !-----
  ! local PSC relevant information:
  !-----
  ! "0" indicates information taken from ECHAM before PSC and Chemistry calc.
  ! "g","l","s","i","n" stand for gas, liquid, solid, ice and nat phase
  ! "t" stands for total
  ! "gl" stands for gas + liquid phase
  !-----
  ! partX: ratio of gas phase amount-of-substance of X 
  !        to gas-and-liquid-phase amount-of-substance of X / (mol/mol)
  ! bX: molality of X in water / (mol/kg)
  ! wX: mass fraction of X in liquid / (kg/kg)
  ! hX: Henry's law constant / (mol/(kg*atm))
  ! N_s: number density of solid aerosol particles in air / (1/m**3)
  ! r_s: radius of mono-modal distribution of solid particles / m
  ! N_l0: number density of liquid aerosol particles in air / (1/m**3)
  ! sigmaaero: width of liquid log-normal distribution / cm
  ! rmean_l: surface median radius of liquid log-normal distribution / cm
  ! A_l: surface area density of liquid particles / (cm**2/cm**3)
  ! ival: liquid parameterisation validity flag
  !-----
  REAL(dp), DIMENSION(nproma,nlev) :: H2SO4_t0
  REAL(dp), DIMENSION(nproma,nlev) :: H2O_t0, H2O_g0, H2O_l0, H2O_i0
  REAL(dp), DIMENSION(nproma,nlev) :: HNO3_t0, HNO3_gl0, HNO3_n0
  REAL(dp), DIMENSION(nproma,nlev) :: HBr_t0
  REAL(dp), DIMENSION(nproma,nlev) :: HOBr_t0, HOBr_g, HOBr_l
  REAL(dp), DIMENSION(nproma,nlev) :: BrNO3_t0
  REAL(dp), DIMENSION(nproma,nlev) :: HCl_t0
  REAL(dp), DIMENSION(nproma,nlev) :: HOCl_t0
  REAL(dp), DIMENSION(nproma,nlev) :: ClNO3_t0
  REAL(dp), DIMENSION(nproma,nlev) :: H2O_g, H2O_i, H2O_l, H2O_gl
  REAL(dp), DIMENSION(nproma,nlev) :: HNO3_g, HNO3_n, HNO3_l, HNO3_gl
  REAL(dp), DIMENSION(nproma,nlev) :: H2SO4_l
  REAL(dp), DIMENSION(nproma,nlev) :: HCl_g,  HCl_l
  REAL(dp), DIMENSION(nproma,nlev) :: HOCl_g, HOCl_l
  REAL(dp), DIMENSION(nproma,nlev) :: HBr_g, HBr_l
  REAL(dp), DIMENSION(nproma,nlev) :: partHNO3, partHBr
  REAL(dp), DIMENSION(nproma,nlev) :: bHNO3, bH2SO4, bHNO3b, bH2SO4b
  REAL(dp), DIMENSION(nproma,nlev) :: wHNO3, wH2SO4, wHCl, &
                                          wHOCl, wHOBr, wnen
  REAL(dp), DIMENSION(nproma,nlev) :: hHCl, hHBr, hHOCl, hHOBr
  REAL(dp), DIMENSION(nproma,nlev) :: N_s
  REAL(dp), PARAMETER :: sigmaaero=1.8_dp
  REAL(dp), DIMENSION(nproma,nlev) :: r_s, rmean_l, A_l, diff, dens
  INTEGER, DIMENSION(nproma,nlev) :: i_val, phase
  LOGICAL, DIMENSION(nproma,nlev) :: LCalcChemArray, val_psc, val_strat
  ! op_rd_20100108+
  REAL(dp), DIMENSION(nproma,nlev) :: r_s_offl
  INTEGER, DIMENSION(nproma,nlev)  :: phase_offl 
  ! op_rd_20100108-
  ! ka_ok_20100118+
  REAL(dp), DIMENSION(nproma,nlev,NSB) :: HNO3_nsarr_0
  REAL(dp), DIMENSION(nproma,nlev,NSB) :: HNO3_nsarr, N_NAT_arr, r_NAT_arr
  REAL(dp), DIMENSION(nproma,nlev)     :: DH2O, DHNO3, vhno3
  REAL(dp), DIMENSION(nproma,nlev)     :: pHNO3, pH2O, pHNO3_over_NAT, &
                                          pH2OS, pHNO3S, SVAPN
  ! radius NAT and ice and number density of NAT and ICE
  REAL(dp), DIMENSION(nproma,nlev)     :: r_N, r_i, N_N, N_i
  REAL(dp), DIMENSION(nproma,nlev,NSB) :: G_arr ! growth factor
  REAL(dp), DIMENSION(nproma,nlev,NSB) :: &   
       v_sed_NAT_arr,   &  ! sedimentation velocity of NAT particles 
                           ! in different sizebins
       SedStep_NAT_arr, &  ! sedimentation step of NAT particles
                           ! in different sizebins           
       HNO3ChangeDueToNATSed_arr
  ! ka_ok_20100118-

  !-----
  ! 2d variables relevant for PSCs
  !-----
  REAL(dp), DIMENSION(nproma, nlev) :: &
!!$    SedStep,                               &   ! ka_ok_20100118
    SedStep_Ice,                               &  ! ka_ok_20100118
    IceChangeDueToSed,                     &
    HNO3ChangeDueToSed

#ifndef _SX
    REAL(DP), PARAMETER :: REALZERO = 0.0_dp
#else
    REAL(DP), PARAMETER :: REALZERO = 1.0E-40_dp
#endif

    INTEGER :: idt, idt2

    ! op_pj_20170215+
    REAL(DP) :: zhumte(kproma,nlev,3)
    REAL(DP) :: zxtte(nproma,nlev,ntrac_gp)
#ifdef MECCA_TAG
    REAL(DP) :: pxtp1(nproma,nlev,ntrac_gp)
#endif
    ! op_pj_20170215-

  INTRINSIC REAL, MAX, INT, MIN, NINT, SUM
  !--------------------------------------------------------------------------

  !-----
  ! Initialise array form of MECCA switch
  !-----
  LCalcChemArray = LCalcChem
  
  !-----
  ! Initialise local flag variables, otherwise certain WHERE constructs
  ! would evaluate undefined values
  !-----
  i_val(:,:) = 0
  phase(:,:) = 0
  val_psc(:,:) = .false. 
  val_strat(:,:) = .false. 
  ! op_pj_20100108+
  IF (l_feedback) THEN
     tloop = 1
  ELSE
     tloop = 2
  END IF
  ! op_pj_20100108-
  ! op_pj_20170215+
  zhumte(:,:,:) = 0.0_dp
  zxtte(:,:,:) = 0.0_dp
#ifdef MECCA_TAG
  pxtp1(:,:,:) = 0.0_dp
  DO jt=1, ntrac_gp
     DO jp=1, kproma
        pxtp1(jp,:,jt) = pxtm1(_RI_X_ZN_(jp,:,jt)) + pxtte(_RI_X_ZN_(jp,:,jt))*time_step_len
     END DO
  END DO
#endif
  ! op_pj_20170215-
  HNO3_nsarr_0(:,:,:) = 0.0_dp
  
  !-----
  ! Local Latitude:
  !-----
  ! jrow = local latitude index (N -> S)

  ! mz_pj_20070209+
  !-----
  ! store information about stratosphere region locally in val_strat
  !-----
  DO jk=1,nlev
     DO jp=1,kproma
        IF (NINT(flt_stratreg(_RI_XYZ__(jp,jrow,jk)))==1) val_strat(jp,jk)=.TRUE.
     ENDDO
  ENDDO
  ! mz_pj_20070209-

  ! op_pj_20100729+
  IF (LCalcChem) THEN
     DO i=1, NKHET
        khet(i)%ptr(_RI_XYZ__(1:kproma,jrow,:)) = 0.0_dp
     END DO
  END IF
  ! op_pj_20100729-

  !-----
  ! store information about psc relevant region locally in val_psc
  !-----
  DO jk=1,nlev
     DO jp=1,kproma
        IF (NINT(flt_pscreg(_RI_XYZ__(jp,jrow,jk)))==1) val_psc(jp,jk)=.TRUE.
     ENDDO
  ENDDO

  !-----
  ! initialise temperature according to previous timestep and other tendencies
  ! and add artificial temperature shift from namelist 
  !-----
  TEMP(1:kproma,:) = tm1_3d(_RI_XYZ__(1:kproma,jrow,:)) &
                    +tte_3d(_RI_XYZ__(1:kproma,jrow,:)) * time_step_len &
                    +TempShift

  level_loop: DO jk=1,nlev
     vector_loop: DO jp=1,kproma
       ! mz_pj_20070209+
       !IF (val_psc(jp,jk)) THEN
       if_val_strat: IF (val_strat(jp,jk)) THEN
       ! mz_pj_20070209-
        ! -------------------------------------------------------
        ! start of block that calculates aerosol properties.
        ! -------------------------------------------------------
        !------------------------------------------------------------------
        ! The following calculations take place inside of the psc relevant
        ! region.
        !------------------------------------------------------------------

        !-----
        ! Initialise pressure acc. to previous timestep and other tendencies:
        !-----
        ! PRESS = pressure / hPa
        !-----
        PRESS(jp,jk)= press_3d(_RI_XYZ__(jp,jrow,jk))/100.0_dp

! op_rd_20100108+
! moved to loop below for consistency (equivalence) with HNO3-climatology
!!$        !-----
!!$        ! H2SO4 initialisation
!!$        !-----
!!$        H2SO4_t0(jp,jk) = H2SO4(jp,jk,jrow)
! op_rd_20100108-

        !-----
        ! H2O initialisation
        !-----
        H2O_i0(jp,jk)  = xim1_3d(_RI_XYZ__(jp,jrow,jk)) &
                        +xite_3d(_RI_XYZ__(jp,jrow,jk))*time_step_len
        H2O_l0(jp,jk)  = xlm1_3d(_RI_XYZ__(jp,jrow,jk)) &
                        +xlte_3d(_RI_XYZ__(jp,jrow,jk))*time_step_len
        H2O_g0(jp,jk)  = qm1_3d(_RI_XYZ__(jp,jrow,jk)) &
                        +qte_3d(_RI_XYZ__(jp,jrow,jk))*time_step_len

        !-----
        ! [kg H2O /kg air] --> [mol H2O/mol air]
        !-----
        H2O_i0(jp,jk)  = H2O_i0(jp,jk) * mmratio
        H2O_l0(jp,jk)  = H2O_l0(jp,jk) * mmratio
! op_ff_20170904+
        !!! ice and liquid water is in kg/kg dry air, 
        !!! but qm1 is in kg/kg_moistair !!!
!!$     H2O_g0(jp,jk)  = H2O_g0(jp,jk) * mmratio
        H2O_g0(jp,jk)  = H2O_g0(jp,jk) * mmratio / (1.0_dp - H2O_g0(jp,jk) )
! op_ff_20170904-

        H2O_t0(jp,jk)  = H2O_g0(jp,jk) + H2O_i0(jp,jk) + H2O_l0(jp,jk)

        tloop_loop: DO tloop_index = 1,tloop             ! op_rd_20100108

           IF (l_feedback .OR. (tloop_index == 2)) THEN  ! op_rd_20100108
              H2SO4_t0(jp,jk) = H2SO4(_RI_XYZ__(jp,jrow,jk))         ! op_pj_20100108
              idt = idt_HNO3_nat
              HNO3_n0(jp,jk) = pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                   +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len
              ! ka_ok_20100118+
              IF (KinPar) THEN
                 DO jt=1, NSB
                    idt=idt_HNO3_natsize(jt)
                    HNO3_nsarr_0(jp,jk,jt) = pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                         +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len 
                 END DO
              END IF
              ! ka_ok_20100118-
              idt = idt_HNO3
              HNO3_gl0(jp,jk)= pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                   +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len
           ELSE                  ! op_rd_20100108
              ! op_rd_20100108+
              H2SO4_t0(jp,jk) = H2SO4(_RI_XYZ__(jp,jrow,jk))
              HNO3_n0(jp,jk)  = 0.0_dp
              HNO3_gl0(jp,jk) = predef_HNO3_tot(_RI_XYZ__(jp,jrow,jk))
           ENDIF
           ! op_rd_20100108-

           HNO3_t0(jp,jk) = HNO3_gl0(jp,jk) + HNO3_n0(jp,jk)

           idt = idt_HBr
           HBr_t0(jp,jk)  = pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len
           
           idt = idt_HOBr
           HOBr_t0(jp,jk) = pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len

           idt = idt_BrNO3
           BrNO3_t0(jp,jk)= pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len

           idt = idt_HCl
           HCl_t0(jp,jk)  = pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len

           idt = idt_HOCl
           HOCl_t0(jp,jk) = pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len

           idt = idt_ClNO3
           ClNO3_t0(jp,jk)= pxtm1(_RI_X_ZN_(jp,jk,idt)) &
                +pxtte(_RI_X_ZN_(jp,jk,idt))*time_step_len

           !-----
           ! set negative values to zero
           !-----
           H2SO4_t0(jp,jk) = MAX(H2SO4_t0(jp,jk),REALZERO)
           H2O_i0(jp,jk)   = MAX(H2O_i0(jp,jk),REALZERO)
           H2O_l0(jp,jk)   = MAX(H2O_l0(jp,jk),REALZERO)
           H2O_g0(jp,jk)   = MAX(H2O_g0(jp,jk),REALZERO)
           H2O_t0(jp,jk)   = MAX(H2O_t0(jp,jk),1.0e-10_dp)
           HNO3_n0(jp,jk)  = MAX(HNO3_n0(jp,jk),REALZERO)
           HNO3_gl0(jp,jk) = MAX(HNO3_gl0(jp,jk),REALZERO)
           HNO3_t0(jp,jk)  = MAX(HNO3_t0(jp,jk),REALZERO)
           HBr_t0(jp,jk)   = MAX(HBr_t0(jp,jk),REALZERO)
           HOBr_t0(jp,jk)  = MAX(HOBr_t0(jp,jk),REALZERO)
           BrNO3_t0(jp,jk) = MAX(BrNO3_t0(jp,jk),REALZERO)
           HCl_t0(jp,jk)   = MAX(HCl_t0(jp,jk),REALZERO)
           HOCl_t0(jp,jk)  = MAX(HOCl_t0(jp,jk),REALZERO)
           ClNO3_t0(jp,jk) = MAX(ClNO3_t0(jp,jk),REALZERO)

           !-----
           ! make actual phase indices (from psc-channel) usable:
           !-----
           phase(jp,jk)=INT(flt_phase(_RI_XYZ__(jp,jrow,jk)))

           !----------------------------
           !Treatment of solid particles
           !----------------------------
           !-----
           ! determine for thermodynamically stable phases
           !-----
           phase(jp,jk) = mz_psc_phase(TEMP(jp,jk),          &
                PRESS(jp,jk), H2O_t0(jp,jk),  &
                HNO3_t0(jp,jk), phase(jp,jk), &
                H2O_i0(jp,jk), HNO3_n0(jp,jk))

           !-----
           ! calculate partitioning between gas phase, ice and NAT
           !-----
           H2O_i(jp,jk)   = mz_psc_ice_H2O(phase(jp,jk),    &
                TEMP(jp,jk), PRESS(jp,jk), &
                H2O_t0(jp,jk))
           H2O_gl(jp,jk)  = H2O_t0(jp,jk)-H2O_i(jp,jk)
           H2O_l(jp,jk)   = 0.0_dp
           H2O_g(jp,jk)   = H2O_gl(jp,jk)
           IF (.NOT.KinPar) THEN   ! ka_ok_20100118
              HNO3_n(jp,jk)  = mz_psc_nat_HNO3(phase(jp,jk),      &
                   TEMP(jp,jk), PRESS(jp,jk),    &
                   H2O_t0(jp,jk), HNO3_t0(jp,jk))
           ! ka_ok_20100118+
           ELSE
              DH2O(jp,jk)=diff_H2O(TEMP(jp,jk),PRESS(jp,jk))
              DHNO3(jp,jk)=diff_HNO3(DH2O(jp,jk))
              vhno3(jp,jk)=speed_HNO3(TEMP(jp,jk))
              pHNO3(jp,jk)=press_HNO3(HNO3_gl0(jp,jk),TEMP(jp,jk),PRESS(jp,jk))
              pH2O(jp,jk)=press_H2O(H2O_g(jp,jk),TEMP(jp,jk),PRESS(jp,jk))
              pHNO3_over_NAT(jp,jk)=press_HNO3_over_NAT(TEMP(jp,jk),pH2O(jp,jk))
              pH2OS(jp,jk)=press_H2OS(H2O_g(jp,jk),TEMP(jp,jk),PRESS(jp,jk))
              pHNO3S(jp,jk)=press_HNO3S(pH2OS(jp,jk),TEMP(jp,jk))
              SVAPN(jp,jk)=svapn_fkt(HNO3_gl0(jp,jk),pH2OS(jp,jk), &
                   PRESS(jp,jk),TEMP(jp,jk))
              CALL K_PARA_NAT(HNO3_nsarr_0(jp,jk,:), HNO3_nsarr(jp,jk,:), &
                   DHNO3(jp,jk), &
                   vhno3(jp,jk), pHNO3(jp,jk), pHNO3_over_NAT(jp,jk), &
                   pHNO3S(jp,jk), &
                   svapn(jp,jk), TEMP(jp,jk), PRESS (jp,jk), time_step_len, &
                   r_N(jp,jk), N_NAT_arr(jp,jk,:), r_NAT_arr(jp,jk,:), &
                   G_arr(jp,jk,:))
              ! convert from 1/cm^3 to 1/m^3
              N_NAT_arr(jp,jk,:)=N_NAT_arr(jp,jk,:)*1.0e6_dp 
              N_N(jp,jk)=sum(N_NAT_arr(jp,jk,:),DIM=1)
              HNO3_n(jp,jk)=SUM(HNO3_nsarr(jp,jk,:),DIM=1)
           ENDIF
           ! ka_ok_20100118-
           ! op_rd_20100108+
           !---------------------------- 
           ! in case of offline mode incorporate changes in 
           ! HNO3_n into clim_HNO3_nat
           !----------------------------
           IF (.NOT. l_feedback .AND. (tloop_index == 1)) THEN 
              clim_HNO3_nat(_RI_XYZ__(jp,jrow,jk)) = HNO3_n(jp,jk)
           END IF
           ! op_rd_20100108-

           HNO3_gl(jp,jk) = HNO3_t0(jp,jk) - HNO3_n(jp,jk)
           HNO3_l(jp,jk)  = 0.0_dp
           HNO3_g(jp,jk)  = HNO3_gl(jp,jk)

           IF (phase(jp,jk)==3 .OR. phase(jp,jk)==2) THEN
              !-----
              ! calculate size of ice and NAT particles
              !-----
              N_s(jp,jk) = mz_psc_N_solid(phase(jp,jk), &
                   TEMP(jp,jk), PRESS(jp,jk), &
                   H2O_i(jp,jk), HNO3_n(jp,jk))
              r_s(jp,jk) = mz_psc_r_solid(phase(jp,jk),    &
                   TEMP(jp,jk), PRESS(jp,jk),    &
                   N_s(jp,jk), H2O_i(jp,jk), &
                   HNO3_n(jp,jk))
           ELSE   ! IF (phase(jp,jk)==3 .OR. phase(jp,jk)==2) THEN ...
              N_s(jp,jk) = 0.0_dp
              r_s(jp,jk) = 0.0_dp
           END IF   ! IF (phase(jp,jk)==3 .OR. phase(jp,jk)==2) THEN ...

           ! ka_ok_20100118+
           IF (KinPar) THEN
              IF (phase(jp,jk)==3) THEN
                 N_i(jp,jk) = mz_psc_N_ice(phase(jp,jk), &
                              TEMP(jp,jk), PRESS(jp,jk), &
                              H2O_i(jp,jk))
                 r_i(jp,jk) = mz_psc_r_ice(phase(jp,jk),    &
                              TEMP(jp,jk), PRESS(jp,jk),    &
                              N_i(jp,jk), H2O_i(jp,jk))
              ELSE
                 N_i(jp,jk) = 0.0_dp
                 r_i(jp,jk) = 0.0_dp
              END IF
           ELSE
              N_i(jp,jk) = N_s(jp,jk)
              r_i(jp,jk) = r_s(jp,jk)
              N_N(jp,jk) = N_s(jp,jk)
              r_N(jp,jk) = r_s(jp,jk)
           END IF
           ! ka_ok_20100118-

           !------------------------------
           ! Treatment of liquid particles
           !------------------------------
           !-----
           ! check conditions whether STS parameterisation holds
           !-----
           i_val(jp,jk)=mz_psc_liq_check(TEMP(jp,jk),       &
                PRESS(jp,jk), H2O_g(jp,jk),    &
                HNO3_g(jp,jk), H2SO4_t0(jp,jk))

           !-----
           ! calculate partitioning (fraction of HX in the gas phase),
           ! henrys law coefficients and aerosolcomposition
           !-----
           bH2SO4b(jp,jk) = mz_psc_liq_bH2SO4b(TEMP(jp,jk),       &
                PRESS(jp,jk), H2O_gl(jp,jk))
           bHNO3b(jp,jk)  = mz_psc_liq_bHNO3b(TEMP(jp,jk),        &
                PRESS(jp,jk), H2O_gl(jp,jk))
           bH2SO4(jp,jk)  = mz_psc_liq_bH2SO4(TEMP(jp,jk),        &
                PRESS(jp,jk), H2SO4_t0(jp,jk),   &
                H2O_gl(jp,jk), HNO3_gl(jp,jk),   &
                bH2SO4b(jp,jk), bHNO3b(jp,jk))
           bHNO3(jp,jk)   = mz_psc_liq_bHNO3(TEMP(jp,jk),         &
                PRESS(jp,jk), H2O_gl(jp,jk),     &
                HNO3_gl(jp,jk), bHNO3b(jp,jk),   &
                bH2SO4b(jp,jk), bH2SO4(jp,jk))
           partHNO3(jp,jk)= mz_psc_liq_partHNO3(TEMP(jp,jk),      &
                PRESS(jp,jk), H2O_gl(jp,jk),     &
                HNO3_gl(jp,jk), bHNO3(jp,jk),    &
                bH2SO4(jp,jk), bH2SO4b(jp,jk))
           wnen(jp,jk)    = mz_psc_liq_wnen(TEMP(jp,jk),          &
                bHNO3(jp,jk), bH2SO4(jp,jk),     &
                bH2SO4b(jp,jk), partHNO3(jp,jk))
           wHNO3(jp,jk)   = mz_psc_liq_wHNO3(TEMP(jp,jk),         &
                bHNO3(jp,jk), bH2SO4(jp,jk),     &
                bH2SO4b(jp,jk), partHNO3(jp,jk), &
                wnen(jp,jk))
           wH2SO4(jp,jk)  = mz_psc_liq_wH2SO4(TEMP(jp,jk),        &
                bH2SO4(jp,jk), bH2SO4b(jp,jk),   &
                partHNO3(jp,jk), wnen(jp,jk))
           hHCl(jp,jk)    = mz_psc_liq_hHCl(TEMP(jp,jk),          &
                wHNO3(jp,jk), wH2SO4(jp,jk),     &
                wnen(jp,jk))
           hHBr(jp,jk)    = mz_psc_liq_hHBr(TEMP(jp,jk),          &
                wHNO3(jp,jk), wH2SO4(jp,jk),     &
                wnen(jp,jk))
           hHOCl(jp,jk)   = mz_psc_liq_hHOCl(TEMP(jp,jk),         &
                bH2SO4(jp,jk), bHNO3(jp,jk))
           hHOBr(jp,jk)   = mz_psc_liq_hHOBr(hHOCl(jp,jk))
           partHBr(jp,jk) = mz_psc_liq_partHBr(PRESS(jp,jk),      &
                HBr_t0(jp,jk), H2SO4_t0(jp,jk),  &
                hHBr(jp,jk), bH2SO4(jp,jk))
           wHCl(jp,jk)    = mz_psc_liq_wHCl(PRESS(jp,jk),         &
                HCl_t0(jp,jk), H2SO4_t0(jp,jk),  &
                hHCl(jp,jk), bH2SO4(jp,jk),      &
                wnen(jp,jk))
           wHOCl(jp,jk)   = mz_psc_liq_wHOCl(PRESS(jp,jk),        &
                HOCl_t0(jp,jk), H2SO4_t0(jp,jk), &
                hHOCl(jp,jk), bH2SO4(jp,jk),     &
                wnen(jp,jk))
           wHOBr(jp,jk)   = mz_psc_liq_wHOBr(PRESS(jp,jk),        &
                HOBr_t0(jp,jk), H2SO4_t0(jp,jk), &
                hHOBr(jp,jk), bH2SO4(jp,jk),     &
                wnen(jp,jk))

           !-----
           ! calculate new gas and liquid phase amount-of-substance ratios
           ! - liquid phase amount-of-substance ratios follow from 
           !   parameterisations
           ! - gas phase amount-of-substance ratios are then calculated as
           !   "gas = gas+liquid - liquid"
           !-----
           H2SO4_l(jp,jk) = H2SO4_t0(jp,jk)
           H2O_l(jp,jk)   = mz_psc_liq_H2O(wH2SO4(jp,jk),         &
                wHNO3(jp,jk), H2SO4_l(jp,jk),    &
                H2O_t0(jp,jk))
           H2O_g(jp,jk)   = H2O_gl(jp,jk)-H2O_l(jp,jk)
           HNO3_l(jp,jk)  = (1.0_dp - partHNO3(jp,jk))            &
                *HNO3_gl(jp,jk)
           HNO3_g(jp,jk)  = HNO3_gl(jp,jk) - HNO3_l(jp,jk)
           HCl_l(jp,jk)   = mz_psc_liq_HCl(wHCl(jp,jk),           &
                H2SO4_l(jp,jk), wH2SO4(jp,jk),   &
                HCl_t0(jp,jk))
           HCl_g(jp,jk)   = HCl_t0(jp,jk) - HCl_l(jp,jk)
           HOCl_l(jp,jk)  = mz_psc_liq_HOCl(wHOCl(jp,jk),         &
                H2SO4_l(jp,jk), wH2SO4(jp,jk),   &
                HOCl_t0(jp,jk))
           HOCl_g(jp,jk)  = HOCl_t0(jp,jk) - HOCl_l(jp,jk)
           HBr_l(jp,jk)   = (1.0_dp - partHBr(jp,jk))             &
                *HBr_t0(jp,jk)
           HBr_g(jp,jk)   = HBr_t0(jp,jk) - HBr_l(jp,jk)
           HOBr_l(jp,jk)  = mz_psc_liq_HOBr(wHOBr(jp,jk),         &
                H2SO4_l(jp,jk), wH2SO4(jp,jk),   &
                HOBr_t0(jp,jk))
           HOBr_g(jp,jk)  = HOBr_t0(jp,jk) - HOBr_l(jp,jk)

           !-----
           ! calculate composition dependent size distribution
           !-----
           A_l(jp,jk)     = mz_psc_surface_liquid(TEMP(jp,jk),  &
                PRESS(jp,jk), H2O_l(jp,jk),         &
                HNO3_l(jp,jk), H2SO4_l(jp,jk))
           rmean_l(jp,jk) = mz_psc_dim_liquid(A_l(jp,jk),       &
                sigmaaero)


           ! op_rd_20100108+
           IF (tloop_index == 1) THEN   ! op_pj_20110304
              !-----
              ! Tendencies in H2O_g, H2O_l, H2O_i:
              ! qte_3d, xlte_3d, xite_3d= .... [kg H2O/kg air]
              !-----
              ! In offline mode update water tendencies only if tloop_index==1 
              ! (climatological impact of HNO3, HNO3_nat on water). 
              ! Water tendencies are
              ! neglected if tloop_index==2 (online HNO3 and HNO3_nat).
              !-----
              zhumte(jp,jk,1) = ((H2O_l(jp,jk)-H2O_l0(jp,jk)) * rmmratio)    &
                   /time_step_len
              zhumte(jp,jk,2) = ((H2O_i(jp,jk)-H2O_i0(jp,jk)) * rmmratio)    &
                   /time_step_len
              zhumte(jp,jk,3) = ((H2O_g(jp,jk)/(mmratio + H2O_g(jp,jk))) - &
                   (H2O_g0(jp,jk)/(mmratio + H2O_g0(jp,jk))) )             &
                   /time_step_len
              !---
              ! r_s_offl needed further down in psc module 
              ! phase_offl needed in psc module during next time step 
              !---
              r_s_offl(jp,jk)   = r_s(jp,jk)
              phase_offl(jp,jk) = phase(jp,jk)

           END IF

        END DO tloop_loop ! DO tloop_index = 1,tloop
        ! op_rd_20100108-

        ! -------------------------------------------------------
        ! end of block that calculates aerosol properties.
        ! -------------------------------------------------------

        ! -------------------------------------------------------
        ! start calculating khet 23...30
        ! -------------------------------------------------------

        IF (LCalcChem) THEN
          !-----
          ! calculate heterogeneous reaction rates
          !-----
          dens(jp,jk)  = mz_psc_density(wH2SO4(jp,jk),         &
                         wHNO3(jp,jk), TEMP(jp,jk))
          diff(jp,jk)  = mz_psc_diff(TEMP(jp,jk),              &
                         dens(jp,jk), bHNO3(jp,jk),            &
                         bH2SO4(jp,jk))

          psc_only_01: IF (val_psc(jp, jk)) THEN ! mz_pj_20070209
          khet(23)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq23(TEMP(jp,jk),    &
                              PRESS(jp,jk), HCl_t0(jp,jk),     &
                              HOCl_t0(jp,jk), bH2SO4(jp,jk),   &
                              bHNO3(jp,jk), wHCl(jp,jk),       &
                              dens(jp,jk), hHOCl(jp,jk),       &
                              A_l(jp,jk), rmean_l(jp,jk),      &
                              diff(jp,jk))
          khet(24)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq24(TEMP(jp,jk),    &
                              PRESS(jp,jk), H2O_g(jp,jk),      &
                              HCl_t0(jp,jk), ClNO3_t0(jp,jk),  &
                              wHCl(jp,jk), dens(jp,jk),        &
                              A_l(jp,jk), rmean_l(jp,jk))
          khet(25)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq25(TEMP(jp,jk),    &
                              PRESS(jp,jk), H2O_g(jp,jk),      &
                              wHCl(jp,jk), dens(jp,jk),        &
                              A_l(jp,jk), rmean_l(jp,jk))
          END IF psc_only_01 ! mz_pj_20070209

          khet(26)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq26(TEMP(jp,jk),    &
                              PRESS(jp,jk), H2O_g(jp,jk),      &
                              A_l(jp,jk))

          psc_only_02: IF (val_psc(jp, jk)) THEN ! mz_pj_20070209
          khet(27)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq27(TEMP(jp,jk),    &
                              PRESS(jp,jk), HCl_t0(jp,jk),     &
                              HOBr_t0(jp,jk), bH2SO4(jp,jk),   &
                              bHNO3(jp,jk), wHCl(jp,jk),       &
                              dens(jp,jk), hHOBr(jp,jk),       &
                              A_l(jp,jk), rmean_l(jp,jk),      &
                              diff(jp,jk))
          khet(28)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq28(TEMP(jp,jk),    &
                              PRESS(jp,jk), HBr_t0(jp,jk),     &
                              HOBr_t0(jp,jk), partHBr(jp,jk),  &
                              bH2SO4(jp,jk), bHNO3(jp,jk),     &
                              wHOBr(jp,jk), dens(jp,jk),       &
                              hHBr(jp,jk), hHOBr(jp,jk),       &
                              A_l(jp,jk), rmean_l(jp,jk),      &
                              diff(jp,jk))
          khet(29)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq29(TEMP(jp,jk),    &
                              PRESS(jp,jk), HBr_t0(jp,jk),     &
                              HOCl_t0(jp,jk), partHBr(jp,jk),  &
                              bH2SO4(jp,jk), bHNO3(jp,jk),     &
                              wHOCl(jp,jk), dens(jp,jk),       &
                              hHOCl(jp,jk), hHBr(jp,jk),       &
                              A_l(jp,jk), rmean_l(jp,jk),      &
                              diff(jp,jk))
          END IF psc_only_02 ! mz_pj_20070209

          khet(30)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq30(TEMP(jp,jk),    &
                              PRESS(jp,jk), H2O_g(jp,jk),      &
                              A_l(jp,jk))

          ! gamma = 0.001
          khet(33)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq_gen(TEMP(jp,jk),    &
                              PRESS(jp,jk), H2O_g(jp,jk), MHg, 0.001_dp, &
                              A_l(jp,jk))
          ! gamma = 0.1
          khet(36)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_liq_gen(TEMP(jp,jk),    &
                              PRESS(jp,jk), H2O_g(jp,jk), MHg, 0.1_dp, &
                              A_l(jp,jk))
        END IF   ! IF (LCalcChem) THEN ...

        ! -------------------------------------------------------
        ! end calculating khet 23...30
        ! -------------------------------------------------------
        psc_only_03: IF (val_psc(jp, jk)) THEN ! mz_pj_20070209
        ! -------------------------------------------------------
        ! start calculating khet 1...11
        ! -------------------------------------------------------

! ka_ok_20100118+
!!$        IF (LCalcChem .AND. phase(jp,jk)==2) THEN
        IF ( ((.NOT.KinPar) .AND. (LCalcChem .AND. phase(jp,jk)==2)) .OR. &
             ( KinPar .AND. LCalcChem  ) ) THEN
! ka_ok_20100118-
          !-----
          ! calculate het. reaction rates on nat particles
          !-----
! r_s -> r_N, N_s -> N_N ! ka_ok_20100118
          khet(1)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat1(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), HCl_g(jp,jk))
          khet(2)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat2(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), H2O_g(jp,jk))
          khet(3)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat3(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), HCl_g(jp,jk))
          khet(4)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat4(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), HCl_g(jp,jk))
          khet(5)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat5(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), H2O_g(jp,jk))
          khet(6)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat6(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), HBr_g(jp,jk))
          khet(7)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat7(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), HCl_g(jp,jk))
          khet(8)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat8(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), HBr_g(jp,jk))
          khet(9)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat9(TEMP(jp,jk),        &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), HCl_g(jp,jk))
          khet(10)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat10(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), HBr_g(jp,jk))
          khet(11)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_nat11(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), H2O_g(jp,jk))
          ! gamma = 0.001
          khet(31)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_sol_gen(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), H2O_g(jp,jk), MHg, 0.001_dp)
          ! gamma = 0.1
          khet(34)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_sol_gen(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_N(jp,jk),        &
                                N_N(jp,jk), H2O_g(jp,jk), MHg, 0.1_dp)
! ka_ok_20100118+
!!$        ELSEIF (LCalcChem .AND. phase(jp,jk)/=2) THEN
        ELSEIF (.NOT.KinPar .AND. (LCalcChem .AND. phase(jp,jk)/=2)) THEN
! ka_ok_20100118-
          !-----
          ! Set reaction rates for reactions on nat to zero if there is no nat,
          ! otherwise values from previous time steps would persist.
          !-----
           khet( 1)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet( 2)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet( 3)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet( 4)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet( 5)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet( 6)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet( 7)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet( 8)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet( 9)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet(10)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet(11)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet(31)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           khet(34)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
!           DO i=1, 11
!              khet(i)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
!           END DO
        END IF   ! IF (LCalcChem .AND. phase(jp,jk)==2) THEN ...
        
        ! -------------------------------------------------------
        ! end calculating khet 1...11
        ! -------------------------------------------------------

        ! -------------------------------------------------------
        ! start calculating khet 12...22
        ! -------------------------------------------------------

        IF (LCalcChem .AND. phase(jp,jk)==3) THEN
          !-----
          ! assume, that even though solid particles are an ice-NAT mixture
          ! the het. react. are similar to those on a pure ice surface
          !-----
! r_s -> r_i, N_s -> N_i ka_ok_20100118
          khet(12)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice12(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), HCl_g(jp,jk))
          khet(13)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice13(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), H2O_g(jp,jk))
          khet(14)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice14(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), HCl_g(jp,jk))
          khet(15)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice15(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), HCl_g(jp,jk))
          khet(16)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice16(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), H2O_g(jp,jk))
          khet(17)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice17(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), HBr_g(jp,jk))
          khet(18)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice18(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), HCl_g(jp,jk))
          khet(19)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice19(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), HBr_g(jp,jk))
          khet(20)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice20(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), HCl_g(jp,jk))
          khet(21)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice21(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), HBr_g(jp,jk))
          khet(22)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_ice22(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), H2O_g(jp,jk))
          ! gamma = 0.001
          khet(32)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_sol_gen(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), H2O_g(jp,jk), MHg, 0.001_dp)
          ! gamma = 0.1
          khet(35)%ptr(_RI_XYZ__(jp,jrow,jk))= mz_psc_het_sol_gen(TEMP(jp,jk),      &
                                PRESS(jp,jk), r_i(jp,jk),        &
                                N_i(jp,jk), H2O_g(jp,jk), MHg, 0.1_dp)
        ELSEIF (LCalcChem .AND. phase(jp,jk)/=3) THEN
          !-----
          ! Set reaction rates for heterogeneous chemical reactions on ice to
          ! zero if there is no ice, otherwise values from previous time
          ! steps would persist.
          ! Note that nullifying reaction rates is unnecessary if
          ! LCalcChem==.false., however, it does not cause problems either.
          !-----
           khet(12)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(13)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(14)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(15)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(16)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(17)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(18)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(19)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(20)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(21)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(22)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(32)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           khet(35)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
!           DO i=12, 22
!              khet(i)%ptr(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
!           END DO
        END IF   ! IF (LCalcChem .AND. phase(jp,jk)==3) THEN ...

        ! -------------------------------------------------------
        ! end calculating khet 12...22
        ! -------------------------------------------------------
        END IF psc_only_03 ! mz_pj_20070209

        !------------------------------------------------------------------
        ! Calculate new tendencies for Tracers
        !------------------------------------------------------------------
        ! Tendencies due to heterogeneous reactions are calculated in mecca
        ! Note that mecca must only "see" gas phase amount-of-substance
        ! fractions
        !-----
        zxtte(jp,jk,idt_HNO3) = zxtte(jp,jk,idt_HNO3) + &
             (HNO3_g(jp,jk)-HNO3_gl0(jp,jk)) / time_step_len
        zxtte(jp,jk,idt_HNO3_nat) = zxtte(jp,jk,idt_HNO3_nat) + &
             (HNO3_n(jp,jk)-HNO3_n0(jp,jk))/time_step_len
        zxtte(jp,jk,idt_HCl) = zxtte(jp,jk,idt_HCl) + &
             (HCl_g(jp,jk)-HCl_t0(jp,jk)) / time_step_len
        zxtte(jp,jk,idt_HOCl) = zxtte(jp,jk,idt_HOCl) + &
             (HOCl_g(jp,jk)-HOCl_t0(jp,jk))/ time_step_len
        zxtte(jp,jk,idt_HBr) = zxtte(jp,jk,idt_HBr) + &
             (HBr_g(jp,jk)-HBr_t0(jp,jk))/time_step_len
        zxtte(jp,jk,idt_HOBr) = zxtte(jp,jk,idt_HOBr) + &
             (HOBr_g(jp,jk)-HOBr_t0(jp,jk))/time_step_len
        IF (KinPar) THEN
           DO jt=1, NSB
              idt = idt_HNO3_natsize(jt)
              zxtte(jp,jk,idt) = zxtte(jp,jk,idt) + &
                   (HNO3_nsarr(jp,jk,jt)-HNO3_nsarr_0(jp,jk,jt)) &
                   /time_step_len
           END DO
        END IF

        !-----------------------------------------
        ! Write results (model output) to channels
        !-----------------------------------------
        HNO3_liq(_RI_XYZ__(jp,jrow,jk)) = HNO3_l(jp,jk)
        HNO3_nat(_RI_XYZ__(jp,jrow,jk)) = HNO3_n(jp,jk)
        HNO3_gas(_RI_XYZ__(jp,jrow,jk)) = HNO3_g(jp,jk)
        HCl_liq(_RI_XYZ__(jp,jrow,jk))  = HCl_l(jp,jk)
        HCl_gas(_RI_XYZ__(jp,jrow,jk))  = HCl_g(jp,jk)
        HOCl_liq(_RI_XYZ__(jp,jrow,jk)) = HOCl_l(jp,jk)
        HOCl_gas(_RI_XYZ__(jp,jrow,jk)) = HOCl_g(jp,jk)
        HBr_liq(_RI_XYZ__(jp,jrow,jk))  = HBr_l(jp,jk)
        HBr_gas(_RI_XYZ__(jp,jrow,jk))  = HBr_g(jp,jk)
        HOBr_liq(_RI_XYZ__(jp,jrow,jk)) = HOBr_l(jp,jk)
        HOBr_gas(_RI_XYZ__(jp,jrow,jk)) = HOBr_g(jp,jk)
        ! ka_ok_20100118+
        IF (KinPar) THEN
           DO jt=1, NSB
              HNO3_natsize(jt)%ptr(_RI_XYZ__(jp,jrow,jk)) = HNO3_nsarr(jp,jk,jt)
           END DO
        END IF
        ! ka_ok_20100118-
! op_rd_20100108+
! Note: in case of l_feedback, phase_offl contains the "online"-value,
!       see loop above ...
        flt_phase(_RI_XYZ__(jp,jrow,jk)) = REAL(phase_offl(jp,jk),KIND=dp)
! op_rd_20100108-
        N_solid(_RI_XYZ__(jp,jrow,jk))   = N_s(jp,jk)
        flt_val(_RI_XYZ__(jp,jrow,jk))   = REAL(i_val(jp,jk),KIND=dp)
        A_liq(_RI_XYZ__(jp,jrow,jk))     = A_l(jp,jk)
        r_SurfMed(_RI_XYZ__(jp,jrow,jk)) = rmean_l(jp,jk)
        ! ka_ok_20100118+
        r_ice(_RI_XYZ__(jp,jrow,jk))     = r_i(jp,jk)
        IF (KinPar) THEN
           N_ice(_RI_XYZ__(jp,jrow,jk))     = N_i(jp,jk)
           N_NAT(_RI_XYZ__(jp,jrow,jk))     = N_N(jp,jk)
           r_NAT(_RI_XYZ__(jp,jrow,jk))     = r_N(jp,jk)
        END IF
        ! ka_ok_20100118-

        ! mz_pj_20070209+
      !ELSE   ! IF (val_psc(jp,jk)) THEN ...
      ELSE   ! IF (val_strat(jp,jk)) THEN ...
        ! mz_pj_20070209-
        !----------------------------------------------------------------------
        ! mz_pj_20070209+
        !! The following calc. take place outside of psc relevant region.
        ! The following calculations take place outside of stratosphere.
        ! mz_pj_20070209-
        !----------------------------------------------------------------------

        ! ka_ok_20100118+
        IF (KinPAr) THEN
           DO jt=1, NSB
              idt = idt_HNO3_natsize(jt)
! op_pj_20170215+
              zxtte(jp,jk,idt) = - pxtte(_RI_X_ZN_(jp,jk,idt)) &
                   - pxtm1(_RI_X_ZN_(jp,jk,idt))/time_step_len  
! op_pj_20170215-
           END DO
        ENDIF
        ! ka_ok_20100118-

        !-----
        ! Set all psc channel objects outside of psc relevant region to zero.
        ! Note that the psc relevant region can vary, therefore, setting the
        ! channel objects zero only before the first time step would not
        ! suffice.
        !-----
        HNO3_liq(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
        HNO3_nat(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
        HNO3_gas(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
        HCl_liq(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
        HCl_gas(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
        HOCl_liq(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
        HOCl_gas(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
        HBr_liq(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
        HBr_gas(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
        HOBr_liq(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
        HOBr_gas(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp

        flt_phase(_RI_XYZ__(jp,jrow,jk))= 0.0_dp
        N_solid(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
! op_rd_20110304+
        IF (.NOT. l_feedback) clim_HNO3_nat(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
! op_rd_20110304-
        flt_val(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
        A_liq(_RI_XYZ__(jp,jrow,jk))    = 0.0_dp
        r_SurfMed(_RI_XYZ__(jp,jrow,jk))= 0.0_dp

        ! ka_ok_20100118+
        r_ice(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp ! op_pj_20100915
        IF (KinPar) THEN
           HNO3_nsarr(jp,jk,:)= 0.0_dp
           DO jt=1, NSB
              HNO3_natsize(jt)%ptr(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
           END DO
           N_ice(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           N_NAT(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
           r_NAT(_RI_XYZ__(jp,jrow,jk))  = 0.0_dp
        END IF
        ! ka_ok_20100118-
      
        IF (LCalcChem) THEN
           DO i=1, NKHET
              khet(i)%ptr(_RI_XYZ__(jp,jrow,jk))   = 0.0_dp
           END DO
        END IF   ! IF (LCalcChem) THEN ...
! mz_pj_20070209+
      !END IF   ! IF (val_psc(jp,jk)) THEN ...
      END IF if_val_strat   ! IF (val_strat(jp,jk)) THEN ...
! mz_pj_20070209-

! mz_pj_20070209+ moved from inside (see above)
      IF (.NOT. val_psc(jp,jk)) THEN
         !-----
         ! Set HNO3_nat outside of psc relevant region to zero. 
         ! The molecules transported as nat-tracer are added to the 
         ! corresponding gas phase tracer. 
         !-----
         idt  = idt_HNO3
         idt2 = idt_HNO3_nat
! op_pj_20170215+
         zxtte(jp,jk,idt) = zxtte(jp,jk,idt) + &
              ( zxtte(jp,jk,idt2) + pxtte(_RI_X_ZN_(jp,jk,idt2)) ) + &
              pxtm1(_RI_X_ZN_(jp,jk,idt2))/time_step_len
         zxtte(jp,jk,idt2) = - pxtte(_RI_X_ZN_(jp,jk,idt2))          &
              - pxtm1(_RI_X_ZN_(jp,jk,idt2))/time_step_len
! op_pj_20170215-
      END IF
! mz_pj_20070209-

    END DO vector_loop
  END DO level_loop

  !---------------------------------------------------------------------
  ! mz_rs_20060123+
  ! define rate coefficients for kpp
  khet_St_3d(ihs_N2O5_H2O)  %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(05)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(16)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(26)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_HOCl_HCl)  %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(03)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(14)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(23)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_ClNO3_HCl) %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(01)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(12)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(24)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_ClNO3_H2O) %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(02)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(13)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(25)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_N2O5_HCl)  %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(04)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(15)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_ClNO3_HBr) %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(06)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(17)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_BrNO3_HCl) %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(07)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(18)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_HOCl_HBr)  %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(08)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(19)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(29)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_HOBr_HCl)  %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(09)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(20)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(27)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_HOBr_HBr)  %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(10)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(21)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(28)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_BrNO3_H2O) %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(11)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(22)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(30)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_Hg) %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(31)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(32)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(33)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  khet_St_3d(ihs_RGM) %ptr(_RI_XYZ__(1:kproma,jrow,:)) = &
    khet(34)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(35)%ptr(_RI_XYZ__(1:kproma,jrow,:)) + &
    khet(36)%ptr(_RI_XYZ__(1:kproma,jrow,:))
  ! mz_rs_20060123-
  !---------------------------------------------------------------------

  ! -------------------------------------------------------
  ! start of PSC Sedimentation (only for PSC region, i.e. val_psc=TRUE)
  ! -------------------------------------------------------

! op_rd_20110304+
  !---
  ! l_feedback==false and KinPar==false:
  !     calculate xite_3d for tloop_index=1 
  !     calculate pxtte(:,:,idt_HNO3_nat) for tloop_index=2
  !
  ! l_feedback==true and KinPar==false:
  !     calculate xite_3d and pxtte(:,:,idt_HNO3_nat)
  !
  ! l_feedback==true and KinPar==true:
  !     separate calculation of pxtte(:,:,idt_HNO3_nat) 
  !
  ! l_feedback==false and KinPar==true:
  !     untested
  !---
  DO tloop_index = 1,tloop    
    IF (.NOT. l_feedback .AND. tloop_index .EQ. 1) THEN
      DO jk=1,nlev
        DO jp=1,kproma
          IF (val_psc(jp,jk)) THEN
            r_solid(_RI_XYZ__(jp,jrow,jk)) = r_s_offl(jp,jk)
          ELSE
            r_solid(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO jk=1,nlev
        DO jp=1,kproma
          IF (val_psc(jp,jk)) THEN
             IF (kinpar) THEN
                r_solid(_RI_XYZ__(jp,jrow,jk)) = r_ice(_RI_XYZ__(jp,jrow,jk))
             ELSE
                r_solid(_RI_XYZ__(jp,jrow,jk)) = r_s(jp,jk)
             END IF
          ELSE
            r_solid(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
          ENDIF
        ENDDO
      ENDDO
    ENDIF
 ! op_rd_20110304-

  !-----
  ! PSC Sedimentation:
! ka_ok_20100118+
!!$  ! 1.) Calculation of sedimentation velocity: v_sed (unit: m/s)
  ! 1.) Calculation of sedimentation velocity: v_sed_ice and v_sed_NAT
  !     (unit: m/s)
! ka_ok_20100118-
  ! 2.) Calculation of sedimentation distance, i.e. the product of
  !   v_sed and time_step_len converted into pressure units: SedStep
  !   (unit: Pa)
! ka_ok_20100118+
  !   -> IF (KinPar) Then seperate for ice and every NAT-sizebin
  ! If KINPAR=.TRUE. then v_sed_ice corresponds to r_ice
  ! and v_sed_NAT to r_NAT; also SedStep_ice to r_ice (v_sed_ice) and
  ! SedStep_NAT to r_NAT (v_sed_NAT).
  ! If KINPAR=.FALSE. then v_sed_ice corresponds to r_solid (r_ice is r_solid)
  ! and v_sed_NAT to r_solid; also SedStep_ice to r_solid (v_sed_ice)
  ! and SedStep_NAT to r_solid (v_sed_NAT).
! ka_ok_20100118-
  !-----
  DO jk=1,nlev
     DO jp=1,kproma
! ka_ok_20100118+
        IF (val_psc(jp,jk)) THEN
           v_sed_ice(_RI_XYZ__(jp,jrow,jk)) = mz_psc_vel(press_3d(_RI_XYZ__(jp,jrow,jk)), &
                TEMP(jp,jk), &
! op_rd_20110304+
                r_solid(_RI_XYZ__(jp,jrow,jk)))
! op_rd_20110304-
           SedStep_ice(jp,jk)    = mz_psc_SedStep(time_step_len, &
                press_3d(_RI_XYZ__(jp,jrow,jk)), TEMP(jp,jk), &
                v_sed_ice(_RI_XYZ__(jp,jrow,jk)))
        ELSE
           v_sed_ice(_RI_XYZ__(jp,jrow,jk)) = 0.0_dp
           SedStep_ice(jp,jk)    = 0.0_dp
        ENDIF
! ka_ok_20100118-
     ENDDO
  ENDDO
! ka_ok_20100118+
  IF (KinPar) THEN
     DO jt=1, NSB
        DO jk=1,nlev
           DO jp=1,kproma
              IF (val_psc(jp,jk)) THEN
                 v_sed_NAT_arr(jp,jk,jt) = &
                      mz_NAT_vel(time_step_len,PRESS(jp,jk), &
                      TEMP(jp,jk),r_NAT_arr(jp,jk,jt), &
                      G_arr(jp,jk,jt), N_NAT_arr(jp,jk,jt))
                 SedStep_NAT_arr(jp,jk,jt) = &
                      mz_psc_SedStep(time_step_len, &
                      press_3d(_RI_XYZ__(jp,jrow,jk)), TEMP(jp,jk), &
                      v_sed_NAT_arr(jp,jk,jt))
              ELSE
                 v_sed_NAT_arr(jp,jk,jt)=0.0_dp
                 SedStep_NAT_arr(jp,jk,jt)=0.0_dp
              ENDIF
           END DO
        END DO
     END DO
     v_sed_NAT(_RI_XYZ__(1:kproma,jrow,:))= &
          SUM(v_sed_NAT_arr(1:kproma,:,:),DIM=3)/REAL(NSB,DP)
  ELSE
     v_sed_NAT(_RI_XYZ__(1:kproma,jrow,:))=v_sed_ice(_RI_XYZ__(1:kproma,jrow,:))
  END IF
! ka_ok_20100118-
  !-----
  ! PSC Sedimentation:
  ! 3.) SedStep must not exceed the height of the next lower grid box and
  !   is, therefore, reduced if necessary.
  !-----
  DO jk=1,nlev-1
     DO jp=1,kproma
        IF( val_psc(jp,jk) ) THEN
           SedStep_ice(jp,jk)=MIN(SedStep_ice(jp,jk), &
                pressi_3d(_RI_XYZ__(jp,jrow,jk+2))-pressi_3d(_RI_XYZ__(jp,jrow,jk+1)))
        ENDIF
     ENDDO
  ENDDO

  IF (KinPar) THEN
     DO jt=1, NSB
        DO jk=1,nlev-1
           DO jp=1,kproma
              IF( val_psc(jp,jk) ) THEN
                 SedStep_NAT_arr(jp,jk,jt)=MIN(SedStep_NAT_arr(jp,jk,jt), &
                      pressi_3d(_RI_XYZ__(jp,jrow,jk+2))-pressi_3d(_RI_XYZ__(jp,jrow,jk+1)))
              END IF
           END DO
        END DO
     END DO
  ENDIF
! ka_ok_20100118-

  !-----
  ! PSC Sedimentation:
  ! 4.) Calculation of change of mass fraction of ice in air due to psc
  !   sedimentation
  ! 5.) Calculation of change of amount-of-substance fraction of HNO3 in air
  !   due to psc sedimentation
  !-----
  IF (l_feedback .OR. tloop_index .EQ. 1) THEN ! op_rd_20110304
  IceChangeDueToSed(1:kproma,:)=mz_psc_sed(                        &
       kproma,nlev,                                                &
       pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)),                               &
       pressi_3d(_RI_XYZ__(1:kproma,jrow,2:nlev+1)),                            &
! ka_ok_20100118+
       SedStep_ice(1:kproma,:),                                    &
! ka_ok_20100118-
! op_ff_20170217+
       xim1_3d(_RI_XYZ__(1:kproma,jrow,:)) + &
       (xite_3d(_RI_XYZ__(1:kproma,jrow,:))+zhumte(:,:,2))*time_step_len, &
! op_ff_20170217-
       val_psc(1:kproma,:))
! op_rd_20110304+
! op_pj_20170215+
      zhumte(:,:,2) = zhumte(:,:,2) + &
           IceChangeDueToSed(1:kproma,:)/time_step_len
! op_pj_20170215-
    ENDIF
! op_rd_20110304-

  !-----
  ! PSC Sedimentation:
  ! 6.) Calculation of change of amount-of-substance fraction of HNO3 in air
  !   due to psc sedimentation
  ! 7.) amount-of-substance fraction of HNO3 in air is modified due to
  !   psc sedimentation
  !-----

  IF (.NOT. KinPar) THEN ! ka_ok_20100118
     IF (l_feedback .OR. tloop_index .EQ. 2) THEN ! op_rd_20110304
     HNO3ChangeDueToSed(1:kproma,:)=mz_psc_sed(                       &
          kproma,nlev,                                                &
          pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)),                                &
          pressi_3d(_RI_XYZ__(1:kproma,jrow,2:nlev+1)),                            &
! ka_ok_20100118+
          SedStep_ice(1:kproma,:),                                    &
! ka_ok_20100118-
          HNO3_nat(_RI_XYZ__(1:kproma,jrow,:)),                                  &
          val_psc(1:kproma,:))
! op_rd_20110211+
      idt = idt_HNO3_nat
! op_pj_20170215+
      zxtte(1:kproma,:,idt) = zxtte(1:kproma,:,idt) + &
           HNO3ChangeDueToSed(1:kproma,:)/time_step_len
! op_pj_20170215-
     END IF
! op_rd_20110211-
! ka_ok_20100118+
  ELSE !IF (.NOT. KinPar)
     DO jt=1, NSB
        HNO3ChangeDueToNATSed_arr(1:kproma,:,jt)=mz_psc_sed(kproma, nlev, &
             pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)),                                &
             pressi_3d(_RI_XYZ__(1:kproma,jrow,2:nlev+1)),                             &
             SedStep_NAT_arr(1:kproma,:,jt),                              &
             HNO3_natsize(jt)%ptr(_RI_XYZ__(1:kproma,jrow,:)),                      &
             val_psc(1:kproma,:))
     END DO
     HNO3ChangeDueToSed(1:kproma,:)=&
          SUM(HNO3ChangeDueToNATSed_arr(1:kproma,:,:),DIM=3)
     DO jt=1, NSB
        idt = idt_HNO3_natsize(jt)
! op_pj_20170215+
        zxtte(1:kproma,:,idt) = zxtte(1:kproma,:,idt) + &
             HNO3ChangeDueToNATSed_arr(1:kproma,:,jt)/time_step_len
! op_pj_20170215-
     END DO
! op_rd_20110304+
     idt = idt_HNO3_nat
! op_pj_20170215+
     zxtte(1:kproma,:,idt) = zxtte(1:kproma,:,idt) + &
          HNO3ChangeDueToSed(1:kproma,:)/time_step_len
! op_pj_20170215-
! op_rd_20110304-
  END IF
! ka_ok_20100118-
  ENDDO !DO tloop_index = 1,tloop ! op_rd_20110304

  ! ka_ok_20100118+
  IF (KinPar) THEN
     HNO3_Change_Sed(_RI_XYZ__(1:kproma,jrow,:)) = HNO3ChangeDueToSed(1:kproma,:)
     SedStep_Ice_noarr(_RI_XYZ__(1:kproma,jrow,:)) = SedStep_ice(1:kproma,:) 
     Ice_Change_Sed_noarr(_RI_XYZ__(1:kproma,jrow,:)) = IceChangeDueToSed(1:kproma,:)
  END IF
  ! ka_ok_20100118-

  !op_kg_20091112+
  level_loop2: DO jk=1,nlev
     vector_loop2: DO jp=1,kproma
        if (HNO3ChangeDueToSed(jp,jk) < 0.0_dp) then
           loss_nat(_RI_XYZ__(jp,jrow,jk))=HNO3ChangeDueToSed(jp,jk)/time_step_len
           prod_nat(_RI_XYZ__(jp,jrow,jk))=0.0_dp
        else
           prod_nat(_RI_XYZ__(jp,jrow,jk))=HNO3ChangeDueToSed(jp,jk)/time_step_len
           loss_nat(_RI_XYZ__(jp,jrow,jk))=0.0_dp
        end if
     END DO vector_loop2
  END DO level_loop2
  !op_kg_20091112-

  ! -------------------------------------------------------
  ! end of PSC Sedimentation (only for PSC region)
  ! -------------------------------------------------------

  ! op_pj_20170215+
#ifndef MESSYTENDENCY
  xlte_3d(_RI_XYZ__(1:kproma,jrow,:)) = &
       xlte_3d(_RI_XYZ__(1:kproma,jrow,:)) + zhumte(:,:,1)
  xite_3d(_RI_XYZ__(1:kproma,jrow,:)) = &
       xite_3d(_RI_XYZ__(1:kproma,jrow,:)) + zhumte(:,:,2)
  qte_3d(_RI_XYZ__(1:kproma,jrow,:))  = &
       qte_3d(_RI_XYZ__(1:kproma,jrow,:))  + zhumte(:,:,3)
#else
  CALL mtend_add_l(my_handle, mtend_id_xl, px=zhumte(:,:,1))
  CALL mtend_add_l(my_handle, mtend_id_xi, px=zhumte(:,:,2))
  CALL mtend_add_l(my_handle, mtend_id_q,  px=zhumte(:,:,3))
#endif

  ! op_pj_20171127+
  ! ADJUST LARGE SCALE CLOUD COVER
#ifndef MESSYTENDENCY
  xlp1(1:kproma,:) = xlm1_3d(_RI_XYZ__(1:kproma,jrow,:))+zhumte(:,:,1)*time_step_len
  xip1(1:kproma,:) = xim1_3d(_RI_XYZ__(1:kproma,jrow,:))+zhumte(:,:,2)*time_step_len
#else
  CALL mtend_get_start_l(mtend_id_xl, v0=xlp1)
  CALL mtend_get_start_l(mtend_id_xi, v0=xip1)
#endif
  ! see cloud_ori
  DO jk=1, nlev
     DO jp=1, kproma
        lo  = (xlp1(jp,jk) < ccwmin)
        lo1 = (xip1(jp,jk) < ccwmin)
        aclc(_RI_XYZ__(jp,jrow,jk)) = MERGE(0.0_dp,aclc(_RI_XYZ__(jp,jrow,jk)),lo.AND.lo1)
     END DO
  END DO
  ! op_pj_20171127-

#ifdef MECCA_TAG
  CALL mecca_tag_calc_xtte4scav(zxtte(:,:,:), 1, pxtp1, kproma)
#endif

#ifndef MESSYTENDENCY
! WITHOUT TENDENCY

#ifdef MECCA_TAG
  ! ... all tracers if MECCA_TAG
  DO idt=1, ntrac_gp
     pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  END DO
#else
  ! ... specific tracers only, without MECCA_TAG
  idt=idt_HNO3
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_HCl
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_HOCl
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_ClNO3
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_HBr
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_HOBr
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_BrNO3
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  !
  idt=idt_HNO3_nat
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  !
  IF (KinPar) THEN
     DO jt=1, NSB
        idt = idt_HNO3_natsize(jt)
        pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
     END DO
  END IF
#endif

#else
! WITH TENDENCY

#ifdef MECCA_TAG
  ! ... all tracers if MECCA_TAG
  CALL mtend_add_l(my_handle, mtend_id_tracer, pxt=zxtte(:,:,:))
#else
  ! ... specific tracers only, without MECCA_TAG
  CALL mtend_add_l(my_handle,idt_HNO3,   px=zxtte(:,:,idt_HNO3))
  CALL mtend_add_l(my_handle, idt_HCl,   px=zxtte(:,:,idt_HCl))
  CALL mtend_add_l(my_handle,idt_HOCl,   px=zxtte(:,:,idt_HOCl))
  CALL mtend_add_l(my_handle, idt_ClNO3, px=zxtte(:,:,idt_ClNO3))
  CALL mtend_add_l(my_handle, idt_HBr,   px=zxtte(:,:,idt_HBr))
  CALL mtend_add_l(my_handle, idt_HOBr , px=zxtte(:,:,idt_HOBr))
  CALL mtend_add_l(my_handle, idt_BrNO3, px=zxtte(:,:,idt_BrNO3))
  !
  CALL mtend_add_l(my_handle, idt_HNO3_nat, px=zxtte(:,:,idt_HNO3_nat))
  !
  IF (KinPar) THEN
     DO jt=1, NSB
        idt = idt_HNO3_natsize(jt)
        CALL mtend_add_l(my_handle, idt, px=zxtte(:,:,idt))
     END DO
  END IF
#endif
  ! MECCA_TAG

#endif
  ! MESSYTENDENCY
  ! op_pj_20170215-

END SUBROUTINE msbm_physc

!=============================================================================

SUBROUTINE msbm_local_end
  !-----
  ! For the next advection step, add HNO3_liq, HCl_liq, HOCl_liq, HBr_liq,
  !  HOBr_liq to gas phase HNO3, HCl, HOCl, HBr, HOBr. 
  !-----

  USE messy_main_grid_def_mem_bi,ONLY: jrow, kproma, nproma, nlev
  USE messy_main_timer,         ONLY:  time_step_len
  USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, ntrac_gp
! op_pj_20170215+
#ifdef MECCA_TAG
  USE messy_main_tracer_mem_bi, ONLY: pxtm1 => qxtm1
#endif
! op_pj_20170215-

  IMPLICIT NONE
  
  INTEGER :: idt
! op_pj_20170215+
  REAL(DP) :: zxtte(nproma,nlev,ntrac_gp)
#ifdef MECCA_TAG
  INTEGER :: jt, jp
  REAL(DP) :: pxtp1(nproma,nlev,ntrac_gp)
#endif

  zxtte(:,:,:) = 0.0_dp
#ifdef MECCA_TAG
  pxtp1(:,:,:) = 0.0_dp
  DO jt=1, ntrac_gp
     DO jp=1, kproma
        pxtp1(jp,:,jt) = pxtm1(_RI_X_ZN_(jp,:,jt)) + pxtte(_RI_X_ZN_(jp,:,jt))*time_step_len
     END DO
  END DO
#endif
! op_pj_20170215-

! op_pj_20170215+
  zxtte(1:kproma,:,idt_HNO3) = HNO3_liq(_RI_XYZ__(1:kproma,jrow,:))/time_step_len
  zxtte(1:kproma,:,idt_HCl)  = HCl_liq(_RI_XYZ__(1:kproma,jrow,:))/time_step_len
  zxtte(1:kproma,:,idt_HOCl) = HOCl_liq(_RI_XYZ__(1:kproma,jrow,:))/time_step_len
  zxtte(1:kproma,:,idt_HBr)  = HBr_liq(_RI_XYZ__(1:kproma,jrow,:))/time_step_len
  zxtte(1:kproma,:,idt_HOBr) = HOBr_liq(_RI_XYZ__(1:kproma,jrow,:))/time_step_len
! op_pj_20170215-

#ifdef MECCA_TAG
  CALL mecca_tag_calc_xtte4scav(zxtte(:,:,:), 1, pxtp1, kproma)
#endif

#ifndef MESSYTENDENCY
! WITHOUT TENDENCY

#ifdef MECCA_TAG
  ! ... all tracers if MECCA_TAG
  DO idt=1, ntrac_gp
     pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  END DO
#else
  ! ... specific tracers only, without MECCA_TAG
  idt=idt_HNO3
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_HCl
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_HOCl
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_HBr
  pxtte(_RI_X_ZN_(1:kproma,:,idt))= pxtte(_RI_X_ZN_(1:kproma,:,idt)) + zxtte(1:kproma,:,idt)
  idt=idt_HOBr
#endif

#else
! WITH TENDENCY

#ifdef MECCA_TAG
  ! ... all tracers if MECCA_TAG
  CALL mtend_add_l(my_handle, mtend_id_tracer, pxt=zxtte(:,:,:))
#else
  ! ... specific tracers only, without MECCA_TAG
  CALL mtend_add_l(my_handle,idt_HNO3, px=zxtte(:,:,idt_HNO3))
  !
  CALL mtend_add_l(my_handle,idt_HCl,  px=zxtte(:,:,idt_HCl))
  !
  CALL mtend_add_l(my_handle,idt_HOCl, px=zxtte(:,:,idt_HOCl))
  !
  CALL mtend_add_l(my_handle,idt_HBr,  px=zxtte(:,:,idt_HBr))
  !
  CALL mtend_add_l(my_handle,idt_HOBr, px=zxtte(:,:,idt_HOBr))
#endif

! TENDENCY
#endif

END SUBROUTINE msbm_local_end

!=============================================================================

SUBROUTINE msbm_read_nml_cpl(status, iou)
  !------------------------------------------
  ! read coupling namelist /CPL/ from msbm.nml
  !------------------------------------------

  ! MESSy
  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  
  IMPLICIT NONE
  
  INTEGER, INTENT(out) :: status
  INTEGER, INTENT(in) :: iou
  
  CHARACTER(len=*), PARAMETER :: substr='msbm_read_nml_cpl'
  LOGICAL :: lex     ! file exists?
  INTEGER :: fstat   ! file status
  
  ! mz_pj_20070209 Tropop_Channel, Tropop_Index added
  ! op_pj_20091013 r_lb added
  NAMELIST /CPL/ LCalcChem, TempShift, Tropop_Channel, Tropop_Index &
       , r_lat, r_lb, r_mb, r_ub  &  ! op_pj_20091013
       , l_feedback    &             ! op_rd_20100108
       , c_H2SO4clim, c_predef_HNO3_tot

  status = 1             ! initialise status flag with error code
  
  CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
  IF (.not.lex) RETURN   ! error: msbm.nml does not exist
  read(iou, nml=cpl, iostat=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
  IF (fstat/=0) RETURN   ! error while reading namelist
  
  ! op_rd_20100108+
  IF (l_feedback) THEN
     write(*,*) 'feedback of PSC on dynamics is ON'
  ELSE
     write(*,*) 'feedback of PSC on dynamics is OFF'
  ENDIF
  ! op_rd_20100108-
  write (*,*) 'switch for calculating psc chemistry: LCalcChem = ', LCalcChem
  write (*,*) 'temperature shift in msbm submodel: TempShift = ', TempShift
  ! mz_pj_20070209+
  write (*,*) 'tropopause channel:',TRIM(Tropop_Channel)
  write (*,*) 'tropopause index  :',TRIM(Tropop_Index)
  ! mz_pj_20070209-

  write (*,*) 'H2SO4 channel / object:',TRIM(c_H2SO4clim%CHA), &
       ' / ',TRIM(c_H2SO4clim%OBJ)
  write (*,*) 'predefined HNO3_NAT channel / object:' &
       ,TRIM(c_predef_HNO3_tot%CHA),' / ',TRIM(c_predef_HNO3_tot%OBJ)
  
  CALL read_nml_close(substr, iou, modstr)
  status = 0   ! no error
  
END SUBROUTINE msbm_read_nml_cpl

!=============================================================================

END MODULE messy_msbm_si

