#include "messy_main_ppd_bi.inc"

! ***********************************************************************
MODULE messy_e4chem_si
  ! ***********************************************************************

  ! ***********************************************************************
  ! MESSy-SMIL FOR SUBMODEL E4CHEM
  !
  ! Author: Andreas Baumgaertner, MPICH, 2009
  ! Based on MESSy-SMIL for MECCA
  ! ***********************************************************************


  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      info_bi, warning_bi, error_bi
  USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp
#ifdef MESSYTENDENCY
  !tendency budget
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
    mtend_get_start_l,      &
    mtend_add_l,            &
    mtend_register,         &
    mtend_id_tracer,        &
    mtend_id_q, mtend_id_t, mtend_id_xi
#endif
  ! MESSy
  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY, PTR_1D_ARRAY, str
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT &
                                    , t_chaobj_cpl ! op_pj_20111201
  USE messy_cmn_photol_mem      ! IP_MAX, ip_*, jname
  ! SUBMODEL E4CHEM
  USE messy_e4chem

  IMPLICIT NONE
  SAVE
  PRIVATE

  INTRINSIC :: NULL

  CHARACTER(LEN=STRLEN_MEDIUM)  :: photrat_channel_gp = ''
  TYPE(T_CHAOBJ_CPL)            :: tropop_idx
  TYPE(T_CHAOBJ_CPL)            :: lcover
  TYPE(T_CHAOBJ_CPL)            :: ratep
  TYPE(T_CHAOBJ_CPL)            :: prec
  TYPE(T_CHAOBJ_CPL)            :: lwc
  TYPE(T_CHAOBJ_CPL)            :: cvprec

  ! op_pj_20160315+
  INTEGER :: i_H2O_feedback = 0 ! no feedback to hydrological cycle
  !                             ! (1: feedback from GP)
  INTEGER :: i_ice_feedback = 1 ! 1: from GP, 0: no feedback
  ! op_pj_20160315-

  ! 2d-field for tropopause index (details via namelist)
  REAL(dp), DIMENSION(:,:), POINTER :: tp_i => NULL()

  ! PSC-region indicator, needed by CLOUD submodel
  REAL(dp), DIMENSION(:,:,:), POINTER :: flt_pscreg => NULL()

  ! ice 
  REAL(dp), DIMENSION(:,:,:), POINTER :: ice => NULL()

  ! cos zenith angle
  TYPE(T_CHAOBJ_CPL)                :: cossza
  REAL(dp), DIMENSION(:,:), POINTER :: cossza_2d     => NULL()

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  NAMELIST /CPL/ photrat_channel_gp, tropop_idx & ! op_pj_20111201
       , lcover, ratep, prec, lwc, cvprec, cossza &
       , i_H2O_feedback, i_ice_feedback ! op_pj_20160315

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: e4chem_initialize
  PUBLIC :: e4chem_new_tracer    ! define tracers
  PUBLIC :: e4chem_init_memory
  PUBLIC :: e4chem_init_tracer   ! initialize tracers
  PUBLIC :: e4chem_init_coupling
  PUBLIC :: e4chem_local_start
  PUBLIC :: e4chem_physc
  PUBLIC :: e4chem_free_memory

  ! pointer to channel objects (photolysis rates)
  TYPE(PTR_3D_ARRAY), PUBLIC, DIMENSION(IP_MAX) :: photrat_gp
  REAL(dp), DIMENSION(:,:), POINTER :: DANI, DANIM
  ! CONVECTIVE CLOUD BASE LEVEL
  REAL(DP) , POINTER :: KCONBOT(:,:)        => NULL()
  ! FOR LARGE SCALE CLOUD COVER
  REAL (DP), POINTER :: PCLCOVER(:,:,:)     => NULL()   
  !  PRECIPITATION FORMATION RATE     [KG KG-1]
  REAL (DP), POINTER :: PMRATEP(:,:,:)      => NULL()   
  ! RECIPITATION FLUX               [KG M-2 S-1]
  REAL (DP), POINTER :: PFPREC(:,:,:)       => NULL()   
  ! LIQUID WATER CONTENT             [KG KG-1]
  REAL (DP), POINTER :: PMLWC(:,:,:)        => NULL()
  ! FRESHLY FORMED PRECIPITATION
  REAL (DP), POINTER :: PCVDPREC(:,:,:)     => NULL()
  ! mz_bk_20101222+
  ! ozone molecules destroyed by BrO+ClO cycle
  REAL(DP), POINTER :: ZDELTAO3_BRV(:,:,:)  => NULL()
  ! ozone budget
  REAL(DP), POINTER :: ZPRODO2(:,:,:)       => NULL()
  REAL(DP), POINTER :: ZPRODCO(:,:,:)       => NULL()
  REAL(DP), POINTER :: ZPRODCH4(:,:,:)      => NULL()
  REAL(DP), POINTER :: ZDESTH12(:,:,:)      => NULL()
  REAL(DP), POINTER :: ZDESTH14(:,:,:)      => NULL()
  REAL(DP), POINTER :: ZDESTN13(:,:,:)      => NULL()
  REAL(DP), POINTER :: ZDESTC1(:,:,:)       => NULL()
  REAL(DP), POINTER :: ZDESTCL2O2(:,:,:)    => NULL()
  REAL(DP), POINTER :: ZDESTCLOH(:,:,:)     => NULL()
  REAL(DP), POINTER :: ZDESTH8(:,:,:)       => NULL()
  REAL(DP), POINTER :: ZDESTH4(:,:,:)       => NULL()
  REAL(DP), POINTER :: ZDESTH11(:,:,:)      => NULL()
  REAL(DP), POINTER :: ZDESTO1(:,:,:)       => NULL()  
  ! mz_bk_20101222-

  INTEGER :: zidt_H2O

  INTEGER :: idt_H          = 0
  INTEGER :: idt_OH         = 0
  INTEGER :: idt_HO2        = 0
  INTEGER :: idt_N          = 0
  INTEGER :: idt_NO         = 0
  INTEGER :: idt_NO2        = 0
  INTEGER :: idt_NO3        = 0
  INTEGER :: idt_N2O5       = 0
  INTEGER :: idt_HNO4       = 0
  INTEGER :: idt_CL         = 0
  INTEGER :: idt_CLO        = 0
  INTEGER :: idt_HOCl       = 0 ! CLOH  --> HOCl
  INTEGER :: idt_CL2O2      = 0
  INTEGER :: idt_CL2        = 0
  INTEGER :: idt_HCHO       = 0 ! CH2O --> HCHO
  INTEGER :: idt_CH3O2      = 0
  INTEGER :: idt_CH4        = 0
  INTEGER :: idt_N2O        = 0
  INTEGER :: idt_H2O2       = 0
  INTEGER :: idt_HCl        = 0
  INTEGER :: idt_CO         = 0
  INTEGER :: idt_CH3OOH     = 0 ! CH3O2H --> CH3OOH
  INTEGER :: idt_ClNO3      = 0
  INTEGER :: idt_CFCl3      = 0
  INTEGER :: idt_CF2Cl2     = 0
  INTEGER :: idt_CH3CL      = 0
  INTEGER :: idt_CCL4       = 0
  INTEGER :: idt_CH3CCL3    = 0
  INTEGER :: idt_H2         = 0
  INTEGER :: idt_HNO3       = 0
  INTEGER :: idt_OHAB       = 0
  INTEGER :: idt_HO2AB      = 0
  INTEGER :: idt_NAT        = 0
  INTEGER :: idt_O3P        = 0
  INTEGER :: idt_O3         = 0
  INTEGER :: idt_O1D        = 0
  INTEGER :: idt_CO2        = 0

CONTAINS

  ! ========================================================================
  SUBROUTINE e4chem_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'e4chem_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    CALL start_message_bi(modstr,'INITIALISATION',modstr)

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL e4chem_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('error in e4chem_read_nml_ctrl',substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast(NSTCHPH, p_io)
    CALL p_bcast(l_fastscav, p_io)
    CALL p_bcast(l_Brparam, p_io) ! mz_ab_20100505

    ! CPL NAMELIST
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL e4chem_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('error in e4chem_read_nml_cpl',substr)
    ENDIF
    CALL p_bcast(i_H2O_feedback, p_io) ! op_pj_20160315
    CALL p_bcast(i_ice_feedback, p_io) ! op_pj_20160315
    CALL p_bcast(photrat_channel_gp, p_io)
    CALL p_bcast(tropop_idx%cha,p_io)
    CALL p_bcast(tropop_idx%obj,p_io)
    CALL p_bcast(ratep%cha, p_io)
    CALL p_bcast(ratep%obj, p_io)
    CALL p_bcast(prec%cha, p_io) 
    CALL p_bcast(prec%obj, p_io) 
    CALL p_bcast(lwc%cha, p_io)
    CALL p_bcast(lwc%obj, p_io)
    CALL p_bcast(cvprec%cha, p_io) 
    CALL p_bcast(cvprec%obj, p_io) 
    CALL p_bcast(cossza%cha, p_io) 
    CALL p_bcast(cossza%obj, p_io) 

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    CALL end_message_bi(modstr,'INITIALISATION',modstr)

  END SUBROUTINE e4chem_initialize
  ! ========================================================================

  !***************************************************************************

  SUBROUTINE e4chem_new_tracer

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,          ONLY: new_tracer, AIR, ON, OFF, &
         set_tracer, &
         AMOUNTFRACTION, &
         I_ADVECT, I_CONVECT, &
         I_VDIFF, &
         I_DRYDEP, &
         I_SCAV, &
         R_MOLARMASS, R_PSS  , &
         R_VINI
    USE messy_main_constants_mem, ONLY: MH, MC, MN, MO, MS, MCl, MF

    IMPLICIT NONE

    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'e4chem_new_tracer'
    CHARACTER(LEN=*), PARAMETER :: setname = GPTRSTR

    CALL new_tracer(status, setname, 'H',          modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_H)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_H, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H, R_molarmass         , MH)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'OH',          modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_OH)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_OH, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OH, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OH, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OH, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OH, R_molarmass         , MO+MH)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OH, R_pss               , 25._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OH, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'HO2',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_HO2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_HO2, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2, R_molarmass         , MH+MO*2._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'N',           modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_N)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_N, I_advect            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N, R_molarmass         , MN)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'NO',          modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_NO)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_NO, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO, R_molarmass         , MN+MO)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO, R_pss               , 2.E-3_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'NO2',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_NO2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_NO2, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO2, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO2, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO2, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO2, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO2, R_molarmass         , MN+MO*2._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO2, R_pss               , 1.0E-2_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO2, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'NO3',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_NO3)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_NO3, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO3, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO3, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO3, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO3, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO3, R_molarmass         , MN+MO*3._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO3, R_pss               , 1.8_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NO3, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'N2O5',        modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_N2O5)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_N2O5, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O5, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O5, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O5, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O5, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O5, R_molarmass         , MN*2.+MO*5._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O5, R_pss               , 1.E4_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O5, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'HNO4',        modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_HNO4)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_HNO4, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO4, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO4, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO4, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO4, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO4, R_molarmass         , MH+MN+MO*4._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO4, R_pss               , 1.E4_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO4, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'Cl',          modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_Cl)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_Cl, I_advect            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl, R_molarmass         , MCl)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'ClO',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_ClO)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_ClO, I_advect            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClO, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClO, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClO, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClO, R_molarmass         , MCl+MO)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClO, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'HOCl',        modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_HOCl)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_HOCl, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HOCl, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HOCl, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HOCl, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HOCl, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HOCl, R_molarmass         , MH+MO+MCl)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HOCl, R_pss               , 6.7e2_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HOCl, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'Cl2O2',       modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_Cl2O2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_Cl2O2, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2O2, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2O2, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2O2, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2O2, R_molarmass         , MCl*2.+MO*2.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2O2, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'Cl2',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_Cl2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_Cl2, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2, R_molarmass         , MCl*2._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2, R_pss               , 7.e-2_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_Cl2, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'HCHO',        modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_HCHO)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_HCHO, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCHO, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCHO, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCHO, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCHO, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCHO, R_molarmass         , MC+MH*2._dp+MO)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCHO, R_pss               , 3.2e3_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCHO, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CH3O2',       modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CH3O2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CH3O2, I_advect            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3O2, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3O2, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3O2, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3O2, R_molarmass         , MC+MH*3.+MO*2._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3O2, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CH4',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CH4)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CH4, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH4, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH4, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH4, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH4, R_molarmass         , MC+MH*4.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH4, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'N2O',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_N2O)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_N2O, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O, R_molarmass         , MN*2.+MO)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_N2O, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'H2O2',        modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_H2O2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_H2O2, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2O2, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2O2, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2O2, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2O2, I_scav              , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2O2, R_molarmass         , MH*2.+MO*2._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2O2, R_pss               , 7.45E4_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2O2, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'HCl',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_HCl)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_HCl, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCl, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCl, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCl, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCl, I_scav              , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCl, R_molarmass         , MH+MCl)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCl, R_pss               , 1.e14_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HCl, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CO',          modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CO)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CO, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO, R_molarmass         , MC+MO)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CH3OOH',      modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CH3OOH)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CH3OOH, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3OOH, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3OOH, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3OOH, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3OOH, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3OOH, R_molarmass         , MC+MH*4.+MO*2._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3OOH, R_pss               , 3.e2_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3OOH, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'ClNO3',       modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_ClNO3)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_ClNO3, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClNO3, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClNO3, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClNO3, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClNO3, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClNO3, R_molarmass         , MCl+MN+MO*3._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClNO3, R_pss               , 1.e30_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_ClNO3, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CFCl3',       modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CFCl3)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CFCl3, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CFCl3, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CFCl3, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CFCl3, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CFCl3, R_molarmass         , MC+MF+MCl*3._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CFCl3, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CF2Cl2',      modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CF2Cl2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CF2Cl2, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CF2Cl2, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CF2Cl2, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CF2Cl2, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CF2Cl2, R_molarmass         , MC+MF*2.+MCl*2.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CF2Cl2, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CH3Cl',       modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CH3Cl)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CH3Cl, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3Cl, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3Cl, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3Cl, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3Cl, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3Cl, R_molarmass         , MC+MH*3.+MCl)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3Cl, R_pss               , 9.4e-2_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3Cl, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CCl4',        modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CCl4)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CCl4, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CCl4, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CCl4, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CCl4, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CCl4, R_molarmass         , MC+MCl*4.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CCl4, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CH3CCl3',     modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CH3CCl3)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CH3CCl3, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3CCl3, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3CCl3, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3CCl3, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3CCl3, R_molarmass         , MC*2.+MH*3.+MCl*3.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CH3CCl3, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'H2',          modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_H2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_H2, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2, R_molarmass         , MH*2.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_H2, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'HNO3',        modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_HNO3)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_HNO3, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO3, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO3, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO3, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO3, I_scav              , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO3, R_molarmass         , MH+MN+MO*3.+MH*7.+MN+MO*6.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO3, R_pss               , 1.E4_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HNO3, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'OHAB',          modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_OHAB)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_OHAB, I_advect            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OHAB, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OHAB, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OHAB, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OHAB, R_molarmass         , MO+MH)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OHAB, R_pss               , 25._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_OHAB, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'HO2AB',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_HO2AB)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_HO2AB, I_advect            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2AB, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2AB, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2AB, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2AB, R_molarmass         , MH+MO*2.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_HO2AB, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'NAT',        modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_NAT)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_NAT, I_advect            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NAT, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NAT, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NAT, I_drydep            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NAT, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NAT, R_molarmass         , MH*7.+MN+MO*6.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_NAT, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'O3',          modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_O3)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_O3, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3, I_drydep            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3, R_molarmass         , MO*3.)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3, R_pss               , 0.01_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'O3P',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_O3P)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_O3P, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3P, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3P, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3P, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3P, R_molarmass         , MO)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O3P, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'O1D',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_O1D)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_O1D, I_advect            , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O1D, I_convect           , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O1D, I_vdiff             , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O1D, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O1D, R_molarmass         , MO)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_O1D, R_VINI              , 1E-24_dp)
    CALL tracer_halt(substr, status)
    CALL new_tracer(status, setname, 'CO2',         modstr, &
         quantity             = AMOUNTFRACTION, &
         unit                 = 'mol/mol', &
         medium               = AIR, &
         idx                  = idt_CO2)
    CALL tracer_halt(substr,status)
    CALL set_tracer(status, setname, idt_CO2, R_vini              , 3.5E-04_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO2, I_advect            , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO2, I_convect           , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO2, I_vdiff             , ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO2, I_scav              , OFF)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, setname, idt_CO2, R_molarmass         , MC+MO*2.)
    CALL tracer_halt(substr, status)

    CALL end_message_bi(modstr, 'TRACER INITIALIZATION', substr)

  END SUBROUTINE e4chem_new_tracer

  !***************************************************************************

  SUBROUTINE e4chem_init_memory

    USE messy_main_grid_def_bi,      ONLY: philat_2d
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, npromz, ngpblks
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_2D_HORIZONTAL
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                           new_attribute
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'e4chem_init_memory'
    INTEGER :: status

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ALLOCATE(RCGAS(NUMTEM, NUMRAT))
    ALLOCATE(sulook(ksul, nproma, ngpblks))

    ! PRECALCULATE REACTION RATES (LOOKUP TABLE)
    CALL INRCGAS

    ! CALCULATE LOOKUP TABLE FOR AEROSOL SURFACE AREA
#ifndef CESM1
    CALL inisulnew(philat_2d(:,:),nproma,ngpblks,npromz=npromz)
#else
    CALL inisulnew(philat_2d(:,:),nproma,ngpblks,vnpromz=npromz)
#endif

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    !-----
    ! psc region indicator
    !-----
    CALL new_channel_object(status, modstr, 'PSC_region', p3=flt_pscreg &
         , lrestreq=.false. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PSC_region' &
         , 'long_name', c='flag indicating PSC relevant region')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'PSC_region', 'units', c=' ')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'ice', p3=ice &
         , lrestreq=.false. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ice', 'units', c='mol/mol')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DANIM', p2=DANIM &
         , reprid=GP_2D_HORIZONTAL, lrestreq=.true. )
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr, 'DANI', p2=DANI &
         , reprid=GP_2D_HORIZONTAL, lrestreq=.false. )
    CALL channel_halt(substr, status)

    ! mz_bk_20101222+
    ! ozone destroyed by bromine
    CALL new_channel_object(status, modstr, 'ZDELTAO3_BRV', p3=ZDELTAO3_BRV &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDELTAO3_BRV' &
         , 'long_name', c='ozone molecules destroyed by BrO+ClO cycle')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDELTAO3_BRV' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    ! ozone budget
    CALL new_channel_object(status, modstr, 'ZPRODO2', p3=ZPRODO2 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZPRODO2' &
         , 'long_name', c='ozone production by photolysis of O2')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZPRODO2' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZPRODCO', p3=ZPRODCO &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZPRODCO' &
         , 'long_name', c='ozone production by oxidation of CO')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZPRODCO' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZPRODCH4', p3=ZPRODCH4 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZPRODCH4' &
         , 'long_name', c='ozone production by oxidation of CH4')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZPRODCH4' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTH12', p3=ZDESTH12 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH12' &
         , 'long_name', c='ozone destruction by O3+HO2')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH12' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTH14', p3=ZDESTH14 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH14' &
         , 'long_name', c='ozone destruction by O3+OH')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH14' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTN13', p3=ZDESTN13 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTN13' &
         , 'long_name', c='ozone destruction by NO2+O(3P)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTN13' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTC1', p3=ZDESTC1 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTC1' &
         , 'long_name', c='ozone destruction by ClO+O(3P)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTC1' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTCL2O2', p3=ZDESTCL2O2 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTCL2O2' &
         , 'long_name', c='ozone destruction by Cl2O2 photolysis')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTCL2O2' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTCLOH', p3=ZDESTCLOH &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTCLOH' &
         , 'long_name', c='ozone destruction by ClOH photolysis')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTCLOH' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTH8', p3=ZDESTH8 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH8' &
         , 'long_name', c='ozone destruction by H2O+O(1D)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH8' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTH4', p3=ZDESTH4 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH4' &
         , 'long_name', c='ozone destruction by OH+O(3P)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH4' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTH11', p3=ZDESTH11 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH11' &
         , 'long_name', c='ozone destruction by O(3P)+HO2')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTH11' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    CALL new_channel_object(status, modstr, 'ZDESTO1', p3=ZDESTO1 &
         , reprid=GP_3D_MID, lrestreq=.false.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTO1' &
         , 'long_name', c='ozone destruction by O(3P)+O3')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ZDESTO1' &
         , 'units', c='molec/s')
    CALL channel_halt(substr, status)
    !
    ! mz_bk_20101222-

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
    IF (i_H2O_feedback > 0) THEN
       CALL mtend_register(my_handle, mtend_id_q)
    ENDIF
    IF (i_ice_feedback > 0) THEN
       CALL mtend_register(my_handle, mtend_id_xi)
    ENDIF
#endif

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE e4chem_init_memory

  !***************************************************************************


  SUBROUTINE e4chem_init_tracer     

    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, xt
    USE messy_main_timer,         ONLY: lstart
    USE messy_main_tracer,        ONLY: get_tracer

    IMPLICIT NONE

    INTRINSIC :: MAX, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'e4chem_init_tracer'
    INTEGER                     :: status   ! status flag
    INTEGER                     :: jt       ! tracer index

    CALL start_message_bi(modstr, 'E4CHEM tracer initialization', substr)

    ! check if H2O tracer exists (if yes, it is probably from H2O submodel)
    CALL get_tracer(status, GPTRSTR, 'H2O',idx=zidt_H2O)
    IF (status == 0) THEN
       CALL info_bi('Using H2O tracer for E4CHEM chemistry.', substr)
       ! op_pj_20160315+
       IF (i_H2O_feedback == 1) THEN
          CALL warning_bi('Feedback to specific humidity is disabled, if H2O tracer is present.', substr)
       ENDIF
       ! op_pj_20160315-
    ELSE
       CALL info_bi('No H2O tracer. Using qm1 and qte instead.', substr)
       zidt_H2O = 0
       ! op_pj_20160315+
       IF (i_H2O_feedback == 1) THEN
          CALL info_bi('Feedback to specific humidity is enabled.',substr)
       END IF
       ! op_pj_20160315-
    ENDIF

    ! op_pj_20160315+
    IF (i_ice_feedback == 1) THEN
       CALL info_bi('Feedback to ice is enabled.',substr)
    END IF
    ! op_pj_20160315-

    ! set negative values of e4chem tracers to zero
    IF (lstart) THEN ! don't do this after a restart !
       gp_tracer_loop: DO jt = 1, ntrac_gp
          IF (TRIM(ti_gp(jt)%tp%ident%submodel) == 'e4chem') THEN
             xt(_RI_XYZN_(:,:,:,jt)) = MAX(xt(_RI_XYZN_(:,:,:,jt)),0._dp)
          ENDIF
       ENDDO gp_tracer_loop
    ENDIF

    CALL end_message_bi(modstr, 'E4CHEM tracer initialization', substr)

  END SUBROUTINE e4chem_init_tracer

  !***************************************************************************

  ! ========================================================================
  SUBROUTINE e4chem_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    ! MESSy
    USE messy_main_channel, ONLY: get_channel_object, get_channel_info

    IMPLICIT NONE
    
    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'e4chem_init_coupling'
    INTEGER          :: status
    INTEGER          :: ip

    CALL start_message_bi(modstr, 'TEST COUPLING',substr) 
    CALL get_channel_info(status, TRIM(photrat_channel_gp))
    IF (status /= 0) THEN
       CALL info_bi('E4CHEM requires the channel '&
            &//TRIM(photrat_channel_gp))
       CALL warning_bi('Switch on the corresponding submodel',substr)
    ENDIF
    DO ip=1, IP_MAX
       CALL get_channel_object(status, TRIM(photrat_channel_gp), &
            'J_'//TRIM(jname(ip)), p3=photrat_gp(ip)%ptr)
       IF (status /= 0) &
            CALL warning_bi( 'J_'//TRIM(jname(ip))//&
            ' not in channel '//TRIM(photrat_channel_gp),substr)
    ENDDO

    ! GET POINTER TO TROPOPAUSE INDEX
    IF (p_parallel_io) THEN
       WRITE(*,*) '  INITIALIZING TROPOPAUSE LEVEL INDEX ...'
    END IF
    !
    CALL get_channel_object(status &
         , TRIM(tropop_idx%cha), TRIM(tropop_idx%obj), p2=tp_i)
    CALL channel_halt(substr, status)
    !
    IF (p_parallel_io) THEN
       WRITE(*,*) '  ... OK: ',&
            TRIM(tropop_idx%cha)//' - '// TRIM(tropop_idx%obj)
    END IF

! op_pj_20111201+
    ! GET POINTER TO COS ZENITH ANGLE
    IF (p_parallel_io) THEN
       WRITE(*,*) '  INITIALIZING COS ZENITH ANGLE ...'
    END IF
    !
    CALL get_channel_object(status, cossza%CHA, cossza%OBJ, p2=cossza_2d )
    CALL channel_halt(substr, status)
    !
    IF (p_parallel_io) THEN
       WRITE(*,*) '  ... OK: ',&
            TRIM(cossza%cha)//' - '// TRIM(cossza%obj)
    END IF

! op_pj_20111201-

    IF (l_fastscav) THEN
       CALL get_channel_object(status, 'cvtrans', 'kbot', p2=kconbot)
       IF (status == 1) &
            call error_bi(substr,&
            'channel object for bottom level of convection not found')
       CALL get_channel_object(status, &
            TRIM(lcover%cha), TRIM(lcover%obj), p3=pclcover)
       IF (status /= 0) CALL error_bi(substr,&
            ' channel object for ls cloud cover not found !')
       CALL get_channel_object(status, &
            TRIM(ratep%cha), TRIM(ratep%obj), p3=pmratep)
       IF (status /= 0) CALL error_bi(substr,&
            ' channel object for ls rain formation rate not found !')
       CALL get_channel_object(status, &
            TRIM(prec%cha), TRIM(prec%obj), p3=pfprec)
       IF (status /= 0) CALL error_bi(substr,&
            ' channel object for ls rain not found !')
       CALL get_channel_object(status, &
            TRIM(lwc%cha), TRIM(lwc%obj), p3=pmlwc)
       IF (status /= 0) CALL error_bi(substr,&
            ' channel object for cloud water content not found !')
       CALL get_channel_object(status, &
            TRIM(cvprec%cha), TRIM(cvprec%obj), p3=pcvdprec)
       IF (status /= 0) CALL error_bi(substr,&
            ' channel object for freshly formed convective precipitation not found !')
    ENDIF

    CALL end_message_bi(modstr, 'TEST COUPLING',substr) 

  END SUBROUTINE e4chem_init_coupling
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE e4chem_local_start

    USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma
    USE messy_main_grid_def_bi,     ONLY: philat_2d
    USE messy_main_data_bi,         ONLY: pmid => press_3d 

    IMPLICIT none
    INTRINSIC :: ABS, REAL

    INTEGER :: jk, jp

    flt_pscreg(_RI_XYZ__(1:kproma,jrow,:))= 0._dp
    
    ! SIMPLE GUESS FOR PSC REGION (NEEDED BY CLOUD)
    DO jk=1,nlev
       DO jp=1,kproma
          IF((REAL(jk,dp)<tp_i(jp,jrow))              &  ! above tropopause
               .and.(abs(philat_2d(jp,jrow))>=40._dp) &  ! abs(latitude) >= 40
               .and.(pmid(_RI_XYZ__(jp,jrow,jk)) >= 2800._dp))   &  ! 2800 Pa as zriceup in E4CHEM core
               flt_pscreg(_RI_XYZ__(jp,jrow,jk)) = 1._dp
       ENDDO
    ENDDO
 
  END SUBROUTINE e4chem_local_start
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE e4chem_physc

    ! This subroutine is the interface between ECHAM5 and CHEM inside the
    ! time loop. Before "CALL CHEMICS", several variables (which
    ! are declared in messy_e4chem_mem.f90) must be updated:
    !
    ! Conc(NSPEC) = concentrations of NSPEC chemical species [mcl/cm3]
    ! temp        = temperature [K]
    ! press       = pressure [Pa]
    ! cair        = c(air) [mcl/cm^3]
    ! JX(ip_*)    = several (1st order) photolysis rate coefficients [1/s]

    ! ECHAM5/MESSy
    USE messy_main_data_bi,   ONLY: &
           pint => pressi_3d  & ! level interface (above) pressure [Pa]
         , pmid => press_3d   & ! mid-level pressures [Pa]
#ifndef MESSYTENDENCY
         , PTM1 => tm1        &
         , PTTE => tte_3d     &
         , PQM1 => qm1        &
         , PQTE => qte_3d     &
         , xite => xite_3d    &
#endif
         , xim1 => xim1_3d
    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    USE messy_main_grid_def_bi,     ONLY: philat_2d
    USE messy_main_timer,           ONLY: MONTH, time_step_len
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte=>qxtte, pxtm1=>qxtm1
#endif
    USE messy_main_constants_mem, ONLY: pi, N_A, R_gas, M_air, M_H2O, vtmpc1

    IMPLICIT NONE

    INTRINSIC :: ACOS, MAX, SIZE, TINY, ASSOCIATED

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'e4chem_physc'
    INTEGER                     :: status   ! error status

    REAL(DP), PARAMETER :: scvmr = M_air/M_H2O ! op_pj_20160315

    ! tracer mixing ratio before chemistry
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zmr
    ! zmr after tracer scaling
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zmrbc
    ! tracer mixing ratio after chemistry
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zmrac

    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: temp
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: press
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: cair

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Conc
    ! specific humidity = m(H2O)/m(air) [kg/kg]
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: sphum

    REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: density

    REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: ZCLCOVER &
         , ZMRATEP, ZFEVAP, ZDP
    ! op_pj_20160315+
    REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: zwvac, zwvbc
    ! op_pj_20160315-
#ifdef MESSYTENDENCY
    REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: ztend
#endif

    ! LOCAL FIELDS
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE  :: JX
    INTEGER :: jk, jp, jt, ip, JL
    INTEGER :: idt ! um_ak_20110722 required for rank identifiers

    REAL(dp), PARAMETER :: U0LIM = pi/2.+0.052335956_dp ! 3deg past sunset

    ALLOCATE(press(kproma,nlev))
    ALLOCATE(cair(kproma,nlev))
    ALLOCATE(temp(kproma,nlev))
    ALLOCATE(sphum(kproma,nlev))
    ALLOCATE(Conc(kproma,nlev,NSPEC))
    ALLOCATE(zmr(kproma,nlev,ntrac_gp))
    ALLOCATE(zmrbc(kproma,nlev,ntrac_gp))
    ALLOCATE(zmrac(kproma,nlev,ntrac_gp))
    ALLOCATE(JX(kproma,nlev,IP_MAX))
    ALLOCATE(density(kproma,nlev))
    ALLOCATE(ZDP(kproma,nlev))
    ALLOCATE(ZCLCOVER(kproma,nlev), ZMRATEP(kproma,nlev))
    ALLOCATE(ZFEVAP(kproma,nlev))
    ! op_pj_20160315+
    ALLOCATE(zwvbc(kproma, nlev))
    ALLOCATE(zwvac(kproma, nlev))
    ! op_pj_20160315-
#ifdef MESSYTENDENCY
    ALLOCATE(ztend(kproma,nlev))
#endif

    ! CALCULATE ICE (VMR) FROM BASEMODEL ICE xi (kg/kg)
#ifndef MESSYTENDENCY
    ice(_RI_XYZ__(1:kproma,jrow,:)) =                     & ! op_pj_20121204
    MAX(( xim1(_RI_XYZ__(1:kproma,jrow,:))                &
         +xite(_RI_XYZ__(1:kproma,jrow,:))*time_step_len)       &
        * 28.84_dp/18.0_dp, tiny(1.0_dp))*flt_pscreg(_RI_XYZ__(1:kproma,jrow,:))
#else
    CALL mtend_get_start_l(mtend_id_xi, v0=ice(_RI_XYZ__(1:kproma,jrow,:)))
    ice(_RI_XYZ__(1:kproma,jrow,:)) = MAX(ice(_RI_XYZ__(1:kproma,jrow,:))*28.84_dp/18.0_dp &
         ,tiny(1.0_dp))*flt_pscreg(_RI_XYZ__(1:kproma,jrow,:))
#endif

#ifdef MESSYTENDENCY
    CALL mtend_get_start_l(mtend_id_t, v0=temp)
    CALL mtend_get_start_l(mtend_id_q, v0=sphum)
    sphum(:,:) = MAX(sphum(:,:), 0.0_dp)
#endif
    DO jk=1,nlev
       DO jp=1,kproma
          press(jp,jk) = pmid(_RI_XYZ__(jp,jrow,jk)) ! pressure
#ifndef MESSYTENDENCY
          temp(jp,jk)  = ptm1(_RI_XYZ__(jp,jrow,jk))  &
               + ptte(_RI_XYZ__(jp,jrow,jk)) * time_step_len
          sphum(jp,jk) = MAX(pqm1(_RI_XYZ__(jp,jrow,jk)) +&
               pqte(_RI_XYZ__(jp,jrow,jk)) * &
               time_step_len, 0._dp)   ! specific humidity
#endif
          cair(jp,jk)  = (N_A/1.E6_dp) * press(jp,jk) &
               / (R_gas*temp(jp,jk)*(1.0_dp+vtmpc1*sphum(jp,jk)))
       ENDDO
    ENDDO

    ! define photolysis rate constants for CHEMICS
    ! If photrat_gp(ip)%ptr is not associated, because
    !  - the tracer does not exist, or
    !  - the photolysis submodel does not provide it
    ! then jx(jp,jk,:) is 0.
    JX(:,:,:) = 0.0_dp
    DO ip=1, IP_MAX
       IF (ASSOCIATED(photrat_gp(ip)%ptr)) THEN 
          DO jk=1,nlev
             DO jp=1,kproma
                JX(jp,jk,ip) = photrat_gp(ip)%ptr(_RI_XYZ__(jp,jrow,jk))
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    DO jt=1, ntrac_gp
#ifdef MESSYTENDENCY
!       CALL mtend_get_start_l(mtend_id_tracer, v0=zmr(1:kproma,:,jt), idt=jt)
       CALL mtend_get_start_l(jt, v0=zmr(1:kproma,:,jt))
#else
      DO jk=1,nlev
          DO jp=1,kproma
             ! estimate tracer mixing ratio (mr) before chemistry (bc)
             zmr(jp,jk,jt) = pxtm1(_RI_X_ZN_(jp,jk,jt)) + &
                  pxtte(_RI_X_ZN_(jp,jk,jt)) *time_step_len
          ENDDO
       ENDDO
#endif
    ENDDO

    zmrbc(:,:,:) = zmr(:,:,:)
    zmrac(:,:,:) = zmr(:,:,:)
    Conc(:,:,:) = 0.0_dp

    CALL mr2c(cair)  ! convert mr to conc

    IF (l_fastscav) THEN
       ZFEVAP(:,:)=0._dp
       ZCLCOVER(:,:)=PCLCOVER(_RI_XYZ__(:,jrow,:))
       ZMRATEP(:,:)=PMRATEP(_RI_XYZ__(:,jrow,:))
       ! CALCULATE DENSITY, convert from g/m3 to g/cm3
       density(1:kproma,:)=pmid(_RI_XYZ__(1:kproma,jrow,:))*M_air/(temp(1:kproma,:)*R_gas)*1e-6_dp
       DO JK = 1, nlev 
          DO JL = 1, kproma 
             ZDP (JL, JK) = pint(_RI_XYZ__(JL,jrow,jk+1)) - pint (_RI_XYZ__(JL,jrow,jk)) 
             IF (ZCLCOVER (JL, JK) .LT.1E-5_dp.OR.PMLWC (_RI_XYZ__(JL,jrow,jk)) &
                  .LT.1E-10_dp) THEN
                ZCLCOVER (JL, JK) = 0._dp
                ZMRATEP (JL, JK) = 0._dp
             ENDIF
          ENDDO
       ENDDO
       CALL CLSCAV (Conc(:,:,:), time_step_len, density(1:kproma,:) &
            , ZDP(1:kproma,:), zmratep(1:kproma,:), pfprec(_RI_XYZ__(1:kproma,jrow,:))  &
            , zfevap(1:kproma,:), zclcover(1:kproma,:) &
            , pmlwc(_RI_XYZ__(1:kproma,jrow,:))    &
            , pcvdprec(_RI_XYZ__(1:kproma,jrow,:)), kconbot(1:kproma,jrow) &
            , temp(1:kproma,:) &
            , kproma, nlev)               
    ENDIF

    ! DANI, DAY NIGHT HOUSEKEEPING FOR CHEMISTRY                    
    ! U0LIM: correction for flat geometry zenith angle (3 deg allowance)
    DO JL = 1, kproma 
       IF (ACOS(cossza_2d(JL,jrow)) <= U0LIM) THEN 
          DANI(JL,jrow) = 1.1_dp
       ELSE 
          DANI(JL,jrow) = -1.1_dp
          IF (DANIM(JL,jrow).GT.1._dp) DANI(JL,jrow) = 0.1_dp 
       ENDIF
    ENDDO
    ! For Debugging+
!!$  DANI = 1.1_dp
!!$  DANIM = 1.1_dp
    ! For Debugging-

    CALL CHEMICS(kproma,    nlev,              jrow,                                &
         temp(1:kproma,:),         sphum(1:kproma,:),        Conc(:,:,:),           & 
         pmid(_RI_XYZ__(1:kproma,jrow,:)),    pint(_RI_XYZ__(1:kproma,jrow,:)),                         &
         DANI(1:kproma,jrow),      DANIM(1:kproma,jrow),                            &
         JX(:,:,ip_O3P),    JX(:,:,ip_O1D),    JX(:,:,ip_NO2),   JX(:,:,ip_HNO3),   &
         JX(:,:,ip_COH2),   JX(:,:,ip_CHOH),   JX(:,:,ip_N2O5),  JX(:,:,ip_HNO4),   &
         JX(:,:,ip_NO2O),   JX(:,:,ip_NOO2),   JX(:,:,ip_H2O2),  JX(:,:,ip_CH3OOH), &
         JX(:,:,ip_O2),     JX(:,:,ip_CFCl3),  JX(:,:,ip_CF2Cl2),JX(:,:,ip_N2O),    &
         JX(:,:,ip_CLONO2), JX(:,:,ip_CL2O2),  JX(:,:,ip_HOCL),  JX(:,:,ip_CCL4),   &
         JX(:,:,ip_CH3CL),  JX(:,:,ip_CH3CCL3),JX(:,:,ip_HCL),   JX(:,:,ip_H2O),    &
         JX(:,:,ip_NO),     JX(:,:,ip_CO2),                                         &  
         philat_2d(1:kproma,jrow), MONTH, time_step_len, tp_i(1:kproma,jrow)  & 
         , ZDELTAO3_BRV(_RI_XYZ__(1:kproma,jrow,:))                                       &
         , ZPRODO2(_RI_XYZ__(1:kproma,jrow,:)),   ZPRODCO(_RI_XYZ__(1:kproma,jrow,:))               &
         , ZPRODCH4(_RI_XYZ__(1:kproma,jrow,:)),  ZDESTH12(_RI_XYZ__(1:kproma,jrow,:))              &
         , ZDESTH14(_RI_XYZ__(1:kproma,jrow,:)),  ZDESTN13(_RI_XYZ__(1:kproma,jrow,:))              &
         , ZDESTC1(_RI_XYZ__(1:kproma,jrow,:)),   ZDESTCL2O2(_RI_XYZ__(1:kproma,jrow,:))            &
         , ZDESTCLOH(_RI_XYZ__(1:kproma,jrow,:)), ZDESTH8(_RI_XYZ__(1:kproma,jrow,:))               &
         , ZDESTH4(_RI_XYZ__(1:kproma,jrow,:)),   ZDESTH11(_RI_XYZ__(1:kproma,jrow,:))              &
         , ZDESTO1(_RI_XYZ__(1:kproma,jrow,:)))

    CALL c2mr(cair) ! convert conc to mr !

    DO idt=1, ntrac_gp
       DO jk=1,nlev
          DO jp=1,kproma
#ifndef MESSYTENDENCY
             pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt)) &
                  + (zmrac(jp,jk,idt) - zmr(jp,jk,idt)) / time_step_len
#else
             ztend(jp,jk) = (zmrac(jp,jk,idt) - zmr(jp,jk,idt)) / time_step_len
#endif
          ENDDO
       ENDDO
#ifdef MESSYTENDENCY
!       CALL mtend_add_l(my_handle, mtend_id_tracer, px=ztend, idt=idt)
       CALL mtend_add_l(my_handle, idt, px=ztend)
#endif
    ENDDO

    ! op_pj_20160315+
    ! NOTE: ALLOW FEEDBACK TO SPECIFIC HUMIDITY ONLY, IF NO H2O TRACER
    !       IS PRESENT. This limitation has been chosen to be downward
    !       compatible with previous versions/setups, in which the H2O
    !       submodel controls the feedback between the H2O tracer
    !       and the specific humidity. The limitation here is to
    !       avoid a double accounting of the feedback.
    !
    IF ( (zidt_H2O == 0) .AND. (i_H2O_feedback == 1) ) THEN
       !
       !
       !        dq    d      H2O           scvmr * d(H2O)/dt
       ! qte = ---- = -- (-----------) = ---------------------
       !        dt    dt  scvmr + H2O      (scvmr + H2O)^2
       !
       !                 q              H2O
       ! H2O = scvmr * -----  => q = -----------
       !               1 - q         scvmr + H2O
       !
#ifndef MESSYTENDENCY
       pqte(_RI_XYZ__(1:kproma,jrow,:)) = &
            pqte(_RI_XYZ__(1:kproma,jrow,:)) + &
            ( zwvac(:,:)/(scvmr+zwvac(:,:)) - zwvbc(:,:)/(scvmr+zwvbc(:,:)) ) &
            / time_step_len
#else
       ztend(:,:) = ( zwvac(:,:)/(scvmr+zwvac(:,:)) - &
            zwvbc(:,:)/(scvmr+zwvbc(:,:)) ) &
            / time_step_len
       CALL mtend_add_l(my_handle, mtend_id_q, px=ztend)
#endif
       !
    END IF
    ! op_pj_20160315-

    IF (I_ice_feedback == 1) THEN ! op_pj_20160315
#ifndef MESSYTENDENCY
    DO jk=1,nlev
       DO jp=1,kproma
          ! multiply with PSC region flag because E4CHEM ice is sensible only 
          ! in this region 
          xite(_RI_XYZ__(jp,jrow,jk))= xite(_RI_XYZ__(jp,jrow,jk))            &
               + (ice(_RI_XYZ__(jp,jrow,jk))*18.0_dp/28.84_dp - xim1(_RI_XYZ__(jp,jrow,jk)))    &
                 /time_step_len*flt_pscreg(_RI_XYZ__(jp,jrow,jk))
       ENDDO
    ENDDO
#else
    DO jk=1,nlev
       DO jp=1,kproma
          ztend(jp,jk) = (ice(_RI_XYZ__(jp,jrow,jk))*18.0_dp/28.84_dp - &
               xim1(_RI_XYZ__(jp,jrow,jk)))    &
               /time_step_len*flt_pscreg(_RI_XYZ__(jp,jrow,jk))
       END DO
    END DO
    CALL mtend_add_l(my_handle, mtend_id_xi, px=ztend)
#endif
    END IF ! op_pj_20160315

    DANIM(:,jrow) = DANI(:,jrow)

    ! DEALLOCATE MEMORY
    DEALLOCATE(press,cair,temp)
    DEALLOCATE(sphum)
    DEALLOCATE(Conc)
    DEALLOCATE(zmr)
    DEALLOCATE(zmrbc)
    DEALLOCATE(zmrac)
    DEALLOCATE(JX)
    DEALLOCATE(density)
    DEALLOCATE(ZDP)
    DEALLOCATE(ZCLCOVER, ZMRATEP)
    DEALLOCATE(ZFEVAP)
    ! op_pj_20160315+
    DEALLOCATE(zwvbc)
    DEALLOCATE(zwvac)
    ! op_pj_20160315-
#ifdef MESSYTENDENCY
    DEALLOCATE(ztend)
#endif

  CONTAINS
    !-------------------------------------------------------------------------

    SUBROUTINE mr2c(c_air) 
      ! convert mixing ratio [mol/mol] to concentration [mcl/cc]

      ! ind_* from messy_e4chem (global in this module)

      REAL(dp), INTENT(IN) :: c_air(:,:)

      ! special case: H2O
      IF (zidt_H2O/=0) THEN
         ! The tracer H2O is defined. It is assumed here that it contains
         ! the correct molar mixing ratio for H2O.
         Conc(:,:,ind_H2O) = c_air(:,:) * zmrbc(:,:,zidt_H2O)
         zwvbc(:,:) = zmrbc(:,:,zidt_H2O) ! op_pj_20160315
      ELSE
         ! op_pj_20160315+
         zwvbc(:,:) = M_air/M_H2O * sphum(:,:) / (1._dp-sphum(:,:))
         ! op_pj_20160315-
         ! The tracer H2O is not defined. Obtain H2O from specific humidity.
         ! op_pj_20160315+
!!$      Conc(:,:,ind_H2O) = c_air(:,:) * M_air/M_H2O * sphum(:,:) / (1._dp-sphum(:,:))
         Conc(:,:,ind_H2O) = c_air(:,:) * zwvbc(:,:)
         ! op_pj_20160315-
      ENDIF

      Conc(:,:,ind_H)       = c_air(:,:) * zmrbc(:,:,idt_H)
      Conc(:,:,ind_OH)      = c_air(:,:) * zmrbc(:,:,idt_OH)
      Conc(:,:,ind_HO2)     = c_air(:,:) * zmrbc(:,:,idt_HO2)
      Conc(:,:,ind_N)       = c_air(:,:) * zmrbc(:,:,idt_N)
      Conc(:,:,ind_NO)      = c_air(:,:) * zmrbc(:,:,idt_NO)
      Conc(:,:,ind_NO2)     = c_air(:,:) * zmrbc(:,:,idt_NO2)
      Conc(:,:,ind_NO3)     = c_air(:,:) * zmrbc(:,:,idt_NO3)
      Conc(:,:,ind_N2O5)    = c_air(:,:) * zmrbc(:,:,idt_N2O5)
      Conc(:,:,ind_HNO4)    = c_air(:,:) * zmrbc(:,:,idt_HNO4)
      Conc(:,:,ind_CL)      = c_air(:,:) * zmrbc(:,:,idt_CL)
      Conc(:,:,ind_CLO)     = c_air(:,:) * zmrbc(:,:,idt_CLO)
      Conc(:,:,ind_HOCl)    = c_air(:,:) * zmrbc(:,:,idt_HOCl)
      Conc(:,:,ind_CL2O2)   = c_air(:,:) * zmrbc(:,:,idt_CL2O2)
      Conc(:,:,ind_CL2)     = c_air(:,:) * zmrbc(:,:,idt_CL2)
      Conc(:,:,ind_HCHO)    = c_air(:,:) * zmrbc(:,:,idt_HCHO)
      Conc(:,:,ind_CH3O2)   = c_air(:,:) * zmrbc(:,:,idt_CH3O2)
      Conc(:,:,ind_CH4)     = c_air(:,:) * zmrbc(:,:,idt_CH4)
      Conc(:,:,ind_N2O)     = c_air(:,:) * zmrbc(:,:,idt_N2O)
      Conc(:,:,ind_H2O2)    = c_air(:,:) * zmrbc(:,:,idt_H2O2)
      Conc(:,:,ind_HCl)     = c_air(:,:) * zmrbc(:,:,idt_HCl)
      Conc(:,:,ind_CO)      = c_air(:,:) * zmrbc(:,:,idt_CO)
      Conc(:,:,ind_CH3OOH)  = c_air(:,:) * zmrbc(:,:,idt_CH3OOH)
      Conc(:,:,ind_ClNO3)   = c_air(:,:) * zmrbc(:,:,idt_ClNO3)
      Conc(:,:,ind_CFCl3)   = c_air(:,:) * zmrbc(:,:,idt_CFCl3)
      Conc(:,:,ind_CF2Cl2)  = c_air(:,:) * zmrbc(:,:,idt_CF2Cl2)
      Conc(:,:,ind_CH3CL)   = c_air(:,:) * zmrbc(:,:,idt_CH3CL)
      Conc(:,:,ind_CCL4)    = c_air(:,:) * zmrbc(:,:,idt_CCL4)
      Conc(:,:,ind_CH3CCL3) = c_air(:,:) * zmrbc(:,:,idt_CH3CCL3)
      Conc(:,:,ind_H2)      = c_air(:,:) * zmrbc(:,:,idt_H2)
      Conc(:,:,ind_HNO3)    = c_air(:,:) * zmrbc(:,:,idt_HNO3)
      Conc(:,:,ind_OHAB)    = c_air(:,:) * zmrbc(:,:,idt_OHAB)
      Conc(:,:,ind_HO2AB)   = c_air(:,:) * zmrbc(:,:,idt_HO2AB)
      Conc(:,:,ind_NAT)     = c_air(:,:) * zmrbc(:,:,idt_NAT)
      Conc(:,:,ind_O3P)     = c_air(:,:) * zmrbc(:,:,idt_O3P)
      Conc(:,:,ind_O3)      = c_air(:,:) * zmrbc(:,:,idt_O3)
      Conc(:,:,ind_O1D)     = c_air(:,:) * zmrbc(:,:,idt_O1D)
      Conc(:,:,ind_CO2)     = c_air(:,:) * zmrbc(:,:,idt_CO2)

! op_pj_20121204+
!!$      Conc(:,:,ind_ICE)     = c_air(:,:) * ice(1:kproma,:,jrow)
      Conc(:,:,ind_ICE)     = c_air(:,:) * ice(_RI_XYZ__(1:kproma,jrow,:))
! op_pj_20121204-

    END SUBROUTINE mr2c

    !-------------------------------------------------------------------------

    SUBROUTINE c2mr(c_air) 
      ! convert concentration [mcl/cc] to mixing ratio [mol/mol]

      ! ind_* from messy_e4chem (global in this module)

      REAL(dp), INTENT(IN) :: c_air(:,:)
      REAL(dp), DIMENSION(kproma,nlev) :: riac ! 1/c(air)

      riac(:,:) = 1._dp/c_air(:,:)

      zmrac(:,:,idt_H)       = riac(:,:) * Conc(:,:,ind_H)
      zmrac(:,:,idt_OH)      = riac(:,:) * Conc(:,:,ind_OH)
      zmrac(:,:,idt_HO2)     = riac(:,:) * Conc(:,:,ind_HO2)
      zmrac(:,:,idt_N)       = riac(:,:) * Conc(:,:,ind_N)
      zmrac(:,:,idt_NO)      = riac(:,:) * Conc(:,:,ind_NO)
      zmrac(:,:,idt_NO2)     = riac(:,:) * Conc(:,:,ind_NO2)
      zmrac(:,:,idt_NO3)     = riac(:,:) * Conc(:,:,ind_NO3)
      zmrac(:,:,idt_N2O5)    = riac(:,:) * Conc(:,:,ind_N2O5)
      zmrac(:,:,idt_HNO4)    = riac(:,:) * Conc(:,:,ind_HNO4)
      zmrac(:,:,idt_CL)      = riac(:,:) * Conc(:,:,ind_CL)
      zmrac(:,:,idt_CLO)     = riac(:,:) * Conc(:,:,ind_CLO)
      zmrac(:,:,idt_HOCl)    = riac(:,:) * Conc(:,:,ind_HOCl)
      zmrac(:,:,idt_CL2O2)   = riac(:,:) * Conc(:,:,ind_CL2O2)
      zmrac(:,:,idt_CL2)     = riac(:,:) * Conc(:,:,ind_CL2)
      zmrac(:,:,idt_HCHO)    = riac(:,:) * Conc(:,:,ind_HCHO)
      zmrac(:,:,idt_CH3O2)   = riac(:,:) * Conc(:,:,ind_CH3O2)
      zmrac(:,:,idt_CH4)     = riac(:,:) * Conc(:,:,ind_CH4)
      zmrac(:,:,idt_N2O)     = riac(:,:) * Conc(:,:,ind_N2O)
      zmrac(:,:,idt_H2O2)    = riac(:,:) * Conc(:,:,ind_H2O2)
      zmrac(:,:,idt_HCl)     = riac(:,:) * Conc(:,:,ind_HCl)
      zmrac(:,:,idt_CO)      = riac(:,:) * Conc(:,:,ind_CO)
      zmrac(:,:,idt_CH3OOH)  = riac(:,:) * Conc(:,:,ind_CH3OOH)
      zmrac(:,:,idt_ClNO3)   = riac(:,:) * Conc(:,:,ind_ClNO3)
      zmrac(:,:,idt_CFCl3)   = riac(:,:) * Conc(:,:,ind_CFCl3)
      zmrac(:,:,idt_CF2Cl2)  = riac(:,:) * Conc(:,:,ind_CF2Cl2)
      zmrac(:,:,idt_CH3CL)   = riac(:,:) * Conc(:,:,ind_CH3CL)
      zmrac(:,:,idt_CCL4)    = riac(:,:) * Conc(:,:,ind_CCL4)
      zmrac(:,:,idt_CH3CCL3) = riac(:,:) * Conc(:,:,ind_CH3CCL3)
      zmrac(:,:,idt_H2)      = riac(:,:) * Conc(:,:,ind_H2)
      zmrac(:,:,idt_HNO3)    = riac(:,:) * Conc(:,:,ind_HNO3)
      zmrac(:,:,idt_NAT)     = riac(:,:) * Conc(:,:,ind_NAT)
      zmrac(:,:,idt_OHAB)    = riac(:,:) * Conc(:,:,ind_OHAB)
      zmrac(:,:,idt_HO2AB)   = riac(:,:) * Conc(:,:,ind_HO2AB)
      zmrac(:,:,idt_O3P)     = riac(:,:) * Conc(:,:,ind_O3P)
      zmrac(:,:,idt_O3)      = riac(:,:) * Conc(:,:,ind_O3)
      zmrac(:,:,idt_O1D)     = riac(:,:) * Conc(:,:,ind_O1D)
      zmrac(:,:,idt_CO2)     = riac(:,:) * Conc(:,:,ind_CO2)

      zwvac(:,:)             = riac(:,:) * Conc(:,:,ind_H2O) ! op_pj_20160315
      ! If the tracer H2O is defined, its chemical tendency must
      ! eventually be put into xtte (via zmrac)
      ! op_pj_20160315+
!!$   IF (zidt_H2O/=0) zmrac(:,:,zidt_H2O) = riac(:,:) * Conc(:,:,ind_H2O)
      IF (zidt_H2O/=0) zmrac(:,:,zidt_H2O) = zwvac(:,:)
      ! op_pj_20160315-

      ice(_RI_XYZ__(1:kproma,jrow,:))     = riac(:,:) * Conc(:,:,ind_ICE)

    END SUBROUTINE c2mr

    !-------------------------------------------------------------------------

  END SUBROUTINE e4chem_physc
  ! ========================================================================



  ! ========================================================================
  SUBROUTINE e4chem_free_memory

    IMPLICIT NONE

    DEALLOCATE(RCGAS, sulook)

  END SUBROUTINE e4chem_free_memory
  ! ========================================================================

  !***************************************************************************

  SUBROUTINE e4chem_read_nml_cpl(status, iou)

    ! E4CHEM MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    ! read namelist for 'coupling' to channel containing reaction rates
    ! Author: Astrid Kerkweg, MPICH, Sep 2004

    ! MESSy
    USE messy_main_tools,          ONLY: read_nml_open, read_nml_check &
                                       , read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'e4chem_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR
   

  END SUBROUTINE e4chem_read_nml_cpl

  ! ***********************************************************************
END MODULE messy_e4chem_si
! ***********************************************************************

