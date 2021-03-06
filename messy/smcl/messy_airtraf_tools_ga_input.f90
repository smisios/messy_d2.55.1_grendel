MODULE messy_airtraf_tools_ga_input

 USE messy_main_constants_mem, ONLY: DP     !HY20140421-3

  IMPLICIT NONE

  !fort.1
  INTEGER, PARAMETER          :: iast_fort_1    = 300 
  INTEGER, PARAMETER          :: iaint_fort_1   = 5
  INTEGER, PARAMETER          :: mcros_fort_1   = 0 
  INTEGER, PARAMETER          :: mmute_fort_1   = 3 
  INTEGER, PARAMETER          :: idummy_fort_1  = 1
  REAL(DP), PARAMETER         :: rave_fort_1    = 1.0_dp
  REAL(DP), PARAMETER         :: rstd_fort_1    = 0.7_dp
  REAL(DP), PARAMETER         :: rout_fort_1    = 0.3_dp
  REAL(DP), PARAMETER         :: rcros_fort_1   = 1.0_dp
  REAL(DP), PARAMETER         :: pcros_fort_1   = 0.2_dp
  REAL(DP), PARAMETER         :: rmute_fort_1   = 0.1_dp
  REAL(DP), PARAMETER         :: pmute_fort_1   = 5.0_dp

  !fort.42
  INTEGER, PARAMETER          :: iter_fort_42    = -1 
  INTEGER, PARAMETER          :: npop_fort_42    = 100
  INTEGER, PARAMETER          :: nisland_fort_42 =  1 
  INTEGER, PARAMETER          :: nobj_fort_42    = -2
  INTEGER, PARAMETER          :: ncon_fort_42    =  0
  INTEGER, PARAMETER          :: ndv_fort_42     = 11
  INTEGER, PARAMETER          :: idum_fort_42    =  1
  INTEGER, PARAMETER          :: iy_fort_42      =  1
  
  !armoga.set
  INTEGER, PARAMETER          :: ndetail_armogaset  = 0     !Default=0, HY20140401-1,2

  !input.d
! ..... General Information ...........................
  CHARACTER(len=4), PARAMETER :: testp_input_d    = 'OUTE'  !Default=OUTE , defined in original input.d 
  INTEGER, PARAMETER          :: imaxmin_input_d  = -1 
  INTEGER, PARAMETER          :: ngend_input_d    = 100     !defined in original input.d 
  INTEGER, PARAMETER          :: iman_input_d     = 0       !Default=0, defined in original input.d  
  INTEGER, PARAMETER          :: idman_input_d    = 1 
  INTEGER, PARAMETER          :: indveva_input_d  = 1       !Default=1, defined in original input.d 
  INTEGER, PARAMETER          :: indvcnr_input_d  = 0
  INTEGER, PARAMETER          :: idebug_input_d   = 0       !defined in original input.d 
  INTEGER, PARAMETER          :: isetdes_input_d  = 1       !Default=1, HY20140331-4, defined in original input.d 
  INTEGER, PARAMETER          :: isubar_input_d   = 0 
  INTEGER, PARAMETER          :: iadinit_input_d  = 0 

! ..... I/O Data .......................................
  INTEGER, PARAMETER          :: mtout_input_d   =  6 
  INTEGER, PARAMETER          :: mtparam_input_d =  1 
  INTEGER, PARAMETER          :: mtgop_input_d   = 31 
  INTEGER, PARAMETER          :: mtgen_input_d   = 42 
  INTEGER, PARAMETER          :: mteva_input_d   = 43 
  INTEGER, PARAMETER          :: ipout_input_d   =  0 
  INTEGER, PARAMETER          :: mtobj_input_d   = 44 
  INTEGER, PARAMETER          :: mtprt_input_d   = 45 
  INTEGER, PARAMETER          :: mtall_input_d   = 46 
  INTEGER, PARAMETER          :: mtdeb_input_d   = 11 
  INTEGER, PARAMETER          :: imtall_input_d  =  1 
  INTEGER, PARAMETER          :: imemory_input_d =  1 
  INTEGER, PARAMETER          :: irank1_input_d  =  0 
  INTEGER, PARAMETER          :: infopt_input_d  =  1       !defined in original input.d 
  INTEGER, PARAMETER          :: mtsys_input_d   = 41 

! ..... Genetic Operator ...............................
  INTEGER, PARAMETER          :: icnrgeo_input_d    = 0     !defined in original input.d 
  INTEGER, PARAMETER          :: icnrrnk_input_d    = 1     !defined in original input.d 
  INTEGER, PARAMETER          :: icnrvio_input_d    = 1     !defined in original input.d 
  INTEGER, PARAMETER          :: mfiteval_input_d   = 1     !defined in original input.d 
  INTEGER, PARAMETER          :: ishare_input_d     = 0     !defined in original input.d 
  INTEGER, PARAMETER          :: msdist_input_d     = 1 
  INTEGER, PARAMETER          :: iallsh_input_d     = 1 
  INTEGER, PARAMETER          :: irksp_input_d      = 0 
  INTEGER, PARAMETER          :: mprank_input_d     = 1 
  INTEGER, PARAMETER          :: ictrlelt_input_d   = 0 
  INTEGER, PARAMETER          :: mselection_input_d = 1     !defined in original input.d  
  INTEGER, PARAMETER          :: ibestn_input_d     = 1 
  INTEGER, PARAMETER          :: isig2n_input_d     = 0 
  INTEGER, PARAMETER          :: ialpnorm_input_d   = 2 
  INTEGER, PARAMETER          :: iarchiv_input_d    = 0 
  INTEGER, PARAMETER          :: icrobd_input_d     = 0 
  INTEGER, PARAMETER          :: ismpop_input_d     = 0     !NEW
  INTEGER, PARAMETER          :: iarbase_input_d    = 0     !NEW 
  INTEGER, PARAMETER          :: mrpair_input_d     = 2     !NEW 
  INTEGER, PARAMETER          :: loopmax_input_d    = 1000  !NEW 
  INTEGER, PARAMETER          :: mfitsel_input_d    = 0     !NEW 
  INTEGER, PARAMETER          :: ibnarc_input_d     = 1 
  INTEGER, PARAMETER          :: idvsm_input_d      = 0     !NEW 

  REAL(DP), PARAMETER         :: shalpha_input_d    = 2.0_dp 
  REAL(DP), PARAMETER         :: rkcoef_input_d     = 0.4_dp 
  REAL(DP), PARAMETER         :: alpdmpar_input_d   = 0.01_dp 
  REAL(DP), PARAMETER         :: arcratio_input_d   = 0.1_dp
  REAL(DP), PARAMETER         :: arclimit_input_d   = 1.0_dp
  REAL(DP), PARAMETER         :: aratio_input_d     = 0.1_dp 
  REAL(DP), PARAMETER         :: xratio_input_d     = 0.5_dp  !NEW 
  REAL(DP), PARAMETER         :: arpopratio_input_d = 0.5_dp  !NEW 
  REAL(DP), PARAMETER         :: ararcratio_input_d = 0.1_dp  !NEW 
  REAL(DP), PARAMETER         :: rnpair_input_d     = 0.3_dp  !NEW 
  REAL(DP), PARAMETER         :: clr_input_d        = 0.08_dp !NEW
  REAL(DP), PARAMETER         :: gvaltol_input_d    = 0.05_dp 

END MODULE messy_airtraf_tools_ga_input 

