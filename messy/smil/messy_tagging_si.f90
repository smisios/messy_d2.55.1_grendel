! **********************************************************************
! INTERFACE TO ECHAM-5 FOR DIAGNOSTIC CHEMICAL
! TRACERS TO ESTIMATE THE IMPACT OF INDIVIDUAL EMISISONS, 
!     Details see messy_tagging.f90 
! Version: see 'modver' in !
!
! Authors: Volker Grewe,    DLR,   August 2008
!          Patrick Joeckel, MPICH, March 2006
!          Eleni Tsati, DLR, August 2013
!          Mariano Mertens, DLR, 2013-2017
!          Vanessa Rieger, DLR, 2017  
!
!       Changelog op_mm_2017
!       - 1. Tagging is now completly namelist controlled. Number of catagories and abbreviations are set in namelist 
!
!
!
!
!         Changelog op_vr_2017
!         - 1. correction of HOx tagging: budget is closed, new reactions for
!              stratosphere are added
!         - 2. new channel output:
!              fraction of odd oxygen per category to total odd oxygen: fracODD = ODDi / ODD
!              CH4 loss per tagging category: ch4loss_tag = frac_OH * ch4loss
!
!
!
!         Changelog op_mm_2015
!         - 1. changed code for MESSy 1 -> MESSy 2 (new tracer, "SMIL"ification)
!         - 2. "smil"ificated code (for use in COSMO)
!         - 3. Fixed various error (divisions by zero etc.) and improved computational performance
!         - 4. changed logic for restarts
!         - 5. tagging tracer now initialised at first timestep
!         - 6. using climatological calculated TP-level instead of fixed tp-level for tracer initialisation
!         - 7. revissed initialisation of tracers; now check if they have been initialised (from rerun file or tracer.nml ; see init_tracer_subroutine)
!         - 8. correction of HOx tagging: budget is closed, new reactions for
!              stratosphere are added
!         - 9. new channel output:
!              fraction of odd oxygen per category to total odd oxygen: fracODD = ODDi / ODD
!              CH4 loss per tagging category: ch4loss_tag = frac_OH * ch4loss
!
!
!
! References:
!     Grewe, 2004, ACP 4.
!     Grewe et al., GMD, 2010 and Grewe, GMD, 2013.
!     Rieger et al., GMD 2018
! 
!
!
!
! References:
!     Grewe, 2004, ACP 4.
!     Grewe et al., GMD, 2010 and Grewe, GMD, 2013.
!
! TO DO
! - New tracer:  use new_tracer instead of new_tracer_old
! - Define prod and loss terms (go through reactions.)
! - Eliminate _i in init_tracer
! - Include C production from CFC loss goes into X_ch4 ! not necessary since not 
! - Some coupling might be better done via namelist: existent
!   Basic idea:
!        Iteration of each tracer follows the idea to distribute loss 
!                                                           and production terms 
!        of the tracer to individual categories. I.e.
!        From Chemistry one gets       d/dt X  = P  - DX
!        This has to be converted into d/dt Xi = Pi -DiXi
!
!        Note that Prod and Loss should be strictly non-negative !!
!                                      
!   Detailed approach:
!        Two types of parameterisations have to be distinguished:
!            (1) Gas phase chemistry
!            (2) all others (wash-out, deposition, etc.) 
!        Add (2)
!        The latter one is treated as a bulk and treated as a process, which acts 
!        uniformly on the individual species. This is correct for dry and wet 
!        deposition, and generally all loss terms, which this approach aims at most.
!        However, it introduces a diffusion in the case of production terms, 
!        i.e. the re-evaporation of HNO3 at lower levels is locally a production term,
!
!        which should have a decomposition into the categories,
!        with a distribution according to the upper level. 
!        Since the current lower level is only known, this
!        is not feasable without large efforts. This is  thought  to be a minor 
!        contributor to the  overall uncertainties and hence ignored.
!        
!        So the tendency of all (2) processes is 
!        
!        | X=Specie from MECCA    Yi=tagged specie SYi= Sum over all tagged 
!        |- Transport !        |- messy_local_start         X0=XTM1+dt*XTE  DX(transprt)=dt*XTE
!        |                            Yi=XTM1=dt*XTE  Y0=Summe Yi
!        |                            EX-trans=X0-Y0  Error due to numerical diffusion
!        |                                            by transport algorithm
!        |                            Yi=Yi*X0/Y0     Correction due to advection 
!        |                                                                   diffusion
!        |- vdiff                     Emissions via offlem
!        |                            Biogenic emissions via onlem. Assumption:
!        |                            All emissions are biogenic = soils
!        |                            Except for other emissions: Lightning in physc
!        |- physc
!           |- processes
!           |- tagging(1) pre-MECCA   update all missing emissions, e.g. lightning
!           |                         X1=XTM1+dt*XTE  DX(processe) = (X1-X0)  
!           |                         Yi=XTM1+dt*XTE  Y1=Summe Yi
!           |                         X1/Y1           Linear factor representing all 
!           |                                                            processes !       
!           |                                         Yi=Yi*X1/Y1  Application of all other 
!           |                                         processes to tagged species
!           |- MECCA
!           |- tagging(2) post-MECCA  X2=XTM1+dt*XTE  DX(MECCA)    = (X2-X1)
!           |                         Distribution of Prod and Loss terms to 
!           |                                   individual categories i
!           |                         i.e. calculate Pi Di consistently with MECCA
!           |                         d/dt Yi = Pi + Di * Yi      
!           |                               Yi=(Yi+Pi*dt)/(1+Di*dt/Yi)  and  Y2+Sum Yi
!           |                         make sure that Pi and Di are calculated to 
!           |                                                            yield X2 = Y2
!           |                         Update tendencies
! **********************************************************************
#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_tagging_si
! **********************************************************************

  ! USE TOOLS
  USE messy_main_blather_bi,   ONLY: start_message_bi, end_message_bi

  USE messy_tagging              ! global data and public 

  USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev,ngpblks
  USE messy_main_tools,       ONLY:PTR_3D_ARRAY

!changed: time_step_len from timer and not main_data_bi
  USE messy_main_timer,       ONLY: time_step_len
  USE messy_main_constants_mem, ONLY: STRLEN_ULONG



IMPLICIT NONE
! PRIVATE

! Global parameters defining the tagging 
! Hard wired since it depends on the number of emissions tagged and the kind of 
!      chemical system u sed. 
!      It would be nice to have a check whether it is consistent.  To Do

!  INTEGER, PARAMETER         ::   N_tag_species=10
!                                        ! Number of categories for tagged species
!                                        ! In principle 8 Emissions + CH4 
                                        !              + O2 photolosis for O3
    INTEGER, PARAMETER         ::   Fam_NOy_species=15 , &   ! Number of N-species in Family 
                                                             ! Family without PAN
                                    Fam_NMHC_species=42      ! Number of NMHC-species 
                                                             !  in Family


  CHARACTER (LEN=10), DIMENSION(Fam_NOy_species), PARAMETER  ::  NOy_species=&
       (/'N         ','NO        ','NO2       ','NO3       ','N2O5      ', &
         'HONO      ','HNO3      ','HNO4      ','NACA      ','MPAN      ', &
         'IC3H7NO3  ','LC4H9NO3  ','ISON      ','ClNO3     ','BrNO3     '/)

! changed - needed to correct families 
! not ProNO2 Definition of moelcules via Namelist?

  REAL(DP), DIMENSION(Fam_NOy_species), PARAMETER  ::   N_in_Fam=  &
                                     (/1._dp,1._dp,1._dp,1._dp,2._dp,1._dp,1._dp,    &
                                       1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,      &
                                         1._dp,1._dp            /)
 
! Check S, Br, Cl, I compounds, how should he feed into NMHC fam.?
!

 CHARACTER (LEN=10), DIMENSION(Fam_NMHC_species), PARAMETER  :: & 
        NMHC_species=     &
        (/ 'CH3OH     ', 'CH3O2     ', 'CH3OOH    ', 'HCHO      ', &
           'HCOOH     ',                                           & ! 5    1C
           'C2H6      ', 'C2H4      ', 'C2H5O2    ', 'C2H5OOH   ', &
           'CH3CHO    ',                                           & ! 
           'CH3CO2H   ', 'CH3CO3    ', 'CH3CO3H   ', 'NACA      ', & ! 9    2C
           'C3H8      ', 'C3H6      ', 'IC3H7O2   ', 'IC3H7OOH  ', &
           'LHOC3H6O2 ',                                           &
           'LHOC3H6OOH', 'CH3COCH3  ', 'CH3COCH2O2', 'HYPERACET ', &
           'ACETOL    ',                                           &   
           'MGLYOX    ', 'MPAN      ', 'IC3H7NO3  ',               & ! 13   3C
           'NC4H10    ', 'LC4H9O2   ', 'LC4H9OOH  ', 'MVK       ', &
           'MVKO2     ',                                           &   
           'MVKOOH    ', 'MEK       ', 'LMEKO2    ', 'LMEKOOH   ', &
           'BIACET    ',                                           & 
           'LC4H9NO3  ',                                           & ! 11   4C
           'C5H8      ', 'ISO2      ', 'ISOOH     ', 'ISON      ' /) ! 4    5C


  REAL(DP), DIMENSION(Fam_NMHC_species), PARAMETER  ::   C_in_Fam =            &
                                     (/1._dp,1._dp,1._dp,1._dp,1._dp,          &
                                       2._dp,2._dp,2._dp,2._dp,2._dp,          &
                                       2._dp,2._dp,2._dp,2._dp,                &
                                       3._dp,3._dp,3._dp,3._dp,3._dp,          &
                                       3._dp,3._dp,3._dp,3._dp,3._dp,          &
                                       3._dp,3._dp,3._dp,                      &
                                       4._dp,4._dp,4._dp,4._dp,4._dp,          &
                                       4._dp,4._dp,4._dp,4._dp,4._dp,          &
                                       4._dp,                                  &
                                       5._dp,5._dp,5._dp,5._dp         /)

  REAL(DP), Parameter                :: C_in_PAN=2.0_dp, factor_small_init=1.e-4

  ! ... STREAM ELEMENTS: (DIAGNOSED FIELDS)
! Input to module 
! chemical species
! Initialization 
! REAL(DP), DIMENSION(:,:),   POINTER :: tp_p_wmo=>NULL()     ! kproma*ngl
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_clim=>NULL()     ! kproma*ngl


! Emissions
  REAL(DP), DIMENSION(:,:,:), POINTER :: telnox=>NULL()       ! kproma*nlev*ngpblks
                                                              ! mol/mol/s

  REAL(DP), DIMENSION(:,:,:), POINTER :: pxtte_pre_emis_NOy=>NULL() ! kproma*nlev*ngpblks
  REAL(DP), DIMENSION(:,:,:), POINTER :: pxtte_pre_emis_NMHC=>NULL()! mol/mol/s

                                          ! kproma*nlev t=tracer d=diagnosed
! Module variables 
! Chemical X0, X1 and X2 values   ! nproma*nlev*ngpblks
  REAL(DP), DIMENSION(:,:,:), POINTER ::   o3_0=>NULL(),   o3_1=>NULL(),  o3_2=>NULL() !, oh_0=>NULL(),  oh_1=>NULL() , oh_2=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER ::  noy_0=>NULL(),  noy_1=>NULL(), noy_2=>NULL() !, ho2_0=>NULL(), ho2_1=>NULL(), ho2_2=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER ::  pan_0=>NULL(),  pan_1=>NULL(), pan_2=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER ::   co_0=>NULL(),   co_1=>NULL(),  co_2=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: nmhc_0=>NULL(), nmhc_1=>NULL(),nmhc_2=>NULL()





! for these only X2 values are needed
  REAL(DP), DIMENSION(:,:,:), POINTER ::                  ch4_1=>NULL(), ch4_2=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER ::                                  oh_2=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER ::                                 ho2_2=>NULL()


! Chemical tagged species advected   nproma*nlev*n_tag_species*ngpblks
! REAL(DP), DIMENSION(:,:,:,:) , POINTER  :: o3_tagging=>NULL(), noy_tagging=>NULL(),   &
!                                            pan_tagging=>NULL(), &
!                                            co_tagging=>NULL(), nmhc_tagging=>NULL()



 REAL(DP), DIMENSION(:,:,:,:) , allocatable  :: o3_tagging, noy_tagging,   &
                                            pan_tagging, &
                                            co_tagging, nmhc_tagging







 ! REAL(DP), DIMENSION(:,:,:,:) , POINTER  :: oh_tagging=>NULL(), ho2_tagging=>NULL() 

!  REAL(DP), DIMENSION(:,:,:)   , POINTER  :: o3_tagging_sum=>NULL(),  &
!                                             noy_tagging_sum=>NULL(), &
!                                             pan_tagging_sum=>NULL(), &
!                                             co_tagging_sum=>NULL(),  &
!                                             nmhc_tagging_sum=>NULL()

  REAL(DP), DIMENSION(:,:,:) , allocatable  :: o3_tagging_sum,  &
                                             noy_tagging_sum, &
                                             pan_tagging_sum, &
                                             co_tagging_sum,  &
                                             nmhc_tagging_sum





 ! REAL(DP), DIMENSION(:,:,:)   , POINTER  :: oh_tagging_sum=>NULL(), ho2_tagging_sum=>NULL() 

  REAL(DP), DIMENSION(:,:,:)   , POINTER  :: do3_err=>NULL(),   &
                                             dnoy_err=>NULL(),  &
                                             dnmhc_err=>NULL(), &
                                             dpan_err=>NULL(),  &
                                             dco_err=>NULL()

   
  REAL(DP), DIMENSION(:,:,:), POINTER :: ch4=>NULL(),co=>NULL(),nmhc=>NULL(),noy=>NULL(),o3=>NULL(),  &
                                          oh=>NULL(),ho2=>NULL(),pan=>NULL() , O1D=>NULL()


! Type definition for Tagging
 TYPE  ::  tagging_specie_global_ptr !(for pointer e.g. channel)
           sequence
           REAL(DP), DIMENSION(:,:,:), POINTER       :: ptr 
           CHARACTER (len=10)                        :: name
           CHARACTER (len=30)                        :: longname
 END TYPE tagging_specie_global_ptr

 TYPE  ::  tagging_specie_global  !(for tracer)
           sequence
           INTEGER                                   :: idx
           CHARACTER (len=10)                        :: name
           CHARACTER (len=30)                        :: longname
 END TYPE tagging_specie_global

 
 !op_mm_16102017+
  TYPE T_CATEGORY_SET
     CHARACTER(LEN=STRLEN_ULONG)    :: abr  = '' ! abreviation
     CHARACTER(LEN=STRLEN_ULONG)    :: name  = '' ! long_name
     INTEGER                        :: ncat
     INTEGER, DIMENSION(:), Pointer  :: category  
  END TYPE T_CATEGORY_SET

  TYPE(T_CATEGORY_SET), DIMENSION(:), POINTER  ::  CATEGORY_SET
 !op_mm_16102017-
  
  

!tagging species
 !TYPE(tagging_specie_global), DIMENSION(N_tag_species)  :: &
 !                                   o3_tag,  &
 !                                   noy_tag, pan_tag, nmhc_tag, co_tag
 
 !op_mm_16102017
TYPE(tagging_specie_global), DIMENSION(:), POINTER  :: &
                                            o3_tag=>NULL(),  noy_tag=>NULL(), &
                                            pan_tag=>NULL(), nmhc_tag=>NULL(), co_tag=>NULL()
                                    
!op_mm_16102017
!TYPE(tagging_specie_global_ptr), DIMENSION(N_tag_species)  :: &
!                                   oh_tag, ho2_tag, ch4loss_tag   !op_vr_20170215 
!TYPE(tagging_specie_global_ptr), DIMENSION(N_tag_species)  ::  frac_ODD  !op_vr_20170215 

TYPE(tagging_specie_global_ptr), DIMENSION(:), POINTER   :: &
                                   oh_tag=>NULL(), ho2_tag=>NULL(), ch4loss_tag=>NULL()   !op_vr_20170215 
TYPE(tagging_specie_global_ptr), DIMENSION(:), POINTER  ::  frac_ODD=>NULL()  !op_vr_20170215 


!op_mm_16102017
! Tendencies of different species
!TYPE(tagging_specie_global_ptr), DIMENSION(:), POINTER  ::           &
!                                  p_ho2_ten, p_ro2_ten, d_oh_ten       & 
!                                 ,d_ho2_ten, d_no_ten, d_ro_ten        &
!                                 ,d_xo_ten
 TYPE(tagging_specie_global_ptr), DIMENSION(:), POINTER::           &
                                  p_ho2_ten, p_ro2_ten, d_oh_ten       & 
                                 ,d_ho2_ten, d_no_ten, d_ro_ten        &
                                 ,d_xo_ten


!------
!debug fields for tropop
!------
!REAL(DP), DIMENSION(:,:,:), POINTER :: debug_press=>NULL()
!REAL(DP), DIMENSION(:,:), POINTER   :: debug_zklim=>NULL()
!REAL(DP), DIMENSION(:,:,:), POINTER :: debug_rts1=>NULL()
!REAL(DP), DIMENSION(:,:,:), POINTER :: debug_rts2=>NULL()
!REAL(DP), DIMENSION(:,:,:), POINTER :: debug_rts3=>NULL()
!REAL(DP), DIMENSION(:,:,:), POINTER :: debug_rts4=>NULL()


!------
! Chemical tendencies for diagnostics
!------

  REAL(DP), DIMENSION(:,:,:), POINTER :: o3prod_ho2=>NULL()  ! kproma*nlev*ngl
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3prod_ro2=>NULL()  ! kproma*nlev*ngl
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3prod_o2=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3prod=>NULL()

  REAL(DP), DIMENSION(:,:,:), POINTER :: o3loss_ho2=>NULL()              
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3loss_ro=>NULL()            
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3loss_oh=>NULL()             
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3loss_no=>NULL()              
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3loss_xo=>NULL()              
  REAL(DP), DIMENSION(:,:,:), POINTER :: o3loss=>NULL()

  REAL(DP), DIMENSION(:,:,:), POINTER :: panprod=>NULL()            
  REAL(DP), DIMENSION(:,:,:), POINTER :: panloss=>NULL()              

  REAL(DP), DIMENSION(:,:,:), POINTER :: coprod=>NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: coloss=>NULL()

  REAL(DP), DIMENSION(:,:,:), POINTER :: ch4loss=>NULL()

  REAL(DP), DIMENSION(:,:,:), POINTER :: ProdtOH=>NULL() 
  REAL(DP), DIMENSION(:,:,:), POINTER :: LosstOH=>NULL() 
  REAL(DP), DIMENSION(:,:,:), POINTER :: ProdtHO2=>NULL() 
  REAL(DP), DIMENSION(:,:,:), POINTER :: LosstHO2=>NULL() 

  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2100=>NULL()           ! H + O2 -> HO2           
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2103=>NULL()           ! OH + O3P -> H + O2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2104=>NULL()           ! OH + O3 -> HO2 + O2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2105=>NULL()           ! OH + H2 -> H2O + H
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2106=>NULL()           ! OH + H2 -> H2O + H
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2107=>NULL()           ! HO2 + O3 -> OH + 2 O2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2109=>NULL()           ! HO2 + OH -> H2O + O2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2110=>NULL()           ! HO2 + HO2 -> H2O2 + O2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2111=>NULL()           ! H2O + O1D -> 2 OH
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG2112=>NULL()           ! H2O2 + OH -> H2O + HO2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG3200=>NULL()           ! NO + OH -> HONO
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG3201=>NULL()           ! NO + HO2 -> NO2 + OH
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG3202=>NULL()           ! NO2 + OH -> HNO3
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG3203=>NULL()           ! NO2 + HO2 -> HNO4
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG3207=>NULL()           ! HNO4 -> NO2 + HO2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG4101=>NULL()           ! CH4 + OH -> CH3O2 + H2O
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG4110=>NULL()           ! CO + OH -> H + CO2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG6203=>NULL()           ! ClO + OH -> .94 Cl + .94 HO2 + .06 HCl + .06 O2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG6204=>NULL()           ! ClO + HO2 -> HOCl + O2
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossG7201=>NULL()           ! BrO + HO2 -> HOBr + O2

  REAL(DP), DIMENSION(:,:,:), POINTER :: LossJ3200=>NULL()           ! HONO + hv -> NO + OH
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossJ3201=>NULL()           ! HNO3 + hv -> NO2 + OH
  REAL(DP), DIMENSION(:,:,:), POINTER :: LossJ4101b=>NULL()          ! HCHO + hv ->  H + CO + HO2

  REAL(DP), DIMENSION(:,:,:), POINTER :: OHlossNMHC=>NULL()          ! NMHC + OH -> NMHC
  REAL(DP), DIMENSION(:,:,:), POINTER :: HO2lossNMHC=>NULL()         ! NMHC + HO2 -> NMHC
  REAL(DP), DIMENSION(:,:,:), POINTER :: OHlossHO2prodNMHC=>NULL()   ! NMHC + OH -> NMHC + HO2
  REAL(DP), DIMENSION(:,:,:), POINTER :: HO2prodNMHCNOy=>NULL()      ! NMHC + NOy -> HO2 + NMHC + NOy
  REAL(DP), DIMENSION(:,:,:), POINTER :: HO2prodNMHCphoto=>NULL()    ! NMHC + hv -> NMHC + HO2

  REAL(DP), DIMENSION(:,:,:), POINTER :: nmhcloss=>NULL()

  REAL(DP), DIMENSION(:,:,:), POINTER :: emis_noy_onlem=>NULL(), &
                                         emis_nmhc_onlem=>NULL()


  !REAL(DP), DIMENSION(:,:,:), POINTER           :: O_tag_lig=>NULL(),O_tag_bio=>NULL()
  !REAL(DP), DIMENSION(:,:,:), POINTER           :: A_tag_lig=>NULL(),A_tag_bio=>NULL()
  !REAL(DP), DIMENSION(:,:,:), POINTER           :: B_tag_lig=>NULL(),B_tag_bio=>NULL()
  !REAL(DP), DIMENSION(:,:,:), POINTER           :: G_tag_lig=>NULL(),G_tag_bio=>NULL()

!------
! Chemical tendencies for tagged species
!------
  REAL(DP), DIMENSION(:,:,:,:),  POINTER ::  do3=>NULL(), dco=>NULL(),  &
                                             dnoy=>NULL(),dnmhc=>NULL(),&
                                             dpan=>NULL() 

!! op_mm_20131206
! unnecessary debug field
!  REAL(DP), DIMENSION(:,:,:),    POINTER ::  do3_out=>NULL(),  &
!                                             dco_out=>NULL(),  &
!                                             dnoy_out=>NULL(), &
!                                             dnmhc_out=>NULL(),&
!                                             dpan_out=>NULL()                         
! op_mm_20131202
!---------------------
!Fields for parsing data to SMCL
!--------------------

REAL(DP), DIMENSION(:,:,:), Pointer ::  oh_tag_data=>NULL(), ho2_tag_data=>NULL(),      &
                                        p_ho2_ten_data=>NULL(), p_ro2_ten_data=>NULL(), &
                                        d_ro2_ten_data=>NULL(), d_ho2_ten_data=>NULL(), &
                                        d_xo_ten_data=>NULL(), d_oh_ten_data=>NULL(),   &
                                        d_no_ten_data => NULL(), d_ro_ten_data=>NULL(), &
                                        ch4loss_tag_data=>NULL(), frac_ODD_data=>NULL()     !op_vr_20170215
                                 



  ! GLOBAL COUPLING SWITCHES
  CHARACTER(len=16),dimension(2) ::  c_lnox !,c_tropop ! op_mm_20140124 removed tropo: now calculated online
  INTEGER                        :: i_diag=1, iout=1   !   (0: kein diagn output, 1: diag outp)
                                                   
  INTEGER :: &
             itrac_NOy_errP,itrac_NMHC_errP,itrac_CO_errP,itrac_PAN_errP,itrac_O3_errP,&!itrac_ohe_errP,itrac_ho2e_errP, &
             itrac_NOy_errN,itrac_NMHC_errN,itrac_CO_errN,itrac_PAN_errN,itrac_O3_errN !itrac_ohe_errN,itrac_ho2e_errN

!Bc dOH=dHO2 + for t(0), 'pxtte_OH'= 'pxtte_HO2' => then the OH_tagging_sum = HO2_tagging_sum  for every t(x). so, we need only one of two.
! 
! PTR TO OZONE TRACER/STREAM ELEMENT
! POINTER to advected chemical species in tracer_gp stream
  INTEGER                               :: idx_o3, idx_CO, idx_PAN, idx_oh, idx_ho2, idx_CH4 , idx_O1D
  INTEGER, DIMENSION(Fam_NOy_species)   :: idx_NOy       ! Pointer Array for NOy  compounds
  INTEGER                               :: ip_HNO3, ip_HCHO     ! needed for init values
  INTEGER, DIMENSION(Fam_NMHC_species)  :: idx_NMHC      ! Pointer Array for NMHC compounds
  
! POINTER to non-advected diagnostics in tracer_gp stream
! ozone:
  INTEGER :: idx_o3prod_ho2, idx_o3prod_ro2, idx_o3prod_MeO2, idx_o3prod_o2
  INTEGER :: idx_o3loss_ho2, idx_o3loss_oh, idx_o3loss_ro,   & 
             idx_o3loss_no, idx_o3loss_xo


  INTEGER :: idx_noprod_N2O

! PAN
  INTEGER :: idx_panprod, idx_panloss

! HOx
  INTEGER :: idx_LossG2100, idx_LossG2103, idx_LossG2104, idx_LossG2105,     &
             idx_LossG2106, idx_LossG2107, idx_LossG2109, idx_LossG2110,     &
             idx_LossG2111, idx_LossG2112, idx_LossG3200, idx_LossG3201,     &
             idx_LossG3202, idx_LossG3203, idx_LossG3207, idx_LossG4101,     &
             idx_LossG4110, idx_LossG6203, idx_LossG6204, idx_LossG7201,     &
             idx_LossJ3200, idx_LossJ3201, idx_LossJ4101b,                   &
             idx_OHlossNMHC, idx_HO2lossNMHC, idx_OHlossHO2prodNMHC,         &
             idx_HO2prodNMHCNOy, idx_HO2prodNMHCphoto


! tagged tracers
!op_mm_16102017
!INTEGER, DIMENSION(N_tag_species)   :: itrac_NOy, itrac_NMHC, itrac_PAN, itrac_CO
!INTEGER, DIMENSION(N_tag_species)   :: itrac_O3  


INTEGER,POINTER, DIMENSION(:)   :: itrac_NOy=>NULL(), itrac_NMHC=>NULL(), itrac_PAN=>NULL() & 
                                  ,itrac_CO=>NULL(), itrac_O3=>NULL()


! Logicals controlling the initialisation
  LOGICAL                             :: linit_tagging_tracer, linit_tagging_tracer_now
  LOGICAL                             :: ltagging ! run tagging method    ! only after initialisation
  LOGICAL                             :: force_reinitialize_tracer ! switch to force initialisation of tagging tracer 

  !op_mm_16102017
  INTEGER, DIMENSION(:), POINTER :: category=>NULL()
  
!day-night_LossO1D_condition
REAL(DP), PARAMETER                             :: limit_o1D= 1.e-21_dp 
REAL(DP), PARAMETER                             :: limit_tagg_pos=1.e-21_dp
REAL(DP), PARAMETER                             :: limit_tagg_neg=-1.e-21_dp
REAL(DP), PARAMETER                             :: tag_limit= 10.e-18_dp
REAL(DP), PARAMETER                             :: pan_limit= 1.e-11_dp
REAL(DP), PARAMETER                             :: pan_limit_1= 1.e-09_dp
  ! -------------------------------------------------------------

  ! PUBLIC ECHAM-5 INTERFACE ROUTINES TO 'messy_SUBMODELS'
  PUBLIC :: tagging_global_start        ! Probably not necessary
  PUBLIC :: tagging_initialize          ! global initialisation of module
  PUBLIC :: tagging_new_tracer          ! define new tracers
  PUBLIC :: tagging_init_tracer         ! initialize tracer fields
 ! PUBLIC :: tagging_init_tracer_online  ! initialize tracer with mecca fields
  PUBLIC :: tagging_init_memory         ! allocate memory and define channels
  PUBLIC :: tagging_init_coupling       ! Initialize ...
  PUBLIC :: tagging_local_start         ! Define X0 values and do online  
                                        !  tracer initialisation in the case of a new start
  PUBLIC :: tagging_physc               ! integrate one time step
  PUBLIC :: tagging_vdiff               ! determine the soil tendencies from emissions
  PUBLIC :: tagging_free_memory         ! free memory

  ! PRIVATE ECHAM-5 INTERFACE ROUTINES
  PRIVATE :: tagging_read_nml_e5     ! initialize 'coupling' to E5/ATTILA
                                       ! ( /CPL/-namelist )
  PRIVATE :: tagging_init_tracer_online  ! initialize tracer with mecca fields
CONTAINS

! ************************************************************************
! PUBLIC ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ------------------------------------------------------------------------
  SUBROUTINE  tagging_initialize

    ! Tagging MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF TAGGING SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT : Not yet done
    !
    ! Author: Patrick Joeckel, MPICH, July 2002
    !         Volker Grewe, DLR, AUGUST 2008

    ! ECHAM5
    !changed finish -> messy_main_blateher_bi error_bi
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY:error_bi, warning_bi
    USE messy_main_tools,   ONLY: find_next_free_unit
    USE messy_main_timer,    ONLY: lstart, lresume

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_initialize'
    INTEGER(I4)         :: iou    ! I/O unit
    INTEGER(I4)         :: status ! error status
    INTEGER             :: ispecie,i,n,ncat,j
    INTEGER, PARAMETER  :: MAXPSTRLEN = 50

    CALL start_message_bi(modstr, 'TAGGING INITIALIZATION', substr)

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL tagging_read_nml(status, iou)
       IF (status /= 0) CALL error_bi('Error in Reading Namelist CTRL', substr)
    END IF


    CALL p_bcast(i_integrate, p_io)
   !CALL p_bcast(i_method, p_io)   ! ET
    CALL p_bcast(i_species  , p_io)
    CALL p_bcast(i_advect_scaling, p_io)
    CALL p_bcast(i_tracer_init, p_io)
    CALL p_bcast(l_adv_err_diag, p_io)

    IF (p_parallel_io) THEN
       n = 0
       DO i=1, N_MAX_CATEGORY
          IF (TRIM(CATEGORY_IN(i)%abr) == '') CYCLE
           IF (TRIM(CATEGORY_IN(i)%long_name) == '') CYCLE
          IF (TRIM(CATEGORY_IN(i)%cat_type) == '') CYCLE
          n=n+1
 
       ENDDO
       WRITE(*,*) '--------------------------------------------'
       WRITE(*,*) 'NUMBER OF TAGGING Catagories :', n
       WRITE(*,*) '--------------------------------------------'
    ENDIF
    
    
    
! op_mm_20171016+
 
   
 CALL p_bcast(n, p_io)
 ALLOCATE(CATEGORY_SET(n))



! do the string cracking (similar as in offemis) 
IF (p_parallel_io) THEN
   n = 0
   DO i=1, N_MAX_CATEGORY
      IF (TRIM(CATEGORY_IN(i)%abr) == '') CYCLE
      IF (TRIM(CATEGORY_IN(i)%long_name) == '') CYCLE
      IF (TRIM(CATEGORY_IN(i)%cat_type) == '') CYCLE
      n=n+1
      ! op_mm_20140116+
      ! crack string category_types
        
      CALL parse_catstr(status, MAXPSTRLEN,CATEGORY_IN(i)%cat_type, category &
                       ,ncat)
      IF (status /= 0) CALL error_bi('parse_catstr reported an error' &
                                     ,substr)

      CATEGORY_SET(n)%ncat = ncat
      ALLOCATE(CATEGORY_SET(n)%category(CATEGORY_SET(n)%ncat))

      ! initial / default values
      CATEGORY_SET(n)%abr = ''
      CATEGORY_SET(n)%name = ''
      CATEGORY_SET(n)%category(:) = 0


      ! loop over categories
      do j=1,  CATEGORY_SET(n)%ncat
         CATEGORY_SET(n)%category(j)  = category(j)
         !save the ids for later use!
         !set istr  
         if (category(j) .eq. 1) then   
            if (istr .ne. 0) then 
               CALL error_bi('istr is used serveral times' &
                            ,substr)
            end if
            istr=n 
         end if
         !ilig
         if (category(j) .eq. 2) then   
            if (ilig .ne. 0) then 
               CALL error_bi('ilig is used serveral times' &
                    ,substr)
            end if
            ilig=n 
         end if

         !iN2O
         if (category(j) .eq. 3) then   
            if (in2o .ne. 0) then 
               CALL error_bi('in2o is used serveral times' &
                    ,substr)
            end if
            in2o=n 
         end if

            !isoi
         if (category(j) .eq. 4) then   
            if (isoi .ne. 0) then 
               CALL error_bi('isoi is used serveral times' &
                    ,substr)
            end if
            isoi=n 
         end if

         !ich4 
         if (category(j) .eq. 5) then   
            if (ich4 .ne. 0) then 
               CALL error_bi('ich4 is used serveral times' &
                    ,substr)
            end if
            ich4=n 
         end if
         

      end do
         
      ! op_mm_20140116-
      CATEGORY_SET(n)%abr     = CATEGORY_IN(i)%abr
      CATEGORY_SET(n)%name   = CATEGORY_IN(i)%long_name
   ENDDO
   
   ! deallocate category (allocated in pars_catstr)
   deallocate(category)
   
   ! op_mm_20171016
   !set number of species 
   N_tag_species=n
   
   ! op_mm_20171016
   ! check if all necessary categiories are defined in namelist 
   if (istr .eq. 0) then
      CALL error_bi('stratospher category is not defined in namelist' &
           ,substr)
   end if

   if (ilig .eq. 0) then
      CALL error_bi('Lightning category is not defined in namelist' &
           ,substr)
   end if

      
   if (iN2O .eq. 0) then
      CALL error_bi('N2O category is not defined in namelist' &
           ,substr)
   end if   


   if (isoi .eq. 0) then
      CALL error_bi('Soil/biogenic category is not defined in namelist' &
           ,substr)
   end if

   if (iCH4 .eq. 0) then
      CALL error_bi('CH4 category is not defined in namelist' &
           ,substr)
   end if   


ENDIF ! end parallel io 





!broadcast of all values 
CALL p_bcast(N_tag_species, p_io)
call p_bcast(istr, p_io)
call p_bcast(ilig,p_io)
call p_bcast(in2o,p_io)
call p_bcast(isoi,p_io)
call p_bcast(ich4,p_io) 


DO n=1,N_tag_species
   CALL p_bcast(CATEGORY_SET(n)%abr, p_io)
   CALL p_bcast(CATEGORY_SET(n)%name, p_io)
   CALL p_bcast(CATEGORY_SET(n)%ncat, p_io)
End DO

DO n=1,N_tag_species
   IF (.NOT. p_parallel_io) THEN
      ALLOCATE(CATEGORY_SET(n)%category(CATEGORY_SET(n)%ncat))
   ENDIF
end do 


DO n=1,N_tag_species
   do j=1,  CATEGORY_SET(n)%ncat
      CALL p_bcast(CATEGORY_SET(n)%category(j) , p_io)
   end do
end do


!allocate tagging types 
!qqx allocate later?? 
! op_mm_20171016
!no some of them are used in tagging_new_tracer! 
ALLOCATE(o3_tag(N_tag_species))
ALLOCATE(noy_tag(N_tag_species))
ALLOCATE(pan_tag(N_tag_species))
ALLOCATE(nmhc_tag(N_tag_species))
ALLOCATE(co_tag(N_tag_species))
ALLOCATE(oh_tag(N_tag_species))
ALLOCATE(ho2_tag(N_tag_species))
ALLOCATE(ch4loss_tag(N_tag_species))
ALLOCATE(frac_odd(N_tag_species))
ALLOCATE(p_ho2_ten(N_tag_species))
ALLOCATE(p_ro2_ten(N_tag_species))
ALLOCATE(d_oh_ten(N_tag_species))
ALLOCATE(d_ho2_ten(N_tag_species))
ALLOCATE(d_no_ten(N_tag_species))
ALLOCATE(d_ro_ten(N_tag_species))
ALLOCATE(d_xo_ten(N_tag_species))
ALLOCATE(itrac_NOy(N_tag_species))
ALLOCATE(itrac_NMHC(N_tag_species))
ALLOCATE(itrac_PAN(N_tag_species))
ALLOCATE(itrac_CO(N_tag_species))
ALLOCATE(itrac_O3(N_tag_species))

 ! op_mm_20171016-
    
    

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL tagging_read_nml_e5(status, iou)
       IF (status /= 0) CALL error_bi('Error in Reading Namelist CPL', substr)
    END IF
   ! CALL p_bcast(c_tropop(1), p_io)
   ! CALL p_bcast(c_tropop(2), p_io)
    CALL p_bcast(c_lnox(1), p_io)
    CALL p_bcast(c_lnox(2), p_io)


    CALL p_bcast(i_diag, p_io)
!   CALL p_bcast(c_a(1), p_io)
!   CALL p_bcast(c_a(2), p_io)
!   CALL p_bcast(c_b(1), p_io)
!   CALL p_bcast(c_b(2), p_io)

!   Initialize variables i_trac_HNO3 and i_trac_HCHO
    ip_HNO3=0
    do ispecie=1,Fam_NOy_species
       if (TRIM(NOy_species(ispecie))=='HNO3') ip_HNO3=ispecie
    enddo
    IF (ip_HNO3==0)call warning_bi('Warning: HNO3 not found in NOy list ',substr)

     ip_HCHO=0
    do ispecie=1,Fam_NMHC_species
       if (TRIM(NMHC_species(ispecie))=='HCHO') ip_HCHO=ispecie
    enddo
    IF (ip_HCHO==0) call warning_bi(' Warning: CHCO not found in NMHC list',substr)

    ltagging=.false.
    force_reinitialize_tracer=.true.

! op_mm_20131126
!completly changed restart algorithm 
! the new part can now be found in init_tracer    
   !! ltagging=.true.                             !  Switch on the tagging algorithm
   !! tagging_lforce_init=(abs(i_tracer_init)==2)
   !! IF (lstart) then 
   !!    linit_tagging_tracer=.true.              !  Do the initialisation
   !!    linit_tagging_tracer_now=.false.         !  But only on second timestep
   !!    ltagging=.false.                         !  Switch off tagging until initialised
   !! else if (lresume.and..not.i_tracer_init==0) then
   !!    linit_tagging_tracer=.true.
   !!    linit_tagging_tracer_now=.true.
   !! else
   !!    linit_tagging_tracer=.false.             ! no initialisation
   !!    linit_tagging_tracer_now=.false.
   !! endif



   ! Write (*,*) 'lforce=',tagging_lforce_init         
   ! write (*,*) 'lstart_2', lstart
   ! write (*,*) 'ltagging2', ltagging
   ! write (*,*) 'linit_tracer_',linit_tagging_tracer
   ! write (*,*) 'linit_tracer_now',linit_tagging_tracer_now 


    CALL end_message_bi(modstr, 'TAGGING INITIALIZATION', substr)

!   !
  END SUBROUTINE tagging_initialize
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE tagging_new_tracer

    ! ECHAM5/MESSy
!changed finish -> error_bi
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, OFF, ON !, LGTRSTR
    USE messy_main_tracer,        ONLY: ON, OFF, SINGLE, AIR

    USE messy_main_constants_mem, ONLY: MO, MN, MC
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, set_tracer,AIR, ON,OFF &
                                        ,I_DRYDEP, I_SEDI,I_SCAV           &
                                        ,R_MOLARMASS, Amountfraction, i_advect 

    IMPLICIT NONE

    ! LOCAL
    INTEGER(I4) :: status, i
    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_new_tracer'
    INTEGER                     :: jt
    !op_mm_16102017
    
    !CHARACTER(LEN=10),Dimension(10):: name
    !CHARACTER(LEN=30),Dimension(10):: long_name
    CHARACTER(LEN=10),Dimension(N_tag_species):: name
    CHARACTER(LEN=30),Dimension(N_tag_species):: long_name
    

    CALL start_message_bi(modstr, 'TRACER REQUEST', substr)
     status=0


!   ! ALLOCATE INDEX ARRAYS
!   IF (ANY(L_GP(:))) THEN
!      ALLOCATE(idx(I_NTRAC))
!      idx(:) = 0
!   END IF
!   IF (ANY(L_LG(:))) THEN
!      ALLOCATE(idx_a(I_NTRAC))
!      idx_a(:) = 0
!   END IF

!   tracer_loop: DO jt=1, I_NTRAC     ! I_NTRAC=5 wo definiert ???VG Vars ueber CPLnml? 
!
!      IF ((.NOT. L_GP(jt)) .AND. (.NOT. L_LG(jt))) CYCLE

!      IF (p_parallel_io) THEN
!         IF (TRIM(C_NAME(jt)%subname) ==  '') THEN
!            WRITE(*,*) 'TRACER: ', TRIM(C_NAME(jt)%name)
!         ELSE
!            WRITE(*,*) 'TRACER: ', &
!                 TRIM(C_NAME(jt)%name)//'_'//TRIM(C_NAME(jt)%subname)
!         END IF
!      END IF




!definition NOy
!op_mm_16102017
!name= (/'NOylig    ','NOybio    ','NOysoi    ','NOyind    ','NOytra    ', & 
!        'NOyshp    ','NOyair    ','NOyn2o    ','NOych4    ','NOystr    ' /) 

!long_name= (/'Lightning NOy                 ',&
!             'Biomass Burning NOy           ',&
!             'Soils NOy                     ',&
!             'Industry NOy                  ',&
!             'Road traffic NOy              ',&
!             'Ships NOy                     ',&
!             'Air traffic NOy               ',&
!             'N2O Deg. + upper boundary NOy ',&
!             'Impact from CH4 degradation   ',&
!             'Impact from Stratosphere      '/) 




noy_tag(:)%idx=0


do i=1, N_tag_species



   noy_tag(i)%name=trim("NOy") // trim(CATEGORY_SET(i)%abr)
   noy_tag(i)%longname=trim("NOy of ") // trim(CATEGORY_SET(i)%name)



   CALL new_tracer(status, GPTRSTR, trim(noy_tag(i)%name), modstr,  &
        idx=noy_tag(i)%idx,             &
        unit='mol/mol', medium=AIR,    &
        quantity = AMOUNTFRACTION,  &
        longname=trim(noy_tag(i)%longname)              )
      
   CALL set_tracer(status, GPTRSTR, noy_tag(i)%idx, I_DRYDEP,OFF)
   CALL set_tracer(status, GPTRSTR, noy_tag(i)%idx, I_SCAV,OFF)
   CALL set_tracer(status, GPTRSTR, noy_tag(i)%idx, I_SEDI,OFF)
   CALL set_tracer(status, GPTRSTR, noy_tag(i)%idx, R_molarmass,MN)
   IF (p_parallel_io) WRITE(*,*) noy_tag(i)%longname,' status',status
   IF (status/=0) CALL error_bi('NOY Tracer not allocated ', substr)
   IF (p_parallel_io) WRITE(*,*) ' ... New tracer', noy_tag(i)%longname, ' defined'   


end do 


! old call 

!  CALL new_tracer_old(status, GPTRSTR,'NOylig',modstr       &
!             ,idx=itrac_noy_lig, unit='mol/mol'                 &
!              ,ndrydep=OFF, nwetdep=OFF, nscav=OFF, nsedi=OFF    &
!               ,medium=AIR, molarmass=MN                          &
!               ,lforce_init=tagging_lforce_init                   &
!               ,longname='Lightning NOy')
!      IF (p_parallel_io) WRITE(*,*) 'NOylig status',status
!       IF (status/=0) CALL FINISH(substr,'NOylig not allocated ')
!       IF (p_parallel_io) WRITE(*,*) ' ... New tracer NOy_lig defined'






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of NMHC tracers 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!name= (/'NMHClig   ','NMHCbio   ','NMHCsoi   ','NMHCind   ','NMHCtra   ', & 
!        'NMHCshp   ','NMHCair   ','NMHCn2o   ','NMHCch4   ','NMHCstr   ' /) 


!long_name= (/'Lightning NMHC                ',&
!             'Biomass Burning NMHC          ',&
!             'Soils NMHC                    ',&
!             'Industry NMHC                 ',&
!             'Road traffic NMHC             ',&
!             'Ships NMHC                    ',&
!             'Air traffic NMHC              ',&
!             'N2O Deg. + upper boundary NMHC',&
!             'Impact from CH4 degradation   ',&
!             'Impact from Stratosphere      '/) 


nmhc_tag(:)%idx=0

do i=1, N_tag_species

  nmhc_tag(i)%name=trim("NMHC") // trim(CATEGORY_SET(i)%abr)
  nmhc_tag(i)%longname=trim("NMHC of ") // trim(CATEGORY_SET(i)%name)




   CALL new_tracer(status, GPTRSTR, trim(nmhc_tag(i)%name), modstr,  &
                  idx=nmhc_tag(i)%idx,                     &
                  quantity= AMOUNTFRACTION, &
                  unit='mol/mol', medium=AIR,    &
                  longname=trim(nmhc_tag(i)%longname)              )
      
   CALL set_tracer(status, GPTRSTR, nmhc_tag(i)%idx, I_DRYDEP,OFF)
   CALL set_tracer(status, GPTRSTR, nmhc_tag(i)%idx, I_SCAV,OFF)
   CALL set_tracer(status, GPTRSTR, nmhc_tag(i)%idx, I_SEDI,OFF)
   CALL set_tracer(status, GPTRSTR, nmhc_tag(i)%idx, R_molarmass,MC)
   IF (p_parallel_io) WRITE(*,*) nmhc_tag(i)%longname,' status',status
   IF (status/=0) CALL error_bi( 'NMHC tracer not allocated ', substr)
   IF (p_parallel_io) WRITE(*,*) ' ... New tracer', nmhc_tag(i)%longname , ' defined'   

end do 

! old call

!       CALL new_tracer_old(status, GPTRSTR,'NMHCbio',modstr           &
!               ,idx=itrac_nmhc_bio, unit='mol/mol'                 &
!              ,ndrydep=OFF, nwetdep=OFF, nscav=OFF, nsedi=OFF      &
!             ,medium=AIR, molarmass=MC                          &
!             ,lforce_init=tagging_lforce_init                   &
!             ,longname='Biomass Burning NMHC')
!      IF (status/=0) CALL FINISH(substr,'NMHC_bio not allocated ')
!       IF (p_parallel_io) WRITE(*,*) ' ... New tracer NMHC_bio defined'



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of CO tracers 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!name= (/'COlig     ','CObio     ','COsoi     ','COind     ','COtra     ', & 
!        'COshp     ','COair     ','COn2o     ','COch4     ','COstr     ' /) 


!long_name= (/'Lightning CO                  ',&
!             'Biomass Burning CO            ',&
!             'Soils CO                      ',&
!             'Industry CO                   ',&
!             'Road traffic CO               ',&
!             'Ships CO                      ',&
!             'Air traffic CO                ',&
!             'N2O Deg. + upper boundary CO  ',&
!             'Impact from CH4 degradation   ',&
!             'Impact from Stratosphere      '/) 


co_tag(:)%idx=0

do i=1, N_tag_species

   co_tag(i)%name=trim("CO") // trim(CATEGORY_SET(i)%abr)
   co_tag(i)%longname=trim("CO of ") // trim(CATEGORY_SET(i)%name)


   CALL new_tracer(status, GPTRSTR, trim(co_tag(i)%name), modstr,  &
                   idx=co_tag(i)%idx,                     &
                   quantity= AMOUNTFRACTION, &
                   unit='mol/mol', medium=AIR,    &
                   longname=trim(co_tag(i)%longname)              )
      
   CALL set_tracer(status, GPTRSTR, co_tag(i)%idx, I_DRYDEP,OFF)
   CALL set_tracer(status, GPTRSTR, co_tag(i)%idx, I_SCAV,OFF)
   CALL set_tracer(status, GPTRSTR, co_tag(i)%idx, I_SEDI,OFF)
   CALL set_tracer(status, GPTRSTR, co_tag(i)%idx, R_molarmass,MC)
   IF (p_parallel_io) WRITE(*,*) co_tag(i)%longname,' status',status
   IF (status/=0) CALL error_bi( ' co not allocated ', substr)
   IF (p_parallel_io) WRITE(*,*) ' ... New tracer', co_tag(i)%longname, ' defined'   

end do 

! old call
 
!       CALL new_tracer_old(status, GPTRSTR,'COlig',modstr           &
 !              ,idx=itrac_co_lig, unit='mol/mol'                 &
  !             ,ndrydep=OFF, nwetdep=OFF, nscav=OFF, nsedi=OFF      &
 !              ,medium=AIR, molarmass=MC                          &
  !             ,lforce_init=tagging_lforce_init                   &
  !             ,longname='Lightning CO')
   !    IF (status/=0) CALL FINISH(substr,'CO_lig not allocated ')
    !   IF (p_parallel_io) WRITE(*,*) ' ... New tracer CO_lig defined'





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of PAN tracers 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

!name= (/'PANlig    ','PANbio    ','PANsoi    ','PANind    ','PANtra    ', & 
!        'PANshp    ','PANair    ','PANn2o    ','PANch4    ','PANstr    ' /) 


!long_name= (/'Lightning PAN                 ',&
!             'Biomass Burning PAN           ',&
!             'Soils PAN                     ',&
!             'Industry PAN                  ',&
!             'Road traffic PAN              ',&
!             'Ships PAN                     ',&
!             'Air traffic PAN               ',&
!             'N2O Deg. + upper boundary PAN ',&
!             'Impact from CH4 degradation   ',&
!             'Impact from Stratosphere      '/) 



pan_tag(:)%idx=0


do i=1, N_tag_species


!   pan_tag(i)%name=name(i)
!   pan_tag(i)%longname=long_name(i)

   pan_tag(i)%name=trim("PAN") // trim(CATEGORY_SET(i)%abr)
   pan_tag(i)%longname=trim("PAN of ") // trim(CATEGORY_SET(i)%name)


   CALL new_tracer(status, GPTRSTR,trim(pan_tag(i)%name), modstr,  &
                   idx=pan_tag(i)%idx,                     &
                   quantity= AMOUNTFRACTION,   &
                   unit='mol/mol', medium=AIR,    &
                   longname=trim(pan_tag(i)%longname)              )
      
   CALL set_tracer(status, GPTRSTR, pan_tag(i)%idx, I_DRYDEP,OFF)
   CALL set_tracer(status, GPTRSTR, pan_tag(i)%idx, I_SCAV,OFF)
   CALL set_tracer(status, GPTRSTR, pan_tag(i)%idx, I_SEDI,OFF)
   CALL set_tracer(status, GPTRSTR, pan_tag(i)%idx, R_molarmass,MC)
   IF (p_parallel_io) WRITE(*,*) pan_tag(i)%longname,' status',status
         IF (status/=0) CALL error_bi('pan not allocated', substr)
   IF (p_parallel_io) WRITE(*,*) ' ... New tracer', pan_tag(i)%longname, ' defined'   

end do 

! old call
 

!       CALL new_tracer_old(status, GPTRSTR,'PANlig',modstr           &
!               ,idx=itrac_pan_lig, unit='mol/mol'                 &
!               ,ndrydep=OFF, nwetdep=OFF, nscav=OFF, nsedi=OFF      &
!               ,medium=AIR, molarmass=MC                          &
!               ,lforce_init=tagging_lforce_init                   &
 !              ,longname='Lightning PAN')
 !      IF (status/=0) CALL FINISH(substr,'PAN_lig not allocated ')
 !      IF (p_parallel_io) WRITE(*,*) ' ... New tracer PAN_lig defined'


     



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of O3 tracers 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!name= (/'O3lig     ','O3bio     ','O3soi     ','O3ind     ','O3tra     ', & 
!        'O3shp     ','O3air     ','O3n2o     ','O3ch4     ','O3str     ' /) 


!long_name= (/'Lightning O3                  ',&
!             'Biomass Burning O3            ',&
!             'Soils O3                      ',&
!             'Industry O3                   ',&
!             'Road traffic O3               ',&
!             'Ships O3                      ',&
!             'Air traffic O3                ',&
!             'N2O Deg. + upper boundary O3  ',&
!             'Impact from CH4 degradation   ',&
!             'Impact from Stratosphere      '/) 


o3_tag(:)%idx=0



do i=1, N_tag_species


!!   o3_tag(i)%name=name(i)
!!   o3_tag(i)%longname=long_name(i)

   o3_tag(i)%name=trim("O3") // trim(CATEGORY_SET(i)%abr)
   o3_tag(i)%longname=trim("O3 of") // trim(CATEGORY_SET(i)%name)


   
   

   CALL new_tracer(status, GPTRSTR, trim(o3_tag(i)%name), modstr,  &
                   idx=o3_tag(i)%idx,                     &
                    quantity= AMOUNTFRACTION, &
                   unit='mol/mol', medium=AIR,    &
                   longname=trim(o3_tag(i)%longname)              )
      
   CALL set_tracer(status, GPTRSTR, o3_tag(i)%idx, I_DRYDEP,OFF)
   CALL set_tracer(status, GPTRSTR, o3_tag(i)%idx, I_SCAV,OFF)
   CALL set_tracer(status, GPTRSTR, o3_tag(i)%idx, I_SEDI,OFF)
   CALL set_tracer(status, GPTRSTR, o3_tag(i)%idx, R_molarmass,3.*MO)
   IF (p_parallel_io) WRITE(*,*) o3_tag(i)%longname,' status',status
   IF (status/=0) CALL error_bi( ' O3 Tracer not allocated ', substr)
   IF (p_parallel_io) WRITE(*,*) ' ... New tracer', o3_tag(i)%longname, ' defined'   

end do 




 
!       CALL new_tracer_old(status, GPTRSTR,'O3lig',modstr           &
!               ,idx=itrac_o3_lig, unit='mol/mol'                 &
!               ,ndrydep=OFF, nwetdep=OFF, nscav=OFF, nsedi=OFF      &
!               ,medium=AIR, molarmass=3.*MO                          &
!               ,lforce_init=tagging_lforce_init                   &
!               ,longname='Lightning O3')
!       IF (status/=0) CALL FINISH(substr,'O3_lig not allocated ')
!       IF (p_parallel_io) WRITE(*,*) ' ... New tracer O3_lig defined'



!changed
!new tracer old -> new tracer for error diagnostics


!-------------------
! Error Diagnostics
!-------------------
IF (l_adv_err_diag) then
!
!Positve Error O3
!
   CALL new_tracer(status, GPTRSTR, 'O3errP', modstr,  &
        idx=itrac_o3_errP,                         &
        unit='mol/mol', medium=AIR,    &
        longname= 'O3 positive error from numerical diffusion by advection'          )
      
   CALL set_tracer(status, GPTRSTR,itrac_o3_errP , I_DRYDEP,OFF)
   CALL set_tracer(status, GPTRSTR,itrac_o3_errP , I_SCAV,OFF)
   CALL set_tracer(status, GPTRSTR,itrac_o3_errP, I_SEDI,OFF)
   CALL set_tracer(status, GPTRSTR,itrac_o3_errP, R_molarmass,3.*MO)
   IF (p_parallel_io) WRITE(*,*) ' O3 positive error from numerical diffusion by advection status',status
   IF (status/=0) CALL error_bi( 'O3 positive error tracer not allocated ', substr)
   IF (p_parallel_io) WRITE(*,*) ' ... New tracer O3 positive error defined'   

!
!Negative Error O3
!
   CALL new_tracer(status, GPTRSTR, 'O3errN', modstr,             &
                   idx=itrac_o3_errN,                             &
                   unit='mol/mol', medium=AIR,                    &
                   longname= 'O3 positive error from numerical diffusion by advection' )
      
   CALL set_tracer(status, GPTRSTR,itrac_o3_errN , I_DRYDEP,OFF)
   CALL set_tracer(status, GPTRSTR,itrac_o3_errN , I_SCAV,OFF)
   CALL set_tracer(status, GPTRSTR,itrac_o3_errN, I_SEDI,OFF)
   CALL set_tracer(status, GPTRSTR,itrac_o3_errN, R_molarmass,3.*MO)
   IF (p_parallel_io) WRITE(*,*) ' O3 negative  error from numerical diffusion by advection status',status
   IF (status/=0) CALL error_bi(' New tracer O3 negative  error from numerical diffusion by advection not allocated ',substr)
   IF (p_parallel_io) WRITE(*,*) ' ... New tracer O3 negative error defined'   
         
!
!POSITIVE ERROR PAN
!
          CALL new_tracer(status, GPTRSTR, 'PANerrP', modstr,             &
                          idx=itrac_pan_errP,                             &
                          unit='mol/mol', medium=AIR,                    &
                           longname= 'PAN positive error from numerical diffusion by advection' )
      
          CALL set_tracer(status, GPTRSTR,itrac_pan_errP , I_DRYDEP,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_pan_errP , I_SCAV,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_pan_errP, I_SEDI,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_pan_errP, R_molarmass,MC)
          IF (p_parallel_io) WRITE(*,*) 'PAN positive error from numerical diffusion by advection  status',status
          IF (status/=0) CALL error_bi('PAN positive error from numerical diffusion by advection not allocated ', substr)
          IF (p_parallel_io) WRITE(*,*) ' ... New tracer PAN positive error  define'   
         


!Negative ERROR PAN

          CALL new_tracer(status, GPTRSTR, 'PANerrN', modstr,             &
                          idx=itrac_pan_errN,                             &
                          unit='mol/mol', medium=AIR,                    &
                          longname= 'PAN negative error from numerical diffusion by advection' )
      
          CALL set_tracer(status, GPTRSTR,itrac_pan_errN , I_DRYDEP,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_pan_errN , I_SCAV,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_pan_errN, I_SEDI,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_pan_errN, R_molarmass,MC)
          IF (p_parallel_io) WRITE(*,*) 'PAN negative error from numerical diffusion by advection  status',status
          IF (status/=0) CALL error_bi('PAN negative error from numerical diffusion by advection not allocated ', substr)
          IF (p_parallel_io) WRITE(*,*) ' ... New tracer PAN negative error  define'   
         



!Positve Error CO

          CALL new_tracer(status, GPTRSTR, 'COerrP', modstr,             &
                          idx=itrac_co_errP,                             &
                          unit='mol/mol', medium=AIR,                    &
                           longname= 'CO positive error from numerical diffusion by advection' )
      
          CALL set_tracer(status, GPTRSTR,itrac_co_errP , I_DRYDEP,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_co_errP , I_SCAV,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_co_errP, I_SEDI,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_co_errP, R_molarmass,MC)
          IF (p_parallel_io) WRITE(*,*) 'CO positive error from numerical diffusion by advection  status',status
          IF (status/=0) CALL error_bi('CO positive error  from numerical diffusion by advection not allocated ', substr)
          IF (p_parallel_io) WRITE(*,*) ' ... New tracer CO positive   define'   
         

!
!Negative Error CO
!
          CALL new_tracer(status, GPTRSTR, 'COerrN', modstr,             &
                          idx=itrac_co_errN,                             &
                          unit='mol/mol', medium=AIR,                    &
                           longname= 'CO negative error from numerical diffusion by advection' )
      
          CALL set_tracer(status, GPTRSTR,itrac_co_errN , I_DRYDEP,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_co_errN , I_SCAV,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_co_errN, I_SEDI,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_co_errN, R_molarmass,MC)
          IF (p_parallel_io) WRITE(*,*) 'CO negative error from numerical diffusion by advection  status',status
          IF (status/=0) CALL error_bi('CO negative error  from numerical diffusion by advection not allocated ', substr)
          IF (p_parallel_io) WRITE(*,*) ' ... New tracer CO negative  error  define' 
  
!         
!Poitive error NMHC
!        

          CALL new_tracer(status, GPTRSTR, 'NMHCerrP', modstr,             &
                          idx=itrac_nmhc_errP,                             &
                          unit='mol/mol', medium=AIR,                    &
                          longname= 'NMHC positive error from numerical diffusion by advection' )
      
          CALL set_tracer(status, GPTRSTR,itrac_nmhc_errP , I_DRYDEP,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_nmhc_errP , I_SCAV,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_nmhc_errP, I_SEDI,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_nmhc_errP, R_molarmass,MC)
          IF (p_parallel_io) WRITE(*,*) 'NMHC positive error from numerical diffusion by advection  status',status
          IF (status/=0) CALL error_bi('NMHC positive error  from numerical diffusion by advection not allocated ', substr)
          IF (p_parallel_io) WRITE(*,*) ' ... New tracer NMHC positive  define'   

!
!Negative error NMHC
!
 
         CALL new_tracer(status, GPTRSTR, 'NMHCerrN', modstr,             &
                          idx=itrac_nmhc_errN,                             &
                          unit='mol/mol', medium=AIR,                    &
                          longname= 'NMHC negative error from numerical diffusion by advection' )
      
          CALL set_tracer(status, GPTRSTR,itrac_nmhc_errN , I_DRYDEP,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_nmhc_errN , I_SCAV,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_nmhc_errN, I_SEDI,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_nmhc_errN, R_molarmass,MC)
          IF (p_parallel_io) WRITE(*,*) 'NMHC negative error from numerical diffusion by advection  status',status
          IF (status/=0) CALL error_bi('NMHC negative error  from numerical diffusion by advection not allocated ', substr)
          IF (p_parallel_io) WRITE(*,*) ' ... New tracer NMHC negative error  define'   

!
!Poitive error NOy
!

          CALL new_tracer(status, GPTRSTR, 'NOyerrP', modstr,             &
                          idx=itrac_noy_errP,                             &
                          unit='mol/mol', medium=AIR,                    &
                          longname= 'NOy positive error from numerical diffusion by advection' )
      
          CALL set_tracer(status, GPTRSTR,itrac_noy_errP , I_DRYDEP,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_noy_errP , I_SCAV,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_noy_errP, I_SEDI,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_noy_errP, R_molarmass,MC)
          IF (p_parallel_io) WRITE(*,*) 'NOy positive error from numerical diffusion by advection  status',status
          IF (status/=0) CALL error_bi('NOy positive error  from numerical diffusion by advection not allocated ', substr)
          IF (p_parallel_io) WRITE(*,*) ' ... New tracer NOy  positive error define'  


!
!Negative  error NOy
!
          CALL new_tracer(status, GPTRSTR, 'NOyerrN', modstr,             &
                          idx=itrac_noy_errN,                             &
                          unit='mol/mol', medium=AIR,                    &
                          longname= 'NOy negative error from numerical diffusion by advection' )
      
          CALL set_tracer(status, GPTRSTR,itrac_noy_errN , I_DRYDEP,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_noy_errN , I_SCAV,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_noy_errN, I_SEDI,OFF)
          CALL set_tracer(status, GPTRSTR,itrac_noy_errN, R_molarmass,MC)
          IF (p_parallel_io) WRITE(*,*) 'NOy negative error from numerical diffusion by advection  status',status
          IF (status/=0) CALL error_bi('NOy negative error  from numerical diffusion by advection not allocated ', substr)
          IF (p_parallel_io) WRITE(*,*) ' ... New tracer NOy  negative error  define'  



ENDIF ! Err_diag
   
CALL end_message_bi(modstr, 'TRACER REQUEST', substr)

END SUBROUTINE tagging_new_tracer


! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE tagging_init_tracer

    ! ECHAM5/MESSy
    USE messy_main_tracer_bi,        ONLY: main_tracer_init_tracer
    USE messy_main_tracer,           ONLY: tracer_iniflag
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR

    USE messy_main_mpi_bi,           ONLY: p_parallel_io

    IMPLICIT NONE

    INTEGER                        :: status, ispec

    LOGICAL                         :: isinit
    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_init_tracer'

    status=0
!    IF (.NOT.linit_tagging_tracer_now) RETURN     ! if initialisation should be done later 
                                           !  i.e. in case of new start
                                           !  Might be wrong if data are read from file
!  linit=.true. 

! old code 
!-------
! Check whether tracer were initialised from re-run file, 
!------
!      CALL tagging_init_tracer_online  
    !print*, "in init_tracer"

   ! SELECT CASE(i_tracer_init)
   !  CASE (0)           ! tracer assumed to be read from re-run file 
   !   Return
   !  CASE (-1,1)          ! discard initialisation by re-run file 
                        ! and do it again


   !   IF (p_parallel_io) WRITE(*,*) modstr,' TRACER INIT: uniformly  done '
   !   RETURN
   !  CASE (-2,2)          ! 
!      CALL tracer_init(modstr)    ! per input files
   !    CALL  main_tracer_init_tracer(4) 
   !   write (*,*)'tagging_init_tracer_2ndcase', (abs(i_tracer_init)==2)
   !  CASE DEFAULT
   !   WRITE(*,*) 'WARNING Tracer initialization for tagging failed'
   !   WRITE(*,*) ' Error in i_tracer_init=',i_tracer_init
   !   status=1
   ! END SELECT

    ! ----
    ! Error tracer are automatically set to zero!
    ! ----

! op_mm_20150204
! Changed for new restart logic 
! First loop over all tracers an see if they have been initialised (from rerun file or tracer.nml)

isinit = .FALSE. 

  DO ispec=1,N_tag_species

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! NOY TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
      call tracer_iniflag(status, GPTRSTR, noy_tag(ispec)%idx, lget=isinit)
      if (status/=0) write(*,*) 'GET tracer_iniflag failed for NOy ',ispec,status
      
      if (.not. isinit)  exit 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! O3 TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     call tracer_iniflag(status, GPTRSTR, o3_tag(ispec)%idx, lget=isinit)
     if (status/=0) write(*,*) 'GET tracer_iniflag failed for O3 ',ispec,status
     
     if (.not. isinit) exit 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PAN TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     call tracer_iniflag(status, GPTRSTR, pan_tag(ispec)%idx, lget=isinit)
     if (status/=0) write(*,*) 'GET tracer_iniflag failed for PAN ',ispec,status
     
     if (.not. isinit) exit 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! NMHC TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     call tracer_iniflag(status, GPTRSTR, nmhc_tag(ispec)%idx, lget=isinit)
     if (status/=0) write(*,*) 'GET tracer_iniflag failed for NMHC ',ispec,status

     if (.not. isinit ) exit 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CO TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     call tracer_iniflag(status, GPTRSTR, co_tag(ispec)%idx, lget=isinit)
     if (status/=0) write(*,*) 'GET tracer_iniflag failed for CO ',ispec,status

     if (.not. isinit) exit 

END DO 

IF (ISINIT) THEN 
   write(*,*) "Tagging Tracers have been initialised via rerun or via tracer namelist" 
   force_reinitialize_tracer = .false. 
   ltagging = .true. 
END IF

!This would mean that tagging tracers are initialised (from rerun file or via tracer namelist)
! However the user might want to force the initialisation 
! this can be done with the namelist option i_tracer_init==0 
! or by not specifying the initialisation via tracer.nml

IF (.not. i_tracer_init==0) then
    force_reinitialize_tracer = .true.
     ltagging=.false.
END IF 

 



    IF (p_parallel_io)  Write(*,*) modstr,' TRACER INITIALIZED_init_tracer'


  END SUBROUTINE tagging_init_tracer
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE tagging_init_tracer_online

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, xt, xtte, xtm1, xtf, ti_gp, ntrac_gp
    USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nlev, nproma,  kproma, npromz       
    USE messy_main_grid_def_bi,   ONLY: philat_2d, hyam, hybm
    USE messy_main_data_bi,       ONLY: aps, press_3d
    USE messy_main_tools,         ONLY: iso2ind          
    USE messy_main_constants_mem, ONLY: pi      
    USE messy_main_timer,         ONLY: lstart, lresume
    USE messy_main_tracer,        ONLY: tracer_iniflag !, lset
    USE messy_main_mpi_bi,        ONLY: p_parallel_io

    IMPLICIT NONE

    INTEGER                     :: jt, i, j, jrow, jp
    INTEGER                     :: ispec, is
    INTEGER                     :: status
   ! INTEGER,  PARAMETER        :: ktp=10      ! tp level now calculated online ! op_mm_20140124
    INTEGER                    :: tp           !variabel for local tp index from climatological tropopause

!for cosmo different order of array necessary
    REAL(DP), DIMENSION(_RI_XYZ__(nproma,ngpblks,nlev) ,2,5) :: rts  ! et put 1 
                                !   factor each specie is initialised
                                !   (:,:,1,:)=trop     (:,:,2,:)= strat


    INTEGER,  DIMENSION(5,N_tag_species)         :: its
                                ! how to initialise the fields as strat or trop component
    INTEGER,  PARAMETER         :: io3=1, inoy=2, ico=3, ipan=4, inmhc=5 !  ,  iohe=6, iho2e=7

    LOGICAL,  PARAMETER         :: linit=.TRUE.
    LOGICAL                     :: isinit            ! op_mm_20150203
    REAL(DP)                    :: rstrat,rtrop      ! number of tracers of 
                                                     ! strat, trops. characteristics
    real(DP),ALLOCATABLE         ::zphi(:)
    INTEGER, ALLOCATABLE         ::TP_LEV(:,:)
    INTEGER                      ::iclim
    real(DP)                     ::zpclim


    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_init_tracer_online'

   
!-------
! Initialise Components characteristic
! Check istr !!! haven't done it yet
!------
   its(:,:)=1         ! everything is troposphere
   its(io3,istr)=2    ! except for istr tracer
   its(inoy,in2o)=2   ! except for in2o
   its(ico,ich4)=2    ! except for methane source
   its(ipan,ich4)=2   ! except for methane source
   its(inmhc,ich4)=2  ! except for methane source
 
!---------------
!added ! op_mm_20131109
!calculation of tropopause
!---------
!calculate climatological tropopause (from messy_tropop_si)
   allocate(zphi(nproma)) 
   allocate(tp_LEV(nproma, ngpblks) )

   zpclim=0.0d0
   zphi=0.0d0
   tp_lev=0


   do jrow=1,ngpblks
#ifndef CESM1
      IF ( jrow   == ngpblks ) THEN
         kproma = npromz
      ELSE
         kproma = nproma
      ENDIF
#else
      kproma = npromz(jrow)
#endif
          
      DO jp=1, kproma
         zphi(jp)    = (philat_2d(jp,jrow)/180.)*pi
         zpclim = 300._dp*100._dp &
              - 215._dp*100._dp*cos(zphi(jp))*cos(zphi(jp))

         CALL iso2ind(press_3d(_RI_XYZ__(jp,jrow,:)), zpclim,   &
                    iclim, lrev=.true.)

         tp_lev(jp, jrow) =iclim

         !Debug channel objects
!!$      debug_press(_RI_XYZ__(jp,jrow,:))=press_3d(_RI_XYZ__(jp,jrow,:))    
!!$      debug_zklim(jp, jrow)=tp_lev(jp, jrow)

      end do
   end do

   deallocate(zphi)

! 
! end calculation tropo
!


!--------
! Initialize mask for initialization  SUM(MASK)=1. for each gridpoint
!--------

 do ispec=1,5
        WRITE(*,*) 'Specie', ispec
       ! Contribution for tropospheric tracers in troposphere
        rstrat=REAL(SUM(its(ispec,:))-N_tag_species,DP)      ! number of stratos tracers
        rtrop=REAL(N_tag_species,DP)-rstrat                 ! number of tropospheric tracers
!changed: differente array order for cosmo

        Write(*,*) 'trop/strat tracer',ispec,rtrop,rstrat
        do i=1, nproma
           do j=1,ngpblks
              tp=tp_lev(i,j)
              rts(_RI_XYZ__(i,j,tp+1:nlev),1,ispec)= (1._dp-factor_small_init)/rtrop
             
              ! Contribution for stratospheric tracer in troposphere
              rts(_RI_XYZ__(i,j,tp+1:nlev),2,ispec)=factor_small_init/rstrat
                  
              ! Contribution for tropospheric tracers in stratosphere
              rts(_RI_XYZ__(i,j,1:tp),1,ispec)= factor_small_init/rtrop
         
              ! Contribution for stratospheric tracers in stratosphere
              rts(_RI_XYZ__(i,j,1:tp),2,ispec)=(1._dp-factor_small_init)/rstrat
 
           end do
        end do
     enddo


     write(*,*) 'MINVAL MAXVAL trop',MINVAL(rts(:,:,:,1,:)),MAXVAL(rts(:,:,:,1,:)), &
                 MINVAL(rts(:,:,:,2,:)),MAXVAL(rts(:,:,:,2,:))

  deallocate (tp_lev) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set negative values of mecca1 tracers to zero
!    IF (lstart) THEN ! don't do this after a restart !
!      tracer_loop: DO jt = 1, ntrac_gp
!        IF (TRIM(ti_gp(jt)%tp%ident%submodel) == 'tagging') THEN
!          xt(:,:,jt,:) = 0._dp ! MAX(xt(:,:,jt,:),0._dp)
!          WRITE (*,*) "Tracer set to 0 : " ,jt
!         ENDIF
!      END DO tracer_loop
!    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------
! Set Initialise tracer mixing ratios and set initialisation flag 
! ------

   ! WRITE (*,*) 'Tracer O3total: idx_o3',idx_o3
    
    !WRITE (*,*) 'Tracer O3total: xtm1',MINVAL(xtm1(:,:,idx_o3,:)),&
      !                '...',MAXVAL(xtm1(:,:,idx_o3,:))

    !WRITE (*,*) 'Tracer O3total: xtf',MINVAL(xtf(:,:,idx_o3,:)),&
     !                 '...',MAXVAL(xtf(:,:,idx_o3,:))


   ! WRITE (*,*) 'Tracer O3total: xt ',MINVAL(xt(:,:,idx_o3,:)),&
    !                  '...',MAXVAL(xt(:,:,idx_o3,:))


! op_mm_20131206
! id of tagging tracer is now noy_tag(i)%idx
! idx_noy not changed yet (id of messy tracer from get tracer          

! changed
! xtf not available in cosmo (no time filter)
! check is xtf is associated


! op_mm_20150203
! changed restart logic 

  DO ispec=1,N_tag_species

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! NOY TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        xtte(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))=0._dp 
        xtm1(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))=0._dp
        xt(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))  =0._dp
        IF (ASSOCIATED(xtf)) THEN 
           xtf(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx)) =0._dp
        end if
     
        DO is=1,Fam_NOy_species
           xtm1(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))= xtm1(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))  &
             + xtm1(_RI_XYZN_(:,:,:,idx_NOy(is)))*rts(:,:,:,its(inoy,ispec),inoy)
           xt(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))  = xt(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))    &
                + xt(_RI_XYZN_(:,:,:,idx_NOy(is)))*rts(:,:,:,its(inoy,ispec),inoy)
           xtte(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx)) = xtte(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))   &
                + xtte(_RI_XYZN_(:,:,:,idx_NOy(is)))*rts(:,:,:,its(inoy,ispec),inoy)
           IF (ASSOCIATED(xtf)) THEN   
              xtf(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))= xtf(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))  &
                   + xtf(_RI_XYZN_(:,:,:,idx_NOy(is)))*rts(:,:,:,its(inoy,ispec),inoy)
           end if
        ENDDO     ! Sum over N, NO, NO2, ...
        ! 0. --> pxtte
          
     ! As tracer is now initialized set iniflag
        call tracer_iniflag(status, GPTRSTR, noy_tag(ispec)%idx, linit)
        if (status/=0) write(*,*) 'tracer_iniflag failed for NOy ',ispec,status
        WRITE (*,*) 'xtm1 Tracer NOy: ',ispec,':',MINVAL(xtm1(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))),&
             '...',MAXVAL(xtm1(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx)))
        WRITE (*,*) 'xt   Tracer NOy: ',ispec,':',MINVAL(xt(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))),&
             '...',MAXVAL(xt(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx)))
        IF (ASSOCIATED(xtf)) THEN        
           WRITE (*,*) 'xtf - xmem Tracer NOy: ',ispec,':',MINVAL(xtf(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))),&
                '...',MAXVAL(xtf(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx)))
        end if
             
        WRITE (*,*) 'xtte Tracer NOy: ',ispec,':',MINVAL(xtte(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx))),&
             '...',MAXVAL(xtte(_RI_XYZN_(:,:,:,noy_tag(ispec)%idx)))

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! O3 TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        xtte(_RI_XYZN_(:,:,:,o3_tag(ispec)%idx))=0._dp  ! pxtte(1:kproma,:,idx_o3)*rts(1:kproma,:,its(io3,ispec),io3)
        xtm1(_RI_XYZN_(:,:,:,o3_tag(ispec)%idx))= xtm1(_RI_XYZN_(:,:,:,idx_o3))*rts(:,:,:,its(io3,ispec),io3)
        xt(_RI_XYZN_(:,:,:,o3_tag(ispec)%idx))= xt(_RI_XYZN_(:,:,:,idx_o3))*rts(:,:,:,its(io3,ispec),io3)
        IF (ASSOCIATED(xtf)) THEN   
           xtf(_RI_XYZN_(:,:,:,o3_tag(ispec)%idx))= xtf(_RI_XYZN_(:,:,:,idx_o3))*rts(:,:,:,its(io3,ispec),io3)
        end if
        call tracer_iniflag(status, GPTRSTR, o3_tag(ispec)%idx, linit )
        if (status/=0) write(*,*) 'tracer_iniflag failed for O3 ',ispec,status
        WRITE (*,*) 'Tracer O3: ',ispec,':',MINVAL(xtm1(_RI_XYZN_(:,:,:,o3_tag(ispec)%idx))),&
             '...',MAXVAL(xtm1(_RI_XYZN_(:,:,:,o3_tag(ispec)%idx)))

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PAN TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        xtte(_RI_XYZN_(:,:,:,pan_tag(ispec)%idx))=0._dp  ! pxtte(1:kproma,:,idx_pan)*rts(1:kproma,:,its(ipan,ispec),ipan)
        xtm1(_RI_XYZN_(:,:,:,pan_tag(ispec)%idx))=xtm1(_RI_XYZN_(:,:,:,idx_pan))*rts(:,:,:,its(ipan,ispec),ipan)
        xt(_RI_XYZN_(:,:,:,pan_tag(ispec)%idx))=xt(_RI_XYZN_(:,:,:,idx_pan))*rts(:,:,:,its(ipan,ispec),ipan)
        IF (ASSOCIATED(xtf)) THEN   
           xtf(_RI_XYZN_(:,:,:,pan_tag(ispec)%idx))=xtf(_RI_XYZN_(:,:,:,idx_pan))*rts(:,:,:,its(ipan,ispec),ipan)
        end if
        call tracer_iniflag(status, GPTRSTR, pan_tag(ispec)%idx, linit )
        if (status/=0) write(*,*) 'tracer_iniflag failed for PAN ',ispec,status
        WRITE (*,*) 'Tracer PAN: ',ispec,':',MINVAL(xtm1(_RI_XYZN_(:,:,:,pan_tag(ispec)%idx))),&
             '...',MAXVAL(xtm1(_RI_XYZN_(:,:,:,pan_tag(ispec)%idx)))
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! NMHC TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     xtte(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))=0._dp
     xtm1(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))=0._dp
     xt(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))=0._dp
     IF (ASSOCIATED(xtf)) THEN   
        xtf(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))=0._dp
     end if
   
     DO is=1,Fam_NMHC_species
        xtm1(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))=xtm1(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx)) &
             + xtm1(_RI_XYZN_(:,:,:,idx_NMHC(is))) * rts(:,:,:,its(inmhc,ispec),inmhc)
        xtte(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))=xtte(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))    &
             + xtte(_RI_XYZN_(:,:,:,idx_NMHC(is))) * rts(:,:,:,its(inmhc,ispec),inmhc)
      xt(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))=xt(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))     &
           + xt(_RI_XYZN_(:,:,:,idx_NMHC(is))) * rts(:,:,:,its(inmhc,ispec),inmhc)
      IF (ASSOCIATED(xtf)) THEN 
         xtf(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))=xtf(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))     &
              + xtf(_RI_XYZN_(:,:,:,idx_NMHC(is))) * rts(:,:,:,its(inmhc,ispec),inmhc)
      end if
   END DO
          
   call tracer_iniflag(status, GPTRSTR, nmhc_tag(ispec)%idx, linit)
   if (status/=0) write(*,*) 'tracer_iniflag failed for NMHC ',ispec,status
   WRITE (*,*) 'Tracer NMHC: ',ispec,':',MINVAL(xtm1(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx))),&
               '...',MAXVAL(xtm1(_RI_XYZN_(:,:,:,nmhc_tag(ispec)%idx)))


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CO TRACER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     xtte(_RI_XYZN_(:,:,:,co_tag(ispec)%idx))=0._dp
     xtm1(_RI_XYZN_(:,:,:,co_tag(ispec)%idx))=xtm1(_RI_XYZN_(:,:,:,idx_co))*rts(:,:,:,its(ico,ispec),ico)
     xt(_RI_XYZN_(:,:,:,co_tag(ispec)%idx))=xt(_RI_XYZN_(:,:,:,idx_co))*rts(:,:,:,its(ico,ispec),ico)
     IF (ASSOCIATED(xtf)) THEN      
        xtf(_RI_XYZN_(:,:,:,co_tag(ispec)%idx))=xtf(_RI_XYZN_(:,:,:,idx_co))*rts(:,:,:,its(ico,ispec),ico)
     end if
     call tracer_iniflag(status, GPTRSTR, co_tag(ispec)%idx, linit)
     if (status/=0) write(*,*) 'tracer_iniflag failed for CO ',ispec,status
     WRITE (*,*) 'Tracer CO: ',ispec,':',MINVAL(xtm1(_RI_XYZN_(:,:,:,co_tag(ispec)%idx))),&
          '...',MAXVAL(xtm1(_RI_XYZN_(:,:,:,co_tag(ispec)%idx)))

ENDDO  ! over all categories, lig, bio, etc


IF (p_parallel_io)  Write(*,*) modstr,' TRACER INITIALIZED_online'
END SUBROUTINE tagging_init_tracer_online




! ------------------------------------------------------------------------
SUBROUTINE tagging_init_memory
! Wird eigentlich nicht gebraucht Loss and Production sowie origin tracers ueber 
! Tracer stream definiert! Doch hier kommen die rein diagnostischen vars!!
    !
    ! define Tagging specific stream(s) and allocate memory for
    ! global fields
    !
    ! Author: 
    !         Volker Grewe, DLR, August 2008

    ! ECHAM5
    !USE mo_memory_base,        ONLY: t_stream, new_stream         &
    !                                ,add_stream_element           &
    !                                ,default_stream_setting       
    USE messy_main_mpi_bi,     ONLY: p_parallel_io

    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngl, ngpblks
    USE messy_main_tools,      ONLY:  PTR_3D_ARRAY
    
    USE messy_main_channel,    ONLY: new_channel, new_channel_object, new_attribute

    USE messy_main_channel_BI,        ONLY: GP_2D_HORIZONTAL
    USE messy_main_channel_error_BI,  ONLY: channel_halt

    USE messy_main_channel_repr,      ONLY: get_representation_id &
                                          , get_representation_info
    IMPLICIT NONE


    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_init_memory'

    
    integer(kind=4)                :: i
    integer                        :: status, repr_idx



    !character (len=10), dimension(n_tag_species)     :: setname, longname


! op_mm_20131202
!special fields for storing data before passing to smcl
allocate(oh_tag_data(nproma, nlev, n_tag_species))
allocate(ho2_tag_data(nproma, nlev, n_tag_species))
allocate(p_ho2_ten_data(nproma, nlev, n_tag_species))
allocate(p_ro2_ten_data(nproma, nlev, n_tag_species))
allocate(d_ro2_ten_data(nproma, nlev, n_tag_species))
allocate(d_ro_ten_data(nproma, nlev, n_tag_species))
allocate(d_oh_ten_data(nproma, nlev, n_tag_species))
allocate(d_no_ten_data(nproma, nlev, n_tag_species))
allocate(d_xo_ten_data(nproma, nlev, n_tag_species))
allocate(d_ho2_ten_data(nproma, nlev, n_tag_species))
allocate(ch4loss_tag_data(nproma, nlev, n_tag_species))
allocate(frac_ODD_data(nproma, nlev, n_tag_species))

!ALLOCATE(spec_field_ls(nproma, nsets, nspec, ngpblks, 1))
ALLOCATE(o3_tagging(nproma,nlev,N_tag_species,ngpblks))
ALLOCATE(noy_tagging(nproma,nlev,N_tag_species,ngpblks))
ALLOCATE(nmhc_tagging(nproma,nlev,N_tag_species,ngpblks))
ALLOCATE(co_tagging(nproma,nlev,N_tag_species,ngpblks))
ALLOCATE(pan_tagging(nproma,nlev,N_tag_species,ngpblks))
 

ALLOCATE(o3_tagging_sum(nproma,nlev,ngpblks))
ALLOCATE(noy_tagging_sum(nproma,nlev,ngpblks))
ALLOCATE(nmhc_tagging_sum(nproma,nlev,ngpblks))
ALLOCATE(co_tagging_sum(nproma,nlev,ngpblks))
ALLOCATE(pan_tagging_sum(nproma,nlev,ngpblks))


ALLOCATE(do3(nproma,nlev,N_tag_species,ngpblks))
ALLOCATE(dnoy(nproma,nlev,N_tag_species,ngpblks))
ALLOCATE(dco(nproma,nlev,N_tag_species,ngpblks))
ALLOCATE(dpan(nproma,nlev,N_tag_species,ngpblks))
ALLOCATE(dnmhc(nproma,nlev,N_tag_species,ngpblks))
 
ALLOCATE(o3_0(nproma,nlev,ngpblks))
ALLOCATE(noy_0(nproma,nlev,ngpblks))
ALLOCATE(pan_0(nproma,nlev,ngpblks))
ALLOCATE(nmhc_0(nproma,nlev,ngpblks))
ALLOCATE(co_0(nproma,nlev,ngpblks))

ALLOCATE(o3_1(nproma,nlev,ngpblks))
ALLOCATE(noy_1(nproma,nlev,ngpblks))
ALLOCATE(pan_1(nproma,nlev,ngpblks))
ALLOCATE(nmhc_1(nproma,nlev,ngpblks))
ALLOCATE(co_1(nproma,nlev,ngpblks))
ALLOCATE(ch4_1(nproma,nlev,ngpblks))

ALLOCATE(o3_2(nproma,nlev,ngpblks))
ALLOCATE(noy_2(nproma,nlev,ngpblks))
ALLOCATE(pan_2(nproma,nlev,ngpblks))
ALLOCATE(nmhc_2(nproma,nlev,ngpblks))
ALLOCATE(co_2(nproma,nlev,ngpblks))
ALLOCATE(ch4_2(nproma,nlev,ngpblks))
ALLOCATE(ho2_2(nproma,nlev,ngpblks))
ALLOCATE(oh_2(nproma,nlev,ngpblks))

! just for out out
! ALLOCATE(do3_out(nproma,nlev,ngpblks))
! ALLOCATE(dnoy_out(nproma,nlev,ngpblks))
! ALLOCATE(dpan_out(nproma,nlev,ngpblks))
! ALLOCATE(dnmhc_out(nproma,nlev,ngpblks))
! ALLOCATE(dco_out(nproma,nlev,ngpblks))

ALLOCATE(do3_err(nproma,nlev,ngpblks))
ALLOCATE(dnoy_err(nproma,nlev,ngpblks))
ALLOCATE(dco_err(nproma,nlev,ngpblks))
ALLOCATE(dpan_err(nproma,nlev,ngpblks))
ALLOCATE(dnmhc_err(nproma,nlev,ngpblks))

ALLOCATE(emis_noy_onlem(nproma,nlev,ngpblks))
ALLOCATE(emis_nmhc_onlem(nproma,nlev,ngpblks))

ALLOCATE(pxtte_pre_emis_NOy(nproma,nlev,ngpblks))
ALLOCATE(pxtte_pre_emis_NMHC(nproma,nlev,ngpblks))

ALLOCATE(o3prod_ho2(nproma,nlev,ngpblks))
ALLOCATE(o3prod_ro2(nproma,nlev,ngpblks))
ALLOCATE(o3prod_o2(nproma,nlev,ngpblks))

ALLOCATE(o3loss_ho2(nproma,nlev,ngpblks))
ALLOCATE(o3loss_ro(nproma,nlev,ngpblks))
ALLOCATE(o3loss_oh(nproma,nlev,ngpblks))
ALLOCATE(o3loss_no(nproma,nlev,ngpblks))
ALLOCATE(o3loss_xo(nproma,nlev,ngpblks))

ALLOCATE(panprod(nproma,nlev,ngpblks))
ALLOCATE(panloss(nproma,nlev,ngpblks))
!ALLOCATE(coprod(nproma,nlev,ngpblks))
!ALLOCATE(coloss(nproma,nlev,ngpblks))
ALLOCATE(nmhcloss(nproma,nlev,ngpblks))
! ALLOCATE(ch4loss(nproma,nlev,ngpblks))

ALLOCATE(HO2(nproma,nlev,ngpblks))
ALLOCATE(OH(nproma,nlev,ngpblks))  
ALLOCATE(ch4(nproma,nlev,ngpblks))   
ALLOCATE(co(nproma,nlev,ngpblks))  
ALLOCATE(NMHC(nproma,nlev,ngpblks)) 
ALLOCATE(NOy(nproma,nlev,ngpblks))
ALLOCATE(O3(nproma,nlev,ngpblks))
ALLOCATE(PAN(nproma,nlev,ngpblks))
ALLOCATE(O1D(nproma,nlev,ngpblks))

ALLOCATE(LossG2100(nproma,nlev,ngpblks))
ALLOCATE(LossG2103(nproma,nlev,ngpblks))
ALLOCATE(LossG2104(nproma,nlev,ngpblks))
ALLOCATE(LossG2105(nproma,nlev,ngpblks))
ALLOCATE(LossG2106(nproma,nlev,ngpblks))
ALLOCATE(LossG2107(nproma,nlev,ngpblks))
ALLOCATE(LossG2109(nproma,nlev,ngpblks))
ALLOCATE(LossG2110(nproma,nlev,ngpblks))
ALLOCATE(LossG2111(nproma,nlev,ngpblks))
ALLOCATE(LossG2112(nproma,nlev,ngpblks))
ALLOCATE(LossG3200(nproma,nlev,ngpblks))
ALLOCATE(LossG3201(nproma,nlev,ngpblks))
ALLOCATE(LossG3202(nproma,nlev,ngpblks))
ALLOCATE(LossG3203(nproma,nlev,ngpblks))
ALLOCATE(LossG3207(nproma,nlev,ngpblks))
ALLOCATE(LossG4101(nproma,nlev,ngpblks))
ALLOCATE(LossG4110(nproma,nlev,ngpblks))
ALLOCATE(LossG6203(nproma,nlev,ngpblks))
ALLOCATE(LossG6204(nproma,nlev,ngpblks))
ALLOCATE(LossG7201(nproma,nlev,ngpblks))

ALLOCATE(LossJ3200(nproma,nlev,ngpblks))
ALLOCATE(LossJ3201(nproma,nlev,ngpblks))
ALLOCATE(LossJ4101b(nproma,nlev,ngpblks))

ALLOCATE(OHlossNMHC(nproma,nlev,ngpblks))
ALLOCATE(HO2lossNMHC(nproma,nlev,ngpblks))
ALLOCATE(OHlossHO2prodNMHC(nproma,nlev,ngpblks))
ALLOCATE(HO2prodNMHCNOy(nproma,nlev,ngpblks))
ALLOCATE(HO2prodNMHCphoto(nproma,nlev,ngpblks))



if (i_diag.eq.0) then
   CALL start_message_bi(modstr, 'NO Channel DEFINITION', substr)
   RETURN
endif

CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

! op_mm_20131206
! stream is now channel
! complete new memory managment
! now using type definitions
! looping over n_tag_species

! define channel
CALL new_channel(status, cname=modstr, lrestreq=.true.)
CALL channel_halt(substr, status)


! ELEMENTS USED FOR GP AND LG
IF (p_parallel_io) WRITE(*,*) ' ... '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!changed
! representation for all fields is GP_3D_MID
   CALL get_representation_id(status, 'GP_3D_MID', repr_idx)
   CALL channel_halt(substr, status)


!Debug channel objects for tropop calculation
!!C$ALL new_channel_object(status, modstr, 'tropop', p2=debug_zklim &
!                         , reprid=GP_2D_HORIZONTAL                          )
!!$CALL channel_halt(substr, status)

!!$CALL new_channel_object(status, modstr, 'rts1', p3=debug_rts1 &
!     , reprid=repr_idx                         )
!!$CALL channel_halt(substr, status)
!!$CALL new_channel_object(status, modstr, 'rts2', p3=debug_rts2 &
!     , reprid=repr_idx                         )
!!$CALL channel_halt(substr, status)
!!$CALL new_channel_object(status, modstr, 'rts3', p3=debug_rts3 &
!     , reprid=repr_idx                         )
!!$CALL channel_halt(substr, status)
!!$CALL new_channel_object(status, modstr, 'rts4', p3=debug_rts4 &
!     , reprid=repr_idx                         )
!!$CALL channel_halt(substr, status)


       

! first define channel objects for global tendencies
! a loop would be nice, but doesnt really make sense here
      !
      !O3 loss
      !
      CALL new_channel_object(status, modstr, 'o3loss', p3=o3loss &
                         , reprid=repr_idx                           )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr,'o3loss' &
           ,'longname' ,c='Ozone loss from tracer')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'o3loss'  &
                      ,'units' ,c='mol/mol/s')
      CALL channel_halt(substr, status)
      IF (p_parallel_io) WRITE(*,*) '.... ozone loss from tracer added to channel ',modstr

      !
      !O3 prod
      !

      CALL new_channel_object(status, modstr, 'o3Prod', p3=o3prod &
                         , reprid=repr_idx                           )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'o3Prod' &
                        ,'longname' ,c='Ozone prod')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'o3Prod'  &
                         ,'units' ,c='mol/mol/s')
      CALL channel_halt(substr, status)
      IF (p_parallel_io) WRITE(*,*) '.... ozone prod added to channel ',modstr

      !
      !CH4loss
      !
      CALL new_channel_object(status, modstr, 'ch4loss', p3=ch4loss &
                              ,reprid=repr_idx                           )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ch4loss' &
                         ,'longname' ,c='CH4 loss')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ch4loss'  &
                      ,'units' ,c='mol/mol/s')
      CALL channel_halt(substr, status)
      IF (p_parallel_io) WRITE(*,*) '.... CH4 loss  added to channel ',modstr

      !
      !CO Prod
      !
      CALL new_channel_object(status, modstr, 'coprod', p3=coprod &
                             , reprid=repr_idx                           )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'coprod' &
                         ,'longname' ,c='CO Prod')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'coprod'  &
                         ,'units' ,c='mol/mol/s')
      CALL channel_halt(substr, status)
      IF (p_parallel_io) WRITE(*,*) '.... CO Prod  added to channel ',modstr

      !
      !CO loss
      !
      CALL new_channel_object(status, modstr, 'coloss', p3=coloss &
                            , reprid=repr_idx                           )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'coloss' &
                         ,'longname' ,c='CO Loss')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'coloss'  &
                         ,'units' ,c='mol/mol/s')
      CALL channel_halt(substr, status)
      IF (p_parallel_io) WRITE(*,*) '.... CO loss  added to channel ',modstr


      !
      !now steady state output for OH and HO2
      !
      
      ! loop over oh
   
!op_mm_20171016   
  !    setname=  (/ 'OHlig    ','OHbio    ','OHsoi    ','OHind    ','OHtra    ',&
  !                 'OHshp    ','OHair    ','OHn2o    ','OHch4    ','OHstr    '&
  !                 /)

  !    longname= (/ 'OH from lightning       ','OH from biomass burning ',&
  !                 'OH from soils           ','OH from industry        ',&
  !                 'OH from road traffic    ','OH from ships           ',&
  !                 'OH from air traffic     ','OH from N2O+ UB         ',&
  !                 'OH from CH4             ','OH from Stratosphere    '&
  !                 /)



      do i=1, N_tag_species
!op_mm_20171016   
!         oh_tag(i)%name=setname(i)
!         oh_tag(i)%longname=longname(i)

        oh_tag(i)%name=trim("OH") // trim(CATEGORY_SET(i)%abr)
        oh_tag(i)%longname=trim("OH of ") // trim(CATEGORY_SET(i)%name)

   
         CALL new_channel_object(status, modstr, TRIM(oh_tag(i)%name), p3=oh_tag(i)%ptr &
                                 , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(oh_tag(i)%name) &
                             ,'longname' ,c=TRIM(oh_tag(i)%longname))
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(oh_tag(i)%name)  &
                            ,'units' ,c='mol/mol')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',TRIM(oh_tag(i)%name),'from  added to channel ',modstr

      end do

     !
     ! loop over ho2
     !
!op_mm_20171016   
  !    setname=  (/     'HO2lig   ','HO2bio   ','HO2soi   ','HO2ind   ','HO2tra   ',&
  !                     'HO2shp   ','HO2air   ','HO2n2o   ','HO2ch4   ','HO2str   '&
  !                     /)

   !   longname= (/     'HO2 from lightning      ','HO2 from biomass burning',&
   !                    'HO2 from soils          ','HO2 from industry       ',&
   !                    'HO2 from road traffic   ','HO2 from shipping       ',&
   !                    'HO2 from air            ','HO2 from N2O + UB       ',&
   !                    'HO2 from CH4            ','HO2 from STR            '&
    !             /)

 
      do i=1, N_tag_species
!op_mm_20171016   
!         ho2_tag(i)%name=setname(i)
!         ho2_tag(i)%longname=longname(i)
 
          ho2_tag(i)%name=trim("HO2") // trim(CATEGORY_SET(i)%abr)
          ho2_tag(i)%longname=trim("HO2 of") // trim(CATEGORY_SET(i)%name)



 
         CALL new_channel_object(status, modstr, TRIM(ho2_tag(i)%name), p3=ho2_tag(i)%ptr &
                                , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(ho2_tag(i)%name) &
                            ,'longname' ,c=TRIM(ho2_tag(i)%longname))
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(ho2_tag(i)%name)  &
                           ,'units' ,c='mol/mol')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',ho2_tag(i)%name,'from  added to channel ',modstr

      end do

!
!now initialise the channel objects for tendencies
!

!       
! loop over p_ho2
!
!op_mm_20171016   
!      setname=      (/                                                            &  
!                      'PO3ho2lig','PO3ho2bio','PO3ho2soi','PO3ho2ind','PO3ho2tra',&
!                     'PO3ho2shp','PO3ho2air','PO3HO2N2O','PO3HO2CH4','PO3HO2str'&
!                      /)


      do i=1, N_tag_species
!op_mm_20171016   
!         p_ho2_ten(i)%name=setname(i)
 
          p_ho2_ten(i)%name=trim("PO3HO2") // trim(CATEGORY_SET(i)%abr)

 
 
         CALL new_channel_object(status, modstr, TRIM(p_ho2_ten(i)%name), p3=p_ho2_ten(i)%ptr &
                                ,reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(p_ho2_ten(i)%name) &
                            ,'longname' ,c='PO3 HO2 rate')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(p_ho2_ten(i)%name)  &
                            ,'units' ,c='mol/mol/s')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',p_ho2_ten(i)%name,'from  added to channel ',modstr

      end do

!
! loop over p_ro2
!
!op_mm_20171016   
!      setname=      (/                                                            &  
!                      'PO3RO2lig','PO3RO2bio','PO3RO2soi','PO3RO2ind','PO3RO2tra',&
!                      'PO3RO2shp','PO3RO2air','PO3RO2N2O','PO3RO2CH4','PO3RO2str'&
!                      /)


      do i=1, N_tag_species
!op_mm_20171016   
!         p_ro2_ten(i)%name=setname(i)
 
        p_ro2_ten(i)%name=trim("PO3RO2") // trim(CATEGORY_SET(i)%abr)

 
        CALL new_channel_object(status, modstr, TRIM(p_ro2_ten(i)%name), p3=p_ro2_ten(i)%ptr &
                               , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(p_ro2_ten(i)%name) &
                           ,'longname' ,c='PO3 RO2 rate')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(p_ro2_ten(i)%name)  &
                           ,'units' ,c='mol/mol/s')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',p_ro2_ten(i)%name,'from  added to channel ',modstr

      end do



!
! loop over DO3 OH
!
!op_mm_20171016   
!      setname=      (/                                                            &  
!                      'DO3OHlig ','DO3OHbio ','DO3OHsoi ','DO3OHind ','DO3OHtra ',&
!                       'DO3OHshp ','DO3OHair ','DO3OHN2O ','DO3OHCH4 ','DO3OHstr '&
!                       /)


       


      do i=1, N_tag_species
!op_mm_20171016            
!         d_oh_ten(i)%name=setname(i)

        d_oh_ten(i)%name=trim("DO3OH") // trim(CATEGORY_SET(i)%abr)  


         CALL new_channel_object(status, modstr, TRIM(d_oh_ten(i)%name), p3=d_oh_ten(i)%ptr &
                                , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_oh_ten(i)%name) &
                           ,'longname' ,c='DO3 OH rate')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_oh_ten(i)%name)  &
                           ,'units' ,c='mol/mol/s')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',d_oh_ten(i)%name,'from  added to channel ',modstr

      end do


! loop over DO3 HO2
!op_mm_20171016   
!      setname=      (/                                                            &  
!                      'DO3HO2lig','DO3HO2bio','DO3HO2soi','DO3HO2ind','DO3HO2tra',&
!                       'DO3HO2shp','DO3HO2air','DO3HO2NO2','DO3HO2CH4','DO3HO2str'&
!                       /)

        

      do i=1, N_tag_species
!op_mm_20171016   
!         d_ho2_ten(i)%name=setname(i)
         d_ho2_ten(i)%name=trim("DO3HO2") // trim(CATEGORY_SET(i)%abr) 



         CALL new_channel_object(status, modstr, TRIM(d_ho2_ten(i)%name), p3=d_ho2_ten(i)%ptr &
                               , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_ho2_ten(i)%name) &
                           ,'longname' ,c='DO3 OH rate')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_ho2_ten(i)%name)  &
                            ,'units' ,c='mol/mol/s')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',d_ho2_ten(i)%name,'from  added to channel ',modstr

      end do

!
! loop over DO3 NO
!
!op_mm_20171016   
!      setname=      (/                                                            &  
!                      'DO3NOlig ','DO3NObio ','DO3NOsoi ','DO3NOind ','DO3NOtra ',&
!                       'DO3NOshp ','DO3NOair ','DO3NON2O ','DO3NOCH4 ','DO3NOstr '&
!                       /)



      do i=1, N_tag_species
!op_mm_20171016   
!         d_no_ten(i)%name=setname(i)
         d_no_ten(i)%name=trim("DO3NO") // trim(CATEGORY_SET(i)%abr)  
   


         CALL new_channel_object(status, modstr, TRIM(d_no_ten(i)%name), p3=d_no_ten(i)%ptr &
                                , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_no_ten(i)%name) &
                          ,'longname' ,c='DO3 NO rate')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_no_ten(i)%name)  &
                           ,'units' ,c='mol/mol/s')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',d_no_ten(i)%name,'from  added to channel ',modstr

      end do

!
! loop over DO3 RO
!
!op_mm_20171016   
!      setname=      (/                                                            &  
!                     'DO3ROlig ','DO3RObio ','DO3ROsoi ','DO3ROind ','DO3ROtra ',&
!                     'DO3ROshp ','DO3ROair ','DO3RON2O ','DO3ROCH4 ','DO3ROstr '&
!                     /)

                      

      do i=1, N_tag_species

!         d_ro_ten(i)%name=setname(i)
         d_ro_ten(i)%name=trim("DO3RO") // trim(CATEGORY_SET(i)%abr)  
   


         CALL new_channel_object(status, modstr, TRIM(d_ro_ten(i)%name), p3=d_ro_ten(i)%ptr &
                                , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_ro_ten(i)%name) &
                           ,'longname' ,c='DO3 RO rate')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_ro_ten(i)%name)  &
                           ,'units' ,c='mol/mol/s')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',d_ro_ten(i)%name,'from  added to channel ',modstr

      end do



! loop over DO3 XO
!op_mm_20171016   
!      setname=      (/                                                            &  
!                     'DO3XOlig ','DO3XObio ','DO3XOsoi ','DO3XOind ','DO3XOtra ',&
!                     'DO3XOshp ','DO3XOair ','DO3XON2O ','DO3XOCH4 ','DO3XOstr ' & 
!                     /)

  
         


      do i=1, N_tag_species
!op_mm_20171016   
!         d_xo_ten(i)%name=setname(i)
         d_xo_ten(i)%name=trim("DO3XO") // trim(CATEGORY_SET(i)%abr) 
  

         CALL new_channel_object(status, modstr, TRIM(d_xo_ten(i)%name), p3=d_xo_ten(i)%ptr &
                                , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_xo_ten(i)%name) &
                           ,'longname' ,c='DO3 XO rate')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(d_xo_ten(i)%name)  &
                            ,'units' ,c='mol/mol/s')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',d_xo_ten(i)%name,'from  added to channel ',modstr

      end do



!op_vr_20170215+
!op_mm_20171016   
!     setname=      (/                                                            &  
!                     'ch4losslig ','ch4lossbio ','ch4losssoi ','ch4lossind ','ch4losstra ',&
!                     'ch4lossshp ','ch4lossair ','ch4lossN2O ','ch4lossCH4 ','ch4lossstr ' & 
!                     /)

!      longname= (/     'ch4loss from lightning      ','ch4loss from biomass burning',&
!                       'ch4loss from soils          ','ch4loss from industry       ',&
!                       'ch4loss from road traffic   ','ch4loss from shipping       ',&
!                       'ch4loss from air            ','ch4loss from N2O + UB       ',&
!                       'ch4loss from CH4            ','ch4loss from STR            '&
!                 /)



      do i=1, N_tag_species
!op_mm_20171016   
!         ch4loss_tag(i)%name=setname(i)
!         ch4loss_tag(i)%longname=longname(i)
 
          ch4loss_tag(i)%name=trim("ch4loss") // trim(CATEGORY_SET(i)%abr)
          ch4loss_tag(i)%longname=trim("CH4loss from ") // trim(CATEGORY_SET(i)%name)

 
         CALL new_channel_object(status, modstr, TRIM(ch4loss_tag(i)%name), p3=ch4loss_tag(i)%ptr &
                                , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(ch4loss_tag(i)%name) &
                           ,'longname' ,c='tagged CH4 loss')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(ch4loss_tag(i)%name)  &
                            ,'units' ,c='mol/mol/s')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',ch4loss_tag(i)%name,'from  added to channel ',modstr

      end do



!op_mm_20171016
  !   setname=      (/                                                            &  
  !                   'fracODDlig ','fracODDbio ','fracODDsoi ','fracODDind ','fracODDtra ',&
  !                   'fracODDshp ','fracODDair ','fracODDN2O ','fracODDCH4 ','fracODDstr ' & 
  !                   /)

  !    longname= (/     'fracODD from lightning      ','fracODD from biomass burning',&
  !                     'fracODD from soils          ','fracODD from industry       ',&
  !                     'fracODD from road traffic   ','fracODD from shipping       ',&
  !                     'fracODD from air            ','fracODD from N2O + UB       ',&
  !                     'fracODD from CH4            ','fracODD from STR            '&
   !              /)



      do i=1, N_tag_species
!op_mm_20171016
!         frac_ODD(i)%name=setname(i)
!         frac_ODD(i)%longname=longname(i)

          frac_ODD(i)%name=trim("fracODD") // trim(CATEGORY_SET(i)%abr)
          frac_ODD(i)%longname=trim("fracODD from ") // trim(CATEGORY_SET(i)%name)


         CALL new_channel_object(status, modstr, TRIM(frac_ODD(i)%name), p3=frac_ODD(i)%ptr &
                                , reprid=repr_idx                           )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(frac_ODD(i)%name) &
                           ,'longname' ,c='fraction of tagged ODD to total ODDD')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr, TRIM(frac_ODD(i)%name)  &
                            ,'units' ,c='')
         CALL channel_halt(substr, status)
         IF (p_parallel_io) WRITE(*,*) ' ...',frac_ODD(i)%name,'from  added to channel ',modstr

      end do

!op_vr_20170215-





      CALL end_message_bi(modstr, 'CHANNEL  DEFINITION', substr)

    END SUBROUTINE tagging_init_memory
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------

  SUBROUTINE tagging_init_coupling
 
    ! ECHAM5
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_bcast 
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer,           ONLY: get_tracer
    USE messy_main_channel,          ONLY: get_channel_object
    USE messy_main_channel_error_BI, ONLY: channel_halt 

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr = 'tagging_init_coupling'
    INTEGER                       :: ierr    ! error flag
    INTEGER                       :: ispec   !


    integer                       :: status


    ! PROCESS NAMELIST

       CALL start_message_bi(modstr, 'COUPLING INITIALIZING', substr)

       ierr=0

!++++++++++++++++++++++++++++++
! Get relevant chemical Tracers
!+++++++++++++++++++++++++++++
! Get NOy Family
       do ispec=1,Fam_NOy_species
          CALL get_tracer(ierr, GPTRSTR, TRIM(NOy_species(ispec)),idx=idx_NOy(ispec))
          if (ierr/=0) CALL error_bi('Tracer in NOy Family not found:'//TRIM(NOy_species(ispec)), substr)
          IF (p_parallel_io) WRITE(*,*) ' ... fechted NOy tracer ',ispec,  &
                             TRIM(NOy_species(ispec)),' from ',GPTRSTR
       enddo

! Get NMHC Family
       do ispec=1,Fam_NMHC_species
         CALL get_tracer(ierr, GPTRSTR, TRIM(NMHC_species(ispec)),idx=idx_NMHC(ispec))
         if (ierr/=0) CALL error_bi('Tracer in NMHC-family not found:'//TRIM(NMHC_species(ispec)), substr)
         IF (p_parallel_io) WRITE(*,*) ' ... fechted NMHC tracer ',ispec,  &
                             TRIM(NMHC_species(ispec)),' from ',GPTRSTR
       enddo

! Get Ozone
       CALL get_tracer(ierr, GPTRSTR, 'O3',idx=idx_O3)
       if (ierr/=0) CALL error_bi('Tracer O3 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted O3 from ',GPTRSTR,'with idx=',idx_O3

       CALL get_tracer(ierr, GPTRSTR, 'O1D',idx=idx_O1D)
       if (ierr/=0) CALL error_bi('Tracer O1D not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted O1D from ',GPTRSTR,'with idx=',idx_O1D

! Get PAN
       CALL get_tracer(ierr, GPTRSTR, 'PAN',idx=idx_PAN)
       if (ierr/=0) CALL error_bi('Tracer PAN not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted PAN from ',GPTRSTR

! Get CO 
       CALL get_tracer(ierr, GPTRSTR, 'CO',idx=idx_CO)
       if (ierr/=0) CALL error_bi('Tracer CO not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted CO from ',GPTRSTR

! Get CH4        
       CALL get_tracer(ierr, GPTRSTR, 'CH4',idx=idx_CH4)
       if (ierr/=0) CALL error_bi('Tracer CH4 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted CH4 from ',GPTRSTR

! Get OH
       CALL get_tracer(ierr, GPTRSTR, 'OH',idx=idx_oh)
       if (ierr/=0) CALL error_bi('Tracer OH not found ', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted OH from ',GPTRSTR

! G1et HO2
       CALL get_tracer(ierr, GPTRSTR, 'HO2',idx=idx_ho2)
       if (ierr/=0) CALL error_bi('Tracer HO2 not found ', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted HO2 from  ',GPTRSTR



! Get diagnostic prod and loss rates  (loss = loss temporary, will be revised)
! ozone loss and prod
       CALL get_tracer(ierr, GPTRSTR, 'o3lossho2',idx=idx_o3loss_ho2) 
       if (ierr/=0) CALL error_bi('Tracer LO3_HO2 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted LO3_HO2'

       CALL get_tracer(ierr, GPTRSTR, 'o3lossoh',idx=idx_o3loss_oh)
       if (ierr/=0) CALL error_bi('Tracer LO3_OH not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted LO3_OH'

       CALL get_tracer(ierr, GPTRSTR, 'o3lossro',idx=idx_o3loss_ro)
       if (ierr/=0) CALL error_bi('Tracer LO3_RO not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted LO3_RO'

       CALL get_tracer(ierr, GPTRSTR, 'o3lossno',idx=idx_o3Loss_no)
       if (ierr/=0) CALL error_bi('Tracer LO3_NO not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted LO3_NO'

       CALL get_tracer(ierr, GPTRSTR, 'o3lossxo',idx=idx_o3loss_xo)
       if (ierr/=0) CALL error_bi('Tracer LO3_XO not found',  substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted LO3_XO'

       CALL get_tracer(ierr, GPTRSTR, 'o3prodro2',idx=idx_o3prod_ro2)
       if (ierr/=0) CALL error_bi('Tracer PO3_RO2 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted PO3_RO2'

       CALL get_tracer(ierr, GPTRSTR, 'o3prodho2',idx=idx_o3prod_ho2)     !ProdHO2 
       if (ierr/=0) CALL error_bi('Tracer PO3_HO2 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted PO3_HO2'
 

       CALL get_tracer(ierr, GPTRSTR, 'o3prodo2',idx=idx_o3prod_o2)
       if (ierr/=0) CALL error_bi('Tracer PO3_O2 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted PO3_O2'

       CALL get_tracer(ierr, GPTRSTR, 'ProdMeO2',idx=idx_o3prod_MeO2)
       if (ierr/=0) CALL error_bi('Tracer ProdMeO2 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted ProdMeO2' 

       CALL get_tracer(ierr, GPTRSTR, 'noprodN2O',idx=idx_noprod_N2O)
       if (ierr/=0) CALL error_bi('Tracer noprodN2O not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted noprodN2O'


! Get  HCs 
! PAN
       CALL get_tracer(ierr, GPTRSTR, 'panprod',idx=idx_panprod)
       if (ierr/=0) CALL error_bi('Tracer PPAN not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted PPAN'

       CALL get_tracer(ierr, GPTRSTR, 'panloss',idx=idx_panloss)
       if (ierr/=0) CALL error_bi('Tracer LPAN not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... fechted LPAN'

! CO
      ! calculated; not from MECCA Diag
! CH4
      ! calculated; not from MECCA Diag


! HOx
       ! get rate (=Loss) of H + O2 -> HO2
       CALL get_tracer(ierr, GPTRSTR, 'LossG2100',idx=idx_LossG2100)
       if (ierr/=0) CALL error_bi('Tracer LossG2100 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2100'

       ! get rate (=Loss) of OH + O3P -> H + O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG2103',idx=idx_LossG2103)
       if (ierr/=0) CALL error_bi('Tracer LossG2103 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2103'

       ! get rate (=Loss) of OH + O3 -> HO2 + O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG2104',idx=idx_LossG2104)
       if (ierr/=0) CALL error_bi('Tracer LossG2104 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2104'

       ! get rate (=Loss) of OH + H2 -> H2O + H
       CALL get_tracer(ierr, GPTRSTR, 'LossG2105',idx=idx_LossG2105)
       if (ierr/=0) CALL error_bi('Tracer LossG2105 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2105'

       ! get rate (=Loss) of HO2 + O3P -> OH + O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG2106',idx=idx_LossG2106)
       if (ierr/=0) CALL error_bi('Tracer LossG2106 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2106'

       ! get rate (=Loss) of HO2 + O3 -> OH + 2 O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG2107',idx=idx_LossG2107)
       if (ierr/=0) CALL error_bi('Tracer LossG2107 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2107'

       ! get rate (=Loss) of HO2 + OH -> H2O + O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG2109',idx=idx_LossG2109)
       if (ierr/=0) CALL error_bi('Tracer LossG2109 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2109'

       ! get rate (=Loss) of HO2 + HO2 -> H2O2 + O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG2110',idx=idx_LossG2110)
       if (ierr/=0) CALL error_bi('Tracer LossG2110 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2110'

       ! get rate (=Loss) of H2O + O1D -> 2 OH
       CALL get_tracer(ierr, GPTRSTR, 'LossG2111',idx=idx_LossG2111)
       if (ierr/=0) CALL error_bi('Tracer LossG2111 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2111'

       ! get rate (=Loss) of H2O2 + OH -> H2O + HO2
       CALL get_tracer(ierr, GPTRSTR, 'LossG2112',idx=idx_LossG2112)
       if (ierr/=0) CALL error_bi('Tracer LossG2112 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG2112'

       ! get rate (=Loss) of NO + OH -> HONO
       CALL get_tracer(ierr, GPTRSTR, 'LossG3200',idx=idx_LossG3200)
       if (ierr/=0) CALL error_bi('Tracer LossG3200 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG3200'

       ! get rate (=Loss) of NO + HO2 -> NO2 + OH
       CALL get_tracer(ierr, GPTRSTR, 'LossG3201',idx=idx_LossG3201)
       if (ierr/=0) CALL error_bi('Tracer LossG3201 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG3201'

       ! get rate (=Loss) of NO2 + OH -> HNO3
       CALL get_tracer(ierr, GPTRSTR, 'LossG3202',idx=idx_LossG3202)
       if (ierr/=0) CALL error_bi('Tracer LossG3202 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG3202'

       ! get rate (=Loss) of NO2 + HO2 -> HNO4
       CALL get_tracer(ierr, GPTRSTR, 'LossG3203',idx=idx_LossG3203)
       if (ierr/=0) CALL error_bi('Tracer LossG3203 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG3203'

       ! get rate (=Loss) of HNO4 -> NO2 + HO2
       CALL get_tracer(ierr, GPTRSTR, 'LossG3207',idx=idx_LossG3207)
       if (ierr/=0) CALL error_bi('Tracer LossG3207 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG3207'

       ! get rate (=Loss) of CH4 + OH -> CH3O2 + H2O
       CALL get_tracer(ierr, GPTRSTR, 'LossG4101',idx=idx_LossG4101)
       if (ierr/=0) CALL error_bi('Tracer LossG4101 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG4101'

       ! get rate (=Loss) of CO + OH -> H + CO2
       CALL get_tracer(ierr, GPTRSTR, 'LossG4110',idx=idx_LossG4110)
       if (ierr/=0) CALL error_bi('Tracer LossG4110 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG4110'

       ! get rate (=Loss) of ClO + OH -> .94 Cl + .94 HO2 + .06 HCl + .06 O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG6203',idx=idx_LossG6203)
       if (ierr/=0) CALL error_bi('Tracer LossG6203 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG6203'

       ! get rate (=Loss) of ClO + HO2 -> HOCl + O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG6204',idx=idx_LossG6204)
       if (ierr/=0) CALL error_bi('Tracer LossG6204 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG6204'

       ! get rate (=Loss) of BrO + HO2 -> HOBr + O2
       CALL get_tracer(ierr, GPTRSTR, 'LossG7201',idx=idx_LossG7201)
       if (ierr/=0) CALL error_bi('Tracer LossG7201 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossG7201'

       ! get rate (=Loss) of HONO + hv -> NO + OH
       CALL get_tracer(ierr, GPTRSTR, 'LossJ3200',idx=idx_LossJ3200)
       if (ierr/=0) CALL error_bi('Tracer LossJ3200 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossJ3200'

       ! get rate (=Loss) of HNO3 + hv -> NO2 + OH
       CALL get_tracer(ierr, GPTRSTR, 'LossJ3201',idx=idx_LossJ3201)
       if (ierr/=0) CALL error_bi('Tracer LossJ3201 not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossJ3201'

       ! get rate (=Loss) of HCHO + hv -> H + CO + HO2
       CALL get_tracer(ierr, GPTRSTR, 'LossJ4101b',idx=idx_LossJ4101b)
       if (ierr/=0) CALL error_bi('Tracer LossJ4101b not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... LossJ4101b'

       ! get rate (=Loss) of NMHC + OH -> NMHC
       CALL get_tracer(ierr, GPTRSTR, 'OHlossNMHC',idx=idx_OHlossNMHC)
       if (ierr/=0) CALL error_bi('Tracer OHlossNMHC not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... OHlossNMHC'

       ! get rate (=Loss) of NMHC + HO2 -> NMHC
       CALL get_tracer(ierr, GPTRSTR, 'HO2lossNMHC',idx=idx_HO2lossNMHC)
       if (ierr/=0) CALL error_bi('Tracer HO2lossNMHC not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... HO2lossNMHC'

       ! get rate (=Loss) of total NMHC + OH -> NMHC + HO2
       CALL get_tracer(ierr, GPTRSTR, 'OHlossHO2prodNMHC',idx=idx_OHlossHO2prodNMHC)
       if (ierr/=0) CALL error_bi('Tracer OHlossHO2prodNMHC not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... OHlossHO2prodNMHC'

       ! get rate (=Loss) of total NMHC + NOy -> HO2 + NMHC + NOy
       CALL get_tracer(ierr, GPTRSTR, 'HO2prodNMHCNOy',idx=idx_HO2prodNMHCNOy)
       if (ierr/=0) CALL error_bi('Tracer HO2prodNMHCNOy not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... HO2prodNMHCNOy'

       ! get rate (=Loss) of NMHC + hv -> NMHC + HO2
       CALL get_tracer(ierr, GPTRSTR, 'HO2prodNMHCphoto',idx=idx_HO2prodNMHCphoto)
       if (ierr/=0) CALL error_bi('Tracer HO2prodNMHCphoto not found', substr)
       IF (p_parallel_io) WRITE(*,*) ' ... HO2prodNMHCphoto'


     
   ! op_mm_20131126
    ! now using climatological tropopause calculated online before first timestep
       ! CALL get_channel_object(status,trim(c_tropop(1)),trim(c_tropop(2)),p2=tp_p_wmo) 
    ! if (ierr/=0) CALL error_bi('Channel not found '//c_tropop(1), substr)
    ! CALL channel_halt(substr, status)  
    ! IF (p_parallel_io) WRITE(*,*) ' ... fechted Tropopause'
     

     

!Get Lightning Emissions

      !old calls
      ! CALL get_stream(lnox,c_lnox(1),ierr)
      ! if (ierr/=0) CALL FINISH(substr,'Stream not found '//c_lnox(1))
      ! CALL get_stream_element(lnox,c_lnox(2),telnox,ierr)
      ! if (ierr/=0) CALL FINISH(substr,'Stream element not found '//c_lnox(2))

       CALL get_channel_object(status,trim(c_lnox(1)),trim(c_lnox(2)),p3=telnox) 
       if (ierr/=0) CALL error_bi('Channel not found '//c_lnox(1), substr)
       CALL channel_halt(substr//' getting lnox', status)  
       IF (p_parallel_io) WRITE(*,*) ' ... fechted Lightning emissions'


       CALL end_message_bi(modstr, 'COUPLING INITIALIZING', substr)

  END SUBROUTINE tagging_init_coupling


!=============================================================================
SUBROUTINE tagging_global_start
!-------------------------------------------------------------------------
! global start routine
!-------------------------------------------------------------------------
  
  USE messy_main_timer,             ONLY: time_step_len, lstart, lresume
  USE messy_main_tracer_mem_bi, ONLY:  xtte


  IMPLICIT none
  CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_global_start'


  if (force_reinitialize_tracer) then 
     call tagging_init_tracer_online 
     ltagging=.true.  
     force_reinitialize_tracer=.false. 
  end if


END SUBROUTINE tagging_global_start



!=============================================================================


SUBROUTINE tagging_local_start
!-------------------------------------------------------------------------
! This subroutine stores the  chemical fields after transport into X0 fields
!    - or should only the tendencies stored?
! And calculates the difference in numerical diffusion due to advection 
!    for the tracer and the sum of the taged tracers 
! And finally scales the taged traces if desired 
!
! Author, V. Grewe, DLR-OP and NCAR Boulder August 2008
!-------------------------------------------------------------------------

  ! ECHAM5/MESSy
  USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1=> qxtm1


  USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma, nproma 
                                    

  !changed: timestep -> messy_timer
  USE messy_main_timer,         ONLY: time_step_len, lstart
  ! MESSy
  USE messy_tagging,            ONLY: i_advect_scaling


  IMPLICIT none
 
 ! REAL(dp), DIMENSION(nproma,nlev) :: H2O_tot, HNO3_tot, Tice, Tnat
 ! INTEGER,  DIMENSION(nproma)      :: index_arr ! mz_bs_20060112
      
 !et to logical
 !LOGICAL,  PARAMETER             :: linit=.TRUE.
  CHARACTER(LEN=*), PARAMETER      :: substr = 'tagging_local_start'
  REAL (DP)                        :: inv_time_step_len  !inverse of time step length ! op_mm_20140124
 
  INTEGER :: ispec, idt, jp,jk

  INTRINSIC MAX, TINY, MINLOC, SUM, MINVAL, MAXVAL, SIZE
 

   !write (*,*) 'ltagging_local1', ltagging

  ! op_mm_20131024+
  !moved init part to tagging global start because it should not be done within the the jrow loop
  ! it should be done within the time loop!  
  if (.not.ltagging) then
     print*,"tagging switched off so far"
     return  
  end if
  ! op_mm_20131024-




!---------
! Set species: X0 values: ozone , carbon monoxide and PAN 
!---------


  idt=idx_O3
  o3_0(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) & 
       + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
   
  idt=idx_CO
  co_0(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &
       + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

  idt=idx_pan
  pan_0(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &
       + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len 
  
    !    oh_0(1:kproma,:,jrow) = pxtm1(1:kproma,:,idx_ohe) &
    !                        + pxtte(1:kproma,:,idx_ohe) * time_step_len 
    !   ho2_0(1:kproma,:,jrow) = pxtm1(1:kproma,:,idx_ho2e) &
    !                        + pxtte(1:kproma,:,idx_ho2e) * time_step_len 


!---------
! Set families NOy and NMHC
!---------  
        
  noy_0(1:kproma,:,jrow) = .0_dp
  do ispec=1,Fam_NOy_species
     idt=idx_NOy(ispec)
     noy_0(1:kproma,:,jrow) = noy_0(1:kproma,:,jrow) + N_in_Fam(ispec) *   &
                            (   pxtm1(_RI_X_ZN_(1:kproma,:,idt))                   &
                             + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len )
  enddo

  nmhc_0(1:kproma,:,jrow) = 0.0_dp
  do ispec=1,Fam_NMHC_species
     idt=idx_NMHC(ispec)
     nmhc_0(1:kproma,:,jrow) = nmhc_0(1:kproma,:,jrow) + C_in_Fam(ispec) * &
                              (pxtm1(_RI_X_ZN_(1:kproma,:,idt))                     &
                              + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len )
  enddo


      
!---------
! Set diagnostic tagged species
!---------


  o3_tagging_sum(1:kproma,:,jrow)=0.0_dp
  noy_tagging_sum(1:kproma,:,jrow)=0.0_dp
  nmhc_tagging_sum(1:kproma,:,jrow)=0.0_dp
  pan_tagging_sum(1:kproma,:,jrow)=0.0_dp
  co_tagging_sum(1:kproma,:, jrow)=0.0_dp

       

  do ispec=1,N_tag_species

     !+++++++
     !O3
     !++++++
     idt=O3_tag(ispec)%idx      
     o3_tagging(1:kproma,:,ispec,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))   &
                                         + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
     o3_tagging_sum(1:kproma,:,jrow) = o3_tagging_sum(1:kproma,:,jrow)+o3_tagging(1:kproma,:,ispec,jrow)

     !+++++++
     !NOy
     !++++++
     idt=noy_tag(ispec)%idx
     noy_tagging(1:kproma,:,ispec,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))   &
                                         + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
     noy_tagging_sum(1:kproma,:,jrow) = noy_tagging_sum(1:kproma,:,jrow)+noy_tagging(1:kproma,:,ispec,jrow)

     !+++++++
     !NMHC
     !++++++
     idt=nmhc_tag(ispec)%idx
     nmhc_tagging(1:kproma,:,ispec,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))   &
                                              + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
     nmhc_tagging_sum(1:kproma,:,jrow) = nmhc_tagging_sum(1:kproma,:,jrow)+nmhc_tagging(1:kproma,:, ispec,jrow)

     !+++++++
     !PAN
     !++++++
     idt=pan_tag(ispec)%idx
     pan_tagging(1:kproma,:,ispec,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))   &
                                             + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
     pan_tagging_sum(1:kproma,:,jrow)= pan_tagging_sum(1:kproma,:,jrow)+pan_tagging(1:kproma,:,ispec,jrow)

     !+++++++
     !CO
     !++++++
     idt=co_tag(ispec)%idx
     co_tagging(1:kproma,:,ispec,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &  
                                            + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
     co_tagging_sum(1:kproma,:,jrow) = co_tagging_sum(1:kproma,:,jrow)+co_tagging(1:kproma,:,ispec,jrow)

  enddo


!---------------------------------------------------------------------
!---------
! Account for tendencies, which are not from chemistry or transport
!  xnew = xold * fac => tendency=(xold*fac-xold)/dt  =  xold*(fac-1)/dt  fac=sum/tag_sum
!---------


      SELECT CASE (I_ADVECT_SCALING)  ! do scaling
          CASE(1)
           IF (l_adv_err_diag) then
              do3_err(1:kproma,:,jrow)=o3_0(1:kproma,:,jrow)-o3_tagging_sum(1:kproma,:,jrow)
              dnoy_err(1:kproma,:,jrow)=noy_0(1:kproma,:,jrow)-noy_tagging_sum(1:kproma,:,jrow)
              dco_err(1:kproma,:,jrow)=co_0(1:kproma,:,jrow)-co_tagging_sum(1:kproma,:,jrow)
              dnmhc_err(1:kproma,:,jrow)=nmhc_0(1:kproma,:,jrow)-nmhc_tagging_sum(1:kproma,:,jrow)
              dpan_err(1:kproma,:,jrow)=pan_0(1:kproma,:,jrow)-pan_tagging_sum(1:kproma,:,jrow)   
              
              ! op_mm_20140123 Added dnmhc_err and dpan_err
           endif

 ! op_mm_20131125
!roll out do-loops to overcome divisions by zero. Advection correction is only done if the tagging sum values larger as a certain threshold. 
!multiply with inverse of time_step_len

           inv_time_step_len=1.0_dp/time_step_len

           do ispec=1,N_tag_species
              !!
              !! O3
              !!
              do jk=1, nlev
                 do jp=1, kproma  
                    idt=o3_tag(ispec)%idx
!                   if (o3_tagging_sum(i,k, jrow) .eq. 0.0_dp) then
                    if (abs(o3_tagging_sum(jp,jk, jrow)) .ge. 1.0e-20_dp) then
                       pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                        &
                                         +o3_tagging(jp,jk,ispec,jrow)                              &
                                         *(o3_0(jp,jk,jrow)/o3_tagging_sum(jp,jk,jrow)-1._dp)       & 
                                         *inv_time_step_len
                      !                   / time_step_len
                    else
                     ! do nothing 
                    end if
                 end do
              end do


              
              !!
              !! NOY
              !!
              do jk=1, nlev
                 do jp=1, kproma
                    idt=noy_tag(ispec)%idx
!                    if (noy_tagging_sum(i,k, jrow) .eq. 0.0_dp) then
                    if (abs(noy_tagging_sum(jp,jk, jrow)) .ge. 1.0e-20) then
                       pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                          &
                                         +noy_tagging(jp,jk,ispec,jrow)                               &
                                         *(noy_0(jp,jk,jrow)/noy_tagging_sum(jp,jk,jrow)-1._dp)       & 
                                         *inv_time_step_len
                       !                  / time_step_len
                    else
                       ! do nothing 
                    end if
                 end do
              end do
             
 
              !!
              !! CO
              !!             

              do jk=1, nlev
                 do jp=1, kproma
                    idt=co_tag(ispec)%idx
!                    if (co_tagging_sum(i,k, jrow) .eq. 0.0_dp) then
                    if (abs(co_tagging_sum(jp,jk, jrow)) .ge. 1.0e-20) then
                       pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                           &
                                         +co_tagging(jp,jk,ispec,jrow)                                 &
                                         * (co_0(jp,jk,jrow)/co_tagging_sum(jp,jk,jrow)-1._dp)         & 
                                         *inv_time_step_len
                       !                  / time_step_len
                    else
                       ! do nothing 
                    end if
                 end do
              end do

              
              !!
              !!NMHC
              !!

              do jk=1, nlev
                 do jp=1, kproma
                    idt=nmhc_tag(ispec)%idx
!                    if (nmhc_tagging_sum(i,k, jrow) .eq. 0.0_dp) then
                    if (abs(nmhc_tagging_sum(jp,jk, jrow)) .ge. 1.0e-20) then
                       pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                             &
                                         +nmhc_tagging(jp,jk,ispec,jrow)                                 &
                                         * (nmhc_0(jp,jk,jrow)/nmhc_tagging_sum(jp,jk,jrow)-1._dp)       & 
                                         *inv_time_step_len
                      !                  / time_step_len
                   else
                        ! do nothing
                   end if
                end do
             end do
     
         

               !!
               !! PAN
               !!
               do  jk = 1, nlev
                  do  jp = 1, kproma
                     idt=pan_tag(ispec)%idx
                     if ( abs(pan_tagging_sum(jp,jk,jrow)) .gt. limit_o1D) then
                          pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                     &
                                          +  pan_tagging(jp,jk,ispec,jrow)                          &
                                          * (pan_0(jp,jk,jrow)/pan_tagging_sum(jp,jk,jrow)-1._dp)   & 
                                          * inv_time_step_len
                        !                   / time_step_len
                     else
          !  pan_tagging_sum(j,k,jrow) = sign(limit_o1D , pan_tagging_sum(j,k,jrow))
                     endif
                  enddo
               enddo

            enddo
   

      !  WRITE(*,*) 'vgvg local start scaling-result: pxtte(4,5,itrac_o3(1))',pxtte(4,5,itrac_o3(1))
       ! Hier noch Platz den Advection fehler zu berechnen 
         CASE(2)  
            ! do nothing
         END SELECT
       END SUBROUTINE tagging_local_start
  !***************************************************************************

  SUBROUTINE tagging_vdiff(i_tagging_ctrl)

    ! ECHAM5
    !changed: finish -> error
    !USE messy_main_mpi_bi,     ONLY: finish
    USE messy_main_blather_bi, ONLY:error_bi
    USE messy_main_grid_def_mem_bi,  ONLY: jrow, kproma, nlev
    USE messy_main_tracer_mem_bi,    ONLY: pxtte => qxtte

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)            :: i_tagging_ctrl             ! 1: pre mecca; 2: post Mecca

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'tagging_vdiff'
    INTEGER(I4)                    :: status, idt
    INTEGER                        :: ispec                      ! species in NOy or NMHC families
    INTRINSIC                      :: LBOUND, UBOUND
    status=0


    if (.not.ltagging) then
       print*,"tagging switched off so faar"
       return  
    end if

!   WRITE(*,*) 'TAGGING VDIFF ',i_tagging_ctrl
!   WRITE(*,*) 'Low/High Bounds pxte_pre.. ',LBOUND(pxtte_pre_emis_NOy),UBOUND(pxtte_pre_emis_NOy)
!   WRITE(*,*) 'Low/High Bounds pxtte ',LBOUND(pxtte(1:kproma,:,idx_NOy(1))), &
!                                UBOUND(pxtte(1:kproma,:,idx_NOy(1)))
!   WRITE(*,*) 'Emis_noy_onlem',LBOUND(emis_noy_onlem),UBound(emis_noy_onlem)


    SELECT CASE(i_tagging_ctrl)
    CASE(1)                    ! Pre-VDIFF Set pre-emission tendencies
        
       pxtte_pre_emis_NOy(1:kproma,:,jrow)=0._dp
       do ispec=1,Fam_NOy_species  
          idt=idx_NOy(ispec)
          pxtte_pre_emis_NOy(1:kproma,1:nlev,jrow)=pxtte_pre_emis_NOy(1:kproma,:,jrow) &
                                                  +pxtte(_RI_X_ZN_(1:kproma,:,idt))
       enddo

       pxtte_pre_emis_NMHC(1:kproma,:,jrow)=0._dp
       do ispec=1,Fam_NMHC_species
           idt=idx_NMHC(ispec)
           pxtte_pre_emis_NMHC(1:kproma,:,jrow)=pxtte_pre_emis_NMHC(1:kproma,:,jrow) &
                                               +pxtte(_RI_X_ZN_(1:kproma,:,idt))
        enddo


   
    CASE(2)                    ! Post onlem diagnose emissions and add to diag tracers
       emis_noy_onlem(1:kproma,:,jrow)=-pxtte_pre_emis_NOy(1:kproma,:,jrow)
        
       do ispec=1,Fam_NOy_species
          idt=idx_NOy(ispec)
          emis_noy_onlem(1:kproma,:,jrow)=emis_noy_onlem(1:kproma,:,jrow)  &
                                          + pxtte(_RI_X_ZN_(1:kproma,:,idt))                     
       enddo
           

       idt=noy_tag(isoi)%idx
       pxtte(_RI_X_ZN_(1:kproma,:,idt))=pxtte(_RI_X_ZN_(1:kproma,:,idt)) &
                                   + emis_noy_onlem(1:kproma,:,jrow)


       emis_nmhc_onlem(1:kproma,:,jrow)=-pxtte_pre_emis_NMHC(1:kproma,:,jrow)

         do ispec=1,Fam_NMHC_species
           idt=idx_NMHC(ispec)
           emis_nmhc_onlem(1:kproma,:,jrow)=emis_nmhc_onlem(1:kproma,:,jrow)  &
                                            + pxtte(_RI_X_ZN_(1:kproma,:,idt))                          
        enddo

         
        idt=nmhc_tag(isoi)%idx
        pxtte(_RI_X_ZN_(1:kproma,:,idt))=pxtte(_RI_X_ZN_(1:kproma,:,idt)) &
                                 +emis_nmhc_onlem(1:kproma,:,jrow)


    CASE DEFAULT
         status=1
    END SELECT
    IF (status/=0) CALL error_bi('incorrect value for i_tagging_ctrl', substr)

  END SUBROUTINE tagging_vdiff

! ---------------------------------------------------------------------------------

  SUBROUTINE tagging_physc(i_tagging_ctrl)

    ! ECHAM5
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_grid_def_mem_bi,  ONLY:  nlev, jrow, kproma, nproma    
     
    !changed -> time_step_len --> timer
    USE messy_main_timer,            ONLY: time_step_len
    USE messy_main_tracer_mem_bi,    ONLY: pxtte => qxtte, pxtm1 => qxtm1, &
                                           ntrac_gp, ti_gp

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)              :: i_tagging_ctrl    ! 1: pre mecca; 2: post Mecca

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'tagging_physc'

! op_mm_20171016
!    REAL(DP),DIMENSION(kproma,nlev) :: alpha,beta,gamma, D1, D2, P
    
    REAL(DP),DIMENSION(kproma,nlev) ::  D1, D2, P
   
    INTEGER(I4) :: status

    INTEGER                         :: ispec        ! tagging categories 
                                                    ! 1,..,N_tag_species

    INTEGER                         :: idt, idt2, jk, jp, j,k
    REAL(DP)                        :: toto3loss, totcoloss, &
                                       totnoyloss, totpanloss, &
                                       totnmhcloss, inv_time_step_len

    !op_mm_20171016
    REAL(DP)                        :: scalefactor 
    status=0


    if (.not.ltagging) then
       print*,"tagging switched off so faar"
       return  
    end if


  !WRITE(*,*) 'START TAGGING_PHYSC:',i_tagging_ctrl

!------------------------------------------------- C  A  S  E  ------------------------------------

      

    SELECT CASE(i_tagging_ctrl)
    CASE(1)                  ! Pre-MECCA: Update emissions       ! Set X1 values and account for loss terms
        

!   
!----------------------------------
! EMission update
!----------------------------------
 
! changed 
!rank identifiere telnox
             ! Lightning
       idt=noy_tag(ilig)%idx
       pxtte(_RI_X_ZN_(1:kproma,:,idt))=pxtte(_RI_X_ZN_(1:kproma,:,idt))  &
                                + telnox(_RI_XYZ__(1:kproma,jrow,:))

     
        
!----------------------------------
! Set X1 values
!----------------------------------
!--------
! Set species: ozone , carbon monoxide and PAN 
!--------
          idt=idx_O3
          o3_1(1:kproma,:,jrow) =  pxtm1(_RI_X_ZN_(1:kproma,:,idt)) & 
                                   + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
        
          idt=idx_CO
          co_1(1:kproma,:,jrow) =  pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &
                                    + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
          idt=idx_pan
          pan_1(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &
                                   + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len 
          idt=idx_ch4
          ch4_1(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt)) &
                                   + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len 
      

        !WRITE (*,*) "Physc pxtm1 03 ", minval(pxtm1(1:kproma,:,idx_O3)),MAXVAL(pxtm1(1:kproma,:,idx_O3))
        !WRITE (*,*) "Physc pxtte 03",  minval(pxtte(1:kproma,:,idx_O3)),MAXVAL(pxtte(1:kproma,:,idx_O3))
        !WRITE (*,*) "Physc o3_0 " ,  minval(o3_0(1:kproma,:,jrow)),MAXVAL(o3_0(1:kproma,:,jrow))
        !WRITE(*,*) 'vgvg physc 1 : o3_0,o3_1', o3_0(4,5,1),o3_1(4,5,1)
!--------
! Set families NOy and NMHC
!--------

!changed
!Rank identifier
!o3tagging(:,:,:,:)->o3tag(:)%ptr(:,:,:)

          noy_1(1:kproma,:,jrow) = 0.0_dp
          do ispec=1,Fam_NOy_species
             idt=idx_NOy(ispec)
             noy_1(1:kproma,:,jrow) = noy_1(1:kproma,:,jrow) + N_in_Fam(ispec) *          &
                                     (pxtm1(_RI_X_ZN_(1:kproma,:,idt))                             &
                                     + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len )
          enddo

      
          nmhc_1(1:kproma,:,jrow) = 0.0_dp
          do ispec=1,Fam_NMHC_species
             idt=idx_NMHC(ispec)
             nmhc_1(1:kproma,:,jrow) = nmhc_1(1:kproma,:,jrow) + C_in_Fam(ispec) *        &
                                      (pxtm1(_RI_X_ZN_(1:kproma,:,idt))                            &
                                        + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len )
          enddo
!--------
! Set diagnostic species: 
!--------
          o3_tagging_sum(1:kproma,:,jrow)=0.0_dp
          noy_tagging_sum(1:kproma,:,jrow)=0.0_dp
          nmhc_tagging_sum(1:kproma,:,jrow)=0.0_dp
          pan_tagging_sum(1:kproma,:,jrow)=0.0_dp
          co_tagging_sum(1:kproma,:,jrow)=0.0_dp


        do ispec=1,N_tag_species


           idt=o3_tag(ispec)%idx
           o3_tagging(1:kproma,:,ispec,jrow) =pxtm1(_RI_X_ZN_(1:kproma,:,idt))                    &
                                              + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
           o3_tagging_sum(1:kproma,:,jrow)   =o3_tagging_sum(1:kproma,:,jrow)              & 
                                              +o3_tagging(1:kproma,:,ispec,jrow)
          

           idt=noy_tag(ispec)%idx
           noy_tagging(1:kproma,:,ispec,jrow)  = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                   &
                                               + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
           noy_tagging_sum(1:kproma,:,jrow)    = noy_tagging_sum(1:kproma,:,jrow)          &
                                                +noy_tagging(1:kproma, :, ispec, jrow)


           idt=NMHC_tag(ispec)%idx
           nmhc_tagging(1:kproma,:,ispec,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                   &
                                                + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
           nmhc_tagging_sum(1:kproma,:,jrow)   =  nmhc_tagging_sum(1:kproma,:,jrow)        & 
                                                   +nmhc_tagging(1:kproma,:,ispec,jrow)
                   
           idt=pan_tag(ispec)%idx
           pan_tagging(1:kproma,:,ispec,jrow)  = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                   &
                                               + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
           pan_tagging_sum(1:kproma,:,jrow)    = pan_tagging_sum(1:kproma,:,jrow)          &
                                               +pan_tagging(1:kproma,:,ispec,jrow)


          idt=co_tag(ispec)%idx
          co_tagging(1:kproma,:,ispec,jrow)   = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                   &
                                              + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
          co_tagging_sum(1:kproma,:,jrow)     = co_tagging_sum(1:kproma,:,jrow)           &
                                               +co_tagging(1:kproma, :, ispec, jrow)


       enddo


        ! Account for tendencies, which are not from chemistry or transport
        !  xnew_tag = xold_tag * fac => tendency=(xold*fac-xold)/dt  
        !                                       =  xold*(fac-1)/dt;  fac=X1/X0 

! op_mm_20131125
! Similar as for the advection correction (local start) the division by a field (here o3_0) can cause divions by zero
! Therefore the tendencies are only calculated if x0 fields are larger then a certain threshold. 
! Also added multiplication of inverse time step instead of division by time step. 



       inv_time_step_len=1.0_dp/time_step_len

       do ispec=1,N_tag_species

           !!
           !! O3
           !!
          do jk=1, nlev
             do jp=1, kproma
                idt=o3_tag(ispec)%idx
!               if (o3_0(i,k,jrow) .eq. 0.0_dp) then
                if (abs(o3_0(jp,jk,jrow)) .ge. 1.0e-20) then
                   pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                            &
                                     +o3_tagging(jp,jk,ispec,jrow)                  &
                                     *(o3_1(jp,jk,jrow)/o3_0(jp,jk,jrow)-1._dp)     & 
                                     *inv_time_step_len
                     !               / time_step_len               
                else
                    ! do nothing 
                end if
             end do
          end do
       


           !!
           !!NOY
           !!
   
          do jk=1, nlev
             do jp=1, kproma
                idt=noy_tag(ispec)%idx
!               if (noy_0(i,k,jrow) .eq. 0.0_dp) then
                if (abs(noy_0(jp,jk,jrow)) .ge. 1.0e-20) then
                   pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                    & 
                                     +noy_tagging(jp,jk,ispec,jrow)                         &
                                     *(noy_1(jp,jk,jrow)/noy_0(jp,jk,jrow)-1._dp)           & 
                                     *inv_time_step_len
                 !                  / time_step_len
                 else
                    ! do nothing
              end if
           end do
        end do
       


        !!
        !!CO
        !!
         do jk=1, nlev
           do jp=1, kproma
              idt=co_tag(ispec)%idx
!              if (co_0(i,k,jrow) .eq. 0.0_dp) then
              if (abs(co_0(jp,jk,jrow)) .ge. 1.0e-20_dp) then
                 pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                 & 
                                  + co_tagging(jp,jk,ispec,jrow)                       &
                                  * (co_1(jp,jk,jrow)/co_0(jp,jk,jrow)-1._dp)          &
                                  *inv_time_step_len
                 !                  / time_step_len
              else
                  ! do nothing 
              end if
           end do
        end do

       
        !!
        !! NMHC
        !!
 
        do jk=1, nlev
           do jp=1, kproma     
              idt=nmhc_tag(ispec)%idx
!             if (nmhc_0(i,k,jrow) .eq. 0.0_dp) then
              if (abs(nmhc_0(jp,jk,jrow)) .ge. 1.0e-20_dp) then
                 pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                   & 
                                   +  nmhc_tagging(jp,jk,ispec,jrow)                     &
                                   * ( nmhc_1(jp,jk,jrow)/nmhc_0(jp,jk,jrow)-1._dp)      & 
                                   *inv_time_step_len
                     !                  / time_step_len
              else
                 ! do nothing 
              end if
           end do
        end do



        !!
        !!PAN
        !!
        do  jk= 1, nlev
           do  jp = 1, kproma  
              if (abs(pan_0(jp,jk,jrow)) .gt. limit_o1D) then 
                 idt=pan_tag(ispec)%idx
                 pxtte(_RI_X_ZN_(jp,jk,idt)) = pxtte(_RI_X_ZN_(jp,jk,idt))                                     &
                                   +  pan_tagging(jp,jk,ispec,jrow)                        &
                                   * ( pan_1(jp,jk,jrow)/pan_0(jp,jk,jrow)-1._dp)          &
                                   *inv_time_step_len
!                  / time_step_len
              else 
         ! pan_0(j,k,jrow)= sing(.... , pan_0(j,k,jrow))
              end if
           enddo
        enddo

     enddo
  
!  write(*,*) 'SHAPE_i_trac', SHAPE(itrac_o3(ispec))
!  write(*,*) 'SHAPE_pxtte', SHAPE(pxtte(1:kproma,:,itrac_o3(ispec)))
!  WRITE (*,*)"Physc o_3 diagnostic ",o3_tagging(1:kproma,:,ispec,jrow)
!  WRITE (*,*)"Physc pxtte diagnostc ", pxtte(1:kproma,:,itrac_o3(ispec))   
!  WRITE(*,*) "Physc_For ispec=3",MINVAL(o3_tagging(1:kproma,:,1:N_tag_species,jrow)),MAXVAL(o3_tagging(1:kproma,:,1:N_tag_species,jrow))
!  WRITE (*,*)"Physc o_3 diagnostic ",MINVAL(o3_tagging(1:kproma,:,:,jrow)),MAXVAL(o3_tagging(1:kproma,:,:,jrow))
!  WRITE (*,*)"Physc pxtte diagnostc ", MINVAL(pxtte(1:kproma,:,itrac_o3(ispec))),MAXVAL(pxtte(1:kproma,:,itrac_o3(ispec)))
!   WRITE (*,*)"Physc pxtte diagnostc OH", MINVAL(pxtte(1:kproma,:,itrac_ohe(ispec))),MAXVAL(pxtte(1:kproma,:,itrac_ohe(ispec)))
!   WRITE (*,*)"Physc pxtte diagnostc HO2 ", MINVAL(pxtte(1:kproma,:,itrac_ho2e(ispec))),MAXVAL(pxtte(1:kproma,:,itrac_ho2e(ispec)))!  WRITE(*,*) 'vgvg physc: pxtte(4,5,itrac_o3(1))',pxtte(4,5,itrac_o3(1))

! op_mm_20140114
!added return; otherwise pxtte of tagging tracer is set two times ech time step (first time with old destruction rates)
     return


!------------------------------------------------  C A S E   2  ------------------------------------------------
       
  CASE(2)                    ! Post-MECCA  calculate chemical changes
!--------
! Update 'Emissions' i.e. N2O -> NO
!       
!--------
! 1.Step set X2 values, note tagged has not to be set since X1=X2 :)      
!--------
! think about it
!
   

     idt=noy_tag(in2o)%idx
     idt2=idx_noprod_N2O
     pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt))+      &
                                    pxtte(_RI_X_ZN_(1:kproma,:,idt2))
        
     idt=idx_noprod_N2O
     noy_1(1:kproma,:,jrow)=noy_1(1:kproma,:,jrow) +         &
          pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len
!--------
! Set species: ozone , carbon monoxide and PAN 
!--------

     idt=idx_o3
     o3_2(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                      & 
                            +pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

     idt=idx_CO
     co_2(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                      &
                             + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

     idt=idx_pan
     pan_2(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                     &
                             + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

     idt=idx_ch4
     ch4_2(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                     &
                             + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len

     idt=idx_oh
     oh_2(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                      &
                            +pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len 

     idt=idx_ho2
     ho2_2(1:kproma,:,jrow) = pxtm1(_RI_X_ZN_(1:kproma,:,idt))                     &
                            + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len 

 
    
!--------
! Set families NOy and NMHC
!--------
        

     noy_2(1:kproma,:,jrow) = 0.0_dp
     do ispec=1,Fam_NOy_species
        idt=idx_NOy(ispec)
        noy_2(1:kproma,:,jrow) = noy_2(1:kproma,:,jrow) + N_in_Fam(ispec) *  &
                                (pxtm1(_RI_X_ZN_(1:kproma,:,idt))                     &
                                 + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len )
     enddo

     nmhc_2(1:kproma,:,jrow) = 0.0_dp
     do ispec=1,Fam_NMHC_species
        idt=idx_NMHC(ispec)
        nmhc_2(1:kproma,:,jrow) = nmhc_2(1:kproma,:,jrow) + C_in_Fam(ispec) * &
                                 (pxtm1(_RI_X_ZN_(1:kproma,:,idt))                     &
                                  + pxtte(_RI_X_ZN_(1:kproma,:,idt)) * time_step_len )
     enddo


!--------
! 2nd step:
! Adapt chemical tendencies to integration scheme 
!--------
! Set tendencies
!--------
!--------
! Ozone          
!--------
        
     idt=idx_o3prod_ho2
     o3prod_ho2(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
         
     idt=idx_o3prod_ro2
     idt2=idx_o3prod_MeO2
     o3prod_ro2(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))  &
                                     +pxtte(_RI_X_ZN_(1:kproma,:,idt2))
     idt=idx_o3prod_o2
     o3prod_o2(1:kproma,:,jrow)     = 2._dp*pxtte(_RI_X_ZN_(1:kproma,:,idt))   
                                           ! 2 O are produced.

        

!--------
! For unkonown reasons Prod might be negative.
!--------
!WRITE(*,*) 'Min/MAX',MINVAL(o3prod_ho2(1:kproma,:,jrow)),MAXVAL(o3prod_ho2(1:kproma,:,jrow))


     o3prod_ho2(1:kproma,:,jrow) = MAX(0._dp,o3prod_ho2(1:kproma,:,jrow))
     o3prod_ro2(1:kproma,:,jrow) = MAX(0._dp,o3prod_ro2(1:kproma,:,jrow))
     o3prod_o2(1:kproma,:,jrow)  = MAX(0._dp,o3prod_o2(1:kproma,:,jrow))
        

!      WRITE(*,*) 'Physc o_3prod Min/MAX',MINVAL(o3prod_ho2(1:kproma,:,jrow)), &
!                                         MAXVAL(o3prod_ho2(1:kproma,:,jrow))

        
     idt=idx_o3loss_ho2
     o3loss_ho2(1:kproma,:,jrow)   = pxtte(_RI_X_ZN_(1:kproma,:,idt)) 
     idt=idx_o3loss_ro
     o3loss_ro(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_o3loss_oh
     o3loss_oh(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_o3loss_no
     o3loss_no(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_o3loss_xo
     o3loss_xo(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))


     ! get reaction rates for HOx tagging
     idt=idx_LossG2100
     LossG2100(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2103
     LossG2103(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2104
     LossG2104(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2105
     LossG2105(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2106
     LossG2106(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2107
     LossG2107(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2109
     LossG2109(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2110
     LossG2110(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2111
     LossG2111(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG2112
     LossG2112(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG3200
     LossG3200(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG3201
     LossG3201(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG3202
     LossG3202(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG3203
     LossG3203(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG3207
     LossG3207(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG4101
     LossG4101(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG4110
     LossG4110(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG6203
     LossG6203(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG6204
     LossG6204(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossG7201
     LossG7201(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossJ3200
     LossJ3200(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossJ3201
     LossJ3201(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_LossJ4101b
     LossJ4101b(1:kproma,:,jrow)   = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_OHlossNMHC
     OHlossNMHC(1:kproma,:,jrow)        = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_HO2lossNMHC
     HO2lossNMHC(1:kproma,:,jrow)       = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_OHlossHO2prodNMHC
     OHlossHO2prodNMHC(1:kproma,:,jrow) = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_HO2prodNMHCNOy
     HO2prodNMHCNOy(1:kproma,:,jrow)    = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_HO2prodNMHCphoto
     HO2prodNMHCphoto(1:kproma,:,jrow)  = pxtte(_RI_X_ZN_(1:kproma,:,idt))



     o3loss(_RI_XYZ__(1:kproma,jrow,:))       = o3loss_ho2(1:kproma,:,jrow)       &
                                     +o3loss_ro(1:kproma,:,jrow)        &
                                     +o3loss_oh(1:kproma,:,jrow)        &
                                     +o3loss_no(1:kproma,:,jrow)        &
                                     +o3loss_xo(1:kproma,:,jrow)

     o3prod(_RI_XYZ__(1:kproma,jrow,:))       = o3prod_ho2(1:kproma,:,jrow)        &
                                      +o3prod_ro2(1:kproma,:,jrow)       &
                                      +o3prod_o2(1:kproma,:,jrow)

!        WRITE (*,*)"Physc O_3 loss/prod", o3loss(1:kproma,:,jrow), &
!                                          o3prod(1:kproma,:,jrow)
!        WRITE (*,*)"Physc O_3 loss", MINVAL(o3loss(1:kproma,:,jrow)),MAXVAL(o3loss(1:kproma,:,jrow))
!        WRITE (*,*)"Physc O_3 prod", MINVAL(o3prod(1:kproma,:,jrow)),MAXVAL(o3prod(1:kproma,:,jrow))
!--------
! PAN  
!--------
  
     idt=idx_panprod
     panprod(1:kproma,:,jrow)   = pxtte(_RI_X_ZN_(1:kproma,:,idt))
     idt=idx_panloss
     panloss(1:kproma,:,jrow)   = pxtte(_RI_X_ZN_(1:kproma,:,idt))

     panprod(1:kproma,:,jrow)   = MAX(0._dp,panprod(1:kproma,:,jrow))
!--------
! HOx 
!-------- 
!op_mm_13102015
! use full tracer values instead of tendency only! 
     idt=idx_OH
    !! OH(1:kproma,:,jrow)=pxtte(_RI_X_ZN_(1:kproma,:,idt))
    OH(1:kproma,:,jrow)= pxtm1(_RI_X_ZN_(1:kproma,:,idt))+pxtte(_RI_X_ZN_(1:kproma,:,idt))*time_step_len
     idt=idx_HO2
    !! HO2(1:kproma,:,jrow)=pxtte(_RI_X_ZN_(1:kproma,:,idt))
    HO2(1:kproma,:,jrow)= pxtm1(_RI_X_ZN_(1:kproma,:,idt))+pxtte(_RI_X_ZN_(1:kproma,:,idt))*time_step_len


!------
! O1D
!-----
     idt=idx_O1D
     O1D(1:kproma,:,jrow)=pxtte(_RI_X_ZN_(1:kproma,:,idt))

     
     
     
!op_mm_20171016+
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!removed all the tagging_tendencies call because the calculated alpha and beta values are not used after they were calculated
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   


     
!--------
! 3rd step Adapt tendencies (calc alpha, beta) obtained from MECCA diagnostics
!--------
!     call tagging_tendencies(kproma,nlev,time_step_len,         &
!             o3_1(1:kproma,:,jrow),o3_2(1:kproma,:,jrow),       &
!             o3prod(_RI_XYZ__(1:kproma,jrow,:)),o3loss(_RI_XYZ__(1:kproma,jrow,:)), & 
!             alpha(1:kproma,:),beta(1:kproma,:))

        
!       Write(*,*) 'i,j,jrow=1,1,2',o3prod(1,1,2),alpha(1,1),o3loss(1,1,2),beta(1,1),o3_2(1,1,2),o3_1(1,1,2)
        !o3loss_ho2(1:kproma,:,jrow)   = beta(1:kproma,:) * o3loss_ho2(1:kproma,:,jrow)
        !o3loss_ro(1:kproma,:,jrow)    = beta(1:kproma,:) * o3loss_ro(1:kproma,:,jrow)
        !o3loss_oh(1:kproma,:,jrow)    = beta(1:kproma,:) * o3loss_oh(1:kproma,:,jrow)
        !o3loss_no(1:kproma,:,jrow)    = beta(1:kproma,:) * o3loss_no(1:kproma,:,jrow)
        !o3loss_xo(1:kproma,:,jrow)    = beta(1:kproma,:) * o3loss_xo(1:kproma,:,jrow)
        !o3loss(1:kproma,:,jrow)       = beta(1:kproma,:) * o3loss(1:kproma,:,jrow)

!        WRITE(*,*) 'vgvg physc-2: a,b',alpha(4,5),beta(4,5)
!        WRITE(*,*) 'vgvg physc-2: o3prod(4,5,1)',o3prod(4,5,1)
!        WRITE(*,*) 'vgvg physc-2: o3loss(4,5,1)',o3loss(4,5,1)


!todo 1:kproma
!     call tagging_tendencies(kproma,nlev,time_step_len,                                             &
!                             pan_1(:,:,jrow),pan_2(:,:,jrow),panprod(:,:,jrow),panloss(:,:,jrow),   &
!                             alpha(1:kproma,:),beta(1:kproma,:))



 !call tagging_tendencies(kproma,nlev,time_step_len,                                             &
 !                            pan_1(1:kproma,:,jrow),pan_2(1:kproma,:,jrow),panprod(1:kproma,:,jrow),panloss(1:kproma,:,jrow),   &
 !                            alpha(1:kproma,:),beta(1:kproma,:))

      
        !panprod(1:kproma,:,jrow) = alpha(1:kproma,:) * panprod(1:kproma,:,jrow)
       ! panloss(1:kproma,:,jrow) = beta(1:kproma,:) * panloss(1:kproma,:,jrow)

        !write(*,*) 'alpha_es',MINVAL(alpha(1:kproma,:)),MAXVAL(alpha(1:kproma,:)) 
        !write(*,*) 'bita_e5',MINVAL(beta(1:kproma,:)),MAXVAL(beta(1:kproma,:))
        !write(*,*) 'panprod_es',MINVAL(panprod(1:kproma,:,jrow)),MAXVAL(panprod(1:kproma,:,jrow)) 
        !write(*,*) 'panloss_es_',MINVAL(panloss(1:kproma,:,jrow)),MAXVAL(panloss(1:kproma,:,jrow)) 
!op_mm_20171016-
        
        
        
!--------
! 4th step Calculated tendencies based on chemical changes
!--------


     call tagging_get_tendencies(kproma,nlev,time_step_len,       &
            ch4_1(1:kproma,:,jrow),ch4_2(1:kproma,:,jrow),        &
            co_1(1:kproma,:,jrow),co_2(1:kproma,:,jrow),          &
            pan_1(1:kproma,:,jrow),pan_2(1:kproma,:,jrow),        &
            nmhc_1(1:kproma,:,jrow),nmhc_2(1:kproma,:,jrow),      &
            ch4loss(_RI_XYZ__(1:kproma,jrow,:)),                            &
            coprod(_RI_XYZ__(1:kproma,jrow,:)),                             &
            coloss(_RI_XYZ__(1:kproma,jrow,:)))   
          
        
!-------
! Calculate Correction factors for HOx chemistry
!-------


     D1(1:kproma,:)= ch4loss(_RI_XYZ__(1:kproma,jrow,:))      &
                    +coprod(_RI_XYZ__(1:kproma,jrow,:))       &
                   +coloss(_RI_XYZ__(1:kproma,jrow,:))
              
     D2(1:kproma,:)= o3prod_ho2(1:kproma,:,jrow)    &
                    +o3prod_ho2(1:kproma,:,jrow)
           
     P (1:kproma,:)= LossG2111(1:kproma,:,jrow)       &                          
                    + o3prod_ho2(1:kproma,:,jrow)   &
                    + o3loss_ho2(1:kproma,:,jrow)
 
! op_vr_20170405 auskommentiert
! because alpha, beta, gamma are not used anymore 
 !    call tagging_HOx_tendencies(kproma,nlev,           &
 !         oh_2(1:kproma,:,jrow),ho2_2(1:kproma,:,jrow), &
 !         LossG2111(1:kproma,:,jrow),                     &
 !         D1(1:kproma,:),D2(1:kproma,:),                &                             
 !         alpha(1:kproma,:),beta(1:kproma,:),gamma(1:kproma,:))
        
      

!--------
! 5th step Calculate tagging chemistry and integrate species = calc tendency
!--------   


!    write(*,*)'O3_tagging_before_chem',  MINVAL(o3_tagging(1:kproma,:,:,jrow)),MAXVAL(o3_tagging(1:kproma,:,:,jrow))
!    write(*,*)'O3_tagging_sum_before_chem',  MINVAL(o3_tagging_sum(1:kproma,:,jrow)),MAXVAL(o3_tagging_sum(1:kproma,:,jrow))
! write(*,*)'pan_tagging_before_chem',  MINVAL(pan_tagging(1:kproma,:,:,jrow)),MAXVAL(pan_tagging(1:kproma,:,:,jrow))
 !   write(*,*)'pan_tagging_sum_before_chem',  MINVAL(pan_tagging_sum(1:kproma,:,jrow)),MAXVAL(pan_tagging_sum(1:kproma,:,jrow))

!op_mm_20131209
!copy data of typre structures to pointer for parsing them to the tagging chemistry routine

         
     oh_tag_data(1:kproma,:,:)=0.0_DP
     ho2_tag_data(1:kproma,:,:)=0.0_DP
     p_ho2_ten_data(1:kproma,:,:)= 0.0_DP
     p_ro2_ten_data(1:kproma,:,:)=0.0_DP
     d_ro_ten_data(1:kproma,:,:)=0.0_DP
     d_ho2_ten_data(1:kproma,:,:)=0.0_DP
     d_xo_ten_data(1:kproma,:,:)=0.0_DP
     d_oh_ten_data(1:kproma,:,:)=0.0_DP
     d_no_ten_data(1:kproma,:,:)=0.0_DP
     ch4loss_tag_data(1:kproma,:,:)=0.0_DP !op_vr_20170215
     frac_ODD_data(1:kproma,:,:)=0.0_DP     !op_vr_20170215



     do ispec=1, n_tag_species
        oh_tag_data(1:kproma,:, ispec)=oh_tag(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))
        ho2_tag_data(1:kproma,:, ispec)=ho2_tag(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))
        p_ho2_ten_data(1:kproma,:, ispec)=P_ho2_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
        p_ro2_ten_data(1:kproma,:, ispec)=P_ro2_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
        d_ro_ten_data(1:kproma,:, ispec)=D_ro_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
        d_ho2_ten_data(1:kproma,:, ispec)=D_ho2_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
        d_xo_ten_data(1:kproma,:, ispec)=D_xo_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
        d_oh_ten_data(1:kproma,:, ispec)=D_oh_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
        d_no_ten_data(1:kproma,:, ispec)=D_no_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
        ch4loss_tag_data(1:kproma,:, ispec)=ch4loss_tag(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
        frac_ODD_data(1:kproma,:, ispec)=frac_ODD(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:)) 
     end do


     call tagging_chemistry(kproma,nlev,time_step_len,                                                          &
        o3_tagging(1:kproma,:,:,jrow),noy_tagging(1:kproma,:,:,jrow), pan_tagging(1:kproma,:,:,jrow),           &
        nmhc_tagging(1:kproma,:,:,jrow),co_tagging(1:kproma,:,:,jrow), co_tagging_sum(1:kproma,:,jrow),         &
        nmhc_tagging_sum(1:kproma,:,jrow),noy_tagging_sum(1:kproma,:,jrow),o3_tagging_sum(1:kproma,:,jrow),     &
        pan_tagging_sum(1:kproma,:,jrow),oh(1:kproma,:,jrow), ho2(1:kproma,:,jrow),do3(1:kproma,:,:,jrow),      &
        dnoy(1:kproma,:,:,jrow), dpan(1:kproma,:,:,jrow),dnmhc(1:kproma,:,:,jrow), dco(1:kproma,:,:,jrow),      &
        oh_tag_data(1:kproma, :,:), ho2_tag_data(1:kproma,:,:), o3prod_ho2(1:kproma,:,jrow),                    &
        o3prod_ro2(1:kproma,:,jrow), o3prod_o2(1:kproma,:,jrow),                                                &
        o3loss_oh(1:kproma,:,jrow),o3loss_ho2(1:kproma,:,jrow),o3loss_ro(1:kproma,:,jrow),                      &
        o3loss_no(1:kproma,:,jrow),o3loss_xo(1:kproma,:,jrow),o3loss(_RI_XYZ__(1:kproma,jrow,:)),                         &
        panprod(1:kproma,:,jrow),panloss(1:kproma,:,jrow),coprod(_RI_XYZ__(1:kproma,jrow,:)),                             &
        coloss(_RI_XYZ__(1:kproma,jrow,:)),ch4loss(_RI_XYZ__(1:kproma,jrow,:)),ch4loss_tag_data(1:kproma,:,:),                      &
        LossG2100(1:kproma,:,jrow),LossG2103(1:kproma,:,jrow),LossG2104(1:kproma,:,jrow),                       &
        LossG2105(1:kproma,:,jrow),LossG2106(1:kproma,:,jrow),LossG2107(1:kproma,:,jrow),                       &
        LossG2109(1:kproma,:,jrow),LossG2110(1:kproma,:,jrow),LossG2111(1:kproma,:,jrow),                       &
        LossG2112(1:kproma,:,jrow),LossG3200(1:kproma,:,jrow),LossG3201(1:kproma,:,jrow),                       &
        LossG3202(1:kproma,:,jrow),LossG3203(1:kproma,:,jrow),LossG3207(1:kproma,:,jrow),                       &
        LossG4101(1:kproma,:,jrow),LossG4110(1:kproma,:,jrow),LossG6203(1:kproma,:,jrow),                       &
        LossG6204(1:kproma,:,jrow),LossG7201(1:kproma,:,jrow),LossJ3200(1:kproma,:,jrow),                       &
        LossJ3201(1:kproma,:,jrow),LossJ4101b(1:kproma,:,jrow),                                                 &
        OHlossNMHC(1:kproma,:,jrow),HO2lossNMHC(1:kproma,:,jrow),                                               &
        OHlossHO2prodNMHC(1:kproma,:,jrow),HO2prodNMHCNOy(1:kproma,:,jrow),                                     &
        HO2prodNMHCphoto(1:kproma,:,jrow),frac_ODD_data(1:kproma,:,:),                                          &
        p_ho2_ten_data(1:kproma, :,:), p_ro2_ten_data(1:kproma,:,:), d_ro_ten_data(1:kproma,:,:),               &
        d_oh_ten_data(1:kproma,:,:), d_ho2_ten_data(1:kproma, :,:), d_xo_ten_data(1:kproma, :,:), d_no_ten_data(1:kproma, :,:), &
        status)


!op_mm_20131209
!copy data back for channel output
     do ispec=1, n_tag_species
        oh_tag(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))=oh_tag_data(1:kproma,:, ispec)
        ho2_tag(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))=ho2_tag_data(1:kproma,:, ispec)
        P_ho2_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))= p_ho2_ten_data(1:kproma,:, ispec)
        P_ro2_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))= p_ro2_ten_data(1:kproma,:, ispec) 
        D_ro_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))=d_ro_ten_data(1:kproma,:, ispec)
        D_ho2_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))=d_ho2_ten_data(1:kproma,:, ispec)
        D_xo_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))=d_xo_ten_data(1:kproma,:, ispec)
        D_oh_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))=d_oh_ten_data(1:kproma,:, ispec) 
        D_no_ten(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))= d_no_ten_data(1:kproma,:, ispec)
        ch4loss_tag(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))=ch4loss_tag_data(1:kproma,:, ispec)
        frac_ODD(ispec)%ptr(_RI_XYZ__(1:kproma,jrow,:))=frac_ODD_data(1:kproma,:, ispec)
     end do

     IF (status/=0) CALL error_bi('Error in Tagging Chemistry ', substr)

 
  CASE DEFAULT
     status=1
  END SELECT
 
!-----
! Update Tendencies
!-----
!! 1/10 loop
!! op_mm_20171016
!change 1/10 with 1/(ntag_catagories) (because number of categories is not fixed anymore) 

scalefactor=1.0_dp/dble(N_tag_species)


  do ispec=1,N_tag_species
     do  k=1, nlev
        do  j=1,kproma
             if ( abs(do3(j,k,ispec,jrow)*time_step_len) .gt. abs(scalefactor*O3_tagging_sum(j,k,jrow)))     then   
    
                do3(j,k,ispec,jrow)=     & 
                     sign(scalefactor*O3_tagging_sum(j,k,jrow)/time_step_len, do3(j,k,ispec,jrow)*O3_tagging_sum(j,k,jrow) )
             end if

             if ( abs(dpan(j,k,ispec,jrow)*time_step_len) .gt. abs(scalefactor*PAN_tagging_sum(j,k,jrow)))   then

                dpan(j,k,ispec,jrow)=    &
                     sign(scalefactor*PAN_tagging_sum(j,k,jrow)/time_step_len, dpan(j,k,ispec,jrow)*PAN_tagging_sum(j,k,jrow) )
             end if

             if ( abs(dNOy(j,k,ispec,jrow)*time_step_len) .gt. abs(scalefactor*NOy_tagging_sum(j,k,jrow)) )     then         
 
                dNOy(j,k,ispec,jrow)=sign(scalefactor*NOy_tagging_sum(j,k,jrow)/time_step_len, &
                     dNOy(j,k,ispec,jrow)*NOy_tagging_sum(j,k,jrow) )
             end if

             if (     abs(dNMHC(j,k,ispec,jrow)*time_step_len) .gt. abs(scalefactor*NMHC_tagging_sum(j,k,jrow)) & 
                  .or. abs(dNMHC(j,k,ispec,jrow)*time_step_len) .lt. tag_limit)     then        
  
                dNMHC(j,k,ispec,jrow)= &
                     sign(scalefactor*NMHC_tagging_sum(j,k,jrow)/time_step_len, dNMHC(j,k,ispec,jrow)*NMHC_tagging_sum(j,k,jrow) )
             end if

             if (      abs(dco(j,k,ispec,jrow)*time_step_len) .gt. abs(scalefactor*co_tagging_sum(j,k,jrow)) & 
                  .or. abs(dco(j,k,ispec,jrow)*time_step_len) .lt. tag_limit )     then       
   
                dco(j,k,ispec,jrow)=sign(scalefactor*CO_tagging_sum(j,k,jrow)/time_step_len, &
                     dco(j,k,ispec,jrow)*CO_tagging_sum(j,k,jrow) )
             end if
             
             if  (O1D(j,1,jrow)  .le.  limit_o1D) then    !day/night indicator
                !do3(j,k,ispec,jrow)= 0.0_dp
                !dpan(j,k,ispec,jrow) = 0.0_dp
                !dnoy(j,k,ispec,jrow)= 0.0_dp
                dnmhc(j,k,ispec,jrow) = 0.0_dp
                dco(j,k,ispec,jrow)=0.0_dp
             end if

             if  (O1D(j,k,jrow)  .le.  limit_o1D) then 
                dnoy(j,k,ispec,jrow)= 0.0_dp
                dpan(j,k,ispec,jrow)= 0.0_dp
             end if
          enddo
       enddo

!!
! Updating Tendencies
!!
       idt=noy_tag(ispec)%idx
       pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt))  + dnoy(1:kproma,:,ispec,jrow)

       idt=o3_tag(ispec)%idx
       pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt))   + do3(1:kproma,:,ispec,jrow)

       idt=pan_tag(ispec)%idx
       pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt))  + dpan(1:kproma,:,ispec,jrow)

       idt=nmhc_tag(ispec)%idx
       pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt)) + dnmhc(1:kproma,:,ispec,jrow)

       idt=co_tag(ispec)%idx
       pxtte(_RI_X_ZN_(1:kproma,:,idt)) = pxtte(_RI_X_ZN_(1:kproma,:,idt))   +  dco(1:kproma,:,ispec,jrow)

   enddo
 
 

! op_mm_20131206
!removed unnecessary output (do3_out etc) 
    !do3_out(1:kproma,:,jrow)=do3(1:kproma,:,iout,jrow)
    !dnoy_out(1:kproma,:,jrow)=dnoy(1:kproma,:,iout,jrow)
    !dpan_out(1:kproma,:,jrow)=-dnoy(1:kproma,:,iout,jrow)   !dpan=-dnoy
    !dnmhc_out(1:kproma,:,jrow)=dnmhc(1:kproma,:,iout,jrow)
    !dco_out(1:kproma,:,jrow)=dco(1:kproma,:,iout,jrow)
  



!----------------------------------------------------------------------------------------------------
! Calculate Error
!------

! op_mm_20131126
! added additinal if check for all tendency updates to overcome divisions by zero 


   IF (l_adv_err_diag) then

      do k=1,nlev
         do j=1,kproma
            jk=k
            toto3loss=MAX(0._dp,(O3_0(j,k,jrow)-O3_1(j,k,jrow))/time_step_len) &
                     +o3loss(_RI_XYZ__(j,jrow,jk))
!!$            if (j.eq.44.and.k.eq.1) write(*,*) 'isnotnan0',j,k,jrow,toto3loss 
!!$        if (isnan(toto3loss)) write(*,*) 'isnan1',j,k,jrow

            toto3loss=toto3loss-MIN(0._dp,do3_err(j,k,jrow))

!!$           if (j.eq.44.and.k.eq.1) write(*,*) 'isnotnan1',toto3loss
!!$        if (isnan(toto3loss)) write(*,*) 'isnan2',j,k,jrow

            if (o3_2(j,k,jrow).gt.0._dp.and.o3_2(j,k,jrow).lt.2.e-5) then
               if ((o3_2(j,k,jrow) .ne. 0.0_dp)) toto3loss=toto3loss/o3_2(j,k,jrow)
            end if

!!$            if (j.eq.44.and.k.eq.1) write(*,*) 'isnotnan2',toto3loss,o3_2(j,k,jrow)
!!$        if (isnan(toto3loss)) write(*,*) 'isnan3',j,k,jrow
!!$            if (j.eq.44.and.k.eq.1) write(*,*) 'isnotnan3',j,k,jrow,toto3loss
        
            idt=itrac_o3_errP
            toto3loss=toto3loss*(pxtm1(_RI_X_ZN_(j,jk,idt))+pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)

!!$        if (isnan(toto3loss)) write(*,*) 'isnan4',j,k,jrow,pxtte(j,k,itrac_o3_errP),do3_err(j,k,jrow), &
!!$                                       pxtm1(j,k,itrac_o3_errP),pxtte(j,k,itrac_o3_errP), &
!!$                                       pxtm1(j,k,itrac_o3_errP)+pxtte(j,k,itrac_o3_errP)*time_step_len

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! D03 ERR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            if (do3_err(j,k,jrow).gt.0.) then
               idt=itrac_o3_errP
               !if o3_1 ) = 0 -> pxtee = pxtte
               if (O3_1(j,k,jrow) .ne. 0.0_dp) then
                  pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                  &
                       + (do3_err(j,k,jrow)-toto3loss)               &
                       /(1._dp+toto3loss*time_step_len/O3_1(j,k,jrow))
               else
                  pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))  
               end if

            else
               idt=itrac_o3_errN
               if (O3_1(j,k,jrow) .ne. 0.0_dp) then
                  pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))           &
                       - (do3_err(j,k,jrow)-toto3loss)        &
                       /(1._dp+toto3loss*time_step_len/O3_1(j,k,jrow))
               else
                  pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))      
               end if
!!$               if (isnan( pxtte(_RI_X_ZN_(j,jk,idt)))) write (*,*) "istnan 2 "
            endif
!!$        if (isnan(pxtte(j,k,itrac_o3_errP))) write(*,*) 'isnan5',j,k,jrow,do3_err(j,k,jrow),toto3loss


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! NOY ERR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
            totnoyloss=(noy_0(j,k,jrow)-noy_1(j,k,jrow))/time_step_len              &
                       +panprod(j,k,jrow)                                           &
                        -telnox(_RI_XYZ__(j,jrow,jk))                                          &
                        -emis_noy_onlem(j,k,jrow)

            totnoyloss=MAX(0._dp,totnoyloss)
            if (dnoy_err(j,k,jrow).gt.0.) then
               idt=itrac_noy_errP
            
             ! if noy_2 = 0 -> pxxte 0 pxtte + dnoy_err
               if (noy_2(j,k,jrow) .ne. 0.0_dp) then
                  pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                    &
                                  + dnoy_err(j,k,jrow)                 &
                                  - (totnoyloss/noy_2(j,k,jrow))       &
                                  *(pxtm1(_RI_X_ZN_(j,jk,idt))                   &
                                  +pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)
               else
                  !!What should happen if noy_2 is to small? pxtte=pxtte or pxtte=pxtte+dnoy_er??
                  pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                    &
                                  +dnoy_err(j,k,jrow)   
               end if
            else
               idt=itrac_noy_errN
               !if noy_2 = 0 -> pxtte = pxtte - dnoy_err
               if (noy_2(j,k,jrow) .ne. 0.0_dp) then 
                  pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                                               &
                                  - dnoy_err(j,k,jrow)                                            &
                                  - ((totnoyloss-dnoy_err(j,k,jrow))/noy_2(j,k,jrow))             &
                                  *(pxtm1(_RI_X_ZN_(j,jk,idt))                                              &
                                  +pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)

               else 
                  pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                                               &
                                  -dnoy_err(j,k,jrow)   
             end if
             
          endif
        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CO ERR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++           

          totcoloss=(co_0(j,k,jrow)-co_1(j,k,jrow))/time_step_len                 &
                    + coloss(_RI_XYZ__(j,jrow,jk))
          totcoloss=MAX(0._dp,totcoloss)

          if (dco_err(j,k,jrow).gt.0.) then
             idt=itrac_co_errP
             ! if co2 = 0 -> pxtte = pxtte + dco_err
             if (co_2(j,k,jrow) .ne. 0.0_dp) then 
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                  &
                                 + dco_err(j,k,jrow)               &
                                 - (totcoloss/co_2(j,k,jrow))      &
                                 *(pxtm1(_RI_X_ZN_(j,jk,idt))                &
                                 +pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)
             else 
                !what to to if co_2 is too small? 
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                  &
                                + dco_err(j,k,jrow)    
             end if
          else
             idt=itrac_co_errN
             ! if co2 = 0 -> pxtte = pxtte - dco_err
             if (co_2(j,k,jrow) .ne. 0.0_dp) then 
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                                      &
                                - dco_err(j,k,jrow)                                    &
                                - ((totcoloss-dco_err(j,k,jrow))/co_2(j,k,jrow))       &
                                *(pxtm1(_RI_X_ZN_(j,jk,idt))                                     &
                                +pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)
             else
                !what to to if co_2 is too small? 
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                  &
                     - dco_err(j,k,jrow)  
             end if
          endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!NMHC ERR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++           


        ! nmhcloss=coprod
          totnmhcloss=(NMHC_0(j,k,jrow)-nmhc_1(j,k,jrow))/time_step_len              &
                      -emis_nmhc_onlem(j,k,jrow)                                     &
                      +coprod(_RI_XYZ__(j,jrow,jk))
          totnmhcloss=MAX(0._dp,totnmhcloss)
          if (dnmhc_err(j,k,jrow).gt.0.) then
             idt=itrac_nmhc_errP
             ! if nmhc_2 = 0 -> pxtte = pxtte + dnmhc_err
             if (nmhc_2(j,k,jrow) .ne. 0.0_dp) then 
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                                    &
                                + dnmhc_err(j,k,jrow)                                &
                                - (totnmhcloss/nmhc_2(j,k,jrow))                     &
                                *(pxtm1(_RI_X_ZN_(j,jk,idt))                                   &
                                +pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)
             else
                !what happens if nmhc_2 is to small?
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))              &
                                + dnmhc_err(j,k,jrow)    
             end if
          else
             idt=itrac_nmhc_errN
             idt2=itrac_o3_errN
             ! if nmhc_2 = 0 -> pxtte = pxtee - dnmhc_err
             if (nmhc_2(j,k,jrow) .ne. 0.0_dp) then
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt2))                                              &
                                - dnmhc_err(j,k,jrow)                                           &
                                - ((totnmhcloss-dnmhc_err(j,k,jrow))/nmhc_2(j,k,jrow))          &
                                *(pxtm1(_RI_X_ZN_(j,jk,idt))                                              &
                                +pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)
             else
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt2))                &
                     -dnmhc_err(j,k,jrow)     
             end if
          endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!PAN ERR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++           


          totpanloss=(pan_0(j,k,jrow)-pan_0(j,k,jrow))/time_step_len              &
               +panloss(j,k,jrow)
          totpanloss=MAX(0._dp,totpanloss)

          if (dpan_err(j,k,jrow).gt.0.) then
             idt=itrac_pan_errP
             ! if pan_2 = 0 -> pxtte = pxtte + dpan_err 
             if (pan_2(j,k,jrow) .ne. 0.0_dp) then  
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                                    &
                                + dpan_err(j,k,jrow)                                 &
                                - (totpanloss/pan_2(j,k,jrow))                       &
                                *(pxtm1(_RI_X_ZN_(j,jk,idt))                                   &
                                +pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)
             else 
                !! What should happen if  pan_2 is too small?
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                &
                     + dpan_err(j,k,jrow)    
             end if

          else
             idt=itrac_pan_errN
             ! if pan_2 = 0 -> pxtte = pxtte + do3err
             if (pan_2(j,k,jrow) .ne. 0.0_dp) then  
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                                               &
                                - do3_err(j,k,jrow)                                             &
                                - ((totpanloss-dpan_err(j,k,jrow))/pan_2(j,k,jrow))             &
                                *(pxtm1(_RI_X_ZN_(j,jk,idt))                                              &
                                 +pxtte(_RI_X_ZN_(j,jk,idt))*time_step_len)
             else
                !! What happen if pan_2 is too small? 
                pxtte(_RI_X_ZN_(j,jk,idt))=pxtte(_RI_X_ZN_(j,jk,idt))                                               &
                                - do3_err(j,k,jrow)                   !think about why d03 err
             end if
          endif
       enddo
    enddo
 ENDIF
 IF (status/=0) CALL error_bi('incorrect value for i_tagging_ctrl', substr)



END SUBROUTINE tagging_physc
! ----------------------------------------------------------------------
 
  SUBROUTINE tagging_free_memory

  ! tagging MODULE ROUTINE
  ! free memory of global fields


  ! Author: Volker Grewe, DLR-OP, NCAR Boulder visit, 2008
  
  IMPLICIT NONE
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_free_memory'
  integer                     :: n

    DEALLOCATE(o3_tagging)
    DEALLOCATE(noy_tagging)
    DEALLOCATE(nmhc_tagging)
    DEALLOCATE(co_tagging)
    DEALLOCATE(pan_tagging)

    !!$DEALLOCATE(oh_tagging)
    !!$DEALLOCATE(ho2_tagging)

    DEALLOCATE(o3_tagging_sum)
    DEALLOCATE(noy_tagging_sum)
    DEALLOCATE(nmhc_tagging_sum)
    DEALLOCATE(co_tagging_sum)
    DEALLOCATE(pan_tagging_sum)

   !!$ DEALLOCATE(oh_tagging_sum)
   !!$ DEALLOCATE(ho2_tagging_sum)

    DEALLOCATE(do3)
    DEALLOCATE(dnoy)
    DEALLOCATE(dco)
    DEALLOCATE(dpan)
    DEALLOCATE(dnmhc)

   !!$DEALLOCATE(do3_out)
   !!$DEALLOCATE(dnoy_out)
   !!$DEALLOCATE(dpan_out)
   !!$DEALLOCATE(dnmhc_out)
   !!$DEALLOCATE(dco_out)

    DEALLOCATE(o3_0)
    DEALLOCATE(noy_0)
    DEALLOCATE(pan_0)
    DEALLOCATE(nmhc_0)
    DEALLOCATE(co_0)

    DEALLOCATE(o3_1)
    DEALLOCATE(noy_1)
    DEALLOCATE(pan_1)
    DEALLOCATE(nmhc_1)
    DEALLOCATE(co_1)
    DEALLOCATE(ch4_1)

    DEALLOCATE(o3_2)
    DEALLOCATE(noy_2)
    DEALLOCATE(pan_2)
    DEALLOCATE(nmhc_2)
    DEALLOCATE(co_2)
    DEALLOCATE(ch4_2)
    DEALLOCATE(ho2_2)
    DEALLOCATE(oh_2)

    DEALLOCATE(do3_err)
    DEALLOCATE(dnoy_err)
    DEALLOCATE(dco_err)
    DEALLOCATE(dpan_err)
    DEALLOCATE(dnmhc_err)

    DEALLOCATE(emis_noy_onlem)
    DEALLOCATE(emis_nmhc_onlem)

    DEALLOCATE(pxtte_pre_emis_NOy)
    DEALLOCATE(pxtte_pre_emis_NMHC)

    DEALLOCATE(o3prod_ho2)
    DEALLOCATE(o3prod_ro2)
    DEALLOCATE(o3prod_o2)
  
    DEALLOCATE(o3loss_ho2)
    DEALLOCATE(o3loss_ro)
    DEALLOCATE(o3loss_oh)
    DEALLOCATE(o3loss_no)
    DEALLOCATE(o3loss_xo)

    DEALLOCATE(panprod)
    DEALLOCATE(panloss)
    !!$DEALLOCATE(coprod)
    !!$DEALLOCATE(coloss)
    DEALLOCATE(nmhcloss)
    !!$DEALLOCATE(ch4loss)

    DEALLOCATE(LossG2100)
    DEALLOCATE(LossG2103)
    DEALLOCATE(LossG2104)
    DEALLOCATE(LossG2105)
    DEALLOCATE(LossG2106)
    DEALLOCATE(LossG2107)
    DEALLOCATE(LossG2109)
    DEALLOCATE(LossG2110)
    DEALLOCATE(LossG2111)
    DEALLOCATE(LossG2112)
    DEALLOCATE(LossG3200)
    DEALLOCATE(LossG3201)
    DEALLOCATE(LossG3202)
    DEALLOCATE(LossG3203)
    DEALLOCATE(LossG3207)
    DEALLOCATE(LossG4101)
    DEALLOCATE(LossG4110)
    DEALLOCATE(LossG6203)
    DEALLOCATE(LossG6204)
    DEALLOCATE(LossG7201)
    DEALLOCATE(LossJ3200)
    DEALLOCATE(LossJ3201)
    DEALLOCATE(LossJ4101b)
    DEALLOCATE(OHlossNMHC)
    DEALLOCATE(HO2lossNMHC)
    DEALLOCATE(OHlossHO2prodNMHC)
    DEALLOCATE(HO2prodNMHCNOy)
    DEALLOCATE(HO2prodNMHCphoto)

    DEALLOCATE(HO2)
    DEALLOCATE(OH)   
    DEALLOCATE(ch4)   
    DEALLOCATE(co)  
    DEALLOCATE(NMHC) 
    DEALLOCATE(NOy)
    DEALLOCATE(O3)
    DEALLOCATE(PAN) 
    DEALLOCATE(O1D)

    
! op_mm_20171016+
    DEALLOCATE(o3_tag) 
    DEALLOCATE(noy_tag)
    DEALLOCATE(pan_tag)
    DEALLOCATE(nmhc_tag)
    DEALLOCATE(co_tag)
    DEALLOCATE(oh_tag)
    DEALLOCATE(ho2_tag)
    DEALLOCATE(ch4loss_tag)
    DEALLOCATE(p_ho2_ten)
    DEALLOCATE(p_ro2_ten)
    DEALLOCATE(d_oh_ten)
    DEALLOCATE(d_ho2_ten)
    DEALLOCATE(d_no_ten)
    DEALLOCATE(d_ro_ten)
    DEALLOCATE(d_xo_ten)
    DEALLOCATE(itrac_NOy)
    DEALLOCATE(itrac_NMHC)
    DEALLOCATE(itrac_PAN)
    DEALLOCATE(itrac_CO)
    DEALLOCATE(itrac_O3)

    do n =1, N_tag_species
       IF (ASSOCIATED(CATEGORY_SET(n)%category)) then
          DEALLOCATE(CATEGORY_SET(n)%category)
       END IF
    END DO

    DEALLOCATE(CATEGORY_SET)

! op_mm_20171016-
  


  END SUBROUTINE tagging_free_memory
! ************************************************************************
! PRIVATE ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE tagging_read_nml_e5(status, iou)

    ! DRADON MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5/ATTILA
    !
    ! Author: Patrick Joeckel, MPICH, Oct 2003

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close
  
    IMPLICIT NONE

    ! I/O
    INTEGER(I4), INTENT(OUT) :: status     ! error status
    INTEGER(I4), INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_read_nml_e5'

    NAMELIST /CPL/  c_lnox, i_diag !c_tropop

    LOGICAL              :: lex      ! file exists ?
    INTEGER(I4)          :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
   ! WRITE(*,*) 'C_tropop: ',C_tropop(1),C_tropop(2)  
    !      WRITE(*,*) 'C_B: ',C_B(1),C_B(2)
    WRITE(*,*) 'C_lnox',C_lnox(1),C_lnox(2)

    SELECT CASE (i_diag)
    CASE (0) 
       write (*,*) 'No additional diagnostic output'
    CASE (1)
       write (*,*) 'Additional diagnostic output for tagging'
    CASE DEFAULT 
       write (*,*) 'i_diag out of range ',i_diag
       RETURN
    END SELECT

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE tagging_read_nml_e5
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_tagging_si
