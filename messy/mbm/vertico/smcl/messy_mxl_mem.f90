MODULE messy_mxl_mem

  !-------------------------------------------------------------------------------------------------------
  !  mxl : MiXed Layer model for the diurnal dynamics of the convective boundary layer  
  !
  !  AUTHOR:  Ruud Janssen, MPIC, Sept 2013
  !           Andrea Pozzer, the other guy did nothing
  !-------------------------------------------------------------------------------------------------------

  ! MESSy
  USE messy_main_constants_mem,  ONLY: DP, SP,  &
                                 STRLEN_MEDIUM,STRLEN_LONG,STRLEN_VLONG,STRLEN_ULONG

  IMPLICIT NONE 

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'mxl'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.1'

  INTEGER, PARAMETER :: NMAX_IMPORT= 100
  INTEGER, PARAMETER :: REPR_SC     = 1
  INTEGER, PARAMETER :: REPR_2D     = 2
  INTEGER, PARAMETER :: REPR_3DMID  = 3
  INTEGER, PARAMETER :: REPR_3DINT  = 4

  PUBLIC :: dp, sp, STRLEN_MEDIUM

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! CTRL NAMELIST 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), PUBLIC :: lat = 0.0_DP
  REAL(DP), PUBLIC :: lon = 0.0_DP

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL_IC NAMELIST => MXL INITIAL CONDITIONS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp) :: hbl_ic       = 0.0_dp   ! initial boundary layer height (m)
  REAL(dp) :: psurf        = 0.0_dp   ! (initial) surface pressure (Pa)
  REAL(dp) :: wthetasmax   = 0.0_dp   ! maximum surface sensible heat flux (K m s-1)
  REAL(dp) :: wqsmax       = 0.0_dp   ! maximum surface specific humidity flux (g kg-1 m s-1)  
  REAL(dp) :: beta         = 0.0_dp   ! ratio of entrainment to surface heat flux (-)  
  REAL(dp) :: omega        = 0.0_dp   ! subsidence rate (s-1)  
  REAL(dp) :: thetam_ic    = 0.0_dp   ! initial mixed layer potential temperature (K)
  REAL(dp) :: dtheta_ic    = 0.0_dp   ! initial potential temperature jump (K)
  REAL(dp) :: gammatheta   = 0.0_dp   ! lapse rate potential temperature in free troposphere (K m-1)
  REAL(dp) :: gammath2     = 0.0_dp   ! lapse rate potential temperature in free troposphere (K m-1)
  LOGICAL  :: l_gamma                  ! switch for second lapse rate
  REAL(dp) :: hcrit        = 0.0_dp   ! crictical boundary layer height (m)
  REAL(dp) :: advtheta     = 0.0_dp   ! advection of heat (K s-1)  
  REAL(dp) :: qm_ic        = 0.0_dp   ! initial mixed layer specific moisture (g kg-1)
  REAL(dp) :: dq_ic        = 0.0_dp   ! initial specific moisture jump (g kg-1)
  REAL(dp) :: gammaq       = 0.0_dp   ! lapse rate specific moisture in free troposphere (g kg-1 m-1)
  REAL(dp) :: advq         = 0.0_dp   ! advection of moisture (g kg-1 s-1)
  CHARACTER(LEN=STRLEN_MEDIUM)  :: f_wthetas    = 'NOFLUX'    ! function heat flux
  CHARACTER(LEN=STRLEN_MEDIUM)  :: f_wqs        = 'NOFLUX'    ! function moisture flux
  REAL(dp) :: starttime_wths = 0.0_dp   ! start time heat flux after start simulation (s) 
  REAL(dp) :: stoptime_wths  = 0.0_dp   ! stop time heat flux after start simulation (s) 
  REAL(dp) :: starttime_wqs  = 0.0_dp   ! start time moisture flux after start simulation (s) 
  REAL(dp) :: stoptime_wqs   = 0.0_dp   ! stop time moisture flux after start simulation (s) 
  REAL(dp) :: starttime_adv  = 0.0_dp   ! start time advection after start simulation (s) 
  REAL(dp) :: stoptime_adv   = 0.0_dp   ! stop time advection after start simulation (s) 

  LOGICAL  :: l_ustconst
  REAL(dp) :: um_ic       = 0.0_dp   ! initial mixed-layer wind velocity in x-direction (m s-1)   
  REAL(dp) :: vm_ic       = 0.0_dp   ! initial mixed-layer wind velocity in y-direction (m s-1)
  REAL(dp) :: ug          = 0.0_dp   ! geostrophic wind velocity in x-direction (m s-1)
  REAL(dp) :: vg          = 0.0_dp   ! geostrophic wind velocity in y-direction (m s-1)
  REAL(dp) :: uws_ic      = 0.0_dp   ! initial surface momentum flux in x-direction (m2 s-2)
  REAL(dp) :: vws_ic      = 0.0_dp   ! initial surface momentum flux in y-direction (m2 s-2)
  REAL(dp) :: gammau      = 0.0_dp   ! wind velocity lapse rate x-direction (m s-1 m-1)
  REAL(dp) :: gammav      = 0.0_dp   ! wind velocity lapse rate y-direction (m s-1 m-1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL_RAD NAMELIST 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL  :: l_radiation
  REAL(dp) :: Cc          = 0.0_dp   ! cloud cover (-)
  REAL(dp) :: salbedo     = 0.0_dp   ! surface albedo (-)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL_SURFLAYER NAMELIST 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL  :: l_surfacelayer
  REAL(dp) :: z0m       = 0.0_dp   ! roughness length for momentum (m)
  REAL(dp) :: z0h       = 0.0_dp   ! roughness length for heat and moisture (m)
  REAL(dp) :: Ch_ic     = 0.0_dp   ! drag coefficient for heat and moisture (-)  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL_LANDSURFACE NAMELIST (see Van Heerwaarden et al (2010), J. Hydrometeorology)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL  :: l_landsurface
  REAL(dp) :: Tsurf     = 0.0_dp   ! surface temperature (K)
  REAL(dp) :: wwilt     = 0.0_dp   ! volumetric soil moisture at wilting point (m3 m-3)     
  REAL(dp) :: w1        = 0.0_dp   ! volumetric soil moisture layer 1 (top) (m3 m-3)   
  REAL(dp) :: w2        = 0.0_dp   ! volumetric soil moisture layer 2 (m3 m-3)   
  REAL(dp) :: wfc       = 0.0_dp   ! volumetric soil moisture at field capacity (m3 m-3)   
  REAL(dp) :: wsat      = 0.0_dp   ! saturated volumetric water content (m3 m-3)    
  REAL(dp) :: CLa       = 0.0_dp   ! Clapp-Hornberger retention curve parameter (-)    
  REAL(dp) :: CLb       = 0.0_dp   ! Clapp-Hornberger retention curve parameter (-)        
  REAL(dp) :: CLc       = 0.0_dp   ! Clapp-Hornberger retention curve parameter (-)        
  REAL(dp) :: C1sat     = 0.0_dp   ! coefficient force term moisture (-)                                      
  REAL(dp) :: C2ref     = 0.0_dp   ! coefficient restore term moisture (-)   
  REAL(dp) :: gD        = 0.0_dp   ! correction factor for vapor pressure deficit (-)    
  REAL(dp) :: rsmin     = 0.0_dp   ! minimum resistance transpiration (s m-1)    
  REAL(dp) :: rssoilmin = 0.0_dp   ! minimum resistance soil evaporation (s m-1)    
  REAL(dp) :: LAI       = 0.0_dp   ! leaf area index of the vegetated fraction (-)     
  REAL(dp) :: cveg      = 0.0_dp   ! vegetation fraction (-)    
  REAL(dp) :: Tsoil1    = 0.0_dp   ! soil temperature layer 1 (top) (K)    
  REAL(dp) :: Tsoil2    = 0.0_dp   ! soil temperature layer 2 (K)    
  REAL(dp) :: Wl        = 0.0_dp   ! equivalent water layer depth for wet vegetation (m)    
  REAL(dp) :: Lambda    = 0.0_dp   ! thermal diffusivity of the skin layer (W m-2 K-1)    
  REAL(dp) :: CGsat     = 0.0_dp   ! saturated soil conductivity for heat (K m-2 J-1)    
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL_DDEP NAMELIST 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp) :: hc       = 0.0_dp   ! canopy height (m)
  REAL(dp) :: drag     = 0.0_dp   ! drag (?)
  REAL(dp) :: soilph   = 0.0_dp   ! soil pH (-)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL_MEGAN NAMELIST 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp) :: laip  = 0.0_dp   ! LAI of previous month (-)
  REAL(dp) :: btr_frac = 0.0_dp   ! broadleaf fraction (-)
  REAL(dp) :: ntr_frac = 0.0_dp   ! needleleaf fraction (-)
  REAL(dp) :: shr_frac = 0.0_dp   ! shrub fraction (-)
  REAL(dp) :: hrb_frac = 0.0_dp   ! herb/grass/crop fraction (-)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL_ONEMIS NAMELIST 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp) :: CH4_conc       = 0.0_dp   ! climatological CH4 concentratoon (mol mol-1 (?))
  REAL(dp) :: NOemisclass1   = 0.0_dp   ! ratio veg./emis class 1 
  REAL(dp) :: NOemisclass2   = 0.0_dp   ! ratio veg./emis class 2
  REAL(dp) :: emis_NO_cult   = 0.0_dp   ! cultivation intensity
  REAL(dp) :: emis_NO_fert   = 0.0_dp   ! fertilizer application

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL_ORACLE NAMELIST 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp) :: OA_bg          = 0.0_dp   ! background organic aerosol (ug m-3)
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MXL DYNAMICS MEMORY
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp), DIMENSION(:,:),   POINTER :: Cc_mem         => NULL()   
  REAL(dp), DIMENSION(:,:),   POINTER :: hbl_mem        => NULL()   
  REAL(dp), DIMENSION(:,:),   POINTER :: hsl_mem        => NULL()   
  REAL(dp), DIMENSION(:,:),   POINTER :: wthetasmax_mem => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wqsmax_mem     => NULL()   
  REAL(dp),                   POINTER :: beta_mem       => NULL()  
  REAL(dp),                   POINTER :: omega_mem      => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: thetam_mem     => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: thetate_mem    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: dtheta_mem     => NULL()
  REAL(dp),                   POINTER :: gammatheta_mem => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: hcrit_mem      => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: advtheta_mem   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: qm_mem         => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: qte_mem        => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: dq_mem         => NULL()
  REAL(dp),                   POINTER :: gammaq_mem     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: advq_mem       => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: rh_mem         => NULL()

  REAL(dp), DIMENSION(:,:),   POINTER :: we_mem         => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: ws_mem         => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: thetavm_mem    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: dthetav_mem    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wthetas_mem    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wqs_mem        => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wthetavs_mem   => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wthetave_mem   => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: wqe_mem        => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: betaq_mem      => NULL() 
  REAL(dp), DIMENSION(:,:),   POINTER :: SH_mem         => NULL() 
  REAL(dp), DIMENSION(:,:),   POINTER :: LE_mem         => NULL() 

  REAL(dp), DIMENSION(:,:,:), POINTER :: um_mem         => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: du_mem         => NULL()
  REAL(dp),                   POINTER :: gammau_mem     => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: vm_mem         => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: dv_mem         => NULL()
  REAL(dp),                   POINTER :: gammav_mem     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: uws_mem        => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: vws_mem        => NULL() 
  REAL(dp), DIMENSION(:,:),   POINTER :: ueff_mem       => NULL()
  REAL(DP), POINTER, DIMENSION(:,:,:) :: pressi_3d_mem  => NULL() ! pressure at interface 
  REAL(DP), POINTER, DIMENSION(:,:,:) :: press_3d_mem   => NULL() ! pressure 
  REAL(DP), POINTER, DIMENSION(:,:,:) :: geopoti_3d_mem => NULL() ! geopotential at interface 
  REAL(DP), POINTER, DIMENSION(:,:,:) :: geopot_3d_mem  => NULL() ! geopotential 
  REAL(DP), POINTER, DIMENSION(:,:)   :: az0_mem        => NULL() ! roughness length 
  REAL(DP), POINTER, DIMENSION(:,:)   :: z0h_mem        => NULL() ! roughness length 
  REAL(DP), POINTER, DIMENSION(:,:)   :: z0m_mem        => NULL() ! roughness length 
  REAL(DP), POINTER, DIMENSION(:,:)   :: Ch_mem         => NULL() ! drag coefficient for heat
  ! mxl_radiation
  REAL(DP), POINTER, DIMENSION(:,:)   :: Swin_mem       => NULL() ! incoming shortwave radiation (W m-2)
  REAL(DP), POINTER, DIMENSION(:,:)   :: Swout_mem      => NULL() ! outgoing shortwave radiation (W m-2)
  REAL(DP), POINTER, DIMENSION(:,:)   :: Lwin_mem       => NULL() ! incoming longwave radiation (W m-2)
  REAL(DP), POINTER, DIMENSION(:,:)   :: Lwout_mem      => NULL() ! outgoing longwave radiation (W m-2)
  REAL(DP), POINTER, DIMENSION(:,:)   :: Rn_mem         => NULL() ! net radiation at the surface (W m-2)
  REAL(DP), POINTER, DIMENSION(:,:)   :: albedo_mem     => NULL() ! surface albedo
  REAL(DP), POINTER                   :: cossza_mem     => NULL() ! cosine of solar zenith angle (-)
  ! mxl_surfacelayer
  REAL(DP), POINTER, DIMENSION(:,:)   :: T2m_mem        => NULL() 
  REAL(DP), POINTER, DIMENSION(:,:)   :: q2m_mem        => NULL() 
  REAL(DP), POINTER, DIMENSION(:,:)   :: u2m_mem        => NULL() 
  REAL(DP), POINTER, DIMENSION(:,:)   :: v2m_mem        => NULL() 
  REAL(DP), POINTER, DIMENSION(:,:)   :: u10m_mem       => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: v10m_mem       => NULL()  
  REAL(DP), POINTER, DIMENSION(:,:)   :: rh2m_mem       => NULL() 
  REAL(DP), POINTER, DIMENSION(:,:)   :: Rib_mem        => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: thetavsurf_mem => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: Cm_mem         => NULL()  
  REAL(DP), POINTER, DIMENSION(:,:)   :: ustar_mem      => NULL()  
  REAL(DP), POINTER, DIMENSION(:,:)   :: thetasurf_mem  => NULL()  
  REAL(DP), POINTER, DIMENSION(:,:)   :: qsurf_mem      => NULL()  
  REAL(DP), POINTER, DIMENSION(:,:)   :: e2m_mem        => NULL()  
  REAL(DP), POINTER, DIMENSION(:,:)   :: esat2m_mem     => NULL()  
  ! mxl_landsurface
  REAL(DP), POINTER, DIMENSION(:,:)   :: ra_mem         => NULL() 
  REAL(DP), POINTER, DIMENSION(:,:)   :: rs_mem         => NULL() 
  REAL(DP), POINTER, DIMENSION(:,:)   :: rssoil_mem     => NULL() 
  REAL(DP), POINTER, DIMENSION(:,:)   :: Tsurf_mem      => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: wwilt_mem      => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: w2_mem         => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: w1_mem         => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: wfc_mem        => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: wsat_mem       => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: CLa_mem        => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: CLb_mem        => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: CLc_mem        => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: C1sat_mem      => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: C2ref_mem      => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: gD_mem         => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: rsmin_mem      => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: rssoilmin_mem  => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: LAI_mem        => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: cveg_mem       => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: Tsoil1_mem     => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: Tsoil2_mem     => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: Wl_mem         => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: Lambda_mem     => NULL()  
  REAL(DP), POINTER, DIMENSION(:,:)   :: CGsat_mem      => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: GR_mem         => NULL()

  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  !! DDEP MEMORY: messy_main_data
  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), POINTER, DIMENSION(:,:)   :: hc_mem         => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: drag_mem       => NULL()
  REAL(DP), POINTER, DIMENSION(:,:,:) :: soilpH_mem     => NULL()

  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  !! MEGAN MEMORY
  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), POINTER, DIMENSION(:,:)   :: laip_mem         => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: btr_frac_mem     => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: ntr_frac_mem     => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: shr_frac_mem     => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: hrb_frac_mem     => NULL()

  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  !! ONEMIS MEMORY
  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), POINTER, DIMENSION(:,:)   :: ch4_conc_mem       => NULL()
  REAL(DP), POINTER, DIMENSION(:,:,:) :: NOemisclass1_mem   => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: NOemisclass2_mem   => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: emis_NO_cult_mem   => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: emis_NO_fert_mem   => NULL()    

  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  !! OFFEMIS MEMORY
  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), POINTER, DIMENSION(:,:)   :: NOemis             => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)   :: O3emis             => NULL()

  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  !! ORACLE MEMORY
  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp), POINTER, DIMENSION(:,:,:) :: OA_bg_mem          => NULL()    

  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  !! JVAL MEMORY
  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), POINTER :: cdisse     => NULL()
  
  ! switch for output
  LOGICAL, PUBLIC :: l_verbose = .true.   ! GLOBAL SWITCH
  LOGICAL, PUBLIC :: l_chem_ft = .false.   ! GLOBAL SWITCH

  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  !! SIMPLE EMISSIONS
  !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  LOGICAL  :: l_emis_simple    
  REAL(dp) :: starttime_emis = 0.0_dp   ! start time heat flux after start simulation (s) 
  REAL(dp) :: stoptime_emis  = 0.0_dp   ! stop time heat flux after start simulation (s) 

END MODULE
