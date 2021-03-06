
!VERSION NUMBER:
!  $Id: donner_types.h,v 13.0 2006/03/28 21:08:59 fms Exp $


!#######################################################################

type donner_initialized_type

!***********************************************************************
!    these variables are defined during initialization, usually based on
!    input namelist values. they may change during model integration.
!    they are preserved across timesteps.
!***********************************************************************

!   do_donner_tracer            tracers are to be transported by the
!                               donner convection scheme ?
!   do_input_cell_liquid_size   cell liquid size to be input ?
!   do_bower_cell_liquid_size   cell liquid size from bower calculation ?
!   do_input_cell_ice_size      cell ice size to be input ?
!   do_default_cell_ice_size    use default cell ice size ?
!   coldstart                   is this a donner_deep coldstart ?
!   total_pts                   total number of points in the processor's
!                               subdomain
!   pts_processed_conv          number of points processed during current
!                               convection calculation
!   conv_alarm                  time remaining until next convection 
!                               calculation [ sec ]


logical  :: do_donner_tracer
logical  :: do_input_cell_liquid_size
logical  :: do_bower_cell_liquid_size
logical  :: do_input_cell_ice_size
logical  :: do_default_cell_ice_size
logical  :: coldstart

integer  :: total_pts
integer  :: pts_processed_conv
integer  :: conv_alarm
integer  :: physics_dt

end type donner_initialized_type


!#######################################################################

type donner_save_type

!***********************************************************************
!    these arrays contain information used to provide the variables
!    needed by moist_processes_mod on every timestep, so they must be
!    preserved across timesteps, in case donner_deep is not called on 
!    every step.
!***********************************************************************
       
!   lag_temp       temperature field to be used in calculating lag
!                  values of cape in order to determine time tendency
!                  of cape [ deg K ]
!   lag_vapor      water vapor mixing ratio field to be used in 
!                  calculating lag values of cape and column vapor
!                  in order to determine their time tendencies
!                  [ kg(h2o) / kg (dry air) ]
!   lag_press      full level pressure field to be used in 
!                  calculating lag values of cape and column vapor
!                  in order to determine their time tendencies
!                  [ Pa ]
!   cememf         normalized moisture forcing, cells+meso 
!                  NOTE: final value of cememf has terms related to 
!                  flux convergence of condensate and mesoscale 
!                  detrainment removed when parameterization is run 
!                  in a model using strat_cloud_mod, since in that 
!                  case those terms will be calculated within that
!                  module.
!                  [ kg(H2O)/kg/sec ]
!   cemetf         normalized thermal forcing, cells+meso [ K/sec ]
!                  (excludes convergence of surface heat flux)
!                  Cumulus thermal forcing defined as in Fig. 3 of 
!                  Donner (1993, JAS).
!                  NOTE: final value of cemetf has terms related to 
!                  flux convergence of condensate and mesoscale 
!                  detrainment removed when parameterization is run 
!                  in a model using strat_cloud_mod, since in that 
!                  case those terms will be calculated within that
!                  module.
!   det_mass_flux  detrained mass flux [ kg(air)/(m**2 sec) ]
!   tprea1         area weighted total normalized precipitation 
!                  [ mm/day ]
!   mass_flux      total mass flux at model levels, convective +
!                  mesoscale  [ kg/(m**2)/sec) ]
!   cell_up_mass_flux
!                  upward cell mass flux at interface levels
!                  [ kg / (m**2 sec) ]
!   dql_strat      increment in cloud liquid over the physics timestep
!                  resulting from donner convection [ kg/kg/sec ]
!   dqi_strat      increment in cloud ice over the physics timestep
!                  resulting from donner convection [ kg/kg/sec ]
!   dqa_strat      increment in cloud area over the physics timestep
!                  resulting  from donner convection [ 1/sec ]
!   humidity_area  fraction of cloud box area taken up by the
!                  cell + meso fraction. it is needed to determine
!                  the large-scale specific humidity field for the
!                  grid box. do not use for radiation calculation,
!                  since this area includes mesoscale downdraft, which
!                  does not contain cloud water.
!                  [ fraction ]
!   humidity_ratio ratio of large-scale specific humidity to specific 
!                  humidity in environment outside convective system
!                  [ dimensionless ]
!   tracer_tends   tendencies to tracer fields due to donner_deep
!                  mod [ kg/kg/sec ]
!   parcel_disp    time-integrated low level displacement [ Pa ]
!   tracername     names of the tracers
!   tracer_units   units associated with each tracer
 

real, dimension(:,:,:),          pointer  ::  cemetf=>NULL()
real, dimension(:,:,:),          pointer  ::  lag_temp=>NULL()
real, dimension(:,:,:),          pointer  ::  lag_vapor=>NULL()
real, dimension(:,:,:),          pointer  ::  lag_press=>NULL()
real, dimension(:,:,:),          pointer  ::  cememf=>NULL()
real, dimension(:,:,:),          pointer  ::  mass_flux=>NULL()
real, dimension(:,:,:),          pointer  ::  cell_up_mass_flux=>NULL()
real, dimension(:,:,:),          pointer  ::  det_mass_flux=>NULL()
real, dimension(:,:,:),          pointer  ::  dql_strat=>NULL()
real, dimension(:,:,:),          pointer  ::  dqi_strat=>NULL()
real, dimension(:,:,:),          pointer  ::  dqa_strat=>NULL()
real, dimension(:,:,:),          pointer  ::  humidity_area=>NULL()
real, dimension(:,:,:),          pointer  ::  humidity_ratio=>NULL()
real, dimension(:,:,:,:),        pointer  ::  tracer_tends=>NULL()
real, dimension(:,:),            pointer  ::  parcel_disp=>NULL()
real, dimension(:,:),            pointer  ::  tprea1=>NULL()
character(len=32), dimension(:), pointer  :: tracername=>NULL()
character(len=32), dimension(:), pointer  :: tracer_units=>NULL()

end type donner_save_type


!#######################################################################

type donner_rad_type

!***********************************************************************
!    these variables are needed to connect the donner deep convection   
!    cloud fields to the radiation package. they are allocated and 
!    deallocated on each deep convection timestep.
!***********************************************************************

!   cell_cloud_frac     fractional area of convective cells in grid
!                       box [ dimensionless ]
!   cell_liquid_amt     liquid water content of convective cells
!                       [ kg(h2o)/kg(air) ]
!   cell_liquid_size    assumed effective size of cell liquid drops
!                       [ microns ]
!   cell_ice_amt        ice water content of cells
!                       [ kg(h2o)/kg(air) ]
!   cell_ice_size       generalized effective diameter for ice in
!                       convective cells [ microns ]
!   meso_cloud_frac     fractional area of mesoscale clouds in grid
!                       box [ dimensionless ]
!   meso_liquid_amt     liquid water content in mesoscale clouds
!                       currently set to 0.0
!                       [ kg(h2o)/kg(air) ]
!   meso_liquid_size    assumed effective size of mesoscale drops
!                       currently set to 0.0 [ microns ]
!   meso_ice_amt        ice water content of mesoscale elements
!                       [ kg(h2o)/kg(air) ]
!   meso_ice_size       generalized ice effective size for anvil ice
!                       [ microns ]
!   nsum                number of time levels of data contained in
!                       the accumulation arrays; needed when time-
!                       averaging of cloud properties to be used in
!                       radiation package is desired


real,    dimension(:,:,:), pointer      :: cell_cloud_frac=>NULL()
real,    dimension(:,:,:), pointer      :: cell_liquid_amt=>NULL()
real,    dimension(:,:,:), pointer      :: cell_liquid_size=>NULL()
real,    dimension(:,:,:), pointer      :: cell_ice_amt=>NULL()
real,    dimension(:,:,:), pointer      :: cell_ice_size=>NULL()
real,    dimension(:,:,:), pointer      :: meso_cloud_frac=>NULL()
real,    dimension(:,:,:), pointer      :: meso_liquid_amt=>NULL()
real,    dimension(:,:,:), pointer      :: meso_liquid_size=>NULL()
real,    dimension(:,:,:), pointer      :: meso_ice_amt=>NULL()
real,    dimension(:,:,:), pointer      :: meso_ice_size=>NULL()
integer, dimension(:,:)  , pointer      :: nsum=>NULL()

end type donner_rad_type


!#######################################################################


type donner_nml_type

!************************************************************************
!    documentation for these variables may be found in the namelist 
!    section of donner_deep.f90. they exist for the duration of the run.
!************************************************************************


integer             ::  model_levels_in_sfcbl
integer             ::  parcel_launch_level  
logical             ::  allow_mesoscale_circulation
integer             ::  donner_deep_freq
character(len=32)   ::  entrainment_constant_source
character(len=16)   ::  cell_liquid_size_type
character(len=16)   ::  cell_ice_size_type 
real                ::  cell_liquid_eff_diam_input
real                ::  cell_ice_geneff_diam_input
real                ::  meso_liquid_eff_diam_input
logical             ::  do_average

end type donner_nml_type


!#######################################################################

type donner_param_type

!***********************************************************************
!    documentation for these variables may be found in the "private data"
!    section of donner_deep.f90.
!***********************************************************************


integer          :: istart
real             :: cp_vapor
real             :: cp_air  
real             :: rdgas     
real             :: rvgas     
real             :: parcel_dp
real             :: upper_limit_for_lfc
real             :: pstop
real             :: grav
real             :: kappa
real             :: dens_h2o
real             :: pie
real             :: seconds_per_day
real             :: hlv
real             :: hls
real             :: hlf 
real             :: kelvin
real             :: cld_base_vert_vel
real             :: dp_of_cloud_model
real             :: cloud_base_radius
real             :: wdet
real             :: rbound
real             :: wbound 
real             :: freeze_fraction 
real             :: virt_mass_co 
real             :: pdeep_mc 
real             :: tr_insert_time 
real             :: autoconv_rate 
real             :: autoconv_threshold 
real             :: tfre 
real             :: dfre
real             :: evap_in_downdrafts 
real             :: evap_in_environ    
real             :: entrained_into_meso 
real             :: d622
real             :: d608
real             :: upper_limit_for_lcl
real             :: tmin
real             :: anvil_precip_efficiency
real             :: meso_lifetime
real             :: meso_ref_omega
real             :: tprime_meso_updrft
real             :: meso_sep 
real             :: ref_press
real             :: meso_down_evap_fraction
real             :: meso_up_evap_fraction
integer          :: kpar
real             :: pdeep_cv
real             :: cdeep_cv
real             :: max_entrainment_constant_gate
real             :: max_entrainment_constant_kep
real             :: r_conv_land
real             :: r_conv_ocean
real             :: N_land
real             :: N_ocean
real             :: delz_land
real             :: delz_ocean
real             :: cell_liquid_eff_diam_def
real             :: cell_ice_geneff_diam_def
integer          :: anvil_levels
  
real, dimension(:), pointer :: arat=>NULL()
real, dimension(:), pointer :: ensemble_entrain_factors_gate=>NULL()
real, dimension(:), pointer :: ensemble_entrain_factors_kep=>NULL()
real, dimension(:), pointer :: dgeice=>NULL()
real, dimension(:), pointer :: relht=>NULL()

end type donner_param_type


!#######################################################################

type donner_column_diag_type

!***********************************************************************
!    these variables are used in defining the column diagnostics that
!    may be requested.
!***********************************************************************
 
!   in_diagnostics_window     are column diagnostics desired anywhere 
!                             in the current window ?
!   num_diag_pts              total number of activated diagnostics 
!                             columns
!   ncols_in_window           number of activated diagnostic columns in 
!                             this window
!   kstart                    output array elements for model levels 
!                             with k index > kstart will be written
!                             to the output file
!   i_dc                      column's i index (window coordinates)
!   j_dc                      column's j index (window coordinates)
!   unit_dc                   column's output file unit number
!   igl_dc                    column's i index (processor coordinates)
!   jgl_dc                    column's j index (processor coordinates)


logical                         :: in_diagnostics_window
integer                         :: num_diag_pts
integer                         :: ncols_in_window
integer                         :: kstart
integer, dimension(:), pointer  :: i_dc=>NULL()
integer, dimension(:), pointer  :: j_dc=>NULL()
integer, dimension(:), pointer  :: unit_dc=>NULL()
integer, dimension(:), pointer  :: igl_dc=>NULL()
integer, dimension(:), pointer  :: jgl_dc=>NULL()
   
end type donner_column_diag_type


!#######################################################################

type donner_conv_type

!***********************************************************************
!    these variables are allocated and deallocated on each donner step.
!    they contain variables needed during the calculation of deep con-
!    vection and its effects.
!**********************************************************************
 
!   cecon          normalized cell condensation/deposition
!                  [ K/sec ]
!   ceefc          normalized cell entropy-flux convergence [ K/sec ]
!                  (excludes convergence of surface flux) Entropy-flux
!                  convergence divided by (p0/p)**(rd/cp).
!   cell_ice_geneff_diam  
!                  cell ice generalized effective diameter.
!                  default is smallest data value given in andge
!                  subroutine (micrometers)
!   cell_liquid_eff_diam  
!                  cell liquid effective diameter. defined 
!                  using Bower parameterization or from namelist
!                  input. (micrometers)
!   cememf_mod     normalized moisture forcing, cells+meso, after 
!                  any modifications needed to prevent the creation
!                  of negative values [ kg(H2O)/kg/sec ]
!   cemfc          normalized cell moisture-flux convergence
!                  (excludes convergence of surface moisture flux)
!                  [ kg(H2O)/kg/sec ]
!   cmus           normalized mesoscale-updraft deposition
!                  [ kg(H2O)/kg/sec]
!   conv_moist_forcing
!                  total cell + meso mixing ratio tendency due to donner
!                  parameterization, including those terms which are
!                  recalculated in the strat_cloud parameterization.
!                  [ kg(h2o) / (kg(dry air) sec) ]
!   conv_temp_forcing
!                  total cell + meso temperature tendency due to donner
!                  parameterization, including those terms which are
!                  recalculated in the strat_cloud parameterization.
!                  [ K / sec ]
!   cual           cloud fraction, cells+meso, normalized by a(1,p_b)
!   cuqi           cell ice content (kg(ice)/kg)
!   cuql           cell liquid content (kg(water)/kg)
!   detmfl         detrained cell mass flux [ kg(air) / (m**2 sec) ]
!   dgeice         mesoscale ice generalized effective size, defined 
!                  as in Fu  (1996, J. Clim.) (micrometers)
!   dmeml          mass flux in mesoscale downdraft [ kg/((m**2) s) ]
!                  (normalized by a(1,p_b)) 
!   ecds           normalized convective downdraft evaporation
!                  [ kg(H2O)/kg/sec ]
!   eces           normalzed convective-updraft evporation/sublimation
!                  [ kg(H2O)/kg/sec ]
!   elt            normalized melting [ K/sec ]
!   emds           normalized mesoscale-downdraft sublimation
!                  [ kg(H2O)/kg/sec ]
!   emes           normalized mesoscale-updraft sublimation
!                  [ kg(H2O)/kg/sec ]
!   fre            normalized freezing [ K/sec ]
!   mrmes           normalized mesoscale moisture-flux convergence
!                  [ kg(H2O)/kg/sec ]
!   tmes           normalized mesoscale entropy-flux convergence
!                  [ K/sec ]
!                  Entropy-flux convergence is mesoscale component
!                  of second term in expression for cumulus thermal
!                  forcing in Fig. 3 of Donner (1993, JAS).
!   uceml          normalized mass fluxes in cell updrafts
!                  [ kg/((m**2)*s ] 
!   umeml          mass flux in mesoscale updraft [ kg/((m**2) s) ]
!                  (normalized by a(1,p_b)) 
!   wmms           normalized mesoscale deposition of water vapor from
!                  cells [ kg(H2O)/kg/sec ]
!   wmps           normalized mesoscale redistribution of water vapor
!                  from cells [ kg(H2O)/kg/sec ]
!   xice           mesoscale ice mass mixing ratio (kg(ice)/kg)
!   xliq           mesoscale liquid mass mixing ratio (kg(ice)/kg)
!   qtceme         tracer tendencies due to donner_deep_mod 
!                  [ kg/kg/s ]
!   qtmes1         tracer time tendency due to mesoscale  motions
!                  [ kg/kg/s ]
!   qtren1         tracer time tendency due to cell scale motions
!                  [ kg/kg/s ]
!   wtp1           redistribution of tracer from cellscale to mesoscale
!                  [ kg/kg/s ]
!   a1             fractional area of index-1 cu subensemble
!   amax           maximum value for a_1(p_b)
!                  See "a Bounds 6/7/97" notes
!   amos           upper limit on cloud fractional area based on
!                  moisture constraint See "Moisture Constraint," 
!                  8/8/97.
!   ampta1         area weighted mesoscale cloud fraction, normal-
!                  ized by a(1,p_b)
!   cell_precip    area weighted convective precipitation rate
!                  [ mm/day ]
!   dcape          time rate of change of cape [ J/(kg s) ]
!   emdi_v         vertical integral of mesoscale-downdraft 
!                  sublimation
!   meso_precip    area weighted mesoscale precipitation rate
!                  [ mm/day ]
!   pb_v           pressure at base of cumulus updrafts (Pa)
!   pmd_v          pressure at top of mesoscale downdraft (Pa)
!   przm           pressure at base of anvil (based on presence of ice)
!   prztm          pressure at top of anvil (based on presence of ice)
!   pzm_v          pressure at base of mesoscale updraft (Pa)
!   pztm_v         pressure at top of mesoscale updraft (Pa)


real, dimension(:,:,:), pointer  ::           &
                 cecon=>NULL(),                    &
                 ceefc=>NULL(),                    &
                 cell_ice_geneff_diam=>NULL(),     &
                 cell_liquid_eff_diam=>NULL(),     &
                 cememf_mod=>NULL(),               &
                 cemfc=>NULL(),                    &
                 cmus=>NULL(),                     &
                 conv_temp_forcing=>NULL(),        &
                 conv_moist_forcing=>NULL(),       &
                 cual=>NULL(),                     &
                 cuqi=>NULL(),                     &
                 cuql=>NULL(),                     &
                 detmfl=>NULL(),                   &
                 dgeice=>NULL(),                   &
                 dmeml=>NULL(),                    &
                 ecds=>NULL(),                     &
                 eces=>NULL(),                     &
                 elt=>NULL(),                      &
                 emds=>NULL(),                     &
                 emes=>NULL(),                     &
                 fre=>NULL(),                      &        
                 mrmes=>NULL(),                    &
                 tmes=>NULL(),                     &
                 uceml=>NULL(),                    &
                 umeml=>NULL(),                    &
                 wmms=>NULL(),                     &
                 wmps=>NULL(),                     &
                 xice=>NULL(),                     &
                 xliq=>NULL()
real, dimension(:,:,:,:), pointer ::          &
                 qtceme=>NULL(),                   &
                 qtmes1=>NULL(),                   &
                 qtren1=>NULL(),                   &
                 wtp1=>NULL()
real, dimension(:,:),   pointer  ::           &
                 a1=>NULL(),                       &
                 amax=>NULL(),                     &
                 amos=>NULL(),                     &
                 ampta1=>NULL(),                   &
                 cell_precip=>NULL(),              &
                 dcape=>NULL(),                    &
                 emdi_v=>NULL(),                   &
                 meso_precip=>NULL(),              &
                 pb_v=>NULL(),                     &
                 pmd_v=>NULL(),                    &
                 przm=>NULL(),                     &
                 prztm=>NULL(),                    &
                 pzm_v=>NULL(),                    &
                 pztm_v=>NULL()

end type donner_conv_type


!#######################################################################

type donner_cape_type

!***********************************************************************
!    these variables are allocated and deallocated on each donner step.
!    they contain variables related to the lifting of a parcel in a 
!    convective environment.
!**********************************************************************
 
!   coin           convective inhibition 
!                  energy required to lift parcel from level istart
!                  to level of free convection. [ J/kg ]
!   plcl           pressure at lifting condensation level [ Pa ]
!   plfc           pressure at level of free convection [ Pa ]
!                  height of plfc .le. height of plcl. if parcel 
!                  becomes buoyant below plcl, cin can be .lt. 0
!   plzb           pressure at level of zero buoyancy [ Pa ]
!   qint           vertically integrated column moisture
!                  [ kg (h20)/(m**2) ]
!   qint_lag       vertically integrated column moisture, calculated
!                  using the lag timestep moisture profile.
!                  [ kg (h20)/(m**2) ]
!   xcape          convective available potential energy. energy 
!                  released as parcel moves from level of free 
!                  convection to level of zero buoyancy, calculated on 
!                  convection step [ J / kg ].
!   xcape_lag      convective available potential energy, calculated
!                  using  the lag timestep sounding. [ J / kg ]

!   the following variables are on the enhanced cape vertical grid
!   with k index 1 closest to the ground:

!   cape_p         pressure levels of cape grid       [ hPa ]
!   env_r          environmental mixing ratio profile [ kg/kg ]
!   env_t          environmental temperature profile  [ K ]
!   parcel_r       parcel mixing ratio                [ kg/kg ]
!   parcel_t       parcel temperature                 [ K ]

!   the following variables are on large-scale model levels:

!   model_p        pressure profile used to define pressure 
!                  profile used in cape calculation [ hPa ]
!   model_r        mixing ratio profile used to define
!                  moisture profile used in cape calculation [ kg/kg ]
!   model_t        temperature profile used to define
!                  temperature profile used in cape calculation [ K ]


real, dimension(:,:), pointer ::          &
                         coin=>NULL(),         &
                         plcl=>NULL(),         &
                         plfc=>NULL(),         &
                         plzb=>NULL(),         &
                         qint=>NULL(),         &
                         qint_lag=>NULL(),     &
                         xcape_lag=>NULL(),    &
                         xcape=>NULL()
real, dimension(:,:,:), pointer ::        &
                         cape_p=>NULL(),       &
                         env_r=>NULL(),        &
                         env_t=>NULL(),        &
                         parcel_r=>NULL(),     &
                         parcel_t=>NULL()
real, dimension (:,:,:), pointer ::       &
                         model_p=>NULL(),      &
                         model_r=>NULL(),      &
                         model_t=>NULL()

end type donner_cape_type


!######################################################################
