
!VERSION NUMBER:
MODULE MESSY_CONVECT_DONNER_DEEP_K
!  $Id: donner_deep_k.F90,v 13.0.4.1.2.1 2006/10/28 15:58:01 rsh Exp $

!module donner_deep_inter_mod

!#include "donner_deep_interfaces.h"

!end module donner_deep_inter_mod

CONTAINS
!######################################################################

subroutine don_d_donner_deep_k   &
         (is, ie, js, je, isize, jsize, nlev_lsm, nlev_hires, ntr, me, &
          cloud_tracers_present,  &
          dt, Param, Nml, temp, mixing_ratio, pfull, phalf, omega, qlin,&
          qiin, qain, land, sfc_sh_flux, sfc_vapor_flux, tr_flux,  &
          tracers, cell_cld_frac, cell_liq_amt, cell_liq_size,  &
          cell_ice_amt, cell_ice_size, meso_cld_frac, meso_liq_amt,  &
          meso_liq_size, meso_ice_amt, meso_ice_size, nsum, & 
          precip, delta_temp, delta_vapor, detf, uceml_inter, mtot,   &
          donner_humidity_area, donner_humidity_ratio, total_precip,  &
          temperature_forcing, moisture_forcing, parcel_rise, &
          delta_ql, delta_qi, delta_qa, qtrtnd, calc_conv_on_this_step, &
          ermesg, Initialized, Col_diag, Don_rad, Don_conv, Don_cape, &
          Don_save, exit_flag)      
                        
!-------------------------------------------------------------------
!    subroutine don_d_donner_deep_k is the primary kernel sub-
!    routine of the donner deep convection parameterization. it receives
!    all input needed from donner_deep_mod and controls the generation
!    of output that is returned to donner_deep_mod, from which it is made
!    accessible to the rest of the model parameterizations, as needed.
!-------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_initialized_type, donner_save_type,&
                             donner_rad_type, donner_nml_type, &
                             donner_param_type, donner_conv_type, &
                             donner_column_diag_type, donner_cape_type

implicit none

!--------------------------------------------------------------------
integer,                 intent(in)     :: is, ie, js, je, isize, jsize,&
                                           nlev_lsm, nlev_hires, ntr, me
logical,                 intent(in)     :: cloud_tracers_present
real,                    intent(in)     :: dt
type(donner_param_type), intent(in)     :: Param
type(donner_nml_type),   intent(in)     :: Nml
real, dimension(isize,jsize,nlev_lsm),                                  &
                         intent(in)     :: temp, mixing_ratio, pfull,   &
                                           omega, qlin, qiin, qain, &
                                           cell_cld_frac,  cell_liq_amt,&
                                           cell_liq_size, cell_ice_amt, &
                                           cell_ice_size, meso_cld_frac,&
                                           meso_liq_amt, meso_liq_size, &
                                           meso_ice_amt, meso_ice_size
real,    dimension(isize,jsize,nlev_lsm+1),                            &
                         intent(in)     :: phalf 
real,    dimension(isize,jsize),                                      &
                         intent(in)     :: land, sfc_sh_flux,       &
                                           sfc_vapor_flux
real,    dimension(isize,jsize,ntr),                                 &
                         intent(in)     :: tr_flux 
real,    dimension(isize,jsize,nlev_lsm,ntr),                         &
                         intent(in)     :: tracers 
integer, dimension(isize,jsize),                                     &
                         intent(in)     :: nsum      
real,    dimension(isize,jsize),                                       & 
                         intent(out)    :: precip      
real, dimension(isize,jsize,nlev_lsm),                                 &
                         intent(out)    :: delta_temp, delta_vapor,&
                                           detf,  mtot, &
                                           donner_humidity_area,&
                                           donner_humidity_ratio, &
                                           temperature_forcing,   &
                                           moisture_forcing, &
                                           delta_ql, delta_qi, &
                                           delta_qa
real, dimension(isize,jsize,nlev_lsm+1),                               &
                         intent(out)    :: uceml_inter
real, dimension(isize,jsize),                                         &
                         intent(out)    :: total_precip, parcel_rise
real,    dimension(isize,jsize,nlev_lsm,ntr),                        &
                         intent(out)    :: qtrtnd 
logical,                 intent(out)    :: calc_conv_on_this_step
character(len=*),        intent(out)    :: ermesg
type(donner_initialized_type),                            &
                         intent(inout)  :: Initialized
type(donner_column_diag_type),                               &
                         intent(inout)  :: Col_diag
type(donner_rad_type),   intent(inout)  :: Don_rad
type(donner_conv_type),  intent(inout)  :: Don_conv
type(donner_cape_type),  intent(inout)  :: Don_cape
type(donner_save_type),  intent(inout)  :: Don_save
!!$type(donner_column_diag_type)                               &
!!$                           :: Col_diag
!!$type(donner_rad_type)     :: Don_rad
!!$type(donner_conv_type)    :: Don_conv
!!$type(donner_cape_type)    :: Don_cape
!!$type(donner_save_type)    :: Don_save

!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     dt             physics time step [ sec ]
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   vapor mixing ratio field at model levels 
!                    [ kg(h20) / kg(dry air) ]
!     pfull          pressure field on model full levels [ Pa ]
!     omega          model omega field at model full levels [ Pa / sec ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qain           large-scale cloud fraction  
!                    [ fraction ]
!     cell_cld_frac  fractional coverage of convective cells in
!                    grid box [ dimensionless ]
!     cell_liq_amt   liquid water content of convective cells
!                    [ kg(h2o) / kg(air) ]
!     cell_liq_size  assumed effective size of cell liquid drops
!                    [ microns ]
!     cell_ice_amt   ice water content of cells
!                    [ kg(h2o) / kg(air) ]
!     cell_ice_size  generalized effective diameter for ice in
!                    convective cells [ microns ]
!     meso_cld_frac  fractional area of mesoscale clouds in grid
!                    box [ dimensionless ]
!     meso_liq_amt   liquid water content in mesoscale clouds
!                    [ kg(h2o) / kg(air) ]
!     meso_liq_size  assumed effective size of mesoscale drops
!                    [ microns ]
!     meso_ice_amt   ice water content of mesoscale elements
!                    [ kg(h2o) / kg(air) ]
!     meso_ice_size  generalized ice effective size for anvil ice
!                    [ microns ]
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     land           fraction of grid box covered by land
!                    [ fraction ]
!     sfc_sh_flux   sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     tr_flux        surface flux of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!     nsum           number of time levels over which the above variables
!                    have so far been summed
!
!   intent(out) variables:
!
!     precip         precipitation generated by deep convection
!                    [ kg(h2o) / m**2 ]
!     delta_temp     temperature increment due to deep convection 
!                    [ deg K ]
!     delta_vapor    water vapor mixing ratio increment due to deep 
!                    convection [ kg(h2o) / kg (dry air) ]
!     detf           detrained cell mass flux at model levels 
!                    [ (kg / (m**2 sec) ) ]
!     mtot           mass flux at model full levels, convective plus 
!                    mesoscale, due to donner_deep_mod 
!                    [ (kg / (m**2 sec) ) ]
!     donner_humidity_area
!                    fraction of grid box in which humidity is affected
!                    by the deep convection, defined as 0.0 below cloud
!                    base and above the mesoscale updraft, and as the
!                    sum of the cell and mesoscale cloud areas in 
!                    between. it is used in strat_cloud_mod to determine
!                    the large-scale specific humidity field for the
!                    grid box. DO NOT use for radiation calculation,
!                    since not all of this area includes condensate.
!                    [ fraction ]
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to specific 
!                    humidity in environment outside convective system
!                    [ dimensionless ]
!     temperature_forcing  
!                    temperature tendency due to donner convection
!                    [ deg K / sec ]
!     moisture_forcing 
!                    vapor mixing ratio tendency due to donner
!                    convection [ kg(h2o) / (kg(dry air) sec ) ]
!     delta_ql       cloud water specific humidity increment due to 
!                    deep convection over the timestep
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qi       cloud ice specific humidity increment due to deep 
!                    convection over the timestep 
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qa       cloud area increment due to deep convection
!                    over the time step [ fraction ]
!     uceml_inter    upward cell mass flux at interface levels 
!                    [ (kg / (m**2 sec) ) ]
!     total_precip   total precipitation rate produced by the
!                    donner parameterization [ mm / day ]
!     parcel_rise    accumulated vertical displacement of a
!                    near-surface parcel as a result of the lowest
!                    model level omega field [ Pa ]
!     qtrtnd         tracer time tendencies due to deep convection
!                    during the time step
!                    [ kg(tracer) / (kg (dry air) sec) ]
!     calc_conv_on_this_step
!                    is this a step on which to calculate
!                    donner convection ?
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!   intent(inout) variables:
!
!     Initialized    donner_initialized_type variable containing
!                    variables which are defiuned during initialization.
!                    these values may be changed during model execution.
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     Don_rad        donner_rad_type derived type variable used to hold 
!                    those fields needed to connect the donner deep 
!                    convection parameterization and the model radiation 
!                    package
!     Don_conv       donner_convection_type derived type variable
!                    containing diagnostics and intermediate results 
!                    describing the nature of the convection produced by
!                    the donner parameterization
!     Don_cape       donner_cape type derived type variable containing 
!                    diagnostics and intermediate results related to the
!                    cape calculation associated with the donner convec-
!                    tion parameterization
!     Don_save       donner_save_type derived type variable containing
!                    those variables which must be preserved across
!                    timesteps
!     
!-----------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real,                                                      &
         dimension (isize, jsize, nlev_lsm) ::  lag_cape_temp,        &
                                                lag_cape_vapor,      &
                                                lag_cape_press, &
                                                dql, dqi, dqa
      real,    dimension (isize, jsize)     ::  current_displ
      logical, dimension (isize, jsize)     ::  exit_flag
      integer                               ::  idiag, jdiag, unitdiag

      integer  :: i, j, k, n   

!--------------------------------------------------------------------
!   local variables:
!
!     lag_cape_temp        temperature field used in lag-time cape 
!                          calculation [ deg K ]
!     lag_cape_vapor       vapor mixing ratio field used in lag-time
!                          cape calculation [ kg(h2o) / kg(dry air) ]
!     lag_cape_press       model full-level pressure field used in 
!                          lag-time cape calculation 
!                          [ kg(h2o) / kg(dry air) ]
!     dql                  tendency of cloud liquid specific humidity
!                          due to donner convection 
!                          [ kg(h2o) / kg(moist air) / sec ]
!     dqi                  tendency of cloud ice specific humidity
!                          due to donner convection 
!                          [ kg(h2o) / kg(moist air) / sec ]
!     dqa                  tendency of large-scale cloud area
!                          due to donner convection 
!                          [ fraction / sec ]
!     current_displ        low-level parcel displacement to use in cape
!                          calculation on this step [ Pa ]
!     exit_flag            logical array indicating whether deep conv-
!                          ection exists in a column
!     idiag                physics window i index of current diagnostic
!                          column
!     jdiag                physics window j index of current diagnostic
!                          column
!     unitdiag             i/o unit assigned to current diagnostic
!                          column
!     i, j, k, n           do-loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' '
!-------------------------------------------------------------------
!    verify that donner_deep_freq is an integral multiple of the model
!    physics time step.
!--------------------------------------------------------------------
      if (MOD (Nml%donner_deep_freq, int(dt)) /= 0) then
        ermesg = 'donner_deep_donner_deep_k: donner_deep timestep NOT &
            &an integral multiple of physics timestep'
        print*, ermesg
        return
      endif

!--------------------------------------------------------------------
!    decrement the time remaining before the convection calculations 
!    on the first entry to this routine on a given timestep. save the
!    current model physics timestep.
!--------------------------------------------------------------------
      if (Initialized%pts_processed_conv == 0) then
        Initialized%conv_alarm  = Initialized%conv_alarm - int(dt)
        Initialized%physics_dt = int(dt)
      endif
!--------------------------------------------------------------------
!    set a flag to indicate whether the convection calculation is to be 
!    done on this timestep. if this is the first call to donner_deep 
!    (i.e., coldstart), convection cannot be calculated because the 
!    lag profiles needed to calculate cape are unavailable, and so
!    a time tendency of cape can not be obtained. otherwise, it is a
!    calculation step or not dependent on whether the convection "alarm"
!    has gone off. 
!---------------------------------------------------------------------
      if (Initialized%coldstart) then
        calc_conv_on_this_step = .false.
      else
        if (Initialized%conv_alarm <= 0) then
          calc_conv_on_this_step = .true.
        else
          calc_conv_on_this_step = .false.
        endif
      endif
!---------------------------------------------------------------------
!    perform the following calculations only if this is a step upon
!    which donner convection is to be calculated.
!---------------------------------------------------------------------
      if (calc_conv_on_this_step) then
 
!-------------------------------------------------------------------
!    call initialize_local_variables_k to allocate and initialize the
!    elements of the donner_conv, donner_cape and donner_rad derived type
!    variables.
!-------------------------------------------------------------------
        call don_d_init_loc_vars_k       &
             (isize, jsize, nlev_lsm, ntr, nlev_hires, cell_cld_frac, &
              cell_liq_amt, cell_liq_size, cell_ice_amt, cell_ice_size, &
              meso_cld_frac, meso_liq_amt, meso_liq_size, meso_ice_amt, &
              meso_ice_size, nsum, Don_conv, Don_cape, Don_rad, ermesg) 

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return
      endif !(calc_conv_on_this_step)
!---------------------------------------------------------------------
!    add the vertical displacement resulting from the current omega
!    field to the accumulated displacement of a parcel which originated
!    at the lowest model full level. prevent the parcel from moving 
!    below its starting point or going out the top of the atmosphere. 
!---------------------------------------------------------------------
      do j=1,jsize       
        do i=1,isize        
          parcel_rise(i,j) = Don_save%parcel_disp(i+is-1,j+js-1) +  &
                             omega(i,j,nlev_lsm)*dt
          parcel_rise(i,j) = MIN (0.0, parcel_rise(i,j))
          parcel_rise(i,j) = MAX (-phalf(i,j,nlev_lsm+1),    &
                                  parcel_rise(i,j))
        end do
      end do
!---------------------------------------------------------------------
!    if there are one or more diagnostic columns in the current physics
!    window, set a flag to so indicate. call don_d_column_input_fields
!    to print out the relevant input fields, location and control 
!    information for these diagnostics columns.   
!---------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then
        if (Col_diag%ncols_in_window > 0) then
          Col_diag%in_diagnostics_window = .true.
          call don_d_column_input_fields_k    &
               (isize, jsize, nlev_lsm, dt, calc_conv_on_this_step, &
                Col_diag, temp, mixing_ratio, pfull, omega, phalf,   &
                parcel_rise, ermesg) 

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    if there are no diagnostic columns in the current physics
!    window, set a flag to so indicate. 
!---------------------------------------------------------------------
        else
          Col_diag%in_diagnostics_window = .false.
        endif

!---------------------------------------------------------------------
!    if column diagnostics are not desired in any model columns, set a 
!    flag to so indicate. 
!---------------------------------------------------------------------
      else
        Col_diag%in_diagnostics_window = .false.
      endif

!---------------------------------------------------------------------
!    perform the following calculations only if this is a step upon 
!    which donner convection is to be calculated.
!---------------------------------------------------------------------
      if (calc_conv_on_this_step) then 

!---------------------------------------------------------------------
!    define the low-level displacement to be used on this step 
!    (current_displ). it is the current time-integrated value, unless 
!    the current lowest-level vertical velocity is downward, in which 
!    case the displacement to be used on the current step is set to 
!    zero, precluding deep convection on this step.
!---------------------------------------------------------------------
        do j=1,jsize       
          do i=1,isize        
            if (omega(i,j,nlev_lsm) > 0.)   then
              current_displ(i,j) = 0. 
            else
              current_displ(i,j) = parcel_rise(i,j)
            endif
          end do
        end do

!---------------------------------------------------------------------
!    define the temperature, vapor mixing ratio and pressure fields to
!    be used in calculating the lag values of cape so that a cape tend-
!    ency due to large-scale forcing may be computed.
!---------------------------------------------------------------------
        lag_cape_temp (:,:,:) = Don_save%lag_temp (is:ie,js:je,:)
        lag_cape_vapor(:,:,:) = Don_save%lag_vapor(is:ie,js:je,:)
        lag_cape_press(:,:,:) = Don_save%lag_press(is:ie,js:je,:)

!---------------------------------------------------------------------
!    call donner_convection_driver to calculate the effects of deep 
!    convection.  
!--------------------------------------------------------------------- 
        call don_d_convection_driver_k   &
             (isize, jsize, nlev_lsm, nlev_hires, ntr, me, &
              cloud_tracers_present, dt, Nml, &
              Initialized, Param, Col_diag, temp, mixing_ratio, pfull, &
              qlin, qiin, qain, lag_cape_temp, lag_cape_vapor,    &
              lag_cape_press, phalf, current_displ, land, sfc_sh_flux, &
              sfc_vapor_flux, tr_flux, tracers, Don_cape, Don_conv, &
              Don_rad, temperature_forcing, moisture_forcing,  &
              total_precip, donner_humidity_ratio, donner_humidity_area,&
              dql, dqi, dqa, exit_flag, ermesg) 
!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return
!--------------------------------------------------------------------- 
!    define the module variables used to preserve output fields across
!    physics timesteps, needed when the donner parameterization is not
!    executed on every physics step.
!
!    1) the variables defining the humidity disturbance caused by
!       donner convection. these variables are needed so that the large-
!       scale humidity may be properly adjusted in strat_cloud_mod.
!--------------------------------------------------------------------- 
        Don_save%humidity_area (is:ie,js:je,:) =    &
                                            donner_humidity_area (:,:,:)
        Don_save%humidity_ratio(is:ie,js:je,:) =    &
                                            donner_humidity_ratio(:,:,:)
!--------------------------------------------------------------------- 
!    2) the total precipitation produced by the donner parameter-
!       ization.
!--------------------------------------------------------------------- 
        Don_save%tprea1(is:ie,js:je) = total_precip(:,:)

!---------------------------------------------------------------------
!    3) the vapor and temperature forcing resulting from the donner
!       deep parameterization that will be output to the calling 
!       routine. NOTE: these values of cemetf and cememf have the terms
!       related to flux convergence of condensate and mesoscale 
!       detrainment removed when parameterization is run in a model 
!       using strat_cloud_mod, since in that case those terms will be 
!       calculated within that module.
!---------------------------------------------------------------------
        Don_save%cememf(is:ie,js:je,:) = moisture_forcing(:,:,:) 
        Don_save%cemetf(is:ie,js:je,:) = temperature_forcing(:,:,:)

!----------------------------------------------------------------------
!    4) the increments which must be applied to the strat_cloud vari-
!       ables as a result of donner convection. the returned tendencies 
!       are converted to increments. 
!----------------------------------------------------------------------
        if (cloud_tracers_present) then
        Don_save%dql_strat(is:ie,js:je,:) = dql(:,:,:)*dt
        Don_save%dqi_strat(is:ie,js:je,:) = dqi(:,:,:)*dt
        Don_save%dqa_strat(is:ie,js:je,:) = dqa(:,:,:)*dt
        endif
!--------------------------------------------------------------------
!    5) the net mass flux and detrained cell mass flux at model full 
!       levels that is associated with the donner deep parameterization. 
!       the net mass flux is needed as input to strat_cloud_mod, while
!       the detrained cell mass flux is needed by cu_mo_trans_mod. 
!--------------------------------------------------------------------
        do k=1,nlev_lsm
          do j=1,jsize
            do i=1,isize
              if ((Don_conv%uceml(i,j,k) <= 1.0e-10) .and.   &
                  (Don_conv%umeml(i,j,k) <= 1.0e-10) .and.   &
                  (Don_conv%dmeml(i,j,k) <= 1.0e-10) ) then
                Don_save%mass_flux(i+is-1,j+js-1,k) = 0.
              else
                Don_save%mass_flux(i+is-1,j+js-1,k) =   &
                                        Don_conv%uceml(i,j,k) + &
                                        Don_conv%umeml(i,j,k) + &
                                        Don_conv%dmeml(i,j,k)
              endif
              if (Don_conv%detmfl(i,j,k) <= 1.0e-10) then
                Don_save%det_mass_flux(i+is-1,j+js-1,k) = 0.
              else
                Don_save%det_mass_flux(i+is-1,j+js-1,k) =      &
                                                  Don_conv%detmfl(i,j,k)
              endif
            end do
          end do
        end do
!--------------------------------------------------------------------
!    6) the upward mass flux at model interface levels that is 
!       associated with the convective cells present in the donner deep 
!       convction parameterization. this is needed by cu_mo_trans_mod.
!       define values at upper and lower boundary to be 0.0.
!--------------------------------------------------------------------
        Don_save%cell_up_mass_flux(:,:,1) = 0.
        Don_save%cell_up_mass_flux(:,:,nlev_lsm+1) = 0.
        do k=2,nlev_lsm
          do j=1,jsize
            do i=1,isize
              Don_save%cell_up_mass_flux(i+is-1,j+js-1,k) =  &
                                    0.5*(Don_Conv%uceml(i,j,k) + &
                                         Don_conv%uceml(i,j,k-1))
            end do
          end do
        end do

!--------------------------------------------------------------------
!    7) the tracer time tendencies resulting from the donner param-
!       eterization. 
!--------------------------------------------------------------------
        if (Initialized%do_donner_tracer) then
          Don_save%tracer_tends(is:ie,js:je,:,:) =      &
                                               Don_conv%qtceme(:,:,:,:)
        endif
!---------------------------------------------------------------------
!    end of if loop for code executed only on steps for which the donner
!    parameterization is requested.
!---------------------------------------------------------------------
      endif ! (calc_conv_on_this_step)
!--------------------------------------------------------------------- 
!    update the module variable containing the total low-level parcel 
!    displacement. this field is updated on every model physics step, 
!    regardless of whether donner convective tendencies are calculated 
!    or not.
!--------------------------------------------------------------------- 
      Don_save%parcel_disp(is:ie,js:je) = parcel_rise(:,:)
!----------------------------------------------------------------------
!    define the output fields to be passed back to the calling routine.
!    these fields are returned on every physics step, regardless of
!    whether or not the donner calculation is done, and so must be 
!    defined from module variables.
!
!    1) the temperature increment due to deep convection
!----------------------------------------------------------------------
      delta_temp(:,:,:) = Don_save%cemetf(is:ie,js:je,:)*dt
!----------------------------------------------------------------------
!    2) the moisture increment due to deep convection. if the moisture
!       tendency results in the production of a negative value, reduce 
!       the tendency to avoid producing the negative mixing ratio.
!----------------------------------------------------------------------
      do k=1,nlev_lsm
        do j=1,jsize      
          do i=1,isize       
            delta_vapor(i,j,k) = Don_save%cememf(i+is-1,j+js-1,k)*dt
            if ((mixing_ratio(i,j,k) + delta_vapor(i,j,k)) < 0.) then
              if (mixing_ratio(i,j,k) > 0.) then
                delta_vapor (i,j,k) = -mixing_ratio(i,j,k)
              else 
                delta_vapor(i,j,k) = 0.0
              endif
            endif
          end do
        end do
      end do
!-------------------------------------------------------------------
!    3) the net mass flux, detrained cell mass flux and upward mass 
!       flux due to convective cells at interface levels resulting from
!       donner convection.
!-------------------------------------------------------------------
      mtot(:,:,:)        = Don_save%mass_flux(is:ie, js:je,:)
      detf(:,:,:)        = Don_save%det_mass_flux(is:ie,js:je,:)
      uceml_inter(:,:,:) = Don_save%cell_up_mass_flux(is:ie,js:je,:)
!-------------------------------------------------------------------
!    4) the increments of the large-scale cloud variables due to deep
!       convection and the variables describing the specific humidity
!       disturbance associated with donner convection.
!-------------------------------------------------------------------
      if (cloud_tracers_present) then
      delta_ql(:,:,:) = Don_save%dql_strat (is:ie, js:je,:)
      delta_qi(:,:,:) = Don_save%dqi_strat (is:ie, js:je,:)
      delta_qa(:,:,:) = Don_save%dqa_strat (is:ie, js:je,:)
      endif
      donner_humidity_area(:,:,:)  =             &
                                   Don_save%humidity_area(is:ie,js:je,:)
      donner_humidity_ratio(:,:,:) =             &
                                   Don_save%humidity_ratio(is:ie,js:je,:)
!----------------------------------------------------------------------
!    5) the precipitation accrued on the current timestep from deep
!       convection. 
!       note: precip    [mm/day] * 1/86400 [day/sec] * 1/1000 [ m/mm] * 
!                  1000 [kg(h2o)/m**3] * dt [sec] = kg/m2, as desired. 
!----------------------------------------------------------------------
      precip(:,:) = Don_save%tprea1(is:ie,js:je)*dt/Param%seconds_per_day
!--------------------------------------------------------------------
!    6) time tendencies of any tracers being transported by donner 
!       convection. if none have been defined, fill the output array
!       with zeroes.
!--------------------------------------------------------------------
      if (Initialized%do_donner_tracer) then
        qtrtnd(:,:,:,:) = Don_save%tracer_tends(is:ie,js:je,:,:)
      else
        qtrtnd(:,:,:,:) = 0.0                            
      endif
!---------------------------------------------------------------------
!    if this is a diagnostics window, output the increments to temper-
!    ature and vapor mixing ratio at levels where donner convection
!    has produced a modification.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          idiag = Col_diag%i_dc(n)
          jdiag = Col_diag%j_dc(n)
          unitdiag = Col_diag%unit_dc(n)
          do k=Col_diag%kstart,nlev_lsm
            if (delta_temp(idiag,jdiag,k) /= 0.0) then
              write (unitdiag, '(a, i4, f20.14, e20.12)') &
                   'in donner_deep: k,ttnd,qtnd',  k,  &
                   delta_temp(idiag,jdiag,k), delta_vapor(idiag,jdiag,k)
            endif
          end do
        end do
      endif
!---------------------------------------------------------------------
!    define the module variables containing the temperature, pressure 
!    and vapor fields that are to be used on the next time step to 
!    calculate a lag-time cape so that the time tendency of cape due
!    to large-scale forcing may be obtained. this field is updated on 
!    every model physics step, so that values are present in case the 
!    next step is a donner calculation step.
!--------------------------------------------------------------------- 
      Don_save%lag_temp (is:ie, js:je,:) = temp + delta_temp
      Don_save%lag_vapor(is:ie, js:je,:) = mixing_ratio + delta_vapor 
      Don_save%lag_press(is:ie, js:je,:) = pfull
!-------------------------------------------------------------------
!    perform the following calculations only if this is a step upon 
!    which donner convection is to be calculated.
!-------------------------------------------------------------------
      if (calc_conv_on_this_step) then
!---------------------------------------------------------------------
!    define the revised moisture tendency produced by donner convection
!    after it has been adjusted to prevent the generation of negative 
!    vapor mixing ratio.
!---------------------------------------------------------------------
        Don_conv%cememf_mod(:,:,:) = delta_vapor(:,:,:)/dt
!---------------------------------------------------------------------
!    if this is a diagnostics window, call donner_column_end_of_step
!    to output various diagnostic fields in the specified diagnostic
!    columns.
!---------------------------------------------------------------------
        if (Col_diag%in_diagnostics_window) then
          call don_d_column_end_of_step_k   &
               (isize, jsize, nlev_lsm, ntr, Col_diag, exit_flag,   &
                total_precip, parcel_rise, temperature_forcing,   &
                moisture_forcing, tracers(:,:,:,:), Don_cape,   &
                Don_conv, ermesg) 

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (trim(ermesg) /= ' ') return
        endif
      endif !(calc_conv_on_this_step) 

!---------------------------------------------------------------------
!    determine if all points on this processor's subdomain have been 
!    processed on this timestep. if so, set the point counter to zero, 
!    and if this was the first time through the parameterization, set
!    the flag so indicating (coldstart) to be .false.. if all points 
!    have been processed and if this was a calculation step, set the 
!    alarm to define the next time at which donner convection is to be 
!    executed.
!----------------------------------------------------------------------
      Initialized%pts_processed_conv = Initialized%pts_processed_conv + &
                                       isize*jsize
      if (Initialized%pts_processed_conv >= Initialized%total_pts) then
        Initialized%pts_processed_conv = 0 
        if (Initialized%coldstart) Initialized%coldstart = .false.
        if (calc_conv_on_this_step) then
          Initialized%conv_alarm = Initialized%conv_alarm +    &
                                   Nml%donner_deep_freq 
        endif 
      endif
!--------------------------------------------------------------------
      RETURN
      
end subroutine don_d_donner_deep_k




!####################################################################

subroutine don_d_init_loc_vars_k      &
         (isize, jsize, nlev_lsm, ntr, nlev_hires, cell_cld_frac,  &
          cell_liq_amt, cell_liq_size, cell_ice_amt, cell_ice_size,   &
          meso_cld_frac, meso_liq_amt, meso_liq_size, meso_ice_amt,   &
          meso_ice_size, nsum, Don_conv, Don_cape, Don_rad, ermesg)     

use messy_convect_donner_types_mod, only : donner_rad_type, donner_conv_type, &
                             donner_cape_type
implicit none

!--------------------------------------------------------------------
!   subroutine don_d_init_loc_vars_k allocates space 
!   for and initializes the array components of the donner_conv_type 
!   variable Don_conv, the donner_cape_type variable Don_cape, and the 
!   donner_rad_type variable Don_rad.
!--------------------------------------------------------------------

integer,                         intent(in)    :: isize, jsize,    &
                                                  nlev_lsm, ntr, &
                                                  nlev_hires
real,dimension(isize,jsize,nlev_lsm),                              &
                                 intent(in)    :: cell_cld_frac,  &
                                                  cell_liq_amt,   &
                                                  cell_liq_size, &
                                                  cell_ice_amt,   &
                                                  cell_ice_size, &
                                                  meso_cld_frac,    &
                                                  meso_liq_amt, &
                                                  meso_liq_size, &
                                                  meso_ice_amt,     &
                                                  meso_ice_size 
integer, dimension(isize,jsize), intent(in)    :: nsum
type(donner_conv_type),          intent(inout) :: Don_conv
type(donner_cape_type),          intent(inout) :: Don_cape
type(donner_rad_type),           intent(inout) :: Don_rad
character(len=*),                intent(out)   :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize         size of x-dimension of physics window
!     jsize         size of y-dimension of physics window
!     nlev_lsm      number of layers in large-scale model
!     ntr           number of tracers to be transported by donner
!                   convection
!     nlev_hires    number of model layers in hi-res cloud model
!                   of the donner deep convection parameterization
!     cell_cld_frac fractional coverage of convective cells in
!                   grid box [ dimensionless ]
!     cell_liq_amt  liquid water content of convective cells
!                   [ kg(h2o) / kg(air) ]
!     cell_liq_size assumed effective size of cell liquid drops
!                   [ microns ]
!     cell_ice_amt  ice water content of cells
!                   [ kg(h2o) / kg(air) ]
!     cell_ice_size generalized effective diameter for ice in
!                   convective cells [ microns ]
!     meso_cld_frac fractional area of mesoscale clouds in grid
!                   box [ dimensionless ]
!     meso_liq_amt  liquid water content in mesoscale clouds
!                   [ kg(h2o) / kg(air) ]
!     meso_liq_size assumed effective size of mesoscale drops
!                   [ microns ]
!     meso_ice_amt  ice water content of mesoscale elements
!                   [ kg(h2o) / kg(air) ]
!     meso_ice_size generalized ice effective size for anvil ice
!                   [ microns ]
!     nsum          number of time levels over which the above variables
!                   have so far been summed
!
!   intent(inout) variables:
!
!     Don_conv     donner_conv_type derived type variable containing 
!                  diagnostics and intermediate results describing the 
!                  nature of the convection produced by the donner 
!                  parameterization
!     Don_cape     donner_cape type derived type variable containing 
!                  diagnostics and intermediate results related to the 
!                  cape calculation associated with the donner 
!                  convection parameterization
!     Don_rad      donner_rad_type derived type variable used to hold 
!                  those fields needed to connect the donner deep 
!                  convection parameterization and the model radiation 
!                  package
!
!   intent(out) variables:
!
!     ermesg        character string containing any error message
!                   that is returned from a kernel subroutine
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg= ' '

!---------------------------------------------------------------------
!    allocate the components of the donner_conv_type variable Don_conv.
!    definitions of these arrays are found in donner_types.h.
!---------------------------------------------------------------------
      allocate (Don_conv%cecon              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%ceefc              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cell_liquid_eff_diam     &
                                            (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cell_ice_geneff_diam     &
                                            (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cememf_mod         (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cemfc              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cmus               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%conv_temp_forcing  (isize, jsize, nlev_lsm) )
      allocate (Don_conv%conv_moist_forcing (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cual               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cuqi               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cuql               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%detmfl             (isize, jsize, nlev_lsm) )
      allocate (Don_conv%dgeice             (isize, jsize, nlev_lsm) )
      allocate (Don_conv%dmeml              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%ecds               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%eces               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%elt                (isize, jsize, nlev_lsm) )
      allocate (Don_conv%emds               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%emes               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%fre                (isize, jsize, nlev_lsm) )
      allocate (Don_conv%mrmes              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%tmes               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%uceml              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%umeml              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%wmms               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%wmps               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%xice               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%xliq               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%a1                 (isize, jsize) )
      allocate (Don_conv%amax               (isize, jsize) )
      allocate (Don_conv%amos               (isize, jsize) )
      allocate (Don_conv%ampta1             (isize, jsize) )
      allocate (Don_conv%cell_precip        (isize, jsize) )
      allocate (Don_conv%dcape              (isize, jsize) )
      allocate (Don_conv%emdi_v             (isize, jsize) )
      allocate (Don_conv%meso_precip        (isize, jsize) )
      allocate (Don_conv%pb_v               (isize, jsize) )
      allocate (Don_conv%pmd_v              (isize, jsize) )
      allocate (Don_conv%przm               (isize, jsize) )
      allocate (Don_conv%prztm              (isize, jsize) )
      allocate (Don_conv%pzm_v              (isize, jsize) )
      allocate (Don_conv%pztm_v             (isize, jsize) )

      allocate (Don_conv%qtceme        (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%qtmes1        (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%qtren1        (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%wtp1          (isize, jsize, nlev_lsm, ntr) )

!---------------------------------------------------------------------
!    initialize the components of the donner_conv_type variable Don_conv.
!---------------------------------------------------------------------
      Don_conv%cecon                = 0.0
      Don_conv%ceefc                = 0.0
      Don_conv%cell_liquid_eff_diam = 0.0
      Don_conv%cell_ice_geneff_diam = 0.0
      Don_conv%cememf_mod           = 0.0
      Don_conv%cemfc                = 0.0
      Don_conv%cmus                 = 0.0
      Don_conv%conv_temp_forcing    = 0.0
      Don_conv%conv_moist_forcing   = 0.0
      Don_conv%cual                 = 0.0
      Don_conv%cuqi                 = 0.0
      Don_conv%cuql                 = 0.0
      Don_conv%detmfl               = 0.0
      Don_conv%dgeice               = 0.0
      Don_conv%dmeml                = 0.0
      Don_conv%ecds                 = 0.0
      Don_conv%eces                 = 0.0
      Don_conv%elt                  = 0.0
      Don_conv%emds                 = 0.0
      Don_conv%emes                 = 0.0
      Don_conv%fre                  = 0.0
      Don_conv%mrmes                = 0.0
      Don_conv%tmes                 = 0.0
      Don_conv%uceml                = 0.0
      Don_conv%umeml                = 0.0
      Don_conv%wmms                 = 0.0
      Don_conv%wmps                 = 0.0
      Don_conv%xice                 = 0.0
      Don_conv%xliq                 = 0.0
      Don_conv%a1                   = 0.0
      Don_conv%amax                 = 0.0
      Don_conv%amos                 = 0.0
      Don_conv%ampta1               = 0.0
      Don_conv%cell_precip          = 0.0
      Don_conv%dcape                = 0.0
      Don_conv%emdi_v               = 0.0
      Don_conv%meso_precip          = 0.0
      Don_conv%pb_v                 = 0.0
      Don_conv%pmd_v                = 0.0
      Don_conv%przm                 = 0.0
      Don_conv%prztm                = 0.0
      Don_conv%pzm_v                = 0.0
      Don_conv%pztm_v               = 0.0
      Don_conv%qtceme               = 0.0
      Don_conv%qtmes1               = 0.0
      Don_conv%qtren1               = 0.0
      Don_conv%wtp1                 = 0.0

!---------------------------------------------------------------------
!    allocate the components of the donner_cape_type variable Don_cape.
!    definitions of these arrays are found in donner_types.h.
!---------------------------------------------------------------------
      allocate (Don_cape%coin       (isize, jsize) )
      allocate (Don_cape%plcl       (isize, jsize) )
      allocate (Don_cape%plfc       (isize, jsize) )
      allocate (Don_cape%plzb       (isize, jsize) )
      allocate (Don_cape%qint_lag   (isize, jsize) )
      allocate (Don_cape%qint       (isize, jsize) )
      allocate (Don_cape%xcape_lag  (isize, jsize) )
      allocate (Don_cape%xcape      (isize, jsize) )
      allocate (Don_cape%cape_p     (isize, jsize, nlev_hires) )
      allocate (Don_cape%env_r      (isize, jsize, nlev_hires) )
      allocate (Don_cape%env_t      (isize, jsize, nlev_hires) )
      allocate (Don_cape%parcel_r   (isize, jsize, nlev_hires) )
      allocate (Don_cape%parcel_t   (isize, jsize, nlev_hires) )
      allocate (Don_cape%model_p    (isize, jsize, nlev_lsm) )
      allocate (Don_cape%model_r    (isize, jsize, nlev_lsm) )
      allocate (Don_cape%model_t    (isize, jsize, nlev_lsm) )

!---------------------------------------------------------------------
!    initialize the components of the donner_cape_type variable Don_cape.
!---------------------------------------------------------------------
      Don_cape%coin        = 0.0
      Don_cape%plcl        = 0.0
      Don_cape%plfc        = 0.0
      Don_cape%plzb        = 0.0   
      Don_cape%qint_lag    = 0.0
      Don_cape%qint        = 0.0
      Don_cape%xcape_lag   = 0.0 
      Don_cape%xcape       = 0.0 
      Don_cape%cape_p      = 0.0
      Don_cape%env_r       = 0.0
      Don_cape%env_t       = 0.0
      Don_cape%parcel_r    = 0.0
      Don_cape%parcel_t    = 0.0
      Don_cape%model_p     = 0.0
      Don_cape%model_r     = 0.0
      Don_cape%model_t     = 0.0

!---------------------------------------------------------------------
!    allocate the components of the donner_rad_type variable Don_rad.
!    definitions of these arrays are found in donner_types.h.
!---------------------------------------------------------------------
      allocate (Don_rad%cell_cloud_frac  (isize, jsize, nlev_lsm) )
      allocate (Don_rad%cell_ice_amt     (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%cell_ice_size    (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%cell_liquid_amt  (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%cell_liquid_size (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_cloud_frac  (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_ice_amt     (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_ice_size    (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_liquid_amt  (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_liquid_size (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%nsum             (isize, jsize) )

!---------------------------------------------------------------------
!    initialize the components of the donner_rad_type variable Don_rad 
!    using the input variables supplied.
!---------------------------------------------------------------------
      Don_rad%cell_cloud_frac  = cell_cld_frac
      Don_rad%cell_ice_amt     = cell_ice_amt
      Don_rad%cell_ice_size    = cell_ice_size
      Don_rad%cell_liquid_amt  = cell_liq_amt
      Don_rad%cell_liquid_size = cell_liq_size
      Don_rad%meso_cloud_frac  = meso_cld_frac
      Don_rad%meso_ice_amt     = meso_ice_amt
      Don_rad%meso_ice_size    = meso_ice_size
      Don_rad%meso_liquid_amt  = meso_liq_amt
      Don_rad%meso_liquid_size = meso_liq_size
      Don_rad%nsum             = nsum

!----------------------------------------------------------------------


end subroutine don_d_init_loc_vars_k



!####################################################################

subroutine don_d_column_input_fields_k  &
         (isize, jsize, nlev_lsm, dt, calc_conv_on_this_step, Col_diag, &
          temp, mixing_ratio, pfull, omega, phalf, parcel_rise, ermesg)

use messy_convect_donner_types_mod, only : donner_column_diag_type     

implicit none

!---------------------------------------------------------------------
!    subroutine don_d_column_input_fields_k outputs the 
!    basic profile information for any diagnostic columns.
!---------------------------------------------------------------------

integer,                            intent(in)  :: isize, jsize, nlev_lsm
real,                               intent(in)  :: dt
logical,                            intent(in)  :: calc_conv_on_this_step
type(donner_column_diag_type),      intent(in)  :: Col_Diag
real, dimension(isize,jsize,nlev_lsm),                            &
                                    intent(in)  :: temp, mixing_ratio, &
                                                   pfull, omega
real, dimension(isize,jsize,nlev_lsm+1),                            &
                                    intent(in)  :: phalf
real, dimension(isize,jsize),       intent(in)  :: parcel_rise              
character(len=*),                   intent(out) :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize         size of x-dimension of physics window
!     jsize         size of y-dimension of physics window
!     nlev_lsm      number of layers in large-scale model
!     dt            physics time step [ sec ]
!     calc_conv_on_this_step
!                   logical indicating whether the deep convection
!                   calculation is to be done on this timestep
!     Col_diag      donner_column_diagtype variable containing the
!                   information defining the columns fro which diagnos-
!                   tics are desired.
!     temp          temperature field at model levels [ deg K ]
!     mixing_ratio  vapor mixing ratio field at model levels 
!                   [ kg(h20) / kg(dry air) ]
!     pfull         pressure field at full-levels 1:nlev_lsm    [ Pa ]
!     omega         model omega field at model full levels 
!                   [ Pa / sec ]
!     phalf         pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     parcel_rise   accumulated vertical displacement of a near-surface 
!                   parcel as a result of the lowest model level omega 
!                   field [ Pa ]
!
!   intent(out) variables:
!
!     ermesg        character string containing any error message
!                   that is returned from a kernel subroutine
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer :: idiag, jdiag, unitdiag
      integer :: n, k       

!----------------------------------------------------------------------
!   local variables:
!
!     idiag         physics window i index of current diagnostic column
!     jdiag         physics window j index of current diagnostic column
!     unitdiag      i/o unit assigned to current diagnostic column
!     n, k          do-loop indices
!
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!    loop over the diagnostic columns in this physics window. output
!    the physics timestep and the window coordinates, and whether the 
!    convection calculation is to be done on this timestep.
!---------------------------------------------------------------------
      do n=1,Col_diag%ncols_in_window
        idiag = Col_diag%i_dc(n)
        jdiag = Col_diag%j_dc(n)
        unitdiag = Col_diag%unit_dc(n)
        write (unitdiag, '(a,f8.1, 2i4)')  &
                ' physics timestep, window i, window j= ',  &
                  dt, idiag, jdiag                         
        write (unitdiag,'(a,l4 )' )  ' conv_calc_on_this_step = ',   &
                  calc_conv_on_this_step

!----------------------------------------------------------------------
!    if the calculation is to be done on this timestep, and convection
!    in the column is not precluded by downward motion at the lowest 
!    level, output the column temperature and mixing ratio profiles
!    over the levels at which output has been requested.
!----------------------------------------------------------------------
        if (calc_conv_on_this_step) then
          if (omega(idiag,jdiag,nlev_lsm) < 0.0) then
            write (unitdiag, '(a)')  & 
                       '                      input profiles'
            write (unitdiag, '(a)')  &
                       '   k   press      temp           mixing ratio '
            write (unitdiag, '(a)')  &
                  '        hPa       deg K    g(h2o) / kg (dry air)  '
            do k=Col_diag%kstart,nlev_lsm
              write (unitdiag, '(i4, 2f10.4, 7x, 1pe13.5)')  &
                   k, 1.0E-02*pfull(idiag,jdiag,k), temp(idiag,jdiag,k),&
                   1.0e03*mixing_ratio(idiag,jdiag,k)
            end do
          endif 
        endif 

!---------------------------------------------------------------------
!    output the surface pressure, omega at the lowest level, and the 
!    accumulated parcel displacement.
!---------------------------------------------------------------------
        write (unitdiag,'(a,f13.4,1pe13.5)')  &
                  ' sfcprs (hPa),  omega_btm (Pa/sec)= ',   &
                  1.0E-02*phalf(idiag,jdiag,nlev_lsm+1),   &
                  omega(idiag,jdiag,nlev_lsm) 
        write (unitdiag,'(a,f13.6)')  ' omint (hPa)= ',   &
                  1.0E-02*parcel_rise(idiag,jdiag)
      end do

!---------------------------------------------------------------------


end subroutine don_d_column_input_fields_k 



!####################################################################

subroutine don_d_convection_driver_k    &
         (isize, jsize, nlev_lsm, nlev_hires, ntr, me,  &
          cloud_tracers_present, dt, Nml,    &
          Initialized, Param, Col_diag, temp, mixing_ratio, pfull, & 
          qlin, qiin, qain, lag_cape_temp, lag_cape_vapor,  &
          lag_cape_press, phalf, current_displ, land, sfc_sh_flux,  &
          sfc_vapor_flux, tr_flux, tracers, Don_cape, Don_conv, &
          Don_rad, temperature_forcing, moisture_forcing, total_precip, &
          donner_humidity_ratio, donner_humidity_area, dql, dqi, dqa, &
          exit_flag, ermesg)

use messy_convect_donner_types_mod, only : donner_initialized_type, donner_rad_type, &
                             donner_param_type, donner_conv_type, &
                             donner_nml_type, &
                             donner_column_diag_type, donner_cape_type
USE MESSY_CONVECT_DONNER_RAD_K,     ONLY: don_r_donner_rad_driver_k
USE MESSY_CONVECT_DONNER_CAPE_K,    ONLY: don_c_def_conv_env_k
USE MESSY_CONVECT_DONNER_MESO_K,    ONLY: don_m_define_anvil_ice_k
USE MESSY_CONVECT_DONNER_LSCLOUD_K, ONLY: don_l_lscloud_driver_k
implicit none

!---------------------------------------------------------------------
!    subroutine don_d_convection_driver_k manages the cal-
!    culation of the effects of deep convection on atmospheric fields 
!    by calling routines to lift a parcel, determine if deep convection 
!    results and, if so, obtain the temperature and moisture forcing and 
!    precipitation produced, and the fields needed to assess the effects
!    of the deep convection on the radiative fluxes and heating and 
!    the large-scale cloud fields of the model.
!---------------------------------------------------------------------

integer,                         intent(in) :: isize, jsize, nlev_lsm,  &
                                               nlev_hires,    ntr, me
logical,                         intent(in) :: cloud_tracers_present
real,                            intent(in) :: dt 
type(donner_nml_type),           intent(in) :: Nml
type(donner_initialized_type),   intent(in) :: Initialized
type(donner_param_type),         intent(in) :: Param
type(donner_column_diag_type),   intent(in) :: Col_diag
real,    dimension(isize,jsize,nlev_lsm),              &
                              intent(in)    :: temp, mixing_ratio,  &
                                               pfull, qlin, &
                                               qiin, qain,   &
                                               lag_cape_temp, &
                                               lag_cape_vapor, &
                                               lag_cape_press
real,    dimension(isize,jsize,nlev_lsm+1),                       &
                              intent(in)    ::  phalf                  
real,    dimension(isize,jsize),                                     &
                              intent(in)    :: current_displ, land, &
                                               sfc_sh_flux,   &
                                               sfc_vapor_flux
real,    dimension(isize,jsize,ntr),                             &
                              intent(in)    :: tr_flux        
real,    dimension(isize,jsize,nlev_lsm,ntr),                      &
                              intent(in)    :: tracers        
type(donner_cape_type),       intent(inout) :: Don_cape
type(donner_conv_type),       intent(inout) :: Don_conv
type(donner_rad_type),        intent(inout) :: Don_rad
real,    dimension(isize,jsize,nlev_lsm),                           &
                              intent(out)   :: temperature_forcing,&
                                               moisture_forcing
real,    dimension(isize,jsize),                                  &
                              intent(out)   :: total_precip
real,    dimension(isize,jsize,nlev_lsm),                              &
                              intent(out)   :: donner_humidity_ratio, &
                                               donner_humidity_area, &
                                               dql, dqi, dqa
character(len=*),             intent(out)   :: ermesg
logical, dimension(isize,jsize),                                &
                              intent(out)   :: exit_flag

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     dt             physics time step [ sec ]
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     Initialized    donner_initialized_type variable containing
!                    variables which are defiuned during initialization.
!                    these values may be changed during model execution.
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!                    tion parameterization
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   vapor mixing ratio field at model levels 
!                    [ kg(h20) / kg(dry air) ]
!     pfull          pressure field on model full levels [ Pa ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qain           large-scale cloud fraction  
!                    [ fraction ]
!     lag_cape_temp  temperature field used in lag-time cape 
!                    calculation [ deg K ]
!     lag_cape_vapor vapor mixing ratio field used in lag-time
!                    cape calculation [ kg(h2o) / kg(dry air) ]
!     lag_cape_press model full-level pressure field used in 
!                    lag-time cape calculation  [ Pa ]
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     current_displ  low-level parcel displacement to use in cape
!                    calculation on this step [ Pa ]
!     land           fraction of grid box covered by land
!                    [ fraction ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     tr_flux        surface flux of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer mixing ratios of tracers transported by the
!                    donner deep convection parameterization
!                    [ kg(tracer) / kg (dry air) ]
!
!   intent(inout) variables:
!
!     Don_cape       donner_cape type derived type variable containing 
!                    diagnostics and intermediate results related to the
!                    cape calculation associated with the donner convec-
!     Don_conv       donner_convection_type derived type variable
!                    containing diagnostics and intermediate results 
!                    describing the nature of the convection produced by
!                    the donner parameterization
!     Don_rad        donner_rad_type derived type variable used to hold 
!                    those fields needed to connect the donner deep 
!                    convection parameterization and the model radiation 
!                    package
!
!   intent(out) variables:
!
!     temperature_forcing  
!                    temperature tendency due to donner convection
!                    [ deg K / sec ]
!     moisture_forcing  
!                    vapor mixing ratio tendency due to donner 
!                    convection [ kg(h2o) / (kg(dry air) sec ) ]
!     total_precip   total precipitation rate produced by the
!                    donner parameterization [ mm / day ]
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to specific 
!                    humidity in environment outside convective system
!                    [ dimensionless ]
!     donner_humidity_area
!                    fraction of grid box in which humidity is affected
!                    by the deep convection, defined as 0.0 below cloud
!                    base and above the mesoscale updraft, and as the
!                    sum of the cell and mesoscale cloud areas in 
!                    between. it is used in strat_cloud_mod to determine
!                    the large-scale specific humidity field for the
!                    grid box. DO NOT use for radiation calculation,
!                    since not all of this area includes condensate.
!                    [ fraction ]
!     dql            tendency of cloud liquid specific humidity
!                    due to donner convection 
!                    [ kg(h2o) / kg(moist air) / sec ]
!     dqi            tendency of cloud ice specific humidity
!                    due to donner convection 
!                    [ kg(h2o) / kg(moist air) / sec ]
!     dqa            tendency of large-scale cloud area
!                    due to donner convection 
!                    [ fraction / sec ]
!     exit_flag      logical array indicating whether deep convection 
!                    exists in a column
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!    call don_c_def_conv_env_k to determine how 
!    a parcel will behave with respect to deep convection in each model
!    column.
!---------------------------------------------------------------------
      call don_c_def_conv_env_k          &
           (isize, jsize, nlev_lsm, nlev_hires, Nml, Param, Col_diag, &
            temp, mixing_ratio, pfull, lag_cape_temp, lag_cape_vapor, &
            lag_cape_press, current_displ, Don_cape, Don_conv, ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    call don_d_cupar to calculate the normalized deep convective 
!    forcing.
!---------------------------------------------------------------------
      call don_d_cupar_k     &
           (isize, jsize, nlev_lsm, nlev_hires, ntr, me, dt, Col_diag, &
            Param, Nml, Initialized, current_displ, sfc_sh_flux,    &
            sfc_vapor_flux, pfull, temp, phalf, tr_flux, tracers,    &
            Don_conv, Don_cape, temperature_forcing, moisture_forcing, &
            total_precip, ermesg, exit_flag)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    call define_donner_anvil_ice to define the ice content profile
!    (Don_conv%xice) and the pressures at top and bottom of mesoscale
!    anvil (Don_conv%prztm, Don_conv%przm).
!----------------------------------------------------------------------
      call don_m_define_anvil_ice_k   &
           (isize, jsize, nlev_lsm, Param, Col_diag, pfull, temp,    &
            exit_flag, Don_conv, ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    call donner_rad_driver to define the cloud ice, cloud liquid and
!    cloud areas of the cell and mesoscale clouds associated with 
!    donner convection so as to make them available to the radiation
!    code.
!---------------------------------------------------------------------
      call don_r_donner_rad_driver_k   &
           (isize, jsize, nlev_lsm, Param, Col_diag, Initialized, &
            pfull, temp, land, exit_flag, Don_conv, Don_rad, Nml, ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    call donner_lscloud_driver to provide the connection between 
!    donner convection and the large-scale cloud scheme. 
!---------------------------------------------------------------------
      call don_l_lscloud_driver_k   &
           (isize, jsize, nlev_lsm, cloud_tracers_present, Param,  &
            Col_diag, pfull, temp,   &
            mixing_ratio, qlin, qiin, qain, phalf, Don_conv, &
            donner_humidity_ratio, donner_humidity_area, dql, dqi,  &
            dqa, ermesg) 

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!!!  QUESTION 3:
!!     A BETTER COMMENT IS NEEDED HERE --
!!     WHY IS THIS REMOVAL NECESSARY ??
!    save the vapor and temperature forcing resulting from the donner
!    deep parameterization. remove from them the contributions due to  
!    vertical transport by the donner mass flux and the mesoscale flux.
!    also remove vapor and temperature tendencies
!    corresponding to these increments from the Donner cumulus
!    thermal forcing and moisture forcing, which included
!    them as evaporatation and/or sublimation in mulsub.
!    assumptions used in strat_cloud_donner_tend to relate detrainment
!    to net mass fluxes differ from those in mulsub, so the
!    increments here do not balance those in mulsub. the difference
!    remains as a phase change.
!    mulsub allowed ice and liquid from convective system to evaporate
!    and/or sublimate as part of thermal and moisture forcing terms
!    remove those tendencies here. different assumptions used to
!    calculate these increments/tendencies here and in mulsub, so
!    some residual phase change will generally remain.
!---------------------------------------------------------------------
      Don_conv%conv_temp_forcing(:,:,:)  = temperature_forcing(:,:,:)
      Don_conv%conv_moist_forcing(:,:,:) = moisture_forcing(:,:,:)

      if (cloud_tracers_present) then
        moisture_forcing(:,:,:) = moisture_forcing(:,:,:) - &
                                            dql(:,:,:) - dqi(:,:,:)
        temperature_forcing(:,:,:) = temperature_forcing(:,:,:) + &
                               dql(:,:,:)*Param%HLV   /(Param%cp_air) + &
                               dqi(:,:,:)*(Param%HLS  )/(Param%cp_air)
      endif

!---------------------------------------------------------------------


end subroutine don_d_convection_driver_k

!######################################################################

subroutine don_d_cupar_k     &
         (isize, jsize, nlev_lsm, nlev_hires, ntr, me, dt, Col_diag, &
          Param, Nml, Initialized, current_displ, sfc_sh_flux,   &
          sfc_vapor_flux, pfull, temp, phalf, tr_flux, tracers,   &
          Don_conv, Don_cape, temperature_forcing, moisture_forcing,  &
          total_precip, ermesg, exit_flag)

!----------------------------------------------------------------------
!    subroutine cupar drives the parameterization for deep cumulus 
!    convection. it returns the temperature and moisture forcing assoc-
!    iated with deep convection, the total convective precipitation
!    and various diagnostics contained in Don_conv and Don_cape to the 
!    calling routine. it is based on (Donner, 1993, J.Atmos.Sci.).
!---------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_initialized_type, donner_nml_type, &
                             donner_param_type, donner_conv_type, &
                             donner_column_diag_type, donner_cape_type
implicit none

!--------------------------------------------------------------------- 
integer,                           intent(in)    :: isize, jsize,    &
                                                    nlev_lsm,    &
                                                    nlev_hires,&
                                                    ntr, me
real,                              intent(in)    :: dt
type(donner_column_diag_type),     intent(in)    :: Col_diag
type(donner_param_type),           intent(in)    :: Param
type(donner_nml_type),             intent(in)    :: Nml
type(donner_initialized_type),     intent(in)    :: Initialized
real,    dimension(isize,jsize),   intent(in)    :: current_displ, &
                                                    sfc_sh_flux,  &
                                                    sfc_vapor_flux
real,    dimension(isize,jsize,nlev_lsm),                      &
                                   intent(in)    :: pfull, temp
real,    dimension(isize,jsize,nlev_lsm+1),                    &  
                                   intent(in)    :: phalf
real,    dimension(isize,jsize,ntr),                           &
                                   intent(in)    :: tr_flux
real,    dimension(isize,jsize,nlev_lsm,ntr),               &
                                   intent(in)    :: tracers
type(donner_conv_type),            intent(inout) :: Don_conv
type(donner_cape_type),            intent(inout) :: Don_cape
real,    dimension(isize,jsize,nlev_lsm),                      &
                                   intent(out)   :: temperature_forcing,&
                                                    moisture_forcing
real,    dimension(isize,jsize),   intent(out)   :: total_precip
character(len=*),                  intent(out)   :: ermesg
logical, dimension(isize,jsize),   intent(out)   :: exit_flag

!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     dt             physics time step [ sec ]
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!                    tion parameterization
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     Initialized    donner_initialized_type variable containing
!                    variables which are defiuned during initialization.
!                    these values may be changed during model execution.
!     current_displ  low-level parcel displacement to use in cape
!                    calculation on this step [ Pa ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     pfull          pressure field at model full levels [ Pa ]
!     temp           temperature field at model full levels [ deg K ]
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     tr_flux        flux across the surface of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer fields that are to be transported by donner
!                    convection [ kg (tracer) / kg (dry air) ]
!
!   intent(out) variables:
!    
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!   intent(inout) variables:
!
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!     Don_cape       donner_cape_type derived type variable containing
!                    fields associated with the calculation of
!                    convective available potential energy (cape).
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:


      real, dimension (isize, jsize, nlev_lsm, ntr) :: xgcm_v
      real, dimension (isize, jsize, ntr)           :: sfc_tracer_flux
      integer                                       :: i, j, k, n 

!---------------------------------------------------------------------
!   local variables:
!
!      xgcm_v              tracer mixing ratio fields transported by 
!                          donner convection, index 1 nearest surface
!                          [ kg(tracer) / kg (dry air) ]
!      sfc_tracer_flux     tracer flux across the surface
!                          [ kg(tracer) / (m**2 sec) ]
!      i, j, k, n          do-loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!    write a message to the output file for each diagnostic column in 
!    this window.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          write (Col_diag%unit_dc(n), '(a, 2i4)')  &
           'cupar: entering cupar with i_dc, j_dc:',    &
                           Col_diag%i_dc(n), Col_diag%j_dc(n)
        end do
      endif

!----------------------------------------------------------------------
!    call donner_deep_check_for_deep_convection_k to determine if deep 
!    convection may at this time be precluded in any of the columns of 
!    this physics window. logical array exit_flag is returned, with a 
!    value of .false. if donner convection is still allowed, a value of 
!    .true. if deep convection is precluded in a particular coluumn.
!----------------------------------------------------------------------
      call don_d_check_for_deep_conv_k   &
           (isize, jsize, nlev_lsm, dt, Param, Col_diag, &
            current_displ,Don_cape, Don_conv, exit_flag, ermesg)
!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!-------------------------------------------------------------------
!    for the tracers that are to be transported by donner_deep_mod,
!    define the tracer input fields that will be needed by the donner
!    cloud model.
!-------------------------------------------------------------------
      if (Initialized%do_donner_tracer) then

!-------------------------------------------------------------------
!    define the tracer fluxes across the surface.
!-------------------------------------------------------------------
        sfc_tracer_flux(:,:,:) = tr_flux(:,:,:)

!-------------------------------------------------------------------
!    define an inverted tracer profile (index 1 nearest ground) for use
!    in the cloud and convection routines.
!------------------------------------------------------------------
        do k=1,nlev_lsm
          xgcm_v(:,:,k,:) = tracers(:,:,nlev_lsm-k+1,:)
        end do

!--------------------------------------------------------------------
!    if tracers are not to be transported by donner_deep_mod, define
!    these tracer input fields to be 0.0.
!--------------------------------------------------------------------
      else
        xgcm_v = 0.
        sfc_tracer_flux = 0.0
      endif 
      
!---------------------------------------------------------------------
!    call subroutine mulsub to calculate normalized (in-cloud) cumulus 
!    forcings, one column at a time. the forcings are normalized by the 
!    cloud area at cloud base level a_1(p_b).
!---------------------------------------------------------------------
      call don_d_mulsub_k   &
           (isize, jsize, nlev_lsm, nlev_hires, ntr, me, dt, Param,   &
!++lwh
            Nml, Col_diag, Initialized, phalf, sfc_vapor_flux, sfc_sh_flux, &
!--lwh
            sfc_tracer_flux, xgcm_v, Don_cape, Don_conv, exit_flag,  &
            total_precip, temperature_forcing, moisture_forcing, &
            ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    call remove_normalization to remove the normalization from the 
!    deep convection diagnostics and forcing terms by multiplying them 
!    by the fractional cloud area. the values thus become grid-box 
!    averages, rather than averages over the cloudy area, and so are 
!    ready to use in the large-scale model equations. all columns in
!    which exit_flag is .true. are given zero values for total_precip,
!    temperature_forcing and moisture_forcing.
!---------------------------------------------------------------------
      call don_d_remove_normalization_k   &
           (isize, jsize, nlev_lsm, exit_flag, Don_conv, total_precip, &
            temperature_forcing, moisture_forcing, ermesg)


!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    in a diagnostics window, output various desired quantities.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          call don_d_output_cupar_diags_k    &
               (isize, jsize, nlev_lsm, Col_diag, n, exit_flag,   &
                total_precip, temperature_forcing, Don_conv, Don_cape, &
                ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (trim(ermesg) /= ' ') return
        end do
      endif  ! (in_diagnostics_window)

!---------------------------------------------------------------------


end subroutine don_d_cupar_k


!#####################################################################


subroutine don_d_check_for_deep_conv_k   &
           (isize, jsize, nlev_lsm, dt, Param, Col_diag, &
            current_displ, Don_cape, Don_conv, exit_flag, ermesg)

!---------------------------------------------------------------------
!    subroutine don_d_check_for_deep_conv_k tests for the 
!    sounding- and upward-motion-based criteria which will prevent deep 
!    convection from occurring in a column. if convection is precluded, 
!    the logical variable exit_flag is set to .true. and additional cal-
!    culations in that column will be skipped. if convection is deter-
!    mined to be possible, exit_flag is set to .false., and additional 
!    calculations in the column will be done.
!---------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_param_type, donner_conv_type, &
                             donner_column_diag_type, donner_cape_type

implicit none

!---------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize, nlev_lsm
real,                             intent(in)    :: dt
type(donner_param_type),          intent(in)    :: Param
type(donner_column_diag_type),    intent(in)    :: Col_diag
real,    dimension(isize,jsize),  intent(in)    :: current_displ
type(donner_cape_type),           intent(inout) :: Don_cape
type(donner_conv_type),           intent(inout) :: Don_conv
logical, dimension (isize,jsize), intent(out)   :: exit_flag 
character(len=*),                 intent(out)   :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     dt             physics time step [ sec ]
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!                    tion parameterization
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     current_displ  low-level parcel displacement to use in cape
!                    calculation on this step [ Pa ]
!
!   intent(inout) variables:
!
!     Don_cape       donner_cape_type derived type variable containing
!                    fields associated with the calculation of
!                    convective available potential energy (cape).
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!
!   intent(out) variables:
!
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
 
      real, dimension (isize, jsize)  ::  pdeet1, pdeet2
      integer                         :: idiag, jdiag, unitdiag
      integer                         :: i, j, k, n 

!---------------------------------------------------------------------
!   local variables:
!
!     pdeet1              pressure depth between the level of free 
!                         convection and the level of zero buoyancy 
!                         [ Pa ]
!     pdeet2              pressure depth between the level of free 
!                         convection and the pressure at lowest 
!                         large-scale model grid level 
!                         [ Pa ]
!     idiag               physics window i index of current diagnostic
!                         column
!     jdiag               physics window j index of current diagnostic
!                         column
!     unitdiag            i/o unit assigned to current diagnostic
!                         column
!     i, j, k, n          do-loop indices
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!    process each column in the physics window.
!---------------------------------------------------------------------
      do j=1,jsize     
        do i=1,isize     

!---------------------------------------------------------------------
!    define the time rates of change of column convective available 
!    potential energy (Don_conv%dcape_v). 
!---------------------------------------------------------------------
          Don_conv%dcape(i,j) = (Don_cape%xcape(i,j) -          &
                                 Don_cape%xcape_lag(i,j))/dt

!---------------------------------------------------------------------
!    define the pressure depth between the level of free convection
!    and the level of zero buoyancy (pdeet1) and the pressure depth 
!    between the level of free convection and the pressure at lowest 
!    large-scale model grid level (pdeet2).
!---------------------------------------------------------------------
          pdeet1(i,j) = Don_cape%plfc(i,j) - Don_cape%plzb(i,j)
          pdeet2(i,j) = Don_cape%plfc(i,j) - Don_cape%model_p(i,j,1)

!---------------------------------------------------------------------
!    check that all criteria for deep convection are satisfied; if so,
!    set exit_flag to be .false., if any of the criteria are not sat-
!    isfied, set exit_flag to .true. the criteria which can be evaluated
!    at this time are:
!       1) cape (Don_Cape%xcape) must be positive;
!       2) cape must be increasing with time (Don_conv%dcape > 0);
!       3) pressure depth between lfc and lzb (pdeet1) must be greater 
!          than pdeep_cv;
!       4) the time-integrated upward displacement of a parcel from the
!          lowest model level (current_displ) must be sufficient to 
!          allow the parcel to have reached the lfc;
!       5) convective inhibition must be less than cdeep_cv.
!---------------------------------------------------------------------
          if ((Don_cape%xcape(i,j) <= 0.)        .or.  &
              (Don_conv%dcape(i,j) <= 0.)        .or. & 
              (pdeet1(i,j) < Param%pdeep_cv )    .or. &
              (pdeet2(i,j) < current_displ(i,j)) .or.  &
              (Don_cape%coin(i,j) > Param%cdeep_cv) )   then
            exit_flag(i,j) = .true.
          else
            exit_flag(i,j) = .false.
          endif
        end do
      end do

!--------------------------------------------------------------------
!    if in diagnostics window, output info concerning the status of
!    deep convection in this column.
!--------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          idiag = Col_diag%i_dc(n)
          jdiag = Col_diag%j_dc(n)
          unitdiag = Col_diag%unit_dc(n)

!--------------------------------------------------------------------
!    for any diagnostic columns in the window for which deep convection
!    is possible, output the integrated upward displacement 
!    (current_displ), the value of cape (Don_Cape%xcape), the time rate
!    of change of cape (Don_conv%dcape) and the logical variable 
!    indicating whether deep convection is precluded in this column 
!    (exit_flag).
!--------------------------------------------------------------------
          if (.not. exit_flag(idiag,jdiag)) then
            write (unitdiag, '(a, f20.12, f20.12, e20.12, l4)')   &
                  'in cupar: omint,cape,dcape, exit_flag',    &
                         current_displ  (idiag,jdiag),   &
                         Don_Cape%xcape (idiag,jdiag),       &
                         Don_conv%dcape (idiag,jdiag),       &
                         exit_flag      (idiag,jdiag)  

!--------------------------------------------------------------------
!    output various thermodynamic parameters (cp_air, cp_vapor, d622,  &
!    rdgas, hlv, rvgas), various sounding levels and features (cape, 
!    convective inhibition, level of zero buoyancy, level of free 
!    convection and the model soundings (p, t, mixing ratio).
!--------------------------------------------------------------------
            write (unitdiag, '(a, 2f12.4)')   &
                   'in cupar: cpi,cpv= ',Param%cp_air, Param%cp_vapor
            write (unitdiag, '(a, 2f12.6, f12.2)')  &
                   'in cupar: rocp,rair,latvap= ',Param%d622, &
                                       Param%rdgas, Param%hlv   
            write (unitdiag, '(a, f12.7)') 'in cupar: rvap= ',Param%rvgas
            write (unitdiag, '(a, 2f14.7, f19.10)')  &
                    'in cupar: cape,cin,plzb= ',  &
                  Don_cape%xcape(idiag,jdiag), &
                  Don_cape%coin(idiag,jdiag), &
                  Don_cape%plzb(idiag,jdiag)
            write (unitdiag, '(a, f19.10)') 'in cupar: plfc= ', &
                  Don_cape%plfc(idiag,jdiag)
            do k=1,nlev_lsm-Col_diag%kstart+1
              write (unitdiag, '(a, i4, f19.10, f20.14, e20.12)') &
                                   'in cupar: k,pr,t,q= ',k,   &
                    Don_cape%model_p(idiag,jdiag,k),   &
                    Don_cape%model_t(idiag,jdiag,k),   &
                    Don_cape%model_r(idiag,jdiag,k)
            end do

!----------------------------------------------------------------------
!    if convection is precluded, output information indicating why.
!----------------------------------------------------------------------
          else
            write (unitdiag, '(a)')   &
               'in cupar: exit_flag is .true., no further calculations&
                              & in this column at this time'
            write (unitdiag, '(a)')   &
                'in cupar: reason(s) for no deep convection:'

!----------------------------------------------------------------------
!    case of no upward motion at lowest level:    
!----------------------------------------------------------------------
            if (current_displ(idiag,jdiag) == 0) then
              write (unitdiag, '(a)')   &
                    'no upward motion at lowest level'
            else 

!----------------------------------------------------------------------
!    case of non-positive cape:    
!----------------------------------------------------------------------
              if (Don_cape%xcape(idiag,jdiag) <= 0.) then      
                write (unitdiag, '(a, f20.12)')   &
                       'non-positive cape, cape = ', &
                           Don_Cape%xcape(idiag,jdiag)      
              endif

!----------------------------------------------------------------------
!    case of non-positive cape time tendency:    
!----------------------------------------------------------------------
              if (Don_conv%dcape(idiag,jdiag) <= 0.) then      
                write (unitdiag, '(a, f20.12)')   &
                      'non-positive cape time tendency, dcape = ', &
                            Don_conv%dcape(idiag,jdiag)      
              endif

              if (Don_cape%plfc(idiag,jdiag) == 0.0 .or.   &
                  Don_cape%plzb(idiag,jdiag) == 0.0) then

!----------------------------------------------------------------------
!    case of sounding not having a level of free convection for 
!    specified parcel:    
!----------------------------------------------------------------------
                if (Don_cape%plfc(idiag,jdiag) == 0.0 ) then
                  write (unitdiag, '(a)')   &
                    'lfc is not definable for parcel used in cape &
                        &calculation'
                endif

!----------------------------------------------------------------------
!    case of sounding not having a level of zero buoyancy for 
!    specified parcel:    
!----------------------------------------------------------------------
                if (Don_cape%plzb(idiag,jdiag) == 0.0) then
                  write (unitdiag, '(a)')   &
                    'lzb is not definable for parcel used in cape &
                        &calculation'
                endif
              else 


!----------------------------------------------------------------------
!    case of sounding not providing a deep enough layer of positive
!    buoyancy:    
!----------------------------------------------------------------------
                if (pdeet1(idiag,jdiag) < Param%pdeep_cv) then      
                  write (unitdiag, '(a, f20.12, a, f20.12,a)')   &
                       'depth of positive buoyancy too shallow, &
                         &plfc - plzb = ',    &
                           pdeet1(idiag,jdiag)*1.0e-02, ' hPa, &
                         & needed depth =', Param%pdeep_cv*1.0e-02, ' hPa'
                endif
                if (Don_cape%plfc(idiag,jdiag) ==  0.0) then 

!----------------------------------------------------------------------
!    case of parcel having insufficient displacement to reach the level
!    of free convection:
!----------------------------------------------------------------------
                else if        &
                  (pdeet2(idiag,jdiag) < current_displ(idiag,jdiag)) then
                  write (unitdiag, '(a, f20.12, a, f20.12, a)')   &
                      'parcel displacement insufficient to reach lfc, &
                       &displacement =',    &
                           current_displ(idiag,jdiag)*1.0e-02, ' hPa, &
                       &needed displacement = ',  &
                            pdeet2(idiag,jdiag)*1.0e-02, ' hPa'
                endif
              endif

!----------------------------------------------------------------------
!    case of sounding having too much convective inhibition:
!----------------------------------------------------------------------
              if (Don_cape%coin(idiag,jdiag) > Param%cdeep_cv) then      
                write (unitdiag, '(a, f20.12, a, f20.12)')   &
                       'convective inhibition too large, cin   = ', &
                             Don_cape%coin(idiag,jdiag), &
                            'max allowable =', Param%cdeep_cv
              endif
            endif
          endif  ! (not exit_flag)
        end do
      endif

!--------------------------------------------------------------------


end subroutine don_d_check_for_deep_conv_k



!#####################################################################

subroutine don_d_mulsub_k   &
         (isize, jsize, nlev_lsm, nlev_hires, ntr, me, dt, Param, Nml, &
!++lwh
          Col_diag, Initialized, phalf, sfc_vapor_flux, sfc_sh_flux, &
!--lwh
          sfc_tracer_flux, xgcm_v, Don_cape, Don_conv, exit_flag, &
          total_precip, temperature_forcing, moisture_forcing, ermesg)

!--------------------------------------------------------------------
!    subroutine don_d_mulsub_k calculates the thermal and moisture
!    forcing produced by an ensemble of cumulus elements and any meso-
!    scale circulation which the ensemble induces, following Donner 
!    (1993, JAS). See also LJD notes, "Cu Closure A," 2/97. calculations
!    at and below this subroutine level are done a column at a time, in 
!    only those columns for which the possibility of deep convection has
!    not yet been ruled out.
!
!                L. Donner  GFDL 27 Apr 97
!---------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_param_type, donner_conv_type, &
                             donner_nml_type, donner_column_diag_type, &
!++lwh
                             donner_cape_type, donner_initialized_type
!--lwh
use messy_convect_donner_meso_k,    only : don_m_meso_effects_k
implicit none

!---------------------------------------------------------------------
integer,                      intent(in)     ::  isize, jsize, nlev_lsm,&
                                                 nlev_hires, ntr, me
real,                         intent(in)     ::  dt
type(donner_param_type),      intent(in)     ::  Param
type(donner_nml_type),        intent(in)     ::  Nml  
type(donner_column_diag_type),                           &
                              intent(in)     ::  Col_diag
!++lwh
type(donner_initialized_type), intent(in)    :: Initialized
!--lwh
real,    dimension(isize,jsize,nlev_lsm+1),                    &
                              intent(in)     ::  phalf
real,    dimension(isize,jsize),                                &
                              intent(in)     ::  sfc_vapor_flux,  &
                                                 sfc_sh_flux
real,    dimension(isize,jsize,ntr),                              &
                              intent(in)     ::  sfc_tracer_flux       
real,    dimension(isize,jsize,nlev_lsm,ntr),               &
                              intent(in)     ::  xgcm_v
type(donner_cape_type),       intent(inout)  ::  Don_cape
type(donner_conv_type),       intent(inout)  ::  Don_conv
logical, dimension(isize,jsize),                            &
                              intent(inout)  ::  exit_flag
real,    dimension(isize,jsize),                          &
                              intent(out)    ::  total_precip        
real,    dimension(isize,jsize,nlev_lsm),                        &
                              intent(out)    ::  temperature_forcing, &
                                                 moisture_forcing
character(len=*),             intent(out)    ::  ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     dt             physics time step [ sec ]
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_tracer_flux 
!                    flux across the surface of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     xgcm_v         tracer fields that are to be transported by donner
!                    convection. index 1 nearest the ground.
!                    [ kg (tracer) / kg (dry air) ]
!
!   intent(inout) variables:
!
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!     Don_cape       donner_cape_type derived type variable containing
!                    fields associated with the calculation of
!                    convective available potential energy (cape).
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!
!   intent(out) variables:
!    
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!---------------------------------------------------------------------

!
!     On Output:
!     
!     ampt             mesoscale cloud fraction, normalized by a(1,p_b)
!     contot           ratio of convective to total precipitation
!     cmui             normalized vertical integral of mesoscale-updraft
!                      deposition (kg(H2O)/((m**2) sec)
!     cmus(nlev)       normalized mesoscale-updraft deposition
!                      (kg(H2O)/kg/sec)
!     cual(nlev)       cloud fraction, cells+meso, normalized by a(1,p_b)
!     cuq(nlev)        ice content in cells, weighted by cell area,
!                      (kg(H2O)/kg)
!                      index 1 at model bottom
!     cuqll(nlev)      liquid content in cells, weighted by cell area,
!                      (kg(H2O)/kg)
!                      index 1 at model bottom
!     ecds(nlev)       normalized convective downdraft evaporation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     eces(nlev)       normalzed convective-updraft evporation/sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emds(nlev)       normalized mesoscale-downdraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emei             normalized vertical integral of mesoscale-updraft
!                      sublimation (kg(h2O)/((m**2) sec)
!     emes(nlev)       normalized mesoscale-updraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     disa(nlev)       normalized thermal forcing, cells+meso (K/sec)
!                      (excludes convergence of surface heat flux)
!                      index 1 at ground. Cumulus thermal forcing defined
!                      as in Fig. 3 of Donner (1993, JAS).
!     disb(nlev)       normalized cell entropy-flux convergence (K/sec)
!                      (excludes convergence of surface flux)
!                      index 1 at ground. Entropy-flux convergence divided
!                      by (p0/p)**(rd/cp).
!     disc(nlev)       normalized cell condensation/deposition
!                      (K/sec)
!                      index 1 at ground.
!     disd(nlev)       normalized cell moisture-flux convergence
!                      (excludes convergence of surface moisture flux)
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dise(nlev)       normalized moisture forcing, cells+meso (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dmeml(nlev)      mass flux in mesoscale downdraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!     elt(nlev)        normalized melting (K/sec)
!                      index 1 at ground.
!     fre(nlev)        normalized freezing (K/sec)
!                      index 1 at ground.
!     pb               pressure at base of cumulus updrafts (Pa)
!     pmd              pressure at top of mesoscale downdraft (Pa)
!     pztm             pressure at top of mesoscale updraft (Pa)
!     mrmes(nlev)       normalized mesoscale moisture-flux convergence
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     qtmes(nlev,ncont)  tracer tendency due to mesoscale tracer-flux
!                        convergence (kg/kg/s) (normalized by a(1,p_b))
!                        index 1 at ground 
!     qtren_v          normalized tracer tendency due to cells...
!                      (lon,lat,vert,tracer index)
!                      Vertical index increases as height increases.
!     sfcq(nlev)       boundary-layer mixing-ratio tendency due to surface
!                      moisture flux (kg(H2O)/kg/sec)
!     sfch(nlev)       boundary-layer heating due to surface heat flux
!                      (K/sec)
!     tmes(nlev)       normalized mesoscale entropy-flux convergence
!                      (K/sec)
!                      Entropy-flux convergence is mesoscale component
!                      of second term in expression for cumulus thermal
!                      forcing in Fig. 3 of Donner (1993, JAS).
!                      index 1 at ground.
!     tpre_v           total normalized precipitation (mm/day)
!     detmfl(nlev)     normalized detrained mass flux from cell
!                      updrafts (kg/((m**2)*s)
!                      (index 1 at atmosphere bottom)
!     uceml(nlev)      normalized mass fluxes in cell updrafts
!                      (kg/((m**2)*s) 
!     umeml(nlev)      mass flux in mesoscale updraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!                      index 1 at ground.
!     wmms(nlev)       normalized mesoscale deposition of water vapor from
!                      cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wmps(nlev)       normalized mesoscale redistribution of water vapor
!                      from cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wtp_v            tracer redistributed by mesoscale processes
!                      (kg/kg/s) (normalized by a(1,p_b))
!                      vertical index increases with increasing height
!                      (lon,lat,vert,tracer index)
!--------------------------------------------------------------------



!!  UNITS
!!     ensmbl_anvil_cond  ! [mm / day ]
!!    ucemh  [kg /sec / m**2 ]
!!    detmfh [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    h1_2   [ deg K / sec ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    enctf  [ deg K / day ]
!!    encmf  [ g(h2o) / kg(air) /day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]
!!    disg   [ deg K / day ]


!        cape     convective available potential energy (J/kg)
!        cin      convective inhibtion (J/kg)
!        cpd      specific heat of dry air at constant pressure (J/(kg K))
!        cpv      specific heat of water vapor [J/(kg K)]
!        dcape    local rate of CAPE change by all processes
!                 other than deep convection [J/(kg s)]
!        dqls     local rate of change in column-integrated vapor
!                 by all processes other than deep convection
!                 {kg(H2O)/[(m**2) s]}
!        epsilo   ratio of molecular weights of water vapor to dry air
!        gravm    gravity constant [m/(s**2)]
!        ilon     longitude index
!        jlat     latitude index
!        mcu      frequency (in time steps) of deep cumulus
!        current_displ  integrated low-level displacement (Pa)
!        cape_p   pressure at Cape.F resolution (Pa)
!                 Index 1 at bottom of model.
!        plfc     pressure at level of free convection (Pa)
!        plzb     pressure at level of zero buoyancy (Pa)
!        pr       pressure at Skyhi vertical resolution (Pa)
!                 Index 1 nearest ground  
!        q        large-scale vapor mixing ratio at Skyhi vertical resolution
!                 [kg(h2O)/kg]
!                 Index 1 nearest ground 
!        qlsd     column-integrated vapor divided by timestep for cumulus
!                 parameterization {kg(H2O)/[(m**2) s]}
!        r        large-scale vapor mixing ratio at Cape.F resolution
!                 [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rpc      parcel vapor mixing ratio from Cape.F [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rd       gas constant for dry air (J/(kg K))
!        rlat     latent heat of vaporization (J/kg)
!        rv       gas constant for water vapor (J/(kg K))
!        t        large-scale temperature at Skyhi vertical resolution (K)
!                 Index 1 nearest ground
!        tcape    large-scale temperature at Cape.F resolution (K)
!                 Index 1 at bottom of model.
!        tpc      parcel temperature from from Cape.F (K)
!                 Index 1 at bottom of model.
!
!     On Input as Parameters:
!
!        kmax     number of vertical levels at Skyhi resolution
!        kpar     number of cumulus sub-ensembles
!        ncap     number of vertical levels in Cape.F resolution
!




!      disa_v              thermal forcing due to deep convection
!                          index 1 nearest surface, normalized by 
!                          cloud area  [ deg K / sec ]
!      dise_v              moisture forcing due to deep convection
!                          index 1 nearest surface, normalized by 
!                          cloud area  [ kg(h2o) / (kg(dry air) *sec ) ]

!----------------------------------------------------------------------
!   local variables:

      real, dimension(nlev_lsm)               ::        &
              ensmbl_cloud_area, cutotal, cmus_tot, cuq, cuql_v, disa, &
              disb, disc, disd, dise, dmeml, uceml, umeml, ecds, eces,  &
              emds, emes, mrmes, tmes, wmms, wmps, detmfl,    &
              meso_cloud_area, disf, disg, disg_2, disn, enctf, encmf, &
              enev, ensmbl_melt, anvil_precip_melt, ensmbl_freeze,&
              temp_tend_freeze, temp_tend_melt, sfcq, sfch

      real, dimension(isize, jsize, nlev_lsm) :: disa_v, dise_v
      real, dimension (nlev_hires)            :: rlsm, emsm, cld_press
      real, dimension( nlev_lsm,ntr)          :: qtmes, qtren, wtp
      real, dimension (nlev_hires,ntr)        :: etsm
      real, dimension (nlev_lsm+1)            :: phalf_c               

      real         ::  al, ampta1, ensmbl_cond, pb, ensmbl_precip,  &
                       pt_ens, &
                       ensmbl_anvil_cond, max_depletion_rate, dqls_v,  &
                       qlsd_v
      logical      ::  lmeso, debug_ijt
      integer      :: diag_unit
      integer      :: kinv
      integer      :: kcont
      integer      :: i, j, k, n


!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' '

!----------------------------------------------------------------------
!    initialize the output arrays.
!----------------------------------------------------------------------
      temperature_forcing = 0.
      moisture_forcing    = 0.
      total_precip        = 0.

!---------------------------------------------------------------------
!    output a message to all diagnostic files indicating entry into
!    subroutine don_d_mulsub_k.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          write(Col_diag%unit_dc(n), '(a, 2i4)')    &
            'in mulsub: i_dc,j_dc= ', Col_diag%i_dc(n), Col_diag%j_dc(n)
        end do
      endif

!--------------------------------------------------------------------
!    LOOP OVER COLUMNS IN CURRENT PHYSICS WINDOW:
!--------------------------------------------------------------------
      do j=1,jsize               
        do i=1,isize                 

!--------------------------------------------------------------------
!    if it is already known that convection is not possible in this 
!    column, cycle to end of this loop and process the next column.
!--------------------------------------------------------------------
          if (exit_flag(i,j)) cycle

!--------------------------------------------------------------------
!    determine if column diagnostics are requested for this column.
!    define the output unit and set debug_ijt to .true. if it is.
!--------------------------------------------------------------------
          debug_ijt = .false.
          diag_unit = -99
          if (Col_diag%in_diagnostics_window ) then
            do n=1,Col_diag%ncols_in_window
              if (j == Col_diag%j_dc(n) .and.      &
                  i == Col_diag%i_dc(n)) then
                debug_ijt = .true.
                diag_unit = Col_diag%unit_dc(n)
                exit
              endif
            end do
          endif

!---------------------------------------------------------------------
!    define an inverted interface level pressure profile phalf_c 
!    (level 1 at the surface).
!---------------------------------------------------------------------
          do k=1,nlev_lsm+1
            phalf_c(k) = phalf(i,j,nlev_lsm+2-k)
          end do

!--------------------------------------------------------------------
!    call don_d_integ_cu_ensemble_k to determine the 
!    characteristics of the clouds in the cumulus ensemble defined in 
!    the current column.
!--------------------------------------------------------------------
          call don_d_integ_cu_ensemble_k             &
               (nlev_lsm, nlev_hires, ntr, me, diag_unit, debug_ijt, &
!++lwh
                Param, Col_diag, Nml, Initialized, Don_cape%model_t(i,j,:), &
!--lwh
                Don_cape%model_r(i,j,:), Don_cape%model_p(i,j,:),  &
                phalf_c, xgcm_v(i,j,:,:), sfc_sh_flux(i,j),     &
                sfc_vapor_flux(i,j), sfc_tracer_flux(i,j,:),   &
                Don_cape%plzb(i,j), exit_flag(i,j),                &
                ensmbl_precip, ensmbl_cond, ensmbl_anvil_cond, pb,  &
                pt_ens, ampta1, Don_conv%amax(i,j), emsm, rlsm,  &
                cld_press, ensmbl_melt, ensmbl_freeze, disb, disc, disd,&
                disg, enctf, encmf, enev, ecds, eces, ensmbl_cloud_area,&
                cuq, cuql_v, detmfl, uceml, qtren, etsm, lmeso, ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (trim(ermesg) /= ' ') return

!--------------------------------------------------------------------
!    if the exit_flag was set within integrate_cumulus_ensemble (due to
!    an ensemble member either a) not reaching an acceptable level of 
!    free convection, b) not producing precipitation, c) having con-
!    densate evaporation within the cloud, or d) not having a net column
!    non-zero moisture forcing (the "moisture constraint") stop the 
!    calculations for this column -- deep convection is turned off here,
!    output fields will reflect the absence of the effects of deep 
!    convection in this column.
!--------------------------------------------------------------------
          if (exit_flag(i,j)) cycle

!--------------------------------------------------------------------
!    if mesoscale circulation is present, call subroutine meso_effects 
!    to obtain full ensemble output fields to be applied to large-scale
!    model fields.
!--------------------------------------------------------------------
          if (lmeso) then
            call don_m_meso_effects_k  &
                 (nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, &
                  Param, Don_cape%model_p(i,j,:),   &
                  Don_cape%model_t(i,j,:), Don_cape%model_r(i,j,:),  &
                  phalf_c, rlsm, emsm, etsm, xgcm_v(i,j,:,:),   &
                  ensmbl_cond, ensmbl_precip, pb, Don_cape%plzb(i,j), &
                  pt_ens, ampta1, ensmbl_anvil_cond, wtp, qtmes,    &
                  anvil_precip_melt, meso_cloud_area, cmus_tot, dmeml, &
                  emds, emes, wmms, wmps, umeml, tmes, mrmes,    &
                  Don_conv%emdi_v(i,j), Don_conv%pmd_v(i,j),   &
                  Don_conv%pztm_v(i,j), Don_conv%pzm_v(i,j),    &
                  Don_conv%meso_precip(i,j), ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
            if (trim(ermesg) /= ' ') return

!--------------------------------------------------------------------
!    define cmus_tot   as the profile of total condensate source to the
!    large-scale flow from the mesoscale circulation; the sum of the 
!    water mass condensed in the mesoscale updraft plus the vapor 
!    transferred from cell to mesoscale and then condensed. 
!--------------------------------------------------------------------
          else
            qtmes = 0.
            wtp = 0.
            umeml = 0.
            dmeml = 0.
            cmus_tot = 0.
            tmes = 0.
            wmms = 0.
            wmps = 0.
            mrmes = 0.
            emds = 0.
            emes = 0.
            anvil_precip_melt = 0.
            meso_cloud_area = 0.
          endif

!---------------------------------------------------------------------
!    if in a diagnostics column, output the profiles of cell-scale 
!    tracer flux convergence (qtren). 
!---------------------------------------------------------------------
          if (debug_ijt) then
            do k=1,nlev_lsm
              do kcont=1,ntr  
                if (qtren(k,kcont) /= 0.00) then
                  write (diag_unit, '(a, 2i4, f19.10, e20.12)')  &
                  'in mulsub: jk, pr,qtren= ', k, kcont,              &
                            Don_cape%model_p(i,j,k), qtren(k,kcont)
                endif
              end do
            end do
          endif

!--------------------------------------------------------------------
!    if in diagnostics column, output the rate of condensate transfer 
!    from cells to anvil (ensmbl_anvil_cond), and the ratio of
!    convective precipitation to total precipitation (contotxx_v).
!--------------------------------------------------------------------
          if (debug_ijt) then
            write (diag_unit, '(a,e20.12, a, e20.12)')  &
              'in mulsub: CATOT= ',ensmbl_anvil_cond,' contot=',  &
                       ensmbl_precip/(ensmbl_precip +    &
                                              Don_conv%meso_precip(i,j))
          endif

!----------------------------------------------------------------------
!    call subroutine define_convective_forcing to combine the cell and
!    mesoscale contributions to the output fields and the time tendency
!    terms that will be returned to the large-scale model. it also
!    call subroutine output_diagnostic_profiles to print various 
!    output fields from the donner_deep parameterization in those 
!    columns for which diagnostics have been requested.
!----------------------------------------------------------------------
          call don_d_def_conv_forcing_k   &
               (nlev_lsm, diag_unit, debug_ijt, lmeso, Param, &
                ensmbl_precip, Don_conv%meso_precip(i,j), &
                meso_cloud_area, anvil_precip_melt, phalf_c, enev,  &
                encmf, ensmbl_freeze, enctf, disg, ecds, eces, emds,   &
                emes, mrmes, tmes, wmps, ensmbl_cloud_area, ensmbl_melt,&
                Don_cape%model_p(i,j,:), Don_cape%model_t(i,j,:),  &
                cmus_tot, wmms, disc, disb, disd, total_precip(i,j), &
                disf, disg_2, disn, dise, disa, cutotal, temp_tend_melt,&
                temp_tend_freeze, ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    call finalize_output_fields to convert to mks units and then store
!    profile arrays into the various components of the donner_conv type
!    derived-type variable Don_conv.
!----------------------------------------------------------------------
          call don_d_finalize_output_fields_k  &
               (nlev_lsm, ntr, i, j, Param, disb, disc,   &
                temp_tend_freeze, temp_tend_melt, tmes, disd, cmus_tot, &
                ecds, eces, emds, emes, wmms, wmps, mrmes, cutotal,   &
                dmeml, detmfl, uceml, umeml, cuq, cuql_v, qtren, qtmes,&
                wtp, Don_conv, ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (trim(ermesg) /= ' ') return

!--------------------------------------------------------------------
!    store some additional output fields in the donner_conv type 
!    variable for later use.
!--------------------------------------------------------------------
          Don_conv%cell_precip(i,j) = ensmbl_precip
          Don_conv%pb_v(i,j) = pb
          Don_conv%ampta1(i,j) = ampta1
          dise_v(i,j,:) = dise(:)/(1.0E03*Param%SECONDS_PER_DAY)
          disa_v(i,j,:) = disa(:)/Param%SECONDS_PER_DAY

          do k=1,nlev_lsm
            kinv = nlev_lsm + 1 - k
            temperature_forcing(i,j,kinv) = disa(k)/Param%SECONDS_PER_DAY
            moisture_forcing(i,j,kinv)    = dise(k)/    &
                                           (1.0e03*Param%SECONDS_PER_DAY)
          end do

!--------------------------------------------------------------------
!    for any diagnostic columns in the window in which deep convection
!    occurred, output the cloud anvil area (Don_conv%ampta1) and the
!    total precipitation produced (total_precip). also output the vert-
!    ical profile of total cloud fraction (Don_conv%cual).
!--------------------------------------------------------------------
          if (debug_ijt) then
            write  (diag_unit, '(a, 2e20.12)')   &
                  'in cupar:  ampt,tpre= ',  &
                            Don_conv%ampta1(i,j), total_precip(i,j)      
            do k=1,nlev_lsm-Col_diag%kstart+1    
              write (diag_unit, '(a, i4, e20.12)')  &
                   'in cupar: k,cual= ',k,  &
                                Don_conv%cual(i,j,nlev_lsm-k+1)
            end do
          endif

!---------------------------------------------------------------------
!    define the time rates of change of column-integrated water vapor
!    (dqls_v) and the time rate of change needed to deplete the column
!    water vapor in a single donner timestep (qlsd_v).
!---------------------------------------------------------------------
          dqls_v = (Don_cape%qint(i,j) - Don_cape%qint_lag(i,j))/dt
          qlsd_v = Don_cape%qint(i,j)/Nml%donner_deep_freq
          max_depletion_rate = dqls_v + qlsd_v

!--------------------------------------------------------------------
!    if in a diagnostic column, output these moisture tendency 
!    variables.
!--------------------------------------------------------------------
          if (debug_ijt) then
            write (diag_unit, '(a, 2e20.12)')   &
                  'in cupar: dqls,qlsd= ', dqls_v, qlsd_v     
          endif

!---------------------------------------------------------------------
!    call determine_cloud_area to define the cloud area of the convect-
!    ive clouds and so close the parameterization. note that exit_flag
!    may be set to .true. within determine_cloud_area, so that the 
!    if (exit_flag) loop must be closed after this call.
!---------------------------------------------------------------------
          if (.not. exit_flag(i,j)) then
            call don_d_determine_cloud_area_k  &
                 (me, nlev_lsm, nlev_hires, diag_unit, debug_ijt, Param,&
                  Nml, max_depletion_rate, Don_conv%dcape(i,j),   &
                  Don_conv%amax(i,j), dise_v(i,j,:), disa_v(i,j,:),     &
                  Don_cape%model_p(i,j,:), Don_cape%model_t(i,j,:), &
                  Don_cape%model_r(i,j,:), Don_cape%env_t(i,j,:), &
                  Don_cape%env_r(i,j,:), Don_cape%parcel_t(i,j,:), &
                  Don_cape%parcel_r(i,j,:), Don_cape%cape_p(i,j,:), &
                  exit_flag(i,j), Don_conv%amos(i,j), Don_conv%a1(i,j), &
                  ermesg)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
            if (trim(ermesg) /= ' ') return
          endif

!--------------------------------------------------------------------
!    if the exit_flag was set within determine_cloud_area (due to
!    not having a net column non-zero moisture forcing (the "moisture 
!    constraint") set a flag to indicate that deep convection is turned
!    off; output fields will be made to reflect the absence of the 
!    effects of deep convection in this column. 
!--------------------------------------------------------------------
        end do
      end do

!--------------------------------------------------------------------



end subroutine don_d_mulsub_k



!######################################################################

subroutine don_d_integ_cu_ensemble_k             &
         (nlev_lsm, nlev_hires, ntr, me, diag_unit, debug_ijt, Param,   &
!++lwh
          Col_diag, Nml, Initialized, temp_c, &
!--lwh
          mixing_ratio_c, pfull_c, phalf_c,   &
          tracers_c, sfc_sh_flux_c, sfc_vapor_flux_c,   &
          sfc_tracer_flux_c, plzb_c, exit_flag_c, ensmbl_precip,    &
          ensmbl_cond, ensmbl_anvil_cond, pb, pt_ens, ampta1, amax, &
          emsm, rlsm, cld_press, ensmbl_melt, ensmbl_freeze, disb, disc,&
          disd, disg, enctf, encmf, enev, ecds, eces, ensmbl_cloud_area,&
          cuq, cuql_v, detmfl, uceml, qtren, etsm, lmeso, ermesg)

!----------------------------------------------------------------------
!    subroutine integrate_cumulus_ensemble works on a single model 
!    column. all profile arrays used in this subroutine and below have 
!    index 1 nearest the surface. it first determines the lifting conden-
!    sation level (if one exists) of a parcel moving from the specified 
!    parcel_launch_level. if an lcl is found, subroutine 
!    donner_cloud_model_cloud_model is called to determine the behavior
!    of each of kpar cloud ensemble mem-
!    bers assumed present in the column (each ensemble member is ass-
!    umed to have a different entrainment rate). if all ensemble members
!    produce deep convection, the ensemble statistics are produced for 
!    use in the large-scale model; otherwise deep convection is not seen
!    in the large-scale model in this grid column. if the ensemble will 
!    support a mesoscale circulation, its impact on the large-scale model
!    fields is also determined. upon completion, the appropriate output 
!    fields needed by the large-scale model are returned to the calling 
!    routine.
!----------------------------------------------------------------------
use messy_convect_donner_types_mod, only : donner_param_type, &
                             donner_nml_type, donner_column_diag_type, &
!++lwh
                             donner_initialized_type
!--lwh
USE MESSY_CONVECT_DONNER_UTIL,       ONLY: don_u_map_hires_i_to_lores_c_k
USE MESSY_CONVECT_DONNER_CLOUDMODEL, ONLY: don_cm_lcl_k,                  &
                                           don_cm_cloud_model_k,          &
                                           don_cm_mesub_k

implicit none 

!----------------------------------------------------------------------
integer,                           intent(in)    :: nlev_lsm,    &
                                                    nlev_hires, ntr, &
                                                    me, diag_unit
logical,                           intent(in)    :: debug_ijt
type(donner_param_type),           intent(in)    :: Param
type(donner_column_diag_type),     intent(in)    :: Col_diag
type(donner_nml_type),             intent(in)    :: Nml   
!++lwh
type(donner_initialized_type),     intent(in)    :: Initialized
!--lwh
real,    dimension(nlev_lsm),      intent(in)    :: temp_c,   &
                                                    mixing_ratio_c,   &
                                                    pfull_c
real,    dimension(nlev_lsm+1),    intent(in)    :: phalf_c
real,    dimension(nlev_lsm,ntr),  intent(in)    :: tracers_c           
real,                              intent(in)    :: sfc_sh_flux_c,   &
                                                    sfc_vapor_flux_c 
real,    dimension(ntr),           intent(in)    :: sfc_tracer_flux_c 
real,                              intent(in)    :: plzb_c
logical,                           intent(inout) :: exit_flag_c  
real,                              intent(out)   :: ensmbl_precip,   &
                                                    ensmbl_cond,&
                                                    ensmbl_anvil_cond, &
                                                    pb, pt_ens, ampta1, &
                                                    amax
real,    dimension(nlev_hires),    intent(out)   :: emsm, rlsm, cld_press
real,    dimension(nlev_lsm),      intent(out)   :: ensmbl_melt,   &
                                                    ensmbl_freeze,&
                                                    disb, disc, disd, &
                                                    disg, enctf, encmf, &
                                                    enev, ecds, eces, &
                                                    ensmbl_cloud_area, &
                                                    cuq, cuql_v, &
                                                    detmfl, uceml
real,    dimension(nlev_lsm,ntr),  intent(out)   :: qtren
real,    dimension(nlev_hires,ntr),intent(out)   :: etsm
logical,                           intent(out)   :: lmeso       
character(len=*),                  intent(out)   :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
! 
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     diag_unit      unit number for column diagnostics output, if 
!                    diagnostics are requested for the current column
!     debug_ijt      logical indicating whether current column requested
!                    column diagnostics
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     temp_c         temperature field at model full levels 
!                    index 1 nearest the surface [ deg K ]
!     mixing_ratio_c        vapor mixing ratio at model full levels 
!                    index 1 nearest the surface
!                    [ kg(h2o) / kg(dry air) ]
!     pfull_c         pressure field at large-scale model full levels 
!                    index 1 nearest the surface [ Pa ]
!     phalf_c        pressure field at large-scale model half-levels 
!                    index 1 nearest the surface [ Pa ]
!     tracers_c      tracer fields that are to be transported by donner
!                    convection.  index 1 nearest the surface 
!                    [ kg (tracer) / kg (dry air) ]
!     sfc_sh_flux_c  sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux_c water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     sfc_tracer_flux_c  
!                    flux across the surface of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     plzb_c         level of zero buoyancy for a parcel lifted from
!                    the parcel_launch_level.  [ Pa ]
!
!   intent(inout) variables:
!
!     exit_flag_c    logical indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    current model column 
!
!   intent(out) variables:
!    
!     ensmbl_precip      sum of precipitation rate over ensemble members,
!                        # 1 to the current, weighted by the area at 
!                        cloud base of each member
!                        [ mm / day ]
!     ensmbl_cond        sum of condensation rate over ensemble members,
!                        # 1 to the current, weighted by the area at 
!                        cloud base of each member
!                        [ mm / day ]
!     ensmbl_anvil_cond  sum of rate of transfer of condensate from cell 
!                        to anvil over ensemble members, # 1 to the c
!                        current, weighted by the area at cloud base of 
!                        each member [ mm / day ]
!     pb                 pressure at cloud base for ensemble (all ensem-
!                        ble members have same base) [ Pa ]
!     pt_ens             pressure at cloud top for the ensemble (top 
!                        pressure of deepest ensemble member) [ Pa ]
!     ampta1             cloudtop anvil area (assumed to be five times
!                        larger than the sum of the cloud top areas of 
!                        the ensemble members, as in Leary and Houze 
!                        (1980).  [ fraction ]
!     amax               maximum allowable area of cloud base that is
!                        allowed; if cloud base area is larger than 
!                        amax, the cloud fractional area somewhere in
!                        the grid box would be greater than one, which 
!                        is non-physical.
!     emsm               vertical profile on the hi-res grid of vertical
!                        moisture flux convergence, summed over ensemble 
!                        members # 1 to the current, each member's cont-
!                        ribution being weighted by its cloud area at 
!                        level k relative to the cloud base area of 
!                        ensemble member #1  
!                        [ kg (h2o) / ( kg(dry air) sec ) ]
!     rlsm               vertical profile on the hi-res grid of conden-
!                        sation rate, summed over ensemble members # 1 to
!                        the current, each member's contribution being 
!                        weighted by its cloud area at level k relative 
!                        to the cloud base area of ensemble member #1
!                        [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!     cld_press          pressures at hi-res model levels [ Pa ]
!     ensmbl_melt        vertical profile on the lo-res grid of ice melt,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!     ensmbl_freeze      vertical profile on the lo-res grid of freezing,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!     disg               vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation
!                        associated with the evaporation of condensate
!                        in the convective downdraft and updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!     enev               vertical profile on the lo-res grid of the      
!                        cloud-area-weighted profile of the potential
!                        cloud water evaporation, summed over ensemble 
!                        members # 1 to the current, each member's con-
!                        tribution being weighted by its cloud area at !
!                        level k relative to the cloud base area of 
!                        ensemble member #1.  this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ ( kg(h2o) ) / ( kg(dry air) sec ) ] 
!     enctf              vertical profile on the lo-res grid of the entr-
!                        opy forcing, consisting of the sum of the
!                        vertical entropy flux convergence and the latent
!                        heat release, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area
!                        of ensemble member #1
!                        [ deg K / day ]                        
!     encmf              vertical profile on the lo-res grid of the      
!                        moisture forcing, consisting of the sum of the
!                        vertical moisture flux convergence and the cond-
!                        ensation, summed over ensemble members # 1 to 
!                        the current, each member's contribution being 
!                        weighted by its cloud area at level k relative 
!                        to the cloud base area of ensemble member #1
!                        [ ( kg(h2o) ) / ( kg( dry air) day ) ] 
!     disb               vertical profile on the lo-res grid of the      
!                        temperature flux convergence, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area 
!                        of ensemble member #1.  
!                        [ deg K / day ] 
!     disc               vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation, 
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!     disd               vertical profile on the lo-res grid of the      
!                        vertical moisture flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        lo-res grid for the current ensemble member 
!                        [  g(h2o) / ( kg(dry air) day ) ]
!     ecds               vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective downdraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!     eces               vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!     ensmbl_cloud_area  total cloud area profile over all ensemble
!                        members on large_scale model grid [ fraction ]
!     cuq                ice water profile on large-scale model grid, 
!                        normalized by ensemble cloud area.
!     cuql_v             liquid water profile on large-scale model grid, 
!                        normalized by ensemble cloud area.
!     uceml              upward mass flux on large_scale model grid     
!                        [ kg (air) / (sec m**2) ]
!     detmfl             detrained mass flux on large-scale model grid
!                        normalized by ensemble cloud area
!                        [ kg (air) / (sec m**2) ]
!     etsm               vertical profile on the hi-res grid of vertical
!                        tracer flux convergence, summed over ensemble 
!                        members # 1 to the current, each member's con-
!                        tribution being weighted by its cloud area at i
!                        level k relative to the cloud base area of 
!                        ensemble member #1 
!                        [ kg (tracer) / ( kg(dry air) sec ) ]
!     qtren              vertical profile on the lo-res grid of the      
!                        vertical tracer flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!     lmeso              logical variable; if .false., then it has been
!                        determined that a mesoscale circulation cannot
!                        exist in the current column. final value not
!                        determined until all ensemble members have been
!                        integrated. 
!     ermesg             character string containing any error message
!                        that is returned from a kernel subroutine
!
!---------------------------------------------------------------------

!     cmui             normalized vertical integral of mesoscale-updraft
!                      deposition (kg(H2O)/((m**2) sec)
!     cmus(nlev)       normalized mesoscale-updraft deposition
!                      (kg(H2O)/kg/sec)
!     emds(nlev)       normalized mesoscale-downdraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emei             normalized vertical integral of mesoscale-updraft
!                      sublimation (kg(h2O)/((m**2) sec)
!     emes(nlev)       normalized mesoscale-updraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     disa(nlev)       normalized thermal forcing, cells+meso (K/sec)
!                      (excludes convergence of surface heat flux)
!                      index 1 at ground. Cumulus thermal forcing defined
!                      as in Fig. 3 of Donner (1993, JAS).
!     disb(nlev)       normalized cell entropy-flux convergence (K/sec)
!                      (excludes convergence of surface flux)
!                      index 1 at ground. Entropy-flux convergence divided
!                      by (p0/p)**(rd/cp).
!     disc(nlev)       normalized cell condensation/deposition
!                      (K/sec)
!                      index 1 at ground.
!     disd(nlev)       normalized cell moisture-flux convergence
!                      (excludes convergence of surface moisture flux)
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dise(nlev)       normalized moisture forcing, cells+meso (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dmeml(nlev)      mass flux in mesoscale downdraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!     elt(nlev)        normalized melting (K/sec)
!                      index 1 at ground.
!     fre(nlev)        normalized freezing (K/sec)
!                      index 1 at ground.
!     pb               pressure at base of cumulus updrafts (Pa)
!     pmd              pressure at top of mesoscale downdraft (Pa)
!     pztm             pressure at top of mesoscale updraft (Pa)
!     mrmes(nlev)       normalized mesoscale moisture-flux convergence
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     qtmes(nlev,ncont)  tracer tendency due to mesoscale tracer-flux
!                        convergence (kg/kg/s) (normalized by a(1,p_b))
!                        index 1 at ground 
!     qtren_v          normalized tracer tendency due to cells...
!                      (lon,lat,vert,tracer index)
!                      Vertical index increases as height increases.
!     sfcq(nlev)       boundary-layer mixing-ratio tendency due to surface
!                      moisture flux (kg(H2O)/kg/sec)
!     sfch(nlev)       boundary-layer heating due to surface heat flux
!                      (K/sec)
!     tmes(nlev)       normalized mesoscale entropy-flux convergence
!                      (K/sec)
!                      Entropy-flux convergence is mesoscale component
!                      of second term in expression for cumulus thermal
!                      forcing in Fig. 3 of Donner (1993, JAS).
!                      index 1 at ground.
!     tpre_v           total normalized precipitation (mm/day)
!     detmfl(nlev)     detrained mass flux from cell updrafts
!                      (normalized by a(1,p_b))
!                      (index 1 near atmosphere bottom)
!                      (kg/((m**2)*s)
!     uceml(nlev)      normalized mass fluxes in cell updrafts
!                      (kg/((m**2)*s) 
!     umeml(nlev)      mass flux in mesoscale updraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!                      index 1 at ground.
!     wmms(nlev)       normalized mesoscale deposition of water vapor from
!                      cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wmps(nlev)       normalized mesoscale redistribution of water vapor
!                      from cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wtp_v            tracer redistributed by mesoscale processes
!                      (kg/kg/s) (normalized by a(1,p_b))
!                      vertical index increases with increasing height
!                      (lon,lat,vert,tracer index)
!--------------------------------------------------------------------


!!  UNITS
!!    ucemh  [kg /sec / m**2 ]
!!    detmfh [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    h1_2   [ deg K / sec ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]




!        cape     convective available potential energy (J/kg)
!        cin      convective inhibtion (J/kg)
!        cpd      specific heat of dry air at constant pressure (J/(kg K))
!        cpv      specific heat of water vapor [J/(kg K)]
!        dcape    local rate of CAPE change by all processes
!                 other than deep convection [J/(kg s)]
!        dqls     local rate of change in column-integrated vapor
!                 by all processes other than deep convection
!                 {kg(H2O)/[(m**2) s]}
!        epsilo   ratio of molecular weights of water vapor to dry air
!        gravm    gravity constant [m/(s**2)]
!        ilon     longitude index
!        jlat     latitude index
!        mcu      frequency (in time steps) of deep cumulus
!        current_displ  integrated low-level displacement (Pa)
!        cape_p   pressure at Cape.F resolution (Pa)
!                 Index 1 at bottom of model.
!        plfc     pressure at level of free convection (Pa)
!        plzb_c   pressure at level of zero buoyancy (Pa)
!        pr       pressure at Skyhi vertical resolution (Pa)
!                 Index 1 nearest ground  
!        q        large-scale vapor mixing ratio at Skyhi vertical resolution
!                 [kg(h2O)/kg]
!                 Index 1 nearest ground 
!        qlsd     column-integrated vapor divided by timestep for cumulus
!                 parameterization {kg(H2O)/[(m**2) s]}
!        r        large-scale vapor mixing ratio at Cape.F resolution
!                 [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rpc      parcel vapor mixing ratio from Cape.F [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rd       gas constant for dry air (J/(kg K))
!        rlat     latent heat of vaporization (J/kg)
!        rv       gas constant for water vapor (J/(kg K))
!        t        large-scale temperature at Skyhi vertical resolution (K)
!                 Index 1 nearest ground
!        tcape    large-scale temperature at Cape.F resolution (K)
!                 Index 1 at bottom of model.
!        tpc      parcel temperature from from Cape.F (K)
!                 Index 1 at bottom of model.
!

!----------------------------------------------------------------------
!   local variables:

      real,    dimension (nlev_hires)     ::                &
              efchr, emfhr, te, mre, rcl, dpf, qlw, dfr, cfracice, alp, &
              cld_evap, flux, ucemh, cuql, cuqli, detmfh

      real,    dimension (nlev_lsm)       ::           &
              h1, q1, pi, em, rlh, cmf, cell_freeze, cell_melt, disf,  &
              meso_melt, meso_freeze, h1_2, disg_2, out, evap_rate, ecd,&
              ece, disl, thlr, qlr, sfcq, sfch

      real,    dimension (nlev_hires,ntr) :: xclo, xtrae, etfhr
      real,    dimension (nlev_lsm,ntr)   :: qtr
      real,    dimension (Param%kpar)     :: cuto, preto, ptma
      integer, dimension (Param%kpar)     :: ncca

      logical ::   lcl_reached                  
      integer ::   ncc_kou, ncc_ens
      integer ::   k,    kou
      integer ::   kc, kcl, kch
      real    ::   al, dints, disga, dp, mrb, pmel, p, sumehf, sumhlr, &
                   summel, pl, dpp, ph, esh, esl, rh, rl, pkc, tveh,   &
                   tvch, dpdzh, ehfh, tvel, tvcl, dpdzl, ehfl, ptt,   &
                   ehf, tve, tvc, dpdz, exf, emfh, emfl, thetf, emff, &
                   sbl, p1, dmela, psmx, esumc, sumf, summ, sumqme,   &
                   sumg, sumn, sumelt, sumfre, summes, esum, sumev,   &
                   esuma, sumlhr, es, etfh, etfl, dint, cu, cell_precip,&
                   precip, conint, ca, apt, qtrsum, qtmesum, rintsum, &
                   rintsum2, intgl_in, intgl_out, alphaw, tb, alpp,   &
                   pcsave, rsc, ensmbl_cld_top_area  

!----------------------------------------------------------------------
!   local variables:
!
!      ensmbl_cld_top_area  
!                       sum of the cloud top areas over ensemble members 
!                       # 1 to the current, normalized by the cloud base
!                       area of ensemble member # 1 [ dimensionless ]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!    if in diagnostics column, output the large-scale model temperature,
!    vapor mixing ratio and full-level pressure profiles (index 1 near-
!    est the surface).
!---------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm-Col_diag%kstart+1
          write (diag_unit, '(a, i4, f20.14, e20.12, f19.10)')&
                'in mulsub: k,T,Q,P= ',k, temp_c(k), mixing_ratio_c(k), pfull_c(k)
        end do
      endif

!--------------------------------------------------------------------
!    call don_cm_lcl_k to calculate the temperature (tb), a
!    pressure (pb) and mixing ratio (mrb) at the lifting condensation 
!    level for a parcel starting from the parcel_launch_level. if a sat-
!    isfactory lcl is not reached for this parcel, the logical variable 
!    lcl_reached will be set to .false..
!--------------------------------------------------------------------
      call don_cm_lcl_k    &
           (Param, temp_c (Nml%parcel_launch_level),    &
            pfull_c       (Nml%parcel_launch_level),    &
            mixing_ratio_c(Nml%parcel_launch_level),   &
            tb, pb, mrb, lcl_reached, ermesg)     

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!--------------------------------------------------------------------
!    if in diagnostics column and an lcl was defined, output the lcl 
!    temperature, pressure and mixing ratio. if an acceptble lcl was 
!    not reached, print a message.
!--------------------------------------------------------------------
      if (debug_ijt) then
        if (lcl_reached) then
          write (diag_unit, '(a, f20.14, f19.10, e20.12)') &
                                'in mulsub: tb,pb,qb= ',tb, pb, mrb  
        else
          write (diag_unit, '(a)') 'in mulsub: lcl not reached'
        endif
      endif

!--------------------------------------------------------------------
!    if an acceptable lcl was not reached, set exit_flag_c so that the
!    remaining computations for this column are bypassed, and return to
!    calling routine. 
!--------------------------------------------------------------------
      if (.not. lcl_reached) then
        exit_flag_c = .true.
        return
      endif
 
!---------------------------------------------------------------------
!    if calculations are continuing, initialize needed variables.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize variables which will accumulate scalar sums over all 
!    ensemble members.
!---------------------------------------------------------------------
      ensmbl_precip       = 0.
      ensmbl_cond         = 0.
      ensmbl_anvil_cond   = 0.
      ensmbl_cld_top_area = 0.

!---------------------------------------------------------------------
!    initialize the variables which will contain the sum over the 
!    ensemble members of the vertical profiles of various quantities 
!    on the cloud-model grid.
!---------------------------------------------------------------------
      do k=1,nlev_hires
        cuql(k)   = 0.
        cuqli(k)  = 0.
        ucemh(k)  = 0.
        detmfh(k) = 0.
        alp(k)    = 0.
        rlsm(k)   = 0.
        emsm(k)   = 0.
        etsm(k,:) = 0.
      end do

!---------------------------------------------------------------------
!    initialize the variables which will contain the sum over the 
!    ensemble members of the vertical profiles of various quantities 
!    on the large-scale model grid.
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        ensmbl_freeze(k)    = 0.
        ensmbl_melt(k)    = 0.
        disb(k)    = 0.
        disc(k)    = 0.
        disd(k)    = 0.
        ecds(k)    = 0.
        eces(k)    = 0.
        enctf(k)   = 0.
        encmf(k)   = 0.
        disg(k)    = 0.
        enev(k)    = 0.
        qtren(k,:) = 0.
      end do

      evap_rate = 0.

!--------------------------------------------------------------------
!    initialize a logical variable which will indicate whether a
!    mesoscale circulation is present in this column. this may be 
!    precluded via the nml variable allow_mesoscale_circulation. if any
!    ensemble members are unable to support a mesoscale circulation, 
!    lmeso will be set to .false. within the following loop over the kpar
!    ensemble members. if the first member of the ensemble (the most 
!    entraining) can, then it is likely (but not guaranteed) that the 
!    ensemble will be able to.
!--------------------------------------------------------------------
      if (Nml%allow_mesoscale_circulation) then
        lmeso = .true.
      else
        lmeso = .false.
      endif

!--------------------------------------------------------------------
!    define the array of cloud model pressure levels (cld_press).
!--------------------------------------------------------------------
      do k=1,nlev_hires
        cld_press(k) = pb + (k-1)*Param%dp_of_cloud_model
      end do

!--------------------------------------------------------------------
!    if this is the first ensemble member, initialize the variables
!    which are defined on this call and will be used by the other
!    ensemble members.
!--------------------------------------------------------------------
      pcsave = phalf_c(1)

!--------------------------------------------------------------------
!    loop over the KPAR members of the cumulus ensemble.
!--------------------------------------------------------------------
      do kou=1,Param%kpar

!-------------------------------------------------------------------
!    define the appropriate entrainment factor (alpp) for this ensemble
!    member using values based on observations either obtained from
!    the GATE or KEP studies.
!-------------------------------------------------------------------
        if (trim(Nml%entrainment_constant_source) == 'gate') then
          alpp = Param%max_entrainment_constant_gate/  &
                           Param%ensemble_entrain_factors_gate(kou)
        else if (trim(Nml%entrainment_constant_source) == 'kep') then
          alpp = Param%max_entrainment_constant_kep/  &
                           Param%ensemble_entrain_factors_kep(kou)
        else
          ermesg = 'invalid entrainment_constant_source'
          return
        endif

!--------------------------------------------------------------------
!    call cloud_model to obtain the in-cloud and environmental profiles
!    and fluxes and column integrals associated with this ensemble 
!    member.
!--------------------------------------------------------------------
        call don_cm_cloud_model_k   &
             (nlev_lsm, nlev_hires, ntr, kou, diag_unit, debug_ijt,   &
!++lwh
!              Param, Col_diag, Initialized, tb, pb, alpp, cld_press, temp_c, &
              Param, Col_diag, tb, pb, alpp, cld_press, temp_c, &
!--lwh
              mixing_ratio_c, pfull_c, phalf_c, tracers_c, pcsave,  &
              exit_flag_c, rcl, dpf, qlw, dfr, flux, ptma(kou), &
              dint, cu, cell_precip, dints, apt, cell_melt, efchr, &
              emfhr, cfracice, etfhr, ncc_kou, ermesg)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return

!--------------------------------------------------------------------
!    if the cloud thickness is less than pdeep_mc, it will not
!    support a mesoscale circulation. set a logical flag to indicate
!    the absence of a mesoscale component for this column's cloud
!    ensemble.
!--------------------------------------------------------------------
        if ((pb - ptma(kou)) < Param%pdeep_mc)  then
          lmeso = .false.
        else
          cell_melt(:) = 0.0
        endif


        if (exit_flag_c) return

!--------------------------------------------------------------------
!    if calculations are continuing, 
!--------------------------------------------------------------------

!----------------------------------------------------------------------
!    if in diagnostics column, output the cloud base (pb) and cloud top
!    (ptma) pressures, and the mesoscale circulation logical variable
!    (lmeso).
!----------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 2f19.10,l4)')    &
                     'in mulsub: PB,PT, lmeso= ', pb, ptma(kou), lmeso
        endif

!---------------------------------------------------------------------
!    define the cloud water from this ensemble member which must be 
!    evaporated if it turns out that there is no mesoscale circulation 
!    associated with the ensemble.
!---------------------------------------------------------------------
        cld_evap(:) = -dpf(:)*(1. - (cell_precip/cu))

!---------------------------------------------------------------------
!    define the pressure one cloud model level above cloud top (ptt).
!---------------------------------------------------------------------
        ptt = ptma(kou) + Param%dp_of_cloud_model

!----------------------------------------------------------------------
!    call define_lo_res_model_profiles to map profiles generated on the
!    cloud model grid to the vertical grid of the large-scale model for
!    this ensemble member.
!----------------------------------------------------------------------
        call don_d_def_lores_model_profs_k        &
             (nlev_lsm, nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt, &
              Param, pb, ptt, sfc_vapor_flux_c, sfc_sh_flux_c,  &
              sfc_tracer_flux_c, pfull_c, phalf_c, cld_press, dpf, dfr, &
              cld_evap, qlw, emfhr, efchr, etfhr, cell_freeze,     &
              evap_rate, h1, h1_2, q1, qtr, ermesg)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    if this member of the ensemble supports a mesoscale circulation,
!    call mesub to obtain various terms related to moving condensate
!    from the convective tower into the mesoscale anvil for this member.
!---------------------------------------------------------------------
        if (lmeso) then
          call don_cm_mesub_k     &
               (nlev_lsm, me, diag_unit, debug_ijt, Param, cu,   &
                cell_precip, dints, plzb_c, pb, ptma(kou), temp_c,  &
                phalf_c, ca, ecd, ece, meso_freeze, meso_melt, ermesg)
        endif

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    call don_d_add_to_ensmbl_sum_hires_k to add this member's 
!    contribution to those fields on the cloud model grid that are being
!    summed over all ensemble members.
!---------------------------------------------------------------------
        call don_d_add_to_ensmbl_sum_hires_k    &
             (nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt, &
              Param%arat(kou), cfracice, rcl, flux, emfhr, dpf, &
              qlw, etfhr, cuql, cuqli, ucemh, alp, rlsm, emsm, detmfh, &
              etsm, ermesg)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    call don_d_add_to_ensmbl_sum_intgl_k to add this member's 
!    contribution to those integrals that are being summed over all 
!    ensemble members.
!---------------------------------------------------------------------
        call don_d_add_to_ensmbl_sum_intgl_k    &
             (diag_unit, debug_ijt, lmeso, Param%arat(kou), ca,  &
              cell_precip, cu, apt, ensmbl_precip, ensmbl_cond,   &
              ensmbl_anvil_cond, ensmbl_cld_top_area, ermesg)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    call don_d_add_to_ensmbl_sum_lores_k to add this member's 
!    contribution to those fields on the lrge-scale model grid that are 
!    being summed over all ensemble members.
!---------------------------------------------------------------------
        call don_d_add_to_ensmbl_sum_lores_k    &
             (nlev_lsm, ntr, diag_unit, debug_ijt, lmeso, Param,   &
              Param%arat(kou), dint, cell_freeze, cell_melt, temp_c,   &
              h1_2, ecd, ece, evap_rate, q1, h1, pfull_c, meso_melt, &
              meso_freeze, phalf_c, qtr, ensmbl_melt, ensmbl_freeze, &
              enctf, encmf, enev, disg, disb, disc, ecds, eces, disd, &
              qtren, ermesg)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return

!--------------------------------------------------------------------
!    save the cloud top (ptma) pressures, the total condensation (cuto),
!    total precpitation (preto) and cloud top index (ncca) from this !
!    ensemble member.
!--------------------------------------------------------------------
        cuto(kou)  = cu
        preto(kou) = cell_precip
        ncca(kou)  = ncc_kou
      end do   ! (kou loop over ensemble members)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! 31   CONTINUE
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!--------------------------------------------------------------------
!    if calculations are continuing: 
!--------------------------------------------------------------------

!----------------------------------------------------------------------
!    define ensemble cloud top pressure (pt_ens) to be the cloud top of 
!    the most penetrative ensemble member. this is frequently, but not 
!    always, the ensemble member with the lowest entrainment rate. 
!    cloud base pressure (pb) is the same for all ensemble members. 
!    define the cloud top index(ncc_ens)  as the highest of any ensemble 
!    member.
!----------------------------------------------------------------------
      pt_ens  = MINVAL (ptma)
      ncc_ens = MAXVAL (ncca)

!----------------------------------------------------------------------
!    divide the ensemble mean ice and liquid condensate terms by the 
!    total cloud area to define the average cloud water and cloud ice 
!    concentrations within the cloudy area, as opposed to averaged over 
!    the entire grid box.
!----------------------------------------------------------------------
      do k=1,ncc_ens
        if (alp(k) > 0.) then
          cuql(k)  = cuql(k)/alp(k)
          cuqli(k) = cuqli(k)/alp(k)
        endif
      end do

!---------------------------------------------------------------------
!    define the cloudtop anvil area (ampta1), assumed to be five times 
!    larger than the sum of the cloud top areas of the ensemble members,
!    as in Leary and Houze (1980), 
!---------------------------------------------------------------------
      ampta1 = 5.*ensmbl_cld_top_area

!---------------------------------------------------------------------
!    if there is no precipitation production in this column, set the 
!    inverse of the max cloud area at any layer in the column to be 0.0.
!---------------------------------------------------------------------
      if (ensmbl_precip == 0.0) then
        amax      = 0.0
      else

!---------------------------------------------------------------------
!    if there is precip in the column, determine the maximum convective 
!    cell area at any level in the column (al). the total normalized 
!    cloud area in the column (cell area + mesoscale area) cannot be 
!    greater than 1.0. this constraint imposes a limit on the cloud area
!    at cloud base (amax). this limit will be imposed in subroutine
!    determine_cloud_area. see "a bounds notes" (7/6/97).
!---------------------------------------------------------------------
        al = MAXVAL (alp)
        amax = 1./(al + ampta1)
      endif

!---------------------------------------------------------------------
!    if in diagnostics column, output the total ensemble condensation,
!    (ensmbl_cond), precipitation (ensmbl_precip), and condensate 
!    transferred into the anvil (ensmbl_anvil_cond). also output 
!    surface pressure (phalf_c(1)), ensemble cloud base nd cloud top 
!    pressures (pb, pt_ens), the flag indicating if a mesoscale circul-
!    ation is present in the grid column (lmeso), and the cloud top anvil
!    area (ampta1).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12, a, e20.12)')  &
                      'in mulsub: CUTOT=', ensmbl_cond, ' PRETOT=', &
                                      ensmbl_precip
        write (diag_unit, '(a, e20.12)') 'in mulsub: CATOT=', &
                                     ensmbl_anvil_cond
        write (diag_unit, '(a, 3f19.10, l4)')  &
              'in mulsub: ps,pb,pt,lmeso= ', phalf_c(1), pb, pt_ens, lmeso
        write (diag_unit, '(a, e20.12)')  &
                                 'in mulsub: ampt= ',ampta1     
      endif

!----------------------------------------------------------------------
!    define the pressure one level above cloud top (ptt).
!----------------------------------------------------------------------
      ptt = pt_ens + Param%dp_of_cloud_model

!--------------------------------------------------------------------
!    call define_ensemble_profiles to produce vertical profiles 
!    representing the ensemble-total cloud area (ensmbl_cloud_area), 
!    cloud liquid (cuql_v), cloud ice (cuq), mass flux(uceml) and
!    detrained mass flux (detmfl).
!--------------------------------------------------------------------
      call don_d_def_ensemble_profs_k    &
           (nlev_lsm, nlev_hires, ncc_ens, diag_unit, debug_ijt, ptt, &
            cld_press, alp, detmfh, ucemh, cuql, cuqli, phalf_c,  &
            ensmbl_cloud_area, cuql_v, cuq, detmfl, uceml, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, call error_mesg.
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    execute the following code if in a diagnostics column.
!---------------------------------------------------------------------
      if (debug_ijt) then

!----------------------------------------------------------------------
!    define the pressure at the large-scale model interface level at or 
!    just above cloud base (psmx).
!----------------------------------------------------------------------
        do k=1,nlev_lsm
          if ((phalf_c(k+1) <= pb) .and. (phalf_c(k) >= pb)) then
            psmx = phalf_c(k+1)
            exit
          endif
        end do

!----------------------------------------------------------------------
!    define the integrated boundary layer heating rate (sbl) due to the 
!    surface heat flux (sfcsf_v). it is defined in units of (deg K)/sec.
!    call don_u_map_hires_i_to_lores_c_k to distribute
!    this heating over the boundary layer.
!---------------------------------------------------------------------
        sbl = Param%grav*sfc_sh_flux_c/((phalf_c(1) - psmx)*Param%cp_air)
        write (diag_unit, '(a, e20.12, 2f19.10)')  &
             'in cm_intgl_to_gcm_col: xav,p1,p2= ',sbl, phalf_c(1), psmx 
        call don_u_map_hires_i_to_lores_c_k   &
             (nlev_lsm, sbl, phalf_c(1), psmx, phalf_c, sfch, ermesg)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return

        do k=1,size(sfch(:))
          if (sfch(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)') &
                            'in cm_intgl_to_gcm_col: k,x= ',k,sfch(k)
          endif
        end do

!----------------------------------------------------------------------
!    define the integrated boundary layer moistening rate (sbl) due to 
!    the surface moisture flux (sfcqf_v), which is defined in units of 
!    kg(h2o) per m**2 per sec. call 
!    don_u_map_hires_i_to_lores_c_k to distribute 
!    this moistening over the boundary layer.
!---------------------------------------------------------------------
        sbl = (sfc_vapor_flux_c*Param%grav)/(phalf_c(1) - psmx)
        write (diag_unit, '(a, e20.12, 2f19.10)')  &
             'in cm_intgl_to_gcm_col: xav,p1,p2= ',sbl, phalf_c(1), psmx 
        call don_u_map_hires_i_to_lores_c_k   &
             (nlev_lsm, sbl, phalf_c(1), psmx, phalf_c, sfcq, ermesg)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (trim(ermesg) /= ' ') return

        do k=1,size(sfcq(:))
          if (sfcq(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)') &
                            'in cm_intgl_to_gcm_col: k,x= ',k,sfcq(k)
          endif
        end do
      endif ! (debug_ijt)

!---------------------------------------------------------------------



end subroutine don_d_integ_cu_ensemble_k 

!#######################################################################

subroutine don_d_column_end_of_step_k  &
         (isize, jsize, nlev_lsm, ntr, Col_diag, exit_flag,   &
          total_precip, parcel_rise, temperature_forcing,&
          moisture_forcing, tracers, Don_cape, Don_conv, ermesg)       

!----------------------------------------------------------------------
!    subroutine don_d_column_end_of_step outputs the final values of
!    significant fields generated by donner_deep_mod in any columns
!    for which column diagnostics were requested, and in which deep
!    convection is present.
!----------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_cape_type, donner_conv_type, &
                             donner_column_diag_type

implicit none

!----------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize,  &
                                                   nlev_lsm, ntr
type(donner_column_diag_type),    intent(in)    :: Col_diag
logical, dimension(isize,jsize),  intent(in)    :: exit_flag
real,    dimension(isize,jsize),  intent(in)    :: total_precip,  &
                                                   parcel_rise
real,    dimension(isize,jsize,nlev_lsm),             &
                                  intent(in)    :: temperature_forcing, &
                                                   moisture_forcing
real,    dimension(isize,jsize,nlev_lsm,ntr),        &
                                  intent(in)    :: tracers
type(donner_cape_type),           intent(inout) :: Don_cape
type(donner_conv_type),           intent(inout) :: Don_conv
character(len=*),                 intent(out)   :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     exit_flag      logical variable indicating whether deep convection
!                    is present or not in each column
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     parcel_rise    accumulated vertical displacement of a 
!                    near-surface parcel as a result of the lowest
!                    model level omega field [ Pa ]
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!
!   intent(inout) variables:
!
!     Don_cape       donner_cape type derived type variable containing
!                    diagnostics related to the cape calculation assoc-
!                    iated with the donner convection parameterization
!     Don_conv       donner_convection_type derived type variable con-
!                    taining diagnostics describing the nature of the 
!                    convection produced by the donner parameterization
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer :: k, n, kcont      ! do-loop indices
      integer :: i, j, unit

      ermesg= ' '

!--------------------------------------------------------------------
!    determine if deep convection exists in any of the columns in the 
!    window for which column diagnostics were requested.
!--------------------------------------------------------------------
      do n=1,Col_diag%ncols_in_window

!--------------------------------------------------------------------
!    determine if deep convection exists in any of the columns in the 
!    window for which column diagnostics were requested. if deep
!    convection is present, output a multitude of values; if deep con-
!    vection is not present, cycle to check the next diagnostics 
!    column in the window.
!--------------------------------------------------------------------
        if (.not. exit_flag(Col_diag%i_dc(n), Col_diag%j_dc(n))) then
          i = Col_diag%i_dc(n)
          j = Col_diag%j_dc(n)
          unit = Col_diag%unit_dc(n)

!---------------------------------------------------------------------
!    output the pressures at lifting condensation level (plcl), at the
!    level of free convection (plfc), and at the level of zero buoyancy
!    (plzb).
!---------------------------------------------------------------------
          write (unit, '(a, e20.12)')  & 
                'in donner_deep: plcl ', Don_cape%plcl(i,j)
          write (unit, '(a, e20.12)')  & 
                 'in donner_deep: plfc ', Don_cape%plfc(i,j)
          write (unit, '(a, e20.12)')  & 
              'in donner_deep: plzb ', Don_cape%plzb(i,j)

!---------------------------------------------------------------------
!    output the lag time value of cape (xcape_lag), the convective 
!    inhibition (coin), the time tendency of cape (dcape) and the lag
!    time column integrated water vapor (qint_lag).
!---------------------------------------------------------------------
          write (unit, '(a, e20.12)')  & 
               'in donner_deep: xcape ',   &
                              Don_cape%xcape_lag(i,j)
          write (unit, '(a, e20.12)')  & 
               'in donner_deep: coin ', Don_cape%coin(i,j)
          write (unit, '(a, e20.12)')  & 
              'in donner_deep: dcape ', Don_conv%dcape(i,j)
          write (unit, '(a, e20.12)')  & 
              'in donner_deep: qint ',  Don_cape%qint_lag(i,j)

!---------------------------------------------------------------------
!    output the total cloud fractional area (a1), the maximum allowed
!    value for a1 (amax), the maximum cloud fractional area based on the
!    moisture constraint (amos), the total precipitation from the col-
!    umn (total_precip), the mesoscale cloud fractional area (ampta1),
!    the displacement of a parcel from its initial location due to 
!    accrued upward motion at the current time (parcel_rise), and the
!    convective precipitation rate (cell_precip).
!---------------------------------------------------------------------
          write (unit, '(a, e20.12)')  & 
               'in donner_deep: a1   ', Don_conv%a1(i,j)
          write (unit, '(a, e20.12)')  & 
            'in donner_deep: amax ', Don_conv%amax(i,j)
          write (unit, '(a, e20.12)')  & 
            'in donner_deep: amos ', Don_conv%amos(i,j)
          write (unit, '(a, e20.12)')  & 
            'in donner_deep: tprea1 ', total_precip(i,j)
          write (unit, '(a, e20.12)')  & 
             'in donner_deep: ampta1 ', Don_conv%ampta1(i,j)
          write (unit, '(a, e20.12)')  & 
              'in donner_deep: omint', parcel_rise(i,j)
          write (unit, '(a, e20.12)')  & 
               'in donner_deep: rcoa1 ', Don_conv%cell_precip(i,j)

!---------------------------------------------------------------------
!    output various 3d fields between the specified highest index at
!    which diagnostics are to be output (kstart) and the nearest 
!    level to the surface (nlev_lsm), provided there has been some effect 
!    of deep convection at the level.
!---------------------------------------------------------------------
          do k=Col_diag%kstart,nlev_lsm
            if (temperature_forcing (i,j,k) == 0.0) cycle
            write (unit, '(a, i4)')'in donner_deep: k = ', k
            write (unit, '(a, e20.12)')  &
                 'in donner_deep: cemetf output to calling routine',  &
                     temperature_forcing(i,j,k)             
            write (unit, '(a, e20.12)')  &
                    'in donner_deep:TOTAL convective cemetf',  &
                     Don_conv%conv_temp_forcing(i,j,k)            
            write (unit, '(a, e20.12)')  &
                      'in donner_deep: ceefc ',     &
                               Don_conv%ceefc(i,j,k)             
            write (unit, '(a, e20.12)')  &
                      'in donner_deep: cecon ',  &
                                Don_conv%cecon(i,j,k)
            write (unit, '(a, e20.12)')  &
                     'in donner_deep: cemfc ',   &
                               Don_conv%cemfc(i,j,k)
            write (unit, '(a, e20.12)')  &
                 'in donner_deep: cememf output to calling routine',  &
                        moisture_forcing(i,j,k)               
            write (unit, '(a, e20.12)')  &
                     'in donner_deep: TOTAL convective cememf',  &
                        Don_conv%conv_moist_forcing (i,j,k)            
            write (unit, '(a, e20.12)')  &
                      'in donner_deep: cememf_mod',  &
                                Don_conv%cememf_mod(i,j,k)            
            write (unit, '(a, e20.12)')  &
                       'in donner_deep: cual  ',  &
                                  Don_conv%cual(i,j,k)              
            write (unit, '(a, e20.12)')  &
                       'in donner_deep: fre   ',   &
                                  Don_conv%fre(i,j,k)             
            write (unit, '(a, e20.12)')  &
                        'in donner_deep: elt   ',  &
                                     Don_conv%elt(i,j,k)            
            write (unit, '(a, e20.12)')  &
                         'in donner_deep: cmus  ',    &
                                     Don_conv%cmus(i,j,k)            
            write (unit, '(a, e20.12)')  &
                         'in donner_deep: ecds ',   &
                                    Don_conv%ecds(i,j,k)             
            write (unit, '(a, e20.12)')  &
                         'in donner_deep: eces  ', &
                                     Don_conv%eces(i,j,k)            
            write (unit, '(a, e20.12)')  &
                         'in donner_deep: emds  ',  &
                                     Don_conv%emds(i,j,k)             
            write (unit, '(a, e20.12)')  &
                          'in donner_deep: emes  ',  &
                                     Don_conv%emes(i,j,k)              
            write (unit, '(a, e20.12)')  &
                          'in donner_deep: qmes  ',  &
                                    Don_conv%mrmes(i,j,k)            
            write (unit, '(a, e20.12)')  &
                          'in donner_deep: wmps  ', &
                                      Don_conv%wmps(i,j,k)            
            write (unit, '(a, e20.12)')  &
                          'in donner_deep: wmms  ',  &
                                       Don_conv%wmms(i,j,k)            
            write (unit, '(a, e20.12)')  &
                           'in donner_deep: tmes  ',  &
                                        Don_conv%tmes(i,j,k)            
            write (unit, '(a, e20.12)')  &
                             'in donner_deep: dmeml ',   &
                                       Don_conv%dmeml(i,j,k)            
            write (unit, '(a, e20.12)')  &
                              'in donner_deep: uceml ',  &
                                       Don_conv%uceml(i,j,k)            
            write (unit, '(a, e20.12)')  &
                              'in donner_deep: detmfl ',  &
                                      Don_conv%detmfl(i,j,k)            
            write (unit, '(a, e20.12)')  &
                               'in donner_deep: umeml ',   &
                                      Don_conv%umeml(i,j,k)            

!---------------------------------------------------------------------
!    output various tracer-related fields for each tracer transported
!    by donner_deep_mod.
!---------------------------------------------------------------------
            do kcont=1,ntr     
              write (unit, '(a, e20.12)')  &
                              'in donner_deep: xgcm1 ',   &
                            tracers(i,j,k,kcont)                 
              write (unit, '(a, e20.12)')  &
                               'in donner_deep: qtren1 ',  &
                              Don_conv%qtren1(i,j,k,kcont)             
              write (unit, '(a, e20.12)')  &
                                'in donner_deep: qtmes1 ',  &
                              Don_conv%qtmes1(i,j,k,kcont)             
              write (unit, '(a, e20.12)')  &
                                   'in donner_deep: qtceme ',   &
                               Don_conv%qtceme(i,j,k,kcont)             
              write (unit, '(a, e20.12)')  &
                                  'in donner_deep: wtp1 ',   &
                                Don_conv%wtp1(i,j,k,kcont)            
            end do
          end do  ! (k loop)
        endif
      end do    ! (n loop)

!--------------------------------------------------------------------


end subroutine don_d_column_end_of_step_k
 



!#####################################################################

subroutine don_d_convert_profile_k     &
         (name_hi, name_lo, n_lo, n_hi, ncc, profile_hi, press_hi, ptop,&
          include_set_value, include_sbl, include_conservation_factor, &
          set_value, sbl, conservation_factor, press_lo, diag_unit,  & 
          debug_ijt, profile_lo, ermesg)

  USE MESSY_CONVECT_DONNER_UTIL,   ONLY: don_u_compare_integrals_k,   &
                                         don_u_set_column_integral_k, & 
                                         don_u_apply_integral_source_k
!----------------------------------------------------------------------
!    subroutine don_d_convert_profile_k takes an input profile 
!    (profile_hi) associated with a character string name_hi on the 
!    hi-res model grid (press_hi) containing ncc_ens levels and extending
!    to a pressure level ptop and maps it to variable profile_lo assoc-
!    iated with character string name_lo on the lo-res model grid defined
!    by press_lo.
!    additonally, if desired, the integral of the profile on the lo-res
!    grid multiplied by conservation_factor may be set to set_value by 
!    modifying the profile below cloud base, or a specified sub-cloud 
!    source (sbl) may be added to the lo-res profile.
!    if column diagnostics are desired (debug_ijt), they are output to
!    diag_unit.
!-----------------------------------------------------------------------
USE MESSY_CONVECT_DONNER_UTIL,   ONLY: don_u_map_hires_c_to_lores_c_k
implicit none

character(len=*),      intent(in)  :: name_hi, name_lo
integer,               intent(in)  :: n_lo, n_hi, ncc
real, dimension(n_hi), intent(in)  :: profile_hi, press_hi
real,                  intent(in)  :: ptop
logical,               intent(in)  :: include_set_value, include_sbl, &
                                      include_conservation_factor
real,                  intent(in)  :: set_value, sbl
real, dimension(n_lo), intent(in)  :: conservation_factor, press_lo
integer,               intent(in)  :: diag_unit
logical,               intent(in)  :: debug_ijt
real, dimension(n_lo), intent(out) :: profile_lo
character(len=*),      intent(out) :: ermesg

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       name_hi       character string associated with input profile
!       name_lo       character string associated with output profile
!       n_lo          number of levels on lo-res grid
!       n_hi          number of levels on hi_res grid
!       ncc           number of layers in input profile that are affected
!                     by presence of cloud; it may be called with 
!                     ncc_kou for each ensemble member 
!                     (from define_lo_res_model_profiles or 
!                     add_to_ensemble_sum_hires), or with ncc_ens
!                     (from define_ensemble_profiles).
!       profile_hi    vertical profile on hi-res model grid
!       press_hi      full pressure levels of hi-res model [ Pa ]
!       ptop          pressure one level above cloud top  [ Pa ]
!       include_set_value
!                     it is desired to force the column integral to a
!                     specified value on the lo-res grid ?
!       include_sbl   it is desired to add a specified value to the
!                     profile in the layers below cloud base ?
!       include_conservation_factor
!                     the integrand which is to be set to set_value 
!                     includes a non-unity factor which multiplies the 
!                     profile ?
!       set_value     value desired for the integral of the
!                     output profile times conservation_factor      
!       sbl           value to be added to the profile in all layers
!                     below cloud base
!       conservation_factor
!                     the column integral of the product of the profile 
!                     and conservation_factor arrays is required to equal
!                     set_value
!       press_lo      interface pressure levels of lo-res model [ Pa ]
!       diag_unit     unit number for column diagnostics file
!       debug_ijt     column diagnostics are desired for this column ?
!
!   intent(out) variables:
!
!       profile_lo    vertical profile on lo-res model grid
!       ermesg        error message produced by any kernel routines
!                     called by this subroutine 
!
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:
  
      real, dimension(n_lo) :: out       ! intermediate lo-res profile 
                                         ! after either setting column
                                         ! integral or adding boundary 
                                         ! layer source
      real, dimension(n_lo) :: conservation_factor_used                
                                         ! conservation_factor array 
                                         ! used in calculation; is array
                                         ! of 1.0 when 
                                         ! include_conservation_factor
                                         ! is .false.
      real                  :: intgl_hi  ! column integral of profile_hi 
      real                  :: intgl_lo  ! column integral of profile_lo
      integer               :: k         ! do-loop index

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  '

!---------------------------------------------------------------------
!    if column diagnostics are desired, output a diagnostic message 
!    indicating the variable that is being processed.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a)')  &
           'in mulsub: map_hi_res_col_to_lo_res_col: ' // trim(name_hi)
      endif

!----------------------------------------------------------------------
!   call don_u_map_hires_c_to_lores_c_k to map the 
!   profile from the hi-res model grid to the lo-res model grid.
!----------------------------------------------------------------------
      call don_u_map_hires_c_to_lores_c_k     &
          (n_lo, ncc+1, profile_hi(1:ncc+1), press_hi(1:ncc+1),  &
           ptop, press_lo, profile_lo, intgl_hi, intgl_lo, ermesg)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') then
        return
      endif

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrals of the
!    profiles on both the hi- and lo-res grids.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
             'in mulsub: rintsum(' // trim(name_lo) // ' ) =',  &
                                              intgl_hi, intgl_lo

!---------------------------------------------------------------------
!    call don_u_compare_integrals_k to assess if the integrals
!    from the two grids are "equal", as they should be.
!---------------------------------------------------------------------
        call don_u_compare_integrals_k    &
                         (intgl_hi, intgl_lo, diag_unit, ermesg)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (trim(ermesg) /= ' ') then
          return
        endif
      endif

!----------------------------------------------------------------------
!    if it is desired to set the integral value for the profile (i.e.,
!    set_value does not equal dummy_set_value), execute the following
!    code.
!----------------------------------------------------------------------
      if (include_set_value) then
        if (include_conservation_factor) then
          conservation_factor_used(:) = conservation_factor(:)
        else
          conservation_factor_used(:) = 1.0                   
        endif
        
!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrands at each
!    level on the lo-res grid.
!---------------------------------------------------------------------
        if (debug_ijt) then
          do k=1,n_lo                   
            if (profile_lo(k) /= 0.0) then
              write (diag_unit, '(a, i4, e20.12)') &
                 'in set_col_integral: k,phr,phr+= ', k, profile_lo(k)* &
                              conservation_factor_used(k)
            endif
          end do
        endif

!-----------------------------------------------------------------------
!    call don_u_set_column_integral_k to adjust the output
!    profile below cloud base so that the desired integral value is
!    obtained.
!-----------------------------------------------------------------------
        call don_u_set_column_integral_k    &
               (n_lo, profile_lo*conservation_factor_used, press_hi(1), &
                press_lo(1), set_value, press_lo, intgl_hi,     &
                intgl_lo, out, ermesg)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (trim(ermesg) /= ' ') then
          return
        endif

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrals and 
!    profiles, both before and after the adjustment to the desired value.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, e20.12)')  &
                           'in set_col_integral: column(in)= ',intgl_hi
          write (diag_unit, '(a, e20.12)')  &
                          'in set_col_integral: column(out)= ',intgl_lo 
          do k=1,n_lo                 
            if (profile_lo(k)*conservation_factor_used(k) /= out(k)) then
              write (diag_unit, '(a, i4, 2e20.12)') &
               'in set_col_integral: k,qtr(in), qtr(out)= ', k,  &
                      profile_lo(k)*conservation_factor_used(k), out(k)
            endif
          end do
        endif

!---------------------------------------------------------------------
!    define the adjusted output profile by removing conservation_factor.
!---------------------------------------------------------------------
        profile_lo(:) = out(:)/conservation_factor_used(:)
      endif !(set_value /= dummy_set_value)

!----------------------------------------------------------------------
!    if a boundary layer source is to be added to the profile, execute
!    the following code.
!----------------------------------------------------------------------
      if (include_sbl .and. sbl /= 0.0) then

!----------------------------------------------------------------------
!    call don_u_apply_integral_source_k to apply the imposed 
!    subcloud source (sbl) to the input profile profile_out, resulting 
!    in the output profile out.  also returned are the column integrals
!    of the input profile (intgl_in) and the integral of the output
!    profile (intgl_out).
!    NOTE: in the original code, the subcloud source was not applied in 
!    the non-entropy case anywhere, and in the entropy case only to the 
!    model layer containing cloud base. 
!    I have MODIFIED THE CODE so that the value is APPLIED FROM SFC TO
!    TOP OF SPECIFIED REGION (CLOUD BASE) IS THIS CORRECT AND WHAT WAS
!    INTENDED ?
!----------------------------------------------------------------------
        call don_u_apply_integral_source_k     &
             (n_lo, profile_lo, press_hi(1), press_lo(1), sbl,  &
              press_lo, intgl_hi, intgl_lo, out, ermesg)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (trim(ermesg) /= ' ') then
          return
        endif

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrals and 
!    profiles, both before and after adding the boundary layer source.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 3e20.12)')  &
             'after apply_subcloud: column(in)= ',   &
                             intgl_hi, press_lo(1), press_hi(1)
          write (diag_unit, '(a, e20.12)')  &
                           'after apply_subcloud: column(out)= ',  &
                                                       intgl_lo 
          do k=1,n_lo                  
            if (profile_lo(k) /= out(k)) then
              write (diag_unit, '(a, i4, 2e20.12)') &
               'in set_col_integral: k,qtr(in), qtr(out)= ', k,  &
                                       profile_lo(k), out(k)
            endif
          end do
        endif

!----------------------------------------------------------------------
!    define the output profile on the lo-res model grid to be returned to
!    the calling routine.
!----------------------------------------------------------------------
        profile_lo(:) = out(:)
      endif  !(sbl /= 0.0)

!---------------------------------------------------------------------


end subroutine don_d_convert_profile_k



!#####################################################################

subroutine don_d_def_ensemble_profs_k    &
         (nlev_lsm, nlev_hires, ncc_ens, diag_unit, debug_ijt, ptt,  &
          cld_press, alp, detmfh, ucemh, cuql, cuqli, phalf_c,  &
          ensmbl_cloud_area, cuql_v, cuq, detmfl, uceml, ermesg)


!---------------------------------------------------------------------
!    subroutine don_d_def_ensemble_profs_k defines vertical 
!    profiles of cloud area, cloud ice, cloud liquid, vertical mass flux
!    and detrained vertical mass flux produced by the entire cumulus 
!    ensemble on the lo-res grid. 
!---------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------
integer,                        intent(in)   :: nlev_lsm, nlev_hires,&
                                                ncc_ens, diag_unit
logical,                        intent(in)   :: debug_ijt
real,                           intent(in)   :: ptt
real,    dimension(nlev_hires), intent(in)   :: cld_press, alp,  &
                                                detmfh, ucemh, cuql,  &
                                                cuqli
real,    dimension(nlev_lsm+1), intent(in)   :: phalf_c
real,    dimension(nlev_lsm),   intent(out)  :: ensmbl_cloud_area,  &
                                                cuql_v, cuq, detmfl, &
                                                uceml     
character(len=*),               intent(out)  :: ermesg

!--------------------------------------------------------------------
!   intent(in) variables:
!
!        ptt        pressure one cloud model level above the ensemble
!                   cloud top [ Pa ]
!        nlev_lsm   number of levels in the low resolution grid
!        nlev_hires       number of levels in the high resolution grid
!        ncc_ens    cloud top index for the ensemble on the hi-res grid
!        cld_press  pressures at hi-res model levels [ Pa ]
!        alp        cloud area profile on hi-res model grid
!                   [ fraction ]
!        detmfh     detrained mass flux (layer above index level)
!                   (on cloud-model grid) (index 1 at cloud base)
!                   [ kg (air) / (sec m**2) ]
!        ucemh      upward mass flux on cloud model grid            
!                   [ kg (air) / (sec m**2) ]
!        cuql       ice water profile on cloud model grid; on input is
!                   normalized by total grid box area, on output is
!                   normalized by ensemble cloud area.
!        cuqli      liquid water profile on cloud model grid; on input 
!                   is normalized by total grid box area, on output is
!                   normalized by ensemble cloud area.
!        phalf_c    pressure at lo-res model half levels [ Pa ]
!        diag_unit  unit for column diagnostics output
!        debug_ijt  are column diagnostics desired in this column ?
!
!   intent(out) variables:
!
!        ensmbl_cloud_area  
!                   total cloud area profile over all ensemble members
!                   on large_scale model grid [ fraction ]
!        cuql_v     liquid water profile on large-scale model grid, 
!                   normalized by ensemble cloud area.
!        cuq        ice water profile on large-scale model grid, 
!                   normalized by ensemble cloud area.
!        detmfl     detrained mass flux on large-scale model grid
!                   normalized by ensemble cloud area
!                   index 1 near large-scale model base
!                   [ kg (air) / (sec m**2) ]
!        uceml      upward mass flux on large_scale model grid       
!                   [ kg (air) / (sec m**2) ]
!        ermesg     error message produced by any kernel routines
!                   called by this subroutine 
!
!--------------------------------------------------------------------

      real, dimension (nlev_lsm) :: conv_fact

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' '

      conv_fact = 0.0

!----------------------------------------------------------------------
!    call don_d_convert_profile_k to map the ensemble-total cloud 
!    area profile from the cloud model grid (alp) to the large-scale 
!    model grid (ensmbl_cloud_area).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('alp', 'cual', nlev_lsm, nlev_hires, ncc_ens, alp, cld_press, &
          ptt, .false., .false., .false.,  0.0, 0.0, conv_fact, &
          phalf_c, diag_unit, debug_ijt, ensmbl_cloud_area, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    call convert_profile to map the ensemble-total condensed ice 
!    profile from the cloud model grid (cuql) to the large-scale model 
!    grid (cuq).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('cuql', 'cuq', nlev_lsm, nlev_hires, ncc_ens, cuql, cld_press,&
          ptt, .false., .false., .false.,  0.0, 0.0, conv_fact, &
          phalf_c, diag_unit, debug_ijt, cuq, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    call convert_profile to map the ensemble-total condensed liquid
!    profile from the cloud model grid (cuql) to the large-scale model 
!    grid (cuq).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('cuqli', 'cuql_v', nlev_lsm, nlev_hires, ncc_ens, cuqli, &
          cld_press, ptt, .false., .false., .false.,  0.0, 0.0,  &
          conv_fact, phalf_c, diag_unit, debug_ijt, cuql_v, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    call convert_profile to map the ensemble-total upward mass flux
!    profile from the cloud model grid (ucemh) to the large-scale model 
!    grid (uceml).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('ucemh', 'uceml', nlev_lsm, nlev_hires, ncc_ens, ucemh, &
          cld_press, ptt, .false., .false., .false.,  0.0, 0.0,   &
          conv_fact, phalf_c, diag_unit, debug_ijt, uceml, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    call convert_profile to map the ensemble-total detrained mass flux
!    profile from the cloud model grid (detmfh) to the large-scale model 
!    grid (detmfl).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('detmfh', 'detmfl', nlev_lsm, nlev_hires, ncc_ens, detmfh, &
          cld_press, ptt, .false., .false., .false.,  0.0, 0.0,  &
          conv_fact, phalf_c, diag_unit, debug_ijt, detmfl, ermesg)
 
!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------

end subroutine don_d_def_ensemble_profs_k
    

!#####################################################################

subroutine don_d_def_lores_model_profs_k               &
         (nlev_lsm, nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt,  &
          Param, pb, ptt, sfc_vapor_flux_c, sfc_sh_flux_c,   &
          sfc_tracer_flux_c, pfull_c, phalf_c, cld_press, dpf, dfr,  &
          cld_evap, qlw, emfhr, efchr, etfhr, cell_freeze, evap_rate, &
          h1, h1_2, q1, qtr, ermesg)
 
!---------------------------------------------------------------------
!    subroutine don_d_def_lores_model_profs_k maps vertical
!    profiles of various fields from the cloud-model grid to the large-
!    scale model grid. also, if desired, the sub-cloud base model levels
!    of the lo-res profiles may be modified so that the column integral 
!    equals a prescribed value (set_value), and / or a given value may
!    be assigned to the sub-cloud base levels.
!    this routine is called for each ensemble member individually, so 
!    that the input and output profiles are weighted by the cloud area 
!    of the ensemble member.
!---------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_param_type

implicit none

!---------------------------------------------------------------------
integer,                            intent(in)    :: nlev_lsm,   &
                                                     nlev_hires,  &
                                                     ntr, ncc_kou, &
                                                     diag_unit
logical,                            intent(in)    :: debug_ijt
type(donner_param_type),            intent(in)    :: Param
real,                               intent(in)    :: pb, ptt,  &
                                                     sfc_vapor_flux_c, &
                                                     sfc_sh_flux_c
real,    dimension(ntr),            intent(in)    :: sfc_tracer_flux_c
real,    dimension(nlev_lsm),       intent(in)    :: pfull_c
real,    dimension(nlev_lsm+1),     intent(in)    :: phalf_c
real,    dimension(nlev_hires),     intent(in)    :: cld_press, dpf, &
                                                     dfr, cld_evap, qlw,&
                                                     emfhr, efchr
real,    dimension(nlev_hires,ntr), intent(in)    :: etfhr
real,    dimension(nlev_lsm),       intent(out)   :: cell_freeze, &
                                                     evap_rate, h1,  &
                                                     h1_2, q1
real,    dimension(nlev_lsm,ntr),   intent(out)   :: qtr
character(len=*),                   intent(out)   :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       nlev_lsm         number of levels on lo-res grid
!       nlev_hires             number of levels on hi_res grid
!       ntr              number of tracers being transported by the 
!                        donner deep convection parameterization
!       ncc_kou          number of layers in hi-res profile that are 
!                        affected by the presence of cloud
!       pb               pressure at cloud base [ Pa ]
!       ptt              pressure one model level above cloud top [ Pa ]
!       sfc_vapor_flux_c flux of water vapor from the surface into the 
!                        sub cloud layer [ kg(h2o) / (m**2 sec) ]
!       sfc_sh_flux_c    flux of sensible heat from the surface into the
!                        sub cloud layer [ W / m**2, or kg / m**3 ] 
!       sfc_tracer_flux_c  flux of tracer from the surface into the sub-
!                        cloud layer  [ kg(tracer) / (m**2 sec) ]
!       pfull_c           pressure at lo-res model full levels [ Pa ]
!       phalf_c         pressure at lo-res model half levels [ Pa ]
!       cld_press        pressure at hi-res model full levels [ Pa ]
!       dpf              condensation rate profile
!                        on hi-res grid  for the current ensemble member
!                        [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!       dfr              profile of                     moisture tendency
!                        due to freezing in the convective updraft 
!                        on hi-res grid  for the current ensemble member
!                        [ (      g(h2o) ) / ( kg(dry air) sec ) ] 
!!!!!!!!======>>>>>>>    NOTE UNITS OF g(h2o). (Verify and change.)
!       cld_evap                             profile of the potential
!                        cloud water evaporation for the curent ensemble
!                        member on th ehi-res grid. this amount of water!
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ (      kg(h2o) ) / ( kg(dry air) sec ) ] 
!       qlw              profile of cloud water for the current ensemble
!                        member [ kg(h2o) / kg(air) ]
!       emfhr            vertical moisture flux convergence profile on 
!                        hi-res grid for the current ensemble member 
!                        [ kg(h2o) / ( kg(dry air) sec ) ]
!       efchr            vertical entropy flux convergence profile on
!                        hi-res grid for the current ensemble member 
!                        [ deg K / sec ]                        
!       etfhr            vertical tracer flux convergence profile on
!                        hi-res grid for the current ensemble member 
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!       diag_unit        unit number of column diagnostics output file
!       debug_ijt        logical indicating whether diagnostics are 
!                        requested for this column 
!
!   intent(out) variables:
!
!       cell_freeze      profile of cloud-area-weighted moisture tendency
!                        due to freezing in the convective updraft 
!                        on lo-res grid for the current ensemble member
!                        [ (      g(h2o) ) / ( kg(dry air) sec ) ] 
!!!!!!!!======>>>>>>>    NOTE UNITS OF g(h2o). (Verify and change.)
!       evap_rate        cloud-area-weighted profile of the potential
!                        cloud water evaporation for the current ensemble
!                        member on the lo-res grid. this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ (      kg(h2o) ) / ( kg(dry air) sec ) ] 
!       h1                                   condensation rate profile
!                        on lo-res grid for the current ensemble member
!                        [ (      kg(h2o) ) / ( kg( dry air) sec ) ] 
!       h1_2             vertical entropy flux convergence profile on
!                        lo-res grid for the current ensemble member 
!                        [ deg K / sec ]                        
!       q1               vertical moisture flux convergence profile on 
!                        lo-res grid for the current ensemble member 
!                        [ kg(h2o) / ( kg(dry air) sec ) ]
!       qtr              vertical tracer flux convergence profile on
!                        lo-res grid for the current ensemble member 
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!       ermesg           error message produced by any kernel routines
!                        called by this subroutine 
!
!---------------------------------------------------------------------


!-----------------------------------------------------------------------
!   local variables:

      real, dimension (nlev_lsm)  ::   pi     ! inverse exner function
                                              ! used for setting column
                                              ! integral value (conserv-
                                              ! ation of theta)
      real, dimension (nlev_lsm)  ::   condensate ! liquid water profile on 
      real, dimension (nlev_lsm)  ::   conv_fact  
      integer                 ::   kcont, k   ! do-loop indices
      real                    ::   sbl        ! value to be used for
                                              ! the profile at levels
                                              ! below cloud base
      real                    ::   set_value  ! desired column integral
                                              ! value 

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' '

      conv_fact = 0.0

!--------------------------------------------------------------------
!    call don_d_convert_profile_k to map the moisture tendency due
!    to freezing from the cloud model grid (dfr) to the large-scale model
!    grid (cell_freeze.
!--------------------------------------------------------------------
      call don_d_convert_profile_k   &
           ('DFR', 'frea', nlev_lsm, nlev_hires, ncc_kou,   &
            dfr(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt,   &
            .false., .false., .false., 0.0, 0.0, conv_fact,   &
            phalf_c, diag_unit, debug_ijt, cell_freeze, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    map the cloud condensate (qlw) from the cloud model to the large-
!    scale model (condensate). this field is only used for diagnostic 
!    purposes.
!----------------------------------------------------------------------
      if (debug_ijt) then
        call don_d_convert_profile_k    &
             ('QLW', 'evap', nlev_lsm, nlev_hires, ncc_kou,   &
              qlw(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt, & 
              .false., .false., .false., 0.0, 0.0, conv_fact,&
              phalf_c, diag_unit, debug_ijt, condensate, ermesg)
      endif
      
!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    map the rate at which condensate which has not precipitated out
!    must evaporate from the cloud model grid (cld_evap) to the lo-res
!    model grid (evap_rate).
!---------------------------------------------------------------------
      call don_d_convert_profile_k    &
           ('QLW', 'evap_rate', nlev_lsm, nlev_hires, ncc_kou,   &
            cld_evap(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt, &
            .false., .false., .false., 0.0, 0.0, conv_fact,&
            phalf_c,   diag_unit,  debug_ijt, evap_rate, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    if in diagnostics column, output profiles of the cloud evaporation
!    rate (cld_evap) and evaporation(qlw) on the hi-res model grid. 
!    cld_evap will be the actual evaporation rate if there turns out to 
!    be  no mesoscale circulation in the column, while qlw is the profile
!    of liquid water produced by the given ensemble member.
!----------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a)') 'in mulsub: P & CLD_EVAP'
        do k=1,ncc_kou-1
            write (diag_unit, '(a, i4, 2e20.12)')  &
                 'in mulsub: k, P & QLW', k, cld_press(k), cld_evap (k)
        end do
        write (diag_unit, '(a)') 'in mulsub: P & QLW'
        do k=1,ncc_kou-1
            write (diag_unit, '(a, i4, 2e20.12)')  &
                 'in mulsub: k, P & QLW', k, cld_press(k), qlw      (k)
        end do
      endif
      
!----------------------------------------------------------------------
!    map the condensation rate from the cloud model (-dpf) to the 
!    large-scale model (h1). h1 is a term appropriate for use in the
!    water vapor equation; i.e., condensation is a loss term.
!----------------------------------------------------------------------
      call don_d_convert_profile_k                 &
           ('RLHR', 'h1', nlev_lsm, nlev_hires, ncc_kou,    &
            -dpf(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt,  &
            .false., .false., .false., 0.0, 0.0, conv_fact,&
            phalf_c, diag_unit, debug_ijt, h1, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    determine the vertical flux convergence of each tracer.
!---------------------------------------------------------------------
      do kcont=1,ntr  

!----------------------------------------------------------------------
!    calculate the imposed subcloud tracer-flux convergence (sbl) in 
!    units of kg(tracer) per kg(dry air) per sec. define set_value so
!    that the column integrated tracer flux convergence will be set to
!    zero.
!----------------------------------------------------------------------
        sbl = (sfc_tracer_flux_c(kcont)*Param%grav)/(phalf_c(1) - pb)
        set_value = 0.0

!----------------------------------------------------------------------
!    call convert_profile to map the vertical tracer flux convergence 
!    from the cloud model (etfhr) to the large-scale model grid (qtr). 
!    force the column integral of the flux convergence to be 0.0; then 
!    add the imposed sub-cloud convergence.
!----------------------------------------------------------------------
        call don_d_convert_profile_k     &
             ('qtrv', 'qtr', nlev_lsm, nlev_hires, ncc_kou,    &
              etfhr(1:ncc_kou+1,kcont), cld_press(1:ncc_kou+1), ptt,   & 
              .true., .true., .false., set_value, sbl, conv_fact, &
              phalf_c, diag_unit, debug_ijt, qtr(:,kcont),ermesg)
      end do

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    define the subcloud moisture flux convergence (sbl) in units of
!    kg(h2o) per kg(air) per sec. sfc_vapor_flux_c is the imposed bound-
!    ary layer moisture source in units of kg(h2o) per m**2 per sec. 
!----------------------------------------------------------------------
      sbl = (sfc_vapor_flux_c*Param%grav)/(phalf_c(1) - pb)
      set_value = 0.0

!----------------------------------------------------------------------
!    call don_d_convert_profile_k to map the vertical moisture 
!    flux convergence from the cloud model (emfhr) to the lo-res model 
!    grid (q1). force the column integral of the flux convergence to be 
!    0.0; then add the imposed sub-cloud convergence.
!----------------------------------------------------------------------
      call don_d_convert_profile_k       &
           ('EMFHR', 'q1', nlev_lsm, nlev_hires, ncc_kou,    &
            emfhr(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt,   &
            .true., .true., .false., set_value, sbl, conv_fact, &
            phalf_c, diag_unit,  debug_ijt, q1, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------
!    calculate the subcloud entropy flux convergence (sbl) in units of
!    deg K per sec. sfc_sh_flux_c is the imposed boundary layer sensible
!    heat source in units of watts per square meter. define the inverse
!    exner function so that an integral constraint on theta may be
!    applied.
!----------------------------------------------------------------------
      sbl = Param%grav*sfc_sh_flux_c/((phalf_c(1) - pb)*Param%cp_air)
      set_value = 0.0
      do k=1,nlev_lsm               
        pi(k) = (1.0e05/pfull_c(k))**Param%kappa
      end do

!----------------------------------------------------------------------
!    map the vertical entropy flux convergence from the cloud model 
!    (efchr) to the large-scale model grid (h1_2). force the column 
!    integral of the flux convergence times inverse exner function (i.e.,
!    theta) to be 0.0; then add the imposed sub-cloud convergence.
!----------------------------------------------------------------------
      call don_d_convert_profile_k       &
           ('EFCHR', 'h1_2', nlev_lsm, nlev_hires, ncc_kou,   &
            efchr(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt, &
            .true., .true., .true., set_value, sbl, pi, &
            phalf_c, diag_unit,  debug_ijt, h1_2, ermesg)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!---------------------------------------------------------------------
!    if in diagnostics column, output the profile of the entropy flux
!    convergence on the lo-res grid, at levels where it is non-zero.
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (debug_ijt) then
          if (h1_2(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)')  &
                                 'in mulsub: JK,H1= ', k, h1_2(k)
          endif
        endif
      end do

!----------------------------------------------------------------------



end subroutine don_d_def_lores_model_profs_k


!#####################################################################

subroutine don_d_add_to_ensmbl_sum_hires_k     &
         (nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt, area_ratio,  &
          cfracice, rcl, flux, emfhr, dpf, qlw, etfhr, cuql, cuqli, &
          ucemh, alp, rlsm, emsm, detmfh, etsm, ermesg)

!-----------------------------------------------------------------------
!    subroutine don_d_add_to_ensmbl_sum_hires_k adds the contrib-
!    utions from this ensemble member to various profiles on the hi-res
!    grid that are being summed over the ensemble.
!-----------------------------------------------------------------------

implicit none

!-----------------------------------------------------------------------
integer,                         intent(in   )  :: nlev_hires, ntr, ncc_kou,&
                                                   diag_unit
logical,                         intent(in   )  :: debug_ijt
real,                            intent(in   )  :: area_ratio
real, dimension(nlev_hires),     intent(in   )  :: cfracice, rcl, flux, &
                                                   emfhr, dpf, qlw     
real, dimension(nlev_hires,ntr), intent(in   )  :: etfhr
real, dimension(nlev_hires),     intent(inout)  :: cuql, cuqli, ucemh, &
                                                   alp, rlsm, emsm,  &
                                                   detmfh
real, dimension(nlev_hires,ntr), intent(inout)  :: etsm         
character(len=*),                intent(  out)  :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      nlev_hires            number of levels on hi_res grid
!      ntr             number of tracers being transported by the 
!                      donner deep convection parameterization
!      ncc_kou         number of layers in hi-res profile that are 
!                      affected by the presence of cloud
!      area_ratio      ratio of cloud base area of this ensemble member
!                      to that of ensemble member # 1. (ensemble member
!                      # 1 assumed to have largest cloud base area)
!                      [ dimensionless ]
!      cfracice        fraction of condensate that is ice [ fraction ]
!      rcl             profile of cloud radius for this ensemble member
!                      [ m ]
!      flux            upward mass flux profile in cloud for this
!                      ensemble member [ kg (air) / sec ]
!      emfhr           vertical moisture flux convergence for this
!                      ensemble member [ kg (h2o) / ( kg(dry air) sec ) ]
!      dpf             condensation rate profile on hi-res grid for the 
!                      current ensemble member
!                      [ kg(h2o) / ( kg( dry air) sec ) ] 
!      qlw             profile of cloud water for the current ensemble
!                      member [ kg(h2o) / kg(air) ]
!      etfhr           vertical tracer flux convergence profile on
!                      hi-res grid for the current ensemble member 
!                      [ kg(tracer) / ( kg(dry air) sec ) ]
!      debug_ijt       logical indicating whether diagnostics are 
!                      requested for this column 
!      diag_unit       unit number of column diagnostics output file
!
!   intent(inout) variables:
!
!      cuql            vertical profile on the hi-res grid of condensed 
!                      ice, summed over ensemble members # 1 to the cur-
!                      rent, each member's contribution being weighted by
!                      its cloud area at level k relative to the cloud 
!                      base area of ensemble member #1
!                      [ kg(h2o) / kg (dry air) ]
!      cuqli           vertical profile on the hi-res grid of condensed 
!                      liquid, summed over ensemble members # 1 to the 
!                      current, each member's contribution being weighted
!                      by its cloud area at level k relative to the cloud
!                      base area of ensemble member #1
!                      [ kg(h2o) / kg (dry air) ]
!      ucemh           vertical profile on the hi-res grid of cell upward
!                      mass flux, summed over ensemble members # 1 to the
!                      current, each member's contribution being weighted
!                      by its cloud area at level k relative to the cloud
!                      base area of ensemble member #1
!                      [ kg (dry air) / ( m**2 sec ) ]
!      alp             vertical profile on the hi-res grid of cloud area
!                      summed over ensemble members # 1 to the current, 
!                      each member's contribution being weighted
!                      by its cloud area at level k relative to the cloud
!                      base area of ensemble member #1
!                      as a result, the cloud area profile is expressed
!                      relative to the cloud base area of ensemble member
!                      # 1. [ dimensionless ]
!      rlsm            vertical profile on the hi-res grid of conden-
!                      sation rate, summed over ensemble members # 1 to
!                      the current, each member's contribution being 
!                      weighted by its cloud area at level k relative to 
!                      the cloud base area of ensemble member #1
!                      [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!      emsm            vertical profile on the hi-res grid of vertical
!                      moisture flux convergence, summed over ensemble 
!                      members # 1 to the current, each member's contrib-
!                      ution being weighted by its cloud area at level k
!                      relative to the cloud base area of ensemble 
!                      member #1  [ kg (h2o) / ( kg(dry air) sec ) ]
!      detmfh          vertical profile on the hi-res grid of detrained
!                      mass flux in the layer above indexed level, summed
!                      over ensemble members # 1 to the current, each 
!                      member's contribution being weighted by its cloud 
!                      area at level k relative to the cloud base area of
!                      ensemble member #1 [ kg (dry air) / ( m**2 sec ) ]
!      etsm            vertical profile on the hi-res grid of vertical
!                      tracer flux convergence, summed over ensemble 
!                      members # 1 to the current, each member's contrib-
!                      ution being weighted by its cloud area at level k
!                      relative to the cloud base area of ensemble 
!                      member #1  [ kg (tracer) / ( kg(dry air) sec ) ]
!
!   intent(out) variables:
!
!      ermesg          error message produced by this subroutine or any
!                      kernel routines called by this subroutine 
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real       :: wt_factor   ! cloud area at level k for current 
                                ! ensemble member, relative to cloud
                                ! base area of ensemble member # 1
     integer     :: k, ktr      ! do-loop indices     

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  '

!--------------------------------------------------------------------
!    add the contributions from this ensemble member to the arrays 
!    accumulating ensemble sums on the cloud model levels.
!--------------------------------------------------------------------
      do k=1,ncc_kou 

!----------------------------------------------------------------------
!    define the factor needed to normalize each ensemble member's con-
!    tribution by the cloud base area of ensemble member #1. wt_factor
!    is the cloud area at level k for ensemble member kou, relative to
!    the cloud area at cloud base (k=1) for ensemble member #1.
!-----------------------------------------------------------------------
        wt_factor = area_ratio*(rcl(k)/rcl(1))**2
        
!----------------------------------------------------------------------
!    add this ensemble member's appropriately weighted contribution to
!    the ensemble-total cloud area (alp), condensed ice (cuql), condensed
!    liquid (cuqli), cell upward mass flux (ucemh), cell detrained mass 
!    flux (detmfh), condensation rate (rlsm), vertical moisture flux 
!    convergence (emsm) and vertical tracer flux convergence (etsm). the
!    weighting factor area_ratio*(rcl(k)/rcl(1))**2 allows the contrib-
!    utions from each member to be added by normalizing each member's 
!    contribution by the cloud base area of ensemble member #1.
!    NOTE: several of the arrays already have some of the normalizing
!    factors already included and so here need only to be multiplied by 
!    a portion of wt_factor.
!----------------------------------------------------------------------
        alp(k)   = alp(k)   + wt_factor                      
        cuql(k)  = cuql(k)  + wt_factor*(cfracice(k)*qlw(k))
        cuqli(k) = cuqli(k) + wt_factor*((1.0 - cfracice(k))*qlw(k))
        ucemh(k) = ucemh(k) + area_ratio*flux(k)/(rcl(1)**2)
        if (k < ncc_kou) then
          if (flux(k+1) < flux(k)) then
            detmfh(k) = detmfh(k) + area_ratio*   &
                        ((flux(k)-flux(k+1))/(rcl(1)**2))
          endif
        else
          if (flux(k) /= 0.)  then
            detmfh(k) = detmfh(k) - area_ratio*flux(k)/(rcl(1)**2)
          endif
        endif
        rlsm(k)   = rlsm(k)   - area_ratio*dpf (k) 
        emsm(k)   = emsm(k)   + area_ratio*emfhr(k)

!----------------------------------------------------------------------
!    if in a diagnostics column, output the total cell upward mass flux 
!    (ucemh) and the cloud area at each level ( alp).
!----------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                    'in mulsub: k,ucemh, alp= ',k,ucemh(k), alp(k)
        endif
      end do

!----------------------------------------------------------------------
!    add this ensemble member's appropriately weighted contribution to
!    the vertical tracer flux convergence (etsm). the weighting factor 
!    area_ratio allows the contributions from each member to be added by
!    normalizing each member's contribution by the cloud base area of 
!    ensemble member #1.
!----------------------------------------------------------------------
      do ktr=1,ntr
        do k=1,ncc_kou 
          etsm(k,ktr) = etsm(k,ktr) + area_ratio*etfhr(k,ktr)
        end do
      end do
      
!---------------------------------------------------------------------


end subroutine don_d_add_to_ensmbl_sum_hires_k 



!#####################################################################

subroutine don_d_add_to_ensmbl_sum_lores_k      &
         (nlev_lsm, ntr, diag_unit, debug_ijt, lmeso, Param,   &
          area_ratio, dint, cell_freeze, cell_melt, temp_c, h1_2, ecd, &
          ece, evap_rate, q1, h1, pfull_c, meso_melt, meso_freeze, &
          phalf_c, qtr, ensmbl_melt, ensmbl_freeze, enctf, encmf, enev,&
          disg, disb, disc, ecds, eces, disd, qtren, ermesg)

!-----------------------------------------------------------------------
!    subroutine don_d_add_to_ensmbl_sum_lores_k adds the contrib-
!    utions from this ensemble member to various profiles on the lo-res
!    grid that are being summed over the ensemble.
!-----------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_param_type

implicit none

!-----------------------------------------------------------------------
integer,                       intent(in   ) :: nlev_lsm, ntr, diag_unit
logical,                       intent(in   ) :: debug_ijt, lmeso
type(donner_param_type),       intent(in)    :: Param
!real,                      intent(in   ) :: seconds_per_day, hls, hlv, &
real,                          intent(in   ) :: area_ratio, dint     
real, dimension(nlev_lsm),     intent(in   ) :: cell_freeze, cell_melt, &
                                                temp_c, h1_2, ecd, ece, &
                                                evap_rate, q1, h1, &
                                                pfull_c, meso_melt,   &
                                                meso_freeze
real, dimension(nlev_lsm+1),   intent(in   ) :: phalf_c  
real, dimension(nlev_lsm,ntr), intent(in   ) :: qtr
real, dimension(nlev_lsm),     intent(inout) :: ensmbl_melt,   &
                                                ensmbl_freeze, enctf, &
                                                encmf, enev, disg, &
                                                disb, disc, ecds, eces, &
                                                disd
real, dimension(nlev_lsm,ntr), intent(inout) :: qtren 
character(len=*),              intent(  out) :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     
!       nlev_lsm         number of levels on lo-res grid
!       ntr              number of tracers being transported by the 
!                        donner deep convection parameterization
!       area_ratio       ratio of cloud base area of this ensemble 
!                        member to that of ensemble member # 1. 
!                        (ensemble member # 1 assumed to have largest 
!                        cloud base area) [ dimensionless ]
!       dint             column sum of moisture tendency due to freezing
!                        in the convective updraft on hi-res grid for the
!                        current ensemble member
!!!!  CHECK ON THESE UNITS !!!!!
!                        [ (      g(h2o) ) / ( kg(dry air) sec ) ] 
!       cell_freeze      profile of cloud-area-weighted moisture tendency
!                        due to freezing in the convective updraft 
!                        on lo-res grid for the current ensemble member
!                        [ (      g(h2o) ) / ( kg(dry air) sec ) ] 
!!!!!!!!======>>>>>>>    NOTE UNITS OF g(h2o). (Verify and change.)
!       cell_melt        in-cloud melting of condensate associated with
!                        convective cells. made up of two parts, 1) that
!                        due to the freezing of liquid carried upwards
!                        in the cell updraft, 2) that due to the melting 
!                        of condensed ice that precipitates out. if meso
!                        circulation is present, this component is zero;
!                        melting will be determined in subroutine mesub.
!                        [ (      g(h2o) ) / ( kg(dry air) day ) ] 
!!   CHECK UNITS HERE !!!!
!!!!!!!!======>>>>>>>    NOTE UNITS OF g(h2o). (Verify and change.)
!       temp_c           temperature at model levels [ deg K ]
!       h1_2             vertical entropy flux convergence profile on
!                        lo-res grid for the current ensemble member 
!                        [ deg K / sec ]                        
!       ecd              profile of condensate evaporated in convective
!                        downdraft on large-scale model grid
!                        [ g(h2o) / kg(air) / day ]
!       ece              profile of condensate evaporated in convective
!                        updraft on large-scale model grid
!                        [ g(h2o) / kg(air) / day ]
!       evap_rate        cloud-area-weighted profile of the potential
!                        cloud water evaporation for the current ensemble
!                        member on the lo-res grid. this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ (      kg(h2o) ) / ( kg(dry air) sec ) ] 
!       phalf_c          pressure at lo-res model interface levels [ Pa ]
!       q1               vertical moisture flux convergence profile on 
!                        lo-res grid for the current ensemble member 
!                        [ kg(h2o) / ( kg(dry air) sec ) ]
!       h1                                   condensation rate profile
!                        on lo-res grid for the current ensemble member
!                        [ (      kg(h2o) ) / ( kg( dry air) sec ) ] 
!       pfull_c          pressure on lo-res model full levels [ Pa ]
!       meso_melt        profile of condensate which is melted in meso-
!                        scale downdraft on large-scale model grid
!                        [ g(h2o) / kg(air) / day ]
!       meso_freeze      profile of condensate which is frozen upon 
!                        entering the anvil on the large-scale grid
!                        [ g(h2o) / kg(air) / day ]
!       qtr              vertical tracer flux convergence profile on
!                        lo-res grid for the current ensemble member 
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!       lmeso            a mesoscale circulation exists in the current
!                        grid box ?
!       debug_ijt        logical indicating whether diagnostics are 
!                        requested for this column 
!       diag_unit        unit number of column diagnostics output file
!
!   intent(inout) variables:
!
!       ensmbl_melt      vertical profile on the lo-res grid of ice melt,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!       ensmbl_freeze    vertical profile on the lo-res grid of freezing,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!       enctf            vertical profile on the lo-res grid of the entr-
!                        opy forcing, consisting of the sum of the
!                        vertical entropy flux convergence and the latent
!                        heat release, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area
!                        of ensemble member #1
!                        [ deg K / day ]                        
!       encmf            vertical profile on the lo-res grid of the      
!                        moisture forcing, consisting of the sum of the
!                        vertical moisture flux convergence and the cond-
!                        ensation, summed over ensemble members # 1 to 
!                        the current, each member's contribution being 
!                        weighted by its cloud area at level k relative 
!                        to the cloud base area of ensemble member #1
!                        [ (      kg(h2o) ) / ( kg( dry air) day ) ] 
!       enev             vertical profile on the lo-res grid of the      
!                        cloud-area-weighted profile of the potential
!                        cloud water evaporation, summed over ensemble 
!                        members # 1 to the current, each member's con-
!                        tribution being weighted by its cloud area at !
!                        level k relative to the cloud base area of 
!                        ensemble member #1.  this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ (      kg(h2o) ) / ( kg(dry air) sec ) ] 
!       disg             vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation
!                        associated with the evaporation of condensate
!                        in the convective downdraft and updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!       disb             vertical profile on the lo-res grid of the      
!                        temperature flux convergence, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area 
!                        of ensemble member #1.  
!                        [ deg K / day ] 
!       disc             vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation, 
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!       ecds             vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective downdraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!       eces             vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!       disd             vertical profile on the lo-res grid of the      
!                        vertical moisture flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        lo-res grid for the current ensemble member 
!                        [  g(h2o) / ( kg(dry air) day ) ]
!       qtren            vertical profile on the lo-res grid of the      
!                        vertical tracer flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!
!    intent(out) variables:
!
!       ermesg           error message produced by any kernel routines
!                        called by this subroutine 
!
!---------------------------------------------------------------------- 

!--------------------------------------------------------------------
!   local variables:


      real, dimension (nlev_lsm) :: rlh  
                                     !  condensation term in temperature
                                     !  equation on lo-res grid for cur-
                                     !  rent ensemble member 
                                     !  [ deg K / day ]
      real, dimension (nlev_lsm) :: cmf 
                                     !  forcing term for moisture 
                                     !  equation on lo-res grid for
                                     !  current ensemble member; sum 
                                     !  of vertical flux convergence 
                                     !  and condensation terms 
                                     !  [ g(h2o) / ( kg(air) day ) ]

     real     :: convrat   !  latent heat factor, appropriate for the 
                           !  temperature at a given model level 
                           !  ( = L / cp ) [ deg K ]
     real     :: wt_factor !  factor used to normalize each ensemble 
                           !  member's contribution to the ensemble sum, 
                           !  based on the cloud area  [ dimensionless ]
     real     :: qtrsum    !  sum of tracer flux convergence over all 
                           !  tracers and all levels and all ensemble 
                           !  members up to the current. used as a diag-
                           !  nostic; sum should be 0.0, if no boundary 
                           !  source term.
                           !  [ kg(tracer) / ( kg(dry air) sec ) ]
     integer  :: k, kcont  !  do-loop indices

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  '

!--------------------------------------------------------------------
!    sum up various cloud-base-area weighted contributions to vertical
!    profiles on the large-scale grid that are being summed over the 
!    ensemble.
!--------------------------------------------------------------------
      do k=1,nlev_lsm       

!---------------------------------------------------------------------
!    define the moisture forcing term (sum of condensation h1 and 
!    vertical flux convergence q1) on the large-scale grid. convert to
!    units of g(h20) per kg(air) per day, requiring multiplication by
!    1.0e3 g(h2o) /kg(h2o) times SECONDS_PER_DAY. add this member's 
!    contribution to the sum over the ensemble (encmf). 
!----------------------------------------------------------------------
        cmf(k) = (-h1(k) + q1(k))*(Param%SECONDS_PER_DAY*1.0e03)
        encmf(k) = encmf(k) + area_ratio*cmf(k)

!----------------------------------------------------------------------
!    define the condensation term in the temperature equation on the 
!    large-scale grid (rlh), using the latent heat of vaporization when 
!    the ambient temperature is above freezing, and the latent heat of 
!    sublimation when ice may be present. add this member's contribution
!    to the sum over the ensemble (disc).
!----------------------------------------------------------------------
        if (temp_c(k) >= Param%tfre) then
          convrat = Param%HLV/Param%CP_AIR
        else
          convrat = Param%HLS/Param%CP_AIR
        endif
        rlh(k) = h1(k)*Param%SECONDS_PER_DAY*convrat 
        disc(k) = disc(k) + area_ratio*rlh(k)

!--------------------------------------------------------------------
!    add this member's weighted contribution to the ensemble's temper-
!    ature flux convergence (disb), the ensemble's water vapor flux 
!    convergence (disd) and the ensemble's entropy flux convergence 
!    (enctf). convert the rates to units of per day, and for disd from
!    kg(h2o) per kg(air) to g(h2o) per kg(air).
!--------------------------------------------------------------------
        disb(k) = disb(k) + area_ratio*(h1_2(k)*Param%SECONDS_PER_DAY)
        disd(k) = disd(k) + area_ratio*(q1(k)*(Param%SECONDS_PER_DAY*1.0e3))
        enctf(k) = enctf(k) + area_ratio*    &
                                     (h1_2(k)*Param%SECONDS_PER_DAY + rlh(k))

!--------------------------------------------------------------------
!    if a mesoscale circulation exists, add this member's contribution
!    to the mesoscale condensate's evaporation associated with convect-
!    ive downdrafts (ecds) and that associated with evaporation into 
!    the environment (eces). if there has been no freezing associated
!    with the mesoscale condensate, define the condensation term for
!    the temperature equation using the latent heat of vaporization
!    (disg). if there has been freezing, then the convective downdraft 
!    heating uses the latent heat of vaporization, whereas the entrain-
!    ment evaporation is of ice and so uses the latent heat of 
!    sublimation.
!--------------------------------------------------------------------
        if (lmeso) then
          ecds(k) = ecds(k) + area_ratio*ecd(k)
          eces(k) = eces(k) + area_ratio*ece(k)
          if (dint == 0.) then
            disg(k) = disg(k) - area_ratio*((ecd(k) + ece(k))*  &
                                Param%hlv/(Param%CP_AIR*1000.))
          else
            disg(k) = disg(k) - area_ratio*(ece(k)*Param%HLS/  &
                      (Param%CP_AIR*1000.))
            disg(k) = disg(k) - area_ratio*(ecd(k)*Param%HLV/  &
                      (Param%CP_AIR*1000.))
          endif
        endif

!---------------------------------------------------------------------
!    add this member's cloud water evaporation rate to the sum over 
!    the ensemble (enev).
!---------------------------------------------------------------------
        enev(k) = enev(k) + area_ratio*evap_rate(k)

!--------------------------------------------------------------------
!    if a mesoscale circulation exists, add the appropriately-weighted
!    anvil freezing and melting terms to the arrays accumulating their 
!    sums over the ensemble (ensmbl_melt, ensmbl_freeze). if in a diag-
!    nostic column, output the anvil (meso_freeze) and ensemble-sum
!    (ensmbl_freeze) freezing profiles.
!--------------------------------------------------------------------
        if (lmeso) then
          ensmbl_melt(k) = ensmbl_melt(k) + area_ratio*meso_melt(k)
          ensmbl_freeze(k) = ensmbl_freeze(k) + area_ratio*meso_freeze(k)
          if (debug_ijt) then
            if (meso_freeze(k) /= 0.0) then
              write (diag_unit, '(a, i4, 2e20.12)')  &
                              'in mulsub: jk,fres,fre= ',   &
                               k, ensmbl_freeze(k), meso_freeze(k)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    add the appropriately-weighted convective cell freezing and 
!    melting terms to the arrays accumulating vertical profiles of 
!    total cloud melting (ensmbl_melt) and freezing (ensmbl_freeze) 
!    over the entire ensemble.  if in diagnostic column, output the 
!    convective cell (cell_freeze) and accumulated (ensmbl_freeze) 
!    freezing profiles.
!--------------------------------------------------------------------
        ensmbl_freeze(k) = ensmbl_freeze(k) + area_ratio*cell_freeze(k)
        ensmbl_melt(k) = ensmbl_melt(k) + area_ratio*cell_melt(k)
        if (debug_ijt) then
          if (cell_freeze(k) /= 0.0) then
            write (diag_unit, '(a, i4, 2e20.12)')  &
                     'in mulsub: jk,fres,frea= ',    &
                                     k, ensmbl_freeze(k), cell_freeze(k)
          endif
        endif
      end do

!---------------------------------------------------------------------
!    if in a diagnostics column, initialize a variable to sum the 
!    pressure-weighted tracer flux convergence, summed over all tracers
!    and all levels, for this ensemble member. upon completion of the 
!    loop, qtrsum should be 0.0, if there are no imposed tracer sources 
!    or sinks.
!---------------------------------------------------------------------
      if (debug_ijt) then
        qtrsum = 0.
        do k=1,nlev_lsm       
          do kcont=1,ntr  
            write (diag_unit, '(a,  i4, 2e20.12)')  &
                      'in mulsub: jk,    qtr,qtren=         ', &
                              k,    qtr(k,kcont),qtren(k,kcont)
            qtrsum = qtrsum + qtr(k,kcont)*(phalf_c(k) - phalf_c(k+1))
            write (diag_unit, '(a,  i4, e20.12)')  &
                         'in mulsub: jk,    qtrsum= ', k,    qtrsum
          end do
        end do
      endif 

!--------------------------------------------------------------------
!   add this ensemble member's appropriately-weighted contributions 
!   to the tracer flux convergence profiles being summed over the 
!   ensemble. 
!--------------------------------------------------------------------
      do k=1,nlev_lsm        
        qtren(k,:) = qtren(k,:) + area_ratio*qtr(k,:)
      end do

!----------------------------------------------------------------------
!    if in diagnostics column, output the profiles of the amount of
!    cloud  water evaporated (if lmeso is true) or the cloud water evap-
!    oration rate (if lmeso is false),  the array ctf (forcing for 
!    entropy eqn) and cmf (forcing for moisture equation) for this 
!    ensemble member.
!----------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm        
          if (rlh(k) /= 0.0 .or. h1_2(k) /= 0.0) then
            write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, p & ctf', k, pfull_c(    k),  &
                                               h1_2(k)*86400. + rlh(k)
          endif
        end do
        do k=1,nlev_lsm        
          if (cmf(k) /= 0.0) then
            write (diag_unit, '(a, i4, f19.10, e20.12)') &
                       'in mulsub: k, p & cmf', k, pfull_c(k), cmf(k)
          endif
        end do
      endif

!---------------------------------------------------------------------


end subroutine don_d_add_to_ensmbl_sum_lores_k 




!#####################################################################

subroutine don_d_add_to_ensmbl_sum_intgl_k        &
         (diag_unit, debug_ijt, lmeso, area_ratio, ca, cell_precip, cu, &
          apt, ensmbl_precip, ensmbl_cond, ensmbl_anvil_cond,  &  
          ensmbl_cld_top_area, ermesg)

!----------------------------------------------------------------------
!    subroutine don_d_add_to_ensmbl_sum_intgl_k adds the contrib-
!    utions from this ensemble member to various global integrals.
!----------------------------------------------------------------------

implicit none

!----------------------------------------------------------------------
integer,          intent(in   ) :: diag_unit
logical,          intent(in   ) :: debug_ijt, lmeso
real,             intent(in   ) :: area_ratio, ca, cell_precip, cu, apt
real,             intent(inout) :: ensmbl_precip, ensmbl_cond, &
                                   ensmbl_anvil_cond, ensmbl_cld_top_area
character(len=*), intent(  out) :: ermesg
 
!----------------------------------------------------------------------
!   intent(in) variables:
!
!      area_ratio       ratio of cloud base area of this ensemble member
!                       to that of ensemble member # 1. (ensemble member
!                       # 1 assumed to have largest cloud base area)
!                       [ dimensionless ]
!      ca               rate of transfer of condensate from cell to 
!                       anvil for this ensemble member 
!                       [ mm / day ]
!      cell_precip      precipitation rate for this ensemble member
!                       [ mm / day ]
!      cu               condensation rate for this ensemble member
!                       [ mm / day ]
!      apt              ratio of cloud top area to cloud base area
!                       for this ensemble member [ dimensionless ]
!      lmeso            logical indicating if mesoscale circulation 
!                       is present
!      debug_ijt        logical indicating whether diagnostics are 
!                       requested for this column 
!      diag_unit        unit number of column diagnostics output file
! 
!   intent(inout) variables:
!
!      ensmbl_precip    sum of precipitation rate over ensemble members,
!                       # 1 to the current, weighted by the area at 
!                       cloud base of each member
!                       [ mm / day ]
!      ensmbl_cond      sum of condensation rate over ensemble members,
!                       # 1 to the current, weighted by the area at 
!                       cloud base of each member
!                       [ mm / day ]
!      ensmbl_anvil_cond 
!                       sum of rate of transfer of condensate from cell 
!                       to anvil over ensemble members, # 1 to the c
!                       current, weighted by the area at cloud base of 
!                       each member [ mm / day ]
!      ensmbl_cld_top_area  
!                       sum of the cloud top areas over ensemble members 
!                       # 1 to the current, normalized by the cloud base
!                       area of ensemble member # 1 [ dimensionless ]
!
!   intent(out) variables:
!
!      ermesg           error message produced by this subroutine or any
!                       kernel routines called by this subroutine 
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  '

!--------------------------------------------------------------------
!    if a mesoscale circulation is present, add this member's cloud-
!    base_area-weighted contribution of condensate transferred to the 
!    anvil (ensmbl_anvil_cond) and cloud top cloud fraction 
!    (ensmbl_cld_top_area) to the arrays accumulating the ensemble sums.
!--------------------------------------------------------------------
      if (lmeso) then
        ensmbl_anvil_cond   = ensmbl_anvil_cond   + area_ratio*ca
        ensmbl_cld_top_area = ensmbl_cld_top_area + area_ratio*apt
      endif

!--------------------------------------------------------------------
!    add this ensemble member's weighted contribution to the total 
!    precipitation (ensmbl_precip) and condensation (ensmbl_cond). 
!--------------------------------------------------------------------
      ensmbl_precip = ensmbl_precip + area_ratio*cell_precip
      ensmbl_cond   = ensmbl_cond   + area_ratio*cu

!---------------------------------------------------------------------



end subroutine don_d_add_to_ensmbl_sum_intgl_k 

!#####################################################################


!######################################################################

subroutine don_d_output_diag_profs_k    &
         (nlev_lsm, diag_unit, pfull_c, disc, disb, disd, disn, encmf,&
          temp_tend_freeze, temp_tend_melt, cmus_tot, emds, emes, wmms, &
          wmps, tmes, mrmes, eces, ecds, disa, dise, disg_2, disf, &
          ermesg)

!---------------------------------------------------------------------
!    subroutine output_diagnostic_profiles prints out vertical profiles
!    of various donner_deep variables in those columns for which 
!    diagnostics have been requested.
!---------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------
integer,                      intent(in)   :: nlev_lsm, diag_unit
real,    dimension(nlev_lsm), intent(in)   :: pfull_c, disc, disb, disd,&
                                              disn, encmf,  &
                                              temp_tend_freeze,    &
                                              temp_tend_melt, cmus_tot, &
                                              emds, emes, wmms, wmps, &
                                              tmes, mrmes, eces, ecds, &
                                              disa, dise, disg_2, disf 
character(len=*),              intent(out) :: ermesg

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     diag_unit
!     pfull_c
!     disc
!     disb
!     disd
!     disn
!     encmf
!     temp_tend_freeze
!     temp_tend_melt
!     cmus_tot
!     emds
!     emes
!     wmms
!     wmps
!     tmes
!     mrmes
!     eces
!     ecds
!     disa
!     dise
!     disg_2
!     disf
!
!----------------------------------------------------------------------

!!  UNITS
!!    ucemh  [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    enctf  [ deg K / day ]
!!    encmf  [ g(h2o) / kg(air) /day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]
!!    disg   [ deg K / day ]



!-------------------------------------------------------------------
!   local variables:

      integer  ::  k      ! do-loop index
     
!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!  disc: cloud ensemble cell condensation heating rate [ deg K / day ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disc(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
              'in mulsub: k, P & LHR =',  k, pfull_c(k),disc(k)
        endif
      end do

!---------------------------------------------------------------------
!  disb  : cloud ensemble cell vertical entropy flux convergence 
!          [ deg K / day ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disb(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)') &
              'in mulsub: k, P & EFC =',  k, pfull_c(k),disb(k)
        endif
      end do

!---------------------------------------------------------------------
!  disd: cloud ensemble cell vertical moisture flux convergence 
!        [ kg(h2o)/ (kg(air) sec) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disd(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
               'in mulsub: k, P & EMF =',  k, pfull_c(k),disd(k)
        endif
      end do

!---------------------------------------------------------------------
!  disn : cloud ensemble cell thermal forcing 
!         [ deg K / sec ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disn(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & cell thermal forcing =',    &
                        k, pfull_c(k),disn(k)
        endif
      end do

!---------------------------------------------------------------------
!  encmf : cloud ensemble cell moisture forcing 
!          [ kg(h2o) / (kg(air) sec) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (encmf(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & cell moisture forcing =',    &
                        k, pfull_c(k),encmf(k)
        endif
      end do

!---------------------------------------------------------------------
!  temp_tend_freeze  : cloud ensemble temperature tendency due to   
!                      freezing [ deg K / sec ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (temp_tend_freeze(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)') &
                      'in mulsub: k, P & meso up freeze        =',    &
                        k, pfull_c(k),temp_tend_freeze(k)
        endif
      end do

!---------------------------------------------------------------------
!  temp_tend_melt  : cloud ensemble plus mesoscale temperature 
!                    tendency due to melting [ deg K / sec ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (temp_tend_melt(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                     'in mulsub: k, P & meso down melt        =',    &
                       k, pfull_c(k),temp_tend_melt(k)
        endif
      end do

!---------------------------------------------------------------------
!  cmus_tot : water mass condensed in mesoscale updraft
!             [ kg(h2o) / (kg(air) sec) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (cmus_tot(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)') &
                      'in mulsub: k, P & meso up con           =',    &
                         k, pfull_c(k),cmus_tot  (    k)
        endif
      end do

!---------------------------------------------------------------------
!  emds : evaporation in mesoscale downdrafts.
!         [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (emds(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                     'in mulsub: k, P & meso down evap        =',    &
                       k, pfull_c(k), emds(k)
        endif
      end do

!---------------------------------------------------------------------
!  emes : evaporation in mesoscale updrafts.
!         [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (emes(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                        'in mulsub: k, P & meso up evap        =',    &
                          k, pfull_c(k),emes(k)
        endif
      end do

!---------------------------------------------------------------------
!  wmms : vapor transferred from cell to mesoscale circulation and 
!         then condensed [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (wmms(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & meso cell con       =',    &
                       k, pfull_c(k),wmms(k)
        endif
      end do

!---------------------------------------------------------------------
!  wmps : vapor transferred from cell to mesoscale circulation  
!         [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (wmps(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & meso vap redist     =',    &
                          k, pfull_c(k),wmps(k)
        endif
      end do

!---------------------------------------------------------------------
!  tmes   : mesoscale temperature flux convergence                
!           [ deg K / day ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (tmes(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & meso efc            =',    &
                        k, pfull_c(k),tmes(k)
        endif
      end do

!---------------------------------------------------------------------
!  mrmes   : mesoscale moisture flux convergence                
!            [ kg(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
!! WAS ORIGINALLY    :    if (tmes  (    k) /= 0.00 ) then
        if (mrmes(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)') &
                        'in mulsub: k, P & meso mfc            =',    &
                         k, pfull_c(k),mrmes(k)
        endif
      end do

!---------------------------------------------------------------------
!  eces   : sublimation in mesoscale updrafts               
!           [ kg(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (eces(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & up con evap         =',    &
                          k, pfull_c(k), eces(k)
        endif
      end do

!---------------------------------------------------------------------
!  ecds_v : sublimation in mesoscale downdrafts             
!           [ kg(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (ecds(k) /= 0.00 ) then
           write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & down con evap         =',    &
                       k, pfull_c(k), ecds(k)
        endif
      end do

!---------------------------------------------------------------------
!  disa   : total temperature tendency due to deep convection
!           [ deg K / day ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disa(k) /= 0.0) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                           'in mulsub: k, p & ens thermal forc', &
                              k, pfull_c(k), disa(k)
         endif
       end do

!---------------------------------------------------------------------
!  dise : total moisture tendency due to deep convection
!         [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (dise(k) /= 0.0) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                          'in mulsub: k, p & ens moisture forc', &
                           k, pfull_c(k), dise(k)
        endif
      end do

!---------------------------------------------------------------------
!  disg2_v : total moisture tendency due to remaining cell condensate
!            [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disg_2(k) /= 0.0) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                         'in mulsub: k, p & thermal modifications', &
                          k, pfull_c(k), disg_2(k)
        endif
      end do

!---------------------------------------------------------------------
!  disf_v : total moisture tendency due to evaporation in cell updrafts
!           [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disf(k) /= 0.0) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, p & moisture modifications', &
                       k, pfull_c(k), disf(k)
        endif
      end do

!---------------------------------------------------------------------



end subroutine don_d_output_diag_profs_k


!#####################################################################

subroutine don_d_def_conv_forcing_k  &
         (nlev_lsm, diag_unit, debug_ijt, lmeso, Param, ensmbl_precip,  &
          meso_precip, meso_cloud_area, anvil_precip_melt, phalf_c,  &
          enev, encmf, ensmbl_freeze, enctf, disg, ecds, eces, emds, &
          emes, mrmes, tmes, wmps, ensmbl_cloud_area, ensmbl_melt,   &
          pfull_c, temp_c, cmus_tot, wmms, disc, disb, disd,  &
          total_precip_c, disf, disg_2, disn, dise, disa, cutotal,   &
          temp_tend_melt, temp_tend_freeze, ermesg)

!---------------------------------------------------------------------
!    subroutine define_convective_forcing produces the effects of
!    the donner_deep parameterization on the large-scale flow, defining
!    the time tendency terms and integral quantities resulting from the
!    parameterization.
!---------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_param_type

implicit none

!---------------------------------------------------------------------
integer,                      intent(in)  :: nlev_lsm, diag_unit
logical,                      intent(in)  :: debug_ijt, lmeso 
type(donner_param_type),      intent(in)  :: Param
real,                         intent(in)  :: ensmbl_precip, meso_precip
real,    dimension(nlev_lsm), intent(in)  :: meso_cloud_area,  &
                                             anvil_precip_melt, phalf_c,&
                                             enev, encmf, ensmbl_freeze,&
                                             enctf, disg, ecds, eces, &
                                             emds, emes, mrmes, tmes,  &
                                             wmps, ensmbl_cloud_area,  &
                                             ensmbl_melt, pfull_c,   &
                                             temp_c, cmus_tot, wmms,   &
                                             disc, disb, disd
real,                         intent(out) :: total_precip_c
real,    dimension(nlev_lsm), intent(out) :: disf, disg_2, disn, dise,  &
                                             disa, cutotal,   &
                                             temp_tend_melt,  &
                                             temp_tend_freeze
character(len=*),             intent(out) :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!        diag_unit         i/o unit for column diagnostics output
!        ensmbl_precip
!        meso_precip       
!        lmeso              a mesoscale circulation is present in this
!                           column ?     
!        debug_ijt          column diagnostics are desired in this 
!                           column ?
!        meso_cloud_area
!        anvil_precip_melt
!        phalf_c            pressure at large-scale model half levels 
!                           [ Pa ]
!        enev
!        encmf
!        ensmbl_freeze
!        enctf
!        disg
!        ecds
!        eces
!        emds
!        emes
!        mrmes
!        tmes
!        wmps
!        ensmbl_cloud_area
!        ensmbl_melt         
!        pfull_c
!        temp_c
!        cmus_tot
!
!    intent(out) variables:
!
!        total_precip_c
!        disf
!        disg_2
!        disn
!        dise
!        disa
!        cutotal
!        temp_tend_melt
!        temp_tend_freeze
!
!---------------------------------------------------------------------



!!  UNITS
!!    ucemh  [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    h1_2   [ deg K / sec ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    enctf  [ deg K / day ]
!!    encmf  [ g(h2o) / kg(air) /day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]
!!    disg   [ deg K / day ]


!--------------------------------------------------------------------
!   local variables:

      real    ::  esum, esuma,        esumc, sumf, summ, sumqme, sumg,&
                  sumn, sumelt, sumfre, summes, disl, disga
      real    ::  dp
      integer ::  k

!--------------------------------------------------------------------
!   local variables:
!
!         esum
!         esuma
!         esumc
!         sumf
!         summ
!         sumqme
!         sumg
!         sumn
!         sumelt
!         sumfre
!         summes
!         disl
!         disga
!         nlev            number of layers in large-scale model
!         k               do-loop index
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  '

!--------------------------------------------------------------------
!    define the total precipitation (total_precip_c) from the parameter-
!    ization as the sum of the convective (ensmbl_precip) and mesoscale
!    (meso_precip) precipitation. 
!--------------------------------------------------------------------
      total_precip_c = ensmbl_precip + meso_precip    

!----------------------------------------------------------------------
!    add the mesoscale cloud area to the cell-ensemble cloud area to
!    obtain the total cloud area profile.    
!----------------------------------------------------------------------
      do k=1,nlev_lsm
        cutotal (k) = ensmbl_cloud_area(k) + meso_cloud_area(k)
      end do

!---------------------------------------------------------------------
!    if in a diagnostics column, output the profiles of ensemble-total 
!    cloud area (ensmbl_cloud_area) and mesoscale cloud area 
!    (meso_cloud_area), total cloud area (cu_total), deposition in 
!    mesoscale updrafts (cmus), evaporation in mesoscale downdrafts 
!    (emds), evaporation from mesoscale updrafts (emes), water vapor 
!    supplied to mesoscale circulation (wmps), melted anvil precip-
!    itation ( anvil_precip_melt), mesoscale temperature flux
!    convergence (tmes) and mesoscale vapor-flux convergence (mrmes).
!---------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm
          write (diag_unit, '(a, i4, f19.10, 3e20.12)') &
                 'in mulsub: jk, pr,cual,cuml, cutot= ', k, pfull_c(k), &
                  ensmbl_cloud_area (k), meso_cloud_area(k), cutotal(k)
          write (diag_unit, '(a, i4, 3e20.12)')  &
                     'in mulsub: jk,cmu,emd,eme= ', k, cmus_tot(k), &
                      emds(k), emes(k)
          write (diag_unit, '(a, i4, 2e20.12)') &
                     'in mulsub: jk,wmm,wmp,elt= ', k,           &
                      wmps(k),  anvil_precip_melt(k)
          write (diag_unit, '(a, i4, f20.14, e20.12)')  &
                      'in mulsub: jk,tmes,qmes= ', k, tmes(k), mrmes(k)
        end do
      endif

!---------------------------------------------------------------------
!    define terms which will appear in the large-scale model equations.
!---------------------------------------------------------------------
      do k=1,nlev_lsm

!----------------------------------------------------------------------
!    combine several of the moisture tendency terms associated with the
!    donner_deep parameterization (disf). if a mesoscale circulation is
!    present, the following terms are included : 1) transfer of vapor 
!    from mesoscale to large-scale flow (cmus_tot), 2) evaporation in 
!    cumulus downdrafts (ecds), 3) evaporation from cumulus updrafts 
!    (eces), 4)  vapor transferred from cells to mesoscale (wmps), 5) 
!    evaporation from mesoscale updrafts (emes), 6) evaporation from 
!    mesoscale downdrafts (emds), and 7) mesoscale moisture-flux 
!    convergence (mrmes). 
!----------------------------------------------------------------------
        if (lmeso) then
          disf(k) = -cmus_tot(k) + ecds(k) + eces(k) + wmps(k) +  &
                    emes(k) + emds(k) + mrmes(k)

!----------------------------------------------------------------------
!    if a mesoscale circulation is not present, disf is simply the 
!    moisture tendency associated with the evaporation of the condensed
!    cloud water that did not precipitate out (enev). convert to units 
!    of g(h2o) per kg(air) per day.
!----------------------------------------------------------------------
        else
          disf(k) = enev(k)*(1.0E03*Param%seconds_per_day)
        endif

!---------------------------------------------------------------------
!    define the sum of disf and the term containing the tendency due 
!    to cell-scale vertical moisture-flux convergence and associated
!    condensation (encmf), and store in array dise.
!---------------------------------------------------------------------
        dise(k) = encmf(k) + disf(k)

!----------------------------------------------------------------------
!    define the temperature tendencies associated with the freezing
!    of updraft liquid (temp_tend_freeze) and the melting of ice
!    falling from the anvil (temp_tend_melt). combine several of the 
!    temperature tendencies associated with the cell component of the
!    donner_deep parameterization (disn). disn is composed of 1) a term
!    combining the vertical flux convergence of temperature and cloud 
!    condensation (enctf), 2) evaporation of liquid in the cell updraft
!    and downdraft (disg), 3) the freezing of updraft liquid 
!    (temp_tend_freeze) and 4) the melting of ice (temp_tend_melt). 
!    separately define the temperature tendency resulting from the 
!    latent heat release associated with sublimation occurring in the 
!    mesoscale updraft and downdraft (disga).
!----------------------------------------------------------------------
        temp_tend_freeze (k) = ensmbl_freeze(k)*Param%hlf/     &
                               (Param%cp_air*1000.)
        temp_tend_melt(k) = -(ensmbl_melt(k) + anvil_precip_melt(k))*  &
                                          Param%hlf/(Param%cp_air*1000.)
        disn(k) = enctf(k) + disg(k) + temp_tend_freeze(k) + &
                  temp_tend_melt(k)
        disga = (emes(k) + emds(k))*Param%hls/(Param%cp_air*1000.)

!--------------------------------------------------------------------
!    if a mesoscale circulation is present, define the heating terms
!    equivalent to the disf array (disg_2). included in disg_2 are 
!    terms associated with 1) transfer of vapor from mesoscale to 
!    large-scale flow (disl), 2) evaporation in cumulus updrafts and 
!    downdrafts (disg), 3) freezing of liquid in the updraft 
!    (temp_tend_freeze), 4) melting of ice (temp_tend_melt), 5) evap-
!    oration in the mesoscale circulation (disga), and 6) mesoscale 
!    temperature flux convergence (tmes).
!--------------------------------------------------------------------
        if (lmeso) then
          disl = cmus_tot(k)*Param%hls/(Param%cp_air*1000.)
          disg_2(k) = disl + disg(k) + temp_tend_freeze (k) +  &
                      temp_tend_melt(k) - disga + tmes(k)

!----------------------------------------------------------------------
!    if a mesoscale circulation is not present, disg_2 is simply the 
!    temperature tendency associated with the evaporation of the 
!    condensed cloud water that did not precipitate out (disf). the
!    latent heat constant is chosen appropriately for the in situ
!    temperature.
!----------------------------------------------------------------------
        else
          if (temp_c(k) .gt. Param%tfre) then
            disg_2(k) = -disf(k)*Param%hlv/(Param%cp_air*1000.)
          else
            disg_2(k) = -disf(k)*Param%hls/(Param%cp_air*1000.)
          endif
        endif

!---------------------------------------------------------------------
!    define the sum of disg_2 and the term containing the tendency due 
!    to cell-scale vertical temperature-flux convergence and associated
!    condensation (enctf), and store in array disa.
!---------------------------------------------------------------------
        disa(k) = enctf(k) + disg_2(k)

!--------------------------------------------------------------------
!    if in a diagnostics column, output the profile of temperature 
!    change associated with evaporation in the mesoscale circulation 
!    (disga).
!--------------------------------------------------------------------
        if (debug_ijt) then
          if (disga /= 0.0) then
            write (diag_unit, '(a, i4, f19.10,  e20.12)')  &
                    'in mulsub: jk,pr,emds,disga= ', k, pfull_c(k), disga
          endif
        endif
      end do

!--------------------------------------------------------------------
!    if in a diagnostics column, compute the column integrals of the 
!    various tendency terms for the vapor and temperature equations. 
!    esum  : total vapor tendency from donner_deep parameterization
!    sumf  : total vapor tendency less the vertical flux convergence and
!            condensation      
!    summ  : vapor tendency due to vertical flux convergence and
!            condensation 
!    sumqme: mesoscale moisture flux convergence
!    esuma : total temperature tendency from donner_deep parameter-
!            ization
!    sumg  : total temperature tendency less the vertical flux conver-
!            gence and condensation
!    esumc : temperature tendency due to vertical flux convergence and
!            condensation  
!    summes: temperature tendency due to mesoscale temperature flux
!            convergence
!    sumelt: temperature tendency due to melting within column
!    sumfre: temperature tendency due to freezing within the column
!    sumn  : temperature tendency associated with the cell component of
!            the donner_deep parameterization
!--------------------------------------------------------------------
      if (debug_ijt) then
        esum  = 0.
        sumf   = 0.
        summ   = 0.
        sumqme = 0.
        esuma = 0.
        sumg = 0.
        esumc  = 0.
        summes = 0.
        sumelt = 0.
        sumfre = 0.
        sumn = 0.
        do k=1,nlev_lsm
          dp = phalf_c(k) - phalf_c(k+1)
          esum   = esum   + dise(k)*dp                     
          sumf   = sumf   + disf(k)*dp
          summ   = summ   + encmf(k)*dp
          sumqme = sumqme + mrmes(k)*dp
          esuma  = esuma  + disa(k)*dp
          sumg   = sumg   + disg_2(k)*dp
          esumc  = esumc  + enctf(k)*dp
          summes = summes + tmes(k)*dp
          sumelt = sumelt + (ensmbl_melt(k) + anvil_precip_melt(k))*dp
          sumfre = sumfre + ensmbl_freeze(k)*dp
          sumn   = sumn   + disn(k)*dp
        end do
!---------------------------------------------------------------------
!    convert the moisture terms to units of mm(h2o) per day.
!---------------------------------------------------------------------
        esum   = esum/(Param%grav*1000.)
        sumf   = sumf/(Param%grav*1000.)
        summ   = summ/(Param%grav*1000.)
        sumqme = sumqme/(Param%grav*1000.)

!---------------------------------------------------------------------
!    convert the temperature terms to units of 
!    (kg(air)/ kg(h2o)) * mm(h2o) per day.
!---------------------------------------------------------------------
        esuma  = (esuma*Param%cp_air)/(Param%grav*Param%hlv)
        sumg   = sumg*Param%cp_air/(Param%grav*Param%hlv)
        esumc  = (esumc*Param%cp_air)/(Param%grav*Param%hlv)
        summes = summes*Param%cp_air/(Param%grav*Param%hlv)
        sumelt = sumelt/(Param%grav*1000.)
        sumfre = sumfre/(Param%grav*1000.)
        sumn   = sumn*Param%cp_air/(Param%grav*Param%hlv)

!--------------------------------------------------------------------
!    output the various column integrals.
!--------------------------------------------------------------------
        write (diag_unit, '(a, e20.12, a)') &
             'in mulsub: ESUM=', esum , ' MM/DAY'
        write (diag_unit, '(a, e20.12, a)') &
              'in mulsub: SUMF= ', sumf, ' MM/DAY'
        write (diag_unit, '(a, e20.12, a)') &
             'in mulsub: SUMM= ', summ, ' MM/DAY'
        write (diag_unit, '(a, e20.12, a)') &
              'in mulsub: sumqme= ', sumqme, ' mm/day'

        write (diag_unit, '(1(a,e20.12))')  &
                  'in mulsub:  ESUMA=',ESUMA
        write (diag_unit, '(a, e20.12,a)')  &
                                'in mulsub: SUMG=',SUMG,' MM/DAY'
        write (diag_unit, '(a, e20.12)') 'in mulsub: ESUMC=',ESUMC
        write (diag_unit, '(a,e20.12,a)')   &
                                'in mulsub: summes=',summes,' mm/day'
        write (diag_unit, '(a, 2e20.12,a)')  &
                     'in mulsub: sumelt,sumfre= ',sumelt,sumfre,  &
                                                ' mm/day'
        write (diag_unit, '(a,e20.12,a)')  &
                                'in mulsub: SUMN= ',SUMN,' MM/DAY'
      endif

!---------------------------------------------------------------------
!    call subroutine output_diagnostic_profiles to print various 
!    output fields from the donner_deep parameterization in those 
!    columns for which diagnostics have been requested.
!---------------------------------------------------------------------
      if (debug_ijt) then
         call don_d_output_diag_profs_k    &
              (nlev_lsm, diag_unit, pfull_c,  disc, disb, disd, disn,  &
               encmf, temp_tend_freeze, temp_tend_melt, cmus_tot, emds, &
               emes, wmms, wmps, tmes, mrmes, eces, ecds, disa, &
               dise, disg_2, disf, ermesg)
      endif

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!----------------------------------------------------------------------


end subroutine don_d_def_conv_forcing_k



!#####################################################################

subroutine don_d_finalize_output_fields_k   &
         (nlev_lsm, ntr, i, j, Param, disb, disc, temp_tend_freeze, &
          temp_tend_melt, tmes, disd, cmus_tot, ecds, eces, emds, emes, &
          wmms, wmps, mrmes, cutotal, dmeml, detmfl, uceml, umeml, cuq, &
          cuql_v, qtren, qtmes, wtp, Don_conv, ermesg)

!----------------------------------------------------------------------
!    subroutine finalize_output_fields stores output variables from 
!    columns with active deep convection into the appropriate elements 
!    of the donner_conv_type variable Don_conv.
!----------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_param_type, donner_conv_type 

implicit none

!---------------------------------------------------------------------
integer,                          intent(in)    :: nlev_lsm, ntr
integer,                          intent(in)    :: i, j
type(donner_param_type),          intent(in)    :: Param
real,    dimension(nlev_lsm),     intent(in)    :: disb, disc,   &
                                                   temp_tend_freeze, &
                                                   temp_tend_melt,  &
                                                   tmes, disd, cmus_tot,&
                                                   ecds, eces, emds,   &
                                                   emes, wmms, wmps,  &
                                                   mrmes, cutotal, &
                                                   dmeml, detmfl, uceml,&
                                                   umeml, cuq, cuql_v
real,    dimension(nlev_lsm,ntr), intent(in)    :: qtren, qtmes, wtp
type(donner_conv_type),           intent(inout) :: Don_conv
character(len=*),                 intent(out)   :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       i, j           i, j indices of the current grid column
!       wmms
!       wmps
!       mrmes
!       emds
!       emes
!       ecds
!       eces
!       disd
!       cmus_tot
!       disb
!       disc
!       temp_tend_melt
!       temp_tend_freeze
!       tmes
!       cutotal
!       cuq
!       cuql_v
!       dmeml
!       detmfl
!       uceml
!       umeml
!       exit_flag
!       total_precip
!       meso_precip
!       qtren
!       qtmes
!       wtp
!
!   intent(inout) variables:
!
!       Don_conv       donner_convection_type derived type variable 
!                      containing fields produced by the donner_deep
!                      convection mod 
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer  :: k, kinv    ! do-loop indices

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!    if deep convection occurred in this column, save various output
!    fields. if it did not, then these components of the Don_conv
!    derived-type variable will retain their default values.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    save the following arrays as elements of the donner_conv type 
!    variable Don_Conv. make sure all arrays are in mks units, requiring
!    conversion of some arrays from per day to per second, g(h2o) to 
!    kg(h2o) and / or mm to m. reverse the vertical index, making these
!    profile arrays compatible with the large-scale model grid ( index 1
!    nearest upper boundary) rather than the cloud model grid (index 1 
!    nearest sfc).
!---------------------------------------------------------------------
      do k=1,nlev_lsm            
        kinv = nlev_lsm + 1 - k
        Don_conv%ceefc (i,j,kinv)   = disb(k)/Param%seconds_per_day
        Don_conv%cecon (i,j,kinv)   = disc(k)/Param%seconds_per_day
        Don_conv%tmes  (i,j,kinv)   = tmes(k)/Param%seconds_per_day
        Don_conv%fre   (i,j,kinv)   = temp_tend_freeze(k)/  &
                                                 Param%seconds_per_day
        Don_conv%elt   (i,j,kinv)   = temp_tend_melt(k)/     &
                                                   Param%seconds_per_day
        Don_conv%cmus  (i,j,kinv)   = cmus_tot(k)/       &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%ecds  (i,j,kinv)   = ecds(k)/           &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%eces  (i,j,kinv)   = eces(k)/            &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%emds  (i,j,kinv)   = emds(k)/            &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%emes  (i,j,kinv)   = emes(k)/            &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%mrmes  (i,j,kinv)   = mrmes(k)/          &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%wmps  (i,j,kinv)   = wmps(k)/             &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%wmms  (i,j,kinv)   = wmms(k)/                 &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%cemfc (i,j,kinv)   = disd(k)/                 &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%cual  (i,j,kinv)   = cutotal(k)
        Don_conv%dmeml (i,j,kinv)   = dmeml(k)
        Don_conv%uceml (i,j,kinv)   = uceml(k)
        if (detmfl(k) <= 1.0e-10) then
          Don_conv%detmfl(i,j,kinv) = 0.
        else
          Don_conv%detmfl(i,j,kinv)   = detmfl(k)
        endif
        Don_conv%umeml (i,j,kinv)   = umeml(k)
        Don_conv%cuqi  (i,j,kinv)   = cuq(k)
        Don_conv%cuql  (i,j,kinv)   = cuql_v(k)
        Don_conv%qtren1(i,j,kinv,:) = qtren(k,:)
        Don_conv%qtmes1(i,j,kinv,:) = qtmes(k,:)
        Don_conv%wtp1  (i,j,kinv,:) = wtp(k,:)
      end do
        

!--------------------------------------------------------------------


end subroutine don_d_finalize_output_fields_k 



!#####################################################################

!#####################################################################

subroutine don_d_determine_cloud_area_k            &
         (me, nlev_lsm, nlev_hires, diag_unit, debug_ijt, Param, Nml,  & 
          max_depletion_rate, dcape, amax, dise_v, disa_v, pfull_c,  &
          temp_c, mixing_ratio_c, env_t, env_r, parcel_t, parcel_r, &
          cape_p, exit_flag, amos, a1, ermesg)

!---------------------------------------------------------------------
!    subroutine determine_cloud_area defines the convective cloud area
!    and so closes the donner_deep parameterization. The arrays 
!    Don_conv%a1 and Don_conv%amos are output by this routine.
!---------------------------------------------------------------------

use messy_convect_donner_types_mod, only: donner_param_type, donner_nml_type 
USE MESSY_CONVECT_DONNER_UTIL,      ONLY: don_u_lo1d_to_hi1d_k
USE MESSY_CONVECT_DONNER_CLOSURE,   ONLY: cu_clo_cumulus_closure_k
implicit none

!-----------------------------------------------------------------------
integer,                      intent(in)    :: me, nlev_lsm,     &
                                               nlev_hires, diag_unit
logical,                      intent(in)    :: debug_ijt
type(donner_param_type),      intent(in)    :: Param
type(donner_nml_type),        intent(in)    :: Nml      
real,                         intent(in)    :: max_depletion_rate,   &
                                               dcape, amax
real, dimension(nlev_lsm),    intent(in)    :: dise_v, disa_v, &
                                               pfull_c, temp_c,  &
                                               mixing_ratio_c 
real, dimension(nlev_hires),  intent(in)    :: env_t, env_r, parcel_t,  &
                                               parcel_r, cape_p
logical,                      intent(inout) :: exit_flag
real,                         intent(out)   :: amos, a1
character(len=*),             intent(out)   :: ermesg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       diag_unit          unit number for column diagnostics file
!       debug_ijt          column_diagnostics are requested in 
!                          current column  ?
!       max_depletion_rate rate of moisture depletion due to convection
!                          that would result in a column without vapor
!                          [ kg(h2o) / ( kg(air) sec ) ]      
!       dcape              time tendency of cape
!       amax
!       dise_v
!       disa_v
!       pfull_c
!       temp_c 
!       mixing_ratio_c
!       env_t
!       env_r
!       parcel_t
!       parcel_r
!       cape_p
!
!   intent(inout) variables:
!
!       exit_flag
!
!   intent(out) variables:
!
!       amos
!       a1
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:
 
      real, dimension (nlev_lsm)         :: a1_vk              
      real, dimension(nlev_hires)        :: qli0_v, qli1_v, qt_v,  &
                                            qr_v, rl_v, ri_v
      real                               :: qtest, tfint, disbar
      integer                            :: k
!----------------------------------------------------------------------
!   local variables:
!
!         a1_vk
!         qli0      normalized component of cumulus condensate forcing
!         qli1      un-normalized component of condensate forcing
!         qt_v      temperature tendency due to deep convection on
!                   cape grid [ deg K / sec ]
!         qr_v      vapor mixing ratio tendency due to deep convection
!                   on cape grid [ kg(h2o) / ( kg(air) sec ]
!         rl_v      large-scale liquid mixing ratio
!         ri_v      large-scale ice mixing ratio 
!         qtest
!         tfint     column integral of moisture time tendency due to
!                   convection  [ mm / sec , or  kg / (m**2 sec ) ]
!         disbar    water vapor time tendency due to deep convection at 
!                   large-scale model interface levels
!                   [ kg(h2o) / ( kg(air) sec ) ]
!         nlev      number of layers in large-scale model
!         k         do-loop index

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!    call map_lo_res_col_to_hi_res_col to interpolate moisture and
!    temperature forcings from large-scale model grid (dise_v, disa_v)
!    to the vertical grid used in the cape calculation (qr_v, qt_v). 
!--------------------------------------------------------------------
      call don_u_lo1d_to_hi1d_k   &
            (nlev_lsm, nlev_hires, disa_v, pfull_c, cape_p, qt_v, ermesg)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

      call don_u_lo1d_to_hi1d_k   &
            (nlev_lsm, nlev_hires, dise_v, pfull_c, cape_p, qr_v, ermesg)


!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!--------------------------------------------------------------------
!    if in a diagnostic column, output the temperature and moisture 
!    forcings on both the cape grid (qt_v, qr_v) and the large-scale
!    model grid (disa_v, dise_v).
!--------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_hires
          if (qr_v(k) /= 0.0 .or. qt_v(k) /= 0.0) then 
            write (diag_unit, '(a, i4, e20.12, f20.14)')  &
                     'in cupar: k,qr,qt= ',k, qr_v(k), qt_v(k)
          endif
        end do
        do k=1,nlev_lsm
          if (dise_v(k) /= 0.0 .or. disa_v(k) /= 0.0) then 
            write (diag_unit, '(a, i4, 2e20.12)')  &
                    'in cupar: k,dise,disa= ',k, dise_v(k), disa_v(k)
          endif
        end do
      endif

!--------------------------------------------------------------------
!   define condensate variables on the cape grid (qli0, qli1, rl_v, 
!   ri_v). these variables are not used in the current version of the
!   cumulus closure scheme implemented in subroutine cumulus_closure, 
!   so they are given values of 0.0.
!--------------------------------------------------------------------
      do k=1,nlev_hires
        qli0_v(k) = 0.
        qli1_v(k) = 0.
        rl_v(k)   = 0.
        ri_v(k)   = 0.
      end do

!--------------------------------------------------------------------
!    call subroutine cumulus_closure to determine cloud base cloud
!    fraction and so close the deep-cumulus parameterization.
!--------------------------------------------------------------------
      call cu_clo_cumulus_closure_k   &
           (nlev_hires, diag_unit, debug_ijt, Param, dcape, &
            cape_p, qli0_v, qli1_v, qr_v, qt_v, env_r, ri_v, &
            rl_v, parcel_r, env_t, parcel_t, a1, ermesg)     

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') return

!--------------------------------------------------------------------
!    calculate the vertical integral of normalized moisture forcing 
!    in the column (tfint) in units of kg (h2o) per m**2 per second, or
!    mm (h2o) per second.
!-------------------------------------------------------------------
      tfint = 0.0
      do k=2,nlev_lsm
        disbar = 0.5*(dise_v(k-1) + dise_v(k))
        tfint = tfint - disbar*(pfull_c(k-1) - pfull_c(k))
      end do
      tfint = tfint/Param%grav

!--------------------------------------------------------------------
!    restrict the cloud-base area fraction produced by subroutine
!    cumulus_closure to be no larger than the cloud base area that 
!    results in total grid box coverage at some higher level (amax). 
!--------------------------------------------------------------------
      a1 = MIN (amax, a1)
!---------------------------------------------------------------------
!    set the cloud-base area fraction to be 0.0 if there is no net
!    column integral of moisture forcing in the column. this is 
!    referred to as the moisture constraint. see "Moisture Constraint",
!    8/8/97. set the exit_flag to .true., turning off convection in
!    this column, output a message, and return to calling subprogram.
!---------------------------------------------------------------------
      if (tfint == 0.) then      
        a1 = 0.
        exit_flag      = .true.
        if (debug_ijt) then
          write (diag_unit, '(a)')  &
                 'convection turned off in column because of moist&
                  &ure constraint; cloud area being set to 0.0'
        endif
        return
      endif

!---------------------------------------------------------------------
!    if in a diagnostic column, output the column integral of the 
!    moisture forcing (tfint) and the fractional cloud area (a1) after
!    assuring that moisture forcing is present in the column.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)')  &
                      'in cupar: tfint= ',tfint       
        write (diag_unit, '(a, e20.12)')  &
                      'in cupar: a1_v = ',a1       
      endif

!---------------------------------------------------------------------
!    restrict cloud fractional area by the moisture constraint. this
!    requirement limits the cloud area so that the moisture tendency 
!    due to the deep convection (tfint - which occurs only within the 
!    cloud fractional area) will not remove more vapor from the column 
!    than is available. here amos is the cloud area over which applic-
!    ation of the convective moisture tendency will result in total
!    vapor depletion in the column.
!---------------------------------------------------------------------
      amos = max_depletion_rate/tfint     
      if (a1 > amos)  then    
        a1 = max(amos, 0.)
      endif 

!---------------------------------------------------------------------
!    for any diagnostic columns in the window in which deep convection
!    was possible, output the column integral of the moisture forcing 
!    (tfint), the max cloud area allowed by the moisture constraint 
!    (amos) and the fractional cloud area after applying the moisture
!    constraint (a1).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')  &
                   'in cupar: tfint,amos,a1= ',  &
                                       tfint, amos, a1  
      endif

!---------------------------------------------------------------------
!    verify that the current value of a1 will not produce negative
!    value of vapor mixing ratio at any level in the column when the
!    convective moisture tendency is applied. determine the large-scale
!    model mixing ratio for the current value of a1 (qtest). if qtest
!    is negative at any level for this value of a1, reset the value 
!    of a1, so that no negative mixing ratios will be produced.
!--------------------------------------------------------------------
      do k=1,nlev_lsm
        qtest = mixing_ratio_c(k) + a1*Nml%donner_deep_freq*dise_v(k)
        if (qtest < 0.) then
          a1_vk(k) = -mixing_ratio_c(k)/(dise_v(k)*Nml%donner_deep_freq)
        else
          a1_vk(k) = a1     
        endif
      end do

!--------------------------------------------------------------------
!    define the a1 for the column as the smallest of those defined
!    in the column. 
!--------------------------------------------------------------------
      a1 = MINVAL (a1_vk)

!---------------------------------------------------------------------
!    if in a diagnostic column, output the final value of a1, after 
!    all necessary constraints have been applied.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)') 'in cupar: a1= ',a1        
      endif


!--------------------------------------------------------------------


end subroutine don_d_determine_cloud_area_k 





!####################################################################

!######################################################################



!######################################################################

subroutine don_d_remove_normalization_k   &
      (isize, jsize, nlev_lsm, exit_flag, Don_conv, total_precip, &
       temperature_forcing, moisture_forcing, ermesg)

!---------------------------------------------------------------------
!    subroutine remove_normalization removes the normalization by the
!    cloud base fractional area from the various convective diagnostics
!    and output fields so that they are ready fro use in the large-scale
!    model equations.
!---------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_conv_type

implicit none 

!---------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize, nlev_lsm
logical, dimension(isize,jsize),  intent(in)    :: exit_flag
type(donner_conv_type),           intent(inout) :: Don_conv
real   , dimension(isize,jsize),  intent(inout) :: total_precip
real   , dimension(isize,jsize,nlev_lsm),                 &
                                  intent(inout) :: temperature_forcing, &
                                                   moisture_forcing
character(len=*),                 intent(out)   :: ermesg

!----------------------------------------------------------------------
!   intent(in) variables:
!
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!
!   intent(inout) variables:
!    
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer :: i, j, k    ! do-loop indices

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' '

!---------------------------------------------------------------------
!    remove normalization from the cumulus diagnostics and forcing terms
!    by multiplying them by the fractional cloud base area. these values
!    thus become grid-box averages, rather than averages over the cloudy
!    area, and so are appropriate to use in the large-scale model
!    equations. 
!---------------------------------------------------------------------
      do j=1,jsize                          
        do i=1,isize

!---------------------------------------------------------------------
!    if deep convection is present in the column, denormalize the 
!    convective fields.
!---------------------------------------------------------------------
          if (.not. exit_flag(i,j)) then
            total_precip(i,j) =  total_precip(i,j)*Don_conv%a1(i,j)
            Don_conv%ampta1(i,j) =  Don_conv%ampta1(i,j)*Don_conv%a1(i,j)
            Don_conv%cell_precip(i,j) =              &
                             Don_conv%cell_precip (i,j)*Don_conv%a1(i,j)
            Don_conv%meso_precip(i,j) =              &
                             Don_conv%meso_precip (i,j)*Don_conv%a1(i,j)
            Don_conv%emdi_v(i,j) = Don_conv%emdi_v(i,j)*Don_conv%a1(i,j)
            do k=1,nlev_lsm                           
              temperature_forcing(i,j,k) =   &
                             temperature_forcing(i,j,k)*Don_conv%a1(i,j)
              Don_conv%ceefc(i,j,k) =   &
                                  Don_conv%ceefc(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cecon(i,j,k) =        &
                                  Don_conv%cecon(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cemfc(i,j,k) =      &
                                  Don_conv%cemfc(i,j,k)*Don_conv%a1(i,j)
              moisture_forcing(i,j,k) =      &
                                moisture_forcing(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cual (i,j,k) =       &
                                   Don_conv%cual(i,j,k)*Don_conv%a1(i,j)
              Don_conv%fre(i,j,k) = Don_conv%fre(i,j,k)*Don_conv%a1(i,j)
              Don_conv%elt(i,j,k) = Don_conv%elt(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cmus(i,j,k) =      &
                                   Don_conv%cmus(i,j,k)*Don_conv%a1(i,j)
              Don_conv%ecds(i,j,k) =      &
                                   Don_conv%ecds(i,j,k)*Don_conv%a1(i,j)
              Don_conv%eces(i,j,k) =      &
                                   Don_conv%eces(i,j,k)*Don_conv%a1(i,j)
              Don_conv%emds(i,j,k) =       &
                                   Don_conv%emds(i,j,k)*Don_conv%a1(i,j)
              Don_conv%emes(i,j,k) =       &
                                   Don_conv%emes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%mrmes(i,j,k) =       &
                                   Don_conv%mrmes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%wmps(i,j,k) =       &
                                   Don_conv%wmps(i,j,k)*Don_conv%a1(i,j)
              Don_conv%wmms(i,j,k) =      &
                                   Don_conv%wmms(i,j,k)*Don_conv%a1(i,j)
              Don_conv%tmes(i,j,k) =      &
                                   Don_conv%tmes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%dmeml(i,j,k) =      &
                                  Don_conv%dmeml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%uceml(i,j,k) =      &
                                  Don_conv%uceml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%detmfl(i,j,k) =      &
                                  Don_conv%detmfl(i,j,k)*Don_conv%a1(i,j)
              Don_conv%umeml(i,j,k) =      &
                                  Don_conv%umeml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%qtren1(i,j,k,:) =     &
                               Don_conv%qtren1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%qtmes1(i,j,k,:) =     &
                               Don_conv%qtmes1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%wtp1(i,j,k,:) =       &
                                 Don_conv%wtp1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%qtceme(i,j,k,:) =   &
                     Don_conv%qtmes1(i,j,k,:) + Don_conv%qtren1(i,j,k,:)
            end do

!---------------------------------------------------------------------
!    if deep convection is not present in the column, define the output
!    fields appropriately.
!---------------------------------------------------------------------
          else
            total_precip(i,j) = 0.
            do k=1,nlev_lsm                           
              temperature_forcing(i,j,k) = 0.
              moisture_forcing(i,j,k) = 0.
            end do
          endif
        end do
      end do

!---------------------------------------------------------------------


end subroutine don_d_remove_normalization_k



!######################################################################

subroutine don_d_output_cupar_diags_k    &
         (isize, jsize, nlev_lsm, Col_diag, n, exit_flag, &
          total_precip, temperature_forcing, Don_conv, Don_cape, ermesg)

!----------------------------------------------------------------------
!----------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_conv_type, donner_cape_type, &
                             donner_column_diag_type

implicit none

!----------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize,  &
                                                   nlev_lsm
type(donner_column_diag_type),    intent(in)    :: Col_diag
integer,                          intent(in)    :: n
logical, dimension(isize,jsize),  intent(in)    :: exit_flag
real, dimension (isize,jsize),    intent(in)    :: total_precip
real, dimension (isize,jsize,nlev_lsm),                      &
                                  intent(in)    :: temperature_forcing
type(donner_conv_type),           intent(inout) :: Don_conv
type(donner_cape_type),           intent(inout) :: Don_cape
character(len=*),                 intent(out)   :: ermesg

      integer  :: idiag, jdiag, unitdiag
      integer  :: i,j,k


!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' '
      idiag = Col_diag%i_dc(n)
      jdiag = Col_diag%j_dc(n)
      unitdiag = Col_diag%unit_dc(n)

!---------------------------------------------------------------------
!    find any other columns in the current physics window in which deep
!    convection has produced a precipitation rate of over 1 mm/day. 
!    output the precipitation rate and cloud areas for each of these 
!    columns.
!---------------------------------------------------------------------
      do j=1,jsize
        do i=1,isize       
          if (.not. exit_flag(i,j) ) then
            if (Don_conv%cell_precip(i,j) > 1.) then
              write (unitdiag, '(a)')  &
                     ' the following columns in the current physics&
                         & window contain deep convection producing&
                         & rainfall rates of over 1mm per day '
              write (unitdiag, '(a, 2i4, 3e20.12)')  &
                        'in cupar: i,j, precip rate, cloud area, &
                                      &anvil area = ,',  &
                               i, j, Don_conv%cell_precip(i,j),    &
                               Don_conv%a1(i,j), Don_conv%ampta1(i,j)
            endif
          endif
        end do
      end do

!---------------------------------------------------------------------
!    if in a diagnostic window, output convection-related upper tropo-
!    spheric heating rates if there is convective precipitation in any
!    of the diagnostic columns.
!---------------------------------------------------------------------
      if (total_precip(idiag,jdiag) /= 0.) then
        do k=Col_diag%kstart,nlev_lsm
          if ((Don_cape%model_p(idiag,jdiag,k) > 100.e02) .and.&
              (Don_cape%model_p(idiag,jdiag,k) < 500.e02)) then 
            if (temperature_forcing(idiag,jdiag,nlev_lsm-k+1) /= 0.) then
              write (unitdiag, '(a, 3i4, f20.14)')    &
                     'in cupar: j_dc,i_dc,k,t= ',  &
                               jdiag, idiag, k,    &
                                    Don_cape%model_t(idiag,jdiag,k)
              write (unitdiag, '(a, e20.12, i4, 2e20.12)')&
                     'in cupar: tprea1,k,pr,cemetf= ',  &
                            total_precip(idiag,jdiag), k,    &
                            Don_cape%model_p(idiag,jdiag,k),   &
                       temperature_forcing(idiag,jdiag,nlev_lsm-k+1 )
            endif
          endif
        end do
      endif

!----------------------------------------------------------------------
!    if in a diagnostic window, output values of convective and total 
!    precipitation and cloud areas,
!----------------------------------------------------------------------
      if (.not. exit_flag(idiag,jdiag) ) then
        write (unitdiag, '(a, 2e20.12)')  &
                     'in cupar: contot,tpre=', &
                 Don_conv%cell_precip(idiag,jdiag) /  &
                                           (total_precip(idiag,jdiag)),&
                        total_precip(idiag,jdiag)
        write (unitdiag, '(a, 2e20.12)') 'in cupar: a1,ampt =',  &
                         Don_conv%a1 (idiag,jdiag), &
                         Don_conv%ampta1(idiag,jdiag)
        write (unitdiag, '(a, e20.12)')  'in cupar: amax= ', &
                          Don_conv%amax(idiag,jdiag)

!----------------------------------------------------------------------
!    if in a diagnostic window, output values of mesoscale and 
!    cell-scale mass fluxes.
!----------------------------------------------------------------------
        do k=Col_diag%kstart,nlev_lsm
          write (unitdiag, '(a, i4, f19.10, 3e20.12)')  &
                 'in cupar: k,pr,uceml,dmeml,umeml= ',  &
                     k,  Don_cape%model_p(idiag,jdiag,nlev_lsm-k+1),  &
                         Don_conv%uceml(idiag,jdiag,k), &
                         Don_conv%dmeml(idiag,jdiag,k),  &
                         Don_conv%umeml(idiag,jdiag,k)
        end do

!----------------------------------------------------------------------
!    if in a diagnostic window, output values of cloud liquid (cuql).
!    at any levels at which heating associated with the donner deep
!    convection is greater than 0.002 deg K / sec, output heating rate,
!    cloud area, cape, cape tendency, and cloud area.
!----------------------------------------------------------------------
        do k=Col_diag%kstart,nlev_lsm
          write (unitdiag, '(a, i4, e20.12)')  &
                              'in donner_deep: k,cuql', &
                              k,Don_conv%cuql (idiag,jdiag    ,k)
          if (ABS(temperature_forcing(idiag,jdiag,k)) > 0.002) then
            write (unitdiag, '(a, i4, e20.12)')  &
                             'in donner_deep: k, cemetf= ',k,   &
                                 temperature_forcing(idiag,jdiag,k)
            write (unitdiag, '(a, i4, e20.12)')  &
                              'in donner_deep: k, cual= ',k,    &
                                    Don_conv%cual(idiag,jdiag,k )
            write (unitdiag, '(a, i4, e20.12)')  &
                            'in donner_deep: k, xcape= ',k,    &
                                  Don_cape%xcape_lag(idiag,jdiag) 
            write (unitdiag, '(a, i4, e20.12)')   &
                              'in donner_deep: k, dcape = ',k,    &
                                  Don_conv%dcape(idiag,jdiag)
            write (unitdiag, '(a, i4, e20.12)')  &
                              'in donner_deep: k,a1    = ',k,    &
                                  Don_conv%a1 (idiag,jdiag)
            write (unitdiag, '(a, i4, e20.12)')   &
                              'in donner_deep: k, amax  = ',k,   &
                                  Don_conv%amax(idiag,jdiag)
          endif
        end do
      endif   ! (not exit_flag)

!--------------------------------------------------------------------



end subroutine don_d_output_cupar_diags_k



!####################################################################

subroutine don_d_dealloc_loc_vars_k   &
         (Don_conv, Don_cape, Don_rad, ermesg)

!----------------------------------------------------------------------
!    subroutine don_d_dealloc_loc_vars_k deallocates the
!    local variables found in subroutine donner_deep of donner_deep_mod.
!    these are limited to the pointer components of the donner_conv_type,
!    donnr_cape_type and donner_rad_type arrays resident there.
!----------------------------------------------------------------------

use messy_convect_donner_types_mod, only : donner_conv_type, donner_cape_type, &
                             donner_rad_type

implicit none

!----------------------------------------------------------------------
type(donner_conv_type),         intent(inout) :: Don_conv
type(donner_cape_type),         intent(inout) :: Don_cape
type(donner_rad_type),          intent(inout) :: Don_rad 
character(len=*),               intent(out)   :: ermesg

!----------------------------------------------------------------------
!   intent(inout) variables:
!
!     Don_conv             donner_convection_type derived type variable
!                          containing diagnostics and intermediate
!                          results describing the nature of the convec-
!                          tion produced by the donner parameterization
!     Don_cape             donner_cape type derived type variable con-
!                          taining diagnostics and intermediate results
!                          related to the cape calculation associated
!                          with the donner convection parameterization
!     Don_rad              donner_rad_type derived type variable used
!                          to hold those fields needed to connect the
!                          donner deep convection parameterization and
!                          the model radiation package
!
!  intent(out) variables:
! 
!     ermesg               character string containing any error message
!                          to be returned to calling routine
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the error message string.
!---------------------------------------------------------------------
      ermesg = ' '

!----------------------------------------------------------------------
!    deallocate the components of the donner_conv_type variable.
!----------------------------------------------------------------------
      deallocate (Don_conv%conv_temp_forcing    )
      deallocate (Don_conv%conv_moist_forcing   )
      deallocate (Don_conv%ceefc                )
      deallocate (Don_conv%cecon                )
      deallocate (Don_conv%cemfc                )
      deallocate (Don_conv%cememf_mod           )
      deallocate (Don_conv%cual                 )
      deallocate (Don_conv%fre                  )
      deallocate (Don_conv%elt                  )
      deallocate (Don_conv%cmus                 )
      deallocate (Don_conv%ecds                 )
      deallocate (Don_conv%eces                 )
      deallocate (Don_conv%emds                 )
      deallocate (Don_conv%emes                 )
      deallocate (Don_conv%mrmes                )
      deallocate (Don_conv%wmps                 )
      deallocate (Don_conv%wmms                 )
      deallocate (Don_conv%tmes                 )
      deallocate (Don_conv%dmeml                )
      deallocate (Don_conv%uceml                )
      deallocate (Don_conv%detmfl               )
      deallocate (Don_conv%umeml                )
      deallocate (Don_conv%xice                 ) 
      deallocate (Don_conv%xliq                 )
      deallocate (Don_conv%qtren1               )
      deallocate (Don_conv%qtceme               )
      deallocate (Don_conv%qtmes1               )
      deallocate (Don_conv%wtp1                 )
      deallocate (Don_conv%dgeice               )
      deallocate (Don_conv%cuqi                 )
      deallocate (Don_conv%cuql                 )
      deallocate (Don_conv%cell_liquid_eff_diam )
      deallocate (Don_conv%cell_ice_geneff_diam )
      deallocate (Don_conv%dcape                )  
      deallocate (Don_conv%a1                   )
      deallocate (Don_conv%amax                 )
      deallocate (Don_conv%amos                 )
      deallocate (Don_conv%ampta1               )
      deallocate (Don_conv%cell_precip          )
      deallocate (Don_conv%meso_precip          )
      deallocate (Don_conv%emdi_v               )
      deallocate (Don_conv%prztm                )
      deallocate (Don_conv%przm                 )
      deallocate (Don_conv%pb_v                 )
      deallocate (Don_conv%pmd_v                )
      deallocate (Don_conv%pztm_v               )
      deallocate (Don_conv%pzm_v                )

!----------------------------------------------------------------------
!    deallocate the components of the donner_cape_type variable.
!----------------------------------------------------------------------
      deallocate (Don_cape%coin       )
      deallocate (Don_cape%plcl       )
      deallocate (Don_cape%plfc       )
      deallocate (Don_cape%plzb       )
      deallocate (Don_cape%xcape      )
      deallocate (Don_cape%xcape_lag  )
      deallocate (Don_cape%parcel_r   )
      deallocate (Don_cape%parcel_t   )
      deallocate (Don_cape%cape_p     )
      deallocate (Don_cape%env_r      )
      deallocate (Don_cape%env_t      )
      deallocate (Don_cape%model_p    )
      deallocate (Don_cape%model_r    )
      deallocate (Don_cape%model_t    )
      deallocate (Don_cape%qint       )     
      deallocate (Don_cape%qint_lag   ) 

!----------------------------------------------------------------------
!    deallocate the components of the donner_rad_type variable.
!----------------------------------------------------------------------
      deallocate (Don_rad%cell_cloud_frac  )
      deallocate (Don_rad%cell_liquid_amt  )
      deallocate (Don_rad%cell_liquid_size )
      deallocate (Don_rad%cell_ice_amt     )
      deallocate (Don_rad%cell_ice_size    )
      deallocate (Don_rad%meso_cloud_frac  )
      deallocate (Don_rad%meso_liquid_amt  )
      deallocate (Don_rad%meso_liquid_size )
      deallocate (Don_rad%meso_ice_amt     )
      deallocate (Don_rad%meso_ice_size    )
      deallocate (Don_rad%nsum             )        

!----------------------------------------------------------------------


end subroutine don_d_dealloc_loc_vars_k 



!######################################################################
END MODULE MESSY_CONVECT_DONNER_DEEP_K
