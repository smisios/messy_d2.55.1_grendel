module messy_clamschem_globalhet

  USE messy_main_constants_mem, ONLY: DP
  implicit none

  ! number of heterogeneous reactions 
  integer,parameter  :: numhet=33


  ! CTRL_HETERO namelist:
!!$    NAMELIST /CTRL_HETERO/ liq_sdist_sigma, densaero_default, aer_h2so4_default, &
!!$                           ciceinit, cnatinit, sat_meltallowed, param_nat_HR, &
!!$                           transform, saturation_criteria, gamma


  real(kind=DP),save :: liq_sdist_sigma, cice, cnat
  real(kind=DP),save :: saturation_criteria(4,4), transform(4,4), gamma(numhet)
  logical,save       :: allice, allnat, param_nat_HR, &
                        sat_meltallowed, cice_from_clim
  integer,save       :: itr0

  ! information on particle densities 
  ! must be preserved between calls to hetero_ken for different trajectories    
  ! (for subroutine hetsolid...)
  !  - note this is analogous to preserving kstate and laststate
  logical,save       :: liquids,satmelting
  real(kind=DP),save :: densnat,densice,denssat,aer_h2so4_default,cnatinit,ciceinit
  real(kind=DP),save :: densaero_default
  

  !  information on the state of the condensed
  !  material which must be stored by hetero.f between 
  !  calls to rufhet for different trajectories. 
  logical,save       :: natcore, mixedpop, liqtest, t_lt_tnat
  integer,save       :: kstate, laststate
  real(kind=DP),save :: told, pressold 
  

  real(kind=DP),allocatable,save :: teold(:),prold(:)
  real(kind=DP),allocatable,save :: densnat_old(:),densice_old(:),denssat_old(:)
  real(kind=DP),allocatable,save :: cnat_old(:),cice_old(:),vliq_save(:)
  integer*1,allocatable,save     :: astate(:),lstate(:)
  logical,allocatable,save       :: log_state(:,:)

end module messy_clamschem_globalhet

