
MODULE messy_clamstracer_global

  USE messy_clams_global,               ONLY: PREC, filenamelen
  USE messy_clamstracer_regional_masks, ONLY: eradata_type

  INTEGER     :: ntracer=0              ! number of tracers

  type :: tracer_type
     character(10) :: name 
     character(10) :: parent_specie_name
     character(10) :: pulse_date               ! pulse startdate (string)
     real(PREC)    :: pulse_time               ! pulse startdate in julian seconds
     integer       :: pulse_length = 0         ! pulse duration in days
     integer       :: pulse_reset_period = 0   ! pulse reset after n years
     real(PREC)    :: pulse_time_chk = 0.      !
     integer       :: overwrite = 0.           !
     integer       :: type_pulsing             ! to determine if pulsing is repeated according to date or according to number of days (continuous pulsing)
     real(PREC)    :: latmin = 0.              ! latitude range (type='grid')
     real(PREC)    :: latmax = 0.              ! 
     real(PREC)    :: lonmin = 0.              ! longitude range (type='grid')
     real(PREC)    :: lonmax = 0.              ! 
     real(PREC)    :: levmin = 0.              ! level range (type='grid')
     real(PREC)    :: levmax = 0.              ! 
     integer       :: seg_no = 0               ! number of map segment (type='map')
     real(PREC)    :: const  = 0.
     real(PREC)    :: lin    = 0.
     real(PREC)    :: quad   = 0.
     real(PREC)    :: amplitude = 0.           ! amplitude of the cycle
     real(PREC)    :: periode   = 0.           ! periode
     real(PREC)    :: phase     = 0.           ! phase of the cycle
     real(PREC)    :: tau    = 0.
     real(PREC), dimension(:), pointer :: values 
  end type tracer_type


  !---------------------------------------------------------------------------
  ! Namelist variables
  !---------------------------------------------------------------------------
  
  INTEGER         :: timestep_tracer=0 ! timestep in hours
  INTEGER         :: step_perp=0       ! perpetuum-step (transient run: 0, else: step in perpetuum loop)

  CHARACTER(20)   :: tracertype=""          ! tracer type: grid/map/region
  CHARACTER(160)  :: tracer_desc_file=""

  ! used for tracertype="region":
  integer         :: set_to_zero = 0 ! set to zero in the level outside of the geographic region where pulsing is done

  REAL(PREC)      :: lev_bound=0.  ! level boundary    
  CHARACTER(filenamelen):: file_lsm=""   ! netcdf-file containing land-sea mask 
  CHARACTER(filenamelen):: file_gph=""   ! netcdf-file containing orographie
  type (eradata_type)   :: lsm     ! land-sea mask 
  type (eradata_type)   :: gph     ! gph (on surface) 


End MODULE messy_clamstracer_global
