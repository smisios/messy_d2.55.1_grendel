Module messy_clamsbmix_global

  USE messy_clams_global,    ONLY: PREC, DP, species_type, &
                                   filenamelen, specnamelen
  USE messy_clamsmix_global, ONLY: gridtype

  TYPE bmixbound_type
     CHARACTER(specnamelen) :: spec
     CHARACTER(filenamelen) :: file
     INTEGER                :: no
     INTEGER                :: add
     INTEGER                :: startyear=1900, endyear=2100
  END type bmixbound_type

  TYPE clamsbound_type
     CHARACTER(specnamelen) :: spec
     CHARACTER(specnamelen) :: spec_bf
     CHARACTER(filenamelen) :: file
     REAL(PREC)             :: lbound, ubound
     INTEGER                :: lno, uno
     INTEGER                :: startyear, endyear
  END type clamsbound_type

  integer,       parameter  :: maxboundspec = 50

  ! aus Namelist CTRL:
  INTEGER        :: timestep_bmix
  INTEGER        :: interpol_from_init=0
  REAL(PREC)     :: lev_in_down, lev_in_up, &
                    lat_in_down, lat_in_up, &
                    delta_lev, &
                    max_dist=1000.
  CHARACTER(filenamelen) :: file_bounds, dir_boundfiles
  CHARACTER(filenamelen) :: bmix_boundlist=''
  CHARACTER(filenamelen) :: clams_boundlist=''
  LOGICAL        :: switch_EMAC_H2O=.FALSE.
  REAL(DP)       :: EMAC_H2O_z=280.
  

  INTEGER              :: nbmixbounds
  TYPE(bmixbound_type), SAVE :: bmixbound(maxboundspec)

  INTEGER              :: nclamsbounds
  TYPE(clamsbound_type), SAVE:: clamsbound(maxboundspec)


  INTEGER        :: nlevs
  REAL(DP)       :: time_init_value
  REAL(PREC)     :: lev_down, lev_up, lat_down, lat_up
  LOGICAL        :: replace_low, replace_up, replace_north, replace_south

  type(gridtype) :: latlon_grid

  real(prec),         dimension(:,:), pointer   :: lev0


end Module messy_clamsbmix_global
