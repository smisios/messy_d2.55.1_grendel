MODULE mo_diagnosis

  USE mo_kind, ONLY: i4, sp, dp, wp
  USE mo_param1, ONLY: ie, ie_g, je, je_g, ke, kep, nbox
  USE mo_planetary_constants, ONLY: g, rocp, rhoref_ice,rhosnic,rhoref_water
  USE mo_commo1, ONLY : area, weto, tho, sao, sictho, ddpo, sicsno, sicomo, &
       zero, lwetol1_g, wetol1_g, lweto, preff, dt, ndtday, almzer, &
       ddue, dduo, dlxp, dlyp, dlyu, dlxv, dzw, tiestu, &
       stabio, uko, wo, tiestw, alat, alat_g, alatpsi_g, alonpsi_g, fswr, dz, &
       lyears, ldays, lmonts, lmont1, lmont2, ldtdayc, &
       eminpo, sicuo, sicve, tauwatu, tauwatv , txo, tye, &
       vke, avo, dvo, amsue, amsuo, lbounds_exch_tp, rhoo, sicsno, sictho, zo, &
       lwith_barotropic_stokes_drift

  USE mo_commoau1
  USE mo_commoau2
  USE mo_mean, ONLY : tmepo,tmcdo,tmceo
  USE mo_boundsexch, ONLY : bounds_exch

  USE mo_parallel, ONLY: gather, have_g_js, scatter, &
       global_sum, p_ioff, p_joff, p_io, p_pe
  USE mo_basin_masks, ONLY: ibek, ibek_g
  USE mo_grid, ONLY: p_suchij, thkcello
  USE mo_units, ONLY: io_stdout, io_ou_difi
  USE mo_grid, ONLY: get_level_index_by_depth

  USE mo_fluxes1,ONLY :aofltxio, aofltyie
  USE mo_levitus, ONLY : tlevi,slevi

  IMPLICIT NONE
  LOGICAL :: lcalc_mixed_layer_thickness=.FALSE.
  LOGICAL :: lcalc_mixed_layer_thickness_sqr=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: zmld(:,:) !< mixed_layer_thickness
  REAL(wp), ALLOCATABLE,TARGET :: zmld_sqr(:,:) !< square_of_mixed_layer_thickness
  REAL(wp), ALLOCATABLE,TARGET :: amld(:,:) !< maximum_mixed_layer_thickness

  LOGICAL :: lcalc_zo_sqr=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: zo_sqr(:,:) !< square of sea surface height

  LOGICAL :: lcalc_rhopoto=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: rhopoto(:,:,:) !< sea_water_potential_density

  LOGICAL :: lcalc_gmsl_diag=.FALSE.
  REAL(wp), TARGET :: gmsl_st !< global_mean_steric_sealevel
  REAL(wp), TARGET :: gmsl_eu !< global_mean_eustatic_sealevel
  REAL(wp) :: reference_volume = 0.0_wp
  REAL(wp) :: reference_area=0.0_wp
  REAL(wp) :: reference_density=0.0_wp

  LOGICAL :: lcalc_sst=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: sst(:,:) !< sea surface temperature

  LOGICAL :: lcalc_sst_sqr=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: sst_sqr(:,:) !< square of sea surface temperature

  LOGICAL :: lcalc_sss=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: sss(:,:) !< sea surface salinity

  LOGICAL :: lcalc_bottom_pressure=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: bottom_pressure(:,:) !< mass_per_unit_area

  LOGICAL :: lcalc_upward_mass_transport=.FALSE.
  LOGICAL :: lcalc_upward_mass_transport_sqr=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: wmo(:,:,:) !< upward_ocean_mass_transport
  REAL(wp), ALLOCATABLE,TARGET :: wmosq(:,:,:) !<square_of_upward_ocean_mass_transport

  LOGICAL :: lcalc_ocean_mass_transport=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: umo(:,:,:) !< ocean_mass_x_transport
  REAL(wp), ALLOCATABLE,TARGET :: vmo(:,:,:) !< ocean_mass_y_transport

  REAL(wp), ALLOCATABLE,TARGET :: secmap(:,:) !< map of the sections

  REAL(wp), ALLOCATABLE,TARGET :: flum(:,:) !< surface_net_downward_heat_flux_where_sea
  REAL(wp), ALLOCATABLE,TARGET :: pem(:,:) !< water_flux_into_ocean

  LOGICAL :: lcalc_sictr=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: sictru(:,:) !< sea_ice_x_transport
  REAL(wp), ALLOCATABLE,TARGET :: sictrv(:,:) !< sea_ice_y_transport

  LOGICAL :: lcalc_ice_transport=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: transix(:,:) !< sea_ice_x_transport
  REAL(wp), ALLOCATABLE,TARGET :: transiy(:,:) !< sea_ice_y_transport

  LOGICAL :: lcalc_ice_stress=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: strairx(:,:) !< zonal_wind_stress_on_ice
  REAL(wp), ALLOCATABLE,TARGET :: strairy(:,:) !< meridional_wind_stress_on_ice
  REAL(wp), ALLOCATABLE,TARGET :: strocx(:,:) !< zonal_ocean_stress_on_ice
  REAL(wp), ALLOCATABLE,TARGET :: strocy(:,:) !< meridional_ocean_stress_on_ice

  LOGICAL :: lcalc_ice_velocities=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: usi(:,:) !< zonal_sea_ice_velity
  REAL(wp), ALLOCATABLE,TARGET :: vsi(:,:) !< meridional_sea_ice_velocity

  REAL(wp), ALLOCATABLE,TARGET :: ice_mask_u(:,:) !< sea_ice_thickness > 1e-3
  REAL(wp), ALLOCATABLE,TARGET :: ice_mask_v(:,:) !< sea_ice_thickness > 1e-3

  LOGICAL :: lcalc_global_mean=.FALSE.
  REAL(wp), TARGET :: global_volume=0.0_dp
  REAL(wp), TARGET :: global_mass=0.0_dp
  REAL(wp), TARGET :: global_mean_temperature=0.0_dp
  REAL(wp), TARGET :: global_mean_salinity=0.0_dp
  REAL(wp), TARGET :: global_sum_temperature=0.0_dp
  REAL(wp), TARGET :: global_sum_salinity=0.0_dp
  REAL(wp), TARGET :: global_salt_content=0.0_dp
  !FIXME: why not -HUGE?
  REAL(wp), TARGET :: const_subsurface_volume=-9e33_wp

  LOGICAL :: lcalc_region=.FALSE.
  LOGICAL :: lcalc_region_kinetic_energy=.FALSE.
  LOGICAL :: lcalc_region_salinity=.FALSE.
  LOGICAL :: lcalc_region_temperature=.FALSE.
  LOGICAL :: lcalc_region_seaice=.FALSE.
  LOGICAL :: lcalc_region_fluxes=.FALSE.



  TYPE ocean_region

    CHARACTER(256) :: RegionName

    REAL(wp) :: eisab
    REAL(wp) :: eiscb

    REAL(wp) :: hflb
    REAL(wp) :: wflb

    REAL(wp) :: area_0m
    REAL(wp) :: salt_0m
    REAL(wp) :: tem_0m
    REAL(wp) :: u2pv2_0m

    REAL(wp) :: area_200m
    REAL(wp) :: salt_200m
    REAL(wp) :: tem_200m
    REAL(wp) :: u2pv2_200m

    REAL(wp) :: area_800m
    REAL(wp) :: salt_800m
    REAL(wp) :: tem_800m
    REAL(wp) :: u2pv2_800m

    REAL(wp) :: area_2000m
    REAL(wp) :: salt_2000m
    REAL(wp) :: tem_2000m
    REAL(wp) :: u2pv2_2000m

  END TYPE ocean_region

  TYPE ocean_section
    INTEGER :: SecStart(2), SecEnd(2), SecLev(2),SecLay1(2),Seclay2(2)
    CHARACTER(256) :: SecName
    CHARACTER(2) :: SecType

    LOGICAL :: register_netheat=.FALSE.
    LOGICAL :: register_netsalt=.FALSE.
    LOGICAL :: register_netwater=.FALSE.
    LOGICAL :: register_sice=.FALSE.
    LOGICAL :: register_layer1=.FALSE.
    LOGICAL :: register_layer2=.FALSE.
    REAL(wp) ,ALLOCATABLE :: HeatTransport(:)
    REAL(wp) ,ALLOCATABLE :: SaltTransport(:)
    REAL(wp) ,ALLOCATABLE :: WaterTransport(:)
    REAL(wp) :: layer1Transport
    REAL(wp) :: layer2transport
    REAL(wp) :: SiceTransport
    REAL(wp) :: netheatTransport
    REAL(wp) :: netsaltTransport
    REAL(wp) :: netwaterTransport

  END TYPE ocean_section

  TYPE(ocean_section),SAVE,TARGET :: section(21)

  TYPE(ocean_region),SAVE,TARGET :: region(9)

  INTEGER :: noffset

  INTEGER,PARAMETER ::     barents_opening              =1
  INTEGER,PARAMETER ::     bering_strait                =2
  INTEGER,PARAMETER ::     canadian_archipelago         =3
  INTEGER,PARAMETER ::     denmark_strait               =4
  INTEGER,PARAMETER ::     drake_passage                =5

  INTEGER,PARAMETER ::     english_channel              =6
  INTEGER,PARAMETER ::     equatorial_undercurrent      =7
  INTEGER,PARAMETER ::     faroe_scotland_channel       =8
  INTEGER,PARAMETER ::     florida_bahamas              =9
  INTEGER,PARAMETER ::     fram_strait                  =10

  INTEGER,PARAMETER ::     iceland_faroe_channel        =11
  INTEGER,PARAMETER ::     indonesian_throughflow       =12
  INTEGER,PARAMETER ::     mozambiqe_channel            =13
  INTEGER,PARAMETER ::     taiwan_and_luzon_straits     =14
  INTEGER,PARAMETER ::     windward_passage             =15
  INTEGER,PARAMETER ::     strait_of_gibraltar          =16
  INTEGER,PARAMETER ::     atlantic_60n                 =17
  INTEGER,PARAMETER ::     atlantic_40n                 =18
  INTEGER,PARAMETER ::     atlantic_30n                 =19
  INTEGER,PARAMETER ::     atlantic_26n                 =20
  INTEGER,PARAMETER ::     atlantic_30s                 =21


  ! regions in the BEK file

  INTEGER,PARAMETER ::  greenland_iceland_norwegian_sea  =1
  INTEGER,PARAMETER ::  arctic_ocean                    =2
  INTEGER,PARAMETER ::  labrador_sea                    =3
  INTEGER,PARAMETER ::  north_atlantic_ocean            =4
  INTEGER,PARAMETER ::  atlantic_ocean                  =5
  INTEGER,PARAMETER ::  southern_ocean                  =6
  INTEGER,PARAMETER ::  indopacific_ocean               =7
  INTEGER,PARAMETER ::  tropical_pacific_ocean          =8
  INTEGER,PARAMETER ::  global_ocean                    =9  ! (sum of 1 to 9 )


  LOGICAL :: lcalc_psi=.FALSE.
  REAL(wp),ALLOCATABLE,TARGET:: psitro(:,:)


  LOGICAL :: lcalc_moc=.FALSE.
  REAL(wp), ALLOCATABLE,TARGET :: atlantic_moc(:,:)
  REAL(wp), ALLOCATABLE,TARGET :: global_moc(:,:)
  REAL(wp), ALLOCATABLE,TARGET :: indopacific_moc(:,:)

  REAL(wp), ALLOCATABLE,TARGET :: atlantic_hfl(:)
  REAL(wp), ALLOCATABLE,TARGET :: global_hfl(:)
  REAL(wp), ALLOCATABLE,TARGET :: indopacific_hfl(:)

  REAL(wp), ALLOCATABLE,TARGET :: atlantic_wfl(:)
  REAL(wp), ALLOCATABLE,TARGET :: global_wfl(:)
  REAL(wp), ALLOCATABLE,TARGET :: indopacific_wfl(:)

  REAL(wp), ALLOCATABLE,TARGET :: atlantic_zonal_mask(:,:)
  REAL(wp), ALLOCATABLE,TARGET :: global_zonal_mask(:,:)
  REAL(wp), ALLOCATABLE,TARGET :: indopacific_zonal_mask(:,:)


  INTEGER :: ltq1,ltq2,ltq3,ltq4,knadw,kaabw
  INTEGER :: ldenmar,lfaroe, lequ, itsdiag

CONTAINS

  SUBROUTINE alloc_mem_diag


    ALLOCATE(zmld(ie,je))
    ALLOCATE(zmld_sqr(ie,je))
    ALLOCATE(amld(ie,je))
    ALLOCATE(secmap(ie,je))
    ALLOCATE(flum(ie,je),pem(ie,je))

    ALLOCATE(ice_mask_u(ie,je),ice_mask_v(ie,je))
    ALLOCATE(sictru(ie,je),sictrv(ie,je))

    ALLOCATE(transix(ie,je),transiy(ie,je))
    ALLOCATE(usi(ie,je),vsi(ie,je))

    ALLOCATE(global_moc(180,kep))
    ALLOCATE(atlantic_moc(180,kep))
    ALLOCATE(indopacific_moc(180,kep))

    ALLOCATE(global_hfl(180))
    ALLOCATE(atlantic_hfl(180))
    ALLOCATE(indopacific_hfl(180))

    ALLOCATE(global_wfl(180))
    ALLOCATE(atlantic_wfl(180))
    ALLOCATE(indopacific_wfl(180))

    ALLOCATE(global_zonal_mask(180,kep))
    ALLOCATE(atlantic_zonal_mask(180,kep))
    ALLOCATE(indopacific_zonal_mask(180,kep))

    ALLOCATE(psitro(ie,je))

    ALLOCATE(zo_sqr(ie,je))

    ALLOCATE(rhopoto(ie,je,ke))

    ALLOCATE(sst(ie,je))
    ALLOCATE(sst_sqr(ie,je))
    ALLOCATE(sss(ie,je))

    ALLOCATE(bottom_pressure(ie,je))

    ALLOCATE(wmo(ie,je,kep))
    ALLOCATE(wmosq(ie,je,kep))
    ALLOCATE(umo(ie,je,ke))
    ALLOCATE(vmo(ie,je,ke))

    ALLOCATE(strairx(ie,je)) !< zonal_wind_stress_on_ice
    ALLOCATE(strairy(ie,je)) !< meridional_wind_stress_on_ice
    ALLOCATE(strocx(ie,je)) !< zonal_ocean_stress_on_ice
    ALLOCATE(strocy(ie,je)) !< meridional_ocean_stress_on_ice


  END SUBROUTINE alloc_mem_diag

  SUBROUTINE diag_zero_init

    INTEGER :: n

    !  define the region name of the BEK file
    region(greenland_iceland_norwegian_sea)%regionName='greenland_iceland_norwegian_sea'
    region(arctic_ocean)%regionName='arctic_ocean'
    region(labrador_sea)%regionName='labrador_sea'
    region(north_atlantic_ocean)%regionName='north_atlantic_ocean'
    region(atlantic_ocean)%regionName='atlantic_ocean'
    region(atlantic_ocean)%regionName='atlantic_ocean'
    region(southern_ocean)%regionName='southern_ocean'
    region(indopacific_ocean)%regionName='indopacific_ocean'
    region(tropical_pacific_ocean)%regionName='tropical_pacific_ocean'
    region(global_ocean)%regionName='global_ocean'




    DO n=1,SIZE(region)

      region(n)%wflb=0.0_dp
      region(n)%hflb=0.0_dp
      region(n)%eiscb=0.0_dp
      region(n)%eisab=0.0_dp

      region(n)%area_0m=0.0_dp
      region(n)%salt_0m=0.0_dp
      region(n)%tem_0m=0.0_dp
      region(n)%u2pv2_0m=0.0_dp

      region(n)%area_200m=0.0_dp
      region(n)%salt_200m=0.0_dp
      region(n)%tem_200m=0.0_dp
      region(n)%u2pv2_200m=0.0_dp

      region(n)%area_800m=0.0_dp
      region(n)%salt_800m=0.0_dp
      region(n)%tem_800m=0.0_dp
      region(n)%u2pv2_800m=0.0_dp

      region(n)%area_2000m=0.0_dp
      region(n)%salt_2000m=0.0_dp
      region(n)%tem_2000m=0.0_dp
      region(n)%u2pv2_2000m=0.0_dp

    ENDDO


    amld(:,:) = zero
    zmld_sqr(:,:) = zero
    psitro(:,:) = zero
    zo_sqr(:,:) = zero
    sst(:,:) = zero
    sst_sqr(:,:) = zero
    sss(:,:) = zero
    sictru(:,:) = zero
    sictrv(:,:) = zero
    transix(:,:) = zero
    transiy(:,:) = zero
    flum(:,:) = zero
    pem(:,:) = zero
    umo = 0.0_dp
    vmo = 0.0_dp
    ice_mask_u(:,:) = zero
    ice_mask_v(:,:) = zero



    DO n=1,SIZE(section)
      ALLOCATE(Section(n)%WaterTransport(ke))
      ALLOCATE(Section(n)%HeatTransport(ke))
      ALLOCATE(Section(n)%SaltTransport(ke))

      Section(n)%WaterTransport(:)=zero
      Section(n)%HeatTransport(:)=zero
      Section(n)%SaltTransport(:)=zero

      Section(n)%SiceTransport=zero
      Section(n)%Layer1Transport=zero
      Section(n)%layer2transport=zero

      Section(n)%NetHeatTransport=zero
      Section(n)%NetSaltTransport=zero
      Section(n)%NetWaterTransport=zero
    ENDDO


  END SUBROUTINE diag_zero_init

  SUBROUTINE diag_ini


    INTEGER :: i1,j1,i2,j2,k,jb
    REAL(wp) dist

    ! initialisations
    !:: common diagnostics

    ! initialise the section map with the land sea mask
    secmap(:,:) = 0._wp



    ! Compute level index for given depths (in m below NN)
    ltq1 = 1 ! surface: equivalent to get_index_level_by_depth(0._dp)
    ltq2 = get_level_index_by_depth(200._dp)
    ltq3 = get_level_index_by_depth(800._dp)
    ltq4 = get_level_index_by_depth(2000._dp)
    ldenmar = get_level_index_by_depth(300._dp)
    lequ = get_level_index_by_depth(350._dp)
    lfaroe = get_level_index_by_depth(400._dp)
    knadw = get_level_index_by_depth(1000._dp)
    kaabw = get_level_index_by_depth(3000._dp)

    CALL ini_area_averages

     ! search for cmip5 sections

     ! Barents opening = (16.8'E, 76.5'N) to (19.2'E, 70.2'N).
     CALL p_suchij(76.5_wp, 16.8_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(70.2_wp, 19.2_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, 'uv', &
          'barents_opening', barents_opening)


     ! Bering Strait = (171'W, 66.2'N) to (166'W, 65'N).
     CALL p_suchij(66.2_wp, -171.0_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(65.0_wp, -166.0_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, &
          'uv', 'bering_strait', bering_strait)

     ! Canadian Archipelego = (128.2'W, 70.6'N) to (59.3'W, 82.1'N).
     CALL p_suchij(70.6_wp, -128.2_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(82.1_wp, -59.3_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, 'uv', &
          'canadian_archipelago', canadian_archipelago)


     ! Denmark Strait = (37'W, 66.1'N) to (22.5'W, 66'N).
     CALL p_suchij(66.1_wp, -37.0_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(66.0_wp, -22.5_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 1, ldenmar, ldenmar, ke, &
          'uv', 'denmark_strait', denmark_strait)


     ! Drake Passage = (68'W, 54'S) to (60'W, 64.7'S).
     CALL p_suchij(-54.0_wp, -68._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(-64.7_wp, -60._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, &
          'uv', 'drake_passage', drake_passage)

     ! English Channel = (1.5'E, 51.1'N to (1.7'E, 51.0'N).
     CALL p_suchij(51.1_wp, 1.5_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(51.0_wp, 1.7_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, 'uv', &
          'english_channel', english_channel)


     ! Equatorial Undercurrent = (155'W, 3'S) to (155'W, 3'N) over the depth range 0-350m.
     CALL p_suchij(-3.0_wp, 155.0_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij( 3.0_wp, 155.0_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, lequ, 0, 0, 0, 0, 'u', &
          'equatorial_undercurrent', equatorial_undercurrent)


     ! Faroe-Scotland Channel = (6.9'W, 62'N) to (5'W, 58.7'N)
     CALL p_suchij(62.0_wp, 6.9_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(58.7_wp, 5.0_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 1, lfaroe, lfaroe, ke, &
          'uv', 'faroe_scotland_channel', faroe_scotland_channel)


     ! Florida-Bahamas Strait = (78.5'W, 26'N) to (80.5'W, 27'N).
     CALL p_suchij(26._wp, -78.5_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(27._wp, -80.5_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, &
          'uv', 'florida_bahamas', florida_bahamas)


     ! Fram Strait = (11.5'W, 81.3'N to (10.5'E, 79.6'N).
     CALL p_suchij(81.3_wp, -11.5_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(79.6_wp,  10.5_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, 'uv', &
          'fram_strait', fram_strait)

     ! Iceland-Faroe Channel = (13.6'W, 64.9'N) to (7.4'W, 62.2'N)
     CALL p_suchij(64.9_wp, -13.6_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(62.2_wp, -7.4_wp, 1, i2, j2, dist, 1._wp)
     CALL define_sectIon(i1, j1, i2, j2, 1, ke, 1, lfaroe, lfaroe, ke, 'uv', &
          'Iceland_faroe_channel', iceland_faroe_channel)


     ! Indonesian Throughflow = (100'E, 6'S) to (140'E, 6'S).
     CALL p_suchij(-6._wp, 100._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(-6._wp, 140._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, 'uv', &
          'indonesian_throughflow', indonesian_throughflow)


     ! Mozambique Channel = (39'E, 16'S) to (45'E, 18'S).
     CALL p_suchij(-16._wp, 39._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(-18._wp, 45._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, 'uv', &
          'mozambiqe_channel', mozambiqe_channel)


     ! Luzon Strait = (121'E, 15'N) to (121'E, 22'N).
     CALL p_suchij(15._wp, 121._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(22._wp, 121._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, 'uv', &
          'taiwan_and_luzon_straits', taiwan_and_luzon_straits)


     ! Windward Passage = (75'W, 20.2'N) to (72.6'W, 19.7'N).
     CALL p_suchij(20.2_wp, -75.0_wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(19.7_wp, -72.6_wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 0, 0, 0, 0, 'uv', &
          'windward_passage', windward_passage)


     ! Strait of Gibraltar = (6'W, 35.5'N) to (6'W, 36.5'N).
     CALL p_suchij(35.6_wp, -6._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(36.5_wp, -6._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, 1, ltq2, ltq2, ke, 'uv', &
          'strait_of_gibraltar', strait_of_gibraltar)


     ! Atlantic 60N = (70'W, 60'N) to (10'E, 60'N).
     CALL p_suchij(60._wp, -70._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(60._wp,  10._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, knadw, kaabw, kaabw, ke, 'uv', &
          'atlantic_60n', atlantic_60n)


     ! Atlantic 40N = (80'W, 40'N) to (0'E, 40'N).
     CALL p_suchij(40._wp, -80._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(40._wp,   0._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, knadw, kaabw, kaabw, ke, 'uv', &
          'atlantic_40n', atlantic_40n)


     ! Atlantic 30N = (90'W, 30'N) to (0'E, 30'N).
     CALL p_suchij(30._wp, -90._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(30._wp,   0._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, knadw, kaabw, kaabw, ke, 'uv', &
          'atlantic_30n', atlantic_30n)


     ! Atlantic 30N = (90'W, 26'N) to (0'E, 26'N).
     CALL p_suchij(26._wp, -90._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(26._wp,   0._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, knadw, kaabw, kaabw, ke, 'uv', &
          'atlantic_30n', atlantic_26n)


     ! Atlantic 30S = (60'W, 30'S) to (20'E, 30'S).
     CALL p_suchij(-30._wp, -60._wp, 1, i1, j1, dist, 1._wp)
     CALL p_suchij(-30._wp,  20._wp, 1, i2, j2, dist, 1._wp)
     CALL define_section(i1, j1, i2, j2, 1, ke, knadw, kaabw, kaabw, ke, 'uv', &
          'atlantic_30s', atlantic_30s)



!!$    IF (ie_g.EQ.182 .AND. je_g.EQ.84) THEN           ! versiongin
!!$
!!$      CALL define_section(123, 22, 123, 55, 1, ke, 'Atlantic 60N', 4)
!!$      CALL define_section(134, 20, 134, 65, 1, ke, 'Atlantic 40N', 5)
!!$      CALL define_section(148, 35, 148, 72, 1, ke, 'Atlantic 30N', 6)
!!$      CALL define_section(172, 35, 172, 72, 1, ke, 'Atlantic 30S', 7)
!!$      CALL define_section(92, 33, 92, 45, 1, ke, 'Fram Strait', 8)
!!$      CALL define_section(115, 40, 115, 50, 1, ke, 'Denmark Strait', 9)
!!$      CALL define_section(121, 33, 124, 33, 1, ke, 'Faroer Bank Channel', 10)
!!$      CALL define_section(172, 35, 172, 72, 1, ke, 'Strait of Gibraltar', 11)
!!$
!!$    ENDIF
!!$
     IF (ie_g.EQ.66 .AND. je_g.EQ.36) THEN         ! versiontoy


       CALL define_section(50, 4, 52, 7, 1, ke, 0, 0, 0, 0, 'uv', &
            'barents_opening', barents_opening)
       CALL define_section(17, 7, 19, 7, 1, ke, 0, 0, 0, 0, 'uv', &
            'bering_strait', bering_strait)
       CALL define_section(42, 3, 45, 3, 1, ke, 0, 0, 0, 0, 'v', &
            'canadian_archipelago', canadian_archipelago)
       CALL define_section(48, 7, 48, 7, 1, ke, 1, ldenmar, ldenmar, ke, 'uv', &
            'denmark_strait', denmark_strait)
       CALL define_section(38, 31, 39, 34, 1, ke, 0, 0, 0, 0, 'uv', &
            'drake_passage', drake_passage)

       CALL define_section(50, 9, 50, 9, 1, ke, 0, 0, 0, 0, 'v', &
            'english_channel', english_channel)
       CALL define_section(14, 18, 14, 20, 1, lequ, 0, 0, 0, 0, 'u', &
            'equatorial_undercurrent', equatorial_undercurrent)
       CALL define_section(49, 8, 49, 8, 1, ke, 1, lfaroe, lfaroe, ke, 'uv', &
            'faroe_scotland_channel', faroe_scotland_channel)
       CALL define_section(37, 12, 38, 12, 1, ke, 0, 0, 0, 0, 'v', &
            'florida_bahamas', florida_bahamas)
       CALL define_section(48, 4, 50, 4, 1, ke, 0, 0, 0, 0, 'uv', &
            'fram_strait', fram_strait)

       CALL define_section(49, 8, 49, 8, 1, ke, 1, lfaroe, lfaroe, ke, 'uv', &
            'iceland_faroe_channel', iceland_faroe_channel)
       CALL define_section(5, 19, 10, 22, 1, ke, 0, 0, 0, 0, 'uv', &
            'indonesian_throughflow', indonesian_throughflow)
       CALL define_section(57, 23, 59, 23, 1, ke, 0, 0, 0, 0, 'v', &
            'mozambiqe_channel', mozambiqe_channel)
       CALL define_section(8, 14, 8, 16, 1, ke, 0, 0, 0, 0, 'u', &
            'taiwan_and_luzon_straits', taiwan_and_luzon_straits)
       CALL define_section(38, 14, 38, 15, 1, ke, 0, 0, 0, 0, 'u', &
            'windward_passage', windward_passage)

       CALL define_section(51, 11, 51, 13, 1, ke, 1, ltq2, ltq2, ke, 'uv', &
            'strait_of_gibraltar', strait_of_gibraltar)
       CALL define_section(42, 7, 52, 7, 1, ke, knadw, kaabw, kaabw, ke, 'uv', &
            'atlantic_60n', atlantic_60n)
       CALL define_section(40, 9, 50, 9, 1, ke, knadw, kaabw, kaabw, ke, 'uv', &
            'atlantic_40n', atlantic_40n)
       CALL define_section(33, 13, 50, 13, 1, ke, knadw, kaabw, kaabw, ke, &
            'uv', 'atlantic_30n', atlantic_30n)
       CALL define_section(33, 13, 50, 13, 1, ke, knadw, kaabw, kaabw, ke, &
            'uv', 'atlantic_26n', atlantic_26n)
       CALL define_section(41, 25, 54, 25, 1, ke, knadw, kaabw, kaabw, ke, &
            'uv', 'atlantic_30s', atlantic_30s)

    ENDIF

    IF (ie_g.EQ.362 .AND. je_g.EQ.192) THEN         ! versiontp10

       CALL define_section(281, 15, 284, 23, 1, ke, 0, 0, 0, 0, 'uv', &
            'barents_opening', barents_opening)
       CALL define_section(100, 28, 103, 28, 1, ke, 0, 0, 0, 0, 'v', &
            'bering_strait', bering_strait)
       CALL define_section(246, 10, 253, 10, 1, ke, 0, 0, 0, 0, 'v', &
            'canadian_archipelago', canadian_archipelago)
       CALL define_section(259, 25, 262, 28, 1, ke, 1, ldenmar, ldenmar, ke, &
            'uv', 'denmark_strait', denmark_strait)
       CALL define_section(211, 157, 220, 166, 1, ke, 0, 0, 0, 0, 'uv', &
            'drake_passage', drake_passage)

       CALL define_section(277, 43, 277, 43, 1, ke, 0, 0, 0, 0, 'v', &
            'english_channel', english_channel)
       CALL define_section(73, 90, 73, 96, 1, lequ, 0, 0, 0, 0, 'u', &
            'equatorial_undercurrent', equatorial_undercurrent)
       CALL define_section(271, 33, 273, 36, 1, ke, 1, lfaroe, lfaroe, ke, &
            'uv', 'faroe_scotland_channel', faroe_scotland_channel)
       CALL define_section(198, 60, 201, 60, 1, ke, 0, 0, 0, 0, 'v', &
            'florida_bahamas', florida_bahamas)
       CALL define_section(269, 13, 278, 13, 1, ke, 0, 0, 0, 0, 'v', &
            'fram_strait', fram_strait)

       CALL define_section(267, 30, 271, 32, 1, ke, 1, lfaroe, lfaroe, ke, &
            'uv', 'iceland_faroe_channel', iceland_faroe_channel)
       CALL define_section(32, 102, 61, 107, 1, ke, 0, 0, 0, 0, 'uv', &
            'indonesian_throughflow', indonesian_throughflow)
       CALL define_section(317, 112, 323, 112, 1, ke, 0, 0, 0, 0, 'v', &
            'mozambiqe_channel', mozambiqe_channel)
       CALL define_section(40, 66, 40, 72, 1, ke, 0, 0, 0, 0, 'u', &
            'taiwan_and_luzon_straits', taiwan_and_luzon_straits)
       CALL define_section(203, 68, 207, 70, 1, ke, 0, 0, 0, 0, 'uv', &
            'windward_passage', windward_passage)

       CALL define_section(273, 56, 273, 58, 1, ke, 1, ltq2, ltq2, ke, 'u', &
            'strait_of_gibraltar', strait_of_gibraltar)
       CALL define_section(232, 30, 283, 30, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_60n', atlantic_60n)
       CALL define_section(207, 43, 278, 43, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_40n', atlantic_40n)
       CALL define_section(198, 54, 270, 54, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_30n', atlantic_30n)
       CALL define_section(198, 60, 272, 60, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_26n', atlantic_26n)
       CALL define_section(224, 133, 298, 133, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_30s', atlantic_30s)


    ENDIF

    IF (ie_g.EQ.802 .AND. je_g.EQ.404) THEN         ! versiontp04

       CALL define_section(623, 30, 629, 48, 1, ke, 0, 0, 0, 0, 'uv', &
            'barents_opening', barents_opening)
       CALL define_section(222, 58, 227, 58, 1, ke, 0, 0, 0, 0, 'v', &
            'bering_strait', bering_strait)
       CALL define_section(546, 19, 562, 19, 1, ke, 0, 0, 0, 0, 'v', &
            'canadian_archipelago', canadian_archipelago)
       CALL define_section(567, 55, 582, 58, 1, ke, 1, ldenmar, ldenmar, ke, &
            'uv', 'denmark_strait', denmark_strait)
       CALL define_section(467, 346, 484, 367, 1, ke, 0, 0, 0, 0, 'uv', &
            'drake_passage', drake_passage)

       CALL define_section(616, 92, 616, 92, 1, ke, 0, 0, 0, 0, 'v', &
            'english_channel', english_channel)
       CALL define_section(161, 195, 161, 210, 1, lequ, 0, 0, 0, 0, 'u', &
            'equatorial_undercurrent', equatorial_undercurrent)
       CALL define_section(602, 69, 602, 77, 1, ke, 1, lfaroe, lfaroe, ke, &
            'uv', 'faroe_scotland_channel', faroe_scotland_channel)
       CALL define_section(439, 128, 445, 128, 1, ke, 0, 0, 0, 0, 'v', &
            'florida_bahamas', florida_bahamas)
       CALL define_section(595, 28, 613, 28, 1, ke, 0, 0, 0, 0, 'v', &
            'fram_strait', fram_strait)

       CALL define_section(592, 62, 602, 68, 1, ke, 1, lfaroe, lfaroe, ke, &
            'uv', 'iceland_faroe_channel', iceland_faroe_channel)
       CALL define_section(52, 218, 133, 230, 1, ke, 0, 0, 0, 0, 'uv', &
            'indonesian_throughflow', indonesian_throughflow)
       CALL define_section(702, 246, 716, 246, 1, ke, 0, 0, 0, 0, 'v', &
            'mozambiqe_channel', mozambiqe_channel)
       CALL define_section(87, 145, 87, 156, 1, ke, 0, 0, 0, 0, 'u', &
            'taiwan_and_luzon_straits', taiwan_and_luzon_straits)
       CALL define_section(452, 148, 456, 150, 1, ke, 0, 0, 0, 0, 'uv', &
            'windward_passage', windward_passage)

       CALL define_section(604, 123, 604, 125, 1, ke, 1, ltq2, ltq2, ke, &
            'u', 'strait_of_gibraltar', strait_of_gibraltar)
       CALL define_section(513, 60, 625, 60, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_60n', atlantic_60n)
       CALL define_section(458, 94, 616, 94, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_40n', atlantic_40n)
       CALL define_section(439, 114, 598, 114, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_30n', atlantic_30n)
       CALL define_section(440, 128, 601, 128, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_26n', atlantic_26n)
       CALL define_section(496, 291, 660, 291, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_30s', atlantic_30s)

    ENDIF

    IF (ie_g.EQ.3602 .AND. je_g.EQ.2394) THEN         ! versiontp6m

      CALL define_section(2777, 139, 2860, 185, 1, ke, 0, 0, 0, 0, 'uv', &
           'barents_opening', barents_opening)
      CALL define_section(981, 267, 1059, 267, 1, ke, 0, 0, 0, 0, 'v', &
           'bering_strait', bering_strait)
      CALL define_section(2434, 88, 2530, 88, 1, ke, 0, 0, 0, 0, 'v', &
           'canadian_archipelago', canadian_archipelago)
      CALL define_section(2603, 219, 2603, 259, 1, ke, 1, ldenmar, ldenmar, &
           ke, 'uv', 'denmark_strait', denmark_strait)
      CALL define_section(2095, 1663, 2210, 1858, 1, ke, 0, 0, 0, 0, 'uv', &
           'drake_passage', drake_passage)

      CALL define_section(2763, 401, 2771, 411, 1, ke, 0, 0, 0, 0, 'v', &
           'english_channel', english_channel)
      CALL define_section(1221, 867, 1221, 939, 1, lequ, 0, 0, 0, 0, 'u', &
           'equatorial_undercurrent', equatorial_undercurrent)
      CALL define_section(2701, 305, 2715, 339, 1, ke, 1, lfaroe, lfaroe, ke, &
           'uv', 'faroe_scotland_channel', faroe_scotland_channel)
      CALL define_section(1975, 567, 1999, 567, 1, ke, 0, 0, 0, 0, 'v', &
           'florida_bahamas', florida_bahamas)
      CALL define_section(2661, 119, 2757, 119, 1, ke, 0, 0, 0, 0, 'v', &
           'fram_strait', fram_strait)

      CALL define_section(2659, 273, 2701, 305, 1, ke, 1, lfaroe, lfaroe, ke, &
           'uv', 'iceland_faroe_channel', iceland_faroe_channel)


      CALL define_section(225, 966, 594, 1036, 1, ke, 0, 0, 0, 0, 'uv', &
           'indonesian_throughflow', indonesian_throughflow)
      CALL define_section(3163, 1090, 3223, 1090, 1, ke, 0, 0, 0, 0, 'v', &
           'mozambiqe_channel', mozambiqe_channel)
      CALL define_section(383, 640, 383, 695, 1, ke, 0, 0, 0, 0, 'u', &
           'taiwan_and_luzon_straits', taiwan_and_luzon_straits)
      CALL define_section(2025, 653, 2047, 667, 1, ke, 0, 0, 0, 0, 'uv', &
           'windward_passage', windward_passage)

      CALL define_section(2713, 545, 2713, 550, 1, ke, 1, ltq2, ltq2, ke, &
           'u', 'strait_of_gibraltar', strait_of_gibraltar)
      CALL define_section(2307, 282, 2804, 282, 1, ke, knadw, kaabw, kaabw, &
           ke, 'v', 'atlantic_60n', atlantic_60n)
      CALL define_section(2034, 434, 2726, 434, 1, ke, knadw, kaabw, kaabw, &
           ke, 'v', 'atlantic_40n', atlantic_40n)
      CALL define_section(1966, 503, 2694, 503, 1, ke, knadw, kaabw, kaabw, &
           ke, 'v', 'atlantic_30n', atlantic_30n)
      CALL define_section(1972, 573, 2692, 573, 1, ke, knadw, kaabw, kaabw, &
           ke, 'v', 'atlantic_26n', atlantic_26n)
      CALL define_section(2223, 1317, 2961, 1317, 1, ke, knadw, kaabw, kaabw, &
           ke, 'v', 'atlantic_30s', atlantic_30s)

    ENDIF
!!$
!!$
!!$    IF(ie_g.EQ.60.AND.je_g.EQ.50)THEN         ! versiongr60
!!$
!!$      CALL define_section(21, 14, 36, 14, 1, ke, 'Atlantic 60N', 4)
!!$      CALL define_section(19, 17, 36, 17, 1, ke, 'Atlantic 40N', 5)
!!$      CALL define_section(17, 21, 37, 21, 1, ke, 'Atlantic 30N', 6)
!!$      CALL define_section(22, 34, 36, 34, 1, ke, 'Atlantic 30S', 7)
!!$      CALL define_section(52, 3, 52, 11, 1, ke, 'Fram Strait', 8)
!!$      CALL define_section(34, 1, 34, 6, 1, ke, 'Denmark Strait', 9)
!!$      CALL define_section(36, 8, 36, 14, 1, ke, 'Faroer Bank Channel', 10)
!!$      CALL define_section(34, 21, 34, 21, 1, ke, 'Strait of Gibraltar', 11)
!!$
!!$    ENDIF
!!$
!!$
    IF (ie_g.EQ.122 .AND. je_g.EQ.101) THEN         ! versiongr30

       CALL define_section(92, 29, 99, 23, 1, ke, 0, 0, 0, 0, 'uv', &
            'barents_opening', barents_opening)
       CALL define_section(1, 39, 3, 39, 1, ke, 0, 0, 0, 0, 'uv', &
            'bering_strait', bering_strait)
       CALL define_section(11, 19, 16, 19, 1, ke, 0, 0, 0, 0, 'uv', &
            'canadian_archipelago', canadian_archipelago)
       CALL define_section(70, 1, 70, 12, 1, ke, 1, ldenmar, ldenmar, ke, &
            'uv', 'denmark_strait', denmark_strait)
       CALL define_section(41, 79, 41, 88, 1, ke, 0, 0, 0, 0, 'uv', &
            'drake_passage', drake_passage)

       CALL define_section(75, 34, 75, 36, 1, ke, 0, 0, 0, 0, 'u', &
            'english_channel', english_channel)
       CALL define_section(119, 65, 119, 65, 1, lequ, 0, 0, 0, 0, 'u', &
            'equatorial_undercurrent', equatorial_undercurrent)
       CALL define_section(75, 28, 77, 25, 1, ke, 1, lfaroe, lfaroe, ke, 'uv', &
            'faroe_scotland_channel', faroe_scotland_channel)
       CALL define_section(36, 47, 37, 47, 1, ke, 0, 0, 0, 0, 'v', &
            'florida_bahamas', florida_bahamas)
       CALL define_section(101, 3, 101, 23, 1, ke, 0, 0, 0, 0, 'uv', &
            'fram_strait', fram_strait)

       CALL define_section(76, 18, 77, 24, 1, ke, 1, lfaroe, lfaroe, ke, 'uv', &
            'iceland_faroe_channel', iceland_faroe_channel)
       CALL define_section(108, 66, 113, 68, 1, ke, 0, 0, 0, 0, 'uv', &
            'indonesian_throughflow', indonesian_throughflow)
       CALL define_section(80, 68, 85, 68, 1, ke, 0, 0, 0, 0, 'uv', &
            'mozambiqe_channel', mozambiqe_channel)
       CALL define_section(109, 57, 109, 60, 1, ke, 0, 0, 0, 0, 'u', &
            'taiwan_and_luzon_straits', taiwan_and_luzon_straits)
       CALL define_section(39, 50, 39, 51, 1, ke, 0, 0, 0, 0, 'u', &
            'windward_passage', windward_passage)

       CALL define_section(69, 43, 69, 43, 1, ke, 1, ltq2, ltq2, ke, 'uv', &
            'strait_of_gibraltar', strait_of_gibraltar)
       CALL define_section(31, 19, 83, 28, 1, ke, knadw, kaabw, kaabw, ke, &
            'uv', 'atlantic_60n', atlantic_60n)
       CALL define_section(40, 9, 50, 9, 1, ke, knadw, kaabw, kaabw, ke, 'v', &
            'atlantic_40n', atlantic_40n)
       CALL define_section(35, 46, 66, 46, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_30n', atlantic_30n)
       CALL define_section(35, 47, 66, 47, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_26n', atlantic_26n)
       CALL define_section(43, 71, 75, 71, 1, ke, knadw, kaabw, kaabw, ke, &
            'v', 'atlantic_30s', atlantic_30s)

    ENDIF


    IF (ie_g.EQ.256 .AND. je_g.EQ.220) THEN         ! versiongr15

       CALL define_section(194, 58, 210, 52, 1, ke, 0, 0, 0, 0, 'uv', &
            'barents_opening', barents_opening)
       CALL define_section(6, 86, 10, 86, 1, ke, 0, 0, 0, 0, 'v', &
            'bering_strait', bering_strait)
       CALL define_section(14, 31, 41, 31, 1, ke, 0, 0, 0, 0, 'v', &
            'canadian_archipelago', canadian_archipelago)
       CALL define_section(119, 15, 152, 24, 1, ke, 1, ldenmar, ldenmar, ke, &
            'uv', 'denmark_strait', denmark_strait)
       CALL define_section(89, 170, 92, 182, 1, ke, 0, 0, 0, 0, 'uv', &
            'drake_passage', drake_passage)

       CALL define_section(159, 76, 159, 74, 1, ke, 0, 0, 0, 0, 'u', &
            'english_channel', english_channel)
       CALL define_section(252, 133, 252, 137, 1, lequ, 0, 0, 0, 0, 'u', &
            'equatorial_undercurrent', equatorial_undercurrent)
       CALL define_section(160, 60, 162, 52, 1, ke, 1, lfaroe, lfaroe, ke, &
            'uv', 'faroe_scotland_channel', faroe_scotland_channel)
       CALL define_section(75, 101, 78, 101, 1, ke, 0, 0, 0, 0, 'v', &
            'florida_bahamas', florida_bahamas)
       CALL define_section(216, 45, 226, 35, 1, ke, 0, 0, 0, 0, 'uv', &
            'fram_strait', fram_strait)

       CALL define_section(160, 39, 162, 52, 1, ke, 1, lfaroe, lfaroe, ke, &
            'uv', 'iceland_faroe_channel', iceland_faroe_channel)
       CALL define_section(222, 140, 245, 143, 1, ke, 0, 0, 0, 0, 'uv', &
            'indonesian_throughflow', indonesian_throughflow)
       CALL define_section(171, 144, 177, 144, 1, ke, 0, 0, 0, 0, 'v', &
            'mozambiqe_channel', mozambiqe_channel)
       CALL define_section(231, 122, 231, 125, 1, ke, 0, 0, 0, 0, 'u', &
            'taiwan_and_luzon_straits', taiwan_and_luzon_straits)
       CALL define_section(78, 105, 84, 107, 1, ke, 0, 0, 0, 0, 'uv', &
            'windward_passage', windward_passage)

       CALL define_section(146, 91, 146, 93, 1, ke, 1, ltq2, ltq2, ke, 'u', &
            'strait_of_gibraltar', strait_of_gibraltar)
       CALL define_section(71, 41, 177, 60, 1, ke, knadw, kaabw, kaabw, ke, 'uv', &
            'atlantic_60n', atlantic_60n)
       CALL define_section(76, 89, 144, 89, 1, ke, knadw, kaabw, kaabw, ke, 'v', 'atlantic_40n', atlantic_40n)
       CALL define_section(75, 99, 138, 99, 1, ke, knadw, kaabw, kaabw, ke, 'v', 'atlantic_30n', atlantic_30n)
       CALL define_section(75, 101, 136, 101, 1, ke, knadw, kaabw, kaabw, ke, 'v', 'atlantic_26n', atlantic_26n)
       CALL define_section(92, 155, 155, 155, 1, ke, knadw, kaabw, kaabw, ke, 'v', 'atlantic_30s', atlantic_30s)

    ENDIF

!#endif /*diag*/



    CALL ini_gmsl_diag




  END SUBROUTINE diag_ini


  SUBROUTINE diagnosis
    INTEGER :: n


!    CALL calc_icecutoff  ! this should go into mo_ocice and should be called in main program once per day

    IF ( lcalc_rhopoto ) THEN
      CALL calc_rhopoto
    END IF


    IF ( lcalc_gmsl_diag ) THEN
      CALL calc_gmsl_diag
    END IF

    IF ( lcalc_psi ) THEN
      CALL calc_psi
    END IF

    IF ( lcalc_zo_sqr ) THEN
      CALL calc_square_of_sea_level
    END IF

    IF ( lcalc_sst ) THEN
      CALL calc_sst
    END IF

    IF ( lcalc_sst_sqr ) THEN
      CALL calc_square_of_sst
    END IF

    IF ( lcalc_sss ) THEN
      CALL calc_sss
    END IF

    IF ( lcalc_moc ) THEN
      CALL calc_moc
    END IF

    CALL calc_flux

    IF ( lcalc_ice_transport .OR. lcalc_ice_stress .OR. lcalc_ice_velocities ) THEN
      CALL calc_ice_mask
    ENDIF

    IF ( lcalc_sictr ) THEN
      CALL calc_sictr
    END IF

    IF ( lcalc_ice_transport ) THEN
      CALL calc_ice_transport
    END IF

    IF ( lcalc_ice_stress ) THEN
      CALL calc_ice_stress
    END IF

    IF ( lcalc_ice_velocities ) THEN
      CALL calc_ice_velocities
    END IF


    IF ( lcalc_mixed_layer_thickness ) THEN
      CALL calc_mixed_layer_thickness
    END IF

    IF ( lcalc_mixed_layer_thickness_sqr ) THEN
      CALL calc_mixed_layer_thickness_sqr
    END IF

    DO n=1,SIZE(section)
      CALL calc_section(n)
    END DO

    IF ( lcalc_region ) THEN
      CALL calc_area_averages
    END IF

    IF ( lcalc_global_mean ) THEN
      CALL calc_global_mean
    END IF

    IF ( lcalc_bottom_pressure ) THEN
      CALL calc_bottom_pressure
    ENDIF

    IF ( lcalc_upward_mass_transport ) THEN
      CALL calc_upward_mass_transport
    ENDIF

    IF ( lcalc_ocean_mass_transport ) THEN
      CALL calc_ocean_mass_transport
    ENDIF

  END SUBROUTINE diagnosis



  ! prepare the static masks for moc diagnostics
  SUBROUTINE ini_moc
    USE mo_parallel, ONLY : have_g_js
    USE mo_commo1, ONLY : lbounds_exch_tp

    REAL(wp) :: zlat

    INTEGER :: i, j, k, jbrei, lbrei, l, jb

    ! mask
    atlantic_zonal_mask(:,:) = 0._wp
    global_zonal_mask(:,:) = 0._wp
    indopacific_zonal_mask(:,:) = 0._wp

    jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

    jbrei=3
    DO i=2,ie-1
      DO j=jb,je-1
        IF (weto(i, j, 1) .GT. 0.5_wp) THEN
          !     1   suedpol
          !     180 nordpol
          lbrei = NINT(90._wp + alat(i,j))
          lbrei=MAX(lbrei,1)
          lbrei=MIN(lbrei,180)

          DO k=1,ke   ! masks at level kep are always zero
            zlat = MIN(dlxp(i, j), dlyp(i, j))/(REAL(2*jbrei, wp) * 111111._wp)
            DO l=-jbrei,jbrei
              lbrei = NINT(90._wp + alat(i,j) + REAL(l, wp) * zlat)
              lbrei=MAX(lbrei,1)
              lbrei=MIN(lbrei,180)
              global_zonal_mask(lbrei, k) = global_zonal_mask(lbrei, k) &
                   + weto(i, j, k) / REAL(2*jbrei + 1, wp)
              IF ((ibek(i, j) .LE. atlantic_ocean)) THEN
                atlantic_zonal_mask(lbrei, k) = atlantic_zonal_mask(lbrei, k) &
                     + weto(i, j, k) / REAL(2*jbrei + 1, wp)
              ELSE
                indopacific_zonal_mask(lbrei, k) = &
                     indopacific_zonal_mask(lbrei, k) &
                     + weto(i, j, k) / REAL(2*jbrei + 1, wp)
              END IF
            END DO
          END DO
        END IF
      END DO
    END DO

    CALL GLOBAL_SUM(atlantic_zonal_mask)
    CALL GLOBAL_SUM(global_zonal_mask)
    CALL GLOBAL_SUM(indopacific_zonal_mask)

  END SUBROUTINE ini_moc

  SUBROUTINE calc_moc
    USE mo_parallel, ONLY : have_g_js
    USE mo_planetary_constants, ONLY : rhoref_water


    REAL(wp) :: zlat
    INTEGER :: i,j,k,jbrei,lbrei,l,lb,jb


    atlantic_moc(:,:) = 0._wp
    global_moc(:,:) = 0._wp
    indopacific_moc(:,:) = 0._wp

    atlantic_hfl(:) = 0._wp
    global_hfl(:) = 0._wp
    indopacific_hfl(:) = 0._wp

    atlantic_wfl(:) = 0._wp
    global_wfl(:) = 0._wp
    indopacific_wfl(:) = 0._wp

    jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

    jbrei=3
    DO i=2,ie-1
      DO j=jb,je-1
        IF (weto(i,j,1) .GT. 0.5_wp) THEN
          !     1   suedpol
          !     180 nordpol
          lbrei = NINT(90._wp + alat(i,j))
          lbrei=MAX(lbrei,1)
          lbrei=MIN(lbrei,180)

          DO k=1,kep
            zlat = MIN(dlxp(i, j), dlyp(i, j)) &
                 / (REAL(2*jbrei, wp) * 111111._wp)
            DO l=-jbrei,jbrei
              lbrei = NINT(90._wp + alat(i,j) + REAL(l, wp) * zlat)
              lbrei=MAX(lbrei,1)
              lbrei=MIN(lbrei,180)

              global_moc(lbrei,k)=global_moc(lbrei,k)-area(i,j)*rhoref_water &
                   * wo(i, j, k) / REAL(2*jbrei + 1, wp)

              IF ( k == 1 ) THEN
                global_hfl(lbrei) = global_hfl(lbrei) &
                     - area(i, j) * flum(i, j) / REAL(2*jbrei + 1, wp)
                global_wfl(lbrei) = global_wfl(lbrei) &
                     - area(i, j) * pem(i, j) / REAL(2*jbrei + 1, wp)
              END IF

              IF((ibek(i,j).LE. atlantic_ocean ))THEN
                atlantic_moc(lbrei, k) = atlantic_moc(lbrei, k) &
                     - area(i, j) * rhoref_water &
                     * wo(i, j, k) / REAL(2*jbrei + 1, wp)

                IF ( k == 1 ) THEN
                  atlantic_hfl(lbrei) = atlantic_hfl(lbrei) &
                       - area(i, j) * flum(i, j) / REAL(2*jbrei + 1, wp)
                  atlantic_wfl(lbrei) = atlantic_wfl(lbrei) &
                       - area(i, j) * pem(i, j) / REAL(2*jbrei + 1, wp)
                END IF

              ELSE
                indopacific_moc(lbrei, k) = indopacific_moc(lbrei, k) &
                     - area(i, j) * rhoref_water &
                     * wo(i, j, k) / REAL(2*jbrei+1, wp)
                IF ( k == 1 ) THEN
                  indopacific_hfl(lbrei) = indopacific_hfl(lbrei) &
                       - area(i, j) * flum(i, j) / REAL(2*jbrei + 1, wp)
                  indopacific_wfl(lbrei) = indopacific_wfl(lbrei) &
                       - area(i, j) * pem(i, j) / REAL(2*jbrei + 1, wp)
                END IF

              END IF
            END DO
          END DO
        END IF
      END DO
    END DO

    CALL GLOBAL_SUM(atlantic_moc)
    CALL GLOBAL_SUM(global_moc)
    CALL GLOBAL_SUM(indopacific_moc)

    CALL GLOBAL_SUM(atlantic_hfl)
    CALL GLOBAL_SUM(global_hfl)
    CALL GLOBAL_SUM(indopacific_hfl)

    CALL GLOBAL_SUM(atlantic_wfl)
    CALL GLOBAL_SUM(global_wfl)
    CALL GLOBAL_SUM(indopacific_wfl)

    IF (p_pe==p_io) THEN
      DO lb=179,1,-1
          global_moc(lb,:)=global_moc(lb+1,:)+global_moc(lb,:)
          atlantic_moc(lb,:)=atlantic_moc(lb+1,:)+atlantic_moc(lb,:)
          indopacific_moc(lb,:)=indopacific_moc(lb+1,:)+indopacific_moc(lb,:)

          global_hfl(lb)=global_hfl(lb+1)+global_hfl(lb)
          atlantic_hfl(lb)=atlantic_hfl(lb+1)+atlantic_hfl(lb)
          indopacific_hfl(lb)=indopacific_hfl(lb+1)+indopacific_hfl(lb)

          global_wfl(lb)=global_wfl(lb+1)+global_wfl(lb)
          atlantic_wfl(lb)=atlantic_wfl(lb+1)+atlantic_wfl(lb)
          indopacific_wfl(lb)=indopacific_wfl(lb+1)+indopacific_wfl(lb)

      END DO

      ! FIXME: tj: what is -9e33? almost -inf?
      global_moc(:,:)=MERGE(global_moc(:,:), -9e33_wp, &
           global_zonal_mask(:,:) > 0.5_wp)
      atlantic_moc(:,:)=MERGE(atlantic_moc(:,:), -9e33_wp, &
           atlantic_zonal_mask(:,:) > 0.5_wp)
      indopacific_moc(:,:)=MERGE(indopacific_moc(:,:),-9e33_wp, &
           indopacific_zonal_mask(:,:) > 0.5_wp)

    END IF ! p_pe==p_pio


  END SUBROUTINE calc_moc



SUBROUTINE calc_flux

  ! Heat flux
  FLUM(:,:)=QSWO(:,:)+QLWO(:,:)+QLAO(:,:)+QSEO(:,:)

  ! Fresh water flux
  PEM(:,:)=PRECH(:,:)+EMINPO(:,:)

END SUBROUTINE calc_flux


SUBROUTINE calc_difi(idate)
  !
  !  calculates 2d-fields helpful for diagnostics of ocean fluxes
  !  e.g. sea ice melting from below etc.
  !  snapshots at the end of months
  !  interface:
  !  idate parameter for header of output
  !
  USE mo_planetary_constants, ONLY : rhosnwa,rhoicwa

  REAL(wp) hice(ie,je),tdz(ie,je),sdz(ie,je)
  INTEGER idate
  INTEGER (kind=i4) i41,i42,i43,i44
  REAL(kind=sp) :: ff_4(ie_g,je_g)
  INTEGER i,j,k
  REAL(wp) :: ff_g(ie_g,je_g)
  !
  !  initialize fields

  hice(:,:) = 0._wp
  tdz(:,:) = 0._wp
  sdz(:,:) = 0._wp

  DO k=1,ke
   DO j=1,je
   DO i=1,ie
   tdz(i,j)=tdz(i,j)+ddpo(i,j,k)*tho(i,j,k)*weto(i,j,k)
   sdz(i,j)=sdz(i,j)+ddpo(i,j,k)*sao(i,j,k)*weto(i,j,k)
   ENDDO
  ENDDO
  ENDDO

  DO j=1,je
   DO i=1,ie
    hice(i,j)=(rhosnwa*sicsno(i,j)+rhoicwa*sictho(i,j))*weto(i,j,1)
    sdz(i,j)=sdz(i,j)+(sao(i,j,1)*(zo(i,j)-rhosnwa*sicsno(i,j))       &
            - (sao(i, j, 1) - 5._wp) * rhoicwa * sictho(i, j)) * weto(i, j, 1)
    tdz(i,j)=tdz(i,j)+tho(i,j,1)*(zo(i,j)-rhosnwa*sicsno(i,j)         &
                                -rhoicwa*sictho(i,j))*weto(i,j,1)
   ENDDO
  ENDDO


!      CALL gather_arr(hice(:,:),ff_g,p_io)
      CALL gather(hice(:,:),ff_g,p_io)
      i41=idate
      i42=201
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.EQ.p_io) THEN
       ff_4 = REAL(ff_g, sp)
       WRITE(io_ou_difi)i41,i42,i43,i44
       WRITE(io_ou_difi)ff_4
      ENDIF

!      CALL gather_arr(tdz,ff_g,p_io)
      CALL gather(tdz,ff_g,p_io)
      i41=idate
      i42=202
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.EQ.p_io) THEN
       ff_4 = REAL(ff_g,sp)
       WRITE(io_ou_difi)i41,i42,i43,i44
       WRITE(io_ou_difi)ff_4
      ENDIF

!      CALL gather_arr(sdz,ff_g,p_io)
      CALL gather(sdz,ff_g,p_io)
      i41=idate
      i42=203
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.EQ.p_io) THEN
       ff_4 = REAL(ff_g, sp)
       WRITE(io_ou_difi)i41,i42,i43,i44
       WRITE(io_ou_difi)ff_4
      ENDIF



  END  SUBROUTINE calc_difi


  SUBROUTINE calc_section(n)

    USE mo_planetary_constants, ONLY : rhosnic

    INTEGER,INTENT(in) :: n
    REAL(wp) :: transport(ke),trmin,trmax
    INTEGER :: ci1,ci2,cj1,cj2
    INTEGER :: i,j,k

    IF (section(n)%register_netheat)    section(n)%HeatTransport(:) = 0.0_wp
    IF (section(n)%register_netsalt)    section(n)%SaltTransport(:) = 0.0_wp
    IF (section(n)%register_netwater)  section(n)%WaterTransport(:) = 0.0_wp
    IF (section(n)%register_sice)       section(n)%SiceTransport = 0.0_wp

    DO k=section(n)%SecLev(2),section(n)%SecLev(1),-1

      IF (section(n)%SecType == 'uv' .OR. section(n)%SecType == 'v') THEN

        IF (section(n)%SecStart(1) > section(n)%SecEnd(1)) THEN
          ci1=section(n)%SecEnd(1)
          ci2=section(n)%SecStart(1)
        ELSE
          ci1=section(n)%SecStart(1)
          ci2=section(n)%SecEnd(1)
        ENDIF

        j=section(n)%SecStart(2)-p_joff

        IF ( j > 2  .AND. j < je-1 ) THEN

          DO i=MAX(2,ci1-p_ioff),MIN(ie-1,ci2-p_ioff)

            IF (section(n)%register_netheat) THEN
              section(n)%HeatTransport(k) = section(n)%HeatTransport(k)     &
                   + vke(i, j, k) * dlxv(i, j) * ddue(i, j, k) &
                   * (tho(i,j,k) + tho(i, j+1, k)) * 0.5_wp * rocp
            ENDIF

            IF (section(n)%register_netsalt) THEN
              section(n)%SaltTransport(k) = section(n)%SaltTransport(k)     &
                   + vke(i, j, k) * dlxv(i, j) * ddue(i, j, k) &
                   * (sao(i, j, k) + sao(i, j+1, k)) * 0.5_wp
            ENDIF

            IF (section(n)%register_netwater .OR. section(n)%register_layer2 &
                 .OR.section(n)%register_layer1 ) THEN
              section(n)%WaterTransport(k) = section(n)%WaterTransport(k)                                           &
                   +vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)*rhoref_water
            ENDIF

            IF (section(n)%register_sice .AND. k == section(n)%SecLev(1)) THEN
              section(n)%SiceTransport=section(n)%SiceTransport+((ABS(sicve(i,j))+sicve(i,j))                    &
                   *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
                   +((sicuo(i,j)-ABS(sicve(i,j))))               &
                   *(sictho(i,j+1)+sicsno(i,j+1)*rhosnic))            &
                   * 0.5_wp * dlxv(i, j) * rhoref_ice
            ENDIF
          ENDDO
        ENDIF
      ENDIF


      IF (section(n)%SecType == 'uv' .OR. section(n)%SecType == 'u') THEN

        IF (section(n)%SecStart(2) > section(n)%SecEnd(2)) THEN
          cj1=section(n)%SecEnd(2)
          cj2=section(n)%SecStart(2)
        ELSE
          cj1=section(n)%SecStart(2)
          cj2=section(n)%SecEnd(2)
        ENDIF

        i=section(n)%SecEnd(1)-p_ioff

        IF ( i > 2  .AND. i < ie-1 ) THEN

          DO j=MAX(2,cj1-p_joff),MIN(je-1,cj2-p_joff)

            IF (section(n)%register_netheat) THEN
              section(n)%HeatTransport(k) = section(n)%HeatTransport(k) &
                   + uko(i, j, k) * dlyu(i, j) * dduo(i, j, k) &
                   * (tho(i, j, k) + tho(i+1, j, k)) * 0.5_wp * rocp
            ENDIF

            IF (section(n)%register_netsalt) THEN
              section(n)%SaltTransport(k) = section(n)%SaltTransport(k) &
                   + uko(i, j, k) * dlyu(i, j) * dduo(i, j, k) &
                   * (sao(i, j, k) + sao(i+1, j, k)) * 0.5_wp
            ENDIF

            IF (section(n)%register_netwater .OR. section(n)%register_layer2 &
                 .OR. section(n)%register_layer1) THEN
              section(n)%WaterTransport(k) = section(n)%WaterTransport(k) &
                   + uko(i, j, k) * dlyu(i, j) * dduo(i, j, k) * rhoref_water
            ENDIF

            IF (section(n)%register_sice .AND.  k == section(n)%SecLev(1) ) THEN
             section(n)%SiceTransport=section(n)%SiceTransport+((ABS(sicuo(i,j))+sicuo(i,j))                    &
                   *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
                   +((sicuo(i,j)-ABS(sicuo(i,j))))               &
                   *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
                   * 0.5_wp * dlyu(i, j) * rhoref_ice
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO

    IF (section(n)%register_netheat) THEN
      CALL global_sum(section(n)%HeatTransport)
      section(n)%NetHeatTransport=SUM(section(n)%HeatTransport(:))
    ENDIF

    IF (section(n)%register_netsalt) THEN
      CALL global_sum(section(n)%SaltTransport)
      section(n)%NetSaltTransport=SUM(section(n)%SaltTransport(:))
    ENDIF

    IF (section(n)%register_netwater .OR. section(n)%register_layer2                                         &
         .OR.section(n)%register_layer1 ) THEN
      CALL global_sum(section(n)%WaterTransport)
      section(n)%NetWaterTransport=SUM(section(n)%WaterTransport(:))
    ENDIF

    IF (section(n)%register_sice) THEN
      CALL global_sum(section(n)%SiceTransport)
      section(n)%SiceTransport=section(n)%SiceTransport
    ENDIF

    IF (section(n)%register_layer2 .OR. section(n)%register_layer1) THEN
      transport(ke)=section(n)%WaterTransport(ke)
      DO k=ke-1,1,-1
        transport(k)=transport(k+1)+section(n)%WaterTransport(k)
      ENDDO
    ENDIF

    IF (section(n)%register_layer2) THEN
!      section(n)%layer2Transport=SUM(section(n)%WaterTransport(section(n)%SecLay2(1):section(n)%SecLay2(2)))
      trmin=MINVAL(transport(section(n)%SecLay2(1):section(n)%SecLay2(2)))
      trmax=MAXVAL(transport(section(n)%SecLay2(1):section(n)%SecLay2(2)))
      section(n)%layer2Transport=MERGE(trmin,trmax,ABS(trmin) > ABS(trmax))
    ENDIF

    IF (section(n)%register_layer1) THEN
!      section(n)%layer1Transport=SUM(section(n)%WaterTransport(section(n)%SecLay1(1):section(n)%SecLay1(2)))
      trmin=MINVAL(transport(section(n)%SecLay1(1):section(n)%SecLay1(2)))
      trmax=MAXVAL(transport(section(n)%SecLay1(1):section(n)%SecLay1(2)))
      section(n)%layer1Transport=MERGE(trmin,trmax,ABS(trmin) > ABS(trmax))
    ENDIF

  END SUBROUTINE calc_section


  SUBROUTINE define_section(i1,j1,i2,j2,k1,k2,b1,b2,b3,b4,stp,name,n)

    USE mo_kind
    USE mo_parallel

    CHARACTER(*) :: name,stp
    INTEGER :: n, i1,j1,i2,j2,k1,k2
    INTEGER :: b1,b2,b3,b4
    INTEGER :: ci1,ci2,cj1,cj2,i,j

    Section(n)%register_netheat   = .FALSE.
    Section(n)%register_netsalt   = .FALSE.
    Section(n)%register_netwater = .FALSE.
    Section(n)%register_sice      = .FALSE.

    Section(n)%register_layer1    = .FALSE.
    Section(n)%register_layer2    = .FALSE.

    Section(n)%SecName=TRIM(name)
    Section(n)%SecType=TRIM(stp)
    Section(n)%SecStart=(/i1,j1/)
    Section(n)%SecEnd=(/i2,j2/)
    Section(n)%SecLev=(/k1,k2/)
    Section(n)%SecLay1=(/b1,b2/)
    Section(n)%SecLay2=(/b3,b4/)


    secmap(:,:) = MERGE(0.0_wp, secmap(:,:), secmap(:,:) == REAL(n, wp))

    ! WRITE SECMAP
    IF (section(n)%SecStart(1) > section(n)%SecEnd(1)) THEN
      ci1=section(n)%SecEnd(1)
      ci2=section(n)%SecStart(1)
    ELSE
      ci1=section(n)%SecStart(1)
      ci2=section(n)%SecEnd(1)
    ENDIF

    IF (section(n)%SecType == 'uv' .OR. section(n)%SecType == 'v') THEN
      j=section(n)%SecStart(2)-p_joff
      IF ( j > 2  .AND. j < je-1 ) THEN
        DO i=MAX(2,ci1-p_ioff),MIN(ie-1,ci2-p_ioff)
          secmap(i,j)=REAL(n,dp)
        ENDDO
      ENDIF
    ENDIF

    IF (section(n)%SecStart(2) > section(n)%SecEnd(2)) THEN
      cj1=section(n)%SecEnd(2)
      cj2=section(n)%SecStart(2)
    ELSE
      cj1=section(n)%SecStart(2)
      cj2=section(n)%SecEnd(2)
    ENDIF

    IF (section(n)%SecType == 'uv' .OR. section(n)%SecType == 'u') THEN
      i=section(n)%SecEnd(1)-p_ioff
      IF ( i > 2  .AND. i < ie-1 ) THEN
        DO j=MAX(2,cj1-p_joff),MIN(je-1,cj2-p_joff)
          secmap(i,j)=REAL(n,dp)
        ENDDO
      ENDIF
    ENDIF


  END SUBROUTINE define_section


  SUBROUTINE calc_mixed_layer_thickness

    USE mo_commo1, ONLY : icontro
    REAL(wp) :: zzzuwe,zzzcou,zzz,sigcrit
    REAL(wp) :: sigh(ie,je)
    INTEGER :: i,j,k

    !   compute mixed-layer depth

!$OMP PARALLEL PRIVATE(i,j,k)

    sigcrit = 0.125_wp

    sigh(:,:)=zero
    zmld(:,:)=zero

    !$OMP DO
    DO j=1,je
       DO i=1,ie
          sigh(i,j)=weto(i,j,1)*sigcrit
          zmld(i,j)=weto(i,j,1)*tiestu(1)
       ENDDO
    ENDDO
    !$OMP END DO

    DO k=2,ke
       !$OMP DO
       DO j=1,je
          DO i=1,ie
             IF (weto(i, j, k) * sigh(i, j) .GT. 1.e-6_wp) THEN
                zzz=MIN(sigh(i,j)/(ABS(stabio(i,j,k))+almzer),dz(k))
                sigh(i,j) = MAX(0._wp, sigh(i,j)-zzz*stabio(i,j,k))
                zmld(i,j)=zmld(i,j)+zzz
             ELSE
                sigh(i,j) = 0._wp
             ENDIF
          ENDDO
       ENDDO
       !$OMP END DO
    ENDDO


    amld(:,:)=MAX(amld(:,:),zmld(:,:))


    IF (icontro /= 0 ) THEN
      zzzuwe = SUM(zmld(2:ie-1,2:je-1)*weto(2:ie-1,2:je-1,1))
      zzzcou = SUM(weto(2:ie-1,2:je-1,1))
      !$OMP SINGLE
      CALL global_sum(zzzuwe,zzzcou)
      !$OMP END SINGLE
      zzzcou = zzzcou + almzer

      WRITE(io_stdout,*)'mean mld: ',zzzuwe/zzzcou
    END IF

    !hh      end mixed layer depths
    !$OMP END PARALLEL

  END SUBROUTINE calc_mixed_layer_thickness

  SUBROUTINE calc_icecutoff
    USE mo_planetary_constants, ONLY :rhoicwa ,rhosnwa
#ifdef PBGC
    USE mo_param1_bgc, ONLY: nocetra
    USE mo_carbch, ONLY: ocetra
    INTEGER l
#endif

    REAL(wp) swmin,eismax,volice,areice,volsno,arges,zmi,zma  &
         ,schwell,ddice

    REAL(wp) svold,svnew,draftold,draftnew

    INTEGER i,j,icou,ierr

    swmin = 0._wp
    icou = 0
    eismax = 0._wp
    volice = 0._wp
    areice = 0._wp
    volsno = 0._wp
    arges = 0._wp
    zmi = 0._wp
    zma = 0._wp

    ierr = 0


    !uwe     ensure that effective thickness of uppermost waterlayer
    !        does not fall below a critical value
    !        in that case ice is converted to water
    !        not heat conserving!!!!
    !

    DO j=1,je
       DO i=1,ie
          IF(lweto(i,j,1))THEN
             schwell = 0.9_wp * (dzw(1) + zo(i,j))

             IF(rhoicwa*sictho(i,j)+rhosnwa*sicsno(i,j).GT.schwell)THEN
                !-et
                WRITE(io_stdout,*)'ice cut-off at ',i,j,sictho(i,j),sao(i,j,1)
                ierr=1

                draftold=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa
                svold=sao(i,j,1)*draftold+sice*sictho(i,j)*rhoicwa

                ddice=MIN((sictho(i,j)*rhoicwa+sicsno(i,j)*rhosnwa     &
                     -schwell)/rhoicwa,sictho(i,j))

                sictho(i,j)=sictho(i,j)-ddice

                draftnew=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa

                !uwe       salt conservation!
!                sao(i,j,1)=(ddice*sice+sao(i,j,1)                      &
!                     *(schwell-ddice*rhoicwa))/schwell

                sao(i,j,1) =(ddice*sice+sao(i,j,1)*draftold)/draftnew

                draftnew=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa
                svnew=sao(i,j,1)*draftnew+sice*sictho(i,j)*rhoicwa

                WRITE(io_stdout,*)'ice cut-off Salt Content Change : ', svnew-svold


#ifdef PBGC
                DO l=1,nocetra
                   ocetra(i,j,1,l)=ocetra(i,j,1,l)*draftold/draftnew
                ENDDO
#endif
                WRITE(io_stdout,*)'eisneu,sneu: ',i,j,sictho(i,j),sao(i,j,1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO


  END SUBROUTINE calc_icecutoff

  SUBROUTINE calc_psi

    USE mo_parallel, ONLY : have_g_js
    USE mo_commo1, ONLY : lbounds_exch_tp
    USE mo_planetary_constants, ONLY : rhoref_water

    INTEGER i,j,k,jb

    REAL(wp) :: uint(ie,je)
    REAL(wp),ALLOCATABLE :: uint_g(:,:)

    IF (p_pe == p_io) THEN
      ALLOCATE(uint_g(ie_g,je_g))
    ELSE
      ALLOCATE(uint_g(0,0))
    ENDIF

    uint(:,:) = 0._wp

    IF (p_pe == p_io) THEN
      uint_g(:,:) = 0.0_wp
    ENDIF

    jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

    !   vertical sum of the x transport
    DO k=1,ke
      uint(:,:)=uint(:,:)-uko(:,:,k)*dlyu(:,:)                       &
               *dduo(:,:,k)*rhoref_water
    ENDDO

    IF ( lwith_barotropic_stokes_drift ) THEN
      DO i=1,ie-1
        uint(i,:)=uint(i,:)-uko(i,:,1)*dlyu(i,:)                     &
             * 0.5_wp * (zo(i,:) + zo(i+1,:)) * rhoref_water
      ENDDO
    ENDIF



    !  now gather uint array, the following loop and sums it up in y direction globally

    CALL gather(uint,uint_g,p_io)

    IF(p_pe==p_io) THEN

      DO j=je_g-1,1,-1
        uint_g(:,j)=uint_g(:,j)+uint_g(:,j+1)
      ENDDO

      ! remove a constant land value
      uint_g(:,:)=uint_g(:,:)-uint_g(2,jb)

    ENDIF

    CALL scatter(uint_g,uint,p_io)

    DO j=2,je-1
      DO i=2,ie-1
        psitro(i,j) = 0.25_wp * (uint(i, j) + uint(i, j-1) &
             +uint(i-1,j)+uint(i-1,j-1))
      ENDDO
    ENDDO

    CALL bounds_exch(1,'s',psitro,'calc_psi')

    DEALLOCATE(uint_g)


  END SUBROUTINE calc_psi

  SUBROUTINE calc_potential_energy_release(is)

    INTEGER :: I,J,K,IS
    REAL(wp) :: swi,thick,sh(ie,je),th(ie,je), rh(ie,je)

    if (is.eq.1) then

       ! mixing diagnostics
       ! potential energy

       tmepo(:,:) = 0._wp

       DO j=1,je
          DO k=1,ke
             DO i=1,ie
                sh(i,j)=sao(i,j,k)
                th(i,j)=tho(i,j,k)
             ENDDO

             CALL adisitj(th,sh,preff(k),j)
             CALL rho1j(th,sh,preff(k),rh,j)

             swi = 0._wp
             IF (k .EQ. 1) swi = 1._wp

             DO i=1,ie
                thick = ddpo(i,j,k)+zo(i,j)*swi
                tmepo(i,j)=tmepo(i,j)-weto(i,j,k)*rh(i,j)*thick*g*thick
             ENDDO
          ENDDO
       ENDDO

    endif

    if (is.eq.2) then

       ! mixing diagnostics
       ! potential energy after convection

       tmcdo(:,:) = 0._wp

       DO j=1,je
           DO k=1,ke
             DO i=1,ie
                sh(i,j)=sao(i,j,k)
                th(i,j)=tho(i,j,k)
             ENDDO

             CALL adisitj(th,sh,preff(k),j)
             CALL rho1j(th,sh,preff(k),rh,j)

             swi = 0._wp
             IF (k .EQ. 1) swi = 1._wp

             DO i=1,ie
                thick = ddpo(i,j,k)+zo(i,j)*swi
                tmcdo(i,j)=tmcdo(i,j)-weto(i,j,k)*rh(i,j)*thick*g*thick
             ENDDO
          ENDDO
       ENDDO

       tmcdo(:,:)=(tmepo(:,:)-tmcdo(:,:))/dt

    endif

    if (is.eq.3) then

       ! mixing diagnostics
       ! potential energy after vertical mixing

       tmceo(:,:) = 0._wp

       DO j=1,je
           DO k=1,ke
             DO i=1,ie
                sh(i,j)=sao(i,j,k)
                th(i,j)=tho(i,j,k)
             ENDDO

             CALL adisitj(th,sh,preff(k),j)
             CALL rho1j(th,sh,preff(k),rh,j)

             swi = 0._wp
             IF (k .EQ. 1) swi = 1._wp

             DO i=1,ie
                thick = ddpo(i,j,k)+zo(i,j)*swi
                tmceo(i,j)=tmceo(i,j)-weto(i,j,k)*rh(i,j)*thick*g*thick
             ENDDO
          ENDDO
       ENDDO
       tmceo(:,:)=(tmepo(:,:)-tmceo(:,:))/dt
    endif


  END SUBROUTINE calc_potential_energy_release


  SUBROUTINE calc_square_of_sea_level

    zo_sqr(:,:)=zo(:,:)*zo(:,:)

  END SUBROUTINE calc_square_of_sea_level

  SUBROUTINE calc_mixed_layer_thickness_sqr

    zmld_sqr(:,:)=zmld(:,:)*zmld(:,:)

  END SUBROUTINE calc_mixed_layer_thickness_sqr

  SUBROUTINE calc_sst

    USE mo_commo1, ONLY : tho
    USE mo_constants, ONLY : tmelt

    sst(:,:)=tho(:,:,1)+tmelt

  END SUBROUTINE calc_sst

  SUBROUTINE calc_square_of_sst

    USE mo_commo1, ONLY : tho
    USE mo_constants, ONLY : tmelt

    sst_sqr(:,:)=(tho(:,:,1)+tmelt)*(tho(:,:,1)+tmelt)

  END SUBROUTINE calc_square_of_sst

  SUBROUTINE calc_sss

    USE mo_commo1, ONLY : sao

    sss(:,:)=sao(:,:,1)

!    WRITE(0,*) 'in calc_sss',MAXVAL(sss(:,:)),MINVAL(sss(:,:))

  END SUBROUTINE calc_sss


  SUBROUTINE ini_area_averages
    !   area averages

    INTEGER :: n,i,j

    DO n=1,SIZE(region)
      region(n)%area_0m = 0.0_dp
      region(n)%area_200m = 0.0_dp
      region(n)%area_800m = 0.0_dp
      region(n)%area_2000m = 0.0_dp
    ENDDO

    !only the 4 levels 0, 200, 800, 2000

    DO j=2,je-1
      DO i=2,ie-1
        n=ibek(i,j)

        IF (lweto(i,j,ltq1)) THEN
          region(n)%area_0m = region(n)%area_0m + area(i,j)
        END IF

        IF (lweto(i,j,ltq2)) THEN
          region(n)%area_200m = region(n)%area_200m + area(i,j)
        END IF

        IF (lweto(i,j,ltq3)) THEN
          region(n)%area_800m = region(n)%area_800m + area(i,j)
        END IF

        IF (lweto(i,j,ltq4)) THEN
          region(n)%area_2000m = region(n)%area_2000m + area(i,j)
        END IF

      END DO
    END DO

    DO n=1,global_ocean-1

    region(global_ocean)%area_0m  = region(global_ocean)%area_0m  &
         + region(n)%area_0m
    region(global_ocean)%area_200m  = region(global_ocean)%area_200m  &
         + region(n)%area_200m
    region(global_ocean)%area_800m  = region(global_ocean)%area_800m  &
         + region(n)%area_800m
    region(global_ocean)%area_2000m  = region(global_ocean)%area_2000m  &
         + region(n)%area_2000m

    END DO


    CALL global_sum(region(:)%area_0m)
    CALL global_sum(region(:)%area_200m)
    CALL global_sum(region(:)%area_800m)
    CALL global_sum(region(:)%area_2000m)


!    WRITE(0,*),'in ini_area_averages',p_pe,region(global_ocean)%area_0m



  END SUBROUTINE ini_area_averages

  SUBROUTINE calc_area_averages

    USE mo_planetary_constants, ONLY : rhosnic
    USE mo_parallel, ONLY : have_g_js
    USE mo_commo1

    !   area averages on the 4 levels 0, 200, 800, 2000

    INTEGER :: i,j,n,jb

    jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

    IF ( lcalc_region_temperature ) THEN

      ! initialsize region
      region(:)%tem_0m = 0.0_dp
      region(:)%tem_200m = 0.0_dp
      region(:)%tem_800m = 0.0_dp
      region(:)%tem_2000m = 0.0_dp

      ! sum up the regions (without halos)
      DO i=2,ie-1
        DO j=jb,je-1

          n=INT(ibek(i,j))

          IF ( n /= 0 ) THEN

            region(n)%tem_0m  = region(n)%tem_0m + area(i,j)*tho(i,j,ltq1)*weto(i,j,ltq1)
            region(n)%tem_200m  = region(n)%tem_200m + area(i,j)*tho(i,j,ltq2)*weto(i,j,ltq2)
            region(n)%tem_800m  = region(n)%tem_800m + area(i,j)*tho(i,j,ltq3)*weto(i,j,ltq3)
            region(n)%tem_2000m  = region(n)%tem_2000m + area(i,j)*tho(i,j,ltq4)*weto(i,j,ltq4)

          END IF

        END DO
      END DO

      ! sum up global_ocean from all regions
      region(global_ocean)%tem_0m  =  SUM(region(:)%tem_0m)
      region(global_ocean)%tem_200m  =  SUM(region(:)%tem_200m)
      region(global_ocean)%tem_800m  = SUM(region(:)%tem_800m)
      region(global_ocean)%tem_2000m  = SUM(region(:)%tem_2000m)

      ! global sum from all pe 's
      CALL global_sum(region(:)%tem_0m)
      CALL global_sum(region(:)%tem_200m)
      CALL global_sum(region(:)%tem_800m)
      CALL global_sum(region(:)%tem_2000m)

      ! divide by the area
      region(:)%tem_0m=region(:)%tem_0m/region(:)%area_0m
      region(:)%tem_200m=region(:)%tem_200m/region(:)%area_200m
      region(:)%tem_800m=region(:)%tem_800m/region(:)%area_800m
      region(:)%tem_2000m=region(:)%tem_2000m/region(:)%area_2000m

    END IF


    IF ( lcalc_region_salinity ) THEN

      ! initialsize region
      region(:)%salt_0m = 0.0_dp
      region(:)%salt_200m = 0.0_dp
      region(:)%salt_800m = 0.0_dp
      region(:)%salt_2000m = 0.0_dp

      ! sum up the regions (without halos)
      DO i=2,ie-1
        DO j=jb,je-1

          n=INT(ibek(i,j))

          IF ( n /= 0 ) THEN

            region(n)%salt_0m = region(n)%salt_0m + area(i,j)*sao(i,j,ltq1)*weto(i,j,ltq1)
            region(n)%salt_200m = region(n)%salt_200m + area(i,j)*sao(i,j,ltq2)*weto(i,j,ltq2)
            region(n)%salt_800m = region(n)%salt_800m + area(i,j)*sao(i,j,ltq3)*weto(i,j,ltq3)
            region(n)%salt_2000m = region(n)%salt_2000m + area(i,j)*sao(i,j,ltq4)*weto(i,j,ltq4)

          END IF

        END DO
      END DO

      ! sum up global_ocean from all regions
      region(global_ocean)%salt_0m  = SUM(region(:)%salt_0m)
      region(global_ocean)%salt_200m  = SUM(region(:)%salt_200m)
      region(global_ocean)%salt_800m  = SUM(region(:)%salt_800m)
      region(global_ocean)%salt_2000m  = SUM(region(:)%salt_2000m)

      ! global sum from all pe 's
      CALL global_sum(region(:)%salt_0m)
      CALL global_sum(region(:)%salt_200m)
      CALL global_sum(region(:)%salt_800m)
      CALL global_sum(region(:)%salt_2000m)

      ! divide by the area
      region(:)%salt_0m=region(:)%salt_0m/region(:)%area_0m
      region(:)%salt_200m=region(:)%salt_200m/region(:)%area_200m
      region(:)%salt_800m=region(:)%salt_800m/region(:)%area_800m
      region(:)%salt_2000m=region(:)%salt_2000m/region(:)%area_2000m

    END IF

    IF ( lcalc_region_kinetic_energy ) THEN

      ! initialsize region
      region(:)%u2pv2_0m = 0.0_dp
      region(:)%u2pv2_200m = 0.0_dp
      region(:)%u2pv2_800m = 0.0_dp
      region(:)%u2pv2_2000m = 0.0_dp

      ! sum up the regions (without halos)
      DO i=2,ie-1
        DO j=jb,je-1

          n=INT(ibek(i,j))

          IF ( n /= 0 ) THEN

            region(n)%u2pv2_0m = region(n)%u2pv2_0m  + area(i,j) * 0.125_dp              &
                 *((uko(i-1,j,ltq1)*amsuo(i-1,j,ltq1)+uko(i,j,ltq1)*amsuo(i,j,ltq1))**2   &
                 +(vke(i,j-1,ltq1)*amsue(i,j-1,ltq1)+vke(i,j,ltq1)*amsue(i,j,ltq1))**2)

            region(n)%u2pv2_200m = region(n)%u2pv2_200m + area(i,j) * 0.125_dp            &
                 *((uko(i-1,j,ltq2)*amsuo(i-1,j,ltq2)+uko(i,j,ltq2)*amsuo(i,j,ltq2))**2   &
                 +(vke(i,j-1,ltq2)*amsue(i,j-1,ltq2)+vke(i,j,ltq2)*amsue(i,j,ltq2))**2)

            region(n)%u2pv2_800m = region(n)%u2pv2_800m + area(i,j) * 0.125_dp            &
                 *((uko(i-1,j,ltq3)*amsuo(i-1,j,ltq3)+uko(i,j,ltq3)*amsuo(i,j,ltq3))**2   &
                 +(vke(i,j-1,ltq3)*amsue(i,j-1,ltq3)+vke(i,j,ltq3)*amsue(i,j,ltq3))**2)

            region(n)%u2pv2_2000m = region(n)%u2pv2_2000m + area(i,j) * 0.125_dp          &
                 *((uko(i-1,j,ltq4)*amsuo(i-1,j,ltq4)+uko(i,j,ltq4)*amsuo(i,j,ltq4))**2   &
                 +(vke(i,j-1,ltq4)*amsue(i,j-1,ltq4)+vke(i,j,ltq4)*amsue(i,j,ltq4))**2)

          END IF

        END DO
      END DO

      ! sum up global_ocean from all regions
      region(global_ocean)%u2pv2_0m  = SUM(region(:)%u2pv2_0m)
      region(global_ocean)%u2pv2_200m  = SUM(region(:)%u2pv2_200m)
      region(global_ocean)%u2pv2_800m  = SUM(region(:)%u2pv2_800m)
      region(global_ocean)%u2pv2_2000m  = SUM(region(:)%u2pv2_2000m)

      ! global sum from all pe 's
      CALL global_sum(region(:)%u2pv2_0m)
      CALL global_sum(region(:)%u2pv2_200m)
      CALL global_sum(region(:)%u2pv2_800m)
      CALL global_sum(region(:)%u2pv2_2000m)

      ! divide by the area
      region(:)%u2pv2_0m=region(:)%u2pv2_0m/region(:)%area_0m
      region(:)%u2pv2_200m=region(:)%u2pv2_200m/region(:)%area_200m
      region(:)%u2pv2_800m=region(:)%u2pv2_800m/region(:)%area_800m
      region(:)%u2pv2_2000m=region(:)%u2pv2_2000m/region(:)%area_2000m

    END IF


    IF ( lcalc_region_fluxes ) THEN

      ! initialsize region
      region(:)%hflb = 0.0_dp
      region(:)%wflb = 0.0_dp

      ! sum up the regions (without halos)
      DO i=2,ie-1
        DO j=jb,je-1

          n=INT(ibek(i,j))

          IF ( n /= 0 ) THEN

            region(n)%hflb=region(n)%hflb + area(i,j)*(qswo(i,j)+qlwo(i,j)+qlao(i,j)+qseo(i,j))
            region(n)%wflb=region(n)%wflb + area(i,j)*(prech(i,j)+eminpo(i,j))

          END IF

        END DO
      END DO

      ! sum up global_ocean from all regions
      region(global_ocean)%hflb  = SUM(region(:)%hflb)
      region(global_ocean)%wflb  = SUM(region(:)%wflb)

      ! global sum from all pe 's
      CALL global_sum(region(:)%hflb)
      CALL global_sum(region(:)%wflb)

    END IF

    IF ( lcalc_region_seaice ) THEN

      ! initialsize region
      region(:)%eiscb = 0.0_dp
      region(:)%eisab = 0.0_dp


      ! sum up the regions (without halos)
      DO i=2,ie-1
        DO j=jb,je-1

          n=INT(ibek(i,j))

          IF ( n /= 0 ) THEN

            region(n)%eiscb=region(n)%eiscb + area(i,j)*(sictho(i,j)+(sicsno(i,j)*rhosnic))
            region(n)%eisab=region(n)%eisab + &
                 MERGE(area(i,j), 0._wp, sicomo(i,j) .GE. 0.15_wp)

          END IF

        END DO
      END DO


      ! sum up global_ocean from all regions
      region(global_ocean)%eisab = SUM(region(:)%eisab)
      region(global_ocean)%eiscb = SUM(region(:)%eiscb)

      ! global sum from all pe 's
      CALL global_sum(region(:)%eiscb)
      CALL global_sum(region(:)%eisab)

    END IF



  END SUBROUTINE calc_area_averages


  SUBROUTINE calc_global_mean



    USE mo_kind, ONLY : dp
    USE mo_parallel, ONLY : have_g_js
    USE mo_planetary_constants, ONLY : rhosnwa,rhoicwa
    USE mo_constants, ONLY : tmelt

    INTEGER :: i,j,k,jb
    REAL(wp) :: vol	

    jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

    global_sum_temperature=0.0_dp
    global_sum_salinity=0.0_dp
    global_mass=0.0_dp
    global_volume=0.0_dp
    global_salt_content=0.0_dp


    DO i=2,ie-1
      DO j=jb,je-1
        DO k=1,ke

          vol= area(i,j)*thkcello(i,j,k)
          global_sum_temperature = global_sum_temperature + tho(i,j,k)*vol
          global_sum_salinity = global_sum_salinity + sao(i,j,k)*vol
          global_volume = global_volume + vol
          global_mass = global_mass + rhoo(i,j,k)*area(i,j)*ddpo(i,j,k)

          IF ( k == 1 ) THEN
            global_salt_content = global_sum_salinity &
                 + 5.0_wp * area(i,j) * rhoicwa * sictho(i, j)
            global_mass = global_mass + rhoo(i,j,1)*area(i,j)*zo(i,j)
          END IF

        END DO
      END DO
    END DO

    CALL global_sum(global_sum_temperature,global_sum_salinity,global_volume,global_mass,global_salt_content)

    IF ( const_subsurface_volume < 0.0_dp ) THEN

      jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

      const_subsurface_volume=0.0_dp
      DO k=2,ke
        const_subsurface_volume=const_subsurface_volume+SUM(area(2:ie-1,jb:je-1)*thkcello(2:ie-1,jb:je-1,k))
      END DO
      CALL global_sum(const_subsurface_volume)
    ENDIF



    global_mean_temperature=(global_sum_temperature/(global_volume+const_subsurface_volume))+tmelt
    global_mean_salinity=global_sum_salinity/(global_volume+const_subsurface_volume)


  END SUBROUTINE calc_global_mean

  SUBROUTINE calc_bottom_pressure

    USE mo_commo1, ONLY : rhoo,ddpo,zo
    USE mo_planetary_constants, ONLY : g, slpref

    INTEGER :: k
    REAL(wp), PARAMETER :: Pa2dbar=1.0e-4_wp

    bottom_pressure(:,:)=g*rhoo(:,:,1)*(ddpo(:,:,1)+zo(:,:))+slpref
    DO k=2,ke
      bottom_pressure(:,:)=bottom_pressure(:,:)+g*rhoo(:,:,k)*ddpo(:,:,k)
    END DO
    bottom_pressure=bottom_pressure*Pa2dbar

  END SUBROUTINE calc_bottom_pressure


  SUBROUTINE calc_rhopoto

    USE mo_param1,ONLY : ie,je,ke
    USE mo_commo1, ONLY : tho,sao
    REAL(wp) :: shelp(ie,je),thelp(ie,je),rhelp(ie,je)
    REAL(wp) :: refpress=0.0_wp
    INTEGER :: j,k

    DO j=1,je
      DO k=1,ke
         thelp(:,j)=tho(:,j,k)
         shelp(:,j)=sao(:,j,k)

         CALL adisitj(thelp,shelp,refpress,j)
         CALL rho1j(thelp,shelp,refpress,rhelp,j)

         rhopoto(:,j,k)=rhelp(:,j)
       ENDDO
     ENDDO

   END SUBROUTINE calc_rhopoto


   SUBROUTINE calc_upward_mass_transport

     USE mo_param1,ONLY : kep
     USE mo_commo1, ONLY : wo, area
     USE mo_planetary_constants, ONLY : rhoref_water

     INTEGER :: k

     DO k=1,kep
       wmo(:,:,k)=wo(:,:,k)*area(:,:)*rhoref_water
     ENDDO

     IF (lcalc_upward_mass_transport_sqr) THEN

       wmosq(:,:,:)=wmo(:,:,:)*wmo(:,:,:)

     ENDIF

   END SUBROUTINE calc_upward_mass_transport

   SUBROUTINE calc_ocean_mass_transport

     USE mo_param1, ONLY : ie,je,ke
     USE mo_planetary_constants, ONLY : rhoref_water
     USE mo_commo1, ONLY : uko,vke,dlxp,dlyp,ddue,dduo,zo

     INTEGER :: i,j,k

     DO k=1,ke
       DO i=1,ie-1
         DO j=1,je-1
           IF (k == 1) THEN
             umo(i,j,k)=rhoref_water*uko(i,j,k)*dlyp(i,j)*dduo(i,j,k)  &
                  + 0.5_wp * (zo(i, j) + zo(i+1, j))
             vmo(i,j,k)=rhoref_water*vke(i,j,k)*dlxp(i,j)*ddue(i,j,k)  &
                  + 0.5_wp * (zo(i, j) + zo(i, j+1))
           ELSE
             umo(i,j,k)=rhoref_water*uko(i,j,k)*dlyp(i,j)*dduo(i,j,k)
             vmo(i,j,k)=rhoref_water*vke(i,j,k)*dlxp(i,j)*ddue(i,j,k)
           ENDIF
         ENDDO
       ENDDO
     ENDDO

     CALL bounds_exch(1,'u',umo,'calc_ocean_mass_transport 1')
     CALL bounds_exch(1,'v',vmo,'calc_ocean_mass_transport 2')


   END SUBROUTINE calc_ocean_mass_transport


   SUBROUTINE calc_ice_mask

     INTEGER :: i,j

     DO i=1,ie1
       ice_mask_u(i,:) = MERGE(0.0_wp, 1.0_wp, &
            (0.5_wp * (sictho(i,:) + sictho(i+1,:)) < 1.e-3_wp))
     ENDDO
     CALL bounds_exch(1,'u+',ice_mask_u,'calc_stress_ice')

     DO j=2,je
       ice_mask_v(:,j) = MERGE(0.0_wp, 1.0_wp, &
            (0.5_wp * (sictho(:,j) + sictho(:,j-1)) < 1.e-3_wp))
     ENDDO
     CALL bounds_exch(1,'v+',ice_mask_v,'calc_stress_ice')

   END SUBROUTINE calc_ice_mask

   SUBROUTINE calc_ice_stress

#ifdef __coupled
     strairx(:,:)=aofltxio(:,:)*ice_mask_u(:,:)
     strairy(:,:)=aofltyie(:,:)*ice_mask_v(:,:)
#else
     strairx(:,:)=txo(:,:)*ice_mask_u(:,:)
     strairy(:,:)=tye(:,:)*ice_mask_v(:,:)
#endif

     strocx(:,:)=tauwatu(:,:)*ice_mask_u(:,:)
     strocy(:,:)=tauwatv(:,:)*ice_mask_v(:,:)

   END SUBROUTINE calc_ice_stress

   SUBROUTINE calc_ice_velocities

     usi(:,:)=sicuo(:,:)*ice_mask_u(:,:)
     vsi(:,:)=sicve(:,:)*ice_mask_v(:,:)

   END SUBROUTINE calc_ice_velocities


   SUBROUTINE calc_ice_transport

     INTEGER :: i,j

     DO i=2,ie1
       DO j=2,je1

         ! sea ice transport (x-direction) in kg s-1 (cmip5)
         transix(i,j)=( (ABS(sicuo(i,j))+sicuo(i,j))                    &
              *(sictho(i,j)+sicsno(i,j)*rhosnic)                       &
              +( sicuo(i,j)-ABS(sicuo(i,j)))                           &
              * (sictho(i+1, j) + sicsno(i+1, j) * rhosnic)) * 0.5_wp  &
              *dlyu(i,j)*rhoref_ice*ice_mask_u(i,j)

         ! sea ice transport (y-direction) in kg s-1 (cmip5)
         transiy(i,j)=( (ABS(sicve(i,j))+sicve(i,j))                    &
              *(sictho(i,j+1)+sicsno(i,j+1)*rhosnic)                   &
              +( sicve(i,j)-ABS(sicve(i,j)))                           &
              * (sictho(i, j) + sicsno(i, j) * rhosnic)) * 0.5_wp      &
              *dlxv(i,j)*rhoref_ice*ice_mask_v(i,j)

       ENDDO
     ENDDO

   END SUBROUTINE calc_ice_transport


   SUBROUTINE calc_sictr

     INTEGER :: i,j

     DO i=1,ie1
       DO j=1,je1

         ! sea ice transport (x-direction) in m2 s-1
         sictru(i,j)=( (ABS(sicuo(i,j))+sicuo(i,j))                    &
              *(sictho(i,j)+sicsno(i,j)*rhosnic)                       &
              +( sicuo(i,j)-ABS(sicuo(i,j)))                           &
              *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic) ) * 0.5_wp

         ! sea ice transport (y-direction) in m2 s-1
         sictrv(i,j)=( (ABS(sicve(i,j))+sicve(i,j))                    &
              *(sictho(i,j+1)+sicsno(i,j+1)*rhosnic)                   &
              +( sicve(i,j)-ABS(sicve(i,j)))                           &
              *(sictho(i,j)+sicsno(i,j)*rhosnic) ) * 0.5_wp


       ENDDO
     ENDDO

   END SUBROUTINE calc_sictr

   SUBROUTINE ini_gmsl_diag

     REAL(wp) :: shelp,thelp
     INTEGER :: k,jb
     REAL(wp) :: reference_rho(ke) !< sea_water_reference_density



     jb=MERGE(3,2,lbounds_exch_tp .AND. have_g_js)

     reference_volume = 0.0_wp
     reference_density = 0.0_wp

     DO k=1,ke

       thelp = 0.0_wp
       shelp = 35.0_wp

       CALL adisit1(thelp,shelp,preff(k))

       CALL rho1(thelp,shelp,preff(k),reference_rho(k))

       reference_volume = reference_volume + SUM(ddpo(2:ie-1,jb:je-1,k) &
            *area(2:ie-1,jb:je-1)*weto(2:ie-1,jb:je-1,k))

     ENDDO

     reference_area=SUM(area(2:ie-1,jb:je-1)*weto(2:ie-1,jb:je-1,1))


     CALL global_sum(reference_volume,reference_area)


     DO k=1,ke
       reference_density = reference_density + SUM( reference_rho(k) * (ddpo(2:ie-1,jb:je-1,k) &
            *area(2:ie-1,jb:je-1)*weto(2:ie-1,jb:je-1,k))/reference_volume)
     ENDDO

     CALL global_sum(reference_density)


   END SUBROUTINE ini_gmsl_diag


   SUBROUTINE calc_gmsl_diag

     REAL(kind=wp) :: density,volume
     INTEGER :: k,jb

     jb=MERGE(3,2,lbounds_exch_tp .AND. have_g_js)

     gmsl_st = 0.0_wp      ! global mean sea level (steric)
     gmsl_eu = 0.0_wp      ! global mean sea level: (eustatic)

     volume = 0.0_wp       ! global volume

     volume = SUM(zo(2:ie-1,jb:je-1)                          &
          *area(2:ie-1,jb:je-1)*weto(2:ie-1,jb:je-1,1))

     DO k =1,ke
       volume = volume + SUM(ddpo(2:ie-1,jb:je-1,k)           &
            *area(2:ie-1,jb:je-1)*weto(2:ie-1,jb:je-1,k))
     ENDDO

     CALL global_sum(volume)

     gmsl_eu = (volume-reference_volume)/reference_area



     density = 0.0_wp
     density = SUM(rhoo(2:ie-1,jb:je-1,1)*(zo(2:ie-1,jb:je-1)     &
          *area(2:ie-1,jb:je-1)*weto(2:ie-1,jb:je-1,1))/reference_volume)
     DO k =1,ke
       density = density + SUM(rhoo(2:ie-1,jb:je-1,k)*(ddpo(2:ie-1,jb:je-1,k)    &
            *area(2:ie-1,jb:je-1)*weto(2:ie-1,jb:je-1,k))/reference_volume)
     ENDDO

     CALL global_sum(density)


     gmsl_st = ((reference_density-density)/reference_density)*reference_volume/reference_area

!     IF (p_pe == p_io) WRITE(0,*) 'gmsl', gmsl_st,gmsl_eu


   END SUBROUTINE calc_gmsl_diag




END MODULE mo_diagnosis
