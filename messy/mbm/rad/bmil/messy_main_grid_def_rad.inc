! -*- f90 -*-

  USE messy_main_grid_def_mem_bi 
  USE messy_main_grid_def,      ONLY: submodstr

  IMPLICIT NONE
  SAVE

  ! NAME AND VERSION OF THE BASEMODEL  
  !
  INTEGER, PARAMETER  :: nlev91=91
  INTEGER, PARAMETER  :: nlev48=48
  INTEGER, PARAMETER  :: nlev40=40
  !
  REAL(dp), DIMENSION(nlev91), TARGET :: hyai_91, hybi_91
  REAL(dp), DIMENSION(nlev48), TARGET :: hyai_48, hybi_48
  REAL(dp), DIMENSION(nlev40), TARGET :: hyai_40, hybi_40
  
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: hyam ! hybrid A coefficient at layer midpoints
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: hybm ! hybrid B coefficient at layer midpoints
  REAL(dp), POINTER,     DIMENSION(:) :: hyai => NULL() ! hybrid A coefficient at layer interfaces
  REAL(dp), POINTER,     DIMENSION(:) :: hybi => NULL() ! hybrid B coefficient at layer interfaces

  REAL(dp), ALLOCATABLE, DIMENSION(:) :: longitude
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: latitude

  REAL(dp), ALLOCATABLE, DIMENSION(:) :: cetah !  *half hybrid vertical levels.
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: ceta  !  *full hybrid vertical levels.

  REAL(dp), ALLOCATABLE, DIMENSION(:) :: zh    !  half levels geometric height
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: zf    !  full levels geometric height
  
  REAL(dp), ALLOCATABLE, TARGET, DIMENSION(:) :: ph  !  half pressure levels.
  REAL(dp), ALLOCATABLE, TARGET, DIMENSION(:) :: pf  !  full pressure levels.

  REAL(dp), ALLOCATABLE, DIMENSION(:) :: coslon
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: sinlon
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: philat
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: philon
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: gl_twomu

  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: coslon_2d
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: sinlon_2d
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: coslat_2d
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: sinlat_2d
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: philat_2d
  ! 
  ! - 2D - fields
  !
  REAL(dp), POINTER, DIMENSION(:,:) :: gboxarea_2d=> NULL() ! gridbox area
  ! 
  ! - 3D - fields
  !
  REAL(dp), POINTER, DIMENSION(:,:,:) :: grmass => NULL() ! grid mass in kg
  REAL(dp), POINTER, DIMENSION(:,:,:) :: grvol  => NULL() ! grid volume in m^3
  REAL(dp), POINTER, DIMENSION(:,:,:) :: deltaz => NULL() ! layerthickness in m

CONTAINS
  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_setup

    IMPLICIT NONE

  END SUBROUTINE main_grid_def_setup
  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_initialize

    USE messy_main_blather,       ONLY: start_message, end_message
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_constants_mem, ONLY: dtr
    USE mo_gaussgrid, ONLY : e5_philat   => philat   &
                           , e5_philon   => philon   &
                           , e5_gl_twomu => gl_twomu &
                           , inigau
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_initialize'
    INTEGER  :: status
    INTEGER  :: i
    REAL(dp) :: zetam, zetap

    CALL start_message(submodstr,'INITIALISATION',substr)

    SELECT CASE(nn)
      CASE(21);  nlon = 64
      CASE(30);  nlon = 96
      CASE(42);  nlon = 128
      CASE(63);  nlon = 192
      CASE(85);  nlon = 256
      CASE(106); nlon = 320
      CASE DEFAULT
        CALL error_bi(': wrong maximal meridional wave number selected (NN).'// &
                     ' Use possible choices: 21, 30, 42, 63, 85, 106.', substr)
    END SELECT
    
    nlat = nlon/2
    
    IF (lzonal_mean) THEN
       nlon = 1
    END IF
    
    nlevp1 = nlev + 1
    nvclev = nlevp1 !   number of levels with vertical coefficients.
    nglat  = nlat   ! GP: NUMBER OF LAT (1 if lcolumn=T)
    nglon  = nlon
    
    ALLOCATE( &
             glat(nlat), glon(nlat) &
           , hyam(nlev), hybm(nlev) &
           , longitude(nlon), latitude(nlat) &
           , cetah(nlevp1), ceta(nlev) &
           , zh(nlevp1), zf(nlev) &
           , ph(nlevp1), pf(nlev) &
           , coslon(nlon*2), sinlon(nlon*2) &
           , philat(nlat), philon(nlon), gl_twomu(nlat) &
           , vct(nvclev*2) &
           , coslon_2d(nlon*2,nlat), sinlon_2d(nlon*2,nlat) &
           , coslat_2d(nlon,nlat), sinlat_2d(nlon,nlat), philat_2d(nlon,nlat) &
           )  
    vct = 0._dp
    cetah = 0._dp
    ceta = 0._dp
    
    nproma  = nlon
    ngpblks = nlat
    kproma  = nproma
    npromz  = nproma
    jrow    = 1
    
    glat(:)  = (/ (i, i=1,nlat)    /)
    glon(:)  = 1
    glats(:) = (/    1 /)
    glate(:) = (/ nlat /)
    glons(:) = (/    1 /)
    glone(:) = (/ nlon /)

    ! -
    ! - initialise the gaussian grid
    ! -
    e5_ngl  = nlat   !   number of gaussian latitudes.
    e5_nhgl = nlat/2 !   (number of gaussian latitudes)/2.
    e5_nlon = nlon
    
    CALL inigau
    
    philat   = e5_philat
    philon   = e5_philon
    gl_twomu = e5_gl_twomu
    
    longitude = philon
    latitude  = philat
    
    ! -
    sinlon(:) = (/ SIN(philon * dtr), SIN(philon * dtr) /)
    coslon(:) = (/ COS(philon * dtr), COS(philon * dtr) /)
    DO i = 1, nlat
       coslon_2d(:,i) = coslon
       sinlon_2d(:,i) = sinlon
    END DO
    DO i = 1, nlon
       philat_2d(i,:) = philat
    END DO
    coslat_2d = COS(philat_2d * dtr)
    sinlat_2d = SIN(philat_2d * dtr)
    !
    ! initialise hybrid coefficients
    !
    CALL ini_hybrid
    !
    ! initialise the vertical grid
    !
    IF (lv_echam) THEN
       !
       ! ECHAM type vertical grid  
       !
       ! initialise the pressure levels
       !
       ! cetah *half hybrid vertical levels.
       ! ceta  *full hybrid vertical levels.
       !
       zetam = vct(1)/apzero + vct(nvclev+1)
       cetah(1) = zetam
       !
       DO i = 1, nlev
          zetap      = vct(i+1)/apzero + vct(nvclev+i+1)
          ceta (i)   = (zetam+zetap)*0.5_dp
          cetah(i+1) = zetap
          zetam      = zetap
       END DO

       ! To get the pressure on full and half model levels
       ! pf and ph are multiplied with aps in global_start.

       pf(:) = ceta(:)
       ph(:) = cetah(:)

       PRINT *,'hybi = ',hybi
    ELSE
      ! vertical coordinate system
      ! -----------------------------------------
      ! pressure level heights
      DO i = 1, nlev
         zf(i) = REAL(nlev - i,dp) ! km
      END DO
      !
      ! in ECHAM the pressure at the full levels is defined in the middle
      ! between the upper and lower layer interfaces: 
      ! pf(jk)=(ph(jk)+ph(jk+1))/2 with ph(1)=0. and ph(nlev+1)=psrf
      !
      ! Here the half level pressure levels have to be constructed for 
      ! the given full level pressure grid defined by the equidistant 
      ! z levels: pf(jk)=p0*exp(-zf/H0)
      !
      ! --> The vertical grid differs by construction from that in 
      !     the GCM, with consequences for the results 
      !
      ! pressure at full levels
      DO i = 1, nlev
         pf(i) = apzero * EXP(-zf(i)/h0)
      END DO

      hyam(:) = pf(:)

      ! pressure at interface levels
      ph(1)=0._dp
      ph(nlevp1)=101325._dp
      DO i = 2, nlev
         zh(i) = ( zf(i-1) + zf(i) )/2._dp
         ph(i) = apzero * EXP(-zh(i)/h0)
      END DO
      PRINT *,'zf = ',zf
      PRINT *,'zh = ',zh
      !
      ceta(:)  = pf(:)/apzero
      cetah(:) = ph(:)/apzero

    END IF
    PRINT *,SHAPE(hyai), ', hyai = ',hyai
    PRINT *,SHAPE(hybm), ', hybm = ',hybm
    PRINT *,SHAPE(hyam), ', hyam = ',hyam
    PRINT *,SHAPE(ceta), ', ceta(rev)  = ',( ceta(i), i=nlev,1,-1 )
    PRINT *,SHAPE(cetah),', cetah(rev) = ',( cetah(i), i=nlevp1,1,-1 )
    PRINT *,SHAPE(ceta), ', ceta       = ',ceta
    PRINT *,SHAPE(cetah),', cetah      = ',cetah
    PRINT *,SHAPE(pf),   ', pf(rev) = ',( pf(i), i=nlev,1,-1 )
    PRINT *,SHAPE(ph),   ', ph(rev) = ',( ph(i), i=nlevp1,1,-1 )
    PRINT *,SHAPE(pf),   ', pf      = ',pf
    PRINT *,SHAPE(ph),   ', ph      = ',ph
    PRINT *,SHAPE(pf),   ', pf*apzero = ',pf*apzero
    PRINT *,SHAPE(ph),   ', ph*apzero = ',ph*apzero
    PRINT *,SHAPE(vct),  ', vct=',vct
    !
    PRINT *,SHAPE(philon),', philon = ',philon
    PRINT *,SHAPE(philat),', philat = ',philat
    PRINT *,SHAPE(glon),', glon = ',glon
    PRINT *,SHAPE(glat),', glat = ',glat
    PRINT *,SHAPE(glons),', glons = ',glons
    PRINT *,SHAPE(glats),', glats = ',glats
    PRINT *,SHAPE(glone),', glone = ',glone
    PRINT *,SHAPE(glate),', glate = ',glate
    !
    CALL end_message(submodstr,'INITIALISATION',substr)
    !
  CONTAINS
    ! ==============================================================================
    SUBROUTINE ini_hybrid
      !
      IF (lv_echam) THEN
         SELECT CASE(nlevp1)
         CASE(91)
            hyai_91 = (/   &
              0., 1.98918199539185, 4.3762001991272, 7.24062299728394,                 &
              10.6779298782349, 14.8027000427246, 19.5461807250977, 25.0011901855469,  &
              31.2744407653809, 38.4886894226074, 46.7850685119629, 56.077018737793,   &
              66.484001159668, 77.9316787719727, 90.5241317749023, 104.375801086426,   &
              119.612701416016, 136.373199462891, 154.809799194336, 175.090103149414,  &
              197.398406982422, 221.714401245117, 248.218902587891, 277.108795166016,  & 
              308.598785400391, 342.922912597656, 380.336090087891, 421.116607666016,  &
              465.567199707031, 514.018493652344, 566.830322265625, 624.059387207031,  &
              686.075073242188, 753.277709960938, 826.101196289062, 905.015686035156,  &
              990.530578613281, 1083.19799804688, 1183.61596679688, 1292.43298339844,  &
              1410.35095214844, 1538.13195800781, 1676.60095214844, 1826.65197753906,  &
              1989.251953125, 2165.36889648438, 2356.15795898438, 2562.85498046875,    &
              2786.81396484375, 3029.52490234375, 3292.64697265625, 3579.71997070312,  &
              3893.169921875, 4234.77685546875, 4608.55810546875, 5015.7158203125,     &
              5456.84912109375, 5932.2880859375, 6443.5322265625, 6987.4560546875,     &
              7559.75, 8156.06396484375, 8775.7900390625, 9416.1591796875,             &
              10073.33984375, 10742.2802734375, 11418.490234375, 12094.240234375,      &
              12760.150390625, 13405.150390625, 14016.9404296875, 14584.3095703125,    &
              15090.9599609375, 15521.2802734375, 15855.01953125, 16071.1396484375,    &
              16145.599609375, 16050.3701171875, 15754.1298828125, 15222.0302734375,   &
              14418.08984375, 13300.7998046875, 11833.7197265625, 9980.625,            &
              7717.98486328125, 5125.14208984375, 2549.96899414062, 783.195007324219,  &
              0., 0., 0. &
              /)
            hybi_91 = (/    &
              0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &
              0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &
              0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &
              0., 2.48594496952137e-05, 0.000107582898635883,                         &
              0.000281307293334976, 0.000592263008002192, 0.00110967794898897,        &
              0.0019136150367558, 0.00309228408150375, 0.00475343177095056,           &
              0.00701756495982409, 0.0100215198472142, 0.0139185795560479,            &
              0.0188941806554794, 0.0251423493027687, 0.0328751802444458,             &
              0.0423218309879303, 0.0537220016121864, 0.0673977434635162,             &
              0.0836439281702042, 0.102857403457165, 0.125381097197533,               &
              0.151691406965256, 0.182314097881317, 0.217832893133163,                &
              0.258867502212524, 0.306137204170227, 0.360260188579559,                &
              0.422180503606796, 0.492536813020706, 0.572195112705231,                &
              0.661917328834534, 0.759698390960693, 0.856437623500824,                &
              0.9287468791008, 0.972985208034515, 0.992281496524811, 1.               &
              /)
            hyai => hyai_91
            hybi => hybi_91
         CASE(48)
            hyai_48 = (/   &
              0., 1.98918533325195, 6.57208967208862, 15.6739025115967, &
              30.6242790222168,  54.5457191467285, 92.5588302612305,    &
              150.504699707031, 235.327453613281, 356.100250244141,     &
              523.919494628906, 751.04296875, 1051.13720703125,         &
              1438.98840332031, 1930.17736816406, 2540.69702148438,     &
              3286.55297851562, 4199.57421875, 5303.95703125,           &
              6624.7041015625, 8187.18505859375, 9976.13671875,         &
              11820.5400390625, 13431.3896484375, 14736.3603515625,     &
              15689.2099609375, 16266.6103515625, 16465.,               &
              16297.6201171875, 15791.599609375, 14985.26953125,        &
              13925.51953125, 12665.2900390625, 11261.23046875,         &
              9771.40625, 8253.2109375, 6761.33984375,                  &
              5345.9140625, 4050.71801757812, 2911.56909179688,         &
              1954.80505371094, 1195.89001464844, 638.14892578125,      &
              271.626495361328, 72.0635986328125, 0., 0., 0.            &
              /)
            hybi_48 = (/   &
              0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,  &
              0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      &
              0.00039999998989515, 0.00289999996311963, 0.00920000020414591,  &
              0.0203000009059906, 0.0370000004768372, 0.0595000013709068,     &
              0.0878999978303909, 0.122000001370907, 0.161400005221367,       &
              0.205699995160103, 0.254200011491776, 0.306199997663498,        &
              0.361099988222122, 0.418199986219406, 0.476700007915497,        &
              0.535899996757507, 0.595099985599518, 0.653599977493286,        &
              0.710600018501282, 0.765399992465973, 0.817200005054474,        &
              0.865000009536743, 0.907700002193451, 0.944199979305267,        &
              0.97299998998642, 0.992299973964691, 1.  &
              /)
            hyai => hyai_48
            hybi => hybi_48
         CASE(40)
            hyai_40 = (/   &
              0., 1.989182, 5.09726, 10.51217, 18.75759, 30.39117, 47.1266,  &
              71.56693, 106.5414, 155.6308, 223.2727, 315.8983, 442.4553,    &
              613.6625, 843.0477, 1147.513, 1547.971, 2070.053, 2744.883,    &
              3609.916, 4709.828, 6097.457, 7834.766, 9959.486, 12317.67,    &
              14597.68, 16441.67, 17538.87, 17708.88, 16970.06, 15587.9,     &
              13621.46, 11101.56, 8127.155, 5125.142, 2549.969, 783.195,     &
              0., 0., 0.  &
                 /)
            hybi_40 = (/   &
              0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                &
              0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                    &
              0.0003389933, 0.003357187, 0.01307004, 0.03407715, 0.07064983, &
              0.1259167,    0.2011954, 0.2955196, 0.4054092, 0.5249322,      &
              0.646108, 0.7596984,  0.8564376, 0.9287469, 0.9729852,         &
              0.9922815, 1.   &
              /)
            hyai => hyai_40
            hybi => hybi_40
         END SELECT

         vct ( 1 : nvclev )  = hyai(:)
         vct ( nvclev+1 : )  = hybi(:)
         !
         hyam(:) = ( vct(1:nvclev-1         ) + vct(2:nvclev)          ) / 2.
         hybm(:) = ( vct(nvclev+1:2*nvclev-1) + vct(nvclev+2:2*nvclev) ) / 2.
         !
      ELSE
         hyai => ph
      END IF
    END SUBROUTINE ini_hybrid
    !
  END SUBROUTINE main_grid_def_initialize
  !----------------------------------------------------------------------------

  SUBROUTINE main_grid_def_init_memory

    !MESSy/BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_mpi_bi,           ONLY: reorder
    ! MESSy/SMCL
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
    USE messy_main_channel_repr,  ONLY: get_representation_id
    USE messy_main_blather,       ONLY: start_message, end_message
    USE mo_gaussgrid,             ONLY : gridarea
    
    IMPLICIT NONE
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_init_memory'
    INTEGER :: status
    INTEGER :: GP_2D_HORIZONTAL_ID, GP_3D_MID_ID
    
    ! create new channel
    CALL start_message(submodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)
    
    CALL get_representation_id(status, 'GP_2D_HORIZONTAL', reprid=GP_2D_HORIZONTAL_ID)
    CALL channel_halt(substr, status)
    
    CALL get_representation_id(status, 'GP_3D_MID', reprid=GP_3D_MID_ID)
    CALL channel_halt(substr, status)
    
    CALL new_channel(status, submodstr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'channel_info', c = 'standard basemodel channel' )
    CALL channel_halt(substr, status)
    !
    ! 2D fields
    !
    CALL new_channel_object(status, submodstr, 'gboxarea_2d', p2=gboxarea_2d, reprid=GP_2D_HORIZONTAL_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'gboxarea_2d', 'long_name', c = 'gridbox area')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'gboxarea_2d', 'units', c = 'm^2')
    CALL channel_halt(substr, status)
    !
    ! 3D fields
    !
    CALL new_channel_object(status, submodstr,  'grmass', p3=grmass, reprid=GP_3D_MID_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grmass', 'long_name', c='grid mass')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grmass', 'units', c='kg')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, submodstr,  'grvol', p3=grvol, reprid=GP_3D_MID_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grvol', 'long_name', c='grid volume')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grvol', 'units', c='m3')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, submodstr,  'deltaz', p3=deltaz, reprid=GP_3D_MID_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'deltaz', 'long_name', c='layerthickness')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'deltaz', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL reorder (gboxarea_2d ,SPREAD(gridarea(:),1,nlon))
    
    CALL end_message(submodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)
    
  END SUBROUTINE main_grid_def_init_memory

  !----------------------------------------------------------------------------

  SUBROUTINE main_grid_def_read_restart
    IMPLICIT NONE
  END SUBROUTINE main_grid_def_read_restart

  !----------------------------------------------------------------------------

  SUBROUTINE main_grid_def_global_start
    IMPLICIT NONE
  END SUBROUTINE main_grid_def_global_start

  !----------------------------------------------------------------------------

  SUBROUTINE main_grid_def_local_start
    IMPLICIT NONE
  END SUBROUTINE main_grid_def_local_start
  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_free_memory
  
    USE mo_gaussgrid, ONLY : cleanup_gaussgrid 
    
    IMPLICIT NONE
    
    DEALLOCATE(glat, glon)
    DEALLOCATE(hyam, hybm)
    DEALLOCATE(longitude, latitude)
    DEALLOCATE(cetah, ceta, zf, zh, pf, ph)
    DEALLOCATE(coslon, sinlon, philat, philon, gl_twomu)
    DEALLOCATE(vct, coslon_2d, sinlon_2d, coslat_2d, sinlat_2d, philat_2d)
    
    CALL cleanup_gaussgrid

  END SUBROUTINE main_grid_def_free_memory
  !----------------------------------------------------------------------------
