MODULE messy_ec2cosmo_e5

  ! MESSY
  USE messy_main_constants_mem, ONLY: dp
  USE messy_ec2cosmo

  IMPLICIT NONE
  PRIVATE

  ! soil surface temp.
  REAL(DP), DIMENSION(:,:), POINTER   :: T_S      => NULL() 
  ! soil temperature
  REAL(DP), DIMENSION(:,:,:), POINTER :: T_SO     => NULL() 
  ! relative soil moisture
  REAL(DP), DIMENSION(:,:,:), POINTER :: W_SO_REL => NULL() 
  ! horizontal wind velocities
  REAL(DP), DIMENSION(:,:,:), POINTER :: Uwind => NULL() 
  REAL(DP), DIMENSION(:,:,:), POINTER :: Vwind => NULL() 

  ! SUBROUTINES
  PUBLIC :: ec2cosmo_init_memory
  PUBLIC :: ec2cosmo_init_coupling
  PUBLIC :: ec2cosmo_global_start

CONTAINS

  SUBROUTINE ec2cosmo_init_memory

    ! MESSY INTERFACE
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL          &
                                         , GP_3D_BELOWSF, GP_3D_MID  &
                                         , DIMID_ILEV, DC_BC       
    USE messy_main_channel_repr,     ONLY: new_representation, AUTO  &
                                         , set_representation_decomp & 
                                         , IRANK, PIOTYPE_COL 
    ! MESSY CORE
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr ='ec2cosmo_init_memory'
    INTEGER                     :: status
    INTEGER                     :: REPR_ICOLUMN
    REAL(dp), DIMENSION(:), POINTER :: icolumn
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg  = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL new_representation(status, REPR_ICOLUMN,   &
         'REPR_ICOLUMN'                             &
         , rank = 1, link = 'x---', dctype = DC_BC  &
         , dimension_ids = (/ DIMID_ILEV /)         &
         , ldimlen       = (/ AUTO /)               &
         , axis = 'Z---'                            &
         )
    CALL channel_halt(substr, status)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))
    
    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1
    
    start(:,1) = 1
    cnt(:,1)   = nlev+1
    meml(:,1)  = 1
    memu(:,1)  = nlev+1
    
    CALL set_representation_decomp(status, REPR_ICOLUMN &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! MAKE OUTPUT CHANNEL FOR COSMO BOUNDARY DATA
    CALL new_channel(status, modstr, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)

    ! WRITE GRID ROTATION INFORMATION INTO ATTRIBUTES (FOR INT2COSMO)
    CALL new_attribute(status, 'grid_north_pole_latitude' &
         , r=90.0_dp)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'grid_north_pole_longitud' &
         , r=180.0_dp)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, 'northpole_grid_longitude' &
         , r=0.0_dp)
    CALL channel_halt(substr, status)

    ! DEFINE  horiz. wind velocities
    CALL new_channel_object(status, modstr  &
         , 'Um1', p3=uwind, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'Um1', 'units', c='m/s')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'Um1', 'long_name', c='horz. wind velocity')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr  &
         , 'Vm1', p3=vwind, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'Vm1', 'units', c='m/s')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'Vm1', 'long_name', c='horz. wind velocity')
    CALL channel_halt(substr, status)
    
   ! DEFINE SOIL SURFACE TEMPERATURE T_S FOR COSMO
    CALL new_channel_object(status, modstr  &
         , 'T_S', p2=T_S, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'T_S', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'T_S', 'long_name', c='soil surface temperature')
    CALL channel_halt(substr, status)

    ! DEFINE SOIL TEMPERATURE T_SO FOR COSMO
    CALL new_channel_object(status, modstr  &
         , 'T_SO', p3=T_SO, reprid=GP_3D_BELOWSF)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'T_SO', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'T_SO', 'long_name', c='soil temperature')
    CALL channel_halt(substr, status)
    
    ! DEFINE RELATIVE (SCALED) SOIL MOISTURE FOR COSMO
    CALL new_channel_object(status, modstr  &
         , 'W_SO_REL', p3=W_SO_REL, reprid=GP_3D_BELOWSF)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'W_SO_REL', 'units', c='1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'W_SO_REL', 'long_name', c='relative (scaled) soil moisture')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'W_SO_REL' &
         , 'comment', c='soil moisture scaled on maximum field capacity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'W_SO_REL', 'comment2'      &
         , c='scaled soil moisture =  MIN(1., w_so(i,j) / scale_factor(i,j))')
    CALL channel_halt(substr, status)

    ! dummy channel object to force interface hybrid coefficients into 
    ! channel output
    CALL new_channel_object(status, modstr  &
         , 'ICOLUMN', p1=icolumn, reprid=REPR_ICOLUMN)
    CALL channel_halt(substr, status)

  END SUBROUTINE ec2cosmo_init_memory
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  SUBROUTINE ec2cosmo_init_coupling
    
    ! MESSY INTERFACE
    USE messy_main_grid_def_mem_bi,  ONLY: nn
    USE messy_main_channel_error_bi, ONLY: channel_halt
    
    ! MESSY CORE
    USE messy_main_channel,          ONLY: new_channel_object_reference &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr ='ec2cosmo_init_coupling'
    INTEGER                     :: status   

    ! SET POINTER FOR VARIABLES NEEDED IN COSMO
    ! AND RENAME VARIABLES TO COSMO NAMES

    ! NOTE: AS WE ARE IN GLOBAL START m1 VALUES ARE IDENTICAL TO
    !       UP-TO-DATE VALUES

    ! TEMPERATURE
    CALL new_channel_object_reference(status, &
         'g1a', 'tm1', modstr, 'T', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'T', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'T', 'long_name', c='air temperature')
    CALL channel_halt(substr, status)

    ! surface air pressure PS
    CALL new_channel_object_reference(status, &
         'g3b', 'aps', modstr, 'PS', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'PS', 'units', c='Pa')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'PS' , 'long_name', c='surface pressure')
    CALL channel_halt(substr, status)

    ! specific humidity QV
    CALL new_channel_object_reference(status, &
         'g1a', 'qm1', modstr, 'QV', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'QV', 'units', c='kg kg-1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'QV' , 'long_name', c='specific humidity')
    CALL channel_halt(substr, status)

    ! specific cloud liquid water content QC
    CALL new_channel_object_reference(status, &
         'g1a', 'xlm1', modstr, 'QC', .FALSE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'QC', 'units', c='kg kg-1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'QC' , 'long_name', c='specific cloud liquid water content')
    CALL channel_halt(substr, status)

    ! specific cloud ice content QI
    CALL new_channel_object_reference(status, &
         'g1a', 'xim1', modstr, 'QI', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'QI', 'units', c='kg kg-1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'QI' , 'long_name', c='specific cloud ice content')
    CALL channel_halt(substr, status)

    ! surface geopotential FIS
    CALL new_channel_object_reference(status, &
         'g3b', 'geosp', modstr, 'geosp', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'geosp', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'geosp' , 'long_name', c='surface geopotential')
    CALL channel_halt(substr, status)

    ! canopy water amount W_I
    CALL new_channel_object_reference(status, &
         'g3b', 'wl', modstr, 'wl', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'wl', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'wl' , 'long_name', c='canopy water amount')
    CALL channel_halt(substr, status)

    ! lwe_thickness_of_surface_snow_amount ! W_SNOW
    CALL new_channel_object_reference(status, &
         'g3b', 'sni', modstr, 'sni', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'sni', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'sni' , 'long_name', c='surface snow amount')
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, &
         'g3b', 'sn', modstr, 'sn', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'sn', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'sn' , 'long_name', c='snow depth')
    CALL channel_halt(substr, status)

    ! land fraction 
    CALL new_channel_object_reference(status, &
         'g3b', 'slf', modstr, 'FR_LAND', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'FR_LAND', 'units', c='1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'FR_LAND' , 'long_name', c='land-sea fraction')
    CALL channel_halt(substr, status)

    ! land sea mask 
    CALL new_channel_object_reference(status, &
         'g3b', 'slm', modstr, 'FR_LAND_MASK', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'FR_LAND_MASK', 'units', c='1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'FR_LAND_MASK' , 'long_name', c='land-sea mask')
    CALL channel_halt(substr, status)

    ! surface temperature 
    CALL new_channel_object_reference(status, &
         'g3b', 'tslm1', modstr, 'tslm1', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'tslm1', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'tslm1' , 'long_name', c='surface temperature of land')
    CALL channel_halt(substr, status)

    ! snow temperature T_SNOW
    CALL new_channel_object_reference(status, &
         'g3b', 'tsi', modstr, 'tsi', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'tsi', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'tsi' , 'long_name', c='surface temperature of ice')
    CALL channel_halt(substr, status)

    ! water temperature
    CALL new_channel_object_reference(status, &
         'g3b', 'tsw', modstr, 'tsw', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'tsw', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'tsw' , 'long_name', c='surface temperature of water')
    CALL channel_halt(substr, status)

    ! soil temperature
    CALL new_channel_object_reference(status, &
         'g3b', 'tsoil', modstr, 'tsoil', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'tsoil', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'tsoil' , 'long_name', c='deep soil temperatures')
    CALL channel_halt(substr, status)

    ! soil wetness
    CALL new_channel_object_reference(status, &
         'g3b', 'ws', modstr, 'ws', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'ws', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'ws' , 'long_name', c='soil wetness')
    CALL channel_halt(substr, status)

    ! field capacity of soil
    CALL new_channel_object_reference(status, &
         'g3b', 'wsmx', modstr, 'wsmx', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'wsmx', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'wsmx' , 'long_name', c='field capacity of soil')
    CALL channel_halt(substr, status)

    ! water temperature
    CALL new_channel_object_reference(status, &
         'g3b', 'seaice', modstr, 'seaice', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'seaice', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'seaice' , 'long_name', c='ice cover (fraction of 1-SLM)')
    CALL channel_halt(substr, status)

    ! specific humidity
    CALL new_channel_object_reference(status, &
         'gl', 'q', modstr, 'q', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'q', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'q' , 'long_name', c='specific humidity')
    CALL channel_halt(substr, status)

    ! cloud water
    CALL new_channel_object_reference(status, &
         'gl', 'xl', modstr, 'xl', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'xl', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'xl' , 'long_name', c='cloud water')
    CALL channel_halt(substr, status)

    ! cloud ice
    CALL new_channel_object_reference(status, &
         'gl', 'xi', modstr, 'xi', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'xi', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'xi' , 'long_name', c='cloud ice')
    CALL channel_halt(substr, status)

    ! mz_pj_20081110+
    ! add dummy field on interface levels
    CALL new_channel_object_reference(status, &
         'ECHAM5', 'pressi', modstr, 'PRESSI', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'PRESSI', 'units', c='Pa')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'PRESSI' , 'long_name', c='pressure at interface levels')
    CALL channel_halt(substr, status)
    ! mz_pj_20081110-

    CALL new_channel_object_reference(status, &
         'sp', 'sd', modstr, 'sd', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'sd', 'units', c='1/s')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'sd' , 'long_name', c='divergence')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'sd' , 'grid_type', c='spectral, triangular truncation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'sd' , 'code', i=155)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'sd' , 'table', i=128)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'sd' , 'axis', c='tz--')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'sd' , 'truncation', i=nn)
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, &
         'sp', 'svo', modstr, 'svo', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'svo', 'units', c='1/s')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'svo' , 'long_name', c='vorticity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'svo' , 'grid_type', c='spectral, triangular truncation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'svo' , 'code', i=138)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'svo' , 'table', i=128)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'svo' , 'axis', c='tz--')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'svo' , 'truncation', i=nn)
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, &
         'sp', 'st', modstr, 'st', .FALSE.) 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'st', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'st' , 'long_name', c='temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'st' , 'grid_type', c='spectral, triangular truncation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'st' , 'code', i=130)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'st' , 'table', i=128)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'st' , 'axis', c='tz--')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'st' , 'truncation', i=nn)
    CALL channel_halt(substr, status)

  END SUBROUTINE ec2cosmo_init_coupling
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  SUBROUTINE ec2cosmo_global_start

    USE messy_main_data_bi,  ONLY: tslm1, tsi, tsw, slf, seaice, um1, vm1   &
                                 , wsmx, ws, tsoil
    USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma, npromz, nlev
    USE messy_main_grid_def_bi,     ONLY: sqcst_2d

    IMPLICIT NONE

    REAL(DP), PARAMETER :: undef = -1.E20_dp ! value for undefined value
    INTEGER             :: jk, jp, jrow, kproma

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr ='ec2cosmo_global_start'

    ! CALCULATE SOIL TEMPERATURE T_S and  RELATIVE SOIL MOISTURE FOR COSMO
    DO jrow=1,ngpblks
       IF ( jrow == ngpblks ) THEN
          kproma = npromz
       ELSE
          kproma = nproma
       END IF
       DO jp = 1, kproma
          IF (slf(jp,jrow) < 0.5_dp ) THEN
             t_s(jp,jrow) = seaice(jp,jrow)*tsi(jp,jrow) &
                  + (1.-seaice(jp,jrow))*tsw(jp,jrow)
            DO jk = 1,5   ! number of soil layers is hardcoded within ECHAM5
                t_so(jp,jk,jrow)      = -1.E20 ! value for undefined value
                w_so_rel(jp,jk,jrow)  = -1.E20 ! value for undefined value
             END DO
          ELSE
             t_s(jp,jrow) = tslm1(jp,jrow)
             DO jk = 1,5   ! number of soil layers is hardcoded within ECHAM5
                t_so(jp,jk,jrow)      = tsoil(jp,jk,jrow) 
                w_so_rel(jp,jk,jrow)  = MIN(1., ws(jp,jrow) / wsmx(jp,jrow))
             END DO
          END IF
          DO jk = 1,nlev
             uwind(jp,jk,jrow) = um1(jp,jk,jrow)  / sqcst_2d(jp,jrow)
             vwind(jp,jk,jrow) = vm1(jp,jk,jrow)  / sqcst_2d(jp,jrow)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE ec2cosmo_global_start
  ! -----------------------------------------------------------------------

END MODULE messy_ec2cosmo_e5
