#include "messy_main_ppd_bi.inc"
#define FAST_BUT_UGLY
! **************************************************************************
MODULE messy_tropop_si
! **************************************************************************

  ! MODULE FOR VARIOUS TROPOPAUSE DIAGNOSTICS
  !
  ! Authors: Patrick Joeckel, MPICH, Sep 2003
  !          - original code
  !          Michael Traub,   MPICH, Nov 2003 
  !          - planetary boundary layer height calculations
  !          Astrid Kerkweg, UniMz, June 2011
  !          - planetary boundary layer height calculations
  !          - ECMWF sea level pressure calculation

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_tools,         ONLY: iso2ind, ind2val
  USE messy_tropop

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  ! NEW CHANNEL OBJECTS: (DIAGNOSED FIELDS)
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_clim   ! climatol. tp pressure
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_clim_i ! level index of climatol. tp
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_clim_f ! fract. of box in troposph.
  REAL(DP), DIMENSION(:,:,:), POINTER :: PV        ! potential vorticity
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_PV     ! PV tropopause (pressure)
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_PV_i   ! PV tropopause (index)
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_PV_f   ! PV fract. of box in tropos.
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_WMO    ! WMO tropopause (pressure)
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_WMO_i  ! WMO tropopause (index)
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_WMO_f  ! WMO fract. of box in trsph.
  REAL(DP), DIMENSION(:,:),   POINTER :: tp       ! PV(< r_lat ) + WMO(> r_lat)
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_i     ! ...
  REAL(DP), DIMENSION(:,:),   POINTER :: tp_f     ! ...
  REAL(DP), DIMENSION(:,:,:), POINTER :: O3_PV    ! parameterized strat. OZONE
  REAL(DP), DIMENSION(:,:,:), POINTER :: N2O_P    ! parameterized strat. N2O
  REAL(DP), DIMENSION(:,:,:), POINTER :: NOy_P    ! parameterized strat. NOy

  ! DIAGNOSTIC OUTPUT FOR BOUNDARY LAYER HEIGHT CALCULATIONS
  REAL(DP), DIMENSION(:,:),   POINTER :: pblh    ! planetary bnd. layer height
  REAL(DP), DIMENSION(:,:),   POINTER :: pblh_i  ! index of pblh
  ! planetary bnd. layer height with bulk richardson number 
  REAL(DP), DIMENSION(:,:),   POINTER :: pblhRi  => NULL()  
  REAL(DP), DIMENSION(:,:),   POINTER :: pblhRi_i=> NULL() ! index of pblh
#ifdef COSMO
  ! bulk richardson number
  REAL(DP), DIMENSION(:,:,:), POINTER :: brn     => NULL() 

  REAL(DP), DIMENSION(:,:,:), POINTER :: u_nrot  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: v_nrot  => NULL()

  REAL(DP), DIMENSION(:,:,:), POINTER :: TI      => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: N2      => NULL()
#endif
  ! sea level pressure (hPa)
  REAL(DP), DIMENSION(:,:),   POINTER :: slp => NULL()  

  ! cold point diagnostics
  REAL(DP), DIMENSION(:,:), POINTER :: cpt    => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: cpt_i  => NULL()

  REAL(dp), DIMENSION(:,:,:), POINTER :: windspeed => NULL() 
  
  PUBLIC :: tropop_initialize
  PUBLIC :: tropop_init_memory
  PUBLIC :: tropop_vdiff

CONTAINS

! -------------------------------------------------------------------------
SUBROUTINE tropop_initialize

  ! ECHAM5/MESSy
  USE messy_main_blather_bi, ONLY: error_bi !um_ak_20110624
  USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
  USE messy_main_tools,      ONLY: find_next_free_unit

  IMPLICIT NONE

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'tropop_initialize'
  INTEGER                     :: iou    ! I/O unit
  INTEGER                     :: status ! error status

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL tropop_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('error in namelist CTRL ',substr)
    END IF

    CALL p_bcast(r_climtp(:), p_io)
    CALL p_bcast(l_wmo_clim_corr, p_io)
    CALL p_bcast(r_lat, p_io)
    CALL p_bcast(r_dyntp_PV, p_io)
    CALL p_bcast(r_press_range_PV(:), p_io)
    CALL p_bcast(l_o3_PV, p_io)
    CALL p_bcast(l_n2o, p_io)
    CALL p_bcast(l_noy, p_io)
    CALL p_bcast(l_pblh, p_io)
    CALL p_bcast(l_slp, p_io)
    CALL p_bcast(l_TI, p_io) !um_ch_20150127
    CALL p_bcast(l_N2, p_io) !um_ch_20150128
    CALL p_bcast(l_cpt, p_io) ! op_sd_201211214
    CALL p_bcast(l_windspeed, p_io)
    
  END SUBROUTINE tropop_initialize
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE tropop_init_memory

  ! TROPOP MODULE ROUTINE (ECHAM-5 INTERFACE)
  !
  ! define tropop specific channel(s) and allocate memory for
  ! global fields
  !
  ! Author: Patrick Joeckel, MPICH, Sep 2003

  ! ECHAM5/MESSy
  USE messy_main_channel_error_bi, ONLY: channel_halt
  USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_2D_HORIZONTAL
  ! MESSy
  USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                       , new_attribute

  IMPLICIT NONE

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'tropop_init_memory'
  INTEGER                     :: status
  
  CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)
  
  ! define new channel
  CALL new_channel(status, modstr, lrestreq=.TRUE.)
  CALL channel_halt(substr, status)
  
  ! TROPOP DIAGNOSTIC OUTPUT
  CALL new_channel_object(status, modstr, 'tp_clim', p2=tp_clim &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_clim', 'long_name', c='clim. tropopause pressure' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_clim', 'units', c='Pa' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_clim_i', p2=tp_clim_i &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_clim_i', 'long_name', c='clim. tropopause level index' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_clim_i', 'units', c='' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_clim_f', p2=tp_clim_f &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_clim_f', 'long_name', c='fract. of clim.-TP box in trop.' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_clim_f', 'units', c='' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'PV', p3=PV &
       , reprid=GP_3D_MID )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'PV', 'long_name', c='potential vorticity' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'PV', 'units', c='PVU' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_PV', p2=tp_PV &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_PV', 'long_name', c='PV tropopause pressure' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_PV', 'units', c='Pa' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_PV_i', p2=tp_PV_i &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_PV_i', 'long_name', c='PV tropopause level index' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_PV_i', 'units', c='' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_PV_f', p2=tp_PV_f &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_PV_f', 'long_name', c='fract. of PV-TP box in trop.' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_PV_f', 'units', c='' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_WMO', p2=tp_WMO &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_WMO', 'long_name', c='WMO tropopause pressure' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_WMO', 'units', c='Pa' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_WMO_i', p2=tp_WMO_i &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_WMO_i', 'long_name', c='WMO tropopause level index' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_WMO_i', 'units', c='' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_WMO_f', p2=tp_WMO_f &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_WMO_f', 'long_name', c='fract. of WMO-TP box in trop.' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_WMO_f', 'units', c='' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp', p2=tp &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp', 'long_name', c='tropopause pressure' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp', 'units', c='Pa' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_i', p2=tp_i &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_i', 'long_name', c='tropopause level index' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_i', 'units', c='' )
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tp_f', p2=tp_f &
       , reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_f', 'long_name', c='fract. of TP box in trop.' )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr &
       , 'tp_f', 'units', c='' )
  CALL channel_halt(substr, status)

  IF (l_O3_PV) THEN
     CALL new_channel_object(status, modstr, 'O3_PV', p3=o3_pv &
          , reprid=GP_3D_MID )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'O3_PV', 'long_name', c='parameterized stratosph. O3(PV)' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'O3_PV', 'units', c='mol/mol' )
     CALL channel_halt(substr, status)
  END IF

  IF (l_N2O) THEN
     CALL new_channel_object(status, modstr, 'N2O_P', p3=n2o_p &
          , reprid=GP_3D_MID )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'N2O_P', 'long_name', c='parameterized stratosph. N2O(O3)' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'N2O_P', 'units', c='mol/mol' )
     CALL channel_halt(substr, status)
  END IF

  IF (l_NOy) THEN
     CALL new_channel_object(status, modstr, 'NOy_P', p3=noy_p &
          , reprid=GP_3D_MID )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'NOy_P', 'long_name', c='parameterized stratosph. NOy(N2O)' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'NOy_P', 'units', c='mol/mol' )
     CALL channel_halt(substr, status)
  END IF

  ! mz_mt_20031119+
  ! BOUNDARY LAYER
  IF (l_PBLH) THEN
#if defined (ECHAM5) || defined(COSMO)
     CALL new_channel_object(status, modstr, 'pblh', p2=pblh &
          , reprid=GP_2D_HORIZONTAL )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'pblh', 'long_name', c='planetary boundary layer' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'pblh', 'units', c='m' )
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'pblh_i', p2=pblh_i &
          , reprid=GP_2D_HORIZONTAL )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'pblh_i', 'long_name', c='planetary boundary layer index' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'pblh_i', 'units', c='' )
     CALL channel_halt(substr, status)

     ! um_ak_20110622+ 
     CALL new_channel_object(status, modstr, 'pblhRi', p2=pblhRi &
          , reprid=GP_2D_HORIZONTAL )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'pblhRi', 'long_name', c='planetary boundary layer' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'pblhRi', 'units', c='m' )
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'pblhRi_i', p2=pblhRi_i &
          , reprid=GP_2D_HORIZONTAL )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'pblhRi_i', 'long_name', c='planetary boundary layer index' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'pblhRi_i', 'units', c='' )
     CALL channel_halt(substr, status)
#endif

#ifdef COSMO
     CALL new_channel_object(status, modstr, 'bulkRi', p3=brn &
          , reprid=GP_3D_MID )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'bulkRi', 'long_name', c='bulk Richardson number' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'bulkRi', 'units', c='-' )
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'u_nrot', p3=u_nrot &
          , reprid=GP_3D_MID )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'u_nrot', 'long_name', c='zonal wind on geogr. grid' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'u_nrot', 'units', c='m/s' )
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'v_nrot', p3=v_nrot &
          , reprid=GP_3D_MID )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'v_nrot', 'long_name', c='meridional wind on geogr. grid' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'v_nrot', 'units', c='m/s' )
     CALL channel_halt(substr, status)

     IF (l_TI) THEN
        CALL new_channel_object(status, modstr, 'TI', p3=TI &
             , reprid=GP_3D_MID )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr &
             , 'TI', 'long_name', c='turbulence index' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr &
             , 'TI', 'units', c=' ' )
        CALL channel_halt(substr, status)
     END IF
     IF (l_N2) THEN
        CALL new_channel_object(status, modstr, 'N2', p3=N2 &
             , reprid=GP_3D_MID )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr &
             , 'N2', 'long_name', c='Brunt-Vaisala frequency' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr &
             , 'N2', 'units', c=' ' )
        CALL channel_halt(substr, status)
     END IF
#endif
  END IF

  IF (l_slp) THEN
     CALL new_channel_object(status, modstr, 'slp', p2=slp &
          , reprid=GP_2D_HORIZONTAL )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'slp', 'long_name', c='sea-level pressure' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'slp', 'units', c='Pa' )
     CALL channel_halt(substr, status)
 END IF

  IF (l_cpt) THEN
     CALL new_channel_object(status, modstr, 'cpt', p2=cpt &
          , reprid=GP_2D_HORIZONTAL )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'cpt', 'long_name', c='pressure at cold point' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'cpt', 'units', c='Pa' )
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'cpt_i', p2=cpt_i &
          , reprid=GP_2D_HORIZONTAL )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr &
          , 'cpt_i', 'long_name', c='cold point level index' )
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'cpt_i', 'units', c='' )
     CALL channel_halt(substr, status)     
   END IF

   IF (l_windspeed) THEN
      CALL new_channel_object(status, modstr, 'windspeed', p3=windspeed &
           , reprid=GP_3D_MID )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr &
           , 'windspeed', 'long_name', c='wind speed at grid midpoints' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'windspeed', 'units', c='m s-1' )
      CALL channel_halt(substr, status)
   END IF
   
 CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

END SUBROUTINE tropop_init_memory
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE tropop_vdiff

  ! BMIL
  USE messy_main_data_bi,     ONLY: press_3d, geopot_3d, um1, vm1 &
                                  , qm1_3d, tm1_3d                &
                                  , tm1, vom1                     &
                                  , aps                           &
                                  , pressi_3d                     &
                                  , tpot_3d, coriol_2d
  USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
  USE messy_main_grid_def_bi,     ONLY: hyam, hybm, grvol, grmass, philat_2d &
                                      , altitudei_msl, altitudei_gnd         &
                                      , altitude_gnd
#if defined(ECHAM5)
  ! variables only required for ECHAM-only PBL calculation
  USE messy_main_data_bi,     ONLY: tpot_3d             &
                                  , zust_2d             &
                                  , zlatkf_2d           &
                                  , zsenkf_2d           &
                                  , rinum_3d
#endif

#if defined(ECHAM5) || defined(CESM1)
   USE messy_main_data_bi,         ONLY: geosp
#endif

#ifdef COSMO
   USE messy_main_data_bi,         ONLY: qv_2m, t_2m                  &
                                       , zsenkf_2d, zlatkf_2d, zust_2d 
   USE messy_main_grid_def_mem_bi, ONLY: pollon, pollat     &
                                       , dlon, dlat, ngpblks
   USE messy_main_grid_def_bi,     ONLY: crlat, philon_2d, hhl, hsurf
   USE messy_main_grid_trafo,      ONLY: uvrot2uv_vec_jloop
#endif

  USE messy_main_timer,         ONLY: MONTH
  USE messy_main_constants_mem, ONLY: pi, g, cp_air, rd

  IMPLICIT NONE
  INTRINSIC :: ABS, COS, MAX, MERGE, MIN, REAL

  ! LOCAL
#ifndef FAST_BUT_UGLY
  INTEGER :: jp
#else
  INTEGER :: jp
#ifdef ECHAM5
  INTEGER :: jk
#endif
#endif
  INTEGER,     DIMENSION(:),   ALLOCATABLE :: tp_k_pv, tp_k_wmo
  REAL(DP),    DIMENSION(:),   ALLOCATABLE :: zcoriol, zphi
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE :: philat_z
#ifndef FAST_BUT_UGLY
  REAL(DP) :: zpclim                   ! pressure of clim. tropopause [Pa]
  INTEGER  :: iclim                    ! level index of clim. tropopause
  INTEGER  :: kmin, kmax               ! level index range for PV tropopause
  INTEGER  :: k
#else
  REAL(DP), DIMENSION(kproma)   :: zpclim  ! pressure of clim. tropopause [Pa]
  INTEGER, DIMENSION(kproma)    :: iclim   ! level index of clim. tropopause
  INTEGER, DIMENSION(kproma)    :: k, kmin, kmax
#endif
  REAL(dp), DIMENSION(nlev)     :: press0
  REAL(dp), DIMENSION(kproma)   :: rho0_ll   ! density at the full model levels(kg/m3)
  REAL(dp), DIMENSION(kproma)   :: dp0_ll    ! pressure thickness of model layers (Pa)
  !
  INTEGER :: i_cpt_i
  
#ifdef ECHAM5 
  REAL(dp)                      :: hhl(kproma,nlev+1)
#endif 

  ! ALLOCATE SPACE FOR TEMPORARY LOCAL FIELDS
  ALLOCATE(zcoriol(kproma))
  ALLOCATE(zphi(kproma))
  ALLOCATE(tp_k_pv(kproma))
  ALLOCATE(tp_k_wmo(kproma))
  ALLOCATE(philat_z(kproma, nlev))

  ! GET VALUES FROM GLOBAL FIELDS
  DO jp=1, kproma
     zcoriol(jp) = coriol_2d(jp,jrow)
     zphi(jp)    = (philat_2d(jp,jrow)/180.)*pi
     philat_z(jp,1:nlev) = philat_2d(jp,jrow)
  END DO

  ! 0) climatological tropopause
#ifndef FAST_BUT_UGLY
  DO jp=1, kproma
     zpclim = r_climtp(1)*100._dp &
          - r_climtp(2)*100._dp*cos(zphi(jp))*cos(zphi(jp))
     tp_clim(jp, jrow) = zpclim
     CALL iso2ind(press_3d(_RI_XYZ__(jp,jrow,:)), zpclim,   &
          iclim, f=tp_clim_f(jp,jrow), lrev=.true.)
     tp_clim_i(jp, jrow) = REAL(iclim,dp)
  END DO
#else
  zpclim = r_climtp(1)*100._dp &
         - r_climtp(2)*100._dp*cos(zphi(1:kproma))*cos(zphi(1:kproma))
  tp_clim(1:kproma, jrow) = zpclim(:)
  CALL iso2ind(kproma, press_3d(_RI_XYZ__(1:kproma,jrow,:)), zpclim,   &
       iclim, f=tp_clim_f(1:kproma,jrow), lrev=.true.)
  tp_clim_i(1:kproma, jrow) = REAL(iclim,dp) 
#endif

  ! 1) PV-tropopause
  ! 1a) PV
#ifndef FAST_BUT_UGLY
  DO jp=1, kproma
     CALL calc_PV(PV(_RI_XYZ__(jp,jrow,:))          &
          ,press_3d(_RI_XYZ__(jp,jrow,:))           &
          ,tm1(_RI_XYZ__(jp,jrow,:))                &
          ,vom1(_RI_XYZ__(jp,jrow,:)), zcoriol(jp)  &
          ,rd, cp_air, g)
  END DO
#else
  CALL calc_PV(kproma, PV(_RI_XYZ__(1:kproma,jrow,:))         &
       ,press_3d(_RI_XYZ__(1:kproma,jrow,:))                  &
       ,tm1(_RI_XYZ__(1:kproma,jrow,:))                       &
       ,vom1(_RI_XYZ__(1:kproma,jrow,:)), zcoriol(1:kproma)   &
       ,rd, cp_air, g)
#endif

  ! USE ABSOLUTE VALUE OF PV 
  PV(_RI_XYZ__(:,jrow,:))=ABS(PV(_RI_XYZ__(:,jrow,:)))

#ifndef FAST_BUT_UGLY
  ! 1b) level index and fraction at PV = const. = r_dyntp_PV
  DO jp=1,kproma
     ! pressure-range to index range
     CALL iso2ind(press_3d(_RI_XYZ__(jp,jrow,:)) &
          , MIN(r_press_range_PV(1),r_press_range_PV(2)), kmin, lrev=.true.)
     CALL iso2ind(press_3d(_RI_XYZ__(jp,jrow,:)) &
          , MAX(r_press_range_PV(1),r_press_range_PV(2)), kmax, lrev=.false.)
     kmin = MAX(kmin, 2)      ! PV is zero at k=1
     kmax = MIN(kmax, nlev-2) ! PV is zero at k=NLEV and k=NLEV-1
     !
     ! PV is zero at K=1, K=NLEV, and K=NLEV-1
     CALL iso2ind(PV(_RI_XYZ__(jp,jrow,kmin:kmax)), r_dyntp_PV, &
          tp_k_pv(jp), f=tp_PV_f(jp,jrow), lrev=.false.)
     ! add offset, which has been skipped in iso2ind
     ! (see array boundaries there)
     tp_k_pv(jp) = tp_k_pv(jp) + kmin - 1
  END DO
  tp_PV_i(1:kproma,jrow) = REAL(tp_k_pv(:),dp) 

  ! 1c) pressure at tropopause
  DO jp=1,kproma
     CALL ind2val(tp_PV(jp,jrow), press_3d(_RI_XYZ__(jp,jrow,:)), &
          tp_k_pv(jp), tp_PV_f(jp, jrow))
  END DO
#else
  ! 1b) level index and fraction at PV = const. = r_dyntp_PV
  ! pressure-range to index range
  zpclim = MIN(r_press_range_PV(1),r_press_range_PV(2))
  CALL iso2ind(kproma, press_3d(_RI_XYZ__(1:kproma,jrow,:)), zpclim, kmin, lrev=.true.)
  zpclim = MAX(r_press_range_PV(1),r_press_range_PV(2))
  CALL iso2ind(kproma, press_3d(_RI_XYZ__(1:kproma,jrow,:)), zpclim, kmax, lrev=.false.)
  kmin = MAX(kmin, 2)      ! PV is zero at k=1
  kmax = MIN(kmax, nlev-2) ! PV is zero at k=NLEV and k=NLEV-1
  !
  DO jp=1, kproma
     CALL iso2ind(PV(_RI_XYZ__(jp,jrow,kmin(jp):kmax(jp))), r_dyntp_PV, &
          tp_k_pv(jp), f=tp_PV_f(jp,jrow), lrev=.false.)
  END DO

  tp_k_pv(1:kproma) = tp_k_pv(1:kproma) + kmin - 1
  tp_pv_i(1:kproma,jrow) = REAL(tp_k_pv(:),dp) 

  ! 1c) pressure at tropopause
  CALL ind2val(kproma, tp_pv(1:kproma,jrow), press_3d(_RI_XYZ__(1:kproma,jrow,:)), &
       tp_k_pv(1:kproma), tp_pv_f(1:kproma, jrow))
#endif

  ! 2) WMO-tropopause (with optional climatological filling)
  ! 2a) pressure at tropopause

  !  calculate standard pressure axis in Pascal:
  press0(:) = hyam(:) + hybm(:) * 101325._dp

  IF (l_wmo_clim_corr) THEN
     CALL WMOtropop(tp_WMO(1:kproma,jrow), kproma, nlev  &
          ,tm1(_RI_XYZ__(1:kproma,jrow,:)), press_3d(_RI_XYZ__(1:kproma,jrow,:)) &
          ,press0, g, rd, cp_air &
          ,zphi(:) ) ! latitude [rad] used for climatolog. fill.
  ELSE
     CALL WMOtropop(tp_WMO(1:kproma,jrow), kproma, nlev  &
          ,tm1(_RI_XYZ__(1:kproma,jrow,:)), press_3d(_RI_XYZ__(1:kproma,jrow,:)) &
          ,press0 ,g, rd, cp_air)            ! un-filled
  END IF

  ! 2b) level index and fraction
#ifndef FAST_BUT_UGLY
  DO jp=1, kproma
     CALL iso2ind(press_3d(_RI_XYZ__(jp,jrow,:)), tp_WMO(jp, jrow), &
          tp_k_wmo(jp), f=tp_WMO_f(jp, jrow))
  END DO
#else
  CALL iso2ind(kproma, press_3d(_RI_XYZ__(1:kproma,jrow,:)), tp_wmo(1:kproma, jrow), &
       tp_k_wmo(1:kproma), f=tp_wmo_f(1:kproma, jrow))
#endif
  tp_WMO_i(1:kproma,jrow) = REAL(tp_k_wmo(:),dp)

  ! 3) combine PV-tropopause and WMO tropopasue (intersect at r_lat)
  tp(1:kproma,jrow) = MERGE(tp_PV(1:kproma,jrow) &
       , tp_WMO(1:kproma,jrow), ABS(zphi(:)) >= (r_lat/180._dp)*pi )
  tp_i(1:kproma,jrow) = MERGE(tp_PV_i(1:kproma,jrow) &
       , tp_WMO_i(1:kproma,jrow), ABS(zphi(:)) >= (r_lat/180._dp)*pi )
  tp_f(1:kproma,jrow) = MERGE(tp_PV_f(1:kproma,jrow) &
       , tp_WMO_f(1:kproma,jrow), ABS(zphi(:)) >= (r_lat/180._dp)*pi )

  IF (L_O3_PV) THEN
     ! 4) calulate parameterized stratospheric ozone
     CALL pv_to_o3(o3_pv(_RI_XYZ__(1:kproma,jrow,:)), pv(_RI_XYZ__(1:kproma,jrow,:)), MONTH)
  END IF

  IF (L_N2O) THEN
     CALL o3_to_n2o(n2o_p(_RI_XYZ__(1:kproma,jrow,:)), o3_pv(_RI_XYZ__(1:kproma,jrow,:))  &
          , philat_z(1:kproma,:))
     ! SET TO ZERO WHERE O3 is ZERO
     n2o_p(_RI_XYZ__(:,jrow,:)) = &
          MERGE(n2o_p(_RI_XYZ__(:,jrow,:)), 0.0_dp, o3_pv(_RI_XYZ__(:,jrow,:)) > 0.0_dp)
  END IF

  IF (L_NOy) THEN
     CALL n2o_to_noy(noy_p(_RI_XYZ__(1:kproma,jrow,:)), n2o_p(_RI_XYZ__(1:kproma,jrow,:)) &
          , philat_z(1:kproma,:))
     ! SET TO ZERO WHERE N2O is ZERO
     noy_p(_RI_XYZ__(:,jrow,:)) = &
          MERGE(noy_p(_RI_XYZ__(:,jrow,:)), 0.0_dp, n2o_p(_RI_XYZ__(:,jrow,:)) > 0.0_dp)
  END IF

  ! CALCULATION OF POTENTIAL TEMPERATURE AT SURFACE
  ! NEEDED FOR PBLH 
  IF (l_pblh) THEN

#if defined(ECHAM5) || defined(COSMO)
     ! 1.) first parameterisation
     ! ALLOCATE SPACE FOR TEMPORARY LOCAL FIELDS
     !ALLOCATE(theta_surf(kproma))
     !
     !theta_surf(:)=tsurf_2d(1:kproma,jrow)  & 
     !     *(1.E5_dp/pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)))**(rd/cp_air)
     !
     CALL pblheight(pblh(1:kproma, jrow), kproma, nlev              &
          , tpot_3d(_RI_XYZ__(1:kproma,jrow,:))                               &
          , qm1_3d(_RI_XYZ__(1:kproma,jrow,:)), geopot_3d(_RI_XYZ__(1:kproma,jrow,:))   &
          , um1(_RI_XYZ__(1:kproma,jrow,:)), vm1(_RI_XYZ__(1:kproma,jrow,:))            &
          !, theta_surf(1:kproma)                                    & 
          , zsenkf_2d(1:kproma, jrow), zlatkf_2d(1:kproma, jrow)    &
          , zust_2d(1:kproma, jrow) &!, cdnl(1:kproma, jrow)           &
          !, zcdh_2d(1:kproma, jrow), cfml(1:kproma, jrow)           &
          , g=g )
     !
     ! CALCULATE PBL INDEX IN ECHAM
     ! THIS INFORMATION IS NEEDED IN ATTILA
     !
#ifndef FAST_BUT_UGLY
     DO jp=1, kproma
        CALL pblindex(k, pblh(jp,jrow), altitudei_gnd(_RI_XYZ__(jp,jrow,:)))
        pblh_i(jp, jrow) = REAL(k,dp)
     END DO
#else
     CALL pblindex(kproma, k, pblh(1:kproma,jrow), &
          altitudei_gnd(_RI_XYZ__(1:kproma,jrow,:)))
     pblh_i(1:kproma, jrow) = REAL(k,dp)
#endif
     !
     ! CLEAN MEMORY
     !DEALLOCATE(theta_surf)
#endif

#ifdef COSMO
     ! 2a). CALCULATE richardson number
     CALL calc_bulk_richardson(kproma, ngpblks, nlev, jrow  &
       , brn(_RI_XYZ__(1:kproma,jrow,:))               &
       , tm1_3d(_RI_XYZ__(1:kproma,jrow,:))            &
       , qm1_3d(_RI_XYZ__(1:kproma,jrow,:))            &
       , um1(1:kproma,:,:), vm1(1:kproma,:,:)          &
       , press_3d(_RI_XYZ__(1:kproma,jrow,:))          &
       , hsurf(1:kproma,jrow), aps(1:kproma,jrow)      &
#ifndef COSMOv509
       , t_2m(1:kproma,jrow), qv_2m(1:kproma,jrow)     &
#else
       , t_2m(1:kproma,jrow,0), qv_2m(1:kproma,jrow,0) &
#endif
       , hhl(1:kproma,jrow,:)           )
#endif

#ifdef ECHAM5  
      ! 2b). CALCULATE richardson number                                                                                                  
     ! calculate level height 
     hhl(1:kproma,1) = 0.0_dp
     do jk=2,nlev+1
        hhl(1:kproma,jk)= altitudei_msl(_RI_XYZ__(1:kproma,jrow,jk))
     end do

     CALL  calc_pbl_brn(kproma, nlev, tm1_3d(_RI_XYZ__(1:kproma,jrow,:))  &
          , qm1_3d(_RI_XYZ__(1:kproma,jrow,:)), press_3d(_RI_XYZ__(1:kproma,jrow,:))   &
          , hhl(1:kproma,:)     &
          , geosp(1:kproma,jrow)/g             &
          , rinum_3d(_RI_XYZ__(1:kproma,jrow,:)) &
          , pblhRi(1:kproma,jrow), pblhRi_i(1:kproma,jrow) )
#endif

#ifdef COSMO
      ! 2b). CALCULATE boundary layer height
   CALL  calc_pbl_brn(kproma, nlev, tm1_3d(_RI_XYZ__(1:kproma,jrow,:))  &
       , qm1_3d(_RI_XYZ__(1:kproma,jrow,:)), press_3d(_RI_XYZ__(1:kproma,jrow,:))   &
       , hhl(1:kproma,jrow,:), hsurf(1:kproma,jrow)             &
       , brn(_RI_XYZ__(1:kproma,jrow,:)) &
       , pblhRi(1:kproma,jrow), pblhRi_i(1:kproma,jrow) )

   u_nrot(_RI_XYZ__(1:kproma,jrow,:))=um1(_RI_XYZ__(1:kproma,jrow,:))
   v_nrot(_RI_XYZ__(1:kproma,jrow,:))=vm1(_RI_XYZ__(1:kproma,jrow,:))

   CALL uvrot2uv_vec_jloop(u_nrot(_RI_XYZ__(1:kproma,jrow,:))          &
        , v_nrot(_RI_XYZ__(1:kproma,jrow,:))                           &
        , philat_2d(1:kproma,jrow), philon_2d(1:kproma,jrow) &
        , pollat, pollon, kproma, nlev)

   IF (l_TI) THEN
      CALL calc_TI(kproma, ngpblks, nlev, jrow, TI(_RI_XYZ__(1:kproma,jrow,:)) &
           , um1(1:kproma,:,:), vm1(1:kproma,:,:), altitude_gnd(1:kproma,:,:) &
           , dlon, dlat, &
           crlat(jrow,1))
   END IF
   IF (l_N2) THEN
      CALL calc_N2(kproma, nlev, N2(_RI_XYZ__(1:kproma,jrow,:)) &
           , tpot_3d(_RI_XYZ__(1:kproma,jrow,:)), altitude_gnd(_RI_XYZ__(1:kproma,jrow,:)))
   END IF
#endif

  END IF

  ! CLEAN MEMORY
  DEALLOCATE(zcoriol)
  DEALLOCATE(zphi)
  DEALLOCATE(tp_k_pv)
  DEALLOCATE(tp_k_wmo)
  DEALLOCATE(philat_z)

  IF (l_slp) THEN
     ! Calculate pressure thickness of model layers(Pa)
     dp0_ll(1:kproma) =  &
          pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) - pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev))
     ! Calculate reference density at the full model levels(kg/m3)
     rho0_ll(1:kproma) = grmass(_RI_XYZ__(1:kproma,jrow,nlev)) / grvol(_RI_XYZ__(1:kproma,jrow,nlev))
     
     CALL calc_slp( slp(1:kproma,jrow), aps(1:kproma,jrow), &
          tm1(_RI_XYZ__(1:kproma,jrow,nlev)),rho0_ll(1:kproma), dp0_ll(1:kproma),  &
          altitudei_msl(_RI_XYZ__(1:kproma,jrow,nlev+1)), kproma, g, rd )
  END IF
  
  IF(l_cpt)THEN

     DO jp=1, kproma

        CALL cptdiag(tm1_3d(_RI_XYZ__(jp,jrow,:)), press_3d(_RI_XYZ__(jp,jrow,:)) &  ! IN
                     , tp(jp,jrow)                            &
                     , cpt(jp,jrow)                           &  ! OUT
                     , i_cpt_i )

        cpt_i(jp,jrow) = REAL(i_cpt_i, DP)

     END DO

  END IF
  
  IF (l_windspeed) THEN
     windspeed(_RI_XYZ__(1:kproma,jrow,:)) = &
          calc_windspeed(um1(_RI_XYZ__(1:kproma,jrow,:)) &
          , vm1(_RI_XYZ__(1:kproma,jrow,:)))
  END IF

END SUBROUTINE tropop_vdiff
! -------------------------------------------------------------------------

! **************************************************************************
END MODULE messy_tropop_si
! **************************************************************************
