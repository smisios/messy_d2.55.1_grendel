!******************************************************************************
! SMIL INTERFACE TO ECHAM5 FOR SURFACE PROCESSES: former routines:
!    surf.f90 lake.f90 icetemp.f90 soiltemp.f90 sicetemp.f90
!
! Authors: see subroutines
!          Sabine Brinkop, 21.9.2012 modularization of surface routines
!
! *****************************************************************************
MODULE messy_surface_si
! *****************************************************************************

#if defined(ECHAM5) || defined (CESM1)
 
  ! BMIL
  USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi &
                                     , error_bi, warning_bi

  USE messy_main_grid_def_mem_bi,ONLY: ngl
  USE messy_main_data_bi,        ONLY: eps,            &
                                       lcouple,        &
                                       wsmx,ws,wl,     & 
                                       tsl,tslm,tslm1, &
#ifdef CESM1
                                       tslnew,         & 
#endif
                                       u10, v10,       & 
                                       snmel,snc,      &
                                       sn,sni,snacl,   &
                                       gld, runoff,    &
                                       rogl, drain,    &
                                       apmegl, orostd, &
                                       rgcgn, grndcapc,&
                                       grndhflx,       &
                                       grndflux,       &
                                       tsoil,          &
                                       grndc, grndd,   &
                                       sodif, glac,    &
                                       ahfres, ahfcon, &
                                       loland_2d,      &
                                       loglac_2d,      &
                                       cvsc,           &
                                       evwsd
  USE messy_main_timer,         ONLY: delta_time, lstart 
  ! SMCL
  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_surface

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: surface_initialize
  PUBLIC :: surface_init_memory
  PUBLIC :: surface_init_coupling
  PUBLIC :: surface_global_start
#ifdef CESM1
  PUBLIC :: surface_local_start
#endif
  PUBLIC :: surface_radiation
  PUBLIC :: surface_mixlo
  PUBLIC :: surface_global_end
  !PRIVATE :: surface_read_nml_cpl

#ifdef ECHAM5
  CHARACTER(LEN=*), PARAMETER :: basemod = 'ECHAM5'
#endif
#ifdef CESM1
  CHARACTER(LEN=*), PARAMETER :: basemod = 'CESM1'
#endif

  ! POINTERS To CHANNEL OBJECTS
  REAL(dp), DIMENSION(:,:), POINTER :: slm       => NULL() ! land mask
  REAL(dp), DIMENSION(:,:), POINTER :: slf       => NULL() ! land fraction
  REAL(dp), DIMENSION(:,:), POINTER :: seaice    => NULL() ! seaice fraction rel to ocean
  REAL(dp), DIMENSION(:,:), POINTER :: siced     => NULL() ! ice depth in m
  REAL(dp), DIMENSION(:,:), POINTER :: icecov    => NULL() ! ice cover (fraction of grid box)
  REAL(dp), DIMENSION(:,:), POINTER :: seacov    => NULL() ! sea cover (fraction of grid box)
  REAL(dp), DIMENSION(:,:), POINTER :: landcov   => NULL() ! land cover (fraction of grid box)
  REAL(dp), DIMENSION(:,:), POINTER :: tsi       => NULL() ! surface temperature of ice in K
  REAL(dp), DIMENSION(:,:), POINTER :: tsw       => NULL() ! surface temperature of water in K
  REAL(dp), DIMENSION(:,:), POINTER :: alake     => NULL() ! lake fraction of grid box
  REAL(dp), DIMENSION(:,:), POINTER :: fluxres   => NULL() ! heat flux residual

  REAL(dp), DIMENSION(:,:), POINTER :: ahflw     => NULL() ! latent heat flux over water in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: ahfsw     => NULL() ! sensible heat flux over water in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: ahfli     => NULL() ! latent heat flux over ice in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: ahfsi     => NULL() ! sensible heat flux over ice in W m-2

  REAL(dp), DIMENSION(:,:), POINTER :: trflw     => NULL() ! LW flux over water in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: trfli     => NULL() ! LW flux over ice in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: soflw     => NULL() ! SW flux over water in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: sofli     => NULL() ! SW flux over ice in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: ahfice    => NULL() ! conductive heat flux in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: evapi     => NULL() ! evaporation over ice

  REAL(dp), DIMENSION(:,:), POINTER :: cvsi      => NULL() ! snow cover over ice (fraction of grid box)
  REAL(dp), DIMENSION(:,:), POINTER :: qres      => NULL() ! res. heat flux for melting sea ice in W m-2
  REAL(dp), DIMENSION(:,:,:),POINTER :: aphm1_3d => NULL() !

  REAL(dp), DIMENSION(:,:), POINTER :: cvs         => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: cvw         => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: rsfl_2d     => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: rsfc_2d     => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: ssfl_2d     => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: ssfc_2d     => NULL() !

  REAL(dp), DIMENSION(:,:), POINTER :: evapl_2d    => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: evapot_2d   => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: aros_2d     => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: apmecal_2d  => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER :: adrain_2d   => NULL() !

  REAL(dp), DIMENSION(:,:), POINTER :: tsurf_2d   => NULL() !

  REAL(dp), POINTER, DIMENSION(:) :: pevwsd! local variable converted to pointer

  ! own channel objects
  REAL(dp), DIMENSION(:,:), POINTER :: wlmx        => NULL() ! ub_ak_20190206

#if defined(CESM1)
  ! define ncregrid event triggers
  REAL(dp), POINTER, DIMENSION(:,:) :: faom      => NULL() ! FAO data set (soil data flags 0...5.)
  REAL(dp), POINTER, DIMENSION(:,:)   :: import_ws => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: import_sn => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: import_tsoil => NULL() 
#endif

  ! ub_ak_20190204+
  TYPE(t_chaobj_cpl)                :: imp_lai
  REAL(dp), POINTER, DIMENSION(:,:) :: lai => NULL() 
  ! ub_ak_20190204-

  LOGICAL :: L_LGMC = .FALSE.      ! submodel LGMC is ON
  LOGICAL :: L_MLOCEAN = .FALSE.   ! submodel MLOCEAN is ON

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

CONTAINS

! ****************************************************************************

!---------------------------------------------------------------------------
  SUBROUTINE surface_initialize

    ! BMIL
    USE messy_main_mpi_bi,         ONLY: p_parallel_io, p_io, p_bcast
    ! SMCL
    USE messy_main_tools,          ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'surface_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    CALL start_message_bi(modstr,'surface_initialize',substr)

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL surface_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('error in surface_read_nml_ctrl', substr)
    END IF

    CALL p_bcast(lice,  p_io)

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL surface_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(imp_lai%CHA,    p_io)
    CALL p_bcast(imp_lai%OBJ,    p_io)

    CALL end_message_bi(modstr,'surface_initialize',substr)

  END SUBROUTINE surface_initialize
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE surface_init_memory

    ! BNIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL
    ! SMCL
    USE messy_main_channel,          ONLY: new_channel, new_channel_object

#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,      ONLY: mtend_get_handle,       &
                                           mtend_register,         &
                                           mtend_id_t
#endif

    IMPLICIT NONE
  
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'surface_init_memory'
    INTEGER :: status

    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)

    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'wlmx', p2=wlmx &
         , reprid = GP_2D_HORIZONTAL, lrestreq=.TRUE. )

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)

#ifdef MESSYTENDENCY
     my_handle=mtend_get_handle(modstr)
     CALL mtend_register(my_handle, mtend_id_t)
#endif

    END SUBROUTINE surface_init_memory
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE surface_init_coupling

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! SMCL
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info
#ifdef CESM1
    USE messy_main_grid_def_mem_bi,  ONLY: ngpblks, npromz
#endif

    IMPLICIT NONE
    INTRINSIC :: TRIM

    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr='surface_init_coupling'
    CHARACTER(LEN=12) :: radcha = ''
    CHARACTER(LEN=12) :: vdicha = ''

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)


#ifdef ECHAM5
    CALL get_channel_info(status, 'e5vdiff')
    IF (status /= 0) THEN
       CALL warning_bi( &
            'channel e5vdiff not available, trying vertex ...' &
            , substr)
       CALL get_channel_info(status, 'vertex')
       IF (status /= 0) THEN
          CALL error_bi(' ... vertex also not available!', substr)
       ELSE
          vdicha = 'vertex'
       ENDIF
    ELSE
       vdicha = 'e5vdiff'
    ENDIF
#endif
#ifdef CESM1
    CALL get_channel_info(status, 'vertdiff')
    IF (status /= 0) THEN
       CALL warning_bi( &
            'channel vertdiff not available, trying e5vdiff ...' &
            , substr)
       CALL get_channel_info(status, 'e5vdiff')
       IF (status /= 0) THEN
          CALL error_bi(' ... vertex e5vdiff not available!', substr)
       ELSE
          vdicha = 'e5vdiff'
       ENDIF
    ELSE
       vdicha = 'vertdiff'
    ENDIF
#endif
    
    ! can currently not be moved from g3b to here (ioinitial.f90)
    CALL get_channel_object(status,basemod,'slm', p2=slm)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'slf', p2=slf)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'seaice', p2=seaice)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'siced', p2=siced)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'icecov', p2=icecov)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,basemod,'seacov', p2=seacov)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,basemod,'landcov', p2=landcov)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'tsi', p2=tsi)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'tsw', p2=tsw)
    CALL channel_halt(substr, status)
  
#ifdef ECHAM5
    ! can currently not be moved from g3b to here (ioinitial.f90)
    CALL get_channel_object(status,basemod,'alake', p2=alake)
#endif
#ifdef CESM1
    CALL get_channel_object(status,'import_grid', 'SURFACEDAT_ALAKE', p2=alake)
#endif  
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,basemod,'fluxres', p2=fluxres)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,TRIM(vdicha),'ahflw', p2=ahflw)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,TRIM(vdicha),'ahfsw', p2=ahfsw)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,TRIM(vdicha),'ahfli', p2=ahfli)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,TRIM(vdicha),'ahfsi', p2=ahfsi)
    CALL channel_halt(substr, status)

    CALL get_channel_info(status, 'rad01')
    IF (status /=0) THEN
       CALL warning_bi( &
            'channel rad01 not available, trying old rad4all ...' &
            , substr)
       CALL get_channel_info(status, 'rad4all')
       IF (status /= 0) THEN
          CALL error_bi(' ... old rad4all also not available!', substr)
       ELSE
          radcha = 'rad4all'
       ENDIF
    ELSE
       radcha = 'rad01'
    ENDIF

    CALL get_channel_object(status, TRIM(radcha), 'trflw', p2=trflw)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status, TRIM(radcha), 'trfli', p2=trfli)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status, TRIM(radcha), 'soflw', p2=soflw)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status, TRIM(radcha), 'sofli', p2=sofli)
    CALL channel_halt(substr, status)
!     
    CALL get_channel_object(status,basemod,'ahfice', p2=ahfice)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,TRIM(vdicha),'evapi', p2=evapi)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'cvsi', p2=cvsi)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'qres', p2=qres)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,basemod,'pressi', p3=aphm1_3d)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'cvs', p2=cvs)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'cvw', p2=cvw)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,TRIM(vdicha),'evapot_2d', p2=evapot_2d)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'rsfl_2d', p2=rsfl_2d)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'rsfc_2d', p2=rsfc_2d)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'ssfl_2d', p2=ssfl_2d)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'ssfc_2d', p2=ssfc_2d)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,TRIM(imp_lai%cha),TRIM(imp_lai%obj) &
         , p2=lai)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,TRIM(vdicha),'evapl_2d', p2=evapl_2d)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'aros', p2=aros_2d)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'adrain', p2=adrain_2d)
    CALL channel_halt(substr, status)
  
    CALL get_channel_object(status,basemod,'apmecal', p2=apmecal_2d)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,basemod,'tsurf', p2=tsurf_2d)
    CALL channel_halt(substr, status)

    ! if LGMC is running, MLOCEAN must be called in global_end instead of physc.
    CALL get_channel_info(status, 'lgmc')
    L_LGMC = (status == 0)
    CALL get_channel_info(status, 'mlocean')
    L_MLOCEAN = (status == 0)
     
#ifdef CESM1
    CALL get_channel_object(status, 'import_grid', 'SURFACEDATIDX_FAO', p2=faom)
    CALL channel_halt(substr//': object FAO in channel import_grid not found!',&
         status)
    CALL get_channel_object(status, 'import_grid', 'SURFACEDAT_WS', p2=import_ws)
    CALL channel_halt(substr//': object WS in channel import_grid not found!',&
         status)    
    CALL get_channel_object(status, 'import_grid', 'SURFACEDAT_SN', p2=import_sn)
    CALL channel_halt(substr//': object SN in channel import_grid not found!',&
         status)
    CALL get_channel_object(status, 'import_grid', 'tsoil_tsoil', p3=import_tsoil)
    CALL channel_halt(substr//': object tsoil in channel import_grid not found!',&
         status)
#endif

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)    

  END SUBROUTINE surface_init_coupling
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE surface_global_start
#ifdef CESM1
    USE messy_main_grid_def_bi,     ONLY: ilat, ilon
    USE messy_main_grid_def_mem_bi, ONLY: npromz, nproma, ngpblks        
    USE messy_main_data_bi,       ONLY: tslclim, sst, seaice &
                                      , siced, tsw, tsi, alake
    USE messy_main_timer,         ONLY: lstart
    USE messy_main_constants_mem, ONLY: tmelt, ctfreez, api=>pi

    CHARACTER(LEN=*), PARAMETER :: substr = 'surface_global_start'
    INTEGER                     :: jp, jk, jl, zkproma, zjrow
    LOGICAL                     :: lok     ! OK?
    REAL(dp) :: ztsl(nproma), zic(nproma), zts(nproma), zdelti


    IF (lstart) THEN
       ! from: echam5/src/ioinitial.f90 ! mz_ab_20160803
       ! Setting of array of variable soil characteristics to be use
       ! in *surf*
       ! Input: FAO soils interpolated from 0.5 degree resolution
       !        to model resolution (simple average).
       ! and setting of array of variable available water storage capacity
       ! to be used in *surf*
       ! Input: Patterson data interpolated from 0.5 degree resolution
       !        to model resolution (simple average).
       DO zjrow =1,ngpblks
          zkproma = npromz(zjrow)
          DO jl=1,zkproma
             IF (NINT(faom(jl,zjrow)) == 1) THEN
                rgcgn(jl,zjrow) = 1.93e+06_dp
                sodif(jl,zjrow) = 8.7e-7_dp
             ELSE IF (NINT(faom(jl,zjrow)) == 2) THEN
                rgcgn(jl,zjrow) = 2.10e+06_dp
                sodif(jl,zjrow) = 8.0e-7_dp
             ELSE IF (NINT(faom(jl,zjrow)) == 3) THEN
                rgcgn(jl,zjrow) = 2.25e+06_dp
                sodif(jl,zjrow) = 7.4e-7_dp
             ELSE IF (NINT(faom(jl,zjrow)) == 4) THEN
                rgcgn(jl,zjrow) = 2.36e+06_dp
                sodif(jl,zjrow) = 7.1e-7_dp
             ELSE IF (NINT(faom(jl,zjrow)) == 5) THEN
                rgcgn(jl,zjrow) = 2.48e+06_dp
                sodif(jl,zjrow) = 6.7e-7_dp
             ELSE
                IF (NINT(faom(jl,zjrow)) == 0) THEN
                   rgcgn(jl,zjrow) = 2.25e+06_dp
                   sodif(jl,zjrow) = 7.4e-7_dp
                ELSE
                   CALL error_bi(' error in FAO data input!',substr)
                   !WRITE (nerr,*) 'faom(',jl,',',zjrow,') = ',faom(jl,zjrow)
                END IF
             END IF
             ws(jl,zjrow) = MIN(import_ws(jl,zjrow),wsmx(jl,zjrow))
          END DO
       END DO
       sn(:,:) = import_sn(:,:)
    ENDIF ! lstart

    ! tsoil is read only for initialisation in first timestep
    ! see echam5/src/initemp.f90 where tsoil is initialised online
    IF (lstart) THEN
       tsoil(:,:,:) = import_tsoil(:,:,:)
       !
       DO zjrow =1,ngpblks

          ! from: echam5/src/initemp.f90
          !  Set all time levels of surface temperature to uppermost soil temp.
          DO jp = 1,nproma
             tsl(jp,zjrow)   = tsoil(jp,1,zjrow)
             tslm(jp,zjrow)  = tsl(jp,zjrow)
             tslm1(jp,zjrow) = tsl(jp,zjrow)
             tslnew(jp,zjrow) = tsl(jp,zjrow)
          END DO
       END DO
       
       ! from: echam5/src/initemp.f90
       DO zjrow =1,ngpblks
          zkproma = npromz(zjrow)
          
          ztsl(1:nproma)=tslclim(1:nproma,zjrow)
          
          DO jp=1,nproma
             zts(jp)=sst(jp,zjrow) ! from CESM
             zic(jp)=seaice(jp,zjrow) ! from CESM
             IF(zic(jp).LE.0.01_dp) zic(jp)=0._dp
          END DO
          !
          !
          !   Initialize temperatures and ice
          !
          !
          DO jp = 1,nproma
             IF (alake(jp,zjrow).GE.0.5_dp) THEN         !  lakes
                tsi(jp,zjrow)=ztsl(jp)
                zdelti=tsi(jp,zjrow)-tmelt
                IF (zdelti.LT.0._dp) THEN
                   siced(jp,zjrow)=1._dp-EXP(-0.005_dp*zdelti**2)
                ELSE
                   siced(jp,zjrow)=0._dp
                END IF
                IF (siced(jp,zjrow).GE.0.1_dp) THEN
                   seaice(jp,zjrow)=1._dp
                   tsw(jp,zjrow)=tmelt
                ELSE
                   siced(jp,zjrow)=0._dp
                   seaice(jp,zjrow)=0._dp
                   tsw(jp,zjrow)=MAX(ztsl(jp),tmelt)
                   tsi(jp,zjrow)=MIN(ztsl(jp),tmelt)
                END IF
             ELSE
                siced(jp,zjrow)=0._dp
                seaice(jp,zjrow)=0._dp
                tsw(jp,zjrow)=MAX(ztsl(jp),tmelt)
                tsi(jp,zjrow)=MIN(ztsl(jp),tmelt)
             END IF
             IF (alake(jp,zjrow).EQ.0.0_dp .AND. slf(jp,zjrow).LE.0.5_dp) THEN          !  ocean
                seaice(jp,zjrow)=MAX(0._dp,MIN(0.99_dp,zic(jp)))
                IF (seaice(jp,zjrow).GT.0._dp) THEN
                   tsi(jp,zjrow)=MIN(ztsl(jp),tmelt)
                   tsw(jp,zjrow)=ctfreez
                   zdelti=tsi(jp,zjrow)-tmelt
                   IF (zdelti.LT.0._dp) THEN
                      siced(jp,zjrow)=2._dp*(1._dp-EXP(-0.005_dp*zdelti**2))+0.1_dp
                   ELSE
                      siced(jp,zjrow)=0.1_dp
                   END IF
                ELSE
                   tsi(jp,zjrow)=tmelt
                   tsw(jp,zjrow)=MAX(ctfreez,zts(jp))
                   siced(jp,zjrow)=0.0_dp
                END IF
             END IF
          END DO
       END DO
    END IF ! lstart
#endif

  END SUBROUTINE surface_global_start
!---------------------------------------------------------------------------

#ifdef CESM1
  SUBROUTINE surface_local_start

    IMPLICIT NONE

    CALL surface_radiation

  END SUBROUTINE surface_local_start
#endif

!---------------------------------------------------------------------------
  SUBROUTINE surface_radiation

    USE messy_main_grid_def_mem_bi, ONLY: jrow, nproma
    USE messy_main_data_bi,       ONLY: lcouple
    USE messy_main_constants_mem, ONLY: tmelt, cwlmax

    IMPLICIT NONE

    INTEGER             :: jp
    REAL(dp)            :: zsn_mm, zsigh
    REAL(dp), PARAMETER :: zepsec  = 1.E-12_dp
    REAL(dp), PARAMETER :: zsigfac = 0.15_dp
    REAL(dp), PARAMETER :: cqsncr = 0.95_dp   

    !*        3.6   COMPUTE LOGICAL MASK FOR LAND AND GLACIER.
    !
    loland_2d(:,jrow) = .FALSE.
    loglac_2d(:,jrow) = .FALSE.
    DO jp=1,nproma
       loland_2d(jp,jrow) = slm(jp,jrow).GT.0._dp
       loglac_2d(jp,jrow) = loland_2d(jp,jrow) .AND. glac(jp,jrow).GT.0._dp
    END DO

#ifndef CESM1
    ! CESM uses its own maps, therefore this block needs to be omitted

    !       3.7 Weighting factors for fractional surface coverage
    !           Accumulate ice portion for diagnostics
    !
    DO jp=1,nproma
       landcov(jp,jrow) = slm(jp,jrow)
       seacov(jp,jrow)  = (1._dp-slm(jp,jrow))*(1._dp-seaice(jp,jrow))
       icecov(jp,jrow)  = 1._dp-landcov(jp,jrow)-seacov(jp,jrow) 
    END DO

    IF (lcouple) THEN
       DO jp=1,nproma
          IF(slf(jp,jrow).GT.1.0_dp-zepsec) THEN
             tsi(jp,jrow)=tmelt
             tsw(jp,jrow)=tmelt
          END IF
       END DO
    ENDIF
#else
    WHERE ((.NOT.loland_2d(:,jrow)).OR.loglac_2d(:,jrow))
       ws(:,jrow) = 0._dp
    END WHERE
#endif
    !
    !      3.8  Skin reservoir, wet skin fraction and snow cover
    !           (bare land, canopy, lake ice)
    !
    wlmx(:,jrow) = 0._dp
    DO jp=1,nproma
       IF (.NOT.loglac_2d(jp,jrow)) THEN
          wlmx(jp,jrow) = cwlmax*(1._dp+lai(jp,jrow))
          cvw(jp,jrow)=MIN(wl(jp,jrow)/wlmx(jp,jrow),1.0_dp)
          zsn_mm=1000._dp*sn(jp,jrow)
          zsigh=SQRT(zsn_mm/(zsn_mm+zepsec+zsigfac*orostd(jp,jrow)))
          cvs(jp,jrow)=cqsncr*TANH(zsn_mm/10._dp)*zsigh
          cvsc(jp,jrow)=MIN(1._dp,snc(jp,jrow)/(wlmx(jp,jrow) &
               - cwlmax+EPSILON(1._dp)))
          IF (cvs(jp,jrow).LT.EPSILON(1._dp) .AND. &
               cvsc(jp,jrow).GE.EPSILON(1._dp)) THEN
             cvs(jp,jrow)=cvsc(jp,jrow)
          END IF
       ELSE
          cvw(jp,jrow)=0._dp
          cvs(jp,jrow)=1._dp
          cvsc(jp,jrow)=0._dp
       END IF
       IF (.NOT. loland_2d(jp,jrow)) THEN
          cvsi(jp,jrow)=TANH(sni(jp,jrow)*100._dp)
       ELSE
          cvsi(jp,jrow)=0._dp
       END IF
    END DO

  END SUBROUTINE surface_radiation
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE surface_mixlo(flag)

    USE messy_main_grid_def_mem_bi,       ONLY: nlev, npromz   &
         ,nproma, ngpblks        & 
         ,jrow 
    USE messy_main_data_bi,       ONLY: tm1              &
         ,qm1, aphm1             &
         ,snc_surf, gld_surf     !op_re_20140904

#ifndef MESSYTENDENCY
#ifdef ECHAM5
    USE messy_main_data_bi,       ONLY: tte_scb
#endif
#ifdef CESM1
    USE messy_main_data_bi,       ONLY: tte_scb=>tte_3d
#endif
#endif

#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,   ONLY: mtend_get_start_l,      &
                                        mtend_add_l,            &
                                        mtend_id_t
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    INTEGER :: zproma
    REAL(DP), DIMENSION(:,:), POINTER :: ztte

#ifdef MESSYTENDENCY
     REAL(dp),DIMENSION(nproma,nlev) :: lo_tte,lo_tm1
#endif

     REAL(dp), POINTER, DIMENSION(:) :: psnc_surf !new pointer for temp result
     REAL(dp), POINTER, DIMENSION(:) :: pgld_surf !new pointer for temp result

     psnc_surf => snc_surf(:,jrow)
     pgld_surf => gld_surf(:,jrow)
     pevwsd    => evwsd(:,jrow)

      IF (L_LGMC) RETURN

#ifdef MESSYTENDENCY
     lo_tte = 0.0_dp
     CALL mtend_get_start_l (mtend_id_t, v0 = lo_tm1)   ! t vorläufig neu
#endif 

#ifdef ECHAM5
     if ( jrow == ngpblks) then
        zproma = npromz
     else
        zproma = nproma
     endif
#endif
#ifdef CESM1
     zproma = npromz(jrow)
#endif

    SELECT CASE(flag)

    CASE(1)

       ALLOCATE(ztte(nproma,nlev))

       ztte(:,:) = 0._dp

       ! Temporary results needed for h2oiso
       psnc_surf(:) = snc(:,jrow)!op_re_20140904
       pgld_surf(:) = gld(:,jrow)!op_re_20140904

#ifdef CESM1
       ! accumulated channel objects (see below):
       grndflux(:,jrow) = 0._dp
#endif
       CALL surface ( zproma,    nproma,            nlev              &
            , ngl,               delta_time,        eps, lstart       &
            , tsl(:,jrow),       tslm(:,jrow),      tslm1(:,jrow)     &  
            , ws(:,jrow),        wl(:,jrow),        wsmx(:,jrow)      &  
            , sn(:,jrow),        snmel(:,jrow),     gld(:,jrow)       &  
            , snc(:,jrow),       u10(:,jrow),       v10(:,jrow)       &  
            , runoff(:,jrow),    rogl(:,jrow),      drain(:,jrow)     &  
            , apmegl(:,jrow),    snacl(:,jrow),     orostd(:,jrow)    &  
            , rgcgn(:,jrow),     sodif(:,jrow),     slm(:,jrow)       &
            , grndcapc(:,jrow),  grndhflx(:,jrow),  grndflux(:,jrow)  &
            , tsoil(:,:,jrow),   grndd(:,:,jrow),   grndc(:,:,jrow)   &
            , tm1(:,:,jrow),     qm1(:,:,jrow),     ztte(:,:)         &
            , aphm1                                                   &
            , cvs(:,jrow),       cvw(:,jrow),       wlmx(:,jrow)      & 
            , evapl_2d(:,jrow),  evapot_2d(:,jrow)                    & 
            , rsfl_2d(:,jrow),   rsfc_2d(:,jrow)                      & 
            , ssfl_2d(:,jrow),   ssfc_2d(:,jrow)                      &
            , aros_2d(:,jrow),   adrain_2d(:,jrow)                    & 
            , apmecal_2d(:,jrow)                                      & 
            , loland_2d(:,jrow), loglac_2d(:,jrow)                    &
            , pevwsd(:)                                               )

#ifndef MESSYTENDENCY
            tte_scb(1:zproma,:,jrow) = tte_scb(1:zproma,:,jrow) +     &
                                       ztte(1:zproma,:)
#else
            lo_tte(1:zproma,:) = ztte(1:zproma,:)
            CALL mtend_add_l (my_handle, mtend_id_t, px = lo_tte)
#endif

       CALL lake ( zproma, delta_time                                    &
            , seaice(:,jrow),    siced(:,jrow),    alake(:,jrow)      &
            , tsi(:,jrow),       tsw(:,jrow)                          &
            , ahflw(:,jrow),     ahfsw(:,jrow),    fluxres(:,jrow)    &
            , trflw(:,jrow),     soflw(:,jrow)                        &
            , evapi(:,jrow),     sni(:,jrow),      cvsi(:,jrow)       &
            , ahfres(:,jrow),    icecov(:,jrow)           )

       CALL licetemp ( zproma, delta_time                                &
            , siced(:,jrow),     sni(:,jrow),      alake(:,jrow)      &
            , tsi(:,jrow),       trfli(:,jrow),    sofli(:,jrow)      &
            , ahfice(:,jrow),    fluxres(:,jrow)                      &
            , ahfcon(:,jrow),    ahfres(:,jrow),   evapi(:,jrow)      &
            , ssfl_2d(:,jrow),   ssfc_2d(:,jrow)                      &
            , ahfsi(:,jrow),     ahfli(:,jrow),    cvsi(:,jrow)       &
            , icecov(:,jrow)     )

#ifdef CESM1
       ! for CESM, averaging is done through CHANNEL 
       ! instead of the echam laccu stream objects routines.
       ! Therefore, accumulation in smcl routines is wrong and corrected here:
       grndflux(:,jrow) = grndflux(:,jrow)/delta_time
#endif
       ! release space
       DEALLOCATE(ztte) 

    CASE (2)                      ! Tendenzen fehlen noch.

#ifdef CESM1
       ahfcon(:,jrow) = 0._dp
       ahfres(:,jrow) = 0._dp
#endif
       CALL sicetemp ( zproma, delta_time, lcouple, L_MLOCEAN         &
            , siced(:,jrow),     sni(:,jrow),      alake(:,jrow)      &
            , slf(:,jrow)                                             &
            , tsi(:,jrow),       trfli(:,jrow),    sofli(:,jrow)      &
            , ahfice(:,jrow),    fluxres(:,jrow),  qres(:,jrow)       &
            , ahfcon(:,jrow),    ahfres(:,jrow)                       &
            , ahfsi(:,jrow),     ahfli(:,jrow)                        &
            , icecov(:,jrow) )

#ifdef CESM1
       ! for CESM, averaging is done through CHANNEL 
       ! instead of the echam laccu stream objects routines.
       ! Therefore, accumulation in smcl routines is wrong and corrected here:
       ahfcon(:,jrow) = ahfcon(:,jrow)/delta_time
       ahfres(:,jrow) = ahfcon(:,jrow)/delta_time
#endif
       tsurf_2d(:,jrow) = landcov(:,jrow) * tslm1(:,jrow) &
                        + icecov(:,jrow)  * tsi(:,jrow)   &
                        + seacov(:,jrow)  * tsw(:,jrow)

    END SELECT

  END SUBROUTINE surface_mixlo
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE surface_global_end(flag)

    USE messy_main_grid_def_mem_bi, ONLY:  npromz &
                                        , nproma, nlev, ngpblks   &
                                        , nlev
    USE messy_main_data_bi,         ONLY: qm1, tm1
#ifndef MESSYTENDENCY
#ifdef ECHAM5
    USE messy_main_data_bi,         ONLY:  tte_scb
#endif
#ifdef CESM1
    USE messy_main_data_bi,         ONLY:  tte_scb=>tte_3d
#endif
#endif

#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,   ONLY:  &
         mtend_get_start_l,              &
         mtend_add_l,                    &
         mtend_id_t
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    INTEGER :: zproma, krow
    REAL(DP), DIMENSION(:,:), POINTER :: ztte

#ifdef MESSYTENDENCY
    REAL(dp), DIMENSION(nproma,nlev)    :: lo_tte,lo_tm1
#endif

    IF (.NOT. L_LGMC) RETURN

#ifdef MESSYTENDENCY       
     lo_tte = 0._dp 
     CALL mtend_get_start_l(mtend_id_t, v0 = lo_tm1)   
#endif

    SELECT CASE(flag)

    CASE(1)

       ALLOCATE(ztte(nproma,nlev))

       ztte(:,:) = 0._dp

       do  krow = 1, ngpblks 
#ifdef ECHAM5 
          if ( krow == ngpblks) then
             zproma = npromz
          else
             zproma = nproma
          endif
#endif
#ifdef CESM1
          zproma = npromz(krow) 
#endif

          CALL surface ( zproma,    nproma,            nlev              &
               , ngl,               delta_time,        eps, lstart       &
               , tsl(:,krow),       tslm(:,krow),      tslm1(:,krow)     &
               , ws(:,krow),        wl(:,krow),        wsmx(:,krow)      &
               , sn(:,krow),        snmel(:,krow),     gld(:,krow)       &
               , snc(:,krow),       u10(:,krow),       v10(:,krow)       &
               , runoff(:,krow),    rogl(:,krow),      drain(:,krow)     &
               , apmegl(:,krow),    snacl(:,krow),     orostd(:,krow)    &
               , rgcgn(:,krow),     sodif(:,krow),     slm(:,krow)       &
               , grndcapc(:,krow),  grndhflx(:,krow),  grndflux(:,krow)  &
               , tsoil(:,:,krow),   grndd(:,:,krow),   grndc(:,:,krow)   &
               , tm1(:,:,krow),     qm1(:,:,krow),     ztte(:,:)         &
               , aphm1_3d(:,:,krow)                                      &  
               , cvs(:,krow),       cvw(:,krow),       wlmx(:,krow)      &
               , evapl_2d(:,krow),  evapot_2d(:,krow)                    &
               , rsfl_2d(:,krow),   rsfc_2d(:,krow)                      & 
               , ssfl_2d(:,krow),   ssfc_2d(:,krow)                      &
               , aros_2d(:,krow),   adrain_2d(:,krow)                    &
               , apmecal_2d(:,krow)                                      &
               , loland_2d(:,krow), loglac_2d(:,krow) &
               , pevwsd(:)                                               )

#ifndef MESSYTENDENCY
          tte_scb(1:zproma,:,krow) = tte_scb(1:zproma,:,krow) + ztte(1:zproma,:)
#else   
          lo_tte(1:zproma,:) = ztte(1:zproma,:)
          CALL mtend_add_l(my_handle, mtend_id_t, px = lo_tte)
#endif                          


          CALL lake ( zproma, delta_time                                 &
               , seaice(:,krow),    siced(:,krow),    alake(:,krow)      &
               , tsi(:,krow),       tsw(:,krow)                          &
               , ahflw(:,krow),     ahfsw(:,krow),    fluxres(:,krow)    &
               , trflw(:,krow),     soflw(:,krow)                        &
               , evapi(:,krow),     sni(:,krow),      cvsi(:,krow)       &
               , ahfres(:,krow),    icecov(:,krow)           )

          CALL licetemp ( zproma, delta_time                             &
               , siced(:,krow),     sni(:,krow),      alake(:,krow)      &
               , tsi(:,krow),       trfli(:,krow),    sofli(:,krow)      &
               , ahfice(:,krow),    fluxres(:,krow)                      &
               , ahfcon(:,krow),    ahfres(:,krow),   evapi(:,krow)      &
               , ssfl_2d(:,krow),   ssfc_2d(:,krow)                      &
               , ahfsi(:,krow),     ahfli(:,krow),    cvsi(:,krow)       &
               , icecov(:,krow)                            )

       ENDDO  ! krow loop

       ! release space
       DEALLOCATE(ztte)     

    CASE (2)   

       do krow = 1, ngpblks   
#ifdef ECHAM5
          if ( krow == ngpblks) then
             zproma = npromz
          else
             zproma = nproma
          endif
#endif
#ifdef CESM1
          zproma = npromz(krow)
#endif

          CALL sicetemp ( zproma, delta_time, lcouple, L_MLOCEAN         &
               , siced(:,krow),     sni(:,krow),      alake(:,krow)      &
               , slf(:,krow)                                             &
               , tsi(:,krow),       trfli(:,krow),    sofli(:,krow)      &
               , ahfice(:,krow),    fluxres(:,krow),  qres(:,krow)       &
               , ahfcon(:,krow),    ahfres(:,krow)                       &
               , ahfsi(:,krow),     ahfli(:,krow)                        &
               , icecov(:,krow) )

       ENDDO  ! krow loop

    END SELECT

  END SUBROUTINE surface_global_end
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE surface_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' 
    
    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ imp_lai

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='surface_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE surface_read_nml_cpl
!---------------------------------------------------------------------------

#endif
! ECHAM5 || CESM1

! **********************************************************************
END MODULE messy_surface_si
! **********************************************************************
