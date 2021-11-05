! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL MLOCEAN 
!
! Author : Markus Kunze, FUB, 2009 - 2011
!
! References: see messy_mlocean.f90
!
! ! op_sb_20121205+
! In case LGMC is running, the mixed layer ocean together with the surface
! routines must be called in global_end. This has to be checked.
! ! op_sb_20121205-
! **********************************************************************

! **********************************************************************
MODULE messy_mlocean_e5
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,        ONLY: start_message_bi, end_message_bi

  ! SMCL
  USE messy_mlocean
  USE messy_mlocean_mlo,     modstr_mlo =>    modstr, modver_mlo    => modver
  USE messy_mlocean_plasim,  modstr_plasim => modstr, modver_plasim => modver

  IMPLICIT NONE
  PRIVATE

  ! CPL-NAMELIST PARAMETERS
  LOGICAL            :: l_diag_out
  CHARACTER(len=256) :: fn_mlo_fluxcorr  ! file name of the mixed layer ocean flux correction

  ! References to messy channel objects:
  REAL(dp), DIMENSION(:,:), POINTER :: slm       => NULL() ! land mask
  REAL(dp), DIMENSION(:,:), POINTER :: slf       => NULL() ! land fraction
  REAL(dp), DIMENSION(:,:), POINTER :: sni       => NULL() ! water equivalent of snow on ice
  REAL(dp), DIMENSION(:,:), POINTER :: seaice    => NULL() ! seaice fraction rel to ocean
  REAL(dp), DIMENSION(:,:), POINTER :: siced     => NULL() ! ice depth in m
  REAL(dp), DIMENSION(:,:), POINTER :: icecov    => NULL() ! ice cover (fraction of grid box)
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
  REAL(dp), DIMENSION(:,:), POINTER :: amlcorr   => NULL() ! flux correction for mixed layer ocean
  REAL(dp), DIMENSION(:,:), POINTER :: amlheat   => NULL() ! ??? ! op_pj_20160617
  REAL(dp), DIMENSION(:,:), POINTER :: evapi     => NULL() ! evaporation over ice
  REAL(dp), DIMENSION(:,:), POINTER :: cvsi      => NULL() ! snow cover over ice (fraction of grid box)
  REAL(dp), DIMENSION(:,:), POINTER :: qres      => NULL() ! residual heat flux for melting sea ice in W m-2
  REAL(dp), DIMENSION(:,:), POINTER :: fsnet     => NULL() ! net surface heat flux
  
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: aflux ! (nlon,ngl,0:13) in global coordinates

  LOGICAL :: L_LGMC = .FALSE.      ! op_sb_20120928 ! submodel LGMC is ON

  ! PUBLIC SUBROUTINES
  PUBLIC :: mlocean_initialize    ! initialize submodel
  PUBLIC :: mlocean_init_memory   ! request memory
  PUBLIC :: mlocean_init_coupling ! set pointers for coupling to BM and other SMs
  PUBLIC :: mlocean_global_start  ! entry point in time loop (all vectors)
  PUBLIC :: mlocean_global_end    ! op_sb_20120921: required, if LGMC is on
  
! op_pj_20110224+
  !PRIVATE :: mlocean_mlayerocean   ! entry point in time loop (current vector)
  !PRIVATE :: mlocean_mlflx         ! entry point in time loop (current vector)
  PUBLIC :: mlocean_mixlo           ! entry point in time loop (current vector)
! op_pj_20110224-  

  PUBLIC :: mlocean_free_memory   ! free allocated memory

  ! PRIVATE SUBROTINES
  !PRIVATE :: mlocean_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE mlocean_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast 
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mlocean_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CTRL namelist
    IF (p_parallel_io) THEN                     ! read only on I/O-PE
       iou = find_next_free_unit(100,200)       ! find free I/O unit
       CALL mlocean_read_nml_ctrl(status, iou)  ! read CTRL-namelist
       IF (status /= 0) CALL error_bi(' ',substr)  ! terminate if error
    END IF
    ! BROADCAST CTRL namelist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(mloswitch,   p_io)
    CALL p_bcast(mldmix,      p_io)
    CALL p_bcast(nflxcorr,    p_io)
    CALL p_bcast(flxscale,    p_io)

       CALL p_bcast(fbase_north, p_io)
       CALL p_bcast(fbase_south, p_io)

       CALL p_bcast(ndiag,       p_io)
       CALL p_bcast(nout,        p_io)
       CALL p_bcast(nocean,      p_io)
       CALL p_bcast(newsurf,     p_io)
       CALL p_bcast(ntspd,       p_io)
       CALL p_bcast(nperpetual_ocean, p_io)
       CALL p_bcast(nprint,      p_io)
       CALL p_bcast(nprhor,      p_io)
       CALL p_bcast(dlayer,      p_io)
       CALL p_bcast(taunc,       p_io)
       CALL p_bcast(vdiffk,      p_io)

    ! READ CPL namelist
    IF (p_parallel_io) THEN                    ! read only on I/O-PE
       iou = find_next_free_unit(100,200)      ! find next free I/O unit
       CALL mlocean_read_nml_cpl(status, iou)  ! read CPL-namelist
       IF (status /= 0) CALL error_bi(' ',substr) ! terminate if error
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(l_diag_out,      p_io)
    CALL p_bcast(fn_mlo_fluxcorr, p_io)

    ! ### PERFORM INITIAL SETUP (CALL RESPECTIVE SMCL ROUTINE(S))
    
    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE mlocean_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE mlocean_init_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel,          ONLY: new_channel , new_channel_object &
                                         , new_channel_object_reference &
                                         , new_attribute 
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mlocean_init_memory'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output
    IF (p_parallel_io) WRITE(*,*) 'add new channel mlocean ...'
    
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, 'g3b', 'tsw', modstr, 'tsw')
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, 'g3b', 'tsi', modstr, 'tsi')
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, 'g3b', 'seaice', modstr, 'seaice')
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, 'g3b', 'siced', modstr, 'siced')
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, 'ECHAM5', 'icecov', modstr, 'icecov')
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, 'ECHAM5', 'fluxres', modstr, 'fluxres')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'amlcorr' &
         , p2 = amlcorr, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'amlcorr', 'long_name' &
         , c='mixed layer ocean flux correction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'amlcorr', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    ! op_pj_20160617+
    CALL new_channel_object(status, modstr, 'amlheat' &
         , p2 = amlheat, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    ! op_pj_20160617-

    CALL new_channel_object(status, modstr, 'fsnet' &
         , p2 = fsnet,reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fsnet', 'long_name' &
         , c='net surface fluxes')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fsnet', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    
    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

  END SUBROUTINE mlocean_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE mlocean_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    ! NOTE: Here, part of ECHAM5 is utilised as MESSy BMIL.
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object,   &
                                           get_channel_info
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_pe
    USE messy_main_blather_bi,       ONLY: info_bi, error_bi, warning_bi
    USE messy_main_mpi_bi,           ONLY: dcl, dcg, scatter_gp
    !

    IMPLICIT NONE     ! op_pj_20130408
    INTRINSIC :: TRIM ! op_pj_20130408

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mlocean_init_coupling'
    INTEGER                     :: i
    INTEGER                     :: status   ! error status
    INTEGER                     :: iou      ! I/O unit
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: fluxcorr
    REAL(dp), DIMENSION(:,:,:),             POINTER :: gl_aflux
    CHARACTER(LEN=12) :: radcha = '' ! op_pj_20130408
    CHARACTER(LEN=12) :: vdicha = '' ! op_pj_20160614

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

    ! ### set pointers to channel objects here

! op_pj_20160614+
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
! op_pj_20160614-

    CALL get_channel_object(status,'ECHAM5','slm', p2=slm)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','slf', p2=slf)
    CALL channel_halt(substr, status)
    
    CALL get_channel_object(status,'ECHAM5','sni', p2=sni)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','seaice', p2=seaice)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','siced', p2=siced)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','tsi', p2=tsi)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','tsw', p2=tsw)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','alake', p2=alake)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','fluxres', p2=fluxres)
    CALL channel_halt(substr, status)

!!$ CALL get_channel_object(status,'ECHAM5','ahflw', p2=ahflw)
    CALL get_channel_object(status,TRIM(vdicha),'ahflw', p2=ahflw)
    CALL channel_halt(substr, status)

!!$ CALL get_channel_object(status,'ECHAM5','ahfsw', p2=ahfsw)
    CALL get_channel_object(status,TRIM(vdicha),'ahfsw', p2=ahfsw)
    CALL channel_halt(substr, status)

!!$ CALL get_channel_object(status,'ECHAM5','ahfli', p2=ahfli)
    CALL get_channel_object(status,TRIM(vdicha),'ahfli', p2=ahfli)
    CALL channel_halt(substr, status)

!!$ CALL get_channel_object(status,'ECHAM5','ahfsi', p2=ahfsi)
    CALL get_channel_object(status,TRIM(vdicha),'ahfsi', p2=ahfsi)
    CALL channel_halt(substr, status)

    ! op_pj_20130408+
    CALL get_channel_info(status, 'rad01')
    IF (status /=0) THEN
       CALL warning_bi( &
            'channel rad01 not available, trying old rad4all ...' &
            , substr)
       CALL get_channel_info(status, 'rad4all')
       IF (status /= 0) THEN
          CALL error_bi(' ... rad4all also not available!', substr)
       ELSE
          radcha = 'rad4all'
       ENDIF
    ELSE
       radcha = 'rad01'
    ENDIF
    ! op_pj_20130408-

    CALL get_channel_object(status, TRIM(radcha), 'trflw', p2=trflw)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, TRIM(radcha), 'trfli', p2=trfli)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, TRIM(radcha), 'soflw', p2=soflw)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, TRIM(radcha), 'sofli', p2=sofli)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','ahfice', p2=ahfice)
    CALL channel_halt(substr, status)

!!$ CALL get_channel_object(status,'ECHAM5','evapi', p2=evapi)
    CALL get_channel_object(status,TRIM(vdicha),'evapi', p2=evapi)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','cvsi', p2=cvsi)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','icecov', p2=icecov)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,'ECHAM5','qres', p2=qres)
    CALL channel_halt(substr, status)

    !
    ! initialise fluxres
    !
    fluxres = 0._dp
    !
    ! READ INPUT DATA FOR FLUX CORRECTION
    !
    IF (p_parallel_io) THEN
       CALL info_bi('This is a MLOCEAN run.', substr)
       !
       ALLOCATE (fluxcorr(dcl%nlon,dcl%nlat,0:13), stat=status)
       IF (status /= 0) CALL error_bi ('Could not allocate fluxcorr.', modstr//' '//substr)
       ! 
       CALL mlocean_readflux(status, fn_mlo_fluxcorr, fluxcorr)
       IF (status /= 0) CALL error_bi ('Could not read flux corr.', modstr//' '//substr)
       !
    END IF
    !
    !  Allocate memory for aflux per PE
    !
    IF (.NOT. ALLOCATED(aflux)) THEN
       ALLOCATE (aflux(dcl%nproma, dcl%ngpblks,0:13), stat = status)
       IF (status /= 0) CALL error_bi ('Could not allocate aflux.', modstr//' '//substr)
    END IF
    !
    NULLIFY (gl_aflux)
    DO i = 0, 13
       IF (p_pe == p_io) gl_aflux => fluxcorr(:,:,i:i)
       CALL scatter_gp (gl_aflux, aflux(:,:,i:i), dcg)
    END DO
    IF (p_parallel_io) THEN
       status = 0
       DEALLOCATE (fluxcorr, stat=status)
       IF (status /= 0) CALL error_bi ('Could not deallocate fluxcorr.', modstr//' '//substr)
    END IF

    ! op_sb_20120928+
    ! if LGMC is running, MLOCEAN must be called in global_end instead of physc.
    CALL get_channel_info(status, 'lgmc')
    L_LGMC = (status == 0)
    ! op_sb_20120928-

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output
    RETURN
  END SUBROUTINE mlocean_init_coupling
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE mlocean_global_start

    ! ------------------------------------------------------------------
    ! This subroutine is called at the beginning of the time loop.
    ! Here, all vectors of the grid-point-fields are accessible.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mlocean_global_start'
    INTEGER                     :: status

  END SUBROUTINE mlocean_global_start
  ! ====================================================================

! op_pj_20110224+
  ! ====================================================================
  SUBROUTINE mlocean_mixlo

    USE messy_main_grid_def_mem_bi,   ONLY : kproma, jrow ! op_pj_20121205

    IMPLICIT NONE

    IF (L_LGMC) RETURN                            ! op_sb_20121205

    CALL mlocean_mlflx(jrow, kproma)              ! op_pj_20121205
    CALL mlocean_mlayerocean(jrow, kproma)        ! op_pj_20121205

  END SUBROUTINE mlocean_mixlo
  ! ====================================================================
! op_pj_20110224-

  ! op_sb_20120921+
  ! ====================================================================
  SUBROUTINE mlocean_global_end

    USE messy_main_grid_def_mem_bi, ONLY : ngpblks  &
                                          ,npromz   &
                                          ,nproma

    IMPLICIT NONE

    INTEGER   :: row
    INTEGER   :: zproma

    ! in case of LGMC (LG moist convection) mlocean
    ! is called in mlocean_global_end instead of physc
    IF (.NOT. L_LGMC) RETURN     

     do row = 1, ngpblks
      if ( row == ngpblks) then
       zproma = npromz
      else
       zproma = nproma
      endif

      CALL mlocean_mlflx(row,zproma)
      CALL mlocean_mlayerocean(row,zproma)
     enddo   

  END SUBROUTINE mlocean_global_end
  ! ====================================================================
  ! op_sb_20120921-

  ! ====================================================================
  SUBROUTINE mlocean_mlflx(jrow, kproma)    ! op_pj_20110224 PRIVATE

    ! The flux correction of the mixed layer ocean is prepared for
    ! the current time step.
    ! Two possibilities to provide a flux correction are possible:
    ! nflxcorr = 1:
    !    A precalculated flux correction can be used. It is only
    !    interpolated in time for the current time step.
    ! nflxcorr = 2:
    !    The net ocean surface energy budget is given with the flux
    !    correction input file. The climatological SST are used to
    !    calculate the final flux correction.
    ! nflxcorr = 0:
    !    Flux correction is set to zero.
    !
    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_data_bi,    ONLY: wgt1, wgt2, nmw1, nmw2 &
                                   , sst, aice

    ! I/O
    INTEGER, INTENT(IN) :: jrow, kproma ! op_pj_20121205

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mlocean_mlflx'

    ! 
    ! Update the flux correction: amlcorr
    !
    SELECT CASE (nflxcorr)
    CASE(1)  ! use prescribed flux correction from input file
       !
       CALL mlocean_mlfluxp (kproma, wgt1, wgt2              &
                            , slf(:,jrow), alake(:,jrow)     &
                            , aflux(:,jrow,nmw1), aflux(:,jrow,nmw2) &
                            , amlcorr(:,jrow))
    CASE(2)  ! use climatological SST and the net ocean surface energy budget,
       !       prescribed by the input file.
       !
       CALL mlocean_mlflux (kproma, nmw1, nmw2, wgt1, wgt2  &
                           , mldmix, slf(:,jrow), amlcorr(:,jrow), alake(:,jrow)   &
                           , sst(:,jrow,nmw1), sst(:,jrow,nmw2)   &
                           , aice(:,jrow,nmw1), aice(:,jrow,nmw2)  &
                           , aflux(:,jrow,nmw1), aflux(:,jrow,nmw2))
       !
    CASE(0)  ! NO FLUX CORRECTION FOR nflxcorr == 0
       !
       amlcorr(:,jrow) = 0 
    END SELECT
    ! 
    RETURN
  END SUBROUTINE mlocean_mlflx
  ! ====================================================================
  SUBROUTINE mlocean_mlayerocean(jrow, kproma)     ! op_pj_20110224 PRIVATE
    !
    ! To change between diffrent realisations of mixed layer oceans,
    ! the integer switch mloswitch (set by the namelist CTRL can be used.
    ! So far only the mixed layer ocean of ECHAM5 is ready to use.
    !
    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_grid_def_bi, ONLY:  philat_2d
    USE messy_main_data_bi,     ONLY:  ahfres &   ! melting of ice in W m-2
! op_pj_20160617+
!!$                                   , amlcorac &
!!$                                   , amlheatac &
! op_pj_20160617-
                                    , ahfcon !&   ! conductive heat flux through ice in W m-2
!!$                                   , ahfres     ! melting of ice in W m-2

    USE messy_main_timer,         ONLY: delta_time 

    ! I/O
    INTEGER, INTENT(IN) :: jrow, kproma ! op_pj_20121205

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mlocean_mlayerocean'
    INTEGER                     :: status
    INTEGER                     :: jl
    LOGICAL, DIMENSION(kproma)  :: lonorth ! .TRUE. for northern latitude
    
    DO jl = 1, kproma
       lonorth(jl) = philat_2d(jl,jrow) > 0.0_dp ! true in northern hemisphere
    END DO

    !  Net ocean surface energy budget:
    !    ->  (latent heat + sensible heat + long wave heating + short wave heating)
    !
    fsnet(:,jrow) = ahflw(:,jrow) + ahfsw(:,jrow) + trflw(:,jrow) + soflw(:,jrow)
    !
    SELECT CASE(mloswitch)
    CASE(1)
       ! 
       ! Calculate surface temperature for open sea gridpoints, ice thickness,..
       !
       CALL mlocean_mlo (kproma, lonorth,  delta_time      &
                       , slm(:,jrow),      alake(:,jrow),    fsnet(:,jrow)  &
                       , evapi(:,jrow),    cvsi(:,jrow),     icecov(:,jrow) &
                       , amlcorr(:,jrow),  seaice(:,jrow),   siced(:,jrow)  &
                       , tsi(:,jrow),      tsw(:,jrow)    &
                       , fluxres(:,jrow),  ahfres(:,jrow),   sni(:,jrow)    &
! op_pj_20160617+
!!$                    , amlcorac(:,jrow), amlheatac(:,jrow))
                       ,                   amlheat(:,jrow))
! op_pj_20160617-
       !
       ! Prognostic calculation of sea-ice temperature
       !
       CALL mlocean_mlicetemp (kproma, delta_time   &
                        , siced(:,jrow),    alake(:,jrow),  slf(:,jrow)     &
                        , trfli(:,jrow),    sofli(:,jrow)    &
                        , ahfsi(:,jrow),    ahfli(:,jrow),  icecov(:,jrow)  &
                        , tsi(:,jrow),      sni(:,jrow),    fluxres(:,jrow) &
                        , ahfres(:,jrow),   qres(:,jrow),   ahfice(:,jrow)  &
                        , ahfcon(:,jrow))
    CASE(2)

    END SELECT
    RETURN
  END SUBROUTINE mlocean_mlayerocean
  ! ====================================================================
  !
  ! ====================================================================
  SUBROUTINE mlocean_free_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to deallocate the memory, which has
    ! been "manually" allocated in mlocean_init_memory.
    ! Note: channel object memory must not be deallocated! This is
    !       performed centrally.
    ! ------------------------------------------------------------------

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mlocean_free_memory'
    INTEGER                     :: status

  END SUBROUTINE mlocean_free_memory
  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE mlocean_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    
    INTRINSIC :: TRIM
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ l_diag_out, fn_mlo_fluxcorr

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mlocean_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE
    WRITE (*,NML=CPL)
    WRITE (*,*) TRIM(substr)//':l_diag_out:     ',l_diag_out
    WRITE (*,*) TRIM(substr)//':fn_mlo_fluxcorr:',TRIM(fn_mlo_fluxcorr)
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE mlocean_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_mlocean_e5
! **********************************************************************
