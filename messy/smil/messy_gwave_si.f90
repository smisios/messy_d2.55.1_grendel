! ***********************************************************************
MODULE messy_gwave_si
! ***********************************************************************

  ! MESSy-SMIL FOR SUBMODEL GWAVE
  !
  ! Author: Andreas Baumgaertner, MPICH, 2010

#if defined(ECHAM5) || defined(BLANK) || defined(CESM1) || defined(MESSYDWARF)

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_constants_mem, ONLY: DP

  USE messy_gwave
#ifdef MESSYTENDENCY
   USE messy_main_tendency_bi,   ONLY: mtend_get_handle,               &
                                       mtend_register,                 &
                                       mtend_get_start_l,              &
                                       mtend_add_l,                    &
                                       mtend_set_sqcst_scal,           &
                                       mtend_id_t,                     &
                                       mtend_id_u,                     &
                                       mtend_id_v

#endif

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: l_gwdrag_u = .TRUE.
  LOGICAL :: l_gwdrag_v = .TRUE.
  LOGICAL :: l_gweddy   = .FALSE.
  LOGICAL :: l_gwheat   = .FALSE.

  REAL(dp), POINTER, DIMENSION(:,:,:) :: dummy

#ifdef MESSYTENDENCY
  INTEGER,SAVE                 :: my_handle
#endif

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: gwave_initialize
#if !defined(MESSYDWARF)
#if defined(ECHAM5) || defined(BLANK) || defined(CESM1)
  PUBLIC :: gwave_init_memory
  PUBLIC :: gwave_physc
#endif
#ifdef ECHAM5
  PUBLIC :: gwave_global_end
#endif
#if defined(ECHAM5) || defined(BLANK) || defined(CESM1)
  PUBLIC :: gwave_free_memory
#endif
#endif
  !PRIVATE :: gwave_read_nml_cpl

CONTAINS

  ! --------------------------------------------------------------------------
  SUBROUTINE gwave_initialize

    ! MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi, warning_bi
    USE messy_main_tools,      ONLY: find_next_free_unit
#if !defined(MESSYDWARF) && !defined(BLANK) || defined(MBM_GWAVE)
    USE messy_main_grid_def_mem_bi, ONLY: lmidatm  & ! ECHAM middle atmosphere switch
         , vct    & ! vertical coefficients table
         , nvclev & ! number of levels with vertical coefficients
         , nn     & ! max meridional wave number for m=0.
         , apzero & ! surface pressure
         , ngl
    USE messy_main_grid_def_bi, ONLY: philat
#endif
    USE messy_gwave_hines,     ONLY: hines_read_nml_ctrl                      &
                                   , lextro, lfront, lozpr, iheatcal, rmscon  &
                                   , kstar, m_min, rms_front, front_thres, pcrit &
                                   , pcons, emiss_lev, emiss_press            &
                                   , lrmscon_lat, rmscon_lo, rmscon_hi        &
                                   , lat_rmscon_lo, lat_rmscon_hi, rmscon_lat &
                                   , icutoff, alt_cutoff

    USE messy_gwave_mk,        ONLY: mk_read_nml_ctrl, gwave_int_fac, p0_gw_mk=>p0_gw
    USE messy_gwave_ym,        ONLY: ym_read_nml_ctrl, p0_gw_ym=>p0_gw


    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gwave_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: jl,jstart,jlast ! sm_gr_20160810

#if defined(MESSYDWARF)
    CALL error_bi('GWAVE needs to be adapted to MESSYDWARF first', substr)
#else
    IF (.NOT.lmidatm) THEN
       CALL warning_bi('No middle atmosphere! - GWAVE has no effect',substr)
       RETURN
    END IF

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL gwave_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('error in gwave_read_nml_ctrl',substr)
    END IF

   CALL p_bcast(gwparam,p_io)
   CALL p_bcast(tf_updates,p_io)

   SELECT CASE (gwparam)
       CASE (1)

          ! set latitude dependent GW source
          !
          ! aim: In particular in coupled (atmosphere/ocean) model runs wave
          ! forcing is sometimes insufficient to produce a QBO.
          ! A general increase of the GW
          ! source (rmscon) may have adverse effects on the stratospheric
          ! mid-latitude circulation. This option allows to increase GW
          ! forcing only in low latitudes. A value of rmscon_lo is used for
          ! latitudes below +/-lat_rmscon_lo, a value of rmscon_hi
          ! above +/- lat_rmscon_hi. Between these latitude boundaries
          ! an interpolation is performed.
          ! Here, default values are given that may be overwritten by
          ! setting parameters in the namelist gwsctl.
          ! This option will only be used if lrmscon_lat=.true.
          ! (default is .false.)

          IF (nn==31)  THEN
             ! note: values have not been sufficiently tuned for T31
             lat_rmscon_lo =  5._dp
             lat_rmscon_hi = 10._dp
             rmscon_lo     = 1.1_dp
             rmscon_hi     = 1.0_dp
          ELSE IF (nn==42)  THEN
             ! note: values have not been sufficiently tuned for T31
             lat_rmscon_lo =  5._dp
             lat_rmscon_hi = 10._dp
             rmscon_lo     = 1.1_dp
             rmscon_hi     = 1.0_dp
          ELSE IF (nn==63)  THEN
             ! note: values have not been sufficiently tuned for T31
             lat_rmscon_lo =  5._dp
             lat_rmscon_hi = 10._dp
             rmscon_lo     = 1.2_dp
             rmscon_hi     = 1.0_dp
          ELSE IF (nn==127)  THEN
             ! note: values have been tuned for use in CMIP5 T127L95 / TP04
             lat_rmscon_lo =  5._dp
             lat_rmscon_hi = 10._dp
             rmscon_lo     = 1.05_dp
             !    IF(lcouple) rmscon_lo     = 1.1_dp
             rmscon_hi     = 1.0_dp
          ELSE
             ! note: these are just default values that have not been tested
             lat_rmscon_lo =  5._dp
             lat_rmscon_hi = 10._dp
             rmscon_lo     = 1.1_dp
             rmscon_hi     = 1.0_dp
          END IF

          IF (p_parallel_io) THEN
             iou = find_next_free_unit(100,200)
             CALL hines_read_nml_ctrl(status, iou, modstr, nn)
             IF (status /= 0) CALL error_bi('hines_read_nml_ctrl',substr)
          ENDIF

          CALL p_bcast (lextro, p_io)
          CALL p_bcast (lfront, p_io)
          CALL p_bcast (lozpr, p_io)
          CALL p_bcast (iheatcal, p_io)
          CALL p_bcast (rmscon, p_io)
          CALL p_bcast (emiss_press, p_io)
          CALL p_bcast (kstar, p_io)
          CALL p_bcast (m_min, p_io)
          CALL p_bcast (rms_front, p_io)
          CALL p_bcast (front_thres, p_io)
          CALL p_bcast (pcrit, p_io)
          CALL p_bcast (pcons, p_io)
          CALL p_bcast (lrmscon_lat,p_io )
          CALL p_bcast (lat_rmscon_lo,p_io )
          CALL p_bcast (lat_rmscon_hi,p_io )
          CALL p_bcast (rmscon_lo,p_io )
          CALL p_bcast (rmscon_hi,p_io )

          CALL p_bcast (icutoff,p_io)
          CALL p_bcast (alt_cutoff,p_io)

          IF ( lrmscon_lat ) THEN
             IF (.NOT. ALLOCATED(rmscon_lat)) ALLOCATE(rmscon_lat(ngl))
             DO jl=1,ngl/2
                IF ( philat(jl) >= lat_rmscon_hi ) THEN
                   rmscon_lat(jl) = rmscon_hi
                ELSE IF ( philat(jl) < lat_rmscon_hi .AND. &
                     philat(jl) > lat_rmscon_lo ) THEN
                   rmscon_lat(jl) = rmscon_hi + (philat(jl)-lat_rmscon_hi)* &
                        (rmscon_lo-rmscon_hi)/(lat_rmscon_lo-lat_rmscon_hi)
                ELSE
                   rmscon_lat(jl) = rmscon_lo
                END IF
                ! copy the values to the other hemisphere
                rmscon_lat(ngl+1-jl)=rmscon_lat(jl)
             END DO

             IF (p_parallel_io) THEN
                WRITE(*,*)"The lat dependent rmscon array is:"
                DO jl = 1,ngl,10
                   jstart = jl
                   jlast  = jl+9
                   IF (jlast > ngl) jlast = ngl
                   WRITE (*, '(10f6.3)') rmscon_lat(jstart:jlast)
                ENDDO
             END IF
          END IF

          ! Note that the function sum is only used to convert an integer
          ! array of one element into an integer.
          emiss_lev = nvclev                              &
               -sum(minloc( vct(1:nvclev)                 &
               +vct(nvclev+1:2*nvclev)*apzero,            &
               vct(1:nvclev)                              &
               +vct(nvclev+1:2*nvclev)*apzero>=emiss_press))
       CASE (2) ! MK
          IF (p_parallel_io) THEN
             iou = find_next_free_unit(100,200)
             CALL mk_read_nml_ctrl(status, iou, modstr)
             IF (status /= 0) CALL error_bi('mk_read_nml_ctrl',substr)
          ENDIF

          CALL p_bcast(gwave_int_fac,p_io)
          CALL p_bcast(p0_gw_mk, p_io)

          emiss_lev = nvclev                              &
               -sum(minloc( vct(1:nvclev)                 &
               +vct(nvclev+1:2*nvclev)*apzero,            &
               vct(1:nvclev)                              &
               +vct(nvclev+1:2*nvclev)*apzero>=p0_gw_mk))
       CASE (3)
          !
       CASE (4)
          !
       CASE (5) ! YM
          IF (p_parallel_io) THEN
             iou = find_next_free_unit(100,200)
             CALL ym_read_nml_ctrl(status, iou, modstr)
             IF (status /= 0) CALL error_bi('ym_read_nml_ctrl',substr)
          ENDIF
          CALL p_bcast (p0_gw_ym, p_io)
    END SELECT

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL gwave_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(substr, 'error in gwave_read_nml_cpl')
    END IF
    CALL p_bcast(l_gwdrag_u,p_io)
    CALL p_bcast(l_gwdrag_v,p_io)
    CALL p_bcast(l_gweddy,p_io)
    CALL p_bcast(l_gwheat,p_io)
#endif

  END SUBROUTINE gwave_initialize
  ! --------------------------------------------------------------------------

#if !defined(MESSYDWARF)
  ! --------------------------------------------------------------------------
  SUBROUTINE gwave_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_blather_bi,       ONLY: warning_bi
    USE messy_main_grid_def_mem_bi,  ONLY: lmidatm, ngl
    USE messy_main_grid_def_bi,      ONLY: philat
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    USE messy_gwave_mk,              ONLY: mkgwinti, fac, usql, nhar, levgws
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gwave_init_memory'
    INTEGER :: status

#ifdef MESSYTENDENCY
    ! THIS MUST BE CALLED BEFORE THE 'RETURN' BELOW, BECAUSE
    ! OTHERWISE, tendency.nml MUST BE ALTERED FOR MA and non-MA SETUPS ...
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register(my_handle, mtend_id_t)
    CALL mtend_register(my_handle, mtend_id_u)
    CALL mtend_register(my_handle, mtend_id_v)
#endif

    IF (.NOT.lmidatm) THEN
       CALL warning_bi('No middle atmosphere! - GWAVE has no effect',substr)
       RETURN
    END IF

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'gwdrag_u', &
         p3=gwdrag_u)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwdrag_u', 'long_name' &
         , c='zonal component of gravity wave drag')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwdrag_u', 'units', c='m/s^2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'dummy', &
         p3=dummy)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dummy', 'long_name' &
         , c='dummy')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dummy', 'units', c='')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'gwdrag_v', &
         p3=gwdrag_v)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwdrag_v', 'long_name' &
         , c='meridional component of gravity wave drag')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwdrag_v', 'units', c='m/s^2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'gwheat', &
         p3=gwheat)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwheat', 'long_name' &
         , c='temperature change due to gravity wave drag')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gwheat', 'units', c='K/s')
    CALL channel_halt(substr, status)

    ! Special channel objects depending on chosen parametrisation
    SELECT CASE (gwparam)
    CASE (1) ! Hines
       CALL new_channel_object(status, modstr,  'gwflux_u', &
            p3=gwflux_u)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwflux_u', 'long_name' &
            , c='zonal component of vertical momentum flux')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwflux_u', 'units', c='Pa')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr,  'gwflux_v', &
            p3=gwflux_v)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwflux_v', 'long_name' &
            , c='meridional component of vertical momentum flux')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwflux_v', 'units', c='Pa')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr,  'gweddy', &
            p3=gweddy)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gweddy', 'long_name' &
            , c='eddy diffusion coefficient due to gravity waves')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gweddy', 'units', c='m2s-1')
       CALL channel_halt(substr, status)

    CASE (2) ! MK
       CALL new_channel_object(status, modstr,  'gweddy', &
            p3=gweddy)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gweddy', 'long_name' &
            , c='eddy diffusion coefficient due to gravity waves')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gweddy', 'units', c='m2s-1')
       CALL channel_halt(substr, status)

    CASE (3)
       !
    CASE (4)
       !
    CASE (5)
       !
       CALL new_channel_object(status, modstr,  'gweddy', &
            p3=gweddy)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gweddy', 'long_name' &
            , c='eddy diffusion coefficient due to gravity waves')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gweddy', 'units', c='m2s-1')
       CALL channel_halt(substr, status)
       !
    END SELECT

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)


    SELECT CASE (gwparam)
    CASE (1) ! Hines
       !
    CASE (2) ! MK
       ALLOCATE(fac(ngl))
       fac(:)=1.
       ALLOCATE(usql(nhar,ngl,2))
       levgws = 1
       CALL mkgwinti(philat, ngl)
    CASE (3)
       !
    CASE (4)
       !
    END SELECT

  END SUBROUTINE gwave_init_memory
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
#if defined(ECHAM5) || defined(BLANK) || defined(CESM1)
  SUBROUTINE gwave_physc

    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_data_bi,       ONLY: temp3d=>tm1, temp_tend=>tte_3d &
                                      , uwind3d=>um1, uwind_tend=>vom_3d  &
                                      , vwind3d=>vm1, vwind_tend=>vol_3d
    USE messy_main_grid_def_mem_bi, ONLY:kproma, nproma, jrow, nlev     &
                                      , lmidatm
    USE messy_main_timer,         ONLY: ztmst=>time_step_len

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER           :: substr='gwave_physc'
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: temp, uwind, vwind
    INTEGER :: jk
    INTEGER :: status ! error status
    REAL(dp) :: zgwdiffco(nproma,nlev)    ! ka_sv_20170406

#ifdef MESSYTENDENCY
    REAL(kind=dp),dimension(nproma,nlev)      :: gwdrag_u_lo
    REAL(kind=dp),dimension(nproma,nlev)      :: gwdrag_v_lo
    REAL(kind=dp),dimension(nproma,nlev)      :: gwheat_lo
#endif

    IF (.NOT. lmidatm) RETURN

#ifdef MESSYTENDENCY
    gwdrag_u_lo(:,:) = 0.0_dp
    gwdrag_v_lo(:,:) = 0.0_dp
    gwheat_lo(:,:)   = 0.0_dp
#endif

    ALLOCATE(temp(nproma,SIZE(temp3d,2))) ! nproma because of hines SMCL
    ALLOCATE(uwind(nproma,SIZE(uwind3d,2)))
    ALLOCATE(vwind(nproma,SIZE(vwind3d,2)))

    ! UPDATE TEMPERATURE AND WINDS
    IF (tf_updates) THEN
#ifndef MESSYTENDENCY
       temp(1:kproma,:)  = temp3d(1:kproma,:,jrow)  + temp_tend(1:kproma,:,jrow)*ztmst
       uwind(1:kproma,:) = uwind3d(1:kproma,:,jrow) + uwind_tend(1:kproma,:,jrow)*ztmst
       vwind(1:kproma,:) = vwind3d(1:kproma,:,jrow) + vwind_tend(1:kproma,:,jrow)*ztmst
#else
       call mtend_get_start_l (mtend_id_t, v0 = temp)
       call mtend_get_start_l (mtend_id_u, v0 = uwind)
       call mtend_get_start_l (mtend_id_v, v0 = vwind)
#endif
    ELSE
       temp(1:kproma,:) = temp3d(1:kproma,:,jrow)
       uwind(1:kproma,:) = uwind3d(1:kproma,:,jrow)
       vwind(1:kproma,:) = vwind3d(1:kproma,:,jrow)
    ENDIF

    SELECT CASE (gwparam)
    CASE(1) ! Hines
       CALL gwave_hines_physc(temp(1:kproma,:)  &
                             ,uwind(1:kproma,:) &
                             ,vwind(1:kproma,:),status &
                             ,zgwdiffco) ! ka_sv_20170406
    CASE(2) ! MK
       CALL gwave_mk_physc(temp(1:kproma,:) &
                          ,uwind(1:kproma,:)&
                          ,vwind(1:kproma,:),status)
    CASE(3) ! HLM
       CALL gwave_hlm_physc(temp(1:kproma,:) &
                           ,uwind(1:kproma,:)&
                           ,vwind(1:kproma,:),status)
    CASE(4) ! AD
       CALL gwave_ad_physc(temp(1:kproma,:) &
                          ,uwind(1:kproma,:)&
                          ,vwind(1:kproma,:),status)
    CASE(5) ! YM
       CALL gwave_ym_physc(temp(1:kproma,:) &
                          ,uwind(1:kproma,:)&
                          ,vwind(1:kproma,:),status)
    END SELECT

    IF (status /= 0) CALL error_bi('gwave_physc',substr)

    IF (l_gwdrag_u) THEN
#ifndef MESSYTENDENCY
       DO jk=1, nlev
          uwind_tend(1:kproma,jk,jrow) = &
               uwind_tend(1:kproma,jk,jrow) + gwdrag_u(1:kproma,jk,jrow)
       END DO
#else
       gwdrag_u_lo(1:kproma,1:nlev) = gwdrag_u(1:kproma,1:nlev,jrow)
       call mtend_add_l (my_handle, mtend_id_u, px = gwdrag_u_lo)
#endif
    ENDIF

    IF (l_gwdrag_v) THEN
#ifndef MESSYTENDENCY
       DO jk=1, nlev
          vwind_tend(1:kproma,jk,jrow) = &
               vwind_tend(1:kproma,jk,jrow) + gwdrag_v(1:kproma,jk,jrow)
       END DO
#else
       gwdrag_v_lo(1:kproma,1:nlev) = gwdrag_v(1:kproma,1:nlev,jrow)
       call mtend_add_l (my_handle, mtend_id_v, px = gwdrag_v_lo)
#endif
    ENDIF

    IF (l_gwheat) THEN
#ifndef MESSYTENDENCY
       DO jk=1, nlev
          temp_tend(1:kproma,jk,jrow) = &
               temp_tend(1:kproma,jk,jrow) + gwheat(1:kproma,jk,jrow)
       END DO
#else
       gwheat_lo(1:kproma,1:nlev) = gwheat(1:kproma,1:nlev,jrow)
       call mtend_add_l (my_handle, mtend_id_t, px = gwheat_lo)
#endif
    ENDIF

    IF (l_gweddy) THEN
       SELECT CASE (gwparam)
       CASE (1)
          gweddy(1:kproma,:,jrow) = zgwdiffco(1:kproma,:)
       CASE (5)
          gweddy(1:kproma,:,jrow) = 0._dp
       END SELECT
    ELSE
       gweddy(1:kproma,:,jrow) = 0._dp
    ENDIF

    DEALLOCATE(temp,uwind,vwind)

  END SUBROUTINE gwave_physc
#endif
  ! --------------------------------------------------------------------------

#ifdef ECHAM5
  ! ------------------------------------------------------------------------
  SUBROUTINE gwave_global_end

  END SUBROUTINE gwave_global_end
  ! ------------------------------------------------------------------------
#endif

  ! --------------------------------------------------------------------------
  SUBROUTINE gwave_free_memory

    USE messy_main_grid_def_mem_bi, ONLY: lmidatm
    USE messy_main_blather_bi,    ONLY: warning_bi
    USE messy_gwave_mk,           ONLY: usql
    USE messy_gwave_hines,        ONLY: rmscon_lat

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER             :: substr='gwave_free_memory'

    IF (.NOT.lmidatm) THEN
       CALL warning_bi('No middle atmosphere! - GWAVE has no effect',substr)
       RETURN
    END IF

    SELECT CASE (gwparam)
    CASE (1)
       !
       IF (ALLOCATED(rmscon_lat)) DEALLOCATE(rmscon_lat)
       !
    CASE (2)
       DEALLOCATE(usql)
    CASE (3)
       !
    CASE (4)
       !
    END SELECT

  END SUBROUTINE gwave_free_memory
 ! --------------------------------------------------------------------------


  ! --------------------------------------------------------------------------
  SUBROUTINE gwave_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (INTERFACE, PRIVATE)
    !

    ! MESSy
    USE messy_main_tools,     ONLY: read_nml_open, read_nml_check &
         , read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'gwave_read_nml_cpl'

    NAMELIST /CPL/  l_gwdrag_u, &
                    l_gwdrag_v, &
                    l_gweddy,   &
                    l_gwheat

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    ! ...


    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE gwave_read_nml_cpl
  ! --------------------------------------------------------------------------

  ! 1: Hines
  ! ========================================================================
  SUBROUTINE gwave_hines_physc(temp, uwind, vwind, status, zgwdiffco)

#if defined(ECHAM5) || defined(CESM1)
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev &
                                        , nvclev, vct
    USE messy_main_grid_def_bi,     ONLY: ilat, gl_twomu, gl_sqcst
    USE messy_main_data_bi,       ONLY: aphm1, apm1, aprflux                  &
#if defined(ECHAM5)
                                      , dalpslm1, dalpsmm1, vom1, dm1, alpsm1 &
                                      , dtlm1, dtmm1, dudlm1, dvdlm1          &
#endif
                                      , aprflux
#endif
    USE messy_gwave_hines

    IMPLICIT NONE

    REAL(dp) ,INTENT(IN) :: temp(:,:)   ! temperature
    REAL(dp) ,INTENT(IN) :: uwind(:,:)  ! meridional wind
    REAL(dp) ,INTENT(IN) :: vwind(:,:)  ! zonal wind

    INTEGER, INTENT(OUT) :: status ! error status

    ! adding output for the turbulence coefficient, HAMMONIA way
    REAL(dp) ,INTENT(OUT) ::  zgwdiffco(:,:)

#if defined(ECHAM5) || defined(CESM1)
    ! Local
    CHARACTER(LEN=*), PARAMETER :: substr='hines_spectrum'


    CALL gwspectrum(jrow, kproma, kproma, nlev          &
            ,aphm1(1:kproma,:),      apm1(1:kproma,:)   &
            ,temp(:,:),  uwind(:,:),  vwind(:,:)        &
            ,aprflux(1:kproma,jrow)                     &
            ,zgwdiffco                                  &
            ,gwdrag_u, gwdrag_v, gwheat, gwflux_u, gwflux_v   &
            ,gl_twomu, gl_sqcst, nvclev, vct            &
#if defined(ECHAM5)
            ,dalpslm1, dalpsmm1, vom1, dm1, alpsm1      &
            ,dtlm1, dtmm1, dudlm1, dvdlm1               &
#else
            ! dummies
            ,aphm1, aphm1, dummy, dummy, aphm1          &
            ,dummy, dummy, dummy,  dummy                &
#endif
            ,ilat, status )
#endif
  END SUBROUTINE gwave_hines_physc
  ! ========================================================================

  ! 2: MK
  ! ****************************************************************
  SUBROUTINE gwave_mk_physc(temp, uwind, vwind, status)
#ifdef ECHAM5
    USE messy_main_constants_mem, ONLY: GRAV0=>g
    USE messy_main_data_bi,       ONLY: press=>press_3d, pressi=>pressi_3d
    USE messy_main_grid_def_mem_bi, ONLY: nlev,kproma,nproma, jrow
    USE messy_main_grid_def_bi,     ONLY: ilat, altitude_gnd
#endif

    USE messy_main_tools,         ONLY: full2half
    USE messy_gwave_mk

    IMPLICIT NONE

    REAL(dp) ,INTENT(IN) :: temp(:,:)   ! temperature
    REAL(dp) ,INTENT(IN) :: uwind(:,:)  ! meridional wind
    REAL(dp) ,INTENT(IN) :: vwind(:,:)  ! zonal wind

    INTEGER, INTENT(OUT) :: status

#ifdef ECHAM5
    ! Or use    pres(:)=ceta(:)*apzero ?

    ! Local
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: u,v,h,t,cp,molvisc &
                                         , dpf,bvfr,pf,ph,th &
                                         , drag_u,drag_v,eddy,gwtemp  &
                                         , su,sv,sh,st,sth,scp,smolvisc&
                                         , sdrag_u,sdrag_v,seddy,sgwtemp &
                                         , placemat, pf_p
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: tempi

    REAL(dp) :: dt
    INTEGER :: jp, jk, plev

    !------------------------
    real(dp) :: lev_step =.1 !25
    integer ::  k, zero_winds

    status = 1 ! ERROR

    plev=63

    ALLOCATE(u(nlev),v(nlev),h(nlev),t(nlev),cp(nlev),molvisc(nlev),dpf(plev),bvfr(plev), pf(plev),ph(plev),th(nlev) &
                 , drag_u(nlev),drag_v(nlev),eddy(nlev),gwtemp(nlev))
    ALLOCATE(tempi(nproma,nlev))

    ALLOCATE(placemat(nlev), pf_p(plev))
    ALLOCATE(su(plev),sv(plev),sh(plev),st(plev),sth(plev),scp(plev),smolvisc(plev) &
         ,sdrag_u(plev),sdrag_v(plev),seddy(plev),sgwtemp(plev))

    ! cp = specific heat capacity
    cp(:) = 1004.64    ! JK-1kg-1
    ! molvisc = molecular viscosity
    molvisc(:) = .5e-4 ! from cmat run

    ! For spline interpolation onto new grid
    pf(1) = p0_gw
    ph(1) = 0.00
    DO k=2,plev
       pf(k) = pf(k-1)*EXP(-1.*lev_step)
       ph(k) = pf(k)*EXP(0.5*lev_step)
       IF(k > 1) dpf(k-1) = pf(k) - pf(k-1)
       pf_p(plev+1-k) = pf(k)
    ENDDO
    pf_p(plev) = pf(1)
    zero_winds = 0    !

    CALL full2half(temp,tempi,press(:,:,jrow),pressi(:,:,jrow))
    DO jp = 1, kproma
       placemat(:)=press(jp,:,jrow)
       DO jk = 1, nlev
          u(nlev+1-jk) = uwind(jp,jk)
          v(nlev+1-jk) = vwind(jp,jk)*(-1.)
          t(nlev+1-jk) = temp(jp,jk)
          h(nlev+1-jk) = altitude_gnd(jp,jk,jrow)
       ENDDO
       DO jk = 1, nlev
          th(jk) = tempi(jp,nlev+1-jk)
       ENDDO

     CALL SplineOn(u, su, 1, zero_winds, placemat, pf_p, nlev, plev, status)
     CALL SplineOn(v, sv, 1, zero_winds, placemat, pf_p, nlev, plev, status)
     CALL SplineOn(t, st, 0, zero_winds, placemat, pf_p, nlev, plev, status)
     CALL SplineOn(h, sh, 2, zero_winds, placemat, pf_p, nlev, plev, status)
     CALL SplineOn(cp, scp, 2, zero_winds, placemat, pf_p, nlev, plev, status)
     CALL SplineOn(molvisc, smolvisc, 2, zero_winds, placemat, pf_p, nlev, plev, status)

   ! Calculate Brunt Vaisailla frequency
    DO k = 1, plev

       ! Evaluate Brunt Vaisala Freq
       IF(k.GT.1.AND.k<plev-1) THEN

          dt=(st(k+1)-st(k-1))/((sh(k+1) - sh(k-1)))
          bvfr(k)=(1./st(k))*(dt + (GRAV0/scp(k)))

          IF(bvfr(k) < 0) THEN
             bvfr(k) = 5.e-3
          ELSE
             bvfr(k)=SQRT((bvfr(k)*GRAV0))
          ENDIF

       ELSE
          bvfr(k)=2.e-2
       ENDIF

       IF(k > 1) sth(k) = 0.5*(st(k) + st(k-1))

    ENDDO

       CALL mkgwintr (su, sv, st, sth, pf, ph, dpf,                             &
         bvfr, sdrag_u, sdrag_v, sgwtemp, ilat(jp,jrow), fac, scp, seddy,           &
         smolvisc, plev)

       CALL SplineOn(sdrag_u,  drag_u,  0, zero_winds, pf_p, placemat, plev, nlev, status)
       CALL SplineOn(sdrag_v,  drag_v,  0, zero_winds, pf_p, placemat, plev, nlev, status)
       CALL SplineOn(seddy,    eddy,    0, zero_winds, pf_p, placemat, plev, nlev, status)
       CALL SplineOn(sgwtemp,  gwtemp,  0, zero_winds, pf_p, placemat, plev, nlev, status)

       ! reverse result!
       DO jk = 1, nlev
         gwdrag_u(jp,jk,jrow) = drag_u(nlev+1-jk)
         gwdrag_v(jp,jk,jrow) = drag_v(nlev+1-jk)*(-1.)
         gwheat  (jp,jk,jrow) = gwtemp(plev+1-jk)
         gweddy  (jp,jk,jrow) = eddy(plev+1-jk)
       ENDDO
    ENDDO

    DEALLOCATE(u,v,h,t,cp,molvisc,dpf,bvfr, pf,ph,th &
               , drag_u,drag_v,eddy,gwtemp)
    DEALLOCATE(tempi)

    status = 0 ! NO ERROR

#endif
  END SUBROUTINE gwave_mk_physc

  ! 3: HLM
  ! ========================================================================
   SUBROUTINE gwave_hlm_physc(temp, uwind, vwind, status)
#ifdef ECHAM5
    USE messy_main_constants_mem,   ONLY: GRAV0=>g, M_air, R_gas
    USE messy_main_data_bi,         ONLY: press=>press_3d
    USE messy_main_grid_def_bi,     ONLY: altitude_gnd
    USE messy_main_grid_def_mem_bi, ONLY: nlev,kproma, jrow
#endif
    USE messy_gwave_hlm

    IMPLICIT NONE

    REAL(dp) ,INTENT(IN) :: temp(:,:)   ! temperature
    REAL(dp) ,INTENT(IN) :: uwind(:,:)  ! meridional wind
    REAL(dp) ,INTENT(IN) :: vwind(:,:)  ! zonal wind

    INTEGER, INTENT(OUT) :: status

#ifdef ECHAM5
    ! Or use    pres(:)=ceta(:)*apzero ?

    ! Local
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: u,v,h,t,cp &
                                         , pf, drag_u,drag_v,eddy,gwtemp &
                                         , gw_eddy_u, gw_eddy_v &
                                         , gwh_u, gwh_v &
                                         , scht

    INTEGER :: jp, jk, plev, direction

    status = 1 ! ERROR

    plev = nlev-30

    ALLOCATE(u(plev),v(plev),h(plev),t(plev),cp(plev) &
            , pf(plev), drag_u(plev),drag_v(plev),eddy(plev),gwtemp(plev))
    ALLOCATE(gw_eddy_u(plev), gw_eddy_v(plev) &
            , gwh_u(plev), gwh_v(plev) &
            , scht(plev))

    ! cp = specific heat capacity
    cp(:) = 1004.64    ! JK-1kg-1

    DO jp = 1, kproma
       DO jk = 1, plev
          u(plev+1-jk)  = uwind(jp,jk)
          v(plev+1-jk)  = vwind(jp,jk)*(-1.)
          t(plev+1-jk)  = temp(jp,jk)
          h(plev+1-jk)  = altitude_gnd(jp,jk,jrow)
          pf(plev+1-jk) = press(jp,jk,jrow)
       ENDDO


       ! scale height (m). M_air only constant up to ~mesopause!
       scht(:) = R_gas*1e3_dp*t/(M_air*GRAV0)

       direction = 1
       CALL HLMGwave(u, scht, t, cp, h,   &
            drag_u, gw_eddy_u, pf, direction,  &
            gwh_u, plev)

       direction = 2
       CALL HLMGwave(v, scht, t, cp, h,   &
            drag_v, gw_eddy_v, pf, direction,  &
            gwh_v, plev)

       ! reverse result!
       DO jk = 1, plev
         gwdrag_u(jp,jk,jrow)=drag_u(plev+1-jk)
         gwdrag_v(jp,jk,jrow)=drag_v(plev+1-jk)*(-1.)
         gwheat(jp,jk,jrow)=(gwh_u(plev+1-jk) + gwh_v(plev+1-jk))*cp(jk)
       !    eddy = gw_eddy_u + gw_eddy_v
       ENDDO

    ENDDO

    DEALLOCATE(u,v,h,t,cp, pf, drag_u,drag_v,eddy,gwtemp &
              ,gw_eddy_u, gw_eddy_v, gwh_u, gwh_v, scht)

    status = 0 ! NO ERROR

    RETURN
#endif
  END SUBROUTINE gwave_hlm_physc

  ! 4: AD
  ! ========================================================================
  SUBROUTINE gwave_ad_physc(temp, uwind, vwind, status)
#ifdef ECHAM5
    USE messy_main_constants_mem,   ONLY: GRAV0=>g, M_air, R_gas
    USE messy_main_data_bi,         ONLY: press=>press_3d
    USE messy_main_grid_def_bi,     ONLY: altitude_gnd
    USE messy_main_grid_def_mem_bi, ONLY: nlev,kproma, jrow
#endif
    USE messy_gwave_ad

    IMPLICIT NONE

    REAL(dp) ,INTENT(IN) :: temp(:,:)   ! temperature
    REAL(dp) ,INTENT(IN) :: uwind(:,:)  ! meridional wind
    REAL(dp) ,INTENT(IN) :: vwind(:,:)  ! zonal wind

    INTEGER, INTENT(OUT) :: status

#ifdef ECHAM5
    ! Or use    pres(:)=ceta(:)*apzero ?

    ! Local
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: u,v,h,t,cp &
                                         , pf, drag_u,drag_v,eddy,gwtemp &
                                         , gw_eddy_u, gw_eddy_v &
                                         , gwh_u, gwh_v &
                                         , scht, density, bruntv
    INTEGER :: jp, jk, plev, direction, flag

    status = 1 ! ERROR

    plev = nlev-10

    ALLOCATE(u(plev),v(plev),h(plev),t(plev),cp(plev) &
            , pf(plev), drag_u(plev),drag_v(plev),eddy(plev),gwtemp(plev))
    ALLOCATE(gw_eddy_u(plev), gw_eddy_v(plev) &
            , gwh_u(plev), gwh_v(plev) &
            , scht(plev), density(plev) &
            , bruntv(plev))

    ! cp = specific heat capacity
    cp(:) = 1004.64    ! JK-1kg-1

    DO jp = 1, kproma
       DO jk = 1, plev
          u(plev+1-jk) = uwind(jp,jk)
          v(plev+1-jk) = vwind(jp,jk)*(-1.)
          t(plev+1-jk) = temp(jp,jk)
          h(plev+1-jk) = altitude_gnd(jp,jk,jrow)
          pf(plev+1-jk) = press(jp,jk,jrow)
       ENDDO


       ! scale height (m). M_air only constant up to ~mesopause!
       scht(:) = R_gas*1e3_dp*t/(M_air*GRAV0)
       ! density (kg/m3)
       density(:) = pf(:)*M_air/(t(:)*R_gas)*1e-3

       flag = 0

       direction = 1
       CALL AlexDGwave(flag, u, scht, t, &
            cp, h, drag_u, gw_eddy_u, pf, density, &
            bruntv, direction, gwh_u, plev)

       direction = 2
       CALL AlexDGwave(flag, v, scht, t, &
            cp, h, drag_v, gw_eddy_v, pf, density, &
            bruntv, direction, gwh_v, plev)

       ! reverse result!
       DO jk = 1, plev
         gwdrag_u(jp,jk,jrow)=drag_u(plev+1-jk)
         gwdrag_v(jp,jk,jrow)=drag_v(plev+1-jk)*(-1.)
         gwheat(jp,jk,jrow)=(gwh_u(plev+1-jk) + gwh_v(plev+1-jk))*cp(jk)
       !    eddy = gw_eddy_u + gw_eddy_v
       ENDDO

    ENDDO

    DEALLOCATE(u,v,h,t,cp, pf, drag_u,drag_v,eddy,gwtemp &
              ,gw_eddy_u, gw_eddy_v, gwh_u, gwh_v &
              ,scht, density, bruntv)

    status = 0 ! NO ERROR

    RETURN
#endif
  END SUBROUTINE gwave_ad_physc

  ! 5: YM
  ! ****************************************************************
  SUBROUTINE gwave_ym_physc(temp, uwind, vwind, status)
#ifdef ECHAM5
    USE messy_main_constants_mem,   ONLY: M_air, R_gas
    USE messy_main_data_bi,         ONLY: press=>press_3d
    USE messy_main_grid_def_bi,     ONLY: altitude_gnd
    USE messy_main_grid_def_mem_bi, ONLY: nlev,kproma,nproma, jrow
#endif
    USE messy_gwave_ym

    IMPLICIT NONE

    REAL(dp) ,INTENT(IN) :: temp(:,:)   ! temperature
    REAL(dp) ,INTENT(IN) :: uwind(:,:)  ! meridional wind
    REAL(dp) ,INTENT(IN) :: vwind(:,:)  ! zonal wind

    INTEGER, INTENT(OUT) :: status

#ifdef ECHAM5
    ! Or use    pres(:)=ceta(:)*apzero ?

    ! Local
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: u,v,h,t,cp,rho &
                                         , dpf,bvfr,pf,ph,th &
                                         , drag_u,drag_v,eddy,gwtemp  &
                                         , su,sv,sh,st,sth,scp,srho&
                                         , sdrag_u,sdrag_v,seddy,sgwtemp &
                                         , placemat, pf_p
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: tempi

    INTEGER :: jp, jk, plev

    real(dp) :: lev_step =.25
    integer ::  k, zero_winds
    status = 1 ! ERROR

    plev = 58   ! in paper z0 (launch level) = 16 km

    ALLOCATE(u(nlev),v(nlev),h(nlev),t(nlev),cp(nlev),rho(nlev),dpf(plev),bvfr(plev), pf(plev),ph(plev),th(nlev) &
                 , drag_u(nlev),drag_v(nlev),eddy(nlev),gwtemp(nlev))
    ALLOCATE(tempi(nproma,nlev))

    ALLOCATE(placemat(nlev), pf_p(plev))
    ALLOCATE(su(plev),sv(plev),sh(plev),st(plev),sth(plev),scp(plev),srho(plev) &
         ,sdrag_u(plev),sdrag_v(plev),seddy(plev),sgwtemp(plev))


    ! For spline interpolation onto new grid
    pf(1) = p0_gw
    DO k=2,plev
       pf(k) = pf(k-1)*EXP(-1.*lev_step)
       IF(k > 1) dpf(k-1) = pf(k) - pf(k-1)
       pf_p(plev+1-k) = pf(k)
    ENDDO
    pf_p(plev) = pf(1)
    zero_winds = 0    !

    DO jp = 1, kproma
       placemat(:)=press(jp,:,jrow)
       DO jk = 1, nlev
          u(nlev+1-jk) = uwind(jp,jk)
          v(nlev+1-jk) = vwind(jp,jk)*(-1.)
          t(nlev+1-jk) = temp(jp,jk)
          h(nlev+1-jk) = altitude_gnd(jp,jk,jrow)
          rho(nlev+1-jk) = press(jp,jk,jrow)*M_air/(temp(jp,jk)*R_gas)*1e-6_dp
       ENDDO

       CALL SplineOn(u, su, 1, zero_winds, placemat, pf_p, nlev, plev, status)
       CALL SplineOn(v, sv, 1, zero_winds, placemat, pf_p, nlev, plev, status)
       CALL SplineOn(t, st, 0, zero_winds, placemat, pf_p, nlev, plev, status)
       CALL SplineOn(h, sh, 2, zero_winds, placemat, pf_p, nlev, plev, status)
       CALL SplineOn(rho, srho, 2, zero_winds, placemat, pf_p, nlev, plev, status)

       CALL EYGwave(su*0., srho, su, sv, st, sh, &
       sdrag_u, sdrag_v, pf, su*0., plev)


       CALL SplineOn(sdrag_u,  drag_u,  0, zero_winds, pf_p, placemat, plev, nlev, status)
       CALL SplineOn(sdrag_v,  drag_v,  0, zero_winds, pf_p, placemat, plev, nlev, status)

       ! reverse result!
       DO jk = 1, nlev
         gwdrag_u(jp,jk,jrow) = drag_u(nlev+1-jk)
         gwdrag_v(jp,jk,jrow) = drag_v(nlev+1-jk)*(-1.)
       ENDDO
    ENDDO

    DEALLOCATE(u,v,h,t,cp,rho,dpf,bvfr, pf,ph,th &
               , drag_u,drag_v,eddy,gwtemp)
    DEALLOCATE(tempi)

    status = 0 ! NO ERROR

#endif
  END SUBROUTINE gwave_ym_physc
  ! ****************************************************************

#endif

!#else from #if defined(ECHAM5) || defined(BLANK) || defined(MESSYDWARF)
#else

IMPLICIT NONE

#endif

! ***********************************************************************
END MODULE messy_gwave_si
! ***********************************************************************
