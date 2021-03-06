!**********************************************************************
MODULE messy_clamscirrus_si
!**********************************************************************
!  Submodel interface for clamscirrus 
!**********************************************************************

  USE messy_clamscirrus
  USE messy_main_blather_bi,  ONLY: start_message_bi, end_message_bi
  USE messy_main_timer_event, ONLY: time_event, io_time_event
  USE messy_clams_global, ONLY: species_type_3d

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  !CLaMS
  REAL(DP), DIMENSION(:), POINTER :: LAT        => NULL()
  REAL(DP), DIMENSION(:), POINTER :: LON        => NULL()
  REAL(DP), DIMENSION(:), POINTER :: LEV        => NULL()
  REAL(DP), DIMENSION(:), POINTER :: TEMP       => NULL()
  REAL(DP), DIMENSION(:), POINTER :: PRESS      => NULL()
  REAL(DP), DIMENSION(:), POINTER :: THETA      => NULL()
  REAL(DP), DIMENSION(:), POINTER :: PV         => NULL()
      
  REAL(DP), DIMENSION(:), POINTER :: H2O        => NULL()
  REAL(DP), DIMENSION(:), POINTER :: H2O_100    => NULL()
  REAL(DP), DIMENSION(:), POINTER :: IWC        => NULL()
  REAL(DP), DIMENSION(:), POINTER :: IWC_100    => NULL()
  REAL(DP), DIMENSION(:), POINTER :: CLWC       => NULL()

#ifdef ECHAM5
  !EMAC
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LAT3D_D   => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LON3D_D   => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_ZETA_D    => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_TEMP_D    => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_PRESS_D   => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_H2O_D     => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_H2O_100_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_IWC_D     => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_IWC_100_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_CLWC_D    => NULL()
  ! op_pj_20170110+
  LOGICAL, SAVE :: L_CLAMSCHEME5 = .FALSE. ! is CLAMSCHEME5 active ?
  ! op_pj_20170110-
#endif

  TYPE(species_type_3d), DIMENSION(:), POINTER :: CIRRUS_tte ! ECHAMTRACER Tendency

! op_pj_20160606+
!!$  TYPE(time_event), PUBLIC :: cirrusevent
!!$  TYPE(io_time_event):: io_cirrusevent
  TYPE(time_event), PUBLIC, SAVE :: cirrusevent
  TYPE(io_time_event), SAVE :: io_cirrusevent
! op_pj_20160606-

  PUBLIC :: clamscirrus_initialize
  PUBLIC :: clamscirrus_init_memory
  PUBLIC :: clamscirrus_init_coupling
  PUBLIC :: clamscirrus_global_end
  PUBLIC :: clamscirrus_free_memory

!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE clamscirrus_initialize

    USE messy_main_tools,      ONLY: find_next_free_unit
    USE messy_clamscirrus_global, ONLY: timestep_cirrus
    USE messy_main_timer,      ONLY: delta_time
    USE messy_main_timer_bi,   ONLY: timer_event_init
    USE messy_main_blather_bi, ONLY: error_bi

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscirrus_initialize'
    INTEGER :: status, iou

    ! Read namelist variables:
    iou = find_next_free_unit(100,200)

    ! Read namelist and set default values:
    CALL clamscirrus_read_nml(status, iou)
    IF (status /= 0) CALL error_bi('Error in clamscirrus_read_nml !',substr)

    if (mod(timestep_cirrus*3600,int(delta_time)) /= 0) then
       call error_bi ("Wrong cirrus timestep !!!",substr)
    endif

    ! Define CIRRUS event:
    io_cirrusevent%counter = timestep_cirrus   
    io_cirrusevent%unit = 'hours'
    io_cirrusevent%adjustment = 'exact'
    io_cirrusevent%offset = -delta_time
    CALL timer_event_init (cirrusevent, io_cirrusevent, 'CIRRUS_Event', 'present')

#ifdef ECHAM5
    ALLOCATE(CIRRUS_tte(3))
    CIRRUS_tte(1)%name = 'H2O'
    CIRRUS_tte(2)%name = 'IWC'
    CIRRUS_tte(3)%name = 'CLWC'
#endif

  END SUBROUTINE clamscirrus_initialize
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamscirrus_init_memory

    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_mpi_bi,           ONLY: p_pe

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscirrus_init_memory'
    INTEGER :: status, i

    IF (p_pe==0) WRITE(*,*) substr

    CALL new_channel(status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt(substr, status)


#ifdef ECHAM5
    DO i = 1,3
       CALL new_channel_object(status, modstr, TRIM(CIRRUS_tte(i)%name)//'_cirtte',&
            p3=CIRRUS_tte(i)%values, reprid=GP_3D_MID, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
    END DO
#endif

  END SUBROUTINE clamscirrus_init_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamscirrus_init_coupling

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: error_bi

    USE messy_main_channel,      ONLY: get_channel_object, get_channel_info
    USE messy_clams_global,      ONLY: init_vertcoorname, nspec, specnames
    USE messy_clams_tools_utils, ONLY: str_found
 
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscirrus_init_coupling'
    INTEGER :: status

    CALL start_message_bi(modstr, 'COUPLING TO DRIVER FIELDS', substr)

    ! Couple positions (LAT, LON, LEV) 
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=LEV)
    CALL channel_halt(substr, status)

    ! Couple parameters (TEMP, PRESS, THETA, PV)
    CALL get_channel_object(status, 'clams', 'TEMP', p1=TEMP)
    CALL channel_halt(substr, status)
    IF (TRIM(init_vertcoorname) == 'press') THEN 
       CALL get_channel_object(status, 'clams', 'LEV', p1=PRESS)
       CALL channel_halt(substr, status)
    ELSE
       CALL get_channel_object(status, 'clams', 'PRESS', p1=PRESS)
       CALL channel_halt(substr, status)
    ENDIF
    CALL get_channel_object(status, 'clams', 'THETA', p1=THETA)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'PV', p1=PV)
    CALL channel_halt(substr, status)

    ! Check, if the species, which will be coupled, exist
    if (.not. str_found (nspec, specnames, 'H2O') .or. &
        .not. str_found (nspec, specnames, 'H2O_100') .or. &
        .not. str_found (nspec, specnames, 'IWC') .or. &
        .not. str_found (nspec, specnames, 'IWC_100') .or. &
        .not. str_found (nspec, specnames, 'CLWC') ) then
       call error_bi &
            ('One or more of the species used for CIRRUS missing: '// &
            'H2O, H2O_100, IWC, IWC_100, CLWC !!!', substr)
    endif

    ! Couple species
    CALL get_channel_object(status, 'clams', 'H2O', p1=H2O)
    CALL channel_halt(substr, status)     
    CALL get_channel_object(status, 'clams', 'H2O_100', p1=H2O_100)
    CALL channel_halt(substr, status)     
    CALL get_channel_object(status, 'clams', 'IWC', p1=IWC)
    CALL channel_halt(substr, status)     
    CALL get_channel_object(status, 'clams', 'IWC_100', p1=IWC_100)
    CALL channel_halt(substr, status)     
    CALL get_channel_object(status, 'clams', 'CLWC', p1=CLWC)
    CALL channel_halt(substr, status)

#ifdef ECHAM5
    CALL get_channel_object(status, 'clams', 'E5_LAT3D', p3=E5_LAT3D_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'E5_LON3D', p3=E5_LON3D_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'E5_ZETA', p3=E5_ZETA_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'ECHAM5', 'tm1', p3=E5_TEMP_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'ECHAM5', 'press', p3=E5_PRESS_D)
    CALL channel_halt(substr, status)
    ! op_pj_20170110+
    ! check, if clamscheme5 is active
    CALL get_channel_info(status,'clamscheme5')
    L_CLAMSCHEME5 = (status == 0)
    ! op_pj_20170110-
#endif

    CALL end_message_bi(modstr, 'COUPLING TO DRIVER FIELDS', substr)

  END SUBROUTINE clamscirrus_init_coupling
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamscirrus_global_end

    USE messy_main_timer,         ONLY: lstart, YEAR, MONTH, &
                                        DAY, HOUR, MINUTE, SECOND, &
                                        delta_time, time_step_len
    USE messy_clamscirrus_global, ONLY: use_traj, timestep_cirrus
    USE messy_main_mpi_bi,        ONLY: p_pe
    USE messy_clams_global,       ONLY: mdi, dnparts_max, dnparts, &
                                        timetype, dates30, irdsec, irdday, lcirrusevent,&
                                        fut_year, fut_month, fut_day, fut_sec, &
                                        H2O_index, H2O_100_index, IWC_index, &
                                        IWC_100_index, CLWC_index, init_h2o_emac, &
                                        ldiagout, lcoupled, &
                                        nparams, paramnames
    USE messy_clams_tools_dateconv, ONLY: incdat_interface
    USE messy_clams_tools_interpolreg, ONLY: interpolate_param
    USE messy_clams_tools_utils,       ONLY: str_pos
!!$ USE messy_main_switch,  ONLY: USE_CLAMSCHEME5 ! op_pj_20170110
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi,   ONLY: xtte, xtm1
    USE messy_main_grid_def_mem_bi, ONLY: nlev, ngpblks, nproma, npromz
#endif

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscirrus_global_end'

    INTEGER  :: status 
    TYPE(timetype) :: itime, ipasttime
#ifdef ECHAM5
    INTEGER :: ilat, ilev, icnt, counter, ncnt
    INTEGER :: E5_nlon_d, E5_nlat_d, E5_nlev_d
    INTEGER :: dngridpoints
    REAL(DP), DIMENSION(:), POINTER :: E5_TEMP_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_PRESS_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_LAT_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_LON_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_ZETA_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_H2O_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_H2O_100_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_IWC_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_IWC_100_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_CLWC_1D
#endif


#ifdef ECHAM5
    IF (lstart) THEN
       IF (init_h2o_emac) THEN
          ! Interpolation of H2O on initial air parcel positions
          itime%year    = YEAR
          itime%month   = MONTH
          itime%day     = DAY
          itime%sec     = HOUR*3600+MINUTE*60+SECOND+delta_time
          
          ! Calculate time of previous windfile
          ipasttime%sec = fut_sec
          ipasttime%day = fut_day
          ipasttime%month = fut_month
          ipasttime%year = fut_year
          CALL incdat_interface( ipasttime%sec,ipasttime%day,ipasttime%month,ipasttime%year, &
            -irdsec,-irdday,0,0, dates30)     
          ! Interpolation of H2O
          H2O = mdi
          CALL interpolate_param(LON, LAT, LEV, itime, ipasttime, &
               str_pos(nparams,paramnames,'H2O'), 'H2O', H2O, lcoupled)
          H2O_100 = H2O
       END IF
       ! IWC set to zero
       IWC = mdi
       IWC(1:dnparts) = 0.
       IWC_100 = IWC
       ! CLWC set to zero
       CLWC = mdi
       CLWC(1:dnparts) = 0.
    END IF
#endif

    IF (lcirrusevent) THEN
       IF (p_pe==0 .and. ldiagout) WRITE(*,*) 'Active CIRRUS event'
       
       use_traj = .TRUE. 
       IF (dnparts>0) THEN
          CALL CIRRUS(LAT, LON, THETA, TEMP, PRESS, PV, H2O, H2O_100, IWC, IWC_100, CLWC, &
               timestep_cirrus*3600._DP, dnparts)
       ENDIF

!!$#ifdef ECHAM5
!!$       ALLOCATE (E5_TEMP_1D(dnparts_max))
!!$       ALLOCATE (E5_PRESS_1D(dnparts_max))
!!$       ALLOCATE (E5_LAT_1D(dnparts_max))
!!$       ALLOCATE (E5_LON_1D(dnparts_max))
!!$       ALLOCATE (E5_ZETA_1D(dnparts_max))
!!$       ALLOCATE (E5_H2O_1D(dnparts_max))
!!$       ALLOCATE (E5_H2O_100_1D(dnparts_max))
!!$       ALLOCATE (E5_IWC_1D(dnparts_max))
!!$       ALLOCATE (E5_IWC_100_1D(dnparts_max))
!!$       ALLOCATE (E5_CLWC_1D(dnparts_max))
!!$       E5_TEMP_1D = mdi
!!$       E5_PRESS_1D = mdi
!!$       E5_LAT_1D = mdi
!!$       E5_LON_1D = mdi
!!$       E5_ZETA_1D = mdi
!!$       E5_H2O_1D = mdi
!!$       E5_H2O_100_1D = mdi
!!$       E5_IWC_1D = mdi
!!$       E5_IWC_100_1D = mdi
!!$       E5_CLWC_1D = mdi
!!$
!!$       E5_nlat_d = SIZE(E5_TEMP_D,3)
!!$       E5_nlon_d = SIZE(E5_TEMP_D,1)
!!$       E5_nlev_d = SIZE(E5_TEMP_D,2)
!!$       dngridpoints = E5_nlat_d*E5_nlon_d*E5_nlev_d
!!$
!!$       counter = 0
!!$       DO ilev = 0, nlev-1
!!$          DO ilat = 0, ngpblks-1
!!$             IF (ilat==ngpblks-1) THEN
!!$                ncnt = npromz
!!$             ELSE
!!$                ncnt = nproma
!!$             END IF
!!$             DO icnt = 0, ncnt-1
!!$                E5_TEMP_1D (counter+1) = E5_TEMP_D (icnt+1,ilev+1,ilat+1)
!!$                E5_PRESS_1D(counter+1) = E5_PRESS_D(icnt+1,ilev+1,ilat+1)/100.
!!$                E5_LAT_1D  (counter+1) = E5_LAT3D_D(icnt+1,ilev+1,ilat+1)
!!$                E5_LON_1D  (counter+1) = E5_LON3D_D(icnt+1,ilev+1,ilat+1)
!!$                E5_ZETA_1D (counter+1) = E5_ZETA_D (icnt+1,ilev+1,ilat+1)
!!$                
!!$                E5_H2O_1D  (counter+1) = xtm1(icnt+1,ilev+1,H2O_index,ilat+1) &
!!$                     + xtte(icnt+1,ilev+1,H2O_index,ilat+1) * time_step_len
!!$                IF (E5_H2O_1D(counter+1).LT.0.) E5_H2O_1D(counter+1) = 0.
!!$                
!!$                E5_H2O_100_1D(counter+1) = xtm1(icnt+1,ilev+1,H2O_100_index,ilat+1) &
!!$                     + xtte(icnt+1,ilev+1,H2O_100_index,ilat+1) * time_step_len
!!$                IF (E5_H2O_100_1D(counter+1).LT.0.) E5_H2O_100_1D(counter+1) = 0.
!!$                
!!$                E5_IWC_1D(counter+1)     = xtm1(icnt+1,ilev+1, IWC_index,ilat+1) &
!!$                     + xtte(icnt+1,ilev+1, IWC_index,ilat+1) * time_step_len
!!$                IF (E5_IWC_1D(counter+1).LT.0.) E5_IWC_1D(counter+1) = 0.
!!$                
!!$                E5_IWC_100_1D(counter+1) = xtm1(icnt+1,ilev+1, IWC_100_index,ilat+1) &
!!$                     + xtte(icnt+1,ilev+1, IWC_100_index,ilat+1) * time_step_len
!!$                IF (E5_IWC_100_1D(counter+1).LT.0.) E5_IWC_100_1D(counter+1) = 0.
!!$                
!!$                E5_CLWC_1D(counter+1)    = xtm1(icnt+1,ilev+1, CLWC_index,ilat+1) &
!!$                     + xtte(icnt+1,ilev+1, CLWC_index,ilat+1) * time_step_len
!!$                IF (E5_CLWC_1D(counter+1).LT.0.) E5_CLWC_1D(counter+1) = 0.
!!$                
!!$                counter = counter + 1
!!$             END DO
!!$          END DO
!!$       END DO
!!$       
!!$       IF (dngridpoints > 0) THEN
!!$          use_traj = .TRUE. 
!!$          CALL CIRRUS(E5_LAT_1D, E5_LON_1D, E5_ZETA_1D, E5_TEMP_1D, E5_PRESS_1D, PV, &
!!$               E5_H2O_1D, E5_H2O_100_1D, E5_IWC_1D, E5_IWC_100_1D, E5_CLWC_1D, &
!!$               timestep_cirrus*3600._DP, dngridpoints)
!!$       END IF
!!$      
!!$       counter = 0
!!$       DO ilev = 0, nlev-1
!!$          DO ilat = 0, ngpblks-1
!!$             IF (ilat==ngpblks-1) THEN
!!$                ncnt = npromz
!!$             ELSE
!!$                ncnt = nproma
!!$             END IF
!!$             DO icnt = 0, ncnt-1
!!$                IF (E5_H2O_1D(counter+1) .GE. 0.) THEN
!!$                   CIRRUS_tte(1)%values(icnt+1,ilev+1,ilat+1) = &
!!$                        (((E5_H2O_1D(counter+1)-xtm1(icnt+1,ilev+1,H2O_index,ilat+1)) / time_step_len) &
!!$                        -  xtte(icnt+1,ilev+1,H2O_index,ilat+1)) / (timestep_cirrus*3600./time_step_len)
!!$                ELSE
!!$                   CIRRUS_tte(1)%values(icnt+1,ilev+1,ilat+1) = 0.
!!$                END IF
!!$                IF (E5_IWC_1D(counter+1) .GE. 0.) THEN
!!$                   CIRRUS_tte(2)%values(icnt+1,ilev+1,ilat+1) = &
!!$                        (((E5_IWC_1D(counter+1)-xtm1(icnt+1,ilev+1,IWC_index,ilat+1)) / time_step_len) &
!!$                        -  xtte(icnt+1,ilev+1,IWC_index,ilat+1)) / (timestep_cirrus*3600./time_step_len)
!!$                ELSE
!!$                   CIRRUS_tte(2)%values(icnt+1,ilev+1,ilat+1)= 0.
!!$                END IF
!!$                IF (E5_CLWC_1D(counter+1) .GE. 0.) THEN
!!$                   CIRRUS_tte(3)%values(icnt+1,ilev+1,ilat+1) = &
!!$                        (((E5_CLWC_1D(counter+1)-xtm1(icnt+1,ilev+1,CLWC_index,ilat+1)) / time_step_len) &
!!$                        -  xtte(icnt+1,ilev+1,CLWC_index,ilat+1)) / (timestep_cirrus*3600./time_step_len)
!!$                ELSE
!!$                   CIRRUS_tte(3)%values(icnt+1,ilev+1,ilat+1) = 0.
!!$                END IF
!!$                counter = counter + 1
!!$             END DO
!!$          END DO
!!$       END DO
!!$       
!!$       DEALLOCATE(E5_TEMP_1D)
!!$       DEALLOCATE(E5_PRESS_1D)
!!$       DEALLOCATE(E5_LAT_1D)
!!$       DEALLOCATE(E5_LON_1D)
!!$       DEALLOCATE(E5_ZETA_1D)
!!$       DEALLOCATE(E5_H2O_1D)
!!$       DEALLOCATE(E5_H2O_100_1D)
!!$       DEALLOCATE(E5_IWC_1D)
!!$       DEALLOCATE(E5_IWC_100_1D)
!!$       DEALLOCATE(E5_CLWC_1D)
!!$#endif
    END IF !lcirrusevent

!!$#ifdef ECHAM5
!!$    ! Add tendencies due to CIRRUS (every timestep):
!!$    ! H2O
!!$    xtte(:,:,H2O_index,:) = xtte(:,:,H2O_index,:) + CIRRUS_tte(1)%values(:,:,:) 
!!$    ! IWC
!!$    xtte(:,:,IWC_index,:) = xtte(:,:,IWC_index,:) + CIRRUS_tte(2)%values(:,:,:) 
!!$    ! CLWC
!!$    xtte(:,:,CLWC_index,:) = xtte(:,:,CLWC_index,:) + CIRRUS_tte(3)%values(:,:,:) 
!!$! op_pj_20170110+
!!$! IF (.NOT. USE_CLAMSCHEME5) THEN
!!$    IF (.NOT. L_CLAMSCHEME5) THEN
!!$! op_pj_20170110-
!!$       WHERE (xtm1(:,:,H2O_index,:)+xtte(:,:,H2O_index,:)*time_step_len .LT. 0.)
!!$          xtte(:,:,H2O_index,:) = - xtm1(:,:,H2O_index,:) / time_step_len   
!!$       END WHERE
!!$       WHERE (xtm1(:,:,IWC_index,:)+xtte(:,:,IWC_index,:)*time_step_len .LT. 0.)
!!$          xtte(:,:,IWC_index,:) = - xtm1(:,:,IWC_index,:) / time_step_len   
!!$       END WHERE
!!$       WHERE (xtm1(:,:,CLWC_index,:)+xtte(:,:,CLWC_index,:)*time_step_len .LT. 0.)
!!$          xtte(:,:,CLWC_index,:) = - xtm1(:,:,CLWC_index,:) / time_step_len   
!!$       END WHERE
!!$    END IF
!!$#endif

  END SUBROUTINE clamscirrus_global_end
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamscirrus_free_memory

    IMPLICIT NONE

#ifdef ECHAM5
    DEALLOCATE (CIRRUS_tte)
#endif

  END SUBROUTINE clamscirrus_free_memory

!--------------------------------------------------------------------
END MODULE messy_clamscirrus_si
