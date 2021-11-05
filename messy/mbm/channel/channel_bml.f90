PROGRAM CHANNEL

  USE channel_mem_bml
  USE messy_main_channel_bi
  USE messy_submodel_si

  IMPLICIT NONE
  INTRINSIC :: MOD

  ! 1 INITIALISATION PHASE
  ! 1.1 read CTRL and CPL namelists from 'channel.nml'
  ! 1.2 setup global attributes
  ! 1.3 setup dimensions, dimension-variables
  ! 1.4 setup representations
  CALL main_channel_initialize

  ! 1.5 create standard basemodel channel(s) / object(s)
  CALL main_channel_init_memory

  ! SI-1 create new channel(s) / object(s)
  CALL submodel_init_memory

  ! 1.6 check, if start or restart (netCDF restart file present)
  DO cy = 9999, 1, -1
     WRITE(cystr,'(i4.4)') cy
     INQUIRE(file = TRIM('restart_'//cystr//'_'//modstr//'.nc') &
          , exist = lresume)
     IF (lresume) EXIT
  END DO
  WRITE(*,*)    '###########################'
  IF (lresume) THEN
     WRITE(*,*) 'RESTART FROM CYCLE '//cystr
  ELSE
     WRITE(*,*) 'START'
  ENDIF
  WRITE(*,*)    '###########################'

  ! 1.7 initialize restart information
  CALL main_channel_setup

  ! 1.8 add channels and object references from CTRL namelist
  ! 1.9 initialize I/O and file timer from CPL namelist
  CALL main_channel_init_coupling(1)
  CALL main_channel_init_coupling(2)

  ! 1.10 read restart information for all available channels
  IF (lresume) CALL main_channel_read_restart

  ! 2. TIME INTEGRATION LOOP
  DO 
     NH = NH + 1
     IF (NH > NH_max) EXIT
     WRITE(*,*) 'TIME STEP: ', NH

     ! SI-2. update channel objects
     CALL submodel_global_end

     ! 2.1 UPDATE TIMER INFORMATION
     CALL main_channel_write_output(1)

     ! 2.2 WRITE OUTPUT
     CALL main_channel_write_output(2)

     ! 2.3 CHECK RESTART
     lresume = MOD(NH, i_restfreq) == 0
     IF (lresume) THEN
        CALL main_channel_write_restart
        EXIT
     ENDIF

  END DO

  ! 3. FINALISING PHASE
  ! SI-3. clean up additional submodel memory
  !       (not contained in channel / objects)
  CALL submodel_free_memory
  ! 3.1 clean up channel / object memory
  CALL main_channel_free_memory

  STOP

END PROGRAM CHANNEL
