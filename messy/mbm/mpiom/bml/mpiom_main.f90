!*******************************************************************************
!
! Authors:
!  Astrid Kerkweg,  UNI-MZ,         2010
!  Patrick Joeckel, DLR,            2010
!
!*****************************************************************************

PROGRAM mpiom

  USE messy_main_channel_bi,    ONLY: main_channel_read_restart
  USE messy_main_blather_bi,    ONLY: messy_blather_endfile_bi, info_bi
  USE messy_main_grid_def_mem_bi, ONLY: ngpblks, jrow
!!$  USE messy_main_tracer_bi,     ONLY: messy_tracer_beforeadv &
!!$                                    , messy_tracer_afteradv
  USE messy_main_timer,         ONLY: lbreak, lstop      &
                                    , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
  USE messy_main_mpi_bi,        ONLY: p_start, p_stop, p_pe

  IMPLICIT NONE

  ! ######################
  ! 1 INITIALISATION PHASE
  ! ######################

  CALL p_start("MPIOM")

  ! 1.1 SETUP OF BASEMODEL and TIMER
  CALL messy_setup(1)           ! setup MESSy TIMER

  ! 1.2 read CTRL and CPL namelists
  ! 1.3 setup global attributes
  ! 1.4 setup dimensions, dimension-variables
  ! 1.5 setup representations
  CALL messy_initialize         ! read submodel namelists

  ! 1.6 define new tracers
  CALL messy_new_tracer         ! define tracers

  ! 1.7 initialise memory (create channel and channel objects)
  CALL messy_init_memory        ! allocate submodel memory

  ! 1.9 couple submodels
  CALL messy_init_coupling      ! couple to other submodels

  ! 1.10 read restart information for all available channels
  CALL messy_read_restart

  ! 1.11 initialise tracer
  CALL messy_init_tracer        ! initialise tracer

  ! 1.12
  ! CALL messy_init_loop

  ! ################################
  ! 2. TIME LOOP / INTEGRATION PHASE
  ! ################################
  time_loop: DO WHILE (.NOT. lbreak)

!!$     CALL messy_tracer_beforeadv
!!$     ! CALL tracer advection here
!!$     CALL messy_tracer_afteradv

     ! READ/UPDATE BOUNDARY CONDITIONS
     CALL messy_time(1)
     !CALL messy_tendency_reset
     CALL messy_global_start   ! first entry point for submodels in time loop

     if (p_pe==0) THEN
        WRITE(*,*) '-------------------------------------------- '
        WRITE(*, '(1x,a,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
       'Time_loop:',YEAR,'-',MONTH,'-',DAY,' ',HOUR,':', MINUTE,':',SECOND
        WRITE(*,*) '-------------------------------------------- '
     endif

     region_loop: DO jrow=1, ngpblks
        !
        CALL messy_local_start ! first entry point for submodels in region loop
        !
        ! CALCULATE PROCESSES AND DIAGNOSTICS
        CALL messy_physc
        !
        CALL messy_local_end   ! last entry point for submodels in region loop
        !
     END DO region_loop

     CALL messy_global_end     ! last entry point for submodels in time loop

     ! WRITE OUTPUT
     CALL messy_write_output

     ! WRITE RESTART FILES
     CALL messy_write_restart

     ! step to next time step / set new dates
     CALL messy_time(2)

  END DO time_loop

  ! #################################
  ! 3. FINALISING PHASE / FREE MEMORY
  ! #################################
  CALL messy_free_memory  ! close output files, deallocate fields

  CALL messy_finalize

  IF (lstop) THEN
     CALL messy_blather_endfile_bi('SIMULATION FINISHED','mpiom_main')
  ELSE
     CALL info_bi('SIMULATION INTERRUPTED','mpiom_main')
  ENDIF

  CALL p_stop

  STOP

END PROGRAM mpiom
