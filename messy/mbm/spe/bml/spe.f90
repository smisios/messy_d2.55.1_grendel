! COLUMN MODEL FOR SOLAR PROTON EVENTS (SPE)
!
! Author:
!    Andreas Baumgaertner, MPICH, 2007-2010
! 
PROGRAM SPE

  ! SMIL
  USE messy_spe_box
  USE messy_main_timer_bi,      ONLY: main_timer_setup  &
                                    , main_timer_initialize &
                                    , main_timer_time
  USE messy_main_timer,         ONLY: lbreak, lresume, lstop, l_rerun      &
                                    , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
  USE messy_main_channel_bi, ONLY: main_channel_setup &
       , main_channel_initialize &
       , main_channel_init_memory &
       , main_channel_init_coupling &
       , main_channel_global_start &
       , main_channel_write_output &
       , main_channel_write_restart &
       , main_channel_free_memory

  IMPLICIT NONE
 
  ! LOCAL
  INTEGER :: nclock = 1

  ! ub_ak_20171109+
  !CALL main_timer_setup(1)
  !CALL main_timer_setup(2)
  !CALL main_channel_setup
  CALL main_timer_setup(1, nclock)
  CALL main_channel_setup
  CALL main_timer_setup(2, nclock)
  ! ub_ak_20171109-

  CALL main_timer_initialize
  CALL main_channel_initialize
  CALL spe_initialize

  CALL main_channel_init_memory  
  CALL spe_init_memory

  CALL main_channel_init_coupling(1)
  CALL main_channel_init_coupling(2)

 ! TIME LOOP 
  DO while (.NOT. lbreak)
     WRITE(*, '(1x,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a3'//&
       &',a4,i4,a8,i5,a3,i5,a10,i5)') &
       YEAR,'-',MONTH,'-',DAY,' ',HOUR,':', MINUTE,':',SECOND  
     ! READ/UPDATE BOUNDARY CONDITIONS
     CALL main_timer_time(1)
     CALL main_channel_global_start
     CALL spe_global_start
     ! CALCULATE PHYSICAL PROCESSES
     CALL spe_physc     

     CALL main_channel_write_output !(1) ! ub_ak_20190110

     !CALL main_channel_write_output(2)  ! ub_ak_20190110

     IF (l_rerun) CALL main_channel_write_restart

     ! step to next time step / set new dates
     CALL main_timer_time(2)
  END DO

  CALL spe_free_memory
  CALL main_channel_free_memory

 END PROGRAM SPE
