MODULE messy_mmd2way_si

 USE messy_mmd2way

IMPLICIT NONE

  INTERFACE mmd2way_setup
     MODULE PROCEDURE mmd2way_setup
  END INTERFACE mmd2way_setup

  INTERFACE mmd2way_initialize
     MODULE PROCEDURE mmd2way_initialize
  END INTERFACE mmd2way_initialize

  INTERFACE mmd2way_init_memory
     MODULE PROCEDURE mmd2way_init_memory
  END INTERFACE mmd2way_init_memory

  INTERFACE mmd2way_init_coupling
     MODULE PROCEDURE mmd2way_init_coupling
  END INTERFACE mmd2way_init_coupling

  INTERFACE mmd2way_init_coupler
     MODULE PROCEDURE mmd2way_init_coupler
  END INTERFACE mmd2way_init_coupler

  INTERFACE mmd2way_init_loop
     MODULE PROCEDURE mmd2way_init_loop
  END INTERFACE mmd2way_init_loop

  INTERFACE mmd2way_global_start
     MODULE PROCEDURE mmd2way_global_start
  END INTERFACE mmd2way_global_start

  INTERFACE mmd2way_global_end
     MODULE PROCEDURE mmd2way_global_end
  END INTERFACE mmd2way_global_end

  INTERFACE mmd2way_write_output
     MODULE PROCEDURE mmd2way_write_output
  END INTERFACE mmd2way_write_output

  INTERFACE mmd2way_write_restart
     MODULE PROCEDURE mmd2way_write_restart
  END INTERFACE mmd2way_write_restart

  INTERFACE mmd2way_read_restart
     MODULE PROCEDURE mmd2way_read_restart
  END INTERFACE mmd2way_read_restart

  INTERFACE mmd2way_free_memory
     MODULE PROCEDURE mmd2way_free_memory
  END INTERFACE mmd2way_free_memory


CONTAINS

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_setup

    USE messy_mmd2way_child_si,   ONLY: mmd2way_child_setup

    IMPLICIT NONE

    CALL mmd2way_child_setup
   
    RETURN

  END SUBROUTINE mmd2way_setup

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_initialize

    USE messy_mmd2way_parent_si,   ONLY: mmd2way_parent_initialize

    IMPLICIT NONE

    CALL mmd2way_parent_initialize
   
    RETURN

  END SUBROUTINE mmd2way_initialize

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_init_memory(flag)

    USE messy_mmd2way_child_si,   ONLY: mmd2way_child_init_memory
    USE messy_mmd2way_parent_si,  ONLY: mmd2way_parent_init_memory

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    SELECT CASE(flag) 

    CASE(1)
       CALL mmd2way_parent_init_memory
    CASE(2)
       CALL mmd2way_child_init_memory
    END SELECT

    RETURN

  END SUBROUTINE mmd2way_init_memory

 !--------------------------------------------------------------

  SUBROUTINE mmd2way_init_coupler
    
    USE messy_mmd2way_child_si,   ONLY: mmd2way_child_init_memory

    IMPLICIT NONE

    CALL mmd2way_child_init_memory

  END SUBROUTINE mmd2way_init_coupler

!--------------------------------------------------------------

  SUBROUTINE mmd2way_init_coupling

    ! MESSy/SMIL
    USE messy_mmd2way_parent_si,   ONLY: mmd2way_parent_init_coupling 

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'mmd2way_init_coupling'

    CALL mmd2way_parent_init_coupling

    RETURN

  END SUBROUTINE mmd2way_init_coupling

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_init_loop

    USE messy_mmd2way_child_si,    ONLY: mmd2way_child_init_loop

    IMPLICIT NONE

    CALL mmd2way_child_init_loop
   
    RETURN

  END SUBROUTINE mmd2way_init_loop

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_global_start

#ifdef MESSYMMD
    ! SMIL
    USE messy_mmd2way_parent_si,   ONLY: mmd2way_parent_global_start
#ifdef ECHAM5
    USE messy_mmd2way_parent_si,   ONLY: mmd2way_parent_global_end &
                                       , lcpl_global_start
    ! MESSy/SMCL
    USE messy_main_timer,          ONLY: lstart
#endif

    IMPLICIT NONE

#ifdef ECHAM5
    IF (lcpl_global_start .AND. .NOT. lstart) THEN
       CALL mmd2way_parent_global_end(1)
       CALL mmd2way_parent_global_end(2)
    ENDIF
#endif

    CALL mmd2way_parent_global_start
#endif
    RETURN

  END SUBROUTINE mmd2way_global_start

  !--------------------------------------------------------------

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_global_end(flag)
#ifdef MESSYMMD
    ! MESSy/SMIL
    USE messy_mmd2way_parent_si,   ONLY: mmd2way_parent_global_end
#ifdef ECHAM5
    USE messy_mmd2way_parent_si,   ONLY: lcpl_global_start
#endif
    USE messy_mmd2way_child_si,    ONLY: mmd2way_child_global_end

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

#ifdef ECHAM5
    IF (.NOT. lcpl_global_start) THEN
#endif
    SELECT CASE(flag)
    CASE(1,3)
       CALL mmd2way_parent_global_end(flag)
    CASE(2)
       CALL mmd2way_child_global_end
       CALL mmd2way_parent_global_end(flag)
    END SELECT
#ifdef ECHAM5
    END IF
#endif
#else
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    IF (flag == 0) RETURN
#endif
    RETURN

  END SUBROUTINE mmd2way_global_end

  !--------------------------------------------------------------

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_write_output(flag)

    USE messy_mmd2way_child_si,    ONLY: mmd2way_child_write_output

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    CALL mmd2way_child_write_output(flag)
   
    RETURN

  END SUBROUTINE mmd2way_write_output

  !--------------------------------------------------------------

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_write_restart(flag)

    USE messy_mmd2way_child_si,    ONLY: mmd2way_child_write_restart

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    CALL mmd2way_child_write_restart(flag)
   
    RETURN

  END SUBROUTINE mmd2way_write_restart

  !--------------------------------------------------------------

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_read_restart(flag)

    USE messy_mmd2way_child_si,    ONLY: mmd2way_child_read_restart

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    CALL mmd2way_child_read_restart(flag)
   
    RETURN

  END SUBROUTINE mmd2way_read_restart

  !--------------------------------------------------------------

  !--------------------------------------------------------------

  SUBROUTINE mmd2way_free_memory

    USE messy_mmd2way_child_si,   ONLY: mmd2way_child_free_memory
    USE messy_mmd2way_parent_si,  ONLY: mmd2way_parent_free_memory

    IMPLICIT NONE

    CALL mmd2way_child_free_memory
    CALL mmd2way_parent_free_memory
   
    RETURN

  END SUBROUTINE mmd2way_free_memory

  !--------------------------------------------------------------

END MODULE messy_mmd2way_si
