MODULE mmd_mpi_wrapper

  USE      mpi
  USE      MMD_handle_communicator, ONLY: m_to_parent_comm, m_to_child_comm, &
                                          m_model_comm, m_model_rank
  USE      MMD_utilities,           ONLY: MMD_DP

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER, PARAMETER :: dp = MMD_DP

! INTERFACE section
  
  INTERFACE MMD_Send_to_Parent
     MODULE procedure MMD_Send_to_Parent_integer
     MODULE procedure MMD_Send_to_Parent_integer_i2
     MODULE procedure MMD_Send_to_Parent_real_r1
     MODULE procedure MMD_Send_to_Parent_real_r2
     MODULE procedure MMD_Send_to_Parent_real_r3
  END INTERFACE

  INTERFACE MMD_Recv_from_Parent
     MODULE procedure MMD_Recv_from_Parent_integer
     MODULE procedure MMD_Recv_from_Parent_integer_i2
     MODULE procedure MMD_Recv_from_Parent_real_r1
     MODULE procedure MMD_Recv_from_Parent_real_r2
     MODULE procedure MMD_Recv_from_Parent_real_r3
  END INTERFACE

  INTERFACE MMD_Send_to_Child
     MODULE procedure MMD_Send_to_Child_integer
     MODULE procedure MMD_Send_to_Child_integer_i2
     MODULE procedure MMD_Send_to_Child_real_r1
     MODULE procedure MMD_Send_to_Child_real_r2
     MODULE procedure MMD_Send_to_Child_real_r3
  END INTERFACE

  INTERFACE MMD_Recv_from_Child
     MODULE procedure MMD_Recv_from_Child_integer
     MODULE procedure MMD_Recv_from_Child_integer_i2
     MODULE procedure MMD_Recv_from_Child_real_r1
     MODULE procedure MMD_Recv_from_Child_real_r2
     MODULE procedure MMD_Recv_from_Child_real_r3
  END INTERFACE

  INTERFACE MMD_Bcast
     MODULE procedure MMD_Bcast_Integer
     MODULE procedure MMD_Bcast_character
  END INTERFACE

  INTERFACE MMD_Inter_Bcast
     MODULE procedure MMD_Inter_Bcast_Integer_1
     MODULE procedure MMD_Inter_Bcast_Character_1 ! um_ak_20150413   
  END INTERFACE

  PUBLIC :: MMD_Send_to_Parent, MMD_Recv_from_Parent
  PUBLIC :: MMD_Send_to_Child,  MMD_Recv_from_Child
  PUBLIC :: MMD_Bcast,          MMD_Inter_Bcast

 contains

  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! SEND TO /RECV FROM PARENT
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Parent_integer (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    INTEGER, DIMENSION(:), INTENT(IN)         :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Parent_rank
    INTEGER, INTENT(IN)                       :: tag 
    INTEGER, INTENT(OUT)                      :: ierr
    
    ierr = 0
    CALL MPI_Send (buf, n, MPI_INTEGER, Parent_rank, tag, m_to_parent_comm &
         , ierr)
    
    RETURN

 END SUBROUTINE MMD_Send_to_Parent_integer
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Parent_integer (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    INTEGER, DIMENSION(:), INTENT(OUT)        :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Parent_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_INTEGER, Parent_rank, tag, m_to_parent_comm, &
         MPI_STATUS_IGNORE, ierr)

    RETURN

 END SUBROUTINE MMD_Recv_from_Parent_integer
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Parent_integer_i2 (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    INTEGER, DIMENSION(:,:), INTENT(IN)       :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Parent_rank
    INTEGER, INTENT(IN)                       :: tag 
    INTEGER, INTENT(OUT)                      :: ierr
    
    ierr = 0
    CALL MPI_Send (buf, n, MPI_INTEGER, Parent_rank, tag, m_to_parent_comm&
         , ierr)
    
    RETURN
    
  END SUBROUTINE MMD_Send_to_Parent_integer_i2
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Parent_integer_i2 (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    INTEGER, DIMENSION(:,:), INTENT(OUT)      :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Parent_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_INTEGER, Parent_rank, tag, m_to_parent_comm, &
         MPI_STATUS_IGNORE, ierr)

    RETURN

  END SUBROUTINE MMD_Recv_from_Parent_integer_i2
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Parent_real_r1 (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    REAL(KIND=DP), DIMENSION(:), INTENT(IN)   :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Parent_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Send (buf, n, MPI_DOUBLE_PRECISION, Parent_rank, tag &
         , m_to_parent_comm, ierr)

    RETURN

 END SUBROUTINE MMD_Send_to_Parent_real_r1
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Parent_real_r1 (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    REAL(KIND=DP), DIMENSION(:), INTENT(OUT)  :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Parent_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_DOUBLE_PRECISION, Parent_rank, tag&
         , m_to_parent_comm, MPI_STATUS_IGNORE, ierr)

    RETURN
  END SUBROUTINE MMD_Recv_from_Parent_real_r1
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Parent_real_r2 (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    REAL(KIND=DP), DIMENSION(:,:), INTENT(IN) :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Parent_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Send (buf, n, MPI_DOUBLE_PRECISION, Parent_rank, tag &
         , m_to_parent_comm, ierr)

    RETURN

 END SUBROUTINE MMD_Send_to_Parent_real_r2
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Parent_real_r2 (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    REAL(KIND=DP), DIMENSION(:,:),INTENT(OUT) :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Parent_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_DOUBLE_PRECISION, Parent_rank, tag &
         , m_to_parent_comm, MPI_STATUS_IGNORE, ierr)

    RETURN

 END SUBROUTINE MMD_Recv_from_Parent_real_r2
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Parent_real_r3 (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: buf
    INTEGER, INTENT(IN)                         :: n
    INTEGER, INTENT(IN)                         :: Parent_rank
    INTEGER, INTENT(IN)                         :: tag
    INTEGER, INTENT(OUT)                        :: ierr

    ierr = 0
    CALL MPI_Send (buf, n, MPI_DOUBLE_PRECISION, Parent_rank, tag &
         , m_to_parent_comm, ierr)

    RETURN

  END SUBROUTINE MMD_Send_to_Parent_real_r3
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Parent_real_r3 (buf, n, Parent_rank, tag, ierr)

    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    REAL(KIND=DP), DIMENSION(:,:,:),INTENT(OUT) :: buf
    INTEGER, INTENT(IN)                         :: n
    INTEGER, INTENT(IN)                         :: Parent_rank
    INTEGER, INTENT(IN)                         :: tag
    INTEGER, INTENT(OUT)                        :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_DOUBLE_PRECISION, Parent_rank, tag &
         , m_to_parent_comm, MPI_STATUS_IGNORE, ierr)

    RETURN

  END SUBROUTINE MMD_Recv_from_Parent_real_r3
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! SENT TO / RECV FROM CHILD
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Child_integer (Child_id, buf, n, Child_rank &
                                       , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    INTEGER, INTENT(IN)                       :: Child_id
    INTEGER, DIMENSION(:), INTENT(IN)         :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Child_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Send (buf, n, MPI_INTEGER, Child_rank, tag &
         , m_to_child_comm(Child_id), ierr)

    RETURN

  END SUBROUTINE MMD_Send_to_Child_integer
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Child_integer (Child_id, buf, n, Child_rank &
                                         , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    INTEGER, INTENT(IN)                       :: Child_id
    INTEGER, DIMENSION(:), INTENT(INOUT)      :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Child_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_INTEGER, Child_rank, tag &
         , m_to_child_comm(Child_id), MPI_STATUS_IGNORE, ierr)

    RETURN

  END SUBROUTINE MMD_Recv_from_Child_integer
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Child_integer_i2 (Child_id, buf, n, Child_rank &
                                           , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    INTEGER, INTENT(IN)                       :: Child_id
    INTEGER, DIMENSION(:,:), INTENT(IN)       :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Child_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Send (buf, n, MPI_INTEGER, Child_rank, tag &
         , m_to_child_comm(Child_id), ierr)

    RETURN

  END SUBROUTINE MMD_Send_to_Child_integer_i2
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Child_integer_i2 (Child_id, buf, n, Child_rank &
                                            , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    INTEGER, INTENT(IN)                       :: Child_id
    INTEGER, DIMENSION(:,:), INTENT(INOUT)    :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Child_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_INTEGER, Child_rank, tag &
         , m_to_child_comm(Child_id), MPI_STATUS_IGNORE, ierr)

    RETURN

  END SUBROUTINE MMD_Recv_from_Child_integer_i2
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Child_real_r1 (Child_id, buf, n, Child_rank &
                                       , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    INTEGER, INTENT(IN)                       :: Child_id
    REAL(KIND=DP), DIMENSION(:), INTENT(IN)   :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Child_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Send (buf, n, MPI_DOUBLE_PRECISION, Child_rank, tag &
         , m_to_child_comm(Child_id), ierr)

    RETURN

  END SUBROUTINE MMD_Send_to_Child_real_r1
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Child_real_r1 (Child_id, buf, n, Child_rank &
                                         , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    INTEGER, INTENT(IN)                       :: Child_id
    REAL(KIND=DP), DIMENSION(:), INTENT(INOUT):: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Child_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_DOUBLE_PRECISION, Child_rank, tag &
         , m_to_child_comm(Child_id), MPI_STATUS_IGNORE, ierr)

    RETURN

  END SUBROUTINE MMD_Recv_from_Child_real_r1
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Child_real_r2 (Child_id, buf, n, Child_rank &
                                       , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    INTEGER, INTENT(IN)                       :: Child_id
    REAL(KIND=DP), DIMENSION(:,:), INTENT(IN) :: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Child_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Send (buf, n, MPI_DOUBLE_PRECISION, Child_rank, tag &
         , m_to_child_comm(Child_id), ierr)

    RETURN

  END SUBROUTINE MMD_Send_to_Child_real_r2
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Child_real_r2 (Child_id, buf, n, Child_rank &
                                         , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    INTEGER, INTENT(IN)                       :: Child_id
    REAL(KIND=DP), DIMENSION(:,:), INTENT(OUT):: buf
    INTEGER, INTENT(IN)                       :: n
    INTEGER, INTENT(IN)                       :: Child_rank
    INTEGER, INTENT(IN)                       :: tag
    INTEGER, INTENT(OUT)                      :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_DOUBLE_PRECISION, Child_rank, tag &
          , m_to_child_comm(Child_id), MPI_STATUS_IGNORE, ierr)

    RETURN

  END SUBROUTINE MMD_Recv_from_Child_real_r2
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Send_to_Child_real_r3 (Child_id, buf, n, Child_rank &
                                       , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Send

    INTEGER, INTENT(IN)                         :: Child_id
    REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: buf
    INTEGER, INTENT(IN)                         :: n
    INTEGER, INTENT(IN)                         :: Child_rank
    INTEGER, INTENT(IN)                         :: tag
    INTEGER, INTENT(OUT)                        :: ierr

    ierr = 0
    CALL MPI_Send (buf, n, MPI_DOUBLE_PRECISION, Child_rank, tag &
         , m_to_child_comm(Child_id),  ierr)

    RETURN
  END SUBROUTINE MMD_Send_to_Child_real_r3
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Recv_from_Child_real_r3 (Child_id, buf, n, Child_rank &
                                         , tag, ierr)
    IMPLICIT     NONE

    EXTERNAL :: MPI_Recv

    INTEGER, INTENT(IN)                         :: Child_id
    REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT):: buf
    INTEGER, INTENT(IN)                         :: n
    INTEGER, INTENT(IN)                         :: Child_rank
    INTEGER, INTENT(IN)                         :: tag
    INTEGER, INTENT(OUT)                        :: ierr

    ierr = 0
    CALL MPI_Recv (buf, n, MPI_DOUBLE_PRECISION, Child_rank, tag   &
         , m_to_child_comm(Child_id), MPI_STATUS_IGNORE, ierr)

    RETURN

  END SUBROUTINE MMD_Recv_from_Child_real_r3
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! INTEGER B_CAST 
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Bcast_Integer (buf, root_pe, comm, ierr)

    IMPLICIT     NONE

    EXTERNAL  :: MPI_Bcast
    INTRINSIC :: PRESENT

    INTEGER, INTENT(INOUT)           :: buf
    INTEGER, INTENT(IN)              :: root_pe
    INTEGER, INTENT(IN), OPTIONAL    :: comm
    INTEGER, INTENT(OUT),OPTIONAL    :: ierr

    !-- local variables
    INTEGER                          :: myComm
    INTEGER                          :: myErr

    IF(PRESENT (comm))  then
       myComm = comm
    ELSE
       myComm = m_model_comm
    END IF

    CALL MPI_Bcast (buf, 1, MPI_Integer, root_pe, myComm, myErr)

    IF(PRESENT (ierr))  then
       ierr = myErr
    END IF

    RETURN

  END SUBROUTINE MMD_Bcast_Integer
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Bcast_character (buf, root_pe, comm, ierr)

    IMPLICIT     NONE

    EXTERNAL  :: MPI_Bcast
    INTRINSIC :: LEN, PRESENT

    CHARACTER(LEN=*), INTENT(INOUT)             :: buf
    INTEGER, INTENT(IN)                         :: root_pe
    INTEGER, INTENT(IN), OPTIONAL               :: comm
    INTEGER, INTENT(OUT),OPTIONAL               :: ierr
    !-- local variables
    INTEGER                                     :: myComm
    INTEGER                                     :: myErr

    IF(PRESENT (comm))  then
       myComm = comm
    ELSE
       myComm = m_model_comm
    END IF
    
    CALL MPI_Bcast (buf, len(buf), MPI_Character, root_pe, myComm, myErr)
    
    IF(PRESENT (ierr))  then
       ierr = myErr
    END IF

    RETURN

 END SUBROUTINE MMD_Bcast_character
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Inter_Bcast_Integer_1 (buf, sender, Child_id, ierr)
    ! um_ak_20150413 sender added

    IMPLICIT     NONE

    EXTERNAL  :: MPI_Bcast
    INTRINSIC :: PRESENT, SIZE
    

    INTEGER, INTENT(INOUT),DIMENSION(:)         :: buf
    LOGICAL, INTENT(IN)                         :: sender ! um_ak_20150413
    INTEGER, INTENT(IN), OPTIONAL               :: Child_id
    INTEGER, INTENT(OUT),OPTIONAL               :: ierr

    !-- local variables
    INTEGER                                     :: myComm
    INTEGER                                     :: myErr
    INTEGER                                     :: root_pe

    ! model PE rank = 0 broadcast to all remote model PEs

    IF(PRESENT (Child_id))  then
       myComm  = m_to_child_comm(Child_id)
    ELSE
       myComm  = m_to_parent_comm
    END IF

    IF (sender) THEN
       IF(m_model_rank == 0)  then
          root_pe = MPI_ROOT
       ELSE
          root_pe = MPI_PROC_NULL
       END IF
    ELSE
       root_pe = 0
    END IF
    
    CALL MPI_Bcast (buf, size(buf), MPI_INTEGER, root_pe, myComm, myErr)
    
    IF(PRESENT (ierr))  then
       ierr = myErr
    END IF

    RETURN

 END SUBROUTINE MMD_Inter_Bcast_Integer_1
  ! --------------------------------------------------------------------------

  ! um_ak_20150413+
  ! --------------------------------------------------------------------------
  SUBROUTINE MMD_Inter_Bcast_Character_1 (buf, sender, Child_id, ierr)

    IMPLICIT     NONE

    EXTERNAL  :: MPI_Bcast
    INTRINSIC :: PRESENT, SIZE
    

    CHARACTER(LEN=*), INTENT(INOUT)         :: buf
    LOGICAL, INTENT(IN)                     :: sender ! um_ak_20150413
    INTEGER, INTENT(IN),OPTIONAL            :: Child_id
    INTEGER, INTENT(OUT),OPTIONAL           :: ierr

    !-- local variables
    INTEGER                                     :: myComm
    INTEGER                                     :: myErr
    INTEGER                                     :: root_pe

    !   PE 0 of the model to all remote model PEs

    IF(PRESENT (Child_id))  then
       myComm  = m_to_child_comm(Child_id)
    ELSE
       myComm  = m_to_parent_comm
    END IF

    IF (sender) THEN
       IF(m_model_rank == 0)  then
          root_pe = MPI_ROOT
       ELSE
          root_pe = MPI_PROC_NULL
       END IF
    ELSE
       root_pe = 0
    END IF
    
    CALL MPI_Bcast (buf, len(buf), MPI_CHARACTER, root_pe, myComm, myErr)
    
    IF(PRESENT (ierr))  then
       ierr = myErr
    END IF

    RETURN

  END SUBROUTINE MMD_Inter_Bcast_Character_1
  ! --------------------------------------------------------------------------

  ! um_ak_20150413-

END MODULE mmd_mpi_wrapper
