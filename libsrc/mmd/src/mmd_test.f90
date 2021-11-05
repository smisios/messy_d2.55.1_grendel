module mmd_test

! Test Data distribution for the Multi Model driver
!
! The parent stores unique information of every grid Point in test_array
! The information will be sent to the child 
! The child compares the information with its local information and aborts
! in case of no match

! Author: Klaus Ketelsen, MPICH, Dec 2008

  USE mpi
  USE MMD_utilities,             ONLY: MMD_STATUS_OK, MMD_MAX_MODEL, MMD_DP
  USE MMD_handle_communicator,   ONLY: m_model_comm, m_model_rank, m_model_npes
  USE MMD_MPI_wrapper,           ONLY: MMD_Send_to_Parent, MMD_Recv_from_Parent&
                                     , MMD_Bcast

  IMPLICIT NONE

  PRIVATE
  SAVE

  INTEGER, PARAMETER :: DP = MMD_DP
  INTEGER, PARAMETER :: MMD_STATUS_TEST_ERROR=99
  

  TYPE ArraySetup
    REAL(KIND=DP),DIMENSION(:,:,:,:,:), POINTER   :: global_array => NULL()
    ! Test array maximum dimension
    REAL(KIND=DP),DIMENSION(:,:,:,:),   POINTER   :: test_array_md => NULL()
    ! Test array with local PE dimension
    REAL(KIND=DP),DIMENSION(:,:,:,:),   POINTER   :: test_array => NULL()
  END TYPE ArraySetup

  TYPE(ArraySetup), DIMENSION(MMD_MAX_MODEL) :: Sv
  TYPE(ArraySetup)                           :: Cl

  INTEGER,PARAMETER                          :: MT = 8   !Number of Test Values
  CHARACTER(LEN=4)                           :: dim_order
  integer                                    :: nxMax,nyMax
  
! Interface section

  INTERFACE MMD_testP_Setup
    MODULE PROCEDURE MMD_testP_Setup
  END INTERFACE MMD_testP_Setup

  INTERFACE MMD_testC_Setup
    MODULE PROCEDURE MMD_testC_Setup
  END INTERFACE MMD_testC_Setup

  INTERFACE MMD_testP_Fill
    MODULE PROCEDURE MMD_testP_Fill
  END INTERFACE MMD_testP_Fill

  INTERFACE MMD_testP_FinishFill
    MODULE PROCEDURE MMD_testP_FinishFill
  END INTERFACE MMD_testP_FinishFill

  INTERFACE MMD_testP_GetTestPtr
    MODULE PROCEDURE MMD_testP_GetTestPtr
  END INTERFACE MMD_testP_GetTestPtr

  INTERFACE MMD_testC_GetTestPtr
    MODULE PROCEDURE MMD_testC_GetTestPtr
  END INTERFACE MMD_testC_GetTestPtr

  INTERFACE MMD_testC_Compare
    MODULE PROCEDURE MMD_testC_Compare
  END INTERFACE MMD_testC_Compare

  INTERFACE MMD_testC_FreeMem
    MODULE PROCEDURE MMD_testC_FreeMem
  END INTERFACE MMD_testC_FreeMem

  INTERFACE MMD_testP_FreeMem
    MODULE PROCEDURE MMD_testP_FreeMem
  END INTERFACE MMD_testP_FreeMem

  PUBLIC :: MMD_testP_Setup, MMD_testC_Setup 
  PUBLIC :: MMD_testP_Fill, MMD_testP_FinishFill
  PUBLIC :: MMD_testP_GetTestPtr,MMD_testC_GetTestPtr
  PUBLIC :: MMD_testC_Compare
  PUBLIC :: MMD_testP_FreeMem, MMD_testC_FreeMem

 CONTAINS

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testP_Setup (Id, nx,ny,cdim_order)

    IMPLICIT NONE

    EXTERNAL  :: MPI_ALLREDUCE

    INTEGER, INTENT(IN)                  :: Id
    INTEGER, INTENT(IN)                  :: nx
    INTEGER, INTENT(IN)                  :: ny
    CHARACTER(LEN=4), INTENT(IN)         :: cdim_order
    ! LOCAL
    INTEGER                              :: istat


    dim_order = cdim_order

!kk In ECHAM5, the number of Blocks (ny) must not be the same on all PE's

    CALL MPI_Allreduce (nx, nxMax, 1, MPI_Integer, MPI_Max, m_model_comm, istat)
    CALL MPI_Allreduce (ny, nyMax, 1, MPI_Integer, MPI_Max, m_model_comm, istat)

    IF(m_model_rank == 0)  THEN
       IF(dim_order == 'XZNY')   THEN
          ALLOCATE(Sv(Id)%global_array(nxMax,1,MT,nyMax,0:m_model_npes-1))
       ELSE
          ALLOCATE(Sv(Id)%global_array(nxMax,nyMax,MT,1,0:m_model_npes-1))
       ENDIF
    ELSE
       ALLOCATE(Sv(Id)%global_array(1,1,1,1,1))
    END IF

    IF(dim_order == 'XZNY')   THEN
      ALLOCATE(Sv(Id)%test_array_md(nxMax,1,MT,nyMax))
      ALLOCATE(Sv(Id)%test_array(nx,1,MT,ny))
    ELSE
      ALLOCATE(Sv(Id)%test_array_md(nxMax,nyMax,MT,1))
      ALLOCATE(Sv(Id)%test_array(nx,ny,MT,1))
    ENDIF

    Sv(Id)%global_array = -999._dp

    RETURN

  END SUBROUTINE MMD_testP_Setup
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testC_Setup (nx,ny)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny

    ALLOCATE(Cl%test_array(nx,ny,MT,1))

    RETURN

  END SUBROUTINE MMD_testC_Setup
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testP_Fill (Id, index_list,lon,lat)

    IMPLICIT NONE

    INTRINSIC :: REAL

    INTEGER,INTENT(IN)              :: Id
    INTEGER,DIMENSION(:),INTENT(IN) :: index_list
    REAL(KIND=DP),INTENT(IN)        :: lon
    REAL(KIND=DP),INTENT(IN)        :: lat

!   Fill Global test array (is called from PE0 only)

    IF(dim_order == 'XZNY')   THEN
      Sv(Id)%global_array(index_list(1),1,1,   index_list(2),index_list(6)) &
           = lon
      Sv(Id)%global_array(index_list(1),1,2,   index_list(2),index_list(6)) &
           = lat
      Sv(Id)%global_array(index_list(1),1,3:MT,index_list(2),index_list(6)) &
           = REAL(index_list(:),dp)
    ELSE
      Sv(Id)%global_array(index_list(1),index_list(2),1,   1,index_list(6)) &
           = lon
      Sv(Id)%global_array(index_list(1),index_list(2),2,   1,index_list(6)) &
           = lat
      Sv(Id)%global_array(index_list(1),index_list(2),3:MT,1,index_list(6)) &
           = REAL(index_list(:),dp)
    END IF

    RETURN

  END SUBROUTINE MMD_testP_Fill
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testP_FinishFill (Id)

    IMPLICIT NONE
    
    EXTERNAL  :: MPI_SCATTER
    INTRINSIC :: SIZE

    INTEGER,INTENT(IN) :: Id

    ! LOCAL
    INTEGER :: istat,isize

    ! Scatter Global array
    isize = SIZE(Sv(Id)%test_array_md) 


    call MPI_Scatter (Sv(Id)%global_array, isize, MPI_DOUBLE_PRECISION,   &
                      Sv(Id)%test_array_md,   isize, MPI_DOUBLE_PRECISION, 0 &
                      , m_model_comm, istat) 
    Sv(Id)%test_array = Sv(Id)%test_array_md(&
         1:SIZE(Sv(Id)%test_array,1),1:SIZE(Sv(Id)%test_array,2),      &
         1:SIZE(Sv(Id)%test_array,3),1:SIZE(Sv(Id)%test_array,4))

    DEALLOCATE (Sv(Id)%global_array)   !Global array is not needed any more
    !test_array_md was only a container with maximum dimension for Scatter
    DEALLOCATE (Sv(Id)%test_array_md)  

    RETURN

  END SUBROUTINE MMD_testP_FinishFill
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testP_GetTestPtr (Id, p, axis, ldim)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    INTEGER,INTENT(IN)                     :: Id  ! Child Id
    CHARACTER(LEN=4), INTENT(OUT)          :: axis
    INTEGER, DIMENSION(4), INTENT(OUT)     :: ldim

    REAL(KIND=DP),DIMENSION(:,:,:,:),POINTER :: p ! Return pointer of test_array
                                                  ! to caller
    axis = dim_order
    ldim = (/SIZE(Sv(Id)%test_array,1),SIZE(Sv(Id)%test_array,2) &
         ,SIZE(Sv(Id)%test_array,3),SIZE(Sv(Id)%test_array,4)/)

    p => Sv(Id)%test_array

    RETURN

  END SUBROUTINE MMD_testP_GetTestPtr
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testC_GetTestPtr (p, axis, ldim)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    REAL(KIND=DP),DIMENSION(:,:,:,:),POINTER  :: p
    CHARACTER(LEN=4), INTENT(OUT)             :: axis
    INTEGER, DIMENSION(4), INTENT(OUT)        :: ldim

    p => Cl%test_array
    axis = 'XYNZ'
    ldim = &
         (/SIZE(Cl%test_array,1),SIZE(Cl%test_array,2) &
         ,SIZE(Cl%test_array,3),SIZE(Cl%test_array,4)/)


    RETURN

  END SUBROUTINE MMD_testC_GetTestPtr
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testC_Compare (lon,lat,istat,PrintUnit)

    IMPLICIT NONE

    INTRINSIC :: ABS, NINT, PRESENT, SIZE

    REAL(KIND=DP),DIMENSION(:,:),INTENT(IN) :: lon  ! original Coarse longitudes
    REAL(KIND=DP),DIMENSION(:,:),INTENT(IN) :: lat  ! original Coarse longitudes
    INTEGER,INTENT(OUT)               :: istat      ! Error Status
    INTEGER,INTENT(IN),OPTIONAL       :: PrintUnit  ! Print Unit (Default = 6)

    !-- Local Variables
    INTEGER                                   :: i,j,nx,ny,ierr_ind,jerr_ind
    INTEGER,DIMENSION(SIZE(Cl%test_array,3))  :: iarray
    REAL(KIND=DP),PARAMETER                   :: eps=0.003_dp
    INTEGER                                   :: iu

    istat = MMD_STATUS_OK
    iu    = 6
    IF(PRESENT(PrintUnit)) iu = PrintUnit

    nx = SIZE(Cl%test_array,1)
    ny = SIZE(Cl%test_array,2)

    DO j=1,ny
      DO i=1,nx
        iarray = nint(Cl%test_array(i,j,:,1))
!       Check coordinates
        IF(ABS(Cl%test_array(i,j,1,1)-lon(i,j)) > eps .OR.            &
           ABS(Cl%test_array(i,j,2,1)-lat(i,j)) > eps) THEN
             istat    = MMD_STATUS_TEST_ERROR
             ierr_ind = i
             jerr_ind = j
             write(0,'(a,i4,8x,a,2i5,10x,2i5,a,4F12.7)')               &
                  'Coordinate Error ',m_model_rank,                   &
                  '   I ',i,j,iarray(3),iarray(4),                    &
                  '   C ',Cl%test_array(i,j,1:2,1),lon(i,j),lat(i,j)
             EXIT
        END IF 

        ! Check indices

        ! In the inner boundary areas (PE decomposition), it may happen that 
        ! the same cell is sent to different PEs. In this case, the indices 
        ! stored in test_array are not unique index checking can only be done,
        ! if child PEnr in test_array equals model rank
        IF(m_model_rank == iarray(7))   THEN
          IF(i /= iarray(5) .OR. j /= iarray(6))   THEN
             istat    = MMD_STATUS_TEST_ERROR
             ierr_ind = i
             jerr_ind = j
             write(0,'(a,i4,8x,a,2i5,10x,2i5,a,4F6.1)') 'Index Error ' &
                  ,m_model_rank,      &
                  '   I ',i,j,iarray(3),iarray(4),                         &
                  '   C ',Cl%test_array(i,j,1:2,1),lon(i,j),lat(i,j)
             EXIT
          END IF
        END IF
      END DO
      IF(istat /= MMD_STATUS_OK) EXIT
    END DO

!   In Case of Error, Print
    IF(istat /= MMD_STATUS_OK) THEN
      DO j=1,ny
        DO i=1,nx
          IF(i /= ierr_ind .AND. j == jerr_ind) write(iu,*) ' '  
          !Blank Line with error cell
          iarray = nint(Cl%test_array(i,j,:,1))
          IF(m_model_rank == iarray(7))   THEN
            write(iu,'(a,3i4,a,6i5,a,4F6.1)') 'compare ',m_model_rank,iarray(7)&
                 ,iarray(8), &
                 '   I ',i,j,iarray(5),iarray(6),iarray(3),iarray(4),      &
                 '   C ',Cl%test_array(i,j,1:2,1),lon(i,j),lat(i,j)
         ELSE
            write(iu,'(a,i4,8x,a,2i5,10x,2i5,a,4F6.1)') 'compare ',m_model_rank&
                 , '   I ',i,j,iarray(3),iarray(4),                          &
                 '   C ',Cl%test_array(i,j,1:2,1),lon(i,j),lat(i,j)
          END IF
          IF(i == ierr_ind .AND. j == jerr_ind) write(iu,*) ' '
        END DO
      END DO
    END IF

    RETURN

  END SUBROUTINE MMD_testC_Compare
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testP_FreeMem(ChildId)

    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: ChildId

    DEALLOCATE(Sv(ChildId)%test_array)

  END SUBROUTINE MMD_testP_FreeMem
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE MMD_testC_FreeMem

    IMPLICIT NONE

    DEALLOCATE(Cl%test_array)

  END SUBROUTINE MMD_testC_FreeMem
  !--------------------------------------------------------------------

END MODULE mmd_test

