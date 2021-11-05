module mmd_handle_communicator

! Handle Communicator for the Multi Model driver
!
! Author: Klaus Ketelsen, MPICH, Oct 2008
! Author: Astrid Kerkweg, UniMz, May 2009/2010/2016 slight modifications

  USE  mpi
  USE  mmd_utilities,     ONLY: MMD_STATUS_OK, MMD_MAX_MODEL,           &
                                MMD_ParentIsECHAM, MMD_ParentIsCOSMO
  IMPLICIT NONE
  PRIVATE
  SAVE

  ! Define Types
  TYPE MMD_layout
    CHARACTER(LEN=7)        :: name
    INTEGER                 :: Parent_Id
    INTEGER                 :: npe
  END TYPE MMD_layout

  ! RETURN status
  PUBLIC ::                            MMD_STATUS_OK
  INTEGER,PARAMETER,PUBLIC          :: MMD_ERROR_NPES = 1
  INTEGER,PARAMETER,PUBLIC          :: MMD_ERROR_MPI  = 2
  INTEGER,PARAMETER,PUBLIC          :: MMD_STATUS_OFF = 3 ! mz_pj_20081216

  ! Coupler Setup
  ! -------------------------------------
  ! Coupler Id of this model 
  INTEGER                                    :: m_my_CPL_Id 
  ! Number of Coupler in layout file (-1= standalone)
  INTEGER                                    :: m_NrOfCpl    
  ! Number of Models not coupled via MMD (e.g OASIS) 
  INTEGER                                    :: m_NrOfNonMMDModels = 0
  !Information of all coupler 
  TYPE(MMD_layout),dimension(MMD_MAX_MODEL)  :: m_couplers   

  ! MPI settings
  ! -------------------------------------
  ! Communicator of this model
  INTEGER,PUBLIC                             :: m_model_comm     = -1     
  ! Communicator to the parent
  INTEGER,PUBLIC                             :: m_to_parent_comm = -1 
  ! Communicator to the child(s)
  INTEGER,dimension(MMD_MAX_MODEL), PUBLIC   :: m_to_child_comm  = -1     
  INTEGER                                    :: m_world_rank = -1 
  INTEGER                                    :: m_world_npes = -1 
  INTEGER,PUBLIC                             :: m_model_rank = -1 
  INTEGER,PUBLIC                             :: m_model_npes = -1 
  ! Parent Type (1 ECHAM, 2 COSMO)
  INTEGER,PUBLIC                             :: m_ParentType = -1       

  ! Indicates this PE is parent for Child model with Id ... 
  INTEGER,DIMENSION(:),POINTER,PUBLIC :: MMD_Parent_for_Child => NULL()

  ! Interface section
  ! -------------------------------------

  INTERFACE MMD_get_model_communicator 
    MODULE PROCEDURE MMD_get_model_communicator
    MODULE PROCEDURE MMD_cag_model_communicator    ! create and get
  END INTERFACE MMD_get_model_communicator 

  INTERFACE MMD_Print_Error_Message
    MODULE PROCEDURE MMD_Print_Error_Message
  END INTERFACE MMD_Print_Error_Message

  INTERFACE MMD_FreeMem_communicator
    MODULE PROCEDURE MMD_FreeMem_communicator
  END INTERFACE MMD_FreeMem_communicator

  PUBLIC :: MMD_get_model_communicator  
  PUBLIC :: MMD_FreeMem_communicator
  PUBLIC :: MMD_Print_Error_Message 

  PUBLIC :: MMD_instance_name               ! ub_ak_20180613
  PUBLIC :: MMD_numberofcoupledmodels       ! ub_ak_20190527
  PUBLIC :: MMD_numberofnonMMDcoupledmodels ! ub_ak_20190527

 CONTAINS

  ! ----------------------------------------------------------------------
   ! ub_ak_20180613+
  !SUBROUTINE MMD_cag_model_communicator (comm, MMD_status)
  SUBROUTINE MMD_cag_model_communicator (comm, MMD_status, myname)
   ! ub_ak_20180613-

    IMPLICIT     NONE

    INTRINSIC :: ADJUSTL, TRIM
    EXTERNAL  :: MPI_Bcast

    INTEGER,INTENT(OUT)                     :: comm
    INTEGER,INTENT(OUT)                     :: MMD_status
    CHARACTER(LEN=7), INTENT(OUT), OPTIONAL :: myname     ! ub_ak_20180613

    ! LOCAL
    INTEGER                             :: i,istat
    INTEGER,DIMENSION(MMD_MAX_MODEL+1)  :: start_PE
    INTEGER                             :: m_my_CPL_rank
    INTEGER                             :: tag, ChildCount
    ! I am active parent for this child ID
    INTEGER,DIMENSION(MMD_MAX_MODEL)    :: activeParent

    MMD_status   = MMD_STATUS_OK
    comm         = -1
    m_my_CPL_Id  = -1
    ChildCount  = 0
    activeParent = -1
    start_PE(:)  = 0     ! mz_pj_20081216

    CALL  MPI_Comm_rank (MPI_COMM_WORLD, m_world_rank, istat)
    IF (istat /=  MPI_SUCCESS) THEN
       istat = MMD_ERROR_MPI
       RETURN
    ENDIF
    CALL  MPI_Comm_size (MPI_COMM_WORLD, m_world_npes, istat)
    IF (istat /=  MPI_SUCCESS) THEN
       istat = MMD_ERROR_MPI
       RETURN
    ENDIF

    IF (m_world_rank == 0) THEN ! only master parent PE 0 reads

      CALL read_coupling_layout (MMD_status)

      IF (MMD_status /= MMD_STATUS_OFF) THEN ! mz_pj_20081216
         ! Compute Start PE of every model
         start_PE(1) = 0
         DO i=2,m_NrOfCpl+1
            start_pe(i) = &
                 start_PE(i-1) + m_couplers(i-1)%npe
         END DO
         IF (start_pe(m_NrOfCpl+1) /= m_world_npes)   THEN
            MMD_status = MMD_ERROR_NPES
            RETURN
         END IF
      END IF ! mz_pj_20081216
   END IF

   CALL MPI_Bcast (m_NrOfCpl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   IF (istat /=  MPI_SUCCESS) THEN
      istat = MMD_ERROR_MPI
      RETURN
   ENDIF
   CALL MPI_Bcast ( start_PE, m_NrOfCpl+1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   IF (istat /=  MPI_SUCCESS) THEN
      istat = MMD_ERROR_MPI
      RETURN
   ENDIF

   ! mz_pj_20081216+
   CALL MPI_Bcast (MMD_status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   IF (istat /=  MPI_SUCCESS) THEN
      istat = MMD_ERROR_MPI
      RETURN
   ENDIF
   IF (MMD_status == MMD_STATUS_OFF) THEN
      comm = MPI_COMM_WORLD
      allocate(MMD_Parent_for_Child(ChildCount))
      MMD_status = MMD_STATUS_OK
      m_NrOfCpl = -1
      m_NrOfNonMMDModels = 0
      RETURN
   END IF
   ! mz_pj_20081216-

   DO i=1,m_NrOfCpl
      IF (m_world_rank >= start_PE(i) .and. m_world_rank < start_PE(i+1) ) THEN
         m_my_CPL_Id = i   
         EXIT
      END IF
   END DO
   m_my_CPL_rank = m_world_rank-start_PE(i)

   ! MPI_COMM_WORLD is the communicator for ALL models (MPI-1 approach)
   ! The communictors for the individual models are created by MPI_Comm_split
   ! The color of the model is represented by the Coupler Id
 
    CALL MPI_Comm_split ( MPI_COMM_WORLD, m_my_CPL_Id, m_my_CPL_rank &
                        , comm, istat)
    IF (istat /= MPI_SUCCESS) THEN
      MMD_status = MMD_ERROR_MPI
      RETURN
    END IF

    ! Get size and rank of the communicator of the model running on this PE
    CALL  MPI_Comm_rank (comm, m_model_rank, istat)
    IF (istat /=  MPI_SUCCESS) THEN
       istat = MMD_ERROR_MPI
       RETURN
    ENDIF
    CALL  MPI_Comm_size (comm, m_model_npes, istat)
    IF (istat /=  MPI_SUCCESS) THEN
       istat = MMD_ERROR_MPI
       RETURN
    ENDIF

    ! Start_Pe of every model brodcasts its Communicator and Parent ID
    ! Pe 0  brodcasts the Parent ID
    DO i=1,m_NrOfCpl
       CALL MPI_Bcast(m_couplers(i)%Parent_Id, 1, MPI_INTEGER,   0   &
            , MPI_COMM_WORLD, istat)
       IF (istat /=  MPI_SUCCESS) THEN
          istat = MMD_ERROR_MPI
          RETURN
       ENDIF
       CALL MPI_Bcast (m_couplers(i)%name,     7, MPI_CHARACTER, 0   &
            , MPI_COMM_WORLD, istat)
       IF (istat /=  MPI_SUCCESS) THEN
          istat = MMD_ERROR_MPI
          RETURN
       ENDIF
    END DO

    m_model_comm = comm

    ! create Intercommunicator to parent and childs
    ! MPI_Intercomm_create creates an intercommunicator 
    ! between 2 groups of different colors
    ! The grouping with done prior with MPI_Comm_split

    m_ParentType = -1                            !Default and NO Parent
    DO i=2,m_NrOfCpl
      IF (m_couplers(i)%Parent_Id == m_my_CPL_Id)   THEN
        tag = 500+i
        CALL MPI_Intercomm_create (comm, 0, MPI_COMM_WORLD, start_pe(i), &
                                 tag, m_to_child_comm(i), istat)
        IF (istat /=  MPI_SUCCESS) THEN
           istat = MMD_ERROR_MPI
           RETURN
        ENDIF
        childCount = childCount+1
        activeParent(i) = 1
      ELSE IF (i == m_my_CPL_Id  &
           ! could be more than one root (e.g. ensembles)
           ! ub_ak_20181004
           .AND. m_couplers(i)%Parent_Id > 0)   THEN       ! 
         ! ub_ak_20180410-
        tag = 500+i
        CALL MPI_Intercomm_create (comm, 0, MPI_COMM_WORLD &
             , start_pe(m_couplers(i)%Parent_Id), tag, m_to_parent_comm, istat)
        IF (istat /=  MPI_SUCCESS) THEN
           istat = MMD_ERROR_MPI
           RETURN
        ENDIF
        IF (m_couplers(i)%Parent_Id == 1)   THEN         
           ! um_ak_20090513+
           !              m_ParentType = MMD_ParentIsECHAM
           ! ub_ak_20190122+
!           IF (ADJUSTL(TRIM(m_couplers(1)%name)) == 'echam') THEN
           IF (m_couplers(1)%name(1:5) == 'echam') THEN
           ! ub_ak_20190122-
              m_ParentType = MMD_ParentIsECHAM             !Parent is ECHAM
           ELSE
              m_ParentType = MMD_ParentIsCOSMO             !Parent is COSMO
           ENDIF
           ! um_ak_20090513-
        ELSE
           m_ParentType = MMD_ParentIsCOSMO               !Parent is COSMO
        END IF
        !kk  write(0,'(a,7i4)') 'child Part',m_world_rank,m_world_npes &
        !kk ,m_model_rank,m_model_npes,tag, start_pe(m_couplers(i)%Parent_Id) &
        !kk ,m_ParentType, m_couplers(i)%Parent_Id

      END IF
      IF (istat /= MPI_SUCCESS) THEN
        MMD_status = MMD_ERROR_MPI
        RETURN
      END IF
    END DO

    ALLOCATE(MMD_Parent_for_Child(ChildCount))
    ChildCount = 0
    DO i=2,m_NrOfCpl
      IF (activeParent(i) == 1)  THEN
        ChildCount = childCount+1
        MMD_Parent_for_Child(ChildCount) = i
      END IF
    END DO

    ! ub_ak_20190527+
    m_NrOfNonMMDModels = 0
    DO i=1,m_NrOfCpl
       IF (m_couplers(i)%Parent_Id < -5) THEN
          m_NrOfNonMMDModels = m_NrOfNonMMDModels + 1
       END IF
    END DO
    ! ub_ak_20190527-

    IF (PRESENT(myname)) myname=m_couplers(m_my_CPL_Id)%name ! ub_ak_20180613

    RETURN

  END SUBROUTINE MMD_cag_model_communicator
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_get_model_communicator (comm)

    IMPLICIT     NONE

    INTEGER,INTENT(OUT)             :: comm

    comm = m_model_comm
    
    RETURN

  END SUBROUTINE MMD_get_model_communicator
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_Print_Error_Message (iu, MMD_status)

    IMPLICIT     NONE

    INTEGER,INTENT(IN)              :: iu
    INTEGER,INTENT(IN)              :: MMD_status

    SELECT CASE (MMD_status)
      CASE (MMD_ERROR_NPES)
        write(iu,*) 'MMD Status Number of PEs Inconsistent'
      CASE (MMD_ERROR_MPI)
        write(iu,*) 'MMD Status Error in MPI CALL'
    END SELECT

    RETURN

  END SUBROUTINE MMD_Print_Error_Message 
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_FreeMem_Communicator

    IMPLICIT NONE

    DEALLOCATE(MMD_Parent_for_Child)

  END SUBROUTINE MMD_FreeMem_Communicator
  !-------------------------------------------------------------------------

  ! Private SUBROUTINEs
  !-------------------------------------------------------------------------
  SUBROUTINE read_coupling_layout (MMD_status)

    IMPLICIT NONE
    
    INTRINSIC :: TRIM

    INTEGER,INTENT(INOUT)           :: MMD_status

    INTEGER                         :: i,iunit,istat
    CHARACTER(LEN=*), PARAMETER     :: fname = 'MMD_layout.nml'
    LOGICAL                         :: lex
    NAMELIST /CPL/ m_couplers

    m_NrOfCpl = 0
    iunit     = 345

    ! mz_pj_20081216+
    INQUIRE(file=TRIM(fname), exist=lex)
    IF (.NOT. lex) THEN
       MMD_status = MMD_STATUS_OFF
       RETURN
    END IF
    ! mz_pj_20081216-

    OPEN(iunit,file=TRIM(fname))
    READ(iunit, NML=CPL, IOSTAT=istat)
    IF (istat /= 0) THEN
       WRITE(*,*) '*** ERROR: READ ERROR in NAMELIST MMD_layout.nml !'
       CLOSE(iunit)
    ELSE
       WRITE(*,*) ' ... OK !'
    END IF

    DO i=1,MMD_MAX_MODEL
       IF (TRIM(m_couplers(i)%name) == '') EXIT
       m_NrOfCpl = i
    END DO
    CLOSE(iunit)

    RETURN

  END SUBROUTINE read_coupling_layout
  !-------------------------------------------------------------------------

  ! ub_ak_20180613+
  CHARACTER(LEN=7) FUNCTION MMD_instance_name()

    IMPLICIT NONE

    MMD_instance_name = m_couplers(m_my_CPL_Id)%name

  END FUNCTION MMD_instance_name
  ! ub_ak_20180613-

  ! ub_ak_20190527+
  INTEGER FUNCTION MMD_numberofcoupledmodels()

    IMPLICIT NONE

    MMD_numberofcoupledmodels = m_NrOfCpl

  END FUNCTION MMD_numberofcoupledmodels

  INTEGER FUNCTION MMD_numberofnonMMDcoupledmodels()

    IMPLICIT NONE

    MMD_numberofnonMMDcoupledmodels = m_NrOfNonMMDModels

  END FUNCTION MMD_numberofnonMMDcoupledmodels
  ! ub_ak_20190527-
END MODULE mmd_handle_communicator
