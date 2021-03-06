MODULE mmd_parent

  USE mpi
  USE MMD_utilities,           ONLY: MMD_DP, ArrayDef_list, ExchDataDef &
                                   , MMD_STATUS_OK
  USE MMD_handle_communicator, ONLY: MMD_Parent_for_Child, m_model_rank
  USE MMD_MPI_wrapper,         ONLY: MMD_Send_to_Child, MMD_Recv_from_Child,&
                                     MMD_Bcast, MMD_Inter_Bcast

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER, PARAMETER                         :: dp = MMD_DP
  TYPE(ExchDataDef), DIMENSION(:), POINTER   :: Child  => NULL()
  TYPE(ExchDataDef), DIMENSION(:), POINTER   :: Parent =>NULL()

  TYPE bufferlen
     INTEGER, DIMENSION(:), POINTER          :: BufLen => NULL()
  END type bufferlen
  TYPE(bufferlen), DIMENSION(:), ALLOCATABLE :: ChldBL ! Child BufLen
  TYPE(bufferlen), DIMENSION(:), ALLOCATABLE :: ParBL! Parent BufLen

  LOGICAL                                    :: ltwoway = .FALSE.

  PUBLIC :: MMD_Parent_for_Child
  PUBLIC :: MMD_STATUS_OK

  ! INTERFACE section
  INTERFACE MMD_P_Allocate_Child
    MODULE PROCEDURE MMD_P_Allocate_Child
  END INTERFACE MMD_P_Allocate_Child

  INTERFACE MMD_P_Init
    MODULE PROCEDURE  MMD_P_Init
  END INTERFACE MMD_P_Init

  INTERFACE MMD_P_Set_Indexlist
    MODULE PROCEDURE MMD_P_Set_Indexlist
  END INTERFACE MMD_P_Set_Indexlist

  INTERFACE MMD_P_Get_ParIndexList
    MODULE PROCEDURE MMD_P_Get_ParIndexList
  END INTERFACE MMD_P_Get_ParIndexList

  INTERFACE MMD_P_Set_ParDataArray_Name 
    MODULE PROCEDURE MMD_P_Set_ParDataArray_Name
  END INTERFACE MMD_P_Set_PArDataArray_Name

  INTERFACE MMD_P_Set_ParDataArray_EndList
    MODULE PROCEDURE MMD_P_Set_ParDataArray_EndList
  END INTERFACE MMD_P_Set_ParDataArray_EndList
  
  INTERFACE MMD_P_Set_ParDataArray
    MODULE PROCEDURE MMD_P_Set_ParDataArray
  END INTERFACE MMD_P_Set_ParDataArray

  INTERFACE MMD_P_GetNextParArray
    MODULE PROCEDURE MMD_P_GetNextParArray
  END INTERFACE MMD_P_GetNextParArray

  INTERFACE MMD_P_Get_DataArray_Name
    MODULE PROCEDURE MMD_P_Get_DataArray_Name
  END INTERFACE MMD_P_Get_DataArray_Name

  INTERFACE MMD_P_GetNextArray
    MODULE PROCEDURE MMD_P_GetNextArray
  END INTERFACE MMD_P_GetNextArray

  INTERFACE MMD_P_Send_Repr
   MODULE PROCEDURE MMD_P_Send_Repr
  END INTERFACE MMD_P_Send_Repr

  INTERFACE MMD_P_Set_DataArray
    MODULE PROCEDURE MMD_P_Set_DataArray
  END INTERFACE MMD_P_Set_DataArray

  INTERFACE MMD_P_setInd_and_AllocMem
    MODULE PROCEDURE MMD_P_setInd_and_AllocMem
  END INTERFACE MMD_P_setInd_and_AllocMem

  INTERFACE MMD_P_FillBuffer
   MODULE PROCEDURE MMD_P_FillBuffer
  END INTERFACE MMD_P_FillBuffer

  INTERFACE MMD_P_GetBuffer
    MODULE PROCEDURE MMD_P_GetBuffer
  END INTERFACE MMD_P_GetBuffer

  INTERFACE MMD_P_FreeMem
     MODULE PROCEDURE MMD_P_FreeMem
  END INTERFACE MMD_P_FreeMem

  ! PUBLIC section
  PUBLIC :: MMD_P_Allocate_Child
  PUBLIC :: MMD_P_Init, MMD_P_Set_Indexlist, MMD_P_GetNextArray
  PUBLIC :: MMD_P_Set_DataArray, MMD_P_Get_DataArray_Name, MMD_P_Send_Repr

  PUBLIC :: MMD_P_Set_ParDataArray_Name, MMD_P_Set_ParDataArray_EndList
  PUBLIC :: MMD_P_GetNextParArray,       MMD_P_Set_ParDataArray
  PUBLIC :: MMD_P_Get_ParIndexList,      MMD_P_GetBuffer

  PUBLIC :: MMD_P_setInd_and_AllocMem,   MMD_P_FillBuffer
  PUBLIC :: MMD_Send_to_Child,           MMD_Recv_from_Child
  PUBLIC :: MMD_Inter_Bcast,             MMD_P_FreeMem

 CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_Allocate_Child(NUMChildren, l2way)

    IMPLICIT NONE

    INTRINSIC :: PRESENT

    INTEGER, INTENT(IN)           :: NumChildren
    LOGICAL, INTENT(IN), OPTIONAL :: l2way
    ! LOCAL
    INTEGER                       :: Id

    IF (PRESENT(l2way)) ltwoway = l2way

    ALLOCATE(Child(NumChildren))
    ALLOCATE(ChldBL (NumChildren))
    
    DO Id=1,NumChildren
       NULLIFY(Child(Id)%Ar)
       NULLIFY(Child(Id)%ArrayStart)
    ENDDO

    IF (ltwoway) THEN
       ALLOCATE(Parent(NumChildren))  
       ALLOCATE(ParBL(NumChildren))     
       DO Id=1,NumChildren
          NULLIFY(Parent(Id)%Ar)        
          NULLIFY(Parent(Id)%ArrayStart)
       ENDDO
    END IF

  END SUBROUTINE MMD_P_Allocate_Child
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_Init (ChildId, Id)

    USE MMD_handle_communicator, ONLY: m_model_comm, m_to_child_comm
 
    IMPLICIT NONE

    EXTERNAL  :: MMDc_P_Init

    INTEGER, INTENT(IN)           :: ChildId
    INTEGER, INTENT(IN)           :: Id

    ! LOCAL
    INTEGER                       :: ip
    
    Child(Id)%RMId =ChildId

    CALL MMDc_P_Init (ChildId, m_model_comm &
         , m_to_child_comm(ChildId), Child(Id)%inter_npes)

    ! allocate structure for Child PE information 
    ALLOCATE(Child(Id)%PEs(0:Child(Id)%inter_npes-1))
    ALLOCATE(ChldBL(Id)%BufLen(0:Child(Id)%inter_npes-1))

    DO ip =0, Child(Id)%inter_npes-1
       Child(Id)%PEs(ip)%NrEle = 0
       NULLIFY(Child(Id)%PEs(ip)%locInd)

       ChldBL(Id)%BufLen(ip)       = 0
    END DO

    ! um_ak_20120413+
    IF (ltwoway) THEN 
       Parent(Id)%RMId =ChildId ! um_ak_20120413
       Parent(Id)%inter_npes = Child(Id)%inter_npes 
    
       ! allocate structure for Child PE information 
       ALLOCATE(Parent(Id)%PEs(0:Parent(Id)%inter_npes-1))
       ALLOCATE(ParBL(Id)%BufLen(0:Parent(Id)%inter_npes-1))

       DO ip =0, Parent(Id)%inter_npes-1
          Parent(Id)%PEs(ip)%NrEle = 0
          NULLIFY(Parent(Id)%PEs(ip)%locInd)
          ParBL(Id)%BufLen(ip) = 0
       END DO
    ENDIF
    ! um_ak_20120413-

    RETURN

   END SUBROUTINE MMD_P_Init
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_Set_Indexlist (Id, index_list)

    USE MMD_utilities,           ONLY: sort
    USE MMD_handle_communicator, ONLY: m_model_comm, m_model_rank, m_model_npes

    IMPLICIT NONE

    EXTERNAL  :: MPI_RECV, MPI_SEND
    INTRINSIC :: ASSOCIATED, SIZE

    ! Index list will be sorted
    INTEGER               , INTENT(IN)    :: Id
    INTEGER,DIMENSION(:,:), INTENT(INOUT) :: index_list  

    ! LOCAL
    INTEGER  :: i,ip,is,ie,ian, ind
    INTEGER  :: status
    INTEGER  :: RemPe
    INTEGER, ALLOCATABLE, DIMENSION(:) :: RemInd
    INTEGER, DIMENSION(1) :: NrEle
    INTEGER  :: num_indpairs

    TYPE indPair
       INTEGER,    DIMENSION(:,:), POINTER :: ij => NULL()
    END TYPE indPair

    TYPE(indPair), DIMENSION(:),   POINTER :: clnt_ind => NULL()


    ! 1. SPLIT INDEX LIST INTO PARENT SPECIFIC LISTS
    ! the PE with model_rank=0 evaluates the index_list and distributes
    ! the parent Pe specific list to each parent PE
    !-----------------------------------------------------------------------
    IF (m_model_rank == 0)   THEN
       CALL sort (index_list, 6)         ! Sort to ascending Parent PE
       is = 1
       DO ip=0,m_model_npes-1
          ! Split into Parent PEs
          ie = is-1                      ! there may my no entry for This PE

          IF (is <= SIZE(index_list,2).AND. is >=0 )  THEN
             DO WHILE ( index_list(6,ie+1) == ip)
                ie = ie+1
                IF ( ie == SIZE(index_list,2)) EXIT
             END DO
             ian = ie-is+1
          ELSE
             is  = -1
             ie  = -2
             ian = 0
          END IF
          
          ! Send data to other parent PEs
          IF (ip == 0)   THEN
             Child(Id)%NrPoints = ian
             IF (ian > 0)   THEN
                ALLOCATE (Child(Id)%index_list_2d(6,ian))
                Child(Id)%index_list_2d(:,1:ian) = index_list(:,is:ie)
             END IF
          ELSE
             CALL MPI_Send (ian, 1, MPI_INTEGER, ip, 1000, m_model_comm,status)
             IF (ian > 0) THEN
                CALL MPI_Send (index_list(1,is), 6*ian, MPI_INTEGER, ip, 1001,&
                     m_model_comm, status)
             END IF
          END IF
          is = ie+1
       END DO
    ELSE
       CALL MPI_Recv (Child(Id)%NrPoints, 1, MPI_INTEGER, 0, 1000 &
            , m_model_comm,   MPI_STATUS_IGNORE, status)
       ian = Child(Id)%NrPoints

       IF (ian > 0) THEN
          ALLOCATE(Child(Id)%index_list_2d(6,ian))
          CALL MPI_RECV (Child(Id)%index_list_2d, 6*ian &
               , MPI_INTEGER, 0, 1001, m_model_comm, MPI_STATUS_IGNORE, status)
       END IF
    END IF

    !-----------------------------------------------------------------------
    ! 2. SORT LIST ALONG Child PE numbers
    !-----------------------------------------------------------------------
    ! a) SET NUMBER OF ELEMENTS ON EACH Child PE TO ZERO
    DO ip=0,Child(Id)%inter_npes-1
       Child(Id)%PEs(ip)%NrEle = 0
    END DO

    ! b) COUNT NUMBER OF ELEMENTS ON EACH Child PE
    DO i=1,Child(Id)%NrPoints
       RemPe = Child(Id)%index_list_2d(5,i)
       Child(Id)%PEs(RemPe)%NrEle = &
            Child(Id)%PEs(RemPe)%NrEle+1
    END DO

    ! c) EVALUATE REMOTE INDEX PAIRs  
    ALLOCATE(RemInd(0:Child(Id)%inter_npes-1))
    ! INITIALIZE remote index
    DO ip=0,Child(Id)%inter_npes-1
       ALLOCATE(Child(Id)%PEs(ip)%locInd(Child(Id)%PEs(ip)%NrEle))
       RemInd(ip) = 0
    END DO

    ! ALLOCATE ARRAY FOR Child i,j index Pairs
    ALLOCATE(clnt_ind(0:Child(Id)%inter_npes-1))
    DO ip=0,Child(Id)%inter_npes-1
       ALLOCATE(clnt_ind(ip)%ij(Child(Id)%PEs(ip)%NrEle,2))
    ENDDO
    ! FIND PARENT and Child i,j Index Pairs per REMOTE PE (RemPe)
    DO i=1,Child(Id)%NrPoints
       RemPe = Child(Id)%index_list_2d(5,i)
       RemInd(RemPe) = RemInd(RemPe) +1 
       ind   = RemInd(RemPe)
       ! PARENT INDICES
       Child(Id)%PEs(RemPe)%locInd(ind)%i = &
            Child(Id)%index_list_2d(1,i)
       Child(Id)%PEs(RemPe)%locInd(ind)%j = &
            Child(Id)%index_list_2d(2,i)
       ! Child INDICES
       clnt_ind(RemPe)%ij(ind,1) = &
            Child(Id)%index_list_2d(3,i) 
       clnt_ind(RemPe)%ij(ind,2) = &
            Child(Id)%index_list_2d(4,i)
    END DO

    ! SEND NUMBER OF ELEMENTS AND Child i,j-INDEX PAIRS 
    ! WHICH ARE SEND FROM THIS PARENT PE TO SPECIFIC Child PE
    DO ip=0,Child(Id)%inter_npes-1
       NrEle(1) = Child(Id)%PEs(ip)%NrEle
       CALL MMD_Send_to_Child( Child(id)%RMId, NrEle, 1 &
                              , ip, 350+ip*2, status)
       num_indpairs = Child(Id)%PEs(ip)%NrEle*2
       CALL MMD_Send_to_Child( Child(id)%RMId, clnt_ind(ip)%ij &
                              , num_indpairs, ip, 351+ip*2, status)
       DEALLOCATE(clnt_ind(ip)%ij)
    END DO
    DEALLOCATE(clnt_ind)

    IF (ASSOCIATED(Child(Id)%index_list_2d))  &
         DEALLOCATE(Child(Id)%index_list_2d)
    DEALLOCATE(RemInd)

    RETURN

  END SUBROUTINE MMD_P_Set_Indexlist
  !-------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_P_Get_ParIndexlist(Id, fractions, wfunc)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    INTEGER,                  INTENT(IN) :: Id
    REAL(dp), DIMENSION(:,:), POINTER    :: fractions
    REAL(dp), DIMENSION(:,:), POINTER    :: wfunc 

    ! LOCAL
    INTEGER                              :: ip, mr, i
    INTEGER                              :: status

    INTEGER, DIMENSION(2)                :: NrEle
    ! NrEle(2) indicates if weighting function wfunc is sent
    !  NrEle(2) == 1 => wfunc is sent
    INTEGER                              :: idx, jdx

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: indfield
    REAL(dp),  DIMENSION(:), ALLOCATABLE :: inddpfield

    ! mr is the equivalent to ip in mmd_client. use here to create same tag     
    mr        = m_model_rank
    fractions = 0._dp
    wfunc     = 0._dp 

    ! GET NUMBER OF ELEMENTS PROVIDED BY EACH PARENT PE
    DO ip=0,Parent(Id)%inter_npes-1
       CALL MMD_Recv_from_Child(Parent(Id)%RMId, NrEle(1:2), 2&
            , ip, 1350+mr*5, status)

       Parent(Id)%PEs(ip)%NrEle = NrEle(1)
       ALLOCATE(Parent(Id)%PEs(ip)%locInd(Parent(Id)%PEs(ip)%NrEle))
       ALLOCATE(indfield(Parent(Id)%PEs(ip)%NrEle,2))
       ALLOCATE(inddpfield(Parent(Id)%PEs(ip)%NrEle))
       CALL MMD_Recv_from_Child(Parent(Id)%RMId, indfield &
            , SIZE(indfield), ip, 1351+mr*5, status) 
       CALL MMD_Recv_from_Child(Parent(Id)%RMId, inddpfield &
            , SIZE(inddpfield), ip, 1352+mr*5, status)
       DO i=1,Parent(Id)%PEs(ip)%NrEle
          Parent(Id)%PEs(ip)%locInd(i)%i      = indfield(i,1)
          Parent(Id)%PEs(ip)%locInd(i)%j      = indfield(i,2)
          Parent(Id)%PEs(ip)%locInd(i)%frac   = inddpfield(i)
          fractions(indfield(i,1),indfield(i,2)) = &
               fractions(indfield(i,1),indfield(i,2)) +  inddpfield(i)
       END DO

       DEALLOCATE(indfield)
       DEALLOCATE(inddpfield)
    ENDDO

    IF (NrEle(2) == 1) THEN
       DO ip=0,Parent(Id)%inter_npes-1
          ALLOCATE(inddpfield(Parent(Id)%PEs(ip)%NrEle))
          CALL MMD_Recv_from_Child(Parent(Id)%RMId, inddpfield &
               , SIZE(inddpfield), ip, 1353+mr*5, status)
          DO i=1,Parent(Id)%PEs(ip)%NrEle
             idx = Parent(Id)%PEs(ip)%locInd(i)%i
             jdx = Parent(Id)%PEs(ip)%locInd(i)%j 
             wfunc(idx,jdx)=  wfunc(idx,jdx) &
                  +  Parent(Id)%PEs(ip)%locInd(i)%frac * inddpfield(i)
          END DO
          DEALLOCATE(inddpfield)
       END DO
    END IF

    RETURN

  END SUBROUTINE MMD_P_Get_ParIndexlist
  ! ------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_Set_ParDataArray_Name(Id &
                                    , par_channel,par_object, chld_channel  &
                                    , chld_object, chld_repr, interpol_method &
                                    , sentunit, istat)

    USE  MMD_utilities,           ONLY: MMD_DA_NAME_ERR               &    
                                      , STRLEN_CHANNEL, STRLEN_OBJECT
    USE  MMD_handle_communicator, ONLY: m_to_child_comm

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, LEN, NULL, TRIM

    INTEGER         ,INTENT(IN)   :: Id
    CHARACTER(LEN=*),INTENT(IN)   :: par_object
    CHARACTER(LEN=*),INTENT(IN)   :: par_channel
    CHARACTER(LEN=*),INTENT(IN)   :: chld_object
    CHARACTER(LEN=*),INTENT(IN)   :: chld_channel
    CHARACTER(LEN=*),INTENT(IN)   :: chld_repr
    INTEGER,         INTENT(IN)   :: interpol_method
    LOGICAL,         INTENT(IN)   :: sentunit ! um_ak_20150413
    INTEGER,         INTENT(OUT)  :: istat

    ! LOCAL
    INTEGER, SAVE                 :: ParentId = 0
    INTEGER, SAVE                 :: myIndex  = 0
    INTEGER                       :: myPe
    TYPE(ArrayDef_list), POINTER  :: ai => NULL()
    TYPE(ArrayDef_list), POINTER  :: an => NULL()
    CHARACTER(LEN=STRLEN_OBJECT)  :: s_object
    CHARACTER(LEN=STRLEN_CHANNEL) :: s_channel
    CHARACTER(LEN=STRLEN_OBJECT)  :: c_repr
    INTEGER                       :: interpM
    INTEGER                       :: isentunit ! um_ak_20150413

    IF (ParentId /= Parent(Id)%RMId) THEN
       ParentId = Parent(Id)%RMId
       myIndex  = 0
    ENDIF

    istat = MMD_STATUS_OK
    ! Name too long
    IF (LEN(trim(par_object)) > STRLEN_OBJECT .OR.            &
       LEN(trim(chld_object)) > STRLEN_OBJECT )  THEN 
      istat = MMD_DA_NAME_ERR
    END IF

    myIndex = myIndex+1

    ! Broadcat to all Parent PEs
    IF (m_model_rank == 0) THEN
      myPE = MPI_ROOT
      s_channel = TRIM(chld_channel)
      s_object  = TRIM(chld_object)
      c_repr    = TRIM(chld_repr)
      interpM   = interpol_method
      IF (sentunit) THEN
         isentunit = 1
      ELSE
         isentunit = 0
      ENDIF
    else
      myPE = MPI_PROC_NULL
    ENDif

    CALL MMD_Bcast ( myIndex,   myPE, comm=m_to_child_comm(ParentId))
    CALL MMD_Bcast ( s_channel, myPE, comm=m_to_child_comm(ParentId))
    CALL MMD_Bcast ( s_object,  myPE, comm=m_to_child_comm(ParentId))
    CALL MMD_Bcast ( c_repr,    myPE, comm=m_to_child_comm(ParentId))
    CALL MMD_Bcast ( interpM,   myPE, comm=m_to_child_comm(ParentId))
    CALL MMD_Bcast ( isentunit, myPE, comm=m_to_child_comm(ParentId))

    ! BUILD LIST DATA ARRAY LIST
    ai => Parent(Id)%ArrayStart
    DO 
       IF (.NOT. ASSOCIATED(ai)) EXIT
       an => ai
       ai => ai%next
    ENDDO
    
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(Parent(id)%ArrayStart)) THEN
       Parent(Id)%ArrayStart => ai  ! SET POINTER TO FIRST OBJECT
    ELSE
       an%next => ai        ! SET NEXT POINTER OF LAST OBJECT
                            ! TO NEW OBJECT
    ENDIF
    
    ai%ArrDef%channel     = TRIM(par_channel)
    ai%ArrDef%object      = TRIM(par_object)
    
    RETURN

  END SUBROUTINE MMD_P_Set_ParDataArray_Name
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_P_Set_ParDataArray_EndList(Id)

    USE  MMD_handle_communicator, ONLY: m_to_child_comm

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Id

    ! LOCAL
    INTEGER                               :: myPe
    INTEGER                               :: myIndex
    INTEGER                               :: ParentId

    ParentId = Parent(Id)%RMId

    ! THIS Routine breaks the list of data arrays
    myIndex  = -1

    IF (m_model_rank == 0) THEN
      myPE = MPI_ROOT
    ELSE
      myPE = MPI_PROC_NULL
    ENDIF
    CALL MMD_Bcast ( myIndex,  myPE, comm=m_to_child_comm(ParentId))

    RETURN

  END SUBROUTINE MMD_P_Set_ParDataArray_EndList
  ! ------------------------------------------------------------------------
  ! um_ak_20120413-

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_Get_DataArray_Name(Id)

    USE mmd_utilities,           ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
    USE MMD_handle_communicator, ONLY: m_to_child_comm

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    INTEGER,INTENT(IN)            :: Id

    ! LOCAL 
    INTEGER                       :: ChildId
    INTEGER                       :: myIndex
    TYPE(ArrayDef_list), POINTER  :: ai => NULL()
    TYPE(ArrayDef_list), POINTER  :: ap => NULL()
    CHARACTER(LEN=STRLEN_OBJECT)  :: par_object
    CHARACTER(LEN=STRLEN_CHANNEL) :: par_channel
    CHARACTER(LEN=STRLEN_OBJECT)  :: chld_repr
    INTEGER                       :: isentunit ! ub_ak_20190614
 
    ChildId  = Child(Id)%RMId

    ! Get Data Array Channel and Name from Child
    DO 
       CALL MMD_Bcast ( myIndex,  0, comm=m_to_child_comm(ChildId))
       IF (myIndex == -1) EXIT

       CALL MMD_Bcast (par_channel, 0, comm=m_to_child_comm(ChildId))
       CALL MMD_Bcast (par_object,  0, comm=m_to_Child_comm(ChildId))
       CALL MMD_Bcast (chld_repr,   0, comm=m_to_Child_comm(ChildId))
       CALL MMD_Bcast (isentunit,   0, comm=m_to_Child_comm(ChildId)) ! ub_ak_20190614

       ! BUILD LIST DATA ARRAY LIST
       ai => Child(Id)%ArrayStart
       DO 
          IF (.NOT. ASSOCIATED(ai)) EXIT
          ap => ai
          ai => ai%next
       ENDDO
          
       ALLOCATE(ai)
       NULLIFY(ai%next)
       IF (.NOT. ASSOCIATED(Child(Id)%ArrayStart)) THEN
          ! SET POINTER TO FIRST OBJECT
          Child(Id)%ArrayStart => ai  
       ELSE
          ap%next => ai    ! SET NEXT POINTER OF LAST OBJECT TO NEW OBJECT
       ENDIF
       
       ai%ArrDef%channel     = par_channel
       ai%ArrDef%object      = par_object 
       ai%ArrDef%repr        = chld_repr
       IF (isentunit == 1) THEN
          ai%ArrDef%l_sentunit  = .TRUE.
       ELSE
          ai%ArrDef%l_sentunit  = .FALSE.
       END IF

    ENDDO
    
    RETURN

   END SUBROUTINE MMD_P_Get_DataArray_Name
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  LOGICAL FUNCTION MMD_P_GetNextArray (Id, MyChannel, myName, repr, l_unit)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM
  
    INTEGER,INTENT(IN)           :: Id
    CHARACTER(LEN=*),INTENT(OUT) :: myChannel
    CHARACTER(LEN=*),INTENT(OUT) :: myName
    CHARACTER(LEN=*),INTENT(OUT) :: repr
    LOGICAL         ,INTENT(OUT) :: l_unit

    ! LOCAL
    LOGICAL, SAVE                :: lfirst =.TRUE.

    IF (lfirst) THEN
       Child(Id)%Ar => Child(Id)%ArrayStart
    ELSE
       Child(Id)%Ar => Child(Id)%Ar%next
       IF (.NOT. ASSOCIATED(Child(Id)%Ar)) THEN
          MMD_P_GetNextArray = .FALSE.
          Child(Id)%Ar => Child(Id)%ArrayStart
          lfirst = .TRUE.
          RETURN
       ENDIF
    ENDIF
    
    MyChannel = TRIM(Child(Id)%Ar%ArrDef%channel)
    myName    = TRIM(Child(Id)%Ar%ArrDef%object)
    repr      = TRIM(Child(Id)%Ar%ArrDef%repr)
    l_unit    = Child(Id)%Ar%ArrDef%l_SentUnit

    lfirst = .FALSE.
    MMD_P_GetNextArray = .TRUE.

    RETURN 

  END FUNCTION MMD_P_GetNextArray
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_Send_Repr(axis, gdimlen, name, att, ChildId)

    USE MMD_utilities,           ONLY: STRLEN_MEDIUM, STRLEN_ATT
    USE MMD_handle_communicator, ONLY: m_model_rank, m_to_Child_comm

    IMPLICIT NONE
    
    CHARACTER(LEN=4),              INTENT(INOUT) :: axis
    INTEGER, DIMENSION(4),         INTENT(INOUT) :: gdimlen
    CHARACTER(LEN=STRLEN_MEDIUM),  INTENT(INOUT) :: name
    CHARACTER(LEN=STRLEN_ATT),     INTENT(INOUT) :: att
    INTEGER,                       INTENT(IN)    :: ChildId
    ! LOCAL
    INTEGER                                      :: myPE
    INTEGER                                      :: i

    IF (m_model_rank == 0) THEN
      myPE = MPI_ROOT
    else
      myPE = MPI_PROC_NULL
    ENDif
    CALL MMD_Bcast ( axis, myPE, comm=m_to_Child_comm(ChildId))
    CALL MMD_Bcast ( name, myPE, comm=m_to_Child_comm(ChildId))
    CALL MMD_Bcast ( att,  myPE, comm=m_to_Child_comm(ChildId))
    DO i=1,4
       CALL MMD_Bcast ( gdimlen(i), myPE, comm=m_to_Child_comm(ChildId))
    ENDDO

  END SUBROUTINE MMD_P_Send_Repr
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_Set_DataArray(Id, status, DimLen, ArrayOrder, p4)

    IMPLICIT NONE

    INTEGER, INTENT(IN)                    :: Id
    INTEGER, INTENT(OUT)                   :: status
    INTEGER, INTENT(IN),DIMENSION(4)       :: DimLen
    CHARACTER(LEN=4), INTENT(IN)           :: ArrayOrder
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: p4

    ! LOCAL
    INTEGER                       :: ip, idim
    INTEGER                       :: ArrFac

    Child(Id)%Ar%Arrdef%dim = DimLen
    Child(Id)%Ar%Arrdef%dim_order   = ArrayOrder
    Child(Id)%Ar%ArrDef%xyzn_dim(:) = -1
    DO idim = 1,4
       IF (ArrayOrder(idim:idim) == 'X') THEN
          Child(Id)%Ar%ArrDef%xyzn_dim(1) = idim
       ELSE IF (ArrayOrder(idim:idim) == 'Y') THEN
          Child(Id)%Ar%ArrDef%xyzn_dim(2) = idim
       ELSE IF (ArrayOrder(idim:idim) == 'Z') THEN
          Child(Id)%Ar%ArrDef%xyzn_dim(3) = idim
       ELSE IF (ArrayOrder(idim:idim) == 'N') THEN
          Child(Id)%Ar%ArrDef%xyzn_dim(4) = idim
       ENDIF
    END DO
    ArrFac = 1
    IF (Child(Id)%Ar%ArrDef%xyzn_dim(3) > 0) THEN
       ArrFac = ArrFac *DimLen(Child(Id)%Ar%ArrDef%xyzn_dim(3))
    ENDIF
    IF (Child(Id)%Ar%ArrDef%xyzn_dim(4) > 0) THEN
       ArrFac = ArrFac *DimLen(Child(Id)%Ar%ArrDef%xyzn_dim(4))
    ENDIF

    ! ALLOCATE ARRAY FOR DIMENSION OF EACH ARRAY TRANSFERED FROM PARENT
    ! TO CHILD
    ALLOCATE(Child(Id)%Ar%Arrdef%ArrLen(0:Child(Id)%inter_npes))
    DO ip = 0,Child(Id)%inter_npes-1
       Child(Id)%Ar%Arrdef%ArrLen(ip) =  &
            ArrFac * Child(Id)%PEs(ip)%NrEle
       ChldBL(Id)%BufLen(ip) = ChldBL(Id)%BufLen(ip) &
            + Child(Id)%Ar%Arrdef%ArrLen(ip)
    END DO

    Child(Id)%Ar%Arrdef%p4 => p4
    
    status = 0

    RETURN

  END SUBROUTINE MMD_P_Set_DataArray
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  LOGICAL FUNCTION MMD_P_GetNextParArray (Id, MyChannel, myName)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM
  
    INTEGER,INTENT(IN)           :: Id
    CHARACTER(LEN=*),INTENT(OUT) :: myChannel
    CHARACTER(LEN=*),INTENT(OUT) :: myName

    ! LOCAL
    LOGICAL, SAVE                :: lfirst =.TRUE.

    IF (lfirst) THEN
       Parent(Id)%Ar => Parent(Id)%ArrayStart
    ELSE
       Parent(Id)%Ar => Parent(Id)%Ar%next
       IF (.NOT. ASSOCIATED(Parent(Id)%Ar)) THEN
          MMD_P_GetNextParArray = .FALSE.
          Parent(Id)%Ar => Parent(Id)%ArrayStart
          lfirst = .TRUE.
          RETURN
       ENDIF
    ENDIF
    
    MyChannel = TRIM(Parent(Id)%Ar%ArrDef%channel)
    myName    = TRIM(Parent(Id)%Ar%ArrDef%object)

    lfirst = .FALSE.
    MMD_P_GetNextParArray = .TRUE.

    RETURN 

  END FUNCTION MMD_P_GetNextParArray
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_Set_ParDataArray (Id, status, DIMLEN, ArrayOrder, p4)

    IMPLICIT NONE
    
    INTRINSIC :: ASSOCIATED

    INTEGER, INTENT(IN)                    :: Id
    INTEGER, INTENT(OUT)                   :: status
    INTEGER, INTENT(IN),DIMENSION(4)       :: DIMLEN
    CHARACTER(LEN=4), INTENT(IN)           :: ArrayOrder
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: p4

    ! LOCAL
    INTEGER                      :: ip, idim
    TYPE(ArrayDef_list), POINTER :: Ar => NULL()
    INTEGER                      :: ArrFac

    Ar => Parent(Id)%Ar
    Ar%ArrDef%dim(:) = 0
    
    DO idim = 1,4
       Ar%ArrDef%dim(idim) = DIMLEN(idim)
    ENDDO
    
    Ar%ArrDef%dim_order = ArrayOrder

    Ar%ArrDef%xyzn_dim(:) = -1
    DO idim = 1,4
       IF (ArrayOrder(idim:idim) == 'X') THEN
          Ar%ArrDef%xyzn_dim(1)    = idim
       ELSE IF (ArrayOrder(idim:idim) == 'Y') THEN
          Ar%ArrDef%xyzn_dim(2) = idim
       ELSE IF (ArrayOrder(idim:idim) == 'Z') THEN
          Ar%ArrDef%xyzn_dim(3) = idim
       ELSE IF (ArrayOrder(idim:idim) == 'N') THEN
          Ar%ArrDef%xyzn_dim(4) = idim
       ENDIF
    END DO

    ! CALCULATE FACTOR FOR ADDITIONAL DIMENSIONS
    ArrFac = 1
    IF (Parent(Id)%Ar%ArrDef%xyzn_dim(3) > 0) THEN
       ArrFac = ArrFac * DIMLEN(Parent(Id)%Ar%ArrDef%xyzn_dim(3))
    END IF
    IF (Parent(Id)%Ar%ArrDef%xyzn_dim(4) > 0) THEN
       ArrFac = ArrFac * DIMLEN(Parent(Id)%Ar%ArrDef%xyzn_dim(4))
    END IF

    ! CALCULATE ArrayLength FOR EACH REMOTE PE
    ALLOCATE(Parent(Id)%Ar%Arrdef%ArrLen(0:Parent(Id)%inter_npes))
    ALLOCATE(Parent(Id)%Ar%Arrdef%ArrIdx(0:Parent(Id)%inter_npes))
    DO ip=0,Parent(Id)%inter_npes-1
       Parent(Id)%Ar%Arrdef%ArrLen(ip) = ArrFac * Parent(Id)%PEs(ip)%NrEle
       Parent(Id)%Ar%Arrdef%ArrIdx(ip) = ParBL(Id)%BufLen(ip) + 1
       ParBL(Id)%BufLen(ip) = ParBL(Id)%BufLen(ip) + Parent(Id)%Ar%Arrdef%ArrLen(ip)
    END DO

    IF (ASSOCIATED(p4)) THEN
       Ar%ArrDef%p4 => p4
    ELSE
       status = 100
       RETURN
    ENDIF

    status = 0

   RETURN

 END SUBROUTINE MMD_P_Set_ParDataArray
  ! ------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_SetInd_and_AllocMem (Id)

    IMPLICIT NONE

    EXTERNAL  :: MMDc_P_SetInd_and_Mem, MMDc_P_SetInd_and_Mem_2way

    INTEGER, INTENT(IN) :: Id

    ! LOCAL

    IF (.NOT. ltwoway) THEN
       CALL MMDc_P_SetInd_and_Mem (Child(Id)%RMId, ChldBL(Id)%BufLen)
    ELSE
       CALL MMDc_P_SetInd_and_Mem_2way ( Child(Id)%RMId &
                                       , ChldBL(Id)%BufLen ,ParBL(Id)%BufLen)
    END IF

    RETURN

  END SUBROUTINE MMD_P_SetInd_and_AllocMem
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_P_FillBuffer(Id, WaitTime)

    USE MMD_utilities,           ONLY: PeDef

    IMPLICIT NONE

    EXTERNAL  :: MMDC_P_FILLBUFFER, MMDC_P_SETBARRIER
    INTRINSIC :: ASSOCIATED, PRESENT
    
    INTEGER, INTENT(IN)            :: Id
    REAL(dp),INTENT(OUT), OPTIONAL :: WaitTime

    ! LOCAL
    TYPE(ArrayDef_list), POINTER   :: Ar  => NULL()
    TYPE(PeDef),         POINTER   :: aPE => NULL()
    INTEGER, DIMENSION(4)          :: aind  ! Index in Pointer
    INTEGER                        :: ij, ip, idim, k, n, ind
    INTEGER                        :: zloopend, nloopend

    REAL(DP), DIMENSION(:), POINTER :: array => NULL()
    REAL(DP)                        :: my_waittime
    LOGICAL, SAVE                   :: lfirst  = .TRUE.

    ! This not only get the WaitTime but also includes on barrier Call
    ! for the 2-way coupling the barrier call is required in the getbuffer
    ! routine
    IF (.NOT. ltwoway .OR. lfirst) THEN
       my_WaitTime = -1.0
       IF (.NOT. ltwoway) &
            CALL MMDc_P_GetWaitTime(Child(Id)%RMId, my_WaitTime)
       IF (PRESENT(WaitTime))  then
          WaitTime = my_WaitTime
       end if
       lfirst = .FALSE.
    ENDIF

    PE_loop: DO ip=0,Child(Id)%inter_npes-1
       Child(Id)%Ar => Child(Id)%ArrayStart
       Ar => Child(Id)%Ar
       aPE =>  Child(Id)%PEs(ip)
       ALLOCATE(array(ChldBL(Id)%BufLen(ip)))
       ind = 0
       
       arrayloop: DO
          !IF (TRIM(Ar%Arrdef%object) == 'ALLTRACER') &
          !   write (0,*) 'PARENT: ',m_model_rank, ip, TRIM(Ar%Arrdef%object) &
          !        , Ar%Arrdef%dim_order, Ar%Arrdef%ArrLen(ip) &
          !        , Ar%Arrdef%xyzn_dim, Ar%Arrdef%dim &
          !        , 'BD ',LBOUND(Ar%Arrdef%p4), UBOUND(Ar%Arrdef%p4)
          IF (Ar%ArrDef%dim_order == 'XY--') THEN
             DO ij = 1, aPE%NrEle
                ind = ind + 1  
                aind(1) = aPE%locInd(ij)%i
                aind(2) = aPE%locInd(ij)%j
                aind(3) = 1
                aind(4) = 1
                array(ind)= Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))
             END DO
         ELSE IF (Ar%ArrDef%dim_order == 'X-Y-') THEN
             DO ij = 1, aPE%NrEle
                ind = ind + 1  
                aind(1) = aPE%locInd(ij)%i
                aind(2) = 1
                aind(3) = aPE%locInd(ij)%j
                aind(4) = 1
                array(ind)= Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))
             END DO
          ELSE IF(Ar%ArrDef%dim_order == 'XYZ-') THEN
             DO ij = 1, aPE%NrEle
                DO k=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(3))
                   ind = ind + 1  
                   aind(1) = aPE%locInd(ij)%i
                   aind(2) = aPE%locInd(ij)%j
                   aind(3) = k
                   aind(4) = 1
                   array(ind)= Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))
                END DO
             END DO
          ELSE IF(Ar%ArrDef%dim_order == 'XZY-') THEN
             DO ij = 1, aPE%NrEle
                DO k=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(3))
                   ind = ind + 1  
                   aind(1) = aPE%locInd(ij)%i
                   aind(2) = k
                   aind(3) = aPE%locInd(ij)%j
                   aind(4) = 1
                   array(ind)= Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))
                END DO
             END DO
          ELSE IF(Ar%ArrDef%dim_order == 'XZNY') THEN
             DO ij = 1, aPE%NrEle
                DO k=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(3))
                   DO n=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(4))
                      ind = ind + 1  
                      aind(1) = aPE%locInd(ij)%i
                      aind(2) = k
                      aind(3) = n
                      aind(4) = aPE%locInd(ij)%j
                      array(ind)= Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))
                   END DO
                END DO
             END DO
          ELSE IF(Ar%ArrDef%dim_order == 'XYNZ') THEN
             DO ij = 1, aPE%NrEle
                DO k=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(3))
                   DO n=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(4))
                      ind = ind + 1  
                      aind(1) = aPE%locInd(ij)%i
                      aind(2) = aPE%locInd(ij)%j
                      aind(3) = n
                      aind(4) = k
                      array(ind)= Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))
                   END DO
                END DO
             END DO
          ELSE IF(Ar%ArrDef%dim_order == 'XYZN') THEN
             DO ij = 1, aPE%NrEle
                DO k=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(3))
                   DO n=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(4))
                      ind = ind + 1  
                      aind(1) = aPE%locInd(ij)%i
                      aind(2) = aPE%locInd(ij)%j
                      aind(3) = k
                      aind(4) = n
                      array(ind)= Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))
                   END DO
                END DO
             END DO
          ELSE
             zloopend = 1
             nloopend = 1
             IF (Ar%ArrDef%xyzn_dim(3)>0) &
                  zloopend=Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(3))
             IF (Ar%ArrDef%xyzn_dim(4)>0) &
                  nloopend=Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(4))

             DO ij = 1, aPE%NrEle
                DO k=1, zloopend  ! Z dimension
                   DO n =1, nloopend ! N dimension
                      ind = ind + 1 
                      aind(:) = 1
                      DO idim=1,4
                         IF ( Ar%Arrdef%dim_order(idim:idim)=='X') THEN
                            aind(idim) = aPE%locInd(ij)%i
                         ELSE IF (Ar%Arrdef%dim_order(idim:idim)=='Y') THEN
                            aind(idim) = aPE%locInd(ij)%j
                         ELSE IF (Ar%Arrdef%dim_order(idim:idim)=='Z') THEN
                            aind(idim) = k
                         ELSE IF (Ar%Arrdef%dim_order(idim:idim)=='N') THEN
                            aind(idim) = n
                         ENDIF
                      END DO
                      array(ind)= Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))
                   END DO
                END DO
             ENDDO
          ENDIF

          Ar => Ar%next
          IF (.NOT. ASSOCIATED(Ar)) EXIT
       END DO arrayloop

       CALL MMDc_P_FillBuffer(Child(Id)%RMId, ip, array)

       DEALLOCATE(array); NULLIFY(array)
       NULLIFY(Ar)
       NULLIFY(aPE)

    END DO PE_loop

    ! Set Barrier when Buffer Full
    CALL MMDc_P_SetBarrier(Child(Id)%RMId)

    RETURN

  END SUBROUTINE MMD_P_FillBuffer
  !-------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_P_GetBuffer(Id, Waittime)

    USE  MMD_utilities,            ONLY: PeDef, ArrayDef_list

    IMPLICIT NONE

    EXTERNAL  :: MMDC_P_GETBUFFER, MMDC_P_GETWAITTIME, MMDC_P_SETBARRIER
    INTRINSIC :: ASSOCIATED, PRESENT

    INTEGER, INTENT(IN)                    :: Id
     REAL(kind=dp),INTENT(INOUT), OPTIONAL :: WaitTime

    ! LOCAL
    TYPE(ArrayDef_list),    POINTER :: Ar  => NULL()
    TYPE(PeDef),            POINTER :: aPE => NULL()
    INTEGER, DIMENSION(4)           :: aind
    INTEGER                         :: ind, ij, ip, k,n, idim
    REAL(DP), DIMENSION(:), POINTER :: array => NULL()
    REAL(dp)                        :: my_waittime
    INTEGER                         :: zloopend, nloopend
    REAL(dp)                        :: frac

    ! This not only get the WaitTime but also includes on barrier Call
    CALL MMDc_P_GetWaitTime(Parent(Id)%RMId, my_WaitTime)
    IF (PRESENT(WaitTime)) WaitTime = my_WaitTime

    PE_loop: DO ip=0,Parent(Id)%inter_npes-1
       IF (.NOT. ASSOCIATED(Parent(Id)%ArrayStart)) EXIT ! um_ak_20150324

       Parent(Id)%Ar => Parent(Id)%ArrayStart
       Ar    => Parent(Id)%Ar
       
       ALLOCATE(array(ParBL(Id)%BufLen(ip)))
       array = 0._dp ! ub_ak_20170711
       aPE =>  Parent(Id)%PEs(ip)

       CALL MMDc_P_GetBuffer(Parent(Id)%RMId,ip, array)

       ind = 0

       arrayloop: DO
          !write (0,*) 'Child: ',ip, TRIM(Ar%Arrdef%object) &
          !     , Ar%Arrdef%dim_order &
          !     , Ar%Arrdef%xyzn_dim, Ar%Arrdef%dim  
          !IF (ASSOCIATED(Ar%Arrdef%p4)) write (0,*) 'Child p4' &
          !     , LBOUND(Ar%Arrdef%p4), UBOUND(Ar%Arrdef%p4)
          ind = Ar%Arrdef%ArrIdx(ip)
          zloopend=1
          nloopend=1
          
          IF (ip==0) Ar%Arrdef%p4(:,:,:,:) = 0._dp
          
          IF (Ar%Arrdef%ArrLen(ip) == 0) THEN
             ! DO NOTHING
          ELSE IF (Ar%Arrdef%dim_order == 'XY--') THEN
             DO ij = 1, aPE%NrEle
                aind(1) = aPE%locInd(ij)%i
                aind(2) = aPE%locInd(ij)%j
                frac    = aPE%locInd(ij)%frac
                Ar%Arrdef%p4(aind(1),aind(2),1,1) =  &
                     Ar%Arrdef%p4(aind(1),aind(2),1,1) + frac * array(ind)
                ind     = ind+1
             END DO
          ELSE IF (Ar%Arrdef%dim_order == 'X-Y-') THEN
             DO ij = 1, aPE%NrEle
                aind(1) = aPE%locInd(ij)%i
                !aind(2) = 1
                aind(3) = aPE%locInd(ij)%j
                frac    = aPE%locInd(ij)%frac
                Ar%Arrdef%p4(aind(1),1,aind(3),1) = &
                     Ar%Arrdef%p4(aind(1),1,aind(3),1) + frac * array(ind)
                ind     = ind+1
             END DO
          ELSE IF (Ar%Arrdef%dim_order == 'XYZ-') THEN
             DO ij = 1, aPE%NrEle
                DO k=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(3))
                   aind(1) = aPE%locInd(ij)%i
                   aind(2) = aPE%locInd(ij)%j
                   aind(3) = k
                   frac     = aPE%locInd(ij)%frac
                   Ar%Arrdef%p4(aind(1),aind(2),aind(3),1) =  &
                        Ar%Arrdef%p4(aind(1),aind(2),aind(3),1) &
                        + frac* array(ind)
                   ind     = ind+1
                END DO
             END DO
          ELSE IF(Ar%ArrDef%dim_order == 'XZY-') THEN
             DO ij = 1, aPE%NrEle
                DO k=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(3))
                   aind(1) = aPE%locInd(ij)%i
                   aind(2) = k
                   aind(3) = aPE%locInd(ij)%j
                   frac     = aPE%locInd(ij)%frac
                   Ar%Arrdef%p4(aind(1),aind(2),aind(3),1) = &
                        Ar%Arrdef%p4(aind(1),aind(2),aind(3),1) &
                        + frac* array(ind)
                   ind = ind + 1  
                END DO
             END DO
          ELSE IF (Ar%Arrdef%dim_order == 'XYNZ') THEN
             test: IF (TRIM(Ar%Arrdef%object) == 'Test_Ar') THEN
                DO ij = 1, aPE%NrEle
                   DO k=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(3))
                      DO n=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(4))
                         aind(1) = aPE%locInd(ij)%i
                         aind(2) = aPE%locInd(ij)%j
                         aind(3) = n
                         aind(4) = k
                         frac     = aPE%locInd(ij)%frac
                         Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) &
                              = array(ind) 
                         ind     = ind + 1
                      END DO
                   END DO
                END DO
             ELSE
                DO ij = 1, aPE%NrEle
                   DO k=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(3))
                      DO n=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(4))
                         aind(1) = aPE%locInd(ij)%i
                         aind(2) = aPE%locInd(ij)%j
                         aind(3) = n
                         aind(4) = k
                         frac     = aPE%locInd(ij)%frac
                         Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) =  &
                              Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) &
                              + frac * array(ind) 
                         ind     = ind + 1
                      END DO
                   END DO
                END DO
             END IF test
          ELSE IF (Ar%Arrdef%dim_order == 'XZNY') THEN
             test2: IF (TRIM(Ar%Arrdef%object) == 'Test_Ar') THEN
                DO ij = 1, aPE%NrEle
                   DO k=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(3))
                      DO n=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(4))
                         aind(1) = aPE%locInd(ij)%i
                         aind(2) = k
                         aind(3) = n
                         aind(4) = aPE%locInd(ij)%j
                         Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4))  &
                              = array(ind) 
                         ind     = ind + 1
                      END DO
                   END DO
                END DO

             ELSE

                DO ij = 1, aPE%NrEle
                   DO k=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(3))
                      DO n=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(4))
                         aind(1) = aPE%locInd(ij)%i
                         aind(2) = k
                         aind(3) = n
                         aind(4) = aPE%locInd(ij)%j
                         frac     = aPE%locInd(ij)%frac
                         Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) =  &
                              Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) &
                              + frac * array(ind) 
                         ind     = ind + 1
                      END DO
                   END DO
                END DO
             END IF test2
          ELSE
             IF (Ar%Arrdef%xyzn_dim(3) > 0)  &
                  zloopend=Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(3))
             IF (Ar%Arrdef%xyzn_dim(4) > 0)  &
                  nloopend=Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(4))
             DO ij = 1, aPE%NrEle
                DO k=1, zloopend   ! Z dimension
                   DO n =1, nloopend ! N dimension
                      aind(:) = 1
                      DO idim=1,4
                         IF ( Ar%Arrdef%dim_order(idim:idim)=='X') THEN
                            aind(idim) = aPE%locInd(ij)%i
                         ELSE IF (Ar%Arrdef%dim_order(idim:idim)=='Y') THEN
                            aind(idim) = aPE%locInd(ij)%j
                         ELSE IF (Ar%Arrdef%dim_order(idim:idim)=='Z') THEN
                            aind(idim) = k
                         ELSE IF (Ar%Arrdef%dim_order(idim:idim)=='N') THEN
                            aind(idim) = n
                         ENDIF
                      END DO
                      frac     = aPE%locInd(ij)%frac
                      Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) =  &
                           Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) &
                           + frac * array(ind) 

                      ind     = ind + 1 
                   END DO
                END DO
             ENDDO
          ENDIF

          Ar => Ar%next
          IF (.NOT. ASSOCIATED(Ar)) EXIT ! exit array loop
       END DO Arrayloop
       DEALLOCATE(array); NULLIFY(array)
       NULLIFY(Ar)
    END DO PE_loop

    ! um_ak_20151208 CALL MMDc_P_SetBarrier(Child(Id)%RMId)

    RETURN

  END SUBROUTINE MMD_P_GetBuffer
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_P_FreeMem(NumChildren)

    IMPLICIT NONE

    EXTERNAL  :: MMDc_P_FreeMem   ! ub_ak_20180202
    INTRINSIC :: ASSOCIATED, NULL

    INTEGER, INTENT(IN)          :: NumChildren

    ! LOCAL
    INTEGER                      :: ip
    INTEGER                      :: Id
    TYPE(ArrayDef_list), POINTER :: ai => NULL()
    TYPE(ArrayDef_list), POINTER :: ae => NULL()

    DO Id = 1, NumChildren
       CALL MMDc_P_FreeMem (Child(Id)%RMId) ! ub_ak_20180202

       DO ip=0,Child(Id)%inter_npes-1
          DEALLOCATE(Child(Id)%PEs(ip)%locInd)
       END DO
       DEALLOCATE(Child(Id)%PEs)
       DEALLOCATE(Child(Id)%Ar%Arrdef%ArrLen)
       
       ai => Child(Id)%ArrayStart

       DO 
          IF (.NOT. ASSOCIATED(ai)) EXIT
          ae => ai
          ai => ai%next
          
          DEALLOCATE(ae)
          NULLIFY(ae)
       END DO

    END DO 
    DEALLOCATE(Child)

    ! PARENT
    IF (ltwoway) THEN

       DO Id = 1, NumChildren
          DO ip=0,Parent(Id)%inter_npes-1
             IF (ASSOCIATED(Parent(Id)%PEs(ip)%locInd)) &
                  DEALLOCATE(Parent(Id)%PEs(ip)%locInd)
          END DO
          DEALLOCATE(Parent(Id)%PEs)
          IF (ASSOCIATED(Parent(Id)%Ar)) THEN
             DEALLOCATE(Parent(Id)%Ar%Arrdef%ArrLen)

             ai => Parent(Id)%ArrayStart

             DO
                IF (.NOT. ASSOCIATED(ai)) EXIT
                ae => ai
                ai => ai%next

                DEALLOCATE(ae)
                NULLIFY(ae)
             END DO
          END IF
       END DO
       DEALLOCATE(Parent)
       DEALLOCATE(ParBL)
    ENDIF 

  END SUBROUTINE MMD_P_FreeMem
  !-------------------------------------------------------------------------

END MODULE mmd_parent
