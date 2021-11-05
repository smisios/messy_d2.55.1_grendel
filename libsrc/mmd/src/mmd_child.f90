MODULE mmd_child

  USE  mpi
  USE  MMD_utilities,            ONLY: MMD_STATUS_OK , MMD_DP                 &
                                     , MMD_ParentIsECHAM, MMD_ParentIsCOSMO   &
                                     , ExchDataDef, ArrayDef_list   
  USE  MMD_handle_communicator,  ONLY: m_model_rank, m_to_parent_comm
  USE  MMD_MPI_wrapper,          ONLY: MMD_Send_to_Parent                     &
                                     , MMD_Recv_from_Parent                   &
                                     , MMD_Bcast, MMD_Inter_Bcast

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER, PARAMETER             :: dp = MMD_DP
  TYPE(ExchDataDef)              :: Me
  TYPE(ExchDataDef)              :: Parent
  INTEGER, DIMENSION(:), POINTER :: BufLen    => NULL()
  INTEGER, DIMENSION(:), POINTER :: ParBufLen => NULL()

  LOGICAL                        :: ltwoway = .FALSE.

  ! INTERFACE section
  INTERFACE MMD_C_Init
    MODULE PROCEDURE MMD_C_Init
  END INTERFACE MMD_C_Init

  INTERFACE MMD_C_Set_DataArray_Name
    MODULE PROCEDURE MMD_C_Set_DataArray_Name
  END INTERFACE MMD_C_Set_DataArray_Name

  INTERFACE MMD_C_Set_DataArray_EndList
    MODULE PROCEDURE MMD_C_Set_DataArray_EndList
  END INTERFACE MMD_C_Set_DataArray_EndList

  INTERFACE MMD_C_Get_Indexlist
   MODULE PROCEDURE MMD_C_Get_Indexlist
  END INTERFACE MMD_C_Get_Indexlist

  INTERFACE MMD_C_GetNextArray
    MODULE PROCEDURE MMD_C_GetNextArray
  END INTERFACE MMD_C_GetNextArray

  INTERFACE MMD_C_get_Repr
     MODULE PROCEDURE MMD_C_get_Repr
  END INTERFACE MMD_C_get_Repr

  INTERFACE MMD_C_Set_DataArray
    MODULE PROCEDURE MMD_C_Set_DataArray
  END INTERFACE MMD_C_Set_DataArray

  INTERFACE MMD_C_Set_ParIndexlist
    MODULE PROCEDURE MMD_C_Set_ParIndexlist
  END INTERFACE MMD_C_Set_ParIndexlist

  INTERFACE MMD_C_Get_ParDataArray_Name
    MODULE PROCEDURE MMD_C_Get_ParDataArray_Name
  END INTERFACE MMD_C_Get_ParDataArray_Name

  INTERFACE MMD_C_GetNextParArray
    MODULE PROCEDURE MMD_C_GetNextParArray
  END INTERFACE MMD_C_GetNextParArray

  INTERFACE MMD_C_Set_ParDataArray
    MODULE PROCEDURE MMD_C_Set_ParDataArray
  END INTERFACE MMD_C_Set_ParDataArray
  
  INTERFACE MMD_C_setInd_and_AllocMem
    MODULE PROCEDURE MMD_C_setInd_and_AllocMem
  END INTERFACE MMD_C_setInd_and_AllocMem

  INTERFACE MMD_C_GetBuffer
    MODULE PROCEDURE MMD_C_GetBuffer
  END INTERFACE MMD_C_GetBuffer

  INTERFACE MMD_C_FillBuffer
    MODULE PROCEDURE MMD_C_FillBuffer
  END INTERFACE MMD_C_FillBuffer

  INTERFACE MMD_C_GetParentType
    MODULE PROCEDURE MMD_C_GetParentType
  END INTERFACE MMD_C_GetParentType

  INTERFACE MMD_C_FreeMem
     MODULE PROCEDURE MMD_C_FreeMem
  END INTERFACE MMD_C_FreeMem

  ! Public section

  PUBLIC :: MMD_C_Init, MMD_C_Set_DataArray_Name, MMD_C_Get_Indexlist 
  PUBLIC :: MMD_C_Set_DataArray_EndList
  PUBLIC :: MMD_C_GetNextArray, MMD_C_Set_DataArray, MMD_C_get_Repr

  PUBLIC :: MMD_C_Get_ParDataArray_Name, MMD_C_Set_ParDataArray
  PUBLIC :: MMD_C_GetNextParArray, MMD_C_Set_ParIndexlist
  PUBLIC :: MMD_C_FillBuffer

  PUBLIC :: MMD_C_setInd_and_AllocMem, MMD_C_GetBuffer, MMD_C_GetParentType
  PUBLIC :: MMD_STATUS_OK, MMD_Send_to_Parent, MMD_Recv_from_Parent
  PUBLIC :: MMD_ParentIsECHAM, MMD_ParentIsCOSMO
  PUBLIC :: MMD_Inter_Bcast, MMD_C_FreeMem

 CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_C_Init(l2way)

    USE  MMD_handle_communicator,  ONLY: m_model_comm
  
    IMPLICIT NONE

    EXTERNAL  :: MMDc_C_init
    INTRINSIC :: PRESENT

    LOGICAL, INTENT(IN), OPTIONAL :: l2way

    ! LOCAL
    INTEGER :: ip

    IF (PRESENT(l2way)) ltwoway = l2way

    ! STORE COMMUNICATOR IN CHILD STRUCT
    CALL MMDc_C_init (m_model_comm, m_to_parent_comm, Me%inter_npes)

    NULLIFY(Me%Ar)               
    NULLIFY(Me%ArrayStart)

    ALLOCATE(Me%PEs(0:Me%inter_npes-1))
    ALLOCATE(BufLen(0:Me%inter_npes-1))

    DO ip=0,Me%inter_npes-1
       BufLen(ip)= 0
       Me%Pes(ip)%NrEle=0
       NULLIFY(Me%Pes(ip)%locInd)
    END DO


    IF (ltwoway) THEN
       NULLIFY(Parent%Ar)         
       NULLIFY(Parent%ArrayStart) 

       Parent%inter_npes = Me%inter_npes
       ALLOCATE(Parent%PEs(0:Parent%inter_npes-1))
       ALLOCATE(ParBufLen(0:Parent%inter_npes-1)) 
       DO ip=0,Parent%inter_npes-1
          ParBufLen(ip)= 0                    
          Parent%Pes(ip)%NrEle=0
          NULLIFY(Parent%Pes(ip)%locInd)
       END DO

    ENDIF

    RETURN

  END SUBROUTINE MMD_C_Init
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE MMD_C_Set_DataArray_Name(par_channel,par_object, chld_channel &
                                     , chld_object, chld_repr, L_sentunit  &
                                     , istat)

    USE  MMD_utilities,           ONLY: MMD_DA_NAME_ERR               &    
                                      , STRLEN_CHANNEL, STRLEN_OBJECT

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, LEN, NULL, TRIM

    CHARACTER(LEN=*),INTENT(IN)   :: par_object
    CHARACTER(LEN=*),INTENT(IN)   :: par_channel
    CHARACTER(LEN=*),INTENT(IN)   :: chld_object
    CHARACTER(LEN=*),INTENT(IN)   :: chld_channel
    CHARACTER(LEN=*),INTENT(IN)   :: chld_repr
    LOGICAL         ,INTENT(IN)   :: l_SentUnit ! ub_ak_20190614
    INTEGER,INTENT(OUT)           :: istat

    ! LOCAL
    INTEGER, SAVE                 :: myIndex = 0
    INTEGER                       :: myPe
    TYPE(ArrayDef_list), POINTER  :: ai => NULL()
    TYPE(ArrayDef_list), POINTER  :: an => NULL()
    CHARACTER(LEN=STRLEN_OBJECT)  :: s_object
    CHARACTER(LEN=STRLEN_CHANNEL) :: s_channel
    CHARACTER(LEN=STRLEN_OBJECT)  :: c_repr
    INTEGER                       :: isentunit ! um_ak_20190614
     

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
      s_channel = TRIM(par_channel)
      s_object  = TRIM(par_object)
      c_repr    = TRIM(chld_repr)
      ! ub_ak_20190614+
      IF (l_sentunit) THEN
         isentunit = 1
      ELSE
         isentunit = 0
      ENDIF
      ! ub_ak_20190614-
    else
      myPE = MPI_PROC_NULL
    ENDif

    CALL MMD_Bcast ( myIndex,   myPE, comm=m_to_parent_comm)
    CALL MMD_Bcast ( s_channel, myPE, comm=m_to_parent_comm)
    CALL MMD_Bcast ( s_object,  myPE, comm=m_to_parent_comm)
    CALL MMD_Bcast ( c_repr,    myPE, comm=m_to_parent_comm)
    CALL MMD_Bcast ( isentunit, myPE, comm=m_to_parent_comm)

    ! BUILD LIST DATA ARRAY LIST
    ai => Me%ArrayStart
    DO 
       IF (.NOT. ASSOCIATED(ai)) EXIT
       an => ai
       ai => ai%next
    ENDDO
    
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(Me%ArrayStart)) THEN
       Me%ArrayStart => ai  ! SET POINTER TO FIRST OBJECT
    ELSE
       an%next => ai        ! SET NEXT POINTER OF LAST OBJECT
                            ! TO NEW OBJECT
    ENDIF
    
    ai%ArrDef%channel     = TRIM(chld_channel)
    ai%ArrDef%object      = TRIM(chld_object)
    
    RETURN

  END SUBROUTINE MMD_C_Set_DataArray_Name
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_C_Set_DataArray_EndList

   IMPLICIT NONE

    ! LOCAL
    INTEGER                               :: myPe
    INTEGER                               :: myIndex

    ! THIS Routine breaks the list of data arrays
    myIndex  = -1

    IF (m_model_rank == 0) THEN
      myPE = MPI_ROOT
    ELSE
      myPE = MPI_PROC_NULL
    ENDIF
    CALL MMD_Bcast ( myIndex,  myPE, comm=m_to_parent_comm)

    RETURN

  END SUBROUTINE MMD_C_Set_DataArray_EndList
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_C_Set_ParIndexlist (index_list, fractions, wfunc)

    USE MMD_utilities,           ONLY: sort
    USE MMD_handle_communicator, ONLY: m_model_comm, m_model_rank, m_model_npes

    IMPLICIT NONE

    EXTERNAL  :: MPI_RECV, MPI_SEND
    INTRINSIC :: ASSOCIATED, SIZE

    ! Index list will be sorted
    INTEGER,  DIMENSION(:,:), INTENT(INOUT)        :: index_list
    REAL(dp), DIMENSION(:,:), INTENT(IN)           :: fractions
    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL :: wfunc  ! weight function
    
    ! LOCAL
    INTEGER  :: i,ip,is,ie,ian, ind
    INTEGER  :: status
    INTEGER  :: RemPe
    INTEGER, ALLOCATABLE, DIMENSION(:) :: RemInd
    INTEGER, DIMENSION(2) :: NrEle  
    ! NrEle(1) gives the number of coupled elements 
    ! NrEle(2) indicates if weighting function wfunc is present
    INTEGER  :: num_indpairs

    TYPE indPair
       INTEGER,  DIMENSION(:,:), POINTER :: ij   => NULL()
       REAL(dp), DIMENSION(:),   POINTER :: frac => NULL()
       REAL(dp), DIMENSION(:),   POINTER :: wf   => NULL()
    END TYPE indPair

    TYPE(indPair), DIMENSION(:), POINTER :: par_ind => NULL()

    IF (PRESENT(wfunc)) THEN
       NrEle(2) = 1
    ELSE
       NrEle(2) = 0
    END IF

    ! 1. SPLIT INDEX LIST INTO CHILD SPECIFIC LISTS
    ! the PE with model_rank=0 evaluates the index_list and distributes
    ! the Child Pe specific list to each Child PE
    !-----------------------------------------------------------------------
    IF (m_model_rank == 0)   THEN
       CALL sort (index_list, 5)         ! Sort to ascending Child PE
       is = 1
       DO ip=0,m_model_npes-1
          ! Split into Child PEs
          ie = is-1                      ! there may be no entry for This PE

          IF (is <= SIZE(index_list,2).AND. is >=0 )  THEN
             DO WHILE ( index_list(5,ie+1) == ip)
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
             Parent%NrPoints = ian
             IF (ian > 0)   THEN
                ALLOCATE (Parent%index_list_2d(6,ian))
                Parent%index_list_2d(:,1:ian) = index_list(:,is:ie)
             END IF
          ELSE
             CALL MPI_Send (ian, 1, MPI_INTEGER, ip, 2000, m_model_comm,status)
             IF (ian > 0) THEN
                CALL MPI_Send (index_list(1,is), 6*ian, MPI_INTEGER, ip, 2001,&
                     m_model_comm, status)
             END IF
          END IF
          is = ie+1
       END DO
    ELSE
       CALL MPI_Recv (ian, 1, MPI_INTEGER, 0, 2000, m_model_comm &
                         , MPI_STATUS_IGNORE, status)
       Parent%NrPoints = ian

       IF (ian > 0) THEN
          ALLOCATE(Parent%index_list_2d(6,ian))
          CALL MPI_RECV (Parent%index_list_2d, 6*ian &
               , MPI_INTEGER, 0, 2001, m_model_comm, MPI_STATUS_IGNORE, status)
       END IF
    END IF

    !-----------------------------------------------------------------------
    ! 2. SORT LIST ALONG Child PE numbers
    !-----------------------------------------------------------------------
    ! a) SET NUMBER OF ELEMENTS ON EACH Parent PE TO ZERO
    DO ip=0,Parent%inter_npes-1
       Parent%PEs(ip)%NrEle = 0
    END DO

    ! b) COUNT NUMBER OF ELEMENTS ON EACH Parent PE
    DO i=1,Parent%NrPoints
       RemPe = Parent%index_list_2d(6,i)
       Parent%PEs(RemPe)%NrEle = Parent%PEs(RemPe)%NrEle+1
    END DO

    ! c) EVALUATE REMOTE INDEX PAIRs  
    ALLOCATE(RemInd(0:Parent%inter_npes-1))
    ! INITIALIZE remote index
    DO ip=0,Parent%inter_npes-1
       ALLOCATE(Parent%PEs(ip)%locInd(Parent%PEs(ip)%NrEle))
       RemInd(ip) = 0
    END DO

    ! ALLOCATE ARRAY FOR Parent i,j index Pairs
    ALLOCATE(par_ind(0:Parent%inter_npes-1))
    DO ip=0,Parent%inter_npes-1
       ALLOCATE(par_ind(ip)%ij(Parent%PEs(ip)%NrEle,2))
       ALLOCATE(par_ind(ip)%frac(Parent%PEs(ip)%NrEle))
    ENDDO

    IF (PRESENT(wfunc)) THEN
       DO ip=0,Parent%inter_npes-1
          ALLOCATE(par_ind(ip)%wf(Parent%PEs(ip)%NrEle))
       END DO
    END IF


    ! FIND PARENT and Child i,j Index Pairs per REMOTE PE (RemPe)
    DO i=1,Parent%NrPoints
       RemPe = Parent%index_list_2d(6,i)
       RemInd(RemPe) = RemInd(RemPe) + 1 
       ind   = RemInd(RemPe)
       ! Child INDICES
       Parent%PEs(RemPe)%locInd(ind)%i = Parent%index_list_2d(3,i)
       Parent%PEs(RemPe)%locInd(ind)%j = Parent%index_list_2d(4,i)
       ! Parent INDICES
       par_ind(RemPe)%ij(ind,1) = Parent%index_list_2d(1,i) 
       par_ind(RemPe)%ij(ind,2) = Parent%index_list_2d(2,i)
       ! set weight fraction
       par_ind(RemPe)%frac(ind) = &
            fractions(Parent%index_list_2d(3,i), Parent%index_list_2d(4,i))
       ! set weight function
       IF (PRESENT(wfunc)) &
       par_ind(RemPe)%wf(ind) = &
            wfunc(Parent%index_list_2d(3,i), Parent%index_list_2d(4,i))
     END DO

    ! SEND NUMBER OF ELEMENTS AND Parent i,j-INDEX PAIRS 
    ! WHICH ARE SEND FROM THIS Child PE TO SPECIFIC Parent PE
    DO ip=0,Parent%inter_npes-1
       NrEle(1) = Parent%PEs(ip)%NrEle
       CALL MMD_Send_to_Parent(NrEle(1:2), 2, ip, 1350+ip*5, status)
       num_indpairs = Parent%PEs(ip)%NrEle*2
       CALL MMD_Send_to_Parent(par_ind(ip)%ij, num_indpairs, ip &
            , 1351+ip*5, status)
       CALL MMD_Send_to_Parent(par_ind(ip)%frac, Parent%PEs(ip)%NrEle , ip &
            , 1352+ip*5, status)
       DEALLOCATE(par_ind(ip)%ij)
       DEALLOCATE(par_ind(ip)%frac)
    END DO

    IF (PRESENT(wfunc)) THEN
       DO ip=0,Parent%inter_npes-1
          CALL MMD_Send_to_Parent(par_ind(ip)%wf, Parent%PEs(ip)%NrEle , ip &
               , 1353+ip*5, status)
          DEALLOCATE(par_ind(ip)%wf)
       END DO
    END IF

    DEALLOCATE(par_ind)

    IF (ASSOCIATED(Parent%index_list_2d)) DEALLOCATE(Parent%index_list_2d)
    DEALLOCATE(RemInd)

    RETURN

  END SUBROUTINE MMD_C_Set_ParIndexlist
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_C_Get_ParDataArray_Name(numFields)

    USE mmd_utilities,           ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
    USE MMD_handle_communicator, ONLY: m_to_parent_comm

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    INTEGER, INTENT(OUT)          :: numFields

    ! LOCAL 
    INTEGER                       :: myIndex
    TYPE(ArrayDef_list), POINTER  :: ai => NULL()
    TYPE(ArrayDef_list), POINTER  :: ap => NULL()
    CHARACTER(LEN=STRLEN_OBJECT)  :: chld_object
    CHARACTER(LEN=STRLEN_CHANNEL) :: chld_channel
    CHARACTER(LEN=STRLEN_OBJECT)  :: chld_repr
    INTEGER                       :: interpM
    INTEGER                       :: isentunit 
 
    ! Get Data Array Channel and Name from Parent
    DO 
       CALL MMD_Bcast ( myIndex,  0, comm=m_to_parent_comm)
       IF (myIndex == -1) EXIT

       numFields = myIndex
       CALL MMD_Bcast (chld_channel, 0, comm=m_to_parent_comm)
       CALL MMD_Bcast (chld_object,  0, comm=m_to_parent_comm)
       CALL MMD_Bcast (chld_repr,    0, comm=m_to_parent_comm)
       CALL MMD_Bcast (interpM,      0, comm=m_to_parent_comm)
       CALL MMD_Bcast (isentunit,    0, comm=m_to_parent_comm) 

       ! BUILD LIST DATA ARRAY LIST
       ai => Parent%ArrayStart
       DO 
          IF (.NOT. ASSOCIATED(ai)) EXIT
          ap => ai
          ai => ai%next
       ENDDO
          
       ALLOCATE(ai)
       NULLIFY(ai%next)
       IF (.NOT. ASSOCIATED(Parent%ArrayStart)) THEN
          ! SET POINTER TO FIRST OBJECT
          Parent%ArrayStart => ai  
       ELSE
          ap%next => ai    ! SET NEXT POINTER OF LAST OBJECT TO NEW OBJECT
       ENDIF
       
       ai%ArrDef%channel     = TRIM(chld_channel)
       ai%ArrDef%object      = TRIM(chld_object) 
       ai%ArrDef%repr        = TRIM(chld_repr)
       ai%ArrDef%interpM     = interpM
       IF (isentunit == 1) THEN
          ai%ArrDef%l_sentunit   = .TRUE.
       ELSE
          ai%ArrDef%l_sentunit   = .FALSE.
       END IF

    ENDDO
    
    RETURN

   END SUBROUTINE MMD_C_Get_ParDataArray_Name
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_C_Get_Indexlist

    IMPLICIT NONE

    INTRINSIC :: SIZE

    INTEGER                              :: ip, mr, i
    INTEGER                              :: status
    INTEGER, DIMENSION(1)                :: NrEle
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: indfield

    ! mr is the equivalent to ip in mmd_server. use here to create same tag
    mr = m_model_rank

    ! GET NUMBER OF ELEMENTS PROVIDED BY EACH PARENT PE
    DO ip=0,Me%inter_npes-1
       CALL MMD_Recv_from_Parent(NrEle, 1, ip, 350+mr*2, status)
       IF (status /= 0) RETURN
       Me%PEs(ip)%NrEle = NrEle(1)
       ALLOCATE(Me%PEs(ip)%locInd(Me%PEs(ip)%NrEle))
       ALLOCATE(indfield(Me%PEs(ip)%NrEle,2))
       CALL MMD_Recv_from_Parent(indfield, SIZE(indfield), ip, 351+mr*2, status)
       IF (status /= 0) RETURN
       DO i=1,Me%PEs(ip)%NrEle
          Me%PEs(ip)%locInd(i)%i = indfield(i,1)
          Me%PEs(ip)%locInd(i)%j = indfield(i,2)
       END DO

       DEALLOCATE(indfield)
    ENDDO

    RETURN

  END SUBROUTINE MMD_C_Get_Indexlist
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  LOGICAL FUNCTION MMD_C_GetNextArray (MyChannel, myName)

    INTRINSIC :: ASSOCIATED, TRIM

    CHARACTER(LEN=*), INTENT(OUT)  :: MyChannel
    CHARACTER(LEN=*), INTENT(OUT)  :: myName

    ! LOCAL
    LOGICAL, SAVE                  :: lfirst =  .TRUE.

    IF (lfirst) THEN
       Me%Ar => Me%ArrayStart
    ELSE
       Me%Ar => Me%Ar%next
       IF (.NOT. ASSOCIATED(Me%Ar)) THEN
          MMD_C_GetNextArray = .FALSE.
          Me%Ar => Me%ArrayStart
          lfirst = .TRUE.
          RETURN
       ENDIF
    ENDIF

    MyChannel = TRIM(Me%Ar%ArrDef%channel)
    myName    = TRIM(Me%Ar%ArrDef%object)

    lfirst = .FALSE.
    MMD_C_GetNextArray = .TRUE.

    RETURN

  END FUNCTION MMD_C_GetNextArray
  ! ------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  LOGICAL FUNCTION MMD_C_GetNextParArray (MyChannel, myName, repr, interpM &
       , l_SentUnit) 

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM
  
    CHARACTER(LEN=*),INTENT(OUT) :: myChannel ! channel name
    CHARACTER(LEN=*),INTENT(OUT) :: myName    ! object name
    CHARACTER(LEN=*),INTENT(OUT) :: repr      ! representation
    INTEGER,         INTENT(OUT) :: interpM   ! interpolation method
    LOGICAL,         INTENT(OUT) :: l_SentUnit 

    ! LOCAL
    LOGICAL, SAVE                :: lfirst =.TRUE.

    IF (lfirst) THEN
       Parent%Ar => Parent%ArrayStart
    ELSE
       Parent%Ar => Parent%Ar%next
       IF (.NOT. ASSOCIATED(Parent%Ar)) THEN
          MMD_C_GetNextParArray = .FALSE.
          Parent%Ar => Parent%ArrayStart
          lfirst = .TRUE.
          RETURN
       ENDIF
    ENDIF
    
    MyChannel = TRIM(Parent%Ar%ArrDef%channel)
    myName    = TRIM(Parent%Ar%ArrDef%object)
    repr      = TRIM(Parent%Ar%ArrDef%repr)
    interpM   = Parent%Ar%ArrDef%interpM
    l_SentUnit= Parent%Ar%ArrDef%l_SentUnit 

    lfirst = .FALSE.
    MMD_C_GetNextParArray = .TRUE.

    RETURN 

  END FUNCTION MMD_C_GetNextParArray
  !-------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_C_get_repr(axis, gdimlen, name, att)

    USE  MMD_utilities,            ONLY: STRLEN_MEDIUM, STRLEN_ATT

    IMPLICIT NONE

    INTRINSIC :: SIZE

    CHARACTER(LEN=4),             INTENT(OUT) :: axis
    INTEGER, DIMENSION(4),        INTENT(OUT) :: gdimlen 
    CHARACTER(LEN=STRLEN_MEDIUM), INTENT(OUT) :: name
    CHARACTER(LEN=STRLEN_ATT),    INTENT(OUT) :: att

    ! LOCAL
    INTEGER  :: i

    CALL MMD_Bcast ( axis,0, comm=m_to_parent_comm)
    CALL MMD_Bcast ( name,0, comm=m_to_parent_comm)
    CALL MMD_Bcast ( att, 0, comm=m_to_parent_comm)
    DO i=1,SIZE(gdimlen)
       CALL MMD_Bcast ( gdimlen(i), 0, comm=m_to_parent_comm)
    ENDDO

  END SUBROUTINE MMD_C_get_repr
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_C_Set_DataArray (status, DIMLEN, ArrayOrder, p4)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    INTEGER, INTENT(OUT)                   :: status
    INTEGER, INTENT(IN),DIMENSION(4)       :: DIMLEN
    CHARACTER(LEN=4), INTENT(IN)           :: ArrayOrder
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: p4

    ! LOCAL
    INTEGER                      :: ip, idim
    TYPE(ArrayDef_list), POINTER :: Ar => NULL()
    INTEGER                      :: ArrFac

    Ar => Me%Ar
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
    IF (Me%Ar%ArrDef%xyzn_dim(3) > 0) THEN
       ArrFac = ArrFac * DIMLEN(Me%Ar%ArrDef%xyzn_dim(3))
    END IF
    IF (Me%Ar%ArrDef%xyzn_dim(4) > 0) THEN
       ArrFac = ArrFac * DIMLEN(Me%Ar%ArrDef%xyzn_dim(4))
    END IF

    ! CALCULATE ArrayLength FOR EACH REMOTE PE
    ALLOCATE(Me%Ar%Arrdef%ArrLen(0:Me%inter_npes))
    ALLOCATE(Me%Ar%Arrdef%ArrIdx(0:Me%inter_npes))
    DO ip=0,Me%inter_npes-1
       Me%Ar%Arrdef%ArrLen(ip) = ArrFac * Me%PEs(ip)%NrEle
       Me%Ar%Arrdef%ArrIdx(ip) = BufLen(ip) + 1
       BufLen(ip) = BufLen(ip) + Me%Ar%Arrdef%ArrLen(ip)
    END DO

    IF (ASSOCIATED(p4)) THEN
       Ar%ArrDef%p4 => p4
    ELSE
       status = 100
       RETURN
    ENDIF

    status = 0

   RETURN

 END SUBROUTINE MMD_C_Set_DataArray
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
 SUBROUTINE MMD_C_Set_ParDataArray(status, DimLen, ArrayOrder, p4)

    IMPLICIT NONE

    INTEGER, INTENT(OUT)                   :: status
    INTEGER, INTENT(IN),DIMENSION(4)       :: DimLen
    CHARACTER(LEN=4), INTENT(IN)           :: ArrayOrder
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: p4

    ! LOCAL
    INTEGER                       :: ip, idim
    INTEGER                       :: ArrFac

    Parent%Ar%Arrdef%dim         = DimLen
    Parent%Ar%Arrdef%dim_order   = ArrayOrder
    Parent%Ar%ArrDef%xyzn_dim(:) = -1
    DO idim = 1,4
       IF (ArrayOrder(idim:idim) == 'X') THEN
          Parent%Ar%ArrDef%xyzn_dim(1) = idim
       ELSE IF (ArrayOrder(idim:idim) == 'Y') THEN
          Parent%Ar%ArrDef%xyzn_dim(2) = idim
       ELSE IF (ArrayOrder(idim:idim) == 'Z') THEN
          Parent%Ar%ArrDef%xyzn_dim(3) = idim
       ELSE IF (ArrayOrder(idim:idim) == 'N') THEN
          Parent%Ar%ArrDef%xyzn_dim(4) = idim
       ENDIF
    END DO
    ArrFac = 1
    IF (Parent%Ar%ArrDef%xyzn_dim(3) > 0) THEN
       ArrFac = ArrFac *DimLen(Parent%Ar%ArrDef%xyzn_dim(3))
    ENDIF
    IF (Parent%Ar%ArrDef%xyzn_dim(4) > 0) THEN
       ArrFac = ArrFac *DimLen(Parent%Ar%ArrDef%xyzn_dim(4))
    ENDIF

    ! ALLOCATE ARRAY FOR DIMENSION OF EACH ARRAY TRANSFERED FROM PARENT
    ! TO CHILD
    ALLOCATE(Parent%Ar%Arrdef%ArrLen(0:Parent%inter_npes))
    DO ip = 0,Parent%inter_npes-1
       Parent%Ar%Arrdef%ArrLen(ip) =  &
            ArrFac * Parent%PEs(ip)%NrEle
       ParBufLen(ip) = ParBufLen(ip) + Parent%Ar%Arrdef%ArrLen(ip)
    END DO

    Parent%Ar%Arrdef%p4 => p4
    
    status = 0

    RETURN

  END SUBROUTINE MMD_C_Set_ParDataArray
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_C_setInd_and_AllocMem 

    IMPLICIT NONE

    EXTERNAL :: MMDc_C_SetInd_and_Mem_2way, MMDc_C_SetInd_and_Mem 

    IF (ltwoway) THEN
       CALL MMDc_C_SetInd_and_Mem_2way(BufLen,ParBufLen)
    ELSE
       CALL MMDc_C_SetInd_and_Mem(BufLen) 
    ENDIF

    RETURN

  END SUBROUTINE MMD_C_setInd_and_AllocMem
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE MMD_C_GetBuffer(WaitTime)

    USE  MMD_utilities,            ONLY: PeDef, ArrayDef_list

    IMPLICIT NONE

    EXTERNAL  :: MMDC_C_GETBUFFER, MMDC_C_GETWAITTIME_2WAY
    INTRINSIC :: ASSOCIATED, PRESENT

    REAL(kind=dp),INTENT(OUT), OPTIONAL        :: WaitTime

    ! LOCAL
    TYPE(ArrayDef_list), POINTER    :: Ar  => NULL()
    TYPE(PeDef),         POINTER    :: aPE => NULL()
    INTEGER, DIMENSION(4)           :: aind
    INTEGER                         :: ind, ij, ip, k,n, idim
    REAL(DP), DIMENSION(:), POINTER :: array => NULL()
    REAL(dp)                        :: my_waittime
    INTEGER                         :: zloopend, nloopend

    ! This not only get the WaitTime but also includes on barrier Call
    IF (ltwoway) THEN
       CALL MMDc_C_GetWaitTime_2way(my_WaitTime)
    ELSE
       CALL MMDc_C_GetWaitTime(my_WaitTime)
    END IF

    IF (PRESENT(WaitTime)) WaitTime = my_WaitTime

    PE_loop: DO ip=0,Me%inter_npes-1
       Me%Ar => Me%ArrayStart
       Ar    => Me%Ar

       ALLOCATE(array(BufLen(ip)))
       array = 0._dp ! ub_ak_20170711

       aPE =>  Me%PEs(ip)
       CALL MMDc_C_GetBuffer(ip, array)

       ind = 0

       arrayloop: DO
          !write (0,*) 'Child: ',ip, TRIM(Ar%Arrdef%object) &
          !     , Ar%Arrdef%dim_order &
          !     , Ar%Arrdef%xyzn_dim, Ar%Arrdef%dim  
          ind = Ar%Arrdef%ArrIdx(ip)
          zloopend=1
          nloopend=1
 
         IF (ip==0) Ar%Arrdef%p4(:,:,:,:) = 0._dp ! ub_ak_20170711

          IF (Ar%Arrdef%dim_order == 'XY--') THEN
             DO ij = 1, aPE%NrEle
                aind(1) = aPE%locInd(ij)%i
                aind(2) = aPE%locInd(ij)%j
                Ar%Arrdef%p4(aind(1),aind(2),1,1) = array(ind)
                ind     = ind+1
             END DO
          ELSE IF (Ar%Arrdef%dim_order == 'X-Y-') THEN
             DO ij = 1, aPE%NrEle
                aind(1) = aPE%locInd(ij)%i
                !aind(2) = 1
                aind(3) = aPE%locInd(ij)%j
                Ar%Arrdef%p4(aind(1),1,aind(3),1) = array(ind)
                ind     = ind+1
             END DO
          ELSE IF (Ar%Arrdef%dim_order == 'XYZ-') THEN
             DO ij = 1, aPE%NrEle
                DO k=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(3))
                   aind(1) = aPE%locInd(ij)%i
                   aind(2) = aPE%locInd(ij)%j
                   aind(3) = k
                   Ar%Arrdef%p4(aind(1),aind(2),aind(3),1) = array(ind)
                   ind     = ind+1
                END DO
             END DO
          ELSE IF(Ar%ArrDef%dim_order == 'XZY-') THEN
             DO ij = 1, aPE%NrEle
                DO k=1,Ar%Arrdef%dim(Ar%ArrDef%xyzn_dim(3))
                   aind(1) = aPE%locInd(ij)%i
                   aind(2) = k
                   aind(3) = aPE%locInd(ij)%j
                   Ar%Arrdef%p4(aind(1),aind(2),aind(3),1)= array(ind)
                   ind = ind + 1  
                END DO
             END DO
          ELSE IF (Ar%Arrdef%dim_order == 'XYNZ') THEN
             DO ij = 1, aPE%NrEle
                DO k=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(3))
                   DO n=1, Ar%Arrdef%dim(Ar%Arrdef%xyzn_dim(4))
                      aind(1) = aPE%locInd(ij)%i
                      aind(2) = aPE%locInd(ij)%j
                      aind(3) = n
                      aind(4) = k
                      Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) = array(ind)
                      ind     = ind+1
                   END DO
                END DO
             END DO
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
                      Ar%Arrdef%p4(aind(1),aind(2),aind(3),aind(4)) = array(ind) 

                      ind     = ind + 1 
                   END DO
                END DO
             ENDDO
          ENDIF

          Ar => Ar%next
          IF (.NOT. ASSOCIATED(Ar)) EXIT ! exit array loop
       END DO Arrayloop
       DEALLOCATE(array)
    END DO PE_loop

    !IF (.NOT. ltwoway) CALL MMDc_C_SetBarrier() 
    RETURN

  END SUBROUTINE MMD_C_GetBuffer
  ! ------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_C_FillBuffer

    USE MMD_utilities,           ONLY: PeDef

    IMPLICIT NONE

    EXTERNAL  :: MMDC_C_FILLBUFFER, MMDC_C_SETBARRIER
    INTRINSIC :: ASSOCIATED
    ! LOCAL
    TYPE(ArrayDef_list), POINTER   :: Ar  => NULL()
    TYPE(PeDef),         POINTER   :: aPE => NULL()
    INTEGER, DIMENSION(4)          :: aind  ! Index in Pointer
    INTEGER                        :: ij, ip, idim, k, n, ind
    INTEGER                        :: zloopend, nloopend

    REAL(DP), DIMENSION(:), POINTER :: array => NULL()
    REAL(DP)                        :: my_waittime

    ! This not only get the WaitTime but also includes on barrier Call
!!$    my_WaitTime = -1.0
!!$    CALL MMDc_C_GetWaitTime(my_WaitTime)
!!$    IF (PRESENT(WaitTime))  then
!!$      WaitTime = my_WaitTime
!!$    end if
!!$    CALL MMDc_C_SetBarrier

    PE_loop: DO ip=0,Parent%inter_npes-1
       IF (.NOT. ASSOCIATED(Parent%ArrayStart)) EXIT ! op_mm_20150326 

       Parent%Ar => Parent%ArrayStart
       Ar => Parent%Ar
       aPE =>  Parent%PEs(ip)
       ALLOCATE(array(ParBufLen(ip)))
       ind = 0

       arrayloop: DO
             !IF (TRIM(Ar%Arrdef%object) == 'ALLTRACER') &
!!$             write (0,*) 'CHILD FILL: ',m_model_rank, ip, TRIM(Ar%Arrdef%object)&
!!$                  , '  Order ',Ar%Arrdef%dim_order  &
!!$                  , ' ArrLen ',Ar%Arrdef%ArrLen(ip)     &
!!$                  , ' XYZN_DIM ',Ar%Arrdef%xyzn_dim &
!!$                  , ' DIM ',Ar%Arrdef%dim &
!!$                  , 'BD ',LBOUND(Ar%Arrdef%p4), UBOUND(Ar%Arrdef%p4)
          IF (Ar%Arrdef%ArrLen(ip) == 0) THEN
             !NOTHING TO DO
          ELSE IF (Ar%ArrDef%dim_order == 'XY--') THEN
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

       CALL MMDc_C_FillBuffer(ip, array)

       DEALLOCATE(array)

    END DO PE_loop

    ! Set Barrier when Buffer Full
    CALL MMDc_C_SetBarrier()

    RETURN

  END SUBROUTINE MMD_C_FillBuffer
  !-------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  INTEGER FUNCTION MMD_C_GetParentType()

    USE  MMD_handle_communicator,  ONLY: m_ParentType

    IMPLICIT NONE

    MMD_C_GetParentType = m_ParentType

    RETURN

  END FUNCTION MMD_C_GetParentType
  ! ------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE MMD_C_FreeMem

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL
    
    ! LOCAL
    INTEGER                      :: ip
    TYPE(ArrayDef_list), POINTER :: ai => NULL()
    TYPE(ArrayDef_list), POINTER :: ae => NULL()

    CALL MMDc_C_FreeMem ! ub_ak_20180202

    DO ip=0,Me%inter_npes-1
       DEALLOCATE(Me%PEs(ip)%locInd)
    END DO
    DEALLOCATE(Me%Ar%Arrdef%ArrLen)
    DEALLOCATE(Me%PEs)

    ai => Me%ArrayStart

    DO 
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    
       DEALLOCATE(ae)
       NULLIFY(ae)
    END DO

    ! Parent 

    IF (ltwoway) THEN

       DO ip=0,Parent%inter_npes-1
          IF (ASSOCIATED(Parent%PEs(ip)%locInd))& 
               DEALLOCATE(Parent%PEs(ip)%locInd)
       END DO

       IF (ASSOCIATED(Parent%Ar))& 
            DEALLOCATE(Parent%Ar%Arrdef%ArrLen)
       DEALLOCATE(Parent%PEs) 

       ai => Parent%ArrayStart

       DO 
          IF (.NOT. ASSOCIATED(ai)) EXIT
          ae => ai
          ai => ai%next
          
          DEALLOCATE(ae)
          NULLIFY(ae)
       END DO
    END IF

  END SUBROUTINE MMD_C_FreeMem
  !-------------------------------------------------------------------------

END MODULE mmd_child
