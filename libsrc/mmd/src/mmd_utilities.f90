MODULE mmd_utilities

  USE      mpi

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: NULL, SELECTED_INT_KIND, SELECTED_REAL_KIND
  
  ! from messy_main_constants_mem.f90
  INTEGER, PARAMETER, PUBLIC  :: MMD_DP = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER, PUBLIC  :: MMD_I8 = SELECTED_INT_KIND(14)

  INTEGER, PARAMETER, PUBLIC  :: DP=MMD_DP

  ! Length of Data Array Name 
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_CHANNEL = 23     
  ! Length of Data Array Name 
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_OBJECT  = 55     
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_MEDIUM  = 24
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_ULONG   = 256
  ! Length of attribute for representation identification
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_ATT     = 800

  ! Definition PARAMETER
  INTEGER, PARAMETER, PUBLIC  :: MMD_ParentIsECHAM = 1
  INTEGER, PARAMETER, PUBLIC  :: MMD_ParentIsCOSMO = 2

  ! return status
  INTEGER, PARAMETER, PUBLIC  :: MMD_STATUS_OK    = 0
  INTEGER, PARAMETER, PUBLIC  :: MMD_DA_NAME_ERR  = 10
  ! Maximum number of coupled models
  INTEGER, PARAMETER, PUBLIC  :: MMD_MAX_MODEL   = 64

  TYPE ArrayDef
     CHARACTER(LEN=STRLEN_CHANNEL)          :: channel = '' ! Name of Channel  
     CHARACTER(LEN=STRLEN_OBJECT)           :: object  = '' ! Name of Object
     CHARACTER(LEN=STRLEN_MEDIUM)           :: repr    = '' ! repr of Object
     ! interpolation Method
     INTEGER                                :: interpM      
     LOGICAL                                :: l_sentunit ! um_ak_20150413
     ! DATA pointer
     REAL(DP), POINTER, DIMENSION(:,:,:,:)  :: p4 => NULL()

     CHARACTER(LEN=4)      :: dim_order = '    ' ! Order of dimensions
     INTEGER, DIMENSION(4) :: xyzn_dim  = 0      ! index of x,y,z,n dimension
     INTEGER, DIMENSION(4) :: dim       = 0      ! Size of Dimensions
     ! ArrLen and ArrIdx are different on each remote PE
     ! Dimension of Array moved
     INTEGER, DIMENSION(:), POINTER :: ArrLen  => NULL() 
     ! Start Index of Array moved
     INTEGER, DIMENSION(:), POINTER :: ArrIdx  => NULL()
  END TYPE ArrayDef

  TYPE ArrayDef_list
     TYPE(ArrayDef)               :: arrdef
     TYPE(ArrayDef_list), POINTER :: next => NULL()
  END TYPE ArrayDef_list

  PUBLIC :: ArrayDef, ArrayDef_list
  
  ! Pair of indices in horizontal plane
  TYPE xy_ind
     INTEGER  :: i
     INTEGER  :: j
     REAL(dp) :: frac    ! weight fraction 
  END TYPE xy_ind
  PUBLIC:: xy_ind

  TYPE PeDef
     INTEGER                             :: NrEle     ! Number of Elements 
     TYPE(xy_ind), POINTER, DIMENSION(:) :: locInd => NULL()
  END TYPE PeDef
  PUBLIC :: PeDef

  TYPE ExchDataDef
    ! NUMBER OF REMOTE MODEL PEs
    INTEGER                            :: inter_npes = 0  
    ! REMOTE MODEL ID
    INTEGER                            :: RMId   = 0
    ! STRUCTURE FOR EACH PE
    TYPE(PeDef), DIMENSION(:), POINTER :: PEs => NULL()
    ! INDEX LIST OF PARENT MODEL POINTS
    INTEGER,   DIMENSION(:,:), POINTER :: index_list_2d => NULL()
    ! Number of Points in index_list
    INTEGER                            :: NrPoints = 0
    ! ARRAY INFORMATION STRUCTURE (SAME ON ALL PEs)
    TYPE(ArrayDef_list), POINTER       :: Ar         => NULL()
    TYPE(ArrayDef_list), POINTER       :: ArrayStart => NULL()
  END TYPE ExchDataDef
  PUBLIC :: ExchDataDef

  INTERFACE sort
    MODULE PROCEDURE sort_2d_i
  END INTERFACE sort

  INTERFACE get_Wtime
    MODULE PROCEDURE get_Wtime
  END INTERFACE get_Wtime

  PUBLIC ::  sort, get_Wtime

 CONTAINS

  !--------------------------------------------------------
  SUBROUTINE sort_2d_i (array,sort_ind)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    INTEGER,DIMENSION(:,:),INTENT(INOUT)         :: array
    INTEGER,INTENT(IN)                           :: sort_ind

    ! LOCAL
    INTEGER                                      :: i,j,n
    INTEGER,DIMENSION(SIZE(array,1))             :: tmp

    n = SIZE(array,2)
    DO j=1,n-1
      DO i=j+1,n
        IF (array(sort_ind,i) < array(sort_ind,j) )  THEN
          tmp = array(:,i)
          array(:,i) = array(:,j)
          array(:,j) = tmp
        END IF
      END DO
    END DO 

    RETURN

  END  SUBROUTINE sort_2d_i
  !--------------------------------------------------------

  !--------------------------------------------------------
  REAL(kind=DP) FUNCTION get_Wtime()

    IMPLICIT NONE

    INTRINSIC :: system_clock, REAL

    ! LOCAL
    INTEGER           :: counter,count_rate
    REAL(KIND=dp)     :: t

    CALL system_clock (count=counter,count_rate=count_rate)

    t         = REAL(counter,dp)
    t         = t/REAL(count_rate,dp)
    get_Wtime = t

    RETURN

  END FUNCTION get_Wtime 

END MODULE mmd_utilities
