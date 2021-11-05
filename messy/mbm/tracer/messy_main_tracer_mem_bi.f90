! **********************************************************************+
MODULE messy_main_tracer_mem_bi
! **********************************************************************+

  ! THIS MODULE PROVIDES THE BML-SPECIFIC AND BML-INDEPENDENT
  ! 'TRACER' INFORMATION TO THE MESSy-SUBMODELS

  USE messy_main_tracer, ONLY: t_trinfo_list, t_trinfo_tp, dp

  IMPLICIT NONE
  INTRINSIC :: NULL
  !PUBLIC is already default
  PRIVATE :: dp
  SAVE

  ! TRACER SET #1 ##################################################
  !
  ! NAME OF TRACER SET
  CHARACTER(LEN=*), PARAMETER              :: S1TRSTR = 's1'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_s1 => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_s1 => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_s1 = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxt   => NULL()
!!$  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtte => NULL()
!!$  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtm1 => NULL()
!!$  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtf  => NULL()
  !
  ! - FOR OUTSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xt    => NULL()
!!$  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtte  => NULL()
!!$  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtm1  => NULL()
!!$  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtf   => NULL()
  !
  ! - FOR INSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxt   => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtte => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtm1 => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtf  => NULL()
  !
  ! #######################################################################

! **********************************************************************+
END MODULE messy_main_tracer_mem_bi
! **********************************************************************+
