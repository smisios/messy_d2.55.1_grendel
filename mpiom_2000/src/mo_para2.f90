MODULE MO_PARA2
  USE mo_kind, ONLY: wp
  USE MO_PARAM1
  IMPLICIT NONE

!UWE
!     AUXILIARY FIELDS FOR ITERATIVE SOLUTION OF BAROTROPIC MATRIX
!     INITIALIZED IN ITPREP

  REAL(wp), ALLOCATABLE :: uf(:,:), vf(:,:), ff(:,:), xx(:,:)

  REAL(wp), ALLOCATABLE :: uf0(:,:),vf0(:,:)

  REAL(wp) :: sorpar

! included in namelist ocectl

  INTEGER :: iter_sor=300 ! default number of iterations in the sor scheme
  REAL(wp) :: rtsorpar=-999._wp ! user specified value for sor parameter , overides sorpar determined in sbr trotest

  !> default number of iterations with rtsorpar_hack , usefull in
  !> higher resolution setup to dampen 2dx noise
  INTEGER :: iter_sor_hack=0
  REAL(wp) :: rtsorpar_hack=-999._wp !user specified value for sor parameter in the last ITER_SOR_HACK steps of the SOR


CONTAINS

  SUBROUTINE alloc_mem_para2

    USE mo_commo1, ONLY : lwith_barotropic_stokes_drift

    ALLOCATE(UF(IE,JE),VF(IE,JE),FF(IE,JE),XX(IE,JE))

    uf=0._wp
    vf=0._wp
    ff=0._wp
    xx=0._wp

    IF ( lwith_barotropic_stokes_drift ) THEN

      ALLOCATE(uf0(ie,je),vf0(ie,je))
      uf0=0._wp
      vf0=0._wp

    ENDIF

  END SUBROUTINE alloc_mem_para2

END MODULE MO_PARA2
