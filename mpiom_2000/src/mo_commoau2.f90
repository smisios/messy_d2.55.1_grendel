!>
!! @ingroup common
!!
      MODULE MO_COMMOAU2
      USE MO_PARAM1
      IMPLICIT NONE
!==>COMMOAU2
     REAL(wp), ALLOCATABLE,TARGET ::PRECO(:,:)
     REAL(wp), ALLOCATABLE,TARGET ::QSWO(:,:)
     REAL(wp), ALLOCATABLE,TARGET ::QLWO(:,:)
     REAL(wp), ALLOCATABLE,TARGET ::QSEO(:,:)
     REAL(wp), ALLOCATABLE,TARGET ::QLAO(:,:)
     REAL(wp), ALLOCATABLE,TARGET ::SICUDO(:,:)
     REAL(wp), ALLOCATABLE,TARGET ::SICVDE(:,:)
     REAL(wp), ALLOCATABLE,TARGET ::PRECH(:,:)
!<==END COMMOAU2

      CONTAINS

      SUBROUTINE alloc_mem_commoau2

      ALLOCATE(                                                         &
     &  PRECO(IE,JE),QSWO(IE,JE),QLWO(IE,JE),QSEO(IE,JE)                &
     & ,QLAO(IE,JE)                                                     &
     & ,SICUDO(IE,JE),SICVDE(IE,JE),PRECH(IE,JE))

      END SUBROUTINE alloc_mem_commoau2
      END MODULE MO_COMMOAU2
