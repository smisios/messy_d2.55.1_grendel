      MODULE MO_ELICOM

      USE mo_kind, ONLY: wp
      IMPLICIT NONE

      REAL(wp), POINTER :: PGL(:,:)

      CONTAINS

      SUBROUTINE alloc_mem_elicom

      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL

      IF (p_pe==p_io) THEN
         ALLOCATE (PGL(IMM,ILL))
      ELSE
         ALLOCATE (PGL(0,0))
      ENDIF

      END SUBROUTINE alloc_mem_elicom
      END MODULE MO_ELICOM
