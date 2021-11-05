MODULE mod_oasis_flush
 
#ifdef NAGFOR
   USE F90_unix,      ONLY: flush
   USE F90_unix_PROC, ONLY: abort
#endif

  USE kinds_mod     ! defines common data types
  !
  IMPLICIT NONE

   private

   public oasis_flush_scrip
   public oasis_abort_scrip
!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------

!==========================================================================
   SUBROUTINE oasis_flush_scrip(nu)


   IMPLICIT NONE

!--------------------------------------------------------------------
   INTEGER(kind=int_kind),INTENT(in) :: nu
!--------------------------------------------------------------------
   character(len=*),parameter :: subname = 'oasis_flush_scrip'
!--------------------------------------------------------------------

   CALL FLUSH(nu)

 END SUBROUTINE oasis_flush_scrip

!==========================================================================

!==========================================================================
   SUBROUTINE oasis_abort_scrip


   IMPLICIT NONE

!--------------------------------------------------------------------
!   character(len=*),parameter :: subname = 'oasis_abort_scrip'
!--------------------------------------------------------------------

   CALL ABORT

 END SUBROUTINE oasis_abort_scrip

!==========================================================================

END MODULE mod_oasis_flush
