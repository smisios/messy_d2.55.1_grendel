!
! Print value of all #define keywords
!

SUBROUTINE print_defines(io_unit)
  IMPLICIT NONE

  INTEGER, INTENT(in) :: io_unit

#ifndef MESSY
#include "defines.list"
#endif

END SUBROUTINE print_defines
