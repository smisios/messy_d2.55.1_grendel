PROGRAM msbm
! DESCRIPTION
! box model version of the ECHAM5/MESSY msbm module
!
! AUTHOR
! J. Buchholz, Max Planck Institute for Chemistry, Mainz, Germany
!
! LAST CHANGES
! 11. August 2004 by J. Buchholz!
  
  USE messy_msbm_box, ONLY: msbm_initialize, msbm_physc, msbm_output
  
  IMPLICIT NONE

  CALL msbm_initialize
  CALL msbm_physc
  CALL msbm_output

END PROGRAM msbm
