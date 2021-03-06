
Module messy_clamsdeepconv_global
  
  USE messy_clams_global,      ONLY: PREC

  IMPLICIT NONE
  
  REAL, PARAMETER :: lv = 2.26e+6
  REAL, PARAMETER :: cp = 1005.




  
  !---------------------------------------------------------------------------
  ! Namelist variables
  !---------------------------------------------------------------------------
  REAL(PREC)       :: cbvf       = 0.
  REAL(PREC)       :: crate      = 1.
  REAL(PREC)       :: min_dtheta = 0.
  REAL(PREC)       :: DC_kind    = 0.

  INTEGER          :: timestep_deepconv=0   ! timestep in hours 
  



  REAL(PREC) :: zeta_min


  LOGICAL :: dc_mix = .false.
  
  
End Module messy_clamsdeepconv_global
