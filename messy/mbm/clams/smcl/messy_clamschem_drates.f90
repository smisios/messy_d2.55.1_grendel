Module messy_clamschem_drates


Contains


! drates.f90 : subroutines that read in reaction rates and rate constants 
!
! 25.08.1999  changed to array calculations instead of single values(jug)


subroutine get_bconst(ibr, values) 

  USE messy_main_constants_mem, ONLY: DP

  use messy_clamschem_global,   only: ntraj
  use messy_clamschem_asad_mod, only: rk, nbrkx

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: ibr 
  real(kind=DP), dimension(ntraj), intent(out) :: values 
!-----------------------------------------------
  values = rk(:ntraj,nbrkx(ibr)) 
  return  
end subroutine get_bconst


 
subroutine get_tconst(itr, values) 

  USE messy_main_constants_mem, ONLY: DP

  use messy_clamschem_global,   only: ntraj
  use messy_clamschem_asad_mod, only: rk, ntrkx

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: itr 
  real(kind=DP), dimension(ntraj), intent(out) :: values 
!-----------------------------------------------
  values = rk(:ntraj,ntrkx(itr)) 
  return  
end subroutine get_tconst 


 
subroutine get_jconst(ijr, values) 

  USE messy_main_constants_mem, ONLY: DP

  use messy_clamschem_global,   only: ntraj
  use messy_clamschem_asad_mod, only: rk, nprkx

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: ijr 
  real(kind=DP), dimension(ntraj), intent(out) :: values 
!-----------------------------------------------
  values = rk(:ntraj,nprkx(ijr)) 
  return  
end subroutine get_jconst


 
subroutine get_hconst(ihr, values) 

  USE messy_main_constants_mem, ONLY: DP

  use messy_clamschem_global,   only: ntraj
  use messy_clamschem_asad_mod, only: rk, nhrkx

  implicit none
       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: ihr 
  real(kind=DP), dimension(ntraj), intent(out) :: values 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!D      value = hk(itrj,ihr)*(abs(y(itrj,scale_index(ihr)))+eps)
!jug    value = rk(itrj,nhrkx(ihr))*(abs(f(itrj,scale_index(ihr)))+eps)
  values = rk(:ntraj,nhrkx(ihr)) 
  return  
end subroutine get_hconst 


 
subroutine get_brates(ibr, values) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:07  11/23/98  
!...Switches: -ymhf               
  USE messy_main_constants_mem, ONLY: DP
  use messy_clamschem_global,   only: ntraj
  use messy_clamschem_asad_mod, only: rk, nbrkx, nspi, y
  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: ibr 
  real(kind=DP), dimension(ntraj), intent(out) :: values 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: ibb 
!-----------------------------------------------
  ibb = nbrkx(ibr) 
  values = rk(:ntraj,nbrkx(ibr))*y(:ntraj,nspi(ibb,1))*y(:ntraj,nspi(ibb,2)) 
  return  
end subroutine get_brates


 
subroutine get_trates(itr, values) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:07  11/23/98  
!...Switches: -ymhf               
!     Purpose: calculate termolecular reaction rates of the reaction e.g. A+B+M->products
!     k*[A]*[B]*[M] for output only
!     06.08.1997    corrected since reaction rate constant tk already includes k*[M]
!                   before there was a factor tnd(itrj)

  USE messy_main_constants_mem, ONLY: DP
  use messy_clamschem_global,   only: ntraj, therm_flag
  use messy_clamschem_asad_mod, only: rk, ntrkx, nspi, y

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: itr 
  real(kind=DP), dimension(ntraj), intent(out) :: values 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: itt 
!-----------------------------------------------
  itt = ntrkx(itr) 
 
  if (therm_flag(itr) == 0) then 
     values = rk(:ntraj,ntrkx(itr))*y(:ntraj,nspi(itt,1))*y(:ntraj,nspi(itt,2)) 
  else if (therm_flag(itr) == 1) then 
     values = rk(:ntraj,ntrkx(itr))*y(:ntraj,nspi(itt,1)) 
  endif
 
  return  
end subroutine get_trates 


 
subroutine get_jrates(ijr, values) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:07  11/23/98  
!...Switches: -ymhf               
  USE messy_main_constants_mem, ONLY: DP
  use messy_clamschem_global,   only: ntraj
  use messy_clamschem_asad_mod, only: rk, nprkx, nspi, y

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: ijr 
  real(kind=DP), dimension(ntraj), intent(out) :: values 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: ijj 
!-----------------------------------------------
  ijj = nprkx(ijr) 
  values = rk(:ntraj,nprkx(ijr))*y(:ntraj,nspi(ijj,1)) 
  return  
end subroutine get_jrates 


 
subroutine get_hrates(ihr, values) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:07  11/23/98  
!...Switches: -ymhf               
  USE messy_main_constants_mem, ONLY: DP
  use messy_clamschem_global,   only: ntraj
  use messy_clamschem_asad_mod, only: rk, nhrkx, nspi, y

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: ihr 
  real(kind=DP), dimension(ntraj), intent(out) :: values 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: ihh 
!-----------------------------------------------
  ihh = nhrkx(ihr) 
  values = rk(:ntraj,nhrkx(ihr))*y(:ntraj,nspi(ihh,1))*y(:ntraj,nspi(ihh,2)) 
  return  
end subroutine get_hrates 
 
End Module messy_clamschem_drates
