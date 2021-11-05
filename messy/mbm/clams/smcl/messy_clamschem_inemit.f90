Module messy_clamschem_inemit

contains

!     inemit   - Initialise emissions scheme
!
!     written by V. Cals 
!     last modification: 25-jun-2001

  

subroutine inemit

  use messy_clamschem_asad_mod,       only: nlemit, lemit, speci
  use messy_clamschem_asad_mod_clams, only: jpspec

  implicit none

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
   
  integer :: i
  
  do i=1, jpspec
     if( trim(speci(i)) == 'NO' ) then
        nlemit(1)=i
        lemit(i)=.true.
     else 
        lemit(i)=.false.
     endif
  end do
  
  return
end subroutine inemit


end Module messy_clamschem_inemit

