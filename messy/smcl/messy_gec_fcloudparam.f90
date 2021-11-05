module messy_gec_fcloudparam

use messy_main_constants_mem, only:r8=>dp

implicit none

public :: f_cloudparam,fprime_cloudparam


contains

real(r8) function f_cloudparam(x,R,d,a,b,pver)
    implicit none
    integer, intent(in) :: pver
    real(r8), intent(in) :: x,R,d(pver),a(pver),b(pver)

    f_cloudparam = R-sum(d/(a+x*b))

  end function f_cloudparam


real(r8) function fprime_cloudparam(x,R,d,a,b,pver)
    implicit none
    integer, intent(in) :: pver
    real(r8), intent(in) :: x,R,d(pver),a(pver),b(pver)
    
    fprime_cloudparam = -sum((-d*b)/((a+x*b)**2))

  end function fprime_cloudparam
end module messy_gec_fcloudparam
