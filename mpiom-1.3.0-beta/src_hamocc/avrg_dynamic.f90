      SUBROUTINE AVRG_DYNAMIC(kpie,kpje,kpke)
!*******************************************************************
!
!**** *AVRG_DYNAMIC* - average dynamic data.
!
!     Patrick Wetzel,        *MPIMet, HH*    30.05.04
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *AVRG_DYNAMIC(kpie,kpje,kpke)*
!
!     *MODULES*     *mo_carbch* - ocean/sediment tracer arrays.
!     *MODULES*     *mo_control_bgc*  - std I/O logical units.
!     *MODULES*     *mo_bgcmean
!     *MODULES*     *mo_param1_bgc 
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************


      USE mo_carbch
      USE mo_control_bgc
      USE mo_bgcmean
      USE mo_dynamic
      use mo_param1_bgc 

implicit none    

      INTEGER :: kpie,kpje,kpke,i,j,k,l


      WRITE(io_stdo_bgc,*) 'AVRG_DYNAMIC mit: ',              &
     &  meantime_2d,meancnt_bgc_2D,kpje,kpie,nbgcm2d

      DO l=1,kdtot
      DO k=1,nbgcdyn
      DO j=1,kpje
      DO i=1,kpie
         bgcdyn(i,j,k,l)=bgcdyn(i,j,k,l)/meancnt_bgc_2D
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO j=1,kpje
      DO i=1,kpie
         bgc_zmld(i,j) = bgc_zmld(i,j)/meancnt_bgc_2D
         bgc_nmld(i,j) = bgc_nmld(i,j)/meancnt_bgc_2D
      ENDDO
      ENDDO


      RETURN
      END
