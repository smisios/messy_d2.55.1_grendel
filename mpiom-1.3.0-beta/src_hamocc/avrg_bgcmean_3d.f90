      SUBROUTINE AVRG_BGCMEAN_3D(kpie,kpje,kpke)
!*******************************************************************
!
!**** *AVRG_BGCMEAN* - average bgcmean data.
!
!     Patrick Wetzel,        *MPIMet, HH*    24.05.04
!
!     Modified
!     --------
!     js 02.08.2005 average sediment tracers
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
!     *CALL*       *AVRG_BGCMEAN_3D(kpie,kpje,kpke)*
!
!     *MODULES*     *mo_carbch* - ocean/sediment tracer arrays.
!     *MODULES*     *mo_control_bgc*  - std I/O logical units.
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
      use mo_param1_bgc 

implicit none    

      INTEGER :: kpie,kpje,kpke,i,j,k,l




      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_3D using:   mean_3D_freq='  &
                                                      ,mean_3D_freq
      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_3D using:    meantime_3d='  &
                                                      ,meantime_3d   
      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_3D using: meancnt_bgc_3D='&
                                                     ,meancnt_bgc_3D




      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_3D at step: ',meantime_3d
      WRITE(io_stdo_bgc,*) 'avrg. total 3D fields with ', meancnt_bgc_3D
        
      DO l=1,nbgct3d
        DO k=1,kpke
        DO j=1,kpje
        DO i=1,kpie      
           bgct3d(i,j,k,l)=bgct3d(i,j,k,l)/meancnt_bgc_3D
        ENDDO
        ENDDO
        ENDDO
        ENDDO

!js average sediment
      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_SED at step: ',meantime_3d
      WRITE(io_stdo_bgc,*) 'avrg. total sed. fields with ', meancnt_bgc_3D

        DO l=1,nbgct_sed
        DO k=1,ks
        DO j=1,kpje
        DO i=1,kpie
           bgct_sed(i,j,k,l)=bgct_sed(i,j,k,l)/meancnt_bgc_3D
        ENDDO
        ENDDO
        ENDDO
        ENDDO

! Reset counter  
    
      meancnt_bgc_3D=0


      RETURN
      END
