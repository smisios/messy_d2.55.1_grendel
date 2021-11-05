      SUBROUTINE AVRG_BGCMEAN_2D(kpie,kpje,kpke)
!*******************************************************************
!
!**** *AVRG_BGCMEAN* - average bgcmean data.
!
!     Patrick Wetzel,        *MPIMet, HH*    24.05.04
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     time-average 2d and bioz fields
!
!     Method
!     -------
!
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *AVRG_BGCMEAN_2D(kpie,kpje,kpke)*
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
      use mo_param1_bgc

implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k,l

      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_2D using:   mean_2D_freq='  &
                                                      ,mean_2D_freq
      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_2D using:    meantime_2d='  &
                                                      ,meantime_2d
      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_2D using: meancnt_bgc_2D='&
                                                     ,meancnt_bgc_2D

      WRITE(io_stdo_bgc,*) 'AVRG_BGCMEAN_2D using: ',              &
     &  meantime_2d,meancnt_bgc_2D,kpje,kpie,nbgcm2d,nbgcm3d,kwrbioz

      DO k=1,nbgcm2d
      DO j=1,kpje
      DO i=1,kpie
        bgcm2d(i, j, k) = bgcm2d(i, j, k) / REAL(meancnt_bgc_2D, wp)
      ENDDO
      ENDDO
      ENDDO

      DO k=1,nbgct2d
      DO j=1,kpje
      DO i=1,kpie
        bgct2d(i, j, k) = bgct2d(i, j, k) / REAL(meancnt_bgc_2D, wp)
      ENDDO
      ENDDO
      ENDDO

      DO l=1,nbgcm3d
      DO k=1,kwrbioz
      DO j=1,kpje
      DO i=1,kpie
        bgcm3d(i, j, k, l) = bgcm3d(i, j, k, l) / REAL(meancnt_bgc_2D, wp)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

! Save no. of steps

!      stepspm(meantime_2d) = meancnt_bgc_2D

! Reset counter

      meancnt_bgc_2D=0

      RETURN
      END
