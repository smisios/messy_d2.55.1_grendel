      SUBROUTINE AVRG_TIMESER_BGC
!*******************************************************************
!
!**** *AVRG_BGCMEAN* - average timeser_bgc data.
!
!     Patrick Wetzel,        *MPIMet, HH*    1.05.03
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
!     *CALL*       *AVRG_TIMESER_BGC*
!
!     *MODULES*     *mo_timeser_bgc* - bgc timeseries parameters.
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



      USE mo_timeser_bgc
      USE mo_control_bgc
implicit none

      INTEGER :: l1,l2
      REAL :: rfreqts1

      WRITE(io_stdo_bgc,*)'Averaging at step',ldtrunbgc,lts1

      rfreqts1=1./FLOAT(nfreqts1)
      DO l2=1,nelets1
      DO l1=1,nvarts1
         ts1(l1,l2,lts1) = ts1(l1,l2,lts1)*rfreqts1
      ENDDO
      ENDDO
!
! Increment timeseries sample counter (initialized to 1 in INI_TIMSER_BGC).
!
      IF ( lts1 .LT. lents1 ) THEN
        lts1 = lts1 + 1
        WRITE(io_stdo_bgc,*)'Sample counter of timeseries 1 is ',lts1
      ENDIF
      RETURN
      END
