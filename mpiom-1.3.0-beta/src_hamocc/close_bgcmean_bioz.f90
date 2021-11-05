      SUBROUTINE CLOSE_BGCMEAN_BIOZ

!****************************************************************
!
!**** *CLOSE_BGCMEAN_BIOZ* - close bgc bioz mean data file.
!
!     Patrick Wetzel,    *MPI-Met, HH*    24.05.04
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     close bgc mean data file.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *CLOSE_BGCMEAN_BIOZ(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
!                  ptiestu) *
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

      USE mo_bgcmean
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_param1_bgc 

      USE mo_control_bgc
      USE mo_parallel
 
      IMPLICIT NONE

      INCLUDE 'netcdf.inc'
      INTEGER :: ncstat

!-----------------------------------------------------------------------
      
      IF (p_pe==p_io) THEN

!         WRITE(0,*) 'CLOSE_BGCMEAN_BIOZ ',nc_bioz_id,ncstat

         !
         ! Close File
         !

         ncstat = NF_CLOSE(nc_bioz_id)

!         WRITE(0,*) 'END OF CLOSE_BGCMEAN_BIOZ ',ncstat
         
         IF ( ncstat .NE. NF_NOERR ) CALL STOP_all('CLOSE_BGCMEAN_BIOZ: Problem with netCDF200'  )
         
!         WRITE(0,*) 'END OF CLOSE_BGCMEAN_BIOZ ',ncstat

      END IF
      !




      RETURN
      END




