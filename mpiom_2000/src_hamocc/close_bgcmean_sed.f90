      SUBROUTINE CLOSE_BGCMEAN_SED
#ifdef PBGC      
!****************************************************************
!
!**** *WRITE_BGCMEAN_SED - calculate and write 3-dimensional bgc mean data.
!
!     Patrick Wetzel,    *MPI-Met, HH*    24.05.04
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Close bgcmean_sed.nc
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *CLOSE_BGCMEAN_SED*
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
      use mo_param1_bgc 

      USE mo_control_bgc
      USE mo_parallel
 
      implicit none


      INCLUDE 'netcdf.inc'
      INTEGER :: ncstat


      
!-----------------------------------------------------------------------

      if (p_pe==p_io) then

      WRITE(io_stdo_bgc,*) 'Close sediment bgcmean file at step:', meantime_3d

!
! Close File
!
!-----------------------------------------------------------------------
!
      ncstat = NF_CLOSE(nc_sed_id)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('CLOSE_BGCMEAN_SED: Problem with netCDF200')

      end if

      RETURN
#endif
      END
