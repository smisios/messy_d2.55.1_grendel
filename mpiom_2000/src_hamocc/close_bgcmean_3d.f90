      SUBROUTINE CLOSE_BGCMEAN_3D
!****************************************************************
!
!**** *WRITE_BGCMEAN_2D* - calculate and write 3-dimensional bgc mean data.
!
!     Patrick Wetzel,    *MPI-Met, HH*    24.05.04
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Write bgc mean data.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *CLOSE_BGCMEAN_3D*
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


      WRITE(io_stdo_bgc,*) 'Close 3D bgcmean data at step:', meantime_3d


!
! Close File
!
!-----------------------------------------------------------------------
!
      ncstat = NF_CLOSE(nc_3d_id)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('CLOSE_BGCMEAN_3D: Problem with netCDF200')

      end if


      RETURN
      END
