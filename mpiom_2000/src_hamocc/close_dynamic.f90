      SUBROUTINE CLOSE_DYNAMIC
!****************************************************************
!
!**** *CLOSE_DYNAMIC* - calculate and write 2-dimensional bgc mean data.
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
!     *CALL*       *CLOSE_DYNAMIC
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

      USE mo_bgcmean
      use mo_param1_bgc
      USE mo_dynamic

      USE mo_control_bgc
      USE mo_parallel

      implicit none

      INCLUDE 'netcdf.inc'
      INTEGER :: ncstat

!
! Close File
!

      if (p_pe==p_io) then

        ncstat = NF_CLOSE(nc_dyn_id)
        IF ( ncstat .NE. NF_NOERR ) STOP 'CLOSE_DYNAMIC: Problem with netCDF200'

      end if

      RETURN
      END
