      SUBROUTINE CLOSE_BGCMEAN_2D
!****************************************************************
!
!**** *CLOSE_BGCMEAN_2D* - calculate and write 2-dimensional bgc mean data.
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
!     *CALL*       *CLOSE_BGCMEAN_2D*
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
      INTEGER :: ncstat,ncoldmod

!-----------------------------------------------------------------------

!
! Set fill mode
!
!-----------------------------------------------------------------------
      IF (p_pe==p_io) THEN
       ncstat = NF_SET_FILL(nc_2d_id,NF_NOFILL, ncoldmod)
       IF ( ncstat .NE. NF_NOERR ) call stop_all('CLOSE_BGCMEAN_2D: Problem with netCDF97')
      END IF

!
! Write bgcmean data : 2-D data
!
!-----------------------------------------------------------------------

!!$      CALL write_netcdf_var(nc_2d_id,'co214flux',bgct2d(1,1,jco214f),1,0)
!!$      CALL write_netcdf_var(nc_2d_id,'co2flux',bgct2d(1,1,jco2flux),1,0)
!!$      CALL write_netcdf_var(nc_2d_id,'n2oflux',bgct2d(1,1,jn2oflux),1,0)
!!$      CALL write_netcdf_var(nc_2d_id,'prorca',bgct2d(1,1,jprorca),1,0)
!!$      CALL write_netcdf_var(nc_2d_id,'prcaca',bgct2d(1,1,jprcaca),1,0)
!!$      CALL write_netcdf_var(nc_2d_id,'silpro',bgct2d(1,1,jsilpro),1,0)
!!$      CALL write_netcdf_var(nc_2d_id,'produs',bgct2d(1,1,jprodus),1,0)

!!$#ifdef PANTHROPOCO2
!!$      CALL write_netcdf_var(nc_2d_id,'co2fant',bgct2d(1,1,jco2fant),1,0)
!!$      CALL write_netcdf_var(nc_2d_id,'cacadiff',bgct2d(1,1,jcacadiff),1,0)
!!$#endif
!!$#ifdef ANTC14
!!$      CALL write_netcdf_var(nc_2d_id,'rantc14',rbomb(1,1),1,0)
!!$#endif

!
! Close File
!
      if (p_pe==p_io) then
        ncstat = NF_CLOSE(nc_2d_id)
        IF ( ncstat .NE. NF_NOERR ) call stop_all('CLOSE_BGCMEAN_2D: Problem with netCDF200'  )
      end if

      RETURN
      END
