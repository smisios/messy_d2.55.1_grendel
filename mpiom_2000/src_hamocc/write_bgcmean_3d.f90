      SUBROUTINE WRITE_BGCMEAN_3D(kpie,kpje,kpke,pddpo)
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
!     *CALL*       *WRITE_BGCMEAN_3D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,ptiestu*
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL :: of model grid.
!     *INTEGER* *kpje*    - 2nd REAL :: of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL :: of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd REAL ::) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st REAL ::) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd REAL ::) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu* - depth of layers [m].
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
      use mo_commo1, only: lyears,lmonts,ldays

      USE mo_control_bgc

      use mo_parallel

      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,l,nk
      REAL(wp) :: pddpo(kpie,kpje,kpke)
      REAL(wp) :: time3d

      INTEGER :: START_3D(4),COUNT_3D(4) ! size of internal array
      INTEGER :: START_1D,COUNT_1D


      INCLUDE 'netcdf.inc'
      INTEGER :: ncvarid,ncstat,ncoldmod


!-----------------------------------------------------------------------

      START_1D =meantime_2d
      COUNT_1D =1

      START_3D(1) =1
      START_3D(2) =1
      START_3D(3) =1
      START_3D(4) =meantime_3d

      COUNT_3D(1) =kpie
      COUNT_3D(2) =kpje
      COUNT_3D(3) =kpke
      COUNT_3D(4) =1

!-----------------------------------------------------------------------


      WRITE(io_stdo_bgc,*) 'Writing 3D bgcmean data at step:', meantime_3d

!
!  Masking bgcmean data.
!
!-----------------------------------------------------------------------

      DO l=1,nbgct3d
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        IF (pddpo(i,j,k) .LT. 0.5_wp) THEN
          bgct3d(i,j,k,l)=rmasko
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO


!
! Set fill mode
!
!-----------------------------------------------------------------------
      if (p_pe==p_io) then

      ncstat = NF_SET_FILL(nc_3d_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('BGCMEAN: Problem with netCDF97')

       time3d=REAL(LDAYS+LMONTS*100+LYEARS*10000, wp)


       ncstat = NF_INQ_VARID(nc_3d_id,'time',ncvarid )
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_3D: Problem with netCDF103d')
       ncstat = NF_PUT_VARA_DOUBLE (nc_3d_id,ncvarid,start_1d,count_1d,time3d)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_3D: Problem with netCDF103d')

      end if
!
! Write bgcmean data : total 3-D bgc data
!
!-----------------------------------------------------------------------
      nk = kpke
      CALL write_netcdf_var(nc_3d_id,'phosph_t',bgct3d(1,1,1,jphosph_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'nitrate_t',bgct3d(1,1,1,jano3_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'silica_t',bgct3d(1,1,1,jsilica_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'iron_t',bgct3d(1,1,1,jiron_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'oxygen_t',bgct3d(1,1,1,joxygen_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'alkali_t',bgct3d(1,1,1,jalkali_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'dic_t',bgct3d(1,1,1,jdic_t),nk,meantime_3d)
#ifdef __c_isotopes
      CALL write_netcdf_var(nc_3d_id,'dic13_t',bgct3d(1,1,1,jdic13_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'dic14_t',bgct3d(1,1,1,jdic14_t),nk,meantime_3d)
#endif
      CALL write_netcdf_var(nc_3d_id,'doc_t',bgct3d(1,1,1,jdoc_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'poc_t',bgct3d(1,1,1,jpoc_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'calc_t',bgct3d(1,1,1,jcalc_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'opal_t',bgct3d(1,1,1,jopal_t),nk,meantime_3d)

#ifdef ANTC14
      CALL write_netcdf_var(nc_3d_id,'ac14_t',bgct3d(1,1,1,jac14_t),nk,meantime_3d)
#endif
#ifdef PCFC
      CALL write_netcdf_var(nc_3d_id,'cfc11_t',bgct3d(1,1,1,jcfc11_t),nk,meantime_3d)
      CALL write_netcdf_var(nc_3d_id,'cfc12_t',bgct3d(1,1,1,jcfc12_t),nk,meantime_3d)
#endif



!
! Reset mean fields
!
!-----------------------------------------------------------------------
!
      DO l=1,nbgct3d
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
          bgct3d(i, j, k, l) = 0._wp
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
