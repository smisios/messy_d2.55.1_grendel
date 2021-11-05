      SUBROUTINE WRITE_BGCMEAN_SED(kpie,kpje,kpke,pddpo)
!****************************************************************
!
!     under construction!!!!!!!!!!!!!!!!
!
!**** *WRITE_BGCMEAN_SED* - calculate and write 3-dimensional bgc sediment mean data.
!
!     Patrick Wetzel,    *MPI-Met, HH*    24.05.04
!
!     Modified
!     --------
!     Joachim Segschneider, *MPI-Met, HH* 02.08.05
!     output of sediment data, derived from WRITE_BGCMEAN_3D
!
!     Purpose
!     -------
!     Write bgc mean sediment data.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *WRITE_BGCMEAN_SEDI(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,ptiestu*
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - zonal dimension of model grid.
!     *INTEGER* *kpje*    - meridional dimension of model grid.
!     *INTEGER* *kpke*    - vertical dimension of model grid.
!     *REAL*    *pddpo*   - depth of grid cell [m].
!     *REAL*    *pdlxp*   - zonal size of scalar grid cell  [m].
!     *REAL*    *pdlyp*   - meridional size of scalar grid cell [m].
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

      USE mo_control_bgc
      use mo_commo1, only: lyears,lmonts,ldays

      use mo_parallel

      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,l,nk
      REAL(wp) :: pddpo(kpie,kpje,kpke)
      REAL(wp) :: timesed
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
      COUNT_3D(3) =ks
      COUNT_3D(4) =1

!-----------------------------------------------------------------------


      WRITE(io_stdo_bgc,*) 'Writing sediment bgcmean data at step:', meantime_3d

!
!  Masking bgcmean data.
!
!-----------------------------------------------------------------------

      DO l=1,nbgct_sed
      DO k=1,ks
      DO j=1,kpje
      DO i=1,kpie
!       IF(bolay(i,j) .LT. 0.5) THEN
        IF (pddpo(i, j, 1) .LT. 0.5_wp) THEN
          bgct_sed(i,j,k,l)=rmasko
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

      ncstat = NF_SET_FILL(nc_sed_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('BGCMEAN: Problem with netCDF97')

       timesed=REAL(LDAYS+LMONTS*100+LYEARS*10000, wp)

       ncstat = NF_INQ_VARID(nc_sed_id,'time',ncvarid )
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_3D: Problem with netCDF10sed')
       ncstat = NF_PUT_VARA_DOUBLE (nc_sed_id,ncvarid,start_1d,count_1d,timesed)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_3D: Problem with netCDF10sed')

      end if

!
! Write bgcmean data : sediment (pore water and solid constituents)
!
!-----------------------------------------------------------------------


      nk = ks
      CALL write_netcdf_var(nc_sed_id,'powaic',bgct_sed(1,1,1,jpowaic),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'powaal',bgct_sed(1,1,1,jpowaal),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'powaph',bgct_sed(1,1,1,jpowaph),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'powaox',bgct_sed(1,1,1,jpowaox),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'pown2',bgct_sed(1,1,1,jpown2),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'powno3',bgct_sed(1,1,1,jpowno3),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'powasi',bgct_sed(1,1,1,jpowasi),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'ssso12',bgct_sed(1,1,1,jssso12),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'ssssil',bgct_sed(1,1,1,jssssil),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'sssc12',bgct_sed(1,1,1,jsssc12),nk,meantime_3d)
      CALL write_netcdf_var(nc_sed_id,'ssster',bgct_sed(1,1,1,jssster),nk,meantime_3d)

!
! Reset mean fields
!
!-----------------------------------------------------------------------
!
      DO l=1,nbgct_sed
      DO k=1,ks
      DO j=1,kpje
      DO i=1,kpie
          bgct_sed(i,j,k,l) = 0._wp
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
