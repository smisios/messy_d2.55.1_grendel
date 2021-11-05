      SUBROUTINE OPEN_DYNAMIC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce)
!****************************************************************
!
!**** *OPEN_DYNAMIC* - calculate and write 2-dimensional bgc mean data.
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
!     *CALL*       *OPEN_DYNAMIC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
!                  pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce) *
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
      USE mo_dynamic

      USE mo_control_bgc

      USE mo_commo1, ONLY: ie_g, je_g
      USE mo_constants, ONLY: aradtogra
      USE mo_parallel

      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,ii,jj
      INTEGER :: kplyear,kplmon,kplday,kpldtoce
      REAL(wp) :: pddpo(kpie,kpje,kpke)
      REAL(wp) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL(wp) :: pgila(kpie*2,kpje*2)
      REAL(wp) :: pgiph(kpie*2,kpje*2)
      REAL(wp) :: ptiestu(kpke)


      INCLUDE 'netcdf.inc'
      INTEGER :: ncvarid,ncstat,ncoldmod,ncdims(4)
      INTEGER :: nclatid,nclonid
      INTEGER :: nstart2(2),ncount2(2),nstride2(2)
      INTEGER :: nbgcdynid, timeid

      INTEGER :: idate(5)
      INTEGER :: chunksize,initialsize
      REAL(wp) :: zfield(kpie,kpje)

!-----------------------------------------------------------------------

      idate(1) = kplyear
      idate(2) = kplmon
      idate(3) = kplday
      idate(4) = kpldtoce
      idate(5) = ldtbgc

      IF (p_pe==p_io) THEN

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Create 2D bgcmean data at date :'    &
     &  ,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_bgc,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_bgc,*) 'Bgc   model step number is ',idate(5)

!      WRITE(io_stdo_bgc,*) 'OPEN_DYNAMIC: ',kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
!     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce

!
! Open netCDF data file
!
!-----------------------------------------------------------------------

      ncstat = NF_CREATE('dynamic_bgc.nc',NF_CLOBBER, nc_dyn_id)
      ncstat = NF_CREATE('dynamic_bgc.nc',IOR(NF_CLOBBER,NF_64BIT_OFFSET), nc_dyn_id)
!      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF1'

      chunksize   = 1024*1024*32 ! 32 MB --> man 3f netcdf
      initialsize = 1024*1024*32 ! 32 MB --> man 3f netcdf

      ncstat = nf__create('dynamic_bgc.nc',NF_CLOBBER,initialsize,chunksize,nc_dyn_id)
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF1'

!
! Define dimension
!
!-----------------------------------------------------------------------

      ncstat = NF_DEF_DIM(nc_dyn_id, 'lon', ie_g, nclonid)
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF2'

      ncstat = NF_DEF_DIM(nc_dyn_id, 'lat', je_g, nclatid)
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF3'

      ncstat = NF_DEF_DIM(nc_dyn_id, 'time',NF_UNLIMITED, timeid)
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF7'

      ncstat = NF_DEF_DIM(nc_dyn_id, 'dyn', nbgcdyn, nbgcdynid)
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF5'

!
! Define global attributes
!
!-----------------------------------------------------------------------

      ncstat = NF_PUT_ATT_TEXT(nc_dyn_id, NF_GLOBAL,'title'             &
     &,42, 'Mean monthly output for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF8'

      ncstat = NF_PUT_ATT_TEXT(nc_dyn_id, NF_GLOBAL,'history'           &
     &,42, 'Mean monthly output for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF9'

      ncstat = NF_PUT_ATT_TEXT(nc_dyn_id, NF_GLOBAL,'conventions'       &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF10'

      ncstat = NF_PUT_ATT_TEXT(nc_dyn_id, NF_GLOBAL,'source'            &
     &,24, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF11'

      ncstat = NF_PUT_ATT_INT(nc_dyn_id, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF11'

!
! Define variables : grid
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid

      CALL NETCDF_DEF_VARSG(nc_dyn_id,8,'scal_lon',2,ncdims,ncvarid,          &
     &   8,'degree E',34,'2-d longitude of scalar grid cells',rmasko,20,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,8,'scal_lat',2,ncdims,ncvarid,          &
     &   8,'degree N',34,'2-d longitude of scalar grid cells',rmasko,21,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,9,'scal_wdep',2,ncdims,ncvarid,          &
     &   5,'meter',37,'2-d water depth at scalar grid points',rmasko,22,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,6,'size_x',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in x-direction',rmasko,23,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,6,'size_y',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in y-direction',rmasko,24,io_stdo_bgc)

!
! Define variables : 2-D mean data
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = timeid

      CALL NETCDF_DEF_VARSG(nc_dyn_id,4,'nmld',3,ncdims,ncvarid,       &
     &   5,'layer',17, 'no. of mld layers',rmasko,40,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,4,'zmld',3,ncdims,ncvarid,         &
     &   5,'meter',17, 'mixed layer depth',rmasko,41,io_stdo_bgc)


!
! Define variables : 3-D mean data
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nbgcdynid
      ncdims(4) = timeid

      CALL NETCDF_DEF_VARSG(nc_dyn_id,3,'adv',4,ncdims,ncvarid,   &
     &   9,'kmol/m**2',9, 'advection',rmasko,50,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,3,'dif',4,ncdims,ncvarid,   &
     &   9,'kmol/m**2',9, 'diffusion',rmasko,50,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,3,'pre',4,ncdims,ncvarid,   &
     &   9,'kmol/m**2',6, 'fwflux',rmasko,50,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,3,'gmp',4,ncdims,ncvarid,   &
     &   9,'kmol/m**2',17, 'Gent & McWilliams',rmasko,50,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(nc_dyn_id,3,'bio',4,ncdims,ncvarid,   &
     &   9,'kmol/m**2',7, 'biology',rmasko,50,io_stdo_bgc)

!
! END Define variables
!
!-----------------------------------------------------------------------

      ncstat = NF_ENDDEF(nc_dyn_id)
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF00'
!
! Set fill mode
!
!-----------------------------------------------------------------------

      ncstat = NF_SET_FILL(nc_dyn_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF97'

      END IF ! p_pe==p_io

!
! Write grid describing data
!
!-----------------------------------------------------------------------

      nstart2(1) = 1
      nstart2(2) = 1
      ncount2(1) = kpie
      ncount2(2) = kpje
      nstride2(1) = 2
      nstride2(2) = 2
      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgila(ii,jj)*aradtogra
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_dyn_id,'scal_lon',zfield(1,1),1,0)

      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgiph(ii,jj)*aradtogra
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_dyn_id,'scal_lat', zfield(1,1),1,0)

      DO j=1,kpje
         DO i=1,kpie
            zfield(i, j) = 0.0_wp
         ENDDO
      ENDDO
      DO k=1,kpke
         DO j=1,kpje
            DO i=1,kpie
               zfield(i,j) = zfield(i,j) + pddpo(i,j,k)
            ENDDO
         ENDDO
      ENDDO

      CALL write_netcdf_var(nc_dyn_id,'size_x',pdlxp(1,1),1,0)
      CALL write_netcdf_var(nc_dyn_id,'size_y',pdlyp(1,1),1,0)

!
! Close File
!
!      ncstat = NF_CLOSE(nc_dyn_id)
!      IF ( ncstat .NE. NF_NOERR ) STOP 'OPEN_DYNAMIC: Problem with netCDF200'


      RETURN
      END
