      SUBROUTINE OPEN_BGCMEAN_BIOZ(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce)
!****************************************************************
!
!**** *OPEN_BGCMEAN_BIOZ* - calculate and write 2-dimensional bgc mean data.
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
!     *CALL*       *OPEN_BGCMEAN_BIOZ(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
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
      INTEGER :: ncvarid,ncstat,ncoldmod,ncdims(4),timeid
      INTEGER :: nclatid,nclonid,nclevid,ncbiozid
      INTEGER :: nstart2(2),ncount2(2),nstride2(2)
      integer :: icode(1)
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

!      WRITE(io_stdo_bgc,*) 'OPEN_BGCMEAN_BIOZ: ',kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
!     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce

!
! Open netCDF data file
!
!-----------------------------------------------------------------------

!      ncstat = NF_CREATE('bgcmean_bioz.nc',NF_CLOBBER, nc_bioz_id)
      ncstat = NF_CREATE('bgcmean_bioz.nc',IOR(NF_CLOBBER,NF_64BIT_OFFSET), nc_bioz_id)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF1')

      chunksize   = 1024*1024*32 ! 32 MB --> man 3f netcdf
      initialsize = 1024*1024*32 ! 32 MB --> man 3f netcdf

!      ncstat = nf__create('bgcmean_bioz.nc',NF_CLOBBER,initialsize,chunksize,nc_bioz_id)
!      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF1')

!
! Define dimension
!
!-----------------------------------------------------------------------

      ncstat = NF_DEF_DIM(nc_bioz_id, 'lon', ie_g, nclonid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF2')

      ncstat = NF_DEF_DIM(nc_bioz_id, 'lat', je_g, nclatid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF3')

      ncstat = NF_DEF_DIM(nc_bioz_id, 'depth', kpke, nclevid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF4')

      ncstat = NF_DEF_DIM(nc_bioz_id, 'kwrbioz', kwrbioz, ncbiozid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF5')

      ncstat = NF_DEF_DIM(nc_bioz_id, 'time',NF_UNLIMITED, timeid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF7'      )

!
! Define global attributes
!
!-----------------------------------------------------------------------

      ncstat = NF_PUT_ATT_TEXT(nc_bioz_id, NF_GLOBAL,'title'             &
     &,42, 'Monthly mean output for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF8')

      ncstat = NF_PUT_ATT_TEXT(nc_bioz_id, NF_GLOBAL,'history'           &
     &,42, 'Monthly mean output for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(nc_bioz_id, NF_GLOBAL,'conventions'       &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF10')

      ncstat = NF_PUT_ATT_TEXT(nc_bioz_id, NF_GLOBAL,'source'            &
     &,24, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF11')

      ncstat = NF_PUT_ATT_INT(nc_bioz_id, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF11')

!
! Define variables : grid
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid

      CALL NETCDF_DEF_VARSG(nc_bioz_id,8,'scal_lon',2,ncdims,ncvarid,          &
     &   8,'degree E',34,'2-d longitude of scalar grid cells',rmasko,20,io_stdo_bgc)
     icode=1
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_bioz_id,8,'scal_lat',2,ncdims,ncvarid,          &
     &   8,'degree N',33,'2-d latitude of scalar grid cells',rmasko,21,io_stdo_bgc)
     icode=2
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_bioz_id,9,'scal_wdep',2,ncdims,ncvarid,          &
     &   5,'meter',37,'2-d water depth at scalar grid points',rmasko,22,io_stdo_bgc)
     icode=3
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'size_x',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in x-direction',rmasko,23,io_stdo_bgc)
     icode=4
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'size_y',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in y-direction',rmasko,24,io_stdo_bgc)
     icode=5
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid

      CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'size_z',3,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in z-direction',rmasko,25,io_stdo_bgc)
     icode=6
      ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

      ncdims(1) = nclevid

      CALL NETCDF_DEF_VARSG(nc_bioz_id,5,'depth',1,ncdims,ncvarid,          &
     &   5,'meter',32,'1-d layer depths of ocean tracer',rmasko,26,io_stdo_bgc)


!
! Define variables : 1-D data
!
!-----------------------------------------------------------------------

      ncdims(1) = timeid

!      CALL NETCDF_DEF_VARSG(nc_bioz_id,9,'steps_p_m',1,ncdims,ncvarid,          &
!     &   5,'steps',25,'model timesteps per month',rmasko,30,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(nc_bioz_id,4,'time',1,ncdims,ncvarid,          &
     &   16,'day as %Y%m%d.%f',4,'time',rmasko,31,io_stdo_bgc)

!
! Define variables : bioz bgc mean data
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = ncbiozid
      ncdims(4) = timeid

       CALL NETCDF_DEF_VARSG(nc_bioz_id,5,'phyto',4,ncdims,ncvarid,        &
     &    10,'kmolP/m**3',27,'Phytoplankton concentration',          &
     &    rmasko,34,io_stdo_bgc)
     icode=22
       ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'grazer',4,ncdims,ncvarid,       &
     &    10,'kmolP/m**3',25,'Zooplankton concentration',            &
     &    rmasko,36,io_stdo_bgc)
     icode=23
       ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

!       CALL NETCDF_DEF_VARSG(nc_bioz_id,3,'poc',4,ncdims,ncvarid,          &
!     &    9,'kmol/m**3',26,'Particulate organic carbon',             &
!     &    rmasko,38,io_stdo_bgc)

       CALL NETCDF_DEF_VARSG(nc_bioz_id,5,'phosy',4,ncdims,ncvarid,        &
     &    14,'kmol/(dt*m**3)',14,'Photosynthesis',                    &
     &    rmasko,40,io_stdo_bgc)
     icode=100
       ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

#ifdef FB_BGC_OCE
       CALL NETCDF_DEF_VARSG(nc_bioz_id,5,'atten',4,ncdims,ncvarid,        &
     &    3,'1/m',34,'short wave attenuation coefficient',           &
     &    rmasko,42,io_stdo_bgc)
     icode=102
       ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)
#endif

       CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'phosph',4,ncdims,ncvarid,       &
     &    9,'kmol/m**3',19,'Dissolved phosphate',                    &
     &    rmasko,44,io_stdo_bgc)
     icode=11
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'oxygen',4,ncdims,ncvarid,       &
     &    9,'kmol/m**3',20,'oxygen concentration',                   &
     &    rmasko,46,io_stdo_bgc)
     icode=12
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_bioz_id,4,'iron',4,ncdims,ncvarid,          &
     &    12,'kmol Fe/m**3',18,'Iron concentration',                 &
     &    rmasko,47,io_stdo_bgc)
     icode=31
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_bioz_id,7,'nitrate',4,ncdims,ncvarid,      &
     &    9,'kmol/m**3',17,'Dissolved nitrate',                      &
     &    rmasko,44,io_stdo_bgc)
     icode=14
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)
       CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'alkali',4,ncdims,ncvarid,       &
     &    9,'kmol/m**3',16,'total alkalinity',                       &
     &    rmasko,44,io_stdo_bgc)
     icode=10
      ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'silica',4,ncdims,ncvarid,       &
     &    9,'kmol/m**3',18,'Dissolved silicate',                     &
     &    rmasko,44,io_stdo_bgc)
     icode=15
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_bioz_id,3,'dic',4,ncdims,ncvarid,          &
     &    9,'kmol/m**3',26,'Dissolved inorganic carbon',             &
     &    rmasko,44,io_stdo_bgc)
     icode=7
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

!      CALL NETCDF_DEF_VARSG(nc_bioz_id,5,'dic13',4,ncdims,ncvarid,          &
!    &    9,'kmol/m**3',28,'Dissolved inorganic carbon13',             &
!    &    rmasko,44,io_stdo_bgc)
!      CALL NETCDF_DEF_VARSG(nc_bioz_id,5,'dic14',4,ncdims,ncvarid,          &
!    &    9,'kmol/m**3',28,'Dissolved inorganic carbon14',             &
!    &    rmasko,44,io_stdo_bgc)


       CALL NETCDF_DEF_VARSG(nc_bioz_id,3,'doc',4,ncdims,ncvarid,          &
     &    9,'kmol/m**3',24,'Dissolved organic carbon',               &
     &    rmasko,44,io_stdo_bgc)
     icode=16
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)

!       CALL NETCDF_DEF_VARSG(nc_bioz_id,3,'dms',4,ncdims,ncvarid,          &
!     &    9,'kmol/m**3',15,'DiMethylSulfate',                        &
!     &    rmasko,44,io_stdo_bgc)
!
!       CALL NETCDF_DEF_VARSG(nc_bioz_id,4,'calc',4,ncdims,ncvarid,         &
!     &    9,'kmol/m**3',17,'calcium carbonate',                      &
!     &    rmasko,44,io_stdo_bgc)
!
!       CALL NETCDF_DEF_VARSG(nc_bioz_id,4,'opal',4,ncdims,ncvarid,         &
!     &    9,'kmol/m**3',15,'biogenic silica',                        &
!     &    rmasko,44,io_stdo_bgc)
!
!       CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'export',4,ncdims,ncvarid,       &
!     &    9,'kmol/m**3',26,'detritus export production',             &
!     &    rmasko,44,io_stdo_bgc)
!
!       CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'expoca',4,ncdims,ncvarid,       &
!     &    9,'kmol/m**3',25,'calcium export production',              &
!     &    rmasko,44,io_stdo_bgc)
!
!       CALL NETCDF_DEF_VARSG(nc_bioz_id,6,'exposi',4,ncdims,ncvarid,       &
!     &    9,'kmol/m**3',22,'opal export production',                 &
!     &    rmasko,44,io_stdo_bgc)

#ifdef AGG
       CALL NETCDF_DEF_VARSG(nc_bioz_id,7,'numbers',4,ncdims,ncvarid,      &
     &    12,'nos/dt*cm**3',32,'Number of marine snow aggregates',    &
     &    rmasko,48,io_stdo_bgc)
     icode=101
     ncstat = NF_PUT_ATT_INT(nc_bioz_id,ncvarid, 'code', NF_INT, 1, icode)
#endif

!
! END Define variables
!
!-----------------------------------------------------------------------

      ncstat = NF_ENDDEF(nc_bioz_id)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF00')
!
! Set fill mode
!
!-----------------------------------------------------------------------

      ncstat = NF_SET_FILL(nc_bioz_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF97')

      END IF ! p_io==p_pe
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
      CALL write_netcdf_var(nc_bioz_id,'scal_lon',zfield(1,1),1,0)

      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgiph(ii,jj)*aradtogra
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_bioz_id,'scal_lat', zfield(1,1),1,0)

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

      CALL write_netcdf_var(nc_bioz_id,'scal_wdep',zfield(1,1),1,0)

      IF(p_pe==p_io) THEN
       ncstat = NF_INQ_VARID(nc_bioz_id,'depth',ncvarid )
       IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF102a')
       ncstat = NF_PUT_VAR_DOUBLE (nc_bioz_id,ncvarid,ptiestu(1) )
       IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF103a')
      END IF

      CALL write_netcdf_var(nc_bioz_id,'size_x',pdlxp(1,1),1,0)
      CALL write_netcdf_var(nc_bioz_id,'size_y',pdlyp(1,1),1,0)
      CALL write_netcdf_var(nc_bioz_id,'size_z',pddpo(1,1,1),kpke,0)

! Close File
!
!      ncstat = NF_CLOSE(nc_bioz_id)
!      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_BIOZ: Problem with netCDF200')


      RETURN
      END
