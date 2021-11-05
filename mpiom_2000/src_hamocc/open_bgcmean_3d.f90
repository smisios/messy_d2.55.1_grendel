      SUBROUTINE OPEN_BGCMEAN_3D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce)

!****************************************************************
!
!**** OPEN_BGCMEAN_3D* - open netcdf files for 3-dimensional bgc mean data.
!
!     Patrick Wetzel,    *MPI-Met, HH*    24.05.04
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     open netcdf files for 3-dimensional bgc mean data
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *OPEN_BGCMEAN_3D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
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
      INTEGER :: ncvarid,ncstat,ncoldmod,ncdims(4)
      INTEGER :: nclatid,nclonid,nclevid
      INTEGER :: timeid,nstart2(2),ncount2(2),nstride2(2)
      INTEGER :: icode(1)
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
      WRITE(io_stdo_bgc,*) 'Creating 3D bgcmean data at date :'    &
     &  ,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_bgc,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_bgc,*) 'Bgc   model step number is ',idate(5)


!
! Open netCDF data file
!
!-----------------------------------------------------------------------

!      ncstat = NF_CREATE('bgcmean_3d.nc',NF_CLOBBER, nc_3d_id)
       ncstat = NF_CREATE('bgcmean_3d.nc',IOR(NF_CLOBBER,NF_64BIT_OFFSET), nc_3d_id)

      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF1')

      chunksize = 1024*1024*32 ! 32 MB --> man 3f netcdf
      initialsize = 1024*1024*32 ! 32 MB --> man 3f netcdf

!      ncstat = nf__create('bgcmean_3d.nc',NF_CLOBBER,initialsize,chunksize,nc_3d_id)
!      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF1')
!
! Define dimension
!
!-----------------------------------------------------------------------

      ncstat = NF_DEF_DIM(nc_3d_id, 'lon', ie_g, nclonid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF2')

      ncstat = NF_DEF_DIM(nc_3d_id, 'lat', je_g, nclatid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF3')

      ncstat = NF_DEF_DIM(nc_3d_id, 'depth', kpke, nclevid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF4')

      ncstat = NF_DEF_DIM(nc_3d_id, 'time',NF_UNLIMITED, timeid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF7'      )
!
! Define global attributes
!
!-----------------------------------------------------------------------

      ncstat = NF_PUT_ATT_TEXT(nc_3d_id, NF_GLOBAL,'title'             &
     &,42, 'Monthly mean output for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF8')

      ncstat = NF_PUT_ATT_TEXT(nc_3d_id, NF_GLOBAL,'history'           &
     &,42, 'Monthly mean output for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(nc_3d_id, NF_GLOBAL,'conventions'       &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF10')

      ncstat = NF_PUT_ATT_TEXT(nc_3d_id, NF_GLOBAL,'source'            &
     &,24, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF11')

      ncstat = NF_PUT_ATT_INT(nc_3d_id, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF11')



!-----------------------------------------------------------------------

      ncdims(1) = timeid

      CALL NETCDF_DEF_VARDB(nc_3d_id,4,'time',1,ncdims,ncvarid,          &
     &   16,'day as %Y%m%d.%f',4,'time',rmasko,31,io_stdo_bgc)

!
! Define variables : grid
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid

      CALL NETCDF_DEF_VARSG(nc_3d_id,8,'scal_lon',2,ncdims,ncvarid,          &
     &   8,'degree E',34,'2-d longitude of scalar grid cells',rmasko,20,io_stdo_bgc)
      icode=1
      ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_3d_id,8,'scal_lat',2,ncdims,ncvarid,          &
     &   8,'degree N',33,'2-d latitude of scalar grid cells',rmasko,21,io_stdo_bgc)
      icode=2
      ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_3d_id,9,'scal_wdep',2,ncdims,ncvarid,          &
     &   5,'meter',37,'2-d water depth at scalar grid points',rmasko,22,io_stdo_bgc)
      icode=3
      ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_3d_id,6,'size_x',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in x-direction',rmasko,23,io_stdo_bgc)
      icode=4
      ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_3d_id,6,'size_y',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in y-direction',rmasko,24,io_stdo_bgc)
      icode=5
      ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid

      CALL NETCDF_DEF_VARSG(nc_3d_id,6,'size_z',3,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in z-direction',rmasko,25,io_stdo_bgc)
      icode=6
      ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

      ncdims(1) = nclevid

      CALL NETCDF_DEF_VARSG(nc_3d_id,5,'depth',1,ncdims,ncvarid,          &
     &   5,'meter',32,'1-d layer depths of ocean tracer',rmasko,26,io_stdo_bgc)


!
! Define variables : total bgc data (js: what means 'total' here?)
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid
      ncdims(4) = timeid




       CALL NETCDF_DEF_VARSG(nc_3d_id,8,'phosph_t',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',19,'Dissolved phosphate',                    &
     &    rmasko,44,io_stdo_bgc)
      icode=11
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,9,'nitrate_t',4,ncdims,ncvarid,    &
     &    9,'kmol/m**3',17,'Dissolved nitrate',                      &
     &    rmasko,44,io_stdo_bgc)
      icode=14
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,8,'silica_t',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',18,'Dissolved silicate',                     &
     &    rmasko,44,io_stdo_bgc)
      icode=15
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

           CALL NETCDF_DEF_VARSG(nc_3d_id,6,'iron_t',4,ncdims,ncvarid,   &
     &    12,'kmol Fe/m**3',18,'Iron concentration',                 &
     &    rmasko,47,io_stdo_bgc)
      icode=31
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

        CALL NETCDF_DEF_VARSG(nc_3d_id,8,'oxygen_t',4,ncdims,ncvarid,    &
     &    9,'kmol/m**3',20,'oxygen concentration',                   &
     &    rmasko,46,io_stdo_bgc)
      icode=12
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,8,'alkali_t',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',16,'total alkalinity',                       &
     &    rmasko,44,io_stdo_bgc)
      icode=10
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,5,'dic_t',4,ncdims,ncvarid,        &
     &    9,'kmol/m**3',26,'dissolved inorganic carbon',             &
     &    rmasko,40,io_stdo_bgc)
      icode=7
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)
#ifdef __c_isotopes
       CALL NETCDF_DEF_VARSG(nc_3d_id,7,'dic13_t',4,ncdims,ncvarid,        &
     &    9,'kmol/m**3',28,'dissolved inorganic carbon13',             &
     &    rmasko,40,io_stdo_bgc)
      icode=8
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,7,'dic14_t',4,ncdims,ncvarid,        &
     &    9,'kmol/m**3',28,'dissolved inorganic carbon14',             &
     &    rmasko,40,io_stdo_bgc)
      icode=9
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)
#endif
       CALL NETCDF_DEF_VARSG(nc_3d_id,5,'doc_t',4,ncdims,ncvarid,        &
     &    9,'kmol/m**3',24,'dissolved organic carbon',               &
     &    rmasko,40,io_stdo_bgc)
      icode=16
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,5,'poc_t',4,ncdims,ncvarid,        &
     &    9,'kmol/m**3',26,'particulate organic carbon',             &
     &    rmasko,40,io_stdo_bgc)
      icode=17
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,6,'calc_t',4,ncdims,ncvarid,       &
     &    9,'kmol/m**3',17,'calcium carbonate',                      &
     &    rmasko,40,io_stdo_bgc)
      icode=24
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,6,'opal_t',4,ncdims,ncvarid,       &
     &    9,'kmol/m**3',15,'biogenic silica',                        &
     &    rmasko,40,io_stdo_bgc)
      icode=27
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

#ifdef ANTC14
       CALL NETCDF_DEF_VARSG(nc_3d_id,6,'ac14_t',4,ncdims,ncvarid,       &
     &    9,'kmol/m**3',17,'anthropogenic C14',                      &
     &    rmasko,42,io_stdo_bgc)
#endif
#ifdef PCFC
       CALL NETCDF_DEF_VARSG(nc_3d_id,7,'cfc11_t',4,ncdims,ncvarid,      &
     &    9,'kmol/m**3',19,'anthropogenic CFC11',                    &
     &    rmasko,43,io_stdo_bgc)
      icode=35
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)

       CALL NETCDF_DEF_VARSG(nc_3d_id,7,'cfc12_t',4,ncdims,ncvarid,      &
     &    9,'kmol/m**3',19,'anthropogenic CFC12',                    &
     &    rmasko,44,io_stdo_bgc)
      icode=36
       ncstat = NF_PUT_ATT_INT(nc_3d_id,ncvarid, 'code', NF_INT, 1, icode)
#endif
!#ifdef AGG     ! js try for output of 'wmass' (settling velocity of aggregation scheme)
!       CALL NETCDF_DEF_VARSG(nc_3d_id,7,'wmass_t',4,ncdims,ncvarid,      &
!     &    5,'m/day',25,'settling velocity of mass',                    &
!     &    rmasko,43,io_stdo_bgc)
!#endif
!
! END Define variables
!
!-----------------------------------------------------------------------

      ncstat = NF_ENDDEF(nc_3d_id)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF00')
!
! Set fill mode
!
!-----------------------------------------------------------------------

      ncstat = NF_SET_FILL(nc_3d_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF97')

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
      CALL write_netcdf_var(nc_3d_id,'scal_lon',zfield(1,1),1,0)

      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgiph(ii,jj)*aradtogra
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_3d_id,'scal_lat', zfield(1,1),1,0)

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
      CALL write_netcdf_var(nc_3d_id,'scal_wdep',zfield(1,1),1,0)

      IF(p_pe==p_io) THEN
        ncstat = NF_INQ_VARID(nc_3d_id,'depth',ncvarid )
        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF102a')
        ncstat = NF_PUT_VAR_DOUBLE (nc_3d_id,ncvarid,ptiestu(1) )
        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF103a')
      ENDIF

      CALL write_netcdf_var(nc_3d_id,'size_x',pdlxp(1,1),1,0)
      CALL write_netcdf_var(nc_3d_id,'size_y',pdlyp(1,1),1,0)
      CALL write_netcdf_var(nc_3d_id,'size_z',pddpo(1,1,1),kpke,0)

!
! Close File
!
!-----------------------------------------------------------------------
!      IF(p_pe==p_io) THEN
!        ncstat = NF_CLOSE(nc_3d_id)
!        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_3D: Problem with netCDF200')
!      ENDIF
!

      RETURN
      END
