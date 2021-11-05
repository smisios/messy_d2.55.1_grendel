      SUBROUTINE OPEN_BGCMEAN_SED(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce)

!****************************************************************
!
!**** *WRITE_BGCMEAN_2D* - calculate and write 3-dimensional bgc mean data.
!
!     Joachim Segschneider  MPI-Met, HH*    03.08.05
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Open bgcmean_sed.nc and write the header
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *OPEN_BGCMEAN_SED(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
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
      USE mo_parallel
 
      implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k,ii,jj
      INTEGER :: kplyear,kplmon,kplday,kpldtoce

      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: pgila(kpie*2,kpje*2)
      REAL :: pgiph(kpie*2,kpje*2)
      REAL :: ptiestu(kpke)


      INCLUDE 'netcdf.inc'
      INTEGER :: ncvarid,ncstat,ncoldmod,ncdims(4)
      INTEGER :: nclatid,nclonid,nclevid
      INTEGER :: timeid,nstart2(2),ncount2(2),nstride2(2)

      INTEGER :: idate(5)
      INTEGER :: chunksize,initialsize

      REAL :: zfield(kpie,kpje)

      REAL :: pi,rad,radi

!-----------------------------------------------------------------------

      pi        = 4.*ATAN(1.)
      rad       = pi/180.
      radi      = 1./rad


      idate(1) = kplyear
      idate(2) = kplmon
      idate(3) = kplday
      idate(4) = kpldtoce
      idate(5) = ldtbgc

      IF (p_pe==p_io) THEN
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Creating sed bgcmean data at date :'   &
     &  ,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_bgc,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_bgc,*) 'Bgc   model step number is ',idate(5)


!
! Open netCDF data file
!
!-----------------------------------------------------------------------

      ncstat = NF_CREATE('bgcmean_sed.nc',NF_CLOBBER, nc_sed_id)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF1')

      chunksize = 1024*1024*32 ! 32 MB --> man 3f netcdf
      initialsize = 1024*1024*32 ! 32 MB --> man 3f netcdf
      
!      ncstat = nf__create('bgcmean_3d.nc',NF_CLOBBER,initialsize,chunksize,nc_sed_id)
!      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF1')
!
! Define dimension
!
!-----------------------------------------------------------------------

      ncstat = NF_DEF_DIM(nc_sed_id, 'lon', ie_g, nclonid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF2')

      ncstat = NF_DEF_DIM(nc_sed_id, 'lat', je_g, nclatid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF3')

      ncstat = NF_DEF_DIM(nc_sed_id, 'depth', ks, nclevid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF4')

      ncstat = NF_DEF_DIM(nc_sed_id, 'time',NF_UNLIMITED, timeid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF7'      )
!
! Define global attributes
!
!-----------------------------------------------------------------------
 
      ncstat = NF_PUT_ATT_TEXT(nc_sed_id, NF_GLOBAL,'title'             &
     &,40, 'annual mean sediment data from HAMOCC5.1')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF8')

      ncstat = NF_PUT_ATT_TEXT(nc_sed_id, NF_GLOBAL,'history'           &
     &,40, 'annual mean sediment data from HAMOCC5.1')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(nc_sed_id, NF_GLOBAL,'conventions'       &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF10')

      ncstat = NF_PUT_ATT_TEXT(nc_sed_id, NF_GLOBAL,'source'            &
     &,16, 'HAMOCC5.1 MPI-OM')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF11')

      ncstat = NF_PUT_ATT_INT(nc_sed_id, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF11')

!
! Define variables : grid
!



!-----------------------------------------------------------------------

      ncdims(1) = timeid

      CALL NETCDF_DEF_VARSG(nc_sed_id,4,'time',1,ncdims,ncvarid,          &
     &   16,'day as %Y%m%d.%f',4,'time',rmasko,31,io_stdo_bgc)

!



!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid

      CALL NETCDF_DEF_VARSG(nc_sed_id,8,'scal_lon',2,ncdims,ncvarid,          &
     &   8,'degree E',34,'2-d longitude of scalar grid cells',rmasko,20,io_stdo_bgc) 

      CALL NETCDF_DEF_VARSG(nc_sed_id,8,'scal_lat',2,ncdims,ncvarid,          &
     &   8,'degree N',33,'2-d latitude of scalar grid cells',rmasko,21,io_stdo_bgc) 

      CALL NETCDF_DEF_VARSG(nc_sed_id,9,'scal_wdep',2,ncdims,ncvarid,          &
     &   5,'meter',37,'2-d water depth at scalar grid points',rmasko,22,io_stdo_bgc) 

      CALL NETCDF_DEF_VARSG(nc_sed_id,6,'size_x',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in x-direction',rmasko,23,io_stdo_bgc) 

      CALL NETCDF_DEF_VARSG(nc_sed_id,6,'size_y',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in y-direction',rmasko,24,io_stdo_bgc) 


      ncdims(1) = nclevid

!     CALL NETCDF_DEF_VARSG(nc_sed_id,6,'size_z',1,ncdims,ncvarid,          &
!    &   5,'meter',24,'depth of sediment layers',rmasko,25,io_stdo_bgc)  
          
      
      ncdims(1) = nclevid

      CALL NETCDF_DEF_VARSG(nc_sed_id,5,'depth',1,ncdims,ncvarid,          &
     &   5,'meter',32,'1-d depths of ocean model layers',rmasko,26,io_stdo_bgc)  
     

!
! Define variables : total bgc data
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid
      ncdims(4) = timeid


       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'powaic',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',27,'Pore water inorganic carbon',               &
     &    rmasko,44,io_stdo_bgc)

       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'powaal',4,ncdims,ncvarid,    &
     &    9,'kmol/m**3',21,'Pore water alkalinity',                    &
     &    rmasko,44,io_stdo_bgc)
     
       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'powaph',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',20,'pore water phosphate',                     &
     &    rmasko,44,io_stdo_bgc)

       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'powaox',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',17,'pore water oxygen',                     &
     &    rmasko,44,io_stdo_bgc)

    
       CALL NETCDF_DEF_VARSG(nc_sed_id,5,'pown2',4,ncdims,ncvarid,   &  
     &    12,'kmol Fe/m**3',19,'pore water nitrogen',                 &  
     &    rmasko,47,io_stdo_bgc)
     
       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'powno3',4,ncdims,ncvarid,    & 
     &    9,'kmol/m**3',18,'pore water nitrate',                        &
     &    rmasko,46,io_stdo_bgc)         
  
       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'powasi',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',19,'pore water silicate',                       &
     &    rmasko,44,io_stdo_bgc)
     
       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'ssso12',4,ncdims,ncvarid,        & 
     &    9,'kmol/m**3',29,'solid sediment organic carbon',            &
     &    rmasko,40,io_stdo_bgc)   
     
       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'ssssil',4,ncdims,ncvarid,        & 
     &    9,'kmol/m**3',19,'solid sediment opal',               &
     &    rmasko,40,io_stdo_bgc)       
     
       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'sssc12',4,ncdims,ncvarid,        & 
     &    9,'kmol/m**3',32,'solid sediment calcium carbonate',             &
     &    rmasko,40,io_stdo_bgc)        
     
       CALL NETCDF_DEF_VARSG(nc_sed_id,6,'ssster',4,ncdims,ncvarid,       & 
     &    9,'kmol/m**3',20,'terrigenous sediment',                       &
     &    rmasko,40,io_stdo_bgc)        
     
!
! END Define variables
!
!-----------------------------------------------------------------------

      ncstat = NF_ENDDEF(nc_sed_id)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF00')
!
! Set fill mode
!
!-----------------------------------------------------------------------

      ncstat = NF_SET_FILL(nc_sed_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF97')

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
            zfield(i,j) = pgila(ii,jj)*radi
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_sed_id,'scal_lon',zfield(1,1),1,0)
               
      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgiph(ii,jj)*radi
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_sed_id,'scal_lat', zfield(1,1),1,0)

      DO j=1,kpje
         DO i=1,kpie
            zfield(i,j) = 0.0
         ENDDO
      ENDDO
      DO k=1,kpke
         DO j=1,kpje
            DO i=1,kpie
               zfield(i,j) = zfield(i,j) + pddpo(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_sed_id,'scal_wdep',zfield(1,1),1,0)

      IF(p_pe==p_io) THEN 
        ncstat = NF_INQ_VARID(nc_sed_id,'depth',ncvarid )
        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF102a')
!       ncstat = NF_PUT_VAR_DOUBLE (nc_sed_id,ncvarid,ptiestu(1) )
        ncstat = NF_PUT_VAR_DOUBLE (nc_sed_id,ncvarid,z_sed(1) )
        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF103a')
      ENDIF

      CALL write_netcdf_var(nc_sed_id,'size_x',pdlxp(1,1),1,0)
      CALL write_netcdf_var(nc_sed_id,'size_y',pdlyp(1,1),1,0)
!
! Close File
!
!-----------------------------------------------------------------------
!      IF(p_pe==p_io) THEN
!        ncstat = NF_CLOSE(nc_sed_id)
!        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_SED: Problem with netCDF200')
!      ENDIF
!

      RETURN
      END
