      SUBROUTINE OPEN_BGCMEAN_2D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce)

!****************************************************************
!
!**** *OPEN_BGCMEAN_2D* - calculate and write 2-dimensional bgc mean data.
!
!     Patrick Wetzel,    *MPI-Met, HH*    24.05.04
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Write bgc mean data.
!     js some fields are defined in open_bgcmean_2d that are not picked up
!        in write_bgcmean. these fields then have no time dimension in bgcmean_2d.nc
!        co214flux, co2flux, n2oflux, prorca, prcaca, silpro, produs
!        -PANTC14  rantc14
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *OPEN_BGCMEAN_2D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
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
      INTEGER :: nclatid,nclonid
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

      IF(p_pe == p_io) THEN

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Create 2D bgcmean data at date :'    &
     &  ,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_bgc,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_bgc,*) 'Bgc   model step number is ',idate(5)

!      WRITE(io_stdo_bgc,*) 'OPEN_BGCMEAN_2D: ',kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
!     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce

!
! Open netCDF data file js: here some date information should be added to name of bgcmean
!                           e.g. from kplyear
!
!-----------------------------------------------------------------------

!      ncstat = NF_CREATE('bgcmean_2d.nc',NF_CLOBBER, nc_2d_id)
      ncstat = NF_CREATE('bgcmean_2d.nc',IOR(NF_CLOBBER,NF_64BIT_OFFSET), nc_2d_id)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF1')

      chunksize   = 1024*1024*32 ! 32 MB --> man 3f netcdf
      initialsize = 1024*1024*32 ! 32 MB --> man 3f netcdf

!      ncstat = nf__create('bgcmean_2d.nc',NF_CLOBBER,initialsize,chunksize,nc_2d_id)
!      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF1')

!
! Define dimension
!
!-----------------------------------------------------------------------

      ncstat = NF_DEF_DIM(nc_2d_id, 'lon', ie_g, nclonid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF2')

      ncstat = NF_DEF_DIM(nc_2d_id, 'lat', je_g, nclatid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF3')

      ncstat = NF_DEF_DIM(nc_2d_id, 'time',NF_UNLIMITED, timeid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF7'      )

!
! Define global attributes
!
!-----------------------------------------------------------------------

      ncstat = NF_PUT_ATT_TEXT(nc_2d_id, NF_GLOBAL,'title'             &
     &,42, 'Monthly mean output for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF8')

      ncstat = NF_PUT_ATT_TEXT(nc_2d_id, NF_GLOBAL,'history'           &
     &,42, 'Monthly mean output for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(nc_2d_id, NF_GLOBAL,'conventions'       &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF10')

      ncstat = NF_PUT_ATT_TEXT(nc_2d_id, NF_GLOBAL,'source'            &
     &,24, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF11')

      ncstat = NF_PUT_ATT_INT(nc_2d_id, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF11')

!
! Define variables : grid
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid

      CALL NETCDF_DEF_VARSG(nc_2d_id,8,'scal_lon',2,ncdims,ncvarid,          &
     &   8,'degree E',34,'2-d longitude of scalar grid cells',rmasko,20,io_stdo_bgc)
      icode=1
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,8,'scal_lat',2,ncdims,ncvarid,          &
     &   8,'degree N',33,'2-d latitude of scalar grid cells',rmasko,21,io_stdo_bgc)
      icode=2
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,9,'scal_wdep',2,ncdims,ncvarid,          &
     &   5,'meter',37,'2-d water depth at scalar grid points',rmasko,22,io_stdo_bgc)
      icode=3
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'size_x',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in x-direction',rmasko,23,io_stdo_bgc)
      icode=4
     ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'size_y',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in y-direction',rmasko,24,io_stdo_bgc)
      icode=5
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)
!
!
! Define variables : 1-D mean data
!
!-----------------------------------------------------------------------

      ncdims(1) = timeid

!      CALL NETCDF_DEF_VARSG(nc_2d_id,9,'steps_p_m',1,ncdims,ncvarid,          &
!     &   5,'steps',25,'model timesteps per month',rmasko,30,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(nc_2d_id,4,'time',1,ncdims,ncvarid,          &
     &   16,'day as %Y%m%d.%f',4,'time',rmasko,31,io_stdo_bgc)

!
! Define variables : 2-D bgc mean data
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = timeid

      CALL NETCDF_DEF_VARSG(nc_2d_id,16,'co2fluxdown_mean',3,ncdims,ncvarid,   &
     &   9,'kmol/m**2',13, 'CO2 Flux Down',rmasko,50,io_stdo_bgc)
      icode=65
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,14,'co2fluxup_mean',3,ncdims,ncvarid,     &
     &   9,'kmol/m**2',11, 'CO2 Flux Up',rmasko,51,io_stdo_bgc)
      icode=66
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,4,'pco2',3,ncdims,ncvarid,                &
     &   3,'ppm',20, 'CO2 partial-pressure',rmasko,52,io_stdo_bgc)
      icode=67
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,3,'dms',3,ncdims,ncvarid,                 &
     &   9,'kmol/m**2',15, 'DiMethylSulfate',rmasko,53,io_stdo_bgc)
      icode=29
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,7,'dmsflux',3,ncdims,ncvarid,             &
     &   9,'kmol/m**2',8, 'DMS Flux',rmasko,54,io_stdo_bgc)
      icode=68
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,7,'dmsprod',3,ncdims,ncvarid,             &
     &   9,'kmol/m**2',14, 'DMS Production',rmasko,55,io_stdo_bgc)
      icode=69
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,7,'dms_bac',3,ncdims,ncvarid,             &
     &   9,'kmol/m**2',25, 'DMS bacterial consumption',rmasko,56,io_stdo_bgc)
      icode=70
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'dms_uv',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',24, 'DMS uv-light destruction',rmasko,57,io_stdo_bgc)
      icode=71
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'oxflux',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',11, 'oxygen flux',rmasko,61,io_stdo_bgc)
      icode=72
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,5,'kwco2',3,ncdims,ncvarid,               &
     &   19,'kmol/(m**2*sec*ppm)',24, 'co2 exchange coefficient',rmasko,62,io_stdo_bgc)
      icode=73
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'niflux',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',13, 'nitrogen flux',rmasko,63,io_stdo_bgc)
      icode=74
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'atmco2',3,ncdims,ncvarid,              &
     &   3,'ppm',15, 'Atmospheric CO2',rmasko,64,io_stdo_bgc)
      icode=61
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,5,'atmo2',3,ncdims,ncvarid,               &
     &   3,'ppm',14, 'Atmospheric O2',rmasko,65,io_stdo_bgc)
      icode=62
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,5,'atmn2',3,ncdims,ncvarid,               &
     &   3,'ppm',14, 'Atmospheric N2',rmasko,66,io_stdo_bgc)
      icode=64
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'opex90',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',17, 'opal flux in 90 m',rmasko,41,io_stdo_bgc)
      icode=75
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,8,'opex1000',3,ncdims,ncvarid,            &
     &   9,'kmol/m**2',19, 'opal flux in 1000 m',rmasko,42,io_stdo_bgc)
      icode=76
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,8,'opex2000',3,ncdims,ncvarid,            &
     &   9,'kmol/m**2',19, 'opal flux in 2000 m',rmasko,43,io_stdo_bgc)
      icode=77
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'caex90',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',17, 'calc flux in 90 m',rmasko,44,io_stdo_bgc)
      icode=78
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,8,'caex1000',3,ncdims,ncvarid,            &
     &   9,'kmol/m**2',19, 'calc flux in 1000 m',rmasko,45,io_stdo_bgc)
      icode=79
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,8,'caex2000',3,ncdims,ncvarid,            &
     &   9,'kmol/m**2',19, 'calc flux in 2000 m',rmasko,46,io_stdo_bgc)
      icode=80
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'coex90',3,ncdims,ncvarid,              &
     &   19,'kmolP/m**2/timestep',17, 'phos flux in 90 m',rmasko,47,io_stdo_bgc)
      icode=81
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,8,'coex1000',3,ncdims,ncvarid,            &
     &   9,'kmol/m**2',19, 'phos flux in 1000 m',rmasko,48,io_stdo_bgc)
      icode=82
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,8,'coex2000',3,ncdims,ncvarid,            &
     &   9,'kmol/m**2',19, 'phos flux in 2000 m',rmasko,49,io_stdo_bgc)
      icode=83
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'export',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',26,'detritus export production',rmasko,47,io_stdo_bgc)
      icode=84
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'expoca',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',25,'calcium export production',rmasko,48,io_stdo_bgc)
      icode=85
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'exposi',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',22,'opal export production',rmasko,49,io_stdo_bgc)
      icode=86
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)
#ifdef ANTC14
      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'ac14fx',3,ncdims,ncvarid,              &
     &   9,'kmol/m**2',22, 'anthropogenic c14 flux',rmasko,49,io_stdo_bgc)
#endif

#ifdef PCFC
      CALL NETCDF_DEF_VARSG(nc_2d_id,7,'cfc11fx',3,ncdims,ncvarid,             &
     &   9,'kmol/m**2',10, 'cfc11 flux',rmasko,50,io_stdo_bgc)
      icode=87
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,7,'cfc12fx',3,ncdims,ncvarid,             &
     &   9,'kmol/m**2',10, 'cfc12 flux',rmasko,51,io_stdo_bgc)
      icode=88
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'pcfc11',3,ncdims,ncvarid,             &
     &   3,'ppm',30, 'cfc11 atmos. partial preassure',rmasko,51,io_stdo_bgc)
      icode=89
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'pcfc12',3,ncdims,ncvarid,             &
     &   3,'ppm',30, 'cfc12 atmos. partial preassure',rmasko,52,io_stdo_bgc)
      icode=90
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)
#endif

      CALL NETCDF_DEF_VARSG(nc_2d_id,9,'co214flux',3,ncdims,ncvarid,       &
     &   9,'kmol/m**2',10, 'CO214 Flux',rmasko,40,io_stdo_bgc)
      icode=91
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,7,'co2flux',3,ncdims,ncvarid,         &
     &   9,'kmol/m**2',8, 'CO2 Flux',rmasko,41,io_stdo_bgc)
      icode=92
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,7,'n2oflux',3,ncdims,ncvarid,    &
     &   9,'kmol/m**2',8, 'N2O Flux',rmasko,42,io_stdo_bgc)
      icode=93
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'prorca',3,ncdims,ncvarid,          &
     &   9,'kmol/m**2',22, 'detritus sediment flux',rmasko,43,io_stdo_bgc)
      icode=94
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'prcaca',3,ncdims,ncvarid,          &
     &   9,'kmol/m**2',19, 'CaCO3 sediment flux',rmasko,44,io_stdo_bgc)
      icode=95
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'silpro',3,ncdims,ncvarid,          &
     &   9,'kmol/m**2',18, 'opal sediment flux',rmasko,45,io_stdo_bgc)
      icode=96
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_2d_id,6,'produs',3,ncdims,ncvarid,          &
     &   7,'kg/m**2',41, 'dust flux to sediment (free + aggregated)',rmasko,46,io_stdo_bgc)
      icode=97
      ncstat = NF_PUT_ATT_INT(nc_2d_id,ncvarid, 'code', NF_INT, 1, icode)

#ifdef ANTC14
      CALL NETCDF_DEF_VARSG(nc_2d_id,7,'rantc14',3,ncdims,ncvarid,          &
     &   5,'ratio',16, 'atmos. C14 ratio',rmasko,46,io_stdo_bgc)
#endif


!
! END Define variables
!
!-----------------------------------------------------------------------

      ncstat = NF_ENDDEF(nc_2d_id)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF00')
!
! Set fill mode
!
!-----------------------------------------------------------------------

      ncstat = NF_SET_FILL(nc_2d_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF97')

      ENDIF ! p_pe == p_io

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
      CALL write_netcdf_var(nc_2d_id,'scal_lon',zfield(1,1),1,0)

      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgiph(ii,jj)*aradtogra
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_2d_id,'scal_lat', zfield(1,1),1,0)

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
!js output of wdep was missing in 2d file
      CALL write_netcdf_var(nc_2d_id,'scal_wdep',zfield(1,1),1,0)

      CALL write_netcdf_var(nc_2d_id,'size_x',pdlxp(1,1),1,0)
      CALL write_netcdf_var(nc_2d_id,'size_y',pdlyp(1,1),1,0)

! Close File
!

!      if (p_pe==p_io) then
!      ncstat = NF_CLOSE(nc_2d_id)
!      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_BGCMEAN_2D: Problem with netCDF200')
!      end if


      RETURN
      END
