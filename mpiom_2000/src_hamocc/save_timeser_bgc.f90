      SUBROUTINE SAVE_TIMESER_BGC
!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/save_timeser_bgc.f90,v $\\
!$Revision: 1.4.2.1.4.1.2.2.4.1.2.2.2.3.2.2 $\\
!$Date: 2006/04/06 10:01:07 $\\
!$Name: mpiom_1_2_0 $\\
!
!****************************************************************
!
!**** *SAVE_TIMSER_BGC*
!
!     S.Legutke,        *MPI-MaD, HH*    03.08.01
!     P.Wetzel          *MPI-Met, HH*    31.01.05
!     - write NetCDF
!
!     Purpose
!     -------
!     - save bgc time series.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *SAVE_TIMESER_BGC*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
!ik added trap depths for time series
!ik changed output to ASCII formatted file
!ik added budget output to be written at the end of bgcout (being too lazy
!ik to write a seperate routine)

      USE mo_timeser_bgc
      USE mo_control_bgc
      USE mo_biomod
      USE mo_carbch
      use mo_param1_bgc

      use mo_mpi
      use mo_parallel

implicit none

      INTEGER :: i,n
      REAL(wp) :: ts1_time(lents1)
      REAL(wp) :: ts1recv(nvarts1,nelets1,lents1)
      CHARACTER*33 ystring

#ifdef PNETCDF
      INCLUDE 'netcdf.inc'

      INTEGER :: ncid,ncstat,ncvarid
      INTEGER :: nvarts1id, nelets1id, nctimeid
      INTEGER :: ncdims(3)

      WRITE(io_stdo_bgc,*)'Save Timeser BGC in NetCDF Format'
#else
      INTEGER :: l,l1,l2,l3
#endif


! For a parallel run, sum up ts1 on p_io
! Note that this code also works for nonparallel runs,
! it just does nothing in this case!

      ! All non-I/O-PEs just send to p_io
      IF(p_pe/=p_io .AND. p_pe<nprocxy) CALL p_send(ts1,p_io,123)

      ! The I/O-PE sums up
      IF(p_pe==p_io) THEN
        DO n=0,nprocxy-1
          IF(n/=p_io) THEN
            CALL p_recv(ts1recv,n,123)
            ts1(:,:,:) = ts1(:,:,:) + ts1recv(:,:,:)
          ENDIF
        ENDDO
      ENDIF

      IF(p_pe==p_io) THEN ! Output is only done by I/O-PE


! divide by depth
!
!      DO l=1,nts
!      DO k=1,lents1
!         ts1(itsphosph,l+1,k) = ts1(itsphosph,l+1,k)
!         ts1(itsopal,l+1,k)   = ts1(itsopal,l+1,k)
!         ts1(itssilica,l+1,k) = ts1(itssilica,l+1,k)
!         ts1(itsphy,l+1,k)    = ts1(itsphy,l+1,k)
!         ts1(itsdet,l+1,k)    = ts1(itsdet,l+1,k)
!#ifdef AGG
!         ts1(itsnos,l+1,k)    = ts1(itsnos,l+1,k)
!#endif
!         ts1(itsphosy,l+1,k)  = ts1(itsphosy,l+1,k)
!         ts1(itszoo,l+1,k)    = ts1(itszoo,l+1,k)
!         ts1(itssco212,l+1,k) = ts1(itssco212,l+1,k)
!         ts1(itsdoc,l+1,k)    = ts1(itsdoc,l+1,k)
!         ts1(itsiron,l+1,k)   = ts1(itsiron,l+1,k)
!      ENDDO
!      ENDDO


#ifdef PNETCDF

!
! Calculate time
!

      ystring(1:33)='seconds since 0001-01-01 23:59:59'

      if (bgcstartyear .ge. 1000) then
        WRITE(ystring(15:18),'(I4)') bgcstartyear
      elseif (bgcstartyear .ge. 100) then
        WRITE(ystring(16:18),'(I3)') bgcstartyear
      elseif (bgcstartyear .ge. 10) then
        WRITE(ystring(17:18),'(I2)') bgcstartyear
      else
        WRITE(ystring(18:18),'(I1)') bgcstartyear
      endif

      if (bgcstartmonth .ge. 10) then
        WRITE(ystring(20:21),'(I2)') bgcstartmonth
      else
        WRITE(ystring(21:21),'(I1)') bgcstartmonth
      endif

      if (bgcstartday .ge. 10) then
        WRITE(ystring(23:24),'(I2)') bgcstartday
      else
        WRITE(ystring(24:24),'(I1)') bgcstartday
      endif

      WRITE(io_stdo_bgc,*) 'Write timeser_bgc in ', ystring

      do i=1,lents1
         ts1_time(i) = dtbgc * REAL(nfreqts1, wp) / 2._wp &
              + dtbgc * REAL(nfreqts1*(i-1), wp)
!      WRITE(io_stdo_bgc,*)'TEST:', ts1_time(i)/86400
      enddo



!
! Open netCDF data file
!
!      ncstat = NF_CREATE('timeser_bgc.nc',NF_CLOBBER, ncid)
      ncstat = NF_CREATE('timeser_bgc.nc',IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('SAVE_TIMESER_BGC: Problem with netCDF1')

      ! Define global attributes
      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'title'             &
     &,21, 'marine bgc timeseries')
      IF ( ncstat .NE. NF_NOERR ) call stop_all ('SAVE_TIMESER_BGC: Problem with netCDF2')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'source'             &
     &,7, 'HAMOCC5')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF3')

!      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'author'             &
!     &,42, 'Patrick Wetzel, <patrick.wetzel@dkrz.de>')
!      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF4')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'institution'             &
     &,45, 'Max-Planck-Institute for Meteorology, Hamburg')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF5')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'Conventions'             &
     &,12, 'non-standard')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('write_all_to_ncdf.f90: Problem with netCDF6')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'History'             &
     &,7, 'created')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('write_all_to_ncdf.f90: Problem with netCDF7')

      WRITE(io_stdo_bgc,*)'Global attributes of Timeser BGC'
!
! Define dimension
! ---------------------------------------------------------------------------------
      ncstat = NF_DEF_DIM(ncid, 'dim_nvarts1',nvarts1,nvarts1id)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF11')

      ncstat = NF_DEF_DIM(ncid, 'dim_nelets1',nelets1,nelets1id)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF12')

!      ncstat = NF_DEF_DIM(ncid, 'dim_lents1',lents1,lents1id)
!      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCD13')


      ncstat = NF_DEF_DIM(ncid, 'time',lents1,nctimeid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF14')

      WRITE(io_stdo_bgc,*)'Define dimensions of Timeser BGC'


! Define grid variables
! ---------------------------------------------------------------------------------


      ncdims(1) = nctimeid

      ncstat = NF_DEF_VAR(ncid,'time',NF_REAL,1,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF21')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name',4, 'time')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF22')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',33,ystring)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF22a')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'calendar',9, 'gregorian')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF22b')


      WRITE(io_stdo_bgc,*)'Define time variable of Timeser BGC'

!               time:long_name = "time" ;
!               time:units = "day since 1999-10-01 00:00:00" ;
!               time:calendar = "gregorian" ;

! Define data variables
! ---------------------------------------------------------------------------------

      ncdims(1) = nvarts1id
      ncdims(2) = nelets1id
      ncdims(3) = nctimeid


      ncstat = NF_DEF_VAR(ncid,'timeseries',NF_REAL,3,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF31')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',4, 'none')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF32')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                  &
     &,15, 'timeseries data')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF33')

      ncstat = NF_PUT_ATT_REAL(ncid,ncvarid,'missing_value',NF_FLOAT,1,99999.)

      WRITE(io_stdo_bgc,*)'Define timeseries variable of Timeser BGC'

! Define data variables
! ---------------------------------------------------------------------------------

      ncdims(1) = nelets1id
      ncdims(2) = nctimeid

      CALL NETCDF_DEF_VARSG(ncid,3,'dic',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',26, 'dissolved inorganic carbon',rmasko,40,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 200)

      CALL NETCDF_DEF_VARSG(ncid,5,'phosy',2,ncdims,ncvarid,                    &
     &   14,'10^6 kmolP/day',14, 'photosynthesis',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 201)

      CALL NETCDF_DEF_VARSG(ncid,6,'phosph',2,ncdims,ncvarid,                   &
     &   9,'kmol/m**3',19,'dissolved phosphate',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 202)

      CALL NETCDF_DEF_VARSG(ncid,6,'oxygen',2,ncdims,ncvarid,                   &
     &   9,'kmol/m**3',16,'dissolved oxygen',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 203)

      CALL NETCDF_DEF_VARSG(ncid,6,'gasnit',2,ncdims,ncvarid,                   &
     &   9,'kmol/m**3',21,'gaseous nitrogen (N2)',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 204)

      CALL NETCDF_DEF_VARSG(ncid,4,'ano3',2,ncdims,ncvarid,                     &
     &   9,'kmol/m**3',17,'dissolved nitrate',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 205)

      CALL NETCDF_DEF_VARSG(ncid,6,'silica',2,ncdims,ncvarid,                   &
     &   9,'kmol/m**3',22,'silicid acid (Si(OH)4)',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 206)

      CALL NETCDF_DEF_VARSG(ncid,3,'doc',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',24,'dissolved organic carbon',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 207)

      CALL NETCDF_DEF_VARSG(ncid,3,'phy',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',27,'phytoplankton concentration',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 208)

      CALL NETCDF_DEF_VARSG(ncid,3,'zoo',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',25,'Zooplankton concentration',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 209)

      CALL NETCDF_DEF_VARSG(ncid,3,'poc',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',25,'particulate organic carbon',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 210)

      CALL NETCDF_DEF_VARSG(ncid,3,'cal',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',17,'Calcium carbonate',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 211)

      CALL NETCDF_DEF_VARSG(ncid,4,'opal',2,ncdims,ncvarid,                     &
     &   9,'kmol/m**3',15,'biogenic silica',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 212)

      CALL NETCDF_DEF_VARSG(ncid,4,'iron',2,ncdims,ncvarid,                     &
     &   9,'kmol/m**3',14,'dissolved iron',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 213)

      CALL NETCDF_DEF_VARSG(ncid,4,'pco2',2,ncdims,ncvarid,                     &
     &   3,'ppm',20,'CO2 partial pressure',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 214)

      CALL NETCDF_DEF_VARSG(ncid,3,'atm',2,ncdims,ncvarid,                      &
     &   3,'ppm',32,'atmospheric CO2 partial pressure',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 215)

      CALL NETCDF_DEF_VARSG(ncid,4,'co2f',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',19,'sea-to-air CO2 flux',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 216)

#ifdef PANTHROPOCO2
      CALL NETCDF_DEF_VARSG(ncid,7,'antpco2',2,ncdims,ncvarid,                     &
     &   9,'kmol/m**2',25,'ant. CO2 partial pressure',rmasko,41,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(ncid,6,'antatm',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',33,'ant. atmospheric partial pressure',rmasko,41,io_stdo_bgc)

      CALL NETCDF_DEF_VARSG(ncid,7,'antco2f',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',24,'ant. sea-to-air CO2 flux',rmasko,41,io_stdo_bgc)
#endif

#ifdef AGG
     CALL NETCDF_DEF_VARSG(ncid,3,'nos',2,ncdims,ncvarid,                       &
     &   9,'kmol/m**2',20,'number of aggregates',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 217)
#endif


      CALL NETCDF_DEF_VARSG(ncid,5,'fdet1',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',24,'detritus flux in depth 1',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 218)

      CALL NETCDF_DEF_VARSG(ncid,5,'fopa1',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',17,'opal flux depth 1',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 219)

      CALL NETCDF_DEF_VARSG(ncid,5,'fcal1',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',30,'calcium carbonate flux depth 1',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 220)

      CALL NETCDF_DEF_VARSG(ncid,5,'fdet2',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',24,'detritus flux in depth 2',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 221)

      CALL NETCDF_DEF_VARSG(ncid,5,'fopa2',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',17,'opal flux depth 2',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 222)

      CALL NETCDF_DEF_VARSG(ncid,5,'fcal2',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',30,'calcium carbonate flux depth 2',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 223)

      CALL NETCDF_DEF_VARSG(ncid,5,'fdet3',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',24,'detritus flux in depth 3',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 224)

      CALL NETCDF_DEF_VARSG(ncid,5,'fopa3',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',17,'opal flux depth 3',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 225)

      CALL NETCDF_DEF_VARSG(ncid,5,'fcal3',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**2',30,'calcium carbonate flux depth 3',rmasko,41,io_stdo_bgc)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 226)

! ---------------------------------------------------------------------------------
! END Define variables
      ncstat = NF_ENDDEF(ncid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF40')

      WRITE(io_stdo_bgc,*)'end define of Timeser BGC'

!
! Set fill mode
! ---------------------------------------------------------------------------------
      ncstat = NF_SET_FILL(ncid,NF_NOFILL,ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF41')

      WRITE(io_stdo_bgc,*)'fill mode of Timeser BGC'

! Save grid variables
! ---------------------------------------------------------------------------------

      ncstat = NF_INQ_VARID(ncid,'time',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF42')

      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1_time)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF43')

      WRITE(io_stdo_bgc,*)'fill time of Timeser BGC'


! Save timeseries data
!---------------------------------------------------------------------
      ncstat = NF_INQ_VARID(ncid,'timeseries',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')

      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      WRITE(io_stdo_bgc,*)'fill timeseries of Timeser BGC'

!
! Write timeser data
!
!-----------------------------------------------------------------------
!      nk=1
!      CALL write_netcdf_var(ncid,'dic',ts1(itssco212,:,:),nk,lents1)
!      CALL write_netcdf_var(ncid,'phosy',ts1(itsphosy,:,:),nk,lents1)

      ncstat = NF_INQ_VARID(ncid,'dic',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itssco212,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'phosy',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsphosy,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'phosph',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsphosph,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'oxygen',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsoxygen,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'gasnit',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsgasnit,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'ano3',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsano3,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'silica',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itssilica,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'doc',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsdoc,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'phy',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsphy,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'zoo',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itszoo,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'poc',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsdet,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'cal',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itscalc,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'opal',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsopal,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'iron',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsiron,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'pco2',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itspco2,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'atm',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsatm,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'co2f',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsco2f,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

#ifdef PANTHROPOCO2
      ncstat = NF_INQ_VARID(ncid,'antpco2',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itspco2a,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'antatm',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsatma,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'antco2f',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsco2fa,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')
#endif

#ifdef AGG
      ncstat = NF_INQ_VARID(ncid,'nos',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsnos,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')
#endif

      ncstat = NF_INQ_VARID(ncid,'fdet1',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its1fdet,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'fopa1',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its1fopa,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'fcal1',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its1fcal,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')


      ncstat = NF_INQ_VARID(ncid,'fdet2',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its2fdet,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'fopa2',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its2fopa,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'fcal2',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its2fcal,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')


      ncstat = NF_INQ_VARID(ncid,'fdet3',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its3fdet,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'fopa3',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its3fopa,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')

      ncstat = NF_INQ_VARID(ncid,'fcal3',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(its3fcal,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF45')




!
! Close File

      ncstat = NF_CLOSE(ncid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF50')

      WRITE(io_stdo_bgc,*)'Timeseries BGC written to "timeser_bgc.nc"'


#else
      IF( lents1 .LT. lts1 ) THEN
         WRITE(io_stdo_bgc,*)                                          &
     &   'Warning : too many samples in time series 1!!!'
      ELSE
         io_timeser_bgc = io_stdi_bgc
         WRITE(io_stdo_bgc,*)'Opening timeser_bgc (',io_timeser_bgc,').'

         CLOSE(io_timeser_bgc)
         OPEN(io_timeser_bgc,FILE='timeser_bgc',STATUS='UNKNOWN'       &
     &       ,FORM='FORMATTED')
         WRITE(io_timeser_bgc,*) 'time steps, number of variables,     &
     &                            number of stations '
         WRITE(io_timeser_bgc,*) lts1,nvarts1,nelets1
         WRITE(io_timeser_bgc,*) (rlonts1(l),l=1,nelets1)
         WRITE(io_timeser_bgc,*) (rlatts1(l),l=1,nelets1)

         DO l2=1,nelets1 !outer vertical loop is station
         DO l3=1,lts1 !inner vertical loop is time
           WRITE(io_timeser_bgc,433)(ts1(l1,l2,l3),l1=1,nvarts1) !horizontal loop is variable type
         ENDDO
         ENDDO

         CLOSE(io_timeser_bgc)
433      FORMAT(19e15.6)

      ENDIF
#endif
      ENDIF ! p_pe==p_io

         WRITE(io_stdo_bgc,*)'Memory deallocation ts1 ...'
         WRITE(io_stdo_bgc,*)'No. of variables: ',nvarts1
         WRITE(io_stdo_bgc,*)'No. of stations : ',nelets1
         WRITE(io_stdo_bgc,*)'No. of time steps: ',lts1
         DEALLOCATE (ts1)

      !CALL p_barrier ! Not really necessary, just to wait for the I/O-PE

      RETURN
      END
