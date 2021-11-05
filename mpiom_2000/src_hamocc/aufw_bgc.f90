      SUBROUTINE AUFW_BGC(kpie,kpje,kpke,pddpo,pgila,pgiph,ptiestu   &
     &                    ,kplyear,kplmon,kplday,kpldtoce)
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/aufw_bgc.f90,v $\\
!$Revision: 1.3.2.1.16.1.4.2 $\\
!$Date: 2006/04/06 10:01:07 $\\

!****************************************************************
!
!**** *AUFW_BGC* - write marine bgc restart data.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - extra SBR for writing bgc data to the restart file.
!     S.Legutke,        *MPI-MaD, HH*    15.08.01
!     - netCDF version (cond.comp. PNETCDF)
!     - chemcm is multiplied with layer-dependent constant in order
!       to be displayable by ncview. It is not read in AUFR_BGC!
!
!     Purpose
!     -------
!     Write restart data for continuation of interrupted integration.
!
!     Method
!     -------
!     The bgc data are written to an extra file, other than the ocean data.
!     The netCDF version also writes grid description variables.
!     The netCDF file is selfdescribing and can be used for
!     visualization (e.g. STOMPP, grads). Before the data are saved, all
!     values at dry cells are set to rmissing (e.g.rmasko).
!     The time stamp of the bgc restart file (idate) is taken from the
!     ocean time stamp through the SBR parameter list. The only time
!     control variable proper to the bgc is the time step number (idate(5)).
!     It can differ from that of the ocean (idate(4)) by the difference
!     of the offsets of restart files.

!       changed : bgc uses ocean time
!
!**   Interface.
!     ----------
!
!     *CALL*       *AUFW_BGC(kpie,kpje,kpke,pddpo,pgila,pgiph,ptiestu)
!
!     *COMMON*     *MO_PARAM1_BGC* - declaration of ocean/sediment tracer.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu* - depth of layers [m].
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
!ik IK introduced ocetra array elements for phyto, grazer, poc (=det), calciu
!ik array indices are: iphy, izoo, idet, icalc
!iktodo IK introduced new variable opal (index iopal)
!ik nocetra is the number of all BGC element (size of ocetra(,,,l))
!ik nosedi is the number of all elements interacting with the sediment


      USE mo_kind, ONLY: wp
      USE mo_biomod, ONLY: rnit, kbo
      USE mo_carbch, ONLY: ocetra, hi, co3, aksp, atm, suppco2
      USE mo_sedmnt, ONLY: powtra, sedlay, sedhpl, burial
      USE mo_control_bgc, ONLY: io_stdo_bgc, ldtbgc, rmasko, rmasks
      USE mo_param1_bgc, ONLY: ipown2, ipowaal, issso12, ipowasi, isssc12, ks, &
           issster, ipowaic, nocetra, ipowaph, issssil, ipowaox, ipowno3

#ifdef OLD_IO
      USE mo_commo1, ONLY: ie_g,je_g,ldays,lmonts,lyears
      USE mo_constants, ONLY: aradtogra
      USE mo_parallel, ONLY : p_io, p_pe, stop_all
#endif
      USE mo_boundsexch, ONLY : bounds_exch

      implicit none
      INTEGER, INTENT(in) :: kpie,kpje,kpke
      REAL(wp), INTENT(in) :: pddpo(kpie,kpje,kpke)
      REAL(wp), INTENT(in) :: pgila(kpie*2,kpje*2)
      REAL(wp), INTENT(in) :: pgiph(kpie*2,kpje*2)
      REAL(wp), INTENT(in) :: ptiestu(kpke)
      INTEGER, INTENT(in) :: kplyear,kplmon,kplday,kpldtoce
#ifdef OLD_IO
!      REAL(wp) chemcm_g(ie_g,je_g,8,12)
      INTEGER jj,ii
      integer icode(1)
#endif
      INTEGER :: i,j,k,l,ll
      INTEGER :: kpicycli
      INTEGER :: idate(5)
      REAL(wp) :: rmissing
#ifdef OLD_IO
#ifdef PNETCDF
      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncvarid,ncstat,ncoldmod,ncdims(4)                  &
     &       ,nclatid,nclonid,nclevid,nclev1id                        &
     &       ,nctraid,nccheid,ncksid,ncsedid                          &
     &       ,ncmonid,nstart2(2),ncount2(2),nstride2(2)

      INTEGER timeid, START_1D,COUNT_1D
      REAL(wp) timeaufw
      REAL(wp) zfield(kpie,kpje)


#endif
#endif
!
! Check of fields written to restart file.
!
      CALL CHCK_BGC(io_stdo_bgc,kpicycli,                                 &
     &'Check values of ocean/sediment tracer written to restart file :',  &
     & kpie,kpje,kpke,pddpo)

      idate(1) = kplyear
      idate(2) = kplmon
      idate(3) = kplday
      idate(4) = kpldtoce
      idate(5) = ldtbgc
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Writing restart file at date :'              &
     &,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_bgc,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_bgc,*) 'Bgc   model step number is ',idate(5)
!
!  Masking aqueous sea water tracer.
!

      rmissing = rmasko
      rmasks   = rmasko


      DO l=1,nocetra
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF (pddpo(i, j, k) .LT. 0.5_wp) THEN
         ocetra(i,j,k,l)=rmasko
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

!
!  Masking other sea water tracer.
!
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF (pddpo(i, j, k) .LT. 0.5_wp) THEN
         hi    (i,j,k)=rmissing
         co3   (i,j,k)=rmissing
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!
!  Masking sediment pore water tracer.
!
      DO  k=1,ks
      DO  j=1,kpje
      DO  i=1,kpie
      IF(kbo(i,j) .LE. 1) THEN
         powtra(i, j, k, ipowno3) = 2.17e-6_wp * rnit
         powtra(i, j, k, ipown2) = 0._wp
         powtra(i, j, k, ipowaic) = 2.27e-3_wp
#ifdef __c_isotopes
         powtra(i, j, k, ipowc13) = 2.27e-3_wp
         powtra(i, j, k, ipowc14) = 2.27e-3_wp
#endif
         powtra(i, j, k, ipowaal) = 2.37e-3_wp
         powtra(i, j, k, ipowaph) = 2.17e-6_wp
         powtra(i, j, k, ipowaox) = 3.e-4_wp
         powtra(i, j, k, ipowasi) = 1.2e-4_wp

         sedlay(i, j, k, issso12) = 1.e-8_wp
#ifdef __c_isotopes
         sedlay(i, j, k, issso13) = 1.e-8_wp
         sedlay(i, j, k, issso14) = 1.e-8_wp
#endif

         sedlay(i, j, k, isssc12) = 1.e-8_wp
#ifdef __c_isotopes
         sedlay(i, j, k, isssc13) = 1.e-8_wp
         sedlay(i, j, k, isssc14) = 1.e-8_wp
#endif

         sedlay(i, j, k, issssil) = 3.e-3_wp
         sedlay(i, j, k, issster) = 30._wp
         burial(i, j, issso12) = 0._wp
         burial(i, j, isssc12) = 0._wp
         burial(i, j, issssil) = 0._wp
         burial(i, j, issster) = 0._wp
         sedhpl(i, j, k) = 0._wp
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!
!  Masking co2.
!
      DO  j=1,kpje
      DO  i=1,kpie
      IF (pddpo(i, j, 1) .LT. 0.5_wp) THEN
      suppco2(i,j)=rmissing
      ENDIF
      ENDDO
      ENDDO

#ifdef PNETCDF

#ifdef OLD_IO /* part 1 */
!
! Open netCDF data file
!
      IF(p_pe==p_io) THEN
!use large files
       ncstat = nf_create('restartw_bgc.nc',IOR(NF_CLOBBER,NF_64BIT_OFFSET),ncid)
!      ncstat = NF_CREATE('restartw_bgc.nc',NF_CLOBBER, ncid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF1')

!
! Define dimension
!
      ncstat = NF_DEF_DIM(ncid, 'lon', ie_g, nclonid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF2')

      ncstat = NF_DEF_DIM(ncid, 'lat', je_g, nclatid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF3')

      ncstat = NF_DEF_DIM(ncid, 'depth', kpke, nclevid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF4')

      ncstat = NF_DEF_DIM(ncid, 'ntra', nocetra, nctraid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF5')

      ncstat = NF_DEF_DIM(ncid, 'nche', 8, nccheid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF6')

      ncstat = NF_DEF_DIM(ncid, 'nks', ks, ncksid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF7')

      ncstat = NF_DEF_DIM(ncid, 'nsed',nsedtra , ncsedid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF8')
! Paddy:
      ncstat = NF_DEF_DIM(ncid, 'nmon',12 , ncmonid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF8a')

      ncstat = NF_DEF_DIM(ncid, 'lev1', 1, nclev1id)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF8b')


      ncstat = NF_DEF_DIM(ncid, 'time',NF_UNLIMITED, timeid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('AUFW: Problem with netCDF7'      )
!

!
! Define global attributes
!
      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'title'               &
     &,35, 'restart data for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'history'             &
     &,35, 'restart data for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'conventions'         &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'source'              &
     &,24, 'Marine bgc model output MPI-OM/grob')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF10')

      ncstat = NF_PUT_ATT_INT(ncid, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF11')

!
! Define variables : grid
!
      ncdims(1) = nclonid
      ncdims(2) = nclatid

      ncstat = NF_DEF_VAR(ncid,'scal_lon',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF13')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',8, 'degree E')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF14')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,34, '2-d longitude of scalar grid cells')
      icode=1
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF15')


      ncstat = NF_DEF_VAR(ncid,'scal_lat',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF16')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',8, 'degree N')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF17')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,33, '2-d latitude of scalar grid cells')
      icode=2
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF18')

      ncstat = NF_DEF_VAR(ncid,'scal_wdep',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF16a')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF17a')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,37, '2-d water depth at scalar grid points')
      icode=3
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF18a')




      ncdims(1) = nclevid

      ncstat = NF_DEF_VAR(ncid,'depth',NF_DOUBLE,1,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF19')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF20')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,32, '1-d layer depths of ocean tracer')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF21')

      ncdims(1) = timeid

      CALL NETCDF_DEF_VARDB(ncid,4,'time',1,ncdims,ncvarid,          &
     &   16,'day as %Y%m%d.%f',4,'time',rmasko,31,io_stdo_bgc)


      ncdims(1) = 1

      ncstat = NF_DEF_VAR(ncid,'calcinpglint',NF_DOUBLE,0,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF22')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',4,'kmol')
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF23')
      ncstat = NF_DEF_VAR(ncid,'orginpglint',NF_DOUBLE,0,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF24')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',4,'kmol')
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF25')
      ncstat = NF_DEF_VAR(ncid,'silinpglint',NF_DOUBLE,0,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF26')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',4,'kmol')
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF27')

!
! Define variables : advected ocean tracer
!
      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid


      CALL NETCDF_DEF_VARDB(ncid,6,'sco212',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',13, 'Dissolved CO2',rmissing,22,io_stdo_bgc)
      icode=7
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF28')

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'sco213',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',15, 'Dissolved CO213',rmissing,22,io_stdo_bgc)
      icode=8
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF28')

      CALL NETCDF_DEF_VARDB(ncid,6,'sco214',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',15, 'Dissolved CO214',rmissing,22,io_stdo_bgc)
      icode=9
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF29')
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'alkali',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',10,'Alkalinity',rmissing,25,io_stdo_bgc)
      icode=10
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF30')

      CALL NETCDF_DEF_VARDB(ncid,6,'phosph',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',19,'Dissolved phosphate',rmissing,28,io_stdo_bgc)
      icode=11
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF31')

      CALL NETCDF_DEF_VARDB(ncid,6,'oxygen',3,ncdims,ncvarid,           &
     &9,'kmol/m**3',16,'Dissolved oxygen',rmissing,31,io_stdo_bgc)
      icode=12
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF32')

      CALL NETCDF_DEF_VARDB(ncid,6,'gasnit',3,ncdims,ncvarid,           &
      9,'kmol/m**3',21,'Gaseous nitrogen (N2)',rmissing,34,io_stdo_bgc)
      icode=13
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF33')

      CALL NETCDF_DEF_VARDB(ncid,4,'ano3',3,ncdims,ncvarid,             &
      9,'kmol/m**3',17,'Dissolved nitrate',rmissing,34,io_stdo_bgc)
      icode=14
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF34')

      CALL NETCDF_DEF_VARDB(ncid,6,'silica',3,ncdims,ncvarid,           &
      9,'kmol/m**3',22,'Silicid acid (Si(OH)4)'                         &
      ,rmissing,40,io_stdo_bgc)
      icode=15
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF35')


      CALL NETCDF_DEF_VARDB(ncid,3,'doc',3,ncdims,ncvarid,              &
      9,'kmol/m**3',24,'Dissolved organic carbon',                      &
      rmissing,40,io_stdo_bgc)
      icode=16
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF36')

      CALL NETCDF_DEF_VARDB(ncid,3,'poc',3,ncdims,ncvarid,              &
      10,'kmolP/m**3',26,'Particulate organic carbon',                    &
      rmissing,46,io_stdo_bgc)
      icode=17
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF37')

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,5,'poc13',3,ncdims,ncvarid,              &
      10,'kmolC/m**3',28,'Particulate organic carbon13',                    &
      rmissing,46,io_stdo_bgc)
      icode=18
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF38')

      CALL NETCDF_DEF_VARDB(ncid,5,'poc14',3,ncdims,ncvarid,              &
      10,'kmolC/m**3',28,'Particulate organic carbon14',                    &
      rmissing,46,io_stdo_bgc)
      icode=19
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode )
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF39')
#endif

      CALL NETCDF_DEF_VARDB(ncid,2,'hi',3,ncdims,ncvarid,               &
      9,'kmol/m**3',26,'Hydrogen ion concentration',                    &
      rmissing,46,io_stdo_bgc)
      icode=20
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF40')

      CALL NETCDF_DEF_VARDB(ncid,3,'co3',3,ncdims,ncvarid,              &
      9,'kmol/m**3',25,'Dissolved carbonate (CO3)',                     &
      rmissing,52,io_stdo_bgc)
      icode=21
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF41')

      CALL NETCDF_DEF_VARDB(ncid,5,'phyto',3,ncdims,ncvarid,            &
     &    10,'kmolP/m**3',27,'Phytoplankton concentration',             &
     &    rmissing,28,io_stdo_bgc)
      icode=22
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF42')

      CALL NETCDF_DEF_VARDB(ncid,6,'grazer',3,ncdims,ncvarid,           &
     &    10,'kmolP/m**3',25,'Zooplankton concentration',               &
     &    rmissing,29,io_stdo_bgc)
      icode=23
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF43')

      CALL NETCDF_DEF_VARDB(ncid,6,'calciu',3,ncdims,ncvarid,           &
     &    10,'kmolP/m**3',17,'Calcium carbonate',                         &
     &    rmissing,30,io_stdo_bgc)
      icode=24
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF44')

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,8,'calciu13',3,ncdims,ncvarid,         &
     &    10,'kmolC/m**3',19,'Calcium carbonate13',                       &
     &    rmissing,30,io_stdo_bgc
      icode=25
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF45')

      CALL NETCDF_DEF_VARDB(ncid,8,'calciu14',3,ncdims,ncvarid,         &
     &    10,'kmolC/m**3',19,'Calcium carbonate14',                       &
     &    rmissing,30,io_stdo_bgc)
      icode=26
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF46')
#endif

      CALL NETCDF_DEF_VARDB(ncid,4,'opal',3,ncdims,ncvarid,             &
     &    9,'kmol/m**3',15,'Biogenic silica',                           &
     &    rmissing,31,io_stdo_bgc)
      icode=27
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF47')

      CALL NETCDF_DEF_VARDB(ncid,3,'n2o',3,ncdims,ncvarid,              &
     &    9,'kmol/m**3',12,'laughing gas',                              &
     &    rmissing,32,io_stdo_bgc)
      icode=28
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF48')

      CALL NETCDF_DEF_VARDB(ncid,3,'dms',3,ncdims,ncvarid,              &
     &    9,'kmol/m**3',15 ,'DiMethylSulfide',                          &
     &    rmissing,33,io_stdo_bgc)
      icode=29
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF49')

      CALL NETCDF_DEF_VARDB(ncid,5,'fdust',3,ncdims,ncvarid,            &
     &    7,'kg/m**3',19,'Non-aggregated dust',                         &
     &    rmissing,34,io_stdo_bgc)
      icode=30
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF50')

      CALL NETCDF_DEF_VARDB(ncid,4,'iron',3,ncdims,ncvarid,             &
     &    9,'kmol/m**3',14,'Dissolved iron',                            &
     &    rmissing,35,io_stdo_bgc)
      icode=31
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF50')

#ifdef AGG
      CALL NETCDF_DEF_VARDB(ncid,4,'snos',3,ncdims,ncvarid,             &
     &    8,'1/cm**3',30,'marine snow aggregates per cm3',              &
     &    rmissing,41,io_stdo_bgc)
      icode=32
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF51')

      CALL NETCDF_DEF_VARDB(ncid,5,'adust',3,ncdims,ncvarid,            &
     &    7,'kg/m**3',15,'Aggregated dust',                             &
     &    rmissing,42,io_stdo_bgc)
      icode=33
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF52')
#endif /*AGG*/

#ifdef ANTC14
      CALL NETCDF_DEF_VARDB(ncid,6,'antc14',3,ncdims,ncvarid,           &
     &    9,'kmol/m**3',17,'anthropogenic C14',                         &
     &    rmissing,41,io_stdo_bgc)
      icode=34
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF53')
#endif
#ifdef PCFC
      CALL NETCDF_DEF_VARDB(ncid,5,'cfc11',3,ncdims,ncvarid,             &
     &    9,'kmol/m**3',5,'CFC11',                                       &
     &    rmissing,41,io_stdo_bgc)
      icode=35
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF54')

      CALL NETCDF_DEF_VARDB(ncid,5,'cfc12',3,ncdims,ncvarid,            &
         9,'kmol/m**3',5,'CFC12',                                       &
         rmissing,41,io_stdo_bgc)
      icode=36
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF55')
#endif

!
! Define variables : aksp
! ----------------------------------------------------------------------
      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid

      CALL NETCDF_DEF_VARDB(ncid,4,'aksp',3,ncdims,ncvarid,             &
         9,'XXXXXXXXX',39 ,'apparent solubility product for calcite',   &
         rmissing,64,io_stdo_bgc)
      icode=37
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF56')

! Define variables : sediment
! ----------------------------------------------------------------------
      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = ncksid
      ncdims(4) = 0

      CALL NETCDF_DEF_VARDB(ncid,6,'ssso12',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',35,'Sediment accumulated organic carbon',      &
     &    rmasks,69,io_stdo_bgc)
      icode=38
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF57')

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'ssso13',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',37,'Sediment accumulated organic carbon13',    &
     &    rmasks,69,io_stdo_bgc)
      icode=39
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF58')

      CALL NETCDF_DEF_VARDB(ncid,6,'ssso14',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',37,'Sediment accumulated organic carbon14',    &
     &    rmasks,69,io_stdo_bgc)
      icode=40
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF59')
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'sssc12',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',38,'Sediment accumulated calcium carbonate',   &
     &    rmasks,69,io_stdo_bgc)
      icode=41
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF60')

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'sssc13',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',40,'Sediment accumulated calcium carbonate13', &
     &    rmasks,69,io_stdo_bgc)
      icode=42
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF61')

      CALL NETCDF_DEF_VARDB(ncid,6,'sssc14',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',40,'Sediment accumulated calcium carbonate14', &
     &    rmasks,69,io_stdo_bgc)
      icode=43
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF62')
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'ssssil',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',25,'Sediment accumulated opal',                &
     &    rmasks,69,io_stdo_bgc)
      icode=44
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF63')

      CALL NETCDF_DEF_VARDB(ncid,6,'ssster',3,ncdims,ncvarid,          &
     &    7,'kg/m**3',25,'Sediment accumulated clay',                &
     &    rmasks,69,io_stdo_bgc)
      icode=45
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF64')

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = 0
      ncdims(4) = 0

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_o12',2,ncdims,ncvarid,         &
     &    9,'kmol/m**2',30,'Burial layer of organic carbon',           &
     &    rmasks,70,io_stdo_bgc)
      icode=46
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF65')

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_c12',2,ncdims,ncvarid,         &
     &    9,'kmol/m**2',33,'Burial layer of calcium carbonate',        &
     &    rmasks,71,io_stdo_bgc)
      icode=47
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF66')

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_sil',2,ncdims,ncvarid,         &
     &    9,'kmol/m**2',20,'Burial layer of opal',                     &
     &    rmasks,72,io_stdo_bgc)
      icode=48
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF67')

      CALL NETCDF_DEF_VARDB(ncid,8,'bur_clay',2,ncdims,ncvarid,        &
     &    7,'kg/m**3',20,'Burial layer of clay',                     &
     &    rmasks,73,io_stdo_bgc)
      icode=49
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF68')

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = ncksid

      CALL NETCDF_DEF_VARDB(ncid,6,'sedhpl',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',34,'Sediment accumulated hydrogen ions',       &
     &    rmasks,74,io_stdo_bgc)
      icode=50
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF69')

      CALL NETCDF_DEF_VARDB(ncid,6,'powaic',3,ncdims,ncvarid,          &
     &    9,'kmol/m**3',23,'Sediment pore water DIC',                  &
     &    rmasks,75,io_stdo_bgc)
      icode=51
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF70')
#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'powc13',3,ncdims,ncvarid,          &
     &    9,'kmol/m**3',25,'Sediment pore water DIC13',                &
     &    rmasks,75,io_stdo_bgc)
      icode=52
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF71')

      CALL NETCDF_DEF_VARDB(ncid,6,'powc14',3,ncdims,ncvarid,          &
     &    9,'kmol/m**3',25,'Sediment pore water DIC13',                &
     &    rmasks,75,io_stdo_bgc)
      icode=53
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF72')
#endif
      CALL NETCDF_DEF_VARDB(ncid,6,'powaal',3,ncdims,ncvarid,       &
     &    7,'eq/m**3',30,'Sediment pore water alkalinity',        &
     &    rmasks,78,io_stdo_bgc)
      icode=54
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF73')

      CALL NETCDF_DEF_VARDB(ncid,6,'powaph',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',29,'Sediment pore water phosphate',         &
     &    rmasks,79,io_stdo_bgc)
      icode=55
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF74')

      CALL NETCDF_DEF_VARDB(ncid,6,'powaox',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',26,'Sediment pore water oxygen',            &
     &    rmasks,84,io_stdo_bgc)
      icode=56
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF75')

      CALL NETCDF_DEF_VARDB(ncid,5,'pown2',3,ncdims,ncvarid,        &
     &    9,'kmol/m**3',36,'Sediment pore water gaseous nitrogen',  &
     &    rmasks,87,io_stdo_bgc)
      icode=57
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF76')

      CALL NETCDF_DEF_VARDB(ncid,6,'powno3',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',33,'Sediment pore water nitrate (NO3)',     &
     &    rmasks,90,io_stdo_bgc)
      icode=58
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF77')

      CALL NETCDF_DEF_VARDB(ncid,6,'powasi',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',42,'Sediment pore water silicid acid (Si(OH)4)',         &
     &    rmasks,91,io_stdo_bgc)
      icode=59
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF78')

!
! Define variables : co2 diffusion
!
      CALL NETCDF_DEF_VARDB(ncid,7,'suppco2',2,ncdims,ncvarid,      &
     &    4,'ppmv',42,'pCO2 from total dissolved inorganic carbon', &
     &    rmissing,92,io_stdo_bgc)
      icode=60
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF79')

      CALL NETCDF_DEF_VARDB(ncid,6,'atmco2',2,ncdims,ncvarid,       &
     &    3,'ppm',15,'atmospheric CO2',                             &
     &    rmissing,93,io_stdo_bgc)
      icode=61
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF80')

      CALL NETCDF_DEF_VARDB(ncid,5,'atmo2',2,ncdims,ncvarid,        &
     &    3,'ppm',14,'atmospheric O2',                              &
     &    rmissing,94,io_stdo_bgc)
      icode=62
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF81')

      CALL NETCDF_DEF_VARDB(ncid,5,'atmn2',2,ncdims,ncvarid,        &
     &    3,'ppm',14,'atmospheric N2',                              &
     &    rmissing,95,io_stdo_bgc)
      icode=64
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF82')

      ncstat = NF_ENDDEF(ncid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF00')


!
! Set fill mode
!
      ncstat = NF_SET_FILL(ncid,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF97')

      ENDIF ! p_pe == p_io

!
! Write grid describing data
!

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

      CALL write_netcdf_var(ncid,'scal_lon',zfield(1,1),1,0)


      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgiph(ii,jj)*aradtogra
         ENDDO
      ENDDO

     CALL write_netcdf_var(ncid,'scal_lat', zfield(1,1),1,0)

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

      CALL write_netcdf_var(ncid,'scal_wdep',zfield(1,1),1,0)



      IF(p_pe==p_io) THEN
        ncstat = NF_INQ_VARID(ncid,'depth',ncvarid )
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all('AUFW: netCDF102')
        ncstat = NF_PUT_VAR_DOUBLE (ncid,ncvarid,ptiestu(1) )
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all('AUFW: netCDF103')




       timeaufw=REAL(LDAYS + LMONTS * 100 + LYEARS * 10000)
       START_1D =1
       COUNT_1D =1


       ncstat = NF_INQ_VARID(ncid,'time',ncvarid )
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF103d')
       ncstat = NF_PUT_VARA_DOUBLE (ncid,ncvarid,start_1d,count_1d,timeaufw)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF103d')


       ncstat = NF_INQ_VARID(ncid,'calcinpglint',ncvarid)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF104a1')
       ncstat = NF_PUT_VARA_DOUBLE (ncid,ncvarid,start_1d,count_1d,calcinpglint)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF104a2')
       ncstat = NF_INQ_VARID(ncid,'orginpglint',ncvarid)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF104b1')
       ncstat = NF_PUT_VARA_DOUBLE (ncid,ncvarid,start_1d,count_1d,orginpglint)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF104b2')
       ncstat = NF_INQ_VARID(ncid,'silinpglint',ncvarid)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF104c1')
       ncstat = NF_PUT_VARA_DOUBLE (ncid,ncvarid,start_1d,count_1d,silinpglint)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF104c2')

      ENDIF







!
! Write restart data : ocean aquateous tracer
!

#ifdef uwe_restart
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie

! replace rmasko by start value so restart file can be used for newly introduced wet points (e.g., lakes)
      IF(pddpo(i,j,k) .LT. 0.5) THEN       ! dry points
         ocetra(i,j,k,isco212)=2.27e-3     ! [kmol/m3]
#ifdef __c_isotopes
         ocetra(i,j,k,isco213)=2.27e-3     ! adjusted to reference ratio 13C/12C=1 (*100)!
         ocetra(i,j,k,isco214)=2.27e-3
#endif /*__c_isotopes*/
         ocetra(i,j,k,ialkali)=2.37e-3
         ocetra(i,j,k,iphosph)=2.17e-6
         ocetra(i,j,k,ioxygen)=3.e-4
         ocetra(i,j,k,igasnit)=0.0
         ocetra(i,j,k,iano3)  =2.17e-6*rnit ! old 32.e-6
         ocetra(i,j,k,isilica)=1.2e-4
         ocetra(i,j,k,idoc)   =1.e-10
         ocetra(i,j,k,iphy)   =1.e-8
         ocetra(i,j,k,izoo)   =1.e-8
         ocetra(i,j,k,idet)   =1.e-8
         ocetra(i,j,k,icalc)  =0.
#ifdef __c_isotopes
         ocetra(i,j,k,idet13) =1.e-8
         ocetra(i,j,k,icalc13)=0.
         ocetra(i,j,k,idet14) =1.e-8
         ocetra(i,j,k,icalc14)=0.
         ocetra(i,j,k,isco214)=0.75*2.27e-3 !Paddy: oldest water: 1600y --> X0.83
#endif   /*__c_isotopes*/
         ocetra(i,j,k,iopal)  =1.e-8
         ocetra(i,j,k,ian2o)  =0.
         ocetra(i,j,k,idms)   =0.
         ocetra(i,j,k,ifdust) =0.
         ocetra(i,j,k,iiron)  =0.6*1.e-9
!        ocetra(i,j,k,ibeten) =0.6*1.e-12
         hi(i,j,k)            =3.e-9
         co3(i,j,k)           =1.e-4
      ENDIF

      ENDDO
      ENDDO
      ENDDO
#endif/*uwe_restart*/
#endif/*def OLD_IO part 1*/

      ll=size(ocetra(1,1,1,:))
      do l = 1, ll
       call bounds_exch(1,'p',ocetra(:,:,:,l),'in aufw 1')
      enddo


#ifdef OLD_IO /* part 2 */
      CALL write_netcdf_var(ncid,'sco212',ocetra(1,1,1,isco212),kpke,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'sco213',ocetra(1,1,1,isco213),kpke,0)
      CALL write_netcdf_var(ncid,'sco214',ocetra(1,1,1,isco214),kpke,0)
#endif
      CALL write_netcdf_var(ncid,'alkali',ocetra(1,1,1,ialkali),kpke,0)
      CALL write_netcdf_var(ncid,'phosph',ocetra(1,1,1,iphosph),kpke,0)
      CALL write_netcdf_var(ncid,'oxygen',ocetra(1,1,1,ioxygen),kpke,0)
      CALL write_netcdf_var(ncid,'gasnit',ocetra(1,1,1,igasnit),kpke,0)
      CALL write_netcdf_var(ncid,'ano3',ocetra(1,1,1,iano3),kpke,0)
      CALL write_netcdf_var(ncid,'silica',ocetra(1,1,1,isilica),kpke,0)
      CALL write_netcdf_var(ncid,'doc',ocetra(1,1,1,idoc),kpke,0)
      CALL write_netcdf_var(ncid,'poc',ocetra(1,1,1,idet),kpke,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'poc13',ocetra(1,1,1,idet13),kpke,0)
      CALL write_netcdf_var(ncid,'poc14',ocetra(1,1,1,idet14),kpke,0)
#endif
      CALL write_netcdf_var(ncid,'phyto',ocetra(1,1,1,iphy),kpke,0)
      CALL write_netcdf_var(ncid,'grazer',ocetra(1,1,1,izoo),kpke,0)
      CALL write_netcdf_var(ncid,'calciu',ocetra(1,1,1,icalc),kpke,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'calciu13',ocetra(1,1,1,icalc13),kpke,0)
      CALL write_netcdf_var(ncid,'calciu14',ocetra(1,1,1,icalc14),kpke,0)
#endif
      CALL write_netcdf_var(ncid,'opal',ocetra(1,1,1,iopal),kpke,0)
      CALL write_netcdf_var(ncid,'n2o',ocetra(1,1,1,ian2o),kpke,0)
      CALL write_netcdf_var(ncid,'dms',ocetra(1,1,1,idms),kpke,0)
      CALL write_netcdf_var(ncid,'fdust',ocetra(1,1,1,ifdust),kpke,0)
      CALL write_netcdf_var(ncid,'iron',ocetra(1,1,1,iiron),kpke,0)
#ifdef AGG
      CALL write_netcdf_var(ncid,'snos',ocetra(1,1,1,inos),kpke,0)
      CALL write_netcdf_var(ncid,'adust',ocetra(1,1,1,iadust),kpke,0)
#endif /*AGG*/
#ifdef ANTC14
      CALL write_netcdf_var(ncid,'antc14',ocetra(1,1,1,iantc14),kpke,0)
#endif
#ifdef PCFC
      CALL write_netcdf_var(ncid,'cfc11',ocetra(1,1,1,icfc11),kpke,0)
      CALL write_netcdf_var(ncid,'cfc12',ocetra(1,1,1,icfc12),kpke,0)
#endif
#endif/*def OLD_IO part 2 */

!
! Write restart data : diagnostic ocean tracer
!
      call bounds_exch(1,'p',hi(:,:,:),'in aufw 3')
      call bounds_exch(1,'p',co3(:,:,:),'in aufw 3')


#ifdef OLD_IO /* part 3 */
      CALL write_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0)
      CALL write_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0)
#endif/*def OLD_IO part 3 */
!
! Write restart data : other fields
!
      call bounds_exch(1,'p',aksp(:,:,:),'in aufw 3')

#ifdef OLD_IO /* part 4 */
      CALL write_netcdf_var(ncid,'aksp',aksp(1,1,1),kpke,0)
#endif/*def OLD_IO part 4 */

      do k=1,ks
         ll=size(sedlay(1,1,1,:))
         do l = 1, ll
            call bounds_exch(1,'p',sedlay(:,:,k,l),'in aufw 1')
         enddo

         ll=size(powtra(1,1,1,:))
         do l = 1, ll
            call bounds_exch(1,'p',powtra(:,:,k,l),'in aufw 2')
         enddo
         call bounds_exch(1,'p',sedhpl(:,:,k),'in aufw 4')
      enddo

      ll=size(burial(1,1,:))
      do l = 1, ll
         call bounds_exch(1,'p',burial(:,:,l),'in aufw 3')
      enddo

      call bounds_exch(1,'p',suppco2(:,:),'in aufw 5')

      ll=size(atm(1,1,:))
      do l = 1, ll
         call bounds_exch(1,'p',atm(:,:,l),'in aufw 6')
      enddo


#ifdef OLD_IO /* part 5 */
! Write restart data : sediment variables.
!
      CALL write_netcdf_var(ncid,'ssso12',sedlay(1,1,1,issso12),ks,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'ssso13',sedlay(1,1,1,issso13),ks,0)
      CALL write_netcdf_var(ncid,'ssso14',sedlay(1,1,1,issso14),ks,0)
#endif
      CALL write_netcdf_var(ncid,'sssc12',sedlay(1,1,1,isssc12),ks,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'sssc13',sedlay(1,1,1,isssc13),ks,0)
      CALL write_netcdf_var(ncid,'sssc14',sedlay(1,1,1,isssc14),ks,0)
#endif
      CALL write_netcdf_var(ncid,'ssssil',sedlay(1,1,1,issssil),ks,0)
      CALL write_netcdf_var(ncid,'ssster',sedlay(1,1,1,issster),ks,0)
      CALL write_netcdf_var(ncid,'bur_o12',burial(1,1,issso12),1,0)
      CALL write_netcdf_var(ncid,'bur_c12',burial(1,1,isssc12),1,0)
      CALL write_netcdf_var(ncid,'bur_sil',burial(1,1,issssil),1,0)
      CALL write_netcdf_var(ncid,'bur_clay',burial(1,1,issster),1,0)
      CALL write_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0)
      CALL write_netcdf_var(ncid,'powaic',powtra(1,1,1,ipowaic),ks,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'powc13',powtra(1,1,1,ipowc13),ks,0)
      CALL write_netcdf_var(ncid,'powc14',powtra(1,1,1,ipowc14),ks,0)
#endif
      CALL write_netcdf_var(ncid,'powaal',powtra(1,1,1,ipowaal),ks,0)
      CALL write_netcdf_var(ncid,'powaph',powtra(1,1,1,ipowaph),ks,0)
      CALL write_netcdf_var(ncid,'powaox',powtra(1,1,1,ipowaox),ks,0)
      CALL write_netcdf_var(ncid,'pown2',powtra(1,1,1,ipown2),ks,0)
      CALL write_netcdf_var(ncid,'powno3',powtra(1,1,1,ipowno3),ks,0)
      CALL write_netcdf_var(ncid,'powasi',powtra(1,1,1,ipowasi),ks,0)
!      CALL write_netcdf_var(ncid,'suppco2',suppco2(1,1),1,0)
      CALL write_netcdf_var(ncid,'atmco2',atm(1,1,iatmco2),1,0)
      CALL write_netcdf_var(ncid,'atmo2',atm(1,1,iatmo2),1,0)
      CALL write_netcdf_var(ncid,'atmn2',atm(1,1,iatmn2),1,0)

      IF(p_pe==p_io) THEN
        ncstat = NF_CLOSE(ncid)
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all('AUFW: netCDF200')
      ENDIF
#endif/*def OLD_IO part 5 */

#else/* def PNETCDF */

#ifdef OLD_IO /* part 6 */
      ! no netcdf
      ! does this work with the MPI anymore - c.f. read_netcdf_var ??
      OPEN(io_rsto_bgc,FILE='restart_bgc',STATUS='UNKNOWN'            &
     &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      WRITE(io_rsto_bgc)                                              &
     &             (((ocetra(i,j,k,iphosph),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,isilica),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,ioxygen),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,iphy   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,izoo   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,idet   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,idoc   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,icalc  ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,iano3  ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,igasnit),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,iopal  ),i=1,kpie),j=1,kpje),k=1,kpke)

      WRITE(io_rsto_bgc) chemcm,hi,co3,aksp                              &
     &            ,(((ocetra(i,j,k,isco212),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocetra(i,j,k,ialkali),i=1,kpie),j=1,kpje),k=1,kpke)

      WRITE(io_rsto_bgc) sedlay,sedhpl

      WRITE(io_rsto_bgc)                                               &
     &             (((powtra(i,j,k,ipowaic),i=1,kpie),j=1,kpje),k=1,ks)&
     &            ,(((powtra(i,j,k,ipowaal),i=1,kpie),j=1,kpje),k=1,ks)&
     &            ,(((powtra(i,j,k,ipowaph),i=1,kpie),j=1,kpje),k=1,ks)&
     &            ,(((powtra(i,j,k,ipowaox),i=1,kpie),j=1,kpje),k=1,ks)&
     &            ,(((powtra(i,j,k,ipowasi),i=1,kpie),j=1,kpje),k=1,ks)&
     &            ,(((powtra(i,j,k,ipowno3),i=1,kpie),j=1,kpje),k=1,ks)&
     &            ,(((powtra(i,j,k,ipown2) ,i=1,kpie),j=1,kpje),k=1,ks)

      REWIND io_rsto_bgc

      CLOSE(io_rsto_bgc)
#endif/*def OLD_IO part 6*/

#endif/* else def PNETCDF */

      WRITE(io_stdo_bgc,*) 'End of AUFW_BGC'
      WRITE(io_stdo_bgc,*) '***************'

      RETURN
      END
