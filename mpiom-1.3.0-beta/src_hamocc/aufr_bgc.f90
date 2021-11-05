      SUBROUTINE AUFR_BGC(kpie,kpje,kpke,pddpo,kplyear,kplmon     &
     &                    ,kplday,kpldtoce)

!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/aufr_bgc.f90,v $\\
!$Revision: 1.3.2.1.16.1.4.2 $\\
!$Date: 2006/04/06 10:01:07 $\\

!****************************************************************
!
!**** *AUFR_BGC* - reads marine bgc restart data.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - SBR for reading bgc data from the restart file.
!     S.Legutke,        *MPI-MaD, HH*    15.08.01
!     - netCDF version (with cond.comp. PNETCDF)
!     - no use of chemc values from netCDF restart
!
!     Patrick Wetzel,    *MPI-Met, HH*    16.04.02
!     - read chemcm(i,j,7,12) from netCDF restart
!
!     Purpose
!     -------
!     Read restart data to continue an interrupted integration.
!
!     Method
!     -------
!     The bgc data are read from a file other than the ocean data.
!     The time stamp of the bgc restart file (idate) is specified from the
!     ocean time stamp through the SBR parameter list of AUFW_BGC. The only 
!     time control variable proper to the bgc is the time step number 
!     (idate(5)). It can differ from that of the ocean (idate(4)) by the 
!     difference of the offsets of restart files.
!
!**   Interface.
!     ----------
!
!     *CALL*       *AUFR_BGC(kpie,kpje,kpke,pddpo
!                            ,kplyear,kplmon,kplday,kpldtoce)*
!
!     *COMMON*     *MO_PARAM1_BGC* - declaration of ocean/sediment tracer.
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
!     *INTEGER* *kplyear*   - year  in ocean restart date
!     *INTEGER* *kplmon*  - month in ocean restart date
!     *INTEGER* *kplday*    - day   in ocean restart date
!     *INTEGER* *kpldtoce*  - step  in ocean restart date
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
!ik IK introduced ocetra array elements for phyto, grazer, poc (=det), calciu 
!ik array indices are: iphy, izoo, idet, icalc
!iktodo IK introduced new variable opal (index iopal)
!ik nocetra is the number of all BGC elements (size of ocetra(,,,l))
!ik nosedi is the number of all elements interacting with the sediment

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc 

      use mo_mpi
      use mo_parallel
      
      implicit none
      
      INTEGER  kpie,kpje,kpke
      INTEGER  kplyear,kplmon,kplday,kpldtoce
      REAL pddpo(kpie,kpje,kpke)
      INTEGER  i,j,k,l,kmon

      INTEGER idate(5)

      INTEGER :: restyear            !  year of restart file
      INTEGER :: restmonth           !  month of restart file
      INTEGER :: restday             !  day of restart file
      INTEGER :: restdtoce           !  time step number from bgc ocean file

#ifdef PNETCDF                
      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncstat,ncvarid
!
! Open netCDF data file
!
!js 2 lines from emr deactivated for production runs
!       write(0,333)kpie,kpje,kpke
!333    format('aufrbgc',3i5)
      IF(p_pe==p_io) THEN
        ncstat = NF_OPEN('restartr_bgc.nc',NF_NOWRITE, ncid)
        IF ( ncstat .NE. NF_NOERR ) CALL STOP_ALL('AUFR: Problem with netCDF1')

!
! Read restart data : date
!

        ncstat = NF_GET_ATT_INT(ncid, NF_GLOBAL,'date', idate)
        IF ( ncstat .NE. NF_NOERR )                                      &
       &CALL STOP_ALL('AUFR: Problem reading date of restart file')
      ENDIF
      CALL p_bcast(idate,p_io)
      restyear  = idate(1)
      restmonth = idate(2)
      restday   = idate(3)
      restdtoce = idate(4)
      ldtbgc = idate(5)
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Date of bgc restart file : '
      WRITE(io_stdo_bgc,*) ' year  = ',restyear
      WRITE(io_stdo_bgc,*) ' month = ',restmonth
      WRITE(io_stdo_bgc,*) ' day   = ',restday
      WRITE(io_stdo_bgc,*) ' oce_step = ',restdtoce
      WRITE(io_stdo_bgc,*) ' bgc_step = ',ldtbgc
      WRITE(io_stdo_bgc,*) ' '

!
! Compare with date read from ocean restart file
!
! As the ocean is already in its first step, its counter has 
! gone up one step already for the year and month. The ocean day
! counter is still at its restart date. Therefore:

!      restmonth = restmonth + 1
!
!      IF (restmonth .GT. 12) THEN
!         restmonth=1
!         restyear=restyear+1
!      ENDIF

! Memorize ocean start date :   
      bgcstartyear  = kplyear
      bgcstartmonth = kplmon
      bgcstartday   = kplday

      IF ( kplyear  .NE. restyear  ) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart years in oce/bgc are not the same : '  &
     &   ,kplyear,'/',restyear,' !!!'
      ENDIF

      IF ( kplmon .NE. restmonth ) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart months in oce/bgc are not the same : '   &
     &   ,kplmon,'/',restmonth,' !!!'
!         STOP 'Stop : restart months in oce/bgc are not the same.'
      ENDIF

      IF ( kplday   .NE. restday   ) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart days in oce/bgc are not the same : '   &
     &   ,kplday,'/',restday,' !!!'
!         STOP 'Stop : restart days in oce/bgc are not the same.'
      ENDIF 

      IF ( kpldtoce .NE. ldtbgc   ) THEN
         WRITE(io_stdo_bgc,*)                                       & 
     &   'WARNING: restart step no.  in oce/bgc are not the same : '&
     &   ,kpldtoce,'/',ldtbgc,' !!!'
      ENDIF

!
! Read restart data : ocean aqueous tracer
!                
      CALL read_netcdf_var(ncid,'sco212',ocetra(1,1,1,isco212),kpke)
! comment out for first run (no 13c in restart file)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sco213',ocetra(1,1,1,isco213),kpke)
      CALL read_netcdf_var(ncid,'sco214',ocetra(1,1,1,isco214),kpke)
#endif

      CALL read_netcdf_var(ncid,'alkali',ocetra(1,1,1,ialkali),kpke)
      CALL read_netcdf_var(ncid,'phosph',ocetra(1,1,1,iphosph),kpke)
      CALL read_netcdf_var(ncid,'oxygen',ocetra(1,1,1,ioxygen),kpke)
      CALL read_netcdf_var(ncid,'gasnit',ocetra(1,1,1,igasnit),kpke)
      CALL read_netcdf_var(ncid,'ano3',ocetra(1,1,1,iano3),kpke)
      CALL read_netcdf_var(ncid,'silica',ocetra(1,1,1,isilica),kpke)
      CALL read_netcdf_var(ncid,'doc',ocetra(1,1,1,idoc),kpke)
      CALL read_netcdf_var(ncid,'phyto',ocetra(1,1,1,iphy),kpke)
      CALL read_netcdf_var(ncid,'grazer',ocetra(1,1,1,izoo),kpke)
      CALL read_netcdf_var(ncid,'poc',ocetra(1,1,1,idet),kpke)

! comment out for first run (no 13c in restart file)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'poc13',ocetra(1,1,1,idet13),kpke)
      CALL read_netcdf_var(ncid,'poc14',ocetra(1,1,1,idet14),kpke)
#endif

      CALL read_netcdf_var(ncid,'calciu',ocetra(1,1,1,icalc),kpke)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'calciu13',ocetra(1,1,1,icalc13),kpke)
      CALL read_netcdf_var(ncid,'calciu14',ocetra(1,1,1,icalc14),kpke)
#endif

      CALL read_netcdf_var(ncid,'opal',ocetra(1,1,1,iopal),kpke)
      CALL read_netcdf_var(ncid,'n2o',ocetra(1,1,1,ian2o),kpke)
      CALL read_netcdf_var(ncid,'dms',ocetra(1,1,1,idms),kpke)
      CALL read_netcdf_var(ncid,'fdust',ocetra(1,1,1,ifdust),kpke)
      CALL read_netcdf_var(ncid,'iron',ocetra(1,1,1,iiron),kpke)

#ifdef AGG
      CALL read_netcdf_var(ncid,'snos',ocetra(1,1,1,inos),kpke)
      CALL read_netcdf_var(ncid,'adust',ocetra(1,1,1,iadust),kpke)
#endif /*AGG*/

#ifdef ANTC14
      CALL read_netcdf_var(ncid,'antc14',ocetra(1,1,1,iantc14),kpke)
#endif
#ifdef PCFC
      CALL read_netcdf_var(ncid,'cfc11',ocetra(1,1,1,icfc11),kpke)
      CALL read_netcdf_var(ncid,'cfc12',ocetra(1,1,1,icfc12),kpke)
#endif

!
!Check aquateous restart data for topography
! 

      DO i    =1,kpie
      DO j    =1,kpje
      DO k    =1,kpke      
      DO l    =1,nocetra 
         IF (pddpo(i,j,k) .le. 0.5 ) THEN
            IF (kchck == 1) THEN
               IF ( ocetra(i,j,k,l) .NE. rmasko ) THEN
                  WRITE(io_stdo_bgc,*) 'ocetra not properly masked at :'  &
     &                              ,i,j,k,l,ocetra(i,j,k,l)
               END IF
            END IF
         ELSE
            IF ( ocetra(i,j,k,l) .EQ. rmasko ) THEN
               WRITE(io_stdo_bgc,*) 'land mask values at wet points :'  &
     &                             ,i,j,k,l,ocetra(i,j,k,l)
!               call STOP_ALL('Stop : Restart file with different topography:     &
!     &         land mask values at wet points')
               IF (pddpo(i-1,j,1) .gt.0.5) ocetra(i,j,k,l)=ocetra(i-1,j,1,l)
               IF (pddpo(i+1,j,1) .gt.0.5) ocetra(i,j,k,l)=ocetra(i+1,j,1,l)
               IF (pddpo(i,j+1,1) .gt.0.5) ocetra(i,j,k,l)=ocetra(i,j+1,1,l)
               IF (pddpo(i,j-1,1) .gt.0.5) ocetra(i,j,k,l)=ocetra(i,j-1,1,l)
            END IF
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

!
! Read restart data : other fields
!
     CALL read_netcdf_var(ncid,'hi',hi(1,1,1),kpke)
     CALL read_netcdf_var(ncid,'co3',co3(1,1,1),kpke)
!     CALL read_netcdf_var(ncid,'aksp',aksp(1,1,1),kpke)

!
! Read restart data : chemical constants
!
!!$      CALL read_netcdf_var(ncid,'chemcm',chemcm(1,1,1,1),8*12)
!!$      DO kmon =1,12
!!$      DO i    =1,kpie
!!$      DO j    =1,kpje
!!$       IF (pddpo(i,j,1) .le. 0.5 ) THEN
!!$          DO k=1,8
!!$             IF ( chemcm(i,j,k,kmon) .NE. rmasko ) THEN
!!$             WRITE(io_stdo_bgc,*) 'chemcm,i,j,k,kmon:',          &
!!$     &                          chemcm(i,j,k,kmon),i,j,k,kmon
!!$                chemcm(i,j,k,kmon)  =   rmasko
!!$             ENDIF
!!$          ENDDO
!!$       ELSE
!!$         chemcm(i,j,1,kmon) = chemcm(i,j,1,kmon)*1.e-09
!!$         chemcm(i,j,2,kmon) = chemcm(i,j,2,kmon)*1.e-13
!!$         chemcm(i,j,3,kmon) = chemcm(i,j,3,kmon)*1.e-09
!!$         chemcm(i,j,4,kmon) = chemcm(i,j,4,kmon)*1.e-06
!!$         chemcm(i,j,5,kmon) = chemcm(i,j,5,kmon)*1.e-08
!!$         chemcm(i,j,6,kmon) = chemcm(i,j,6,kmon)*1.e-04
!!$         chemcm(i,j,7,kmon) = chemcm(i,j,7,kmon)*1.e-04
!!$         chemcm(i,j,8,kmon) = chemcm(i,j,8,kmon)*1.e-04
!!$       ENDIF
!!$      ENDDO
!!$      ENDDO
!!$      ENDDO
!!$!
! Read restart data : sediment variables.
! js: reading of terrigenous sediment was missing until 02.08.2005
      CALL read_netcdf_var(ncid,'ssso12',sedlay(1,1,1,issso12),ks)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'ssso13',sedlay(1,1,1,issso13),ks)
      CALL read_netcdf_var(ncid,'ssso14',sedlay(1,1,1,issso14),ks)
#endif

      CALL read_netcdf_var(ncid,'sssc12',sedlay(1,1,1,isssc12),ks)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sssc13',sedlay(1,1,1,isssc13),ks)
      CALL read_netcdf_var(ncid,'sssc14',sedlay(1,1,1,isssc14),ks)
#endif

      CALL read_netcdf_var(ncid,'ssssil',sedlay(1,1,1,issssil),ks)
      CALL read_netcdf_var(ncid,'ssster',sedlay(1,1,1,issster),ks)
      CALL read_netcdf_var(ncid,'bur_o12',burial(1,1,issso12),1)
      CALL read_netcdf_var(ncid,'bur_c12',burial(1,1,isssc12),1)
      CALL read_netcdf_var(ncid,'bur_sil',burial(1,1,issssil),1)
      CALL read_netcdf_var(ncid,'bur_clay',burial(1,1,issster),1)
! js: need to add 13C 14C in burial layer (also aufw_bgc.f90)
      CALL read_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks)
      CALL read_netcdf_var(ncid,'powaic',powtra(1,1,1,ipowaic),ks)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'powc13',powtra(1,1,1,ipowc13),ks)
      CALL read_netcdf_var(ncid,'powc14',powtra(1,1,1,ipowc14),ks)
#endif

      CALL read_netcdf_var(ncid,'powaal',powtra(1,1,1,ipowaal),ks)
      CALL read_netcdf_var(ncid,'powaph',powtra(1,1,1,ipowaph),ks)
      CALL read_netcdf_var(ncid,'powaox',powtra(1,1,1,ipowaox),ks)
      CALL read_netcdf_var(ncid,'pown2',powtra(1,1,1,ipown2),ks)
      CALL read_netcdf_var(ncid,'powno3',powtra(1,1,1,ipowno3),ks)
      CALL read_netcdf_var(ncid,'powasi',powtra(1,1,1,ipowasi),ks)

#ifdef DIFFAT 
!
! Read restart data : co2 diffusion
!
      CALL read_netcdf_var(ncid,'atmco2',atm(1,1,iatmco2),1)
      CALL read_netcdf_var(ncid,'atmo2',atm(1,1,iatmo2),1)
      CALL read_netcdf_var(ncid,'atmn2',atm(1,1,iatmn2),1)
#endif

      if(p_pe==p_io) ncstat = NF_CLOSE(ncid)

#else  !! Attention !! Update required !!!!!

     ! Attention - this doesn't work any more for MPI
     ! One big FORTRAN READ like here is a very, very bad idea
     ! for parallel programs !!!!!

#ifndef NOMPI
     call stop_all("we can't read in bgc restart data with MPI and no NETCDF now")
! would need to read these in on p_io and scatter if we were actually i
! going to do this - c.f. read_netcdf_var.F90
#else

      OPEN(io_rsti_bgc,FILE='restart_bgc',STATUS='UNKNOWN'              &
     &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      READ(io_rsti_bgc)                                                 &
     &            (((ocetra(i,j,k,iphosph),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocetra(i,j,k,isilica),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocetra(i,j,k,ioxygen),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocetra(i,j,k,iphy   ),i=1,kpie),j=1,kpje),k=1,kpke)& 
     &           ,(((ocetra(i,j,k,izoo   ),i=1,kpie),j=1,kpje),k=1,kpke)& 
     &           ,(((ocetra(i,j,k,idet   ),i=1,kpie),j=1,kpje),k=1,kpke)& 
     &           ,(((ocetra(i,j,k,idoc   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocetra(i,j,k,icalc  ),i=1,kpie),j=1,kpje),k=1,kpke)& 
     &           ,(((ocetra(i,j,k,iano3  ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocetra(i,j,k,igasnit),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocetra(i,j,k,iopal  ),i=1,kpie),j=1,kpje),k=1,kpke)


      READ(io_rsti_bgc) chemcm,hi,co3,aksp                              &
     &           ,(((ocetra(i,j,k,isco212),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocetra(i,j,k,ialkali),i=1,kpie),j=1,kpje),k=1,kpke)

      READ(io_rsti_bgc) sedlay,sedhpl

      READ(io_rsti_bgc)                                                 &
     &            (((powtra(i,j,k,ipowaic),i=1,kpie),j=1,kpje),k=1,ks)  &
     &           ,(((powtra(i,j,k,ipowaal),i=1,kpie),j=1,kpje),k=1,ks)  &
     &           ,(((powtra(i,j,k,ipowaph),i=1,kpie),j=1,kpje),k=1,ks)  &
     &           ,(((powtra(i,j,k,ipowaox),i=1,kpie),j=1,kpje),k=1,ks)  &
     &           ,(((powtra(i,j,k,ipowasi),i=1,kpie),j=1,kpje),k=1,ks)  &
     &           ,(((powtra(i,j,k,ipowno3),i=1,kpie),j=1,kpje),k=1,ks)  &
     &           ,(((powtra(i,j,k,ipown2) ,i=1,kpie),j=1,kpje),k=1,ks)

      CLOSE (io_rsti_bgc)
#endif

#endif

!     
!  Masking aqueous sea water tracer.
!
      DO l=1,nocetra
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(pddpo(i,j,k) .LT. 0.5) THEN
         ocetra(i,j,k,l)=rmasko
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
! construction site!!!!!!----------------------------------------------------------<<<<<
!     DO k=1,kpke
!     DO j=1,kpje
!     DO i=1,kpie
! below line only for initialisation from scratch       
!     ocetra(i,j,k,isco214)=ocetra(i,j,k,isco212)*(0.96-0.004*k)
! insert 14C correction here?? 
! shifted to ini_bgc (dlxp, dlyp not in aufr)

! initialisation of 13C using phosphate (use only once! after introduction of 13C)
!aufsetz! (comment out if not)
!     ocetra(i,j,k,isco213)=ocetra(i,j,k,isco212)*                    &
!    &                  (1+0.001*(2.7 - 1.1e6*ocetra(i,j,k,iphosph)))

! set poc pool to zero
!     if (ocetra(i,j,k,idet).ne.rmasko) ocetra(i,j,k,idet)=0.


!     ENDDO
!     ENDDO
!     ENDDO      
!
!  Masking other sea water tracers.
!
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(pddpo(i,j,k) .LT. 0.5) THEN
!ik      phyto(i,j,k)=rmasko
!ik      grazer(i,j,k)=rmasko
!ik      poc(i,j,k)=rmasko
         hi(i,j,k)=rmasko
!ik      calciu(i,j,k)=rmasko
         co3(i,j,k)=rmasko
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!
!  Masking sediment pore water tracers.
!
      DO  k=1,ks
      DO  j=1,kpje
      DO  i=1,kpie 
      IF(kbo(i,j) .LE. 1) THEN
         powtra(i,j,k,ipowaic)=rmasks
#ifdef __c_isotopes
         powtra(i,j,k,ipowc13)=rmasks
         powtra(i,j,k,ipowc14)=rmasks
#endif
         powtra(i,j,k,ipowaal)=rmasks
         powtra(i,j,k,ipowaph)=rmasks
         powtra(i,j,k,ipowaox)=rmasks
         powtra(i,j,k,ipown2)=rmasks
         powtra(i,j,k,ipowno3)=rmasks
         powtra(i,j,k,ipowasi)=rmasks
         sedlay(i,j,k,issso12)=rmasks
#ifdef __c_isotopes
         sedlay(i,j,k,issso13)=rmasks
         sedlay(i,j,k,issso14)=rmasks
#endif
         sedlay(i,j,k,isssc12)=rmasks
#ifdef __c_isotopes
         sedlay(i,j,k,isssc13)=rmasks
         sedlay(i,j,k,isssc14)=rmasks
#endif
         sedlay(i,j,k,issssil)=rmasks
         sedlay(i,j,k,issster)=rmasks
       burial(i,j,issso12)=rmasks
       burial(i,j,isssc12)=rmasks
       burial(i,j,issssil)=rmasks
         burial(i,j,issster)=rmasks
         sedhpl(i,j,k)=rmasks
      else
!aufsetz! comment out if not
!        sedlay(i,j,k,issster)=min(sedlay(i,j,k,issster),1000.)  ! only for initialization!!!!!!!!!
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      return
!     
!  Restrict to positive values (until error is found only !!!!) js: which error?
!
      DO l=1,nocetra
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(pddpo(i,j,k) .GT. 0.5) THEN
       ocetra(i,j,k,l)=MAX(ocetra(i,j,k,l),0.)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(pddpo(i,j,k) .GT. 0.5) THEN
!ik         phyto (i,j,k)=MAX(phyto (i,j,k),0.)
!ik         grazer(i,j,k)=MAX(grazer(i,j,k),0.)
!ik         poc   (i,j,k)=MAX(poc   (i,j,k),0.)
         hi    (i,j,k)=MAX(hi    (i,j,k),1.e-12)
!ik         calciu(i,j,k)=MAX(calciu(i,j,k),0.)
      ENDIF
      ENDDO
      ENDDO
      ENDDO


      DO  k=1,ks
      DO  j=1,kpje
      DO  i=1,kpie 
      IF(kbo(i,j) .GT. 1) THEN
         powtra(i,j,k,ipowaic)=MAX(powtra(i,j,k,ipowaic),0.)
#ifdef __c_isotopes
         powtra(i,j,k,ipowc13)=MAX(powtra(i,j,k,ipowc13),0.)
         powtra(i,j,k,ipowc14)=MAX(powtra(i,j,k,ipowc14),0.)
#endif
         powtra(i,j,k,ipowaal)=MAX(powtra(i,j,k,ipowaal),0.)
         powtra(i,j,k,ipowaph)=MAX(powtra(i,j,k,ipowaph),0.)
         powtra(i,j,k,ipowaox)=MAX(powtra(i,j,k,ipowaox),0.)
         powtra(i,j,k,ipown2) =MAX(powtra(i,j,k,ipown2) ,0.)
         powtra(i,j,k,ipowno3)=MAX(powtra(i,j,k,ipowno3),0.)
         powtra(i,j,k,ipowasi)=MAX(powtra(i,j,k,ipowasi),0.)
         sedlay(i,j,k,issso12)=MAX(sedlay(i,j,k,issso12),0.)
         sedlay(i,j,k,isssc12)=MAX(sedlay(i,j,k,isssc12),0.)
#ifdef __c_isotopes
         sedlay(i,j,k,isssc13)=MAX(sedlay(i,j,k,isssc13),0.)
         sedlay(i,j,k,isssc14)=MAX(sedlay(i,j,k,isssc14),0.)
#endif
         sedlay(i,j,k,issssil)=MAX(sedlay(i,j,k,issssil),0.)
         sedlay(i,j,k,issster)=MAX(sedlay(i,j,k,issster),0.)
         sedhpl(i,j,k)        =MAX(sedhpl(i,j,k)    ,1.e-12)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

! aufsetz? (new ernst)
!     call sedshi(kpie,kpje)
!     call maschk(kpie,kpje,kpke,100)

      RETURN
      END
