
      SUBROUTINE INI_BGC(kpaufr,kpicycli,pdt,kpndtrun,kpie,kpje,kpke   &
     &           ,pddpo,ptho,psao,pdlxp,pdlyp,ptiestu,ptiestw          &
     &           ,kplyear,kplmonth,kplday,kpldtoce,pyears,pmonts       &
     &           ,pgila,pgiph)
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/ini_bgc.f90,v $\\
!$Revision: 1.2.2.1.20.2 $\\
!$Date: 2006/04/03 11:27:49 $\\

!****************************************************************
!
!**** *INI_BGC* - initialize marine bio-geo-chemistry module.
!
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     - initialize sediment layering
!     - initialize bgc variables
!     - calculation of chemical constants
!     - read restart fields
!     - calculate budgets
!
!     Method
!     -------
!     - Noch zu tun : Biharmonic tracer diffusion (mit AULAPTS > ALMZER)
!                     Convective adjustment
!                     sea ice growth influence on tracer concentration
!
!**   Interface.
!     ----------
!
!     *CALL*       *INI_BGC(kpaufr,kpicycli,pdt,kpndtrun,kpie,kpje,kpke
!                          ,pddpo,ptho,psao,pdlxp,pdlyp,ptiestu
!                          ,kplyear,kplmonth,kplday,kpldtoce)*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpaufr*   - 1/0 for read / do not read restart file
!     *INTEGER* *kpicycli* - flag for cyclicity.
!     *REAL*    *pdt*      - ocean model time step [sec].
!     *INTEGER* *kpndtrun* - total no. of time steps of run.
!     *INTEGER* *kpie*     - zonal dimension of model grid.
!     *INTEGER* *kpje*     - meridional dimension of model grid.
!     *INTEGER* *kpke*     - vertical dimension of model grid.
!     *REAL*    *pddpo*    - size of scalar grid cell (3rd REAL) [m].
!     *REAL*    *ptho*     - potential temperature [deg C].
!     *REAL*    *psao*     - salinity [psu.].
!     *REAL*    *pdlxp*    - size of scalar grid cell (zonal) [m].
!     *REAL*    *pdlyp*    - size of scalar grid cell (meridional) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu*  - depth of level [m].
!     *INTEGER* *kplyear*  - year  in ocean restart date
!     *INTEGER* *kplmonth* - month in ocean restart date
!     *INTEGER* *kplday*   - day   in ocean restart date
!     *INTEGER* *kpldtoce* - step  in ocean restart date
!
!     Externals
!     ---------
!     READ_NAMELIST CALL INI_TIMESER_BGC BODENSED BELEG_BGC AUFR_BGC
!     CHEMIN CHCK_BGC PRINT_BGC
!
!**********************************************************************

      USE mo_kind, ONLY: dp
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
      use mo_param1_bgc
      use mo_bgc_diagnostic
      use mo_mpi
      USE mo_parallel, ONLY: stop_all
#ifdef PDYNAMIC_BGC
      use mo_dynamic
#endif /* PDYNAMIC_BGC */
#ifdef AVFLUX
      use mo_avflux
#endif
#ifndef NO_NEW_IO
      USE mo_bgc_varlist,ONLY : build_bgc_varlist, bgc_varlist
      USE mo_bgc_iolist,ONLY : read_namelist_bgc, build_iolist_bgc, bgc_iolist
      USE mo_iolist,ONLY : readRestartCDI
#endif
      USE mo_bgc_diagnostic, ONLY: bgc_diagnostic_init
      USE mo_grid, ONLY: get_level_index_by_depth

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif


      implicit none
      INTEGER :: pyears,pmonts,kpie,kpje,kpke
      INTEGER :: kplyear,kplmonth,kplday,kpldtoce
      INTEGER :: kpaufr,kpicycli,kpndtrun

      REAL(wp) :: pddpo(kpie,kpje,kpke)
      REAL(wp) :: ptho (kpie,kpje,kpke)
      REAL(wp) :: psao (kpie,kpje,kpke)
      REAL(wp) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL(wp) :: pgila(kpie*2,kpje*2)
      REAL(wp) :: pgiph(kpie*2,kpje*2)
      REAL(wp) :: ptiestu(kpke+1),ptiestw(kpke+1)

      REAL(wp) :: pdt

#ifdef DMSASSIM
      REAL(wp) :: dmsrms_ref
      REAL(wp) :: assim_fac(2)
      INTEGER :: ndmsdr,idmsdr,i,dmssuc(5),dmsdir(5)
#endif
      INTEGER :: ierror

      !
! Set control constants ( mo_control_bgc )
!
      dtbgc = pdt                   !  time step length [sec]
      ndtdaybgc = NINT(86400._wp / dtbgc)  !  time steps per day [no.]
      dtb = 1._wp / REAL(ndtdaybgc, wp)        !  time step length [days]

      icyclibgc = kpicycli
      ndtrunbgc = kpndtrun

!
! Initialize some namelist parameters
!
      isac = 1

      mean_2D_freq = 2
      mean_3D_freq = 3
!
! Initialize time step counter of run.
!
      ldtrunbgc = 0

!
! Set namelist defaults and read namelist

      CALL READ_NAMELIST

#ifndef NO_NEW_IO
      ierror = 0
      CALL read_namelist_bgc(ierror)

      IF ( ierror .NE. 0 ) THEN
        IF (p_pe == p_io) CALL stop_all('error in reading namelist - run aborted' )
      ELSE
        IF (p_pe == p_io) WRITE(0,*)'read BGCCTL successfully'
      ENDIF
#endif/*ndef NO_NEW_IO */

!
! Initialize global settings.
!
      CALL param1_bgc_init
!
! Initialize bgc time series.
!
      CALL INI_TIMESER_BGC(kpke,ptiestw)

!
! Initialize bgcmean
!
      n90depth   = get_level_index_by_depth(90._dp)
      n1000depth = get_level_index_by_depth(1000._dp)
      n2000depth = get_level_index_by_depth(2000._dp)
!e    io_stdo_bgc=0
      write(io_stdo_bgc,*)'fluxes will be sampled at levels :', n90depth, n1000depth, n2000depth
!e    write(0,322)io_stdo_bgc
!  322 format('filezuordnung',i5)
      meantime_2d = 0
      meantime_3d = 0

!!$      if (mean_2D_day.eq.1) then
!!$        nmeantime_2d = ndtrunbgc*dtb
!!$      else
!!$        nmeantime_2d = 12*pyears + pmonts
!!$      endif
!!$
!!$      if (mean_3D_month.eq.1) then
!!$        nmeantime_3d = 12*pyears + pmonts
!!$      else
!!$        nmeantime_3d = 1
!!$     meantime_3d  = 1
!!$      endif

      IF (mean_2D_freq .EQ. 1) nmeantime_2d = INT(REAL(ndtrunbgc, wp) * dtb)
      IF (mean_3D_freq .EQ. 1) nmeantime_3d = INT(REAL(ndtrunbgc, wp) * dtb)
      if (mean_2D_freq.eq.2) nmeantime_2d = 12*pyears + pmonts
      if (mean_3D_freq.eq.2) nmeantime_3d = 12*pyears + pmonts
      if (mean_2D_freq.eq.3) nmeantime_2d = 1
      if (mean_3D_freq.eq.3) nmeantime_3d = 1

      if (mean_2D_freq.eq.4) nmeantime_2d = 1
      if (mean_3D_freq.eq.4) nmeantime_3d = 1


      CALL ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

!
! Allocate memory : biology
!
      CALL ALLOC_MEM_BIOMOD(kpie,kpje,kpke)

!
! Allocate memory : sediment
!
      CALL ALLOC_MEM_SEDMNT(kpie,kpje)

!
! Allocate memory : inorganic carbon cycle
!
      CALL ALLOC_MEM_CARBCH(kpie,kpje,kpke)

#ifdef PDYNAMIC_BGC
      CALL ALLOC_MEM_DYNAMIC(kpie,kpje,kpke)
#endif /* PDYNAMIC_BGC */

!
! Initialize diagnostic variables
!
      CALL bgc_diagnostic_init(io_stdo_bgc, kpie, kpje, kpke, kwrbioz)

!
! Initialize sediment layering
!
      CALL BODENSED(kpie,kpje,kpke,pddpo)

!
! Initialize sediment and ocean tracers.
!
      CALL BELEG_BGC(kpie,kpje,kpke,pddpo)

!
! Initialize chemical constants: this routine needs 3D temp./sal.
! input from the ocean model. For MPI-OM this means, that it must not be
! called before SBR AUFR.
!
! If kpaufr.eq.1 the initial chem. constants are read from the restart file.
! If kpaufr.eq.0 they have to be calculated now. !js is this still so?
! Parameter "-13" is initializing all months.

      CALL CHEMCON(-13,kpie,kpje,kpke,psao,ptho,                &
     &     pddpo,pdlxp,pdlyp,ptiestu,kplmonth)

!

#ifndef NO_NEW_IO
!      WRITE(0,*) 'before build_varlist'
      CALL build_bgc_varlist
!      WRITE(0,*) 'before build_iolist'
      CALL build_iolist_bgc
#endif/*ndef NO_NEW_IO */

! Read restart fields from restart file
!
      IF(kpaufr.eq.1) THEN
#ifdef _PROFILE
         CALL trace_start ('restartreadbgc', 12)
#endif

#ifdef OLD_IO
         ! Old I/O is overridden by new I/O, but we chose to run both as a test.
         CALL AUFR_BGC(kpie,kpje,kpke,pddpo,kplyear,kplmonth,   &
     &                 kplday,kpldtoce)
#endif
#ifndef NO_NEW_IO
         ! May override old I/O (see above).
         CALL readRestartCDI(bgc_iolist,bgc_varlist)
         ! Old I/O is overridden by new I/O, but we chose to run both as a test.
         CALL AUFR_MASK_BGC(kpie,kpje,kpke,pddpo,kplyear,kplmonth,   &
     &                 kplday,kpldtoce)
#endif/*ndef NO_NEW_IO */

#ifdef _PROFILE
         CALL trace_stop ('restartreadbgc', 12)
#endif
      ENDIF

! aufsetz! (for initialization of 14C)
!     call c14_correction(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,psao,ptho)

#ifdef OLD_IO
! Open bgc output files
      CALL OPEN_BGCMEAN_2D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,pgila,pgiph     &
     &                    ,ptiestu,kplyear,kplmonth,kplday,kpldtoce)
      CALL OPEN_BGCMEAN_BIOZ(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,pgila,pgiph   &
     &                    ,ptiestu,kplyear,kplmonth,kplday,kpldtoce)
      CALL OPEN_BGCMEAN_3D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,pgila,pgiph     &
     &                    ,ptiestu,kplyear,kplmonth,kplday,kpldtoce)
      CALL OPEN_BGCMEAN_SED(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,pgila,pgiph    &
     &                    ,ptiestu,kplyear,kplmonth,kplday,kpldtoce)
#ifdef PDYNAMIC_BGC
      CALL OPEN_DYNAMIC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,pgila,pgiph     &
     &                    ,ptiestu,kplyear,kplmonth,kplday,kpldtoce)
#endif /* PDYNAMIC_BGC */
#endif/*def OLD_IO */


! store starting values to calculated rate of change of various tracers within euphotic zone
! CMIP5 requirement

        call store_tracer(pddpo)

#ifdef DMSASSIM

!       ndmsdr = 5
!       idmsdr = 1
!       dmsrms_ref = 1.e99
!       assim_fac(1) = 0.9
!       assim_fac(2) = 1.1
!       dmssuc   = 1
!       dmsdir   = 1
!
!       OPEN(321,FILE='SUCHWURZ',ACCESS='SEQUENTIAL',FORM='FORMATTED')
!       write(321,3210) ndmsdr,idmsdr,dmssuc,dmsdir,dmsrms_ref,assim_fac,dmspar
!       CLOSE(321)

     IF(p_pe == p_io) THEN

       OPEN(321,FILE='SUCHWURZ',ACCESS='SEQUENTIAL',FORM='FORMATTED')
       do i=1,100000
       read(321,3210,err=1466,end=1466) ndmsdr,idmsdr,dmssuc,dmsdir,      &
     &                                  dmsrms_ref,assim_fac,dmspar
       write(io_stdo_bgc,*)'INI_BGC :',i,ndmsdr,idmsdr,dmssuc,dmsdir,      &
     &                                  dmsrms_ref,assim_fac,dmspar
       enddo

1466   continue
3210  format(12i4,8e15.9)
      CLOSE(321)

     ENDIF
!
!      call DMS_ASSIM(kpie,kpje)
!
#endif  /*DMSASSIM*/


      if (lspinbgc) then
      CALL SPINUP_BGC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp)
      endif

!
! Global inventory of all tracers
!

      CALL INVENTORY_BGC(kpie,kpje,kpke)

!     CALL CHCK_BGC(io_stdo_bgc,kpicycli,                         &
!    &'Check values of ocean tracer at exit from SBR INI_BGC :',  &
!    & kpie,kpje,kpke,pddpo)


#ifdef __cpl_co2
#ifdef AVFLUX
!     initializes fields used for redistribution of co2 fluxes
      call avflux_ini
#endif
#endif


      RETURN
      END
