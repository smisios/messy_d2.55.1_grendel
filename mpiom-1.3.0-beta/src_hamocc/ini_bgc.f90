
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

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
      use mo_param1_bgc 
      use mo_mpi
#ifdef PDYNAMIC_BGC
      use mo_dynamic
#endif /* PDYNAMIC_BGC */    

      implicit none
      INTEGER :: pyears,pmonts,kpie,kpje,kpke
      INTEGER :: kplyear,kplmonth,kplday,kpldtoce
      INTEGER :: kpaufr,kpicycli,kpndtrun,k
      INTEGER :: ii,jj,kk
      
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: psao (kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: pgila(kpie*2,kpje*2)
      REAL :: pgiph(kpie*2,kpje*2)
      REAL :: ptiestu(kpke+1),ptiestw(kpke+1)

      REAL :: pdt

#ifdef DMSASSIM
      REAL :: dmsrms_ref
      REAL :: assim_fac(2)
      INTEGER :: ndmsdr,idmsdr,i,dmssuc(5),dmsdir(5)
#endif      
      
      !                    
! Set control constants ( mo_control_bgc )
!
      dtbgc = pdt                   !  time step length [sec]
      ndtdaybgc=NINT(86400./dtbgc)  !  time steps per day [no.]
      dtb=1./ndtdaybgc              !  time step length [days]
      
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
!
      CALL READ_NAMELIST

!                        
! Initialize bgc time series.
!
      CALL INI_TIMESER_BGC(kpke,ptiestw)

!                        
! Initialize bgcmean
!
      
      n90depth   = 1
      n1000depth = 1
      n2000depth = 1
      DO k = 1, kpke
        if (ptiestw(k).lt.90)   n90depth   = k
        if (ptiestw(k).lt.1000) n1000depth = k
        if (ptiestw(k).lt.2000) n2000depth = k
      ENDDO
!e    io_stdo_bgc=0
      write(io_stdo_bgc,*)'fluxes will be sampled at levels :', n90depth, n1000depth, n2000depth
!e    write(0,322)io_stdo_bgc
  322 format('filezuordnung',i5)
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
!!$      meantime_3d  = 1
!!$      endif     

      if (mean_2D_freq.eq.1) nmeantime_2d = ndtrunbgc*dtb
      if (mean_3D_freq.eq.1) nmeantime_3d = ndtrunbgc*dtb
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
! Read restart fields from restart file
!
      IF(kpaufr.eq.1) THEN
         CALL AUFR_BGC(kpie,kpje,kpke,pddpo,kplyear,kplmonth,   &
     &                 kplday,kpldtoce)
      ENDIF

! aufsetz! (for initialization of 14C)
!     call c14_correction(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,psao,ptho)

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


#ifdef DMSASSIM
       
!       ndmsdr = 5
!       idmsdr = 1
!       dmsrms_ref = 1.e99
!       assim_fac(1) = 0.9
!       assim_fac(2) = 1.1
!       dmssuc       = 1
!       dmsdir       = 1
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
   
      
!#ifdef DIFFAT
! correction of alkalinity during spin-up of pco2 
!
! temporary test js

      if (lspinbgc) then
      CALL SPINUP_BGC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp)
      endif
!#endif

!
! Global inventory of all tracers
!

      CALL INVENTORY_BGC(kpie,kpje,kpke)

!     CALL CHCK_BGC(io_stdo_bgc,kpicycli,                         &
!    &'Check values of ocean tracer at exit from SBR INI_BGC :',  &
!    & kpie,kpje,kpke,pddpo) 

      RETURN
      END
