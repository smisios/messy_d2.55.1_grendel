MODULE MESSY_SURFACE

  
  USE messy_main_constants_mem,     ONLY: DP,g,alf,        &
                                          cpd=>cp_air,     &
                                          vtmpc2, tmelt,   &
                                          rhoh2o=>rho_H2O, &
                                          cwlmax

  IMPLICIT NONE
  PRIVATE
  SAVE

  CHARACTER(len=*), PUBLIC, PARAMETER :: modstr = 'surface'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver = '1.1'

  ! efficency of interception of precipitation as rain
  Real(dp), Parameter :: cvinter = 0.25_dp  ! see mo_vegetation and iniphy

  ! hard coded number of soil layers
  INTEGER,  Parameter :: jpgrnd  = 5        ! see mo_parameters


  ! control namelist parameters now in surface.nml, see mo_param_switches
  LOGICAL, PUBLIC :: lice 

  ! maximum moisture content of the skin reservoir
      
  PUBLIC  :: surface_read_nml_ctrl
  PUBLIC  :: surface
  PUBLIC  :: lake
! PRIVATE :: soiltemp
  PUBLIC  :: licetemp
  PUBLIC  :: sicetemp 

CONTAINS

!===============================================================================

  SUBROUTINE surface_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

   ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status

    NAMELIST /CTRL/ lice                  ! switches now in surface.nml
                                                 ! former mo_param_switches
   ! LOCAL USE 
    CHARACTER(LEN=*), PARAMETER :: substr = 'surface_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    lice   = .false.

    status = 1            ! error

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR
  
  END SUBROUTINE surface_read_nml_ctrl

! *******************************************************************

! *******************************************************************

  SUBROUTINE surface ( kproma, kbdim, klev          &
                      ,ngl, zdtime, eps, lstart     &
!
! - 1D from mo_memory_g3b
          , ptsl,       ptslm,     ptslm1           &
          , pws,        pwl,       pwsmx            &
          , psn,        psnmel,    pgld             &
          , psnc,       pu10,      pv10             &
          , prunoff,    progl,     pdrain           &
          , papmegl,    psnacl,    porostd          &
          , prgcgn,     psodif,    pslm             &
          , pgrndcapc,  pgrndhflx, pgrndflux        &
! - 2D from mo_memory_g3b (soil variables)
          , ptsoil,     pgrndd,    pgrndc           &
          , ptm1,       pqm1,      ztte             &
          , paphm1                                  &
! - variables internal to physics
          , pcvs,       pcvw,      pwlmx            &
          , pevapl,     pevapot                     &
          , prsfl,      prsfc                       &
          , pssfl,      pssfc                       &
          , pros_hd,    pdrain_hd                   &
          , palac                                   &
          , lpland,     lpglac                      &
          , zevwsd                                  )  !op_re_20140904

 IMPLICIT NONE

     INTEGER,  INTENT(IN)    :: kproma
     INTEGER,  INTENT(IN)    :: kbdim
     INTEGER,  INTENT(IN)    :: klev
     INTEGER,  INTENT(IN)    :: ngl
     REAL(dp), INTENT(IN)    :: zdtime
     REAL(dp), INTENT(IN)    :: eps
     LOGICAL,  INTENT(IN)    :: lstart

     REAL(dp), INTENT(INOUT) :: ptsl(:)      ! Sfc temp [K] at timestep t+dt (unfiltered)
     REAL(dp), INTENT(INOUT) :: ptslm(:)     ! Sfc temp [K] at timestep t (unfiltered)
     REAL(dp), INTENT(INOUT) :: ptslm1(:)    ! Sfc temp [K] at timestep t-dt (filtered)
     REAL(dp), INTENT(INOUT) :: pws(:)       ! Soil water content [m]
     REAL(dp), INTENT(INOUT) :: pwl(:)       ! Water content [m] in skin reservoir (vegetation and bare soil)
     REAL(dp), INTENT(IN)    :: pwsmx(:)     ! Water holding capacity [m] of the soil
     REAL(dp), INTENT(INOUT) :: psn(:)       ! Snow depth [m water equivalent] at the ground
     REAL(dp), INTENT(INOUT) :: psnmel(:)    ! Snow melt [kg/m**2] (accumulated for diagnostics)    
     REAL(dp), INTENT(INOUT) :: pgld(:)      ! Glacier depth (including snow) [m water equivalent]
     REAL(dp), INTENT(INOUT) :: psnc(:)      ! Snow depth [m water equivalent] at the canopy
     REAL(dp), INTENT(IN)    :: pu10(:)      ! Wind (u-component) at 10m height [m/s] ... from 'vdiff'
     REAL(dp), INTENT(IN)    :: pv10(:)      ! Wind (v-component) at 10m height [m/s] ... from 'vdiff'
     REAL(dp), INTENT(INOUT) :: prunoff(:)   ! Total runoff [kg/m**2] at non-glacier points (accumul.)
     REAL(dp), INTENT(INOUT) :: progl(:)     ! Glacier runoff (rain+snow/ice melt) [kg/m**2] (accumul.)
     REAL(dp), INTENT(INOUT) :: pdrain(:)    ! Drainage at non-glacier points [kg/m**2] (accumul.)
     REAL(dp), INTENT(INOUT) :: papmegl(:)   ! Precip-Evap at glacier points [kg/m**2] (accumul.)
     REAL(dp), INTENT(INOUT) :: psnacl(:)    ! Snow budget at non-glacier points [kg/m**2] (accumul.)
     REAL(dp), INTENT(IN)    :: porostd(:)   ! Subgrid standard diviation [m] used in runoff scheme
     REAL(dp), INTENT(IN)    :: prgcgn(:)    ! Volumetric heat capacity of the soil [j/m**3/K]
                                             ! set in ionitial
     REAL(dp), INTENT(IN)    :: psodif(:)    ! Thermal diffusivity of the soil [m**2/s]
     REAL(dp), INTENT(IN)    :: pslm(:)      ! Land sea mask
     REAL(dp), INTENT(INOUT) :: pgrndcapc(:) ! Heat capacity of the uppermost soil layer [j/m**2/K]
     REAL(dp), INTENT(INOUT) :: pgrndhflx(:) ! Soil heat flux at the surface [W/m**2]
     REAL(dp), INTENT(INOUT) :: pgrndflux(:)
   
     REAL(dp), INTENT(INOUT) :: ptsoil(:,:)  ! Temperature [K] in the five soil layers
     REAL(dp), INTENT(INOUT) :: pgrndc(:,:)    ! Coefficients used in ptsoil calculation (*soiltemp*)
     REAL(dp), INTENT(INOUT) :: pgrndd(:,:)    ! Coefficients used in ptsoil calculation (*soiltemp*)
     REAL(dp), INTENT(IN)    :: ptm1(:,:)    ! Air temp [K] at timestep t-dt (filtered)
     REAL(dp), INTENT(IN)    :: pqm1(:,:) 
     REAL(dp), INTENT(OUT)   :: ztte(:,:)
     REAL(dp), INTENT(IN)    :: paphm1(:,:)

     REAL(dp), INTENT(IN)    :: pcvs(:)      ! Fractional snow cover (function of psn in *physc*)
     REAL(dp), INTENT(IN)    :: pcvw(:)      ! Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
     REAL(dp), INTENT(IN)    :: pwlmx(:)     ! Skin reservoir [m] (calculated in *vdiff* as a function
!                                                of vegetation index and leaf area index)
     REAL(dp), INTENT(IN)    :: pevapl(:)    ! Total evaporation, including sublimation [kg/m**2/s]
     REAL(dp), INTENT(IN)    :: pevapot(:)   ! Potential evaporation/sublimation [kg/m**2/s]
     REAL(dp), INTENT(IN)    :: prsfl(:)     ! Large scale rainfall [kg/m**2/s]
     REAL(dp), INTENT(IN)    :: prsfc(:)     ! Convective rainfall [kg/m**2/s]
     REAL(dp), INTENT(IN)    :: pssfl(:)     ! Large scale snowfall [kg/m**2/s]
     REAL(dp), INTENT(IN)    :: pssfc(:)     ! Convective snowfall [kg/m**2/s]
     REAL(dp), INTENT(OUT)   :: pros_hd(:)   ! Runoff for HD-Model (does NOT include drainage) [m]
     REAL(dp), INTENT(OUT)   :: pdrain_hd(:) ! Drainage for HD-Model [m]
     REAL(dp), INTENT(OUT)   :: palac(:)     ! Precipitation minus sublimation at glacier points
     LOGICAL,  INTENT(IN)    :: lpland(:)    ! logical mask land
     LOGICAL,  INTENT(IN)    :: lpglac(:)    ! logical mask glacier
     REAL(dp), INTENT(INOUT) :: zevwsd(:)    ! New pointer needed for h2oiso

! The following local variables represent the respective fluxes
! integrated over one timestep (delta_time) and divided by the density of
! water (rhoh2o). Units are m water equivalent.
!
! zraind        Total rain
! zsnowd        Total snow
! zevttd        Total evaporation
! zevsnd        Sublimation from snow
! zevwld        Evaporation from the skin reservoir
! zevwsd        Evaporation from the soil and from the skin reservoir
! zros          Total runoff (including drainage) at non-glacier points
! zdrain        Drainage at non-glacier points
! zrogl         Runoff at glacier points (rain and melt, but no calving)
! zsnmel        Snow/ice melt at land and glacier points
! zsncmelt      Snow melt in the canopy
! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
! zmlres        Residual melt water available for infiltration into the
!               non-frozen soil after being intercepted by the
!               skin reservoir
!
!       Rest folgt spaeter ....
!
!     *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
!              CALCULATE FLUXES OF TOTAL RAIN, TOTAL SNOW AND EVAPO-
!              RATION FROM THE THREE RESERVOIRS (SN, WS, WL)
!              CONVERT FLUXES (KG/M**2*S) TO CHANGES OF WATER LEVELS (M)
!              DURING TIMESTEP DELTA_TIME.
!
!     J.F.GELEYN     E.C.M.W.F.     08/06/82.
!     MODIFIED BY
!     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
!     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
!     J.-P. SCHULZ   MPI - 1997 : IMPLEMENTATION OF IMPLICIT
!                                 COUPLING BETWEEN LAND SURFACE
!                                 AND ATMOSPHERE.
!     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
!     MODIFIED BY M. ESCH        MPI - APR  1999
!     MODIFIED BY E. ROECKNER    MPI - JAN  2001
!     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
!     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir 
!                                                for snow changed
!     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
!
!     MODIFICATION
!
!     PURPOSE
!
!     INTERFACE.
!
!          *SURF* IS CALLED FROM *PHYSC*.
!
!     METHOD.
!
!     EXTERNALS.
!
!          NONE.
!
!     REFERENCE.
!
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
!
 
  INTEGER :: jl
!
!  local arrays
!
  REAL(dp) ::                                                         &
       zraind (kbdim),       zsnowd(kbdim),         zevttd(kbdim)     &
     , zevsnd(kbdim),                               zevwld(kbdim)     &
     , zros(kbdim),          zdrain(kbdim),         zrogl(kbdim)      &
     , zsnmel(kbdim),        zsn(kbdim)                               &
     , zmlres(kbdim),        zsncmelt(kbdim)                          &
     , zdp(kbdim),           zlfdcp(kbdim)!

  ! local scalars
  REAL(dp) ::                                                         &
       zorvari, zorvars, zdrmin, zdrmax, zdrexp, zsmelt, zsnmlt,      &
       zmprcp, zwlp, zwdtr, zwslim, zconw2, zconw3,            &
       zroeff, zbws, zb1, zbm, zconw1, zlyeps, zvol, zinfil, zprfl,   &
       zlysic, zwsup, zc2, zc3, zsncp, zexpt, zexpw, zsncmax, zsncwind
  REAL(dp) :: zrcp, zsncfac

  zorvari=100._dp
  zorvars=1000._dp*64._dp/REAL(ngl)
  zdrmin=0.001_dp/(3600._dp*1000._dp)
  zdrmax=0.1_dp/(3600._dp*1000._dp)
  zdrexp=1.5_dp
  zc2=1.87E5_dp
  zc3=1.56E5_dp
  zsncfac=rhoh2o*g/zdtime

!     ------------------------------------------------------------------
!
!*    1.     Convert water fluxes to [m water equivalent * timestep]
!
  DO 110 jl=1,kproma
     zrcp=1._dp/(cpd*(1._dp+vtmpc2*MAX(0.0_dp,pqm1(jl,klev))))
     zlfdcp(jl)=alf*zrcp
     zdp(jl)=paphm1(jl,klev+1)-paphm1(jl,klev)
     zsnmel(jl)=0._dp
     zsncmelt(jl)=0._dp
     zros(jl)=0._dp
     zdrain(jl)=0._dp
     zsn(jl)=0._dp
     zmlres(jl)=0._dp
     palac(jl)=0._dp
     zrogl(jl)=0._dp
     zraind(jl)=(prsfl(jl)+prsfc(jl))             *zdtime/rhoh2o
     zsnowd(jl)=(pssfl(jl)+pssfc(jl))             *zdtime/rhoh2o
     zevttd(jl)=pevapl(jl)                        *zdtime/rhoh2o
     zevsnd(jl)=pcvs(jl)*pevapot(jl)              *zdtime/rhoh2o
     zevwld(jl)=(1._dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o
     zevwsd(jl)=zevttd(jl)-zevsnd(jl)
     ztte(jl,:) = 0._dp
 110 END DO
!
!     ------------------------------------------------------------------
!
!*    2.     Budgets of snow (canopy and ground) and glaciers
!
!*    2.1    Snow changes in the canopy (interception of snowfall,
!            sublimation, melting, unloading due to wind)
!
      DO 210 jl=1,kproma
         if (lpland(jl) .and. (.not. lpglac(jl))) then
            zsn(jl)=zsnowd(jl)+zevsnd(jl)
            zsncmax=MAX(0.0_dp,pwlmx(jl)-cwlmax)
            zmprcp=MIN(zsnowd(jl)*cvinter,zsncmax-psnc(jl))
            zsncp=psnc(jl)+zmprcp
            zsnowd(jl)=zsnowd(jl)-zmprcp
            psnc(jl)=MIN(MAX(0._dp,zsncp+zevsnd(jl)),zsncmax)
            zevsnd(jl)=zevsnd(jl)-(psnc(jl)-zsncp)

            zexpt=MAX(0._dp,ptm1(jl,klev)+3._dp-tmelt)*zdtime/zc2

            zexpw=SQRT(pu10(jl)**2+pv10(jl)**2)*zdtime/zc3
            zsncmelt(jl)=psnc(jl)*(1._dp-EXP(-zexpt))
            psnc(jl)=psnc(jl)-zsncmelt(jl)
            zsncwind=psnc(jl)*(1._dp-EXP(-zexpw))
            psnc(jl)=psnc(jl)-zsncwind
            zsnowd(jl)=zsnowd(jl)+zsncwind
            ztte(jl,klev) = - (zsncmelt(jl)*zsncfac*zlfdcp(jl)/zdp(jl))
!
!   pwl(jl)=pwl(jl)+zsncmelt(jl) see section 2.5
!
         ELSE
            psnc(jl)=0._dp
         END IF
          
 210  END DO

!
!*    2.2    Snowfall and sublimation on land (excluding glaciers)
!
      DO 220 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            psn(jl)=psn(jl)+zsnowd(jl)+zevsnd(jl)
            IF (psn(jl).LT.0._dp) THEN
               zevwsd(jl)=zevwsd(jl)+psn(jl)
               psn(jl)=0._dp
            END IF
         ELSE
            psn(jl)=0._dp
         END IF
 220  END DO
!
!*    2.3    Snowfall and sublimation on glaciers and diagnostics
!
      DO 230 jl=1,kproma
         IF (lpglac(jl)) THEN
            pgld(jl)=pgld(jl)+zsnowd(jl)+zevsnd(jl)
            palac(jl)=zraind(jl)+zsnowd(jl)+zevttd(jl)
            zrogl(jl)=zraind(jl)
         END IF
 230  END DO
!
!*    2.4    Snow and glacier melt
!
   IF (.NOT. lstart) THEN
      DO 240 jl=1,kproma
         IF (lpland(jl).AND.ptsl(jl).GT.tmelt) THEN
            IF (lpglac(jl)) THEN
               zsnmel(jl)=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
               pgld(jl)=pgld(jl)-zsnmel(jl)
               zrogl(jl)=zrogl(jl)+zsnmel(jl)
               ptsl(jl)=tmelt
            ELSE IF (psn(jl).GT.0._dp) THEN
               zsmelt=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
               zsnmel(jl)=MIN(psn(jl),zsmelt)
               ptsl(jl)=ptsl(jl)-zsnmel(jl)*alf*rhoh2o/pgrndcapc(jl)
               psn(jl)=MAX(psn(jl)-zsnmel(jl),0._dp)
             END IF
         END IF
 240  END DO
   END IF
!
!*    2.5    Snow budget and meltwater (glacier-free land only)
!
      DO 250 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            pwl(jl)=pwl(jl)+zsncmelt(jl)
            zsnmlt=zsnmel(jl)+MAX(0._dp,pwl(jl)-pwlmx(jl))
            pwl(jl)=MIN(pwlmx(jl),pwl(jl))
            zmlres(jl)=zsnmlt
            zsnmel(jl)=zsnmel(jl)+zsncmelt(jl)
            zsn(jl)=zsn(jl)-zsnmel(jl)
         END IF
 250  END DO
!
!     ------------------------------------------------------------------
!*    3.     Soil temperatures

      CALL soiltemp (   kproma,     kbdim,     zdtime, lstart  &
                       ,ptsl,       ptsoil,    psn             &
                       ,pgrndc,     pgrndd,    pgrndcapc       &
                       ,pgrndhflx,  psodif,    prgcgn          &
                       ,lpland,     lpglac   )
!     ------------------------------------------------------------------
!
!*    4.     Water budget
!
!*    4.1    Skin reservoir (vegetation and bare soil)
!
      DO 410 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
!
!*    4.1.1  Interception of rain
!
            zmprcp=MIN(zraind(jl)*cvinter,pwlmx(jl)-pwl(jl))
            zwlp=pwl(jl)+zmprcp
            zraind(jl)=zraind(jl)-zmprcp
!
!*    4.1.2  Evaporation or dew collection
!
            pwl(jl)=MIN(MAX(0._dp,zwlp+zevwld(jl)),pwlmx(jl))
            zevwsd(jl)=zevwsd(jl)-(pwl(jl)-zwlp)
          ELSE
           pwl(jl)=0._dp
          END IF
 410  END DO
!
!*    4.2    Soil reservoir
!
      DO 420 jl=1,kproma
          IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
             !zwptr=0.90_dp*pwsmx(jl)
             zwdtr=0.90_dp*pwsmx(jl)
             zwslim=0.05_dp*pwsmx(jl)
             zconw2=pwsmx(jl)-zwdtr
             zconw3=zdrmax-zdrmin
             zroeff=MAX(0._dp, porostd(jl)-zorvari)   &
                          /(porostd(jl)+zorvars)
             zbws=MAX(MIN(zroeff,0.5_dp),0.01_dp)
             zb1=1._dp+zbws
             zbm=1._dp/zb1
             zconw1=pwsmx(jl)*zb1
             zlyeps=0._dp
             zvol=0._dp
             zinfil=0._dp
!
!*    4.2.1  Surface runoff, infiltration and evaporation from soil
!
             IF (zevwsd(jl) >= 0.0_dp) THEN
               zprfl=zmlres(jl)+zraind(jl)+zevwsd(jl)
             ELSE
                pws(jl)=pws(jl)+zevwsd(jl)
               zprfl=zmlres(jl)+zraind(jl)
             END IF
             IF (ptsoil(jl,1).LT.tmelt) THEN
                zros(jl)=zprfl
             ELSE
                IF (zprfl.GT.0._dp.AND.pws(jl).GT.zwslim) THEN
                   IF (pws(jl).GT.pwsmx(jl)) THEN
                      zlyeps=pws(jl)-pwsmx(jl)
                   ELSE
                      zlyeps=0._dp
                   END IF
                   zlysic=(pws(jl)-zlyeps)/pwsmx(jl)
                   zlysic=MIN(zlysic,1._dp)
                   zvol=(1._dp-zlysic)**zbm-zprfl/zconw1
                   zros(jl)=zprfl-(pwsmx(jl)-pws(jl))
                   IF (zvol.GT.0._dp) THEN
                      zros(jl)=zros(jl)+pwsmx(jl)*zvol**zb1
                   END IF
                   zros(jl)=MAX(zros(jl),0._dp)
                   zinfil=zprfl-zros(jl)
                ELSE
                   zros(jl)=0._dp
                   zinfil=zprfl
                END IF
                pws(jl)=pws(jl)+zinfil
             END IF
!
!*    4.2.2  Drainage and total runoff
!
             IF (pws(jl).LE.zwslim) THEN
                zdrain(jl)=0._dp
             ELSE
                IF (ptsoil(jl,1).GT.tmelt) THEN
                   zdrain(jl)=zdrmin*pws(jl)/pwsmx(jl)
                   IF (pws(jl).GT.zwdtr) THEN
                      zdrain(jl)=zdrain(jl)+zconw3*                  &
                                ((pws(jl)-zwdtr)/zconw2)**zdrexp
                   END IF
                   zdrain(jl)=zdrain(jl)*zdtime
                   zdrain(jl)=MIN(zdrain(jl),pws(jl)-zwslim)
                   pws(jl)=pws(jl)-zdrain(jl)
                ELSE
                   zdrain(jl)=0._dp
                END IF
             END IF
             zwsup=MAX(pws(jl)-pwsmx(jl),0._dp)
             pws(jl)=pws(jl)-zwsup
             zros(jl)=zros(jl)+zdrain(jl)+zwsup
          ELSE
             pws(jl)=0._dp
          END IF
 420   END DO
!
!*     4.2.3  Runoff and drainage for the HD-Model
!
      pros_hd(:) = 0.0_dp
      pdrain_hd(:) = 0.0_dp
      DO 423 jl=1,kproma
        pros_hd(jl)=zros(jl)-zdrain(jl)
        pdrain_hd(jl)=zdrain(jl)
 423  END DO
!
!     ------------------------------------------------------------------
!
!*    5.     Time filter for surface temperature
!
   IF (.NOT.lstart) THEN
      DO 510 jl=1,kproma
         ptslm1(jl)=ptslm(jl)+eps*(ptslm1(jl)-2._dp*ptslm(jl)+ptsl(jl))
         ptslm(jl)=ptsl(jl)
  510  END DO
   ELSE
      DO 511 jl=1,kproma
         ptslm1(jl)=ptslm(jl)
 511  END DO
   END IF
!
!     ------------------------------------------------------------------
!
!*    6.     Accumulate fluxes for diagnostics
!
!     6.1    Water fluxes
!
      DO 601 jl=1,kproma
         prunoff(jl)= prunoff(jl) +zros(jl)   *rhoh2o*pslm(jl)
         psnmel(jl) = psnmel(jl)  +zsnmel(jl) *rhoh2o*pslm(jl)
         papmegl(jl)= papmegl(jl) +palac(jl)  *rhoh2o*pslm(jl)
         pdrain(jl) = pdrain(jl)  +zdrain(jl) *rhoh2o*pslm(jl)
         psnacl(jl) = psnacl(jl)  +zsn(jl)    *rhoh2o*pslm(jl)
         progl(jl)  = progl(jl)   +zrogl(jl)  *rhoh2o*pslm(jl)
 601  END DO
!
!     6.2     Ground heat flux
!
      DO 602 jl=1,kproma
         pgrndflux(jl)=pgrndflux(jl)+pgrndhflx(jl)*pslm(jl)*zdtime
 602  END DO
!
  RETURN

END SUBROUTINE surface

!***********************************************************************

!***********************************************************************

  SUBROUTINE lake ( kproma, zdtime                                     &
                  , pseaice,    psiced,    palake                      &
                  , ptsi,       ptsw                                   &
                  , pahflw,     pahfsw,    pfluxres                    &
                  , ptrflw,     psoflw                                 &
                  , pevapi,     psni,      pcvsi                       &
                  , pahfres,    pfri                      )  

!
!  ---------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER,  INTENT(IN)    :: kproma
  REAL(dp), INTENT(IN)    :: zdtime
  REAL(dp), INTENT(INOUT) :: pseaice(:)
  REAL(dp), INTENT(INOUT) :: psiced(:)
  REAL(dp), INTENT(IN)    :: palake(:)
  REAL(dp), INTENT(INOUT) :: ptsi(:)
  REAL(dp), INTENT(INOUT) :: ptsw(:)
  REAL(dp), INTENT(IN)    :: pahflw(:)
  REAL(dp), INTENT(IN)    :: pahfsw(:)
  REAL(dp), INTENT(INOUT) :: pfluxres(:)
  REAL(dp), INTENT(IN)    :: ptrflw(:)
  REAL(dp), INTENT(IN)    :: psoflw(:)
  REAL(dp), INTENT(IN)    :: pevapi(:)
  REAL(dp), INTENT(IN)    :: psni(:)
  REAL(dp), INTENT(IN)    :: pcvsi(:)
  REAL(dp), INTENT(INOUT) :: pahfres(:)
  REAL(dp), INTENT(IN)    :: pfri(:)

!
! local fields
!
  INTEGER :: jl
  REAL(dp):: zalpha, zalphas, zrho_sn, ziscond, zrhoice        &
           , zdice, zrhoilf, zdtrilf, zfreez, zdmix            &
           , zcpwater, zmixcap, zmcapdt, zmcaprilf, zfluxw, zts        &
           , zfres, zconhflx, zsubice, zhi

! Executable statements
!
! 1. Set up constants
!
  zalpha=2.1656_dp
  zalphas=0.31_dp
  zrho_sn=330._dp
  ziscond=zalpha/zalphas*rhoh2o/zrho_sn
  zrhoice=910._dp
  zdice=0.10_dp
  zrhoilf=zrhoice*alf
  zdtrilf=zdtime/zrhoilf
  zfreez=-zdice/zdtrilf
  zdmix=10._dp
  zcpwater=4218._dp
  zmixcap=rhoh2o*zcpwater*zdmix
  zmcapdt=zdtime/zmixcap
  zmcaprilf=zmixcap/zrhoilf
!
! 2. Lake temperature and ice thickness
!
  DO jl=1,kproma
!
  IF (palake(jl).GE.0.5_dp) THEN                           ! lake points
!
     IF (pseaice(jl).LT.0.5_dp) THEN                       ! open water
!
        zfluxw=pahflw(jl)+pahfsw(jl)+ptrflw(jl)+psoflw(jl)
!
!       Lake temperature (ptsw)
!
        zts=ptsw(jl)+zmcapdt*(zfluxw+pfluxres(jl))
        ptsi(jl)=tmelt
        pfluxres(jl)=0._dp
        psiced(jl)=0._dp
        IF (zts.GE.tmelt) THEN                  ! open water (unchanged)
           ptsw(jl)=zts
        ELSE                                    ! check ice formation
           ptsw(jl)=tmelt
           zfres=(zts-tmelt)/zmcapdt            ! < 0.
           IF (zfres.LE.zfreez) THEN            ! ice formation
              psiced(jl)=zmcaprilf*(tmelt-zts)  ! > zdice
              pseaice(jl)=1._dp
           ELSE
              pfluxres(jl)=zfres
           END IF
        END IF
!  ---------------------------------------------------------------------
     ELSE IF (psiced(jl).GE.zdice) THEN
!
!       Ice thickness (psiced)
!
        zconhflx=zalpha*(ptsi(jl)-tmelt)/(psiced(jl)+ziscond*psni(jl))
        zsubice=(1._dp-pcvsi(jl))*pevapi(jl)*zdtime/zrhoice
        zhi=psiced(jl)-zdtrilf*(zconhflx+pfluxres(jl))+zsubice
        ptsw(jl)=tmelt
        IF (zhi.GE.zdice) THEN
           psiced(jl)=zhi
           pseaice(jl)=1._dp
           pfluxres(jl)=0._dp
        ELSE IF (zhi.LE.0._dp) THEN               ! complete melting
           ptsw(jl)=tmelt-zhi/zmcaprilf        ! ptsw > tmelt
           psiced(jl)=0._dp
           pseaice(jl)=0._dp
           pfluxres(jl)=0._dp
        ELSE                                   ! incomplete melting
           psiced(jl)=zdice
           pseaice(jl)=1._dp
           pfluxres(jl)=(zdice-zhi)/zdtrilf
           pahfres(jl)=pahfres(jl)-zdtime*pfri(jl)*pfluxres(jl)
        END IF
     END IF
   END IF
  END DO

     RETURN
  END SUBROUTINE lake

! ************************************************************************

!*************************************************************************
  SUBROUTINE licetemp (kproma, zdtime                                  &
                  , psiced,     psni,      palake                      &
                  , ptsi,       ptrfli,    psofli                      &
                  , pahfice,    pfluxres                               &
                  , pahfcon,    pahfres,   pevapi                      &
                  , pssfl,      pssfc                                  &
                  , pahfsi,     pahfli,     pcvsi                      &
                  , pfri   )

  ! Description:
  !
  ! Prognostic calculation of lake-ice temperature
  !
  ! Method:
  !
  ! *licetemp* called from physc
  ! *physc* called gpc
  !
  ! Authors:
  !
  ! E. Roeckner, MPI, June 2000
  !
  ! for more details see file AUTHORS
  !
!
  IMPLICIT NONE
!
  INTEGER,  INTENT(IN)    :: kproma
  REAL(dp), INTENT(IN)    :: zdtime
  REAL(dp), INTENT(IN)    :: psiced(:)
  REAL(dp), INTENT(INOUT) :: psni(:)
  REAL(dp), INTENT(IN)    :: palake(:)
  REAL(dp), INTENT(INOUT) :: ptsi(:)
  REAL(dp), INTENT(IN)    :: ptrfli(:)
  REAL(dp), INTENT(IN)    :: psofli(:)
  REAL(dp), INTENT(INOUT) :: pahfice(:) 
  REAL(dp), INTENT(INOUT) :: pfluxres(:)  
  REAL(dp), INTENT(INOUT) :: pahfcon(:)
  REAL(dp), INTENT(INOUT) :: pahfres(:)
  REAL(dp), INTENT(IN)    :: pevapi(:)
  REAL(dp), INTENT(IN)    :: pssfl(:)
  REAL(dp), INTENT(IN)    :: pssfc(:)
  REAL(dp), INTENT(IN)    :: pahfsi(:)
  REAL(dp), INTENT(IN)    :: pahfli(:)
  REAL(dp), INTENT(IN)    :: pcvsi(:)
  REAL(dp), INTENT(IN)    :: pfri(:)
!
! Arguments

  INTEGER :: jl
  REAL(dp):: zalpha, zalphas, zrho_sn, ziscond, zcpice, zrhoice        &
           , zdice, zcpcon, zcpdt, zsnowd, zevsnd, zsniced, zicefl     &
           , zsflx, zmelfac, zsmelt

!  Executable statements
!
! 1. Set up constants
!
  zalpha=2.1656_dp
  zalphas=0.31_dp
  zrho_sn=330._dp
  ziscond=zalpha/zalphas*rhoh2o/zrho_sn
  zcpice=2106._dp
  zrhoice=910._dp
  zdice=0.10_dp
  zcpcon=zrhoice*zcpice*zdice
  zcpdt=zcpcon/zdtime

! 2. Compute new skin-temperature

   IF (lice) THEN

      DO jl=1,kproma
      IF (palake(jl).GE.0.5_dp) THEN
         IF (psiced(jl).GE.zdice) THEN                         ! ice
            zsnowd=(pssfl(jl)+pssfc(jl))*zdtime/rhoh2o
            zevsnd=pcvsi(jl)*pevapi(jl)*zdtime/rhoh2o
            psni(jl)=MAX(psni(jl)+zsnowd+zevsnd,0._dp)
            zsniced=psiced(jl)+ziscond*psni(jl)
            zicefl=zalpha*tmelt/zsniced
            zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)          &
                  +pfluxres(jl)
            pfluxres(jl)=0._dp
            ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/                    &
                                         (zcpdt+zalpha/zsniced)
            IF (ptsi(jl).GT.tmelt) THEN
               zmelfac=(zcpdt+zalpha/zsniced)*zdtime/(alf*rhoh2o)
               zsmelt=MIN(zmelfac*(ptsi(jl)-tmelt),psni(jl))
               psni(jl)=psni(jl)-zsmelt
               zsniced=psiced(jl)+ziscond*psni(jl)
               ptsi(jl)=ptsi(jl)-zsmelt/zmelfac
               pahfres(jl)=pahfres(jl)+zsmelt*alf*rhoh2o*pfri(jl)
            END IF
            IF (ptsi(jl).GT.tmelt) THEN
               pfluxres(jl)=(zcpdt+zalpha/zsniced)*(ptsi(jl)-tmelt)
               ptsi(jl)=tmelt
            END IF
            pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pfluxres(jl)
            pahfice(jl)=zalpha*(ptsi(jl)-tmelt)/zsniced
         ELSE                                                 ! water
            pahfice(jl)=0._dp
            ptsi(jl)=tmelt
            psni(jl)=0._dp
         END IF
         pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
      END IF
      END DO

!        Necessary computations if subroutine is bypassed
   ELSE
         DO jl = 1, kproma
         ptsi(jl)=tmelt
         psni(jl)=0._dp
         END DO
   END IF

  RETURN
END SUBROUTINE licetemp

! *******************************************************************

! *******************************************************************

SUBROUTINE soiltemp (  kproma,     kbdim,     zdtime, lstart        &
                      ,pts,        ptsoil,    psn                   &
                      ,pgrndc,     pgrndd,    pgrndcapc             &
                      ,pgrndhflx,  psodif,    prgcgn                &
                      ,ldland,     ldglac   )
!
!   AUTHOR:  FREDERIC HOURDIN     30/01/92
!
!            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
!            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
!
!            J.-P. SCHULZ   MPI - OCTOBER 1997 :
!               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
!               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
!               ECHAM4 GCM.
!            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
!            U.Schlese DKRZ - February 2000  new soil temperatures
!            L Kornblueh, MPI, January 2003, removed MERGE
!
!
!   OBJECTIVE:  COMPUTATION OF:
!               THE GROUND TEMPERATURE EVOLUTION
!               THE GROUND SPECIFIC HEAT "CAPCAL"
!               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
!
!
!   METHOD:  IMPLICIT TIME INTEGRATION
!
!   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
!           T(K+1) = C(K) + D(K)*T(K)  (1)
!   THE COEFFICIENTS C (=GRNDC) AND D (=GRNDD) ARE COMPUTED AT THE
!   T-DT TIME-STEP.
!   ROUTINE STRUCTURE:
!   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
!   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
!     PROFILE FOR THE T+DT TIME-STEP
!   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
!     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
!            FDIFF = A + B TS(T+DT)
!     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
!            WITH F0 = A + B (TS(T))
!                 CAPCAL = B*DT
!
!   INTERFACE:
!
!   ARGUMENTS:
!
!   INPUT:
!
!   PTIMESTEP         TIME-STEP (S)
!   PTSOL             INITIAL TEMPERATURE AT SOIL SURFACE
!   ZSO_TDIF          SOIL TEMPERATURE DIFFUSIVITY [M**2/S] (FAO)
!   ZSO_CAPA          SOIL VOL. HEAT CAPACITY    [J/M**3/K] (FAO)
!   SNOW              SNOW DEPTH (MM LIQUID WATER EQUIVALENT)
!
!   OUTPUT:
!
!   PTN(GRID_SIZE,NGRNDMX)   GROUND TEMPERATURES OF NGRNDMX LAYERS
!   CGRND(GRID_SIZE,NGRNDMX) COEFFICIENT OF THE SOIL TEMPERATURE SCHEME
!   DGRND(GRID_SIZE,NGRNDMX) COEFFICIENT OF THE SOIL TEMPERATURE SCHEME
!
!     ------------------------------------------------------------------
!
!   DECLARATIOION                      
!
  IMPLICIT NONE

  INTEGER,  INTENT(IN)    :: kproma
  INTEGER,  INTENT(IN)    :: kbdim
  REAL(dp), INTENT(IN)    :: zdtime
  LOGICAL,  INTENT(IN)    :: lstart
  REAL(dp), INTENT(IN)    :: pts(:)
  REAL(dp), INTENT(INOUT) :: ptsoil(:,:)
  REAL(dp), INTENT(IN)    :: psn(:)
  REAL(dp), INTENT(INOUT) :: pgrndc(:,:)
  REAL(dp), INTENT(INOUT) :: pgrndd(:,:)
  REAL(dp), INTENT(INOUT) :: pgrndcapc(:)
  REAL(dp), INTENT(INOUT) :: pgrndhflx(:)
  REAL(dp), INTENT(IN)    :: psodif(:)
  REAL(dp), INTENT(IN)    :: prgcgn(:)
  LOGICAL,  INTENT(IN)    :: ldland(:)
  LOGICAL,  INTENT(IN)    :: ldglac(:) 
!
!  local Variables
!
  INTEGER :: jl, jk
  REAL(dp):: zso_cond(kbdim), zso_capa(kbdim)
  REAL(dp):: z1(kbdim)
  REAL(dp):: zd1(jpgrnd)
  REAL(dp):: zdz1(kbdim,jpgrnd),   zdz2(kbdim,jpgrnd)
  REAL(dp):: zkappa(kbdim,jpgrnd), zcapa(kbdim,jpgrnd)
  REAL(dp):: zsnow_h, zx1, zx2
  REAL(dp):: zrici, zdifiz, zsn_cond, zsn_dens, zsn_capa

  REAL(dp), DIMENSION(:), POINTER :: cdel, cmid    ! from iniphy
!
!     ------------------------------------------------------------------
!
!*    1.  SPECIFYING THE DEPTHS OF THE TEMPERATURE LEVELS.
!
!*    1.1 SOME CONSTANTS USED IN THE TEMPERATURE SCHEME.
!
  zrici = 2.09e+06_dp        ! volumetric heat capacity of ice [j/m**3/k]
  zdifiz = 12.e-07_dp        ! temperature diffusivity of ice  [m**2/s]
  zsn_cond = 0.31_dp         ! snow thermal conductivity [j/s/m/k]
  zsn_dens = 330.0_dp        ! snow density              [kg/m**3]
  zsn_capa = 634500.0_dp     ! snow vol. heat capacity   [j/m**3/k]

!
!*    1.2 COMPUTING SOME USEFUL CONSTANTS.
!
! declaration imported from iniphy, only used here
!
   allocate(cdel(jpgrnd))
   allocate(cmid(jpgrnd))

!  THICKNESS OF SOIL LAYERS
!
      CDEL(1)=0.065_dp
      CDEL(2)=0.254_dp
      CDEL(3)=0.913_dp
      CDEL(4)=2.902_dp
      CDEL(5)=5.700_dp
!
!  DEPTH OF MIDS OF SOIL LAYERS
!
      CMID(1)=CDEL(1)*0.5_dp
      CMID(2)=CDEL(1)+CDEL(2)*0.5_dp
      CMID(3)=CDEL(1)+CDEL(2)+CDEL(3)*0.5_dp
      CMID(4)=CDEL(1)+CDEL(2)+CDEL(3)+CDEL(4)*0.5_dp
      CMID(5)=CDEL(1)+CDEL(2)+CDEL(3)+CDEL(4)+CDEL(5)*0.5_dp

  DO jk = 1,jpgrnd-1
     zd1(jk) = 1._dp/(cmid(jk+1)-cmid(jk))
  END DO

!
!*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
!*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
!
  DO jl = 1,kproma
    IF (ldglac(jl)) THEN
      zso_capa(jl) = zrici
      zso_cond(jl) = zso_capa(jl)*zdifiz
    ELSE
      zso_capa(jl) = prgcgn(jl)
      zso_cond(jl) = zso_capa(jl)*psodif(jl)
    END IF
  END DO
!
!*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
!
  DO jk = 1,jpgrnd
     DO jl = 1,kproma
        zkappa(jl,jk) = zso_cond(jl)
        zcapa(jl,jk)  = zso_capa(jl)
     END DO
  END DO
!
!   --------------------------------------------------------------
!   COMPUTATION OF THE GROUND TEMPERATURES USING THE CGRD AND DGRD
!   COEFFICIENTS COMPUTED AT THE PREVIOUS TIME-STEP
!   --------------------------------------------------------------
!
  IF(.NOT.lstart) THEN

! NOTE: At lstart, tsoil (g3b) is initialised in initemp (called from scan1).
!       This should be moved to surface_global_start or surface_initialize
!       and removed from initemp!

!
!   Upper layer
!
    DO jl = 1,kproma
       IF (ldland(jl)) THEN
             ptsoil(jl,1)=pts(jl)
       END IF
    END DO
!
!   Deeper layers
!
    DO jk = 1,jpgrnd-1
       DO jl = 1,kproma
          IF (ldland(jl)) THEN
           ptsoil(jl,jk+1)=pgrndc(jl,jk)+pgrndd(jl,jk)*ptsoil(jl,jk)
          END IF
       END DO
    END DO
  END IF
!
!   ---------------------------------------------------------------
!   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
!   ---------------------------------------------------------------
!
  DO jl = 1,kproma
     IF (ldland(jl)) THEN
           zsnow_h = psn(jl)*rhoh2o / zsn_dens
!
!*       Special treatment for first layer
!
        IF ( zsnow_h .GT. cmid(2) ) THEN
           zcapa(jl,1) = zsn_capa
           zkappa(jl,1) = zsn_cond
        ELSE IF ( zsnow_h .GT. 0.0_dp .AND. zsnow_h .LE. cmid(2) ) THEN
           zx1 = zsnow_h / cmid(2)
           zx2 = ( cmid(2) - zsnow_h) / cmid(2)
           zcapa(jl,1) = zx1 * zsn_capa + zx2 * zso_capa(jl)
           zkappa(jl,1) = 1.0_dp / ( zx1 / zsn_cond +                  &
                                                  zx2 / zso_cond(jl) )
        ELSE
           zcapa(jl,1) = zso_capa(jl)
           zkappa(jl,1) = zso_cond(jl)
        ENDIF
!
        DO jk = 2, jpgrnd - 2
           IF ( zsnow_h .GT. cmid(jk+1) ) THEN
              zcapa(jl,jk) = zsn_capa
              zkappa(jl,jk) = zsn_cond
           ELSE IF ( zsnow_h .GT. cmid(jk) .AND.                       &
                                    zsnow_h .LE. cmid(jk+1) ) THEN
              zx1 = (zsnow_h - cmid(jk)) * zd1(jk)
              zx2 = ( cmid(jk+1) - zsnow_h) * zd1(jk)
              zcapa(jl,jk) = zx1*zsn_capa + zx2*zso_capa(jl)
              zkappa(jl,jk) = 1.0_dp / ( zx1 / zsn_cond +              &
                                                  zx2 / zso_cond(jl) )
           ELSE
              zcapa(jl,jk) = zso_capa(jl)
              zkappa(jl,jk) = zso_cond(jl)
           END IF
        END DO
     END IF
  END DO
!
  DO jk=1,jpgrnd
     DO jl=1,kproma
        IF (ldland(jl)) THEN
           zdz2(jl,jk)=zcapa(jl,jk)*cdel(jk)/zdtime
        END IF
     END DO
  END DO
!
  DO jk=1,jpgrnd-1
     DO jl=1,kproma
        IF (ldland(jl)) THEN
           zdz1(jl,jk)=zd1(jk)*zkappa(jl,jk)
        END IF
     END DO
  END DO
!
  DO jl=1,kproma
     IF (ldland(jl)) THEN
        z1(jl)=zdz2(jl,jpgrnd)+zdz1(jl,jpgrnd-1)
        pgrndc(jl,jpgrnd-1)=zdz2(jl,jpgrnd)*ptsoil(jl,jpgrnd)/z1(jl)
        pgrndd(jl,jpgrnd-1)=zdz1(jl,jpgrnd-1)/z1(jl)
     END IF
  END DO
!
  DO jk=jpgrnd-1,2,-1
     DO jl=1,kproma
        IF (ldland(jl)) THEN
           z1(jl)=1._dp/(zdz2(jl,jk)+zdz1(jl,jk-1) +                   &
                                 zdz1(jl,jk)*(1._dp-pgrndd(jl,jk)))
           pgrndc(jl,jk-1)=(ptsoil(jl,jk)*zdz2(jl,jk) +                &
                                 zdz1(jl,jk)*pgrndc(jl,jk))*z1(jl)
           pgrndd(jl,jk-1)=zdz1(jl,jk-1)*z1(jl)
        END IF
     END DO
  END DO
!
!   ---------------------------------------------------------
!   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
!   CALORIFIC CAPACITY OF THE GROUND:
!   ---------------------------------------------------------
!
  DO jl=1,kproma
     IF (ldland(jl)) THEN
        pgrndhflx(jl)=zdz1(jl,1)*(pgrndc(jl,1)                         &
                                   +(pgrndd(jl,1)-1._dp)*ptsoil(jl,1))
        pgrndcapc(jl)=(zdz2(jl,1)*zdtime+                              &
                           zdtime * (1._dp-pgrndd(jl,1)) * zdz1(jl,1))
     END IF
  END DO
!
!     ------------------------------------------------------------------

  ! release space
  deallocate(cmid)
  deallocate(cdel)
!
  RETURN
END SUBROUTINE soiltemp

! *******************************************************************

! *******************************************************************

SUBROUTINE sicetemp (kproma, zdtime, lcouple, lmlo                  &
                   , psiced,     psni,      palake                  &
                   , pslf                                           &
                   , ptsi,       ptrfli,    psofli                  &
                   , pahfice,    pfluxres,  pqres                   &
                   , pahfcon,    pahfres                            &
                   , pahfsi,     pahfli                             &
                   , pfri                           )

  ! Description:
  !
  ! Prognostic calculation of sea-ice temperature
  !
  ! Method:
  !
  ! *sicetemp* called from physc
  ! *physc*    called from gpc
  !
  ! Authors:
  !
  ! F. Lunkeit, MI, April 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A, Rhodin, MPI, Jan 1999, argument list added
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! R. Voss, U.Schlese, December 1999, modifications for coupling
  ! I. Kirchner, MPI, December 2000, time control
  ! U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
  ! U. Schlese, MPI December 2002, ice thickness over ocean
  !
  ! for more details see file AUTHORS
  !

  USE messy_main_constants_mem,   ONLY: ctfreez  

  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: kproma
  REAL(dp), INTENT(IN)    :: zdtime
  LOGICAL,  INTENT(IN)    :: lcouple
  LOGICAL,  INTENT(IN)    :: lmlo
  REAL(dp), INTENT(IN)    :: psiced(:)
  REAL(dp), INTENT(INOUT) :: psni(:)
  REAL(dp), INTENT(IN)    :: palake(:)
  REAL(dp), INTENT(IN)    :: pslf(:)
  REAL(dp), INTENT(INOUT) :: ptsi(:)
  REAL(dp), INTENT(IN)    :: ptrfli(:)
  REAL(dp), INTENT(IN)    :: psofli(:)
  REAL(dp), INTENT(INOUT) :: pahfice(:)    
  REAL(dp), INTENT(INOUT) :: pfluxres(:)
  REAL(dp), INTENT(INOUT) :: pqres(:)
  REAL(dp), INTENT(INOUT) :: pahfcon(:)

  REAL(dp), INTENT(INOUT) :: pahfres(:)
  REAL(dp), INTENT(IN)    :: pahfsi(:)
  REAL(dp), INTENT(IN)    :: pahfli(:)
  REAL(dp), INTENT(IN)    :: pfri(:)

  REAL(dp):: zalpha, zalphas, zrho_sea, zrho_sn, ziscond, zcpice       &
           , zrhoice, zdice, zcpcon, zcpdt, zsniced, zicefl, zsflx
  INTEGER :: jl
  
!  Executable statements
!
!-- 1. Set up constants

  zalpha = 2.1656_dp
  zalphas=0.31_dp
  zrho_sea=1025._dp
  zrho_sn=330._dp
  ziscond=zalpha/zalphas*zrho_sea/zrho_sn
  zcpice = 2106._dp
  zrhoice = 910._dp
  zdice = 0.10_dp
  zcpcon = zrhoice*zcpice*zdice
  zcpdt = zcpcon/zdtime
!
!-- 2. Compute new skin-temperature
!    
   IF (lice) THEN
!
    IF (.NOT.lmlo .AND. .NOT.lcouple) THEN     ! ECHAM5

      DO jl=1,kproma
        IF (palake(jl).EQ.0._dp) THEN
          IF (psiced(jl).GE.zdice) THEN
            zsniced=psiced(jl)+ziscond*psni(jl)
            zicefl=zalpha*ctfreez/zsniced
            zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
            ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/                    &
                                        (zcpdt+zalpha/zsniced)
            IF (ptsi(jl).GT.tmelt) THEN
               pqres(jl)=(zalpha/zsniced+zcpdt)*(ptsi(jl)-tmelt)
               ptsi(jl)=tmelt
            ELSE
               pqres(jl)=0._dp
            END IF
            pahfice(jl)=zalpha*(ptsi(jl)-ctfreez)/zsniced
          ELSE
            pqres(jl)=0._dp
            pahfice(jl)=0._dp
            ptsi(jl)=tmelt
          END IF
          pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
          pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pqres(jl)
       END IF
     END DO
    ELSE IF (lcouple) THEN                      ! ECHAM5-IPCC
      DO jl=1,kproma
         IF (pslf(jl) .LT. 1.0_dp) THEN
            zsniced=MAX(psiced(jl),zdice)+ziscond*psni(jl)
            zicefl=zalpha*ctfreez/zsniced
            zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
            ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/                    &
                                        (zcpdt+zalpha/zsniced)
            IF (ptsi(jl).GT.tmelt) THEN
               pqres(jl)=(zalpha/zsniced+zcpdt)*(ptsi(jl)-tmelt)
               ptsi(jl)=tmelt
            ELSE
               pqres(jl)=0._dp
            END IF
            pahfice(jl)=zalpha*(ptsi(jl)-ctfreez)/zsniced
            pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
            pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pqres(jl)
        ELSE
            pqres(jl)=0._dp
            pahfice(jl)=0._dp
            ptsi(jl)=tmelt
        END IF
     END DO

   END IF
!
!       Necessary computations if subroutine is bypassed
!
      ELSE
         DO jl = 1, kproma
         ptsi(jl)=tmelt
         psni(jl)=0._dp
         END DO
      END IF
!
  RETURN
END SUBROUTINE sicetemp

! ***********************************************************************

END MODULE messy_surface
