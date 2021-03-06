MODULE messy_main_tools_wiso

  USE messy_main_constants_mem, ONLY: dp, M_H2O, M_HH18O, M_HDO, rho_H2O, tmelt
  
  IMPLICIT NONE
  PUBLIC
  SAVE

  PRIVATE :: dp
  
  ! ----------------------------------
  ! determine the specific phase
  ! ----------------------------------
  INTEGER, PARAMETER :: i_vap = 1
  INTEGER, PARAMETER :: i_liq = 2
  INTEGER, PARAMETER :: i_ice = 3
  INTEGER, PARAMETER :: kphase = 3  ! max. number of phases
  ! strings for loop calling tracers
  CHARACTER(LEN=3), DIMENSION(kphase), PARAMETER :: &
       phasetag = (/'vap', 'liq', 'ice'/)

  ! -----------------------------------
  ! determine the specific isotopologue
  ! -----------------------------------
  INTEGER, PARAMETER :: i_HHO   = 1
  INTEGER, PARAMETER :: i_HH18O = 2
  INTEGER, PARAMETER :: i_HDO   = 3
  INTEGER, PARAMETER :: mwiso   = 3  ! max. number of isotopologues
  ! strings for loops calling tracers
  CHARACTER(LEN=5), DIMENSION(mwiso), PARAMETER :: &
       isotag  = (/'HHO  ', 'HH18O', 'HDO  '/)
  ! molar masses of tracers in the order: HHO, HH18O, HDO
  REAL(DP), DIMENSION(3)   :: molarmasses = (/ M_H2O, M_HH18O, M_HDO /)
  
  ! natural isotope-ratio
  REAL(dp), PARAMETER :: tnat(mwiso) = (/ 1.0_dp, 2005.2e-6_dp, 155.76e-6_dp /)
  ! density of isotopic water
  REAL(dp), PARAMETER :: twisorhoh2o(mwiso) = (/ rho_H2O, rho_H2O, rho_H2O /)
  ! mean ocean concentration of the tracer
  REAL(dp), PARAMETER :: toce(mwiso) = (/ 1.0_dp, tnat(2), tnat(3) /)

  ! initial isotope deviation from SMOW in the atmosphere
  REAL(dp), PARAMETER :: twisoatm(mwiso) = &
       (/ 0.0_dp, -20.0_dp/1000.0_dp, -150.0_dp/1000.0_dp/)
  ! relation of the diffusivities rel. water
  REAL(dp), PARAMETER :: tdifrel(mwiso) = &
       (/ 1.0_dp, 1.0_dp/0.9723_dp, 1.0_dp/0.9755_dp/)
  ! factors for initial deviation from SMOW in the surf. layers
  REAL(dp), PARAMETER :: twisosur1(mwiso) = &
       (/ 0.0_dp, 0.69_dp/1000.0_dp, 5.6_dp/1000.0_dp /)
  ! factors for initial deviation from SMOW in the surf. layers
  REAL(dp), PARAMETER :: twisosur2(mwiso) = &
       (/ 0.0_dp, 13.6_dp/1000.0_dp, 100.0_dp/1000.0_dp/)
  ! factors for the kinetic fractionation over ocean
  REAL(dp), PARAMETER :: tkinsl(mwiso) = (/ 0.0_dp, 0.006_dp, 0.00528_dp /)
  ! factors for the kinetic fractionation over ocean
  REAL(dp), PARAMETER :: tkinfa1(mwiso) = &
       (/ 0.0_dp, 0.000285_dp, 0.0002508_dp/)
  ! factors for the kinetic fractionation over ocean
  REAL(dp), PARAMETER :: tkinfa2(mwiso) = &
       (/ 0.0_dp, 0.00082_dp, 0.0007216_dp/)
  ! factors for the computation of the fractionation over liquid
  REAL(dp), PARAMETER :: talphal1(mwiso) = &
       (/ 0._dp, 1137._dp, 24844._dp /)
  ! factors for the computation of the fractionation over liquid
  REAL(dp), PARAMETER :: talphal2(mwiso) = &
       (/ 0._dp, -0.4156_dp, -76.248_dp /)
  ! factors for the computation of the fractionation over liquid
  REAL(dp), PARAMETER :: talphal3(mwiso) = &
       (/ 0._dp, -2.0667e-3_dp, 52.612e-3_dp /)

  ! factors for the computation of the fractionation over ice
  REAL(dp), PARAMETER :: talphas1(mwiso) = (/ 0._dp, 0._dp, 16288._dp /)
  ! by M. Ellehoj:
  !REAL(dp), PARAMETER :: talphas1(mwiso) = (/ 0._dp,  8312.5_dp, 48888._dp /)

  ! factors for the computation of the fractionation over ice 
  REAL(dp), PARAMETER :: talphas2(mwiso) = (/ 0._dp, 11.839_dp, 0._dp /)
  ! by M. Ellehoj:
  !REAL(dp), PARAMETER :: talphas2(mwiso) = (/ 0._dp,  -49.192_dp, -203.10_dp /)

  ! factors for the computation of the fractionation over ice 
  REAL(dp), PARAMETER :: talphas3(mwiso) = (/ 0._dp, -0.028244_dp, -0.0934_dp /)
  ! by M. Ellehoj:
  !REAL(dp), PARAMETER :: talphas3(mwiso) = (/ 0._dp,  0.0831_dp, 0.2133_dp /)
  
  ! factors for the computation of the fractionation
  REAL(dp), PARAMETER :: talph1(mwiso) = (/ 0.0_dp, 1137.0_dp, 24844.0_dp /)
  ! factors for the computation of the fractionation
  REAL(dp), PARAMETER :: talph2(mwiso) = (/ 0.0_dp, -0.4156_dp, -76.248_dp /)
  ! factors for the computation of the fractionation
  REAL(dp), PARAMETER :: talph3(mwiso) = &
       (/ 0.0_dp, -2.0667e-3_dp, 52.612e-3_dp /)
  
  ! factors for the computation of the fractionation over ice
  REAL(dp), PARAMETER :: talps1(mwiso) = (/ 0.0_dp, 11.839_dp, 16288.0_dp /)
  ! factors for the computation of the fractionation over ice 
  REAL(dp), PARAMETER :: talps2(mwiso) = (/ 0.0_dp, -0.028244_dp, -0.0934_dp /)

  ! base value for temperatur-dependency of supersaturation
  ! during kinetic effects
  REAL(dp), PARAMETER :: tsatbase=1.01_dp
  !REAL(dp), PARAMETER :: tsatbase=1.0_dp   ! M. Ellehoj
  ! temperatur-dependency of supersaturation during kinetic effects
  REAL(dp), PARAMETER :: tsatfac=0.0045_dp 
  !REAL(dp), PARAMETER :: tsatfac=0.002_dp  ! M. Ellehoj

  ! a constant for alpha_eff for equilibrium below cloud base
  REAL(dp), PARAMETER :: tdifexp=0.58_dp
  ! coefficient for effective fractionation below cumulus clouds
  REAL(dp), PARAMETER :: thumwiso1=0.75_dp
  ! coefficient for effective fractionation below cumulus clouds
  REAL(dp), PARAMETER :: thumwiso2=0.25_dp
  ! coefficient for raindrop equilibrium below cumulus clouds
  REAL(dp), PARAMETER :: twisoeqcu=0.5_dp
  ! coefficient for raindrop equilibrium below large scale clouds
  REAL(dp), PARAMETER :: twisoeqls=0.9_dp
  
  ! temperature limit for homogenous ice cumulus clouds
  ! (set to same default Echam5 value (-35C) as for large-scale clouds)
  !!$REAL(dp), PARAMETER :: twisoice=cthomi
  REAL(dp), PARAMETER :: cthomi  = tmelt - 35.0_dp
  
  ! minimum threshold value for the snow layer depth on glaciers 
  REAL(dp), PARAMETER :: snglacmn=0.001_dp
  ! maximum threshold value for the snow layer depth on glaciers 
  REAL(dp), PARAMETER :: snglacmx=0.02_dp

  INTEGER,  PARAMETER :: lwisofracl=0 ! switch for evapotranspiration over land:
                                      ! 0: no fractionation over land surface
                                      ! 1: fractionation for soil evaporation,
                                      !    only
                                      ! 2: fractionation for soil evaporation
                                      !    and plant transpiration

  INTEGER,  PARAMETER :: lwisokinl=0  ! switch for kinetic fractionation over
                                      ! land:
                                      ! 0: no kinetic fractionation over land
                                      !    surface
                                      ! 1: kinetic fractionation according
                                      !    to Brutsaert (1975) as for
                                      !    frac. over water
                                      ! 2: kinetic fractionation according
                                      !    to Mathieu & Bariac (1996)

  ! mininum value for calculating different water isotopic ratios
  REAL(dp), PARAMETER :: cwisomin= 1.0e-15_dp
  ! security parameter for calculation of delta vaues
  REAL(dp), PARAMETER :: cwisosec= 1.0e-12_dp
  
CONTAINS

  ! ----------------------------------------------------------------------
  SUBROUTINE wiso_frac_liq(kproma, kbdim, kwiso, pt, pwisofracliq)

    ! wiso_frac_liq calculates fractionation coefficients for vapour-liquid
    ! for a band of longitudinal grid points, simultaneously 

    IMPLICIT NONE

    ! input arguments
    INTEGER,  INTENT(IN) :: kproma, kbdim, kwiso
    REAL(dp), INTENT(IN) :: pt(kbdim) ! temperature
    ! input/output arguments
    REAL(dp), INTENT(INOUT) :: pwisofracliq(kbdim,kwiso)
  
    ! local variables
    INTEGER     :: jl, jt
    
    ! fractionation over liquid water
    DO jt=1,kwiso
       DO jl=1,kproma
          pwisofracliq(jl,jt) = exp( talphal1(jt)/(pt(jl)**2) &
               + talphal2(jt)/pt(jl)+talphal3(jt) )
       END DO
    END DO

  END SUBROUTINE wiso_frac_liq
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE wiso_frac_liq_ice(kproma,kbdim,kwiso,pt,pwisofracliq,pwisofracice)

    ! wiso_frac_ice calculates fractionation coefficients for
    ! vapour-liquid and vapour-ice
    ! for a band of longitudinal grid points, simultaneously 

    USE messy_main_constants_mem, ONLY: tmelt
  
    IMPLICIT NONE

    ! input arguments
    INTEGER, INTENT(IN)     :: kproma, kbdim, kwiso
    REAL(dp), INTENT(IN)    :: pt(kbdim)
    ! input/output arguments
    REAL(dp), INTENT(INOUT) :: pwisofracliq(kbdim,kwiso) &
                             , pwisofracice(kbdim,kwiso)
  
    ! local variables
    INTEGER     :: jl,jt
    REAL(dp)    :: zsatval
    
    ! fractionation over liquid water
    DO jt=1,kwiso
       DO jl=1,kproma
          pwisofracliq(jl,jt) = &
               exp(talphal1(jt)/(pt(jl)**2)+talphal2(jt)/pt(jl)+talphal3(jt))
       END DO
    END DO
  
    ! fractionation over ice
    DO jt=1,kwiso
       DO jl=1,kproma
          pwisofracice(jl,jt) = &
               exp(talphas1(jt)/(pt(jl)**2)+talphas2(jt)/pt(jl)+talphas3(jt))
!!$          IF (nwisotyp(jt).ne.1.and.pt(jl).lt.tmelt) THEN
          IF ( (jt /= i_HHO) .AND. (pt(jl) < tmelt)) THEN
             ! effective fractionation over ice if necessary
             zsatval=tsatbase-tsatfac*(pt(jl)-tmelt)
             pwisofracice(jl,jt) = &
                  pwisofracice(jl,jt) * &
                  (zsatval/(1._dp+pwisofracice(jl,jt) * &
                  (zsatval-1._dp)*tdifrel(jt)))
          ENDIF
       END DO
    END DO

  END SUBROUTINE wiso_frac_liq_ice
  ! ----------------------------------------------------------------------
  
  ! ----------------------------------------------------------------------
  SUBROUTINE wiso_cloudadj(kproma,kbdim,klev,kwiso,kk,           &
                        ptmst,pcons2,pdp,                        &
                        ptpone,pqpone,pqspone,pxlpone,           &
                        paclc,                                   &
                        pwisoqm1,  pwisoqte,                     &
                        pwisoxlm1, pwisoxlte,                    &
                        prfl, pwisorfl)   
!
!      G.HOFFMANN          MPI MET, HAMBURG      1992 
!      M. WERNER,          MPI BGC, JENA,        2004
!      M. WERNER           AWI, BREMERHAVEN      2009
!
!      PURPOSE
!      -------
!      CALCULATES THE FRACTIONATION FACTOR (LIQUID/VAPOUR),
!      GETS CLOUD LIQUID ISOTOPE AND VAPOUR AND THEN
!      LIQUID PRECIPITATION AND VAPOUR IN EQUILIBRIUM
!
!      INTERFACE
!      ---------
!      THIS SUBROUTINE IS CALLED FROM
!        *CLOUD*
!
!      INPUT  TEMPERATURE , 
!             CLOUD LIQUID WATER,
!             VAPOUR AND SATURATION MIXING RATIO OF VAPOUR,
!             PRECIPITATION AND FRACTIONAL CLOUD COVER,
!             OLD ISOTOPE MIXING RATIO AND TENDENCY,
!             ISOTOPE PRECIPITATION
!      OUTPUT NEW ISOTOPE TENDENCIES 
!
!      EXTERNALS
!      ---------
!      NONE

    IMPLICIT NONE

    ! input arguments  
    INTEGER, INTENT(IN) :: kproma,kbdim,klev,kwiso,kk
    REAL(dp), INTENT(IN) :: ptmst,pcons2
    REAL(dp), INTENT(IN) :: pdp(kbdim),ptpone(kbdim),pqpone(kbdim),  &
         pqspone(kbdim),pxlpone(kbdim),           &
         pwisoqm1(kbdim,klev,kwiso),              &
         pwisoxlm1(kbdim,klev,kwiso),             &
         prfl(kbdim),paclc(kbdim,klev)
    ! input/output arguments  
    REAL(dp), INTENT(INOUT) :: pwisoqte(kbdim,klev,kwiso),   &
         pwisoxlte(kbdim,klev,kwiso),  &
         pwisorfl(kbdim,kwiso)
    ! local variables
    REAL(dp) :: zwisofracliq(kbdim,kwiso)
    REAL(dp) :: zdelta, zwisoqpone,zwisoxlpone,                           &
         zrafac,zqrain,zwisoqrain,                                 &
         zqinc,zqoutc,zwisoqinc,zwisoqoutc,zdeltab,zdeltamb,zdenom
    INTEGER :: jl,jt

    CALL wiso_frac_liq(kproma,kbdim,kwiso,ptpone,zwisofracliq)

    DO jt=1,kwiso        
       DO jl=1,kproma
          ! Equilibration Part I: bring liquid cloud water in
          ! isotopic equilibrium with surrounding vapour
          zwisoqpone =pwisoqm1(jl,kk,jt) +pwisoqte(jl,kk,jt) *ptmst 
          zwisoxlpone=pwisoxlm1(jl,kk,jt)+pwisoxlte(jl,kk,jt)*ptmst
          IF (pqpone(jl).ge.0.0_dp .and. pxlpone(jl).ge.0.0_dp .and. zwisoqpone.ge.0.0_dp .and. zwisoxlpone.ge.0.0_dp) THEN

             IF (pxlpone(jl).ne.0.0_dp) THEN
                ! if some liquid cloud water will still exist in next time step: 
                ! assume fractionation between liquid and vapour (closed system)
                IF (pqpone(jl).gt.cwisomin) THEN
                   ! apply fractionation for positive vapour values, only
                   IF ((pqpone(jl)+zwisofracliq(jl,jt)*pxlpone(jl)).lt.cwisomin) THEN
                      zdelta=0.0_dp
                   ELSE
                      zdelta=(zwisofracliq(jl,jt)*pxlpone(jl)*zwisoqpone-zwisoxlpone*pqpone(jl))         &
                           /(pqpone(jl)+zwisofracliq(jl,jt)*pxlpone(jl))
                      IF (jt == i_HHO) zdelta = 0.0_dp
                   ENDIF

                   ! limit changes to absolut amount of available vapour
                   ! and liquid water
                   IF (zwisoqpone.GE.0.0_dp) THEN
                      zdelta=MIN(zdelta,zwisoqpone)
                      zdelta=MAX(zdelta,-zwisoxlpone)
                   ELSE
                      zdelta=MAX(zdelta,zwisoqpone)
                      zdelta=MIN(zdelta,-zwisoxlpone)
                   ENDIF
                   pwisoqte(jl,kk,jt) =pwisoqte(jl,kk,jt) -zdelta/ptmst
                   pwisoxlte(jl,kk,jt)=pwisoxlte(jl,kk,jt)+zdelta/ptmst
                ENDIF
             ELSE
                ! if no default cloud water will exist in next time step:
                ! adjust isotope tendencies, accordingly
                zdelta=pwisoxlm1(jl,kk,jt)+pwisoxlte(jl,kk,jt)*ptmst
                pwisoqte(jl,kk,jt) =pwisoqte(jl,kk,jt) +zdelta/ptmst
                pwisoxlte(jl,kk,jt)=pwisoxlte(jl,kk,jt)-zdelta/ptmst
             ENDIF
          ENDIF

          ! Equilibration Part II: bring precipitation in isotopic
          ! equilibrium with surrounding vapour,
          ! distinguish here between rain still within a cloud
          ! versus rain already outside a cloud

          ! convert rain to same units as water vapour
          zrafac=pcons2*pdp(jl)
          zqrain=prfl(jl)/zrafac
          zwisoqrain=pwisorfl(jl,jt)/zrafac
          zwisoqpone=pwisoqm1(jl,kk,jt)+pwisoqte(jl,kk,jt)*ptmst

          IF (pqpone(jl).ge.0.0_dp .and. zqrain.ge.0.0_dp .and. zwisoqpone.ge.0.0_dp .and. zqrain.ge.0.0_dp) THEN

             ! if some rain will exist in next time step:
             ! assume fractionation between rain and vapour (closed system)

             IF (zqrain.ne.0.0_dp) THEN

                IF (pqpone(jl).gt.cwisomin) THEN
                   ! apply fractionation for positive vapour values, only
                   ! if a cloud exist: assume that vapour in cloud
                   ! (zqinc,zwisoinc) is at saturation value with
                   ! a delta value equal to the value of the whole grid cell
                   IF (paclc(jl,kk).gt.1.0e-10_dp) THEN
                      zqinc=pqspone(jl)
                      zdelta=tnat(jt)
                      IF (pqpone(jl).GT.cwisomin) &
                           zdelta=MIN(zwisoqpone/pqpone(jl),1.0_dp)!remin
                      IF (ABS(1.0_dp-zdelta).LT.cwisosec) zdelta=1.0_dp
                      IF (jt == i_HHO) zdelta = 1.0_dp
                      zwisoqinc=pqspone(jl)*zdelta
                   ELSE
                      zqinc=pqpone(jl) 
                      zwisoqinc=zwisoqpone
                   ENDIF

                   ! if some cloud-free areas exist: calculate normal
                   ! vapour and isotope value outside cloud (zqoutc, zwisoqoutc)
                   ! as the residual of the whole grid box - cloud values
                   IF (paclc(jl,kk).lt.(1.0_dp-1.0e-10_dp)) THEN
                      zqoutc=(pqpone(jl)-paclc(jl,kk)*zqinc)/(1.0_dp-paclc(jl,kk))
                      zwisoqoutc=(zwisoqpone-paclc(jl,kk)*zwisoqinc)/&
                           (1.0_dp-paclc(jl,kk))
                   ELSE
                      zqoutc=pqpone(jl)
                      zwisoqoutc=zwisoqpone
                   ENDIF

                   ! correct for negative values
                   IF ((zwisoqoutc).lt.0.0_dp) THEN
                      zwisoqoutc=zwisoqpone
                      zwisoqinc=zwisoqpone
                      zqinc=pqpone(jl)
                      zqoutc=pqpone(jl)
                   ENDIF

                   ! calculate isotope equilibration of precipitation
                   ! in the cloudy part of the grid box
                   zdenom=zwisofracliq(jl,jt)*twisoeqls*zqrain+zqinc
                   IF (zdenom.lt.cwisomin) THEN
                      zdeltab=0.0_dp
                   ELSE
                      zdeltab=twisoeqls*(zwisoqinc*zwisofracliq(jl,jt)*zqrain-zqinc*zwisoqrain)/zdenom
                      IF (jt == i_HHO) zdeltab =0.0_dp
                   ENDIF

                   ! calculate isotope equilibration of precipitation
                   ! in the cloud-free part of the grid box
                   zdenom=zwisofracliq(jl,jt)*twisoeqls*zqrain+zqoutc
                   IF (zdenom.lt.cwisomin) THEN
                      zdeltamb=0.0_dp
                   ELSE
                      zdeltamb=twisoeqls*(zwisoqoutc*zwisofracliq(jl,jt)*zqrain-zqoutc*zwisoqrain)/zdenom
                      IF (jt == i_HHO) zdeltamb = 0.0_dp
                   ENDIF

                   ! calculate isotope changes of the whole grid box
                   ! (cloudy and cloud-free part)
                   zdelta=paclc(jl,kk)*zdeltab+(1.0_dp-paclc(jl,kk))*zdeltamb

                ELSE
                   ! nothing happens if no positive vapour exist in next
                   ! time step (pqpone(jl).le.cwisomin) 
                   zdelta=0.0_dp
                ENDIF

                ! if no rain will exist in next time step (=all rain evaporates):
                ! - add isotope rain tendency to vapour without any
                !   fractionation  
             ELSE
                zdelta=-pwisorfl(jl,jt)/zrafac
             ENDIF

             ! add up isotopic changes caused by equilibration of rain with vapour

             pwisorfl(jl,jt)=pwisorfl(jl,jt)+zdelta*zrafac
             pwisoqte(jl,kk,jt)=pwisoqte(jl,kk,jt)-zdelta/ptmst

          ENDIF
       END DO
    END DO

  END SUBROUTINE wiso_cloudadj
  ! ----------------------------------------------------------------------
  
END MODULE messy_main_tools_wiso
