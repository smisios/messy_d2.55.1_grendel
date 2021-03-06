! ***********************************************************************
MODULE messy_gec
  ! ***********************************************************************
  ! Authors: Andreas Baumgaertner, 2014

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE
  INTRINSIC :: NULL

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'gec'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'

  ! GLOBAL CTRL NAMELIST VARIABLES
  REAL(dp),     PUBLIC :: sa  =  0.5_dp ! solar activity [0..1]
  INTEGER,      PUBLIC :: ionradon_method = 1
  INTEGER,      PUBLIC :: iongcr_method = 1
  INTEGER,      PUBLIC :: ioncloud_method = 1
  INTEGER,      PUBLIC :: ionaerosol_method = 1
  ! mz_se_20180530+
  ! if true then calculated ion concentration in this submodel
  LOGICAL,      PUBLIC :: lssion = .TRUE.
  LOGICAL,      PUBLIC :: ltotalipr = .TRUE.
  ! mz_se_20180530-


  INTEGER :: lev65km = 1 ! level above 65km -- no conductivity calculation above 
  INTEGER :: lev15km = 1 ! level above 15km -- no cloud conductivity calculation above 
  INTEGER :: current_lev_top = 1 ! vertical integration limits for current parametrization
  INTEGER :: current_lev_bottom =1 
  ! The regression numbers are obtained from GECsources/regression_v2.m
  REAL(dp) :: regress_land_a = 0.!0.025
  REAL(dp) :: regress_land_b = 22.
  REAL(dp) :: regress_ocean_a = 0.!0.026
  REAL(dp) :: regress_ocean_b = 22.

  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: carma_processed1a
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: carma_processed1b
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: carma_processed2

  PUBLIC :: GEC_init
  PUBLIC :: GEC_conductivity
  PUBLIC :: GEC_sourcecurrent
  PUBLIC :: GEC_solver
  PUBLIC :: GEC_read_nml_ctrl
  PUBLIC :: gec_clean

  REAL(dp), PARAMETER   :: el_charge = 1.60217646e-19 ! elementary (electron) charge (C)
  INTEGER, PARAMETER                 :: nCltypes=5, nClsizes=9, xsize=10000
  REAL(dp)    :: Rcolhat_to_Rcolnocld_eta50(nCltypes,nClsizes) &
               , Rcolhat_to_Rcolnocld_eta25(nCltypes,nClsizes) &
               , Rcolhat_to_Rcolnocld_eta10(nCltypes,nClsizes) &
               , Rcoldashcld_Rcolnocld_ratio(nCltypes,nClsizes) 
  REAL(dp)    :: x(xsize)            &
               , L(xsize)            &
               , contribution(xsize) &
               , clcontrib(nClsizes) &
               , cldfrac(nCltypes,nClsizes)
  REAL(dp), DIMENSION(nClsizes)      :: CLhoriz_all=(/0.5,1.,2.,6.,15.,40.,100.,300.,2000./) ! km

CONTAINS

  !==========================================================================
  SUBROUTINE gec_init(pref_mid, nlev)
    
    REAL(dp), INTENT(IN) :: pref_mid(:)
    INTEGER, INTENT(IN) :: nlev
    
    INTEGER :: cl, cli, clim1, k
    
    ! no calculation over 65 km / 10 Pa (0.1 hPa)
    DO k=1,nlev
       IF (pref_mid(k) > 10.) EXIT
       lev65km = k
    END DO
    
    ! no calculation over 15 km / 10000 Pa (100 hPa)
    DO k=1,nlev
       IF (pref_mid(k) > 10000.) EXIT
       lev15km = k
    END DO
    
    ! find heights for current parametrization
    DO k=1,nlev
       IF (pref_mid(k) > 20000.) EXIT
       current_lev_top = k
    END DO
    DO k=1,nlev
       IF (pref_mid(k) > 80000.) EXIT
       current_lev_bottom = k
    END DO
    
    Rcolhat_to_Rcolnocld_eta50 = &
         transpose(reshape( &
         (/147., 172., 218., 601., 1043., 1144., 1161., 1168., 1172., &
         110., 125., 144., 207., 445., 774., 892., 944., 970., &
         119., 138., 160., 274., 615., 812., 866., 890., 902., &
         124., 147., 174., 309., 863., 1313., 1451., 1512., 1543., &
         100., 102., 106., 120., 146., 210., 268., 296., 311. /) /100., &
         (/size(Rcolhat_to_Rcolnocld_eta50,2), size(Rcolhat_to_Rcolnocld_eta50,1) /) ))
    Rcolhat_to_Rcolnocld_eta25 = &
         transpose(reshape( &
         (/142., 166., 205., 431., 599., 626., 625., 625., 625., &
         109., 124., 142., 195., 336., 461., 500., 517., 521., &
         118., 135., 156., 241., 401., 467., 483., 489., 491., &
         122., 144., 170., 278., 565., 730., 776., 796., 800., &
         100., 102., 106., 118., 137., 170., 189., 199., 200. /) /100., &
         (/size(Rcolhat_to_Rcolnocld_eta50,2), size(Rcolhat_to_Rcolnocld_eta50,1) /) ))
    Rcolhat_to_Rcolnocld_eta10 = &
         transpose(reshape( &
         (/134., 153., 178., 262., 301., 303., 299., 298., 297., &
         107., 120., 135., 169., 220., 248., 255., 258., 259., &
         114., 129., 145., 188., 232., 245., 246., 247., 247., &
         118., 138., 158., 222., 315., 352., 360., 363., 364., &
         100., 101., 105., 114., 123., 133., 136., 138., 138.  /) /100., &
         (/size(Rcolhat_to_Rcolnocld_eta50,2), size(Rcolhat_to_Rcolnocld_eta50,1) /) ))
    
    Rcoldashcld_Rcolnocld_ratio = &
         transpose(reshape( &
         (/1172.,1172.,1172.,1172.,1172.,1172.,1172.,1172.,1172., &
         970.,970.,970.,970.,970.,970.,970.,970.,970., &
         902.,902.,902.,902.,902.,902.,902.,902.,902., &
         1543.,1543.,1543.,1543.,1543.,1543.,1543.,1543.,1543., &
         311.,311.,311.,311.,311.,311.,311.,311.,311. /) /100., &
         (/size(Rcoldashcld_Rcolnocld_ratio,2), size(Rcoldashcld_Rcolnocld_ratio,1) /) ))
    
    ! total cloud cover contribution of all cloud sizes in CLhoriz_all
    x=logspace(-6._dp,0._dp,10000)
    L=x*2000. ! km
    contribution=1.-x**(2.-1.7)
    DO cl=1,nClsizes
       cli=find(L>=CLhoriz_all(cl))
       IF (cl>1) then
          clim1=find(L>=CLhoriz_all(cl-1))        
       ELSE
          clim1=find(L>=0.)        
       END IF
       clcontrib(cl)=contribution(clim1)-contribution(cli)
    END DO
    
  END SUBROUTINE gec_init
  !==========================================================================

  ! =========================================================================
  SUBROUTINE GEC_conductivity ( pmid, pint              & ! INPUT
                              , tm1                     & ! INPUT
                              , density                 & ! INPUT
                              , numdensity              & ! INPUT
                              , mobility                & ! INPUT
                              , landfrac                & ! INPUT
                              , zm, zi, zmsea           & ! INPUT
                              , Rn222                   & ! INPUT
                              , cldnodeepcu             & ! INPUT
                              , philat, philon          & ! INPUT
                              , ilat, ilon              & ! INPUT
                              , nproma, kproma, nlev, jrow & ! INPUT
                              , dt                      & ! INPUT
                              , lsa                     & ! INPUT!mz_se_20170206
                              , con                     & ! OUTPUT
                              , con_nocld               & ! OUTPUT
                              , Rcoltilde               & ! OUTPUT
                              , Rcol_nocld              & ! OUTPUT
                              , ionrate_radon           & ! OUTPUT
                              , ionrate_direct          & ! OUTPUT
                              , ionrate_gcr             & ! OUTPUT
                              , ionrate_spe             & ! OUTPUT
                              , ionrate                 & ! OUTPUT
                              , ioncon                  & ! OUTPUT
                              , status                  & ! OUTPUT
                              )

    USE messy_gec_fcloudparam, only: f_cloudparam, fprime_cloudparam

    IMPLICIT NONE

    ! Input PARAMETERs
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: pmid, pint    ! 
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: tm1           ! 
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: density       ! air density midlevel [g cm-3]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: numdensity    ! number density midlevel [g cm-3]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: mobility      ! [cm2/V/s]
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: landfrac      ! 
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: zm, zi, zmsea ! geopotential heights [m]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: Rn222         ! 
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: cldnodeepcu   ! 
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: philat        ! latitude philat(ilat) [deg]
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: philon        ! longitude philon(ilon)
    INTEGER,     DIMENSION(:,:), INTENT(IN) :: ilat          ! latitude index (1:kproma,jrow)
    INTEGER,     DIMENSION(:,:), INTENT(IN) :: ilon          ! longitude index (1:kproma,jrow)
    INTEGER,                     INTENT(IN) :: nproma 
    INTEGER,                     INTENT(IN) :: kproma 
    INTEGER,                     INTENT(IN) :: nlev
    INTEGER,                     INTENT(IN) :: jrow
    REAL(dp),                    INTENT(IN) :: dt
    REAL(dp),                    INTENT(IN) :: lsa           ! solar cycle parameter !mz_se_20170206 
    ! Output parameters
    REAL(DP), DIMENSION(:,:),   INTENT(OUT) :: con       
    REAL(DP), DIMENSION(:,:),   INTENT(OUT) :: con_nocld       
    REAL(DP), DIMENSION(:),     INTENT(OUT) :: Rcoltilde       
    REAL(DP), DIMENSION(:),     INTENT(OUT) :: Rcol_nocld       
    REAL(DP), DIMENSION(:,:),   INTENT(INOUT) :: ionrate_radon ! made INOUT mz_se 
    REAL(DP), DIMENSION(:,:),   INTENT(OUT) :: ionrate_direct       
    REAL(DP), DIMENSION(:,:),   INTENT(INOUT) :: ionrate_gcr   ! made INOUT mz_se
    REAL(DP), DIMENSION(:,:),   INTENT(OUT) :: ionrate_spe       
    REAL(DP), DIMENSION(:,:),   INTENT(INOUT) :: ionrate       ! made INOUT mz_se
    REAL(DP), DIMENSION(:,:),   INTENT(INOUT) :: ioncon        ! made INOUT mz_se
    INTEGER,                    INTENT(OUT) :: status   ! error status

    ! Local
    real(dp), dimension(nproma)         :: cldtot, cldlow, cldmed, cldhgh
    REAL(dp) :: Rcolhat_cl_eta50(nCltypes,nClsizes) &
              , Rcolhat_cl_eta25(nCltypes,nClsizes) &
              , Rcolhat_cl_eta10(nCltypes,nClsizes) &
              , Rcoldash_cl(nCltypes,nClsizes)
    REAL(dp) :: Rcoltildetmp_eta50(nCltypes,nClsizes) &
              , Rcoltildetmp_eta25(nCltypes,nClsizes) &
              , Rcoltildetmp_eta10(nCltypes,nClsizes) &
              , Rcolsimpletmp(nCltypes)
    REAL(dp) :: gamma, gamma0
    REAL(DP) :: alpha(nproma,nlev)
    REAL(DP) :: ion_aer(nproma,nlev)
    REAL(DP) :: ion_cld(nproma,nlev)
    INTEGER :: i, k, level, jp, cl, cltype, j ! counters 


    INTRINSIC ABS

    status=1
!!$ ioncon =0._dp ! mz_se_20180530
    con = 0._dp
    con_nocld = 0._dp
    Rcoltilde = 0._dp
    Rcol_nocld = 0._dp

     ! Ionization rate from Radon 222Rn
    SELECT CASE (ionradon_method)
    CASE(0)
       ionrate_radon = 0._dp
    CASE(1)  
       CALL ionrate_radon_joeckel(pmid, pint, tm1, Rn222, nproma, kproma, nlev, jrow, dt, ionrate_radon)
    CASE(2)
       ionrate_radon(:,nlev) = landfrac(:)*10._dp * exp(-zm(:,nlev)/200._dp) 
    CASE(3) ! mz_se_20170206
       ! ionrate_radon not overwriten 
    END SELECT

    ! Ionization from alpha, beta, gamma, Rn220 (direct)
    IF(ionradon_method /= 3) & ! don't do this unless original GEC settings are used !mz_se_20170206
         & CALL ionrate_alpha_beta_gamma_Rn220(zm, landfrac, philat, ilat &
         , kproma, jrow, nlev, ionrate_direct)

    ! Ionization rate from GCR
    SELECT CASE (iongcr_method)
    CASE(0)
       ionrate_gcr = 0._dp
    CASE(1)
       CALL ionrate_gcr_tinsley(zmsea, philat, philon, ilat, ilon &
                              , lsa, nproma, kproma, nlev, jrow, ionrate_gcr)
       ! changed sa to lsa above to allow coupling to GCR modulation !mz_se_20170206
    CASE(2)
!       CALL ionrate_gcr_nairas(state,kproma,ionrate_gcr)
    CASE(3) ! added mz_se 20170206
      ! iongcr_method not overwrite 
    END SELECT


    ! Ion-ion recombination rate, Eq. 6a-c TZ06
    ! Could also do this with reference atmosphere!
    ! mz_se_20180530+
    ! calculate alpha only if the total steady state ions are calculeted 
    ! in this submodel ...
    IF(lssion) THEN
       ioncon =0._dp
       ! mz_se_20180530+
       DO i = 1, kproma
          DO k = lev65km, nlev  
             SELECT CASE (int(zmsea(i,k)))
                ! < 10 km
             CASE (:10000) 
                alpha(i,k) = 6e-8*(300./tm1(i,k))**.5 &
                     +1.702e-6*(300./tm1(i,k))**(-1.984) &
                     *(numdensity(i,k)/2.69e19)**(-.451)
                ! 10-20 km
             CASE (10001 : 20000) 
                alpha(i,k) = 6e-8*(300./tm1(i,k))**.5 &
                     +1.035e-6*(300./tm1(i,k))**4.374 &
                     *(numdensity(i,k)/2.69e19)**.769
                ! >= 20 km
             CASE (20001:) 
                alpha(i,k) = 6e-8*(300./tm1(i,k))**.5 &
                     +6.471e-6*(300./tm1(i,k))**(-.191) &
                     *(numdensity(i,k)/2.69e19)**(.901)
             END SELECT
          END DO
       END DO
    END IF! mz_se_20180530

    ! Ion-aerosol attachment 
    SELECT CASE (ionaerosol_method) 
    CASE(0) 
       ion_aer      = 0._dp
    CASE(1)
       ion_aer(:kproma,:) =  carma_processed1a(:kproma,:,jrow) & ! total_ion_aerosol_volc
                            +carma_processed1b(:kproma,:,jrow) & ! total_ion_aerosol_base
                            +carma_processed2(:kproma,:,jrow)
    END SELECT

!!$    ! Ion loss through clouds
!!$    SELECT CASE (ioncloud_method)
!!$    CASE(0)
!!$       ! ion_cld = 0._dp
!!$    CASE(1) 
!!$       ! nothing to do here
!!$    END SELECT

    ionrate_spe    = 0._dp

    ! Scale ionisation rate from STP to actual density
    DO i = 1, kproma  
       DO k = lev65km, nlev 
          IF(ionradon_method /= 3) THEN ! mz_se_20180530
             ionrate_radon(i,k)  = ionrate_radon(i,k)   *pmid(i,k)/100000._dp
             ionrate_direct(i,k) = ionrate_direct(i,k)  *pmid(i,k)/100000._dp
          END IF ! mz_se_20180530
          IF(iongcr_method /= 3) THEN ! mz_se_20180530
             ionrate_gcr(i,k)    = ionrate_gcr(i,k)     *pmid(i,k)/100000._dp
          END IF ! mz_se_20180530
          ! ionrate_spe is already correct

          ! sum all types
          IF(ltotalipr) & ! mz_se_20180530
          ionrate(i,k) = ionrate_radon(i,k)  &
                       + ionrate_direct(i,k) &
                       + ionrate_gcr(i,k)    &
                       + ionrate_spe(i,k)

       END DO
    END DO

    Rcoltilde = 0._dp
    Rcol_nocld = 0._dp

    ! cldtot, cldlow, cldmed, cldhgh
    CALL cldsav(jrow, kproma, nproma, nlev, nlev+1, cldnodeepcu, pmid, pint &
               , cldtot, cldlow, cldmed, cldhgh)

    ! Calculate ionconcentration, conductivity, column resistance
    DO i = 1, kproma  
       DO k = lev65km, nlev

          ! Ionconcentration: solution to eq. 3 from TZ06, measured in 1/cm3
          ! Conducivity: Eq. 4 TZ06,  measured in S/m

          ! no clouds
          IF(lssion) & ! mz_se_20180530
          ioncon(i,k) = (sqrt(4.*alpha(i,k)*ionrate(i,k) & 
                      + (ion_aer(i,k))**2.) -  ion_aer(i,k)) / (2.*alpha(i,k))          
          con_nocld(i,k) = ioncon(i,k) * el_charge * mobility(i,k) * 100._dp
          Rcol_nocld(i) = Rcol_nocld(i) + (zi(i,k)-zi(i,k+1))/con_nocld(i,k)


          ! cloud cover for cloud size bins for every cloud type
          ! for cumulus+stratocumulus, altostratus, altocumulus, nimbostratus, cirrus    
          ! paper: g = cldfrac; f(type) = cldlow,...; C_h(hi) = clcontrib(cl)
          DO cl=1,nClsizes
             cldfrac(1,cl) = (cldlow(i)*(11.5+12.1)/(11.5+12.1+1.8))*clcontrib(cl)
             cldfrac(2,cl) = (cldmed(i)*7.8/(8.6+7.8+2.4))*clcontrib(cl)
             cldfrac(3,cl) = (cldmed(i)*8.6/(8.6+7.8+2.4))*clcontrib(cl)
             cldfrac(4,cl) = (cldmed(i)*2.4/(8.6+7.8+2.4))*clcontrib(cl)
             cldfrac(5,cl) = (cldhgh(i))*clcontrib(cl)
          END DO

          Rcolhat_cl_eta50 = Rcolhat_to_Rcolnocld_eta50 *Rcol_nocld(i)
          Rcoldash_cl= Rcoldashcld_Rcolnocld_ratio*Rcol_nocld(i)
          DO cltype=1,nCltypes
             DO  cl=1,nClsizes        
                Rcoltildetmp_eta50(cltype,cl)=1./Rcolhat_cl_eta50(cltype,cl) *cldfrac(cltype,cl)
             END DO
          END DO
          Rcoltilde(i)=1./(sum(sum(Rcoltildetmp_eta50,2),1)+1./Rcol_nocld(i)*(1.-cldtot(i)))

          ! parametrize for every layer
          con(i,lev65km:nlev) = con_nocld(i,lev65km:nlev)
          IF (cldtot(i)>0.01) then
             DO j=3,1,-1
                gamma0=1.-1./maxval(cldnodeepcu(i,lev15km:nlev))
                gamma0=gamma0+1./(10.**REAL(j))
                CALL newtonsolve(f_cloudparam,fprime_cloudparam,gamma0,gamma,&
                              Rcoltilde(i),zi(i,lev15km:nlev)-zi(i,lev15km+1:nlev+1), &
                              (1.-cldnodeepcu(i,lev15km:nlev))*con_nocld(i,lev15km:nlev), &
                                  cldnodeepcu(i,lev15km:nlev) *con_nocld(i,lev15km:nlev),nlev-lev15km+1)
                IF (abs(gamma)<1e8) exit
             END DO
             IF (abs(gamma)>=1e8) gamma=gamma0
             con(i,lev15km:nlev)=    (1.-cldnodeepcu(i,lev15km:nlev))*con_nocld(i,lev15km:nlev) &
                                  +gamma*cldnodeepcu(i,lev15km:nlev) *con_nocld(i,lev15km:nlev)             
          END IF

       END DO ! k
    END DO ! i

    ! No error
    status=0

  END SUBROUTINE GEC_conductivity
  !================================================================================================

  !================================================================================================
  SUBROUTINE ionrate_alpha_beta_gamma_Rn220 (zm, landfrac, philat, ilat, kproma, jrow, nlev, ionpairs)

    REAL(dp), DIMENSION(:,:),   INTENT(IN)  :: zm             ! geopotential height [m]
    REAL(dp), DIMENSION(:),     INTENT(IN)  :: landfrac       ! land fraction
    REAL(dp), DIMENSION(:),     INTENT(IN)  :: philat         ! latitude [deg]
    INTEGER,  DIMENSION(:,:),   INTENT(IN)  :: ilat           ! 
    INTEGER,                    INTENT(IN)  :: kproma, jrow, nlev
    REAL(dp), DIMENSION(:,:),   INTENT(OUT) :: ionpairs       ! radon ion pair production rate [1/cm3/s]

    ! local
    INTEGER :: i,k

    ionpairs(:,:) = 0._dp

    ! alpha, beta, gamma, Rn220, see TZ06 [23]
    DO k=nlev-5,nlev ! only affects lowest 5 levels
       DO i=1, kproma
          IF (abs(philat(ilat(i,jrow))) < 60._dp) then
             ! 10 ionpairs cm-3 s-1, scale height 200m
             ionpairs(i,k) = landfrac(i) * 10._dp * exp(-zm(i,k)/200._dp) 
          ELSEIF (philat(ilat(i,jrow)) > 60._dp .and. philat(ilat(i,jrow)) < 70._dp) then
             ! 5 ionpairs cm-3 s-1
             ionpairs(i,k) = landfrac(i) * 5._dp * exp(-zm(i,k)/200._dp) 
          END IF
       END DO
    END DO

  END SUBROUTINE ionrate_alpha_beta_gamma_Rn220
  !================================================================================================

  !================================================================================================
  SUBROUTINE ionrate_gcr_tinsley ( zmsea              & ! INPUT
                                 , philat, philon     & ! INPUT
                                 , ilat, ilon         & ! INPUT
                                 , lsa                & ! INPUT
                                 , nproma, kproma, nlev, jrow & ! INPUT
                                 , ionpairs           & ! OUTPUT 
                                 )

    ! Galactic cosmic rate ionization rates parametrized as in 
    ! Tinsley, B. A., and L. Zhou (2006), Initial results of a global 
    ! circuit model with variable stratospheric and tropospheric
    ! aerosols, J. Geophys. Res., 111 (D10), D16205, doi:10.1029/2005JD006988.
    ! referred  to as TZ06
    
    USE messy_main_constants_mem, ONLY: pi, rearth=>radius_earth, r2d=>rtd, d2r=>dtr

    IMPLICIT NONE

    ! Input parameters
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: zmsea
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: philat ! latitude philat(ilat) [deg]
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: philon ! longitude philon(ilon)
    INTEGER,     DIMENSION(:,:), INTENT(IN) :: ilat   ! latitude index (1:kproma,jrow)
    INTEGER,     DIMENSION(:,:), INTENT(IN) :: ilon   ! longitude index (1:kproma,jrow)
    REAL(dp),                    INTENT(IN) :: lsa    ! local sa
    INTEGER,                     INTENT(IN) :: nproma, kproma 
    INTEGER,                     INTENT(IN) :: nlev
    INTEGER,                     INTENT(IN) :: jrow
    ! Output parameters
    REAL(dp),  DIMENSION(:,:),   INTENT(OUT) :: ionpairs   ! GCR ion pair production rate [1/cm3/s]

    ! Local
    REAL(dp)                    :: sol_cycle_factor,X1, X2, X3, s
    REAL(dp)                    :: thetar_min, thetar_max
    REAL(dp), DIMENSION(6)      :: z_equator, z_knee, theta_knee &
                                 , rate_equator, rate_knee
    REAL(dp), DIMENSION(kproma,6) :: z, q
    INTEGER                     :: i, k, l
    REAL(dp) :: alatm(nproma) ! magnetic latitude/longitude grid (radians)
    REAL(dp) :: magplat= 78.8_dp/180._dp*pi    ! latitude of geomagnetic pole [rad]
    REAL(dp) :: magplon=289.1_dp/180._dp*pi    ! longitude of geomagnetic pole [rad]
    REAL(dp) :: coslat, coslon, sinlat, sinlon
    ! matrices for transformation to geomagnetic latitude
    REAL(dp), DIMENSION(3,3) :: geo2mag1, geo2mag2, geo2magf 
    ! vectors to grid point in geographic and magnetic coordinates
    REAL(dp), DIMENSION(3)   :: vgeo, vmag 


    DO i=1, kproma
       ! GEO to MAG, e.g. Hapgood (1992,1997)
       ! Initialize coordinate transformation variables and matrices
       coslat=COS(philat(ilat(i,jrow))/180._dp*pi)
       coslon=COS(philon(ilon(i,jrow))/180._dp*pi)
       sinlat=SIN(philat(ilat(i,jrow))/180._dp*pi)
       sinlon=SIN(philon(ilon(i,jrow))/180._dp*pi)
       geo2mag1=RESHAPE( (/ COS(magplat-pi/2._dp), 0._dp, -SIN(magplat-pi/2._dp),  &  ! column 1
                                0._dp,             1._dp, 0._dp,                   &  ! column 2
                          SIN(magplat-pi/2._dp), 0._dp, COS(magplat-pi/2._dp)/),   &  ! column 3 
            (/3,3/) )
       geo2mag2=RESHAPE( (/ COS(magplon), -SIN(magplon), 0._dp,    &  ! column 1
                            SIN(magplon), COS(magplon),  0._dp,    &  ! column 2
                            0._dp        ,       0._dp,  1._dp /), &  ! column 3 
                            (/3,3/) )
       geo2magf=MATMUL(geo2mag1,geo2mag2)
       vgeo(1)=rearth*coslat*coslon
       vgeo(2)=rearth*coslat*sinlon
       vgeo(3)=rearth*sinlat
       vmag=MATMUL(geo2magf,vgeo)
       alatm(i)=ATAN2(vmag(3),SQRT(vmag(1)**2+vmag(2)**2))
       ! maglon(i)=ACOS(vmag(1)/SQRT(vmag(1)**2+vmag(2)**2)) ! magn. long. not needed
       ! IF (vmag(2)<0._dp) maglon(i)=360._dp-maglon(i) ! magn. long. not needed
    END DO

    ionpairs = 0._dp 

 
    ! Table 1 from TZ06

    ! Column 2: Altitude of Level at Equator (km)
    z_equator(1) = 0._dp      ! sea level                       
    z_equator(2) = 6.5_dp     ! intermediate altitude a         
    z_equator(3) = 9.0_dp     ! intermediate altitude b         
    z_equator(4) = 11.0_dp    ! intermediate altitude c         
    z_equator(5) = 16.0_dp    ! height of maximum production    
    z_equator(6) = 1000._dp   ! above                           
    z_equator(:) = z_equator(:)*1000.  ! km --> m

    ! Column 3: Altitude of Level at Knee (km)
    z_knee(1) = 0._dp
    z_knee(2) = 9.8_dp
    z_knee(3) = 16.0_dp
    z_knee(4) = 21.0_dp
    z_knee(5) = 32.0_dp
    z_knee(6) = 1000._dp 
    z_knee(:) = z_knee(:)*1000. ! km --> m

    ! Column 4+5: Knee latitude 
    theta_knee(1) = 50.-(50.-49.)*lsa 
    theta_knee(2) = 52.-(52.-50.)*lsa  
    theta_knee(3) = 58.-(58.-55.)*lsa
    theta_knee(4) = 60.-(60.-56.)*lsa
    theta_knee(5) = 63.-(63.-56.)*lsa
    theta_knee(6) = 63.-(63.-56.)*lsa
    theta_knee(:) = theta_knee(:)*d2r

    ! Column 7+8: Ion production at Equator [1/cm3/s]
    rate_equator(1) = 1.4 *(1.-0.040*lsa) 
    rate_equator(2) = 13.5*(1.-0.042*lsa) 
    rate_equator(3) = 34. *(1.-0.044*lsa) 
    rate_equator(4) = 64. *(1.-0.047*lsa) 
    rate_equator(5) = 98. *(1.-0.050*lsa) 
    rate_equator(6) = 45. *(1.-0.060*lsa) 

    ! Column 9+10: Ion production at knee
    rate_knee(1) = 2.2 *(1.-0.10*lsa) 
    rate_knee(2) = 145.*(1.-0.20*lsa) 
    rate_knee(3) = 325.*(1.-0.30*lsa) 
    rate_knee(4) = 435.*(1.-0.45*lsa) 
    rate_knee(5) = 500.*(1.-0.50*lsa) 
    rate_knee(6) = 550.*(1.-0.60*lsa) 

    ! GCR Ionization Reference Altitudes

    ! loop over all rows
    DO i=1, kproma   
       DO l=1, 6 ! 6 reference altitudes
          IF (abs(alatm(i))<theta_knee(l)) then
             z(i,l)=z_equator(l)+(z_knee(l)-z_equator(l)) &
                  *sin(alatm(i))**4._dp/sin(theta_knee(l))**4._dp
             q(i,l) = rate_equator(l)+(rate_knee(l)-rate_equator(l)) &
                  *sin(alatm(i))**4._dp/sin(theta_knee(l))**4._dp
          ELSE 
             z(i,l) = z_knee(l)
             q(i,l) = rate_knee(l);
          END IF
       END DO
    END DO

    ! loop over all levels
    DO k=nlev,lev65km,-1 
       ! loop over all rows
       DO i=1, kproma   
          IF (zmsea(i,k) < z(i,2)) then 
             s = (z(i,2)-z(i,1))/log(q(i,2)/q(i,1)) ! scale height
             ionpairs(i,k) = q(i,2)*exp((zmsea(i,k)-z(i,2))/s) 
          ELSEIF (zmsea(i,k) < z(i,3)) then
             s = (z(i,3)-z(i,2))/log(q(i,3)/q(i,2))
             ionpairs(i,k) = q(i,3)*exp((zmsea(i,k)-z(i,3))/s) 
          ELSEIF (zmsea(i,k) < z(i,4)) then
             s = (z(i,4)-z(i,3))/log(q(i,4)/q(i,3))
             ionpairs(i,k) = q(i,4)*exp((zmsea(i,k)-z(i,4))/s)
          ELSEIF (zmsea(i,k) < z(i,5)) then
             s = (z(i,5)-z(i,4))/(log(q(i,5)/q(i,4)))**.5
             ionpairs(i,k) = q(i,5)*exp(-((zmsea(i,k)-z(i,5))/s)**2.)
          ELSE            
             s = (z(i,5)-z(i,4))/(log(q(i,5)/q(i,4)))**.5
             ionpairs(i,k) = q(i,6)+(q(i,5)-q(i,6))*exp(-((zmsea(i,k)-z(i,5))/s)**2.)
          END IF
       END DO
    END DO

  END SUBROUTINE ionrate_gcr_tinsley
  !================================================================================================

  !================================================================================================
  SUBROUTINE ionrate_radon_joeckel (pmid, pint, tm1, Rn222, nproma, kproma, nlev, jrow, dt, ionpairs)

    use messy_main_constants_mem,    ONLY: boltz=>k_B    & ! Boltzmann constant J/K/molecule
                                         , mwdry=>M_air    !  mean molar mass of air

    REAL(dp), DIMENSION(:,:), INTENT(IN)    :: pmid, pint, tm1, Rn222     ! 
    INTEGER,                  INTENT(IN)    :: nproma, kproma, jrow, nlev
    REAL(dp),                 INTENT(IN)    :: dt           ! time step
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: ionpairs ! radon ion pair production rate [1/cm3/s]


    ! local
!!$    REAL(dp)              :: dp_sfc(nproma)      ! delta-p at surface
    REAL(dp)              :: decay(nproma,nlev),  Rn222flux_mmrs(nproma)
    INTEGER               :: rn_flux
    INTEGER               :: i, k
    REAL(dp), PARAMETER   :: mwrn222 = 222.018  ! molecular weight of Rn222

    decay(:,:) = 0._dp
    ionpairs(:,:) = 0._dp
    Rn222flux_mmrs(:) = 0._dp

!!$    dp_sfc(1:kproma) = abs( pint(1:kproma,nlev) &
!!$         - pint(1:kproma,nlev+1) )

!!$    CALL flux2mmrs(Rn222flux_mmrs(1:kproma), Rn222flux(1:kproma,jrow,mon), & 
!!$         dp_sfc(1:kproma),mwrn222)


    ! Ionisation rate: just decay
    DO i = 1, kproma
       DO k = lev65km, nlev
          CALL int_radon(decay(i,k), Rn222(i,k), dt)
       END DO
    END DO
    ! convert decay from kg/kg/s to 1/cm3/s
    DO k = lev65km, nlev
       decay(1:kproma,k) = decay(1:kproma,k)*mwdry/mwrn222 * pmid(1:kproma,k) / boltz / tm1(1:kproma,k) *1.E-06
       ! convert Bq m^(-2) s^(-1) to ion pairs m^(-2) s^(-1)
       ionpairs(1:kproma,k) = (-1.)*decay(1:kproma,k) * 3.5E+05
    END DO

  contains

    ! ----------------------------------------------------------------------
    elemental SUBROUTINE int_radon(y, x, dt, q)

      ! Radon MODULE ROUTINE (CORE)
      ! MAIN INTERFACE ROUTINE OF Radon DIAGNOSTIC
      !
      ! INPUT:
      !        dt         time step length
      !        q          Radon source strength (mol/mol/s)
      !        x          tracer, e.g., 222-Rn at t=t0
      !
      ! OUTPUT:
      !        y          Radon tENDency for t0-:-t0+dt  (?/s)
      !
      ! Note:
      !       [q] and [x] must be in compatible units,
      !           e.g. mol/mol/s and mol/mol
      !
      ! Author: Patrick Joeckel, MPICH, Jul 2003

      implicit none

      ! I/O
      REAL(dp), INTENT(OUT)             :: y   ! Radon tENDency
      REAL(dp), INTENT(IN)              :: x   ! Radon
      REAL(dp), INTENT(IN)              :: dt  ! time step in seconds
      REAL(dp), INTENT(IN), optional    :: q   ! Radon source

      ! LOCAL
      ! HALF-LIFE
      ! e.g., 222Rn-half-life 3.8d
      REAL(dp), PARAMETER :: t2_Rn =  &
           3.8_dp*24.0_dp*3600.0_dp             ! 222Rn -> 218Po

      ! PRE-CALCULTED FOR CHAIN INTEGRATION
      REAL(dp) :: zq  ! source term
      REAL(dp) :: tau ! decay time

      IF (present(q)) then
         zq = q
      ELSE
         zq = 0.0
      END IF

      ! TENDENCY:
      ! DIFFERENTIAL EQUATION (HERE ANALYTICALLY SOLVED)
      ! dc/dt = -(1/tau)*c + q ; c(0) = c0
      !    with  tau = (t2_Rn / ln(2))
      !         --> c(t) = exp(-t/tau)*(c0-q*tau+q*tau*exp(t/tau))
      !             c(t) = c0*(exp(-t/tau)-1) + q*tau*(1-exp(-t/tau))
      !             --------------------------------------
      !         ==> c(t) - c0 = (c0-q*tau)*(exp(-t/tau)-1)
      !             --------------------------------------

      tau = t2_Rn / log(2.0_dp)

      y = ( x - zq * tau ) * ( exp(-dt/tau ) - 1._dp ) / dt

    END SUBROUTINE int_radon
    ! ----------------------------------------------------------------------

  END SUBROUTINE ionrate_radon_joeckel
  !================================================================================================

  ! =========================================================================
  SUBROUTINE GEC_sourcecurrent(landfrac, zi, zmmu, kproma, nlev, hourUT, current, currentcorr)
    ! source current parametrization based on updraft mass flux

    REAL(dp), DIMENSION(:),     INTENT(IN)  :: landfrac       ! land fraction
    REAL(dp), DIMENSION(:,:),   INTENT(IN)  :: zmmu           ! updraft mass flux
    REAL(dp), DIMENSION(:,:),   INTENT(IN)  :: zi
    INTEGER,                    INTENT(IN)  :: kproma, nlev
    INTEGER,                    INTENT(IN)  :: hourUT
    REAL(dp), DIMENSION(:),     INTENT(OUT) :: current        ! upward current
    REAL(dp), DIMENSION(:),     INTENT(OUT) :: currentcorr    ! 

    ! The following is obtained from GECsources/carnegie_correct_CESM.m
    REAL(dp) :: Carnegie_correction_land(24)  = (/0.7601, 0.7279, 0.6490, &
         0.6334, 0.6552, 0.7254, 0.8204, 0.8918, 0.8950, 0.8834, 0.9133, &
         1.0153, 1.1656, 1.2728, 1.2715, 1.1650, 1.1432, 1.1823, 1.2892, &
         1.3766, 1.3670, 1.2038, 1.0292, 0.8690 /)
    REAL(dp) :: Carnegie_correction_ocean(24) = (/0.9584, 0.9619, 0.9568, &
         0.9782, 0.9808, 1.0037, 1.0030, 1.0133, 1.0182, 1.0207, 1.0203, &
         1.0338, 1.0171, 1.0139, 1.0134, 1.0017, 0.9879, 0.9935, 0.9769, &
         0.9843, 0.9797, 0.9833, 0.9797, 0.9617 /)
    REAL(dp) :: weights(nlev)     
    INTEGER  :: i, k, itim, cl, Cltype
    REAL(dp) :: zmmu_int


    DO i=1,kproma
       weights(:)=0._dp
       zmmu_int=0._dp
       DO k=current_lev_top,current_lev_bottom
          weights(k)=zi(i,k)-zi(i,k+1)
          zmmu_int=zmmu_int+zmmu(i,k)*weights(k)
       END DO
       zmmu_int=zmmu_int/sum(weights)
       IF (landfrac(i)>0.5_dp) then
          current(i)=regress_land_a+regress_land_b*zmmu_int
          currentcorr(i)=current(i)*Carnegie_correction_land(int(hourUT)+1)
       ELSE
          current(i)=regress_ocean_a+regress_ocean_b*zmmu_int
          currentcorr(i)=current(i)*Carnegie_correction_ocean(int(hourUT)+1)
       END IF
    END DO


  END SUBROUTINE GEC_sourcecurrent
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE GEC_solver
    IMPLICIT NONE
  END SUBROUTINE GEC_solver
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE GEC_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ sa &
         , ionradon_method, iongcr_method, ioncloud_method & ! mz_se_20180530
         , ionaerosol_method, lssion, ltotalipr              ! mz_se_20180530

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='gec_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE GEC_read_nml_ctrl
  ! =========================================================================

  !--------------------------------------------------------------------------
  SUBROUTINE gec_clean

    IMPLICIT NONE

    INTRINSIC ASSOCIATED

    ! DATA
!!$    IF (ASSOCIATED(Ap_data)) DEALLOCATE(Ap_data)
!!$    NULLIFY(Ap_data)

  END SUBROUTINE gec_clean
  !--------------------------------------------------------------------------

  !===============================================================================
  subroutine cldsav(lchnk, ncol, pcols, pver, pverp, cld, pmid, pint &
                    , cldtot, cldlow, cldmed, cldhgh)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute total & 5 levels of cloud fraction assuming maximum-random overlap.
    ! Pressure ranges for the 5 cloud levels are specified.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: W. Collins
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: lchnk                ! chunk identifier
    integer, intent(in) :: ncol                 ! number of atmospheric columns
    integer, intent(in) :: pcols, pver, pverp

    real(dp), intent(in) :: cld(:,:)     ! Cloud fraction
    real(dp), intent(in) :: pmid(:,:)    ! Level pressures
    real(dp), intent(in) :: pint(:,:)    ! Level pressures
    real(dp)  :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
    !    maximally overlapped region.
    !    0->pmxrgn(i,1) is range of pressure for
    !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
    !    2nd region, etc

    integer :: nmxrgn(pcols)        ! Number of maximally overlapped regions
    !
    ! Output arguments
    !
    real(dp), intent(out) :: cldtot(:)       ! Total random overlap cloud cover
    real(dp), intent(out) :: cldlow(:)      ! Low random overlap cloud cover
!!$    real(dp), intent(out) :: cldlow2(:)      ! Low random overlap cloud cover
    real(dp), intent(out) :: cldmed(:)      ! Middle random overlap cloud cover
!!$    real(dp), intent(out) :: cldmed2(:)      ! Middle random overlap cloud cover
    real(dp), intent(out) :: cldhgh(:)       ! High random overlap cloud cover

    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: i,k, ityp                  ! Longitude,level indices
    integer :: irgn(pcols)          ! Max-overlap region index
    integer :: max_nmxrgn           ! maximum value of nmxrgn over columns
    real(dp) :: clrsky(pcols)       ! Max-random clear sky fraction
    real(dp) :: clrskymax(pcols)    ! Maximum overlap clear sky fraction
    !------------------------------Parameters-------------------------------
    real(dp) :: plowmax             ! Max prs for low cloud cover range
    real(dp) :: plowmin             ! Min prs for low cloud cover range
!!$    real(dp) :: plow2max             ! Max prs for low cloud cover range
!!$    real(dp) :: plow2min             ! Min prs for low cloud cover range
    real(dp) :: pmedmax             ! Max prs for mid cloud cover range
    real(dp) :: pmedmin             ! Min prs for mid cloud cover range
!!$    real(dp) :: pmed2max             ! Max prs for mid cloud cover range
!!$    real(dp) :: pmed2min             ! Min prs for mid cloud cover range
    real(dp) :: phghmax             ! Max prs for hgh cloud cover range
    real(dp) :: phghmin             ! Min prs for hgh cloud cover range
    !
    !CLvert_all1    =[ 1 3   2   2    8  ]*1e3; % cloud bottom
    !CLvert_all2    =[ 2 5   3   5    9.5]*1e3; % cloud top
    !1=900 2=800 3=700 5=500 8=360 9.5=280
!!$    parameter (plowmax = 120000._dp,plowmin = 80000._dp, & ! 1-2
!!$               plow2max =  60000._dp,plow2min = 45000._dp, & ! 4-5
!!$               pmedmax =  80000._dp,pmedmin = 70000._dp, & ! 2-3
!!$               pmed2max =  70000._dp,pmed2min = 60000._dp, & ! 3-4
!!$               phghmax =   45000._dp,phghmin =  20000._dp)   ! 8-9.5
    parameter (plowmax = 120000._dp,plowmin = 68000._dp, & !
               pmedmax =  68000._dp,pmedmin = 40000._dp, & !
               phghmax =  40000._dp,phghmin = 10000._dp)   !


!!$    real(dp) ptypmin(6)
!!$    real(dp) ptypmax(6)
!!$    data ptypmin /phghmin,  plowmin, plow2min, pmedmin, pmed2min, phghmin/
!!$    data ptypmax /plowmax, plowmax, plow2max, pmedmax, pmed2max, phghmax/
    real(dp) :: ptypmin(4)
    real(dp) :: ptypmax(4)
    data ptypmin /phghmin,  plowmin, pmedmin, phghmin/
    data ptypmax /plowmax,  plowmax, pmedmax, phghmax/

    cldtot = 0._dp
    cldlow = 0._dp
    cldmed = 0._dp
    cldhgh = 0._dp

    call cldovrlap(lchnk, ncol, pint, cld, nmxrgn, pmxrgn)
    !
    !-----------------------------------------------------------------------
    !
    ! Initialize region number
    !

    max_nmxrgn = -1
    do i=1,ncol
       max_nmxrgn = max(max_nmxrgn,nmxrgn(i))
    end do

    do ityp = 1, 4
       irgn(1:ncol) = 1
       do k =1,max_nmxrgn-1
          do i=1,ncol
             if (pmxrgn(i,irgn(i)) < ptypmin(ityp) .and. irgn(i) < nmxrgn(i)) then
                irgn(i) = irgn(i) + 1
             end if
          end do
       end do
       !
       ! Compute cloud amount by estimating clear-sky amounts
       !
       clrsky(1:ncol)    = 1.0_dp
       clrskymax(1:ncol) = 1.0_dp
       do k = 1, pver
          do i=1,ncol
             if (pmid(i,k) >= ptypmin(ityp) .and. pmid(i,k) <= ptypmax(ityp)) then
                if (pmxrgn(i,irgn(i)) < pmid(i,k) .and. irgn(i) < nmxrgn(i)) then
                   irgn(i) = irgn(i) + 1
                   clrsky(i) = clrsky(i) * clrskymax(i)
                   clrskymax(i) = 1.0_dp
                endif
                clrskymax(i) = min(clrskymax(i),1.0_dp-cld(i,k))
             endif
          end do
       end do
       if (ityp == 1) cldtot(1:ncol)  = 1.0_dp - (clrsky(1:ncol) * clrskymax(1:ncol))
       if (ityp == 2) cldlow(1:ncol)  = 1.0_dp - (clrsky(1:ncol) * clrskymax(1:ncol))
!!$       if (ityp == 3) cldlow2(1:ncol) = 1.0_dp - (clrsky(1:ncol) * clrskymax(1:ncol))
       if (ityp == 3) cldmed(1:ncol)  = 1.0_dp - (clrsky(1:ncol) * clrskymax(1:ncol))
!!$       if (ityp == 5) cldmed2(1:ncol) = 1.0_dp - (clrsky(1:ncol) * clrskymax(1:ncol))
       if (ityp == 4) cldhgh(1:ncol)  = 1.0_dp - (clrsky(1:ncol) * clrskymax(1:ncol))
    end do

    contains
!===============================================================================
  subroutine cldovrlap(lchnk   ,ncol    ,pint    ,cld     ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Partitions each column into regions with clouds in neighboring layers.
! This information is used to implement maximum overlap in these regions
! with random overlap between them.
! On output,
!    nmxrgn contains the number of regions in each column
!    pmxrgn contains the interface pressures for the lower boundaries of
!           each region! 
! Method: 

! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------

!
! Input arguments
!
    integer, intent(in) :: lchnk                ! chunk identifier
    integer, intent(in) :: ncol                 ! number of atmospheric columns

    real(dp), intent(in) :: pint(:,:)   ! Interface pressure
    real(dp), intent(in) :: cld(:,:)     ! Fractional cloud cover
!
! Output arguments
!
    integer,  intent(out) :: nmxrgn(:)    ! Number of maximally overlapped regions
    real(dp), intent(out) :: pmxrgn(:,:)  ! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
!
!---------------------------Local variables-----------------------------
!
    integer i                    ! Longitude index
    integer k                    ! Level index
    integer n                    ! Max-overlap region counter

    real(dp) pnm(pcols,pverp)    ! Interface pressure

    logical cld_found            ! Flag for detection of cloud
    logical cld_layer(pver)      ! Flag for cloud in layer
!
!------------------------------------------------------------------------
!
    nmxrgn=0._dp
    pmxrgn=0._dp

    do i = 1, ncol
       cld_found = .false.
       cld_layer(:) = cld(i,:) > 0.0_dp
       pmxrgn(i,:) = 0.0_dp
       pnm(i,:)=pint(i,:)*10._dp
       n = 1
       do k = 1, pver
          if (cld_layer(k) .and.  .not. cld_found) then
             cld_found = .true.
          else if ( .not. cld_layer(k) .and. cld_found) then
             cld_found = .false.
             if (count(cld_layer(k:pver)) == 0) then
                exit
             endif
             pmxrgn(i,n) = pnm(i,k)
             n = n + 1
          endif
       end do
       pmxrgn(i,n) = pnm(i,pverp)
       nmxrgn(i) = n
    end do

    return
  end subroutine cldovrlap

  end subroutine cldsav

 FUNCTION LOGSPACE (A, B, N)
    integer, intent(in) :: N
    real(dp) :: LOGSPACE(N)
    real(dp), intent(in) ::  A, B
    real(dp) :: dindgen(N)
    integer :: i
    do i=1,N
       dindgen(i)=real(i-1,dp)
    end do    
    LOGSPACE = dindgen/(N-1.0) * (B-A) + A
    LOGSPACE = 10.**LOGSPACE        
  END FUNCTION LOGSPACE

  function find(array)
    logical, dimension(:), intent(in) :: array
    integer :: find
    do find=1,size(array)
       if (array(find)) return
    end do
  end function find

  subroutine newtonsolve(f, fp, x0, x, R,d,a,b,pver)

    ! Estimate the zero of f(x) using Newton's method. 
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!) 
    !   the number of iterations iters

    implicit none
    integer, intent(in) :: pver
    real(dp), intent(in) :: x0
    real(dp), intent(in) :: R,d(pver),a(pver),b(pver)
    real(dp), external :: f, fp
    logical :: debug=.FALSE.
    real(dp), intent(out) :: x
!!$    integer, intent(out) :: iters

    integer, parameter :: maxiter = 20
    real(dp), parameter :: tol = 1.e15
    real(dp), parameter :: toolarge = 1.e8

    ! Declare any local variables:
    real(dp) :: deltax, fx, fxprime
    integer :: k


    ! initial guess
    x = x0

!!$    if (debug) then
!!$        print 11, x
!!$ 11     format('Initial guess: x = ', e22.15)
!!$        endif

    ! Newton iteration to find a zero of f(x) 
    do k=1,maxiter

       ! evaluate function and its derivative:
       fx = f(x,R,d,a,b,pver)
       fxprime = fp(x,R,d,a,b,pver)

        if (debug) then
            write(*,*) 'newtonsolve, after',k,x,fx
! 12         format('newtonsolve, after', i3, ' iterations, x = ', e22.15,', ',e22.15)
         endif

       if (abs(fx) < tol) then
          exit  ! jump out of do loop
       endif
       if (abs(x) > toolarge) then
          exit  ! jump out of do loop
       endif

       ! compute Newton increment x:
       deltax = fx/fxprime

       ! update x:
       x = x - deltax

    enddo


    if (k > maxiter) then
       ! might not have converged

       fx = f(x,R,d,a,b,pver)
       if (abs(fx) > tol) then
          !write(*,*) '*** Warning: newtonsolve has not yet converged'
          x=1e12
       endif
    endif

!!$    ! number of iterations taken:
!!$    iters = k-1


  end subroutine newtonsolve


  ! ***********************************************************************
END MODULE messy_gec
! ***********************************************************************
