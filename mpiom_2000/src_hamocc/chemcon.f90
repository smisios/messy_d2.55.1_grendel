      SUBROUTINE CHEMCON(modal,kpie,kpje,kpke,psao,ptho,       &
     &           pddpo,pdlxp,pdlyp,ptiestu,kplmonth)
!**********************************************************************
!
!**** *CHEMCON* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - rename: tracer(i,j,k)=ocetra(i,j,k,itracer)
!     - interfacing with ocean model.
!     - use of zonal mean temperature/salinity.
!     akw3(kpje,kpke),akb3(kpje,kpke),ak13(kpje,kpke),ak23(kpje,kpke)
!     aksp(kpje,kpke) again 3d fields! (EMR suggestion)
!
!     Patrick Wetzel,    *MPI-Met, HH*    15.04.02
!     - rename to CHEMCON
!     - chemc(i,j,7) to chemcm(i,j,8,12)
!
!     S.Lorenz/JO.Beismann, OpenMP parallel    *MPI-Met, HH*  24.08.07
!
!     Purpose
!     -------
!     Computes chemical constants in the surface layer (chemcm)
!     and in the water column (ak13, ak23, akb3, aksp)
!
!     Method
!     -------
!     .
!
!     *CALL*       *CHEMCON(modal,kpie,kpje,kpke,psao,ptho,
!                          pddpo,pdlxp,pdlyp,ptiestu,kplmonth)*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *modal*   - .
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *ptiestu* - level depths [m].
!     *INTEGER* *kplmonth*  - month in ocean model
!
!     Externals
!     ---------
!     .
!**********************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      use mo_param1_bgc

      USE mo_control_bgc

!      USE mo_param1, ONLY: ie_g
      USE mo_commoau1, ONLY: tmelt
      USE mo_parallel

implicit none

      INTEGER :: i,j,k,kplmonth,kpie,kpje,kpke,modal,kchem,kmon

      REAL(wp) :: pddpo(kpie,kpje,kpke)
      REAL(wp) :: ptho (kpie,kpje,kpke)
      REAL(wp) :: psao (kpie,kpje,kpke)
      REAL(wp) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL(wp) :: ptiestu(kpke+1)

      REAL(wp) :: ZERO,TENM7,SMICR,THOUSI,PERC,FOURTH,THIRD
      REAL(wp) :: HALF,ONE,TWO,TEN
      REAL(wp) :: devkb,devk1t,devk2t,devkbt,devkst,devks
      REAL(wp) :: rgas,bor1,bor2,oxyco
      REAL(wp) :: c00,c01,c02,c03,c04,c05,c10,c11,c12,c13,c20,c21,c22,c23
      REAL(wp) :: cb0,cb1,cb2,cb3,cw0,cw1,cw2
      REAL(wp) :: ox0,ox1,ox2,ox3,ox4,ox5,ox6
      REAL(wp) :: an0,an1,an2,an3,an4,an5,an6
      REAL(wp) :: akcc1, akcc2, akcc3, akcc4
      REAL(wp) :: salchl, arafra, calfra,sucall,aracal
      REAL(wp) :: devk1, devk2, t,q,s,cl
      REAL(wp) :: cek0, ckb, ck1, ck2, ckw, oxy, ani
      REAL(wp) :: ak1, ak2, akb, akw, ak0, aksp0
      REAL(wp) :: p,cp,tc!,bor, rho, rrr
      REAL(wp) :: a1,a2,a3,b1,b2,b3,atn2o,rs

       WRITE(io_stdo_bgc,*)  ' CHEMCON is called ',                     &
      &               'with modal,kpie,kpje,kpke,kplmonth,ptiestu:  ',  &
      &                modal,kpie,kpje,kpke,kplmonth,ptiestu

! kplmonth is 0 when chemcon is called from INI_BGC
      if (kplmonth .eq. 0) then
         kplmonth = 1
      endif
!
!     -----------------------------------------------------------------
!*         1. SET HALF PRECISION CONSTANTS
!             --- ---- --------- ---------
!
      ZERO=0._wp
      TENM7=10._wp**(-7.0_wp)
      SMICR=1.E-6_wp
      THOUSI=1._wp/1000._wp
      PERC=0.01_wp
      FOURTH=0.25_wp
      THIRD=1._wp/3._wp
      HALF=0.5_wp
      ONE=1._wp
      TWO=2._wp
      TEN=10._wp
!
!     -----------------------------------------------------------------
!*         3. SET CONVERSION FACTOR SALINITY -> CHLORINITY
!             ------ ------- -- ---- ------- -- ----------
!             (AFTER WOOSTER ET AL., 1969)
!
      SALCHL=1._wp/1.80655_wp
!
!     -----------------------------------------------------------------
!*         4. SET ZERO DEG CENTIGRADE AT KELVIN SCALE
!             --- ---- --- ---------- -- ------ -----
!
!
!     -----------------------------------------------------------------
!*         5. SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG)
!             (SEE BROECKER A. PENG, 1982, P. 26)
!             ([CA++](MOLES/KG)=1.026E-2*(S/35.) AFTER
!             CULKIN(1965), CF. BROECKER ET AL. 1982)
!             ------------- --- -------- -- --- -----
!
      calcon = 1.03e-2_wp
!
!     -----------------------------------------------------------------
!*         6. SET COEFFICIENTS FOR APPARENT SOLUBILITY EQUILIBRIUM
!             OF CALCITE (INGLE, 1800, EQ. 6)
!             -- ------- ------- ----- --- ----------- -----------
!
       akcc1 = -34.452_wp
       akcc2 = -39.866_wp
       akcc3 = 110.21_wp
       akcc4 = -7.5752e-6_wp
!
!     -----------------------------------------------------------------
!*        6A. SET FRACTION OF ARAGONITE IN BIOGENIC CACO3 PARTICLES
!             --- -------- --------- -- -------- ----- ---------
!
       arafra = 0._wp   ! js: obsolete
!
!     -----------------------------------------------------------------
!*        6B. FRACTION OF CALCITE IN BIOGENIC CACO3 PARTICLES
!             -------- ------- -- -------- ----- ---------
!
       calfra = 1._wp - arafra
       SUCALL=ARAFRA+CALFRA

!
!     -----------------------------------------------------------------
!*         7. FACTOR TO GET APPARENT SOLUBILITY PRODUCT FOR
!             ARAGONIT BY MULTIPLICATION WITH APPARENT SOLUBILITY
!             PRODUCT FOR CALCIT (CF. BERNER, 1976,
!             OR BROECKER ET AL., 1982)
!             -- -------- -- ---- ----- -- ------ --- ------- -----
!
      aracal = arafra * 1.45_wp + calfra
!     WRITE(io_stdo_bgc,*) 'ARACAL=',ARACAL
!
!     -----------------------------------------------------------------
!*         8. SET COEFFICIENTS FOR SEAWATER PRESSURE CORRECTION OF

!             (AFTER CULBERSON AND PYTKOWICZ, 1968, CF. BROECKER
!             ET AL., 1982)
!             ------- -------- --- ---------- ----- --- -------- -
!
      devk1 = 24.2_wp
      devk2 = 16.4_wp
      devkb = 27.5_wp
      devk1t = 0.085_wp
      devk2t = 0.04_wp
      devkbt = 0.095_wp
!
!     -----------------------------------------------------------------
!*         9. SET COEFFICIENTS FOR PRESSURE CORRECTION OF SOLUBILITY
!             PRODUCT OF CACO3 CORRESPONDING TO ARAGONITE/CALCITE
!             RATIO AFTER EDMOND AND GIESKES (1970),
!             P. 1285
!             -- ---- -- --------- ----- ------ --- ------- -------
!
      devkst = 0.23_wp
      devks = 32.8_wp * arafra + 35.4_wp * calfra
!     WRITE(io_stdo_bgc,*) '***DEVKS=',DEVKS,' DEVKST=',DEVKST,TEN
!
!     -----------------------------------------------------------------
!*        11. SET UNIVERSAL GAS CONSTANT
!             --- --------- --- --------
!
         rgas = 83.143_wp
!
!     -----------------------------------------------------------------
!*        12. SET BORON CONCENTRATION IN SEA WATER
!*            IN G/KG PER O/OO CL ACCORDING
!             TO RILEY AND SKIRROW, 1965 (P.250)
!             -- ----- --- -------- ---- ------- -
!
      bor1 = 0.00023_wp
!
!     -----------------------------------------------------------------
!*        13. SET INVERSE OF ATOMIC WEIGHT OF BORON [G**-1]
!             (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
!             ----- -- ------- -------- ----- ----- ---- ---------------
!
      bor2 = 1._wp/10.82_wp
!
!     -----------------------------------------------------------------
!*        14. SET INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS
!             [CM**3]
!             ---------- -- ------ ----- ------ -- -- ----- ---
!
      oxyco = 1._wp/22414.4_wp
!
!     -----------------------------------------------------------------
!*        15. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR CO2 IN ML/L
!             WEISS, R. F. (1974)
!             CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A
!             NON IDEAL GAS. MARINE CHEMISTRY, VOL. 2, 203-215.
!     -----------------------------------------------------------------

      c00 = -58.0931_wp          ! C null null
      c01 = 90.5069_wp
      c02 = 22.2940_wp
      c03 = 0.027766_wp
      c04 = -0.025888_wp
      c05 = 0.0050578_wp
!
!     -----------------------------------------------------------------
!*        16. SET COEFF. FOR 1. DISSOC. OF CARBONIC ACID
!             (EDMOND AND GIESKES, 1970)
!             ------- --- -------- ------- -------- ----
!
      c10 = 812.27_wp
      c11 = 3.356_wp
      c12 = -0.00171_wp
      c13 =  0.000091_wp

!     -----------------------------------------------------------------
!*        17. SET COEFF. FOR 2. DISSOC. OF CARBONIC ACID
!             (EDMOND AND GIESKES, 1970)
!             ------- --- -------- ------- -------- ----
!
      c20 = 1450.87_wp
      c21 = 4.604_wp
      c22 = -0.00385_wp
      c23 =  0.000182_wp
!
!     -----------------------------------------------------------------
!*        18. SET COEFF. FOR 1. DISSOC. OF BORIC ACID
!             (EDMOND AND GIESKES, 1970)
!             ------- --- -------- ------- ----- ----
!
      cb0 = 2291.90_wp
      cb1 = 0.01756_wp
      cb2 = -3.385_wp
      cb3 = -0.32051_wp
!
!     -----------------------------------------------------------------
!*        19. SET COEFF. FOR DISSOC. OF WATER
!             (DICKSON AND RILEY, 1979, EQ. 7, COEFFICIENT
!             CW2 CORRECTED FROM 0.9415 TO 0.09415 AFTER
!             PERS. COMMUN. TO B. BACASTOW, 1988)
!             ----- ------- -- -- --------- ------ -------
!
      cw0 = 3441._wp
      cw1 = 2.241_wp
      cw2 = -0.09415_wp
!
!     -----------------------------------------------------------------
!*        20. SET VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L
!             (WEISS, 1970)
!             ------- ------ --------- --------- --- -- -- ----
!
      ox0 = -173.4292_wp
      ox1 = 249.6339_wp
      ox2 = 143.3483_wp
      ox3 = -21.8492_wp
      ox4 = -0.033096_wp
      ox5 = 0.014259_wp
      ox6 = -0.0017_wp

!     -----------------------------------------------------------------
!*            SET VOLUMETRIC SOLUBILITY CONSTANTS FOR N2 IN ML/L
!             WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN
!             OXYGEN AND ARGON IN WATER AND SEAWATER.
!             DEEP-SEA RESEARCH, VOL. 17, 721-735.
!     -----------------------------------------------------------------

       an0 = -172.4965_wp
       an1 = 248.4262_wp
       an2 = 143.0738_wp
       an3 = -21.7120_wp
       an4 = -0.049781_wp
       an5 = 0.025018_wp
       an6 = -0.0034861_wp

!      Constants for laughing gas solubility
!      (WEISS, 1974, MARINE CHEMISTRY)
!      --------------------------------------
       a1 = -62.7062_wp
       a2 = 97.3066_wp
       a3 = 24.1406_wp
       b1 = -0.058420_wp
       b2 = 0.033193_wp
       b3 = -0.0051313_wp
       atn2o = 3.e-7_wp

      rrrcl = salchl * 1.025_wp * bor1 * bor2
!
!     -----------------------------------------------------------------
!*        21. CHEMICAL CONSTANTS - SURFACE LAYER
!             -------- --------- - ------- -----
!
!     call bounds_exch(1,'p+',ptho,'chemcon 1')
!     call bounds_exch(1,'p+',psao,'chemcon 2')
!     call bounds_exch(1,'p+',pddpo,'chemcon 3')


!$OMP PARALLEL PRIVATE(t,q,s,rs,cl,cek0,ckb,ck1,ck2,ckw,oxy,ani, &
!$OMP                  ak1,ak2,akb,akw,ak0,aksp0,cp,tc,p)
      IF(modal.LE.0) THEN

!$OMP DO
      DO 1 j=1,kpje
      DO 1 i=1,kpie
      IF (pddpo(i, j, 1) .GT. 0.5_wp) THEN           ! wet cell

!
!*        21.1 SET ABSOLUTE TEMPERATURE
!              ------------------------
      T = ptho(i,j,1) + tmelt                  ! degC to K
      Q=T*PERC                                 ! perc=0.01
      s = MAX(25._wp, psao(i, j, 1))           ! minimum salinity 25

!      Laughing gas solubility (WEISS, 1974)
!      --------------------------------------
      rs = a1 + a2 * (100._wp / t) + a3 * LOG(t/100._wp)                  &
     &    +s*( b1 +b2*(t/100._wp) + b3*(t/100._wp)**2)

       satn2o(i,j)=atn2o*exp(rs)

!
!*        21.2 CHLORINITY (WOOSTER ET AL., 1969)
!              ---------------------------------
!
      CL=S*SALCHL
!
!*        21.3 LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
!              -------------------------------------------------

      cek0 = c00 + c01 / Q + c02 * LOG(q) + s * (c03 + c04 * q + c05 * Q**2)

!
!*        21.4 PK1, PK2 OF CARB. ACID, PKB OF BORIC ACID
!              -----------------------------------------
!*             AFTER EDMOND AND GIESKES (1970)
!              -------------------------------

      CKB=CB0/T+CB1*T+CB2+CB3*CL**THIRD
      ck1 = c10 / t + c11 + c12 * s * LOG(t) + c13 * s**2
      ck2 = c20 / t + c21 + c22 * s * LOG(t) + c23 * s**2

!
!*        21.5 CKW (H2O) (DICKSON AND RILEY, 1979)
!              ------------------------------------

      CKW=CW0/T+CW1+CW2*SQRT(S)

!
!*****CKW COULD ADDITIONALLY BE EXPRESSED SALIN. DEPENDENT *********
!

!
!*        21.6 LN(K0) OF SOLUBILITY OF O2 (EQ. 4, WEISS, 1970)
!              -----------------------------------------------

      oxy = ox0 + ox1 / q + ox2 * LOG(q) + ox3 * q &
           + s * (ox4 + ox5 * q + ox6 * q**2)

!*       SOLUBILITY OF N2
!        WEISS, R. F. (1970), DEEP-SEA RESEARCH, VOL. 17, 721-735.
!              -----------------------------------------------

      ani = an0 + an1 / q + an2 * LOG(q) + an3 * q &
           + s * (an4 + an5 * q + an6 * q**2)

!
!*        21.7 K1, K2 OF CARB. ACID, KB OF BORIC ACID (EDMOND AND GIESKES,1970)
!              ----------------------------------------------------------------
      AK1=TEN**(-CK1)
      AK2=TEN**(-CK2)
      AKB=TEN**(-CKB)
!
!*        21.8 IONIC PRODUCT OF WATER KW (H2O) (DICKSON AND RILEY, 1979)
!              ----------------------------------------------------------------
      AKW=TEN**(-CKW)
      AKW3(I,J,1)=AKW
!
!*       21.9 CO2 SOLUBILITY IN SEAWATER (WEISS, 1974, CF. EQ. 12)
!              ----------------------------------------------------------------
      AK0=EXP(CEK0)*SMICR
!
!*       21.10 DENSITY OF SEAWATER AND TOTAL BORATE IN MOLES/L

!      RRR=RHO(S,ptho(i,j,1),ZERO) *THOUSI
!      BOR=BOR1*RRR*CL*BOR2

!      reformulation after R. Bacastow
!      BOR=1.22e-5*S

!
!*       21.11 SET CHEMICAL CONSTANTS
      CHEMCM(i,j,5,kplmonth)=AK0
      CHEMCM(i,j,4,kplmonth)=ak1
      CHEMCM(i,j,3,kplmonth)=ak2
      CHEMCM(i,j,1,kplmonth)=akb
      CHEMCM(i,j,2,kplmonth)=AKW
      CHEMCM(i,j,6,kplmonth)=1.22e-5_wp * S    !  = BOR
!
!*       21.12 O2/N2 SOLUBILITY IN SEAWATER (WEISS, 1970)
!              ----------------------------------------------------------------
      CHEMCM(i,j,7,kplmonth)=EXP(OXY)*OXYCO
      CHEMCM(i,j,8,kplmonth)=EXP(ANI)*OXYCO

!      CHEMC(I,J,7,LAUMON)=EXP(OXY)*OXYCO/196800.
!      CHEMC(I,J,8,LAUMON)=EXP(ani)*OXYCO/802000.

      ENDIF

1     CONTINUE

      ENDIF
!Paddy:

      IF(modal.EQ.-13) THEN          ! call from ini_bgc
!$OMP SINGLE
         WRITE(io_stdo_bgc,*) 'CHEMCON: initializing all month '
!$OMP END SINGLE
!$OMP DO
         DO kmon=1,12
         DO kchem=1,8
         DO j=1,kpje
         DO i=1,kpie
            CHEMCM(i,j,kchem,kmon) = CHEMCM(i,j,kchem,kplmonth)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF


!js   IF ( kchck .NE. 0 ) THEN
!js      DO kchem=1,8
!js         WRITE(io_stdo_bgc,*) 'chemcm:k=',kchem
!js         CALL EXTR(kpie,kpje,chemcm(1,1,kchem,kplmonth),pddpo(1,1,1),   &
!js  &          rmasko,io_stdo_bgc)
!js      ENDDO

!js    i = ie_g/2 - p_ioff
!js      IF(i>=1 .AND. i<=kpie) THEN
!js        WRITE(io_stdo_bgc,*)                                              &
!js  &       'Matrix with chemical constants (CHEMCM) at i=kpie/2; k=1->8:'
!js        WRITE(io_stdo_bgc,6003)                                           &
!js  &       (j,INT(ptiestu(1)),(chemcm(i,j,k,kplmonth),k=1,8),j=1,kpje)
!js      ENDIF
!js   ENDIF

! 6003 FORMAT(I3.0,'(',I5.0,'m):',7e12.5)


!
!     -----------------------------------------------------------------
!*        22. CHEMICAL CONSTANTS - DEEP OCEAN
!             ----------------------------------------------------------------

!$OMP DO
      DO 2 k=1,kpke
!
!*        22.1 APPROX. SEAWATER PRESSURE AT U-POINT DEPTH (BAR)
!              ----------------------------------------------------------------

      p = 1.025e-1_wp * ptiestu(k)
!      WRITE(io_stdo_bgc,*) 'CHEMCON: P=1.025E-1*ptiestu(k)', P,ptiestu(k),k

      DO 2 i=1,kpie
      DO 2 j=1,kpje

!
!*        22.1.1 Zonal mean temperature/salinity
!                -------------------------------
!
!      tzsum = 0.0
!      szsum = 0.0
!      vzsum = 0.0
!      DO i=1,kpie
!         IF(pddpo(i,j,k).GT.0.5) THEN
!            tzsum = tzsum+ptho(i,j,k)*pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)
!            szsum = szsum+psao(i,j,k)*pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)
!            vzsum = vzsum+            pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)
!         ENDIF
!      ENDDO
!      IF(vzsum.GT.0.5) THEN
!         tzmean = tzsum/vzsum
!         szmean = szsum/vzsum
!      ELSE
!         tzmean = 0.0
!         szmean = 0.0
!      ENDIF
!
!
!*        22.2 SET LIMITS FOR SEAWATER TEMP. AND SALINITY
!              ----------------------------------------------------------------
!              (THIS IS DONE TO AVOID COMPUTATIONAL CRASH AT DRY
!               POINTS DURING CALCULATION OF CHEMICAL CONSTANTS)
!
!*        22.3 SET [H+] (FIRST GUESS)
!              ----------------------------------------------------------------

!
!*        22.4 SET ABSOLUTE TEMPERATURE
!              ----------------------------------------------------------------
      t = ptho(i,j,k) + tmelt
      q=t*perc
      s = MAX(34._wp, psao(i, j, k))

!
!*        22.5 CHLORINITY (WOOSTER ET AL., 1969)
!              ----------------------------------------------------------------
      CL=S*SALCHL
!
!*        22.6 LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
!              ----------------------------------------------------------------
      cek0 = c00 + c01 / q + c02 * LOG(q) + s * (c03 + c04 * q + c05 * q**2)
!
!*        22.7 PK1, PK2 OF CARBONIC ACID, PKB OF BORIC ACID
!              ----------------------------------------------------------------
!              AFTER EDMOND AND GIESKES (1970)

      CKB=CB0/T+CB1*T+CB2+CB3*CL**THIRD
!
      ck1 = c10 / t + c11 + c12 * s * LOG(t) + c13 * s**2
      ck2 = c20 / t + c21 + c22 * s * LOG(t) + c23 * s**2

!*        22.8 LN(K0) OF SOLUBILITY OF O2 (EQ. 4, WEISS, 1970)
!              ----------------------------------------------------------------
      oxy = ox0 + ox1 / q + ox2 * LOG(q) + ox3 * q &
           + s * (ox4 + ox5 * q + ox6 * q**2)

      satoxy(i,j,k)=exp(oxy)*oxyco
!
!*        22.9 K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
      AK1=TEN**(-CK1)
      AK2=TEN**(-CK2)
      AKB=TEN**(-CKB)
!
!*       22.10 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
!              ----------------------------------------------------------------
!              (S=27-43, T=2-25 DEG C) AT P=0 (ATMOSPH. PRESSURE)
!              (INGLE, 1800, EQ. 6)

       aksp0 = 1.E-7_wp * (akcc1 + akcc2 * s**(THIRD) &
            + akcc3 * LOG10(s) + akcc4 * t**2)
       CKW=CW0/T+CW1+CW2*SQRT(S)
       AKW3(I,J,K)=TEN**(-CKW)
!
!*       22.11 FORMULA FOR CP AFTER EDMOND AND GIESKES (1970)
!              ----------------------------------------------------------------

!              (REFERENCE TO CULBERSON AND PYTKOQICZ (1968) AS MADE
!              IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
!              TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF P IN
!              DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CP IS
!              MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
!              WITH BASIS E IN THE FORMULA FOR AKSP (CF. EDMOND
!              AND GIESKES (1970), P. 1285 AND P. 1286 (THE SMALL
!              FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
!              SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON
!              P. 1285))
!      WRITE(io_stdo_bgc,*)  'CHEMCON: CP=P/(RGAS*T)', CP,P,RGAS,T
      CP=P/(RGAS*T)
!
!*       22.12 KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
!              CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
!              (CF. BROECKER ET AL., 1982)

      TC=ptho(i,j,k)

!      WRITE(io_stdo_bgc,*)  ' CHEMCON: modal&(CP*(DEVKB-DEVKBT*TC)) ',modal,CP,DEVKB,DEVKBT,TC

      AKB3(I,J,K)=AKB*EXP(CP*(DEVKB-DEVKBT*TC))
      AK13(I,J,K)=AK1*EXP(CP*(DEVK1-DEVK1T*TC))
      AK23(I,J,K)=AK2*EXP(CP*(DEVK2-DEVK2T*TC))
!
!        22.13 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE (OR ARAGONITE)
!              ----------------------------------------------------------------
!              AS FUNCTION OF PRESSURE FOLLOWING EDMOND AND GIESKES (1970)
!              (P. 1285) AND BERNER (1976)
      AKSP(I,J,K)=ARACAL*AKSP0*EXP(CP*(DEVKS-DEVKST*TC))



!
!*       22.14 DENSITY OF SEAWATER AND TOTAL BORATE CONCENTR. [MOLES/L]
!              ----------------------------------------------------------------

!      rrr=rho(s,tc,p)*thousi
!      bor=bor1*rrr*cl*bor2

!     reformulation after R. Bacastow
!     BOR=1.22e-5*S

!     rrrcl=salchl*1.025*bor1*bor2
! #slo# moved outside loops, s.a.



!
!     -----------------------------------------------------------------
!*        23. INITIATE [H+] AND [CO3--]
!             -------- ---- --- -------
!
!*       23.1  FIRST GUESSES FOR [CO3--] AND [H+]
!              ----------------------------------------------------------------

!      CARALK=OCETRA(i,j,k,IALKALI)-BOR/(ONE+TENM7/AKB3(I,J,K))

!
!*       20.15 DENSITY OF SEAWATER AND TOTAL BORATE IN MOLES/L
!              ----------------------------------------------------------------
!      RRR=RHO(S,TC,P)*THOUSI
!      BOR=BOR1*RRR*CL*BOR2

!     reformulation after R. Bacastow
!      BOR=1.22e-5*S
!
!*       20.16 [CO3--], FIRST ESTIMATE
!              ----------------------------------------------------------------

2     CONTINUE
!$OMP END PARALLEL

!     -----------------------------------------------------------------
!*        21. ITERATION TO INITIATE [H+] AND [CO3--]
!             --------- -- -------- ---- --- -------
!      BT=BOR
!      BREMS=1.5
!      DO 43 KI=1,30
!         zer(ki)=0.
!         DO 43 k=1,kpke
!         DO 43 j=1,kpje
!         DO 43 i=1,kpie
!
!            IF(pddpo(i,j,k).GT.0.5) THEN
!
!            ak1=ak13(i,j,k)
!            ak2=ak23(i,j,k)
!            H=HI(i,j,k)
!            R=CO3(i,j,k)
!            ALKA=OCETRA(i,j,k,IALKALI)
!            C=OCETRA(i,j,k,ISCO212)
!            T1=h/ak13(i,j,k)
!            T2=h/ak23(i,j,k)
!            AKW=AKW3(J,K)
!            BT=rrrcl*psao(i,j,k)
!            AKB=AKB3(J,K)
!            alk=ocetra(i,j,k,ialkali)
!            A=!*(2.+t2)/(1.+t2+t2*t1)  +AKW/H-H+BT/(1.+H/AKB)-ALK
!            A=!*(2.+t2)/(1.+t2+t2*t1)  +AKW/H-H+BT/(1.+H/AKB)-ALK
!            zer(ki)=zer(ki)+a**2
!            DADH=!*(1./(AK2*(1.+T2+T2*T1))-(2.+T2)*(1./AK2+2.*T1/AK2)/
!     1          (1.+T2+T2*T1)**2)
!     1          -AKW/H**2-1.-(BT/AKB)/(1.+H/AKB)**2
!            dddhhh=a/dadh
!            reduk=MAX(1.,1.2*abs(dddhhh/h))
!            H=H-dddhhh/reduk
!
!            HI(i,j,k)=HI(i,j,k)-dddhhh
!            co3(i,j,k)
!     1      =c/(1.+hi(i,j,k)*(1.+hi(i,j,k)/ak13(i,j,k))/ak23(i,j,k))
!
!           ENDIF
!
!43      CONTINUE

!      WRITE(io_stdo_bgc,*)  ' '
!      WRITE(io_stdo_bgc,*)  'CHEMCON: convergence ', zer

      RETURN
      END
