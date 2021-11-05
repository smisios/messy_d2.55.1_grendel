      SUBROUTINE CARCHM_ANT                                           &
     &           (kpie,kpje,kpke,pddpo,psao,ptho,psicomo,             &
     &            pfu10,kplmon,kplday,kmonlen)
!**********************************************************************
!
!**** *CARCHM* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - rename: ssso12(i,j,k)=sedlay(i,j,k,issso12 ) etc.; no equivalences
!     - rename: powasi(i,j,k )=powtra(i,j,1,ipowasi) etc.; no equivalences
!     - interfacing with ocean model
!
!     Purpose
!     -------
!     Inorganic carbon cycle.
!
!     Method
!     -------
!     Surface fluxes of CO2 / N2O / dms
!     Dissolution of calcium
!     Note: O2 solubility in seawater is calculated in chemin
!     kchck=1 can be used to check max/min of bgc arrays on wet/dry cells.
!
!     To do :
!     -------
!
!     Surface fluxes should not go through sea ice !
!
!
!     *CALL* *CARCHM(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,psao,ptho,psicomo,
!                                 pfu10,kplyear,kplmon,kplday,kmonlen)*
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL of model grid.
!     *INTEGER* *kpje*    - 2nd REAL of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd REAL) [m].
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *ptho*    - potential temperature (3rd REAL).
!     *REAL*    *ppao*    - sea level pressure [Pascal].
!     *REAL*    *psicomo* - sea ice (2nd REAL).
!     *REAL*    *pfu10*   - forcing field wind speed (2nd REAL).
!
!     Externals
!     ---------
!     .
!
!**********************************************************************

#ifdef PANTHROPOCO2
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_timeser_bgc
      USE mo_control_bgc
      USE mo_bgcmean
      use mo_param1_bgc 
      use mo_parallel

      implicit none
      INTEGER :: i,j,k,l,kpie,kpje,kpke
      INTEGER :: kplmon,kplday,kmonlen,laumo1
      
      REAL psao(kpie,kpje,kpke)
      REAL pddpo(kpie,kpje,kpke)
      REAL psicomo(kpie,kpje)
      REAL pfu10(kpie,kpje)
      REAL ptho(kpie,kpje,kpke)

      REAL :: supsat, undsa, dissol
      REAL :: fluxd,fluxu
      REAL :: dddhhh,dadh,a,h,c,alk,t1,t2
      REAL :: akbi,ak2i,ak1i
      REAL :: kwco2
      REAL :: scco2
      REAL :: contppm, Xconvxa
      REAL :: atco2,pco2
      REAL :: AHI, ABE,RMONLEN,RPLDAY
      REAL :: AK0,AK1,AK2,AKB,AKW,BT,oxysa,anisa
      REAL :: co2_start, co2_end
!hh#ifdef PCFC
      REAL :: kwcfc11, kwcfc12, sccfc11, sccfc12
      REAL :: sol_cfc11, sol_cfc12
      REAL :: flux_cfc11,flux_cfc12
      REAL :: ta,d 
      REAL :: cfc11_atm_1,cfc11_atm_2,cfc11_atm_3
      REAL :: cfc12_atm_1,cfc12_atm_2,cfc12_atm_3
      REAL :: pcfc11,cfc11_start,cfc11_end
      REAL :: pcfc12,cfc12_start,cfc12_end
!hh#endif


      laumo1=kplmon+1
      IF(laumo1.GT.12) laumo1=1
      
      RMONLEN=kmonlen
      RPLDAY=kplday
      AHI=RPLDAY/RMONLEN
      ABE=1.-AHI

!
!     -----------------------------------------------------------------
!*         1. SET HALF PRECISION CONSTANTS
!             --- ---- --------- ---------
!
!
      k=1

      DO 1 j=1,kpje
      DO 1 i=1,kpie
!PW : coupled model, only diag:
#ifdef __cpl_co2
     atm(i,j,iantco2)=co2conc(i,j)*(28.970/44.011)*1.e6
#endif

      IF(pddpo(i,j,1).GT.0.5) THEN   ! js: wet cell (this might not be useful for antco2)

!*       21.11 SET CHEMICAL CONSTANTS
      AK0  =AHI*CHEMCM(i,j,5,LAUMO1)+ABE*CHEMCM(i,j,5,kplmon)
      AK1  =AHI*CHEMCM(i,j,4,LAUMO1)+ABE*CHEMCM(i,j,4,kplmon)
      AK2  =AHI*CHEMCM(i,j,3,LAUMO1)+ABE*CHEMCM(i,j,3,kplmon)
      AKB  =AHI*CHEMCM(i,j,1,LAUMO1)+ABE*CHEMCM(i,j,1,kplmon)
      AKW  =AHI*CHEMCM(i,j,2,LAUMO1)+ABE*CHEMCM(i,j,2,kplmon)
      BT   =AHI*CHEMCM(i,j,6,LAUMO1)+ABE*CHEMCM(i,j,6,kplmon)
      oxysa=AHI*CHEMCM(i,j,7,LAUMO1)+ABE*CHEMCM(i,j,7,kplmon)
      anisa=AHI*CHEMCM(i,j,8,LAUMO1)+ABE*CHEMCM(i,j,8,kplmon)

      ak1i=1./ak1
      ak2i=1./ak2
      akbi=1./akb
      
      contppm=1./0.35e-3

!
!*********************************************************************
!
!  Compute the Schmidt number of CO2 in seawater and the transfer 
!  (piston) velocity using the formulation presented
!  by Wanninkhof (1992, J. Geophys. Res., 97, 7373-7382).
!  Input is temperature in deg C.
!  CO2 Schmidt number after Wanninkhof (1992, J. Geophys. Res., 97,
!  CFC Schmidt number ref:  Zheng et al (1998), JGR, vol 103,No C1 
!  ---> www.ipsl.jussieu.fr/OCMIP/
!*********************************************************************
!
       scco2   = 2073.1 - 125.62*ptho(i,j,1) + 3.6276*ptho(i,j,1)**2  &
     &             - 0.043219*ptho(i,j,1)**3

#ifdef PCFC
       sccfc12 = 3845.4 - 228.95*ptho(i,j,1) + 6.1908*ptho(i,j,1)**2  &  
     &         - 0.067430*ptho(i,j,1)**3

       sccfc11 = 3501.8 - 210.31*ptho(i,j,1) + 6.1851*ptho(i,j,1)**2  &
     &             - 0.07513*ptho(i,j,1)**3
#endif
!
!  Compute the transfer (piston) velocity in m/s
!
!     a : coefficient for a piston velocity in cm/hr 
!     Xconv : conversion to m/sec
!     Xconvxa = Xconv * a = 1/3.6e+05 * 0.337
      Xconvxa = 9.3611e-07      
      
         kwco2   = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2        &
     &         * (660/scco2)**0.5 

#ifdef PCFC
         kwcfc12 = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2        &
     &         * (660/sccfc12)**0.5 

         kwcfc11 = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2        &
     &         * (660/sccfc11)**0.5 
#endif    
 
!js: from cpl_co2
!PW : coupled model kg/kg --> ppm:
! molecular weight dry air : 28.970
! molecular weight CO2     : 44.011

#ifdef __cpl_co2 
     atco2=co2conc(i,j)*(28.970/44.011)*1.e6
#endif

#ifdef DIFFAT 
#ifndef __cpl_co2
      atco2=atm(i,j,iantco2)
#endif
#endif

#ifdef DIFFAT
#ifndef __cpl_co2 
! same switches as above, does not make sense - check (25.09.2006)
!co2_atm_1, _2, _3 is the prescribed atmopheric pCO2 at beginning, middle, end of year
      if (kplmon .le. 6) then
         co2_start = co2_atm_1+(kplmon-1)*(co2_atm_2-co2_atm_1)/6
       co2_end   = co2_atm_1+(kplmon  )*(co2_atm_2-co2_atm_1)/6
       atco2     = AHI*co2_end+ABE*co2_start
      else
         co2_start = co2_atm_2+(kplmon-7)*(co2_atm_3-co2_atm_2)/6
       co2_end   = co2_atm_2+(kplmon-6)*(co2_atm_3-co2_atm_2)/6
       atco2     = AHI*co2_end+ABE*co2_start     
      endif

      bgcm2d(i,j,jantco2)=bgcm2d(i,j,jantco2)    &    ! js: why in ifdef??
     &                                     +atco2      
#endif
#endif
       
! Calculate new hi_ant concentration

         h=hi_ant(i,j,k)
         c=ocetra(i,j,k,isco2_ant)
         akw=akw3(i,j,k)
         bt=rrrcl*psao(i,j,k)
         akb=akb3(i,j,k)
         alk=ocetra(i,j,k,ialk_ant)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2

         dddhhh=a/dadh
         hi_ant(i,j,k)=hi_ant(i,j,k)-dddhhh
         h=hi_ant(i,j,k)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2 
         dddhhh=a/dadh
         hi_ant(i,j,k)=hi_ant(i,j,k)-dddhhh
         h=hi_ant(i,j,k)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2
         dddhhh=a/dadh
         hi_ant(i,j,k)=hi_ant(i,j,k)-dddhhh
         h=hi_ant(i,j,k)

!
!
! Calculate pCO2 [ppmv] from total dissolved inorganic carbon (DIC: SCO212)
! the calculation also includes solubility
!
         pco2=c/((1.+ak1*(1.+ak2/h)/h)*ak0)

         fluxd=atco2*kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.
         fluxu=pco2 *kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.       

#ifdef __cpl_co2
!PW : coupled model, flux in kg/m^2/sec for exchange with ECHAM
          co2flux(i,j) = 44.011*(fluxu-fluxd)/dtbgc  ! used in mo_couple.f90 of C4MIP setup
#endif

#ifdef DIFFAT 
#ifndef __cpl_co2
         atm(i,j,iantco2)=atm(i,j,iantco2)+(fluxu-fluxd)*contppm
#endif
#endif
       ocetra(i,j,1,isco2_ant)=                                     &
     &   ocetra(i,j,1,isco2_ant)+(fluxd-fluxu)/pddpo(i,j,1)

!
! write output for bgcmean      
!
        bgct2d(i,j,jco2fant)  =  bgct2d(i,j,jco2fant)+(fluxu-fluxd)

        bgcm2d(i,j,jco2antdn)=bgcm2d(i,j,jco2antdn)    &
     &                                  +fluxd
        bgcm2d(i,j,jco2antup)=bgcm2d(i,j,jco2antup)    &
     &                                  +fluxu
        bgcm2d(i,j,jpco2_ant)=bgcm2d(i,j,jpco2_ant)    &
     &                                  +pco2
     

#ifdef PCFC

!     CFC 11 and 12 Solubilities in seawater
!     ref: Warner & Weiss (1985) , Deep Sea Research, vol32
!     coefficient for solubility in  mol/l/atm
!  ----------------------------------------
!     konstants given for mol/(l * atm) or kmol/(m^3 * atm)
!
!     for CFC 11
!     ----------

      ta = ( ptho(i,j,1) + 273.16)* 0.01
      
      d  = ( -0.0157274 * ta + 0.091459)* ta - 0.142382
      sol_cfc11 = exp ( - 229.9261                        &
     &                  + 319.6552 / ta                   &
     &                  + 119.4471 * alog ( ta )          &
     &                  - 1.39165  * ta * ta  + psao(i,j,1)* d )
!    
!     for CFC/12
!     ----------

      d    = ( -0.0153924 * ta + 0.091015)* ta - 0.143566
      sol_cfc12 = exp ( - 218.0971                        &
     &                  + 298.9702 / ta                   &
     &                  + 113.8049 * alog ( ta )          &
     &                  - 1.39165  * ta * ta  + psao(i,j,1)* d )

!
!     conversion from kmol/(m^3 * atm) to kmol/(m3 * pptv) 
!     --------------------------------------------------
      sol_cfc11 = 1.0e-12 * sol_cfc11
      sol_cfc12 = 1.0e-12 * sol_cfc12


!       F = Kw (Csat - Csurf)
!       Csat = alpha * pCFC * P/Po

      cfc11_atm_1 = cfc11_atm_1s*(1.-cfc_int(i,j))+cfc11_atm_1n*cfc_int(i,j)
      cfc11_atm_2 = cfc11_atm_2s*(1.-cfc_int(i,j))+cfc11_atm_2n*cfc_int(i,j)
      cfc11_atm_3 = cfc11_atm_3s*(1.-cfc_int(i,j))+cfc11_atm_3n*cfc_int(i,j)

      cfc12_atm_1 = cfc12_atm_1s*(1.-cfc_int(i,j))+cfc12_atm_1n*cfc_int(i,j)
      cfc12_atm_2 = cfc12_atm_2s*(1.-cfc_int(i,j))+cfc12_atm_2n*cfc_int(i,j)
      cfc12_atm_3 = cfc12_atm_3s*(1.-cfc_int(i,j))+cfc12_atm_3n*cfc_int(i,j)
      
      if (kplmon .le. 6) then
         cfc11_start = cfc11_atm_1+(kplmon+5)*(cfc11_atm_2-cfc11_atm_1)/12
       cfc11_end   = cfc11_atm_1+(kplmon+6)*(cfc11_atm_2-cfc11_atm_1)/12
       pcfc11      = AHI*cfc11_end + ABE*cfc11_start
      else
         cfc11_start = cfc11_atm_2+(kplmon-7)*(cfc11_atm_3-cfc11_atm_2)/12
       cfc11_end   = cfc11_atm_2+(kplmon-8)*(cfc11_atm_3-cfc11_atm_2)/12
       pcfc11      = AHI*cfc11_end + ABE*cfc11_start
      endif
      
      if (kplmon .le. 6) then
         cfc12_start = cfc12_atm_1+(kplmon+5)*(cfc12_atm_2-cfc12_atm_1)/12
       cfc12_end   = cfc12_atm_1+(kplmon+6)*(cfc12_atm_2-cfc12_atm_1)/12
       pcfc12      = AHI*cfc12_end + ABE*cfc12_start
      else
         cfc12_start = cfc12_atm_2+(kplmon-7)*(cfc12_atm_3-cfc12_atm_2)/12
       cfc12_end   = cfc12_atm_2+(kplmon-8)*(cfc12_atm_3-cfc12_atm_2)/12
       pcfc12      = AHI*cfc12_end + ABE*cfc12_start
      endif
            
       flux_cfc11 = kwcfc11*dtbgc *(sol_cfc11*pcfc11 - ocetra(i,j,1,icfc11))
       flux_cfc12 = kwcfc12*dtbgc *(sol_cfc12*pcfc12 - ocetra(i,j,1,icfc12))

       ocetra(i,j,1,icfc11) = ocetra(i,j,1,icfc11) + flux_cfc11/pddpo(i,j,1)
       ocetra(i,j,1,icfc12) = ocetra(i,j,1,icfc12) + flux_cfc12/pddpo(i,j,1)


!
! write output for bgcmean      
!
        bgcm2d(i,j,jcfc11fx) = bgcm2d(i,j,jcfc11fx)    &
     &                                         + flux_cfc11
        bgcm2d(i,j,jcfc12fx) = bgcm2d(i,j,jcfc12fx)    &
     &                                         + flux_cfc12
        bgcm2d(i,j,jpcfc11 ) = bgcm2d(i,j,jpcfc11)     &
     &                                         + pcfc11
        bgcm2d(i,j,jpcfc12 ) = bgcm2d(i,j,jpcfc12)     &
     &                                         + pcfc12
#endif

      ENDIF      ! wet cell
1     CONTINUE

!
!     -----------------------------------------------------------------
!*        22. CHEMICAL CONSTANTS - DEEP OCEAN

      DO 43 k=1,kpke
      DO 43 j=1,kpje
      DO 43 i=1,kpie
         IF(pddpo(i,j,k).GT.0.5) THEN
            h=hi_ant(i,j,k)
            c=ocetra(i,j,k,isco2_ant)
            t1=h/ak13(i,j,k)
            t2=h/ak23(i,j,k)
            ak2=ak23(i,j,k)
            akw=akw3(i,j,k)
            bt=rrrcl*psao(i,j,k)
            akb=akb3(i,j,k)
            alk=ocetra(i,j,k,ialk_ant)
            a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h/akb)-alk
            dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*(1./ak2+2.*t1/ak2)/ &
     &          (1.+t2+t2*t1)**2)                                      &
     &          -akw/h**2-1.-(bt/akb)/(1.+h/akb)**2
            dddhhh=a/dadh
            h=h-dddhhh              
            hi_ant(i,j,k)=hi_ant(i,j,k)-dddhhh
            co3_ant(i,j,k)=                                                &
     &      c/(1.+hi_ant(i,j,k)*(1.+hi_ant(i,j,k)/ak13(i,j,k))/ak23(i,j,k))
     
         ENDIF
43      CONTINUE

!
! Dissolution of calcium
! Note : mixed layer (k=1) is assumed to be always supersaturated
!        (in reality saturation depends on temperature/DIC/alkalinity)
!
        DO 11 k=2,kpke
        DO 11 j=1,kpje
        DO 11 i=1,kpie
         IF(pddpo(i,j,k).GT.0.5) THEN
           supsat=co3_ant(i,j,k)-97.*aksp(i,j,k)   ! -97 =1./calcon=1.03e-2 
           undsa=MAX(0.,-supsat)
           dissol=MIN(undsa,dremcalc*ocetra(i,j,k,icalc_ant))
           ocetra(i,j,k,icalc_ant)=ocetra(i,j,k,icalc_ant)-dissol
           ocetra(i,j,k,ialk_ant)=ocetra(i,j,k,ialk_ant)+2.*dissol
           ocetra(i,j,k,isco2_ant)=ocetra(i,j,k,isco2_ant)+dissol
         ENDIF
11    CONTINUE

!
! Sampling timeseries-1 : concentrations at specific positions
!
      DO l=1,nts
         i = its1(l)-p_ioff
         j = jts1(l)-p_joff
         k = 1
       
       
         IF(i<=1 .OR. i>=kpie .OR. j<=1 .OR. j>=kpje) CYCLE
         AK0  =AHI*CHEMCM(i,j,5,LAUMO1)+ABE*CHEMCM(i,j,5,kplmon)
         AK1  =AHI*CHEMCM(i,j,4,LAUMO1)+ABE*CHEMCM(i,j,4,kplmon)
         AK2  =AHI*CHEMCM(i,j,3,LAUMO1)+ABE*CHEMCM(i,j,3,kplmon)
         
       h=hi_ant(i,j,k)
         c=ocetra(i,j,k,isco2_ant)
       
       pco2=c/((1.+ak1*(1.+ak2/h)/h)*ak0)
       
#ifdef DIFFAT              
       atco2 = atm(i,j,iantco2)
#endif              
       
         fluxd=atco2*kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.
         fluxu=pco2 *kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.

       ts1(itsco2fa,l+1,lts1) = ts1(itsco2fa,l+1,lts1)               &
     &                          + fluxu - fluxd
       ts1(itspco2a,l+1,lts1) = ts1(itspco2a,l+1,lts1)               &
     &                          + pco2 
       ts1(itsatma,l+1,lts1)  = ts1(itsatma,l+1,lts1)                & 
     &                          + atco2  

      ENDDO

!
! Check maximum/minimum values on wet/dry cells.
! Check for negative values of hi_ant.
!
      IF( kchck .EQ. 1) THEN
         CALL CHCK_BGC(io_stdo_bgc,icyclibgc,                          &
     &       'Check values of ocean tracer at exit from CARCHM_ANT :', &
     &       kpie,kpje,kpke,pddpo)

         DO  k=1,kpke
            DO  j=1,kpje
               DO  i=1,kpie
                  IF( hi_ant(i,j,k) .LT. 0.0 ) THEN
                     WRITE(io_stdo_bgc,*)                              &
     &                   'CARCHM_ANT: invalid values of hi_ant at i,j,k=',     &
     &                   i,j,k,hi_ant(i,j,k)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

#endif /*PANTHROPOCO2*/
      RETURN
      END
