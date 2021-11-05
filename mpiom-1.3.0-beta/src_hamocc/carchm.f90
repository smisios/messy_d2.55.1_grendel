      SUBROUTINE CARCHM                                               &
     &           (kpie,kpje,kpke,pddpo,psao,ptho,psicomo,             &
     &            pfu10,kplmon,kplday,kmonlen,pdlxp,pdlyp)
!**********************************************************************
!
!**** *CARCHM* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - rename: ssso12(i,j,k)=sedlay(i,j,k,issso12 ) etc.; no equivalence statements
!     - rename: powasi(i,j,k )=powtra(i,j,1,ipowasi) etc.; no equivalence statements
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
!SL" the transfer velocity (piston velocity) is frequently expressed as
!SL  3 m/day or 12.5 cm/hour for typical moderate wind conditions
!SL
!SL  We set Vst = 3m/86400 s (to be modified in the coupled mode with
!SL  the actual scalar windspeed according to Wanninkhof).
!SL  Vst/Dzw(1) then is the time constant for the surface layer"
!SL
!SL  da ist dann klar, dasz noch mit dt multipliziert werden musz.
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
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st REAL) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd REAL) [m].
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

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_timeser_bgc
      USE mo_control_bgc
      USE mo_bgcmean
      use mo_param1_bgc 

      USE MO_COMMO1
      USE MO_COMMOAU1

      use mo_parallel

      implicit none
      INTEGER :: i,j,k,l,kpie,kpje,kpke, kk,ii,jj
      INTEGER :: kplmon,kplday,kmonlen,laumo1
      
      REAL psao(kpie,kpje,kpke)
      REAL pddpo(kpie,kpje,kpke)
      REAL psicomo(kpie,kpje)
      REAL pfu10(kpie,kpje)
      REAL ptho(kpie,kpje,kpke)
      REAL pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      
      REAL :: supsat, undsa, dissol
      REAL :: fluxd,fluxu
      REAL :: dddhhh,dadh,a,h,c,alk,t1,t2
      REAL :: akbi,ak2i,ak1i
      REAL :: kwco2,kwo2,kwdms
      REAL :: scco2,sco2,scdms,atmload12
      REAL :: contppm, Xconvxa
      REAL :: oxflux,niflux,dmsflux,n2oflux
      REAL :: ato2, atn2, atco2,pco2
      REAL :: AHI, ABE,RMONLEN,RPLDAY
      REAL :: AK0,AK1,AK2,AKB,AKW,BT,oxysa,anisa
#ifdef __c_isotopes
      REAL :: r13,r14,rat13,rat14
      REAL :: flux14d,flux14u,flux13d,flux13u
      REAL :: atc13,atc14
#endif
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

      REAL :: thickness

#ifdef ANTC14
      REAL :: fantc14d,fantc14u
#endif

!      WRITE(*,*) 'CARCHM called with :',                             &
!     &           kpie,kpje,kpke,pddpo(50,50,1),psao(50,50,1),        &
!     &           ptho(50,50,1),psicomo(50,50),pfu10(50,50),          &
!     &           kplyear,kplmon,kplday,kmonlen

!      call bounds_exch('p+',pddpo,'carchm 1')      
!      call bounds_exch('p+',psao,'carchm 2')      
!      call bounds_exch('p+',ptho,'carchm 3')      
!      call bounds_exch('p+',psicomo,'carchm 4')      
!      call bounds_exch('p+',pfu10,'carchm 5')      
!      call bounds_exch('p+',atm,'carchm 6')      
!      call bounds_exch('p+',zo,'carchm 7')      
!      call bounds_exch('p+',sictho,'carchm 7')      
!      call bounds_exch('p+',sicsno,'carchm 8')
  
!      call bounds_exch('p+',ddpo,'carchm 9')      
!      call bounds_exch('p+',pdlxp,'carchm 10')      
!      call bounds_exch('p+',pdlyp,'carchm 11')        

!      call bounds_exch('p+',akb3,'carchm 12')        
!      call bounds_exch('p+',akw3,'carchm 13')        
!      call bounds_exch('p+',hi,'carcahm 14')        
   
!      do kk=1,12
!      call bounds_exch('p+',chemcm(:,:,:,kk),'carchm 15')      
!      enddo
!      do kk=1,NOCETRA
!      call bounds_exch('p+',ocetra(:,:,:,kk),'carchm 16')      
!      enddo



      laumo1=kplmon+1
      IF(laumo1.GT.12) laumo1=1
        
      RMONLEN=kmonlen
      RPLDAY=kplday
      AHI=RPLDAY/RMONLEN
      ABE=1.-AHI
      atmload12=atm_co2

#ifdef __c_isotopes
      write(io_stdo_bgc,*)'eingang carchm',atm_co2,atm_c13,atm_c14 !e
#endif

#ifndef DIFFAT
       atco2 = atm_co2
       ato2  = atm_o2 
       atn2  = atm_n2 


#ifdef __c_isotopes
         atc13=  atm_c13
         atc14=  atm_c14
#endif /*__c_isotopes*/
#endif /*DIFFAT*/ 


      k=1      ! surface layer



      DO 1 j=1,kpje
!      DO 1 i=2,kpie
      DO 1 i=1,kpie

!#ifdef __cpl_co2
!     atm(i,j,iantco2)=co2conc(i,j)*(28.970/44.011)*1.e6
!#endif

      IF(pddpo(i,j,1).GT.0.5) THEN

!*       21.11 SET CHEMICAL CONSTANTS

!!$      AK0  =AHI*CHEMCM(i,j,5,LAUMO1)+ABE*CHEMCM(i,j,5,kplmon)
!!$      AK1  =AHI*CHEMCM(i,j,4,LAUMO1)+ABE*CHEMCM(i,j,4,kplmon)
!!$      AK2  =AHI*CHEMCM(i,j,3,LAUMO1)+ABE*CHEMCM(i,j,3,kplmon)
!!$      AKB  =AHI*CHEMCM(i,j,1,LAUMO1)+ABE*CHEMCM(i,j,1,kplmon)
!!$      AKW  =AHI*CHEMCM(i,j,2,LAUMO1)+ABE*CHEMCM(i,j,2,kplmon)
!!$      BT   =AHI*CHEMCM(i,j,6,LAUMO1)+ABE*CHEMCM(i,j,6,kplmon)
!!$      oxysa=AHI*CHEMCM(i,j,7,LAUMO1)+ABE*CHEMCM(i,j,7,kplmon)
!!$      anisa=AHI*CHEMCM(i,j,8,LAUMO1)+ABE*CHEMCM(i,j,8,kplmon)

      AK0  =CHEMCM(i,j,5,kplmon)
      AK1  =CHEMCM(i,j,4,kplmon)
      AK2  =CHEMCM(i,j,3,kplmon)
      AKB  =CHEMCM(i,j,1,kplmon)
      AKW  =CHEMCM(i,j,2,kplmon)
      BT   =CHEMCM(i,j,6,kplmon)
      oxysa=CHEMCM(i,j,7,kplmon)
      anisa=CHEMCM(i,j,8,kplmon)


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
!                                                          7373-7382)
!  DMS Schmidt number after Saltzmann et al. (1993, J. Geophys. Res. 98,
!                                                      16,481-16,486)
!  O2 Schmidt number after Keeling et al. (1998, Global Biogeochem.
!                                                Cycles, 12, 141-163)
!  CFC Schmidt number ref:  Zheng et al (1998), JGR, vol 103,No C1
!  ---> www.ipsl.jussieu.fr/OCMIP/
!*********************************************************************
!
       scco2 = 2073.1 - 125.62*ptho(i,j,1) + 3.6276*ptho(i,j,1)**2  &
     &           - 0.043219*ptho(i,j,1)**3

       scdms = 2674.0-147.12*ptho(i,j,1)+3.726*ptho(i,j,1)**2       &
     &       - 0.038*ptho(i,j,1)**3

       sco2  = 1638.0 - 81.83*ptho(i,j,1) + 1.483*ptho(i,j,1)**2    & 
     &       - 0.008004*ptho(i,j,1)**3  
#ifdef PCFC
       sccfc12 = 3845.4 - 228.95*ptho(i,j,1) + 6.1908*ptho(i,j,1)**2  &
     &         - 0.067430*ptho(i,j,1)**3
       sccfc11 = 3501.8 - 210.31*ptho(i,j,1) + 6.1851*ptho(i,j,1)**2  &
     &               - 0.07513*ptho(i,j,1)**3
#endif
!
!  Compute the transfer (piston) velocity in m/s
!
!     a : coefficient for a piston velocity in cm/hr 
!     Xconv : conversion to m/sec
!     Xconvxa = Xconv * a = 1/3.6e+05 * 0.337
      Xconvxa = 9.3611e-07      
      
!     660 = Schmidt number of CO2 @ 20 degC in seawater
         kwco2 = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2        &
     &         * (660/scco2)**0.5 

         kwdms = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2        & 
     &         * (660/scdms)**0.5 

         kwo2  = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2        & 
     &         * (660/sco2)**0.5 
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

! checking/work needed
#ifdef DIFFAT              
#ifndef __cpl_co2
       atco2 = atm(i,j,iatmco2)
#endif
       ato2  = atm(i,j,iatmo2)
       atn2  = atm(i,j,iatmn2)
#ifdef __c_isotopes
       atc13 = atm(i,j,iatmc13)
       atc14 = atm(i,j,iatmc14)
#endif
#endif


#ifdef PANTHROPOCO2
#ifndef __cpl_co2
#ifndef DIFFAT
! introduced key NCEP, 18.10.06
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

#endif
#endif
#endif /*PANTHROPOCO2*/



! Surface flux of oxygen

       oxflux=kwo2*dtbgc*(ocetra(i,j,1,ioxygen)                      &
     &                          -oxysa*(ato2/196800)) ! *ppao(i,j)/101300. ! sea level pressure normalization


         ocetra(i,j,1,ioxygen)=ocetra(i,j,1,ioxygen)                   &
     &                               -oxflux/pddpo(i,j,1)

! Surface flux of gaseous nitrogen (same piston velocity as for O2)
         
       niflux=kwo2*dtbgc*(ocetra(i,j,1,igasnit)                      &
     &                     -anisa*(atn2/802000)) ! *ppao(i,j)/101300.

         ocetra(i,j,1,igasnit)=ocetra(i,j,1,igasnit)                   &
     &                               -niflux/pddpo(i,j,1)

! Surface flux of laughing gas (same piston velocity as for O2 and N2)
         
       n2oflux=kwo2*dtbgc*(ocetra(i,j,1,ian2o)                       &
     &                      -satn2o(i,j)) ! *ppao(i,j)/101300.

         ocetra(i,j,1,ian2o)=ocetra(i,j,1,ian2o)                       &
     &                               -n2oflux/pddpo(i,j,1)

#ifdef DIFFAT              
         atm(i,j,iatmo2)=atm(i,j,iatmo2) + oxflux *contppm
         atm(i,j,iatmn2)=atm(i,j,iatmn2) + niflux *contppm
         atm(i,j,iatmn2)=atm(i,j,iatmn2) + n2oflux*contppm !closing mass balance
#endif       


! Surface flux of dms
        
         dmsflux = kwdms*dtbgc*ocetra(i,j,1,idms)  
         ocetra(i,j,1,idms)=ocetra(i,j,1,idms)-dmsflux/pddpo(i,j,1)
    
!
! Write output for bgcmean
!     
         bgct2d(i,j,jo2flux)  =bgct2d(i,j,jo2flux) +oxflux
         bgct2d(i,j,jn2flux)  =bgct2d(i,j,jn2flux) +niflux
         bgct2d(i,j,jn2oflux) =bgct2d(i,j,jn2oflux)+n2oflux
       
       bgcm2d(i,j,joxflux)=                             &
     &          bgcm2d(i,j,joxflux)+oxflux
       bgcm2d(i,j,jniflux)=                             &
     &          bgcm2d(i,j,jniflux)+niflux
       bgcm2d(i,j,jdmsflux)=                            &
     &              bgcm2d(i,j,jdmsflux) + dmsflux
       bgcm2d(i,j,jdms)=                                &
     &              bgcm2d(i,j,jdms) + ocetra(i,j,1,idms)
       
! back to CO2:
! Calculate new hi concentration  (four iterations)


         h=hi(i,j,k)
         c=ocetra(i,j,k,isco212)
         akw=akw3(i,j,k)                             ! IONIC PRODUCT OF WATER
         bt=rrrcl*psao(i,j,k)                        ! salinity
         akb=akb3(i,j,k)
         alk=ocetra(i,j,k,ialkali)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2

         dddhhh=a/dadh
         hi(i,j,k)=hi(i,j,k)-dddhhh
         h=hi(i,j,k)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2 
         dddhhh=a/dadh
         hi(i,j,k)=hi(i,j,k)-dddhhh
         h=hi(i,j,k)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2
         dddhhh=a/dadh
         hi(i,j,k)=hi(i,j,k)-dddhhh
         h=hi(i,j,k)


!
! Calculate pCO2 [ppmv] from total dissolved inorganic carbon (DIC: SCO212)
! the calculation also includes solubility
!
!            co212
         pco2=  c  /((1.+ak1*(1.+ak2/h)/h)*ak0)
       suppco2(i,j)=pco2                     !

         fluxd=atco2*kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.
         fluxu=pco2 *kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.

#ifdef __cpl_co2
!PW : coupled model, flux in kg/m^2/sec for exchange with ECHAM
         co2flux(i,j) = 44.011*(fluxu-fluxd)/dtbgc  ! used in mo_couple.f90 of C4MIP setup
#endif
       
!proxies d13C stuff
#ifdef __c_isotopes
       Roc13=ocetra(i,j,1,isco213)/(ocetra(i,j,1,isco212)+1e-25)
       Roc14=ocetra(i,j,1,isco214)/(ocetra(i,j,1,isco212)+1e-25)

#ifdef DIFFAT
         rat13=atm(i,j,iatmc13)/atm(i,j,iatmco2)
         rat14=atm(i,j,iatmc14)/atm(i,j,iatmco2)
#else
         rat13=atc13/atco2
         rat14=atc14/atco2
#endif

       flux13d=atc13*kwco2*dtbgc*ak0              ! *ppao(i,j)/101300. (atm to ocean)
         flux13u=pco2 *kwco2*dtbgc*ak0*Roc13*0.9935 ! *ppao(i,j)/101300. (ocean to atm,
                                                                        !twofold fractionation through evaporation)

       flux14d=atc14*kwco2*dtbgc*ak0              ! *ppao(i,j)/101300.
         flux14u=pco2 *kwco2*dtbgc*ak0*Roc14*0.987  ! *ppao(i,j)/101300. (twofold fractionation through evaporation)
#endif /*__c_isotopes*/

#ifdef ANTC14
       Roc14=ocetra(i,j,1,iantc14)/(ocetra(i,j,1,isco212)+1e-15)

       fantc14d=atco2*kwco2*dtbgc*ak0*Rbomb(i,j)  ! *ppao(i,j)/101300.
         fantc14u=pco2 *kwco2*dtbgc*ak0*Roc14       ! *ppao(i,j)/101300.

         ocetra(i,j,1,iantc14)=                                     &
     &   ocetra(i,j,1,iantc14)+(fantc14d-fantc14u)/pddpo(i,j,1)

       bgcm2d(i,j,jac14fx) = bgcm2d(i,j,jac14fx)     &
     &                                         + fantc14d-fantc14u
#endif

       
#ifdef DIFFAT                
#ifndef __cpl_co2
         atm(i,j,iatmco2)=atm(i,j,iatmco2)+(fluxu-fluxd)*contppm
#endif
#ifdef __c_isotopes
         atm(i,j,iatmc13)=atm(i,j,iatmc13)+(flux13u-flux13d)*contppm
         atm(i,j,iatmc14)=atm(i,j,iatmc14)+(flux14u-flux14d)*contppm
#endif
#else

!!$#ifdef __c_isotopes
!!$        atm_c13=atm_c13 + pdlxp(i,j)*pdlyp(i,j)              &
!!$     &         *(flux13u-flux13d)/atmacmol
!!$        atm_c14=atm_c14 + pdlxp(i,j)*pdlyp(i,j)              &
!!$     &         *(flux14u-flux14d)/atmacmol
!!$#endif /*__c_isotopes*/
#endif /*DIFFAT*/       

! new concentrations ocean (kmol/m3 -->ppm)
          thickness    = (DDPO(i,j,k)                               &    ! still k=1
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA               &
     &                   -SICSNO(I,J)*RHOSNWA)

!       ocetra(i,j,1,isco212)=                                     &
!    &   ocetra(i,j,1,isco212)+(fluxd-fluxu)/ddpo(i,j,k)

       ocetra(i,j,1,isco212)=                                     &
     &   ocetra(i,j,1,isco212)+(fluxd-fluxu)/thickness   ! test 19.10.06  conserves mass but check conflict with dilute_bgc/dilcor

#ifdef __c_isotopes
         ocetra(i,j,1,isco213)=                                     &
     &   ocetra(i,j,1,isco213)+(flux13d-flux13u)/thickness
         ocetra(i,j,1,isco214)=                                     &
     &   ocetra(i,j,1,isco214)+(flux14d-flux14u)/thickness
#endif
     
!
! write output for bgcmean      
!
        bgct2d(i,j,jco2flux)=bgct2d(i,j,jco2flux)+(fluxu-fluxd)
#ifdef __c_isotopes
        bgct2d(i,j,jc13flux)=bgct2d(i,j,jc13flux)+(flux13u-flux13d)
        bgct2d(i,j,jc14flux)=bgct2d(i,j,jc14flux)+(flux14u-flux14d)
#endif

        bgcm2d(i,j,jco2fxd)=bgcm2d(i,j,jco2fxd)    &
     &                                  +fluxd
        bgcm2d(i,j,jco2fxu)=bgcm2d(i,j,jco2fxu)    &
     &                                  +fluxu
        bgcm2d(i,j,jpco2)  =bgcm2d(i,j,jpco2)      &
     &                                  +pco2
      bgcm2d(i,j,jkwco2) =bgcm2d(i,j,jkwco2)     &
     &                                  +kwco2*ak0

 1000 CONTINUE 

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
     &                                     + flux_cfc11
        bgcm2d(i,j,jcfc12fx) = bgcm2d(i,j,jcfc12fx)    &
     &                                     + flux_cfc12
        bgcm2d(i,j,jpcfc11 ) = bgcm2d(i,j,jpcfc11)     &
     &                                     + pcfc11
        bgcm2d(i,j,jpcfc12 ) = bgcm2d(i,j,jpcfc12)     &
     &                                     + pcfc12
#endif /*CFC*/
  
      ENDIF ! wet cell
    1 CONTINUE    ! i,j loop



      call maschk(kpie,kpje,kpke,76)

!
!     -----------------------------------------------------------------
!*        22. CHEMICAL CONSTANTS - water column
   
      DO 43 k=1,kpke
      DO 43 j=1,kpje
      DO 43 i=1,kpie
         IF(pddpo(i,j,k).GT.0.5) THEN
            h=hi(i,j,k)
            c=ocetra(i,j,k,isco212)
            t1=h/ak13(i,j,k)
            t2=h/ak23(i,j,k)
            ak2=ak23(i,j,k)
            akw=akw3(i,j,k)
            bt=rrrcl*psao(i,j,k)
            akb=akb3(i,j,k)
            alk=ocetra(i,j,k,ialkali)
            a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h/akb)-alk
            dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*(1./ak2+2.*t1/ak2)/ &
     &          (1.+t2+t2*t1)**2)                                      &
     &          -akw/h**2-1.-(bt/akb)/(1.+h/akb)**2
            dddhhh=a/dadh
            h=h-dddhhh              
            hi(i,j,k)=hi(i,j,k)-dddhhh
            co3(i,j,k)=                                                &
     &      c / (1.+ hi(i,j,k)*(1.+ hi(i,j,k)/ak13(i,j,k)) / ak23(i,j,k))

! C14 decay, c14dec is set in BELEG_BGC , c14ret = 1-c14dec 
#ifdef __c_isotopes
            ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)*c14ret
            ocetra(i,j,k,idet14) =ocetra(i,j,k,idet14) *c14ret
            ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14)*c14ret
#endif
          ENDIF ! wet cell
   43  CONTINUE ! i,j,k loop

! C14 decay in the sediment
#ifdef __c_isotopes
        do k=1,ks
        do j=1,kpje
        do i=1,kpie
        if(bolay(i,j).gt.0.) then
        sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)*c14ret
        sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc14)*c14ret
        powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)*c14ret
        endif
        enddo
        enddo
        enddo
        atm_c14=atm_c14*c14ret+c14prod*c14ret
#endif

!
! Dissolution of calcium
! Note : mixed layer (k=1) is assumed to be always supersaturated
!        (saturation in reality depends on temperature/DIC/alkalinity)
!
        call maschk(kpie,kpje,kpke,77)
        DO 11 k=2,kpke
        DO 11 j=1,kpje
        DO 11 i=1,kpie
         IF(pddpo(i,j,k).GT.0.5) THEN
           supsat=co3(i,j,k)-97.*aksp(i,j,k)                       ! 97. = 1./1.03e-2 (MEAN TOTAL [CA++] IN SEAWATER [kmol/m3])
           undsa=MAX(0.,-supsat)
           dissol=MIN(undsa,dremcalc*ocetra(i,j,k,icalc))              
#ifdef __c_isotopes
           r13=dissol*ocetra(i,j,k,icalc13)                        &
      &             /(ocetra(i,j,k,icalc)+1.e-25)
           r14=dissol*ocetra(i,j,k,icalc14)                        &
      &             /(ocetra(i,j,k,icalc)+1.e-25)
#endif
           ocetra(i,j,k,icalc)=ocetra(i,j,k,icalc)-dissol
           ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)+2.*dissol

           ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+dissol
#ifdef __c_isotopes
           ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc13)-r13
           ocetra(i,j,k,isco213)=ocetra(i,j,k,isco213)+r13               ! remineralized calcite shells
           ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14)-r14
           ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)+r14
#endif
         ENDIF   ! wet cell
   11 CONTINUE   ! i,j,k loop
      call maschk(kpie,kpje,kpke,78)
!
! Sampling timeseries-1 : concentrations at specific positions
!
      DO l=1,nts               ! no of positions for which time series are produced
         i = its1(l)-p_ioff
         j = jts1(l)-p_joff
         k = 1
       
         IF(i<=1 .OR. i>=kpie .OR. j<=1 .OR. j>=kpje) CYCLE

!!$         AK0  =AHI*CHEMCM(i,j,5,LAUMO1)+ABE*CHEMCM(i,j,5,kplmon)
!!$         AK1  =AHI*CHEMCM(i,j,4,LAUMO1)+ABE*CHEMCM(i,j,4,kplmon)
!!$         AK2  =AHI*CHEMCM(i,j,3,LAUMO1)+ABE*CHEMCM(i,j,3,kplmon)


         AK0  =CHEMCM(i,j,5,kplmon)
         AK1  =CHEMCM(i,j,4,kplmon)
         AK2  =CHEMCM(i,j,3,kplmon)

       c=ocetra(i,j,k,isco212)
       h=hi(i,j,k)
       pco2=c/((1.+ak1*(1.+ak2/h)/h)*ak0)
       
       atco2 = atm_co2

#ifdef __cpl_co2
         atco2=co2conc(i,j)*(28.970/44.011)*1.e6
#endif

#ifdef DIFFAT
#ifndef __cpl_co2
         atco2=atm(i,j,iatmco2)
#endif
#endif

#ifndef DIFFAT
#ifndef __cpl_co2
! same switches as above, does not make sense - check (25.09.2006) should probably be ifndef DIFFAT changed 241006
!co2_atm_1, _2, _3 is the prescribed atmopheric pCO2 at beginning, middle, end of year

!!$      if (kplmon .le. 6) then
!!$         co2_start = co2_atm_1+(kplmon-1)*(co2_atm_2-co2_atm_1)/6
!!$         co2_end   = co2_atm_1+(kplmon  )*(co2_atm_2-co2_atm_1)/6
!!$         atco2     = AHI*co2_end+ABE*co2_start
!!$      else
!!$         co2_start = co2_atm_2+(kplmon-7)*(co2_atm_3-co2_atm_2)/6
!!$         co2_end   = co2_atm_2+(kplmon-6)*(co2_atm_3-co2_atm_2)/6
!!$         atco2     = AHI*co2_end+ABE*co2_start
!!$      endif

#endif
#endif

       
         fluxd=atco2*kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.
         fluxu=pco2 *kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.

       ts1(itsco2f,l+1,lts1) = ts1(itsco2f,l+1,lts1)                 &
     &                         + fluxu - fluxd
       ts1(itspco2,l+1,lts1) = ts1(itspco2,l+1,lts1)                 &
     &                         + pco2 
       ts1(itsatm,l+1,lts1)  = ts1(itsatm,l+1,lts1)                  & 
     &                         + atco2  

      ENDDO
!
! Check maximum/minimum values on wet/dry cells.
! Check for negative values of hi.
!
      IF( kchck .EQ. 1) THEN
         CALL CHCK_BGC(io_stdo_bgc,icyclibgc,                          &
     &       'Check values of ocean tracer at exit from SBR CARCHM :', &
     &       kpie,kpje,kpke,pddpo)

         DO  k=1,kpke
            DO  j=1,kpje
               DO  i=1,kpie
                  IF( hi(i,j,k) .LT. 0.0 ) THEN
                     WRITE(io_stdo_bgc,*)                              &
     &                   'CARCHM: invalid values of hi at i,j,k=',     &
     &                   i,j,k,hi(i,j,k)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END
