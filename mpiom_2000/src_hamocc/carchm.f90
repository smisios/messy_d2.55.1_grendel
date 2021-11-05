      SUBROUTINE CARCHM                                               &
     &           (kpie,kpje,kpke,pddpo,psao,ptho,psicomo,             &
     &            pfu10,kplmon,kplday,kmonlen,pdlxp,pdlyp,ptiestu)
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
!                                 pfu10,kplyear,kplmon,kplday,kmonlen,ptiestu)*
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
!     *REAL*    *ptiestu* - depth of scalar grid cell [m].
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

      USE mo_bgc_diagnostic              !CDI tinka

      USE MO_COMMO1
      USE mo_commoau1, ONLY: tmelt
      USE mo_planetary_constants, ONLY : rhoicwa, rhosnwa

      use mo_parallel
#ifdef AVFLUX
      use mo_avflux
#endif

      implicit none
      INTEGER :: i,j,k,l,kpie,kpje,kpke, iter!, kk
      INTEGER :: kplmon,kplday,kmonlen,laumo1

      REAL(wp) psao(kpie,kpje,kpke)
      REAL(wp) pddpo(kpie,kpje,kpke)
      REAL(wp) psicomo(kpie,kpje)
      REAL(wp) pfu10(kpie,kpje)
      REAL(wp) ptho(kpie,kpje,kpke)
      REAL(wp) pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL(wp),intent(in) :: ptiestu(kpke+1)
      INTEGER iflag(kpie,kpje)

      REAL(wp) :: supsat, undsa, dissol
      REAL(wp) :: fluxd,fluxu
      REAL(wp) :: dddhhh,dadh,a,h,c,alk,t1,t2
      REAL(wp) :: akbi,ak2i,ak1i
      REAL(wp) :: kwco2,kwo2,kwdms

      REAL(wp) c_int(kpie,kpje)  !help fields needed for integrating carbon flux from river input (millennium)
      REAL(wp) s_int(kpie,kpje)  !help fields needed for integrating carbon flux from river input (millennium)
      REAL(wp) o_int(kpie,kpje)  !help fields needed for integrating carbon flux from river input (millennium)
      REAL(wp) gs_c,gs_o,gs_s    !help fields needed for integrating carbon flux from river input (millennium)

#ifdef __cpl_co2
      REAL(wp) :: kwco2_cpl
#endif
      REAL(wp) :: scco2,sco2,scdms
      REAL(wp) :: Xconvxa
      REAL(wp) :: oxflux,niflux,dmsflux,nlaughflux
      REAL(wp) :: ato2, atn2, atco2,pco2
      REAL(wp) :: AHI, ABE,RMONLEN,RPLDAY
      REAL(wp) :: AK0,AK1,AK2,AKB,AKW,BT,oxysa,anisa
#ifdef __c_isotopes
      REAL(wp) :: r13,r14,rat13,rat14
      REAL(wp) :: flux14d,flux14u,flux13d,flux13u
      REAL(wp) :: atc13,atc14
#endif
#ifdef PCFC
      REAL(wp) :: kwcfc11, kwcfc12, sccfc11, sccfc12
      REAL(wp) :: sol_cfc11, sol_cfc12
      REAL(wp) :: flux_cfc11,flux_cfc12
      REAL(wp) :: ta,d
      REAL(wp) :: cfc11_atm_1,cfc11_atm_2,cfc11_atm_3
      REAL(wp) :: cfc12_atm_1,cfc12_atm_2,cfc12_atm_3
      REAL(wp) :: pcfc11,cfc11_start,cfc11_end
      REAL(wp) :: pcfc12,cfc12_start,cfc12_end
#endif

      REAL(wp) :: thickness
      REAL(wp) :: supsatup,satdiff,depthdiff   ! needed to calculate depth of lysokline

#ifdef ANTC14
      REAL(wp) :: fantc14d,fantc14u
#endif

!      WRITE(*,*) 'CARCHM called with :',                             &
!     &           kpie,kpje,kpke,pddpo(50,50,1),psao(50,50,1),        &
!     &           ptho(50,50,1),psicomo(50,50),pfu10(50,50),          &
!     &           kplyear,kplmon,kplday,kmonlen

!      call bounds_exch(1,'p+',pddpo,'carchm 1')
!      call bounds_exch(1,'p+',psao,'carchm 2')
!      call bounds_exch(1,'p+',ptho,'carchm 3')
!      call bounds_exch(1,'p+',psicomo,'carchm 4')
!      call bounds_exch(1,'p+',pfu10,'carchm 5')
!      call bounds_exch(1,'p+',atm,'carchm 6')
!      call bounds_exch(1,'p+',zo,'carchm 7')
!      call bounds_exch(1,'p+',sictho,'carchm 7')
!      call bounds_exch(1,'p+',sicsno,'carchm 8')

!      call bounds_exch(1,'p+',ddpo,'carchm 9')
!      call bounds_exch(1,'p+',pdlxp,'carchm 10')
!      call bounds_exch(1,'p+',pdlyp,'carchm 11')

!      call bounds_exch(1,'p+',akb3,'carchm 12')
!      call bounds_exch(1,'p+',akw3,'carchm 13')
!      call bounds_exch(1,'p+',hi,'carcahm 14')

!      do kk=1,12
!      call bounds_exch(1,'p+',chemcm(:,:,:,kk),'carchm 15')
!      enddo
!      do kk=1,NOCETRA
!      call bounds_exch(1,'p+',ocetra(:,:,:,kk),'carchm 16')
!      enddo



      laumo1=kplmon+1
      IF(laumo1.GT.12) laumo1=1

      rmonlen = REAL(kmonlen, wp)
      rplday = REAL(kplday, wp)
      AHI=RPLDAY/RMONLEN
      abe = 1._wp - ahi


      k=1      ! surface layer

#ifdef __cpl_co2
#ifdef AVFLUX
      co2flux_act = 0.0_wp
#endif
#endif

      DO 1 j=1,kpje
      DO 1 i=1,kpie



      IF (pddpo(i, j, 1) .GT. 0.5_wp) THEN

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


      ak1i = 1._wp / ak1
      ak2i = 1._wp / ak2
      akbi = 1._wp / akb

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
!
!  NB: most values (except for O2) were replaced by a new fit provided by
!      Matthias Groeger, MPI-M (personal communication 14.04.2010).
!      This fit is valid over a larger temperature range.
!      For reference, the previous values are provided as a comment.
!
!*********************************************************************

!        scco2 = 2073.1 - 125.62*ptho(i,j,1) + 3.6276*ptho(i,j,1)**2  &
!      &       - 0.043219*ptho(i,j,1)**3
       scco2 = 142.653_wp + 1955.9_wp *  exp(-0.0663147_wp*ptho(i,j,1))

!        scdms = 2674.0-147.12*ptho(i,j,1)+3.726*ptho(i,j,1)**2       &
!      &       - 0.038*ptho(i,j,1)**3
       scdms = 186.560_wp + 2506.78_wp * exp(-0.0618603_wp*ptho(i,j,1))

       sco2  = 1638.0_wp - 81.83_wp*ptho(i,j,1) + 1.483_wp*ptho(i,j,1)**2    &
     &       - 0.008004_wp*ptho(i,j,1)**3
#ifdef PCFC
!        sccfc12 = 3845.4_wp - 228.95_wp*ptho(i,j,1) + 6.1908_wp*ptho(i,j,1)**2  &
!      &         - 0.067430_wp*ptho(i,j,1)**3
       sccfc12 = 16.491_wp + 3631.41_wp * exp(-0.0651186_wp*ptho(i,j,1))

!        sccfc11 = 3501.8_wp - 210.31_wp*ptho(i,j,1) + 6.1851_wp*ptho(i,j,1)**2  &
!      &               - 0.07513_wp*ptho(i,j,1)**3
       sccfc11 = 316.571_wp + 3176.76_wp * exp(-0.066488_wp*ptho(i,j,1))
#endif
!
!  Compute the transfer (piston) velocity in m/s
!
!     a : coefficient for a piston velocity in cm/hr
!     Xconv : conversion to m/sec
!     Xconvxa = Xconv * a = 1/3.6e+05_wp * 0.337_wp
      Xconvxa = 9.3611e-07_wp

!     660 = Schmidt number of CO2 @ 20 degC in seawater
#ifdef __cpl_co2
      !FIXME: this should probably use SQRT
         kwco2_cpl = (1._wp - psicomo(i,j)) * Xconvxa                    &
     &           * (660._wp/scco2)**0.5_wp
#endif
         kwco2 = (1._wp - psicomo(i, j)) * Xconvxa * pfu10(i, j)**2        &
     &           * (660._wp / scco2)**0.5_wp

         kwdms = (1._wp - psicomo(i, j)) * Xconvxa * pfu10(i, j)**2        &
     &           * (660._wp / scdms)**0.5_wp

         kwo2  = (1._wp - psicomo(i, j)) * Xconvxa * pfu10(i, j)**2        &
     &           * (660._wp / sco2)**0.5_wp
#ifdef PCFC
         kwcfc12 = (1._wp - psicomo(i, j)) * Xconvxa * pfu10(i, j)**2        &
     &         * (660._wp / sccfc12)**0.5_wp

         kwcfc11 = (1._wp - psicomo(i, j)) * Xconvxa * pfu10(i, j)**2        &
     &         * (660._wp / sccfc11)**0.5_wp
#endif
!js: from cpl_co2
!PW : coupled model kg/kg --> ppm:
! molecular weight dry air : 28.970
! molecular weight CO2     : 44.011

#ifdef __cpl_co2
          atco2 = atm(i,j,iatmco2)
!         KEEP O2 AND N2 FIXED (VALUES SHOULD BE TAKEN FROM BELEG_BGC)
          ato2 = 196800._wp
          atn2 = 802000._wp
          atm(i,j,iatmo2)=ato2
          atm(i,j,iatmn2)=atn2
#else
          atco2 = atm(i,j,iatmco2)
          ato2  = atm(i,j,iatmo2)
          atn2  = atm(i,j,iatmn2)
#ifdef __c_isotopes
          atc13 = atm(i,j,iatmc13)
          atc14 = atm(i,j,iatmc14)
#endif
#endif


! Surface flux of oxygen

         oxflux=kwo2*dtbgc*(ocetra(i,j,1,ioxygen)                      &
     &                      -oxysa*(ato2 / 196800._wp)) ! *ppao(i,j)/101300. ! sea level pressure normalization


         ocetra(i,j,1,ioxygen)=ocetra(i,j,1,ioxygen)                   &
     &                               -oxflux/pddpo(i,j,1)

        bgcflux(i,j,kdpo2)  =196800._wp * (ocetra(i,j,1,ioxygen)/oxysa - 1._wp)          !CDI tinka

! only for mass balance, a neg. oxflux means a gain for the ocean
         o2flux(i,j)=o2flux(i,j)+oxflux
! Surface flux of gaseous nitrogen (same piston velocity as for O2)

         niflux=kwo2*dtbgc*(ocetra(i,j,1,igasnit)                      &
     &                     -anisa*(atn2/802000._wp)) ! *ppao(i,j)/101300.

         ocetra(i,j,1,igasnit)=ocetra(i,j,1,igasnit)                   &
     &                               -niflux/pddpo(i,j,1)
         n2flux(i,j)=n2flux(i,j)+niflux
! Surface flux of laughing gas (same piston velocity as for O2 and N2)

         nlaughflux=kwo2*dtbgc*(ocetra(i,j,1,ian2o)                       &
     &                      -satn2o(i,j)) ! *ppao(i,j)/101300.

         ocetra(i,j,1,ian2o)=ocetra(i,j,1,ian2o)                       &
     &                               -nlaughflux/pddpo(i,j,1)
         n2oflux(i,j)=n2oflux(i,j)+nlaughflux

#ifndef __cpl_co2
         atm(i,j,iatmo2)=atm(i,j,iatmo2) + oxflux *contppm ! only used for diagnostic, not for mass balance
         atm(i,j,iatmn2)=atm(i,j,iatmn2) + niflux *contppm
         atm(i,j,iatmn2)=atm(i,j,iatmn2) + nlaughflux*contppm
#endif

! Surface flux of dms

         dmsflux = kwdms*dtbgc*ocetra(i,j,1,idms)
         ocetra(i,j,1,idms)=ocetra(i,j,1,idms)-dmsflux/pddpo(i,j,1)

! Write output for bgcflux (CDI)
         bgcflux(i,j,ko2flux)  =oxflux/dtbgc               ! CDI tinka
         bgcflux(i,j,kn2flux)  =niflux/dtbgc               ! CDI tinka
         bgcflux(i,j,kn2oflux) =nlaughflux/dtbgc              ! CDI tinka
         bgcflux(i,j,kdmsflux) =dmsflux/dtbgc              ! CDI tinka
!
! Write output for bgcmean
!
         bgct2d(i,j,jo2flux)  =bgct2d(i,j,jo2flux) +oxflux
         bgct2d(i,j,jn2flux)  =bgct2d(i,j,jn2flux) +niflux
         bgct2d(i,j,jn2oflux) =bgct2d(i,j,jn2oflux)+nlaughflux

         bgcm2d(i,j,joxflux)=                             &
     &          bgcm2d(i,j,joxflux)+oxflux
         bgcm2d(i,j,jniflux)=                             &
     &          bgcm2d(i,j,jniflux)+niflux
         bgcm2d(i,j,jdmsflux)=                            &
     &          bgcm2d(i,j,jdmsflux) + dmsflux
         bgcm2d(i,j,jdms)=                                &
     &          bgcm2d(i,j,jdms) + ocetra(i,j,1,idms)


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
         a=c*(2._wp + t2)/(1._wp + t2 + t2 * t1) + akw/h - h &
              + bt/(1._wp + h * akbi) - alk
         dadh=c*(1._wp / (ak2 * (1._wp + t2 + t2 * t1)) &
              - (2._wp + t2) * ((1._wp + 2._wp * t1) * ak2i)/  &
     &       (1._wp + t2 + t2 * t1)**2)                                      &
     &       -akw/h**2 - 1._wp - (bt * akbi) / (1._wp + h * akbi)**2

         dddhhh=a/dadh
         h = MAX(h - dddhhh, 1.e-10_wp) ! Prevent overshooting to negative values at start of iteration

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2._wp + t2)/(1._wp + t2 + t2 * t1) + akw/h - h &
              + bt/(1._wp + h * akbi) - alk
         dadh=c*(1._wp/(ak2*(1._wp + t2 + t2 * t1)) &
              - (2._wp + t2) * ((1._wp + 2._wp * t1) * ak2i)/  &
     &       (1._wp + t2 + t2 * t1)**2)                                      &
     &       - akw / h**2 - 1._wp - (bt * akbi) / (1._wp + h * akbi)**2
         dddhhh=a/dadh
         h = MAX(h - dddhhh, 1.E-10_wp) ! Prevent overshooting to negative values at start of iteration

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2._wp + t2)/(1._wp + t2 + t2 * t1) + akw/h - h &
              + bt/(1._wp + h * akbi) - alk
         dadh = c * (1._wp / (ak2 * (1._wp + t2 + t2 * t1)) &
              - (2._wp + t2) * ((1._wp + 2._wp * t1) * ak2i) /  &
     &       (1._wp + t2 + t2 * t1)**2)                                      &
     &       - akw / h**2 - 1._wp - (bt * akbi) / (1._wp + h * akbi)**2
         dddhhh=a/dadh
         h = MAX(h - dddhhh, 1.e-10_wp) ! Prevent overshooting to negative values at start of iteration
         hi(i,j,k)=h


!
! Calculate pCO2 [ppmv] from total dissolved inorganic carbon (DIC: SCO212)
! the calculation also includes solubility
!
!            co212
         pco2=  c  /((1._wp + ak1 * (1._wp + ak2/h)/h) * ak0)
         suppco2(i,j)=pco2                     !

         fluxd=atm(i,j,iatmco2)*kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.
         fluxu=pco2 *kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.

#ifdef __cpl_co2
#ifdef AVFLUX
! store fluxd and fluxu in case of coupling for correcting ECHAM flux pattern
! with oceanic small scale delta pco2 pattern
         ! FIXME: replace 44.011 with parameter
         co2flux_act(i,j)=(fluxu-fluxd) * 44.011_wp/dtbgc  ! convert to seconds and
                                                      ! from gC to gCO2, unit as co2flux_cpl
#endif
!--      coupled model, flux in kg/m^2/sec from atmospheric model (ECHAM5)
         co2trans(i,j) = kwco2_cpl * ak0
! only for mass balance, a neg. co2flux means a gain for the ocean,
! ECHAM flux, converted to kmol m-2 dt-1
         co2flux(i,j) = co2flux(i,j)+co2flux_cpl(i,j)*dtbgc/44.011_wp
#else
         atm(i,j,iatmco2)=atm(i,j,iatmco2)+(fluxu-fluxd)*contppm        ! only for diagnotic
#endif  /*__cpl_co2*/

! new concentrations ocean (kmol/m3 -->ppm)
          thickness    = (DDPO(i,j,k)                               &    ! still k=1
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA               &
     &                   -SICSNO(I,J)*RHOSNWA)

#ifdef __cpl_co2
#ifndef AVFLUX
         ocetra(i,j,1,isco212)=                                     &
     &   ocetra(i,j,1,isco212)-((co2flux_cpl(i,j)*dtbgc)/(44.011_wp*thickness))

#else
! ECHAM co2flux modified with oceanic small scale structure after do loop
#endif
#else
         ocetra(i,j,1,isco212)=                                     &
     &   ocetra(i,j,1,isco212)+(fluxd-fluxu)/thickness
! only for mass balance, a neg. co2flux means a gain for the ocean
         co2flux(i,j) = co2flux(i,j)+fluxu-fluxd              ! kmol m-2 dt-1
#endif

!proxies d13C stuff
#ifdef __c_isotopes
         Roc13=ocetra(i,j,1,isco213)/(ocetra(i,j,1,isco212)+1e-25)
         Roc14=ocetra(i,j,1,isco214)/(ocetra(i,j,1,isco212)+1e-25)

         rat13=atm(i,j,iatmc13)/atco2
         rat14=atm(i,j,iatmc14)/atco2

         flux13d=atm(i,j,iatmc13)*kwco2*dtbgc*ak0              ! *ppao(i,j)/101300. (atm to ocean)
         flux13u=pco2 *kwco2*dtbgc*ak0*Roc13*0.9935_wp ! *ppao(i,j)/101300. (ocean to atm,
                                                                        !twofold fractionation through evaporation)

         flux14d=atm(i,j,iatmc14)*kwco2*dtbgc*ak0              ! *ppao(i,j)/101300.
         flux14u=pco2 *kwco2*dtbgc*ak0*Roc14*0.987_wp  ! *ppao(i,j)/101300. (twofold fractionation through evaporation)

         atm(i,j,iatmc13)=atm(i,j,iatmc13)+(flux13u-flux13d)*contppm
         atm(i,j,iatmc14)=atm(i,j,iatmc14)+(flux14u-flux14d)*contppm
         ocetra(i,j,1,isco213)=                                     &
     &   ocetra(i,j,1,isco213)+(flux13d-flux13u)/thickness
         ocetra(i,j,1,isco214)=                                     &
     &   ocetra(i,j,1,isco214)+(flux14d-flux14u)/thickness
         c13flux(i,j) = c13flux(i,j)-(flux13d-flux13u)
         c14flux(i,j) = c14flux(i,j)-(flux14d-flux14u)
#endif /*__c_isotopes*/

#ifdef ANTC14
         Roc14=ocetra(i,j,1,iantc14)/(ocetra(i,j,1,isco212)+1e-15)

         fantc14d=atco2*kwco2*dtbgc*ak0*Rbomb(i,j)  ! *ppao(i,j)/101300.
         fantc14u=pco2 *kwco2*dtbgc*ak0*Roc14       ! *ppao(i,j)/101300.

         ocetra(i,j,1,iantc14)=                                     &
     &   ocetra(i,j,1,iantc14)+(fantc14d-fantc14u)/pddpo(i,j,1)

         bgcm2d(i,j,jac14fx) = bgcm2d(i,j,jac14fx)     &
     &                         + fantc14d-fantc14u
         bgcflux(i,j,kac14fx) = (fantc14d-fantc14u)/dtbgc     !CDI tinka
#endif

!
! write output for bgcflux !CDI tinka
!

#ifdef __cpl_co2
        bgcflux(i,j,kco2flux)=co2flux_cpl(i,j)/44.011_wp        !CDI tinka co2flux_cpl in kgCO2/m2/s
                                                             !    converted to kmolC/m2/s
#else
        bgcflux(i,j,kco2flux)=(fluxu-fluxd)/dtbgc             !CDI tinka

        bgcflux(i,j,kco2fxd)=fluxd/dtbgc                      !CDI tinka
#ifdef __c_isotopes
        bgcflux(i,j,kc13flux)=(flux13u-flux13d)/dtbgc               !CDI tinka
        bgcflux(i,j,kc14flux)=(flux14u-flux14d)/dtbgc               !CDI tinka
#endif
#endif
        bgcflux(i,j,kco2fxu)=fluxu/dtbgc                           !CDI tinka
        bgcflux(i,j,kpco2)  =pco2                                  !CDI tinka
        bgcflux(i,j,kdpco2) =pco2-atco2                            !CDI tinka
        bgcflux(i,j,kkwco2) =kwco2*ak0                             !CDI tinka kwco2 in m/s

!
! write output for bgcmean
!
#ifdef __cpl_co2
        bgct2d(i,j,jco2flux)=bgct2d(i,j,jco2flux)+co2flux_cpl(i,j)*dtbgc/44.011_wp
        bgcm2d(i,j,jatmco2)=bgcm2d(i,j,jatmco2)+atm(i,j,iatmco2)
        bgcm2d(i,j,jatmn2)=bgcm2d(i,j,jatmn2)+atm(i,j,iatmn2)
        bgcm2d(i,j,jatmo2)=bgcm2d(i,j,jatmo2)+atm(i,j,iatmo2)
!     bgcm2d(i,j,jco2fxd) is used in avflux_redis
#else
        bgct2d(i,j,jco2flux)=bgct2d(i,j,jco2flux)+(fluxu-fluxd)
        bgcm2d(i,j,jatmco2)=bgcm2d(i,j,jatmco2)+atm(i,j,iatmco2)
        bgcm2d(i,j,jatmn2)=bgcm2d(i,j,jatmn2)+atm(i,j,iatmn2)
        bgcm2d(i,j,jatmo2)=bgcm2d(i,j,jatmo2)+atm(i,j,iatmo2)
        bgcm2d(i,j,jco2fxd)=bgcm2d(i,j,jco2fxd)  +fluxd
#ifdef __c_isotopes
        bgct2d(i,j,jc13flux)=bgct2d(i,j,jc13flux)+(flux13u-flux13d)
        bgct2d(i,j,jc14flux)=bgct2d(i,j,jc14flux)+(flux14u-flux14d)
#endif
#endif
        bgcm2d(i,j,jco2fxu)=bgcm2d(i,j,jco2fxu)    &
     &                              +fluxu
        bgcm2d(i,j,jpco2)  =bgcm2d(i,j,jpco2)      &
     &                              +pco2
        bgcm2d(i,j,jkwco2) =bgcm2d(i,j,jkwco2)     &
     &                              +kwco2*ak0
!end of atm.flux calculation

! for stationary state over millennia
         ocetra(i,j,1,isco212) = ocetra(i,j,1,isco212) + calcinp / thickness
         ocetra(i,j,1,ialkali) = ocetra(i,j,1,ialkali) + 2._wp * calcinp / thickness
         ocetra(i,j,1,idoc)    = ocetra(i,j,1,idoc)    + orginp  / thickness
         ocetra(i,j,1,isilica) = ocetra(i,j,1,isilica) + silinp  / thickness

#ifdef PCFC

!     CFC 11 and 12 Solubilities in seawater
!     ref: Warner & Weiss (1985) , Deep Sea Research, vol32
!     coefficient for solubility in  mol/l/atm
!  ----------------------------------------
!     konstants given for mol/(l * atm) or kmol/(m^3 * atm)
!
!     for CFC 11
!     ----------

      ta = ( ptho(i,j,1) + tmelt)* 0.01_wp

      d  = ( -0.0157274_wp * ta + 0.091459_wp)* ta - 0.142382_wp
      sol_cfc11 = exp ( - 229.9261_wp                        &
     &                  + 319.6552_wp / ta                   &
     &                  + 119.4471_wp * alog ( ta )          &
     &                  - 1.39165_wp  * ta * ta  + psao(i,j,1)* d )
!
!     for CFC/12
!     ----------

      d    = ( -0.0153924_wp * ta + 0.091015_wp)* ta - 0.143566_wp
      sol_cfc12 = exp ( - 218.0971_wp                        &
     &                  + 298.9702_wp / ta                   &
     &                  + 113.8049_wp * alog ( ta )          &
     &                  - 1.39165_wp  * ta * ta  + psao(i,j,1)* d )

!
!     conversion from kmol/(m^3 * atm) to kmol/(m3 * pptv)
!     --------------------------------------------------
      sol_cfc11 = 1.0e-12_wp * sol_cfc11
      sol_cfc12 = 1.0e-12_wp * sol_cfc12


!       F = Kw (Csat - Csurf)
!       Csat = alpha * pCFC * P/Po

      cfc11_atm_1 = cfc11_atm_1s*(1._wp-cfc_int(i,j))+cfc11_atm_1n*cfc_int(i,j)
      cfc11_atm_2 = cfc11_atm_2s*(1._wp-cfc_int(i,j))+cfc11_atm_2n*cfc_int(i,j)
      cfc11_atm_3 = cfc11_atm_3s*(1._wp-cfc_int(i,j))+cfc11_atm_3n*cfc_int(i,j)

      cfc12_atm_1 = cfc12_atm_1s*(1._wp-cfc_int(i,j))+cfc12_atm_1n*cfc_int(i,j)
      cfc12_atm_2 = cfc12_atm_2s*(1._wp-cfc_int(i,j))+cfc12_atm_2n*cfc_int(i,j)
      cfc12_atm_3 = cfc12_atm_3s*(1._wp-cfc_int(i,j))+cfc12_atm_3n*cfc_int(i,j)
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
! write output for CDI
        bgcflux(i,j,kcfc11fx) =flux_cfc11/dtbgc           !CDI tinka
        bgcflux(i,j,kcfc12fx) =flux_cfc12/dtbgc           !CDI tinka
        bgcflux(i,j,kpcfc11 ) =pcfc11              !CDI tinka
        bgcflux(i,j,kpcfc12 ) =pcfc12              !CDI tinka
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


#ifdef __cpl_co2
!     redistribution of echam co2 fluxes by KD SIX (NOV 2007);
!     sbr updates co2flux_cpl, ocetra(i,j,1,isco212) and bgcm2d(i,j,jco2fxd)
#ifdef AVFLUX
      call avflux_redis
#endif
#endif

!     atmospheric change due to fluxes

      globalmean_co2 = 0._wp

      globalmean_n2  = 0._wp
      globalmean_o2  = 0._wp
      CALL global_mean_2d(atm(:,:,iatmco2), globalmean_co2)
      CALL global_mean_2d(atm(:,:,iatmn2 ), globalmean_n2)
      CALL global_mean_2d(atm(:,:,iatmo2 ), globalmean_o2)
      IF (p_pe==p_io .AND. ldtdayc==1) THEN
         WRITE(0,*) 'co2 = ', globalmean_co2, ', n2 = ', globalmean_n2, ', o2 = ', globalmean_o2
      END IF

#ifndef __cpl_co2
      IF(diffat) THEN
        CALL ATMOTR(kpie,kpje,kpke,kplmon,pdlxp,pdlyp)   ! distribution of CO2 emissions to atm(iantco2)
      WRITE(io_stdo_bgc,*) ' in diffat '
      ELSE
        atm(:,:,iatmco2)=globalmean_co2    ! instantaneous mixing to homogenity
        atm(:,:,iatmo2)=globalmean_o2
        atm(:,:,iatmn2)=globalmean_n2
      ENDIF
#endif

! fixme : find a better place
! millennium :  sum over carbon input from sediment
      do j=1,kpje
         do i=1,kpie
               c_int(i,j)=weto(i,j,1)*dlxp(i,j)*dlyp(i,j)*calcinp
               o_int(i,j)=weto(i,j,1)*dlxp(i,j)*dlyp(i,j)*orginp
               s_int(i,j)=weto(i,j,1)*dlxp(i,j)*dlyp(i,j)*silinp
         end do
      end do

      call global_sum_2d_pio(c_int,gs_c)
      call global_sum_2d_pio(o_int,gs_o)
      call global_sum_2d_pio(s_int,gs_s)
      calcinpglint=calcinpglint+gs_c
      orginpglint=orginpglint+gs_o
      silinpglint=silinpglint+gs_s
      WRITE(io_stdo_bgc,*) ' sed reflux ',orginpglint,calcinpglint,silinpglint

      call maschk(kpie,kpje,kpke,76)

!
!     -----------------------------------------------------------------
!*        22. CHEMICAL CONSTANTS - water column

      DO 43 iter=1,3
      DO 43 k=1,kpke
      DO 43 j=1,kpje
      DO 43 i=1,kpie
         IF(pddpo(i,j,k).GT.0.5_wp) THEN
            h=hi(i,j,k)
            c=ocetra(i,j,k,isco212)
            t1=h/ak13(i,j,k)
            t2=h/ak23(i,j,k)
            ak2=ak23(i,j,k)
            akw=akw3(i,j,k)
            bt=rrrcl*psao(i,j,k)
            akb=akb3(i,j,k)
            alk=ocetra(i,j,k,ialkali)
            ! Determine hydrogen ion HI so that ALK(DIC,BT,HI) matches given alk by Newton iteration
            ! Actual mismatch
            a = c * (2._wp + t2) / (1._wp+t2+t2*t1)+akw/h-h+bt/(1._wp+h/akb)-alk
            ! Derivative
            dadh = c * (1._wp / (ak2 * (1._wp + t2 + t2 * t1)) &
                 - (2._wp + t2) * (1._wp / ak2 + 2._wp * t1 / ak2)/ &
     &          (1._wp + t2 + t2 * t1)**2)                          &
     &          - akw / h**2 - 1._wp - (bt / akb) / (1._wp + h / akb)**2
            dddhhh=a/dadh
            h = MAX(h - dddhhh, 1.e-10_wp) ! Prevent overshooting to negative values at start of iteration
            hi(i,j,k) = h
            co3(i,j,k) = c/(1._wp+h*(1._wp+h/ak13(i,j,k))/ak23(i,j,k))

          ENDIF ! wet cell
   43  CONTINUE ! i,j,k,iter loop

#ifdef __c_isotopes
            ! C14 decay, c14dec is set in BELEG_BGC , c14ret = 1-c14dec
            ocetra(:,:,:,isco214)=ocetra(:,:,:,isco214)*c14ret
            ocetra(:,:,:,idet14) =ocetra(:,:,:,idet14) *c14ret
            ocetra(:,:,:,icalc14)=ocetra(:,:,:,icalc14)*c14ret


! C14 decay in the sediment

        do k=1,ks
        do j=1,kpje
        do i=1,kpie
        if(bolay(i,j).gt.0._wp) then
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
        do j=1,kpje
        do i=1,kpie
         iflag(i,j)=0
         bgcflux(i,j,klysokl)=0._wp
         IF(pddpo(i,j,1).GT.0.5_wp) THEN
            supsat=co3(i,j,1)-97._wp*aksp(i,j,1)
            IF (supsat .LT. 0._wp) then
               bgcflux(i,j,klysokl)=ptiestu(1)
               iflag(i,j)=1
            endif
         END IF
        end do
        end do

        DO 11 k=2,kpke
        DO 11 j=1,kpje
        DO 11 i=1,kpie
         IF(pddpo(i,j,k).GT.0.5_wp) THEN
           supsat=co3(i,j,k)-97._wp*aksp(i,j,k)                       ! 97. = 1./1.03e-2 (MEAN TOTAL [CA++] IN SEAWATER [kmol/m3])
           undsa = MAX(0._wp, -supsat)
           if(supsat .lt.0._wp .and. iflag(i,j).eq.0) then
            iflag(i,j)=1
            supsatup=co3(i,j,k-1)-97._wp*aksp(i,j,k-1)
            depthdiff=.5_wp * (pddpo(i,j,k)+pddpo(i,j,k-1))
            satdiff=supsatup-supsat
            bgcflux(i,j,klysokl)=ptiestu(k-1)+depthdiff*(supsatup/satdiff)  ! depth of lysokline
           end if

           dissol=MIN(undsa,dremcalc*ocetra(i,j,k,icalc))
#ifdef __c_isotopes
           r13=dissol*ocetra(i,j,k,icalc13)                        &
      &             /(ocetra(i,j,k,icalc)+1.e-25_wp)
           r14=dissol*ocetra(i,j,k,icalc14)                        &
      &             /(ocetra(i,j,k,icalc)+1.e-25_wp)
#endif
           ocetra(i,j,k,icalc)=ocetra(i,j,k,icalc)-dissol
           ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)+2._wp*dissol

           ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+dissol

           bgc_o_pro(i,j,k,kdissol)  =  dissol/dtbgc                      !CDI tinka
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
         pco2=c/((1._wp+ak1*(1._wp+ak2/h)/h)*ak0)


         fluxd=atm(i,j,iatmco2)*kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.
         fluxu=pco2 *kwco2*dtbgc*ak0 ! *ppao(i,j)/101300.

#ifdef __cpl_co2
         ts1(itsco2f,l+1,lts1) = ts1(itsco2f,l+1,lts1)                 &
     &                         + co2flux_cpl(i,j)*dtbgc/44.011_wp
#else
         ts1(itsco2f,l+1,lts1) = ts1(itsco2f,l+1,lts1)                 &
     &                         + fluxu - fluxd
#endif
         ts1(itspco2,l+1,lts1) = ts1(itspco2,l+1,lts1)                 &
     &                         + pco2
         ts1(itsatm,l+1,lts1)  = ts1(itsatm,l+1,lts1)                  &
     &                         + atm(i,j,iatmco2)

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
                  IF( hi(i,j,k) .LT. 0.0_wp ) THEN
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
