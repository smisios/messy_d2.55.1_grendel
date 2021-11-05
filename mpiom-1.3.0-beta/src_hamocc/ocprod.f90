      SUBROUTINE OCPROD(kpie,kpje,kpke,ptho,pddpo,                     &
     &                    pdlxp,pdlyp,pdpio,kplmon)

!$Source: /scratch/local1/m212047/patrick/SRC_MPI/src_hamocc/RCS/ocprod.f90,v $\\
!$Revision: 1.1 $\\
!$Date: 2005/01/28 08:37:45 $\\

!**********************************************************************
!
!**** *OCPROD* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Purpose
!     -------
!     compute biological production, settling of debris, and related biogeochemistry
!
!     Method:
!     ------
!     kchck=1 can be used to check max/min of bgc arrays on wet/dry cells.
!     Note: prosil is used only for k=1,2. It is adressed, however, for
!           k=1,4 in loop 100. To save memory,  js:note seems dated???
!
!     _ant fields are natural PLUS anthropogenic (not anthropogenic only!!!)
!
!     *CALL*       *OCPROD*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *pdpio*   - inverse thickness of grid cell (3rd dimension)[m].
!
!     Externals
!     ---------
!     .
!**********************************************************************
! nocetra is the number of all BGC elements (size of ocetra(,,,l))

      USE mo_timeser_bgc
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      use mo_param1_bgc 

      USE mo_control_bgc
      USE mo_bgcmean

      use mo_parallel

implicit none

      INTEGER :: kplmon,kpie,kpje,kpke
      INTEGER :: i,j,k,l
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdpio(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: abs_bgc(kpie,kpje,kpke)                       !js name implies "absorbtion bgc", but it is the light fraction 
                                                            !                                  (surface layer = 1.) 
      REAL :: dmsp1,dmsp2,dmsp3,dmsp4,dmsp5,dmsp6
      REAL :: atten,avphy,avanut,avanfe,pho,xa,xn,ya,yn,phosy,         &
     &        volcell,avgra,grazing,avsil,graton,                      &
     &        gratpoc,grawa,bacfra,phymor,zoomor,excdoc,exud,          &
     &        export, delsil, delcar, sterph, sterzo, remin,           &
     &        docrem, opalrem, remin2o, aou, refra
      REAL :: epfac, tloc, eppz15    ! not used
      
      REAL :: dustinp, dustinp_volcano 
      REAL :: gutscale
      REAL :: fopa, fdet, fcal
      REAL :: absorption
      REAL :: dmsprod,dms_bac,dms_uv 
      REAL :: detref, detrl
      INTEGER :: volcano
#ifdef __c_isotopes
      REAL :: rem13,rem14
      REAL :: rl13, rl14
      REAL :: rocean13, rocean14, flui13, flui14
      REAL :: fcal13,fcal14
      REAL :: d13C, d14C
#endif
#ifdef AGG
      REAL :: fphy                               ! sinking speed wphy * concentration of phytoplankton (for time series only?)
      REAL :: wmass(kpie,kpje,kpke)              ! sinking speed for 'mass of aggregates'
      REAL :: wnumb(kpie,kpje,kpke)              ! sinking speed for 'numbers of aggregates'
      REAL :: aggregate(kpie,kpje,kpke)          ! aggregation (should be renamed)
      REAL :: dustagg(kpie,kpje,kpke)            ! aggregation of dust

      REAL :: avmass, avsize, avsmin1, avsmax1, avnos, anosloss     
      INTEGER :: nosin1, nosde1, nosin2, nosde2, nosin3, nosde3     ! nosin = nos increase, nosde = nos decrease
      REAL :: zmornos, avm, avn, eps, e1,e2,e3,e4,es1,es3
      REAL :: TopM, TopF, snow,fshear,sagg1,sagg2,sagg4
      REAL :: sett_agg,shear_agg,effsti,dfirst,dshagg,dsett
      REAL :: checksize,nacheck,flar1,flar2,flar3
      REAL :: fTSFac,fTMFac,fTopF,fTopM,wphy,wphyup,wnosup,wnos

#endif 

      volcano=0           ! volcano=0 no dust from volcano, =1 ashes from volcano in dustinp_volcano
!
! Constant parameters
!
! parameter definition in beleg_bgc.f90
      call maschk(kpie,kpje,kpke,0)

      dmsp6=dmspar(6)
      dmsp5=dmspar(5)
      dmsp4=dmspar(4)
      dmsp3=dmspar(3)
      dmsp2=dmspar(2)
      dmsp1=dmspar(1)

! Calculate swr absorption by water and phytoplankton

! Almost half of the SWR is absorbed in the surface layer (0.04*10m). --> *0.6 is transmitted (js: for atten_c=0!)
! 1) Upper 2 layers
      DO j=1,kpje
      DO i=1,kpie
        abs_bgc(i,j,1)=1.                                 ! surface layer has full light available
#ifdef FB_BGC_OCE     
          abs_oce(i,j,1)=1.
#endif       
        IF(pddpo(i,j,1).GT.0.5) THEN


        atten=atten_w+atten_c*ocetra(i,j,1,iphy)        ! atten_w = 0.04 m^-1, atten_c=7.32e5 (defined in beleg_bgc)
                                                          ! for ocetra(iphy)=5e-7 attenuation is 0.36/m : seems too much
          absorption= exp(-atten*pddpo(i,j,1))            ! absorption in surface layer: 0.67 for atten_c=0, sfc layer 10m
        abs_bgc(i,j,2)=absorption                       ! second layer 

! (Implicit would be faster:)
!        abs_bgc(i,j,2)=1./(1+atten*pddpo(i,j,1))

#ifdef FB_BGC_OCE
          abs_oce(i,j,2)=atten_f*absorption               ! atten_f =0.4 (beleg_bgc.f90)
!            abs_oce(i,j,2)=0.5/(1+atten*pddpo(i,j,1))
#endif 

!            bgcm3d(i,j,1,jatten)=                            &
!     &                           bgcm3d(i,j,1,jatten) +atten
     
        ENDIF                  
      ENDDO
      ENDDO
      call maschk(kpie,kpje,kpke,-2)       ! entrance of ocprod
      
! 2) water column
      DO k=3,kpke           ! (js) why run over whole water colum? kwrbioz sufficient
      DO j=1,kpje
      DO i=1,kpie
        IF(pddpo(i,j,k-1).GT.0.5) THEN  ! wet point

!        atten=0.04+0.03*2.44e7*ocetra(i,j,k-1,iphy)
!        abs_bgc(i,j,k)=abs_bgc(i,j,k-1)/(1+atten*pddpo(i,j,k-1))

        atten=atten_w+atten_c*ocetra(i,j,k-1,iphy)
          absorption= exp(-atten*pddpo(i,j,k-1))
        abs_bgc(i,j,k)=abs_bgc(i,j,k-1)*absorption
#ifdef FB_BGC_OCE     
!            abs_oce(i,j,k)=abs_oce(i,j,k-1)/(1+atten*pddpo(i,j,k-1))
            abs_oce(i,j,k)=abs_oce(i,j,k-1)*absorption
#endif 
          IF(k-1.LE.kwrbioz) THEN
!            bgcm3d(i,j,k-1,jatten)=                            &
!     &                           bgcm3d(i,j,k-1,jatten) +atten
          ENDIF
        ENDIF                  
      ENDDO
      ENDDO
      ENDDO
      call maschk(kpie,kpje,kpke,-3)

! dust flux from the atmosphere to the surface layer; 
! dust fields are monthly mean values in units of kg/m2/year
! dissolved iron is a fixed fraction (3.5%), and immediately released
! 1% of the iron input is dissolved [='bio-available'] (see perc_diron in beleg_bgc)

      dustinp_volcano=0.

      do j=1,kpje
      do i=1,kpie
       if(pddpo(i,j,1).gt.0.5) then   ! wet point
#ifdef __cpl_dust
        dustinp=dustdep(i,j)/365.*dtb*pdpio(i,j,1)

        if (volcano.eq.1) then
          dustinp_volcano=dustdep(i,j)/365.*dtb*pdpio(i,j,1)
          dustinp=dusty(i,j,kplmon)/365.*dtb*pdpio(i,j,1)
        endif
#else
        dustinp=dusty(i,j,kplmon)/365.*dtb*pdpio(i,j,1)
#endif
        ocetra(i,j,1,ifdust)=ocetra(i,j,1,ifdust)+dustinp + dustinp_volcano 
        ocetra(i,j,1,iiron)=ocetra(i,j,1,iiron)+dustinp*perc_diron  &! perc_diron =.035*.01 /55.85 ! why *0.8 in IPCC_HAM?
    &                                          +dustinp_volcano*perc_diron*0.8/3.5  ! volcanic ash has lower iron content

       endif      
      enddo
      enddo
!
! sum-up 3d fields in euphotic zone (averaged monthly in avrg_bgcmean_2d.f90) 
! these go into bgcmean_bioz files
! 
      DO k=1,kwrbioz
      DO j=1,kpje
      DO i=1,kpie
         IF(pddpo(i,j,k).GT.0.5) THEN

            bgcm3d(i,j,k,jphyto)  =                   &
     &      bgcm3d(i,j,k,jphyto)  + ocetra(i,j,k,iphy)
            bgcm3d(i,j,k,jgrazer) =                   &
     &      bgcm3d(i,j,k,jgrazer) + ocetra(i,j,k,izoo)
            bgcm3d(i,j,k,jphosph) =                   &
     &      bgcm3d(i,j,k,jphosph) + ocetra(i,j,k,iphosph)
            bgcm3d(i,j,k,joxygen) =                   &
     &      bgcm3d(i,j,k,joxygen) + ocetra(i,j,k,ioxygen)
            bgcm3d(i,j,k,jiron)   =                   &
     &      bgcm3d(i,j,k,jiron)   + ocetra(i,j,k,iiron)
            bgcm3d(i,j,k,jano3)   =                   &
     &      bgcm3d(i,j,k,jano3)   + ocetra(i,j,k,iano3)
            bgcm3d(i,j,k,jalkali) =                   &
     &      bgcm3d(i,j,k,jalkali) + ocetra(i,j,k,ialkali)        
            bgcm3d(i,j,k,jsilica) =                   &
     &      bgcm3d(i,j,k,jsilica) + ocetra(i,j,k,isilica)      
            bgcm3d(i,j,k,jdic)    =                   &
     &      bgcm3d(i,j,k,jdic)    + ocetra(i,j,k,isco212)

! delta notation for 13C, 14C (js 6.3.2006)
#ifdef __c_isotopes
            d13C = ((ocetra(i,j,k,isco213)/ocetra(i,j,k,isco212)) -1.) *1000.
            bgcm3d(i,j,k,jdic13_t)  = bgcm3d(i,j,k,jdic13_t) +d13C

!           bgcm3d(i,j,k,jdic13)  =                     &
!    &      bgcm3d(i,j,k,jdic13)  + ocetra(i,j,k,isco213)

            d14C = ((ocetra(i,j,k,isco214)/ocetra(i,j,k,isco212)) -1.) *1000.
            D14C = d14C - 2.* (d13C + 25.) * (1.+ d14C *1.e-3)
            bgcm3d(i,j,k,jdic14_t)  = bgcm3d(i,j,k,jdic14_t)  + D14C

!           bgcm3d(i,j,k,jdic14)  =                     &
!    &      bgcm3d(i,j,k,jdic14)  + ocetra(i,j,k,isco214)
#endif

            bgcm3d(i,j,k,jdoc)    =                   &
     &      bgcm3d(i,j,k,jdoc)    + ocetra(i,j,k,idoc)        
!            bgcm3d(i,j,k,jdms)    =                   &      ! these are not in write_bgcmean etc.
!     &      bgcm3d(i,j,k,jdms)    + ocetra(i,j,k,idms) 
!            bgcm3d(i,j,k,jpoc)    =                   &
!     &      bgcm3d(i,j,k,jpoc)    + ocetra(i,j,k,idet)
!            bgcm3d(i,j,k,jcalc)    =                   &
!     &      bgcm3d(i,j,k,jcalc)    + ocetra(i,j,k,icalc)   
!            bgcm3d(i,j,k,jopal)    =                   &
!     &      bgcm3d(i,j,k,jopal)    + ocetra(i,j,k,iopal)   
#ifdef AGG
            bgcm3d(i,j,k,jnos)    =                   &
          bgcm3d(i,j,k,jnos)    + ocetra(i,j,k,inos)
#endif     
         ENDIF
      ENDDO
      ENDDO
      ENDDO
          
! sum-up 3d tracers, averaged annually in avrg_bgcmean_3d.f90
! these go into bgcmean_3d files
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
         IF(pddpo(i,j,k).GT.0.5) THEN
       
             bgct3d(i,j,k,jphosph_t) =                                     &
     &       bgct3d(i,j,k,jphosph_t) + ocetra(i,j,k,iphosph)
             bgct3d(i,j,k,joxygen_t) =                                     &
     &       bgct3d(i,j,k,joxygen_t) + ocetra(i,j,k,ioxygen)
             bgct3d(i,j,k,jiron_t)   =                                     &
     &       bgct3d(i,j,k,jiron_t)   + ocetra(i,j,k,iiron)
             bgct3d(i,j,k,jano3_t)   =                                     &
     &       bgct3d(i,j,k,jano3_t)   + ocetra(i,j,k,iano3)
             bgct3d(i,j,k,jalkali_t) =                                     &
     &       bgct3d(i,j,k,jalkali_t) + ocetra(i,j,k,ialkali)       
             bgct3d(i,j,k,jsilica_t) =                                     &
     &       bgct3d(i,j,k,jsilica_t) + ocetra(i,j,k,isilica)   
             bgct3d(i,j,k,jdic_t)    =                                     &
     &       bgct3d(i,j,k,jdic_t)    + ocetra(i,j,k,isco212)     

! 13C and 14C in delta notation

#ifdef __c_isotopes
             d13C =  ((ocetra(i,j,k,isco213)/ocetra(i,j,k,isco212)) -1.) *1000.
             bgct3d(i,j,k,jdic13_t)  = bgct3d(i,j,k,jdic13_t) + d13C

!            bgct3d(i,j,k,jdic13_t)  =                                         &
!           &bgct3d (i, j, k, jdic13_t) + ocetra (i, j, k, isco213)

!            d14C =  ((ocetra(i,j,k,isco214)/ocetra(i,j,k,isco212))/bifr14 -1.) *1000.
             d14C =  ((ocetra(i,j,k,isco214)/ocetra(i,j,k,isco212)) -1.) *1000.    ! 10.4.06
             D14C =  d14C - 2* (d13C + 25) * (1+ d14C *1.e-3)
             bgct3d(i,j,k,jdic14_t)  = bgct3d(i,j,k,jdic14_t)  + d14C

!            bgct3d(i,j,k,jdic14_t)  =                                         &
!    &       bgct3d(i,j,k,jdic14_t)  + ocetra(i,j,k,isco214)
#endif

             bgct3d(i,j,k,jdoc_t)    =                                     &
     &       bgct3d(i,j,k,jdoc_t)    + ocetra(i,j,k,idoc)      

             bgct3d(i,j,k,jpoc_t)    =                                     &
     &       bgct3d(i,j,k,jpoc_t)    + ocetra(i,j,k,idet)    
! add poc_13, poc14 (also needed in mo_bgcmean.f90)

             bgct3d(i,j,k,jcalc_t)   =                                     &
     &       bgct3d(i,j,k,jcalc_t)   + ocetra(i,j,k,icalc)  
             bgct3d(i,j,k,jopal_t)   =                                     &
     &       bgct3d(i,j,k,jopal_t)   + ocetra(i,j,k,iopal)  

#ifdef ANTC14
           bgct3d(i,j,k,jac14_t) =                                       &
     &           bgct3d(i,j,k,jac14_t) + ocetra(i,j,k,iantc14)
#endif
#ifdef PCFC

             bgct3d(i,j,k,jcfc11_t)    =                                   &
     &           bgct3d(i,j,k,jcfc11_t)    + ocetra(i,j,k,icfc11)
             bgct3d(i,j,k,jcfc12_t)    =                                   &
     &           bgct3d(i,j,k,jcfc12_t)    + ocetra(i,j,k,icfc12)
#endif         

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      call maschk(kpie,kpje,kpke,-4)

!js sediment
      DO k=1,ks
      DO j=1,kpje
      DO i=1,kpie
         IF(bolay(i,j).GT.0.) THEN 
             bgct_sed(i,j,k,jpowaic)  =                      &
     &       bgct_sed(i,j,k,jpowaic)  + powtra(i,j,k,ipowaic)
             bgct_sed(i,j,k,jpowaal)  =                      &
     &       bgct_sed(i,j,k,jpowaal)  + powtra(i,j,k,ipowaal)
             bgct_sed(i,j,k,jpowaph)  =                      &
     &       bgct_sed(i,j,k,jpowaph)  + powtra(i,j,k,ipowaph)
             bgct_sed(i,j,k,jpowaox)  =                      &
     &       bgct_sed(i,j,k,jpowaox)  + powtra(i,j,k,ipowaox)
             bgct_sed(i,j,k,jpown2)   =                      &
     &       bgct_sed(i,j,k,jpown2)   + powtra(i,j,k,ipown2)
             bgct_sed(i,j,k,jpowno3)  =                      &
     &       bgct_sed(i,j,k,jpowno3)  + powtra(i,j,k,ipowno3)
             bgct_sed(i,j,k,jpowasi)  =                      &
     &       bgct_sed(i,j,k,jpowasi)  + powtra(i,j,k,ipowasi)
             bgct_sed(i,j,k,jssso12)  =                      &
     &       bgct_sed(i,j,k,jssso12)  + sedlay(i,j,k,issso12)
             bgct_sed(i,j,k,jssssil)  =                      &
     &       bgct_sed(i,j,k,jssssil)  + sedlay(i,j,k,issssil)
             bgct_sed(i,j,k,jsssc12)  =                      &
     &       bgct_sed(i,j,k,jsssc12)  + sedlay(i,j,k,isssc12)
             bgct_sed(i,j,k,jssster)  =                      &
     &       bgct_sed(i,j,k,jssster)  + sedlay(i,j,k,issster)


         ENDIF
      ENDDO
      ENDDO
      ENDDO

      call maschk(kpie,kpje,kpke,-5)

!mz_ap_20070626-
#ifndef MESSY
!
! Sampling timeseries-1 : global inventory (this will be stored in same file as station data with index 0)
!
      DO k=1,kpke
      DO j=2,kpje-1
         DO i=2,kpie-1
            IF(pddpo(i,j,k).GT.0.5) THEN
               pho=pi_alpha*strahl(i,j)*abs_bgc(i,j,k)        &        ! strahl =swr*(1-icecover_fraction) pi_alpha=0.02*dtb
     & *(1.+0.06*ptho(i,j,k)*(1.+0.03*ptho(i,j,k)))        !17022006 ! temperature dependency from ernst, 0.06^T -> ln exp...
                                                                     ! emr: "optimum" 0.8, radiation ~400Wm^-2 ->.002*400=.8
               avphy=MAX(phytomi,ocetra(i,j,k,iphy)) 
               avanut=MAX(0.,MIN(ocetra(i,j,k,iphosph),                &  ! available nutrients [kmol P/m3] (phosphate
     &             rnoi*ocetra(i,j,k,iano3)))                             ! rnoi=1./16.                     + nitrate
               avanfe=MAX(0.,MIN(avanut,ocetra(i,j,k,iiron)/riron))       ! riron=6.1e-4                    iron limitation
               xa=avanfe                                             ! try 10% iron-independent plankton: .9*avanfe+.1*avanut???
               xn=xa/(1.+pho*avphy/(xa+bkphy))
               phosy=MAX(0.,xa-xn)

               volcell = pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)*1.e-6          ! why *1.e-6? pdlxp [m], pdlyp [m], pddpo [m]
                                                                           ! does not fit unit in output file [kmol/m3]
                                                                           ! unit in output file should be [10^6 kmol P /day]

               ts1(itssco212,1,lts1) = ts1(itssco212,1,lts1)           &   ! lts1 counts time steps for sampling
     &                               + ocetra(i,j,k,isco212)*volcell
               ts1(itsphosy,1,lts1)  = ts1(itsphosy,1,lts1)            & 
     &                               + phosy*volcell
#ifdef __c_isotopes
               ts1(itssco213,1,lts1) = ts1(itssco213,1,lts1)           &
     &                               + ocetra(i,j,k,isco213)*volcell
               ts1(itssco214,1,lts1) = ts1(itssco214,1,lts1)           &
     &                               + ocetra(i,j,k,isco214)*volcell
#endif
!js: was itsphosy instead of itsphosph changed 6.10.2005
               ts1(itsphosph,1,lts1) = ts1(itsphosph,1,lts1)           & 
     &                               + ocetra(i,j,k,iphosph)*volcell
               ts1(itsoxygen,1,lts1) = ts1(itsoxygen,1,lts1)           & 
     &                               + ocetra(i,j,k,ioxygen)*volcell
               ts1(itsgasnit,1,lts1) = ts1(itsgasnit,1,lts1)           & 
     &                               + ocetra(i,j,k,igasnit)*volcell
               ts1(itsano3,1,lts1)   = ts1(itsano3,1,lts1)             &
     &                               + ocetra(i,j,k,iano3)*volcell
               ts1(itssilica,1,lts1) = ts1(itssilica,1,lts1)           & 
     &                               + ocetra(i,j,k,isilica)*volcell
               ts1(itsdoc,1,lts1)    = ts1(itsdoc,1,lts1)              & 
     &                               + ocetra(i,j,k,idoc)*volcell
               ts1(itsphy,1,lts1)    = ts1(itsphy,1,lts1)              & 
     &                               + ocetra(i,j,k,iphy)*volcell
               ts1(itszoo,1,lts1)    = ts1(itszoo,1,lts1)              & 
     &                               + ocetra(i,j,k,izoo)*volcell
               ts1(itsdet,1,lts1)    = ts1(itsdet,1,lts1)              & 
     &                               + ocetra(i,j,k,idet)*volcell
               ts1(itscalc,1,lts1)   = ts1(itscalc,1,lts1)             & 
     &                               + ocetra(i,j,k,icalc)*volcell
               ts1(itsopal,1,lts1)   = ts1(itsopal,1,lts1)             & 
     &                               + ocetra(i,j,k,iopal)*volcell
               ts1(itsiron,1,lts1)   = ts1(itsiron,1,lts1)             & 
     &                               + ocetra(i,j,k,iiron)*volcell
            ENDIF
         ENDDO
      ENDDO
      ENDDO
!
      call maschk(kpie,kpje,kpke,-6)
! Sampling timeseries-1 : concentrations at specific positions as sum over euphotic layer
!
      DO k = 1,kwrbioz               ! js: was 4, changed to kwrbioz 9.2.2006
      DO l=1,nts               ! no. of stations
         i = its1(l)-p_ioff
         j = jts1(l)-p_joff
         IF(i<=1 .OR. i>=kpie .OR. j<=1 .OR. j>=kpje) CYCLE

         pho=pi_alpha*strahl(i,j)*abs_bgc(i,j,k)             & 
     & *(1.+0.06*ptho(i,j,k)*(1.+0.03*ptho(i,j,k)))                  ! temperature dependency from ernst, 0.06^T -> ln exp...

         avphy=MAX(phytomi,ocetra(i,j,k,iphy))  
         avgra=MAX(grami,ocetra(i,j,k,izoo))  
         avanut=MAX(0.,MIN(ocetra(i,j,k,iphosph),                     &
     &          rnoi*ocetra(i,j,k,iano3)))
         avanfe=MAX(0.,MIN(avanut,ocetra(i,j,k,iiron)/riron))
         xa=avanfe
         xn=xa/(1.+pho*avphy/(xa+bkphy)) 
         phosy=MAX(0.,xa-xn)

         ya=avphy+phosy
         yn=(ya+grazra*avgra*phytomi/(avphy+bkzoo))                   &
     &            /(1.+grazra*avgra/(avphy+bkzoo)) 
         grazing=MAX(0.,ya-yn)

         ts1(itsphosph,l+1,lts1) = ts1(itsphosph,l+1,lts1)            &
     &                           + ocetra(i,j,k,iphosph)*pddpo(i,j,k)
         ts1(itsopal,l+1,lts1)   = ts1(itsopal,l+1,lts1)              &
     &                           + ocetra(i,j,k,iopal)*pddpo(i,j,k)
         ts1(itssilica,l+1,lts1) = ts1(itssilica,l+1,lts1)            &
     &                           + ocetra(i,j,k,isilica)*pddpo(i,j,k)
         ts1(itsphy,l+1,lts1)    = ts1(itsphy,l+1,lts1)               &
     &                           + ocetra(i,j,k,iphy)*pddpo(i,j,k)
         ts1(itsdet,l+1,lts1)    = ts1(itsdet,l+1,lts1)               &
     &                           + ocetra(i,j,k,idet)*pddpo(i,j,k) 
#ifdef AGG
         ts1(itsnos,l+1,lts1)    = ts1(itsnos,l+1,lts1)               &
     &                           + ocetra(i,j,k,inos)*pddpo(i,j,k)
#endif
         ts1(itsphosy,l+1,lts1)  = ts1(itsphosy,l+1,lts1)             & 
                               + phosy*pddpo(i,j,k)
         ts1(itszoo,l+1,lts1)    = ts1(itszoo,l+1,lts1)               &
     &                           + ocetra(i,j,k,izoo)*pddpo(i,j,k)
         ts1(itssco212,l+1,lts1) = ts1(itssco212,l+1,lts1)            & 
     &                           + ocetra(i,j,k,isco212)*pddpo(i,j,k)
! 13C, 14C in delta notation?
#ifdef __c_isotopes
         ts1(itssco213,l+1,lts1) = ts1(itssco213,l+1,lts1)            &
     &                           + ocetra(i,j,k,isco213)*pddpo(i,j,k)
         ts1(itssco214,l+1,lts1) = ts1(itssco214,l+1,lts1)            &
     &                           + ocetra(i,j,k,isco214)*pddpo(i,j,k)
#endif
         ts1(itsdoc,l+1,lts1)    = ts1(itsdoc,l+1,lts1)               & 
     &                           + ocetra(i,j,k,idoc)*pddpo(i,j,k)
         ts1(itsiron,l+1,lts1)   = ts1(itsiron,l+1,lts1)              & 
     &                           + ocetra(i,j,k,iiron)*pddpo(i,j,k)     
      ENDDO
      ENDDO
! end time series (output only), below is the relevant part
#endif

#ifdef AGG
!***********************************************************************
!
!  special resetting for particle numbers, that sets their concentration
!  (particles per volume, ocetra(inos)) depending on the mass of marine snow:
!
!  Compartments have already been set to 0 in 
!  ADVECTION_BGC.h and OCTDIFF_BGC.h. js: ???
!
!  Ensure that if there is no mass, there are no particles, and 
!  that the number of particles is in the right range (this is crude, but
!  is supposed to be needed only due to numerical errors such as truncation or 
!  overshoots during advection)
!
! (1) avnos<<avmass, such that eps = FractDim + 1: increase numbers
!     such that eps = FractDim + 1 + safe (currently set to 1.e-6 in BELEG_BGC) 
!
! (2) avnos>>avmass, such that  Nbar (=Mass/Nos/cellmass) <=1: decrease numbers
!     such that Nbar=1.1 (i.e. 1.1 cells per aggregate, set in BELEG_BGC) 

       DO  k=1,kpke
         DO j=1,kpje
           DO i=1,kpie

            IF(pddpo(i,j,k).GT.0.5) THEN               ! wet cell

             avmass = ocetra(i,j,k,iphy) + ocetra(i,j,k,idet)
             snow = avmass*1.e+6                                     ! why *1.e6??

! look for max. and min average size = Nbar             
             if(ocetra(i,j,k,inos).gt.0.) then

                avsize=snow/ocetra(i,j,k,inos)/cellmass              ! js: not further used other than below 2 lines
                avsmin1=MIN(avsize,avsmin1)                          ! js: not used
                avsmax1=MAX(avsize,avsmax1)                          ! js: not used

             endif

! check whether the numbers have to be decreased or increased     

             if (snow*pupper.gt.ocetra(i,j,k,inos)) then
               nosin1 = nosin1 + 1                                  ! counter 'nos increase' (not further used), not set to 0.
             endif

             if (snow*plower.lt.ocetra(i,j,k,inos)) then
               nosde1 = nosde1 + 1                                  ! counter 'nos decrease' (not further used)
             endif

             ocetra(i,j,k,inos) = MAX( snow*pupper, ocetra(i,j,k,inos)) 
             ocetra(i,j,k,inos) = MIN( snow*plower, ocetra(i,j,k,inos))    !js (MAX/MIN correct?)

            ENDIF     ! endif wet cell
           ENDDO
       ENDDO
      ENDDO
         
#endif  /*AGG*/

      call maschk(kpie,kpje,kpke,1)

!
! Biological productivity in the euphotic zone (upper 90m)
!
      DO 100 K=1,kwrbioz

      DO 1 j=1,kpje
      DO 1 i=1,kpie

      IF(pddpo(i,j,k).GT.0.5) THEN

#ifdef AGG
         avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/


#ifdef __c_isotopes
         rocean13=ocetra(i,j,k,isco213)/ocetra(i,j,k,isco212)                ! "ratio oceanic 13C/12C" (at least in carchm)
         rocean14=ocetra(i,j,k,isco214)/ocetra(i,j,k,isco212)                ! ratio oceanic 14C/12C
#endif

         avphy=MAX(phytomi,ocetra(i,j,k,iphy))                            ! 'available' phytoplankton
         avgra=MAX(grami,ocetra(i,j,k,izoo))                              ! 'available' zooplankton
         avsil=MAX(0.,ocetra(i,j,k,isilica))                              ! available silicate
         avanut=MAX(0.,MIN(ocetra(i,j,k,iphosph),                      &  ! available nutrients (phosphate   [kmol P /m3]
     &          rnoi*ocetra(i,j,k,iano3)))                                !                     + nitrate)
         avanfe=MAX(0.,MIN(avanut,ocetra(i,j,k,iiron)/riron))             ! available iron

         pho=pi_alpha*strahl(i,j)*abs_bgc(i,j,k)       &                  ! biological production
     & *(1.+0.06*ptho(i,j,k)*(1.+0.03*ptho(i,j,k)))      !js17022006   ! temperature dependency from ernst, 0.06^T -> ln exp...
                                                                          ! * abs_bgc: absorbtion coefficient
         xa=avanfe
         xn = xa / (1.+ pho*avphy / (xa+bkphy) )                          ! bkphy = half saturation constant
         phosy=MAX(0.,xa-xn)                                              ! photo synthesis
         xn=MAX(xn,1.e-10)
         ya=avphy+phosy                                                   ! new phytoplankton concentration before grazing
         yn=(ya+grazra*avgra*phytomi/(avphy+bkzoo))                    &  ! grazing
     &            /(1.+grazra*avgra/(avphy+bkzoo))
         grazing=MAX(0.,ya-yn)                                            ! what about grazing below euphotic zone?
         graton=epsher*(1.-zinges)*grazing                                ! "grazing to (re-dissolved) nutrients" 
         gratpoc=(1.-epsher)*grazing                                      ! epsher=0.8 "grazing to POC"
         grawa=epsher*zinges*grazing                                      ! grazer 'wachstum(?)'

         bacfra=remido*ocetra(i,j,k,idoc)                                 ! remido = remineralization rate of DOM
         phymor=dyphy*MAX(0.,(ocetra(i,j,k,iphy)-2.*phytomi))             ! phytoplankton mortality dyphy=.008*dt
         zoomor=spemor*MAX(0.,(ocetra(i,j,k,izoo)-2.*grami))              ! zooplankton mortality  
         excdoc=gammaz*MAX(0.,(ocetra(i,j,k,izoo)-2.*grami))              ! excretion to DOC (zooplankton)
         exud=gammap*MAX(0.,(ocetra(i,j,k,iphy)-2.*phytomi))              ! exudation to DOC (phytoplankton)

         ocetra(i,j,k,iphosph)=                                        &
     &       ocetra(i,j,k,iphosph) + bacfra- phosy+graton+ecan*zoomor

         ocetra(i,j,k,iano3)=                                          &
     &       ocetra(i,j,k,iano3)+(bacfra-phosy+graton+ecan*zoomor)*rnit

         export= zoomor*(1.-ecan) + phymor + gratpoc                   ! ecan=.95, gratpoc= .2*grazing [P-units]

         ocetra(i,j,k,idet)=ocetra(i,j,k,idet) + export                ! k=1,8

! new from emr version 4.5.06
#ifdef __c_isotopes
         ocetra(i,j,k,idet13)=ocetra(i,j,k,idet13) + export*rcar*bifr13
         ocetra(i,j,k,idet14)=ocetra(i,j,k,idet14) + export*rcar*bifr14
         ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)-export*rcar*bifr14

! 13C, 14C 'removal' by exportproduction

         flui13=max(rocean13-1.,0.)                                    ! assumes rocean >1 in euphotic layer 
         flui14=max(rocean14-1.,0.)                                    ! flui means what? (flui14 never used) flux...?

! are these two needed? ('delcar' is removed later) (also there is no equivalent for isco212)
! short term removal from dissolved pool "for gas exchange" (in carchm), effect is fairly small
! 12C in P-units, 13C/14C in C-units (*122)  

         ocetra(i,j,k,isco213)=ocetra(i,j,k,isco213)                    &
     &                        -rcar*export* bifr13 * flui13              !bifr = biogenic fractionation (0.98)  ! *flui13
         ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)                    &
     &                        -rcar*export* bifr14 * flui14              ! *rocean14 (--> flui14?) 4.5.06

         ocetra(i,j,k,idet13)=ocetra(i,j,k,idet13)                      &
     &                       +rcar*export*flui13* bifr13                 ! det in P-units, det13 in C-units (*122)
         ocetra(i,j,k,idet14)=ocetra(i,j,k,idet14)                      &! (kwrbioz loop)
     &                       +rcar*export*flui14* bifr14
#endif /*__c_isotopes*/

#ifdef AGG       
         delsil=MIN(ropal*phosy*avsil/(avsil+bkopal),0.5*avsil) 
       delcar=rcalc*MIN(calmax*phosy,(phosy-delsil/ropal))            ! 
#else
         delsil=MIN(ropal*export*avsil/(avsil+bkopal),0.5*avsil)      
         delcar=rcalc * export * bkopal/(avsil+bkopal)                  ! 'detritus linked calcium carbonate ' ?P units
#endif

! DMS (js: slightly out of place)
         dmsprod = (dmsp5*delsil+dmsp4*delcar)                        &
     &            *(1.+1./(ptho(i,j,k)+dmsp1)**2)         
       dms_bac = dmsp3*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)       &   ! bacterial consumption
     &             *(ocetra(i,j,k,idms)/(dmsp6+ocetra(i,j,k,idms)))       
       dms_uv  = dmsp2*4.*pho*ocetra(i,j,k,idms)                        ! decay due to UV-radiation
       
         ocetra(i,j,k,idms) = ocetra(i,j,k,idms)                      &
     &                          + dmsprod - dms_bac - dms_uv

! end DMS
       
         ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)-delcar +          &    ! - CACO3 production
     &              rcar*( bacfra - phosy + graton + ecan*zoomor)           ! + remineralization C-units

#ifdef __c_isotopes
         ocetra(i,j,k,isco213)=ocetra(i,j,k,isco213)                   &
     &                         -delcar*rocean13                        &
     &                         + bifr13*rcar*( bacfra - phosy + graton + ecan*zoomor) ! bifr=0.98

         ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)-delcar*rocean14
! js: for efficency below line (which should in principle be there) is neglected (additional tracer field would be needed 
!     to account for radioactive decay of 14C in particles)
!     &      + bifr14*rcar*( bacfra - phosy + graton + ecan*zoomor)
#endif /*__c_isotopes*/

         ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)-2.*delcar -       &
     &              rnit*( bacfra - phosy + graton + ecan*zoomor)

         ocetra(i,j,k,iphy)=ocetra(i,j,k,iphy) + phosy - grazing       &
     &                     - phymor - exud  
         ocetra(i,j,k,ioxygen)=ocetra(i,j,k,ioxygen)                   &
     &                        +ro2ut*(phosy-bacfra)                    &
     &                        -(graton+ecan*zoomor)*ro2ut

         ocetra(i,j,k,izoo)=ocetra(i,j,k,izoo)+grawa-excdoc-zoomor

         ocetra(i,j,k,idoc)=ocetra(i,j,k,idoc)-bacfra+excdoc+exud

         ocetra(i,j,k,icalc)=ocetra(i,j,k,icalc)+delcar 

#ifdef __c_isotopes
         ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc13) + delcar *rocean13

         ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14) + delcar *rocean14
#endif

         ocetra(i,j,k,isilica)=ocetra(i,j,k,isilica)-delsil            &
     &                        +dremopal*ocetra(i,j,k,iopal)  

         ocetra(i,j,k,iopal)=ocetra(i,j,k,iopal)+delsil                &
     &                      -dremopal*ocetra(i,j,k,iopal)  

         ocetra(i,j,k,iiron)=ocetra(i,j,k,iiron)                       &
     &          +(bacfra-phosy+graton+ecan*zoomor)*riron               &
     &          - relaxfe*MAX(ocetra(i,j,k,iiron)-fesoly,0.)           


#ifdef AGG
!***********************************************************************
! effects of biological processes on number of particles:
! photosynthesis creates POM
! exudation removes POM
! grazing removes POM; but only the fraction that is not egested as 
! fecal pellets again (grawa remains in zoo, graton goes to po4)     
! none of the processes at the current time is assumed to change
! the size distribution (subject to change)
! NOTE that phosy, exud etc. are in kmol/m3! 
! Thus divide by avmass (kmol/m3)
!**********************************************************************

        if(avmass.gt.0.) then
           avnos = ocetra(i,j,k,inos) 
           anosloss = (phosy-exud-graton-grawa)*avnos/avmass
           ocetra(i,j,k,inos) = ocetra(i,j,k,inos)+anosloss
        endif  

!***********************************************************************
! dead zooplankton corpses come with their own, flat distribution
! this flow even takes place if there is neither nos nor mass
! NOTE: zoomor is in kmol/m3!! Thus multiply flow by 1.e+6
!***********************************************************************

      zmornos = zoomor * (1.-ecan) * zdis * 1.e+6
      ocetra(i,j,k,inos) = ocetra(i,j,k,inos)+zmornos

#endif /*AGG*/

!
! add up over kwrbioz for total inventory (used in inventory_bgc.f90)
!
         expoor(i,j)=expoor(i,j)+pddpo(i,j,k)*export*rcar      
         expoca(i,j)=expoca(i,j)+pddpo(i,j,k)*delcar
         exposi(i,j)=exposi(i,j)+pddpo(i,j,k)*delsil
!
! write output for bgcmean (sum over kwrbioz for 2d fields)
!                          (should that be at the bottom of euphotic zone? (at least for export))
!

         bgcm2d(i,j,jdmsprod)  = bgcm2d(i,j,jdmsprod)             &
     &                                   + dmsprod*pddpo(i,j,k)
         bgcm2d(i,j,jdms_bac)  = bgcm2d(i,j,jdms_bac)             &
     &                                     + dms_bac*pddpo(i,j,k)
         bgcm2d(i,j,jdms_uv)   = bgcm2d(i,j,jdms_uv)              &
     &                                     + dms_uv*pddpo(i,j,k)
         bgcm2d(i,j,jexport)   = bgcm2d(i,j,jexport)              &
     &                                   + export*rcar*pddpo(i,j,k)
         bgcm2d(i,j,jexpoca)   = bgcm2d(i,j,jexpoca)              &
     &                                   + delcar*pddpo(i,j,k)
         bgcm2d(i,j,jexposi)   = bgcm2d(i,j,jexposi)              &
     &                                   + delsil*pddpo(i,j,k)


         bgcm3d(i,j,k,jphosy)  = bgcm3d(i,j,k,jphosy) + phosy 


      ENDIF      ! pddpo(i,j,k).GT.0.5
1     CONTINUE   ! kpie, kpje
100   CONTINUE   ! kwrbioz

      call maschk(kpie,kpje,kpke,2)

#ifdef AGG
       DO  k=1,kpke
         DO j=1,kpje
           DO i=1,kpie
            IF(pddpo(i,j,k).GT.0.5) THEN
             avmass = ocetra(i,j,k,iphy) + ocetra(i,j,k,idet)
             snow = avmass*1.e+6
!  check whether the numbers had to be decreased or increased
             if (snow*pupper.gt.ocetra(i,j,k,inos)) then
               nosin2 = nosin2 + 1
             endif
             if (snow/cellmass.lt.ocetra(i,j,k,inos)) then
               nosde2 = nosde2 + 1
             endif  
            ENDIF
           ENDDO
       ENDDO
      ENDDO
#endif /*AGG*/

      CALL contro(148)


!-----below euphotic zone
      DO 20 k=kwrbioz+1,kpke
         DO 201 j=1,kpje
         DO 201 i=1,kpie
            IF(pddpo(i,j,k).GT.0.5) THEN
#ifdef AGG
            avmass=ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/          
            sterph=dphymor*MAX(0.,ocetra(i,j,k,iphy)-phytomi)                ! 'sterberate' phytoplankton
            sterzo=dzoomor*MAX(0.,ocetra(i,j,k,izoo)-grami)                  ! 'sterberate' zooplankton
                 ocetra(i,j,k,iphy)=ocetra(i,j,k,iphy)-sterph   
            ocetra(i,j,k,izoo)=ocetra(i,j,k,izoo)-sterzo   

            IF(ocetra(i,j,k,ioxygen).gt.5.e-8) THEN                          ! remineralization of poc12 using oxygen
               remin =MIN(drempoc*ocetra(i,j,k,idet),                  &
     &                   0.5*ocetra(i,j,k,ioxygen)/ro2ut)

               detref=remin/(ocetra(i,j,k,idet)+1.e-20)                      ! 'detritus remineralized fraction' (?)
#ifdef __c_isotopes
               rem13=detref*ocetra(i,j,k,idet13)                             ! remineralization of poc13
               rem14=detref*ocetra(i,j,k,idet14)                             !                     poc14
#endif

               docrem=MIN(dremdoc *ocetra(i,j,k,idoc),                 &     ! remineralization of doc
     &                   0.5*(ocetra(i,j,k,ioxygen)-5.e-8)/ro2ut)
            else                   !changed from =max(remin,0.) etc. to =0. js3.5.2006
               remin =0.
#ifdef __c_isotopes
               rem13 =0.
               rem14 =0.
#endif
               docrem=0.
            endif 
          
            ocetra(i,j,k,idet)=ocetra(i,j,k,idet)-remin+sterph+sterzo  
#ifdef __c_isotopes
            ocetra(i,j,k,idet13)=ocetra(i,j,k,idet13)                  &
     &                          +rcar*bifr13*(sterph+sterzo) - rem13
            ocetra(i,j,k,idet14)=ocetra(i,j,k,idet14)-rem14
#endif

            ocetra(i,j,k,ialkali)= ocetra(i,j,k,ialkali)               &
     &                            -rnit*(remin+docrem)

            ocetra(i,j,k,isco212)= ocetra(i,j,k,isco212)               &
     &                            +rcar*(remin+docrem)
#ifdef __c_isotopes
            ocetra(i,j,k,isco213)= ocetra(i,j,k,isco213)               &
     &                            +rcar*docrem*bifr13 + rem13                     ! rem13 not *rcar as idet13 is in C-units
!           ocetra(i,j,k,isco214)= ocetra(i,j,k,isco214)                          ! no biogenic component for 14C
#endif

            ocetra(i,j,k,idoc)=ocetra(i,j,k,idoc)-docrem
            ocetra(i,j,k,ioxygen)= ocetra(i,j,k,ioxygen)               &
     &                            -ro2ut*(remin+docrem)
            ocetra(i,j,k,iphosph)=ocetra(i,j,k,iphosph)+remin+docrem
            ocetra(i,j,k,iano3)=ocetra(i,j,k,iano3)+(remin+docrem)*rnit
            ocetra(i,j,k,iiron)=ocetra(i,j,k,iiron)                    &
     &                       +(remin+docrem)*riron                     &
     &                       -relaxfe*MAX(ocetra(i,j,k,iiron)-fesoly,0.)
!***********************************************************************
! as ragueneau (2000) notes, Si(OH)4sat is about 1000 umol, but
! Si(OH)4 varies only between 0-100 umol
! so the expression dremopal*(Si(OH)4sat-Si(OH)4) would change the 
! rate only from 0 to 100%     
!***********************************************************************
            opalrem=dremopal*ocetra(i,j,k,iopal)
            ocetra(i,j,k,iopal)=ocetra(i,j,k,iopal)-opalrem
            ocetra(i,j,k,isilica)=ocetra(i,j,k,isilica)+opalrem

!***********************************************************************
!           There is about 1.e4 O2 on 1 N2O molecule (Broecker&Peng)
!           refra : Tim Rixen, pers. communication
!***********************************************************************
            aou=satoxy(i,j,k)-ocetra(i,j,k,ioxygen)
            refra=1.+3.*(0.5+sign(0.5,aou-1.97e-4))
            ocetra(i,j,k,ian2o)  =ocetra(i,j,k,ian2o)                  &
     &                                +(remin+docrem)*1.e-4*ro2ut*refra
            ocetra(i,j,k,igasnit)=ocetra(i,j,k,igasnit)                &
     &                                -(remin+docrem)*1.e-4*ro2ut*refra
            ocetra(i,j,k,ioxygen)= ocetra(i,j,k,ioxygen)               &
     &                            -(remin+docrem)*1.e-4*ro2ut*refra*0.5

!careful, pho is very small at large depths (js: why careful? photolysis of DMS will be small, which it should)
!            ocetra(i,j,k,idms)=ocetra(i,j,k,idms)                      &
!     &                        -dmsp2*8.*pho*ocetra(i,j,k,idms)         &
!     &           -dmsp3*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)

       dms_bac = dmsp3*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)       
       
         ocetra(i,j,k,idms) = ocetra(i,j,k,idms)- dms_bac
       

#ifdef AGG
!***********************************************************************
! loss of snow aggregates (by numbers) due to remineralization of poc
! gain of snow aggregates (by numbers) due to zooplankton mortality
! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
!***********************************************************************
           if(avmass.gt.0.) then  
              avnos = ocetra(i,j,k,inos)
              ocetra(i,j,k,inos) = ocetra(i,j,k,inos)                  & 
     &                           - remin * avnos/avmass
           endif
!***********************************************************************
! dead zooplankton corpses come with their own, flat distribution
! this flow even takes place if there is neither nos nor mass
! NOTE: zoomor is in kmol/m3!! Thus multiply flow by 1.e+6
!***********************************************************************
           zmornos = sterzo * zdis * 1.e+6
           ocetra(i,j,k,inos) = ocetra(i,j,k,inos) + zmornos 
#endif /*AGG*/

            ENDIF
 201     CONTINUE
 20   CONTINUE     

!-----below euphotic zone

      call maschk(kpie,kpje,kpke,3)

!---------------------------------------------------------------------
      DO 30 k=kwrbioz+1,kpke
         DO 30 j=1,kpje
         DO 30 i=1,kpie
            IF(ocetra(i,j,k,ioxygen).LT.5.e-7) THEN                            ! denitrification
#ifdef AGG
               avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/          
          
               remin=0.5*drempoc*MIN(ocetra(i,j,k,idet),               &       ! remineralization using NO3
     &                           0.5*ocetra(i,j,k,iano3)/rnit23)

               detref=remin/(ocetra(i,j,k,idet)+1.e-60)                        ! P-units
#ifdef __c_isotopes
               rem13=detref*ocetra(i,j,k,idet13)                               ! C-units
               rem14=detref*ocetra(i,j,k,idet14)                               ! C-units
#endif

               remin2o=dremn2o*MIN(ocetra(i,j,k,idet),                    &    ! remineralization using N2O
     &                              0.003*ocetra(i,j,k,ian2o)/(2*ro2ut))

               detrl=remin2o/(ocetra(i,j,k,idet)+1.e-60)                       ! detrl?
#ifdef __c_isotopes
               rl13=detrl*ocetra(i,j,k,idet13)                                 ! C-units
               rl14=detrl*ocetra(i,j,k,idet14)                                 ! C-units
#endif

               ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)-rnit*(remin + remin2o)

               ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+rcar*(remin + remin2o)
! proxies 13C, 14C 
#ifdef __c_isotopes
               ocetra(i,j,k,isco213)= ocetra(i,j,k,isco213)               &
     &                               + (rem13+rl13  )    
!    &                               +rcar* (rem13+rl13  ) ! changed 3.5.2006
               ocetra(i,j,k,isco214)= ocetra(i,j,k,isco214)               &
     &                               +(rem14+rl14    )
#endif

               ocetra(i,j,k,idet)   =ocetra(i,j,k,idet)   -     (remin + remin2o)
! proxies
#ifdef __c_isotopes
               ocetra(i,j,k,idet13) = ocetra(i,j,k,idet13)                 &
     &                              - (rem13+rl13  )
               ocetra(i,j,k,idet14) = ocetra(i,j,k,idet14)                 &
     &                              - (rem14+rl14    )
#endif

               ocetra(i,j,k,iphosph)=ocetra(i,j,k,iphosph)+     (remin + remin2o)
               ocetra(i,j,k,iano3)  =ocetra(i,j,k,iano3)  -rnit23*remin + rnit*(remin + remin2o)
               ocetra(i,j,k,igasnit)=ocetra(i,j,k,igasnit)+rnit13*remin + 2*ro2ut*remin2o
               ocetra(i,j,k,iiron)  =ocetra(i,j,k,iiron)  +riron*(remin + remin2o)
             ocetra(i,j,k,ian2o)  =ocetra(i,j,k,ian2o)  -2*ro2ut*remin2o  

#ifdef AGG
!***********************************************************************
! loss of snow aggregates (numbers) due to remineralization of poc
! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
!***********************************************************************
           if(avmass.gt.0.) then  
              avnos = ocetra(i,j,k,inos)
              ocetra(i,j,k,inos) = ocetra(i,j,k,inos)         &
     &                - (remin+remin2o)*avnos/avmass
           endif
#endif /*AGG*/

            ENDIF
30    CONTINUE

      call maschk(kpie,kpje,kpke,4)

#ifdef AGG
      DO  k=1,kpke
        DO j=1,kpje
          DO i=1,kpie

           IF(pddpo(i,j,k).GT.0.5) THEN

            avmass = ocetra(i,j,k,iphy) + ocetra(i,j,k,idet)
            snow = avmass*1.e+6

!           check whether the numbers had to be decreased or increased
            if (snow*pupper.gt.ocetra(i,j,k,inos)) then
              nosin3 = nosin3 + 1
            endif

            if (snow/cellmass.lt.ocetra(i,j,k,inos)) then
              nosde3 = nosde3 + 1
            endif

           ENDIF
          ENDDO
        ENDDO
      ENDDO
#endif /*AGG*/

!js why twice ifdef AGG?
#ifdef AGG

! **********************AGGREGATION*****(by Iris Kriest)***************
! General:
! Sinking speed, size distribution and aggregation are calculated 
! as in Kriest and Evans, 2000.
! I assume that opal and calcium carbonate sink at the same speed as P (mass).
!
! Sinking speed and aggregation: I assume that if there is no phosphorous mass,
! the sinking speed is the maximal sinking speed of aggregates. I further
! assume that then there are no particles, and that the rate of aggregation
! is 0. This scheme removes no P in the absence of P, but still opal and/or
! calcium carbonate.
! This could or should be changed, because silica as well as carbonate
! shells will add to the aggregate mass, and should be considered.
! Puh. Does anyone know functional relationships between
! size and Si or CaCO3? Perhaps in a later version, I have to
! take the relationship between mass and size (i.e., density)?
!
! 1. Size distribution and loss of marine snow aggregates due to aggregation 
! (aggregate(i,j,k)) and sinking speed of mass and numbers (wmass(i,j,k)
! and wnumb(i,j,k) are calculated in a loop over 2-kpke. 
!
! 2. The depth of the first layer may change due to ice drift, etc.
! This puts a restriction onto the maximum sinking speed.
! I currently set the max. size for size dependent sinking onto
! one appropriate for this depth.
!
! 3. The fluxes out of the bottom layer are calculated from sinking speed
! and mass concentration, and form the boundary condition for the sediment.
!
! 4. The fluxes in layer kpke->2 are calculated from sinking speed and mass
! concentration, sinking speed, aggregation and number concentration. 

! 5. The fluxes in layer 1 are calculated from sinking speed and mass
! concentration, sinking speed, aggregation and number concentration.  (??)
!************************************************************************

      do k=2,kpke

      do i=1,kpie
      do j=1,kpje
        if(pddpo(i,j,k).gt.0.5) then
          avm = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
          if(avm.gt.0.) then
          snow = avm*1.e+6
          avn = ocetra(i,j,k,inos)
          eps = ((1.+ FractDim)*snow-avn*cellmass) /                   &
     &           (snow-avn*cellmass)

! prevent epsilon from becoming exactly one of the values which are 
! needed for the division 
          if (abs(eps-3.).lt.1.e-15) eps=3.+ vsmall
          if (abs(eps-4.).lt.1.e-15) eps=4.+ vsmall
          if (abs(eps-3.-SinkExp).lt.1.e-15) eps=3.+SinkExp+vsmall
          if (abs(eps-1.-SinkExp-FractDim).lt.1.e-15)                  &
     &        eps=1.+SinkExp+FractDim+vsmall

          e1 = 1. - eps
          e2 = 2. - eps
          e3 = 3. - eps
          e4 = 4. - eps
          es1 = e1 + SinkExp
          es3 = e3 + SinkExp
          TopF = (alar1/alow1)**e1     ! alar1 'largest diameter', alow1 'smallest diameter'
          TopM = TopF*TMFac

! SINKING SPEED FOR THIS LAYER
          wmass(i,j,k) = cellsink * ( (FractDim+e1)/ (FractDim+es1)    &
     &         +TopM*TSFac*SinkExp/ (FractDim+es1))
          wnumb(i,j,k) = cellsink * (e1/es1+TopF*TSFac*SinkExp/es1)

! AGGREGATION

! As a first step, assume that shear in the upper 4 layers is high and 
! zero below. Subject to change. js: this is for 20 layer version.
!                                    include 40 layer version
!                                  should be replaced by check for depth. done 29072005js 
          if (k.le. n90depth) then
            fshear = fsh
          else
            fshear = 0.
          endif     


! shear kernel:
      sagg1 = (TopF-1.)*(TopF*alar3-alow3)*e1/e4                       &
     &   + 3.*(TopF*alar1-alow1)*(TopF*alar2-alow2)*e1*e1/(e2*e3)
      sagg2 = TopF*(                                                   &
     &    (alar3+3.*(alar2*alow1*e1/e2+alar1*alow2*e1/e3)+alow3*e1/e4) &
     &   - TopF*alar3*(1.+3*(       e1/e2+       e1/e3)+     e1/e4))
      sagg4 = TopF*TopF*4.*alar3
      shear_agg = (sagg1+sagg2+sagg4)*fshear

! settlement kernel:
      sagg1 = (TopF * TopF * alar2 * TSFac - alow2)                    &
     &   * SinkExp / (es3 * e3 * (es3 + e1))                           &
     &   + alow2 * ((1. - TopF * TSFac) / (e3 * es1)                   &
     &   - (1. - TopF) / (es3*e1))
      sagg2 = TopF * e1 * (TSFac * ( alow2 - TopF * alar2) / e3        &
     &   - (alow2 - TopF * alar2 * TSFac) / es3)
      sett_agg =  (e1*e1*sagg1+sagg2)*fse
      
      effsti=Stick*(ocetra(i,j,k,iopal)*1.e+6/ropal)/                  &
     &  ((ocetra(i,j,k,iopal)*1.e+6/ropal)+snow)

      aggregate(i,j,k) = (shear_agg+sett_agg)*effsti*avn*avn

! dust aggregation:
! shear kernel:
      dfirst=dustd3+3.*dustd2*alar1+3.*dustd1*alar2+alar3                ! dustd3: dust diameter**3, d2: **2
      dshagg=e1*fsh*(dfirst*TopF/e1-(                                  &
     &  (TopF-1.)/e1*dustd3+3.*(TopF*alar1-alow1)/e2*dustd2            &
     &   +3.*(TopF*alar2-alow2)/e3*dustd1+(TopF*alar3-alow3)/e4))

! settlement kernel:
      dsett=fse*dustd2*((e1+SinkExp*TopF*TSFac)/es1-dustsink/cellsink)
      
      dustagg(i,j,k) = effsti*avn*ocetra(i,j,k,ifdust)                 &
     &                *(dshagg+dsett)

      else    ! available mass le 0
        wmass(i,j,k)=TSFac*cellsink
        wnumb(i,j,k)=0.
        aggregate(i,j,k)=0.
        dustagg(i,j,k)=0.
        ocetra(i,j,k,inos)=0.
      endif

      endif   ! wet cell

      enddo   ! je
      enddo   ! ie
      enddo   ! ke

      call maschk(kpie,kpje,kpke,5)

! EVALUATE SINKING RATE AND AGGREGATION FOR SURFACE LAYER, WHICH MAY BE
! LESS DEEP THAN INITIALLY SET BECAUSE OF EVAPORATION, ICE ETC.

      DO j=1,kpje
      DO i=1,kpie
         if(pddpo(i,j,1).gt.0.5) then

!ik evaluate safe length scale for size dependent sinking and
!ik aggregation, and the resulting sinking rate and aggregation rate.
!ik zo may reduce the first layer depth to values that are small and
!ik may cause the sinking length to exceed the layers depth.
!ik to be safe, for this upper layer set the upper size such that
!ik loss due to sinking is at max the whole inventory of this box.
!ik aggregation will be calculated accordingly.

         checksize = (pddpo(i,j,1)/cellsink)**(1./SinkExp)*alow1
         if(alar1.gt.checksize) then
           nacheck=nacheck+1              ! js: seems not to be used
         endif
         flar1 = MIN(alar1,checksize)     ! reduce diameter of largest particle
         flar2 = flar1 * flar1
         flar3 = flar2 * flar1
         fTSFac = (flar1/alow1)**SinkExp
         fTMFac = (flar1/alow1)**FractDim

! SIZE DITRIBUTION
         avm = ocetra(i,j,1,iphy)+ocetra(i,j,1,idet)   ! available mass   (js: add dust here to account for ballast effect?)
         if(avm.gt.0.) then
         snow = avm*1.e+6
         avn = ocetra(i,j,1,inos)                      ! available numbers
         eps = ((1.+ FractDim)*snow-avn*cellmass) /                    &  ! exponential coefficient of size distribution
     &           (snow-avn*cellmass)
     
         if (abs(eps-3.).lt.1.e-15) eps=3.+ vsmall
         if (abs(eps-4.).lt.1.e-15) eps=4.+ vsmall
         if (abs(eps-3.-SinkExp).lt.1.e-15) eps=3.+SinkExp+vsmall
         if (abs(eps-1.-SinkExp-FractDim).lt.1.e-15)                   &
     &        eps=1.+SinkExp+FractDim+vsmall

         e1 = 1. - eps
         e2 = 2. - eps
         e3 = 3. - eps
         e4 = 4. - eps
         es1 = e1 + SinkExp
         es3 = e3 + SinkExp

         fTopF = (flar1/alow1)**e1
         fTopM = fTopF*fTMFac

! SINKING SPEEDS
         wmass(i,j,1) = cellsink * ( (FractDim+e1)/ (FractDim+es1)     &
     &          +fTopM*fTSFac*SinkExp/ (FractDim+es1))
         wnumb(i,j,1) = cellsink * (e1/es1+fTopF*fTSFac*SinkExp/es1)

! AGGREGATION
      sagg1 = (fTopF-1.)*(fTopF*flar3-alow3)*e1/e4                     &
     &   + 3.*(fTopF*flar1-alow1)*(fTopF*flar2-alow2)*e1*e1/(e2*e3)
      sagg2 = fTopF*(                                                  &
     &    (flar3+3.*(flar2*alow1*e1/e2+flar1*alow2*e1/e3)+alow3*e1/e4) &
     &   - fTopF*flar3*(1.+3*(       e1/e2+       e1/e3)+     e1/e4))
      sagg4 = fTopF*fTopF*4.*flar3
      shear_agg = (sagg1+sagg2+sagg4)*fsh

      sagg1 = (fTopF * fTopF * flar2 * fTSFac - alow2)                 &
     &   * SinkExp / (es3 * e3 * (es3 + e1))                           &
     &   + alow2 * ((1. - fTopF * fTSFac) / (e3 * es1)                 &
     &   - (1. - fTopF) / (es3*e1))
      sagg2 = fTopF * e1 * (fTSFac * ( alow2 - fTopF * flar2) / e3     &
     &   - (alow2 - fTopF * flar2 * fTSFac) / es3)
      sett_agg =  (e1*e1*sagg1+sagg2)*fse


      effsti=Stick*(ocetra(i,j,1,iopal)*1.e+6/ropal)/                  &
     &  ((ocetra(i,j,1,iopal)*1.e+6/ropal)+snow)

      aggregate(i,j,1) = (shear_agg+sett_agg)*effsti*avn*avn

! dust aggregation:
! shear kernel:
      dfirst=dustd3+3.*dustd2*flar1+3.*dustd1*flar2+flar3
      dshagg=e1*fsh*(dfirst*fTopF/e1-(                                 &
     &  (fTopF-1.)/e1*dustd3+3.*(fTopF*flar1-alow1)/e2*dustd2          &
     &   +3.*(fTopF*flar2-alow2)/e3*dustd1+(fTopF*flar3-alow3)/e4))

! settlement kernel:
      dsett=fse*dustd2*((e1+SinkExp*fTopF*fTSFac)/es1-dustsink/cellsink)
      
      dustagg(i,j,1) = effsti*avn*ocetra(i,j,1,ifdust)                 &
     &                *(dshagg+dsett)

      else                            ! available mass le 0.
        wmass(i,j,1)=fTSFac*cellsink
        wnumb(i,j,1)=0.
        aggregate(i,j,1)=0.
        dustagg(i,j,1)=0.
        ocetra(i,j,1,inos)=0.
      endif

      endif    ! wet cell

      enddo
      enddo
      
! EVALUATE SINKING RATE AND AGGREGATION FOR BOTTOM LAYER, WHICH MAY BE
! LESS THICK THAN THE MINIMUM LAYER THICKNESS

      DO j=1,kpje
      DO i=1,kpie
         if(pddpo(i,j,1).gt.0.5) then
         if(alar1max(i,j).lt.alar1) then

!ik take safe length scale for size dependent sinking and
!ik aggregation, and the resulting sinking rate and aggregation rate.

         flar1 = alar1max(i,j)
         flar2 = flar1 * flar1
         flar3 = flar2 * flar1
         fTSFac = TSFmax(i,j)
         fTMFac = TMFmax(i,j)

! SIZE DITRIBUTION
         avm = ocetra(i,j,kbo(i,j),iphy)+ocetra(i,j,kbo(i,j),idet)
         if(avm.gt.0.) then
         snow = avm*1.e+6       ! why *1.e6?
         avn = ocetra(i,j,kbo(i,j),inos)
         eps = ((1.+ FractDim)*snow-avn*cellmass) /                    &
     &           (snow-avn*cellmass)
     
         if (abs(eps-3.).lt.1.e-15) eps=3.+ vsmall
         if (abs(eps-4.).lt.1.e-15) eps=4.+ vsmall
         if (abs(eps-3.-SinkExp).lt.1.e-15) eps=3.+SinkExp+vsmall
         if (abs(eps-1.-SinkExp-FractDim).lt.1.e-15)                   &
     &        eps=1.+SinkExp+FractDim+vsmall

         e1 = 1. - eps
         e2 = 2. - eps
         e3 = 3. - eps
         e4 = 4. - eps
         es1 = e1 + SinkExp
         es3 = e3 + SinkExp

         fTopF = (flar1/alow1)**e1
         fTopM = fTopF*fTMFac

! SINKING SPEEDS
         wmass(i,j,kbo(i,j)) = cellsink *                              &
     &        ( (FractDim+e1)/ (FractDim+es1)                          &
     &          +fTopM*fTSFac*SinkExp/ (FractDim+es1))
         wnumb(i,j,kbo(i,j)) = cellsink *                              &
     &          (e1/es1+fTopF*fTSFac*SinkExp/es1)

! AGGREGATION
      sagg1 = (fTopF-1.)*(fTopF*flar3-alow3)*e1/e4                     &
     &   + 3.*(fTopF*flar1-alow1)*(fTopF*flar2-alow2)*e1*e1/(e2*e3)
      sagg2 = fTopF*(                                                  &
     &    (flar3+3.*(flar2*alow1*e1/e2+flar1*alow2*e1/e3)+alow3*e1/e4) &
     &   - fTopF*flar3*(1.+3*(       e1/e2+       e1/e3)+     e1/e4))
      sagg4 = fTopF*fTopF*4.*flar3
      shear_agg = (sagg1+sagg2+sagg4)*fsh

      sagg1 = (fTopF * fTopF * flar2 * fTSFac - alow2)                 &
     &   * SinkExp / (es3 * e3 * (es3 + e1))                           &
     &   + alow2 * ((1. - fTopF * fTSFac) / (e3 * es1)                 &
     &   - (1. - fTopF) / (es3*e1))
      sagg2 = fTopF * e1 * (fTSFac * ( alow2 - fTopF * flar2) / e3     &
     &   - (alow2 - fTopF * flar2 * fTSFac) / es3)
      sett_agg =  (e1*e1*sagg1+sagg2)*fse


      effsti=Stick*(ocetra(i,j,kbo(i,j),iopal)*1.e+6/ropal)/           &
     &  ((ocetra(i,j,kbo(i,j),iopal)*1.e+6/ropal)+snow)

      aggregate(i,j,kbo(i,j)) = (shear_agg+sett_agg)*effsti*avn*avn

! dust aggregation:
! shear kernel:
      dfirst=dustd3+3.*dustd2*flar1+3.*dustd1*flar2+flar3
      dshagg=e1*fsh*(dfirst*fTopF/e1-(                                 &
     &  (fTopF-1.)/e1*dustd3+3.*(fTopF*flar1-alow1)/e2*dustd2          &
     &   +3.*(fTopF*flar2-alow2)/e3*dustd1+(fTopF*flar3-alow3)/e4))

! settlement kernel:
      dsett=fse*dustd2*((e1+SinkExp*fTopF*fTSFac)/es1-dustsink/cellsink)
      
      dustagg(i,j,kbo(i,j)) = effsti*avn*ocetra(i,j,kbo(i,j),ifdust)   &
     &                *(dshagg+dsett)

      else
        wmass(i,j,kbo(i,j))=fTSFac*cellsink
        wnumb(i,j,kbo(i,j))=0.
        aggregate(i,j,kbo(i,j))=0.
        dustagg(i,j,kbo(i,j))=0.
        ocetra(i,j,kbo(i,j),inos)=0.
      endif ! avm

      endif ! alar1max

      endif ! pddpo

      enddo
      enddo

      call maschk(kpie,kpje,kpke,6)
      
!
! Sampling timeseries-1 : sedimentation at specific positions (no global value!)
!                         particle flux
!
      DO l=1,nts
         i = its1(l)-p_ioff
         j = jts1(l)-p_joff
         IF(i<=1 .OR. i>=kpie .OR. j<=1 .OR. j>=kpje) CYCLE

! first depth 

         if(k1ts1(l).gt.0) then
             wphy = wmass(i,j,k1ts1(l))
           
             fphy=wphy*ocetra(i,j,k1ts1(l),iphy)  
             fopa=wphy*ocetra(i,j,k1ts1(l),iopal)  
             fdet=wphy*ocetra(i,j,k1ts1(l),idet)  
             fcal=wphy*ocetra(i,j,k1ts1(l),icalc) 
             ts1(its1fdet,l+1,lts1) = ts1(its1fdet,l+1,lts1) + fphy + fdet
             ts1(its1fopa,l+1,lts1) = ts1(its1fopa,l+1,lts1) + fopa
             ts1(its1fcal,l+1,lts1) = ts1(its1fcal,l+1,lts1) + fcal
          else
             ts1(its1fdet,l+1,lts1) = ts1(its1fdet,l+1,lts1) -9999.
             ts1(its1fopa,l+1,lts1) = ts1(its1fopa,l+1,lts1) -9999.
             ts1(its1fcal,l+1,lts1) = ts1(its1fcal,l+1,lts1) -9999.
          endif          

! second depth 

         if(k2ts1(l).gt.0) then
             wphy = wmass(i,j,k2ts1(l))
           
             fphy=wphy*ocetra(i,j,k2ts1(l),iphy)  
             fopa=wphy*ocetra(i,j,k2ts1(l),iopal)  
             fdet=wphy*ocetra(i,j,k2ts1(l),idet)  
             fcal=wphy*ocetra(i,j,k2ts1(l),icalc) 
             ts1(its2fdet,l+1,lts1) = ts1(its2fdet,l+1,lts1) + fphy + fdet
             ts1(its2fopa,l+1,lts1) = ts1(its2fopa,l+1,lts1) + fopa
             ts1(its2fcal,l+1,lts1) = ts1(its2fcal,l+1,lts1) + fcal
          else
             ts1(its2fdet,l+1,lts1) = ts1(its2fdet,l+1,lts1) -9999.
             ts1(its2fopa,l+1,lts1) = ts1(its2fopa,l+1,lts1) -9999.
             ts1(its2fcal,l+1,lts1) = ts1(its2fcal,l+1,lts1) -9999.
          endif
! third depth 
         if(k3ts1(l).gt.0) then
             wphy = wmass(i,j,k3ts1(l))
           
             fphy=wphy*ocetra(i,j,k3ts1(l),iphy)  
             fopa=wphy*ocetra(i,j,k3ts1(l),iopal)  
             fdet=wphy*ocetra(i,j,k3ts1(l),idet)  
             fcal=wphy*ocetra(i,j,k2ts1(l),icalc) 
             ts1(its3fdet,l+1,lts1) = ts1(its3fdet,l+1,lts1) + fphy + fdet
             ts1(its3fopa,l+1,lts1) = ts1(its3fopa,l+1,lts1) + fopa
             ts1(its3fcal,l+1,lts1) = ts1(its3fcal,l+1,lts1) + fcal
          else
             ts1(its3fdet,l+1,lts1) = ts1(its3fdet,l+1,lts1) -9999.
             ts1(its3fopa,l+1,lts1) = ts1(its3fopa,l+1,lts1) -9999.
             ts1(its3fcal,l+1,lts1) = ts1(its3fcal,l+1,lts1) -9999.
         endif
      ENDDO
!
! prepare output for bgcmean files
!
!js: here the check should be depth >90 m etc (13.09.2005)
!    (the original approach worked here (ifdef AGG) because wmass
!    is zero for depth <kbo. it does not work below (ifndef AGG)
      DO j=1,kpje
      DO i=1,kpie
!js      if(pddpo(i,j,kbo(i,j)).gt.0.5) then

! fluxes at 90 m (aggregation part)
         if (kbo(i,j).gt.n90depth) then
         bgcm2d(i,j,jcoex90) = bgcm2d(i,j,jcoex90)              &
     &   +(ocetra(i,j,n90depth,iphy)                            &
     &   +ocetra(i,j,n90depth,idet))*wmass(i,j,n90depth)
         bgcm2d(i,j,jopex90) = bgcm2d(i,j,jopex90) +            &
     &   ocetra(i,j,n90depth,iopal)*wmass(i,j,n90depth)
         bgcm2d(i,j,jcaex90) = bgcm2d(i,j,jcaex90) +            &
     &   ocetra(i,j,n90depth,icalc)*wmass(i,j,n90depth)

! move masking from write_bgcmean_2d to here?
!        elseif
!         bgcm2d(i,j,jcoex90) = rmasko etc...

         
         endif
! fluxes at about 1000 m
         if (kbo(i,j).gt.n1000depth) then
         bgcm2d(i,j,jcoex1000) = bgcm2d(i,j,jcoex1000) +        &
     &   (ocetra(i,j,n1000depth,iphy)                           &
     &   +ocetra(i,j,n1000depth,idet))*wmass(i,j,n1000depth)
         bgcm2d(i,j,jopex1000) = bgcm2d(i,j,jopex1000) +        &
     &   ocetra(i,j,n1000depth,iopal)*wmass(i,j,n1000depth)
         bgcm2d(i,j,jcaex1000) = bgcm2d(i,j,jcaex1000) +        &
     &   ocetra(i,j,n1000depth,icalc)*wmass(i,j,n1000depth)
         endif
! fluxes at about 1950 m
         if (kbo(i,j).gt.n2000depth) then
         bgcm2d(i,j,jcoex2000) = bgcm2d(i,j,jcoex2000) +        &
     &   (ocetra(i,j,n2000depth,iphy)                           &
     &   +ocetra(i,j,n2000depth,idet))*wmass(i,j,n2000depth)
         bgcm2d(i,j,jopex2000) = bgcm2d(i,j,jopex2000) +        &
     &   ocetra(i,j,n2000depth,iopal)*wmass(i,j,n2000depth)
         bgcm2d(i,j,jcaex2000) = bgcm2d(i,j,jcaex2000) +        &
     &   ocetra(i,j,n2000depth,icalc)*wmass(i,j,n2000depth)
         endif

!js      endif
      ENDDO
      ENDDO

!IK COMPUTE FLUXES FOR BOUNDARY CONDITION/BOTTOM LAYER

!js fluxes to sediment   (still AGG)
      
      DO 36 j=1,kpje
      DO 36 i=1,kpie
         if(pddpo(i,j,kbo(i,j)).gt.0.5) then
         wphy = wmass(i,j,kbo(i,j))

         prorca(i,j) = ocetra(i,j,kbo(i,j),iphy)  *wphy                &
     &               + ocetra(i,j,kbo(i,j),idet)  *wphy
         prcaca(i,j) = ocetra(i,j,kbo(i,j),icalc) *wphy
         silpro(i,j) = ocetra(i,j,kbo(i,j),iopal) *wphy
         produs(i,j) = ocetra(i,j,kbo(i,j),ifdust)*dustsink            &
     &               + ocetra(i,j,kbo(i,j),iadust)*wphy    

!
! prepare output for bgcmean files 
!
         bgct2d(i,j,jprorca)=bgct2d(i,j,jprorca) +prorca(i,j)
         bgct2d(i,j,jprcaca)=bgct2d(i,j,jprcaca) +prcaca(i,j)
         bgct2d(i,j,jsilpro)=bgct2d(i,j,jsilpro) +silpro(i,j)
         bgct2d(i,j,jprodus)=bgct2d(i,j,jprodus) +produs(i,j)


        endif
36      CONTINUE

! COMPUTE FLUXES FOR LAYERS 2 TO kpke
      DO 2 K=kpke,2,-1
      DO 34 j=1,kpje
      DO 34 i=1,kpie
         if(pddpo(i,j,k).gt.0.5) then

! SINKING SPEED FOR UPPER LAYER
          wphyup = wmass(i,j,k-1)    ! settling velocity of mass
          wnosup = wnumb(i,j,k-1)    ! settling velocity of number of marine snow aggregates

! SINKING SPEED FOR ACTUAL LAYER
          wphy = wmass(i,j,k)
          wnos = wnumb(i,j,k)

! SUM-UP FLUXES (compute new concentrations)
        ocetra(i,j,k,iphy) =ocetra(i,j,k,iphy) +                       &
     &      (ocetra(i,j,k-1,iphy)*wphyup-ocetra(i,j,k,iphy)*wphy)      &
     &       *pdpio(i,j,k)
        ocetra(i,j,k,idet) =ocetra(i,j,k,idet) +                       &
     &      (ocetra(i,j,k-1,idet)*wphyup-ocetra(i,j,k,idet)*wphy)      &
     &       *pdpio(i,j,k)
#ifdef __c_isotopes
        ocetra(i,j,k,idet13) =ocetra(i,j,k,idet13) +                   &
     &      (ocetra(i,j,k-1,idet13)*wphyup-ocetra(i,j,k,idet13)*wphy)  &
     &       *pdpio(i,j,k)
        ocetra(i,j,k,idet14) =ocetra(i,j,k,idet14) +                   &
     &      (ocetra(i,j,k-1,idet14)*wphyup-ocetra(i,j,k,idet14)*wphy)  &
     &       *pdpio(i,j,k)

        ocetra(i,j,k,icalc13) =ocetra(i,j,k,icalc13) +                 &
     &      (ocetra(i,j,k-1,icalc13)*wphyup-ocetra(i,j,k,icalc13)*wphy)&
     &       *pdpio(i,j,k)
        ocetra(i,j,k,icalc14) =ocetra(i,j,k,icalc14) +                 &
     &      (ocetra(i,j,k-1,icalc14)*wphyup-ocetra(i,j,k,icalc14)*wphy)&
     &       *pdpio(i,j,k)
#endif

        ocetra(i,j,k,icalc) =ocetra(i,j,k,icalc) +                     &
     &      (ocetra(i,j,k-1,icalc)*wphyup-ocetra(i,j,k,icalc)*wphy)    &
     &       *pdpio(i,j,k)
        ocetra(i,j,k,iopal) =ocetra(i,j,k,iopal) +                     &
     &      (ocetra(i,j,k-1,iopal)*wphyup-ocetra(i,j,k,iopal)*wphy)    &
     &       *pdpio(i,j,k)
        ocetra(i,j,k,inos) =ocetra(i,j,k,inos) - aggregate(i,j,k) +    &
     &      (ocetra(i,j,k-1,inos)*wnosup-ocetra(i,j,k,inos)*wnos)      &
     &       *pdpio(i,j,k)
! sinking of free dust and loss due to attachment to aggregated dust
        ocetra(i,j,k,ifdust) =ocetra(i,j,k,ifdust) - dustagg(i,j,k) +  &
     &      (ocetra(i,j,k-1,ifdust)-ocetra(i,j,k,ifdust))*dustsink     &
     &       *pdpio(i,j,k)
! sinking of aggregated dust and gain due to attachment of free dust to aggregates
        ocetra(i,j,k,iadust) =ocetra(i,j,k,iadust) + dustagg(i,j,k) +  &
     &      (ocetra(i,j,k-1,iadust)*wphyup-ocetra(i,j,k,iadust)*wphy)  &
     &       *pdpio(i,j,k)
      endif

   34 CONTINUE   ! end i,j-loop

    2 CONTINUE   ! end k-loop



!IK  COMPUTE FLUXES FOR SURFACE LAYER

      DO 35 j=1,kpje
      DO 35 i=1,kpie
         if(pddpo(i,j,1).gt.0.) then

         wphy = wmass(i,j,1)
         wnos = wnumb(i,j,1)

! SUM-UP FLUXES
        ocetra(i,j,1,iphy) =ocetra(i,j,1,iphy)                          &
     &          -ocetra(i,j,1,iphy)*wphy*pdpio(i,j,1) 
        ocetra(i,j,1,idet) =ocetra(i,j,1,idet)                          &
     &          -ocetra(i,j,1,idet)*wphy*pdpio(i,j,1)
        ocetra(i,j,1,icalc) =ocetra(i,j,1,icalc)                        &
     &          -ocetra(i,j,1,icalc)*wphy*pdpio(i,j,1)
        ocetra(i,j,1,iopal) =ocetra(i,j,1,iopal)                        &
     &         -ocetra(i,j,1,iopal)*wphy*pdpio(i,j,1)
        ocetra(i,j,1,inos) =ocetra(i,j,1,inos) - aggregate(i,j,1)       &
     &         -ocetra(i,j,1,inos)*wnos*pdpio(i,j,1)

! sinking of free dust and loss of free dust to aggregated dust
        ocetra(i,j,1,ifdust) =ocetra(i,j,1,ifdust) - dustagg(i,j,1)     &
     &         -ocetra(i,j,1,ifdust)*dustsink*pdpio(i,j,1)

! sinking of aggregated dust and gain due to attachment of free dust to aggregates
        ocetra(i,j,1,iadust) =ocetra(i,j,1,iadust) + dustagg(i,j,1)     &
     &         -ocetra(i,j,1,iadust)*wphy*pdpio(i,j,1)

          endif
35    CONTINUE 

      call maschk(kpie,kpje,kpke,8)


#endif /*AGG*/
!-------------------------------------------end aggregation part

#ifndef AGG
!
! Sampling timeseries-1 : 'iris'-sedimentation at specific positions (i.e., downward flux of matter)
!

      call maschk(kpie,kpje,kpke,-7)

      DO l=1,nts
         i = its1(l)-p_ioff
         j = jts1(l)-p_joff
         IF(i<=1 .OR. i>=kpie .OR. j<=1 .OR. j>=kpje) CYCLE

! first sampling depth of station in time series

         if(k1ts1(l).gt.0) then             
             fopa=wopal*ocetra(i,j,k1ts1(l),iopal)  
             fdet=wpoc *ocetra(i,j,k1ts1(l),idet)  
             fcal=wcal *ocetra(i,j,k1ts1(l),icalc) 
             ts1(its1fdet,l+1,lts1) = ts1(its1fdet,l+1,lts1) + fdet
             ts1(its1fopa,l+1,lts1) = ts1(its1fopa,l+1,lts1) + fopa
             ts1(its1fcal,l+1,lts1) = ts1(its1fcal,l+1,lts1) + fcal
          else
             ts1(its1fdet,l+1,lts1) = ts1(its1fdet,l+1,lts1) -9999.
             ts1(its1fopa,l+1,lts1) = ts1(its1fopa,l+1,lts1) -9999.
             ts1(its1fcal,l+1,lts1) = ts1(its1fcal,l+1,lts1) -9999.
          endif          

! second depth 

         if(k2ts1(l).gt.0) then
             fopa=wopal*ocetra(i,j,k2ts1(l),iopal)  
             fdet=wpoc *ocetra(i,j,k2ts1(l),idet)  
             fcal=wcal *ocetra(i,j,k2ts1(l),icalc) 
             ts1(its2fdet,l+1,lts1) = ts1(its2fdet,l+1,lts1) + fdet
             ts1(its2fopa,l+1,lts1) = ts1(its2fopa,l+1,lts1) + fopa
             ts1(its2fcal,l+1,lts1) = ts1(its2fcal,l+1,lts1) + fcal
          else
             ts1(its2fdet,l+1,lts1) = ts1(its2fdet,l+1,lts1) -9999.
             ts1(its2fopa,l+1,lts1) = ts1(its2fopa,l+1,lts1) -9999.
             ts1(its2fcal,l+1,lts1) = ts1(its2fcal,l+1,lts1) -9999.
          endif  

! third depth 

         if(k3ts1(l).gt.0) then
             fopa=wopal*ocetra(i,j,k3ts1(l),iopal)  
             fdet=wpoc *ocetra(i,j,k3ts1(l),idet)  
             fcal=wcal *ocetra(i,j,k3ts1(l),icalc) 
             ts1(its3fdet,l+1,lts1) = ts1(its3fdet,l+1,lts1) + fdet
             ts1(its3fopa,l+1,lts1) = ts1(its3fopa,l+1,lts1) + fopa
             ts1(its3fcal,l+1,lts1) = ts1(its3fcal,l+1,lts1) + fcal
          else
             ts1(its3fdet,l+1,lts1) = ts1(its3fdet,l+1,lts1) -9999.
             ts1(its3fopa,l+1,lts1) = ts1(its3fopa,l+1,lts1) -9999.
             ts1(its3fcal,l+1,lts1) = ts1(its3fcal,l+1,lts1) -9999.
         endif
      ENDDO
!
! write output for bgcmean       
!
!js  check against 90, 1000, 2000 m otherwise output is corrupted
!    were water depth is >0 but less than export depth
      DO j=1,kpje
      DO i=1,kpie
!js      if(pddpo(i,j,kbo(i,j)).gt.0.5) then
! write out fluxes at 90 m (fixed settling velocity part)
         if (kbo(i,j).gt.n90depth) then
         bgcm2d(i,j,jcoex90) = bgcm2d(i,j,jcoex90) +     &
     &   ocetra(i,j,n90depth,idet)*wpoc
         bgcm2d(i,j,jopex90) = bgcm2d(i,j,jopex90) +     &
     &   ocetra(i,j,n90depth,iopal)*wopal
         bgcm2d(i,j,jcaex90) = bgcm2d(i,j,jcaex90) +     &
     &   ocetra(i,j,n90depth,icalc)*wcal
         endif
! write out fluxes at about 1000 m
         if (kbo(i,j).gt.n1000depth) then
         bgcm2d(i,j,jcoex1000) = bgcm2d(i,j,jcoex1000) + &
     &   ocetra(i,j,n1000depth,idet)*wpoc
         bgcm2d(i,j,jopex1000) = bgcm2d(i,j,jopex1000) + &
     &   ocetra(i,j,n1000depth,iopal)*wopal
         bgcm2d(i,j,jcaex1000) = bgcm2d(i,j,jcaex1000) + &
     &   ocetra(i,j,n1000depth,icalc)*wcal
         endif
! write out fluxes at about 1950 m
         if (kbo(i,j).gt.n2000depth) then
         bgcm2d(i,j,jcoex2000) = bgcm2d(i,j,jcoex2000) + &
     &   ocetra(i,j,n2000depth,idet)*wpoc
         bgcm2d(i,j,jopex2000) = bgcm2d(i,j,jopex2000) + &
     &   ocetra(i,j,n2000depth,iopal)*wopal
         bgcm2d(i,j,jcaex2000) = bgcm2d(i,j,jcaex2000) + &
     &   ocetra(i,j,n2000depth,icalc)*wcal
         endif

!js      endif
      ENDDO
      ENDDO
! 23  -> 15  (k=, now n1000depth)
! 29  -> 18  (now n2000depth)
!

! implicit method:
! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
! -->           
! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
! sedimentation=w*dt*C(ks,T+dt)
!

      call maschk(kpie,kpje,kpke,-8)

       k=1                         ! -----------surface layer
       DO j=1,kpje
       DO i=1,kpie
          IF(pddpo(i,j,k).GT.0.5) THEN

             ocetra(i,j,k,idet) =(ocetra(i,j,k,idet)*pddpo(i,j,k))/    &
     &                           (pddpo(i,j,k)+wpoc)
!proxies
#ifdef __c_isotopes
         ocetra(i,j,k,idet13) =(ocetra(i,j,k,idet13)*pddpo(i,j,k))/    &
     &                           (pddpo(i,j,k)+wpoc)
         ocetra(i,j,k,idet14) =(ocetra(i,j,k,idet14)*pddpo(i,j,k))/    &
     &                           (pddpo(i,j,k)+wpoc)
#endif

             ocetra(i,j,k,icalc)=(ocetra(i,j,k,icalc)*pddpo(i,j,k))/   &
     &                           (pddpo(i,j,k)+wcal)

!proxies
#ifdef __c_isotopes
         ocetra(i,j,k,icalc13)=(ocetra(i,j,k,icalc13)*pddpo(i,j,k))/   &
     &                           (pddpo(i,j,k)+wcal)
         ocetra(i,j,k,icalc14)=(ocetra(i,j,k,icalc14)*pddpo(i,j,k))/   &
     &                           (pddpo(i,j,k)+wcal)
#endif
 
             ocetra(i,j,k,iopal)=(ocetra(i,j,k,iopal)*pddpo(i,j,k))/   &
     &                           (pddpo(i,j,k)+wopal)    
 
             ocetra(i,j,k,ifdust)=(ocetra(i,j,k,ifdust)*pddpo(i,j,k))/ &
     &                            (pddpo(i,j,k)+wdust)    
     
          ENDIF
      enddo
      enddo
    
      DO 10 k=2,kpke                    ! ---------------------- water column
        DO 12 j=1,kpje
        DO 12 i=1,kpie
          IF(pddpo(i,j,k).GT.0.5) THEN

             ocetra(i,j,k,idet)=(ocetra(i,j,k  ,idet)*pddpo(i,j,k)     &
     &                              +ocetra(i,j,k-1,idet)*wpoc)/           &
     &                          (pddpo(i,j,k)+wpoc)
#ifdef __c_isotopes
             ocetra(i,j,k,idet13)=(ocetra(i,j,k  ,idet13)*pddpo(i,j,k) &
     &                          +ocetra(i,j,k-1,idet13)*wpoc)/         &
     &                          (pddpo(i,j,k)+wpoc)
             ocetra(i,j,k,idet14)=(ocetra(i,j,k  ,idet14)*pddpo(i,j,k) &
     &                          +ocetra(i,j,k-1,idet14)*wpoc)/         &
     &                          (pddpo(i,j,k)+wpoc)
#endif

             ocetra(i,j,k,icalc)=(ocetra(i,j,k  ,icalc)*pddpo(i,j,k)   &
     &                               +ocetra(i,j,k-1,icalc)*wcal)/         &
     &                           (pddpo(i,j,k)+wcal)
#ifdef __c_isotopes
         ocetra(i,j,k,icalc13)=(ocetra(i,j,k  ,icalc13)*pddpo(i,j,k)   &
     &                         +ocetra(i,j,k-1,icalc13)*wcal)/         &
     &                           (pddpo(i,j,k)+wcal)
         ocetra(i,j,k,icalc14)=(ocetra(i,j,k  ,icalc14)*pddpo(i,j,k)   &
     &                         +ocetra(i,j,k-1,icalc14)*wcal)/         &
     &                         (pddpo(i,j,k)+wcal)
#endif
 
             ocetra(i,j,k,iopal)=(ocetra(i,j,k  ,iopal)*pddpo(i,j,k)   &
     &                               +ocetra(i,j,k-1,iopal)*wopal)/        &
     &                           (pddpo(i,j,k)+wopal)        
 
             ocetra(i,j,k,ifdust)=(ocetra(i,j,k ,ifdust)*pddpo(i,j,k)  &
     &                                +ocetra(i,j,k-1,ifdust)*wdust)/      &
     &                            (pddpo(i,j,k)+wdust)        

          ENDIF
12      CONTINUE
10    CONTINUE

      call maschk(kpie,kpje,kpke,9)

! ------------------------------flux to sediment
      DO 33 j=1,kpje
      DO 33 i=1,kpie
       IF(pddpo(i,j,1).GT.0.5) THEN
         prorca(i,j)=ocetra(i,j,kbo(i,j),idet )*wpoc
         prcaca(i,j)=ocetra(i,j,kbo(i,j),icalc)*wcal
#ifdef __c_isotopes
         pror13(i,j)=ocetra(i,j,kbo(i,j),idet13 )*wpoc
         prca13(i,j)=ocetra(i,j,kbo(i,j),icalc13)*wcal
         pror14(i,j)=ocetra(i,j,kbo(i,j),idet14 )*wpoc
         prca14(i,j)=ocetra(i,j,kbo(i,j),icalc14)*wcal
#endif
         silpro(i,j)=ocetra(i,j,kbo(i,j),iopal)*wopal
         produs(i,j)=ocetra(i,j,kbo(i,j),ifdust)*wdust

!
! Paddy: write output for bgcmean       
!
         bgct2d(i,j,jprorca)=bgct2d(i,j,jprorca) + prorca(i,j)
         bgct2d(i,j,jprcaca)=bgct2d(i,j,jprcaca) + prcaca(i,j)
         bgct2d(i,j,jsilpro)=bgct2d(i,j,jsilpro) + silpro(i,j)
         bgct2d(i,j,jprodus)=bgct2d(i,j,jprodus) + produs(i,j)

       ENDIF
33    CONTINUE

      call maschk(kpie,kpje,kpke,10)

#endif /* not AGG*/



!
! Check maximum/minimum values in wet/dry cells.
!
      IF( kchck .EQ. 1) CALL CHCK_BGC(io_stdo_bgc,icyclibgc,           &
     &'Check values of ocean tracer at exit from SBR OCPROD :',        &
     & kpie,kpje,kpke,pddpo)


      RETURN
      END
