      SUBROUTINE POWACH(kpie,kpje,kpke,pdlxp,pdlyp,psao,pwo)
!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/powach.f90,v $\\
!$Revision: 1.2.10.1.2.2.4.1.2.2.2.3.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!$Name: mpiom_1_2_0 $\\
!
!**********************************************************************
!
!**** *POWACH* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     S.Lorenz/JO.Beismann, OpenMP parallel    *MPI-Met, HH*  24.08.07
!
!     Purpose
!     -------
!     .
!
!     Method
!     -------
!     .
!
!**   Interface.
!     ----------
!
!     *CALL*       *POWACH*
!
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL :: of model grid.
!     *INTEGER* *kpje*    - 2nd REAL :: of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL :: of model grid.
!     *REAL*    *psao*    - potential temperature [deg C].
!     *REAL*    *pwo*     - vertical velocity in scalar points [m/s].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!
!     Externals
!     ---------
!     none.
!ssso12 refersto the phosphorus of organic debris
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      use mo_param1_bgc

implicit none

      INTEGER :: i,j,k, iter
      INTEGER :: kpie,kpje,kpke

      REAL(wp) :: psao(kpie,kpje,kpke)

      REAL(wp) :: sedb1(kpie,0:ks),sediso(kpie,0:ks)
      REAL(wp) :: solrat(kpie,ks),powcar(kpie,ks)
      REAL(wp) :: aerob(kpie,ks),anaerob(kpie,ks),ansulf(kpie,ks)
      REAL(wp) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)

#ifdef __c_isotopes
      REAL(wp) :: ratc13, ratc14, rato13, rato14, poso13, poso14
#endif
      REAL(wp) :: disso, dissot, undsa, silsat, posol, dissot1, dissot2
      REAL(wp) :: umfa,denit,bt,alk,c
      REAL(wp) :: ak1,ak2,akb,akw
      REAL(wp) :: h,t1,t2,a,dadh,dddhhh,satlev

! *****************************************************************
! accelerated sediment
! needed for boundary layer ventilation in fast sediment routine

      REAL(wp) :: pwo(kpie,kpje,kpke+1)
      REAL(wp) :: bolven(kpie)

! A LOOP OVER J
! RJ: This loop must go from 1 to kpje in the parallel version,
!     otherways we had to do a boundary exchange
!e    write(0,737)nsedtra
!737   format('powach',i6)

! Silicate saturation concentration is 1 mol/m3

      silsat=0.001_wp

! Dissolution rate constant of opal (disso) [1/(kmol Si(OH)4/m3)*1/sec]

!      disso=1.e-8_wp
      disso=1.e-6_wp ! test vom 03.03.04 half live sil ca. 20.000 yr
      dissot=disso*dtbgc

! Degradation rate constant of POP (disso) [1/(kmol O2/m3)*1/sec]

      disso=0.01_wp/86400._wp  !  disso=3.e-5 was quite high
      dissot1=disso*dtbgc

! Dissolution rate constant of CaCO3 (disso) [1/(kmol CO3--/m3)*1/sec]
      disso=1.e-7_wp
      dissot2=disso*dtbgc

!ik      denit = 1.e-6*dtbgc
      denit = 0.01_wp / 86400._wp *dtbgc

!$OMP PARALLEL PRIVATE (                                                &
!$OMP          bolven,undsa,posol,umfa,bt,alk,c, ak1,ak2,akb,akw,       &
!$OMP          ratc13,ratc14,rato13,rato14,poso13,poso14,               &
!$OMP          h,t1,t2,a,dadh,dddhhh,satlev,                            &
!$OMP          sedb1,sediso,solrat,powcar,aerob,anaerob)

!$OMP DO
      DO 8888 j=1,kpje

      DO 1189 k=1,ks
      DO 1189 i=1,kpie
         solrat(i,k) =0._wp
         powcar(i,k) =0._wp
         anaerob(i,k)=0._wp
         aerob(i,k)  =0._wp
         ansulf(i,k) =0._wp
1189  CONTINUE

! calculate bottom ventilation rate for scaling of sediment-water exchange
      do 1170 i=1,kpie
      bolven(i) = 1._wp
1170  continue

      DO 1171 k=0,ks
      DO 1171 i=1,kpie
        sedb1(i,k)=0._wp
        sediso(i,k)=0._wp
1171  CONTINUE



! CALCULATE SILICATE-OPAL CYCLE AND SIMULTANEOUS SILICATE DIFFUSION
!******************************************************************

! Evaluate boundary conditions for sediment-water column exchange.
! Current undersaturation of bottom water: sedb(i,0) and
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

      DO 3 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            undsa=silsat-powtra(i,j,1,ipowasi)
            sedb1(i,0)=bolay(i,j)*(silsat-ocetra(i,j,kbo(i,j),isilica)) &
     &                 *bolven(i)
            solrat(i,1)=                                                &
     &      (sedlay(i,j,1,issssil)+silpro(i,j)/(porsol(1)*seddw(1)))    &
     &      *dissot/(1._wp + dissot*undsa)*porsol(1)/porwat(1)
         ENDIF
3     CONTINUE

! Evaluate sediment undersaturation and degradation.
! Current undersaturation in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k)

      DO 2 k=1,ks
      DO 2 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            undsa=silsat-powtra(i,j,k,ipowasi)
            sedb1(i,k)=seddw(k)*porwat(k)*(silsat-powtra(i,j,k,ipowasi))
            IF(k.GT.1)solrat(i,k)=sedlay(i,j,k,issssil)                 &
     &                 *dissot/(1._wp+dissot*undsa)*porsol(k)/porwat(k)
         ENDIF
2     CONTINUE

! Solve for new undersaturation sediso, from current undersaturation sedb1,
! and first guess of new solid sediment solrat.


      CALL powadi(j,kpie,solrat,sedb1,sediso,bolven)

! Update water column silicate, and store the flux for budget.
! Add biogenic opal flux to top sediment layer.

      DO 4 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
         sedfluxo(i,j,ipowasi)=sedfluxo(i,j,ipowasi) +                &
     &  (silsat-sediso(i,0)-ocetra(i,j,kbo(i,j),isilica))*bolay(i,j)

            ocetra(i,j,kbo(i,j),isilica)=silsat-sediso(i,0)
            sedlay(i,j,1,issssil)=                                    &
     &        sedlay(i,j,1,issssil)+silpro(i,j)/(porsol(1)*seddw(1))
            silpro(i,j)=0._wp
         ENDIF
4     CONTINUE

! Calculate updated degradation rate from updated undersaturation.
! Calculate new solid sediment.
! Update pore water concentration from new undersaturation.

      DO 5 k=1,ks
      DO 5 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            solrat(i,k)=sedlay(i,j,k,issssil)                        &
     &                  *dissot/(1._wp+dissot*sediso(i,k))
            posol=sediso(i,k)*solrat(i,k)
            sedlay(i,j,k,issssil)=                                   &
     &          sedlay(i,j,k,issssil)-posol
            powtra(i,j,k,ipowasi)=silsat-sediso(i,k)
         ENDIF
5     CONTINUE


! CALCULATE OXYGEN-POC CYCLE AND SIMULTANEOUS OXYGEN DIFFUSION
!*************************************************************

! This scheme is not based on undersaturation, but on O2 itself

! Evaluate boundary conditions for sediment-water column exchange.
! Current concentration of bottom water: sedb(i,0) and
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

      DO 13 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            undsa=powtra(i,j,1,ipowaox)
            sedb1(i,0)=bolay(i,j)*ocetra(i,j,kbo(i,j),ioxygen)         &
     &                 *bolven(i)
            solrat(i,1)=                                               &
     &       (sedlay(i,j,1,issso12)+prorca(i,j)/(porsol(1)*seddw(1)))  &
     &          *ro2ut*dissot1/(1._wp + dissot1*undsa)*porsol(1)/porwat(1)
         ENDIF
13    CONTINUE

! Evaluate sediment concentration and degradation.
! Current concentration in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k)

      DO 12 k=1,ks
      DO 12 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            undsa=powtra(i,j,k,ipowaox)
            sedb1(i,k)=seddw(k)*porwat(k)*powtra(i,j,k,ipowaox)
            IF(k.GT.1)solrat(i,k)=sedlay(i,j,k,issso12)               &
     &         *ro2ut*dissot1/(1._wp+dissot1*undsa)*porsol(k)/porwat(k)
         ENDIF
12     CONTINUE

! Solve for new O2 concentration sediso, from current concentration sedb1,
! and first guess of new solid sediment solrat.

      CALL powadi(j,kpie,solrat,sedb1,sediso,bolven)

! Update water column oxygen, and store the flux for budget (opwflux). ! js: opwflux not in present model code
! Add organic carbon flux 'prorca' to top sediment layer.

      DO 14 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            ocetra(i,j,kbo(i,j),ioxygen)=sediso(i,0)
            sedlay(i,j,1,issso12)                                     &
     &      =sedlay(i,j,1,issso12)+prorca(i,j)/(porsol(1)*seddw(1))
#ifdef __c_isotopes
            sedlay(i,j,1,issso13)                                     &
     &      =sedlay(i,j,1,issso13)+pror13(i,j)/(porsol(1)*seddw(1))
            sedlay(i,j,1,issso14)                                     &
     &      =sedlay(i,j,1,issso14)+pror14(i,j)/(porsol(1)*seddw(1))
#endif

         prorca(i,j)=0._wp
#ifdef __c_isotopes
         pror13(i,j)=0._wp
         pror14(i,j)=0._wp
#endif
         ENDIF
14    CONTINUE


! Calculate updated degradation rate from updated concentration.
! Calculate new solid sediment.
! Update pore water concentration.
! Store flux in array aerob, for later computation of DIC and alkalinity.
      DO 15 k=1,ks
            umfa=porsol(k)/porwat(k)
      DO 15 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            solrat(i,k)=sedlay(i,j,k,issso12)                         &
     &                 *dissot1/(1._wp+dissot1*sediso(i,k))
            posol=sediso(i,k)*solrat(i,k)
#ifdef __c_isotopes
            rato13=sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+1.e-24_wp)
            rato14=sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+1.e-24_wp)
            poso13=posol*rato13
            poso14=posol*rato14
#endif
            aerob(i,k)=posol*umfa !this has P units: kmol P/m3 of pore water
            sedlay(i,j,k,issso12)=sedlay(i,j,k,issso12)-posol
            powtra(i,j,k,ipowaph)=powtra(i,j,k,ipowaph)+posol*umfa
            powtra(i,j,k,ipowno3)=powtra(i,j,k,ipowno3)+posol*rnit*umfa
            powtra(i,j,k,ipowaox)=sediso(i,k)
#ifdef __c_isotopes
            sedlay(i,j,k,issso13)=sedlay(i,j,k,issso13)-poso13
            sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)-poso14
            powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowc13)+poso13*umfa
            powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)+poso14*umfa
#endif
         ENDIF
15    CONTINUE

! CALCULATE NITRATE REDUCTION UNDER ANAEROBIC CONDITIONS EXPLICITELY
!*******************************************************************

! Denitrification rate constant of POP (disso) [1/sec]
! Store flux in array anaerob, for later computation of DIC and alkalinity.

      DO 124 k=1,ks
         umfa=porsol(k)/porwat(k)
      DO 124 i=1,kpie
         IF (bolay(i, j) .GT. 0._wp) THEN
         IF (powtra(i, j, k, ipowaox) .LT. 1.e-6_wp) THEN
           posol = denit * MIN(0.5_wp * powtra(i, j, k, ipowno3)/nitdem, &
     &                          sedlay(i,j,k,issso12))
#ifdef __c_isotopes
          rato13=sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+1.e-24_wp)
          rato14=sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+1.e-24_wp)
#endif
            anaerob(i,k)=posol*umfa !this has P units: kmol P/m3 of pore water
            sedlay(i,j,k,issso12)=sedlay(i,j,k,issso12)-posol
            powtra(i,j,k,ipowaph)=powtra(i,j,k,ipowaph)+posol*umfa
!            powtra(i,j,k,ipowno3)=powtra(i,j,k,ipowno3)-98.*posol*umfa
!            powtra(i,j,k,ipown2)=powtra(i,j,k,ipown2)+57.*posol*umfa
!tk27012010 changed no3 use and n2 production in denitrification
            powtra(i,j,k,ipowno3)=powtra(i,j,k,ipowno3)-nitdem*posol*umfa
            powtra(i,j,k,ipown2)=powtra(i,j,k,ipown2)+n2prod*posol*umfa
            powh2obud(i,j,k)=powh2obud(i,j,k)+0.5_wp*n2prod*posol*umfa

#ifdef __c_isotopes
            poso13=posol*rato13
            poso14=posol*rato14
            sedlay(i,j,k,issso13)=sedlay(i,j,k,issso13)-poso13
            sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)-poso14
            powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowc13)+poso13*umfa
            powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)+poso14*umfa
#endif

         ENDIF   ! oxygen <1.e-6
         endif   ! bolay
124   CONTINUE

!    sulphate reduction in sediments
      DO 125 k=1,ks
         umfa=porsol(k)/porwat(k)
      DO 125 i=1,kpie
         if(bolay(i,j).gt.0._wp) then
         IF(powtra(i,j,k,ipowaox).LT.1.e-6_wp) THEN
            posol=denit*0.01_wp * sedlay(i,j,k,issso12)
#ifdef __c_isotopes
            rato13=sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+1.e-24_wp)
            rato14=sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+1.e-24_wp)
#endif
            ansulf(i,k)=posol*umfa !this has P units: kmol P/m3 of pore water
            sedlay(i,j,k,issso12)=sedlay(i,j,k,issso12)-posol
            powtra(i,j,k,ipowaph)=powtra(i,j,k,ipowaph)+posol*umfa
            powtra(i,j,k,ipowno3)=powtra(i,j,k,ipowno3)+posol*umfa*rno3
            powh2obud(i,j,k)=powh2obud(i,j,k)-ro2ut*posol*umfa
#ifdef __c_isotopes
            poso13=posol*rato13
            poso14=posol*rato14
            sedlay(i,j,k,issso13)=sedlay(i,j,k,issso13)-poso13
            sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)-poso14
            powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowc13)+poso13*umfa
            powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)+poso14*umfa
#endif

         endif
         ENDIF
125   CONTINUE

!        if(j.eq.kpje)       call maschk(kpie,kpje,kpke,12)

! CALCULATE CaCO3-CO3 CYCLE AND SIMULTANEOUS CO3-UNDERSATURATION DIFFUSION
!*************************************************************************


! COMPUTE NEW POWCAR=CARBONATE ION CONCENTRATION IN THE SEDIMENT
! FROM CHANGED ALKALINITY (NITRATE PRODUCTION DURING REMINERALISATION)
! AND DIC GAIN. ITERATE 5 TIMES. THIS CHANGES PH (SEDHPL) OF SEDIMENT.

      DO 10 ITER=1,5

      DO 1 K=1,KS
      DO 1 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            bt=rrrcl*psao(i,j,kbo(i,j))
!alkalinity is increased during denitrification due to consumption of H+ via NO3 (see Wolf-Gladrow etal,2007)
            alk=powtra(i,j,k,ipowaal)-(ansulf(i,k)+aerob(i,k))*rnit + nitdem*anaerob(i,k)
            c=powtra(i,j,k,ipowaic)+(anaerob(i,k)+aerob(i,k)+ansulf(i,k))*rcar
            ak1=ak13(i,j,kbo(i,j))
            ak2=ak23(i,j,kbo(i,j))
            akb=akb3(i,j,kbo(i,j))
            akw=akw3(i,j,kbo(i,j))
            h=sedhpl(i,j,k)
            t1=h/ak1
            t2=h/ak2
            a=c*(2._wp+t2)/(1._wp+t2+t2*t1)  +akw/h-h+bt/(1._wp + h/akb)-alk
            dadh=c*( 1._wp/(ak2*(1._wp+t2+t2*t1))-(2._wp+t2)*(1._wp/ak2+2._wp*t1/ak2)/  &
     &          (1._wp+t2+t2*t1)**2)                                        &
     &          -akw/h**2-1._wp-(bt/akb)/(1._wp+h/akb)**2
            dddhhh=a/dadh
            sedhpl(i,j,k) = MAX(h - dddhhh, 1.e-11_wp)
            powcar(i, k) = c / (1._wp + t2 * (1._wp + t1))
         ENDIF
1     CONTINUE

10    CONTINUE

! Evaluate boundary conditions for sediment-water column exchange.
! Current undersaturation of bottom water: sedb(i,0) and
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

! CO3 saturation concentration is aksp/calcon as in CARCHM
! (calcon defined in BELEG_BGC with 1.03e-2; 1/calcon =~ 97.)

      DO 23 i=1,kpie
         IF(bolay(i,j) .GT. 0._wp) THEN
            satlev=aksp(i,j,kbo(i,j))/calcon+2.e-5_wp
            undsa = MAX(satlev - powcar(i, 1), 0._wp)
            sedb1(i,0)=bolay(i,j)*(satlev-co3(i,j,kbo(i,j)))             &
     &                 *bolven(i)
            solrat(i,1)=                                                 &
     &         (sedlay(i,j,1,isssc12)+prcaca(i,j)/(porsol(1)*seddw(1)))  &
     &         *dissot2/(1._wp + dissot2*undsa)*porsol(1)/porwat(1)
         ENDIF
23     CONTINUE

! Evaluate sediment undersaturation and degradation.
! Current undersaturation in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k)

      DO 22 k=1,ks
      DO 22 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            undsa = MAX(aksp(i, j, kbo(i, j)) / calcon - powcar(i, k), 0._wp)
            sedb1(i,k)=seddw(k)*porwat(k)*undsa
            IF(k.GT.1)solrat(i,k)=sedlay(i,j,k,isssc12)                 &
     &          *dissot2/(1._wp+dissot2*undsa)*porsol(k)/porwat(k)
            IF(undsa.LE.0._wp) solrat(i,k)=0._wp
         ENDIF
22     CONTINUE

! Solve for new undersaturation sediso, from current undersaturation sedb1,
! and first guess of new solid sediment solrat.

      CALL powadi(j,kpie,solrat,sedb1,sediso,bolven)

! There is no exchange between water and sediment with respect to co3 so far.
! Add calcite flux 'prcaca' to uppermost sediment layer.
      DO 24 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
            sedlay(i,j,1,isssc12)=                                     &
     &      sedlay(i,j,1,isssc12)+prcaca(i,j)/(porsol(1)*seddw(1))
#ifdef __c_isotopes
            sedlay(i,j,1,isssc13)=                                     &
     &      sedlay(i,j,1,isssc13)+prca13(i,j)/(porsol(1)*seddw(1))
            sedlay(i,j,1,isssc14)=                                     &
     &      sedlay(i,j,1,isssc14)+prca14(i,j)/(porsol(1)*seddw(1))
#endif
         prcaca(i,j)=0._wp
#ifdef __c_isotopes
         prca13(i,j)=0._wp
         prca14(i,j)=0._wp
#endif
         ENDIF
24    CONTINUE

! Calculate updated degradation rate from updated undersaturation.
! Calculate new solid sediment.
! No update of powcar pore water concentration from new undersaturation so far.
! Instead, only update DIC, and, of course, alkalinity.
! This also includes gains from aerobic and anaerobic decomposition.

      DO 25 k=1,ks
           umfa=porsol(k)/porwat(k)
      DO 25 i=1,kpie
         IF(bolay(i,j).GT.0._wp) THEN
#ifdef __c_isotopes
          ratc13=sedlay(i,j,k,isssc13)/(sedlay(i,j,k,isssc12)+1.e-24_wp)
          ratc14=sedlay(i,j,k,isssc14)/(sedlay(i,j,k,isssc12)+1.e-24_wp)
#endif
           solrat(i,k)=sedlay(i,j,k,isssc12)                           &
     &                 *dissot2/(1._wp+dissot2*sediso(i,k))
           posol=sediso(i,k)*solrat(i,k)
           sedlay(i,j,k,isssc12)=sedlay(i,j,k,isssc12)-posol
           powtra(i,j,k,ipowaic)=powtra(i,j,k,ipowaic)                 &
     &        +posol*umfa+(aerob(i,k)+anaerob(i,k)+ansulf(i,k))*rcar
           powtra(i,j,k,ipowaal)=powtra(i,j,k,ipowaal)                 &
     &        +2._wp*posol*umfa-rnit*(aerob(i,k)+ansulf(i,k))+nitdem*anaerob(i,k)
           pown2bud(i,j,k)=pown2bud(i,j,k)+2._wp*n2prod*anaerob(i,k)
#ifdef __c_isotopes
            poso13=posol*ratc13
            poso14=posol*ratc14
       sedlay(i,j,k,isssc13)=sedlay(i,j,k,isssc13)-poso13
       sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc14)-poso14
       powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowc13)+poso13*umfa
       powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)+poso14*umfa
#endif

         ENDIF
25    CONTINUE


8888  CONTINUE
!$OMP END DO
!$OMP END PARALLEL

#ifdef CHK_OCPROD
       call maschk(kpie,kpje,kpke,13)
#endif
      CALL DIPOWA(kpie,kpje,kpke,pdlxp,pdlyp,pwo)
#ifdef CHK_OCPROD
       call maschk(kpie,kpje,kpke,14)
#endif


!ik add clay sedimentation
!ik this is currently assumed to depend on total and corg sedimentation:
!ik f(POC) [kg C] / f(total) [kg] = 0.05
!ik thus it is
       do j=1,kpje
       do i=1,kpie
       sedlay(i,j,1,issster) = sedlay(i,j,1,issster)                 &
     &                       + produs(i,j)/(porsol(1)*seddw(1))
       enddo
       enddo


       DO 91 j=1,kpje
       DO 91 i=1,kpie
         silpro(i,j)=0._wp
         prorca(i,j)=0._wp
#ifdef __c_isotopes
         pror13(i,j)=0._wp
         pror14(i,j)=0._wp
         prca13(i,j)=0._wp
         prca14(i,j)=0._wp
#endif
         prcaca(i,j)=0._wp
         produs(i,j)=0._wp
91     CONTINUE

      RETURN
      END
