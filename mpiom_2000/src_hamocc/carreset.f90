      SUBROUTINE carreset                                              &
     &           (kpie,kpje,kpke,pddpo,psao,ptho,psicomo,             &
     &            pfu10,kplmon,kplday,kmonlen,pdlxp,pdlyp)
!**********************************************************************
#ifdef BPGC
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_timeser_bgc
      USE mo_control_bgc
      USE mo_bgcmean
      USE mo_param1_bgc

      USE mo_parallel

      USE MO_COMMO1, only: weto,zo,sictho,sicsno
      USE mo_planetary_constants, ONLY: rhoicwa, rhosnwa

      IMPLICIT none

      INTEGER :: i,j,k,l,kpie,kpje,kpke,itit
      INTEGER :: kplmon,kplday,kmonlen,laumo1

      REAL(wp) psao(kpie,kpje,kpke)
      REAL(wp) pddpo(kpie,kpje,kpke),a0(kpie,kpje,kpke)
      REAL(wp) psicomo(kpie,kpje)
      REAL(wp) pfu10(kpie,kpje)
      REAL(wp) ptho(kpie,kpje,kpke)
      REAL(wp) pdlxp(kpie,kpje),pdlyp(kpie,kpje)


      REAL(wp) :: supsat, undsa, dissol,r13,r14,rat13,rat14
      REAL(wp) :: fluxd,fluxu,flux14d,flux14u,flux13d,flux13u
      REAL(wp) :: dddhhh,dadh,a,h,c,alk,t1,t2
      REAL(wp) :: akbi,ak2i,ak1i,dtja,aa,bb,cc
      REAL(wp) :: kwco2,kwo2,kwdms,vol,volfak,carsum,vorfak,deldic,dicinv
      REAL(wp) :: scco2,sco2,scdms
      REAL(wp) :: contppm, Xconvxa
      REAL(wp) :: oxflux,niflux,dmsflux,n2oflux
      REAL(wp) :: ato2, atn2, atco2,pco2,atc13,atc14
      REAL(wp) :: AHI, ABE,RMONLEN,RPLDAY
      REAL(wp) :: AK0,AK1,AK2,AKB,AKW,BT,oxysa,anisa

#ifdef ANTC14
      REAL(wp) :: fantc14d,fantc14u
#endif

!      WRITE(*,*) 'CARCHM called with :',                             &
!     &           kpie,kpje,kpke,pddpo(50,50,1),psao(50,50,1),        &
!     &           ptho(50,50,1),psicomo(50,50),pfu10(50,50),          &
!     &           kplyear,kplmon,kplday,kmonlen
      laumo1=kplmon+1
      IF(laumo1.GT.12) laumo1=1
      kplday=0
      RMONLEN=kmonlen
      RPLDAY=kplday
      AHI=RPLDAY/RMONLEN
      ABE=1.-AHI


      volfak=0.
      carsum=0.
      do k=1,kpke
      do j=2,kpje-1
      do i=2,kpie-1
         vol   = pdlxp(i,j) * pdlyp(i,j) * pddpo(i,j,k)
         volfak=volfak + vol
         carsum=carsum + vol * weto(i,j,k) * ocetra(i,j,k,isco212)
      enddo
      enddo
      enddo
      CALL global_sum(volfak)
      CALL global_sum(carsum)

      vorfak=12.e-12

      IF (p_pe == p_io) WRITE(0,*) 'gigatons?', vorfak*carsum

      dicinv=vorfak*carsum
      deldic=(38500.-dicinv)/(volfak*vorfak)

!     -----------------------------------------------------------------
!*        22. CHEMICAL CONSTANTS - DEEP OCEAN
      do 43 itit=1,6
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
            hi(i,j,k)=max(hi(i,j,k)-dddhhh,1.e-10)
            co3(i,j,k)=                                                &
     &      c/(1.+hi(i,j,k)*(1.+hi(i,j,k)/ak13(i,j,k))/ak23(i,j,k))
! C14 decay, c14dec is set in BELEG_BGC
            a0(i,j,k)=c*h**2/(h**2+h*ak13(i,j,k)+ak13(i,j,k)*ak23(i,j,k))
         ENDIF
43    CONTINUE
!
      ocetra(:,:,:,isco212)=ocetra(:,:,:,isco212)+deldic

      carsum=0.
      do k=1,kpke
      do j=2,kpje-1
      do i=2,kpie-1
         vol   = pdlxp(i,j) * pdlyp(i,j) * pddpo(i,j,k)
         carsum=carsum + vol * weto(i,j,k) * ocetra(i,j,k,isco212)
      enddo
      enddo
      enddo
      CALL global_sum(carsum)

      IF (p_pe == p_io) WRITE(0,*) 'gigatons?', vorfak*carsum

      DO 143 k=1,kpke
      DO 143 j=1,kpje
      DO 143 i=1,kpie
         IF(pddpo(i,j,k).GT.0.5) THEN
            h=hi(i,j,k)
            c=ocetra(i,j,k,isco212)
            ak2=ak23(i,j,k)
            akw=akw3(i,j,k)
            bt=rrrcl*psao(i,j,k)
            akb=akb3(i,j,k)
            aa=a0(i,j,k)-c
            bb=a0(i,j,k)*ak13(i,j,k)
            cc=a0(i,j,k)*ak13(i,j,k)*ak23(i,j,k)
            h=-0.5*bb/aa+sqrt((0.5*bb/aa)**2-cc/aa)
            !  das naechste nur um confusion zu vermeiden
            hi(i,j,k)=h
            t1=h/ak13(i,j,k)
            t2=h/ak23(i,j,k)
            a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h/akb)
            ocetra(i,j,k,ialkali)=a
         ENDIF
143   CONTINUE

      carsum=0.
      do k=1,kpke
      do j=2,kpje-1
      do i=2,kpie-1
         vol   = pdlxp(i,j) * pdlyp(i,j) * pddpo(i,j,k)
         carsum=carsum + vol * weto(i,j,k) * ocetra(i,j,k,ialkali)
      enddo
      enddo
      enddo
      CALL global_sum(carsum)

      IF (p_pe == p_io) WRITE(0,*) 'gigatons alk?', vorfak*carsum

      RETURN
#endif
      END
