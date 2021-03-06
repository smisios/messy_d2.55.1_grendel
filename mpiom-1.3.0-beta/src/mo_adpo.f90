MODULE MO_ADPO

  USE MO_PARAM1

  REAL ,ALLOCATABLE :: TRP(:,:,:),TRM(:,:,:)

#ifdef bounds_exch_tp
  REAL ,ALLOCATABLE :: TRPHELP(:,:),TRMHELP(:,:)
#endif 

  REAL ,ALLOCATABLE :: WTP(:,:,:),WTM(:,:,:)

CONTAINS

SUBROUTINE ALLOC_MEM_ADPO

  ALLOCATE(TRP(IE,JE,KEP),TRM(IE,JE,0:KE))


#ifdef bounds_exch_tp
   ALLOCATE(TRPHELP(IE,KE),TRMHELP(IE,KE))
#endif 

  ALLOCATE(WTP(IE,JE,KEP),WTM(IE,JE,0:KE))

END SUBROUTINE ALLOC_MEM_ADPO



SUBROUTINE OCADPO(TRF)

  USE MO_PARAM1
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_COMMOAU1
  USE MO_COMMOBBL
  USE MO_UNITS

  REAL :: TRF(IE,JE,KE)

  !     ROUTINE OCADPO
  !
  !     COMPUTES ADVECTION OF ONE TRACERFIELD
  !
  !     BY ERNST MAIER-REIMER 1/1999
  !     MODIFIED 1/2000 UWE MIKOLAJEWICZ
  !     MAKE MASS CONSERVING WITH FREE SURFACE LAYER
  !     MODIFIED 3/2000 UWE MIKOLAJEWICZ
  !     INCLUDE BOTTOM BOUNDARY LAYER TRANSPORTS IN
  !     ADVECTION SCHEME (IOCAD.EQ.4 :: ADPO + SLOPECON_ADPO)
  !     NEEDS TO BE RUN TOGETHER WITH SR SLOPETRANS
  !
  !     INPUT/OUTPUT  
  !     TRF(IE,JE,KE)     TRACER FIELD 
  !
  !     USES VELOCITIES (UKO,VKE,WO)

!#ifdef bounds_exch_save
  CALL bounds_exch('p',trf,'mo_adpo 1') 
!#endif
!$OMP PARALLEL PRIVATE(i,j,k,klo,kup,suminf,abl,zwabl,up,um,scha,sch,schi,rn,zzsurf)

if(IOCAD.eq.4)then

  DO K=1,KEP
!$OMP DO
    DO J=1,JE
      DO I=1,IE
        WOBACK(I,J,K)=WO(I,J,K)
      END DO
    END DO
!$OMP END DO
  END DO

!     COMPUTE NEW VERTICAL VELOCITIES INCLUDING BBL TRANSPORTS

!$OMP DO
  DO J=2,JE1
    DO K=KE,1,-1
      DO I=2,IE1
        SUMINF=0.
        
        IF(KDWUBBL(I-1,J).EQ.K)SUMINF=SUMINF+UBBL(I-1,J)
        IF(KUPUBBL(I-1,J).EQ.K)SUMINF=SUMINF-UBBL(I-1,J)
        IF(KDWUBBL(I,J).EQ.K)SUMINF=SUMINF-UBBL(I,J)
        IF(KUPUBBL(I,J).EQ.K)SUMINF=SUMINF+UBBL(I,J)
        
        IF(KDWVBBL(I,J).EQ.K)SUMINF=SUMINF+VBBL(I,J)
        IF(KUPVBBL(I,J).EQ.K)SUMINF=SUMINF-VBBL(I,J)
        IF(KDWVBBL(I,J-1).EQ.K)SUMINF=SUMINF-VBBL(I,J-1)
        IF(KUPVBBL(I,J-1).EQ.K)SUMINF=SUMINF+VBBL(I,J-1)
        
        WO(I,J,K)=WO(I,J,K+1)+SUMINF/(DLXP(I,J)*DLYP(I,J))
        
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
 

  DO K=1,KE
!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        WO(I,J,K)=WOBACK(I,J,K)+WO(I,J,K)
      ENDDO
    ENDDO
!$OMP END DO
  ENDDO


!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',WO,'mo_adpo 2')
!$OMP END SINGLE
!#endif
 
endif    ! iocad=4   


  DO K=1,KE
!$OMP DO
    DO J=1,JE
      DO I=1,IE
        S1O(I,J,K)=WETO(I,J,K)*DLXP(I,J)*DLYP(I,J)*DDPO(I,J,K)
        T1O(I,J,K)=TRF(I,J,K)
        TRP(I,J,K)=0.
        TRM(I,J,K)=0.
        WTP(I,J,K)=0.
        WTM(I,J,K)=0.
      END DO
    END DO
!$OMP END DO
  END DO
  
!$OMP DO
  DO J=1,JE
    DO I=1,IE
      TRM(I,J,0)=0.
      WTM(I,J,0)=0.
      TRP(I,J,KEP)=0.
      WTP(I,J,KEP)=0.
      S1O(I,J,1)=S1O(I,J,1)+DLXP(I,J)*DLYP(I,J)*WETO(I,J,1)*              &
           (ZO(I,J)-WO(I,J,1)*DT-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)
!      S1O(I,J,1)=S1O(I,J,1)+DLXP(I,J)*DLYP(I,J)*WETO(I,J,1)*              &
!           (ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)
    ENDDO
  ENDDO
!$OMP END DO

  
!      VERTICAL TRANSPORTS
  

  DO K=1,KE
    zzsurf=1.
    if(k.eq.1)zzsurf=0.
    KLO=MIN(K+1,KE)
    KUP=MAX(K-1,1)
!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I,J,KUP)-TRF(I,J,KLO))
        ZWABL=ABS(TRF(I,J,KUP)+TRF(I,J,KLO)-2.*TRF(I,J,K))
        UP=0.5*DT*DLYP(I,J)*DLXP(I,J)*(WO(I,J,K)+ABS(WO(I,J,K)))*zzsurf
        UM=0.5*DT*DLYP(I,J)*DLXP(I,J)*(ABS(WO(I,J,K+1))-WO(I,J,K+1))
        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))

!CUWE       SCH=MIN(1.,SCHA*S1O(I,J,K)/(UP+UM+1.E-20))

#ifdef SMOADV
        SCH=MIN(WETO(I,J,KLO),SCHA)
#else
        SCH=MIN(WETO(I,J,KLO),SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
#endif
        IF(K.EQ.1.OR.K.EQ.KE)SCH=0.

        IF(IOCAD.EQ.1)SCH=0.  ! PURE UPWIND

        SCHI=1.-SCH  

        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KUP)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KLO)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
!$OMP END DO
  END DO


  
!$OMP DO
  DO J=1,JE
    DO I=1,IE
      TRP(I,J,KEP)=0.
      WTP(I,J,KEP)=0.
      TRM(I,J,0)=0.
      WTM(I,J,0)=0.
    END DO
  END DO
!$OMP END DO

! RJ: no boundary exchange necessary here, since TRP, WTP, TRM, WTM
!     are used only in the inner domain in the following loop


  DO K=1,KE
!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        IF(WETO(I,J,K).GT.0.5) THEN
          RN=S1O(I,J,K)+WTP(I,J,K+1)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J,K-1)
          T1O(I,J,K)=(TRF(I,J,K)*S1O(I,J,K)+TRP(I,J,K+1)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J,K-1))/RN
          S1O(I,J,K)=RN
        ENDIF
      END DO
    END DO
!$OMP END DO

!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        TRF(I,J,K)=T1O(I,J,K)
      END DO
    END DO
!$OMP END DO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',TRF,'mo_adpo 3')
  CALL bounds_exch('p',S1O,'mo_adpo 4')
!$OMP END SINGLE
!#endif

!      TRANSPORTS IN X-DIRECTION


  DO K=1,KE
!$OMP DO
    DO J=1,JE
      DO I=1,IE
        TRP(I,J,K)=0.
        TRM(I,J,K)=0.
        WTP(I,J,K)=0.
        WTM(I,J,K)=0.
      END DO
    END DO
!$OMP END DO
  END DO

 

  DO K=1,KE
     if(iocad.eq.4)then 
!$OMP DO
        DO J=1,JE
           DO I=1,IE
              BBU(I,J)=0.
              IF(KDWUBBL(I,J).EQ.K)BBU(I,J)=UBBL(I,J)
              IF(KUPUBBL(I,J).EQ.K)BBU(I,J)=-UBBL(I,J)
           END DO
        END DO
!$OMP END DO
     endif
    
!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I+1,J,K)-TRF(I-1,J,K))
        ZWABL=ABS(TRF(I+1,J,K)+TRF(I-1,J,K)-2.*TRF(I,J,K))

        if(iocad.eq.4)then
           UP=0.5*DT*DDUO(I,J,K)*DLYU(I,J)*(UKO(I,J,K)+ABS(UKO(I,J,K)))     &
                +0.5*DT*(BBU(I,J)+ABS(BBU(I,J)))
        else
           UP=0.5*DT*DDUO(I,J,K)*DLYU(I,J)*(UKO(I,J,K)+ABS(UKO(I,J,K)))
        endif

        if(iocad.eq.4)then
           UM=0.5*DT*DDUO(I-1,J,K)*(ABS(UKO(I-1,J,K))-UKO(I-1,J,K))         &
                *DLYU(I-1,J)                                                &
                +0.5*DT*(ABS(BBU(I-1,J))-BBU(I-1,J))
        else
           UM=0.5*DT*DDUO(I-1,J,K)*(ABS(UKO(I-1,J,K))-UKO(I-1,J,K))         &
                *DLYU(I-1,J)
        endif

        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO(I,J,K)                &
             *WETO(I+1,J,K)*WETO(I-1,J,K)
#ifdef SMOADH
        SCH=MIN(1.,SCHA)
#else
        SCH=MIN(1.,SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
#endif

        IF(IOCAD.EQ.1)SCH=0.   !PURE UPWIND

        SCHI=1.-SCH

        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I+1,J,K)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I-1,J,K)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
!$OMP END DO
  END DO
  
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('u+',TRP,'mo_adpo 5')
  CALL bounds_exch('u+',TRM,'mo_adpo 6')
  CALL bounds_exch('u+',WTP,'mo_adpo 7')
  CALL bounds_exch('u+',WTM,'mo_adpo 8')
!$OMP END SINGLE
!#endif  

  DO K=1,KE
!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        IF(WETO(I,J,K).GT.0.5) THEN   
          RN=S1O(I,J,K)+WTP(I-1,J,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I+1,J,K)
          T1O(I,J,K)=(TRF(I,J,K)*S1O(I,J,K)+TRP(I-1,J,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I+1,J,K))/RN
          S1O(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        TRF(I,J,K)=T1O(I,J,K)
      END DO
    END DO
!$OMP END DO
  END DO
  

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',TRF,'mo_adpo 9')
  CALL bounds_exch('p',S1O,'mo_adpo 10')
!$OMP END SINGLE
!#endif 
 
!      TRANSPORTS IN Y-DIRECTION



  DO K=1,KE
!$OMP DO
    DO J=1,JE
      DO I=1,IE
        TRP(I,J,K)=0.
        TRM(I,J,K)=0.
        WTP(I,J,K)=0.
        WTM(I,J,K)=0.
      END DO
    END DO
!$OMP END DO
  END DO
 



  DO K=1,KE

     if(iocad.eq.4)then
! INITIALIZE BBL TRANSPORTS FOR LEVEL
!$OMP DO
        DO J=1,JE
           DO I=1,IE
              BBV(I,J)=0.
              IF(K.EQ.KDWVBBL(I,J))BBV(I,J)=VBBL(I,J)
              IF(K.EQ.KUPVBBL(I,J))BBV(I,J)=-VBBL(I,J)
           ENDDO
        ENDDO
!$OMP END DO
     endif

!$OMP DO    
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I,J-1,K)-TRF(I,J+1,K))
        ZWABL=ABS(TRF(I,J-1,K)+TRF(I,J+1,K)-2.*TRF(I,J,K))
        if(iocad.eq.4)then
           UP=0.5*DT*DDUE(I,J-1,K)*DLXV(I,J-1)                              &
                *(VKE(I,J-1,K)+ABS(VKE(I,J-1,K)))                           &
                +0.5*DT*(BBV(I,J-1)+ABS(BBV(I,J-1)))
        else
           UP=0.5*DT*DDUE(I,J-1,K)*DLXV(I,J-1)                              &
                *(VKE(I,J-1,K)+ABS(VKE(I,J-1,K)))
        endif
        if(iocad.eq.4)then
           UM=0.5*DT*DDUE(I,J,K)*DLXV(I,J)*(ABS(VKE(I,J,K))-VKE(I,J,K))     &
                +0.5*DT*(ABS(BBV(I,J))-BBV(I,J))
        else
           UM=0.5*DT*DDUE(I,J,K)*DLXV(I,J)*(ABS(VKE(I,J,K))-VKE(I,J,K))
        endif
        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO(I,J,K)                &
             *WETO(I,J-1,K)*WETO(I,J+1,K)
#ifdef SMOADH
        SCH=MIN(1.,SCHA)
#else
        SCH=MIN(1.,SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
#endif

        IF(IOCAD.EQ.1)SCH=0.   !PURE UPWIND

        SCHI=1.-SCH
        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J-1,K)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J+1,K)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
!$OMP end DO
  END DO

!#ifdef bounds_exch_save  
!$OMP SINGLE
  CALL bounds_exch('v+',TRP,'mo_adpo 11')
  CALL bounds_exch('v+',TRM,'mo_adpo 12')
  CALL bounds_exch('v+',WTP,'mo_adpo 13')
  CALL bounds_exch('v+',WTM,'mo_adpo 14')
!$OMP END SINGLE
!#endif

#ifdef bounds_exch_tp
      do k=1,ke
      do i=2,ie-1
!      trm(i,1,k)=trp(ie+2-i,3,k)
!      wtm(i,1,k)=wtp(ie+2-i,3,k)
         trm(i,1,k)=trp(i,1,k)    
         wtm(i,1,k)=wtp(i,1,k)
      enddo
      enddo
#endif

  

  DO K=1,KE
!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        IF(WETO(I,J,K).GT.0.5) THEN
          
          RN=S1O(I,J,K)+WTP(I,J+1,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J-1,K)
          T1O(I,J,K)=(TRF(I,J,K)*S1O(I,J,K)+TRP(I,J+1,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J-1,K))/RN
          S1O(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        TRF(I,J,K)=T1O(I,J,K)
      END DO
    END DO
!$OMP END DO
  END DO

!#ifdef bounds_exch_save  
!$OMP SINGLE
  CALL bounds_exch('p',TRF,'mo_adpo 15')
!$OMP END SINGLE
!#endif 
 
!     RESET VERTICAL VELOCITIES TO INITIAL VALUES WITHOUT BBL TRANSPORTS
  
  if(iocad.eq.4)then
     DO K=1,KEP
!$OMP DO
        DO J=1,JE
           DO I=1,IE
              WO(I,J,K)=WOBACK(I,J,K)
           END DO
        END DO
!$OMP END DO
     END DO
  endif

!$OMP END PARALLEL  

END SUBROUTINE OCADPO

END MODULE MO_ADPO
