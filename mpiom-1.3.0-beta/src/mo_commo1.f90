      MODULE MO_COMMO1

      USE mo_param1

      IMPLICIT NONE

      INTEGER IMITTEL,IAUFR,IAUFW,JAHR,IKLIM,NNNDT,IMMA,JMMA
      INTEGER ISTART,I3DREST
      REAL DT,ROCON

      REAL ZERO,ONE,TWO,TENM4,FOUR,HALF,FOURTH,EIGHTH
      REAL G,AH,ROCP,ROCD,PI,RADIUS,OMEGA,RELTEM,RELSAL,DH,DTI
      REAL TIPHAS,COPHAS,SIPHAS
      INTEGER LITER &
           ,NDUMMY,N1,M1,N2,M2,N3,M3,N4,M4,KB,KBM,KM,MATR,NMAX &
           ,IUNIT,IFLAG,imean &
           ,LYEARS,LDTYEAR, LYEAR1,LYEAR2,LYEAR,LMONTH &
           ,LMONTS,LDTMONTH,LMONT1,LMONT2,MONLEN(12),NDTMONTH &
           ,LDAYS,LDAY,LDTDAYC,LDT &
           ,NYEARS,NMONTS,NDAYS,LY_END,LY_START,LM_START,IOASISFLUX &
           ,IMOCDIAG 

      LOGICAL LFORCEDIAG,LHFLDIAG,LCONVDIAG,LDIFFDIAG,LGMDIAG,LGRIDINFO &
             ,LCALCDIFI

      INTEGER IOCAD,ICONTRO

      REAL ALMZER,STABN,STABO,CONO,CONN, &
     &   GHO,GHN,DTDECT,CLICN,CLICO,AULAPTS,AULAPUV,AUS, &
     &   AH00,CAH00,AV0,DV0,WT,CSTABEPS,CAVOCON,CDVOCON
      INTEGER NTIPHA
!
!:: GRID DESCRIPTORS AND RELATED FIELDS
      REAL, ALLOCATABLE :: &
     &   DLXP(:,:),DLYP(:,:),DLXU(:,:),DLYV(:,:), &
     &   DLXV(:,:),DLYU(:,:),DTDXUO(:,:), &
     &   DTDXUE(:,:),DTDYO(:,:),DTDXPO(:,:),DTDXPE(:,:), &
     &   DPYE(:,:),DPYO(:,:), &
     &   DDUE(:,:,:),DDUO(:,:,:),DDPO(:,:,:),DPIO(:,:,:), &
     &   DDPSIO(:,:,:),DWI(:), &
     &   DEUTE(:,:),DEUTO(:,:),DEUTIE(:,:),DEUTIO(:,:), &
     &   DEPTO(:,:)

      REAL, ALLOCATABLE :: DZ(:),DZW(:),TIESTU(:),TIESTW(:),DI(:)

!:: GEOGRAPHIC POSITION:
      REAL, ALLOCATABLE :: GILA(:,:),GIPH(:,:),ALAT(:,:)
      REAL, ALLOCATABLE :: GILA_G(:,:),GIPH_G(:,:)
!:: LAND SEA MASKS ON SKALAR/VECTOR POINTS
      REAL, ALLOCATABLE :: WETO(:,:,:),AMSUE(:,:,:),AMSUO(:,:,:)

      REAL, ALLOCATABLE :: WETO_G(:,:,:)
!
!:: ONE-DIMENSIONAL FIELDS
      REAL, ALLOCATABLE :: &
     &   DOTRA(:),DOTRT(:),UPTRT(:),WINTUR(:), &
     &   DBACKV(:),ABACKV(:),TAF(:),SAF(:),PREFF(:), &
     &   UPTRA(:),ROREF(:),SKAL(:),B(:),X(:)
!
!:: THREE-DIMENSIONAL FIELDS
      REAL, ALLOCATABLE :: &
     &   RHOO(:,:,:),UKO(:,:,:),VKE(:,:,:), &
     &   THO(:,:,:),SAO(:,:,:),WO(:,:,:), &
     &   PO(:,:,:),VK1E(:,:,:), &
     &   UK1O(:,:,:),VOE(:,:,:),UOO(:,:,:), &
     &   T1O(:,:,:),S1O(:,:,:),DVO(:,:,:),AVO(:,:,:), &
     &   STABIO(:,:,:)
!
!:: TWO-DIMENSIONAL FIELDS
      REAL, ALLOCATABLE :: &
     &   ZO(:,:),Z1E(:,:),Z1O(:,:), &
     &   U1E(:,:),V1E(:,:),U1O(:,:),V1O(:,:), &
     &   UROT(:,:),UDIV(:,:),USNO(:,:),VSNE(:,:), &
     &   USO(:,:),VSE(:,:),UCOS(:,:),VCOS(:,:), &
     &   SHELP(:,:),THELP(:,:),RHELP(:,:), &
     &   TLOW(:,:),SLOW(:,:),TSUP(:,:),SSUP(:,:), &
     &   PXOIN(:,:),PYEIN(:,:),UZO(:,:),VZE(:,:), &
     &   FTWOU(:,:),FTWOV(:,:),CURVAV(:,:),B1E(:,:), &
     &   B1O(:,:)
      INTEGER, ALLOCATABLE :: NUM(:,:),KCONDEP(:,:),KBOT(:,:)
      INTEGER, ALLOCATABLE :: NUM_G(:,:)
!
!:: ICE-MODEL FIELDS
      REAL, ALLOCATABLE :: &
     &   SICTHO(:,:),SICOMO(:,:),SICDIO(:,:),SICSHO(:,:), &
     &   SICPO(:,:), SICOMP(:,:),                                 &
     &   SICUO(:,:),SICVE(:,:), TICE(:,:), &
     &   HIBDELO(:,:),HIBZETO(:,:),HIBETO(:,:), &
     &   HIBDELE(:,:),HIBZETE(:,:),HIBETE(:,:), &
     &   TAUWATU(:,:),TAUWATV(:,:), &
     &   SICSNO(:,:),FWO(:,:),SHO(:,:)


!  FIELDS FOR COMMON DIAGNOSTICS
      REAL, ALLOCATABLE :: &
     &   HFLM(:,:), &
     &   PMEM(:,:),EISTRX(:,:),EISTRY(:,:),AMLD(:,:,:),PSIUWE(:,:)
    
!
!:: FORCING FIELDS
      REAL, ALLOCATABLE :: &
     &   TXO(:,:),TYE(:,:),TAFO(:,:), &
     &   HEATO(:,:),EMINPO(:,:),RELSAO(:,:),RELTHO(:,:), &
     &   FCLOU(:,:),FSWR(:,:),FU10(:,:), &
     &   FPREC(:,:),FTDEW(:,:),RIVRUN(:,:)
#ifndef RIVER_GIRIV
      REAL RIVAL(NUMRIV,12),FRIV(NUMRIV)
      REAL RIVLON(NUMRIV),RIVLAT(NUMRIV),DDRIV(NUMRIV)
      INTEGER IRIVI(NUMRIV),IRIVJ(NUMRIV)
#endif /* RIVER_GIRIV*/

#ifdef GLACCALV
      REAL GLACLAT(NUMGLAC),GLACLON(NUMGLAC),GLACVAL(NUMGLAC)
      INTEGER IGLAC(NUMGLAC),JGLAC(NUMGLAC)
#endif

      REAL, ALLOCATABLE :: GIRIV(:,:),FSLP(:,:)


!
      REAL, ALLOCATABLE :: &
     &   CNUMWE(:),CNUMWO(:),CNUMUE(:),CNUMUO(:), &
     &   CMERWE(:,:),CMERWO(:,:),CMERUE(:,:),CMERUO(:,:)
!


!mz_ap_20071126+
#ifdef MESSY
! for chemistry coupling. It is changed in messy_mpiom_e5.f90!
      REAL, ALLOCATABLE :: swsum(:,:),swsumi(:,:), swrab(:,:,:)
#endif
!mz_ap_20071126-


#ifdef GMBOLUS
      REAL, ALLOCATABLE :: WGO(:,:,:),BOLX(:,:),BOLY(:,:)
#endif
      CONTAINS

      SUBROUTINE alloc_mem_commo1

      USE MO_PARAM1

      ALLOCATE (DZ(KE+1),DZW(KE),TIESTU(KE+1),TIESTW(KE+1),DI(KE+1))
!
!:: GRID DESCRIPTORS AND RELATED FIELDS
      ALLOCATE( DLXP(IE,JE),DLYP(IE,JE),DLXU(IE,JE),DLYV(IE,JE),        &
     &  DLXV(IE,JE),DLYU(IE,JE),DTDXUO(IE,JE),                          &
     &  DTDXUE(IE,JE),DTDYO(IE,JE),DTDXPO(IE,JE),DTDXPE(IE,JE),         &
     &  DPYE(IE,JE),DPYO(IE,JE),                                        &
     &  DDUE(IE,JE,KE),DDUO(IE,JE,KE),DDPO(IE,JE,KE),DPIO(IE,JE,KE),    &
     &  DDPSIO(IE,JE,KE),DWI(KE),                                       &
     &  DEUTE(IE,JE),DEUTO(IE,JE),DEUTIE(IE,JE),DEUTIO(IE,JE),          &
     &  DEPTO(IE,JE) )
!
!:: GEOGRAPHIC POSITION:
      ALLOCATE(GILA(ITO,JTO),GIPH(ITO,JTO),ALAT(IE,JE))
      ALLOCATE(GILA_G(2*IE_G,2*JE_G),GIPH_G(2*IE_G,2*JE_G))
!:: LAND SEA MASKS ON SKALAR/VECTOR POINTS
      ALLOCATE(WETO(IE,JE,KE),AMSUE(IE,JE,KE),AMSUO(IE,JE,KE))

      ALLOCATE(WETO_G(IE_G,JE_G,KE))
!
!:: ONE-DIMENSIONAL FIELDS
      ALLOCATE( DOTRA(KE),DOTRT(KE),UPTRT(KE),WINTUR(KEP),              &
     &   DBACKV(KEP),ABACKV(KEP),TAF(KE),SAF(KE),PREFF(KE),             &
     &   UPTRA(KE),                                                     &
     &   ROREF(KE),SKAL(ILL),B(ILT),X(ILT))
!
!:: THREE-DIMENSIONAL FIELDS
      ALLOCATE( RHOO(IE,JE,KE),UKO(IE,JE,KE),VKE(IE,JE,KE),             &
     &   THO(IE,JE,KE),SAO(IE,JE,KE),WO(IE,JE,KEP),                     &
     &   PO(IE,JE,KE),VK1E(IE,JE,KE),                                   &
     &   UK1O(IE,JE,KE),VOE(IE,JE,KE),UOO(IE,JE,KE),                    &
     &   T1O(IE,JE,KE),S1O(IE,JE,KE),DVO(IE,JE,KEP),AVO(IE,JE,KEP),     &
     &   STABIO(IE,JE,KE))
!
!:: TWO-DIMENSIONAL FIELDS
      ALLOCATE( ZO(IE,JE),Z1E(IE,JE),Z1O(IE,JE),                        &
     & U1E(IE,JE),V1E(IE,JE),U1O(IE,JE),V1O(IE,JE),        &
     & UROT(IE,JE),UDIV(IE,JE),USNO(IE,JE),VSNE(IE,JE),                 &
     & USO(IE,JE),VSE(IE,JE),UCOS(IE,JE),VCOS(IE,JE),                   &
     & SHELP(IE,JE),THELP(IE,JE),RHELP(IE,JE),                          &
     & TLOW(IE,JE),SLOW(IE,JE),TSUP(IE,JE),SSUP(IE,JE),                 &
     & PXOIN(IE,JE),PYEIN(IE,JE),UZO(IE,JE),VZE(IE,JE),                 &
     & FTWOU(IE,JE),FTWOV(IE,JE),CURVAV(IE,JE),B1E(IE,JE),              &
     & B1O(IE,JE),NUM(IE,JE)                                &
     &,KCONDEP(IE,JE),KBOT(IE,JE))

      ALLOCATE(NUM_G(IE_G,JE_G))
!
!:: ICE-MODEL FIELDS
       ALLOCATE(SICTHO(IE,JE),SICOMO(IE,JE),SICDIO(IE,JE),SICSHO(IE,JE) &
               ,SICPO(IE,JE),  SICOMP(IE,JE)                            &
               ,SICUO(IE,JE),SICVE(IE,JE),TICE(IE,JE)                   &
               ,HIBDELO(IE,JE),HIBZETO(IE,JE),HIBETO(IE,JE)             &
               ,HIBDELE(IE,JE),HIBZETE(IE,JE),HIBETE(IE,JE)             &
               ,TAUWATU(IE,JE),TAUWATV(IE,JE)                           &
               ,SICSNO(IE,JE),FWO(IE,JE),SHO(IE,JE))
               

!  FIELDS FOR COMMON DIAGNOSTICS
      ALLOCATE( HFLM(IE,JE)                 &
     &,PMEM(IE,JE),EISTRX(IE,JE),EISTRY(IE,JE)                          &
     &,AMLD(IE,JE,12),PSIUWE(IE,JE))
!
!:: FORCING FIELDS
      ALLOCATE( TXO(IE,JE),TYE(IE,JE),TAFO(IE,JE),                      &
     & HEATO(IE,JE),EMINPO(IE,JE),RELSAO(IE,JE),RELTHO(IE,JE),          &
     & FCLOU(IE,JE),FSWR(IE,JE),FU10(IE,JE)                             &
     &,FPREC(IE,JE),FTDEW(IE,JE),RIVRUN(IE,JE))
#ifdef RIVER_GIRIV
      ALLOCATE(GIRIV(IE,JE))
#endif /* RIVER_GIRIV*/

!
      ALLOCATE( CNUMWE(KE),CNUMWO(KE),CNUMUE(KE),CNUMUO(KE),            &
     & CMERWE(JE,KE),CMERWO(JE,KE),CMERUE(JE,KE),CMERUO(JE,KE))
!
!mz_ap_20071126+
#ifdef MESSY
      ALLOCATE(swsum(ie,je),swsumi(ie,je), swrab(ie,je,ke) )
#endif

!mz_ap_20071126-
#ifdef GMBOLUS
      ALLOCATE(WGO(IE,JE,KEP),BOLX(IE,JE),BOLY(IE,JE))
#endif
!<===END COMMO1
      END SUBROUTINE alloc_mem_commo1

      END MODULE MO_COMMO1
