  SUBROUTINE ml_ocean ( kproma                                         &
                  , pslm                                               &
                  , lonorth                                            &
                  , pseaice,    psiced,    palake                      &
                  , ptsi,       ptsw                                   &
                  , pahflw,     pahfsw,    pfluxres                    &
                  , ptrflw,     psoflw                                 &
                  , pamlcorr,   pamlcorac, pamlheatac                  &
                  , pevapi,     psni,      pcvsi                       &
                  , pahfres,    pfri                      )  

!
!  ---------------------------------------------------------------------
!
  USE mo_kind,           ONLY: dp                   ! op_pj_20110224
  USE mo_constants,      ONLY: stbo, rhoh2o, alf
  USE mo_physc2,         ONLY: ctfreez
  USE mo_radiation
  USE mo_time_control,   ONLY: delta_time
!
! Arguments
!
  INTEGER, INTENT(in) :: kproma                     ! op_pj_20110224
  REAL(dp) ::                                                                &
          pseaice(kproma),     psiced(kproma),     palake(kproma)            &
         ,ptsi(kproma),        ptsw(kproma)                                  &
         ,pahflw(kproma),      pahfsw(kproma),     pfluxres(kproma)          &
         ,ptrflw(kproma),      psoflw(kproma)                                &
         ,pamlcorr(kproma),    pamlcorac(kproma),  pamlheatac(kproma)        &
         ,pevapi(kproma),      psni(kproma),       pcvsi(kproma)             &
         ,pahfres(kproma),     pfri(kproma),       pslm(kproma)
  LOGICAL :: lonorth(kproma)      ! .true. for northern latitude

! Local variables

  REAL    :: zdtime, zfw, zheat, zfbase

! Executable statements
!
! 1. Set up constants
!
  zdtime = delta_time
  zalpha=2.1656_dp
  zalphas=0.31_dp
  zrho_sea=1025._dp
  zrho_sn=330._dp
  ziscond=zalpha/zalphas*zrho_sea/zrho_sn
  zcpice=2106._dp
  zrhoice=910._dp
  zdice=0.10_dp
  zrhoilf=zrhoice*alf
  zicecap=zrhoice*zcpice*zdice
  zdtrilf=zdtime/zrhoilf
  zfreez=-zdice/zdtrilf
  zdmix=50._dp
  zcpwater=3994._dp
  zrho_sea=1025._dp
  zmixcap=zrho_sea*zcpwater*zdmix
  zmcapdt=zdtime/zmixcap
  zmcaprilf=zmixcap/zrhoilf
!
! 2. Mixed layer ocean temperature and ice thickness
!
  DO jl=1,kproma
!
  IF (palake(jl).EQ.0. .AND. pslm(jl).lt.0.5) THEN      ! no lake points
!
     IF (pseaice(jl).LT.0.5) THEN                      ! open water
!
        zfluxw=pahflw(jl)+pahfsw(jl)+ptrflw(jl)+psoflw(jl)-pamlcorr(jl)
!
!       Water temperature (ptsw)
!
        zts=ptsw(jl)+zmcapdt*(zfluxw+pfluxres(jl))
        IF(zts.LT.ctfreez) THEN
          pamlcorr(jl)=MIN(0.0,pamlcorr(jl))
          zfluxw=pahflw(jl)+pahfsw(jl)+ptrflw(jl)+psoflw(jl)-pamlcorr(jl)
          zts=ptsw(jl)+zmcapdt*(zfluxw+pfluxres(jl))
        END IF
        ptsi(jl)=ctfreez
        pfluxres(jl)=0.
        psiced(jl)=0.
        IF (zts.GE.ctfreez) THEN                  ! open water (unchanged)
           ptsw(jl)=zts
        ELSE                                    ! check ice formation
           ptsw(jl)=ctfreez
           zfres=(zts-ctfreez)/zmcapdt            ! < 0.
           IF (zfres.LE.zfreez) THEN            ! ice formation
              psiced(jl)=zmcaprilf*(ctfreez-zts)  ! > zdice
              pseaice(jl)=1.
           ELSE
              pfluxres(jl)=zfres
           END IF
        END IF
!  ---------------------------------------------------------------------
     ELSE IF (psiced(jl).GE.zdice) THEN
!
        IF (lonorth(jl)) THEN
          zfbase=0.0
        ELSE
          zfbase=20.0
        END IF
!
!       Ice thickness (psiced)
!
        zconhflx=zalpha*(ptsi(jl)-ctfreez)/(psiced(jl)+ziscond*psni(jl))
        zsubice=(1.-pcvsi(jl))*pevapi(jl)*zdtime/zrhoice
        pamlcorr(jl)=MIN(0.0,pamlcorr(jl))-zfbase
        zhi=psiced(jl)-zdtrilf*(zconhflx+pfluxres(jl)-pamlcorr(jl))+zsubice
        ptsw(jl)=ctfreez
        IF (zhi.GE.zdice) THEN
           psiced(jl)=zhi
           pseaice(jl)=1.
           pfluxres(jl)=0.
        ELSE IF (zhi.LE.0.) THEN               ! complete melting
           ptsw(jl)=ctfreez-zhi/zmcaprilf        ! ptsw > ctfreez
           pfluxres(jl)=-rhoh2o*alf*psni(jl)/zdtime
           psiced(jl)=0.
           psni(jl)=0.
           pseaice(jl)=0.
        ELSE                                   ! incomplete melting
           psiced(jl)=zdice
           pseaice(jl)=1.
           pfluxres(jl)=(zdice-zhi)/zdtrilf
           pahfres(jl)=pahfres(jl)-zdtime*pfri(jl)*pfluxres(jl)
        END IF
     END IF
     zfw=1.-pfri(jl)-pslm(jl)
! accumulate mixed layer variables
     zheat=zmixcap*zfw*ptsw(jl)-zrhoilf*pfri(jl)*psiced(jl)
     pamlheatac(jl)=pamlheatac(jl)+zheat*zdtime
     pamlcorac(jl)=pamlcorac(jl)+pamlcorr(jl)*zdtime
   END IF
  END DO
!  ---------------------------------------------------------------------
     RETURN
  END SUBROUTINE ml_ocean
