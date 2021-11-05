MODULE mo_ocice
 !     hiblers ice model as described in
 !     w.d.hibler iii, a dynamic thermodynamic sea ice model. j. phys. oceanogr.
 !     9, 815 - 846, 1979

 !
 !     modifications:
 !     uwe 2.2.00
 !       include advection of ice velocities
 !       correct tau-w term using mixture of old and new velocities
 !        new fileds included:
 !         speedu    difference in speed between water and ice (old)
 !         speedv
 !     uwe 8.3.00
 !        include tauwat, proper determination of water ice stress
 !        to be used in ocwind
 !     uwe 28.6.00
 !        make sw-penetration applicable under sea ice too
 !     uwe 10/00
 !        include implicite treatment of mass advection in iteration
 !         of ice velocities
  USE mo_param1
  USE mo_planetary_constants, ONLY: g, rhoref_water, rhoref_ice            &
       , rhoref_snow , inv_rhoref_water, rhoicwa, rhosnic, rhosnwa, slpref
  USE mo_parallel
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1
  USE mo_commoau1
  USE mo_commoau2
  USE mo_units
  USE mo_swr_absorption, ONLY: dynamic_swr_absorption, heatabs, &
       lfb_bgc_oce, swsum, subsurface_swr_absorption, swr_frac_from_chl
  USE mo_tro
  USE mo_tidal, ONLY : sal, tipoto
  USE mo_contro, ONLY : contro
#ifdef PBGC
  USE mo_param1_bgc , only: nocetra
  USE mo_carbch , only: ocetra
#endif

#ifdef __coupled
  USE mo_fluxes1
#endif

  IMPLICIT NONE
  ! compute salt transport correction
  LOGICAL :: lsaoclose = .FALSE.
  REAL(wp), PARAMETER :: cw = 0.0045_wp
  !> the coefficients e and p from (7), (8)
  REAL(wp), PARAMETER :: hicce = 2._wp, hiccp = 20._wp

# ifdef MESSY
INTERFACE ocice
  MODULE PROCEDURE ocice
  MODULE PROCEDURE ocice_coupled
END INTERFACE
INTERFACE ice_thermodynamics
  MODULE PROCEDURE ice_thermodynamics
  MODULE PROCEDURE ice_thermodynamics_coupled
END INTERFACE
#endif /*MESSY*/ 


CONTAINS

SUBROUTINE ocice

  INTEGER i,j
  REAL(wp) sicomo1
  REAL(wp) siotho(ie,je),sioomo(ie,je),siosno(ie,je)
  !     taf atmospheric temperature
  !     sicth  ice thickness
  !     hiccp is divided by mean density of 1000 kg / m**3
  !     sicom  ice compactness
  !     sicdi  ice flow divergence
  !     sicsh  ice flow shear
  !     hibzet,hibdel,hibet are the fields of (7) - (9)
  !
  !      part 4: increase of existing ice
  ! in theory a distinction should be made bewteen the ice-covered part of a
  ! grid-cell for which this do-loop applies and the ice free part, which
  ! should be treated as do loops 20 and 25. such a distinction enhances the
  ! overall ice-growth which results in irrealistically high values in some
  ! grid points and subsequent numerical instabilities. inclusion of diffusion
  ! on icethickness and compactness also gave rise to instabilities.
  !
  !      cutoff value for prognostic calculation of ice velocities
  !

  ! Need halos for sea ice variables for copying
  CALL bounds_exch(1,'p',sicomo,'mo_ocice 100')
  CALL bounds_exch(1,'p',sicsno,'mo_ocice 101')
  CALL bounds_exch(1,'p',sictho,'mo_ocice 102')

!$OMP PARALLEL  private(i,j,sicomo1)

!$OMP DO
  ! Normalize sea ice variables
  DO j=1, je
    DO i=1, ie
      sicomo1 = sicomo(i, j)
      sictho(i, j) = MAX(0._wp, sictho(i, j) * weto(i, j, 1))
      IF(sictho(i, j) <= 0._wp) sicomo1 = 0._wp
      sicomo1 = MAX(0._wp, sicomo1*weto(i, j, 1))
      sicsno(i, j) = MAX(0._wp, sicsno(i, j) * weto(i, j, 1))
      sicomo(i, j) = MIN(1._wp, sicomo1)
    END DO
  END DO
!$OMP END DO

!$OMP DO
  ! store old values for sr growth
  DO j=1,je
    DO i=1,ie
      siotho(i,j)=sictho(i,j)
      sioomo(i,j)=sicomo(i,j)
      siosno(i,j)=sicsno(i,j)
    END DO
  END DO
!$OMP END DO

! @remove Moved bounds exchange into ice_dynamics

!$OMP END PARALLEL

  ! Run sea-ice related computations
  IF (icontro .NE. 0) CALL contro(-17)
  CALL ice_dynamics
  IF (icontro .NE. 0) CALL contro(-18)
  ! Store advection values into new model.
  CALL ice_advection
  IF (icontro .NE. 0) CALL contro(-19)
  ! Compute ice thermodynamics with old model values,
  ! use new model values for budget computations.
  CALL ice_thermodynamics(siotho, sioomo, siosno)
  IF (icontro .NE. 0) CALL contro(-20)
END SUBROUTINE ocice

#ifdef MESSY
SUBROUTINE ocice_coupled(AOFLFRIO,AOFLFRWO,AOFLRHIO,AOFLCHIO,AOFLNHWO,AOFLSHWO,AOFLWSVO)
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLFRIO!solid freshwater flux (over ice only) 
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLFRWO!liquid freshwater flux (water and ice)
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLRHIO!residual heat flux (sea-ice topmelt heat flux)
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLCHIO!conductive heat flux through ice 
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLNHWO!net heat flux over water
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLSHWO!downwelling solar radiation
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLWSVO!wind stress velocity 

  INTEGER i,j
  REAL(wp) sicomo1
  REAL(wp) siotho(ie,je),sioomo(ie,je),siosno(ie,je)
  !     taf atmospheric temperature
  !     sicth  ice thickness
  !     hiccp is divided by mean density of 1000 kg / m**3
  !     sicom  ice compactness
  !     sicdi  ice flow divergence
  !     sicsh  ice flow shear
  !     hibzet,hibdel,hibet are the fields of (7) - (9)
  !
  !      part 4: increase of existing ice
  ! in theory a distinction should be made bewteen the ice-covered part of a
  ! grid-cell for which this do-loop applies and the ice free part, which
  ! should be treated as do loops 20 and 25. such a distinction enhances the
  ! overall ice-growth which results in irrealistically high values in some
  ! grid points and subsequent numerical instabilities. inclusion of diffusion
  ! on icethickness and compactness also gave rise to instabilities.
  !
  !      cutoff value for prognostic calculation of ice velocities
  !

  ! Need halos for sea ice variables for copying
  CALL bounds_exch(1,'p',sicomo,'mo_ocice 100')
  CALL bounds_exch(1,'p',sicsno,'mo_ocice 101')
  CALL bounds_exch(1,'p',sictho,'mo_ocice 102')

!$OMP PARALLEL  private(i,j,sicomo1)

!$OMP DO
  ! Normalize sea ice variables
  DO j=1, je
    DO i=1, ie
      sicomo1 = sicomo(i, j)
      sictho(i, j) = MAX(0._wp, sictho(i, j) * weto(i, j, 1))
      IF(sictho(i, j) <= 0._wp) sicomo1 = 0._wp
      sicomo1 = MAX(0._wp, sicomo1*weto(i, j, 1))
      sicsno(i, j) = MAX(0._wp, sicsno(i, j) * weto(i, j, 1))
      sicomo(i, j) = MIN(1._wp, sicomo1)
    END DO
  END DO
!$OMP END DO

!$OMP DO
  ! store old values for sr growth
  DO j=1,je
    DO i=1,ie
      siotho(i,j)=sictho(i,j)
      sioomo(i,j)=sicomo(i,j)
      siosno(i,j)=sicsno(i,j)
    END DO
  END DO
!$OMP END DO

! @remove Moved bounds exchange into ice_dynamics

!$OMP END PARALLEL

  ! Run sea-ice related computations
  IF (icontro .NE. 0) CALL contro(-17)
  CALL ice_dynamics
  IF (icontro .NE. 0) CALL contro(-18)
  ! Store advection values into new model.
  CALL ice_advection
  IF (icontro .NE. 0) CALL contro(-19)
  ! Compute ice thermodynamics with old model values,
  ! use new model values for budget computations.
  CALL ice_thermodynamics(siotho, sioomo, siosno, AOFLFRIO,AOFLFRWO,AOFLRHIO,AOFLCHIO,AOFLNHWO,AOFLSHWO,AOFLWSVO)
  IF (icontro .NE. 0) CALL contro(-20)
END SUBROUTINE ocice_coupled
#endif

SUBROUTINE ice_dynamics

  INTEGER i,j,iter

  REAL(wp) rhoicsn,rhowaic,uepsi
  REAL(wp) eps11(ie,je),eps22(ie,je),eps12(ie,je)
  REAL(wp) uh(ie,je),vh(ie,je)
  REAL(wp) ux,uy,vx,vy
  REAL(wp) e12,e11,e22,hibcc,pst,rst,sicm,sicom,dxdx,dydy,dxdy
  REAL(wp) alpalt,alpneu
  REAL(wp) siouo(ie,je),siove(ie,je)
  REAL(wp) effico(ie,je),effice(ie,je)
  REAL(wp) cweffo(ie,je),cweffe(ie,je)
  REAL(wp) zsic(ie,je)
  REAL(wp) zschaltd
  REAL(wp) zschalt
  REAL(wp) ucor(ie,je),vcor(ie,je)
  REAL(wp) zdiff(ie,je)
  REAL(wp) zzl(ie,je),zzr(ie,je),zzo(ie,je),zzu(ie,je)
  REAL(wp) rs,rv,speedmi,reval
  REAL(wp) speedu(ie,je),speedv(ie,je)

  reval = 0.01_wp
  speedmi = 0.01_wp
  rhowaic=rhoref_water/rhoref_ice
  rhoicsn=rhoref_ice/rhoref_snow
  uepsi = 1.e-8_wp

!$OMP PARALLEL  private(i,j,ux,uy,vx,vy,e12,e11,e22,hibcc, &
!$OMP   pst,rst,sicm,sicom,dxdx,dydy,dxdy,rs,rv)

!$OMP SINGLE
  CALL bounds_exch(1,'u',sicuo,'mo_ocice 103')
  CALL bounds_exch(1,'vf',sicve,'mo_ocice 104')
!$OMP END SINGLE

! DN: Commenting out this block gives binary different results. Why??

  ! Also copy boundaries -> no bounds exchange needed.
!$OMP DO
  DO j=1,je
     DO i=1,ie
        sicomp(i,j)=sicomo(i,j)
    END DO
  END DO
!$OMP END DO

!$OMP DO
      DO J=1,JE
       DO I=1,IE
         SICUDO(I,J)=SICUO(I,J)
         SICVDE(I,J)=SICVE(I,J)
       ENDDO
      ENDDO
!$OMP END DO



!      ice dynamics

!$OMP DO
  DO j=1,je
    DO i=1,ie
      eps11(i, j) = 0._wp
      eps12(i, j) = 0._wp
      eps22(i, j) = 0._wp
      uh(i, j) = 0._wp
      vh(i, j) = 0._wp
      zdiff(i, j) = 0._wp
      effico(i, j) = 0._wp
      effice(i, j) = 0._wp
      speedu(i, j) = 0._wp
      speedv(i, j) = 0._wp
      zzl(i, j) = 0._wp
      zzr(i, j) = 0._wp
      zzo(i, j) = 0._wp
      zzu(i, j) = 0._wp
      ucor(i, j) = 0._wp
      vcor(i, j) = 0._wp
      siouo(i, j) = 0._wp
      siove(i, j) = 0._wp
   END DO
  END DO
!$OMP END DO

!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      ux=(sicuo(i,j)-sicuo(i-1,j))/dlxp(i,j)
      vy=(sicve(i,j-1)-sicve(i,j))/dlyp(i,j)
      vx=(sicve(i+1,j)-sicve(i,j))/dlxp(i,j)
      uy=(sicuo(i,j)-sicuo(i,j+1))/dlyp(i,j)
      eps11(i,j)=ux
      eps22(i,j)=vy
      eps12(i, j) = 0.5_wp * (vx + uy)
      IF(amsuo(i, j, 1) .GT. 0.5_wp)                                         &
           speedu(i, j)                                                      &
           = MAX(SQRT((uepsi + sicuo(i, j) - uko(i, j, 1))**2                &
           &          + (0.25_wp * (sicve(i, j) + sicve(i+1, j)              &
           &                        + sicve(i, j-1) + sicve(i+1, j-1)        &
           &                        - vke(i, j, 1) - vke(i+1, j, 1)          &
           &                        - vke(i,j-1,1) - vke(i+1, j-1, 1)))**2)  &
           , speedmi)
      IF (amsue(i, j, 1) .GT. 0.5_wp)                                        &
           speedv(i, j)                                                      &
           = MAX(SQRT((uepsi + sicve(i, j) - vke(i, j, 1))**2                &
           &          + (0.25_wp * (sicuo(i, j) + sicuo(i-1, j)              &
           &                        + sicuo(i, j+1) + sicuo(i-1,j+1)         &
           &                        - uko(i, j, 1) - uko(i-1, j, 1)          &
           &                        - uko(i, j+1, 1) - uko(i-1, j+1, 1)))**2)&
           ,speedmi)
    END DO
  END DO
!$OMP END DO


!$OMP SINGLE
!#ifdef bounds_exch_save
  CALL bounds_exch(1,'p',eps11,'mo_ocice 106')
  CALL bounds_exch(1,'p',eps22,'mo_ocice 107')
  CALL bounds_exch(1,'s',eps12,'mo_ocice 108')
!#endif
  CALL bounds_exch(1,'u',speedu,'mo_ocice 109')
  CALL bounds_exch(1,'v',speedv,'mo_ocice 110')
!$OMP END SINGLE

!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      e12 = 0.25_wp * (eps12(i-1, j-1) + eps12(i, j-1) &
           &           + eps12(i-1, j) + eps12(i, j))
      !      argu=((eps11(i,j)-eps22(i,j))**2+4.*e12**2)/hicce**2             &
      !    & +(eps11(i,j)+eps22(i,j))**2
      hibdelo(i, j) &
           = SQRT(((eps11(i, j) - eps22(i, j))**2 + 4._wp * e12**2)/hicce**2  &
           &      + (eps11(i, j) + eps22(i, j))**2)
    END DO
  END DO
!$OMP END DO

!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      e11 = 0.25_wp * (eps11(i, j) + eps11(i+1, j) + eps11(i, j+1) &
           &           + eps11(i+1,j+1))
      e22 = 0.25_wp * (eps22(i, j) + eps22(i+1, j) + eps22(i,j+1) &
           &           + eps22(i+1, j+1))
      hibdele(i, j) = SQRT(((e11 - e22)**2 + 4._wp * eps12(i, j)**2)/hicce**2 &
           &               + (e11 + e22)**2)
    END DO
  END DO
!$OMP END DO
!#ifdef bounds_exch_save
!$OMP SINGLE
   CALL bounds_exch(1,'p',hibdelo,'mo_ocice 111')
   CALL bounds_exch(1,'s',hibdele,'mo_ocice 112a')
!$OMP END SINGLE
!#endif
!uwe introduce hibcc, corresponds to c IN hibler, jpo 9,825

  hibcc = 20._wp

  alpneu = 1._wp

  alpalt = 1._wp - alpneu

!$OMP DO
  DO j=2,je1
    DO i=2,ie1

!uwe    hibler (17),
      ! FIXME: why not almzer instead of 1.e-20?
      pst = hiccp * sictho(i, j) * EXP(-hibcc * (1._wp - sicomo(i, j))) &
           * weto(i, j, 1)
      sicm=(sictho(i,j)+sictho(i+1,j)+sictho(i,j+1)+sictho(i+1,j+1))            &
           /(weto(i,j,1)+weto(i+1,j,1)+weto(i,j+1,1)+weto(i+1,j+1,1)+1.e-20_wp)
      sicom=(sicomo(i,j)+sicomo(i+1,j)+sicomo(i,j+1)+sicomo(i+1,j+1))           &
           /(weto(i,j,1)+weto(i+1,j,1)+weto(i,j+1,1)+weto(i+1,j+1,1)+1.e-20_wp)
      sicom = MIN(sicom, 1._wp)
      sicom = MAX(sicom, 0._wp)
      sicm = MAX(sicm, 0._wp)
      rst = hiccp * sicm * EXP(-hibcc*(1._wp - sicom))
      uh(i,j)=sicuo(i,j)
      vh(i,j)=sicve(i,j)
      siouo(i,j)=sicuo(i,j)
      siove(i,j)=sicve(i,j)
      ! FIXME: this makes not sense, hibdelo(i,j) is already >almzer
      ! because of the predicate
      IF(hibdelo(i,j).GT.almzer)THEN
        hibzeto(i, j) = alpneu * pst * 0.5_wp / MAX(hibdelo(i, j), almzer) &
             +alpalt*hibzeto(i,j)
        hibeto(i,j)=alpneu*hibzeto(i,j)/hicce**2+alpalt*hibeto(i,j)
      ELSE
        hibzeto(i, j) = 0._wp
        hibeto(i, j) = 0._wp
      ENDIF
      !FIXME: this makes no sense, hibdele(i,j) is already >almzer
      IF(hibdele(i,j).GT.almzer) THEN
        hibzete(i, j) = alpneu * rst * 0.5_wp/MAX(hibdele(i,j),almzer)           &
             +alpalt*hibzete(i,j)
        hibete(i,j)=alpneu*hibzete(i,j)/hicce**2+alpalt*hibete(i,j)
      ELSE
        hibzete(i,j) = 0._wp
        hibete(i,j) = 0._wp
      ENDIF
      !c       spur=eps11(i,j)+eps22(i,j)
      !
      effico(i, j) = 0.5_wp * (sictho(i, j) + sictho(i+1, j)           &
           +rhoicsn*(sicsno(i,j)+sicsno(i+1,j)))*amsuo(i,j,1)
      effice(i, j) = 0.5_wp * (sictho(i, j) + sictho(i, j+1)           &
           +rhoicsn*(sicsno(i,j)+sicsno(i,j+1)))*amsue(i,j,1)
    ENDDO
  ENDDO
!$OMP END DO

!$OMP SINGLE
!#ifdef bounds_exch_save
  CALL bounds_exch(1,'u',uh,'mo_ocice 112')
  CALL bounds_exch(1,'v',vh,'mo_ocice 113')
  CALL bounds_exch(1,'u',siouo,'mo_ocice 114')
  CALL bounds_exch(1,'v',siove,'mo_ocice 115')
  CALL bounds_exch(1,'p',hibzeto,'mo_ocice 116')
  CALL bounds_exch(1,'p',hibeto,'mo_ocice 117')
  CALL bounds_exch(1,'s',hibzete,'mo_ocice 118')
  CALL bounds_exch(1,'s',hibete,'mo_ocice 119')
  CALL bounds_exch(1,'u',effico,'mo_ocice 120')
  CALL bounds_exch(1,'v',effice,'mo_ocice 121')
!#endif
!$OMP END SINGLE


!      enhance friction coefficient for very thin ice
!      ==> smoothe transition towards water velocities for effic ==> 0.

!$OMP DO
  DO j=1,je
    DO i=1,ie
      !c         cweffo(i,j)=cw*MAX(1.,1./(25.*(effico(i,j)+1.e-6)**2))
      !c         cweffe(i,j)=cw*MAX(1.,1./(25.*(effice(i,j)+1.e-6)**2))
      cweffo(i, j) = cw &
           * MAX(1._wp, 1._wp / (10._wp * effico(i, j) + 1.e-6_wp))**1
      cweffe(i, j) = cw &
           * MAX(1._wp, 1._wp / (10._wp * effice(i, j) + 1.e-6_wp))**1
      sicudo(i, j) &
           = (0.8_wp * sicudo(i, j) + 0.2_wp * sicuo(i, j)) * amsuo(i, j, 1)
      sicvde(i, j) &
           = (0.8_wp * sicvde(i, j) + 0.2_wp * sicve(i, j)) * amsue(i, j, 1)
      zsic(i, j) = sictho(i, j)+rhosnic*sicsno(i, j)
    END DO
  END DO
!$OMP END DO

  !      old ice velocities, u on v-point and v on u-point
  !
  !      zschaltd: upwind type diffusion

  zschaltd = 1._wp
!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      ucor(i, j) = 0.25_wp * (uh(i, j) + uh(i-1, j) + uh(i, j+1) + uh(i-1, j+1))
      vcor(i, j) = 0.25_wp * (vh(i, j) + vh(i+1, j) + vh(i, j-1) + vh(i+1, j-1))
      !
      !     sea level change due to diffusion of sea ice
      !
      zdiff(i, j) = zschaltd * 0.5_wp * dt * rhoicwa *                   &
           (ABS(sicudo(i-1,j))*dlyu(i-1,j)*(zsic(i-1,j)-zsic(i,j))       &
           +ABS(sicudo(i,j))*dlyu(i,j)*(zsic(i+1,j)-zsic(i,j))           &
           +ABS(sicvde(i,j))*dlxv(i,j)*(zsic(i,j+1)-zsic(i,j))           &
           +ABS(sicvde(i,j-1))*dlxv(i,j-1)*(zsic(i,j-1)-zsic(i,j)))      &
           /(dlxp(i,j)*dlyp(i,j))

      IF ( ltidal ) THEN
        zdiff(i,j)=zdiff(i,j)       +sal*tipoto(i,j)/g
      ENDIF

      zdiff(i,j)=zdiff(i,j) + (fslp(i,j)-slpref)/(rhoref_water*g)

!     sea ice thickness protection to avaid negative layer thickness
      zdiff(i,j) = zdiff(i,j) + MAX(0._wp, &
           zsic(i,j) - 0.75_wp * (dzw(1) + zo(i,j)))**2



      !
      !     thickness protection, avoid negative layer thickness
      !       induce divergent ice transports
      !
      !c     x  +  MAX(0.,zsic(i,j)-0.6667*(dzw(1)+zo(i,j)))
      !c     x     /MAX(0.01,dzw(1)+zo(i,j)-zsic(i,j))
      !
    END DO
  END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'u',ucor,'mo_ocice 122')
  CALL bounds_exch(1,'v',vcor,'mo_ocice 123')
  CALL bounds_exch(1,'p',zdiff,'mo_ocice 124')
!$OMP END SINGLE
!#endif

  zschalt = 1._wp
  !      main iteration

  !      iteration_loop:

DO iter = 1, 40

!$OMP DO
    DO j=1,je
      DO i=1,ie
        uh(i, j) = 0._wp
        vh(i, j) = 0._wp
      END DO
    END DO
!$OMP END DO

!$OMP DO
    DO j=2,je1
      DO i=2,ie1
        zzl(i,j)=dt*rhoicwa*(                                            &
             dlyu(i-1,j)*sicuo(i-1,j)*effico(i-1,j)                       &
             +dlxv(i,j)*sicve(i,j)*effice(i,j)                             &
             -dlxv(i,j-1)*sicve(i,j-1)*effice(i,j-1))                      &
             /(dlxp(i,j)*dlyp(i,j))
        zzr(i,j)=dt*rhoicwa*(                                            &
             -dlyu(i+1,j)*sicuo(i+1,j)*effico(i+1,j)                       &
             +dlxv(i+1,j)*sicve(i+1,j)*effice(i+1,j)                       &
             -dlxv(i+1,j-1)*sicve(i+1,j-1)*effice(i+1,j-1))                &
             /(dlxp(i+1,j)*dlyp(i+1,j))
        speedu(i, j) = 0.5_wp * (MAX(SQRT((uepsi + sicuo(i, j) - uko(i, j, 1))**2 &
             +(0.25_wp * (sicve(i,j)+sicve(i+1,j)+sicve(i,j-1)+sicve(i+1,j-1) &
             -vke(i,j,1)-vke(i+1,j,1)-vke(i,j-1,1)-vke(i+1,j-1,1)))**2)   &
             ,speedmi)+speedu(i,j))
        speedv(i, j) = 0.5_wp * (MAX(SQRT((uepsi + sicve(i, j) - vke(i, j, 1))**2 &
             +(0.25_wp*(sicuo(i,j)+sicuo(i-1,j)+sicuo(i,j+1)+sicuo(i-1,j+1) &
             -uko(i,j,1)-uko(i-1,j,1)-uko(i,j+1,1)-uko(i-1,j+1,1)))**2)   &
             ,speedmi)+speedv(i,j))
      END DO
    END DO
!$OMP END DO

!$OMP SINGLE
!#ifdef bounds_exch_save
    CALL bounds_exch(1,'p',zzl,'mo_ocice 124a')
    CALL bounds_exch(1,'p',zzr,'mo_ocice 125')
!#endif
    CALL bounds_exch(1,'u',speedu,'mo_ocice 126')
    CALL bounds_exch(1,'v',speedv,'mo_ocice 127')
!$OMP END SINGLE

!$OMP DO
    DO j=2,je1
      DO i=2,ie1

        IF (effico(i,j).GT.reval)THEN
          dxdx=dlxp(i,j)**2
          dydy=dlyp(i,j)**2
          dxdy=dlxp(i,j)*dlyp(i,j)
          rs = effico(i, j) * (siouo(i, j) &
               + 0.5_wp * dt * ftwou(i, j) * (vcor(i,j)          &
               + 0.25_wp * (sicve(i,j)+sicve(i+1,j)+sicve(i,j-1)+sicve(i+1,j-1))))     &
               +cweffo(i,j)*dt*uko(i,j,1)*speedu(i,j)                           &
               + dt*hiccp*(sictho(i,j)*EXP(-hibcc*(1._wp-sicomo(i,j)))              &
               -sictho(i+1,j)*EXP(-hibcc*(1._wp-sicomo(i+1,j))))/dlxu(i,j)       &
               + dt*effico(i,j)*g*(zo(i,j)-zo(i+1,j)+zdiff(i,j)-zdiff(i+1,j)     &
               +zschalt*(zzl(i,j)-zzr(i,j)))/dlxu(i,j)                       &
               +dt*rhowaic &
#ifdef __coupled
               !sv 31.08.99 included fluxes option
               !   (i.e. distinguish bewteen tau over water and ice)
               !sv here only tau over ice is considered (i.e. txo is replaced by aofltxio)
               !sv tau over water is used IN sbr ocwind
               *aofltxio(i,j)*inv_rhoref_water &
#else
               *txo(i,j) &
#endif
               * (sicomo(i, j) + sicomo(i+1, j)) * 0.5_wp
               !svx 31.08.99

          rv=                                                              &
               +(hibzeto(i,j)+hibeto(i,j))    *sicuo(i-1,j)/dxdx                 &
               +(hibzeto(i+1,j)+hibeto(i+1,j))*sicuo(i+1,j)/dxdx                 &
               +(hibzeto(i+1,j)-hibeto(i+1,j))*(sicve(i+1,j-1)-sicve(i+1,j))     &
               /dxdy                                                            &
               -(hibzeto(i,j)-hibeto(i,j))*(sicve(i,j-1)-sicve(i,j))/dxdy        &
               +(hibete(i,j-1)*sicuo(i,j-1)+hibete(i,j)*sicuo(i,j+1))/dydy       &
               +hibete(i,j-1)*(sicve(i+1,j-1)-sicve(i,j-1))/dxdy                 &
               -hibete(i,j)*(sicve(i+1,j)-sicve(i,j))/dxdy

          uh(i,j)=(rs+rv*dt)/(effico(i,j)                                  &
               +cweffo(i,j)*dt*speedu(i,j)                                 &
               +dt*effico(i,j)*g*dt*rhoicwa*effico(i,j)*zschalt*              &
               (1._wp / (dlxp(i, j) * dlyp(i, j)) &
               & + 1._wp/(dlxp(i+1, j) * dlyp(i+1, j))) * dlyu(i, j)/dlxu(i, j)&
               +dt*(hibzeto(i,j)+hibeto(i,j)+hibzeto(i+1,j)+hibeto(i+1,j))/dxdx &
               +dt*(hibete(i,j-1)+hibete(i,j))/dydy)
        ENDIF

      END DO
    END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch(1,'u',uh,'mo_ocice 128')
!$OMP END SINGLE
!#endif

!$OMP DO
    DO j=2,je1
      DO i=2,ie1
        zzo(i,j)=dt*rhoicwa*(                                           &
             dlyu(i-1,j)*sicuo(i-1,j)*effico(i-1,j)                      &
             -dlyu(i,j)*sicuo(i,j)*effico(i,j)                            &
             -dlxv(i,j-1)*sicve(i,j-1)*effice(i,j-1))                     &
             /(dlxp(i,j)*dlyp(i,j))
        zzu(i,j)=dt*rhoicwa*(                                           &
             dlyu(i-1,j+1)*sicuo(i-1,j+1)*effico(i-1,j+1)                &
             -dlyu(i,j+1)*sicuo(i,j+1)*effico(i,j+1)                      &
             +dlxv(i,j+1)*sicve(i,j+1)*effice(i,j+1))                     &
             /(dlxp(i,j+1)*dlyp(i,j+1))
      END DO
    END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch(1,'p',zzo,'mo_ocice 129')
    CALL bounds_exch(1,'p',zzu,'mo_ocice 130')
!$OMP END SINGLE
!#endif

!$OMP DO
    DO j=2,je1
      DO i=2,ie1
        IF(effice(i,j).GT.reval) THEN
          dxdx=dlxp(i,j)**2
          dydy=dlyp(i,j)**2
          dxdy=dlxp(i,j)*dlyp(i,j)
          rs = effice(i, j) * (siove(i, j) &
               - 0.5_wp * dt * ftwov(i, j) *(ucor(i,j)           &
               + 0.25_wp * (sicuo(i, j) + sicuo(i-1, j) + sicuo(i, j+1) &
               &            + sicuo(i-1, j+1))))     &
               +cweffe(i,j)*dt*vke(i,j,1)*speedv(i,j)                           &
               + dt * hiccp * (sictho(i, j+1) &
               &    * EXP(-hibcc*(1._wp - sicomo(i, j+1))) &
               &    - sictho(i, j) * EXP(-hibcc * (1._wp - sicomo(i, j)))) &
               / dlyv(i, j)            &
               +effice(i,j)*g*dt*(zo(i,j+1)-zo(i,j)+zdiff(i,j+1)-zdiff(i,j)     &
               +zschalt*(zzu(i,j)-zzo(i,j)))/dlyv(i,j)                      &
               +rhowaic &
#ifdef __coupled
               !sv 31.08.99 included fluxes option
               !    (i.e. distinguish between tau over water and ice)
               !sv  here only tau over ice is considered (i.e. tye is replaced by aofltyie)
               !sv         tau over water is used IN sbr ocwind
               *aofltyie(i,j)*inv_rhoref_water &
#else
               *tye(i,j) &
#endif
               * dt * (sicomo(i, j) + sicomo(i, j+1)) * 0.5_wp
               !svx 31.08.99
          rv=(hibzeto(i,j)+hibeto(i,j))*sicve(i,j-1)/dydy                      &
               +   (hibzeto(i,j+1)+hibeto(i,j+1))*sicve(i,j+1)/dydy            &
               +(hibete(i,j)*sicve(i+1,j)+hibete(i-1,j)*sicve(i-1,j))/dxdx     &
               +(hibzeto(i,j)-hibeto(i,j))*(sicuo(i,j)-sicuo(i-1,j))/dxdy      &
               -(hibzeto(i,j+1)-hibeto(i,j+1))*(sicuo(i,j+1)-sicuo(i-1,j+1))   &
               /dxdy                                                           &
               +hibete(i,j)*(sicuo(i,j)-sicuo(i,j+1))/dxdy                     &
               -hibete(i-1,j)*(sicuo(i-1,j)-sicuo(i-1,j+1))/dxdy
          !
          vh(i,j)=(rs+rv*dt)/(effice(i,j)                                      &
               +cweffe(i,j)*dt*speedv(i,j)                                     &
               +dt*effice(i,j)*g*dt*zschalt*effice(i,j)*rhoicwa*               &
               (1._wp/(dlxp(i,j)*dlyp(i,j))+1._wp/(dlxp(i,j+1)*dlyp(i,j+1)))   &
               *dlxv(i,j)/dlyv(i,j)                                            &
               +dt*(hibzeto(i,j)+hibeto(i,j)+hibzeto(i,j+1)+hibeto(i,j+1))/dydy&
               +dt*(hibete(i,j)+hibete(i-1,j))/dxdx)
          !
        ENDIF
      END DO
    END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch(1,'v',vh,'mo_ocice 131')
!$OMP END SINGLE
!#endif

!$OMP DO
    DO j=2,je1
      DO i=2,ie1
        IF(effico(i,j).GT.reval)THEN
          sicuo(i, j) = 0.5_wp * (uh(i, j) + sicuo(i, j)) * amsuo(i, j, 1)
        ELSE
          sicuo(i,j)=uko(i,j,1)
        ENDIF
        IF(effice(i,j).GT.reval)THEN
          sicve(i, j) = 0.5_wp * (vh(i, j) + sicve(i, j)) * amsue(i, j, 1)
        ELSE
          sicve(i,j)=vke(i,j,1)
        ENDIF
      END DO
    END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch(1,'u',sicuo,'mo_ocice 132')
    CALL bounds_exch(1,'vf',sicve,'mo_ocice 133')
!$OMP END SINGLE
!#endif

    !      write(io_stdout,*) ' iteration   u ',iter
    !      write(io_stdout,681)((amsuo(i,j,1)*sicuo(i,j),j=1,25),i=1,ie)
    !      write(io_stdout,*) ' iteration   v ',iter
    !      write(io_stdout,681)((amsue(i,j,1)*sicve(i,j),j=1,25),i=1,ie)


  END DO ! iteration_loop



!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      tauwatu(i,j)=cw*rhoref_ice*speedu(i,j)*(sicuo(i,j)-uko(i,j,1))
      tauwatv(i,j)=cw*rhoref_ice*speedv(i,j)*(sicve(i,j)-vke(i,j,1))
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL


  IF( have_g_js .and.  .not. lbounds_exch_tp ) THEN
     tauwatu(:,1) = 0._wp
     tauwatv(:,1) = 0._wp
  ENDIF
  IF( have_g_je.and.  .not. lbounds_exch_tp ) THEN
     tauwatu(:,je) = 0._wp
     tauwatv(:,je) = 0._wp
  ENDIF

!#ifdef bounds_exch_save
  CALL bounds_exch(1,'u',tauwatu,'mo_ocice 133')
  CALL bounds_exch(1,'v',tauwatv,'mo_ocice 134')
!#endif

END SUBROUTINE ice_dynamics



!>
!! Compute advection of sea ice related parameters
!!
!! Relates to global variables.
!!
!! @par Side-effects
!! - mo_commo1::zo
!!
SUBROUTINE ice_advection

  REAL(wp) uh(ie,je),vh(ie,je),wh(ie,je),zschalt

  INTEGER i,j

  REAL(wp) uwin,uwou,uein,ueou,vsin,vsou,vnin,vnou

  zschalt = 1._wp

!$OMP PARALLEL private(i,j,uwin,uwou,uein,ueou,vsin,vsou,vnin,vnou)

!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      uwin=0.5_wp*(sicuo(i-1,j)+ABS(sicuo(i-1,j)))*dt*dlyu(i-1,j)
      uwou=0.5_wp*(ABS(sicuo(i-1,j))-sicuo(i-1,j))*dt*dlyu(i-1,j)
      uein=0.5_wp*(ABS(sicuo(i,j))-sicuo(i,j))*dt*dlyu(i,j)
      ueou=0.5_wp*(ABS(sicuo(i,j))+sicuo(i,j))*dt*dlyu(i,j)
      vsin=0.5_wp*(ABS(sicve(i,j))+sicve(i,j))*dt*dlxv(i,j)
      vsou=0.5_wp*(ABS(sicve(i,j))-sicve(i,j))*dt*dlxv(i,j)
      vnin=0.5_wp*(ABS(sicve(i,j-1))-sicve(i,j-1))*dt*dlxv(i,j-1)
      vnou=0.5_wp*(ABS(sicve(i,j-1))+sicve(i,j-1))*dt*dlxv(i,j-1)

      uh(i,j)=(sictho(i,j)*(area(i,j)-uwou-ueou-vsou-vnou)               &
             +sictho(i-1,j)*uwin+sictho(i+1,j)*uein+                  &
              sictho(i,j+1)*vsin+sictho(i,j-1)*vnin)*areain(i,j)

      vh(i,j)=(sicomo(i,j)*(area(i,j)-uwou-ueou-vsou-vnou)               &
             +sicomo(i-1,j)*uwin+sicomo(i+1,j)*uein+                  &
             sicomo(i,j+1)*vsin+sicomo(i,j-1)*vnin)*areain(i,j)

      wh(i,j)=(sicsno(i,j)*(area(i,j)-uwou-ueou-vsou-vnou)               &
             +sicsno(i-1,j)*uwin+sicsno(i+1,j)*uein+                  &
             sicsno(i,j+1)*vsin+sicsno(i,j-1)*vnin)*areain(i,j)

      ! Include effect on sealevel zo (by Uwe M.)

      zo(i,j)=zo(i,j)+(uh(i,j)-sictho(i,j))*rhoicwa*zschalt &
                     +(wh(i,j)-sicsno(i,j))*rhosnwa*zschalt

    ENDDO
  ENDDO
!$OMP END DO

! Update sea ice variables only after advection is completed!

!$OMP DO
  DO j=2,je1
    DO i=2,ie1

      sictho(i,j)=uh(i,j)
      sicomo(i,j)=vh(i,j)
      sicsno(i,j)=wh(i,j)

    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL


END SUBROUTINE ice_advection

SUBROUTINE ice_thermodynamics(siotho, sioomo, siosno)

  REAL(wp), INTENT (INOUT) :: siotho(ie,je),sioomo(ie,je),siosno(ie,je)

#ifdef PBGC
  INTEGER :: l
#endif

  CALL bounds_exch(1,'p',zo,'mo_ocice 135/155')
  CALL bounds_exch(1,'p',sictho,'mo_ocice 136/143')
  CALL bounds_exch(1,'p',sicomo,'mo_ocice 137/144')
  CALL bounds_exch(1,'p',sicsno,'mo_ocice 138/145')


#ifdef PBGC
  CALL swr_frac_from_chl ! calculate available sw radiation fraction ; used in hamocc

  IF (LFB_BGC_OCE) THEN
    CALL dynamic_swr_absorption  ! calculate dynamic absortion factor from swr_frac
  ENDIF
#endif

  IF (lsaoclose) CALL dilcor_gtrf

#ifdef PBGC
  CALL dilcor_gtrf2
#endif

#ifdef __coupled

  heatabs(:,:)=0.  ! can be done in growth

  CALL growth(siotho,sioomo,siosno,sictho,sicomo,sicsno,            &
       aoflfrwo,aoflfrio,aoflnhwo,aoflchio,aoflrhio,       &
       aoflshwo,                                                    &
       sao(:,:,1),tho(:,:,1),dzw(1),zo,dt,fwo,weto(:,:, 1),                &
       qswo,qlwo,qseo,qlao,preco,prech,heatabs,swsum)


#else/*def __coupled*/
!#ifdef  bounds_exch_save
  CALL bounds_exch(1,'p',fprec,'mo_ocice 146')
  CALL bounds_exch(1,'p',tafo,'mo_ocice 147')
  CALL bounds_exch(1,'p',ftdew,'mo_ocice 148')
  CALL bounds_exch(1,'p',fclou,'mo_ocice 149')
  CALL bounds_exch(1,'p',fslp,'mo_ocice 150')
  CALL bounds_exch(1,'p',fu10,'mo_ocice 151')
  CALL bounds_exch(1,'p',fswr,'mo_ocice 152')
  CALL bounds_exch(1,'p',sao,'mo_ocice 153')
  CALL bounds_exch(1,'p',tho,'mo_ocice 154')
  CALL bounds_exch(1,'p',fwo,'mo_ocice 156')
  CALL bounds_exch(1,'p',weto,'mo_ocice 158')
  CALL bounds_exch(1,'p',qswo,'mo_ocice 159')
  CALL bounds_exch(1,'p',qlwo,'mo_ocice 160')
  CALL bounds_exch(1,'p',qseo,'mo_ocice 161')
  CALL bounds_exch(1,'p',qlao,'mo_ocice 162')
  CALL bounds_exch(1,'p',heatabs,'mo_ocice 165')
  CALL bounds_exch(1,'p',preco,'mo_ocice 163')
  CALL bounds_exch(1,'p',prech,'mo_ocice 164')
!#endif


  heatabs(:,:) = 0._wp  ! can be done in growth

  CALL growth(siotho,sioomo,siosno,sictho,sicomo,sicsno,        &
       fprec, & ! used to be extended by: + giriv(i,j)
       tafo,ftdew,fclou,fslp,fu10,fswr,                         &
       sao(:,:,1),tho(:,:,1),dzw(1),zo,dt,fwo,weto(:,:,1),         &
       qswo,qlwo,qseo,qlao,preco,prech,heatabs,swsum)


#endif/*else def __coupled*/
!#ifdef bounds_exch_save
  CALL bounds_exch(1,'p',siotho,'mo_ocice 166')
  CALL bounds_exch(1,'p',sioomo,'mo_ocice 167')
  CALL bounds_exch(1,'p',siosno,'mo_ocice 168')
  CALL bounds_exch(1,'p',sictho,'mo_ocice 169')
  CALL bounds_exch(1,'p',sicomo,'mo_ocice 170')
  CALL bounds_exch(1,'p',sicsno,'mo_ocice 171')
  CALL bounds_exch(1,'p',fprec,'mo_ocice 172')
  CALL bounds_exch(1,'p',tafo,'mo_ocice 173')
  CALL bounds_exch(1,'p',ftdew,'mo_ocice 174')
  CALL bounds_exch(1,'p',fclou,'mo_ocice 175')
  CALL bounds_exch(1,'p',fslp,'mo_ocice 176')
  CALL bounds_exch(1,'p',fu10,'mo_ocice 177')
  CALL bounds_exch(1,'p',fswr,'mo_ocice 178')
  CALL bounds_exch(1,'p',sao,'mo_ocice 179')
  CALL bounds_exch(1,'p',tho,'mo_ocice 180')
  CALL bounds_exch(1,'p',zo,'mo_ocice 181')
  CALL bounds_exch(1,'p',fwo,'mo_ocice 182')
  CALL bounds_exch(1,'p',weto,'mo_ocice 184')
  CALL bounds_exch(1,'p',qswo,'mo_ocice 185')
  CALL bounds_exch(1,'p',qlwo,'mo_ocice 186')
  CALL bounds_exch(1,'p',qseo,'mo_ocice 187')
  CALL bounds_exch(1,'p',qlao,'mo_ocice 188')
  CALL bounds_exch(1,'p',preco,'mo_ocice 189')
  CALL bounds_exch(1,'p',prech,'mo_ocice 190')
  CALL bounds_exch(1,'p',heatabs,'mo_ocice 191')
!#endif

#ifdef PBGC
  DO l=1,nocetra
    CALL dilcor_ptrf2(ocetra(1,1,1,l))
  END DO
#endif

  IF (lsaoclose) CALL dilcor_ptrf

  CALL subsurface_swr_absorption ! this can go into growths


END SUBROUTINE ice_thermodynamics

#ifdef MESSY
SUBROUTINE ice_thermodynamics_coupled(siotho, sioomo, siosno,      &
           AOFLFRIO,AOFLFRWO,AOFLRHIO,AOFLCHIO,AOFLNHWO,AOFLSHWO,AOFLWSVO)

  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLFRIO!solid freshwater flux (over ice only) 
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLFRWO!liquid freshwater flux (water and ice)
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLRHIO!solid freshwater flux (over ice only)
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLCHIO!conductive heat flux through ice 
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLNHWO!net heat flux over water
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLSHWO!downwelling solar radiation
  REAL(dp),DIMENSION(:,:),INTENT(IN):: AOFLWSVO!wind stress velocity 

  REAL(wp), INTENT (INOUT) :: siotho(ie,je),sioomo(ie,je),siosno(ie,je)

#ifdef PBGC
  INTEGER :: l
#endif

  CALL bounds_exch(1,'p',zo,'mo_ocice 135/155')
  CALL bounds_exch(1,'p',sictho,'mo_ocice 136/143')
  CALL bounds_exch(1,'p',sicomo,'mo_ocice 137/144')
  CALL bounds_exch(1,'p',sicsno,'mo_ocice 138/145')


#ifdef PBGC
  CALL swr_frac_from_chl ! calculate available sw radiation fraction ; used in hamocc

  IF (LFB_BGC_OCE) THEN
    CALL dynamic_swr_absorption  ! calculate dynamic absortion factor from swr_frac
  ENDIF
#endif

  IF (lsaoclose) CALL dilcor_gtrf

#ifdef PBGC
  CALL dilcor_gtrf2
#endif

  heatabs(:,:)=0.  ! can be done in growth

  CALL growth_coupled(siotho,sioomo,siosno,sictho,sicomo,sicsno,            &
       aoflfrwo,aoflfrio,aoflnhwo,aoflchio,aoflrhio,       &
       aoflshwo,                                                    &
       sao(:,:,1),tho(:,:,1),dzw(1),zo,dt,fwo,weto(:,:, 1),                &
       qswo,qlwo,qseo,qlao,preco,prech,heatabs,swsum)


  !#ifdef bounds_exch_save
  CALL bounds_exch(1,'p',siotho,'mo_ocice 166')
  CALL bounds_exch(1,'p',sioomo,'mo_ocice 167')
  CALL bounds_exch(1,'p',siosno,'mo_ocice 168')
  CALL bounds_exch(1,'p',sictho,'mo_ocice 169')
  CALL bounds_exch(1,'p',sicomo,'mo_ocice 170')
  CALL bounds_exch(1,'p',sicsno,'mo_ocice 171')
  CALL bounds_exch(1,'p',fprec,'mo_ocice 172')
  CALL bounds_exch(1,'p',tafo,'mo_ocice 173')
  CALL bounds_exch(1,'p',ftdew,'mo_ocice 174')
  CALL bounds_exch(1,'p',fclou,'mo_ocice 175')
  CALL bounds_exch(1,'p',fslp,'mo_ocice 176')
  CALL bounds_exch(1,'p',fu10,'mo_ocice 177')
  CALL bounds_exch(1,'p',fswr,'mo_ocice 178')
  CALL bounds_exch(1,'p',sao,'mo_ocice 179')
  CALL bounds_exch(1,'p',tho,'mo_ocice 180')
  CALL bounds_exch(1,'p',zo,'mo_ocice 181')
  CALL bounds_exch(1,'p',fwo,'mo_ocice 182')
  CALL bounds_exch(1,'p',weto,'mo_ocice 184')
  CALL bounds_exch(1,'p',qswo,'mo_ocice 185')
  CALL bounds_exch(1,'p',qlwo,'mo_ocice 186')
  CALL bounds_exch(1,'p',qseo,'mo_ocice 187')
  CALL bounds_exch(1,'p',qlao,'mo_ocice 188')
  CALL bounds_exch(1,'p',preco,'mo_ocice 189')
  CALL bounds_exch(1,'p',prech,'mo_ocice 190')
  CALL bounds_exch(1,'p',heatabs,'mo_ocice 191')
!#endif

#ifdef PBGC
  DO l=1,nocetra
    CALL dilcor_ptrf2(ocetra(1,1,1,l))
  END DO
#endif

  IF (lsaoclose) CALL dilcor_ptrf

  CALL subsurface_swr_absorption ! this can go into growths

END SUBROUTINE ice_thermodynamics_coupled
#endif



END MODULE mo_ocice

