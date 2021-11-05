!> Calulation of Advection
!> @author Ernst Maier-Reimer, Uwe Mikolajewicz, Stephan Lorenz, ...
!> @date Last modified on 11.12.2008 by Nils Fischer and Mario Krapp
!> @todo Documentation of BBL part of this routine (Johann?)
!> @todo "old" TO DO (still needed???): BBL and ocadpo are not default.
!>  Is it necessary to compute the new vertical velocities in ocadpo (20 times with BGC)?
!!
!> This module calculates the advection of tracer fields. Corresponding
!! to the namelist parameter 'iocad' (1,2,3,8) different schemes are used:
!! iocad=1   upwind advection scheme (TVD, 1st order)
!! iocad=2   central difference scheme
!! iocad=3   weighted central difference/upwind scheme (2nd order AND TVD)
!! iocad=8   as iocad=3 but vertical tracer transport with one third of time step
!!
!!First, new vertical velocities are calculated, including the BBL transport velocities.
!!Second, advection of scalar traces is computed with a second order total variation
!!diminishing (TVD) scheme (Sweby,1984).
!!The total variation of a solution is defined as:
!!@f$ TV_{(u^{n+1})} = \sum{\mid u^{n+1}_{k+1}-u^{n+1}_{k}\mid} @f$
!!A difference scheme is defined as being total variation
!!diminishing (TVD) if:
!!@f$ TV_{(u^{n+1})} \le TV_{(u^{n})} @f$
!!Momentum advection of tracers is by a mixed scheme that
!!employs a weighted average of both central-difference and upstream methods.
!!The weights are chosen in a two step process. First, according to the
!!ratio of the first minus the second spatial derivative
!!over the first spatial derivative of the advected quantity T,
!!@f$ r=max\left\{0,{|T^{\prime}|-|T^{\prime\prime}| \over |T^{\prime}|}\right\} @f$
!!with @f$ T^{\prime}=T_{k-1}-T_{k+1}@f$ and @f$T^{\prime\prime}=T_{k-1}+T_{k+1}-2\cdot T_{k} @f$.
!!If the second derivative is small usage of central-differencing is save and therefor favorably.
!!In contrast, if the first derivative is small and the second derivative is large there
!!is an extrema in the middle and an upstream scheme is preferred.
!!
!!In a second step the ratio is weighted with strength of the flow
!!(the time it takes to fully ventilate the grid-box).
!!Water transport in and out of a grid-box is given as:
!!@f{eqnarray*} U_{in}&=&0.5\cdot\delta t \cdot\delta x \cdot\delta y \cdot(w_k +|w_{k-1}|) \nonumber \\\
!! U_{out}&=&0.5\cdot\delta t \cdot\delta x \cdot\delta y \cdot(|w_{k+1}|-w_{k+1}) @f}
!!The total weight @f$R@f$ is defined as:
!!@f$ R=min\left\{1,{\delta x \cdot\delta y \cdot\delta z \over U_{in}+U_{out}}\cdot r\right\} @f$
!!If the flow is weak and @f$ r @f$ is small, the magnitude of this ratio
!!is less than 1 and the weights favor usage of central-differencing.
!!With a stronger flow or a larger @f$ r @f$ the upstream scheme is preferred.
!!The idea is to incorporate the benefit of positive-definiteness of the upstream scheme
!!(and thus limit numerically spurious tracer sources and sinks),
!!while avoiding large implicit numerical diffusion in regions where strong gradients
!!exist in the tracer field.
!!
!!Advection is computed as follows. Tracer transport in and out of a grid-box is defined as:
!!@f{eqnarray*}
!!  T_{in} & = & U_{in}\cdot(1-R)\cdot T_k + R\cdot 0.5 \cdot (T_k+T_{k-1}) \nonumber \\\
!! T_{out}&=& U_{out}\cdot(1-R)\cdot T_k + R\cdot 0.5 \cdot (T_k+T_{k+1})
!! @f}
!!
!!The new tracer concentration @f$ T_k^{n+1} @f$ is given by the old concentration @f$ T_k^n @f$
!!plus tracer in-, and out-fluxes,
!!normalized by the "new" volume of the grid-box (volume plus in-, and out-fluxes of water).
!!@f{eqnarray*}
!!  T_k^{n+1}&=&  \left(T_k^n \cdot\delta x \cdot\delta y \cdot\delta z + T_{in \, k+1} - T_{in \, k} \
!! - T_{out \, k} + T_{out \, k-1}\right)\cdot{1\over V_{new}} \nonumber \\\
!! V_{new}&=& \delta x \cdot\delta y \cdot\delta z + U_{in \, k+1} - U_{in \, k} - U_{out \, k} + U_{out \, k-1} @f}
!!
!! References:
!! High Resolution Schemes Using Flux Limiters for Hyperbolic Conservation Laws.
!! P.K. Sweby, SIAM Journal on Numerical Analysis, Vol. 21, No. 5. 1984.
MODULE MO_ADPO
  USE mo_kind, ONLY: wp
  USE MO_PARAM1
  USE mo_boundsexch, ONLY : bounds_exch
  IMPLICIT NONE
!  REAL(wp) ,ALLOCATABLE :: TRPHELP(:,:),TRMHELP(:,:)  ! help fields needed for tripolar setups

  REAL(wp) ,ALLOCATABLE :: woh1(:,:,:) !< array for vertical velocity(BBL corrected)
  REAL(wp), ALLOCATABLE :: svo(:,:,:)  !< array for the volume of the gridbox
  REAL(wp), ALLOCATABLE :: dlxypa(:,:) !< array for the value for the x-y plane

CONTAINS
!>
!! Allocating memory for the different fields
!!
SUBROUTINE ALLOC_MEM_ADPO

!   if ( lbounds_exch_tp ) then
!      ALLOCATE(TRPHELP(IE,KE),TRMHELP(IE,KE))
!   endif
   allocate(woh1(ie,je,kep))
   allocate(dlxypa(ie,je),svo(ie,je,ke))

END SUBROUTINE ALLOC_MEM_ADPO

!> Computes bottom boundary layer velocities (ibbl_transport=1)
!! for each timestep to be used in sbrt ocadpo_trf
!!
!! By Stephan Lorenz 08/2007
!!

SUBROUTINE ocadpo_base

! mak: Which variables are used from other modules? Suggestion: -> USE <module>, ONLY: <variables>
  USE mo_param1   ! ie,ie1,je,je1,ke
  USE mo_commo1   ! dlxp,dlyp,wo,weto,ddpo,iocad,zo,sictho,sicsno,dt
  USE mo_planetary_constants, ONLY: rhoicwa,rhosnwa
  USE mo_commobbl ! ubbl,vbbl,kupubbl,kupvbbl,kdwubbl,kdwvbbl

  INTEGER :: i, j, k
  REAL(wp) :: suminf
!$OMP PARALLEL PRIVATE(i,j,k,suminf)


!$OMP WORKSHARE
! Help field initialised with wo (upward_sea_water_velocity)
   woh1(:,:,:) = wo
!$OMP END WORKSHARE

! Initializing the field which holds the volume of every grid box (svo).
! Uses: land-sea mask for all levels (weto), the area and the thickness of the grid box.

!$OMP DO
  DO k=1,ke
    DO j=1,je
      DO i=1,ie
        svo(i,j,k)=weto(i,j,k)* area(i,j) *ddpo(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!     COMPUTE NEW VERTICAL VELOCITIES INCLUDING BBL TRANSPORTS

  if (ibbl_transport==1) then

!$OMP DO
! BBL shift is calculated.

    DO K=KE,1,-1
      DO J=2,JE1
        DO I=2,IE1
          suminf = 0._wp

          IF (KDWUBBL(I-1,J).EQ.K) SUMINF=SUMINF+UBBL(I-1,J)
          IF (KUPUBBL(I-1,J).EQ.K) SUMINF=SUMINF-UBBL(I-1,J)
          IF (KDWUBBL(I,J).EQ.K) SUMINF=SUMINF-UBBL(I,J)
          IF (KUPUBBL(I,J).EQ.K) SUMINF=SUMINF+UBBL(I,J)

          IF (KDWVBBL(I,J).EQ.K) SUMINF=SUMINF+VBBL(I,J)
          IF (KUPVBBL(I,J).EQ.K) SUMINF=SUMINF-VBBL(I,J)
          IF (KDWVBBL(I,J-1).EQ.K) SUMINF=SUMINF-VBBL(I,J-1)
          IF (KUPVBBL(I,J-1).EQ.K) SUMINF=SUMINF+VBBL(I,J-1)

! new vertical velocity from bbl transports is computed. wo/woh1 at
! bottom are zero, woh1 initiali values are overwritten by bbl
! velocity.

          woh1(I,J,K)=woh1(I,J,K+1)+SUMINF*areain(I,J)

        ENDDO
      ENDDO
    ENDDO

!$OMP DO
    DO K=1,KE
      DO J=2,JE1
        DO I=2,IE1
! new vertical velocity =  bbl vertical velocity + non bbl vertical velocity
          woh1(I,J,K)=woh1(I,J,K)+wo(I,J,K)
        ENDDO
      ENDDO
    ENDDO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch(1,'p',woh1,'ocadpo_base 1')
!$OMP END SINGLE
!#endif

  endif ! ibbl_transport=1

!     Compute old volume of surface layer according to old sealevel
!       before update of ZO in update_zo

!$OMP DO
  DO J=1,JE
    DO I=1,IE
      ! svo(surface volume new) = svo(old) + dlxypa (gridbox area) * &
      ! weto (land/sea mask) * ZO (sea surface height old)- woh1(vertical velocity old + bbl velocity) * &
      ! DT - SICTHO (sea ics thickness) * RHOICWA(ICe/WAter volume ratio)
      !    - SICSNO (SNOw thickness) * RHOSNWA (SNOw/WAter volume ratio)
       ! The surface volume is calculated by adding/substracting  the contribution of
       ! sea level (zo) , sea ice (sictho), snow on sea ice (sicsno).

      svo(I,J,1)=svo(I,J,1)+ area(i,j) *WETO(I,J,1)*              &
           (ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)
    ENDDO
  ENDDO

!$OMP END PARALLEL

END SUBROUTINE ocadpo_base
!> ROUTINE OCADPO
!!
!! COMPUTES ADVECTION OF ONE TRACERFIELD
!!
!! BY ERNST MAIER-REIMER 1/1999
!! MODIFIED 1/2000 UWE MIKOLAJEWICZ
!! MAKE MASS CONSERVING WITH FREE SURFACE LAYER
!! MODIFIED 3/2000 UWE MIKOLAJEWICZ
!! INCLUDE BOTTOM BOUNDARY LAYER TRANSPORTS IN
!! ADVECTION SCHEME (IBBL_TRANSPORT=1)
!! NEEDS TO BE RUN TOGETHER WITH SR SLOPETRANS
!!
!! INPUT/OUTPUT
!! TRF(IE,JE,KE)     TRACER FIELD
!!
!! USES VELOCITIES (UKO,VKE,WO)
!!
!! Changes by S. Lorenz, J.-O. Biesmann, 2007-08 optimised
!! vertical velocities updated in ocadpo_base
SUBROUTINE ocadpo_trf(trf)

  USE mo_param1   ! je1,ke,ie1
  USE mo_commo1   ! iocad,dt,ddue,dduo,dlxv,vke,uko
  USE mo_commobbl ! kdwubbl,kupubbl,kdwvbbl,kupvbbl
!  USE mo_parallel, ONLY : p_joff

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif

  IMPLICIT NONE

  REAL(wp) :: trf(ie,je,ke)    !< array for a tracer field (e.g. salinity, temperature)
  REAL(wp) :: s2o(ie,je,ke)    !< array for new gridbox volume
  REAL(wp) :: trp(ie,je,kep)   !< array for tracer transport into a grid box
  REAL(wp) :: trm(ie,je,0:ke)  !< array for tracer transport out of a grid box
  REAL(wp) :: wtp(ie,je,kep)   !< array for water transport into a grid box
  REAL(wp) :: wtm(ie,je,0:ke)  !< array for water transport out of a grid box
  REAL(wp) :: b1u, b1um, b1v, b1vm !< Something for BBL-transport (TO DO! Johann)
  REAL(wp) :: up, um           !< water transport into/out of a grid box
  REAL(wp) :: abl, zwabl       !< first and second spatial derivatives of tracer
  REAL(wp) :: scha, sch, schi  !< weighting factors r, R (see documentation)
  REAL(wp) :: rn               !< updated grid box volume

  INTEGER :: i, j, k

#ifdef _PROFILE
  CALL trace_start ('ocadpo_trf 1', 6)
#endif

#ifndef TRACER_OMP
!$OMP PARALLEL PRIVATE(i,j,k,abl,zwabl,up,um,scha,sch,schi,rn,   &
!$OMP                  b1u,b1um,b1v,b1vm)
#endif


!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'p',trf,'ocadpo_trf 1')
!$OMP END SINGLE
!!#endif

! initialise local variables
#ifndef TRACER_OMP
!$OMP WORKSHARE
#endif
  s2o(:,:,:)=svo(:,:,:)
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif


!!
!! Vertical, Zonal, and Meridional Transports
!!

!! VERTICAL TRANSPORTS

  CALL ocadpo_trf_z(trf, s2o, wtm, wtp, trm, trp)


!! TRANSPORTS IN X-DIRECTION
!! comments same as in vertical transport (s.a.) except for BBL-stuff (to do!)


#ifndef TRACER_OMP
!$OMP DO
#endif
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        b1u  = 0._wp
        b1um = 0._wp
        IF ( ibbl_transport == 1 ) THEN
          IF (KDWUBBL(I,J).EQ.K) b1u=UBBL(I,J)
          IF (KUPUBBL(I,J).EQ.K) b1u=-UBBL(I,J)
          IF (KDWUBBL(I-1,J).EQ.K) b1um=UBBL(I-1,J)
          IF (KUPUBBL(I-1,J).EQ.K) b1um=-UBBL(I-1,J)
        ENDIF
        ABL=ABS(TRF(I+1,J,K)-TRF(I-1,J,K))
        zwabl = ABS(trf(i+1, j, k) + trf(i-1, j, k) - 2._wp * trf(i, j, k))

        up = 0.5_wp * dt * dduo(i, j, k) * dlyu(i, j) &
             * (uko(i, j, k) + ABS(uko(i, j, k))) &
             + 0.5_wp * dt * (b1u + ABS(b1u))

        um = 0.5_wp * dt * dduo(i-1, j, k) &
             * (ABS(uko(i-1, j, k)) - uko(i-1, j, k)) * dlyu(i-1, j) &
             + 0.5_wp * dt * (ABS(b1um) - b1um)

        SELECT CASE (iocad)
        CASE (1)     ! pure upwind scheme
          sch = 0._wp
        CASE (2)     ! central-difference scheme
          sch = 1._wp
        CASE DEFAULT
          scha = MAX(0._wp, (abl - zwabl) /(abl + 1.E-20_wp)) * weto(i, j, k) &
               * weto(i+1, j, k) * weto(i-1, j, k)
#ifdef SMOADH
          sch = MIN(1._wp, scha)
#else
          sch = MIN(1._wp, scha*s2o(i, j, k) / (up + um + 1.E-20_wp))
#endif
        END SELECT

        schi = 1._wp - sch

        trp(i, j, k) = up * (schi * trf(i, j, k) + sch * 0.5_wp &
             * (trf(i, j, k) + trf(i+1, j, k)))
        trm(i, j, k) = um * (schi * trf(i, j, k) + sch * 0.5_wp &
             * (trf(i, j, k) + trf(i-1, j, k)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'u+',TRP,'ocadpo_trf 5')
  CALL bounds_exch(1,'u+',TRM,'ocadpo_trf 6')
  CALL bounds_exch(1,'u+',WTP,'ocadpo_trf 7')
  CALL bounds_exch(1,'u+',WTM,'ocadpo_trf 8')
!$OMP END SINGLE
!#endif

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        IF (weto(i, j, k) .GT. 0.5_wp) THEN
          RN=s2o(I,J,K)+WTP(I-1,J,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I+1,J,K)
          TRF(I,J,K)=(TRF(I,J,K)*s2o(I,J,K)+TRP(I-1,J,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I+1,J,K))/RN
          s2o(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'p',TRF,'ocadpo_trf 9')
  CALL bounds_exch(1,'p',s2o,'ocadpo_trf 10')
!$OMP END SINGLE
!#endif

!! VERTICAL TRANSPORTS with one third of time step
  IF (iocad==8) CALL ocadpo_trf_z(trf, s2o, wtm, wtp, trm, trp)


!! TRANSPORTS IN Y-DIRECTION
!! comments same as in vertical transport (s.a.) except for BBL-stuff (to do!)

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO K=1,KE
    DO J=2,JE1

!     hh: is this really needed ?
!      IF ( lbounds_exch_tp .AND. p_joff.EQ.0 ) THEN
!        DO i=1,ie
!          trf(i,1,k)=2.*trf(i,2,k)-trf(i,3,k)
!        ENDDO
!      ENDIF

      DO I=2,IE1
        b1v = 0._wp
        b1vm = 0._wp
        IF ( ibbl_transport == 1 ) THEN
          IF (K.EQ.KDWVBBL(I,J)) b1v=VBBL(I,J)
          IF (K.EQ.KUPVBBL(I,J)) b1v=-VBBL(I,J)
          IF (K.EQ.KDWVBBL(I,J-1)) b1vm=VBBL(I,J-1)
          IF (K.EQ.KUPVBBL(I,J-1)) b1vm=-VBBL(I,J-1)
        ENDIF

        ABL=ABS(TRF(I,J-1,K)-TRF(I,J+1,K))
        zwabl = ABS(trf(i, j-1, k) + trf(i, j+1, k) - 2._wp * trf(i, j, k))
        up = 0.5_wp * dt * ddue(i, j-1, k) * dlxv(i, j-1) &
             &  * (vke(i, j-1, k) + ABS(vke(i, j-1, k))) &
             + 0.5_wp * dt * (b1vm + ABS(b1vm))
        um = 0.5_wp * dt*ddue(i, j, k) * dlxv(i, j) * (ABS(vke(i, j, k))   &
             - vke(i, j, k)) + 0.5_wp * dt * (ABS(b1v) - b1v)

        SELECT CASE (iocad)
        CASE (1)     ! pure upwind scheme
          sch = 0._wp
        CASE (2)     ! central-difference scheme
          sch = 1._wp
        CASE DEFAULT
          scha = MAX(0._wp, (abl - zwabl) / (abl + 1.E-20_wp)) * weto(i, j, k) &
               * weto(I,J-1,K)*WETO(I,J+1,K)
#ifdef SMOADH
          sch = MIN(1._wp, scha)
#else
          sch = MIN(1._wp, scha * s2o(i, j, k)/(up + um + 1.E-20_wp))
#endif
        END SELECT

        schi = 1._wp - sch
        trp(i, j, k) = up * (schi * trf(i, j, k) + sch * 0.5_wp &
             * (trf(i, j, k) + trf(i, j-1, k)))
        trm(i, j, k) = um &
             * (schi * trf(i, j, k) &
             &  + sch * 0.5_wp * (trf(i, j, k) + trf(i, j+1, k)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'vv',TRP,'ocadpo_trf 11')
  CALL bounds_exch(1,'vv',TRM,'ocadpo_trf 12')
  CALL bounds_exch(1,'vv',WTP,'ocadpo_trf 13')
  CALL bounds_exch(1,'vv',WTM,'ocadpo_trf 14')
!$OMP END SINGLE
!#endif


#ifndef TRACER_OMP
!$OMP DO
#endif
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        IF (lweto(i,j,k)) THEN

          RN=s2o(I,J,K)+WTP(I,J+1,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J-1,K)
          TRF(I,J,K)=(TRF(I,J,K)*s2o(I,J,K)+TRP(I,J+1,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J-1,K))/RN
          s2o(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch(1,'p',TRF,'ocadpo_trf 15')
!$OMP END SINGLE
!#endif


!! VERTICAL TRANSPORTS with one third of time step
  IF (iocad==8) CALL ocadpo_trf_z(trf, s2o, wtm, wtp, trm, trp)

#ifndef TRACER_OMP
!$OMP END PARALLEL
#endif

#ifdef _PROFILE
  CALL trace_stop ('ocadpo_trf 1', 6)
#endif
END SUBROUTINE ocadpo_trf


SUBROUTINE ocadpo_trf_z(trf, s2o, wtm, wtp, trm, trp)

  USE mo_param1,   ONLY: ie, je, ke, ke1, kep
  USE mo_commo1,   ONLY: dt, area, weto, lweto, iocad

  IMPLICIT NONE

  REAL(wp), INTENT(INOUT) :: trf(ie,je,ke)    !< array for a tracer field (e.g. salinity, temperature)
  REAL(wp), INTENT(INOUT) :: s2o(ie,je,ke)    !< array for new gridbox volume

  REAL(wp), INTENT(INOUT) :: wtm(ie,je,0:ke)  !< array for "out" water transport
  REAL(wp), INTENT(INOUT) :: wtp(ie,je,kep)   !< array for "in"  water transport
  REAL(wp), INTENT(INOUT) :: trm(ie,je,0:ke)  !< array for "out" tracer transport
  REAL(wp), INTENT(INOUT) :: trp(ie,je,kep)   !< array for "in"  tracer transport
  REAL(wp)                :: um, up           !  out/in water transport
  REAL(wp)                :: abl, zwabl       !  first and second spatial derivatives of tracer
  REAL(wp)                :: scha, sch, schi  !  weighting factors r, R (s. docu)
  REAL(wp)                :: rn               !  updated grid box volume
  REAL(wp)                :: dtfrac           !  time step fraction for vertical transports

  INTEGER             :: i, j, k, klo, kup


!! VERTICAL TRANSPORTS with fractional time step defined by dtfrac

#ifndef TRACER_OMP
!$OMP PARALLEL PRIVATE(i,j,k,klo,kup,abl,zwabl,up,um,scha,sch,schi,rn)
#endif

  IF (iocad==8) THEN
    dtfrac = 1._wp/3._wp! vertical tracer transport with one third of time step
  ELSE
    dtfrac = 1._wp
  ENDIF

#ifndef TRACER_OMP
!$OMP DO
#endif
  ! Special treatment of top (0 or 1) and bottom (ke or kep) layers
  DO J=1,JE
    DO I=1,IE
      ! weighted water transport downward (i.e. "out") UM dependent on direction water flux woh1 from lower box (k=2)
      um = 0.5_wp * dt * area(i, j) &
           * (ABS(woh1(i, j, 2)) - woh1(i, j, 2)) * dtfrac
      ! Tracer transport up/downward (i.e. "in"/"out") TRP/TRM
      trp(i, j, 1) = 0._wp
      trm(i, j, 0) = 0._wp
      TRM(I,J,1)=UM*TRF(I,J,1)
      ! assign water transport
      wtp(i, j, 1) = 0._wp
      wtm(i, j, 0) = 0._wp
      WTM(I,J,1)=UM

      ! weighted water transport up (i.e. "in") UP dependent on direction water flux woh1 from lowermst box to box above
      up = 0.5_wp * dt * area(i, j) &
           * (woh1(i, j, ke) + ABS(woh1(i, j, ke))) * dtfrac
      um = 0.5_wp * dt * area(i, j) &
           * (ABS(woh1(i, j, ke+1)) - woh1(i, j, ke+1)) * dtfrac
      ! Tracer transport up/downward (i.e. "in"/"out") TRP/TRM
      TRP(I,J,KE)=UP*TRF(I,J,KE)
      TRM(I,J,KE)=UM*TRF(I,J,KE)
      ! assign water transport
      WTP(I,J,KE)=UP
      WTM(I,J,KE)=UM
      ! no transport in the bottom
      wtp(i, j, kep) = 0._wp
      trp(i, j, kep) = 0._wp
    END DO
  ENDDO

  DO K=2,KE1

    KLO=MIN(K+1,KE)   ! grid-box below
    KUP=MAX(K-1,1)    ! grid-box above
    DO J=1,JE
      DO I=1,IE

        ! weighted water transport up/downwards UP/UM dependent on direction water
        ! flux woh1 from lowermst box to box above
        up = 0.5_wp * dt * area(i, j) &
             * (woh1(i, j, k) + ABS(woh1(i, j, k))) * dtfrac
        um = 0.5_wp * dt * area(i, j) &
             * (ABS(woh1(i, j, k+1)) - woh1(i, j, k+1)) * dtfrac

        SELECT CASE (iocad)
        CASE (1)     ! pure upwind scheme
          sch = 0._wp
        CASE (2)     ! central-difference scheme
          sch = 1._wp
        CASE DEFAULT
          ! Calculating first and second derivatives (ABL, ZWABL), weighting SCHA
          ! (corresponding to "r" in documentation)
          ABL=ABS(TRF(I,J,KUP)-TRF(I,J,KLO))
          zwabl = ABS(trf(i, j, kup) + trf(i, j, klo) - 2._wp * trf(i, j, k))
          scha = MAX(0._wp, (abl - zwabl) / (abl + 1.E-20_wp))

#ifdef SMOADV
          ! SCH is corresponding to "R" in documentation
          SCH=MIN(WETO(I,J,KLO),SCHA)
#else
          sch = MIN(weto(i, j, klo), &
               scha * s2o(i, j, k) / (up + um + 1.E-20_wp))
#endif
        END SELECT

        ! SCHI introduced for faster computation
        schi = 1._wp - sch
        ! Tracer transport up/downward (TRP/TRM)
        trp(i, j, k) = up * (schi * trf(i, j, k) &
             &               + sch * 0.5_wp * (trf(i, j, k) + trf(i, j, kup)))
        trm(i, j, k) = um * (schi * trf(i, j, k) &
             &               + sch * 0.5_wp * (trf(i, j, k) + trf(i, j, klo)))
        ! assign water transport
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
  END DO

  DO k=1,ke
    DO J=1,JE
      DO I=1,IE
        IF (lweto(i,j,k)) THEN
          ! calculate new gridbox volume
          RN=s2o(I,J,k)+WTP(I,J,K+1)-WTP(I,J,k)-WTM(I,J,k)+WTM(I,J,k-1)
          ! update/calculate new tracer field
          TRF(I,J,K)=(TRF(I,J,K)*s2o(I,J,k)+TRP(I,J,k+1)-TRP(I,J,k)      &
                      -TRM(I,J,k)+TRM(I,J,K-1))/RN
          ! assign gridbox volume to s2o
          s2o(I,J,K)=RN
        ENDIF
      END DO
    END DO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
!   CALL bounds_exch(1,'p',TRF,'ocadpo_trf 3')   ! not necessary
!   CALL bounds_exch(1,'p',s2o,'ocadpo_trf 4')   ! not necessary
!$OMP END SINGLE
!#endif

#ifndef TRACER_OMP
!$OMP END PARALLEL
#endif

END SUBROUTINE ocadpo_trf_z


!!!! ATTENTION: SBR ocadpo_trf was modified! The changes are not implemented
!!!! in the SBR ocadpo_trf_halo below
!!$SUBROUTINE ocadpo_trf_halo(trf)
!!$
!!$  !     ROUTINE OCADPO
!!$  !
!!$  !     COMPUTES ADVECTION OF ONE TRACERFIELD
!!$  !
!!$  !     BY ERNST MAIER-REIMER 1/1999
!!$  !     MODIFIED 1/2000 UWE MIKOLAJEWICZ
!!$  !     MAKE MASS CONSERVING WITH FREE SURFACE LAYER
!!$  !     MODIFIED 3/2000 UWE MIKOLAJEWICZ
!!$  !     INCLUDE BOTTOM BOUNDARY LAYER TRANSPORTS IN
!!$  !     ADVECTION SCHEME (IOCAD.EQ.4 :: ADPO + SLOPECON_ADPO)
!!$  !     NEEDS TO BE RUN TOGETHER WITH SR SLOPETRANS
!!$  !
!!$  !     INPUT/OUTPUT
!!$  !     TRF(IE,JE,KE)     TRACER FIELD
!!$  !
!!$  !     USES VELOCITIES (UKO,VKE,WO)
!!$  !
!!$  !     Changes by S. Lorenz, J.-O. Biesmann, 2007-08 optimised
!!$  !       vertical velocities updated in ocadpo_base
!!$
!!$  USE mo_param1
!!$  USE mo_parallel
!!$  USE mo_commo1
!!$  USE mo_commoau1
!!$  USE mo_commobbl
!!$#ifdef _PROFILE
!!$  USE mo_profile,      ONLY: trace_start, trace_stop
!!$#endif
!!$
!!$  REAL(wp) :: trf(ie,je,ke)
!!$
!!$  REAL(wp) :: aux(ie,je),aux_h(-1:ie+2,-1:je+2)
!!$
!!$  REAL(wp) :: trf_h(-1:ie+2,-1:je+2,ke)
!!$  REAL(wp) :: woh2_h(-1:ie+2,-1:je+2,kep)
!!$  REAL(wp) :: area_h(-1:ie+2,-1:je+2)
!!$  REAL(wp) :: weto_h(-1:ie+2,-1:je+2,ke)
!!$  REAL(wp) :: svo_h(-1:ie+2,-1:je+2,ke)
!!$  REAL(wp) :: dduo_h(-1:ie+2,-1:je+2,ke)
!!$  REAL(wp) :: uko_h(-1:ie+2,-1:je+2,ke)
!!$  REAL(wp) :: dlyu_h(-1:ie+2,-1:je+2)
!!$  REAL(wp) :: ddue_h(-1:ie+2,-1:je+2,ke)
!!$  REAL(wp) :: vke_h(-1:ie+2,-1:je+2,ke)
!!$  REAL(wp) :: dlxv_h(-1:ie+2,-1:je+2)
!!$  integer :: kdwubbl_h(-1:ie+2,-1:je+2)
!!$  integer :: kupubbl_h(-1:ie+2,-1:je+2)
!!$  REAL(wp) :: ubbl_h(-1:ie+2,-1:je+2)
!!$  integer :: kdwvbbl_h(-1:ie+2,-1:je+2)
!!$  integer :: kupvbbl_h(-1:ie+2,-1:je+2)
!!$  REAL(wp) :: vbbl_h(-1:ie+2,-1:je+2)
!!$
!!$
!!$
!!$  REAL(wp) :: s2o_h(-1:ie+2,-1:je+2,ke)
!!$  REAL(wp) :: trp_h(-1:ie+2,-1:je+2,kep), trm_h(-1:ie+2,-1:je+2,0:ke)
!!$  REAL(wp) :: wtp_H(-1:ie+2,-1:je+2,kep), wtm_h(-1:ie+2,-1:je+2,0:ke)
!!$
!!$  REAL(wp) :: b1u, b1um, b1v, b1vm, zzsurf
!!$
!!$#ifdef _PROFILE
!!$  CALL trace_start ('mo_adpotrf 1', 6)
!!$#endif
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP PARALLEL PRIVATE(i,j,k,klo,kup,suminf,abl,zwabl,up,um,scha,sch,schi,rn,   &
!!$!$OMP                  zzsurf,b1u,b1um,b1v,b1vm)
!!$#endif
!!$
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$ CALL bounds_exch(1,'p',trf,'ocadpo_trf 1')
!!$!$OMP END SINGLE
!!$!!#endif
!!$
!!$ call sethalon_new('p',area,area_h,2)
!!$ call sethalon_new('u+',dlyu,dlyu_h, 2)
!!$ call sethalon_new('v+',dlxv,dlxv_h, 2)
!!$ call sethalon_new('u',ubbl,ubbl_h, 2)
!!$ call sethalon_new('v',vbbl,vbbl_h, 2)
!!$
!!$ aux=real(kdwubbl)
!!$ call sethalon_new('u+',aux,aux_h, 2)
!!$ kdwubbl_h=int(aux_h)
!!$
!!$ aux=real(kupubbl)
!!$ call sethalon_new('u+',aux,aux_h, 2)
!!$ kupubbl_h=int(aux_h)
!!$
!!$ aux=real(kupvbbl)
!!$ call sethalon_new('v+',aux,aux_h, 2)
!!$ kupvbbl_h=int(aux_h)
!!$
!!$ aux=real(kdwvbbl)
!!$ call sethalon_new('v+',aux,aux_h, 2)
!!$ kdwvbbl_h=int(aux_h)
!!$
!!$
!!$
!!$ do k=1,ke
!!$ call sethalon_new('p',trf(:,:,k),trf_h(:,:,k),2)
!!$ call sethalon_new('p',weto(:,:,k),weto_h(:,:,k),2)
!!$ call sethalon_new('p',svo(:,:,k),svo_h(:,:,k),2)
!!$ call sethalon_new('u+',dduo(:,:,k),dduo_h(:,:,k), 2)
!!$ call sethalon_new('u',uko(:,:,k),uko_h(:,:,k), 2)
!!$ call sethalon_new('v+',ddue(:,:,k),ddue_h(:,:,k), 2)
!!$ call sethalon_new('v',vke(:,:,k),vke_h(:,:,k), 2)
!!$ enddo
!!$
!!$ do k=1,kep
!!$ call sethalon_new('p',woh2(:,:,k),woh2_h(:,:,k),2)
!!$ enddo
!!$
!!$if ( iocad == 1 ) then
!!$
!!$! #slo# - initialise local variables
!!$#ifndef TRACER_OMP
!!$!$OMP WORKSHARE
!!$#endif
!!$   s2o_h(:,:,:)=0.
!!$   trp_h(:,:,:)=0.
!!$   trm_h(:,:,:)=0.
!!$   wtp_h(:,:,:)=0.
!!$   wtm_h(:,:,:)=0.
!!$#ifndef TRACER_OMP
!!$!$OMP END WORKSHARE
!!$#endif
!!$
!!$!
!!$! Vertical, Zonal, and Meridional Transports for iocad == 1 - Pure Upwind
!!$!
!!$
!!$!      VERTICAL TRANSPORTS
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$    DO k=1,ke
!!$  DO J=-1,JE+2
!!$
!!$      zzsurf = 1.
!!$      if (k == 1) zzsurf=0.
!!$      DO I=-1,IE+2
!!$        UP=0.5*DT * area_h(i,j) * (woh2_H(I,J,K)+ABS(woh2_h(I,J,K)))*zzsurf
!!$        UM=0.5*DT * area_h(i,j) * (ABS(woh2_h(I,J,K+1))-woh2_h(I,J,K+1))
!!$        TRP_h(I,J,K)=UP*TRF_h(I,J,K)
!!$        TRM_h(I,J,K)=UM*TRF_h(I,J,K)
!!$        WTP_h(I,J,K)=UP
!!$        WTM_h(I,J,K)=UM
!!$      END DO
!!$    END DO
!!$  END DO
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$    DO K=1,KE
!!$  DO J=-1,JE+2
!!$
!!$      DO I=-1,IE+2
!!$        IF(WETO_h(I,J,K).GT.0.5) THEN
!!$          RN=svo_h(I,J,K)+WTP_h(I,J,K+1)-WTP_h(I,J,K)-WTM_h(I,J,K)+WTM_h(I,J,K-1)
!!$          TRF_h(I,J,K)=(TRF_h(I,J,K)*svo_h(I,J,K)+TRP_h(I,J,K+1)-TRP_h(I,J,K)      &
!!$               -TRM_h(I,J,K)+TRM_h(I,J,K-1))/RN
!!$          s2o_h(I,J,K)=RN
!!$        ENDIF
!!$      END DO
!!$    END DO
!!$  END DO
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$!  CALL bounds_exch(1,'p',TRF,'ocadpo_trf 3')   ! #slo#: not necessary?
!!$!  CALL bounds_exch(1,'p',s2o,'ocadpo_trf 4')   ! #slo#: not necessary?
!!$!$OMP END SINGLE
!!$!#endif
!!$
!!$!      TRANSPORTS IN X-DIRECTION
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$    DO K=1,KE
!!$     DO J=1,JE+2
!!$      DO I=0,IE+1
!!$
!!$        UP=0.5*DT*DDUO_h(I,J,K)*DLYU_h(I,J)*(UKO_h(I,J,K)+ABS(UKO_h(I,J,K)))
!!$
!!$        UM=0.5*DT*DDUO_h(I-1,J,K)*(ABS(UKO_h(I-1,J,K))-UKO_h(I-1,J,K))         &
!!$                *DLYU_h(I-1,J)
!!$
!!$        TRP_h(I,J,K)=UP*TRF_h(I,J,K)
!!$        TRM_h(I,J,K)=UM*TRF_h(I,J,K)
!!$        WTP_h(I,J,K)=UP
!!$        WTM_h(I,J,K)=UM
!!$      END DO
!!$    END DO
!!$  END DO
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$!  CALL bounds_exch(1,'u+',TRP,'ocadpo_trf 5')
!!$!  CALL bounds_exch(1,'u+',TRM,'ocadpo_trf 6')
!!$!  CALL bounds_exch(1,'u+',WTP,'ocadpo_trf 7')
!!$!  CALL bounds_exch(1,'u+',WTM,'ocadpo_trf 8')
!!$!$OMP END SINGLE
!!$!#endif
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$    DO K=1,KE
!!$     DO J=-1,JE+2
!!$      DO I=0,IE+1
!!$        IF(WETO_h(I,J,K).GT.0.5) THEN
!!$          RN=s2o_h(I,J,K)+WTP_h(I-1,J,K)-WTP_h(I,J,K)-WTM_h(I,J,K)+WTM_h(I+1,J,K)
!!$          TRF_h(I,J,K)=(TRF_h(I,J,K)*s2o_h(I,J,K)+TRP_h(I-1,J,K)-TRP_h(I,J,K)      &
!!$               -TRM_h(I,J,K)+TRM_h(I+1,J,K))/RN
!!$          s2o_h(I,J,K)=RN
!!$        ENDIF
!!$      ENDDO
!!$    ENDDO
!!$  END DO
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$!  CALL bounds_exch(1,'p',TRF,'ocadpo_trf 9')    ! #slo#: not necessary?
!!$!  CALL bounds_exch(1,'p',s2o,'ocadpo_trf 10')   ! #slo#: not necessary?
!!$!$OMP END SINGLE
!!$!#endif
!!$
!!$!      TRANSPORTS IN Y-DIRECTION
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$    DO K=1,KE
!!$     DO J=0,JE+1
!!$      DO I=-1,IE+2
!!$        UP=0.5*DT*DDUE_h(I,J-1,K)*DLXV_h(I,J-1)                              &
!!$             *(VKE_h(I,J-1,K)+ABS(VKE_h(I,J-1,K)))
!!$        UM=0.5*DT*DDUE_h(I,J,K)*DLXV_h(I,J)*(ABS(VKE_h(I,J,K))-VKE_h(I,J,K))
!!$
!!$        TRP_h(I,J,K)=UP*TRF_h(I,J,K)
!!$        TRM_h(I,J,K)=UM*TRF_h(I,J,K)
!!$        WTP_h(I,J,K)=UP
!!$        WTM_h(I,J,K)=UM
!!$      END DO
!!$    END DO
!!$  END DO
!!$
!!$!
!!$! Vertical, Zonal, and Meridional Transports for iocad == 3, 4
!!$!
!!$
!!$ELSE     ! i.e. iocad==3 or iocad==4
!!$
!!$!      VERTICAL TRANSPORTS
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP WORKSHARE
!!$#endif
!!$   s2o_h(:,:,:)=0.
!!$#ifndef TRACER_OMP
!!$!$OMP END WORKSHARE
!!$#endif
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$!here JBlock
!!$    DO K=1,KE
!!$       IF (k.EQ.2)THEN
!!$          DO J=-1,JE+2
!!$             DO I=-1,IE+2
!!$                UM=0.5*DT * area_h(i,j) * (ABS(woh2_h(I,J,2))-woh2_h(I,J,2))
!!$                TRP_h(I,J,1)=0.
!!$                TRM_h(I,J,0)=0.
!!$                TRM_h(I,J,1)=UM*TRF_h(I,J,1)
!!$                WTP_h(I,J,1)=0.
!!$                WTM_h(I,J,0)=0.
!!$                WTM_h(I,J,1)=UM
!!$
!!$                UP=0.5*DT * area_h(i,j) * (woh2_h(I,J,KE)+ABS(woh2_h(I,J,KE)))
!!$                UM=0.5*DT * area_h(i,j) * (ABS(woh2_h(I,J,KE+1))-woh2_h(I,J,KE+1))
!!$                TRP_h(I,J,KE)=UP*TRF_h(I,J,KE)
!!$                TRM_h(I,J,KE)=UM*TRF_h(I,J,KE)
!!$                WTP_h(I,J,KE)=UP
!!$                WTM_h(I,J,KE)=UM
!!$                WTP_h(I,J,KEP)=0.
!!$                TRP_h(I,J,KEP)=0.
!!$             END DO
!!$          ENDDO
!!$       ENDIF
!!$
!!$       IF (k.ge.2 .AND. K.Lt.KE) THEN
!!$
!!$          KLO=MIN(K+1,KE)
!!$          KUP=MAX(K-1,1)
!!$          DO J=-1,JE+2
!!$             DO I=-1,IE+2
!!$
!!$                ABL=ABS(TRF_h(I,J,KUP)-TRF_h(I,J,KLO))
!!$                ZWABL=ABS(TRF_h(I,J,KUP)+TRF_h(I,J,KLO)-2.*TRF_h(I,J,K))
!!$                UP=0.5*DT * area_h(i,j) * (woh2_h(I,J,K)+ABS(woh2_h(I,J,K)))
!!$                UM=0.5*DT * area_h(i,j) * (ABS(woh2_h(I,J,K+1))-woh2_h(I,J,K+1))
!!$                SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
!!$
!!$#ifdef SMOADV
!!$                SCH=MIN(WETO_h(I,J,KLO),SCHA)
!!$#else
!!$                SCH=MIN(WETO_h(I,J,KLO),SCHA*svo_h(I,J,K)/(UP+UM+1.E-20))
!!$#endif
!!$
!!$                SCHI=1.-SCH
!!$
!!$                TRP_h(I,J,K)=UP*(SCHI*TRF_h(I,J,K)+SCH*0.5*(TRF_h(I,J,K)+TRF_h(I,J,KUP)))
!!$                TRM_h(I,J,K)=UM*(SCHI*TRF_h(I,J,K)+SCH*0.5*(TRF_h(I,J,K)+TRF_h(I,J,KLO)))
!!$                WTP_h(I,J,K)=UP
!!$                WTM_h(I,J,K)=UM
!!$             END DO
!!$          END DO
!!$       ENDIF
!!$    END DO
!!$    do k=1,ke
!!$       DO J=-1,JE+2
!!$          DO I=-1,IE+2
!!$             IF(WETO_h(I,J,K).GT.0.5) THEN
!!$                RN=svo_h(I,J,k)+WTP_h(I,J,K+1)-WTP_h(I,J,k)-WTM_h(I,J,k)+WTM_h(I,J,k-1)
!!$                TRF_h(I,J,K)=(TRF_h(I,J,K)*svo_h(I,J,k)+TRP_h(I,J,k+1)-TRP_h(I,J,k)      &
!!$                     -TRM_h(I,J,k)+TRM_h(I,J,K-1))/RN
!!$                s2o_h(I,J,K)=RN
!!$             ENDIF
!!$          END DO
!!$       END DO
!!$    END DO
!!$    !here jblock
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$!  CALL bounds_exch(1,'p',TRF,'ocadpo_trf 3')   ! #slo#: not necessary?
!!$!  CALL bounds_exch(1,'p',s2o,'ocadpo_trf 4')   ! #slo#: not necessary?
!!$!$OMP END SINGLE
!!$!#endif
!!$#ifdef _PROFILE
!!$  CALL trace_stop ('mo_adpotrf 1', 6)
!!$#endif
!!$!      TRANSPORTS IN X-DIRECTION
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$    DO K=1,KE
!!$  DO J=0,JE+1
!!$      DO I=0,IE+1
!!$        b1u  = 0.
!!$        b1um = 0.
!!$         if ( iocad == 4 ) then
!!$           IF (KDWUBBL_h(I,J).EQ.K) b1u=UBBL_h(I,J)
!!$           IF (KUPUBBL_h(I,J).EQ.K) b1u=-UBBL_h(I,J)
!!$           IF (KDWUBBL_h(I-1,J).EQ.K) b1um=UBBL_h(I-1,J)
!!$           IF (KUPUBBL_h(I-1,J).EQ.K) b1um=-UBBL_h(I-1,J)
!!$         endif
!!$        ABL=ABS(TRF_h(I+1,J,K)-TRF_h(I-1,J,K))
!!$        ZWABL=ABS(TRF_h(I+1,J,K)+TRF_h(I-1,J,K)-2.*TRF_h(I,J,K))
!!$
!!$        UP=0.5*DT*DDUO_h(I,J,K)*DLYU_h(I,J)*(UKO_h(I,J,K)+ABS(UKO_h(I,J,K)))     &
!!$             + 0.5*DT*(b1u+ABS(b1u))
!!$
!!$        UM=0.5*DT*DDUO_h(I-1,J,K)*(ABS(UKO_h(I-1,J,K))-UKO_h(I-1,J,K))         &
!!$             *DLYU_h(I-1,J)                                                &
!!$             + 0.5*DT*(ABS(b1um)-b1um)
!!$
!!$        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO_h(I,J,K)                &
!!$             *WETO_h(I+1,J,K)*WETO_h(I-1,J,K)
!!$#ifdef SMOADH
!!$        SCH=MIN(1.,SCHA)
!!$#else
!!$        SCH=MIN(1.,SCHA*s2o_h(I,J,K)/(UP+UM+1.E-20))
!!$#endif
!!$
!!$        SCHI=1.-SCH
!!$
!!$        TRP_h(I,J,K)=UP*(SCHI*TRF_h(I,J,K)+SCH*0.5*(TRF_h(I,J,K)+TRF_h(I+1,J,K)))
!!$        TRM_h(I,J,K)=UM*(SCHI*TRF_h(I,J,K)+SCH*0.5*(TRF_h(I,J,K)+TRF_h(I-1,J,K)))
!!$        WTP_h(I,J,K)=UP
!!$        WTM_h(I,J,K)=UM
!!$      END DO
!!$    END DO
!!$
!!$  END DO
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$!  CALL bounds_exch(1,'u+',TRP,'ocadpo_trf 5')
!!$!  CALL bounds_exch(1,'u+',TRM,'ocadpo_trf 6')
!!$!  CALL bounds_exch(1,'u+',WTP,'ocadpo_trf 7')
!!$!  CALL bounds_exch(1,'u+',WTM,'ocadpo_trf 8')
!!$!$OMP END SINGLE
!!$!#endif
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$    DO K=1,KE
!!$  DO J=0,JE+1
!!$
!!$      DO I=0,IE+1
!!$        IF(WETO_h(I,J,K).GT.0.5) THEN
!!$          RN=s2o_h(I,J,K)+WTP_h(I-1,J,K)-WTP_H(I,J,K)-WTM_h(I,J,K)+WTM_h(I+1,J,K)
!!$          TRF_h(I,J,K)=(TRF_h(I,J,K)*s2o_h(I,J,K)+TRP_H(I-1,J,K)-TRP_h(I,J,K)      &
!!$               -TRM_h(I,J,K)+TRM_h(I+1,J,K))/RN
!!$          s2o_h(I,J,K)=RN
!!$        ENDIF
!!$      ENDDO
!!$    ENDDO
!!$  END DO
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$!  CALL bounds_exch(1,'p',TRF,'ocadpo_trf 9')    ! #slo#: not necessary?
!!$!  CALL bounds_exch(1,'p',s2o,'ocadpo_trf 10')   ! #slo#: not necessary?
!!$!$OMP END SINGLE
!!$!#endif
!!$
!!$!      TRANSPORTS IN Y-DIRECTION
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$  DO K=1,KE
!!$     DO J=0,JE+1
!!$
!!$!    hh: is this really needed ?
!!$        if ( lbounds_exch_tp .and. p_joff.eq.0 ) then
!!$           do i=-1,ie+2
!!$              trf_h(i,1,k)=2.*trf_h(i,2,k)-trf_h(i,3,k)
!!$           enddo
!!$        endif
!!$
!!$        DO I=0,IE+1
!!$           b1v = 0.
!!$           b1vm = 0.
!!$           if ( iocad == 4 ) then
!!$              IF (K.EQ.KDWVBBL_h(I,J)) b1v=VBBL_h(I,J)
!!$              IF (K.EQ.KUPVBBL_h(I,J)) b1v=-VBBL_h(I,J)
!!$              IF (K.EQ.KDWVBBL_h(I,J-1)) b1vm=VBBL_h(I,J-1)
!!$              IF (K.EQ.KUPVBBL_h(I,J-1)) b1vm=-VBBL_h(I,J-1)
!!$           endif
!!$
!!$        ABL=ABS(TRF_h(I,J-1,K)-TRF_h(I,J+1,K))
!!$        ZWABL=ABS(TRF_h(I,J-1,K)+TRF_h(I,J+1,K)-2.*TRF_h(I,J,K))
!!$        UP=0.5*DT*DDUE_h(I,J-1,K)*DLXV_h(I,J-1)                              &
!!$             *(VKE_h(I,J-1,K)+ABS(VKE_h(I,J-1,K)))                           &
!!$                + 0.5*DT*(b1vm+ABS(b1vm))
!!$        UM=0.5*DT*DDUE_h(I,J,K)*DLXV_h(I,J)*(ABS(VKE_h(I,J,K))-VKE_h(I,J,K))     &
!!$                + 0.5*DT*(ABS(b1v)-b1v)
!!$        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO_h(I,J,K)                &
!!$             *WETO_h(I,J-1,K)*WETO_h(I,J+1,K)
!!$#ifdef SMOADH
!!$        SCH=MIN(1.,SCHA)
!!$#else
!!$        SCH=MIN(1.,SCHA*s2o_h(I,J,K)/(UP+UM+1.E-20))
!!$#endif
!!$
!!$        SCHI=1.-SCH
!!$        TRP_h(I,J,K)=UP*(SCHI*TRF_h(I,J,K)+SCH*0.5*(TRF_h(I,J,K)+TRF_h(I,J-1,K)))
!!$        TRM_h(I,J,K)=UM*(SCHI*TRF_h(I,J,K)+SCH*0.5*(TRF_h(I,J,K)+TRF_h(I,J+1,K)))
!!$        WTP_h(I,J,K)=UP
!!$        WTM_h(I,J,K)=UM
!!$      END DO
!!$    END DO
!!$  END DO
!!$
!!$ENDIF    ! iocad
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$!  CALL bounds_exch(1,'vv',TRP,'ocadpo_trf 11')
!!$!  CALL bounds_exch(1,'vv',TRM,'ocadpo_trf 12')
!!$!  CALL bounds_exch(1,'vv',WTP,'ocadpo_trf 13')
!!$!  CALL bounds_exch(1,'vv',WTM,'ocadpo_trf 14')
!!$!$OMP END SINGLE
!!$!#endif
!!$
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP DO
!!$#endif
!!$    DO K=1,KE
!!$  DO J=0,JE+1
!!$
!!$      DO I=0,IE+1
!!$        IF(WETO_h(I,J,K).GT.0.5) THEN
!!$
!!$          RN=s2o_h(I,J,K)+WTP_h(I,J+1,K)-WTP_h(I,J,K)-WTM_h(I,J,K)+WTM_h(I,J-1,K)
!!$          TRF_h(I,J,K)=(TRF_h(I,J,K)*s2o_h(I,J,K)+TRP_h(I,J+1,K)-TRP_h(I,J,K)      &
!!$               -TRM_h(I,J,K)+TRM_h(I,J-1,K))/RN
!!$          s2o_h(I,J,K)=RN
!!$        ENDIF
!!$      ENDDO
!!$    ENDDO
!!$  END DO
!!$
!!$!#ifdef bounds_exch_save
!!$!$OMP SINGLE
!!$!  CALL bounds_exch(1,'p',TRF,'ocadpo_trf 15')
!!$!$OMP END SINGLE
!!$!#endif
!!$
!!$ trf(:,:,:)=trf_h(1:ie,1:je,:)
!!$
!!$
!!$#ifndef TRACER_OMP
!!$!$OMP END PARALLEL
!!$#endif
!!$
!!$END SUBROUTINE ocadpo_trf_halo



END MODULE MO_ADPO
