!> Calulation of horizontal and vertical diffusion of tracers.
!> @author Ernst Maier-Reimer, Uwe Mikolajewicz, ... 1999
!> @date Last modified on 09.03.2009 by Helmuth Haak
!> @todo Documentation
!> This module holds subroutines for horizontal and vertical diffusion of tracers.
!! Included subroutines are:
!! octdiff_base  : Calulation of tracer-independent matrices for horizontal and vertical diffusion
!! octdiff_trf   : Calulation of horizontal and vertical diffusion


MODULE mo_diffusion

  USE mo_kind, ONLY: dp, wp
  USE mo_param1, ONLY: ie, ie1, je, je1, ke
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1, ONLY: lisopyc, almzer, area, amsue, amsuo, ddpo, ddue, dduo, di, &
       dlxp, dlxu, dlxv, dlyp, dlyu, dlyv, dlxui, dlyvi, dt, dvo, dzw, rhoo, &
       sicsno, sictho, stabio, weto, zo
!  USE mo_levitus, ONLY:
  USE mo_planetary_constants, ONLY: rhoicwa, rhosnwa


  IMPLICIT NONE

  REAL(wp), ALLOCATABLE :: vol_term(:,:,:)     ! Local Active Volume of Cell

  REAL(wp), ALLOCATABLE :: sloploy(:,:,:)     &! Isopycnal Density Slope (upper left)
                      ,slopluy(:,:,:)     &! Isopycnal Density Slope (lower left)
                      ,sloproy(:,:,:)     &! Isopycnal Density Slope (upper right)
                      ,slopruy(:,:,:)     &! Isopycnal Density Slope (lower right)
                      ,sloplox(:,:,:)     &! Isopycnal Density Slope (upper left)
                      ,sloplux(:,:,:)     &! Isopycnal Density Slope (lower left)
                      ,sloprox(:,:,:)     &! Isopycnal Density Slope (upper right)
                      ,sloprux(:,:,:)      ! Isopycnal Density Slope (lower right)

  REAL(wp), ALLOCATABLE :: xcoeff_lo(:,:,:),zcoeff_lox(:,:,:)   &!constrained slope
                      ,xcoeff_lu(:,:,:),zcoeff_lux(:,:,:)   &!constrained slope
                      ,xcoeff_ro(:,:,:),zcoeff_rox(:,:,:)   &!constrained slope
                      ,xcoeff_ru(:,:,:),zcoeff_rux(:,:,:)   &!constrained slope
                      ,ycoeff_lo(:,:,:),zcoeff_loy(:,:,:)   &!constrained slope
                      ,ycoeff_lu(:,:,:),zcoeff_luy(:,:,:)   &!constrained slope
                      ,ycoeff_ro(:,:,:),zcoeff_roy(:,:,:)   &!constrained slope
                      ,ycoeff_ru(:,:,:),zcoeff_ruy(:,:,:)    !constrained slope

  REAL(wp), ALLOCATABLE :: tridsy(:,:,:,:) !  tridiagonal matrix used in the implicit vertical diffusion scheme

  REAL(wp), ALLOCATABLE :: zbott(:,:),zsurf(:,:) ! help fields

  INTEGER, ALLOCATABLE :: ko(:,:),ku(:,:) ! help fields

  REAL(wp), ALLOCATABLE :: xflux(:,:,:),yflux(:,:,:) ! harmonic horizontal diffusion fluxes

  REAL(wp), ALLOCATABLE :: aulx(:,:),auly(:,:) ! biharmonic horizontal diffusion coefficients


  REAL(wp) :: aulapts, ah00

CONTAINS

!> Calulation of tracer-independent matrices for horizontal and vertical diffusion
!> @author Ernst Maier-Reimer, Uwe Mikolajewicz, ... 1999
!> @date Modified by R. Smith, 2004-09-27: tracer-independent matrices calculated separately
!> @date Last modified 2009-02-23 by Stephan Lorenz
!> @todo Specific documentation of ...
!> @todo The Gent and McWilliams (1990) eddy induced tracer transport (routine ocjitr) and
!> the isoneutral diffusion (octdiff_base and octdiff_trf) could be combined following the skew
!> flux method by Griffies (1998), which is not implemented in MPIOM
!>
!!
!> This module calculates the tracer-independent matrices for the calculation of
!! horizontal and vertical diffusion of tracer fields, done in routine octdiff_trf.f90.
!! Isopycnal density slopes for isoneutral diffusion (cpp-key ISOPYK) or local coefficient
!! for explicit harmonic and biharmonic diffusion are calculated.
!!
!! The implementation follows Griffies et al. (1998). They utilized the concept of coordinate rotation
!! along the isoneutral direction of the symmetric isoneutral diffusion tensor of Redi (1982)
!! by applying the small slope approximation tensor of Gent and McWilliams (1990).
!!
!! Vertical tracer diffusion follows this equation,
!! \f[
!!\frac{\partial T}{\partial t}=\frac{\partial}{\partial z}\big( D_v
!!\frac{\partial T}{\partial z}  \big)
!!\f]
!! \f$T\f$ is the tracer, \f$D_v\f$ is the diffusivity.
!! The discretized form looks as follows:
!! \f[
!!T^t_k - T^{t-1}_k = \tau \big( D_{vk} \frac{T^n_{k-1} - T^t_k} {\Delta
!!z_{k-1}\Delta z} \big) - \tau\big(
!!D_{v(k+1)}\frac{T^t_k-T^t_{k+1}}{\Delta z_{k+1}\Delta z } \big)
!!\f]
!!Here, \f$t\f$ is time increment, \f\tau\f$ is the time step.
!!We can rewrite this as
!!\f[
!!T^{t-1}_k = \underbrace {\big( 1 + \frac{\tau D_{vk} }{\Delta z_{k-1}\Delta z}
!!+ \frac{\tau D_{v(k+1)} }{\Delta z_{k+1}\Delta z}}_{\alpha_k} \big) T^t_k
!!\underbrace {-\frac{\tau D_{vk} }{\Delta z_{k-1}\Delta z}}_{\beta_k}
!!T^t_{k-1}
!!\underbrace{-\frac{\tau D_{v(k+1)} }{\Delta z_{k+1}\Delta
!!z}}_{\gamma_k} T^t_{k+1}
!!\f]
!!Where
!!\f[
!!\alpha_k=1-\beta_k-\gamma_k
!!\f]
!!This takes a form of a matrix k*k that can be solved from bottom to top. The matrix is
!!tridiagonal, and the coefficients of the diagones are \f$\beta\f$, \f$\alpha\f$,
!!and \f$\gamma\f$, which are described above. After some algebra
!!with the matrix, the \f$\beta\f$ are eliminated, by multiplying
!!with
!!\f[
!!\frac{\beta _k}{1-\beta _{k-1} - \gamma _{k-1}}
!!\f]
!!After that the \f$T_k^t\f$ can be found.
!!
!! References:
!!  Griffies et al. 1998; JPO 28:805-830
!!  Griffies 1998; JPO 28:831-83?
!!  Gent and McWilliams 1990; JPO 20:150-155
!!  Redi 1982; JPO 12:1154-1158
!!
SUBROUTINE octdiff_base
  !     Old comments:
  !     -------------
  !
  !     COMPUTES DIFFUSION OF TEMPERATURE THO AND SALINITY SAO
  !
  !     UWE MIKOLAJEWICZ 12/99
  !
  !     VERTICAL DIFFUSION IMPLICITE
  !     HORIZONTAL DIFFUSION BOTH HARMONIC AND BIHARMONIC, BOTH EXPLICITE
  !                                AH00         AULAPTS
  !     ACTUAL COEFFICIENTS SPATIALLY VARIABLE WITH RESOLUTION
  !                               AHX,AHY      AULX,AULY
  !
  !     VARIABLE THICKNESS OF SURFACE LAYER NOW INCLUDED
  !
  !     Changes R. Johanni, 2003-11-10:
  !     DSLOPX, DSLOPY, DVSLOP were calculated but nowhere used -> removed
  !
  !     Changes R. Smith, 2004-09-27
  !     tracer-independent matrices calculated separately first for use
  !     with multiple tracers

  REAL(dp) :: ahx(ie,je), ahy(ie,je), ahxd(ie,je), ahyd(ie,je), zzsurf(ke)

  INTEGER :: i, j, k
  REAL(dp) :: otmp
  REAL(dp) :: scale, slcut, stabmin, cutdel
  REAL(dp) ::rho_diff, xterm, zterm, yterm
  !
  DO k = 1, ke
    zzsurf(k) = 0.0_dp
  ENDDO
  zzsurf(1) = 1.0_dp
  !
  !     Minimum for Stability
  !

!!      Vertical stability STABIO>0. - the vertical density gradient -
!!      is calculated in routine octher.f90
!!      In case of instability no isopycnal slopes are calculated

  stabmin = 1.e-7_dp   !  minimum for STABIO - the vertical stability
  cutdel  = 0.02_dp    !  factor to constrain slopes not to rise above one layer (?)

!$OMP PARALLEL PRIVATE(i,j,k,scale,slcut,rho_diff,xterm,zterm,yterm,otmp)
  !
  IF(ah00 > almzer)THEN

  !     Local Diffusion Coefficients:
  !      - AHX  and AHY  are for harmonic diffusion
  !      - AHXD and AHYD are for isopycnal diffusion
  !
!$OMP DO
     DO j = 1, je
        DO i = 1, ie
           ahx(i,j)  = ah00*MAX(dlxu(i,j),dlyu(i,j))
           ahy(i,j)  = ah00*MAX(dlxv(i,j),dlyv(i,j))      !emr vpoints
           ahxd(i,j) = ahx(i,j)*dt*dlyu(i,j)
           ahyd(i,j) = ahy(i,j)*dt*dlxv(i,j)
        ENDDO
     ENDDO

!!$!$OMP SINGLE
!!$  call bounds_exch(1,'p',zo)
!!$  call bounds_exch(1,'p',sictho)
!!$  call bounds_exch(1,'p',sicsno)
!!$  call bounds_exch(1,'p',weto)
!!$!$OMP END SINGLE


!$OMP DO
     DO k = 1, ke

       DO j = 1, je
         ko(j,k) = MAX(k-1,1)
         ku(j,k) = MIN(k+1,ke)
         zsurf(j,k) = REAL(k-ko(j,k),dp)
         zbott(j,k) = REAL(ku(j,k)-k,dp)
       ENDDO


!!      Local Active Volume of Cell
!!       - subtract ice and snow in surface layer (ZZSURF)
!!       - ZO is surface layer including surface elevation

       DO j = 1, je
         DO i = 1, ie
!             vol_term(i,j,k) = dlxp(i,j)*dlyp(i,j)*(zzsurf(k)*(zo(i,j) &
!                -rhoicwa*sictho(i,j)-rhosnwa*sicsno(i,j))            &
!                +ddpo(i,j,k)+(1.0_dp-weto(i,j,k)))
           vol_term(i,j,k) = area(i,j)*(zzsurf(k)*(zo(i,j) &
                -rhoicwa*sictho(i,j)-rhosnwa*sicsno(i,j))            &
                +ddpo(i,j,k)+(1.0_dp-weto(i,j,k)))

         ENDDO
       ENDDO


!!      Calculation of Isopycnal Density Slopes
!!       - tracer independent values
!!       - discretized slopes (SLOPLOX/SLOPLUX/SLOPROX/SLOPRUX) are calculated
!!         following Eq. (31) of Griffies et al.(1998)
!!       - 4 contributions of different density slopes at the boundary faces are considered,
!!         see x-z grid arrangement (Fig. 3) and grid stencil (Fig. 9) of Griffies et al. (1998)
!!       - constrain slopes to a maximum of one layer through SCALE and SLCUT (?)


!!      Calculation in X Direction - tracer independent values

       IF (lisopyc) THEN

         DO j = 2, je-1
           DO i = 1, ie1
             IF (amsuo(i,j,k) > 0.0_dp) THEN

                slcut = cutdel*dzw(k)**2/(ahx(i,j)*dt)
                rho_diff = (rhoo(i+1,j,k)-rhoo(i,j,k))/dlxu(i,j)
                xterm = 0.25_dp*ahxd(i,j)*dduo(i,j,k)
                zterm = 0.25_dp*dt*ahx(i,j)

                !  contribution from left/up
                sloplox(i,j,k) = rho_diff/MAX(stabio(i,j,k),stabmin)
                scale = 1.0_dp/MAX(1.0_dp,(sloplox(i,j,k)/slcut)**2)
                xcoeff_lo(i,j,k)  = xterm*zsurf(j,k)*scale
                zcoeff_lox(i,j,k) = zterm*zsurf(j,k)*scale*sloplox(i,j,k)

                !  contribution from left/down
                sloplux(i,j,k) = rho_diff/MAX(stabio(i,j,ku(j,k)),stabmin)
                scale = 1.0_dp/MAX(1.0_dp,(sloplux(i,j,k)/slcut)**2)
                xcoeff_lu(i,j,k)  = xterm*zbott(j,k)*weto(i,j,ku(j,k))*scale
                zcoeff_lux(i,j,k) = zterm*zbott(j,k)*scale*sloplux(i,j,k)*      &
                     weto(i,j,ku(j,k))

                !  contribution from right/up
                sloprox(i,j,k) = rho_diff/MAX(stabio(i+1,j,k),stabmin)
                scale = 1.0_dp/MAX(1.0_dp,(sloprox(i,j,k)/slcut)**2)
                xcoeff_ro(i,j,k)  = xterm*zsurf(j,k)*scale
                zcoeff_rox(i,j,k) = zterm*zsurf(j,k)*scale*sloprox(i,j,k)
                !
                !  contribution from right/down
                sloprux(i,j,k) = rho_diff/MAX(stabio(i+1,j,ku(j,k)),stabmin)
                scale = 1.0_dp/MAX(1.0_dp,(sloprux(i,j,k)/slcut)**2)
                xcoeff_ru(i,j,k)  = xterm*zbott(j,k)*weto(i+1,j,ku(j,k))*scale
                zcoeff_rux(i,j,k) = zterm*zbott(j,k)*scale*sloprux(i,j,k)*     &
                     weto(i+1,j,ku(j,k))
             ENDIF
          ENDDO
       ENDDO

     ELSE   !.NOT. lisopyc

       DO j = 1, je
         DO i = 1, ie
           IF (amsuo(i,j,k) > 0.0_dp) THEN
             !!      Calculation of Harmonic Horizontal Diffusion Flux
             xflux(i,j,k) = ahxd(i,j)*dduo(i,j,k)
           ENDIF
         ENDDO
       ENDDO

     ENDIF ! lisopyc

!!      Calculation in Y Direction - tracer independent values (see X Direction)

     IF (lisopyc) THEN
       DO j = 1, je1
          DO i = 2, ie1
             IF (amsue(i,j,k) > 0.0_dp) THEN

                slcut = cutdel*dzw(k)**2/(ahy(i,j)*dt)
                rho_diff = (rhoo(i,j+1,k)-rhoo(i,j,k))/dlyv(i,j)
                yterm = 0.25_dp*ahyd(i,j)*ddue(i,j,k)
                zterm = 0.25_dp*dt*ahy(i,j)
                !
                !  contribution from left/up
                sloploy(i,j,k) = rho_diff/MAX(stabio(i,j,k),stabmin)
                scale = 1.0_dp/MAX(1.0_dp,(sloploy(i,j,k)/slcut)**2)
                ycoeff_lo(i,j,k)  = yterm*zsurf(j,k)*scale
                zcoeff_loy(i,j,k) = zterm*zsurf(j,k)*scale*sloploy(i,j,k)
                !
                !  contribution from left/down
                slopluy(i,j,k) = rho_diff/MAX(stabio(i,j,ku(j,k)),stabmin)
                scale = 1.0_dp/MAX(1.0_dp,(slopluy(i,j,k)/slcut)**2)
                ycoeff_lu(i,j,k)  = yterm*zbott(j,k)*weto(i,j,ku(j,k))*scale
                zcoeff_luy(i,j,k) = zterm*zbott(j,k)*scale*slopluy(i,j,k)*     &
                     weto(i,j,ku(j,k))
                !
                !  contribution from right/up
                sloproy(i,j,k) = rho_diff/MAX(stabio(i,j+1,k),stabmin)
                scale = 1.0_dp/MAX(1.0_dp,(sloproy(i,j,k)/slcut)**2)
                ycoeff_ro(i,j,k)  = yterm*zsurf(j,k)*scale
                zcoeff_roy(i,j,k) = zterm*zsurf(j,k)*scale*sloproy(i,j,k)

                !
                !  contribution from right/down
                slopruy(i,j,k) = rho_diff/MAX(stabio(i,j+1,ku(j,k)),stabmin)
                scale = 1.0_dp/MAX(1.0_dp,(slopruy(i,j,k)/slcut)**2)
                ycoeff_ru(i,j,k)  = yterm*zbott(j,k)*weto(i,j+1,ku(j,k))*scale
                zcoeff_ruy(i,j,k) = zterm*zbott(j,k)*scale*slopruy(i,j,k)*     &
                 weto(i,j+1,ku(j,k))
             ENDIF
          ENDDO
       ENDDO

     ELSE !.NOT. lisopyc

       DO j = 1, je
         DO i = 1, ie
           IF (amsue(i,j,k) > 0.0_dp) THEN
             !!      Calculation of Harmonic Horizontal Diffusion Flux
             yflux(i,j,k) = ahyd(i,j)*ddue(i,j,k)
           ENDIF
         ENDDO
       ENDDO

     ENDIF ! lisopyc

   ENDDO   ! k-loop

   IF (.NOT. LISOPYC) THEN
     !$!$OMP SINGLE
     CALL bounds_exch(1,'u+',xflux)
     CALL bounds_exch(1,'v+',yflux)
     !$!$OMP END SINGLE
   ENDIF

 ENDIF ! ah00.gt.almzer

  ! VERTICAL DIFFUSION - tracer independent values
  !     IMPLICIT VERTICAL DIFFUSION

!!! The tridsy are the elements of the tridiagonal matrix, described
  !!above. tridsy(:,:,:,2) is the alpha, trdsy(:,:,:,1) is beta and
  !tridsy(:,:,:,3) is gamma.

!CDIR NOLOOPCHG
  DO k = 1,ke
    ! INCLUDE ACTUAL LEVEL THICKNESS
!$OMP DO
    DO j = 2,je1
      tridsy(:,j,k,1) = - dt*dvo(:,j,k)*weto(:,j,k)*di(k)         &
           /(ddpo(:,j,k)+zzsurf(k)*(zo(:,j)-rhoicwa*sictho(:,j)   &
           -rhosnwa*sicsno(:,j))+almzer)
      tridsy(:,j,k,3) = - dt*dvo(:,j,k+1) * di(k+1)               &
           /(ddpo(:,j,k)+zzsurf(k)*(zo(:,j)-rhoicwa*sictho(:,j)   &
           -rhosnwa*sicsno(:,j))+almzer)
      tridsy(:,j,k,2) = 1.0_dp - tridsy(:,j,k,1) - tridsy(:,j,k,3)
    END DO
  ENDDO
  !

!!$!CDIR NOLOOPCHG
!!$  DO k = 2,ke
!!$!$OMP DO
!!$    DO j = 2,je1
!!$      tridsy(:,j,k-1,1) = tridsy(:,j,k,1) / tridsy(:,j,k-1,2)
!!$      tridsy(:,j,k,2)   = tridsy(:,j,k,2) - tridsy(:,j,k-1,3)     &
!!$                         *tridsy(:,j,k,1) / tridsy(:,j,k-1,2)
!!$    END DO
!!$  ENDDO


!HH This workaround is suggested by NEC
!   it is needed to avoid different results
!   for openmp and non-openmp optimisation

!!!This is the elimination of the lower diagonal of the matrix.

!CDIR NOLOOPCHG
   DO k = 2,ke
!$OMP DO
     DO j = 2,je-1
       DO i=1,ie
       otmp = 1.0_wp / tridsy(i,j,k-1,2)
       tridsy(i,j,k-1,1) = tridsy(i,j,k,1) * otmp
       tridsy(i,j,k,2)   = tridsy(i,j,k,2) - tridsy(i,j,k-1,3)     &
                         * tridsy(i,j,k,1) * otmp
       END DO
     END DO
   END DO


  !
  ! BIHARMONIC DIFFUSION - tracer independent values
  !
  IF (aulapts > almzer) THEN
!$OMP DO
#ifndef AULREDSC
!CDIR NOUNROLL
    DO j = 1, je
      DO i = 1, ie
        aulx(i,j) = aulapts*dlxu(i,j)**4
        auly(i,j) = aulapts*dlyv(i,j)**4
      ENDDO
    ENDDO
#else
!CDIR NOUNROLL
    DO j = 1, je
      DO i = 1, ie
        aulx(i,j) = 1.e5_dp*aulapts*dlxu(i,j)**3
        auly(i,j) = 1.e5_dp*aulapts*dlyv(i,j)**3
      ENDDO
    ENDDO
#endif
  ENDIF
  !
!$OMP END PARALLEL
  !
END SUBROUTINE octdiff_base

!> Calulation of horizontal and vertical diffusion
!> @author Ernst Maier-Reimer, Uwe Mikolajewicz, ... 1999
!> @date Modified by R. Smith, 2004-09-27: tracer-independent matrices calculated separately
!> @date Last modified 2009-02-23 by Stephan Lorenz
!> @todo The Gent and McWilliams (1990) eddy induced tracer transport (routine ocjitr) and
!> the isoneutral diffusion (octdiff_base and octdiff_trf) could be combined following the skew
!> flux method by Griffies (1998), which is not implemented in MPIOM
!>
!!
!> This module calculates the components of the small slope approximation tensor for the calculation of
!! horizontal diffusion of tracer fields, using the density slopes calculated in routine octdiff_base.f90.
!!
!! The implementation follows Griffies et al. (1998). They utilized the concept of coordinate rotation
!! along the isoneutral direction of the symmetric isoneutral diffusion tensor of Redi (1982)
!! by applying the small slope approximation tensor of Gent and McWilliams (1990).
!!
!!
!! References:
!!  Griffies et al. 1998; JPO 28:805-830
!!  Griffies 1998; JPO 28:831-83?
!!  Gent and McWilliams 1990; JPO 20:150-155
!!  Redi 1982; JPO 12:1154-1158
SUBROUTINE octdiff_trf(trf)



  !     Old comments:
  !     -------------
  !
  !     COMPUTES DIFFUSION OF TEMPERATURE THO AND SALINITY SAO
  !
  !     UWE MIKOLAJEWICZ 12/99
  !
  !     VERTICAL DIFFUSION IMPLICITE
  !     HORIZONTAL DIFFUSION BOTH HARMONIC AND BIHARMONIC, BOTH EXPLICITE
  !                                AH00         AULAPTS
  !     ACTUAL COEFFICIENTS SPATIALLY VARIABLE WITH RESOLUTION
  !                                AHX,AHY      AULX,AULY
  !
  !     VARIABLE THICKNESS OF SURFACE LAYER NOW INCLUDED
  !
  !     Changes R. Johanni, 2003-11-10:
  !     DSLOPX, DSLOPY, DVSLOP were calculated but nowhere used -> removed
  !
  !     Changes R. Smith, 2004-09-27
  !     tracer-independent matrices calculated separately first for use
  !     with multiple tracers
  !
  !     Changes S. J. Lorenz, J.-O. Biesmann, 08/2007, optimisation for NEC

  INTEGER i, j, k

  REAL(wp) tv_diff, tfluz, zzsur0
  REAL(wp) tflux(ie,je,ke), tfluy(ie,je,ke)
  REAL(wp) trf(ie,je,ke), ten(ie,je,ke), t2o(ie,je,ke)


!#ifdef bounds_exch_save
!#slo#: Here - no bounds_exch necessary
!! !$OMP SINGLE
!!   CALL bounds_exch(1,'p',trf,'octdiff_trf 1')
!! !$OMP END SINGLE
!#endif

#ifndef TRACER_OMP
!$OMP PARALLEL PRIVATE(i,j,k,tv_diff,tfluz,zzsur0)

!$OMP WORKSHARE
#endif
  t2o(:,:,:)=trf(:,:,:)*vol_term(:,:,:)
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif

  IF (ah00.GT.almzer) THEN

#ifndef TRACER_OMP
!$OMP WORKSHARE
#endif
     tflux(:,:,:) = 0.0_wp
     tfluy(:,:,:) = 0.0_wp
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif

        IF (lisopyc) THEN

!!      Calculation of Isoneutral Diffusion Flux
!!       - tracer dependent values
!!       - discretized horizontal (TFLUX/TFLUY) and vertical (TFLUZ) flux components of small slope
!!         approximation tensor are calculated following Eqs. (30) and (33) of Griffies et al.(1998)
!!       - horizontal tracer gradient from i to i+1 is considered only
!!       - 4 contributions of vertical tracer gradients along isopycnal density slopes
!!         (calculated in octdiff_base) at the boundary faces are considered,
!!         see x-z grid arrangement (Fig. 3) and grid stencil (Fig. 9) of Griffies et al. (1998)

#ifndef TRACER_OMP
!$OMP DO
#endif

!!      Calculation in X Direction - tracer dependent values

        DO k=1,ke
           DO j=1,je
           DO i=1,ie1
              IF (amsuo(i, j, k) .GT. 0._wp) THEN
                 tv_diff=(trf(i+1,j,k)-trf(i,j,k))*dlxui(i,j)

                 !  flux contributions from left/up
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_lo(i,j,k)*(tv_diff  &
                      +sloplox(i,j,k)*                                &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*di(k))
                 ! FIXME: why multiply with 1.0?
                 tfluz=zcoeff_lox(i,j,k)*(sloplox(i,j,k)*             &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*                  &
                      area(i, j) * di(k) + 1._wp * area(i, j) * tv_diff)

                 t2o(i,j,ko(j,k))=t2o(i,j,ko(j,k))-tfluz
                 t2o(i,j,k)=t2o(i,j,k)+tfluz

                 !  flux contributions from left/down
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_lu(i,j,k)*(tv_diff  &
                      +sloplux(i,j,k)*                                &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*di(ku(j,k)))
                 ! FIXME: why multiply with 1.0?
                 tfluz=zcoeff_lux(i,j,k)*(sloplux(i,j,k)*             &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*                  &
                      area(i, j) * di(ku(j, k)) + 1._wp * area(i, j) * tv_diff)

                 t2o(i,j,k)=t2o(i,j,k)-tfluz
                 t2o(i,j,ku(j,k))=t2o(i,j,ku(j,k))+tfluz

              ENDIF
           END DO
        END DO
     END DO


#ifndef TRACER_OMP
!$OMP DO
#endif

!!      Calculation in X Direction - tracer dependent values, continued

        DO k=1,ke
           DO j=1,je
           DO i=1,ie1
              IF (amsuo(i, j, k) .GT. 0._wp) THEN
                 tv_diff=(trf(i+1,j,k)-trf(i,j,k))*dlxui(i,j)

                 !  flux contributions from right/up
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_ro(i,j,k)*(tv_diff &
                      +sloprox(i,j,k)*                               &
                      (trf(i+1,j,ko(j,k))-trf(i+1,j,k))*di(k))
                 ! FIXME: why multiply with 1.0?
                 tfluz=zcoeff_rox(i,j,k)*(sloprox(i,j,k)*            &
                      (trf(i+1,j,ko(j,k))-trf(i+1,j,k))*             &
                      area(i+1, j) * di(k) + 1._wp * area(i+1, j) * tv_diff)

                 t2o(i+1,j,ko(j,k))=t2o(i+1,j,ko(j,k))-tfluz
                 t2o(i+1,j,k)=t2o(i+1,j,k)+tfluz

                 !  flux contributions from right/down
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_ru(i,j,k)*(tv_diff &
                      +sloprux(i,j,k)*                               &
                      (trf(i+1,j,k)-trf(i+1,j,ku(j,k)))*di(ku(j,k)))
                 ! FIXME: why multiply with 1.0?
                 tfluz=zcoeff_rux(i,j,k)*(sloprux(i,j,k)*            &
                      (trf(i+1,j,k)-trf(i+1,j,ku(j,k)))*             &
                      area(i+1, j) * di(ku(j, k)) + 1._wp * area(i+1, j) * tv_diff)

                 t2o(i+1,j,k)=t2o(i+1,j,k)-tfluz
                 t2o(i+1,j,ku(j,k))=t2o(i+1,j,ku(j,k))+tfluz

              ENDIF
           END DO
        END DO
      END DO

    ELSE ! .NOT. lisopyc

!!      Calculation of Harmonic Horizontal Tracer Diffusion in X Direction

#ifndef TRACER_OMP
!$OMP DO
#endif
      DO k=1,ke
        DO j=1,je
          DO i=1,ie1
            IF (amsuo(i, j, k) .GT. 0._wp) THEN
              tv_diff=(trf(i+1,j,k)-trf(i,j,k))*dlxui(i,j)
              tflux(i,j,k)=tv_diff*xflux(i,j,k)
            ENDIF
          END DO
        END DO
      END DO

    ENDIF ! lisopyc

!$OMP SINGLE
     call bounds_exch(1,'u',tflux,'octdiff_trf 2 ')
!$OMP END SINGLE

#ifndef TRACER_OMP
!$OMP DO
#endif
     DO k = 1, ke
        DO j = 1, je
!hh not needed due to following bounds_exch
!              t2o(1,j,k)   = t2o(1,j,k)   + tflux(1,j,k)
!              t2o(ie,j,k)  = t2o(ie,j,k)  - tflux(ie1,j,k)
           DO i = 2, ie1
              t2o(i,j,k)   = t2o(i,j,k)   + tflux(i,j,k) - tflux(i-1,j,k)
           END DO
        END DO
     END DO

!$OMP SINGLE
     call bounds_exch(1,'p',t2o,'octdiff_trf 3 ')
!$OMP END SINGLE


     IF (lisopyc) THEN

!!      Calculation of Isopycnal Horizontal Diffusion in Y Direction

#ifndef TRACER_OMP
!$OMP DO
#endif
        DO k=1,ke
           DO j=1,je1

           DO i=1,ie
              IF (amsue(i, j, k) .GT. 0._wp) THEN
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))*dlyvi(i,j)

                 !  flux contributions from left/up
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_lo(i,j,k)*(tv_diff &
                      +sloploy(i,j,k)*                               &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*di(k))
                 ! FIXME: why multiply with 1.0?
                 tfluz=zcoeff_loy(i,j,k)*(sloploy(i,j,k)*            &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*                 &
                      area(i, j) * di(k) + 1._wp * area(i, j) * tv_diff)

                 t2o(i,j,ko(j,k))=t2o(i,j,ko(j,k))-tfluz
                 t2o(i,j,k)=t2o(i,j,k)+tfluz

                 !  flux contributions from left/down
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_lu(i,j,k)*(tv_diff &
                      +slopluy(i,j,k)*                               &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*di(ku(j,k)))
                 ! FIXME: why multiply with 1.0?
                 tfluz=zcoeff_luy(i,j,k)*(slopluy(i,j,k)*            &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*                 &
                      area(i, j) * di(ku(j, k)) + 1._wp * area(i, j) * tv_diff)

                 t2o(i,j,k)=t2o(i,j,k)-tfluz
                 t2o(i,j,ku(j,k))=t2o(i,j,ku(j,k))+tfluz

              END IF
           END DO
        END DO
     END DO


#ifndef TRACER_OMP
!$OMP DO
#endif
     DO k=1,ke
        DO j=1,je1

          DO i=1,ie
             IF (amsue(i, j, k) .GT. 0._wp) THEN
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))*dlyvi(i,j)

                 !  flux contributions from right/up
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_ro(i,j,k)*(tv_diff &
                      +sloproy(i,j,k)*                               &
                      (trf(i,j+1,ko(j,k))-trf(i,j+1,k))*di(k))
                 ! FIXME: why multiply with 1.0?
                 tfluz=zcoeff_roy(i,j,k)*(sloproy(i,j,k)*            &
                      (trf(i,j+1,ko(j,k))-trf(i,j+1,k))*             &
                      area(i, j+1) * di(k) + 1._wp * area(i, j+1) * tv_diff)

                 t2o(i,j+1,ko(j,k))=t2o(i,j+1,ko(j,k))-tfluz
                 t2o(i,j+1,k)=t2o(i,j+1,k)+tfluz

                 !  flux contributions from right/down
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_ru(i,j,k)*(tv_diff &
                      +slopruy(i,j,k)*                               &
                      (trf(i,j+1,k)-trf(i,j+1,ku(j,k)))*di(ku(j,k)))
                 ! FIXME: why multiply with 1.0?
                 tfluz=zcoeff_ruy(i,j,k)*(slopruy(i,j,k)*            &
                      (trf(i,j+1,k)-trf(i,j+1,ku(j,k)))*             &
                      area(i, j+1) * di(ku(j, k)) + 1._wp * area(i, j+1) * tv_diff)

                 t2o(i,j+1,k)=t2o(i,j+1,k)-tfluz
                 t2o(i,j+1,ku(j,k))=t2o(i,j+1,ku(j,k))+tfluz

              END IF
           END DO
        END DO
     END DO

   ELSE ! .NOT. lisopyc

!!      Calculation of Harmonic Horizontal Tracer Diffusion in Y Direction

#ifndef TRACER_OMP
!$OMP DO
#endif
     DO k=1,ke
     DO j=1,je1

           DO i=1,ie
              IF (amsue(i, j, k) .GT. 0._wp) THEN
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))*dlyvi(i,j)
                 tfluy(i,j,k)=yflux(i,j,k)*tv_diff
              ENDIF
           END DO
        END DO
     END DO

   ENDIF ! lisopyc

!        SUMMATION
!
!$OMP SINGLE
     CALL bounds_exch(1,'v',tfluy,'octdiff_trf 4')
!$OMP END SINGLE

#ifndef TRACER_OMP
! #slo# this loop cannot be parallelised in k - due to outerunroll ??
! !!CDIR OUTERUNROLL=20
!$OMP DO
#endif
     DO k=1,ke
        DO j=2,je1
           DO i=1,ie
              t2o(i,j,k)=t2o(i,j,k)+tfluy(i,j,k)-tfluy(i,j-1,k)
!              t2o(i,j+1,k)=t2o(i,j+1,k)-tfluy(i,j,k)
           END DO
        END DO
     END DO

!$OMP SINGLE
     CALL bounds_exch(1,'p',t2o,'octdiff trf 5')
!$OMP END SINGLE

#ifndef TRACER_OMP
!$OMP DO
#endif
     DO j=2,je1
        DO k=1,ke
           DO i=1,ie
              ten(i,j,k)=t2o(i,j,k)/vol_term(i,j,k)-trf(i,j,k)
           END DO
        END DO
     END DO

  ELSE ! ah00.gt.almzer

#ifndef TRACER_OMP
!$OMP WORKSHARE
#endif
    ten(:,:,:) = 0._wp
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif

  END IF ! ah00.gt.almzer

  ! VERTICAL DIFFUSION -tracer dependent

  !!The tracer value at old time step is trf

#ifndef TRACER_OMP
!$OMP DO
#endif
   DO k=1,ke
     DO j=1,je
       DO i=1,ie
         t2o(i,j,k)=trf(i,j,k)
       END DO
     END DO
   END DO

  !     IMPLICIT VERTICAL DIFFUSION

   !!This is the right hand side matrix operation, when eliminating the
   !lower diagonal of the tridiagonal matrix.

#ifndef TRACER_OMP
!$OMP DO
#endif
    DO k=2,ke
       DO j=2,je-1
      DO i=1,ie
        t2o(i,j,k)=t2o(i,j,k)-tridsy(i,j,k-1,1)*t2o(i,j,k-1)
      END DO
    END DO
  END DO

  !!This is the solution for the last ke tracer.

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO j = 2,je-1
    DO i = 1,ie
      trf(i,j,ke)=t2o(i,j,ke)/tridsy(i,j,ke,2)
    END DO
  END DO

  !!This is the solution for tracers in all other levels.
  !!System is solved from lower level ke toward the top.

#ifndef TRACER_OMP
!$OMP DO
#endif
    DO k = ke-1, 1, -1
  DO j = 2, je-1

      DO i = 1,ie
        trf(i,j,k)=(t2o(i,j,k)-tridsy(i,j,k,3)*trf(i,j,k+1))              &
                  /tridsy(i,j,k,2)
      END DO
    END DO
  END DO


!$OMP SINGLE
  CALL bounds_exch(1,'p',trf,'octdiff_trf 6')
!$OMP END SINGLE


#ifndef TRACER_OMP
!$OMP DO
#endif
    DO k = 1,ke
  DO j = 2,je1

      DO i = 1,ie
        trf(i,j,k)=trf(i,j,k)+ten(i,j,k)
      END DO
    END DO
  END DO

!$OMP SINGLE
  CALL bounds_exch(1,'p',trf,'octdiff_trf 7')
!$OMP END SINGLE

  !
  ! BIHARMONIC tracer DIFFUSION - tracer dependent
  !
  IF (aulapts.GT.almzer) THEN

#ifndef TRACER_OMP
!$OMP DO
#endif
    DO k=1,ke
      DO j=2,je1
        DO i=2,ie1
          t2o(i,j,k)=weto(i,j,k)*(                                           &
              (weto(i-1,j,k)*(trf(i-1,j,k)-trf(i,j,k))*dlxui(i-1,j)          &
              +weto(i+1,j,k)*(trf(i+1,j,k)-trf(i,j,k))*dlxui(i,j))/dlxp(i,j) &
              +(weto(i,j-1,k)*(trf(i,j-1,k)-trf(i,j,k))*dlyvi(i,j-1)         &
              +weto(i,j+1,k)*(trf(i,j+1,k)-trf(i,j,k))*dlyvi(i,j))/dlyp(i,j))
        END DO
      END DO
    END DO

!$OMP SINGLE
    CALL bounds_exch(1,'p',t2o,'octdiff_trf 8')
!$OMP END SINGLE

#ifndef TRACER_OMP
!$OMP DO
#endif
    DO k=1,ke
      zzsur0 = 0._wp
      if (k == 1) zzsur0 = 1._wp
      DO j=2,je1
        DO i=2,ie1
          trf(i,j,k)=trf(i,j,k)                                              &
               - (weto(i,j,k)/((ddpo(i,j,k)+almzer                           &
               +zzsur0*(zo(i,j)-rhoicwa*sictho(i,j)-rhosnwa*sicsno(i,j)))    &
               *area(i,j)))*dzw(k)*(                                        &
               weto(i-1,j,k)*aulx(i-1,j)*dlyu(i-1,j)*                        &
               (t2o(i-1,j,k)-t2o(i,j,k))*dlxui(i-1,j)                        &
               +weto(i+1,j,k)*aulx(i,j)*dlyu(i,j)*                           &
               (t2o(i+1,j,k)-t2o(i,j,k))*dlxui(i,j)                          &
               +weto(i,j-1,k)*auly(i,j-1)*dlxv(i,j-1)*                       &
               (t2o(i,j-1,k)-t2o(i,j,k))*dlyvi(i,j-1)                        &
               +weto(i,j+1,k)*auly(i,j)*dlxv(i,j)*                           &
               (t2o(i,j+1,k)-t2o(i,j,k))*dlyvi(i,j))
        END DO
      END DO
    END DO

!$OMP SINGLE
    CALL bounds_exch(1,'p',trf,'octdiff_trf 9')
!$OMP END SINGLE

  ENDIF ! aulapts > almzer

#ifndef TRACER_OMP
!$OMP END PARALLEL
#endif

END SUBROUTINE octdiff_trf

SUBROUTINE alloc_mem_octdiff

  ALLOCATE(vol_term(ie,je,ke))
  vol_term(:,:,:) = 0._wp

  IF (lisopyc) THEN
  !:: three-dimensional fields
  ALLOCATE(sloploy(ie,je,ke), slopluy(ie,je,ke), sloproy(ie,je,ke),       &
           slopruy(ie,je,ke), sloplox(ie,je,ke), sloplux(ie,je,ke),       &
           sloprox(ie,je,ke), sloprux(ie,je,ke),                          &
           xcoeff_lo(ie,je,ke),xcoeff_lu(ie,je,ke),xcoeff_ro(ie,je,ke),   &
           xcoeff_ru(ie,je,ke),zcoeff_loy(ie,je,ke),zcoeff_luy(ie,je,ke), &
           zcoeff_roy(ie,je,ke),zcoeff_ruy(ie,je,ke),zcoeff_lox(ie,je,ke),&
           zcoeff_lux(ie,je,ke),zcoeff_rox(ie,je,ke),zcoeff_rux(ie,je,ke),&
           ycoeff_lo(ie,je,ke),ycoeff_lu(ie,je,ke),ycoeff_ro(ie,je,ke),   &
           ycoeff_ru(ie,je,ke))

  sloploy(:,:,:) = 0._wp
  slopluy(:,:,:) = 0._wp
  sloproy(:,:,:) = 0._wp
  slopruy(:,:,:) = 0._wp
  sloplox(:,:,:) = 0._wp
  sloplux(:,:,:) = 0._wp
  sloprox(:,:,:) = 0._wp
  sloprux(:,:,:) = 0._wp
  xcoeff_lo(:,:,:) = 0._wp
  xcoeff_lu(:,:,:) = 0._wp
  xcoeff_ro(:,:,:) = 0._wp
  xcoeff_ru(:,:,:) = 0._wp
  zcoeff_loy(:,:,:) = 0._wp
  zcoeff_luy(:,:,:) = 0._wp
  zcoeff_roy(:,:,:) = 0._wp
  zcoeff_ruy(:,:,:) = 0._wp
  zcoeff_lox(:,:,:) = 0._wp
  zcoeff_lux(:,:,:) = 0._wp
  zcoeff_rox(:,:,:) = 0._wp
  zcoeff_rux(:,:,:) = 0._wp
  ycoeff_lo(:,:,:) = 0._wp
  ycoeff_lu(:,:,:) = 0._wp
  ycoeff_ro(:,:,:) = 0._wp
  ycoeff_ru(:,:,:) = 0._wp

ELSE ! .NOT. lisopyc

  ALLOCATE(xflux(ie,je,ke),yflux(ie,je,ke))

  xflux(:,:,:) = 0._wp
  yflux(:,:,:) = 0._wp

ENDIF ! lisopyc

  WRITE(0,*)aulapts ,almzer
  IF (aulapts.GT.almzer) THEN
    ALLOCATE(aulx(ie,je),auly(ie,je))
    aulx(:,:) = 0._wp
    auly(:,:) = 0._wp
  ENDIF

  ALLOCATE(zbott(je,ke),zsurf(je,ke),       &
       ko(je,ke),ku(je,ke))

  ALLOCATE( tridsy(ie,je,ke,3))
  tridsy(:,:,:,:) = 0._wp



      END SUBROUTINE alloc_mem_octdiff

end module mo_diffusion
