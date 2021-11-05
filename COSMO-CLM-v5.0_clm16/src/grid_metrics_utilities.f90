!+ Module for metric utilities
!------------------------------------------------------------------------------

MODULE grid_metrics_utilities

!------------------------------------------------------------------------------
!
! Description:
!   This module generates all the data describing the grid and the 
!   coordinate transformation.
!   The module contains only subroutines which don't need any 
!   other module (only exception: data_parameters).
!   These routines are called by 'init_grid_metrics' (possible
!   boundary exchange, periodic BC's or 2D-simulations
!   are also treated there)
!
! Method:
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8236 1493
!  email:  michael.baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V4_23        2012/05/10 Michael Baldauf
!  Initial version: only some subroutines are included now, not to 
!  shuffle too much code around in one version change
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,  ONLY:  &
  ireals,    & ! KIND-type parameter for real variables
  iintegers    ! KIND-type parameter for standard integer variables

USE data_runcontrol,  ONLY:  &
  l2tls,            & ! forecast with 2-TL integration scheme
  l_dzeta_d_needed, &
  l2dim,            &
  lperi_x,          & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                              ! or with Davies conditions (.FALSE.)
  lperi_y             ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                      ! or with Davies conditions (.FALSE.)

USE data_modelconfig, ONLY:  &
  ie,           & ! number of grid points in zonal direction
  je,           & ! number of grid points in meridional direction
  ke,           & ! number of grid points in vertical direction
  ke1,          & ! KE+1
  jstart,       & ! start index for the forecast of w, t, qv, qc and pp
  jend,         & ! end index for the forecast of w, t, qv, qc and pp
  jstartu,      & ! start index for the forecast of u
  jendu,        & ! end index for the forecast of u
  jstartv,      & ! start index for the forecast of v
  jendv,        & ! end index for the forecast of v
  jstartpar,    & ! start index for computations in the parallel program
  jendpar,      & ! end index for computations in the parallel program
  dlon,         & ! grid point distance in zonal direction (in degrees)
  dlat,         & ! grid point distance in meridional direction (in degrees)
  eddlon,       &
  eddlat,       &
  dt

USE data_fields,      ONLY:  &
  hhl               ! geometical height of half model levels        ( m )

USE data_parallel,    ONLY:  &
  my_cart_id,       & ! rank of this subdomain in the cartesian communicator
  my_cart_neigh,    & ! neighbors
  num_compute,      & ! number of compute PEs
  nboundlines,      & ! number of boundary lines of the domain for which
                      ! no forecast is computed = overlapping boundary
                      ! lines of the subdomains
  icomm_cart,       & ! communicator that belongs to the cartesian grid
  imp_reals,        & ! determines the correct REAL type used in the model
                      ! for MPI
  ncomm_type,       & ! type of communication
  sendbuf,          & ! sending buffer for boundary exchange:
                      ! 1-4 are used for sending, 5-8 are used for receiving
  isendbuflen         ! length of one column of sendbuf

USE environment,      ONLY:  exchg_boundaries

!==============================================================================

IMPLICIT NONE

!==============================================================================

  REAL  (KIND=ireals), ALLOCATABLE ::           &
    sqrtg_r_s (:,:,:),   & ! 1 / square root of G at scalar points       (1/m)
    sqrtg_r_u (:,:,:),   & ! 1 / square root of G at u points            (1/m)
    sqrtg_r_v (:,:,:),   & ! 1 / square root of G at v points            (1/m)
    sqrtg_r_w (:,:,:),   & ! 1 / square root of G at w points            (1/m)
    dzeta_dlam(:,:,:),   & ! d zeta / d lambda (for constant phi,    z)
                           ! at the scalar position                      ( 1 )
    dzeta_dphi(:,:,:)      ! d zeta / d phi    (for constant lambda, z)
                           ! at the scalar position                      ( 1 )

  REAL  (KIND=ireals), ALLOCATABLE ::  &
    wgtfac   (:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfac_u (:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfac_v (:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfacq  (:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfacq_u(:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfacq_v(:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgt_one_sided(:,:,:)   ! weighting factor for vertical interpolation ( 1 )

!==============================================================================

CONTAINS

!==============================================================================
!==============================================================================

SUBROUTINE init_grid_metrics (ierrorstat)

!------------------------------------------------------------------------------
!
! Description:
!   call of the routines which compute the necessary metric variables
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT) :: ierrorstat

  INTEGER :: j, k
  INTEGER :: istat

  INTEGER (KIND=iintegers) :: kzdims(24)
  INTEGER (KIND=iintegers) :: izerror
  CHARACTER(100)           :: yzerrmsg

!------------------------------------------------------------------------------

  ierrorstat = 0

  istat = 0

  IF ( my_cart_id == 0 ) WRITE(*,*) "Subr.[init_grid_metrics] ..."

! Allocation of fields
! Note: there is no corresponding deallocate at the moment
  ALLOCATE ( sqrtg_r_s(ie,je,ke), STAT=istat ); sqrtg_r_s=0.0_ireals
  ALLOCATE ( sqrtg_r_u(ie,je,ke), STAT=istat ); sqrtg_r_u=0.0_ireals
  ALLOCATE ( sqrtg_r_v(ie,je,ke), STAT=istat ); sqrtg_r_v=0.0_ireals
  ALLOCATE ( sqrtg_r_w(ie,je,ke1),STAT=istat ); sqrtg_r_w=0.0_ireals
  IF ( l_dzeta_d_needed ) THEN
    ALLOCATE ( dzeta_dlam(ie,je,ke), STAT=istat ); dzeta_dlam=0.0_ireals
    ALLOCATE ( dzeta_dphi(ie,je,ke), STAT=istat ); dzeta_dphi=0.0_ireals
  END IF

  ALLOCATE ( wgtfac    (ie,je,ke1), STAT=istat ); wgtfac    = 0.0_ireals
  ALLOCATE ( wgtfac_u  (ie,je,ke1), STAT=istat ); wgtfac_u  = 0.0_ireals
  ALLOCATE ( wgtfac_v  (ie,je,ke1), STAT=istat ); wgtfac_v  = 0.0_ireals
  ALLOCATE ( wgtfacq   (ie,je,3),   STAT=istat ); wgtfacq   = 0.0_ireals
  ALLOCATE ( wgtfacq_u (ie,je,3),   STAT=istat ); wgtfacq_u = 0.0_ireals
  ALLOCATE ( wgtfacq_v (ie,je,3),   STAT=istat ); wgtfacq_v = 0.0_ireals

!------------------------------------------------------------------------------

  IF ( l2tls ) THEN
    CALL weighting_factors_full2half( hhl,                                &
              wgtfac, wgtfac_u, wgtfac_v, wgtfacq, wgtfacq_u, wgtfacq_v,  &
              ie, je, ke )
  END IF


  IF ( l2dim ) THEN
    DO j=1, nboundlines
      wgtfac(:,jstart-j,:) = wgtfac(:,jstart,:)
      wgtfac(:,jend  +j,:) = wgtfac(:,jstart,:)

      wgtfac_u(:,jstartu-j,:) = wgtfac_u(:,jstartu,:)
      wgtfac_u(:,jendu  +j,:) = wgtfac_u(:,jstartu,:)

      wgtfac_v(:,jstartv-j,:) = wgtfac_v(:,jstartv,:)
      wgtfac_v(:,jendv  +j,:) = wgtfac_v(:,jstartv,:)

      wgtfacq(:,jstart-j,:) = wgtfacq(:,jstart,:)
      wgtfacq(:,jend  +j,:) = wgtfacq(:,jstart,:)

      wgtfacq_u(:,jstartu-j,:) = wgtfacq_u(:,jstartu,:)
      wgtfacq_u(:,jendu  +j,:) = wgtfacq_u(:,jstartu,:)

      wgtfacq_v(:,jstartv-j,:) = wgtfacq_v(:,jstartv,:)
      wgtfacq_v(:,jendv  +j,:) = wgtfacq_v(:,jstartv,:)

    END DO
    wgtfac_v (:,jendv,:) = wgtfac_v (:,jstartv,:)   ! notwendig ???????
    wgtfacq_v(:,jendv,:) = wgtfacq_v(:,jstartv,:)   ! notwendig ???????
  END IF


#ifndef NOMPI
  kzdims(1:24) =                           &
       (/ ke+1, ke+1, ke+1,   3,   3,      &
             3,    0,    0,   0,   0,      &
             0,    0,    0,   0,   0,      &
             0,    0,    0,   0,   0,      &
             0,    0,    0,   0  /)

  ! nicht optimierte Version (icase=0 <--> ldatatypes=.FALSE.)
  CALL exchg_boundaries                                                &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,  &
        kzdims, jstartpar, jendpar, nboundlines, nboundlines, my_cart_neigh,     &
        lperi_x, lperi_y, l2dim, &
        0, .FALSE., ncomm_type, izerror, yzerrmsg,       &
        wgtfac   (:,:,:),     &
        wgtfac_u (:,:,:),     &
        wgtfac_v (:,:,:),     &
        wgtfacq  (:,:,:),     &
        wgtfacq_u(:,:,:),     &
        wgtfacq_v(:,:,:) )

  IF ( izerror /= 0_iintegers ) THEN
    ierrorstat = 20005
  END IF
#endif


  ! --- vertical weighting factors (2) ------------------------------

  ALLOCATE( wgt_one_sided (ie,je,3), STAT=istat )

  CALL weighting_factors_one_sided( hhl, wgt_one_sided, ie, je, ke )


  ! --- metric coeeficients of the terrain following coordinate -----

  ! initialize sqrt(G) for all time level schemes
  CALL calc_sqrtg_r (hhl, sqrtg_r_s, sqrtg_r_u, sqrtg_r_v, sqrtg_r_w, &
                     ie, je, ke )


  IF ( l2dim ) THEN
    DO j=1, nboundlines
      sqrtg_r_s(:,jstart-j,:) = sqrtg_r_s(:,jstart,:)
      sqrtg_r_s(:,jend  +j,:) = sqrtg_r_s(:,jstart,:)

      sqrtg_r_w(:,jstart-j,:) = sqrtg_r_w(:,jstart,:)
      sqrtg_r_w(:,jend  +j,:) = sqrtg_r_w(:,jstart,:)

      sqrtg_r_u(:,jstartu-j,:) = sqrtg_r_u(:,jstartu,:)
      sqrtg_r_u(:,jendu  +j,:) = sqrtg_r_u(:,jstartu,:)

      sqrtg_r_v(:,jstartv-j,:) = sqrtg_r_v(:,jstartv,:)
      sqrtg_r_v(:,jendv  +j,:) = sqrtg_r_v(:,jstartv,:)
    END DO
    sqrtg_r_v(:,jendv,:) = sqrtg_r_v(:,jstartv,:)   ! notwendig ???????
  END IF

#ifndef NOMPI
  kzdims(1:24) =                         &
       (/ ke , ke ,  ke, ke+1,   0,      &
            0,   0,   0,    0,   0,      &
            0,   0,   0,    0,   0,      &
            0,   0,   0,    0,   0,      &
            0,   0,   0,    0  /)

  ! nicht optimierte Version (icase=0 <--> ldatatypes=.FALSE.)
  CALL exchg_boundaries                                                &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,  &
        kzdims, jstartpar, jendpar, nboundlines, nboundlines, my_cart_neigh,     &
        lperi_x, lperi_y, l2dim, &
        0, .FALSE., ncomm_type, izerror, yzerrmsg,       &
        sqrtg_r_s(:,:,:),     &
        sqrtg_r_u(:,:,:),     &
        sqrtg_r_v(:,:,:),     &
        sqrtg_r_w(:,:,:) )

  IF ( izerror /= 0_iintegers ) THEN
    ierrorstat = 20005
  END IF
#endif


  IF ( l_dzeta_d_needed ) THEN

    CALL metric_coeffs ( hhl, sqrtg_r_s, eddlon, eddlat,       &
      dzeta_dlam, dzeta_dphi, ie, je, ke )


    IF ( l2dim ) THEN
      DO j=1, nboundlines
        dzeta_dlam(:,jstart-j,:) = dzeta_dlam(:,jstart,:)
        dzeta_dlam(:,jend  +j,:) = dzeta_dlam(:,jstart,:)

        dzeta_dphi(:,jstart-j,:) = dzeta_dphi(:,jstart,:)
        dzeta_dphi(:,jend  +j,:) = dzeta_dphi(:,jstart,:)

      END DO
    END IF


#ifndef NOMPI
    kzdims(1:24) =                        &
         (/ ke , ke ,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0  /)

    ! nicht optimierte Version (icase=0 <--> ldatatypes=.FALSE.)
    CALL exchg_boundaries                                                &
         (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,  &
          kzdims, jstartpar, jendpar, nboundlines, nboundlines, my_cart_neigh,     &
          lperi_x, lperi_y, l2dim, &
          0, .FALSE., ncomm_type, izerror, yzerrmsg,       &
          dzeta_dlam(:,:,:),     &
          dzeta_dphi(:,:,:)    )

    IF ( izerror /= 0_iintegers ) THEN
      ierrorstat = 20005
    END IF
#endif

  END IF

END SUBROUTINE init_grid_metrics

!==============================================================================
!==============================================================================

SUBROUTINE calc_sqrtg_r (hhl, sqrtg_r_s, sqrtg_r_u, sqrtg_r_v, sqrtg_r_w, &
                         ie, je, ke)

!------------------------------------------------------------------------------
!
! Description:
!   calculate the reciprocal square root of G (Dim.: [1/m]) at several grid positions
!   e.g. as a preliminary for Subr. 'integral_3D'.
!   This routine is only needed, if the 3-timelevel scheme is used
!   (because sqrtg_r_s, ... will not be calculated there)
!   A good place for calling this routine is e.g.
!   in Subr. 'init_dynamics' (file oranize_dynamics.f90)
!
! Method:
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN)  ::         &
    ie, je, ke              ! dimensions

  REAL (KIND = ireals),     INTENT(IN)  ::         &
    hhl       (ie,je,ke+1)  ! height of half levels

  REAL (KIND = ireals),     INTENT(OUT) ::         &
    sqrtg_r_s (ie,je,ke), & ! for scalar grid points
    sqrtg_r_u (ie,je,ke), & ! for u-grid points
    sqrtg_r_v (ie,je,ke), & ! for v-grid points
    sqrtg_r_w (ie,je,ke+1)  ! for w-grid points

  INTEGER (KIND=iintegers) :: i, j, k

!------------------------------------------------------------------------------

  ! for scalar , u- and v- points
  DO  k = 1, ke

    DO j = 1, je
      DO i = 1, ie
        sqrtg_r_s(i,j,k) = 1.0_ireals / ( hhl(i,j,k)-hhl(i,j,k+1) )
      ENDDO
    ENDDO
    ! (--> no exchange necessary)

    DO j = 1, je
      DO i = 1, ie-1
        sqrtg_r_u(i,j,k) = 2.0_ireals / ( hhl(i,j,k  )+hhl(i+1,j,k  )    &
                                         -hhl(i,j,k+1)-hhl(i+1,j,k+1) )
      ENDDO
      ! approx. for lateral boundaries:
      sqrtg_r_u(ie,j,k) = 1.0_ireals / ( hhl(ie,j,k  )-hhl(ie,j,k+1) ) 
    ENDDO

    DO j = 1, je-1
      DO i = 1, ie
        sqrtg_r_v(i,j,k) = 2.0_ireals / ( hhl(i,j,k  )+hhl(i,j+1,k  )    &
                                         -hhl(i,j,k+1)-hhl(i,j+1,k+1) )
      ENDDO
    ENDDO
    ! approx. for lateral boundaries:
    DO i = 1, ie
      sqrtg_r_v(i,je,k) = 1.0_ireals / ( hhl(i,je,k)-hhl(i,je,k+1) )
    ENDDO

  ENDDO

  ! and for the w-points
  DO  k = 2, ke
    DO j = 1, je
      DO i = 1, ie
        ! version A:
         sqrtg_r_w(i,j,k) = 2.0_ireals / ( hhl(i,j,k-1)-hhl(i,j,k+1) )

        ! version B:
        !sqrtg_r_w(i,j,k) =    wgtfac(i,j,k)   / ( hhl(i,j,k)   - hhl(i,j,k+1) )  &
        !      + (1.0_ireals - wgtfac(i,j,k) ) / ( hhl(i,j,k-1) - hhl(i,j,k)   )
      ENDDO
    ENDDO
  ENDDO

  ! one-sided differences in k=1 and k=ke+1
  DO j = 1, je
    DO i = 1, ie
      sqrtg_r_w(i,j,1)    = 1.0_ireals / ( hhl(i,j,1 )-hhl(i,j,2   ) )
      sqrtg_r_w(i,j,ke+1) = 1.0_ireals / ( hhl(i,j,ke)-hhl(i,j,ke+1) )
    ENDDO
  ENDDO

END SUBROUTINE calc_sqrtg_r

!==============================================================================
!==============================================================================

SUBROUTINE metric_coeffs( hhl, sqrtg_r_s, eddlon, eddlat,       &
                          dzeta_dlam, dzeta_dphi, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   calculate the metric coefficients (see Doms, Schaettler (2002) section 3):
!     d zeta / d lambda (at constant phi,    z) = - dzeta/dz * dz/dlambda
!     d zeta / d phi    (at constant lambda, z) = - dzeta/dz * dz/dphi
!   at the scalar grid positions.
!   (remarks:
!     dzeta/dz = 1/( dz/dzeta ) = -1/sqrt(G) (see subroutine calc_sqrtg_r)
!     dz/dzeta = -sqrt(G)  )
!   This subroutine is closely related to calc_sqrtg_r and should be
!   called after that.
!
! Method:
!   centered differences using hhl for calculating dz/d...
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN)  ::         &
    ie, je, ke              ! dimensions

  REAL (KIND = ireals),     INTENT(IN)  ::         &
    hhl       (ie,je,ke+1)  ! height of half levels

  REAL (KIND = ireals),     INTENT(IN) ::          &
    sqrtg_r_s (ie,je,ke)    ! = 1 / sqrt(G)  (at scalar grid point position)

  REAL (KIND = ireals),     INTENT(IN) ::          &
    eddlon,   & ! 1 / dlon, dlon in degrees (!)    &
    eddlat      ! 1 / dlat, dlat in degrees (!)

  REAL (KIND=ireals), DIMENSION(1:ie,1:je,1:ke), INTENT(OUT) ::   &
    dzeta_dlam, dzeta_dphi

  INTEGER (KIND=iintegers) :: i, j, k

!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = 1, je
      DO  i = 2, ie-1

        dzeta_dlam(i,j,k) = sqrtg_r_s(i,j,k) *           &
          ( ( hhl(i+1,j,k  ) - hhl(i-1,j,k  ) )      &
          + ( hhl(i+1,j,k+1) - hhl(i-1,j,k+1) ) )    &
          * 0.25_ireals * eddlon

      END DO
    END DO
  END DO
  ! approximation for lateral boundaries:
  dzeta_dlam(1, :,:) = 0.0_ireals
  dzeta_dlam(ie,:,:) = 0.0_ireals


  DO k = 1, ke
    DO j = 2, je-1
      DO  i = 1, ie

        dzeta_dphi(i,j,k) = sqrtg_r_s(i,j,k) *           &
          ( ( hhl(i,j+1,k  ) - hhl(i,j-1,k  ) )      &
          + ( hhl(i,j+1,k+1) - hhl(i,j-1,k+1) ) )    &
          * 0.25_ireals * eddlat

      END DO
    END DO
  END DO
  ! approximation for lateral boundaries:
  dzeta_dphi(:,1 ,:) = 0.0_ireals
  dzeta_dphi(:,je,:) = 0.0_ireals

END SUBROUTINE metric_coeffs

!==============================================================================
!==============================================================================

SUBROUTINE weighting_factors_full2half( hhl,                        &
        wgtfac, wgtfac_u, wgtfac_v, wgtfacq, wgtfacq_u, wgtfacq_v,  &
        ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
! Weighting factors for linear interpolations from full levels 
! to half levels (for metric terms) are calculated. All fields are 
! set on the full (sub)-domain. In parallel programs the neighboring
! values have to be exchanged after calling this routine.
!
! Input:
!   hhl
!
! Output:
!   wgtfac, wgtfac_u, wgtfac_v, wgtfacq, wgtfacq_u, wgtfacq_v
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN) :: ie, je, ke

  REAL (KIND=ireals), INTENT(IN)  :: hhl       (ie,je,ke+1)

  REAL (KIND=ireals), INTENT(OUT) :: wgtfac    (ie,je,ke+1)
  REAL (KIND=ireals), INTENT(OUT) :: wgtfac_u  (ie,je,ke+1)
  REAL (KIND=ireals), INTENT(OUT) :: wgtfac_v  (ie,je,ke+1)
  REAL (KIND=ireals), INTENT(OUT) :: wgtfacq   (ie,je,   3)
  REAL (KIND=ireals), INTENT(OUT) :: wgtfacq_u (ie,je,   3)
  REAL (KIND=ireals), INTENT(OUT) :: wgtfacq_v (ie,je,   3)

  INTEGER            :: i, j, k
  INTEGER            :: istat
  REAL (KIND=ireals) :: z1, z2, z3

!------------------------------------------------------------------------------

  DO  k = 2, ke
    DO  j = 1, je
      !CDIR ON_ADB(hhl)
      DO  i = 1, ie
        wgtfac(i,j,k) = ( hhl(i,j,k-1) - hhl(i,j,k  ) )  &
          &           / ( hhl(i,j,k-1) - hhl(i,j,k+1) )
      ENDDO
    ENDDO
  ENDDO

  DO  j = 1, je
    !CDIR ON_ADB(hhl)
    !CDIR ON_ADB(wgtfac)
    DO  i = 1, ie
      wgtfac(i,j,1)    = ( hhl(i,j,2) - hhl(i,j,1)   )     &
        &              / ( hhl(i,j,3) - hhl(i,j,1)   )
      wgtfac(i,j,ke+1) = ( hhl(i,j,ke  ) - hhl(i,j,ke+1) )  &
        &              / ( hhl(i,j,ke-1) - hhl(i,j,ke+1) )
    ENDDO
  ENDDO

  DO  k = 1, ke+1
    DO  j = 1, je
      !CDIR ON_ADB(wgtfac)
      DO  i = 1, ie-1
        wgtfac_u(i,j,k) = 0.5_ireals*(wgtfac(i+1,j,k)+wgtfac(i,j,k))
      ENDDO
    ENDDO
  ENDDO
  wgtfac_u(ie,:,:) = 0.5_ireals  ! approximation for lateral boundaries

  DO  k = 1, ke+1
    DO  j = 1, je-1
      !CDIR ON_ADB(wgtfac)
      DO  i = 1, ie
        wgtfac_v(i,j,k) = 0.5_ireals*(wgtfac(i,j+1,k)+wgtfac(i,j,k))
      ENDDO
    ENDDO
  ENDDO
  wgtfac_v(:,je,:) = 0.5_ireals  ! approximation for lateral boundaries

  ! Factors for quadratic extrapolation of surface pressure 
  ! (used if ldyn_bbc = .FALSE.)
  DO  j = 1, je
    !CDIR ON_ADB(wgtfacq)
    DO  i = 1, ie
      z1 = 0.5_ireals * ( hhl(i,j,ke  ) - hhl(i,j,ke+1 ) )
      z2 = 0.5_ireals * ( hhl(i,j,ke  ) + hhl(i,j,ke-1) ) - hhl(i,j,ke+1)
      z3 = 0.5_ireals * ( hhl(i,j,ke-1) + hhl(i,j,ke-2) ) - hhl(i,j,ke+1)
      wgtfacq(i,j,3) = z1*z2/(z2-z3)/(z1-z3)
      wgtfacq(i,j,2) = (z1-wgtfacq(i,j,3)*(z1-z3))/(z1-z2)
      wgtfacq(i,j,1) = 1._ireals - (wgtfacq(i,j,2) + wgtfacq(i,j,3))
    ENDDO
  ENDDO

  DO k = 1, 3
    DO  j = 1, je
      !CDIR ON_ADB(wgtfacq)
      DO  i = 1, ie-1
        wgtfacq_u(i,j,k) = 0.5_ireals*(wgtfacq(i+1,j,k)+wgtfacq(i,j,k))
      ENDDO
    ENDDO
  ENDDO
  wgtfacq_u(ie,:,:) = 0.5_ireals   ! approximation for lateral boundaries

  DO k = 1, 3
    DO  j = 1, je-1
      !CDIR ON_ADB(wgtfacq)
      DO  i = 1, ie
        wgtfacq_v(i,j,k) = 0.5_ireals*(wgtfacq(i,j+1,k)+wgtfacq(i,j,k))
      ENDDO
    ENDDO
  ENDDO
  wgtfacq_v(:,je,:) = 0.5_ireals   ! approximation for lateral boundaries

END SUBROUTINE weighting_factors_full2half

!==============================================================================
!==============================================================================

SUBROUTINE weighting_factors_one_sided ( hhl, wgt_one_sided, &
                 ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
! calculate weighting factors for vertical 2nd order one sided 
! finite difference formula
! (valid for scalar variables like pressure p)
! dp/dz = wgt_one_sided(1) * p(ke)
!       + wgt_one_sided(2) * p(ke-1)
!       + wgt_one_sided(3) * p(ke-2)
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN)  :: ie, je, ke
  REAL    (KIND=ireals),    INTENT(IN)  :: hhl          (ie,je,ke+1)

  REAL    (KIND=ireals),    INTENT(OUT) :: wgt_one_sided(ie,je,   3)

  REAL    (KIND=ireals) :: delta1, delta2

  INTEGER :: i, j

!------------------------------------------------------------------------------

  DO j=1, je
    DO i=1, ie

      ! distance between main levels at ke-2 and ke:
      delta2 = 0.5_ireals * ( ( hhl(i,j,ke-2) + hhl(i,j,ke-1) )   &
        &                   - ( hhl(i,j,ke  ) + hhl(i,j,ke+1) ) )

      ! distance between main levels at ke-1 and ke:
      delta1 = 0.5_ireals * ( ( hhl(i,j,ke-1) + hhl(i,j,ke)   )   &
        &                   - ( hhl(i,j,ke  ) + hhl(i,j,ke+1) ) )

      ! quadratic weights:
      wgt_one_sided(i,j,3) = delta1 / delta2 / ( delta1 - delta2)
      wgt_one_sided(i,j,2) = ( 1.0_ireals - wgt_one_sided(i,j,3) * delta2) / delta1
      wgt_one_sided(i,j,1) = -( wgt_one_sided(i,j,2) + wgt_one_sided(i,j,3) )

      ! linear weights:
      !wgt_one_sided(i,j,3) =  0.0_ireals
      !wgt_one_sided(i,j,2) =  1.0_ireals / delta1
      !wgt_one_sided(i,j,1) = -1.0_ireals / delta1

    END DO
  END DO

  ! i=4
  ! j=4
  ! WRITE(*,'(A,3I4,3E14.5)') "wgt ", my_cart_id, i, j,          &
  !   & wgt_one_sided(i, j,1), wgt_one_sided(i, j,2), wgt_one_sided(i, j,3)

END SUBROUTINE weighting_factors_one_sided

!==============================================================================

END MODULE grid_metrics_utilities
