SUBROUTINE si1_extended_const

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc => local_decomposition
  USE mo_semi_impl,     ONLY: betazq
  USE mo_control,       ONLY: nlev
  USE mo_time_control,  ONLY: time_step_len
  USE mo_scan_buffer,   ONLY: u0
  USE mo_geoloc,        ONLY: rz1u0
  USE mo_gaussgrid,     ONLY: gl_racst
  USE mo_transpose,     ONLY: reorder

  IMPLICIT NONE

  INTEGER :: jlev
  INTEGER :: jlat, lnlat, jglat
  REAL(dp):: zdt
  REAL(dp):: zrz1u0(ldc% nglon, nlev, ldc% nglat)

  zdt = 0.5_dp*time_step_len

  lnlat = ldc%nglat

  DO jlat = 1, lnlat
    jglat = ldc% glat(jlat) ! global latitude index N -> S

    DO jlev = 1, nlev
      zrz1u0(:,jlev,jlat)=(betazq*zdt*gl_racst(jglat))*u0(jlev,jlat)
    END DO
  END DO

  CALL reorder (rz1u0, zrz1u0)

END SUBROUTINE si1_extended_const

!OCL NOALIAS

SUBROUTINE si1(krow)

  ! Description:
  !
  ! 1st part of the semi-implicit scheme (done in grid point space). 
  !
  ! Method:
  !
  ! This subroutine computes the contribution in
  ! grid points to the semi-implicit scheme.
  !
  ! *si1* is called from *gpc*.
  !
  ! Reference:
  ! 1-appendix *b1:organisation of the spectral model.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! I. Kirchner, MPI, November 2000, date/time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! I. Kirchner, MPI, August 2002, bugfix tendency diagnostics
  ! 
  ! for more details see file AUTHORS
  !
  USE mo_kind,          ONLY: dp
#ifndef MESSY
  USE mo_scan_buffer,   ONLY: vol, vom,  qte, rh, xlte, xite, tte, xtte, &
                              alps, alpste, d, t, vo, dm
  USE mo_memory_g1a,    ONLY: vom1, qm1, xlm1, xim1, xtm1, dm1, tm1,   &
                              alpsm1
#else
  USE mo_scan_buffer,   ONLY: vol, vom,  qte, rh, xlte, xite, tte, xtte, &
                              alps, alpste, d, t, vo, dm,                &
                              xtte_a, & ! mz_pj_20030930
                              xtte_c ! op_sb_20191009
  USE mo_memory_g1a,    ONLY: vom1, qm1, xlm1, xim1, xtm1, dm1, tm1,   &
                              alpsm1, xtm1_a, &
                              xtm1_c   ! op_sb_20191009
#endif
  USE mo_geoloc,        ONLY: rz1u0, racst_2d, sqcst_2d
  USE mo_control,       ONLY: ltdiag, nlev
  USE mo_time_control,  ONLY: l_putdata, time_step_len
  USE mo_constants,     ONLY: a
#ifndef MESSY
  USE mo_tracer,        ONLY: trlist
#else
  ! mz_rs_20040329+
  USE messy_main_tracer_mem_bi, ONLY: ntrac_lg, ti_lg, ntrac_gp, ti_gp, ON &
                                    , I_integrate                          &
                                    , ti_cl,ntrac_cl,xt_c ! op_sb_20191009
  ! mz_rs_20040329-
#endif
  USE mo_diag_tendency, ONLY: pdiga, pdsga, pdigs, pdsgs, dio_index
  USE mo_decomposition, ONLY: ldc=>local_decomposition

  IMPLICIT NONE

  INTEGER :: krow

  !  Local scalars: 
  REAL(dp) :: z1u0, z2, z3, zdt, ztwodt
  INTEGER :: jlev, jl, jrow, jt
#ifndef MESSY
  INTEGER :: nglon, nproma, nbdim
#else
  INTEGER :: nproma, nbdim
#endif
  LOGICAL :: loperm

  !  Local arrays: 
  REAL(dp) :: zscale(ldc% nproma)
  REAL(dp), TARGET  :: zd(ldc% nproma, nlev)
  REAL(dp)          :: zp(ldc% nproma), zt(ldc% nproma, nlev)
  REAL(dp), POINTER :: zr(:,:)

  !  External subroutines 

  EXTERNAL conteq, pgrad


  !  Executable statements 

  !-- 1. Locate and allocate storage

  jrow  = krow        ! local  continuous index
#ifndef MESSY
  nglon = ldc% nglon  ! local  number of latitudes
#endif

  nbdim = ldc% nproma

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

  !   zp(:)   = 0.0_dp
  !   zt(:,:) = 0.0_dp
  !   zd(:,:) = 0.0_dp

  !-- 1.2 Equivalence arrays *zd* and *zr*

  zr => zd

  !-- 2. Skip over *si1* during initialisation iterations
  !      or prepare some temporary constants.

  !-- 2.1 Compute temporary constants

  zdt    = 0.5_dp*time_step_len
  ztwodt =        time_step_len

  z2 = ztwodt/a

  IF (ltdiag) THEN
    IF (l_putdata(dio_index)) THEN
      ! rescale vom/vol terms calculated in physc
      zscale(1:nproma) = sqcst_2d(1:nproma,jrow)
      pdiga(:,:,3: 5,jrow) = pdiga(:,:,3: 5,jrow)*SPREAD(SPREAD(zscale(:),2,nlev),3,3)
      pdiga(:,:,8:10,jrow) = pdiga(:,:,8:10,jrow)*SPREAD(SPREAD(zscale(:),2,nlev),3,3)
      ! correction of pdiga only at accumulation time steps
      ! rescaling with 
      !   pdiga(*,1...10,*) with 2*dt/a
      !   pdiga(*,11..23,*) with 2*dt
      pdiga(:,:, 1:10,jrow) = pdiga(:,:, 1:10,jrow)*z2
      pdiga(:,:,11:20,jrow) = pdiga(:,:,11:20,jrow)*ztwodt
      pdiga(:,:,24:25,jrow) = pdiga(:,:,24:25,jrow)*ztwodt
#ifndef MESSY
      pdsga(:,      1,jrow) = pdsga(:,      1,jrow)*ztwodt
#else
      ! mz_jb_20040829+
      pdsga(:,1,    1,jrow) = pdsga(:,1,    1,jrow)*ztwodt
      ! mz_jb_20040829-
#endif
    ENDIF
  END IF

  !-- 3. Semi implicit computations

  !-- 3.1 Vorticity and humidity equations
  !       and term dmu of divergence equation.
  DO jlev = 1, nlev

     ! get explicit part of time schema, (2*dt) included
     ! explicit part of vorticity
     IF (ltdiag) THEN
!DIR$ CONCURRENT
       DO jl = 1, nproma
         z1u0 = rz1u0(jl,jlev,jrow)
         pdigs(jl,jlev,1,jrow) =                                       &
              -z1u0*(vom1(jl,jlev,jrow)-2._dp*vo(jl,jlev,jrow))
       END DO
     END IF

!DIR$ IVDEP
     DO jl = 1, nproma

        z1u0 = rz1u0(jl,jlev,jrow)
        z3   = ztwodt*racst_2d(jl,jrow)

        dm(jl,jlev,jrow)  =  z2 * vol (jl,jlev,jrow)

        vol(jl,jlev,jrow) =  z3 * vol (jl,jlev,jrow)                  &
                         -z1u0*(vom1(jl,jlev,jrow) - 2._dp*vo(jl,jlev,jrow))
        vom(jl,jlev,jrow) = -z2  * vom (jl,jlev,jrow)

        qm1(jl,jlev,jrow)  = qm1 (jl,jlev,jrow) + ztwodt*qte (jl,jlev,jrow)
        xlm1(jl,jlev,jrow) = xlm1(jl,jlev,jrow) + ztwodt*xlte(jl,jlev,jrow)
        xim1(jl,jlev,jrow) = xim1(jl,jlev,jrow) + ztwodt*xite(jl,jlev,jrow)
     END DO
  END DO

#ifndef MESSY
  DO jt = 1, trlist% ntrac
     if (trlist% ti(jt)% nint /= 1) CYCLE
     DO jlev = 1, nlev
        DO jl = 1, nproma
           xtm1(jl,jlev,jt,jrow) = xtm1(jl,jlev,jt,jrow) +             &
                                                ztwodt*xtte(jl,jlev,jt,jrow)
        END DO
     END DO
  END DO
#else
  ! mz_rs_20040329+
  DO jt = 1, ntrac_gp
    IF (ti_gp(jt)%tp%meta%cask_i(I_integrate) == ON) THEN
      DO jlev = 1, nlev
        DO jl = 1, nproma
          xtm1(jl,jlev,jt,jrow) = &
            xtm1(jl,jlev,jt,jrow) + ztwodt*xtte(jl,jlev,jt,jrow)
        END DO
      END DO
    END IF
  END DO
  ! mz_rs_20040329-

  ! mz_pj_20030930+
  ! JROW = ldc%ngpblks -> TRANSFER AT END OF LOOP
  IF (jrow == ldc%ngpblks) THEN
    DO jt = 1, ntrac_lg
      IF (ti_lg(jt)%tp%meta%cask_i(I_integrate) == ON) THEN
        xtm1_a(:,:,jt,1) = xtm1_a(:,:,jt,1) + ztwodt*xtte_a(:,:,jt,1)
      END IF
    END DO
  END IF
  ! mz_pj_20030930-
  ! op_sb_20191007+
  IF (jrow == ldc%ngpblks) THEN
    DO jt = 1, ntrac_cl
      IF (ti_cl(jt)%tp%meta%cask_i(I_integrate) == ON) THEN
!!$      xtm1_c(:,:,jt,1) = xt_c(:,:,jt,1) + time_step_len*xtte_c(:,:,jt,1)
         xtm1_c(:,:,jt,1) = xt_c(:,:,jt,1) + zdt*xtte_c(:,:,jt,1)
      END IF
    END DO
   END IF
   ! op_sb_20191007-
#endif

  !-- 3.2 Compute implicit contribution of divergence
  !       to temperature and surface equations.

  DO jlev = 1, nlev
     DO jl = 1, nproma
        zd(jl,jlev) = .5_dp*dm1(jl,jlev,jrow) - d(jl,jlev,jrow)
     END DO
  END DO
  loperm = .FALSE.
  CALL conteq(zt,zp,zd,nbdim,nproma,loperm)

  IF (ltdiag) THEN
#ifndef MESSY
      pdigs(1:nglon,:,3,jrow) = 2._dp*zt(:,:)   ! explicit part of temperature
      pdsgs(1:nglon    ,jrow) = 2._dp*zp(:)     ! explicit part of pressure
#else
      ! mz_jb_20040829+
      ! explicit part of temperature: 
      pdigs(1:nproma,:,3,jrow) = 2.0_dp*zt(1:nproma,:)
      ! explicit part of pressure: 
      pdsgs(1:nproma,1  ,jrow) = 2.0_dp*zp(1:nproma)
      ! mz_jb_20040829-
#endif
  ENDIF

  !-- 3.3 Update *zt* and *zp* to compute the contribution
  !       of temperature and surface pressure to the
  !       divergence equation.

  DO jlev = 1, nlev
     DO jl = 1, nproma
        zt(jl,jlev) = zt(jl,jlev) + tm1(jl,jlev,jrow) - t(jl,jlev,jrow) +   &
                      zdt*tte(jl,jlev,jrow)
     END DO
  END DO
  DO jl = 1, nproma
     zp(jl) = zp(jl) + alpsm1(jl,jrow) - alps(jl,jrow) + zdt*alpste(jl,jrow)
  END DO
  CALL pgrad(zr,zt,nbdim,zp,nproma)

  ! explicit part of divergence
#ifndef MESSY
!DIR$ CONCURRENT
  IF (ltdiag) pdigs(1:nglon,:,2,jrow) = -zr(:,:)*ztwodt
#else
  ! mz_jb_20040505+: nproma introduced
  IF (ltdiag) pdigs(1:nproma,:,2,jrow) = -zr(1:nproma,:)*ztwodt
  ! mz_jb_20040505-
#endif

  !-- 3.4 Complete computation of the terms to be
  !       passed to *fftd*.

  DO jlev = 1, nlev
!DIR$ CONCURRENT
     DO jl = 1, nproma
        tm1(jl,jlev,jrow) = 2._dp*zt(jl,jlev) - tm1(jl,jlev,jrow) +    &
                                               2._dp*t(jl,jlev,jrow)
        rh(jl,jlev,jrow) = -ztwodt*(rh(jl,jlev,jrow)+zr(jl,jlev))
     END DO
  END DO

  DO jl = 1, nproma
     alpsm1(jl,jrow) = 2._dp*zp(jl) - alpsm1(jl,jrow) + 2._dp*alps(jl,jrow)
  END DO

  RETURN
END SUBROUTINE si1
