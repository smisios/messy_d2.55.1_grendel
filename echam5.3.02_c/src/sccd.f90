SUBROUTINE sccd

  ! Description:
  !
  ! This subroutine computes the final value of divergence
  !
  ! Method:
  !
  ! *sccd* is called from *stepon*
  !
  ! Results:
  ! The implicit equation is inverted with the help of
  ! temporary array *zd*
  !
  ! Reference:
  ! See echam3 manual  eq. 2.5.44
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, August 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! T. Diehl, DKRZ, July 1999, parallel version 
  !
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_memory_sp,     ONLY: sd
  USE mo_tmp_buffer,    ONLY: cn
  USE mo_control,       ONLY: ltdiag, nlev
  USE mo_diag_tendency, ONLY: pddiv

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: ic, is, jl, snsp
#if defined(__uxp__) || defined(__SX__) || defined(ES) || defined(_UNICOSMP)
  INTEGER :: i
#endif

  !  Local arrays: 
  REAL(dp) :: zd(nlev,2)
  INTEGER :: np1(lc%snsp)

  !  Intrinsic functions 
  INTRINSIC DOT_PRODUCT

  !  Executable statements 

  snsp = lc%snsp
  np1 = lc%np1

!-- 1. Invert divergence equation

!$OMP PARALLEL PRIVATE(is,ic,jl,zd)
!$OMP DO
  DO  is = 1, snsp

     ic = np1(is)
!CDIR SHORTLOOP
      DO jl = 1, nlev
#if defined(__uxp__) || defined(__SX__) || defined(ES) || defined(_UNICOSMP)
        zd(jl,1)=0.0_dp
        zd(jl,2)=0.0_dp
!CDIR SHORTLOOP
        DO i=1,nlev
          zd(jl,1) = zd(jl,1) + cn(i,jl,ic)*sd(i,1,is)
          zd(jl,2) = zd(jl,2) + cn(i,jl,ic)*sd(i,2,is)
        END DO
#else
        zd(jl,1) = DOT_PRODUCT(cn(1:nlev,jl,ic),sd(1:nlev,1,is))
        zd(jl,2) = DOT_PRODUCT(cn(1:nlev,jl,ic),sd(1:nlev,2,is))
#endif
      END DO

      ! store implicit part of divergence equation
      IF (ltdiag) pddiv(:,:,is,7) = pddiv(:,:,is,7) + zd(:,:) - sd(:,:,is)

      sd(:,:,is) = zd(:,:)

  END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE sccd
