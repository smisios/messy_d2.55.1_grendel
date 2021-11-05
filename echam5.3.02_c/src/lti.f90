!OCL NOALIAS

#ifdef CRAY
#define dgemm sgemm
#endif

#if defined(__uxp__) || defined(__SX__) || defined (ES) || defined(_UNICOSMP)
#define FAST_AND_DIRTY 1
#endif

#if defined (ES)
#define _DGEMM_MACRO_TEMPLATE 
#endif

SUBROUTINE lti

  ! Description:
  !
  ! Inverse Legendre transforms
  !
  ! Method:
  !
  ! This subroutine performs inverse *legendre transforms
  !
  ! *lti* is called from *scan2*
  !
  ! Results:
  ! *lti* computes the *fourier components
  ! for the current latitude line.
  !
  ! Authors:
  !
  ! D. W. Dent, ECMWF, December 1984, original source
  ! U. Schlese, DKRZ, December 1994, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! L. Kornblueh, MPI, November 2002, optimization
  ! L. Kornblueh, MPI, february 2004, optimization
  ! R. Smith, MPI, February 2004, optimization
  ! S. Shingu, NEC, March 2004, optimization (ES) 
  !
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp 
  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_memory_f,      ONLY: fad, fadu0, fatp, fatpm, fau, fau0, fav, favo, &
                              fsd, fsdu0, fstp, fstpm, fsu, fsu0, fsv, fsvo
  USE mo_legendre,      ONLY: leginv
  USE mo_decomposition, ONLY: lc => local_decomposition

  IMPLICIT NONE

  !  Local loop bounds

  INTEGER          :: nllev, nllevp1, nlmp1, nlnm0
  INTEGER ,POINTER :: nlmp(:), nlnp(:)

  !  Global bounds
  INTEGER          :: nhgl

  !  Local scalars:

  LOGICAL, SAVE :: ini_flag = .FALSE.

  INTEGER :: j, j0, ik, ims, inn, inp, ins, irow, iwrk, jl, jm
  LOGICAL :: lotypa

  !  Local arrays:

  REAL(dp) :: pnmit   (lc% nlat/2 ,lc% lnsp) ,&
              anmit   (lc% nlat/2 ,lc% lnsp) ,&
              pnmiuvt (lc% nlat/2 ,lc% lnsp) ,&
              anmiuvt (lc% nlat/2 ,lc% lnsp)

  REAL(dp) :: ftmp1 ( lc% nlat /2 , lc% nllev  , 2 ), &
              ftmp2 ( lc% nlat /2 , lc% nllev  , 2 ), &
              ftmp3 ( lc% nlat    , lc% nllev  , 2 ), &
              ftmp4 ( lc% nlat    , lc% nllev  , 2 ), &
              ftmpp1( lc% nlat    , lc% nllevp1, 2 )

  REAL(dp), POINTER :: fd (:,:,:,:), fdu0(:,:),     &
                       ftp(:,:,:,:), ftpm(:,:,:,:), &
                       fu0(:,:), fvo (:,:,:,:)

  REAL(dp), SAVE, ALLOCATABLE :: pnmil(:,:),pnmiuvl(:,:)
  REAL(dp), SAVE, ALLOCATABLE :: faud(:,:,:,:,:), favr(:,:,:,:,:), &
                                 fsur(:,:,:,:,:), fsvd(:,:,:,:,:)

  !  External subroutines
#ifndef _DGEMM_MACRO_TEMPLATE 
  EXTERNAL :: dgemm
#else
  EXTERNAL :: dgemmnt
#endif

  !  Intrinsic functions

  INTRINSIC MOD

  !  Executable statements

  !-- Set local loop bounds

  nllev   =  lc% nllev    ! number of levels
  nllevp1 =  lc% nllevp1  ! number of levels + 1
  nlmp1   =  lc% nlm      ! number of m wave numbers
  nlnm0   =  lc% nlnm0    ! number of coefficients for m=0

  nlmp    => lc% nlmp     ! displacement of the first point of columns
  nlnp    => lc% nlnp     ! number of points on each column
  nhgl    =  lc% nlat/2   ! global halv number of gaussian latitudes

  !-- Set initial values for transforms (zero coefficients (1,m=0 n=0)))

  IF (nlnm0 > 0) lvo(:,1,1) = 0.0_dp
  IF (nlnm0 > 0) ld (:,1,1) = 0.0_dp

  IF ( .NOT. ini_flag ) THEN

    ALLOCATE( faud ( lc% nlat /2 ,lc% nllev, 2, lc% nlm, 2))
    ALLOCATE( favr ( lc% nlat /2 ,lc% nllev, 2, lc% nlm, 2)) 
    ALLOCATE( fsur ( lc% nlat /2 ,lc% nllev, 2, lc% nlm, 2))
    ALLOCATE( fsvd ( lc% nlat /2 ,lc% nllev, 2, lc% nlm, 2))
    ALLOCATE( pnmil   (2*nhgl,lc% lnsp))
    ALLOCATE( pnmiuvl (2*nhgl,lc% lnsp))

    CALL leginv(pnmit=pnmit ,anmit=anmit ,pnmiuvt=pnmiuvt ,anmiuvt=anmiuvt)
!$OMP PARALLEL PRIVATE(irow)
!$OMP DO
    DO irow = 1, nhgl
      pnmil  (     irow,:) = pnmit (irow,:)
      pnmil  (nhgl+irow,:) = anmit (irow,:)
      pnmiuvl(     irow,:) = pnmiuvt(irow,:)
      pnmiuvl(nhgl+irow,:) = anmiuvt(irow,:)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ini_flag = .TRUE.

  ENDIF

  !-- Inverse *Legendre transforms

!$OMP PARALLEL PRIVATE(j0,j,iwrk,lotypa,ik,fd,ftp,ftpm,fvo,ims,inp,ins,jl, &
!$OMP                  ftmp1, ftmp2, ftmpp1, ftmp3, ftmp4)
!$OMP DO
  DO j0 = 1, nlmp1
    j = (nlmp1 - j0/2) * MOD(j0,2) + j0/2 * MOD(j0+1,2)

    DO iwrk = 0, 1
      lotypa = MOD(iwrk,2) == 0
      IF (lotypa) THEN
        ik = 1
        fd   => fsd
        ftp  => fstp
        ftpm => fatpm
        fvo  => fsvo
      ELSE
        ik = 2
        fd   => fad
        ftp  => fatp
        ftpm => fstpm
        fvo  => favo
      END IF

      ims  = nlmp(j) + 1         ! offset to local m columns  (spectral coef.)
      inp  = nlnp(j)             ! column length

      IF (lotypa) THEN
        ins  = (inp+1)/2
      ELSE
        ins  = inp/2
        ims  = ims  + 1
      END IF

      IF (ins>0) THEN
#ifdef _DGEMM_MACRO_TEMPLATE 
        CALL dgemmnt(  nhgl, 2*nllev,  ins, pnmil(1,ims), nhgl*2*2, ld (1,1,ims), nllev*2*2,  ftmp1 (1,1,1),   nhgl)
        CALL dgemmnt(  nhgl, 2*nllev,  ins, pnmil(1,ims), nhgl*2*2, lvo(1,1,ims), nllev*2*2,  ftmp2 (1,1,1),   nhgl)
        CALL dgemmnt(2*nhgl, 2*nllevp1,ins, pnmil(1,ims), nhgl*2*2, ltp(1,1,ims), nllevp1*2*2,ftmpp1(1,1,1), 2*nhgl)
#else
        CALL dgemm('N','T',  nhgl,2*nllev  ,ins,1.0_dp,pnmil(1,ims),nhgl*2*2,ld (1,1,ims),nllev*2*2  ,0.0_dp,ftmp1 (1,1,1),  nhgl)
        CALL dgemm('N','T',  nhgl,2*nllev  ,ins,1.0_dp,pnmil(1,ims),nhgl*2*2,lvo(1,1,ims),nllev*2*2  ,0.0_dp,ftmp2 (1,1,1),  nhgl)
        CALL dgemm('N','T',2*nhgl,2*nllevp1,ins,1.0_dp,pnmil(1,ims),nhgl*2*2,ltp(1,1,ims),nllevp1*2*2,0.0_dp,ftmpp1(1,1,1),2*nhgl)
#endif

#ifdef FAST_AND_DIRTY
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
        DO jl = 1, 2*nllev
          fd (jl,1,j,:) = ftmp1(:,jl,1)
          fvo(jl,1,j,:) = ftmp2(:,jl,1)
        ENDDO
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
        DO jl = 1, 2*nllevp1
          ftp (jl,1,j,:) = ftmpp1(     1:  nhgl,jl,1)
          ftpm(jl,1,j,:) = ftmpp1(nhgl+1:2*nhgl,jl,1)
        ENDDO
#else 
      DO jl = 1, nllev
        DO jm = 1,2
          fd (jl,jm,j,:) = ftmp1(:,jl,jm)
          fvo(jl,jm,j,:) = ftmp2(:,jl,jm)
        ENDDO
      ENDDO
      DO jl = 1, nllevp1
        DO jm = 1,2
          ftp (jl,jm,j,:) = ftmpp1(     1:  nhgl,jl,jm)
          ftpm(jl,jm,j,:) = ftmpp1(nhgl+1:2*nhgl,jl,jm)
        ENDDO
      ENDDO
#endif
      ELSE
        fd  (:,:,j,:) = 0.0_dp
        fvo (:,:,j,:) = 0.0_dp
        ftp (:,:,j,:) = 0.0_dp
        ftpm(:,:,j,:) = 0.0_dp
      END IF

      ims  = nlmp(j) + 1         ! offset to local m columns  (spectral coef.)
      inp  = nlnp(j)             ! column length

      IF (lotypa) THEN
        ins  = (inp+1)/2
      ELSE
        ins  = inp/2
        ims  = ims  + 1
      END IF

      IF (ins>0) THEN
#ifdef _DGEMM_MACRO_TEMPLATE 
        CALL dgemmnt(2*nhgl,2*nllev,ins,pnmiuvl(1,ims),nhgl*2*2,ld (1,1,ims),nllev*2*2,ftmp3,2*nhgl)
        CALL dgemmnt(2*nhgl,2*nllev,ins,pnmiuvl(1,ims),nhgl*2*2,lvo(1,1,ims),nllev*2*2,ftmp4,2*nhgl)
#else
        CALL dgemm('N','T',2*nhgl,2*nllev,ins,1.0_dp,pnmiuvl(1,ims),nhgl*2*2,ld (1,1,ims),nllev*2*2,0.0_dp,ftmp3,2*nhgl)
        CALL dgemm('N','T',2*nhgl,2*nllev,ins,1.0_dp,pnmiuvl(1,ims),nhgl*2*2,lvo(1,1,ims),nllev*2*2,0.0_dp,ftmp4,2*nhgl)
#endif
        faud(:,:,:,j,ik) = ftmp3(     1:  nhgl,:,:) ! Pnm * Dnm
        fsvd(:,:,:,j,ik) = ftmp3(nhgl+1:2*nhgl,:,:) ! d(Pnm)/d(mu) * Dnm
        favr(:,:,:,j,ik) = ftmp4(     1:  nhgl,:,:) ! Pnm * Znm
        fsur(:,:,:,j,ik) = ftmp4(nhgl+1:2*nhgl,:,:) ! d(Pnm)/d(mu) * Znm
      ELSE
        faud(:,:,:,j,ik) = 0.0_dp
        fsvd(:,:,:,j,ik) = 0.0_dp
        favr(:,:,:,j,ik) = 0.0_dp
        fsur(:,:,:,j,ik) = 0.0_dp
      END IF
    END DO

  END DO
!$OMP END DO
!$OMP END PARALLEL

  !-- Combine rotational and divergent parts of u and v

!DIR$ CONCURRENT
!DIR$ PREFERSTREAM
!$OMP PARALLEL PRIVATE(j,jl,irow)
!$OMP DO
  DO j = 1, nlmp1

#ifdef FAST_AND_DIRTY
    DO jl = 1, nllev
!DIR$ IVDEP
      DO irow = 1, nhgl
        fsu(jl,1,j,irow) =  fsur(irow,jl,1,j,2) + faud(irow,jl,2,j,1)
        fsu(nllev+jl,1,j,irow) =  fsur(irow,jl,2,j,2) - faud(irow,jl,1,j,1)
        fau(jl,1,j,irow) =  fsur(irow,jl,1,j,1) + faud(irow,jl,2,j,2)
        fau(nllev+jl,1,j,irow) =  fsur(irow,jl,2,j,1) - faud(irow,jl,1,j,2)
        fsv(jl,1,j,irow) =  favr(irow,jl,2,j,1) - fsvd(irow,jl,1,j,2)
        fsv(nllev+jl,1,j,irow) = -favr(irow,jl,1,j,1) - fsvd(irow,jl,2,j,2)
        fav(jl,1,j,irow) =  favr(irow,jl,2,j,2) - fsvd(irow,jl,1,j,1)
        fav(nllev+jl,1,j,irow) = -favr(irow,jl,1,j,2) - fsvd(irow,jl,2,j,1)
      END DO
    END DO
#else
    DO jl = 1, nllev
!DIR$ IVDEP
      DO irow = 1, nhgl
        ! Re{fs_u} = Re{A: d(Pnm)/d(mu) * Znm} + Im{S: Pnm * Dnm}
        fsu(jl,1,j,irow) =  fsur(irow,jl,1,j,2) + faud(irow,jl,2,j,1)

        ! Im{fs_u} = Im{A: d(Pnm)/d(mu) * Znm} - Re{S: Pnm * Dnm}
        fsu(jl,2,j,irow) =  fsur(irow,jl,2,j,2) - faud(irow,jl,1,j,1)

        !    1  2                      1  2                  1  2
        !    Re,Im                     Re,Im                 Re,Im
        !      |                         |                     |
        !      |                         |  1 2                |  1 2
        !      |                         |  S,A                |  S,A
        !      | m                       | m |                 | m |
        fau(jl,1,j,irow) =  fsur(irow,jl,1,j,1) + faud(irow,jl,2,j,2)
        fau(jl,2,j,irow) =  fsur(irow,jl,2,j,1) - faud(irow,jl,1,j,2)
        fsv(jl,1,j,irow) =  favr(irow,jl,2,j,1) - fsvd(irow,jl,1,j,2)
        fsv(jl,2,j,irow) = -favr(irow,jl,1,j,1) - fsvd(irow,jl,2,j,2)
        fav(jl,1,j,irow) =  favr(irow,jl,2,j,2) - fsvd(irow,jl,1,j,1)
        fav(jl,2,j,irow) = -favr(irow,jl,1,j,2) - fsvd(irow,jl,2,j,1)
      END DO
    END DO
#endif
  END DO
!$OMP END DO
!$OMP END PARALLEL
  
  !-- Transform mean wind

  IF ( nlnm0 > 0 ) THEN
!$OMP PARALLEL PRIVATE(ik,fu0,fdu0,inn)
!$OMP SECTIONS
!$OMP SECTION

    ik = 1
    fu0  => fsu0 
    fdu0 => fadu0
    inn  =  (nlnm0+1)/2
#ifdef _DGEMM_MACRO_TEMPLATE 
    CALL dgemmnt(nllev,nhgl,inn,lu0(1,ik),nllev*2,pnmil(1,ik),nhgl*2*2,fu0(1,1),nllev)
#else
    CALL dgemm('N','T',nllev,nhgl,inn,1.0_dp,lu0(1,ik),nllev*2,pnmil(     1,ik),nhgl*2*2,0.0_dp,fu0 (1,1),nllev)
#endif
!$OMP SECTION
    ik = 1
    fu0  => fsu0 
    fdu0 => fadu0
    inn  =  (nlnm0+1)/2
#ifdef _DGEMM_MACRO_TEMPLATE 
    CALL dgemmnt(nllev,nhgl,inn,lu0(1,ik),nllev*2,pnmil(nhgl+1,ik),nhgl*2*2,fdu0(1,1),nllev)
#else
    CALL dgemm('N','T',nllev,nhgl,inn,1.0_dp,lu0(1,ik),nllev*2,pnmil(nhgl+1,ik),nhgl*2*2,0.0_dp,fdu0(1,1),nllev)
#endif
!$OMP SECTION
    ik = 2
    fu0  => fau0 
    fdu0 => fsdu0
    inn  =  nlnm0/2
#ifdef _DGEMM_MACRO_TEMPLATE 
    CALL dgemmnt(nllev,nhgl,inn,lu0(1,ik),nllev*2,pnmil(1,ik),nhgl*2*2,fu0(1,1),nllev)
#else
    CALL dgemm('N','T',nllev,nhgl,inn,1.0_dp,lu0(1,ik),nllev*2,pnmil(     1,ik),nhgl*2*2,0.0_dp,fu0 (1,1),nllev)
#endif
!$OMP SECTION
    ik = 2
    fu0  => fau0 
    fdu0 => fsdu0
    inn  =  nlnm0/2
#ifdef _DGEMM_MACRO_TEMPLATE 
    CALL dgemmnt(nllev,nhgl,inn,lu0(1,ik),nllev*2,pnmil(nhgl+1,ik),nhgl*2*2,fdu0(1,1),nllev)
#else
    CALL dgemm('N','T',nllev,nhgl,inn,1.0_dp,lu0(1,ik),nllev*2,pnmil(nhgl+1,ik),nhgl*2*2,0.0_dp,fdu0(1,1),nllev)
#endif
!$OMP END SECTIONS
!$OMP END PARALLEL
  END IF

END SUBROUTINE lti

#ifdef _DGEMM_MACRO_TEMPLATE

! keep this to allow compiler to use intrinsic macro expansion

SUBROUTINE DGEMMNT ( m, n, k, a, lda, b, ldb, c, ldc )

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  !     .. Scalar Arguments ..
  INTEGER :: m, n, k, lda, ldb, ldc
  !     .. Array Arguments ..
  REAL(dp) ::  a(lda,*), b(ldb,*), c(ldc,*)
  !     .. Local Scalars ..
  INTEGER :: i, j, l
  !     .. Parameters ..
  REAL(dp), PARAMETER ::  one  = 1.0_dp, &
                          zero = 0.0_dp
  !
  !     Form  C = C+A*B'
  !
  DO j = 1, n
    DO i = 1, m
      c( i, j ) = zero
    ENDDO
  ENDDO
  DO j = 1, n
    DO l = 1, k
      IF( b( j, l ) /= zero )THEN
        DO i = 1, m
          c( i, j ) = c( i, j ) + b( j, l )*a( i, l )
        ENDDO
      END IF
    ENDDO
  ENDDO
  
  RETURN
  
END SUBROUTINE dgemmnt
#endif



