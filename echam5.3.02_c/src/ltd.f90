!OCL NOALIAS

#ifdef CRAY
#define dgemv sgemv
#define dgemm sgemm
#endif

#if defined(__uxp__) || defined(__SX__)  || defined(ES) || defined(_UNICOSMP)
#define FAST_AND_DIRTY 1
#endif

SUBROUTINE ltd

  ! Description:
  !
  ! Direct Legendre transforms.
  !
  ! Method:
  !
  ! This subroutine performs direct *Legendre transforms for
  ! the divergence equation,
  ! the temperature and surface pressure equations,
  ! the vorticity equation,
  ! the mean wind.
  !
  ! *ltd* is called from *scan1sl*
  !
  ! Results:
  ! *ltd* adds in the spectral arrays:-
  !      *ld*  the contribution of the current latitude line
  !      *ltp* the contribution of the current latitude line
  !      *lvo* the contribution of the current latitude line
  !      *lu0* the contribution of the current latitude line
  !
  ! Authors:
  !
  ! D. W. Dent, ECMWF, December 1984, original source
  ! U. Schlese, DKRZ, in 1991, and 1994, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! L. Kornblueh, MPI, November 2002, optimization
  !
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_memory_f,      ONLY: fadl, fadm, far, fatp1, faul, fazl, fazm, fsdl, &
                              fsdm, fsr, fstp1, fsul, fszl, fszm
  USE mo_legendre,      ONLY: legmod
  USE mo_decomposition, ONLY: lc => local_decomposition

  IMPLICIT NONE

  !  Local loop bounds

  INTEGER          :: nllev, nllevp1, nlmp1, nlnm0, lnsp
  INTEGER ,POINTER :: nlmp(:), nlnp(:)

  !  Local scalars:

  INTEGER :: ims, ins, inu, is, iu, j, jh, jn, jhr, nhgl, kk, j0
#ifndef FAST_AND_DIRTY
  INTEGER :: kf, km
#endif  

  REAL(dp), POINTER :: ful (:,:)

  REAL(dp), POINTER :: fdl2 (:,:,:,:)
  REAL(dp), POINTER :: fdm2 (:,:,:,:)
  REAL(dp), POINTER :: fr2  (:,:,:,:)
  REAL(dp), POINTER :: ftp12(:,:,:,:)
  REAL(dp), POINTER :: fzl2 (:,:,:,:)
  REAL(dp), POINTER :: fzm2 (:,:,:,:)

  REAL(dp), SAVE, ALLOCATABLE :: pnml(:,:)
  REAL(dp), SAVE, ALLOCATABLE :: anml(:,:)
  REAL(dp), SAVE, ALLOCATABLE :: rnml(:,:)

  REAL(dp), SAVE, ALLOCATABLE :: fpstack(:,:,:)
  REAL(dp), SAVE, ALLOCATABLE :: fastack(:,:,:)
  REAL(dp), SAVE, ALLOCATABLE :: lpstack(:,:)
  REAL(dp), SAVE, ALLOCATABLE :: lastack(:,:)
  REAL(dp), SAVE, ALLOCATABLE :: lrstack(:,:)

  !  External subroutines 
  EXTERNAL dgemv, dgemm

  !  Executable statements

  !-- Set local loop bounds

  nllev   =  lc% nllev    ! number of levels
  nllevp1 =  lc% nllevp1  ! number of levels + 1
  nlmp1   =  lc% nlm      ! number of m wave numbers
  nlnm0   =  lc% nlnm0    ! number of coefficients for m=0
  lnsp    =  lc% lnsp     ! number of complex spectral coefficients on this pe
  nlmp    => lc% nlmp     ! displacement of the first point of columns
  nlnp    => lc% nlnp     ! number of points on each column
  nhgl    =  lc% nlat/2   ! global halv number of gaussian latitudes

  IF ( .NOT. ALLOCATED(fpstack) ) THEN

    ALLOCATE(fpstack((2*nllev+nllevp1)*2,nlmp1,nhgl))
    ALLOCATE(fastack(2*nllev*2,nlmp1,nhgl))
    ALLOCATE(lpstack((2*nllev+nllevp1)*2+1,lnsp))
    ALLOCATE(lastack(2*nllev*2+1,lnsp))
    ALLOCATE(lrstack(nllev*2+1,lnsp))

    ALLOCATE( anml( lc% nlat/2, lc% lnsp))
    ALLOCATE( pnml( lc% nlat/2, lc% lnsp))
    ALLOCATE( rnml( lc% nlat/2, lc% lnsp))

    ! derive local wavenumber index
    ! calculate legendre coefficents for each latitude

    CALL legmod(pnmt=pnml,anmt=anml,rnmt=rnml)
  END IF

  !-- 1. Legendre transforms

  ! Loop over hemispheres (1: north, 2: south)

  DO jh = 1, 2
    iu = 2 - jh

    IF (jh==1) THEN
      fdl2  => fsdl (:,:,:,:)
      fdm2  => fadm (:,:,:,:)
      fr2   => fsr  (:,:,:,:)
      fzl2  => fszl (:,:,:,:)
      fzm2  => fazm (:,:,:,:)
      ftp12 => fstp1(:,:,:,:)
    ELSE
      fdl2  => fadl (:,:,:,:)
      fdm2  => fsdm (:,:,:,:)
      fr2   => far  (:,:,:,:)
      fzl2  => fazl (:,:,:,:)
      fzm2  => fszm (:,:,:,:)
      ftp12 => fatp1(:,:,:,:)
    END IF

    !-- 1.1 Transforms for d, vo, t and p

    ! Loop over latitudes

    DO jhr = 1, nhgl
      DO j=1,nlmp1
!CDIR NODEP
#ifdef FAST_AND_DIRTY        
        DO kk=1,nllev*2
          fpstack(kk,j,jhr)             = fdl2(kk,1,j,jhr)
          fpstack(kk+nllev*2,j,jhr)     = fzl2(kk,1,j,jhr)
          fastack(kk,j,jhr)             = fdm2(kk,1,j,jhr)
          fastack(kk+nllev*2,j,jhr)     = fzm2(kk,1,j,jhr)
        ENDDO
        DO kk=1,nllevp1*2
          fpstack(kk+nllev*4,j,jhr)     = ftp12(kk,1,j,jhr)
        ENDDO
#else
        DO kf = 1,2
!DIR$ CONCURRENT
          DO kk = 1, nllev
            km = (kf-1)*nllev+kk 
            fpstack(km,j,jhr)             = fdl2(kk,kf,j,jhr)
            fpstack(km+nllev*2,j,jhr)     = fzl2(kk,kf,j,jhr)
            fastack(km,j,jhr)             = fdm2(kk,kf,j,jhr)
            fastack(km+nllev*2,j,jhr)     = fzm2(kk,kf,j,jhr)
          ENDDO
        ENDDO
        DO kf = 1,2
          DO kk=1,nllevp1
            km = (kf-1)*nllevp1+kk 
            fpstack(km+nllev*4,j,jhr)     = ftp12(kk,kf,j,jhr)
          ENDDO
        ENDDO
#endif
      ENDDO

    END DO

!CSD$ PARALLEL DO PRIVATE(j,ims,ins,jn,is)
!$OMP PARALLEL PRIVATE(j,ims,ins,jn,is,j0)
!$OMP DO
    DO j0 = 1, nlmp1
      j = (nlmp1 - j0/2) * MOD(j0,2) + j0/2 * MOD(j0+1,2)
      ims  = nlmp(j) - iu   ! offset to local m columns  (spectral coef.)
      ins  = nlnp(j) + iu   ! column length

      DO jn = 2, ins, 2
        is  = ims  + jn
        CALL dgemv('N',4*nllev+2*nllevp1,nhgl,  1._dp,fpstack(1,j,1),(4*nllev+2*nllevp1)*nlmp1,pnml(1,is),1,0._dp,lpstack(1,is),1)
        CALL dgemv('N',4*nllev          ,nhgl, -1._dp,fastack(1,j,1),4*nllev*nlmp1,            anml(1,is),1,0._dp,lastack(1,is),1)
        CALL dgemv('N',2*nllev          ,nhgl,  1._dp,fr2(1,1,j,1),  2*nllev*nlmp1,            rnml(1,is),1,0._dp,lrstack(1,is),1)
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO

  END DO

  !-- 1.2 Transforms for the mean wind

  IF ( nlnm0 > 0 ) THEN
!CSD$ PARALLEL DO PRIVATE(iu,inu,ful)
!$OMP PARALLEL PRIVATE(iu,inu,ful,jh)
!$OMP DO
    DO jh = 1,2
      iu = 2 - jh
      
      inu = (nlnm0+iu)/2
      
      IF (jh==1) THEN
        ful  => fsul 
      ELSE
        ful  => faul
      END IF

      CALL dgemm("N","N",nllev,inu,nhgl,1.0_dp,ful(1,1),nllev,pnml(1,2-iu),2*nhgl,1.0_dp,lu0(1,2-iu),2*nllev)
      
    END DO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO
  END IF

#ifdef FAST_AND_DIRTY        
!DIR$ CONCURRENT 
!DIR$ PREFERVECTOR
  DO kk=1,nllev*2
    ld(kk,1,:)=lpstack(kk,:)+lastack(kk,:)+lrstack(kk,:)
    lvo(kk,1,:)=lpstack(2*nllev+kk,:)+lastack(2*nllev+kk,:)
  ENDDO
!DIR$ CONCURRENT 
!DIR$ PREFERVECTOR
  DO kk=1,nllevp1*2
    ltp(kk,1,:)=lpstack(4*nllev+kk,:)
  ENDDO
#else
  DO kf = 1,2
!DIR$ CONCURRENT
    DO kk = 1, nllev
      km = (kf-1)*nllev+kk 
      ld(kk,kf,:)=lpstack(km,:)+lastack(km,:)+lrstack(km,:)
      lvo(kk,kf,:)=lpstack(2*nllev+km,:)+lastack(2*nllev+km,:)
    ENDDO
  ENDDO
  DO kf = 1,2
    DO kk=1,nllevp1
      km = (kf-1)*nllevp1+kk 
      ltp(kk,kf,:)=lpstack(4*nllev+km,:)
    ENDDO
  ENDDO
#endif

END SUBROUTINE ltd
