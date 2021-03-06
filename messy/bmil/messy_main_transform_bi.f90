! **************************************************************************
MODULE messy_main_transform_bi
! **************************************************************************
#if defined(ECHAM5) || defined(MBM_CLAMS)
!#define MPIPERF_BAD
!#define _MPI_ALLGATHER
#define _MPI_PSUM
  !
  ! This module contains high-level routines for the
  ! transformation and transposition
  ! of decomposed fields between the different
  ! representations / decompositions.
  !
  ! The transposition routines
  !    trp_gp_fs_3d, trp_fs_ls_3d, trp_ls_sp_3d
  ! are based on the respective routines
  ! in mo_transpose.f90 (by Andreas Rhodin, DWD/MPIM, April     2002)
  !
  ! Transpositions in grid-point space
  !   indx, reorder, gather_gp, scatter_gp
  ! are used from mo_transpose.f90.
  !
  ! TRANSFORMATIONS BETWEEN SPECTRAL- AND GRIDPOINT-SPACE IN ONE STEP
  ! (FROM/TO RESPECTIVE DECOMPOSITION):
  !
  !                               sp2gp
  !   ---> trp_ls_sp_3d(-1, ...)    |   transpos. spec. space -> Leg. space
  !   ---> trf_lti                  |   inverse Legendre-Trafo
  !   ---> trp_fs_as(-1, ...)       |   antisym. + sym. -> Fourier
  !   ---> trp_fs_ls_3d(-1, ...)    |   transpos. Leg. space -> Four. space
  !   ---> trf_ffti                 |   inverse Fourier-Trafo
  !   ---> trp_gp_fs_3d(-1, ...)    V   transpos. Four. space -> gridpoint
  !
  !   ---> trp_ls_sp_3d(+1, ...)    ^   transpos. Leg. space -> spec. space
  !   ---> trf_ltd                  |   direct Legendre-Trafo
  !   ---> trp_fs_as(+1, ...)       |   Fourier -> antisym. + sym.
  !   ---> trp_fs_ls_3d(+1, ...)    |   transpos. Four. space -> Leg. space
  !   ---> trf_fftd                 |   direct Fourier-Trafo
  !   ---> trp_gp_fs_3d(+1, ...)    |   transpos. gridpoint -> Four. space
  !                               gp2sp
  !
  ! NEW: (1) TRANSPOSITION BETWEEN DECOMPOSED GRIDPOINT-SPACE AND
  !          GLOBAL GRIDPOINT-SPACE (VISIBLE ON ALL PEs):
  !           --> trp_gpdc_gpgl
  !
  !      (2) INDEX-DECOMPOSITION AND RECOMPOSITION,
  !          I.E., - DISTRIBUTION OF SUBSETS OF GLOBAL FIELDS FROM THE
  !                  I/O-PE TO DIFFERENT PEs
  !                - COMBINATION OF SUBSETS FROM DIFFERENT PEs TO GLOBAL
  !                  FIELDS ON THE I/O-PE
  !           --> get_dc_index, scatter_glix, gather_glix
  !
  !      Those can be used for processes, where a different decomposition
  !      makes more sense, e.g.,
  !       - parallelization of advection algorithms
  !         using the number of tracers
  !       - parallelization of ATTILA
  !         using the number of air parcels
  !
  ! AUTHORS:
  !     Michael Traub and Patrick J?ckel, MPICH, Jul 2003
  !     Patrick Joeckel,                  MPICH, Sep 2003
  !     Patrick Joeckel,                  MPICH, Feb 2004
  !     Patrick Joeckel,                  MPICH, Jul 2005
  !
  ! -----------------------------------------------------------------------
  ! declare explicit shape arguments in routine (un)pack_fs_buf
  ! on vector machines
  !
#if ((defined CRAY)&&(! defined _CRAYMPP)) || (defined __SX__) || (defined __uxp__)
#define EXPLICIT
#else
#undef EXPLICIT
#endif
  !
#ifdef CRAY
#define dgemm sgemm
#endif
  !
#if defined(__uxp__) || defined(__SX__)
#define FAST_AND_DIRTY 1
#endif
  ! -----------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather_bi,    ONLY: error_bi
  USE messy_main_grid_def_mem_bi
  USE messy_main_mpi_bi,        ONLY: p_send, p_recv, p_sendrecv &
                                    , p_pe, p_io                 &
                                    , p_bcast, p_nprocs          &
                                    , reorder
#ifndef MBM_CLAMS
  USE messy_main_grid_def_bi,   ONLY: philon, philat
  USE messy_main_mpi_bi,        ONLY: dcg, dcl, indx
#endif


  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  ! FLAGS FOR trp_gpdc_gpgl(-1, lf, gf, FLAG)
  INTEGER, PARAMETER, PUBLIC :: M_SUM = 0  ! SUM OVER ALL PEs
  INTEGER, PARAMETER, PUBLIC :: M_AVE = 1  ! AVERAGE OVER ALL PEs
  INTEGER, PARAMETER, PUBLIC :: M_STD = 2  ! STANDARD DEV. BETWEEN ALL PEs
  INTEGER, PARAMETER, PUBLIC :: M_LOC = 3  ! EXTRACT ONLY LOCAL

#if defined (ECHAM5) || defined (MBM_CLAMS)

#ifdef MBM_CLAMS
  INTEGER, PARAMETER :: tag_scatter_ix   = 310
  INTEGER, PARAMETER :: tag_gather_ix    = 320
  PUBLIC :: get_dc_index   ! get decomposition along index
  PUBLIC :: scatter_glix   ! scatter field along index between different PEs
  PUBLIC :: gather_glix    ! collect field along index from different PEs
  INTERFACE scatter_glix
     MODULE PROCEDURE scatter_glix_1d
     MODULE PROCEDURE scatter_glix_4d
  END INTERFACE
  !
  INTERFACE gather_glix
     MODULE PROCEDURE gather_glix_4d
     MODULE PROCEDURE gather_glix_1d
  END INTERFACE
#else

  ! DECOMPOSITION INFORMATION
  PUBLIC :: get_mn_ls      ! get (m(k),n(k)) - lists of spectral coefficients
  PUBLIC :: get_dc_index   ! get decomposition along index

  ! TRANSFORMATIONS USING TRANSPOSITIONS
  PUBLIC :: sp2gp          ! spectral to gridpoint
  PUBLIC :: gp2sp          ! gridpoint to spectral
  PUBLIC :: trf_ls_dz2uv   ! transform divergence + vorticity to (u,v)*cos(lat)
  !                        ! ... in legendre space
  !
  INTERFACE trf_ls_dz2uv
     MODULE PROCEDURE trf_ls_dz2uv_0
     MODULE PROCEDURE trf_ls_dz2uv_1
  END INTERFACE

  ! TRANSPOSITIONS
  INTEGER, PARAMETER :: tag_tr_fs_ls     = 210
  INTEGER, PARAMETER :: tag_tr_ls_sp     = 220
  INTEGER, PARAMETER :: tag_tr_gpdc_gpgl = 300
  INTEGER, PARAMETER :: tag_scatter_ix   = 310
  INTEGER, PARAMETER :: tag_gather_ix    = 320
  !
  PUBLIC :: trp_gp_fs_3d   ! grid point space <-> fourier space
  PUBLIC :: trp_fs_ls_3d   ! fourier space    <-> legendre space
  PUBLIC :: trp_ls_sp_3d   ! legendre space   <-> spectral space
  PUBLIC :: trp_fs_as      ! fourier space    <-> asymmetric + symmetric comp.
  PUBLIC :: trp_gpdc_gpgl  ! decomposed gridpoint space  <->
                           !                        global gp-space on each PE
  PUBLIC :: scatter_glix   ! scatter field along index between different PEs
  PUBLIC :: gather_glix    ! collect field along index from different PEs
  !
  INTERFACE trp_fs_as
     MODULE PROCEDURE trp_fs_as_0
     MODULE PROCEDURE trp_fs_as_1
  END INTERFACE
  !
  INTERFACE trp_gpdc_gpgl
     MODULE PROCEDURE trp_gpdc_gpgl_4d
     MODULE PROCEDURE trp_gpdc_gpgl_3d
     MODULE PROCEDURE trp_gpdc_gpgl_2d
  END INTERFACE

   INTERFACE scatter_glix
     MODULE PROCEDURE scatter_glix_1d
     MODULE PROCEDURE scatter_glix_4d
  END INTERFACE
  !
  INTERFACE gather_glix
     MODULE PROCEDURE gather_glix_4d
     MODULE PROCEDURE gather_glix_1d
  END INTERFACE

  ! TRANSFORMATIONS
  PUBLIC :: trf_lti        ! inverse Legendre transformation
  PUBLIC :: trf_ltd        ! direct Legendre transformation
  PUBLIC :: trf_ffti       ! inverse Fourier transformation
  PUBLIC :: trf_fftd       ! direct Fourier transformation
  !
  INTERFACE trf_lti
     MODULE PROCEDURE trf_lti_0
     MODULE PROCEDURE trf_lti_1
     MODULE PROCEDURE trf_lti_2
  END INTERFACE
  !
  INTERFACE trf_ltd
     MODULE PROCEDURE trf_ltd_0
     MODULE PROCEDURE trf_ltd_1
     MODULE PROCEDURE trf_ltd_2
  END INTERFACE

  !PRIVATE :: set_ilon_ilat
  INTEGER, DIMENSION(:,:), POINTER, SAVE :: ilon => NULL()
  INTEGER, DIMENSION(:,:), POINTER, SAVE :: ilat => NULL()

  PUBLIC :: trf_ltd_a
  INTERFACE trf_ltd_a
     MODULE PROCEDURE trf_ltd_a_0
     MODULE PROCEDURE trf_ltd_a_1
     MODULE PROCEDURE trf_ltd_a_2
  END INTERFACE
  PUBLIC :: trf_ltd_r
  INTERFACE trf_ltd_r
     MODULE PROCEDURE trf_ltd_r_0
     MODULE PROCEDURE trf_ltd_r_1
     MODULE PROCEDURE trf_ltd_r_2
  END INTERFACE

  INTERFACE zonal_average
     ! CALCULATES ZONAL AVERAGE IN DECOMPOSITION
     MODULE PROCEDURE zonal_average_3d
     MODULE PROCEDURE zonal_average_2d
  END INTERFACE
  PUBLIC :: zonal_average       ! calculate zonal average (2D, 3D)
  !                             ! for any point in decomposition (GP)

  INTERFACE locate_in_decomp
     MODULE PROCEDURE locate_in_decomp_nrst
     MODULE PROCEDURE locate_in_decomp_4
  END INTERFACE
  PUBLIC :: locate_in_decomp
! ##########################################################################
! END #ifdef MBM_CLAMS ... #else
! ##########################################################################
#endif
#endif

CONTAINS

#if defined (ECHAM5) || defined(MBM_CLAMS)
  !==========================================================================

#ifndef MBM_CLAMS
  SUBROUTINE set_ilon_ilat

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NINT, REAL, SPREAD

    ! LOCAL
    INTEGER                           :: jlat, jlon
    REAL(dp)  :: tmp (dcl% nglon, dcl% nglat)
    REAL(dp)  :: tmp2(dcl% nproma, dcl% ngpblks)

    IF (ASSOCIATED(ilon) .AND. ASSOCIATED(ilat)) RETURN
    IF (ASSOCIATED(ilon)) THEN
       DEALLOCATE(ilon); NULLIFY(ilon)
    END IF
    IF (ASSOCIATED(ilat)) THEN
       DEALLOCATE(ilat)
       NULLIFY(ilat)
    END IF

    ALLOCATE (ilon (dcl% nproma, dcl% ngpblks))
    ALLOCATE (ilat (dcl% nproma, dcl% ngpblks))

    tmp = REAL(SPREAD (dcl% glat(:), 1, dcl%nglon),dp)
    CALL reorder (tmp2 ,tmp)
    ilat = NINT(tmp2)

    DO jlat = 1, dcl% nglat     ! local  latitude index N -> S
      DO jlon = 1, dcl% nglon   ! local  longitude index
        tmp(jlon,jlat) = REAL(jlon + dcl% glon(jlat),dp)
      END DO
    END DO
    CALL reorder (tmp2 ,tmp)
    ilon = NINT(tmp2)

  END SUBROUTINE set_ilon_ilat
  !==========================================================================

  !==========================================================================
  SUBROUTINE get_mn_ls(m, n)

    ! GENERATE LIST OF (m, n) IN LEGENDRE SPACE

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER, DIMENSION(:), POINTER :: m, n   ! INTENT(OUT)

    ! LOCAL
    INTEGER :: jsp    ! counter for spectral coeff. pairs (m, n)
    INTEGER :: jm, jn, jjn

    IF (ASSOCIATED(m)) THEN
       DEALLOCATE(m) ; NULLIFY(m)
    END IF
    IF (ASSOCIATED(n)) THEN
       DEALLOCATE(n) ; NULLIFY(m)
    END IF

    ! (LEGENDRE SPACE)
    ALLOCATE(m(dcl%lnsp))
    ALLOCATE(n(dcl%lnsp))

    jsp = 0
    DO jm=1, dcl%nlm
       jjn = 0
       DO jn=1, dcl%nlnp(jm)
          jsp = jsp + 1
          m(jsp) = dcl%lm(jm)
          n(jsp) = dcl%lm(jm) + jjn
          jjn = jjn + 1
       END DO
    END DO

  END SUBROUTINE get_mn_ls
  !==========================================================================

  !==========================================================================
  SUBROUTINE trp_gp_fs_3d (sign, pgp, pfs)
    !
    ! transpose
    !   sign= 1 : grid point space  -> Fourier space
    !   sign=-1 : grid point space <-  Fourier space
    !

    IMPLICIT NONE

    INTRINSIC :: MOD, ASSOCIATED, SIZE

    ! I/O
    INTEGER              ,INTENT(in)     :: sign        ! 1:gp>fs; -1:gp<fs
    REAL(dp)             ,INTENT(inout)  :: pgp(:,:,:)  ! gridpoint space 3d
    REAL(dp)             ,INTENT(inout)  :: pfs(:,:,:)  ! Fourier space

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'trp_gp_fs_3d'
    INTEGER             :: i, j, k, l
    INTEGER             :: imype                  ! index of this pe
    REAL(dp),POINTER    :: gp(:,:,:)    => NULL()  ! gridpoint
    REAL(dp),POINTER    :: fs(:,:,:)    => NULL()  ! fourier
    REAL(dp),POINTER    :: bufs (:,:,:) => NULL()  ! send buffer
    REAL(dp),POINTER    :: bufr(:,:,:)  => NULL()  ! receive buffer
    INTEGER             :: seta          ! my set A
    INTEGER             :: nprocb        ! number of PEs in set A
    INTEGER,ALLOCATABLE :: idx_com (:)   ! PEs to communicate with
    INTEGER             :: ks, ke, nk    ! vertical range of buffer
    INTEGER             :: nglat, nglon  ! gridpoint space no. lons,lats
    INTEGER             :: glons(2)      ! first longitudes in gridspace
    INTEGER             :: glone(2)      ! last  longitudes in gridspace
    INTEGER             :: nglh(2)       ! number of lats in each domains
    INTEGER             :: nlon          ! global number of longitudes
    LOGICAL             :: lreg          ! regular lon/lat ordering

    ALLOCATE( gp( SIZE(pgp,1), dcl%nlev +1, SIZE(pgp,3))); gp(:,:,:) = 0.0_dp
    ALLOCATE( fs( SIZE(pfs,1), dcl%nflevp1, SIZE(pfs,3))); fs(:,:,:) = 0.0_dp
    !
    ! loop 1: send if (dest_pe > src_pe); else recv
    !
    imype  = indx (p_pe, dcg)
    seta   = dcg(imype)% set_a
    nprocb = dcg(imype)% nprocb
    lreg   = dcg(imype)% lreg
    IF (dcg(imype)% col_1d) RETURN
    ALLOCATE (idx_com (0:nprocb-1))
    CALL plan_b_sr (idx_com)
    i = imype         ! PE index (my)
    SELECT CASE (sign)
    CASE (1)
       !
       ! grid point -> Fourier
       !
       gp(:,:size(pgp,2),:) = pgp(:,:,:)
       DO k = 0, nprocb-1
          j = idx_com (           k        ) ! PE index (send)
          l = idx_com (MOD(nprocb-k,nprocb)) ! PE index (recv)
          CALL alloc_gp_fs_buf (gp_i=i, fs_i=j, zbuf=bufs)
          CALL pack_gp_buf
          IF(i/=j) CALL alloc_gp_fs_buf (gp_i=l, fs_i=i, zbuf=bufr)
          CALL send_recv
          CALL unpack_buf_fs
          IF (ASSOCIATED(bufr)) THEN
             DEALLOCATE(bufr) ; NULLIFY(bufr)
          END IF
       END DO
       !
       ! zero latitudes > nlat
       !
       fs (dcg(imype)% nlon+1:,:,:) = 0.
       pfs(:,:,:) = fs(:,:size(pfs,2),:)
    CASE (-1)
       !
       ! Fourier -> grid point
       !
       fs(:,:size(pfs,2),:) = pfs(:,:,:)
       DO k = 0, nprocb-1
          j = idx_com (           k        ) ! PE index (send)
          l = idx_com (MOD(nprocb-k,nprocb)) ! PE index (recv)
          CALL alloc_gp_fs_buf (gp_i=j, fs_i=i, zbuf=bufs)
          CALL pack_fs_buf
          IF(i/=j) CALL alloc_gp_fs_buf (gp_i=i, fs_i=l, zbuf=bufr)
          CALL send_recv
          CALL unpack_buf_gp
          IF (ASSOCIATED(bufr)) THEN
             DEALLOCATE(bufr) ; NULLIFY(bufr)
          END IF
       END DO
       pgp(:,:,:) = gp(:,:size(pgp,2),:)
    CASE default
       CALL error_bi ('invalid SIGN parameter (not 1,-1)', substr)
    END SELECT
    !
    ! frt at ECMWF doesnt yet deallocate
    !
    DEALLOCATE(idx_com)
    DEALLOCATE(gp,fs)

  CONTAINS
    !---------------------------------------------------------------------
    SUBROUTINE alloc_gp_fs_buf (gp_i, fs_i, zbuf)
      INTEGER       :: gp_i, fs_i ! indices of grid point ...
      ! ... and Fourier space pe
      REAL(dp),POINTER :: zbuf(:,:,:)
      !
      ! derive bounds and allocate buffer
      !
      nk    = dcg(fs_i)%nflevp1
      ks    = dcg(fs_i)% flevs
      ke    = dcg(fs_i)% fleve
      nglat = dcg(gp_i)% nglat
      nglon = dcg(gp_i)% nglon
      glons = dcg(gp_i)% glons
      glone = dcg(gp_i)% glone
      nlon  = dcg(gp_i)% nlon
      nglh  = dcg(gp_i)% nglh
      ALLOCATE (zbuf (nglon, nk, nglat) )
    END SUBROUTINE alloc_gp_fs_buf
    !--------------------------------------------------------------------
    SUBROUTINE pack_gp_buf
      !
      ! pack message to send/recv buffer buf
      !
      ! pack 3d arrays
      !
      IF(nk > 0) THEN
         IF (lreg) THEN
            bufs (:,:nk,:) = gp (:nglon,ks:ke,:)
         ELSE
            CALL reorder (bufs (:,:nk,:), gp (:,ks:ke,:))
         ENDIF
      ENDIF
    END SUBROUTINE pack_gp_buf
    !-------------------------------------------------------------------
    SUBROUTINE unpack_buf_gp
      !
      ! unpack grid point space from send/recv buffer buf
      !
      ! unpack 3d arrays
      !
      IF (ke == dcg(imype)% nlev+1) THEN
        ke = ke - 1
        nk = nk - 1
      ENDIF
      !
      IF(nk > 0) THEN
         IF (lreg) THEN
            gp (:,ks:ke,:) = bufr (:,:nk,:)
         ELSE
            CALL reorder (gp (:,ks:ke,:), bufr (:,:nk,:))
         ENDIF
      ENDIF
      !
    END SUBROUTINE unpack_buf_gp
    !-------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs
      !
      ! unpack message to fourier buffer fs
      !
      !
      ! unpack first segment
      !
      fs(glons(1):glone(1),:,:nglat/2) = bufr(:,:,:nglat/2)
      !
      ! unpack second segment
      !
      IF (glone(2)>glons(2)) THEN
         fs(glons(2):glone(2),:,nglat/2+1:) = bufr(:,:,nglat/2+1:)
      ELSE
         !
         ! unpack second segment, split into longitudes
         !
         fs    (glons(2)        :nlon           ,:,nglat/2+1:) = &
              bufr(                :nlon-glons(2)+1,:,nglat/2+1:)

         fs    (1               :glone(2)       ,:,nglat/2+1:) = &
              bufr(nglon-glone(2)+1:               ,:,nglat/2+1:)
      ENDIF
    END SUBROUTINE unpack_buf_fs
    !------------------------------------------------------------------
    SUBROUTINE pack_fs_buf
      !
      ! pack fourier buffer fs to buffer
      !
      !
      ! pack first segment
      !
      bufs(:,:,:nglh(1)) = fs(glons(1):glone(1),:,:nglh(1))
      !
      ! pack second segment
      !
      IF(nglh(1)>0) THEN
         IF (glone(2)>glons(2)) THEN
            bufs(:,:,nglh(1)+1:) = fs(glons(2):glone(2),:,nglh(1)+1:)
         ELSE
            !
            ! pack second segment, split into longitudes
            !
            bufs  (                :nlon-glons(2)+1,:,nglh(1)+1:) = &
                 fs (glons(2)        :nlon           ,:,nglh(1)+1:)
            bufs  (nglon-glone(2)+1:               ,:,nglh(1)+1:) = &
                 fs (1               :glone(2)       ,:,nglh(1)+1:)
         ENDIF
      ENDIF
    END SUBROUTINE pack_fs_buf
    !-----------------------------------------------------------------
    SUBROUTINE send_recv
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      IF(i/=j) THEN
         CALL p_sendrecv (bufs,  dcg(j)% pe, &
              bufr, dcg(l)% pe, tag_tr_ls_sp)
         DEALLOCATE (bufs) ; NULLIFY(bufs)
      ELSE
         bufr => bufs
      ENDIF
    END SUBROUTINE send_recv
    !----------------------------------------------------------------
    SUBROUTINE plan_b_sr (idx_com)
      !
      ! set up plan of PEs to communicate with
      !
      INTRINSIC :: CSHIFT
      INTEGER :: idx_com (0:nprocb-1)
      INTEGER :: i,k,n
      !
      ! get PE identification number
      !
      k = 0
      DO i = dcl%spe, dcl%epe
         IF (dcg(i)% set_a /= seta) CYCLE
         idx_com (k) = i ! dcg(i)% pe
         IF (i == imype) n = k
         k = k + 1
      END DO
      idx_com = CSHIFT (idx_com,n)
    END SUBROUTINE plan_b_sr
    !---------------------------------------------------------------------
  END SUBROUTINE trp_gp_fs_3d
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trp_fs_ls_3d (sign, fs, ls)
    !
    ! transpose
    !   sign= 1 : Fourier space  -> Legendre space
    !   sign=-1 : Fourier space <-  Legendre space

    IMPLICIT NONE

    INTRINSIC :: MOD, SIZE, ASSOCIATED

    ! I/O
    INTEGER              ,INTENT(in)     :: sign ! 1:fs>ls; -1:gs<ls
    REAL(dp)             ,INTENT(inout)  :: fs   (:,:,:)   ! fs
    REAL(dp)             ,INTENT(inout)  :: ls   (:,:,:)   ! ls

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'trp_fs_ls_3d'
    INTEGER              :: i, j, k, l   ! loop indices
    INTEGER              :: imype        ! index of this pe
    INTEGER              :: nlev         ! number of levels
    INTEGER              :: nflat        ! number of latitudes (2*nhgl)
    INTEGER              :: flats(2)     ! first latitude in Fourier space
    INTEGER              :: flate(2)     ! last  latitude in Fourier space
    INTEGER              :: nlm          ! number of m coeff. in Leg. space
    INTEGER ,POINTER     :: intr(:)=>NULL()      ! index array
    INTEGER              :: setb         ! set B
    INTEGER              :: nproca       ! number of PEs in set A
    INTEGER              :: n2mp1        ! total number of coeff. from lgti
    REAL(dp),POINTER     :: bufs(:,:,:)=>NULL()  ! send buffer
    REAL(dp),POINTER     :: bufr(:,:,:)=>NULL()  ! receive buffer
    INTEGER ,ALLOCATABLE :: idx_com (:)  ! PEs to communicate with
    INTEGER              :: nb,nf,n2     ! explicit shapes
    !
    ! invariant local variables
    !
    nlev   = SIZE (fs,2)
    imype  = indx (p_pe, dcg)
    setb   =  dcg(imype)% set_b
    nproca =  dcg(imype)% nproca
    n2mp1  = (dcg(imype)% nm + 1) * 2
    ALLOCATE (idx_com (0:nproca-1))
    CALL plan_a_sr (idx_com)
    SELECT CASE (sign)
    CASE (1)
       !
       ! Fourier -> Legendre space
       !
       i = imype         ! PE index (my)
       DO k = 0, nproca-1
          j = idx_com (           k        ) ! PE index (send)
          l = idx_com (MOD(nproca-k,nproca)) ! PE index (recv)
          CALL alloc_fs_ls_buf (fs_i=i, ls_i=j, zbuf=bufs)
          CALL pack_fs_buf
          IF(i/=j) CALL alloc_fs_ls_buf (fs_i=l, ls_i=i, zbuf=bufr)
          CALL send_recv
          CALL unpack_buf_ls
          IF (ASSOCIATED(bufr)) THEN
             DEALLOCATE(bufr) ; NULLIFY(bufr)
          END IF
       END DO
    CASE (-1)
       !
       ! Legendre -> Fourier
       !
       i = imype         ! PE index (my)
       DO k = 0, nproca-1
          j = idx_com (           k        ) ! PE index (send)
          l = idx_com (MOD(nproca-k,nproca)) ! PE index (recv)
          CALL alloc_fs_ls_buf (fs_i=j, ls_i=i, zbuf=bufs)
          CALL pack_ls_buf
          IF(i/=j) CALL alloc_fs_ls_buf (fs_i=i, ls_i=l, zbuf=bufr)
          CALL send_recv
          CALL unpack_buf_fs
          IF (ASSOCIATED(bufr)) THEN
             DEALLOCATE(bufr) ; NULLIFY(bufr)
          END IF
       END DO
       ! set coefficients, which are not provided by inv. Legendre
       ! transformation, to zero
       fs (n2mp1+1:,:,:) = 0.0_dp
    CASE default
       CALL error_bi ('invalid SIGN parameter (not 1,-1)', substr)
    END SELECT
    !
    ! frt at ECMWF doesnt yet deallocate
    !
    DEALLOCATE (idx_com)

  CONTAINS
    !------------------------------------------------------------------
    SUBROUTINE send_recv
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      IF(i/=j) THEN
         CALL p_sendrecv (bufs,  dcg(j)% pe, &
              bufr, dcg(l)% pe, tag_tr_fs_ls)
         DEALLOCATE (bufs) ; NULLIFY(bufs)
      ELSE
         bufr => bufs
      ENDIF
    END SUBROUTINE send_recv
    !------------------------------------------------------------------
    SUBROUTINE alloc_fs_ls_buf (fs_i, ls_i, zbuf)
      INTEGER       :: fs_i, ls_i  ! indices of Fourier and Lagendre space pe
      REAL(dp),POINTER :: zbuf(:,:,:)
      !
      ! derive bounds and allocate buffer
      !
      nflat =  dcg(fs_i)% nflat
      flats =  dcg(fs_i)% flats
      flate =  dcg(fs_i)% flate
      nlm   =  dcg(ls_i)% nlm
      intr  => dcg(ls_i)% intr
      ALLOCATE (zbuf(2*nlm,nlev,nflat))
      !
      ! determine sizes of buffer and fourier space variable
      ! for explicit shape dummy variables (to enforce vectorization)
      !
      nb = 2*nlm
      nf = SIZE (fs,1)
      n2 = nlev * nflat
    END SUBROUTINE alloc_fs_ls_buf
    !----------------------------------------------------------------
    SUBROUTINE pack_fs_buf
#ifdef EXPLICIT
      CALL pack_fs_buf_ex (fs(1,1,1), bufs(1,1,1), intr, nf,nb,n2)
#else
      bufs(:,:,:) = fs (intr,:,:)
#endif
    END SUBROUTINE pack_fs_buf
    !----------------------------------------------------------------
    !----------------------------------------------------------------
    SUBROUTINE unpack_buf_fs
#ifdef EXPLICIT
      CALL unpack_buf_fs_ex (bufr(1,1,1), fs(1,1,1), intr, nb,nf,n2)
#else
      fs (intr,:,:) = bufr(:,:,:)
#endif
    END SUBROUTINE unpack_buf_fs
    !----------------------------------------------------------------
    SUBROUTINE unpack_buf_ls
#ifdef EXPLICIT
      CALL unpack_buf_ls_ex(ls(1,1,1),bufr(1,1,1),2*nlm*nlev, &
           size(LS,3),nflat,1,flats(1),flate(1),1,nflat/2)
      CALL unpack_buf_ls_ex(ls(1,1,1),bufr(1,1,1),2*nlm*nlev, &
           size(LS,3),nflat,1,flats(2),flate(2),nflat/2+1,nflat)
#else
      ls(:,:,flats(1):flate(1)) = bufr (:,:,         :nflat/2)
      ls(:,:,flats(2):flate(2)) = bufr (:,:,nflat/2+1:       )
#endif
    END SUBROUTINE unpack_buf_ls
    !----------------------------------------------------------------

    !--------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf
#ifdef EXPLICIT
      CALL pack_ls_buf_ex(bufs(1,1,1),ls(1,1,1),2*nlm*nlev, &
           size(LS,3),nflat,1,flats(1),flate(1),1,nflat/2)
      CALL pack_ls_buf_ex(bufs(1,1,1),ls(1,1,1),2*nlm*nlev, &
           size(LS,3),nflat,1,flats(2),flate(2),nflat/2+1,nflat)
#else
      bufs (:,:,         :nflat/2) = ls(:,:,flats(1):flate(1))
      bufs (:,:,nflat/2+1:       ) = ls(:,:,flats(2):flate(2))
#endif
    END SUBROUTINE pack_ls_buf
    !----------------------------------------------------------------
    SUBROUTINE plan_a_sr (idx_com)
      !
      ! set up plan of PEs to communicate with
      !
      INTRINSIC :: CSHIFT
      INTEGER :: idx_com (0:nproca-1)
      INTEGER :: i,k,n
      !
      ! get PE identification number
      !
      k = 0
      DO i = dcl%spe, dcl%epe
         IF (dcg(i)% set_b /= setb) CYCLE
         idx_com (k) = i ! dcg(i)% pe
         IF (i == imype) n = k
         k = k + 1
      END DO
      idx_com = CSHIFT (idx_com,n)
    END SUBROUTINE plan_a_sr
    !------------------------------------------------------------------
  END SUBROUTINE trp_fs_ls_3d
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trp_ls_sp_3d (sign, pls, psp)
    !
    ! transpose
    !   sign= 1 : Legendre space  -> spectral space
    !   sign=-1 : Legendre space <-  spectral space

    IMPLICIT NONE

    INTRINSIC :: MOD, ASSOCIATED, SIZE

    ! I/O
    INTEGER              ,INTENT(in)     :: sign       ! 1:ls>sp; -1:ls<sp
    REAL(dp)             ,INTENT(inout)  :: pls(:,:,:) ! Legendre space
    REAL(dp)             ,INTENT(inout)  :: psp(:,:,:) ! spectral space

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'trp_ls_sp_3d'
    REAL(dp),POINTER    :: ls(:,:,:)=>NULL()   ! legendre
    REAL(dp),POINTER    :: sp(:,:,:)=>NULL()   ! spectral
    REAL(dp),POINTER    :: bufs(:,:,:)=>NULL()   ! send buffer
    REAL(dp),POINTER    :: bufr(:,:,:)=>NULL()   ! receive buffer
    INTEGER             :: seta          ! this set A
    INTEGER             :: k, i, j, l    ! loop indices
    INTEGER             :: imype         ! decomp. table index of this pe
    INTEGER             :: nllevp1       ! number of levels in Legendre space
    INTEGER             :: llevs         ! first level in Legendre space
    INTEGER             :: lleve         ! last level in Legendre space
    INTEGER             :: snsp          ! number of coeff. in sp. space
    INTEGER             :: ssps          ! first coefficients in spec. space
    INTEGER             :: sspe          ! last coefficients in spec. space
    INTEGER             :: ke, nk        ! actual last level, no. of levels
    INTEGER             :: nprocb        ! number of PEs in set A
    INTEGER,ALLOCATABLE :: idx_com (:)   ! PEs to communicate with

    ALLOCATE (ls(dcl%nllevp1, SIZE(pls,2), SIZE(pls,3))); ls(:,:,:) = 0.0_dp
    ALLOCATE (sp(dcl%nlev+1, SIZE(psp,2), SIZE(psp,3))); sp(:,:,:) = 0.0_dp

    imype  = indx (p_pe, dcg)
    seta   = dcg(imype)% set_a
    nprocb = dcg(imype)% nprocb
    ALLOCATE (idx_com (0:nprocb-1))
    CALL plan_b_sr (idx_com)
    i = imype         ! PE index (my)
    SELECT CASE (sign)
    CASE (1)
       !
       ! Legendre space -> spectral space
       !
       ls(:SIZE(pls,1),:,:) = pls(:,:,:)
       DO k = 0, nprocb-1
          j = idx_com (           k        ) ! PE index (send)
          l = idx_com (MOD(nprocb-k,nprocb)) ! PE index (recv)
          CALL alloc_ls_sp_buf (ls_i=i, sp_i=j, zbuf=bufs)
          CALL pack_ls_buf
          IF(i/=j) CALL alloc_ls_sp_buf (ls_i=l, sp_i=i, zbuf=bufr)
          CALL send_recv
          CALL unpack_buf_sp
          IF (ASSOCIATED(bufr)) THEN
             DEALLOCATE(bufr) ; NULLIFY(bufr)
          END IF
       END DO
       psp(:,:,:) = sp(:SIZE(psp,1),:,:)
    CASE (-1)
       !
       ! Legendre space <- spectral space
       !
       sp(:SIZE(psp,1),:,:) = psp(:,:,:)
       DO k = 0, nprocb-1
          j = idx_com (           k        ) ! PE index (send)
          l = idx_com (MOD(nprocb-k,nprocb)) ! PE index (recv)
          CALL alloc_ls_sp_buf (ls_i=j, sp_i=i, zbuf=bufs)
          CALL pack_sp_buf
          IF(i/=j) CALL alloc_ls_sp_buf (ls_i=i, sp_i=l, zbuf=bufr)
          CALL send_recv
          CALL unpack_buf_ls
          IF (ASSOCIATED(bufr)) THEN
             DEALLOCATE(bufr) ; NULLIFY(bufr)
          END IF
       END DO
       pls(:,:,:) = ls(:SIZE(pls,1),:,:)
    CASE default
       CALL error_bi ('invalid SIGN parameter (not 1,-1)', substr)
    END SELECT
    !
    ! frt at ECMWF doesnt yet deallocate
    !
    DEALLOCATE (idx_com)
    DEALLOCATE (ls, sp)

  CONTAINS
    !-------------------------------------------------------------------
    SUBROUTINE alloc_ls_sp_buf (ls_i, sp_i, zbuf)
      INTEGER :: ls_i, sp_i   ! Legendre and spectral space pe indices
      REAL(dp),POINTER :: zbuf(:,:,:)
      !
      ! derive bounds and allocate buffer
      !
      nllevp1 = dcg(ls_i)% nllevp1 ! number of levels in Legendre space
      llevs   = dcg(ls_i)% llevs  ! first level in Legendre space
      lleve   = dcg(ls_i)% lleve  ! last level in Legendre space
      snsp    = dcg(sp_i)% snsp   ! number of coefficients in sp. space
      ssps    = dcg(sp_i)% ssps   ! first coefficients in spectral space
      sspe    = dcg(sp_i)% sspe   ! last coefficients in spectral space

      ALLOCATE (zbuf (nllevp1, 2, snsp))

    END SUBROUTINE alloc_ls_sp_buf
    !------------------------------------------------------------------
    SUBROUTINE pack_ls_buf
      INTRINSIC :: SIZE
      bufs (:SIZE(ls,1),:,:)     = ls (:,:,ssps:sspe)
    END SUBROUTINE pack_ls_buf
    !------------------------------------------------------------------
    SUBROUTINE unpack_buf_sp
      INTRINSIC :: MIN, SIZE
      ke = MIN(SIZE(sp,1), lleve); nk = ke-llevs+1
      IF(nk>0) sp (llevs:ke, :, :) = bufr (:nk,:,:)
    END SUBROUTINE unpack_buf_sp
    !------------------------------------------------------------------
    SUBROUTINE pack_sp_buf
      INTRINSIC :: MAX, MIN, SIZE
      ke = MIN(SIZE(sp,1), lleve); nk = MAX(0, ke-llevs+1)
      bufs (:nk,:,:) = sp (llevs:ke, :, :); bufs(nk+1:,:,:) = 0.
    END SUBROUTINE pack_sp_buf
    !------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls
      INTRINSIC :: SIZE
      ls (:,:,ssps:sspe) = bufr (:SIZE(ls,1),:,:)
    END SUBROUTINE unpack_buf_ls
    !------------------------------------------------------------------
    SUBROUTINE send_recv
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      IF(i/=j) THEN
         CALL p_sendrecv (bufs,  dcg(j)% pe, &
              bufr, dcg(l)% pe, tag_tr_ls_sp)
         DEALLOCATE (bufs) ; NULLIFY(bufs)
      ELSE
         bufr => bufs
      ENDIF
    END SUBROUTINE send_recv
    !-----------------------------------------------------------------
    SUBROUTINE plan_b_sr (idx_com)
      !
      ! set up plan of PEs to communicate with
      !
      INTRINSIC :: CSHIFT
      INTEGER :: idx_com (0:nprocb-1)
      INTEGER :: i,k,n
      !
      ! get PE identification number
      !
      k = 0
      DO i = dcl%spe, dcl%epe
         IF (dcg(i)% set_a /= seta) CYCLE
         idx_com (k) = i ! dcg(i)% pe
         IF (i == imype) n = k
         k = k + 1
      END DO
      idx_com = CSHIFT (idx_com,n)
    END SUBROUTINE plan_b_sr
    !-----------------------------------------------------------------
  END SUBROUTINE trp_ls_sp_3d
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_lti_0(fa, fs, ls)

    ! INVERSE LEGENDRE TRANSFORMATION (0 PARAMETER)
    !
    !
    ! Author: P. Joeckel, MPICH, Jan 2004
    !
    !
    !
    ! antisymmetric:
    !                    T
    !                   ---
    !                   \
    !      FA (m, mu) =  >    LS([m,n]) * P([m,n], mu)
    !                   /
    !                   ---
    !                   n=m, n:even
    !
    ! symmetric:
    !                    T
    !                   ---
    !                   \
    !      FS (m, mu) =  >    LS([m,n]) * P([m,n], mu)
    !                   /
    !                   ---
    !                   n=m, n:odd
    !
    !
    !
    !
    ! Notes:
    !      -  P([m,n], mu) : associated Legendre Polynomial at position mu
    !      -  T            : (triangular) truncation (= MAX(m))
    !

    ! ECHAM5
    USE mo_legendre,      ONLY: leginv

    IMPLICIT NONE

    ! I/O
    REAL(dp), DIMENSION(:,:), TARGET, INTENT(OUT) :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:), TARGET, INTENT(OUT) :: fs  ! symmetric
    REAL(dp), DIMENSION(:),           INTENT(IN)  :: ls

    ! LOCAL
    REAL(dp), DIMENSION(:,:), POINTER  :: fd=>NULL()
    ! LOOP VARAIABLES
    INTEGER                        :: nlnsp  ! number of (m,n) pairs
    INTEGER                        :: nhgl   ! 1/2 number of latitudes
    INTEGER                        :: nlm    ! number of wave numbers m
    ! ECHAM5 SPECIFIC DECOMPOSITION INFORMATION
    INTEGER ,POINTER :: nlmp(:)=>NULL() ! offset to local m columns
    INTEGER ,POINTER :: nlnp(:)=>NULL() ! column length (No. of n's)
    ! LOCAL SCALARS
    LOGICAL, SAVE                  :: ini_flag = .FALSE. ! P(m,n,mu) calculated

    INTEGER :: jm           ! LOOP OVER WAVE NUMBERS m
    INTEGER :: i            ! LOOP FOR 'odd', 'even' n
    INTEGER :: knm          ! INDEX FOR (m,n=m) OR (m,n=m+1)
    INTEGER :: knT          ! INDEX FOR (m,n=T) OR (m,n=T-1)
    INTEGER :: nn, nnh      ! NUMBER OF n FOR THIS m

    !  LOCAL ARRAYS
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE :: pnmit

    !  INTRINSIC FUNCTIONS
    INTRINSIC :: MATMUL, SIZE

    ! INIT
    nlnsp = SIZE(ls)    ! dcl% lnsp
    nlm   = SIZE(fa,1)  ! dcl% nlm
    nhgl  = SIZE(fa,2)  ! dcl% nlat/2
    !

    !-- ECHAM5 specific decomposition information
    nlmp    => dcl% nlmp     ! displacement of the first point of columns
    nlnp    => dcl% nlnp     ! number of points on each column


    IF ( .NOT. ini_flag ) THEN

       ALLOCATE(pnmit(nhgl, nlnsp))

       CALL leginv(pnmit=pnmit)

       ini_flag = .TRUE.

    ENDIF

    DO jm = 1, nlm

       DO i = 1, 2

          knm  = nlmp(jm) + 1 ! offset to local m columns  (spectral coef.)
          nn   = nlnp(jm)     ! column length (No. of n for current m)

          IF (i == 1) THEN
             fd   => fs      ! summation over n for odd n -> symmetric
             nnh  = (nn+1)/2
          ELSE
             fd   => fa      ! summation over n for even n -> antisymmetric
             knm  = knm + 1
             nnh  = nn/2
          END IF

          IF (nnh>0) THEN

             knT = knm + 2*nnh - 2

             fd(jm, :) = MATMUL(pnmit(:,knm:knT:2),ls(knm:knT:2))

             !WRITE(*,*) 'TRF_LTI (pe =',p_pe,'): ',knm,knT,jm,' --> ',i

          ELSE

             fd(jm,:) = 0.0_dp

          END IF

       END DO

    END DO

  END SUBROUTINE trf_lti_0
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_lti_1(fa, fs, ls)

    ! INVERSE LEGENDRE TRANSFORMATION (1 PARAMETER SUCH AS Re,Im OR level)
    !
    !
    ! Author: P. Joeckel, MPICH, Jan 2004
    !
    !

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(OUT) :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(OUT) :: fs  ! symmetric
    REAL(dp), DIMENSION(:,:),           INTENT(IN)  :: ls

    ! LOCAL
    INTEGER :: n1, i1

    n1 = SIZE(ls,1)

    DO i1 = 1, n1
       CALL trf_lti_0(fa(i1,:,:),fs(i1,:,:),ls(i1,:))
    END DO

  END SUBROUTINE trf_lti_1
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_lti_2(fa, fs, ls)

    ! INVERSE LEGENDRE TRANSFORMATION (2 PARAMETER SUCH AS Re,Im AND level)
    !
    !
    ! Author: P. Joeckel, MPICH, Jan 2004
    !
    !

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(OUT) :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(OUT) :: fs  ! symmetric
    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: ls

    ! LOCAL
    INTEGER  :: n1, n2, i1, i2

    n1 = SIZE(ls,1)
    n2 = SIZE(ls,2)

    DO i2 = 1, n2
       DO i1 = 1, n1
          CALL trf_lti_0(fa(i1,i2,:,:),fs(i1,i2,:,:),ls(i1,i2,:))
       END DO
    END DO

  END SUBROUTINE trf_lti_2
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ltd_0(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (0 PARAMETER)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !
    !
    !
    !
    !                    ---
    !                    \
    !      LS ([m, n]) =  >    FA(m, mu) * Q([m,n], mu)   ; n: even
    !                    /
    !                    ---
    !                    mu
    !
    !
    !                    ---
    !                    \
    !      LS ([m, n]) =  >    FS(m, mu) * Q([m,n], mu)   ; n: odd
    !                    /
    !                    ---
    !                    mu
    !

    ! ECHAM5
    USE mo_legendre,      ONLY: legmod

    IMPLICIT NONE

    ! I/O
    REAL(dp), DIMENSION(:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:), TARGET, INTENT(IN)  :: fs  ! symmetric


    ! LOCAL
    REAL(dp), DIMENSION(:,:), POINTER  :: fd=>NULL()
    ! LOOP VARAIABLES
    INTEGER                        :: nlnsp  ! number of (m,n) pairs
    INTEGER                        :: nhgl   ! 1/2 number of latitudes
    INTEGER                        :: nlm    ! number of wave numbers m
    ! ECHAM5 SPECIFIC DECOMPOSITION INFORMATION
    INTEGER ,POINTER :: nlmp(:)=>NULL() ! offset to local m columns
    INTEGER ,POINTER :: nlnp(:)=>NULL() ! column length (No. of n's)
    ! LOCAL SCALARS
    LOGICAL, SAVE                  :: ini_flag = .FALSE. ! P(m,n,mu) calculated

    INTEGER :: jh           ! LOOP OVER HEMISPHERES
    INTEGER :: iu
    INTEGER :: jm           ! LOOP OVER WAVE NUMBERS m
    INTEGER :: knm          ! INDEX FOR (m,n=m) OR (m,n=m+1)
    INTEGER :: knT          ! INDEX FOR (m,n=T) OR (m,n=T-1)

    !  LOCAL ARRAYS
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE :: pnmt

    !  INTRINSIC FUNCTIONS
    INTRINSIC :: MATMUL, MOD, SIZE

    ! INIT
    nlnsp = SIZE(ls)    ! dcl% lnsp
    nlm   = SIZE(fa,1)  ! dcl% nlm
    nhgl  = SIZE(fa,2)  ! dcl% nlat/2
    !

    !-- ECHAM5 specific decomposition information
    nlmp    => dcl% nlmp     ! displacement of the first point of columns
    nlnp    => dcl% nlnp     ! number of points on each column


    IF ( .NOT. ini_flag ) THEN

       ALLOCATE(pnmt(nhgl, nlnsp))

       CALL legmod(pnmt=pnmt)

       ini_flag = .TRUE.

    ENDIF

    DO jh = 1, 2
       iu = 2 - jh

       IF (jh==1) THEN
          fd  => fs (:,:)
       ELSE
          fd  => fa (:,:)
       END IF

       DO jm = 1, nlm

          knm = nlmp(jm) - iu + 2
          knT = nlmp(jm) + nlnp(jm)
          knT = knT - MOD(knT-knm,2)

          IF (knm > nlnsp) CYCLE

          ls(knm:knT:2) = MATMUL(fd(jm,:), pnmt(:,knm:knT:2))

          !WRITE(*,*) 'TRF_LTD (pe =',p_pe,'): ',knm, knT, jm, ' --> ', jh

       END DO

    END DO

  END SUBROUTINE trf_ltd_0
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ltd_1(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (1 PARAMETER SUCH AS Re,Im OR level)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(IN)  :: fs  ! symmetric

    ! LOCAL
    INTEGER :: n1, i1

    n1 = SIZE(ls,1)

    DO i1 = 1, n1
       CALL trf_ltd_0(ls(i1,:),fa(i1,:,:),fs(i1,:,:))
    END DO

  END SUBROUTINE trf_ltd_1
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ltd_2(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (2 PARAMETER SUCH AS Re,Im AND level)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(IN)  :: fs  ! symmetric

    ! LOCAL
    INTEGER :: n1, i1, n2, i2

    n1 = SIZE(ls,1)
    n2 = SIZE(ls,2)

    DO i2 = 1, n2
       DO i1 = 1, n1
          CALL trf_ltd_0(ls(i1,i2,:),fa(i1,i2,:,:),fs(i1,i2,:,:))
       END DO
    END DO

  END SUBROUTINE trf_ltd_2
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trp_fs_as_0(sign, f, fa, fs)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    !   sign= 1 : fourier -> antisymmteric + symmetric
    !   sign=-1 : fourier <- antisymmteric + symmetric

    ! I/O
    INTEGER,                    INTENT(IN)    :: sign
    REAL(dp), DIMENSION(:,:),   INTENT(INOUT) :: f
    REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: fa, fs

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'trp_fs_as_0'
    INTEGER :: nlm  ! NUMBER OF WAVE NUMBERS m
    INTEGER :: nhgl ! 1/2 NUMBER OF LATITUDES

    !     = SIZE(fa,1)   ! =2
    nlm   = SIZE(fa,2)   ! dcl% nlm
    nhgl  = SIZE(fa,3)   ! dcl% nlat / 2

    SELECT CASE (sign)
    CASE (1)
       !
       ! fourier -> antisymmteric + symmetric
       !
       fa(:,:,:) = 0.0_dp
       fs(:,:,:) = 0.0_dp
       !
       fa(1,:,:) = ( f(1:2*nlm-1:2 ,1:nhgl)                    &
                 -   f(1:2*nlm-1:2, 2*nhgl:nhgl+1:-1) ) * 0.5_dp
       fa(2,:,:) = ( f(2:2*nlm:2   ,1:nhgl)                    &
                 -   f(2:2*nlm:2,   2*nhgl:nhgl+1:-1) ) * 0.5_dp
       !
       fs(1,:,:) = ( f(1:2*nlm-1:2 ,1:nhgl)                    &
                 +   f(1:2*nlm-1:2, 2*nhgl:nhgl+1:-1) ) * 0.5_dp
       fs(2,:,:) = ( f(2:2*nlm:2   ,1:nhgl)                    &
                 +   f(2:2*nlm:2,   2*nhgl:nhgl+1:-1) ) * 0.5_dp
       !
    CASE (-1)
       !
       ! fourier <- antisymmteric + symmetric
       !
       f(:,:) = 0.0_dp
       !
       f(1:2*nlm-1:2 ,1:nhgl) = fs(1,:,1:nhgl) + fa(1,:,1:nhgl)
       f(2:2*nlm:2   ,1:nhgl) = fs(2,:,1:nhgl) + fa(2,:,1:nhgl)
       !
       f(1:2*nlm-1:2, 2*nhgl:nhgl+1:-1) = fs(1,:,1:nhgl) - fa(1,:,1:nhgl)
       f(2:2*nlm:2,   2*nhgl:nhgl+1:-1) = fs(2,:,1:nhgl) - fa(2,:,1:nhgl)
       !
       !
    CASE default
       CALL error_bi ('invalid SIGN parameter (not 1,-1)', substr)
    END SELECT

  END SUBROUTINE trp_fs_as_0
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trp_fs_as_1(sign, f, fa, fs)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    !   sign= 1 : fourier -> antisymmteric + symmetric
    !   sign=-1 : fourier <- antisymmteric + symmetric
    !
    ! 1 PARAMETER (SUCH AS level)

    ! I/O
    INTEGER,                      INTENT(IN)    :: sign
    REAL(dp), DIMENSION(:,:,:),   INTENT(INOUT) :: f
    REAL(dp), DIMENSION(:,:,:,:), INTENT(INOUT) :: fa, fs

    ! LOCAL
    INTEGER :: n1, i1

    n1 = SIZE(fa,1)

    DO i1=1,n1
       CALL trp_fs_as_0(sign, f(:,i1,:), fa(i1,:,:,:), fs(i1,:,:,:))
    END DO

  END SUBROUTINE trp_fs_as_1
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ffti(ffs)

    ! ECHAM5
#if defined(ECHAM5)
#ifdef FFT991
    USE mo_fft991, ONLY: fft991cy
#else
    USE mo_fft992, ONLY: fft992
#endif
#endif

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: ffs

    ! LOCAL
    INTEGER                              :: inc, isign
    INTEGER                              :: nlon, nlp2, nlev, nlat
    REAL(dp), DIMENSION(:), ALLOCATABLE  :: zwork

    ! TRANSFORMATION: INVERSE FOURIER
    inc    = 1
    isign  = 1
    nlon   = SIZE(ffs,1) - 2   ! dcl% nlon
    nlp2   = nlon + 2
    nlev   = SIZE(ffs,2)
    nlat   = SIZE(ffs,3)       ! dcl% nflat

    ALLOCATE(zwork(nlp2 * nlev * nlat))

#if defined(ECHAM5)
#ifdef FFT991
    CALL fft991cy(ffs, zwork, inc, nlp2, nlon, nlev*nlat, isign)
#else
    CALL fft992(ffs, inc, nlp2, nlon, nlev*nlat, isign)
#endif
#endif

    DEALLOCATE(zwork)

  END SUBROUTINE trf_ffti
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_fftd(ffs)

    ! ECHAM5
#if defined(ECHAM5)
#ifdef FFT991
    USE mo_fft991, ONLY: fft991cy
#else
    USE mo_fft992, ONLY: fft992
#endif
#endif

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: ffs

    ! LOCAL
    INTEGER                              :: inc, isign
    INTEGER                              :: nlon, nlp2, nlev, nlat
    REAL(dp), DIMENSION(:), ALLOCATABLE  :: zwork

    ! TRANSFORMATION: INVERSE FOURIER
    inc    = 1
    isign  = -1
    nlon   = SIZE(ffs,1) - 2  ! dcl% nlon
    nlp2   = nlon + 2
    nlev   = SIZE(ffs,2)
    nlat   = SIZE(ffs,3)  ! dcl% nflat

    ALLOCATE(zwork(nlp2 * nlev * nlat))

#if defined(ECHAM5)
#ifdef FFT991
    CALL fft991cy(ffs, zwork, inc, nlp2, nlon, nlev*nlat, isign)
#else
    CALL fft992(ffs,inc,nlp2,nlon,nlev*nlat,isign)
#endif
#endif

    DEALLOCATE(zwork)

  END SUBROUTINE trf_fftd
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE sp2gp(sp, gp, lzm0)

    IMPLICIT NONE

    INTRINSIC :: MIN, SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:), INTENT(IN)  :: sp
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT) :: gp
    LOGICAL,                    INTENT(IN)  :: lzm0

    ! LOCAL
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zsp
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zls
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET  :: zfsa
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET  :: zfss
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zlfs
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zffs
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zgp
    !
    INTEGER :: nproma, ngpblks, nlev, nllev, nflev
    INTEGER :: snsp, lnsp, nhgl, nlat, nflat, nlm, nlon

    ! INIT
    nlev    = SIZE(sp,1)
    !2      = SIZE(sp,2)
    snsp    = SIZE(sp,3)
    !
    nproma  = SIZE(gp,1)
    nlev    = SIZE(gp,2)
    ngpblks = SIZE(gp,3)

    nllev = MIN(nlev, dcl%nllev)
    nflev = MIN(nlev, dcl%nflev)
    lnsp  = dcl%lnsp
    nhgl  = dcl%nlat/2
    nlat  = dcl%nlat
    nflat = dcl%nflat
    nlm   = dcl%nlm
    nlon  = dcl%nlon

    ! WORKSPACE
    ALLOCATE(zsp(nlev, 2, snsp))          ! COPY NEEDED FOR INTENT(INOUT)
    !
    ALLOCATE(zls (nllev, 2, lnsp))
    ALLOCATE(zfsa(nllev, 2, nlm, nhgl))
    ALLOCATE(zfss(nllev, 2, nlm, nhgl))
    !
    ALLOCATE(zlfs(nlm*2, nflev, nlat))
    ALLOCATE(zffs(nlon+2, nflev, nflat))
    !
    ALLOCATE(zgp(nproma, nlev, ngpblks))  ! COPY NEEDED FOR INTENT(INOUT)

    zsp(:,:,:) = sp(:,:,:)
    zgp(:,:,:) = 0.0_dp
    gp(:,:,:)  = 0.0_dp

    ! TRANSPOSITION: SPECTRAL SPACE -> LEGENDRE SPACE
    zls(:,:,:) = 0.0_dp
    CALL trp_ls_sp_3d(-1, zls, zsp)

    IF (lzm0) THEN
       IF (dcl% nlnm0 > 0) zls(:,1,1) = 0.0_dp
    END IF

    ! TRANSFORMATION: INVERSE LEGENDRE
    zfsa(:,:,:,:) = 0.0_dp
    zfss(:,:,:,:) = 0.0_dp
    CALL trf_lti(zfsa, zfss, zls)
    !
    ! ANTISYMMETRIC + SYMMTERIC -> FOURIER
    zlfs(:,:,:) = 0.0_dp
    CALL trp_fs_as(-1, zlfs, zfsa, zfss)

    ! TRANSPOSITION: LEGENDRE SPACE -> FOURIER SPACE
    zffs(:,:,:) = 0.0_dp
    CALL trp_fs_ls_3d(-1, zffs, zlfs)

    ! TRANSFORMATION: INVERSE FOURIER
    CALL trf_ffti(zffs)

    ! TRANSPOSITION: FOURIER SPACE -> GRIDPOINT SPACE
    CALL trp_gp_fs_3d(-1, zgp, zffs)
    gp(:,:,:) = zgp(:,:,:)

    ! FREE MEMORY
    DEALLOCATE(zsp)
    !
    DEALLOCATE(zls)
    DEALLOCATE(zfsa)
    DEALLOCATE(zfss)
    !
    DEALLOCATE(zlfs)
    DEALLOCATE(zffs)
    !
    DEALLOCATE(zgp)

  END SUBROUTINE sp2gp
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE gp2sp(gp, sp, lzm0)

    IMPLICIT NONE

    INTRINSIC :: MIN, SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:), INTENT(IN)  :: gp
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT) :: sp
    LOGICAL,                    INTENT(IN)  :: lzm0

    ! LOCAL
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zsp
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zls
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET  :: zfsa
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET  :: zfss
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zlfs
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zffs
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE          :: zgp
    !
    INTEGER :: nproma, ngpblks, nlev, nllev, nflev
    INTEGER :: snsp, lnsp, nhgl, nlat, nflat, nlm, nlon

    ! INIT
    nlev    = SIZE(sp,1)
    !2      = SIZE(sp,2)
    snsp    = SIZE(sp,3)
    !
    nproma  = SIZE(gp,1)
    nlev    = SIZE(gp,2)
    ngpblks = SIZE(gp,3)

    nllev = MIN(nlev, dcl%nllev)
    nflev = MIN(nlev, dcl%nflev)
    lnsp  = dcl%lnsp
    nhgl  = dcl%nlat/2
    nlat  = dcl%nlat
    nflat = dcl%nflat
    nlm   = dcl%nlm
    nlon  = dcl%nlon

    ! WORKSPACE
    ALLOCATE(zsp(nlev, 2, snsp))          ! COPY NEEDED FOR INTENT(INOUT)
    !
    ALLOCATE(zls (nllev, 2, lnsp))
    ALLOCATE(zfsa(nllev, 2, nlm, nhgl))
    ALLOCATE(zfss(nllev, 2, nlm, nhgl))
    !
    ALLOCATE(zlfs(nlm*2, nflev, nlat))
    ALLOCATE(zffs(nlon+2, nflev, nflat))
    !
    ALLOCATE(zgp(nproma, nlev, ngpblks))  ! COPY NEEDED FOR INTENT(INOUT)

    zgp(:,:,:) = gp(:,:,:)
    zsp(:,:,:) = 0.0_dp
    sp(:,:,:)  = 0.0_dp

    ! TRANSPOSITION: GRIDPOINT SPACE -> FOURIER SPACE
    zffs = 0.0_dp
    CALL trp_gp_fs_3d(1, zgp, zffs)

    ! TRANSFORMATION: DIRECT FOURIER
    CALL trf_fftd(zffs)

    ! TRANSPOSITION: FOURIER SPACE -> LEGENDRE SPACE
    zlfs(:,:,:) = 0.0_dp
    CALL trp_fs_ls_3d(1, zffs, zlfs)

    !  FOURIER -> ANTISYMMETRIC + SYMMTERIC
    zfsa(:,:,:,:) = 0.0_dp
    zfss(:,:,:,:) = 0.0_dp
    CALL trp_fs_as(1, zlfs, zfsa, zfss)
    !
    ! TRANSFORMATION: DIRECT LEGENDRE
    zls(:,:,:) = 0.0_dp
    CALL trf_ltd(zls, zfsa, zfss)

    IF (lzm0) THEN
       IF (dcl% nlnm0 > 0) zls(:,1,1) = 0.0_dp
    END IF

    ! TRANSPOSITION: LEGENDRE SPACE -> SPECTRAL SPACE
    CALL trp_ls_sp_3d(1, zls, zsp)

    ! COPY RESULT
    sp(:,:,:) = zsp(:,:,:)

    ! FREE MEMORY
    DEALLOCATE(zsp)
    !
    DEALLOCATE(zls)
    DEALLOCATE(zfsa)
    DEALLOCATE(zfss)
    !
    DEALLOCATE(zlfs)
    DEALLOCATE(zffs)
    !
    DEALLOCATE(zgp)

  END SUBROUTINE gp2sp
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ls_dz2uv_0(d, z, u, v)

    ! CALCULATE (U,V)*cos(lat) BY SPECTRAL TRANSFORMATION
    ! (O PARAMETER)

    ! ECHAM5
    USE mo_legendre,      ONLY: leginv

    IMPLICIT NONE

    ! I/O
    REAL(dp), DIMENSION(:,:),  INTENT(IN)  :: d  ! divergence
    REAL(dp), DIMENSION(:,:),  INTENT(IN)  :: z  ! vorticity
    REAL(dp), DIMENSION(:,:),  INTENT(OUT) :: u  ! u * cos(lat)
    REAL(dp), DIMENSION(:,:),  INTENT(OUT) :: v  ! v * cos(lat)

    ! LOCAL
    ! LOOP VARAIABLES
    INTEGER                        :: nlnsp  ! number of (m,n) pairs
    INTEGER                        :: nhgl   ! 1/2 number of latitudes
    INTEGER                        :: nlm    ! number of wave numbers m
    ! ECHAM5 SPECIFIC DECOMPOSITION INFORMATION
    INTEGER ,POINTER :: nlmp(:)=>NULL() ! offset to local m columns
    INTEGER ,POINTER :: nlnp(:)=>NULL() ! column length (No. of n's)
    ! LOCAL SCALARS
    LOGICAL, SAVE                  :: ini_flag = .FALSE. ! P(m,n,mu) calculated

    INTEGER :: jm           ! LOOP OVER WAVE NUMBERS m
    INTEGER :: i            ! LOOP FOR 'symmetric', 'asymmetric'
    INTEGER :: knm          ! INDEX FOR (m,n=m) OR (m,n=m+1)
    INTEGER :: knT          ! INDEX FOR (m,n=T) OR (m,n=T-1)
    INTEGER :: nn, nnh      ! NUMBER OF n FOR THIS m

    ! LOCAL ARRAYS
    !   Pnm * a * m/(n(n+1))
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE     :: pnmi_uv
    !   d(Pnm)/d(mu) * a * (1.-mu**2)/(n(n+1))
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE     :: anmi_uv
    !
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE         :: fas_d,  fas_z
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE         :: fas_dd, fas_dz
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE, TARGET :: fau, fsu, fav, fsv

    !  INTRINSIC FUNCTIONS
    INTRINSIC :: MATMUL, SIZE, TRANSPOSE

    ! INIT
    ! 2   = SIZE(d,1)     ! Re,Im
    nlnsp = SIZE(d,2)     ! dcl% lnsp

    !-- ECHAM5 specific decomposition information
    nlm   =  dcl%nlm
    nhgl  =  dcl%nlat/2
    nlmp  => dcl% nlmp     ! displacement of the first point of columns
    nlnp  => dcl% nlnp     ! number of points on each column

    !                              1 2
    !                      1 2     S A
    !                     Re,Im     |
    !                       |       |
    ALLOCATE( fas_d ( nhgl, 2, nlm, 2))
    ALLOCATE( fas_z ( nhgl, 2, nlm, 2))
    ALLOCATE( fas_dd( nhgl, 2, nlm, 2))
    ALLOCATE( fas_dz( nhgl, 2, nlm, 2))

    ALLOCATE(fau(2, nlm, nhgl))
    ALLOCATE(fsu(2, nlm, nhgl))
    ALLOCATE(fav(2, nlm, nhgl))
    ALLOCATE(fsv(2, nlm, nhgl))

    IF ( .NOT. ini_flag ) THEN
       ALLOCATE(pnmi_uv(nhgl, nlnsp))
       ALLOCATE(anmi_uv(nhgl, nlnsp))
       CALL leginv(pnmiuvt=pnmi_uv ,anmiuvt=anmi_uv)
       ini_flag = .TRUE.
    ENDIF

    DO jm = 1, nlm

       DO i = 1, 2

          knm  = nlmp(jm) + 1 ! offset to local m columns  (spectral coef.)
          nn   = nlnp(jm)     ! column length (No. of n for current m)

          IF (i == 1) THEN
             nnh  = (nn+1)/2
          ELSE
             knm  = knm + 1
             nnh  = nn/2
          END IF

          IF (nnh>0) THEN

             knT = knm + 2*nnh - 2

             !     1:nhgl
             !      |     S A
             !      |Re,Im |                                Re,Im
             !      | |  m |                                  |
             fas_d (:,1,jm,i) = MATMUL(pnmi_uv(:,knm:knT:2),d(1,knm:knT:2))
             fas_d (:,2,jm,i) = MATMUL(pnmi_uv(:,knm:knT:2),d(2,knm:knT:2))

             fas_z (:,1,jm,i) = MATMUL(pnmi_uv(:,knm:knT:2),z(1,knm:knT:2))
             fas_z (:,2,jm,i) = MATMUL(pnmi_uv(:,knm:knT:2),z(2,knm:knT:2))

             fas_dd(:,1,jm,i) = MATMUL(anmi_uv(:,knm:knT:2),d(1,knm:knT:2))
             fas_dd(:,2,jm,i) = MATMUL(anmi_uv(:,knm:knT:2),d(2,knm:knT:2))

             fas_dz(:,1,jm,i) = MATMUL(anmi_uv(:,knm:knT:2),z(1,knm:knT:2))
             fas_dz(:,2,jm,i) = MATMUL(anmi_uv(:,knm:knT:2),z(2,knm:knT:2))

          ELSE

             fas_d (:,:,jm,i) = 0.0_dp
             fas_z (:,:,jm,i) = 0.0_dp
             fas_dd(:,:,jm,i) = 0.0_dp
             fas_dz(:,:,jm,i) = 0.0_dp

          END IF

       END DO

    END DO

    ! COMBINE DIVERGENT AND ROTATIONAL PARTS OF U AND V
    !
    !  1 2                        1 2   1 2            1 2      1 2
    ! Re,Im                      Re,Im  S A           Re,Im     S A
    !   |   1:nhgl          1:nhgl |     |       1:nhgl |        |
    !   | m |                      |  |  m  |           |  |  m  |
    fsu(1,:,:) = TRANSPOSE( fas_dz(:, 1, :, 2) + fas_d (:, 2, :, 1) )
    fsu(2,:,:) = TRANSPOSE( fas_dz(:, 2, :, 2) - fas_d (:, 1, :, 1) )
    fau(1,:,:) = TRANSPOSE( fas_dz(:, 1, :, 1) + fas_d (:, 2, :, 2) )
    fau(2,:,:) = TRANSPOSE( fas_dz(:, 2, :, 1) - fas_d (:, 1, :, 2) )
    !
    fsv(1,:,:) = TRANSPOSE( fas_z (:, 2, :, 1) - fas_dd(:, 1, :, 2) )
    fsv(2,:,:) = TRANSPOSE(-fas_z (:, 1, :, 1) - fas_dd(:, 2, :, 2) )
    fav(1,:,:) = TRANSPOSE( fas_z (:, 2, :, 2) - fas_dd(:, 1, :, 1) )
    fav(2,:,:) = TRANSPOSE(-fas_z (:, 1, :, 2) - fas_dd(:, 2, :, 1) )


    ! DIRECT LEGENDRE TRANSFORMATION
    CALL trf_ltd_1(u, fau, fsu)
    CALL trf_ltd_1(v, fav, fsv)

    ! FREE MEMORY
    DEALLOCATE( fas_d )
    DEALLOCATE( fas_z )
    DEALLOCATE( fas_dd)
    DEALLOCATE( fas_dz)
    !
    DEALLOCATE(fau)
    DEALLOCATE(fsu)
    DEALLOCATE(fav)
    DEALLOCATE(fsv)

  END SUBROUTINE trf_ls_dz2uv_0
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ls_dz2uv_1(d, z, u, v)

    ! CALCULATE (U,V)*cos(lat) BY SPECTRAL TRANSFORMATION
    ! (1 PARAMETER, SUCH AS level)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    !             Re,Im lnsp
    !                 | |
    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: d  ! divergence
    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: z  ! vorticity
    REAL(dp), DIMENSION(:,:,:),           INTENT(OUT) :: u  ! u * cos(lat)
    REAL(dp), DIMENSION(:,:,:),           INTENT(OUT) :: v  ! v * cos(lat)

    ! LOCAL
    INTEGER :: n1, i1

    n1 = SIZE(d,1)

    DO i1 = 1, n1
       CALL trf_ls_dz2uv_0(d(i1,:,:), z(i1,:,:), u(i1,:,:), v(i1,:,:))
    END DO

  END SUBROUTINE trf_ls_dz2uv_1
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE trp_gpdc_gpgl_4d(sign, lf, gf, method)

    ! ECHAM5/MESSy
#ifdef MPIPERF_BAD
    USE messy_main_mpi_bi, ONLY: gather_field, scatter_gp
#else
    USE messy_main_mpi_bi, ONLY: gather_field
#endif
#ifdef _MPI_ALLGATHER
    USE messy_main_mpi_bi, ONLY: allgather_field
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT, SIZE
#ifdef MPIPERF_BAD
    INTRINSIC :: ABS, REAL, SQRT
#endif

    ! I/O
    INTEGER,                  INTENT(IN)           :: sign
    REAL(dp), DIMENSION(:,:,:,:), POINTER          :: lf    ! decomposed field
    REAL(dp), DIMENSION(:,:,:,:), POINTER          :: gf    ! global field
    INTEGER,                  INTENT(IN), OPTIONAL :: method

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr = 'trp_gpdc_gpgl_4d'
    INTEGER                               :: i
    INTEGER                               :: zmethod
    INTEGER                               :: kproma, jp
    INTEGER                               :: jrow
#ifdef MPIPERF_BAD
    REAL(dp)                              :: qave, qstd  ! quotients
    INTEGER                               :: pe
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: zgf=>NULL()         ! result
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: zgft=>NULL()        ! for transfer
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: zgfh=>NULL()        ! for std
#else
    REAL(dp), DIMENSION(:,:,:), POINTER :: ltmp   => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: gtmp   => NULL()
#endif
    INTEGER, DIMENSION(4) :: gls

    SELECT CASE (sign)
    CASE (1)
       !
       ! decomposed field -> global field
       !
       ! INIT
       IF (ASSOCIATED(gf)) THEN
          DEALLOCATE(gf)
          NULLIFY(gf)
       END IF
       !
#ifndef _MPI_ALLGATHER
       ! GATHER INFO ON IO-CPU
       IF (p_pe == p_io) THEN
          ALLOCATE(gf(dcl%nlon,SIZE(lf, 2),SIZE(lf,3),dcl%nlat))
       ELSE
          NULLIFY(gf)
       END IF
       gls(:) = (/ dcl%nlon, SIZE(lf, 2), SIZE(lf,3), dcl%nlat /)
       CALL gather_field(gf, gls, lf)
       !
       ! DISTRIBUTE GLOBAL FIELD TO ALL CPUs
       IF (p_pe /= p_io) THEN
          ALLOCATE(gf(dcl%nlon,SIZE(lf, 2),SIZE(lf,3),dcl%nlat))
       END IF
       CALL p_bcast(gf, p_io)
#else
       ALLOCATE(gf(dcl%nlon,SIZE(lf, 2),SIZE(lf,3),dcl%nlat))
       CALL allgather_field(gf,lf)
#endif
       !
    CASE (-1)
       !
       ! decomposed field <- global field
       !
       IF (.NOT.PRESENT(method)) THEN
          zmethod = M_LOC          ! default
       ELSE
          zmethod = method
       END IF
       !
       IF (zmethod == M_LOC) THEN
          ! PICK OUT SUBSET ON LOCAL CPU
          CALL set_ilon_ilat
          lf(:,:,:,:) = 0.0_dp
          DO jrow = 1, dcl%ngpblks
             IF ( jrow == dcl%ngpblks ) THEN
                kproma = dcl%npromz
             ELSE
                kproma = dcl%nproma
             END IF
             DO jp = 1, kproma
                lf(jp,:,:,jrow) = gf(ilon(jp, jrow),:,:,ilat(jp, jrow))
             END DO
          END DO
          RETURN
       END IF
       !
#ifdef MPIPERF_BAD
       ! I/O-PE
       IF (p_pe == p_io) THEN
          ALLOCATE(zgft(SIZE(gf,1),SIZE(gf,2),SIZE(gf,3),SIZE(gf,4)))
          ALLOCATE(zgf (SIZE(gf,1),SIZE(gf,2),SIZE(gf,3),SIZE(gf,4)))
          ALLOCATE(zgfh(SIZE(gf,1),SIZE(gf,2),SIZE(gf,3),SIZE(gf,4)))
          SELECT CASE (zmethod)
          CASE (M_SUM,M_AVE)
             zgf(:,:,:,:)  = gf(:,:,:,:)
          CASE (M_STD)
             zgf (:,:,:,:) = gf(:,:,:,:)
             zgfh(:,:,:,:) = gf(:,:,:,:)**2
          CASE (M_LOC)
             ! DO NOTHING ; SEE ABOVE
          CASE DEFAULT
             CALL error_bi('invalid METHOD parameter', substr)
          END SELECT
       ELSE
          NULLIFY(zgft)
          NULLIFY(zgf)
          NULLIFY(zgfh)
       END IF

       ! LOOP OVER ALL NON-I/O PEs
       IF (p_pe /= p_io) THEN
          CALL p_send (gf, p_io, tag_tr_gpdc_gpgl)
       ELSE ! p_pe = p_io
          DO i=1, p_nprocs
             pe = dcg(i)%pe
             IF (pe /= p_pe) THEN
                CALL p_recv( zgft, pe, tag_tr_gpdc_gpgl)
                SELECT CASE (zmethod)
                CASE (M_SUM,M_AVE)
                   zgf(:,:,:,:)  = zgf(:,:,:,:)  + zgft(:,:,:,:)
                CASE (M_STD)
                   zgf(:,:,:,:)  = zgf(:,:,:,:)  + zgft(:,:,:,:)
                   zgfh(:,:,:,:) = zgfh(:,:,:,:) + zgft(:,:,:,:)**2
                CASE (M_LOC)
                   ! DO NOTHING ; SEE ABOVE
                CASE DEFAULT
                   CALL error_bi('invalid METHOD parameter', substr)
                END SELECT
             END IF
          END DO
       END IF
       !
       ! RESULT
       IF (p_pe == p_io) THEN
          SELECT CASE (zmethod)
          CASE (M_SUM)
          CASE (M_AVE)
             IF (p_nprocs > 0) THEN
                qave = 1.0_dp/REAL(p_nprocs,dp)
             ELSE
                qave = 1.0_dp
             END IF
             zgf(:,:,:,:) = zgf(:,:,:,:)*qave
          CASE (M_STD)
             IF (p_nprocs > 0) THEN
                qave = 1.0_dp/REAL(p_nprocs,dp)
             ELSE
                qave = 1.0_dp
             END IF
             IF ((p_nprocs-1) > 0) THEN
                qstd = 1.0_dp/REAL(p_nprocs-1,dp)
             ELSE
                qstd = 0.0_dp
             END IF
             !
             zgf(:,:,:,:) = SQRT( ABS( (zgfh(:,:,:,:) - &
             REAL(p_nprocs,dp)*(zgf(:,:,:,:)*qave)**2)*qstd ) )
          CASE (M_LOC)
             !
          CASE DEFAULT
             CALL error_bi('invalid METHOD parameter', substr)
          END SELECT
       END IF
       !
       ! SCATTER RESULT INTO DECOMPOSED FIELDS
       CALL scatter_gp(zgf, lf, dcg)
       !
       ! CLEAN MEMORY
       IF (ASSOCIATED(zgft)) THEN
          DEALLOCATE(zgft) ; NULLIFY(zgft)
       END IF
       IF (ASSOCIATED(zgf)) THEN
          DEALLOCATE(zgf) ; NULLIFY(zgf)
       END IF
       IF (ASSOCIATED(zgfh)) THEN
          DEALLOCATE(zgfh) ; NULLIFY(zghf)
       END IF
       !
#else
       DO i=1, SIZE(gf,3)
          gtmp => gf(:,:,i,:)
          ltmp => lf(:,:,i,:)
          CALL trp_gpdc_gpgl_3d(sign, ltmp, gtmp, method)
       END DO
#endif

    CASE default
       CALL error_bi ('invalid SIGN parameter (not 1,-1)', substr)
    END SELECT

  END SUBROUTINE trp_gpdc_gpgl_4d
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE trp_gpdc_gpgl_3d(sign, lf, gf, method)

    ! ECHAM5/MESSy
#ifdef MPIPERF_BAD
    USE messy_main_mpi_bi, ONLY: gather_field, scatter_gp
#else
    USE messy_main_mpi_bi, ONLY: gather_field, p_sum
#endif
#ifdef _MPI_ALLGATHER
    USE messy_main_mpi_bi, ONLY: allgather_field
#endif

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, PRESENT, REAL, SIZE, SQRT

    ! I/O
    INTEGER,                INTENT(IN)           :: sign
    REAL(dp), DIMENSION(:,:,:), POINTER          :: lf    ! decomposed field
    REAL(dp), DIMENSION(:,:,:), POINTER          :: gf    ! global field
    INTEGER,                INTENT(IN), OPTIONAL :: method

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr = 'trp_gpdc_gpgl_3d'
    REAL(dp)                            :: qave, qstd  ! quotients
    INTEGER                             :: zmethod
    INTEGER                             :: kproma, jp
    INTEGER                             :: jrow
    REAL(dp), DIMENSION(:,:,:), POINTER :: zgf  => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: zgft => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: zgfh => NULL()
#ifndef MPIPERF_BAD
    REAL(dp), DIMENSION(:,:,:), POINTER :: zlf   => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: zlfh  => NULL()
#else
    INTEGER                             :: pe, i
#endif
    INTEGER, DIMENSION(3) :: gls

    SELECT CASE (sign)
    CASE (1)
       !
       ! decomposed field -> global field
       !
       ! INIT
       IF (ASSOCIATED(gf)) THEN
          DEALLOCATE(gf)
          NULLIFY(gf)
       END IF
       !
#ifndef _MPI_ALLGATHER
       ! GATHER INFO ON IO-CPU
       IF (p_pe == p_io) THEN
          ALLOCATE(gf(dcl%nlon,SIZE(lf, 2),dcl%nlat))
       ELSE
          NULLIFY(gf)
       END IF
       gls(:) = (/ dcl%nlon, SIZE(lf, 2), dcl%nlat /)
       CALL gather_field(gf, gls, lf)
       !
       ! DISTRIBUTE GLOBAL FIELD TO ALL CPUs
       IF (p_pe /= p_io) THEN
          ALLOCATE(gf(dcl%nlon,SIZE(lf, 2),dcl%nlat))
       END IF
       CALL p_bcast(gf, p_io)
#else
       ALLOCATE(gf(dcl%nlon,SIZE(lf, 2),dcl%nlat))
       CALL allgather_field(gf, lf)
#endif
       !
    CASE (-1)
       !
       ! decomposed field <- global field
       !
       IF (.NOT.PRESENT(method)) THEN
          zmethod = M_LOC          ! default
       ELSE
          zmethod = method
       END IF
       !
       IF (zmethod == M_LOC) THEN
          ! PICK OUT SUBSET ON LOCAL CPU
          CALL set_ilon_ilat
          lf(:,:,:) = 0.0_dp
          DO jrow = 1, dcl%ngpblks
             IF ( jrow == dcl%ngpblks ) THEN
                kproma = dcl%npromz
             ELSE
                kproma = dcl%nproma
             END IF
             DO jp = 1, kproma
                lf(jp,:,jrow) = gf(ilon(jp, jrow),:,ilat(jp, jrow))
             END DO
          END DO
          RETURN
       END IF
       !
#ifdef MPIPERF_BAD
       ! I/O-PE
       IF (p_pe == p_io) THEN
          ALLOCATE(zgft(SIZE(gf,1),SIZE(gf,2),SIZE(gf,3)))
          ALLOCATE(zgf (SIZE(gf,1),SIZE(gf,2),SIZE(gf,3)))
          ALLOCATE(zgfh(SIZE(gf,1),SIZE(gf,2),SIZE(gf,3)))
          SELECT CASE (zmethod)
          CASE (M_SUM,M_AVE)
             zgf(:,:,:)  = gf(:,:,:)
          CASE (M_STD)
             zgf (:,:,:) = gf(:,:,:)
             zgfh(:,:,:) = gf(:,:,:)**2
          CASE (M_LOC)
             ! DO NOTHING ; SEE ABOVE
          CASE DEFAULT
             CALL error_bi('invalid METHOD parameter', substr)
          END SELECT
       ELSE
          NULLIFY(zgft)
          NULLIFY(zgf)
          NULLIFY(zgfh)
       END IF

       ! LOOP OVER ALL NON-I/O PEs
       IF (p_pe /= p_io) THEN
          CALL p_send (gf, p_io, tag_tr_gpdc_gpgl)
       ELSE ! p_pe = p_io
          DO i=1, p_nprocs
             pe = dcg(i)%pe
             IF (pe /= p_pe) THEN
                CALL p_recv( zgft, pe, tag_tr_gpdc_gpgl)
                SELECT CASE (zmethod)
                CASE (M_SUM,M_AVE)
                   zgf(:,:,:)  = zgf(:,:,:)  + zgft(:,:,:)
                CASE (M_STD)
                   zgf(:,:,:)  = zgf(:,:,:)  + zgft(:,:,:)
                   zgfh(:,:,:) = zgfh(:,:,:) + zgft(:,:,:)**2
                CASE (M_LOC)
                   ! DO NOTHING ; SEE ABOVE
                CASE DEFAULT
                   CALL error_bi('invalid METHOD parameter', substr)
                END SELECT
             END IF
          END DO
       END IF
       !
       ! RESULT
       IF (p_pe == p_io) THEN
          SELECT CASE (zmethod)
          CASE (M_SUM)
          CASE (M_AVE)
             IF (p_nprocs > 0) THEN
                qave = 1.0_dp/REAL(p_nprocs,dp)
             ELSE
                qave = 1.0_dp
             END IF
             zgf(:,:,:) = zgf(:,:,:)*qave
          CASE (M_STD)
             IF (p_nprocs > 0) THEN
                qave = 1.0_dp/REAL(p_nprocs,dp)
             ELSE
                qave = 1.0_dp
             END IF
             IF ((p_nprocs-1) > 0) THEN
                qstd = 1.0_dp/REAL(p_nprocs-1,dp)
             ELSE
                qstd = 0.0_dp
             END IF
             !
             zgf(:,:,:) = SQRT( ABS( (zgfh(:,:,:) - &
               REAL(p_nprocs,dp)*(zgf(:,:,:)*qave)**2)*qstd ) )
          CASE (M_LOC)
             !
          CASE DEFAULT
             CALL error_bi('invalid METHOD parameter', substr)
          END SELECT
       END IF
       !
       ! SCATTER RESULT INTO DECOMPOSED FIELDS
       CALL scatter_gp(zgf, lf, dcg)
       !
       ! CLEAN MEMORY
       IF (ASSOCIATED(zgft)) THEN
          DEALLOCATE(zgft) ; NULLIFY(zgft)
       END IF
       IF (ASSOCIATED(zgf))  THEN
          DEALLOCATE(zgf) ; NULLIFY(zgf)
       END IF
       IF (ASSOCIATED(zgfh)) THEN
          DEALLOCATE(zgfh) ; NULLIFY(zgfh)
       END IF
       !
#else
       ALLOCATE(zgf(SIZE(gf,1),SIZE(gf,2),SIZE(gf,3)))
       zgf(:,:,:) = p_sum(gf)

       ! PICK OUT SUBSET ON LOCAL CPU
       CALL set_ilon_ilat
       ALLOCATE(zlf(SIZE(lf,1),SIZE(lf,2),SIZE(lf,3)))
       zlf(:,:,:) = 0.0_dp
       DO jrow = 1, dcl%ngpblks
          IF ( jrow == dcl%ngpblks ) THEN
             kproma = dcl%npromz
          ELSE
             kproma = dcl%nproma
          END IF
          DO jp = 1, kproma
             zlf(jp,:,jrow) = zgf(ilon(jp, jrow),:,ilat(jp, jrow))
          END DO
       END DO

       IF (zmethod == M_STD) THEN
          ALLOCATE(zgfh(SIZE(gf,1),SIZE(gf,2),SIZE(gf,3)))
          ALLOCATE(zgft(SIZE(gf,1),SIZE(gf,2),SIZE(gf,3)))
          zgft(:,:,:) = gf(:,:,:)**2
          zgfh(:,:,:) = p_sum(zgft)
          ALLOCATE(zlfh(SIZE(lf,1),SIZE(lf,2),SIZE(lf,3)))
          zlfh(:,:,:) = 0.0_dp
          DO jrow = 1, dcl%ngpblks
             IF ( jrow == dcl%ngpblks ) THEN
                kproma = dcl%npromz
             ELSE
                kproma = dcl%nproma
             END IF
             DO jp = 1, kproma
                zlfh(jp,:,jrow) = zgfh(ilon(jp, jrow),:,ilat(jp, jrow))
             END DO
          END DO
       END IF

       SELECT CASE (zmethod)
       CASE (M_SUM)
          lf(:,:,:) = zlf(:,:,:)
       CASE (M_AVE)
          IF (p_nprocs > 0) THEN
             qave = 1.0_dp/REAL(p_nprocs,dp)
          ELSE
             qave = 1.0_dp
          END IF
          lf(:,:,:) = zlf(:,:,:)*qave
       CASE (M_STD)
          IF (p_nprocs > 0) THEN
             qave = 1.0_dp/REAL(p_nprocs,dp)
          ELSE
             qave = 1.0_dp
          END IF
          IF ((p_nprocs-1) > 0) THEN
             qstd = 1.0_dp/REAL(p_nprocs-1,dp)
          ELSE
             qstd = 0.0_dp
          END IF
          !
          zlf(:,:,:) = SQRT( ABS( (zlfh(:,:,:) - &
               REAL(p_nprocs,dp)*(zlf(:,:,:)*qave)**2)*qstd ) )
       CASE (M_LOC)
          !
       CASE DEFAULT
          CALL error_bi('invalid METHOD parameter', substr)
       END SELECT

       ! CLEAN
       IF (ASSOCIATED(zgf))  THEN
          DEALLOCATE(zgf) ; NULLIFY(zgf)
       END IF
       IF (ASSOCIATED(zgfh)) THEN
          DEALLOCATE(zgfh) ; NULLIFY(zgfh)
       END IF
       IF (ASSOCIATED(zgft)) THEN
          DEALLOCATE(zgft) ; NULLIFY(zgft)
       END IF
       IF (ASSOCIATED(zlf))  THEN
          DEALLOCATE(zlf) ; NULLIFY(zlf)
       END IF
       IF (ASSOCIATED(zlfh)) THEN
          DEALLOCATE(zlfh) ; NULLIFY(zlfh)
       END IF
#endif

    CASE default
       CALL error_bi ('invalid SIGN parameter (not 1,-1)', substr)
    END SELECT

  END SUBROUTINE trp_gpdc_gpgl_3d
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE trp_gpdc_gpgl_2d(sign, lf, gf, method)

    ! ECHAM5/MESSy
#ifdef MPIPERF_BAD
    USE messy_main_mpi_bi, ONLY: gather_field, scatter_gp
#else
    USE messy_main_mpi_bi, ONLY: gather_field, p_sum
#endif
#ifdef _MPI_ALLGATHER
    USE messy_main_mpi_bi, ONLY: allgather_field
#endif

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, PRESENT, REAL, SIZE, SQRT

    ! I/O
    INTEGER,              INTENT(IN)           :: sign
    REAL(dp), DIMENSION(:,:), POINTER          :: lf    ! decomposed field
    REAL(dp), DIMENSION(:,:), POINTER          :: gf    ! global field
    INTEGER,              INTENT(IN), OPTIONAL :: method

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'trp_gpdc_gpgl_2d'
    REAL(dp)                          :: qave, qstd  ! quotients
    INTEGER                           :: zmethod
    INTEGER                           :: kproma, jp
    INTEGER                           :: jrow
    REAL(dp), DIMENSION(:,:), POINTER :: zgf  => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: zgft => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: zgfh => NULL()
#ifndef MPIPERF_BAD
    REAL(dp), DIMENSION(:,:), POINTER :: zlf   => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: zlfh  => NULL()
#else
    INTEGER                           :: pe, i
#endif

    SELECT CASE (sign)
    CASE (1)
       !
       ! decomposed field -> global field
       !
       ! INIT
       IF (ASSOCIATED(gf)) THEN
          DEALLOCATE(gf)
          NULLIFY(gf)
       END IF
       !
#ifndef _MPI_ALLGATHER
       ! GATHER INFO ON IO-CPU
       IF (p_pe == p_io) THEN
          ALLOCATE(gf(dcl%nlon,dcl%nlat))
       ELSE
          NULLIFY(gf)
       END IF
       CALL gather_field(gf, lf)
       !
       ! DISTRIBUTE GLOBAL FIELD TO ALL CPUs
       IF (p_pe /= p_io) THEN
          ALLOCATE(gf(dcl%nlon,dcl%nlat))
       END IF
       CALL p_bcast(gf, p_io)
#else
       ALLOCATE(gf(dcl%nlon,dcl%nlat))
       CALL allgather_field(gf, lf)
#endif
       !
    CASE (-1)
       !
       ! decomposed field <- global field
       !
       IF (.NOT.PRESENT(method)) THEN
          zmethod = M_LOC          ! default
       ELSE
          zmethod = method
       END IF
       !
       IF (zmethod == M_LOC) THEN
          ! PICK OUT SUBSET ON LOCAL CPU
          CALL set_ilon_ilat
          lf(:,:) = 0.0_dp
          DO jrow = 1, dcl%ngpblks
             IF ( jrow == dcl%ngpblks ) THEN
                kproma = dcl%npromz
             ELSE
                kproma = dcl%nproma
             END IF
             DO jp = 1, kproma
                lf(jp,jrow) = gf(ilon(jp, jrow),ilat(jp, jrow))
             END DO
          END DO
          RETURN
       END IF
       !
#ifdef MPIPERF_BAD
       ! I/O-PE
       IF (p_pe == p_io) THEN
          ALLOCATE(zgft(SIZE(gf,1),SIZE(gf,2)))
          ALLOCATE(zgf (SIZE(gf,1),SIZE(gf,2)))
          ALLOCATE(zgfh(SIZE(gf,1),SIZE(gf,2)))
          SELECT CASE (zmethod)
          CASE (M_SUM,M_AVE)
             zgf(:,:)  = gf(:,:)
          CASE (M_STD)
             zgf (:,:) = gf(:,:)
             zgfh(:,:) = gf(:,:)**2
          CASE (M_LOC)
             ! DO NOTHING ; SEE ABOVE
          CASE DEFAULT
             CALL error_bi('invalid METHOD parameter', substr)
          END SELECT
       ELSE
          NULLIFY(zgft)
          NULLIFY(zgf)
          NULLIFY(zgfh)
       END IF

       ! LOOP OVER ALL NON-I/O PEs
       IF (p_pe /= p_io) THEN
          CALL p_send (gf, p_io, tag_tr_gpdc_gpgl)
       ELSE ! p_pe = p_io
          DO i=1, p_nprocs
             pe = dcg(i)%pe
             IF (pe /= p_pe) THEN
                CALL p_recv( zgft, pe, tag_tr_gpdc_gpgl)
                SELECT CASE (zmethod)
                CASE (M_SUM,M_AVE)
                   zgf(:,:)  = zgf(:,:)  + zgft(:,:)
                CASE (M_STD)
                   zgf(:,:)  = zgf(:,:)  + zgft(:,:)
                   zgfh(:,:) = zgfh(:,:) + zgft(:,:)**2
                CASE (M_LOC)
                   ! DO NOTHING ; SEE ABOVE
                CASE DEFAULT
                   CALL error_bi('invalid METHOD parameter', substr)
                END SELECT
             END IF
          END DO
       END IF
       !
       ! RESULT
       IF (p_pe == p_io) THEN
          SELECT CASE (zmethod)
          CASE (M_SUM)
          CASE (M_AVE)
             IF (p_nprocs > 0) THEN
                qave = 1.0_dp/REAL(p_nprocs,dp)
             ELSE
                qave = 1.0_dp
             END IF
             zgf(:,:) = zgf(:,:)*qave
          CASE (M_STD)
             IF (p_nprocs > 0) THEN
                qave = 1.0_dp/REAL(p_nprocs,dp)
             ELSE
                qave = 1.0_dp
             END IF
             IF ((p_nprocs-1) > 0) THEN
                qstd = 1.0_dp/REAL(p_nprocs-1,dp)
             ELSE
                qstd = 0.0_dp
             END IF
             !
             zgf(:,:) = SQRT( ABS( (zgfh(:,:) - &
               REAL(p_nprocs,dp)*(zgf(:,:)*qave)**2)*qstd ) )
          CASE (M_LOC)
             !
          CASE DEFAULT
             CALL error_bi('invalid METHOD parameter', substr)
          END SELECT
       END IF
       !
       ! SCATTER RESULT INTO DECOMPOSED FIELDS
       CALL scatter_gp(zgf, lf, dcg)
       !
       ! CLEAN MEMORY
       IF (ASSOCIATED(zgft)) THEN
          DEALLOCATE(zgft) ; NULLIFY(zgft)
       END IF
       IF (ASSOCIATED(zgf))  THEN
          DEALLOCATE(zgf) ; NULLIFY(zgf)
       END IF
       IF (ASSOCIATED(zgfh)) THEN
          DEALLOCATE(zgfh) ; NULLIFY(zgfh)
       END IF
       !
#else
       ALLOCATE(zgf(SIZE(gf,1),SIZE(gf,2)))
       zgf(:,:) = p_sum(gf)

       ! PICK OUT SUBSET ON LOCAL CPU
       CALL set_ilon_ilat
       ALLOCATE(zlf(SIZE(lf,1),SIZE(lf,2)))
       zlf(:,:) = 0.0_dp
       DO jrow = 1, dcl%ngpblks
          IF ( jrow == dcl%ngpblks ) THEN
             kproma = dcl%npromz
          ELSE
             kproma = dcl%nproma
          END IF
          DO jp = 1, kproma
             zlf(jp,jrow) = zgf(ilon(jp, jrow),ilat(jp, jrow))
          END DO
       END DO

       IF (zmethod == M_STD) THEN
          ALLOCATE(zgfh(SIZE(gf,1),SIZE(gf,2)))
          ALLOCATE(zgft(SIZE(gf,1),SIZE(gf,2)))
          zgft(:,:) = gf(:,:)**2
          zgfh(:,:) = p_sum(zgft)
          ALLOCATE(zlfh(SIZE(lf,1),SIZE(lf,2)))
          zlfh(:,:) = 0.0_dp
          DO jrow = 1, dcl%ngpblks
             IF ( jrow == dcl%ngpblks ) THEN
                kproma = dcl%npromz
             ELSE
                kproma = dcl%nproma
             END IF
             DO jp = 1, kproma
                zlfh(jp,jrow) = zgfh(ilon(jp, jrow),ilat(jp, jrow))
             END DO
          END DO
       END IF

       SELECT CASE (zmethod)
       CASE (M_SUM)
          lf(:,:) = zlf(:,:)
       CASE (M_AVE)
          IF (p_nprocs > 0) THEN
             qave = 1.0_dp/REAL(p_nprocs,dp)
          ELSE
             qave = 1.0_dp
          END IF
          lf(:,:) = zlf(:,:)*qave
       CASE (M_STD)
          IF (p_nprocs > 0) THEN
             qave = 1.0_dp/REAL(p_nprocs,dp)
          ELSE
             qave = 1.0_dp
          END IF
          IF ((p_nprocs-1) > 0) THEN
             qstd = 1.0_dp/REAL(p_nprocs-1,dp)
          ELSE
             qstd = 0.0_dp
          END IF
          !
          zlf(:,:) = SQRT( ABS( (zlfh(:,:) - &
               REAL(p_nprocs,dp)*(zlf(:,:)*qave)**2)*qstd ) )
       CASE (M_LOC)
          !
       CASE DEFAULT
          CALL error_bi('invalid METHOD parameter', substr)
       END SELECT

       ! CLEAN
       IF (ASSOCIATED(zgf))  THEN
          DEALLOCATE(zgf) ; NULLIFY(zgf)
       END IF
       IF (ASSOCIATED(zgfh)) THEN
          DEALLOCATE(zgfh) ; NULLIFY(zgfh)
       END IF
       IF (ASSOCIATED(zgft)) THEN
          DEALLOCATE(zgft) ; NULLIFY(zgft)
       END IF
       IF (ASSOCIATED(zlf))  THEN
          DEALLOCATE(zlf) ; NULLIFY(zlf)
       END IF
       IF (ASSOCIATED(zlfh)) THEN
          DEALLOCATE(zlfh) ; NULLIFY(zlfh)
       END IF
#endif

    CASE default
       CALL error_bi ('invalid SIGN parameter (not 1,-1)', substr)
    END SELECT

  END SUBROUTINE trp_gpdc_gpgl_2d
  ! =========================================================================

! ##########################################################################
! END #ifndef (MBM_CLAMS)
! ##########################################################################
#endif

  ! =========================================================================
  SUBROUTINE scatter_glix_1d(gl, lc, xpe, XNG)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    REAL(dp), DIMENSION(:), POINTER  :: gl
    REAL(dp), DIMENSION(:), POINTER  :: lc
    INTEGER, INTENT(IN), OPTIONAL    :: xpe
    INTEGER, INTENT(IN), OPTIONAL    :: XNG

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr='scatter_glix_1d'
    INTEGER, DIMENSION(:,:), POINTER :: idx=>NULL()
    LOGICAL, DIMENSION(:),   POINTER :: ldo=>NULL()
    INTEGER                          :: NG
    INTEGER                          :: NL
    INTEGER                          :: pe, zxpe
    INTEGER                          :: i

    IF (PRESENT(xpe)) THEN
       zxpe = xpe
    ELSE
       zxpe = p_io
    END IF

    IF (PRESENT(XNG)) THEN
       NG = XNG
    ELSE
       IF (p_pe == zxpe) THEN
          NG = SIZE(gl)
       END IF
       CALL p_bcast(NG, zxpe)
    ENDIF

    CALL get_dc_index(NG, idx, ldo)
    NL = idx(p_pe,2)-idx(p_pe,1) + 1

    IF (ldo(p_pe)) THEN
       IF (.NOT.ASSOCIATED(lc)) THEN
          ALLOCATE(lc(NL))
       ELSE
          IF (SIZE(lc) /= NL) THEN
             CALL error_bi('wrong size of local array', substr)
          END IF
       END IF
    END IF

    !
    ! send if p_pe = zxpe
    !
    IF (p_pe == zxpe) THEN
       DO i = 1, p_nprocs
#if defined(MBM_CLAMS)
          pe = i-1
#else
          pe = dcg(i)% pe
#endif
          IF (pe /= p_pe) THEN
             IF (ldo(pe)) &
                  CALL p_send(gl(idx(pe,1):idx(pe,2)), pe, tag_scatter_ix)
          ELSE
             IF (ldo(pe)) &
                  lc(1:NL) = gl(idx(pe,1):idx(pe,2))
          END IF
       END DO
    ELSE
       ! receive if p_pe /= zxpe
       IF (ldo(p_pe)) CALL p_recv(lc(1:NL), zxpe, tag_scatter_ix)
    END IF

    ! CLEAN UP
    DEALLOCATE(idx)
    DEALLOCATE(ldo)

  END SUBROUTINE scatter_glix_1d
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE scatter_glix_4d(gl, lc, index, xpe, xishpg)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRODUCT, SIZE, UBOUND

    ! I/O
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: gl    ! global field
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: lc    ! local field
    INTEGER, INTENT(IN)                   :: index ! decomposition index
    INTEGER, INTENT(IN), OPTIONAL         :: xpe
    INTEGER, DIMENSION(4), INTENT(IN), OPTIONAL :: XISHPG

    ! LOCAL
    INTEGER, PARAMETER               :: nrank = 4
    CHARACTER(LEN=*), PARAMETER      :: substr='scatter_glix_4d'
    INTEGER, DIMENSION(:,:), POINTER :: idx=>NULL()
    LOGICAL, DIMENSION(:),   POINTER :: ldo=>NULL()
    INTEGER                          :: NG
    INTEGER                          :: pe, zxpe
    INTEGER                          :: i
    INTEGER                          :: ishpg(4) ! global shape
    INTEGER                          :: ishpl(4) ! local shape
    INTEGER                          :: ixg(4,2) ! global start/stop index
    INTEGER                          :: ixl(4,2) ! local start/stop index

    IF ((index < 1).OR.(index > nrank)) THEN
       CALL error_bi('index out of range', substr)
    END IF

    IF (PRESENT(xpe)) THEN
       zxpe = xpe
    ELSE
       zxpe = p_io
    END IF

    ! SHAPE OF 4-DARRAY
    IF (PRESENT(XISHPG)) THEN
       ishpg(:) = XISHPG(:)
    ELSE
       IF (p_pe == zxpe) THEN
          DO i=1, nrank
             ishpg(i) = SIZE(gl,i)
          END DO
       END IF
       CALL p_bcast(ishpg, zxpe)
    ENDIF
    NG = ishpg(index)                ! size at decomp. index of global field

    CALL get_dc_index(NG, idx, ldo)
    ishpl(:) = ishpg(:)              ! size of local field is size of global
                                     ! field, except at decomp. index
    ishpl(index) = idx(p_pe,2)-idx(p_pe,1) + 1

    DO i=1, nrank
       ixl(i,1) = 1                  ! start/stop index in local field ...
       ixl(i,2) = ishpl(i)           ! ... is given by shape
       ixg(i,1) = 1                  ! start/stop index in global field
       ixg(i,2) = ishpg(i)           ! is size of global field ...
                                     ! ... except at decomp. index
                                     ! -> set below
    END DO

    IF (ldo(p_pe)) THEN
       IF (.NOT.ASSOCIATED(lc)) THEN
          ALLOCATE(lc(ishpl(1),ishpl(2),ishpl(3),ishpl(4)))
       ELSE
          IF (SIZE(lc) /= PRODUCT(ishpl)) THEN
             WRITE(*,*) 'p_pe=',p_pe,': ',ishpl,' =/= ',UBOUND(lc)
             CALL error_bi('wrong size of local array', substr)
          END IF
       END IF
    END IF

    !
    ! send if p_pe = zxpe
    !
    IF (p_pe == zxpe) THEN
       DO i = 1, p_nprocs
#if defined(MBM_CLAMS)
          pe = i-1
#else
          pe = dcg(i)% pe
#endif
          ! set global start/stop indices at decomp. index
          ixg(index, 1) = idx(pe,1)
          ixg(index, 2) = idx(pe,2)
          IF (pe /= p_pe) THEN
             IF (ldo(pe)) &
                  CALL p_send(gl(ixg(1,1):ixg(1,2)  &
                  ,   ixg(2,1):ixg(2,2)  &
                  ,   ixg(3,1):ixg(3,2)  &
                  ,   ixg(4,1):ixg(4,2)) &
                  , pe, tag_scatter_ix)
          ELSE
             IF (ldo(pe)) &
                  lc(ixl(1,1):ixl(1,2)  &
                  ,   ixl(2,1):ixl(2,2)  &
                  ,   ixl(3,1):ixl(3,2)  &
                  ,   ixl(4,1):ixl(4,2)) = &
                  gl(ixg(1,1):ixg(1,2)   &
                  ,   ixg(2,1):ixg(2,2)  &
                  ,   ixg(3,1):ixg(3,2)  &
                  ,   ixg(4,1):ixg(4,2))
          END IF
       END DO
    ELSE
       ! reveive if p_pe /= zxpe
       IF (ldo(p_pe)) CALL p_recv(lc(ixl(1,1):ixl(1,2)  &
            ,   ixl(2,1):ixl(2,2)  &
            ,   ixl(3,1):ixl(3,2)  &
            ,   ixl(4,1):ixl(4,2)) &
            , zxpe, tag_scatter_ix)
    END IF

    ! CLEAN UP
    DEALLOCATE(idx)
    DEALLOCATE(ldo)

  END SUBROUTINE scatter_glix_4d
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE gather_glix_1d(gl, lc, xpe, XNG)

#ifdef _MPI_PSUM
    USE messy_main_mpi_bi, ONLY: p_sum
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    REAL(dp), DIMENSION(:), POINTER  :: gl
    REAL(dp), DIMENSION(:), POINTER  :: lc
    INTEGER, INTENT(IN), OPTIONAL    :: xpe
    INTEGER, INTENT(IN), OPTIONAL    :: XNG

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr='gather_glix_1d'
    INTEGER, DIMENSION(:,:), POINTER :: idx=>NULL()
    LOGICAL, DIMENSION(:),   POINTER :: ldo=>NULL()
    INTEGER                          :: NG
    INTEGER                          :: NL
    INTEGER                          :: pe, zxpe
    INTEGER                          :: i
#ifndef _MPI_PSUM
    INTEGER                          :: NLs
#endif

    IF (PRESENT(xpe)) THEN
       zxpe = xpe
    ELSE
       zxpe = p_io
    END IF


    IF (PRESENT(XNG)) THEN
       NG = XNG
    ELSE
       NL = SIZE(lc)      ! size of local field

#ifndef _MPI_PSUM
       IF (p_pe == zxpe) THEN
          NG = NL         ! size of global field
       END IF

       ! SUM OVER ALL NON-zxpe p_pe
       IF (p_pe /= zxpe) THEN
          CALL p_send (NL, zxpe, tag_gather_ix)
       ELSE ! p_pe = zxpe
          DO i=1, p_nprocs
#if defined(MBM_CLAMS)
             pe = i-1
#else
             pe = dcg(i)%pe
#endif
             IF (pe /= p_pe) THEN
                CALL p_recv( NLs, pe, tag_gather_ix)
                NG = NG + NLs
             END IF
          END DO
       END IF

       CALL p_bcast(NG, zxpe)
#else
       NG = p_sum(NL)
#endif
    ENDIF

    CALL get_dc_index(NG, idx, ldo)

    IF (p_pe == zxpe) THEN
       IF (ASSOCIATED(gl)) THEN
          IF (SIZE(gl) /= NG) THEN
             CALL error_bi('wrong size of global array', substr)
          END IF
       ELSE
          ALLOCATE(gl(NG))
       END IF
    ELSE
       IF (ASSOCIATED(gl)) DEALLOCATE(gl)
       NULLIFY(gl)
    ENDIF

    ! send if pe /= zxpe
    !
    IF (p_pe /= zxpe) THEN
       IF (ldo(p_pe)) CALL p_send (lc(:), zxpe, tag_gather_ix)
    ELSE ! p_pe = zxpe
       DO i=1, p_nprocs
#if defined(MBM_CLAMS)
          pe = i-1
#else
          pe = dcg(i)%pe
#endif
          IF (pe /= p_pe) THEN
             IF (ldo(pe)) &
                  CALL p_recv( gl(idx(pe,1):idx(pe,2)), pe, tag_gather_ix)
          ELSE
             IF (ldo(pe)) &
                  gl(idx(pe,1):idx(pe,2)) = lc(:)
          END IF
       END DO
    END IF

    ! CLEAN UP
    DEALLOCATE(idx)
    DEALLOCATE(ldo)

  END SUBROUTINE gather_glix_1d
  ! =========================================================================
  ! =========================================================================
  SUBROUTINE gather_glix_4d(gl, lc, index, xpe, XNG)

#ifdef _MPI_PSUM
    USE messy_main_mpi_bi, ONLY: p_sum
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRODUCT, SIZE, UBOUND

    ! I/O
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: gl     ! global field
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: lc     ! local field
    INTEGER, INTENT(IN)                   :: index  ! decomposition index
    INTEGER, INTENT(IN), OPTIONAL         :: xpe
    INTEGER, INTENT(IN), OPTIONAL         :: XNG

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr='gather_glix_4d'
    INTEGER, PARAMETER               :: nrank = 4
    INTEGER, DIMENSION(:,:), POINTER :: idx=>NULL()
    LOGICAL, DIMENSION(:),   POINTER :: ldo=>NULL()
    INTEGER                          :: NG, NL
    INTEGER                          :: pe, zxpe
    INTEGER                          :: i
    INTEGER                          :: ishpg(4) ! global shape
    INTEGER                          :: ishpl(4) ! local shape
    INTEGER                          :: ixg(4,2) ! global start/stop index
    INTEGER                          :: ixl(4,2) ! local start/stop index
#ifndef _MPI_PSUM
    INTEGER                          :: NLs
#endif

    IF ((index < 1).OR.(index > nrank)) THEN
       CALL error_bi('index out of range', substr)
    END IF

    IF (.NOT. ASSOCIATED(lc)) THEN
       CALL error_bi('pointer to local field is not associated', substr)
    END IF

    IF (PRESENT(xpe)) THEN
       zxpe = xpe
    ELSE
       zxpe = p_io
    END IF

    ! SHAPE OF 4-DARRAY
    DO i=1, nrank
       ishpl(i) = SIZE(lc,i)
    END DO

    IF (PRESENT(XNG)) THEN
       NG = XNG
    ELSE
       NL = ishpl(index)         ! size at decomp. index of local field

#ifndef _MPI_PSUM
       IF (p_pe == zxpe) THEN
          NG = NL        ! size at decomp. index of global field (summation)
       END IF
       !
       ! SUM OVER ALL NON-zxpe p_pe
       IF (p_pe /= zxpe) THEN
          CALL p_send (NL, zxpe, tag_gather_ix)
       ELSE ! p_pe = zxpe
          DO i=1, p_nprocs
#if defined(MBM_CLAMS)
             pe = i-1
#else
             pe = dcg(i)%pe
#endif
             IF (pe /= p_pe) THEN
                CALL p_recv( NLs, pe, tag_gather_ix)
                NG = NG + NLs
             END IF
          END DO
       END IF
#else
       NG = p_sum(NL)
#endif
    ENDIF

    ! lcoal field start/stop indices given by size
    DO i=1, nrank
       ixl(i,1) = 1
       ixl(i,2) = SIZE(lc, i)
    END DO

#ifndef _MPI_PSUM
    CALL p_bcast(NG, zxpe)     ! global field dim. at decomp. index
#endif
    ishpg(:) = ishpl(:)        ! global shape is local shape ...
    ishpg(index) = NG          ! ... except at decomposition index

    ! global field start/stop indices
    ixg(:,:) = ixl(:,:)        ! global start/stop index is local start/stop
    ! index, except at decomposition index
    CALL get_dc_index(NG, idx, ldo)
    ! -> set below

    IF (p_pe == zxpe) THEN
       IF (ASSOCIATED(gl)) THEN
          IF (SIZE(gl) /= PRODUCT(ishpg)) THEN
             WRITE(*,*) 'SIZE(gl) = ',SIZE(gl),' =/= ',PRODUCT(ishpg)
             WRITE(*,*) 'SHAPE OF GLOBAL ARRAY: ',UBOUND(gl)
             WRITE(*,*) 'EXPECTED SHAPE       : ',ishpg(:)
             CALL error_bi('wrong size of global array', substr)
          END IF
       ELSE
          ALLOCATE(gl(ishpg(1),ishpg(2),ishpg(3),ishpg(4)))
       END IF
    ELSE
       IF (ASSOCIATED(gl)) DEALLOCATE(gl)
       NULLIFY(gl)
    ENDIF

    ! send if pe /= zxpe
    !
    IF (p_pe /= zxpe) THEN
       IF (ldo(p_pe)) CALL p_send (lc(ixl(1,1):ixl(1,2)  &
            ,  ixl(2,1):ixl(2,2)  &
            ,  ixl(3,1):ixl(3,2)  &
            ,  ixl(4,1):ixl(4,2)) &
            , zxpe, tag_gather_ix)
    ELSE ! p_pe = zxpe
       DO i=1, p_nprocs
#if defined(MBM_CLAMS)
          pe = i-1
#else
          pe = dcg(i)%pe
#endif
          ! set global start/stop indices at decomp. index
          ixg(index,1) = idx(pe,1)
          ixg(index,2) = idx(pe,2)
          IF (pe /= p_pe) THEN
             IF (ldo(pe)) &
                  CALL p_recv( gl(ixg(1,1):ixg(1,2)  &
                  ,   ixg(2,1):ixg(2,2)  &
                  ,   ixg(3,1):ixg(3,2)  &
                  ,   ixg(4,1):ixg(4,2)) &
                  , pe, tag_gather_ix)
          ELSE
             IF (ldo(pe)) &
                  gl(ixg(1,1):ixg(1,2)     &
                  ,   ixg(2,1):ixg(2,2)    &
                  ,   ixg(3,1):ixg(3,2)    &
                  ,   ixg(4,1):ixg(4,2)) = &
                  lc(ixl(1,1):ixl(1,2)  &
                  ,  ixl(2,1):ixl(2,2)  &
                  ,  ixl(3,1):ixl(3,2)  &
                  ,  ixl(4,1):ixl(4,2))
          END IF
       END DO
    END IF

    ! CLEAN UP
    DEALLOCATE(idx)
    DEALLOCATE(ldo)

  END SUBROUTINE gather_glix_4d
  ! =========================================================================

#ifndef MBM_CLAMS
    ! ========================================================================
  SUBROUTINE trf_ltd_a_0(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (0 PARAMETER)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !
    !
    !
    !
    !                    ---
    !                    \
    !      LS ([m, n]) =  >    FA(m, mu) * Q([m,n], mu)   ; n: even
    !                    /
    !                    ---
    !                    mu
    !
    !
    !                    ---
    !                    \
    !      LS ([m, n]) =  >    FS(m, mu) * Q([m,n], mu)   ; n: odd
    !                    /
    !                    ---
    !                    mu
    !

    !  same as trf_ltd_0, but with the parameter anmt insted of pnmt

    ! ECHAM5
    USE mo_legendre,      ONLY: legmod

    IMPLICIT NONE

    ! I/O
    REAL(dp), DIMENSION(:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:), TARGET, INTENT(IN)  :: fs  ! symmetric


    ! LOCAL
    REAL(dp), DIMENSION(:,:), POINTER  :: fd=>NULL()
    ! LOOP VARAIABLES
    INTEGER                        :: nlnsp  ! number of (m,n) pairs
    INTEGER                        :: nhgl   ! 1/2 number of latitudes
    INTEGER                        :: nlm    ! number of wave numbers m
    ! ECHAM5 SPECIFIC DECOMPOSITION INFORMATION
    INTEGER ,POINTER :: nlmp(:)=>NULL() ! offset to local m columns
    INTEGER ,POINTER :: nlnp(:)=>NULL() ! column length (No. of n's)
    ! LOCAL SCALARS
    LOGICAL, SAVE                  :: ini_flag = .FALSE. ! P(m,n,mu) calculated

    INTEGER :: jh           ! LOOP OVER HEMISPHERES
    INTEGER :: iu
    INTEGER :: jm           ! LOOP OVER WAVE NUMBERS m
    INTEGER :: knm          ! INDEX FOR (m,n=m) OR (m,n=m+1)
    INTEGER :: knT          ! INDEX FOR (m,n=T) OR (m,n=T-1)

    !  LOCAL ARRAYS
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE :: anmt

    !  INTRINSIC FUNCTIONS
    INTRINSIC :: MATMUL, MOD, SIZE

    ! INIT
    nlnsp = SIZE(ls)    ! dcl% lnsp
    nlm   = SIZE(fa,1)  ! dcl% nlm
    nhgl  = SIZE(fa,2)  ! dcl% nlat/2
    !

    !-- ECHAM5 specific decomposition information
    nlmp    => dcl% nlmp     ! displacement of the first point of columns
    nlnp    => dcl% nlnp     ! number of points on each column

    IF ( .NOT. ini_flag ) THEN

       ALLOCATE(anmt(nhgl, nlnsp))

       CALL legmod(anmt=anmt)

       ini_flag = .TRUE.

    ENDIF

    DO jh = 1, 2
       iu = 2 - jh

       IF (jh==1) THEN
          fd  => fs (:,:)
       ELSE
          fd  => fa (:,:)
       END IF

       DO jm = 1, nlm

          knm = nlmp(jm) - iu + 2
          knT = nlmp(jm) + nlnp(jm)
          knT = knT - MOD(knT-knm,2)

          IF (knm > nlnsp) CYCLE

          ls(knm:knT:2) = MATMUL(fd(jm,:), anmt(:,knm:knT:2))

          !WRITE(*,*) 'TRF_LTD (pe =',p_pe,'): ',knm, knT, jm, ' --> ', jh

       END DO

    END DO

  END SUBROUTINE trf_ltd_a_0
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ltd_a_1(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (1 PARAMETER SUCH AS Re,Im OR level)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(IN)  :: fs  ! symmetric

    ! LOCAL
    INTEGER :: n1, i1

    n1 = SIZE(ls,1)

    DO i1 = 1, n1
       CALL trf_ltd_a_0(ls(i1,:),fa(i1,:,:),fs(i1,:,:))
    END DO

  END SUBROUTINE trf_ltd_a_1
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ltd_a_2(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (2 PARAMETER SUCH AS Re,Im AND level)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(IN)  :: fs  ! symmetric

    ! LOCAL
    INTEGER :: n1, i1, n2, i2

    n1 = SIZE(ls,1)
    n2 = SIZE(ls,2)

    DO i2 = 1, n2
       DO i1 = 1, n1
          CALL trf_ltd_a_0(ls(i1,i2,:),fa(i1,i2,:,:),fs(i1,i2,:,:))
       END DO
    END DO

  END SUBROUTINE trf_ltd_a_2
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ltd_r_0(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (0 PARAMETER)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !
    !
    !
    !
    !                    ---
    !                    \
    !      LS ([m, n]) =  >    FA(m, mu) * Q([m,n], mu)   ; n: even
    !                    /
    !                    ---
    !                    mu
    !
    !
    !                    ---
    !                    \
    !      LS ([m, n]) =  >    FS(m, mu) * Q([m,n], mu)   ; n: odd
    !                    /
    !                    ---
    !                    mu
    !

    !  same as trf_ltd_0, but with the parameter rnmt insted of pnmt

    ! ECHAM5
    USE mo_legendre,      ONLY: legmod

    IMPLICIT NONE

    ! I/O
    REAL(dp), DIMENSION(:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:), TARGET, INTENT(IN)  :: fs  ! symmetric


    ! LOCAL
    REAL(dp), DIMENSION(:,:), POINTER  :: fd=>NULL()
    ! LOOP VARAIABLES
    INTEGER                        :: nlnsp  ! number of (m,n) pairs
    INTEGER                        :: nhgl   ! 1/2 number of latitudes
    INTEGER                        :: nlm    ! number of wave numbers m
    ! ECHAM5 SPECIFIC DECOMPOSITION INFORMATION
    INTEGER ,POINTER :: nlmp(:)=>NULL() ! offset to local m columns
    INTEGER ,POINTER :: nlnp(:)=>NULL() ! column length (No. of n's)
    ! LOCAL SCALARS
    LOGICAL, SAVE                  :: ini_flag = .FALSE. ! P(m,n,mu) calculated

    INTEGER :: jh           ! LOOP OVER HEMISPHERES
    INTEGER :: iu
    INTEGER :: jm           ! LOOP OVER WAVE NUMBERS m
    INTEGER :: knm          ! INDEX FOR (m,n=m) OR (m,n=m+1)
    INTEGER :: knT          ! INDEX FOR (m,n=T) OR (m,n=T-1)

    !  LOCAL ARRAYS
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE :: rnmt
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE, SAVE :: pnmt

    !  INTRINSIC FUNCTIONS
    INTRINSIC :: MATMUL, MOD, SIZE

    ! INIT
    nlnsp = SIZE(ls)    ! dcl% lnsp
    nlm   = SIZE(fa,1)  ! dcl% nlm
    nhgl  = SIZE(fa,2)  ! dcl% nlat/2
    !

    !-- ECHAM5 specific decomposition information
    nlmp    => dcl% nlmp     ! displacement of the first point of columns
    nlnp    => dcl% nlnp     ! number of points on each column

    IF ( .NOT. ini_flag ) THEN

       ALLOCATE(rnmt(nhgl, nlnsp))
       ALLOCATE(pnmt(nhgl, nlnsp))

       CALL legmod(pnmt=pnmt, rnmt=rnmt)

       ini_flag = .TRUE.

    ENDIF

    DO jh = 1, 2
       iu = 2 - jh

       IF (jh==1) THEN
          fd  => fs (:,:)
       ELSE
          fd  => fa (:,:)
       END IF

       DO jm = 1, nlm

          knm = nlmp(jm) - iu + 2
          knT = nlmp(jm) + nlnp(jm)
          knT = knT - MOD(knT-knm,2)

          IF (knm > nlnsp) CYCLE

          ls(knm:knT:2) = MATMUL(fd(jm,:), rnmt(:,knm:knT:2))

          !WRITE(*,*) 'TRF_LTD (pe =',p_pe,'): ',knm, knT, jm, ' --> ', jh

       END DO

    END DO

  END SUBROUTINE trf_ltd_r_0
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ltd_r_1(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (1 PARAMETER SUCH AS Re,Im OR level)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(IN)  :: fs  ! symmetric

    ! LOCAL
    INTEGER :: n1, i1

    n1 = SIZE(ls,1)

    DO i1 = 1, n1
       CALL trf_ltd_r_0(ls(i1,:),fa(i1,:,:),fs(i1,:,:))
    END DO

  END SUBROUTINE trf_ltd_r_1
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE trf_ltd_r_2(ls, fa, fs)

    ! DIRECT LEGENDRE TRANSFORMATION (2 PARAMETER SUCH AS Re,Im AND level)
    !
    ! Authors: P. Joeckel, MPICH, Jan 2004
    !

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:),           INTENT(OUT) :: ls
    REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(IN)  :: fa  ! anti-symmetric
    REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(IN)  :: fs  ! symmetric

    ! LOCAL
    INTEGER :: n1, i1, n2, i2

    n1 = SIZE(ls,1)
    n2 = SIZE(ls,2)

    DO i2 = 1, n2
       DO i1 = 1, n1
          CALL trf_ltd_r_0(ls(i1,i2,:),fa(i1,i2,:,:),fs(i1,i2,:,:))
       END DO
    END DO

  END SUBROUTINE trf_ltd_r_2
  ! ========================================================================

! --------------------------------------------------------------------
  SUBROUTINE zonal_average_3d(f, fza)

    ! CALCULATION OF ZONAL AVERAGE OF A 3D FIELD IN
    ! GRIDPOINT DECOMPOSITION
    !
    ! THE RESULT IS A FIELD OF EQUAL RANK/SIZE WITH THE ZONAL
    ! AVERAGE VALUE AT EACH POINT
    !
    ! Author: Patrick Joeckel, MPICH, December 2003

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, ngpblks, nproma, npromz
    USE messy_main_grid_def_bi,      ONLY: ilat
    USE messy_main_blather_bi,       ONLY: error_bi

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, REAL, SIZE, SUM

    ! I/O
    REAL(dp), DIMENSION(:,:,:), POINTER :: f   ! 3D-field (decomp.) (IN)
    REAL(dp), DIMENSION(:,:,:), POINTER :: fza ! zonal average (str. dec.) (OUT)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER        :: substr = 'zonal_average_3d'
    ! 'global' fields on all CPUs
    REAL(dp), DIMENSION(:,:,:), POINTER    :: fg   => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER    :: fgza => NULL()
    INTEGER :: jk, jp, zjrow, zkproma

    IF (ASSOCIATED(fza)) THEN
       IF (SIZE(fza,1) /= SIZE(f,1) ) &
            CALL error_bi('WRONG DIMENSION LENGTH (1)', substr)
       IF (SIZE(fza,2) /= SIZE(f,2) ) &
            CALL error_bi('WRONG DIMENSION LENGTH (2)', substr)
       IF (SIZE(fza,3) /= SIZE(f,3) ) &
            CALL error_bi('WRONG DIMENSION LENGTH (3)', substr)
    ELSE
       ALLOCATE(fza(SIZE(f,1),SIZE(f,2),SIZE(f,3)))
    END IF
    fza(:,:,:) = 0.0

    CALL trp_gpdc_gpgl(1, f, fg)           ! fg   -> nlon x nlev x nlat
    ALLOCATE(fgza(SIZE(fg,2),SIZE(fg,3)))  ! fgza ->        nlev x nlat

    ! UNWEIGHTED ZONAL AVERAGE
    fgza(:,:) = SUM(fg(:,:,:),1)/REAL(SIZE(fg,1),dp)

    ! RE-DISTRIBUTE INTO DECOMPOSITION
    DO zjrow=1, ngpblks

       IF ( zjrow == ngpblks ) THEN
          zkproma = npromz
       ELSE
          zkproma = nproma
       END IF

       DO jk=1, nlev
          DO jp=1, zkproma

             fza(jp, jk, zjrow) = fgza(jk, ilat(jp, zjrow))

          END DO
       END DO

    END DO

    IF (ASSOCIATED(fgza)) THEN
       DEALLOCATE(fgza) ; NULLIFY(fgza)
    END IF
    IF (ASSOCIATED(fg)) THEN
       DEALLOCATE(fg) ; NULLIFY(fg)
    END IF

  END SUBROUTINE zonal_average_3d
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE zonal_average_2d(f, fza)

    ! CALCULATION OF ZONAL AVERAGE OF A 2D FIELD IN
    ! GRIDPOINT DECOMPOSITION
    !
    ! THE RESULT IS A FIELD OF EQUAL RANK/SIZE WITH THE ZONAL
    ! AVERAGE VALUE AT EACH POINT
    !
    ! Author: Patrick Joeckel, MPICH, December 2003

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi,  ONLY: ngpblks, nproma, npromz
    USE messy_main_grid_def_bi,      ONLY: ilat
    USE messy_main_blather_bi,       ONLY: error_bi

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, REAL, SIZE, SUM

    ! I/O
    REAL(dp), DIMENSION(:,:), POINTER :: f    ! 3D-field (decomp.) (IN)
    REAL(dp), DIMENSION(:,:), POINTER :: fza  ! zonal average (str. dec.) (OUT)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER        :: substr = 'zonal_average_2d'
    ! 'global' fields on all CPUs
    REAL(dp), DIMENSION(:,:), POINTER    :: fg   => NULL()
    REAL(dp), DIMENSION(:),   POINTER    :: fgza => NULL()
    INTEGER :: jp, zjrow, zkproma

    IF (ASSOCIATED(fza)) THEN
       IF (SIZE(fza,1) /= SIZE(f,1) ) &
            CALL error_bi('WRONG DIMENSION LENGTH (1)', substr)
       IF (SIZE(fza,2) /= SIZE(f,2) ) &
            CALL error_bi('WRONG DIMENSION LENGTH (2)',substr)
    ELSE
       ALLOCATE(fza(SIZE(f,1),SIZE(f,2)))
    END IF
    fza(:,:) = 0.0

    CALL trp_gpdc_gpgl(1, f, fg)            ! fg   -> nlon x nlat
    ALLOCATE(fgza(SIZE(fg,2)))              ! fgza ->        nlat

    ! UNWEIGHTED ZONAL AVERAGE
    fgza(:) = SUM(fg(:,:),1)/REAL(SIZE(fg,1),dp)

    ! RE-DISTRIBUTE INTO DECOMPOSITION
    DO zjrow=1, ngpblks

       IF ( zjrow == ngpblks ) THEN
          zkproma = npromz
       ELSE
          zkproma = nproma
       END IF

       DO jp=1, zkproma

          fza(jp, zjrow) = fgza(ilat(jp, zjrow))

       END DO

    END DO

    IF (ASSOCIATED(fgza)) THEN
       DEALLOCATE(fgza) ; NULLIFY(fgza)
    END IF
    IF (ASSOCIATED(fg)) THEN
       DEALLOCATE(fg) ; NULLIFY(fg)
    END IF

  END SUBROUTINE zonal_average_2d
  ! --------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE locate_in_decomp_nrst(status, lon, lat, pe, jp, jrow, jgx, jgy)

    ! LOCATES GEOGRAPHICAL POSITION IN VECTOR-PARALLEL-DECOMPOSITION
    !
    ! Author: Patrick Joeckel, MPICH, December 2004

    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_grid_def_mem_bi,  ONLY: ngpblks, nproma
    USE messy_main_grid_def_bi,      ONLY: philat, philon
    USE messy_main_tools,            ONLY: nn_index

    IMPLICIT NONE
    INTRINSIC :: PRESENT, REAL

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    REAL(DP), INTENT(IN)  :: lon     ! LONGITUDE [-180 ... 180]
    REAL(DP), INTENT(IN)  :: lat     ! LATITUDE  [-90 ... 90]
    INTEGER,  INTENT(OUT) :: pe      ! PROCESS ID
    INTEGER,  INTENT(OUT) :: jp      ! COLUMN
    INTEGER,  INTENT(OUT) :: jrow    ! ROW
    INTEGER,  INTENT(OUT), OPTIONAL :: jgx  ! lon index in global field
    INTEGER,  INTENT(OUT), OPTIONAL :: jgy  ! lat index in global field

    ! LOCAL
    LOGICAL,  SAVE                          :: linit = .FALSE.
    REAL(DP), DIMENSION(:,:), POINTER, SAVE :: gpe   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER, SAVE :: gjp   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER, SAVE :: gjrow => NULL()
    REAL(DP), DIMENSION(:,:), POINTER       :: lpe   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER       :: ljp   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER       :: ljrow => NULL()
    INTEGER                                 :: zjrow, zkproma
    REAL(DP)                                :: zlon
    INTEGER                                 :: jjgx, jjgy

    INTRINSIC :: NINT

    ! INITIALISE
    status = 0
    IF (Present(jgx))  jgx = -1
    IF (Present(jgy))  jgy = -1
    PE   = -1
    jp   = -1
    jrow = -1

    ! CHECKS
    IF ((lon < -180.0_DP) .OR. (lon > 360.0_DP)) THEN
       status = 1  ! LONGITUDE OUT OF RANGE
       RETURN
    ELSE
       ! CHANGE INTERVAL [-180 ... 180] -> [0 ... 360]
       IF (lon < 0.0_DP) THEN
          zlon = lon + 360.0_DP
       ELSE
          zlon = lon
       END IF
    END IF
    !
    IF ((lat < -90.0_DP) .OR. (lat > 90.0_DP)) THEN
       status = 2  ! LATITUDE OUT OF RANGE
       RETURN
    END IF

    ! INITIALIZATION
    IF (.NOT. linit) THEN
       ! FIELDS NOT YET INITIALIZED
       ALLOCATE(lpe  (nproma, ngpblks))
       ALLOCATE(ljp  (nproma, ngpblks))
       ALLOCATE(ljrow(nproma, ngpblks))
       lpe(:,:) = REAL(p_pe, DP)
       DO zjrow=1, ngpblks
          DO zkproma=1, nproma
             ljp  (zkproma,zjrow) = REAL(zkproma, DP)
             ljrow(zkproma,zjrow) = REAL(zjrow,   DP)
          END DO
       END DO
       !
       CALL trp_gpdc_gpgl(1, lpe,   gpe)
       CALL trp_gpdc_gpgl(1, ljp,   gjp)
       CALL trp_gpdc_gpgl(1, ljrow, gjrow)
       !
       DEALLOCATE(lpe) ; NULLIFY(lpe)
       DEALLOCATE(ljp) ; NULLIFY(ljp)
       DEALLOCATE(ljrow) ; NULLIFY(ljrow)
       !
       linit = .TRUE.
       !
    END IF

    ! SET INDICES
    CALL nn_index(philon, zlon, jjgx)
    CALL nn_index(philat,  lat, jjgy)
    pe     = NINT(gpe  (jjgx,jjgy))
    jp     = NINT(gjp  (jjgx,jjgy))
    jrow   = NINT(gjrow(jjgx,jjgy))

    IF (PRESENT(jgx)) jgx = jjgx
    IF (PRESENT(jgy)) jgy = jjgy

    status = 0

  END SUBROUTINE locate_in_decomp_nrst
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE locate_in_decomp_4(status, lon, lat, pe, jp, jrow, w, gjx, gjy)

    ! LOCATES GEOGRAPHICAL POSITION IN VECTOR-PARALLEL-DECOMPOSITION
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006

    USE messy_main_mpi_bi,          ONLY: p_pe
    USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma
    USE messy_main_grid_def_bi,     ONLY: philat, philon
    USE messy_main_tools,           ONLY: nn_index, bilin_weight

    IMPLICIT NONE
    INTRINSIC :: PRESENT, REAL, MAX, MIN

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    REAL(DP),               INTENT(IN)  :: lon     ! LONGITUDE [-180 ... 180]
    REAL(DP),               INTENT(IN)  :: lat     ! LATITUDE  [-90 ... 90]
    INTEGER,  DIMENSION(4), INTENT(OUT) :: pe      ! PROCESS ID
    INTEGER,  DIMENSION(4), INTENT(OUT) :: jp      ! COLUMN
    INTEGER,  DIMENSION(4), INTENT(OUT) :: jrow    ! ROW
    REAL(DP), DIMENSION(4), INTENT(OUT) :: w       ! WEIGHT
    ! lon index in global field
    INTEGER,  DIMENSION(4), INTENT(OUT), OPTIONAL :: gjx
    ! lat index in global field
    INTEGER,  DIMENSION(4), INTENT(OUT), OPTIONAL :: gjy

    ! LOCAL
    LOGICAL,  SAVE                          :: linit = .FALSE.
    REAL(DP), DIMENSION(:,:), POINTER, SAVE :: gpe   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER, SAVE :: gjp   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER, SAVE :: gjrow => NULL()
    REAL(DP), DIMENSION(:,:), POINTER       :: lpe   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER       :: ljp   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER       :: ljrow => NULL()
    INTEGER                                 :: zjrow, zkproma
    REAL(DP)                                :: zlon
    INTEGER                                 :: jgx1, jgy1, jgx2, jgy2
    INTEGER                                 :: jgx, jgy
    REAL(DP), DIMENSION(4,2)                :: vn
    REAL(DP), DIMENSION(2)                  :: v

    INTRINSIC :: NINT

    ! INITIALISE
    status = 0
    IF (PRESENT(gjx))  gjx = -1
    IF (PRESENT(gjy))  gjy = -1
    PE   = -1
    jp   = -1
    jrow = -1
    w(:) = 0.0_dp

    ! CHECKS
    IF ((lon < -180.0_DP) .OR. (lon > 360.0_DP)) THEN
       status = 1  ! LONGITUDE OUT OF RANGE
       RETURN
    ELSE
       ! CHANGE INTERVAL [-180 ... 180] -> [0 ... 360]
       IF (lon <= 0.0_DP) THEN
          zlon = lon + 360.0_DP
       ELSE
          zlon = lon
       END IF
    END IF
    !
    IF ((lat < -90.0_DP) .OR. (lat > 90.0_DP)) THEN
       status = 2  ! LATITUDE OUT OF RANGE
       RETURN
    END IF

    ! INITIALIZATION
    IF (.NOT. linit) THEN
       ! FIELDS NOT YET INITIALIZED
       ALLOCATE(lpe  (nproma, ngpblks))
       ALLOCATE(ljp  (nproma, ngpblks))
       ALLOCATE(ljrow(nproma, ngpblks))
       lpe(:,:) = REAL(p_pe, DP)
       DO zjrow=1, ngpblks
          DO zkproma=1, nproma
             ljp  (zkproma,zjrow) = REAL(zkproma, DP)
             ljrow(zkproma,zjrow) = REAL(zjrow,   DP)
          END DO
       END DO
       !
       CALL trp_gpdc_gpgl(1, lpe,   gpe)
       CALL trp_gpdc_gpgl(1, ljp,   gjp)
       CALL trp_gpdc_gpgl(1, ljrow, gjrow)
       !
       DEALLOCATE(lpe) ; NULLIFY(lpe)
       DEALLOCATE(ljp) ; NULLIFY(ljp)
       DEALLOCATE(ljrow) ; NULLIFY(ljrow)
       !
       linit = .TRUE.
       !
    END IF

    ! SET INDICES
    CALL nn_index(philon, zlon, jgx1, jgx2)
    CALL nn_index(philat,  lat, jgy1, jgy2)

    v(1) = zlon
    v(2) = lat

    jgx = MIN(jgx1, jgx2)
    jgy = MIN(jgy1, jgy2)
    vn(1,1) = philon(jgx)
    vn(1,2) = philat(jgy)
    pe(1)     = NINT(gpe  (jgx,jgy))
    jp(1)     = NINT(gjp  (jgx,jgy))
    jrow(1)   = NINT(gjrow(jgx,jgy))

    IF (PRESENT(gjx)) gjx(1) = jgx
    IF (PRESENT(gjy)) gjy(1) = jgy

    jgx = MAX(jgx1, jgx2)
    jgy = MIN(jgy1, jgy2)
    vn(2,1) = philon(jgx)
    vn(2,2) = philat(jgy)
    pe(2)     = NINT(gpe  (jgx,jgy))
    jp(2)     = NINT(gjp  (jgx,jgy))
    jrow(2)   = NINT(gjrow(jgx,jgy))

    IF (PRESENT(gjx)) gjx(2) = jgx
    IF (PRESENT(gjy)) gjy(2) = jgy

    jgx = MAX(jgx1, jgx2)
    jgy = MAX(jgy1, jgy2)
    vn(3,1) = philon(jgx)
    vn(3,2) = philat(jgy)
    pe(3)     = NINT(gpe  (jgx,jgy))
    jp(3)     = NINT(gjp  (jgx,jgy))
    jrow(3)   = NINT(gjrow(jgx,jgy))

    IF (PRESENT(gjx)) gjx(3) = jgx
    IF (PRESENT(gjy)) gjy(3) = jgy

    jgx = MIN(jgx1, jgx2)
    jgy = MAX(jgy1, jgy2)
    vn(4,1) = philon(jgx)
    vn(4,2) = philat(jgy)
    pe(4)     = NINT(gpe  (jgx,jgy))
    jp(4)     = NINT(gjp  (jgx,jgy))
    jrow(4)   = NINT(gjrow(jgx,jgy))

    IF (PRESENT(gjx)) gjx(4) = jgx
    IF (PRESENT(gjy)) gjy(4) = jgy

    CALL bilin_weight(vn, v, w)

    status = 0

  END SUBROUTINE locate_in_decomp_4
  !-------------------------------------------------------------------------
! ##########################################################################
! END #ifndef MBM_CLAMS
! ##########################################################################
#endif
#endif

! ##########################################################################
! END #if defined(ECHAM5) || defined(MBM_CLAMS)
! ##########################################################################
#endif

! ============================================================================
! ============================================================================
! ============================================================================
! ============================================================================
! ============================================================================
#ifdef COSMO

  USE messy_main_constants_mem,  ONLY: dp, pi
  USE messy_main_grid_trafo,     ONLY: phi2phirot, rla2rlarot

  IMPLICIT NONE

  INTERFACE locate_in_decomp
     MODULE PROCEDURE locate_in_decomp_nrst
     MODULE PROCEDURE locate_in_decomp_4
  END INTERFACE
  PUBLIC :: locate_in_decomp

CONTAINS

  !--------------------------------------------------------------------------
  SUBROUTINE locate_in_decomp_nrst(status, lon, lat, pe, jp, jrow, jgx, jgy &
       , lrot, linclbnd, lprinterr, grlon, grlat)

    USE messy_main_mpi_bi,          ONLY: p_nprocs, isubpos, nboundlines &
                                        , nprocy
    USE messy_main_grid_def_mem_bi, ONLY: startlon_tot, startlat_tot &
                                        , dlon, dlat, ie_tot, je_tot &
                                        , polgam, pollon, pollat
    USE messy_main_constants_mem,   ONLY: iouerr

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    REAL(DP), INTENT(IN)  :: lon     ! LONGITUDE
    REAL(DP), INTENT(IN)  :: lat     ! LATITUDE
    INTEGER,  INTENT(OUT) :: pe      ! PROCESS ID
    INTEGER,  INTENT(OUT) :: jp      ! COLUMN
    INTEGER,  INTENT(OUT) :: jrow    ! ROW
    INTEGER,  INTENT(OUT), OPTIONAL :: jgx  ! lon index in global field
    INTEGER,  INTENT(OUT), OPTIONAL :: jgy  ! lat index in global field
    LOGICAL,  INTENT(IN),  OPTIONAL :: lrot ! .TRUE. if incoming point already
                                            ! on rotated grid
    LOGICAL,  INTENT(IN),  OPTIONAL :: linclbnd !.TRUE. if boundary should be
                                                ! included into location search
    LOGICAL,  INTENT(IN),  OPTIONAL :: lprinterr ! .TRUE. for error output
    REAL(DP), INTENT(OUT), OPTIONAL :: grlon     ! LONGITUDE
    REAL(DP), INTENT(OUT), OPTIONAL :: grlat     ! LATITUDE

    ! LOCAL
    INTEGER                               :: i,j,ix,jy,p
    LOGICAL                               :: lrotated
    LOGICAL,  SAVE                        :: linit =.FALSE.
    REAL(DP), DIMENSION(:), POINTER, SAVE :: philat => NULL()
    REAL(DP), DIMENSION(:), POINTER, SAVE :: philon => NULL()
    INTEGER                               :: istartlon
    INTEGER                               :: istartlat
    INTEGER                               :: idlon, idlat
    REAL(dp), PARAMETER                   :: intfac = 100000._dp
    REAL(DP)                              :: zlon  ! LONGITUDE
    ! POINT in Rotated grid:
    REAL(DP)                              :: rlon     ! LONGITUDE
    REAL(DP)                              :: rlat     ! LATITUDE

    ! INITIALISE
    status = 0
    IF (Present(jgx))  jgx = -1
    IF (Present(jgy))  jgy = -1
    PE   = -1
    jp   = -1
    jrow = -1

    ! CHECKS
    IF ((lon < -180.0_DP) .OR. (lon > 360.0_DP)) THEN
       status = 1  ! LONGITUDE OUT OF RANGE
       RETURN
    ELSE
       ! CHANGE INTERVAL [0 ... 360] -> [-180 ... 180]
       IF (lon > 180.0_DP) THEN
          zlon = lon - 360.0_DP
       ELSE
          zlon = lon
       END IF
    END IF
    !
    IF ((lat < -90.0_DP) .OR. (lat > 90.0_DP)) THEN
       status = 2  ! LATITUDE OUT OF RANGE
       RETURN
    END IF

    IF (PRESENT(lrot)) THEN
       lrotated = lrot
    ELSE
       lrotated = .FALSE.
    ENDIF

    IF (lrotated) THEN
       rlat = lat
       rlon = zlon
    ELSE
       ! Transform POINT into rotated grid
       rlat =   phi2phirot ( lat , zlon , pollat, pollon)
       rlon =   rla2rlarot ( lat , zlon , pollat, pollon, polgam)
    ENDIF

    IF (.NOT. linit) THEN
       istartlon = INT (startlon_tot * intfac)
       istartlat = INT (startlat_tot * intfac)
       idlon = INT(dlon * intfac)
       idlat = INT(dlat * intfac)

       ! NOTE: startlon_tot, startlat_tot are the coordinates of the lower left
       !       gridpoint mids. i.e. we need to shift it by 0.5 dlon or dlat
       !       to get the interfaces of the grid boxes to check, if a
       !       coordinate pair is located in one specific grid box
       ! Longitude (rotated)
       ALLOCATE(philon(ie_tot+1))
       DO i=1, ie_tot+1
          philon(i) = (REAL(istartlon+(i-1)*idlon,dp)-0.5_dp * &
               REAL(idlon,dp))/intfac
       END DO
       ! Latitude (rotated)
       ALLOCATE(philat(je_tot+1))
       DO j=1, je_tot+1
          philat(j) = (REAL(istartlat+(j-1)*idlat,dp)-0.5_dp * &
               REAL(idlat,dp))/intfac
       END Do
       !
       linit = .TRUE.
       !
    END IF

    ! find global longitude
    ix = -1
    DO i=1,ie_tot
       IF (philon(i) <= rlon .AND. philon(i+1) >= rlon)  THEN
          ix = i
          EXIT
       END IF
    END DO

    ! find global latitude
    jy = -1
    do j=1,je_tot
       IF (philat(j) <= rlat .AND. philat(j+1) >= rlat)  THEN
          jy = j
          EXIT
       END IF
    END DO

    IF (PRESENT(jgx)) jgx = ix
    IF (PRESENT(jgy)) jgy = jy

    ix_jy_undef: IF (jy /= -1 .AND. ix /= -1) THEN

       ! find PE and local indices
       DO p = p_nprocs-1, 0, -1
          ! NOTE: isubpos must not be <=1 but always >1 (due to halo)
          IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
             pe   = p
             jp   = ix - isubpos(p,1) + nboundlines + 1
             jrow = jy - isubpos(p,2) + nboundlines + 1
             EXIT
          END IF
       END DO

       ! geographical point not located in inner model domain
       IF ( (pe == -1) .OR. (jp == -1) .or. (jrow == -1) )   THEN
          if_pres: IF (PRESENT(linclbnd)) THEN
             if_lincb: IF (linclbnd) THEN
                ! 0 < ix <  isubpos(p,1)
                IF (ix <= isubpos(0,1)) THEN
                   ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
                   IF (jy <= isubpos(0,2)) THEN
                      pe   = 0
                      jp   = ix-isubpos(0,1) + nboundlines + 1
                      jrow = jy-isubpos(0,2) + nboundlines + 1
                   ELSE
                      ! .AND. jy > isubpos(p,2)
                      ! => search Pe on left domain bound.
                      DO p=nprocy-1,0,-1
                         IF (jy >= isubpos(p,2)) THEN
                            pe   = p
                            jp   = ix-isubpos(p,1) + nboundlines + 1
                            jrow = jy-isubpos(p,2) + nboundlines + 1
                            EXIT
                         END IF
                      END DO
                   END IF
                ELSE IF (jy <= isubpos(0,2)) THEN
                   ! case ix <= isubpos(p,1) already done
                   DO p=p_nprocs-nprocy,0,-4
                      IF (ix >= isubpos(p,1)) THEN
                         pe   = p
                         jp   = ix-isubpos(p,1) + nboundlines + 1
                         jrow = jy-isubpos(p,2) + nboundlines + 1
                         EXIT
                      END IF
                   END DO
                END IF
             END IF if_lincb
          END IF if_pres

       ENDIF

    END IF ix_jy_undef

    ! geographical point not located in model domain
    IF ( (pe == -1) .OR. (jp == -1) .or. (jrow == -1) ) THEN
       IF (PRESENT(lprinterr)) THEN
          IF (lprinterr) &
               write(iouerr,*) 'locate_in_decomp ' &
               , lrotated, lon,lat, rlon, rlat &
               , ' IX ',ix,jy, pe, jp, jrow &
               , 'BOUND', philon(1),philon(ie_tot),philat(1),philat(je_tot)
       ENDIF
       status = 555
       RETURN
    END IF

    IF (PRESENT(grlon)) grlon = philon(ix)+0.5_dp*dlon
    IF (PRESENT(grlat)) grlat = philat(jy)+0.5_dp*dlat

  END SUBROUTINE locate_in_decomp_nrst
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE locate_in_decomp_4(status, lon, lat, pe, jp, jrow, w, gjx, gjy &
       , lrot, linclbnd, lprinterr)

    ! LOCATES GEOGRAPHICAL POSITION IN VECTOR-PARALLEL-DECOMPOSITION
    !
    ! Author: Patrick Joeckel, DLR, May 2011

    USE messy_main_mpi_bi,       ONLY: p_pe, p_nprocs, isubpos, nboundlines &
                                     , nprocy
    USE messy_main_grid_def_mem_bi, ONLY: startlon_tot, startlat_tot &
                                     , dlon, dlat, ie_tot, je_tot &
                                     , polgam, pollon, pollat
    USE messy_main_tools,        ONLY: nn_index, bilin_weight

    IMPLICIT NONE
    INTRINSIC :: PRESENT, REAL, MAX, MIN

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    REAL(DP),               INTENT(IN)  :: lon     ! LONGITUDE [-180 ... 180]
    REAL(DP),               INTENT(IN)  :: lat     ! LATITUDE  [-90 ... 90]
    INTEGER,  DIMENSION(4), INTENT(OUT) :: pe      ! PROCESS ID
    INTEGER,  DIMENSION(4), INTENT(OUT) :: jp      ! COLUMN
    INTEGER,  DIMENSION(4), INTENT(OUT) :: jrow    ! ROW
    REAL(DP), DIMENSION(4), INTENT(OUT) :: w       ! WEIGHT
    ! lon index in global field
    INTEGER,  DIMENSION(4), INTENT(OUT), OPTIONAL :: gjx
    ! lat index in global field
    INTEGER,  DIMENSION(4), INTENT(OUT), OPTIONAL :: gjy
    ! .TRUE. if incoming point already on rotated grid
    LOGICAL,                INTENT(IN),  OPTIONAL :: lrot
    !.TRUE. if boundary should be included into location search
    LOGICAL,                INTENT(IN),  OPTIONAL :: linclbnd
    ! .TRUE. for error output
    LOGICAL,                INTENT(IN),  OPTIONAL :: lprinterr

    ! LOCAL
    LOGICAL,  SAVE                        :: linit = .FALSE.
    REAL(DP), DIMENSION(:), POINTER, SAVE :: philon => NULL()
    REAL(DP), DIMENSION(:), POINTER, SAVE :: philat => NULL()
    REAL(DP)              :: zlon
    ! POINT in Rotated grid:
    REAL(DP)              :: rlon     ! LONGITUDE
    REAL(DP)              :: rlat     ! LATITUDE
    LOGICAL               :: lrotated
    !
    INTEGER                                 :: i, j, p
    INTEGER                                 :: jgx1, jgy1, jgx2, jgy2
    INTEGER                                 :: ix, jy
    REAL(DP), DIMENSION(4,2)                :: vn
    REAL(DP), DIMENSION(2)                  :: v
    INTEGER                                 :: istartlon
    INTEGER                                 :: istartlat
    INTEGER                                 :: idlon, idlat
    REAL(dp), PARAMETER                     :: intfac = 100000._dp

    ! INIT
    status = 0
    IF (Present(gjx))  gjx(:) = -1
    IF (Present(gjy))  gjy(:) = -1
    PE(:)   = -1
    jp(:)   = -1
    jrow(:) = -1
    w(:) = 0.0_dp

    ! CHECKS
    IF ((lon < -180.0_DP) .OR. (lon > 360.0_DP)) THEN
       status = 1  ! LONGITUDE OUT OF RANGE
       RETURN
    ELSE
       ! CHANGE INTERVAL [0 ... 360] -> [-180 ... 180] (for COSMO!!!)
       IF (lon >= 180.0_DP) THEN
          zlon = lon - 360.0_DP
       ELSE
          zlon = lon
       END IF
    END IF
    !
    IF ((lat < -90.0_DP) .OR. (lat > 90.0_DP)) THEN
       status = 2  ! LATITUDE OUT OF RANGE
       RETURN
    END IF

    IF (PRESENT(lrot)) THEN
       lrotated = lrot
    ELSE
       lrotated = .FALSE.
    ENDIF

    IF (lrotated) THEN
       rlat = lat
       rlon = zlon
    ELSE
       ! Transform POINT into rotated grid
       rlat =   phi2phirot ( lat , zlon , pollat, pollon)
       rlon =   rla2rlarot ( lat , zlon , pollat, pollon, polgam)
    ENDIF

    ! check, if point is in range
    IF ( (rlat < startlat_tot-0.5_dp*dlat) .OR. &
         (rlat > startlat_tot+(je_tot-0.5_dp)*dlat) .OR. &
         (rlon < startlon_tot-0.5_dp*dlon) .OR. &
         (rlon > startlon_tot+(ie_tot-0.5_dp)*dlon) ) THEN
       ! geographical point not located in model domain
       status = 555
       RETURN
    ENDIF

    IF (.NOT. linit) THEN
       istartlon = INT (startlon_tot * intfac)
       istartlat = INT (startlat_tot * intfac)
       idlon = INT(dlon * intfac)
       idlat = INT(dlat * intfac)

       ! Longitude (rotated)
       ALLOCATE(philon(ie_tot))
       DO i=1,ie_tot
          philon(i) = REAL(istartlon+(i-1)*idlon,dp) / intfac
       END DO
       ! Latitude (rotated)
       ALLOCATE(philat(je_tot))
       DO j=1,je_tot
          philat(j) = REAL(istartlat+(j-1)*idlat,dp) / intfac
       END Do
       !
       linit = .TRUE.
       !
    END IF

    ! SET INDICES
    CALL nn_index(philon, rlon, jgx1, jgx2, 2.0_dp, status)
    CALL nn_index(philat, rlat, jgy1, jgy2, 2.0_dp, status)
    ! OUT OF RANGE
    IF (status /= 0) THEN
       status = 555
       RETURN
    END IF

    v(1) = rlon
    v(2) = rlat

    ix = MIN(jgx1, jgx2)
    jy = MIN(jgy1, jgy2)
    vn(1,1) = philon(ix)
    vn(1,2) = philat(jy)
    DO p=p_nprocs-1,0,-1
       ! NOTE: isubpos must not be <=1, always >1 (due to halo)
       IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
          pe(1)   = p
          jp(1)   = ix-isubpos(p,1) + nboundlines +1
          jrow(1) = jy-isubpos(p,2) + nboundlines +1
          EXIT
       END IF
    END DO

    IF (PRESENT(gjx)) gjx(1) = ix
    IF (PRESENT(gjy)) gjy(1) = jy

    ! geographical point not located in inner model domain
    IF ( (pe(1)== -1) .OR. (jp(1) == -1) .or. (jrow(1) == -1) )   THEN
       if_pres: IF (PRESENT(linclbnd)) THEN
          if_lincb: IF (linclbnd) THEN
             ! 0 < ix <  isubpos(p,1)
             IF (ix <= isubpos(0,1)) THEN
                ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
                IF (jy <= isubpos(0,2)) THEN
                   pe(1)   = 0
                   jp(1)   = ix-isubpos(0,1) + nboundlines + 1
                   jrow(1) = jy-isubpos(0,2) + nboundlines + 1
                ELSE
                   ! .AND. jy > isubpos(p,2)
                   ! => search Pe on left domain bound.
                   DO p=nprocy-1,0,-1
                      IF (jy >= isubpos(p,2)) THEN
                         pe(1)   = p
                         jp(1)   = ix-isubpos(p,1) + nboundlines + 1
                         jrow(1) = jy-isubpos(p,2) + nboundlines + 1
                         EXIT
                      END IF
                   END DO
                END IF
             ELSE IF (jy <= isubpos(0,2)) THEN
                ! case ix <= isubpos(p,1) already done
                DO p=p_nprocs-nprocy,0,-4
                   IF (ix >= isubpos(p,1)) THEN
                      pe(1)   = p
                      jp(1)   = ix-isubpos(p,1) + nboundlines + 1
                      jrow(1) = jy-isubpos(p,2) + nboundlines + 1
                      EXIT
                   END IF
                END DO
             END IF
          END IF if_lincb
       END IF if_pres

    ENDIF

    ix = MAX(jgx1, jgx2)
    jy = MIN(jgy1, jgy2)
    vn(2,1) = philon(ix)
    vn(2,2) = philat(jy)
    DO p=p_nprocs-1,0,-1
       ! NOTE: isubpos must not be <=1, always >1 (due to halo)
       IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
          pe(2)   = p
          jp(2)   = ix-isubpos(p,1) + nboundlines +1
          jrow(2) = jy-isubpos(p,2) + nboundlines +1
          EXIT
       END IF
    END DO

    IF (PRESENT(gjx)) gjx(2) = ix
    IF (PRESENT(gjy)) gjy(2) = jy

    ! geographical point not located in inner model domain
    IF ( (pe(2)== -1) .OR. (jp(2) == -1) .or. (jrow(2) == -1) )   THEN
       if_pres2: IF (PRESENT(linclbnd)) THEN
          if_lincb2: IF (linclbnd) THEN
             ! 0 < ix <  isubpos(p,1)
             IF (ix <= isubpos(0,1)) THEN
                ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
                IF (jy <= isubpos(0,2)) THEN
                   pe(2)   = 0
                   jp(2)   = ix-isubpos(0,1) + nboundlines + 1
                   jrow(2) = jy-isubpos(0,2) + nboundlines + 1
                ELSE
                   ! .AND. jy > isubpos(p,2)
                   ! => search Pe on left domain bound.
                   DO p=nprocy-1,0,-1
                      IF (jy >= isubpos(p,2)) THEN
                         pe(2)   = p
                         jp(2)   = ix-isubpos(p,1) + nboundlines + 1
                         jrow(2) = jy-isubpos(p,2) + nboundlines + 1
                         EXIT
                      END IF
                   END DO
                END IF
             ELSE IF (jy <= isubpos(0,2)) THEN
                ! case ix <= isubpos(p,1) already done
                DO p=p_nprocs-nprocy,0,-4
                   IF (ix >= isubpos(p,1)) THEN
                      pe(2)   = p
                      jp(2)   = ix-isubpos(p,1) + nboundlines + 1
                      jrow(2) = jy-isubpos(p,2) + nboundlines + 1
                      EXIT
                   END IF
                END DO
             END IF
          END IF if_lincb2
       END IF if_pres2

    ENDIF

    ix = MAX(jgx1, jgx2)
    jy = MAX(jgy1, jgy2)
    vn(3,1) = philon(ix)
    vn(3,2) = philat(jy)
    DO p=p_nprocs-1,0,-1
       ! NOTE: isubpos must not be <=1, always >1 (due to halo)
       IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
          pe(3)   = p
          jp(3)   = ix-isubpos(p,1) + nboundlines +1
          jrow(3) = jy-isubpos(p,2) + nboundlines +1
          EXIT
       END IF
    END DO

    IF (PRESENT(gjx)) gjx(3) = ix
    IF (PRESENT(gjy)) gjy(3) = jy

    ! geographical point not located in inner model domain
    IF ( (pe(3)== -1) .OR. (jp(3) == -1) .or. (jrow(3) == -1) )   THEN
       if_pres3: IF (PRESENT(linclbnd)) THEN
          if_lincb3: IF (linclbnd) THEN
             ! 0 < ix <  isubpos(p,1)
             IF (ix <= isubpos(0,1)) THEN
                ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
                IF (jy <= isubpos(0,2)) THEN
                   pe(3)   = 0
                   jp(3)   = ix-isubpos(0,1) + nboundlines + 1
                   jrow(3) = jy-isubpos(0,2) + nboundlines + 1
                ELSE
                   ! .AND. jy > isubpos(p,2)
                   ! => search Pe on left domain bound.
                   DO p=nprocy-1,0,-1
                      IF (jy >= isubpos(p,2)) THEN
                         pe(3)   = p
                         jp(3)   = ix-isubpos(p,1) + nboundlines + 1
                         jrow(3) = jy-isubpos(p,2) + nboundlines + 1
                         EXIT
                      END IF
                   END DO
                END IF
             ELSE IF (jy <= isubpos(0,2)) THEN
                ! case ix <= isubpos(p,1) already done
                DO p=p_nprocs-nprocy,0,-4
                   IF (ix >= isubpos(p,1)) THEN
                      pe(3)   = p
                      jp(3)   = ix-isubpos(p,1) + nboundlines + 1
                      jrow(3) = jy-isubpos(p,2) + nboundlines + 1
                      EXIT
                   END IF
                END DO
             END IF
          END IF if_lincb3
       END IF if_pres3

    ENDIF

    ix = MIN(jgx1, jgx2)
    jy = MAX(jgy1, jgy2)
    vn(4,1) = philon(ix)
    vn(4,2) = philat(jy)
    DO p=p_nprocs-1,0,-1
       ! NOTE: isubpos must not be <=1, always >1 (due to halo)
       IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
          pe(4)   = p
          jp(4)   = ix-isubpos(p,1) + nboundlines +1
          jrow(4) = jy-isubpos(p,2) + nboundlines +1
          EXIT
       END IF
    END DO

    IF (PRESENT(gjx)) gjx(4) = ix
    IF (PRESENT(gjy)) gjy(4) = jy

    ! geographical point not located in inner model domain
    IF ( (pe(4)== -1) .OR. (jp(4) == -1) .or. (jrow(4) == -1) )   THEN
       if_pres4: IF (PRESENT(linclbnd)) THEN
          if_lincb4: IF (linclbnd) THEN
             ! 0 < ix <  isubpos(p,1)
             IF (ix <= isubpos(0,1)) THEN
                ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
                IF (jy <= isubpos(0,2)) THEN
                   pe(4)   = 0
                   jp(4)   = ix-isubpos(0,1) + nboundlines + 1
                   jrow(4) = jy-isubpos(0,2) + nboundlines + 1
                ELSE
                   ! .AND. jy > isubpos(p,2)
                   ! => search Pe on left domain bound.
                   DO p=nprocy-1,0,-1
                      IF (jy >= isubpos(p,2)) THEN
                         pe(4)   = p
                         jp(4)   = ix-isubpos(p,1) + nboundlines + 1
                         jrow(4) = jy-isubpos(p,2) + nboundlines + 1
                         EXIT
                      END IF
                   END DO
                END IF
             ELSE IF (jy <= isubpos(0,2)) THEN
                ! case ix <= isubpos(p,1) already done
                DO p=p_nprocs-nprocy,0,-4
                   IF (ix >= isubpos(p,1)) THEN
                      pe(4)   = p
                      jp(4)   = ix-isubpos(p,1) + nboundlines + 1
                      jrow(4) = jy-isubpos(p,2) + nboundlines + 1
                      EXIT
                   END IF
                END DO
             END IF
          END IF if_lincb4
       END IF if_pres4

    ENDIF

    CALL bilin_weight(vn, v, w)

    status = 0

  END SUBROUTINE locate_in_decomp_4
  !-------------------------------------------------------------------------

! ##########################################################################
! END #if defined(COSMO)
! ##########################################################################
#endif

#if defined(CESM1)

  use messy_main_constants_mem, only: dp

  ! FLAGS FOR trp_gpdc_gpgl(-1, lf, gf, FLAG)
  INTEGER, PARAMETER, PUBLIC :: M_SUM = 0  ! SUM OVER ALL PEs
  INTEGER, PARAMETER, PUBLIC :: M_AVE = 1  ! AVERAGE OVER ALL PEs
  INTEGER, PARAMETER, PUBLIC :: M_STD = 2  ! STANDARD DEV. BETWEEN ALL PEs
  INTEGER, PARAMETER, PUBLIC :: M_LOC = 3  ! EXTRACT ONLY LOCAL

  INTERFACE locate_in_decomp
     MODULE PROCEDURE locate_in_decomp_nrst
     MODULE PROCEDURE locate_in_decomp_4
  END INTERFACE
  PUBLIC :: locate_in_decomp

  INTERFACE trp_gpdc_gpgl
     MODULE PROCEDURE trp_gpdc_gpgl_4d
     MODULE PROCEDURE trp_gpdc_gpgl_3d
     MODULE PROCEDURE trp_gpdc_gpgl_2d
  END INTERFACE
  PUBLIC :: trp_gpdc_gpgl  ! decomposed gridpoint space  <->
                           !                        global gp-space on each PE

contains

  ! only dummies
  SUBROUTINE locate_in_decomp_nrst(status, lon, lat, pe, jp, jrow, jgx, jgy &
       , lrot, linclbnd, lprinterr)

     IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    REAL(DP), INTENT(IN)  :: lon     ! LONGITUDE
    REAL(DP), INTENT(IN)  :: lat     ! LATITUDE
    INTEGER,  INTENT(OUT) :: pe      ! PROCESS ID
    INTEGER,  INTENT(OUT) :: jp      ! COLUMN
    INTEGER,  INTENT(OUT) :: jrow    ! ROW
    INTEGER,  INTENT(OUT), OPTIONAL :: jgx  ! lon index in global field
    INTEGER,  INTENT(OUT), OPTIONAL :: jgy  ! lat index in global field
    LOGICAL,  INTENT(IN),  OPTIONAL :: lrot ! .TRUE. if incoming point already
                                            ! on rotated grid
    LOGICAL,  INTENT(IN),  OPTIONAL :: linclbnd !.TRUE. if boundary should be
                                                ! included into location search
    LOGICAL,  INTENT(IN),  OPTIONAL :: lprinterr ! .TRUE. for error output

  end SUBROUTINE locate_in_decomp_nrst

  SUBROUTINE locate_in_decomp_4(status, lon, lat, pe, jp, jrow, w, gjx, gjy &
       , lrot, linclbnd, lprinterr)

    IMPLICIT NONE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    REAL(DP),               INTENT(IN)  :: lon     ! LONGITUDE [-180 ... 180]
    REAL(DP),               INTENT(IN)  :: lat     ! LATITUDE  [-90 ... 90]
    INTEGER,  DIMENSION(4), INTENT(OUT) :: pe      ! PROCESS ID
    INTEGER,  DIMENSION(4), INTENT(OUT) :: jp      ! COLUMN
    INTEGER,  DIMENSION(4), INTENT(OUT) :: jrow    ! ROW
    REAL(DP), DIMENSION(4), INTENT(OUT) :: w       ! WEIGHT
    ! lon index in global field
    INTEGER,  DIMENSION(4), INTENT(OUT), OPTIONAL :: gjx
    ! lat index in global field
    INTEGER,  DIMENSION(4), INTENT(OUT), OPTIONAL :: gjy
    ! .TRUE. if incoming point already on rotated grid
    LOGICAL,                INTENT(IN),  OPTIONAL :: lrot
    !.TRUE. if boundary should be included into location search
    LOGICAL,                INTENT(IN),  OPTIONAL :: linclbnd
    ! .TRUE. for error output
    LOGICAL,                INTENT(IN),  OPTIONAL :: lprinterr

  end SUBROUTINE locate_in_decomp_4


  SUBROUTINE trp_gpdc_gpgl_4d(sign, lf, gf, method)
    IMPLICIT NONE
    INTRINSIC :: ABS, ASSOCIATED, PRESENT, REAL, SIZE, SQRT
    INTEGER,                  INTENT(IN)           :: sign
    REAL(dp), DIMENSION(:,:,:,:), POINTER          :: lf    ! decomposed field
    REAL(dp), DIMENSION(:,:,:,:), POINTER          :: gf    ! global field
    INTEGER,                  INTENT(IN), OPTIONAL :: method
  END SUBROUTINE trp_gpdc_gpgl_4d

  SUBROUTINE trp_gpdc_gpgl_3d(sign, lf, gf, method)
    IMPLICIT NONE
    INTRINSIC :: ABS, ASSOCIATED, PRESENT, REAL, SIZE, SQRT
    INTEGER,                  INTENT(IN)           :: sign
    REAL(dp), DIMENSION(:,:,:), POINTER          :: lf    ! decomposed field
    REAL(dp), DIMENSION(:,:,:), POINTER          :: gf    ! global field
    INTEGER,                  INTENT(IN), OPTIONAL :: method
  END SUBROUTINE trp_gpdc_gpgl_3d

  SUBROUTINE trp_gpdc_gpgl_2d(sign, lf, gf, method)
    IMPLICIT NONE
    INTRINSIC :: ABS, ASSOCIATED, PRESENT, REAL, SIZE, SQRT
    INTEGER,                  INTENT(IN)           :: sign
    REAL(dp), DIMENSION(:,:), POINTER          :: lf    ! decomposed field
    REAL(dp), DIMENSION(:,:), POINTER          :: gf    ! global field
    INTEGER,                  INTENT(IN), OPTIONAL :: method
  END SUBROUTINE trp_gpdc_gpgl_2d

! ##########################################################################
! END #if defined(CESM1)
! ##########################################################################
#endif

! ##########################################################################
! #if defined(ICON)
! ##########################################################################
#if defined(ICON)
#define _MPI_PSUM

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather_bi, ONLY: error_bi
  USE messy_main_mpi_bi,     ONLY: dcg, p_io0, p_pe, p_nprocs &
     &                           , p_send, p_recv, p_sendrecv &
     &                           , p_bcast, p_comm_work, p_max

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  INTEGER, PARAMETER :: tag_gather_ix    = 320
  INTEGER, PARAMETER :: tag_scatter_ix   = 321
  PUBLIC :: get_dc_index   ! get decomposition along index
  PUBLIC :: gather_trix
  PUBLIC :: gather_trix_1d
  PUBLIC :: scatter_trix

  TYPE t_cartesian_coord
     REAL(DP) :: x(3)
  END TYPE t_cartesian_coord

  INTERFACE locate_in_decomp
     MODULE PROCEDURE locate_in_decomp_nrst
     MODULE PROCEDURE locate_in_decomp_stencil
  END INTERFACE
  PUBLIC :: locate_in_decomp

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE locate_in_decomp_nrst(status, lon, lat, pe, jp, jrow, glbidx)

    ! LOCATES GEOGRAPHICAL POSITION IN VECTOR-PARALLEL-DECOMPOSITION

    USE messy_main_mpi_bi,        ONLY: p_pe
    USE messy_main_constants_mem, ONLY: DTR_icon, PI_2
    USE messy_main_channel_mem,   ONLY: dom_current
    USE messy_main_bmluse_bi,     ONLY: p_patch, grf_bdywidth_c &
                                      , min_rlcell_int          &
                                      , idx_1d, get_indices_c

    IMPLICIT NONE
    INTRINSIC :: PRESENT, SIN, COS, ACOS

    ! I/O
    INTEGER,  INTENT(OUT)           :: status
    REAL(DP), INTENT(IN)            :: lon     ! LONGITUDE [-180 ... 180]
    REAL(DP), INTENT(IN)            :: lat     ! LATITUDE  [-90 ... 90]
    INTEGER,  INTENT(OUT)           :: pe      ! PROCESS ID
    INTEGER,  INTENT(OUT)           :: jp      ! COLUMN (index)
    INTEGER,  INTENT(OUT)           :: jrow    ! ROW    (block)
    INTEGER,  INTENT(OUT), OPTIONAL :: glbidx  ! index in global field

    ! LOCAL
    REAL(DP)                        :: zlon, zlat, plon, plat

    TYPE(t_cartesian_coord)         :: p_x
    TYPE(t_cartesian_coord)         :: c_x(3)
    REAL(DP)                        :: clon, clat                             &
                                     , p_sinlon, p_coslon, p_sinlat, p_coslat &
                                     , x_sinlon, x_coslon, x_sinlat, x_coslat &
                                     , px_coslon
    INTEGER                         :: v_idx, v_blk
    LOGICAL                         :: ccw1, ccw2, ccw3
    INTEGER                         :: i_startblk, i_endblk, i_startidx &
                                     , i_endidx, rl_start, rl_end
    INTEGER                         :: jc, jb, zglbidx
    INTEGER                         :: ii
    INTEGER                         :: host_pe

    ! INITIALISE
    status = 0

    IF (PRESENT(glbidx)) zglbidx = -1
    pe   = -1
    jp   = -1
    jrow = -1

    ! CHECKS
    IF ((lon < -180.0_DP) .OR. (lon > 360.0_DP)) THEN
       status = 1  ! LONGITUDE OUT OF RANGE
       RETURN
    ELSE
       ! grid coordinates in ICON in radians range from -180 .. 180 !!
       ! interval change has to be adapted for lon > 180 ... !!
       ! calculation of ccw is not affected, as handling with cartesian
       ! coordinates
       ! (transformed with trigonometric functions cos/sin)
       ! CHANGE INTERVAL [0 ... 360] -> [-180 ... 180]
       ! and transform to radians
       IF (lon > 180.0_dp) THEN
          zlon = (lon - 360.0_dp) * DTR_icon
       ELSE
          zlon = lon * DTR_icon
       END IF
    END IF
    !
    IF ((lat < -90.0_DP) .OR. (lat > 90.0_DP)) THEN
       status = 2  ! LATITUDE OUT OF RANGE
       RETURN
    END IF
    ! transform to radians
    zlat = lat * DTR_icon

    ! calculate cartesian coordinates from geographical coordinates
    p_sinlon = SIN(zlon)
    p_coslon = COS(zlon)
    p_sinlat = SIN(zlat)
    p_coslat = COS(zlat)
    p_x%x(1) = p_coslon*p_coslat
    p_x%x(2) = p_sinlon*p_coslat
    p_x%x(3) = p_sinlat

    ! loop over prognostic points on PE
    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch(dom_current)%cells%start_block(rl_start)
    i_endblk   = p_patch(dom_current)%cells%end_block(rl_end)

    DO jb = i_startblk, i_endblk
       CALL get_indices_c(p_patch(dom_current), jb, i_startblk, i_endblk &
            &           , i_startidx, i_endidx, rl_start, rl_end)
       DO jc = i_startidx, i_endidx
          ! get cell's coordinate
          plon = p_patch(dom_current)%cells%center(jc,jb)%lon
          plat = p_patch(dom_current)%cells%center(jc,jb)%lat
          ! for "high" latitudes (>45 degree), do a simple test,
          ! if point and cell are on same hemisphere, so avoid calculation
          ! of trigonometric functions in some cases
          IF (ABS(zlat) > (PI_2 * 0.5_dp)) THEN
             IF (SIGN(1._dp,zlat) /= SIGN(1._dp,plat)) CYCLE
          ELSE
             ! additional test, to avoid selecting cells on the
             ! direct opposite of the sphere
             ! use arc distance between the actual point and the centre
             ! of the cell, if it is larger than 90 degrees (PI/2),
             ! it is not the cell we want...
             x_sinlon = SIN(plon)
             x_coslon = COS(plon)
             x_sinlat = SIN(plat)
             x_coslat = COS(plat)
             px_coslon = COS(zlon - plon)
             IF (ACOS(p_sinlat * x_sinlat + p_coslat * x_coslat * px_coslon) &
                  > PI_2) CYCLE
          END IF
          ! get the coordinates of the three points sourrounding the
          ! current cell
          DO ii=1,3
             v_idx  = p_patch(dom_current)%cells%vertex_idx(jc,jb,ii)
             v_blk  = p_patch(dom_current)%cells%vertex_blk(jc,jb,ii)
             clon = p_patch(dom_current)%verts%vertex(v_idx,v_blk)%lon
             clat = p_patch(dom_current)%verts%vertex(v_idx,v_blk)%lat
             ! calculate cartesian coordinates from geographical coordinates
             x_sinlon = SIN(clon)
             x_coslon = COS(clon)
             x_sinlat = SIN(clat)
             x_coslat = COS(clat)
             c_x(ii)%x(1) = x_coslon*x_coslat
             c_x(ii)%x(2) = x_sinlon*x_coslat
             c_x(ii)%x(3) = x_sinlat
          END DO

          ! test the combination of two vertices and actual point
          ! form a triangle with vectors going counter-clockwise
          ccw1 = ccw(c_x(1),c_x(2),p_x)
          ccw2 = ccw(c_x(2),c_x(3),p_x)
          ccw3 = ccw(c_x(3),c_x(1),p_x)

          ! if all combinations result in same direction of vectors,
          ! the point is inside
          IF ( ((      ccw1) .AND. (      ccw2) .AND. (      ccw3)) .OR. &
             & ((.NOT. ccw1) .AND. (.NOT. ccw2) .AND. (.NOT. ccw3)) ) THEN
             pe     = p_pe
             jp     = jc
             jrow   = jb
             IF (PRESENT(glbidx)) &
               glbidx = &
               p_patch(dom_current)%cells%decomp_info%glb_index(idx_1d(jc,jb))
             EXIT
          END IF
       END DO
    END DO

    host_pe = p_max(pe, p_comm_work)
    IF (host_pe < 0) THEN ! all are -1
       ! point is outside domain
       status = 555
       RETURN
    ELSE
       ! make sure that all PEs return the same information
       pe   = host_pe
       jp   = p_max(jp, p_comm_work)
       jrow = p_max(jrow, p_comm_work)
       IF (PRESENT(glbidx)) glbidx = p_max(glbidx, p_comm_work)
       status = 0
    END IF

  END SUBROUTINE locate_in_decomp_nrst
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE locate_in_decomp_stencil(status, lon, lat, pe, jp, jrow, w &
       , glbidx, stencil)

    ! LOCATES GEOGRAPHICAL POSITION IN VECTOR-PARALLEL-DECOMPOSITION

    USE messy_main_mpi_bi,        ONLY: p_pe
    USE messy_main_channel_mem,   ONLY: dom_current
    USE messy_main_constants_mem, ONLY: DTR_icon
    USE messy_main_tools,         ONLY: inv_dist_weight
    USE messy_main_bmluse_bi,     ONLY: p_patch, idx_1d

    IMPLICIT NONE
    INTRINSIC :: PRESENT, SIZE

    ! I/O
    INTEGER,                INTENT(OUT)           :: status
    REAL(DP),               INTENT(IN)            :: lon     ! LONGITUDE [-180 ... 180]
    REAL(DP),               INTENT(IN)            :: lat     ! LATITUDE  [-90 ... 90]
    INTEGER,  DIMENSION(:), INTENT(OUT)           :: pe      ! PROCESS ID
    INTEGER,  DIMENSION(:), INTENT(OUT)           :: jp      ! COLUMN (index)
    INTEGER,  DIMENSION(:), INTENT(OUT)           :: jrow    ! ROW    (block)
    REAL(DP), DIMENSION(:), INTENT(OUT)           :: w       ! WEIGHT
    INTEGER,  DIMENSION(:), INTENT(OUT), OPTIONAL :: glbidx ! index in global field
    INTEGER,                 INTENT(IN), OPTIONAL :: stencil

    ! LOCAL
    INTEGER, DIMENSION(3)                         :: v_blk, v_idx
    INTEGER, DIMENSION(3,6)                       :: c_blk, c_idx
    REAL(DP), DIMENSION(:,:), ALLOCATABLE         :: neighbours
    REAL(DP), DIMENSION(2)                        :: position
    INTEGER                                       :: zstencil, zglbidx
    INTEGER                                       :: ii, jj, kk

    ! INITIALISE
    status = 0
    IF (PRESENT(glbidx))  glbidx(:) = -1
    pe(:)   = -1
    jp(:)   = -1
    jrow(:) = -1
    w(:)    = 0.0_dp

    IF (PRESENT(stencil)) THEN
       zstencil = stencil
    ELSE
       zstencil = SIZE(jp)
    END IF
    IF ((zstencil /= 4) .AND. (zstencil /= 13)) THEN
       status = 3 ! STENCIL HAS TO BE OF {4,13}
       RETURN
    END IF

    CALL locate_in_decomp_nrst(status, lon, lat, pe(1), jp(1), jrow(1), zglbidx)
    IF (status /= 0) RETURN ! 1, 2, 555

    IF (PRESENT(glbidx)) glbidx(1) = zglbidx

    ! find the neighbours and fill the pe, jp, jrow arrays
    ! at least the halo cells are on this pe
    pe(:) = pe(1)

    only_on_valid_pe: IF (p_pe == pe(1)) THEN

       ! the current triangle is already added to the stencil,
       ! so we go on with position 2
       kk = 2
       ! 1. add the 3 direct neighbors to the stencil
       DO ii = 1, 3
          ! get line and block indices of direct neighbors
          jp(kk)     = p_patch(dom_current)%cells%neighbor_idx(jp(1),jrow(1),ii)
          jrow(kk)   = p_patch(dom_current)%cells%neighbor_blk(jp(1),jrow(1),ii)
          IF (PRESENT(glbidx)) &
               glbidx(kk) = &
               p_patch(dom_current)%cells%decomp_info%glb_index( &
               idx_1d(jp(kk),jrow(kk)) )
          kk = kk + 1
       ENDDO

       ! 2. loop over the vertices and add all the cells
       !    that are no direct neighbors and not our Control volume.
       !    This should only be done, if we use a 13-stencil,
       !    otherwise we found everything we need!
       IF (zstencil == 13) THEN
          ! get line and block indices of cell vertices
          v_idx(1:3) = p_patch(dom_current)%cells%vertex_idx(jp(1),jrow(1),1:3)
          v_blk(1:3) = p_patch(dom_current)%cells%vertex_blk(jp(1),jrow(1),1:3)
          ! for each vertex: get all the cells which share this vertex
          DO jj = 1,3
             c_idx(jj,:) = &
                  p_patch(dom_current)%verts%cell_idx(v_idx(jj),v_blk(jj),:)
             c_blk(jj,:) = &
                  p_patch(dom_current)%verts%cell_blk(v_idx(jj),v_blk(jj),:)
          ENDDO

          DO ii = 1,3   ! loop over vertices
             ! loop over cells around each vertex
             DO jj = 1,p_patch(dom_current)%verts%num_edges(v_idx(ii),v_blk(ii))

                IF (.NOT.(   (c_idx(ii,jj) == jp(1) .AND. c_blk(ii,jj) == jrow(1)) &
                     &  .OR. (c_idx(ii,jj) == jp(2) .AND. c_blk(ii,jj) == jrow(2)) &
                     &  .OR. (c_idx(ii,jj) == jp(3) .AND. c_blk(ii,jj) == jrow(3)) &
                     &  .OR. (c_idx(ii,jj) == jp(4) .AND. c_blk(ii,jj) == jrow(4)) &
                     &  .OR. (c_idx(ii,jj) == 0     .AND. c_blk(ii,jj) == 0      ) ) ) THEN

                   jp(kk)   = c_idx(ii,jj)
                   jrow(kk) = c_blk(ii,jj)
                   IF (PRESENT(glbidx)) &
                        & glbidx(kk) = &
                        p_patch(dom_current)%cells%decomp_info%glb_index( &
                        idx_1d(jp(kk),jrow(kk)) )
                   kk = kk + 1
                ENDIF
             ENDDO
          ENDDO
       ENDIF ! (zstencil == 13)

       ! return number of actual available cells in this stencil
       ! at the given point
       ! (could be one less for 13 stencil, when on pentagon)
       zstencil = kk - 1

       ALLOCATE(neighbours(zstencil,2))
       DO kk = 1, zstencil
          neighbours(kk,1) = p_patch(dom_current)%cells%center( &
               jp(kk),jrow(kk) )%lon
          neighbours(kk,2) = p_patch(dom_current)%cells%center( &
               jp(kk),jrow(kk) )%lat
       END DO

       ! CHANGE INTERVAL [-180 ... 180] -> [0 ... 360]
       ! and transform to radians
       IF (lon < 0.0_DP) THEN
          position(1) = (lon + 360.0_DP) * DTR_icon
       ELSE
          position(1) = lon * DTR_icon
       END IF
       !
       ! transform to radians
       position(2) = lat * DTR_icon

       CALL inv_dist_weight(neighbours,position,w,status=status)
       DEALLOCATE(neighbours)

    END IF only_on_valid_pe
    CALL p_bcast(status, pe(1), p_comm_work)
    IF (status /= 0) RETURN

    ! make sure that all PEs return the same information
    jp(:)   = p_max(jp(:),   p_comm_work)
    jrow(:) = p_max(jrow(:), p_comm_work)
    w(:)    = p_max(w(:),    comm=p_comm_work)
    IF (PRESENT(glbidx)) glbidx(:) = p_max(glbidx(:), p_comm_work)

  END SUBROUTINE locate_in_decomp_stencil
  !-------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  ! cf. FUNCTION ccw_spherical in shr_horizontal/mo_delauny_types.f90
  ! The ICON implementation works on TYPE t_point, to avoid using this here,
  ! this is a slightly modified copy, working with a simplier TYPE, which should
  ! meet our needs.
  !> Locates a point relative to a directed arc. Let v1, v2, and v3 be
  !  distinct points on the sphere, and denote by v1 -> v2 the
  !  geodesic connecting v1 and v2 and directed toward v2. This test
  !  determines which of the two hemispheres defined by v1 -> v2
  !  contains v3, or, in other words, if we "turn left" when going
  !  from v1 to v2 to v3.
  ! -------------------------------------------------------------------------
  PURE FUNCTION ccw(v1,v2,v3)
     LOGICAL :: ccw
     TYPE(t_cartesian_coord), INTENT(IN) :: v1, v2, v3
     REAL(DP) ccw_r

    ! det(v1,v2,v3) = <v1 x v2, v3> = | v1 x v2 | cos(a)
    !
    ! where a is the angle between v3 and the normal to the plane
    ! defined by v1 and v2.
     ccw_r = v3%x(1)*(v1%x(2)*v2%x(3) - v2%x(2)*v1%x(3)) &
           - v3%x(2)*(v1%x(1)*v2%x(3) - v2%x(1)*v1%x(3)) &
           + v3%x(3)*(v1%x(1)*v2%x(2) - v2%x(1)*v1%x(2))

    ! we apply a static error of
    !   | e - e'| <= 3*2^-48
    ! to decide if a floating-point evaluation e' of an expression e
    ! has the correct sign, see Section 2.2 of
    !
    ! Burnikel, C.; Funke, S. & Seel, M.
    ! "Exact geometric computation using Cascading"
    ! International Journal of Computational Geometry & Applications,
    ! World Scientific, 2001, 11, 245-266
     IF (ABS(ccw_r) <= 1.1e-14_dp) THEN
        ccw = ccw_q128(v1,v2,v3)
     ELSE
        ccw = ccw_r <= 0._dp
     END IF

  ! ------------------------------------------------------------------------
  END FUNCTION ccw
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  !> Locates a point relative to a directed arc. Let v1, v2, and v3 be
  !  distinct points on the sphere, and denote by v1 -> v2 the
  !  geodesic connecting v1 and v2 and directed toward v2. This test
  !  determines which of the two hemispheres defined by v1 -> v2
  !  contains v3, or, in other words, if we "turn left" when going
  !  from v1 to v2 to v3.
  PURE FUNCTION ccw_q128(v1,v2,v3)
    INTEGER, PARAMETER :: QR_K = SELECTED_REAL_KIND (2*precision(1.0_dp))
    LOGICAL :: ccw_q128
    TYPE (t_cartesian_coord), INTENT(IN)  :: v1,v2,v3
    REAL(QR_K) :: v1_x, v1_y, v1_z,v2_x, v2_y, v2_z,v3_x, v3_y, v3_z

    ! det(v1,v2,v3) = <v1 x v2, v3> = | v1 x v2 | cos(a)
    !
    ! where a is the angle between v3 and the normal to the plane
    ! defined by v1 and v2.

    v1_x = v1%x(1)
    v1_y = v1%x(2)
    v1_z = v1%x(3)
    v2_x = v2%x(1)
    v2_y = v2%x(2)
    v2_z = v2%x(3)
    v3_x = v3%x(1)
    v3_y = v3%x(2)
    v3_z = v3%x(3)

    ccw_q128 = v3_x*(v1_y*v2_z - v2_y*v1_z) &
      &   -    v3_y*(v1_x*v2_z - v2_x*v1_z) &
      &   +    v3_z*(v1_x*v2_y - v2_x*v1_y)  <= 0._dp
  END FUNCTION ccw_q128
  !-------------------------------------------------------------------------

  ! =========================================================================
  SUBROUTINE gather_trix(gl, lc, index, xpe, XNG)

#ifdef _MPI_PSUM
    USE messy_main_mpi_bi, ONLY: p_sum
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRODUCT, SIZE, UBOUND

    ! I/O
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: gl     ! global field
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: lc     ! local field
    INTEGER, INTENT(IN)                   :: index  ! decomposition index
    INTEGER, INTENT(IN), OPTIONAL         :: xpe
    INTEGER, INTENT(IN), OPTIONAL         :: XNG

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr='gather_glix_4d'
    INTEGER, PARAMETER               :: nrank = 4
    INTEGER, DIMENSION(:,:), POINTER :: idx=>NULL()
    LOGICAL, DIMENSION(:),   POINTER :: ldo=>NULL()
    INTEGER                          :: NG, NL
    INTEGER                          :: NLs
    INTEGER                          :: pe, zxpe
    INTEGER                          :: i
    INTEGER                          :: ishpg(4) ! global shape
    INTEGER                          :: ishpl(4) ! local shape
    INTEGER                          :: ixg(4,2) ! global start/stop index
    INTEGER                          :: ixl(4,2) ! local start/stop index

    IF ((index < 1).OR.(index > nrank)) THEN
       CALL error_bi('index out of range', substr)
    END IF

    IF (.NOT. ASSOCIATED(lc)) THEN
       CALL error_bi('pointer to local field is not associated', substr)
    END IF

    IF (PRESENT(xpe)) THEN
       zxpe = xpe
    ELSE
       zxpe = p_io0
    END IF

    ! SHAPE OF 4-DARRAY
    DO i=1, nrank
       ishpl(i) = SIZE(lc,i)
    END DO

    IF (PRESENT(XNG)) THEN
       NG = XNG
    ELSE
       NL = ishpl(index)         ! size at decomp. index of local field

#ifndef _MPI_PSUM
       IF (p_pe == zxpe) THEN
          NG = NL        ! size at decomp. index of global field (summation)
       END IF
       !
       ! SUM OVER ALL NON-zxpe p_pe
       IF (p_pe /= zxpe) THEN
          CALL p_send (NL, zxpe, tag_gather_ix)
       ELSE ! p_pe = zxpe
          DO i=1, p_nprocs
             pe = i-1
             IF (pe /= p_pe) THEN
                CALL p_recv( NLs, pe, tag_gather_ix)
                NG = NG + NLs
             END IF
          END DO
       END IF
#else
       NG = p_sum(NL)
#endif
    ENDIF

    ! lcoal field start/stop indices given by size
    DO i=1, nrank
       ixl(i,1) = 1
       ixl(i,2) = SIZE(lc, i)
    END DO

#ifndef _MPI_PSUM
    CALL p_bcast(NG, zxpe)     ! global field dim. at decomp. index
#endif
    ishpg(:) = ishpl(:)        ! global shape is local shape ...
    ishpg(index) = NG          ! ... except at decomposition index

    ! global field start/stop indices
    ixg(:,:) = ixl(:,:)        ! global start/stop index is local start/stop
    ! index, except at decomposition index
    CALL get_dc_index(NG, idx, ldo)
    ! -> set below

    IF (p_pe == zxpe) THEN
       IF (ASSOCIATED(gl)) THEN
          IF (SIZE(gl) /= PRODUCT(ishpg)) THEN
             WRITE(*,*) 'SIZE(gl) = ',SIZE(gl),' =/= ',PRODUCT(ishpg)
             WRITE(*,*) 'SHAPE OF GLOBAL ARRAY: ',UBOUND(gl)
             WRITE(*,*) 'EXPECTED SHAPE       : ',ishpg(:)
             CALL error_bi('wrong size of global array', substr)
          END IF
       ELSE
          ALLOCATE(gl(ishpg(1),ishpg(2),ishpg(3),ishpg(4)))
       END IF
    ELSE
       IF (ASSOCIATED(gl)) DEALLOCATE(gl)
       NULLIFY(gl)
    ENDIF

    ! send if pe /= zxpe
    !
    IF (p_pe /= zxpe) THEN
       IF (ldo(p_pe)) CALL p_send (lc(ixl(1,1):ixl(1,2)  &
            ,  ixl(2,1):ixl(2,2)  &
            ,  ixl(3,1):ixl(3,2)  &
            ,  ixl(4,1):ixl(4,2)) &
            , zxpe, tag_gather_ix)
    ELSE ! p_pe = zxpe
       DO i=1, p_nprocs
          pe = i-1
          ! set global start/stop indices at decomp. index
          ixg(index,1) = idx(pe,1)
          ixg(index,2) = idx(pe,2)
          IF (pe /= p_pe) THEN
             IF (ldo(pe)) THEN
                CALL p_recv( gl(ixg(1,1):ixg(1,2)  &
                  ,   ixg(2,1):ixg(2,2)  &
                  ,   ixg(3,1):ixg(3,2)  &
                  ,   ixg(4,1):ixg(4,2)) &
                  , pe, tag_gather_ix)
             END IF
          ELSE
             IF (ldo(pe)) &
                  gl(ixg(1,1):ixg(1,2)     &
                  ,   ixg(2,1):ixg(2,2)    &
                  ,   ixg(3,1):ixg(3,2)    &
                  ,   ixg(4,1):ixg(4,2)) = &
                  lc(ixl(1,1):ixl(1,2)  &
                  ,  ixl(2,1):ixl(2,2)  &
                  ,  ixl(3,1):ixl(3,2)  &
                  ,  ixl(4,1):ixl(4,2))
          END IF
       END DO
    END IF

    ! CLEAN UP
    DEALLOCATE(idx) ; NULLIFY(idx)
    DEALLOCATE(ldo) ; NULLIFY(ldo)

  END SUBROUTINE gather_trix
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE scatter_trix(gl, lc, index, xpe, xishpg)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRODUCT, SIZE, UBOUND

    ! I/O
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: gl    ! global field
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: lc    ! local field
    INTEGER, INTENT(IN)                   :: index ! decomposition index
    INTEGER, INTENT(IN), OPTIONAL         :: xpe
    INTEGER, DIMENSION(4), INTENT(IN), OPTIONAL :: XISHPG

    ! LOCAL
    INTEGER, PARAMETER               :: nrank = 4
    CHARACTER(LEN=*), PARAMETER      :: substr='scatter_trix'
    INTEGER, DIMENSION(:,:), POINTER :: idx=>NULL()
    LOGICAL, DIMENSION(:),   POINTER :: ldo=>NULL()
    INTEGER                          :: NG
    INTEGER                          :: pe, zxpe
    INTEGER                          :: i
    INTEGER                          :: ishpg(4) ! global shape
    INTEGER                          :: ishpl(4) ! local shape
    INTEGER                          :: ixg(4,2) ! global start/stop index
    INTEGER                          :: ixl(4,2) ! local start/stop index

    IF ((index < 1).OR.(index > nrank)) THEN
       CALL error_bi('index out of range', substr)
    END IF

    IF (PRESENT(xpe)) THEN
       zxpe = xpe
    ELSE
       zxpe = p_io0
    END IF

    ! SHAPE OF 4-DARRAY
    IF (PRESENT(XISHPG)) THEN
       ishpg(:) = XISHPG(:)
    ELSE
       IF (p_pe == zxpe) THEN
          DO i=1, nrank
             ishpg(i) = SIZE(gl,i)
          END DO
       END IF
       CALL p_bcast(ishpg, zxpe)
    ENDIF
    NG = ishpg(index)                ! size at decomp. index of global field

    CALL get_dc_index(NG, idx, ldo)
    ishpl(:) = ishpg(:)              ! size of local field is size of global
                                     ! field, except at decomp. index
    ishpl(index) = idx(p_pe,2)-idx(p_pe,1) + 1

    DO i=1, nrank
       ixl(i,1) = 1                  ! start/stop index in local field ...
       ixl(i,2) = ishpl(i)           ! ... is given by shape
       ixg(i,1) = 1                  ! start/stop index in global field
       ixg(i,2) = ishpg(i)           ! is size of global field ...
                                     ! ... except at decomp. index
                                     ! -> set below
    END DO

    IF (ldo(p_pe)) THEN
       IF (.NOT.ASSOCIATED(lc)) THEN
          ALLOCATE(lc(ishpl(1),ishpl(2),ishpl(3),ishpl(4)))
       ELSE
          IF (SIZE(lc) /= PRODUCT(ishpl)) THEN
             WRITE(*,*) 'p_pe=',p_pe,': ',ishpl,' =/= ',UBOUND(lc)
             CALL error_bi('wrong size of local array', substr)
          END IF
       END IF
    END IF

    !
    ! send if p_pe = zxpe
    !
    IF (p_pe == zxpe) THEN
       DO i = 1, p_nprocs
          pe = i-1
          ! set global start/stop indices at decomp. index
          ixg(index, 1) = idx(pe,1)
          ixg(index, 2) = idx(pe,2)
          IF (pe /= p_pe) THEN
             IF (ldo(pe)) &
                  CALL p_send(gl(ixg(1,1):ixg(1,2)  &
                  ,   ixg(2,1):ixg(2,2)  &
                  ,   ixg(3,1):ixg(3,2)  &
                  ,   ixg(4,1):ixg(4,2)) &
                  , pe, tag_scatter_ix)
          ELSE
             IF (ldo(pe)) &
                  lc(ixl(1,1):ixl(1,2)  &
                  ,   ixl(2,1):ixl(2,2)  &
                  ,   ixl(3,1):ixl(3,2)  &
                  ,   ixl(4,1):ixl(4,2)) = &
                  gl(ixg(1,1):ixg(1,2)   &
                  ,   ixg(2,1):ixg(2,2)  &
                  ,   ixg(3,1):ixg(3,2)  &
                  ,   ixg(4,1):ixg(4,2))
          END IF
       END DO
    ELSE
       ! reveive if p_pe /= zxpe
       IF (ldo(p_pe)) CALL p_recv(lc(ixl(1,1):ixl(1,2)  &
            ,   ixl(2,1):ixl(2,2)  &
            ,   ixl(3,1):ixl(3,2)  &
            ,   ixl(4,1):ixl(4,2)) &
            , zxpe, tag_scatter_ix)
    END IF

    ! CLEAN UP
    DEALLOCATE(idx) ; NULLIFY(idx)
    DEALLOCATE(ldo) ; NULLIFY(ldo)

  END SUBROUTINE scatter_trix
  ! =========================================================================

  !-------------------------------------------------------------------------
  SUBROUTINE gather_trix_1d(gl, lc, xpe)


    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    REAL(dp), DIMENSION(:), POINTER  :: gl
    REAL(dp), DIMENSION(:), POINTER  :: lc
    INTEGER, INTENT(IN), OPTIONAL    :: xpe

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr='gather_trix'
    INTEGER, DIMENSION(:,:), POINTER :: idx=>NULL()
    LOGICAL, DIMENSION(:),   POINTER :: ldo=>NULL()
    INTEGER                          :: NG
    INTEGER                          :: NL, NLs
    INTEGER                          :: pe, zxpe
    INTEGER                          :: i

    IF (PRESENT(xpe)) THEN
       zxpe = xpe
    ELSE
       zxpe = p_io0
    END IF

    NL = SIZE(lc)      ! size of local field

    IF (p_pe == zxpe) THEN
       NG = NL         ! size of global field
    END IF

    ALLOCATE(idx(0:p_nprocs-1,2))
    !Start and end index of all p_pe's
    ALLOCATE(ldo(0:p_nprocs-1))


    ! SUM OVER ALL NON-zxpe p_pe
    IF (p_pe /= zxpe) THEN
       CALL p_send (NL, zxpe, tag_gather_ix)
    ELSE ! p_pe = zxpe

       DO i=1, p_nprocs
          pe = i-1
          IF (pe /= p_pe) THEN
             CALL p_recv( NLs, pe, tag_gather_ix)
             NG = NG + NLs
          ELSE
             NLs=NL
          END IF

          IF (i==1) THEN
             idx(pe,1)=0
          ELSE
             idx(pe,1)=idx(pe-1,1)
          END IF
          idx(pe,2)=idx(pe,1)+NLs-1
          IF (idx(pe,1)>idx(pe,2)) THEN
                ldo(pe)=.false.
          ELSE
                ldo(pe)=.true.
          END IF
       END DO
    END IF

    IF (p_pe == zxpe) THEN
       IF (ASSOCIATED(gl)) THEN
          IF (SIZE(gl) /= NG) THEN
             CALL error_bi('wrong size of global array', substr)
          END IF
       ELSE
          ALLOCATE(gl(NG))
       END IF
    ELSE
       IF (ASSOCIATED(gl)) DEALLOCATE(gl)
       NULLIFY(gl)
    ENDIF

    ! send if pe /= zxpe
    !
    IF (p_pe /= zxpe) THEN
       IF (ldo(p_pe)) CALL p_send (lc(:), zxpe, tag_gather_ix)
    ELSE ! p_pe = zxpe
       DO i=1, p_nprocs
          pe = i-1
!          pe = dcg(i)%pe
          IF (pe /= p_pe) THEN
             IF (ldo(pe)) &
                  CALL p_recv( gl(idx(pe,1):idx(pe,2)), pe, tag_gather_ix)
          ELSE
             IF (ldo(pe)) &
                  gl(idx(pe,1):idx(pe,2)) = lc(:)
          END IF
       END DO
    END IF

    ! CLEAN UP
    DEALLOCATE(idx) ; NULLIFY(idx)
    DEALLOCATE(ldo) ; NULLIFY(ldo)

  END SUBROUTINE gather_trix_1d
  !-------------------------------------------------------------------------

#endif
! ##########################################################################
! END #if defined(ICON)
! ##########################################################################

! ##########################################################################
! #if defined(MESSYDWARF)
! ##########################################################################
#if defined(MESSYDWARF)

  USE messy_main_constants_mem, ONLY: dp

  INTERFACE locate_in_decomp
     MODULE PROCEDURE locate_in_decomp_nrst
     MODULE PROCEDURE locate_in_decomp_4
  END INTERFACE
  PUBLIC :: locate_in_decomp
  PUBLIC :: get_dc_index

CONTAINS

  !--------------------------------------------------------------------------
  SUBROUTINE locate_in_decomp_nrst(status, lon, lat, pe, jp, jrow, jgx, jgy &
       , lprinterr)

    USE messy_main_mpi_bi,       ONLY: p_pe!, p_nprocs, isubpos, nboundlines &
                                    ! , nprocy
    USE messy_main_grid_def_mem_bi, ONLY: longitude_i, latitude_i    &
                                     , nlon, nlat
    USE messy_main_constants_mem, ONLY: iouerr

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    REAL(DP), INTENT(IN)  :: lon     ! LONGITUDE
    REAL(DP), INTENT(IN)  :: lat     ! LATITUDE
    INTEGER,  INTENT(OUT) :: pe      ! PROCESS ID
    INTEGER,  INTENT(OUT) :: jp      ! COLUMN
    INTEGER,  INTENT(OUT) :: jrow    ! ROW
    INTEGER,  INTENT(OUT), OPTIONAL :: jgx  ! lon index in global field
    INTEGER,  INTENT(OUT), OPTIONAL :: jgy  ! lat index in global field
    LOGICAL,  INTENT(IN),  OPTIONAL :: lprinterr ! .TRUE. for error output

    ! LOCAL
    INTEGER                               :: i,j,ix,jy,p

    ! INITIALISE
    status = 0
    IF (Present(jgx))  jgx = -1
    IF (Present(jgy))  jgy = -1
    PE   = -1
    jp   = -1
    jrow = -1

    ! CHECKS
    IF ((lon < -180.0_DP) .OR. (lon > 360.0_DP)) THEN
       status = 1  ! LONGITUDE OUT OF RANGE
       RETURN
    END IF
    !
    IF ((lat < -90.0_DP) .OR. (lat > 90.0_DP)) THEN
       status = 2  ! LATITUDE OUT OF RANGE
       RETURN
    END IF

    IF (lon <= longitude_i(1) .OR. lon > longitude_i(nlon+1) .OR. &
          lat < latitude_i(1) .OR. lat > latitude_i(nlat+1)) THEN
       ! geographical point not located in model domain
       status = 555
       RETURN
    ENDIF

    ! find global longitude
    ix = -1
    DO i=1 ,nlon
       IF (longitude_i(i) <= lon .AND. longitude_i(i+1) >= lon)  THEN
          ix = i
          EXIT
       END IF
    END DO

    ! find global latitude
    jy = -1
    do j=1, nlat
       IF (latitude_i(j) <= lat .AND. latitude_i(j+1) >= lat)  THEN
          jy = j
          EXIT
       END IF
    END DO

    IF (PRESENT(jgx)) jgx = ix
    IF (PRESENT(jgy)) jgy = jy

!!$    ix_jy_undef: IF (jy /= -1 .AND. ix /= -1) THEN
!!$
!!$       ! find PE and local indices
!!$       DO p = p_nprocs-1, 0, -1
!!$          ! NOTE: isubpos must not be <=1 but always >1 (due to halo)
!!$          IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
!!$             pe   = p
!!$             jp   = ix - isubpos(p,1) + nboundlines + 1
!!$             jrow = jy - isubpos(p,2) + nboundlines + 1
!!$             EXIT
!!$          END IF
!!$       END DO
!!$
!!$       ! geographical point not located in inner model domain
!!$       IF ( (pe == -1) .OR. (jp == -1) .or. (jrow == -1) )   THEN
!!$          if_pres: IF (PRESENT(linclbnd)) THEN
!!$             if_lincb: IF (linclbnd) THEN
!!$                ! 0 < ix <  isubpos(p,1)
!!$                IF (ix <= isubpos(0,1)) THEN
!!$                   ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
!!$                   IF (jy <= isubpos(0,2)) THEN
!!$                      pe   = 0
!!$                      jp   = ix-isubpos(0,1) + nboundlines + 1
!!$                      jrow = jy-isubpos(0,2) + nboundlines + 1
!!$                   ELSE
!!$                      ! .AND. jy > isubpos(p,2)
!!$                      ! => search Pe on left domain bound.
!!$                      DO p=nprocy-1,0,-1
!!$                         IF (jy >= isubpos(p,2)) THEN
!!$                            pe   = p
!!$                            jp   = ix-isubpos(p,1) + nboundlines + 1
!!$                            jrow = jy-isubpos(p,2) + nboundlines + 1
!!$                            EXIT
!!$                         END IF
!!$                      END DO
!!$                   END IF
!!$                ELSE IF (jy <= isubpos(0,2)) THEN
!!$                   ! case ix <= isubpos(p,1) already done
!!$                   DO p=p_nprocs-nprocy,0,-4
!!$                      IF (ix >= isubpos(p,1)) THEN
!!$                         pe   = p
!!$                         jp   = ix-isubpos(p,1) + nboundlines + 1
!!$                         jrow = jy-isubpos(p,2) + nboundlines + 1
!!$                         EXIT
!!$                      END IF
!!$                   END DO
!!$                END IF
!!$             END IF if_lincb
!!$          END IF if_pres
!!$
!!$       ENDIF
!!$
!!$    END IF ix_jy_undef
    IF (jy > 0 .AND. ix > 0) THEN
       pe = 1
       jp = ix
       jrow = jy
    END IF

    ! geographical point not located in model domain
    IF ( (pe == -1) .OR. (jp == -1) .or. (jrow == -1) ) THEN
       IF (PRESENT(lprinterr)) THEN
          IF (lprinterr) &
               write(iouerr,*) 'locate_in_decomp ', lon,lat &
               , ' IX ',ix,jy, pe, jp, jrow &
               , 'BOUND', longitude_i(1), longitude_i(nlon+1) &
                        ,  latitude_i(1),  latitude_i(nlat+1)
       ENDIF
       status = 555
       RETURN
    END IF

  END SUBROUTINE locate_in_decomp_nrst
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE locate_in_decomp_4(status, lon, lat, pe, jp, jrow, w, gjx, gjy &
       , lprinterr)

    ! LOCATES GEOGRAPHICAL POSITION IN VECTOR-PARALLEL-DECOMPOSITION
    !
    ! Author: Patrick Joeckel, DLR, May 2011

    USE messy_main_mpi_bi,       ONLY: p_pe!, p_nprocs, isubpos, nboundlines &
                                     !, nprocy
    USE messy_maingrid_def_mem_bi, ONLY: longitude_i, latitude_i, nlon, nlat &
                                     , longitude  , latitude
    USE messy_main_tools,        ONLY: nn_index, bilin_weight

    IMPLICIT NONE
    INTRINSIC :: PRESENT, REAL, MAX, MIN

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    REAL(DP),               INTENT(IN)  :: lon     ! LONGITUDE [-180 ... 180]
    REAL(DP),               INTENT(IN)  :: lat     ! LATITUDE  [-90 ... 90]
    INTEGER,  DIMENSION(4), INTENT(OUT) :: pe      ! PROCESS ID
    INTEGER,  DIMENSION(4), INTENT(OUT) :: jp      ! COLUMN
    INTEGER,  DIMENSION(4), INTENT(OUT) :: jrow    ! ROW
    REAL(DP), DIMENSION(4), INTENT(OUT) :: w       ! WEIGHT
    ! lon index in global field
    INTEGER,  DIMENSION(4), INTENT(OUT), OPTIONAL :: gjx
    ! lat index in global field
    INTEGER,  DIMENSION(4), INTENT(OUT), OPTIONAL :: gjy
    ! .TRUE. for error output
    LOGICAL,                INTENT(IN),  OPTIONAL :: lprinterr

    ! LOCAL
    LOGICAL,  SAVE                        :: linit = .FALSE.
    REAL(DP), DIMENSION(:), POINTER, SAVE :: philon => NULL()
    REAL(DP), DIMENSION(:), POINTER, SAVE :: philat => NULL()
    !
    REAL(DP)                                :: dlon, dlat
    !
    INTEGER                                 :: i, j, p
    INTEGER                                 :: jgx1, jgy1, jgx2, jgy2
    INTEGER                                 :: ix, jy
    REAL(DP), DIMENSION(4,2)                :: vn
    REAL(DP), DIMENSION(2)                  :: v

    ! INIT
    status = 0
    IF (Present(gjx))  gjx(:) = -1
    IF (Present(gjy))  gjy(:) = -1
    PE(:)   = -1
    jp(:)   = -1
    jrow(:) = -1
    w(:) = 0.0_dp

    ! CHECKS
    IF ((lon < -180.0_DP) .OR. (lon > 360.0_DP)) THEN
       status = 1  ! LONGITUDE OUT OF RANGE
       RETURN
    END IF
    !
    IF ((lat < -90.0_DP) .OR. (lat > 90.0_DP)) THEN
       status = 2  ! LATITUDE OUT OF RANGE
       RETURN
    END IF

    ! check, if point is in range
    IF (lon < longitude_i(1) .OR. lon > longitude_i(nlon+1) .OR. &
         lat < latitude_i(1) .OR. lat > latitude_i(nlat+1)) THEN
       ! geographical point not located in model domain
       status = 555
       RETURN
    ENDIF

    ! SET INDICES
    CALL nn_index(longitude, lon, jgx1, jgx2, 2.0_dp, status)
    CALL nn_index(latitude, lat, jgy1, jgy2, 2.0_dp, status)
    ! OUT OF RANGE
    IF (status /= 0) THEN
       status = 555
       RETURN
    END IF

    v(1) = lon
    v(2) = lat

    ix = MIN(jgx1, jgx2)
    jy = MIN(jgy1, jgy2)
    vn(1,1) = longitude(ix)
    vn(1,2) = latitude(jy)
!!$    DO p=p_nprocs-1,0,-1
!!$       ! NOTE: isubpos must not be <=1, always >1 (due to halo)
!!$       IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
!!$          pe(1)   = p
!!$          jp(1)   = ix-isubpos(p,1) + nboundlines +1
!!$          jrow(1) = jy-isubpos(p,2) + nboundlines +1
!!$          EXIT
!!$       END IF
!!$    END DO
    pe(1)   = p_pe
    jp(1)   = ix
    jrow(1) = jy

    IF (PRESENT(gjx)) gjx(1) = ix
    IF (PRESENT(gjy)) gjy(1) = jy

!!$    ! geographical point not located in inner model domain
!!$    IF ( (pe(1)== -1) .OR. (jp(1) == -1) .or. (jrow(1) == -1) )   THEN
!!$       if_pres: IF (PRESENT(linclbnd)) THEN
!!$          if_lincb: IF (linclbnd) THEN
!!$             ! 0 < ix <  isubpos(p,1)
!!$             IF (ix <= isubpos(0,1)) THEN
!!$                ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
!!$                IF (jy <= isubpos(0,2)) THEN
!!$                   pe(1)   = 0
!!$                   jp(1)   = ix-isubpos(0,1) + nboundlines + 1
!!$                   jrow(1) = jy-isubpos(0,2) + nboundlines + 1
!!$                ELSE
!!$                   ! .AND. jy > isubpos(p,2)
!!$                   ! => search Pe on left domain bound.
!!$                   DO p=nprocy-1,0,-1
!!$                      IF (jy >= isubpos(p,2)) THEN
!!$                         pe(1)   = p
!!$                         jp(1)   = ix-isubpos(p,1) + nboundlines + 1
!!$                         jrow(1) = jy-isubpos(p,2) + nboundlines + 1
!!$                         EXIT
!!$                      END IF
!!$                   END DO
!!$                END IF
!!$             ELSE IF (jy <= isubpos(0,2)) THEN
!!$                ! case ix <= isubpos(p,1) already done
!!$                DO p=p_nprocs-nprocy,0,-4
!!$                   IF (ix >= isubpos(p,1)) THEN
!!$                      pe(1)   = p
!!$                      jp(1)   = ix-isubpos(p,1) + nboundlines + 1
!!$                      jrow(1) = jy-isubpos(p,2) + nboundlines + 1
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$             END IF
!!$          END IF if_lincb
!!$       END IF if_pres
!!$
!!$    ENDIF

    ix = MAX(jgx1, jgx2)
    jy = MIN(jgy1, jgy2)
    vn(2,1) = longitude(ix)
    vn(2,2) = latitude(jy)
!!$    DO p=p_nprocs-1,0,-1
!!$       ! NOTE: isubpos must not be <=1, always >1 (due to halo)
!!$       IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
!!$          pe(2)   = p
!!$          jp(2)   = ix-isubpos(p,1) + nboundlines +1
!!$          jrow(2) = jy-isubpos(p,2) + nboundlines +1
!!$          EXIT
!!$       END IF
!!$    END DO
    pe(2)   = p_pe
    jp(2)   = ix
    jrow(2) = jy

    IF (PRESENT(gjx)) gjx(2) = ix
    IF (PRESENT(gjy)) gjy(2) = jy

!!$    ! geographical point not located in inner model domain
!!$    IF ( (pe(2)== -1) .OR. (jp(2) == -1) .or. (jrow(2) == -1) )   THEN
!!$       if_pres2: IF (PRESENT(linclbnd)) THEN
!!$          if_lincb2: IF (linclbnd) THEN
!!$             ! 0 < ix <  isubpos(p,1)
!!$             IF (ix <= isubpos(0,1)) THEN
!!$                ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
!!$                IF (jy <= isubpos(0,2)) THEN
!!$                   pe(2)   = 0
!!$                   jp(2)   = ix-isubpos(0,1) + nboundlines + 1
!!$                   jrow(2) = jy-isubpos(0,2) + nboundlines + 1
!!$                ELSE
!!$                   ! .AND. jy > isubpos(p,2)
!!$                   ! => search Pe on left domain bound.
!!$                   DO p=nprocy-1,0,-1
!!$                      IF (jy >= isubpos(p,2)) THEN
!!$                         pe(2)   = p
!!$                         jp(2)   = ix-isubpos(p,1) + nboundlines + 1
!!$                         jrow(2) = jy-isubpos(p,2) + nboundlines + 1
!!$                         EXIT
!!$                      END IF
!!$                   END DO
!!$                END IF
!!$             ELSE IF (jy <= isubpos(0,2)) THEN
!!$                ! case ix <= isubpos(p,1) already done
!!$                DO p=p_nprocs-nprocy,0,-4
!!$                   IF (ix >= isubpos(p,1)) THEN
!!$                      pe(2)   = p
!!$                      jp(2)   = ix-isubpos(p,1) + nboundlines + 1
!!$                      jrow(2) = jy-isubpos(p,2) + nboundlines + 1
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$             END IF
!!$          END IF if_lincb2
!!$       END IF if_pres2
!!$
!!$    ENDIF

    ix = MAX(jgx1, jgx2)
    jy = MAX(jgy1, jgy2)
    vn(3,1) = longitude(ix)
    vn(3,2) = latitude(jy)
!!$    DO p=p_nprocs-1,0,-1
!!$       ! NOTE: isubpos must not be <=1, always >1 (due to halo)
!!$       IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
!!$          pe(3)   = p
!!$          jp(3)   = ix-isubpos(p,1) + nboundlines +1
!!$          jrow(3) = jy-isubpos(p,2) + nboundlines +1
!!$          EXIT
!!$       END IF
!!$    END DO
    pe(3)   = p_pe
    jp(3)   = ix
    jrow(3) = jy

    IF (PRESENT(gjx)) gjx(3) = ix
    IF (PRESENT(gjy)) gjy(3) = jy

!!$    ! geographical point not located in inner model domain
!!$    IF ( (pe(3)== -1) .OR. (jp(3) == -1) .or. (jrow(3) == -1) )   THEN
!!$       if_pres3: IF (PRESENT(linclbnd)) THEN
!!$          if_lincb3: IF (linclbnd) THEN
!!$             ! 0 < ix <  isubpos(p,1)
!!$             IF (ix <= isubpos(0,1)) THEN
!!$                ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
!!$                IF (jy <= isubpos(0,2)) THEN
!!$                   pe(3)   = 0
!!$                   jp(3)   = ix-isubpos(0,1) + nboundlines + 1
!!$                   jrow(3) = jy-isubpos(0,2) + nboundlines + 1
!!$                ELSE
!!$                   ! .AND. jy > isubpos(p,2)
!!$                   ! => search Pe on left domain bound.
!!$                   DO p=nprocy-1,0,-1
!!$                      IF (jy >= isubpos(p,2)) THEN
!!$                         pe(3)   = p
!!$                         jp(3)   = ix-isubpos(p,1) + nboundlines + 1
!!$                         jrow(3) = jy-isubpos(p,2) + nboundlines + 1
!!$                         EXIT
!!$                      END IF
!!$                   END DO
!!$                END IF
!!$             ELSE IF (jy <= isubpos(0,2)) THEN
!!$                ! case ix <= isubpos(p,1) already done
!!$                DO p=p_nprocs-nprocy,0,-4
!!$                   IF (ix >= isubpos(p,1)) THEN
!!$                      pe(3)   = p
!!$                      jp(3)   = ix-isubpos(p,1) + nboundlines + 1
!!$                      jrow(3) = jy-isubpos(p,2) + nboundlines + 1
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$             END IF
!!$          END IF if_lincb3
!!$       END IF if_pres3
!!$
!!$    ENDIF

    ix = MIN(jgx1, jgx2)
    jy = MAX(jgy1, jgy2)
    vn(4,1) = longitude(ix)
    vn(4,2) = latitude(jy)
!!$    DO p=p_nprocs-1,0,-1
!!$       ! NOTE: isubpos must not be <=1, always >1 (due to halo)
!!$       IF (ix >= isubpos(p,1) .AND. jy >= isubpos(p,2) )   THEN
!!$          pe(4)   = p
!!$          jp(4)   = ix-isubpos(p,1) + nboundlines +1
!!$          jrow(4) = jy-isubpos(p,2) + nboundlines +1
!!$          EXIT
!!$       END IF
!!$    END DO
    pe(4)   = p_pe
    jp(4)   = ix
    jrow(4) = jy

    IF (PRESENT(gjx)) gjx(4) = ix
    IF (PRESENT(gjy)) gjy(4) = jy

!!$    ! geographical point not located in inner model domain
!!$    IF ( (pe(4)== -1) .OR. (jp(4) == -1) .or. (jrow(4) == -1) )   THEN
!!$       if_pres4: IF (PRESENT(linclbnd)) THEN
!!$          if_lincb4: IF (linclbnd) THEN
!!$             ! 0 < ix <  isubpos(p,1)
!!$             IF (ix <= isubpos(0,1)) THEN
!!$                ! .AND. 0 < jy < isubpos(p,2) => left lower corner (pe == 0)
!!$                IF (jy <= isubpos(0,2)) THEN
!!$                   pe(4)   = 0
!!$                   jp(4)   = ix-isubpos(0,1) + nboundlines + 1
!!$                   jrow(4) = jy-isubpos(0,2) + nboundlines + 1
!!$                ELSE
!!$                   ! .AND. jy > isubpos(p,2)
!!$                   ! => search Pe on left domain bound.
!!$                   DO p=nprocy-1,0,-1
!!$                      IF (jy >= isubpos(p,2)) THEN
!!$                         pe(4)   = p
!!$                         jp(4)   = ix-isubpos(p,1) + nboundlines + 1
!!$                         jrow(4) = jy-isubpos(p,2) + nboundlines + 1
!!$                         EXIT
!!$                      END IF
!!$                   END DO
!!$                END IF
!!$             ELSE IF (jy <= isubpos(0,2)) THEN
!!$                ! case ix <= isubpos(p,1) already done
!!$                DO p=p_nprocs-nprocy,0,-4
!!$                   IF (ix >= isubpos(p,1)) THEN
!!$                      pe(4)   = p
!!$                      jp(4)   = ix-isubpos(p,1) + nboundlines + 1
!!$                      jrow(4) = jy-isubpos(p,2) + nboundlines + 1
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$             END IF
!!$          END IF if_lincb4
!!$       END IF if_pres4
!!$
!!$    ENDIF

    CALL bilin_weight(vn, v, w)

    status = 0

  END SUBROUTINE locate_in_decomp_4
  !-------------------------------------------------------------------------

#endif

! ##########################################################################
! END #if defined(MESSYDWARF)
! ##########################################################################

#ifdef EXPLICIT
    SUBROUTINE pack_fs_buf_ex (fs, buf, intr, nf,nb,n2)
      !
      ! explicit shape array argument version of pack_fs_buf
      ! (enforces vectorization)
      !
      REAL(dp):: fs  (nf,n2)
      REAL(dp):: buf (nb,n2)
      INTEGER :: intr(nb)
      INTEGER :: nf,nb,n2
      INTEGER :: i
      DO i=1,nb
         !CDIR SELECT(VECTOR)
         buf(i,:) = fs (intr(i),:)
      END DO
    END SUBROUTINE pack_fs_buf_ex
    !----------------------------------------------------------------

    !-----------------------------------------------------------------
    SUBROUTINE unpack_buf_fs_ex (buf, fs, intr, nb,nf,n2)
      !
      ! explicit shape size array argument version of unpack_buf_fs
      ! (enforces vectorization)
      !
      REAL(dp):: buf (nb,n2)
      REAL(dp):: fs  (nf,n2)
      INTEGER :: intr(nb)
      INTEGER :: nf,nb,n2
      INTEGER :: i
      DO i=1,nb
         !CDIR SELECT(VECTOR)
         fs (intr(i),:) = buf(i,:)
      END DO
    END SUBROUTINE unpack_buf_fs_ex
    !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls_ex(ls,buf,n1,nfls,nf,nv,i1,i2,j1,j2)
      REAL(dp):: ls(n1,nfls,nv)
      REAL(dp):: buf (n1,nf,nv)
      INTEGER :: i1,i2,j1,j2
      INTEGER :: i,j,k

      DO k=1,nv
         DO j=i1,i2
            DO i=1,n1
               ls(i,j,k)=buf(i,j-i1+j1,k)
            END DO
         END DO
      END DO
    END SUBROUTINE unpack_buf_ls_ex
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf_ex(buf,ls,n1,nfls,nf,nv,i1,i2,j1,j2)
      REAL(dp):: ls(n1,nfls,nv)
      REAL(dp):: buf (n1,nf,nv)
      INTEGER :: i1,i2,j1,j2
      INTEGER :: i,j,k

      DO k=1,nv
         DO j=i1,i2
            DO i=1,n1
               buf(i,j-i1+j1,k)=ls(i,j,k)
            END DO
         END DO
      END DO
    END SUBROUTINE pack_ls_buf_ex
    !----------------------------------------------------------------
#endif

  ! =========================================================================
  SUBROUTINE get_dc_index(N, IDX, LDO)

    USE messy_main_mpi_bi,         ONLY: p_nprocs

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT

    ! I/O
    INTEGER,                                   INTENT(IN) :: N
    INTEGER, DIMENSION(:,:), POINTER                      :: IDX
    LOGICAL, DIMENSION(:),   POINTER, OPTIONAL            :: LDO

    ! LOCAL
    INTEGER :: DN, i

    IF (ASSOCIATED(IDX)) THEN
       DEALLOCATE(IDX)
       NULLIFY(IDX)
    END IF
    ALLOCATE(IDX(0:p_nprocs-1,2))
    ! START INDEX: IDX(p_pe,1)
    ! STOP  INDEX: IDX(p_pe,2)
    IF (PRESENT(LDO)) THEN
       IF (ASSOCIATED(LDO)) THEN
          DEALLOCATE(LDO)
          NULLIFY(LDO)
       END IF
       ALLOCATE(LDO(0:p_nprocs-1))
    END IF

    DN = N/p_nprocs
    DO i=0, p_nprocs-1
       IDX(i,1) = i*DN + 1
       IDX(i,2) = IDX(i,1) + DN - 1
    END DO
    ! CORRECT FOR ODD Ns
    IDX(p_nprocs-1,2) = N

    IF (PRESENT(LDO)) THEN
       DO i=0, p_nprocs-1
          IF (IDX(i,1) > IDX(i,2)) THEN
             LDO(i) = .false.
          ELSE
             LDO(i) = .true.
          END IF
       END DO
    END IF

  END SUBROUTINE get_dc_index
  ! =========================================================================
! **************************************************************************
END MODULE messy_main_transform_bi
! **************************************************************************
