!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide YAXT support.
!!
!! @remarks
!!   This module contains routines for description of domain
!!   decomposition for parallel output with CDI-PIO and global data
!!   transpositions using YAXT.
!!
!! @author Jörg Behrens, DKRZ, Hamburg (2013-06-20): Initial version
!!         Irina Fast,   DKRZ, Hamburg (2013-07-18): Revised version
!!
!! (jb, 2015): modified version - adapted to messy_d2.52c/echam5.3.02

!#ifndef DEBUG
!#define DEBUG
!#endif

!#define HAVE_YAXT_0_4_5

MODULE mo_echam_yaxt
  USE mo_kind,          ONLY: dp
  USE mo_control,       ONLY: lyaxt_transposition, ltimer
  USE mo_exception,     ONLY: message, message_text, finish
#ifdef HAVE_YAXT
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_loc, c_null_ptr, c_ptr, c_int
  USE mo_mpi,           ONLY: p_real_dp, p_pe, p_io
  USE mo_decomposition, ONLY: ldc => local_decomposition, debug_parallel, &
       & gdc => global_decomposition
  USE mo_control,       ONLY: control_nsp => nsp
#ifndef MESSY
    USE mo_tracer,        ONLY: advection_jps => jps, trlist, ntrac
#else
    USE mo_parameters,    ONLY: advection_jps => jps
    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp, on_switch => ON, I_ADVECT
#endif
  USE yaxt,             ONLY: xt_initialize, xt_finalize, xt_abort,       &
                              xt_idxlist, xt_idxvec_new, xt_idxempty_new, &
                              xt_xmap, xt_xmap_all2all_new, xt_redist, &
                              xt_redist_p2p_off_new, xt_redist_s_exchange1, &
                              xt_redist_collection_new, xt_redist_s_exchange, &
                              xt_stripe, xt_idxstripes_new, xt_is_null, &
                              xt_idxlist_c2f, xt_xmap_delete, xt_redist_repeat_new
  USE mo_real_timer, ONLY: new_timer, timer_start, timer_stop
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: yaxt_initialize, yaxt_finalize, setup_yaxt_decomposition, &
            generate_yaxt_redist, add_yaxt_gp_nlevs, yaxt_init_gp_coords

  PUBLIC :: yaxt_tr_gp_to_ffsl, yaxt_tr_ffsl_to_gp, &
       &    yaxt_tr_gp_to_fs, yaxt_tr_fs_to_gp, &
       &    yaxt_tr_fs_to_ls, yaxt_tr_ls_to_fs, &
       &    yaxt_tr_ls_to_sp, yaxt_tr_sp_to_ls
  PUBLIC :: yaxt_tr_gp_to_ffsl_2d, yaxt_tr_gp_to_ffsl_3d, yaxt_tr_ffsl_to_gp_3d, &
       &    yaxt_tr_gp_to_ffsl_3d_repeated, yaxt_ffsl_ghost_update_3d_repeated

  PUBLIC :: yaxt_simple_gather_gp_3d
  PUBLIC :: yaxt_simple_cross_gather_gp_3d
  PUBLIC :: yaxt_simple_gather_gp_2d
  PUBLIC :: yaxt_simple_cross_gather_gp_2d

  PUBLIC :: yaxt_simple_scatter_gp_2d
  PUBLIC :: yaxt_simple_scatter_gp_3d

  PUBLIC :: yaxt_ffsl_ghost_update_2d, yaxt_ffsl_ghost_update_3d, yaxt_ffsl_ghost_update_4d
  PUBLIC :: yaxt_max_tracer_aggregate_size

#ifdef HAVE_YAXT
  PUBLIC :: ilon_map, ilat_map, nonid_ilon_map, nonid_ilat_map

  PUBLIC :: gp_gdeco, sp_sdeco
#endif

! check memory layout expectations of ECHAM data
#define CHECK_MEM_STRIDES
#define CHECK_SHAPES

  ! methods for generation of index lists
  INTEGER, PARAMETER :: idx_descr_vector  = 1
  INTEGER, PARAMETER :: idx_descr_stripes = 2
  INTEGER, PARAMETER :: idx_descr_section = 3

#ifdef HAVE_YAXT

  ! agregation switches (aggregation leads to bigger and fewer messages):
  LOGICAL, PARAMETER :: aggregate_default    = .TRUE.
  LOGICAL, PARAMETER :: aggregate_gp_to_fs   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_fs_to_gp   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_fs_to_ls   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_ls_to_fs   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_ls_to_sp   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_sp_to_ls   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_gp_to_ffsl = aggregate_default
#ifdef HAVE_YAXT_0_4_5
  LOGICAL, PARAMETER :: aggregate_ffsl_to_gp = aggregate_default
#else
  LOGICAL, PARAMETER :: aggregate_ffsl_to_gp = .FALSE.
#endif

  ! tie each decomposition (idxlist) in echam to a fixed data layout (offsets):
  TYPE yaxt_deco
    TYPE(xt_idxlist)     :: idxlist
    INTEGER, ALLOCATABLE :: offset(:)
  END TYPE yaxt_deco

  ! Note:
  ! We distinguish 5 kinds of index spaces: gdeco, sdeco, zmdeco, m0deco and hdeco.
  ! The general decomposition naming is: {short decomposition name}_{index space name}
  ! Only decompositions belonging to the same index space are combined within a redist object.
  ! The decompositions are tied to offsets as used in ECHAM.

  ! Decos depending on various level counts:
  TYPE(yaxt_deco), ALLOCATABLE, SAVE :: gp_gdeco(:) ! standard gridpoint decomposition of grid point data
  TYPE(yaxt_deco), ALLOCATABLE, SAVE :: sp_sdeco(:) ! fourier-decomposition of grid point data
  TYPE(yaxt_deco), ALLOCATABLE, SAVE :: fs_gdeco(:) ! fourier-decomposition of grid point data
  TYPE(yaxt_deco), ALLOCATABLE, SAVE :: ls_sdeco(:) ! legendre decomposition of spectral point data

  ! 1-dim index space permutations due to nlonstrip and nlatstrip:
  ! nlonstrip=/0 or  nlatstrip/=0 is not supported yet:
  INTEGER,PARAMETER :: nlonstrip = 0
  INTEGER,PARAMETER :: nlatstrip = 0
  INTEGER, SAVE, ALLOCATABLE :: ilon_map(:), ilat_map(:)
  LOGICAL, PARAMETER :: preserve_latsym =.TRUE.
  LOGICAL, SAVE :: nonid_ilon_map = .FALSE.
  LOGICAL, SAVE :: nonid_ilat_map = .FALSE.

  ! global gp decos:
  TYPE(yaxt_deco), SAVE:: global_gp_2d_gdeco
  TYPE(yaxt_deco), SAVE:: global_gp_3d_gdeco

  ! empty deco:
  TYPE(yaxt_deco), SAVE :: empty_deco

  ! standard gridpoint decomposition for tracers:
  TYPE(yaxt_deco), SAVE:: gp_3d_fb_single_tracer_gdeco ! offsets as in mo_tpcore::fb_gp for first cnst
  TYPE(yaxt_deco), SAVE:: gp_3d_xtm1_single_tracer_gdeco ! offsets as in mo_tpcore::xtm1 for first tracer

  ! FFSL decompositions:
  TYPE(yaxt_deco), SAVE :: ffsl_2d_gdeco ! decomposition description for FFSL 2d
  TYPE(yaxt_deco), SAVE :: ffsl_3d_gdeco ! decomposition description for FFSL 3d

  INTEGER, PARAMETER :: max_ffsl_ghost_size = 3
  ! ffsl_2d_gz_gdeco(istate,ghost_size)
  !    with ghost_size = ghost data latitude size,
  !    and istate = {0:src state; 1:target state}
  TYPE(yaxt_deco), SAVE :: ffsl_2d_gz_gdeco(0:1,3) ! FFSL 2d with ghost zones
  TYPE(yaxt_deco), SAVE :: ffsl_3d_gz_gdeco(0:1,3) ! FFSL 3d with ghost zones

  ! subset of fs_gdeco(nlev+1), used for fs <--> ls:
  TYPE(yaxt_deco), SAVE :: fs_subset_gdeco

  ! 2d FS decomposition:
  TYPE(yaxt_deco), SAVE :: fs_2d_gdeco

  ! decomposition of GP zonal means for all levels:
  TYPE(yaxt_deco), SAVE :: gp_zmdeco

  ! decomposition of GP horizontal means for all levels:
  TYPE(yaxt_deco), SAVE :: gp_hmdeco

  ! decomposition of GP zonal means for selected levels (fs levels) only
  TYPE(yaxt_deco), SAVE :: gp_sel_zmdeco

  ! decomposition of FS zonal means:
  TYPE(yaxt_deco), SAVE :: fs_zmdeco

  ! decomposition of FS horizontal means:
  TYPE(yaxt_deco), SAVE :: fs_hmdeco

  ! 3d LS decomposition:
  TYPE(yaxt_deco), SAVE :: ls_3d_gdeco

  ! decomposition of LS zonal means:
  TYPE(yaxt_deco), SAVE :: ls_zmdeco

  ! decomposition of LS coeff. with m = 0:
  TYPE(yaxt_deco), SAVE :: ls_m0deco

  ! decomposition of SP coeff. with m = 0:
  TYPE(yaxt_deco), SAVE :: sp_m0deco

#endif

  ! global dimensions
  INTEGER :: nlon, nlat, nlev

  ! gp-deco info:
  INTEGER :: nproca, nprocb, nproma, npromz, ngpblks
  INTEGER :: glons(2), glone(2), glats(2), glate(2)

  ! ffsl-deco info:
  INTEGER :: ffsl_nlat, ffsl_nlev
  INTEGER :: ffsl_gp_lat1, ffsl_gp_lat2 ! ffsl lat-borders in gp-oriented coords

  ! fs-deco info:
  INTEGER :: flats(2), flate(2), nflat, flevs, fleve, nflev, nflevp1

  ! tracer info:
  INTEGER, PARAMETER :: yaxt_max_tracer_aggregate_size = 32
  INTEGER :: pcnst, trdim_size, active_trnum
  INTEGER, ALLOCATABLE :: active_trpos(:) ! position of active tracer within trlist

  ! mtracer_fb_block_X_size, with X=max,tail: these are tracer aggregation sizes used as
  ! temp. workaround until yaxt releases more flexible redist_repeat constructors.
  INTEGER, PARAMETER :: ffsl_4d_max_block_size = 4
  INTEGER :: ffsl_4d_block_size

#ifdef HAVE_YAXT
  ! single var redists:

  ! gp -> ffsl:
  TYPE(xt_redist), SAVE :: gp2ffsl_4d_xtm1_redist
  TYPE(xt_redist), SAVE :: gp2ffsl_3d_xtm1_single_tracer_redist
  TYPE(xt_redist), SAVE :: gp2ffsl_3d_redist
  TYPE(xt_redist), SAVE :: gp2ffsl_2d_redist
  TYPE(xt_redist), SAVE :: gp2ffsl_all_redist

  ! ffsl -> gp:
#ifdef HAVE_YAXT_0_4_5
  TYPE(xt_redist), SAVE :: ffsl2gp_4d_fb_redist
#endif
  TYPE(xt_redist), SAVE :: ffsl2gp_3d_fb_single_tracer_redist
  TYPE(xt_redist), SAVE :: ffsl2gp_3d_redist
  TYPE(xt_redist), SAVE :: ffsl2gp_2d_redist
  TYPE(xt_redist), SAVE :: ffsl2gp_all_redist

  ! ffsl ghost update:
  TYPE(xt_redist), SAVE :: ffsl_update_2d_redist(max_ffsl_ghost_size)
  TYPE(xt_redist), SAVE :: ffsl_update_3d_redist(max_ffsl_ghost_size)
  TYPE(xt_redist), SAVE :: ffsl_update_4d_redist(max_ffsl_ghost_size)

  ! gp -> fs:
  TYPE(xt_redist), SAVE :: gp2fs_3d_redist
  TYPE(xt_redist), SAVE :: gp2fs_2d_redist
  TYPE(xt_redist), SAVE :: gp2fs_zm_redist
  TYPE(xt_redist), SAVE :: gp2fs_all_redist

  ! fs -> gp:
  TYPE(xt_redist), SAVE :: fs2gp_3d_redist
  TYPE(xt_redist), SAVE :: fs2gp_2d_redist
  TYPE(xt_redist), SAVE :: fs2gp_zm_redist
  TYPE(xt_redist), SAVE :: fs2gp_hm_redist
  TYPE(xt_redist), SAVE :: fs2gp_all_redist

  ! fs -> ls:
  TYPE(xt_redist), SAVE :: fs2ls_3d_redist
  TYPE(xt_redist), SAVE :: fs2ls_zm_redist
  TYPE(xt_redist), SAVE :: fs2ls_all_redist

  ! ls -> fs:
  TYPE(xt_redist), SAVE :: ls2fs_3d_redist
  TYPE(xt_redist), SAVE :: ls2fs_zm_redist
  TYPE(xt_redist), SAVE :: ls2fs_all_redist

  ! ls -> sp:
  TYPE(xt_redist), ALLOCATABLE, SAVE :: ls2sp_3d_redist(:)
  TYPE(xt_redist), SAVE :: ls2sp_m0_redist
  TYPE(xt_redist), SAVE :: ls2sp_all_redist

  ! sp -> ls:
  TYPE(xt_redist), ALLOCATABLE, SAVE :: sp2ls_3d_redist(:)
  TYPE(xt_redist), SAVE :: sp2ls_m0_redist
  TYPE(xt_redist), SAVE :: sp2ls_all_redist

  ! gather gp within nproca x nprocb process space:
  TYPE(xt_redist), SAVE :: simple_gather_gp_2d_redist
  TYPE(xt_redist), SAVE :: simple_gather_gp_3d_redist

  ! scatter gp within nproca x nprocb process space:
  TYPE(xt_redist), SAVE :: simple_scatter_gp_2d_redist
  TYPE(xt_redist), SAVE :: simple_scatter_gp_3d_redist

  ! cross gather from compute instance to debug instance
  TYPE(xt_redist), SAVE :: simple_cross_gather_gp_2d_redist
  TYPE(xt_redist), SAVE :: simple_cross_gather_gp_3d_redist

  ! repeated  3d-tracer communication:
   TYPE(xt_redist), SAVE :: gp2ffsl_m3d_redist(yaxt_max_tracer_aggregate_size)
   TYPE(xt_redist), SAVE :: ffsl_update_m3d_redist(yaxt_max_tracer_aggregate_size,max_ffsl_ghost_size)
#ifdef CHECK_MEM_STRIDES
  INTERFACE check_mem_strides
    MODULE PROCEDURE check_mem_strides_4d
    MODULE PROCEDURE check_mem_strides_3d
    MODULE PROCEDURE check_mem_strides_2d
  END INTERFACE check_mem_strides
#endif

#endif /* HAVE_YAXT */

  ! nproca x nprocb domain:
  INTEGER :: ab_comm, ab_rank, ab_size

  LOGICAL, SAVE :: decos_are_finalized = .FALSE.
  LOGICAL, SAVE, ALLOCATABLE :: gp_nlevs(:)

  ! gather:
  LOGICAL, SAVE :: is_gather_target = .FALSE.
  LOGICAL, SAVE :: is_cross_gather_target = .FALSE.
  LOGICAL, SAVE :: require_cross_gather = .FALSE.

  ! timer
  INTEGER, SAVE :: it_gp2ffsl = 0
  INTEGER, SAVE :: it_ffsl2gp = 0

CONTAINS

  SUBROUTINE yaxt_initialize(comm)
    INTEGER, INTENT(IN) :: comm     ! mpi communicator
#ifndef HAVE_YAXT
  IF (lyaxt_transposition) THEN
    WRITE (message_text, *) 'For array transpositions with YAXT '// &
       &'specification of -DHAVE_YAXT is required. Please, recompile '// &
       &'and link with YAXT libary OR set lyaxt_transposition=.FALSE. '// &
       &'in namelist RUNCTL.'
    CALL finish('yaxt_initialize', message_text)
  ENDIF
#else
  CALL xt_initialize(comm)
  it_gp2ffsl = new_timer('yaxt_gp2ffsl')
  it_ffsl2gp = new_timer('yaxt_ffsl2gp')
#endif
  END SUBROUTINE yaxt_initialize

!-------------------------------------------------------------------------------

  SUBROUTINE yaxt_finalize()
#ifdef HAVE_YAXT
    CALL xt_finalize()
#endif
  END SUBROUTINE yaxt_finalize

!-------------------------------------------------------------------------------
  SUBROUTINE yaxt_init_gp_coords
#ifdef HAVE_YAXT
    ! Gaussian grid decomposition data
    nproca  = ldc%nproca
    nprocb  = ldc%nprocb
    nproma  = ldc%nproma
    npromz  = ldc%npromz
    ngpblks = ldc%ngpblks

    nlon  = ldc%nlon
    nlat  = ldc%nlat
    nlev  = ldc%nlev

    glons = ldc%glons
    glone = ldc%glone
    glats = ldc%glats
    glate = ldc%glate

    ! generate index space permutations:
    CALL gen_index_maps

    !IF (nonid_ilon_map .OR. nonid_ilat_map) CALL prepare_new_gp_coords()
    CALL prepare_new_gp_coords
#endif
  END SUBROUTINE yaxt_init_gp_coords

  SUBROUTINE add_yaxt_gp_nlevs(levels)
    INTEGER, INTENT(in) :: levels(:)
#ifdef HAVE_YAXT
    INTEGER :: min_lev, nlevs_size_new, nlevs_size_old
    LOGICAL, ALLOCATABLE :: tmp(:)

    ! For every 3d-array that is to be written via mo_output, this subroutine must
    ! have been called with the size of third dimension as element of the argument vector.
    ! Multiple calls are possible but must precede the call of setup_yaxt_decomposition.
    ! This routine has to be called at least once if setup_yaxt_decomposition is called

    IF (decos_are_finalized) THEN
      WRITE (message_text, *) 'Late call of add_yaxt_gp_nlevs. Yaxt decos are already finalized.'
      CALL finish('mo_echam_yaxt: add_yaxt_gp_nlevs', 'Decos are already finalized.')
    ENDIF

    min_lev = MINVAL(levels)
    IF (min_lev<1) THEN
      WRITE (message_text, *) 'invalid level',MINVAL(levels)
      CALL finish('mo_echam_yaxt: add_yaxt_gp_nlevs', message_text)
    ENDIF

    nlevs_size_old = SIZE(gp_nlevs)
    nlevs_size_new = MAXVAL(levels)

    IF (nlevs_size_new > nlevs_size_old) THEN
      IF (ALLOCATED(gp_nlevs)) THEN
        ! realloc:
        ALLOCATE(tmp(nlevs_size_old))
        tmp = gp_nlevs
        DEALLOCATE(gp_nlevs)
        ALLOCATE(gp_nlevs(nlevs_size_new))
        gp_nlevs(1:nlevs_size_old) = tmp
        gp_nlevs(nlevs_size_old+1:) = .FALSE.
      ELSE
        ! alloc:
        ALLOCATE(gp_nlevs(nlevs_size_new))
        gp_nlevs = .FALSE.
      ENDIF
    END IF

    gp_nlevs(levels(:)) = .TRUE.
#else

#endif
  END SUBROUTINE add_yaxt_gp_nlevs

  SUBROUTINE setup_yaxt_decomposition()
#ifdef HAVE_YAXT
    INTEGER :: itracer, nt

    ! FFSL decomposition data
    ! translate ffsl region boundaries to gp latitudes
    ffsl_gp_lat1 = nlat + 1 - ldc%ffsl%latn
    ffsl_gp_lat2 = nlat + 1 - ldc%ffsl%lats
    ffsl_nlat = ffsl_gp_lat2 - ffsl_gp_lat1 + 1

    ! FS decomposition data
    flats = ldc%flats
    flate = ldc%flate
    nflat = ldc%nflat
    flevs = ldc%flevs
    fleve = ldc%fleve
    nflev = ldc%nflev
    nflevp1 = ldc%nflevp1

    ! tracer data:
    trdim_size = get_trdim_size()
    active_trnum = get_active_trnum()

    pcnst = advection_jps + active_trnum
    ALLOCATE(active_trpos(active_trnum))
#ifndef MESSY
    nt = 0
    DO itracer = 1, tracdef_trlist%ntrac
      IF (tracdef_trlist%ti(itracer)%ntran /= 0) THEN
        nt = nt + 1
        active_trpos(nt) = itracer
      ENDIF
    ENDDO
#else
    nt = 0
    DO itracer=1, ntrac_gp
      IF (ti_gp(itracer)%tp%meta%cask_i(I_advect) == on_switch) THEN
        nt = nt + 1
        active_trpos(nt) = itracer
      END IF
    END DO
#endif

    ! gather targets:
    is_gather_target = (ldc%spe==1 .AND. p_pe==0) .OR. (ldc%spe==2 .AND. p_pe==1)
    require_cross_gather = (debug_parallel >= 0)
    is_cross_gather_target = require_cross_gather .AND. (ldc%spe==1 .AND. p_pe==0)

    ! empty deco:
    CALL setup_empty_deco

    ! define index lists and offsets for Gaussian grid
    CALL setup_gp_gdeco(idx_descr_stripes)
    !CALL setup_gp_gdeco(idx_descr_vector)

    ! spectral coeff. space decompositions:
    CALL setup_sp_sdeco(idx_descr_stripes)
    CALL setup_ls_sdeco

    CALL setup_fs_2d_gdeco
    CALL setup_fs_3d_gdeco
    CALL setup_ls_3d_gdeco

    ! define index lists and offsets for FFSL 2D data
    CALL setup_ffsl_2d_gdeco(idx_descr_stripes)
    ! define index lists and offsets for FFSL 3D data
    CALL setup_ffsl_3d_gdeco(idx_descr_stripes)

    ! zonal-means space decompositions:
    CALL setup_gp_zmdeco
    CALL setup_gp_hmdeco
    CALL setup_fs_zmdeco
    CALL setup_fs_hmdeco
    CALL setup_ls_zmdeco

    ! spectral m==0 coeff. space decompositions:
    CALL setup_sp_m0deco
    CALL setup_ls_m0deco

    decos_are_finalized = .TRUE.
#endif /* HAVE_YAXT */
  END SUBROUTINE setup_yaxt_decomposition

  SUBROUTINE generate_yaxt_redist(model_comm)
    INTEGER, INTENT(IN) :: model_comm  ! communicator of application without I/O servers
#ifdef HAVE_YAXT
    INTEGER, PARAMETER :: max_num_redists = 15
    TYPE(xt_redist) :: redists(max_num_redists)
    INTEGER, PARAMETER :: redist_cs = 1 ! best setting for static aggregation
    INTEGER :: num_redists, ir, n
    CALL set_comms(model_comm)

    ! generate gp <--> ffsl redistributions:

    gp2ffsl_2d_redist = new_redist(gp_gdeco(1), ffsl_2d_gdeco)
    ffsl2gp_2d_redist = new_redist(ffsl_2d_gdeco, gp_gdeco(1))
    gp2ffsl_3d_redist = new_redist(gp_gdeco(nlev), ffsl_3d_gdeco)
    ffsl2gp_3d_redist = new_redist(ffsl_3d_gdeco, gp_gdeco(nlev))
    gp2ffsl_3d_redist = new_redist(gp_gdeco(nlev), ffsl_3d_gdeco)
    ffsl2gp_3d_redist = new_redist(ffsl_3d_gdeco, gp_gdeco(nlev))
    ffsl2gp_3d_fb_single_tracer_redist = new_redist(ffsl_3d_gdeco, gp_3d_fb_single_tracer_gdeco)
    gp2ffsl_m3d_redist(1) = gp2ffsl_3d_redist

    IF (active_trnum>0) gp2ffsl_3d_xtm1_single_tracer_redist = &
         & new_redist(gp_3d_xtm1_single_tracer_gdeco, ffsl_3d_gdeco)

    IF (aggregate_gp_to_ffsl) THEN
      redists(1:5) = gp2ffsl_3d_redist
      redists(6:7) = gp2ffsl_2d_redist
      num_redists = 7
      gp2ffsl_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    IF (aggregate_ffsl_to_gp) THEN
#ifdef HAVE_YAXT_0_4_5
      CALL my_ffsl2gp_4d_redist
      redists(1) = ffsl2gp_4d_fb_redist
      redists(2) = ffsl2gp_2d_redist
      redists(3) = ffsl2gp_3d_redist
      num_redists = 3
      ffsl2gp_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
#else
      CALL die("ffsl2gp_4d_fb_redist not supported",__LINE__)
#endif
    ENDIF

    ! ffsl ghost update:
    ffsl_update_2d_redist(1) = new_redist(ffsl_2d_gz_gdeco(0,1), ffsl_2d_gz_gdeco(1,1))
    ffsl_update_2d_redist(3) = new_redist(ffsl_2d_gz_gdeco(0,3), ffsl_2d_gz_gdeco(1,3))
    ffsl_update_3d_redist(1) = new_redist(ffsl_3d_gz_gdeco(0,1), ffsl_3d_gz_gdeco(1,1))
    ffsl_update_3d_redist(3) = new_redist(ffsl_3d_gz_gdeco(0,3), ffsl_3d_gz_gdeco(1,3))
    !CALL my_big_ffsl_redists(1)
    CALL my_big_ffsl_update_redists(3)


    ! generate gp <--> fs redistributions:
    gp2fs_3d_redist = new_redist(gp_gdeco(nlev), fs_gdeco(nlev))
    fs2gp_3d_redist = new_redist(fs_gdeco(nlev), gp_gdeco(nlev))

    gp2fs_2d_redist = new_redist(gp_gdeco(1), fs_2d_gdeco)
    fs2gp_2d_redist = new_redist(fs_2d_gdeco, gp_gdeco(1))

    gp2fs_zm_redist = new_redist(gp_sel_zmdeco, fs_zmdeco)
    fs2gp_zm_redist = new_redist(fs_zmdeco, gp_zmdeco)

    ! horizontal average: only one direction
    fs2gp_hm_redist = new_redist(fs_hmdeco, gp_hmdeco)

    IF (aggregate_gp_to_fs) THEN
      DO ir = 1, 7
        redists(ir) = gp2fs_3d_redist
      ENDDO
      redists(8) = gp2fs_2d_redist
      DO ir = 9, 11
        redists(ir) = gp2fs_zm_redist
      ENDDO
      num_redists = 11
      gp2fs_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    IF (aggregate_fs_to_gp) THEN
      DO ir = 1, 9
        redists(ir) = fs2gp_3d_redist
      ENDDO
      DO ir = 10, 12
        redists(ir) = fs2gp_2d_redist
      ENDDO
      DO ir = 13, 15
        redists(ir) = fs2gp_zm_redist
      ENDDO
      num_redists = 15
      fs2gp_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    ! generate fs <--> ls redistributions:
    fs2ls_3d_redist = new_redist(fs_subset_gdeco, ls_3d_gdeco)
    ls2fs_3d_redist = new_redist(ls_3d_gdeco, fs_subset_gdeco)

    fs2ls_zm_redist = new_redist(fs_zmdeco, ls_zmdeco)
    ls2fs_zm_redist = new_redist(ls_zmdeco, fs_zmdeco)

    IF (aggregate_fs_to_ls) THEN
      DO ir = 1, 6
        redists(ir) = fs2ls_3d_redist
      ENDDO
      redists(7) = fs2ls_zm_redist
      num_redists = 7
      fs2ls_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    IF (aggregate_ls_to_fs) THEN
      DO ir = 1, 9
        redists(ir) = ls2fs_3d_redist
      ENDDO
      DO ir = 10, 12
        redists(ir) = ls2fs_zm_redist
      ENDDO
      num_redists = 12
      ls2fs_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    ! generate ls <--> sp redistributions:
    ALLOCATE(ls2sp_3d_redist(nlev:nlev+1), sp2ls_3d_redist(nlev:nlev+1))
    DO n = nlev, nlev+1
      ls2sp_3d_redist(n) = new_redist(ls_sdeco(n), sp_sdeco(n))
      sp2ls_3d_redist(n) = new_redist(sp_sdeco(n), ls_sdeco(n))
    ENDDO
    ls2sp_m0_redist = new_redist(ls_m0deco, sp_m0deco)
    sp2ls_m0_redist = new_redist(sp_m0deco, ls_m0deco)


    IF (aggregate_ls_to_sp) THEN
      redists(1) = ls2sp_3d_redist(nlev)
      redists(2) = ls2sp_3d_redist(nlev)
      redists(3) = ls2sp_3d_redist(nlev+1)
      redists(4) = ls2sp_m0_redist
      num_redists = 4
      ls2sp_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    IF (aggregate_sp_to_ls) THEN
      redists(1) = sp2ls_3d_redist(nlev)
      redists(2) = sp2ls_3d_redist(nlev)
      redists(3) = sp2ls_3d_redist(nlev+1)
      redists(4) = sp2ls_m0_redist
      num_redists = 4
      sp2ls_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    ! gp gather
    simple_gather_gp_2d_redist = new_redist(gp_gdeco(1), global_gp_2d_gdeco)
    simple_gather_gp_3d_redist = new_redist(gp_gdeco(nlev), global_gp_3d_gdeco)

    ! gp scatter
    simple_scatter_gp_2d_redist = new_redist(global_gp_2d_gdeco, gp_gdeco(1))
    simple_scatter_gp_3d_redist = new_redist(global_gp_3d_gdeco, gp_gdeco(nlev))

    IF (require_cross_gather) THEN
      IF (is_cross_gather_target) THEN
        simple_cross_gather_gp_2d_redist = new_redist(empty_deco, global_gp_2d_gdeco, model_comm)
        simple_cross_gather_gp_3d_redist = new_redist(empty_deco, global_gp_3d_gdeco, model_comm)
      ELSE
        simple_cross_gather_gp_2d_redist = new_redist(gp_gdeco(1), empty_deco, model_comm)
        simple_cross_gather_gp_3d_redist = new_redist(gp_gdeco(nlev), empty_deco, model_comm)
      ENDIF
    ENDIF

  CONTAINS

    TYPE(xt_redist) FUNCTION new_redist(src, dst, xmap_comm)
      TYPE(yaxt_deco), INTENT(in) :: src, dst
      INTEGER, INTENT(in), OPTIONAL :: xmap_comm
      TYPE(xt_xmap) :: x

      IF (PRESENT(xmap_comm)) THEN
        x = xt_xmap_all2all_new(src%idxlist, dst%idxlist, xmap_comm)
      ELSE
        x = xt_xmap_all2all_new(src%idxlist, dst%idxlist, ab_comm)
      ENDIF

      new_redist = xt_redist_p2p_off_new(x, src%offset, dst%offset, p_real_dp)

      CALL xt_xmap_delete(x)
    END FUNCTION new_redist

    SUBROUTINE my_big_ffsl_update_redists(nghost)
      USE mpi, ONLY: mpi_address_kind, MPI_SUCCESS
      INTEGER, INTENT(in) :: nghost
      INTEGER(mpi_address_kind) :: src_extent, dst_extent
      INTEGER(mpi_address_kind) :: base_address, temp_address
      INTEGER(mpi_address_kind) :: lb, extent
      INTEGER :: ierror, nrep, i
      INTEGER(c_int) :: displ(pcnst)

      CALL mpi_type_get_extent(p_real_dp, lb, extent, ierror)
      IF (ierror /= MPI_SUCCESS) CALL die('mpi_type_get_extent failed', __LINE__)
      src_extent =  ( nlon * (ffsl_gp_lat2-ffsl_gp_lat1+1 + 2*nghost) * ffsl_nlev ) *extent
      dst_extent =  src_extent
      nrep = pcnst
      displ = (/ (i, i=0,nrep-1) /)
      ffsl_update_4d_redist(nghost) = xt_redist_repeat_new(ffsl_update_3d_redist(nghost), src_extent, dst_extent, &
           & nrep, displ)
    END SUBROUTINE my_big_ffsl_update_redists

#if HAVE_YAXT_0_4_5
    ! requires at least yaxt version 0.4.5
    SUBROUTINE my_ffsl2gp_4d_redist
      USE mpi, ONLY: mpi_address_kind, MPI_SUCCESS
      INTEGER(mpi_address_kind) :: src_extent, dst_extent
      INTEGER(mpi_address_kind) :: base_address, temp_address
      INTEGER(mpi_address_kind) :: lb, dp_extent
      INTEGER :: ierror, nrep, i
      INTEGER(c_int) :: displ(pcnst)

      CALL mpi_type_get_extent(p_real_dp, lb, dp_extent, ierror)
      IF (ierror /= MPI_SUCCESS) CALL die('mpi_type_get_extent failed', __LINE__)
      src_extent =  ( nlon * (ffsl_gp_lat2-ffsl_gp_lat1+1) * ffsl_nlev ) * dp_extent
      dst_extent =  nproma * nlev * dp_extent
      nrep = pcnst
      displ = (/ (i, i=0,nrep-1) /)
      ffsl2gp_4d_fb_redist = xt_redist_repeat_new(ffsl2gp_3d_fb_single_tracer_redist, src_extent, dst_extent, &
           & nrep, displ)
    END SUBROUTINE my_ffsl2gp_4d_redist
#endif


#endif /* HAVE_YAXT */
  END SUBROUTINE generate_yaxt_redist

  SUBROUTINE yaxt_ffsl_ghost_update_2d(ffsl2d, nghost)
    REAL(dp), INTENT(inout) :: ffsl2d(:,:)
    INTEGER, INTENT(in) :: nghost
#ifdef HAVE_YAXT

    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_ffsl_ghost_update_2d'

#ifdef CHECK_SHAPES
    INTEGER :: sh(2)
    sh(:) = (/ nlon, ffsl_gp_lat2-ffsl_gp_lat1+1 + 2*nghost /)
    IF (ANY(SHAPE(ffsl2d) /= sh)) THEN
      WRITE(0,*) 'shape(ffsl2d)=',SHAPE(ffsl2d)
      WRITE(0,*) 'sh=',sh
      WRITE(0,*) 'nghost=',nghost
      WRITE(0,*) 'ffsl_gp_lat1,ffsl_gp_lat2=',ffsl_gp_lat1,ffsl_gp_lat2
      CALL die(subname//": unexpected data shape",__LINE__)
    ENDIF
#endif

    IF ( nghost /= 1 .AND. nghost /= 3 ) &
         & CALL die(subname//": bad value for nghost",__LINE__)
    CALL cont_s_exchange1_inplace(ffsl_update_2d_redist(nghost),ffsl2d)

#endif
  END SUBROUTINE yaxt_ffsl_ghost_update_2d

  SUBROUTINE yaxt_ffsl_ghost_update_3d(ffsl3d, nghost)
    REAL(dp), INTENT(inout) :: ffsl3d(:,:,:)
    INTEGER, INTENT(in) :: nghost
#ifdef HAVE_YAXT

    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_ffsl_ghost_update_3d'

#ifdef CHECK_SHAPES
    INTEGER :: sh(3)
    sh(:) = (/ nlon, ffsl_gp_lat2-ffsl_gp_lat1+1 + 2*nghost, ffsl_nlev /)
    IF (ANY(SHAPE(ffsl3d) /= sh)) THEN
      WRITE(0,*) 'shape(ffsl3d)=',SHAPE(ffsl3d)
      WRITE(0,*) 'sh=',sh
      WRITE(0,*) 'nghost=',nghost
      WRITE(0,*) 'ffsl_gp_lat1,ffsl_gp_lat2=',ffsl_gp_lat1,ffsl_gp_lat2
      CALL die(subname//": unexpected data shape",__LINE__)
    ENDIF
#endif

    IF ( nghost /= 1 .AND. nghost /= 3 ) &
         & CALL die(subname//": bad value for nghost",__LINE__)

    CALL cont_s_exchange1_inplace(ffsl_update_3d_redist(nghost),ffsl3d)

#endif
  END SUBROUTINE yaxt_ffsl_ghost_update_3d

  SUBROUTINE yaxt_ffsl_ghost_update_3d_repeated(m_ffsl3d, nghost)
    REAL(dp), INTENT(inout) :: m_ffsl3d(:,:,:,:)
    INTEGER, INTENT(in) :: nghost
#ifdef HAVE_YAXT

    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_ffsl_ghost_update_3d'
    INTEGER :: nt
#ifdef CHECK_SHAPES
    INTEGER :: sh(3)
    sh(:) = (/ nlon, ffsl_gp_lat2-ffsl_gp_lat1+1 + 2*nghost, ffsl_nlev /)
    IF (ANY(SHAPE(m_ffsl3d(:,:,:,1)) /= sh)) THEN
      WRITE(0,*) 'shape(ffsl3d)=',SHAPE(m_ffsl3d)
      WRITE(0,*) 'sh=',sh
      WRITE(0,*) 'nghost=',nghost
      WRITE(0,*) 'ffsl_gp_lat1,ffsl_gp_lat2=',ffsl_gp_lat1,ffsl_gp_lat2
      CALL die(subname//": unexpected data shape",__LINE__)
    ENDIF
#endif

    nt = SIZE(m_ffsl3d,4)
    IF (nt < 1 .OR. nt > yaxt_max_tracer_aggregate_size) CALL die(subname//": bad value for nt",__LINE__)
    IF (nghost /= 1 .AND. nghost /= 3 ) &
         & CALL die(subname//": bad value for nghost",__LINE__)

    IF (xt_is_null(ffsl_update_m3d_redist(nt, nghost))) CALL my_init

    CALL cont_s_exchange1_inplace(ffsl_update_m3d_redist(nt, nghost),m_ffsl3d)

  CONTAINS

    SUBROUTINE my_init
      USE mpi, ONLY: mpi_address_kind, MPI_SUCCESS
      INTEGER(mpi_address_kind) :: src_extent, dst_extent
      INTEGER(mpi_address_kind) :: lb, dp_extent
      INTEGER :: i, ierror
      INTEGER(c_int) :: displ(nt)
 
      CALL mpi_type_get_extent(p_real_dp, lb, dp_extent, ierror)
      IF (ierror /= MPI_SUCCESS) CALL die('mpi_type_get_extent failed', __LINE__)
      src_extent = ( nlon * (ffsl_gp_lat2-ffsl_gp_lat1+1 + 2*nghost) * ffsl_nlev ) * dp_extent
      dst_extent = src_extent
      displ = (/ (i, i=0,nt-1) /)
      ffsl_update_m3d_redist(nt, nghost) =  xt_redist_repeat_new(ffsl_update_3d_redist(nghost), src_extent, dst_extent, &
           & nt, displ)
    END SUBROUTINE my_init

#endif
  END SUBROUTINE yaxt_ffsl_ghost_update_3d_repeated


  SUBROUTINE yaxt_ffsl_ghost_update_4d(ffsl4d, nghost)
    REAL(dp), INTENT(inout) :: ffsl4d(:,:,:,:)
    INTEGER, INTENT(in) :: nghost
#ifdef HAVE_YAXT

    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_ffsl_ghost_update_4d'

#ifdef CHECK_SHAPES
    INTEGER :: sh(4)
    sh(:) = (/ nlon, ffsl_gp_lat2-ffsl_gp_lat1+1 + 2*nghost, ffsl_nlev, pcnst /)
    IF (ANY(SHAPE(ffsl4d) /= sh)) THEN
      WRITE(0,*) 'shape(ffsl4d)=',SHAPE(ffsl4d)
      WRITE(0,*) 'sh=',sh
      WRITE(0,*) 'nghost=',nghost
      WRITE(0,*) 'ffsl_gp_lat1,ffsl_gp_lat2=',ffsl_gp_lat1,ffsl_gp_lat2
      CALL die(subname//": unexpected data shape",__LINE__)
    ENDIF
#endif

    IF ( nghost /= 1 .AND. nghost /= 3 ) &
         & CALL die(subname//": bad value for nghost",__LINE__)

    CALL cont_s_exchange1_inplace(ffsl_update_4d_redist(nghost),ffsl4d)

#endif
  END SUBROUTINE yaxt_ffsl_ghost_update_4d

#ifdef HAVE_YAXT
  SUBROUTINE cont_s_exchange1(redist, src_field, dst_field)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(dp), TARGET, INTENT(in) :: src_field(*)
    REAL(dp), TARGET, INTENT(inout) :: dst_field(*)

    CALL xt_redist_s_exchange1(redist, C_LOC(src_field), C_LOC(dst_field))

  END SUBROUTINE cont_s_exchange1

  SUBROUTINE cont_s_exchange1_inplace(redist, field)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(dp), TARGET, INTENT(inout) :: field(*)
    TYPE(c_ptr) :: p
    p = C_LOC(field)
    CALL xt_redist_s_exchange1(redist, p, p)

  END SUBROUTINE cont_s_exchange1_inplace
#endif

  SUBROUTINE yaxt_simple_gather_gp_2d(global_gp, local_gp)
#ifndef LF
    REAL(dp), POINTER, INTENT(inout) :: global_gp(:,:)
#else
    REAL(dp), POINTER :: global_gp(:,:)
#endif
    REAL(dp), TARGET, INTENT(in) :: local_gp(:,:)

#ifdef HAVE_YAXT
    REAL(dp), SAVE :: dummy(1) = 0.0_dp
    INTEGER :: gshape(2), mshape(2)
    LOGICAL :: lsrc

    lsrc = (ab_rank == 0)
    mshape = SHAPE(local_gp)
    IF (ASSOCIATED(global_gp)) THEN
      gshape = SHAPE(global_gp)
    ELSE
      gshape = 0
    ENDIF
    IF (ANY(mshape /= (/ nproma, ngpblks /) )) &
         & CALL die('yaxt_simple_gather_gp_2d: unexpected shape of local array', __LINE__)

    IF (lsrc) THEN
      IF (ANY(gshape /= (/ nlon, nlat /) )) &
           & CALL die('yaxt_simple_gather_gp_2d: unexpected shape of global array', __LINE__)
      CALL cont_s_exchange1(simple_gather_gp_2d_redist, local_gp, global_gp)
    ELSE
      CALL cont_s_exchange1(simple_gather_gp_2d_redist, local_gp, dummy)
    ENDIF
#endif
  END SUBROUTINE yaxt_simple_gather_gp_2d

  SUBROUTINE yaxt_simple_scatter_gp_2d(global_gp, local_gp)
#ifndef LF
    REAL(dp), POINTER, INTENT(in) :: global_gp(:,:)
#else
    REAL(dp), POINTER :: global_gp(:,:)
#endif
    REAL(dp), TARGET, INTENT(out) :: local_gp(:,:)
#ifdef HAVE_YAXT
    REAL(dp), SAVE, TARGET :: dummy(1) = 0.0_dp
    INTEGER :: gshape(2), mshape(2)
    LOGICAL :: lsrc

    lsrc = (ab_rank == 0)
    mshape = SHAPE(local_gp)

    IF (ANY(mshape /= (/ nproma, ngpblks /) )) &
         & CALL die('yaxt_simple_scatter_gp_2d: unexpected shape of local array', __LINE__)

    IF (ASSOCIATED(global_gp)) THEN
      gshape = SHAPE(global_gp)
    ELSE
      gshape = 0
    ENDIF

    IF (lsrc) THEN
      IF (ANY(gshape /= (/ nlon, nlat /) )) &
           & CALL die('yaxt_simple_scatter_gp_2d: unexpected shape of global array', __LINE__)
      CALL cont_s_exchange1(simple_scatter_gp_2d_redist, global_gp, local_gp)
    ELSE
      CALL cont_s_exchange1(simple_scatter_gp_2d_redist, dummy, local_gp)
    ENDIF

    IF (nproma /= npromz) CALL zero_fractional_gp_block_2d(local_gp)
#endif
  END SUBROUTINE yaxt_simple_scatter_gp_2d

  SUBROUTINE yaxt_simple_gather_gp_3d(global_gp, local_gp)
#ifndef LF
    REAL(dp), POINTER, INTENT(inout) :: global_gp(:,:,:)
#else
    REAL(dp), POINTER :: global_gp(:,:,:)
#endif
    REAL(dp), TARGET, INTENT(in) :: local_gp(:,:,:)
#ifdef HAVE_YAXT
    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_simple_gather_gp_3d'
    REAL(dp), SAVE :: dummy(1) = 0.0_dp
    INTEGER :: gshape(3), mshape(3)
    LOGICAL :: lsrc

    lsrc = (ab_rank == 0)
    mshape = SHAPE(local_gp)
    IF (mshape(2) /= nlev) THEN
      CALL simple_fallback
      RETURN
    ENDIF
    IF (ANY(mshape /= (/ nproma, nlev, ngpblks /) )) &
         & CALL die(subname//': unexpected shape of local array', __LINE__)
    IF (ASSOCIATED(global_gp)) THEN
      gshape = SHAPE(global_gp)
    ELSE
      gshape = 0
    ENDIF


    IF (lsrc) THEN
      IF (ANY(gshape /= (/ nlon, nlev, nlat /) )) &
           & CALL die(subname//': unexpected shape of global array', __LINE__)
      CALL cont_s_exchange1(simple_gather_gp_3d_redist, local_gp, global_gp)
    ELSE
      CALL cont_s_exchange1(simple_gather_gp_3d_redist, local_gp, dummy)
    ENDIF

  CONTAINS

    SUBROUTINE simple_fallback
      REAL(dp), ALLOCATABLE :: gtmp(:,:)
      REAL(dp) :: mtmp(nproma,ngpblks)
      INTEGER :: k

      IF (lsrc) THEN
        ALLOCATE(gtmp(nlon,nlat))
      ELSE
        ALLOCATE(gtmp(1,1))
      ENDIF
      DO k = 1, mshape(2)
        IF (lsrc) mtmp = local_gp(:,k,:)
        CALL cont_s_exchange1(simple_gather_gp_2d_redist, mtmp, gtmp)
        IF (lsrc) global_gp(:,k,:) = gtmp
      ENDDO
    END SUBROUTINE simple_fallback
#endif
  END SUBROUTINE yaxt_simple_gather_gp_3d

  SUBROUTINE yaxt_simple_scatter_gp_3d(global_gp, local_gp)
#ifndef LF
    REAL(dp), POINTER, INTENT(in) :: global_gp(:,:,:)
#else
    REAL(dp), POINTER :: global_gp(:,:,:)
#endif
    REAL(dp), TARGET, INTENT(out) :: local_gp(:,:,:)
#ifdef HAVE_YAXT
    REAL(dp), SAVE :: dummy(1) = 0.0_dp
    INTEGER :: gshape(3), mshape(3)
    LOGICAL :: lsrc

    lsrc = (ab_rank == 0)

    mshape = SHAPE(local_gp)
    IF (mshape(2) /= nlev) THEN
      CALL simple_fallback
      RETURN
    ENDIF

    IF (ANY(mshape /= (/ nproma, nlev, ngpblks /) )) &
         & CALL die('yaxt_simple_scatter_gp_3d: unexpected shape of local array', __LINE__)

    IF (ASSOCIATED(global_gp)) THEN
      gshape = SHAPE(global_gp)
    ELSE
      gshape = 0
    ENDIF

    IF (lsrc) THEN
      IF (ANY(gshape /= (/ nlon, nlev, nlat /) )) &
           & CALL die('yaxt_simple_scatter_gp_3denerate_yaxt_redist: unexpected shape of global array', __LINE__)
      CALL cont_s_exchange1(simple_scatter_gp_3d_redist, global_gp, local_gp)
    ELSE
      CALL cont_s_exchange1(simple_scatter_gp_3d_redist, dummy, local_gp)
    ENDIF

    IF (nproma /= npromz) CALL zero_fractional_gp_block_3d(local_gp)

  CONTAINS

    SUBROUTINE simple_fallback
      REAL(dp), ALLOCATABLE :: gtmp(:,:)
      REAL(dp) :: mtmp(nproma,ngpblks)
      INTEGER :: k

      IF (lsrc) THEN
        ALLOCATE(gtmp(nlon,nlat))
      ELSE
        ALLOCATE(gtmp(1,1))
        gtmp = 0.0_dp
      ENDIF
      DO k = 1, mshape(2)
        IF (lsrc) gtmp = global_gp(:,k,:)
        CALL cont_s_exchange1(simple_scatter_gp_2d_redist, gtmp, mtmp)
        local_gp(:,k,:) = mtmp
      ENDDO
    END SUBROUTINE simple_fallback
#endif
  END SUBROUTINE yaxt_simple_scatter_gp_3d

  SUBROUTINE yaxt_simple_cross_gather_gp_2d(global_gp, local_gp)
#ifndef LF
    REAL(dp), POINTER, INTENT(inout) :: global_gp(:,:)
#else
    REAL(dp), POINTER :: global_gp(:,:)
#endif
    REAL(dp), TARGET, INTENT(in) :: local_gp(:,:)
    INTEGER :: i,j
#ifdef HAVE_YAXT
    REAL(dp), SAVE, TARGET :: dummy(1) = 0.0_dp
    INTEGER :: gshape(2), mshape(2)
    LOGICAL :: lsrc

    lsrc = (ab_rank == 0)
    mshape = SHAPE(local_gp)
    IF (ASSOCIATED(global_gp)) THEN
      gshape = SHAPE(global_gp)
    ELSE
      gshape = 0
    ENDIF
    IF (ANY(mshape /= (/ nproma, ngpblks /) )) &
         & CALL die('yaxt_simple_cross_gather_gp_2d: unexpected shape of local array', __LINE__)

    IF (lsrc) THEN
      IF (ANY(gshape /= (/ nlon, nlat /) )) &
           & CALL die('yaxt_simple_cross_gather_gp_2d: unexpected shape of global array', __LINE__)
      CALL cont_s_exchange1(simple_cross_gather_gp_2d_redist, local_gp, global_gp)
    ELSE
      CALL cont_s_exchange1(simple_cross_gather_gp_2d_redist, dummy, global_gp)
    ENDIF
#endif
  END SUBROUTINE yaxt_simple_cross_gather_gp_2d

  SUBROUTINE yaxt_simple_cross_gather_gp_3d(global_gp, local_gp)
#ifndef LF
    REAL(dp), POINTER, INTENT(inout) :: global_gp(:,:,:)
#else
    REAL(dp), POINTER :: global_gp(:,:,:)
#endif
    REAL(dp), TARGET, INTENT(in) :: local_gp(:,:,:)
#ifdef HAVE_YAXT
    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_cross_simple_gather_gp_3d'
    REAL(dp), SAVE :: dummy(1) = 0.0_dp
    INTEGER :: gshape(3), mshape(3)
    LOGICAL :: lsrc

    lsrc = (ab_rank == 0)
    mshape = SHAPE(local_gp)
    IF (mshape(2) /= nlev) THEN
      CALL simple_fallback
      RETURN
    ENDIF
    IF (ANY(mshape /= (/ nproma, nlev, ngpblks /) )) &
         & CALL die(subname//': unexpected shape of local array', __LINE__)
    IF (ASSOCIATED(global_gp)) THEN
      gshape = SHAPE(global_gp)
    ELSE
      gshape = 0
    ENDIF

    IF (lsrc) THEN
      IF (ANY(gshape /= (/ nlon, nlev, nlat /) )) &
           & CALL die(subname//': unexpected shape of global array', __LINE__)
      CALL cont_s_exchange1(simple_cross_gather_gp_3d_redist, local_gp, global_gp)
    ELSE
      CALL cont_s_exchange1(simple_cross_gather_gp_3d_redist, local_gp, dummy)
    ENDIF

  CONTAINS

    SUBROUTINE simple_fallback
      REAL(dp), ALLOCATABLE :: gtmp(:,:)
      REAL(dp) :: mtmp(nproma,ngpblks)
      INTEGER :: k

      IF (lsrc) THEN
        ALLOCATE(gtmp(nlon,nlat))
      ELSE
        ALLOCATE(gtmp(1,1))
        gtmp = 0.0_dp
      ENDIF
      mtmp = 0.0_dp ! todo: we only need to zero the fractional last block
      DO k = 1, mshape(2)
        IF (lsrc) mtmp = local_gp(:,k,:)
        CALL cont_s_exchange1(simple_cross_gather_gp_2d_redist, mtmp, gtmp)
        IF (lsrc) global_gp(:,k,:) = gtmp
      ENDDO
    END SUBROUTINE simple_fallback
#endif
  END SUBROUTINE yaxt_simple_cross_gather_gp_3d

  SUBROUTINE yaxt_tr_gp_to_ffsl_2d(gp2d1, ffsl2d1)
    REAL(dp), TARGET, INTENT(in)  :: gp2d1(:,:)     ! decomposed with 2d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(out) :: ffsl2d1(:,:)   ! decomposed with 2d ffsl gridpoint deco.
#ifdef HAVE_YAXT
    CALL timer_start(it_gp2ffsl)
    IF (   ANY(SHAPE(gp2d1)   /= (/ nproma, ngpblks /) ) .OR. &
         & ANY(SHAPE(ffsl2d1) /= (/ nlon, ffsl_nlat /) ) ) &
         & CALL die("yaxt_tr_gp_to_ffsl_2d: shape mismatch",__LINE__)

    CALL cont_s_exchange1(gp2ffsl_2d_redist, gp2d1, ffsl2d1)
    CALL timer_stop(it_gp2ffsl)
#endif
  END SUBROUTINE yaxt_tr_gp_to_ffsl_2d

  SUBROUTINE yaxt_tr_gp_to_ffsl_3d(gp3d, ffsl3d)
    REAL(dp), TARGET, INTENT(in)  :: gp3d(:,:,:)   ! decomposed with 3d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(out) :: ffsl3d(:,:,:) ! decomposed with 3d ffsl deco.
#ifdef HAVE_YAXT
    CALL timer_start(it_gp2ffsl)
    IF (   ANY(SHAPE(gp3d)   /= (/ nproma, nlev, ngpblks /) ) .OR. &
         & ANY(SHAPE(ffsl3d) /= (/ nlon, ffsl_nlat, ffsl_nlev /) ) ) &
         & CALL die("yaxt_tr_gp_to_ffsl_3d: shape mismatch",__LINE__)

    CALL cont_s_exchange1(gp2ffsl_3d_redist, gp3d, ffsl3d)
    CALL timer_stop(it_gp2ffsl)
#endif
  END SUBROUTINE yaxt_tr_gp_to_ffsl_3d

  SUBROUTINE yaxt_tr_gp_to_ffsl_3d_repeated(m_gp3d, m_ffsl3d)
    REAL(dp), TARGET, INTENT(in)  :: m_gp3d(:,:,:,:)   ! decomposed with 3d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(out) :: m_ffsl3d(:,:,:,:) ! decomposed with 3d ffsl deco.
#ifdef HAVE_YAXT
    INTEGER :: it, nt
    CALL timer_start(it_gp2ffsl)
    nt = SIZE(m_gp3d,4)

    IF (SIZE(m_ffsl3d,4) /= nt) CALL die("yaxt_tr_gp_to_ffsl_3d_repeated: size mismatch",__LINE__)
    IF (   ANY(SHAPE(m_gp3d(:,:,:,1))   /= (/ nproma, nlev, ngpblks /) ) .OR. &
         & ANY(SHAPE(m_ffsl3d(:,:,:,1)) /= (/ nlon, ffsl_nlat, ffsl_nlev /) ) ) &
         & CALL die("yaxt_tr_gp_to_ffsl_3d_repeated: shape mismatch",__LINE__)
    IF (nt > yaxt_max_tracer_aggregate_size) &
         & CALL die("yaxt_tr_gp_to_ffsl_3d_repeated: (nt > yaxt_max_tracer_aggregate_size)",__LINE__)
    IF (xt_is_null(gp2ffsl_m3d_redist(nt))) CALL my_init

    CALL cont_s_exchange1(gp2ffsl_m3d_redist(nt), m_gp3d(:,:,:,:), m_ffsl3d(:,:,:,:))
    CALL timer_stop(it_gp2ffsl)
  CONTAINS

    SUBROUTINE my_init
      USE mpi, ONLY: mpi_address_kind, MPI_SUCCESS
      INTEGER(mpi_address_kind) :: src_extent, dst_extent
      INTEGER(mpi_address_kind) :: lb, dp_extent
      INTEGER :: i, nrep, ierror
      INTEGER(c_int) :: displ(nt)

      CALL mpi_type_get_extent(p_real_dp, lb, dp_extent, ierror)
      IF (ierror /= MPI_SUCCESS) CALL die('mpi_type_get_extent failed', __LINE__)
      src_extent =  nproma * nlev * ngpblks * dp_extent
      dst_extent = nlon * ffsl_nlat * ffsl_nlev * dp_extent
      nrep = nt
      displ = (/ (i, i=0,nrep-1) /)
      gp2ffsl_m3d_redist(nrep) =  xt_redist_repeat_new(gp2ffsl_3d_redist, src_extent, dst_extent, &
           & nrep, displ)
    END SUBROUTINE my_init

#endif
  END SUBROUTINE yaxt_tr_gp_to_ffsl_3d_repeated

#ifdef HAVE_YAXT
  SUBROUTINE cont_tr_gp_to_ffsl(gp3d1, gp3d2, gp3d3, gp3d4, gp3d5, gp2d1, gp2d2, &
       &          ffsl3d1, ffsl3d2, ffsl3d3, ffsl3d4, ffsl3d5, ffsl2d1, ffsl2d2)
    REAL(dp), TARGET, INTENT(in)  :: gp3d1(*)
    REAL(dp), TARGET, INTENT(in)  :: gp3d2(*)
    REAL(dp), TARGET, INTENT(in)  :: gp3d3(*)
    REAL(dp), TARGET, INTENT(in)  :: gp3d4(*)
    REAL(dp), TARGET, INTENT(in)  :: gp3d5(*)
    REAL(dp), TARGET, INTENT(in)  :: gp2d1(*)
    REAL(dp), TARGET, INTENT(in)  :: gp2d2(*)
    REAL(dp), TARGET, INTENT(out) :: ffsl3d1(*)
    REAL(dp), TARGET, INTENT(out) :: ffsl3d2(*)
    REAL(dp), TARGET, INTENT(out) :: ffsl3d3(*)
    REAL(dp), TARGET, INTENT(out) :: ffsl3d4(*)
    REAL(dp), TARGET, INTENT(out) :: ffsl3d5(*)
    REAL(dp), TARGET, INTENT(out) :: ffsl2d1(*)
    REAL(dp), TARGET, INTENT(out) :: ffsl2d2(*)
    TYPE(c_ptr) :: gp_c_ptr(7), ffsl_c_ptr(7)

    IF (aggregate_gp_to_ffsl) THEN
      gp_c_ptr = (/ C_LOC(gp3d1), C_LOC(gp3d2), C_LOC(gp3d3), &
           &        C_LOC(gp3d4), C_LOC(gp3d5), &
           &        C_LOC(gp2d1), C_LOC(gp2d2) /)
      ffsl_c_ptr = (/ C_LOC(ffsl3d1), C_LOC(ffsl3d2), C_LOC(ffsl3d3), &
           &          C_LOC(ffsl3d4), C_LOC(ffsl3d5), &
           &          C_LOC(ffsl2d1), C_LOC(ffsl2d2) /)
      CALL xt_redist_s_exchange(gp2ffsl_all_redist, gp_c_ptr, ffsl_c_ptr)
    ELSE
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d1), C_LOC(ffsl3d1) )
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d2), C_LOC(ffsl3d2) )
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d3), C_LOC(ffsl3d3) )
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d4), C_LOC(ffsl3d4) )
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d5), C_LOC(ffsl3d5) )
      CALL xt_redist_s_exchange1(gp2ffsl_2d_redist, C_LOC(gp2d1), C_LOC(ffsl2d1) )
      CALL xt_redist_s_exchange1(gp2ffsl_2d_redist, C_LOC(gp2d2), C_LOC(ffsl2d2) )
    ENDIF

  END SUBROUTINE cont_tr_gp_to_ffsl
#endif

  SUBROUTINE yaxt_tr_gp_to_ffsl(gp3d1, gp3d2, gp3d3, gp3d4, gp3d5, gp2d1, gp2d2, &
       &                        ffsl3d1, ffsl3d2, ffsl3d3, ffsl3d4, ffsl3d5, ffsl2d1, ffsl2d2)
    REAL(dp), INTENT(in)  :: gp3d1(:,:,:)   ! decomposed with 3d standard gridpoint deco.
    REAL(dp), INTENT(in)  :: gp3d2(:,:,:)   ! ""
    REAL(dp), INTENT(in)  :: gp3d3(:,:,:)   ! ""
    REAL(dp), INTENT(in)  :: gp3d4(:,:,:)   ! ""
    REAL(dp), INTENT(in)  :: gp3d5(:,:,:)   ! ""
    REAL(dp), INTENT(in)  :: gp2d1(:,:)     ! decomposed with 2d standard gridpoint deco.
    REAL(dp), INTENT(in)  :: gp2d2(:,:)     ! ""
    REAL(dp), INTENT(out) :: ffsl3d1(:,:,:) ! decomposed with 3d ffsl deco.
    REAL(dp), INTENT(out) :: ffsl3d2(:,:,:) ! ""
    REAL(dp), INTENT(out) :: ffsl3d3(:,:,:) ! ""
    REAL(dp), INTENT(out) :: ffsl3d4(:,:,:) ! ""
    REAL(dp), INTENT(out) :: ffsl3d5(:,:,:) ! ""
    REAL(dp), INTENT(out) :: ffsl2d1(:,:)   ! decomposed with 2d ffsl gridpoint deco.
    REAL(dp), INTENT(out) :: ffsl2d2(:,:)   ! ""
#ifdef HAVE_YAXT
    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_tr_gp_to_ffsl'

#ifdef CHECK_SHAPES
    INTEGER :: gp3_shape(3), gp2_shape(2), ffsl3_shape(3), ffsl2_shape(2)
#endif
    CALL timer_start(it_gp2ffsl)
#ifdef CHECK_SHAPES
    gp3_shape = (/ nproma, nlev, ngpblks /)
    gp2_shape = (/ nproma, ngpblks /)
    ffsl3_shape = (/ nlon, ffsl_nlat, ffsl_nlev /)
    ffsl2_shape = (/ nlon, ffsl_nlat /)
    IF (nproma * nlev * ngpblks * nlon * ffsl_nlat * ffsl_nlev == 0) &
         & CALL die(subname//": found unexpected zero shape product",__LINE__)
    IF (   ANY(SHAPE(gp3d1) /= gp3_shape) .OR. &
         & ANY(SHAPE(gp3d2) /= gp3_shape) .OR. &
         & ANY(SHAPE(gp3d3) /= gp3_shape) .OR. &
         & ANY(SHAPE(gp3d4) /= gp3_shape) .OR. &
         & ANY(SHAPE(gp3d5) /= gp3_shape) ) &
         & CALL die(subname//": gp3_shape mismatch",__LINE__)
    IF (   ANY(SHAPE(gp2d1) /= gp2_shape) .OR. &
         & ANY(SHAPE(gp2d2) /= gp2_shape) ) &
         & CALL die(subname//": gp2_shape mismatch",__LINE__)
    IF (   ANY(SHAPE(ffsl3d1) /= ffsl3_shape) .OR. &
         & ANY(SHAPE(ffsl3d2) /= ffsl3_shape) .OR. &
         & ANY(SHAPE(ffsl3d3) /= ffsl3_shape) .OR. &
         & ANY(SHAPE(ffsl3d4) /= ffsl3_shape) .OR. &
         & ANY(SHAPE(ffsl3d5) /= ffsl3_shape) ) &
         & CALL die(subname//": ffsl3_shape mismatch",__LINE__)
    IF (   ANY(SHAPE(ffsl2d1) /= ffsl2_shape) .OR. &
         & ANY(SHAPE(ffsl2d2) /= ffsl2_shape) ) &
         & CALL die(subname//": ffsl3_shape mismatch",__LINE__)
#endif

    CALL cont_tr_gp_to_ffsl(gp3d1, gp3d2, gp3d3, gp3d4, gp3d5, gp2d1, gp2d2, &
         &  ffsl3d1, ffsl3d2, ffsl3d3, ffsl3d4, ffsl3d5, ffsl2d1, ffsl2d2)
    CALL timer_stop(it_gp2ffsl)
#endif
  END SUBROUTINE yaxt_tr_gp_to_ffsl

  SUBROUTINE yaxt_tr_ffsl_to_gp_3d(gp3d, ffsl3d)
    REAL(dp), TARGET, INTENT(out) :: gp3d(:,:,:)     ! decomposed with 3d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(in)  :: ffsl3d(:,:,:)   ! decomposed with 3d ffsl deco.

#ifdef HAVE_YAXT
    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_tr_ffsl_to_gp_3d'

    INTEGER :: gp3d_shape(3), ffsl3d_shape(3)
    CALL timer_start(it_ffsl2gp)
#ifdef CHECK_SHAPES
    gp3d_shape = (/ nproma, nlev, ngpblks /)
    ffsl3d_shape = (/ nlon, ffsl_nlat, ffsl_nlev /)
    IF (nproma * nlev * ngpblks * nlon * ffsl_nlat * ffsl_nlev == 0) &
         & CALL die(subname//": found unexpected zero shape product",__LINE__)
    IF (ANY(SHAPE(gp3d) /= gp3d_shape)) CALL die(subname//": gp3d_shape mismatch",__LINE__)
    IF (ANY(SHAPE(ffsl3d) /= ffsl3d_shape)) CALL die(subname//": ffsl3d_shape mismatch",__LINE__)
#endif

    CALL cont_s_exchange1(ffsl2gp_3d_redist, ffsl3d, gp3d)
    IF (npromz /= nproma) CALL zero_fractional_gp_block_3d(gp3d)
    CALL timer_stop(it_ffsl2gp)
#endif
  END SUBROUTINE yaxt_tr_ffsl_to_gp_3d

#ifdef HAVE_YAXT
  SUBROUTINE cont_tr_ffsl_to_gp(gp4d, gp2d, gp3d, ffsl4d, ffsl2d, ffsl3d)
    REAL(dp), TARGET, INTENT(out) :: gp4d(nproma, nlev, pcnst, ngpblks)
    REAL(dp), TARGET, INTENT(out) :: gp2d(nproma, ngpblks)
    REAL(dp), TARGET, INTENT(out) :: gp3d(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(in)  :: ffsl4d(nlon, ffsl_nlat, ffsl_nlev, pcnst)
    REAL(dp), TARGET, INTENT(in)  :: ffsl2d(nlon, ffsl_nlat)
    REAL(dp), TARGET, INTENT(in)  :: ffsl3d(nlon, ffsl_nlat, ffsl_nlev)

    TYPE(c_ptr) :: gp_c_ptr(3), ffsl_c_ptr(3)
    INTEGER :: ib, it

    IF (aggregate_ffsl_to_gp) THEN
      ffsl_c_ptr = (/ C_LOC(ffsl4d), C_LOC(ffsl2d), C_LOC(ffsl3d) /)
      gp_c_ptr = (/  C_LOC(gp4d), C_LOC(gp2d), C_LOC(gp3d) /)
      CALL xt_redist_s_exchange(ffsl2gp_all_redist, ffsl_c_ptr, gp_c_ptr)
    ELSE
      DO ib = 1, pcnst
        it = ib
        CALL xt_redist_s_exchange1(ffsl2gp_3d_fb_single_tracer_redist,C_LOC(ffsl4d(1,1,1,ib)),C_LOC(gp4d(1,1,it,1)))
      ENDDO

      CALL xt_redist_s_exchange1(ffsl2gp_2d_redist, C_LOC(ffsl2d(1,1)), C_LOC(gp2d(1,1)))
      CALL xt_redist_s_exchange1(ffsl2gp_3d_redist, C_LOC(ffsl3d(1,1,1)), C_LOC(gp3d(1,1,1)))
    ENDIF

  END SUBROUTINE cont_tr_ffsl_to_gp
#endif

  SUBROUTINE yaxt_tr_ffsl_to_gp(gp4d, gp2d, gp3d, ffsl4d, ffsl2d, ffsl3d)
    REAL(dp), TARGET, INTENT(out) :: gp4d(:,:,:,:)   ! decomposed with 4d standard (gridpoint&tracer) deco (tpcore::fb_gp).
    REAL(dp), TARGET, INTENT(out) :: gp2d(:,:)       ! decomposed with 2d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(out) :: gp3d(:,:,:)     ! decomposed with 3d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(in)  :: ffsl4d(:,:,:,:) ! decomposed with 4d ffsl deco (tpcore::fb).
    REAL(dp), TARGET, INTENT(in)  :: ffsl2d(:,:)     ! decomposed with 2d ffsl deco.
    REAL(dp), TARGET, INTENT(in)  :: ffsl3d(:,:,:)   ! decomposed with 3d ffsl deco.
#ifdef HAVE_YAXT
    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_tr_ffsl_to_gp'

    INTEGER :: gp4d_shape(4), gp3d_shape(3), gp2d_shape(2)
    INTEGER :: ffsl4d_shape(4), ffsl3d_shape(3), ffsl2d_shape(2)
    CALL timer_start(it_ffsl2gp)
#ifdef CHECK_SHAPES
    gp4d_shape = (/ nproma, nlev, pcnst, ngpblks /)
    gp3d_shape = (/ nproma, nlev, ngpblks /)
    gp2d_shape = (/ nproma, ngpblks /)
    ffsl4d_shape = (/ nlon, ffsl_nlat, ffsl_nlev, pcnst /)
    ffsl3d_shape = (/  nlon, ffsl_nlat, ffsl_nlev /)
    ffsl2d_shape = (/ nlon, ffsl_nlat /)
    IF (ANY(SHAPE(gp4d) /= gp4d_shape)) THEN
      WRITE(0,*) 'gp4d_shape=',gp4d_shape
      WRITE(0,*) 'shape(gp4d)=',SHAPE(gp4d)
      CALL die(subname//": gp4d shape mismatch",__LINE__)
    ENDIF
    IF (ANY(SHAPE(gp3d) /= gp3d_shape)) &
         & CALL die(subname//": gp3d shape mismatch",__LINE__)
    IF (ANY(SHAPE(gp2d) /= gp2d_shape)) &
         & CALL die(subname//": gp2d shape mismatch",__LINE__)
    IF (ANY(SHAPE(ffsl4d) /= ffsl4d_shape)) &
         & CALL die(subname//": ffsl4d shape mismatch",__LINE__)
    IF (ANY(SHAPE(ffsl4d) /= ffsl4d_shape)) &
         & CALL die(subname//": ffsl3d shape mismatch",__LINE__)
    IF (ANY(SHAPE(ffsl2d) /= ffsl2d_shape)) &
         & CALL die(subname//": ffsl2d shape mismatch",__LINE__)
#endif

    CALL cont_tr_ffsl_to_gp(gp4d, gp2d, gp3d, ffsl4d, ffsl2d, ffsl3d)

    IF (npromz /= nproma) THEN
      CALL zero_fractional_gp_block_4d(gp4d)
      CALL zero_fractional_gp_block_3d(gp3d)
      CALL zero_fractional_gp_block_2d(gp2d)
    ENDIF
    CALL timer_stop(it_ffsl2gp)
#endif
  END SUBROUTINE yaxt_tr_ffsl_to_gp

#ifdef HAVE_YAXT
  SUBROUTINE cont_tr_gp_to_fs(gp1, gp2, gp3, gp4, gp5, gp6, gp7, &
       &                      sf3, zm1, zm2, zm3, fs, fs0)
    !
    !   grid point space  -> Fourier space
    !
    REAL(dp), TARGET, INTENT(in)    :: gp1(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(in)    :: gp2(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(in)    :: gp3(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(in)    :: gp4(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(in)    :: gp5(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(in)    :: gp6(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(in)    :: gp7(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(in)    :: sf3(nproma, ngpblks)
    REAL(dp), TARGET, INTENT(in)    :: zm1(nlev,*)
    REAL(dp), TARGET, INTENT(in)    :: zm2(nlev,*)
    REAL(dp), TARGET, INTENT(in)    :: zm3(nlev,*)
    REAL(dp), TARGET, INTENT(inout) :: fs(nlon+2, nflevp1, nflat,*)
    REAL(dp), TARGET, INTENT(inout) :: fs0(nflev, nflat, *)

    TYPE(c_ptr) :: gp_c_ptr(7)
    TYPE(c_ptr) :: zm_c_ptr(3)
    TYPE(c_ptr) :: all_src_c_ptr(11), all_dst_c_ptr(11)
    INTEGER :: m

    CALL zero_fs

    gp_c_ptr = (/ C_LOC(gp1), C_LOC(gp2), C_LOC(gp3), &
         &        C_LOC(gp4), C_LOC(gp5), C_LOC(gp6), &
         &        C_LOC(gp7)  /)

    zm_c_ptr = (/ C_LOC(zm1), C_LOC(zm2), C_LOC(zm3) /)

    IF (aggregate_gp_to_fs) THEN
      DO m = 1,7
        all_src_c_ptr(m) =  gp_c_ptr(m)
        all_dst_c_ptr(m) =  C_LOC(fs(1,1,1,m))
      ENDDO
      all_src_c_ptr(8) =  C_LOC(sf3)
      all_dst_c_ptr(8) =  C_LOC(fs(1,nflevp1,1,3))
      DO m = 1,3
        all_src_c_ptr(8+m) =  zm_c_ptr(m)
        all_dst_c_ptr(8+m) =  C_LOC(fs0(1,1,m))
      ENDDO
      CALL xt_redist_s_exchange(gp2fs_all_redist, all_src_c_ptr, all_dst_c_ptr)
    ELSE
      DO m = 1, 7
        CALL xt_redist_s_exchange1(gp2fs_3d_redist, gp_c_ptr(m), C_LOC(fs(1,1,1,m)));
      ENDDO
      CALL xt_redist_s_exchange1(gp2fs_2d_redist, C_LOC(sf3), C_LOC(fs(1,nflevp1,1,3)))
      DO m = 1, 3
        CALL xt_redist_s_exchange1(gp2fs_zm_redist, zm_c_ptr(m), C_LOC(fs0(1,1,m)))
      ENDDO
    ENDIF

  CONTAINS

    ! zero selected 4d fs data
    SUBROUTINE zero_fs
      INTEGER :: i, jj, kk, p

      DO p = 1, 7
        DO jj = 1, nflat
          DO kk = 1, nflevp1
            DO i = nlon+1, nlon+2
              fs(i,kk,jj,p) = 0.0_dp
            ENDDO
          ENDDO
          DO kk = nflev+1,nflevp1
            DO i = 1, nlon+2
              fs(i,kk,jj,p) = 0.0_dp
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE zero_fs

  END SUBROUTINE cont_tr_gp_to_fs
#endif

  SUBROUTINE yaxt_tr_gp_to_fs(gp1, gp2, gp3, gp4, gp5, gp6, gp7, &
       &                      sf3, zm1, zm2, zm3, fs, fs0)
    !
    !   grid point space  -> Fourier space
    !
    REAL(dp), TARGET, INTENT(in)    :: gp1(:,:,:)   ! gridpoint space 3d
    REAL(dp), TARGET, INTENT(in)    :: gp2(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp3(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp4(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp5(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp6(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp7(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: sf3(:,:)     ! gridpoint space 2d
    REAL(dp), TARGET, INTENT(in)    :: zm1(:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(in)    :: zm2(:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(in)    :: zm3(:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(inout) :: fs (:,:,:,:) ! Fourier space
    REAL(dp), TARGET, INTENT(inout) :: fs0(:,:,:)   ! zonal mean, Four.
#ifdef HAVE_YAXT
    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_tr_gp_to_fs'

#ifdef CHECK_SHAPES
    INTEGER :: nglat
    INTEGER :: gp3d_shape(3), gp2d_shape(2), zm_shape(2), fs3d_shape(3), fs0_2d_shape(2)
    INTEGER :: tmp_4d_shape(4), tmp_3d_shape(3)
    nglat = ldc%nglat
    gp3d_shape = (/ nproma, nlev, ngpblks /)
    gp2d_shape = (/ nproma, ngpblks /)
    zm_shape = (/ nlev, nglat /)
    fs3d_shape = (/ nlon+2, nflevp1, nflat /)
    fs0_2d_shape = (/ nflev, nflat /)
    IF (nproma * nlev * ngpblks * (nlon+2) * nflev * nflat == 0) &
         & CALL die(subname//": found unexpected zero shape product",__LINE__)
    IF (   ANY(SHAPE(gp1) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp2) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp3) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp4) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp5) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp6) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp7) /= gp3d_shape) ) &
         & CALL die(subname//": gp3d_shape mismatch",__LINE__)
    IF (   ANY(SHAPE(sf3) /= gp2d_shape) ) &
         & CALL die(subname//": gp2d_shape mismatch",__LINE__)
    IF (   ANY(SHAPE(zm1) /= zm_shape) .OR. &
         & ANY(SHAPE(zm2) /= zm_shape) .OR. &
         & ANY(SHAPE(zm3) /= zm_shape) ) &
         & CALL die(subname//": zm_shape mismatch",__LINE__)
    tmp_4d_shape = SHAPE(fs)
    IF (   ANY(tmp_4d_shape(1:3) /= fs3d_shape) ) &
         & CALL die(subname//": fs3d_shape mismatch",__LINE__)
    tmp_3d_shape = SHAPE(fs0)
    IF (   ANY(tmp_3d_shape(1:2) /= fs0_2d_shape) ) &
         & CALL die(subname//": fs0_2d_shape mismatch",__LINE__)
#endif

    CALL cont_tr_gp_to_fs(gp1, gp2, gp3, gp4, gp5, gp6, gp7, &
         &                      sf3, zm1, zm2, zm3, fs, fs0)
#endif
  END SUBROUTINE yaxt_tr_gp_to_fs

#ifdef HAVE_YAXT
  SUBROUTINE cont_tr_fs_to_gp(gp1, gp2, gp3, gp4, gp5, gp6, gp7, gp8, gp9, &
       &                      sf1, sf2, sf3, zm1, zm2, zm3, &
       &                      fs, fs0)
    !
    !   grid point space  -> Fourier space
    !
    REAL(dp), TARGET, INTENT(inout) :: gp1(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: gp2(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: gp3(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: gp4(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: gp5(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: gp6(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: gp7(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: gp8(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: gp9(nproma, nlev, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: sf1(nproma, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: sf2(nproma, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: sf3(nproma, ngpblks)
    REAL(dp), TARGET, INTENT(inout) :: zm1(nlev, *)
    REAL(dp), TARGET, INTENT(inout) :: zm2(nlev, *)
    REAL(dp), TARGET, INTENT(inout) :: zm3(nlev, *)
    REAL(dp), TARGET, INTENT(in)    :: fs(nlon+2, nflevp1, nflat, *)
    REAL(dp), TARGET, INTENT(in)    :: fs0(nflev, nflat, *)

    TYPE(c_ptr) :: gp_c_ptr(9)
    TYPE(c_ptr) :: fs_c_ptr(9)
    TYPE(c_ptr) :: sf_c_ptr(3)
    TYPE(c_ptr) :: fs2d_c_ptr(3)
    TYPE(c_ptr) :: zm_c_ptr(3)
    TYPE(c_ptr) :: fs0_c_ptr(3)
    TYPE(c_ptr) :: all_src_c_ptr(15), all_dst_c_ptr(15)
    INTEGER :: m

    ! fs2gp_3d_redist:
    gp_c_ptr = (/ C_LOC(gp1), C_LOC(gp2), C_LOC(gp3), &
         &        C_LOC(gp4), C_LOC(gp5), C_LOC(gp6), &
         &        C_LOC(gp7), C_LOC(gp8), C_LOC(gp9)  /)
    DO m = 1, 9
      fs_c_ptr(m) = C_LOC(fs(1,1,1,m))
    ENDDO

    !fs2gp_2d_redist:
    sf_c_ptr = (/ C_LOC(sf1), C_LOC(sf2), C_LOC(sf3) /)
    DO m = 1, 3
      fs2d_c_ptr(m) = C_LOC(fs(1,nflevp1,1,m))
    ENDDO

    !fs2gp_zm_redist:
    zm_c_ptr = (/ C_LOC(zm1), C_LOC(zm2), C_LOC(zm3) /)
    DO m = 1, 3
      fs0_c_ptr(m) = C_LOC(fs0(1,1,m))
    ENDDO

    IF (aggregate_fs_to_gp) THEN
      DO m = 1,9
        all_src_c_ptr(m) =  fs_c_ptr(m)
        all_dst_c_ptr(m) =  gp_c_ptr(m)
      ENDDO
      DO m = 1,3
        all_src_c_ptr(9+m) =  fs2d_c_ptr(m)
        all_dst_c_ptr(9+m) =  sf_c_ptr(m)
      ENDDO
      DO m = 1,3
        all_src_c_ptr(12+m) =  fs0_c_ptr(m)
        all_dst_c_ptr(12+m) =  zm_c_ptr(m)
      ENDDO

      CALL xt_redist_s_exchange(fs2gp_all_redist, all_src_c_ptr, all_dst_c_ptr)

    ELSE
      DO m = 1, 9
        CALL xt_redist_s_exchange1(fs2gp_3d_redist, fs_c_ptr(m), gp_c_ptr(m))
      ENDDO
      DO m = 1, 3
        CALL xt_redist_s_exchange1(fs2gp_2d_redist, fs2d_c_ptr(m), sf_c_ptr(m) )
      ENDDO
      DO m = 1, 3
        CALL xt_redist_s_exchange1(fs2gp_zm_redist, fs0_c_ptr(m), zm_c_ptr(m))
      ENDDO

    ENDIF

    IF (npromz /= nproma) THEN
      CALL zero_fractional_gp_block_3d(gp1)
      CALL zero_fractional_gp_block_3d(gp2)
      CALL zero_fractional_gp_block_3d(gp3)
      CALL zero_fractional_gp_block_3d(gp4)
      CALL zero_fractional_gp_block_3d(gp5)
      CALL zero_fractional_gp_block_3d(gp6)
      CALL zero_fractional_gp_block_3d(gp7)
      CALL zero_fractional_gp_block_3d(gp8)
      CALL zero_fractional_gp_block_3d(gp9)
      CALL zero_fractional_gp_block_2d(sf1)
      CALL zero_fractional_gp_block_2d(sf2)
      CALL zero_fractional_gp_block_2d(sf3)
    ENDIF

  END SUBROUTINE cont_tr_fs_to_gp
#endif

  SUBROUTINE yaxt_tr_fs_to_gp(gp1, gp2, gp3, gp4, gp5, gp6, gp7, gp8, gp9, &
       &                      sf1, sf2, sf3, zm1, zm2, zm3, &
       &                      fs, fs0)
    !
    !   grid point space  -> Fourier space
    !
    REAL(dp), TARGET, INTENT(inout)  :: gp1    (:,:,:)   ! gridpoint space 3d
    REAL(dp), TARGET, INTENT(inout)  :: gp2    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp3    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp4    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp5    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp6    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp7    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp8    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp9    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: sf1    (:,:)     ! gridpoint space 2d
    REAL(dp), TARGET, INTENT(inout)  :: sf2    (:,:)     !
    REAL(dp), TARGET, INTENT(inout)  :: sf3    (:,:)     !
    REAL(dp), TARGET, INTENT(inout)  :: zm1    (:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(inout)  :: zm2    (:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(inout)  :: zm3    (:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(in)     :: fs     (:,:,:,:) ! Fourier space
    REAL(dp), TARGET, INTENT(in)     :: fs0    (:,:,:)   ! zonal mean, Four.
#ifdef HAVE_YAXT
    CHARACTER(len=*), PARAMETER :: subname = 'yaxt_tr_fs_to_gp'
#ifdef CHECK_SHAPES
    INTEGER :: nglat
    INTEGER :: gp3d_shape(3), gp2d_shape(2), zm_shape(2), fs3d_shape(3), fs0_2d_shape(2)
    INTEGER :: tmp_4d_shape(4), tmp_3d_shape(3), tmp_2d_shape(2)
    nglat = ldc%nglat
    gp3d_shape = (/ nproma, nlev, ngpblks /)
    gp2d_shape = (/ nproma, ngpblks /)
    zm_shape = (/ nlev, nglat /)
    fs3d_shape = (/ nlon+2, nflevp1, nflat /)
    fs0_2d_shape = (/ nflev, nflat /)
    IF (nproma * nlev * ngpblks * (nlon+2) * nflev * nflat == 0) &
         & CALL die(subname//": found unexpected zero shape product",__LINE__)
    IF (   ANY(SHAPE(gp1) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp2) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp3) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp4) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp5) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp6) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp7) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp8) /= gp3d_shape) .OR. &
         & ANY(SHAPE(gp9) /= gp3d_shape) ) &
         & CALL die(subname//": gp3d_shape mismatch",__LINE__)
    IF (   ANY(SHAPE(sf1) /= gp2d_shape) .OR. &
         & ANY(SHAPE(sf2) /= gp2d_shape) .OR. &
         & ANY(SHAPE(sf3) /= gp2d_shape) ) &
         & CALL die(subname//": gp2d_shape mismatch",__LINE__)
    IF (   ANY(SHAPE(zm1) /= zm_shape) .OR. &
         & ANY(SHAPE(zm2) /= zm_shape) .OR. &
         & ANY(SHAPE(zm3) /= zm_shape)    ) &
         & CALL die(subname//": zm_shape mismatch",__LINE__)
    tmp_4d_shape = SHAPE(fs)
    IF (   ANY(tmp_4d_shape(1:3) /= fs3d_shape) ) &
         & CALL die(subname//": fs3d_shape mismatch",__LINE__)
    tmp_3d_shape = SHAPE(fs0)
    IF (   ANY(tmp_3d_shape(1:2) /= fs0_2d_shape) ) &
         & CALL die(subname//": fs0_2d_shape mismatch",__LINE__)
#endif

    CALL cont_tr_fs_to_gp(gp1, gp2, gp3, gp4, gp5, gp6, gp7, gp8, gp9, &
         &                sf1, sf2, sf3, zm1, zm2, zm3, &
         &                fs, fs0)

#endif
  END SUBROUTINE yaxt_tr_fs_to_gp

  SUBROUTINE yaxt_tr_fs_to_ls(fs, ls, fs0, ls0)
    REAL(dp), TARGET, INTENT(in)  :: fs(:,:,:,:)
    REAL(dp), TARGET, INTENT(out) :: ls(:,:,:,:)
    REAL(dp), TARGET, INTENT(in)  :: fs0(:,:,:)
    REAL(dp), TARGET, INTENT(out) :: ls0(:,:,:)
#ifdef HAVE_YAXT
    INTEGER, PARAMETER :: n4_expected = 6
    INTEGER :: m, n4
    TYPE(c_ptr) :: fs_c_ptr(7)
    TYPE(c_ptr) :: ls_c_ptr(7)

#ifdef CHECK_MEM_STRIDES
    CALL check_mem_strides(fs, (/ 1, nlon+2, (nlon+2)*nflevp1, (nlon+2)*nflevp1*nflat /), 8, __LINE__)
    CALL check_mem_strides(fs0, (/ 1, nflev, nflev*nflat /), 8, __LINE__)
    CALL check_mem_strides(ls, (/ 1, 2*ldc%nlm, 2*ldc%nlm*nflevp1, 2*ldc%nlm*nflevp1*nlat /), 8, __LINE__)
    CALL check_mem_strides(ls0, (/ 1, nflev, nflev*nlat /), 8, __LINE__)
#endif

    n4 = MIN(SIZE(fs,4), SIZE(ls,4))
    IF (aggregate_fs_to_ls .AND. n4 == n4_expected) THEN
      DO m = 1, 6
        fs_c_ptr(m) =  C_LOC(fs(1,1,1,m))
        ls_c_ptr(m) =  C_LOC(ls(1,1,1,m))
      ENDDO
      m = 1
      fs_c_ptr(7) =  C_LOC(fs0(1,1,m))
      ls_c_ptr(7) =  C_LOC(ls0(1,1,m))
      CALL xt_redist_s_exchange(fs2ls_all_redist, fs_c_ptr, ls_c_ptr)
    ELSE
      DO m = 1, n4
        CALL xt_redist_s_exchange1(fs2ls_3d_redist, C_LOC(fs(1,1,1,m)), C_LOC(ls(1,1,1,m)) )
      ENDDO
      m=1
      CALL xt_redist_s_exchange1(fs2ls_zm_redist, C_LOC(fs0(1,1,m)), C_LOC(ls0(1,1,m)) )
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_fs_to_ls

  SUBROUTINE yaxt_tr_ls_to_fs(fs, ls, fs0, ls0)
    REAL(dp), TARGET, INTENT(out)  :: fs(:,:,:,:)
    REAL(dp), TARGET, INTENT(in)   :: ls(:,:,:,:)
    REAL(dp), TARGET, INTENT(out)  :: fs0(:,:,:)
    REAL(dp), TARGET, INTENT(in)   :: ls0(:,:,:)
#ifdef HAVE_YAXT
    TYPE(c_ptr) :: fs_c_ptr(12)
    TYPE(c_ptr) :: ls_c_ptr(12)
    INTEGER, PARAMETER :: n1_expected = 9
    INTEGER, PARAMETER :: n2_expected = 3
    INTEGER :: m, n1, n2

#ifdef CHECK_MEM_STRIDES
    CALL check_mem_strides(fs, (/ 1, nlon+2, (nlon+2)*nflevp1, (nlon+2)*nflevp1*nflat /), 8, __LINE__)
    CALL check_mem_strides(fs0, (/ 1, nflev, nflev*nflat /), 8, __LINE__)
    CALL check_mem_strides(ls, (/ 1, 2*ldc%nlm, 2*ldc%nlm*nflevp1, 2*ldc%nlm*nflevp1*nlat /), 8, __LINE__)
    CALL check_mem_strides(ls0, (/ 1, nflev, nflev*nlat /), 8, __LINE__)
#endif

    n1 = MIN(SIZE(fs,4), SIZE(ls,4))
    n2 = MIN(SIZE(fs0,3), SIZE(ls0,3))

    IF (aggregate_ls_to_fs .AND. n1==n1_expected .AND. n2==n2_expected) THEN
      DO m = 1, n1
        fs_c_ptr(m) =  C_LOC(fs(1,1,1,m))
        ls_c_ptr(m) =  C_LOC(ls(1,1,1,m))
      ENDDO
      DO m = 1, n2
        fs_c_ptr(n1+m) =  C_LOC(fs0(1,1,m))
        ls_c_ptr(n1+m) =  C_LOC(ls0(1,1,m))
      ENDDO
      CALL xt_redist_s_exchange(ls2fs_all_redist, ls_c_ptr, fs_c_ptr)
      DO m = 1, n1
        CALL zero_fs_3d(fs(:,:,:,m))
      ENDDO
    ELSE
      DO m = 1, n1
        CALL xt_redist_s_exchange1(ls2fs_3d_redist, C_LOC(ls(1,1,1,m)), C_LOC(fs(1,1,1,m)) )
        CALL zero_fs_3d(fs(:,:,:,m))
      ENDDO
      DO m = 1, n2
        CALL xt_redist_s_exchange1(ls2fs_zm_redist, C_LOC(ls0(1,1,m)), C_LOC(fs0(1,1,m)) )
      ENDDO
    ENDIF

  CONTAINS

    ! 3d fs post processing
    SUBROUTINE zero_fs_3d(f)
      REAL(dp), INTENT(inout) :: f(:,:,:)
      INTEGER :: i, jj, kk

      DO jj = 1, nflat
        DO kk = 1, nflevp1
          DO i = 2*(ldc%nm+1)+1, nlon+2
            f(i,kk,jj) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE zero_fs_3d
#endif
  END SUBROUTINE yaxt_tr_ls_to_fs

  SUBROUTINE yaxt_tr_ls_to_sp(ls1, sp1, ls2, sp2, ls3, sp3, ls0, sp0)
    REAL(dp), TARGET, INTENT(in)    :: ls1(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(inout) :: sp1(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(in)    :: ls2(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(inout) :: sp2(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(in)    :: ls3(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(inout) :: sp3(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(in)    :: ls0(:,:)   ! Legendre (m=0 only)
    REAL(dp), TARGET, INTENT(inout) :: sp0(:,:)   ! spectral (m=0 only)
#ifdef HAVE_YAXT
    TYPE(c_ptr) :: ls_c_ptr(4), sp_c_ptr(4), ls0_c_ptr, sp0_c_ptr
    ! use to prevent operations on null pointer
    REAL(dp), TARGET, SAVE :: dummy(1)

#ifdef CHECK_MEM_STRIDES
    INTEGER :: ls_nk(nlev:nlev+1)
    ls_nk(nlev)   = MIN(ldc%lleve,nlev)   - ldc%llevs + 1
    ls_nk(nlev+1) = MIN(ldc%lleve,nlev+1) - ldc%llevs + 1
    CALL check_mem_strides(ls1, (/ 1, ls_nk(nlev), 2*ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(ls3, (/ 1, ls_nk(nlev), 2*ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(ls2, (/ 1, ls_nk(nlev+1), 2*ls_nk(nlev+1) /), 8, __LINE__)
    CALL check_mem_strides(sp1, (/ 1, nlev, 2*nlev /), 8, __LINE__)
    CALL check_mem_strides(sp3, (/ 1, nlev, 2*nlev /), 8, __LINE__)
    CALL check_mem_strides(sp2, (/ 1, nlev+1, 2*(nlev+1)/), 8, __LINE__)
    CALL check_mem_strides(ls0, (/ 1, ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(sp0, (/ 1, nlev /), 8, __LINE__)
#endif

    !todo: use call to assumed size/explicit shape variant to save check
    IF (SIZE(ls0)>0) THEN
      ls0_c_ptr = C_LOC(ls0(1,1))
    ELSE
      ls0_c_ptr = C_LOC(dummy)
    ENDIF
    IF (SIZE(sp0)>0) THEN
      sp0_c_ptr = C_LOC(sp0(1,1))
    ELSE
      sp0_c_ptr = C_LOC(dummy)
    ENDIF

    IF (aggregate_ls_to_sp) THEN
      ls_c_ptr(1) = C_LOC(ls1(1,1,1))
      sp_c_ptr(1) = C_LOC(sp1(1,1,1))

      ls_c_ptr(2) = C_LOC(ls3(1,1,1))
      sp_c_ptr(2) = C_LOC(sp3(1,1,1))

      ls_c_ptr(3) = C_LOC(ls2(1,1,1))
      sp_c_ptr(3) = C_LOC(sp2(1,1,1))

      ls_c_ptr(4) = ls0_c_ptr
      sp_c_ptr(4) = sp0_c_ptr
      CALL xt_redist_s_exchange(ls2sp_all_redist, ls_c_ptr, sp_c_ptr)
    ELSE
      CALL xt_redist_s_exchange1(ls2sp_3d_redist(nlev), C_LOC(ls1(1,1,1)), C_LOC(sp1(1,1,1)) )
      CALL xt_redist_s_exchange1(ls2sp_3d_redist(nlev), C_LOC(ls3(1,1,1)), C_LOC(sp3(1,1,1)) )
      CALL xt_redist_s_exchange1(ls2sp_3d_redist(nlev+1), C_LOC(ls2(1,1,1)), C_LOC(sp2(1,1,1)) )
      CALL xt_redist_s_exchange1(ls2sp_m0_redist, ls0_c_ptr, sp0_c_ptr)
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_ls_to_sp

  SUBROUTINE yaxt_tr_sp_to_ls(ls1, sp1, ls2, sp2, ls3, sp3, ls0, sp0)
    REAL(dp), TARGET, INTENT(inout) :: ls1(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(in)    :: sp1(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(inout) :: ls2(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(in)    :: sp2(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(inout) :: ls3(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(in)    :: sp3(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(inout) :: ls0(:,:)   ! Legendre (m=0 only)
    REAL(dp), TARGET, INTENT(in)    :: sp0(:,:)   ! spectral (m=0 only)
#ifdef HAVE_YAXT
    TYPE(c_ptr) :: ls_c_ptr(4), sp_c_ptr(4), ls0_c_ptr, sp0_c_ptr
    ! use to prevent operations on null pointer
    REAL(dp), TARGET, SAVE :: dummy(1)

#ifdef CHECK_MEM_STRIDES
    INTEGER :: ls_nk(nlev:nlev+1)
    ls_nk(nlev)   = MIN(ldc%lleve,nlev)   - ldc%llevs + 1
    ls_nk(nlev+1) = MIN(ldc%lleve,nlev+1) - ldc%llevs + 1
    CALL check_mem_strides(ls1, (/ 1, ls_nk(nlev), 2*ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(ls3, (/ 1, ls_nk(nlev), 2*ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(ls2, (/ 1, ls_nk(nlev+1), 2*ls_nk(nlev+1) /), 8, __LINE__)
    CALL check_mem_strides(sp1, (/ 1, nlev, 2*nlev /), 8, __LINE__)
    CALL check_mem_strides(sp3, (/ 1, nlev, 2*nlev /), 8, __LINE__)
    CALL check_mem_strides(sp2, (/ 1, nlev+1, 2*(nlev+1)/), 8, __LINE__)
    CALL check_mem_strides(ls0, (/ 1, ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(sp0, (/ 1, nlev /), 8, __LINE__)
#endif

    IF (SIZE(ls0)>0) THEN
      ls0_c_ptr = C_LOC(ls0(1,1))
    ELSE
      ls0_c_ptr = C_LOC(dummy)
    ENDIF
    IF (SIZE(sp0)>0) THEN
      sp0_c_ptr = C_LOC(sp0(1,1))
    ELSE
      sp0_c_ptr = C_LOC(dummy)
    ENDIF

    IF (aggregate_sp_to_ls) THEN
      ls_c_ptr(1) = C_LOC(ls1(1,1,1))
      sp_c_ptr(1) = C_LOC(sp1(1,1,1))

      ls_c_ptr(2) = C_LOC(ls3(1,1,1))
      sp_c_ptr(2) = C_LOC(sp3(1,1,1))

      ls_c_ptr(3) = C_LOC(ls2(1,1,1))
      sp_c_ptr(3) = C_LOC(sp2(1,1,1))

      ls_c_ptr(4) = ls0_c_ptr
      sp_c_ptr(4) = sp0_c_ptr
      CALL xt_redist_s_exchange(sp2ls_all_redist, sp_c_ptr, ls_c_ptr)
    ELSE
      CALL xt_redist_s_exchange1(sp2ls_3d_redist(nlev), C_LOC(sp1(1,1,1)), C_LOC(ls1(1,1,1)) )
      CALL xt_redist_s_exchange1(sp2ls_3d_redist(nlev), C_LOC(sp3(1,1,1)), C_LOC(ls3(1,1,1)) )
      CALL xt_redist_s_exchange1(sp2ls_3d_redist(nlev+1), C_LOC(sp2(1,1,1)), C_LOC(ls2(1,1,1)) )
      CALL xt_redist_s_exchange1(sp2ls_m0_redist, sp0_c_ptr, ls0_c_ptr)
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_sp_to_ls

#ifdef HAVE_YAXT
  !
  ! start of yaxt-only section until end of module
  !
  SUBROUTINE prepare_new_gp_coords
    INTEGER :: ia, ib, ir, ig, jg, i, j
#ifdef DYNAMIC_DECOMPOSITION
    ia = nproma
    ib = 0
    DO ir = 1, 2
      DO jg = glats(ir), glate(ir)
        j = ilat_map(jg)
        DO ig = glons(ir), glone(ir)
          i = ilon_map(ig)
          ia = ia + 1
          IF (ia > nproma) THEN
            ib = ib + 1
            ia = ia - nproma
          ENDIF
          ldc%ilon(ia,ib) = i
          ldc%ilat(ia,ib) = j
        ENDDO
      ENDDO
    ENDDO
#else
    IF (nonid_ilon_map .OR. nonid_ilat_map) THEN
      CALL die("Define DYNAMIC_DECOMPOSITION in order to have ldc%ilon, ldc%ilat.",__LINE__)
    ENDIF
#endif

  END SUBROUTINE prepare_new_gp_coords

  SUBROUTINE gen_index_maps
    LOGICAL :: is_id
    INTEGER :: nhl, j

    ALLOCATE(ilon_map(nlon), ilat_map(nlat))
    IF (nlonstrip <= 0) THEN
      CALL gen_id_map(ilon_map)
      nonid_ilon_map = .FALSE.
    ELSE
      CALL gen_round_robin(ilon_map, nprocb, nlonstrip, is_id)
      nonid_ilon_map = (.NOT. is_id)
    ENDIF

    IF (nlatstrip <= 0) THEN
      CALL gen_id_map(ilat_map)
    ELSE
      IF (preserve_latsym) THEN
        nhl = nlat/2
        IF (nhl*2  /= nlat) CALL die('preserve_latsym requires even nlat',__LINE__)
        CALL gen_round_robin(ilat_map(1:nhl), nproca, nlatstrip, is_id)
        DO j = 1, nhl
          ilat_map(nhl+j) = nlon - ilat_map(nhl-j+1) + 1
        ENDDO
      ELSE
        CALL gen_round_robin(ilat_map, 2*nproca, nlatstrip, is_id)
        nonid_ilat_map = (.NOT. is_id)
      ENDIF
    ENDIF

  END SUBROUTINE gen_index_maps

  SUBROUTINE gen_id_map(map)
    INTEGER, INTENT(out) :: map(:)
    INTEGER :: i
    DO i = 1, SIZE(map)
      map(i) = i
    ENDDO
  END SUBROUTINE gen_id_map

  SUBROUTINE gen_round_robin(map,ntiles,nstrip, is_id)
    INTEGER, INTENT(out) :: map(:)
    INTEGER, INTENT(in) :: ntiles ! number of tiles
    INTEGER, INTENT(in) :: nstrip ! width of stripe
    LOGICAL, intent(out) :: is_id

    INTEGER :: h, npt, i, j, p
    INTEGER :: ts(ntiles), te(ntiles), tn(ntiles)
    INTEGER :: simple_npt
    INTEGER :: simple_map(SIZE(map))

    npt = SIZE(map)
    ! h = modified nstrip, such that h devides npt:
    ! h is groupsize of points that stay together
    h = MAX(1,MIN(npt,nstrip))
    DO WHILE (MOD(npt,h) /= 0)
      h = h-1
    ENDDO
    simple_npt = npt/h
    CALL  gen_simple_round_robin(simple_map(1:simple_npt),ntiles)
    is_id = .TRUE.
    DO i = 1, simple_npt
      DO j=1,h
        p = (i-1)*h+j
        map(p)=h*(simple_map(i)-1)+j
        IF (map(p) /= p) is_id = .FALSE.
      ENDDO
    ENDDO
    IF (ANY(map == 0) .OR. 2*SUM(map) /= npt*(npt+1)) CALL die('gen_round_robin: internal error',__LINE__)
  END SUBROUTINE gen_round_robin

  SUBROUTINE gen_simple_round_robin(map,ntiles)
    INTEGER, INTENT(out) :: map(:)
    INTEGER, INTENT(in) :: ntiles ! number of tiles

    CHARACTER(len=*), PARAMETER :: context = 'gen_idx_maps: simple_round_robin: '
    INTEGER :: i,m,n,npt,p, pstep, idx, it, j
    INTEGER :: ts(ntiles), te(ntiles), tn(ntiles)

    npt = SIZE(map)

    ! tiles:
    IF (ntiles == npt) THEN
      CALL gen_id_map(map)
      RETURN
    ELSEIF (ntiles > npt) THEN
      CALL die(context//'( ntiles > npt )',__LINE__)
    ENDIF

    n = npt/ntiles
    m = npt-n*ntiles
    DO i = 1, m
      tn(i) = n+1
    ENDDO
    DO i = m+1, ntiles
      tn(i) = n
    ENDDO

    ts(1) = 1
    te(1) = tn(1)
    DO i = 2, ntiles
      ts(i) = te(i-1) + 1
      te(i) = ts(i) + tn(i) -1
    ENDDO
    IF (te(ntiles) /= npt) CALL die(context//'internal error',__LINE__)

    ! round roubin:
    map(:) = 0
    pstep = 0
    idx = 1
    DO WHILE (idx <= npt)
      DO it = 1, ntiles

        p = ts(it) + pstep
        IF (p <= te(it) .AND. map(p) == 0) THEN
          map(p) = idx
          idx = idx+1
        ENDIF

      ENDDO
      pstep = pstep+1
    ENDDO
    IF (ANY(map == 0) .or. 2*SUM(map) /= npt*(npt+1)) CALL die(context//'internal error',__LINE__)

  END SUBROUTINE gen_simple_round_robin

  SUBROUTINE zero_fractional_gp_block_3d(g)
    REAL(dp), INTENT(inout) :: g(:,:,:)
    INTEGER :: i, k
    DO k=1, nlev
      DO i = npromz+1, nproma
        g(i,k,ngpblks) = 0.0_dp
      ENDDO
    ENDDO
  END SUBROUTINE zero_fractional_gp_block_3d

  SUBROUTINE zero_fractional_gp_block_2d(g)
    REAL(dp), INTENT(inout) :: g(:,:)
    INTEGER :: i
    DO i = npromz+1, nproma
      g(i,ngpblks) = 0.0_dp
    ENDDO
  END SUBROUTINE zero_fractional_gp_block_2d

  INTEGER FUNCTION get_active_trnum() RESULT(active_trnum)
#ifndef MESSY
    active_trnum == COUNT(tracdef_trlist%ti(1:tracdef_trlist%ntrac)%ntran /= 0)
#else
    INTEGER :: itracer
    active_trnum = 0
    DO itracer=1, ntrac_gp
      IF (ti_gp(itracer)%tp%meta%cask_i(I_advect) == on_switch) active_trnum = active_trnum + 1
    END DO
#endif
  END FUNCTION get_active_trnum

  INTEGER FUNCTION get_trdim_size() RESULT(trdim_size)
#ifndef MESSY
    trdim_size = tracdef_ntrac
#else
    USE mo_memory_g1a, ONLY: xtm1
    trdim_size = SIZE(xtm1,3)
#endif
  END FUNCTION get_trdim_size

  LOGICAL FUNCTION match_active_trpos() RESULT(m)
    INTEGER :: iit, it
    LOGICAL :: is_active
    m = ( active_trnum == get_active_trnum() )
    IF (.NOT. m) RETURN


    DO iit = 1, active_trnum
      it = active_trpos(iit)
#ifndef MESSY
      is_active = (tracdef_trlist%ti(it)%ntran /= 0)
#else
      is_active = (ti_gp(it)%tp%meta%cask_i(I_advect) == on_switch)
#endif
      IF (.NOT. is_active) THEN
        m = .FALSE.
        RETURN
      ENDIF
    ENDDO

  END FUNCTION match_active_trpos

  LOGICAL FUNCTION next_active_tracer(current_it) RESULT(found)
    INTEGER, INTENT(inout) :: current_it
    INTEGER :: it

    found = .false.
#ifndef MESSY
    DO it = current_it+1, tracdef_trlist%ntrac
      IF (tracdef_trlist%ti(it)%ntran /= 0) THEN
        current_it = it
        found = .TRUE.
        RETURN
      ENDIF
    ENDDO
#else
    DO it=current_it+1, ntrac_gp
      IF (ti_gp(it)%tp%meta%cask_i(I_advect) == on_switch) THEN
        current_it = it
        found = .TRUE.
        RETURN
      ENDIF
    END DO
#endif
  END FUNCTION next_active_tracer

#if 0
  ! zero selected 4d fs data
  SUBROUTINE zero_fs_4d(f,m1,m2)
    REAL(dp), INTENT(inout) :: f(:,:,:,:)
    INTEGER, INTENT(in) :: m1, m2
    INTEGER :: i, jj, kk, m, f_shape(4)

    f_shape = SHAPE(f)
    DO m = m1, m2
      DO jj = 1, nflat
        DO kk = 1, nflevp1
          DO i = nlon+1, nlon+2
            f(i,kk,jj,m) = 0.0_dp
          ENDDO
        ENDDO
        DO kk = nflev+1,nflevp1
          DO i = 1, nlon+2
            f(i,kk,jj,m) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE zero_fs_4d
#endif

  ! zero tail of fractional nproma block for 4d arrays
  SUBROUTINE zero_fractional_gp_block_4d(g)
    REAL(dp), INTENT(inout) :: g(:,:,:,:)
    INTEGER :: i, k, it, nt
    nt = SIZE(g,3)
    DO it = 1, nt
      DO k=1, nlev
        DO i = npromz+1, nproma
          g(i,k,it,ngpblks) = 0.0_dp
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE zero_fractional_gp_block_4d

  SUBROUTINE setup_empty_deco
    empty_deco%idxlist =  xt_idxempty_new()
    ALLOCATE(empty_deco%offset(1))
    empty_deco%offset(1) = 0
  END SUBROUTINE setup_empty_deco

  SUBROUTINE setup_gp_gdeco(idx_description_type)
    INTEGER, INTENT(IN) :: idx_description_type
    INTEGER :: nglon, nglat

    nglon = ldc%nglon
    nglat = ldc%nglat

    CALL setup_gp

  CONTAINS

    !
    ! Define global index list and offsets for Gaussian grid
    !
    SUBROUTINE setup_gp
      INTEGER :: nidx, nstr, nk, i
      INTEGER :: trpos(pcnst)

      ALLOCATE(gp_gdeco(SIZE(gp_nlevs)))
      ! the next line is only needed for some buggy compilers
      DO nk = 1, SIZE(gp_gdeco)
        gp_gdeco(nk)%idxlist = xt_idxlist_c2f(c_null_ptr)
      ENDDO

      DO nk = 1, SIZE(gp_nlevs)
        ! generate index list
        IF (gp_nlevs(nk)) THEN
          IF (.NOT. xt_is_null(gp_gdeco(nk)%idxlist)) &
               & CALL finish('mo_echam_yaxt: setup_gp', 'Internal error.')
          nidx = nglon * nk * nglat
          SELECT CASE (idx_description_type)
            ! index list as vector
            CASE (idx_descr_vector)
              gp_gdeco(nk)%idxlist = new_gp_idxvec(nidx, nk)
            ! index list by stripes
            CASE (idx_descr_stripes)
              nstr = nk * nglat
              gp_gdeco(nk)%idxlist = new_gp_idxstripes(nstr, nk)
            CASE DEFAULT
              WRITE (message_text, '(a,i0,a)') &
              'Unsupported global index desciption type: ', &
              idx_description_type, ' (Gaussain grid)'
              CALL finish('setup_yaxt_decomposition', message_text)
          END SELECT

          ! generate offsets
          ALLOCATE(gp_gdeco(nk)%offset(nidx))
          CALL set_gp_offset(nidx, nk, gp_gdeco(nk)%offset)
        ENDIF
      ENDDO

      gp_3d_fb_single_tracer_gdeco%idxlist =  gp_gdeco(nlev)%idxlist
      nidx = nglon * nlev * nglat
      ALLOCATE(gp_3d_fb_single_tracer_gdeco%offset(nidx))
      CALL set_gp_3d_single_tracer_offset(nidx, nlev, pcnst, gp_3d_fb_single_tracer_gdeco%offset)

      gp_3d_xtm1_single_tracer_gdeco%idxlist =  gp_gdeco(nlev)%idxlist
      nidx = nglon * nlev * nglat
      ALLOCATE(gp_3d_xtm1_single_tracer_gdeco%offset(nidx))
      CALL set_gp_3d_single_tracer_offset(nidx, nlev, trdim_size, gp_3d_xtm1_single_tracer_gdeco%offset)

      IF (is_gather_target) THEN
        ! gather target
        nidx = nlon * nlat
        global_gp_2d_gdeco%idxlist = new_global_gp_3d_idxvec(nidx, 1)
        ALLOCATE(global_gp_2d_gdeco%offset(nidx))
        CALL set_global_gp_3d_offset(nidx, 1, global_gp_2d_gdeco%offset)
      ELSE
        global_gp_2d_gdeco%idxlist = xt_idxempty_new()
        ALLOCATE(global_gp_2d_gdeco%offset(1))
        global_gp_2d_gdeco%offset = 0
      ENDIF

      ! gather 3d:
      IF (is_gather_target) THEN
        ! gather target
        nidx = nlon * nlev * nlat
        global_gp_3d_gdeco%idxlist = new_global_gp_3d_idxvec(nidx, nlev)
        ALLOCATE(global_gp_3d_gdeco%offset(nidx))
        CALL set_global_gp_3d_offset(nidx, nlev, global_gp_3d_gdeco%offset)
      ELSE
        global_gp_3d_gdeco%idxlist = xt_idxempty_new()
        ALLOCATE(global_gp_3d_gdeco%offset(1))
        global_gp_3d_gdeco%offset = 0
      ENDIF

    END SUBROUTINE setup_gp

    TYPE(xt_idxlist) FUNCTION new_global_gp_3d_idxvec(nidx, kmax)
      INTEGER, INTENT(in) :: nidx, kmax
      INTEGER :: idx(nidx)
      INTEGER :: ig, jg
      INTEGER :: i, j, k, p

      p = 0
      DO jg = 1, nlat
        j = jg
        DO k = 1, kmax
          DO ig = 1, nlon
            i = ig
            p = p + 1
            idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
          ENDDO !i
        ENDDO !j
      ENDDO !k
      IF (p /= nidx) CALL die("new_global_gp_3d_idxvec: internal error",__LINE__)
      new_global_gp_3d_idxvec = xt_idxvec_new(idx, nidx)

    END FUNCTION new_global_gp_3d_idxvec

    SUBROUTINE set_global_gp_3d_offset(nidx, kmax, offset)
      INTEGER, INTENT(IN)  :: nidx, kmax
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: i, j, k, p, r, ib, ia

      p = 0
      DO j = 1, nlat
        DO k = 1, kmax
          DO i = 1, nlon
            p = p + 1
            offset(p) = p-1
          ENDDO !i
        ENDDO !k
      ENDDO !j

    END SUBROUTINE set_global_gp_3d_offset

    !
    !
    ! Gaussian grid
    ! use index vector to describe decomposition
    TYPE(xt_idxlist) FUNCTION new_gp_idxvec(nidx, nlv)
      INTEGER, INTENT(IN) :: nidx, nlv
      INTEGER :: idx(nidx)
      INTEGER :: i, ig, j, jg, k, p, r
      p = 0
      DO k = 1, nlv
        DO r = 1, 2
          DO jg = glats(r), glate(r)
            j = ilat_map(jg)
            DO ig = glons(r), glone(r)
              i = ilon_map(ig)
              p = p + 1
              idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
            ENDDO !i
          ENDDO !j
        ENDDO !r
      ENDDO !k
      new_gp_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_gp_idxvec

    !
    ! Gaussian grid
    ! use index stripes to describe decomposition
    TYPE(xt_idxlist) FUNCTION new_gp_idxstripes(nstr, nlv)
      INTEGER, INTENT(IN) :: nstr, nlv
      TYPE(xt_stripe) :: s(nstr)
      INTEGER :: p, k, r, j
      IF (nonid_ilon_map .OR. nonid_ilat_map) &
           &       CALL die("new_gp_idxstripes(nstr, nlv): (nonid_ilon_map .OR. nonid_ilat_map) ",__LINE__)

      s(:)%stride = 1
      p = 0 ! strides index
      DO k = 1, nlv
        DO r = 1, 2
          DO j = glats(r), glate(r)
            p = p + 1
            s(p)%start = glons(r) + nlon * ( (j-1) + nlat * (k-1) ) - 1
            s(p)%nstrides = nglon
          ENDDO !j
        ENDDO !r
      ENDDO !k
      new_gp_idxstripes = xt_idxstripes_new(s, nstr)
    END FUNCTION new_gp_idxstripes

    !
    ! Gaussian grid
    ! set offsets
    SUBROUTINE set_gp_offset(nidx, nlv, offset)
      INTEGER, INTENT(IN)  :: nidx, nlv
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: i, j, k, p, r, ib, ia

      p = 0
      DO k = 1, nlv
        ib = 1
        ia = 0
        DO r = 1, 2
          DO j = glats(r), glate(r)
            DO i = glons(r), glone(r)
              p = p + 1
              ia = ia +1
              IF (ia > nproma) THEN
                ib = ib + 1
                ia = 1
              ENDIF
              offset(p) = ia + nproma * (k - 1) + nproma * nlv * (ib - 1) - 1
            ENDDO !i
          ENDDO !j
        ENDDO !r
      ENDDO !k
    END SUBROUTINE set_gp_offset

    !
    ! Gaussian grid with single tracer
    ! -- just for testing -- will be replaced by std-3d-gp-deco
    TYPE(xt_idxlist) FUNCTION new_gp_3d_single_tracer_idxvec(nidx, nlv, nt)
      INTEGER, INTENT(IN) :: nidx, nlv, nt
      INTEGER :: idx(nidx)
      INTEGER :: i, ig, j, jg, k, it, p, r
      p = 0
      DO k = 1, nlv
        DO r = 1, 2
          DO jg = glats(r), glate(r)
            j = ilat_map(jg)
            DO ig = glons(r), glone(r)
              i = ilon_map(ig)
              p = p + 1
              idx(p) = (i-1) + nlon * ( (j-1) + nlat * (k-1) )
            ENDDO !i
          ENDDO !j
        ENDDO !r
      ENDDO !k
      IF (p /= nidx) CALL die('bad case', __LINE__)
      new_gp_3d_single_tracer_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_gp_3d_single_tracer_idxvec

    !
    ! Gaussian grid with tracer dimension
    ! set offsets
    SUBROUTINE set_gp_3d_single_tracer_offset(nidx, nlv, trsize, offset)
      INTEGER, INTENT(IN)  :: nidx, nlv, trsize
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: i, j, k, it, p, r, ib, ia
      INTEGER :: w_ia, w_ib, w_k, w_it

      ! offsets weights for different loops:
      w_ia = 1
      w_k  = w_ia * nproma
      w_it = w_k * nlv
      w_ib = w_it * trsize

      it = 1
      p = 0
      DO k = 1, nlv
        ib = 1
        ia = 0
        DO r = 1, 2
          DO j = glats(r), glate(r)
            DO i = glons(r), glone(r)
              p = p + 1
              ia = ia +1
              IF (ia > nproma) THEN
                ib = ib + 1
                ia = 1
              ENDIF
              IF (p > nidx) CALL die("set_gp_3d_single_tracer_offset: bad case",__LINE__)
              offset(p) = (ia-1)*w_ia + (ib-1)*w_ib + (k-1)*w_k !always zero: + (it-1)*w_it
            ENDDO !i
          ENDDO !j
        ENDDO !r
      ENDDO !k
      IF (p /= nidx) CALL die("set_gp_3d_single_tracer_offset: bad case",__LINE__)
    END SUBROUTINE set_gp_3d_single_tracer_offset

  END SUBROUTINE setup_gp_gdeco

  ! decomposition & offsets of horizontal means
  SUBROUTINE setup_gp_hmdeco
    INTEGER :: num_indices

    num_indices = nlev
    gp_hmdeco%idxlist = new_gp_hm_idxvec(num_indices)

    ALLOCATE(gp_hmdeco%offset(num_indices))
    CALL set_gp_hm_offset(num_indices, gp_hmdeco%offset)

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_gp_hm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: k
      DO k = 1, nlev
        idx(k) = k - 1
      ENDDO
      new_gp_hm_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_gp_hm_idxvec

    SUBROUTINE set_gp_hm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: k
      DO k = 1, nlev
        offset(k) = k-1
      ENDDO
    END SUBROUTINE set_gp_hm_offset

  END SUBROUTINE setup_gp_hmdeco

  ! decomposition & offsets of zonal means
  SUBROUTINE setup_gp_zmdeco
    INTEGER :: nglat
    nglat = ldc%nglat

    CALL setup_gp_zm

  CONTAINS

    SUBROUTINE setup_gp_zm
      INTEGER :: nidx
      ! zonal means for all levels:
      nidx = nlev*nglat

      gp_zmdeco%idxlist = new_gp_zm_idxvec(nidx)
      ALLOCATE(gp_zmdeco%offset(nidx))
      CALL set_gp_zm_offset(nidx, gp_zmdeco%offset)

      ! zonal means for fs levels:
      nidx = nflev*nglat
      gp_sel_zmdeco%idxlist = new_gp_sel_zm_idxvec(nidx)
      ALLOCATE(gp_sel_zmdeco%offset(nidx))
      CALL set_gp_sel_zm_offset(nidx, gp_sel_zmdeco%offset)

    END SUBROUTINE setup_gp_zm

    ! indices for grid point zonal means for all levels
    TYPE(xt_idxlist) FUNCTION new_gp_zm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: k, j, jg, r, p

      p = 0
      DO r = 1, 2
        DO jg = glats(r), glate(r)
          j = jg
          DO k = 1, nlev
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = k + nlev * (j-1)  - 1
          ENDDO
        ENDDO
      ENDDO

      new_gp_zm_idxvec = xt_idxvec_new(idx, nidx)

    END FUNCTION new_gp_zm_idxvec

    ! offsets for grid point zonal means for all levels
    SUBROUTINE set_gp_zm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: jj, k, p

      ! data shape: (nlev, nglat)
      ! see: variables ul, u0, du0 in module mo_scan_buffer
      p = 0
      DO jj = 1, nglat
        DO k = 1, nlev
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          offset(p) = k - 1 + nlev*(jj-1)
        ENDDO
      ENDDO
    END SUBROUTINE set_gp_zm_offset

    ! indices for grid point zonal means within fs levels
    TYPE(xt_idxlist) FUNCTION new_gp_sel_zm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: k, j, jg, r, p, confined_fleve

      ! indices belong to zonal index space
      confined_fleve = MIN(fleve, nlev)

      p = 0
      DO r = 1, 2
        DO jg = glats(r), glate(r)
          j = jg!ilat_map(jg)
          DO k = flevs, confined_fleve
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = k + nlev * (j-1)  - 1
          ENDDO
        ENDDO
      ENDDO

      new_gp_sel_zm_idxvec = xt_idxvec_new(idx, nidx)

    END FUNCTION new_gp_sel_zm_idxvec

    ! offsets for grid point zonal means within fs levels
    SUBROUTINE set_gp_sel_zm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: jj, k, p, confined_fleve

      ! data shape: (nlev, nglat)
      ! see: variables ul, u0, du0 in module mo_scan_buffer
      confined_fleve = MIN(fleve, nlev)
      p = 0
      DO jj = 1, nglat
        DO k = flevs, confined_fleve
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          offset(p) = k - 1 + nlev*(jj-1)
        ENDDO
      ENDDO
    END SUBROUTINE set_gp_sel_zm_offset

  END SUBROUTINE setup_gp_zmdeco


  SUBROUTINE setup_fs_2d_gdeco
    INTEGER :: num_indices

    ! the surface case is only used within the extra fs layer:
    IF (nflevp1 > nflev) THEN
      num_indices = nflat*nlon
      fs_2d_gdeco%idxlist = new_fs_2d_idxvec(num_indices)
      ALLOCATE(fs_2d_gdeco%offset(num_indices))
      CALL set_fs_2d_offset(num_indices, fs_2d_gdeco%offset)
    ELSE
      ! this pe has no extra layer
      fs_2d_gdeco%idxlist = xt_idxempty_new()
      ! we still need a valid address for offsets
      ALLOCATE(fs_2d_gdeco%offset(1))
      fs_2d_gdeco%offset(1) = 0
    ENDIF

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_fs_2d_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)

      INTEGER :: i, j, r, p

      p = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          DO i = 1, nlon
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = i + nlon * (j-1)  - 1
          ENDDO
        ENDDO
      ENDDO

      new_fs_2d_idxvec = xt_idxvec_new( idx, nidx)
    END FUNCTION new_fs_2d_idxvec

    SUBROUTINE set_fs_2d_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: r, i, j, jj, p

      ! fs 2d-data: (1:nlon+2, nflevp1, 1:nflat)

      p = 0
      jj = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          jj = jj + 1
          DO i = 1, nlon
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            offset(p) = i + (nlon+2) * nflevp1 * (jj-1) - 1
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE set_fs_2d_offset

  END SUBROUTINE setup_fs_2d_gdeco

  SUBROUTINE setup_fs_3d_gdeco

    CALL setup_fs_3d

  CONTAINS

    SUBROUTINE setup_fs_3d
      INTEGER :: nidx, confined_fleve, global_kmax, nk(nlev:nlev+1)

      ! std fs deco:
      nk = (/ nflev, nflevp1 /)
      ALLOCATE(fs_gdeco(nlev:nlev+1))
      DO global_kmax = nlev, nlev+1
        confined_fleve = MIN(fleve, global_kmax)
        nidx = nk(global_kmax)*nflat*nlon
        fs_gdeco(global_kmax)%idxlist = new_fs_3d_idxvec(nidx, nlon, confined_fleve)
        ALLOCATE(fs_gdeco(global_kmax)%offset(nidx))
        CALL set_fs_3d_offset(nidx, fs_gdeco(global_kmax)%offset, nlon, confined_fleve)
      ENDDO

      ! special subset deco for transpositions between fs and ls:
      nidx = nflevp1*nflat*2*(ldc%nm+1)
      fs_subset_gdeco%idxlist = new_fs_3d_idxvec(nidx, 2*(ldc%nm+1), fleve)
      ALLOCATE(fs_subset_gdeco%offset(nidx))
      CALL set_fs_3d_offset(nidx, fs_subset_gdeco%offset, 2*(ldc%nm+1), fleve)

    END SUBROUTINE setup_fs_3d

    ! fs 3d indices
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_fs_3d_idxvec(nidx,active_nlon,active_fleve)
      INTEGER, INTENT(IN) :: nidx, active_nlon, active_fleve
      INTEGER :: idx(nidx)

      INTEGER :: i, j, k, r, p

      ! assumptions:
      !
      ! in fs-space echam distributes nlev+1 levels in block form,
      ! these blocks are described by (flevs:fleve), we have nflevp1 = fleve-flevs+1
      ! and nflev= min(fleve,nlev)-flevs, empty blocks have nlevp1 <= 0,
      ! the transpositions from/to gp-deco only touch levels <= nlev,
      ! the level nlev+1 carries extra data (outside the normal grid point space)

      IF (fleve-flevs+1 /= nflevp1) CALL die('bad case',__LINE__)

      p = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          DO k = flevs, active_fleve
            DO i = 1, active_nlon
              p = p + 1
              IF (p > nidx) CALL die('bad case', __LINE__)
              idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      new_fs_3d_idxvec = xt_idxvec_new( idx, nidx)
    END FUNCTION new_fs_3d_idxvec

    SUBROUTINE set_fs_3d_offset(nidx, offset, active_nlon, active_fleve)
      INTEGER, INTENT(IN) :: nidx, active_nlon, active_fleve
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: r, i, j, k, jc, kc, p

      ! assumptions:
      !
      ! fs data layout:
      ! shape(fftz) =  (dc% nlon+2, dc% nflevp1, dc% nflat, nvar)

      p = 0
      jc = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          jc = jc + 1
          DO k = flevs, active_fleve
            kc = k - flevs + 1
            DO i = 1, active_nlon
              p = p + 1
              IF (p > nidx) CALL die('bad case', __LINE__)
              offset(p) = i-1 + (nlon+2) * ( (kc-1)  + nflevp1 * (jc-1))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE set_fs_3d_offset

  END SUBROUTINE setup_fs_3d_gdeco

  ! fs horizontal means:
  SUBROUTINE setup_fs_hmdeco
    INTEGER :: num_indices, confined_fleve

    confined_fleve = MIN(fleve, nlev)
    num_indices = nflev
    fs_hmdeco%idxlist = new_fs_hm_idxvec(num_indices)
    ALLOCATE(fs_hmdeco%offset(num_indices))
    CALL set_fs_hm_offset(num_indices, fs_hmdeco%offset)

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_fs_hm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: p,r,k
      p = 0
      DO k = flevs, confined_fleve
        p = p + 1
        IF (p > nidx) CALL die('bad case', __LINE__)
        idx(p) = k - 1
      ENDDO
      IF (p /= nidx) CALL die('bad case', __LINE__)
      new_fs_hm_idxvec = xt_idxvec_new(idx, nidx)

    END FUNCTION new_fs_hm_idxvec

    SUBROUTINE set_fs_hm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: k, p

      ! offsets for horizontal means

      ! data shape: (nflev)
      p = 0
      DO k = 1, nflev
        p = p + 1
        IF (p > nidx) CALL die('bad case', __LINE__)
        offset(p) = p -1
      ENDDO
    END SUBROUTINE set_fs_hm_offset


  END SUBROUTINE setup_fs_hmdeco

  ! fs zonal means:
  SUBROUTINE setup_fs_zmdeco
    INTEGER :: num_indices

    num_indices = nflev*nflat
    fs_zmdeco%idxlist = new_fs_zm_idxvec(num_indices)
    ALLOCATE(fs_zmdeco%offset(num_indices))
    CALL set_fs_zm_offset(num_indices, fs_zmdeco%offset)

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_fs_zm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)

      INTEGER :: k, j, r, p, confined_fleve

      confined_fleve = MIN(fleve, nlev)
      p = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          DO k = flevs, confined_fleve
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = k + nlev * (j-1)  - 1
          ENDDO
        ENDDO
      ENDDO
      new_fs_zm_idxvec = xt_idxvec_new( idx, nidx)
    END FUNCTION new_fs_zm_idxvec

    SUBROUTINE set_fs_zm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: k, jj, p

      ! offsets for zonal means

      ! data shape: (nflev, nflat)
      ! see: first two dimensions of fbm0 in mo_buffer_fft

      p = 0
      DO jj = 1, nflat
        DO k = 1, nflev
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          offset(p) = p -1
        ENDDO
      ENDDO
    END SUBROUTINE set_fs_zm_offset

  END SUBROUTINE setup_fs_zmdeco

  !
  ! Define global index list and offsets for spectral coefficients
  !
  SUBROUTINE setup_sp_sdeco(idx_description_type)
    INTEGER, INTENT(IN) :: idx_description_type
    INTEGER :: nm, nsm, snsp
    INTEGER, ALLOCATABLE :: sm(:), nmp(:), snn0(:), snnp(:)

    ! Spectral decomposition data
    nm  = ldc%nm   ! spectral truncation
    nsm = ldc%nsm  ! number of m wave numbers on PE
    snsp= ldc%snsp ! number of spectral (complex) coefficients on PE
    ALLOCATE(sm(nsm), snn0(nsm), snnp(nsm), nmp(nm+2))
    sm  = ldc%sm
    nmp = ldc%nmp  ! displacement of the first point of
                   ! m-columns with respect to the first point
                   ! of the first m-column
    snn0 = ldc%snn0
    snnp = ldc%snnp

    CALL setup_sp

  CONTAINS

    SUBROUTINE setup_sp
      INTEGER :: nlev_list(3), nidx, nstr, i

      ! array with numbers of vertical levels / tiles
      nlev_list = (/1, nlev, nlev+1/)

      ALLOCATE(sp_sdeco(MAXVAL(nlev_list)))
      ! the next line is only needed for some buggy compilers
      DO i = 1, SIZE(sp_sdeco)
        sp_sdeco(i)%idxlist = xt_idxlist_c2f(c_null_ptr)
      ENDDO

      DO i = 1, SIZE(nlev_list)
        ! generate index list
        IF (xt_is_null(sp_sdeco(nlev_list(i))%idxlist)) THEN
          nidx = 2 * snsp * nlev_list(i)
          SELECT CASE (idx_description_type)
            ! index list as vector
          CASE (idx_descr_vector)
            sp_sdeco(nlev_list(i))%idxlist = new_sp_idxvec(nidx, nlev_list(i))
            ! index list by stripes
          CASE (idx_descr_stripes)
            nstr = nsm * nlev_list(i)
            sp_sdeco(nlev_list(i))%idxlist = new_sp_idxstripes(nstr, nlev_list(i))
          CASE DEFAULT
            WRITE (message_text, '(a,i0,a)') &
                 'Unsupported global index desciption type: ', &
                 idx_description_type, ' (spectral data)'
            CALL finish('setup_yaxt_decomposition', message_text)
          END SELECT

          ! generate offsets
          ALLOCATE(sp_sdeco(nlev_list(i))%offset(nidx))
          CALL set_sp_offset(nidx, nlev_list(i), sp_sdeco(nlev_list(i))%offset)
        ENDIF
      ENDDO

    END SUBROUTINE setup_sp

    !
    ! spectral data
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_sp_idxvec(nidx, nlv)
      INTEGER, INTENT(IN) :: nidx, nlv
      INTEGER :: idx(nidx)
      INTEGER :: nsp, mp1, spgls, spgle
      INTEGER :: p, k, im, n, c

      ! global 3d spectral index space: sp3_index(1:2, 0:nsp-1, 1:nlev)
      ! sp3_index(c,n,k) := c-1 + 2*n + 2*nsp*(k-1)

      nsp = (nm + 1) * ( nm + 2) / 2

      p = 0
      DO k = 1, nlv
        DO im = 1, nsm
          mp1  = sm(im) + 1
          spgls = nmp(mp1) + snn0(im)
          spgle = spgls + snnp(im) - 1
          DO n = spgls, spgle
            DO c = 1, 2
              p = p + 1
              idx(p) = c + 2 * n + 2 * nsp * (k - 1) - 1
            ENDDO
          ENDDO
        END DO
      END DO
      new_sp_idxvec = xt_idxvec_new(idx, nidx)

    END FUNCTION new_sp_idxvec

    !
    ! spectral data
    ! describe decomposition by index stripes
    TYPE(xt_idxlist) FUNCTION new_sp_idxstripes(nstr, nlv)
      INTEGER, INTENT(IN) :: nstr, nlv
      TYPE(xt_stripe) :: s(nstr)
      INTEGER :: nsp, p, k, im, mp1, spgls

      nsp = (nm + 1) * ( nm + 2) / 2
      s(:)%stride = 1
      p = 0 ! strides index
      DO k = 1, nlv
        DO im = 1, nsm
          mp1  = sm(im) + 1
          spgls = nmp(mp1) + snn0(im)
          p = p + 1
          s(p)%start = 2 * spgls + 2 * nsp * (k - 1)
          s(p)%nstrides = 2 * snnp(im)
        END DO
      END DO
      new_sp_idxstripes = xt_idxstripes_new(s, nstr)

    END FUNCTION new_sp_idxstripes
    !
    ! spectral data
    ! set offsets
    SUBROUTINE set_sp_offset(nidx, nlv, offset)
      INTEGER, INTENT(IN)  :: nidx, nlv
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: p, k, n, c
      INTEGER :: im, mp1, spgls, spgle

      !
      ! local echam 3d spectral offsets: sp3_offset(1:nlv,1:2,1:snsp)
      ! sp3_offset(k,c,n) = k-1 + nlv*(c-1) + nlv*2*(n-1)
      !

      p = 0
      DO k = 1, nlv
        DO n = 1, snsp
          DO c = 1, 2
            p = p + 1
            offset(p) = k + nlv * (c - 1) + nlv * 2 * (n - 1) - 1
          ENDDO !k
        ENDDO !n
      ENDDO !c

    END SUBROUTINE set_sp_offset

  END SUBROUTINE setup_sp_sdeco

  !
  ! decomposition & offsets for spectral coeff. with m==0
  SUBROUTINE setup_sp_m0deco()
    INTEGER :: num_indices

    IF (ldc%nsnm0 == 0) THEN
      sp_m0deco = empty_deco
      RETURN
    ENDIF

    num_indices = nlev*ldc%nsnm0
    sp_m0deco%idxlist = new_sp_m0_idxvec(num_indices, nlev)
    ALLOCATE(sp_m0deco%offset(num_indices))
    CALL set_sp_m0_offset(num_indices, nlev, sp_m0deco%offset)

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_sp_m0_idxvec(nidx, global_kmax)
      INTEGER, INTENT(IN) :: nidx, global_kmax
      INTEGER :: idx(nidx)
      INTEGER :: n, p, k, nlev, nsnm0, snn0

      ! see scatter_sp0 for details

      nlev   = global_kmax
      nsnm0  = ldc%nsnm0
      snn0   = ldc%snn0(1)

      p = 0
      DO n = snn0+1,snn0+nsnm0
        DO k = 1, nlev
          p = p + 1
          idx(p) = k + nlev * (n-1)
        ENDDO
      END DO

      new_sp_m0_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_sp_m0_idxvec

    SUBROUTINE set_sp_m0_offset(nidx, global_kmax, offset)
      INTEGER, INTENT(IN)  :: nidx, global_kmax
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: p, k, n

      p = 0
      DO n = 1, ldc%nsnm0
        DO k = 1, global_kmax
          p = p + 1
          offset(p) = p-1
        ENDDO
      ENDDO

    END SUBROUTINE set_sp_m0_offset

  END SUBROUTINE setup_sp_m0deco

  !
  ! Define global index list and offsets for FFSL 2D-grid
  !
  SUBROUTINE setup_ffsl_2d_gdeco(idx_description_type)
    INTEGER, INTENT(IN) :: idx_description_type

    CALL setup_ffsl_2d

  CONTAINS

    SUBROUTINE setup_ffsl_2d
      INTEGER :: i, istate, nghost
      INTEGER,PARAMETER :: nghost_list(2) = (/ 1, 3 /)

      ffsl_2d_gdeco%idxlist  = my_idxlist(0)
      CALL my_offset(ffsl_2d_gdeco%offset, 0, 0)

      ffsl_2d_gz_gdeco(:,:)%idxlist = xt_idxlist_c2f(c_null_ptr)

      DO i = 1, SIZE(nghost_list)
        nghost = nghost_list(i)

        ! start state:
        istate = 0
        ffsl_2d_gz_gdeco(istate,nghost)%idxlist  = ffsl_2d_gdeco%idxlist
        CALL my_offset(ffsl_2d_gz_gdeco(istate,nghost)%offset, 0, nghost)

        ! end state:
        istate = 1
        ffsl_2d_gz_gdeco(istate,nghost)%idxlist = my_idxlist(nghost)
        CALL my_offset(ffsl_2d_gz_gdeco(istate,nghost)%offset, nghost, nghost)

      ENDDO


    END SUBROUTINE setup_ffsl_2d

    TYPE(xt_idxlist) FUNCTION my_idxlist(nghost) RESULT(y)
      INTEGER, INTENT(in) :: nghost
      INTEGER :: nidx, nstr, jlat1, jlat2

      jlat1 = MAX(1,    ffsl_gp_lat1 - nghost)
      jlat2 = MIN(nlat, ffsl_gp_lat2 + nghost)
      ! generate index list
      SELECT CASE(idx_description_type)
      CASE (idx_descr_vector)
        nidx = nlon * (jlat2-jlat1+1)
        y  = new_ffsl_2d_idxvec(nidx, jlat1, jlat2)
      CASE (idx_descr_stripes)
        nstr = (jlat2-jlat1+1)
        y  = new_ffsl_2d_idxstripes(nstr, jlat1, jlat2)
      CASE DEFAULT
        CALL finish('setup_yaxt_decomposition', &
             'Index list description type is not valid')
      END SELECT
    END FUNCTION my_idxlist

    !
    ! FFSL 2D data
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_ffsl_2d_idxvec(nidx, jlat1, jlat2)
      INTEGER, INTENT(IN) :: nidx, jlat1, jlat2
      INTEGER :: idx(nidx)
      INTEGER :: i, j, p

      p = 0
      DO j = jlat1, jlat2
        DO i = 1, nlon
          p = p + 1
          idx(p) =i + nlon * (j-1) - 1
        ENDDO
      ENDDO
      new_ffsl_2d_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ffsl_2d_idxvec

    !
    ! FFSL 2D data
    ! describe decomposition by index stripes
    TYPE(xt_idxlist) FUNCTION new_ffsl_2d_idxstripes(nstr, jlat1, jlat2)
      INTEGER, INTENT(IN) :: nstr, jlat1, jlat2
      TYPE(xt_stripe)     :: s(nstr)
      INTEGER :: j, p


      s(:)%stride = 1
      s(:)%nstrides = nlon
      p = 0
      DO j = jlat1, jlat2
        p = p + 1
        s(p)%start = 1 + nlon * (j-1) - 1
      ENDDO
      new_ffsl_2d_idxstripes = xt_idxstripes_new(s, nstr)
    END FUNCTION new_ffsl_2d_idxstripes

    SUBROUTINE my_offset(offset, nghost, nghost_data)
      INTEGER, ALLOCATABLE, INTENT(out) :: offset(:)
      INTEGER, INTENT(in) :: nghost, nghost_data
      INTEGER :: nidx, jlat1, jlat2, jlat1_data, jlat2_data

      jlat1 = MAX(1,    ffsl_gp_lat1 - nghost)
      jlat2 = MIN(nlat, ffsl_gp_lat2 + nghost)
      jlat1_data = ffsl_gp_lat1 - nghost_data
      jlat2_data = ffsl_gp_lat2 + nghost_data
      nidx = nlon * (jlat2-jlat1+1)
      ! generate offsets:
      ALLOCATE(offset(nidx))
      CALL set_ffsl_2d_offset(nidx, offset, jlat1, jlat2, jlat1_data, jlat2_data)

    END SUBROUTINE my_offset

    !
    ! FFSL 2D data
    ! set offsets
    SUBROUTINE set_ffsl_2d_offset(nidx, offset, jlat1, jlat2, jlat1_data, jlat2_data)
      INTEGER, INTENT(IN)  :: nidx, jlat1, jlat2, jlat1_data, jlat2_data
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: coords_off(nlon, jlat1_data:jlat2_data)
      INTEGER :: i, j, ic, jc, p

      IF (jlat1 < jlat1_data .OR. jlat2 > jlat2_data) CALL die('bad case', __LINE__)

      ! data offsets - here we take care of the reversed latitute iteration in ffsl space
      ! layout: ffsl_field(1:nlon,jlat1_data:jlat2_data)

      p = 0
      DO jc = jlat2_data, jlat1_data, -1 ! change in j-orientation
        DO ic = 1, nlon
          coords_off(ic, jc) = p
          p = p + 1
        ENDDO
      ENDDO

      p = 0
      DO j = jlat1, jlat2
        jc = j
        DO i = 1, nlon
          ic = i
          p = p + 1
          offset(p) = coords_off(ic, jc)
        ENDDO
      ENDDO

    END SUBROUTINE set_ffsl_2d_offset

  END SUBROUTINE setup_ffsl_2d_gdeco

  !
  ! Define global index list and offsets for FFSL 3D-grid
  !
  SUBROUTINE setup_ffsl_3d_gdeco(idx_description_type)
    INTEGER, INTENT(IN) :: idx_description_type
    INTEGER, ALLOCATABLE :: ffsl_kstack(:)

    CALL setup_ffsl_3d

  CONTAINS

    SUBROUTINE setup_ffsl_3d
      INTEGER,PARAMETER :: nghost_list(2) = (/ 1, 3 /)
      INTEGER :: i, istate, nghost
      INTEGER :: kstack(ldc%nlev)

      CALL get_ffsl_kstack(kstack, ffsl_nlev)
      ALLOCATE(ffsl_kstack(ffsl_nlev))
      ffsl_kstack = kstack(1:ffsl_nlev)

      ffsl_3d_gdeco%idxlist =  my_idxlist(0)
      CALL my_offset(ffsl_3d_gdeco%offset, 0, 0)

      ffsl_3d_gz_gdeco(:,:)%idxlist = xt_idxlist_c2f(c_null_ptr)

      DO i = 1, SIZE(nghost_list)
        nghost = nghost_list(i)

        ! start state:
        istate = 0
        ffsl_3d_gz_gdeco(istate,nghost)%idxlist  = ffsl_3d_gdeco%idxlist
        CALL my_offset(ffsl_3d_gz_gdeco(istate,nghost)%offset, 0, nghost)

        ! end state:
        istate = 1
        ffsl_3d_gz_gdeco(istate,nghost)%idxlist = my_idxlist(nghost)
        CALL my_offset(ffsl_3d_gz_gdeco(istate,nghost)%offset, nghost, nghost)

      ENDDO

#if 0
      INTEGER :: nidx, nstr
      INTEGER :: kstack(ldc%nlev)

      !CALL get_ffsl_kstack(kstack, ffsl_nlev)
      !ALLOCATE(ffsl_kstack(ffsl_nlev))
      !ffsl_kstack = kstack(1:ffsl_nlev)
      nidx = nlon * ffsl_nlat * ffsl_nlev

      ! generate index list
      SELECT CASE(idx_description_type)
        CASE (idx_descr_vector)
          ffsl_3d_gdeco%idxlist = new_ffsl_3d_idxvec(nidx)
        CASE (idx_descr_stripes)
          nstr = ffsl_nlev * ffsl_nlat
          ffsl_3d_gdeco%idxlist = new_ffsl_3d_idxstripes(nstr)
        CASE DEFAULT
          CALL finish('setup_yaxt_decomposition', &
              'Index list description type is not valid')
      END SELECT

      ! generate offsets
      ALLOCATE(ffsl_3d_gdeco%offset(nidx))
      CALL set_ffsl_3d_offset(nidx, ffsl_3d_gdeco%offset)
#endif
    END SUBROUTINE setup_ffsl_3d

    TYPE(xt_idxlist) FUNCTION my_idxlist(nghost) RESULT(y)
      INTEGER, INTENT(in) :: nghost
      INTEGER :: nidx, nstr, jlat1, jlat2
      INTEGER :: kstack(ldc%nlev)

      jlat1 = MAX(1,    ffsl_gp_lat1 - nghost)
      jlat2 = MIN(nlat, ffsl_gp_lat2 + nghost)
      ! generate index list
      SELECT CASE(idx_description_type)
      CASE (idx_descr_vector)
        nidx = nlon * (jlat2-jlat1+1) * ffsl_nlev
        y  = new_ffsl_3d_idxvec(nidx, jlat1, jlat2)
      CASE (idx_descr_stripes)
        nstr = (jlat2-jlat1+1) * ffsl_nlev
        y  = new_ffsl_3d_idxstripes(nstr, jlat1, jlat2)
      CASE DEFAULT
        CALL finish('setup_yaxt_decomposition', &
             'Index list description type is not valid')
      END SELECT
    END FUNCTION my_idxlist

    SUBROUTINE my_offset(offset, nghost, nghost_data)
      INTEGER, ALLOCATABLE, INTENT(out) :: offset(:)
      INTEGER, INTENT(in) :: nghost, nghost_data
      INTEGER :: nidx, jlat1, jlat2, jlat1_data, jlat2_data

      jlat1 = MAX(1,    ffsl_gp_lat1 - nghost)
      jlat2 = MIN(nlat, ffsl_gp_lat2 + nghost)
      jlat1_data = ffsl_gp_lat1 - nghost_data
      jlat2_data = ffsl_gp_lat2 + nghost_data
      nidx = nlon * (jlat2-jlat1+1) * ffsl_nlev
      ! generate offsets:
      ALLOCATE(offset(nidx))
      CALL set_ffsl_3d_offset(nidx, offset, jlat1, jlat2, jlat1_data, jlat2_data)

    END SUBROUTINE my_offset

    !
    ! FFSL 3D data
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_ffsl_3d_idxvec(nidx, jlat1, jlat2)
      INTEGER, INTENT(IN) :: nidx, jlat1, jlat2
      INTEGER :: idx(nidx)

      INTEGER :: i, j, k, kk, p

      p = 0
      DO kk = 1, ffsl_nlev
        k = ffsl_kstack(kk)
        DO j = jlat1, jlat2
          DO i = 1, nlon
            p = p + 1
            idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
          ENDDO
        ENDDO
      ENDDO
      new_ffsl_3d_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ffsl_3d_idxvec
    !
    ! FFSL 3D data
    ! describe decomposition by index stripes
    TYPE(xt_idxlist) FUNCTION new_ffsl_3d_idxstripes(nstr, jlat1, jlat2)
      INTEGER, INTENT(IN) :: nstr, jlat1, jlat2
      TYPE(xt_stripe) :: s(nstr)

      INTEGER :: j, k, kk, p

      s(:)%stride = 1
      s(:)%nstrides = nlon
      p = 0
      DO kk = 1, ffsl_nlev
        k = ffsl_kstack(kk)
        DO j = jlat1, jlat2
          p = p + 1
          s(p)%start = 1 + nlon * ( (j-1) + nlat * (k-1) ) - 1
        ENDDO
      ENDDO
      new_ffsl_3d_idxstripes = xt_idxstripes_new(s, nstr)
    END FUNCTION new_ffsl_3d_idxstripes
    !
    ! FFSL 3D data
    ! set offsets
    SUBROUTINE set_ffsl_3d_offset(nidx, offset, jlat1, jlat2, jlat1_data, jlat2_data)
      INTEGER, INTENT(IN) :: nidx, jlat1, jlat2, jlat1_data, jlat2_data
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: coords_off(nlon, jlat1_data:jlat2_data, ffsl_nlev)
      INTEGER :: i, j, kk, ic, jc, kc, p

      IF (jlat1 < jlat1_data .OR. jlat2 > jlat2_data) CALL die('bad case', __LINE__)

      ! access to data offsets via coords
      p = 0
      DO kc = 1, ffsl_nlev
        DO jc = jlat2_data, jlat1_data, -1 ! change in j-orientation
          DO ic = 1, nlon
            coords_off(ic, jc, kc) = p
            p = p + 1
          ENDDO
        ENDDO
      ENDDO

      ! offsets in index order (must match order in idxlist definition)
      p = 0
      DO kk = 1, ffsl_nlev
        kc = kk
        DO j = jlat1, jlat2
          jc = j
          DO i = 1, nlon
            ic = i
            p = p + 1
            offset(p) = coords_off(ic, jc, kc)
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE set_ffsl_3d_offset

  END SUBROUTINE setup_ffsl_3d_gdeco


  !
  ! Define global index list and offsets for LS 3D
  !
  SUBROUTINE setup_ls_3d_gdeco()
    INTEGER :: num_indices

    num_indices = nlat * nflevp1 * 2 * ldc%nlm

    ! generate index list
    ls_3d_gdeco%idxlist = new_ls_3d_idxvec(num_indices)

    ! generate offsets
    ALLOCATE(ls_3d_gdeco%offset(num_indices))
    CALL set_ls_3d_offset(num_indices, ls_3d_gdeco%offset)

  CONTAINS

    ! Legendre space indices for 3d data,
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_ls_3d_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)

      INTEGER :: i, j, k, ii, p
      INTEGER, POINTER :: intr(:) ! index array

      intr => ldc%intr
      p = 0
      DO j = 1, nlat
        DO k = flevs, fleve
          DO ii = 1, 2*ldc%nlm
            i = intr(ii)
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
            IF (idx(p) <0 .OR. idx(p) > nlon*nlat*(nlev+1)) CALL die("idx out of range",__LINE__)
          ENDDO
        ENDDO
      ENDDO
      new_ls_3d_idxvec = xt_idxvec_new( idx, nidx)

    END FUNCTION new_ls_3d_idxvec

    !
    ! Legendre space 3D data
    ! set offsets
    SUBROUTINE set_ls_3d_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: ii, j, k, p

      ! assumption:
      ! fftl shape = (nlm*2, nflevp1, nlat , nvar)

      p = 0
      DO j = 1, nlat
        DO k = flevs, fleve
          DO ii = 1, 2*ldc%nlm
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            offset(p) = p - 1
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE set_ls_3d_offset

  END SUBROUTINE setup_ls_3d_gdeco

  !
  ! Define global index list and offsets for LS zm
  !
  SUBROUTINE setup_ls_zmdeco()
    INTEGER :: num_indices

    IF (ldc%nlnm0 > 0) THEN
      num_indices = nflev*nlat
      ls_zmdeco%idxlist = new_ls_zm_idxvec(num_indices)
      ALLOCATE(ls_zmdeco%offset(num_indices))
      CALL set_ls_zm_offset(num_indices, ls_zmdeco%offset)
    ELSE
      ls_zmdeco%idxlist = xt_idxempty_new()
      ALLOCATE(ls_zmdeco%offset(1))
      ls_zmdeco%offset(1)=0
    ENDIF

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_ls_zm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: j, k, p, confined_fleve

      confined_fleve = MIN(fleve, nlev)

      p = 0
      DO j = 1, nlat
        DO k = flevs, confined_fleve
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          idx(p) =  k + nlev * (j-1)  - 1
        ENDDO
      ENDDO

      new_ls_zm_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ls_zm_idxvec

    SUBROUTINE set_ls_zm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: j, k, p

      ! offsets for zonal means
      ! data shape: (nflev, nlat)
      ! see: first two dimensions of lbm0 in mo_buffer_fft

      p = 0
      DO j = 1, nlat
        DO k = 1, nflev
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          offset(p) = p -1
        ENDDO
      ENDDO
    END SUBROUTINE set_ls_zm_offset

  END SUBROUTINE setup_ls_zmdeco

  ! legendre decomposition of spectral coeff.:
  SUBROUTINE setup_ls_sdeco()

    CALL setup_ls

  CONTAINS

    SUBROUTINE setup_ls
      INTEGER :: nidx, global_kmax, nk, ke

      ALLOCATE(ls_sdeco(nlev:nlev+1))

      DO global_kmax = nlev, nlev+1

        ke     = MIN(ldc%lleve, global_kmax)
        nk     = ke - ldc%llevs + 1
        nidx = nk * 2 * ldc%lnsp

        ! generate index list
        ls_sdeco(global_kmax)%idxlist = new_ls3_idxvec(nidx, global_kmax)

        ! generate offsets
        ALLOCATE(ls_sdeco(global_kmax)%offset(nidx))
        CALL set_ls3_offset(nidx, ls_sdeco(global_kmax)%offset, global_kmax)

      ENDDO
    END SUBROUTINE setup_ls

    TYPE(xt_idxlist) FUNCTION new_ls3_idxvec(nidx, global_kmax)
      INTEGER, INTENT(IN) :: nidx, global_kmax
      INTEGER :: idx(nidx)

      ! global 3d legendre index space: ls3_index(1:2, 0:nsp-1, 1:nlev)
      ! ls3_index(c,n,k) := c-1 + 2*n + 2*nsp*(k-1)
      ! Note: sp3 and ls3 use the same global index space

      INTEGER :: mp, nsp
      INTEGER :: c, p, k, ke, im, mp1, mpgl, np, n

      nsp = control_nsp
      ke     = MIN(ldc%lleve, global_kmax)

      p = 0
      DO k = ldc%llevs, ke
        mp=0
        DO im=1,ldc%nlm
          mp1  = ldc%lm(im)+1
          np   = ldc%nnp(mp1)
          mpgl = ldc%nmp(mp1)
          DO n = mpgl+1, mpgl+np
            DO c = 1, 2
              p = p + 1
              idx(p)= c-1 + 2 * (n-1) + 2 * nsp * (k-1)
            ENDDO
          ENDDO
          mp = mp + np
        ENDDO
      ENDDO

      new_ls3_idxvec = xt_idxvec_new( idx, nidx)
    END FUNCTION new_ls3_idxvec

    SUBROUTINE set_ls3_offset(nidx, offset, global_kmax)
      INTEGER, INTENT(IN) :: nidx, global_kmax
      INTEGER, INTENT(out) :: offset(nidx)
      INTEGER :: mp1, llevs, lleve, np, ke, nk
      INTEGER :: p, k, im, n, c

      llevs  = ldc%llevs
      lleve  = ldc%lleve
      ke = MIN (lleve,global_kmax)
      nk = ke - llevs + 1

      p = 0
      DO k = 1, nk
        DO n = 0, ldc%lnsp-1
          DO c = 1, 2
            p = p + 1
            offset(p) = k + nk*(c-1) + (n-0)*2*nk - 1
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE set_ls3_offset

  END SUBROUTINE setup_ls_sdeco

  SUBROUTINE setup_ls_m0deco()
    INTEGER :: num_indices, global_kmax, nk, ke

    global_kmax = nlev
    ke     = MIN(ldc%lleve, global_kmax)
    nk     = ke - ldc%llevs + 1
    num_indices = nk * ldc%nlnm0
    ls_m0deco%idxlist = new_ls_m0_idxvec(num_indices)
    ALLOCATE(ls_m0deco%offset(num_indices))
    CALL set_ls_m0_offset(num_indices,ls_m0deco%offset)

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_ls_m0_idxvec(nidx)
      INTEGER, INTENT(in) :: nidx
      INTEGER :: idx(nidx), p, k, n
      p = 0
      DO n = 1, ldc%nlnm0
        DO k = ldc%llevs,ke
          p = p + 1
          idx(p) = k + global_kmax*(n-1)
        ENDDO
      ENDDO

      new_ls_m0_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ls_m0_idxvec

    SUBROUTINE set_ls_m0_offset(nidx, offset)
      INTEGER, INTENT(in) :: nidx
      INTEGER, INTENT(out) :: offset(nidx)
      INTEGER :: p, k, n
      p = 0
      DO n = 1, ldc%nlnm0
        DO k = ldc%llevs,ke
          p = p + 1
          offset(p) = p-1
        ENDDO
      ENDDO
    END SUBROUTINE set_ls_m0_offset

  END SUBROUTINE setup_ls_m0deco


    SUBROUTINE get_ffsl_kstack(kstack, kstack_n)
      INTEGER, INTENT(OUT) :: kstack(:)
      INTEGER, INTENT(OUT) :: kstack_n
      INTEGER :: dest_idx(2), dest_set_b
      INTEGER :: k, ia, my_idx

      my_idx = ldc%pe + 1

      ! echam::mo_transpose send-logic
      kstack_n = 0
      kstack = 0
      DO k = 1, ldc%nlev
        !
        ! for PEs with odd set_a:
        !   Northern Hemisphere sent to same set_a
        !   Southern Hemisphere sent to set_a+1 (unless set_a  ==  nproca)
        ! for PEs with even set_a:
        !   Northern Hemisphere sent to set_a-1
        !   Southern Hemisphere sent to same set_a
        !
        dest_set_b  =  MOD(k-1,nprocb)+1 ! set_b of receiving PE for this k
        IF( MOD(ldc%set_a,2)  ==  1 ) THEN
          dest_idx(1) = ldc%mapmesh(dest_set_b,ldc%set_a) ! N-region goes to my nprocb proc-family
          ia = MIN(ldc%set_a+1,ldc%nproca)
          dest_idx(2) = ldc%mapmesh(dest_set_b,ia) ! S-region goes preferably to northern nprocb proc-family
        ELSE
          dest_idx(1) = ldc%mapmesh(dest_set_b,ldc%set_a-1) ! N-region goes to southern nprocb proc-family
          dest_idx(2) = ldc%mapmesh(dest_set_b,ldc%set_a) ! S-region goes to my nprocb proc-family
        ENDIF

        IF (dest_idx(1) == my_idx .OR. dest_idx(2) == my_idx) THEN
          kstack_n = kstack_n + 1
          kstack(kstack_n) = k
        ENDIF

      ENDDO
    END SUBROUTINE get_ffsl_kstack

  ! set some required sub communicators of the given model communicator
  SUBROUTINE set_comms(model_comm)
    USE mpi
    INTEGER, INTENT(IN) :: model_comm  ! communicator of application without I/O servers
    INTEGER :: ierror, color, key
    IF (debug_parallel < 0) THEN
      ! use own echam communicator to avoid tag collision:
      CALL mpi_comm_dup(model_comm, ab_comm, ierror)
      IF (ierror /= MPI_SUCCESS) CALL die('mpi_comm_dup failed', __LINE__)
    ELSE
      ! If echam runs in debug_parallel mode then the process space size is nproca x nprocb + 1.
      ! In that case we want to split the process space into two spaces of sizes nproca x nprocb and 1 x 1.
      color = ldc%spe
      key = ldc%pe
      CALL mpi_comm_split(model_comm, color, key, ab_comm, ierror)
      IF (ierror /= MPI_SUCCESS) CALL die('mpi_comm_split failed', __LINE__)
    ENDIF

    ! nproca x nprocb space:
    CALL mpi_comm_size(ab_comm, ab_size, ierror)
    CALL mpi_comm_rank(ab_comm, ab_rank, ierror)
    IF (ab_size /= nproca * nprocb) CALL die('set_comms: internal error (1)',__LINE__)

  END SUBROUTINE set_comms

  SUBROUTINE exchange(redist, f, g)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(dp), TARGET, INTENT(in) :: f(*)
    REAL(dp), TARGET, INTENT(out) :: g(*)
    CALL xt_redist_s_exchange1(redist, C_LOC(f), C_LOC(g));
  END SUBROUTINE exchange

  SUBROUTINE set_gp_2d_testdata(f)
    REAL(dp), INTENT(OUT) :: f(:,:) !f(ia,ib)
    INTEGER :: i, ig, j, jg, ia, ib, r
    INTEGER :: ival

    f = 0.0_dp

    ib = 1
    ia = 0
    DO r = 1, 2
      DO jg = glats(r), glate(r)
        j = ilat_map(jg)
        DO ig = glons(r), glone(r)
          i = ilon_map(ig)
          ia = ia +1
          IF (ia > nproma) THEN
            ib = ib + 1
            ia = 1
          ENDIF
          ival = i + (j-1)*1000
          IF (ia > SIZE(f,1) .OR. ib > SIZE(f,2)) CALL die("bad case", __LINE__)
          f(ia,ib) = REAL(ival, dp)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE set_gp_2d_testdata

  SUBROUTINE set_gp_3d_testdata(f)
    REAL(dp), INTENT(OUT) :: f(:,:,:) !f(ia,k,ib)
    INTEGER :: i, ig, j, jg, k, ia, ib, r
    INTEGER :: ival

    DO k = 1, nlev
      ib = 1
      ia = 0
      DO r = 1, 2
        DO jg = glats(r), glate(r)
          j = ilat_map(jg)
          DO ig = glons(r), glone(r)
            i = ilon_map(ig)
            ia = ia +1
            IF (ia > nproma) THEN
              ib = ib + 1
              ia = 1
            ENDIF
            ival = i + ( (j-1) + (k-1)*1000)*1000
            f(ia,k,ib) = REAL(ival, dp)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE set_gp_3d_testdata

  ! module specific termination subroutine,
  ! simplifies code
  SUBROUTINE die(reason, line)
    USE mo_exception,  ONLY: finish
    USE yaxt,          ONLY: xt_finalize, xt_abort
    CHARACTER(len=*), INTENT(in) :: reason
    INTEGER, INTENT(in) :: line
    CHARACTER(len=*), PARAMETER :: file = __FILE__

    WRITE(message_text,*) reason,', file ',file,', line',line
    CALL xt_finalize()
    CALL xt_abort(TRIM(reason), __FILE__, line)
    CALL finish('mo_echam_yaxt', message_text)
    WRITE(0,*) 'model_abort: fallback stop'
    STOP
  END SUBROUTINE die

#ifdef CHECK_MEM_STRIDES
  SUBROUTINE check_mem_strides_4d(f, memstride, element_size, line)
    USE mpi
    CHARACTER(len=*), PARAMETER :: prefix = 'check_mem_strides_4d: '
    REAL(dp), INTENT(in) :: f(:,:,:,:)
    INTEGER, INTENT(in) :: memstride(:) ! vector of memory strides in units of element_size
    INTEGER, INTENT(in) :: element_size, line
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a, b
    INTEGER :: ierror

    IF (SIZE(memstride) /= 4) CALL die(prefix//'bad usage',__LINE__)

    IF (SIZE(f)<2) RETURN
    CALL mpi_get_address(f(1,1,1,1), a, ierror)

    IF (SIZE(f,1) >= 2) THEN
      CALL mpi_get_address(f(2,1,1,1), b, ierror)
      IF ( b-a /= INT(memstride(1)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 1', line)
    ENDIF

    IF (SIZE(f,2) >= 2) THEN
      CALL mpi_get_address(f(1,2,1,1), b, ierror)
      IF ( b-a /= INT(memstride(2)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 2', line)
    ENDIF

    IF (SIZE(f,3) >= 2) THEN
      CALL mpi_get_address(f(1,1,2,1), b, ierror)
      IF ( b-a /= INT(memstride(3)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 3', line)
    ENDIF

    IF (SIZE(f,4) >= 2) THEN
      CALL mpi_get_address(f(1,1,1,2), b, ierror)
      IF ( b-a /= INT(memstride(4)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 4', line)
    ENDIF

  END SUBROUTINE check_mem_strides_4d

  SUBROUTINE check_mem_strides_3d(f, memstride, element_size, line)
    USE mpi
    CHARACTER(len=*), PARAMETER :: prefix = 'check_mem_strides_3d: '
    REAL(dp), INTENT(in) :: f(:,:,:)
    INTEGER, INTENT(in) :: memstride(:) ! vector of memory strides in units of element_size
    INTEGER, INTENT(in) :: element_size, line
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a, b
    INTEGER :: ierror

    IF (SIZE(memstride) /= 3) CALL die(prefix//'bad usage',__LINE__)

    IF (SIZE(f)<2) RETURN
    CALL mpi_get_address(f(1,1,1), a, ierror)

    IF (SIZE(f,1) >= 2) THEN
      CALL mpi_get_address(f(2,1,1), b, ierror)
      IF ( b-a /= INT(memstride(1)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 1', line)
    ENDIF

    IF (SIZE(f,2) >= 2) THEN
      CALL mpi_get_address(f(1,2,1), b, ierror)
      IF ( b-a /= INT(memstride(2)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 2', line)
    ENDIF

    IF (SIZE(f,3) >= 2) THEN
      CALL mpi_get_address(f(1,1,2), b, ierror)
      IF ( b-a /= INT(memstride(3)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 3', line)
    ENDIF
  END SUBROUTINE check_mem_strides_3d

  SUBROUTINE check_mem_strides_2d(f, memstride, element_size, line)
    USE mpi
    CHARACTER(len=*), PARAMETER :: prefix = 'check_mem_strides_2d: '
    REAL(dp), INTENT(in) :: f(:,:)
    INTEGER, INTENT(in) :: memstride(:) ! vector of memory strides in units of element_size
    INTEGER, INTENT(in) :: element_size, line
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a, b
    INTEGER :: ierror

    IF (SIZE(memstride) /= 2) CALL die(prefix//'bad usage',__LINE__)

    IF (SIZE(f)<2) RETURN
    CALL mpi_get_address(f(1,1), a, ierror)

    IF (SIZE(f,1) >= 2) THEN
      CALL mpi_get_address(f(2,1), b, ierror)
      IF ( b-a /= INT(memstride(1)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 1', line)
    ENDIF

    IF (SIZE(f,2) >= 2) THEN
      CALL mpi_get_address(f(1,2), b, ierror)
      IF ( b-a /= INT(memstride(2)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 2', line)
    ENDIF
  END SUBROUTINE check_mem_strides_2d

#endif /* CHECK_MEM_STRIDES */

#endif /* HAVE_YAXT */

END MODULE mo_echam_yaxt
