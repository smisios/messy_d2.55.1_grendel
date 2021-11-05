#define _DIAGONALEX_
! **************************************************************************
MODULE messy_main_decomp_bi
  ! **************************************************************************
  ! Author: Andreas Baumgaertner,  MPI-C, UCL, 2010
  ! TO DO:
  ! - move DC_xx from main_channel_bi to here
  ! - can we make ECHAM gp etc instead of total e5?
  ! - move indx0 from main_mpi_bi to here?
  ! - several dc variables rather than struct? array not possible...
  ! - move structs to core

#if defined(ECHAM5)
  USE mo_decomposition, ONLY: decomp_e5 => pe_decomposed
#endif
  USE messy_main_mpi_bi, ONLY: p_nprocs
  USE messy_main_blather_bi, ONLY: error_bi

  ! SMCL
  USE messy_main_decomp

  IMPLICIT NONE

  PRIVATE

  ! -------------------------------------------------------------------
  ! DECOMPOSITION TYPES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! DC_C1  -- CESM1 GRIDPOINT DECOMPOSITION
  TYPE decomp_c1
     INTEGER          :: nlon      ! global number of longitudes = nlon
     INTEGER          :: nlat      ! global number of latitudes  = nlat
     INTEGER          :: nglat     ! number of latitudes  on PE
     INTEGER          :: nglon     ! number of longitudes on PE
     !INTEGER          :: glats     ! start value of latitudes
     !INTEGER          :: glate     ! end value of latitudes
     !INTEGER          :: glons     ! start value of longitudes
     !INTEGER          :: glone     ! end value of longitudes
     INTEGER ,POINTER :: glat(:)   ! global latitude index N->S   ! dummy
     INTEGER ,POINTER :: glon(:)   ! offset to global longitude   ! dummy
     !INTEGER          :: nprocb    ! #of PEs for dimension that counts longitudes
     !INTEGER          :: nproca    ! #of PEs for dimension that counts latitudes
     INTEGER          :: pe        ! PE id
     !INTEGER          :: set_b  = 1   ! PE id in direction of longitudes
     !INTEGER          :: set_a  = 1   ! PE id in direction of latitudes
     !INTEGER ,POINTER :: mapmesh(:,:) ! indirection array mapping from a
                                      ! logical 2-d mesh to the processor index
                                      ! numbers
     ! Irregular grid
     INTEGER          :: ngpblks      ! number of rows
     INTEGER          :: nproma       ! number of columns
     INTEGER,POINTER  :: kproma(:)    ! number of columns in this row
     LOGICAL          :: lreg=.TRUE.  ! flag for regular grid
  END TYPE decomp_c1
  PUBLIC :: decomp_c1
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! DC_GP  -- ECHAM5 GRIDPOINT
  !TYPE decomp_gp
  !END TYPE decomp_gp
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! DUMMY FOR FURTHER DECOMPOSITION
  TYPE decomp_nn
     INTEGER         :: scales
     INTEGER         :: merx
  END TYPE decomp_nn
  PUBLIC :: decomp_nn
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! DUMMY FOR FURTHER DECOMPOSITION
  TYPE decomp_cosmo
     INTEGER         :: NOT
     INTEGER         :: REQUIRED
  END TYPE decomp_cosmo
  PUBLIC :: decomp_cosmo
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! TYPE FOR ALL DECOMPOSITIONS
  ! -------------------------------------------------------------------
  TYPE decomp
     LOGICAL                   :: lreg = .TRUE.
     INTEGER                   :: pe     ! PE id
#if defined(ECHAM5)
     TYPE (decomp_e5)          :: e5     ! ECHAM5 (gp etc.)
     TYPE (decomp_e5), POINTER :: e5g(:) ! global
#endif
#if defined(COSMO)
     TYPE (decomp_cosmo)       :: cosmo     ! COSMO
#endif
#if defined(CESM1)
     TYPE (decomp_c1)          :: c1     ! CESM1 gridpoint DC_C1
     TYPE (decomp_c1), POINTER :: c1g(:) ! global
#endif
     TYPE (decomp_nn) :: nn ! Clever and new
  END TYPE decomp
  PUBLIC :: decomp

  TYPE(decomp), PUBLIC, SAVE :: dc
  ! -------------------------------------------------------------------


  ! CALLED FROM initialize.f90 in ECHAM5, NOT messy_main_control_bi
  PUBLIC :: main_decomp_setup

  !PRIVATE :: decompose_e5
  !PRIVATE :: decompose_cg

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE main_decomp_setup

!    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
#if defined(CESM1)
    USE messy_main_tools,      ONLY: find_next_free_unit
#endif
#ifdef CESM1
    USE pmgrid,    ONLY: plat, plon
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_decomp_setup'

#if defined(CESM1)
    INTEGER :: iou, status

!!$    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_decomp_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
!!$    END IF

!!$    ! BROADCAST RESULTS
!!$    CALL p_bcast(NPX,     p_io)
!!$    CALL p_bcast(NPY,     p_io)
#endif

#if defined(ECHAM5)
    ALLOCATE(dc%e5g(p_nprocs))
    CALL decompose_e5(dc%e5, dc%e5g)

    NPY = dc%e5%nproca
    NPX = dc%e5%nprocb
    dc%pe = dc%e5%pe
#endif
#if defined(COSMO)
!    CALL decompose_c4(dc%c4)
#endif
#if defined(CESM1)
   ! IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_decomp_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
   ! END IF

    ! BROADCAST RESULTS ! can't broadcast yet
   ! CALL p_bcast(NLAT,     p_io)
  !  CALL p_bcast(NLON,     p_io)
#elif ECHAM5
    ! TAKE NLAT & NLON FROM ECHAM
    NLAT = dc%e5%nlat
    NLON = dc%e5%nlon
#endif
#ifdef CESM1
    ! TAKE NLAT & NLON FROM CESM (overwrite namelist)
    nlat = PLAT
    nlon = PLON
#endif

#ifdef CESM1
    ALLOCATE(dc%c1g(p_nprocs))
    CALL decompose_c1(dc%c1, dc%c1g)
    dc%pe = dc%c1%pe
#endif

  END SUBROUTINE main_decomp_setup
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

#if defined(ECHAM5)
   ! -------------------------------------------------------------------
  SUBROUTINE decompose_e5(dcl, dcg)

    USE mo_decomposition, ONLY: decomp_e5 => pe_decomposed      &
                              , dcg_e5 => global_decomposition  &
                              , dcl_e5 => local_decomposition
    IMPLICIT NONE

    TYPE(decomp_e5), INTENT(INOUT)               :: dcl
    TYPE(decomp_e5), DIMENSION(:), INTENT(INOUT) :: dcg

    ! echam5/src/init_decomposition.f90, normally called from initialize.f90
    CALL init_decomposition

    ! copy local ECHAM decomposition
    dc%e5 = dcl_e5
    ! copy global ECHAM decomposition
    dcg = dcg_e5

  END SUBROUTINE decompose_e5
  ! -------------------------------------------------------------------
#endif

#if defined(CESM1)
  ! -------------------------------------------------------------------
  SUBROUTINE main_decomp_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (SMIL)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CPL/ NLAT, NLON

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_decomp_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! CHECK NAMELIST
    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_decomp_read_nml_cpl
  ! -------------------------------------------------------------------
#endif

#ifdef CESM1
  SUBROUTINE decompose_c1(dcl, dcg)
    !
    USE messy_main_mpi_bi,  ONLY: p_nprocs, p_pe, p_parallel_io, p_bcast, p_io
    USE ppgrid,             ONLY: begchunk, endchunk
    USE phys_grid,          ONLY: lchunks
    USE pmgrid,    ONLY: plat, plon

    IMPLICIT NONE

    TYPE(decomp_c1), INTENT(INOUT)               :: dcl
    TYPE(decomp_c1), DIMENSION(:), INTENT(INOUT) :: dcg
    !
    ! local variables
    !
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_decomp_decompose_c1'
    INTEGER :: p, ngl, i, i2, inp, inprest
    INTEGER :: nptrlat(NPY+1), nptrlon(NPX+1)

    real(dp) :: test_dp=0._dp
    real(dp) :: test_dp5(5) = 0._dp
    logical :: test_l=.FALSE.
    CHARACTER(LEN=4) :: test_c = 'test'

!!$    ! broadcast test
!!$    p = 0
!!$    IF (p_parallel_io) then
!!$       p = 3
!!$    endif
!!$    CALL p_bcast(p,p_io)
!!$    write(*,*) 'broadcast test: p_parallel_io,p_io, p_pe,p ',p_parallel_io,p_io, p_pe,p
!!$    IF (p_parallel_io) then
!!$       test_dp = 2
!!$    endif
!!$    CALL p_bcast(test_dp,p_io)
!!$    write(*,*) 'broadcast test: p_parallel_io,p_io, p_pe,test_dp ',p_parallel_io,p_io, p_pe,test_dp
!!$    IF (p_parallel_io) then
!!$       test_dp5(4) = 2
!!$    endif
!!$    CALL p_bcast(test_dp5,p_io)
!!$    write(*,*) 'broadcast test: p_parallel_io,p_io, p_pe,test_dp5 ',p_parallel_io,p_io, p_pe,test_dp5
!!$    IF (p_parallel_io) then
!!$       test_l = .TRUE.
!!$    endif
!!$    CALL p_bcast(test_l,p_io)
!!$    write(*,*) 'broadcast test: p_parallel_io,p_io, p_pe,test_l ',p_parallel_io,p_io, p_pe,test_l
!!$    IF (p_parallel_io) then
!!$       test_c = 'ney1'
!!$    endif
!!$    CALL p_bcast(test_c,p_io)
!!$    write(*,*) 'broadcast test: p_parallel_io,p_io, p_pe,test_c ',p_parallel_io,p_io, p_pe,test_c



    ! First, populate dcg (copy to dcl later)
    p = 0
    DO i = 1, p_nprocs
       dcg(i)%pe    = p                 ! local PE id, 0:p_nprocs-1

       dcg(i)%nlon = nlon
       dcg(i)%nlat = nlat
       !dcg(i)%nglon = ! dummy
       !dcg(i)%nglat = ! dummy
       !dcg(i)%glat  = ! dummy
       !dcg(i)%glon  = ! dummy
       !dcg(i)%nproca = NPY ! dummy
       !dcg(i)%nprocb = NPX ! dummy
       dcg(i)%nproma = PCOLS
       dcg(i)%ngpblks = endchunk-begchunk+1
       ALLOCATE(dcg(i)%kproma(dcg(i)%ngpblks))
       dcg(i)%kproma(:) = lchunks(:)%ncols
!!$       write(*,*) 'messy_main_decomp_bi: i,plon,PLON,nlon,nlat,p,p_nprocs,ngpblks,nproma,p_parallel_io,p_pe ', &
!!$                                         i,plon,PLON,nlon,nlat,p,p_nprocs,dcg(i)%ngpblks,dcg(i)%nproma,p_parallel_io,p_pe
       p = p + 1
    END DO

    ! copy global decomposition table entry to local decomposition
    DO i = 1, p_nprocs
       IF (dcg(i)%pe == p_pe) THEN
          dcl%pe = dcg(i)%pe
          dcl = dcg(i)
       END IF
    END DO
  END SUBROUTINE decompose_c1
#endif

! **************************************************************************
END MODULE messy_main_decomp_bi
! **************************************************************************
