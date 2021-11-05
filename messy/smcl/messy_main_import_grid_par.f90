MODULE MESSY_MAIN_IMPORT_GRID_PAR

IMPLICIT NONE

  INTEGER, PARAMETER, PUBLIC :: RGSTAT_NULL = -1

  TYPE t_mpi_def
     INTEGER :: rank  = -1 
     INTEGER :: nproc = -1 
     INTEGER :: comm  = -1 
  END TYPE t_mpi_def
  PUBLIC :: t_mpi_def

  TYPE(t_mpi_def), PUBLIC, SAVE :: my_mpi

  INTEGER,  PUBLIC,SAVE         :: me = -1
  INTEGER,  PUBLIC,SAVE         :: nproc, comm
  INTEGER,  PUBLIC,PARAMETER    :: MAX_working_PEs=16
  INTEGER,  PUBLIC,SAVE         :: nproc_work  = 1 ! um_ak_20140327 =1 added
  LOGICAL,  PUBLIC,SAVE         :: i_am_worker = .FALSE.
  INTEGER,  PUBLIC,DIMENSION(:), ALLOCATABLE, SAVE :: pe_list

  PUBLIC  :: EXPAND_FILENAME
  PUBLIC  :: DISTRIBUTE_VARS_ON_PES
  PUBLIC  :: INIT_PIMPGRID
  PUBLIC  :: INIT_PARALLEL 

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE DISTRIBUTE_VARS_ON_PES(ivar, ovar, scl, namatts, RGT, RGTstr, RGrange)

  USE messy_main_grid_netcdf,   ONLY: GRD_MAXSTRLEN, t_ncatt_array &
                                    , INIT_NCATT, COPY_NCATT
  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: DP

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ivar(:)   ! input variable names
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ovar(:)   ! output variable names
  REAL(dp)                    , POINTER :: scl(:)    ! scaling
  TYPE (t_ncatt_array)        , POINTER :: namatts(:) ! list of attributes
  INTEGER                     , POINTER :: RGT(:)    ! regridding type
  CHARACTER (LEN=3)           , POINTER :: RGTstr(:) ! regridding type string
  INTEGER                     , POINTER :: RGrange(:,:) ! range for ixf regridding

  ! LOCAL
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: zivar(:)   ! input variable names
  CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: zovar(:)   ! output variable names
  REAL(dp)                    , POINTER :: zscl(:)    ! scaling
  TYPE (t_ncatt_array)        , POINTER :: znamatts(:) ! list of attributes
  INTEGER                     , POINTER :: zRGT(:)    ! regridding type
  CHARACTER (LEN=3)           , POINTER :: zRGTstr(:) ! regridding type string
  INTEGER                     , POINTER :: zRGrange(:,:) ! range for ixf regridding
  !
  INTEGER :: nv, i
  INTEGER :: pe_inc, active_proc
  INTEGER :: n, no
  INTEGER :: ix,jx ! ub_ak_20190724
  INTEGER, DIMENSION(:), ALLOCATABLE :: znv

  ! SAVE INPUT
  nv = SIZE(ivar)
  ALLOCATE(zivar(nv))
  ALLOCATE(zovar(nv))
  ALLOCATE(zscl(nv))
  ALLOCATE(znamatts(nv)) ! ub_ak_20190724
  ALLOCATE(zRGT(nv))
  ALLOCATE(zRGTstr(nv))
  ALLOCATE(zRGrange(nv,2))
  !
  zivar(:) = ivar(:)
  zovar(:) = ovar(:)
  zscl(:)  = scl(:)
  zRGrange(:,1:2) = RGrange(:,1:2)

  ! ub_ak_20190724+
  DO ix = 1, SIZE(namatts)
     znamatts(ix)%natts = namatts(ix)%natts
     ALLOCATE(znamatts(ix)%atts( namatts(ix)%natts))
     DO jx = 1, namatts(ix)%natts
        CALL COPY_NCATT(znamatts(ix)%atts(jx), namatts(ix)%atts(jx))
     END DO
  END DO
  ! ub_ak_20190724-
  zRGT(:)   = RGT(:)
  zRGTstr(:) = RGTstr(:)
  
  ! RE-SET OUTPUT
  DEALLOCATE(ivar);   NULLIFY(ivar)
  DEALLOCATE(ovar);   NULLIFY(ovar)
  DEALLOCATE(scl);    NULLIFY(scl)
  DEALLOCATE(RGT);    NULLIFY(RGT)
  DEALLOCATE(RGTstr); NULLIFY(RGTstr)
  DEALLOCATE(RGrange); NULLIFY(RGrange)
  ! ub_ak_20190724+
  IF (ASSOCIATED(namatts)) THEN
     DO ix = 1, SIZE(namatts)
        DO jx = 1, namatts(ix)%natts
           CALL INIT_NCATT(namatts(ix)%atts(jx))
        END DO
        IF (ASSOCIATED(namatts(ix)%atts)) THEN
           DEALLOCATE(namatts(ix)%atts); NULLIFY(namatts(ix)%atts)
           namatts(ix)%natts = 0 ! initialize
        END IF
     END DO
     DEALLOCATE(namatts)
  END IF
  NULLIFY(namatts)
  ! ub_ak_20190724-

  ! determine working PEs
  i_am_worker = .FALSE.
  nproc_work = MIN(nproc, nv, MAX_working_PEs)
  nproc_work = MAX(nproc_work,1)
  IF (ALLOCATED(pe_list)) DEALLOCATE(pe_list)
  ALLOCATE(pe_list(0:nproc_work-1))
  pe_inc = nproc/nproc_work
  DO i=0,nproc_work-1
     pe_list(i) = 0+i*pe_inc
  END DO

  ! distribute variables across working PEs
  ALLOCATE(znv(0:nproc_work-1))
  znv(:) = 0
  active_proc = 0
  DO i=1, nv
     znv(active_proc) = znv(active_proc) + 1
     i_am_worker = i_am_worker .OR. (pe_list(active_proc) == me)
     active_proc = active_proc+1
     active_proc = MOD(active_proc,nproc_work)
  END DO

  no = 0
  DO i=0, nproc_work-1
     n = znv(i)
     IF (pe_list(i) == me) THEN
        ALLOCATE(ivar(n))
        ivar(:) = zivar(no+1:no+n)
        ALLOCATE(ovar(n))
        ovar(:) = zovar(no+1:no+n)
        ALLOCATE(scl(n))
        scl(:) = zscl(no+1:no+n)
        ALLOCATE(namatts(n))
        DO ix = 1,n
           namatts(ix)%natts = znamatts(no+ix)%natts
           DO jx = 1, namatts(ix)%natts
              CALL COPY_NCATT( namatts(ix)%atts(jx),znamatts(no+ix)%atts(jx))
          END DO
        END DO
        ALLOCATE(RGT(n))
        RGT(:) = zRGT(no+1:no+n)
        ALLOCATE(RGTstr(n))
        RGTstr(:) = zRGTstr(no+1:no+n)
        ALLOCATE(RGrange(n,2))
        RGrange(:,1:2) = zRGrange(no+1:no+n,1:2)
     END IF
     no = no + n
  END DO

  ! CLEAN MEMORY
  DEALLOCATE(zivar);   NULLIFY(zivar)
  DEALLOCATE(zovar);   NULLIFY(zovar)
  DEALLOCATE(zscl);    NULLIFY(zscl)
  DEALLOCATE(zRGT);    NULLIFY(zRGT)
  DEALLOCATE(zRGTstr); NULLIFY(zRGTstr)
  DEALLOCATE(zRGrange); NULLIFY(zRGrange)
  ! ub_ak_20190724+
  IF (ASSOCIATED(znamatts)) THEN
     DO ix = 1, SIZE(znamatts)
        DO jx = 1, znamatts(ix)%natts
           CALL INIT_NCATT(znamatts(ix)%atts(jx))
        END DO
        DEALLOCATE(znamatts(ix)%atts); NULLIFY(znamatts(ix)%atts)
     END DO
     DEALLOCATE(znamatts)
  END IF
  NULLIFY(znamatts)
  ! ub_ak_20190724-

END SUBROUTINE DISTRIBUTE_VARS_ON_PES
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPAND_FILENAME(outfile, me)

  IMPLICIT NONE

  ! I/O
  CHARACTER(len=*),INTENT(INOUT)  :: outfile  
  INTEGER,         INTENT(IN)     :: me

  ! LOCAL
  CHARACTER(len=4)                :: cpe

  IF(me /= -1 .AND. outfile /= "") THEN 
     WRITE(cpe,'(i3.3,''_'')') me
     outfile = cpe//TRIM(outfile)
  END IF

END SUBROUTINE EXPAND_FILENAME
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_PIMPGRID

  IMPLICIT NONE

  IF (ALLOCATED(pe_list)) DEALLOCATE(pe_list)
  ALLOCATE(pe_list(0:0))
  pe_list(0) = 0

END SUBROUTINE INIT_PIMPGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_PARALLEL(p_pe, p_nprocs, p_all_comm)

  IMPLICIT NONE

  ! I/O
  INTEGER, INTENT(IN) :: p_pe
  INTEGER, INTENT(IN) :: p_nprocs
  INTEGER, INTENT(IN) :: p_all_comm

  my_mpi%rank  = p_pe
  my_mpi%nproc = p_nprocs
  my_mpi%comm  = p_all_comm

END SUBROUTINE INIT_PARALLEL
! ------------------------------------------------------------------

END MODULE messy_main_import_grid_par
