
MODULE mo_tr_allgather

  ! Patrick Joeckel, DLR, August 2016

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc => local_decomposition,  &
                              gdc => global_decomposition, &
                              debug_parallel
  USE mo_mpi,           ONLY: p_nprocs, p_allgather
  USE mo_transpose,     ONLY: reorder

  IMPLICIT NONE

  PRIVATE
  SAVE

#ifdef _USE_MGATHER
  INCLUDE 'mpif.h'
#endif

  INTERFACE allgather_field
     MODULE PROCEDURE allgather_gp432
     MODULE PROCEDURE allgather_gp32
     MODULE PROCEDURE allgather_gp2
  END INTERFACE allgather_field

  PUBLIC :: allgather_field

  ! module state:
  LOGICAL :: require_init = .TRUE.

  ! field dimensions
  INTEGER :: nlon_max = 0
  INTEGER :: nlat_max = 0
  INTEGER :: nlon = 0
  INTEGER :: nlat = 0

CONTAINS

  SUBROUTINE init_tr_allgather
    CHARACTER(len=*), PARAMETER :: &
         context = 'mo_tr_allgather::init_tr_allgather'
    INTEGER :: i

    IF (.NOT. require_init) RETURN
    require_init = .FALSE.
  
    DO i=1, p_nprocs
       IF (gdc(i)%nglon > nlon_max) nlon_max = gdc(i)%nglon
       IF (gdc(i)%nglat > nlat_max) nlat_max = gdc(i)%nglat
    END DO

    nlon = gdc(1)%nlon
    nlat = gdc(1)%nlat

write(*,*) 'QQQ ',context,'',nlon_max,nlat_max,nlon,nlat
    
  END SUBROUTINE init_tr_allgather
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE allgather_gp432(global, local)

    ! nlon x nlev x nt x nlat
    ! nlon x nlev x nlat x 1

    REAL(dp), POINTER                :: global(:,:,:,:)
    REAL(dp), TARGET,  INTENT(in)    :: local(:,:,:,:)

    REAL(dp), POINTER :: global3d(:,:,:)
    REAL(dp), POINTER :: local3d(:,:,:)
    
    INTEGER :: d4, d3, d2, n

    IF (require_init) CALL init_tr_allgather

    d4 = SIZE(global,4)
    d2 = SIZE(global,2)

    IF (d4 == 1) THEN
       NULLIFY(global3d)
       global3d => global(:,:,:,1)
       CALL allgather_gp32 (global3d, local(:,:,:,1))
    ELSE
       d3 = SIZE(global,3)
       NULLIFY(global3d) ; ALLOCATE(global3d(nlon,d2,nlat))
       NULLIFY(local3d)
       ALLOCATE(local3d(SIZE(local,1),SIZE(local,2),SIZE(local,3)))
       DO n=1, d3
          local3d(:,:,:) = local(:,:,n,:)
          CALL allgather_gp32(global3d, local3d)
          global(:,:,n,:) = global3d(:,:,:)
       END DO
       DEALLOCATE(global3d) ; NULLIFY(global3d)
       DEALLOCATE(local3d) ; NULLIFY(local3d)
    END IF

  END SUBROUTINE allgather_gp432
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE allgather_gp32(global, local)

    ! nlon x nlev x nlat
    ! nlon x nlat x 1

    REAL(dp), POINTER                :: global(:,:,:)
    REAL(dp), TARGET,  INTENT(in)    :: local(:,:,:)

    REAL(dp), POINTER :: global2d(:,:)
    REAL(dp), POINTER :: local2d(:,:)
    
    INTEGER :: d2, d3, n

    IF (require_init) CALL init_tr_allgather

    d3 = SIZE(global,3)
    
    IF (d3 == 1) THEN
       NULLIFY(global2d)
       global2d => global(:,:,1)
       CALL allgather_gp2(global2d, local(:,:,1))
    ELSE
       d2 = SIZE(global,2)
       NULLIFY(global2d) ; ALLOCATE(global2d(nlon,nlat))
       NULLIFY(local2d) ; ALLOCATE(local2d(SIZE(local,1),SIZE(local,3)))
       DO n=1, d2
          local2d(:,:) = local(:,n,:)
          CALL allgather_gp2(global2d, local2d)
          global(:,n,:) = global2d(:,:)
       END DO
       DEALLOCATE(global2d) ; NULLIFY(global2d)
       DEALLOCATE(local2d) ; NULLIFY(local2d)
    ENDIF

  END SUBROUTINE allgather_gp32
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE allgather_gp2(global, local)

    ! nlon x nlat

    REAL(dp), POINTER                :: global(:,:)
    REAL(dp), TARGET,  INTENT(in)    :: local(:,:)

    REAL(dp), POINTER, DIMENSION(:,:)   :: sendbuf
    REAL(dp), POINTER, DIMENSION(:,:,:) :: recvbuf

    INTEGER :: i, pe, nlon
    LOGICAL :: lreg
    INTEGER :: gx1, gx2, gy1, gy2
    INTEGER :: lx1, lx2, ly1, ly2

    IF (require_init) CALL init_tr_allgather

    NULLIFY(recvbuf)
    ALLOCATE(recvbuf(nlon_max,nlat_max,0:p_nprocs-1))
    recvbuf(:,:,:) = 0.0_dp

    nlon = ldc%nlon
    lreg = ldc%lreg

    ALLOCATE(sendbuf(nlon_max,nlat_max))
    sendbuf(:,:) = 0.0_dp

    IF (lreg) THEN
       sendbuf(:,:) = local(:,:)
    ELSE
       CALL reorder(sendbuf(1:ldc%nglon,1:ldc%nglat),local)
    END IF

!write(*,*) 'QQQ calling p_allgather ...',lreg,UBOUND(sendbuf),UBOUND(recvbuf)
    CALL p_allgather(recvbuf, sendbuf)
!write(*,*) 'QQQ ... called  p_allgather',UBOUND(global)

    DO i = lbound(gdc,1), ubound(gdc,1)
!!$    DO i = 1, p_nprocs
       pe    = gdc(i)%pe
       
       ! unpack first segment
       gx1 = gdc(i)%glons(1)
       gx2 = gdc(i)%glone(1)
       gy1 = gdc(i)%glats(1)
       gy2 = gdc(i)%glate(1)

       lx1 = LBOUND(recvbuf,1)
       lx2 = gdc(i)%nglon ! UBOUND(recvbuf,1)
       ly1 = LBOUND(recvbuf,2)
       ly2 = gdc(i)%nglh(1)

!write(*,*) 'QQQa ',i,pe,'[',gx1,':',gx2,',',gy1,':',gy2,'] <<-- [',&
!                            lx1,':',lx2,',',ly1,':',ly2,':',pe,']'

       global(gx1:gx2, gy1:gy2) = recvbuf(lx1:lx2, ly1:ly2, pe)

       ! unpack second segment
       IF (gdc(i)%nglh(2) > 0) THEN
          IF (gdc(i)%glone(2) > gdc(i)%glons(2)) THEN

             gx1 = gdc(i)%glons(2)
             gx2 = gdc(i)%glone(2)
             gy1 = gdc(i)%glats(2)
             gy2 = gdc(i)%glate(2)

             lx1 = LBOUND(recvbuf,1)
             lx2 = gdc(i)%nglon ! UBOUND(recvbuf,1)
             ly1 = gdc(i)%nglat/2+1
             ly2 = gdc(i)%nglat ! UBOUND(recvbuf,2)

!write(*,*) 'QQQb ',i,pe,'[',gx1,':',gx2,',',gy1,':',gy2,'] <<-- [',&
!                            lx1,':',lx2,',',ly1,':',ly2,':',pe,']'

             global(gx1:gx2, gy1:gy2) = recvbuf(lx1:lx2, ly1:ly2, pe)

          ELSE
             ! unpacking second segment, split in longitudes

             gx1 = gdc(i)%glons(2)
             gx2 = nlon
             gy1 = gdc(i)%glats(2)
             gy2 = gdc(i)%glate(2)

             lx1 = LBOUND(recvbuf,1)
             lx2 = nlon-gdc(i)%glons(2)+1
             ly1 = gdc(i)%nglat/2+1
             ly2 = gdc(i)%nglat ! UBOUND(recvbuf,2)             

!write(*,*) 'QQQc ',i,pe,'[',gx1,':',gx2,',',gy1,':',gy2,'] <<-- [',&
!                            lx1,':',lx2,',',ly1,':',ly2,':',pe,']'

             global(gx1:gx2, gy1:gy2) = recvbuf(lx1:lx2, ly1:ly2, pe)

             gx1 = LBOUND(global,1)
             gx2 = gdc(i)%glone(2)
             gy1 = gdc(i)%glats(2)
             gy2 = gdc(i)%glate(2)

             lx1 = gdc(i)%nglon-gdc(i)%glone(2)+1
             lx2 = nlon ! UBOUND(recvbuf,1)
             ly1 = gdc(i)%nglat/2+1
             ly2 = gdc(i)%nglat ! UBOUND(recvbuf,2)      

!write(*,*) 'QQQd ',i,pe,'[',gx1,':',gx2,',',gy1,':',gy2,'] <<-- [',&
!                            lx1,':',lx2,',',ly1,':',ly2,':',pe,']'

             global(gx1:gx2, gy1:gy2) = recvbuf(lx1:lx2, ly1:ly2, pe)

          ENDIF
       ENDIF
    END DO

    DEALLOCATE(sendbuf)
    NULLIFY(sendbuf)
    
    DEALLOCATE(recvbuf)
    NULLIFY(recvbuf)

  END SUBROUTINE allgather_gp2
!-----------------------------------------------------------------------------

END MODULE mo_tr_allgather
