! **********************************************************************
MODULE messy_submod3
! **********************************************************************

  ! SMCL
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'submod3'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'

  PUBLIC :: dp
  PUBLIC :: osm_fill
  PUBLIC :: osm_diff
  PUBLIC :: osm_decay

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE osm_fill(xt, scale)

    IMPLICIT NONE
    INTRINSIC :: REAL, SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: xt    ! tracer data
    REAL(DP),                   INTENT(IN)    :: scale ! arbitrary scale
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr = 'osm_fill'
    INTEGER                             :: d1,d2,d3  ! dimensions
    INTEGER                             :: j1,j2,j3  ! counter
    REAL(dp)                            :: rd1, rd2, rd3
    REAL(dp)                            :: rj1, rj2, rj3

    ! FILL TRACER WITH SOME ARBITRARY VALUES
    d1 = SIZE(xt,1)
    d2 = SIZE(xt,2)
    d3 = SIZE(xt,3)

    rd1 = REAL(d1,dp)
    rd2 = REAL(d2,dp)
    rd3 = REAL(d3,dp)

    DO j1=1,d1
       rj1 = REAL(j1,dp)
       DO j2=1,d2
          rj2 = REAL(j2,dp)
          DO j3=1,d3
             rj3 = REAL(j3,dp)

             xt(j1,j2,j3) = ( &
                  (rj1 - rd1*0.5_dp)**2 + (rj2 - rd2*0.5_dp)**2 &
                  )**0.5_dp * rj3/rd3 * scale
          END DO
       END DO
    END DO

  END SUBROUTINE osm_fill
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE osm_diff(xt)

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: xt    ! tracer data
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr = 'osm_diff'
    INTEGER                             :: d1,d2,d3     ! dimensions
    INTEGER                             :: j1,j2,j3     ! counter
    INTEGER                             :: j1p, j1m
    INTEGER                             :: j2p, j2m
    INTEGER                             :: j3p, j3m
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: tmp

    ! 'DIFFUSE' TRACERS
    d1 = SIZE(xt,1)
    d2 = SIZE(xt,2)
    d3 = SIZE(xt,3)

    ALLOCATE(tmp(d1,d2,d3))

    DO j1=1, d1
       j1m = j1-1
       j1p = j1+1
       IF (j1p > d1) j1p = 1
       IF (j1m <  1) j1m = d1

       DO j2=1, d2
          j2m = j2-1
          j2p = j2+1
          IF (j2p > d2) j2p = 1
          IF (j2m <  1) j2m = d2

          DO j3=1, d3
             j3m = j3-1
             j3p = j3+1
             IF (j3p > d3) j3p = 1
             IF (j3m <  1) j3m = d3

             tmp(j1,j2,j3) = ( &
                    xt(j1p,j2,j3) + xt(j1m,j2,j3) &
                  + xt(j1,j2p,j3) + xt(j1,j2m,j3) &
                  + xt(j1,j2,j3p) + xt(j1,j2,j3p) &
                  ) / 6.0_dp
          END DO
       END DO
    END DO

    xt(:,:,:) = tmp(:,:,:)
    DEALLOCATE(tmp)

  END SUBROUTINE osm_diff
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE osm_decay(qxt, decay)

    IMPLICIT NONE
    INTRINSIC :: REAL, SIZE

    ! I/O
    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: qxt   ! tracer data
    REAL(dp),                 INTENT(IN)    :: decay

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr = 'osm_decay'

    ! SIMULATE SIMPLE DECAY
    IF (decay <= 0.0_dp) RETURN
    qxt(:,:) = qxt(:,:) / decay

  END SUBROUTINE osm_decay
  ! --------------------------------------------------------------------

! **********************************************************************
END MODULE messy_submod3
! **********************************************************************
