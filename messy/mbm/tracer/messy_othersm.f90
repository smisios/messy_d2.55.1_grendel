! **********************************************************************
MODULE messy_othersm
! **********************************************************************

  ! SMCL
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'othersm'
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

    ! FILL TRACER WITH SOME ARBITRARY VALUES
    d1 = SIZE(xt,1)
    d2 = SIZE(xt,2)
    d3 = SIZE(xt,3)

    DO j1=1,d1
       DO j2=1,d2
          DO j3=1,d3
             xt(j1,j2,j3) = REAL(j1*(d3-j3+1), dp) &
                  * REAl(j3,dp)/REAL(d3,dp) * scale
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
    INTEGER                             :: jj1,jj2,jj3  ! counter

    ! 'DIFFUSE' TRACERS
    d1 = SIZE(xt,1)
    d2 = SIZE(xt,2)
    d3 = SIZE(xt,3)

    DO j1=1, d1
       DO j2=1, d2
          DO j3=1, d3
             jj1 = j1+1
             IF (jj1 > d1) jj1 = 1
             jj2 = j2+1
             IF (jj2 > d2) jj1 = 1
             jj3 = j3+1
             IF (jj3 > d3) jj3 = 1
             xt(j1,j2,j3) = ( &
                  xt(j1,j2,j3)  + xt(jj1,j2,j3) + &
                  xt(j1,jj2,j3) + xt(j1,j2,jj3) ) / 4.0_dp
          END DO
       END DO
    END DO

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
END MODULE messy_othersm
! **********************************************************************
