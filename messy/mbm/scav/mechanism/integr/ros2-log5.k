#INTEGRATOR ros2

#INLINE F95_DECL
  ! from ros2-log5.k:
  LOGICAL,  PARAMETER :: logsteps           = .TRUE.
  INTEGER,  PARAMETER :: nsubsteps          = 5
  ! mz_cb_20030607+
  REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0.01 ,0.04, 0.16, 0.32, 0.47 /)
  !REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0.02 ,0.12 ,0.22 ,0.30, 0.34 /)
  ! mz_cb_20030607-
#ENDINLINE
