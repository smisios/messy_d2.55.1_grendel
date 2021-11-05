#INTEGRATOR ros2

#INLINE F95_DECL
  ! from ros2-log6.k:
  LOGICAL,  PARAMETER :: logsteps           = .TRUE.
  INTEGER,  PARAMETER :: nsubsteps          = 6
  REAL(dp), PARAMETER :: substep(nsubsteps) = &
    (/ 0.01, 0.04, 0.12, 0.2, 0.29, 0.34 /)
#ENDINLINE

