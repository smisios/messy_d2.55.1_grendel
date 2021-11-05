#INTEGRATOR ros2

#INLINE F95_DECL
  ! from ros2-log3.k:
  LOGICAL,  PARAMETER :: logsteps           = .TRUE.
  INTEGER,  PARAMETER :: nsubsteps          = 3
  REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0.04 ,0.35 ,0.61 /)
#ENDINLINE
