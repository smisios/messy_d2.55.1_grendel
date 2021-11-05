#INTEGRATOR ros2-auto

#INLINE F95_DECL
  ! from ros2-auto.k:
  LOGICAL,  PARAMETER :: logsteps           = .FALSE.
  INTEGER,  PARAMETER :: nsubsteps          = 1
  REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0. /)
#ENDINLINE
