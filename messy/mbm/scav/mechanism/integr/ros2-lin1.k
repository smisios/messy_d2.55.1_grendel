#INTEGRATOR ros2

#INLINE F95_DECL
  ! from ros2-lin.k:
  LOGICAL,  PARAMETER :: logsteps           = .FALSE.
  INTEGER,  PARAMETER :: nsubsteps          = 1
  REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0. /)
#ENDINLINE
