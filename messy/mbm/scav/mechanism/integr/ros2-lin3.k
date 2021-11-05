#INTEGRATOR ros2

#INLINE F95_DECL
  ! from ros2-lin3.k:
  LOGICAL,  PARAMETER :: logsteps           = .FALSE.
  INTEGER,  PARAMETER :: nsubsteps          = 3
  REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0., 0., 0. /)
#ENDINLINE
