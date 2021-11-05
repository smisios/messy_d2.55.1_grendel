#INTEGRATOR ros2

#INLINE F95_DECL
  ! from ros2-lin9.k:
  LOGICAL,  PARAMETER :: logsteps           = .FALSE.
  INTEGER,  PARAMETER :: nsubsteps          = 9
  REAL(dp), PARAMETER :: substep(nsubsteps) = &
    (/ 0., 0., 0., 0., 0., 0., 0., 0., 0. /)
#ENDINLINE
