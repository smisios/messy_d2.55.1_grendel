#INTEGRATOR ros2

#INLINE F95_DECL
  ! from ros2-log10.k:
  LOGICAL,  PARAMETER :: logsteps           = .TRUE.
  INTEGER,  PARAMETER :: nsubsteps          = 10
  REAL(dp), PARAMETER :: substep(nsubsteps) = &
    (/ 0.0005, 0.0015, 0.005, 0.02, 0.06, 0.1, &
       0.15, 0.1875, 0.23775, 0.23775 /)
#ENDINLINE

