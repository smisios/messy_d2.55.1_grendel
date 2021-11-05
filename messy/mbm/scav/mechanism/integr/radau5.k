#INTEGRATOR radau5

#INLINE F95_DECL
  ! from radau5.k:
  LOGICAL,  PARAMETER :: logsteps           = .FALSE.
  INTEGER,  PARAMETER :: nsubsteps          = 1
  REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0. /)
#ENDINLINE
