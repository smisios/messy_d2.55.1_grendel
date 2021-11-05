#INTEGRATOR ros3-ice

#INLINE F95_DECL
  ! from ros3.k:
  LOGICAL,  PARAMETER :: logsteps           = .FALSE.
  INTEGER,  PARAMETER :: nsubsteps          = 1
  REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0. /)
#ENDINLINE

#INLINE F95_INIT
        ! from ros3.k:
        ! STEPMIN=0.0001
        STEPMIN=1.E-37 ! mz_ak_20040707
        STEPMAX=4800.
        STEPSTART=STEPMIN
#ENDINLINE


