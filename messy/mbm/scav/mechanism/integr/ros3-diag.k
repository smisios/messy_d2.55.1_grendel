{ros3 integrator with diagnostic output}
#INTEGRATOR ./ros3-diag

#INLINE F95_DECL
  ! from ros3-diag.k
  LOGICAL,  PARAMETER :: logsteps           = .FALSE.
  INTEGER,  PARAMETER :: nsubsteps          = 1
  REAL(dp), PARAMETER :: substep(nsubsteps) = (/ 0. /)
  !mz_rs_20031025+ diagnostic info
  INTEGER Nfcn, Njac, Naccept, Nreject, Nier, Ndecr, Nsame, Nincr
  REAL(dp) :: Hmin, Hstart, Hlow, Hmax
  INTEGER :: Nros3call=0
  !mz_rs_20031025-

#ENDINLINE

#INLINE F95_INIT
        ! from ros3-diag.k
        STEPMIN=0.0001
        STEPSTART=100.*STEPMIN
        STEPMAX=4800.
#ENDINLINE
