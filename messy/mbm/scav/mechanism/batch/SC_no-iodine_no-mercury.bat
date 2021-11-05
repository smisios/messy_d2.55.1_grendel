# -*- Shell-script -*-

# Author: Rolf Sander (2014),
# - based on CCMI-base-01.bat from Patrick Joeckel (2013)
# - create SCAV-L-mechanism with: (Sc && !I && !Hg)
# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

 set apn          = 1                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set ignoremassbalance                # MIM1 violates the mass balance
 set wanted       = "(Sc && \!I && \!Hg)"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                  # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = n                 # perform (isotope) tagging?
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
 set vlen         = 16                # only for kp4 and integr=rosenbrock_vec
 set decomp       = 0                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set latex        = y                 # latex list of reactions?
 set graphviz     = n                 # graphviz plots?
 set deltmp       = y                 # delete temporary xmecca files?
