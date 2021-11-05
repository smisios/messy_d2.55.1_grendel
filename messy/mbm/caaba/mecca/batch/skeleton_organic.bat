# -*- Shell-script -*-

# Author: Rolf Sander (2016 - 2018)
# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      = skeleton_organic  # reduce to skeletal mechanism
 set wanted       = "(((Tr && (G || Het) && \!I) || St) && \!Hg)" # from CCMI-base-02
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = n                 # perform (isotope) tagging?
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing kp4: 0/1/2/3/q; kpp: y/n/q
 set latex        = y                 # latex list of reactions?
 set graphviz     = y                 # graphviz plots?
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set deltmp       = y                 # delete temporary xmecca files?
