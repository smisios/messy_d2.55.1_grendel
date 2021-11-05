# -*- Shell-script -*-

# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gasspcfile   = eqn/mcm/terpenes.spc
 set gaseqnfile   = eqn/mcm/terpenes.eqn
 set ignoremassbalance                # MCM violates the mass balance
 set setfixlist   = "O2;"
 set rplfile      =                   # no replacements
 set wanted       = "1"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = n                 # perform (isotope) tagging?
 set kppoption    = k                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_posdef # integrator
 set decomp       = n                 # remove indirect indexing
 set latex        = y                 # latex list of reactions?
 set graphviz     = y                 # graphviz plots?
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set deltmp       = n                 # delete temporary xmecca files?
