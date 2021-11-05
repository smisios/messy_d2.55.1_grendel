# -*- Shell-script -*-

# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

 set apn          = 2                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      =                   # no replacements
 set wanted       = "Tr && (G || (Aa && Mbl)) && \!I && \!Hg"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = n                 # perform (isotope) tagging?
 set setfixlist   = "O2; N2; N2O;"    # chemically fixed species
 set kppoption    = k                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_posdef # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = n                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set latex        = n                 # latex list of reactions?
 set graphviz     = n                 # graphviz plots?
 set deltmpkp4    = y                 # delete temporary kp4 files? (used in kp4.tcsh)
 set deltmp       = n                 # delete temporary xmecca files?
