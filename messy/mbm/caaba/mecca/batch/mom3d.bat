# -*- Shell-script -*-

# MOM for 3D model runs (with oracle and tagging)
# Author: Andrea Pozzer (2017)

# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      = mom3d             # few changes, mainly sulfur chemistry
 set ignoremassbalance                 
 set ignorecarboncount                # ORACLE miss tracer specifics
 set wanted       = "(((Tr && (G || Het) && \!I) || St) && \!Hg)"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?

 set tag          = y                 # perform (isotope) tagging?
 set tagcfg       = E                 # empty configurations for tagging
 set embud        = p                 # extended budgeting : OH reactivity calculation

 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_posdef # integrator
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set latex        = n                 # latex list of reactions?
 set graphviz     = n                 # graphviz plots?
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set deltmp       = y                 # delete temporary xmecca files?
