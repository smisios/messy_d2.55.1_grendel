# -*- Shell-script -*-

# Author: Patrick Joeckel (2008)
# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

# Mod: Sergey Gromov (2014)
# exb_example - an example of using the extended budgeting tools for MECCA

 set gaseqnfile   = gas.eqn
 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set rplfile      =                   #
 set wanted       = "((Tr && (G || Het) && \!Cl && \!Br && \!I) || St) && \!Hg"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = y                 # perform (isotope) tagging?
 set tagcfg       = E                 # empty cfg
 set embud        = e                 # ext. budgeting
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
#set kppoption    = k                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
#set integr       = rosenbrock_vec    # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set latex        = n                 # latex list of reactions?
 set graphviz     = n                 # graphviz plots?
 set deltmptex    = y                 # delete temporary LaTeX files?
 set deltmp       = y                 # delete temporary xmecca files?
