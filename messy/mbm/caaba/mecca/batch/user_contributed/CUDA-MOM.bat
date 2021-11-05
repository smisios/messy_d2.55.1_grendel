# -*- Shell-script -*-

# Author: Sergey Gromov (2018),
# - based on CCMI-base-02.bat from Patrick Joeckel (2013)

# batch file for MECCA/KPP-CUDA enabled intergation tests

# tests with MOM

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn

 set rplfile      =                   # no replacements
#set ignoremassbalance                # MIM1 violates the mass balance
 set wanted       = "(((Tr && (G || Het) && \!I) || St) && \!Hg)"

 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?

 set tag          = n                 # perform (isotope) tagging?
#set tagcfg       = E                 # empty configurations for tagging
#set embud        = p                 # extended budgeting : OH reactivity calculation

 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing kp4: 0/1/2/3/q; kpp: y/n/q

 set deltmpkp4    = y                 # delete temporary kp4 files?
 set latex        = n                 # latex list of reactions?
 set graphviz     = n                 # graphviz plots?
 set deltmp       = y                 # delete temporary xmecca files?
