# -*- Shell-script -*-

# Author: Sergey Gromov [2015]
# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

 ln -i -s messy_mecca_kpp-PM_tans.kpp messy_mecca_kpp.kpp

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas-PM_tans.eqn
 set gasspcfile   = gas-PM_tans.spc
 set gastexfile   = gas-PM_tans.tex
 set rplfile      =                   # no replacements
 set wanted       = "\!Ara"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = y                 # perform (isotope) tagging?
 set tagcfg       = chr               # isocarbon, isohydrogen, radiocarbon
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_posdef # integrator
 set decomp       = n                 # remove indirect indexing
 set latex        = n                 # latex list of reactions
 set graphviz     = n                 # graphviz plots?
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set deltmp       = y                 # delete temporary xmecca files?
