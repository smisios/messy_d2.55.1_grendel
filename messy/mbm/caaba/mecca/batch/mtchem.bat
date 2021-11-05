# -*- Shell-script -*-

# Author: Andreas Baumgaertner (2011)
# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

# Activate submodel="mtchem" to create an additional submodel that can be
# used in parallel with the MECCA submodel for 3D model runs. However,
# deactivate the line for creating the box model chemistry mechanism:
 set submodel     = "mtchem"          # submodel name

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      =                   # no replacements
 set wanted       = "Up"
 set enthalpy     = y                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = n                 # perform (isotope) tagging?
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
#set integr       = rosenbrock_vec    # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set latex        = y                 # latex list of reactions?
 set graphviz     = n                 # graphviz plots?
 set deltmp       = y                 # delete temporary xmecca files?
