# -*- Shell-script -*-

# Author: S.Gromov [2016]
# setup: DIUMOD - simulation of CH4 using average/modulated r-n rates

# The shell variables defined here will be used by xmecca
# when it is run in batch mode (i.e. not interactive).

 set ignoremassbalance                # the reaction scheme "violates" element balance
 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      = tag-clean
#eval:
 set wanted       = "\!Ara"           #
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = y                 # perform (isotope) tagging?
 set tagcfg       = rch               # radiocarbon, isocarbon, isohydrogen
 set embud        = d                 # produce DIUMOD chemistry
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
#set integr       = rosenbrock_posdef_h211b_qssa # integrator
#set integr       = rosenbrock_posdef # integrator
 set integr       = rosenbrock_mz     # integrator
#set integr       = rosenbrock_vec    # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files
 set latex        = n                 # latex list of reactions
 set graphviz     = n                 # graphviz plots
 set deltmp       = y                 # delete temporary xmecca files?

#set setfixlist   = "CO2; N2;"
 set setfixlist   = "CO2; O2; N2; O1Da; O3Pa; Cl; OHa; O1Dm; O3Pm; OHm; jCH4a; jCH4m;"
