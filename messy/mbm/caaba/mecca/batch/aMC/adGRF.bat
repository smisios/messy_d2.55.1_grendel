# -*- Shell-script -*-

# Sergey Gromov (2019)
# atmospheric dynamics & global radioactive fallout

 set gaseqnfile   = eqn/tag/adGRF.eqn
 set gasspcfile   = eqn/tag/adGRF.spc
#set gastblfile   = gas.tbl

# features:
 set ignoremassbalance

 set rplfile      =                   #
 set wanted       = "\!Ara"
 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?

 set tag          = n                 # perform (isotope) tagging?
#set tagcfg       =                   # tagging cfg(s)
#set embud        = c                 # extended budgeting

 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
#set kppoption    = k                 # k=kpp, 4=kp4, q=quit

# integrator
#set integr       = rosenbrock_posdef_h211b_qssa
 set integr       = rosenbrock_mz
#set integr       = rosenbrock_vec
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files
 set latex        = n                 # latex list of reactions
 set deltmptex    = y                 # delete temporary LaTeX files?
 set deltmp       = y                 # delete temporary xmecca files?
 set graphviz     = n                 # graphviz plots

 set setfixlist   = 'H2O;'            # H2O is required for CAABA functioning
