# -*- Shell-script -*-

# author: S.Gromov (2015)
# setup: fully resolved clumped O isotopes based on Yeung et al. mechanism

 set ignoremassbalance
 set gaseqnfile   = gas.eqn
 set rplfile      = tag-clean
#eval:
 set wanted       = "\!Ara"           #
 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?

 set tag          = y                 # perform (isotope) tagging?
 set tagcfg       = x                 # tagging cfg
 set embud        = yx                # extended budgeting

#set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set kppoption    = k                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
#set integr       = rosenbrock_vec    # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
#set decomp       = y                 # remove indirect indexing
 set deltmpkp4    = y                 # delete temporary kp4 files
 set latex        = n                 # latex list of reactions
 set deltmptex    = y                 # delete temporary LaTeX files?
 set deltmp       = y                 # delete temporary xmecca files?
 set graphviz     = n                 # graphviz plots

 set setfixlist   = 'Dummy;'
