# -*- Shell-script -*-

# Sergey Gromov (2015-2016)
# isotope CO2 Monte-Carlo box-model study for CARIBIC-2 data

# features:
 set ignoremassbalance
 set gaseqnfile   = gas.eqn
 set rplfile      = tag-clean
 set wanted       = "\!Ara"

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?

 set tag          = y                 # perform (isotope) tagging?
 set tagcfg       = bB                # tagging cfg(s)
#set embud        = c                 # extended budgeting

 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
#set kppoption    = k                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
#set integr       = rosenbrock_vec    # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files
 set latex        = n                 # latex list of reactions
 set deltmptex    = y                 # delete temporary LaTeX files?
 set deltmp       = y                 # delete temporary xmecca files?
 set graphviz     = n                 # graphviz plots

 set setfixlist   = 'Dummy; EMISTN; I1EMISTN; I2EMISTN; I12EMISTN; I13EMISTN; EMISTS; I1EMISTS; I2EMISTS; I12EMISTS; I13EMISTS;'
