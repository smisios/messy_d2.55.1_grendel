# -*- Shell-script -*-

# author: S. Gromov [MPIC, 2020]
# setup: stratospheric isotopic composition of carbon monoxide (SICM) tagging

 set ignoremassbalance                # the reaction scheme "violates" element balance
 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      = mim1-PalMod_nmo
#eval:
 set wanted       = "(((Tr && (G || Het) && \!I) || St) && \!Hg)"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = y                 # perform (isotope) tagging?
 set tagcfg       = kmOsx             # m|O|s|k|x: methane carbon | O tagging | strat. CO | CO kie | CIO
 set embud        = px                # ext. budgeting form PalMod & CIO
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

 set setfixlist   = "N2;"

