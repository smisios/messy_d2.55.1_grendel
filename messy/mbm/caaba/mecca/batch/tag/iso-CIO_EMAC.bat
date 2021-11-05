# -*- Shell-script -*-

# author: S.Gromov / 2015-2020
# clumped O2 isotope chemistry setup

 set ignoremassbalance

# earlier CH and TR setups
#set wanted       = "(((Tr && (G || Het) && \!I) || St) && \!Hg)"  # CCMI-like setup
# latest CH setup based on SICM
 set gaseqnfile   = gas.eqn
 set rplfile      = mim1-PalMod_nmo
 set wanted       = "(((Tr && (G || Het) && \!I) || St) && \!Hg)"
#
 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?

 set tag          = y                 # perform (isotope) tagging?
 set tagcfg       = xm#sk             # x:CIO | m:methane carbon | s:strat.CO | k:CO KIE
 set embud        = px                # p:PalMod | x:CIO

 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
#set kppoption    = k                 # k=kpp, 4=kp4, q=quit
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

 set setfixlist   = 'Dummy; N2;'      # CO2; H2;
