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
 set tagcfg       = aA                # tagging cfg(s)
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

 set setfixlist   = 'Dummy; \
  HL; I12HL; I13HL; I16HL; I17HL; I18HL; \
  LL; I12LL; I13LL; I16LL; I17LL; I18LL; \
  RC; I12RC; I13RC; I16RC; I17RC; I18RC; \
  CM; I12CM; I13CM; I16CM; I17CM; I18CM; \
  HL_CH4; HL_N2O; HL_SF6; HL_CO; HL_O3; \
  LL_CH4; LL_N2O; LL_SF6; LL_CO; LL_O3; \
  BL_CH4; BL_N2O; BL_SF6; BL_CO; BL_O3; \
  CM_CH4; CM_N2O; CM_SF6; CM_CO; CM_O3; '

# xx; I12xx; I13xx; I16xx; I17xx; I18xx;
# HL_xx; LL_xx; BL_xx; 
# PC; I16RC; I17RC; I18RC;
