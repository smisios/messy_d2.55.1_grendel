# -*- Shell-script -*-

# Author: Christopher Kaiser
# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

#  mz_rs_20151118+
# WARNING: DMS reactions are now correctly labeled with "C", like all
# reactions that contain species with >= 2 C atoms. To ensure that these
# reactions are still included, the wanted string had to be adjusted to
# "S || \!C".
#  mz_rs_20151118-

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      = mim1-simple_MADE  # back to MIM1, then other replacements
 set ignoremassbalance                # MIM1 violates the mass balance
 set wanted       = "Tr && G && (S || \!C) && \!Cl && \!Br && \!I && \!Hg"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = n                 # perform (isotope) tagging?
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_posdef # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = y                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set latex        = y                 # latex list of reactions?
 set graphviz     = n                 # graphviz plots?
 set deltmp       = y                 # delete temporary xmecca files?
