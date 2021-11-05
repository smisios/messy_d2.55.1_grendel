# -*- Shell-script -*-

# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

#set apn          = 2                 # number of aerosol phases [0...99, default=0]
set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      =                   # no replacements
#set wanted        = "Tr && (G || Aa) && \!S && \!Hg && \!Cl && \!I && \!Br"  # tropo, except..., + aerosol
#set wanted        = "Tr && G && \!S && \!Hg && \!Cl && \!I && \!Br && \!C"  # tropo, except ... (no NMHCs!)
set wanted        = "Tr && G && \!S && \!Hg && \!Cl && \!I && \!Br"  # tropo, except ... (since no input or rxn in tropo)
#set wanted       = "\!Ara && \!I && \!Br"  # all except rain chemistry, I, Br
#set wanted       = "\!Ara"           # all except rain chemistry
#set wanted       = "Tr && (G || (Aa && Mbl)) && \!I && \!Hg"
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers? don't put anything if 'no'
 set rxnrates     = y                 # calculate accumulated reaction rates?
 set tag          = n                 # perform (isotope) tagging?
 set kppoption    = k                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_posdef # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = n                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set latex        = y                 # latex list of reactions?
 set graphviz     = y                 # graphviz plots?
 set deltmp       = n                 # delete temporary xmecca files?
