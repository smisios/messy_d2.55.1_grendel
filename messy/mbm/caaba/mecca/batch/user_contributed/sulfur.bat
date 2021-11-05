# -*- Shell-script -*-

# create detailed sulfur chemistry using sulfur.rpl from J. Schallock

# Author: Rolf Sander (2016)

 set apn          = 0                 # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 set rplfile      = sulfur            # full sulfur chemistry
 set wanted       = "(((Tr && (G || Het) && \!I && \!C) || St) && \!Hg)"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile =                   # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = n                 # perform (isotope) tagging?
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set latex        = y                 # latex list of reactions?
 set graphviz     = y                 # graphviz plots?
 set deltmp       = y                 # delete temporary xmecca files?

##############################################################################
echo "THE FILE sulfur.rpl IS UNDER CONSTRUCTION! DO NOT USE IT (YET)!"
exit 1
##############################################################################
