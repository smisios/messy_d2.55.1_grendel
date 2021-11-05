# -*- Shell-script -*-

# Author: Patrick Joeckel, 2020
# - based on CCMI-base-02.bat from Patrick Joeckel (2013)
# - create SCAV-L-mechanism with: (Sc && !I && !Hg)
# - H tagging and D/H isotope-inclusive setups (Winterstein, Gromov)
# - corrections to kinetics for 3D runs (Winterstein, Gromov)
#
# The shell variables defined here will be used by xmecca 
# when it is run in batch mode (i.e. not interactive).

 set apn          = 0            # number of aerosol phases [0...99, default=0]
 set gaseqnfile   = gas.eqn
 ### back to MIM1, then other replacements
 set rplfile      = mim1-CCMI2-base-01-ISO-nmo 
 set ignoremassbalance                # MIM1 violates the mass balance
 set wanted       = "(((Tr && (G || Het) && \!I) || St) && \!Hg)"
 set enthalpy     = n                 # activate enthalpy in kJ/mol?
 set mcfct        = n                 # Monte-Carlo factor?
 set diagtracfile = CCMI2-base-01-ISO-nmo # diagnostic tracers?
 set rxnrates     = n                 # calculate accumulated reaction rates?
 set tag          = y                 # perform (isotope) tagging?
  set tagcfg       = j
#  set embud        = H
 ### chemically fixed species
#  set setfixlist   = "O2; O3; N2; N2O; NO2; CO2; CH4; Cl; Br;"    
 set kppoption    = 4                 # k=kpp, 4=kp4, q=quit
 set integr       = rosenbrock_mz     # integrator
 set vlen         = 256               # only for kp4 and integr=rosenbrock_vec
 set decomp       = 1                 # remove indirect indexing
                                      # kp4: 0/1/2/3/q; kpp: y/n/q
 set deltmpkp4    = y                 # delete temporary kp4 files?
 set latex        = y                 # latex list of reactions?
 set graphviz     = y                 # graphviz plots?
 set deltmp       = y                 # delete temporary xmecca files?
