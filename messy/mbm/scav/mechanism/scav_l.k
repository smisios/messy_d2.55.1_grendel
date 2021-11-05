{define chemical species}
#INCLUDE ./scav.spc

{select which species are fixed}
#SETFIX CO2; O2; N2; H2O_l; 

{define chemical equations}
#INCLUDE ./scav.eqn

{select programming language}
#USE FORTRAN95

{set initial values}
#INITVALUES
ALL_SPEC = 0.;

{choose an integrator and define substeps}
#INCLUDE ./integr.k

{the driver is just a dummy file}
#DRIVER ./driver
