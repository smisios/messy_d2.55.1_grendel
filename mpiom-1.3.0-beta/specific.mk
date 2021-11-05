### THIS FILE CONTAINS SPECIFIC RULES AND IS INCLUDED IN THE MPIOM Makefile

#### NO  ARCHITECTURE SPECIFIC
#EXAMPLE: # NEC-SX6 / SUPER-UX
#ifeq ($(strip $(ARCH)), SX)
#mo_spitfire.o: mo_spitfire.f90
#	$(F90) $(F90ALL) -pi auto exp=minmod exp=medan exp=putyslice exp=cfdot1dp2 exp=cfint1x2 line=2000 -c $< -o $@
#mo_tpcore.o: mo_tpcore.f90
#	$(F90) $(F90ALL) -pi auto exp=xmist exp=fxppm exp=kmppm exp=lmppm exp=xtp noexp=map1_ppm_gp,ppm2m,steepz nest=3 line=1000 -c $< -o $@
#mo_transpose.o: mo_transpose.f90
#	$(F90) $(F90ALL) -Npi -c $< -o $@
#lti.o: lti.f90
#	$(F90) $(F90ALL) -Npi -c $< -o $@
#endif
############################################################################

# IBM-p690 / AIX #############################################################
ifeq ($(strip $(ARCH)), rs6000)
## NOTE: Version 13.01.0000.0008 obviously has problems with run-time checks,
##       when nudging is switched on ...
#mo_util_string.o: mo_util_string.f90
#	$(F90) $(F90DEFSMPIOM) $(F90NOR8) $(INCLUDES) -qnocheck -qnoflttrap -c $< -o $@
#mo_mpi.o: mo_mpi.f90
#	$(F90) $(F90DEFSMPIOM) $(F90NOR8) $(INCLUDES) -qnocheck -qnoflttrap -c $< -o $@
endif
# ############################################################################
