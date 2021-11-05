### THIS FILE CONTAINS SPECIFIC RULES AND IS INCLUDED IN THE ECHAM5 Makefile
mo_mpi.o: mo_mpi.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<

### ECHAM5 - ARCHITECTURE SPECIFIC ###########################################

# NEC-SX6 / SUPER-UX #########################################################
ifeq ($(strip $(ARCH)), SX)
#mo_spitfire.o: mo_spitfire.f90
#	$(F90) $(F90ALL) -pi auto exp=minmod exp=medan exp=putyslice exp=cfdot1dp2 exp=cfint1x2 line=2000 -c $< -o $@
#mo_tpcore.o: mo_tpcore.f90
#	$(F90) $(F90ALL) -pi auto exp=xmist exp=fxppm exp=kmppm exp=lmppm exp=xtp noexp=map1_ppm_gp,ppm2m,steepz nest=3 line=1000 -c $< -o $@
mo_transpose.o: mo_transpose.f90
	$(F90) $(F90ALL) -Npi -c $< -o $@
lti.o: lti.f90
	$(F90) $(F90ALL) -Npi -c $< -o $@
mo_tr_gather.o: mo_tr_gather.f90
	$(F90) $(F90ALL) -f2003 nocbind -c $< -o $@
endif
##############################################################################

# EARTH SIMULATOR ############################################################
ifeq ($(strip $(ARCH)), ES)
mo_spitfire.o: mo_spitfire.f90
	$(F90) $(F90ALL) -pi auto exp=minmod exp=medan exp=putyslice exp=cfdot1dp2 exp=cfint1x2 line=2000 -c $< -o $@
mo_tpcore.o: mo_tpcore.f90
#warning: Don't change line to 2000. This will give wrong code on the SX!!!
	$(F90) $(F90ALL) -pi auto exp=xmist,fxppm,kmppm,lmppm,xtp noexp=map1_ppm_gp,ppm2m,steepz nest=3 line=1000 -c $< -o $@
mo_transpose.o: mo_transpose.f90
	$(F90) $(F90ALL) -Npi -c $< -o $@
lti.o: lti.f90
	$(F90) $(F90ALL) -Npi -c $< -o $@
endif
##############################################################################

# CRAY #######################################################################
ifeq ($(strip $(ARCH)), CRAY_PVP)
mo_buffer_fft.o: mo_buffer_fft.f90
	$(F90) $(F90ALL) -Ovector1 -c $< -o $@ 
mo_grib.o: mo_grib.f90
	$(F90) $(F90ALL) -Ovector1 -c $< -o $@ 
endif
##############################################################################

# Compaq Alpha / OSF1 ########################################################
ifeq ($(strip $(ARCH)), alpha)
mo_memory_base.o: mo_memory_base.f90
	$(F90) $(F90ALL) -O0 -c $< -o $@

mo_tracer.o: mo_tracer.f90
	$(F90) $(F90ALL) -O0 -c $< -o $@
endif
##############################################################################

# IBM-p690 / AIX #############################################################
ifeq ($(strip $(ARCH)), rs6000)
mo_tracer.o: mo_tracer.f90
	$(F90) $(F90ALL) -qnoopt -c $< -o $@
# NOTE: just exclude -qcheck
#       For some reason the compiler (11.1.0.4) reports 
#       an error in line 267 of the code. This is probably due to
#       optimisation, since with introducing "write(*,*) jn" the error
#       disappears.
mo_legendre.o: mo_legendre.f90
	$(F90) $(F90ALL) -qnocheck -c $< -o $@
## NOTE: Version 13.01.0000.0008 obviously has problems with run-time checks ...
##       when nudging is switched on ...
# mo_nmi.o: mo_nmi.f90
# 	$(F90) $(F90ALL) -qnocheck -qnoflttrap -c $< -o $@
# mo_exception.o: mo_exception.f90
# 	$(F90) $(F90ALL) -qnocheck -qnoflttrap -c $< -o $@
# mo_namelist.o: mo_namelist.f90
# 	$(F90) $(F90ALL) -qnocheck -qnoflttrap -c $< -o $@
# control.o: control.f90
# 	$(F90) $(F90ALL) -qnocheck -qnoflttrap -c $< -o $@
#mo_util_string.o: mo_util_string.f90
#	$(F90) $(F90ALL) -qnocheck -qnoflttrap -c $< -o $@
#inipost.o: inipost.f90
#	$(F90) $(F90ALL) -qnocheck -qnoflttrap -c $< -o $@
#mo_mpi.o: mo_mpi.f90
#	$(F90) $(F90ALL) -qnocheck -qnoflttrap -c $< -o $@
endif
##############################################################################

# IBM-Darwin / MAC ###########################################################
ifeq ($(strip $(ARCH)), darwin)
mo_real_timer.o: mo_real_timer.f90
	$(F90) $(F90ALL) -O0 -c $< -o $@ 

mo_couple.o: mo_couple.f90
	$(F90) $(F90ALL) -qnoextname -c $< -o $@

mo_column.o: mo_column.f90
	$(F90) $(F90ALL) -qnoextname -c $< -o $@

physc.o: physc.f90
	$(F90) $(F90ALL) -O0 -c $< -o $@ 

vdiff.o: vdiff.f90
	$(F90) $(F90ALL) -O0 -c $< -o $@ 

endif
############################################################################

# NAG compiler ###############################################################
ifeq ($(COMPILER), NAG)
# refers to compile error:
# 1) Error: mo_mpi.f90, line 1230: Inconsistent data type INTEGER
#      (previously REAL(real32)) for argument 1 in reference to MPI_SEND
# and downgrades it to a warning;
# and to runtime error:
# 2) Runtime Error: ../../messy/bmil/messy_main_control_echam5.inc, line 10:
#      Invalid reference to procedure MESSY_SETUP - 1 arguments expected but
#      it was called with 2 arguments
#    called in mo_mpi.f90, line 942: CALL messy_setup(0, p_all_comm)
mo_mpi.o: mo_mpi.f90
	$(F90) $(F90ALL) -mismatch -c $<
# refers to runtime error:
# 3) Runtime Error: ../src/pres.f90, line 1: Invalid reference to procedure
#      PRES - Actual argument for dummy argument PS (number 3) is not an array
physc.o: physc.f90
	$(F90) $(F90ALL) -mismatch -c $<
endif
############################################################################
