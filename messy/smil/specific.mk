### THIS FILE CONTAINS SPECIFIC RULES AND IS INCLUDED IN THE BASEMODEL Makefile

# ### COMPILE WITHOUT -r8, --dbl or similar ..., i.e. in standard real precision
# messy_jval_si.o: messy_jval_si.f90
# 	$(F90) $(F90NOR8) -c $<
# messy_jvst_si.o: messy_jvst_si.f90
# 	$(F90) $(F90NOR8) -c $<

# messy_cloud_si.o: messy_cloud_si.f90
# 	$(F90) $(F90NOR8) -c $<

# Linux64 / Intel Compiler ###################################################
ifeq ($(ARCH), LINUX64)
ifeq ($(COMPILER), INTEL)

# NOTE: compiler error workaround (v14.1, v15.0), reduced optimisation required
messy_msbm_si.o: messy_msbm_si.f90
	$(F90) $(F90ALL) -O0 -no-vec -c $<  -o $@

#BULL: internal compiler error with intel v15.0.2
#messy_made3_si.o: messy_made3_si.f90
#	$(F90) $(F90ALL) -O0 -no-vec -c $<  -o $@

endif
endif
##############################################################################

# IBM-p690 / AIX #############################################################
ifeq ($(strip $(ARCH)), rs6000)

messy_scav_si.o: messy_scav_si.f90
#	$(F90) $(F90ALL) -qnoopt -c $<  -o $@
	$(F90) -qrealsize=8 $(INCLUDES) -q64 -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -O2 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qfree=f90 -qspillsize=32648 -Q $(F90DEFS) -c $<  -o $@

#messy_onemis_si.o: messy_onemis_si.f90
#	$(F90) $(F90ALL) -qnoopt -c $<  -o $@

### NOTE: JVAL with full optimisation produces wrong results ...
messy_jval_si.o: messy_jval_si.f90
	$(F90) $(F90NOR8) -qnoopt -c $<  -o $@
messy_jvst_si.o: messy_jvst_si.f90
	$(F90) $(F90NOR8) -qnoopt -c $<  -o $@

endif
##############################################################################

# NEC-SX6 / SUPER-UX #########################################################
ifeq ($(strip $(ARCH)), SX)

LIST=-Wf,-L noinclist fmtlist objlist source mrgmsg summary transform map,-pvctl fullmsg,-msg o
INLINE=-pi auto nest=10 line=1000 fullmsg
OPT=-Chopt -Pstack -Wf,-pvctl vwork=stack
INIT=-Wf,-init stack=zero heap=zero -Wf,-K a
# -Wl,-Z 1000000
# -Ep

messy_s4d_si.o: messy_s4d_si.f90
	$(F90) -Wf,-A dbl4 $(INCLUDES) -Wf,-L fmtlist transform map summary,-pvctl fullmsg vwork=stack vworksz=4M $(F90DEFS) -Cvsafe -c $< -o $@ 

messy_mecca_aero_si.o: messy_mecca_aero_si.f90
	$(F90)  -Wf,-A dbl4 $(INCLUDES) $(F90DEFS) $(LIST) $(INIT) -Cvsafe -Pstack -Wf,-pvctl vwork=stack -c $< -o $@ 

messy_dradon_si.o: messy_dradon_si.f90
	$(F90) -Wf,-A idbl4 $(INCLUDES) -Wf,-L fmtlist transform map summary,-pvctl fullmsg vwork=stack vworksz=4M -Chopt $(F90DEFS) -c $< -o $@ 

messy_sedi_si.o: messy_sedi_si.f90
	$(F90) -Wf,-A idbl4 $(INCLUDES) -Wf,-L fmtlist transform map summary,-pvctl fullmsg vwork=stack vworksz=4M -pi auto expin=../../messy/smcl/messy_sedi.f90 nest=3 line=1000 -Chopt $(F90DEFS) -c $< -o $@

messy_scav_si.o: messy_scav_si.f90
	$(F90) -Wf,-A idbl4 $(INCLUDES) -Wf,-L fmtlist transform map summary,-pvctl fullmsg vwork=stack vworksz=4M -pi auto expin=../../messy/smcl/messy_scav.f90,../../messy/smcl/messy_scav_aer.f90 nest=3 line=10000 -Chopt $(F90DEFS) -c $< -o $@

endif
##############################################################################

# IBM-Darwin / MAC ###########################################################
ifeq ($(strip $(ARCH)), darwin)

messy_jval_si.o: messy_jval_si.f90
	$(F90) $(F90NOR8) -O0 -c $< -o $@ 

endif
############################################################################
