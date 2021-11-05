### THIS FILE CONTAINS SPECIFIC RULES AND IS INCLUDED IN THE ECHAM5 Makefile

### ECHAM5/MESSy - ARCHITECTURE SPECIFIC #####################################

# Compaq Alpha / OSF1 ########################################################
ifeq ($(strip $(ARCH)), alpha)
messy_ncregrid_tools_bi.o: messy_ncregrid_tools_bi.f90
	$(F90) $(F90ALL) -O0 -c $< -o $@

endif
##############################################################################

# Linux64 / Intel Compiler ###################################################
# ifeq ($(ARCH), LINUX64)
# ifeq ($(COMPILER), INTEL)

# endif
# endif
##############################################################################

# IBM-p690 / AIX #############################################################
ifeq ($(strip $(ARCH)), rs6000)
#messy_ncregrid_tools_bi.o: messy_ncregrid_tools_bi.f90
#	$(F90) $(F90ALL) -qnoopt -c $<  -o $@
### NOTE: PHOTO with full optimisation produces wrong results ...
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

messy_convect_e5.o: messy_convect_e5.f90
	$(F90)  -Wf,-A dbl4 $(INCLUDES) -Wf,-L fmtlist transform map summary,-pvctl fullmsg vwork=stack vworksz=4M -pi auto expin=../../messy/smcl/messy_convect.f90 nest=3 line=1000 -Chopt $(F90DEFS) -c $< -o $@

messy_psc_e5.o: messy_psc_e5.f90
	$(F90) -Wf,-A idbl4 $(INCLUDES) $(LIST) $(INIT) $(OPT) -pi auto expin=../../messy/smcl/messy_psc.f90 nest=3 line=1000 -Chopt $(F90DEFS) -c $< -o $@

endif
##############################################################################

# Linux compiler #############################################################
ifneq ($(findstring $(strip $(ARCH)), LINUX LINUX64), )
### default:
FOPT=$(F90FLAGS)
### compile always without checks ...
ifneq ($(findstring $(strip $(COMPILER)), LF95), )
   FOPT=--ap -O3 -Cpp
endif
ifneq ($(findstring $(strip $(COMPILER)), G95), )
   FOPT=-cpp -O3 -fno-second-underscore -ffree-line-length-huge -fno-backslash
endif
ifneq ($(findstring $(strip $(COMPILER)), INTEL), )
   FOPT=-fpp -O2 -fp-model strict -fno-alias -no-ansi-alias -lpthread -save-temps
endif
ifneq ($(findstring $(strip $(COMPILER)), GFORTRAN), )
   FOPT=-cpp -fno-second-underscore -ffree-line-length-none -fno-range-check -O3
endif
messy_ncregrid_interface.o: messy_ncregrid_interface.f90
	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
messy_ncregrid_tools_bi.o: messy_ncregrid_tools_bi.f90
	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
endif
##############################################################################

# IBM-Darwin / MAC ###########################################################
ifeq ($(strip $(ARCH)), darwin)
messy_ncregrid_tools_bi.o: messy_ncregrid_tools_bi.f90
	$(F90) $(F90ALL) -O0 -c $< -o $@ 

endif
############################################################################
