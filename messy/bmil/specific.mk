### THIS FILE CONTAINS SPECIFIC RULES AND IS INCLUDED IN THE ECHAM5 Makefile

# NEC-SX6 / SUPER-UX #########################################################
ifeq ($(strip $(ARCH)), SX)
# -Wl,-Z 1000000
messy_main_channel_bi.o: messy_main_channel_bi.f90
	$(F90) -Wf,-A idbl4 $(INCLUDES) -Wf,-L fmtlist transform map summary -Chopt $(F90DEFS) -c $< -o $@ 
endif
##############################################################################

# Linux compiler #############################################################
ifneq ($(findstring $(strip $(ARCH)), LINUX LINUX64), )

# default:
FOPT=$(F90FLAGS)
# compile always without checks ...
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

# messy_main_grid_bi.o: messy_main_grid_bi.f90
# 	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_netcdf_bi.o: messy_main_grid_netcdf_bi.f90
# 	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_bi.o: messy_main_import_bi.f90
# 	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_grid_bi.o: messy_main_import_grid_bi.f90
# 	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_grid_tools_bi.o: messy_main_import_grid_tools_bi.f90
# 	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_rgt_bi.o: messy_main_import_rgt_bi.f90
# 	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_ts_bi.o: messy_main_import_ts_bi.f90
# 	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<

endif
##############################################################################
