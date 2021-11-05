## NEC-SX6 / SUPER-UX ########################################################
ifeq ($(ARCH), SX)
# -Wl,-Z 1000000
messy_main_grid_netcdf.o: messy_main_grid_netcdf.f90
	$(F90) $(INCLUDES) -Wf,-L fmtlist transform map summary -Cvsafe $(F90DEFS) -c $<

messy_main_grid.o: messy_main_grid.f90
	$(F90) $(INCLUDES) -Wf,-L fmtlist transform map summary -Chopt $(F90DEFS) -c $<
endif
##############################################################################
