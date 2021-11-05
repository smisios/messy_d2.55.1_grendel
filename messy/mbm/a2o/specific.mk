## NEC-SX6 / SUPER-UX ########################################################
ifeq ($(ARCH), SX)
messy_ncregrid_netcdf.o: messy_ncregrid_netcdf.f90
	$(F90) $(INCLUDES) -Wf,-L fmtlist transform map summary -Wl,-Z 1000000 -Cvsafe $(F90DEFS) -c $<

messy_ncregrid_geohyb.o: messy_ncregrid_geohyb.f90
	$(F90) $(INCLUDES) -Wf,-L fmtlist transform map summary -Wl,-Z 1000000 -Chopt $(F90DEFS) -c $<
endif
##############################################################################
