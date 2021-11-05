### THIS FILE CONTAINS SPECIFIC RULES AND IS INCLUDED IN THE MESSy Makefile

## DO COMPILE WITH THE -r8 (or equivalent) OPTION
messy_convect_donner_additions.o: messy_convect_donner_additions.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_deep_k.o: messy_convect_donner_deep_k.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_types.o: messy_convect_donner_types.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_cape_k.o: messy_convect_donner_cape_k.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_lscloud_k.o: messy_convect_donner_lscloud_k.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_util.o: messy_convect_donner_util.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_cloudmodel.o: messy_convect_donner_cloudmodel.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_meso_k.o: messy_convect_donner_meso_k.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_closure.o: messy_convect_donner_closure.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect_donner_rad_k.o: messy_convect_donner_rad_k.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<
messy_convect.o: messy_convect.f90
	$(F90) $(F90R8) $(F90NOR8) -c $<

## ARCHITECTURE SPECIFIC

# Compaq Alpha / OSF1 ########################################################
ifeq ($(ARCH), alpha)
messy_main_tools.o: messy_main_tools.f90
	$(F90) $(F90NOR8) -O0 -c $<
endif
##############################################################################

# Linux64 / Intel Compiler ###################################################
ifeq ($(ARCH), LINUX64)
ifeq ($(COMPILER), INTEL)

# NOTE: just exclude -check
messy_mecca_aero.o: messy_mecca_aero.f90
	$(F90) $(F90NOR8) -WB -c $<
messy_mecca.o: messy_mecca.f90
	$(F90) $(F90NOR8) -WB -c $<
messy_mecca_kpp.o: messy_mecca_kpp.f90
	$(F90) $(F90NOR8) -WB -c $<
messy_mtchem_kpp.o: messy_mtchem_kpp.f90
	$(F90) $(F90NOR8) -WB -c $<
## mz_ab_20150420
#messy_main_import_grid.o: messy_main_import_grid.f90
#	$(F90) $(F90NOR8) -O0 -debug all -check all -traceback -c $<

endif
endif
##############################################################################

# IBM-p690 / AIX #############################################################
ifeq ($(ARCH), rs6000)
### NOTE: PHOTO, JVAL with full optimisation produce wrong results ...
messy_jval.o: messy_jval.f90
	$(F90) $(F90NOR8) -qnoopt -c $<
messy_jvst.o: messy_jvst.f90
	$(F90) $(F90NOR8) -qnoopt -c $<
messy_mecca_kpp.o: messy_mecca_kpp.f90
	$(F90) $(F90NOR8) -Q-kppdecomp -c $<
#messy_mecca_kpp.o: messy_mecca_kpp.f90
#	$(F90) $(F90NOR8) -qnoopt -c $<

# NOTE: just exclude -qcheck
messy_mecca.o: messy_mecca.f90
	$(F90) $(F90NOR8) -qnocheck -c $<
messy_mecca_kpp.o: messy_mecca_kpp.f90
	$(F90) $(F90NOR8) -qnocheck -c $<
messy_mecca%_kpp.o: messy_mecca%_kpp.f90
	$(F90) $(F90NOR8) -qnocheck -c $<
messy_mtchem_kpp.o: messy_mtchem_kpp.f90
	$(F90) $(F90NOR8) -qnocheck -c $<

# NOTE: exclude -qhot
messy_convect_donner_deep_k.o: messy_convect_donner_deep_k.f90
	$(F90) $(F90R8) $(F90NOR8) -qnohot -c $<
messy_cloud_ori.o: messy_cloud_ori.f90
	$(F90) $(F90NOR8) -qnohot -c $<


### compile always optimized and without checks ...
LF90=$(F90)
LFLAGS=-q64 -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -O3 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qspillsize=32648 -qarch=auto -qtune=auto -Q -qhot -qxlf90=nosignedzero -bdatapsize:64k -bstackpsize:64k -btextpsize:64k
#
messy_main_grid_trafo_nrgd.o: messy_main_grid_trafo_nrgd.f90
	$(LF90) $(LFLAGS) $(F90DEFS) $(INCLUDES) -c $<
messy_main_grid_trafo_ncrg_base.o: messy_main_grid_trafo_ncrg_base.f90
	$(LF90) $(LFLAGS) $(F90DEFS) $(INCLUDES) -c $<
messy_main_grid_trafo_scrp.o: messy_main_grid_trafo_scrp.f90
	$(LF90) $(LFLAGS) $(F90DEFS) $(INCLUDES) -c $<
messy_main_grid_trafo_scrp_base.o: messy_main_grid_trafo_scrp_base.f90
	$(LF90) $(LFLAGS) $(F90DEFS) $(INCLUDES) -c $<

messy_main_grid.o: messy_main_grid.f90
	$(LF90) $(LFLAGS) $(F90DEFS) $(INCLUDES) -c $<
messy_main_grid_netcdf.o: messy_main_grid_netcdf.f90
	$(LF90) $(LFLAGS) $(F90DEFS) $(INCLUDES) -c $<
messy_main_grid_trafo.o: messy_main_grid_trafo.f90
	$(LF90) $(LFLAGS) $(F90DEFS) $(INCLUDES) -c $<
messy_main_grid_tools.o: messy_main_grid_tools.f90
	$(LF90) $(LFLAGS) $(F90DEFS) $(INCLUDES) -c $<
endif
##############################################################################

# NEC-SX6 / SUPER-UX #########################################################
ifeq ($(ARCH), SX)

LIST=-Wf,-L noinclist fmtlist objlist source mrgmsg summary transform map,-pvctl fullmsg,-msg o
INLINE=-pi auto nest=10 line=1000 fullmsg
OPT=-Chopt -Pstack -Wf,-pvctl vwork=stack
INIT=-Wf,-init stack=zero heap=zero -Wf,-K a
# -Ep
# -Wl,-Z 1000000

messy_main_compilerinfo_mem.o: messy_main_compilerinfo_mem.f90
	$(F90) $(F90NOR8) -f5 -c $< -o $@
messy_main_tools.o: messy_main_tools.f90
	$(F90) $(F90NOR8) -f2003 nocbind -c $< -o $@

messy_scav_l_kpp.o: messy_scav_l_kpp.f90
	$(F90) $(INCLUDES) $(PERF) $(LIST) -pi noauto $(OPT) $(INIT) $(F90DEFS) -c $<

messy_scav_inter.o: messy_scav_inter.f90
	$(F90) $(INCLUDES) $(PERF) $(LIST) -pi noauto $(OPT) $(INIT) $(F90DEFS) -c $<

messy_scav_i_kpp.o: messy_scav_i_kpp.f90
	$(F90) $(INCLUDES) $(PERF) $(LIST) -pi noauto $(OPT) $(INIT) $(F90DEFS) -c $<

endif
##############################################################################

# Linux compiler #############################################################
ifneq ($(findstring $(strip $(ARCH)), LINUX LINUX64), )
#messy_mecca_kpp.o: messy_mecca_kpp.f90
#	$(F90) $(F90NOR8) -O0 -c $<
#messy_scav_l_kpp.o: messy_scav_l_kpp.f90
#	$(F90) $(F90NOR8) -O0 -c $<

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

# messy_main_grid.o: messy_main_grid.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_mpi.o: messy_main_grid_mpi.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_netcdf.o: messy_main_grid_netcdf.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_tools.o: messy_main_grid_tools.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_trafo.o: messy_main_grid_trafo.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_trafo_nrgd_base.o: messy_main_grid_trafo_nrgd_base.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_trafo_nrgd.o: messy_main_grid_trafo_nrgd.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_trafo_scrp_base.o: messy_main_grid_trafo_scrp_base.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_grid_trafo_scrp.o: messy_main_grid_trafo_scrp.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import.o: messy_main_import.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_grid.o: messy_main_import_grid.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_grid_par.o: messy_main_import_grid_par.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_rgt.o: messy_main_import_rgt.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<
# messy_main_import_ts.o: messy_main_import_ts.f90
#	$(F90) $(FOPT) $(F90DEFS) $(INCLUDES) -c $<

endif
##############################################################################
