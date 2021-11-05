AM_FCFLAGS = $(MPI) $(PRISM) \
	$(FPP_DEFOPT)SOR \
	$(FPP_DEFOPT)CLOCK \
	$(FPP_DEFOPT)QLOBERL $(FPP_DEFOPT)BULK_KARA $(FORCING_FC_DEFS) \
	$(FPP_DEFOPT)EISREST $(FPP_DEFOPT)REDWMICE \
	$(FPP_DEFOPT)DIAG $(FPP_DEFOPT)zzKONVDIAG \
	$(FPP_DEFOPT)zzAMLDDIAG \
	$(FPP_DEFOPT)zzONESIDED \
	$(HAMOCC) \
	$(OLD_IO_FCDEFS) \
	$(NEW_IO_FCDEFS) \
	$(CDI) $(CDI_INCLUDE) \
	$(FC_MOD_FLAG). $(FPP_INCOPT).\
	$(FC_MOD_FLAG)@top_srcdir@/src \
	$(FPP_INCOPT)@top_srcdir@/src \
	$(FC_MOD_FLAG)@top_srcdir@/support \
	$(FPP_INCOPT)@top_srcdir@/support \
	$(FC_MOD_FLAG)@top_builddir@/support \
	$(PRISM_INCLUDE) $(NETCDF_INCLUDE) $(MPI_INCLUDE) $(INCLUDES)

bin_PROGRAMS=mpiom.x

mpiom_x_SOURCES= ../src/absturz.f90					\
	../src/adisit1.f90						\
	../src/adisitj.f90						\
	../src/bartim.f90						\
	../src/beleg.f90						\
	../src/beleg_zero.f90						\
	../src/calcgmvis.f90						\
	../src/mo_contro.f90						\
	../src/cppoutput.f90						\
	../src/print_defines.f90					\
	../src/findbot.f90						\
	../src/growth.f90						\
	../src/mo_adpo.f90						\
	../src/mo_average_list.f90					\
	../src/mo_basin_masks.f90					\
	../src/mo_boundsexch.f90					\
	../src/mo_varlist.f90						\
	../src/mo_iolist.f90						\
	../src/mo_commo1.f90						\
	../src/mo_commoau1.f90						\
	../src/mo_commoau2.f90						\
	../src/mo_commoau3.f90						\
	../src/mo_commobbl.f90						\
	../src/mo_constants.f90						\
	../src/mo_couple.f90						\
	../src/mo_diagnosis.f90						\
	../src/mo_diffusion.f90						\
	../src/mo_eddydiag.f90						\
	../src/mo_elicom.f90						\
	../src/mo_file_list.f90						\
	../src/mo_fluxes1.f90						\
	../src/mo_forcing.f90						\
	../src/mo_grid.f90						\
	../src/mo_grid_elementals.f90					\
	../src/mo_kind.f90						\
	../src/mo_levitus.f90						\
	../src/mo_mean.f90						\
	../src/mo_model_time.f90                                        \
	../src/mo_mpi.f90						\
	../src/mo_ncar_ocean_fluxes.f90					\
	../src/mo_nudge_ts.f90						\
	../src/mo_number_map.f90					\
	../src/mo_ocice.f90						\
	../src/mo_octher.f90						\
	../src/mo_ocean_vertical_mixing.f90 ../src/mo_omip.f90		\
	../src/mo_para2.f90						\
	../src/mo_parallel.f90						\
	../src/mo_parallel_diags.f90					\
	../src/mo_param1.f90						\
	../src/mo_param3.f90						\
	../src/mo_planetary_constants.f90				\
	../src/mo_postprocess.f90 \
	../src/mo_profile.f90						\
	../src/mo_range_map.f90 \
	../src/mo_restart.f90						\
	../src/mo_rotation.f90					        \
	../src/mo_runoff.f90						\
	../src/mo_swr_absorption.f90					\
	../src/mo_tidal.f90						\
	../src/mo_tro.f90 ../src/mo_units.f90 ../src/mo_io_config.f90	\
	../src/mo_util_string.f90                      			\
	../src/mpiom.f90 ../src/mo_convection.f90 ../src/ocadfs.f90	\
	../src/ocbarp.f90 ../src/occlit.f90 ../src/occlit2.f90		\
	../src/ocjitr.f90 ../src/ocmodmom.f90 ../src/ocschep.f90	\
	../src/octimf.f90 ../src/ocuad.f90 ../src/ocvad.f90		\
	../src/mo_ocvisc.f90 ../src/ocvtot.f90 ../src/ocwind.f90		\
	../src/open_stdout.f90						\
	../src/rho.f90							\
	../src/rho1j.f90						\
	../src/rho1.f90						\
	../src/rho2.f90							\
	../src/rotate_2_north_east.f90 ../src/rotate_u.f90		\
	../src/rotate_v.f90						\
	../src/trian.f90						\
	../src/vcheck.f90						\
	../src/wrte_debug.f90						\
	../src/mo_wncdf.f90						\
	../src/wrte_konvdiag.f90					\
	../src/wrte_mean.f90						\
	../src/wrte_mfl.f90						\
									\
	../src_hamocc/mo_avflux.f90

nodist_mpiom_x_SOURCES=mo_iso_c_kinds.f90 mo_getenv.f90 version.c

EXTRA_mpiom_x_SOURCES=../src/mo_cpl_oasis4.f90				\
	../src/mo_io_backend_cdi.f90					\
	../src/mo_linked_list.inc					\
	../src/portability/mo_iso_c_kinds.f90.c				\
	../src/portability/mo_getenv.f90.c				\
									\
	../src_hamocc/aufr_bgc.f90					\
	../src_hamocc/mo_bgc_varlist.f90				\
	../src_hamocc/mo_bgc_iolist.f90 				\
	../src_hamocc/atmotr.f90					\
	../src_hamocc/aufw_bgc.f90					\
	../src_hamocc/avrg_bgcmean_2d.f90				\
	../src_hamocc/avrg_bgcmean_3d.f90				\
	../src_hamocc/avrg_dynamic.f90					\
	../src_hamocc/avrg_timeser_bgc.f90				\
	../src_hamocc/beleg_bgc.f90					\
	../src_hamocc/bgc.f90						\
	../src_hamocc/bodensed.f90					\
	../src_hamocc/c14_correction.f90				\
	../src_hamocc/carchm.f90					\
	../src_hamocc/carreset.f90					\
	../src_hamocc/chck_bgc.f90					\
	../src_hamocc/chemcon.f90					\
	../src_hamocc/close_bgcmean_2d.f90				\
	../src_hamocc/close_bgcmean_3d.f90				\
	../src_hamocc/close_bgcmean_bioz.f90				\
	../src_hamocc/close_bgcmean_sed.f90				\
	../src_hamocc/close_dynamic.f90					\
	../src_hamocc/compute_dyn_diff.f90				\
	../src_hamocc/compute_mld.f90					\
	../src_hamocc/compute_pco2.f90					\
	../src_hamocc/cyano.f90						\
	../src_hamocc/dilute_bgc.f90					\
	../src_hamocc/dipowa.f90					\
	../src_hamocc/end_bgc.f90					\
	../src_hamocc/extr.f90						\
	../src_hamocc/get_dust.f90					\
	../src_hamocc/ini_bgc.f90					\
	../src_hamocc/ini_timeser_bgc.f90				\
	../src_hamocc/inventory_bgc.f90					\
	../src_hamocc/inventory_bgc2.f90				\
	../src_hamocc/maschk.f90					\
	../src_hamocc/mo_bgc_diagnostic.f90				\
	../src_hamocc/mo_bgcmean.f90					\
	../src_hamocc/mo_biomod.f90					\
	../src_hamocc/mo_carbch.f90					\
	../src_hamocc/mo_control_bgc.f90				\
	../src_hamocc/mo_dynamic.f90					\
	../src_hamocc/mo_param1_bgc.f90					\
	../src_hamocc/mo_sedmnt.f90					\
	../src_hamocc/mo_timeser_bgc.f90				\
	../src_hamocc/netcdf_def_vardb.f90				\
	../src_hamocc/netcdf_def_varsg.f90				\
	../src_hamocc/ocprod.f90					\
	../src_hamocc/open_bgcmean_2d.f90				\
	../src_hamocc/open_bgcmean_3d.f90				\
	../src_hamocc/open_bgcmean_bioz.f90				\
	../src_hamocc/open_bgcmean_sed.f90				\
	../src_hamocc/open_dynamic.f90					\
	../src_hamocc/powach.f90					\
	../src_hamocc/powadi.f90					\
	../src_hamocc/read_namelist.f90					\
	../src_hamocc/read_netcdf_var.f90				\
	../src_hamocc/save_timeser_bgc.f90				\
	../src_hamocc/sedshi.f90					\
	../src_hamocc/spinup_bgc.f90					\
	../src_hamocc/write_bgcmean_2d.f90				\
	../src_hamocc/write_bgcmean_3d.f90				\
	../src_hamocc/write_bgcmean_bioz.f90				\
	../src_hamocc/write_bgcmean_sed.f90				\
	../src_hamocc/write_dynamic.f90					\
	../src_hamocc/write_netcdf_var.f90				\
	../src_hamocc/wrte_meanbgc.f90

BUILT_SOURCES=mo_iso_c_kinds.f90 mo_getenv.f90 defines.list version.c
CLEANFILES=mo_iso_c_kinds.f90 mo_getenv.f90

mo_iso_c_kinds.f90: ../src/portability/mo_iso_c_kinds.f90.c
	$(CPP) $(DEFS) $(CPPFLAGS) $(DEFAULT_INCLUDES) `test -f ../src/portability/mo_iso_c_kinds.f90.c || echo $(srcdir)/`../src/portability/mo_iso_c_kinds.f90.c | @GREP@ -v '^#' >$@

mo_getenv.f90: ../src/portability/mo_getenv.f90.c
	$(CPP) $(DEFS) $(CPPFLAGS) $(DEFAULT_INCLUDES) `test -f ../src/portability/mo_getenv.f90.c || echo $(srcdir)/`../src/portability/mo_getenv.f90.c | @GREP@ -v '^#' >$@

mpiom_x_LDADD=../support/libsupport.a \
	@PRISM_LIB@ @CDI_LIB@ @NETCDF_LIB@ @MPI_LIB@ @LIBS@

if COUPLE_OASIS4
mpiom_x_SOURCES += ../src/mo_cpl_oasis4.f90
endif

if USEPBGC
HAMOCC_SOURCES = ../src_hamocc/aufr_bgc.f90 ../src_hamocc/atmotr.f90	\
	../src_hamocc/aufw_bgc.f90 ../src_hamocc/avrg_bgcmean_2d.f90	\
	../src_hamocc/mo_bgc_varlist.f90				\
	../src_hamocc/mo_bgc_iolist.f90 				\
	../src_hamocc/avrg_bgcmean_3d.f90				\
	../src_hamocc/avrg_dynamic.f90					\
	../src_hamocc/avrg_timeser_bgc.f90 ../src_hamocc/beleg_bgc.f90	\
	../src_hamocc/bgc.f90 ../src_hamocc/bodensed.f90		\
	../src_hamocc/c14_correction.f90 ../src_hamocc/carchm.f90	\
	../src_hamocc/carreset.f90		\
	../src_hamocc/chck_bgc.f90 ../src_hamocc/chemcon.f90		\
	../src_hamocc/close_bgcmean_2d.f90				\
	../src_hamocc/close_bgcmean_3d.f90				\
	../src_hamocc/close_bgcmean_bioz.f90				\
	../src_hamocc/close_bgcmean_sed.f90				\
	../src_hamocc/close_dynamic.f90					\
	../src_hamocc/compute_dyn_diff.f90				\
	../src_hamocc/compute_mld.f90 ../src_hamocc/compute_pco2.f90	\
	../src_hamocc/cyano.f90	\
	../src_hamocc/dilute_bgc.f90 ../src_hamocc/dipowa.f90		\
	../src_hamocc/end_bgc.f90 ../src_hamocc/extr.f90		\
	../src_hamocc/get_dust.f90 ../src_hamocc/ini_bgc.f90		\
	../src_hamocc/ini_timeser_bgc.f90				\
	../src_hamocc/inventory_bgc.f90					\
	../src_hamocc/inventory_bgc2.f90 ../src_hamocc/maschk.f90	\
	../src_hamocc/mo_bgc_diagnostic.f90				\
	../src_hamocc/mo_bgcmean.f90 ../src_hamocc/mo_biomod.f90	\
	../src_hamocc/mo_carbch.f90 ../src_hamocc/mo_control_bgc.f90	\
	../src_hamocc/mo_dynamic.f90 ../src_hamocc/mo_param1_bgc.f90	\
	../src_hamocc/mo_sedmnt.f90 ../src_hamocc/mo_timeser_bgc.f90	\
	../src_hamocc/netcdf_def_vardb.f90				\
	../src_hamocc/netcdf_def_varsg.f90 ../src_hamocc/ocprod.f90	\
	../src_hamocc/open_bgcmean_2d.f90				\
	../src_hamocc/open_bgcmean_3d.f90				\
	../src_hamocc/open_bgcmean_bioz.f90				\
	../src_hamocc/open_bgcmean_sed.f90				\
	../src_hamocc/open_dynamic.f90 ../src_hamocc/powach.f90		\
	../src_hamocc/powadi.f90 ../src_hamocc/read_namelist.f90	\
	../src_hamocc/read_netcdf_var.f90				\
	../src_hamocc/save_timeser_bgc.f90 ../src_hamocc/sedshi.f90	\
	../src_hamocc/spinup_bgc.f90					\
	../src_hamocc/write_bgcmean_2d.f90				\
	../src_hamocc/write_bgcmean_3d.f90				\
	../src_hamocc/write_bgcmean_bioz.f90				\
	../src_hamocc/write_bgcmean_sed.f90				\
	../src_hamocc/write_dynamic.f90					\
	../src_hamocc/write_netcdf_var.f90				\
	../src_hamocc/wrte_meanbgc.f90
mpiom_x_SOURCES += $(HAMOCC_SOURCES)
endif
if USE_CDI_IO
mpiom_x_SOURCES += \
	../src/mo_io_backend_cdi.f90
endif

clean-hamocc:
	@hamocc_obj=`for i in $(HAMOCC_SOURCES) \
	  ../src_hamocc/mo_avflux.f90 ; do \
	  echo "$$i" | sed -e 's:.*/::;s/.f90$$/.$(OBJEXT)/' ; \
	  echo "$$i" | sed -e 's:.*/::;s/.f90$$/.$(FCMODEXT)/' ; done`; \
	  echo rm -f $$hamocc_obj ; \
	  rm -f $$hamocc_obj

#
# Retrieve version control information for inclusion into executable.
#
version.c: $(mpiom_x_SOURCES) Makefile
	@echo retrieving version control information
	@echo "char const version_control_information[] = \"\$$Revision: `$(top_srcdir)/util/vc_info $(top_srcdir)` \$$\";" > $@

all-local: ./$(DEPDIR)/FC.deps

# This ultra-evil hack allows us to insert a file depending on all
# other sources into the list of sources without causing circular
# dependencies.
defines.list: cpp_keys_list
	@echo post-processing pre-processor information
	@sort -f  cpp_keys_list | \
	  $(top_srcdir)/util/expand_template >defines.list
cpp_keys_list: $(mpiom_x_SOURCES) $(nodist_mpiom_x_SOURCES)
	@echo retrieving pre-processor information
	@unique=`for i in $(mpiom_x_SOURCES) $(nodist_mpiom_x_SOURCES); do \
	  if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | grep '\.f90$$' | grep -v '^print_defines.f90$$'` ; \
	  FC="$(FC)" FPP="$(FPP)" \
	  @top_builddir@/config/makef90depends --ignore-module iso_c_binding \
	  --noprint-obj-deps \
	  --dump-cpp-keys=cpp_keys_list -- \
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) $(AM_FPPFLAGS) -- $$unique >/dev/null


./$(DEPDIR)/FC.deps: $(mpiom_x_SOURCES) $(nodist_mpiom_x_SOURCES) Makefile
	@echo retrieving dependency and pre-processor information
	@unique=`for i in $(mpiom_x_SOURCES) $(nodist_mpiom_x_SOURCES); do \
	  if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | grep '\.f90$$' ` ; \
	  FC="$(FC)" FPP="$(FPP)" \
	  @top_builddir@/config/makef90depends --ignore-module iso_c_binding \
	  --ignore-module mo_fpe \
	  --strip-obj-dirname \
	  --dump-cpp-keys=cpp_keys_list -- \
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) $(AM_FPPFLAGS) -- $$unique >./$(DEPDIR)/FC.deps

@AMDEP_TRUE@@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@

clean-local:
	-rm -f *.$(FCMODEXT) i.*.L
