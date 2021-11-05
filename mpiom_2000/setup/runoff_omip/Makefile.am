noinst_PROGRAMS=runoff_omip

runoff_omip_SOURCES=runoff_omip.f90 ../../src/mo_constants.f90 \
	../../src/mo_kind.f90 ../../src/mo_planetary_constants.f90 \
	../../src/mo_grid_elementals.f90

./$(DEPDIR)/FC.deps: $(runoff_omip_SOURCES) Makefile
	@unique=`for i in $(runoff_omip_SOURCES); do \
	  if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | grep '\.f90$$' ` ; \
	  FC="$(FC)" FCFLAGS="$(FCFLAGS)" FPPFLAGS="$(FPPFLAGS)" \
	  @top_builddir@/config/makef90depends --ignore-module iso_c_binding \
	  --strip-obj-dirname \
	  --dump-cpp-keys=cpp_keys_list -- \
	  $(AM_CPPFLAGS) $(CPPFLAGS) $(FCFLAGS) -- $$unique >./$(DEPDIR)/FC.deps

@AMDEP_TRUE@@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@

clean-local:
	-rm -f *.@FCMODEXT@ i.*.L
