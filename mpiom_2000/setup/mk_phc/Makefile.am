noinst_PROGRAMS=forcing_phc

forcing_phc_SOURCES=forcing_phc.f90 \
	../../support/iso_varying_string.f90 \
	../../src/mo_kind.f90 ../../src/mo_constants.f90

./$(DEPDIR)/FC.deps: $(forcing_phc_SOURCES) Makefile
	@unique=`for i in $(forcing_phc_SOURCES); do \
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
