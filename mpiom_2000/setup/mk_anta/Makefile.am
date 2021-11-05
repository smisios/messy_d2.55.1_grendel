noinst_PROGRAMS=gridgen zeko

gridgen_SOURCES=gridgen.f90 ../../src/mo_kind.f90 \
	../../src/mo_constants.f90 ../../support/iso_varying_string.f90

zeko_SOURCES=zeko.f90 ../../src/mo_kind.f90 \
	../../src/mo_constants.f90 ../../support/iso_varying_string.f90

EXTRA_DIST = readme.mk_anta

./$(DEPDIR)/FC.deps: $(gridgen_SOURCES) $(zeko_SOURCES) Makefile
	@unique=`for i in $(gridgen_SOURCES) $(zeko_SOURCES); do \
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


