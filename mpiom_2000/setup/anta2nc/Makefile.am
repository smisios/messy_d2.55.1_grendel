noinst_PROGRAMS=anta2nc

AM_FPPFLAGS=$(NETCDF_INCLUDE)

anta2nc_SOURCES=anta2nc.f90 ../../support/iso_varying_string.f90 \
	../../src/mo_kind.f90

anta2nc_LDADD=$(NETCDF_LIB)

./$(DEPDIR)/FC.deps: $(anta2nc_SOURCES) Makefile
	@unique=`for i in $(anta2nc_SOURCES); do \
	  if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | grep '\.f90$$' ` ; \
	  FC="$(FC)" FCFLAGS="$(FCFLAGS)" FPPFLAGS="$(FPPFLAGS)" \
	  @top_builddir@/config/makef90depends --ignore-module iso_c_binding \
	  --strip-obj-dirname \
	  --dump-cpp-keys=cpp_keys_list -- \
	  $(AM_CPPFLAGS) $(CPPFLAGS) $(FCFLAGS) -- $$unique >./$(DEPDIR)/FC.deps

@AMDEP_TRUE@@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@

clean-local:
	-rm -f *.mod i.*.L
