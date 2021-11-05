# -*- Makefile -*-

PREFIX          = ../..

LIB		= $(PREFIX)/lib/libnetcdf90.a
F90NOR8		= $(INCLUDES) $(F90FLAGS) $(F90DEFS)

OBJS	        = netcdf.o typeSizes.o
MODS      := $(OBJS:.o=.mod)
MODS_INST := $(addprefix $(PREFIX)/include/, $(MODS))

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)

typeSizes.o:	typeSizes.f90
	$(F90) $(F90NOR8) -c typeSizes.f90

netcdf.o:	\
	    netcdf.f90 typeSizes.o netcdf_constants.f90 netcdf_externals.f90 \
	    netcdf_overloads.f90 netcdf_visibility.f90 netcdf_file.f90 \
	    netcdf_dims.f90 netcdf_attributes.f90 netcdf_variables.f90 \
	    netcdf_text_variables.f90 netcdf_expanded.f90
	$(F90) $(F90NOR8) -c netcdf.f90

install: all
	cp -pf *.mod $(PREFIX)/include/.

clean:
	rm -f $(OBJS) *.mod

distclean: clean
	rm -f $(MODS_INST) $(PREFIX)/include/typesizes.mod
	rm -f $(LIB)
