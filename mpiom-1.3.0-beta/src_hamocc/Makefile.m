# -*- Makefile -*-
# ----------------------------------------------

PREFIX       = ../..

LIB          = $(PREFIX)/lib/libhamocc.a
MAKEFILE_INC = depend.mk

### INCLUDES = $(addprefix $(MODOPT)../../, $(SRCDIRS)) \
### 	   $(addprefix $(MODOPT)../../, $(LIBSRCS)) \
###            $(ECHAM5_INCLUDE) -I../../include

INCLUDES     = $(NETCDF_INCLUDE) $(MPI2_INCLUDE) $(PNETCDF_INCLUDE) $(MPIOM_INCLUDE) 
F90NOR8      = $(F90FLAGS) $(F90DEFS) $(INCLUDES) -I../../messy/smcl  -I../src
# ----------------------------------------------------------------------
SRCS0 := $(wildcard *.f90)
SRCS  := $(filter-out F%.f90, $(SRCS0))
OBJS  := $(SRCS:.f90=.o)
MODS  := $(SRCS:.f90=.mod)
MODS_INST := $(addprefix $(PREFIX)/include/, $(SRCS))

F_makedepend =  ../../messy/util/sfmakedepend.pl --file=$(MAKEFILE_INC)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90R8) $(F90DEFSMPIOM) $(F90NOR8) $(INCLUDES) -c $< -o $@

# ----------------------------------------------------------------------

all: $(LIB)

$(LIB): depend $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)

# check files
list:
	@echo "------------------------------------------------"
	@echo "SRCS = $(SRCS)"
	@echo "------------------------------------------------"
	@echo
	@echo "------------------------------------------------"
	@echo "OBJS = $(OBJS)"
	@echo "------------------------------------------------"

depend $(MAKEFILE_INC): $(SRCS)
	$(F_makedepend)

install: all

clean:
	-rm -f $(OBJS)
	-rm -f *.mod
	-rm -f depend.mk
	-rm -f depend.mk.old

distclean: clean
	-rm -f $(LIB)
	-rm -f Makefile

# ----------------------------------------------------------------------
include $(MAKEFILE_INC)
#include specific.mk
# ----------------------------------------------------------------------
