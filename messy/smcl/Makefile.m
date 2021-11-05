# -*- Makefile -*-
# ----------------------------------------------

PREFIX       = ../..

LIB          = $(PREFIX)/lib/libmessy.a
MAKEFILE_INC = depend.mk

INCLUDES     = $(NETCDF_INCLUDE) $(MPI2_INCLUDE) $(PNETCDF_INCLUDE) $(CDI_INCLUDE) $(POSIX90_INCLUDE) $(T8CODE_INCLUDE) $(FORPY_INCLUDE) $(ASYNCF_INCLUDE)
F90NOR8      = $(F90FLAGS) $(F90DEFS) $(INCLUDES) $(SPEC_LIB)
# ----------------------------------------------------------------------
SRCS0 := $(wildcard *.f90)
SRCS  := $(filter-out F%.f90, $(SRCS0))
OBJS  := $(SRCS:.f90=.o)
MODS  := $(SRCS:.f90=.mod)
MODS_INST := $(addprefix $(PREFIX)/include/, $(SRCS))

F_makedepend =  ../util/sfmakedepend.pl --file=$(MAKEFILE_INC)

.SUFFIXES: $(SUFFIXES) .f90 .md5

# optional MECCA/KPP-CUDA integration
-include messy_mecca_kpp-cuda.mk

%.o: %.f90
	$(F90) $(F90NOR8) -c $< -o $@
# ----------------------------------------------------------------------

all: messy_main_compilerinfo_mem.f90 messy_main_tracer_chemprop.inc $(LIB)

$(LIB): depend $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	$(AR) -dv $(LIB) messy_ncregrid_interface.o

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
	-rm -f messy_main_compilerinfo_mem.f90
	-rm -f messy_main_compilerinfo_mem.f90.md5

# ----------------------------------------------------------------------
# THIS IS ALL TO GET THE COMPILER INFORMATION INTO THE EXECUTABLE(S)
# ----------------------------------------------------------------------
ifeq ($(strip $(MD5SUM_PRESENT)), yes)

%.md5: FORCE
	$(if $(filter-out $(shell cat $@ 2>/dev/null),$(shell md5sum $*)),md5sum $* > $@)

else

%.md5: FORCE
	$(if $(filter-out $(shell cat $@ 2>/dev/null),$(shell sum $*)),sum $* > $@)

endif

messy_main_compilerinfo_mem.f90.md5: messy_main_compilerinfo_mem.f90
messy_main_constants_mem.o: messy_main_compilerinfo_mem.f90.md5 messy_main_compilerinfo_mem.o
messy_main_compilerinfo_mem.o: messy_main_compilerinfo_mem.f90.md5 messy_main_compilerinfo_mem.f90

messy_main_compilerinfo_mem.f90: messy_main_compilerinfo_mem.f90.in FORCE
	@echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++...'
	@echo "COMPILER: $(F90VERS)"
	@echo "F90     : $(F90)"
	@echo "F90FLAGS: $(F90FLAGS)"
	@echo "F90DEFS : $(F90DEFS)"
	@echo "INCLUDES: $(INCLUDES)"
	@echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++...'
	-@rm -f messy_main_compilerinfo_mem.f90
	@sed -e "s|@F90VERS@|$(F90VERS)|g" \
	     -e "s|@F90@|$(F90)|g" \
	     -e "s|@F90FLAGS@|$(F90FLAGS)|g" \
	     -e "s|@F90DEFS@|$(F90DEFS)|g" \
	     -e "s|@INCLUDES@|$(INCLUDES)|g" \
	     -e '/.\{133\}/{' -e ': toolong' \
	     -e 'h' \
	     -e 's/\(.\{131\}\)\(.*\)/\1\&/' -e p -e g \
	     -e 's/\(.\{131\}\)\(.*\)/\&\2/' \
	     -e 't toolong' -e d \
	     -e '}' \
	     messy_main_compilerinfo_mem.f90.in \
	     > messy_main_compilerinfo_mem.f90
	@touch -r messy_main_compilerinfo_mem.f90.in messy_main_compilerinfo_mem.f90

messy_main_tracer_chemprop.inc: ../mbm/tracer/chemprop/messy_main_tracer_chemprop.tbl
	(cd ../mbm/tracer/chemprop; ./xchemprop)


FORCE:
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
include $(MAKEFILE_INC)
include specific.mk
# ----------------------------------------------------------------------
