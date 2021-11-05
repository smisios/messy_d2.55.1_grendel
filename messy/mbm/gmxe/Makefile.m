# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .
PROG       = gmxe.exe

# ----------------------------------------------
SRCS =  messy_main_constants_mem.f90 \
 messy_main_tools.f90 \
 messy_main_blather.f90 \
 messy_gmxe_mem.f90 \
 messy_main_tools_kinetics.f90 \
 messy_gmxe_aerchem_inp_kpp.f90 \
 messy_gmxe_aerchem_kpp.f90 \
 messy_gmxe_aerchem_liq.f90 \
 messy_gmxe_aerchem.f90 \
 messy_gmxe_kappa.f90 \
 messy_gmxe_eqsam4clim.f90 \
 messy_gmxe_isorropia2.f90 \
 messy_gmxe_oc_aging.f90 \
 messy_gmxe_soa.f90 \
 messy_gmxe.f90 \
 messy_gmxe_box.f90

OBJS := $(SRCS:.f90=.o)

GMXE_LIB     = $(NETCDF_LIB) -L../../../lib/ -lisorropia
GMXE_INCLUDE = $(NETCDF_INCLUDE)

ifeq ($(strip $(POSIX90_DEF)), HAVE_POSIX90)
GMXE_INCLUDE += $(MODOPT) ../../../libsrc/posix90/src
GMXE_LIB     += -L ../../../lib -lposix90
endif

all: $(PROG)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*~
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f $(srcdir)/fort.[0-9]*

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(F90DEFS) $(OBJS) $(GMXE_LIB) -o $@

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90DEFS) $(F90FLAGS) $(GMXE_INCLUDE) -c $<

# ----------------------------------------------
messy_gmxe_aerchem.o : messy_gmxe_aerchem.f90 messy_main_tools.o messy_gmxe_aerchem_inp_kpp.o messy_gmxe_aerchem_kpp.o messy_gmxe_aerchem_liq.o messy_main_constants_mem.o messy_gmxe_mem.o
messy_gmxe_aerchem_inp_kpp.o : messy_gmxe_aerchem_inp_kpp.f90 messy_main_tools_kinetics.o messy_main_constants_mem.o
messy_gmxe_aerchem_kpp.o : messy_gmxe_aerchem_kpp.f90 messy_main_tools.o messy_gmxe_aerchem_inp_kpp.o messy_main_constants_mem.o
messy_gmxe_aerchem_liq.o : messy_gmxe_aerchem_liq.f90 messy_main_constants_mem.o messy_main_tools.o messy_gmxe_aerchem_kpp.o
messy_gmxe_box.o : messy_gmxe_box.f90 messy_gmxe_kappa.o messy_gmxe_aerchem_kpp.o messy_gmxe_aerchem_liq.o messy_gmxe_oc_aging.o messy_gmxe_eqsam4clim.o messy_gmxe_isorropia2.o messy_gmxe_soa.o messy_gmxe_aerchem.o messy_main_tools.o messy_gmxe.o messy_gmxe_mem.o messy_main_constants_mem.o
messy_gmxe_eqsam4clim.o : messy_gmxe_eqsam4clim.f90 messy_gmxe_kappa.o messy_gmxe_mem.o messy_main_constants_mem.o
messy_gmxe.o : messy_gmxe.f90 messy_gmxe_eqsam4clim.o messy_gmxe_isorropia2.o messy_gmxe_kappa.o messy_gmxe_aerchem.o messy_gmxe_oc_aging.o messy_main_tools.o messy_main_constants_mem.o messy_gmxe_soa.o messy_gmxe_mem.o messy_main_blather.o
messy_gmxe_isorropia2.o : messy_gmxe_isorropia2.f90 messy_gmxe_kappa.o messy_gmxe_mem.o messy_main_constants_mem.o
messy_gmxe_kappa.o : messy_gmxe_kappa.f90 messy_gmxe_mem.o messy_main_constants_mem.o
messy_gmxe_mem.o : messy_gmxe_mem.f90 messy_main_constants_mem.o
messy_gmxe_oc_aging.o : messy_gmxe_oc_aging.f90 messy_gmxe_mem.o
messy_gmxe_soa.o : messy_gmxe_soa.f90 messy_main_constants_mem.o
messy_main_blather.o : messy_main_blather.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
messy_main_tools_kinetics.o : messy_main_tools_kinetics.f90 messy_main_constants_mem.o
# ----------------------------------------------
