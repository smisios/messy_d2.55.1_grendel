# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .
PROG       = scav.exe

# ----------------------------------------------
SRCS =  messy_main_constants_mem.f90 \
 messy_main_tools.f90 \
 messy_main_blather.f90 \
 messy_scav_mem.f90 \
 messy_main_tools_kinetics.f90 \
 messy_scav_inp_kpp.f90 \
 messy_scav_liq.f90 \
 messy_scav_ice.f90 \
 messy_scav_aer.f90 \
 messy_scav_l_kpp.f90 \
 messy_scav_i_kpp.f90 \
 messy_scav.f90 \
 messy_scav_inter.f90 \
 messy_scav_box.f90

OBJS := $(SRCS:.f90=.o)

all: $(PROG)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*~
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f $(srcdir)/fort.*

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) -o $@ $(OBJS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

# ----------------------------------------------
messy_main_blather.o : messy_main_blather.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
messy_main_tools_kinetics.o : messy_main_tools_kinetics.f90 messy_main_constants_mem.o
messy_scav_aer.o : messy_scav_aer.f90 messy_scav_mem.o messy_main_constants_mem.o
messy_scav_box.o : messy_scav_box.f90 messy_main_tools.o messy_scav_liq.o messy_scav_inter.o messy_scav_mem.o messy_scav.o messy_scav_l_kpp.o messy_main_constants_mem.o
messy_scav.o : messy_scav.f90 messy_scav_inp_kpp.o messy_scav_ice.o messy_main_constants_mem.o messy_main_tools.o messy_scav_i_kpp.o messy_scav_liq.o messy_scav_aer.o messy_scav_mem.o messy_scav_l_kpp.o
messy_scav_ice.o : messy_scav_ice.f90 messy_scav_liq.o messy_scav_aer.o messy_main_constants_mem.o messy_scav_inp_kpp.o messy_scav_mem.o messy_scav_i_kpp.o
messy_scav_i_kpp.o : messy_scav_i_kpp.f90 messy_main_tools.o messy_scav_inp_kpp.o messy_main_constants_mem.o
messy_scav_inp_kpp.o : messy_scav_inp_kpp.f90 messy_main_tools_kinetics.o messy_main_constants_mem.o
messy_scav_inter.o : messy_scav_inter.f90 messy_scav_aer.o messy_scav.o messy_main_tools.o messy_main_constants_mem.o messy_scav_i_kpp.o messy_scav_l_kpp.o messy_scav_mem.o
messy_scav_liq.o : messy_scav_liq.f90 messy_main_constants_mem.o messy_scav_aer.o messy_main_tools_kinetics.o messy_scav_inp_kpp.o messy_scav_mem.o messy_scav_l_kpp.o
messy_scav_l_kpp.o : messy_scav_l_kpp.f90 messy_main_tools.o messy_scav_inp_kpp.o messy_main_constants_mem.o
messy_scav_mem.o : messy_scav_mem.f90 messy_main_constants_mem.o messy_scav_l_kpp.o messy_scav_i_kpp.o
# ----------------------------------------------
