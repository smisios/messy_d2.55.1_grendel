# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .
PROG       = msbm.exe

# ----------------------------------------------
SRCS =  messy_main_constants_mem.f90 \
 messy_main_tools.f90 \
 messy_main_blather.f90 \
 messy_msbm.f90 \
 messy_msbm_box.f90 \
 msbm.f90


OBJS := $(SRCS:.f90=.o)

all: $(PROG)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*~
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod

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
messy_msbm_box.o : messy_msbm_box.f90 messy_msbm.o
messy_msbm.o : messy_msbm.f90 messy_main_tools.o messy_main_constants_mem.o
msbm.o : msbm.f90 messy_msbm_box.o
# ----------------------------------------------
