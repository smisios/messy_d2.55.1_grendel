# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .
PROG       = dradon.exe

# ----------------------------------------------
SRCS =  messy_main_constants_mem.f90 \
	messy_main_tools.f90 \
	messy_main_blather.f90 \
	messy_dradon.f90 \
	messy_dradon_box.f90 \
	dradon.f90


OBJS := $(SRCS:.f90=.o)

all: $(PROG)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f $(srcdir)/*.out

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) -o $@ $(OBJS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

# ----------------------------------------------
dradon.o : dradon.f90 messy_dradon_box.o
messy_dradon_box.o : messy_dradon_box.f90 messy_main_tools.o messy_dradon.o
messy_dradon.o : messy_dradon.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_blather.o : messy_main_blather.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
# ----------------------------------------------
