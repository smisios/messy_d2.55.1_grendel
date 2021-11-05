# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .
PROG       = random.exe

# ----------------------------------------------
SRCS =  messy_main_constants_mem.f90 \
	messy_main_tools.f90 \
        messy_main_blather.f90 \
        messy_main_rnd.f90 \
        messy_main_rnd_lux.f90 \
        messy_main_rnd_mtw.f90 \
        messy_main_rnd_mtw_ja.f90 \
	random.f90


OBJS := $(SRCS:.f90=.o)

all: $(PROG)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f $(srcdir)/*.dat

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
messy_main_rnd.o : messy_main_rnd.f90 messy_main_rnd_lux.o messy_main_rnd_mtw.o messy_main_rnd_mtw_ja.o messy_main_tools.o messy_main_constants_mem.o
messy_main_rnd_lux.o : messy_main_rnd_lux.f90
messy_main_rnd_mtw.o : messy_main_rnd_mtw.f90
messy_main_rnd_mtw_ja.o : messy_main_rnd_mtw_ja.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
random.o : random.f90 messy_main_tools.o messy_main_rnd.o messy_main_constants_mem.o
# ----------------------------------------------
