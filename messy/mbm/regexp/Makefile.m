# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .
PROG       = regexp.exe

# ----------------------------------------------
SRCS =  messy_main_constants_mem.f90 \
	messy_main_tools.f90 \
	messy_main_blather.f90 \
	regexp.f90

ifeq ($(strip $(POSIX90_DEF)), HAVE_POSIX90)
REGEXP_INCLUDE += $(MODOPT) ../../../libsrc/posix90/src
REGEXP_LIB     += -L ../../../lib -lposix90
endif

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
	$(F90) $(F90FLAGS) $(F90DEFS) -o $@ $(OBJS) $(REGEXP_LIB)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(F90DEFS) $(REGEXP_INCLUDE) $(DEFS) -c $<

# ----------------------------------------------
messy_main_blather.o : messy_main_blather.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
regexp.o : regexp.f90 messy_main_tools.o
# ----------------------------------------------
