# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

# overwrite compiler with non-MPI wrapped
F90 := $(shell (. ../../util/locate_f90.sh; echo $$F90))

LIBS      = $(NETCDF_LIB)
INCLUDES  = $(NETCDF_INCLUDE)

# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
PROG       = ncdx.exe

# ----------------------------------------------

SRCS = f2kcli.f90 ncdx.f90

OBJS := $(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJS) $(LIBS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -c $<

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.


clean:
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

# ----------------------------------------------
f2kcli.o : f2kcli.f90
ncdx.o : ncdx.f90 f2kcli.o
# ----------------------------------------------
