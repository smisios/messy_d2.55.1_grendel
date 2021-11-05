# -*- Makefile -*-
# ----------------------------------------------
PROG = import_ts.exe

SRCS0 = $(wildcard *.f90)
SRCS  = $(filter-out F%.f90, $(SRCS0))

OBJS := $(SRCS:.f90=.o)

MAKEFILE_INC = depend.mk
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

IMPORT_TS_LIB     = $(CDI_LIB) $(NETCDF_LIB) $(PNETCDF_LIB)
IMPORT_TS_INCLUDE = $(CDI_INCLUDE) $(NETCDF_INCLUDE) $(PNETCDF_INCLUDE)

ifeq ($(strip $(POSIX90_DEF)), HAVE_POSIX90)
IMPORT_TS_INCLUDE += $(MODOPT) ../../../libsrc/posix90/src
IMPORT_TS_LIB     += -L ../../../lib -lposix90
endif

all: $(PROG)

.PHONY:depend

# update file dependencies
depend $(MAKEFILE_INC): 
	$(F_makedepend) $(SRCS)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*.f
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f *.log
	rm -f *.dat
	rm -f *~

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

.PHONY: list
list:
	@echo '------------------------------------------------------'
	@echo "SRCS = $(SRCS)"
	@echo '------------------------------------------------------'

$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(F90DEFS) -o $@ $(OBJS) $(IMPORT_TS_LIB)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(F90DEFS) $(IMPORT_TS_INCLUDE) -c $<

# ----------------------------------------------
