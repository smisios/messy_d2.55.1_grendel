# Makefile syntax info:
# When a line starts with '@', the echoing of that line is suppressed.
# When a line starts with '-', errors in the command line are ignored.

### name of the executable that will be produced
PROG = $(INSTALLDIR)/dissoc.exe

# complete list of all f90 source files (alphabetic order)
SRCS0 =	$(wildcard *.f90)
SRCS  = $(filter-out F%.f90, $(SRCS0))

# the object files are the same as the source files but with suffix ".o"
OBJS := $(SRCS:.f90=.o)

MAKEFILE_INC = depend.mk

# If you don't have sfmakedepend.pl, get it from:
# http://people.arsc.edu/~kate/Perl
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

.PHONY: all
all: $(PROG)

# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) $(NETCDF_LIB) -o $@

# update file dependencies
depend $(MAKEFILE_INC): $(SRCS)
	$(F_makedepend) $(SRCS)

# list the configuration:
.PHONY: list
list:
	@echo "------------------------------------------------"
	@echo "SRCS           = $(SRCS)"
	@echo "------------------------------------------------"
	@echo "OBJS           = $(OBJS)"
	@echo "------------------------------------------------"
	@echo "OUTPUT         = $(OUTPUT)"
	@echo "SYSTEM         = $(SYSTEM)"
	@echo "HOST           = $(HOST)"
	@echo "COMPILER       = $(COMPILER)"
	@echo "F90            = $(F90)"
	@echo "BITS           = $(BITS)"
	@echo "F90FLAGS       = $(F90FLAGS)"
	@echo "NETCDF_INCLUDE = $(NETCDF_INCLUDE)"
	@echo "NETCDF_LIB     = $(NETCDF_LIB)"
	@echo "------------------------------------------------"

.PHONY: install
install: all

.PHONY: clean
clean:
	-rm -f depend.mk.old $(OBJS) *.mod *.log *~

.PHONY: distclean
distclean: clean
	-rm -f $(PROG)
	-rm -f depend.mk*
	-rm -f *.nc
	-rm -f *.dat
	-rm -f fort.*
	-rm -f *.exe

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(F90) $(F90FLAGS) $(NETCDF_INCLUDE) -c $<
