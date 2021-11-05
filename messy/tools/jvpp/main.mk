### name of the executable that will be produced
PROG       = $(INSTALLDIR)/jvpp.exe

# complete list of all f90 source files (alphabetic order)
SRCS0 =	$(wildcard *.f90)
SRCS  = $(filter-out F%.f90, $(SRCS0))

# the object files are the same as the source files but with suffix ".o"
OBJS := $(SRCS:.f90=.o)

MAKEFILE_INC = depend.mk

# If you don't have sfmakedepend.pl, get it from:
# http://people.arsc.edu/~kate/Perl
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

all: $(PROG)

# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) -o $@

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
	@echo "HOST           = $(HOST)"
	@echo "COMPILER       = $(COMPILER)"
	@echo "F90            = $(F90)"
	@echo "BITS           = $(BITS)"
	@echo "F90FLAGS       = $(F90FLAGS)"
	@echo "------------------------------------------------"

# this command is executed by "gmake clean"
clean:
	rm -f $(OBJS) *.mod *.old *~

distclean: clean
	rm -f $(PROG)
	rm -f depend.mk*
	rm -f workdir_176/* workdir_eff/* workdir_f90/* workdir_jnl/*

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(F90) $(F90FLAGS) -c $<
