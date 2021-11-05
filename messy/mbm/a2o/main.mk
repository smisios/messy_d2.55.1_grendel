### name of the executable that will be produced
PROG = $(INSTALLDIR)/a2o.exe

# complete list of all f90 source files
SRCS0 = $(wildcard *.f90)
SRCS  = $(filter-out F%.f90, $(SRCS0))

# the object files are the same as the source files but with suffix ".o"
OBJS := $(SRCS:.f90=.o)

MAKEFILE_INC = depend.mk

# If you don't have the perl script sfmakedepend, get it from:
# http://people.arsc.edu/~kate/Perl
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

all: $(PROG)

# the dependencies depend on the link
# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) $(NETCDF_LIB) -o $@

# update file dependencies
depend $(MAKEFILE_INC): $(SRCS)
	$(F_makedepend) $(SRCS)


# forcheck
.PHONY: check
check:
	../../util/mfchk a2o "_A2O" "." "." "$(F90)" "$(FCKLIBS)"

# put the same find commands here as in cmg command:
.PHONY: TAGS
TAGS:
	@F90FILES=`find . -name "*.f90" -type f` ;\
	 INCFILES=`find . -name "*.inc" -type f` ;\
	 NMLFILES=`find . -name "*.nml" -type f` ;\
	 XSCRIPTS=`find . -name "x*" -type f -perm -100` ;\
	 etags -l fortran $$F90FILES $$INCFILES $$NMLFILES \
           -lnone $$XSCRIPTS

# check files
list:
	@echo "SRCS           = $(SRCS)"
	@echo "OUTPUT         = $(OUTPUT)"
	@echo "SYSTEM         = $(SYSTEM)"
	@echo "HOST           = $(HOST)"
	@echo "COMPILER       = $(COMPILER)"
	@echo "BITS           = $(BITS)"
	@echo "F90FLAGS       = $(F90FLAGS)"
	@echo "NETCDF_INCLUDE = $(NETCDF_INCLUDE)"
	@echo "NETCDF_LIB     = $(NETCDF_LIB)"

# this command is executed by "gmake clean"
clean:
	rm -f depend.mk.old *.o *.mod *.log *~

distclean: clean
	rm -f $(PROG)
	rm -f depend.mk* 
	rm -f *.nc
	rm -f *.dat
	rm -f *.exe

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(F90) $(DEFOPT)$(MPI) $(DEFOPT)_A2O $(F90FLAGS) $(NETCDF_INCLUDE) -c $<
