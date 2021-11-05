### name of the executable that will be produced
PROG = $(INSTALLDIR)/meteodiag.exe

# complete list of all f90 source files
SRCS0 = $(wildcard *.f90)
SRCS  = $(filter-out F%.f90, $(SRCS0))

# the object files are the same as the source files but with suffix ".o"
OBJS := $(SRCS:.f90=.o)

F90FLAGS += ${DEFOPT}NOMPI

MAKEFILE_INC = depend.mk

# If you don't have the perl script sfmakedepend, get it from:
# http://people.arsc.edu/~kate/Perl
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

all: $(PROG)

# the dependencies depend on the link
# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) $(LIBS) -o $@

# update file dependencies
depend $(MAKEFILE_INC): $(SRCS)
	$(F_makedepend) $(SRCS)

.PHONY: zip
zip:
	./zipmeteodiag.tcsh zip

.PHONY: zipall
zipall:
	./zipmeteodiag.tcsh zipall

# check files
list:
	@echo "SRCS = $(SRCS)"

# this command is executed by "gmake clean"
clean:
	rm -f $(OBJS) *.mod *.log F*.f *.i90 *.optrpt *~ rm -f $(PROG) 

distclean: clean
	rm -f $(PROG)
	rm -f depend.mk* 
	rm -f *.nc
	rm -f *.dat

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) -c $<
