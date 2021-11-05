##############################################################################
# ENTER YOUR FORTRAN 90 COMPILER SETTINGS HERE:
##############################################################################

# settings for compiler g95 from www.g95.org:
F90  = g95
F90FLAGS = -cpp -fimplicit-none -fno-second-underscore -ftrace=full -fno-backslash -fbounds-check -g -O0 -Wall

##############################################################################
# NO NEED TO CHANGE ANYTHING BELOW THIS LINE:
##############################################################################

# name of the executable that will be produced:
PROG = mecca.exe

# complete list of all f90 source files
SRCS = $(wildcard *.f90)

# the object files are the same as the source files but with suffix ".o"
OBJS := $(SRCS:.f90=.o)

##############################################################################

$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) -o $@

# Update file dependencies
depend depend.mk: $(SRCS)
	./sfmakedepend --file=depend.mk $(SRCS)
# If you don't have sfmakedepend.pl, get it from:
# http://people.arsc.edu/~kate/Perl

.PHONY: clean
clean:
	-rm -f depend.mk.old $(OBJS) *.mod *.log *~

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

# list of dependencies (via USE statements)
include depend.mk

##############################################################################
