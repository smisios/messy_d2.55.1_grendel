# -*- Makefile -*-
# ----------------------------------------------
PROG = mecca.exe

#SRCS0 = $(wildcard *.f90)
#SRCS  = $(filter-out F%.f90, $(SRCS0))
SRCS0Q = $(wildcard *.f90)
SRCSL = $(shell find -L * -type l -print)
SRCS0 = $(filter-out $(SRCSL), $(SRCS0Q))
SRCS  = $(filter-out F%.f90, $(SRCS0))

OBJS := $(SRCS:.f90=.o)

MAKEFILE_INC = depend.mk
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

ifeq ($(strip $(POSIX90_DEF)), HAVE_POSIX90)
MECCA_INCLUDE += $(MODOPT) ../../../libsrc/posix90/src
MECCA_LIB     += -L ../../../lib -lposix90
endif

all: $(PROG)

.PHONY:depend

# update file dependencies
depend $(MAKEFILE_INC): 
	$(F_makedepend) $(SRCS)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f *.log
	rm -f *~

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

.PHONY: run
run: all 
	$(srcdir)/$(PROG)

.PHONY: list
list:
	@echo '------------------------------------------------------'
	@echo "SRCS = $(SRCS)"
	@echo '------------------------------------------------------'

$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(F90DEFS) -o $@ $(OBJS) $(MECCA_LIB)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(F90DEFS) $(MECCA_INCLUDE) -c $<

# ----------------------------------------------
