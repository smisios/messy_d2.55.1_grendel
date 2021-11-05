# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

INSTALLDIR = ../../../bin

all: $(PROG)

install: all

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) -o $@ $(OBJS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

include main.mk

include depend.mk
