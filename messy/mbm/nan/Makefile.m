# -*- Makefile -*-

##############################################################################

### F90      = 
### F90FLAGS =

### name of the executable that will be produced
PROG = ../../../bin/nan.exe

##############################################################################

include main.mk

install: all
#	cp -pf $(PROG) ../../../bin/.

# list of dependencies (via USE statements)
include depend.mk
