# -*- Makefile -*-

##############################################################################

# F90      = 
# F90FLAGS = 
F90  = gfortran
F90FLAGS = -cpp -std=f2003 -Wall -O0 -g -fno-backslash -fbounds-check -fimplicit-none
#F90  = ifort
#F90FLAGS = -fpp -O0 -g 

### name of the executable that will be produced
PROG = nan.exe

##############################################################################

include main.mk

install: all
	cp -pf $(PROG) ../../../bin/.

# list of dependencies (via USE statements)
include depend.mk
