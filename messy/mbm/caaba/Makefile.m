# -*- makefile-gmake -*-

##############################################################################

srcdir     = .
PROG       = caaba.exe

INSTALLDIR  = ../../../bin

##############################################################################

OUTPUT = NETCDF

ADDEFS = E4CHEM

# F90NOR8 is needed for specific.mk:
F90NOR8 = $(F90FLAGS)

include main.mk

install: all

# list of dependencies (via USE statements)
include depend.mk
-include ../../smcl/specific.mk
