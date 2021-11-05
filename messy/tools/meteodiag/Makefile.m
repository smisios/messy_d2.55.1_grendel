# -*- Makefile -*-

##############################################################################

srcdir     = .
PROG       = meteodiag.exe

# overwrite compiler with non-MPI wrapped
F90 := $(shell (. ../../util/locate_f90.sh; echo $$F90))

INSTALLDIR  = ../../../bin

INCLUDES   = $(NETCDF_INCLUDE)
LIBS       = $(NETCDF_LIB)

##############################################################################

#ADDEFS = $(DEFOPT)E4CHEM

include main.mk

install: all

# list of dependencies (via USE statements)
include depend.mk
