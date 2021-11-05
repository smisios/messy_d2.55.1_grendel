# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

# overwrite compiler with non-MPI wrapped
F90 := $(shell (. ../../util/locate_f90.sh; echo $$F90))

LIBS      = $(NETCDF_LIB)
INCLUDES  = $(NETCDF_INCLUDE)

INSTALLDIR = ../../../bin/.
# ----------------------------------------------
include main.mk

install: all
# ----------------------------------------------
