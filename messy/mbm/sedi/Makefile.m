# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

#F90      = @F90@
#F90FLAGS = @F90FLAGS@

LIBS      = $(NETCDF_LIB)
INCLUDES  = $(NETCDF_INCLUDE)

INSTALLDIR = ../../../bin/.
# ----------------------------------------------
include main.mk

install: all

# ----------------------------------------------
