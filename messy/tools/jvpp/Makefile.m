# -*- makefile-gmake -*-

##############################################################################

srcdir     = .
PROG       = jvpp.exe

# overwrite compiler with non-MPI wrapped
F90 := $(shell (. ../../util/locate_f90.sh; echo $$F90))

INSTALLDIR  = ../../../bin

##############################################################################

include main.mk

install: all

# list of dependencies (via USE statements)
include depend.mk
