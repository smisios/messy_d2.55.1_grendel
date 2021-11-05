# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .

## overwrite compiler with non-MPI wrapped
#F90 := $(shell (. ../../util/locate_f90.sh; echo $$F90))

include main.mk

include depend.mk
# ----------------------------------------------
