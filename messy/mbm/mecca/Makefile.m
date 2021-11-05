# -*- makefile-gmake -*-
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .

# F90NOR8 is needed for specific.mk:
F90NOR8 = $(F90FLAGS)

include main.mk

include depend.mk
-include ../../smcl/specific.mk
# ----------------------------------------------
