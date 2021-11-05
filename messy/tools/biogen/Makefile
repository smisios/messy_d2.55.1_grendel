# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

### FOR GRAND:
F90       = lf95
F90FLAGS  = --f95 --chk a,e,s,u -Cpp
LIBS      = -L/soft/netcdf/v3.6.0b6_lf/lib -lnetcdf
INCLUDES  = -I/soft/netcdf/v3.6.0b6_lf/include

### FOR COMPAQ/ALPHA TRUE64 OSF1:
#F90       = f90
#LIBS      = -L/soft/netcdf/v3.5.0/lib -lnetcdf
#INCLUDES  = -I/soft/netcdf/v3.5.0/include
###
###     Also, because of the many large arrays, the datasize and stacksize
###     have to be increased (tcsh):
###	- limit datasize 1000000
###     - limit stacksize 200000
###     - limit memoryuse unlimited

INSTALLDIR = .
# ----------------------------------------------
include main.mk
# ----------------------------------------------
