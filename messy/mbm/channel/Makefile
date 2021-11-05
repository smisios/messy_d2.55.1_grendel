# -*- Makefile -*-
##############################################################################
srcdir     = .
bindir     = .
includedir = .

### choose Fortran90 compiler and options

SYSTEM := $(shell uname)
ifeq ($(SYSTEM),OSF1)
  ### f90 on alphas
  F90      = f90
  F90FLAGS = -O -cpp -D_alpha
  DEFOPT   = -D
  #
  NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.0b6/include
  NETCDF_LIB      = -L/soft/netcdf/v3.6.0b6/lib -lnetcdf
  #
else
  ifeq ($(SYSTEM),Linux)

    F90 = gfortran
    F90FLAGS = -cpp -D__linux__ -fno-second-underscore -ffree-line-length-none -fno-range-check -O3 
    DEFOPT = -D
    #
    NETCDF_INCLUDE  = -I/export/home/joec_pa/Programs/netcdf/3.6.3/gfortran/7.3.1/include
    NETCDF_LIB      = -L/export/home/joec_pa/Programs/netcdf/3.6.3/gfortran/7.3.1/lib -lnetcdf

  else
	ERROR
  endif
endif

##############################################################################

### The above block defines the variables F90, F90FLAGS, INCLUDES, and
### LIBS. The if...else...endif constructs try to find suitable values
### for different architectures and machines. To implement the
### boxmodel on a new machine, you can simply overwrite the above
### definitions by activating the following block and entering
### appropriate values:
### F90      = 
### F90FLAGS = 

##############################################################################

include main.mk

include depend.mk
