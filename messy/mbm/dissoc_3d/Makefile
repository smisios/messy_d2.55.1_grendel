# -*- Makefile -*-
##############################################################################

SHELL = /bin/sh

### choose Fortran90 compiler and options

SYSTEM := $(shell uname)
NODENAME := $(shell uname -n)
NODENAME := $(subst ., ,$(NODENAME))
NODENAME := $(word 1, $(NODENAME))
NODENAME := $(subst jr,jr ,$(NODENAME))
NODENAME := $(word 1, $(NODENAME))

ifeq ($(SYSTEM),OSF1)
  ### f90 on alphas
  F90      = f90
  F90FLAGS = -O -cpp -D_alpha
  DEFOPT   = -D
  MODOPT = -I
  EXTMODCASE = asis
  #
  NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.0b6/include
  NETCDF_LIB      = -L/soft/netcdf/v3.6.0b6/lib -lnetcdf
  #
else
  ifeq ($(SYSTEM),Linux)

    # F90 = ifort
    F90 = mpif90
    ### -fpp2  = run preprocessor before compiling
    ###         (necessary because of some compiler directives in the code)
    F90FLAGS =  -fpp -O3 #-check all -fpe0 -traceback #-warn unused  #-g #-pg -zero #-CB  #k #-O0 -g
    DEFOPT = -D
    MODOPT = -I
    EXTMODCASE = asis

    ### NetCDF-Library:
    # JURECA:
    ifeq ($(NODENAME),jr)
	NETCDF_INCLUDE  = 
	NETCDF_LIB      = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
    # IEK-7 Cluster:
    else	
       # get version number of ifc
       version := $(shell ifort --version | grep ifort)
       version := $(word 3, $(version))
       NETCDF_INCLUDE  = -I/usr/nfs/software/icg1/mod/mod_linux_ifc64_$(version)
       NETCDF_LIB      = -L/usr/nfs/software/icg1/lib/lib_linux_ifc64_$(version) -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
   endif

    ASAD_INCLUDE = -I../include_asad
    CHEM_INCLUDE = -I../include_chem

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
