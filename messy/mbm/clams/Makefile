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


TMPF90 := $(shell mpif90 -show)
TMPF90 := $(word 1, $(TMPF90))


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

    # JURECA:
    ifeq ($(NODENAME),jr)

	ifeq ($(TMPF90),gfortran)
           F90    = mpif90
           F90FLAGS = -cpp -O2 -fno-second-underscore -ffree-line-length-none -fno-range-check
           DEFOPT = -D
           MODOPT = -I
           EXTMODCASE = asis
           NETCDF_INCLUDE  = -I $(EBROOTNETCDFMINFORTRAN)/include
           NETCDF_LIB      = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
#           QHULLLIB = -L/homec/jicg11/jicg1108/clams_lib -lqhull
	else	
           F90    = mpif90
           F90FLAGS = -fpp -O2 -qopt-report-phase=vec
#           F90FLAGS = -fpp -O2 -sox -cpp -D__INTEL_v11 -fp-model strict -align all 
	   DEFOPT = -D
           MODOPT = -I
           EXTMODCASE = asis
           NETCDF_INCLUDE  = 
           NETCDF_LIB      = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
#           QHULLLIB = -L/homec/jicg11/jicg1108/clams_lib -lqhull
	endif


    # IEK-7 Cluster:
    else	

       # F90 = ifort
       F90 = mpif90
       ### -fpp  = run preprocessor before compiling
       ###         (necessary because of some compiler directives in the code)
       F90FLAGS = -fpp -O2
       #F90FLAGS = -fpp -g -check bounds -check all -fpe0 -traceback #-warn unused -pg -zero -CB  
       DEFOPT = -D
       MODOPT = -I
       EXTMODCASE = asis

       # get version number of ifc
       version := $(shell ifort --version | grep ifort)
       version := $(word 3, $(version))

       NETCDF_INCLUDE  = -I/usr/nfs/software/icg1/mod/mod_linux_ifc64_$(version)
       NETCDF_LIB      = -L/usr/nfs/software/icg1/lib/lib_linux_ifc64_$(version) -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
#       QHULLLIB = -L/usr/nfs/software/icg1/lib/lib_linux_ifc64_$(version) -lqhull

   endif

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
