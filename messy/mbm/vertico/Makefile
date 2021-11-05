# -*- Makefile -*-
##############################################################################

SHELL = /bin/sh

### choose Fortran90 compiler and options

SYSTEM := $(shell uname)

ifeq ($(SYSTEM),OSF1)
  ### f90 on alphas
  F90      = f90
  F90FLAGS = -O -cpp -D_alpha
  DEFOPT   = -D
  MODOPT = -I
  EXTMODCASE = asis
  #
  NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.0b6/include
  NETCDF_LIB      = -L/soft/netcdf/v3.6.0b6/lib -lnetcdff
  #
else
  ifeq ($(SYSTEM),Linux)

    ### lf95 = Lahey
    ### -Cpp  = run C-preprocessor before compiling
    ###         (necessary because of some compiler directives in the code)
    ### --pca = protect constant argument (necessary for unknown reason)
    #F90FLAGS = -Cpp --pca
    #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0
    #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0 -g --trap --verbose
    ############
    # INTEL    #
    ############
    #F90 = mpiifort    
    # or add -DNOMPI option in the preprocessor option
    F90 = ifort
    F90FLAGS = -cpp -g -autodouble -debug full -traceback -O1 -D__INTEL_v11 -DNOMPI #-Wa
    NETCDF_INCLUDE  = -I/gpfs/airsat/airfs/software/x86_64-linux/netcdf/v4.3.0_i/include
    NETCDF_LIB      = -L/gpfs/airsat/airfs/software/x86_64-linux/netcdf/v4.3.0_i/lib -lnetcdff
    ### ############
    ### # GFORTRAN #
    ### ############
    #F90 = mpif90 
    # or add -DNOMPI option in the preprocessor option
    #F90 = gfortran
    ###### for debugging:
    #F90FLAGS = -cpp -fno-second-underscore -ffree-line-length-none -fno-range-check -fbacktrace -fbounds-check -g -O0 -Wall -fdump-core -DNOMPI
    #NETCDF_INCLUDE  = -I/gpfs/airsat/airfs/software/x86_64-linux/netcdf/v4.3.0_gf/include
    #NETCDF_LIB      = -L/gpfs/airsat/airfs/software/x86_64-linux/netcdf/v4.3.0_gf/lib -lnetcdff
    ### fast, optimized model run:
    ### # F90FLAGS = -cpp -fno-second-underscore -ffree-line-length-none
    ################
    DEFOPT = -D
    MODOPT = -I
    EXTMODCASE = asis
    ### #
    ### #NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.2_lf/include
    ### #NETCDF_LIB      = -L/soft/netcdf/v3.6.2_lf/lib -lnetcdf
    ### #NETCDF_INCLUDE  = -I/athome/joec_pa/_data/soft/netcdf/v3.6.3_lf62d/include
    ### #NETCDF_LIB      = -L/athome/joec_pa/_data/soft/netcdf/v3.6.3_lf62d/lib -lnetcdf
    ### NETCDF_INCLUDE  = -I/gpfs/airsat/airfs/software/x86_64-linux/netcdf/v4.2_gf_for/include
    ### NETCDF_LIB      = -L/gpfs/airsat/airfs/software/x86_64-linux/netcdf/v4.2_gf_for/lib -lnetcdff
    ### # FORCHECK
    ### #FCKLIBS  = /data/joec_pa/soft/forcheck/v14/lib/netcdf90
    ### #FCKLIBS  = /afs/ipp/mpc/soft/@sys/forcheck/v14/netcdf90

#     ### g95 (64-bit)
#     F90 = g95_64
#     F90FLAGS = -fno-second-underscore -cpp -O0 -fbounds-check -ftrace=full
#     #F90FLAGS = -fno-second-underscore -cpp
#     DEFOPT = -D
#     MODOPT = -I
#     EXTMODCASE = asis
#     #
#     NETCDF_INCLUDE  = -I/soft/netcdf_64/v4.0_g95_64/include
#     NETCDF_LIB      = -L/soft/netcdf_64/v4.0_g95_64/lib -lnetcdf

  else
	ERROR
  endif
endif

# activate following lines for intel compiler (8.0.039) to overwrite
# previous settings
#F90      = ifort
#F90FLAGS = -cpp -O0
#DEFOPT   = -D
#MODOPT = -I
#EXTMODCASE = asis

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
