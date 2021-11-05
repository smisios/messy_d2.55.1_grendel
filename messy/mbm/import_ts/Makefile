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

    ### lf95 = Lahey
    F90 = lfc
    ### -Cpp  = run C-preprocessor before compiling
    ###         (necessary because of some compiler directives in the code)
    ### --pca = protect constant argument (necessary for unknown reason)
    #F90FLAGS = -Cpp --pca
    #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0
    #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0 -g --trap --verbose
    F90FLAGS = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap
    DEFOPT = -D
    #
    #NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.2_lf/include
    #NETCDF_LIB      = -L/soft/netcdf/v3.6.2_lf/lib -lnetcdf
#    NETCDF_INCLUDE  = -I/athome/joec_pa/_data/soft/netcdf/v3.6.3_lf62d/include
#    NETCDF_LIB      = -L/athome/joec_pa/_data/soft/netcdf/v3.6.3_lf62d/lib -lnetcdf
    NETCDF_INCLUDE  = -I/afs/ipp/mpc/soft/@sys/netcdf/v3.6.2/include
    NETCDF_LIB      = -L/afs/ipp/mpc/soft/@sys/netcdf/v3.6.2/lib -lnetcdf

#     ### g95 (64-bit)
#     F90 = g95_64
#     F90FLAGS = -fno-second-underscore -cpp -O0 -fbounds-check -ftrace=full
#     #F90FLAGS = -fno-second-underscore -cpp
#     DEFOPT = -D
#     #
#     NETCDF_INCLUDE  = -I/soft/netcdf_64/v4.0_g95_64/include
#     NETCDF_LIB      = -L/soft/netcdf_64/v4.0_g95_64/lib -lnetcdf
  else ifeq ($(OSTYPE),aix)

    F90 = xlf95_r
    ### -Cpp  = run C-preprocessor before compiling
    ###         (necessary because of some compiler directives in the code)
    F90FLAGS = -q64 -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -O3 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qfree=f90 -qspillsize=32648 -qarch=pwr4 -qtune=pwr4 -Q -qhot
    DEFOPT = -WF,-D
    NETCDF_LIB      = -L/pf/b/b302019/soft/netcdf/v3.6.3/lib -lnetcdf
    NETCDF_INCLUDE  = -I/pf/b/b302019/soft/netcdf/v3.6.3/include -I/pf/b/b302019/soft/netcdf/v3.6.3/include

  else
	ERROR
  endif
endif

# activate following lines for intel compiler (8.0.039) to overwrite
# previous settings
#F90      = ifort
#F90FLAGS = -cpp -O0
#DEFOPT   = -D

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
F90DEFS0 = 
F90DEFS  = $(addprefix $(DEFOPT), $(F90DEFS0))

include main.mk

include depend.mk
