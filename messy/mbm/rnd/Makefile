##############################################################################

### directory for executable:
INSTALLDIR = .

##############################################################################

### Try to find the compiler for current machine automatically based on
### operating system:

SYSTEM := $(shell uname)
HOST := $(shell hostname)

DEFOPT = -D

### 32 or 64-bit architecture?
ifeq "$(findstring 64,$(MACHTYPE))" "64"
  BITS = 64
else
  BITS = 32
endif

ifeq ($(SYSTEM),Linux)
  #COMPILER = lahey
  COMPILER = g95
endif

ifeq ($(SYSTEM),OSF1)
  COMPILER = f90_alpha
endif

### If this automatic mechanim doesn't work for you, define the compiler here:
#COMPILER = 

### If this doesn't work for you either, enter the correct values for
### F90, F90FLAGS, etc. in the block near the end of this file

##############################################################################

ifeq ($(COMPILER),lahey)
  ### Lahey compiler:
  F90        = lf95
  F90R8      = --dbl
  ### Choose compiler options:
  ### For lf95-6.2a, --staticlink must not be used!
  #F90FLAGS  = -Cpp --pca --ap
  F90FLAGS   = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap
  ### --staticlink is necessary for lf95-6.2d because then the
  ### undefined reference to `__divdi3' will be found in
  ### /usr/local/lf9562/lib/libfst.a
  #F90FLAGS   = -Cpp --staticlink --pca --ap
  #F90FLAGS  = -Cpp --staticlink --chk a,e,s,u --pca --ap -O0 -g --trap
endif

# ----------------------------------------------------------------------------

ifeq ($(COMPILER),g95)
  ### www.g95.org:
  #F90  = g95
  ### F90 is defined below, depending on 32/64 bit and version 0.9/0.91

  ### Choose compiler options:
  ### -Wall           = Enable most warning messages
  ### -cpp            = run the C preprocessor
  ### -O<n>           = optimization level <n>=0...3, zero = no optimization
  ### -O              = same as -O1
  ### -g              = debug option
  ### -fbounds-check  = check array and substring bounds at runtime
  ### -fimplicit-none = IMPLICIT NONE, unless overridden by IMPLICIT statements
  ### -ftrace=full    = show line number when arithmetic exceptions occur
  #F90FLAGS  = -cpp
  #F90FLAGS = -cpp -O0 -pg -fbounds-check -fimplicit-none -Wall -ftrace=full
  F90FLAGS = -cpp -O0 -pg -fbounds-check -fimplicit-none -ftrace=full

  F90R8      = -r8
  ifeq ($(BITS),64)
    ### g95 version 0.9:
    #F90  = /soft/g95_64/v0.9/bin/g95
    #NETCDF_INCLUDE = -I/soft/netcdf_64/v3.6.2_g_64/include
    #NETCDF_LIB     = -L/soft/netcdf_64/v3.6.2_g_64/lib -lnetcdf
    ### g95 version 0.91:
    F90  = /soft/g95_64/v0.91/bin/g95
    NETCDF_INCLUDE = -I/soft/netcdf_64/v3.6.2_g95v091_64/include
    NETCDF_LIB     = -L/soft/netcdf_64/v3.6.2_g95v091_64/lib -lnetcdf
  else
    ### g95 version 0.91:
    #F90  = /soft/g95/v0.91/bin/i686-suse-linux-gnu-g95
    #NETCDF_INCLUDE = -I/soft/netcdf/v3.6.2_g95v091/include
    #NETCDF_LIB     = -L/soft/netcdf/v3.6.2_g95v091/lib -lnetcdf
    ### g95 version 0.92:
    F90  = /soft/g95/v0.92/bin/i686-unknown-linux-gnu-g95
    NETCDF_INCLUDE = -I/soft/netcdf/v4.1.1_g95v092/include
    NETCDF_LIB     = -L/soft/netcdf/v4.1.1_g95v092/lib -lnetcdf
    ### g95 version 0.93: (netcdf doesn't work)
    #F90  = /soft/g95/v0.93/bin/i686-pc-linux-gnu-g95
    #NETCDF_INCLUDE = -I/soft/netcdf/v4.1.2_g95v093/include
    #NETCDF_LIB     = -L/soft/netcdf/v4.1.2_g95v093/lib -lnetcdf
  endif
endif

# ----------------------------------------------------------------------------

ifeq ($(COMPILER),f90_alpha)
  ### UNIX alpha:
  F90        = f90
  F90FLAGS   = -O -cpp
  F90R8      = -r8
endif

##############################################################################

### To implement the boxmodel on a new machine, you can simply overwrite
### the above definitions by activating the following block and entering
### appropriate values for your system:
### F90:       Fortran90 compiler
### F90FLAGS:  Fortran90 compiler options

#F90       = 
#F90FLAGS  = 

##############################################################################

# targets
include main.mk

# list of dependencies (via USE statements)
include depend.mk

##############################################################################
