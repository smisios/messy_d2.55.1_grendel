##############################################################################
INSTALLDIR = .

### Try to find compiler settings for current machine automatically based on
### operating system and host name:

SYSTEM := $(shell uname)

ifeq ($(SYSTEM),OSF1)
  compiler = f90_alpha
endif

ifeq ($(SYSTEM),Linux)
  compiler = lahey
endif

ifeq ($(SYSTEM),AIX)
  compiler = xlf95
endif

ifeq ($(SYSTEM),Darwin)
  compiler = xlf95
endif

### If this automatic mechanim doesn't work for you, activate one of
### the next lines:
#compiler = f90_alpha
#compiler = lahey
#compiler = xlf95
#compiler = g95

### If this doesn't work for you either, enter the correct values for
### F90, F90FLAGS, etc. in the block near the end of this file

##############################################################################

### lf95 = Lahey compiler

ifeq ($(compiler),lahey)
  F90        = lf95
  #F90FLAGS  = -Cpp --pca
  F90FLAGS   = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap --verbose
  F90R8      = --dbl
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),f90_alpha)
  F90        = f90
  F90FLAGS   = -O -cpp 
  F90R8      = -r8
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),g95)
  F90        = g95
  F90FLAGS   = -cpp -O0
  F90R8      = -r8
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),xlf95)
  # AIX
  F90      = xlf95_r
  F90FLAGS = -q64 -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -O3 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qfree=f90 -qspillsize=32648 -qarch=pwr4 -qtune=pwr4 -Q -qhot
  F90R8    = -qrealsize=8
endif

##############################################################################

### To implement the boxmodel on a new machine, you can simply overwrite
### the above definitions by activating the following block and entering
### appropriate values for your system:
### F90:       Fortran90 compiler
### F90FLAGS:  Fortran90 compiler options
### F90R8:     Fortran90 compiler option to set double precision

#F90       = 
#F90FLAGS  = 
#F90R8     =

##############################################################################

# targets
include main.mk

# list of dependencies (via USE statements)
include depend.mk
