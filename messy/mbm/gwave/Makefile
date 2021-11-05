##############################################################################

INSTALLDIR = .
##############################################################################

### Try to find the compiler for current machine automatically based on
### operating system:

SYSTEM := $(shell uname)

ifeq ($(SYSTEM),OSF1)
  compiler = f90_alpha
endif

ifeq ($(SYSTEM),Linux)
  compiler = lahey
  #compiler = intel
  #compiler = g95
  MODOPT = -I
  DEFOPT   = -D
endif

ifeq ($(SYSTEM),AIX)
  compiler = xlf95
   DEFOPT = -WF,-D
   MODOPT = -I
   EXTMODCASE = asis
endif

### If this automatic mechanim doesn't work for you, define the compiler here:
#compiler = 

### If this doesn't work for you either, enter the correct values for
### F90, F90FLAGS, etc. in the block near the end of this file

##############################################################################

ifeq ($(compiler),lahey)
  ### lf95 = Lahey compiler
  ###
  ### --staticlink is necessary since lf95-6.2d because then the
  ### undefined reference to `__divdi3' will be found in
  ### /usr/local/lf9562/lib/libfst.a
  ### However, for lf95-6.2a, --staticlink must not be used!
  F90        = lf95
  #F90FLAGS  = -Cpp --pca
  F90FLAGS   = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap
  #F90FLAGS  = -Cpp --staticlink --pca
  #F90FLAGS  = -Cpp --staticlink --chk a,e,s,u --pca --ap -O0 -g --trap
  NETCDF_INCLUDE   = -I/soft/netcdf/v3.6.0b6_lf/include
  NETCDF_LIB       = -L/soft/netcdf/v3.6.0b6_lf/lib -lnetcdf
  #NETCDF_INCLUDE = -I/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.2/include
  #NETCDF_LIB     = -L/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.2/lib -lnetcdf
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),f90_alpha)
  ### UNIX alpha
  F90        = f90
  F90FLAGS   = -O -cpp 
  INCLUDES   = -I/soft/netcdf/v3.6.0b6/include
  LIBS       = -L/soft/netcdf/v3.6.0b6/lib -lnetcdf
  #INCLUDES  = -I/usr/local/netcdf/v3.6.0b6/include
  #LIBS      = -L/usr/local/netcdf/v3.6.0b6/lib -lnetcdf
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),intel)
  ### intel compiler (8.0.039):
  F90        = ifort
  F90FLAGS   = -cpp -O0
  NETCDF_INCLUDE = -I/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.3/include
  NETCDF_LIB     = -L/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.3/lib -lnetcdf
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),g95)
  ### gnu compiler www.g95.org
  F90        = g95
  F90FLAGS   = -cpp -O0
  INCLUDES   = -I/soft/netcdf/v3.6.0b6_g/include
  LIBS       = -L/soft/netcdf/v3.6.0b6_g/lib -lnetcdf
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),xlf95)
  # AIX
  F90       = xlf95_r
#  F90FLAGS  = -O3 -q64 -g -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qfree=f90 -qspillsize=32648 -qarch=pwr6 -qtune=pwr6 -Q -qhot
  F90FLAGS  = -O0 -qcheck -qfullpath -q64 -g -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qfree=f90 -qspillsize=32648 
  NETCDF_INCLUDE = -I/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.3/include
  NETCDF_LIB     = -L/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.3/lib -lnetcdf
endif

##############################################################################

### To implement the boxmodel on a new machine, you can simply overwrite
### the above definitions by activating the following block and entering
### appropriate values for your system:
### F90:       Fortran90 compiler
### F90FLAGS:  Fortran90 compiler options
### INCLUDES:  path to include files (for netcdf)
### LIBS:      path to library files and library name (for netcdf)

#F90       = 
#F90FLAGS  = 
#INCLUDES  = 
#LIBS      = 

##############################################################################


##############################################################################

# targets
include main.mk

##############################################################################
