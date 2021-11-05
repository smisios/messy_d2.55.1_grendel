##############################################################################

### select either ASCII or NETCDF output:
#OUTPUT     = ASCII
OUTPUT     = NETCDF

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
  ### Choose compiler:
  #COMPILER = lahey
  #COMPILER = intel
  COMPILER = g95
  #COMPILER = gfortran
  #COMPILER = pgf90
endif

ifeq ($(SYSTEM),OSF1)
  COMPILER = f90_alpha
endif

ifeq ($(SYSTEM),AIX)
  COMPILER = xlf95
  DEFOPT = -WF,-D
endif

### If this automatic mechanim doesn't work for you, define the compiler here:
#COMPILER = 

### If this doesn't work for you either, enter the correct values for
### F90, F90FLAGS, etc. in the block near the end of this file

##############################################################################

ifeq ($(COMPILER),lahey)
  ### Lahey compiler:
  F90        = lf95

  ### Choose compiler options:
  ### For lf95-6.2a, --staticlink must not be used!
  #F90FLAGS  = -Cpp --pca --ap
  F90FLAGS   = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap
  ### --staticlink is necessary for lf95-6.2d because then the
  ### undefined reference to `__divdi3' will be found in
  ### /usr/local/lf9562/lib/libfst.a
  #F90FLAGS   = -Cpp --staticlink --pca --ap
  #F90FLAGS  = -Cpp --staticlink --chk a,e,s,u --pca --ap -O0 -g --trap

  ifeq ($(BITS),64)
    NETCDF_INCLUDE  = -I/soft/netcdf_64/v3.6.2_lf_64/include
    NETCDF_LIB      = -L/soft/netcdf_64/v3.6.2_lf_64/lib -lnetcdf
  else
    NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.2_lf/include
    NETCDF_LIB      = -L/soft/netcdf/v3.6.2_lf/lib -lnetcdf
  endif
  #NETCDF_INCLUDE = -I/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.2/include
  #NETCDF_LIB     = -L/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.2/lib -lnetcdf

endif

# ----------------------------------------------------------------------------

ifeq ($(COMPILER),g95)
  ### www.g95.org:
  #F90  = g95
  ### F90 is defined below, depending on 32/64 bit and version

  ### Choose compiler options:
  ### -cpp                    = run the C preprocessor
  ### -O<n>                   = optimization level <n>=0...3, zero = no optimization
  ### -O                      = same as -O1
  ### -g                      = debug option
  ### -fbounds-check          = check array and substring bounds at runtime
  ### -fimplicit-none         = IMPLICIT NONE, unless overridden by IMPLICIT statements
  ### -ftrace=full            = show line number when arithmetic exceptions occur
  ### -fno-backslash          = treat backslashes literally
  ### -fno-second-underscore  = don't append a second trailing underscore in names having an underscore
  ### -ffree-line-length-huge = allow very large source lines (10k)
  ### -Wall                   = Enable most warning messages
  #F90FLAGS  = -cpp -fimplicit-none -fno-second-underscore -ftrace=full -fno-backslash -fbounds-check -g -O0
  F90FLAGS   = -cpp -fimplicit-none -fno-second-underscore -ftrace=full -fno-backslash -fbounds-check -g -O0 -Wall

  ifeq ($(BITS),64)
    ### g95 version 0.93:
    F90  = /soft/g95_64/v0.93/bin/g95
    NETCDF_INCLUDE = -I/soft/netcdf_64/v3.6.3_g95_dint32/include
    NETCDF_LIB     = -L/soft/netcdf_64/v3.6.3_g95_dint32/lib -lnetcdf
    # DLR pa1:
    # F90 = g95
    # NETCDF_INCLUDE = -I/export/opt/PA/prgs/netcdf/3.6.3/g95_dint32/0.93/include
    # NETCDF_LIB     = -L/export/opt/PA/prgs/netcdf/3.6.3/g95_dint32/0.93/lib -lnetcdf
  else
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

# from AB:
#F90  = /soft/g95/v0.93/bin/g95
#NETCDF_INCLUDE = -I/soft/netcdf/v4.1.1_g95_64/include
#NETCDF_LIB     = -L/soft/netcdf/v4.1.1_g95_64/lib -lnetcdf

# ----------------------------------------------------------------------------

ifeq ($(COMPILER),f90_alpha)
  ### UNIX alpha:
  F90             = f90
  F90FLAGS        = -O -cpp
  NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.0b6/include
  NETCDF_LIB      = -L/soft/netcdf/v3.6.0b6/lib -lnetcdf
endif

# ----------------------------------------------------------------------------

### Note that the intel compiler is neither recommended nor supported
### for JVAL. Use it at your own risk.
ifeq ($(COMPILER),intel)
  ### intel compiler (8.0.039):
  F90             = ifort
  F90FLAGS        = -cpp -O0
  # SOFT:
  NETCDF_INCLUDE  = -I/soft/netcdf/v3.5.1b11_i/include
  NETCDF_LIB      = -L/soft/netcdf/v3.5.1b11_i/lib -lnetcdf
  # DRESDEN:
  #NETCDF_INCLUDE = -I/gpfs/airsat/airfs/software/x86_64-linux/netcdf/v4.3.0_i/include
  #NETCDF_LIB     = -L/gpfs/airsat/airfs/software/x86_64-linux/netcdf/v4.3.0_i/lib -lnetcdf -lnetcdff 
  # MPC RZG:
  NETCDF_INCLUDE  = -I/mpcdata/software/x86_64-suse-linux/netcdf/v3.6.2_i/include
  NETCDF_LIB      = -L/mpcdata/software/x86_64-suse-linux/netcdf/v3.6.2_i/lib -lnetcdf 
  #NETCDF_INCLUDE  = -I/mpcdata/software/x86_64-suse-linux/netcdf/v4.0_i/include
  #NETCDF_LIB      = -L/mpcdata/software/x86_64-suse-linux/netcdf/v4.0_i/lib -lnetcdf 
endif

# ----------------------------------------------------------------------------

ifeq ($(COMPILER),gfortran)
  ### gnu compiler: http://gcc.gnu.org/fortran/index.html
  ### see compiler options: gfortran --help -v
  ### gfortran settings from Hella Riede (2013)
  ### don't use '-x f95-cpp-input' in F90FLAGS -> linker error:
  ###   not a bug, see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=29497
  F90      = gfortran
  ### for production: -O3: optimize, -cpp: enable pre-processing,
  ###   -ffree-line-length-none: no line width limit
  F90FLAGS = -O3 -cpp -ffree-line-length-none
  ###   -fdefault-real-8 = double precision
  ### for debug: -O0: no optimization, -Wall: all error messages,
  ###   -fbounds-check: check bounds before indexing arrays
  #F90FLAGS = -O0 -cpp -fbounds-check -Wall
  ### for debug: do not compile, check syntax
  #F90FLAGS = -cpp -fsyntax-only -fdump-parse-tree 
  ### other flags to try while troubleshooting:
  ###   -std=legacy, -fimplicit-none, -finit-local-zero, -fall-intrinsics,
  ###   -ffpe-trap=invalid,zero,overflow,underflow,denormal
  ifeq ($(BITS),64)
    #NETCDF_INCLUDE = -I/soft/netcdf_64/v3.6.2_gf_64/include
    #NETCDF_LIB     = -L/soft/netcdf_64/v3.6.2_gf_64/lib -lnetcdf
    NETCDF_INCLUDE = -I/usr/include
    NETCDF_LIB     = -L/usr/lib -lnetcdf -lnetcdff
  else
    NETCDF_INCLUDE = -I/soft/netcdf/v3.6.2_gf/include
    NETCDF_LIB     = -L/soft/netcdf/v3.6.2_gf/lib -lnetcdf
  endif
endif

# ----------------------------------------------------------------------------

ifeq ($(COMPILER),pgf90)
  ### portland group f90 (needs version 5.2-4 or newer!)
  ### (input provided by Tom Sobieraj and Jacek Kaminski)
  F90            = pgf90
  F90FLAGS       = -Mpreprocess -O
  ### next line suggested by Mike Long (email 2008-10-14):
  #F90FLAGS      = -Mpreprocess -O3 -Munroll -Mnoframe -Mlre
  NETCDF_INCLUDE = -I/usr/local/include
  NETCDF_LIB     = -L/usr/local/lib -lnetcdf
endif

# ----------------------------------------------------------------------------

ifeq ($(COMPILER),xlf95)
  ### AIX:
  F90       = xlf95_r
  F90FLAGS  = -q64 -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -O3 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qfree=f90 -qspillsize=32648 -qarch=auto -qtune=auto -Q -qhot
  NETCDF_INCLUDE = -I/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.0b6/include
  NETCDF_LIB     = -L/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.0b6/lib -lnetcdf
endif

##############################################################################

### To implement the boxmodel on a new machine, you can simply overwrite
### the above definitions by activating the following block and entering
### appropriate values for your system:
### F90:            Fortran90 compiler
### F90FLAGS:       Fortran90 compiler options (including the cpp option!)
### NETCDF_INCLUDE: path to netcdf include files
### NETCDF_LIB:     path to netcdf library and library name

#F90            =
#F90FLAGS       =
#NETCDF_INCLUDE =
#NETCDF_LIB     =

##############################################################################

### if ASCII output is selected, it is likely that there is no netcdf library
### available. Therefore, the coresponding variables are set to "nothing":

ifeq ($(OUTPUT),ASCII)
  NETCDF_INCLUDE  =
  NETCDF_LIB      =
endif

##############################################################################

# targets
include main.mk

# list of dependencies (via USE statements)
include depend.mk
