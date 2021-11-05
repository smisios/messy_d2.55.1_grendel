##############################################################################

### Choose a compiler:
COMPILER = g95
#COMPILER = ifort
#COMPILER = lahey

### 32 or 64-bit architecture?
ifeq "$(findstring 64,$(MACHTYPE))" "64"
  BITS = 64
else
  BITS = 32
endif

### directory for executable:
INSTALLDIR = .

##############################################################################

ifeq ($(COMPILER),lahey)
  ### Lahey compiler:
  F90        = lf95
  ### Choose compiler options:
  #F90FLAGS   = -Cpp
  F90FLAGS   = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap
  F90R8      = --dbl
endif

# ----------------------------------------------------------------------------

ifeq ($(COMPILER),ifort)
  ### ifort compiler:
  F90        = ifort
  ### Choose compiler options:
  F90FLAGS   = -cpp -O0 -g
  F90R8      = -r8
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
  F90R8    = -r8
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

##############################################################################

# targets
include main.mk

# list of dependencies (via USE statements)
include depend.mk

##############################################################################
