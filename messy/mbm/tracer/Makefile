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
endif
ifeq ($(SYSTEM),Linux)
  ### Choose compiler:
  #COMPILER = lahey
  COMPILER = g95
  ifeq ($(COMPILER),lahey)
    ### lf95 = Lahey
    F90 = lf95
    ### -Cpp  = run C-preprocessor before compiling
    ###         (necessary because of some compiler directives in the code)
    ### --pca = protect constant argument (necessary for unknown reason)
    #F90FLAGS = -Cpp --pca
    #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0
    #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0 -g --trap --verbose
    F90FLAGS = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap
    DEFOPT = -D
  endif
  ifeq ($(COMPILER),g95)
    ### g95 version 0.93:
    # F90 = g95
    F90  = /soft/g95_64/v0.93/bin/g95
    ### Choose compiler options:
    ### -cpp                    = run the C preprocessor
    ### -fimplicit-none         = IMPLICIT NONE, unless overridden by IMPLICIT statements
    ### -fno-second-underscore  = don't append a second trailing underscore in
    ###                           names having an underscore
    ### -ftrace=full            = show line number when arithmetic exceptions occur
    ### -fno-backslash          = treat backslashes literally
    ### -fbounds-check          = check array and substring bounds at runtime
    ### -g                      = debug option
    ### -O<n>                   = optimization level <n>=0...3, zero = no optimization
    ### -O                      = same as -O1
    ### -Wall                   = Enable most warning messages
    F90FLAGS   = -cpp -fimplicit-none -fno-second-underscore -ftrace=full -fno-backslash -fbounds-check -g -O0 -Wall
    DEFOPT = -D
  endif
else
  ERROR
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

include main.mk

include depend.mk
