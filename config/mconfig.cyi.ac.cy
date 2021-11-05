# -*- Makefile -*- Time-stamp: <2019-02-11 14:46:45 joec_pa>
###############################################################################
##### PLATFORM SPECIFIC SETTINGS (LINUX64) ####################################
###############################################################################
ARCH       = LINUX64

###########################################################################
### DEFAULTS
###########################################################################
#### DEFAULT COMPILE MODE
if test -z "$RUNMODE" ; then
   RUNMODE=PRODUCTION
fi

#### SET DEFAULT PROFILING MODE
if test -z "$PROFILEMODE" ; then
   PROFILEMODE=NONE
fi

##### TOOLS
#AR       = ar
#ARFLAGS  = cvr
#NMFLAGS  = 

##### F95 EXTENSIONS
DEFOPT   = -D
MODOPT   = -I

###########################################################################
##### SYSTEM / HOST SPECIFIC NON-DEFAULTS
###########################################################################

  #########################################################################
  ### cyclone.hpcf.cyi.ac.cy - Cyprus Institute HPC Facility 
  #########################################################################

     CXX      = 
     CC       = mpiicc 
     CFLAGS   = -O -fp-model strict -Df2cFortran
     F90      = mpiifort 
     F90VERS  = `$F90 --version | awk '{if (NR == 1) print}'`
     ##### F95 EXTENSIONS
     DEFOPT   = -D
     MODOPT   = -I
     ##### F95 COMPILER FLAG TO OBTAIN DOUBLE PRECISION
     F90R8    = -autodouble 
     ##### F95 COMPILER FLAGS
     FFLAGS = -fpp -O2 -fp-model strict
     ###
     F90ADDOPT = -xAVX 
     ###
     case $RUNMODE in
        DEBUGOPT*)
           F90FLAGS = -fpp -g -check bounds -O2 -traceback -fp-model strict -fno-alias -no-ansi-alias -lpthread -save-temps
           ;;
        DEBUG*)
           F90FLAGS = -fpp -O0 $F90ADDOPT -debug all -check all -fpe0 -traceback -fp-model strict -fno-alias -no-ansi-alias -lpthread -save-temps
           ;;
        PRODUCTION*)
           F90FLAGS = -fpp -O3 $F90ADDOPT -fp-model strict -fno-alias -no-ansi-alias -lpthread -save-temps
           ;;
     esac

     ### (1) MESSAGE PASSING INTERFACE (options a and b are exclusive!)
     ####    a) use mpi-compiler wrappers (preferable!) and keep
     ####       MPIROOT and MPI_LIB unset.
     ####    b) set MPIROOT and MPI_LIB (e.g. for MPI-1 
     ###        MPI_LIB = mpichf90nc mpich)
     MPIROOT    = $EBROOTIMPI
     MPI_LIB    = 
     ### (2) NETCDF LIBRARY (options a and b are exclusive!)
     ###     a) set NETCDFROOT (must contain lib/libnetcdf.a and
     ###        include/netcdf.inc) (for necdf-3)
     ###     b) set SPEC_NETCDF_INC to include path and
     ###        SPEC_NETCDF_LIB to ld-options (-L ...  -l...)
     NETCDFROOT = $EBROOTNETCDFMINFORTRAN
     SPEC_NETCDF_INC =
     SPEC_NETCDF_LIB =
     ###
     ### (3) BLAS and LAPACK LIBRARIES (options a, b and are exclusive)
     ###     a) keep all entries empty -> blas/lapack contained in 
     ###        distribution will be compiled and used
     ###     b) specify *ROOT (path to library) and
     ###        and *_LIB (name of library without "-llib" and ".a",
     ###        e.g., BLASROOT=/usr/lib and BLAS_LIB=blas for
     ###        /usr/lib/libblas.a)
     ###     c) specifiy SPEC_LABLA_LIB with full ld options)
     BLASROOT   =
     BLAS_LIB   =
     LAPACKROOT =
     LAPACK_LIB =
     SPEC_LABLA_LIB =
     ### (4) EMOS and SPHERE LIBRARIES (for INTERA only); 
     ###     a) keep empty, if not available
     ###     b) similar as option b) for BLAS/LAPACK
     EMOSROOT   =
     EMOS_LIB   =
     SPHEREROOT =
     SPHERE_LIB = 
     ### (5) ADDITONAL LIBRARIES REQUIRED FOR LINKING (full ld options)
     SPEC_LIB   =
     ### (6) SPECIAL LIBRARIES FOR FORCHECK
     FCKLIBS    =     

###########################################################################
###########################################################################
###########################################################################
