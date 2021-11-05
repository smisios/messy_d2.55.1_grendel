# -*- Makefile -*- ! mz_sg_20210410(@buran)
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
#HOSTNAME=`hostname`
#HOST=`host $HOSTNAME`
#if (( $? )) ; then
#HOST=$HOSTNAME
#else
#HOST=`host $HOSTNAME | awk '{print $1}'`
#fi

##### C-COMPILER AND FLAGS
case $HOST in

  #########################################################################
  ### Buran Cluster @IGCE Moscow
  #########################################################################
  buran*)
     ### available compiler
     if test -z "$COMPILER" ; then
        COMPILER=`which mpif90 >/dev/null && mpif90 -show || echo "<none>"`
        COMPILER=`echo ${COMPILER} | awk '{print $1}'`
        echo -e "\ndetected compiler: ${COMPILER}"
     fi

   # using tools from custom/OSL software stack (module netcdf-fortran) with the fallback option of $NETCDFROOT
    #NETCDFROOT = ${NETCDF_DIR}
     NETCDF_CFG=`which nf-config 2>/dev/null || which nf-fortran-config 2>/dev/null`
     SPEC_NETCDF_INC = `$NETCDF_CFG --includedir || echo "${NETCDFROOT}/include"`
     SPEC_NETCDF_LIB = `$NETCDF_CFG --flibs      || echo "-L${NETCDFROOT}/lib -lnetcdf -lnetcdff"`
   # MPI_DIR is set via custom/OSL software stack (modules openmpi/impi)
     MPIROOT    = ${MPI_DIR}
    #MPI_LIB    =
    #INCLUDE    =
    #PNETCDF_DIR is set via custom/OSL software stack (modules pnetcdf)
     PNETCDFROOT = ${PNETCDF_DIR}
    #YAXT_DIR is set via custom software stack @buran (modules yaxt)
     YAXTROOT   = ${YAXT_DIR}
     BLASROOT   =
     BLAS_LIB   =
     LAPACKROOT =
     LAPACK_LIB =

   # compiler-specific settings
     case `basename $COMPILER` in
      ### GNU: requires 'ml restore messy-gnu' or 'module load gnu openmpi netcdf-fortran [pnetcdf] [yaxt]'
      gfortran*)
         # currently @ gnu mpif90
         ### ### gcc/gcc481
         ### ### openmpi/1.6.5/gfortran/4.8.1 
         CXX      = g++
         CC       = mpicc
         CFLAGS   = -g -O -Df2cFortran -mavx
         F90      = `which \mpif90`
         F90VERS  = `$F90 --version | awk '{if (NR == 1) print}'`
         ##### F95 EXTENSIONS
         DEFOPT   = -D
         MODOPT   = -I
         ##### F95 COMPILER FLAG TO OBTAIN DOUBLE PRECISION
         F90R8    = -fdefault-real-8
         ##### F77 OPTIONS
         FFLAGS   = -fno-second-underscore -ffree-line-length-none -fno-range-check
         ##### F90 COMMON OPT.
         F90COMOPT = $FFLAGS -cpp -traditional-cpp -D__linux__
         ###
         case $RUNMODE in
            DEBUG)
               ### for debugging:
               F90FLAGS = $F90COMOPT -O0 -fbacktrace -fbounds-check -g -Wall -fdump-core
               # -pedantic
               # -DNOENDIANCONVERT
               # -std=f95
               ;;
            DEBUGOPT)
               ### for debugging:
               F90FLAGS = $F90COMOPT -O3 -fbacktrace -g -fdump-core
               ;;
            PRODUCTION)
               ### fast, optimized model run:
               F90FLAGS = $F90COMOPT -O3
               ;;
            *)
               ERRSTR='Error: wrong run-time mode selected! To solve ./configure [options] RUNMODE=DEBUG|PRODUCTION'
               ;;
         esac
         ;;

      ### INTEL: requires 'ml restore messy-intel' or 'module load intel impi netcdf-fortran [pnetcdf] [yaxt]'
       ifort*)
         CXX      = mpiicpc
         CC       = mpiicc
         CFLAGS   = -g -O -Df2cFortran -xHost
         F90      = `which \mpif90`
         F90VERS  = `$F90 --version | awk '{if (NR == 1) print}'`
         ##### F95 EXTENSIONS
         DEFOPT   = -D
         MODOPT   = -I
         ##### F95 COMPILER FLAG TO OBTAIN DOUBLE PRECISION
         F90R8    = -autodouble
         ##### F77 OPTIONS
         FFLAGS   = -g -fpp -O2 -fp-model strict -xHost
         ##### F90 COMMON OPT.
         F90COMOPT = -sox -cpp -D__INTEL_v11 -fpp -fp-model strict -align all -fno-alias -lpthread -save-temps
         ### -mp1 -mieee-fp
         FORT_INTEGER_LEN=4
         FORT_REAL_LEN=4
         # the following solves KP4 compilation problem with older intel/16 not supporting newer gnu/7 C/C++ headers
         CC       = $CC -D__PURE_INTEL_C99_HEADERS__ -D_Float32=float -D_Float64=double -D_Float32x=_Float64 -D_Float64x=_Float128
         CXX      = $CXX '-D__builtin_addressof(obj)=reinterpret_cast<_Tp*>(&const_cast<char&>(reinterpret_cast<const volatile char&>(__r)))'
         ###
         case $RUNMODE in
            DEBUG)
               F90FLAGS = $F90COMOPT -traceback -O0 -g -debug all -check noarg_temp_created
              #-ftrapuv
               ;;
            DEBUGOPT)
              #F90FLAGS = $F90COMOPT -traceback -O3 -g -debug all -check noarg_temp_created
               F90FLAGS = $F90COMOPT -traceback -O2
               ;;
            PRODUCTION)
               F90FLAGS = $F90COMOPT -O2
               # -Bstatic
               # -no-prec-sqrt
               # -qlargepage -blpdata
               # -fno-alias -no-ansi-alias
               # -static
               ;;
            *)
               ERRSTR='Error: wrong run-time mode selected! To solve ./configure [options] RUNMODE=DEBUG|PRODUCTION'
               ;;
         esac
         ;;

       *)
         ERRSTR='Error: no or wrong compiler suite selected! To solve ./configure [options] COMPILER=GFORTRAN'
         ;;

     esac
     ;;

esac

###########################################################################
###########################################################################
###########################################################################
