# -*- Makefile -*-
###############################################################################
##### PLATFORM SPECIFIC SETTINGSFOR mistral @ DKRZ ############################
###############################################################################
ARCH       = LINUX64

###########################################################################
### DEFAULTS
###########################################################################
#### DEFAULT COMPILE MODE
if test -z "$RUNMODE" ; then
   RUNMODE=PRODUCTION
fi

#### DEFAULT COMPILER
if test -z "$COMPILER" ; then
   COMPILER=INTEL
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

  #########################################################################
  ### cartesius.surfsara.nl
  #########################################################################

     # TMPF90=`mpif90 -show | awk '{print $1}'`
     # case $TMPF90 in
     #      ifort*)
     #        COMPILER=INTEL
     #        ;;
     #      gfortran*)
     #        COMPILER=GFORTRAN
     #        ;;
     #      *)
     #        ERRSTR='Error: no openmpi module loaded!'
     #        ;;
     # esac

     case $COMPILER in

      INTEL)
         ### ### intel/2016 (16.0.3) mpi/openmpi-1.8.8-intel
         ### special flag to run optimized
         ### - on HASWELL processors only (requires
         ###   #SBATCH --constraint=haswell
         ###   in run-script
	 #SPECFLAG = -xCORE-AVX2
         ### - on IVY/SANDY BRIDGE and HASWELL processors
         #SPECFLAG = -xAVX -axCORE-AVX2
         ###
         SPECFLAG =
         ###
         CXX      = g++
         CC       = mpicc
         CFLAGS   = -O -Df2cFortran $SPECFLAG `nf-config --cflags`
         #F90 = mpif90
         F90 = `which \mpiifort`
         F90VERS  = `$F90 --version | awk '{if (NR == 1) print}'`
         ##### F95 EXTENSIONS
         DEFOPT   = -D
         MODOPT   = -I
         ##### F95 COMPILER FLAG TO OBTAIN DOUBLE PRECISION
         F90R8    = -autodouble
         ##### F77 OPTIONS
         FFLAGS   = -fpp -O2 -fp-model strict $SPECFLAG
         ##### OpenMP options
         F90OMP    = -qopenmp -heap-arrays
         F90OMPLIB = -qopenmp -heap-arrays
         ##### MPI C++ compiler, flags and libs (for GUESS)
         MPICXX  = $zMPICXX
         ###
         MPIROOT    =
         MPI_LIB    =
         #PNETCDFROOT = /export/opt/PA/prgs/pnetcdf/1.6.0/ifort/14.0-4
         NETCDFROOT =
         #SPEC_NETCDF_INC = /home/fyin/netcdf/intel
#         SPEC_NETCDF_INC = /home/fyin/netcdf
#         SPEC_NETCDF_LIB = -L/usr/lib64 -lnetcdf -lnetcdff
         SPEC_NETCDF_INC = `nf-config --fflags | sed 's|-I||g'`
         SPEC_NETCDF_LIB = `nf-config --flibs`
         BLASROOT   =
         BLAS_LIB   =
         LAPACKROOT =
         LAPACK_LIB =
         SPEC_LABLA_LIB =
         SPEC_LIB   =

         case $RUNMODE in
            DEBUG)
               F90FLAGS = -sox -g $SPECFLAG -traceback -debug all -check all,noarg_temp_created -fpp -O0 -fp-model strict -align all -save-temps -no-wrap-margin
# -align all
               MPICXXFLAGS = -std=c++98 -g $SPECFLAG -traceback -O0 -DDEBUG -fp-model strict
               ;;
            DEBUGOPT)
               F90FLAGS = -sox -g -traceback -debug all -check all -fpp -O2 $SPECFLAG -fp-model strict -align all -save-temps -no-wrap-margin
               MPICXXFLAGS = -std=c++98 -O2 $SPECFLAG -DDEBUG -fp-model strict -fp-model source -fp-speculation=off -fno-alias -no-fma -fimf-arch-consistency=true -align -fno-alias -lpthread
               ;;
            PRODUCTION)
               F90FLAGS = -sox -fpp -g -O2 $SPECFLAG -fp-model strict -align all -save-temps -no-wrap-margin
               MPICXXFLAGS = -std=c++98 -O2 $SPECFLAG -DNODEBUG -fp-model strict -fp-model source -fp-speculation=off -fno-alias -no-fma -fimf-arch-consistency=true -align -fno-alias -lpthread
               ;;
            *)
               ERRSTR='Error: wrong run-time mode selected! To solve ./configure [options] RUNMODE=DEBUG|PRODUCTION'
               ;;
         esac
         ;;

      GFORTRAN)
         ### ### gcc/gcc483
         ### ### mpi/openmpi-1.8.8-gnu
         CXX      = g++
         CC       = mpicc
         CFLAGS   = -g -O -Df2cFortran `nf-config --cflags`
         F90      = `which \mpif90`
         F90VERS  = `$F90 --version | awk '{if (NR == 1) print}'`
         ##### F95 EXTENSIONS
         DEFOPT   = -D
         MODOPT   = -I
         ##### F95 COMPILER FLAG TO OBTAIN DOUBLE PRECISION
         F90R8    = -fdefault-real-8
         ##### F77 OPTIONS
         FFLAGS   = -fno-second-underscore -ffree-line-length-none -fno-range-check
         ##### OpenMP options
         F90OMP    =
         F90OMPLIB =
         ##### MPI C++ compiler, flags and libs (for GUESS)
         MPICXX = $zMPICXX -std=c++1y
         #####
         MPIROOT    =
         MPI_LIB    =
         #PNETCDFROOT = /export/opt/PA/prgs/pnetcdf/1.6.0/gfortran/4.8.1
         NETCDFROOT =
         #SPEC_NETCDF_INC = /home/fyin/netcdf/gnu
         #SPEC_NETCDF_INC = /home/fyin/netcdf
         SPEC_NETCDF_INC = `nf-config --fflags | sed 's|-I||g'`
         SPEC_NETCDF_LIB = `nf-config --flibs`
         BLASROOT   =
         BLAS_LIB   =
         LAPACKROOT =
         LAPACK_LIB =
         SPEC_LABLA_LIB =
         SPEC_LIB   = -I/usr/include

         case $RUNMODE in
            DEBUG)
               ### for debugging:
               F90FLAGS = -cpp -D__linux__ -fno-second-underscore -ffree-line-length-none -fno-range-check -fbacktrace -fbounds-check -g -O0 -Wall -fdump-core
               # -pedantic
               # -DNOENDIANCONVERT
               # -std=f95
               MPICXXFLAGS = -O0 -DNDEBUG
               ;;
            DEBUGOPT)
               ### for debugging:
               F90FLAGS = -cpp -D__linux__ -fno-second-underscore -ffree-line-length-none -fno-range-check -O2 -fbacktrace -fbounds-check -g -Wall -fdump-core
               MPICXXFLAGS = -O2 -DNDEBUG -Wuninitialized -fno-common
               ;;
            PRODUCTION)
               ### fast, optimized model run:
               F90FLAGS = -cpp -D__linux__ -fno-second-underscore -ffree-line-length-none -fno-range-check -O3 -g
               MPICXXFLAGS = -O2 -DNDEBUG -Wuninitialized -fno-common
               ;;
            *)
               ERRSTR='Error: wrong run-time mode selected! To solve ./configure [options] RUNMODE=DEBUG|PRODUCTION'
               ;;
         esac
         ;;

      *)
         ERRSTR='Error: no or wrong compiler suite selected! To solve ./configure [options] COMPILER=INTEL|GFORTRAN'
         ;;

     esac

     FCKLIBS    =

##########################################################################
##########################################################################
##########################################################################
