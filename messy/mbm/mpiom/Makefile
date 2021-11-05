# -*- Makefile -*-
##############################################################################

SHELL = /bin/sh

#SYSTEM := $(shell uname)
SYSTEM := LINUX64
BIT = 64
MPIOMVERS = MPIOM_13B

ifeq ($(SYSTEM),OSF1)
  ### f90 on alphas
  F90      = f90
  F90FLAGS = -O -cpp -D_alpha
  DEFOPT   = -D
  MODOPT = -I
  EXTMODCASE = asis
  #
  NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.0b6/include
  NETCDF_LIB      = -L/soft/netcdf/v3.6.0b6/lib -lnetcdf
  #
else
  ifeq ($(SYSTEM),LINUX64)
  
    ### choose Fortran90 compiler and options

    #compiler = lahey
    compiler = intel
    #compiler = gfortran
    #compiler = sun
    #compiler = g95
    
    MPI = MPI
    ENDIAN = LITTLE_ENDIAN

    ifeq ($(compiler),lahey)
      ### lf95 = Lahey
      F90 = /opt/mpi/bullxmpi/1.2.8.4/bin/mpif90
      ### -Cpp  = run C-preprocessor before compiling
      ###         (necessary because of some compiler directives in the code)
      ### --pca = protect constant argument (necessary for unknown reason)
      #F90FLAGS = -Cpp --pca
      #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0
      #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0 -g --trap --verbose
      F90      = /opt/mpi/bullxmpi/1.2.8.4/bin/mpif90
      F90R8    = -autodouble
      F90VERS  = ifort (IFORT) 14.0.3 20140422
      F90FLAGS = -sox -g -traceback -debug all -check all -fpp -O0 -fp-model strict -align all -save-temps  -no-wrap-margin
      DEFOPT   = -D
      MODOPT   = -I
      EXTMODCASE = asis
      #
      #NETCDF_INCLUDE  = -I/soft/netcdf/v3.6.2_lf/include
      #NETCDF_LIB      = -L/soft/netcdf/v3.6.2_lf/lib -lnetcdf
      #NETCDF_INCLUDE  = -I/athome/joec_pa/_data/soft/netcdf/v3.6.3_lf62d/include
      #NETCDF_LIB      = -L/athome/joec_pa/_data/soft/netcdf/v3.6.3_lf62d/lib -lnetcdf
      LIBSRCS         =  messy/smcl mpiom/src libsrc/cmatutils libsrc/isorropia libsrc/qhull echam5/support echam5/lapack echam5/blas
      NETCDF_INCLUDE  = -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/include -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/include
      NETCDF_LIB      = -L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/lib -lnetcdff -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/lib -L/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-gcc48/lib -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-gcc48/lib -L/sw/rhel6-x64/hdf5/hdf5-1.8.14-threadsafe-gcc48/lib -Wl,-rpath,/sw/rhel6-x64/hdf5/hdf5-1.8.14-threadsafe-gcc48/lib -L/sw/rhel6-x64/sys/libaec-0.3.2-gcc48/lib -Wl,-rpath,/sw/rhel6-x64/sys/libaec-0.3.2-gcc48/lib -lnetcdf -lhdf5_hl -lhdf5 -lsz -lcurl -lz -lnetcdf
      # FORCHECK
      #FCKLIBS  = /data/joec_pa/soft/forcheck/v14/lib/netcdf90
      FCKLIBS  = /afs/ipp/mpc/soft/@sys/forcheck/v14/netcdf90
    endif

#     ### g95 (64-bit)
#     F90 = g95_64
#     F90FLAGS = -fno-second-underscore -cpp -O0 -fbounds-check -ftrace=full
#     #F90FLAGS = -fno-second-underscore -cpp
#     DEFOPT = -D
#     MODOPT = -I
#     EXTMODCASE = asis
#     #
#     NETCDF_INCLUDE  = -I/soft/netcdf_64/v4.0_g95_64/include
#     NETCDF_LIB      = -L/soft/netcdf_64/v4.0_g95_64/lib -lnetcdf

    ifeq ($(compiler),intel)
      ### intel
      F90 = /sw/rhel6-x64/mpi/openmpi-2.0.2p1_hpcx-intel14/bin/mpif90
      F90FLAGS = -sox -fpp -g -O2 -xCORE-AVX2 -fp-model strict -align all -save-temps -DBULL -I/sw/rhel6-x64/sys/bullxlib-1.0.0/include -L/sw/rhel6-x64/sys/bullxlib-1.0.0/lib -Wl,-rpath,/sw/rhel6-x64/sys/bullxlib-1.0.0/lib -lbullxMATH  -no-wrap-margin
      DEFOPT = -D
      MODOPT = -I
      EXTMODCASE = asis
      
      #F90FLAGS   = -cpp -Bdynamic -O0 -g -check bounds -fltconsistency -fp-speculation off -prec-sqrt -align all -fno-alias
      #F90FLAGS   = -cpp -Bdynamic -O3 -fltconsistency -fp-speculation off -prec-sqrt -align all -fno-alias

      #NETCDF_ROOT    = /pf/b/b324024/local/netcdf/netcdf_ifc_64
      #NETCDF         = $(NETCDF_ROOT)
      #NETCDF_INCLUDE = -I${NETCDF_ROOT}include -I${NETCDF_GF_ROOT}include 
      #NETCDF_LIB     = -L${NETCDF_ROOT}lib -lnetcdf
      
      NETCDF_LIB      = -L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/lib -lnetcdff -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/lib -L/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-openmpi2-intel14/lib -L/sw/rhel6-x64/hdf5/hdf5-1.8.18-parallel-openmpi2-intel14/lib -L/sw/rhel6-x64/sys/libaec-0.3.4-gcc48/lib -L/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.1-openmpi2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-openmpi2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/hdf5/hdf5-1.8.18-parallel-openmpi2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/sys/libaec-0.3.4-gcc48/lib -Wl,-rpath,/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.1-openmpi2-intel14/lib -lnetcdf -lhdf5 -lhdf5_hl -lsz -lz -lcurl -lnetcdf
      NETCDF_INCLUDE  = -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/include -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/include
                       
      ifeq ($(BIT),64)
         MPI_DIR         = ${MPICH_ROOT}
      else
         MPI_DIR         =  /home/aktuell/local/mpich2/mpich2_ifc_${BIT}
      endif
      MPI_INCLUDE     = -I${MPI_DIR}/include
      MPI_LIB         = -L${MPI_DIR}/lib -lmpichf90 -lmpich -lmpl
    endif
    ifeq ($(compiler),gfortran)
      ### GNU gfortran
      MPICH_ROOT = ${HOME}/local/mpich2/mpich2_gf95_64
      
      F90 = ${MPICH_ROOT}/bin/mpif90
      F90FLAGS   = -O0 -g -rdynamic  -x f95-cpp-input -fno-range-check -ffixed-line-length-none
      F90FLAGSLD = -O0 -g -rdynamic  -fno-range-check -ffixed-line-length-none
      F90R8    = -autodouble
      DEFOPT = -D
      MODOPT = -I
      EXTMODCASE = asis
      
      LIBSRCS         =  messy/smcl mpiom/src libsrc/cmatutils libsrc/isorropia libsrc/qhull echam5/support echam5/lapack echam5/blas
      #F90FLAGS   = -cpp -Bdynamic -O0 -g -check bounds -fltconsistency -fp-speculation off -prec-sqrt -align all -fno-alias
      #F90FLAGS   = -cpp -Bdynamic -O3 -fltconsistency -fp-speculation off -prec-sqrt -align all -fno-alias

      NETCDF_ROOT=/usr/local/netcdf/netcdf-4.4.1.1/
      NETCDF_GF_ROOT=/usr/local/netcdf/netcdf-fortran-4.4.4/
      HDF5_ROOT=/usr/
      ZLIB_ROOT=/usr/
      NETCDF_INCLUDE = -I${NETCDF_ROOT}include -I${NETCDF_GF_ROOT}include 
      NETCDF_LIB     = -L${NETCDF_ROOT}lib64 -L${NETCDF_GF_ROOT}lib64 -lnetcdff -lnetcdf \
		       -L${HDF5_ROOT}lib64 -lhdf5_hl -lhdf5  -L${ZLIB_ROOT}lib64 -lz -lm
      ifeq ($(BIT),64)
         MPI_DIR         = ${MPICH_ROOT}
      else
         MPI_DIR         =  /home/aktuell/local/mpich2/mpich2_ifc_${BIT}
      endif
      MPI_INCLUDE     = -I${MPI_DIR}/include
      MPI_LIB         = -L${MPI_DIR}/lib64 -lmpichf90 -lmpich -lmpl
    endif

  else
	ERROR
  endif
endif

# activate following lines for intel compiler (8.0.039) to overwrite
# previous settings
#F90      = ifort
#F90FLAGS = -cpp -O0
#DEFOPT   = -D
#MODOPT = -I
#EXTMODCASE = asis

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
