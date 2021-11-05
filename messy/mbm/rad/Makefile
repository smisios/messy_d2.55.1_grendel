##############################################################################

INSTALLDIR = .
##############################################################################

### Try to find the compiler for current machine automatically based on
### operating system:

SYSTEM := $(shell uname)
MYBIT  := $(shell uname -m)
HOST   := $(shell hostname)
ifeq ($(findstring 64,$(MYBIT)),64)
  BIT = 64
else
  BIT = 32
endif

ifeq ($(SYSTEM),OSF1)
  compiler = f90_alpha
endif

ifeq ($(SYSTEM),Linux)
  #compiler = lahey
  compiler = intel
  #compiler = gfortran
  #compiler = sun
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
  #F90        = ifort
  F90        	  = /sw/rhel6-x64/mpi/openmpi-2.0.2p1_hpcx-intel14/bin/mpif90
#  F90FLAGS   = -cpp -O0 -autodouble -fp-model precise
#  F90FLAGS   = -cpp -Bdynamic -O0 -g -check bounds -fltconsistency -fp-speculation off -prec-sqrt -align all -fno-alias
  F90FLAGS   = -sox -fpp -g -O2 -xCORE-AVX2 -fp-model strict -align all -save-temps -DBULL -I/sw/rhel6-x64/sys/bullxlib-1.0.0/include -L/sw/rhel6-x64/sys/bullxlib-1.0.0/lib -Wl,-rpath,/sw/rhel6-x64/sys/bullxlib-1.0.0/lib -lbullxMATH  -no-wrap-margin
#  F90FLAGS   = -cpp -O0 -g -fltconsistency -fp-speculation off -prec-sqrt -align all -fno-alias
#  F90FLAGS   = -cpp -O3 -align all -fno-alias
#  NETCDF_INCLUDE = -I/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.3/include
#  NETCDF_LIB     = -L/afs/ipp-garching.mpg.de/home/j/joeckel/extra_space/@sys/netcdf/v3.6.3/lib -lnetcdf
#  NETCDF_INCLUDE  = -I/home/aktuell/local/netcdf/netcdf_ifc_${BIT}/include
#  NETCDF_LIB      = -L/home/aktuell/local/netcdf/netcdf_ifc_${BIT}/lib -lnetcdf
#  NETCDF = $(NETCDF_ROOT)
#  NETCDF_INCLUDE = -I${NETCDF_ROOT}include -I${NETCDF_GF_ROOT}include 
#  NETCDF_LIB     = -L${NETCDF_ROOT}lib -L${NETCDF_GF_ROOT}lib -lnetcdff -lnetcdf \
                   -L${HDF5_ROOT}lib -lhdf5_hl -lhdf5  -L${ZLIB_ROOT}lib -lz -lm \
                   -L/net/opt/szip/2.1/gcc472/7/x86_64/lib -lsz
  
  NETCDF_LIB      = -L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/lib -lnetcdff -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/lib -L/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-openmpi2-intel14/lib -L/sw/rhel6-x64/hdf5/hdf5-1.8.18-parallel-openmpi2-intel14/lib -L/sw/rhel6-x64/sys/libaec-0.3.4-gcc48/lib -L/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.1-openmpi2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-openmpi2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/hdf5/hdf5-1.8.18-parallel-openmpi2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/sys/libaec-0.3.4-gcc48/lib -Wl,-rpath,/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.1-openmpi2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/hdf4/hdf-4.2.10-gcc48/lib -lnetcdf -lhdf5 -lhdf5_hl -lsz -lz -lcurl -lnetcdf
  NETCDF_INCLUDE  = -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/include -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/include
  ifeq ($(BIT),64)
    MPI_DIR         = ${MPICH_ROOT}
  else
    MPI_DIR         =  /home/aktuell/local/mpich2/mpich2_ifc_${BIT}
  endif
  MPI_INCLUDE     = -I${MPI_DIR}/include
  MPI_LIB         = -L${MPI_DIR}/lib -lmpichf90 -lmpich -lmpl
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),sun)
  ### SUN f90 compiler:
  F90        = sunf95
  F90FLAGS   = -fpp -xpp=fpp -O0 -g
  ifeq ($(HOST),strat06)
    NETCDF_INCLUDE  = -I/usr/local/netcdf/netcdf_sun_${BIT}/include
    NETCDF_LIB      = -L/usr/local/netcdf/netcdf_sun_${BIT}/lib -lnetcdf
    MPI_DIR         = /home/kunze/local/mpich2/mpich2_sun_$(BIT)
  else
    NETCDF_INCLUDE  = -I/home/aktuell/local/netcdf/netcdf_sun_${BIT}/include
    NETCDF_LIB      = -L/home/aktuell/local/netcdf/netcdf_sun_${BIT}/lib -lnetcdf
    MPI_DIR         = /home/aktuell/local/mpich2/mpich2_sun_$(BIT)
  endif
  MPI_INCLUDE     = -I${MPI_DIR}/include
  MPI_LIB         = -L${MPI_DIR}/lib64 -lmpichf90 -lmpich -lmpl
endif

# ----------------------------------------------------------------------------

ifeq ($(compiler),gfortran)
  ### SUN f90 compiler:
  F90        = gfortran
  F90FLAGS   = -O0 -g -static -static-libgfortran -x f95-cpp-input -fno-range-check
  NETCDF_INCLUDE  = -I${HOME}/local/netcdf/netcdf_gf95_${BIT}/include
  NETCDF_LIB      = -L${HOME}/local/netcdf/netcdf_gf95_${BIT}/lib -lnetcdf
  MPI_DIR         := $(shell dirname ${MPICH2_LIB})
  MPI_DIR         =  ${HOME}/local/mpich2/mpich2_gf95_${BIT}
#  MPI_DIR         = ${MPI_HOME}
  MPI_INCLUDE     = -I${MPI_DIR}/include
  MPI_LIB         = -L${MPI_DIR}/lib -lmpichf90 -lmpich -lmpl
#  MPI_LIB         = -L${MPI_DIR}/lib -lmpi_f90 -lvt-mpi

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
