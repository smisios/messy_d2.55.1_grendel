# ----------------------------------------------
SRCS =  mo_f2kcli.f90 \
	messy_ncregrid_base.f90 \
	messy_ncregrid_netcdf.f90 \
	messy_ncregrid_geohyb.f90 \
	messy_ncregrid_control.f90 \
	messy_ncregrid_interface.f90 \
	messy_ncregrid_diag.f90 \
	ncregrid.f90

OBJS := $(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJS) $(LIBS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -c $<

# ----------------------------------------------
mo_f2kcli.o : mo_f2kcli.f90
messy_ncregrid_base.o : messy_ncregrid_base.f90 messy_ncregrid_mpi.o
messy_ncregrid_control.o : messy_ncregrid_control.f90 \
	messy_ncregrid_diag.o \
	messy_ncregrid_interface.o \
	messy_ncregrid_geohyb.o \
	messy_ncregrid_netcdf.o \
	messy_ncregrid_base.o
messy_ncregrid_diag.o : messy_ncregrid_diag.f90 \
	messy_ncregrid_geohyb.o \
	messy_ncregrid_netcdf.o \
	messy_ncregrid_base.o
messy_ncregrid_geohyb.o : messy_ncregrid_geohyb.f90 \
	messy_ncregrid_netcdf.o \
	messy_ncregrid_base.o
messy_ncregrid_interface.o : messy_ncregrid_interface.f90 \
	messy_ncregrid_geohyb.o \
	messy_ncregrid_base.o
messy_ncregrid_netcdf.o : messy_ncregrid_netcdf.f90 \
	messy_ncregrid_base.o
messy_ncregrid_mpi.o : messy_ncregrid_mpi.f90
ncregrid.o : ncregrid.f90 \
	messy_ncregrid_control.o \
	mo_f2kcli.o
# ----------------------------------------------
include ../specific.mk
