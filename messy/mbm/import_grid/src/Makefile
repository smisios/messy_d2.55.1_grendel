# ----------------------------------------------
SRCS =  mo_f2kcli.f90 \
        messy_main_constants_mem.f90 \
	messy_main_blather.f90 \
	messy_main_import_grid.f90 \
	messy_main_import_grid_par.f90 \
	messy_main_grid.f90 \
        messy_main_grid_mpi.f90 \
	messy_main_grid_netcdf.f90 \
	messy_main_grid_tools.f90 \
	messy_main_grid_trafo.f90 \
	messy_main_grid_trafo_nrgd.f90 \
	messy_main_grid_trafo_nrgd_base.f90 \
	messy_main_grid_trafo_scrp.f90 \
	messy_main_grid_trafo_scrp_base.f90 \
	messy_main_tools.f90 \
	import_grid.f90

OBJS := $(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJS) $(LIBS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -c $<

# ----------------------------------------------
import_grid.o : import_grid.f90 messy_main_constants_mem.o messy_main_tools.o messy_main_grid_tools.o messy_main_grid_trafo_scrp.o messy_main_grid_trafo_nrgd_base.o messy_main_grid_trafo_nrgd.o messy_main_grid_trafo.o messy_main_grid_netcdf.o messy_main_grid.o messy_main_import_grid.o mo_f2kcli.o
messy_main_blather.o : messy_main_blather.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_grid.o : messy_main_grid.f90 messy_main_tools.o messy_main_grid_netcdf.o messy_main_constants_mem.o
messy_main_grid_mpi.o : messy_main_grid_mpi.f90
messy_main_grid_netcdf.o : messy_main_grid_netcdf.f90 messy_main_tools.o messy_main_grid_mpi.o messy_main_constants_mem.o
messy_main_grid_tools.o : messy_main_grid_tools.f90 messy_main_grid.o messy_main_grid_netcdf.o messy_main_constants_mem.o
messy_main_grid_trafo.o : messy_main_grid_trafo.f90 messy_main_grid_netcdf.o messy_main_grid.o messy_main_constants_mem.o messy_main_tools.o
messy_main_grid_trafo_nrgd_base.o : messy_main_grid_trafo_nrgd_base.f90 messy_main_grid_trafo.o messy_main_grid_netcdf.o messy_main_constants_mem.o messy_main_tools.o
messy_main_grid_trafo_nrgd.o : messy_main_grid_trafo_nrgd.f90 messy_main_tools.o messy_main_grid.o messy_main_grid_netcdf.o messy_main_grid_trafo.o messy_main_grid_trafo_nrgd_base.o messy_main_constants_mem.o
messy_main_grid_trafo_scrp_base.o : messy_main_grid_trafo_scrp_base.f90 messy_main_constants_mem.o
messy_main_grid_trafo_scrp.o : messy_main_grid_trafo_scrp.f90 messy_main_grid_tools.o messy_main_tools.o messy_main_grid.o messy_main_grid_trafo.o messy_main_grid_trafo_scrp_base.o messy_main_constants_mem.o messy_main_grid_netcdf.o
messy_main_import_grid.o : messy_main_import_grid.f90 messy_main_grid_trafo_scrp.o messy_main_grid_trafo_nrgd_base.o messy_main_grid_trafo_nrgd.o messy_main_grid_tools.o messy_main_grid.o messy_main_grid_netcdf.o messy_main_grid_trafo.o messy_main_tools.o messy_main_import_grid_par.o messy_main_constants_mem.o
messy_main_import_grid_par.o : messy_main_import_grid_par.f90 messy_main_constants_mem.o messy_main_grid_netcdf.o
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
mo_f2kcli.o : mo_f2kcli.f90
