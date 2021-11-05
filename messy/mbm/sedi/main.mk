### name of the executable that will be produced
PROG       = $(INSTALLDIR)/sedi.exe

# complete list of all f90 source files
SRCS  = sedi_column.f90 \
        messy_main_constants_mem.f90 \
        messy_main_tools.f90 \
        messy_main_blather.f90 \
        messy_sedi_column.f90 \
        messy_sedi.f90 \
        sedi_column_netcdf.f90

# --------------------------------------------------------------------
OBJS := $(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) $(LIBS) -o $@

list:
	@echo "------------------------------------------------"
	@echo "SRCS = $(SRCS)"
	@echo "------------------------------------------------"
	@echo
	@echo "------------------------------------------------"
	@echo "OBJS = $(OBJS)"
	@echo "------------------------------------------------"

clean:
	rm -f $(OBJS) *.mod *.log

distclean: clean
	rm -f $(PROG)
	rm -f *.nc

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) -c $<

# ------------------------------------------------------------------
messy_main_blather.o : messy_main_blather.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
messy_sedi_column.o : messy_sedi_column.f90 messy_main_tools.o messy_sedi.o sedi_column_netcdf.o messy_main_constants_mem.o
messy_sedi.o : messy_sedi.f90 messy_main_tools.o messy_main_constants_mem.o
sedi_column.o : sedi_column.f90 messy_sedi_column.o messy_main_constants_mem.o
sedi_column_netcdf.o : sedi_column_netcdf.f90 messy_main_constants_mem.o
# ------------------------------------------------------------------
