### name of the executable that will be produced
PROG       = $(INSTALLDIR)/mpiom_grid.exe

# complete list of all f90 source files
SRCS  = mpiom_grid.f90 

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

zip:
	zip -r mpiom_grid.zip *
clean:
	rm -f *~
	rm -f $(OBJS) *.mod


distclean: clean
	rm -f $(PROG)
	rm -f ./output/*
	rm -f *.nc

%.o: %.f90
	$(F90) $(F90FLAGS) $(LIBS) $(INCLUDES) -c $<

# ------------------------------------------------------------------
# makdepf90 *.f90
mpiom_grid.o : mpiom_grid.f90  
# ------------------------------------------------------------------
