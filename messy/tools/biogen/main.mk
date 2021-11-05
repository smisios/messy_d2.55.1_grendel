### name of the executable that will be produced
PROG       = $(INSTALLDIR)/biogen.exe

# complete list of all f90 source files
SRCS  = biogen.f90 \
        tools.f90

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
	zip -r biogen.zip *
clean:
	rm -f *~
	rm -f $(OBJS) *.mod

distclean: clean
	rm -f $(PROG)
	rm -f *.nc
	rm -f output/*.nc

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) -c $<

# ------------------------------------------------------------------
# makdepf90 *.f90
biogen.o : biogen.f90 tools.o 
tools.o : tools.f90
# ------------------------------------------------------------------
