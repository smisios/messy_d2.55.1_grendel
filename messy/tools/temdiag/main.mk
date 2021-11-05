### name of the executable that will be produced
PROG       = $(INSTALLDIR)/temdiag.exe

# complete list of all f90 source files
SRCS  = $(wildcard *.f90)

# --------------------------------------------------------------------
OBJS := $(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJS) $(LIBS)

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

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) -c $<

# ------------------------------------------------------------------
messy_main_blather.o : messy_main_blather.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
temdiag.o : temdiag.f90 messy_main_tools.o messy_main_constants_mem.o
# ------------------------------------------------------------------
