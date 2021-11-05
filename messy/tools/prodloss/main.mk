### name of the executable that will be produced
PROG       = $(INSTALLDIR)/prodloss.exe

# complete list of all f90 source files
SRCS  = prodloss.f90

# --------------------------------------------------------------------
OBJS := $(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) -o $@

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
	$(F90) $(F90FLAGS) -c $<

# ------------------------------------------------------------------
prodloss.o: 
# ------------------------------------------------------------------
