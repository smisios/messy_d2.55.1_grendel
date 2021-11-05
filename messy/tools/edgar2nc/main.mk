### name of the executable that will be produced
PROG       = $(INSTALLDIR)/edgar2nc.exe

# complete list of all f90 source files
SRCS  = mo_f2kcli.f90 \
        edgar2nc.f90

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
	rm -f *~
	rm -f $(OBJS) *.mod

distclean: clean
	rm -f $(PROG)
	rm -f *.nc

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) -c $<

# ------------------------------------------------------------------
edgar2nc.o: edgar2nc.f90 mo_f2kcli.o
mo_f2kcli.o: mo_f2kcli.f90
# ------------------------------------------------------------------
