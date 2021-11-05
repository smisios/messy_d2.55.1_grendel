### name of the executable that will be produced
PROG       = $(INSTALLDIR)/strato3bud.exe

# complete list of all f90 source files
SRCS  = readData.f90 \
        strato3bud.f90 \
        writeData.f90

NTAG := $(shell grep -E -i 'imdouble|imtag' ../../mbm/caaba/mecca/mecca.eqn | wc -l )

# --------------------------------------------------------------------
OBJS := $(SRCS:.f90=.o)

all: $(PROG)

strato3bud.o: rate_coeff.inc

rate_coeff.inc: ../../mbm/caaba/mecca/mecca.eqn
	@rm -f rate_coeff.inc
	./eqn2f90.tcsh > rate_coeff.inc
	@echo "------------------------------------------------"
	@echo "rate coefficients from mecca.eqn:"
	@echo "------------------------------------------------"
	@cat rate_coeff.inc
	@echo "------------------------------------------------"

$(PROG): rate_coeff.inc $(OBJS)
ifeq ($(NTAG), 0)
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJS) $(LIBS)
else
	@echo "------------------------------------------------"
	@cat rate_coeff.inc
	@echo "------------------------------------------------"
	@echo '###########################################################'
	@echo 'stratobud cannot be compiled for tagged/doubled mechanisms.'
	@echo '###########################################################'
	@echo
endif

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
	rm -f rate_coeff.inc

%.o: %.f90
ifeq ($(NTAG), 0)
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -c $<
else
	@echo
endif

# ------------------------------------------------------------------

# ------------------------------------------------------------------
