# -*- Makefile -*-
# ----------------------------------------------
PROG = channel.exe

SRCS0 = $(wildcard *.f90)
SRCS  = $(filter-out F%.f90, $(SRCS0))

OBJS := $(SRCS:.f90=.o)

F90DEFS0 = NOMPI MBM_CHANNEL
F90DEFS  = $(addprefix $(DEFOPT), $(F90DEFS0))

MAKEFILE_INC = depend.mk
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

CHANNEL_LIB     = $(CDI_LIB) $(NETCDF_LIB) $(PNETCDF_LIB)
CHANNEL_INCLUDE = $(CDI_INCLUDE) $(NETCDF_INCLUDE) $(PNETCDF_INCLUDE)

ifeq ($(strip $(POSIX90_DEF)), HAVE_POSIX90)
CHANNEL_INCLUDE += $(MODOPT) ../../../libsrc/posix90/src
CHANNEL_LIB     += -L../../../lib -lposix90
endif

ifeq ($(strip $(FORPY_DEF)), HAVE_FORPY)
CHANNEL_INCLUDE += $(MODOPT) ../../../libsrc/forpy
CHANNEL_LIB     += -L../../../lib -lforpy $(PYTHON_LIB)
endif

all: $(PROG)

.PHONY:depend

# update file dependencies
depend $(MAKEFILE_INC):
	$(F_makedepend) $(SRCS)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f *.log
	rm -f *~
	rm -f *.nc
	rm -f *.mc

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

.PHONY: run
run: all
	$(srcdir)/$(PROG)

.PHONY: zip
zip:
	zip -or messy_channel_mbm.zip * -x '*.o' -x '*.mod' -x '*.log' -x '*~' -x '*.old' -x Makefile.m -x '*.zip' -x '*.exe'

.PHONY: list
list:
	@echo '------------------------------------------------------'
	@echo "SRCS = $(SRCS)"
	@echo '------------------------------------------------------'

.PHONY: check
check:
	-forchk -define BOX -rigor -cond -f95 -obs -ff -decl -ext -intr -spec -ancmpl -anprg -anref -shcom -shinc -shmod -shprg -shref -shsrc -shsub -inf -plen 25 -pwid 132 -l fchk_tmp.lst -rep fchk_report.log *.f90 >& forcheck.log
	@echo 'CREATING LOG-FILES ...........................................'
	@rm -f *.log
	@gawk -f ./lst2log.gawk fchk_tmp.lst
	@rm -f fchk_tmp.lst
	@ls -l *.log
	@echo '..............................................................'

$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(F90DEFS) -o $@ $(OBJS) $(CHANNEL_LIB)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(F90DEFS) $(CHANNEL_INCLUDE) -c $<

# ----------------------------------------------
