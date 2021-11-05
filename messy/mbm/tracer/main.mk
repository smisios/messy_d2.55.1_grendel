# -*- Makefile -*-
# ----------------------------------------------
PROG = tracer.exe

SRCS0 = $(wildcard *.f90)
SRCS  = $(filter-out F%.f90, $(SRCS0))

OBJS := $(SRCS:.f90=.o)

MAKEFILE_INC = depend.mk
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

ifeq ($(strip $(POSIX90_DEF)), HAVE_POSIX90)
TRACER_INCLUDE += $(MODOPT) ../../../libsrc/posix90/src
TRACER_LIB     += -L ../../../lib -lposix90
endif

all: $(PROG)

.PHONY: depend

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

distclean: clean
#	rm -f chemprop/messy_main_tracer_chemprop.inc
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

.PHONY: run
run: all ptrac.nml-raw
	rm -fr ptrac.nml
	./copynml.sh ptrac.nml-raw ptrac.nml
	$(srcdir)/$(PROG)

.PHONY: zip
zip:
	zip -or messy_tracer_mbm.zip * -x '*.o' -x '*.mod' -x '*.log' -x '*~' -x '*.old' -x Makefile.m -x '*.zip' -x '*.exe'

.PHONY: list
list:
	@echo '------------------------------------------------------'
	@echo "SRCS = $(SRCS)"
	@echo '------------------------------------------------------'

.PHONY: check
check:
	-forchk -define MBM_TRACER -rigor -cond -f95 -obs -ff -decl -ext -intr -spec -ancmpl -anprg -anref -shcom -shinc -shmod -shprg -shref -shsrc -shsub -inf -plen 25 -pwid 132 -l fchk_tmp.lst -rep fchk_report.log *.f90 >& forcheck.log
	@echo 'CREATING LOG-FILES ...........................................'
	@rm -f *.log
	@gawk -f ./lst2log.gawk fchk_tmp.lst
	@rm -f fchk_tmp.lst
	@ls -l *.log
	@echo '..............................................................'

$(PROG): depend messy_main_tracer_chemprop.inc $(OBJS)
	$(F90) $(F90FLAGS) $(F90DEFS) $(DEFOPT)MBM_TRACER -o $@ $(OBJS) $(TRACER_LIB)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(F90DEFS) $(F90DEFS) $(DEFOPT)MBM_TRACER $(TRACER_INCLUDE) -c $<

messy_main_tracer_chemprop.inc: chemprop/messy_main_tracer_chemprop.tbl
	(cd chemprop; ./xchemprop)

# ----------------------------------------------
