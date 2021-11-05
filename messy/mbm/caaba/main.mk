# Makefile syntax info:
# When a line starts with '@', the echoing of that line is suppressed.
# When a line starts with '-', errors in the command line are ignored.

### name of the executable that will be produced
PROG = $(INSTALLDIR)/caaba.exe

# complete list of all f90 source files (alphabetic order)
SRCS0 =	$(wildcard *.f90)
SRCS  = $(filter-out F%.f90, $(SRCS0))

# the object files are the same as the source files but with suffix ".o"
OBJS := $(SRCS:.f90=.o)

MAKEFILE_INC = depend.mk

# If you don't have sfmakedepend.pl, get it from:
# http://people.arsc.edu/~kate/Perl
F_makedepend = ./sfmakedepend --file=$(MAKEFILE_INC)

# do not delete $(ADDEFS), since it is set by Makefile.m
F90DEFS0 = $(OUTPUT) $(ADDEFS)
# mz_rs_20160219 additional CPP defs for tagging:
-include mecca/tag/f90defs0_mecca_tag.mk
F90DEFS  = $(addprefix $(DEFOPT), $(F90DEFS0))

.PHONY: all
all: $(PROG)

# the dependencies depend on the link
# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): depend $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) $(NETCDF_LIB) -o $@

# update file dependencies
depend $(MAKEFILE_INC): $(SRCS)
	$(F_makedepend) $(SRCS)

.PHONY: zip
zip:
	@./_zipcaaba.tcsh zip verbose

.PHONY: zipall
zipall:
	@./_zipcaaba.tcsh zipall verbose

# validate Fortran code (links, TABs and long lines):
.PHONY: validate
validate:
	@./_validate.tcsh

# forcheck:
.PHONY: forcheck
forcheck:
	@forchk -rigor -cond -f95 -obs -ff -decl -ext -intr -spec -ancmpl -anprg -anref -shcom -shinc -shmod -shprg -shref -shsrc -shsub -inf -plen 25 -pwid 132 *.f90 /soft/ECHAM5/lib/netcdf90.flb > forcheck.log 2>&1
	@echo "(forcheck found: 8=errors, 4=warnings, 2=infos)"

# put the same find commands here as in cmg command:
.PHONY: TAGS
TAGS:
	@F90FILES=`find . -name "*.f90" -type f` ;\
	 INCFILES=`find . -name "*.inc" -type f` ;\
	 NMLFILES=`find . -name "*.nml" -type f` ;\
	 TEXFILES=`find . -name "*.tex" -type f` ;\
	 KPPFILES1=`find . -name "*.eqn" -type f` ;\
	 KPPFILES2=`find . -name "*.spc" -type f` ;\
	 KPPFILES3=`find . -name "*.kpp" -type f` ;\
	 AWKFILES=`find . -name "*.awk" -type f` ;\
	 BATFILES=`find . -name "*.bat" -type f` ;\
	 RPLFILES=`find . -name "*.rpl" -type f` ;\
	 TBLFILES=`find . -name "*.tbl" -type f` ;\
	 PYFILES=`find . -name "*.py" -type f` ;\
	 JNLFILES=`find . -name "*.jnl" -type f` ;\
	 SCRIPTS=`find . -type f -perm -100 -not -regex ".*\.exe"` ;\
	 etags -l fortran $$F90FILES $$INCFILES $$NMLFILES \
           -lnone $$TEXFILES $$KPPFILES1 $$KPPFILES2 $$KPPFILES3 \
	   $$AWKFILES $$BATFILES $$RPLFILES $$TBLFILES $$PYFILES $$JNLFILES $$SCRIPTS

# list the configuration:
.PHONY: list
list:
	@echo "------------------------------------------------"
	@echo "SRCS           = $(SRCS)"
	@echo "------------------------------------------------"
	@echo "OBJS           = $(OBJS)"
	@echo "------------------------------------------------"
	@echo "OUTPUT         = $(OUTPUT)"
	@echo "SYSTEM         = $(SYSTEM)"
	@echo "HOST           = $(HOST)"
	@echo "COMPILER       = $(COMPILER)"
	@echo "F90            = $(F90)"
	@echo "BITS           = $(BITS)"
	@echo "F90FLAGS       = $(F90FLAGS)"
	@echo "NETCDF_INCLUDE = $(NETCDF_INCLUDE)"
	@echo "NETCDF_LIB     = $(NETCDF_LIB)"
	@echo "------------------------------------------------"

.PHONY: clean
clean:
	-rm -f depend.mk.old $(OBJS) *.mod *.log *.s *.i90 *.o *.optrpt *.pyc *~

.PHONY: distclean
distclean: clean
	-rm -f $(PROG)
	-rm -f depend.mk*
	-rm -f *.nc
	-rm -f *.dat
	-rm -f mecca/latex/*.aux
	-rm -f mecca/latex/*.blg
	-rm -f mecca/latex/*.dvi
	-rm -f mecca/latex/*.log
	-rm -f *.exe

messy_main_tracer_chemprop.inc: mecca/tracer/chemprop/messy_main_tracer_chemprop.tbl
	(cd mecca/tracer/chemprop; ./xchemprop)

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(F90) $(F90DEFS) $(F90FLAGS) $(NETCDF_INCLUDE) -c $<
