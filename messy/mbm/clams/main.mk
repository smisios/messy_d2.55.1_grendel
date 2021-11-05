# -*- Makefile -*-
##############################################################################

###srcdir = .
###top_srcdir = ..
### VPATH = .
# prefix = .
# exec_prefix = ${prefix}

# bindir = ${exec_prefix}/bin
# sbindir = ${exec_prefix}/sbin
# libexecdir = ${exec_prefix}/libexec
# datadir = ${prefix}/share
# sysconfdir = ${prefix}/etc
# libdir = ${exec_prefix}/lib
# includedir = ${prefix}/include
# oldincludedir = /usr/include
# infodir = ${prefix}/share/info
# mandir = ${prefix}/share/man

# sharedstatedir = ${prefix}/com
# localstatedir = ${prefix}/var

# program_transform_name = s,x,x,

# ---------------------------------------------------------------
ifeq ($(strip $(ARCH)), )
  ARCH = $(SYSTEM)
  mfile = Makefile
else
  mfile = Makefile.m
endif
# ---------------------------------------------------------------
PROG     = ../../../bin/clams.exe

QHULLLIB=-L../../../../lib -lqhull

CLAMS_LIB     = $(ICGLIB) $(NETCDF_LIB) $(PNETCDF_LIB) $(QHULLLIB)
CLAMS_INCLUDE = $(NETCDF_INCLUDE) $(PNETCDF_INCLUDE) $(QHULLLIB_INCLUDE)

ifeq ($(strip $(POSIX90_DEF)), HAVE_POSIX90)
CLAMS_INCLUDE += $(MODOPT) ../../../../libsrc/posix90/src
CLAMS_LIB     += -L../../../../lib -lposix90
endif

ifeq ($(strip $(FORPY_DEF)), HAVE_FORPY)
CLAMS_INCLUDE += $(MODOPT) ../../../../libsrc/forpy
CLAMS_LIB     += -L../../../../lib -lforpy $(PYTHON_LIB)
endif


# ---------------------------------------------------------------
SRCDIRS  = smcl smil bmil bml
LIBSRCS  = smcl

INCLUDES = $(MODOPT)./ $(addprefix $(MODOPT)../, $(SRCDIRS)) \
	   $(addprefix $(MODOPT)../, $(LIBSRCS)) \
	   $(CLAMS_INCLUDE)

# ---------------------------------------------------------------

# ---------------------------------------------------------------
LIBS     = $(CLAMS_LIB)
#F90DEFS0 = MBM_CLAMS MESSY MPI
### mo MMD
ZF90DEFS0 = $(filter-out MESSYMMD,$(F90DEFS0)) MBM_CLAMS MESSY MPI

F90DEFS  = $(addprefix $(DEFOPT), $(ZF90DEFS0))
F90NOR8  = $(INCLUDES) $(F90FLAGS) $(F90DEFS)
F90ALL   = $(F90NOR8)
LDFLAGS  = $(F90NOR8)
# ---------------------------------------------------------------

# ----------------------------------------------------------------------
 ifeq (,$(filter __$(ARCH),$(notdir $(CURDIR))))
# FROM HERE UP TO 'else' RULES AND TARGETS ARE VALID IN THE
# BASE DIRECTORY OF THIS DISTRIBUTION
# ----------------------------------------------------------------------
.SUFFIXES:

OBJDIR := __$(ARCH)

MAKETARGET = --no-print-directory -C $(OBJDIR) -f $(CURDIR)/$(mfile) \
		  SRCDIR=$(CURDIR) $(MAKECMDGOALS)

#.PHONY: $(OBJDIR)
# $(OBJDIR):
#	+@[ -d $@ ] || mkdir -p $@
#	+@$(MAKE) $(MAKETARGET)

NOOBJTARG = clean distclean

.PHONY: force
force:
ifneq (,$(strip $(findstring $(MAKECMDGOALS), $(NOOBJTARG))))
	@echo 'Do not need (OBJDIR) for $(MAKECMDGOALS)'
else
	@echo 'Need (OBJDIR) for $(MAKECMDGOALS)'
	@if ! test -r Makefile ; then \
	  echo "Makefile does not exist!" ; \
	  echo "Recursion not possible!"; \
	  exit 1 ; \
	fi
	@echo $(OBJDIR)
	@pwd
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	@cd $(OBJDIR) ; [ -f main.mk ] || ln -s ../main.mk main.mk ; cd ..
	@touch $(mfile)
	@$(MAKE) $(MAKETARGET)
endif

#% :: force ; @:

#Makefile : ;
#%.mk :: ;
.PHONY: install
install: force

.PHONY: clean
clean:
	-find $(OBJDIR) -type f -name '*' -exec rm -f '{}' \;
	-rm -fr $(OBJDIR)
	-find . -type f -name '*.log' -exec rm -f '{}' \;
	-find . -type f -name '*.nc' -exec rm -f '{}' \;
	-find . -type f -name '*.mc' -exec rm -f '{}' \;

# Note: "-rm -fr $(OBJDIR)" must be the last command
.PHONY: distclean
distclean: clean
	rm -fr $(OBJDIR)
	rm -f  $(PROG)

# ----------------------------------------------------------------------
# MISC
# ----------------------------------------------------------------------
.PHONY: check
check:
	../../util/mfchk clams "$(ZF90DEFS0)" "$(SRCDIRS) $(LIBSRCS)" "$(ECHAM5_INCLUDE)" "$(F90)" "$(FCKLIBS)"

# ----------------------------------------------------------------------
 else
# FROM HERE UP TO 'endif' RULES AND TARGETS ARE VALID IF THIS Makefile
# IS CALLED BY ITSELF FROM WITHIN THE (ARCHITECTIRE SPECIFIC) 'OBJDIR'
# ----------------------------------------------------------------------

REDIR    = ..

### SEARCH PATHs
vpath %.f90 $(REDIR)/smcl
vpath %.f90 $(REDIR)/smil
vpath %.f90 $(REDIR)/bmil
vpath %.f90 $(REDIR)/bml

SRCS0   := $(foreach DIR, $(SRCDIRS), $(wildcard $(REDIR)/$(DIR)/*.f90))
SRCS    := $(sort $(notdir $(SRCS0)))
OBJS    := $(SRCS:.f90=.o)

MAKEFILE_INC = depend.mk
F_makedepend = $(REDIR)/sfmakedepend --modcase $(EXTMODCASE) \
	       --depend=obj --libdeps --file=$(MAKEFILE_INC) \
	       $(addprefix --srcdir $(REDIR)/, $(SRCDIRS)) \
	       $(addprefix -I$(REDIR)/, $(SRCDIRS)) -I$(REDIR)/include \
	       $(addprefix --modsearch $(REDIR)/, $(LIBSRCS))

# ----------------------------------------------------------------------

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90ALL) -c $<

# ----------------------------------------------------------------------

.PHONY: all
all: $(REDIR)/$(PROG)
	@echo ''
	@echo 'The compilation has been successful.'
	@echo ''

.PHONY: install
install: all

$(REDIR)/$(PROG): depend $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

# update file dependencies
depend $(MAKEFILE_INC): $(SRCS)
	$(F_makedepend)

# check files
.PHONY: list
list:
	@echo
	@echo "------------------------------------------------"
	@echo "SRCS = $(SRCS)"
	@echo "------------------------------------------------"
	@echo
	@echo "------------------------------------------------"
	@echo "OBJS = $(OBJS)"
	@echo "------------------------------------------------"
	@echo

# .PHONY: DEFAULT
# .DEFAULT:
#	@echo '  -------------------------------'
#	@echo '  UNKNOWN TARGET '$(MAKECMDGOALS)
#	@echo '  -------------------------------'

# ----------------------------------------------------------------------

include $(MAKEFILE_INC)
-include $(REDIR)/bml/specific.mk
-include $(REDIR)/bmil/specific.mk
-include $(REDIR)/smil/specific.mk
-include $(REDIR)/smcl/specific.mk

# ----------------------------------------------------------------------
endif
# ----------------------------------------------------------------------
