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
PROG = ../../../bin/spe.exe

BLANK_LIB     = $(NETCDF_LIB) $(PNETCDF_LIB)
BLANK_INCLUDE = $(NETCDF_INCLUDE) $(PNETCDF_INCLUDE)

# ---------------------------------------------------------------
SRCDIRS  = smcl smil bmil bml
LIBSRCS  = smcl

INCLUDES = $(MODOPT)./ $(addprefix $(MODOPT)../, $(SRCDIRS)) \
	   $(addprefix $(MODOPT)../, $(LIBSRCS)) \
           $(BLANK_INCLUDE)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
LIBS     = $(BLANK_LIB)
F90DEFS0 = BLANK MBM_SPE NOMPI
F90DEFS  = $(addprefix $(DEFOPT), $(F90DEFS0))
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
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	@cd $(OBJDIR) ; [ -f main.mk ] || ln -s ../main.mk main.mk ; cd ..
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
	-find . -type f -name '*.log' -exec rm -f '{}' \;
	@echo ''
	@echo 'Distribution is clean now!'

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
	../../util/mfchk spe "$(F90DEFS0)" "$(SRCDIRS) $(LIBSRCS)" "$(ECHAM5_INCLUDE)" "$(F90)" "$(FCKLIBS)"

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
# 	@echo '  -------------------------------'
# 	@echo '  UNKNOWN TARGET '$(MAKECMDGOALS)
# 	@echo '  -------------------------------'

# ----------------------------------------------------------------------

include $(MAKEFILE_INC)
-include $(REDIR)/bml/specific.mk
-include $(REDIR)/bmil/specific.mk
-include $(REDIR)/smil/specific.mk
-include $(REDIR)/smcl/specific.mk

# ----------------------------------------------------------------------
endif
# ----------------------------------------------------------------------
