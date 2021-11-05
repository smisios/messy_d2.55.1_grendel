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
PROG = ../../../bin/vertico.exe

VERTICO_LIB     = $(NETCDF_LIB) $(PNETCDF_LIB)
VERTICO_INCLUDE = $(NETCDF_INCLUDE) $(PNETCDF_INCLUDE)

# ---------------------------------------------------------------
SRCDIRS  = smcl smil bmil bml
LIBSRCS  = smcl

INCLUDES = $(MODOPT)./ $(addprefix $(MODOPT)../, $(SRCDIRS)) \
	   $(addprefix $(MODOPT)../, $(LIBSRCS)) \
           $(VERTICO_INCLUDE)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
LIBS     = $(VERTICO_LIB)
F90DEFS0 = VERTICO NOMPI
# mz_rs_20160219 additional CPP defs for tagging:
-include ../caaba/mecca/tag/f90defs0_mecca_tag.mk
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

# op_pj_20150507+
MECCATAG_SMCL = $(shell find ../caaba/mecca/smcl -name 'messy_mecca_tag_*.f90' -print -or -name 'messy_mecca_tag_*.inc' -print)
MECCATAG_SMIL = $(shell find ../../smil -name 'messy_mecca_tag*_si.f90' -print)
# op_pj_20150507-
# op_pj_20150811+
MECCAPOLY_SMCL = $(shell find ../caaba/mecca/smcl -name 'messy_mecca[0-9][0-9][0-9]_kpp.f90' -print)
MECCAPOLY_SMIL = $(shell find ../../smil -name 'messy_mecca[0-9][0-9][0-9]_*_si.inc' -print)
# op_pj_20150811-

OBJDIR := __$(ARCH)

MAKETARGET = $(MAKE) --no-print-directory -C $(OBJDIR) -f $(CURDIR)/$(mfile) \
                  SRCDIR=$(CURDIR) $(MAKECMDGOALS)

#.PHONY: $(OBJDIR)
# $(OBJDIR):
#	+@[ -d $@ ] || mkdir -p $@
#	+@$(MAKETARGET)

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
# op_pj_20150507+
	@echo '-----------------------------------------------------------'
	@echo updating links to mecca tag files ...
	@find -L smcl -type l -print | xargs rm -f
	@find -L smil -type l -print | xargs rm -f
	@echo "SMCL: "$(MECCATAG_SMCL)
ifeq ($(MECCATAG_SMCL), )
	@echo " ... nothing required for SMCL"
else
	@cd smcl ;\
	 $(foreach f,$(MECCATAG_SMCL),ln -fs ../$(f) . ;)\
	cd -
endif
	@echo "SMIL: "$(MECCATAG_SMIL)
ifeq ($(MECCATAG_SMIL), )
	@echo " ... nothing required for SMIL"
else
	@cd smil ; \
	 $(foreach f,$(MECCATAG_SMIL),ln -s ../$(f) . ;)\
	cd -
endif
	@find smcl -type l -name 'messy_mecca_tag_*.f90' -print
	@find smcl -type l -name 'messy_mecca_tag_*.inc' -print
	@find smil -type l -name 'messy_mecca_tag_*_si.f90' -print
	@echo '-----------------------------------------------------------'
# op_pj_20150507-
# op_pj_20150811+
	@echo '-----------------------------------------------------------'
	@echo updating links to polymecca files ...
	@find -L smcl -type l -print | xargs rm -f
	@find -L smil -type l -print | xargs rm -f
	@echo "SMCL: "$(MECCAPOLY_SMCL)
ifeq ($(MECCAPOLY_SMCL), )
	@echo " ... nothing required for SMCL"
else
	@cd smcl ;\
	 $(foreach f,$(MECCAPOLY_SMCL),ln -fs ../$(f) . ;)\
	cd -
endif
	@echo "SMIL: "$(MECCAPOLY_SMIL)
ifeq ($(MECCAPOLY_SMIL), )
	@echo " ... nothing required for SMIL"
else
	@cd smil ; \
	 $(foreach f,$(MECCAPOLY_SMIL),ln -s ../$(f) . ;)\
	cd -
endif
	@find smcl -type l -name 'messy_mecca[0-9][0-9][0-9]_kpp.f90' -print
	@find smil -type l -name 'messy_mecca[0-9][0-9][0-9]_*_si.inc' -print
	@echo '-----------------------------------------------------------'
# op_pj_20150811-
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	@cd $(OBJDIR) ; [ -f main.mk ] || ln -s ../main.mk main.mk ; cd ..
	@$(MAKETARGET)
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
# op_pj_20150507+
	-find . -type l -name 'messy_mecca_tag_*.f90' -exec rm -f '{}' \;
	-find . -type l -name 'messy_mecca_tag_*.inc' -exec rm -f '{}' \;
	-find . -type l -name 'messy_mecca_tag_*_si.f90' -exec rm -f '{}' \;
# op_pj_20150507-
# op_pj_20150811+
	-find . -type l -name 'messy_mecca[0-9][0-9][0-9]_*.*' -exec rm -f '{}' \;	
# op_pj_20150811-


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
	../../util/mfchk vertico "$(F90DEFS0)" "$(SRCDIRS) $(LIBSRCS)" "$(ECHAM5_INCLUDE)" "$(F90)" "$(FCKLIBS)"

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
#	$(F_makedepend)

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
