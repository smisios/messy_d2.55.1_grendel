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
PROG = ../../../bin/mpiom.exe

ifeq ($(strip $(MESSY)), MESSY)
MPIOM_LIB     = $(SPEC_LIB) $(NETCDF_LIB) $(PNETCDF_LIB) -L../../../../lib -lmpiom 
MPIOM_INCLUDE = $(NETCDF_INCLUDE) $(PNETCDF_INCLUDE)

ifeq ($(strip $(MESSYMMD)), MESSYMMD)
MPIOM_LIB     += -L../../../../lib -lmmd
MPIOM_INCLUDE += -I../../../../libsrc/mmd/src
endif

# ---------------------------------------------------------------
SRCDIRS  = smcl smil bmil bml 
LIBSRCS  = smcl ../../../mpiom/src

INCLUDES = $(MODOPT)./ $(addprefix $(MODOPT)../, $(SRCDIRS)) \
	   $(addprefix $(MODOPT)../, $(LIBSRCS)) \
           $(MPIOM_INCLUDE)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
LIBS     = $(MPIOM_LIB)
#F90DEFS0 = MPI MBM_MPIOM MPIOM_13B MESSY LITTLE_ENDIAN 
F90DEFS0 = $(MPI) MBM_MPIOM $(MPIOMVERS) MESSY $(MESSYMMD) $(ENDIAN)
F90DEFS  = $(addprefix $(DEFOPT), $(F90DEFS0))
F90NOR8  = $(INCLUDES) $(F90FLAGS) $(F90DEFS)
F90ALL   = $(F90NOR8)
LDFLAGS  = $(F90NOR8)
# ---------------------------------------------------------------
else
define err_message

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mpiom-1.3.0-beta cannot be compiled with --disable-MESSY anymore !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

endef
endif
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
ifeq ($(strip $(MESSY)), MESSY)
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
else
	$(warning $(err_message))
endif

#Makefile : ;
#%.mk :: ;
.PHONY: install
install: force

.PHONY: clean
clean:
ifeq ($(strip $(MESSY)), MESSY)
	-find $(OBJDIR) -type f -name '*' -exec rm -f '{}' \;
	-rm -fr $(OBJDIR)
	-find . -type f -name '*.log' -exec rm -f '{}' \;
	-find . -type f -name '*.nc' -exec rm -f '{}' \;
	-find . -type f -name '*.mc' -exec rm -f '{}' \;
else
	$(warning $(err_message))
endif

# Note: "-rm -fr $(OBJDIR)" must be the last command
.PHONY: distclean
distclean: clean
ifeq ($(strip $(MESSY)), MESSY)
	rm -fr $(OBJDIR)
	rm -f  $(PROG)
else
	$(warning $(err_message))
endif

# ----------------------------------------------------------------------
# MISC
# ----------------------------------------------------------------------
.PHONY: check
check:
ifeq ($(strip $(MESSY)), MESSY)
	../../util/mfchk mpiom "$(F90DEFS0)" "$(SRCDIRS) $(LIBSRCS)" "$(ECHAM5_INCLUDE)" "$(F90)" "$(FCKLIBS)"
else
	$(warning $(err_message))
endif

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

messy_main_tracer_chemprop.inc: ../tracer/chemprop/messy_main_tracer_chemprop.tbl
	(cd ../tracer/chemprop; ./xchemprop)
# ----------------------------------------------------------------------

.PHONY: all
ifeq ($(strip $(MESSY)), MESSY)
all: $(REDIR)/$(PROG)
	@echo ''
	@echo 'The compilation has been successful.'
	@echo ''
else
	$(warning $(err_message))
endif

.PHONY: install
install: all

ifeq ($(strip $(MESSY)), MESSY)
ifeq ($(findstring $(strip $(COMPILER)), LF), )
$(REDIR)/$(PROG): depend $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
else
### for some reason, LF95 requires the re-build of mpiom ...
$(REDIR)/$(PROG): mpiom depend $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.PHONY: mpiom
mpiom:
	@echo '######################################################'
	@echo COMPILER IS $(COMPILER) - REBUILD OF mpiom REQUIRED ###
	@echo '######################################################'
	cd ../../../../mpiom/src ;\
	$(MAKE) -f Makefile.m clean ;\
	$(MAKE) -f Makefile.m ;\
	cd -
endif
else
	$(warning $(err_message))
endif

# update file dependencies
depend $(MAKEFILE_INC): $(SRCS)
ifeq ($(strip $(MESSY)), MESSY)
	$(F_makedepend)
else
	$(warning $(err_message))
endif

# check files
.PHONY: list
list:
ifeq ($(strip $(MESSY)), MESSY)
	@echo
	@echo "------------------------------------------------"
	@echo "SRCS = $(SRCS)"
	@echo "------------------------------------------------"
	@echo
	@echo "------------------------------------------------"
	@echo "OBJS = $(OBJS)"
	@echo "------------------------------------------------"
	@echo
else
	$(warning $(err_message))
endif

# .PHONY: DEFAULT
# .DEFAULT:
# 	@echo '  -------------------------------'
# 	@echo '  UNKNOWN TARGET '$(MAKECMDGOALS)
# 	@echo '  -------------------------------'

# ----------------------------------------------------------------------

include $(MAKEFILE_INC)
messy_mpiom_mem_e5.o: mo_mpi.o

-include $(REDIR)/bml/specific.mk
-include $(REDIR)/bmil/specific.mk
-include $(REDIR)/smil/specific.mk
-include $(REDIR)/smcl/specific.mk

# ----------------------------------------------------------------------
endif
# ----------------------------------------------------------------------
