# -*- Makefile -*-
# ----------------------------------------------

export

SHELL = /bin/sh

# ---------------------------------------------------------------

###srcdir = .
###top_srcdir = .
### VPATH = .
prefix = .
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${datarootdir}
sysconfdir = ${prefix}/etc
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include
oldincludedir = /usr/include
infodir = ${datarootdir}/info
mandir = ${datarootdir}/man
datarootdir = ${prefix}/share

sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var

program_transform_name = s,x,x,

# ---------------------------------------------------------------
# ARCHITECTURE
ARCH     = LINUX64
ENDIAN   = LITTLE_ENDIAN
LINKMODE = 

# TOOLS
AR       = ar
ARFLAGS  = crv
NMFLAGS  = 
AS       = as
MD5SUM_PRESENT = yes

# VERSION CONTROL
MESSYVCS = 

# COMPILERS AND FLAGS
COMPILER = 
CPP      = mpicc -E
CPPFLAGS = 
DEFS     = -DHAVE_CONFIG_H

CXX      = g++
MPICXX   = 
MPICXXFLAGS = 

CC       = mpicc
CFLAGS   = -O -fp-model strict -Df2cFortran

FC       = mpifort
FFLAGS   = -O3 -fpp -heap-arrays

F90      = mpifort
F90R8    = -autodouble
F90VERS  = ifort (IFORT) 19.1.1.217 20200306
F90FLAGS = -O3 -fpp -heap-arrays -fp-model strict -lpthread -save-temps -fno-alias -align all
DEFOPT   = -D
MODOPT   = -I
### mz_rs_20160217 THESE F90DEFS ARE ONLY FOR THE SMCL, NOT FOR THE SMIL:
F90DEFS0 =  MESSY LITTLE_ENDIAN _LINUX64  PNCREGRID MESSYMMD  MPIOM_13B         
### mz_rs_20160217 additional CPP defs for tagging:
-include messy/mbm/caaba/mecca/tag/f90defs0_mecca_tag.mk
F90DEFS  = $(addprefix $(DEFOPT), $(F90DEFS0))

# BOX MODEL OF int2cosmo no MMD and no I2CINC enabled.
F90DEFS0I2C =  MESSY LITTLE_ENDIAN _LINUX64  PNCREGRID MPIOM_13B 
ifneq (,$(filter $(strip $(COMPILER)),G95 LF))
   F90DEFS0I2C += IOMSG_NOT_AVAIL
endif
F90DEFSI2C  = $(addprefix $(DEFOPT), $(F90DEFS0I2C))

# LIBRARIES
LIBSRCS         = echam5/blas echam5/lapack echam5/support libsrc/crm libsrc/dwdlibgrib1/source libsrc/ifs libsrc/isorropia libsrc/mmd/src libsrc/qhull libsrc/rttov7_synsat_vector/src messy/smcl mpiom/src 
NETCDF_LIB      = -L/home/stergios/sw/netcdf-3.6.3-intel/lib -lnetcdf
NETCDF_INCLUDE  = -I/home/stergios/sw/netcdf-3.6.3-intel/include -I/home/stergios/sw/netcdf-3.6.3-intel/include
PNETCDF_LIB     = 
PNETCDF_INCLUDE = 
MPI             = 
MPI2_LIB        = 
MPI2_INCLUDE    = 
YAXTROOT        = 
YAXT_INCLUDE    = 
YAXT_LIB        = 
CDI_INCLUDE     = 
CDI_LIB         = 
CUDA_LIB        = 
POSIX90_INCLUDE = 
POSIX90_LIB     = 
POSIX90_DEF     = 
T8CODEROOT      = 
T8CODE_INCLUDE  = 
T8CODE_LIB      = 
T8CODE_DEF      = 
FORPY_INCLUDE   = 
FORPY_LIB       = 
FORPY_DEF       = 
PYTHON_LIB      = 
ASYNCF_INCLUDE  = 
ASYNCF_LIB      = 

# TARGETS
MESSY          = MESSY
MESSYMMD       = MESSYMMD
ECHAM5_LIB     =  -L../../lib -lmpiom -L../../lib -llapack -L../../lib -lblas -L../../lib -lsupport -L../../lib -lmessy -lcrm -lqhull -lisorropia  -L../../lib -lmmd -L/home/stergios/sw/netcdf-3.6.3-intel/lib -lnetcdf      
ECHAM5_INCLUDE =    -I../../libsrc/mmd/src -I/home/stergios/sw/netcdf-3.6.3-intel/include -I/home/stergios/sw/netcdf-3.6.3-intel/include     
CESM1_LIB      = 
CESM1_INCLUDE  = 
CLM_LIB        = 
CLM_INCLUDE    = 
MPIOMVERS      = MPIOM_13B
COSMO          = COSMO

### REVISION NUMBER (used in ./messy/smcl/messy_main_constants_mem.f90)
VCSREV=
ifeq ($(strip $(MESSYVCS)), GIT)
 GIT_VERSION := $(shell git --no-pager describe --tags --always --dirty --broken)
 GIT_BRANCH  := $(shell git branch | grep \* | cut -d ' ' -f2)
 GIT_COMMIT  := $(shell git rev-parse --verify HEAD)
 GIT_DATE    := $(firstword $(shell git --no-pager show --date=iso-strict --format="%ad" --name-only))
 BUILD_DATE  := $(shell date --iso=seconds)

 VCSREV := $(GIT_VERSION)_$(GIT_COMMIT)_$(GIT_DATE)_$(BUILD_DATE)
endif

ifneq ($(strip $(VCSREV)),)
#F90DEFS  += $(DEFOPT)_VCSREV_="'$(VCSREV)'"
 F90DEFS  += $(DEFOPT)_VCSREV_=\''$(VCSREV)'\'
endif
# ---------------------------------------------------------------

# ---------------------------------------------------------------
### MESSy TOOLS
# for testing purposes, it can be useful to create just a few tools:
#TOOLS    = kpp kp4 kpp1 jvpp edgar2nc ncdx biogen
#TOOLDIRS := $(addprefix messy/tools/, $(TOOLS))
TOOLDIRS  := $(sort $(wildcard messy/tools/*))
TOOLS     := $(subst messy/tools/,,$(TOOLDIRS))

### MESSy BASEMODELS
# for testing purposes, it can be useful to create just a few MBMs:
#MBMS     = random caaba
#MBMDIRS  := $(addprefix messy/mbm/, $(MBMS))
MBMDIRS   := $(sort $(wildcard messy/mbm/*))

### OTHER BASEMODELS
#BASEMODELS = cosmo echam5 mpiom cesm1
BASEMODELS = cosmo echam5 cesm1 icon clm fesom2
BASEDIRS := $(addprefix ./, $(BASEMODELS))

### MESSy documentation
DOCUDIRS := $(sort $(wildcard messy/docu/*))

comma:= ,
empty:=
space:= $(empty) $(empty)
#---------------------------------------------------------------

# ---------------------------------------------------------------
.SUFFIXES:

.PHONY: all
all: modlog basemodels

.PHONY: dist
dist: modlog tools mbm basemodels

.PHONY: forceall
forceall:
	@echo ' -------------------------------------------------------------'
	@echo ' TOUCHING messy_main_constants_mem.f90 TO FORCE RECOMPILATION '
	@echo ' -------------------------------------------------------------'
	@touch messy/smcl/messy_main_constants_mem.f90
	$(MAKE)

.PHONY: help
help:
	@echo ''
	@echo '  posible targets:'
	@echo '  -------------------------------------------------------'
	@echo '   $(MAKE) [t]                : build all libs/mod/obj'
	@echo '   $(MAKE) all                : (= $(MAKE))'
	@echo '   $(MAKE) forceall           : rebuild all'
	@echo ''
	@echo '   $(MAKE) dist               : build tools, mbm, and basemodels'
	@echo ''
	@echo '   $(MAKE) tools [t]          : build all tools'
	@echo '   $(MAKE) toolsclean [t]     : clean all tools'
	@echo '   $(MAKE) toolsdistclean [t] : distclean all tools'
	@echo ''
	@echo '   $(MAKE) docu [t]           : build all docus'
	@echo '   $(MAKE) docuclean [t]      : clean all docus'
	@echo '   $(MAKE) docudistclean [t]  : distclean all docus'
	@echo ''
	@echo '   $(MAKE) mbm [t]            : build all messy basemodels'
	@echo '   $(MAKE) mbmclean [t]       : clean all messy basemodels'
	@echo '   $(MAKE) mbmdistclean [t]   : distclean all messy basemodels'
	@echo ''
	@echo '   Note: [t] indicates the additional option'
	@echo '           target=... '
	@echo '         where ... is a comma-separated list of'
	@echo '         tools, mbms, basemodels, or docus respectively,'
	@echo '         the make process can be restricted to.'
	@echo ''
	@echo '   $(MAKE) libs               : build all libraries'
	@echo '   $(MAKE) libsclean          : clean all libraries'
	@echo '   $(MAKE) libsdistclean      : distclean all libraries'
	@echo ''
	@echo '   $(MAKE) purge              : delete all *~ files'
	@echo '   $(MAKE) clean              : delete libs/mod/obj'
	@echo '   $(MAKE) veryclean          : distclean except for workdir/*'
	@echo '   $(MAKE) distclean          : clean up distribution'
	@echo '   $(MAKE) meccalinks         : reset mecca link structure'
	@echo ''
	@echo '   $(MAKE) help               : output this help'
	@echo '   $(MAKE) list               : list configuration'
	@echo '   $(MAKE) TAGS               : create new TAGS / tags file for emacs / vi'
	@echo ''
	@echo '   $(MAKE) check              : run forcheck'
	@echo '   $(MAKE) messycheck         : check MESSy standard conformity'
	@echo '   $(MAKE) validate           : validate code'
	@echo ''
	@echo '   $(MAKE) zip                : zip source-code'
	@echo '   $(MAKE) zip1r              : zip 1D and rerun files'
	@echo '   $(MAKE) zipall             : zip complete directory'
	@echo '   $(MAKE) tarall             : tar complete directory'
	@echo '   $(MAKE) rar                : rar source-code'
	@echo '  -------------------------------------------------------'
	@echo ''

.PHONY: VCSTAG
VCSTAG:
	@echo $(VCSREV)

# list the configuration:
.PHONY: list
list:
	@echo "--- ARCHITECTURE ----------------------------------------------"
	@echo "MACHINE        = `uname -a`"
	@echo "ARCH           = $(ARCH)"
	@echo "ENDIAN         = $(ENDIAN)"
	@echo "LINKMODE       = $(LINKMODE)"
	@echo "--- TOOLS -----------------------------------------------------"
	@echo "AR             = $(AR)"
	@echo "ARFLAGS        = $(ARFLAGS)"
	@echo "NMFLAGS        = $(NMFLAGS)"
	@echo "AS             = $(AS)"
	@echo "MD5SUM_PRESENT = $(MD5SUM_PRESENT)"
	@echo "--- VERSION CONTROL -------------------------------------------"
	@echo "MESSYVCS       = $(MESSYVCS)"
	@echo "VCSREV         = $(VCSREV)"
	@echo "--- COMPILERS AND FLAGS ---------------------------------------"
	@echo "COMPILER       = $(COMPILER)"
	@echo "CPP            = $(CPP)"
	@echo "CPPFLAGS       = $(CPPFLAGS)"
	@echo "DEFS           = $(DEFS)"
	@echo "CXX            = $(CXX)"
	@echo "CC             = $(CC)"
	@echo "CFLAGS         = $(CFLAGS)"
	@echo "FC             = $(FC)"
	@echo "FFLAGS         = $(FFLAGS)"
	@echo "F90            = $(F90)"
	@echo "F90VERS        = $(F90VERS)"
	@echo "F90R8          = $(F90R8)"
	@echo "F90FLAGS       = $(F90FLAGS)"
	@echo "DEFOPT         = $(DEFOPT)"
	@echo "MODOPT         = $(MODOPT)"
	@echo "F90DEFS0       = $(F90DEFS0)"
	@echo "F90DEFS        = $(F90DEFS)"
	@echo "F90DEFS0I2C    = $(F90DEFS0I2C)"
	@echo "F90DEFSI2C     = $(F90DEFSI2C)"
	@echo "--- LIBRARIES -------------------------------------------------"
	@echo "LIBSRCS        = $(LIBSRCS)"
	@echo "NETCDF_LIB     = $(NETCDF_LIB)"
	@echo "NETCDF_INCLUDE = $(NETCDF_INCLUDE)"
	@echo "CDI_LIB        = $(CDI_LIB)"
	@echo "CDI_INCLUDE    = $(CDI_INCLUDE)"
	@echo "POSIX90_LIB    = $(POSIX90_LIB)"
	@echo "POSIX90_INCLUDE= $(POSIX90_INCLUDE)"
	@echo "MPI            = $(MPI)"
	@echo "MPI2_LIB       = $(MPI2_LIB)"
	@echo "MPI2_INCLUDE   = $(MPI2_INCLUDE)"
	@echo "YAXTROOT       = $(YAXTROOT)"
	@echo "YAXT_LIB       = $(YAXT_LIB)"
	@echo "YAXT_INCLUDE   = $(YAXT_INCLUDE)"
	@echo "OASIS3MCT_LIB  = $(OASIS3MCT_LIB)"
	@echo "OASIS3MCT_INCLUDE = $(OASIS3MCT_INCLUDE)"
	@echo "--- TARGETS --------------------------------------------------"
	@echo "MESSY          = $(MESSY)"
	@echo "ECHAM5_LIB     = $(ECHAM5_LIB)"
	@echo "ECHAM5_INCLUDE = $(ECHAM5_INCLUDE)"
	@echo "CESM1_LIB      = $(CESM1_LIB)"
	@echo "CESM1_INCLUDE  = $(CESM1_INCLUDE)"
	@echo "COSMO          = $(COSMO)"
	@echo "--------------------------------------------------------------"
	@echo "ENVIRONMENT:"
	@env
	@echo "--------------------------------------------------------------"
	@echo "--------------------------------------------------------------"
	@echo "SETTINGS:"
	@set
	@echo "--------------------------------------------------------------"
	@echo ''
	@for DIR in $(BASEDIRS) ;\
	  do \
	    if test -r $$DIR/Makefile ; then \
	      back=`pwd`; \
	      cd $$DIR ;\
	      echo "======================================================" ;\
	      echo "CONFIGURATION OF $$DIR" ;\
	      echo "======================================================" ;\
	      $(MAKE) -f Makefile list ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	      cd $$back ; \
	    fi ; \
	  done
	@echo ''

.PHONY: modlog
modlog:
	@if ! test -z $$MODULESHOME ; then \
	    if test -r $$MODULESHOME/init/sh ; then \
	       . $$MODULESHOME/init/sh ;\
	       timestamp=`date +"%Y%m%d_%H%M%S"` ;\
	       module list 2> \
		 module_compile.$$timestamp.log 1>&2 ; \
	    fi \
	else \
	   echo "no \$$MODULESHOME" ; \
	fi

.PHONY: compinfo
compinfo:
	 @if test -r messy/smcl/Makefile.m ; then \
	    back=`pwd`; \
	    cd messy/smcl ;\
	    $(MAKE) -f Makefile.m messy_main_compilerinfo_mem.f90 ;\
	   cd $$back ;\
	 fi

# ---------------------------------------------------------------------
# LIBRARIES
# ---------------------------------------------------------------------
# make libraries
.PHONY: libs
libs:
	@for DIR in $(LIBSRCS) ;\
	  do \
	    echo '' ;\
	    echo "### libs: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m install ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    cd $$back ; \
	  done
	@echo ''
	@echo 'The compilation has been successful.'
	@echo 'NOTES:'
	@echo '- The libraries are now in lib.'
	@echo ''

.PHONY: libsclean
libsclean:
	@for DIR in $(LIBSRCS) ;\
	  do \
	    echo '' ;\
	    echo "### libs: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m clean ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    find . -name 'i.*.L' -exec rm -f '{}' \; ; \
	    find . -name 'i.*.O' -exec rm -f '{}' \; ; \
	    find . -name 'F*.f' -exec rm -f '{}' \; ; \
	    find . -name 'F*.f90' -exec rm -f '{}' \; ; \
	    find . -name '*.i90' -exec rm -f '{}' \; ; \
	    find . -name '*.optrpt' -exec rm -f '{}' \; ; \
	    find . -name '*.d' -exec rm -f '{}' \; ; \
	    cd $$back ; \
	  done

.PHONY: libsdistclean
libsdistclean:
	@for DIR in $(LIBSRCS) ;\
	  do \
	    echo '' ;\
	    echo "### libs: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m distclean ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    if test -r Makefile.conf ; then \
	      rm -f Makefile.conf ; \
	    fi ; \
	    cd $$back ; \
	  done

# ---------------------------------------------------------------------
# TOOLS
# ---------------------------------------------------------------------
# make messy tools
.PHONY: tools
tools:
	@if [ "$$target" != "" ] ; then \
	    ZTOOLDIRS="$(addprefix messy/tools/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "tools: $$ZTOOLDIRS" ;\
	else \
	    ZTOOLDIRS="$(TOOLDIRS)" ;\
	fi ;\
	for DIR in $$ZTOOLDIRS ;\
	  do \
	    echo '' ;\
	    echo "### tools: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m install ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    cd $$back ; \
	  done
	@echo ''
	@echo 'The compilation has been successful.'
	@echo 'NOTES:'
	@echo '- Copies of tools executables are now in bin.'
	@echo '- Some tools need to be run from within the respective'
	@echo '  subdirectory (messy/tools/*).'
	@echo ''

# clean messy tools
.PHONY: toolsclean
toolsclean:
	@if [ "$$target" != "" ] ; then \
	    ZTOOLDIRS="$(addprefix messy/tools/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "tools: $$ZTOOLDIRS" ;\
	else \
	    ZTOOLDIRS="$(TOOLDIRS)" ;\
	fi ;\
	for DIR in $$ZTOOLDIRS ;\
	  do \
	    echo '' ;\
	    echo "### tools: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m clean ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    find . -name 'i.*.L' -exec rm -f '{}' \; ; \
	    find . -name 'i.*.O' -exec rm -f '{}' \; ; \
	    find . -name 'F*.f' -exec rm -f '{}' \; ; \
	    find . -name 'F*.f90' -exec rm -f '{}' \; ; \
	    find . -name '*.i90' -exec rm -f '{}' \; ; \
	    find . -name '*.optrpt' -exec rm -f '{}' \; ; \
	    find . -name '*.d' -exec rm -f '{}' \; ; \
	    cd $$back ; \
	  done

# distclean messy tools
.PHONY: toolsdistclean
toolsdistclean:
	@if [ "$$target" != "" ] ; then \
	    ZTOOLDIRS="$(addprefix messy/tools/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "tools: $$ZTOOLDIRS" ;\
	else \
	    ZTOOLDIRS="$(TOOLDIRS)" ;\
	fi ;\
	for DIR in $$ZTOOLDIRS ;\
	  do \
	    echo '' ;\
	    echo "### tools: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m distclean ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    cd $$back ; \
	  done

# ---------------------------------------------------------------------
# MBM MODELS
# ---------------------------------------------------------------------
# make mbm models
.PHONY: mbm
mbm: compinfo libs
	@if [ "$$target" != "" ] ; then \
	    ZMBMDIRS="$(addprefix messy/mbm/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "mbm: $$ZMBMDIRS" ;\
	else \
	    ZMBMDIRS="$(MBMDIRS)" ;\
	fi ;\
	for DIR in $$ZMBMDIRS ;\
	  do \
	    echo '' ;\
	    echo "### mbm: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m install ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    cd $$back ; \
	  done
	@echo ''
	@echo 'The compilation has been successful.'
	@echo 'NOTES:'
	@echo '- Copies of messy basemodel executables are now in bin.'
	@echo '- Some mbm models need to be run from within the respective'
	@echo '  subdirectory (messy/mbm/*).'
	@echo ''

# clean mbm models
.PHONY: mbmclean
mbmclean:
	@if [ "$$target" != "" ] ; then \
	    ZMBMDIRS="$(addprefix messy/mbm/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "mbm: $$ZMBMDIRS" ;\
	else \
	    ZMBMDIRS="$(MBMDIRS)" ;\
	fi ;\
	for DIR in $$ZMBMDIRS ;\
	  do \
	    echo '' ;\
	    echo "### mbm: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m clean ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    find \. -name 'i.*.L' -exec rm -f '{}' \; ; \
	    find \. -name 'i.*.O' -exec rm -f '{}' \; ; \
	    find \. -name 'F*.f' -exec rm -f '{}' \; ; \
	    find \. -name 'F*.f90' -exec rm -f '{}' \; ; \
	    find \. -name '*.i90' -exec rm -f '{}' \; ; \
	    find \. -name '*.optrpt' -exec rm -f '{}' \; ; \
	    find \. -name '*.d' -exec rm -f '{}' \; ; \
	    cd $$back ; \
	  done

# distclean mbm models
.PHONY: mbmdistclean
mbmdistclean:
	@if [ "$$target" != "" ] ; then \
	    ZMBMDIRS="$(addprefix messy/mbm/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "mbm: $$ZMBMDIRS" ;\
	else \
	    ZMBMDIRS="$(MBMDIRS)" ;\
	fi ;\
	for DIR in $$ZMBMDIRS ;\
	  do \
	    echo '' ;\
	    echo "### mbm: $$DIR ###" ;\
	    echo '' ;\
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m distclean ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    fi ; \
	    cd $$back ; \
	  done

# ---------------------------------------------------------------------
# DOCUS
# ---------------------------------------------------------------------
# make messy docu
.PHONY: docu
docu:
	@if [ "$$target" != "" ] ; then \
	    ZDOCUDIRS="$(addprefix messy/docu/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "docus: $$ZDOCUDIRS" ;\
	else \
	    ZDOCUDIRS="$(DOCUDIRS)" ;\
	fi ;\
	for DIR in $$ZDOCUDIRS ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	       echo '' ;\
	       echo "### docu: $$DIR ###" ;\
	       echo '' ;\
	       $(MAKE) -f Makefile.m docu ; status=$$? ; \
	       if [ $$status != 0 ] ; then \
		 echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	       fi ; \
	    fi ; \
	    cd $$back ; \
	  done
	@echo ''
	@echo 'The creation of the documentation has been successful.'
	@echo 'NOTES:'
	@echo '- Copies of pdfs are now in messy/docu/pdf.'
	@echo ''

# clean messy docu
.PHONY: docuclean
docuclean:
	@if [ "$$target" != "" ] ; then \
	    ZDOCUDIRS="$(addprefix messy/docu/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "docus: $$ZDOCUDIRS" ;\
	else \
	    ZDOCUDIRS="$(DOCUDIRS)" ;\
	fi ;\
	for DIR in $$ZDOCUDIRS ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	       echo '' ;\
	       echo "### docu: $$DIR ###" ;\
	       echo '' ;\
	       $(MAKE) -f Makefile.m clean ; status=$$? ; \
	       if [ $$status != 0 ] ; then \
		 echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	       fi ; \
	    fi ; \
	    cd $$back ; \
	  done

# distclean messy docu
.PHONY: docudistclean
docudistclean:
	@if [ "$$target" != "" ] ; then \
	    ZDOCUDIRS="$(addprefix messy/docu/, $(subst $(comma),$(space),$(target)))" ;\
	    echo "docus: $$ZDOCUDIRS" ;\
	else \
	    ZDOCUDIRS="$(DOCUDIRS)" ;\
	fi ;\
	for DIR in $$ZDOCUDIRS ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	       echo '' ;\
	       echo "### docu: $$DIR ###" ;\
	       echo '' ;\
	       $(MAKE) -f Makefile.m distclean ; status=$$? ; \
	       if [ $$status != 0 ] ; then \
		 echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	       fi ; \
	    fi ; \
	    cd $$back ; \
	  done

# ---------------------------------------------------------------------
# BASE MODELS
# ---------------------------------------------------------------------
# make libraries
.PHONY: basemodels
basemodels: libs
	@if [ "$$target" != "" ] ; then \
	    ZBASEDIRS="$(subst $(comma),$(space),$(target))" ;\
	    echo "basemodels: $$ZBASEDIRS" ;\
	else \
	    ZBASEDIRS="$(BASEDIRS)" ;\
	fi ;\
	for DIR in $$ZBASEDIRS ;\
	  do \
	    echo '' ;\
	    echo "### basemodels: $$DIR ###" ;\
	    echo '' ;\
	    if test -r $$DIR/Makefile ; then \
	      back=`pwd`; \
	      cd $$DIR ;\
	      $(MAKE) -f Makefile ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	      cd $$back ; \
	    else \
	      echo "$$DIR/Makefile not present" ; \
	    fi ; \
	  done
	@echo ''
	@echo 'The compilation has been successful.'
	@echo 'NOTES:'
	@echo '- Copies of the basemodel executables are now in bin.'
	@echo ''

.PHONY: basemodelsclean
basemodelsclean:
	@if [ "$$target" != "" ] ; then \
	    ZBASEDIRS="$(subst $(comma),$(space),$(target))" ;\
	    echo "basemodels: $$ZBASEDIRS" ;\
	else \
	    ZBASEDIRS="$(BASEDIRS)" ;\
	fi ;\
	for DIR in $$ZBASEDIRS ;\
	  do \
	    echo '' ;\
	    echo "### basemodels: $$DIR ###" ;\
	    echo '' ;\
	    if test -r $$DIR/Makefile ; then \
	      back=`pwd`; \
	      cd $$DIR ;\
	      $(MAKE) -f Makefile clean ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	      cd $$back ; \
	    fi ; \
	  done

.PHONY: basemodelsdistclean
basemodelsdistclean:
	@if [ "$$target" != "" ] ; then \
	    ZBASEDIRS="$(subst $(comma),$(space),$(target))" ;\
	    echo "basemodls: $$ZBASEDIRS" ;\
	else \
	    ZBASEDIRS="$(BASEDIRS)" ;\
	fi ;\
	for DIR in $$ZBASEDIRS ;\
	  do \
	    echo '' ;\
	    echo "### basemodels: $$DIR ###" ;\
	    echo '' ;\
	    if test -r $$DIR/Makefile ; then \
	      back=`pwd`; \
	      cd $$DIR ;\
	      $(MAKE) -f Makefile distclean ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	      cd $$back ; \
	    fi ; \
	  done

# ---------------------------------------------------------------------
# MISC
# ---------------------------------------------------------------------
# purge deletes all tilde files
.PHONY: purge
purge:
	-find . -name '*~' -exec rm -f '{}' \;

# clean deletes all
# - tilde files
# - files produced by the compiler, except for executables and libraries
# - log files
.PHONY: clean
clean: purge libsclean toolsclean mbmclean basemodelsclean docuclean
	-find . -name '*.log' -exec rm -f '{}' \;
	-find . -name 'fort.[0-9]*' -exec rm -f '{}' \;
	-find . -name '*.pyc' -exec rm -f '{}' \;

# veryclean makes distclean but leaves the working-directory untouched
.PHONY: veryclean
veryclean: savedata distclean restoredata

.PHONY: meccalinks
meccalinks:
	messy/util/mecca1_link.sh
	messy/util/mecca_link.sh

.PHONY: savedata
savedata:
	mv -f workdir workdir.save

.PHONY: restoredata
restoredata:
	mv -f workdir.save workdir

# distclean
# - includes clean; it further deletes
# - all executables
# - files created by configure
# - the content of the working-directory
.PHONY: distclean
distclean: clean libsdistclean toolsdistclean mbmdistclean \
	   basemodelsdistclean docudistclean
	@echo ''
	@echo '### final distclean ###'
	@echo ''
	-find . -name '*.exe' -exec rm -f '{}' \;
	-find . -name '*.a' -exec rm -f '{}' \;
	-find . -name 'F*.f' -exec rm -f '{}' \;
	-find . -name '*.i90' -exec rm -f '{}' \;
	-find . -name '*.optrpt' -exec rm -f '{}' \;
	-find . -name '*.i'   -exec rm -f '{}' \;
	-find . -name '*.lst' -exec rm -f '{}' \;
	-find . -name 'i.*.L' -exec rm -f '{}' \;
	-find . -name 'i.*.O' -exec rm -f '{}' \;
	-find . -name 'tmp_*' -exec rm -f '{}' \;
	-find . -name '*.d' -exec rm -f '{}' \;
	-rm -f config.cache config.status config/config.h
	-rm -f Makefile
	-rm -fr workdir/*
	-rm -f ./messy/util/*.pl
	-rm -f ./echam5*/Makefile
	-rm -f ./cesm1*/Makefile
	# -rm -f ./mpiom*/Makefile
	-rm -f ./cosmo*/Makefile
	-rm -f ./icon*/Makefile
	-rm -f ./fesom2*/Makefile
	-rm -f ./clm*/Makefile
	-rm -f ./messy/util/locate_f90.sh
	-rm -f ./messy/libsrc/mct/Makefile.conf
	@echo ''
	@echo 'Distribution is clean now!'
	@echo 'Restart with ./configure [options]'
	@echo ''

.PHONY: zip
zip:
	messy/util/zip.sh

.PHONY: zip1r
zip1r:
	messy/util/zip1r.sh

.PHONY: zipall
zipall:
	messy/util/zipall.sh

.PHONY: tarall
tarall:
	@DIR=`pwd | sed 's|/.*/||'` ; \
	if test -r $$DIR.tar ; then \
	 echo "tar file exists ! Please (re)move it first:" ; exit 1 ;\
	 echo "rm $$DIR.tar" ; exit 1 ;\
	fi ; \
	back=`pwd`; \
	cd .. ; \
	if test -r $$DIR.tar ; then \
	 echo "tar file exists ! Please (re)move it first:" ; exit 1 ;\
	 echo "rm ../$$DIR.tar" ; exit 1 ;\
	fi ; \
	tar cvf $$DIR.tar $$DIR ; \
	mv -f $$DIR.tar $$DIR/. ; \
	cd $$back

.PHONY: rar
rar:
	messy/util/rar.sh

.PHONY: check
check:
	@for DIR in $(BASEDIRS) ;\
	  do \
	   if test -r $$DIR/Makefile ; then \
	    back=`pwd`; \
	    cd $$DIR ;\
	      $(MAKE) -f Makefile check ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	    cd $$back ; \
	   fi ; \
	  done

.PHONY: messycheck
messycheck:
	@for MODEL in $(BASEMODELS) ;\
	  do \
	    messy/util/messy_check_all $$MODEL verbose ; \
	  done

.PHONY: validate
validate:
	@for MODEL in $(BASEMODELS) ;\
	  do \
	   messy/util/validate $$MODEL ;\
	  done

# put the same find commands here as in messy/util/mmg!
.PHONY: TAGS
TAGS:
	@rm -f TAGS ctags
	@find . -name "*.f90"         -type f | xargs etags -a -l fortran ;\
	 find . -name "*.f90-*"       -type f | xargs etags -a -l fortran ;\
	 find . -name "*.inc"         -type f | xargs etags -a -l fortran ;\
	 find . -name "*.nml"         -type f | xargs etags -a -l fortran ;\
	 find . -name "*.f"           -type f | xargs etags -a -l fortran ;\
	 find . -name "*.F"           -type f | xargs etags -a -l fortran ;\
	 find . -name "*.F90"         -type f | xargs etags -a -l fortran ;\
	 find . -name "*.awk"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.bash"        -type f | xargs etags -a -l none    ;\
	 find . -name "*.bat"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.sh"          -type f | xargs etags -a -l none    ;\
	 find . -name "*.tcsh"        -type f | xargs etags -a -l none    ;\
	 find . -name "x*" -perm -100 -type f | xargs etags -a -l none    ;\
	 find . -name "*.eqn"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.k"           -type f | xargs etags -a -l none    ;\
	 find . -name "*.kpp"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.spc"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.rpl"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.mk"          -type f | xargs etags -a -l none    ;\
	 find . -name "Makefile.*"    -type f | xargs etags -a -l none    ;\
	 find . -name "*.jnl"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.tbl"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.tex"         -type f | xargs etags -a -l none    ;\
	 find . -name "*.f90"         -type f | xargs ctags    --language-force=fortran ;\
	 find . -name "*.f90-*"       -type f | xargs ctags -a --language-force=fortran ;\
	 find . -name "*.inc"         -type f | xargs ctags -a --language-force=fortran ;\
	 find . -name "*.nml"         -type f | xargs ctags -a --language-force=fortran ;\
	 find . -name "*.f"           -type f | xargs ctags -a --language-force=fortran ;\
	 find . -name "*.F"           -type f | xargs ctags -a --language-force=fortran ;\

.PHONY: depend
depend:
	@for DIR in $(LIBSRCS) ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    if test -r Makefile.m ; then \
	      $(MAKE) -f Makefile.m depend ; \
	    fi ; \
	    cd $$back ; \
	  done
	@for DIR in $(BASEDIRS) ;\
	  do \
	    if test -r $$DIR/Makefile ; then \
	      back=`pwd`; \
	      cd $$DIR ;\
	      $(MAKE) -f Makefile depend ; status=$$? ; \
	      if [ $$status != 0 ] ; then \
		echo "Exit status from $$MAKE was $$status" ; exit $$status ; \
	      fi ; \
	      cd $$back ; \
	    else \
	      echo "$$DIR/Makefile does not exist" ; \
	    fi ; \
	  done

# ---------------------------------------------------------------------
# END
# ---------------------------------------------------------------------
