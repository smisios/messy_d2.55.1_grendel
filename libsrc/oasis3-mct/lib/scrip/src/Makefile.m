# -*- Makefile -*-
# NAG compiler ###############################################################
ifeq ($(COMPILER), NAG)
 F90FLAGS := $(F90FLAGS) -mismatch
 FFLAGS := $(FFLAGS) -mismatch
endif

SRCS1 =  $(shell ls *.F90)
OBJS1 =  $(patsubst %.F90, %.o, $(SRCS1))

SRCS3   = $(shell ls *.F)
OBJS3   = $(patsubst %.F, %.o, $(SRCS3))

SRCS4 =  $(shell ls *.f)
OBJS4 =  $(patsubst %.f, %.o, $(SRCS4))

VPATH =	.:$(COUPLE)/lib/scrip/src:$(COUPLE)/src

LIBRARY	= ../../../../../lib/libscrip.a

clean:
	rm -f i.* *.o *.mod

distclean: clean
	rm -f $(LIBRARY)

install: all

all:	$(LIBRARY)

$(LIBRARY): $(OBJS1)  $(OBJS3) $(OBJS4)
	$(AR) $(ARFLAGS) $(LIBRARY) $(OBJS1) $(OBJS3) $(OBJS4)

.SUFFIXES:
.SUFFIXES: .o .f .F .f90 .F90 .c

%.o: %.F90
	$(F90) $(F90FLAGS) $(NETCDF_INCLUDE) -c   $<

%.o: %.f90
	$(F90) $(F90FLAGS) -c   $<

%.o: %.F
	$(F90) $(FFLAGS) $(NETCDF_INCLUDE) -c   $<

%.o: %.f
	$(F90) $(FFLAGS) -c   $<

%.o: %.c
	$(CC) $(CCFLAGS) -c   $<

#
# ALL dependencies ...
#
constants.o : constants.f kinds_mod.o 
distance.o : distance.f kinds_mod.o constants.o 
fracnnei.o : fracnnei.f mod_oasis_flush.o remap_vars.o grids.o constants.o kinds_mod.o 
fracnnei_vmm.o : fracnnei_vmm.f mod_oasis_flush.o remap_vars.o constants.o kinds_mod.o 
grids.o : grids.f mod_oasis_flush.o iounits.o constants.o kinds_mod.o 
iounits.o : iounits.f mod_oasis_flush.o kinds_mod.o 
kinds_mod.o : kinds_mod.f 
remap_bicubic.o : remap_bicubic.f mod_oasis_flush.o remap_vars.o grids.o constants.o kinds_mod.o 
remap_bilinear.o : remap_bilinear.f mod_oasis_flush.o remap_vars.o grids.o constants.o kinds_mod.o 
remap_bilinear_reduced.o : remap_bilinear_reduced.f mod_oasis_flush.o remap_vars.o grids.o constants.o kinds_mod.o 
remap_gauswgt.o : remap_gauswgt.f mod_oasis_flush.o remap_vars.o grids.o constants.o kinds_mod.o 
remap_write.o : remap_write.f mod_oasis_flush.o netcdf.o remap_vars.o grids.o constants.o kinds_mod.o 
timers.o : timers.f kinds_mod.o 
mod_oasis_flush.o : mod_oasis_flush.F90 kinds_mod.o 
remap_bicubic_reduced.o : remap_bicubic_reduced.F90 mod_oasis_flush.o remap_vars.o grids.o constants.o kinds_mod.o 
netcdf.o : netcdf.F constants.o kinds_mod.o 
remap_conserv.o : remap_conserv.F mod_oasis_flush.o remap_vars.o grids.o timers.o constants.o kinds_mod.o 
remap_distwgt.o : remap_distwgt.F mod_oasis_flush.o remap_vars.o grids.o constants.o kinds_mod.o 
remap_vars.o : remap_vars.F mod_oasis_flush.o grids.o constants.o kinds_mod.o 
scrip.o : scrip.F mod_oasis_flush.o remap_write.o remap_bicubic_reduced.o remap_bilinear_reduced.o remap_bicubic.o remap_bilinear.o remap_gauswgt.o remap_distwgt.o remap_conserv.o remap_vars.o grids.o timers.o iounits.o constants.o kinds_mod.o 

