# -*- makefile-gmake -*-

OSRCS1  := $(wildcard ../src/*.F90)

PF := oas_

CHAN=MPI1
CPPDEF    = $(DEFOPT)use_comm_$(CHAN) $(DEFOPT)__VERBOSE $(DEFOPT)FLD3D

#not used:
#$(DEFOPT)TREAT_OVERLAY

#SRCS1	= $(shell ls *.F90)
SRCS1	:= $(subst ../src/,,${OSRCS1})

OBJS0	:= $(patsubst %.F90, %.o, $(SRCS1))
OBJS1   := $(OBJS0) GPTLget_memusage.o

VPATH 	= .:/../include:\
          ../include:\
          ../../scrip/src:\
          ../../mct/build/mct:\
          ../../mct/build/mpeu:\
	  $(NETCDF_INCLUDE):$(MPI_INCLUDE)

LIBRARY = ../../../../../lib/libpsmile.a

default: all

clean:
	rm -f i.* *.o *.mod

distclean:
	rm -f $(LIBRARY)
	rm -f *.F90 *.c

install: all

all:	$(LIBRARY)

#all:
#	@echo $(OSRCS1)
#	@echo $(SRCS1)

$(LIBRARY): $(SRCS1) $(OBJS1) 
	$(AR) $(ARFLAGS) $(LIBRARY) $(OBJS1) 

INCLS = -I. \
        -I../include \
        -I../../scrip/src \
        -I../../mct/build/mct \
        -I../../mct/build/mpeu \
        $(NETCDF_INCLUDE) $(MPI2_INCLUDE)
#        -I$(ARCHDIR)/build/lib/pio

INCLSC = -I../include

.SUFFIXES:
.SUFFIXES: .o .F90 .c

%.F90 :: ../src/%.F90
	sed "s/\(${PF}\)*mct_/${PF}mct_/g" ../src/$< > $@

%.o: %.F90
	$(F90) $(F90FLAGS) $(CPPDEF) $(INCLS)  -c   $<

%.o: %.c
	$(CC) $(CCFLAGS) $(CPPDEF) $(INCLSC) -c   $<

GPTLget_memusage.c: ../src/GPTLget_memusage.c
	sed "s/\(${PF}\)*mct_/${PF}mct_/g" $< > $@

GPTLget_memusage.o: GPTLget_memusage.c
	$(CC) $(CCFLAGS) $(CPPDEF) -DHAVE_SLASHPROC -c   $<

mod_psmile_io.o: mod_psmile_io.F90
	$(F90) $(F90FLAGS) $(INCLS)  -c   $<

#$(PF)mct_mod.o: ../../mct/build/mct/$(PF)mct_mod.o
#	cp -p ../../mct/build/mct/$(PF)mct_mod.o .

#
# ALL dependencies ...
#

GPTLget_memusage.o:
mod_oasis_kinds.o:

# mod_oasis_parameters.o: mod_oasis_kinds.o
# mod_oasis_data.o: mod_oasis_kinds.o 
# mod_oasis_sys.o: mod_oasis_kinds.o mod_oasis_data.o
# mod_oasis_mem.o: mod_oasis_kinds.o GPTLget_memusage.o mod_oasis_sys.o
# mod_oasis_timer.o: mod_oasis_kinds.o mod_oasis_sys.o mod_oasis_data.o
# mod_oasis_mpi.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o mod_oasis_timer.o
# mod_oasis_string.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o
# mod_oasis_namcouple.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o mod_oasis_mpi.o mod_oasis_string.o
# mod_oasis_part.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o mod_oasis_timer.o mod_oasis_mpi.o $(PF)mct_mod.o
# mod_oasis_var.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o mod_oasis_timer.o mod_oasis_mpi.o mod_oasis_part.o
# mod_oasis_ioshr.o:  mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_mpi.o mod_oasis_string.o $(PF)mct_mod.o
# mod_oasis_io.o:  mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o mod_oasis_ioshr.o mod_oasis_mpi.o $(PF)mct_mod.o
# mod_oasis_grid.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_io.o mod_oasis_part.o mod_oasis_mpi.o mod_oasis_timer.o $(PF)mct_mod.o
# mod_oasis_map.o: mod_oasis_kinds.o  mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o mod_oasis_mpi.o mod_oasis_var.o mod_oasis_part.o  \
#         mod_oasis_string.o mod_oasis_namcouple.o mod_oasis_timer.o mod_oasis_io.o  $(PF)mct_mod.o
# mod_oasis_coupler.o: mod_oasis_kinds.o  mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_map.o mod_oasis_parameters.o mod_oasis_mpi.o mod_oasis_var.o mod_oasis_part.o  \
#         mod_oasis_string.o mod_oasis_namcouple.o mod_oasis_timer.o mod_oasis_io.o \
#         mod_oasis_mem.o $(PF)mct_mod.o
# mod_oasis_advance.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o mod_oasis_var.o mod_oasis_part.o mod_oasis_mpi.o \
#         mod_oasis_coupler.o mod_oasis_timer.o mod_oasis_io.o mod_oasis_mem.o $(PF)mct_mod.o \
#         mod_oasis_map.o
# mod_oasis_method.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o mod_oasis_coupler.o mod_oasis_namcouple.o \
#         mod_oasis_timer.o mod_oasis_ioshr.o mod_oasis_advance.o mod_oasis_grid.o \
#         mod_oasis_mpi.o mod_oasis_part.o mod_oasis_var.o mod_oasis_mem.o $(PF)mct_mod.o
# mod_oasis_getput_interface.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
#         mod_oasis_parameters.o mod_oasis_var.o mod_oasis_advance.o $(PF)mct_mod.o
# mod_oasis_auxiliary_routines.o: mod_oasis_kinds.o mod_oasis_data.o mod_oasis_sys.o \
# 	mod_oasis_parameters.o mod_oasis_var.o mod_oasis_mpi.o \
# 	mod_oasis_coupler.o mod_oasis_timer.o mod_oasis_io.o $(PF)mct_mod.o
# mod_prism.o: mod_oasis_kinds.o mod_oasis_part.o mod_oasis_sys.o \
# 	mod_oasis_getput_interface.o mod_oasis_parameters.o \
# 	mod_oasis_grid.o mod_oasis_method.o mod_oasis_var.o
# mod_oasis.o: mod_oasis_kinds.o mod_oasis_part.o mod_oasis_sys.o \
# 	mod_oasis_getput_interface.o mod_oasis_parameters.o mod_oasis_auxiliary_routines.o \
# 	mod_oasis_grid.o mod_oasis_method.o mod_oasis_var.o

mod_oasis_advance.o : mod_oasis_advance.F90 mod_oasis_mpi.o mod_oasis_io.o mod_oasis_sys.o mod_oasis_var.o mod_oasis_timer.o mod_oasis_part.o mod_oasis_map.o mod_oasis_coupler.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_auxiliary_routines.o : mod_oasis_auxiliary_routines.F90 mod_oasis_mpi.o mod_oasis_io.o mod_oasis_sys.o mod_oasis_var.o mod_oasis_timer.o mod_oasis_coupler.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_coupler.o : mod_oasis_coupler.F90 mod_oasis_timer.o mod_oasis_io.o mod_oasis_string.o mod_oasis_mpi.o mod_oasis_var.o mod_oasis_part.o mod_oasis_map.o mod_oasis_sys.o mod_oasis_namcouple.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_data.o : mod_oasis_data.F90 mod_oasis_kinds.o 
mod_oasis.o : mod_oasis.F90 mod_oasis_sys.o mod_oasis_auxiliary_routines.o mod_oasis_grid.o mod_oasis_getput_interface.o mod_oasis_var.o mod_oasis_part.o mod_oasis_method.o mod_oasis_namcouple.o mod_oasis_parameters.o mod_oasis_kinds.o 
mod_oasis_getput_interface.o : mod_oasis_getput_interface.F90 mod_oasis_sys.o mod_oasis_var.o mod_oasis_advance.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_grid.o : mod_oasis_grid.F90 mod_oasis_timer.o mod_oasis_mpi.o mod_oasis_part.o mod_oasis_sys.o mod_oasis_io.o mod_oasis_data.o 
mod_oasis_io.o : mod_oasis_io.F90 mod_oasis_ioshr.o mod_oasis_sys.o mod_oasis_mpi.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_ioshr.o : mod_oasis_ioshr.F90 mod_oasis_mpi.o mod_oasis_string.o mod_oasis_sys.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_kinds.o : mod_oasis_kinds.F90 
mod_oasis_map.o : mod_oasis_map.F90 mod_oasis_timer.o mod_oasis_io.o mod_oasis_string.o mod_oasis_mpi.o mod_oasis_var.o mod_oasis_part.o mod_oasis_sys.o mod_oasis_namcouple.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_mem.o : mod_oasis_mem.F90 mod_oasis_sys.o mod_oasis_kinds.o 
mod_oasis_method.o : mod_oasis_method.F90 mod_oasis_string.o mod_oasis_mpi.o mod_oasis_grid.o mod_oasis_ioshr.o mod_oasis_timer.o mod_oasis_advance.o mod_oasis_coupler.o mod_oasis_var.o mod_oasis_part.o mod_oasis_namcouple.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_sys.o mod_oasis_mem.o mod_oasis_kinds.o 
mod_oasis_mpi.o : mod_oasis_mpi.F90 mod_oasis_timer.o mod_oasis_sys.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_namcouple.o : mod_oasis_namcouple.F90 mod_oasis_string.o mod_oasis_mpi.o mod_oasis_sys.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_parameters.o : mod_oasis_parameters.F90 mod_oasis_kinds.o 
mod_oasis_part.o : mod_oasis_part.F90 mod_oasis_timer.o mod_oasis_mpi.o mod_oasis_sys.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_string.o : mod_oasis_string.F90 mod_oasis_timer.o mod_oasis_sys.o mod_oasis_data.o mod_oasis_parameters.o mod_oasis_kinds.o 
mod_oasis_sys.o : mod_oasis_sys.F90 mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_timer.o : mod_oasis_timer.F90 mod_oasis_sys.o mod_oasis_data.o mod_oasis_kinds.o 
mod_oasis_var.o : mod_oasis_var.F90 mod_oasis_part.o mod_oasis_timer.o mod_oasis_mpi.o mod_oasis_sys.o mod_oasis_parameters.o mod_oasis_data.o mod_oasis_kinds.o 
mod_prism.o : mod_prism.F90 mod_oasis_sys.o mod_oasis_auxiliary_routines.o mod_oasis_grid.o mod_oasis_getput_interface.o mod_oasis_var.o mod_oasis_part.o mod_oasis_method.o mod_oasis_namcouple.o mod_oasis_parameters.o mod_oasis_kinds.o 
