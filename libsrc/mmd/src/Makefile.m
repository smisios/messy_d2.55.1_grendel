# -*- Makefile -*-

PREFIX = ../../..

LIB =	$(PREFIX)/lib/libmmd.a

SRCSF90 = mmd_child.f90 mmd_handle_communicator.f90 mmd_mpi_wrapper.f90\
          mmd_parent.f90 mmd_utilities.f90 mmd_test.f90

SRCSC   = mmdc_child.c mmdc_parent.c mmdc_util.c

OBJS := $(SRCSF90:.f90=.o) $(SRCSC:.c=.o)
MODS := $(SRCSF90:.f90=.mod)

.SUFFIXES:
.SUFFIXES: .c .o .f90

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

%.o: %.f90
#	$(F90) $(FFLAGS) -c $<
	$(F90) $(F90FLAGS) -c $<

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS) 

install: all

clean:
	rm -f $(OBJS)
	rm -f $(MODS)

distclean: clean
	rm -f $(LIB)

mmd_handle_communicator.o: mmd_utilities.o
mmd_mpi_wrapper.o: mmd_handle_communicator.o
mmd_child.o: mmd_mpi_wrapper.o
mmd_parent.o: mmd_mpi_wrapper.o
mmd_test.o: mmd_mpi_wrapper.o
mmdc_util.o: mmdc_util.h
mmdc_child.o: mmdc_util.h
mmdc_parent.o: mmdc_util.h


mmd_child.o: mmd_child.f90 mmd_mpi_wrapper.o mmd_handle_communicator.o mmd_utilities.o 
mmd_handle_communicator.o: mmd_handle_communicator.f90 mmd_utilities.o 
mmd_mpi_wrapper.o: mmd_mpi_wrapper.f90 mmd_utilities.o mmd_handle_communicator.o 
mmd_parent.o: mmd_parent.f90 mmd_mpi_wrapper.o mmd_handle_communicator.o mmd_utilities.o 
mmd_utilities.o: mmd_utilities.f90 
mmd_test.o: mmd_mpi_wrapper.o  mmd_test.f90
mmd_partest.o: mmd_mpi_wrapper.o  mmd_partest.f90


