
.SUFFIXES: .o .c .F .f .f90 .F90 .a

F90=mpxlf95
FFLAGS=-O3 -q64 -qsuffix=cpp=f90 -qstrict
CC=xlc
CFLAGS=-O3 -q64 
AR=ar -X64


            LIB = ../../lib/libmmd.a

        MEMBERS = \
                  $(LIB)(MMD_handle_communicator.o) \
                  $(LIB)(MMD_MPI_wrapper.o)         \
                  $(LIB)(mmd_server.o)              \
                  $(LIB)(mmd_client.o)
	
.c.a:
	@echo $<
	$(CC) $(CFLAGS) -c $(DEF) $<
	$(AR) r $@ $%
	rm -f $%

.f.a:
	@echo $<
	$(F90) $(FFLAGS) -c $<
	$(AR) r $@ $%
	rm -f $%

.f90.a:
	@echo $<
	$(F90) $(FFLAGS) -c $<
	$(AR) r $@ $%
	rm -f $%

$(LIB): $(MEMBERS)

clean: 
	rm -f *.mod $(LIB)

$(LIB)(MMD_MPI_wrapper.o):$(LIB)(MMD_handle_communicator.o)
$(LIB)(mmd_client.o): $(LIB)(MMD_MPI_wrapper.o)
$(LIB)(mmd_server.o): $(LIB)(MMD_MPI_wrapper.o)
