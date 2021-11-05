# -*- makefile-gmake -*-

# ----------------------------------------------
####################################

# ----------------------------------------------

srcdir     = .
LIB        = ../../lib/libcrm.a
# ----------------------------------------------

# ----------------------------------------------

SOURCES   := $(shell cat Srcfiles)

CURDIR    :=$(shell pwd)

$(CURDIR)/Depends: $(CURDIR)/Srcfiles 
	./mkDepends -t . . Srcfiles > $@

$(CURDIR)/Srcfiles: 
	./mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES)))


install: all
all: $(LIB)

distclean: clean

$(LIB):	$(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)


clean: 
	rm -f ./$(LIB); rm -rf ./*.mod ./*.o ./Depends ./Srcfiles

.SUFFIXES:
.SUFFIXES: .F .F90 .c .s .o

.F.o:
	$(F77) -c $(FFLAGS) $(F90R8) $<

.F90.o:
	$(F90) -c $(F90FLAGS) $(F90R8) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

include $(CURDIR)/Depends


