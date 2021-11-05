# -*- Makefile -*-

PREFIX = ../..

LIB =	$(PREFIX)/lib/libforpy.a

SRCSF90 = forpy_mod.F90

OBJS := $(SRCSF90:.F90=.o)
MODS := $(SRCSF90:.F90=.mod)

.SUFFIXES:
.SUFFIXES: .o .mod .F90

%.o: %.F90
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
