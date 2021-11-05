# -*- Makefile -*-

PREFIX          = ../..

LIB =	$(PREFIX)/lib/libblas.a

SRCS =  $(shell ls *.f)

OBJS =  $(SRCS:.f=.o)

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS) 

%.o: %.f
	$(FC) $(FFLAGS) -c $<

install: all

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(LIB)
