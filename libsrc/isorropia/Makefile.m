# -*- Makefile -*-

# ----------------------------------------------
####################################

# ----------------------------------------------

srcdir     = .
LIB        = ../../lib/libisorropia.a

# ----------------------------------------------
SRCS =  isocom.f \
	isofwd.f \
	isorev.f

OBJS := $(SRCS:.f=.o)

# ----------------------------------------------

all:	$(LIB)

install: all

$(LIB):	$(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
.SUFFIXES: $(SUFFIXES) .f

%.o: %.f
	$(FC) $(FFLAGS)  -c $<

clean:
	rm -f ./*.o ./*~ ./*.mod

distclean: clean
	rm -f ./$(LIB)

# ----------------------------------------------
