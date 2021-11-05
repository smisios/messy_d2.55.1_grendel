# -*- Makefile -*-

# ----------------------------------------------
export
SHELL     = sh

F90       = $(FC)
F90FLAGS  = $(FCFLAGS)

####################################

# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
PROG       = isorropia.exe

# ----------------------------------------------
SRCS =  isocom.f \
	isofwd.f \
	isorev.f \
        main.f

OBJS := $(SRCS:.f=.o)

# ----------------------------------------------

all:	$(PROG)

install:	all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

$(PROG):  $(OBJS)
	$(F90) $(F90LAGS) $(LIB) -o $@ $(OBJS)
.SUFFIXES: $(SUFFIXES) .f

%.o: %.f
	$(F90) $(F90FLAGS)  -c $<

clean:
	rm -f ./*.o ./*~ ./*.mod ./$(PROG)

distclean: clean
	rm -fr ./config.log ./config.status ./Makefile
	rm -f $(bindir)/$(PROG)

# ----------------------------------------------
