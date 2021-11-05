# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

LIBS      = $(NETCDF_LIB)
INCLUDES  = $(NETCDF_INCLUDE)

# ----------------------------------------------

srcdir     = ./src
bindir     = ../../../bin
includedir = ./include
PROG       = import_grid.exe

# ----------------------------------------------

all: $(PROG)

$(PROG):
	cd $(srcdir) ; \
	$(MAKE)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.
	if test ! -d $(includedir) ; then mkdir $(includedir) ; fi
	cp -pf $(srcdir)/*.mod $(includedir)/.

clean:
	rm -f $(srcdir)/*.~
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)
	rm -fr $(includedir)

# ----------------------------------------------
