# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

# ----------------------------------------------

srcdir     = ./src
bindir     = ../../../bin
includedir = ./include
PROG       = kp4.exe

# ----------------------------------------------

all: $(PROG)

$(PROG):
	cd $(srcdir) ; \
	$(MAKE)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*~
	rm -f $(srcdir)/*.o

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/kp4.exe

# ----------------------------------------------
