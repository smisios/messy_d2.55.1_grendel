#
# Posix90 top-level Makefile
#

PREFIX := /usr/local
MODDIR := $(PREFIX)/f90/mod
LIBDIR := $(PREFIX)/lib
INFODIR := $(PREFIX)/info

all : lib tests doc

install : 
	install -d $(MODDIR)
	install src/*.mod $(MODDIR)
	install -d $(LIBDIR)
	install src/libposix90.a  $(LIBDIR)
	install -d $(INFODIR)
	install doc/posix90.info $(INFODIR)

lib :
	(cd src ; make -k libposix90)

tests :
	(cd src/test ; make -k)

doc :
	( cd doc ; make -k)

clean :
	(cd doc ; make clean)
	(cd src ; make clean)
	(cd src/test ; make clean)
	rm *~

.PHONY : doc