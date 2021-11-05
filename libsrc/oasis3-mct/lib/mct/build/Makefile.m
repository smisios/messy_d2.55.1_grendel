# -*- makefile-gmake -*-

SHELL = /bin/sh

include Makefile.conf

#libdir = $(abs_top_builddir)/../../../../../lib # op_pj_20140409
libdir = ../../../../../lib

SUBDIRS = $(MPISERPATH) $(MPEUPATH) $(MCTPATH)

# TARGETS
subdirs: copy
	@set -e; for dir in $(SUBDIRS); do \
	  cd $$dir;                \
	  $(MAKE) -f Makefile.m;    \
	  cd $(abs_top_builddir);  \
	done

clean:
	@set -e; for dir in $(SUBDIRS); do \
	  if test -d $$dir ; then \
	  cd $$dir;                \
	  $(MAKE) -f Makefile.m clean; \
	  cd $(abs_top_builddir);  \
	  fi ; \
	done

install: subdirs
	@set -e; for dir in $(SUBDIRS); do \
	  cd $$dir;                \
	  $(MAKE) -f Makefile.m install; \
	  cd $(abs_top_builddir);  \
	done

distclean: clean
	rm -fr mct
	rm -fr mpeu
	rm -fr mpi-serial
	rm -f Makefile.conf.in
	rm -f Makefile.conf
	rm -f configure
	rm -f config.status
	rm -f config.log
	rm -f config.h.in
	rm -f config.h
	rm -f $(libdir)/libmpeu.a
	@rm -f $(libdir)/libmpi-serial.a
	rm -f $(libdir)/liboasmct.a

examples: subdirs
	@cd $(EXAMPLEPATH) && $(MAKE)

copy:
	@if ! test -d ./mpeu ; then \
	   mkdir ./mpeu ;\
	fi ;
	@cp -pf ../mpeu/Makefile.m mpeu/.
	@if ! test -d ./mpi-serial ; then \
	   mkdir mpi-serial ;\
	fi ;
	@cp -pf ../mpi-serial/Makefile.m mpi-serial/.
	@if ! test -d ./mct ; then \
	   mkdir mct ; \
	fi
	@cp -pf ../mct/Makefile.m mct/.
