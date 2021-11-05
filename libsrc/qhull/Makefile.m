# -*- makefile-gmake -*-

PREFIX          = ../..

LIB            = $(PREFIX)/lib/libqhull.a

# OBJS in execution frequency order.  CFILES after qhull.c are alphabetical
OBJS = user.o global.o stat.o io.o geom2.o poly2.o \
       merge.o qhull.o geom.o poly.o qset.o mem.o unix.o

CFILES= unix.c qhull.c geom.c geom2.c global.c io.c mem.c merge.c poly.c \
        poly2.c qset.c stat.c user.c
HFILES= user.h qhull.h qhull_a.h geom.h io.h mem.h merge.h poly.h qset.h stat.h
DFILES= Announce.txt REGISTER.txt COPYING.txt README.txt Changes.txt
TFILES= rbox.txt qhull.txt
FILES=  Makefile rbox.c user_eg.c q_test q_egtest q_eg
MFILES= qhull.man rbox.man qh-man.htm qh-faq.htm \
        qh-rbox.htm qh-impre.htm qh-opt.htm qh-eg.htm qh-c.htm

ADDOPT= 

ifeq (,$(filter $(strip $(COMPILER)), PGI))
   ADDOPT= -fno-strict-aliasing
endif

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)

install: all

unix.o:   qhull.h user.h mem.h
qhull.o:  $(HFILES)
geom.o:   $(HFILES)
geom2.o:  $(HFILES)
global.o: $(HFILES)
io.o:     $(HFILES)
mem.o:    mem.h 
merge.o:  $(HFILES)
poly.o:   $(HFILES)
poly2.o:  $(HFILES)
qset.o:   qset.h mem.h 
stat.o:   $(HFILES)
user.o:   $(HFILES)

.c.o:
	$(CC) -c $(CCOPTS1) $(ADDOPT) $<

clean:
	rm -f $(OBJS) core

distclean: clean
	rm -f $(LIB)
