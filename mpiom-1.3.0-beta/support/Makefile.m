# -*- Makefile -*-

PREFIX = ../..

LIB =	$(PREFIX)/lib/libsupport.a

SRCS =	gribex.c  util_convert.c  util_pbio.c  util_sysinfo.c  util_system.c  util_timer.c

INCLUDES = -I../../config

ifeq ($(ARCH), SX)
SRCS =	gribex.c util_convert.c util_pbio.c util_sysinfo.c util_system.c util_timer.c
endif

OBJS := $(SRCS:.c=.o)

.SUFFIXES:
.SUFFIXES: .c .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS) 

install: all

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(LIB)

ifeq ($(ARCH), SX)
gribex.o: gribex.c
	$(CC) $(CFLAGS) $(INCLUDES) -pvctl,fullmsg,noassume,loopcnt=1000000 -Orestrict=all -Onooverlap -c gribex.c
rtc_sx.o: rtc_sx.s
	$(AS) -c rtc_sx.s
endif

ifeq ($(ARCH), ES)
gribex.o: gribex.c
	$(CC) $(CFLAGS) $(INCLUDES) -pvctl,fullmsg,noassume,loopcnt=1000000 -Orestrict=all -Onooverlap -c gribex.c
rtc_sx.o: rtc_sx.s
	$(AS) -c rtc_sx.s
endif
