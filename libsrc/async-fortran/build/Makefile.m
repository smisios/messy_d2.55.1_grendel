# -*- makefile-gmake -*-

export

SHELL = /bin/sh

REDIR = ..
SRCDIRS = src src/async_threads_cpp
ADDFLAGS =

PREFIX = ../../..

LIB =   $(PREFIX)/lib/libasyncf.a

vpath %.cpp $(REDIR)/src/async_threads_cpp
vpath %.F90 $(REDIR)/src

INCLUDES = -I. -I$(REDIR)/src/async_threads_cpp

.SUFFIXES:
.SUFFIXES: .cpp .cpp.o .o .F90

%.cpp.o: %.cpp
	$(MPICXX) $(MPICXXFLAGS) $(ADDFLAGS) $(CFLAGS) $(INCLUDES) -c -o $(@F) $<

%.o: %.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $(@F) $<

SRCS0   := $(foreach DIR, $(SRCDIRS), $(wildcard $(REDIR)/$(DIR)/*.cpp))
SRCS    := $(sort $(notdir $(SRCS0)))
OBJS    := $(SRCS:.cpp=.cpp.o)

SRCSF900 := $(foreach DIR, $(SRCDIRS), $(wildcard $(REDIR)/$(DIR)/*.F90))
SRCSF90  := $(sort $(notdir $(SRCSF900)))
OBJSF90  := $(SRCSF90:.F90=.o)

all: $(LIB)

$(LIB): $(OBJS) $(OBJSF90)
	echo $(SRCS) $(SRCSF90)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS) $(OBJSF90)

install: all

clean:
	rm -f $(OBJS) $(OBJSF90)

distclean: clean
	rm -f $(LIB)
	rm -f *.mod
