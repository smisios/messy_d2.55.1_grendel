# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = ../../../bin
includedir = .
PROG       = cloud.exe

# ----------------------------------------------
SRCS =  messy_main_constants_mem.f90 \
 messy_main_tools.f90 \
 messy_main_tools_wiso.f90 \
 messy_main_blather.f90 \
 messy_cloud_ori.f90 \
 messy_cloud_droplet.f90 \
 messy_cloud_box.f90 \

OBJS := $(SRCS:.f90=.o)

all: $(PROG)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*~
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f $(srcdir)/fort.[0-9]*

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) -o $@ $(OBJS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

# ----------------------------------------------
messy_cloud_box.o : messy_cloud_box.f90 messy_cloud_droplet.o messy_main_tools.o messy_main_constants_mem.o
messy_cloud_droplet.o : messy_cloud_droplet.f90 messy_cloud_ori.o messy_main_tools.o messy_main_constants_mem.o
messy_cloud_ori.o : messy_cloud_ori.f90 messy_main_tools_wiso.o messy_main_tools.o messy_main_constants_mem.o
messy_main_blather.o : messy_main_blather.f90 messy_main_tools.o messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_constants_mem.o
messy_main_tools_wiso.o : messy_main_tools_wiso.f90 messy_main_constants_mem.o
# ----------------------------------------------
