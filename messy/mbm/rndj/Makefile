# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh
# ----------------------------------------------

srcdir     = .
bindir     = .
includedir = .
PROG       = rndj.exe

F90        = g95
#F90FLAGS   = -r8 -cpp -ftrace=full
F90FLAGS   = -r8 -cpp -ftrace=full -fbounds-check

#F90        = lf95
#F90FLAGS   = --dbl --chk a,e,s --ap -g --trap -O0 -Cpp
##F90FLAGS   = --dbl --ap -O3 -Cpp

#F90        = xlf95_r
#F90FLAGS   = -g -qfullpath -qcheck -qflttrap=ov:zero:inv:en:nanq -qrealsize=8 -q64 -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -O3 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qspillsize=32648 -qarch=auto -qtune=auto -Q -qhot -qxlf90=nosignedzero -bdatapsize:64k -bstackpsize:64k -btextpsize:64k
#F90FLAGS   = -g -qfullpath -qrealsize=8 -q64 -qsuppress=1518-061:1518-128 -qsuppress=1500-036 -O3 -qstrict -qMAXMEM=-1 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qspillsize=32648 -qarch=auto -qtune=auto -Q -qhot -qxlf90=nosignedzero -bdatapsize:64k -bstackpsize:64k -btextpsize:64k

# ----------------------------------------------
SRCS =  $(wildcard *.f90)

OBJS := $(SRCS:.f90=.o)

all: $(PROG)

install: all
	cp -pf $(srcdir)/$(PROG) $(bindir)/.

clean:
	rm -f $(srcdir)/*.o
	rm -f $(srcdir)/*.mod
	rm -f $(srcdir)/*.dat

distclean: clean
	rm -f $(srcdir)/$(PROG)
	rm -f $(bindir)/$(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(F90R8) -o $@ $(OBJS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(F90R8) -c $<

messy_main_blather.o : messy_main_blather.f90 messy_main_constants_mem.o
messy_main_constants_mem.o : messy_main_constants_mem.f90
messy_main_rnd.o : messy_main_rnd.f90 messy_main_rnd_mtw_ja.o messy_main_rnd_lux.o messy_main_rnd_mtw.o messy_main_tools.o messy_main_constants_mem.o
messy_main_rnd_lux.o : messy_main_rnd_lux.f90
messy_main_rnd_mtw.o : messy_main_rnd_mtw.f90
messy_main_rnd_mtw_ja.o : messy_main_rnd_mtw_ja.f90
messy_main_tools.o : messy_main_tools.f90 messy_main_blather.o messy_main_constants_mem.o
rndj.o : rndj.f90 messy_main_tools.o messy_main_rnd.o messy_main_constants_mem.o
rnd_jump.o : rnd_jump.f90 messy_main_tools.o messy_main_rnd.o messy_main_constants_mem.o
