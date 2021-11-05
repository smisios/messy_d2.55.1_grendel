#
#
#

FC := gfortran
FFLAGS = -g -I ..
CFLAGS = -g -fno-leading-underscore

FSRCS := $(wildcard *.f90)
OBJS := $(patsubst %.f90,%.o,$(FSRCS)) 

all : tests

tests : f90_unix_dir_tst f90_unix_dirent_tst f90_unix_env_tst f90_unix_errno_tst f90_unix_file_tst f90_unix_io_tst f90_unix_proc_tst f90_unix_regexp_tst f90_unix_signal_tst f90_unix_time_tst

f90_unix_dir_tst : f90_unix_dir_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_dirent_tst : f90_unix_dirent_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_env_tst : f90_unix_env_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_errno_tst : f90_unix_errno_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_file_tst : f90_unix_file_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_io_tst : f90_unix_io_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_proc_tst : f90_unix_proc_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_regexp_tst : f90_unix_regexp_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_signal_tst : f90_unix_signal_tst.o f90_unix_signal_tst_ccode.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

f90_unix_time_tst : f90_unix_time_tst.o ../libposix90.a
	$(FC) $^ -lposix90 -L.. -o $@

clean :
	rm -f *.o *.mod *~ *_tst

dirty :
	touch *.f90

%.o : %.f90
	$(FC) $(FFLAGS) -c $<
