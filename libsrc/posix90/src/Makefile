#
#
#

FC := gfortran
CFLAGS := -g -fno-leading-underscore
FFLAGS := -g -W -Wall -fbounds-check # -pedantic # -std=f95

FSRCS := $(wildcard *.f90)
CSRCS := $(wildcard *ccode.c)
OBJS := $(patsubst %.f90,%.o,$(FSRCS))  $(patsubst %.c,%.o,$(CSRCS))

all : $(OBJS)

libposix90: libposix90.a($(OBJS))

libposix90.a($(OBJS)): $(OBJS)
	$(AR) r $@ $*.o

test :
	(cd test ; $(MAKE) -k)

f90_unix_dir.f90 :  f90_unix_dir_const.inc f90_unix_errno.o

f90_unix_dir_const.inc : f90_unix_dir_const
	./f90_unix_dir_const > f90_unix_dir_const.inc

f90_unix_dir_const: f90_unix_dir_const.o
	$(CC) f90_unix_dir_const.o -o $@



f90_unix_dirent.f90 :  f90_unix_dirent_const.inc f90_unix_errno.o

f90_unix_dirent_const.inc : f90_unix_dirent_const
	./f90_unix_dirent_const > f90_unix_dirent_const.inc

f90_unix_dirent_const: f90_unix_dirent_const.o
	$(CC) f90_unix_dirent_const.o -o $@



f90_unix_errno.f90: f90_unix_errno_const.inc f90_unix_tools.mod

f90_unix_errno_const.inc : f90_unix_errno_const
	./f90_unix_errno_const > f90_unix_errno_const.inc

f90_unix_errno_const: f90_unix_errno_const.o
	$(CC) f90_unix_errno_const.o -o $@

f90_unix_errno_tst.o : f90_unix_errno.mod

f90_unix_errno_tst : f90_unix_errno_tst.o f90_unix_errno_ccode.o \
		f90_unix_errno.o f90_unix_tools.o 
	$(FC) $^ -o $@


f90_unix_env.f90: f90_unix_env_const.inc f90_unix_errno.mod f90_unix_time.mod f90_unix_env_ccode.o

f90_unix_env_const.inc : f90_unix_env_const
	./f90_unix_env_const > f90_unix_env_const.inc

f90_unix_env_const: f90_unix_env_const.o
	$(CC) f90_unix_env_const.o -o $@


f90_unix_file.f90: f90_unix_file_const.inc f90_unix_errno.mod f90_unix_file_ccode.o

f90_unix_file_const.inc : f90_unix_file_const
	./f90_unix_file_const > f90_unix_file_const.inc

f90_unix_file_const: f90_unix_file_const.o
	$(CC) f90_unix_file_const.o -o $@


f90_unix_io.f90 :  f90_unix_env.mod f90_unix_io_const.inc

f90_unix_io_const.inc : f90_unix_io_const
	./f90_unix_io_const > f90_unix_io_const.inc

f90_unix_io_const: f90_unix_io_const.o
	$(CC) f90_unix_io_const.o -o $@


f90_unix_proc.f90: f90_unix_proc_const.inc f90_unix_errno.mod f90_unix_proc_ccode.o

f90_unix_proc_const.inc : f90_unix_proc_const
	./f90_unix_proc_const > f90_unix_proc_const.inc

f90_unix_proc_const: f90_unix_proc_const.o
	$(CC) f90_unix_proc_const.o -o $@


f90_unix_regexp.f90: f90_unix_regexp_const.inc f90_unix_errno.mod f90_unix_regexp_ccode.o

f90_unix_regexp_const.inc : f90_unix_regexp_const
	./f90_unix_regexp_const > f90_unix_regexp_const.inc

f90_unix_regexp_const: f90_unix_regexp_const.o
	$(CC) f90_unix_regexp_const.o -o $@


f90_unix_signal.f90: f90_unix_signal_const.inc f90_unix_errno.mod f90_unix_signal_ccode.o

f90_unix_signal_const.inc : f90_unix_signal_const
	./f90_unix_signal_const > f90_unix_signal_const.inc

f90_unix_signal_const: f90_unix_signal_const.o
	$(CC) f90_unix_signal_const.o -o $@


f90_unix_time.f90: f90_unix_time_const.inc f90_unix_errno.mod f90_unix_time_ccode.o

f90_unix_time_const.inc : f90_unix_time_const
	./f90_unix_time_const > f90_unix_time_const.inc

f90_unix_time_const: f90_unix_time_const.o
	$(CC) f90_unix_time_const.o -o $@


clean :
	rm -f *.o *.mod *.inc libposix90.a *_const *f90~ *c~

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<
