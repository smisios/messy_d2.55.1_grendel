noinst_PROGRAMS =
noinst_LIBRARIES = libsupport.a
EXTRA_libsupport_a_SOURCES = 						\
	cdilib.c							\
	rtc_sx.s							\
	aix-fpe.c							\
	cfortran.h							\
	gribex.c

EXTRA_DIST = 								\
	cfortran.doc							\
	iso_varying_string.f90						\
	aix_fptrap.inc

libsupport_a_SOURCES=util_reshape.c util_convert.c util_pbio.c		\
	util_sysinfo.c util_system.c util_timer.c util_backtrace.c \
	mo_fpe.f90

AM_CPPFLAGS=$(NETCDF_INCLUDE) \
	-DHAVE_CF_INTERFACE $(TWS_LL_INCLUDE) \
	$(INTERNAL_CDI_CPPFLAGS) \
        -I$(top_srcdir)/support
libsupport_a_LIBADD = $(NEC_ES_SX_HIGH_RESOLUTION_TIMER_OBJ)	\
	$(FPTRAP_WRAPPER_OBJ) \
	$(INTERNAL_CDILIB) $(LIBOBJS)

if NEC_ES_SX_HIGH_RESOLUTION_TIMER
NEC_ES_SX_HIGH_RESOLUTION_TIMER_OBJ=rtc_sx.o
AM_CFLAGS=-pvctl,fullmsg,noassume,loopcnt=1000000 \
	  -Orestrict=all -Onooverlap
rtc_sx.o: rtc_sx.s
	$(CCAS) -c `test -f 'rtc_sx.s' || echo '$(srcdir)/'`rtc_sx.s
endif

if INTERNAL_CDI
INTERNAL_CDILIB=cdilib.o
endif

if HAVE_IBM_TWS_LOAD_LEVELER
noinst_PROGRAMS += node_usage
endif

if USE_AIX_FPTRAP
FPTRAP_WRAPPER_OBJ=aix-fpe.o
endif

node_usage_SOURCES = node_usage.c

node_usage_LDADD = @TWS_LL_LIB@

clean-local:
	-rm -f *.@FCMODEXT@ i.*.L
