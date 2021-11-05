dnl Process this file with autoconf to produce a configure script.
dnl
dnl Thomas Jahns, January 30 2009
dnl  - added automake support and standard Fortran support (FC macros
dnl    instead of home-grown F90)
dnl Luis Kornblueh August 31 1998
dnl Uwe Schulzweida November 4 2002
dnl  -add cross compiling for NEC (use --host=sx6)
dnl  -use autoconf version 2.52
dnl  -add cross compiling for Earth Simulator, Yokohama (use --host=es)
dnl
AC_PREREQ([2.64])
AC_INIT([mpiom],[1.2.0],[haak@dkrz.de],[mpiom.1.2.0])
AC_CONFIG_SRCDIR([src/mpiom.f90])
AC_CONFIG_AUX_DIR([config])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([1.11 foreign])
AM_MAINTAINER_MODE
dnl
dnl Checks for package options.
dnl ######################################################################


dnl ######################################################################

AC_CONFIG_HEADER([config/config.h])

AC_ARG_ENABLE([mh-file],
  [AC_HELP_STRING([--enable-mh-file=<name>],dnl
   [select specific build and host configuration @<:@default: auto@:>@])],dnl
   ,[enable_mh_file=auto])
AC_ARG_WITH([setup],
  [AC_HELP_STRING([--with-setup=<name>],dnl
     [use user specific setup (mh_setup)])],,)
AS_IF([test x"$enable_mh_file" = xauto],
  [host_frag=""
   AS_CASE([$host],
     [alpha-*-*], [host_frag=config/mh-alpha],
     [5000-fujitsu-*], [host_frag=config/mh-fujitsu-5000],
     [*-fujitsu-*], [host_frag=config/mh-fujitsu],
     [sr8000-hitachi-*], [host_frag=config/mh-hitachi-sr8000],
     [*90-cray-*], [host_frag=config/mh-cray-c90],
     [sv1-cray-*], [host_frag=config/mh-cray-c90],
     [nv1-cray-*], [host_frag=config/mh-cray-x1],
     [t3e-cray-*], [host_frag=config/mh-cray-t3e],
     [xd1-cray-*], [host_frag=config/mh-cray-xd1],
     [xt3-cray-*], [host_frag=config/mh-cray-xt3],
     [mips*-sgi-irix6.*], [host_frag=config/mh-irix6],
     [*-*-solaris2*], [host_frag=config/mh-solaris],
     [sx*-*-*], [host_frag=config/mh-nec],
     [es*-*-*], [host_frag=config/mh-es],
     [*-ibm-aix*], [host_frag=config/mh-aix],
     [ia64-*-linux-*], [host_frag=config/mh-linux-ia64 ],
     [x86_64-*-linux-*], [host_frag=config/mh-linux-x64],
     [i*86-*-linux-*], [host_frag=config/mh-linux],
     [*-apple-darwin*], [host_frag=config/mh-darwin],
     [*-*-hpux*], [host_frag=config/mh-hpux],
     [AC_MSG_WARN([This configuration is not supported.
Please create a valid config/mh-* file.])])],dnl
  [host_frag="config/$enable_mh_file"])
dnl
AS_IF([test x"$enable_mh_file" != xno -a x"$enable_mh_file" != xoff],
  [AS_IF([test x"$host_frag" != x -a ! -f $srcdir/$host_frag],
     [AC_MSG_FAILURE([machine dependent configuration file $host_frag \
does not exist!], 1)])

   AS_IF([test "$ARCH" = CRAY_XT3],
     [cross_compiling=yes])

   mh_setup=""
   AS_IF([test x"$with_setup" != x],
     [AC_MSG_CHECKING([for user specific setup])
      AC_MSG_RESULT($with_setup)
      mh_setup=$with_setup
   ])
   AC_GET_MH([$srcdir/$host_frag])])
AC_ARG_WITH([fc-override],dnl
[AC_HELP_STRING([--with-fc-override=fcprg],dnl
   [override Fortran compiler selected by mh-file
    (i.e. to use instrumenting variant)])],
   [FC="$with_fc_override"])
AC_ARG_WITH([cc-override],dnl
[AC_HELP_STRING([--with-cc-override=ccprg],dnl
   [override C compiler selected by mh-file
    (i.e. to use instrumenting variant)])],
   [CC="$with_cc_override"])
dnl ######################################################################
dnl Checks for programs.
dnl ######################################################################
CHECK_GNU_MAKE
dnl we prefer GNU make because the System V make on many commercial
dnl platforms is slow
AC_CHECK_PROGS([MAKE], [$_cv_gnu_make_command gmake gnumake make])
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MKDIR_P
AC_PROG_CC
AC_PROG_CPP
AC_PROG_GREP
AC_PROG_FC
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_LANG_PUSH([Fortran])
AC_PROG_FPP
AC_LANG_POP([Fortran])
ACX_FC_POINTER_INTENT
AM_PROG_AS
AS_CASE([$host],
  [sx*-*-*], NEC_ES_SX_HIGH_RESOLUTION_TIMER=yes,
  [es*-*-*], NEC_ES_SX_HIGH_RESOLUTION_TIMER=yes)
AM_CONDITIONAL([NEC_ES_SX_HIGH_RESOLUTION_TIMER],
  [test "x$NEC_ES_SX_HIGH_RESOLUTION_TIMER" = xyes])
AC_PROG_RANLIB
AX_WITH_PERL([/bin/false], $PATH:$prefix/bin:/usr/bin:/usr/local/bin)
AX_PROG_PERL_VERSION([5.005],[],
        AC_MSG_ERROR([Perl (version>=5.0005) is required]))
dnl ######################################################################
dnl Find climate data operators version 1.x.
dnl ######################################################################
AC_CACHE_CHECK([for climate data operators binary], [ac_cv_path_CDO],
  [AC_PATH_PROGS_FEATURE_CHECK([CDO],[cdo],[$ac_path_CDO -V 2>&1 | \
      grep '^Climate Data Operators version 1' >conftest.cdo.log 2>&1 \
      && ac_cv_path_CDO=$ac_path_CDO ac_path_CDO_found=:],dnl
      AC_MSG_WARN([No functional version of cdo found.
      Setup-scripts will expect plain 'cdo' binary to be present.
      Running the setup-scripts might fail therefor.])
     [ac_cv_path_CDO=cdo])])
AC_SUBST([CDO], [$ac_cv_path_CDO])
dnl ######################################################################
dnl Find generic mapping tools (GMT) package.
dnl
dnl We require installed versions of pscoast and psxy.
dnl ######################################################################
AC_CACHE_CHECK([for gmt package binaries], [ac_cv_path_GMT],
  [AS_IF([test x${GMT_BINDIR+set} != xset -a x${GMT_ROOT+set} = xset],dnl
     [GMT_BINDIR="${GMT_ROOT}/bin"])
   AC_PATH_PROGS_FEATURE_CHECK([GMT],[GMT],dnl
     [GMT_EXECUTABLES=`PAGER=cat $ac_path_GMT 2>&1 | \
      grep '^@<:@ 	@:>@*executables' | \
      sed 's/^@<:@ 	@:>@*executables@<:@ 	@:>@*//'`
      AS_IF([test -x $GMT_EXECUTABLES/psxy],dnl
        [ac_cv_path_GMT=$GMT_EXECUTABLES/GMT ac_path_GMT_found=:],dnl
        [AS_IF([test -x `AS_DIRNAME([$ac_path_GMT])`/psxy],dnl
           [ac_cv_path_GMT=$ac_path_GMT ac_path_GMT_found=:])])],,dnl
[$GMT_BINDIR:$PATH])])
AS_IF([test x${ac_cv_path_GMT+set} = xset],dnl
  [GMT_BINDIR=`AS_DIRNAME([$ac_cv_path_GMT])`],dnl
  [AC_MSG_WARN([No functional install of GMT programs found.
      Some diagnostics during the setup scripts might not be produced.])])
AC_SUBST([GMT_BINDIR])
dnl ######################################################################
dnl Determine queuing system
dnl ######################################################################
AC_CHECK_PROGS([QUEUE_SUBMIT], [qsub llsubmit], [])
AM_SUBST_NOTMAKE(QUEUE_SUBMIT)
AS_CASE([$QUEUE_SUBMIT],
  [qsub],
    [AC_CACHE_CHECK([for queuing system type], [ac_cv_path_QUEUE_TYPE],
      [AC_PATH_PROGS_FEATURE_CHECK([QUEUE_TYPE], [$QUEUE_SUBMIT],
        [AS_IF(
          [$QUEUE_SUBMIT -help 2> /dev/null | head -1 | grep 'GE' > /dev/null],
          [ac_cv_path_QUEUE_TYPE=SGE],
          [ac_cv_path_QUEUE_TYPE=NQS]
         )]
       )]
     )],
  [llsubmit],
    [ac_cv_path_QUEUE_TYPE=TWS],
  dnl *default*
    [ac_cv_path_QUEUE_TYPE=none]
)
AC_SUBST([QUEUE_TYPE], [$ac_cv_path_QUEUE_TYPE])
AM_SUBST_NOTMAKE(QUEUE_TYPE)
dnl ######################################################################


dnl ######################################################################
dnl Properties of the compiler
dnl ######################################################################
AM_PROG_CC_C_O
AC_PROG_FC_C_O
AC_GET_FORTNAMES
ACX_SL_FC_MOD_SUFFIX(,[F90MODEXT=mod])
AC_LANG_PUSH([Fortran])
ACX_SL_FC_CHECK_MOD_PATH_FLAG
ACX_FC_FLUSH_CALL
AC_LANG_POP([Fortran])

dnl ######################################################################
dnl Properties of the compile environment checked unconditionally
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h fortran.h limits.h malloc.h netdb.h pwd.h \
 unistd.h sys/param.h time.h sys/time.h sys/unistd.h sys/utsname.h \
 inttypes.h sys/types.h])

AS_IF([test x${PLATFORM_OPT_LIBM+set} = xset],
  [ACX_ALTERNATIVE_LIBM([[$PLATFORM_OPT_LIBM]])],
  [ACX_ALTERNATIVE_LIBM])
dnl ######################################################################
dnl Optional feature arguments
dnl ######################################################################
dnl
dnl ######################################################################
dnl Checks for MPI.
dnl ######################################################################

AC_ARG_ENABLE([MPI],
  [AC_HELP_STRING([--enable-MPI],dnl
   [compile for and use MPI environment @<:@default: yes@:>@])])

AS_IF([test x"$enable_MPI" = x"no"],
  [MPI="${FPP_DEFOPT}NOMPI"
   MPI_LIB=""],
  [ACX_FORTRAN_PACKAGE([MPI], [mpif.h],,,dnl
    [AC_MSG_FAILURE([Required include mpif.h not found or not compilable.])],dnl
[mpi_waitall],[mpi mpi_f90 mpi_f77 mpich],dnl
ACX_M4_GENERATE_SUBSETS([[-lmpi_f77],[-lmpi]],[ ]),,dnl
    [AC_MSG_FAILURE([Cannot link MPI programs.])],[])
   ACX_C_PACKAGE([MPI], [mpi.h],,,dnl
      [],[MPI_Waitall],[mpi mpich],,,dnl
      [AC_MSG_NOTICE([Cannot link MPI programs.])],[])
   AS_IF([test $have_MPI_c_bindings = no],
     [AC_MSG_NOTICE([MPI C bindings unavailable, this might affect parallel libraries.])])
  ])
dnl
dnl
dnl ######################################################################
dnl  check wether we are building with coupling code
dnl ######################################################################
dnl  OASIS4 support:
dnl ######################################################################
AC_ARG_WITH([coupler-oasis4],
  [AS_HELP_STRING([--with-coupler-oasis4],
    [support for coupled models with the OASIS4 coupler, @<:@default=no@:>@])],
  [with_coupler_oasis4=yes],
  [with_coupler_oasis4=no])
AM_CONDITIONAL([COUPLE_OASIS4], [test x$with_coupler_oasis4 = xyes])
dnl ######################################################################
dnl  OASIS3 support:
dnl ######################################################################
AC_ARG_ENABLE([PRISM],
  [AC_HELP_STRING([--enable-PRISM],dnl
   [use PRISM for coupled model run @<:@default: no@:>@])])

PRISM_INCLUDE=""

AS_IF([test "$enable_PRISM" = "yes"],
  [PRISM='$(FPP_DEFOPT)use_comm_MPI1 $(FPP_DEFOPT)__prism \
	$(FPP_DEFOPT)__coupled $(FPP_DEFOPT)__CLIM_Box'
   AS_IF([test -n "$PRISMROOT"],
     [PRISM_LIB='-L$(PRISMROOT)/lib -lpsmile.MPI1 -lmpp_io '$PRISM_LIB
#      PRISM_INCLUDE='-I$(PRISMROOT)/include '
      PRISM_INCLUDE='$(FC_MOD_FLAG)$(PRISMROOT)/build/lib/psmile.MPI1'],
     [PRISM=""
      PRISM_LIB=""])],
  [PRISMROOT=""
   PRISM=""
   PRISM_LIB=""
   PRISM_INCLUDE=""])
AC_SUBST([PRISMROOT])
AC_SUBST([PRISM_LIB])
AC_SUBST([PRISM_INCLUDE])
AC_SUBST([PRISM])
dnl
dnl ######################################################################
dnl
dnl
dnl ######################################################################
dnl Choose method for surface forcing input (other than coupling)
dnl ######################################################################
AS_IF([test x$with_coupler_oasis4 != xyes],
  [AC_ARG_WITH([forcing-method],
     [AS_HELP_STRING([--with-forcing-method],
       [enable support for given forcing method (omip|core), @<:@default=omip@:>@])],,
     [with_forcing_method=omip])
   AS_CASE([$with_forcing_method],
     [omip],[FORCING_FC_DEFS=""],
     [ncep],[FORCING_FC_DEFS="\$(FPP_DEFOPT)CORE"],
     [AC_MSG_FAILURE([Unsupported surface forcing method requested.])])
   AC_MSG_NOTICE([Using $with_forcing_method surface forcing.])
])
AC_SUBST([FORCING_FC_DEFS])
dnl ######################################################################
dnl Checks for NETCDF.
dnl ######################################################################
dnl   NetCDF might use HDF5 on some platforms, which in turn might
dnl   require szip and almost definitely requires zlib
dnl ----------------------------------------------------------------------
dnl   test for ZLIB compression library
dnl
ACX_C_PACKAGE([zlib],[zlib.h],,,
  [AC_MSG_NOTICE([zlib library header not found, zlib compression not supported])],[deflate],[z],,,dnl
  [AC_MSG_NOTICE([cannot link zlib library, zlib compression not supported])])
AS_IF([test $have_zlib_c_bindings = no],
  [ASX_VAR_UNSET([ZLIB_C_INCLUDE])
   ASX_VAR_UNSET([ZLIB_C_LIB])])
dnl ----------------------------------------------------------------------
dnl   test for SZIP compression library
dnl
ACX_C_PACKAGE([szip],[szlib.h],[#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif],,
  [AC_MSG_NOTICE([szip library header not found, szip compression not supported])],[SZ_Decompress],[sz],,,dnl
  [AC_MSG_NOTICE([cannot link szip library, szip compression not supported])])
AS_IF([test $have_szip_c_bindings = no],
  [ASX_VAR_UNSET([SZIP_C_INCLUDE])
   ASX_VAR_UNSET([SZIP_C_LIB])])
dnl
dnl ----------------------------------------------------------------------
dnl   test for HDF5 file format library, netCDF 4.x usually requires it
dnl
AC_ARG_ENABLE([hdf5],
  [AC_HELP_STRING([--disable-hdf5],dnl
     [disable checks for HDF5 library @<:@default: auto@:>@])],,
  [enable_hdf5=auto])

AS_IF([test x"$enable_hdf5" != xno],
  [ACX_C_PACKAGE([hdf5],[hdf5.h],,ACX_M4_GENERATE_SUBSETS([[$SZIP_C_INCLUDE],[$ZLIB_C_INCLUDE],[$MPI_C_INCLUDE]],[ ]),dnl
     [AC_MSG_NOTICE([HDF5 library header not found, will not support netCDF 4.x])],[H5open],[hdf5],ACX_M4_GENERATE_SUBSETS([[$SZIP_C_LIB],[$ZLIB_C_LIB],[$LIBM]],[ ]),,dnl
     [AC_MSG_NOTICE([cannot link HDF5 library])])
   AS_IF([test $have_hdf5_c_bindings = yes],dnl
     [AC_DEFINE([HAVE_LIBHDF5],[1],[HDF 5 library is present if defined to 1])],dnl
     [ASX_VAR_UNSET([HDF5_C_INCLUDE])
      ASX_VAR_UNSET([HDF5_C_LIB])
      AS_IF([test x"$enable_hdf5" = xyes],
        [AC_MSG_FAILURE([Requested package HDF5 unavailable.])])
     ])
  ])
dnl
dnl ----------------------------------------------------------------------
dnl   test for netCDF file format library
dnl
ACX_FORTRAN_PACKAGE([NetCDF],[netcdf.inc],,,dnl
  [AC_MSG_FAILURE([Required include netcdf.inc not found or not compilable.])],dnl
[nf_open],[netcdf netcdff],ACX_M4_GENERATE_SUBSETS([[-lnetcdf],[-lhdf5_hl],[$HDF5_C_LIB]],[ ]),,dnl
  [AC_MSG_FAILURE([Cannot link with NetCDF library.])],[])


dnl szlib is required for some GRIB files:

dnl ######################################################################

AC_ARG_ENABLE([CDI],
  [AC_HELP_STRING([--disable-CDI],dnl
   [make use of CDI library for I/O @<:@default: enabled@:>@])])

enable_internal_cdi=no
AS_IF([test x"$enable_CDI" != x"no"],
  [CDI=""
   dnl ######################################################################
   dnl Prerequisites for File-I/O
   dnl ######################################################################

   dnl ----------------------------------------------------------------------
   dnl  test for JASPER compression library
   dnl
   AC_ARG_WITH([jasper],[AS_HELP_STRING([--without-jasper],dnl
      [disable support for JasPer JPEG2000 file compression @<:@default: auto@:>@])])
   AS_IF([test x"$with_jasper" = xno],
     [have_jasper_c_bindings=no],
     [ACX_C_PACKAGE([JasPer],[jasper/jasper.h],,,
        [AC_MSG_NOTICE([JasPer library header not found, JPEG2000 compression not supported])],[jas_stream_memopen],[jasper],,,dnl
        [AC_MSG_NOTICE([cannot link JasPer library, JPEG2000 compression not supported])])])
   AS_IF([test x$have_JasPer_c_bindings = xyes],
     [AC_DEFINE([HAVE_LIBJASPER],[1],[Use the JasPer JPEG2000 compression library.])
      AC_MSG_NOTICE([Compiling with JPEG2000 compression support via JasPer.])],
     [AC_MSG_NOTICE([Compiling without JPEG2000 compression support via JasPer.])
      ASX_VAR_UNSET([JASPER_C_INCLUDE])
      ASX_VAR_UNSET([JASPER_C_LIB])]
     )

   dnl ----------------------------------------------------------------------
   dnl  test wether to use internal or external cdi library
   ACX_FORTRAN_PACKAGE([CDI], [cdi.inc],,,dnl
     [ACX_FORTRAN_CHECK_INCLUDE_PATHS_IFELSE([cdi.inc],dnl
       [`cd $srcdir ; pwd`/support],dnl
       [enable_internal_cdi=yes],dnl
       [AC_MSG_FAILURE([Required include cdi.inc not found or not compilable.])],,,[internal])],dnl
[streamopenwrite],[cdi],dnl
ACX_M4_GENERATE_SUBSETS([[$ZLIB_C_LIB],[$SZIP_C_LIB],[$NETCDF_LIB],[$LIBM],[-lcpp]],[ ]),,dnl
    [AS_IF([test x$enable_internal_cdi = xno],
       [AC_MSG_NOTICE([Cannot link with external CDI library.])
        enable_internal_cdi=no])],[])],
  [CDI="${FPP_DEFOPT}NOCDI"
   CDI_LIB=""
   CDI_INCLUDE=""
   CDIROOT=""])

INTERNAL_CDI_CPPFLAGS=
AS_IF([test x$enable_internal_cdi = xyes],
  [AC_MSG_NOTICE([Using internal CDI library and headers.])
   dnl ----------------------------------------------------------------------
   dnl  use ZLIB compression library in internal CDI?
   dnl (this has already been probed for NetCDF)
   dnl
   AC_ARG_WITH([cdi-zlib],[AS_HELP_STRING([--without-cdi-zlib],dnl
     [disable support for zlib file compression @<:@default: auto@:>@])])
   AS_IF([test x"$with_cdi_zlib" != xno -a x$have_zlib_c_bindings = xyes],
     [AC_DEFINE([HAVE_LIBZ],[1],[Use the zlib compression library.])
      INTERNAL_CDI_CPPFLAGS="${INTERNAL_CDI_CPPFLAGS+$INTERNAL_CDI_CPPFLAGS }$ZLIB_C_INCLUDE"
      AC_MSG_NOTICE([Compiling with ZLIB compression support.])],
     [AC_MSG_NOTICE([Compiling without ZLIB compression support.])])
   dnl ----------------------------------------------------------------------
   dnl use for SZIP compression library for internal CDI?
   dnl (this has already been probed for NetCDF)
   dnl
   AC_ARG_WITH([cdi-szip],[AS_HELP_STRING([--without-cdi-szip],dnl
     [disable internal CDI support for szip file compression @<:@default: auto@:>@])])
   AS_IF([test x"$with_szip" != xno -a x$have_szip_c_bindings = xyes],
     [AC_DEFINE([HAVE_LIBSZ],[1],[Use the szip compression library.])
      INTERNAL_CDI_CPPFLAGS="${INTERNAL_CDI_CPPFLAGS+$INTERNAL_CDI_CPPFLAGS }$SZIP_C_INCLUDE"
      AC_MSG_NOTICE([Compiling with SZIP compression support.])],
     [AC_MSG_NOTICE([Compiling without SZIP compression support.])]
     )
   dnl ----------------------------------------------------------------------
   dnl build internal CDI internal GRIB I/O layer?
   dnl
   AC_ARG_ENABLE([cdi-grib],[AS_HELP_STRING([--disable-cdi-grib],dnl
      [disable internal CDI support for GRIB file format @<:@default: enabled@:>@])])
   AS_IF([test x"$enable_cdi_grib" = xno],
     [AC_MSG_NOTICE([Compiling without support for GRIB file format.])],
     [AC_MSG_NOTICE([Compiling internal CDI with support for GRIB file format.])
      AC_DEFINE([HAVE_LIBGRIB],[1],[Include GRIB file format support in internal CDI library.])])
   dnl ----------------------------------------------------------------------
   dnl build internal CDI internal CGRIBEX I/O layer?
   dnl
   AC_ARG_ENABLE([cdi-cgribex],[AS_HELP_STRING([--disable-cdi-cgribex],dnl
      [disable internal CDI support for CGRIBEX @<:@default: enabled@:>@])])
   AS_IF([test x"$enable_cdi_cgribex" = xno],
     [AC_MSG_NOTICE([Compiling without CGRIBEX support.])],
     [AC_MSG_NOTICE([Compiling internal CDI with CGRIBEX support.])
      AC_DEFINE([HAVE_LIBCGRIBEX],[1],[Include CGRIBEX support in internal CDI library.])])
   dnl ----------------------------------------------------------------------
   dnl build internal CDI internal EXTRA I/O layer?
   dnl
   AC_ARG_ENABLE([cdi-extra],[AS_HELP_STRING([--disable-cdi-extra],dnl
      [disable internal CDI support for EXTRA file format @<:@default: enabled@:>@])])
    AS_IF([test x"$enable_extra" = xno],
      [AC_MSG_NOTICE([Compiling without support for EXTRA file format.])],
      [AC_MSG_NOTICE([Compiling with support for EXTRA file format.])
       AC_DEFINE([HAVE_LIBEXTRA],[1],[Include EXTRA file format support in internal CDI library.])])
   dnl ----------------------------------------------------------------------
   dnl build internal CDI internal SERVICE I/O layer?
   dnl
   AC_ARG_ENABLE([cdi-service],[AS_HELP_STRING([--disable-cdi-service],dnl
     [disable support for SERVICE file format @<:@default: enabled@:>@])])
   AS_IF([test x"$enable_cdi_service" = xno],
     [AC_MSG_NOTICE([Compiling without support for SERVICE file format.])],
     [AC_MSG_NOTICE([Compiling with support for SERVICE file format.])
      AC_DEFINE([HAVE_LIBSERVICE],[1],[Include SERVICE file format support in internal CDI library.])])
   dnl ----------------------------------------------------------------------
   dnl build internal CDI internal IEG I/O layer?
   dnl
   AC_ARG_ENABLE([cdi-ieg],[AS_HELP_STRING([--disable-cdi-ieg],dnl
     [disable support for IEG file format @<:@default: enabled@:>@])])
   AS_IF([test x"$enable_ieg" = xno],
     [AC_MSG_NOTICE([Compiling without support for IEG file format.])],
     [AC_MSG_NOTICE([Compiling with support for IEG file format.])
      AC_DEFINE([HAVE_LIBIEG],[1],[Include IEG file format support in internal CDI library.])])
   dnl ----------------------------------------------------------------------
   dnl build internal CDI netCDF I/O layer?
   dnl
   AC_ARG_ENABLE([cdi_netcdf],[AS_HELP_STRING([--disable-cdi-netcdf],dnl
     [disable internal CDI support for netCDF file format @<:@default: enabled@:>@])])
   AS_IF([test x"$enable_cdi_netcdf" = xno],
     [AC_MSG_NOTICE([Compiling without support for netCDF file format.])],
     [AC_MSG_NOTICE([Compiling with support for netCDF file format.])
      AC_DEFINE([HAVE_LIBNETCDF],[1],[Include netCDF file format support in internal CDI library.])
      INTERNAL_CDI_CPPFLAGS="${INTERNAL_CDI_CPPFLAGS+$INTERNAL_CDI_CPPFLAGS }$NETCDF_INCLUDE"])
   dnl
   dnl ----------------------------------------------------------------------
   dnl enable diagnostics with mallinfo?
   AC_CHECK_HEADERS([malloc.h],dnl
     [AC_SEARCH_LIBS([mallinfo],[malloc],dnl
        [AC_DEFINE([HAVE_MALLINFO],[1],[Include mallinfo diagnostics])])])
   dnl
   dnl ----------------------------------------------------------------------
   dnl
   CDI=""
   CDI_LIB="$SZIP_C_LIB $ZLIB_C_LIB $JASPER_C_LIB $HDF5_C_LIB"
   CDI_INCLUDE="$FPP_INCOPT`cd $srcdir ; pwd`/support"],dnl
  [AS_IF([test x$enable_internal_cdi = xyes],dnl
     [AC_MSG_NOTICE([Using external CDI library and headers.])])])
AM_CONDITIONAL([USE_CDI_IO], [test x"$enable_CDI" != x"no"])
AM_CONDITIONAL([INTERNAL_CDI], [test x"$enable_internal_cdi" = xyes])
AC_SUBST([INTERNAL_CDI_CPPFLAGS])

AC_ARG_ENABLE([HAMOCC],
  [AC_HELP_STRING([--enable-HAMOCC],dnl
   [compile with ocean biogeochemistry model HAMOCC @<:@default: no@:>@])])

HAMOCC=""

AS_IF([test x"$enable_HAMOCC" = xyes],
  [HAMOCC="${FPP_DEFOPT}PBGC ${FPP_DEFOPT}PNETCDF ${FPP_DEFOPT}zzDIFFAT ${FPP_DEFOPT}zzFB_BGC_OCE"
   ENABLE_HAMOCC_RUN=true],
  [HAMOCC="${FPP_DEFOPT}zzPBGC ${FPP_DEFOPT}zzPNETCDF ${FPP_DEFOPT}zzDIFFAT ${FPP_DEFOPT}zzFB_BGC_OCE"
   ENABLE_HAMOCC_RUN=false])
AM_CONDITIONAL([USEPBGC],[test x"$enable_HAMOCC" = x"yes"])
dnl # Set default for run script preparation script.
dnl # Avoid setting this in Makefiles.
AC_SUBST(ENABLE_HAMOCC_RUN)
AM_SUBST_NOTMAKE(ENABLE_HAMOCC_RUN)

dnl ######################################################################

AC_ARG_ENABLE([low-mem],
  [AC_HELP_STRING([--enable-low-mem],dnl
   [save memory by reducing precomputation tables @<:@default: off@:>@])])
AS_CASE(["$enable_low_mem"],
  [yes|on], [FCFLAGS="$FCFLAGS ${FPP_DEFOPT}WETO_STORE_LOGICAL"],
  [off|no|''], [:],
  [*,*|bitvect], [flags="$enable_low_mem"
     while test -n "$flags" ; do
       flag=`echo "$flags" | sed 's/,.*//'`
       flags=`echo "$flags" | sed 's/@<:@^,@:>@*\(,\|$\)//'`
       AS_CASE(["$flag"],
         [bitvect], [FCFLAGS="$FCFLAGS ${FPP_DEFOPT}WETO_STORE_BITVECTOR"],
         [*], [AC_MSG_ERROR([Unrecognized option passed to --enable-low-mem])])
     done],
  [*], [AC_MSG_ERROR([Unrecognized option passed to --enable-low-mem])])

dnl ######################################################################

AC_ARG_ENABLE([setup-programs],
  [AC_HELP_STRING([--enable-setup-programs],dnl
   [build programs to compute initialization data @<:@default: off@:>@])],
  [enable_setup_programs=yes],[enable_setup_programs=no])
AM_CONDITIONAL([SETUPPROGS],[test x"$enable_setup_programs" = x"yes"])

dnl ######################################################################

AC_ARG_ENABLE([old-io],
  [AC_HELP_STRING([--enable-old-io],dnl
   [old I/O mechanism @<:@default: disabled@:>@])])
AS_IF([test x"$enable_old_io" = xyes],
  [OLD_IO_FCDEFS=${FPP_DEFOPT}OLD_IO],
  [OLD_IO_FCDEFS=])
AC_SUBST([OLD_IO_FCDEFS])

dnl ######################################################################

AC_ARG_ENABLE([new-io],
  [AC_HELP_STRING([--disable-new-io],dnl
   [new I/O mechanism @<:@default: enabled@:>@])])
AS_IF([test x"$enable_new_io" = xno],
  [NEW_IO_FCDEFS=${FPP_DEFOPT}NO_NEW_IO],
  [NEW_IO_FCDEFS=])
AC_SUBST([NEW_IO_FCDEFS])

dnl ######################################################################
dnl Let's see if the IBM Tivoli Workload Scheduler (TWS) LoadLeveler
dnl can be queried
dnl ######################################################################
ACX_TWS_LOADLEVELER
dnl Checks for header files.


dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
dnl
dnl Checks for library functions.
AC_CHECK_FUNCS([getrusage gettimeofday sysconf uname])
dnl   check method to trigger floating-point exceptions
AS_CASE([$host],
  [*-ibm-aix*], [
    AC_CHECK_HEADER([fptrap.h],[acx_aix_fptrap=true],[acx_aix_fptrap=false])
    AS_IF([test x$acx_aix_fptrap = xtrue],
      [AC_LANG_PUSH([C])
       for func in fp_any_enable fp_is_enabled fp_enable fp_disable \
           fp_enable_all fp_disable_all ; do
         AC_SEARCH_LIBS([$func],,,[acx_aix_fptrap=false
            break])
       done
       AC_LANG_POP([C])])
    ])
AM_CONDITIONAL([USE_AIX_FPTRAP],[test x$acx_aix_fptrap = xtrue])

dnl
dnl Check for byte ordering
AC_C_BIGENDIAN
dnl
dnl Check wether char is unsigned by default
AC_C_CHAR_UNSIGNED
AS_IF([test x$ac_cv_c_char_unsigned = xyes],
  [C_CHAR_IS_UNSIGNED=1],
  [C_CHAR_IS_UNSIGNED=0])
AC_SUBST([C_CHAR_IS_UNSIGNED])
dnl
dnl check for C datatype sizes
AC_CHECK_SIZEOF([char])
AC_CHECK_SIZEOF([short])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])
AC_CHECK_SIZEOF([float])
AC_CHECK_SIZEOF([double])
AC_CHECK_SIZEOF([long double])
AC_CHECK_SIZEOF([int *])
dnl
dnl type conversion depends on these types being same size
dnl or sizeof(long) > sizeof(int)
ACX_CHECK_SIZEOF_RELATION([long], [int])
ACX_CHECK_SIZEOF_RELATION([long long], [long])
dnl Check for Fortran types:
dnl first find wether we have iso_c_binding module
ACX_FORTRAN_CHECK_MODULE([ISO_C_BINDING], [integer(c_int) :: a])
dnl Defines FORT_INT1=1, FORT_INT2=1, FORT_INT4=1, FORT_INT8=1, FORT_INT16=1
dnl         FORT_REAL4=1, FORT_REAL8=1, FORT_REAL16=1
dnl         FORT_COMPLEX8=1, FORT_COMPLEX16=1, FORT_COMPLEX32=1
dnl         FORT_DOUBLE_COMPLEX=1 if the types are available.
ACX_FORTRAN_USUAL_TYPE_KINDS
dnl ACX_FORTRAN_C_INT
ACX_FORTRAN_C_INT64_T
ACX_FORTRAN_C_FLOAT
dnl For each defined type, find the corresponding C datatype.
dnl its a bit complicated, but necessary.
dnl If none can be found, use 0.

AS_IF([test -z "$FORT_INTEGER_LEN"],
  [ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE([integer])
   AS_IF([test -z "$acx_cv_fortran_sizeof_integer"],
     [echo "Unable to determine the length of a FORTRAN integer"
      echo "You can set the environment variable FORT_INTEGER_LEN"
      echo "with the number of bytes in a FORTRAN integer and rerun"
      echo "configure."],
     [fort_integer_len=$acx_cv_fortran_sizeof_integer])],
  [fort_integer_len=$FORT_INTEGER_LEN])

AS_IF([test -z "$FORT_REAL_LEN"],
  [ACX_FORTRAN_RUN_CHECK_SIZEOF([real])
   AS_IF([test -z "$acx_cv_fortran_sizeof_real"],
     [echo "Unable to determine the length of a FORTRAN real"
      echo "You can set the environment variable FORT_REAL_LEN"
      echo "with the number of bytes in a FORTRAN real and rerun"
      echo "configure."],
     [fort_real_len="$acx_cv_fortran_sizeof_real"])],
  [fort_real_len=$FORT_REAL_LEN])

AS_IF([test x$fort_integer_len != x4],
  [AC_MSG_NOTICE([INTEGER length should be 4 but is ${fort_integer_len}.])
   AC_MSG_FAILURE([Default type size assumption violated.], 1)
  ])

AC_MSG_CHECKING(for C type corresponding to Fortran INTEGER)
ac_lo=
AS_IF([test "$fort_integer_len" -gt 1],
  [for c_type in long_long long int short char ; do
     eval typelen=\$"ac_cv_sizeof_${c_type}"
     AS_IF([test $fort_integer_len = "$typelen"],
      [# Fortran type fort_integer_len corresponds to C type $ctype
       FORT_INT=${c_type}
       AC_MSG_RESULT([$c_type])
       ac_lo=1
       break])
   done])
AS_IF([test x$ac_lo = x],
  [eval FORT_INT=0
   AC_MSG_RESULT([not supported])])

AC_MSG_CHECKING([for C type corresponding to Fortran REAL])
ac_lo=
AS_IF([test "$fort_real_len" -gt 1],
  [for c_type in float double long_double  ; do
     eval typelen=\$"ac_cv_sizeof_${c_type}"
     AS_IF([test $fort_real_len = "$typelen"],
       [# Fortran type fort_real_len corresponds to C type $ctype
        FORT_REAL=${c_type}
        AC_MSG_RESULT([$c_type])
        ac_lo=1
        break])
   done])
AS_IF([test x$ac_lo = x],
  [eval FORT_REAL=0
   AC_MSG_RESULT([not supported])])
AC_DEFINE_UNQUOTED(FORT_INT, $FORT_INT,
  [Defined to C type of same size as Fortran INTEGER])
AC_DEFINE_UNQUOTED(FORT_REAL, $FORT_REAL,
  [Defined to C type of same size as Fortran REAL])
AC_DEFINE_UNQUOTED($FORTRANNAMES,,
  [Mnemonic describing type of name mangling used in Fortran objects])

dnl checking for library routines and intrinsics that might require
dnl emulation
ACX_FORTRAN_CHECK_LIBFUNC([get_environment_variable],[
  integer :: home_length
  character(len=10), parameter :: home_varname = 'HOME'
  call get_environment_variable(home_varname, length=home_length)
],[], [AC_LIBOBJ([get_environment_variable])])
AM_CONDITIONAL([GETENV_EMULATION],
  [test x$acx_fortran_routine_get_environment_variable = xno])


AR="${AR:-ar}"
ARFLAGS="${ARFLAGS:-cru}"
AC_SUBST([AR])
AC_SUBST([ARFLAGS])
AC_SUBST([ARCH])
# translate cpp DEFS into corresponding directives for the Fortran
# compiler to pass to the preprocessor
AC_CONFIG_COMMANDS_PRE([FCDEFS=`for ac_def in $DEFS ; do echo $ac_def ; done |
  sed 's/^-D/'"$FPP_DEFOPT"'/'`])
AC_SUBST([FCDEFS])

AC_SUBST([FCFLAGS])
AC_SUBST([FCLDFLAGS])
AC_SUBST([INCLUDES])
AC_SUBST([MPI])
AC_SUBST([CDI])
AC_SUBST([HAMOCC])
dnl
dnl Defaults for setup scripts
AC_ARG_WITH([pool-root],
  [AC_HELP_STRING([--with-pool-root=<DIR>],dnl
     [use model setup and initialization data in DIR @<:@default=/pool/data@:>@])],dnl
     [POOL_ROOT=$withval],dnl
     [POOL_ROOT=${POOL_ROOT-/pool/data}])
AC_SUBST([POOL_ROOT])
dnl
AC_CONFIG_FILES([Makefile])

AS_IF([test x = x"$PERL"],
  [AC_MSG_ERROR([PERL is undefined!])])

AC_CONFIG_FILES([util/f2html.pl util/move_to_lowercase_f.pl \
   util/move_to_uppercase_f.pl \
   contrib/aix/geometry.pl \
   run/prepare_sge_mpiom_hamocc_omip \
   run/prepare_sge_mpiom_hamocc \
   plots/2dplot/plot2d \
   setup/build-setup.sh \
   setup/mk_anta/mk_anta-etopo2.job \
   setup/mk_anta/tp_gridgen.job \
   setup/anta2nc/anta2nc.job \
   setup/mk_arcgri/mk_arcgri.job \
   setup/mk_phc/mk_phc.job \
   setup/forcing_omip/forcing_omip.job \
   setup/mk_topo/mk_topo.job \
   setup/mk_BEK/mk_BEK.job \
   setup/runoff_omip/runoff_omip.job \
  ], [chmod a+x $ac_file])
dnl
ACX_CHECK_DOXYGEN_LAYOUT
AC_SUBST_FILE(DOXYFILE_HEADER)
SVN_VERSION=`svnversion`
AC_SUBST(SVN_VERSION)
AC_CONFIG_FILES([Doxyfile])
AC_CONFIG_FILES([config/makef90depends], [chmod a+x config/makef90depends])
m4_define([f90makefiles],dnl
[bin/Makefile \
   support/Makefile \
   setup/mk_BEK/Makefile \
   setup/anta2nc/Makefile \
   setup/forcing_era40/Makefile \
   setup/forcing_ncep/Makefile \
   setup/forcing_omip/Makefile \
   setup/mk_anta/Makefile \
   setup/mk_arcgri/Makefile \
   setup/mk_phc/Makefile \
   setup/mk_topo/Makefile \
   setup/runoff_omip/Makefile
])dnl
AC_CONFIG_FILES(f90makefiles,dnl
  [sed '
/FCCOMPILE =/s/(FCFLAGS)/(FCFLAGS) \$(FCDEFS) \$(CPPFLAGS) \$(AM_FPPFLAGS)/
/FCLINK =/s/FCFLAGS/FCLDFLAGS/g
/\$(FCLINK)/s/(LIBS)/(LIBS) $(FCLIBS)/g
/	\$(FC) \$(AM_FCFLAGS)/s/(FCFLAGS)/(FCFLAGS) \$(FCDEFS) \$(FPPFLAGS) \$(AM_FPPFLAGS)/' \
   $ac_file  \
    >conftest.mk.fclinkrewrite;
   cmp $ac_file conftest.mk.fclinkrewrite >/dev/null \
    || mv conftest.mk.fclinkrewrite $ac_file; \
    rm -f conftest.mk.fclinkrewrite;
  ])
AC_CONFIG_FILES([\
   doc/Makefile \
   setup/Makefile \
   setup/setup-common.inc \
  ])
# force make to regenerate FC.deps file (because it depends on Makefile)
AC_CONFIG_COMMANDS([force_dep_regeneration],[sleep 1
touch f90makefiles])
ACX_USE_ALTERNATIVE_LIBM
AC_OUTPUT

echo ""
echo "Configuration done. Run $MAKE to build MPIOM."
echo ""
dnl remove a.out file that some tests leave lying around
/bin/rm -f a.out
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
