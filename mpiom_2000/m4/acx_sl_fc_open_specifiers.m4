# ACX_SL_FC_OPEN_SPECIFIERS(specifier ...)
# ------------------------------------
#
# The Fortran OPEN statement is a rich source of portability problems,
# since there are numerous common extensions which consiste of extra
# specifiers, several of which are useful when they are available.
# For each of the specifiers in the (whitespace-separated) argument
# list, define HAVE_FC_OPEN_mungedspecifier if the specifier may be
# given as  argument to the OPEN statement.  The `mungedspecifier' is the
# `specifier' converted to uppercase and with all characters outside
# [a-zA-Z0-9_] deleted.  Note that this may include `specifiers' such
# as "access='append'" and "[access='sequential',recl=1]" (note quoting of
# comma) to check combinations of specifiers.  You may not include a
# space in the `specifier', even quoted.  Each argument must be a
# maximum of 65 characters in length (to abide by Fortran 77
# line-length limits).
#
dnl Multiple m4_quote instances are necessary in case specifier includes comma.
dnl In the Fortran OPEN line, include status='scratch' unless status=???
dnl is in the specifier being tested.
dnl Put specifier on continuation line, in case it's long.
AC_DEFUN([ACX_SL_FC_OPEN_SPECIFIERS],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_LANG_PUSH([Fortran])
   AC_FOREACH([Specifier],
     m4_quote(m4_toupper([$1])),
     [m4_define([mungedspec],
      m4_bpatsubst(m4_quote(Specifier), [[^a-zA-Z0-9_]], []))
      AC_CACHE_CHECK([whether ${FC} supports OPEN specifier ]m4_quote(Specifier),
        [ac_cv_fc_spec_]mungedspec,
        [AC_COMPILE_IFELSE(AC_LANG_PROGRAM([],
           [      OPEN(UNIT=99,]m4_if(m4_bregexp(m4_quote(Specifier), [\<STATUS *=]), -1, [STATUS='SCRATCH'[,]], [])
     :m4_quote(Specifier)[)]),
          [ac_cv_fc_spec_]mungedspec=yes,
          [ac_cv_fc_spec_]mungedspec=no)])
        AS_IF([test $ac_cv_fc_spec_[]mungedspec = yes],
          AC_DEFINE([HAVE_FC_OPEN_]mungedspec, 1,
            [Define to 1 if the Fortran compiler supports OPEN specifier ]m4_quote(Specifier))
        )])
   AC_LANG_POP([Fortran])
])# ACX_SL_FC_OPEN_SPECIFIERS
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
