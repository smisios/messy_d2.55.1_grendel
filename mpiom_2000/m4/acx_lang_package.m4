dnl
dnl ACX_GENERIC_PACKAGE(PACKAGE, [INCLUDE], INC_FLAG, [EXTRA-INCLUDES],
dnl [EXTRA-INCLUDEFLAGS], [ACTION-IF_HEADER-NOT-FOUND], [FUNCTION],
dnl LIBFLAG, [LIB-CANDIDATES], [EXTRA-LIBS], [EXTRA-LIBFLAGS],
dnl [ACTION-IF-LIB-NOT-FOUND], [DEFAULT-ROOT])
dnl -------------------------------------------------------------------
dnl Check wether INCLUDE can be compiled and FUNCTION is found in
dnl LIB-CANDIDATES with current language compiler. Sets PACKAGE_LANG_LIB
dnl and PACKAGE_LANG_INCLUDE variables to necessary compiler
dnl switches and arguments such that inclusion/compilation succeed for
dnl program using INCLUDE or FUNCTION respectively.  Also defines
dnl configure --with arguments for PACKAGEROOT, PACKAGE-LIB and
dnl PACKAGE-INCLUDE.
dnl The code sets the shell variable have_PACKAGE__AC_LANG_ABBREV_bindings
dnl to yes if all requested tests succeeded, to no if any test failed.
dnl The library check will not run if the header check failed.
AC_DEFUN([ACX_GENERIC_PACKAGE],
  [AC_REQUIRE([_ASX_TR_ARG_PREPARE])
   AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])=yes
   AC_SUBST(AS_TR_CPP([$1][root]))
   AC_ARG_WITH(ASX_TR_ARG([$1-root]),
     [AS_HELP_STRING([--with-]ASX_TR_ARG([$1])[-root],
        [set directory to search for $1 headers and library]m4_ifval([$13],[, @<:@default=$13@:>@]))],
        [AS_TR_CPP([$1][root])="$AS_TR_SH([with_]ASX_TR_ARG([$1])[_root])"],
        m4_ifval([$13], [AS_TR_CPP([$1][root])=$13]))
   AS_VAR_SET_IF(AS_TR_CPP([$1][root]),
     [AS_VAR_SET_IF([AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])],,
        [AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])="$8$[]AS_TR_CPP([$1][root])/lib"])
      AS_VAR_SET_IF([AS_TR_CPP([$1_]_AC_LANG_ABBREV[_INCLUDE])],,dnl
        [AS_TR_CPP([$1_]_AC_LANG_ABBREV[_INCLUDE])="$3$[]AS_TR_CPP([$1][root])/include"])])
   m4_ifval([$2],
     [AC_ARG_WITH(ASX_TR_ARG([$1-include]),
        [AS_HELP_STRING([--with-[]ASX_TR_ARG([$1])[]-include],
           [specifically set directory to search for $1 headers, ]dnl
[@<:@default=$]AS_TR_SH(ASX_TR_ARG([with_$1_root]))[/include@:>@])],
        AS_TR_CPP([$1_]_AC_LANG_ABBREV[_INCLUDE])[="$3$AS_TR_SH(ASX_TR_ARG([with_$1_include]))"],
        [])
      AC_ARG_VAR(AS_TR_CPP([$1_]_AC_LANG_ABBREV[_INCLUDE]),dnl
[specifically set flags to use when compiling sources
using $1 includes.])
      ACX_LANG_CHECK_INCLUDE_PATHS_IFELSE([$2],[],
        [AS_TR_CPP([$1_]_AC_LANG_ABBREV[_INCLUDE])="$[]AS_TR_SH([acx_cv_]_AC_LANG_ABBREV[_include_]$2)"],dnl
        [AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])=no
         $6],dnl
        [$4],[$[]AS_TR_CPP([$1_]_AC_LANG_ABBREV[_INCLUDE])],m4_ifval([$5],[[$5]],[[[]]]))])
   m4_ifval([$7],
     [AC_ARG_WITH(ASX_TR_ARG([$1-lib]),
        [AS_HELP_STRING([--with-]ASX_TR_ARG([$1])[-lib],
        [specifically set directory to search for $1 library, ]dnl
[@<:@default=$]AS_TR_SH(ASX_TR_ARG([with_$1_root]))[/lib@:>@])],
        AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])[="$8$AS_TR_SH(ASX_TR_ARG([with_$1_lib]))"],
        [])
      AC_ARG_VAR(AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB]),dnl
[specifically set flags to use when linking $1.])
      AC_SUBST(AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB]))
      AS_IF([test x$]AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])[ = xyes],
        [AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_$7_]_AC_LANG_ABBREV)dnl
         ACX_OPTION_SEARCH_LIBS_MULTI([$7],[$9],,dnl
           [AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])=no
            $12],[$10],m4_ifval([$11],[$11 ])$[]AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB]))
         AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])="AS_VAR_GET([ac_Search])"
         AS_VAR_POPDEF([ac_Search])dnl
        ])
     ])dnl
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
