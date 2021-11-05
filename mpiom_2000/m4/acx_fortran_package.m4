dnl
dnl ACX_FORTRAN_PACKAGE(PACKAGE,
dnl    [INCLUDE], [EXTRA-INCLUDES], [EXTRA-INCLUDEFLAGS],
dnl       [ACTION-IF_HEADER-NOT-FOUND],
dnl    [FUNCTION], [LIB-CANDIDATES], [EXTRA-LIBS], [EXTRA-LIBFLAGS],
dnl      [ACTION-IF-LIB-NOT-FOUND],
dnl    [DEFAULT-ROOT])
dnl -------------------------------------------------------------------
dnl Check wether INCLUDE can be compiled and FUNCTION is found in
dnl LIB-CANDIDATES. Sets PACKAGE_LIB and PACKAGE_INCLUDE variables to
dnl necessary Fortran compiler switches and arguments such that
dnl inclusion/compilation succeed for program using INCLUDE or
dnl FUNCTION respectively.
dnl Also defines configure --with arguments for PACKAGEROOT,
dnl PACKAGE-LIB and PACKAGE-INCLUDE.
AC_DEFUN([ACX_FORTRAN_PACKAGE],
  [AC_LANG_PUSH([Fortran])
   AC_REQUIRE([AC_PROG_FPP])
   AC_REQUIRE([_ASX_TR_ARG_PREPARE])
   AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])=yes
   AC_SUBST(AS_TR_CPP([$1][root]))
   AC_ARG_WITH(ASX_TR_ARG([$1-root]),
     [AS_HELP_STRING([--with-]ASX_TR_ARG([$1])[-root],
        [set directory to search for $1 headers and library]m4_ifval([$11],[, @<:@default=$11@:>@]))],
        [AS_TR_CPP([$1][root])="$AS_TR_SH([with_]ASX_TR_ARG([$1])[_root])"],
        m4_ifval([$11], [AS_TR_CPP([$1][root])=$11]))
   AS_VAR_SET_IF(AS_TR_CPP([$1][root]),
     [AS_TR_CPP([$1_LIB])="-L$[]AS_TR_CPP([$1][root])/lib $[]AS_TR_CPP([$1_LIB])"
      AS_TR_CPP([$1_INCLUDE])="$FPP_INCOPT$[]AS_TR_CPP([$1][root])/include $[]AS_TR_CPP([$1_INCLUDE])"])
   m4_ifval([$2],
     [AC_ARG_WITH(ASX_TR_ARG([$1-include]),
        [AS_HELP_STRING([--with-[]ASX_TR_ARG([$1])[]-include],
           [specifically set directory to search for $1 headers, ]dnl
[@<:@default=$]AS_TR_SH(ASX_TR_ARG([with_$1_root]))[/include@:>@])],
        AS_TR_CPP([$1_INCLUDE])[="$FPP_INCOPT$AS_TR_SH(ASX_TR_ARG([with_$1_include]))"],
        [])
      AC_SUBST(AS_TR_CPP([$1_INCLUDE]))
      ACX_FORTRAN_CHECK_INCLUDE_PATHS_IFELSE([$2],[],
        [],dnl
        [AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])=no
         $5],dnl
        [$3],m4_ifval([$4],[$4 ])[$[]AS_TR_CPP([$1_INCLUDE])])])
   m4_ifval([$6],
     [AC_ARG_WITH(ASX_TR_ARG([$1-lib]),
        [AS_HELP_STRING([--with-]ASX_TR_ARG([$1])[-lib],
        [specifically set directory to search for $1 library, ]dnl
[@<:@default=$]AS_TR_SH(ASX_TR_ARG([with_$1_root]))[/lib@:>@])],
        AS_TR_CPP([$1_LIB])[="-L$AS_TR_SH(ASX_TR_ARG([with_$1_lib]))"],
        [])
      AC_SUBST(AS_TR_CPP([$1_LIB]))
      AS_IF([test x$]AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])[ = xyes],
        [AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_$6_]_AC_LANG_ABBREV)dnl
         ACX_OPTION_SEARCH_LIBS_MULTI([$6],[$7],,dnl
           [$10],[$8],m4_ifval([$9],[$9 ])$[]AS_TR_CPP([$1_LIB]))
         AS_TR_CPP([$1_LIB])="AS_VAR_GET([ac_Search])"
         AS_VAR_POPDEF([ac_Search])])
     ])dnl
   AC_LANG_POP([Fortran])dnl
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
