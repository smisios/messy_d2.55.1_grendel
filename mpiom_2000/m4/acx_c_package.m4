dnl
dnl ACX_C_PACKAGE(PACKAGE,
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
AC_DEFUN([ACX_C_PACKAGE],
  [AC_LANG_PUSH([C])
   ACX_GENERIC_PACKAGE([$1],[$2],[-I],[$3],[$4],[$5],[$6],[-L],[$7],[$8],[$9],[$10],[$11],[$12])
   AC_LANG_POP([C])dnl
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
