# ACX_ASSERT_LANG_IS_FORTRAN_VARIANT
# Current language must be Fortran or Fortran 77
AC_DEFUN([ACX_ASSERT_LANG_IS_FORTRAN_VARIANT],
  [m4_case(_AC_LANG, [Fortran], [],
     [Fortran 77], [],
     [m4_fatal([$0: current language is not Fortran: ] _AC_LANG)])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
