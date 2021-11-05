# ACX_LANG_CHECK_INCLUDE_PATHS_IFELSE(C)(HEADER-FILE, PATH...
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [INCLUDES], [INCLUDE-FLAGS])
# ------------------------------------------------------------------
# Check the compiler accepts HEADER-FILE. The INCLUDES might be defaulted.
AC_DEFUN([ACX_LANG_CHECK_INCLUDE_PATHS_IFELSE(C)],dnl
  [ACX_GENERIC_CHECK_INCLUDE_PATHS_IFELSE(CPPFLAGS,-I,$@)])
#
m4_define([ACX_LANG_INCLUDE_PROGRAM(C)],dnl
  [AC_LANG_SOURCE([AC_INCLUDES_DEFAULT([$2])
@%:@include <$1>])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
