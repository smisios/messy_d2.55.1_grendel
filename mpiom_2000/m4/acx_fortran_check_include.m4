# _ACX_FORTRAN_CHECK_INCLUDE_IFELSE(HEADER-FILE,
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [INCLUDES])
# ------------------------------------------------------------------
# Check the compiler accepts HEADER-FILE. The INCLUDES might be defaulted.
m4_defun([_ACX_FORTRAN_CHECK_INCLUDE_IFELSE],
  [AC_LANG_PUSH([Fortran])
   AC_COMPILE_IFELSE([AC_LANG_SOURCE([program conftest
$4
@%:@include <$1>
end program conftest])],
     [AS_VAR_SET([acx_Include], [yes])])
   AC_LANG_POP([Fortran])
   AS_VAR_SET_IF([acx_Include],[$2],[$3])dnl
])# _AC_CHECK_HEADER_NEW
# ACX_FORTRAN_CHECK_INCLUDE_IFELSE(HEADER-FILE,
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [INCLUDES])
# ------------------------------------------------------------------
# Check the compiler accepts HEADER-FILE. The INCLUDES might be defaulted.
AC_DEFUN([ACX_FORTRAN_CHECK_INCLUDE_IFELSE],
  [AS_VAR_PUSHDEF([acx_Include], [acx_cv_fortran_include_$1])dnl
   AC_CACHE_CHECK([for $1],[acx_Include],dnl
     [_ACX_FORTRAN_CHECK_INCLUDE_IFELSE([$1],[$2],[$3],[$4])])
   AS_VAR_POPDEF([acx_Include])])
# ACX_FORTRAN_CHECK_INCLUDE_PATHS_IFELSE(HEADER-FILE, PATH...
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [INCLUDES], [INCLUDE-FLAGS], [TAG])
# ------------------------------------------------------------------
# Check the compiler accepts HEADER-FILE. The INCLUDES might be defaulted.
# TAG defaults to extra
AC_DEFUN([ACX_FORTRAN_CHECK_INCLUDE_PATHS_IFELSE],
  [AC_REQUIRE([AC_PROG_FPP])
   AS_VAR_PUSHDEF([acx_Include], [acx_cv_fortran_include_$1])dnl
   AC_MSG_CHECKING([for $1 m4_default([$7],[extra]) include path])
   AC_CACHE_VAL([acx_Include],dnl
     [ac_include_search_FCFLAGS_SAVE="$FCFLAGS"
      for ac_incdir in '' $2; do
        AS_IF([test -z "$ac_incdir"],
          [ac_res="none required"
           FCFLAGS="$6 $ac_include_search_FCFLAGS_SAVE"],
          [ac_res="$FPP_INCOPT$ac_incdir"
           FCFLAGS="$6 $ac_res $ac_include_search_FCFLAGS_SAVE"])
        _ACX_FORTRAN_CHECK_INCLUDE_IFELSE([$1],dnl
          [AS_IF([test -z "$ac_incdir"],dnl
             [AS_VAR_SET([acx_Include],["$6"])],dnl
             [AS_VAR_SET([acx_Include],["]m4_ifval([$6],[$6 ])[$FPP_INCOPT$ac_incdir"])])],,[$5])
        AS_VAR_SET_IF([acx_Include], [break])dnl
      done
      FCFLAGS="$ac_include_search_FCFLAGS_SAVE"])
   AS_VAR_SET_IF([acx_Include],dnl
     [AS_IF([test x"AS_VAR_GET([acx_Include])" = x],dnl
        [AC_MSG_RESULT([(none required)])],dnl
        [AC_MSG_RESULT([AS_VAR_GET([acx_Include])])])],
     [AC_MSG_RESULT([not found])])
   AS_VAR_SET_IF([acx_Include], [$3], [$4])[]dnl
   AS_VAR_POPDEF([acx_Include])dnl
])# ACX_FORTRAN_CHECK_INCLUDE_PATHS_IFELSE
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
