dnl ACX_FORTRAN_CHECK_MODULE(mod_name, user code)
dnl   defines HAVE_FORTRAN_MODULE_mod_name if module is present and
dnl   user code enclosed in program/end program compiles on use
AC_DEFUN([ACX_FORTRAN_CHECK_MODULE],
  [AS_VAR_PUSHDEF([have_fortran_module],
    [AS_TR_SH([acx_cv_fortran_have_module_$1])])dnl
   AC_CACHE_CHECK([Fortran module $1], [have_fortran_module],[
     AC_LANG_PUSH([Fortran])
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM(, [
       use $1
       $2
       ])],
       [AS_VAR_SET([have_fortran_module], [yes])],
       [AS_VAR_SET([have_fortran_module], [no])])
     AC_LANG_POP([Fortran])])
   AS_IF([test x"AS_VAR_GET([have_fortran_module])" = xyes],dnl
     [AC_DEFINE(AS_TR_CPP([HAVE_FORTRAN_MODULE_$1]), 1,dnl
        [Defined if Fortran module $1 is available])])
   AS_VAR_POPDEF([have_fortran_module])dnl
])
dnl Local Variables:
dnl mode: autoconf
dnl End:
