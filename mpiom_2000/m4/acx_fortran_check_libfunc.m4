dnl ACX_FORTRAN_CHECK_LIBFUNC(SUBROUTINE, ROUTINE-CALL,
dnl   ACTION-IF-TRUE, ACTION-IF-FALSE, [OPTIONAL-PROLOGUE])
dnl   Defines HAVE_FORTRAN_ROUTINE_SUBROUTINE if SUBROUTINE is present
dnl   and a fortran program with ROUTINE-CALL compiles and links
dnl   (after the optional PROLOGUE).
dnl   Also runs shell commands ACTION-IF-* accordingly.
AC_DEFUN([ACX_FORTRAN_CHECK_LIBFUNC],dnl
  [AS_VAR_PUSHDEF([have_fortran_routine],
    [AS_TR_SH([acx_cv_fortran_have_routine_$1])])dnl
   AC_CACHE_CHECK([Fortran subroutine $1], [have_fortran_routine],[
     AC_LANG_PUSH([Fortran])
     AC_LINK_IFELSE([AC_LANG_PROGRAM(, [
       $5
       $2
       ])],
       [AS_VAR_SET([have_fortran_routine], [yes])],
       [AS_VAR_SET([have_fortran_routine], [no])])
     AC_LANG_POP([Fortran])])
   AS_IF([test x"AS_VAR_GET([have_fortran_routine])" = xyes],dnl
     [AC_DEFINE(AS_TR_CPP([HAVE_FORTRAN_ROUTINE_$1]), 1,dnl
        [Defined if Fortran routine $1 is available])
        $3],
     [$4])
   AS_VAR_POPDEF([have_fortran_routine])dnl
  ])
dnl Local Variables:
dnl mode: autoconf
dnl End:
