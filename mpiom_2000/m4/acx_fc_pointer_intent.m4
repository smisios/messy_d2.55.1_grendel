dnl test wether the fortran compiler supports specifying intent on pointer
dnl arguments
AC_DEFUN([ACX_FC_POINTER_INTENT],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_LANG_PUSH([Fortran])
   AC_MSG_CHECKING([wether $FC supports pointer arguments with intent attribute])
   AC_COMPILE_IFELSE([SUBROUTINE f(p)
     REAL, POINTER, INTENT(IN) :: p
END SUBROUTINE f
],dnl
     [FPPFLAGS="${FPPFLAGS} ${FPP_DEFOPT}HAVE_POINTER_WITH_ATTR_INTENT"
      AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])])
   AC_LANG_POP([Fortran])])