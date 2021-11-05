dnl
dnl ACX_FC_FLUSH_CALL determines wether Fortran flush can be invoked
dnl as
dnl a. intrinsic (Fortran 2003 mode),
dnl b. a subroutine call or
dnl c. a subroutine call after a 'USE f90_unix_io'-statement.
AC_DEFUN([ACX_FC_FLUSH_CALL],
  [AC_REQUIRE([AC_PROG_FPP])
   AC_CACHE_CHECK([for Fortran FLUSH variation],[acx_cv_fc_flush],
     [AC_LANG_PUSH([Fortran])
      AC_MSG_CHECKING([])
      while true ; do
        AC_LINK_IFELSE(dnl
[AC_LANG_PROGRAM([], [dnl
      OPEN(10,file='conftest.outflush')
      WRITE (10, *) 'conftest for flush'
      FLUSH(10)])],
          [acx_cv_fc_flush='intrinsic'
           break])
        AC_LINK_IFELSE(dnl
[AC_LANG_PROGRAM([], [dnl
      OPEN(10,file='conftest.outflush')
      WRITE (10, *) 'conftest for flush'
      CALL FLUSH(10)])],
          [acx_cv_fc_flush='call'
           break])
        AC_LINK_IFELSE(dnl
[AC_LANG_PROGRAM([], [dnl
      USE f90_unix_io, ONLY: flush
      OPEN(10,file='conftest.outflush')
      WRITE (10, *) 'conftest for flush'
      CALL FLUSH(10)])],
          [acx_cv_fc_flush='call after use f90_unix_io'
           break])
        acx_cv_fc_flush='unknown'
        break
      done])
   AS_CASE(["$acx_cv_fc_flush"],dnl
     [intrinsic],,
     [call],[FPPFLAGS="$FPPFLAGS ${FPP_DEFOPT}'FLUSH(fd)=CALL flush(fd)'"],
     ['call after use f90_unix_io'],dnl
     [FPPFLAGS="$FPPFLAGS ${FPP_DEFOPT}FLUSH_NEEDS_F90_UNIX_IO ${FPP_DEFOPT}'FLUSH=CALL flush'"])
   AC_LANG_POP([Fortran])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
