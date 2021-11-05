# ACX_SL_FC_HAVE_PERCENTLOC
# ---------------------
# Test whether the FC compiler has the %LOC extension.  If so, define
# the preprocessor variable HAVE_PERCENTLOC to be 1.
AC_DEFUN([ACX_SL_FC_HAVE_PERCENTLOC],
         [AC_REQUIRE([AC_PROG_FC])dnl
          AC_CACHE_CHECK([whether ${FC} has the %LOC extension],
                         [ac_cv_fc_have_percentloc],
                         [AC_LANG_PUSH([Fortran])
                          AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
      INTEGER I, ADDR
      I = 1
      ADDR = %LOC( I )
])],
                          ac_cv_fc_have_percentloc=yes,
                          ac_cv_fc_have_percentloc=no)
                          AC_LANG_POP([Fortran])])
          if test $ac_cv_fc_have_percentloc = yes; then
              AC_DEFINE([HAVE_PERCENTLOC], 1,
                        [Define to 1 if the Fortran compiler supports the VAX %LOC extension])
          fi
])# ACX_SL_FC_HAVE_PERCENTLOC
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:


