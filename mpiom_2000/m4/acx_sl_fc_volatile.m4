# ACX_SL_FC_HAVE_VOLATILE
# -------------------
# Test whether the FC compiler supports the VOLATILE statement. VOLATILE
# is used to stop the optimisation of a variable, so that it can be modified
# outside of the program itself. If supported set HAVE_VOLATILE to be 1.
AC_DEFUN([ACX_SL_FC_HAVE_VOLATILE],
         [AC_REQUIRE([AC_PROG_FC])dnl
          AC_CACHE_CHECK([whether ${FC} supports VOLATILE],
                         [ac_cv_fc_have_volatile],
                         [AC_LANG_PUSH([Fortran])
                          AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
      INTEGER I
      VOLATILE I
])],
                          ac_cv_fc_have_volatile=yes,
                          ac_cv_fc_have_volatile=no)
                          AC_LANG_POP([Fortran])])
          if test $ac_cv_fc_have_volatile = yes; then
              AC_DEFINE([HAVE_VOLATILE], 1,
                        [Define to 1 if the Fortran compiler supports VOLATILE])
          fi
])# ACX_SL_FC_HAVE_VOLATILE
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
