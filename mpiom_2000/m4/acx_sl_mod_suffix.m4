# ACX_SL_FC_MOD_SUFFIX([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
# -----------------
# Determines the form of the filename of modules produced
# by the Fortran compiler.
# Tests for all forms of file extension I've (TOHW) found in the
# wild. Note that at least one compiler (PGI??) changes the
# case of the basename as well. Whether this happens is
# encoded in the variable ac_fc_mod_uppercase.
#
# This macro depends, of course, on the Fortran compiler producing
# module files. See comment to AC_FC_MOD_PATH_FLAG.
#
# FIXME: This will fail if an F77-only compiler is used.
# Currently we warn and continue. We should maybe error out.
#
AC_DEFUN([ACX_SL_FC_MOD_SUFFIX],
  [AC_MSG_CHECKING([for suffix of module files])
   AC_ARG_VAR([FCMODEXT], [file extension of compiled Fortran module files])
   ac_fc_mod_uppercase=no
   AC_LANG_PUSH([Fortran])
   AC_COMPILE_IFELSE([
      module conftest
       implicit none
       integer :: i
      end module conftest
     ])
   while :; do
     acx_fc_mod_name=
   m4_foreach([acx_fc_mod_name],dnl
     [[conftest.$FCMODEXT], [conftest.mod], [conftest.MOD], [conftest.M],
      [CONFTEST.MOD]],dnl
     [AS_IF([test -n "acx_fc_mod_name" -a -f "acx_fc_mod_name"],dnl
        [[acx_fc_mod_name]="acx_fc_mod_name" ; break])
     ])
     break
   done
   rm -f conftest*
   AC_LANG_POP([Fortran])
dnl
   AS_CASE(["$acx_fc_mod_name"],dnl
     [conftest.$FCMODEXT], [:],
     [CONFTEST.$FCMODEXT], [ac_fc_mod_uppercase=yes],
     [conftest.mod], [FCMODEXT=mod],
     [conftest.MOD], [FCMODEXT=MOD],
     [conftest.M], [FCMODEXT=M],
     [CONFTEST.MOD], [FCMODEXT=MOD
        ac_fc_mod_uppercase=yes])
   AC_MSG_RESULT([${FCMODEXT-not found}])
   AS_VAR_TEST_SET([FCMODEXT], [$1], [m4_ifval([$2],[$2],dnl
     [AC_MSG_WARN([Could not find Fortran module file extension.])])])
dnl
   AS_IF([test $ac_fc_mod_uppercase = yes],
     [FCMODCASE=uc
      AC_MSG_NOTICE([Fortran module filenames are uppercase.])],
     [FCMODCASE=lc])
   AC_SUBST([FCMODEXT])
   AC_SUBST([FCMODCASE])
])dnl _ACX_SL_FC_MOD_SUFFIX
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
