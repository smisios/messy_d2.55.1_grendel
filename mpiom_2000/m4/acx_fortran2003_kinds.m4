dnl
dnl
dnl ACX_FORTRAN_TYPE_KIND(FORTRAN_TYPE, FORTRAN_KIND, C_CORRESPONDENCE,
dnl                       [PROLOGUE], [C_PROLOGUE])
dnl
dnl First tries wether a fortran declaration
dnl FORTRAN_TYPE(kind=FORTRAN_KIND) compiles and if not tries to find
dnl a fortran kind constant that corresponds to C type
dnl C_CORRESPONDENCE.
dnl PROLOGUE could be something like 'use iso_c_binding'
dnl
AC_DEFUN([ACX_FORTRAN_TYPE_KIND],
 [ACX_FORTRAN_CHECK_TYPE([$1(kind=$2)],
    [],
    [m4_ifval([$5],[AC_CHECK_SIZEOF([$3],,[$5])],[AC_CHECK_SIZEOF([$3])])
     AS_VAR_PUSHDEF([acx_c_size], [ac_cv_sizeof_$3])dnl
     m4_case([$1],[integer],dnl
       [ac_kind_lo=1 ac_kind_hi=32
        while :; do
          ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE([$1(kind=$ac_kind_lo)])
          AS_VAR_PUSHDEF([acx_fortran_ikind_Type],dnl
            [acx_cv_fortran_type_$1(kind=$ac_kind_lo)])dnl
          AS_VAR_PUSHDEF([acx_fortran_Sizeof],dnl
            [acx_cv_fortran_sizeof_$1(kind=$ac_kind_lo)])dnl
          AS_IF([test -n AS_VAR_GET([acx_fortran_ikind_Type]) -a "AS_VAR_GET([acx_fortran_Sizeof])" = "AS_VAR_GET([acx_c_size])"],
            [acx_fortran_kind_subst="$ac_kind_lo" ac_kind_lo=
             ac_kind_hi= ; break],
            [ac_kind_lo=`expr $ac_kind_lo + 1`
             AS_IF([test $ac_kind_lo -gt $ac_kind_hi],
               [ac_kind_lo= ac_kind_hi= ; break])])
          AS_VAR_POPDEF([acx_fortran_Sizeof])
          AS_VAR_POPDEF([acx_fortran_ikind_Type])
        done],[real],dnl
       [ac_kind_lo=1 ac_kind_hi=32
        while :; do
          ACX_FORTRAN_RUN_CHECK_SIZEOF([$1(kind=$ac_kind_lo)])
          AS_VAR_PUSHDEF([acx_fortran_rkind_Type],dnl
            [acx_cv_fortran_type_$1(kind=$ac_kind_lo)])dnl
          AS_VAR_PUSHDEF([acx_fortran_Sizeof],dnl
            [acx_cv_fortran_sizeof_$1(kind=$ac_kind_lo)])dnl
          AS_IF([test -n AS_VAR_GET([acx_fortran_rkind_Type]) -a "AS_VAR_GET([acx_fortran_Sizeof])" = "AS_VAR_GET([acx_c_size])"],
            [acx_fortran_kind_subst="$ac_kind_lo" ac_kind_lo=
             ac_kind_hi= ; break],
            [ac_kind_lo=`expr $ac_kind_lo + 1`
             AS_IF([test $ac_kind_lo -gt $ac_kind_hi],
               [ac_kind_lo= ac_kind_hi= ; break])])
          AS_VAR_POPDEF([acx_fortran_Sizeof])
          AS_VAR_POPDEF([acx_fortran_rkind_Type])
        done],
       [AC_MSG_WARN([Cannot derive C type correspondence for Fortran type $1(kind=$2).])])
     AS_VAR_POPDEF([acx_c_size])
   ], [$4])
  AS_IF([test x$acx_fortran_kind_subst != x],
    [AC_MSG_NOTICE([Substituting $1 kind $acx_fortran_kind_subst for $3])
     AC_DEFINE_UNQUOTED(AS_TR_CPP($2), $acx_fortran_kind_subst, [type kind override])
     AC_DEFINE_UNQUOTED(HAVE_FORTRAN_ISO_[]AS_TR_CPP($2), $acx_fortran_kind_subst, [type kind override])])
  ASX_VAR_UNSET([acx_fortran_kind_subst])
])
dnl
AC_DEFUN([ACX_FORTRAN_C_INT],dnl
  [ACX_FORTRAN_TYPE_KIND([integer], [c_int], [int],[use iso_c_binding])])
AC_DEFUN([ACX_FORTRAN_C_INT64_T],dnl
  [AC_REQUIRE([AC_TYPE_INT64_T])
   ACX_FORTRAN_TYPE_KIND([integer], [c_int64_t], [int64_t],dnl
  [use iso_c_binding])])
AC_DEFUN([ACX_FORTRAN_C_FLOAT],dnl
  [ACX_FORTRAN_TYPE_KIND([real], [c_float], [float],[use iso_c_binding])]))
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
