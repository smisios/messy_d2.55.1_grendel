dnl ACX_CHECK_SIZEOF_RELATION(TYPE_A, TYPE_B)
dnl sets C_SIZEOF_TYPE_A_IS_GREATER_THAN_SIZEOF_TYPE_B to 1 or 0
dnl depending on actual relation
AC_DEFUN([ACX_CHECK_SIZEOF_RELATION],
  [AS_VAR_PUSHDEF([sizeof_type_a], [ac_cv_sizeof_$1])
   AS_VAR_PUSHDEF([sizeof_type_b], [ac_cv_sizeof_$2])
   AC_MSG_CHECKING([if sizeof($1) is greater than sizeof($2)])
   AS_IF([test AS_VAR_GET([sizeof_type_a]) -gt AS_VAR_GET([sizeof_type_b])],
     [AS_TR_CPP([C_$1_IS_LARGER_THAN_$2])=1
      AC_MSG_RESULT([yes])],
     [AS_TR_CPP([C_$1_IS_LARGER_THAN_$2])=0
      AC_MSG_RESULT([no])])
   AC_SUBST(AS_TR_CPP([C_$1_IS_LARGER_THAN_$2]))])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
