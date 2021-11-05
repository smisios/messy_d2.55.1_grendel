# ACX_FORTRAN_CHECK_TYPE(TYPE,
#		     [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#		     [USEINCLUDE])
# ------------------------------------------------------------
# Check whether the type TYPE is supported by the system, maybe via
# the provided includes and module use statements.
#
# Analogous to C specific macro AC_CHECK_TYPE in the newer syntax version
AC_DEFUN([ACX_FORTRAN_CHECK_TYPE],
[
  AS_VAR_PUSHDEF([acx_fortran_Type], [acx_cv_fortran_type_$1])dnl
  AC_CACHE_CHECK([for Fortran type $1], [acx_fortran_Type],
    [AC_LANG_PUSH([Fortran])
     AC_COMPILE_IFELSE([
         AC_LANG_PROGRAM(,[$4
       $1 a
       a = max(a,  a)])
       ],
       [AS_VAR_SET([acx_fortran_Type], [yes])],
       [AS_VAR_SET([acx_fortran_Type], [no])]
     )
     AC_LANG_POP([Fortran])])
  m4_ifval([$2$3],dnl
    [AS_IF([test x]AS_VAR_GET([acx_fortran_Type])[ = xyes], [$2], [$3])])dnl
  AS_VAR_POPDEF([acx_fortran_Type])dnl
])
# ACX_FORTRAN_CHECK_TYPES(TYPES,
#		         [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#		         [INCLUDES = DEFAULT-INCLUDES])
# --------------------------------------------------------
# TYPES is an m4 list.  There are no ambiguities here, we mean the newer
# AC_CHECK_TYPE.
AC_DEFUN([ACX_FORTRAN_CHECK_TYPES],
[m4_foreach([AC_Type], [$1],
  [ACX_FORTRAN_CHECK_TYPE(AC_Type,
    [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_FORTRAN_[]AC_Type), 1,
			[Define to 1 if the system has the type `]AC_Type['.])
$2],
    [$3],
    [$4])])])
dnl
dnl _ACX_CHECK_TYPE_DEFINE_HELPER
dnl
m4_define([_ACX_CHECK_TYPE_DEFINE_HELPER],
  [AS_VAR_PUSHDEF([shellVar], [$2])dnl
   ACX_FORTRAN_CHECK_TYPES([$1],
     [AS_VAR_SET(shellVar, 1)])
   AS_VAR_POPDEF([shellVar])])
dnl ACX_FORTRAN_TYPE_KINDS tests to see if the following fortran datatypes are
dnl supported: INTEGER(1), INTEGER(2), INTEGER(4), INTEGER(8),
dnl            REAL(4), REAL(8), REAL(16),
dnl            DOUBLE_COMPLEX, COMPLEX(4), COMPLEX(8), COMPLEX(16)
dnl
AC_DEFUN([ACX_FORTRAN_USUAL_TYPE_KINDS],dnl
[
_ACX_CHECK_TYPE_DEFINE_HELPER([integer(kind=1)], [FORT_INT1])
_ACX_CHECK_TYPE_DEFINE_HELPER([integer(kind=2)], [FORT_INT2])
_ACX_CHECK_TYPE_DEFINE_HELPER([integer(kind=4)], [FORT_INT4])
_ACX_CHECK_TYPE_DEFINE_HELPER([integer(kind=8)], [FORT_INT8])
dnl
_ACX_CHECK_TYPE_DEFINE_HELPER([real(kind=4)], [FORT_REAL4])
_ACX_CHECK_TYPE_DEFINE_HELPER([real(kind=8)], [FORT_REAL8])
_ACX_CHECK_TYPE_DEFINE_HELPER([real(kind=16)], [FORT_REAL16])
dnl
_ACX_CHECK_TYPE_DEFINE_HELPER([double complex], [DOUBLE_COMPLEX])
_ACX_CHECK_TYPE_DEFINE_HELPER([complex(kind=4)], [FORT_COMPLEX8])
_ACX_CHECK_TYPE_DEFINE_HELPER([complex(kind=8)], [FORT_COMPLEX16])
_ACX_CHECK_TYPE_DEFINE_HELPER([complex(kind=16)], [FORT_COMPLEX32])
])dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
