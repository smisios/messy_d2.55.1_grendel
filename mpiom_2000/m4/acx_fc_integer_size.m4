dnl
dnl _ACX_CHECK_SIZEOF_PROGRAM(Fortran)(type, size_to_test,
dnl                           [ACTION-IF-TRUE], [ACTION-IF-FALSE])
dnl
dnl  produce program that fails to compile if size_to_test equals
dnl  the size of type
dnl  
m4_define([_ACX_LANG_CHECK_SIZEOF_INTEGRAL_TYPE_PROGRAM(Fortran)],
[AC_LANG_PROGRAM([], [
   contains
     $3
     integer function callee(a)
       integer a(:)
       integer i
       $1 c
       do i=1,1024,(bit_size(c) - $2)
         a(i) = a(i)
       end do
       callee = a(1)
     end function callee
     subroutine caller
       integer b(1024)
       $1 c
       b(1) = callee(b(::(bit_size(c) - $2)))
     end subroutine
  ])
])
dnl
dnl ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE(type, [optinal-prologue])
dnl
dnl sets shell variable acx_cv_fortran_sizeof_mangled_type
dnl to size of type in bytes and defines cpp macro named like
dnl FORTRAN_SIZEOF_type
dnl where type is a mangled version of the fortran type declaration
dnl
AC_DEFUN([ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE],
[
  AC_REQUIRE([ACX_C_CHAR_BITS])
  AS_VAR_PUSHDEF([acx_fortran_Sizeof], [acx_cv_fortran_sizeof_$1])dnl
  AS_VAR_PUSHDEF([acx_fortran_Type], [acx_cv_fortran_type_$1])dnl
  ACX_FORTRAN_CHECK_TYPE([$1])
  AS_IF([test x]AS_VAR_GET([acx_fortran_Type])[ = xyes],[
     AC_CACHE_CHECK([size of Fortran type $1],[acx_fortran_Sizeof],[
       AC_LANG_PUSH([Fortran])
       ac_lo=$acx_cv_c_char_bits ac_hi=`expr 32 \* $acx_cv_c_char_bits`
       while :; do
       AC_COMPILE_IFELSE([_AC_LANG_DISPATCH(_ACX_LANG_CHECK_SIZEOF_INTEGRAL_TYPE_PROGRAM,
           _AC_LANG, [$1], $ac_lo, $2)],
 	 [ac_lo=`expr $ac_lo + $acx_cv_c_char_bits`
          AS_IF([test $ac_lo -gt $ac_hi],
            [ac_lo= ac_hi= ; break])],
         [ac_lo=`expr $ac_lo / $acx_cv_c_char_bits` ; break])
       done
       AC_LANG_POP([Fortran])
       AS_IF([test -z "$ac_lo"],
         [AC_MSG_FAILURE([cannot compute sizeof ($1)], 77)],
         [AS_VAR_SET([acx_fortran_Sizeof], [$ac_lo])])])
    AC_DEFINE_UNQUOTED(AS_TR_CPP(fortran_sizeof_$1),dnl
       AS_VAR_GET([acx_fortran_Sizeof]),dnl
       [The size of `$1', as computed by sizeof.])])
  AS_VAR_POPDEF([acx_fortran_Type])dnl
  AS_VAR_POPDEF([acx_fortran_Sizeof])dnl
])
dnl Local Variables:
dnl mode: autoconf
dnl End:
