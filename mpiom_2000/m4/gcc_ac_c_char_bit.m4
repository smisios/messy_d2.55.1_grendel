dnl Probe number of bits in a byte.
dnl Derived from gcc acinclude.m4 macro gcc_AC_C_CHAR_BIT.
dnl Note C89 requires CHAR_BIT >= 8.
dnl Defines CHAR_BIT if not found and sets shell variable acx_cv_c_char_bits
dnl to number of bits per character if determinable.
dnl
AC_DEFUN([ACX_C_CHAR_BITS],
[AC_CHECK_DECL([CHAR_BIT],dnl
[ac_cv_have_decl_CHAR_BIT=yes],[ac_cv_have_decl_CHAR_BIT=no],dnl
[[#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif]])
AC_CACHE_CHECK(number of bits in a byte, acx_cv_c_char_bits,
  [i=8
   acx_cv_c_char_bits=
   while test $i -lt 65; do
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,
       [switch(0) {
        case (unsigned char)((unsigned long)1 << $i)
             == ((unsigned long)1 << $i):
        case (unsigned char)((unsigned long)1<<($i-1))
             == ((unsigned long)1<<($i-1)):
        ; }])],
        [acx_cv_c_char_bits=$i; break])
     i=`expr $i + 1`
   done
   test -z "$acx_cv_c_char_bits" && acx_cv_c_char_bits=failed
  ])
AS_IF([test $acx_cv_c_char_bits = failed],
  [AC_MSG_ERROR(cannot determine number of bits in a byte)],
  [AS_IF([test $ac_cv_have_decl_CHAR_BIT = no],
      [AC_DEFINE_UNQUOTED(CHAR_BIT, $acx_cv_c_char_bits,
      [Define as the number of bits in a byte, if `limits.h' doesn't.])
  ])])])
dnl Local Variables:
dnl mode: autoconf
dnl End:
