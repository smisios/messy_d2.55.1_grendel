dnl ACX_LANG_OTHER_SUFFIX_CONFTEST(SUFFIX, SOURCE)
dnl produce source in same language environment as AC_LANG_CONFTEST uses
dnl but change ac_ext to SUFFIX, so that a file conftest.SUFFIX will be
dnl written instead (useful e.g. to produce a header file).
AC_DEFUN([ACX_LANG_OTHER_SUFFIX_CONFTEST],
  [m4_define([_ACX_LANG_OLD],_AC_LANG)
   AC_LANG_PUSH(_AC_LANG)
   ac_ext="$1"
   AC_LANG_CONFTEST([$2])
   AC_LANG_POP(_ACX_LANG_OLD)
   m4_undefine([_ACX_LANG_OLD])
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
