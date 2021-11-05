# ACX_SL_FC_LITERAL_BACKSLASH
# -----------------------
#
# Check whether the compiler regards the backslash character as an escape
# character.  The Standard doesn't say anything about this, but many Unix
# Fortran compilers interpret '\n', for example, as a newline, and '\\' as
# a single backslash.
#
# Test the behaviour of the currently selected compiler, and define
# FC_LITERAL_BACKSLASH to 1 if backslashes are treated literally -- that is
# if '\\' is interpreted as a _pair_ of backslashes and thus that '\n' is
# interpreted as a pair of characters.
AC_DEFUN([ACX_SL_FC_LITERAL_BACKSLASH],
   [AC_REQUIRE([AC_PROG_FC])dnl
    AC_CACHE_CHECK([whether ${FC} interprets backslashes literally],
        ac_cv_fc_literal_backslash,
       [AC_LANG_PUSH([Fortran])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
dnl Four backslashes here -- this is expanded by the shell in writing
dnl the text to the file.  We want to end up with TEST\\TEST in the source.
[
      write(*,'("TEST\\\\TEST")')
])],
## 'strings' is portable, yes?
           [if strings conftest.$ac_objext | grep 'TEST\\\\TEST' >/dev/null; then
                ac_cv_fc_literal_backslash=yes
            else
                ac_cv_fc_literal_backslash=no
            fi],
           [AC_MSG_WARN([cannot compile a program with backslashes!])
            ac_cv_fc_literal_backslash=unknown])
        AC_LANG_POP([Fortran])])
    if test $ac_cv_fc_literal_backslash = yes; then
        AC_DEFINE([FC_LITERAL_BACKSLASH], 1,
                  [Define to 1 if the Fortran compiler interprets '\\' as a pair of characters, not one])
    fi
])# ACX_SL_FC_LITERAL_BACKSLASH
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:

