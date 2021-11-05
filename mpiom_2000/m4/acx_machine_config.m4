dnl
dnl AC_GET_MH
dnl
dnl Load machine specific compiler options
dnl
AC_DEFUN([AC_GET_MH],
[AC_MSG_CHECKING([for machine dependent configuration])
changequote(,)
cat > confsed <<EOF
:redo
/\\\\\$/ {
    s/\\\\\$//g
    h
    s/.*//
    n
    s/^[ 	]*//
    H
    x
    s/\\n//
    bredo
}
/^[ 	]*[A-Za-z_][0-9A-Za-z_]*[	 ]*=/ {
    s/[ 	]*=[ 	]*/="/
    /="/s/$/"/
    p
    s/^\\([ 	]*\\)\\([A-Za-z_][0-9A-Za-z_]*\\)[	 ]*=.*/\\1export \\2/
}
/^ *#/{
    d
}
p

:next
EOF
changequote([,])
sed -n -f confsed $1 > conftest
. ./conftest
/bin/rm -f confsed conftest
if test "$1" != 0 ; then
    AC_MSG_RESULT($1)
else
    AC_MSG_RESULT(unavailable)
fi
])dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
