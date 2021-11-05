# _AS_TR_ARG_PREPARE
# -----------------
AC_DEFUN([_ASX_TR_ARG_PREPARE],
[AS_REQUIRE([_AS_CR_PREPARE])dnl
# Sed expression to map a string onto a valid argument string part.
asx_tr_arg="eval sed 'y%*+%pp%;s%[[^-$as_cr_alnum]]%-%g'"
])


# AS_TR_ARG(EXPRESSION)
# --------------------
# Transform EXPRESSION into a nice-looking argument suffix, i.e.
# for AC_ARG_WITH or AC_ARG_ENABLE.
# sh/m4 polymorphic.
AC_DEFUN([ASX_TR_ARG],
[AS_REQUIRE([_$0_PREPARE])dnl
AS_LITERAL_IF([$1],
	      [m4_bpatsubst(m4_translit([m4_translit([$1], [A-Z], [a-z])], [*+], [pp]),
			    [[^a-zA-Z0-9-]], [-])],
	      [`echo "$1" | $asx_tr_arg`])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
