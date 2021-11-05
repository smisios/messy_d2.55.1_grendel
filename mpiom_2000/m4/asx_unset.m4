# ASX_VAR_UNSET(VARIABLE)
# ---------------------------
# Unset shell VARIABLE.
# If the variable contains indirections (e.g. `ac_cv_func_$ac_func')
# perform whenever possible at m4 level, otherwise sh level.
AC_DEFUN([ASX_VAR_UNSET],
[AS_LITERAL_IF([$1],
	       [unset $1],
	       [unset `eval "$1"`])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
