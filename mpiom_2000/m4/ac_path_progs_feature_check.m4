# _AC_PATH_PROGS_FEATURE_CHECK(VARIABLE, PROGNAME-LIST, FEATURE-TEST,
#                              [ACTION-IF-NOT-FOUND], [PATH=$PATH])
# -------------------------------------------------------------------
# FEATURE-TEST is called repeatedly with $ac_path_VARIABLE set to the
# name of a program in PROGNAME-LIST found in PATH.  FEATURE-TEST must set
# $ac_cv_path_VARIABLE to the path of an acceptable program, or else
# ACTION-IF-NOT-FOUND is executed; the default action (for internal use
# only) issues a fatal error message.  If a suitable $ac_path_VARIABLE is
# found in the FEATURE-TEST macro, it can set $ac_path_VARIABLE_found=':'
# to accept that value without any further checks.
m4_define([_AC_PATH_PROGS_FEATURE_CHECK],
[if test -z "$$1"; then
  ac_path_$1_found=false
  # Loop through the user's path and test for each of PROGNAME-LIST
  _AS_PATH_WALK([$5],
  [for ac_prog in $2; do
    for ac_exec_ext in '' $ac_executable_extensions; do
      ac_path_$1="$as_dir/$ac_prog$ac_exec_ext"
      AS_EXECUTABLE_P(["$ac_path_$1"]) || continue
$3
      $ac_path_$1_found && break 3
    done
  done])dnl
  if test -z "$ac_cv_path_$1"; then
    m4_default([$4],
      [AC_MSG_ERROR([no acceptable m4_bpatsubst([$2], [ .*]) could be dnl
found in m4_default([$5], [\$PATH])])])
  fi
else
  ac_cv_path_$1=$$1
fi
])


# AC_PATH_PROGS_FEATURE_CHECK(VARIABLE, PROGNAME-LIST,
#                             FEATURE-TEST, [ACTION-IF-NOT-FOUND=:],
#                             [PATH=$PATH])
# ----------------------------------------------------------------
# Designed to be used inside AC_CACHE_VAL.  It is recommended,
# but not required, that the user also use AC_ARG_VAR([VARIABLE]).
# If VARIABLE is not empty, set the cache variable
# $ac_cv_path_VARIABLE to VARIABLE without any further tests.
# Otherwise, call FEATURE_TEST repeatedly with $ac_path_VARIABLE
# set to the name of a program in PROGNAME-LIST found in PATH.  If
# no invocation of FEATURE-TEST sets $ac_cv_path_VARIABLE to the
# path of an acceptable program, ACTION-IF-NOT-FOUND is executed.
# FEATURE-TEST is invoked even when $ac_cv_path_VARIABLE is set,
# in case a better candidate occurs later in PATH; to accept the
# current setting and bypass further checks, FEATURE-TEST can set
# $ac_path_VARIABLE_found=':'.  Note that, unlike AC_CHECK_PROGS,
# this macro does not have any side effect on the current value
# of VARIABLE.
m4_define([AC_PATH_PROGS_FEATURE_CHECK],
[_$0([$1], [$2], [$3], m4_default([$4], [:]), [$5])dnl
])
