## -------------------------------------------- ##
##    Looking for Fortran Preprocessor.         ##
## -------------------------------------------- ##

# ACX_SL_FC_FIXEDFORM([SRCEXT], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# -------------------------------------------------------------------
# Look for compiler flags to make the Fortran (FC) compiler accept
# fixed-format source code, with a source extension of SRCEXT,
# and puts any necessary flags in FCFLAGS_fixed_<SRCEXT>.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile fixed-format code using new extension) and ACTION-IF-FAILURE
# (defaults to failing with an error message) if not.
#
# The known flags are:
#                -FI: Intel compiler (icc, ecc)
#            -qfixed: IBM compiler (xlf)
#             -fixed: NAG compiler
#           -Mnofree: PGI compiler
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
AC_DEFUN([ACX_SL_FC_FIXEDFORM],
  [AC_REQUIRE([AC_PROG_FC])
   AC_LANG_PUSH([Fortran])dnl
   AS_VAR_PUSHDEF([acx_sl_fc_fixed], [acx_sl_cv_fc_fixedform_$1])
   AC_CACHE_CHECK([for Fortran flag needed to allow fixed-form source for .$1 suffix],
     AS_VAR_PUSHDEF([], [ac_cv_fc_fixedform_$1]),
[ac_cv_fc_fixedform_$1=unknown
ac_ext=$1
ac_fc_fixedform_FCFLAGS_save=$FCFLAGS
for ac_flag in none -FI "-qfixed -qsuffix=$ac_ext" -fixed --fix -Mnofree
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_fixedform_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([
      PROGRAM FIXEDFORM
C THIS COMMENT SHOULD CONFUSE FREEFORM COMPILERS
      PRI  NT*, 'HELLO '//
     .      'WORLD.'
      ENDP ROGRAM
],
                    [ac_cv_fc_fixedform_$1=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_fixedform_FCFLAGS_save
])
if test "x$ac_cv_fc_fixedform_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Cannot compile fixed-form source with .$1 suffix], 77)])
else
  if test "x$ac_cv_fc_fixedform_$1" != xnone; then
    AC_SUBST(FCFLAGS_fixed_[]$1, "$ac_cv_fc_fixedform_$1")
  fi
  $2
fi
AC_LANG_POP([Fortran])dnl
])# ACX_SL_FC_FIXEDFORM

# AC_LANG_PP_MANDATORY_SOURCE
# ----------------------
# expands to program text that will compile only iff the source is
# preprocessed prior to compilation
AC_DEFUN([AC_LANG_PP_MANDATORY_SOURCE],
[_AC_LANG_DISPATCH([$0], _AC_LANG, $@)])

m4_define([AC_LANG_PP_MANDATORY_SOURCE(Fortran 77)],dnl
[AC_LANG_SOURCE([      PROGRAM FIXEDFORM
C THIS COMMENT SHOULD CONFUSE FREEFORM COMPILERS
      PRI  NT*, 'HELLO '//
     .      'WORLD.'
@%:@ifndef OK
      ENDP ROGRAM
@%:@endif
])])

m4_define([AC_LANG_PP_MANDATORY_SOURCE(Fortran)],
[AC_LANG_SOURCE([
program freeform
! FIXME: how to best confuse non-freeform compilers?
print *, 'Hello ', &
'world.'
@%:@ifndef OK
end program
@%:@endif
])])


# _ACX_SL_FPP_FIXEDFORM_F
# -------------------
# Related to ACX_SL_FPP_FIXEDFORM, but used only from _ACX_SL_PROG_FC_FPP.
# How do we directly compile a preprocessable .F file?
# This should be a no-op on all systems except those with case-sensitive
# filenames, and those which can't do direct compilation anyway.
# Do not put result into cached variable if it fails.
AC_DEFUN([_ACX_SL_FPP_FIXEDFORM_F],[
ac_ext=F
ac_fpp_fixedform_FCFLAGS_save=$FCFLAGS
for ac_flag in none "-x f77-cpp-input" "-FI -cpp" "-qfixed -qsuffix=cpp=F" \
    "-fixed -fpp" "-lfe \"-Cpp\" --fix"
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fpp_fixedform_FCFLAGS_save $ac_flag"
    AC_COMPILE_IFELSE([AC_LANG_PP_MANDATORY_SOURCE],
  [ac_cv_fpp_fixedform_F=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fpp_fixedform_FCFLAGS_save
if test "x$ac_cv_fpp_fixedform_F" = x; then
  AC_MSG_WARN([Cannot compile fixed-form preprocessable Fortran with a .F extension.])
else
  if test "$ac_cv_fpp_fixedform_F" != none; then
    AC_SUBST(FPPFLAGS_fixed_F, "$ac_cv_fpp_fixedform_F")
  fi
fi
])# _ACX_SL_FPP_FIXEDFORM_F


# ACX_SL_FPP_FIXEDFORM([SRCEXT], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# --------------------------------------------------------------------
# Look for compiler flags to make the Fortran (FC) compiler accept
# preprocessed fixed-format source code, with a source extension of
# SRCEXT, and puts any necessary flags in FPPFLAGS_fixed_<SRCEXT>.
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile fixed-format code using new extension) and ACTION-IF-FAILURE
# (defaults to failing with an error message) if not.
#
# Mostly, this is applicable only when using direct compilation.
# However the macro also sets FPP_PREPROCESS_EXT and FPP_COMPILE_EXT,
# based on SRCEXT.  SRCEXT can be either 'EXT' or 'EXT1:ext2'; in the
# first case, the preprocessor extension is 'EXT', and the compile
# extension 'ext' (ie, the preprocessor extension, lowercased); in the
# second, the preprocessor extension is 'EXT1' and the compile
# extension 'ext2'.
#
# The known flags are:
#              -x f77-cpp-input: g77
#                      -FI -cpp: Intel compiler (ifort)
# -qfixed -qsuffix=cpp=<SRCEXT>: IBM compiler (xlf)
#                   -fixed -fpp: NAG compiler
#             -lfe "-Cpp" --fix: Lahey compiler
#                      -Mnofree: PGI (no flag for preprocessing available)
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
# NB when updating this list of flags, also update those of the previous
# macro.
AC_DEFUN([ACX_SL_FPP_FIXEDFORM],
  [AC_REQUIRE([AC_PROG_FPP])
   AC_LANG_PUSH([Fortran])dnl
dnl Extract preprocessor extension _ac_ppext from $1, part preceding any ':'
   m4_define([_ac_ppext],  m4_bpatsubst([$1], [:.*]))dnl
   AC_CACHE_CHECK([for Fortran flag needed to allow preprocessed fixed-form source for ._ac_ppext suffix],
                ac_cv_fpp_fixedform_[]_ac_ppext,
[if test $ac_cv_fpp_build_rule = direct; then
  ac_cv_fpp_fixedform_[]_ac_ppext=unknown
  ac_ext=_ac_ppext
  ac_fpp_fixedform_FCFLAGS_save=$FCFLAGS
  for ac_flag in none "-x f77-cpp-input" "-FI -cpp" "-qfixed -qsuffix=cpp=_ac_ppext" "-fixed -fpp" "-lfe \"-Cpp\" --fix"
  do
    test "x$ac_flag" != xnone && FCFLAGS="$ac_fpp_fixedform_FCFLAGS_save $ac_flag"
    AC_COMPILE_IFELSE([AC_LANG_PP_MANDATORY_SOURCE],
                 [ac_cv_fpp_fixedform_[]_ac_ppext=$ac_flag; break])
  done
  rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
  FCFLAGS=$ac_fpp_fixedform_FCFLAGS_save
else
  ac_cv_fpp_fixedform_[]_ac_ppext=none
fi # test $ac_cv_fpp_build_rule = direct
])
if test "x$ac_cv_fpp_fixedform_[]_ac_ppext" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Cannot compile fixed-form source with ._ac_ppext suffix], 77)])
else
  if test "x$ac_cv_fpp_fixedform_[]_ac_ppext" != xnone; then
    AC_SUBST(FPPFLAGS_fixed_[]_ac_ppext, "$ac_cv_fpp_fixedform_[]_ac_ppext")
  fi
  $2
fi

FPP_PREPROCESS_EXT=_ac_ppext
FPP_COMPILE_EXT=m4_if(m4_index([$1], :), -1,
                      m4_tolower([$1]),
                      m4_bpatsubst([$1], [.*:]))

AC_LANG_POP([Fortran])dnl
])# ACX_SL_FPP_FIXEDFORM


# ACX_SL_FC_FREEFORM([SRCEXT], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# ------------------------------------------------------------------
# Look for compiler flags to make the Fortran (FC) compiler accept
# free-format source code, with a source extension of SRCEXT,
# and puts any necessary flags in FCFLAGS_free_<SRCEXT>.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile fixed-format code using new extension) and ACTION-IF-FAILURE
# (defaults to failing with an error message) if not.
#
# For backwards compatibility, this macro may be called without
# specifying SRCEXT, in which case, a default extension of f90
# is used.  This usage is deprecated.
#
# The known flags are:
#                        -ffree-form: GNU g77
#                                -FR: Intel compiler (icc, ecc)
#                              -free: Compaq compiler (fort), NAG compiler
#         -qfree -qsuffix=f=<SRCEXT>: IBM compiler (xlf) (generates a warning
#                                         with recent versions)
#     -qfree=f90 -qsuffix=f=<SRCEXT>: Newer xlf versions
#                             --nfix: Lahey compiler
#                 -Mfree, -Mfreeform: Portland Group compiler
#                          -freeform: SGI compiler
#                            -f free: Absoft Fortran
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
AC_DEFUN([ACX_SL_FC_FREEFORM],
[AC_REQUIRE([AC_PROG_FC])
AC_LANG_PUSH([Fortran])dnl
dnl default _AC_EXT to 'f90', if no argument is given.
m4_define([_AC_EXT], m4_if($1, [], f90, $1))dnl
AC_CACHE_CHECK([for Fortran flag needed to allow free-form source for .]_AC_EXT[ suffix],
                ac_cv_fc_freeform_[]_AC_EXT,
[ac_cv_fc_freeform_[]_AC_EXT=unknown
ac_ext=_AC_EXT
ac_fc_freeform_FCFLAGS_save=$FCFLAGS
for ac_flag in none -ffree-form -FR -free "-qfree=f90" "-qfree=f90 -qsuffix=f=$ac_ext"\
               -qfree "-qfree -qsuffix=f=$ac_ext" -Mfree -Mfreeform \
               -freeform "-f free" --nfix
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_freeform_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([
program freeform
! FIXME: how to best confuse non-freeform compilers?
print *, 'Hello ', &
'world.'
end program],
                    [ac_cv_fc_freeform_[]_AC_EXT=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_freeform_FCFLAGS_save
])
if test "x$ac_cv_fc_freeform_[]_AC_EXT" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Cannot compile free-form source with .]_AC_EXT[ suffix], 77)])
else
  if test "x$ac_cv_fc_freeform_[]_AC_EXT" != xnone; then
dnl  if the first argument was absent, then implement the old behaviour,
dnl  and simply append to variable FCFLAGS
    m4_if($1, [],
      [FCFLAGS="$FCFLAGS $ac_cv_fc_freeform_[]_AC_EXT"],
      [AC_SUBST(FCFLAGS_free_[]_AC_EXT, "$ac_cv_fc_freeform_[]_AC_EXT")])
  fi
  $2
fi
AC_LANG_POP([Fortran])dnl
])# ACX_SL_FC_FREEFORM


# ACX_SL_FPP_FREEFORM([SRCEXT], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# ------------------------------------------------------------------
# Look for compiler flags to make the Fortran (FC) compiler accept
# preprocessed free-format source code, with a source extension of SRCEXT,
# and puts any necessary flags in FPPFLAGS_free_<SRCEXT>.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile fixed-format code using new extension) and ACTION-IF-FAILURE
# (defaults to failing with an error message) if not.
#
# Mostly, this is applicable only when using direct compilation.
# However the macro also sets FPP_PREPROCESS_EXT and FPP_COMPILE_EXT,
# based on SRCEXT.  SRCEXT can be either 'EXT' or 'EXT1:ext2'; in the
# first case, the preprocessor extension is 'EXT', and the compile
# extension 'ext' (ie, the preprocessor extension, lowercased); in the
# second, the preprocessor extension is 'EXT1' and the compile
# extension 'ext2'.
#
# The known flags are:
#     -ffree-form -x f77-cpp-input: GNU g77
#                         -FR -cpp: Intel compiler (ifort)
#                       -free -cpp: Compaq compiler (fort), NAG compiler
#     -qfree -qsuffix=cpp=<SRCEXT>: IBM compiler (xlf) (generates a warning
#                                       with recent versions)
# -qfree=f90 -qsuffix=cpp=<SRCEXT>: Newer xlf versions
#                --nfix -lfe="Cpp": Lahey compiler
#               -Mfree, -Mfreeform: PGI (no flag for preprocessing available)
#                        -freeform: SGI compiler
#                          -f free: Absoft Fortran
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
AC_DEFUN([ACX_SL_FPP_FREEFORM],
[AC_REQUIRE([AC_PROG_FPP])
AC_LANG_PUSH([Fortran])dnl
dnl Extract preprocessor extension _ac_ppext from $1, part preceding any ':'
m4_define([_ac_ppext],  m4_bpatsubst([$1], [:.*]))dnl
AC_CACHE_CHECK([for Fortran flag needed to allow free-form source for ._ac_ppext suffix],
                ac_cv_fpp_freeform_[]_ac_ppext,
[if test $ac_cv_fpp_build_rule = direct; then
   ac_cv_fpp_freeform_[]_ac_ppext=unknown
   ac_ext=_ac_ppext
   ac_fpp_freeform_FCFLAGS_save=$FCFLAGS
   for ac_flag in none -ffree-form -FR -free "-qfree=f90" "-qfree=f90 -qsuffix=cpp=_ac_ppext"\
                  -qfree "-qfree -qsuffix=cpp=_ac_ppext" -Mfree -Mfreeform \
                  -freeform "-f free" --nfix
   do
      test "x$ac_flag" != xnone && FCFLAGS="$ac_fpp_freeform_FCFLAGS_save $ac_flag"
      AC_COMPILE_IFELSE([AC_LANG_PP_MANDATORY_SOURCE],
                        [ac_cv_fpp_freeform_[]_ac_ppext=$ac_flag; break])
   done
   rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
   FCFLAGS=$ac_fpp_freeform_FCFLAGS_save
else
   ac_cv_fpp_freeform_[]_ac_ppext=none
fi # test $ac_cv_fpp_build_rule = direct
])
if test "x$ac_cv_fpp_freeform_[]_ac_ppext" = xunknown; then
   m4_default([$3],
              [AC_MSG_ERROR([Cannot compile free-form source with ._ac_ppext suffix], 77)])
else
   if test "x$ac_cv_fpp_freeform_[]_ac_ppext" != xnone; then
      AC_SUBST(FPPFLAGS_free_[]_ac_ppext, "$ac_cv_fpp_freeform_[]_ac_ppext")
   fi
   $2
fi

FPP_PREPROCESS_EXT=_ac_ppext
FPP_COMPILE_EXT=m4_if(m4_index([$1], :), -1,
                      m4_tolower([$1]),
                      m4_bpatsubst([$1], [.*:]))

AC_LANG_POP([Fortran])dnl
])# ACX_SL_FPP_FREEFORM


# AC_LANG_TEST_PREPROC_COMPILE
# --------------------
# try to build a program in the current language that will only
# compile when preprocessed.
AC_DEFUN([AC_LANG_TEST_PREPROC_COMPILE],
[_AC_LANG_DISPATCH([$0], _AC_LANG, $@)])


dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
