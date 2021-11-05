dnl ACX_ALTERNATIVE_LIBM([CANDIDATES])
dnl tries to find a suitable -lm substitute
dnl defines an output variable LIBM to use in tests potentially adding -lm
dnl otherwise
AC_DEFUN([ACX_ALTERNATIVE_LIBM],
  [AC_REQUIRE([AX_WITH_PERL])
   AC_LANG_PUSH([C])
   AC_MSG_CHECKING([for alternative math library])
   AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_acosh_]_AC_LANG_ABBREV)dnl
   while :; do
     m4_foreach([libmSubst], [$1],
       [# try to remove -lm from LIBS first!
        LIBS_save=$LIBS
        LIBS=`echo "$LIBS" | $PERL -e 'while(<>) { s{\s\K-lm(?=\s|$)}{'"libmSubst"'}g ; print }'`
        _ACX_OPTION_SEARCH_LIBS([acosh],,[libmSubst])
        LIBS=$LIBS_save
        AS_VAR_SET_IF([ac_Search],[break])
       ])
     break
   done
   AS_VAR_SET_IF([ac_Search],
     [LIBM=AS_VAR_GET([ac_Search])
      LIBS=`echo "$LIBS" |  $PERL -e 'while(<>) { s{\s\K-lm(?=\s|$)}{'"$LIBM"'}g ; print }'`
      AC_MSG_RESULT([$LIBM])],
     [LIBM=-lm
      AC_MSG_RESULT([default -lm])])
   AS_VAR_POPDEF([ac_Search])dnl        _ACX_OPTION_SEARCH_LIBS([acosh],,
   AC_SUBST([LIBM])
   AC_LANG_POP([C])
  ])

dnl USE_ALTERNATIVE_LIBM([LIBVARS])
dnl substitute a libm alternative in LIBS FCLIBS output variables and
dnl in each element of [LIBVARS]
AC_DEFUN([ACX_USE_ALTERNATIVE_LIBM],
  [m4_foreach([ACX_LibSet],[[LIBS],[FCLIBS]],
     [ACX_LibSet=`echo "$ACX_LibSet" | $PERL -e 'while(<>) { s{\s\K-lm(?=\s|$)}{'"$LIBM"'}g ; print }'`
     ])
   m4_ifval($1,[m4_foreach([ACX_LibSet],[$1],
        [ACX_LibSet=`echo "$ACX_LibSet" | $PERL -e 'while(<>) { s{\s\K-lm(?=\s|$)}{'"$LIBM"'}g ; print }'`
     ])])
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
