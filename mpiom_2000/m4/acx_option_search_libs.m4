dnl _ACX_OPTION_SEARCH_LIBS(FUNCTION, SEARCH-LIBS,
dnl   [OTHER-LIBRARIES], [EXTRA-FLAGS])
dnl same as ACX_OPTION_SEARCH_LIBS but does not cache result
AC_DEFUN([_ACX_OPTION_SEARCH_LIBS],
  [acx_option_func_search_save_LIBS="$LIBS"
   AC_LANG_CONFTEST([AC_LANG_CALL([], [$1])])
   for ac_lib in '' $2; do
     AS_IF([test -z "$ac_lib"],
       [ac_res="none required"
        LIBS="m4_ifval([$4],[$4 ])m4_ifnblank($3,[$3 ])$acx_option_func_search_save_LIBS"],
       [ac_res="-l$ac_lib"
        LIBS="m4_ifval([$4],[$4 ])$ac_res m4_ifnblank($3,[$3 ])$acx_option_func_search_save_LIBS"])
     AC_LINK_IFELSE([], [AS_IF([test x"$ac_res" = x"none required"],dnl
        [AS_VAR_SET([ac_Search],["]m4_ifval([$4],[$4 ])[$3"])],dnl
        [AS_VAR_SET([ac_Search],["]m4_ifval([$4],[$4 ])[-l$ac_lib $3"])])])
     AS_VAR_SET_IF([ac_Search], [break])
   done
   rm conftest.$ac_ext
   LIBS="$acx_option_func_search_save_LIBS"])
dnl ACX_OPTION_SEARCH_LIBS(FUNCTION, SEARCH-LIBS,
dnl   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND], [OTHER-LIBRARIES],
dnl   [EXTRA-FLAGS])
dnl --------------------------------------------------------
dnl Search for a library defining FUNC, if it's not already available.
dnl Do not add said library to default LIBS output variable.
dnl Use $ac_lib or $ac_res in ACTION if needed. Uses OTHER-LIBRARIES
dnl unconditionally, which might provoke linker errors.
AC_DEFUN([ACX_OPTION_SEARCH_LIBS],
  [AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_$1_]_AC_LANG_ABBREV)dnl
   AC_CACHE_CHECK([for library containing $1], [ac_Search],
     [_ACX_OPTION_SEARCH_LIBS([$1],[$2],[$5],[$6])])
   ac_res=AS_VAR_GET([ac_Search])
   AS_VAR_SET_IF([ac_Search],
     [$3],
     [$4])
   AS_VAR_POPDEF([ac_Search])dnl
])
dnl ACX_OPTION_SEARCH_LIBS_MULTI(FUNCTION, SEARCH-LIBS,
dnl   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND], [OTHER-LIBRARIES],
dnl   [EXTRA-FLAGS])
dnl --------------------------------------------------------
dnl Search for a library defining FUNC, if it's not already available.
dnl Do not add said library to default LIBS output variable.
dnl Use $ac_lib or $ac_res in ACTION if needed.
dnl Tries first to link without any OTHER-LIBRARY, then successively
dnl tries each library in the list.
AC_DEFUN([ACX_OPTION_SEARCH_LIBS_MULTI],
  [AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_$1_]_AC_LANG_ABBREV)dnl
   AC_MSG_CHECKING([for library containing $1])
   AC_CACHE_VAL([ac_Search],dnl
     [while :; do
        m4_if(m4_car($5),[[]],,[_ACX_OPTION_SEARCH_LIBS([$1],[$2],,[$6])
        AS_VAR_SET_IF([ac_Search], [break])])
        m4_foreach([ACX_LibSet], [$5],dnl
          [_ACX_OPTION_SEARCH_LIBS([$1],[$2],[ACX_LibSet],[$6])
           AS_VAR_SET_IF([ac_Search],[break])
          ])
        break
      done])
   AS_VAR_SET_IF([ac_Search],dnl
     [AS_IF([test x"AS_VAR_GET([ac_Search])" = x],dnl
        [AC_MSG_RESULT([(none required)])],dnl
        [AC_MSG_RESULT([AS_VAR_GET([ac_Search])])])],dnl
     [AC_MSG_RESULT([not found])])
   AS_VAR_SET_IF([ac_Search],[$3],[$4])
   AS_VAR_POPDEF([ac_Search])dnl
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
