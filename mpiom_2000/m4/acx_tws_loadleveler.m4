AC_DEFUN([ACX_TWS_LOADLEVELER],
  [ibm_tws_load_leveler=no
   AC_SUBST([TWS_LL_LIB])
   AC_SUBST([TWS_LL_INCLUDE])
   AC_ARG_WITH([loadleveler-root],
     [AS_HELP_STRING([--with-loadleveler-root],
        [specify installation directory of IBM TWS LoadLeveler])],
        [TWS_LL_ROOT=$with_loadleveler_root])
   AC_ARG_WITH([loadleveler-include],
     [AS_HELP_STRING([--with-loadleveler-include],
        [specify installation directory of IBM TWS LoadLeveler llapi.h header])],
        [TWS_LL_INCLUDE=-I$with_loadleveler_include],
        [AS_VAR_SET_IF([TWS_LL_INCLUDE],,dnl
           [AS_VAR_SET_IF([TWS_LL_ROOT],dnl
              [TWS_LL_INCLUDE="-I$TWS_LL_ROOT/include"])])])
   AC_ARG_WITH([loadleveler-lib],
     [AS_HELP_STRING([--with-loadleveler-lib],
        [specify installation directory of IBM TWS LoadLeveler llapi.h header])],
        [TWS_LL_LIB=-L$with_loadleveler_lib],
        [AS_VAR_SET_IF([TWS_LL_LIB],,dnl
           [AS_VAR_SET_IF([TWS_LL_ROOT], [TWS_LL_LIB="-L$TWS_LL_ROOT/lib"])])])

   ac_save_LIBS=$LIBS
   ac_save_CPPFLAGS=$CPPFLAGS
   AS_VAR_SET_IF([TWS_LL_INCLUDE],
     [CPPFLAGS="$CPPFLAGS $TWS_LL_INCLUDE"])
   AS_VAR_SET_IF([TWS_LL_LIB],
     [LIBS="$TWS_LL_LIB $LIBS"])
   AC_CHECK_HEADER([llapi.h],dnl
     [AC_CHECK_LIB([llapi],[ll_query],dnl
        [ibm_tws_load_leveler=yes
         TWS_LL_LIB="$TWS_LL_LIB -lllapi"])])
   LIBS=$ac_save_LIBS
   CPPFLAGS=$ac_save_CPPFLAGS
   AM_CONDITIONAL([HAVE_IBM_TWS_LOAD_LEVELER],dnl
     [test x"$ibm_tws_load_leveler" = xyes])dnl
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
