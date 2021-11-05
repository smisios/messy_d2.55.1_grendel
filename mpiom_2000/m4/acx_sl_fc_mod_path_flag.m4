# ACX_SL_FC_MOD_PATH_FLAG([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
# -------------------------
# Check which flag is necessary to alter the compiler's search path
# for module files.
# This obviously requires that the compiler has some notion of
# module files as separate from object files and some sensible
# method of altering its search path. This will therefore not work
# on early Cray F90 compilers, or on v5 (and 6?) of ifc.
#
# Nearly every compiler I have found uses -Ipath for this purpose;
# Sun F95 v7.1 (at least), uses -Mpath
#
AC_DEFUN([ACX_SL_FC_CHECK_MOD_PATH_FLAG],dnl
  [ACX_ASSERT_LANG_IS_FORTRAN_VARIANT
   AC_REQUIRE([AC_PROG_FC])
   AS_VAR_PUSHDEF([mod_flag],[acx_sl_cv_fc_mod_path_flag_]_AC_LANG_ABBREV)dnl
   ASX_VAR_UNSET([mod_flag])
   AC_CACHE_CHECK([for flag to alter module search path],[mod_flag],dnl
     [mkdir conftestdir
      cd conftestdir
      AC_COMPILE_IFELSE([AC_LANG_SOURCE([      module cnftst
       implicit none
       integer :: i
      end module cnftst])],,
        [AC_MSG_ERROR([Cannot compile fortran modules])])
      cd ..
      for i in -I -M -module -p; do
        FCFLAGS_save=$FCFLAGS
        FCFLAGS="$FCFLAGS ${i}conftestdir"
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[dnl
       use cnftst
       implicit none
       i = 0
])
             ],[AS_VAR_SET([mod_flag],[$i]) ; FCFLAGS=$FCFLAGS_save ; break],
             [:])
        FCFLAGS=$FCFLAGS_save
      done
      FCFLAGS=$FCFLAGS_save
      rm -rf conftestdir
      AS_VAR_SET_IF([mod_flag],dnl
        [m4_default([$1],[:])],dnl
        [m4_default([$2],[AC_MSG_ERROR([Cannot find flag to alter module search path])])])])
   FC_MOD_FLAG=AS_VAR_GET([mod_flag])
   AC_SUBST([FC_MOD_FLAG])
   AS_VAR_POPDEF([mod_flag])dnl
])# ACX_SL_FC_MOD_PATH_FLAG
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
