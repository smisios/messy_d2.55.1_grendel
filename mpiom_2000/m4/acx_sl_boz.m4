# ACX_SL_FC_HAVE_BOZ
# --------------
# Test whether the FC compiler supports BOZ constants in the Fortran
# 95 style. These are integer constants written in the format
# B'xxx', O'xxx' and Z'xxx'. If so set the preprocessor variable
# HAVE_BOZ to be 1.
AC_DEFUN([ACX_SL_FC_HAVE_BOZ],
  [AC_REQUIRE([AC_PROG_FC])dnl
     AC_CACHE_CHECK([whether ${FC} supports F95 BOZ constants],
       [ac_cv_fc_have_boz],
       [AC_LANG_PUSH([Fortran])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
      INTEGER J, K, L
      PARAMETER ( J = B'1111111111111111' )
      PARAMETER ( K = O'7777' )
      PARAMETER ( L = Z'FF' )
          ])],dnl
          [ac_cv_fc_have_boz=yes],
          [ac_cv_fc_have_boz=no])
        AC_LANG_POP([Fortran])])
     AS_IF([test x$ac_cv_fc_have_boz = xyes],
       [AC_DEFINE([HAVE_BOZ], 1,dnl
          [Define to 1 if the Fortran compiler supports F95 boz constants])])dnl
])# ACX_SL_FC_HAVE_BOZ


# ACX_SL_FC_HAVE_TYPELESS_BOZ
# ---------------------------
# Test whether the FC compiler supports typeless BOZ constants in the Fortran
# 95 style. These are (usually integer) constants written in the format
# X'xxx', but may also be typeless, which allows the initialisation of any type
# to a specific bit pattern. If so set the preprocessor variable
# HAVE_TYPELESS_BOZ to be 1.
#
# A problem with this test is that it may compile but the assignments are not
# actually typeless and have been done as integer casts.  To stop this we need
# to run the program and check if an integer equals a floating point value, by
# value, they shouldn't for a bit pattern assignment. Uses a "EXIT(1)" to
# signal a problem. This is non-standard, so the test may fail for that reason.
AC_DEFUN([ACX_SL_FC_HAVE_TYPELESS_BOZ],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether ${FC} supports F95 typeless BOZ constants],
     [ac_cv_fc_have_typeless_boz],
     [AC_LANG_PUSH([Fortran])
      AC_RUN_IFELSE([AC_LANG_SOURCE([
      PROGRAM TMP
      INTEGER I
      PARAMETER ( I = X'FF7FFFFF' )
      REAL D
      PARAMETER ( D = X'FF7FFFFF' )
      LOGICAL L
      PARAMETER ( L = X'A55A5AA5' )
      IF ( D .EQ. I ) THEN
         CALL EXIT( 1 )
      END IF
      END
       ])],
        [ac_cv_fc_have_typeless_boz=yes],
        [ac_cv_fc_have_typeless_boz=no])
      AC_LANG_POP([Fortran])])
  AS_IF([test x$ac_cv_fc_have_typeless_boz = xyes],
    [AC_DEFINE([HAVE_TYPELESS_BOZ], 1,dnl
       [Define to 1 if the Fortran compiler supports F95 typeless boz constants])])
])# ACX_SL_FC_HAVE_TYPELESS_BOZ


# ACX_SL_FC_HAVE_OLD_TYPELESS_BOZ
# ---------------------------
# Test whether the FC compiler supports typeless BOZ constants in the OLD (VMS
# and g77) Fortran style. These are constants written in the format 'xxx'X,
# which allows the initialisation of any type to a specific bit pattern. If so
# set the preprocessor variable HAVE_OLD_TYPELESS_BOZ to be 1.
AC_DEFUN([ACX_SL_FC_HAVE_OLD_TYPELESS_BOZ],
         [AC_REQUIRE([AC_PROG_FC])dnl
          AC_CACHE_CHECK([whether ${FC} supports OLD style typeless BOZ constants],
                         [ac_cv_fc_have_old_typeless_boz],
                         [AC_LANG_PUSH([Fortran])
                          AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
      INTEGER I
      PARAMETER ( I = 'FF'X )
      DOUBLE PRECISION D
      PARAMETER ( D = 'FFEFFFFFFFFFFFFF'X )
])],
                          ac_cv_fc_have_old_typeless_boz=yes,
                          ac_cv_fc_have_old_typeless_boz=no)
                          AC_LANG_POP([Fortran])])
          if test $ac_cv_fc_have_old_typeless_boz = yes; then
              AC_DEFINE([HAVE_OLD_TYPELESS_BOZ], 1,
                        [Define to 1 if the Fortran compiler supports OLD style typeless boz constants])
          fi
])# ACX_SL_FC_HAVE_OLD_TYPELESS_BOZ
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
