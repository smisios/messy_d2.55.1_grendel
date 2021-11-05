dnl
dnl ACX_FORTRAN_RUN_CHECK_SIZEOF(var_for_size, )
dnl
dnl sets var_for_size to the size.  Ignores if the size cannot be determined
dnl
AC_DEFUN([ACX_FORTRAN_RUN_CHECK_SIZEOF],
  [AC_REQUIRE([ACX_C_CHAR_BITS])
   AS_VAR_PUSHDEF([acx_fortran_Type], [acx_cv_fortran_type_$1])dnl
   ACX_FORTRAN_CHECK_TYPE([$1])
   AS_IF([test x]AS_VAR_GET([acx_fortran_Type])[ = xyes],[
     AS_VAR_PUSHDEF([acx_fortran_Sizeof],[acx_cv_fortran_sizeof_$1])dnl
     AC_CACHE_CHECK([size of Fortran type $1], [acx_fortran_Sizeof],dnl
       [AC_LANG_PUSH([Fortran])
        AC_RUN_IFELSE(AC_LANG_PROGRAM([], [
    integer :: itest
    real    :: rtest
    integer :: integer_io_length
    integer :: real_io_length
    integer :: integer_bits
    integer :: real_bits
    inquire(iolength=integer_io_length) itest
    inquire(iolength=real_io_length) rtest
    integer_bits=bit_size(itest)
    real_bits = real_io_length * integer_bits / integer_io_length
    open(10,file="conftestval")
    write(10,*)real_bits
    close(10)
]),dnl
          [AS_VAR_SET([acx_fortran_Sizeof], `cat conftestval`)
           AS_VAR_SET([acx_fortran_Sizeof],[`expr ]AS_VAR_GET([acx_fortran_Sizeof])[ / $acx_cv_c_char_bits`])],
          [AS_VAR_SET([acx_fortran_Sizeof], [unavailable])])
        /bin/rm -f conftestval
        AC_LANG_POP([Fortran])
       ])
     AS_IF([test x"]AS_VAR_GET([acx_fortran_Sizeof])[" = xunavailable], [$2], [$3])
     m4_ifval([$2], $2=AS_VAR_GET([acx_fortran_Sizeof]))
     AS_VAR_POPDEF([acx_fortran_Sizeof])dnl
  ])dnl
  AS_VAR_POPDEF([acx_fortran_Type])dnl
])
dnl Local Variables:
dnl mode: autoconf
dnl End:
