dnl
dnl Check for doxygen properties
dnl 
AC_DEFUN([ACX_CHECK_DOXYGEN_LAYOUT],
[AC_CHECK_PROGS([DOXYGEN], [doxygen], [:])
AS_IF([test x"$DOXYGEN" = "x:"],[DOXYFILE_HEADER=Doxyfile-1.5.6],
  [AC_MSG_CHECKING([if doxygen supports layouting])
   AS_IF(["$DOXYGEN" -l /dev/null > /dev/null 2>&1],
     [AC_MSG_RESULT([yes])
      DOXYFILE_HEADER=Doxyfile-1.5.8],
     [AC_MSG_RESULT([no])
      DOXYFILE_HEADER=Doxyfile-1.5.6])])
])dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
