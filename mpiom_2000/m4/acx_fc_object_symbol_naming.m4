dnl
dnl Get the format of Fortran names.  Uses F77, FFLAGS, and sets WDEF.
dnl If the test fails, sets NOF77 to 1, HAS_FORTRAN to 0
dnl
AC_DEFUN([AC_GET_FORTNAMES],[
dnl Check for strange behavior of Fortran.  For example, some FreeBSD
dnl systems use f2c to implement f77, and the version of f2c that they
dnl use generates TWO (!!!) trailing underscores
dnl Currently, WDEF is not used but could be...
dnl
dnl Eventually, we want to be able to override the choices here and
dnl force a particular form.  This is particularly useful in systems
dnl where a Fortran compiler option is used to force a particular
dnl external name format (rs6000 xlf, for example).
AC_FC_WRAPPERS
dnl cat > conftest.f <<EOF
dnl       subroutine mpir_init_fop( a )
dnl       integer a
dnl       a = 1
dnl       return
dnl       end
dnl EOF
dnl $F90 $FFLAGS -c conftest.f > /dev/null 2>&1
dnl if test ! -s conftest.o ; then
dnl     echo "Unable to test Fortran compiler"
dnl     echo "(compiling a test program failed to produce an "
dnl     echo "object file)."
dnl     HAS_FORTRAN=0
dnl elif test -z "$FORTRANNAMES" ; then
dnl dnl We have to be careful here, since the name may occur in several
dnl dnl forms.  We try to handle this by testing for several forms
dnl dnl directly.
dnl     if test $arch_CRAY ; then
dnl dnl Cray doesn't accept -a ...
dnl         nameform1=`strings conftest.o | grep mpir_init_fop_  | sed -n -e '1p'`
dnl         nameform2=`strings conftest.o | grep MPIR_INIT_FOP   | sed -n -e '1p'`
dnl         nameform3=`strings conftest.o | grep mpir_init_fop   | sed -n -e '1p'`
dnl         nameform4=`strings conftest.o | grep mpir_init_fop__ | sed -n -e '1p'`
dnl     else
dnl         nameform1=`strings -a conftest.o | grep mpir_init_fop_  | sed -n -e '1p'`
dnl         nameform2=`strings -a conftest.o | grep MPIR_INIT_FOP   | sed -n -e '1p'`
dnl         nameform3=`strings -a conftest.o | grep mpir_init_fop   | sed -n -e '1p'`
dnl         nameform4=`strings -a conftest.o | grep mpir_init_fop__ | sed -n -e '1p'`
dnl     fi
dnl     /bin/rm -f conftest.f conftest.o
dnl     if test -n "$nameform4" ; then
dnl 	echo "Fortran externals are lower case and have 1 or 2 trailing underscores"
dnl 	FORTRANNAMES="FORTRANDOUBLEUNDERSCORE"
dnl     elif test -n "$nameform1" ; then
dnl         echo "Fortran externals have a trailing underscore and are lowercase"
dnl 	FORTRANNAMES="FORTRANUNDERSCORE"
dnl     elif test -n "$nameform2" ; then
dnl 	echo "Fortran externals are uppercase"
dnl 	FORTRANNAMES="FORTRANCAPS"
dnl     elif test -n "$nameform3" ; then
dnl 	echo "Fortran externals are lower case"
dnl 	FORTRANNAMES="FORTRANNOUNDERSCORE"
dnl     else
dnl 	echo "Unable to determine the form of Fortran external names"
dnl 	echo "Make sure that the compiler $F90 can be run on this system"
dnl         HAS_FORTRAN=0
dnl     fi
dnl fi
])dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
