# _ACX_SL_FPP_BUILD_RULE
# ------------------
# Figure out how to build from cpp/Fortran sources
#
# If we need to use a separate preprocessor, we must override make's
# `direct' .F.o rule in order to do `indirect' compilation
# (.F -> .f then .f -> .o).
#
# Configure variables set here are as follows.  The items in this list
# are suffixed with `[direct]', `[indirect]' or `[both]'.  In the
# first two cases, the variable has a useful value only in the given
# mode, and an unspecified, and therefore unreliable, value in the
# other; in the last, it has a value in both modes.
#
#   FPP [indirect]
#     The name of a suitable preprocessor.
#
#   FPP_COMPILE_EXT [both]
#     This contains the file extension which the Fortran compiler will
#     accept as containing source not to be preprocessed.  It is most
#     typically 'f' (the default), but could be different if set by a
#     call to ACX_SL_FPP_(FIXED|FREE)FORM.
#
#   FPP_PREPROCESS_EXT [both]
#     The partner of FPP_COMPILE_EXT, containing the file extension
#     which is taken to indicate Fortran source to be preprocessed.
#     The default is 'F', but could be different if set by a
#     call to ACX_SL_FPP_(FIXED|FREE)FORM.
#
#   FPP_MAKE_FLAGS [direct]
#     This is used to include CPP/FPP related flags into the compiler
#     call if we compile directly, and leave them out otherwise.
#
#   FPP_OUTPUT [both]
#     This is used to redirect FPP output to the .f file in those
#     cases where FPP writes to stdout rather than to a file.  It is
#     defined as either "" or ">$@".
#
#   FPPDIRECT_TRUE, FPPDIRECT_FALSE [both]
#     If the macro decides that we must use `direct' mode, then it
#     sets FPPDIRECT_TRUE to be blank, and FPPDIRECT_FALSE to be '#',
#     or vice versa if we are to use `indirect' mode.  These may be
#     used within a Makefile.in as follows:
#       @FPPDIRECT_TRUE@.@FPP_PREPROCESS_EXT@.o:
#       @FPPDIRECT_TRUE@        $(PPFCCOMPILE) -c -o $@ $<
#       @FPPDIRECT_FALSE@.@FPP_PREPROCESS_EXT@.@FPP_COMPILE_EXT:
#       @FPPDIRECT_FALSE@        $(FPP) $(DEFS) ... $< @FPP_OUTPUT@
#     If you use automake, then you may possibly recognise that as an
#     automake conditional (which is predeclared, so you do not need
#     to include AM_CONDITIONAL(FPPDIRECT, ???) in your configure.ac),
#     which might be used more straightforwardly in your Makefile.am
#     as follows:
#       if FPPDIRECT
#       .@FPP_PREPROCESS_EXT@.o:
#               $(PPFCCOMPILE) -c -o $@ $<
#       else !FPPDIRECT
#       .@FPP_PREPROCESS_EXT@.@FPP_COMPILE_EXT:
#               $(FPP) $(DEFS) ... $< @FPP_OUTPUT@
#       endif !FPPDIRECT
#
# These are used in Automake's lang_ppfc_finish subroutine.
#
# NOTE 1: There would seem to be a problem here with the use of .F as
# the extension for preprocessed files.  On case-insensitive
# filesystems such as HFS+, as used on MacOS X, foo.F and foo.f are
# the same file.  This means that indirect compilation would lose badly, since
# converting foo.F to foo.f would clobber the original.  This is
# probably not a problem in practice, since the compilers (g77, gfortran,
# nag, and xlf) actually likely to be used on OS X -- which is a
# recent platform, and thus with only recent Fortrans on it -- can all
# do direct compilation of preprocessable Fortran.  Just in case, we
# check below whether we are in this fatal situation, and collapse
# noisily if necessary.
#
# NOTE 2: Martin Wilck's original version of these macros noted that it
# was necessary to generate explicit rules for .F -> .o compilations
# in order to override make's builtin rules in a portable manner
# (i.e. without using make extensions).  Not all makes do chains of
# implicit rules, so we cannot depend on .F.f, .f.o rules generating
# a .f file.  We need unified .F.o and .F.lo rules, but that's
# complicated, an alternative is to name the intermediary .f files in
# the Makefiles.  Again, this may not be much of a problem in fact,
# since the main culprit seems to be Solaris make, but Solaris f77
# can do direct compilation, so that the issue of chaining rules
# doesn't arise.
#
# NOTE 3: POSIX/Single-Unix states that inference rules can be
# redefined, and there's no warning against this in Autoconf's section
# on `Limitations of Make'.
#
# NOTE 4: FPP_OUTPUT is set to either "" or ">$@".  The latter is OK
# in an implicit rule, but will potentially lose in an explicit rule,
# since POSIX does not require that $@ is defined in such a rule, and
# there are still a few makes which do not define it in that context.
# As with Note 1, however, this is probably more a theoretical problem
# than a practical one.
#
AC_DEFUN([_ACX_SL_FPP_BUILD_RULE],
[# FPP is defined by this stage.  If the processing mode is 'direct', then
# this will almost certainly be defined as blank, but we should make no
# committments to this in the documentation, in case we want to change
# our minds about that in future.
AC_SUBST(FPP)

# Default the FPP_PREPROCESS_EXT and FPP_COMPILE_EXT to the most usual ones
FPP_PREPROCESS_EXT=F
FPP_COMPILE_EXT=f

# Switch on the processing mode, direct/indirect, which has been determined
# in AC_PROG_FPP before this macro is called.  The FPPDIRECT_(TRUE|FALSE)
# variables implement an automake (configure-time) conditional, which is
# created, not through an invocation of AM_CONDITIONAL, but implicitly
# within automake.in (qv).
if test $ac_cv_fpp_build_rule = direct; then
   # The simple case: the chosen Fortran compiler can handle preprocessing,
   # so we don't need a separate preprocessing stage.
   FPPDIRECT_TRUE=
   FPPDIRECT_FALSE='#'
   # The flags here are those included in the 'compile' field of the
   # 'ppfc' language in automake.in, minus the {AM_,}FCFLAGS variables.
   # It's not _absolutely_ guaranteed that these are the correct ones,
   # and I (NG) would be open to argument about adding both {AM_,}CPPFLAGS and
   # {AM_,}FCFLAGS, but this set appears to work.
   FPP_MAKE_FLAGS='$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS)'
else
   FPPDIRECT_TRUE='#'
   FPPDIRECT_FALSE=
   FPP_MAKE_FLAGS=
fi

if test -z "$ac_fpp_out"; then
   FPP_OUTPUT=" "
else
   FPP_OUTPUT=">\[$]@"
fi

AC_SUBST(FPPDIRECT_TRUE)
AC_SUBST(FPPDIRECT_FALSE)
AC_SUBST(FPP_MAKE_FLAGS)
AC_SUBST(FPP_PREPROCESS_EXT)
AC_SUBST(FPP_COMPILE_EXT)
AC_SUBST(FPP_OUTPUT)
])# _ACX_SL_FPP_BUILD_RULE
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End: