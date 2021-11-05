!+ Data module for global KIND type parameters
!-------------------------------------------------------------------------------

MODULE data_parameters

!-------------------------------------------------------------------------------
!
! Description:
!  Global parameters for defining the KIND types of the real- and integer-
!  variables are defined.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.8        1998/08/03 Ulrich Schaettler
!  Eliminated intgribf, intgribc, irealgrib, iwlength and put it to data_io.
! 1.10       1998/09/29 Ulrich Schaettler
!  Eliminated parameters for grid point and diagnostic calculations.
! 3.13       2004/12/03 Ulrich Schaettler
!  Introduced intgribf, intgribc, irealgrib, iwlength (again)
! 3.18       2006/03/03 Ulrich Schaettler
!  Introduced KIND parameters idouble, isingle for generic formulation of
!  some utility routines
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced: use MESSy KIND definitions
! V4_28        2013/07/12 Ulrich Schaettler
!  Implemented KIND parameters int_ga for grib_api interface (number of bytes,
!   which could be 4 or 8 byte integers)
!  Implemented global KIND parameter int_dp for 8 byte integers
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure:
!   can remove all ifdef MESSY here
! V5_00_clm8   2016/01/29 Ulrich Schaettler, H.-J. Panitz (HJP)
!  Increase value of SELECTED_INT_KIND from 12 to 18 
!   (suggested by U. Schaettler) in order to avoid a
!    possible crash of the model after a total simulation time of 68 years
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

#ifdef GRIBAPI
USE grib_api
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! 1. KIND-Parameters for the Program:
! -----------------------------------

  INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers

! 2. KIND-Parameters for the variables in the GRIB-library
! --------------------------------------------------------

  INTEGER, PARAMETER       ::                                         &
    intgribf  = KIND(1),                                              &
!   intgribf  = 4,      &  ! (if using libgrib1 on the T3E)
       ! Kind type for Fortran integer variables used in the GRIB library
       ! this normally is the Standard integer with the exception of using
       ! "libgrib1" (former supplib) on a machine with 8 byte INTEGER default
       ! (like Cray-machines; then intgribf has be set to 32-bit INTEGER).

    intgribc  = KIND(1),                                              &
       ! Kind type for C integer variables used in the GRIB library
       ! this always is the Standard integer

    irealgrib = KIND(1.0)
       ! Kind type for Fortran real variables used in the GRIB library
       ! this is the Standard real of the machine


  INTEGER                  ::                                         &
    ! this variable has to be set at the beginning of the program
    ! (at the beginning of organize_data)
    iwlength   ! length of integers used in the griblib in byte
               ! 8: for dwdlib on Cray PVP and T3E systems
               ! 4: for dwdlib on SGI systems and griblib on all systems

! 3. KIND-Parameters for the generic formulation of some utility routines:
! ------------------------------------------------------------------------

    ! The distinction between ireals (working precision) and irealgrib
    ! is not enough, because it could be possible that these KIND parameters
    ! are the same. Compilers could get in trouble then, because they could
    ! not decide, which routine to take then.
    ! Therefore we define the KIND parameters idouble (for double precision
    ! or 8 byte reals) and isinge (for single precision, or 4 byte reals)

  INTEGER, PARAMETER       ::                                         &
       idouble   = KIND (1.0D0),                                      &
       isingle   = KIND (1.0)

! 4. KIND-Parameters for different INTEGER precision:
! ---------------------------------------------------

  INTEGER, PARAMETER       ::                                         &
!
!HJP 2016-01-29 Begin
!
!  Increase value of SELECTED_INT_KIND from 12 to 18 
!   (suggested by U. Schaettler) in order to avoid a
!    possible crash of the model after a total simulation time of 68 years
!      int_dp    = SELECTED_INT_KIND (12),                            &
               ! should represent integers up to 10**12
               ! which should be a INTEGER*8 (in the old notation)
       int_dp    = SELECTED_INT_KIND (18),                            &
!HJP 2016-01-29 End

! integer precision necessary for grib_api in interfaces where length of
! message in bytes is involved
#ifdef GRIBAPI
       int_ga    = kindOfSize              ! should be INTEGER *8 where necessary
#else
       int_ga    = SELECTED_INT_KIND (8)   ! INTEGER *4
#endif

!==============================================================================

END MODULE data_parameters
