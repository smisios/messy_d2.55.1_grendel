!+ External procedure organize_eps for organizing the ensemble prediction mode
!------------------------------------------------------------------------------

SUBROUTINE organize_eps(yaction, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!
!   This procedure is the driving routine for EPS mode set-up
!
!   Internal Routines/Functions currently contained:
!
!    - subroutine input_eps
!      This subroutine reads, checks and distributes the NAMELIST input
!      associated with the EPS
!
! Method:
!   Determine whether a certain action has to be performed in this time step
!
! Current Code Owner: DWD/FE15, Susanne Theis
!  phone:  +49  069 8062 2741
!  fax:    +49  69  8062 3721
!  email:  susanne.theis@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.21       2006/12/04 Susanne Theis
!  Initial release
! V4_5         2008/09/10 Christoph Gebhardt
!  Add namelist parameters for modifying values of lai, plcov, rootdp
! V4_8         2009/02/16 Ulrich Schaettler
!  Corrections in the settings of default values
! V4_10        2009/09/11 Ulrich Schaettler
!  Corrections in the settings of rmin_plcov
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V5_00_clm9   2016/05/11 H.-J. Panitz, IMK/KIT
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!   adapted from COSMO_5.1
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,    ONLY :  &
    iintegers,        & ! KIND-type parameter for integer variables
    ireals              ! KIND-type parameters for real variables

USE parallel_utilities, ONLY :  &
    distribute_values     ! subroutine for MPI value distribution

USE data_parallel,      ONLY :  &
    my_world_id,   & ! rank of this subdomain in the global communicator
    icomm_world,   & ! communicator that belongs to igroup_world, i.e.
                     ! = MPI_COMM_WORLD
    nproc,         & ! total number of processors
    intbuf,        & ! buffers for distributing the Namelist
    realbuf,       & !
    logbuf,        & !
    imp_reals,     & !correct type for MPI
    imp_integers,  & !
    imp_logical

USE data_runcontrol,    ONLY :  &
    iepsmem, iepstyp, iepstot, & ! EPS Member-Id, Typ-Id
                                 ! and total number of members
    fac_plcov,                 & ! modification factor for PLCOV
    rmin_plcov,                & ! lower limit of PLCOV
    rmax_plcov,                & ! upper limit of PLCOV
    fac_rootdp,                & ! modification factor for ROOTDP
    rmin_rootdp,               & ! lower limit of ROOTDP
    rmax_rootdp,               & ! upper limit of ROOTDP
    fac_lai,                   & ! modification factor for LAI
    rmin_lai,                  & ! lower limit of LAI
    rmax_lai,                  & ! upper limit of LAI
    nuspecif                     ! file number of YUSPECIF

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Parameter list:
CHARACTER (LEN= *),       INTENT(IN)            ::                      &
  yaction      ! action to be performed

INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
  ierror       ! error status

CHARACTER (LEN= *),       INTENT(OUT)           ::                      &
  yerrmsg      ! error message


! local variables:

  INTEGER (KIND=iintegers) ::                       &
    nuin,                     &! unit-number of INPUT-File
    izerrstat                  ! error control


  CHARACTER (LEN= 9)       ::                       &
    yinput                     ! Namelist INPUT file

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine organize_eps
!------------------------------------------------------------------------------

ierror  = 0_iintegers
yerrmsg = '   '

!------------------------------------------------------------------------------
! Section 1: Input of the Namelist
!------------------------------------------------------------------------------

IF (yaction == 'input') THEN

  ! Open NAMELIST-INPUT file
  IF (my_world_id == 0) THEN
    PRINT *,'    INPUT OF THE NAMELISTS FOR ENSEMBLE'
    yinput   = 'INPUT_EPS'
    nuin     =  1
    OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=izerrstat)
    IF(izerrstat /= 0) THEN
      yerrmsg  = ' ERROR    *** Error while opening file INPUT_EPS *** '
      ierror   = 2001
      RETURN
    ENDIF
  ENDIF

  ! Read the NAMELIST-group
  CALL input_epsctl (nuspecif, nuin, izerrstat)

  IF (my_world_id == 0) THEN
    ! Close file for input of the NAMELISTS
    CLOSE (nuin    , STATUS='KEEP')
  ENDIF

  IF (izerrstat < 0) THEN
    yerrmsg = 'ERROR *** while reading NAMELIST Group /EPSCTL/ ***'
    ierror  = 2003
  ELSEIF (izerrstat > 0) THEN
    yerrmsg = 'ERROR *** Wrong values occured in NAMELIST INPUT_EPS ***'
    ierror  = 2004
  ENDIF


ENDIF

!------------------------------------------------------------------------------
! Internal Procedures
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Internal procedure in "organize_diagnosis" for the input of NAMELIST diactl
!------------------------------------------------------------------------------

SUBROUTINE input_epsctl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine checks and distributes the NAMELIST input
!   associated with the EPS
!
! Method:
!   copied from organize_physics.f90 (input_phyctl)
!   and modified according to the current needs
!
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program
!   after reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------

! Parameter list:
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! local variables:
!
! Variables for default values
  INTEGER (KIND=iintegers) ::                       &
    iepsmem_d,               & ! Ensemble member-ID
    iepstot_d,               & ! total number of ensemble members
    iepstyp_d                  ! Ensemble typ-ID

  REAL    (KIND=ireals)    ::                       &
    fac_plcov_d,          & ! modification factor for PLCOV
    rmin_plcov_d,         & ! lower limit of PLCOV
    rmax_plcov_d,         & ! upper limit of PLCOV
    fac_rootdp_d,         & ! modification factor for ROOTDP
    rmin_rootdp_d,        & ! lower limit of ROOTDP
    rmax_rootdp_d,        & ! upper limit of ROOTDP
    fac_lai_d,            & ! modification factor for LAI
    rmin_lai_d,           & ! lower limit of LAI
    rmax_lai_d              ! upper limit of LAI

! Miscellaneous variables
  INTEGER (KIND=iintegers) ::                       &
    ierr                       ! error variable for distribute_values
                               ! in this form useless (no input/output)
  CHARACTER (LEN= 8)       ::                       &
    yinput                     ! Namelist INPUT file

!HJP Begin 2016-05-11
  CHARACTER(LEN=250)         :: iomsg_str
!HJP End   2016-05-11

! Define the namelist group
  NAMELIST /epsctl/ iepsmem, iepstot, iepstyp, fac_lai, rmin_lai, rmax_lai,  &
                    fac_plcov, rmin_plcov, rmax_plcov,                       &
                    fac_rootdp, rmin_rootdp, rmax_rootdp

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_epsctl
!------------------------------------------------------------------------------

ierrstat = 0_iintegers

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  iepsmem_d    = -(1_iintegers)
  iepstot_d    = -(1_iintegers)
  iepstyp_d    = -(1_iintegers)

  fac_plcov_d    = 1._ireals
  rmin_plcov_d   = 0._ireals
  rmax_plcov_d   = 1._ireals
  fac_rootdp_d   = 1._ireals
  rmin_rootdp_d  = 0._ireals
  rmax_rootdp_d  = 2._ireals
  fac_lai_d      = 1._ireals
  rmin_lai_d     = 0._ireals
  rmax_lai_d     = 8._ireals

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!------------------------------------------------------------------------------

  iepsmem = iepsmem_d
  iepstot = iepstot_d
  iepstyp = iepstyp_d

  fac_plcov      = fac_plcov_d
  rmin_plcov     = rmin_plcov_d
  rmax_plcov     = rmax_plcov_d
  fac_rootdp     = fac_rootdp_d
  rmin_rootdp    = rmin_rootdp_d
  rmax_rootdp    = rmax_rootdp_d
  fac_lai        = fac_lai_d
  rmin_lai       = rmin_lai_d
  rmax_lai       = rmax_lai_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

!HJP Begin 2016-05-11
! READ (nuin, epsctl)
  iomsg_str(:) = ' '
  READ (nuin, epsctl, IOSTAT=ierr, IOMSG=iomsg_str)

  IF (ierr /= 0) WRITE (*,'(A,A)') 'Namelist-ERROR EPSCTL: ', TRIM(iomsg_str)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (ierr, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (ierr /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!HJP END   2016-05-11

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

! Check whether values are OK
  IF ( iepsmem < 0 .OR. iepstyp < 0 .OR. iepstot < 0 .OR. &
       iepsmem > iepstot ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** iepsmem, iepstot, iepstyp *** ',iepsmem,iepstot,iepstyp
    ierrstat = 1002
    RETURN
  ENDIF

  IF ( fac_plcov < 0.0_ireals .OR. rmin_plcov > rmax_plcov ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** fac_plcov, rmin_plcov, rmax_plcov *** ',  &
            fac_plcov, rmin_plcov, rmax_plcov
    ierrstat = 1002
    RETURN
  ENDIF

  IF ( rmax_plcov > 1.0_ireals ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** rmax_plcov >1.0 is not allowed *** ', rmax_plcov
    ierrstat = 1002
    RETURN
  ENDIF

  IF ( fac_rootdp < 0. .OR. rmin_rootdp > rmax_rootdp ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** fac_rootdp, rmin_rootdp, rmax_rootdp *** ',  &
            fac_rootdp, rmin_rootdp, rmax_rootdp
    ierrstat = 1002
    RETURN
  ENDIF

  IF ( fac_lai < 0. .OR. rmin_lai > rmax_lai ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** fac_lai, rmin_lai, rmax_lai *** ',  &
            fac_lai, rmin_lai, rmax_lai
    ierrstat = 1002
    RETURN
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = iepsmem
    intbuf  ( 2) = iepstot
    intbuf  ( 3) = iepstyp

    realbuf  ( 1) = fac_plcov
    realbuf  ( 2) = rmin_plcov
    realbuf  ( 3) = rmax_plcov
    realbuf  ( 4) = fac_rootdp
    realbuf  ( 5) = rmin_rootdp
    realbuf  ( 6) = rmax_rootdp
    realbuf  ( 7) = fac_lai
    realbuf  ( 8) = rmin_lai
    realbuf  ( 9) = rmax_lai
  ENDIF

  CALL distribute_values (intbuf,  3, 0, imp_integers,  icomm_world, ierr)
  CALL distribute_values (realbuf, 9, 0, imp_reals,  icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    iepsmem  = intbuf  ( 1)
    iepstot  = intbuf  ( 2)
    iepstyp  = intbuf  ( 3)

    fac_plcov   = realbuf ( 1)
    rmin_plcov  = realbuf ( 2)
    rmax_plcov  = realbuf ( 3)
    fac_rootdp  = realbuf ( 4)
    rmin_rootdp = realbuf ( 5)
    rmax_rootdp = realbuf ( 6)
    fac_lai     = realbuf ( 7)
    rmin_lai    = realbuf ( 8)
    rmax_lai    = realbuf ( 9)
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A23)') '0     NAMELIST:  epsctl'
  WRITE (nuspecif, '(A23)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value', &
                                               'Default Value', 'Format'
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12  ,T59,A3)')                      &
                               'iepsmem',iepsmem, iepsmem_d  ,' I '
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12  ,T59,A3)')                      &
                               'iepstot',iepstot, iepstot_d  ,' I '
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12  ,T59,A3)')                      &
                               'iepstyp',iepstyp, iepstyp_d  ,' I '
  WRITE (nuspecif, '(T8,A,T21,F12.2,T40,F12.2  ,T59,A3)')                  &
                               'fac_plcov',fac_plcov, fac_plcov_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,F12.2,T40,F12.2  ,T59,A3)')                  &
                            'rmin_plcov',rmin_plcov, rmin_plcov_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,F12.2,T40,F12.2  ,T59,A3)')                  &
                            'rmax_plcov',rmax_plcov, rmax_plcov_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,F12.3,T40,F12.3  ,T59,A3)')                  &
                            'fac_rootdp',fac_rootdp, fac_rootdp_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,F12.3,T40,F12.3  ,T59,A3)')                  &
                         'rmin_rootdp',rmin_rootdp, rmin_rootdp_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,F12.3,T40,F12.3  ,T59,A3)')                  &
                         'rmax_rootdp',rmax_rootdp, rmax_rootdp_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,F12.2,T40,F12.2  ,T59,A3)')                  &
                                     'fac_lai',fac_lai, fac_lai_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,F12.2,T40,F12.2  ,T59,A3)')                  &
                                  'rmin_lai',rmin_lai, rmin_lai_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,F12.2,T40,F12.2  ,T59,A3)')                  &
                                  'rmax_lai',rmax_lai, rmax_lai_d  ,' R '
  WRITE (nuspecif, '(A2)')  '  '

ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_epsctl

!==============================================================================

!------------------------------------------------------------------------------
! End of external procedure organize_eps
!------------------------------------------------------------------------------

END SUBROUTINE organize_eps
