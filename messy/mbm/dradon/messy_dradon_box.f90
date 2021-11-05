! ***************************************************************************
MODULE messy_dradon_box
! ***************************************************************************

  USE messy_dradon, ONLY: DP, NI, modstr

  IMPLICIT NONE

  ! CPL_BOX NAMELIST
  REAL(DP) :: delta_time = 1800   ! time step in seconds
  INTEGER  :: n_steps    = 1000   ! number of time steps to integrate
  INTEGER  :: n_out      = 5      ! output step interval
  REAL(DP) :: r_222Rn_0  = 0.0_DP ! start value for 222Rn

  ! I/O
  INTEGER, PARAMETER :: unit = 21

  ! VARIABLES
  ! CASE 1: source + decay WITHOUT operator splitting
  REAL(DP) :: Rn_222    = 0.0_DP
  REAL(DP) :: Rn_222_te = 0.0_DP
  !
  ! CASE 2: decay chain: Rn source and decay WITH operator splitting
  REAL(DP), DIMENSION(NI) :: Rn    ! isotopes at t=0
  REAL(DP), DIMENSION(NI) :: Rn_te ! isotopes tendencies

  ! SUBROUTINES
  PUBLIC :: initialize
  PUBLIC :: integrate
  PUBLIC :: finalize
  !PRIVATE :: open_output
  !PRIVATE :: close_output
  !PRIVATE :: dradon_read_nml_cpl

CONTAINS

  ! ================================================================
  ! PUBLIC ROUTINES
  ! ================================================================

  ! ----------------------------------------------------------------
  SUBROUTINE initialize

    USE messy_dradon, ONLY: dradon_read_nml_ctrl, init_radon_chain &
                          , R_Rn_cflux_land

    IMPLICIT NONE

    INTEGER :: status

    CALL dradon_read_nml_ctrl(status, unit)
    IF (status /= 0) STOP

    WRITE(*,*) '**********************************************'
    WRITE(*,*) '*** NOTE: Parameters'
    WRITE(*,*) '***       - I_Rn_flux_method'
    WRITE(*,*) '***       - R_Rn_cflux_ocean'
    WRITE(*,*) '***       ignored in BOX-model calculation!'
    WRITE(*,*) '*** 222Rn emission [atoms m^(-2) s^(-1)] is '
    WRITE(*,*) '***       R_Rn_cflux_land = ',R_Rn_cflux_land
    WRITE(*,*) '**********************************************'

    CALL dradon_read_nml_cpl(status, unit)
    IF (status /= 0) STOP

    CALL open_output

    CALL init_radon_chain

  END SUBROUTINE initialize
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  SUBROUTINE integrate

    USE messy_dradon, ONLY: int_radon_chain, int_radon, R_Rn_cflux_land

    IMPLICIT NONE

    INTEGER :: jstep

    ! INITIALIZE
    Rn_222 = r_222Rn_0
    Rn(:)  = 0.0_DP
    Rn(1)  = r_222Rn_0

    ! TIME LOOP
    time_loop: DO jstep = 0, n_steps

       ! CHECK OUTPUT
       IF (MOD(jstep, n_out) == 0) &
            WRITE(unit,'(7(e14.6,1x))') jstep*delta_time, Rn_222, Rn(:)
       
       ! CALCULATE NEW TENDENCIES AND UPDATE NEW VALUE
       ! CASE 1:
       Rn_222_te = 0.0_DP
       CALL int_radon(Rn_222_te, Rn_222, delta_time, R_Rn_cflux_land)
       Rn_222 = Rn_222 + delta_time * Rn_222_te

       ! CASE 2:
       Rn_te(:) = 0.0_DP
       ! SOURCE
       Rn_te(1) = R_Rn_cflux_land
       Rn(1) = Rn(1) + delta_time * Rn_te(1)
       ! CHAIN DECAY
       Rn_te(:) = 0.0_DP
       CALL int_radon_chain(Rn_te(:), Rn(:), delta_time)
       Rn(:) = Rn(:) + delta_time * Rn_te(:)

    END DO time_loop

  END SUBROUTINE integrate
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  SUBROUTINE finalize

    IMPLICIT NONE

    CALL close_output

  END SUBROUTINE finalize
  ! ----------------------------------------------------------------

  ! ================================================================
  ! PRIVATE ROUTINES
  ! ================================================================

  ! ----------------------------------------------------------------------
  SUBROUTINE dradon_read_nml_cpl(status, iou)

    ! DRADON MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5/ATTILA
    !
    ! Author: Patrick Joeckel, MPICH, Oct 2003

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_read_nml_cpl'

    NAMELIST /CPL_BOX/ delta_time, n_steps, n_out, r_222Rn_0

    ! LOCAL
    LOGICAL :: lex      ! file exists ?
    INTEGER :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL_BOX', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_BOX, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_BOX', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
    WRITE(*,*) 'TIME STEP            :', delta_time,' s'
    WRITE(*,*) 'NUMBER OF TIME STEPS :', n_steps
    WRITE(*,*) 'OUTPUT STEP INTERVAL :', n_out

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE dradon_read_nml_cpl
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE open_output

    IMPLICIT NONE

    OPEN(unit=unit, STATUS='unknown', file='dradon.dat')
    WRITE(unit,'(7(a14,1x))') &
         'time','222Rn','222Rn','218Po','214Pb','214Bi','210Pb'
    WRITE(unit,'(7(a14,1x))') &
         '--------------', &
         '--------------', &
         '--------------', &
         '--------------', &
         '--------------', &
         '--------------', &
         '--------------'

  END SUBROUTINE open_output
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE close_output

    IMPLICIT NONE

    CLOSE(unit)

  END SUBROUTINE close_output
  ! ----------------------------------------------------------------------

! ***************************************************************************
END MODULE messy_dradon_box
! ***************************************************************************
