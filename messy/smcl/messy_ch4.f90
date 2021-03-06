! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL CH4
! (SIMPLIFIED METHANE CHEMISTRY WITH PRESCRIBED EDUCTS)
!
! Author : Patrick Joeckel, DLR-IPA, May 2012
!          (based on original coding in submodel BUFLY by Martin Weismann)
!
! References:
!
! **********************************************************************
!>
!> \mainpage Introduction
!> 
!> This MESSy submodel calculates the methane rections with
!> OH, Cl, and O<SUP>1</SUP>D. The educts need to be prescribed either
!> on-line by running a full chemistry module such as MECCA, or imported 
!> from external files (via IMPORT).
!>
!> \authors Patrick J?ckel, DLR-IPA
!>  - May 2012: original code based on submodel BUFLY by Martin Weismann
!>  - June 2013: Lagrangian extension
!>  - July 2014: methane age classes
!> \authors Franziska Frank, DLR-IPA
!>  - October 2014: methane isotopologues and kinetic fractionation
!> 

! **********************************************************************
MODULE messy_ch4
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'ch4'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.3'
 
  ! CTRL-NAMELIST PARAMETERS
  ! ISO++
  REAL(DP), DIMENSION(2), PUBLIC :: KIE_CH4_13C_OH  = (/1, 0/)
  REAL(DP), DIMENSION(2), PUBLIC :: KIE_CH4_13C_O1D = (/1, 0/)
  REAL(DP), DIMENSION(2), PUBLIC :: KIE_CH4_13C_CL  = (/1, 0/)
  REAL(DP), DIMENSION(2), PUBLIC :: KIE_CH4_13C_jval  = (/1, 0/)
  REAL(DP), DIMENSION(2), PUBLIC :: KIE_CH4_D1_OH   = (/1, 0/)
  REAL(DP), DIMENSION(2), PUBLIC :: KIE_CH4_D1_O1D  = (/1, 0/)
  REAL(DP), DIMENSION(2), PUBLIC :: KIE_CH4_D1_CL   = (/1, 0/)
  REAL(DP), DIMENSION(2), PUBLIC :: KIE_CH4_D1_jval   = (/1, 0/)

  INTEGER, PARAMETER, PUBLIC :: iso_12C = 1 
  INTEGER, PARAMETER, PUBLIC :: iso_13C = 2
  INTEGER, PARAMETER, PUBLIC :: iso_D0  = 3
  INTEGER, PARAMETER, PUBLIC :: iso_D1  = 4
  INTEGER, DIMENSION(4), PUBLIC :: iso_id = &
       (/ iso_12C, iso_13C, iso_D0, iso_D1 /)
  CHARACTER(LEN=3), DIMENSION(4), PARAMETER, PUBLIC :: iso_name = &
       (/'12C','13C','D0 ','D1 '/)
  ! ISO--

  ! PUBLIC SUBROUTINES (to be called from messy_ch4_e5.f90)
  PUBLIC :: ch4_read_nml_ctrl
  ! ### add your own public subroutines here
  PUBLIC :: ch4_integrate

  PUBLIC :: sca_tend
  PUBLIC :: adj_tend  

  PRIVATE :: calc_KIE

CONTAINS

  ! =========================================================================
  ELEMENTAL SUBROUTINE ch4_integrate(CH4_te, CH4, OH, O1D, Cl, j_CH4 &
       , temp, press, sphum, iso_id)

    USE messy_main_constants_mem,    ONLY: N_A, R_gas, M_air, M_H2O

    IMPLICIT NONE
    INTRINSIC :: exp, log

    ! I/O
    REAL(dp), INTENT(OUT) :: CH4_te   !< methane tendency [mol/mol/s]
    REAL(dp), INTENT(IN)  :: CH4      !< methane [mol/mol]
    REAL(dp), INTENT(IN)  :: OH       !< hydroxyl radical [mol/mol]
    REAL(dp), INTENT(IN)  :: O1D      !< excited oxygen [mol/mol]
    REAL(dp), INTENT(IN)  :: Cl       !< chlorine [mol/mol]
    REAL(dp), INTENT(IN)  :: j_CH4    !< photolysis rate [1/s]
    !
    REAL(dp), INTENT(IN)  :: temp     !< temperature [K]
    REAL(dp), INTENT(IN)  :: press    !< pressure [Pa]
    REAL(dp), INTENT(IN)  :: sphum    !< specific humidity [kg/kg]
    ! ISO++
    INTEGER, INTENT(IN), OPTIONAL :: iso_id !< id of isotopologue
    ! ISO--

    ! LOCAL
    REAL(dp), PARAMETER   :: vtmpc1 = M_air / M_H2O - 1._dp
    ! ... reaction coefficients (1st order reactions)
    REAL(dp)              :: k_O1D               ! [cm^3/s]
    REAL(dp)              :: k_OH                ! [cm^3/s]
    REAL(dp)              :: k_Cl                ! [cm^3/s]
    REAL(dp)              :: cair     ! concentration of air [1/cm^3]
    REAL(dp)              :: k_jval              ! [1/s]
    ! ... reaction rates 
    REAL(dp)              :: r_O1D               ! [1/s]
    REAL(dp)              :: r_OH                ! [1/s]
    REAL(dp)              :: r_Cl                ! [1/s]
    ! ISO++
    ! kinetic isotope effect
    REAL(dp)              :: KIE                 
    ! ISO--

    cair = (N_A/1.E6_dp) * press / (R_gas*temp *(1.0_dp + vtmpc1*sphum))

    ! standard reaction rates (apply if no, or no valid iso_id is provided)
    k_OH = 1.85E-20_dp * exp(2.82_dp*log(temp) - 987._dp/temp)
    k_Cl = 6.6E-12_dp  * exp(-1240._dp/temp)
    k_O1D = 1.75E-10_dp
    k_jval = j_CH4
    
    ! ISO++
    IF (PRESENT(iso_id)) THEN
       SELECT CASE(iso_id)
          CASE(iso_12C)
             ! standard reaction rates
          CASE(iso_13C)
             ! For isotopologues: KIE = k_X / k_X_iso
             !                  =>k_X_iso = KIE^-1 * k_X
             !
             CALL calc_KIE(KIE_CH4_13C_OH, temp, KIE)
             k_OH = 1.0_dp/KIE * k_OH
             CALL calc_KIE(KIE_CH4_13C_Cl, temp, KIE)
             k_Cl = 1.0_dp/KIE * k_Cl
             CALL calc_KIE(KIE_CH4_13C_O1D, temp, KIE)
             k_O1D = 1.0_dp/KIE * k_O1D
             CALL calc_KIE(KIE_CH4_13C_jval, temp, KIE)
             k_jval = 1.0_dp / KIE * k_jval
          CASE(iso_D0)
             ! standard reaction rates
          CASE(iso_D1)
             ! For isotopologues: KIE = k_X / k_X_iso
             !                  =>k_X_iso = KIE^-1 * k_X
             CALL calc_KIE(KIE_CH4_D1_OH, temp, KIE)
             k_OH = 1.0_dp/KIE * k_OH
             CALL calc_KIE(KIE_CH4_D1_Cl, temp, KIE)
             k_Cl = 1.0_dp/KIE * k_Cl
             CALL calc_KIE(KIE_CH4_D1_O1D, temp, KIE)
             k_O1D = 1.0_dp/KIE * k_O1D  
             CALL calc_KIE(KIE_CH4_D1_jval, temp, KIE)
             k_jval = 1.0_dp / KIE * k_jval
          END SELECT
    END IF
    ! ISO--

    r_OH  = k_OH  * (cair * OH)   ! 1/s
    r_Cl  = k_Cl  * (cair * Cl)   ! 1/s
    r_O1D = k_O1D * (cair * O1D)  ! 1/s

    ! linearised solution
    CH4_te = -1._dp * CH4 * (r_OH + r_Cl + r_O1D + k_jval) 

  END SUBROUTINE ch4_integrate
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE sca_tend(m, mte, s, ste, dt, a)

    ! CH4 MODULE ROUTINE (SMCL)
    !
    ! LINEARIZE MULTI-COMPONENT 'TAGGED' TRACER TENDENCIES:
    ! TRACER TENDENCIES ARE ADJUSTET TO FORCE:
    ! m + mte*dt = (f1 + t1*dt) + (f2 + t2*dt) + ...
    ! - m is the 'master' tracer, f1, f2 the 'fractional' tracers,
    ! - s = f1 + f2 + ... is the sum of the fractional tracers,
    ! - mte, t1, and t2 are the respective tendencies, 
    ! - ste = t1 + t2 + ...
    !
    ! This can be used to correct for non-linearities resulting from
    ! tracer gradient dependent processes, such as, e.g., the advection.
    !
    ! Author: Patrick Joeckel, DLR, Jun 2014

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(IN)  :: m    ! master tracer
    REAL(DP), INTENT(IN)  :: mte  ! tendency of sum tracer
    REAL(DP), INTENT(IN)  :: s    ! sum of fractional tracers
    REAL(DP), INTENT(IN)  :: ste  ! sum of fract. tendencies
    REAL(DP), INTENT(IN)  :: dt   ! time step length
    REAL(DP), INTENT(OUT) :: a    ! resulting weight (correction factor)

    ! LOCAL
    REAL(DP),         PARAMETER :: EPSILON = TINY(0.0_DP)    ! zero
    REAL(DP) :: sum

    ! m + mte*dt =!= (f1 + t1*dt) + (f2 + t2*dt) + ...
    ! a := (m+mte*dt) / ((f1+t1*dt)+(f2+t2*dt)+ ...)
    ! (f1+t1*dt)' := a*(f1+t1*dt)
    ! (f2+t2*dt)' := a*(f2+t2*dt)
    ! ...
    ! => (f1+t1*dt)' + (f2+t2*dt)' + ...' = (s+t*dt)   O.K.
    !
    ! do not change tracers, rather adjust tendencies
    ! cond.: (f2' == f2) AND (f1' == f1) AND ...
    ! => f1 + t1'*dt = a*f1 + a*t1*dt ; f2 + t2'*dt = a*f2 + a*t2*dt ; ...
    ! => t1' = f1*(a-1)/dt + a*t1  ; t2' = f2*(a-1)/dt + a*t2; ...

    sum = s + ste*dt
    IF (ABS(sum) >= EPSILON) THEN
       a = (m + mte*dt)/sum
    ELSE
       a = 1.0_DP
    ENDIF

  END SUBROUTINE sca_tend
! =========================================================================

! =========================================================================
  ELEMENTAL SUBROUTINE adj_tend(f, t, a, dt, tadj)
    
    ! I/O
    REAL(DP), INTENT(IN)    :: f     ! fractional tracer
    REAL(DP), INTENT(IN)    :: t     ! tendency of fractional tracer
    REAL(DP), INTENT(IN)    :: a     ! correction (scaling) factor
    REAL(DP), INTENT(IN)    :: dt    ! time step length
    REAL(DP), INTENT(OUT)   :: tadj  ! additional tendency for adjustment 

    ! LOCAL
    REAL(DP) :: ta ! adjusted tendency

    ta = t*a + f*(a-1.0_dp)/dt
    tadj = ta - t                ! ta = t + tadj

  END SUBROUTINE adj_tend
! =========================================================================

! ISO++
! =========================================================================
  PURE SUBROUTINE calc_KIE(KIE_AB_val, temp_t, KIE_t)

    ! I/O
    REAL(DP), DIMENSION(2), INTENT(IN)  :: KIE_AB_val
    REAL(DP),               INTENT(IN)  :: temp_t
    REAL(DP),               INTENT(OUT) :: KIE_t

    ! KIE(T) = A * exp(B/T) with T := temperature
    KIE_t = KIE_AB_val(1) * exp(KIE_AB_val(2) / temp_t)

  END SUBROUTINE calc_KIE
! =========================================================================
! ISO--

  ! =========================================================================
  SUBROUTINE ch4_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ KIE_CH4_13C_OH, KIE_CH4_13C_O1D, KIE_CH4_13C_CL, &
                    KIE_CH4_13C_jval, KIE_CH4_D1_OH, KIE_CH4_D1_O1D, &
                    KIE_CH4_D1_CL, KIE_CH4_D1_jval

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='ch4_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    KIE_CH4_13C_OH  = (/1, 0/)
    KIE_CH4_13C_O1D = (/1, 0/)
    KIE_CH4_13C_CL  = (/1, 0/)
    KIE_CH4_13C_jval= (/1, 0/)
    KIE_CH4_D1_OH   = (/1, 0/)
    KIE_CH4_D1_O1D  = (/1, 0/)
    KIE_CH4_D1_CL   = (/1, 0/)
    KIE_CH4_D1_jval = (/1, 0/)

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR
    
  END SUBROUTINE ch4_read_nml_ctrl
  ! =========================================================================

! **********************************************************************
END MODULE messy_ch4
! **********************************************************************

