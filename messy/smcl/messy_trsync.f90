! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL TRSYNC
! (SYNCHRONIZATION OF ISOTOPOLOGICAL TRACERS IN H2OISO AND CH4/MECCA-TAG)
!
! Author : Franziska Frank, DLR-IPA, April 2016
!          
!
! References:
!
! **********************************************************************
!>
!> \mainpage Introduction
!> 
!> This MESSy submodel converts the tracer from kg/kg_moistair to 
!> mol/mol_dryair and back again
!>
!> \authors Franziska Frank, DLR-IPA
!>  - April 2016: subroutines convert_to_kgkg and convert_to_molmol
!> 

! **********************************************************************
MODULE messy_trsync
  ! **********************************************************************

  ! ----------- >
  USE messy_main_constants_mem, ONLY: DP, M_air

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'trsync'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'

  ! ### add your own public subroutines here
  PUBLIC  :: convert_unit

  ! ### add your own private subroutines here
  PRIVATE :: convert_to_kgkg
  PRIVATE :: convert_to_molmol
  PRIVATE :: convert_to_kgkg_te
  PRIVATE :: convert_to_molmol_te

CONTAINS

  ! =========================================================================
  ELEMENTAL SUBROUTINE convert_unit(traten, case, type, molarmass, spechum, spechum_te, tracer)

    ! ------------------------------------------------------------------
    ! Conversion from one unit to the other 
    ! type specifies if a tracer or a tendency will be converted
    !
    ! case = 1 : kg/kg_moistair -> mol/mol_dryair
    !      = 2 : mol/mol_dryair -> kg/kg_moistair
    !      > 3 : not implemented
    ! type = 1 : tracer conversion
    !      = 2 : tendency conversion
    !      > 3 : not implemented
    ! ------------------------------------------------------------------

  !  USE  messy_main_blather,    ONLY: warning

    IMPLICIT NONE 

    ! I/O
    REAL(dp), INTENT(INOUT)        :: traten        ! tracer or tendency
    INTEGER, INTENT(IN)            :: case
    INTEGER, INTENT(IN)            :: type
    REAL(dp), INTENT(IN)           :: molarmass
    REAL(dp), INTENT(IN)           :: spechum
    REAL(dp), INTENT(IN), OPTIONAL :: spechum_te    ! optional but necessary for
    REAL(dp), INTENT(IN), OPTIONAL :: tracer        ! tendency conversion (i.e. 
                                                    ! type = 2

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'convert_unit'

    IF ((type == 2) &
         .AND. (( .NOT. PRESENT(spechum_te) ) &
         .OR. ( .NOT. PRESENT(tracer)))) THEN
   !    CALL warning("The function convert_unit cannot be used like this!", substr)
       RETURN
    END IF

    SELECT CASE (type)
    ! Convert tracer
    CASE(1)
       SELECT CASE (case)
       CASE(1)
          CALL convert_to_molmol(traten, molarmass, spechum)
       CASE(2)
          CALL convert_to_kgkg(traten, molarmass, spechum)
       END SELECT
    ! Convert tendency
    CASE(2)
       SELECT CASE (case)
       CASE(1)
!          CALL convert_to_molmol(traten, molarmass, spechum)
          CALL convert_to_molmol_te(traten, tracer, molarmass, spechum, spechum_te)
       CASE(2)
!          CALL convert_to_kgkg(traten, molarmass, spechum)
          CALL convert_to_kgkg_te(traten, tracer, molarmass, spechum, spechum_te)
       END SELECT

    END SELECT
    
  END SUBROUTINE convert_unit
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE convert_to_kgkg(tr_a, molarmass, spechum)

    ! ------------------------------------------------------------------
    ! Conversion from mol/mol_dryair to kg/kg_moistair
    ! ------------------------------------------------------------------

    ! I/O
    REAL(dp), INTENT(INOUT) :: tr_a          ! [mol/mol_dryair]
    REAL(dp), INTENT(IN)    :: molarmass
    REAL(dp), INTENT(IN)    :: spechum

    ! LOCAL
    REAL(DP)                :: scvmr 
    REAL(dp)                :: mue
    
    ! Molarmass ratio
    scvmr = M_air/molarmass 
    
    ! Conversion factor
    mue = scvmr * (1.0_dp / (1.0_dp - spechum))

    ! Convert tracer
    tr_a = tr_a / mue

  END SUBROUTINE convert_to_kgkg
  ! =========================================================================
 
 ! =========================================================================
  ELEMENTAL SUBROUTINE convert_to_kgkg_te(tr_a_te, tr_a, molarmass, spechum, spechum_te)

    ! ------------------------------------------------------------------
    ! Conversion from mol/mol_dryair to kg/kg_moistair of a tendency
    ! ------------------------------------------------------------------

    ! I/O
    REAL(dp), INTENT(INOUT) :: tr_a_te       ! [mol/mol_dryair]
    REAL(dp), INTENT(IN)    :: tr_a          ! [mol/mol_dryair]
    REAL(dp), INTENT(IN)     :: molarmass
    REAL(dp), INTENT(IN)    :: spechum
    REAL(dp), INTENT(IN)    :: spechum_te

    ! LOCAL
    REAL(DP)                :: scvmr
    REAL(dp)                :: mue
    REAL(dp)                :: eta
    
    ! Molarmass ratio
    scvmr = M_air/molarmass 
   
    ! Conversion factors
    mue = scvmr * (1.0_dp / (1.0_dp - spechum))
    eta = spechum_te / scvmr

    ! Convert tendency
    tr_a_te = tr_a_te / mue  - tr_a * eta

  END SUBROUTINE convert_to_kgkg_te
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE convert_to_molmol(tr_b, molarmass, spechum)

    ! ------------------------------------------------------------------
    ! Conversion from kg/kg_moistair to mol/mol_dryair
    ! ------------------------------------------------------------------

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(INOUT) :: tr_b      ! [kg/kg_moistair]
    REAL(dp), INTENT(IN)    :: molarmass
    REAL(dp), INTENT(IN)    :: spechum

    ! LOCAL
    REAL(DP)                :: scvmr
    REAL(dp)                :: mue
    
    ! Molarmass ratio
    scvmr = M_air/molarmass 
   
    ! Conversion factor
    mue = scvmr * (1.0_dp / (1.0_dp - spechum))

    ! Convert tracer
    tr_b = tr_b * mue

  END SUBROUTINE convert_to_molmol
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE convert_to_molmol_te(tr_b_te, tr_b, molarmass, spechum, spechum_te)

    ! ------------------------------------------------------------------
    ! Conversion from kg/kg_moistair to mol/mol_dryair of a tendency
    ! ------------------------------------------------------------------

!    USE messy_main_constants_mem,    ONLY: M_air, M_H2O

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(INOUT) :: tr_b_te ! [kg/kg_moistair]
    REAL(dp), INTENT(IN)    :: tr_b    ! [kg/kg_moistair]
    REAL(dp), INTENT(IN)    :: molarmass
    REAL(dp), INTENT(IN)    :: spechum
    REAL(dp), INTENT(IN)    :: spechum_te

    ! LOCAL
    REAL(DP)                :: scvmr
    REAL(dp)                :: mue
    REAL(dp)                :: eta    

    ! Molarmass ratio
    scvmr = M_air/molarmass

    ! Conversion factor
    mue = scvmr * (1.0_dp / (1.0_dp - spechum))
    eta = spechum_te / scvmr

    ! Convert tendency
    tr_b_te = tr_b_te * mue + tr_b * mue**2 * eta

  END SUBROUTINE convert_to_molmol_te
  ! =========================================================================

  ! **********************************************************************
END MODULE messy_trsync
! **********************************************************************
