MODULE MESSY_CLAMSCHEM_DEFS_MOD

   USE messy_clams_global, ONLY: prec

IMPLICIT NONE

! Describes the chemical species used in the model
TYPE CHCH_T
  INTEGER           :: idum      ! dummy integer
  CHARACTER(LEN=10) :: speci     ! species name
  INTEGER           :: nodd      ! no of odd atoms
  CHARACTER(LEN=10) :: ctype     ! Species type
  CHARACTER(LEN=10) :: family    ! Family
  INTEGER           :: switch1   ! 1 if dry deposits
  INTEGER           :: switch2   ! 1 if wet deposits
  INTEGER           :: switch3   ! > 0 if an emitted species
ENDTYPE CHCH_T
PUBLIC CHCH_T

! Describes the bimolecular reaction rates
TYPE RATB_T
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL(PREC)        :: K0        ! rate coeff param K0
  REAL(PREC)        :: alpha     ! rate coeff param alpha
  REAL(PREC)        :: beta      ! rate coeff param beta
  REAL(PREC)        :: pyield1   ! product1 fractional yield
  REAL(PREC)        :: pyield2   ! product2 fractional yield
  REAL(PREC)        :: pyield3   ! product3 fractional yield
  REAL(PREC)        :: pyield4   ! product4 fractional yield
ENDTYPE RATB_T
PUBLIC RATB_T

! Describes heterogenous reactions
TYPE RATH_T
  CHARACTER(LEN=10) :: react1    ! reactant1 name
  CHARACTER(LEN=10) :: react2    ! reactant2 name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL(PREC)        :: pyield1   ! product yield
  REAL(PREC)        :: pyield2   ! product yield
  REAL(PREC)        :: pyield3   ! product yield
  REAL(PREC)        :: pyield4   ! product yield
ENDTYPE RATH_T
PUBLIC RATH_T

! describes photolytic reactions
TYPE RATJ_T
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL(PREC)        :: pyield1   ! product yield
  REAL(PREC)        :: pyield2   ! product yield
  REAL(PREC)        :: pyield3   ! product yield
  REAL(PREC)        :: pyield4   ! product yield
  REAL(PREC)        :: jfacta    ! quantum yield
  CHARACTER(LEN=10) :: fname     ! file name/label
ENDTYPE RATJ_T
PUBLIC RATJ_T

! Describes termolecular reactions
TYPE RATT_T
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  REAL(PREC)        :: F         ! rate coeff param F
  REAL(PREC)        :: K1        ! rate coeff param K1
  REAL(PREC)        :: alpha1    ! rate coeff param alpha1
  REAL(PREC)        :: beta1     ! rate coeff param beta1
  REAL(PREC)        :: K2        ! rate coeff param K2
  REAL(PREC)        :: alpha2    ! rate coeff param alpha2
  REAL(PREC)        :: beta2     ! rate coeff param beta2
  REAL(PREC)        :: pyield1   ! product1 fractional yield
  REAL(PREC)        :: pyield2   ! product2 fractional yield
ENDTYPE RATT_T
PUBLIC RATT_T

TYPE(CHCH_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: chch_defs
TYPE(RATB_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: ratb_defs
TYPE(RATH_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: rath_defs
TYPE(RATJ_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: ratj_defs
TYPE(RATT_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: ratt_defs

END MODULE MESSY_CLAMSCHEM_DEFS_MOD
