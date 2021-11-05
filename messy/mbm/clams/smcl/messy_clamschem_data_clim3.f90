MODULE MESSY_CLAMSCHEM_DATA_CLIM3

  ! update notes: 04/2019 update Cl+CH4, OH+CO to JPL2015 (jug)
  USE messy_clamschem_defs_mod,       ONLY: chch_t, ratb_t, ratj_t, ratt_t, rath_t
  USE messy_clamschem_asad_mod_clams, ONLY: jpspec, jpbk, jptk, jphk,  &
                                            jppj, chemdata_type

  IMPLICIT NONE

  TYPE(CHCH_T), DIMENSION(:), ALLOCATABLE :: chch_defs_clim3
  TYPE(RATB_T), DIMENSION(:), ALLOCATABLE :: ratb_defs_clim3
  TYPE(RATT_T), DIMENSION(:), ALLOCATABLE :: ratt_defs_clim3
  TYPE(RATJ_T), DIMENSION(:), ALLOCATABLE :: ratj_defs_clim3
  TYPE(RATH_T), DIMENSION(:), ALLOCATABLE :: rath_defs_clim3
  
CONTAINS

SUBROUTINE CLAMS_CHEM_INIT_CLIM3()

  IMPLICIT NONE
  
  INTEGER            :: ierr
  CHARACTER (LEN=72) :: cmessage 

  ALLOCATE(chch_defs_clim3(jpspec))
  ALLOCATE(ratb_defs_clim3(jpbk))
  ALLOCATE(ratj_defs_clim3(jppj))
  ALLOCATE(ratt_defs_clim3(jptk))
  ALLOCATE(rath_defs_clim3(jphk))
    
! TYPE CHCH_T
!   INTEGER           :: idum      ! dummy integer
!   CHARACTER(LEN=10) :: speci     ! species name
!   INTEGER           :: nodd      ! no of odd atoms
!   CHARACTER(LEN=10) :: ctype     ! Species type
!   CHARACTER(LEN=10) :: family    ! Family
!   INTEGER           :: switch1   ! 1 if dry deposits
!   INTEGER           :: switch2   ! 1 if wet deposits
!   INTEGER           :: switch3   ! > 0 if an emitted species
! ENDTYPE CHCH_T
chch_defs_clim3=(/ &
chch_t(  1,'CH4       ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t(  2,'CO        ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t(  3,'CO2       ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t(  4,'O3        ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t(  5,'N2O       ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t(  6,'F11       ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t(  7,'F12       ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t(  8,'H2O       ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t(  9,'HCl       ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t( 10,'HCN       ',  1,'TR        ','          ',  0,  0,  0),  &
chch_t( 11,'O2        ',  1,'CT        ','          ',  0,  0,  0),  &
chch_t( 12,'OH        ',  1,'CT        ','          ',  0,  0,  0),  &
chch_t( 13,'O(1D)     ',  1,'CT        ','          ',  0,  0,  0),  &
chch_t( 14,'Cl        ',  1,'CT        ','          ',  0,  0,  0),  &
chch_t( 15,'HO2       ',  1,'CT        ','          ',  0,  0,  0)   &
/)

! Describes the bimolecular reaction rates
! TYPE RATB_T
!   CHARACTER(LEN=10) :: react1    ! reactant name
!   CHARACTER(LEN=10) :: react2    ! reactant name
!   CHARACTER(LEN=10) :: prod1     ! product1 name
!   CHARACTER(LEN=10) :: prod2     ! product2 name
!   CHARACTER(LEN=10) :: prod3     ! product3 name
!   CHARACTER(LEN=10) :: prod4     ! product4 name
!   REAL              :: K0        ! rate coeff param K0
!   REAL              :: alpha     ! rate coeff param alpha
!   REAL              :: beta      ! rate coeff param beta
!   REAL              :: pyield1   ! product1 fractional yield
!   REAL              :: pyield2   ! product2 fractional yield
!   REAL              :: pyield3   ! product3 fractional yield
!   REAL              :: pyield4   ! product4 fractional yield
! ENDTYPE RATB_T
ratb_defs_clim3=(/              & 
ratb_t('OH        ','CH4       ','H2O       ','H2O       ','CO        ',   &
'          ',2.45E-12,  0.00,    1775.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ','CH4       ','H2O       ','H2O       ','CO        ',   &
'          ',1.75E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Cl        ','CH4       ','H2O       ','H2O       ','CO        ',   &
'          ',7.10E-12,  0.00,   1270.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('N2O       ','O(1D)     ','N2        ','O2        ','          ',   &
'          ',1.19E-10,  0.00,    -20.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','O3        ','HO2       ','O2        ','          ',   &
'          ',0.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','O3        ','OH        ','O2        ','O2        ',   &
'          ',2.00E-14,  0.00,    490.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HCN       ','OH        ','          ','          ','          ',   &
'          ',1.20E-13,  0.00,    400.00, 0.000, 0.000, 0.000, 0.000)    &
  /)




! describes photolytic reactions
! TYPE RATJ_T
!   CHARACTER(LEN=10) :: react1    ! reactant name
!   CHARACTER(LEN=10) :: react2    ! reactant name
!   CHARACTER(LEN=10) :: prod1     ! product1 name
!   CHARACTER(LEN=10) :: prod2     ! product2 name
!   CHARACTER(LEN=10) :: prod3     ! product3 name
!   CHARACTER(LEN=10) :: prod4     ! product4 name
!   REAL              :: pyield1   ! product yield
!   REAL              :: pyield2   ! product yield
!   REAL              :: pyield3   ! product yield
!   REAL              :: pyield4   ! product yield
!   REAL              :: jfacta    ! quantum yield
!   CHARACTER(LEN=10) :: fname     ! file name/label
! ENDTYPE RATJ_T
ratj_defs_clim3 = (/                  &
ratj_t('O2       ','PHOTON    ','O3         ','O3        ','          ', &
     '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('N2O      ','PHOTON    ','O(1D)      ','N2        ','          ', &
     '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('F11      ','PHOTON    ','Cl         ','Cl2       ','          ', &
     '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('F12      ','PHOTON    ','Cl         ','Cl        ','          ',&
     '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       ')      &
  /)
   
! Describes termolecular reactions
! TYPE RATT_T
!   CHARACTER(LEN=10) :: react1    ! reactant name
!   CHARACTER(LEN=10) :: react2    ! reactant name
!   CHARACTER(LEN=10) :: prod1     ! product1 name
!   CHARACTER(LEN=10) :: prod2     ! product2 name
!   REAL              :: F         ! rate coeff param F
!   REAL              :: K1        ! rate coeff param K1
!   REAL              :: alpha1    ! rate coeff param alpha1
!   REAL              :: beta1     ! rate coeff param beta1
!   REAL              :: K2        ! rate coeff param K2
!   REAL              :: alpha2    ! rate coeff param alpha2
!   REAL              :: beta2     ! rate coeff param beta2
!   REAL              :: pyield1   ! product1 fractional yield
!   REAL              :: pyield2   ! product2 fractional yield
! ENDTYPE RATT_T
ratt_defs_clim3=(/ &
ratt_t('O(3P)     ','O2        ','O3        ','m         ', &
     0.0, 6.00E-34, -2.40,     0.0, 0.00E+00,  0.00,       0.0, 0.000, 0.000) , &            
ratt_t('OH        ','CO        ','CO2       ','m         ', &
     0.6, 5.90E-33, -1.00,     0.0, 1.10E-12,  1.30,       0.0, 0.000, 0.000)   &            
  /)
  
! Describes heterogenous reactions
! TYPE RATH_T
!   CHARACTER(LEN=10) :: react1    ! reactant1 name
!   CHARACTER(LEN=10) :: react2    ! reactant2 name
!   CHARACTER(LEN=10) :: prod1     ! product1 name
!   CHARACTER(LEN=10) :: prod2     ! product2 name
!   CHARACTER(LEN=10) :: prod3     ! product3 name
!   CHARACTER(LEN=10) :: prod4     ! product4 name
!   REAL              :: pyield1   ! product yield
!   REAL              :: pyield2   ! product yield
!   REAL              :: pyield3   ! product yield
!   REAL              :: pyield4   ! product yield
! ENDTYPE RATH_T
rath_defs_clim3=(/ &
rath_t('ClONO2    ','H2O       ','HOCl      ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000), &
rath_t('ClONO2    ','HCl       ','Cl2       ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000), &
rath_t('HOCl      ','HCl       ','Cl2       ','H2O       ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000), &
rath_t('N2O5      ','H2O       ','HONO2     ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000), &
rath_t('N2O5      ','HCl       ','ClNO2     ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000),  &
rath_t('ClONO2    ','HBr       ','BrCl      ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000),  &
rath_t('BrONO2    ','HCl       ','BrCl      ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000),  &
rath_t('HOCl      ','HBr       ','BrCl      ','H2O       ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000),  &
rath_t('HOBr      ','HCl       ','BrCl      ','H2O       ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000),  &
rath_t('HOBr      ','HBr       ','Br2       ','H2O       ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000),  &
rath_t('BrONO2    ','H2O       ','HOBr      ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000)  &
  /)


  END SUBROUTINE CLAMS_CHEM_INIT_CLIM3


  SUBROUTINE CLAMS_CHEM_CLEAN_CLIM3

    implicit none

    DEALLOCATE(chch_defs_clim3)
    DEALLOCATE(ratb_defs_clim3)
    DEALLOCATE(ratj_defs_clim3)
    DEALLOCATE(ratt_defs_clim3)

  END SUBROUTINE CLAMS_CHEM_CLEAN_CLIM3


END MODULE MESSY_CLAMSCHEM_DATA_CLIM3
