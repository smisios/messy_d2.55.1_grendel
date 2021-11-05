MODULE MESSY_CLAMSCHEM_DATA_STANDARD
  !ju_jg_20190425
  ! update notes: 04/2019 update MeCl + OH, Cl+CH4,
  !               ClO+MeOO, HO2+NO2, NO2+NO3, ClO+ClO,
  !               BrO+NO2, OH+CO  to JPL2015 (jug)

  USE messy_clamschem_defs_mod,       ONLY: chch_t, ratb_t, rath_t, ratj_t, ratt_t
  USE messy_clamschem_asad_mod_clams, ONLY: jpspec, jpbk, jptk, jphk,  &
                                            jppj, chemdata_type

  IMPLICIT NONE

  TYPE(CHCH_T), DIMENSION(:), ALLOCATABLE :: chch_defs_std
  TYPE(RATB_T), DIMENSION(:), ALLOCATABLE :: ratb_defs_std
  TYPE(RATT_T), DIMENSION(:), ALLOCATABLE :: ratt_defs_std
  TYPE(RATJ_T), DIMENSION(:), ALLOCATABLE :: ratj_defs_std
  TYPE(RATH_T), DIMENSION(:), ALLOCATABLE :: rath_defs_std
  
CONTAINS

SUBROUTINE CLAMS_CHEM_INIT_STANDARD()

  IMPLICIT NONE
  
  INTEGER            :: ierr
  CHARACTER (LEN=72) :: cmessage 

  ALLOCATE(chch_defs_std(jpspec))
  ALLOCATE(ratb_defs_std(jpbk))
  ALLOCATE(ratj_defs_std(jppj))
  ALLOCATE(ratt_defs_std(jptk))
  ALLOCATE(rath_defs_std(jphk))
    
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
chch_defs_std=(/ &
chch_t( 1, 'O(1D)     ',  1, 'FM        ','Ox        ',  0,  0,  0),  &
chch_t( 2, 'O(3P)     ',  1, 'FM        ','Ox        ',  0,  0,  0),  &
chch_t( 3, 'O3        ',  1, 'FM        ','Ox        ',  0,  0,  0),  &
chch_t( 4, 'NO        ',  1, 'FM        ','NOx       ',  0,  0,  0),  &
chch_t( 5, 'NO3       ',  1, 'FM        ','NOx       ',  0,  0,  0),  &
chch_t( 6, 'NO2       ',  1, 'FM        ','NOx       ',  0,  0,  0),  &
chch_t( 7, 'N2O5      ',  2, 'TR        ','          ',  0,  0,  0),  &
chch_t( 8, 'HO2NO2    ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t( 9, 'HONO2     ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(10, 'OH        ',  1, 'SS        ','          ',  0,  0,  0),  &
chch_t(11, 'HO2       ',  1, 'SS        ','          ',  0,  0,  0),  &
chch_t(12, 'H2O2      ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(13, 'H2        ',  1, 'CT        ','          ',  0,  0,  0),  &
chch_t(14, 'CH4       ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(15, 'CO        ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(16, 'CO2       ',  1, 'CT        ','          ',  0,  0,  0),  &
chch_t(17, 'HCHO      ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(18, 'MeOO      ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(19, 'MeOOH     ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(20, 'MeOH      ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(21, 'MeO2NO2   ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(22, 'H2O       ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(23, 'O2        ',  1, 'CT        ','          ',  0,  0,  0),  &
chch_t(24, 'N2        ',  1, 'CT        ','          ',  0,  0,  0),  &
chch_t(25, 'Cl2       ',  2, 'TR        ','          ',  0,  0,  0),  &
chch_t(26, 'Cl        ',  1, 'FM        ','ClOx      ',  0,  0,  0),  &
chch_t(27, 'Cl2O2     ',  2, 'FM        ','ClOx      ',  0,  0,  0),  &
chch_t(28, 'ClO       ',  1, 'FM        ','ClOx      ',  0,  0,  0),  &
chch_t(29, 'OClO      ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(30, 'HOCl      ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(31, 'HCl       ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(32, 'ClONO2    ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(33, 'ClNO2     ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(34, 'HBr       ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(35, 'Br        ',  1, 'FM        ','BrOx      ',  0,  0,  0),  &
chch_t(36, 'HOBr      ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(37, 'BrO       ',  1, 'FM        ','BrOx      ',  0,  0,  0),  &
chch_t(38, 'BrONO2    ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(39, 'BrCl      ',  1, 'TR        ','          ',  0,  0,  0),  &
chch_t(40, 'Br2       ',  1, 'TR        ','          ',  0,  0,  0)  &
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
ratb_defs_std=(/              & 
ratb_t('O(3P)     ', 'O3        ', 'O2        ', 'O2        ', '          ',  &
'          ', 8.00E-12,   .00,   2060.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ', 'O2        ', 'O(3P)     ', 'O2        ', '          ',  &
'          ', 3.30E-11,   .00,    -55.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ', 'H2O       ', 'OH        ', 'OH        ', '          ',  &
'          ', 1.63E-10,   .00,    -60.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ', 'H2        ', 'OH        ', 'HO2       ', '          ',  &
'          ', 1.20E-10,   .00,       .0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ', 'N2        ', 'O(3P)     ', 'N2        ', '          ',  &
'          ', 2.15E-11,   .00,   -110.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ', 'CH4       ', 'OH        ', 'MeOO      ', '          ',  &
'          ', 1.75E-10,   .00,       .0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'O3        ', 'HO2       ', 'O2        ', '          ',  &
'          ', 1.70E-12,   .00,    940.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'HO2       ', 'H2O       ', 'O2        ', '          ',  &
'          ', 4.80E-11,   .00,   -250.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'H2O2      ', 'H2O       ', 'HO2       ', '          ',  &
'          ', 1.80E-12,   .00,       .0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'O3        ', 'OH        ', 'O2        ', 'O2        ',  &
'          ', 1.00E-14,   .00,    490.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'HO2       ', 'H2O2      ', 'O2        ', '          ',  &
'          ', 3.00E-13,   .00,   -460.0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('O(3P)     ', 'NO2       ', 'NO        ', 'O2        ', '          ',  &
'          ', 5.10E-12,   .00,   -210.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'NO3       ', 'HO2       ', 'NO2       ', '          ',  &
'          ', 2.20E-11,   .00,       .0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'HONO2     ', 'H2O       ', 'NO3       ', '          ',  &
'          ', 1.50E-13,   .00,       .0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('OH        ', 'HO2NO2    ', 'NO2       ', 'H2O       ', 'O2        ',  &
'          ', 1.30E-12,  0.00,   -380.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'NO        ', 'NO2       ', 'OH        ', '          ',  &
'          ', 3.30E-12,   .00,   -270.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'NO3       ', 'NO2       ', 'OH        ', 'O2        ',  &
'          ', 3.50E-12,   .00,       .0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO        ', 'O3        ', 'NO2       ', 'O2        ', '          ',  &
'          ', 3.00E-12,   .00,   1500.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO        ', 'NO3       ', 'NO2       ', 'NO2       ', '          ',  &
'          ', 1.50E-11,   .00,   -170.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO2       ', 'O3        ', 'NO3       ', 'O2        ', '          ',  &
'          ', 1.20E-13,   .00,   2450.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'CH4       ', 'MeOO      ', 'H2O       ', '          ',  &
'          ', 2.45E-12,   .00,   1775.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'HCHO      ', 'H2O       ', 'CO        ', 'HO2       ',  &
'          ', 5.50E-12,   .00,   -125.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'MeOH      ', 'H2O       ', 'HCHO      ', 'HO2       ',  &
'          ', 2.90E-12,   .00,    345.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'MeOOH     ', 'MeOO      ', 'H2O       ', '          ',  &
'          ', 2.70E-12,   .00,   -200.0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('OH        ', 'MeOOH     ', 'HCHO      ', 'OH        ', 'H2O       ',  &
'          ', 1.10E-12,   .00,   -200.0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('HO2       ', 'MeOO      ', 'MeOOH     ', 'O2        ', '          ',  &
'          ', 4.10E-13,   .00,   -750.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOO      ', 'MeOO      ', 'MeOH      ', 'HCHO      ', 'O2        ',  &
'          ', 9.50E-14,   .00,   -390.0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('MeOO      ', 'NO        ', 'NO2       ', 'HCHO      ', 'HO2       ',  &
'          ', 2.80E-12,   .00,   -300.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOOH     ', 'Cl        ', 'HCl       ', 'HCHO      ', 'OH        ',  &
'          ', 5.70E-11,   .00,      0.0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('O(3P)     ', 'ClO       ', 'Cl        ', 'O2        ', '          ',  &
'          ', 2.80E-11,   .00,    -85.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'Cl2       ', 'HOCl      ', 'Cl        ', '          ',  &
'          ', 2.60E-12,   .00,   1100.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'ClO       ', 'HO2       ', 'Cl        ', '          ',  &
'          ', 7.40E-12,   .00,   -270.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'ClO       ', 'HCl       ', 'O2        ', '          ',  &
'          ', 6.00E-13,   .00,   -230.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'HCl       ', 'H2O       ', 'Cl        ', '          ',  &
'          ', 1.80E-12,   .00,    250.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'HOCl      ', 'H2O       ', 'ClO       ', '          ',  &
'          ', 3.00E-12,   .00,    500.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'Cl        ', 'HCl       ', 'O2        ', '          ',  &
'          ', 1.40E-11,   .00,   -270.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'Cl        ', 'OH        ', 'ClO       ', '          ',  &
'          ', 3.60E-11,   .00,    375.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'ClO       ', 'HOCl      ', 'O2        ', '          ',  &
'          ', 2.60E-12,   .00,   -290.0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('HO2       ', 'ClO       ', 'HCl       ', 'O3        ', '          ',  &
'          ', 0.00E+00,   .00,       .0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('Cl        ', 'O3        ', 'ClO       ', 'O2        ', '          ',  &
'          ', 2.30E-11,   .00,    200.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Cl        ', 'H2        ', 'HCl       ', 'HO2       ', '          ',  &
'          ', 3.05E-11,   .00,   2270.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Cl        ', 'CH4       ', 'HCl       ', 'MeOO      ', '          ',  &
'          ', 7.10E-12,   .00,   1270.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Cl        ', 'HCHO      ', 'HCl       ', 'CO        ', 'HO2       ',  &
'          ', 8.10E-11,   .00,     30.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Cl        ', 'MeOH      ', 'HCl       ', 'HCHO      ', 'HO2       ',  &
'          ', 5.50E-11,   .00,       .0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Cl        ', 'OClO      ', 'ClO       ', 'ClO       ', '          ',  &
'          ', 3.40E-11,   .00,   -160.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Cl        ', 'HOCl      ', 'Cl2       ', 'OH        ', '          ',  &
'          ', 3.09E-12,   .00,    130.0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('Cl        ', 'HOCl      ', 'ClO       ', 'HCl       ', '          ',  &
'          ', 3.10E-13,   .00,    130.0, 0.000, 0.000, 0.000, 0.000) ,   & 
ratb_t('Cl        ', 'ClONO2    ', 'Cl2       ', 'NO3       ', '          ',  &
'          ', 6.50E-12,   .00,   -135.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('ClO       ', 'NO        ', 'NO2       ', 'Cl        ', '          ',  &
'          ', 6.40E-12,   .00,   -290.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('ClO       ', 'MeOO      ', 'Cl        ', 'HCHO      ', 'HO2       ',  &
'          ', 1.80E-12,   .00,    600.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(3P)     ', 'BrO       ', 'Br        ', 'O2        ', '          ',  &
'          ', 1.90E-11,   .00,   -230.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'HBr       ', 'H2O       ', 'Br        ', '          ',  &
'          ', 5.50E-12,   .00,   -200.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'Br        ', 'HBr       ', 'O2        ', '          ',  &
'          ', 4.80E-12,   .00,    310.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ', 'BrO       ', 'HOBr      ', 'O2        ', '          ',  &
'          ', 4.50E-12,   .00,   -460.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Br        ', 'O3        ', 'BrO       ', 'O2        ', '          ',  &
'          ', 1.60E-11,   .00,    780.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Br        ', 'HCHO      ', 'HBr       ', 'CO        ', 'HO2       ',  &
'          ', 1.70E-11,   .00,    800.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('BrO       ', 'NO        ', 'NO2       ', 'Br        ', '          ',  &
'          ', 8.80E-12,   .00,   -260.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('BrO       ', 'ClO       ', 'Br        ', 'OClO      ', '          ',  &
'          ', 9.50E-13,   .00,   -550.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('BrO       ', 'ClO       ', 'Br        ', 'Cl        ', 'O2        ',  &
'          ', 2.30E-12,   .00,   -260.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('BrO       ', 'ClO       ', 'BrCl      ', 'O2        ', '          ',  &
'          ', 4.10E-13,   .00,   -290.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('BrO       ', 'BrO       ', 'Br        ', 'Br        ', 'O2        ',  &
'          ', 2.40E-12,   .00,    -40.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('BrO       ', 'BrO       ', 'Br2       ', 'O2        ', '          ',  &
'          ', 2.80E-14,   .00,   -860.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(3P)     ', 'HOBr      ', 'OH        ', 'BrO       ', '          ',  &
'          ', 1.20E-10,   .00,    430.0, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ', 'Br2       ', 'HOBr      ', 'Br        ', '          ',  &
'          ', 2.10E-11,   .00,   -240.0, 0.000, 0.000, 0.000, 0.000)     &
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
rath_defs_std=(/ &
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



! ! describes photolytic reactions
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
ratj_defs_std = (/                  &
ratj_t('BrONO2    ','PHOTON    ','BrO       ','NO2       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'(?)       '),   &
ratj_t('BrONO2    ','PHOTON    ','Br        ','NO3       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'(?)       '),   &
ratj_t('BrCl      ','PHOTON    ','Br        ','Cl        ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('Cl2       ','PHOTON    ','Cl        ','Cl        ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('Cl2O2     ','PHOTON    ','Cl        ','Cl        ','O2         ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('ClNO2     ','PHOTON    ','Cl        ','NO2       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('ClONO2    ','PHOTON    ','Cl        ','NO3       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('ClONO2    ','PHOTON    ','ClO       ','NO2       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('H2O2      ','PHOTON    ','OH        ','OH        ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'          '),   &
ratj_t('HCHO      ','PHOTON    ','CO        ','HO2       ','HO2        ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'(?)       '),   &
ratj_t('HCHO      ','PHOTON    ','H2        ','CO        ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'(?)       '),   &
ratj_t('HO2NO2    ','PHOTON    ','HO2       ','NO2       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('HO2NO2    ','PHOTON    ','OH        ','NO3       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('HOBr      ','PHOTON    ','OH        ','Br        ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'(?)       '),   &
ratj_t('HOCl      ','PHOTON    ','OH        ','Cl        ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('HONO2     ','PHOTON    ','OH        ','NO2       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('MeO2NO2   ','PHOTON    ','MeOO      ','NO2       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'(?)       '),   &
ratj_t('MeOOH     ','PHOTON    ','HCHO      ','OH        ','HO2        ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('N2O5      ','PHOTON    ','NO3       ','NO2       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('NO2       ','PHOTON    ','NO        ','O(3P)     ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('NO3       ','PHOTON    ','NO        ','O2        ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('NO3       ','PHOTON    ','NO2       ','O(3P)     ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('O2        ','PHOTON    ','O(3P)     ','O(3P)     ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('O3        ','PHOTON    ','O2        ','O(3P)     ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('O3        ','PHOTON    ','O2        ','O(1D)     ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('OClO      ','PHOTON    ','O(3P)     ','ClO       ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       '),   &
ratj_t('Br2       ','PHOTON    ','Br        ','Br        ','           ', &
'          ',   0.0,   0.0,   0.0,   0.0, 100.0,'JPL       ')   &
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
ratt_defs_std=(/ &
ratt_t('O(3P)     ','O2        ','O3        ','m          ',  &
.00, 6.00E-34, -2.40,       .0, 0.00E+00,   .00,       .0, 0.000, 0.000) , &  
ratt_t('OH        ','NO2       ','HONO2     ','m          ',  &
.60, 1.80E-30, -3.00,       .0, 2.80E-11,   .00,       .0, 0.000, 0.000) , &  
ratt_t('HO2       ','NO2       ','HO2NO2    ','m          ',  &
.60, 1.90E-31, -3.40,       .0, 4.00E-12, -0.50,       .0, 0.000, 0.000) , &  
ratt_t('HO2NO2    ','m         ','HO2       ','NO2        ',  &
.60, 9.05E-05, -3.40,  10900.0, 1.90E+15, -0.50,  10900.0, 0.000, 0.000) , &  
ratt_t('NO2       ','NO3       ','N2O5      ','m          ',  &
.60, 2.40E-30, -3.00,       .0, 1.60E-12,  -.10,       .0, 0.000, 0.000) , &  
ratt_t('N2O5      ','m         ','NO2       ','NO3        ',  &
.60, 4.14E-04, -3.00,  10840.0, 2.76E+14,  -.10,  10840.0, 0.000, 0.000) , &  
ratt_t('MeOO      ','NO2       ','MeO2NO2   ','m          ',  &
.60, 1.00E-30, -4.80,       .0, 7.20E-12, -2.10,       .0, 0.000, 0.000) , &  
ratt_t('MeO2NO2   ','m         ','MeOO      ','NO2        ',  &
.60, 1.05E-02, -4.80,  11234.0, 7.58E+16, -2.10,  11234.0, 0.000, 0.000) , &  
ratt_t('ClO       ','NO2       ','ClONO2    ','m          ',  &
.60, 1.80E-31, -3.40,       .0, 1.50E-11, -1.90,       .0, 0.000, 0.000) , &  
ratt_t('ClO       ','ClO       ','Cl2O2     ','m          ',  &
.60, 1.90E-32, -3.60,       .0, 3.70E-12, -1.60,       .0, 0.000, 0.000) , &  
ratt_t('Cl2O2     ','m         ','ClO       ','ClO        ',  &
.60, 8.80E-06, -3.60,   8537.0, 1.71E+15, -1.60,   8537.0, 0.000, 0.000) , &  
ratt_t('BrO       ','NO2       ','BrONO2    ','m          ',  &
.60, 5.40E-31, -3.10,       .0, 6.50E-12, -2.90,       .0, 0.000, 0.000) , &  
ratt_t('OH        ','CO        ','CO2       ','HO2        ',  &
.60, 5.90E-33, -1.00,       .0, 1.10E-12,  1.30,       .0, 0.000, 0.000)  &  
  /)

   END SUBROUTINE CLAMS_CHEM_INIT_STANDARD

   SUBROUTINE CLAMS_CHEM_CLEAN_STANDARD

     implicit none

     DEALLOCATE(chch_defs_std)
     DEALLOCATE(ratb_defs_std)
     DEALLOCATE(ratj_defs_std)
     DEALLOCATE(ratt_defs_std)
     DEALLOCATE(rath_defs_std)

   END SUBROUTINE CLAMS_CHEM_CLEAN_STANDARD

END MODULE MESSY_CLAMSCHEM_DATA_STANDARD
