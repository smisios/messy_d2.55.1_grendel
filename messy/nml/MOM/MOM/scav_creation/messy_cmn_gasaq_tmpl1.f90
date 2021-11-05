!*****************************************************************************
!                Time-stamp: <2018-07-03 19:13:24 joec_pa>
!*****************************************************************************

! This file contains physicochemical data of species that exist in the
! GAS phase as well as in the AQueous phase. Currently, it contains only
! the soluble species from the MECCA gas.spc file. However, further
! species could be added here as long as the name does not conflict with
! any MECCA names. Several subroutines are provided that can supply
! these data to any MESSy submodel that needs them, e.g. MECCA or SCAV.

! Usage:
! - step 1: During the initialization phase, "CALL cmn_gasaq_initialize"
!           from the base model.
! - step 2: Later, but also during the initialization phase, each
!           submodel that needs the data should create its own data
!           array with appropriate indices (e.g. ind_*) and fill it
!           using the function get_gasaq() which is provided here.
! - step 3: During the time loop, the submodels should access the
!           submodel-specific data array. get_gasaq() should not be used
!           during the time loop because it needs to find the values by
!           comparing the names.

MODULE messy_cmn_gasaq

  USE messy_main_constants_mem, ONLY: DP, STRLEN_KPPSPECIES, &
    MH, MC, MN, MF, MNa, MO, MS, MCl, MBr, MI, MHg
  USE messy_main_tools,         ONLY: str
  USE messy_main_blather,       ONLY: warning

  IMPLICIT NONE

  PRIVATE
  REAL(DP), PARAMETER :: DUMMY   = -999.999_dp
  INTEGER,  PARAMETER :: MAXSIZE = 1000 ! ! mz_sg_20160806: increased to fit MOM species
  INTEGER :: n_gasaq

  TYPE GASAQ_TYPE
    CHARACTER(STRLEN_KPPSPECIES) :: name ! species name
    REAL(DP) :: Henry_T0   ! Henry constant at T0 = 298 K [M/atm]
    REAL(DP) :: Henry_Tdep ! its temperature dependence [K]
    REAL(DP) :: alpha_T0   ! accommodation coefficient alpha at T0 = 298 K [1]
    REAL(DP) :: alpha_Tdep ! its temperature dependence [K]
    REAL(DP) :: dryreac    ! dryreac
    REAL(DP) :: pss        ! pseudo-soil-solubility
    REAL(DP) :: M          ! molar mass [kg]
  END TYPE GASAQ_TYPE
  TYPE(GASAQ_TYPE), DIMENSION(MAXSIZE) :: gasaq

  ! public subroutines and functions:
  PUBLIC :: cmn_gasaq_initialize, get_gasaq

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE cmn_gasaq_initialize(status)

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: status

    ! LOCAL
    LOGICAL, SAVE :: ldone = .FALSE.

    status = 0 ! status = okay
    IF (ldone) RETURN

    ! set default values:
    gasaq(:) = GASAQ_TYPE('               ', DUMMY, DUMMY, DUMMY, DUMMY, DUMMY,DUMMY, DUMMY)

    CALL def_all_species(status) ! define all species
!mz_ap_20190311+ create Henry/alpha for all species
    !IF (status/=0) RETURN
    CALL add_all_henry(status)   ! add Henry's law coefficients
    CALL add_all_alpha(status)   ! add accommodation coefficients
    CALL add_all_dryreac(status) ! add dryreac and pseudo-soil-solubility
    CALL final_check(status)
    ldone = .TRUE.

  END SUBROUTINE cmn_gasaq_initialize

  ! --------------------------------------------------------------------------

  SUBROUTINE def_all_species(status)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'def_all_species'

    n_gasaq = 0
    ! O:
    CALL add_species('O2',       MO*2.)
    CALL add_species('O3',       MO*3.)
    ! H:
    CALL add_species('OH',       MO+MH)
    CALL add_species('HO2',      MH+MO*2.)
    CALL add_species('H2O2',     MH*2.+MO*2.)
    CALL add_species('H2O',      MH*2.+MO)
    ! N:
    CALL add_species('NH3',      MN+MH*3.)
    CALL add_species('NO',       MN+MO)
    CALL add_species('NO2',      MN+MO*2.)
    CALL add_species('NO3',      MN+MO*3.)
    CALL add_species('N2O5',     MN*2.+MO*5.)
    CALL add_species('HONO',     MH+MO+MN+MO)
    CALL add_species('HNO3',     MH+MN+MO*3.)
    CALL add_species('HNO4',     MH+MN+MO*4.)
    ! C:
    CALL add_species('CH3OH',    MC+MH*4.+MO)
    CALL add_species('CH3O2',    MC+MH*3.+MO*2.)
    CALL add_species('CH3OOH',   MC+MH*4.+MO*2.)
    CALL add_species('CO2',      MC+MO*2.)
    CALL add_species('HCHO',     MC+MH*2.+MO)
    CALL add_species('HCOOH',    MC+MH*2.+MO*2.)

    CALL add_species('CH3CO2H',  MC*2.+MH*4.+MO*2.)
    CALL add_species('PAN',      MC*2.+MH*3.+MO*5.+MN)
    CALL add_species('C2H5O2',   MC*2.+MH*5.+MO*2.)
    CALL add_species('CH3CHO',   MC*2.+MH*4.+MO)
    ! mz_ht_20130510+
    CALL add_species('OXL',      MC*2.+MH*2.+MO*4.) ! only SCAV
    CALL add_species('GLYOX',    MC*2.+MH*2.+MO*2.)
    CALL add_species('HOCH2CO2H',MC*2.+MH*4.+MO*3.)
    CALL add_species('HOCH2CHO', MC*2.+MH*4.+MO*2.)
    ! mz_ht_20130510-

    CALL add_species('CH3COCH3', MC*3.+MH*6.+MO)
    ! mz_ht_20130510+
    CALL add_species('CH3COCO2H', MC*3.+MH*4.+MO*3.) ! pyruvic acid only SCAV
    CALL add_species('MGLYOX',   MC*3.+MH*4.+MO*2.)
    ! mz_ht_20130510-
    ! Cl:
    CALL add_species('Cl2',      MCl*2.)
    CALL add_species('HCl',      MH+MCl)
    CALL add_species('HOCl',     MH+MO+MCl)
    CALL add_species('ClNO3',    MCl+MN+MO*3.)
    ! Br:
    CALL add_species('Br2',      MBr*2.)
    CALL add_species('HBr',      MH+MBr)
    CALL add_species('HOBr',     MH+MO+MBr)
    CALL add_species('BrNO3',    MBr+MN+MO*3.)
    CALL add_species('BrCl',     MBr+MCl)
    ! I:
    CALL add_species('I2',       MI*2.)
    CALL add_species('IO',       MI+MO)
    CALL add_species('OIO',      MI+MO*2.)
    CALL add_species('I2O2',     MI*2.+MO*2.)
    CALL add_species('HI',       MH+MI)
    CALL add_species('HOI',      MH+MO+MI)
    CALL add_species('HIO3',     MH+MI+MO*3.)
    CALL add_species('INO2',     MI+MN+MO*2.)
    CALL add_species('INO3',     MI+MN+MO*3.)
    CALL add_species('ICl',      MI+MCl)
    CALL add_species('IBr',      MI+MBr)
    ! S:
    CALL add_species('SO2',      MS+MO*2.)
    CALL add_species('H2SO4',    MH*2.+MS+MO*4.)
    CALL add_species('CH3SO3H',  MC+MH*4.+MS+MO*3.)
    CALL add_species('DMS',      MC*2.+MH*6.+MS)
    CALL add_species('DMSO',     MC*2.+MH*6.+MS+MO)
    ! Hg:
    CALL add_species('Hg',       MHg)
    CALL add_species('HgO',      MHg+MO)
    CALL add_species('HgCl2',    MHg+MCl*2.)
    CALL add_species('HgBr2',    MHg+MBr*2.)
    CALL add_species('HgCl',     MHg+MCl)
    CALL add_species('HgBr',     MHg+MBr)
    CALL add_species('ClHgBr',   MHg+MCl+MBr)
    CALL add_species('BrHgOBr',  MHg+MO+MBr*2.)
    CALL add_species('ClHgOBr',  MHg+MO+MCl+MBr)
    ! mz_ht_20130510+
    ! Passive
    CALL add_species('SO2t',     MS + MO*2.)      ! only SCAV
    CALL add_species('NH50W',    MH + MN + MO*3.) ! only SCAV
    ! mz_ht_20130510-
