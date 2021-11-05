
  CONTAINS

    SUBROUTINE add_alpha(name, alpha_T0, alpha_Tdep)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)    :: name
      REAL(DP),             INTENT(IN)    :: alpha_T0
      REAL(DP),             INTENT(IN)    :: alpha_Tdep
      INTEGER :: i
      LOGICAL :: l_found
      l_found = .FALSE.
      DO i = 1, n_gasaq
        IF (TRIM(gasaq(i)%name)==name) THEN
          l_found = .TRUE.
          IF ((ABS(gasaq(i)%alpha_T0-DUMMY)>TINY(0._dp)).OR. &
            (ABS(gasaq(i)%alpha_Tdep-DUMMY)>TINY(0._dp))) THEN
            CALL warning("alpha for "//name//" has been added already.", substr)
            status = -1
          ELSE
            gasaq(i)%alpha_T0   = alpha_T0
            gasaq(i)%alpha_Tdep = alpha_Tdep
          ENDIF
        ENDIF
      ENDDO
      IF (.NOT.l_found) THEN
        CALL warning(name//" does not exist, cannot add alpha.", substr)
        status = -1
      ENDIF
    END SUBROUTINE add_alpha

  END SUBROUTINE add_all_alpha

  ! --------------------------------------------------------------------------

  SUBROUTINE final_check(status)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'final_check'
    INTEGER :: i

    DO i = 1, n_gasaq
      IF ((ABS(gasaq(i)%Henry_T0-DUMMY)<TINY(0._dp)).OR. &
        (ABS(gasaq(i)%Henry_Tdep-DUMMY)<TINY(0._dp))) THEN
        CALL warning("Henry for "//TRIM(gasaq(i)%name)// &
          " was never added.", substr)
        status = -1
      ENDIF
      IF ((ABS(gasaq(i)%alpha_T0-DUMMY)<TINY(0._dp)).OR. &
        (ABS(gasaq(i)%alpha_Tdep-DUMMY)<TINY(0._dp))) THEN
        CALL warning("alpha for "//TRIM(gasaq(i)%name)// &
          " was never added.", substr)
        status = -1
      ENDIF
    ENDDO

  END SUBROUTINE final_check

  ! --------------------------------------------------------------------------

  INTEGER FUNCTION get_gasaq(name, &
    Henry_T0, Henry_Tdep, alpha_T0, alpha_Tdep, M, pss, dryreac)

    IMPLICIT NONE

    CHARACTER(LEN=*),   INTENT(IN)  :: name
    REAL(DP), OPTIONAL, INTENT(OUT) :: Henry_T0
    REAL(DP), OPTIONAL, INTENT(OUT) :: Henry_Tdep
    REAL(DP), OPTIONAL, INTENT(OUT) :: alpha_T0
    REAL(DP), OPTIONAL, INTENT(OUT) :: alpha_Tdep
    REAL(DP), OPTIONAL, INTENT(OUT) :: M
    REAL(DP), OPTIONAL, INTENT(OUT) :: pss
    REAL(DP), OPTIONAL, INTENT(OUT) :: dryreac

    INTEGER :: i

    get_gasaq = 1 ! set status to error until species "name" is found
    DO i = 1, MAXSIZE
      IF (TRIM(gasaq(i)%name)==name) THEN
        IF(PRESENT(Henry_T0))   Henry_T0   = gasaq(i)%Henry_T0
        IF(PRESENT(Henry_Tdep)) Henry_Tdep = gasaq(i)%Henry_Tdep
        IF(PRESENT(alpha_T0))   alpha_T0   = gasaq(i)%alpha_T0
        IF(PRESENT(alpha_Tdep)) alpha_Tdep = gasaq(i)%alpha_Tdep
        IF(PRESENT(M))          M          = gasaq(i)%M
        IF(PRESENT(pss))        pss        = gasaq(i)%pss
        IF(PRESENT(dryreac))    dryreac    = gasaq(i)%dryreac
        get_gasaq = 0 ! status = okay
      ENDIF
    ENDDO

  END FUNCTION get_gasaq

!-------------------------------------------------------------------------

 SUBROUTINE add_all_dryreac(status)

   ! setting up pseudo soil solubility and dryreac_sf, 
   ! both required by drydep

    USE messy_main_constants_mem, ONLY: HUGE_DP
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'add_all_dryreac'

    ! O:
    CALL add_dryreac('O2',       0.0_dp,             0.0_dp) 
    CALL add_dryreac('O3',       1.0_dp,             0.01_dp)
    ! H:
    CALL add_dryreac('OH',       0.0_dp,             25._dp)
    CALL add_dryreac('HO2',      0.0_dp,             0.0_dp)
    CALL add_dryreac('H2O2',     1.0_dp,             7.45e4_dp)
    CALL add_dryreac('H2O',      0.0_dp,             0.0_dp)
    ! N:
    CALL add_dryreac('NH3',      1.0_dp,             2.e4_dp) 
    CALL add_dryreac('NO',       1.0_dp,             2.e-3_dp)
    CALL add_dryreac('NO2',      1.0_dp,             1.0e-2_dp)
    CALL add_dryreac('NO3',      1.0_dp,             1.8_dp)
    CALL add_dryreac('N2O5',     1.0_dp,             1.e4_dp)
    CALL add_dryreac('HONO',     0.1_dp,             4.9E1_dp)
    CALL add_dryreac('HNO3',     1.0_dp,             1.e4_dp)
    CALL add_dryreac('HNO4',     1.0_dp,             1.e4_dp)
    ! C:
    CALL add_dryreac('CH3OH',    0.1_dp,             2.2e2_dp)
    CALL add_dryreac('CH3O2',    0.0_dp,             0.0_dp)
    CALL add_dryreac('CH3OOH',   0.1_dp,             3.e2_dp)
    CALL add_dryreac('CO2',      0.0_dp,             0.0_dp)
    CALL add_dryreac('HCHO',     0.1_dp,             3.2e3_dp)
    CALL add_dryreac('HCOOH',    0.0_dp,             4.e6_dp)

    CALL add_dryreac('CH3CO2H',  0.1_dp,             4.1e3_dp)
    CALL add_dryreac('PAN',      0.1_dp,             2.8_dp)
    CALL add_dryreac('C2H5O2',   0.0_dp,             0.0_dp)
    CALL add_dryreac('CH3CHO',   0.1_dp,             1.3e1_dp)
    ! mz_ht_20130510+
    CALL add_dryreac('OXL',      0.1_dp,             3.26E6_dp) ! only SCAV
    CALL add_dryreac('GLYOX',    0.1_dp,             4.2E5_dp)
    CALL add_dryreac('HOCH2CO2H',0.1_dp,             1.1E4_dp)
    CALL add_dryreac('HOCH2CHO', 0.1_dp,             4.10E4_dp)
    ! mz_ht_20130510-

    CALL add_dryreac('CH3COCH3', 0.1_dp,             3.e1_dp)
    ! mz_ht_20130510+
    CALL add_dryreac('CH3COCO2H',0.1_dp,             3.10E5_dp) ! only SCAV
    CALL add_dryreac('MGLYOX',   0.1_dp,             3.7E3_dp)
    ! mz_ht_20130510-
    ! Cl:
    CALL add_dryreac('Cl2',      0.1_dp,             7.e-2_dp)
    CALL add_dryreac('HCl',      1.0_dp,             1e14_dp)
    CALL add_dryreac('HOCl',     1.0_dp,             6.7e2_dp)
    CALL add_dryreac('ClNO3',    1.0_dp,             1.e30_dp)
    ! Br:
    CALL add_dryreac('Br2',      0.1_dp,             0.7_dp)
    CALL add_dryreac('HBr',      1.0_dp,             1.3e17_dp)
    CALL add_dryreac('HOBr',     1.0_dp,             9.1e1_dp)
    CALL add_dryreac('BrNO3',    1.0_dp,             1.0e30_dp)
    CALL add_dryreac('BrCl',     0.1_dp,             1.0_dp)
    ! I:
    CALL add_dryreac('I2',       0.1_dp,             3.0_dp)
    CALL add_dryreac('IO',       0.0_dp,             0.0_dp)
    CALL add_dryreac('OIO',      0.0_dp,             0.0_dp)
    CALL add_dryreac('I2O2',     0.0_dp,             0.0_dp)
    CALL add_dryreac('HI',       1.0_dp,             1.e30_dp)
    CALL add_dryreac('HOI',      1.0_dp,             1.e30_dp)
    CALL add_dryreac('HIO3',     1.0_dp,             1.e30_dp)
    CALL add_dryreac('INO2',     0.1_dp,             4.5_dp)
    CALL add_dryreac('INO3',     1.0_dp,             1.e30_dp)
    CALL add_dryreac('ICl',      0.1_dp,             1.1e2_dp)
    CALL add_dryreac('IBr',      0.0_dp,             0.0_dp)
    ! S:
    CALL add_dryreac('SO2',      0.0_dp,             1.2_dp)
    CALL add_dryreac('H2SO4',    0.0_dp,             0.0_dp)

    CALL add_dryreac('CH3SO3H',  1.0_dp,             1.e30_dp)
    CALL add_dryreac('DMS',      0.0_dp,             0.0_dp)
    CALL add_dryreac('DMSO',     0.1_dp,             5.e4_dp)
    ! Hg:
    CALL add_dryreac('Hg',       0.1_dp,             0.13_dp)
    CALL add_dryreac('HgO',      1.0_dp,             2.4E7_dp)
    CALL add_dryreac('HgCl2',    1.0_dp,             2.4E7_dp)
    CALL add_dryreac('HgBr2',    1.0_dp,             2.4E7_dp)
    CALL add_dryreac('HgCl',     1.0_dp,             2.4E7_dp)
    CALL add_dryreac('HgBr',     1.0_dp,             2.4E7_dp)
    CALL add_dryreac('ClHgBr',   1.0_dp,             2.4E7_dp)
    CALL add_dryreac('BrHgOBr',  1.0_dp,             2.4E7_dp)
    CALL add_dryreac('ClHgOBr',  1.0_dp,             2.4E7_dp)
    ! mz_ht_20130510+
    ! Passive
    CALL add_dryreac('SO2t',     0.0_dp,             1.2_dp)   ! only SCAV
    CALL add_dryreac('NH50W',    1.0_dp,             1.0E4_dp) ! only SCAV
    ! mz_ht_20130510-

  CONTAINS

    SUBROUTINE add_dryreac(name, dryreac, pss)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)    :: name
      REAL(DP),         INTENT(IN)    :: dryreac
      REAL(DP),         INTENT(IN)    :: pss
      INTEGER :: i
      LOGICAL :: l_found
      l_found = .FALSE.
      DO i = 1, n_gasaq
        IF (TRIM(gasaq(i)%name)==name) THEN
          l_found = .TRUE.
          IF ((ABS(gasaq(i)%dryreac-DUMMY)>TINY(0._dp)).OR. &
            (ABS(gasaq(i)%pss-DUMMY)>TINY(0._dp))) THEN
            CALL warning("Drydep values for "//name//" have been added already.", substr)
            status = -1
          ELSE
            gasaq(i)%dryreac = dryreac
            gasaq(i)%pss     = pss
          ENDIF
        ENDIF
      ENDDO
      IF (.NOT.l_found) THEN
        CALL warning(name//" does not exist, cannot add drydep values.", substr)
        status = -1
      ENDIF
    END SUBROUTINE add_dryreac

  END SUBROUTINE add_all_dryreac

!*****************************************************************************
END MODULE messy_cmn_gasaq
!*****************************************************************************

