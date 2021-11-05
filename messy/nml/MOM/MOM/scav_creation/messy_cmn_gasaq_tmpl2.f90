    IF (n_gasaq>MAXSIZE) THEN
      CALL warning("MAXSIZE too small, set to >= "//str(n_gasaq), substr)
      status = -1
    ENDIF

  CONTAINS

    SUBROUTINE add_species(name, M)

      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL(DP),         INTENT(IN) :: M     ! molar mass [g/mol]
      INTEGER :: i

      n_gasaq = n_gasaq + 1
      IF (n_gasaq>MAXSIZE) RETURN

      ! loop over previously defined species (loop is skipped completely
      ! if this is the first call of add_species):
      DO i = 1, n_gasaq
        IF (TRIM(gasaq(i)%name)==name) THEN
          CALL warning("Species "//name//" has been defined already.", substr)
          status = -1
        ENDIF
      ENDDO

      gasaq(n_gasaq)%name = name
      gasaq(n_gasaq)%M    = M / 1E3_dp ! converted to [kg/mol]
    END SUBROUTINE add_species

  END SUBROUTINE def_all_species

  ! --------------------------------------------------------------------------

  SUBROUTINE add_all_henry(status)

    USE messy_main_constants_mem, ONLY: BIG_DP
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'add_all_henry'

    ! The following definitions are read by henry2tex.awk which
    ! transforms the data into a LaTeX table. Therefore, the syntax must be:
    ! "CALL add_henry('XYZ', KH, minDHR) ! {&REF}"
    ! O:
    CALL add_henry('O2',       1.3E-3_dp,          1500._dp) ! {&190}
    CALL add_henry('O3',       1.2E-2_dp,          2560._dp) ! {&87}
    ! H:
    CALL add_henry('OH',       3.0E1_dp,           4300._dp) ! {&515}
    CALL add_henry('HO2',      3.9E3_dp,           5900._dp) ! {&515}
    CALL add_henry('H2O2',     1.E5_dp,            6338._dp) ! {&311}
    CALL add_henry('H2O',      BIG_DP,                0._dp) ! {&&}
    ! N:
    CALL add_henry('NH3',      58._dp,             4085._dp) ! {&87}
    CALL add_henry('NO',       1.9E-3_dp,          1480._dp) ! {&449}
    CALL add_henry('NO2',      7.0E-3_dp,          2500._dp) ! {&&59}
    CALL add_henry('NO3',      2._dp,              2000._dp) ! {&219}
    CALL add_henry('N2O5',     BIG_DP,                0._dp) ! {&&}
    CALL add_henry('HONO',     4.9E1_dp,           4780._dp) ! {&449}
    CALL add_henry('HNO3',     2.45E6_dp/1.5E1_dp, 8694._dp) ! {&&530}
    CALL add_henry('HNO4',     1.2E4_dp,           6900._dp) ! {&797}
    ! C1:
    CALL add_henry('CH3OH',    2.20E2_dp,          5200._dp) ! {&483}
    CALL add_henry('CH3O2',    6._dp,              5600._dp) ! {&&46}
    CALL add_henry('CH3OOH',   3.0E2_dp,           5322._dp) ! {&311}
    CALL add_henry('CO2',      3.1E-2_dp,          2423._dp) ! {&87}
    CALL add_henry('HCHO',     7.0E3_dp,           6425._dp) ! {&87}
    CALL add_henry('HCOOH',    3.7E3_dp,           5700._dp) ! {&87}
    ! C2:
    CALL add_henry('CH3CO2H',  4.1E3_dp,           6200._dp) ! {&1945}
    CALL add_henry('PAN',      2.8_dp,             5730._dp) ! {&1945}
    CALL add_henry('C2H5O2',   6._dp,              5600._dp) ! {&&}
    CALL add_henry('CH3CHO',   1.29E1_dp,          5890._dp) ! {&1945}
    ! mz_ht_20130510+
    CALL add_henry('OXL',      3.26E6_dp,          7285._dp) ! {&3062} only SCAV
    CALL add_henry('GLYOX',    4.19E5_dp,          7481._dp) ! {&2626}
    CALL add_henry('HOCH2CO2H',2.83E4_dp,          4029._dp) ! {&2401}
    CALL add_henry('HOCH2CHO', 4.10E4_dp,          4600._dp) ! {&484}
    ! mz_ht_20130510-
    ! C3:
    CALL add_henry('CH3COCH3', 28.1_dp,            5050._dp) ! {&1945}
    ! mz_ht_20130510+
    CALL add_henry('CH3COCO2H',3.11E5_dp,          5090._dp) ! {&2626} only SCAV
    CALL add_henry('MGLYOX',   3.70E3_dp,          7500._dp) ! {&484}
    ! mz_ht_20130510-
    ! Cl:
    CALL add_henry('Cl2',      9.2E-2_dp,          2081._dp) ! {&1038}
    CALL add_henry('HCl',      2./1.7_dp,          9001._dp) ! {&530}
    CALL add_henry('HOCl',     6.6E2_dp,           5862._dp) ! {&315}
    CALL add_henry('ClNO3',    BIG_DP,                0._dp) ! {&&}
    ! Br:
    CALL add_henry('Br2',      7.7E-1_dp,          3837._dp) ! {&1038}
    CALL add_henry('HBr',      1.3_dp,            10239._dp) ! {&&530}
    CALL add_henry('HOBr',     1.3E3_dp,           5862._dp) ! {&&288}
    CALL add_henry('BrNO3',    BIG_DP,                0._dp) ! {&&}
    CALL add_henry('BrCl',     9.4E-1_dp,          5600._dp) ! {&1038}
    ! I:
    CALL add_henry('I2',       3._dp,              4431._dp) ! {&582}
    CALL add_henry('IO',       4.5E2_dp,           5862._dp) ! {&&}
    CALL add_henry('OIO',      BIG_DP,                0._dp) ! {&&}
    CALL add_henry('I2O2',     BIG_DP,                0._dp) ! {&&}
    CALL add_henry('HI',       BIG_DP,                0._dp) ! {&&}
    CALL add_henry('HOI',      4.5E2_dp,           5862._dp) ! {&&162}
    CALL add_henry('HIO3',     BIG_DP,                0._dp) ! {&&}
    CALL add_henry('INO2',     BIG_DP,                0._dp) ! {&&}
    CALL add_henry('INO3',     BIG_DP,                0._dp) ! {&&}
    CALL add_henry('ICl',      1.1E2_dp,           5600._dp) ! {&&}
    CALL add_henry('IBr',      2.4E1_dp,           5600._dp) ! {&&}
    ! S:
    CALL add_henry('SO2',      1.2_dp,             3120._dp) ! {&87}
    CALL add_henry('H2SO4',    1.E11_dp,              0._dp) ! {&&}
    CALL add_henry('CH3SO3H',  BIG_DP,                0._dp) ! {&&}
    CALL add_henry('DMS',      5.4E-1_dp,          3500._dp) ! {&1525}
    CALL add_henry('DMSO',     5.E4_dp,            6425._dp) ! {&&389}
    ! Hg:
    CALL add_henry('Hg',       0.13_dp,               0._dp) ! {&2171}
    CALL add_henry('HgO',      3.2E6_dp,              0._dp) ! {&2285}
    CALL add_henry('HgCl2',    2.4E7_dp,              0._dp) ! {&2285}
    CALL add_henry('HgBr2',    2.4E7_dp,              0._dp) ! {&&}
    CALL add_henry('HgCl',     2.4E7_dp,              0._dp) ! {&2285}
    CALL add_henry('HgBr',     2.4E7_dp,              0._dp) ! {&&}
    CALL add_henry('ClHgBr',   2.4E7_dp,              0._dp) ! {&&}
    CALL add_henry('BrHgOBr',  2.4E7_dp,              0._dp) ! {&&}
    CALL add_henry('ClHgOBr',  2.4E7_dp,              0._dp) ! {&&}
    ! mz_ht_20130510+
    ! Passive
    CALL add_henry('SO2t',     4.2e3_dp,           3120._dp) ! {&&} only SCAV
    CALL add_henry('NH50W',    3.5E12_dp,          8694._dp) ! {&&} only SCAV
    ! mz_ht_20130510-

