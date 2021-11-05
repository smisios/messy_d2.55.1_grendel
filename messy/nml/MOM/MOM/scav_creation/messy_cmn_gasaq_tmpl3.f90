
  CONTAINS

    SUBROUTINE add_henry(name, Henry_T0, Henry_Tdep)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)    :: name
      REAL(DP),         INTENT(IN)    :: Henry_T0
      REAL(DP),         INTENT(IN)    :: Henry_Tdep
      INTEGER :: i
      LOGICAL :: l_found
      l_found = .FALSE.
      DO i = 1, n_gasaq
        IF (TRIM(gasaq(i)%name)==name) THEN
          l_found = .TRUE.
          IF ((ABS(gasaq(i)%Henry_T0-DUMMY)>TINY(0._dp)).OR. &
            (ABS(gasaq(i)%Henry_Tdep-DUMMY)>TINY(0._dp))) THEN
            CALL warning("Henry for "//name//" has been added already.", substr)
            status = -1
          ELSE
            gasaq(i)%Henry_T0   = Henry_T0
            gasaq(i)%Henry_Tdep = Henry_Tdep
          ENDIF
        ENDIF
      ENDDO
      IF (.NOT.l_found) THEN
        CALL warning(name//" does not exist, cannot add Henry.", substr)
        status = -1
      ENDIF
    END SUBROUTINE add_henry

  END SUBROUTINE add_all_henry

  ! --------------------------------------------------------------------------

  SUBROUTINE add_all_alpha(status)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'add_all_alpha'

    ! default values:
    REAL(DP), PARAMETER :: alpha_T0   = 0.1_dp
    REAL(DP), PARAMETER :: alpha_Tdep = 0._dp

    ! The following definitions are read by alpha2tex.awk which
    ! transforms the data into a LaTeX table. Therefore, the syntax must be:
    ! "CALL add_alpha('XYZ', alpha0, minDHR) ! {&REF}"
    ! O:
    CALL add_alpha('O2',       0.01_dp,         2000._dp) ! {&&}
    CALL add_alpha('O3',       0.002_dp,      alpha_Tdep) ! {&&826}
    ! H:
    CALL add_alpha('OH',       0.01_dp,       alpha_Tdep) ! {&&1047}
    CALL add_alpha('HO2',      0.5_dp,        alpha_Tdep) ! {&1864}
    CALL add_alpha('H2O2',     0.077_dp,        3127._dp) ! {&32}
    CALL add_alpha('H2O',      0.0_dp,        alpha_Tdep) ! {&&}
    ! N:
    CALL add_alpha('NH3',      0.06_dp,       alpha_Tdep) ! {&&826}
    CALL add_alpha('NO',       5.0E-5_dp,     alpha_Tdep) ! {&&448}
    CALL add_alpha('NO2',      0.0015_dp,     alpha_Tdep) ! {&&176}
    CALL add_alpha('NO3',      0.04_dp,       alpha_Tdep) ! {&&1048}
    CALL add_alpha('N2O5',     alpha_T0,      alpha_Tdep) ! {&&826}
    CALL add_alpha('HONO',     0.04_dp,       alpha_Tdep) ! {&&826}
    CALL add_alpha('HNO3',     0.5_dp,        alpha_Tdep) ! {&&930}
    CALL add_alpha('HNO4',     alpha_T0,      alpha_Tdep) ! {&&826}
    ! C1:
    CALL add_alpha('CH3OH',    alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('CH3O2',    0.01_dp,         2000._dp) ! {&&}
    CALL add_alpha('CH3OOH',   0.0046_dp,       3273._dp) ! {&844}
    CALL add_alpha('CO2',      0.01_dp,         2000._dp) ! {&&}
    CALL add_alpha('HCHO',     0.04_dp,       alpha_Tdep) ! {&&826}
    CALL add_alpha('HCOOH',    0.014_dp,        3978._dp) ! {&826}
    ! C2:
    CALL add_alpha('CH3CO2H',  2.0E-2_dp,       4079._dp) ! {&2574}
    CALL add_alpha('PAN',      alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('C2H5O2',   alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('CH3CHO',   3.0E-2_dp,     alpha_Tdep) ! {&&}
    ! mz_ht_20130510+
    CALL add_alpha('OXL',      alpha_T0,      alpha_Tdep) ! {&&} only SCAV
    CALL add_alpha('GLYOX',    alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('HOCH2CO2H',alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('HOCH2CHO', alpha_T0,      alpha_Tdep) ! {&&}
    ! mz_ht_20130510-
    ! C3:
    CALL add_alpha('CH3COCH3', 3.72E-3_dp,      6395._dp) ! {&2574}
    ! mz_ht_20130510+
    CALL add_alpha('CH3COCO2H',alpha_T0,      alpha_Tdep) ! {&&} only SCAV
    CALL add_alpha('MGLYOX',   alpha_T0,      alpha_Tdep) ! {&&}
    ! mz_ht_20130510-

    ! Cl:
    CALL add_alpha('Cl2',      0.038_dp,        6546._dp) ! {&380}
    CALL add_alpha('HCl',      0.074_dp,        3072._dp) ! {&&1161}
    CALL add_alpha('HOCl',     0.5_dp,        alpha_Tdep) ! {&&}
    CALL add_alpha('ClNO3',    0.108_dp,      alpha_Tdep) ! {&&1647}
    ! Br:
    CALL add_alpha('Br2',      0.038_dp,        6546._dp) ! {&380}
    CALL add_alpha('HBr',      0.032_dp,        3940._dp) ! {&&1161}
    CALL add_alpha('HOBr',     0.5_dp,        alpha_Tdep) ! {&&930}
    CALL add_alpha('BrNO3',    0.063_dp,      alpha_Tdep) ! {&&1647}
    CALL add_alpha('BrCl',     0.038_dp,        6546._dp) ! {&&}
    ! I:
    CALL add_alpha('I2',       0.01_dp,         2000._dp) ! {&&}
    CALL add_alpha('IO',       0.5_dp,          2000._dp) ! {&&}
    CALL add_alpha('OIO',      0.01_dp,       alpha_Tdep) ! {&&}
    CALL add_alpha('I2O2',     alpha_T0,        2000._dp) ! {&&}
    CALL add_alpha('HI',       0.036_dp,        4130._dp) ! {&&1161}
    CALL add_alpha('HOI',      0.5_dp,        alpha_Tdep) ! {&&}
    CALL add_alpha('HIO3',     0.01_dp,       alpha_Tdep) ! {&&}
    CALL add_alpha('INO2',     alpha_T0,        2000._dp) ! {&&}
    CALL add_alpha('INO3',     alpha_T0,        2000._dp) ! {&&}
    CALL add_alpha('ICl',      0.018_dp,        2000._dp) ! {&2159}
    CALL add_alpha('IBr',      0.018_dp,        2000._dp) ! {&&}
    ! S:
    CALL add_alpha('SO2',      0.11_dp,       alpha_Tdep) ! {&826}
    CALL add_alpha('H2SO4',    0.65_dp,       alpha_Tdep) ! {&&1205}
    CALL add_alpha('CH3SO3H',  0.076_dp,        1762._dp) ! {&389}
    CALL add_alpha('DMS',      alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('DMSO',     0.048_dp,        2578._dp) ! {&389}
    ! Hg:
    CALL add_alpha('Hg',       alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('HgO',      alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('HgCl2',    alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('HgBr2',    alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('HgCl',     alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('HgBr',     alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('ClHgBr',   alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('BrHgOBr',  alpha_T0,      alpha_Tdep) ! {&&}
    CALL add_alpha('ClHgOBr',  alpha_T0,      alpha_Tdep) ! {&&}
    ! mz_ht_20130510+
    ! Passive
    CALL add_alpha('SO2t',     0.11_dp,       alpha_Tdep) ! {&&} only SCAV
    CALL add_alpha('NH50W',    0.5_dp,        alpha_Tdep) ! {&&} only SCAV
    ! mz_ht_20130510-
    ! mz_sg_20160930+: MOM species

