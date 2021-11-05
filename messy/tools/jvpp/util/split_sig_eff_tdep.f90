PROGRAM split_sig_eff_tdep

IMPLICIT NONE
CHARACTER(LEN=80) :: oneline, oneline2, outputfilename, inputfilename
INTEGER :: iostatus, ji, jt, maxv3, maxv2
CHARACTER(LEN=1) :: ji_string
CHARACTER(LEN=2) :: temp_string, temp_string2

DO jt = 1,15 ! temp
  WRITE(temp_string, '(I2)')   jt ! temperature number as string
  WRITE(temp_string2,'(I2.2)') jt ! temperature number as string (2 digits)
  DO ji = 1,8 ! interval
    WRITE(ji_string,'(I1.1)') ji ! interval number as string

    inputfilename = '../../result-ok-tdep/'//ji_string// &
      '/sig_temp.'//TRIM(ADJUSTL(temp_string))
    PRINT *, "inputfilename = ", TRIM(inputfilename)
    OPEN(UNIT=10, FILE=TRIM(inputfilename), STATUS='OLD')

    IF (ji==1) THEN ! Schumann-Runge bands
      outputfilename = '../spectra_eff_static_tdep_old/v3_du.ef'//ji_string// &
        '_T'//TRIM(ADJUSTL(temp_string2))
      PRINT *, "outputfilename = ", TRIM(outputfilename)
      OPEN(UNIT=20, FILE=TRIM(outputfilename), STATUS='UNKNOWN')
      READ(10, '(A)', IOSTAT=iostatus) oneline
      WRITE (20,'(A)') TRIM(oneline) ! maxv3 = 4 values in one line
      CLOSE(20)
      outputfilename = '../spectra_eff_static_tdep_old/v2.ef'//ji_string// &
        '_T'//TRIM(ADJUSTL(temp_string2))
      PRINT *, "outputfilename = ", TRIM(outputfilename)
      OPEN(UNIT=20, FILE=TRIM(outputfilename), STATUS='UNKNOWN') ! maxv2 = 231
    ELSE ! interval 2-8
      outputfilename = '../spectra_eff_static_tdep_old/v3_du.ef'//ji_string// &
        '_T'//TRIM(ADJUSTL(temp_string2))
      PRINT *, "outputfilename = ", TRIM(outputfilename)
      OPEN(UNIT=20, FILE=TRIM(outputfilename), STATUS='UNKNOWN') ! maxv3 = 600
    ENDIF

    DO
      READ(10, '(A)', IOSTAT=iostatus) oneline
      IF (iostatus < 0) EXIT ! exit do loop at end of file
      oneline2 = ADJUSTL(oneline) ! remove leading white space
      IF (VERIFY(oneline2(1:1),'0123456789-')/=0) THEN
        CLOSE(20)
        SELECT CASE(TRIM(oneline2))
        CASE("CFC_12")   ; oneline = "CF2Cl2"
        CASE("CFC-11")   ; oneline = "CFCl3"
        CASE("CH2Br2")   ; oneline = "CH2Br2"
        CASE("CH3CCL3")  ; oneline = "CH3CCl3"
        CASE("CH3CL")    ; oneline = "CH3Cl"
        CASE("CH3COCH3") ; oneline = "CH3COCH3"
        CASE("CHBr3")    ; oneline = "CHBr3"
        CASE("CHOH")     ; oneline = "CHOH"
        CASE("CLONO2")   ; oneline = "ClNO3"
        CASE("COH2")     ; oneline = "COH2"
        CASE("H2O2")     ; oneline = "H2O2"
        CASE("HNO3")     ; oneline = "HNO3"
        CASE("HNO4")     ; oneline = "HNO4"
        CASE("N2O")      ; oneline = "N2O"
        CASE("N2O5")     ; oneline = "N2O5"
        CASE("NO2")      ; oneline = "NO2"
        CASE("NO2O")     ; oneline = "NO2O"
        CASE("NOO2")     ; oneline = "NOO2"
        CASE("O1D")      ; oneline = "O1D"
        CASE("O2")       ; oneline = "O2"
        CASE("O3P")      ; oneline = "O3P"
        CASE("PAN")      ; oneline = "PAN"
        !CASE DEFAULT ; oneline = "notused_"//TRIM(oneline2)
        CASE DEFAULT ; oneline = "dummy"
        END SELECT
        outputfilename = '../spectra_eff_static_tdep_old/'// &
          TRIM(oneline(1:INDEX(oneline,' ')-1))// &
          '.ef'//ji_string//'_T'//TRIM(ADJUSTL(temp_string2))
        PRINT *, "outputfilename = ", TRIM(outputfilename)
        OPEN(UNIT=20, FILE= TRIM(outputfilename), STATUS='UNKNOWN')
      ENDIF
      WRITE (20,'(A)') TRIM(oneline)
    ENDDO
    CLOSE(20)
    CLOSE(10)

  ENDDO
ENDDO

END PROGRAM split_sig_eff_tdep
