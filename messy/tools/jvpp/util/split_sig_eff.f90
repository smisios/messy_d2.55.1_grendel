! lf95 -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap split_sig_eff.f90

PROGRAM split_sig_eff

IMPLICIT NONE
CHARACTER(LEN=80) :: oneline, oneline2, outputfilename
INTEGER :: iostatus, ji
LOGICAL :: l_firstline = .FALSE.
CHARACTER(LEN=1) :: ji_string

OPEN(UNIT=30, FILE='xcheck', STATUS='UNKNOWN')
WRITE (30,'(A)') '#'

DO ji = 1,8

  WRITE(ji_string,'(I1.1)') ji ! interval number as string
  OPEN(UNIT=10, FILE='../../result-ok/sig_eff.'//ji_string, STATUS='OLD')

  WRITE (30,'(A)',ADVANCE='NO') 'cat '

  DO
    READ(10, '(A)', IOSTAT=iostatus) oneline
    IF (iostatus < 0) EXIT ! exit do loop at end of file
    oneline2 = ADJUSTL(oneline) ! remove leading white space
    IF (VERIFY(oneline2(1:1),'0123456789-')/=0) THEN
      IF (.NOT.l_firstline) CLOSE(20)
      SELECT CASE(TRIM(oneline2))
      CASE("br2")      ; oneline = "Br2"
      CASE("brcl_not") ; oneline = "BrCl"
      CASE("brno2")    ; oneline = "BrNO2"
      CASE("brno3")    ; oneline = "BrNO3"
      CASE("bro_not")  ; oneline = "BrO"
      CASE("c3h7i")    ; oneline = "C3H7I"
      CASE("ccl4")     ; oneline = "CCl4"
      CASE("cf2cl2")   ; oneline = "CF2Cl2"
      CASE("cf2clbr")  ; oneline = "CF2ClBr"
      CASE("cfcl3")    ; oneline = "CFCl3"
      CASE("ch2br2")   ; oneline = "CH2Br2"
      CASE("ch2clbr")  ; oneline = "CH2ClBr"
      CASE("ch2cli")   ; oneline = "CH2ClI"
      CASE("ch2i2")    ; oneline = "CH2I2"
      CASE("ch3br")    ; oneline = "CH3Br"
      CASE("ch3ccl3")  ; oneline = "CH3CCl3"
      CASE("ch3cho")   ; oneline = "CH3CHO"
      CASE("ch3cl")    ; oneline = "CH3Cl"
      CASE("ch3coch3") ; oneline = "CH3COCH3"
      CASE("ch3cocho") ; oneline = "MGLYOX"
      CASE("ch3i")     ; oneline = "CH3I"
      CASE("ch3ooh")   ; oneline = "CH3OOH"
      CASE("chbr3")    ; oneline = "CHBr3"
      CASE("chcl2br")  ; oneline = "CHCl2Br"
      CASE("chclbr2")  ; oneline = "CHClBr2"
      CASE("choh")     ; oneline = "CHOH"
      CASE("cl2")      ; oneline = "Cl2"
      CASE("cl2o2")    ; oneline = "Cl2O2"
      CASE("clno2")    ; oneline = "ClNO2"
      CASE("clno3")    ; oneline = "ClNO3"
      CASE("co2")      ; oneline = "CO2"
      CASE("coh2")     ; oneline = "COH2"
      CASE("h2o2")     ; oneline = "H2O2"
      CASE("h2o")      ; oneline = "H2O"
      CASE("hcl")      ; oneline = "HCl"
      CASE("h_no2")    ; oneline = "h_NO2"
      CASE("hno3")     ; oneline = "HNO3"
      CASE("hno4")     ; oneline = "HNO4"
      CASE("h_o2")     ; oneline = "h_O2"
      CASE("h_o3")     ; oneline = "h_O3"
      CASE("hobr")     ; oneline = "HOBr"
      CASE("hocl")     ; oneline = "HOCl"
      CASE("hoi")      ; oneline = "HOI"
      CASE("hono")     ; oneline = "HONO"
      CASE("i2")       ; oneline = "I2"
      CASE("ibr")      ; oneline = "IBr"
      CASE("icl")      ; oneline = "ICl"
      CASE("ino2")     ; oneline = "INO2"
      CASE("ino3")     ; oneline = "INO3"
      CASE("io")       ; oneline = "IO"
      CASE("n2o5")     ; oneline = "N2O5"
      CASE("n2o")      ; oneline = "N2O"
      CASE("no2")      ; oneline = "NO2"
      CASE("no2o")     ; oneline = "NO2O"
      CASE("noo2")     ; oneline = "NOO2"
      CASE("o1d")      ; oneline = "O1D"
      CASE("o2")       ; oneline = "O2"
      CASE("o3p")      ; oneline = "O3P"
      CASE("oclo_not") ; oneline = "OClO"
      CASE("paa")      ; oneline = "CH3CO3H"
      CASE("pan")      ; oneline = "PAN"
      !CASE DEFAULT ; oneline = "notused_"//TRIM(oneline2)
      CASE DEFAULT ; oneline = "dummy"
      END SELECT
      outputfilename = TRIM(oneline(1:INDEX(oneline,' ')-1))//'.ef'//ji_string
      print *, TRIM(outputfilename)
      WRITE (30,'(A)',ADVANCE='NO') TRIM(outputfilename)//' '
      OPEN(UNIT=20, FILE='../spectra_eff_static_old/'//TRIM(outputfilename), STATUS='UNKNOWN')
      WRITE (20,'(A)') '# '//TRIM(oneline)
    ELSE
      WRITE (20,'(A)') TRIM(oneline)
    ENDIF
  ENDDO
  CLOSE(20)
  CLOSE(10)

  WRITE (30,'(A)') '> sig_eff.'//ji_string//'_repr'
  WRITE (30,'(A)') 'diff sig_eff.'//ji_string//'_repr sig_eff.'//ji_string

ENDDO

CLOSE(30)

END PROGRAM split_sig_eff
