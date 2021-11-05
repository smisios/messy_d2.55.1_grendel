! ****************************************************************************
!                Time-stamp: <2020-09-11 23:20:54 sander>
!     Author: Rolf Sander, Max-Planck Institute, Mainz, Germany, 2009-2014
!*****************************************************************************

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

! ****************************************************************************

MODULE jvpp_step3

  USE messy_main_constants_mem, ONLY: dp
  USE jvpp_mem, ONLY: NBIN, HLINE1, IO_TMP, IO_JTC, IO_JTD, IO_ERR, IO_REC, &
                      IO_NML, IO_F90, IO_CAT, IO_CAL, IO_LOG, IO_TEX, &
                      DIR_EFF, DIR_F90, DIR_JNL, l_species, l_tdep, inputdir

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: INCFILE = '../messy_jval_jvpp.inc'
  INTEGER, PARAMETER :: DIM55   = 55
  INTEGER, PARAMETER :: DIM58   = 58
  INTEGER, PARAMETER :: NTEMP   = 15 ! number of temperatures
  INTEGER :: ji, ji0, maxv2, maxv3
  CHARACTER(LEN=1) :: ji_string, ji0_string
  CHARACTER(LEN=80) :: header
  CHARACTER(LEN=80) :: v2_filename, v3_du_filename, xxx_eff_filename
  CHARACTER(LEN=80) :: recalc_filename, error_filename
  CHARACTER(LEN=5) :: suffix
  REAL(DP) :: T_ref      ! reference temperature [K]

  ! errors for deg_tconst and deg_tdep:
  INTEGER, PARAMETER :: FIT_PROBLEM = -111
  INTEGER, PARAMETER :: ALL_ZERO    = -222
  INTEGER, PARAMETER :: NO_TDEP     = -333

  ! parameters for the 8 intervals (T-const):
  REAL(DP), DIMENSION(DIM58) :: a0_1_1_xxx, a0_2_1_xxx, a0_3_1_xxx, a0_4_1_xxx
  REAL(DP), DIMENSION(DIM58) :: a0_1_2_xxx, a0_2_2_xxx, a0_3_2_xxx, a0_4_2_xxx
  REAL(DP), DIMENSION(DIM55) :: a1_1_xxx, a1_2_xxx, a1_3_xxx
  REAL(DP), DIMENSION(DIM55) :: b1_1_xxx, b1_2_xxx, b1_3_xxx
  REAL(DP), DIMENSION(2:3,DIM55) :: a23_xxx, b23_xxx
  REAL(DP), DIMENSION(4:7,0:3) :: a47_xxx
  ! parameters for the 8 intervals (T-dep):
  REAL(DP), DIMENSION(DIM58)   :: c0_1_1_xxx, c0_2_1_xxx
  REAL(DP), DIMENSION(DIM58)   :: c0_1_2_xxx, c0_2_2_xxx
  REAL(DP), DIMENSION(0:3)     :: c1_xxx
  REAL(DP), DIMENSION(2:7,0:3) :: c27_xxx

  ! NAMELIST /JVPP/:
  LOGICAL :: l_hardcoded
  CHARACTER(LEN=80)  :: lya_ir ! Lyman-alpha or IR addition
  INTEGER :: deg_tconst(8) ! T-const: degree of fit
  ! 1,2,3   for interval    1
  ! 1,2     for interval    2
  ! 1       for intervals 3-4
  ! 0,1,2,3 for intervals 5-8
  INTEGER :: deg_tdep(8) ! T-dep: degree of fit
  ! 1,2     for interval    1
  ! 1,2,3   for intervals 2-8
  INTEGER :: fj_corr
  CHARACTER(LEN=300) :: texrxn
  CHARACTER(LEN=20)  :: eqntag
  INTEGER, PARAMETER :: N_texnote_lines = 10
  CHARACTER(LEN=80), DIMENSION(N_texnote_lines) :: texnote

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE process_species(species, jp)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: species
    INTEGER, INTENT(IN) :: jp
    INTEGER :: zdeg_tdep, i
    CHARACTER(LEN=*), PARAMETER :: spaces = '             '

    WRITE(IO_LOG,*) HLINE1
    WRITE(IO_LOG,*) '***** SPECIES: '//species
    WRITE(IO_LOG,*) HLINE1
    WRITE(IO_LOG,*)

    ! write info to LaTeX file:
    IF (TRIM(texrxn)=="") THEN
      print *, '*** WARNING: empty texrxn for ', species
      texrxn = species
    ENDIF
    WRITE(IO_TEX,'(A)') &
      '\myhline '//TRIM(eqntag)//' & \chem{'//TRIM(texrxn)//'} & '
    IF (l_hardcoded) THEN
      WRITE(IO_TEX,'(A)') "Hardcoded from old JVAL code."
    ENDIF
    DO i = 1, N_texnote_lines
      IF (TRIM(texnote(i))/='') &
        WRITE(IO_TEX,'(A)') TRIM(texnote(i))
    ENDDO
    WRITE(IO_TEX,'(A)') '\\'

    WRITE(IO_CAL,'(A)') '    IF (lp(ip_'//species//'))'// &
      spaces(1:13-LEN(species))//'CALL jval_cal_'//species

    IF (l_hardcoded) THEN
      WRITE(*,'(A10,A)') species, '   hardcoded'
      ! copy hardcoded jval subroutine:
      WRITE(IO_CAT,'(A)') 'echo "  ! NOTE: this subroutine was not ' &
        //'calculated by JVPP but simply copied" >> '//INCFILE
      WRITE(IO_CAT,'(A)') 'echo "  ! from jval_cal_' &
        //species//'.f90" >> '//INCFILE
      WRITE(IO_CAT,'(A)') 'echo "" >> '//INCFILE
      WRITE(IO_CAT,'(A)') &
        'cat ../'//inputdir//'/hardcoded/jval_cal_'//species//'.f90 >> '//INCFILE
      WRITE(IO_CAT,'(A)') 'echo "" >> '//INCFILE
      WRITE(IO_CAT,'(A)') 'echo "  ! **'//HLINE1//'" >> '//INCFILE
      WRITE(IO_CAT,'(A)') 'echo "" >> '//INCFILE
      RETURN
    ENDIF

    IF (.NOT.l_species(jp)) THEN
      WRITE(*,'(A10,40X,3A)') species, '"', TRIM(lya_ir), '"'
      ! create jval code (only lya_ir) for current species:
      WRITE(IO_CAT,'(A)') 'cat jval_cal_'//species//'.f90 >> '//INCFILE
      OPEN(IO_F90, FILE=DIR_F90//'/jval_cal_'//species//'.f90', &
        STATUS='UNKNOWN')
      WRITE(IO_F90,'(A)') '  SUBROUTINE jval_cal_'//species
      WRITE(IO_F90,'(A)')
      WRITE(IO_F90,'(A)') '    ! '//TRIM(texrxn)
      WRITE(IO_F90,'(A)')
      WRITE(IO_F90,'(A)') '    INTEGER :: j, k'
      WRITE(IO_F90,'(A)') '    REAL    :: dj'
      ! calculations of sig_xxx:
      WRITE(IO_F90,'(A)')
      WRITE(IO_F90,'(4X,A)') 'jval_2d(ip_'//species//')%ptr(:,:) = 0.0_dp'
      WRITE(IO_F90,'(4X,A)') 'DO k = 1,klev'
      WRITE(IO_F90,'(4X,A)') '  DO j = 1,kproma_day'
      ! define dj and jval_gp(ip_xxx):
      WRITE(IO_F90,'(8X,A)',ADVANCE='NO') 'dj = 0.'
      IF (lya_ir/="") THEN
        WRITE(IO_F90,'(A)') ' &'
        WRITE(IO_F90,'(10X,A)',ADVANCE='NO') '+ '//TRIM(lya_ir)
      ENDIF
      WRITE(IO_F90,'(A)')
      WRITE(IO_F90,'(8X,A,I1,A)') 'jval_2d(ip_'//species// &
        ')%ptr(iu0(j),k) = REAL(MAX(0.0, dj*fj_corr(j,', fj_corr, ')),dp)'
      WRITE(IO_F90,'(6X,A)') 'ENDDO'
      WRITE(IO_F90,'(4X,A)') 'ENDDO'
      WRITE(IO_F90,'(A)')
      WRITE(IO_F90,'(2X,A)') 'END SUBROUTINE jval_cal_'//species
      WRITE(IO_F90,'(A)')
      WRITE(IO_F90,'(A)') '  ! **'//HLINE1
      WRITE(IO_F90,'(A)')
      CLOSE(IO_F90)
      RETURN
    ENDIF

    ! change deg_tdep if temperature dependence is not available:
    IF (.NOT.l_tdep(jp)) deg_tdep(:) = NO_TDEP

    ! T-const:
    DO ji=1, NBIN
      ji0 = ji - 1
      WRITE(ji_string,'(I1.1)') ji
      WRITE(IO_JTC,'(A)') '! SPECIES: '//species//',  INTERVAL: '//ji_string
      WRITE(IO_LOG,*) '*********************** SPECIES: '//species// &
        ' *** INTERVAL: '//ji_string//'***********'
      IF (deg_tconst(ji)>=0) THEN
        IF (ji==1) THEN ! Schumann-Runge bands
          maxv3=127
          maxv2=231
        ELSE ! interval 2-8
          maxv3=600
          maxv2=50
        ENDIF

        SELECT CASE(ji)
        CASE (1)   ; CALL process_interval_1(species)
        CASE (2)   ; CALL process_interval_2(species)
        CASE (3:4) ; CALL process_interval_3_4(species)
        CASE (5:8) ; CALL process_interval_5_8(species)
        END SELECT
      ENDIF
    ENDDO
    WRITE(IO_LOG,*)
    WRITE(IO_JTC,'(A)') 'GO newpage'
    WRITE(IO_JTC,'(A)') '! '//HLINE1

    ! T-dep:
    IF (l_tdep(jp)) THEN
      DO ji=1, NBIN
        zdeg_tdep = deg_tdep(ji)
        IF (zdeg_tdep>=0) THEN
          WRITE(ji_string,'(I1.1)') ji
          ji0 = ji - 1
          WRITE(IO_JTD,'(A)') '! SPECIES: '//species//',  INTERVAL: '//ji_string
          WRITE(IO_LOG,*) '*********************** SPECIES: '//species// &
            ' *** INTERVAL: '//ji_string//'* (T-dep) *'
          IF (ji==1) THEN ! Schumann-Runge bands:
            T_ref = 240._dp
            maxv3 =   4
            maxv2 = 231
          ELSE ! interval 2-8:
            T_ref = 250._dp
            maxv3 = 600
            maxv2 =   2
          ENDIF
          SELECT CASE(ji)
          CASE (1)   ; CALL process_tdep_interval_1(species)
          CASE (2)   ; CALL process_tdep_interval_2(species)
          CASE (3:8) ; CALL process_tdep_interval_3_8(species)
          END SELECT
        ENDIF
      ENDDO
      WRITE(IO_LOG,*)
      WRITE(IO_JTD,'(A)') 'GO newpage'
    ENDIF

    ! write info about degrees of fits:
    WRITE(*,'(A10,A)',ADVANCE='NO') species, '  '
    DO ji = 1, NBIN
      IF (deg_tconst(ji)<0) THEN
        SELECT CASE(deg_tconst(ji))
        CASE (FIT_PROBLEM) ; WRITE(*,'(A)',ADVANCE='NO') ' _'
        CASE (ALL_ZERO)    ; WRITE(*,'(A)',ADVANCE='NO') ' _'
        CASE DEFAULT       ; WRITE(*,'(A)',ADVANCE='NO') ' -'
        END SELECT
      ELSE
        WRITE(*,'(I2)',ADVANCE='NO') deg_tconst(ji)
      ENDIF
    ENDDO
    WRITE(*,'(A)',ADVANCE='NO') '   '
    DO ji = 1, NBIN
      IF (deg_tdep(ji)<0) THEN
        SELECT CASE(deg_tdep(ji))
        CASE (FIT_PROBLEM) ; WRITE(*,'(A)',ADVANCE='NO') ' _'
        CASE (ALL_ZERO)    ; WRITE(*,'(A)',ADVANCE='NO') ' _'
        CASE (NO_TDEP)     ; WRITE(*,'(A)',ADVANCE='NO') '  '
        CASE DEFAULT       ; WRITE(*,'(A)',ADVANCE='NO') ' -'
        END SELECT
      ELSE
        WRITE(*,'(I2)',ADVANCE='NO') deg_tdep(ji)
      ENDIF
    ENDDO
    WRITE(*,'(3A)') '   "', TRIM(lya_ir), '"'

    ! create jval code for current species:
    WRITE(IO_CAT,'(A)') 'cat jval_cal_'//species//'.f90 >> '//INCFILE
    OPEN(IO_F90, FILE=DIR_F90//'/jval_cal_'//species//'.f90', &
      STATUS='UNKNOWN')
    WRITE(IO_F90,'(A)') '  SUBROUTINE jval_cal_'//species
    WRITE(IO_F90,'(A)')
    WRITE(IO_F90,'(A)') '    ! '//TRIM(texrxn)
    WRITE(IO_F90,'(A)')
    WRITE(IO_F90,'(A)') '    INTEGER :: j, k'
    WRITE(IO_F90,'(A)') '    REAL    :: dj'
    WRITE(IO_F90,'(A)') '    REAL, DIMENSION(0:MAXWAV) :: sig_'//species//''
    WRITE(IO_F90,'(A)')

    ! T-const parameter definitions:
    WRITE(IO_F90,'(4X,A)') '! T-const parameters:'
    DO ji=1, NBIN
      ji0 = ji - 1
      WRITE(ji0_string,'(I1.1)') ji0
      IF (deg_tconst(ji)>=0) THEN
        SELECT CASE(ji)
        CASE (1)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: a0_1_1_'//species//'(dim58) = '
          CALL write_array(a0_1_1_xxx)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: a0_1_2_'//species//'(dim58) = '
          CALL write_array(a0_1_2_xxx)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: a0_2_1_'//species//'(dim58) = '
          CALL write_array(a0_2_1_xxx)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: a0_2_2_'//species//'(dim58) = '
          CALL write_array(a0_2_2_xxx)
          IF (deg_tconst(1)>=2) THEN
            WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
              'REAL, PARAMETER :: a0_3_1_'//species//'(dim58) = '
            CALL write_array(a0_3_1_xxx)
            WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
              'REAL, PARAMETER :: a0_3_2_'//species//'(dim58) = '
            CALL write_array(a0_3_2_xxx)
          ENDIF
          IF (deg_tconst(1)==3) THEN
            WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
              'REAL, PARAMETER :: a0_4_1_'//species//'(dim58) = '
            CALL write_array(a0_4_1_xxx)
            WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
              'REAL, PARAMETER :: a0_4_2_'//species//'(dim58) = '
            CALL write_array(a0_4_2_xxx)
          ENDIF
        CASE (2)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: a1_1_'//species//'(dim55) = '
          CALL write_array(a1_1_xxx)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: b1_1_'//species//'(dim55) = '
          CALL write_array(b1_1_xxx)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: a1_2_'//species//'(dim55) = '
          CALL write_array(a1_2_xxx)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: b1_2_'//species//'(dim55) = '
          CALL write_array(b1_2_xxx)
          IF (deg_tconst(2)==2) THEN
            WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
              'REAL, PARAMETER :: a1_3_'//species//'(dim55) = '
            CALL write_array(a1_3_xxx)
            WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
              'REAL, PARAMETER :: b1_3_'//species//'(dim55) = '
            CALL write_array(b1_3_xxx)
          ENDIF
        CASE (3:4)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: a'//ji0_string//'_'//species//'(dim55) = '
          CALL write_array(a23_xxx(ji0,:))
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: b'//ji0_string//'_'//species//'(dim55) = '
          CALL write_array(b23_xxx(ji0,:))
        CASE (5:8)
          WRITE(IO_F90,'(4X,A,I1.1,A)',ADVANCE='NO') 'REAL, PARAMETER :: a'// &
            ji0_string//'_'//species//'(', deg_tconst(ji)+1, ') = '
          CALL write_array(a47_xxx(ji0,0:deg_tconst(ji)))
        END SELECT
      ENDIF
    ENDDO

    ! T-dep parameter definitions:
    IF (l_tdep(jp)) THEN
      WRITE(IO_F90,'(4X,A)') '! T-dep parameters:'
    ENDIF
    DO ji=1, NBIN
      zdeg_tdep = deg_tdep(ji)
      ji0 = ji - 1
      WRITE(ji0_string,'(I1.1)') ji0
      IF (zdeg_tdep>=0) THEN
        SELECT CASE(ji)
        CASE (1)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: c0_1_1_'//species//'(dim58) = '
          CALL write_array(c0_1_1_xxx)
          WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
            'REAL, PARAMETER :: c0_1_2_'//species//'(dim58) = '
          CALL write_array(c0_1_2_xxx)
          IF (zdeg_tdep==2) THEN
            WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
              'REAL, PARAMETER :: c0_2_1_'//species//'(dim58) = '
            CALL write_array(c0_2_1_xxx)
            WRITE(IO_F90,'(4X,A)',ADVANCE='NO') &
              'REAL, PARAMETER :: c0_2_2_'//species//'(dim58) = '
            CALL write_array(c0_2_2_xxx)
          ENDIF
        CASE (2)
          WRITE(IO_F90,'(4X,A,I1,A)',ADVANCE='NO') 'REAL, PARAMETER :: c1_'// &
            species//'(', zdeg_tdep+1, ') = '
          CALL write_array(c1_xxx(0:zdeg_tdep))
        CASE (3:8)
          WRITE(IO_F90,'(4X,A,I1,A)',ADVANCE='NO') 'REAL, PARAMETER :: c'// &
            ji0_string//'_'//species//'(', zdeg_tdep+1, ') = '
          CALL write_array(c27_xxx(ji0,0:zdeg_tdep))
        END SELECT
      ENDIF
    ENDDO

    ! calculations of sig_xxx:
    WRITE(IO_F90,'(A)')
    WRITE(IO_F90,'(4X,A)') 'jval_2d(ip_'//species//')%ptr(:,:) = 0.0_dp'
    WRITE(IO_F90,'(4X,A)') 'DO k = 1,klev'
    WRITE(IO_F90,'(4X,A)') '  DO j = 1,kproma_day'

    ! 1:
    IF (deg_tconst(1)>=0) THEN
      WRITE(IO_F90,'(8X,A)') 'sig_'//species//'(0) = &'
      CALL write_sig_tdep_1(deg_tdep(1))
      SELECT CASE(deg_tconst(1))
      CASE (1)
        WRITE(IO_F90,'(8X,A)') '  p1(a0_1_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_1_2_'//species//'(i0(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  a0_2_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_2_2_'//species//'(i0(j,k)),v3_du1(j,k))'
      CASE (2)
        WRITE(IO_F90,'(8X,A)') '  p2(a0_1_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_1_2_'//species//'(i0(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  a0_2_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_2_2_'//species//'(i0(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  a0_3_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_3_2_'//species//'(i0(j,k)), v3_du1(j,k))'
      CASE (3)
        WRITE(IO_F90,'(8X,A)') '  p3(a0_1_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_1_2_'//species//'(i0(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  a0_2_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_2_2_'//species//'(i0(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  a0_3_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_3_2_'//species//'(i0(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  a0_4_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  a0_4_2_'//species//'(i0(j,k)), v3_du1(j,k))'
      END SELECT
    ENDIF
    ! 2:
    IF (deg_tconst(2)>=0) THEN
      WRITE(IO_F90,'(8X,A)') 'sig_'//species//'(1) = &'
      CALL write_sig_tdep_2(deg_tdep(2))
      SELECT CASE(deg_tconst(2))
      CASE (1)
        WRITE(IO_F90,'(8X,A)') '  (p1(b1_1_'//species//'(i1(j,k)),a1_1_'// &
          species//'(i1(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  v3_du1(j,k)) + &'
        WRITE(IO_F90,'(8X,A)') '  p1(b1_2_'//species//'(i1(j,k)),a1_2_'// &
          species//'(i1(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  v3_du1(j,k)) * v2s_m(j,k))'
      CASE (2)
        WRITE(IO_F90,'(8X,A)') '  (p1(b1_1_'//species//'(i1(j,k)),a1_1_'// &
          species//'(i1(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  v3_du1(j,k)) + &'
        WRITE(IO_F90,'(8X,A)') '  p1(b1_2_'//species//'(i1(j,k)),a1_2_'// &
          species//'(i1(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  v3_du1(j,k)) * v2s_m(j,k) + &'
        WRITE(IO_F90,'(8X,A)') '  p1(b1_3_'//species//'(i1(j,k)),a1_3_'// &
          species//'(i1(j,k)), &'
        WRITE(IO_F90,'(8X,A)') '  v3_du1(j,k)) * v2s_m(j,k)**2)'
      END SELECT
    ENDIF
    ! 3:
    WRITE(ji0_string,'(I1.1)') 2
    IF (deg_tconst(3)>=0) THEN
      WRITE(IO_F90,'(8X,A)') 'sig_'//species//'(2) = &'
      CALL write_sig_tdep_3_8(deg_tdep(3))
      WRITE(IO_F90,'(8X,A)') '  p1(b2_'//species// &
        '(i2(j,k)), a2_'//species//'(i2(j,k)), &'
      WRITE(IO_F90,'(8X,A)') '  v3_du1(j,k))'
    ENDIF
    ! 4:
    WRITE(ji0_string,'(I1.1)') 3
    IF (deg_tconst(4)>=0) THEN
      WRITE(IO_F90,'(8X,A)') 'sig_'//species//'(3) = &'
      CALL write_sig_tdep_3_8(deg_tdep(4))
      WRITE(IO_F90,'(8X,A)') '  p1(b3_'//species// &
        '(i3(j,k)), a3_'//species//'(i3(j,k)), &'
      WRITE(IO_F90,'(8X,A)') '  v3_du2(j,k))'
    ENDIF
    ! 5-8:
    DO ji=5,8
      ji0 = ji - 1
      WRITE(ji0_string,'(I1.1)') ji0
      IF (deg_tconst(ji)>=0) THEN
        WRITE(IO_F90,'(8X,A)') 'sig_'//species//'('//ji0_string//') = &'
        CALL write_sig_tdep_3_8(deg_tdep(ji))
        SELECT CASE(deg_tconst(ji))
        CASE (0)
          WRITE(IO_F90,'(8X,A)') '  a'//ji0_string//'_'//species//'(1)'
        CASE (1)
          WRITE(IO_F90,'(8X,A)') '  p1(a'//ji0_string//'_'//species// &
            '(1),a'//ji0_string//'_'//species//'(2),v3_du2(j,k))'
        CASE (2)
          WRITE(IO_F90,'(8X,A)') '  p2(a'//ji0_string//'_'//species//'(1),a' &
            //ji0_string//'_'//species//'(2),a'//ji0_string//'_'//species// &
            '(3),v3_du2(j,k))'
        CASE (3)
          WRITE(IO_F90,'(8X,A)') '  p3(a'//ji0_string//'_'//species//'(1),a' &
            //ji0_string//'_'//species//'(2),a'//ji0_string//'_'//species// &
            '(3), &'
          WRITE(IO_F90,'(8X,A)') '  a'//ji0_string//'_'//species//'(4),v3_du2(j,k))'
        END SELECT
      ENDIF
    ENDDO
    ! define dj and jval_gp(ip_xxx):
    WRITE(IO_F90,'(8X,A)',ADVANCE='NO') 'dj = 0.'
    IF (lya_ir/="") THEN
      WRITE(IO_F90,'(A)') ' &'
      WRITE(IO_F90,'(10X,A)',ADVANCE='NO') '+ '//TRIM(lya_ir)
    ENDIF
    DO ji=1, NBIN
      ji0 = ji - 1
      WRITE(ji0_string,'(I1.1)') ji0
      IF (deg_tconst(ji)>=0) THEN
        WRITE(IO_F90,'(A)') ' &'
        WRITE(IO_F90,'(10X,A)',ADVANCE='NO') &
          '+ sig_'//species//'('//ji0_string//') * fint(j,k,'//ji0_string//')'
      ENDIF
    ENDDO
    WRITE(IO_F90,'(A)')
    WRITE(IO_F90,'(A)')
    WRITE(IO_F90,'(8X,A,I1,A)') 'jval_2d(ip_'//species// &
      ')%ptr(iu0(j),k) = REAL(MAX(0.0, dj*fj_corr(j,', fj_corr, ')),dp)'
    WRITE(IO_F90,'(6X,A)') 'ENDDO'
    WRITE(IO_F90,'(4X,A)') 'ENDDO'
    WRITE(IO_F90,'(A)')
    WRITE(IO_F90,'(2X,A)') 'END SUBROUTINE jval_cal_'//species
    WRITE(IO_F90,'(A)')
    WRITE(IO_F90,'(A)') '  ! **'//HLINE1
    WRITE(IO_F90,'(A)')
    CLOSE(IO_F90)

    CONTAINS

      ! ----------------------

      SUBROUTINE write_sig_tdep_1(zdeg_tdep)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: zdeg_tdep ! deg_tdep as a local variable
        SELECT CASE(zdeg_tdep)
        CASE (1)
          WRITE(IO_F90,'(10X,A)') 'p1(1.,c0_1_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
          WRITE(IO_F90,'(10X,A)') 'c0_1_2_'//species//'(i0(j,k)),tnorm_sr(j,k)) * &'
        CASE (2)
          WRITE(IO_F90,'(10X,A)') 'p2(1.,c0_1_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
          WRITE(IO_F90,'(10X,A)') 'c0_1_2_'//species//'(i0(j,k)), &'
          WRITE(IO_F90,'(10X,A)') 'c0_2_1_'//species//'(i0(j,k))*dlv2(j,k) + &'
          WRITE(IO_F90,'(10X,A)') 'c0_2_2_'//species//'(i0(j,k)), tnorm_sr(j,k)) * &'
        END SELECT
      END SUBROUTINE write_sig_tdep_1

      ! ----------------------

      SUBROUTINE write_sig_tdep_2(zdeg_tdep)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: zdeg_tdep ! deg_tdep as a local variable
        SELECT CASE(zdeg_tdep)
        CASE (1)
          WRITE(IO_F90,'(10X,A)') 'p1(c1_'//species//'(1),c1_'//species// &
            '(2),tnorm(j,k)) * &'
        CASE (2)
          WRITE(IO_F90,'(10X,A)') 'p2(c1_'//species//'(1),c1_'//species//'(2), &'
          WRITE(IO_F90,'(10X,A)') 'c1_'//species//'(3),tnorm(j,k)) * &'
        CASE (3)
          WRITE(IO_F90,'(10X,A)') 'p3(c1_'//species//'(1),c1_'//species//'(2), &'
          WRITE(IO_F90,'(10X,A)') 'c1_'//species//'(3),c1_'//species// &
            '(4),tnorm(j,k)) * &'
        END SELECT
      END SUBROUTINE write_sig_tdep_2

      ! ----------------------

      SUBROUTINE write_sig_tdep_3_8(zdeg_tdep)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: zdeg_tdep ! deg_tdep as a local variable
        SELECT CASE(zdeg_tdep)
        CASE (1)
          WRITE(IO_F90,'(10X,A)') 'p1(c'//ji0_string//'_'//species//'(1),c'// &
            ji0_string//'_'//species//'(2),tnorm(j,k)) * &'
        CASE (2)
          WRITE(IO_F90,'(10X,A)') 'p2(c'//ji0_string//'_'//species//'(1),c'// &
            ji0_string//'_'//species//'(2),c'//ji0_string//'_'//species// &
            '(3),tnorm(j,k))* &'
        CASE (3)
          WRITE(IO_F90,'(10X,A)') 'p3(c'//ji0_string//'_'//species//'(1),c'// &
            ji0_string//'_'//species//'(2),c'//ji0_string//'_'//species// &
            '(3),c'//ji0_string//'_'//species//'(4),tnorm(j,k))* &'
        END SELECT
      END SUBROUTINE write_sig_tdep_3_8

      ! ----------------------

  END SUBROUTINE process_species

  ! --------------------------------------------------------------------------

  SUBROUTINE write_array(array)

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER :: i, n

    ! set small numbers to zero:
    WHERE (ABS(array(:)) < 1.0E-38_dp) array(:) = 0._dp

    n = SIZE(array)
    IF (n>4) THEN
      WRITE(IO_F90,'(A)') '(/ &'
      WRITE(IO_F90,'(A)',ADVANCE='NO') '      '
    ELSE
      WRITE(IO_F90,'(A)',ADVANCE='NO') '(/ '
    ENDIF
    DO i=1,n
      WRITE(IO_F90,'(1P,6E11.4)',ADVANCE='NO') array(i)
      IF (i<n) WRITE(IO_F90,'(A)',ADVANCE='NO') ', ' ! no comma after last value
      IF (MODULO(i,6)==0) THEN ! line break after 6 numbers
        WRITE(IO_F90,'(A)') '&'
        WRITE(IO_F90,'(A)',ADVANCE='NO') '      '
      ENDIF
    ENDDO
    WRITE(IO_F90,'(A)') ' /)'

  END SUBROUTINE write_array

  ! --------------------------------------------------------------------------

  SUBROUTINE read_file_1d(filename, nheader, DATA)

    IMPLICIT NONE

    CHARACTER(LEN=*),       INTENT(IN)  :: filename
    INTEGER,                INTENT(OUT) :: nheader
    REAL(DP), DIMENSION(:), INTENT(OUT) :: DATA

    WRITE(IO_LOG,*) 'Reading file: '//TRIM(filename)
    OPEN(IO_TMP, FILE=TRIM(filename), STATUS='OLD')
    nheader = 0
    DO
      READ(IO_TMP,'(A)') header
      IF (header(1:1)/='#') EXIT ! exit do loop if not a header line
      nheader = nheader + 1 ! count header lines
    ENDDO
    BACKSPACE(IO_TMP) ! back to previous line
    READ(IO_TMP,*) DATA
    CLOSE(IO_TMP)

  END SUBROUTINE read_file_1d

  ! --------------------------------------------------------------------------

  SUBROUTINE read_file_2d(filename, nheader, DATA)

    IMPLICIT NONE

    CHARACTER(LEN=*),         INTENT(IN)  :: filename
    INTEGER,                  INTENT(OUT) :: nheader
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: DATA

    INTEGER :: k2, k3

    WRITE(IO_LOG,*) 'Reading file: '//TRIM(filename)
    OPEN(IO_TMP, FILE=TRIM(filename), STATUS='OLD')
    nheader = 0
    DO
      READ(IO_TMP,'(A)') header
      IF (header(1:1)/='#') EXIT ! exit do loop if not a header line
      nheader = nheader + 1 ! count header lines
    ENDDO
    BACKSPACE(IO_TMP) ! back to previous line
    DO k2=1,SIZE(DATA,1)
      DO k3=1,SIZE(DATA,2)
        READ(IO_TMP,*) DATA(k2,k3)
      ENDDO
    ENDDO
    CLOSE(IO_TMP)

  END SUBROUTINE read_file_2d

  ! --------------------------------------------------------------------------

  SUBROUTINE process_interval_1(species)

    USE messy_main_math_lsq, ONLY: poly_fit

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: species

    REAL(DP) :: r2
    INTEGER  :: nheader_v2, nheader_v3_du, nheader_xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: v2, v3_du
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: e, recalc
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: c0, c1, c2, c3
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: c0_recalc, c1_recalc, c2_recalc, c3_recalc
    REAL(DP) :: coeff(0:3)
    REAL(DP) :: v3_lim, mini, maxi
    INTEGER :: k2, k3_lim, k3, k, zdeg_tconst

    ALLOCATE(recalc(maxv2,maxv3), e(maxv2,maxv3))
    ALLOCATE(v2(maxv2), v3_du(maxv3), xxx_eff(maxv2,maxv3))
    ALLOCATE(c0(maxv2), c1(maxv2), c2(maxv2), c3(maxv2))
    ALLOCATE(c0_recalc(maxv2), c1_recalc(maxv2), c2_recalc(maxv2), c3_recalc(maxv2))

    ! read input files:
    v3_du_filename = DIR_EFF//'/v3_du.ef'//ji_string
    CALL read_file_1d(v3_du_filename, nheader_v3_du, v3_du)
    v2_filename = DIR_EFF//'/v2.ef'//ji_string
    CALL read_file_1d(v2_filename, nheader_v2, v2)
    xxx_eff_filename = DIR_EFF//'/'//species//'.ef'//ji_string
    CALL read_file_2d(xxx_eff_filename, nheader_xxx_eff, xxx_eff)
    WRITE(IO_LOG,*) 'xxx_eff: ', MINVAL(xxx_eff), '...', MAXVAL(xxx_eff)

    ! calculate parameters:
    zdeg_tconst = deg_tconst(1)
    DO k2 = 1,maxv2
      ! IDL syntax was: v3_lim = (40.*LOG(v2(k2))-1660.) > 200.
      v3_lim = MAX(40._dp*LOG(v2(k2))-1660._dp,200._dp)
      k3_lim = maxv3
      DO k3 = 1,maxv3
        IF (v3_du(k3)<=v3_lim) k3_lim = k3
      ENDDO
      mini = MINVAL(xxx_eff(k2,1:k3_lim))
      maxi = MAXVAL(xxx_eff(k2,1:k3_lim))
      WRITE(IO_LOG,'(1P,A,I3,A,I3,A,E11.4,A,E11.4,A)', ADVANCE='NO') &
        ' xxx_eff(',k2,',1:',k3_lim,'): ', mini, ' ...', maxi, ', '
      IF (maxi-mini>TINY(0._dp)) THEN
        CALL poly_fit(v3_du(1:k3_lim), xxx_eff(k2,1:k3_lim), &
          zdeg_tconst, coeff(0:zdeg_tconst), r2, IO_LOG)
      ELSE
        coeff(:) = 0._dp
        coeff(0) = maxi
        r2 = 1._dp
        WRITE(IO_LOG,'(A,G20.12)') ' R^2 = ', r2
      ENDIF
      ! Do not use this interval if fit didn't work:
      IF (ABS(r2)<=TINY(0._dp)) deg_tconst(1) = FIT_PROBLEM
      ! Do not use this interval if all values are zero:
      IF (ABS(maxi)<=TINY(0._dp)) deg_tconst(1) = ALL_ZERO

      c0(k2) = coeff(0)
      c1(k2) = coeff(1)
      IF (zdeg_tconst>=2) c2(k2) = coeff(2)
      IF (zdeg_tconst>=3) c3(k2) = coeff(3)
    ENDDO

    DO k = 1, DIM58-1
      a0_1_1_xxx(k)  = (c0(4*k-3)-c0(4*k-3+4)) / (LOG(v2(4*k-3))-LOG(v2(4*k-3+4)))
      a0_1_2_xxx(k)  = c0(4*k-3) - a0_1_1_xxx(k)*LOG(v2(4*k-3))
      a0_2_1_xxx(k)  = (c1(4*k-3)-c1(4*k-3+4)) / (LOG(v2(4*k-3))-LOG(v2(4*k-3+4)))
      a0_2_2_xxx(k)  = c1(4*k-3) - a0_2_1_xxx(k)*LOG(v2(4*k-3))
      c0_recalc(4*k-3) = a0_1_1_xxx(k)*LOG(v2(4*k-3)) + a0_1_2_xxx(k)
      c0_recalc(4*k-2) = a0_1_1_xxx(k)*LOG(v2(4*k-2)) + a0_1_2_xxx(k)
      c0_recalc(4*k-1) = a0_1_1_xxx(k)*LOG(v2(4*k-1)) + a0_1_2_xxx(k)
      c0_recalc(4*k)   = a0_1_1_xxx(k)*LOG(v2(4*k))   + a0_1_2_xxx(k)
      c1_recalc(4*k-3) = a0_2_1_xxx(k)*LOG(v2(4*k-3)) + a0_2_2_xxx(k)
      c1_recalc(4*k-2) = a0_2_1_xxx(k)*LOG(v2(4*k-2)) + a0_2_2_xxx(k)
      c1_recalc(4*k-1) = a0_2_1_xxx(k)*LOG(v2(4*k-1)) + a0_2_2_xxx(k)
      c1_recalc(4*k)   = a0_2_1_xxx(k)*LOG(v2(4*k))   + a0_2_2_xxx(k)
      IF (zdeg_tconst>=2) THEN
        a0_3_1_xxx(k)  = (c2(4*k-3)-c2(4*k-3+4)) / (LOG(v2(4*k-3))-LOG(v2(4*k-3+4)))
        a0_3_2_xxx(k)  = c2(4*k-3) - a0_3_1_xxx(k)*LOG(v2(4*k-3))
        c2_recalc(4*k-3) = a0_3_1_xxx(k)*LOG(v2(4*k-3)) + a0_3_2_xxx(k)
        c2_recalc(4*k-2) = a0_3_1_xxx(k)*LOG(v2(4*k-2)) + a0_3_2_xxx(k)
        c2_recalc(4*k-1) = a0_3_1_xxx(k)*LOG(v2(4*k-1)) + a0_3_2_xxx(k)
        c2_recalc(4*k)   = a0_3_1_xxx(k)*LOG(v2(4*k))   + a0_3_2_xxx(k)
      ENDIF
      IF (zdeg_tconst>=3) THEN
        a0_4_1_xxx(k)  = (c3(4*k-3)-c3(4*k-3+4)) / (LOG(v2(4*k-3))-LOG(v2(4*k-3+4)))
        a0_4_2_xxx(k)  = c3(4*k-3) - a0_4_1_xxx(k)*LOG(v2(4*k-3))
        c3_recalc(4*k-3) = a0_4_1_xxx(k)*LOG(v2(4*k-3)) + a0_4_2_xxx(k)
        c3_recalc(4*k-2) = a0_4_1_xxx(k)*LOG(v2(4*k-2)) + a0_4_2_xxx(k)
        c3_recalc(4*k-1) = a0_4_1_xxx(k)*LOG(v2(4*k-1)) + a0_4_2_xxx(k)
        c3_recalc(4*k)   = a0_4_1_xxx(k)*LOG(v2(4*k))   + a0_4_2_xxx(k)
      ENDIF
    ENDDO
    k = DIM58
    a0_1_1_xxx(k) = (c0(229)-c0(231)) / (LOG(v2(229))-LOG(v2(231)))
    a0_1_2_xxx(k) = c0(229) - a0_1_1_xxx(DIM58)*LOG(v2(229))
    a0_2_1_xxx(k) = (c1(229)-c1(231)) / (LOG(v2(229))-LOG(v2(231)))
    a0_2_2_xxx(k) = c1(229) - a0_2_1_xxx(DIM58)*LOG(v2(229))
    c0_recalc(229)  = a0_1_1_xxx(k)*LOG(v2(229)) + a0_1_2_xxx(k)
    c0_recalc(230)  = a0_1_1_xxx(k)*LOG(v2(230)) + a0_1_2_xxx(k)
    c0_recalc(231)  = a0_1_1_xxx(k)*LOG(v2(231)) + a0_1_2_xxx(k)
    c1_recalc(229)  = a0_2_1_xxx(k)*LOG(v2(229)) + a0_2_2_xxx(k)
    c1_recalc(230)  = a0_2_1_xxx(k)*LOG(v2(230)) + a0_2_2_xxx(k)
    c1_recalc(231)  = a0_2_1_xxx(k)*LOG(v2(231)) + a0_2_2_xxx(k)
    IF (zdeg_tconst>=2) THEN
      a0_3_1_xxx(DIM58) = (c2(229)-c2(231)) / (LOG(v2(229))-LOG(v2(231)))
      a0_3_2_xxx(DIM58) = c2(229) - a0_3_1_xxx(DIM58)*LOG(v2(229))
      c2_recalc(229)      = a0_3_1_xxx(DIM58)*LOG(v2(229)) + a0_3_2_xxx(DIM58)
      c2_recalc(230)      = a0_3_1_xxx(DIM58)*LOG(v2(230)) + a0_3_2_xxx(DIM58)
      c2_recalc(231)      = a0_3_1_xxx(DIM58)*LOG(v2(231)) + a0_3_2_xxx(DIM58)
    ENDIF
    IF (zdeg_tconst>=3) THEN
      a0_4_1_xxx(DIM58) = (c3(229)-c3(231)) / (LOG(v2(229))-LOG(v2(231)))
      a0_4_2_xxx(DIM58) = c3(229) - a0_4_1_xxx(DIM58)*LOG(v2(229))
      c3_recalc(229)      = a0_4_1_xxx(DIM58)*LOG(v2(229)) + a0_4_2_xxx(DIM58)
      c3_recalc(230)      = a0_4_1_xxx(DIM58)*LOG(v2(230)) + a0_4_2_xxx(DIM58)
      c3_recalc(231)      = a0_4_1_xxx(DIM58)*LOG(v2(231)) + a0_4_2_xxx(DIM58)
    ENDIF

    ! calculate error (in percent):
    DO k3 = 1,maxv3
      recalc(:,k3) = c0_recalc(:) + c1_recalc(:) * v3_du(k3)
      IF (zdeg_tconst>=2) recalc(:,k3) = recalc(:,k3) + c2_recalc(:)*v3_du(k3)**2
      IF (zdeg_tconst>=3) recalc(:,k3) = recalc(:,k3) + c3_recalc(:)*v3_du(k3)**3
    ENDDO
    IF (deg_tconst(1)>=0) THEN
      e = (recalc/xxx_eff-1._dp) * 100._dp
    ELSE
      e = 0.0_dp ! dummy value
    ENDIF
    WRITE(IO_LOG,*) 'Maximum error: ', MAXVAL(ABS(E)), ' %'
    error_filename = DIR_JNL//'/'//species//'.er'//ji_string
    recalc_filename = DIR_JNL//'/'//species//'.rc'//ji_string
    OPEN(IO_ERR, FILE=error_filename, STATUS='UNKNOWN')
    OPEN(IO_REC, FILE=recalc_filename, STATUS='UNKNOWN')
    WRITE(IO_ERR,'(A)') '# '//species//' error in %'
    WRITE(IO_REC,'(A)') '# '//species//' recalc'
    DO K2 = 1, maxv2
      DO K3 = 1, maxv3
        WRITE(IO_ERR,'(E12.4)') E(K2,K3)
        WRITE(IO_REC,'(E12.4)') recalc(K2,K3)
      ENDDO
    ENDDO
    CLOSE(IO_ERR)
    CLOSE(IO_REC)

    ! write to ferret jnl file:
    WRITE(IO_JTC,'(A)') 'CANCEL DATA/ALL'
    WRITE(IO_JTC,'(A,I0,A)') &
      'DEFINE AXIS/UNITS="1...maxv2"/X=1:', maxv2, ':1 v2axis'
    WRITE(IO_JTC,'(A,I0,A)') &
      'DEFINE AXIS/UNITS="1...maxv3"/Y=1:', maxv3, ':1 v3axis'
    WRITE(IO_JTC,'(A)') 'DEFINE GRID/X=v3axis/Y=v2axis mygrid'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_xxx_eff, &
      '/VAR="'//species//'"/GRID=mygrid "../'//TRIM(xxx_eff_filename)//'"'
    WRITE(IO_JTC,'(A)') 'FILE/SKIP=1/VAR=recalc/GRID=mygrid "../'// &
      TRIM(recalc_filename)//'"'
    WRITE(IO_JTC,'(A)') 'FILE/SKIP=1/VAR=error/GRID=mygrid "../'// &
      TRIM(error_filename)//'"'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_v2, &
      '/VAR="v2" "../'//TRIM(v2_filename)//'"'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_v3_du, &
      '/VAR="v3_du" "../'//TRIM(v3_du_filename)//'"'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A,I0,A)') 'SHADE/SET/TITLE="'//species// &
      ' (int. 1)" '//species//'[d=1]'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL SHADE'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'SHADE/SET/TITLE="'//species// &
      ' (int. 1) recalc" recalc[d=2]'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL SHADE'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'SHADE/SET/TITLE="'//species// &
      ' (int. 1) error" error[d=3]'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL SHADE'
    WRITE(IO_JTC,'(A)')

    DEALLOCATE(v2, v3_du, xxx_eff)
    DEALLOCATE(recalc, e)
    DEALLOCATE(c0, c1, c2, c3)
    DEALLOCATE(c0_recalc, c1_recalc, c2_recalc, c3_recalc)

  END SUBROUTINE process_interval_1

  ! --------------------------------------------------------------------------

  SUBROUTINE process_interval_2(species)

    USE messy_main_math_lsq, ONLY: poly_fit

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: species

    REAL(DP) :: r2
    INTEGER  :: i, nheader_v2, nheader_v3_du, nheader_xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: v2, v3_du
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: recalc, e
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: recalc0, recalc1, recalc2, c0, c1, c2
    REAL(DP) :: coeff(0:2)
    REAL(DP) :: v2_lim, mini, maxi
    INTEGER :: k2_lim, k2, kg, kp, km, k3, zdeg_tconst

    ALLOCATE(v2(maxv2), v3_du(maxv3), xxx_eff(maxv2,maxv3))
    ALLOCATE(recalc0(maxv3), recalc1(maxv3), recalc2(maxv3))
    ALLOCATE(c0(maxv3), c1(maxv3), c2(maxv3))
    ALLOCATE(recalc(maxv2,maxv3), e(maxv2,maxv3))

    ! read input files:
    v3_du_filename = DIR_EFF//'/v3_du.ef'//ji_string
    CALL read_file_1d(v3_du_filename, nheader_v3_du, v3_du)
    v2_filename = DIR_EFF//'/v2.ef'//ji_string
    CALL read_file_1d(v2_filename, nheader_v2, v2)
    xxx_eff_filename = DIR_EFF//'/'//species//'.ef'//ji_string
    CALL read_file_2d(xxx_eff_filename, nheader_xxx_eff, xxx_eff)
    WRITE(IO_LOG,*) 'xxx_eff: ', MINVAL(xxx_eff), '...', MAXVAL(xxx_eff)

    ! calculate parameters:
    zdeg_tconst = deg_tconst(2)
    v2 = v2 * 3.767e-20_dp * 0.01_dp
    DO k3 = 1, maxv3
      v2_lim = 2._dp*v3_du(k3) + 100._dp
      k2_lim = maxv2
      DO k2 = 1, maxv2
        IF(v2(k2)<=v2_lim) k2_lim = k2
      ENDDO
      mini = MINVAL(xxx_eff(1:k2_lim,k3))
      maxi = MAXVAL(xxx_eff(1:k2_lim,k3))
      WRITE(IO_LOG,'(1P,A,I3,A,I3,A,E11.4,A,E11.4,A)', ADVANCE='NO') &
        ' xxx_eff(',k3,',1:',k2_lim,'): ', mini, ' ...', maxi, ', '
      IF (maxi-mini>TINY(0._dp)) THEN
        CALL poly_fit(v2(1:k2_lim), xxx_eff(1:k2_lim,k3), &
          zdeg_tconst, coeff(0:zdeg_tconst), r2, IO_LOG)
      ELSE
        coeff(:) = 0._dp
        coeff(0) = mini
        r2 = 1._dp
        WRITE(IO_LOG,'(A,G20.12)') ' R^2 = ', r2
      ENDIF
      ! Do not use this interval if fit didn't work:
      IF (ABS(r2)<=TINY(0._dp)) deg_tconst(2) = FIT_PROBLEM
      ! Do not use this interval if all values are zero:
      IF (ABS(maxi)<=TINY(0._dp)) deg_tconst(2) = ALL_ZERO

      c0(k3) = coeff(0)
      c1(k3) = coeff(1)
      IF (zdeg_tconst>=2) c2(k3) = coeff(2)
    ENDDO

    DO kg = 1, 40
      a1_1_xxx(kg) = (c0(kg*5-4)-c0(kg*5+1)) / (v3_du(kg*5-4)-v3_du(kg*5+1))
      b1_1_xxx(kg) = c0(kg*5-4) - a1_1_xxx(kg)*v3_du(kg*5-4)
      a1_2_xxx(kg) = (c1(kg*5-4)-c1(kg*5+1)) / (v3_du(kg*5-4)-v3_du(kg*5+1))
      b1_2_xxx(kg) = c1(kg*5-4) - a1_2_xxx(kg)*v3_du(kg*5-4)
      IF (zdeg_tconst>=2) THEN
      a1_3_xxx(kg) = (c2(kg*5-4)-c2(kg*5+1)) / (v3_du(kg*5-4)-v3_du(kg*5+1))
      b1_3_xxx(kg) = c2(kg*5-4) - a1_3_xxx(kg)*v3_du(kg*5-4)
      ENDIF
      DO k3 = 1, 10
        recalc0(kg*5+k3-5) = a1_1_xxx(kg)*v3_du(kg*5+k3-5) + b1_1_xxx(kg)
        recalc1(kg*5+k3-5) = a1_2_xxx(kg)*v3_du(kg*5+k3-5) + b1_2_xxx(kg)
        IF (zdeg_tconst>=2) &
        recalc2(kg*5+k3-5) = a1_3_xxx(kg)*v3_du(kg*5+k3-5) + b1_3_xxx(kg)
      ENDDO
    ENDDO

    DO i = 0, 14
      kg = i+40
      kp = i*25+225
      km = i*25+200
      a1_1_xxx(kg+1) = (c0(km+1)-c0(kp+1)) / (v3_du(km+1)-v3_du(kp+1))
      b1_1_xxx(kg+1) = c0(km+1) - a1_1_xxx(kg+1)*v3_du(km+1)
      a1_2_xxx(kg+1) = (c1(km+1)-c1(kp+1)) / (v3_du(km+1)-v3_du(kp+1))
      b1_2_xxx(kg+1) = c1(km+1) - a1_2_xxx(kg+1)*v3_du(km+1)
      IF (zdeg_tconst>=2) THEN
      a1_3_xxx(kg+1) = (c2(km+1)-c2(kp+1)) / (v3_du(km+1)-v3_du(kp+1))
      b1_3_xxx(kg+1) =  c2(km+1) - a1_3_xxx(kg+1)*v3_du(km+1)
      ENDIF
      DO k3 = 0, 24
        recalc0(km+k3+1) = a1_1_xxx(kg+1)*v3_du(km+k3+1) + b1_1_xxx(kg+1)
        recalc1(km+k3+1) = a1_2_xxx(kg+1)*v3_du(km+k3+1) + b1_2_xxx(kg+1)
        IF (zdeg_tconst>=2) &
        recalc2(km+k3+1) = a1_3_xxx(kg+1)*v3_du(km+k3+1) + b1_3_xxx(kg+1)
      ENDDO
    ENDDO

    km = 575
    DO k3 = 0, 24
      recalc0(km+k3+1) = a1_1_xxx(DIM55)*v3_du(km+k3+1) + b1_1_xxx(DIM55)
      recalc1(km+k3+1) = a1_2_xxx(DIM55)*v3_du(km+k3+1) + b1_2_xxx(DIM55)
      IF (zdeg_tconst>=2) &
      recalc2(km+k3+1) = a1_3_xxx(DIM55)*v3_du(km+k3+1) + b1_3_xxx(DIM55)
    ENDDO

    DO k2 = 1, maxv2
      DO k3 = 1, maxv3
        recalc(k2,k3) = recalc0(k3) + recalc1(k3)*v2(k2)
        IF (zdeg_tconst>=2) recalc(k2,k3) = &
          recalc0(k3) + recalc1(k3)*v2(k2) + recalc2(k3)*v2(k2)**2
      ENDDO
    ENDDO

    IF (deg_tconst(2)>=0) THEN
      e = (recalc/xxx_eff-1._dp) * 100._dp
    ELSE
      e = 0.0_dp ! dummy value
    ENDIF
    WRITE(IO_LOG,*) 'Maximum error: ', MAXVAL(ABS(E)), ' %'
    error_filename = DIR_JNL//'/'//species//'.er'//ji_string
    recalc_filename = DIR_JNL//'/'//species//'.rc'//ji_string
    OPEN(IO_ERR, FILE=error_filename, STATUS='UNKNOWN')
    OPEN(IO_REC, FILE=recalc_filename, STATUS='UNKNOWN')
    WRITE(IO_ERR,'(A)') '# '//species//' error in %'
    WRITE(IO_REC,'(A)') '# '//species//' recalc'
    DO K2 = 1, maxv2
      DO K3 = 1, maxv3
        WRITE(IO_ERR,'(E12.4)') E(K2,K3)
        WRITE(IO_REC,'(E12.4)') recalc(K2,K3)
      ENDDO
    ENDDO
    CLOSE(IO_ERR)
    CLOSE(IO_REC)

    ! write to ferret jnl file:
    WRITE(IO_JTC,'(A)') 'CANCEL DATA/ALL'
    WRITE(IO_JTC,'(A,I0,A)') &
      'DEFINE AXIS/UNITS="1...maxv2"/X=1:', maxv2, ':1 v2axis'
    WRITE(IO_JTC,'(A,I0,A)') &
      'DEFINE AXIS/UNITS="1...maxv3"/Y=1:', maxv3, ':1 v3axis'
    WRITE(IO_JTC,'(A)') 'DEFINE GRID/X=v3axis/Y=v2axis mygrid'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_xxx_eff, &
      '/VAR="'//species//'"/GRID=mygrid "../'//TRIM(xxx_eff_filename)//'"'
    WRITE(IO_JTC,'(A)') 'FILE/SKIP=1/VAR=recalc/GRID=mygrid "../'// &
      TRIM(recalc_filename)//'"'
    WRITE(IO_JTC,'(A)') 'FILE/SKIP=1/VAR=error/GRID=mygrid "../'// &
      TRIM(error_filename)//'"'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_v2, &
      '/VAR="v2" "../'//TRIM(v2_filename)//'"'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_v3_du, &
      '/VAR="v3_du" "../'//TRIM(v3_du_filename)//'"'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'SHADE/SET/TITLE="'//species// &
      ' (int. 2)" '//species//'[d=1]'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL SHADE'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'SHADE/SET/TITLE="'//species// &
      ' (int. 2) recalc" recalc[d=2]'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL SHADE'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'SHADE/SET/TITLE="'//species// &
      ' (int. 2) error" error[d=3]'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL SHADE'
    WRITE(IO_JTC,'(A)')

    DEALLOCATE(v2, v3_du, xxx_eff)
    DEALLOCATE(recalc0, recalc1, recalc2)
    DEALLOCATE(c0, c1, c2)
    DEALLOCATE(recalc, e)

  END SUBROUTINE process_interval_2

  ! --------------------------------------------------------------------------

  SUBROUTINE process_interval_3_4(species)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: species

    INTEGER  :: i, nheader_v3_du, nheader_xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: v3_du, xxx_eff, recalc, E
    INTEGER :: K, Kg, KP, KM, K3

    ALLOCATE(v3_du(maxv3), xxx_eff(maxv3), recalc(maxv3), E(maxv3))

    ! read input files:
    v3_du_filename = DIR_EFF//'/v3_du.ef'//ji_string
    CALL read_file_1d(v3_du_filename, nheader_v3_du, v3_du)
    xxx_eff_filename = DIR_EFF//'/'//species//'.ef'//ji_string
    CALL read_file_1d(xxx_eff_filename, nheader_xxx_eff, xxx_eff)
    WRITE(IO_LOG,*) 'xxx_eff: ', MINVAL(xxx_eff), '...', MAXVAL(xxx_eff)

    ! Do not use this interval if all values are zero:
    IF (ABS(MAXVAL(xxx_eff))<=TINY(0._dp)) deg_tconst(ji) = ALL_ZERO

    ! calculate parameters:
    DO K=1,40
      a23_xxx(ji0,K) = (xxx_eff(5*K-4)-xxx_eff(5*K+1)) / &
        (v3_du(5*K-4)-v3_du(5*K+1))
      b23_xxx(ji0,K) =  xxx_eff(5*K-4) - a23_xxx(ji0,K)*v3_du(5*K-4)
      DO I=0,4
        recalc(5*K+I-4) = a23_xxx(ji0,K) * v3_du(5*K+I-4) + b23_xxx(ji0,K)
      ENDDO
    ENDDO
    DO I=0,14
      Kg = I+41
      KP = I*25+226
      KM = I*25+201
      a23_xxx(ji0,Kg) = (xxx_eff(KM)-xxx_eff(KP)) / (v3_du(KM)-v3_du(KP))
      b23_xxx(ji0,Kg) = xxx_eff(KM) - a23_xxx(ji0,Kg)*v3_du(KM)
      DO K=0,24
        recalc(KM+K) = a23_xxx(ji0,Kg)*v3_du(KM+K) + b23_xxx(ji0,Kg)
      ENDDO
    ENDDO
    KM = 575
    DO K=1,25
      recalc(KM+K) = a23_xxx(ji0,DIM55)*v3_du(KM+K) + b23_xxx(ji0,DIM55)
    ENDDO

    ! calculate error (in percent):
    DO k3 = 1,maxv3
      IF (ABS(xxx_eff(k3))>TINY(0._dp)) THEN
        E(k3) = (recalc(k3)/xxx_eff(k3)-1._dp) * 100._dp
      ELSE
        E(k3) = 0._dp
      ENDIF
    ENDDO
    WRITE(IO_LOG,*) 'Maximum error: ', MAXVAL(ABS(E)), ' %'
    error_filename = DIR_JNL//'/'//species//'.er'//ji_string
    recalc_filename = DIR_JNL//'/'//species//'.rc'//ji_string
    OPEN(IO_ERR, FILE=error_filename, STATUS='UNKNOWN')
    OPEN(IO_REC, FILE=recalc_filename, STATUS='UNKNOWN')
    WRITE(IO_ERR,'(A)') '# '//species//' error in %'
    WRITE(IO_REC,'(A)') '# '//species//' recalc'
    DO K3 = 1, maxv3
      WRITE(IO_ERR,'(E12.4)') E(K3)
      WRITE(IO_REC,'(E12.4)') recalc(K3)
    ENDDO
    CLOSE(IO_ERR)
    CLOSE(IO_REC)

    ! write to ferret jnl file:
    WRITE(IO_JTC,'(A)') 'CANCEL DATA/ALL'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_xxx_eff, &
      '/VAR="'//species//'" "../'//TRIM(xxx_eff_filename)//'"'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=1/VAR="recalc" "../'// &
      TRIM(recalc_filename)//'"'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=1/VAR="error" "../'// &
      TRIM(error_filename)//'"'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_v3_du, &
      '/VAR="v3_du" "../'//TRIM(v3_du_filename)//'"'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'PLOT/SET/LINE/TITLE="'//species// &
      ' vs v3 (int. '//ji_string//')"/VS/COLOR=2 v3_du, '// &
      species//'[d=1]'
    WRITE(IO_JTC,'(A)') 'PPL XLAB "" ; PPL YLAB "1/cm2"'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL PLOT'
    WRITE(IO_JTC,'(A)') &
      'PLOT/LINE/OVER/DASH=(0.1,0.1,0.1,0.1)/TITLE="'//species// &
      ' vs v3 (int. '//ji_string//') Recalc"/VS/COLOR=4 v3_du, recalc[d=2]'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'PLOT/SET/LINE/TITLE="'//species// &
      ' vs v3 (int. '//ji_string//') Error"/VS v3_du, error[d=3]'
    WRITE(IO_JTC,'(A)') 'PPL XLAB "" ; PPL YLAB "error"'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL PLOT'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'LET abs_error = '//species//'[d=1]-recalc[d=2]'
    WRITE(IO_JTC,'(A)') 'PLOT/SET/LINE/TITLE="'//species// &
      ' vs v3 (int. '//ji_string//') abs error"/VS v3_du, abs_error'
    WRITE(IO_JTC,'(A)') 'PPL XLAB "" ; PPL YLAB "1/cm2"'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL PLOT'
    WRITE(IO_JTC,'(A)')

  END SUBROUTINE process_interval_3_4

  ! --------------------------------------------------------------------------

  SUBROUTINE process_interval_5_8(species)

    USE messy_main_math_lsq, ONLY: poly_fit

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: species

    REAL(DP) :: r2
    INTEGER  :: nheader_v3_du, nheader_xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: v3_du, xxx_eff

    ALLOCATE(v3_du(maxv3), xxx_eff(maxv3))

    ! read input files:
    v3_du_filename = DIR_EFF//'/v3_du.ef'//ji_string
    CALL read_file_1d(v3_du_filename, nheader_v3_du, v3_du)
    xxx_eff_filename = DIR_EFF//'/'//species//'.ef'//ji_string
    CALL read_file_1d(xxx_eff_filename, nheader_xxx_eff, xxx_eff)
    WRITE(IO_LOG,*) 'xxx_eff: ', MINVAL(xxx_eff), '...', MAXVAL(xxx_eff)

    ! calculate parameters:
    a47_xxx(ji-1,:) = 0._dp
    CALL poly_fit(v3_du, xxx_eff, deg_tconst(ji), a47_xxx(ji-1,:), r2, IO_LOG)
    ! Do not use this interval if fit didn't work:
    IF ((ABS(r2)<=TINY(0._dp)).AND.(deg_tconst(ji)/=0)) deg_tconst(ji) = FIT_PROBLEM
    ! Note that for deg_tconst=0, the result is okay even though r2=0
    DEALLOCATE(v3_du, xxx_eff)

    ! write to ferret jnl file:
    WRITE(IO_JTC,'(A)') 'CANCEL DATA/ALL'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_xxx_eff, &
      '/VAR="'//species//'" "../'//TRIM(xxx_eff_filename)//'"'
    WRITE(IO_JTC,'(A,I0,A)') 'FILE/SKIP=', nheader_v3_du, &
      '/VAR="v3_du" "../'//TRIM(v3_du_filename)//'"'
    WRITE(IO_JTC,'(A)') 'GO nextviewport'
    WRITE(IO_JTC,'(A)') 'PLOT/SET/LINE/TITLE="'//species// &
      ' vs v3 (int. '//ji_string//')"/VS/COLOR=2 v3_du, '// &
      species//'[d=1]'
    WRITE(IO_JTC,'(A)') 'PPL XLAB "" ; PPL YLAB "1/cm2"'
    WRITE(IO_JTC,'(A)') 'GO ppl_mylayout ; PPL PLOT'
    WRITE(IO_JTC,'(1P,A,E12.4,A,E12.4,A,E12.4,A,E12.4,A)') 'let fitfunc = (', &
      a47_xxx(ji-1,0), ') + (', a47_xxx(ji-1,1), ')*v3_du + (', &
      a47_xxx(ji-1,2), ')*v3_du^2 + (', a47_xxx(ji-1,3), ')*v3_du^3'
    WRITE(IO_JTC,'(A)') 'PLOT/LINE/OVER/DASH=(0.1,0.1,0.1,0.1)/VS'// &
      '/COLOR=4 v3_du, fitfunc'
    WRITE(IO_JTC,*)

  END SUBROUTINE process_interval_5_8

  ! --------------------------------------------------------------------------

  SUBROUTINE process_tdep_interval_1(species)

    USE messy_main_math_lsq, ONLY: poly_fit

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: species
    !INTEGER, PARAMETER :: nTmin=1, nTmax=NTEMP
    INTEGER, PARAMETER :: nTmin=2, nTmax=12 ! [190,200,...,290]
    INTEGER, PARAMETER :: KT = 2     ! arbitrary value as in temp1.pro
    INTEGER, PARAMETER :: KSTART = 3 ! arbitrary value as in temp1.pro
    REAL(DP) :: r2
    REAL(DP) :: coeff(0:2)
    INTEGER  :: jtemp, nheader_v2, nheader_v3_du, nheader_xxx_eff, k2, k
    REAL(DP), DIMENSION(NTEMP) :: T_norm ! normalized temperatures [K]
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: v2, v3_du
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: v2_tdep, v3_du_tdep, &
      xxx_eff_tdep, xxx_eff_tdep_rel
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: c0, c1, c2

    ALLOCATE(v2(maxv2), v3_du(maxv3), xxx_eff(maxv2,maxv3))
    ALLOCATE(v2_tdep(maxv2,NTEMP), v3_du_tdep(maxv3,NTEMP), &
      xxx_eff_tdep(maxv2,NTEMP), xxx_eff_tdep_rel(maxv2,NTEMP))
    ALLOCATE(c0(maxv2), c1(maxv2), c2(maxv2))

    ! read input files:
    DO jtemp=nTmin,nTmax
      WRITE(suffix,'(A,I3.3,A)') '_', 10*jtemp+170, 'K'
      T_norm(jtemp) = (REAL(10*jtemp+170,dp)-T_ref)/T_ref
      v3_du_filename = DIR_EFF//'/v3_du.ef'//ji_string//suffix
      CALL read_file_1d(v3_du_filename, nheader_v3_du, v3_du)
      v3_du_tdep(:,jtemp) = v3_du(:)
      v2_filename = DIR_EFF//'/v2.ef'//ji_string//suffix
      CALL read_file_1d(v2_filename, nheader_v2, v2)
      v2_tdep(:,jtemp) = v2(:)
      xxx_eff_filename = DIR_EFF//'/'//species//'.ef'//ji_string//suffix
      CALL read_file_2d(xxx_eff_filename, nheader_xxx_eff, xxx_eff)
      xxx_eff_tdep(:,jtemp) = xxx_eff(:,KT)
      WRITE(IO_LOG,*) 'xxx_eff for T-dep: ', MINVAL(xxx_eff), '...', &
        MAXVAL(xxx_eff)
    ENDDO

    ! divide by data at the reference temperature (index 7 = 240 K):
    DO jtemp=nTmin,nTmax
      DO k2 = 1, maxv2
        IF (xxx_eff_tdep(k2,7) /= 0._dp) THEN
          xxx_eff_tdep_rel(k2,jtemp) = xxx_eff_tdep(k2,jtemp)/xxx_eff_tdep(k2,7)
        ELSE
          xxx_eff_tdep_rel(k2,jtemp) = 0._dp
        ENDIF
      ENDDO
    ENDDO

    ! calculate parameters:
    DO k2 = 1, maxv2
      CALL poly_fit(T_norm(nTmin:nTmax), xxx_eff_tdep_rel(k2,nTmin:nTmax), &
        deg_tdep(1), coeff(:), r2, IO_LOG)
      WRITE (IO_LOG,*) 'poly_fit in process_tdep_interval_1:'
      WRITE (IO_LOG,*) 'species     = ', species
      WRITE (IO_LOG,*) 'deg_tdep(1) = ', deg_tdep(1)
      WRITE (IO_LOG,*) 'coeff       = ', coeff(:)
      WRITE (IO_LOG,*) 'r2          = ', r2
      c0(k2) = coeff(0)
      c1(k2) = coeff(1)
      IF (deg_tdep(1)==2) c2(k2) = coeff(2)
    ENDDO

    ! Do not use this interval if fit didn't work:
    IF (ABS(r2)<=TINY(0._dp)) THEN
      deg_tdep(1) = FIT_PROBLEM
      WRITE (IO_LOG,*) 'deg_tdep(1) was changed to: ', deg_tdep(1)
    ENDIF

    DO k = kstart, 57
      c0_1_1_xxx(k) = (c1(4*k-3) - c1(4*k+1)) / (LOG(v2(4*k-3))-LOG(v2(4*k+1)))
      c0_1_2_xxx(k) =  c1(4*k-3) - c0_1_1_xxx(k)*LOG(v2(4*k-3))
      IF (deg_tdep(1)==2) THEN
        c0_2_1_xxx(k) = (c2(4*k-3) - c2(4*k+1)) / (LOG(v2(4*k-3))-LOG(v2(4*k+1)))
        c0_2_2_xxx(k) =  c2(4*k-3) - c0_2_1_xxx(k)*LOG(v2(4*k-3))
      ENDIF
    ENDDO
    k = 58
    c0_1_1_xxx(k) = (c1(4*k-3) - c1(maxv2)) / (LOG(v2(4*k-3))-LOG(v2(maxv2)))
    c0_1_2_xxx(k) =  c1(4*k-3) - c0_1_1_xxx(k)*LOG(v2(4*k-3))
    IF (deg_tdep(1)==2) THEN
      c0_2_1_xxx(k) = (c2(4*k-3) - c2(maxv2)) / (LOG(v2(4*k-3))-LOG(v2(maxv2)))
      c0_2_2_xxx(k) =  c2(4*k-3) - c0_2_1_xxx(k)*LOG(v2(4*k-3))
    ENDIF
    DO k = 1, kstart-1
      c0_1_1_xxx(k) =  c0_1_1_xxx(kstart)
      c0_1_2_xxx(k) =  c0_1_2_xxx(kstart)
      IF (deg_tdep(1)==2) THEN
        c0_2_1_xxx(k) =  c0_2_1_xxx(kstart)
        c0_2_2_xxx(k) =  c0_2_2_xxx(kstart)
      ENDIF
    ENDDO

    c0_1_1_xxx(:) = c0_1_1_xxx(:)
    c0_1_2_xxx(:) = c0_1_2_xxx(:)
    IF (deg_tdep(1)==2) THEN
      c0_2_1_xxx(:) = c0_2_1_xxx(:)
      c0_2_2_xxx(:) = c0_2_2_xxx(:)
    ENDIF

    DEALLOCATE(v2, v3_du, xxx_eff, v2_tdep, v3_du_tdep, xxx_eff_tdep, &
      xxx_eff_tdep_rel, c0, c1, c2)

  END SUBROUTINE process_tdep_interval_1

  ! --------------------------------------------------------------------------

  SUBROUTINE process_tdep_interval_2(species)

    USE messy_main_math_lsq, ONLY: poly_fit

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: species

    INTEGER, PARAMETER :: KF  = 200 ! arbitrary value as in temp2.pro
    INTEGER, PARAMETER :: KV2 =   1 ! arbitrary value as in temp2.pro
    REAL(DP) :: r2
    INTEGER  :: jtemp, nheader_v2, nheader_v3_du, nheader_xxx_eff, k3
    REAL(DP), DIMENSION(NTEMP) :: T_norm ! normalized temperatures [K]
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: v2, v3_du
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: v2_tdep, v3_du_tdep, &
      xxx_eff_tdep, xxx_eff_tdep_rel

    ALLOCATE(v2(maxv2), v3_du(maxv3), xxx_eff(maxv2,maxv3))
    ALLOCATE(v2_tdep(maxv2,NTEMP), v3_du_tdep(maxv3,NTEMP), &
      xxx_eff_tdep(maxv3,NTEMP), xxx_eff_tdep_rel(maxv3,NTEMP))

    ! read input files:
    DO jtemp=1,NTEMP
      WRITE(suffix,'(A,I3.3,A)') '_', 10*jtemp+170, 'K'
      T_norm(jtemp) = (REAL(10*jtemp+170,dp)-T_ref)/T_ref
      v3_du_filename = DIR_EFF//'/v3_du.ef'//ji_string//suffix
      CALL read_file_1d(v3_du_filename, nheader_v3_du, v3_du)
      v3_du_tdep(:,jtemp) = v3_du(:)
      v2_filename = DIR_EFF//'/v2.ef'//ji_string//suffix
      CALL read_file_1d(v2_filename, nheader_v2, v2)
      v2_tdep(:,jtemp) = v2(:)
      xxx_eff_filename = DIR_EFF//'/'//species//'.ef'//ji_string//suffix
      CALL read_file_2d(xxx_eff_filename, nheader_xxx_eff, xxx_eff)
      xxx_eff_tdep(:,jtemp) = xxx_eff(KV2,:)
      WRITE(IO_LOG,*) 'xxx_eff for T-dep: ', MINVAL(xxx_eff), '...', MAXVAL(xxx_eff)
    ENDDO

    ! divide by data at the reference temperature (index 8 = 250 K):
    DO jtemp=1,NTEMP
      DO K3 = 1, maxv3
        IF (xxx_eff_tdep(k3,8) /= 0._dp) THEN
          xxx_eff_tdep_rel(k3,jtemp) = xxx_eff_tdep(k3,jtemp)/xxx_eff_tdep(k3,8)
        ELSE
          xxx_eff_tdep_rel(k3,jtemp) = 0._dp
        ENDIF
      ENDDO
    ENDDO

    ! calculate parameters:
    c1_xxx(:) = 0._dp
    CALL poly_fit(T_norm(:), xxx_eff_tdep_rel(KF,:), &
      deg_tdep(2), c1_xxx(0:deg_tdep(2)), r2, IO_LOG)
    WRITE (IO_LOG,*) 'poly_fit in process_tdep_interval_2:'
    WRITE (IO_LOG,*) 'species     = ', species
    WRITE (IO_LOG,*) 'deg_tdep(2) = ', deg_tdep(2)
    WRITE (IO_LOG,*) 'c1_xxx      = ', c1_xxx(0:deg_tdep(2))
    WRITE (IO_LOG,*) 'r2          = ', r2
    ! Do not use this interval if fit didn't work:
    IF (ABS(r2)<=TINY(0._dp)) THEN
      deg_tdep(2) = FIT_PROBLEM
      WRITE (IO_LOG,*) 'deg_tdep(2) was changed to: ', deg_tdep(2)
    ENDIF
    ! write fitted data to file:
    recalc_filename = DIR_JNL//'/'//species//'_'//ji_string//'.dat'
    OPEN(IO_REC, FILE=recalc_filename, STATUS='UNKNOWN')
    WRITE(IO_REC,'(A)') '# T_norm, xxx_eff(KF), xxx_eff(20), xxx_eff(500)'
    DO jtemp=1,NTEMP
      WRITE(IO_REC,*) T_norm(jtemp), xxx_eff_tdep_rel(KF,jtemp), &
        xxx_eff_tdep_rel(20,jtemp), xxx_eff_tdep_rel(500,jtemp)
    ENDDO
    CLOSE(IO_REC)

    ! write to ferret jnl file:
    WRITE(IO_JTD,'(A)') 'CANCEL DATA/ALL'
    WRITE(IO_JTD,'(A)') 'FILE/SKIP=1/VAR="T_norm, eff59, eff20, eff500" "../'// &
      TRIM(recalc_filename)//'"'
    WRITE(IO_JTD,'(A)') 'GO nextviewport'
    WRITE(IO_JTD,'(1P,A,E12.4,A,E12.4,A,E12.4,A,E12.4,A,E12.4,A)') &
      'let fitfunc = (', c1_xxx(0), ') + (', c1_xxx(1), ')*T_norm + (', &
      c1_xxx(2), ')*T_norm + (', c1_xxx(3), ')*T_norm^3'
    WRITE(IO_JTD,'(A)') 'PLOT/LINE/SET/TITLE="'//species// &
      ' vs T_norm (int. '//ji_string//')"/VS T_norm, fitfunc'
    WRITE(IO_JTD,'(A)') 'PPL XLAB "" ; PPL YLAB "T corr. factor"'
    WRITE(IO_JTD,'(A)') 'GO ppl_mylayout ; PPL PLOT'
    WRITE(IO_JTD,'(A)') 'PLOT/OVER/LINE/TITLE="eff20"/VS'// &
      '/DASH=(0.1,0.1,0.1,0.1)/COLOR=2 T_norm, eff20'
    WRITE(IO_JTD,'(A)') 'PLOT/OVER/LINE/TITLE="eff59"/VS'// &
      '/DASH=(0.1,0.1,0.1,0.1)/COLOR=3 T_norm, eff59'
    WRITE(IO_JTD,'(A)') 'PLOT/OVER/LINE/TITLE="eff500"/VS'// &
      '/DASH=(0.1,0.1,0.1,0.1)/COLOR=4 T_norm, eff500'
    WRITE(IO_JTD,*)

    DEALLOCATE(v3_du, xxx_eff, v2_tdep, v3_du_tdep, xxx_eff_tdep, &
      xxx_eff_tdep_rel)

  END SUBROUTINE process_tdep_interval_2

  ! --------------------------------------------------------------------------

  SUBROUTINE process_tdep_interval_3_8(species)

    USE messy_main_math_lsq, ONLY: poly_fit

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: species

    INTEGER, PARAMETER :: KF=59 ! arbitrary value as in temp3_8.pro
    REAL(DP) :: r2
    REAL(DP), DIMENSION(NTEMP) :: T_norm ! normalized temperatures [K]
    INTEGER  :: jtemp, nheader_v3_du, nheader_xxx_eff, k3
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: v3_du, xxx_eff
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: v3_du_tdep, xxx_eff_tdep, xxx_eff_tdep_rel

    ALLOCATE(v3_du(maxv3), xxx_eff(maxv3))
    ALLOCATE(v3_du_tdep(maxv3,NTEMP), xxx_eff_tdep(maxv3,NTEMP), &
      xxx_eff_tdep_rel(maxv3,NTEMP))

    ! read input files:
    DO jtemp=1,NTEMP
      WRITE(suffix,'(A,I3.3,A)') '_', 10*jtemp+170, 'K'
      T_norm(jtemp) = (REAL(10*jtemp+170,dp)-T_ref)/T_ref
      v3_du_filename = DIR_EFF//'/v3_du.ef'//ji_string//suffix
      CALL read_file_1d(v3_du_filename, nheader_v3_du, v3_du)
      v3_du_tdep(:,jtemp) = v3_du(:)
      xxx_eff_filename = DIR_EFF//'/'//species//'.ef'//ji_string//suffix
      CALL read_file_1d(xxx_eff_filename, nheader_xxx_eff, xxx_eff)
      xxx_eff_tdep(:,jtemp) = xxx_eff(:)
      WRITE(IO_LOG,*) 'xxx_eff for T-dep: ', MINVAL(xxx_eff), '...', &
        MAXVAL(xxx_eff)
    ENDDO

    ! divide by data at the reference temperature (index 8 = 250 K):
    DO jtemp=1,NTEMP
      DO K3 = 1, maxv3
        IF (xxx_eff_tdep(k3,8) /= 0._dp) THEN
          xxx_eff_tdep_rel(k3,jtemp) = xxx_eff_tdep(k3,jtemp)/xxx_eff_tdep(k3,8)
        ELSE
          xxx_eff_tdep_rel(k3,jtemp) = 0._dp
        ENDIF
      ENDDO
    ENDDO

    ! calculate parameters:
    c27_xxx(ji-1,:) = 0._dp
    CALL poly_fit(T_norm(:), xxx_eff_tdep_rel(KF,:), &
      deg_tdep(ji), c27_xxx(ji-1,0:deg_tdep(ji)), r2, IO_LOG)
    WRITE (IO_LOG,*) 'poly_fit in process_tdep_interval_3_8:'
    WRITE (IO_LOG,*) 'species     = ', species
    WRITE (IO_LOG,*) 'deg_tdep(', ji, ') = ', deg_tdep(ji)
    WRITE (IO_LOG,*) 'c27_xxx     = ', c27_xxx(ji-1,0:deg_tdep(ji))
    WRITE (IO_LOG,*) 'r2          = ', r2
    ! Do not use this interval if fit didn't work:
    IF (ABS(r2)<=TINY(0._dp)) THEN
      deg_tdep(ji) = FIT_PROBLEM
      WRITE (IO_LOG,*) 'deg_tdep(', ji, ') was changed to: ', deg_tdep(ji)
    ENDIF
    ! write fitted data to file:
    recalc_filename = DIR_JNL//'/'//species//'_'//ji_string//'.dat'
    OPEN(IO_REC, FILE=recalc_filename, STATUS='UNKNOWN')
    WRITE(IO_REC,'(A)') '# T_norm, xxx_eff(KF), xxx_eff(20), xxx_eff(500)'
    DO jtemp=1,NTEMP
      WRITE(IO_REC,*) T_norm(jtemp), xxx_eff_tdep_rel(KF,jtemp), &
        xxx_eff_tdep_rel(20,jtemp), xxx_eff_tdep_rel(500,jtemp)
    ENDDO
    CLOSE(IO_REC)

    ! write to ferret jnl file:
    WRITE(IO_JTD,'(A)') 'CANCEL DATA/ALL'
    WRITE(IO_JTD,'(A)') 'FILE/SKIP=1/VAR="T_norm, eff59, eff20, eff500" "../'// &
      TRIM(recalc_filename)//'"'
    WRITE(IO_JTD,'(A)') 'GO nextviewport'
    WRITE(IO_JTD,'(1P,A,E12.4,A,E12.4,A,E12.4,A,E12.4,A)') 'let fitfunc = (', &
      c27_xxx(ji-1,0), ') + (', c27_xxx(ji-1,1), ')*T_norm + (', &
      c27_xxx(ji-1,2), ')*T_norm^2 + (', c27_xxx(ji-1,3), ')*T_norm^3'
    WRITE(IO_JTD,'(A)') 'PLOT/LINE/SET/TITLE="'//species// &
      ' vs T_norm (int. '//ji_string//')"/VS T_norm, fitfunc'
    WRITE(IO_JTD,'(A)') 'PPL XLAB "" ; PPL YLAB "T corr. factor"'
    WRITE(IO_JTD,'(A)') 'GO ppl_mylayout ; PPL PLOT'
    WRITE(IO_JTD,'(A)') 'PLOT/OVER/LINE/TITLE="eff20"/VS'// &
      '/DASH=(0.1,0.1,0.1,0.1)/COLOR=2 T_norm, eff20'
    WRITE(IO_JTD,'(A)') 'PLOT/OVER/LINE/TITLE="eff59"/VS'// &
      '/DASH=(0.1,0.1,0.1,0.1)/COLOR=3 T_norm, eff59'
    WRITE(IO_JTD,'(A)') 'PLOT/OVER/LINE/TITLE="eff500"/VS'// &
      '/DASH=(0.1,0.1,0.1,0.1)/COLOR=4 T_norm, eff500'
    WRITE(IO_JTD,*)

    DEALLOCATE(v3_du, xxx_eff, v3_du_tdep, xxx_eff_tdep, xxx_eff_tdep_rel)

  END SUBROUTINE process_tdep_interval_3_8

  ! --------------------------------------------------------------------------

  SUBROUTINE conv_eff_poly

    USE jvpp_mem, ONLY: HLINE1, HLINE2b
    USE messy_cmn_photol_mem ! IP_MAX, ip_*, jname
    IMPLICIT NONE
    LOGICAL :: lex ! file exists?
    INTEGER :: jp

    WRITE(*,'(A)') HLINE1
    WRITE(*,'(A)') '*** STEP 3'
    WRITE(*,'(A)') HLINE1
    WRITE(*,*)

    ! create script that concatenates all jval_cal_*.f90 subroutines:
    OPEN(IO_CAT, FILE=DIR_F90//'/cat_jval.tcsh', STATUS='UNKNOWN')
    WRITE(IO_CAT,'(A)') '#! /bin/tcsh -f'
    WRITE(IO_CAT,'(A)') 'echo "  ! -*- f90 -*-"            > '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo ""                          >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  ! **'//HLINE1//'"        >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo ""                          >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  ! This file contains '// &
      'the subroutines jval_cal and jval_cal_*. They are" >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  ! created '// &
      'automatically by jvpp_step3.f90, do not edit!"     >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo ""                          >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  ! Based on f90 files:"   >> '//INCFILE
    WRITE(IO_CAT,'(A)') '(cd .. ; ls -l *.f90)'// &
      '| sed "s/^/  ! /"                                  >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  ! Input is from:"        >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  !'// &
      ' `( cd ../'//inputdir//'/spectra ; pwd )`"         >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'ls -l ../'//inputdir//'/spectra '// &
      '| sed "s/^/  ! /"                                  >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  !'// &
      ' `( cd ../'//inputdir//'/hardcoded ; pwd )`"       >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'ls -l ../'//inputdir//'/hardcoded '// &
      '| sed "s/^/  ! /"                                  >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  ! Contents of jvpp.log:" >> '//INCFILE
    WRITE(IO_CAT,'(A)') '(cd .. ; cat jvpp.log)'// &
      ' | sed "s/^/  ! /"                                 >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo ""                          >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo "  ! **'//HLINE1//'"        >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'echo ""                          >> '//INCFILE
    WRITE(IO_CAT,'(A)') 'cat jval_cal.f90                 >> '//INCFILE

    ! create ferret jnl file for T-const:
    OPEN(IO_JTC, FILE=DIR_JNL//'/jvpp_step3_tconst.jnl', STATUS='UNKNOWN')
    WRITE(IO_JTC,'(A)') '\CANCEL MODE verify ! like "unset echo" under unix'
    WRITE(IO_JTC,'(A)') 'CANCEL VARIABLE/ALL'
    WRITE(IO_JTC,'(A)') 'CANCEL SYMBOL/ALL'
    WRITE(IO_JTC,'(A)') 'GO initviewport 3 3 noheader 0.6'
    WRITE(IO_JTC,*)

    ! create ferret jnl file for T-dep:
    OPEN(IO_JTD, FILE=DIR_JNL//'/jvpp_step3_tdep.jnl', STATUS='UNKNOWN')
    WRITE(IO_JTD,'(A)') '\CANCEL MODE verify ! like "unset echo" under unix'
    WRITE(IO_JTD,'(A)') 'CANCEL VARIABLE/ALL'
    WRITE(IO_JTD,'(A)') 'CANCEL SYMBOL/ALL'
    WRITE(IO_JTD,'(A)') 'GO initviewport 3 3 noheader 0.6'
    WRITE(IO_JTD,*)

    ! create file for SUBROUTINE jval_cal:
    OPEN(IO_CAL, FILE=DIR_F90//'/jval_cal.f90', STATUS='UNKNOWN')
    WRITE(IO_CAL,'(A)') '  SUBROUTINE jval_cal'
    WRITE(IO_CAL,'(A)')
    WRITE(IO_CAL,'(A)') '    ! CALLs to individual jval_cal_* subroutines:'

    ! create LaTeX file:
    OPEN(IO_TEX, FILE='references_jvpp.tex', STATUS='UNKNOWN')

    WRITE(*,'(A)') '------------------------------------------------------------------------'
    WRITE(*,'(A)') '   species        deg_tconst           deg_tdep   Lyman-alpha or IR'
    WRITE(*,'(A)') '             1 2 3 4 5 6 7 8    1 2 3 4 5 6 7 8   addition'
    WRITE(*,'(A)') '------------------------------------------------------------------------'
    DO jp = 1, IP_MAX
      CALL read_nml(lex, TRIM(jname(jp)))
      IF (lex) THEN
        CALL process_species(TRIM(jname(jp)), jp)
      ELSE
        WRITE(*,'(A10,A)') TRIM(jname(jp)), '   no nml file'
      ENDIF
    ENDDO
    WRITE(*,'(A)') '------------------------------------------------------------------------'

    ! finish ferret jnl file for T-const:
    WRITE(IO_JTC,'(A)') 'GO exitviewport'
    WRITE(IO_JTC,'(A)') 'exit'
    CLOSE(IO_JTC)

    ! finish ferret jnl file for T-dep:
    WRITE(IO_JTD,'(A)') 'GO exitviewport'
    WRITE(IO_JTD,'(A)') 'exit'
    CLOSE(IO_JTD)

    ! finish file for SUBROUTINE jval_cal:
    WRITE(IO_CAL,'(A)')
    WRITE(IO_CAL,'(A)') '  END SUBROUTINE jval_cal'
    WRITE(IO_CAL,'(A)')
    WRITE(IO_CAL,'(A)') '  ! **'//HLINE1
    WRITE(IO_CAL,'(A)')
    CLOSE(IO_CAL)

    CLOSE(IO_CAT)
    CLOSE(IO_TEX)

  END SUBROUTINE conv_eff_poly

  ! --------------------------------------------------------------------------

  SUBROUTINE read_nml(lex, species)

    IMPLICIT NONE

    LOGICAL,          INTENT(OUT) :: lex ! file exists?
    CHARACTER(LEN=*), INTENT(IN)  :: species
    INTEGER :: fstat ! file status
    CHARACTER(LEN=80) :: nml_filename

    NAMELIST /JVPP/ l_hardcoded, lya_ir, deg_tconst, deg_tdep, fj_corr, &
      texrxn, texnote, eqntag

    nml_filename = inputdir//'/spectra/'//species//'.nml'

    ! check if file exists:
    INQUIRE(FILE=nml_filename, EXIST=lex)
    IF (.NOT.lex) RETURN

    ! default values:
    l_hardcoded   = .FALSE.
    lya_ir        = ""
    ! mz_rs_20130404+
    !deg_tconst(:) = (/ 1, 1, 1, 1, 3, 3, 3, 3 /)
    ! changed to 2 to fix problem with J_H2O
    deg_tconst(:) = (/ 2, 1, 1, 1, 3, 3, 3, 3 /)
    ! mz_rs_20130404-
    deg_tdep(:) = 2
    fj_corr     = 7
    texrxn      = ""
    texnote(:)  = ""
    eqntag      = ""

    OPEN(IO_NML, FILE=nml_filename)
    READ(IO_NML, NML=JVPP, IOSTAT=fstat)
    IF (fstat /= 0) THEN
      PRINT *, 'error while reading namelist for ', species
      STOP
    ENDIF
    CLOSE(IO_NML)

  END SUBROUTINE read_nml

  ! --------------------------------------------------------------------------

END MODULE jvpp_step3

! ****************************************************************************
