! ****************************************************************************
!                Time-stamp: <2020-09-15 16:37:38 sander>
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

MODULE jvpp_step1

  USE messy_main_constants_mem, ONLY: DP
  USE jvpp_mem, ONLY: HLINE2b, DIR_176, &
                      l_sig, l_s2t, l_s3t, l_phi, l_tc1, l_tc2, l_tc3, &
                      IO_SIG, IO_176, NWAV, wave_nm, inputdir
  IMPLICIT NONE
  REAL(DP) :: wave_nm_min(NWAV), wave_nm_max(NWAV)

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE process_species(species, filetype, l_ex, n, maxcso, mini, maxi)

    USE messy_main_math_spline
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: species
    CHARACTER(LEN=*), INTENT(IN)  :: filetype
    LOGICAL,          INTENT(OUT) :: l_ex ! file exists ?
    INTEGER,          INTENT(OUT) :: n
    REAL(DP),         INTENT(OUT) :: maxcso, mini, maxi

    INTEGER, PARAMETER :: NMAXWAV = 5000
    INTEGER :: jw, iostatus, nheader, N_tdep, jtdep
    REAL(DP), DIMENSION(NMAXWAV) :: waveo
    REAL(DP), DIMENSION(NMAXWAV,3) :: cso
    REAL(DP), DIMENSION(3) :: csn
    CHARACTER(LEN=80) :: header, infile, outfile, spline_method
    REAL(DP) :: dummy

    infile  = inputdir//'/spectra/'//species//filetype
    outfile = DIR_176//'/'//species//filetype//'_176'
    ! check if file exists:
    INQUIRE(file=TRIM(infile), exist=l_ex)
    IF (.NOT.l_ex) THEN
      RETURN
    ENDIF

    N_tdep = 1
    IF (filetype=='.s2t') N_tdep = 2
    IF (filetype=='.s3t') N_tdep = 3
    IF (filetype=='.tc2') N_tdep = 2
    IF (filetype=='.tc3') N_tdep = 3

    OPEN(IO_SIG, FILE=TRIM(infile),  STATUS='OLD')
    OPEN(IO_176, FILE=TRIM(outfile), STATUS='UNKNOWN')

    nheader = 0
    DO
      READ(IO_SIG,'(A)') header
      IF (header(1:1)/='#') EXIT ! exit do loop if not a header line
      nheader = nheader + 1 ! count header lines
    ENDDO
    BACKSPACE(IO_SIG) ! back to previous line

    IF ((filetype=='.s2t').OR.(filetype=='.s3t')) THEN
      ! read 2 lines with temperature infos and copy them to outputfile:
      READ(IO_SIG,'(A)') header
      WRITE(IO_176,'(A)') TRIM(header)
      READ(IO_SIG,'(A)') header
      WRITE(IO_176,'(A)') TRIM(header)
    ENDIF

    n = 1
    DO
      READ(IO_SIG,*,IOSTAT=iostatus) waveo(n), cso(n,1:N_tdep)
      IF (iostatus < 0) EXIT ! exit do loop at end of file
      ! check for monotonicity of input wavelength:
      IF (n>1) THEN
        IF (waveo(n)<=waveo(n-1)) THEN
          PRINT *, species, waveo(n-1), waveo(n)
          STOP "ERROR: wavelength in input file does not increase."
        ENDIF
      ENDIF
      n = n + 1
      IF (n > NMAXWAV) THEN
         print *, 'ERROR: NMAXWAV limit reached (',NMAXWAV,')!'
         STOP
      ENDIF
    ENDDO
    n = n - 1
    CLOSE(IO_SIG)

    DO jw = 1,NWAV
      DO jtdep=1, N_tdep
        IF ((wave_nm(jw)>=waveo(1)).AND.(wave_nm(jw)<=waveo(n))) THEN
          ! select default interpolation method:
          spline_method = "integration"
          !spline_method = "linear_val"
          ! select another method, if necessary:
          ! IF (species=="XYZ") spline_method = "constant_val"
          SELECT CASE(TRIM(spline_method))
          CASE ("integration")   ; CALL spline_integration   &
            (waveo(1:n), cso(1:n,jtdep), wave_nm_min(jw), wave_nm_max(jw), csn(jtdep))
          CASE ("constant_val")  ; CALL spline_constant_val  &
            (n, waveo(1:n), cso(1:n,jtdep), wave_nm(jw), csn(jtdep))
          CASE ("linear_val")    ; CALL spline_linear_val    &
            (n, waveo(1:n), cso(1:n,jtdep), wave_nm(jw), csn(jtdep), dummy)
          CASE ("quadratic_val") ; CALL spline_quadratic_val &
            (n, waveo(1:n), cso(1:n,jtdep), wave_nm(jw), csn(jtdep), dummy) ! n must be odd
          CASE ("b_val")         ; CALL spline_b_val         &
            (n, waveo(1:n), cso(1:n,jtdep), wave_nm(jw), csn(jtdep))
          CASE DEFAULT ; STOP "ERROR: unknown spline method"
          END SELECT
        ELSE
          ! outside of the measured wavelength range, assume zero:
          csn(jtdep) = 0._dp
        ENDIF
      ENDDO
      WRITE(IO_176,'(3(1PE13.5))') csn(1:N_tdep)
    ENDDO
    CLOSE(IO_176)

    maxcso = MAXVAL(cso(1:n,N_tdep))
    mini   = MINVAL(waveo(1:n))
    maxi   = MAXVAL(waveo(1:n))

    ! comparison plots cannot be made for all species:
    IF (species=="O3P")     RETURN ! dat_m17 uses s3t file
    IF (species=="O1D")     RETURN ! dat_m17 uses s3t file
    IF (species=="NO2O")    RETURN ! dat_m17 uses s2t file
    IF (species=="NOO2")    RETURN ! dat_m17 uses s2t file
    IF (species=="COH2")    RETURN ! dat_m17 uses s2t file
    IF (species=="CHOH")    RETURN ! dat_m17 uses s2t file
    IF (species=="ClNO3")   RETURN ! dat_m17 uses s3t file
    IF (species=="CH3Cl")   RETURN ! dat_m17 uses s3t file
    IF (species=="CH3CCl3") RETURN ! dat_m17 uses s3t file
    IF (species=="ClONO2")  RETURN ! not available in dat_m17

  END SUBROUTINE process_species

  ! --------------------------------------------------------------------------

  SUBROUTINE spline_integration(xdata, ydata, xmin, xmax, yval)

    IMPLICIT NONE
    REAL(DP), INTENT(IN)  :: xdata(:) ! waveo(1:n)
    REAL(DP), INTENT(IN)  :: ydata(:) ! cso(1:n,jtdep)
    REAL(DP), INTENT(IN)  :: xmin     ! wave_nm_min(jw)
    REAL(DP), INTENT(IN)  :: xmax     ! wave_nm_max(jw)
    REAL(DP), INTENT(OUT) :: yval     ! csn(jtdep)
    INTEGER :: i,n
    REAL(DP) :: xstart, xend

    ! Integrate the spectrum from xmin to xmax. The spectrum is
    ! given by ydata(:) at the wavelengths xdata(:). Between the data
    ! points, linear interpolation is used.
    n = SIZE(xdata)
    yval = 0._dp
    DO i = 1,n-1
      xstart = MAX(xdata(i),xmin)
      xend   = MIN(xdata(i+1),xmax)
      IF (xend>xstart) THEN
        yval = yval + (xend-xstart) * &
          (ydata(i)+(ydata(i+1)-ydata(i))* &
          (((xend+xstart)/2._dp-xdata(i))/(xdata(i+1)-xdata(i))))
      ENDIF
    ENDDO
    yval = yval / (xmax-xmin)

  END SUBROUTINE spline_integration

  ! --------------------------------------------------------------------------

  SUBROUTINE conv_sig_176

    USE messy_cmn_photol_mem ! IP_MAX, ip_*, jname

    USE jvpp_mem, ONLY: HLINE1, IO_LOG
    IMPLICIT NONE
    INTEGER  :: jp, n, jw
    REAL(DP) :: maxcso, mini, maxi

    WRITE(*,'(A)') HLINE1
    WRITE(*,'(A)') '*** STEP 1'
    WRITE(*,'(A)') HLINE1
    WRITE(*,*)

    ! write wave.dat (only for ferret plots):
    OPEN(IO_176, FILE=DIR_176//'/wave.dat', STATUS='UNKNOWN')
    DO jw = 1,NWAV
      WRITE(IO_176,'(F8.2)') wave_nm(jw)
      ! activate next line to print sigma(Cl2O2) for lambda> 420nm, see
      ! JPL2011 for details: 
      ! PRINT *, wave_nm(jw), 9.5E-16 * exp(-0.0281*wave_nm(jw))
    ENDDO
    CLOSE(IO_176)

    WRITE(IO_LOG,*) HLINE1
    ! calculate lower and upper borders of wavelength range:
    DO jw=1,NWAV
      IF (jw==1) THEN
        wave_nm_min(jw) = wave_nm(1) - (wave_nm(2)-wave_nm(1)) / 2._dp
      ELSE
        wave_nm_min(jw) = (wave_nm(jw) + wave_nm(jw-1)) / 2._dp
      ENDIF
      IF (jw==NWAV) THEN
        wave_nm_max(jw) = wave_nm(jw) + (wave_nm(jw)-wave_nm(jw-1)) / 2._dp
      ELSE
        wave_nm_max(jw) = (wave_nm(jw) + wave_nm(jw+1)) / 2._dp
      ENDIF
        WRITE(IO_LOG,'(I5,3F10.3)') &
          jw, wave_nm_min(jw), wave_nm(jw), wave_nm_max(jw)
    ENDDO

    ! loop over all photolysis reactions:
    DO jp=1,IP_MAX
      ! --------------------------------------------
      ! spectra (sigma):
      CALL process_species(TRIM(jname(jp)), '.sig', l_sig(jp), n, &
        maxcso, mini, maxi)
      CALL process_species(TRIM(jname(jp)), '.s2t', l_s2t(jp), n, &
        maxcso, mini, maxi)
      CALL process_species(TRIM(jname(jp)), '.s3t', l_s3t(jp), n, &
        maxcso, mini, maxi)
      IF (l_sig(jp)) WRITE(*,'(A)', ADVANCE='NO') 'sigma '
      IF (l_s2t(jp)) WRITE(*,'(A)', ADVANCE='NO') 's2t   '
      IF (l_s3t(jp)) WRITE(*,'(A)', ADVANCE='NO') 's3t   '
      IF (l_sig(jp).OR.l_s2t(jp).OR.l_s3t(jp)) THEN
        WRITE(*,'(A,1PE8.2,A,I5,A,0PF6.2,A,0PF6.2,2A)') 'up to ', &
          maxcso, ' in', n, ' data points from ', &
          mini, '-', maxi, ' nm for '//TRIM(jname(jp))
      ELSE
        WRITE(*,'(A)') HLINE2b//' '//TRIM(jname(jp))
      ENDIF
      ! --------------------------------------------
      ! quantum yields (phi):
      CALL process_species(TRIM(jname(jp)), '.phi', l_phi(jp), n, &
        maxcso, mini, maxi)
      IF (l_phi(jp)) THEN
        WRITE(*,'(A,I4,A,0PF6.2,A,0PF6.2,A)') '   - quantum yields '// &
          'phi:', n, ' data points from ', mini, '-', maxi, ' nm'
      ENDIF
      ! --------------------------------------------
      ! temperature coefficients:
      CALL process_species(TRIM(jname(jp)), '.tc1', l_tc1(jp), n, &
        maxcso, mini, maxi)
      CALL process_species(TRIM(jname(jp)), '.tc2', l_tc2(jp), n, &
        maxcso, mini, maxi)
      CALL process_species(TRIM(jname(jp)), '.tc3', l_tc3(jp), n, &
        maxcso, mini, maxi)
      IF (l_tc1(jp).OR.l_tc2(jp).OR.l_tc3(jp)) THEN
        WRITE(*,'(A,I4,A)', ADVANCE='NO') '   - T-dep with '
        IF (l_tc1(jp)) WRITE(*,'(A)', ADVANCE='NO') '1'
        IF (l_tc2(jp)) WRITE(*,'(A)', ADVANCE='NO') '2'
        IF (l_tc3(jp)) WRITE(*,'(A)', ADVANCE='NO') '3'
        WRITE(*,'(A,I4,A,0PF6.2,A,0PF6.2,A)') ' coeff:', &
          n, ' data points from ', mini, '-', maxi, ' nm'
      ENDIF
      ! --------------------------------------------
    ENDDO

    PRINT *

  END SUBROUTINE conv_sig_176

  ! --------------------------------------------------------------------------

END MODULE jvpp_step1

! ****************************************************************************
