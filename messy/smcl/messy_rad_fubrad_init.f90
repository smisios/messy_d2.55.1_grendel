! ***************************************************************************
MODULE  messy_rad_fubrad_init


  USE messy_main_constants_mem, ONLY : dp
  USE messy_rad_fubrad_mem
  USE messy_rad_fubrad_srb_kck, ONLY : fubrad_srb_kck_mem_ini   &
                                     , fubrad_srb_kck_mem_clean
  USE messy_rad_fubrad_srb_km,  ONLY : fubrad_srb_km_mem_ini    &
                                     , fubrad_srb_km_mem_clean
  
  IMPLICIT NONE
  PRIVATE
  SAVE
  
  PUBLIC :: fubrad_initialize
  PUBLIC :: rad_fubrad_read_nml_ctrl
  PUBLIC :: fubrad_solar_time_control
  PUBLIC :: fubrad_clean_memory
  PUBLIC :: fubrad_ini_param
  PUBLIC :: fubrad_initialize_fluxes
  PUBLIC :: fubrad_global_flux_ini
  PUBLIC :: fubrad_initialize_cross_sec
  PUBLIC :: fubrad_initialize_cross_sec_dyn
  PUBLIC :: fubrad_initialize_fluxes_dyn
  PUBLIC :: fubrad_ini_param_dyn
  
  ! Interface of public subroutines
  
  INTERFACE fubrad_initialize
     MODULE PROCEDURE fubrad_initialize
  END INTERFACE

  INTERFACE rad_fubrad_read_nml_ctrl
     MODULE PROCEDURE rad_fubrad_read_nml_ctrl
  END INTERFACE

  INTERFACE fubrad_solar_time_control
     MODULE PROCEDURE fubrad_solar_time_control
  END INTERFACE
 
  INTERFACE fubrad_clean_memory
     MODULE PROCEDURE fubrad_clean_memory
  END INTERFACE

  INTERFACE fubrad_ini_param
     MODULE PROCEDURE fubrad_ini_param
  END INTERFACE

  INTERFACE fubrad_initialize_fluxes
     MODULE PROCEDURE fubrad_initialize_fluxes
  END INTERFACE

  INTERFACE fubrad_global_flux_ini
     MODULE PROCEDURE fubrad_global_flux_ini
  END INTERFACE

  INTERFACE fubrad_initialize_cross_sec
     MODULE PROCEDURE fubrad_initialize_cross_sec
  END INTERFACE
  
  INTERFACE fubrad_initialize_cross_sec_dyn
     MODULE PROCEDURE fubrad_initialize_cross_sec_dyn
  END INTERFACE
  
  INTERFACE fubrad_initialize_fluxes_dyn
     MODULE PROCEDURE fubrad_initialize_fluxes_dyn
  END INTERFACE
  
  INTERFACE fubrad_ini_param_dyn
     MODULE PROCEDURE fubrad_ini_param_dyn
  END INTERFACE
  
CONTAINS

  SUBROUTINE  rad_fubrad_read_nml_ctrl(status, iou)
    !
    ! read CTRL namelist, check it, and initialize global switches/variables
    ! Author: Patrick Joeckel, (MPICH), Nov 2007
    !
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close, ucase
    USE messy_main_constants_mem, ONLY: HLINE2
    
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: modstr = 'rad'
    CHARACTER(LEN=*), PARAMETER  :: substr = 'rad_fubrad_read_nml_ctrl'
    LOGICAL                      :: lex          ! file exists?
    INTEGER                      :: fstat        ! file status

    NAMELIST /CTRL_FUBRAD/ solfac, nbands     &
                                 , sr         &
                                 , resolution &
                                 , blev

    ! initialize
    status = 1 ! error

    ! Default resolution of FUBRad:
    nbands = 55

    ! By reducing the truncation of the Chebyshev polynomials of the
    ! Koppers and Murtagh (KM) Schumann-Runge parametrization, computational
    ! efficiency is increased. A truncation of 12, according to KM,
    ! already gives a very good fit. Maximum truncation is 20.
    sr % type = 'STR'
    sr % ntr  = 12

    ! type of resolution in subroutine o3fluxes: 
    !         DEFAULT - standard resolution in o3fluxes
    !         DYNAMIC - for tests with finer or coarser resolution
    resolution = 'DEFAULT'

    ! base level for FUBRAD, i.e. the lowest (w.r.t height) 
    ! pressure level in Pa where FUBRAD is operating
    blev = 7000._dp
    
    CALL read_nml_open(lex, substr, iou, 'CTRL_FUBRAD', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML = CTRL_FUBRAD, IOSTAT = fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_FUBRAD', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    CALL ucase(sr%type)
    
    SELECT CASE(TRIM(sr%type))
    CASE('STR','STROBEL','STR-HR22','STR-HR44','STR-FLX','STROBEL-FLX', &
         'KCK','KOCKARTS','KM','KOPPERS-MURTAGH','ZHU'); CONTINUE
    CASE DEFAULT
      CALL nml_error
      RETURN
    END SELECT
    !
    IF (TRIM(resolution) == 'DEFAULT') THEN ! op_mk_20170831
      SELECT CASE(TRIM(sr%type))
      CASE('STR','STROBEL','ZHU')
        IF (nbands /= 55 .AND. nbands /= 106) THEN
           CALL nml_error
           RETURN
        END IF
      CASE('STR-FLX','STROBEL-FLX','STR-HR44') 
        IF (nbands /= 146) THEN
           CALL nml_error
           RETURN
        END IF
      CASE('STR-HR22') 
        IF (nbands /= 81 .AND. nbands /= 124) THEN
           CALL nml_error
           RETURN
        END IF  
      CASE('KCK','KOCKARTS')
        IF (nbands /= 143) THEN
           CALL nml_error
           RETURN
        END IF
      CASE('KM','KOPPERS-MURTAGH')
        IF (nbands /= 144) THEN
           CALL nml_error
           RETURN
        END IF
      END SELECT
    END IF
    WRITE(*,*) HLINE2
    WRITE(*,*) 'SOLAR CYCLE PARAMETER solfac = ',solfac
    WRITE(*,*) 'NOTE: THIS IS POSSIBLY OBSOLETE DEPENDING ON CPL NAMELIST.'

    WRITE(*,*) ''
    WRITE(*,*) 'FUBRad BAND RESOLUTION :',nbands
    WRITE(*,*) ''
    WRITE(*,*) 'SCHUMANN-RUNGE BANDS:',TRIM(sr%type)
    WRITE(*,*) 'TRUNCATION:',sr%ntr,' (ONLY IMPORTANT FOR SR % TYPE: KM)'

    WRITE(*,*) HLINE2
    ! end check namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! no error
    RETURN
  CONTAINS
    SUBROUTINE nml_error
      WRITE(*,*) TRIM(sr%type),' WRONG TYPE OF SCHUMANN-RUNGE PARAMETRIZATION'
      WRITE(*,*) 'OR '//TRIM(sr%type)//', ',nbands,' COMBINATION.'
      WRITE(*,*) 'POSSIBLE ARE ONLY:'
      WRITE(*,*) 'STR, STR-FLX, STR-HR22, STR-HR44, KCK, KM, ZHU'
      WRITE(*,*) 
      WRITE(*,*) 'STR      -> STROBEL         : NBANDS =  55, 106'
      WRITE(*,*) 'STR-HR22 -> STROBEL HEATR.  : NBANDS =  81'
      WRITE(*,*) 'STR-HR44 -> STROBEL HEATR.  : NBANDS = 146'
      WRITE(*,*) 'STR-FLX  -> STROBEL FLUX    : NBANDS = 146'
      WRITE(*,*) 'KCK      -> KOCKARTS        : NBANDS = 143'
      WRITE(*,*) 'KM       -> KOPPERS-MURTAGH : NBANDS = 144'
      RETURN
    END SUBROUTINE nml_error
  END SUBROUTINE rad_fubrad_read_nml_ctrl
  !
  SUBROUTINE fubrad_solar_time_control(status, cdisse, zsolc_fubrad, zval)
    !
    ! I/O
    INTEGER,  INTENT(OUT)    :: status
    REAL(DP), INTENT(IN)     :: cdisse
    REAL(DP), INTENT(OUT)    :: zsolc_fubrad  ! total solar irradiance (1 AU)
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: zval

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'fubrad_solar_time_control'
    !
    ! SOLAR CONSTANT FOR SOLAR MAXIMUM AND MINIMUM
    !   CAUTION: CORRESPONDS TO MEAN SOLAR CONSTANT OF 1366
    !   BUT VALUE IN messy_main_constants_mem IS 1365
    REAL(dp), PARAMETER    :: solc_max=1366.60_dp, solc_min=1365.65_dp

    status = 1 ! ERROR
    !
    ! SAVE cdisse for internal use in FUBRAD
    cdisse_fubrad = cdisse
    !
    external_sol: IF (PRESENT(zval)) THEN

       IF (.NOT. ASSOCIATED(val)) ALLOCATE(val(SIZE(zval)))
       val(:) = zval(:)
       !
       SELECT CASE(SIZE(val))
       CASE(1)
          !
          lscale = .TRUE.
          !
          ! The F10.7cm (in sfu) used here must be adjusted to the Sun-Earth
          ! distance of 1 AU.
          ! In order to obtain the total solar insolation from the total
          ! solar irradiance (the latter is at 1 AU), the factor cdisse
          ! is multiplied (see messy_rad_e5.f90).
          !
          ! NOTE: expert guess: solfac = 0 for 70 sfu ; = 1 for 270 sfu
          solfac = (val(1) - 70.0) / 200.0
          !
          ! calculate the total solar irradiance (i.e. at 1 AU)
          zsolc_fubrad = solc_max*solfac + solc_min*(1.-solfac)
          !
       !CASE(50, 56, 107)  ! fb_mk_20120119: added resolution 107 and 56
       CASE(56, 82, 107, 125, 144, 145, 147)  ! fb_mk_20150115: 50 disabled, 125 added
          !
          lscale = .FALSE.
          !
          ! The spectral data used here must be adjusted to the Sun-Earth
          ! distance of 1 AU. 
          ! In order to obtain the total solar insolation from the total
          ! solar irradiance (the latter is at 1 AU), the factor cdisse
          ! is multiplied (see messy_rad_e5.f90).
          !
          zsolc_fubrad = val(IDX_TSI) ! total solar irradiance (at 1 AU) [W/m2]
          !
       CASE DEFAULT
          !
          IF (TRIM(resolution) == 'DEFAULT') THEN ! op_mk_20170831
             WRITE(*,*) substr, ': ERROR IN CHOICE OF fubrad_solar (CPL)'//&
                 &'! Number of parameters is ',SIZE(val),&
                 ' but only 1, 56, 82, or 107 are currently implemented.'  ! fb_mk_20131010: 50 disabled
             RETURN
          ELSE
            lscale = .FALSE.
            zsolc_fubrad = val(IDX_TSI) ! total solar irradiance (at 1 AU) [W/m2]
          END IF
       END SELECT

    ELSE
       ! constant solfac (CTRL namelist)

       ! calculate the total solar irradiance (i.e. at 1 AU)
       zsolc_fubrad = solc_max*solfac + solc_min*(1.-solfac)
       !
    END IF external_sol

    status = 0

  END SUBROUTINE fubrad_solar_time_control
  !
  ! ---------------------------------------------------------------------------
  !
  SUBROUTINE fubrad_ini_param(status)
    !
    ! Purpose: Initialize the indices of the Chappuis fluxes in input data array.
    ! --------
    ! Author: Markus Kunze, FU-Berlin, January, 2012.
    ! -------
    INTEGER, INTENT(out) :: status
    !
    ! LOCAL
    INTEGER,               PARAMETER :: nchap_49    = 1
    INTEGER,               PARAMETER :: nchap_55    = 6
    INTEGER,               PARAMETER :: nchap_106   = 57
    INTEGER,               PARAMETER :: nchap_81    = 14
    INTEGER,               PARAMETER :: nherz_def   = 15
    INTEGER,               PARAMETER :: nhart_def   = 10
    INTEGER,               PARAMETER :: nhug_def    = 18

    !
    INTEGER, DIMENSION(2), PARAMETER :: IDX_CHA_49  = (/50,50/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_CHA_55  = (/51,56/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_CHA_106 = (/51,107/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_SRB_81 = (/ 6,24/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HRZ_81 = (/25,39/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HAR_81 = (/40,49/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HUG_81 = (/50,67/)
    INTEGER,               SAVE      :: IDX_ADD_81 = 68
    INTEGER, DIMENSION(2), PARAMETER :: IDX_CHA_81 = (/69,82/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_SRB_124 = (/ 6,24/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HRZ_124 = (/25,39/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HAR_124 = (/40,49/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HUG_124 = (/50,67/)
    INTEGER,               SAVE      :: IDX_ADD_124 = 68
    INTEGER, DIMENSION(2), PARAMETER :: IDX_CHA_124 = (/69,125/)
    !
    INTEGER, DIMENSION(2), PARAMETER :: IDX_SRB_143 = (/28,43/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HRZ_143 = (/44,58/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HAR_143 = (/59,68/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HUG_143 = (/69,86/)
    INTEGER,               SAVE      :: IDX_ADD_143 = 87
    INTEGER, DIMENSION(2), PARAMETER :: IDX_CHA_143 = (/88,144/)
    !
    INTEGER, DIMENSION(2), PARAMETER :: IDX_SRB_144 = (/28,44/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HRZ_144 = (/45,59/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HAR_144 = (/60,69/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HUG_144 = (/70,87/)
    INTEGER,               SAVE      :: IDX_ADD_144 = 88
    INTEGER, DIMENSION(2), PARAMETER :: IDX_CHA_144 = (/89,145/)
    !
    INTEGER, DIMENSION(2), PARAMETER :: IDX_SRB_146 = (/28,46/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HRZ_146 = (/47,61/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HAR_146 = (/62,71/)
    INTEGER, DIMENSION(2), PARAMETER :: IDX_HUG_146 = (/72,89/)
    INTEGER,               SAVE      :: IDX_ADD_146 = 90
    INTEGER, DIMENSION(2), PARAMETER :: IDX_CHA_146 = (/91,147/)
    ! 
    ! ntr_km  -  maximum truncation for Koppers and Murtagh (1996)
    ! ntr_k   -  size of parameter arrays for Kockarts (1994)
    INTEGER,               PARAMETER :: ntr_km = 20
    INTEGER,               PARAMETER :: ntr_k  = 12

    CHARACTER(LEN=*),      PARAMETER :: substr = 'fubrad_ini_param'
    status = 1 ! ERROR
    !
    nherz = nherz_def
    nhart = nhart_def
    nhug  = nhug_def
    SELECT CASE(nbands)
    CASE(49)
       IDX_CHA = IDX_CHA_49
       nchap   = nchap_49
    CASE(55)
       IDX_CHA = IDX_CHA_55
       nchap   = nchap_55
    CASE(81)
       IDX_SRB = IDX_SRB_81
       IDX_HRZ = IDX_HRZ_81
       IDX_HAR = IDX_HAR_81
       IDX_HUG = IDX_HUG_81
       IDX_ADD = IDX_ADD_81
       IDX_CHA = IDX_CHA_81
       nsrb    = nsrb_st
       nchap   = nchap_81
    CASE(106)
       IDX_CHA = IDX_CHA_106
       nchap   = nchap_106
    CASE(124)
       IDX_SRB = IDX_SRB_124
       IDX_HRZ = IDX_HRZ_124
       IDX_HAR = IDX_HAR_124
       IDX_HUG = IDX_HUG_124
       IDX_ADD = IDX_ADD_124
       IDX_CHA = IDX_CHA_124
       nsrb    = nsrb_st
       nchap   = nchap_106
    CASE(143)
       IDX_SRB = IDX_SRB_143
       IDX_HRZ = IDX_HRZ_143
       IDX_HAR = IDX_HAR_143
       IDX_HUG = IDX_HUG_143
       IDX_ADD = IDX_ADD_143
       IDX_CHA = IDX_CHA_143
       nsrb    = nsrb_k
       nchap   = nchap_106   
       ntr_max = ntr_k
    CASE(144)
       IDX_SRB = IDX_SRB_144
       IDX_HRZ = IDX_HRZ_144
       IDX_HAR = IDX_HAR_144
       IDX_HUG = IDX_HUG_144
       IDX_ADD = IDX_ADD_144
       IDX_CHA = IDX_CHA_144
       nsrb    = nsrb_km
       nchap   = nchap_106
       ntr_max = ntr_km
    CASE(146)
       IDX_SRB = IDX_SRB_146
       IDX_HRZ = IDX_HRZ_146
       IDX_HAR = IDX_HAR_146
       IDX_HUG = IDX_HUG_146
       IDX_ADD = IDX_ADD_146
       IDX_CHA = IDX_CHA_146
       nsrb    = nsrb_st
       nchap   = nchap_106  
    CASE DEFAULT
       WRITE(*,*) substr, ': ERROR IN CHOICE OF NBANDS '// &
                          ' ! NBANDS = ',nbands,' NOT VALID; ONLY 49, 55, OR 106.'
       RETURN
    END SELECT
    !
    status = 0
    RETURN
  END SUBROUTINE fubrad_ini_param
  ! 
  ! ---------------------------------------------------------------------------
#include "messy_rad_fubrad_ini_param_dyn.inc"

  SUBROUTINE fubrad_initialize (nbdim,klev)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbdim,klev
    !

    !   allocate module local memory
    
    ALLOCATE (my_pph(nbdim,klev+1))
    ALLOCATE (my_ppf(nbdim,klev))
    ALLOCATE (prmu0(nbdim))
    
    ! SCHEME FIELDS
    ALLOCATE (zFherz(nherz), zFherz_max(nherz), zFherz_min(nherz))
    ALLOCATE (zsigherz2(nherz), zsigherz3(nherz))
    ALLOCATE (zFhart(nhart), zFhart_max(nhart), zFhart_min(nhart))
    ALLOCATE (zsighart(nhart))
    ALLOCATE (zFhug(nhug),   zFhug_max(nhug),   zFhug_min(nhug))
    ALLOCATE (zsighug(nhug))
    ALLOCATE (zFchap(nchap), zFchap_max(nchap), zFchap_min(nchap))
    ALLOCATE (zsigchap(nchap))

    ALLOCATE (zFsrb(nsrb), zFsrb_max(nsrb), zFsrb_min(nsrb))
    ALLOCATE (zFsrc(nsrc), zFsrc_max(nsrc), zFsrc_min(nsrc))
    ALLOCATE (zsigsrc(nsrc))
    ALLOCATE (przsec(nbdim))
    SELECT CASE(TRIM(sr%type))
    CASE('KM', 'KOPPERS-MURTAGH'); CALL fubrad_srb_km_mem_ini(nbdim,klev)
    CASE('kck','KOCKARTS');        CALL fubrad_srb_kck_mem_ini(nbdim,klev)
    END SELECT
    return
  END SUBROUTINE fubrad_initialize
  ! -
  SUBROUTINE fubrad_clean_memory

    IF (ASSOCIATED(val)) DEALLOCATE(val)
    NULLIFY(val)

    DEALLOCATE (my_pph)
    DEALLOCATE (prmu0, przsec)
    
    ! SCHEME FIELDS
    DEALLOCATE (zFherz, zFherz_max, zFherz_min)
    DEALLOCATE (zsigherz2, zsigherz3)
    DEALLOCATE (zFhart, zFhart_max, zFhart_min)
    DEALLOCATE (zsighart)
    DEALLOCATE (zFhug,   zFhug_max,   zFhug_min)
    DEALLOCATE (zsighug)
    DEALLOCATE (zFchap, zFchap_max, zFchap_min)
    DEALLOCATE (zsigchap)
    DEALLOCATE (zFsrb, zFsrb_max, zFsrb_min)
    DEALLOCATE (zFsrc, zFsrc_max, zFsrc_min)
    DEALLOCATE (zsigsrc)
    DEALLOCATE (przsec)
    DEALLOCATE (my_ppf)
    SELECT CASE(TRIM(sr%type))
    CASE('KM', 'KOPPERS-MURTAGH'); CALL fubrad_srb_km_mem_clean
    CASE('KCK','KOCKARTS');        CALL fubrad_srb_kck_mem_clean
    END SELECT

  END SUBROUTINE fubrad_clean_memory

  ! =====================================================================
  !
  SUBROUTINE fubrad_initialize_fluxes(status)
    ! 
    ! Purpose:  Initialize FUBRad local arrays with the predefined
    ! --------  solar fluxes for maximum and minimum conditions.
    !
    ! References:  
    ! ----------- from J. Lean (2007, pers. comm.) integrated over WMO bands
    ! http://www.geo.fu-berlin.de/en/met/ag/strat/research/SOLARIS/Input_data/index.html
    !               J. Lean (Lean et al. 1997, Lean, 2001, Lean et al., 2005)
    !
    INTEGER, INTENT(out) :: status
    !
    status = 0
    !
    ! SOLAR MINIMUM VALUES: September 1986
    ! SOLAR MAXIMUM VALUES: November 1989
    !
    ! fb_mk_20150113+
    !
    ! Flux in Schumann-Runge bands (175-0 - 202.5 nm)
    SELECT CASE(nsrb)
    CASE(nsrb_st)
       zFsrb_max = (/ 8.05984630E-04, 1.18707882E-03, 1.56492603E-03, &
                      1.09544822E-03, 1.55191674E-03, 1.90348575E-03, &
                      2.42908721E-03, 3.56904107E-03, 3.64587082E-03, &
                      3.81670198E-03, 5.13323570E-03, 6.87343083E-03, &
                      8.60657775E-03, 8.98061340E-03, 1.31016332E-02, &
                      1.23855224E-02, 9.50047766E-03, 1.12412954E-02, &
                      1.76614784E-02 /)
       zFsrb_min = (/ 7.20947458E-04, 1.06782775E-03, 1.39397970E-03, &
                      9.75785790E-04, 1.38660860E-03, 1.69673305E-03, &
                      2.14378124E-03, 3.09457770E-03, 3.28200823E-03, &
                      3.47182235E-03, 4.65831139E-03, 6.22634161E-03, &
                      7.82121108E-03, 8.18495666E-03, 1.19256116E-02, &
                      1.13316445E-02, 8.74120835E-03, 1.03579606E-02, &
                      1.63000345E-02 /)
    END SELECT
    
    ! Flux in Herzberg continuum (206.186-243.902 nm, WMO band 18-32)
    ! with the following subinterals:
    !  1. 206.186-208.333 (wavelength range for solar variability: 206.5-208.5nm)
    !  2. 208.333-210.526 (208.5-210.5nm)
    !  3. 210.526-212.766 (210.5-212.5nm)
    !  4. 212.766-215.054 (212.5-215.5nm)
    !  5. 215.054-217.391 (215.5-217.5nm)
    !  6. 217.391-219.780 (217.5-219.5nm)
    !  7. 219.780-222.222 (219.5-222.5nm)
    !  8. 222.222-224.719 (222.5-224.5nm)
    !  9. 224.719-227.273 (224.5-227.5nm)
    ! 10. 227.273-229.885 (227.5-229.5nm)
    ! 11. 229.885-232.558 (229.5-232.5nm)
    ! 12. 232.558-235.294 (232.5-235.5nm)
    ! 13. 235.294-238.095 (235.5-238.5nm)
    ! 14. 238.095-240.964 (238.5-240.5nm) 
    ! 15. 240.964-243.902 (240.5-243.5nm)
    
    !
    zFherz_max = (/ 0.0273643, 0.0445575, 0.0685641, 0.111326,  0.0685956,  & 
                    0.0902246, 0.134914 , 0.125568 , 0.141672,  0.101853,   & 
                    0.155390 , 0.131883,  0.148747,  0.0864101, 0.176693 /)
    zFherz_min = (/ 0.0251796, 0.0423016, 0.0658109, 0.107089,  0.0658054,  & 
                    0.0867229, 0.129816,  0.121350,  0.136316,  0.0982412,  & 
                    0.150068,  0.127245,  0.143429,  0.0831662, 0.171376 /)
    !
    ! Flux in Hartley bands (243.902-277.778 nm, WMO band 33-42) 
    ! with the following subintervals:
    !  1. 243.902-246.914 (243.5-246.5nm)
    !  2. 246.914-250.000 (246.5-249.5nm)
    !  3. 250.000-253.165 (249.5-253.5nm)
    !  4. 253.165-256.410 (253.5-256.5nm)
    !  5. 256.410-259.740 (256.5-259.5nm)
    !  6. 259.740-263.158 (259.5-263.5nm)
    !  7. 263.158-266.667 (263.5-266.5nm)
    !  8. 266.667-270.270 (266.5-270.5nm)
    !  9. 270.270-273.973 (270.5-273.5nm)
    ! 10. 273.973-277.778 (273.5-277.5nm) 
    !
    !
    zFhart_max = (/ 0.170959, 0.155456, 0.200624, 0.212513, 0.365152,  & 
                    0.408812, 0.750635, 0.983705, 0.668155, 0.762717 /)
    zFhart_min = (/ 0.165423, 0.150055, 0.193896, 0.207455, 0.359191, & 
                    0.399182, 0.742960, 0.976196, 0.662157, 0.755110 /)
    !
    ! Flux in Huggins bands (277.778-362.500 nm, WMO band 43-60) 
    ! with the following subintervals:
    !  1. 277.778-281.690 (277.5-281.5nm)
    !  2. 281.690-285.714 (281.5-285.5nm)
    !  3. 285.714-289.855 (285.5-289.5nm)
    !  4. 289.855-294.118 (289.5-294.5nm)
    !  5. 294.118-298.507 (294.5-298.5nm)
    !  6. 298.507-303.030 (298.5-302.5nm)
    !  7. 303.030-307.692 (302.5-307.5nm)
    !  8. 307.692-312.500 (307.5-312.5nm)
    !  9. 312.500-317.500 (312.5-317.5nm)
    ! 10. 317.500-322.500 (317.5-322.5nm)
    ! 11. 322.500-327.500 (322.5-327.5nm)
    ! 12. 327.500-332.500 (327.5-332.5nm)
    ! 13. 332.500-337.500 (332.5-337.5nm)
    ! 14. 337.500-342.500 (337.5-342.5nm)
    ! 15. 342.500-347.500 (342.5-347.5nm)
    ! 16. 347.500-352.500 (347.5-352.5nm)
    ! 17. 352.500-357.500 (352.5-357.5nm)
    ! 18. 357.500-362.500 (357.5-362.5nm)
    !
    zFhug_max = (/ 0.534589, 1.07536, 1.31067, 2.77411, 2.13764, & 
                   1.92990,  3.01134, 3.27898, 3.44617, 3.73954, & 
                   4.27710,  5.06794, 4.66309, 4.84510, 4.69533, & 
                   4.87496,  5.14510, 4.56993  /)
    zFhug_min = (/ 0.518797, 1.06331, 1.29965, 2.76701, 2.12876, & 
                   1.91880,  3.00146, 3.26644, 3.43556, 3.73176, & 
                   4.27059,  5.06376, 4.65588, 4.83874, 4.68749, & 
                   4.86631,  5.13192, 4.55463 /)
    !
    ! Only the Chappuis band changes with spectral resolution.
    !
    SELECT CASE(nchap)
    CASE(1)
       !
       ! CHAPPUIS FOR THE INTERVAL 407.5-682.5nm (AS MORCRETTE INFRARED STARTS AT 680nm) 
       ! Flux Top of Atmosphere (WMO): 492.6 W/m2
       ! Flux Top of Atmosphere (Atlas3): 496.9 W/m2
       ! Best fit to WMO: 322 W/m2
       !
       zFchap_max = (/ 322.193_dp /)
       zFchap_min = (/ 321.875_dp /)
       !
    CASE(6)
       !
       !  Intervals according to WMO(1986) Table 7-8:
       !
       !  425.000 -    475.000    450.000 NM CHAPPUIS band
       !  475.000 -    525.000    500.000 NM CHAPPUIS band
       !  525.000 -    575.000    550.000 NM CHAPPUIS band
       !  575.000 -    625.000    600.000 NM CHAPPUIS band
       !  625.000 -    675.000    650.000 NM CHAPPUIS band
       !  675.000 -    690.000    682.500 NM CHAPPUIS band
       !
       zFchap_max =(/ 9.63409355E+01, 9.70892952E+01, 9.22633696E+01, &
                      8.68381210E+01, 7.87469930E+01, 2.22164418E+01 /)
       !
       zFchap_min =(/ 9.62362517E+01, 9.69963518E+01, 9.21742312E+01, &
                      8.67603030E+01, 7.86763093E+01, 2.21970821E+01 /)
       !
    CASE(14)
       !
       ! 407.500 -  427.679 nm CHAPPUIS band ( 1)
       ! 427.679 -  447.857 nm CHAPPUIS band ( 2)
       ! 447.857 -  468.036 nm CHAPPUIS band ( 3)
       ! 468.036 -  488.214 nm CHAPPUIS band ( 4)
       ! 488.214 -  508.393 nm CHAPPUIS band ( 5)
       ! 508.393 -  528.571 nm CHAPPUIS band ( 6)
       ! 528.571 -  548.750 nm CHAPPUIS band ( 7)
       ! 548.750 -  568.929 nm CHAPPUIS band ( 8)
       ! 568.929 -  589.107 nm CHAPPUIS band ( 9)
       ! 589.107 -  609.286 nm CHAPPUIS band (10)
       ! 609.286 -  629.464 nm CHAPPUIS band (11)
       ! 629.464 -  649.643 nm CHAPPUIS band (12)
       ! 649.643 -  669.821 nm CHAPPUIS band (13)
       ! 669.821 -  690.000 nm CHAPPUIS band (14)
       !
       zFchap_max =(/ 3.49590204E+01, 3.55037154E+01, 4.18772836E+01, &
                      4.10843663E+01, 3.97251964E+01, 3.70621641E+01, &
                      3.78872201E+01, 3.69192701E+01, 3.63705880E+01, &
                      3.52508219E+01, 3.36416883E+01, 3.26547881E+01, &
                      3.08794918E+01, 3.00390112E+01 /)
       zFchap_min =(/ 3.49157151E+01, 3.54579482E+01, 4.18361025E+01, &
                      4.10452625E+01, 3.96844341E+01, 3.70296990E+01, &
                      3.78492023E+01, 3.68835681E+01, 3.63382512E+01, &
                      3.52199916E+01, 3.36104175E+01, 3.26249712E+01, &
                      3.08521285E+01, 3.00128256E+01 /)

    CASE(57)
       !
       !  Intervals according to WMO(1986) Table 7-4
       !  intervals 70 - 125 and (687.5 - 690.0 nm):
       !
       zFchap_max =(/ 8.38023090E+00, 8.84969924E+00, 8.82348939E+00, &
                      8.62602673E+00, 7.58919248E+00, 8.62354606E+00, &
                      9.09461718E+00, 9.74504723E+00, 1.04172103E+01, &
                      1.03503691E+01, 1.04510711E+01, 1.02976645E+01, &
                      1.01567530E+01, 1.04501740E+01, 1.05442585E+01, &
                      9.63371623E+00, 9.86067818E+00, 9.98649227E+00, &
                      9.61737179E+00, 9.85266048E+00, 9.68630639E+00, &
                      9.08395390E+00, 8.93512634E+00, 9.17870593E+00, &
                      9.45238714E+00, 9.42373031E+00, 9.23953472E+00, &
                      9.40423082E+00, 9.39937560E+00, 9.32216943E+00, &
                      8.97250877E+00, 8.98102765E+00, 8.96034554E+00, &
                      9.11598332E+00, 9.02230194E+00, 9.02901041E+00, &
                      8.67777950E+00, 8.86597015E+00, 8.71690018E+00, &
                      8.73061145E+00, 8.52591381E+00, 8.26997562E+00, &
                      8.36369691E+00, 8.25210406E+00, 8.27356254E+00, &
                      8.18377032E+00, 8.09246998E+00, 7.99439162E+00, &
                      7.94865294E+00, 7.38910865E+00, 7.69918205E+00, &
                      7.67550839E+00, 7.61577411E+00, 7.49669016E+00, &
                      7.42217856E+00, 7.37055896E+00, 3.68616883E+00 /)
       !
       zFchap_min =(/ 8.36906226E+00, 8.83813912E+00, 8.81121606E+00, &
                      8.61785533E+00, 7.57606934E+00, 8.61341706E+00, &
                      9.08439706E+00, 9.73335261E+00, 1.04068336E+01, &
                      1.03400894E+01, 1.04407866E+01, 1.02876462E+01, &
                      1.01465919E+01, 1.04408738E+01, 1.05347647E+01, &
                      9.62395199E+00, 9.85050971E+00, 9.97673358E+00, &
                      9.60664312E+00, 9.84281909E+00, 9.67695763E+00, &
                      9.07617550E+00, 8.92836401E+00, 9.17011470E+00, &
                      9.44334927E+00, 9.41398363E+00, 9.22997094E+00, &
                      9.39513668E+00, 9.39024960E+00, 9.31321340E+00, &
                      8.96364451E+00, 8.97248122E+00, 8.95189069E+00, &
                      9.10805935E+00, 9.01434573E+00, 9.02122425E+00, &
                      8.66970407E+00, 8.85832001E+00, 8.70919152E+00, &
                      8.72319644E+00, 8.51819875E+00, 8.26189374E+00, &
                      8.35620500E+00, 8.24432013E+00, 8.26612649E+00, &
                      8.17636815E+00, 8.08513923E+00, 7.98706293E+00, &
                      7.94133635E+00, 7.38280210E+00, 7.69211001E+00, &
                      7.66874609E+00, 7.60913138E+00, 7.49018916E+00, &
                      7.41570301E+00, 7.36408865E+00, 3.68299880E+00 /)
    CASE DEFAULT
       status = 1
    END SELECT
    !
    ! Additional flux in bands (362.5 - 407.5 nm, WMO band 61-69)
    ! (non absorbing spectral regions)
    !
    SELECT CASE(nbands)
    CASE(49)
       zFadd_max = 54.25_dp
       zFadd_min = 54.25_dp
    CASE(55)
       zFadd_max = 8.52516110E+01_dp
       zFadd_min = 8.50803922E+01_dp
    CASE(81,106,124,143,144,146)
       zFadd_max = 5.48537541E+01_dp
       zFadd_min = 5.47216557E+01_dp
    CASE DEFAULT
       status = 1
    END SELECT
    RETURN
  END SUBROUTINE fubrad_initialize_fluxes

  SUBROUTINE fubrad_initialize_fluxes_dyn(status)
    ! 
    ! Purpose:  Initialize FUBRad local arrays with the predefined
    ! --------  solar fluxes for maximum and minimum conditions.
    !
    INTEGER, INTENT(out) :: status
    !
    status = 0
    !
    ! SOLAR MINIMUM VALUES: September 1986
    ! SOLAR MAXIMUM VALUES: November 1989
    !
    ! Flux in Schumann-Runge bands (175-0 - 202.5 nm)
    SELECT CASE(nsrb)
    CASE(nsrb_st)
       zFsrb_max = (/ 8.05984630E-04, 1.18707882E-03, 7.86775613E-04, &
                      1.09544822E-03, 1.55191674E-03, 1.81728849E-03, &
                      2.42908721E-03, 3.49670241E-03, 3.71820949E-03, &
                      3.81670198E-03, 5.13323570E-03, 6.87343083E-03, &
                      8.60657775E-03, 9.79054973E-03, 1.31016332E-02, &
                      1.23855224E-02, 9.50047766E-03, 1.53692742E-02, &
                      2.24848271E-02 /)
       zFsrb_min = (/ 7.20947458E-04, 1.06782775E-03, 7.06996980E-04, &
                      9.75785790E-04, 1.38660860E-03, 1.61996065E-03, &
                      2.14378124E-03, 3.02947248E-03, 3.34711345E-03, &
                      3.47182235E-03, 4.65831139E-03, 6.22634161E-03, &
                      7.82121108E-03, 8.92851837E-03, 1.19256116E-02, &
                      1.13316445E-02, 8.74120835E-03, 1.41599244E-02, &
                      2.07382260E-02 /)
    END SELECT
    ! 
    ! NOTE:
    ! -----
    ! The values in the following bands are only dummy values and not
    ! intended to by used by FUBRAD. In this case the model needs input
    ! from the IMPORT channel.
    !
    ! Flux in Herzberg continuum (206.186-243.902 nm, WMO band 18-32)
    !
    zFherz_max = 0._dp
    zFherz_min = 0._dp
    !
    ! Flux in Hartley bands (243.902-277.778 nm, WMO band 33-42) 
    ! 
    zFhart_max = 0._dp
    zFhart_min = 0._dp
    !
    ! Flux in Huggins bands (277.778-362.500 nm, WMO band 43-60) 
    !
    zFhug_max = 0._dp
    zFhug_min = 0._dp
    !
    ! Chappuis bands.
    !
    zFchap_max = 0._dp
    zFchap_min = 0._dp
    !
    ! Additional flux in bands (362.5 - 407.5 nm, WMO band 61-69)
    ! (non absorbing spectral regions)
    !
    zFadd_max = 0._dp
    zFadd_min = 0._dp
   
    RETURN
  END SUBROUTINE fubrad_initialize_fluxes_dyn
  ! 
  ! =====================================================================
  !
  SUBROUTINE fubrad_initialize_cross_sec(status)
    ! 
    ! Purpose:  Initialize FUBRad local arrays with the absorbtion
    ! --------  cross sections for O3 and O2
    !
    ! References: WMO report (1986)
    ! ----------- Molina and Molina (1986)
    !
    INTEGER, INTENT(out) :: status
    !
    status = 0
    !
    ! cross sections for Schumann-Runge continuum taken from
    ! Strobel (1978) table 1, central wavelength 126.0 - 174.0 nm
    !
    zsigsrc = (/ 4.3E-19_dp, 2.8E-19_dp, 5.0E-19_dp, 1.4E-18_dp, 2.3E-18_dp, &
                 8.0E-18_dp, 1.3E-17_dp, 1.4E-17_dp, 1.5E-17_dp, 1.5E-17_dp, &
                 1.3E-17_dp, 1.2E-17_dp, 1.1E-17_dp, 1.0E-17_dp, 8.5E-18_dp, &
                 7.3E-18_dp, 6.0E-18_dp, 4.7E-18_dp, 3.4E-18_dp, 2.5E-18_dp, &
                 1.8E-18_dp, 1.2E-18_dp, 8.5E-19_dp, 5.9E-19_dp, 3.7E-19_dp /)
    !
    ! cross sections from WMO report 1986
    !
    zsigherz2 = (/ 7.33E-24_dp, 6.99E-24_dp, 6.45E-24_dp, 5.81E-24_dp, &
                   5.23E-24_dp, 4.71E-24_dp, 4.26E-24_dp, 3.80E-24_dp, &
                   3.35E-24_dp, 2.90E-24_dp, 2.45E-24_dp, 2.05E-24_dp, &
                   1.69E-24_dp, 1.30E-24_dp, 0.93E-24_dp /)
    !
    ! cross sections from Molina&Molina 1986 (263K where available 298K elsewhere)    
    !                                                                                  
    zsigherz3 = (/ 4.325E-19_dp, 53.87E-20_dp, 69.32E-20_dp, 90.25E-20_dp, &
                   118.0E-20_dp, 153.7E-20_dp, 198.8E-20_dp, 254.7E-20_dp, &
                   322.1E-20_dp, 401.0E-20_dp, 490.1E-20_dp, 589.7E-20_dp, &
                   696.5E-20_dp, 806.7E-20_dp, 914.6E-20_dp /)
    !
    ! cross sections from Molina&Molina 1986 (263K where available 298K elsewhere)    
    ! 
    zsighart = (/ 1007E-20_dp,  1088E-20_dp, 1132E-20_dp,  1155E-20_dp,  &
                  1129E-20_dp,  1071E-20_dp, 973.5E-20_dp, 844.2E-20_dp, &
                  699.8E-20_dp, 548.4E-20_dp /)
    !
    ! cross sections from Molina&Molina 1986 (263K where available 298K elsewhere)    
    ! 
    zsighug = (/ 405.5E-20_dp,  280.1E-20_dp,  181.2E-20_dp,  109.8E-20_dp,  &
                 62.69E-20_dp,  34.31E-20_dp,  18.38E-20_dp,  9.659E-20_dp,  &
                 4.916E-20_dp,  2.457E-20_dp,  1.175E-20_dp,  0.5988E-20_dp, &
                 0.2627E-20_dp, 0.1116E-20_dp, 0.0586E-20_dp, 2.66E-22_dp,   &
                 1.09E-22_dp,   5.49E-23_dp /)
    !
    ! Only the Chappuis band changes with spectral resolution.
    !
    SELECT CASE(nchap)
    CASE(1)
       !
       ! CHAPPUIS FOR THE INTERVAL 407.5-682.5nm (AS MORCRETTE INFRARED STARTS AT 680nm) 
       ! Best fit to WMO: 322 W/m2
       !
       zsigchap = (/ 3.157E-21_dp /)
       !
    CASE(6)
       !
       ! cross sections from WMO report 1986 Table 7-8 Band (450(+/-25) - 650.0(+/-25) nm)
       !
       ! K. Bogumil, J. Orphal, and J. P. Burrows,
       ! University of Bremen - Institute of Environmental Physics
       ! SCIAMACHY PFM Satellite Spectrometer:
       !    value taken at 682.5954 nm 243K for Band (675.0 - 690.0 nm)
       !
       zsigchap = (/ 2.34E-22, 1.25E-21, 3.39E-21, &
                     4.46E-21, 2.47E-21, 1.32E-21 /)
       !
    ! op_mk_20170831+
    CASE(14)
       !
       ! cross sections Brion et al. (1998):
       !   O3_CRS_BDM_243K.dat (407.50-519.01 nm); 
       !   O3_CRS_BDM_295K.dat (519.02-690.00 nm)
       !
       zsigchap = (/  3.8600E-23, 1.1161E-22, 2.7951E-22, 6.1046E-22, &
                      1.1516E-21, 1.8443E-21, 2.9661E-21, 3.9149E-21, &   
                      4.6335E-21, 4.9137E-21, 4.1407E-21, 3.0127E-21, &
                      2.1129E-21, 1.3912E-21 /)

    CASE(57)
       
       !
       ! cross sections from WMO report 1986 Band (70-125) (407.5 - 687.5 nm)
       !
       ! K. Bogumil, J. Orphal, and J. P. Burrows,
       ! University of Bremen - Institute of Environmental Physics
       ! SCIAMACHY PFM Satellite Spectrometer:
       !
       ! value taken at 688.6679 nm 243K for Band (687.5 - 690.0 nm)
       !
       zsigchap = (/  2.91E-23_dp, 3.14E-23_dp, 3.99E-23_dp, 6.54E-23_dp, &
                      6.83E-23_dp, 8.66E-23_dp, 1.25E-22_dp, 1.49E-22_dp, &
                      1.71E-22_dp, 2.12E-22_dp, 3.57E-22_dp, 3.68E-22_dp, &
                      4.06E-22_dp, 4.89E-22_dp, 7.11E-22_dp, 8.43E-22_dp, &
                      8.28E-22_dp, 9.09E-22_dp, 1.22E-21_dp, 1.62E-21_dp, &
                      1.58E-21_dp, 1.60E-21_dp, 1.78E-21_dp, 2.07E-21_dp, &
                      2.55E-21_dp, 2.74E-21_dp, 2.88E-21_dp, 3.07E-21_dp, &
                      3.17E-21_dp, 3.36E-21_dp, 3.88E-21_dp, 4.31E-21_dp, &
                      4.67E-21_dp, 4.75E-21_dp, 4.55E-21_dp, 4.35E-21_dp, &
                      4.42E-21_dp, 4.61E-21_dp, 4.89E-21_dp, 4.84E-21_dp, &
                      4.54E-21_dp, 4.24E-21_dp, 3.90E-21_dp, 3.60E-21_dp, &
                      3.43E-21_dp, 3.17E-21_dp, 2.74E-21_dp, 2.61E-21_dp, &
                      2.42E-21_dp, 2.20E-21_dp, 2.02E-21_dp, 1.85E-21_dp, &
                      1.67E-21_dp, 1.54E-21_dp, 1.42E-21_dp, 1.25E-21_dp, &
                      1.70E-21_dp /)
    CASE DEFAULT
       status = 1
    END SELECT
    RETURN
  END SUBROUTINE fubrad_initialize_cross_sec
 
#include "messy_rad_fubrad_initialize_cross_sec_dyn.inc"

  !
  ! =====================================================================
  !
  SUBROUTINE fubrad_global_flux_ini
    ! 
    ! Purpose:  Initialize fluxes for FUBRad for the next time step.
    ! --------  Called from rad_fubrad_global_start.
    !
    ! The code originally distributed in several subroutines is 
    ! gathered here to initialize the fluxes only once.
    !
    ! Author: Markus Kunze, FU-Berlin, January, 2012.
    ! -------
    ! zFlya    - Lyman-Alpha             (121.5 nm)  in W m^-2
    ! zFsch1   - Schuman-Runge Continuum (125-152nm) in mW m^-2
    ! zFs      - Schuman-Runge Continuum (152-166nm) in mW m^-2
    ! zFl      - Schuman-Runge Continuum (166-175nm) in mW m^-2
    ! zFd      - difference, zFs - zFl
    ! soscale  - Schumann-Runge Bands    (175-205nm) dimensionless
    !            solar scaling factormean (=1), max (>1), min (<1)
    ! 
    REAL(dp), PARAMETER :: zFlya_max = 0.00853464_dp
    REAL(dp), PARAMETER :: zFlya_min = 0.00548156_dp
    REAL(dp), PARAMETER :: zFsch1_max = 1.93832_dp
    REAL(dp), PARAMETER :: zFsch1_min = 1.55666_dp
    REAL(dp), PARAMETER :: zFs_max = 1.86478_dp
    REAL(dp), PARAMETER :: zFs_min = 1.61346_dp
    REAL(dp), PARAMETER :: zFl_max = 3.90776_dp 
    REAL(dp), PARAMETER :: zFl_min = 3.48729_dp
    REAL(dp), PARAMETER :: soscale_max = 1.05129_dp
    REAL(dp), PARAMETER :: soscale_min = 0.960729_dp
    
    !* Calculate Flux at top of atmosphere from max and min values 
    !  and state of solar cycle:
    !
    IF (lscale) THEN
       zFsrb(:)  = zFsrb_min(:) *(1.-solfac) + zFsrb_max(:) *solfac
       zFherz(:) = zFherz_min(:)*(1.-solfac) + zFherz_max(:)*solfac
       zFhart(:) = zFhart_min(:)*(1.-solfac) + zFhart_max(:)*solfac
       zFhug(:)  = zFhug_min(:) *(1.-solfac) + zFhug_max(:) *solfac
       zFchap(:) = zFchap_min(:)*(1.-solfac) + zFchap_max(:)*solfac
       fladd     = zFadd_min    *(1.-solfac) + zFadd_max    *solfac
       !
       ! Lyman-Alpha             (126.5 nm)  in W m^-2
       zFlya = zFlya_max*solfac + zFlya_min*(1.-solfac)
       !
       ! Schuman-Runge Continuum (125-175nm) in mW m^-2:
       !
       zFsch1 = zFsch1_max * solfac + zFsch1_min * (1.-solfac)
       !
       zFs = zFs_max*solfac + zFs_min*(1.-solfac)
       zFl = zFl_max*solfac + zFl_min*(1.-solfac)
       !
       ! scaling factor for Schumann-Runge bands:
       !
       soscale = soscale_max * solfac + soscale_min * (1.-solfac)
    ELSE
       zFlya   = val(IDX_LYA)
       zFsch1  = val(IDX_SR1)
       zFs     = val(IDX_SRS)
       zFl     = val(IDX_SRL)
       soscale = val(IDX_SRB(1))
       zFsrc(:) = val(IDX_SRC(1):IDX_SRC(2))
       zFsrb(:) = val(IDX_SRB(1):IDX_SRB(2))
       zFherz(:)= val(IDX_HRZ(1):IDX_HRZ(2))
       zFhart(:)= val(IDX_HAR(1):IDX_HAR(2))
       zFhug(:) = val(IDX_HUG(1):IDX_HUG(2))
       zFchap(:)= val(IDX_CHA(1):IDX_CHA(2))
       fladd    = (zFadd_min + zFadd_max)*0.5_dp
       IF (nbands > 49 &
           .OR. TRIM(resolution)=='DYNAMIC') &
           fladd = val(IDX_ADD)
    END IF
    !
    zFd = zFs - zFl  ! d: difference, zFs-zFl (Schuman-Runge Continuum)
    !
    RETURN
  END SUBROUTINE fubrad_global_flux_ini
  !
END MODULE  messy_rad_fubrad_init
