! ***********************************************************************
MODULE messy_spacenox
  ! ***********************************************************************
  ! Parameterization of NOx sources related to solar, interplanetary, and
  ! geomagnetic activity 
  !
  ! Currently implemented:
  ! 1. Energetic Particle Precipitation Indirect Effect (EPP IE), also called
  !    thermospheric NOx. Low-energy electrons, related to auroral activity,
  !    produce NO in the thermosphere at about 110 km. Downward transport 
  !    (diffusion, advection) can transport the NO down into the mesosphere.
  ! 
  ! 2. Galactic Cosmic Rays (GCR). Production of N, NO and OH according
  !    to Brasseur and Solomon, Aeronomy of the Middle Atmosphere, p329-330
  ! 
  ! Authors: Andreas Baumgaertner,
  ! MPI-CH Mainz (abaumg@mpch-mainz.mpg.de), Sept 2007

  ! TODO:
  !    see TODO list in messy_spacenox_e5.f90

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'spacenox'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.4'

  ! GLOBAL CTRL NAMELIST VARIABLES
  REAL(dp),     PUBLIC :: EPPIE_latN =  55.0_dp
  REAL(dp),     PUBLIC :: EPPIE_latS = -55.0_dp
  ! coefficient for scaling of EPPIE in NH/SH (read from namelist)
  REAL(dp),     PUBLIC :: EPPIE_coeffN = 1.0_dp 
  REAL(dp),     PUBLIC :: EPPIE_coeffS = 1.0_dp ! 
  ! coefficient for scaling of EPPIE in NH/SH according to time of year
  REAL(dp),     PUBLIC :: EPPIE_seasoncoeffN = 1.0_dp 
  REAL(dp),     PUBLIC :: EPPIE_seasoncoeffS = 1.0_dp 

  ! Ap DATA LINEARLY INTERPOLATED IN TIME
  REAL(DP), DIMENSION(:), POINTER, PUBLIC :: Ap_data => NULL()


  PUBLIC :: SPACENOX_EPPIE
  PUBLIC :: SPACENOX_GCR
  PUBLIC :: SPACENOX_read_nml_ctrl
  PUBLIC :: spacenox_clean

CONTAINS

  ! =========================================================================
  ! *************************************************************************
  ! NOx production from EPP IE (Energetic particle precipitation 
  ! Indirect Effect) - downward transport from the thermosphere
  ! *************************************************************************
  SUBROUTINE SPACENOX_EPPIE (grheight,density      & ! INPUT
       ,philat ,ilat                               & ! INPUT
       ,jrow, kproma                               & ! INPUT
       ,te_no                                      & ! OUTPUT 
       ,status                                     & ! OUTPUT
       )

    USE messy_main_constants_mem, ONLY: M_air & ! mean molar mass of air
                                       ,N_A ! Avogadro constant [1/mol]

    IMPLICIT NONE

    ! Input Parameters
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: density  ! air density midlevel [g cm-3]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: grheight ! height of box [m]
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: philat ! latitude philat(ilat)
    INTEGER,     DIMENSION(:,:), INTENT(IN) :: ilat  ! latitude index (1:kproma,jrow)
    INTEGER,                     INTENT(IN) :: jrow
    INTEGER,                     INTENT(IN) :: kproma 
    ! Output Parameters
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: te_no  ! SPACENOX NO tendency [mol/mol/s]
    INTEGER,                    INTENT(OUT) :: status   ! error status

    ! Local
    REAL(DP) :: NO_prod ! additional NO

    INTEGER :: level,jp ! counters 


    INTRINSIC ABS

    status=1

    ! initialize arrays
    te_no(:,:,jrow)=0._dp

    ! loop over levels
    DO level=1, 2
       ! loop over all rows
       DO jp=1, kproma
          IF (philat(ilat(jp,jrow))>EPPIE_latN) THEN
             NO_prod=(Ap_data(1)**2.5_dp)*2.5E5_dp/(grheight(jp,level)*100._dp)*EPPIE_coeffN*EPPIE_seasoncoeffN ! #/cm3/s
             ! convert #/cm3/s to mole fraction, volume mixing ratio [mol/mol/s]
             NO_prod=(NO_prod/N_A) / (density(jp,level)/M_air)
             ! add tendencies at topmost two layers
             te_no(jp, level, jrow)=.5_dp*NO_prod 
          END IF
          IF (philat(ilat(jp,jrow))<EPPIE_latS) THEN
             NO_prod=(Ap_data(1)**2.5_dp)*2.5E5_dp/(grheight(jp,level)*100._dp)*EPPIE_coeffS*EPPIE_seasoncoeffS ! #/cm3/s
             ! convert #/cm3/s to mole fraction, volume mixing ratio [mol/mol/s]
             NO_prod=(NO_prod/N_A) / (density(jp,level)/M_air)
             ! add tendencies at topmost two layers
             te_no(jp, level, jrow)=.5_dp*NO_prod 
          END IF
       END DO ! rows
    END DO ! levels
    ! No error
    status=0

  END SUBROUTINE SPACENOX_EPPIE


  ! =========================================================================
  ! ****************************************************************
  ! NOx production from GCRs
  ! ****************************************************************
  SUBROUTINE SPACENOX_GCR (density,numdensity      & ! INPUT
       ,philat ,ilat                               & ! INPUT
       ,jrow, kproma, nlev                         & ! INPUT
       ,year                                       & ! INPUT 
       ,te_n, te_no, te_oh                         & ! OUTPUT --> streams
       ,status                                     & ! OUTPUT
       )

    USE messy_main_constants_mem, ONLY: pi  & 
         ,M_air & ! mean molar mass of air
         ,N_A ! Avogadro constant [1/mol]

    IMPLICIT NONE

    ! Input Parameters
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: density  ! air density midlevel [g cm-3]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: numdensity  ! total number density midlevel [# cm-3]
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: philat ! latitude philat(ilat)
    INTEGER,     DIMENSION(:,:), INTENT(IN) :: ilat   ! latitude index (1:kproma,jrow)
    INTEGER,                     INTENT(IN) :: jrow
    INTEGER,                     INTENT(IN) :: kproma 
    INTEGER,                     INTENT(IN) :: nlev     ! number of levels
    INTEGER,                     INTENT(IN) :: year
    ! Output Parameters
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: te_n   ! GCR N tendency [mol/mol/s]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: te_no  ! GCR NO tendency [mol/mol/s]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: te_oh  ! GCR OH tendency [mol/mol/s]
    INTEGER,                    INTENT(OUT) :: status ! error status

    ! Local
    REAL(DP) :: N_prod  ! additional N
    REAL(DP) :: NO_prod ! additional NO
    REAL(DP) :: OH_prod ! additional OH
    REAL(DP) :: sol_cycle_factor,X1, X2, X3, phi, ionpairs

    INTEGER :: level,jp ! counters

    INTRINSIC ABS, SIN, COS

    status=1

    ! initialize arrays
    te_n (:,:,jrow)=0._dp
    te_no(:,:,jrow)=0._dp
    te_oh(:,:,jrow)=0._dp

    sol_cycle_factor=COS(2._dp*pi*(REAL(year,dp)-5.2_dp)/11._dp)
    X1=1.74E-18_dp
    ! X2 is described by the linear expression a*solcyclestate+b :
    X2=4.55E-18_dp*sol_cycle_factor+2.385E-17_dp

    ! loop over all levels
    DO level=1, nlev
       ! loop over all rows
       DO jp=1, kproma
          phi=philat(ilat(jp,jrow))/180._dp*pi ! latitude in rad
          X3=0.6_dp+0.8_dp*ABS(COS(phi))
          IF (ABS(philat(ilat(jp,jrow))) <= 53._dp) THEN

             IF (numdensity(jp,level) >= 3E17_dp) THEN
                ionpairs = (X1+X2*ABS(SIN(phi))**4) &
                     * (3E17_dp)**(1._dp-X3)*(numdensity(jp,level)**X3)
             ELSE
                ionpairs = (X1+X2*ABS(SIN(phi))**4) &
                     * (numdensity(jp,level))
             END IF
             N_prod =0.55_dp*ionpairs ! #/cm3/s
             NO_prod=0.7_dp *ionpairs ! #/cm3/s
             OH_prod=2._dp  *ionpairs ! #/cm3/s

             ! convert #/cm3/s to mole fraction, volume mixing ratio [mol/mol/s]
             te_n (jp, level, jrow)=N_prod /(N_A*density(jp,level))*M_air
             te_no(jp, level, jrow)=NO_prod/(N_A*density(jp,level))*M_air
             te_oh(jp, level, jrow)=OH_prod/(N_A*density(jp,level))*M_air

          ELSE ! latitudes greater than 53 deg
             ! X3=X1 for solar max, X3=X1+X2 for solar min, expressed as linear expression
             X3=2.46E-18_dp*sol_cycle_factor+1.686E-17_dp

             ! shortcut: Because numdensity = density*N_A/M_air there is no altitude dependence
             te_n (jp, level, jrow)=0.55_dp*X3 
             te_no(jp, level, jrow)=0.7_dp *X3 
             te_oh(jp, level, jrow)=2._dp  *X3 

             ! Write(*,*) 'greater53 ', level, ',',jp,',',philat(ilat(jp,jrow)), ',', &
             ! ionpairs, ',',numdensity(jp,level),',',te_oh(jp, level, jrow)
          END IF ! latitude dependence


       END DO ! rows
    END DO ! levels
    ! No error
    status=0

  END SUBROUTINE SPACENOX_GCR
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE SPACENOX_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/  EPPIE_latN, EPPIE_latS &
         ,EPPIE_coeffN, EPPIE_coeffS 

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='spacenox_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE SPACENOX_read_nml_ctrl
  ! =========================================================================

  !--------------------------------------------------------------------------
  SUBROUTINE spacenox_clean

    IMPLICIT NONE

    INTRINSIC ASSOCIATED

    ! DATA
    IF (ASSOCIATED(Ap_data)) DEALLOCATE(Ap_data)
    NULLIFY(Ap_data)

  END SUBROUTINE spacenox_clean
  !--------------------------------------------------------------------------

  ! ***********************************************************************
END MODULE messy_spacenox
! ***********************************************************************
