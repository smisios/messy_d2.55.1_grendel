MODULE messy_m7

  ! MODULE FOR M7 CORE
  !
  ! Swen Metzger    (metzger@mpch-mainz.mpg.de), MPI-CHEM, Dec 2003
  ! Astrid Kerkweg (akerkweg@mpch-mainz.mpg.de), MPI-CHEM, Dec 2003
  ! M7 adopted to the structure of the Modular Earth Submodel System (MESSy)
  ! M7 was originally implemented by Philip Stier, MPI-Met, Hamburg, 2001-2003
  ! Original M7 source code (box model) by J. Wilson, 
  !   E. Vignati, JRC/EI, 09/2000

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, SP, I4, I8, avo => N_A

  IMPLICIT NONE
  PRIVATE
  SAVE
  
  PUBLIC :: DP, SP, I4, I8

  ! ----------- <

  INTRINSIC :: EXP, LOG, MAX, MIN, SQRT, EPSILON, ABS, AINT
  PRIVATE   :: EXP, LOG, MAX, MIN, SQRT, EPSILON, ABS, AINT

  ! GLOBAL PARAMETER
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'm7'    ! name of module
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.1'   ! module version

  ! CTRL-NAMELIST PARAMETERS
  LOGICAL, PUBLIC :: lm7        = .FALSE.  ! Aerosol dynamics and 
                                   ! thermodynamics scheme M7
  LOGICAL, PUBLIC :: lmass_diag = .FALSE.  ! Mass balance check in m7_interface 
  LOGICAL, PUBLIC :: lcdnc      = .FALSE.  ! Interactive calculation of the Cloud
                                           !  Droplet Number Concentration
  LOGICAL, PUBLIC :: licnc      = .FALSE.  ! Interactive calculation of the Ice 
                                           !  Crystal   Number Concentration
  LOGICAL, PUBLIC :: lsnucl     = .FALSE.  ! Nucleation
  LOGICAL, PUBLIC :: lscoag     = .FALSE.  ! Coagulation
  LOGICAL, PUBLIC :: lscond     = .FALSE.  ! Condensation of H2SO4
  INTEGER(i4), PUBLIC :: nnucl  = 1        ! Choice of the nucleation scheme:
                                           !    nnucl = 1  Vehkamaeki (2002)
                                           !          = 2  Kulmala (1998)

  ! GLOBAL PARAMETERS ========================================================
  !
  !--- 1) Numbers of compounds and modes of m7:
  INTEGER, PUBLIC,PARAMETER :: naermod=18,        & !number of all compounds
                            nmod=7,            & !number of modes
                            nss=2,             & !number of sea salt compounds 
                            nsol=4,            & !number of soluble  compounds 
                            ngas=3,            & !number of gaseous  compounds 
                            nsulf=4              !number of sulfate  compounds 
  !--- 2) List of indeces corresponding to the compound masses and
  !       mode numbers:
  !--- 2.1) Mass index (in array aerml and ttn): 
  !
  !         Attention:
  !         The mass of sulfate compounds is always given in [molec. cm-3] 
  !         whilst the mass of other compounds is given in [ug cm-3].
  !
  !         Compounds:
  !
  !           so4 = sulphate
  !           bc  = black carbon
  !           oc  = organic carbon, 
  !           ss  = sea salt
  !           du  = dust 
  !
  !         Modes:
  !
  !           n   = nucleation mode
  !           k   = aitken mode 
  !           a   = accumulation mode
  !           c   = coarse mode
  !
  !         Type:
  !
  !           s   = soluble mode
  !           i   = insoluble mode
  !                                                                  
  !  COMPOUND:
  INTEGER(i4), PUBLIC, PARAMETER ::                                   &
       iso4ns=1, iso4ks=2, iso4as=3, iso4cs=4,  & !- Sulfate
       ibcks =5, ibcas =6, ibccs =7, ibcki =8,  & !- Black Carbon
       iocks =9, iocas=10, ioccs=11, iocki=12,  & !- Organic Carbon
       issas=13, isscs=14,                      & !- Sea Salt
       iduas=15, iducs=16, iduai=17, iduci=18     !- Dust  
  ! MODE:      |         |         |         |         |
  !    nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !    soluble | soluble | soluble | soluble | insol.  | insol.  | insol  .|
  !
  !
  !--- 2.2) Number index (in array aernl):
  !
  INTEGER(i4), PUBLIC, PARAMETER ::                                      &
       inucs=1,  iaits=2,  iaccs=3,  icoas=4,  iaiti=5,  iacci=6,  icoai=7    
  ! MODE:     |         |         |         |         |
  !   nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !   soluble | soluble | soluble | soluble | insol.  | insol.  | insol.  |
  !
  !
  !--- 3) Definition of the modes of M7:
  !--- 3.1) Threshold radii between the different modes [cm]:
  !         Used for the repartititioning in m7_dconc.
  !
  REAL(dp) :: crdiv(4)=(/ 0.0005e-4_dp, 0.005e-4_dp, 0.05e-4_dp, 0.5e-4_dp /)
  !                                   |            |           |      
  !                                   |            |           |
  !                       nucleation - - aitken   - - accum   - - coarse mode
  !---3.2) Standard deviation for the modes:
  !
  REAL(dp), PUBLIC, PARAMETER :: sigma(nmod)=&
       (/ 1.59_dp, 1.59_dp, 1.59_dp, 2.00_dp, 1.59_dp, 1.59_dp, 2.00_dp /)
  !
  !
  !--- 4.) Conversion factors for lognormal particle size distributions:
  !        Calulated in m7_initialize_core. 
  !
  ! Conversion factor: count median radius to radius of average surface
  REAL(dp)            :: cmr2ras(nmod) 
  ! Conversion factor: count median radius to radius of average mass
  REAL(dp), PUBLIC     :: cmr2ram(nmod) 
  ! Conversion factor: radius of average mass to count median radius
  REAL(dp)            :: ram2cmr(nmod) 
  !
  !--- 5) Assumed thresholds for occurence of specific quantities:
  REAL(dp), PARAMETER :: cmin_aerml     = 1.E-12_dp  ! Aerosol mass
  REAL(dp), PARAMETER :: cmin_aernl     = 1.E-7_dp   ! Aerosol number minimum
  !                     
  !--- 6) Chemical constants:
  !
  !--- 7) Accomodation coefficient of H2SO4 on aerosols:
  !       (reduced for insoluble modes)
  REAL(dp), PARAMETER :: caccso4(nmod) = &
       (/ 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.3_dp, 0.3_dp, 0.3_dp /)
  !
  !--- 8) Critical relative humidity:
  !       Assumed relative humidity for the 
  !       Na2SO4 / NaCl system below which crystalization occurs.
  !       (estimated from Tang, I.N.; JGR 102, D2 1883-1893)
  REAL(dp), PARAMETER :: crh    = 0.45_dp
  !
  !--- 9) Physical constants: -------------------------------------------------
  !
  !--- 9.1) General physical constants: 
  !
  REAL(dp), PARAMETER :: &
       bk      = 1.38e-16_dp,     & ! Bolzman constant []
       rerg    = 8.314E+7_dp,     & ! Ideal gas constant [erg.K-1.mole-1]
       r_kcal  = 1.986E-3_dp        ! Ideal gas constant [kcal K-1.mole-1]
  !
  !--- 9.2) Type specific physical constants:
  !
  REAL(dp), PARAMETER ::  &
       dh2so4  = 1.841_dp,      & ! Density          H2SO4  [g cm-3]
       ddust   = 2.650_dp,      & ! Density          du     [g cm-3]
       dbc     = 2._dp,         & ! Density          bc     [g cm-3]
       doc     = 2._dp,         & ! Density          oc     [g cm-3]
       dnacl   = 2.165_dp,      & ! Density          NaCl   [g cm-3]
       dna2so4 = 2.68_dp,       & ! Density          Na2SO4 [g cm-3]
       dnahso4 = 2.435_dp,      & ! Density          NaHSO4 [g cm-3]
       dh2o    = 1.0_dp,        & ! Density          H2O    [g cm-3]
       wh2so4  = 98.0734_dp,    & ! Molecular weight H2SO4  [g mol-1]
       wnacl   = 58.443_dp,     & ! Molecular weight NaCl   [g mol-1]
       wna2so4 = 142.0376_dp,   & ! Molecular weight Na2SO4 [g mol-1]
       wnahso4 = 120.0555_dp      ! Molecular weight NaHSO4 [g mol-1]
  !
  !--- 10) Assumed parameters:
  !        Assumed mass of an nucleated sulfate particle [molecules]
  REAL(dp), PARAMETER :: critn=100._dp
  !        Factor that limits the condensation of sulfate to fmax
  !        times the available sulfate in the gas phase [1]. (m7_dgas)
  REAL(dp), PARAMETER :: fmax=0.95_dp
  !         Assumed required layer thickness of
  !         sulfate to transfer an insoluble 
  !         particle to a soluble mode. It is
  !         given in units of layers of 
  !         monomolecular sulfate. Determines the
  !         transfer rate from insoluble to soluble modes. 
  REAL(dp), PARAMETER :: cLayerThickness = 1.0_dp
  !
  !--- 11) Computational constants: -------------------------------------------
  REAL(dp), PARAMETER :: sqrt2=1.4142136_dp,    & 
       pi=3.141592654_dp
  !
  !--- 12) Data used for the calculation of the aerosol properties ------------
  !        under ambient conditions:
  !       (Included the conversion from Pa to hPa in the first parameter.)
  REAL(dp), PARAMETER :: wvb(17)=                                          &
       (/   95.80188_dp,  -28.5257_dp,  -1.082153_dp,    0.1466501_dp,     &
       -20627.51_dp,    0.0461242_dp,   -0.003935_dp,     -3.36115_dp,     &
       -0.00024137_dp,  0.067938345_dp, 0.00000649899_dp,  8616124.373_dp, &
       1.168155578_dp, -0.021317481_dp,  0.000270358_dp, -1353332314.0_dp, &
       -0.002403805_dp /)
  !
  REAL(dp), PARAMETER :: gmb(9)=                                   &
       (/ 1.036391467_dp, 0.00728531_dp, -0.011013887_dp, -0.068887407_dp, &
       0.001047842_dp, 0.001049607_dp, 0.000740534_dp, -1.081202685_dp,    &
       -0.0000029113_dp /)

  !---    Logical mask for coagulation kernel: -------------------------------------
  !       (The coagulation kernel mask is symmetric and not all 
  !        values are used for physical considerations. As its 
  !        calculation is very expensive, a mask is used to 
  !        calculate only the necessarey elements.)
  LOGICAL :: locoagmask(nmod,nmod)

  DATA locoagmask(1:nmod,1) / .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.  /

  DATA locoagmask(1:nmod,2) / .FALSE., .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.  /

  DATA locoagmask(1:nmod,3) / .FALSE., .FALSE., .TRUE.,  .FALSE., .TRUE.,  .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,4) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,5) / .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.,  .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,6) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,7) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /

  ! END GLOBAL PARAMETERS ====================================================

  ! SUBROUTINES
  PUBLIC :: m7_read_nml_ctrl
  PUBLIC :: m7_initialize_core
  PUBLIC :: m7_main
  !
  !PRIVATE m7_averageproperties
  !PRIVATE m7_equiz
  !PRIVATE m7_equimix
  !PRIVATE m7_equil 
  !PRIVATE m7_dgas 
  !PRIVATE m7_coaset 
  !PRIVATE m7_nuck
  !PRIVATE m7_delcoa
  !PRIVATE m7_concoag
  !PRIVATE m7_coat
  !PRIVATE m7_dnum
  !PRIVATE m7_cumnor
  !PRIVATE m7_dconc
  !PRIVATE nucl_kulmala
  !PRIVATE nucl_vehkamaeki

CONTAINS

! -------------------------------------------------------------------------
 SUBROUTINE m7_read_nml_ctrl(status, iou)

    ! M7  MODULE ROUTINE (CORE)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Sep 2003
    ! Author: Swen Metzger,    MPICH, Oct 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER(i4), INTENT(OUT) :: status ! error status
    INTEGER,     INTENT(IN)  :: iou    ! logical I/O unit

    !--- Local variables:
    NAMELIST /CTRL/ lm7,lmass_diag,lcdnc,licnc,lscoag,lscond,lsnucl,nnucl

    CHARACTER(LEN=*), PARAMETER :: substr = 'm7_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE

    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    !--- 1) Read namelist:
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! consistency checks and diagnostic outpur
    WRITE(*,*) ''
    WRITE(*,*) ''
    WRITE(*,*) '----------------------------------------------------------'
    WRITE(*,*) '----------------------------------------------------------'
    WRITE(*,*) '--- Initialization of settings for aerosol module M7   ---'
    WRITE(*,*) '---'
    WRITE(*,*) '---    Default values of aeroctl modified by setaero:'
    WRITE(*,*) '---'
    WRITE(*,*) '---    New settings: lm7    = ', lm7
    WRITE(*,*) '---                  M7     = ', modstr,modver
    WRITE(*,*) '---                            => Wilson et al. (2001)'
    IF (lmass_diag) THEN
       WRITE(*,*) '---    Mass balance check in m7_interface activated'
    ELSE
       WRITE(*,*) '---    Mass balance check in m7_interface deactivated'
    END IF
    WRITE(*,*) '---'
    WRITE(*,*) '---                  lcdnc = ', lcdnc
    WRITE(*,*) '---                  licnc = ', licnc
    WRITE(*,*) '---'
    WRITE(*,*) '---                  lscoag = ', lscoag
    WRITE(*,*) '---                  lscond = ', lscond
    WRITE(*,*) '---                  lsnucl = ', lsnucl
    IF (nnucl==1) THEN
       WRITE(*,*) '---                  nnucl  = ', nnucl
       WRITE(*,*) '---                            => Vehkamaeki et al., 2002'
    ELSE IF (nnucl==2)  THEN
       WRITE(*,*) '---                  nnucl  = ', nnucl
       WRITE(*,*) '---                            => Kulmala et al., 1998'
    ELSE IF (lsnucl .AND. (nnucl/=1 .OR. nnucl/=2)) THEN
       lsnucl=.FALSE.
       WRITE(*,*) 'm7_read_nml_ctrl: nucleation requested but no scheme selected'
       WRITE(*,*) 'nucleation switched off !'
    END IF
    WRITE(*,*) '----------------------------------------------------------'
    WRITE(*,*) ''
 
    !--- Set the lower mode boundary for the nucleation mode:
    !       (Depends on the choice of the nucleation scheme.)

    SELECT CASE (nnucl) 
    CASE(1)
       crdiv(1)=0.0005e-4_dp

    CASE(2)
       crdiv(1)=( critn * wh2so4/avo/dh2so4*0.75_dp/pi )**(1._dp/3._dp)
    END SELECT

    CALL read_nml_close(substr, iou, modstr)
 
    status = 0  ! no ERROR

  END SUBROUTINE m7_read_nml_ctrl
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE m7_initialize_core

    ! Purpose:
    ! ---------
    ! Initializes constants and parameters 
    ! used in the m7 aerosol model.
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                          may 2001
    !
    ! Interface:
    ! ---------
    ! *m7_initialize_core* is called from *m7_init* in *messy_m7_e5*
    !

    IMPLICIT NONE

    INTEGER(i4) :: jmod
      
    DO jmod=1, nmod

       !--- 1) Calculate conversion factors for lognormal distributions:----
       !       Radius of average mass (ram) to count median radius (cmr) and 
       !       vice versa. Count median radius to radius of average 
       !       mass (ram).
       !       These factors depend on the standard deviation (sigma)
       !       of the lognormal distribution.
       !       (Based on the Hatch-Choate Conversins Equations;
       !        see Hinds, Chapter 4.5, 4.6 for more details.
       !        In particular equation 4.53.)

       !--- Count Median Radius to Radius of Average Mass:

       cmr2ram(jmod) = EXP(1.5_dp*(LOG(sigma(jmod)))**2)

       !--- Radius of Average Mass to Count Median Radius:

       ram2cmr(jmod) = 1._dp / cmr2ram(jmod)

       !--- Count Median Radius to Radius of Average Surface:

       cmr2ras(jmod) = EXP(1.0_dp*(LOG(sigma(jmod)))**2)


    END DO


  END SUBROUTINE m7_initialize_core
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE m7_main(kproma, kbdim,   klev, dt,   &  ! ECHAM indices
              papp1, prelhum, ptp1,           &  !   "   thermodynamics
              pso4g,          paerml, paernl, &  !  M7   tracers
              pm6rp, pm6dry,  prhop,  pww      ) !   "   aerosol properties
  !
  !   ****m7* Aerosol model for the system so4,bc,oc,ss,dust in 7 modes.
  !
  !   Authors:
  !   ---------
  !   E. Vignati, JRC/EI     (original source)                   2000
  !   P. Stier, MPI          (f90-version, changes, comments)    2001 
  !
  !   Purpose
  !   ---------
  !   Aerosol model for the system so4,bc,oc,ss,dust in 7 modes.
  !
  !   Interface:
  !   ----------
  !   m7_main is called from *messy_m7_box or messy_m7_e5*
  !
  !   ! private routines
  !   *m7_averageproperties* 
  !       calculates the average mass for all modes and the particle 
  !       dry radius and density for the insoluble modes.
  !    
  !   *m7_equiz*   
  !       calculates the ambient radius of sulphate particles
  !
  !   *m7_equimix* 
  !       calculates the ambient radius of so4,bc,oc (dust) particles
  !
  !   *m7_equil* 
  !       calculates the ambient radius of so4,ss particles 
  !
  !   *m7_dgas*    
  !       calculates the sulfate condensation on existing particles
  !
  !   *m7_dnum*    
  !       calculates new gas phase sulfate and aerosol numbers and masses
  !       after condensation, nucleation and coagulation over one timestep
  !
  !   *m7_dconc*   
  !       repartitions aerosol number and mass between the
  !       the modes to account for condensational growth and the formation
  !       of an accumulation mode from the upper tail of the aitken mode and 
  !       of a coarse mode from the upper tail of the accumulation mode
  !
  !--- Parameter list:
  !
  !  papp1      = atmospheric pressure at time t+1 [Pa]
  !  prelhum    = atmospheric relative humidity [% (0-1)]
  !  ptp1       = atmospheric temperature at time t+1 [K]
  !  pso4g      = mass of gas phase sulfate [molec. cm-3]
  !  paerml     = total aerosol mass for each compound 
  !               [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and dust]
  !  paernl     = aerosol number for each mode [cm-3]
  !  prhop      = mean mode particle density [g cm-3]
  !  pm6rp      = mean mode actual radius (wet radius for soluble modes 
  !               and dry radius for insoluble modes) [cm]
  !  pm6dry     = dry radius for soluble modes [cm]
  !  pww        = aerosol water content for each mode [kg(water) m-3(air)]
  !
  !--- Local variables:
  !
  !  zttn       = average mass for single compound in each mode 
  !               [in molec. for sulphate and in ug for bc, oc, ss, and dust]
  !  zhplus     = number of h+ in mole [???] (kg water)-1
  !  zso4_x     = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  !               (calculated in dgas used in concoag)

  ! Parameters:

    IMPLICIT NONE

  INTEGER  :: kproma, kbdim, klev
  REAL(dp) :: dt
  REAL(dp) :: prelhum(kbdim,klev), papp1(kbdim,klev),           &
              ptp1(kbdim,klev), pso4g(kbdim,klev)
 
  REAL(dp)    :: paerml(kbdim,klev,naermod), paernl(kbdim,klev,nmod),     &
             pm6rp(kbdim,klev,nmod),     pm6dry(kbdim,klev,nsol),     &
             prhop(kbdim,klev,nmod),     pww(kbdim,klev,nmod)

  ! Local variables:

  REAL(dp)    :: zso4_5(kbdim,klev),         zso4_6(kbdim,klev),          &
             zso4_7(kbdim,klev)

  REAL(dp)    :: zhplus(kbdim,klev,nss)

  REAL(dp)    :: zttn(kbdim,klev,naermod)

  !
  !--- 0) Initialisations: -------------------------------------------------
  !
  zhplus(:,:,:) = 0._dp
  pm6dry(:,:,:) = 0._dp
  pm6rp(:,:,:)  = 0._dp
  zttn(:,:,:)   = 0._dp
  prhop(:,:,:)  = 0._dp
  pww(:,:,:)    = 0._dp 
  zso4_5(:,:)   = 0._dp
  zso4_6(:,:)   = 0._dp
  zso4_7(:,:)   = 0._dp
  !
  !
  !--- 1) Calculation of particle properties under ambient conditions: -----
  !
  !--- 1.1) Calculate mean particle mass for all modes 
  !         and dry radius and density for the insoluble modes.
  !
!CDIR NOIEXPAND
  CALL m7_averageproperties(kproma, kbdim, klev, paernl, paerml, zttn, &
       pm6rp, prhop)
  !
  !--- 1.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  !
  !         Sulfate particles:
  !
!CDIR NOIEXPAND
  CALL m7_equiz(kproma,  kbdim, klev,   &
                papp1,   zttn,  ptp1,   &
                prelhum, pm6rp, pm6dry, &
                prhop,   pww,   paernl  )
  !         
  !         Mixed particles with sulfate, b/o carbon and dust: 
  !
!CDIR NOIEXPAND
  CALL m7_equimix(kproma,  kbdim, klev,   &
                  papp1,   zttn,  ptp1,   &
                  prelhum, pm6rp, pm6dry, &
                  prhop,   pww,   paernl  )
  !
  !         Accumulation and coarse mode particles in presence of
  !         sea salt particles:
  !
!CDIR NOIEXPAND
  CALL m7_equil(kproma, kbdim,  klev,   prelhum, paerml, paernl, &
                pm6rp,  pm6dry, zhplus, pww,     prhop,  ptp1    )
  !
  !--- 2) Calculate changes in aerosol mass and gas phase sulfate ----------
  !       due to sulfate condensation:
  !       No change in particle mass/number relationships.
  !
!CDIR NOIEXPAND
  IF (lscond) CALL m7_dgas(kproma, kbdim,  klev,  pso4g,  paerml, paernl, & 
                           ptp1,   papp1,  pm6rp, zso4_5, zso4_6, zso4_7, &
                           dt                                  )
  !
  !
  !--- 3) Calculate change in particle number concentrations ---------------
  !       due to nucleation and coagulation:
  !       Change particle mass/number relationships.
  !
!CDIR NOIEXPAND
  IF (lsnucl.OR.lscoag) CALL m7_dnum(kproma, kbdim,   klev,          &
                                     pso4g,  paerml,  paernl, ptp1,  &
                                     papp1,  prelhum, pm6rp,  prhop, &
                                     zso4_5, zso4_6,  zso4_7, dt) !qqq
  !
  !
  !--- 4) Recalculation of particle properties under ambient conditions: ---
  !
  !--- 4.1) Recalculate mean masses for all modes 
  !         and dry radius and density for the insoluble modes.
  !

!CDIR NOIEXPAND
  CALL m7_averageproperties(kproma, kbdim, klev, paernl, paerml, zttn, pm6rp, &
       prhop)
  !
  !--- 4.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  !
  !         Sulfate particles:
  !
!CDIR NOIEXPAND
  CALL m7_equiz(kproma,  kbdim, klev,   &
                papp1,   zttn,  ptp1,   &
                prelhum, pm6rp, pm6dry, &
                prhop,   pww,   paernl  )
  !
  !         Mixed particles with sulfate, b/o carbon and dust:
  !
!CDIR NOIEXPAND
  CALL m7_equimix(kproma,  kbdim, klev,   &
                  papp1,   zttn,  ptp1,   &
                  prelhum, pm6rp, pm6dry, &
                  prhop,   pww,   paernl  )
  !
  !         Accumulation and coarse mode particles in presence of
  !         sea salt particles:  
  ! 
!CDIR NOIEXPAND
  CALL m7_equil(kproma, kbdim,  klev,   prelhum, paerml, paernl,  &
                pm6rp,  pm6dry, zhplus, pww,     prhop,  ptp1     )
  !
  !--- 5) Repartitition particles among the modes: -------------------------
  !
  IF (lscond.OR.lscoag) THEN
     
!CDIR NOIEXPAND
     CALL m7_dconc(kproma, kbdim, klev, paerml, paernl, pm6dry) 
     
  END IF
  !
  !--- 6) Recalculation of particle properties under ambient conditions: ---
  !
  !--- 6.1) Calculate mean particle mass for all modes 
  !         and dry radius and density for the insoluble modes:
  !
!CDIR NOIEXPAND
  CALL m7_averageproperties(kproma, kbdim, klev, paernl, paerml, zttn, pm6rp, &
       prhop)
  !
  !--- 6.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  !
  !         Sulfate particles:
  !
!CDIR NOIEXPAND
  CALL m7_equiz(kproma,  kbdim, klev,   &
                papp1,   zttn,  ptp1,   &
                prelhum, pm6rp, pm6dry, &
                prhop,   pww,   paernl  )
  !
  !         Mixed particles with sulfate, b/o carbon and dust: 
  !
!CDIR NOIEXPAND
  CALL m7_equimix(kproma,  kbdim, klev,   &
                  papp1,   zttn,  ptp1,   &
                  prelhum, pm6rp, pm6dry, &
                  prhop,   pww,   paernl  )
  !
  !         Accumulation and coarse mode particles in presence of
  !         sea salt particles:
  !
!CDIR NOIEXPAND
  CALL m7_equil(kproma, kbdim,  klev, prelhum, paerml, paernl, &
                pm6rp,  pm6dry, zhplus, pww,   prhop,  ptp1    )
  !
END SUBROUTINE m7_main
! -------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! PRIVATE ROUTINES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

SUBROUTINE m7_averageproperties(kproma, kbdim, klev, paernl, paerml, pttn, &
     pm6rp, prhop)
  !
  !  Author:
  !  --------
  !  E. Vignati, JRC/EI (original source)                10/2000
  !  P. Stier, MPI      (f90-version, changes, comments)    2001
  !
  !  Purpose:                                                           
  !  ---------                                                           
  !  Calculation of the mean particle mass (pttn).
  !     [molecules cm-3] for the sulphate mass
  !     [ug m-3] for the other compounds
  !
  !  Calculation of the (dry) radius and the density 
  !  of the particles of the insoluble modes.
  !
  !  Interface:
  !  ----------
  !  m7_averageproperties is called from m7
  !
  !  Externals:
  !  ----------
  !  none
  !
  !--- Parameter list:
  !
  ! paerml(kbdim,klev,naermod) = total aerosol mass for each compound 
  !                             [molec. cm-3 for sulfate and ug m-3 for others]
  ! paernl(kbdim,klev,nmod)    = aerosol number for each mode [cm-3]
  ! pttn(kbdim,klev,naermod)   = average mass for single compound in each mode 
  !                             [in molec. for sulphate and in ug for others]
  ! pm6rp(kbdim,klev,nmod)     = mean mode actual radius (wet radius for soluble
  !                             modes and dry radius for insoluble modes) [cm]
  ! prhop(kbdim,klev,nmod)     = mean mode particle density [g cm-3]
  !
  !--- Local variables:
  !
  ! zinsvol                   = average volume for single particle in the 
  !                             insolulbe mode [cm3]
  ! zinsmas                   = average mass for single particle in the 
  !                             insolulbe mode [g]

  !--- Parameters:

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: paerml(kbdim,klev,naermod), paernl(kbdim,klev,nmod), & 
             pttn(kbdim,klev,naermod),   pm6rp(kbdim,klev,nmod),  & 
             prhop(kbdim,klev,nmod)
  
  !--- Local variables:

  INTEGER(i4) :: jmod,         jk,          jl
  
  REAL(dp)    :: zinsvol,                   zinsmas

  REAL(dp)    :: zeps
  
  !--- 0) Initialization:

  zeps=EPSILON(1.0_dp)

  !
  !--- 1) Calculate mean particle masses at start of timestep: ----------------
  !
  ! To be able to compute a intra-modal coagulation coefficient 
  ! for the nucleation
  ! mode for the case of no pre-existing particles but coagulation of freshly 
  ! formed particles during the timestep, pttn is set to the mass of the 
  ! critical cluster
  ! for this case. This allows to calculate an ambient radius of the 
  ! freshly formed particles and subsequently the calculation of the coagulation
  ! coefficient. This mass is "virtual" as it is not added to the mode but used 
  ! only for the described computation of the coagulation coefficient. 
  ! !@@@ Check whether this is always fulfilled. 
  
  DO jmod=1,nsol
     DO jk=1,klev
        DO jl=1,kproma
           IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,jmod)&
                .GT. cmin_aerml) THEN

              pttn(jl,jk,jmod)=paerml(jl,jk,jmod)/paernl(jl,jk,jmod)

           ELSE IF (jmod == 1 .AND. paernl(jl,jk,jmod) <= cmin_aernl .AND.&
                paerml(jl,jk,jmod) <= cmin_aerml) THEN

              pttn(jl,jk,jmod)=critn

           END IF
        END DO
     END DO
  END DO
  !
  !-- 3) Calculation of the mean mass pttn [ug] for each compound in the modes:
  !       [Factor 1.E-6 to convert(ug m-3)/cm-3 into ug]
  !
  DO jmod=2,nmod
     DO jk=1,klev
        DO jl=1,kproma
           IF (jmod.EQ.2) THEN
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,ibcks)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,ibcks)=paerml(jl,jk,ibcks)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ibcks)=0._dp
              END IF
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,iocks) &
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,iocks)=paerml(jl,jk,iocks)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iocks)=0._dp
              END IF
           END IF
           !        
           IF (jmod.EQ.3) THEN
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,ibcas)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,ibcas)=paerml(jl,jk,ibcas)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ibcas)=0._dp
              END IF
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,iocas)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,iocas)=paerml(jl,jk,iocas)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iocas)=0._dp
              END IF
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl.AND. paerml(jl,jk,issas) &
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,issas)=paerml(jl,jk,issas)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,issas)=0._dp
              END IF
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,iduas) &
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,iduas)=paerml(jl,jk,iduas)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iduas)=0._dp
              END IF
           END IF
           !        
           IF (jmod.EQ.4) THEN
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,ibccs)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,ibccs)=paerml(jl,jk,ibccs)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ibccs)=0._dp
              END IF
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,ioccs)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,ioccs)=paerml(jl,jk,ioccs)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ioccs)=0._dp
              END IF
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,isscs)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,isscs)=paerml(jl,jk,isscs)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,isscs)=0._dp
              END IF
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,iducs)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,iducs)=paerml(jl,jk,iducs)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iducs)=0._dp
              END IF
           END IF
           !        
           IF (jmod.EQ.5) THEN
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,ibcki)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,ibcki)=paerml(jl,jk,ibcki)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ibcki)=0._dp
              END IF
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,iocki)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,iocki)=paerml(jl,jk,iocki)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iocki)=0._dp
              END IF
           END IF
           !        
           IF (jmod.EQ.6) THEN
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,iduai)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,iduai)=paerml(jl,jk,iduai)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iduai)=0._dp
              END IF
           END IF
           !
           IF (jmod.EQ.7) THEN
              IF (paernl(jl,jk,jmod) .GT. cmin_aernl .AND. paerml(jl,jk,iduci)&
                   .GT. cmin_aerml) THEN
                 pttn(jl,jk,iduci)=paerml(jl,jk,iduci)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iduci)=0._dp
              END IF
           END IF
        END DO
     END DO
  END DO

  !
  !--- 4) Calculate count median radii for lognormal distribution from --------
  !       mass for insoluble modes:

  DO jk=1,klev
     DO jl=1,kproma

        !--- 4.1) Aitken mode insoluble:

        zinsmas=1.e-6*(pttn(jl,jk,ibcki)+pttn(jl,jk,iocki))
        zinsvol=1.e-6*(pttn(jl,jk,ibcki)/dbc+pttn(jl,jk,iocki)/doc)
        IF (zinsvol > zeps) THEN
           prhop(jl,jk,iaiti)=zinsmas/zinsvol
           pm6rp(jl,jk,iaiti)=(0.75_dp/pi*1.e-6*    &
                (pttn(jl,jk,ibcki)/dbc+pttn(jl,jk,iocki)/doc))**(1._dp/3.)&
                *ram2cmr(iaiti)
        ELSE
           prhop(jl,jk,iaiti)=0._dp
           pm6rp(jl,jk,iaiti)=0._dp
        END IF

        !--- 4.2) Accumulation mode insoluble:

        IF (pttn(jl,jk,iduai) > zeps) THEN 
           prhop(jl,jk,iacci)=ddust
           pm6rp(jl,jk,iacci)=(0.75_dp/pi*1.e-6&
                *pttn(jl,jk,iduai)/ddust)**(1._dp/3.)*ram2cmr(iacci)
        ELSE
           prhop(jl,jk,iacci)= 0._dp
           pm6rp(jl,jk,iacci)= 0._dp
        END IF

        !--- 4.3) Coarse mode insoluble:

        IF (pttn(jl,jk,iduci) > zeps) THEN
           prhop(jl,jk,icoai)=ddust
           pm6rp(jl,jk,icoai)=(0.75_dp/pi*1.e-6 &
                *pttn(jl,jk,iduci)/ddust)**(1._dp/3._dp)*ram2cmr(icoai)
        ELSE
           prhop(jl,jk,icoai)=0._dp
           pm6rp(jl,jk,icoai)=0._dp
        END IF

     END DO
  END DO

END SUBROUTINE m7_averageproperties
!------------------------------------------------------------------------------!
SUBROUTINE m7_equiz(kproma,  kbdim, klev,   &
                    papp1,   pttn,  ptp1,   &
                    prelhum, pm6rp, pm6dry, &
                    prhop,   pww,   paernl  )
  !
  !   *m7_equiz*   calculates the ambient radii of the sulphate particles
  !
  !    Authors:
  !    --------
  !    J. Wilson, E. Vignati, JRC/EI (original source)                05/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)    2001
  !
  !    Purpose:
  !    --------
  !    This routine calculates the ambient radii for sulfate particles
  !    with mass of ttn molecules, converts them to count mean radii and
  !    stores them in the array with address pm6rp.
  !    It additionally calculates the ambient particle density.
  !
  !    Method:
  !    -------
  !    The calculations of the ambient particle properties are based on 
  !    parameterisations of the mass of sulfate and density derived 
  !    by Julian Wilson from a regression analysis of results of solving 
  !    the generalised Kelvin equation using (F. J. Zeleznik, J. Phys. Chem.
  !    Ref. Data 20, 1157, 1991), for an H2SO4-H2O mixture, in the 
  !    following parameter ranges: 

  !       1e2   < pttn    < 1E11  [molecules]
  !       0.2   < prelhum < 0.9   [1]
  !       240   < ptp1    < 330   [K]
  !       10000 < papp1   < 100000  [Pa]
  !
  !    Due to the limitations of the parametrisation, the ambient temperature
  !    is restricted to a minimum of 240 K within this subroutine. 
  !      
  !    Interface:
  !    ----------
  !    *m7_equiz* is called from *m7* and *m7_dconc*
  !
  !    Externals:
  !    ----------
  !    none
  !
  !
  ! pttn      = average mass for single compound in each mode 
  !             [in molec. for sulphate and in ug for bc, oc, ss, and dust]
  ! pm6rp     = count mean radius under ambient conditions [cm]
  ! pm6dry    = count mean radius under dry conditions [cm]
  ! paernl    = aerosol number for each mode [cm-3]
  ! pww       = aerosol water content for each mode [kg(water) m-3(air)]
  ! zwso4     = percentage by mass of sulfate in a H2O-H2SO4 particle 
  !             containing pttn molecules of sulfate under ambient conditions
  ! zvso4     = volume of pttn molecules of sulfate [cm3]
  ! zmso4     = mass of pttn molecules of sulfate [g]
  ! zdso4h2o  = density of sulfate-h2o fraction of a particle with average 
  !             mass [g.cm-3]
  ! zmso4h2o  = mass of sulfate-h2o fraction of a particle with average mass [g]
  ! zvso4h2o  = volume of sulfate-h2o fraction of a particle with average 
  !             mass [cm3]

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: papp1(kbdim,klev),        ptp1(kbdim,klev),       &
             prelhum(kbdim,klev)

  REAL(dp)    :: pttn(kbdim,klev,naermod), prhop(kbdim,klev,nmod), &
             pm6dry(kbdim,klev,nsol),  pm6rp(kbdim,klev,nmod), &
             pww(kbdim,klev,nmod),     paernl(kbdim,klev,nmod)

  !--- Local variables:

  INTEGER(i4) :: jk, jl, jmod

  REAL(dp)    :: zaerelse,                                         &
             zwso4,       zvso4,       zmso4,                  &
             zvso4h2o,    zmso4h2o,    zdso4h2o,               &
             zapp1,       ztk,         zrh

  REAL(dp)    :: ztk2,        zln3,        zln32,                  &
             zlnm,        zss2,        zlnm2

!CDIR unroll=5
  DO 100 jmod=1,nsol
     DO 90 jk=1,klev
        DO 80 jl=1,kproma

           !--- 1) Determine mass of non sulfate compounds in a particle: ----

           SELECT CASE (jmod)

              CASE (1)
                 zaerelse=0._dp 
              CASE (2)
                 zaerelse = pttn(jl,jk,ibcks)+pttn(jl,jk,iocks)
              CASE (3)
                 zaerelse = pttn(jl,jk,issas)+pttn(jl,jk,ibcas)+  &
                            pttn(jl,jk,iocas)+pttn(jl,jk,iduas)
              CASE (4)
                 zaerelse = pttn(jl,jk,isscs)+pttn(jl,jk,ibccs)+  &
                            pttn(jl,jk,ioccs)+pttn(jl,jk,iducs)
           END SELECT

           !--- 2) Calculation of the particle properties in the absense of ---
           !       other compounds than sulfate:

           IF (pttn(jl,jk,jmod) > 0.0_dp .AND. zaerelse < cmin_aerml) THEN

              !
              !--- 2.1) Calculation of the ambient particle properties: -------
              !      
              !--- Constrain ambient temperature to conditions for which the 
              !    parametrisation of the liquid water content works:
              
              ! Temperature:
              ztk = ptp1(jl,jk)
              ztk = MAX(ztk , 240._dp)

              ! Relative Humidity:
              zrh = prelhum(jl,jk)
              zrh = MAX(zrh , 0.05_dp)
              zrh = MIN(zrh , 0.90_dp)

              !--- Assign auxiliary variables:

              zapp1=papp1(jl,jk)
              zlnm = LOG(pttn(jl,jk,jmod))
              zlnm2 = zlnm*zlnm
              zss2 = zrh**2
              ztk2 = ztk*ztk
              zln3 = zlnm/3.0
              zln32 = zln3*zln3
              !
              !--- Percentage by weight of sulfate in the particle [%]:
              !    (Here we ignore any insoluble mass.)
              !
              zwso4 =wvb(1) + wvb(2)*zlnm + wvb(3)*zrh*zlnm + wvb(4)*ztk*zlnm +&
                     wvb(5)*zrh/ztk + wvb(6)*zlnm2*zrh + wvb(7)*zlnm2*ztk +    &
                     wvb(8)*zlnm*zss2 + wvb(9)*zlnm*ztk2 + wvb(10)*zlnm2*zss2 +&
                     wvb(11)*zlnm2*ztk2 + wvb(12)*zss2/ztk2 + wvb(13)*zlnm2 +  &
                     wvb(14)*zlnm2*zlnm + wvb(15)*zlnm2*zlnm2 +                &
                     wvb(16)*zss2*zrh/(ztk2*ztk) + wvb(17)*LOG(zrh*ztk/zapp1)

              !--- Dry mass of sulfate in an average particle [g]:

              zmso4 = pttn(jl,jk,jmod)*wh2so4/avo

              !--- Dry volume of sulfate in an average particle [cm3]:
              !    Any temperature or pressure dependency of the 
              !    sulfate density is ingored.

              zvso4 = zmso4/dh2so4

              !--- Mass of sulfate + water in an average particle [g]:

              zmso4h2o = zmso4/(zwso4/100.0)

              !--- Density of the sulfate-water fraction of an average 
              !particle [g cm-3]:
              !@@@ Check: changed zwvso4 into zwso4 (now the mass!)

              zdso4h2o = gmb(1) + gmb(2)*zwso4 + gmb(3)*zln3 + gmb(4)*zrh +   &
                         gmb(5)*ztk + gmb(6)*zln32 + gmb(7)*zln3/zrh +        &
                         gmb(8)*zln3/ztk + gmb(9)*ztk2
         
              !--- Limits for zdso4h2o: H2O(0.99) and pure H2SO4 (1.841):
              !
              zdso4h2o=MAX(zdso4h2o,dh2o)
              zdso4h2o=MIN(zdso4h2o,dh2so4)
              
              !-- Volume of sulfate-water fraction of an average particle [cm3]:

              zvso4h2o = zmso4h2o/zdso4h2o

              !--- 2.2) Calculatiion of the particle radii: --------------------

              !--- 2.2.1) Dry count mean radius [cm]:
 
              pm6dry(jl,jk,jmod)=((zvso4)*0.75_dp/pi)**(1._dp/3.)*ram2cmr(jmod)

              !--- 2.2.2) Equilibrium wet count mean radius [cm]:
             
              pm6rp(jl,jk,jmod) =&
                   ((zvso4h2o)*0.75_dp/pi)**(1._dp/3.)*ram2cmr(jmod)
              
              !--- 2.3) Assignment of the particle density [g cm-3]: -----------
 
              prhop(jl,jk,jmod)=zdso4h2o

              !--- 2.4) Store aerosol water for each mode [kg(water) m-3(air)]:

              pww(jl,jk,jmod)=(zmso4h2o-zmso4)*paernl(jl,jk,jmod)*1.E3

           END IF
80      END DO
90   END DO
100 END DO

END SUBROUTINE m7_equiz
!------------------------------------------------------------------------------!
SUBROUTINE m7_equimix(kproma,  kbdim, klev,   &
                      papp1,   pttn,  ptp1,   &
                      prelhum, pm6rp, pm6dry, &
                      prhop,   pww,   paernl  )
  !
  !    *m7_equimix*   calculates the ambient radii of the particles with
  !                   sulphate, b/o carbon and dust.
  !
  !    Authors:
  !    --------
  !    J. Wilson and E. Vignati, JRC (original source)               05/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)   2001
  !
  !    Purpose:
  !    --------
  !    This routine calculates the ambient radii for mixed particles without
  !    sea salt with mass of ttn molecules, converts them to count mean radii 
  !    and stores them in the array with address pm6rp.
  !    It additionally calculates the ambient particle density.
  !
  !    Method:
  !    -------
  !    The calculations of the ambient particle properties are based on 
  !    parameterisations of the mass of sulfate and density derived 
  !    by Julian Wilson from a regression analysis of results of solving 
  !    the generalised Kelvin equation using (F. J. Zeleznik, J. Phys. Chem.
  !    Ref. Data 20, 1157, 1991), for an H2SO4-H2O mixture, in the 
  !    following parameter ranges: 

  !       1e2   < pttn    < 1E11  [molecules]
  !       0.2   < prelhum < 0.9   [1]
  !       240   < ptp1    < 330   [K]
  !       10000 < papp1   < 100000  [Pa]
  !
  !    Due to the limitations of the parametrisation, the ambient temperature
  !    is restricted to a minimum of 240 K within this subroutine. 
  !
  !    For this application to mixed aerosols with an insoluble core we
  !    assume the H2O uptake by the particle to be that of a pure 
  !    H2SO4 / H2O particle containing the H2SO4 mass of the mixed aerosol.
  !
  !    Interface:
  !    ----------
  !    *m7_equimix* is called from *m7*
  !
  !    Externals:
  !    ----------
  !    none
  !
  ! pttn      = average mass for single compound in each mode 
  !             [in molec. for sulphate and in ug for bc, oc, ss, and dust]
  ! pm6rp     = count mean radius under ambient conditions
  ! pm6dry    = count mean radius under dry conditions
  ! paernl    = aerosol number for each mode [cm-3]
  ! pww       = aerosol water content for each mode [kg(water) m-3(air)]
  ! zwso4     = percentage by mass of sulfate in a H2O-H2SO4 particle 
  !             containing pttn molecules of sulfate under ambient conditions
  ! zvso4     = volume of pttn molecules of sulfate [cm3]
  ! zmso4     = mass of pttn molecules of sulfate [g]
  ! zdso4h2o  = density of sulfate-h2o fraction of a particle with average 
  !             mass [g.cm-3]
  ! zmso4h2o  = mass of sulfate-h2o fraction of a particle with average mass [g]
  ! zvso4h2o  = volume of sulfate-h2o fraction of a particle with average 
  !             mass [cm3]
  ! zinsvol   = total volume of insoluble compounds in a single particle of 
  !             average mass [cm3]
  ! zinsmass  = total mass of insoluble compounds in a single 
  !             particle of average mass [cm3]

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: papp1(kbdim,klev),        ptp1(kbdim,klev),       &
             prelhum(kbdim,klev)

  REAL(dp)    :: pttn(kbdim,klev,naermod), prhop(kbdim,klev,nmod), &
             pm6dry(kbdim,klev,nsol),  pm6rp(kbdim,klev,nmod), &
             pww(kbdim,klev,nmod),     paernl(kbdim,klev,nmod)
  !
  ! Local variables:
  !
  INTEGER(i4) :: jk, jl, jmod
  
  REAL(dp)    :: zseasalt,    zinsvol,     zinsmas,                &
             zwso4,       zvso4,       zmso4,                  &
             zvso4h2o,    zmso4h2o,    zdso4h2o,               &
             zapp1,       zrh,         ztk
                    
  REAL(dp)    :: zlnm2,       zln3,        zln32,      ztk2,       &
             zlnm,        zss2       
  !
  !
!CDIR unroll=5
  DO 100 jmod=2,nsol
     DO 90 jk=1,klev
        DO 80 jl=1,kproma

           !--- 1) Split particle quantities into soluble (sea salt) and -------
           !   non soluble (organic carbon + black carbon + dust) parts:
           !  (N.B. densities are assumed independent of temperature & pressure)

           SELECT CASE (jmod)

              CASE (2)
                 zseasalt=0._dp
                 zinsvol=1.e-6*(pttn(jl,jk,ibcks)/dbc+pttn(jl,jk,iocks)/doc)
                 zinsmas=1.e-6*(pttn(jl,jk,ibcks)+pttn(jl,jk,iocks))
              CASE (3)
                 zseasalt=pttn(jl,jk,issas)
                 zinsvol=1.e-6*(pttn(jl,jk,ibcas)/dbc+pttn(jl,jk,iocas)/doc+ &
                                pttn(jl,jk,iduas)/ddust)
                 zinsmas=1.e-6*(pttn(jl,jk,ibcas)+pttn(jl,jk,iocas)+         &
                                pttn(jl,jk,iduas))
              CASE (4)
                 zseasalt=pttn(jl,jk,isscs)
                 zinsvol=1.e-6*(pttn(jl,jk,ibccs)/dbc+pttn(jl,jk,ioccs)/doc+ &
                                pttn(jl,jk,iducs)/ddust)
                 zinsmas=1.e-6*(pttn(jl,jk,ibccs)+pttn(jl,jk,ioccs)+         &
                                pttn(jl,jk,iducs))
           END SELECT

   !--- 2) Calculation of the particle properties in the absense of sea salt: --
   !  

           IF (pttn(jl,jk,jmod) > 0.0_dp .AND. zinsvol > 0.0_dp .AND. &
                zseasalt <  cmin_aerml ) THEN

              !--- 2.1) Calculation of the ambient particle properties: --------
              !         
              !--- Constrain ambient temperature and relative humidity to 
              !    conditions for which the parametrisation of the liquid 
              !    water content works:
              !
              ! Temperature:
              ztk = ptp1(jl,jk)
              ztk = MAX(ztk , 240._dp)

              ! Relative Humidity:
              zrh = prelhum(jl,jk)
              zrh = MAX(zrh , 0.05_dp)
              zrh = MIN(zrh , 0.90_dp)
              
              !--- Assign auxiliary variables:

              zapp1=papp1(jl,jk)
              zlnm = LOG(pttn(jl,jk,jmod))
              zlnm2 = zlnm*zlnm
              zss2 = zrh**2
              ztk2 = ztk*ztk
              zln3 = zlnm/3.0
              zln32 = zln3*zln3

              !--- Percentage by weight of sulfate in the particle [%]:
              !    (Here we ignore any insoluble mass.)

              zwso4 =wvb(1) + wvb(2)*zlnm + wvb(3)*zrh*zlnm + wvb(4)*ztk*zlnm +&
                     wvb(5)*zrh/ztk + wvb(6)*zlnm2*zrh + wvb(7)*zlnm2*ztk +    &
                     wvb(8)*zlnm*zss2 + wvb(9)*zlnm*ztk2 + wvb(10)*zlnm2*zss2 +&
                     wvb(11)*zlnm2*ztk2 + wvb(12)*zss2/ztk2 + wvb(13)*zlnm2 +  &
                     wvb(14)*zlnm2*zlnm + wvb(15)*zlnm2*zlnm2 +                &
                     wvb(16)*zss2*zrh/(ztk2*ztk) + wvb(17)*LOG(zrh*ztk/zapp1)

              !--- Dry mass of sulfate in an average particle [g]:

              zmso4 = pttn(jl,jk,jmod)*wh2so4/avo

              !--- Dry volume of sulfate in an average particle [cm3]:
              !    Any temperature or pressure dependency of the 
              !    sulfate density is ingored.

              zvso4 = zmso4/dh2so4

              !--- Mass of sulfate + water in an average particle [g]:

              zmso4h2o = zmso4/(zwso4/100.0)

              !--- Density of the sulfate-water fraction of an average 
              !particle [g cm-3]:
              !@@@ Check: changed zwvso4 into zwso4 (now the mass!)

              zdso4h2o = gmb(1) + gmb(2)*zwso4 + gmb(3)*zln3 + gmb(4)*zrh +  &
                         gmb(5)*ztk + gmb(6)*zln32 + gmb(7)*zln3/zrh +       &
                         gmb(8)*zln3/ztk + gmb(9)*ztk2
         
              !--- Limits for zdso4h2o: H2O(0.99) and pure H2SO4 (1.841):
              !
              zdso4h2o=MAX(zdso4h2o,dh2o)
              zdso4h2o=MIN(zdso4h2o,dh2so4)
              
              !---Volume of sulfate-water fraction of an average particle [cm3]:

              zvso4h2o = zmso4h2o/zdso4h2o

              !--- 2.2) Calculatiion of the particle radii: --------------------

              !--- 2.2.1) Dry count mean radius [cm]:
 
              pm6dry(jl,jk,jmod)=((zvso4+zinsvol)*0.75_dp/pi)**(1._dp/3._dp)*&
                                  ram2cmr(jmod)

              !--- 2.2.2) Equilibrium wet count mean radius [cm]:
             
              pm6rp(jl,jk,jmod) =((zvso4h2o+zinsvol)*0.75_dp/&
                   pi)**(1._dp/3._dp)* ram2cmr(jmod)
              
              !--- 2.3) Calculation of the particle density [g cm-3]:----------
 
              prhop(jl,jk,jmod)=(zmso4h2o+zinsmas)/                  &
                                (zvso4h2o+zinsvol)                 

              !--- 2.4) Store aerosol water for each mode [kg(water) m-3(air)]:

              pww(jl,jk,jmod)=(zmso4h2o-zmso4)*paernl(jl,jk,jmod)*1.e3

           END IF

80      END DO
90   END DO
100 END DO
  
END SUBROUTINE m7_equimix
!------------------------------------------------------------------------------!
SUBROUTINE m7_equil (kproma, kbdim,  klev,   prelhum, paerml, paernl, &
                     pm6rp,  pm6dry, phplus, pww,     prhop,  ptp1    )
  !
  !**** *m7_equil* calculates the ambient radii of accumulation 
  !                and coarse mode particles in presence of sea salt
  !     
  !
  !  Authors:
  !  --------
  !  E. Vignati, JRC/EI (original source)                    01/2000
  !  P. Stier, MPI      (f90-version, changes, comments)        2001
  !
  !  Purpose:
  !  --------
  !  This routine calculates the equilibrium radius of sea salt
  !  particles, and of sea salt particles mixed with sulphate for 
  !  accumulation and coarse modes. 
  !
  !  Method:
  !  -------
  !  from the mass of sea salt and sulphate (in ug/m3),
  !   
  !
  !**Interface:
  !  ----------
  !  *m7_equil* is called from *m7*
  !
  !  Externals:
  !  ----------
  !  none
  !
  !  References:
  !  -----------
  !  Jacobson, M.Z., Tabazadeh, A., Turco, R.P., (1996). Simulating
  !     equilibrium within aerosols and nonequilibrium between gases
  !     and aerosols. 
  !  Tang, I.N., (1997). Thermodynamic and optical properties of 
  !     mixed-salt aerosols of atmospheric importance. 
  !
  !@@@  ToDo:  Rewrite with identifiers of ions and ion pairs!!!!!
  !
  !--- Local variables: to be completed
  !
  !--- Parameter list:
  !
  ! phplus = molality of H+ [mol/kg(water)]
  !
  !
  !--- Local variables:
  ! 
  !  zna(:,:,i) = sodium concentration [ug m-3]
  !               i=1: accumulation mode; i=2: coarse mode
  !  zcl(:,:,i) = chlorine concentration [ug m-3]
  !               i=1: accumulation mode; i=2: coarse mode
  !
  ! !@@@ To be completed!
  !
  ! Indices for the arrays zmm and zmmr:
  !  -----------------------------------
  ! |     ions               ion pair   |
  ! |  (mole m-3)            (mole m-3) |
  ! | i   zmm(i)             zmmr(i)    |
  ! | 1   Na+                NaCl       |
  ! | 2   Cl-                NaHSO4     |
  ! | 3   SO4--              Na2SO4     |
  ! | 4   HSO4-              H2-SO4     |
  ! | 5   H+                            |
  !  -----------------------------------

  !--- Parameters:

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: prelhum(kbdim,klev),        ptp1(kbdim,klev)
  
  REAL(dp)    :: paerml(kbdim,klev,naermod), paernl(kbdim,klev,nmod), &
             pm6rp(kbdim,klev,nmod),     pm6dry(kbdim,klev,nsol), &
             prhop(kbdim,klev,nmod),     pww(kbdim,klev,nmod)

  REAL(dp)    :: phplus(kbdim,klev,nss)
             
  !--- Local variables:

  INTEGER(i4) :: i,           jl,            jk,           jmod

  REAL(dp)    :: zaw,         zvolw,         zdryvol,      zdrymass,  &
             zkw,         zdryvol_mean,  zambvol_mean

  REAL(dp)    :: zmm(5),      zmmr(5),       zmol(5),      zmo(4),    &
             zmmt(4)
  !   
  !
  DO 50 jmod=1,nss
     DO 60 jk=1,klev
        DO 70 jl=1,kproma

           IF ((paerml(jl,jk,jmod+12) > cmin_aerml) .AND. &
                (paernl(jl,jk,jmod+2)>cmin_aernl))  THEN
              !
              !- 1) Dry Calculations: ---------------------------------------
              !
              !- 1.1) Calculate initial concentrations of the compounds 
              !         in mole/m+3:
              !
              !--- Na, Cl:
              !
              zmm(1)=paerml(jl,jk,jmod+12)*1.E-6 / wnacl   ! n(Na)    [mole m-3]
              zmm(2)=paerml(jl,jk,jmod+12)*1.E-6 / wnacl   ! n(Cl)    [mole m-3]
              !      [         g m-3           ] / [g mole-1]       = [mole m-3]
              !
              !--- SO4--:
              !
              zmm(3)=paerml(jl,jk,jmod+2) *1.E+6 / avo     ! n(H2SO4) [mole m-3]
              !      [      m-3                ] / [mole-1]         = [mole m-3]
              !
              !--- HSO4-:
              !
              zmm(4)=0._dp                                ! n(HSO4)  [mole m-3]
              !
              !- 1.2) Calculation of the concentration of the different species:
              !         The ions are supposed to be arranged such that
              !         sodium is associated to sulphate first in
              !         Na2SO4, the remaining sodium is associated to Cl
              !         in form of NaCl.
              !
              zmmt(1)=zmm(1)                ! n(Na)  
              zmmt(3)=zmm(3)                ! n(SO4)    
              !
              zmmr(3)=MIN(zmmt(1)/2._dp , zmmt(3))     ! n(Na2SO4) 
              zmmt(1)=zmmt(1)-2._dp*zmmr(3) ! Remaining n(Na) after association
                                         ! with Na2SO4: n(Na)=n(Na)-2*n(Na2SO4)
              !
              zmmr(1)=MIN(zmm(2),zmmt(1))    ! n(NaCl) 
              !
              zmm(2)=zmmr(1)                 ! n(Cl) bound in NaCl, the rest is
                                             ! assumed to evaporate in presence
                                             !  of SO4
              !
              zmmr(2)=0._dp                  ! n(NaHSO4)
              !
              zmmr(4)=zmm(3)-zmmr(2)-zmmr(3) ! n(H2-SO4)(t)=n(H2SO4)(t0)-
                                             ! n(NaHSO4)-n(Na2SO4)
              !                              ! as n(H2SO4)(t0)=n(SO4--)(t0)
              !
              !--- 1.3) Total aerosol dry volume [cm3/m3]:
              !
              zdryvol= zmmr(1)*wnacl/dnacl               + &
                       zmmr(2)*wnahso4/dnahso4           + &
                       zmmr(3)*wna2so4/dna2so4           + &
                       zmmr(4)*wh2so4/dh2so4             + &
                       paerml(jl,jk,jmod+14)*1.e-6/ddust + &
                       paerml(jl,jk,jmod+5) *1.e-6/dbc   + &
                       paerml(jl,jk,jmod+9) *1.e-6/doc
              !
              !--- 1.4) Mean aerosol dry volume [cm+3]:
              zdryvol_mean = zdryvol    / (paernl(jl,jk,jmod+2)*1.E6)
              ! [cm+3]     = [cm+3/m+3] /  [         m-3            ]
              !
              !--- 1.5) Dry radius [cm]:
              !
              pm6dry(jl,jk,jmod+2)=((3._dp/(4._dp*pi))&
                   *zdryvol_mean)**(1._dp/3._dp) * ram2cmr(jmod+2)
              !
              !--- 1.6) Total aerosol dry mass [gr/m3]:
              !
              zdrymass= zmmr(1)*wnacl                    + &
                        zmmr(2)*wnahso4                  + &
                        zmmr(3)*wna2so4                  + &
                        zmmr(4)*wh2so4                   + &
                        paerml(jl,jk,jmod+14)*1.e-6      + & ! Dust
                        paerml(jl,jk,jmod+5)*1.e-6       + & ! Black Carbon
                        paerml(jl,jk,jmod+9)*1.e-6           ! Organic Carbon
              !
              !
              !--- 2) Wet calculations: ----------------------------------------
              !
              !--- Set threshold for relative humidity:
              !    If RH is smaller than the Critical Relative Humidity 
              !    (currently crh=0.45) the equilibrium radius is set to the dry
              !    radius:

              IF (prelhum(jl,jk) < crh) THEN
                 !
                 phplus(jl,jk,jmod)  = 0._dp 
                 !
                 pww(jl,jk,jmod+2)     = 0._dp
                 !
                 pm6rp(jl,jk,jmod+2) = pm6dry(jl,jk,jmod+2)
                 !
                 prhop(jl,jk,jmod+2) = zdrymass/zdryvol
                 !
              ELSE
                 !
                 !- 2.1) Calculate thermodynamic properties under ambient 
                 !       conditions
                 !
                 !- 2.1.1) Water activity:
                 !
                 zaw=prelhum(jl,jk)
                 !
                 !- 2.1.2) Molality as function of the water activity:
                 !         Currently sulfate is assumed to be fully dissociated,
                 !         i.e. zmmr(2)=0. and zmo(2) is not calculated.
                 !
                 !           Changed reference to Jacobson et al. (1996):
                 !
                 !--- NaCl:
                 
                 zmo(1)=(-1.918004E2_dp+2.001540E3*zaw-8.557205E3*zaw**2       &
                        +1.987670E4*zaw**3-2.717192E4*zaw**4+2.187103E4*zaw**5 &
                         -9.591577E3*zaw**6+1.763672E3*zaw**7            )**2 

                 !--- NaHSO4:

                 zmo(2)=(+4.662777E0_dp-1.128472E1_dp*zaw+7.049464E1*zaw**2    &
                        -2.788050E2*zaw**3+6.103105E2*zaw**4-7.409417E2*zaw**5 &
                         +4.614577E2*zaw**6-1.150735E2*zaw**7       )**2 

                 !--- Na2SO4:

                 zmo(3)=(-3.295311E3_dp+3.188349E4*zaw-1.305168E5*zaw**2       &
                        +2.935608E5*zaw**3-3.920423E5*zaw**4+3.109519E5*zaw**5 &
                         -1.356439E5*zaw**6+2.510249E4*zaw**7              )**2

                 !--- H2-SO4:
                 
                 zmo(4)=(+5.611895_dp-1.387446E1*zaw+1.750682E1*zaw**2        &
                        +7.138146E1*zaw**3-3.109173E2*zaw**4+4.662288E2*zaw**5 &
                         -3.128612E2*zaw**6+7.76097E1*zaw**7              )**2
                 !
                 !
                 !- 2.2) Calculation of the water content in kg water/m3 air:
                 !       (zmmr[mole/m3(air)]/zmo[mole/kg(water)] =
                 !       zww[kg(water)/m3(air)]
                 !
                 pww(jl,jk,jmod+2)=zmmr(1)/zmo(1)+zmmr(2)/zmo(2)+&
                                            zmmr(3)/zmo(3)+zmmr(4)/zmo(4)
                 !
                 !--- 2.3) Calculate the molality of the ions
                 !
                 !--- 2.3.1) For Na+, Cl-, SO4--, HSO4- :
                 !
                 DO i=1,4
                    zmol(i)=zmm(i)/pww(jl,jk,jmod+2)
                 END DO
                 !
                 !--- 2.3.2) For h+ :
                 !
                 !    [H+] = -[Na+] + [Cl-] + [OH-] + 2[SO4--] +[HSO4-]
                 !
                 !    with [OH-] = kw/[H+]
                 !
                 !--- Calculate autodissociation constant (kw) for water:
                 !    (Seinfeld &Pandis (1998): Eq. (6.5) + Table 6.5)

                 zkw=1.0e-14*exp( (13.35_dp/r_kcal) &
                      * (1._dp/298._dp - 1._dp/ptp1(jl,jk)) )

                 !--- Calculate molality of H+:

                 zmol(5)=( (-zmol(1)+zmol(2)+2._dp*zmol(3)+zmol(4)) +    &
                 SQRT( (zmol(1)-zmol(2)-2._dp*zmol(3)-zmol(4))**2 +      &
                 4._dp*zkw ) ) / 2._dp
                 !
                 !
                 zmm(5)=zmol(5)*pww(jl,jk,jmod+2)
                 !
                 phplus(jl,jk,jmod)=zmol(5)
                 !
                 !
                 !--- 2.4) Calculation of the wet radius
                 !
                 !--- 2.4.1) Total water volume  [cm3/m3]:
                 !
                 zvolw=pww(jl,jk,jmod+2)/(dh2o*1.e-3)  
                                       ![cm3/m3]=[kg/m3]/([g/cm3]*[1.E-3 kg/g])
                 !
                 !
                 !--- 2.4.2) Mean aerosol ambient volume:
                 zambvol_mean = (zdryvol+zvolw) / (paernl(jl,jk,jmod+2)*1.e6)
                 ! [cm+3]     = [  cm+3/m+3   ] / [         m-3             ]  
                 
                 !--- 2.4.3) Equilibrium wet count mean radius [cm]:
                 !
                 pm6rp(jl,jk,jmod+2)=((3._dp/&
                      (4._dp*pi))*zambvol_mean)**(1._dp/3._dp) &
                                      * ram2cmr(jmod+2)
                 !
                 !--- 2.4.4) Calculation of the particle density (g cm-3):
                 !
                 prhop(jl,jk,jmod+2)=(zdrymass+zvolw*dh2o)/(zvolw+zdryvol)
                 !

              END IF !(prelhum(jl,jk) < crh)

           END IF !((paerml(jl,jk,jmod+12)> 1.e-15)
                  !.AND.(paernl(jl,jk,jmod+2)> 1.e-10)) 

70      END DO
60   END DO
50 END DO
  !
END SUBROUTINE m7_equil
!------------------------------------------------------------------------------!
SUBROUTINE m7_dgas(kproma, kbdim, klev,  pso4g,  paerml, paernl, &
                   ptp1,   papp1, pm6rp, pso4_5, pso4_6, pso4_7, &
                   ztmst                                         ) 
  !
  !**** *m7_dgas*  calculates the transfer of mass due to 
  !                sulfate condensation on soluble and 
  !                insoluble modes.
  !
  !    Authors:
  !    -----------
  !    J. Wilson, E. Vignati, JRC/EI (original source)                05/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)    2001   
  !
  !    Purpose:
  !    -----------
  !    This routine calculates the changes in aerosol mass and
  !    gas phase sulfate due to sulfate condensation.
  !
  !**  Interface:
  !    -----------
  !    *m7_dgas* is called from *m7*
  !
  !    Method:
  !    -----------------
  !    The transfer of sulfate to the particles is based on
  !    Fuchs (1959). Soluble and insoluble particles are distinguished
  !    by a different accomodation coefficient "caccso4" defined above.
  !    (Currently 1.0 for soluble and 0.3 for insoluble modes). 
  !
  !    Externals
  !    -----------
  !    none
  !
  !    References:
  !    -----------
  !    Fuchs, N.A. (1959). Evaporation and droplet growth in gaseous media;
  !       Pergamon, New York, pp72.
  !
  !--- Parameter list:
  !
  ! pso4g     = mass of gas phase sulfate [molec. cm-3]
  ! pm6rp     = mean mode actual radius (wet radius for soluble modes 
  !             and dry radius for insoluble modes) [cm]
  ! pso4_x    = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  ! 
  !--- Local Variables:
  !
  ! zde       = molecular diffusion []
  ! zvelb     = velocity []
  ! zcondo    = condensation coefficient []
  ! zc2(nmod) = flux of sulfate condensing on the respective mode 
  !             per sulfate gas phase concentration []
  ! zcondo    = total flux of condensing sulfate 
  !             per sulfate gas phase concentration []
  ! zfcond    = total mass of condensing sulfate for one timestep []

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: ptp1(kbdim,klev),        papp1(kbdim,klev),                 &
                 pso4g(kbdim,klev),       pso4_5(kbdim,klev),                &
                 pso4_6(kbdim,klev),      pso4_7(kbdim,klev)
  !
  REAL(dp)    :: paernl(kbdim,klev,nmod)
  REAL(dp)    :: paerml(kbdim,klev,naermod)
  REAL(dp)    :: pm6rp(kbdim,klev,nmod)
  !
  ! Local variables:

  INTEGER(i4) :: jl,         jk,         jmod
  
  REAL(dp)    :: zfcond,     zftot,      zpbyone,    zde2,                   & 
             zvelb,      zxibc,      zm6rp,      zf1,                    &
             ztmst

  REAL(dp)    :: zcondo(kbdim,klev)

  REAL(dp)    :: zc2(kbdim,klev,nmod)

  zcondo(:,:)=0.0_dp
      
  zc2(:,:,:) = 0.0_dp  
  
  !--- 1) Calculate condensation rate for cm diameter sulphate aerosols: ---
  !
  DO jmod=1,nmod
     DO jk=1,klev
        DO jl=1,kproma
           IF (pm6rp(jl,jk,jmod).GT.1.e-15_dp) THEN

              !--- Diffusion coefficient (Reference???):

              zpbyone=1000.0 / (papp1(jl,jk)/100.0)

              zde2=0.073 * zpbyone * (ptp1(jl,jk) / 298.15)**1.5  

              !--- Mean molecule velocity (Moore, 1962 (S+P equ. 8.2)):

              zvelb=SQRT(8.0 * rerg * ptp1(jl,jk) / pi / wh2so4)  

              !--- ???Fuchs???

              zxibc=8.0 * zde2 / pi / zvelb
              !
              ! Use count median radius:

              zm6rp=pm6rp(jl,jk,jmod)

              !-- Distance from particle up to which the kinetic regime applies:

              zf1=( (zm6rp + zxibc)**3.0 - (zm6rp**2.0 + zxibc**2.0)**1.5 ) / &
                  (3.0 * zm6rp * zxibc) - zm6rp

              !--- Diffusive flux to single particle surface:
              !    (Elisabetta's thesis: fraction in equ. 2.26)

              zc2(jl,jk,jmod)=(4.0 * pi * zde2 * zm6rp ) /                    &
                              ((4.0 * zde2) / (zvelb * zm6rp * caccso4(jmod)) +&
                               (zm6rp/(zm6rp+zf1))                            )

              !--- Total diffusive flux to all particles in the respective mode:
              !    (per concentration of gas phase sulfate)

              zc2(jl,jk,jmod)=zc2(jl,jk,jmod) * paernl(jl,jk,jmod)

              !--- Total diffusive flux to all particles in all modes:
              !    (per concentration of gas phase sulfate)

              zcondo(jl,jk)=zcondo(jl,jk)+ zc2(jl,jk,jmod)  

           END IF
        END DO
     END DO
  END DO
  !
  !--- 2) Calculation of the new sulfate aerosol masses and of the ---------
  !       mass of sulfate condensing on the respective modes:
  !
  DO jk=1,klev
     DO jl=1,kproma
        IF(zcondo(jl,jk).GT.0._dp.AND.pso4g(jl,jk).GT.1.e-10_dp) THEN

           !--- Total diffusive flux to all particles in all modes:

           zfcond=zcondo(jl,jk)*pso4g(jl,jk)

           !--- Total mass of sulfate condensing during 1 timestep:

           zftot=zfcond*ztmst
           
           !--- Limit condensing sulfate to 
           !    fmax[%] x (available gas-phase sulfate) :

           zfcond=MIN(zftot,(pso4g(jl,jk)*fmax)) 

           !--- Remaining gas phase sulfate:
           
           pso4g(jl,jk)=pso4g(jl,jk)-zfcond
           
           !--- Add mass to sulfate compounds:
           !    zc2(:,:,jmod)/zcondo = fraction of total sulfate flux
           !                           that condenses on mode jmod
           !    => (   "    )*zfcond = mass of sulfate condensing on
           !                           the respective mode

           paerml(jl,jk,iso4ns)=paerml(jl,jk,iso4ns)+ &
                                zc2(jl,jk,iso4ns)/zcondo(jl,jk)*zfcond
           paerml(jl,jk,iso4ks)=paerml(jl,jk,iso4ks)+ &
                                zc2(jl,jk,iso4ks)/zcondo(jl,jk)*zfcond
           paerml(jl,jk,iso4as)=paerml(jl,jk,iso4as)+ &
                                zc2(jl,jk,iso4as)/zcondo(jl,jk)*zfcond
           paerml(jl,jk,iso4cs)=paerml(jl,jk,iso4cs)+ &
                                zc2(jl,jk,iso4cs)/zcondo(jl,jk)*zfcond

           !--- Mass of sulphate condensing on the insoluble modes:
           !    (Transfer from insoluble to soluble modes 
           !     calculated in m7_concoag.)
           
           pso4_5(jl,jk)=zc2(jl,jk,5)/zcondo(jl,jk)*zfcond
           pso4_6(jl,jk)=zc2(jl,jk,6)/zcondo(jl,jk)*zfcond
           pso4_7(jl,jk)=zc2(jl,jk,7)/zcondo(jl,jk)*zfcond

        ELSE 

           pso4_5(jl,jk)=0._dp
           pso4_6(jl,jk)=0._dp
           pso4_7(jl,jk)=0._dp

        END IF
     END DO
  END DO

END SUBROUTINE m7_dgas
!------------------------------------------------------------------------------!
SUBROUTINE m7_dnum(kproma, kbdim,   klev,          &
                   pso4g,  paerml,  paernl, ptp1,  &
                   papp1,  prelhum, pm6rp,  prhop, &
                   pso4_5, pso4_6,  pso4_7, ztmst) 
  !
  !  *m7_dnum*  changes gas phase sulfate, aerosol numbers and masses
  !             due to nucleation and coagulation
  !
  !  Authors:
  !  ---------
  !  J. Wilson  and E. Vignati, JRC/EI (original source)                09/2000
  !  P. Stier, MPI                     (f90-version, changes, comments)    2001 
  !
  !  Version: 
  !  --------- 
  !  This version is equivalent to the version dnum2 of the m7 boxmodel. 
  !
  !  Purpose
  !  ---------
  !  This routine calculates new gas phase sulfate and aerosol
  !  numbers and masses after timestep ztmst.
  !
  !  Interface
  !  -----------
  !  *m7_dnum* is called from *m7*
  !
  !  Externals
  !  -----------
  !  *m7_coaset*   calculates the coagulation kernels
  !  *m7_nuck*     calculates the integral mass nucleated sulfate and
  !                the number of nucleated particles over one timestep
  !  *m7_delcoa*   integrates equations for the changes in aerosol numbers
  !                dn/dt=c -an2 -bn over one timestep and calculates the 
  !                corresponding changes in aerosol masses
  !  *m7_concoag*  calculates particle numbers and mass moved from the 
  !                insoluble to the mixed modes due to coagulation with
  !                smaller mixed particles and condensation of H2SO4.
  !
  !  Warning:
  !  --------
  !  For optimization purposes currently only "physically reasonable" elements 
  !  of the coagulation kernel zcom are calculated in m7_concoag.
  !  These elements are specified in the matrix locoagmask in mo_aero_m7.
  !  Check carefully and adapt locoagmask accordingly before changes in
  !  the code below.
  !
  !--- Parameters:
  !
  !  pso4g      = mass of gas phase sulfate [molec. cm-3]
  !  paerml     = total aerosol mass for each compound 
  !               [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and dust]
  !  paernl     = aerosol number for each mode [cm-3]
  !  ptp1       = atmospheric temperature at time t+1 [K]
  !  papp1      = atmospheric pressure at time t+1 [Pa]
  !  prelhum    = atmospheric relative humidity [%]
  !  pm6rp      = mean mode actual radius (wet radius for soluble modes 
  !               and dry radius for insoluble modes) [cm]
  !  prhop      = mean mode particle density [g cm-3]
  !  pso4_x     = mass of sulphate condensed on insoluble mode x []
  !
  !--- Local variables:
  !
  ! zcom            = general coagulation coefficient []
  !                   (calculated in m7_coaset)
  ! za              = effectively used coagulation coefficient for 
  !                   unimodal coagulation []
  ! zb              = effectively used coagulation coefficient for 
  !                   inter-modal coagulation []
  ! zbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x (finally calculated in m7_delcoa)
  !                   that is moved to mode y / y+1 (modes 5,6,7 and 
  !                   mode 1,2 resp.) [1]
  !                   !@@@ Careful! Inconsistent usage!!!
  ! za4delt(:,:,:)  = change in H2SO4 mass of the respective mode over 
  !                   one timstep 
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in m7_nuck)
  !                      - coagulation (calculated in m7_concoag)
  ! zanew           = number of nucleated particles
  !                   zanew=za4delt/critn i.e. mass of formed sulfate 
  !                   divided by an assumed mass of a nucleus. []
  !                   (calculated in m7_nuck)
  ! zxxavy          = average mass of species xx in mode y []
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  !                   [molecules for sulfate and ug for others]
  !

  !--- Parameters:

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: pso4g(kbdim,klev),          ptp1(kbdim,klev),             &
             papp1(kbdim,klev),          prelhum(kbdim,klev),              &
             pso4_5(kbdim,klev),         pso4_6(kbdim,klev),               &
             pso4_7(kbdim,klev) 

  REAL(dp)    :: paerml(kbdim,klev,naermod), paernl(kbdim,klev,nmod),      &
             pm6rp(kbdim,klev,nmod),     prhop(kbdim,klev,nmod)

             
  ! Local Variables:

  INTEGER(i4) :: jl, jk

  REAL(dp):: ztmst ! mz_ak_20041014 !qqq

  REAL(dp)    :: zanew(kbdim,klev)

  REAL(dp)    :: za(kbdim,klev,nmod),        zb(kbdim,klev,nmod),              &
             za4delt(kbdim,klev,naermod)

  REAL(dp)    :: zbfract1(kbdim,klev,nmod-1),zbfract2(kbdim,klev,nmod-1),      &
             zbfract5(kbdim,klev,3),     zbfract6(kbdim,klev,2),           &
             zbfract7(kbdim,klev,2)
             

  REAL(dp)    :: zcom(kbdim,klev,nmod,nmod)


  !--- 0) Initialisations: ----------------------------------------------

  za4delt(:,:,:) = 0._dp  ! Mode 1 changed by m7_nuck if lsnucl=TRUE .
                          ! Has to be initialized for the other modes.
  zanew(:,:)     = 0._dp  ! Changed by m7_nuck if lsnucl=TRUE .

  !--- 1) Calculate  coagulation coefficients: --------------------------
  !
!CDIR NOIEXPAND
  IF (lscoag) CALL m7_coaset(kproma, kbdim, klev,  paernl, ptp1, &
                             papp1,  pm6rp, prhop, zcom          )
  !
  !--- 2) Calculate nucleation rate, number of nucleated particles ------
  !       and changes in gas phase sulfate mass.
  !
!CDIR NOIEXPAND
  IF (lsnucl) CALL m7_nuck(kproma, kbdim,   klev,  pso4g,   &
                           ptp1,   prelhum, zanew, za4delt, &
                           ztmst                             )

  !
  !--- 3) Assign coagulation coefficients (unimodal coagulation)---------
  !       and the normalised fraction of particle numbers moved 
  !       between the modes (inter-modal coagulation):
  !
  !       The general equation for dn/dt for each mode is:
  !
  !       dn/dt=-za*n^2 - zb*n + zc
  !
  !       where za=unimodal coagulation coefficient (zcom(mod))
  !             zb=inter-modal coagulation with higher modes
  !                (zcom(mod) * n(jmod+1))
  !             zc=particle formation rate 
  !                (=zanew/ztmst if jmod=1, zero for higher modes) 
  !
  !             zb is zero when n(jmod+1)=zero, or jmod=naermod
  !
  IF (lscoag) THEN

     DO jk=1,klev
        DO jl=1,kproma 

           !---  3.1) Unimodal coagulation coefficients:
           !@@@Coag:
           za(jl,jk,1)=zcom(jl,jk,1,1)/2.    ! Unimodal coagulation
           za(jl,jk,2)=zcom(jl,jk,2,2)/2.    ! only allowed for modes
           za(jl,jk,3)=zcom(jl,jk,3,3)/2.    ! 1,2,3 and 5.
           za(jl,jk,4)=0._dp
           za(jl,jk,5)=zcom(jl,jk,5,5)/2.
           za(jl,jk,6)=0._dp
           za(jl,jk,7)=0._dp
           !
           !---  3.2) Inter-modal coagulation - soluble modes
           ! 
           !--- Sum all higher mode coagulation terms for 
           !    soluble modes 1,2,3,4:
           !
           !--- 3.2.1) Mode 1:
           
           !--- Number of particles (zbfract1(:,:,x)) that are moved 
           !    from mode 1 to the mode x+1 :
           !    !@@@ Clumsy! Change to x - also in concoag!!!
           zbfract1(jl,jk,1)=zcom(jl,jk,2,1)*paernl(jl,jk,2)
           zbfract1(jl,jk,2)=zcom(jl,jk,3,1)*paernl(jl,jk,3)
           zbfract1(jl,jk,3)=zcom(jl,jk,4,1)*paernl(jl,jk,4)
           zbfract1(jl,jk,4)=zcom(jl,jk,5,1)*paernl(jl,jk,5)
           zbfract1(jl,jk,5)=zcom(jl,jk,6,1)*paernl(jl,jk,6)
           zbfract1(jl,jk,6)=zcom(jl,jk,7,1)*paernl(jl,jk,7)
           !
           !--- Sum of all particles that are moved from mode 1:
           !
           zb(jl,jk,1)=zbfract1(jl,jk,1)+zbfract1(jl,jk,2)+            &
                       zbfract1(jl,jk,3)+zbfract1(jl,jk,4)+            &
                       zbfract1(jl,jk,5)+zbfract1(jl,jk,6)
           !
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 1:
           !
           IF (zb(jl,jk,1).GT.0.0_dp) THEN
              zbfract1(jl,jk,1)=zbfract1(jl,jk,1)/zb(jl,jk,1)
              zbfract1(jl,jk,2)=zbfract1(jl,jk,2)/zb(jl,jk,1)
              zbfract1(jl,jk,3)=zbfract1(jl,jk,3)/zb(jl,jk,1)
              zbfract1(jl,jk,4)=zbfract1(jl,jk,4)/zb(jl,jk,1)
              zbfract1(jl,jk,5)=zbfract1(jl,jk,5)/zb(jl,jk,1)
              zbfract1(jl,jk,6)=zbfract1(jl,jk,6)/zb(jl,jk,1)
           END IF
           !
           !--- 3.2.2) Mode 2:
           !
           !--- Number of particles (zbfract1(:,:,x)) that are moved 
           !    from mode 2 to the mode x+1 :
           !
           zbfract2(jl,jk,2)=zcom(jl,jk,3,2)*paernl(jl,jk,3)
           zbfract2(jl,jk,3)=zcom(jl,jk,4,2)*paernl(jl,jk,4)
           zbfract2(jl,jk,4)=0._dp
           zbfract2(jl,jk,5)=zcom(jl,jk,6,2)*paernl(jl,jk,6)
           zbfract2(jl,jk,6)=zcom(jl,jk,7,2)*paernl(jl,jk,7)
           
           !--- Sum of all particles that are moved from mode 2:
           
           zb(jl,jk,2)=zbfract2(jl,jk,2)+zbfract2(jl,jk,3)+            &
                       zbfract2(jl,jk,4)+zbfract2(jl,jk,5)+            &
                       zbfract2(jl,jk,6)

           !--- Normalize particle numbers by the total number of 
           !    particles moved from mode 2:

           IF (zb(jl,jk,2).GT.0.0_dp) THEN
              zbfract2(jl,jk,2)=zbfract2(jl,jk,2)/zb(jl,jk,2)
              zbfract2(jl,jk,3)=zbfract2(jl,jk,3)/zb(jl,jk,2)
              zbfract2(jl,jk,4)=zbfract2(jl,jk,4)/zb(jl,jk,2)
              zbfract2(jl,jk,5)=zbfract2(jl,jk,5)/zb(jl,jk,2)
              zbfract2(jl,jk,6)=zbfract2(jl,jk,6)/zb(jl,jk,2)
           END IF
        
           !--- 3.2.3) Mode 3 and Mode 4 - considered ineffective:

           zb(jl,jk,3)=0.0_dp

           zb(jl,jk,4)=0.0_dp

           !
           !--- 3.3) Inter-modal coagulation - insoluble modes
           !
           !         For the insoluble modes coagulation with soluble modes
           !         is a sink as they are transfered to the corresponding
           !         mixed/solublemode. Therefore, terms with a lower mode 
           !         number or the same mode number are included. 
           !         (!@@@ There are still some switches for testing.)
           !
           !--- 3.3.1) Mode 5:
           !
           !--- Number of particles (zbfract5(:,:,x)) that are moved 
           !    from mode 5 to the mode x:
           
           zbfract5(jl,jk,1)=0._dp
           zbfract5(jl,jk,2)=zcom(jl,jk,2,5)*paernl(jl,jk,2)
           zbfract5(jl,jk,3)=zcom(jl,jk,3,5)*paernl(jl,jk,3)
           
           !--- Sum of all particles that are moved from mode 5:
           
           zb(jl,jk,5)=zbfract5(jl,jk,1)+zbfract5(jl,jk,2)+zbfract5(jl,jk,3)
           
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 5:
           
           IF (zb(jl,jk,5).GT.0.0_dp) THEN
              zbfract5(jl,jk,1)=zbfract5(jl,jk,1)/zb(jl,jk,5)
              zbfract5(jl,jk,2)=zbfract5(jl,jk,2)/zb(jl,jk,5)
              zbfract5(jl,jk,3)=zbfract5(jl,jk,3)/zb(jl,jk,5)
           END IF
           
           !--- 3.3.2) Mode 6:
           
           !--- Number of particles (zbfract6(:,:,x)) that are moved 
           !    from mode 6 to the mode x:
           
           !@@@zbfract6(jl,jk,1)=zcom(jl,jk,1,6)*paernl(jl,jk,1)
           !@@@zbfract6(jl,jk,2)=zcom(jl,jk,2,6)*paernl(jl,jk,2)
           zbfract6(jl,jk,1)=0._dp
           zbfract6(jl,jk,2)=0._dp
           
           !--- Sum of all particles that are moved from mode 6:
           
           zb(jl,jk,6)=zbfract6(jl,jk,1)+zbfract6(jl,jk,2)
           
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 6:
           
           IF (zb(jl,jk,6).GT.0.0_dp) THEN
              zbfract6(jl,jk,1)=zbfract6(jl,jk,1)/zb(jl,jk,6)
              zbfract6(jl,jk,2)=zbfract6(jl,jk,2)/zb(jl,jk,6)
           END IF
           
           !--- 3.3.3) Mode 7:
           
           !--- Number of particles (zbfract7(:,:,x)) that are moved 
           !    from mode 7 to the mode x:
           
           !@@@ zbfract7(jl,jk,1)=zcom(jl,jk,1,7)*paernl(jl,jk,1)
           !@@@ zbfract7(jl,jk,2)=zcom(jl,jk,2,7)*paernl(jl,jk,2)
           zbfract7(jl,jk,1)=0._dp
           zbfract7(jl,jk,2)=0._dp  
           
           !--- Sum of all particles that are moved from mode 7:
           
           zb(jl,jk,7)=zbfract7(jl,jk,1)+zbfract7(jl,jk,2)
           
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 7:
           
           IF (zb(jl,jk,7).GT.0.0_dp) THEN
              zbfract7(jl,jk,1)=zbfract7(jl,jk,1)/zb(jl,jk,7)
              zbfract7(jl,jk,2)=zbfract7(jl,jk,2)/zb(jl,jk,7)
           END IF

        END DO
     END DO
     !
  ELSE

     za(:,:,:)       = 0._dp
     zb(:,:,:)       = 0._dp
     zbfract1(:,:,:) = 0._dp
     zbfract2(:,:,:) = 0._dp
     zbfract5(:,:,:) = 0._dp
     zbfract6(:,:,:) = 0._dp
     zbfract7(:,:,:) = 0._dp

  END IF !(lscoag)
  !
  !
!CDIR NOIEXPAND
  CALL m7_delcoa(kproma,   kbdim,  klev,     paerml,   &
                 paernl,   pm6rp,  za4delt,  zanew,    &
                 za,       zb,     zbfract1, zbfract2, &
                 zbfract5, pso4_5, pso4_6,   pso4_7  ,   &
                  ztmst                                   ) !qqq

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

CONTAINS

SUBROUTINE m7_coaset(kproma, kbdim, klev,  paernl, ptp1, &
                     papp1,  pm6rp, prhop, pcom          )
  !
  ! *m7_coaset*  calculates the coagulation kernels between the modes
  !
  ! Authors:
  ! ---------
  ! J. Wilson  and E. Vignati, JRC/EI (original source)                09/2000
  ! P. Stier, MPI                     (f90-version, changes, comments)    2001 
  !
  ! Modifications:
  ! --------------
  ! Philip Stier, MPI                             2001
  !
  ! Purpose
  ! ---------
  ! This routine calculates the coaglation kernels between particles
  ! with the count median radii of the three modes.
  ! Coagulation allowed between the following modes:
  ! soluble modes:   1+1=1, 2+2=2, 1+2=2, 1+3=3, 1+4=4, 2+3=3, 2+4=4
  ! insoluble modes: 2i+2i=2i
  ! mixed modes:     1+2i=2, 1+4i=4, 2+2i=2, 2+4i=4, 3+2i=3, 4+2i=4
  !
  ! Interface:
  ! -----------
  ! *m7_coaset* is called from *m7_dnum*
  !
  ! Externals:
  ! -----------
  ! none
  !
  ! Reference:
  ! -----------
  ! The calculations are based on:
  ! Fuchs, N.A. (1964). The Mechanics of Aerosols. Pergamon Press. Oxford. 
  ! (Chapter VII, Section 49)
  !
  ! Warning:
  ! --------
  ! For optimization purposes currently only "physically reasonable" elements 
  ! of the coagulation kernel pcom are calculated in m7_concoag. These elements
  ! are specified in the matrix locoagmask in messy_m7. Check carefully and 
  ! adapt locoagmask accordingly  before changes in the code below.
  !
  !--- Parameter list:
  !
  !  paernl            = aerosol number for each mode [cm-3]
  !  ptp1              = atmospheric temperature at time t+1 [K]
  !  papp1             = atmospheric pressure at time t+1 [Pa]
  !  pm6rp             = mean mode actual radius (wet radius for soluble modes 
  !                      and dry radius for insoluble modes) [cm]
  !  prhop             = mean mode particle density [g cm-3]
  !  pcom(:,:,jm1,jm2) = Coagulation coefficient for modes jm1 and jm2 []
  !
  !--- List of local variables:
  !
  ! zwlc              = mean free pathlength []
  ! zairvisc          = air viscosity []
  ! zrpav             = average radius of the two interacting modes
  ! zpvx              = volume of the xth interacting mode
  ! zpmx              = mass of the xth interacting mode
  !
  ! zrknudx           = knudsen number of the xth interacting mode
  ! zpd2x             = particle diffusion of the xth interacting mode
  !
  !--- Parameters:
  !
    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: ptp1(kbdim,klev),          papp1(kbdim,klev)

  REAL(dp)    :: pm6rp(kbdim,klev,nmod),    prhop(kbdim,klev,nmod),                 &
             paernl(kbdim,klev,nmod)

  REAL(dp)    :: pcom(kbdim,klev,nmod,nmod)

  !--- Local variables:
  
  INTEGER(i4) :: jm2, jm1, jl, jk
  !
  REAL(dp)    :: zpbyone,     zwlc,        zairvisc,    zbtketc,   &
             zrpvm1,      zrpvm2,      zrpav,       zpv1,      &
             zpm1,        zcv21,       zpv2,        zpm2,      &
             zcv22,       zcv2av,      zrknud1,     zrknud2,   &
             ze1,         ze2,         zpd21,       zpd22,     &
             zpdav,       zxd1,        zh21,        zxd2,      &
             zh22,        zh2av,       zcoc,        zhu1,      &
             zhu2,        zhu

  !--- 1) Calculation of the coagulation coefficient: -------------------------
  !
  DO jm2=1,nmod
     DO jm1=jm2,nmod

        ! mz_ak_20041014+
        IF (locoagmask(jm1,jm2)) THEN
           ! mz_ak_20041014-

        DO jk=1,klev
           DO jl=1,kproma
         !@@@ Included additional security checks to account for inconsistencies
              !    between mass, numbers, densities and radii:
              IF (paernl(jl,jk,jm1) > cmin_aernl  .AND. &
                   paernl(jl,jk,jm2) > cmin_aernl .AND. &
                   pm6rp(jl,jk,jm1)  > 1.e-10_dp .AND. &
                   pm6rp(jl,jk,jm2)  > 1.e-10_dp .AND. &
                   prhop(jl,jk,jm1)  > 1.e-10_dp .AND. &
                   prhop(jl,jk,jm2)  > 1.e-10_dp        ) THEN        
           
                 !--- 1.1) Calculate ambient properties:

                 !--- Mean free pathlength ? (from Knudsen Number below):
                 !    Parametrisation?
                 zpbyone=1000.0 / (papp1(jl,jk)/100.0)
                 zwlc=6.6e-6 * ptp1(jl,jk) / 293.15 * zpbyone

                 !--- Viscosity:
                 zairvisc=1.827e-4 * (ptp1(jl,jk) / 291.15)**0.74

                 !--- 
                 zbtketc=bk * ptp1(jl,jk) / 6.0 / pi / zairvisc

                 !--- Count median radii of the respective modes:
                 zrpvm1=pm6rp(jl,jk,jm1)
                 zrpvm2=pm6rp(jl,jk,jm2)

                 !--- Average radius of the modes:
                 zrpav=(zrpvm1 + zrpvm2) / 2.0

                 !--- Volume and mass of mode jm1:
                 zpv1=4.0 / 3.0 * pi * zrpvm1**3.0
                 zpm1=zpv1 * prhop(jl,jk,jm1)

                 !--- Squared mean particle velocity of mode jm1:
                 zcv21=8.0 * bk * ptp1(jl,jk) / (pi * zpm1)

                 !--- Volume and mass of particles in mode jm2:
                 zpv2=4.0 / 3.0 * pi * zrpvm2**3.0
                 zpm2=zpv2 * prhop(jl,jk,jm2)

                 !--- Squared mean particle velocity of mode jm2:
                 zcv22=8.0 * bk * ptp1(jl,jk) / (pi * zpm2)

                 !--- Fuchs: G_r (below Eq. 49.27):
                 zcv2av=SQRT(zcv21 + zcv22)

                 !--- Knudsen numbers of the modes:
                 !@@@ Check: is zwlc mean free path then zrknud=zwlc/zrpvm!
                 zrknud1=zwlc/zrpvm1/2.0
                 zrknud2=zwlc/zrpvm2/2.0

                 !--- Diffusivities of the respective modes:
                 !    !@@@ Parameterisation?
                 ze1=EXP(-0.43/zrknud1)
                 ze2=EXP(-0.43/zrknud2)
                 zpd21=zbtketc * (1.0_dp + zrknud1*2.492 + zrknud1*0.84*ze1) &
                                      / zrpvm1
                 zpd22=zbtketc * (1.0_dp + zrknud2*2.492 + zrknud1*0.84*ze2) &
                                      / zrpvm2

                 !--- Average diffusivity of the modes:
                 zpdav=(zpd21 + zpd22) / 2.0

                 !--- Average mean free path of particles in jm1:
                 zxd1=8.0 * zpd21 / (pi*SQRT(zcv21))

                 !- Mean distance from surface after mean free path (Eq. 49.13):
                 zh21=(((2.0_dp*zrpvm1 + zxd1)**3.0 -                        &
                       (4.0_dp*zrpvm1*zrpvm1 + zxd1*zxd1)**1.5) /            &
                       (6.0_dp*zrpvm1*zxd1) - 2._dp*zrpvm1   ) * sqrt2

                 !- Average mean free path of particles in jm2:
                 zxd2=8.0 * zpd22 / (pi*SQRT(zcv22))

                 !- Mean distance from surface after mean free path (Eq. 49.13):

                 zh22=(((2.0_dp*zrpvm2 + zxd2)**3.0 -                        &
                       (4.0_dp*zrpvm2*zrpvm2 + zxd2*zxd2)**1.5) /            &
                       (6.0_dp*zrpvm2*zxd2) - 2._dp*zrpvm2  ) * sqrt2

                 !--- Fuchs: delta_r !@@@ (why division by sqrt2?)
                 zh2av=SQRT(zh21*zh21 + zh22*zh22) / sqrt2
                 
                 !- 1.2) Calculation of the coagulation coefficient 
                 !   pcom (Eq. 49.26):
                 !   Factor 16 instead factor 8 as in Fuchs as his formulation
                 !   applies for the inter-modal coagulation. This is taken into
                 !   account in the assignment of the inter-modal coagulation
                 !    coefficient.

                 zcoc=16.0_dp * pi * zpdav * zrpav

                 !--- Calculation of beta=1/zhu (Eq. 49.27):
                 zhu1=4.0_dp * zpdav / (zcv2av * zrpav)
                 zhu2=zrpav / (zrpav + zh2av /2.0_dp)
                 zhu=zhu1 +  zhu2

                 !--- Coagulation coefficient following (Eq.49.26):
                 pcom(jl,jk,jm1,jm2)=zcoc / zhu

              ELSE
                 pcom(jl,jk,jm1,jm2)=0._dp
              END IF

              !--- 2) Mirror the coagulation matrix (symmetric): --------------

              pcom(jl,jk,jm2,jm1)=pcom(jl,jk,jm1,jm2)

           END DO
        END DO

        ! mz_ak_20041014+
        ELSE

           pcom(1:kproma,:,jm1,jm2)=0.0_dp

        END IF ! locoagmask
        ! mz_ak_20041014-

     END DO
  END DO

END SUBROUTINE m7_coaset
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
SUBROUTINE m7_nuck(kproma, kbdim,   klev,  pso4g,   &
                   ptp1,   prelhum, panew, pa4delt, &
                   ztmst                            )
  !      
  !  Authors:
  !  --------
  !  J. Wilson, E. Vignati, JRC/EI, original source                 09/2000
  !  P. Stier, MPI                  f90-version, 
  !                                 changes, 
  !                                 comments,
  !                                 modularisation and 
  !                                 implementation of 
  !                                 Vehkamaeki (2002)             2001-2003
  !
  !  Purpose:                                                           
  !  --------
  !  This routine calls the routines for the computation of the 
  !  nucleation rate znucrate [molec. cm-3 s-1] and, for 
  !  Vehkamaeki (2002), of the number of molecules in the critical
  !  cluster from a given gas phase H2SO4 concentration 
  !  pso4g [molec. cm-3]. It also calculates the integrated change of 
  !  H2SO4 gas phase mass over one timestep due to nucleation 
  !  pa4delt(:,:,1) [molec. cm-3] as well as the number of nucleated 
  !  particles panew [1] during the timestep. Whilst this is done 
  !  analytically for the old Kulmala (1998) parameterization, it has
  !  to be done numerically for the new Vehkamaeki (2002) scheme.
  !
  ! mz_ak_20041014+
    ! The old parameterization of Kulmala (1998), that apparently contains
  ! some errors is kept for consistency. It is recommended to use the 
  ! new scheme of Vehkamaeki et al. (2002). However, this parameterization
  ! is no longer analytically integrable, the effect of this is to be 
  ! scrutinized.
  ! The modularized version of the Kulmala parameterization has been tested
  ! to give bit-identical results with the old version.
  ! 
  ! References:
  ! -----------
  ! Vehkamaeki et al. (2002), An improved parameterization for sulfuric
  !    acid/water nucleation rates for tropospheric and stratospheric
  !    conditions, J. Geophys. Res, 107, D22, 4622
  ! Kulmala et al. (1998), Parameterizations for sulfuric acid/water
  !    nucleation rates. J. Geophys. Res., 103, No D7, 8301-8307
  ! Vignatti, E. (1999), Modelling Interactions between Aerosols and 
  !    Gaseous Compounds in the Polluted Marine Atmosphere. PhD-Thesis,
  !    RISO National Laborartory Copenhagen, Riso-R-1163(EN)
  !

  ! mz_ak_20041014- !qqq
  !  Interface:
  !  ----------
  !  *m7_nuck* is called from *m7_dnum*
  !
  !  Method:
  !  -------
  !
  !  1) Kulmala et al. (1998):
  !
  !  For the Kulmala et al. (1998) scheme the formula for binary
  !  nucleation is rewritten to take the form 
  !  znucrate = exp[zalpha+ln(pso4g)*beta]. 
  !  Equation numbers are taken from Kulmala et al. (1998).
  !  After the calculation of the nucleation rate znucrate, it is 
  !  integrated in 2) analytically over one timestep, i.e.:
  !
  !  Integration of:
  ! 
  !  znucrate=d(critn*znav)/dt=exp[zalpha + zbeta*ln(znav)]
  !
  !  gives znav(t0+dt, znav(t0)) where znav(t0) is pso4g.
  !  znav is temporarily stored in zso4g_new and confined between
  !  0 and pso4g. 
  !  The number of nucleated particles is then calculated by 
  !  assuming a fixed critical mass critn of newly nucleated 
  !  particles and dividing the total mass of nucleated sulfate 
  !  by this critical mass of one nucleated particle. 
  !
  !
  !  2) Vehkamaeki et al. (2002):
  !  
  !  An analytical integration of the nucleation equation is not
  !  possible, therefor the nucleation rate is simply multiplied
  !  by dt, implying a fixed gas-phase H2SO4 concentration over 
  !  the timestep.
  !@@@ The dependency of the results on the used timestep 
  !@@@ should be checked carefully in a sensitivity study! 
  !  The number of nucleated particles is then calculated by 
  !  taking the calculated critical mass critn of newly nucleated 
  !  particles and dividing the total mass of nucleated sulfate 
  !  by this critical mass of one nucleated particle. 
  !  
  !  Externals:
  !  ----------
  !  nucl_kulmala
  !  nucl_vehkamaeki
  !
  !  References:
  !  -----------
  !  Vehkamaeki et al. (2002), An improved parameterization for sulfuric
  !     acid/water nucleation rates for tropospheric and stratospheric
  !     conditions, J. Geophys. Res, 107, D22, 4622
  !  Kulmala et al. (1998), Parameterizations for sulfuric acid/water
  !     nucleation rates. J. Geophys. Res., 103, No D7, 8301-8307
  !  Vignatti, E. (1999), Modelling Interactions between Aerosols and 
  !     Gaseous Compounds in the Polluted Marine Atmosphere. PhD-Thesis,
  !     RISO National Laborartory Copenhagen, Riso-R-1163(EN)
  !
  !
  !--- Parameters:
  !
  ! pso4g          = mass of gas phase sulfate [molec. cm-3]
  ! ptp1           = atmospheric temperature at time t+1 [K]
  ! prelhum        = atmospheric relative humidity [%]
  ! pa4delt(:,:,1) = mass of H2SO4 added to the nucleation mode due 
  !                  to nucleation of H2SO4 over ztmst. 
  !                  Equilvalent to the integral of H2SO4 gas loss 
  !                  due to nucleation over timestep ztmst. [molec. cm-3]
  ! panew          = number of nucleated particles over timestep ztmst
  !                  panew=pa4delt/critn i.e. mass of formed sulfate 
  !                  divided by an assumed mass of a nucleus. [cm-3]
  !        
  !--- Local variables:
  !
  ! zso4g_new        = temporay storage of gas phase sulfate [molec. cm-3]
  !
  ! See comments!

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: pso4g(kbdim,klev),          ptp1(kbdim,klev),     &
             prelhum(kbdim,klev),        panew(kbdim,klev)

  REAL(dp)    :: pa4delt(kbdim,klev,naermod)

  ! Local variables:

  INTEGER(i4) :: jk,          jl
  REAL(dp)    :: ztmst,          zf1,          zeps !qqq

  REAL(dp)    :: znucrate(kbdim,klev),  & ! nucleation rate [m-3 s-1] 
                                          !@@@ Make consistent with Kulmala
             zncrit(kbdim,klev),    & ! number of molecules in the 
                                      ! critical cluster [1]
             zso4g_new(kbdim,klev), & ! new gas phase sulfate 
                                      ! concentration [molec. cm-3]
             zalpha(kbdim,klev),    & ! auxiliary coefficient for the analytical
                                      ! integration of Kulmala
             zbeta(kbdim,klev)        ! auxiliary coefficient for the analytical
                                      ! integration of Kulmala

  !--- 0) Initialisations: ----------------------------------------------------
  zeps=EPSILON(1.0_dp)

  !--- 1) Calculate nucleation rate:

  IF(nnucl==1) THEN

!CDIR NOIEXPAND
     CALL nucl_vehkamaeki(kproma,   kbdim,   klev,  & ! ECHAM5 dimensions
                          ptp1,     prelhum, pso4g, & ! ECHAM5 temperature, 
                                                      ! relative humidity
                          znucrate, zncrit          )  


     !--- Calculate updated gas phase concentration:
     !
     !    N(t)   = N(0)   - znucrate   * zncrit * dt
     !    [cm-3] = [cm-3] - [cm-3 s-1] * [1]    * [s] 

     DO jk=1, klev
        DO jl=1, kproma

           zso4g_new(jl,jk)=pso4g(jl,jk)-(znucrate(jl,jk)*zncrit(jl,jk)*ztmst)

        END DO
     END DO

  ELSE IF (nnucl==2) THEN

     zncrit(:,:)=critn

!CDIR NOIEXPAND
     CALL nucl_kulmala(kproma,  kbdim,   klev,    &
                       pso4g,   ptp1,    prelhum, &
                       znucrate, zalpha, zbeta    )

     !--- 2) Analytical integration of the nucleation rate (Eq. 19) -----------
     !       over timestep ztmst assuming no new H2SO4(g) production:
     !
     !       d(N_av/critn)/dt=exp(alpha + beta*ln(N_av)) => ... =>
     !           
     !       N_av(t0+dt)=
     !       [N_av(t0)**(1-beta) + critn*(1-beta)exp(alpha)*dt]**(1/(1-beta)
     !

     DO jk=1, klev
        DO jl=1, kproma

           IF (znucrate(jl,jk) .GT. 1e-10_dp) THEN
              zf1 = pso4g(jl,jk)**(1.0_dp-zbeta(jl,jk))-&
                          critn*EXP(zalpha(jl,jk))*(1.0_dp-zbeta(jl,jk))*ztmst
              zso4g_new(jl,jk) = EXP(LOG(zf1)/(1.0_dp - zbeta(jl,jk)))
           ELSE
              zso4g_new(jl,jk) = pso4g(jl,jk)
           END IF

        END DO
     END DO

  END IF

  !--- 3) Calculate changes in gas-phase and aerosol mass and aerosol numbers: -

  DO jk=1, klev
     DO jl=1, kproma

        ! mz_ak_20041014+
        IF( znucrate(jl,jk) > zeps  ) THEN
        ! mz_ak_20041014-

  !--- 3.1) Security check:

        zso4g_new(jl,jk) = MAX(zso4g_new(jl,jk), 0.0_dp)
        zso4g_new(jl,jk) = MIN(zso4g_new(jl,jk), pso4g(jl,jk))
   
  !--- 3.2) Calculate mass of nucleated H2SO4 (equals the net 
  !         gas phase H2SO4 loss):
  ! 
        pa4delt(jl,jk,1) = (pso4g(jl,jk)-zso4g_new(jl,jk))
  !
  !--- 3.3) Calculate the number of nucleated particles (nucleated mass 
  !         divided by the assumed mass of a critical cluster critn):

        panew(jl,jk)=pa4delt(jl,jk,1)/zncrit(jl,jk)      
  !
  !--- 3.4) Calculate changes in gas phase H2SO4 due to nucleation:
  !
        pso4g(jl,jk)=pso4g(jl,jk)-pa4delt(jl,jk,1)
     END IF !qqq
  END DO
END DO

END SUBROUTINE m7_nuck
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~?
SUBROUTINE m7_delcoa(kproma,   kbdim,  klev,     paerml,   &
                     paernl,   pm6rp,  pa4delt,  panew,    &
                     pa,       pb,     pbfract1, pbfract2, &
                     pbfract5, pso4_5, pso4_6,   pso4_7    ,   &
                      ztmst                                 ) !qqq
  ! 
  !    Authors: 
  !    --------- 
  !    E. Vignati and J. Wilson, JRC/EI (original source)                09/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)    2001 
  ! 
  !    Version: 
  !    --------- 
  !    This version is equivalent to the version delco_n2 of the m7 boxmodel. 
  !    + use of the analytical solution 
  ! 
  !    Purpose 
  !    --------- 
  !    This routine calculates changes in number concentration of 
  !    each aerosol mode over the time step, due to coagulation with 
  !    the current mode and all higher ones. 
  ! 
  !    Method: 
  !    ----------- 
  !    *delcoa*  integrates for each mode dn/dt=c -a*n^2 -b*n  over ztmst  
  ! 
  !    The resulting particles are assumed to reside in the 
  !    mode of highest mode of the pair of particles colliding. 
  !    1+1=>1 1+2=>2 1+3=>3, 2+2=>2 2+3=>3, 3+3=>3. 
  !    zc is now non zero for mode 1 only (nucleation). 
  !    All formation of higher mode particles is handled in dconc. 
  ! 
  !    For climatological studies, 5 day accumulation mode concs are  
  !    within a factor of 2 of the full model.  
  ! 
  !    Interface 
  !    ----------- 
  !    *m7_delcoa* is called from *m7_dnum* 
  ! 
  !    Externals 
  !    ----------- 
  !    none 
  ! 
  !--- Parameter list:
  !
  ! paerml          = total aerosol mass for each compound 
  !                   [molec. cm-3 for sulphate; ug m-3 for others]
  ! paernl          = aerosol number for each mode [cm-3]
  ! pm6rp           = mean mode actual radius (wet radius for soluble modes 
  !                   and dry radius for insoluble modes) [cm]
  ! pa4delt(:,:,:)  = change in H2SO4 mass of the respective mode over 
  !                   one timstep 
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in m7_nuck)
  !                      - coagulation (calculated in m7_concoag)
  ! panew           = number of nucleated particles (during 1 timestep) [1] 
  ! pa              = unimodal coagulation coefficient (zcom(mod)) []
  ! pb              = inter-modal coagulation with higher modes
  !                   (zcom(mod) * n(jmod+1)) []
  ! pbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x that is moved to mode y+1 [1]
  ! pso4_x          = mass of sulphate condensed on insoluble 
  !                   mode x [molec. cm-3]
  ! 
  !--- Local variables:
  !
  ! zansum          = aerosol number in the respective mode [cm-3]
  ! zxxsum          = aerosol mass for compound xx in the respective
  !                   mode, e.g. xx = bc, oc, a4 (sulfate) 
  !                   [g cm-3 for bc,oc and molec. cm-3 for sulfate]
  ! zxxav           = average mass of a sulfate particle in the respective
  !                   mode [molecules]
  ! zxxavy          = average mass of species xx in mode y []
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  !                   [molecules for sulfate and ug for others]
  ! zanli(:,:,:)    = Number of particles moved by the inter-modal 
  !                   coagulation []
  ! zansq(:,:,:)    = Number of particles moved by the intra-modal 
  !                   coagulation []
  ! zaernt(:,:,:)   = New particle number n(t+dt) after the integration
  !                   of the aerosol dynamics equation [cm-3]

  !--- Parameters:

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev
 
  REAL(dp)    :: paerml(kbdim,klev,naermod),   paernl(kbdim,klev,nmod),     & 
                 pa4delt(kbdim,klev,naermod),  panew(kbdim,klev),           &
                 pm6rp(kbdim,klev,nmod)
 
  REAL(dp)    :: pbfract1(kbdim,klev,nmod-1),  pbfract2(kbdim,klev,nmod-1),  & 
                 pbfract5(kbdim,klev,3)
 
  REAL(dp)    :: pa(kbdim,klev,nmod),          pb(kbdim,klev,nmod)
 
  REAL(dp)    :: pso4_5(kbdim,klev),           pso4_6(kbdim,klev),           & 
                 pso4_7(kbdim,klev)
 
  ! Local variables: 
 
  INTEGER(i4) :: jmod, jl, jk , kmod
  REAL(dp)    :: ztmst
 
  REAL(dp)    :: za4av1(kbdim,klev),           za4av2(kbdim,klev),          & 
             zbcav2(kbdim,klev),           zocav2(kbdim,klev),              &
             zbcav3(kbdim,klev),           zocav3(kbdim,klev),              & 
             zbcav4(kbdim,klev),           zocav4(kbdim,klev),              & 
             zbcav5(kbdim,klev),           zocav5(kbdim,klev),              & 
             zduav6(kbdim,klev),           zduav7(kbdim,klev) 
 
  REAL(dp):: zanli(kbdim,klev,nmod),       zansq(kbdim,klev,nmod),          &
             zaernt(kbdim,klev,nmod)

  REAL(dp):: zm6rp(nmod),                  zcrtcst(nmod) ! mz_ak_20041014

  ! Auxiliary variables: 
 
  REAL(dp)    :: zansum, zbcsum, zocsum, za4sum, za4av,                     & 
             ztop,   zbot,   zatot,  zanse,  zanle,                         & 
             ze1,    zf1,    zf2,    zf3,    zf4,    zr1,                   &
             zc
 
  ! mz_ak_20041014+
   REAL(dp):: zamloss5,      zamloss6,      zamloss7,      zanloss5,         &
             zanloss6,      zanloss7,      ztloss,        zbfofac,           &
             zbfnfac,       zbftot,        zaerni,        zanytnt,          &
             zanytni,       zanytns,       zanytnm,       ztotn,            &
             zaerns
 
  ! mz_ak_20041014-

  za4av       = 0._dp 
  za4av1(:,:) = 0._dp
  za4av2(:,:) = 0._dp
  zbcav2(:,:) = 0._dp 
  zocav2(:,:) = 0._dp 
  zbcav3(:,:) = 0._dp
  zocav3(:,:) = 0._dp
  zbcav4(:,:) = 0._dp
  zocav4(:,:) = 0._dp
  zbcav5(:,:) = 0._dp
  zocav5(:,:) = 0._dp
 
  !--- 1) Insoluble modes 
  ! 
  DO jk=1,klev 
     DO jl=1,kproma 

        !--- Dust modes insoluble: 
        !    Calculate average mass (used in concoag)

        IF(paernl(jl,jk,iacci) > 0.0_dp) THEN
           zduav6(jl,jk) = (paerml(jl,jk,iduai)*1.E-6) / paernl(jl,jk,iacci)
        ELSE
           zduav6(jl,jk)=0._dp
        END IF
        IF(paernl(jl,jk,icoai) > 0.0_dp) THEN
           zduav7(jl,jk) = (paerml(jl,jk,iduci)*1.E-6) / paernl(jl,jk,icoai)
        ELSE
           zduav7(jl,jk)=0._dp
        END IF

        !--- Aitken mode insoluble:
        !    Only considered process: 
        !    Coagulation and transfer from the insoluble aitken mode 
        !    to the soluble aitken and accumulation modes.

        zansum=paernl(jl,jk,iaiti) 
        zbcsum=paerml(jl,jk,ibcki)*1.e-6 
        zocsum=paerml(jl,jk,iocki)*1.e-6 
        zaernt(jl,jk,iaiti)=zansum 
        zanli(jl,jk,iaiti)=0.0_dp 
        zansq(jl,jk,iaiti)=0.0_dp 

        !--- Calculations only in presence of sufficient particles:

        IF (zansum > cmin_aernl) THEN 
 
           ! --- Average mass for a bc/oc particle in the aitken mode [ug]:

           zbcav5(jl,jk)=zbcsum/zansum 
           zocav5(jl,jk)=zocsum/zansum 

           !--- 1.1) Case of no coagulation:
           !
           !         (pa(jl,jk,iaiti) <  cmin_aerml .AND. pb(jl,jk,iaiti) < cmin_aerml) 
           !
           !         => Nothing to be done
              
           !--- 1.2) Case with coagulation:

           IF (pa(jl,jk,iaiti) >= 1.e-15_dp .OR. &
                                      pb(jl,jk,iaiti) >=1.e-15_dp ) THEN
                 
              !--- 1.2.1) Case of no inter-modal coagulation:
              !           dn/dt = -a*n**2 => 
              !           n(t)  = n0/(1 + n0*a*(t-t0))

              IF (pb(jl,jk,iaiti) < 1.e-15_dp) THEN 
                 zaernt(jl,jk,iaiti)=&
                      zansum/(1.0_dp+zansum*pa(jl,jk,iaiti)*ztmst) 
                 zanli(jl,jk,iaiti)=0.0_dp 
                 zansq(jl,jk,iaiti)=zansum-zaernt(jl,jk,iaiti) 

              !--- 1.2.2) Case with inter- and intra-modal coagulation:
              !         dn/dt = -a*n**2 - b*n => 
              !         n(t)  = (b*n0*exp(-b(t-t0)))/((n0*a)(1-exp(-b(t-t0)))+b)
              ELSE   

                 !--- Calculate n(t+dt):

                 ze1=EXP(-pb(jl,jk,iaiti)*ztmst) 
                 ztop=pb(jl,jk,iaiti)*zansum*ze1 
                 zbot=zansum*pa(jl,jk,iaiti)*(1.0_dp-ze1)+pb(jl,jk,iaiti) 
                 zaernt(jl,jk,iaiti)=ztop/zbot 
                 !--- Limit n(t+dt) to available particle in the mode:
                 zaernt(jl,jk,iaiti)=MIN(zaernt(jl,jk,iaiti), zansum) 

                 !--- Total change in particle numbers of the mode due to 
                 !    coagulation:
                 zatot=zansum-zaernt(jl,jk,iaiti) 
                 !--- Contribution of the intra-modal coagulation:
                 zanse=zansum*zansum*pa(jl,jk,iaiti)
                 !--- Contribution of the inter-modal coagulation:
                 zanle=zansum*pb(jl,jk,iaiti) 
                 !--- Number of particles moved by the inter-modal coagulation:
                 zanli(jl,jk,iaiti)=zatot*zanle/(zanse+zanle) 
                 !--- Number of particles moved by the intra-modal coagulation:
                 zansq(jl,jk,iaiti)=zatot*zanse/(zanse+zanle) 

              END IF 
              !--- 1.2.3) Change masses of the insoluble aitken mode due to 
              !           intra-modal coagulation and the coagulation with the
              !           nucleation mode (transfers to the soluble modes
              !           of the particles coagulating with the nucleation mode 
              !           are done in m7_concoag):
 
              paerml(jl,jk,ibcki)=( zaernt(jl,jk,iaiti) +                    &
                                    zansq(jl,jk,iaiti)  +                    &
                                    pbfract5(jl,jk,1)*zanli(jl,jk,iaiti) ) * & 
                                    zbcav5(jl,jk)*1.e6 

              paerml(jl,jk,iocki)=( zaernt(jl,jk,iaiti) +    &
                                    zansq(jl,jk,iaiti)  +                    &
                                    pbfract5(jl,jk,1)*zanli(jl,jk,iaiti) ) * & 
                                    zocav5(jl,jk)*1.e6 

              !-- 1.2.4) Change the numbers of the insoluble aitken mode due to 
              !           intra-modal coagulation:
 
              paernl(jl,jk,iaiti)=zaernt(jl,jk,iaiti)   +       &
                                  pbfract5(jl,jk,1)*zanli(jl,jk,iaiti)
 
              !--- 1.2.5) Store changes in masses of compounds in the insoluble 
              !           aitken mode due to inter-modal coagulation:
              !           (zanli(:,:,x)   = total number of particles moved 
              !                               from mode x
              !            pbfract5(:,:,x)= fraction of the total number of
              !                              particles
              !                             moved from mode 5 that is moved 
              !                              to mode x   )
 
              pa4delt(jl,jk,ibcks)=pbfract5(jl,jk,2)*zanli(jl,jk,iaiti)*&
                                        zbcav5(jl,jk)*1.e6 
              pa4delt(jl,jk,iocks)=pbfract5(jl,jk,2)*zanli(jl,jk,iaiti)*&
                                        zocav5(jl,jk)*1.e6 
              pa4delt(jl,jk,ibcas)=pbfract5(jl,jk,3)*zanli(jl,jk,iaiti)*&
                                        zbcav5(jl,jk)*1.e6 
              pa4delt(jl,jk,iocas)=pbfract5(jl,jk,3)*zanli(jl,jk,iaiti)*&
                                        zocav5(jl,jk)*1.e6 
 
           END IF
        END IF
     END DO
  END DO
  ! 
  !--- 2) Soluble modes: --------------------------------------------------
  ! 
!CDIR unroll=5
  mode : DO jmod=1,nsol 
     level : DO jk=1,klev 
        longitude : DO jl=1,kproma 

           !--- Nucleation mode:
           IF (jmod .EQ. 1) THEN 
              zansum=paernl(jl,jk,jmod)+panew(jl,jk) 
              za4sum=paerml(jl,jk,jmod)+pa4delt(jl,jk,1)
           !--- Others:
           ELSE 
              zansum=paernl(jl,jk,jmod) 
              za4sum=paerml(jl,jk,jmod)
           END IF 

           IF (jmod.EQ.2) THEN 
              zbcsum=paerml(jl,jk,ibcks)*1.e-6 
              zocsum=paerml(jl,jk,iocks)*1.e-6 
           END IF 
           IF (jmod.EQ.3) THEN 
              zbcsum=paerml(jl,jk,ibcas)*1.e-6 
              zocsum=paerml(jl,jk,iocas)*1.e-6 
           END IF 
           IF (jmod.EQ.4) THEN 
              zbcsum=paerml(jl,jk,ibccs)*1.e-6 
              zocsum=paerml(jl,jk,ioccs)*1.e-6 
           END IF 

           zaernt(jl,jk,jmod)=zansum 
           zanli(jl,jk,jmod)=0.0_dp 
           zansq(jl,jk,jmod)=0.0_dp  

           !--- Calculations only in presence of sufficient particles:

           IF (zansum > cmin_aernl) THEN 

              za4av=za4sum/zansum 

              IF (jmod.EQ.1) THEN 
                 za4av1(jl,jk)=za4av 
              ELSE IF (jmod.EQ.2) THEN 
                 zbcav2(jl,jk)=zbcsum/zansum 
                 zocav2(jl,jk)=zocsum/zansum 
                 za4av2(jl,jk)=za4av 
              ELSE IF (jmod.EQ.3) THEN 
                 zbcav3(jl,jk)=zbcsum/zansum 
                 zocav3(jl,jk)=zocsum/zansum                  
              ELSE IF (jmod.EQ.4) THEN 
                 zbcav4(jl,jk)=zbcsum/zansum 
                 zocav4(jl,jk)=zocsum/zansum 
              END IF

              !--- 2.1) Case of no coagulation:
              !         
              IF (pa(jl,jk,jmod) <  1.e-15_dp .AND. &
                                   pb(jl,jk,jmod) < 1.e-15_dp) THEN 
 
                 !--- Nucleation in mode 1 only. 
                 !    Nothing to be done for other modes.

                 IF(jmod.EQ.1) THEN 
                    paerml(jl,jk,jmod)=za4sum
                    paernl(jl,jk,jmod)=zansum
                 END IF 
 
              !--- 2.2) Case with coagulation:

              ELSE 

                 !--- 2.2.1) Case of no nucleation:

                 !--- Not Mode 1 or Nucleation rate below 1/s:

                 IF ( (jmod .NE. 1) .OR.&
                              (panew(jl,jk)/ztmst < 1.0_dp) ) THEN

                    paernl(jl,jk,jmod)=zansum 

                    !--- 2.2.1a) Case of no inter-modal coagulation:
                    !            dn/dt = -a*n**2 => 
                    !            n(t)  = n0/(1 + n0*a*(t-t0))

                    IF (pb(jl,jk,jmod) < 1.e-15_dp) THEN 
                       zaernt(jl,jk,jmod)=zansum/&
                                      (1.0_dp+zansum*pa(jl,jk,jmod)*ztmst) 
                       zanli(jl,jk,jmod)=0.0_dp 
                       zansq(jl,jk,jmod)=zansum-zaernt(jl,jk,jmod) 

                    !--- 2.2.1b) Case with inter- and intra-modal coagulation:
                    !   dn/dt = -a*n**2 - b*n => 
                    !   n(t)  = (b*n0*exp(-b(t-t0)))/((n0*a)(1-exp(-b(t-t0)))+b)

                    ELSE            
                       !--- Calculate n(t+dt):
                       ze1=EXP(-pb(jl,jk,jmod)*ztmst) 
                       ztop=pb(jl,jk,jmod)*zansum*ze1 
                       zbot=zansum*pa(jl,jk,jmod)*(1.0_dp-ze1)+pb(jl,jk,jmod) 
                       zaernt(jl,jk,jmod)=ztop/zbot 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum) 
                       !--- Total change in particle numbers of the mode 
                       !    due to coagulation:
                       zatot=zansum-zaernt(jl,jk,jmod) 
                       !--- Contribution of the intra-modal coagulation:
                       zanse=zansum*zansum*pa(jl,jk,jmod)
                       !--- Contribution of the inter-modal coagulation:
                       zanle=zansum*pb(jl,jk,jmod) 
                       !--- Number of particles moved by the inter-modal 
                       !    coagulation:
                       zanli(jl,jk,jmod)=zatot*zanle/(zanse+zanle) 
                       !--- Number of particles moved by the intra-modal 
                       ! coagulation:
                       zansq(jl,jk,jmod)=zatot*zanse/(zanse+zanle) 
                    END IF 
                 
                 !--- 2.2.2) Case with nucleation:

                 ELSE IF ( (jmod .EQ. 1) .AND. &
                            (panew(jl,jk)/ztmst >= 1.0_dp) ) THEN

                 !--- 2.2.2a) Nucleation, inter- and intra-modal coagulation:
                 ! dn/dt = -a*n**2 - b*n + c => 
                 ! n(t)  = -(b/(2a)) + 
                 ! R/2a * [ ((1 - (-2ax0-b+R)/(+2ax0+b+R))exp(-Rt)) /
                 !           ((1 + (-2ax0-b+R)/(+2ax0+b+R))exp(-Rt))  ]
                 ! where:  R=SQRT(b**2+4ac)
                 !
                 ! If b/=0 then always a/=0. The only case where a would be 0
                 ! and b unequal zero is the case of no pre-existing particles 
                 ! in the nucleation mode but pre-existing particles in other
                 ! modes. For this case a is calculated for an assumed radius
                 ! of a critical cluster in m7_coaset. 

                    IF (pb(jl,jk,jmod) >= 1.e-15_dp) THEN
                       !--- Calculate n(t):
                       !--- c:
                       zc=panew(jl,jk)/ztmst
                       !--- R:
                       zf1=pb(jl,jk,jmod)*pb(jl,jk,jmod)+4.0*pa(jl,jk,jmod)*zc 
                       zr1=SQRT(zf1) 
                       !--- exp(-Rt):
                       ze1=EXP(-zr1*ztmst) 
                       !--- 2ax0+b:
                       zf2=2.0*pa(jl,jk,jmod)*paernl(jl,jk,jmod)+pb(jl,jk,jmod) 
                       !--- Term in squared bracket:
                       zf3=ze1*(zr1-zf2)/(zr1+zf2) 
                       zf4=(1.0_dp-zf3)/(1.0_dp+zf3) 
                       !--- n(t):
                       zaernt(jl,jk,jmod)=(zr1*zf4-pb(jl,jk,jmod))/&
                                                      2.0/pa(jl,jk,jmod) 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum) 
                       !--- Total change in particle numbers of the mode 
                       !     due to coagulation:
                       zatot=zansum-zaernt(jl,jk,jmod) 
                       !--- Contribution of the intra-modal coagulation:
                       zanse=zansum*zansum*pa(jl,jk,jmod)
                       !--- Contribution of the inter-modal coagulation:
                       zanle=zansum*pb(jl,jk,jmod) 
                       !--- Number of particles moved by the inter-modal 
                       !    coagulation:
                       zanli(jl,jk,jmod)=zatot*zanle/(zanse+zanle) 
                       !--- Number of particles moved by the intra-modal 
                       !    coagulation:
                       zansq(jl,jk,jmod)=zatot*zanse/(zanse+zanle) 

                 !--- 2.2.2b) Nucleation and intra-modal coagulation:
                 !        dn/dt = -a*n**2 - b*n + c with b=0 =>
                 !        dn/dt = -a*n**2 + c => 
                 !        n(t)  = R/2a * [ ((1 - (-2ax0+R)/(+2ax0+R))exp(-Rt)) /
                 !                       ((1 + (-2ax0+R)/(+2ax0+R))exp(-Rt))  ]
                 !        where:  R=SQRT(4ac)
                 !        Can be shown to be equivalent to:
                 !
                 !        n(t)  = R1*((x0+R1)/(x0-R1)+exp(-SQRT(-4ac)t)) / 
                 !                    ((x0+R1)/(x0-R1)-exp(-SQRT(-4ac)t))
                 !        where R1=SQRT(c/a)

                    ELSE IF (pb(jl,jk,jmod) < 1.e-15_dp) THEN 
                       !--- c:
                       zc=panew(jl,jk)/ztmst
                       !--- R1:
                       zr1=SQRT(zc/pa(jl,jk,jmod)) 
                       !--- exp(-Rt):
                       ze1=EXP(-zr1*2.0*pa(jl,jk,jmod)*ztmst)
                       !--- n(t):
                       zf1=(paernl(jl,jk,jmod)+zr1)/(paernl(jl,jk,jmod)-zr1) 
                       ztop=zr1*(zf1+ze1) 
                       zbot=zf1-ze1 
                       IF (zbot < 1.e-15_dp) THEN 
                          zaernt(jl,jk,jmod)=zansum 
                       ELSE 
                          zaernt(jl,jk,jmod)=ztop/zbot 
                       END IF 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum) 
                       !--- Number of particles moved by the inter-modal 
                       !    coagulation:
                       zanli(jl,jk,jmod)=0.0_dp 
                       !--- Number of particles moved by the intra-modal
                       !    coagulation:
                       zansq(jl,jk,jmod)=zansum-zaernt(jl,jk,jmod) 
                    END IF 
                 END IF 
                 ! mz_ak_20041014+
                  !---2.2.3 New bit for insoluble/souble coagulation
                 !--- sum total insoluble+soluble paticles in mode jmod JJNW
                 IF (jmod .EQ. 1 .AND. zanli(jl,jk,jmod)>0.0_dp) THEN
                 zaerni=paernl(jl,jk,iaiti)+paernl(jl,jk,iacci)+&
                      paernl(jl,jk,icoai)
                 zaerns=zansum+paernl(jl,jk,iaits)
!                 zaerns=zansum
                 ztotn=zaerns+zaerni
                 IF (zaerns .gt. zaerni .and. zaerni .gt. 0.0_dp) THEN
                    ! calculate analytical solution no of mixed particles for 
                    ! coagulation between paernl(jl,jk,jmod) soluble particles 
                    ! and zaerni insouble of the same dimensions
                    IF (zaerni .gt. 1.0_dp) then
                       zanytni=4.0_dp*zaerni/&
                            ((2.0_dp+pa(jl,jk,jmod)*ztmst*ztotn)*          &
                            (2.0_dp+pa(jl,jk,jmod)*ztmst*(ztotn-zaerni)))
                    ELSE
                      zanytni = 0.0_dp
                    ENDIF
                    zanytnt=2.0_dp*ztotn/(2.0_dp+pa(jl,jk,jmod)*ztmst*ztotn)
                    zanytns=4.0_dp*zaerns/((2.0_dp+pa(jl,jk,jmod)*ztmst*ztotn)*&
                            (2.0_dp+pa(jl,jk,jmod)*ztmst*(ztotn-zaerns)))
                    zanytnm=zanytnt-(zanytni+zanytns)

                    !scale analytical solution to real aernt
                    zanytnm=min(zanytnm,zaerni)
                    zanytni=zaerni-zanytnm
                    zanytns=zaernt(jl,jk,jmod)
!CDIR UNROLL=7
                    DO kmod=1,nmod
                       zm6rp(kmod)=pm6rp(jl,jk,kmod)
                    END DO
                    CALL m7_coat(zm6rp,zcrtcst)
                    zamloss5=paernl(jl,jk,5)/zaerni*zanytnm*zcrtcst(5)
                    zanloss5=zamloss5/za4av
                    zamloss6=paernl(jl,jk,6)/zaerni*zanytnm*zcrtcst(6)
                    zanloss6=zamloss6/za4av
                    zamloss7=paernl(jl,jk,7)/zaerni*zanytnm*zcrtcst(7)
                    zanloss7=zamloss7/za4av
                    ztloss=zanloss5+zanloss6+zanloss7

                    ztloss=min(ztloss,zansq(jl,jk,jmod)*0.95_dp)           
                    zbfofac=zanli(jl,jk,jmod)/(zanli(jl,jk,jmod)+ztloss)
                    zbfnfac=ztloss/(zanli(jl,jk,jmod)+ztloss)
                    zanli(jl,jk,jmod)=zanli(jl,jk,jmod)+ztloss
                    zansq(jl,jk,jmod)=zansq(jl,jk,jmod)-ztloss
                    zbftot=0.0_dp
!CDIR UNROLL=7
                    DO kmod=1,nmod
                       IF(kmod>jmod) THEN
                          pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)*zbfofac
                          IF (kmod.GE.5) THEN
                             pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)+                  &
                                                       zbfnfac*paernl(jl,jk,kmod)/zaerni
                          END IF
                          zbftot=zbftot+pbfract1(jl,jk,kmod-jmod)
                       END IF
                    END DO       
!CDIR UNROLL=7
                    DO kmod=1,nmod
                       IF (kmod>jmod) THEN
                          pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)/zbftot
                       END IF
                    END DO
                 END IF
                 END IF
                 !---- End of new inslouble/soluble caogulation routine JJNW
                
                 !--- 2.3) Change masses and numbers of the respective modes 
                 !           to account-----------
                 !      for intra-modal coagulation (zansq) and coagulation with
                 !      higher modes (zaernt):
                 !
                 !--- 2.3.1) Change mass of the sulfur compounds:

                 paerml(jl,jk,jmod)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*za4av
 
                 !--- 2.3.2) Change mass of the carbon compounds:

                 IF (jmod.EQ.2) THEN 
                    paerml(jl,jk,ibcks)=(zaernt(jl,jk,jmod)+&
                                          zansq(jl,jk,jmod))*zbcav2(jl,jk)*1.e6 
                    paerml(jl,jk,iocks)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*&
                                            zocav2(jl,jk)*1.e6 
                 ELSE IF (jmod.EQ.3) THEN 
                    paerml(jl,jk,ibcas)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*&
                                            zbcav3(jl,jk)*1.e6 
                    paerml(jl,jk,iocas)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*&
                                            zocav3(jl,jk)*1.e6 
                 ELSE IF (jmod.EQ.4) THEN 
                    paerml(jl,jk,ibccs)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*&
                                            zbcav4(jl,jk)*1.e6 
                    paerml(jl,jk,ioccs)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*&
                                            zocav4(jl,jk)*1.e6 
                 END IF 

                 !-- 2.3.3) Particle numbers:
 
                 paernl(jl,jk,jmod)=zaernt(jl,jk,jmod)!@@+zansq(jl,jk,jmod)/2.0 
 
                 !- 2.4) Calculate changes in particle masses due to inter-modal
                 !         coagulation:

                 !--- 2.4.1) Transfer of mass from mode 1 to other modes:

                 IF (jmod .EQ. 1) THEN 
                     
                    ! Mass from 1 to 2: 
 
                    pa4delt(jl,jk,2)=pbfract1(jl,jk,1)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 2 due to coag. with 5:
 
                    pa4delt(jl,jk,2)=pa4delt(jl,jk,2)+&
                                        pbfract1(jl,jk,4)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 3:
 
                    pa4delt(jl,jk,3)=pbfract1(jl,jk,2)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 3 due to coag. with 6:
 
                    pa4delt(jl,jk,3)=pa4delt(jl,jk,3)+pbfract1(jl,jk,5)*&
                                        zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 4: 
 
                    pa4delt(jl,jk,4)=pbfract1(jl,jk,3)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 4 due to coag. with 7:
 
                    pa4delt(jl,jk,4)=pa4delt(jl,jk,4)+pbfract1(jl,jk,6)*&
                                      zanli(jl,jk,1)*za4av 

                 !---  2.4.2) Transfer of mass from mode 2 to other modes:
 
                 ELSE IF (jmod .EQ. 2) THEN 
 
                    ! Mass from 2 to 3: 
                     
                    pa4delt(jl,jk,3)=pa4delt(jl,jk,3)+pbfract2(jl,jk,2)*&
                                                       zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibcas)=pa4delt(jl,jk,ibcas)+                & 
                            pbfract2(jl,jk,2)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.e6 
                    pa4delt(jl,jk,iocas)=pa4delt(jl,jk,iocas)+                & 
                            pbfract2(jl,jk,2)*zanli(jl,jk,2)*zocav2(jl,jk)*1.e6 
 
                    ! Mass from 2 to 3 due to coag. with 6: 
 
                    pa4delt(jl,jk,3)=pa4delt(jl,jk,3)+pbfract2(jl,jk,5)*&
                                         zanli(jl,jk,2)*za4av 
 
                    ! mz_ak_20041102+
                    pa4delt(jl,jk,ibcas)=pa4delt(jl,jk,ibcas)+       & 
                         pbfract2(jl,jk,5)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.e6 
                    pa4delt(jl,jk,iocas)=pa4delt(jl,jk,iocas)+       & 
                         pbfract2(jl,jk,5)*zanli(jl,jk,2)*zocav2(jl,jk)*1.e6 
                    ! mz_ak_20041102-

                   ! Mass from 2 to 4: 
 
                    pa4delt(jl,jk,4)=pa4delt(jl,jk,4)+pbfract2(jl,jk,3)*&
                                              zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibccs)=pa4delt(jl,jk,ibccs)+                & 
                        pbfract2(jl,jk,3)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.E6 
                    pa4delt(jl,jk,ioccs)=pa4delt(jl,jk,ioccs)+                & 
                        pbfract2(jl,jk,3)*zanli(jl,jk,2)*zocav2(jl,jk)*1.E6 

                    ! Mass from 2 to 4 due to coag. with 7: 
 
                    pa4delt(jl,jk,4)=pa4delt(jl,jk,4)+&
                                         pbfract2(jl,jk,6)*zanli(jl,jk,2)*za4av 
                    ! mz_ak_20041102+
                    pa4delt(jl,jk,ibccs)=pa4delt(jl,jk,ibccs)+          & 
                         pbfract2(jl,jk,6)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.E6 
                    pa4delt(jl,jk,ioccs)=pa4delt(jl,jk,ioccs)+          & 
                         pbfract2(jl,jk,6)*zanli(jl,jk,2)*zocav2(jl,jk)*1.E6 
                    ! mz_ak_20041102-

                    ! Mass from 2 due to coagulation of 2 with 5 remains in 2:
                    !
                    !@@@ No effect as pbfract2(:,:,4)=0._dp!
                    !@@@ (Not needed as it is assumed that 5 coagulates with 2 
                    !@@@ and therefor the masses in 2 remain unchanged!)

                    pa4delt(jl,jk,2)=pa4delt(jl,jk,2)+&
                             pbfract2(jl,jk,4)*zanli(jl,jk,2)*za4av
                    pa4delt(jl,jk,ibcks)=pa4delt(jl,jk,ibcks)+                &
                            pbfract2(jl,jk,4)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.E6
                    pa4delt(jl,jk,iocks)=pa4delt(jl,jk,iocks)+                &
                            pbfract2(jl,jk,4)*zanli(jl,jk,2)*zocav2(jl,jk)*1.E6
 
                 END IF
              END IF
           END IF
        END DO longitude 
     END DO level 
  END DO mode 

  !--- 3) Calculate transfer from the insoluble to the soluble modes: ---------
 

!CDIR NOIEXPAND
  CALL m7_concoag (kproma,   kbdim,   klev,                      &
                   paerml,   paernl,  pm6rp,  pa4delt, zanli,    & 
                   za4av1,   za4av2,  zbcav5, zocav5,  zduav6,   & 
                   zduav7,   pso4_5,  pso4_6, pso4_7,            &
                   pbfract1, pbfract2                            )
 
 
  !--- 4) Final change of the aerosol masses due to nucleation, ----------------
  !       inter-modal coagulation and condensation on the insoluble modes:
  !       (Nucleation mode already done above.)

  DO jmod=2,naermod 
     paerml(1:kproma,:,jmod)=paerml(1:kproma,:,jmod)+pa4delt(1:kproma,:,jmod) 
  END DO 
   
END SUBROUTINE m7_delcoa 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
SUBROUTINE m7_concoag (kproma,   kbdim,   klev,                     &
                       paerml,   paernl,  pm6rp,  pa4delt, panli,   &
                       pa4av1,   pa4av2,  pbcav5, pocav5,  pduav6,  &
                       pduav7,   pso4_5,  pso4_6, pso4_7,           &
                       pbfract1, pbfract2                           )
  !
  !   *m7_concoag*
  !
  !   Author:
  !   ----------
  !   E. Vignati, JRC/EI     (original source)                09/2000
  !   P. Stier, MPI          (f90-version, changes, comments)    2001 

  !   Version:
  !   ----------
  !   This version is equivalent to the version concoa_n of the boxmodel. 
  !
  !   Purpose
  !   ----------
  !   m7_concoag transfers aerosol mass and numbers from the insoluble
  !   to the soluble modes.
  !
  !   Interface:
  !   ----------
  !   *m7_concoag* is called from *m7_delcoa*
  !
  !   Externals
  !   ----------
  !   none
  !
  !--- Parameters:
  !
  ! paerml          = total aerosol mass for each compound 
  !                   [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and 
  !                    dust]
  ! paernl          = aerosol number for each mode [cm-3]
  ! pm6rp           = mean mode actual radius (wet radius for soluble modes 
  !                   and dry radius for insoluble modes) [cm]
  ! pa4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one 
  !                   timstep 
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in m7_nuck)
  !                      - coagulation (calculated here in m7_concoag)
  ! pxxavy          = average mass of species xx in mode y []!@@@
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  ! panli(:,:,x)    = total number of particles moved by inter-modal 
  !                   coagulation from mode x [cm-3]
  ! pbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x that is moved to mode y+1 [1] 
  !                   !@@@ Clumsy notation! Should be moved to mode y !!!
  ! pso4_x          = mass of sulphate condensed on insoluble
  !                    mode x [molec. cm-3]
  !
  !--- Local variables / Constants:
  !
  ! zso4x    = available mass of sulfate from mode 1 and 2 
  !            condensing and coagulating on mode x (x = insoluble modes 5,6,7).
  !
  ! zcrtcst  = Critical constant, i.e. number of sulfate molecules to cover 
  !            an average particle of the mode with a layer of the thickness
  !            determined by cLayerThickness in mo_aero_m7. Calculated by
  !            m7_coat.
  !
  ! =>        zso4x/zcrtcst is the total number of particles that could be moved
  !            from insoluble mode x to soluble modes.
  !
  ! zcrit_x  = total available number of particles in mode x that are moved from
  !            insoluble mode x to the corresponding soluble mode.

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp)    :: pso4_5(kbdim,klev),          pso4_6(kbdim,klev),             &
             pso4_7(kbdim,klev),                                              &
             pa4av1(kbdim,klev),          pa4av2(kbdim,klev),                 &
             pbcav5(kbdim,klev),          pocav5(kbdim,klev),                 &
             pduav6(kbdim,klev),          pduav7(kbdim,klev)
              
  REAL(dp)    :: paerml(kbdim,klev,naermod),  paernl(kbdim,klev,nmod),        &
             pbfract1(kbdim,klev,nmod-1), pbfract2(kbdim,klev,nmod-1),        &
             panli(kbdim,klev,nmod),      pa4delt(kbdim,klev,naermod),        &
             pm6rp(kbdim,klev,nmod)


  ! Local variables:

  INTEGER(i4) :: jl,          jk,            jmod   

  REAL(dp)    :: zcrit_5,     zcrit_6,       zcrit_7,                         &
                 zso45,       zso46,         zso47,                           &
                 zeps

  REAL(dp)    :: zm6rp(nmod), zcrtcst(nmod)

  !--- 0) Initializations:

  zeps=EPSILON(1._dp)


  !--- 1) Redistribution of mass and numbers after nucleation, coagulation ----
  !       and coagulation calculated in the preceeding subroutines:

  DO jk=1,klev
     DO jl=1,kproma

        !--- 1.1) Sum mass of sulphate added to modes 5, 6, and 7 due to 
        !         coagulation with modes 1 and 2 (1st term) and the mass
        !         of sulfate condensed on the insoluble mode x (pso4_x):
        
        zso45=panli(jl,jk,1)*pbfract1(jl,jk,4)*pa4av1(jl,jk)+pso4_5(jl,jk)

        zso46=panli(jl,jk,1)*pbfract1(jl,jk,5)*pa4av1(jl,jk)+                &
              panli(jl,jk,2)*pbfract2(jl,jk,5)*pa4av2(jl,jk)+pso4_6(jl,jk)

        zso47=panli(jl,jk,1)*pbfract1(jl,jk,6)*pa4av1(jl,jk)+                &
              panli(jl,jk,2)*pbfract2(jl,jk,6)*pa4av2(jl,jk)+pso4_7(jl,jk)

        !--- 1.2) Determine number of particles that can be sufficiently coated
        !       by the available sulfate to be transfered to the soluble modes:

        !    Optimization of the call of m7_coat to allow for unroll and 
        !    subsequent vectorization.

!CDIR UNROLL=7
        DO jmod = 1, nmod
          zm6rp(jmod) = pm6rp(jl,jk,jmod)
        END DO

        CALL m7_coat(zm6rp,zcrtcst)

        !@@@ Changed security check to allow for inconsistent radii:

        IF(paernl(jl,jk,iaiti) >= 1.E-5_dp .AND. zcrtcst(5)>zeps) THEN
           zcrit_5=MIN(paernl(jl,jk,iaiti), zso45/zcrtcst(5))
        ELSE
           zcrit_5=0._dp
        END IF
        IF(paernl(jl,jk,iacci) >= 1.E-5_dp .AND. zcrtcst(6)>zeps) THEN
           zcrit_6=MIN(paernl(jl,jk,iacci), zso46/zcrtcst(6))
        ELSE
           zcrit_6=0._dp
        END IF
        IF(paernl(jl,jk,icoai) >= 1.E-5_dp .AND. zcrtcst(7)>zeps) THEN
           zcrit_7=MIN(paernl(jl,jk,icoai), zso47/zcrtcst(7))
        ELSE
           zcrit_7=0._dp
        END IF

        !--- 1.3) Number of particles moved from the mode 5 to 2 due to
        !         interaction with 1 and due to condensation:
        
        paernl(jl,jk,iaits)=paernl(jl,jk,iaits)+zcrit_5
        paernl(jl,jk,iaiti)=paernl(jl,jk,iaiti)-zcrit_5
        
        !--- 1.4) Mass moved from mode 5 to 2:
        
        pa4delt(jl,jk,2)=pa4delt(jl,jk,2)+pso4_5(jl,jk)
        pa4delt(jl,jk,ibcks)=pa4delt(jl,jk,ibcks)+zcrit_5*pbcav5(jl,jk)*1.e6
        pa4delt(jl,jk,iocks)=pa4delt(jl,jk,iocks)+zcrit_5*pocav5(jl,jk)*1.e6

        !--- 1.5) Mass remaining in mode 5:
        
        paerml(jl,jk,ibcki)=paerml(jl,jk,ibcki)-zcrit_5*pbcav5(jl,jk)*1.e6
        paerml(jl,jk,iocki)=paerml(jl,jk,iocki)-zcrit_5*pocav5(jl,jk)*1.e6

        !--- 1.6) Number of particles moved from the mode 6 to 3:
        
        paernl(jl,jk,iaccs)=paernl(jl,jk,iaccs)+zcrit_6
        paernl(jl,jk,iacci)=paernl(jl,jk,iacci)-zcrit_6
        
        !--- 1.7) Mass moved from mode 6 to 3:

        pa4delt(jl,jk,3)=pa4delt(jl,jk,3)+pso4_6(jl,jk)
        pa4delt(jl,jk,iduas)=pa4delt(jl,jk,iduas)+zcrit_6*pduav6(jl,jk)*1.e6

        !--- 1.8) Mass remaining in mode 6:

        paerml(jl,jk,iduai)=paerml(jl,jk,iduai)-zcrit_6*pduav6(jl,jk)*1.e6
        
        !--- 1.9) Number of particles moved from the mode 7 to 4:

        paernl(jl,jk,icoas)=paernl(jl,jk,icoas)+zcrit_7
        paernl(jl,jk,icoai)=paernl(jl,jk,icoai)-zcrit_7

        !--- 1.10) Mass moved from mode 7 to 4:
        
        pa4delt(jl,jk,4)=pa4delt(jl,jk,4)+pso4_7(jl,jk)
        pa4delt(jl,jk,iducs)=pa4delt(jl,jk,iducs)+zcrit_7*pduav7(jl,jk)*1.e6

        !--- 1.11) Mass remaining in mode 7:

        paerml(jl,jk,iduci)=paerml(jl,jk,iduci)-zcrit_7*pduav7(jl,jk)*1.e6

     END DO
  END DO

END SUBROUTINE m7_concoag
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  SUBROUTINE m7_coat(pm6rp_lon_lev, pcrtcst)
        
    ! Purpose:
    ! ---------
    ! *m7_coat* calculates the number of sulfate 
    !           molecules required to coat a particle
    !           with cLayerThickness of sulfate
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                          2001
    !
    ! Interface:
    ! ---------
    ! *m7_coat* is called from *m7_concoag*
    !

    IMPLICIT NONE

    INTEGER(i4)         :: jmod 

    REAL(dp)            :: pm6rp_lon_lev(nmod)! Ambient radii for current
                                              ! longitude and level [cm]
    REAL(dp)            :: pcrtcst(nmod)  ! Critical constant, i.e. number of
                                      ! sulfate to cover an average particle
                                      ! of the mode with a layer of the 
                                      ! thickness determined by cLayerThickness.
    REAL(dp)            :: zras(nmod)     ! Radius of average surface 
                                          ! for a single particle [cm]
    REAL(dp)            :: zas(nmod)      ! Average surface 
                                          ! for single particle [cm+2]
        
    REAL(dp), PARAMETER :: csurf_molec = 2.39E-15_dp ! Average cross-section 
                                          ! of a single H2SO4 molecule [cm+2]

    !--- 1) Calculate the radii of average surface for modes 5-7:

    zras(5) = pm6rp_lon_lev(5) * cmr2ras(5)
    zras(6) = pm6rp_lon_lev(6) * cmr2ras(6)
    zras(7) = pm6rp_lon_lev(7) * cmr2ras(7)

    DO jmod=5, 7
           
       !--- 2) Calculate the average surface of an particle for modes 5-7:

       zas(jmod)    = 4._dp * zras(jmod)**2 * pi

       !--- 3) Determine the number of sulfate molecules needed to form
       !       a cLayerThickness thick layer of sulfate on the particles
       !       in modes 5-7:
           
       pcrtcst(jmod) = (zas(jmod) / csurf_molec) * cLayerThickness

    END DO

  END SUBROUTINE m7_coat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
END SUBROUTINE m7_dnum
!------------------------------------------------------------------------------!
SUBROUTINE m7_dconc (kproma, kbdim, klev, paerml, paernl, pm6dry)
  !
  !    *m7_dconc*  changes aerosol numbers and masses to account for
  !                condensational growth of the mode mean radii
  !  
  !    Authors:
  !    --------
  !    J. Wilson and E. Vignati, JRC (original source)            May 2000
  !    P. Stier, MPI-MET (f90 version, changes, comments)             2001
  !
  !    Purpose:
  !    --------
  !    This routine repartitions aerosol number and mass between the
  !    the modes to account for condensational growth and the formation
  !    of an accumulation mode from the upper tail of the aitken mode.
  !
  !    Interface:
  !    ----------
  !    *m7_dconc* is called from *m7*
  !
  !    Method:
  !    -------
  !    The routine calculates the cumulativ number and mass distribution of the
  !    modes up to the respective mode boundary:
  !
  !                        / x                              _
  !                 N      |       1           1   ln(R)-ln(R)  2
  !    N(0,x) = ---------  |   --------  exp(- - ( ----------- )   ) d ln(R) 
  !             ln(sigma)  |   sqrt(2PI)       2    ln(sigma)
  !                        / 0 
  !                         
  !                         /tx                   2
  !                        |        1            t
  !           =     N      |     --------  exp(- - ) d t 
  !                        |     sqrt(2PI)       2 
  !                        /-inf
  ! 
  !    where:                   
  !
  !                        _
  !               ln(R)-ln(R)
  !    t      =   -----------
  !                ln(sigma)
  !
  !    and:
  !                        _
  !               ln(x)-ln(R)
  !    tx     =   -----------
  !                ln(sigma)
  !    _
  !    R is the Count Mean Radius or the Mass Mean Radius.
  !
  !    Practically, the routine m7_cumnor calculates the fraction 
  !    of the number and
  !    mass distribution for each mode lying below the respective upper mode 
  !    boundary (1).
  !    In a next step the net fraction of each mode lying between the upper and
  !    lower
  !    mode boundaries are summed up (2) and the numbers and masses exceeding 
  !    the mode
  !    boundaries are transfered to the neighboring larger mode (3). 
  !    Finally, these quantities are stored in the respective arrays
  !    paernl and paerml (4).
  !    The repartititioning is currently only done for the soluble modes as it
  !     is
  !    assumed that insoluble modes are rather transfered to the soluble modes 
  !    and grow as soluble particles.
  !
  !    Externals:
  !    ----------
  !    None
  !
  !--- Parameter list:
  !
  !    paerml(kbdim,klev,naermod) = total aerosol mass for each compound 
  !                              [molec. cm-3 for sulfate and ug m-3 for others]
  !    paernl(kbdim,klev,nmod) = aerosol number for each mode [cm-3]  !
  !    sigma(jmod)            = standard deviation of mode jmod [1]
  !    crdiv                  = threshold radii between the different modes [cm]
  !                        crdiv(jmod) is the lower bound and crdiv(jmod+1) is 
  !                                the upper bound of the respective mode
  !
  !--- Local Variables:
  !
  !    zfconn(:,:,jnum,jmod)  = absolute fraction of the number of particles
  !                       in mode jmod, 
  !                       with CMD=2*pm6dry(jmod) and a geometric standard 
  !                       deviation zrcsig, that are smaller than crdiv(jnum+1).
  !                         I.e. with 0 < CMD < crdiv(jnum+1) [1]
  !    zfconm(:,:,jnum,jmod) = absolute fraction of the mass in mode jmod, 
  !                         with CMD=2*pm6dry(jmod) and a geometric standard 
  !                         deviation zrcsig, that are smaller than crdiv(jnum+1).
  !                         I.e. with 0 < CMD < crdiv(jnum+1) [1]
  !

    IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL(dp):: paerml(kbdim,klev,naermod),   paernl(kbdim,klev,nmod),   &
             pm6dry(kbdim,klev,nsol)

  ! Local variables:
  !
  INTEGER :: jmod,          jnum,          jl,           jk

  REAL(dp):: zrcsig,     &!zarg1, 
             zarg2,      zdpcm,   &
             zarg3,      zdpam,         zcongn,        zdummy,  &
             zr1,        zr2,           zttnj,        zavnj,         zmrj,    &
             zmt,        znt,           zavmt,        zmcr,          zfconmj, &
             zntnew,     zmtnew,        zdm,          zeps
  
  REAL(dp)    :: zambc2(4),     zambc3(4),     zambc4(4),          &
                 zamoc2(4),     zamoc3(4),     zamoc4(4),          &
                                zamss3(4),     zamss4(4),          &  
                                zamdu3(4),     zamdu4(4)

  REAL(dp):: ztotmass(kbdim,klev)

  REAL(dp)    :: zsumn(kbdim,klev,nmod),       zsumm(kbdim,klev,nmod),         &
                 zsumbc(kbdim,klev,3),         zsumoc(kbdim,klev,3),           &
                 zsumss(kbdim,klev,2),         zsumdu(kbdim,klev,2)
 
  REAL(dp)    :: zfconn(kbdim,klev,nsol,nsol), zfconm(kbdim,klev,nsol,nsol)


  !--- 0) Initialisations: -----------------------------------------------------
  zeps=EPSILON(1.0_dp) 

  zsumn(:,:,:)    = 0._dp
  zsumm(:,:,:)    = 0._dp
  zsumbc(:,:,:)   = 0._dp
  zsumoc(:,:,:)   = 0._dp
  zsumss(:,:,:)   = 0._dp
  zsumdu(:,:,:)   = 0._dp

  zfconm(:,:,:,:) = 0._dp
  zfconn(:,:,:,:) = 0._dp

  !
  !--- 1) Identify how much the mode jmod has grown into the next higher mode --
  !
  DO jmod=1,nsol-1

     !--- Total mass of the mode in equivalent molecules of sulfate:

     SELECT CASE(jmod)
     CASE(1)
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4ns)
     CASE(2)
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4ks) + &
                               (paerml(1:kproma,:,ibcks)/dbc+paerml(1:kproma,:,iocks)/doc)  &
                               *dh2so4/wh2so4*avo*1.E-12
     CASE(3)
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4as) + &
                               (paerml(1:kproma,:,ibcas)/dbc+paerml(1:kproma,:,iocas)/doc+  &
                                paerml(1:kproma,:,issas)/dnacl+paerml(1:kproma,:,iduas)/ddust) &
                               *dh2so4/wh2so4*avo*1.E-12
     END SELECT

     DO jnum=jmod,nsol-1
        DO jk=1,klev
           DO jl=1,kproma
              IF (paernl(jl,jk,jmod) .GT. zeps .AND. & ! mz_ak_20041014
                       pm6dry(jl,jk,jmod) .GT. 0.0_dp) THEN

                 !--- 1.1) Calculate necessary parameters:                 

                 !--- Geometric Standard Deviation:
                 zrcsig=LOG(sigma(jmod))

                 !--- Mass Median Radius:
!                 zarg1=pm6dry(jl,jk,jmod)*cmedr2mmedr(jmod)

                 !--- Count Median Radius:
                 zarg2=pm6dry(jl,jk,jmod)

                 !--- Threshold radius between the modes:
                 zarg3=crdiv(jnum+1)

                 !--- Transfer to logarithmic scale:
!                 zdpmm=LOG(zarg1)
                 zdpcm=LOG(zarg2)
                 zdpam=LOG(zarg3)

                 !--- Distance of the CMD of the mode from the threshold mode  
                 !    diameter in terms of geometric standard deviations:

                 zcongn=(zdpam-zdpcm)/zrcsig

                 !--- Distance of the MMD of the mode from the threshold mode  
                 !    diameter in terms of geometric standard deviations (t):

!                 zcongm=(zdpam-zdpmm)/zrcsig

                 !--- Calculate the cumulative of the log-normal number 
                 !    distribution:

                 CALL m7_cumnor(zcongn,zfconn(jl,jk,jnum,jmod),zdummy)

                 ! mz_ak_20041014+
                 !--- Limit transfer only to adjacent modes:

                 IF (jnum .GT. jmod) THEN
                    zfconn(jl,jk,jnum,jmod)= 1.0_dp
                 END IF

                 !--- Set minimum radius and maximum radius:

                 zr1 = crdiv(jmod)
                 zr2 = crdiv(jmod+1)

                 !--- Radius of average mass for a lognormal distribution

                 zdm = EXP((LOG(zr1)+LOG(zr2))/2.0_dp)*cmr2ram(jmod)

                 !--- Average mass contained in the mode

                 zttnj = ztotmass(jl,jk)/paernl(jl,jk,jmod)

                 !--- Average number of sulfate molecules or equivalent for mixed modes,
                 !    for a particle with radius zdm

                 zavnj=zdm**3.0*pi*avo*dh2so4/wh2so4/0.75_dp

                 !- If the average mass contained in the mode is larger than
                 !  the average mass for a particle with radius zdm, the ! 
                 !  transfer of number and mass is done,
                 !   else there is no transfer

                 IF (zttnj .GT. zavnj .AND. jnum .EQ. jmod) THEN

                    !--- Mass remaining in the mode

                    zmrj=zfconn(jl,jk,jnum,jmod)*paernl(jl,jk,jmod)*zavnj

                    !--- Mass transferred

                    zmt=ztotmass(jl,jk)-zmrj

                    !--- Numbers transferred 

                    znt=(1.0_dp-zfconn(jl,jk,jnum,jmod))*paernl(jl,jk,jmod)

                    !--- Average mass of particles transferred

                    IF(znt>zeps) THEN
                       zavmt=zmt/znt
                    ELSE
                       zavmt=0.0_dp
                    END IF

                    !--- Average mass of particles of radius zr2

                    zmcr=(zr2*cmr2ram(jmod))**3.0*pi*avo*dh2so4/wh2so4/0.75_dp

                    !- If the average mass of particle transferred is smaller
                    !  than the average mass of particles with radius zr2 then
                    !  reduce the particles transferred so that zavmt=zmcr, else
                    !  calculate the mass fraction transferred zfconmj

                    IF (zavmt .GE. zmcr) THEN
                       zfconmj=zmrj/ztotmass(jl,jk)
                    ELSE
                       zntnew = znt/(1.0_dp + (zmcr-zavmt)/(zavmt-zavnj))
                       zmtnew = zntnew*zmcr
                       zfconmj = 1.0_dp - zmtnew/ztotmass(jl,jk)
                       zfconn(jl,jk,jnum,jmod) = 1.0_dp - &
                            zntnew/paernl(jl,jk,jmod)
                    END IF
                    zfconm(jl,jk,jnum,jmod)=zfconmj
                 ELSE
                    zfconn(jl,jk,jnum,jmod)=1._dp
                    zfconm(jl,jk,jnum,jmod)=1._dp
                 END IF
              ELSE
                 zfconn(jl,jk,jnum,jmod)=1._dp
                 zfconm(jl,jk,jnum,jmod)=1._dp
              END IF
           END DO
        END DO
     END DO
  END DO

  DO jmod=1,nsol ! mz_ak_20041014
     DO jk=1,klev
        DO jl=1,kproma

!        DO jmod=1,nsol
           
           !--- 2) Calculate the net fraction of mode jmod that is transfered -
           !       to the mode jnew zfconn(:,:,jnew,jmod) :
           
          !--- Numbers:
           zfconn(jl,jk,4,jmod)=1.0_dp              -zfconn(jl,jk,3,jmod)
           zfconn(jl,jk,3,jmod)=zfconn(jl,jk,3,jmod)-zfconn(jl,jk,2,jmod)
           zfconn(jl,jk,2,jmod)=zfconn(jl,jk,2,jmod)-zfconn(jl,jk,1,jmod)

           !--- Mass:
           zfconm(jl,jk,4,jmod)=1.0_dp             -zfconm(jl,jk,3,jmod)
           zfconm(jl,jk,3,jmod)=zfconm(jl,jk,3,jmod)-zfconm(jl,jk,2,jmod)
           zfconm(jl,jk,2,jmod)=zfconm(jl,jk,2,jmod)-zfconm(jl,jk,1,jmod)

           !--- 3) Sum the net masses and numbers transfered between the modes: 

           !--- 3.1) Soluble mode numbers and sulfate mass:

           zsumn(jl,jk,1)=zsumn(jl,jk,1)+paernl(jl,jk,jmod)*zfconn(jl,jk,1,jmod)
           zsumm(jl,jk,1)=zsumm(jl,jk,1)+paerml(jl,jk,jmod)*zfconm(jl,jk,1,jmod)
           zsumn(jl,jk,2)=zsumn(jl,jk,2)+paernl(jl,jk,jmod)*zfconn(jl,jk,2,jmod)
           zsumm(jl,jk,2)=zsumm(jl,jk,2)+paerml(jl,jk,jmod)*zfconm(jl,jk,2,jmod)
           zsumn(jl,jk,3)=zsumn(jl,jk,3)+paernl(jl,jk,jmod)*zfconn(jl,jk,3,jmod)
           zsumm(jl,jk,3)=zsumm(jl,jk,3)+paerml(jl,jk,jmod)*zfconm(jl,jk,3,jmod)
           zsumn(jl,jk,4)=zsumn(jl,jk,4)+paernl(jl,jk,jmod)*zfconn(jl,jk,4,jmod)
           zsumm(jl,jk,4)=zsumm(jl,jk,4)+paerml(jl,jk,jmod)*zfconm(jl,jk,4,jmod)

        END DO
        ! mz_ak_20041014+
     END DO
  END DO
  
  DO jk=1,klev
     DO jl=1,kproma
        
        ! mz_ak_20041014-

        !--- 3.2) Non-sulfate masses:

        zambc2(2)=paerml(jl,jk,ibcks)*zfconm(jl,jk,2,2)
        zambc2(3)=paerml(jl,jk,ibcks)*zfconm(jl,jk,3,2)
        zambc2(4)=paerml(jl,jk,ibcks)*zfconm(jl,jk,4,2)
        zamoc2(2)=paerml(jl,jk,iocks)*zfconm(jl,jk,2,2)
        zamoc2(3)=paerml(jl,jk,iocks)*zfconm(jl,jk,3,2)
        zamoc2(4)=paerml(jl,jk,iocks)*zfconm(jl,jk,4,2)
        zambc3(3)=paerml(jl,jk,ibcas)*zfconm(jl,jk,3,3)
        zambc3(4)=paerml(jl,jk,ibcas)*zfconm(jl,jk,4,3)
        zamoc3(3)=paerml(jl,jk,iocas)*zfconm(jl,jk,3,3)
        zamoc3(4)=paerml(jl,jk,iocas)*zfconm(jl,jk,4,3)
        zambc4(4)=paerml(jl,jk,ibccs)
        zamoc4(4)=paerml(jl,jk,ioccs)
        zamss3(3)=paerml(jl,jk,issas)*zfconm(jl,jk,3,3)
        zamss3(4)=paerml(jl,jk,issas)*zfconm(jl,jk,4,3)
        zamdu3(3)=paerml(jl,jk,iduas)*zfconm(jl,jk,3,3)
        zamdu3(4)=paerml(jl,jk,iduas)*zfconm(jl,jk,4,3)
        zamss4(4)=paerml(jl,jk,isscs)
        zamdu4(4)=paerml(jl,jk,iducs)
        
        zsumbc(jl,jk,1)=zambc2(2)
        zsumbc(jl,jk,2)=zambc2(3)+zambc3(3)
        zsumbc(jl,jk,3)=zambc2(4)+zambc3(4)+zambc4(4)
        zsumoc(jl,jk,1)=zamoc2(2)
        zsumoc(jl,jk,2)=zamoc2(3)+zamoc3(3)
        zsumoc(jl,jk,3)=zamoc2(4)+zamoc3(4)+zamoc4(4)
        zsumss(jl,jk,1)=zamss3(3)
        zsumss(jl,jk,2)=zamss3(4)+zamss4(4)
        zsumdu(jl,jk,1)=zamdu3(3)
        zsumdu(jl,jk,2)=zamdu3(4)+zamdu4(4)

     END DO
  END DO

  !--- 4) Store final masses and numbers of the modes: ------------------------

  DO jmod=1,nsol
     DO jk=1,klev
        DO jl=1,kproma

           !--- Particle numbers:
           
           paernl(jl,jk,jmod)=zsumn(jl,jk,jmod)

           !--- Sulfate mass:

           paerml(jl,jk,jmod)=zsumm(jl,jk,jmod)

        END DO
     END DO
  END DO

  !--- Non sulfate masses:

  DO jk=1,klev
     DO jl=1,kproma

        paerml(jl,jk,ibcks)=zsumbc(jl,jk,1)
        paerml(jl,jk,ibcas)=zsumbc(jl,jk,2)
        paerml(jl,jk,ibccs)=zsumbc(jl,jk,3)
        paerml(jl,jk,iocks)=zsumoc(jl,jk,1)
        paerml(jl,jk,iocas)=zsumoc(jl,jk,2)
        paerml(jl,jk,ioccs)=zsumoc(jl,jk,3)
        paerml(jl,jk,issas)=zsumss(jl,jk,1)
        paerml(jl,jk,isscs)=zsumss(jl,jk,2)
        paerml(jl,jk,iduas)=zsumdu(jl,jk,1)
        paerml(jl,jk,iducs)=zsumdu(jl,jk,2)

     END DO
  END DO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
CONTAINS
SUBROUTINE m7_cumnor ( arg, RESULT, ccum )
  !
  !****************************************************************************
  !
  !! CUMNOR computes the cumulative normal distribution.
  !
  !
  !     the integral from -infinity to x of
  !          (1/sqrt(2*pi)) exp(-u*u/2) du
  !
  !  Author:
  !  -------
  !  Original source:
  !
  !    W. J. Cody    Mathematics and Computer Science Division
  !                  Argonne National Laboratory
  !                  Argonne, IL 60439
  !
  !    DCDFLIB is attributed to Barry Brown, James Lovato, and Kathy Russell
  !            bwb@odin.mda.uth.tmc.edu.
  !
  !    Adopted to ECHAM/M7:
  !
  !    Philip Stier  (MPI-MET)                    2001
  !
  !
  !  Reference:
  !  ----------
  !
  !    W D Cody, 
  !    "ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special 
  !    Function Routines and Test Drivers"
  !    ACM Transactions on Mathematical Software,
  !    Volume 19, 1993, pages 22-32.
  !
  !  Parameters:
  !
  !     ARG --> Upper limit of integration.
  !                                        X is REAL(dp)
  !
  !     RESULT <-- Cumulative normal distribution.
  !                                        RESULT is REAL(dp)
  !
  !     CCUM <-- Complement of Cumulative normal distribution.
  !                                        CCUM is REAL(dp)
  !
  !
  ! Original Comments:
  !
  !
  ! This function evaluates the normal distribution function:
  !
  !                              / x
  !                     1       |       -t*t/2
  !          P(x) = ----------- |      e       dt
  !                 sqrt(2 pi)  |
  !                             /-oo
  !
  !   The main computation evaluates near-minimax approximations
  !   derived from those in "Rational Chebyshev approximations for
  !   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
  !   This transportable program uses rational functions that
  !   theoretically approximate the normal distribution function to
  !   at least 18 significant decimal digits.  The accuracy achieved
  !   depends on the arithmetic system, the compiler, the intrinsic
  !   functions, and proper selection of the machine-dependent
  !   constants.
  !
  !  Explanation of machine-dependent constants.
  !
  !   MIN   = smallest machine representable number.
  !
  !   EPS   = argument below which anorm(x) may be represented by
  !           0.5  and above which  x*x  will not underflow.
  !           A conservative value is the largest machine number X
  !           such that   1.0 + X = 1.0   to machine precision.
  !
  !  Error returns
  !
  !  The program returns  ANORM = 0     for  ARG .LE. XLOW.
  !
  !  Author: 
  !
  !    W. J. Cody
  !    Mathematics and Computer Science Division
  !    Argonne National Laboratory
  !    Argonne, IL 60439
  !
  !  Latest modification: March 15, 1992
  !
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: a = (/ &
       2.2352520354606839287d00, &
       1.6102823106855587881d02, &
       1.0676894854603709582d03, &
       1.8154981253343561249d04, &
       6.5682337918207449113d-2 /)
  REAL(dp) arg
  REAL(dp), PARAMETER, DIMENSION ( 4 ) :: b = (/ &
       4.7202581904688241870d01, &
       9.7609855173777669322d02, &
       1.0260932208618978205d04, &
       4.5507789335026729956d04 /)
  REAL(dp), PARAMETER, DIMENSION ( 9 ) :: c = (/ &
       3.9894151208813466764d-1, &
       8.8831497943883759412d00, &
       9.3506656132177855979d01, &
       5.9727027639480026226d02, &
       2.4945375852903726711d03, &
       6.8481904505362823326d03, &
       1.1602651437647350124d04, &
       9.8427148383839780218d03, &
       1.0765576773720192317d-8 /)
  REAL(dp) ccum
  REAL(dp), PARAMETER, DIMENSION ( 8 ) :: d = (/ &
       2.2266688044328115691d01, &
       2.3538790178262499861d02, &
       1.5193775994075548050d03, &
       6.4855582982667607550d03, &
       1.8615571640885098091d04, &
       3.4900952721145977266d04, &
       3.8912003286093271411d04, &
       1.9685429676859990727d04 /)
  REAL(dp) del
!@@@ REAL(dp) dpmpar
  REAL(dp) eps
  INTEGER i
  REAL(dp) min
  REAL(dp), PARAMETER, DIMENSION ( 6 ) :: p = (/ &
       2.1589853405795699d-1, &
       1.274011611602473639d-1, &
       2.2235277870649807d-2, &
       1.421619193227893466d-3, &
       2.9112874951168792d-5, &
       2.307344176494017303d-2 /)
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: q = (/ &
       1.28426009614491121d00, &
       4.68238212480865118d-1, &
       6.59881378689285515d-2, &
       3.78239633202758244d-3, &
       7.29751555083966205d-5 /)
  REAL(dp) RESULT
  REAL(DP), PARAMETER :: root32 = 5.656854248E0_dp
  REAL(DP), PARAMETER :: sixten = 16.0_dp
  REAL(DP) temp
  REAL(DP), PARAMETER :: sqrpi = 3.9894228040143267794E-1_dp
  REAL(DP), PARAMETER :: thrsh = 0.66291E0_dp
  REAL(DP) x
  REAL(DP) xden
  REAL(DP) xnum
  REAL(DP) y
  REAL(DP) xsq
  !
  !  Machine dependent constants
  !
  eps = EPSILON ( 1.0E0_dp ) * 0.5E0_dp
  !
  !@@@ Simplified calculation of the smallest machine representable number
  !    (Higher accuracy than needed!)
  !
  !@@@ min = dpmpar(2)

  min = epsilon ( 1.0E0_dp)

  x = arg
  y = ABS ( x )

  IF ( y <= thrsh ) THEN
     !
     !  Evaluate  anorm  for  |X| <= 0.66291
     !
     IF ( y > eps ) THEN
        xsq = x * x
     ELSE
        xsq = 0.0_dp
     END IF

     xnum = a(5) * xsq
     xden = xsq
     DO i = 1, 3
        xnum = ( xnum + a(i) ) * xsq
        xden = ( xden + b(i) ) * xsq
     END DO
     RESULT = x * ( xnum + a(4) ) / ( xden + b(4) )
     temp = RESULT
     RESULT = 0.5_dp + temp
     ccum = 0.5_dp - temp
     !
     !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
     !
  ELSE IF ( y <= root32 ) THEN

     xnum = c(9) * y
     xden = y
!CDIR UNROLL=7
     DO i = 1, 7
        xnum = ( xnum + c(i) ) * y
        xden = ( xden + d(i) ) * y
     END DO
     RESULT = ( xnum + c(8) ) / ( xden + d(8) )
     xsq = AINT ( y * sixten ) / sixten
     del = ( y - xsq ) * ( y + xsq )
     RESULT = EXP(-xsq*xsq*0.5) * EXP(-del*0.5) * RESULT
     ccum = 1.0_dp - RESULT

     IF ( x > 0.0_dp  ) THEN
        temp = RESULT
        RESULT = ccum
        ccum = temp
     END IF
     !
     !  Evaluate  anorm  for |X| > sqrt(32).
     !
  ELSE

     RESULT = 0.0_dp
     xsq = 1.0 / ( x * x )
     xnum = p(6) * xsq
     xden = xsq
     DO i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
     END DO

     RESULT = xsq * ( xnum + p(5) ) / ( xden + q(5) )
     RESULT = ( sqrpi - RESULT ) / y
     xsq = AINT ( x * sixten ) / sixten
     del = ( x - xsq ) * ( x + xsq )
     RESULT = EXP ( - xsq * xsq * 0.5 ) * EXP ( - del * 0.5 ) * RESULT
     ccum = 1.0_dp - RESULT

     IF ( x > 0.0_dp ) THEN
        temp = RESULT
        RESULT = ccum
        ccum = temp
     END IF

  END IF

  IF ( RESULT < min ) THEN
     RESULT = 0.0_dp
  END IF

  IF ( ccum < min ) THEN
     ccum = 0.0_dp
  END IF

END SUBROUTINE m7_cumnor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
END SUBROUTINE m7_dconc
!-----------------------------------------------------------------------------!
  SUBROUTINE nucl_kulmala(kproma,  kbdim,  klev,    &  ! ECHAM5 dims
                          pso4g,   ptp1,   prelhum, & 
                          pbnrate, palpha, pbeta    )

    
    !  Authors:
    !  --------
    !  P. Stier, MPI-Met, Hamburg,    from the original f77 code
    !                                 in the routine m7_nuck        2001-2003
    !  J. Wilson, E. Vignati, JRC/EI, original source                 09/2000
    !
    !  Purpose:                                                           
    !  --------                                                           
    !  This routine calculates the instananeous nucleation rate            
    !  znucrate [molec. cm-3 s-1] from a given gas phase H2SO4 concentration    
    !  pso4g [molec. cm-3]. It also calculates the integrated change of 
    !  H2SO4 gas phase mass over one timestep due to nucleation 
    !  pa4delt(:,:,1) [molec. cm-3] as well as the number of nucleated 
    !  particles panew [1] during the timestep.
    !
    !  Interface:
    !  ----------
    !  *nucl_kulmala* is called from *m7_nuck*
    !
    !  Method:
    !  -------
    !  Kulmala et al. (1998) 's formula for binary nucleation is 
    !  rewritten to take the form znucrate = exp[zalpha+ln(pso4g)*beta]. 
    !  Equation numbers are taken from Kulmala et al. (1998).
    !  After the calculation of the nucleation rate znucrate, it is 
    !  integrated in 2) analytically over one timestep, i.e.:
    !
    !  Integration of:
    ! 
    !  znucrate=d(critn*znav)/dt=exp[zalpha + zbeta*ln(znav)]
    !
    !  gives znav(t0+dt, znav(t0)) where znav(t0) is pso4g.
    !  znav is temporarily stored in zso4g_new and confined between
    !  0 and pso4g. 
    !  The number of nucleated particles is then calculated by 
    !  assuming a fixed critical mass critn of newly nucleated 
    !  particles and dividing the total mass of nucleated sulfate 
    !  by this critical mass of one nucleated particle. 
    !
    !  Externals:
    !  ----------
    !  None

    IMPLICIT NONE

    INTEGER     :: kproma,  kbdim, klev

    REAL(dp)    :: pso4g(kbdim,klev),         ptp1(kbdim,klev),      &
                   prelhum(kbdim,klev),       pbnrate(kbdim,klev),   &
                   palpha(kbdim,klev),        pbeta(kbdim,klev)

    INTEGER(i4) :: jl, jk

    REAL(dp)    :: znwv,        zln_nac,      ztk,         zsupsat,  &
                   zpeh2o,      zpeh2so4,     zra,         zxal,     &
                   ztkn,        zssn,         zdelta

    !---1) Calculation of the nucleation rate: ----------------------------
    
    DO jk=1,klev
       DO jl=1,kproma

          IF (pso4g(jl,jk) .GT. 1e-5_dp) THEN
             ztk=ptp1(jl,jk)
             zsupsat=prelhum(jl,jk)
             !
             !- 1.1) Restrict t, and rh to limits where the parameterization ok:
             !
             ztkn = MAX(ztk, 220.0_dp)
             zssn = MIN(zsupsat, 0.90_dp)

             !
             !--- 1.2) Equlibrium vapour pressures (Jaeker-Mirabel (1995), JGR):
             !
             !--- H2O equlibrium vapour pressure (Tabata):
             !
             zpeh2o=0.750064*(10.**(8.42926609_dp-1827.17843/          &
                  ztkn-71208.271/ztkn/ztkn))*1333./bk/ztkn
             !
             !--- H2SO4 equlibrium vapour pressure at 360 
             !@@@ Check source: which Ayers?
             !
             zpeh2so4=EXP(-10156./ztkn+16.259_dp)*7.6e2*1333./bk/ztkn
             !
             !--- H2SO4 equlibrium vapour pressure - correction of ayers
             !    by kulmala - currently not used
             !
             !     payers=exp(-10156/360+16.259)*7.6e2
             !     zpeh2so4=exp(log(payers)+10156*(-1./ztkn+1./360.+0.38/
             !              (905-360) * 
             !              (1+log(360./ztkn)-360./ztkn)))*1333/bk/ztkn
             !
             !--- 1.3) Relative acidity (0.0 -1.0):
             ! 
             zra=pso4g(jl,jk)/zpeh2so4
             !
             !--- 1.4) Water vapour molecule concentration [cm-3]:
             !
             znwv=zsupsat*zpeh2o
             !
             !--- 1.5) Factor delta in Eq. 22:
             ! 
             zdelta=1.0_dp+(ztkn-273.15_dp)/273.15
             !
             !--- 1.6) Molefraction of H2SO4 in the critical cluster 
             !         minus the H2SO4(g) term in Eq. 17:
             !
             zxal =1.2233_dp-0.0154*zra/(zra+zssn)-0.0415*LOG(znwv)+ 0.0016*ztkn
             !
             !--- 1.7) Exponent of the critical cluster (Eq. 18):
             !
             zln_nac = -14.5125_dp+0.1335*ztkn-10.5462_dp*zssn &
                  +1958.4_dp*zssn/ztkn
             !
             !--- 1.8) Sum of all terms in Eq. 20 containing H2SO4(g):
             !
             pbeta(jl,jk) = 25.1289_dp - 4890.8_dp/ztkn &
                  + 7643.4_dp*0.0102_dp/ztkn - &
                  2.2479_dp*zdelta*zssn - 1.9712_dp*0.0102_dp*zdelta/zssn
             ! 
             !--- 1.9) Sum all terms in Eq. 20 not containing H2SO4(g):
             !
             palpha(jl,jk) = zln_nac*(-25.1289_dp + 4890.8_dp/ztkn + &
                             2.2479_dp*zdelta*zssn) - &
                             1743.3_dp/ztkn + zxal*(7643.4_dp/ztkn - &
                             1.9712_dp*zdelta/zssn)
             !
             !--- 1.10) Nucleation rate [cm-3 s-1] (Kulmala et al., 1998):
             !
             pbnrate(jl,jk) = EXP(palpha(jl,jk)+LOG(pso4g(jl,jk))*pbeta(jl,jk))

          ELSE

             palpha(jl,jk) =0._dp 
             pbeta(jl,jk)  =0._dp
             pbnrate(jl,jk)=0._dp

          END IF ! pso4g(jl,jk) .GT. 1e-5
       END DO ! kproma
    END DO !klev


  END SUBROUTINE nucl_kulmala

!===========================================================================!

  SUBROUTINE nucl_vehkamaeki(kproma,   kbdim, klev,       &  ! ECHAM5 dimensions
                             ptp1,     prhd,  pmolecH2SO4, & ! ECHAM5 
                                  ! temperature, relative humidity
                             pxtrnucr, pntot               )  
                                  ! nucleation rate, number of molecules in the
                                  ! critical cluster
                         
    !
    !   Authors:
    !   ---------
    !   C. TIMMRECK, MPI HAMBURG                                           2002
    !
    !   Purpose
    !   ---------
    !   Calculation of classical nucleation rate
    !               
    !   calculation of the nucleation rate after Vehkamaeki et al. (2002)
    !   The calculation of the nucrate ZKNH2SO4 is in cm^-3 s^-1
    !   and a coarse approxmation for the first class
    !
    !   Modifications:
    !   --------------
    !   R. Hommel; rewrite in f90, adopted to ECHAM5; MPI HAMBURG;    Dec. 2002
    !   P. Stier; bugfixes, modularisation and optimization; MPI HAMBURG; 
    !                                                                 2003-2004
    !
    !
    !   H2SO4 still fixed to xxx molc/cm3, no sulfur cycle coupling yet
    !
    !   References:
    !   -----------
    !   Vehkamaeki et al. (2002), An improved parameterization for sulfuric
    !      acid/water nucleation rates for tropospheric and stratospheric
    !      conditions, J. Geophys. Res, 107, D22, 4622
    !
    !   Parameters
    !   ----------
    !   prho = prhop_neu in *sam*
    !
    !   prhd = relative humidity in sam_aeroprop & sam_nucl
    !
    !   pxtrnucr = nucleation rate in [1/m3s]
    !   xrhoc    = density of the critical nucleus in kg/m^3
    !
    !----------------------------------------------------
    
    IMPLICIT NONE

    INTEGER     :: kproma, kbdim, klev
    INTEGER(i4) :: jk, jl
    
    !----------------------------------------------------
    !
    
    REAL(dp)    ::   ptp1(kbdim,klev), prhd(kbdim,klev), &
                     pxtrnucr(kbdim,klev),  &
                     pmolecH2SO4(kbdim,klev), &            ! revisited, ok
                     pntot(kbdim,klev)
    
    !----------------------------------------------------  
    ! Local Arrays
    
    REAL(dp)::   zrhoa, zrh, zt, x, zjnuc, &! zrc, &
                 zntot, zxmole

    !--- 0) Initializations:
    
    DO jk=1, klev
       DO jl=1,kproma
          
  !----1.) Parameterization of  nucleation rate after Vehkamaeki et al. (2002) 
        ! t: temperature in K (190.15-300.15K)                              
        ! zrh: saturatio ratio of water (0.0001-1)   
        ! zrhoa: sulfuric acid concentration in 1/cm3 (10^4-10^11 1/cm3) 
        ! jnuc: nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)        
        ! ntot: total number of molecules in the critical cluster (ntot>4) 
        ! x: molefraction of H2SO4 in the critical cluster                 
        ! rc: radius of the critical cluster in nm                    

        ! Calculate nucleation only for valid thermodynamic conditions:

        IF( (pmolecH2SO4(jl,jk)>=1.e+4_dp)                      .AND. &
            (prhd(jl,jk) >=1.e-4_dp)                            .AND. &
            (ptp1(jl,jk)>=190.15_dp .AND. ptp1(jl,jk)<=300.15_dp)  ) THEN 

        zrhoa=MIN(pmolecH2SO4(jl,jk),1.e11_dp)
        zrh=MIN(prhd(jl,jk),1.0_dp)
        zt=ptp1(jl,jk)
          
        ! Equation (11) - molefraction of H2SO4 in the critical cluster

        x=0.7409967177282139_dp - 0.002663785665140117*zt   &
          + 0.002010478847383187*LOG(zrh)    &
          - 0.0001832894131464668*zt*LOG(zrh)    &
          + 0.001574072538464286*LOG(zrh)**2        &  
          - 0.00001790589121766952*zt*LOG(zrh)**2    &
          + 0.0001844027436573778*LOG(zrh)**3     &
          -  1.503452308794887e-6*zt*LOG(zrh)**3    &
          - 0.003499978417957668*LOG(zrhoa)   &
          + 0.0000504021689382576*zt*LOG(zrhoa)   

        zxmole=x !qqq

        ! Equation (12) - nucleation rate in 1/cm3s
          
        zjnuc =0.1430901615568665_dp + 2.219563673425199*zt -   &
              0.02739106114964264*zt**2 +     &
              0.00007228107239317088*zt**3 + 5.91822263375044/x +     &
              0.1174886643003278*LOG(zrh) + 0.4625315047693772*zt*LOG(zrh) -   &
              0.01180591129059253*zt**2*LOG(zrh) +     &
              0.0000404196487152575*zt**3*LOG(zrh) +    &
              (15.79628615047088*LOG(zrh))/x -     &
              0.215553951893509*LOG(zrh)**2 -    &
              0.0810269192332194*zt*LOG(zrh)**2 +     &
              0.001435808434184642*zt**2*LOG(zrh)**2 -    & 
              4.775796947178588e-6*zt**3*LOG(zrh)**2 -     &
              (2.912974063702185*LOG(zrh)**2)/x -   &
              3.588557942822751*LOG(zrh)**3 +     &
              0.04950795302831703*zt*LOG(zrh)**3 -     &
              0.0002138195118737068*zt**2*LOG(zrh)**3 +    & 
              3.108005107949533e-7*zt**3*LOG(zrh)**3 -     &
              (0.02933332747098296*LOG(zrh)**3)/x +     &
              1.145983818561277*LOG(zrhoa) -    &
              0.6007956227856778*zt*LOG(zrhoa) +    &
              0.00864244733283759*zt**2*LOG(zrhoa) -    &
              0.00002289467254710888*zt**3*LOG(zrhoa)! -    &

        zjnuc =zjnuc - &
             (8.44984513869014*LOG(zrhoa))/x +    &
              2.158548369286559*LOG(zrh)*LOG(zrhoa) +   & 
              0.0808121412840917*zt*LOG(zrh)*LOG(zrhoa) -    &
              0.0004073815255395214*zt**2*LOG(zrh)*LOG(zrhoa) - &
              4.019572560156515e-7*zt**3*LOG(zrh)*LOG(zrhoa) +    &
              (0.7213255852557236*LOG(zrh)*LOG(zrhoa))/x +    &
              1.62409850488771*LOG(zrh)**2*LOG(zrhoa) -    &
              0.01601062035325362*zt*LOG(zrh)**2*LOG(zrhoa) +   & 
              0.00003771238979714162*zt**2*LOG(zrh)**2*LOG(zrhoa) +    &
              3.217942606371182e-8*zt**3*LOG(zrh)**2*LOG(zrhoa) -    &
              (0.01132550810022116*LOG(zrh)**2*LOG(zrhoa))/x +    &
              9.71681713056504*LOG(zrhoa)**2 -    &
              0.1150478558347306*zt*LOG(zrhoa)**2 +    &
              0.0001570982486038294*zt**2*LOG(zrhoa)**2 +    &
              4.009144680125015e-7*zt**3*LOG(zrhoa)**2 +    &
              (0.7118597859976135*LOG(zrhoa)**2)/x -    &
              1.056105824379897*LOG(zrh)*LOG(zrhoa)**2 +    &
              0.00903377584628419*zt*LOG(zrh)*LOG(zrhoa)**2 -    &
              0.00001984167387090606*zt**2*LOG(zrh)*LOG(zrhoa)**2 +    &
              2.460478196482179e-8*zt**3*LOG(zrh)*LOG(zrhoa)**2 -    &
              (0.05790872906645181*LOG(zrh)*LOG(zrhoa)**2)/x -    &
              0.1487119673397459*LOG(zrhoa)**3 +    &
              0.002835082097822667*zt*LOG(zrhoa)**3 -    &
              9.24618825471694e-6*zt**2*LOG(zrhoa)**3 +    &
              5.004267665960894e-9*zt**3*LOG(zrhoa)**3 -    &
              (0.01270805101481648*LOG(zrhoa)**3)/x
        
        zjnuc=EXP(zjnuc)      !   add. Eq. (12) [1/(cm^3s)]      


        ! Equation (13) - total number of molecules in the critical cluster

        zntot =-0.002954125078716302_dp - 0.0976834264241286*zt +   &
               0.001024847927067835*zt**2 - 2.186459697726116e-6*zt**3 -    &
               0.1017165718716887/x - 0.002050640345231486*LOG(zrh) -   &
               0.007585041382707174*zt*LOG(zrh) +    &
               0.0001926539658089536*zt**2*LOG(zrh) -   &
               6.70429719683894e-7*zt**3*LOG(zrh) -    &
               (0.2557744774673163*LOG(zrh))/x +   &
               0.003223076552477191*LOG(zrh)**2 +   &
               0.000852636632240633*zt*LOG(zrh)**2 -    &
               0.00001547571354871789*zt**2*LOG(zrh)**2 +   &
               5.666608424980593e-8*zt**3*LOG(zrh)**2 +    &
               (0.03384437400744206*LOG(zrh)**2)/x +   &
               0.04743226764572505*LOG(zrh)**3 -    &
               0.0006251042204583412*zt*LOG(zrh)**3 +   &
               2.650663328519478e-6*zt**2*LOG(zrh)**3 -    &
               3.674710848763778e-9*zt**3*LOG(zrh)**3 -   &
               (0.0002672510825259393*LOG(zrh)**3)/x -    &
               0.01252108546759328*LOG(zrhoa) !+   &

        zntot =zntot + &
               0.005806550506277202*zt*LOG(zrhoa) -    &
               0.0001016735312443444*zt**2*LOG(zrhoa) +   &
               2.881946187214505e-7*zt**3*LOG(zrhoa) +    &
               (0.0942243379396279*LOG(zrhoa))/x -   &
               0.0385459592773097*LOG(zrh)*LOG(zrhoa) -   & 
               0.0006723156277391984*zt*LOG(zrh)*LOG(zrhoa) + &
               2.602884877659698e-6*zt**2*LOG(zrh)*LOG(zrhoa) +    &
               1.194163699688297e-8*zt**3*LOG(zrh)*LOG(zrhoa) -   &
               (0.00851515345806281*LOG(zrh)*LOG(zrhoa))/x -    &
               0.01837488495738111*LOG(zrh)**2*LOG(zrhoa) +   &
               0.0001720723574407498*zt*LOG(zrh)**2*LOG(zrhoa) -   & 
               3.717657974086814e-7*zt**2*LOG(zrh)**2*LOG(zrhoa) -    &
               5.148746022615196e-10*zt**3*LOG(zrh)**2*LOG(zrhoa) +    &
               (0.0002686602132926594*LOG(zrh)**2*LOG(zrhoa))/x -   &
               0.06199739728812199*LOG(zrhoa)**2 +    &
               0.000906958053583576*zt*LOG(zrhoa)**2 -   &
               9.11727926129757e-7*zt**2*LOG(zrhoa)**2 -    &
               5.367963396508457e-9*zt**3*LOG(zrhoa)**2 -   &
               (0.007742343393937707*LOG(zrhoa)**2)/x +    &
               0.0121827103101659*LOG(zrh)*LOG(zrhoa)**2 -   &
               0.0001066499571188091*zt*LOG(zrh)*LOG(zrhoa)**2 +    &
               2.534598655067518e-7*zt**2*LOG(zrh)*LOG(zrhoa)**2 -    &
               3.635186504599571e-10*zt**3*LOG(zrh)*LOG(zrhoa)**2 +    &
               (0.0006100650851863252*LOG(zrh)*LOG(zrhoa)**2)/x +   &
               0.0003201836700403512*LOG(zrhoa)**3 -    &
               0.0000174761713262546*zt*LOG(zrhoa)**3 +   &
               6.065037668052182e-8*zt**2*LOG(zrhoa)**3 -    &
               1.421771723004557e-11*zt**3*LOG(zrhoa)**3 +   &
               (0.0001357509859501723*LOG(zrhoa)**3)/x
               
        zntot = EXP(zntot)  !  add. Eq. (13)

          
        ! Equation (14) - radius of the critical cluster in nm

!        zrc=EXP(-1.6524245_dp+0.42316402*x+0.33466487*LOG(zntot)) ! [nm]

        ! Conversion [nm -> m]

!        zrxc(jl)=zrc*1e-9

        !----1.2) Limiter

         IF(zjnuc<1.e-7_dp .OR. zntot<4.0_dp) zjnuc=0.0

        ! limitation to 1E+10 [1/cm3s]
      
        zjnuc=MIN(zjnuc,1.e10_dp)

        pxtrnucr(jl,jk) = zjnuc

        ! mz_ak_20041102+
        ! convert total number of molecules in the critical cluster
        ! to number of sulfate molecules:

        pntot(jl,jk)=zntot*zxmole
        ! mz_ak_20041102-

        ELSE ! pmolecH2SO4, ptp1 , prhd out of range

        pntot(jl,jk)   =0.0_dp
        pxtrnucr(jl,jk)=0.0_dp

        END IF

     END DO ! kproma
  END DO ! klev

END SUBROUTINE nucl_vehkamaeki

END MODULE messy_m7
