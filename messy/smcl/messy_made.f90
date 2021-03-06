!        1         2         3         4         5         6         7         8
!2345678901234567890123456789012345678901234567890123456789012345678901234567890

MODULE messy_made

  ! MODULE FOR MADE CORE
  !
  ! MADE adopted to the structure of the Modular Earth Submodel System (MESSy).
  ! MADE was originally implemented by A. Lauer, DLR Oberpfaffenhofen, 2001-2003
  ! Original MADE source code (box model) by I. Ackermann, et al.,
  ! University Cologne, Germany, 1998.

  ! Version history:
  ! v1.1 by Axel Lauer, DLR
  !      Modification by Mattia Righi, DLR, 2011:
  !      - new dust size distribution parameters according to
  !        Dentener et al., 2006 and Pringle et al., 2010
  ! v1.2 by Christopher Kaiser, DLR, 2013 (christopher.kaiser@dlr.de)
  !      - removed cloud-related calculations (now done by CLOUD)
  !      - several mainly cosmetic changes after running MESSy's `gmake check'
  !        and `gmake messycheck'

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, SP, I4, I8, &
      pi, GRAV => g, AVO => N_A, RGASUNIV => R_gas,   &
      BOLTZ => k_B, MWAIR => M_air, MWH2O => M_H2O,   &
      MWC => MC, MWCl => MCl, MWN => MN, MWNa => MNa, &
      MWS => MS, RHOH2O => rho_h2o

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! ----------- <

  ! GLOBAL PARAMETER
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'made'   ! name of module
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.2'    ! module version

  ! CTRL-NAMELIST PARAMETERS (and defaults)
!DEBUG+
  LOGICAL          :: l_aeroproc = .true.
  LOGICAL          :: l_eqsam    = .true.
  LOGICAL          :: l_coag     = .true.
  LOGICAL          :: l_nuclcond = .true.
  LOGICAL          :: l_rename   = .true.
!DEBUG-
  LOGICAL, PUBLIC  :: lmade   = .false.     ! MADE on/off switch
  REAL(dp), PUBLIC :: SGININ  = 1.70_dp     ! geometric std dev., Aitken mode
  REAL(dp), PUBLIC :: SGINIA  = 2.00_dp     ! geometric std. dev., acc. mode
  REAL(dp), PUBLIC :: SGINIC  = 2.20_dp     ! geometric std. dev., coarse mode
  REAL(dp), PUBLIC :: BCTIME  = 59888.0_dp  ! half-life, BC (soot-aging) [s]
  REAL(dp), PUBLIC :: POMTIME = 59888.0_dp  ! half-life, POM (soot aging) [s]
  LOGICAL          :: l_nuc   = .true.      ! switch for nucleation

  ! GLOBAL PARAMETERS ========================================================

  !--- 1.) Numbers of compounds and modes of MADE:

  INTEGER, PUBLIC, PARAMETER :: nmod    = 3   ! number of aerosol modes
  INTEGER, PUBLIC, PARAMETER :: nspcsda = 27  ! number of all compounds
  INTEGER, PUBLIC, PARAMETER :: akn     = 1   ! index Aitken mode
  INTEGER, PUBLIC, PARAMETER :: acc     = 2   ! index acc. mode
  INTEGER, PUBLIC, PARAMETER :: cor     = 3   ! index coarse mode

  !--- 2.) List of indeces corresponding to the compound masses and
  !        mode numbers (in MADE internal array 'cblk'):

   INTEGER, PARAMETER, PUBLIC :: VSO4AJ    =  1  ! accumulation mode sulfate
   INTEGER, PARAMETER, PUBLIC :: VSO4AI    =  2  ! Aitken mode sulfate
!  INTEGER, PARAMETER, PUBLIC :: VSO4AC    =  3  ! coarse mode sulfate (fly ash)
   INTEGER, PARAMETER, PUBLIC :: VNH4AJ    =  3  ! accumulation mode ammonium
   INTEGER, PARAMETER, PUBLIC :: VNH4AI    =  4  ! Aitken mode ammonium
   INTEGER, PARAMETER, PUBLIC :: VNO3AJ    =  5  ! accumulation mode nitrate
   INTEGER, PARAMETER, PUBLIC :: VNO3AI    =  6  ! Aitken mode nitrate
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGARO1J =  7  ! accumulation mode
!sorgam                                             ! anthropogenic organic
!sorgam                                             ! aerosol from aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGARO1I =  8  ! Aitken mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGARO2J =  9  ! acc. mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGARO2I = 10  ! Aitken mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGALK1J = 11  ! acc. mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! alkanes and others except
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGALK1I = 12  ! Aitken mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! alkanes and others except
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGOLE1J = 13  ! acc. mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! alkanes and others except
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGOLE1I = 14  ! Aitken mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! alkanes and others except
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGBA1J  = 15  ! accumulation mode biogenic
!sorgam                                             ! aerosol from a-pinene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGBA1I  = 16  ! Aitken mode biogenic
!sorgam                                             ! aerosol from a-pinene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGBA2J  = 17  ! accumulation mode biogenic
!sorgam                                             ! aerosol from a-pinene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGBA2I  = 18  ! Aitken mode biogenic
!sorgam                                             ! aerosol from a-pinene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGBA3J  = 19  ! accumulation mode biogenic
!sorgam                                             ! aerosol from limonene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGBA3I  = 20  ! Aitken mode biogenic
!sorgam                                             ! aerosol from limonene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGBA4J  = 21  ! accumulation mode biogenic
!sorgam                                             ! aerosol from limonene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VORGBA4I  = 22  ! Aitken mode biogenic
!sorgam                                             ! aerosol from limonene
   INTEGER, PARAMETER, PUBLIC :: VORGPAJ   =  7  ! acc. mode (primary) organic
                                             ! aerosol
   INTEGER, PARAMETER, PUBLIC :: VORGPAI   =  8  ! Aitken mode (primary) organic
                                             ! aerosol
   INTEGER, PARAMETER, PUBLIC :: VECJ      =  9  ! acc. mode aerosol elemental
                                             ! carbon
   INTEGER, PARAMETER, PUBLIC :: VECI      = 10  ! Aitken mode elemental carbon
!unused   INTEGER, PARAMETER, PUBLIC :: VP25AJ    = 27  ! accumulation mode primary
!unused                                             ! PM2.5
!unused   INTEGER, PARAMETER, PUBLIC :: VP25AI    = 28  ! Aitken mode primary PM2.5
!unused   INTEGER, PARAMETER, PUBLIC :: VANTHA    = 29  ! coarse mode anthropogenic
!unused                                             ! aerosol
!unused   INTEGER, PARAMETER, PUBLIC :: VSEAS     = 30  ! coarse mode marine aerosol
!unused   INTEGER, PARAMETER, PUBLIC :: VSOILA    = 31  ! coarse mode soil-derived
!unused                                             ! aerosol
   INTEGER, PARAMETER, PUBLIC :: VNU0      = 11  ! Aitken mode number
   INTEGER, PARAMETER, PUBLIC :: VAC0      = 12  ! accumulation mode number
   INTEGER, PARAMETER, PUBLIC :: VCORN     = 13  ! coarse mode number
   INTEGER, PARAMETER, PUBLIC :: VH2OAJ    = 14  ! accumulation mode aerosol water
   INTEGER, PARAMETER, PUBLIC :: VH2OAI    = 15  ! Aitken mode aerosol water
   INTEGER, PARAMETER, PUBLIC :: VH2OAC    = 16  ! coarse mode aerosol water
   INTEGER, PARAMETER, PUBLIC :: VNU3      = 17  ! Aitken mode 3'rd moment
   INTEGER, PARAMETER, PUBLIC :: VAC3      = 18  ! accumulation mode 3'rd moment
   INTEGER, PARAMETER, PUBLIC :: VCOR3     = 19  ! coarse mode 3rd moment
!sorgam   INTEGER, PARAMETER, PUBLIC :: VCVARO1   = 40  ! cond. vapor from aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VCVARO2   = 41  ! cond. vapor from aromatics
!sorgam   INTEGER, PARAMETER, PUBLIC :: VCVALK1   = 42  ! cond. vapor from anth.
!sorgam                                             ! alkanes
!sorgam   INTEGER, PARAMETER, PUBLIC :: VCVOLE1   = 43  ! cond. vapor from anth.
!sorgam                                             ! olefines
!sorgam   INTEGER, PARAMETER, PUBLIC :: VCVAPI1   = 44  ! cond. vapor from a-pinene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VCVAPI2   = 45  ! cond. vapor from a-pinene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VCVLIM1   = 46  ! cond. vapor from limonene
!sorgam   INTEGER, PARAMETER, PUBLIC :: VCVLIM2   = 47  ! cond. vapor from limonene
   INTEGER, PARAMETER, PUBLIC :: VSULF     = 20  ! sulfuric acid vapor
   INTEGER, PARAMETER, PUBLIC :: VHNO3     = 21  ! nitric acid vapor
   INTEGER, PARAMETER, PUBLIC :: VNH3      = 22  ! ammonia gas
!   INTEGER, PARAMETER, PUBLIC :: VHCl      = 23  ! hydrochloric acid (gas phase)
   INTEGER, PARAMETER, PUBLIC :: VSEASI    = 23  ! Aitken mode sea salt
   INTEGER, PARAMETER, PUBLIC :: VSEASJ    = 24  ! accumulation mode sea salt
   INTEGER, PARAMETER, PUBLIC :: VSEASC    = 25  ! coarse mode sea salt
   INTEGER, PARAMETER, PUBLIC :: VDUSTJ    = 26  ! accumulation mode dust
   INTEGER, PARAMETER, PUBLIC :: VDUSTC    = 27  ! coarse mode dust

  !--- 3.) Definition of the modes of MADE:

  ! *** Standard deviation for the modes (set via namelist):

  REAL(dp), PUBLIC :: sigma(nmod) = (/ 1.7_dp, 2.0_dp, 2.2_dp /)

  ! *** Initial geometric mean diameter:

  REAL(dp), PARAMETER, PUBLIC :: dgini(nmod) = &
       (/ 0.01e-6_dp, 0.10e-6_dp, 1.00e-6_dp /)

  ! *** Minimum particle number concentrations [1/m3]:

  REAL(dp), PARAMETER, PUBLIC :: nummin(nmod) = (/ 1000.0_dp, 10.0_dp, 0.1_dp /)

  !--- 4.) Parameters for lognormal particle size distributions
  !        (calulated in made_initialize_core):

  REAL(dp), PUBLIC :: massmin(nmod)   ! minimum mass (total)  (ug/m3)
  REAL(dp), PUBLIC :: mminso4(nmod)   ! minimum mass SO4      (ug/m3)
  REAL(dp), PUBLIC :: mminnh4(nmod)   ! minimum mass NH4      (ug/m3)
  REAL(dp), PUBLIC :: mminno3(nmod)   ! minimum mass NO3      (ug/m3)
  REAL(dp), PUBLIC :: mminec(nmod)    ! minimum mass EC       (ug/m3)
  REAL(dp), PUBLIC :: mminorg(nmod)   ! minimum mass ORG      (ug/m3)
  REAL(dp), PUBLIC :: mminseas(nmod)  ! minimum mass sea salt (ug/m3)
  REAL(dp), PUBLIC :: mmindust(nmod)  ! minimum mass dust     (ug/m3)

  !--- 5.) Constants for calculating corresponding particle number
  !        contrations from mass emitted (dust).
  !        (Defined in subroutine made_initialize_core)

  REAL(dp), PUBLIC :: massfrac_du_j   ! mass fraction dust, acc. mode
  REAL(dp), PUBLIC :: massfrac_du_c   ! mass fraction dust, coarse mode
  REAL(dp), PUBLIC :: DU2AC0          ! acc. mode
  REAL(dp), PUBLIC :: DU2CORN         ! coarse mode

  ! *** Scalar variables for fixed standard deviations.

  REAL(dp) :: e1(nmod)        ! exp( log^2( sigmag )/8 )
  REAL(dp) :: es04(nmod)      !              " **4
  REAL(dp) :: es05(nmod)      !              " **5
  REAL(dp) :: es08(nmod)      !              " **8
  REAL(dp) :: es09(nmod)      !              " **9
  REAL(dp) :: es12(nmod)      !              " **12
  REAL(dp) :: es16(nmod)      !              " **16
  REAL(dp) :: es20(nmod)      !              " **20
  REAL(dp) :: es25(nmod)      !              " **25
  REAL(dp) :: es32(nmod)      !              " **32
  REAL(dp), PUBLIC :: es36(nmod)      !              " **36
  REAL(dp) :: es49(nmod)      !              " **49
  REAL(dp) :: es64(nmod)      !              " **64
  REAL(dp) :: es100(nmod)     !              " **100
  REAL(dp) :: xxlsg(nmod)     ! log(sigma)
  REAL(dp) :: l2sigma(nmod)   ! log(sigma) ** 2

  !--- 5.) Physical constants: -------------------------------------------------

  ! *** component densities [kg/m3]:

  REAL(dp), PARAMETER, PUBLIC :: RHOSO4  = 1.8e3_dp ! bulk dens. aero. sulfate
  REAL(dp), PARAMETER, PUBLIC :: RHONH4  = 1.8e3_dp ! bulk dens. aero. ammonium
  REAL(dp), PARAMETER, PUBLIC :: RHONO3  = 1.8e3_dp ! bulk dens. aero. nitrate
  REAL(dp), PARAMETER, PUBLIC :: RHOORG  = 1.0e3_dp ! bulk dens. aero. organics
  REAL(dp), PARAMETER, PUBLIC :: RHOSOIL = 2.5e3_dp ! bulk dens. aero. soil dust
                                                    ! (AeroCom 2000)
! REAL(dp), PARAMETER, PUBLIC :: RHOSOIL = 2.6e3_dp  ! bulk dens. aero. soil dust
                                                    ! (original MADE)
  REAL(dp), PARAMETER, PUBLIC :: RHOSEAS = 2.2e3_dp ! bulk dens. marine aerosol
  REAL(dp), PARAMETER, PUBLIC :: RHOANTH = 2.2e3_dp ! bulk dens. anth. aerosol
  PUBLIC :: RHOH2O

  ! *** molecular weights [g/mol]

  REAL(dp), PARAMETER, PUBLIC :: MWSO4   = 96.0576_dp  ! molecular weight for SO4
  REAL(dp), PARAMETER, PUBLIC :: MWHNO3  = 63.01287_dp ! molecular weight for HNO3
  REAL(dp), PARAMETER, PUBLIC :: MWNO3   = 62.0649_dp  ! molecular weight for NO3
  REAL(dp), PARAMETER, PUBLIC :: MWNH3   = 17.03061_dp ! molecular weight for NH3
  REAL(dp), PARAMETER, PUBLIC :: MWORG   = 180.0_dp    ! molec. weight for organic
                                                       ! species, first guess
  REAL(dp), PARAMETER, PUBLIC :: MWSEAS  = 58.443_dp   ! molecular weight for NaCl
  REAL(dp), PARAMETER, PUBLIC :: MWSOIL  = 40.08_dp    ! mol. w. for mineral dust
                                                       ! (adopted from M7)
  REAL(dp), PARAMETER, PUBLIC :: MWEC    = 12.011_dp   ! molec. w. for elemental
                                                       ! carbon
  REAL(dp), PARAMETER, PUBLIC :: MWH2SO4 = 98.07948_dp ! molec. w. for sulfuric
                                                       ! acid
  REAL(dp), PARAMETER, PUBLIC :: MWNH4   = 18.03858_dp ! molecular weight for NH4
!  REAL(dp), PARAMETER, PUBLIC :: MWHCl   = 36.4609_dp  ! molecular weight for
                                                       ! hydrochloric acid
  REAL(dp), PARAMETER :: MWCa    = 40.078_dp   ! molecular weight for elemental
                                               ! calcium
  REAL(dp), PARAMETER :: MWK     = 39.0983_dp  ! molecular weight for elemental
                                               ! potassium
  REAL(dp), PARAMETER :: MWMg    = 24.305_dp   ! molecular weight for elemental
                                               ! magnesium
  PUBLIC :: MWH2O, MWAIR

  ! *** mathematical constants:

  REAL(dp), PARAMETER :: TWOPI   = 2.0_dp * pi         ! 2*pi
  REAL(dp), PARAMETER :: THREEPI = 3.0_dp * pi         ! 3*pi
  REAL(dp), PARAMETER :: F6DPI   = 6.0_dp / pi         ! 6/pi
  REAL(dp), PARAMETER :: F6DPI9  = 1.0e9_dp  * F6DPI   ! 1.0e+9 * 6/pi
  REAL(dp), PARAMETER, PUBLIC :: F6DPIM9 = 1.0e-9_dp * F6DPI   ! 1.0e-9 * 6/pi
  REAL(dp), PARAMETER :: SQRT2   = 1.4142135623731_dp  ! sqrt(2)
  REAL(dp), PARAMETER :: LGSQT2  = 0.34657359027997_dp ! ln(sqrt(2))
  REAL(dp), PARAMETER :: DLGSQT2 = 1.0_dp / LGSQT2     ! 1/ln(sqrt(2))
  REAL(dp), PARAMETER :: ONE3    = 1.0_dp / 3.0_dp     ! 1/3
  REAL(dp), PARAMETER :: TWO3    = 2.0_dp / 3.0_dp     ! 2/3

  ! *** Factors for converting aerosol mass concentration [ug m**-3] to
  !     to 3rd moment concentration [m3 m^-3]

  REAL(dp), PARAMETER, PUBLIC :: SO4FAC  = F6DPIM9 / RHOSO4
  REAL(dp), PARAMETER, PUBLIC :: NH4FAC  = F6DPIM9 / RHONH4
  REAL(dp), PARAMETER, PUBLIC :: H2OFAC  = F6DPIM9 / RHOH2O
  REAL(dp), PARAMETER, PUBLIC :: NO3FAC  = F6DPIM9 / RHONO3
  REAL(dp), PARAMETER, PUBLIC :: ORGFAC  = F6DPIM9 / RHOORG
  REAL(dp), PARAMETER, PUBLIC :: SOILFAC = F6DPIM9 / RHOSOIL
  REAL(dp), PARAMETER, PUBLIC :: SEASFAC = F6DPIM9 / RHOSEAS
  REAL(dp), PARAMETER, PUBLIC :: ANTHFAC = F6DPIM9 / RHOANTH

  ! *** misc. constants:

  REAL(dp), PARAMETER :: P0  = 101325.0_dp ! starting std. surf. pressure [Pa]
  REAL(dp), PARAMETER :: T0  = 273.15_dp   ! standard temperature [K]
  REAL(dp), PARAMETER :: TS0 = 288.15_dp   ! starting std. surf. temperature [K]
  PUBLIC :: GRAV, AVO

  !--- 6.) Parameters for splitting of hydrophobic/hydrophilic soot (BC, OC).

  INTEGER, PARAMETER :: NUM_SPLIT_FACS = 8  ! number of splitting factors
                                                ! (for array dimensioning)
  INTEGER, PARAMETER, PUBLIC :: BC_I  = 1           ! BC, hydrophilic (Aitken mode)
  INTEGER, PARAMETER, PUBLIC :: BC_J  = 2           ! BC, hydrophilic (acc. mode)
  INTEGER, PARAMETER, PUBLIC :: BC2_I = 3           ! BC, hydrophobic (Aitken mode)
  INTEGER, PARAMETER, PUBLIC :: BC2_J = 4           ! BC, hydrophobic (acc. mode)
  INTEGER, PARAMETER, PUBLIC :: OC_I  = 5           ! OC, hydrophilic (Aitken mode)
  INTEGER, PARAMETER, PUBLIC :: OC_J  = 6           ! OC, hydrophilic (acc. mode)
  INTEGER, PARAMETER, PUBLIC :: OC2_I = 7           ! OC, hydrophobic (Aitken mode)
  INTEGER, PARAMETER, PUBLIC :: OC2_J = 8           ! OC, hydrophobic (acc. mode)

  ! END GLOBAL PARAMETERS ====================================================

  ! SUBROUTINES
  PUBLIC :: made_read_nml
  PUBLIC :: made_initialize_core
  PUBLIC :: made_main

   CONTAINS

! -------------------------------------------------------------------------

SUBROUTINE made_read_nml(status, iou)

    ! MADE MODULE ROUTINE (CORE)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    !--- Local variables:
!DEBUG+
!!$    NAMELIST /CTRL/ LMADE, SGININ, SGINIA, SGINIC, BCTIME, POMTIME
    NAMELIST /CTRL/ L_AEROPROC, L_EQSAM, L_COAG, L_NUCLCOND, L_RENAME, &
                    LMADE, SGININ, SGINIA, SGINIC, BCTIME, POMTIME, L_NUC
!DEBUG-

    CHARACTER(LEN=*), PARAMETER :: substr = 'made_read_nml'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE

    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! ---> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    !--- 1.) Read namelist:
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! consistency checks and diagnostic output
    WRITE(*,*) ''
    WRITE(*,*) ''
    WRITE(*,*) '------------------------------------------------------------'
    WRITE(*,*) '------------------------------------------------------------'
    WRITE(*,*) '---  Initialization of settings for aerosol module MADE  ---'
    WRITE(*,*) '---'
    WRITE(*,'(1X,A,L1)')        '---  lmade (enable/disable MADE)         = ', &
         lmade
    IF (lmade) THEN
       WRITE(*,'(1X,A,F4.2)')   '---  sginin (sigma Aitken mode)          = ', &
            sginin
       WRITE(*,'(1X,A,F4.2)')   '---  sginia (sigma acc. mode)            = ', &
            sginia
       WRITE(*,'(1X,A,F4.2)')   '---  sginic (sigma coarse mode)          = ', &
            sginic
       WRITE(*,'(1X,A,F7.1,A)') '---  bctime (half-life hydrophobic BC)   = ', &
            bctime, ' s'
       WRITE(*,'(1X,A,F7.1,A)') '---  pomtime (half-life hydrophobic POM) = ', &
            pomtime, ' s'
       WRITE(*,'(1X,A,L1)')     '---  l_nuc (enable/disable nucleation)   = ', &
            l_nuc
    ELSE
       WRITE(*,'(1X,A)') '---  ---> MADE switched off'
    END IF
    WRITE(*,*) '------------------------------------------------------------'
    WRITE(*,*) '------------------------------------------------------------'
    WRITE(*,*) ''
    WRITE(*,*) ''

    CALL read_nml_close(substr, iou, modstr)
 
    status = 0  ! no ERROR

END SUBROUTINE made_read_nml

! -------------------------------------------------------------------------

SUBROUTINE made_initialize_core

   ! Purpose:
   ! ---------
   ! Initializes constants and parameters used in the MADE aerosol model.
   !
   ! Author:
   ! ---------
   ! Axel Lauer, DLR
   !
   ! Interface:
   ! ---------
   ! *made_initialize_core* is called from *made_init* in *messy_made_e5*

   IMPLICIT NONE
   INTRINSIC LOG, EXP

   REAL(dp) :: rho_av   ! average density of a 'minimum'-particle
   REAL(dp) :: rho_sum  ! sum of the densities of all species forming a
                        ! 'minimum'-particle
   REAL(dp) :: mw_sum   ! sum of the molar weights of all species forming a
                        ! 'minimum'-particle
   REAL(dp) :: r        ! radius of log-normal distribution [m]
   REAL(dp) :: zsigma   ! standard deviation of log-normal distribution

   INTEGER  :: i        ! loop index

   ! ------- code starts here ------------

   SIGMA(akn) = SGININ
   SIGMA(acc) = SGINIA
   SIGMA(cor) = SGINIC

   ! *** Compute these once and they will all be saved in "COMMON"

   DO i = 1, nmod
      xxlsg(i)   = LOG(sigma(i))
      l2sigma(i) = xxlsg(i) ** 2
      e1(i)      = EXP(0.125_dp * l2sigma(i))
      es04(i)    = e1(i) ** 4
      es05(i)    = es04(i) * e1(i)
      es08(i)    = es04(i) * es04(i)
      es09(i)    = es04(i) * es05(i)
      es12(i)    = es04(i) * es04(i) * es04(i)
      es16(i)    = es08(i) * es08(i)
      es20(i)    = es16(i) * es04(i)
      es25(i)    = es16(i) * es09(i)
      es32(i)    = es16(i) * es16(i)
      es36(i)    = es16(i) * es20(i)
      es49(i)    = es25(i) * es20(i) * es04(i)
      es64(i)    = es32(i) * es32(i)
      es100(i)   = es36(i) * es64(i)
   END DO

   ! Calculate minimum mass conc. [ug/m3] allowed for each species,
   !
   ! assume the following molar fractions (Aitken + acc. mode):
   !
   ! 1 SO4 + 1 NO3 + 3 NH4 + 0 EC + 0 OC
   !
   ! ---> fully neutralized NH4(2)SO4 / NH4NO3 aerosol, no EC/OC

   rho_sum = rhoso4 + rhono3 + 3.0_dp * rhonh4
   rho_av  = rho_sum / 5.0_dp
   mw_sum  = mwso4 + mwno3 + 3.0_dp * mwnh4

   massmin(akn)  = 1.0e9_dp * rho_av * pi * nummin(akn) &
                   * dgini(akn)**3 * es36(akn) / 6.0_dp   ! [ug/m3]

   mminso4(akn)  =          mwso4 / mw_sum * massmin(akn) ! [ug/m3]
   mminnh4(akn)  = 3.0_dp * mwnh4 / mw_sum * massmin(akn) ! [ug/m3]
   mminno3(akn)  =          mwno3 / mw_sum * massmin(akn) ! [ug/m3]
   mminec(akn)   = 1.0e-30_dp                             ! [ug/m3]
   mminorg(akn)  = 1.0e-30_dp                             ! [ug/m3]
   mmindust(akn) = 0.0_dp                                 ! [ug/m3]
   mminseas(akn) = 1.0e-30_dp                             ! [ug/m3]

   massmin(acc)  = 1.0e9_dp * rho_av * pi * nummin(acc) &
                   * dgini(acc)**3 * es36(acc) / 6.0_dp   ! [ug/m3]

   mminso4(acc)  =          mwso4 / mw_sum * massmin(acc) ! [ug/m3]
   mminnh4(acc)  = 3.0_dp * mwnh4 / mw_sum * massmin(acc) ! [ug/m3]
   mminno3(acc)  =          mwno3 / mw_sum * massmin(acc) ! [ug/m3]
   mminec(acc)   = 1.0e-30_dp                             ! [ug/m3]
   mminorg(acc)  = 1.0e-30_dp                             ! [ug/m3]
   mmindust(acc) = 1.0e-30_dp                             ! [ug/m3]
   mminseas(acc) = 1.0e-30_dp                             ! [ug/m3]

   ! Coarse mode: 1 Seasalt + 1 Dust

   rho_sum = rhoseas + rhosoil
   rho_av  = rho_sum / 2.0_dp

   massmin(cor) = 1.0e9_dp * rho_av * pi * nummin(cor) &
                  * dgini(cor)**3 * es36(cor) / 6.0_dp    ! [ug/m3]

   mminso4(cor)  = 0.0_dp                                 ! [ug/m3]
   mminnh4(cor)  = 0.0_dp                                 ! [ug/m3]
   mminno3(cor)  = 0.0_dp                                 ! [ug/m3]
   mminec(cor)   = 0.0_dp                                 ! [ug/m3]
   mminorg(cor)  = 0.0_dp                                 ! [ug/m3]
   mminseas(cor) = 0.5_dp * massmin(cor)                  ! [ug/m3]
   mmindust(cor) = 0.5_dp * massmin(cor)                  ! [ug/m3]

   ! Factors for calculating particle number concentration from mass
   ! emitted for mineral dust.

   ! accumulation mode, AeroCom 2000 (global annual mean)

   r             = 0.21e-6_dp   ! radius [m]
   zsigma        = 1.59_dp          ! sigma
!   massfrac_du_j = 0.1095_dp        ! mass fraction ("dust_small_ncf")
   massfrac_du_j = 0.014_dp         ! mass fraction (Dentener et al., 2006)

   DU2AC0        = massfrac_du_j &
                   * 3.0_dp/(4.0_dp*pi*rhosoil*r**3*exp(4.5_dp*log(zsigma)**2))
   ! additional factor: #/kg ---> #/mol
   DU2AC0        = DU2AC0 * MWAIR * 1.0e-3_dp

   ! coarse mode, AeroCom 2000 (global annual mean)

   r             = 0.65e-6_dp   ! radius [m]
   zsigma        = 2.00_dp          ! sigma
!   massfrac_du_c = 0.8905_dp        ! mass fraction ("dust_small_ncf")
   massfrac_du_c = 0.986_dp         ! mass fraction (Dentener et al., 2006)

   DU2CORN       = massfrac_du_c &
                   * 3.0_dp/(4.0_dp*pi*rhosoil*r**3*exp(4.5_dp*log(zsigma)**2))
   ! additional factor: #/kg ---> #/mol
   DU2CORN       = DU2CORN * MWAIR * 1.0e-3_dp

END SUBROUTINE made_initialize_core

! -------------------------------------------------------------------------

SUBROUTINE made_main( BLKSIZE, NUMCELLS, PRESSURE, TEMPERATURE, RELHUM,  &
                      PTMST, PSO4RAT, PSOA, CLOUDCOVER, CBLK,            &
                      RH_HIST_AKN, RH_HIST_ACC, RH_HIST_COR,             &
                      DGNUC, DGACC, DGCOR, DGDRYNUC, DGDRYACC, DGDRYCOR, &
                      PDENSN, PDENSA, PDENSC )

   !   *** *MADE* Modal Aerosol Dynamics Model (for Europe).
   !
   !   Authors:
   !   ---------
   !   I. Ackermann, University of Cologne, Germany (original source)  1998
   !   A. Lauer, DLR (ECHAM-/Messy-version)                            2005
   !
   !   Purpose
   !   ---------
   !   Aerosol dynamics model (SO4,NH4,NO3,BC,POM,SS,DU,H2O), 3 modes.
   !
   !   Interface:
   !   ----------
   !   made_main is called from *messy_made_box* or *messy_made_e5*
   !
   !   Externals:
   !   ----------
   !
   !   *made_aeroproc* AEROSOL DYNAMICS DRIVER ROUTINE
   !
   !   Parameter list (i=input, o=output):
   !   -----------------------------------
   !
   !   i-  *BLKSIZE*     size of arrays
   !   i-  *NUMCELLS*    number of gridcells in arrays
   !   i-  *PRESSURE*    full level pressures [Pa]
   !   i-  *TEMPERATURE* temperature [K]
   !   i-  *RELHUM*      relative humidity [frac] (0-1)
   !   i-  *PTMST*       time step [s]
   !   i-  *PSO4RAT*     H2SO4(g), rate of change [ug/m3/s]
   !   i-  *PSOA*        SOA (gas phase) "emissions" [ug/m3/s]
   !   i-  *CLOUDCOVER*  fractional cloud cover (0-1)
   !   io  *CBLK*        tracer concentration array (gas+aerosol)
   !   io  *RH_HIST_AKN* rel. hum. history (---> hysteresis)
   !   io  *RH_HIST_ACC* rel. hum. history (---> hysteresis)
   !   io  *RH_HIST_COR* rel. hum. history (---> hysteresis)
   !   -o  *DGNUC*       modal (wet) diameter, Aitken mode [m]
   !   -o  *DGACC*       modal (wet) diameter, accumulation mode [m]
   !   -o  *DGCOR*       modal (wet) diameter, coarse mode [m]
   !   -o  *DGDRYNUC*    modal (dry) diameter, Aitken mode [m]
   !   -o  *DGDRYACC*    modal (dry) diameter, accumulation mode [m]
   !   -o  *DGDRYCOR*    modal (dry) diameter, coarse mode [m]
   !   -o  *PDENSN*      average density (wet), Aitken mode [kg/m3]
   !   -o  *PDENSA*      average density (wet), accumulation mode [kg/m3]
   !   -o  *PDENSC*      average density (wet), coarse mode [kg/m3]

   IMPLICIT NONE
   INTRINSIC MAX

   ! input/output parameters

   INTEGER,  INTENT(in)    :: BLKSIZE               ! size of arrays
   INTEGER,  INTENT(in)    :: NUMCELLS              ! actual number of grid
                                                    ! cells
   REAL(dp), INTENT(in)    :: PTMST                 ! time step [s]
   REAL(dp), INTENT(in)    :: CLOUDCOVER(BLKSIZE)   ! frac. cloud cover (0-1)
   REAL(dp), INTENT(in)    :: PSO4RAT(BLKSIZE)      ! H2SO4(g), rate of change
                                                    ! [ug/m3/s]
   REAL(dp), INTENT(in)    :: PRESSURE(BLKSIZE)     ! pressure[Pa]
   REAL(dp), INTENT(in)    :: PSOA(BLKSIZE)         ! SOA "emissions"
                                                    ! [ug/m3/s]
   REAL(dp), INTENT(in)    :: TEMPERATURE(BLKSIZE)  ! temperature [K]

   REAL(dp), INTENT(inout) :: RH_HIST_ACC(BLKSIZE)  ! rel. hum. history
                                                    ! (---> hysteresis)
   REAL(dp), INTENT(inout) :: RH_HIST_AKN(BLKSIZE)  ! rel. hum. history
                                                    ! (---> hysteresis)
   REAL(dp), INTENT(inout) :: RH_HIST_COR(BLKSIZE)  ! rel. hum. history
                                                    ! (---> hysteresis)
   REAL(dp), INTENT(inout) :: CBLK(BLKSIZE,NSPCSDA) ! tracer conc.
   REAL(dp), INTENT(inout) :: RELHUM(BLKSIZE)       ! rel. humidity [frac]
                                                    ! (0-1)

   ! modal (wet) diameters [m]

   REAL(dp), INTENT(out) :: DGNUC(BLKSIZE)          ! Aitken mode
   REAL(dp), INTENT(out) :: DGACC(BLKSIZE)          ! accumulation mode
   REAL(dp), INTENT(out) :: DGCOR(BLKSIZE)          ! coarse mode

   ! modal (dry) diameters [m]

   REAL(dp), INTENT(out) :: DGDRYNUC(BLKSIZE)       ! Aitken mode
   REAL(dp), INTENT(out) :: DGDRYACC(BLKSIZE)       ! accumulation mode
   REAL(dp), INTENT(out) :: DGDRYCOR(BLKSIZE)       ! coarse mode

   ! average modal density (wet) [kg/m3]

   REAL(dp), INTENT(out) :: PDENSN(BLKSIZE)         ! Aitken mode
   REAL(dp), INTENT(out) :: PDENSA(BLKSIZE)         ! accumulation mode
   REAL(dp), INTENT(out) :: PDENSC(BLKSIZE)         ! coarse mode

   ! local variables

   INTEGER :: NCELL                    ! loop index (gridcells)
   INTEGER :: LDUM                     ! aux. loop index

   REAL(dp) :: SO4RAT_IN(BLKSIZE)      ! H2SO4(g) production rate [ug/m3/s]
   REAL(dp) :: CLOUDFREE               ! cloud free fraction of grid box
   REAL(dp) :: CBLK_M1(BLKSIZE,NSPCSDA)! copy of 'CBLK' array before MADE;
                                       ! used to update tracer changes

!sorgam   REAL(dp) :: drog_in(BLKSIZE,LDROG)     ! accumulated production of
!sorgam                                          ! anthropogenic AND biogenic
!sorgam                                          ! organic aerosol precursor
!sorgam                                          ! [ug m-3]

   REAL(dp) :: XLM(BLKSIZE)       ! atmospheric mean free path [m]
   REAL(dp) :: AMU(BLKSIZE)       ! atmospheric dynamic viscosity [kg/m/s]

   ! modal mass concentrations [ug/m3]

   REAL(dp) :: PMASSN(BLKSIZE)    ! Aitken mode
   REAL(dp) :: PMASSA(BLKSIZE)    ! accumulation mode
   REAL(dp) :: PMASSC(BLKSIZE)    ! coarse mode

   REAL(dp) :: KNNUC(BLKSIZE)     ! nuclei mode  Knudsen number
   REAL(dp) :: KNACC(BLKSIZE)     ! accumulation Knudsen number
   REAL(dp) :: KNCOR(BLKSIZE)     ! coarse mode Knudsen number

   REAL(dp) :: FCONCN(BLKSIZE)      ! condensation rate (SO4), Aitken mode
   REAL(dp) :: FCONCA(BLKSIZE)      ! condensation rate (SO4), acc. mode
   REAL(dp) :: FCONCN_ORG (BLKSIZE) ! condensation rate (SOA), Aitken mode
   REAL(dp) :: FCONCA_ORG (BLKSIZE) ! condensation rate (SOA), acc. mode

   REAL(dp) :: DMDT(BLKSIZE)      ! rate of production of new so4 mass
                                  ! concentration by particle formation

   REAL(dp) :: DNDT(BLKSIZE)      ! rate of production of new particle number
                                  ! concentration by particle formation

   REAL(dp) :: CGRN3(BLKSIZE)     ! Aitken mode
   REAL(dp) :: CGRA3(BLKSIZE)     ! accumulation mode

   REAL(dp) :: URN00(BLKSIZE)     ! Aitken mode 0th moment
                                  ! self-coagulation rate
   REAL(dp) :: URA00(BLKSIZE)     ! accumulation mode 0th moment
                                  ! self-coagulation rate

   REAL(dp) :: BRNA01(BLKSIZE)    ! rate for 0th moment

   REAL(dp) :: BRNA31(BLKSIZE)    ! intermodal 3rd moment transfer rate by
                                  ! by intermodal coagulation

   REAL(dp) :: DELTASO4A(BLKSIZE) ! increment of concentration added to
                                  ! sulfate aerosol by condensation [ug/m3]

   REAL(dp) :: MERGENUM(BLKSIZE)  ! mode merging, number [#/m3]
                                  ! (for diagnostics only)
   REAL(dp) :: MERGEM3(BLKSIZE)   ! mode merging, 3rd moment [m3/m3]
                                  ! (for diagnostics only)

   REAL(dp) :: nu3dry             ! 3rd moments, dry (aux. variables)
   REAL(dp) :: ac3dry
   REAL(dp) :: cor3dry

   ! *** parameter ***

   REAL(dp), PARAMETER :: CONMIN  = 1.0E-30_dp  ! conc. lower limit [ug/m3]
   REAL(dp), PARAMETER :: DGMIN   = 1.0E-9_dp   ! lowest particle diameter [m]

   ! -------------------------------- Start code ------------------------------

   ! Validate input parameters and make a copy of current tracer
   ! concentrations (gas+aerosol) for calculating the updated
   ! tracer concentrations after aerosol dynamics calculations.

   DO NCELL = 1, NUMCELLS      ! loop over gridcells

      DO LDUM = 1, NSPCSDA
         CBLK_M1(NCELL,LDUM) = CBLK(NCELL,LDUM)
      END DO

      ! Reset diagnostic arrays
      MERGENUM(NCELL) = 0.0_dp
      MERGEM3(NCELL) = 0.0_dp

      ! copy H2SO4(g) production rate [ug/m3/s]
      SO4RAT_IN(NCELL) = PSO4RAT(NCELL)

   END DO    ! loop over gridcells (NCELL)

   ! *** CALL AEROSOL DYNAMICS DRIVER ***

   IF (l_aeroproc) THEN
   CALL MADE_AEROPROC( BLKSIZE, NUMCELLS, CBLK, PTMST,           &
                       TEMPERATURE, PRESSURE, RELHUM,            &
                       SO4RAT_IN, PSOA,                          &
!sorgam                ORGARO1RAT, ORGARO2RAT,                   &
!sorgam                ORGALK1RAT, ORGOLE1RAT,                   &
!sorgam                ORGBIO1RAT, ORGBIO2RAT,                   &
!sorgam                ORGBIO3RAT, ORGBIO4RAT,                   &
!sorgam                DROG, LDROG, NCV, NACV,                   &
                       XLM, AMU, DGNUC, DGACC, DGCOR,            &
                       PMASSN, PMASSA, PMASSC,                   &
                       PDENSN, PDENSA, PDENSC,                   &
                       KNNUC, KNACC, KNCOR, FCONCN, FCONCA,      &
                       FCONCN_ORG, FCONCA_ORG,                   &
                       DMDT, DNDT, CGRN3, CGRA3,                 &
                       URN00, URA00, BRNA01, BRNA31,             &
                       DELTASO4A,                                &
                       RH_HIST_AKN, RH_HIST_ACC, RH_HIST_COR,    &
                       MERGENUM, MERGEM3 )

   ! Aerosol dynamics are only calculated for cloud free conditions.
   ! Within a cloud covered area, all aerosol dynamics are suspended.
   ! Thus, the tendencies calculated by SUBROUTINE AERO_DRIVER are
   ! scaled with (1 - cloud cover) to get a new mean value for the
   ! whole grid box.

   ! cloudy area of grid box = cloud cover
   ! cloud free area         = 1 - cloud cover

   DO NCELL = 1, NUMCELLS

      ! CBLK(new) - CBLK_M1 = (CBLK(current)-CBLK_M1) * (1-cloud cover)
      CLOUDFREE = 1.0_dp - CLOUDCOVER(NCELL)

      DO LDUM = 1, NSPCSDA
         CBLK(NCELL,LDUM) = (CBLK(NCELL,LDUM) - CBLK_M1(NCELL,LDUM)) &
                            * CLOUDFREE + CBLK_M1(NCELL,LDUM)
      END DO
   END DO

   ! *** Get new distribution information:

   CALL MADE_MODPAR( BLKSIZE, NUMCELLS,           &
                     CBLK, TEMPERATURE, PRESSURE, &
                     PMASSN, PMASSA, PMASSC,      &
                     PDENSN, PDENSA, PDENSC,      &
                     XLM, AMU,                    &
                     DGNUC, DGACC, DGCOR,         &
                     KNNUC, KNACC, KNCOR )

   ! Calculate modal diameters (dry) from 3rd moment.
   ELSE
      XLM        = 6.e-8_dp
      AMU        = 2.e-5_dp
      DGNUC      = 1.e-8_dp
      DGACC      = 1.e-7_dp
      DGCOR      = 1.e-6_dp
      PMASSN     = 2._dp
      PMASSA     = 20._dp
      PMASSC     = 40._dp
      PDENSN     = (rhoso4 + rhono3 + 3._dp * rhonh4) / 5._dp
      PDENSA     = PDENSN
      PDENSC     = (rhoseas + rhosoil) / 2._dp
      KNNUC      = 2._dp * XLM / DGNUC
      KNACC      = 2._dp * XLM / DGACC
      KNCOR      = 2._dp * XLM / DGCOR
      FCONCN     = 0._dp
      FCONCA     = 0._dp
      FCONCN_ORG = 0._dp
      FCONCA_ORG = 0._dp
      DMDT       = 0._dp
      DNDT       = 0._dp
      CGRN3      = 0._dp
      CGRA3      = 0._dp
      URN00      = 0._dp
      URA00      = 0._dp
      BRNA01     = 0._dp
      BRNA31     = 0._dp
      DELTASO4A  = 0._dp
      MERGENUM   = 0._dp
      MERGEM3    = 0._dp
   END IF

   DO ncell = 1, numcells
      nu3dry  = MAX(conmin, cblk(ncell,vnu3)  - h2ofac * cblk(ncell,vh2oai))
      ac3dry  = MAX(conmin, cblk(ncell,vac3)  - h2ofac * cblk(ncell,vh2oaj))
      cor3dry = MAX(conmin, cblk(ncell,vcor3) - h2ofac * cblk(ncell,vh2oac))

      dgdrynuc(ncell) = MAX(dgmin, (nu3dry  / (cblk(ncell,vnu0)  &
                                    * es36(akn)))**one3)
      dgdryacc(ncell) = MAX(dgmin, (ac3dry  / (cblk(ncell,vac0)  &
                                    * es36(acc)))**one3)
      dgdrycor(ncell) = MAX(dgmin, (cor3dry / (cblk(ncell,vcorn) &
                                    * es36(cor)))**one3)
   END DO

END SUBROUTINE made_main

! -------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! PRIVATE ROUTINES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

SUBROUTINE MADE_AEROPROC ( BLKSIZE, NUMCELLS,                     &
                           CBLK, DT, BLKTA, BLKPRS,               &
                           BLKRH, SO4RAT, SOA_MADE,               &
!sorgam                    ORGARO1RAT, ORGARO2RAT,                &
!sorgam                    ORGALK1RAT, ORGOLE1RAT,                &
!sorgam                    ORGBIO1RAT, ORGBIO2RAT,                &
!sorgam                    ORGBIO3RAT, ORGBIO4RAT,                &
!sorgam                    DROG, LDROG, NCV, NACV,                &
                           XLM, AMU, DGNUC, DGACC, DGCOR,         &
                           PMASSN, PMASSA, PMASSC,                &
                           PDENSN, PDENSA, PDENSC,                &
                           KNNUC, KNACC, KNCOR,                   &
                           FCONCN, FCONCA,                        &
                           FCONCN_ORG, FCONCA_ORG,                &
                           DMDT, DNDT, CGRN3, CGRA3,              &
                           URN00, URA00, BRNA01, C30,             &
                           DELTASO4A,                             &
                           RH_HIST_AKN, RH_HIST_ACC, RH_HIST_COR, &
                           MERGENUM, MERGEM3 )

   !ia*********************************************************************
   !ia                                                                    *
   !ia     AEROSOL DYNAMICS DRIVER ROUTINE                                *
   !ia     based on MODELS3 formulation by FZB                            *
   !ia     Modified by IA in November 97                                  *
   !ia                                                                    *
   !ia     Revision history                                               *
   !ia     When    WHO     WHAT                                           *
   !ia     ----    ----    ----                                           *
   !ia     ????    FZB     BEGIN                                          *
   !ia     05/97   IA      Adapted for use in CTM2-S                      *
   !ia     11/97   IA      Modified for new model version                 *
   !ia                     see comments under iarev02                     *
   !ia     01/98   IA      Adapted for 3-D use                            *
   !ia                     see comments under ia3d                        *
   !ia                                                                    *
   !ia     Called BY:      RPMMOD3                                        *
   !ia                                                                    *
   !ia     Calls to:       EQL3, MODPAR, COAGRATE, NUCLCOND, AEROSTEP     *
   !ia                     GETVSED                                        *
   !ia                                                                    *
   !ia*********************************************************************

   IMPLICIT NONE

   ! *** INPUT ***

   INTEGER,  INTENT(in)    :: BLKSIZE               ! dimension of arrays
   INTEGER,  INTENT(in)    :: NUMCELLS              ! actual number of cells
                                                    ! in arrays
   REAL(dp), INTENT(inout) :: CBLK(BLKSIZE,NSPCSDA) ! main array of variables
   REAL(dp), INTENT(in)    :: DT                    ! time step [s]

   ! Meteorological information:

   REAL(dp), INTENT(in)    :: BLKTA(BLKSIZE)        ! air temperature [K]
   REAL(dp), INTENT(in)    :: BLKPRS(BLKSIZE)       ! air pressure in [Pa]
   REAL(dp), INTENT(in)    :: BLKRH(BLKSIZE)        ! frac. rel. humidity (0-1)
   REAL(dp), INTENT(inout) :: RH_HIST_AKN(BLKSIZE)  ! rel. hum. history
                                                    ! (---> hysteresis)
   REAL(dp), INTENT(inout) :: RH_HIST_ACC(BLKSIZE)  ! rel. hum. history
                                                    ! (---> hysteresis)
   REAL(dp), INTENT(inout) :: RH_HIST_COR(BLKSIZE)  ! rel. hum. history
                                                    ! (---> hysteresis)

   ! Chemical production rates: [ug/m3/s]

   REAL(dp), INTENT(inout) :: SO4RAT(BLKSIZE)       ! sulfate gas-phase
                                                    ! production rate
   REAL(dp), INTENT(in)    :: SOA_MADE(BLKSIZE)     ! SOA (gas phase)
                                                    ! "emissions" [ug/m3/s]

!sorgam   INTEGER, INTENT(in)    :: LDROG           ! # of organic aerosol
!sorgam                                             ! precursor`
!sorgam   INTEGER, INTENT(in)    :: NCV             ! total # of cond. vapors
!sorgam                                             ! & SOA species
!sorgam   INTEGER, INTENT(in)    :: NACV            ! # of anthrop. cond.
!sorgam                                             ! vapors & SOA species
!sorgam
!sorgam   !bs * organic condensable vapor production rate
!sorgam   REAL(dp) :: DROG(BLKSIZE,LDROG)    ! Delta ROG conc. [ppm]
!sorgam
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     aromatics
!sorgam   REAL(dp) :: ORGARO1RAT(BLKSIZE)
!sorgam 
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     aromatics
!sorgam   REAL(dp) :: ORGARO2RAT(BLKSIZE)
!sorgam 
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     alkanes & others
!sorgam   REAL(dp) :: ORGALK1RAT(BLKSIZE)
!sorgam 
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     alkenes & others
!sorgam   REAL(dp) :: ORGOLE1RAT(BLKSIZE)
!sorgam 
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO1RAT(BLKSIZE)
!sorgam 
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO2RAT(BLKSIZE)
!sorgam 
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO3RAT(BLKSIZE)
!sorgam 
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO4RAT(BLKSIZE)

   ! *** OUTPUT ***

   ! atmospheric properties:

   REAL(dp), INTENT(out)   :: XLM(BLKSIZE)        ! atmospheric mean free
                                                  ! path [m]
   REAL(dp), INTENT(out)   :: AMU(BLKSIZE)        ! atmospheric dynamic
                                                  ! viscosity [kg/m/s]

   ! modal diameters [m]:

   REAL(dp), INTENT(out)   :: DGNUC(BLKSIZE)      ! Aitken mode
   REAL(dp), INTENT(out)   :: DGACC(BLKSIZE)      ! accumulation mode
   REAL(dp), INTENT(out)   :: DGCOR(BLKSIZE)      ! coarse mode

   ! aerosol properties:

   ! modal mass concentrations [ug/m3]

   REAL(dp), INTENT(out)   :: PMASSN(BLKSIZE)     ! Aitken mode
   REAL(dp), INTENT(out)   :: PMASSA(BLKSIZE)     ! accumulation mode
   REAL(dp), INTENT(out)   :: PMASSC(BLKSIZE)     ! coarse mode

   ! average modal particle densities [kg/m3]

   REAL(dp), INTENT(out)   :: PDENSN(BLKSIZE)     ! Aitken mode
   REAL(dp), INTENT(out)   :: PDENSA(BLKSIZE)     ! accumulation mode
   REAL(dp), INTENT(out)   :: PDENSC(BLKSIZE)     ! coarse mode

   ! average modal Knudsen numbers

   REAL(dp), INTENT(out)   :: KNNUC(BLKSIZE)      ! Aitken mode
   REAL(dp), INTENT(out)   :: KNACC(BLKSIZE)      ! accumulation mode
   REAL(dp), INTENT(out)   :: KNCOR(BLKSIZE)      ! coarse mode

   ! modal condensation factors (see comments in NUCLCOND)

   REAL(dp), INTENT(out)   :: FCONCN(BLKSIZE)     ! SO4, Aitken mode
   REAL(dp), INTENT(out)   :: FCONCA(BLKSIZE)     ! SO4, accumulation mode
   REAL(dp), INTENT(out)   :: FCONCN_ORG(BLKSIZE) ! organics, Aitken mode
   REAL(dp), INTENT(out)   :: FCONCA_ORG(BLKSIZE) ! organics, accumulation mode

   ! rates for secondary particle formation

   ! production of new mass concentration [ug/m3/s]

   REAL(dp), INTENT(out)   :: DMDT(BLKSIZE)       ! rate of production of new
                                                  ! so4 mass by particle
                                                  ! formation

   ! production of new number concentration [number/m3/s]

   REAL(dp), INTENT(out)   :: DNDT(BLKSIZE)       ! rate of producton of new
                                                  ! particle number by
                                                  ! particle formation

   ! growth rate for third moment by condensation of precursor
   ! vapor on existing particles [3rd mom/m3/s]

   REAL(dp), INTENT(out)   :: CGRN3(BLKSIZE)      ! Aitken mode
   REAL(dp), INTENT(out)   :: CGRA3(BLKSIZE)      ! accumulation mode

   ! rates for coaglulation [m3/s]:

   ! unimodal rates

   REAL(dp), INTENT(out)   :: URN00(BLKSIZE)      ! Aitken mode 0th moment
                                                  ! self-coagulation rate
   REAL(dp), INTENT(out)   :: URA00(BLKSIZE)      ! acc. mode 0th moment
                                                  ! self-coagulation rate

   ! bimodal rates: Aitken mode with accumulation mode

   REAL(dp), INTENT(out)   :: BRNA01(BLKSIZE)     ! rate for 0th moment

   ! 3rd moment intermodal transfer rate replaces coagulation rate
   ! (FSB 1/12/98)

   REAL(dp), INTENT(out)   :: C30(BLKSIZE)        ! intermodal 3rd moment
                                                  ! transfer rate by
                                                  ! intermodal coagulation

   REAL(dp), INTENT(out)   :: DELTASO4A(BLKSIZE)  ! increment of conc. added to
                                                  ! sulfate aerosol by
                                                  ! condensation [ug/m3]
   REAL(dp), INTENT(out)   :: MERGENUM(BLKSIZE)   ! mode merging, number [#/m3]
                                                  ! (for diagnostics only)
   REAL(dp), INTENT(out)   :: MERGEM3(BLKSIZE)    ! mode merging, 3rd moment
                                                  ! [3rd mom./m3]
                                                  ! (for diagnostics only)

   ! ///////////////////// Begin code ///////////////////////////////////

   ! *** Get water, ammonium and nitrate content:

   IF (l_eqsam) THEN
   CALL MADE_EQL3X(    BLKSIZE, NUMCELLS,           &
                       CBLK, BLKTA, BLKRH,          &
                       RH_HIST_AKN, RH_HIST_ACC,    &
                       RH_HIST_COR )

   ! *** Get size distribution information:

   CALL MADE_MODPAR(   BLKSIZE, NUMCELLS,           &
                       CBLK, BLKTA, BLKPRS,         &
                       PMASSN, PMASSA, PMASSC,      &
                       PDENSN, PDENSA, PDENSC,      &
                       XLM, AMU,                    &
                       DGNUC, DGACC, DGCOR,         &
                       KNNUC, KNACC, KNCOR )
   ELSE
      XLM        = 6.e-8_dp
      AMU        = 2.e-5_dp
      DGNUC      = 1.e-8_dp
      DGACC      = 1.e-7_dp
      DGCOR      = 1.e-6_dp
      PMASSN     = 2._dp
      PMASSA     = 20._dp
      PMASSC     = 40._dp
      PDENSN     = (rhoso4 + rhono3 + 3._dp * rhonh4) / 5._dp
      PDENSA     = PDENSN
      PDENSC     = (rhoseas + rhosoil) / 2._dp
      KNNUC      = 2._dp * XLM / DGNUC
      KNACC      = 2._dp * XLM / DGACC
      KNCOR      = 2._dp * XLM / DGCOR
   END IF

   ! *** Calculate coagulation rates for fine particles:

   IF (l_coag) THEN
   CALL MADE_COAGRATE( BLKSIZE, NUMCELLS,           &
                       CBLK, BLKTA,                 &
                       PDENSN, PDENSA,              &
                       AMU,                         &
                       DGNUC, DGACC,                &
                       KNNUC, KNACC,                &
                       URN00, URA00,                &
                       BRNA01, C30 )
   ELSE
      URN00  = 0._dp
      URA00  = 0._dp
      BRNA01 = 0._dp
      C30    = 0._dp
   END IF

   ! *** Get condensation and particle formation (nucleation) rates:

   IF (l_nuclcond) THEN
   CALL MADE_NUCLCOND( BLKSIZE, NUMCELLS,          &
                       CBLK, DT,                   &
                       BLKTA, BLKPRS, BLKRH,       &
                       SO4RAT, SOA_MADE,           &
!sorgam                ORGARO1RAT, ORGARO2RAT,     &
!sorgam                ORGALK1RAT, ORGOLE1RAT,     &
!sorgam                ORGBIO1RAT, ORGBIO2RAT,     &
!sorgam                ORGBIO3RAT, ORGBIO4RAT,     &
!sorgam                DROG, LDROG, NCV, NACV,     &
                       DGNUC, DGACC,               &
                       FCONCN, FCONCA,             &
                       FCONCN_ORG, FCONCA_ORG,     &
                       DMDT, DNDT, DELTASO4A,      &
                       CGRN3, CGRA3 )
   ELSE
      FCONCN     = 0._dp
      FCONCA     = 0._dp
      FCONCN_ORG = 0._dp
      FCONCA_ORG = 0._dp
      DMDT       = 0._dp
      DNDT       = 0._dp
      DELTASO4A  = 0._dp
      CGRN3      = 0._dp
      CGRA3      = 0._dp
   END IF

   ! *** Advance forward in time  DT seconds:

   CALL MADE_AEROSTEP( BLKSIZE, NUMCELLS,          &
                       CBLK, DT,                   &
!sorgam                ORGARO1RAT, ORGARO2RAT,     &
!sorgam                ORGALK1RAT, ORGOLE1RAT,     &
!sorgam                ORGBIO1RAT, ORGBIO2RAT,     &
!sorgam                ORGBIO3RAT, ORGBIO4RAT,     &
                       DGNUC, DGACC,               &
                       FCONCN, FCONCA,             &
                       FCONCN_ORG, FCONCA_ORG,     &
                       DMDT, DNDT, DELTASO4A,      &
                       SOA_MADE,                   &
                       URN00, URA00,               &
                       BRNA01, C30,                &
                       CGRN3, CGRA3,               &
                       MERGENUM, MERGEM3 )

END SUBROUTINE MADE_AEROPROC

!------------------------------------------------------------------------------!

SUBROUTINE MADE_EQL3X( BLKSIZE, NUMCELLS, CBLK, BLKTA, BLKRH, &
                       RH_HIST_AKN, RH_HIST_ACC, RH_HIST_COR )

   !***********************************************************************
   !**  DESCRIPTION:
   !    Calculates the distribution of ammonia/ammonium, nitric
   !    acid/nitrate, and water between the gas and aerosol phases as the
   !    total sulfate, ammonia, and nitrate concentrations, relative
   !    humidity and temperature change. The evolution of the aerosol mass
   !    concentration due to the change in aerosol chemical composition is
   !    calculated.
   !**  REVISION HISTORY:
   !    NEW VERSION using EQSAM (Metzger et al.) 11/2002 by Axel Lauer
   !    prototype 1/95 by Uma and Carlie
   !    Revised   8/95 by US to calculate air density in stmt func
   !              and collect met variable stmt funcs in one include file
   !    Revised 7/26/96 by FSB to use block concept.
   !    Revise 12/1896 to do do i-mode calculation.
   !***********************************************************************

   IMPLICIT NONE

   ! *** INPUT ***

   INTEGER, INTENT(in)     :: BLKSIZE               ! dimension of arrays
   INTEGER, INTENT(in)     :: NUMCELLS              ! actual number of cells
                                                    ! in arrays

   REAL(dp), INTENT(inout) :: CBLK(BLKSIZE,NSPCSDA) ! main array of variables

   ! Meteorological information in blocked arays:

   REAL(dp), INTENT(in)    :: BLKTA(BLKSIZE)        ! air temperature [K]
   REAL(dp), INTENT(in)    :: BLKRH(BLKSIZE)        ! fractional rel. humidity
   REAL(dp), INTENT(inout) :: RH_HIST_AKN(BLKSIZE)  ! rel. humidity history,
                                                    ! Aitken mode
                                                    ! (---> hysteresis)
   REAL(dp), INTENT(inout) :: RH_HIST_ACC(BLKSIZE)  ! rel. humidity history,
                                                    ! acc. mode
                                                    ! (---> hysteresis)
   REAL(dp), INTENT(inout) :: RH_HIST_COR(BLKSIZE)  ! rel. humidity history,
                                                    ! coarse mode
                                                    ! (---> hysteresis)

   ! *** LOCAL ***

   INTEGER :: LCELL                   ! loop counter

   ! input arrays for EQSAM

   REAL(dp) :: inSO4(BLKSIZE)         ! SO4 (w/o seasalt associated SO4) [ug/m3]
   REAL(dp) :: inNO3(BLKSIZE)         ! aerosol NO3 [ug/m3]
   REAL(dp) :: inNH3(BLKSIZE)         ! NH3 (gas phase) [ug/m3]
   REAL(dp) :: inNH4(BLKSIZE)         ! NH4 [ug/m3]
!   REAL(dp) :: inH2SO4(BLKSIZE)       ! H2SO4 (gas phase) [ug/m3]
   REAL(dp) :: inHNO3(BLKSIZE)        ! HNO3 (gas phase) [ug/m3]
!   REAL(dp) :: inHCl(BLKSIZE)         ! HCl (gas phase) [ug/m3]
   REAL(dp) :: inSS(BLKSIZE)          ! sea salt [ug/m3]

   ! output from EQSAM

   REAL(dp) :: PNO3(BLKSIZE)          ! NO3 [ug/m3]
   REAL(dp) :: WH2O(BLKSIZE)          ! aerosol water [ug/m3]
   REAL(dp) :: PNH4(BLKSIZE)          ! NH4 [ug/m3]
   REAL(dp) :: GNH3(BLKSIZE)          ! NH3 (gas phase) [ug/m3]
   REAL(dp) :: GNO3(BLKSIZE)          ! HNO3 (gas phase) [ug/m3]
!   REAL(dp) :: GHCl(BLKSIZE)          ! HCl (gase phase) [ug/m3]

   ! ....................................................................

   ! *** AITKEN MODE ***

   DO LCELL = 1, NUMCELLS     !  loop on cells

      ! aerosol SO4
      inSO4(LCELL) = CBLK(LCELL,VSO4AI)

      ! aerosol NO3
      inNO3(LCELL) = CBLK(LCELL,VNO3AI)

      ! aerosol NH4
      inNH4(LCELL)  = CBLK(LCELL,VNH4AI)

      ! sea salt
      inSS(LCELL)   = CBLK(LCELL,VSEASI)

      ! gasphase concentrations
!      inH2SO4(LCELL) = CBLK(LCELL,VSULF)
      inNH3(LCELL)  = CBLK(LCELL,VNH3)
      inHNO3(LCELL) = CBLK(LCELL,VHNO3)
!      inHCl(LCELL)  = CBLK(LCELL,VHCl)

   END DO ! end loop over all cells

   CALL MADE_EQSAM( BLKTA, BLKRH, RH_HIST_AKN,                       &
                    inSO4, inNH4, inNH3, inNO3, inHNO3, inSS,        &
!                    inHCl,                                           &
                    PNH4, PNO3, WH2O, GNH3, GNO3,                    &
!                    GHCl,                                            &
                    blksize, numcells )

   ! update modal concs.

   DO LCELL=1,NUMCELLS
      CBLK(LCELL,VH2OAI) = WH2O(LCELL)
      CBLK(LCELL,VNH4AI) = PNH4(LCELL)
      CBLK(LCELL,VNO3AI) = PNO3(LCELL)

      CBLK(LCELL,VNH3)   = GNH3(LCELL)
      CBLK(LCELL,VHNO3)  = GNO3(LCELL)
!      CBLK(LCELL,VHCl)   = GHCl(LCELL)
   END DO

   ! *** ACCUMULATION MODE ***

   DO LCELL = 1, NUMCELLS     !  loop on cells

      ! aerosol SO4
      inSO4(LCELL)  = CBLK(LCELL,VSO4AJ)

      ! aerosol NO3
      inNO3(LCELL)  = CBLK(LCELL,VNO3AJ)

      ! aerosol NH4
      inNH4(LCELL)  = CBLK(LCELL,VNH4AJ)

      ! sea salt
      inSS(LCELL)   = CBLK(LCELL,VSEASJ)

      ! gasphase concentrations
!      inH2SO4(LCELL) = CBLK(LCELL,VSULF)
      inNH3(LCELL)  = CBLK(LCELL,VNH3)
      inHNO3(LCELL) = CBLK(LCELL,VHNO3)
!      inHCl(LCELL)  = CBLK(LCELL,VHCl)

   END DO ! end loop over all cells (=longitudes)

   CALL MADE_EQSAM( BLKTA, BLKRH, RH_HIST_ACC,                       &
                    inSO4, inNH4, inNH3, inNO3, inHNO3, inSS,        &
!                    inHCl,                                           &
                    PNH4, PNO3, WH2O, GNH3, GNO3,                    &
!                    GHCl,                                            &
                    blksize, numcells )

   ! update modal concs.

   DO LCELL=1,NUMCELLS
      CBLK(LCELL,VH2OAJ) = WH2O(LCELL)
      CBLK(LCELL,VNH4AJ) = PNH4(LCELL)
      CBLK(LCELL,VNO3AJ) = PNO3(LCELL)

      CBLK(LCELL,VNH3)   = GNH3(LCELL)
      CBLK(LCELL,VHNO3)  = GNO3(LCELL)
!      CBLK(LCELL,VHCl)   = GHCl(LCELL)
   END DO

   ! *** COARSE MODE ***

   DO LCELL = 1, NUMCELLS     !  loop on cells

      ! aerosol SO4 - not yet implemented
      inSO4(LCELL)  = 0.0_dp

      ! aerosol NO3 - not yet implemented
      inNO3(LCELL)  = 0.0_dp

      ! aerosol NH4 - not yet implemented
      inNH4(LCELL)  = 0.0_dp

      ! sea salt
      inSS(LCELL)   = CBLK(LCELL,VSEASC)

      ! gasphase concentrations
!      inH2SO4(LCELL) = 0.0_dp
      inNH3(LCELL)  = 0.0_dp ! - no NH4 implemented!
      inHNO3(LCELL) = 0.0_dp ! - no NO3 implemented!
!      inHCl(LCELL)  = CBLK(LCELL,VHCl)

   END DO ! end loop over all cells (=longitudes)

   CALL MADE_EQSAM( BLKTA, BLKRH, RH_HIST_COR,                       &
                    inSO4, inNH4, inNH3, inNO3, inHNO3, inSS,        &
!                    inHCl,                                           &
                    PNH4, PNO3, WH2O, GNH3, GNO3,                    &
!                    GHCl,                                            &
                    blksize, numcells )

   ! update modal concs.

   DO LCELL=1,NUMCELLS
      CBLK(LCELL,VH2OAC) = WH2O(LCELL)
!      CBLK(LCELL,VHCl)   = GHCl(LCELL)
   END DO

END SUBROUTINE MADE_EQL3X

!------------------------------------------------------------------------------!

SUBROUTINE MADE_EQSAM( temp, rh, rh_hist,                               &
                       inSO4, inNH4, inNH3, inNO3, inHNO3, inSS,        &
!                       inHCl,                                           &
                       PNH4, PNO3, WH2O, GNH3, GNO3,                    &
!                       GHCl,                                            &
                       blksize, numcells )

   !_______________________________________________________________________
   !      Written by Swen Metzger 3/11/99. Modified 2002, 2003.
   !      Modified for use with MADE by Axel Lauer, DLR (2004).
   !
   !      Department of Atmospheric Chemistry, Max-Planck-Institute
   !      for Chemistry.
   !
   !      email: metzger@mpch-mainz.mpg.de
   !      http://www.mpch-mainz.mpg.de/~metzger
   !
   !      COPYRIGHT 1999-2003
   !
   !      purpose
   !      -------
   !      EQSAM (EQuilibrium Simplified Aerosol Model) is a new and
   !      computationally efficient thermodynamic aerosol composition
   !      model that allows to calculate the gas/aerosol equilibrium
   !      partitioning, including aerosol water, sufficiently fast and
   !      accurate for global (or even regional) modeling. EQSAM is based
   !      on a number of parameterizations, including single solute
   !      molalities and activity coefficients (AC). The thermodynamic
   !      framework (domains and subdomains, internally mixed aerosols) 
   !      is the same as of more sophisticated thermodynamic equilibrium
   !      models (EQMs), e.g. of ISORROPIA (Nenes et al., 1998). Details
   !      are given in the references below (and the references therein).
   !
   !      The main assumption on which EQSAM/EQMs are based is
   !      thermodynamical and chemical equilibrium. From this assumption
   !      it directly follows that the aerosol water activity (aw) equals
   !      the ambient relative humidity (RH), if the water vapor pressure
   !      is sufficiently larger than the partial vapor pressure of the
   !      aerosol compounds. This is approximately true for tropospheric
   !      aerosols. Given the large amount of water vapor present, water
   !      vapor and aerosol water equilibrate relatively faster compared
   !      to all other aerosol compounds. This is subsequently also true
   !      for single aerosol compounds. The water activity of single
   !      solutes must also equal RH under this assumption. Therefore, the
   !      so called ZSR-relation is (and can be) used to calculate the
   !      aerosol associated water mass (simply from the sum of all water
   !      mass fractions that are derived from measured single solute
   !      molalities).
   !
   !      In contrast to other EQMs, EQSAM utilizes the fact that the RH
   !      fixes the water activity (under the above assumptions) and the
   !      consequence that any changes in RH also causes changes in the
   !      aerosol water mass and, hence, aerosol activity (including
   !      activity coefficients). Thus, an decrease (increase) in RH
   !      decrease (increases) the aerosol water mass (and water activity).
   !      This can change the aerosol composition, e.g. due to
   !      condensation (evaporation/crystallization), because the vapor
   !      pressure above the aerosol reduces (increases). In turn, a vapor
   !      pressure reduction (increase) due to changes in the aerosol
   !      composition is compensated by an associated condensation
   !      (evaporation) of water vapor to maintain the aerosol molality to
   !      remain constant (because aw=RH). Furthermore, the aerosol water 
   !      mainly depends on the aerosol mass and the type of solute, so
   !      that parameterizations of single solute molalities and activity
   !      coefficients can be defined, only depending on the type of
   !      solute and RH. The advantage of using such parameterizations is
   !      that the entire aerosol equilibrium composition can be solved
   !      analytically, i.e. non-iteratively, which considerably reduces
   !      the amount of CPU time that is usually need for aerosol
   !      thermodynamic calculations (especially if an EQM is incorporated
   !      in an aerosol dynamical model that is in turn embedded in a high
   !      resolution regional or global model).
   !
   !      However, EQSAM should still be regarded as a starting point for
   !      further developments. There is still room for improvements. For
   !      instance, this code is not yet numerically optimized (vectorized)
   !      and a number of improvements with respect to an explicit
   !      treatment of additional equilibrium reactions, missing (or only
   !      implicit) dissociation, and a basic parameterization of the
   !      water uptake.
   !
   !      Note that EQSAM was originally developed to calculate the
   !      gas/aerosol equilibrium partitioning of the
   !      ammonium-sulfate-nitrate-water system for climate models,
   !      excluding solid compounds. This version (eqsam_v03d.f90) is
   !      extended with respect to sea salt. Solids/hysteresis are treated
   !      in a simplified manner. Results of a box model comparison with
   !      ISORROPIA will be available from the web page. Please also note
   !      that the water uptake is based on additional (unpublished)
   !      parameterizations for single solute molalities, which are
   !      derived from tabulated measurements used in ISORROPIA. Note
   !      further that this extended version (eqsam_v03d.f90) is not yet
   !      published. A publication is in progress.
   !
   ! ToDo:
   !     Split ion-pairs into ions for water parameterizations (since info
   !     is actually available)
   !     Include uptake/dissociation of NH3, HNO3, HCl (mainly to get pH
   !     right at near neutral conditions)
   !     Extension to K+,Ca++,Mg++, CO2/(CO3)2--/HCO3-,SOA,etc..
   !     (maybe not)
   !     Vectorization. Translation of hardcoded formulas in array syntax.
   !     I/O Interface and program structure clean up.
   !     EQSAM info webpage.
   !
   ! Version History:
   !
   !  eqsam_v03d.f90 (MPI-CH, June 2003):
   !   - gama parameterizations now according to Metzger 2002 (JGR
   !     Appendix)
   !   - improved pH calculations (still restricted to strong acids)
   !   - removed bug that lead to too high nitrate formation at dry and
   !     cold regions (UT/LS) 
   !   - removed bug in solid/hysteresis calculations 
   !     (both bugs introduced in eqsam_v03b.f90 by cleaning up
   !     eqsam_v02a.f90)
   !   
   !  eqsam_v03c.f90 (MPI-CH, April 2003):
   !   - more accurate paramterizations of single solute molalities
   !     (Na, Cl species)
   !   - cleanded up RHD subdomain structure
   !   - improved water uptake (Na, Cl species)
   !
   !  eqsam_v03b.f90 (MPI-CH, March 2003):
   !                 System extended to HCl,Cl-/Na+.
   !                 Parameterization (fit) of additional HNO3 uptake
   !                 removed. Instead, complete analytical solution of
   !                 equilibrium reactions, based on the AC-RH
   !                 relationship.
   !  eqsam_v03.f90  (IMAU, October 1999):
   !                 Test version (included in TM3).
   !  eqsam_v02a.f90 (IMAU, April 2000):
   !                 Box model version.
   !  eqsam_v02.f90  (IMAU, October 1999):
   !                 TM3 version.
   !                 Version including solids and additional HNO3 uptake
   !                 on acidic aerosols (parameterized).
   !  eqsam_v01b.f90 (MPI-CH, January 2003):
   !                 Same as eqsam_v01a.f90 (additional lines though
   !                 uncommented for test purposes only).
   !  eqsam_v01a.f90 (IMAU, April 2000):
   !                 Box model version.
   !  eqsam_v01.f90  (IMAU, October 1999):
   !                 TM3 version.
   !                 First and most basic version (without solids) for
   !                 better vectorization (for global modeling).
   !                 System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, H2O
   !                 based on equilibrium / internal mixture assumption /
   !                 aw=rh / ZSR-relation parameterization of activcity
   !                 coefficients (AC), i.e. an AC-RH relationship
   !
   !
   !      interface
   !      ---------
   !      call made_eqsam( temp, rh, rh_hist, inSO4, inNH4, inNH3, &
   !                       inNO3, inHNO3, inSS, PNH4, PNO3,        &
   !                       WH2O, GNH3, GNO3,                       &
   !                       numcells )
   !
   !      temp     = temperature                                  [K]
   !      rh       = relativ humidity                             [0-1]
   !      rh_hist  = old relative humidity to calculate aerosol hysteresis
   !                 1: dry history (--> solids)
   !                 2: wet history
   !      inSO4    = input  sulfate mass (w/o seasalt-so4)        [ug/m3]
   !      inNH4    = input  ammonium mass                         [ug/m3]
   !      inNH3    = input  ammonia mass (gas phase)              [ug/m3]
   !      inNO3    = input  nitrate mass                          [ug/m3]
   !      inHNO3   = input  nitric acid (gas phase)               [ug/m3]
   !      inSS     = input  seasalt mass                          [ug/m3]
!   !      inHCl    = input  hydrochloric acid (gas phase)         [ug/m3]
   !      PNH4     = output ammonium mass                         [ug/m3]
   !      PNO3     = output nitrate mass                          [ug/m3]
   !      WH2O     = output aerosol water                         [ug/m3]
   !      GNH3     = output ammonia (gas phase)                   [ug/m3]
   !      GNO3     = output nitric acid (gas phase)               [ug/m3]
!   !      GHCl     = output hydrochloric acid (gas phase)         [ug/m3]
   !      blksize  = size of arrays
   !      numcells = actual number of cells in arrays
   !
   !      method
   !      ------
   !      equilibrium / internal mixture assumption / aw=rh
   !      System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, HCl,Cl-/Na+, H2O
   !              (K+,Ca++,Mg++)
   !      external
   !      --------
   !
   !      references
   !      ----------
   !      Swen Metzger Ph.D Thesis, University Utrecht, 2000.
   !       http://www.library.uu.nl/digiarchief/dip/diss/1930853/inhoud.htm
   !
   !      Metzger, S. M., F. J. Dentener, J. Lelieveld, and S. N. Pandis,
   !       GAS/AEROSOL PARTITIONING I: A COMPUTATIONALLY EFFICIENT MODEL,
   !       J Geophys. Res., 107, D16, 10.1029/2001JD001102, 2002
   !       http://www.agu.org/journals/jd/jd0216/2001JD001102/index.html
   !      Metzger, S. M., F. J. Dentener, A. Jeuken, and M. Krol,
   !       J. Lelieveld, GAS/AEROSOL PARTITIONING II: GLOBAL MODELING
   !       RESULTS, J Geophys. Res., 107, D16, 10.1029/2001JD001103, 2002.
   !       http://www.agu.org/journals/jd/jd0216/2001JD001103/index.html
   !_______________________________________________________________________

   implicit none
   intrinsic abs, exp, log, max, min, nint, real, sqrt

   real(dp), parameter :: RH_HIST_DW = 1.50_dp ! mean value for mixture of wet
                                               ! (2) and dry (1) gridboxes
                                               ! (needed for HYSTERESIS)
   real(dp), parameter :: TT0 = 298.15_dp      ! = T0 in EQSAM, avoid conflicts
                                               ! with global parameter T0
!   real(dp), parameter :: TT1 = 298.0_dp
   real(dp), parameter :: R  = 82.0567e-6_dp  ! in cu.m*atm/deg/mole

   real(dp), parameter :: RHMAX = 0.99_dp     ! restrict to max / min RH
   real(dp), parameter :: RHMIN = 0.0001_dp

   real(dp), parameter :: MWH20 = 55.51_dp*18.01_dp ! "H-2-Null"
   real(dp), parameter :: ZERO  = 0.0_dp

   real(dp), parameter :: ZEPS  = 1.0e-19_dp  ! epsilon

   real(dp) :: dum                            ! auxiliary variable

   ! exponents of AC-RH functions
   real(dp), parameter :: GF1 = 0.25_dp
   real(dp), parameter :: K = 2._dp
   !______________________________________________
   integer, parameter :: NPAIR = 10

   integer :: il
   integer :: IFLAG

   real(dp) :: RHL
   real(dp) :: T0T,TT,RHD,KAN,KAC,GAMA,GG,GF,GFN
   real(dp) :: X00,X01,X02,X03,X04,X05,X08,X09,X10,X11
   real(dp) :: X0,X1,X2,X3,X4,X5,X6,XK10,XK6
   real(dp) :: ZFLAG,ZKAN,ZKAC,COEF
   real(dp) :: TNH4,TSO4,TNO3,TNa,TCl,TPo,TCa,TMg
   real(dp) :: GSO4,PCl,PNa,PSO4
   real(dp) :: ASO4,ANO3,ANH4,ACl,ANa,SNH4,SSO4,SNO3,SCl,SNa

   ! input for EQSAM

   integer, intent(in)     :: blksize            ! dimension of arrays
   integer, intent(in)     :: numcells           ! number of fields in arrays

   real(dp), intent(in)    :: temp(blksize)     ! temperature [K]
   real(dp), intent(in)    :: rh(blksize)       ! relative humidity (0-1)
   real(dp), intent(inout) :: rh_hist(blksize)  ! history of rel. humidity
                                                 ! --> hysteresis
   real(dp), intent(in)    :: inSO4(blksize)    ! SO4(p) (w/o seasalt-so4)
                                                 ! [ug/m3]
   real(dp), intent(in)    :: inNH4(blksize)    ! NH4(p)  [ug/m3]
   real(dp), intent(in)    :: inNH3(blksize)    ! NH3(g)  [ug/m3]
   real(dp), intent(in)    :: inNO3(blksize)    ! NO3(p)  [ug/m3]
   real(dp), intent(in)    :: inHNO3(blksize)   ! HNO3(g) [ug/m3]
   real(dp), intent(in)    :: inSS(blksize)     ! Seasalt [ug/m3]
!   real(dp), intent(in)    :: inHCl(blksize)    ! HCl(g)  [ug/m3]

   ! output from EQSAM

   real(dp), intent(out)   :: PNH4(blksize)     ! NH4(p)  [ug/m3]
   real(dp), intent(out)   :: PNO3(blksize)     ! NO3(p)  [ug/m3]
   real(dp), intent(out)   :: WH2O(blksize)     ! H2O(p) total [ug/m3]
   real(dp), intent(out)   :: GNH3(blksize)     ! NH3(g)  [ug/m3]
   real(dp), intent(out)   :: GNO3(blksize)     ! HNO3(g) [ug/m3]
!   real(dp), intent(out)   :: GHCl(blksize)     ! HCl(g)  [ug/m3]

   !_______________________________________________

   real(dp) :: w1(8), w2(8)
   ! RHD / MRHD arrays for different aerosol types
   real(dp) :: RHDA(8),RHDE(8),RHDX(8),RHDZ(8)
   ! arrays of ion pairs
   real(dp) :: M0(NPAIR),MW(NPAIR),NW(NPAIR),ZW(NPAIR)
   !
   ! salt solutes:
   !   1 = NACl,    2 = (NA)2SO4,     3 = NANO3,  4 = (NH4)2SO4,
   !   5 = NH4NO3,  6 = NH4CL,        7 = 2H-SO4, 8 = NH4HSO4,
   !   9 = NAHSO4, 10 = (NH4)3H(SO4)2
   !
   ! mole mass of the salt solute
   DATA MW / 58.5_dp, 142.0_dp,  88.0_dp, 132.0_dp, 80.0_dp, 53.5_dp, 98.0_dp,&
            115.0_dp, 120.0_dp, 247.0_dp /
   ! square of max. dissocation number (not consistent)
   DATA NW / 2.0_dp, 2.5_dp, 2.5_dp, 2.5_dp, 3.5_dp, 1.0_dp, 4.5_dp, 2.0_dp,  &
             2.0_dp, 2.5_dp /
   ! exponents of water activity functions
   DATA ZW / 0.67_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, &
             1.0_dp,  1.0_dp /
   !
   ! RHD / MRHD values as of ISORROPIA / SCAPE (T=298.15K)
   DATA RHDA / 0.32840_dp, 0.4906_dp, 0.6183_dp, 0.7997_dp, 0.67500_dp,       &
               0.5000_dp, 0.4000_dp, 0.0000_dp /
   ! Temp. coeff.
   DATA RHDE / -1860.0_dp, -431.0_dp, 852.00_dp, 80.000_dp, 262.000_dp,       &
                3951.0_dp,  384.00_dp,  0.0000_dp /
   !_______________________________________________________________________

   do il = 1, numcells

      TT = temp(il)           ! temperature [K]
      RHL = RH(il)            ! relative humidity [0-1]
      !______________________________________________

!      zflag=1._dp
      iflag=1

      ! Na+ (ss + xsod) (a) [mol/m^3 air]
      w1(1) = 1.0e-6_dp * (0.3061_dp * inSS(il) / MWNa)
      ! H2SO4    + SO4-- (p)  [mol/m^3 air]
      w1(2) = 1.0e-6_dp * ((inSO4(il) + 0.0768_dp * inSS(il)) / MWSO4)
      ! NH3 (g) + NH4+ (p)  [mol/m^3 air]
      w1(3) = 1.0e-6_dp * (inNH3(il) / MWNH3 + inNH4(il) / MWNH4)
      ! HNO3 (g) + NO3- (p) [umol/m^3 air]
      w1(4) = 1.0e-6_dp * (inNO3(il) / MWNO3 + inHNO3(il) / MWHNO3)
      ! HCl(g)+Cl-(p)[mol/m^3 air]
!      w1(5) = 1.0e-6_dp * (inHCl(il) / MWHCl + 0.5504_dp * inSS(il) / MWCl)
      w1(5) = 1.0e-6_dp * 0.5504_dp * inSS(il) / MWCl
      ! K+ (p)       [mol/m^3 air]
      w1(6) = 1.0e-6_dp * (0.0110_dp * inSS(il) / MWK)
      ! Ca++ (p)     [mol/m^3 air]
      w1(7) = 1.0e-6_dp * (0.0116_dp * inSS(il) / MWCa)
      ! Mg++ (p)     [mol/m^3 air]
      w1(8) = 1.0e-6_dp * (0.0369_dp * inSS(il) / MWMg)

      TNa   = w1(1)           ! total input sodium   (g+p) 
      TSO4  = w1(2)           ! total input sulfate  (g+p) 
      TNH4  = w1(3)           ! total input ammonium (g+p)
      TNO3  = w1(4)           ! total input nitrate  (g+p) 
      TCl   = w1(5)           ! total input chloride (g+p) 
      TPo   = w1(6)           ! total input potasium (g+p) 
      TCa   = w1(7)           ! total input calcium  (g+p)
      TMg   = w1(8)           ! total input magnesium(g+p)

      ! SULFATE RICH

      if ((w1(1)+w1(3)+w1(6)+2._dp*(w1(7)+w1(8))) <= (2._dp*w1(2))) &
         iflag=3 !zflag=3._dp

      ! SULFATE VERY RICH CASE if (NH4+Na+K+2(Ca+Mg))/SO4 < 1

      if ((w1(1)+w1(3)+w1(6)+2._dp*(w1(7)+w1(8))) <= w1(2))         &
         iflag=4 !zflag=4._dp

      ! SULFATE NEUTRAL CASE

      if ((w1(1)+w1(3)+w1(6)+2._dp*(w1(7)+w1(8))) > (2._dp*w1(2)))  &
         iflag=2 !zflag=2._dp

      ! SULFATE POOR AND CATION POOR CASE

      if ((w1(1)+w1(6)+2._dp*(w1(7)+w1(8))) > (2._dp*w1(2)))        &
         iflag=1 !zflag=1._dp

      RHL = MAX(RHMIN, RHL)
      RHL = MIN(RHL, RHMAX)

      ! CALCULATE TEMPERATURE DEPENDENCY FOR SOME RHDs

      DUM = (1._dp/TT-1._dp/TT0)

      RHDX(1)=RHDA(1)*exp(RHDE(1)*DUM)
      RHDZ(1)=RHDX(1)
      RHDX(2)=RHDA(2)*exp(RHDE(2)*DUM)
      RHDZ(2)=RHDX(2)
      RHDX(3)=RHDA(3)*exp(RHDE(3)*DUM)
      RHDZ(3)=RHDX(3)
      RHDX(4)=RHDA(4)*exp(RHDE(4)*DUM)
      RHDZ(4)=RHDX(4)
      RHDX(5)=RHDA(5)*exp(RHDE(5)*DUM)
      RHDZ(5)=RHDX(5)
      RHDX(6)=RHDA(6)*exp(RHDE(6)*DUM)
      RHDZ(6)=RHDX(6)
      RHDX(7)=RHDA(7)*exp(RHDE(7)*DUM)
      RHDZ(7)=RHDX(7)
      RHDX(8)=RHDA(8)*exp(RHDE(8)*DUM)
      RHDZ(8)=RHDX(8)
      
      ! ACCOUNT FOR VARIOUS AMMOMIUM/SODIUM SULFATE SALTS ACCORDING
      ! TO MEAN VALUE AS OF ISORROPIA

      GG=2.0_dp    ! (Na)2SO4 / (NH4)2SO4 IS THE PREFFERED SPECIES
                   ! FOR SULFATE DEFICIENT CASES
!      IF (ZFLAG == 3._dp) THEN
      IF (IFLAG == 3) THEN
         IF (RHL <= RHDZ(7)) THEN      ! ACCOUNT FOR MIXTURE OF
                                       ! (NH4)2SO4(s) & NH4HSO4(s) &
                                       ! (NH4)3H(SO4)2(s) 
             GG=1.677_dp               ! (Na)2SO4 & NaHSO4
             !GG=1.5
          ELSEIF (RHL > RHDZ(7).AND.RH(il) <= RHDZ(5)) THEN
             ! MAINLY (Na)2SO4 / (NH4)2SO4(s) & (NH4)3H(SO4)2(s)
             GG=1.75_dp
             !GG=1.5
          ELSEIF (RHL >= RHDZ(5)) THEN ! (NH4)2SO4(S) & NH4HSO4(S) &
                                       ! SO4-- & HSO4-
             GG=1.5_dp                    ! (Na)2SO4 & NaHSO4
          ENDIF
      ENDIF

!      IF (ZFLAG == 4._dp) GG=1.0_dp    ! IF SO4 NEUTRALIZED, THEN ONLY
      IF (IFLAG == 4) GG=1.0_dp        ! IF SO4 NEUTRALIZED, THEN ONLY
                                       ! AS NaHSO4 / NH4HSO4(S) OR
                                       ! HSO4- / H2SO4
      RHD=RHL

      IF (RH_HIST(il) < RH_HIST_DW) THEN ! GET RHD FOR
                                         ! SOLIDS/HYSTERESIS
         !
         ! GET LOWEST DELIQUESCENCE RELATIVE HUMIDITIES ACCORDING
         ! TO THE CONCENTRATION DOMAIN (APROXIMATION) BASED ON
         ! RHD / MRHD ISORROPIA/SCAPE
         !
         w2(1)=1._dp
         w2(2)=1._dp
         w2(3)=1._dp
         w2(4)=1._dp
         w2(5)=1._dp
         w2(6)=1._dp
         w2(7)=1._dp
         w2(8)=1._dp

         ! skip compound in RHD calculation if value is concentration
         ! is zero or rather small

         if (w1(1) <= 1.e-12_dp) w2(1)=0._dp
         if (w1(2) <= 1.e-12_dp) w2(2)=0._dp
         if (w1(3) <= 1.e-12_dp) w2(3)=0._dp
         if (w1(4) <= 1.e-12_dp) w2(4)=0._dp
         if (w1(5) <= 1.e-12_dp) w2(5)=0._dp
         if (w1(6) <= 1.e-12_dp) w2(6)=0._dp
         if (w1(7) <= 1.e-12_dp) w2(7)=0._dp
         if (w1(8) <= 1.e-12_dp) w2(8)=0._dp

         ! GET LOWEST RHD ACCORDING TO THE CONCENTRATION DOMAIN

         ! iflag=1 (cation rich)  ...
         ! 1. sea salt      aerosol          : RHDX(1)=MgCl2
         ! 2. mineral dust  aerosol          : RHDX(2)=Ca(NO3)2
         !
         ! iflag=2 (sulfate neutral) ...
         ! 3. ammonium + nitrate             : RHDX(3)= NH4NO3
         ! 4. ammonium + sulfate             : RHDX(4)=(NH4)2SO4
         ! 5. ammonium + sulfate mixed salt  : RHDX(5)=
         !                                   (NH4)3H(SO4)2, (NH4)2SO4
         ! 6. ammonium + nitrate  + sulfate  : RHDX(6)=(NH4)2SO4,
         !                                   NH4NO3,NA2SO4, NH4Cl
         !
         ! iflag=3 (sulfate rich) ...
         ! 7. ammonium + sulfate  (1:1,1.5)  : RHDX(7)= NH4HSO4
         !
         ! iflag=4 (sulfate very rich) ...
         ! 8. sulfuric acid                  : RHDX(8)= H2SO4

!         IF (ZFLAG == 1._dp) THEN
         IF (IFLAG == 1) THEN

            RHD=W2(1)+W2(5)              ! Na+  dependency
!            IF(RHD == 0._dp) RHDX(1)=1._dp 
            IF(RHD < ZEPS) RHDX(1)=1._dp 
            RHD=W2(6)+W2(7)+W2(8)        ! K+/Ca++/Mg++ dependency (incl. ss)
!            IF(RHD == 0._dp) RHDX(2)=1._dp
            IF(RHD < ZEPS) RHDX(2)=1._dp

            RHD=MIN(RHDX(1),RHDX(2))

!         ELSEIF (ZFLAG == 2._dp) THEN
         ELSEIF (IFLAG == 2) THEN

            RHD=W2(3)*W2(4)              ! NH4+ & NO3- dependency
!            IF(RHD == 0._dp) RHDX(3)=1._dp 
            IF(RHD < ZEPS) RHDX(3)=1._dp 
            RHD=W2(2)+W2(3)              ! NH4+ & SO4-- dependency
!            IF(GG  /= 2._dp) RHD=0._dp   ! account only for pure (NH4)2SO4
            IF(ABS(GG - 2._dp) > ZEPS) RHD=0._dp  ! account only for pure
                                                  ! (NH4)2SO4
!            IF(RHD == 0._dp) RHDX(4)=1._dp
            IF(RHD < ZEPS) RHDX(4)=1._dp
            RHD=W2(2)+W2(3)              ! NH4+ & SO4-- dependency
!            IF(RHD == 0._dp) RHDX(5)=1._dp
            IF(RHD < ZEPS) RHDX(5)=1._dp
            RHD=W2(2)+W2(3)+W2(4)+W2(5)  ! (NH4)2SO4, NH4NO3, Na2SO4
                                         ! NH4Cl dependency
!            IF(RHD == 0._dp) RHDX(6)=1._dp
            IF(RHD < ZEPS) RHDX(6)=1._dp

            RHD=MIN(RHDX(3),RHDX(4),RHDX(5),RHDX(6))

!         ELSEIF (ZFLAG == 3._dp) THEN
         ELSEIF (IFLAG == 3) THEN

            RHD=W2(2)+W2(3)              ! NH4+ & SO4-- dependency
!            IF(RHD == 0._dp) RHDX(7)=1._dp
            IF(RHD < ZEPS) RHDX(7)=1._dp
            RHD=RHDX(7)

!         ELSEIF (ZFLAG == 4._dp) THEN
         ELSEIF (IFLAG == 4) THEN

            RHD=W2(2)                    ! H2SO4 dependency (assume no dry
                                         ! aerosol)
!            IF(RHD == 0._dp) RHDX(8)=1._dp
            IF(RHD < ZEPS) RHDX(8)=1._dp

            RHD=RHDX(8)

         ENDIF ! IFLAG
      ENDIF ! SOLIDS

      ! GET WATER ACTIVITIES ACCORDING TO METZGER, 2000.
      ! FUNCTION DERIVED FROM ZSR RELATIONSHIP DATA (AS USED IN ISORROPIA)

      DUM = (1._dp/RHL-1._dp)
      M0(1)  = ((NW(1) *MWH20/MW(1) *DUM))**ZW(1)
      M0(2)  = ((NW(2) *MWH20/MW(2) *DUM))!**ZW(2)
      M0(3)  = ((NW(3) *MWH20/MW(3) *DUM))!**ZW(3)
      M0(4)  = ((NW(4) *MWH20/MW(4) *DUM))!**ZW(4)
      M0(5)  = ((NW(5) *MWH20/MW(5) *DUM))!**ZW(5)
      M0(6)  = ((NW(6) *MWH20/MW(6) *DUM))!**ZW(6)
      M0(7)  = ((NW(7) *MWH20/MW(7) *DUM))**ZW(7)
      M0(8)  = ((NW(8) *MWH20/MW(8) *DUM))!**ZW(8)
      M0(9)  = ((NW(9) *MWH20/MW(9) *DUM))!**ZW(9)
      M0(10) = ((NW(10)*MWH20/MW(10)*DUM))!**ZW(10)

      ! CALCULATE TEMPERATURE DEPENDENT EQUILIBRIUM CONSTANTS

      T0T=TT0/TT
      COEF=1.0_dp+LOG(T0T)-T0T

      ! EQUILIBRIUM CONSTANT NH4NO3(s) <==> NH3(g) + HNO3(g)
      ! [atm^2] (ISORROPIA)

      XK10 = 5.746e-17_dp
      XK10= XK10 * EXP(-74.38_dp*(T0T-1.0_dp) + 6.120_dp*COEF)
      KAN = XK10/(R*TT)/(R*TT)

      ! EQUILIBRIUM CONSTANT  NH4Cl(s) <==> NH3(g) + HCl(g)
      ! [atm^2] (ISORROPIA)

      XK6  = 1.086e-16_dp
      XK6 = XK6 * EXP(-71.00*(T0T-1.0_dp) + 2.400_dp*COEF)
      KAC = XK6/(R*TT)/(R*TT)

      ! GET MEAN MOLAL IONIC ACTIVITY COEFF ACCORDING TO METZGER, 2002.

      GAMA=0.0_dp
      ZFLAG=REAL(IFLAG)
      IF (RHL >= RHD) GAMA=(RHL**IFLAG/(1000._dp/ZFLAG*(1._dp-RHL)+ZFLAG))
      !                          ^^^^^
      ! using IFLAG instead of ZFLAG might give better performance
      GAMA = GAMA**GF1      ! ONLY GAMA TYPE OF NH4NO3, NaCl, etc.
                            ! NEEDED SO FAR

      GAMA=0.0_dp
      GFN=K*K               ! K=2, i.e. condensation of 2 water
                            ! molecules per 1 mole ion pair
      GF=GFN*GF1            ! = GFN[=Nw=4] * GF1[=(1*1^1+1*1^1)/2/Nw=1/4] = 1
                            ! ONLY GAMA TYPE OF NH4NO3, NH4Cl,
                            ! etc. needed so far

      IF (RHL >= RHD) GAMA=RHL**GF/((GFN*MWH20*(1._dp/RHL-1._dp)))**GF1

      GAMA = MIN(GAMA,1.0_dp)  ! FOCUS ON 0-1 SCALE
      GAMA = MAX(GAMA,0.0_dp)
      GAMA = (1._dp-GAMA)**K   ! transplate into aqueous phase
                            ! equillibrium and account for 
                            ! enhanced uptake of aerosol precursor
                            ! gases with increasing RH
                            ! (to match the results of ISORROPIA)

      ! CALCULATE RHD DEPENDENT EQ: IF RH <  RHD =>
      ! NH4NO3(s) <==> NH3 (g) + HNO3(g) (ISORROPIA)
      ! IF RH >> RHD => HNO3  (g)   -> NO3 (aq)

      X00  = MAX(ZERO,MIN(TNa,GG*TSO4))       ! MAX SODIUM SULFATE
      X0   = MAX(ZERO,MIN(TNH4,GG*TSO4-X00))  ! MAX AMMOMIUM SULFATE
      X01  = MAX(ZERO,MIN(TNa-X00, TNO3))     ! MAX SODIUM NITRATE
      X1   = MAX(ZERO,MIN(TNH4-X0,TNO3-X01))  ! MAX AMMOMIUM NITRATE
      !
      X02  = MAX(ZERO,MIN(TNa-X01-X00,TCl))   ! MAX SODIUM CHLORIDE
      X03  = MAX(ZERO,MIN(TNH4-X0-X1,TCl-X02))! MAX AMMOMIUM CHLORIDE

      X2   = MAX(TNH4-X1-X0-X03,ZERO)        ! INTERIM RESIDUAL NH3
      X3   = MAX(TNO3-X1-X01,ZERO)           ! INTERIM RESIDUAL HNO3
      X04  = MAX(TSO4-(X0+X00)/GG,ZERO)      ! INTERIM RESIDUAL H2SO4
      X05  = MAX(TCl-X03-X02,ZERO)           ! INTERIM RESIDUAL HCl
!      X06  = MAX(TNa-X02-X01-X00,ZERO)       ! INTERIM RESIDUAL Na
      ! (should be zero for electro-neutrality in input data)
      !
      ZKAN=2._dp
      IF (RHL >= RHD) ZKAN=ZKAN*GAMA

      X4   = X2 + X3
      X5   = SQRT(X4*X4+KAN*ZKAN*ZKAN)
      X6   = 0.5_dp*(-X4+X5)
      X6   = MIN(X1,X6)
      
!      GHCl(il) = X05                   ! INTERIM RESIDUAl HCl
      GNH3(il) = X2 + X6               ! INTERIM RESIDUAl NH3
      GNO3(il) = X3 + X6               ! RESIDUAl HNO3
      GSO4 = X04                       ! RESIDUAl H2SO4
      PNa  = X02 + X01 + X00           ! RESIDUAl Na (neutralized)
      
      ZKAC=2._dp
      IF(RHL >= RHD) ZKAC=ZKAC*GAMA

!      X08   = GNH3(il) + GHCl(il)
      X08   = GNH3(il)
      X09   = SQRT(X08*X08+KAC*ZKAC*ZKAC)
      X10   = 0.5_dp*(-X08+X09)
      X11   = MIN(X03,X10)

!      GHCl(il) = GHCl(il) + X11        ! RESIDUAL HCl
      GNH3(il) = GNH3(il) + X11        ! RESIDUAL NH3

      ! GO SAVE ...

!      IF (GHCl(il) < 0._dp) GHCl(il)=0._dp
      IF (GSO4     < 0._dp) GSO4=0._dp
      IF (GNH3(il) < 0._dp) GNH3(il)=0._dp
      IF (GNO3(il) < 0._dp) GNO3(il)=0._dp
      IF (PNa      < 0._dp) PNa=0._dp
      IF (GSO4     > TSO4)  GSO4=TSO4
      IF (GNH3(il) > TNH4)  GNH3(il)=TNH4
      IF (GNO3(il) > TNO3)  GNO3(il)=TNO3
!      IF (GHCl(il) > TCl)   GHCl(il)=TCl
      IF (PNa      > TNa)   PNa=TNa

      ! DEFINE AQUEOUSE PHASE (NO SOLID NH4NO3 IF NO3/SO4>1,
      ! TEN BRINK, ET AL., 1996, ATMOS ENV, 24, 4251-4261)

!      IF (RH_HIST(il) == 1._dp.AND.RHL < RHD) THEN  ! SOLIDS/HYSTERESIS
      IF ((NINT(RH_HIST(il)) == 1).AND.(RHL < RHD)) THEN  ! SOLIDS/HYSTERESIS

         ! EVERYTHING DRY, ONLY H2SO4 (GSO4) REMAINS IN THE AQUEOUSE PHASE

         ANH4 = 0._dp
         ASO4 = 0._dp
         ANO3 = 0._dp
         ACl  = 0._dp
         ANa  = 0._dp

      ELSE  !  SUPERSATURATED SOLUTIONS NO SOLID FORMATION

         ASO4 = TSO4 - GSO4
         ANH4 = TNH4 - GNH3(il)
         ANO3 = TNO3 - GNO3(il)
!         ACl  = TCl  - GHCl(il)
         ACl  = TCl
         ANa  = PNa

      ENDIF ! SOLIDS/HYSTERESIS

      ! CALCULATE AEROSOL WATER [kg/m^3(air)]
      !
      ! salt solutes:
      !   1 = NaCl,    2 = (Na)2SO4, 3 = NaNO3,  4 = (NH4)2SO4,
      !   5 = NH4NO3,  6 = NH4Cl,    7 = 2H-SO4, 8 = NH4HSO4,
      !   9 = NaHSO4, 10 = (NH4)3H(SO4)2
      !
!      IF (ZFLAG == 1._dp) THEN
      IF (IFLAG == 1) THEN
         WH2O(il)    = ASO4/M0(2) + ANO3/M0(3) + ACl/M0(6)
!      ELSE IF (ZFLAG == 2._dp) THEN
      ELSE IF (IFLAG == 2) THEN
         WH2O(il)    = ASO4/M0(9) + ANO3/M0(5) + ACl/M0(6)
!      ELSE IF (ZFLAG == 3._dp) THEN
      ELSE IF (IFLAG == 3) THEN
         WH2O(il)    = ASO4/M0(8) + ANO3/M0(5) + ACl/M0(6)
!      ELSE IF (ZFLAG == 4._dp) THEN
      ELSE IF (IFLAG == 4) THEN
         WH2O(il)    = ASO4/M0(8) + GSO4/M0(7)
      END IF

      !
      !-------------------------------------------------------
      ! calculate diagnostic output consistent with other EQMs ...
      !
      ASO4 = ASO4 + GSO4  ! assuming H2SO4 remains aqueous

      TNa   = TNa  * 1.e6_dp ! total input sodium   (g+p)  [umol/m^3]
      TSO4  = TSO4 * 1.e6_dp ! total input sulfate  (g+p)  [umol/m^3]
      TNH4  = TNH4 * 1.e6_dp ! total input ammonium (g+p)  [umol/m^3]
      TNO3  = TNO3 * 1.e6_dp ! total input nitrate  (g+p)  [umol/m^3]
      TCl   = TCl  * 1.e6_dp ! total input chloride (g+p)  [umol/m^3]
      TPo   = TPo  * 1.e6_dp ! total input potasium (g+p)  [umol/m^3]
      TCa   = TCa  * 1.e6_dp ! total input calcium  (g+p)  [umol/m^3]
      TMg   = TMg  * 1.e6_dp ! total input magnesium(g+p)  [umol/m^3]
      !
      ! residual gas:
      GNH3(il) = GNH3(il) * 1.e6_dp  ! residual NH3 [ug/m3]
      GSO4     = GSO4     * 1.e6_dp  ! residual H2SO4
      GNO3(il) = GNO3(il) * 1.e6_dp  ! residual HNO3
!      GHCl(il) = GHCl(il) * 1.e6_dp  ! residual HCl

      ! total particulate matter (neutralized)
      PNH4(il)=TNH4-GNH3(il)      ! particulate ammonium   [umol/m^3]
      PNO3(il)=TNO3-GNO3(il)      ! particulate nitrate    [umol/m^3]
!      PCl     =TCl -GHCl(il)      ! particulate chloride   [umol/m^3]
      PCl     =TCl                ! particulate chloride   [umol/m^3]
      PNa     =TNa                ! particulate sodium     [umol/m^3]
      PSO4    =TSO4               ! particulate sulfate    [umol/m^3]

      ! liquid matter
      ASO4 = ASO4 * 1.e6_dp  ! aqueous phase sulfate  [umol/m^3] 
      ANH4 = ANH4 * 1.e6_dp  ! aqueous phase ammonium [umol/m^3]
      ANO3 = ANO3 * 1.e6_dp  ! aqueous phase nitrate  [umol/m^3]
      ACl  = ACl  * 1.e6_dp  ! aqueous phase chloride [umol/m^3]
      ANa  = ANa  * 1.e6_dp  ! aqueous phase sodium   [umol/m^3]

      ! solid matter
      SNH4=PNH4(il)-ANH4  ! solid phase ammonium   [umol/m^3]
      SSO4=PSO4-ASO4      ! solid phase sulfate    [umol/m^3]
      SNO3=PNO3(il)-ANO3  ! solid phase nitrate    [umol/m^3]
      SCl =PCl -ACl       ! solid phase chloride   [umol/m^3]
      SNa =PNa -ANa       ! solid phase sodium     [umol/m^3]

      ! GO SAVE ...

      IF (SNH4 < 0._dp) SNH4=0._dp
      IF (SSO4 < 0._dp) SSO4=0._dp
      IF (SNO3 < 0._dp) SNO3=0._dp
      IF (SCl  < 0._dp) SCl =0._dp
      IF (SNa  < 0._dp) SNa =0._dp

      ! convert aerosol water from [kg/m^3] to [ug/m^3]
      WH2O(il)    = WH2O(il)    * 1.0e9_dp
      IF (WH2O(il) < 1.0e-3_dp) WH2O(il)    = 0.0_dp

      ! UPDATE HISTORY RH FOR HYSTERESIS (ONLINE CALCULATIONS ONLY)

      RH_HIST(il) = 2._dp                                       ! wet
!      IF(WH2O(il) == 0._dp) RH_HIST(il)=1._dp                   ! dry
      IF(WH2O(il) < ZEPS) RH_HIST(il)=1._dp                     ! dry

      !
      ! store aerosol species for diagnostic output:
      !___________________________________________________________
      ! Output values:

      ! convert from [umol/m3] to [ug/m3]
      GNH3(il)    = GNH3(il)    * MWNH3
      GNO3(il)    = GNO3(il)    * MWHNO3
!      GHCl(il)    = GHCl(il)    * MWHCl
      PNH4(il)    = PNH4(il)    * MWNH4
      PNO3(il)    = PNO3(il)    * MWNO3
   enddo

END SUBROUTINE MADE_EQSAM

!------------------------------------------------------------------------------!

SUBROUTINE MADE_MODPAR( BLKSIZE, NUMCELLS,          &
                        CBLK,                       &
                        BLKTA, BLKPRS,              &
                        PMASSN, PMASSA, PMASSC,     &
                        PDENSN, PDENSA, PDENSC,     &
                        XLM, AMU,                   &
                        DGNUC, DGACC, DGCOR,        &
                        KNNUC, KNACC, KNCOR )

   !***********************************************************************
   !
   !**    DESCRIPTION:
   !       Calculates modal parameters and derived variables,
   !       log-squared of std deviation, mode mean size, Knudsen number)
   !       based on current values of moments for the modes.
   ! FSB   Now calculates the 3rd moment, mass, and density in all 3 modes.
   !**
   !**    Revision history:
   !       Adapted 3/95 by US and CJC from EAM2's MODPAR and INIT3
   !       Revised  7/23/96 by FSB to use COMMON blocks and small blocks
   !       instead of large 3-d arrays, and to assume a fixed std.
   !       Revised 12/06/96 by FSB to include coarse mode
   !       Revised 1/10/97 by FSB to have arrays passed in call vector
   !**********************************************************************

   IMPLICIT NONE
   INTRINSIC MAX, SQRT

   ! *** input ***

   INTEGER, INTENT(in)     :: BLKSIZE               ! dimension of arrays
   INTEGER, INTENT(in)     :: NUMCELLS              ! actual number of cells in
                                                    ! arrays
   REAL(dp), INTENT(inout) :: CBLK(BLKSIZE,NSPCSDA) ! main array of variables
   REAL(dp), INTENT(in)    :: BLKTA(BLKSIZE)        ! Air temperature [K]
   REAL(dp), INTENT(in)    :: BLKPRS(BLKSIZE)       ! Air pressure in [Pa]

   ! *** output ***

   REAL(dp), INTENT(out)   :: PMASSN(BLKSIZE)       ! mass concentration in
                                                    ! Aitken mode [ug/m3]
   REAL(dp), INTENT(out)   :: PMASSA(BLKSIZE)       ! mass concentration in
                                                    ! acc. mode [ug/m3]
   REAL(dp), INTENT(out)   :: PMASSC(BLKSIZE)       ! mass concentration in
                                                    ! coarse mode [ug/m3]
   REAL(dp), INTENT(out)   :: PDENSN(BLKSIZE)       ! average particle density
                                                    ! in Aitken mode [kg/m3]
   REAL(dp), INTENT(out)   :: PDENSA(BLKSIZE)       ! average particle density
                                                    ! in acc. mode [kg/m3]
   REAL(dp), INTENT(out)   :: PDENSC(BLKSIZE)       ! average particle density
                                                    ! in coarse mode [kg/m3]
   REAL(dp), INTENT(out)   :: XLM(BLKSIZE)          ! atm. mean free path [m]
   REAL(dp), INTENT(out)   :: AMU(BLKSIZE)          ! atm. dynamic viscosity
                                                    ! [kg m-1 s-1]
   REAL(dp), INTENT(out)   :: DGNUC(BLKSIZE)        ! mean diameter, Aitken
                                                    ! mode [m]
   REAL(dp), INTENT(out)   :: DGACC(BLKSIZE)        ! d acc. mode
   REAL(dp), INTENT(out)   :: DGCOR(BLKSIZE)        ! d coarse mode
   REAL(dp), INTENT(out)   :: KNNUC(BLKSIZE)        ! Knudsen number, Aitken
                                                    ! mode
   REAL(dp), INTENT(out)   :: KNACC(BLKSIZE)        ! Kn acc. mode
   REAL(dp), INTENT(out)   :: KNCOR(BLKSIZE)        ! Kn coarse mode

   ! *** parameter ***

   REAL(dp), PARAMETER :: CONMIN  = 1.0E-30_dp  ! conc. lower limit [ug/m3]
   REAL(dp), PARAMETER :: DGMIN   = 1.0E-9_dp   ! lowest particle diameter [m]
   REAL(dp), PARAMETER :: DENSMIN = 1.0E03_dp   ! lowest particle dens. [kg/m3]

   ! *** local ***

   INTEGER :: LCELL                ! loop counter

   ! *** set up  aerosol  3rd moment, mass, density

   DO  LCELL = 1, NUMCELLS

      ! *** Aitken-mode

      CBLK(LCELL,VNU3) = (   SO4FAC  * CBLK(LCELL,VSO4AI)    &
                           + NH4FAC  * CBLK(LCELL,VNH4AI)    &
                           + H2OFAC  * CBLK(LCELL,VH2OAI)    &
                           + NO3FAC  * CBLK(LCELL,VNO3AI)    &
!sorgam                    + ORGFAC  * CBLK(LCELL,VORGARO1I) &
!sorgam                    + ORGFAC  * CBLK(LCELL,VORGARO2I) &
!sorgam                    + ORGFAC  * CBLK(LCELL,VORGALK1I) &
!sorgam                    + ORGFAC  * CBLK(LCELL,VORGOLE1I) &
!sorgam                    + ORGFAC  * CBLK(LCELL,VORGBA1I)  &
!sorgam                    + ORGFAC  * CBLK(LCELL,VORGBA2I)  &
                           + ORGFAC  * CBLK(LCELL,VORGPAI)   &
!unused                    + ANTHFAC * CBLK(LCELL,VP25AI)    &
                           + ANTHFAC * CBLK(LCELL,VECI)      &
                           + SEASFAC * CBLK(LCELL,VSEASI) )

      CBLK(LCELL,VNU3) = MAX(CONMIN, CBLK(LCELL, VNU3))

      ! *** accumulation mode

      CBLK(LCELL,VAC3) = (   SO4FAC *  CBLK(LCELL,VSO4AJ)    &
                           + NH4FAC *  CBLK(LCELL,VNH4AJ)    &
                           + H2OFAC *  CBLK(LCELL,VH2OAJ)    &
                           + NO3FAC *  CBLK(LCELL,VNO3AJ)    &
!sorgam                    + ORGFAC *  CBLK(LCELL,VORGARO1J) &
!sorgam                    + ORGFAC *  CBLK(LCELL,VORGARO2J) &
!sorgam                    + ORGFAC *  CBLK(LCELL,VORGALK1J) &
!sorgam                    + ORGFAC *  CBLK(LCELL,VORGOLE1J) &
!sorgam                    + ORGFAC *  CBLK(LCELL,VORGBA1J)  &
!sorgam                    + ORGFAC *  CBLK(LCELL,VORGBA2J)  &
                           + ORGFAC *  CBLK(LCELL,VORGPAJ)   &
!unused                    + ANTHFAC * CBLK(LCELL,VP25AJ)    &
                           + ANTHFAC * CBLK(LCELL,VECJ)      &
                           + SEASFAC * CBLK(LCELL,VSEASJ)    &
                           + SOILFAC * CBLK(LCELL,VDUSTJ) )

      CBLK(LCELL,VAC3) = MAX(CONMIN, CBLK(LCELL,VAC3))

      ! *** coarse mode

      CBLK(LCELL,VCOR3) = (   SOILFAC * CBLK(LCELL,VDUSTC)   &
                            + SEASFAC * CBLK(LCELL,VSEASC)   &
                            + H2OFAC  * CBLK(LCELL,VH2OAC) )

      CBLK(LCELL,VCOR3) = MAX(CONMIN, CBLK(LCELL,VCOR3))

      ! *** now get particle mass and density

      ! *** Aitken-mode:

      PMASSN(LCELL) = (   CBLK(LCELL,VSO4AI)    &
                        + CBLK(LCELL,VNH4AI)    &
                        + CBLK(LCELL,VH2OAI)    &
                        + CBLK(LCELL,VNO3AI)    &
!sorgam                 + CBLK(LCELL,VORGARO1I) &
!sorgam                 + CBLK(LCELL,VORGARO2I) &
!sorgam                 + CBLK(LCELL,VORGALK1I) &
!sorgam                 + CBLK(LCELL,VORGOLE1I) &
!sorgam                 + CBLK(LCELL,VORGBA1I)  &
!sorgam                 + CBLK(LCELL,VORGBA2I)  &
                        + CBLK(LCELL,VORGPAI)   &
!unused                 + CBLK(LCELL,VP25AI)    &
                        + CBLK(LCELL,VSEASI)    &
                        + CBLK(LCELL,VECI) )

      PMASSN(LCELL) = MAX(CONMIN, PMASSN(LCELL))

      ! *** accumulation mode:

      PMASSA(LCELL) = (   CBLK(LCELL,VSO4AJ)    &
                        + CBLK(LCELL,VNH4AJ)    &
                        + CBLK(LCELL,VH2OAJ)    &
                        + CBLK(LCELL,VNO3AJ)    &
!sorgam                 + CBLK(LCELL,VORGARO1J) &
!sorgam                 + CBLK(LCELL,VORGARO2J) &
!sorgam                 + CBLK(LCELL,VORGALK1J) &
!sorgam                 + CBLK(LCELL,VORGOLE1J) &
!sorgam                 + CBLK(LCELL,VORGBA1J)  &
!sorgam                 + CBLK(LCELL,VORGBA2J)  &
                        + CBLK(LCELL,VORGPAJ)   &
!unused                 + CBLK(LCELL,VP25AJ)    &
                        + CBLK(LCELL,VECJ)      &
                        + CBLK(LCELL,VSEASJ)    &
                        + CBLK(LCELL,VDUSTJ) )

      PMASSA(LCELL) = MAX(CONMIN, PMASSA(LCELL))

      ! *** coarse mode:

      PMASSC(LCELL) = (   CBLK(LCELL,VSEASC)    &
                        + CBLK(LCELL,VDUSTC)    &
!unused                 + CBLK(LCELL,VANTHA)    &
                        + CBLK(LCELL,VH2OAC) )

      PMASSC(LCELL) = MAX(CONMIN, PMASSC(LCELL))

      ! *** now get particle density, mean free path, and dynamic viscosity

      ! Density and mean free path
      ! *** density in [kg m-3]

      PDENSN(LCELL) = MAX(DENSMIN, (F6DPIM9*PMASSN(LCELL) / CBLK(LCELL,VNU3)))
      PDENSA(LCELL) = MAX(DENSMIN, (F6DPIM9*PMASSA(LCELL) / CBLK(LCELL,VAC3)))
      PDENSC(LCELL) = MAX(DENSMIN, (F6DPIM9*PMASSC(LCELL) / CBLK(LCELL,VCOR3)))

      ! *** Calculate mean free path [m]:

      XLM(LCELL) = 6.6328E-8_dp*P0*BLKTA(LCELL) / (TS0*BLKPRS(LCELL))

      ! *** 6.6328E-8 is the sea level values given in Table I.2.8
      ! *** on page 10 of U.S. Standard Atmosphere 1962

      ! ***     Calcualte dynamic viscosity [ kg m-1 s-1]:

      ! *** U.S. Standard Atmosphere 1962 page 14 expression
      !     for dynamic viscosity is:
      !     dynamic viscosity =  beta * T * sqrt(T) / (T + S)
      !     where beta = 1.458e-6 [kg sec^-1 K**-0.5], S = 110.4 [K].

      AMU(LCELL) = 1.458E-6_dp * BLKTA(LCELL) * SQRT(BLKTA(LCELL)) / &
                   (BLKTA(LCELL) + 110.4_dp)

      ! Standard deviation fixed in both modes, so diagnose diameter
      ! from 3rd moment and number concentrations:

      ! calculate diameters

      DGNUC(LCELL) = MAX(DGMIN, (CBLK(LCELL,VNU3)  / (CBLK(LCELL,VNU0)  &
                                 * ES36(akn)))**ONE3)
      DGACC(LCELL) = MAX(DGMIN, (CBLK(LCELL,VAC3)  / (CBLK(LCELL,VAC0)  &
                                 * ES36(acc)))**ONE3)
      DGCOR(LCELL) = MAX(DGMIN, (CBLK(LCELL,VCOR3) / (CBLK(LCELL,VCORN) &
                                 * ES36(cor)))**ONE3)

      ! calculate Knudsen numbers

      KNNUC(LCELL) = 2.0_dp * XLM(LCELL) / DGNUC(LCELL)
      KNACC(LCELL) = 2.0_dp * XLM(LCELL) / DGACC(LCELL)
      KNCOR(LCELL) = 2.0_dp * XLM(LCELL) / DGCOR(LCELL)

   END DO ! end loop over cells

END SUBROUTINE MADE_MODPAR

!------------------------------------------------------------------------------!

SUBROUTINE MADE_COAGRATE( BLKSIZE, NUMCELLS,          &
                          CBLK,                       &
                          BLKTA,                      &
                          PDENSN, PDENSA,             &
                          AMU,                        &
                          DGNUC, DGACC,               &
                          KNNUC, KNACC,               &
                          URN00, URA00,               &
                          BRNA01, C30 )

   !***********************************************************************
   !**    DESCRIPTION:  calculates aerosol coagulation rates for unimodal
   !       and bimodal coagulation using E. Whitby 1990's prescription.
   !
   !.......   Rates for coaglulation:
   !.......   Unimodal Rates:
   !.......   URN00:  nuclei       mode 0th moment self-coagulation rate
   !.......   URA00:  accumulation mode 0th moment self-coagulation rate
   !
   !.......   Bimodal Rates:  (only 1st order coeffs appear)
   !.......   NA-- nuclei  with accumulation coagulation rates,
   !.......   AN-- accumulation with nuclei coagulation rates
   !.......   BRNA01:  rate for 0th moment ( d(nuclei mode 0) / dt term)
   !.......   BRNA31:        "  3rd        ( d(nuclei mode 3) / dt term)
   !**
   !**
   !**    Revision history:
   !       prototype 1/95 by Uma and Carlie
   !       Revised   8/95 by US for calculation of density from stmt func
   !                 and collect met variable stmt funcs in one include
   !                 file
   !      REVISED 7/25/96 by FSB to use block structure
   !      REVISED 9/13/96 BY FSB for Uma's FIXEDBOTH case only.
   !      REVISED 11/08/96 BY FSB the Whitby Shankar convention on signs
   !                              changed. All coagulation coefficients
   !                              returned with positive signs. Their
   !                              linearization is also abandoned.
   !                              Fixed values are used for the corrections
   !                              to the free-molecular coagulation
   !                              integrals.
   !                              The code forces the harmonic means to be
   !                              evaluated in 64 bit arithmetic on 32 bit
   !                              machines.
   !      REVISED 11/14/96 BY FSB Internal units are now MKS, moment /
   !                              unit-volume
   !
   !      REVISED 1/12/98 by FSB  C30 replaces BRNA31 as an array. This
   !                              was done
   !                              because BRNA31 can become zero on a
   !                              workstation because of limited precision.
   !                              With the change in aerostep to omit
   !                              update of the 3rd moment C30 is the only
   !                              variable now needed.
   !                              The logic using ONE88 to force REAL*8
   !                              arithmetic has been removed and all
   !                              intermediates are now REAL*8.
   !
   !**********************************************************************

   IMPLICIT NONE
   INTRINSIC SQRT

   ! *** input ***

   INTEGER,  INTENT(in)  :: BLKSIZE               ! dimension of arrays
   INTEGER,  INTENT(in)  :: NUMCELLS              ! actual number of cells in
                                                  ! arrays
   REAL(dp), INTENT(in)  :: CBLK(BLKSIZE,NSPCSDA) ! main array of variables
   REAL(dp), INTENT(in)  :: BLKTA(BLKSIZE)        ! Air temperature [K]
   REAL(dp), INTENT(in)  :: PDENSN(BLKSIZE)       ! average particle density in
                                                  ! Aitken mode [kg/m3]
   REAL(dp), INTENT(in)  :: PDENSA(BLKSIZE)       ! average particle density in
                                                  ! accumulation mode [kg/m3]
   REAL(dp), INTENT(in)  :: AMU(BLKSIZE)          ! atm. dynamic viscosity
                                                  ! [kg m^-1 s^-1]
   REAL(dp), INTENT(in)  :: DGNUC(BLKSIZE)        ! Aitken mode mean
                                                  ! diameter [m]
   REAL(dp), INTENT(in)  :: DGACC(BLKSIZE)        ! acc. mode mean diameter [m]
   REAL(dp), INTENT(in)  :: KNNUC(BLKSIZE)        ! Aitken mode Knudsen number
   REAL(dp), INTENT(in)  :: KNACC(BLKSIZE)        ! acc. mode Knudsen number

   ! *** output ***

   REAL(dp), INTENT(out) :: URN00(BLKSIZE)        ! intramodal coagulation rate
                                                  ! (Aitken mode)
   REAL(dp), INTENT(out) :: URA00(BLKSIZE)        ! intramodal coagulation rate
                                                  ! (acc. mode)
   REAL(dp), INTENT(out) :: BRNA01(BLKSIZE)       ! intermodal coagulaton rate
                                                  ! (Aitken <--> acc., number)
   REAL(dp), INTENT(out) :: C30(BLKSIZE)          ! intermodal 3rd moment
                                                  ! transfer rate by
                                                  ! intermodal coagulation
                                                  ! (Aitken <--> acc.)

   ! *** local variables (must be 64 bit ---> REAL*8) ***

   REAL(dp) :: KNCNUC, KNCACC         ! coeffs for unimodal NC coag rate
   REAL(dp) :: KFMNUC, KFMACC         ! coeffs for unimodal FM coag rate
   REAL(dp) :: KNC, KFM               ! coeffs for bimodal NC, FM coag rate
   REAL(dp) :: BENCNN, BENCNA         ! NC 0th moment coag rate (both modes)
   REAL(dp) :: BENCM3N                ! NC 3rd moment coag rate (Aitken mode)
   REAL(dp) :: BEFMNN, BEFMNA         ! FM 0th moment coag rate (both modes)
   REAL(dp) :: BEFM3N                 ! FM 3rd moment coag rate (Aitken mode)
   REAL(dp) :: BETANN, BETANA         ! composite coag rates, mom 0 (both modes)
   REAL(dp) :: BRNA31                 ! intermodal coagulation rate for
                                      ! 3rd moments
   REAL(dp) :: S1                     ! scratch subexpression
   REAL(dp) :: T1, T2                 ! scratch subexpressions
   REAL(dp) :: T16 !, T26             ! T1**6, T2**6
   REAL(dp) :: RAT, RIN               ! ratio of acc to nuc size and its invers
   REAL(dp) :: RSQT, RSQ4             ! sqrt( rat ), rsqt**4
   REAL(dp) :: RSQTI, RSQI3           ! sqrt( 1/rat ), sqrt( 1/rat**3 )
   REAL(dp) :: DGN3                   ! dgnuc**3
   REAL(dp) :: DGA3                   ! dgacc**3

   ! *** Fixed values for correctionss to coagulation integrals
   !     for free-molecular (FM) case (must be 64 bit ---> REAL*8). ***

   REAL(dp), PARAMETER :: BM0  = 0.8_dp   ! 0.8D0
   REAL(dp), PARAMETER :: BM0I = 0.9_dp   ! 0.9D0
   REAL(dp), PARAMETER :: BM3I = 0.9_dp   ! 0.9D0
   REAL(dp), PARAMETER :: A    = 1.246_dp ! approx Cunningham corr. factor
                                          ! 1.246D0

   INTEGER :: LCELL                         ! loop counter

   !.......................................................................
   !   begin body of subroutine  COAGRATE

   !...........   Main computational grid-traversal loops
   !...........   for computing coagulation rates.

   ! *** Both modes have fixed std devs. ***

   DO LCELL = 1, NUMCELLS     !  loop on LCELL

      ! *** moment independent factors

      S1 = TWO3 * BOLTZ *  BLKTA(LCELL) / AMU(LCELL)

      ! For unimodal coagualtion:

      KNCNUC = S1
      KNCACC = S1

      KFMNUC = SQRT(3.0_dp * BOLTZ * BLKTA(LCELL) / PDENSN(LCELL))
      KFMACC = SQRT(3.0_dp * BOLTZ * BLKTA(LCELL) / PDENSA(LCELL))

      ! For bimodal coagulation:

      KNC  = S1
      KFM  = SQRT(6.0_dp * BOLTZ * BLKTA(LCELL) &
                  / (PDENSN(LCELL) + PDENSA(LCELL)))

      ! Begin unimodal coagulation rate calculations:

      ! Near-continuum (NC) regime.

      DGN3 = DGNUC(LCELL)**3
      DGA3 = DGACC(LCELL)**3

      T1  = SQRT(DGNUC(LCELL))
      T2  = SQRT(DGACC(LCELL))
      T16 = DGN3          ! = T1**6
!      T26 = DGA3         ! = T2**6

      ! Note rationalization of fractions and subsequent cancellations
      ! from the formulation in  Whitby et al. (1990)

      BENCNN = KNCNUC * (1.0_dp + ES08(akn) &
                         + A * KNNUC(LCELL) * (ES04(akn) + ES20(akn)))
      BENCNA = KNCACC * (1.0_dp + ES08(acc) &
                         + A * KNACC(LCELL) * (ES04(acc) + ES20(acc)))

      ! Free molecular (FM) regime. Uses fixed value for correction factor BM0

      BEFMNN = KFMNUC * T1 * (E1(akn) + ES25(akn) + 2.0_dp * ES05(akn)) * BM0
      BEFMNA = KFMACC * T2 * (E1(acc) + ES25(acc) + 2.0_dp * ES05(acc)) * BM0

      ! Calculate half the harmonic mean between unimodal rates
      ! free molecular (FM) and near-continuum (NC) regimes.

      ! FSB    64 bit evaluation

      BETANN = BENCNN * BEFMNN / (BENCNN + BEFMNN)
      BETANA = BENCNA * BEFMNA / (BENCNA + BEFMNA)

      URN00(LCELL) = BETANN
      URA00(LCELL) = BETANA

      ! *** End of unimodal coagulation calculations.

      ! Begin bimodal coagulation rate calculations:

      RAT  = DGACC(LCELL) / DGNUC(LCELL)
      RIN  = 1.0_dp / RAT   ! 1.0D0 / RAT
      RSQT = SQRT(RAT)
      RSQ4 = RAT**2

      RSQTI = 1.0_dp / RSQT  ! 1.0D0_dp / RSQT
      RSQI3 = RIN * RSQTI

      ! Near-continuum coeffs:
      ! 0th moment Aitken mode bimodal coag coefficient

      BENCNN = KNC * (2.0_dp + A * KNNUC(LCELL)                             &
                               * (ES04(akn) + RAT * ES16(akn) * ES04(acc))  &
                             + A * KNACC(LCELL)                             &
                               * (ES04(acc) + RIN * ES16(acc) * ES04(akn))  &
                             + (RAT + RIN) * ES04(akn) * ES04(acc))

      ! 3rd moment Aitken mode bimodal coag coefficient

      BENCM3N = KNC * DGN3 * (2.0_dp * ES36(akn)                            &
                + A * KNNUC(LCELL)                                          &
                    * (ES16(akn) + RAT * ES04(akn) * ES04(acc))             &
                + A * KNACC(LCELL)                                          &
                    * (ES36(akn) * ES04(acc) + RIN * ES64(akn) * ES16(acc)) &
                + RAT * ES16(akn) * ES04(acc) + RIN * ES64(akn) * ES04(acc))

      ! Free molecular regime coefficients:
      ! Uses fixed value for correction factor BM0I, BM3I

      ! 0th moment nuc mode coeff

      BEFMNN = KFM * BM0I * T1 * (E1(akn) + RSQT * E1(acc)                  &
                                + 2.0_dp * RAT   * E1(akn)     * ES04(acc)  &
                                +          RSQ4  * ES09(akn)   * ES16(acc)  &
                                +          RSQI3 * ES16(akn)   * ES09(acc)  &
                                + 2.0_dp * RSQTI * ES04(akn)   * E1(acc))

      ! 3rd moment nuc mode coeff

      BEFM3N = KFM * BM3I * T1 * T16 * (ES49(akn)                           &
                                +          RSQT  * ES36(akn)  * E1(acc)     &
                                + 2.0_dp * RAT   * ES25(akn)  * ES04(acc)   &
                                +          RSQ4  * ES09(akn)  * ES16(acc)   &
                                +          RSQI3 * ES100(akn) * ES09(acc)   &
                                + 2.0_dp * RSQTI * ES64(akn)  * E1(acc))

      ! Calculate half the harmonic mean between bimodal rates
      ! free molecular and near-continuum regimes.

      ! FSB    Force 64 bit evaluation

      BRNA01(LCELL) = BENCNN  * BEFMNN / (BENCNN  + BEFMNN)
      ! BRNA31 now is a scalar
      BRNA31        = BENCM3N * BEFM3N / (BENCM3N + BEFM3N)
      ! 3rd moment transfer by intermodal coagulation
      C30(LCELL)    = BRNA31  * CBLK(LCELL,VAC0) * CBLK(LCELL,VNU0)

      ! End bimodal coagulation rate.

   END DO ! end of main lop over cells

END SUBROUTINE MADE_COAGRATE

!------------------------------------------------------------------------------!

SUBROUTINE MADE_NUCLCOND( BLKSIZE, NUMCELLS,            &
                          CBLK, DT,                     &
                          BLKTA, BLKPRS, BLKRH,         &
                          SO4RAT, SOA_MADE,             &
!sorgam                   ORGARO1RAT, ORGARO2RAT,       &
!sorgam                   ORGALK1RAT, ORGOLE1RAT,       &
!sorgam                   ORGBIO1RAT, ORGBIO2RAT,       &
!sorgam                   ORGBIO3RAT, ORGBIO4RAT,       &
!sorgam                   DROG, LDROG, NCV, NACV,       &
                          DGNUC, DGACC,                 &
                          FCONCN, FCONCA,               &
                          FCONCN_ORG, FCONCA_ORG,       &
                          DMDT, DNDT, DELTASO4A,        &
                          CGRN3, CGRA3 )

   !***********************************************************************
   !**    DESCRIPTION:  calculates aerosol nucleation and condensational
   !**    growth rates using Binkowski and Shankar (1995) method.
   !
   ! *** In this version, the method od RPM is followed where
   !     the diffusivity, the average molecular ve3locity, and
   !     the accomodation coefficient for sulfuric acid are used for
   !     the organics. This is for consistency.
   !     Future versions will use the correct values.  FSB 12/12/96
   !
   !
   !**
   !**    Revision history:
   !       prototype 1/95 by Uma and Carlie
   !       Corrected 7/95 by Uma for condensation of mass not nucleated
   !       and mass conservation check
   !       Revised   8/95 by US to calculate air density in stmt function
   !                 and collect met variable stmt funcs in one include
   !                 file
   !       Revised 7/25/96 by FSB to use block structure.
   !       Revised 9/17/96 by FSB to use Y&K or K&W Nucleation mechanism
   !       Revised 11/15/96 by FSB to use MKS,  and mom m^-3 units.
   !       Revised 1/13/97 by FSB to pass arrays and simplify code.
   !       Added   23/03/99 by BS growth factors for organics
   !**********************************************************************

   IMPLICIT NONE
   INTRINSIC EXP, MAX, MIN, SQRT

   ! *** arguments ***

   ! *** input ***

   INTEGER,  INTENT(in)    :: BLKSIZE               ! dimension of arrays
   INTEGER,  INTENT(in)    :: NUMCELLS              ! actual number of cells
                                                    ! in arrays
!sorgam   INTEGER,  INTENT(in)    :: LDROG          ! # of organic aerosol
!sorgam                                             ! precursor

   REAL(dp), INTENT(in)    :: CBLK(BLKSIZE,NSPCSDA) ! main array of variables
   REAL(dp), INTENT(in)    :: DT                    ! model time step [s]
   REAL(dp), INTENT(in)    :: BLKTA(BLKSIZE)        ! air temperature [K]
   REAL(dp), INTENT(in)    :: BLKPRS (BLKSIZE)      ! air pressure in [Pa]
   REAL(dp), INTENT(in)    :: BLKRH(BLKSIZE)        ! frac. rel. humidity (0-1)
   REAL(dp), INTENT(inout) :: SO4RAT(BLKSIZE)       ! sulfate gas-phase
                                                    ! production rate [ug/m3/s]
   REAL(dp), INTENT(in)    :: SOA_MADE(BLKSIZE)     ! SOA (gas phase)
                                                    ! "emissions" [ug/m3/s]
!sorgam
!sorgam   INTEGER :: NCV             ! total # of cond. vapors & SOA species
!sorgam   INTEGER :: NACV            ! # of anthrop. cond. vapors & SOA species
!sorgam
!sorgam   !bs * anthropogenic organic condensable vapor production rate
!sorgam   REAL(dp) :: DROG(BLKSIZE,LDROG)  ! Delta ROG conc. [ppm]
!sorgam
!sorgam   REAL(dp) :: ORGARO1RAT(BLKSIZE)  ! anthropogenic organic aerosol mass
!sorgam                                    ! production rate from aromatics
!sorgam                                    ! [ug/m3/s]
!sorgam   REAL(dp) :: ORGARO2RAT(BLKSIZE)  ! anthropogenic organic aerosol mass
!sorgam                                    ! production rate from aromatics
!sorgam                                    ! [ug/m3/s]
!sorgam   REAL(dp) :: ORGALK1RAT(BLKSIZE)  ! anthropogenic organic aerosol mass
!sorgam                                    ! production rate from alkanes &
!sorgam                                    ! others [ug/m3/s]
!sorgam   REAL(dp) :: ORGOLE1RAT(BLKSIZE)  ! anthropogenic organic aerosol mass
!sorgam                                    ! production rate from alkenes &
!sorgam                                    ! others [ug/m3/s]
!sorgam   !bs * biogenic organic condensable vapor production rate
!sorgam   REAL(dp) :: ORGBIO1RAT(BLKSIZE)  ! biogenic organic aerosol
!sorgam                                    ! production rate [ug/m3/s]
!sorgam   REAL(dp) :: ORGBIO2RAT(BLKSIZE)  ! biogenic organic aerosol
!sorgam                                    ! production rate [ug/m3/s]
!sorgam   REAL(dp) :: ORGBIO3RAT(BLKSIZE)  ! biogenic organic aerosol
!sorgam                                    ! production rate [ug/m3/s]
!sorgam   REAL(dp) :: ORGBIO4RAT(BLKSIZE)  ! biogenic organic aerosol
!sorgam                                    ! production rate [ug/m3/s]

   REAL(dp), INTENT(in)    :: DGNUC(BLKSIZE)        ! Aitken mode gemetric
                                                    ! mean diameter [m]
   REAL(dp), INTENT(in)    :: DGACC(BLKSIZE)        ! acc. mode geometric mean
                                                    ! diameter [m]

   ! *** output ***

   REAL(dp), INTENT(out)   :: FCONCN(BLKSIZE)       ! reciprocal condensation
                                                    ! rate, Aitken mode [1/s]
   REAL(dp), INTENT(out)   :: FCONCA(BLKSIZE)       ! reciprocal condensation
                                                    ! rate, acc. mode [1/s]
   REAL(dp), INTENT(out)   :: FCONCN_ORG(BLKSIZE)   ! reciprocal condensation
                                                    ! rate, Aitken mode [1/s]
   REAL(dp), INTENT(out)   :: FCONCA_ORG(BLKSIZE)   ! reciprocal condensation
                                                    ! rate, acc. mode [1/s]
   REAL(dp), INTENT(out)   :: DMDT(BLKSIZE)         ! rate of production of
                                                    ! new so4 mass by particle
                                                    ! formation [ug/m3/s]
   REAL(dp), INTENT(out)   :: DNDT(BLKSIZE)         ! rate of producton of new
                                                    ! particle number by
                                                    ! particle formation
                                                    ! [#/m3/s]
   REAL(dp), INTENT(out)   :: DELTASO4A(BLKSIZE)    ! increment of conc. added
                                                    ! to sulfate aerosol by
                                                    ! condensation [ug/m3]
   REAL(dp), INTENT(out)   :: CGRN3(BLKSIZE)        ! growth rate for 3rd
                                                    ! moment for Aitken mode
                                                    ! [3rd mom/m3/s]
   REAL(dp), INTENT(out)   :: CGRA3(BLKSIZE)        ! growth rate for 3rd
                                                    ! moment for acc. mode
                                                    ! [3rd mom/m3/s]

   ! SCRATCH local variables and their descriptions:

   INTEGER :: LCELL                ! LOOP INDEX

   REAL(dp) :: CONDRATE            ! condensation rate [mom-3/g/s]
   REAL(dp) :: CHEMRAT_ORG         ! conv rate for organics [mom-3/g/s]
   REAL(dp) :: AM1N, AM1A          ! 1st mom density (Aitken, acc. modes)
                                   ! [mom_1/g-air]
   REAL(dp) :: AM2N, AM2A          ! 2nd mom density (Aitken, acc. modes)
                                   ! [mom_2/g-air]
   REAL(dp) :: GNC3N, GNC3A        ! near-cont fns (Aitken, acc) for mom-3
                                   ! density
   REAL(dp) :: GFM3N, GFM3A        ! free-mol  fns (Aitken, acc) for mom-3
                                   ! density
   REAL(dp) :: FCONC               ! total reciprocal condensation rate
   REAL(dp) :: TD                  ! d * tinf (cgs)

   ! Constant to force 64 bit evaluation of an expression
   ! (must be 64 bit ---> REAL*8)
   REAL(dp), PARAMETER :: ONE88 = 1.0_dp ! 1.0D0

   ! *** variables to set up sulfate and organic condensation rates

!   REAL(dp) :: VAPOR1(BLKSIZE)     ! sulfuric acid vapor at current time step
   REAL(dp) :: VAPOR2(BLKSIZE)     ! Sulfuric acid vapor prior to addition
                                   ! from chemistry and emissions
   REAL(dp) :: OLDSULF(BLKSIZE)    ! "old" conc. of sulfuric acid vapor [ug/m3]
   REAL(dp) :: DELTAVAP            ! change to vapor at previous time step
                                   ! including condensation and chem. production

   !bs * start update

   REAL(dp) :: DIFFCORR

   !bs     REAL ALPHSULF ! Accommodation coefficient for sulfuric acid
   !bs     PARAMETER ( ALPHSULF = 0.05 ) ! my be set to one in future
   !bs
   !bs     REAL DIFFSULF ! molecular diffusivity for sulfuric acid
   !bs                   ! [ m**2 /sec ]
   !bs     PARAMETER( DIFFSULF = 0.08E-4 ) ! may be changed in future
   !bs
   !bs * 23/03/99 updates of ALPHSULF and DIFFSULF adopted fro new code
   !bs * from FSB
   !bs * DIFFSULF is calculated from Reid, Prausnitz, and Poling, The
   !bs * properties
   !bs * of gases and liquids, 4th edition, McGraw-Hill, 1987, pp 587-588.
   !bs * Equation (11-4.4) was used.
   !bs * The value is at T = 273.16 K and P = 1.01325E05 Pa
   !bs * Temperature dependence is included for DIFFSULF via DIFFCORR
   !bs * (see below).

   REAL(dp), PARAMETER :: ALPHSULF = 1.0_dp ! accommodation coefficient for
                                            ! sulfuric acid
                                            !bs updated from code of FSB

   REAL(dp), PARAMETER :: DIFFSULF = 9.362223E-06_dp ! molecular diffusivity for
                                                     ! sulfuric acid [m2/s]
                                                     !bs updated from code of
                                                     !bs FSB

   REAL(dp), PARAMETER :: ALPHAORG = 1.0_dp          !bs Accomodation coef.
                                                     !bs for organics
                                                     !bs Bowman et al. '97 use
                                                     !bs alpha = 1.
                                                     !bs Kleeman et al. '99
                                                     !bs propose alpha = 0.1

   !bs * DIFFSULF is calculated from Reid, Prausnitz, and Poling, The
   !bs * properties of gases and liquids, 4th edition, McGraw-Hill, 1987,
   !bs * pp 587-588.
   !bs * Equation (11-4.4) was used.
   !bs * The value is at T = 273.16 K and P = 1.01325E05 Pa

   REAL(dp), PARAMETER :: DIFFORG = 5.151174E-06_dp !bs molecular diff. for
                                                    !bs organics [m2/s]

   REAL(dp), PARAMETER :: CCONC_ORG = 2.0_dp * PI * DIFFORG !bs factor for NC
                                                            !bs condensation for
                                                            !bs organics

   !bs analogue to CCOFM but for organics
   REAL(dp) :: CCOFM_ORG
   REAL(dp) :: CSQT_ORG 

   !bs * end update

   ! FSB  CCOFM is  the accommodation coefficient
   !      times the mean molecular velocity for h2so4 without the
   !      temperature after some algebra

   REAL(dp) :: CCOFM
   REAL(dp) :: CSQT        ! ccofm * sqrt(ta)

   ! *** CCONC is the factor for near-continuum condensation.
   !      times the mean molecular velocity for h2so4 without the
   !      temperature after some algebra

   REAL(dp), PARAMETER :: CCONC = 2.0_dp * PI * DIFFSULF ! [m2/s]

   LOGICAL :: FIRSTIME
   DATA FIRSTIME / .TRUE. /
   SAVE FIRSTIME

   SAVE CCOFM
   SAVE CCOFM_ORG

   !.......................................................................
   ! begin body of subroutine MADE_NUCLCOND

   IF (FIRSTIME) THEN
      ! *** set factor for free molecular condensation
      CCOFM     = ALPHSULF * SQRT(PI*RGASUNIV/(2.0_dp*MWH2SO4*1.0E-3_dp))
      CCOFM_ORG = ALPHAORG * SQRT(PI*RGASUNIV/(2.0_dp*MWORG*1.0E-3_dp))
      FIRSTIME  = .FALSE.
   END IF

   ! Array initializations
   DMDT = 0.0_dp
   DNDT = 0.0_dp

   ! Main computational grid-traversal loop nest
   ! for computing condensation and nucleation:

   DO LCELL = 1, NUMCELLS ! 1st loop over NUMCELLS

      ! *** First moment:
      AM1N = CBLK(LCELL,VNU0) * DGNUC(LCELL) * ES04(akn)
      AM1A = CBLK(LCELL,VAC0) * DGACC(LCELL) * ES04(acc)

      ! near-continuum factors [1/s]

      !bs * adopted from code of FSB
      !bs * correction to DIFFSULF for temperature and pressure

      DIFFCORR = (P0 / BLKPRS(LCELL)) * (BLKTA(LCELL) / 273.16_dp)**1.75_dp

      GNC3N = CCONC * AM1N * DIFFCORR
      GNC3A = CCONC * AM1A * DIFFCORR

      ! *** Second moment:

      AM2N = CBLK(LCELL,VNU0) * DGNUC(LCELL) * DGNUC(LCELL) * ES16(akn)
      AM2A = CBLK(LCELL,VAC0) * DGACC(LCELL) * DGACC(LCELL) * ES16(acc)

      CSQT  = CCOFM * SQRT(BLKTA(LCELL)) ! put in temperature factor

      ! free molecular factors [1/s]

      GFM3N = CSQT * AM2N
      GFM3A = CSQT * AM2A

      ! *** Condensation factors in [s-1] for h2so4
      ! *** In the future, separate factors for condensing organics will
      !     be included. In this version, the h2so4 values are used.

      ! Twice the harmonic mean of fm, nc functions:

      ! *** Force 64 bit evaluation:

      FCONCN(LCELL) = ONE88 * GNC3N * GFM3N / (GNC3N + GFM3N)
      FCONCA(LCELL) = ONE88 * GNC3A * GFM3A / (GNC3A + GFM3A)
      FCONC = FCONCN(LCELL) + FCONCA(LCELL)

      ! *** NOTE: FCONCN and FCONCA will be redefined below <<<<<<

      !bs * start modifications for organcis

      GNC3N = CCONC_ORG * AM1N * DIFFCORR
      GNC3A = CCONC_ORG * AM1A * DIFFCORR

      CSQT_ORG  = CCOFM_ORG * SQRT(BLKTA(LCELL))
      GFM3N = CSQT_ORG * AM2N
      GFM3A = CSQT_ORG * AM2A

      FCONCN_ORG(LCELL) = ONE88 * GNC3N * GFM3N / (GNC3N + GFM3N)
      FCONCA_ORG(LCELL) = ONE88 * GNC3A * GFM3A / (GNC3A + GFM3A)

      !bs * end modifications for organics

      ! *** calculate the total change to sulfuric acid vapor from production
      !     and condensation

      ! VAPOR1: curent sulfuric acid vapor
!      VAPOR1(LCELL) = CBLK(LCELL,VSULF)
      ! OLDSULF: vapor at prev. time step
      OLDSULF(LCELL) = MAX(0.0_dp, CBLK(LCELL,VSULF) - SO4RAT(LCELL)*DT)

      DELTAVAP = (SO4RAT(LCELL) / FCONC - OLDSULF(LCELL)) * &
                 (1.0_dp - EXP(-FCONC * DT))

      ! VAPOR2: H2SO4(g) conc. taking into account changes due to
      !         condensation and chemical production
      VAPOR2(LCELL) = MAX(0.0_dp, OLDSULF(LCELL) + DELTAVAP)

      ! *** Calculate increment in total sufate aerosol mass concentration
      !     (by condensation)

      ! *** This follows the method of Youngblood & Kreidenweis.

      DELTASO4A(LCELL) = MIN(SO4RAT(LCELL) * DT - DELTAVAP, CBLK(LCELL,VSULF))
      DELTASO4A(LCELL) = MAX(0.0_dp, DELTASO4A(LCELL) * MWSO4/MWH2SO4)

      ! *** zero out growth coefficients

      CGRN3(LCELL) = 0.0_dp
      CGRA3(LCELL) = 0.0_dp

      !al   Redefine SO4RAT to be maximum mass production rate allowed
      !al   by nucleation:
      !al   = (available H2SO4 vapor - amount condensed) / DT
      !al   This ensures no more H2SO4 vapor to be nucleated than available.

      SO4RAT(LCELL) = (CBLK(LCELL,VSULF) * MWSO4 / MWH2SO4 - &
                       DELTASO4A(LCELL)) / DT

   ENDDO                     ! End 1st loop over NUMCELLS

   IF (l_nuc) THEN
      ! *** Do Vehkam?ki et al., 2002 (J. Geophys. Res.) nucleation
      CALL MADE_NUKLEFIT( BLKTA, BLKRH, VAPOR2, SO4RAT, BLKSIZE, NUMCELLS, &
           DNDT, DMDT )
   END IF

!sorgam   !bs * Secondary organic aerosol module (SORGAM)
!sorgam 
!sorgam   CALL SORGAM( BLKTA, BLKPRS,
!sorgam                ORGARO1RAT, ORGARO2RAT,
!sorgam                ORGALK1RAT, ORGOLE1RAT,
!sorgam                ORGBIO1RAT, ORGBIO2RAT,
!sorgam                ORGBIO3RAT, ORGBIO4RAT,
!sorgam                DROG, LDROG, NCV, NACV,
!sorgam                CBLK, BLKSIZE, NUMCELLS,
!sorgam                DT )
!sorgam
!sorgam   ! *  Secondary organic aerosol module (SORGAM)

   DO LCELL = 1, NUMCELLS    !2nd loop over numcells

      ! *** redefine FCONCN & FCONCA to be the nondimensional fractional
      !     condensation factors

      TD = 1.0_dp / (FCONCN(LCELL) + FCONCA(LCELL))
      FCONCN(LCELL) = TD * FCONCN(LCELL)
      FCONCA(LCELL) = TD * FCONCA(LCELL)
      !bs
      TD = 1.0_dp / (FCONCN_ORG(LCELL) + FCONCA_ORG(LCELL))
      FCONCN_ORG(LCELL) = TD * FCONCN_ORG(LCELL)
      FCONCA_ORG(LCELL) = TD * FCONCA_ORG(LCELL)
      !bs

      !*** note CHEMRAT includes species other than sulfate.

!      CHEMRAT     = SO4FAC * SO4RAT(LCELL)         ! [mom3 m-3 s-1]
!sorgam      CHEMRAT_ORG = ORGFAC * ( ORGARO1RAT(LCELL)
!sorgam                               + ORGARO2RAT(LCELL)
!sorgam                               + ORGALK1RAT(LCELL)
!sorgam                               + ORGOLE1RAT(LCELL)
!sorgam                               + ORGBIO1RAT(LCELL)
!sorgam                               + ORGBIO2RAT(LCELL)
!sorgam                               + ORGBIO3RAT(LCELL)
!sorgam                               + ORGBIO4RAT(LCELL) )
!sorgam                                                ! [mom3 m-3 s-1]
      CHEMRAT_ORG = ORGFAC * SOA_MADE(LCELL)    ! [mom3 m-3 s-1]

      ! *** Calculate the production rates for new particles

      CGRN3(LCELL) = SO4FAC * DMDT(LCELL) ! Rate of increase of 3rd moment

!      CHEMRAT = CHEMRAT - CGRN3(LCELL ) !bs Remove the rate of new
!                                        !bs particle 3rd moment
!                                        !bs production from CHEMRAT.
!
!      CHEMRAT= MAX(CHEMRAT, 0.0_dp) ! Prevent CHEMRAT from being negative
!
      ! *** Now calculate the rate of condensation on existing particles.

!      CGRN3(LCELL) =  CGRN3(LCELL) + CHEMRAT     * FCONCN(LCELL)  &
!sorgam                             + CHEMRAT_ORG * FCONCN_ORG(LCELL)
!
!      CGRA3(LCELL) =    CHEMRAT     * FCONCA(LCELL)     &
!                      + CHEMRAT_ORG * FCONCA_ORG(LCELL)

      CONDRATE = SO4FAC * DELTASO4A(LCELL) / DT ! [mom3 m**-3 s-1]

      CGRN3(LCELL) = CGRN3(LCELL) + CONDRATE * FCONCN(LCELL) &
                     + CHEMRAT_ORG * FCONCN_ORG(LCELL)
      CGRA3(LCELL) = CONDRATE * FCONCA(LCELL) + CHEMRAT_ORG * FCONCA_ORG(LCELL)

   END DO  !  end 2nd loop over NUMCELLS

END SUBROUTINE MADE_NUCLCOND

!------------------------------------------------------------------------------!

SUBROUTINE made_nuklefit( blkta, blkrh, h2so4, so4rat, blksize, numcells, &
                          ndot1, mdot1 )

   ! NUKLEFIT:
   ! ---------
   !
   ! Calculates binary nucleation rate using revised theory,
   ! stauffer+binder&stauffer kinetics and noppel hydrate correction.
   !
   ! Reference:
   ! ----------
   !
   ! Vehkam?ki, H.; Kulmala, M.; Napari, I.; Lehtinen, K. E. J.;
   ! Timmreck, C.; Noppel, M.; Laaksonen, A.: An improved parameterization
   ! for sulfuric acid-water nucleation rates for tropospheric and
   ! stratospheric conditions, J. Geophys. Res. 10.1029/2002JD002184,
   ! 19 November 2002.
   !
   !   Authors:
   !   ---------
   !   C. TIMMRECK, MPI HAMBURG                                           2002
   !   A. Lauer, DLR - modifications for improved performance and use
   !                   with ECHAM/MADE
   !
   ! Input:
   ! ------
   !
   ! t     temperature [K]
   ! rh    relative humidity %/100 (100% relative humidity  means rh=1)
   ! rhoa  concentration of h2so4 vapour [1/m3]
   !
   ! Output:
   ! -------
   !
   ! x     mole fraction in the core of the critical cluster
   ! nwtot total number of water molecules in the critical cluster
   ! natot total number of h2so4 molecules in the critical cluster
   ! rc    radius of the critical cluster core [m]
   ! jnuc  nucleation rate [1/m3/s]
   !
   ! ... validity:
   !
   !     t    = 190.15K - 300.15K
   !     rh   = 0.0001 - 1 (0.01% - 100%)
   !     rhoa = 1e10 - 1e17 /m3 (1e4 - 1e11/cm3)
   !     jnuc = 1e-1 - 1e16/(m3s)  1e-7 - 1e10 /(cm3s) 
   !     ntot >= 4 ! make sure that you check also this, otherwise results
   !               ! can be rubbish!!!

   IMPLICIT NONE
   INTRINSIC EXP, LOG, MAX, MIN

   ! *** input ***

   integer,  intent(in)  :: blksize         ! dimension of arrays
   integer,  intent(in)  :: numcells        ! # of cells

   real(dp), intent(in)  :: blkta(blksize)  ! temperature [K]
   real(dp), intent(in)  :: blkrh(blksize)  ! rel. humidity [frac]
   real(dp), intent(in)  :: h2so4(blksize)  ! H2SO4(g) [ug/m3]
   real(dp), intent(in)  :: so4rat(blksize) ! maximum mass production rate
                                            ! allowed [ug(SO4)/m3/s]

   ! *** output ***

   real(dp), intent(out) :: ndot1(blksize)  ! nucleation rate [#/m3/s]
   real(dp), intent(out) :: mdot1(blksize)  ! so4 mass production rate
                                            ! [ug/m3/s]

   ! *** local ***

   double precision rhoa  ! H2SO4(g) [molecules/m**3]
   double precision t     ! temperature [k]
   double precision rh    ! relative humidity [frac]
!   double precision rhotres
!   double precision nwc, rc
   double precision x, nac, jnuc
   double precision ntot
   double precision logrh, logrh2, logrh3
   double precision t2, t3
   double precision logrhoa, logrhoa2, logrhoa3

   real(dp) :: mp         ! mass of sulfate in a 3.5 nm particle

   integer :: ncell       ! loop index

   !*** constants ***

   ! conversion factor molecules ---> ug
   double precision, parameter :: molec2ug = mwh2so4 * 1.0e6_dp / avo 
   ! conversion factor ug ---> molecules
   double precision, parameter :: ug2molec = avo * 1.0e-6_dp / mwh2so4

   real(dp) :: mp35       ! arithmetic statement function to compute
                          ! the mass of sulfate in a 3.5 nm particle
   ! coefficients for cubic in mp35
   real(dp), parameter :: a0 =  1.961385e2_dp
   real(dp), parameter :: a1 = -5.564447e2_dp
   real(dp), parameter :: a2 =  8.828801e2_dp
   real(dp), parameter :: a3 = -5.231409e2_dp

   ! *** functions ***

   ! *** function for the mass of sulfate in a 3.5 nm sphere
   ! *** obtained from a fit to the number of sulfate monomers in
   !     a 3.5 nm particle. Uses data from Nair & Vohra.

   double precision rr ! dummy variable for statement function
   mp35(rr) = molec2ug * (a0 + rr * (a1 + rr * (a2 + rr * a3)))

   ! ... code start ...

   DO ncell = 1, numcells
      ndot1(ncell) = 0.0
      mdot1(ncell) = 0.0
   END DO

   DO ncell = 1, numcells

      ! convert H2SO4(g) conc. from [ug/m3] to [molec/m3]

      rhoa = h2so4(ncell) * ug2molec

      ! ensure valid parameter range:

      ! temperature
      t = max(190.15_dp, blkta(ncell))
      t = min(t, 300.15_dp)
      ! relative humidity
      rh = max(0.0001_dp, blkrh(ncell))
      rh = min(rh, 1.0_dp)
      ! h2so4 concentration
      rhoa = max(1.0e10_dp, rhoa)
      rhoa = min(rhoa, 1.0e17_dp)

      ! calculate variables used several times

      logrh  = log(rh)
      logrh2 = logrh  * logrh
      logrh3 = logrh2 * logrh

      t2 = t * t
      t3 = t2 * t

      logrhoa  = log(rhoa/1.0e6_dp)
      logrhoa2 = logrhoa  * logrhoa
      logrhoa3 = logrhoa2 * logrhoa

      ! start nucleation calculations

      ! mole fraction of sulfuric acid in the critical cluster

      x    =  0.7409967177282139_dp                       - 0.002663785665140117_dp*t                  &
            + 0.002010478847383187_dp*logrh               - 0.0001832894131464668_dp*t*logrh           &
            + 0.001574072538464286_dp*logrh2              - 0.00001790589121766952_dp*t*logrh2         &
            + 0.0001844027436573778_dp*logrh3             - 1.503452308794887e-6_dp*t*logrh3           &
            - 0.003499978417957668_dp*logrhoa             + 0.0000504021689382576_dp*t*logrhoa

      ! nucleation rate

      jnuc =  0.1430901615568665_dp                       + 2.219563673425199_dp*t                     &
            - 0.02739106114964264_dp*t2                   + 0.00007228107239317088_dp*t3               &
            + 5.91822263375044_dp/x                       + 0.1174886643003278_dp*logrh                &
            + 0.4625315047693772_dp*t*logrh               - 0.01180591129059253_dp*t2*logrh            &
            + 0.0000404196487152575_dp*t3*logrh           + (15.79628615047088_dp*logrh)/x             &
            - 0.215553951893509_dp*logrh2                 - 0.0810269192332194_dp*t*logrh2             &
            + 0.001435808434184642_dp*t2*logrh2           - 4.775796947178588e-6_dp*t3*logrh2          &
            - (2.912974063702185_dp*logrh2)/x             - 3.588557942822751_dp*logrh3                &
            + 0.04950795302831703_dp*t*logrh3             - 0.0002138195118737068_dp*t2*logrh3         &
            + 3.108005107949533e-7_dp*t3*logrh3           - (0.02933332747098296_dp*logrh3)/x          &
            + 1.145983818561277_dp*logrhoa                - 0.6007956227856778_dp*t*logrhoa            &
            + 0.00864244733283759_dp*t2*logrhoa           - 0.00002289467254710888_dp*t3*logrhoa       &
            - (8.44984513869014_dp*logrhoa)/x             + 2.158548369286559_dp*logrh*logrhoa         &
            + 0.0808121412840917_dp*t*logrh*logrhoa       - 0.0004073815255395214_dp*t2*logrh*logrhoa  &
            - 4.019572560156515e-7_dp*t3*logrh*logrhoa    + (0.7213255852557236_dp*logrh*logrhoa)/x    &
            + 1.62409850488771_dp*logrh2*logrhoa          - 0.01601062035325362_dp*t*logrh2*logrhoa    &
            + 0.00003771238979714162_dp*t2*logrh2*logrhoa + 3.217942606371182e-8_dp*t3*logrh2*logrhoa  &
            - (0.01132550810022116_dp*logrh2*logrhoa)/x   + 9.71681713056504_dp*logrhoa2               &
            - 0.1150478558347306_dp*t*logrhoa2            + 0.0001570982486038294_dp*t2*logrhoa2       &
            + 4.009144680125015e-7_dp*t3*logrhoa2         + (0.7118597859976135_dp*logrhoa2)/x         &
            - 1.056105824379897_dp*logrh*logrhoa2         + 0.00903377584628419_dp*t*logrh*logrhoa2    &
            - 0.00001984167387090606_dp*t2*logrh*logrhoa2 + 2.460478196482179e-8_dp*t3*logrh*logrhoa2  &
            - (0.05790872906645181_dp*logrh*logrhoa2)/x   - 0.1487119673397459_dp*logrhoa3             &
            + 0.002835082097822667_dp*t*logrhoa3          - 9.24618825471694e-6_dp*t2*logrhoa3         &
            + 5.004267665960894e-9_dp*t3*logrhoa3         - (0.01270805101481648_dp*logrhoa3)/x

      ! total number of molecules in the critical cluster

      ntot =-0.002954125078716302_dp                      - 0.0976834264241286_dp*t                    &
            + 0.001024847927067835_dp*t2                  - 2.186459697726116e-6_dp*t3                 &
            - 0.1017165718716887_dp/x                     - 0.002050640345231486_dp*logrh              &
            - 0.007585041382707174_dp*t*logrh             + 0.0001926539658089536_dp*t2*logrh          &
            - 6.70429719683894e-7_dp*t3*logrh             - (0.2557744774673163_dp*logrh)/x            &
            + 0.003223076552477191_dp*logrh2              + 0.000852636632240633_dp*t*logrh2           &
            - 0.00001547571354871789_dp*t2*logrh2         + 5.666608424980593e-8_dp*t3*logrh2          &
            + (0.03384437400744206_dp*logrh2)/x           + 0.04743226764572505_dp*logrh3              &
            - 0.0006251042204583412_dp*t*logrh3           + 2.650663328519478e-6_dp*t2*logrh3          &
            - 3.674710848763778e-9_dp*t3*logrh3           - (0.0002672510825259393_dp*logrh3)/x        &
            - 0.01252108546759328_dp*logrhoa              + 0.005806550506277202_dp*t*logrhoa          &
            - 0.0001016735312443444_dp*t2*logrhoa         + 2.881946187214505e-7_dp*t3*logrhoa         &
            + (0.0942243379396279_dp*logrhoa)/x           - 0.0385459592773097_dp*logrh*logrhoa        &
            - 0.0006723156277391984_dp*t*logrh*logrhoa    + 2.602884877659698e-6_dp*t2*logrh*logrhoa   &
            + 1.194163699688297e-8_dp*t3*logrh*logrhoa    - (0.00851515345806281_dp*logrh*logrhoa)/x   &
            - 0.01837488495738111_dp*logrh2*logrhoa       + 0.0001720723574407498_dp*t*logrh2*logrhoa  &
            - 3.717657974086814e-7_dp*t2*logrh2*logrhoa   - 5.148746022615196e-10_dp*t3*logrh2*logrhoa &
            + (0.0002686602132926594_dp*logrh2*logrhoa)/x - 0.06199739728812199_dp*logrhoa2            &
            + 0.000906958053583576_dp*t*logrhoa2          - 9.11727926129757e-7_dp*t2*logrhoa2         &
            - 5.367963396508457e-9_dp*t3*logrhoa2         - (0.007742343393937707_dp*logrhoa2)/x       &
            + 0.0121827103101659_dp*logrh*logrhoa2        - 0.0001066499571188091_dp*t*logrh*logrhoa2  &
            + 2.534598655067518e-7_dp*t2*logrh*logrhoa2   - 3.635186504599571e-10_dp*t3*logrh*logrhoa2 &
            + (0.0006100650851863252_dp*logrh*logrhoa2)/x + 0.0003201836700403512_dp*logrhoa3          &
            - 0.0000174761713262546_dp*t*logrhoa3         + 6.065037668052182e-8_dp*t2*logrhoa3        &
            - 1.421771723004557e-11_dp*t3*logrhoa3        + (0.0001357509859501723_dp*logrhoa3)/x

      ntot = exp(ntot)

      mdot1(ncell) = 0.0_dp
      ndot1(ncell) = 0.0_dp

      IF (ntot >= 4.0_dp) THEN
         ! number of sulfuric acid molecules in the cluster
         nac = x * ntot
!         ! number of water molecules in the cluster
!         ! nwc =(1.-x)*ntot
!          ! radius of the cluster in nm
!          rc=exp(-1.6524245+0.42316402*x+0.33466487*log(ntot)) ! nm
!          rc=rc*1.e-9 ! m
          ndot1(ncell) = exp(jnuc) * 1.0e6_dp           ! #/m3/s
!          mdot1(ncell) = ndot1(ncell) * nac * molec2ug ! ug/m3/s

         ! *** get the mass of sulfate in a 3.5 nm particle

         mp = mp35(rh) ! number of micrograms of sulfate
                       ! in a 3.5 nm particle at ambient RH

         ! assume log-normal distribution with d=3.5 nm,
         ! sigma = sigma(Aitken mode)
         mdot1(ncell) = ndot1(ncell) * mp * es36(akn)

         IF (mdot1(ncell) > so4rat(ncell)) THEN
            mdot1(ncell) = so4rat(ncell)     ! limit nucleated mass by
                                             ! available mass
            ndot1(ncell) = mdot1(ncell) / &  ! adjust DNDT to this
!                           (nac * molec2ug)
                           (mp * es36(akn))
         ENDIF

         IF (mdot1(ncell) == 0.0_dp) ndot1(ncell) = 0.0_dp
      END IF ! if (ntot >= 4.0_dp)

!      ! treshold concentration of h2so4 (1/cm^3) which produces
!      ! nucleation rate 1/(cm^3 s) as a function of rh and t
!
!      rhotres=exp( -279.2430007512709_dp             &
!                   + 11.73439886096903_dp*rh         &
!                   + 22700.92970508331_dp/t          &
!                   - (1088.644983466801_dp*rh)/t     &
!                   + 1.144362942094912_dp*t          &
!                   - 0.03023314602163684_dp*rh*t     &
!                   - 0.001302541390154324_dp*t2      &
!                   - 6.386965238433532_dp*logrh      &
!                   + (854.980361026715_dp*logrh)/t   &
!                   + 0.00879662256826497_dp*t*logrh)

   END DO

END SUBROUTINE made_nuklefit

!------------------------------------------------------------------------------!

SUBROUTINE MADE_AEROSTEP ( BLKSIZE, NUMCELLS,                       &
                           CBLK, DT,                                &
!sorgam                    ORGARO1RAT, ORGARO2RAT,                  &
!sorgam                    ORGALK1RAT, ORGOLE1RAT,                  &
!sorgam                    ORGBIO1RAT, ORGBIO2RAT,                  &
!sorgam                    ORGBIO3RAT, ORGBIO4RAT,                  &
!unused                    EPM25I, EPM25J, EORGI,EORGJ, EECI, EECJ, &
!unused                    ESOIL, ESEAS, EPMCOARSE,                 &
                           DGNUC, DGACC,                            &
                           FCONCN, FCONCA,                          &
                           FCONCN_ORG, FCONCA_ORG,                  &
!unused                    PMASSN, PMASSA, PMASSC,                  &
                           DMDT, DNDT, DELTASO4A, SOA_MADE,         &
                           URN00, URA00,                            &
                           BRNA01, C30,                             &
                           CGRN3, CGRA3,                            &
                           MERGENUM, MERGEM3 )

   !***********************************************************************
   !
   !      NOTE:
   !
   ! ***  DESCRIPTION: Integrate the Number and Mass equations
   !                   for each mode over the time interval DT.
   !
   !      PRECONDITIONS:
   !       AEROSTEP() must follow calls to all other dynamics routines.
   !
   ! ***   Revision history:
   !       Adapted  3/95 by UAS and CJC from EAM2's code.
   !       Revised  7/29/96 by FSB to use block structure
   !       Revised 11/15/96 by FSB dropped flow-through and cast
   !                        number solver into Riccati equation form.
   !       Revised  8/8/97  by FSB to have mass in Aitken and accumulation
   !                        modes each predicted rather than total mass and
   !                        Aitken mode mass. Also used a local
   !                        approximation for the error function. Also
   !                        added coarse mode.
   !       Revised  9/18/97 by FSB to fix mass transfer from Aitken to
   !                        accumulation mode by coagulation
   !       Revised 10/27/97 by FSB to modify code to use primay emissions
   !                        and to correct 3rd moment updates.
   !                        Also added coarse mode.
   !       Revised 11/4/97  by FSB to fix error in other anthropogenic PM2.5
   !       Revised 11/5/97  by FSB to fix error in MSTRNSFR
   !       Revised 11/6/97  FSB to correct the expression for FACTRANS to
   !                        remove the 6/pi coefficient. UAS found this.
   !       Revised 12/15/97 by FSB to change equations for mass
   !                        concentration to a chemical production form
   !                        with analytic solutions for the Aitken mode
   !                        and to remove time stepping of the 3rd moments.
   !                        The mass conc. in the accumulation mode is
   !                        updated with a forward Euler step.
   !       Revised 1/6/98   by FSB Lowered minimum concentration for
   !                        sulfate aerosol to 0.1 [ ng / m**3 ].
   !       Revised 1/12/98  C30 replaces BRNA31 as a variable. C30
   !                        represents intermodal transfer rate of 3rd
   !                        moment in place of 3rd moment coagulation rate.
   !       Revised 6/5/98   added new renaming criterion
   !       Added   7/6/99   by BS condensational groth factors for organics
   !
   !**********************************************************************

   IMPLICIT NONE
   INTRINSIC EXP, LOG, MAX, SIGN, SQRT

   ! *** ARGUMENTS ***

   INTEGER,  INTENT(in)    :: BLKSIZE               ! dimension of arrays
   INTEGER,  INTENT(in)    :: NUMCELLS              ! actual number of cells
                                                    ! in arrays
   REAL(dp), INTENT(inout) :: CBLK(BLKSIZE,NSPCSDA) ! main array of variables
   REAL(dp), INTENT(in)    :: DT                    ! time step [s]

   ! *** Chemical production rates: [ug/m3/s]

!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     aromatics
!sorgam   REAL(dp) :: ORGARO1RAT(BLKSIZE)
!sorgam   REAL(dp) :: ORGARO2RAT(BLKSIZE)
!sorgam
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     alkanes & others
!sorgam   REAL(dp) :: ORGALK1RAT(BLKSIZE)
!sorgam   REAL(dp) :: ORGOLE1RAT(BLKSIZE)
!sorgam
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO1RAT(BLKSIZE)
!sorgam   REAL(dp) :: ORGBIO2RAT(BLKSIZE)
!sorgam   REAL(dp) :: ORGBIO3RAT(BLKSIZE)
!sorgam   REAL(dp) :: ORGBIO4RAT(BLKSIZE)

!unused   ! *** Primary emissions rates: [ug/m3/s]
!unused
!unused   ! *** emissions rates for unidentified PM2.5 mass
!unused   REAL(dp) :: EPM25I(BLKSIZE)   ! Aitken mode
!unused   REAL(dp) :: EPM25J(BLKSIZE)   ! accumululaton mode
!unused
!unused   ! *** emissions rates for primary organic aerosol
!unused   REAL(dp) :: EORGI(BLKSIZE)    ! Aitken mode
!unused   REAL(dp) :: EORGJ(BLKSIZE)    ! accumululaton mode
!unused
!unused   ! *** emissions rates for elemental carbon
!unused   REAL(dp) :: EECI(BLKSIZE)     ! Aitken mode
!unused   REAL(dp) :: EECJ(BLKSIZE)     ! accumululaton mode
!unused
!unused   ! *** emissions rates for coarse mode particles
!unused   REAL(dp) :: ESOIL(BLKSIZE)           ! soil derived coarse aerosols
!unused   REAL(dp) :: ESEAS(BLKSIZE)           ! marine coarse aerosols
!unused   REAL(dp) :: EPMCOARSE(BLKSIZE)       ! anthropogenic coarse aerosols

   REAL(dp), INTENT(in)    :: DGNUC(BLKSIZE)      ! Aitken mode mean
                                                  ! diameter [m]
   REAL(dp), INTENT(in)    :: DGACC(BLKSIZE)      ! acc. mode mean diameter [m]
   REAL(dp), INTENT(in)    :: FCONCN(BLKSIZE)     ! reciprocal condensation
                                                  ! rate, Aitken mode [1/s]
   REAL(dp), INTENT(in)    :: FCONCA(BLKSIZE)     ! reciprocal condensation
                                                  ! rate, acc. mode [1/s]
   REAL(dp), INTENT(in)    :: FCONCN_ORG(BLKSIZE) ! reciprocal condensation
                                                  ! rate for organics,
                                                  ! Aitken mode [1/s]
   REAL(dp), INTENT(in)    :: FCONCA_ORG(BLKSIZE) ! reciprocal condensation
                                                  ! rate for organics,
                                                  ! acc. mode [1/s]
   REAL(dp), INTENT(in)    :: DMDT(BLKSIZE)       ! rate of production of new
                                                  ! so4 mass by particle
                                                  ! formation [ug/m3/s]
   REAL(dp), INTENT(in)    :: DNDT(BLKSIZE)       ! rate of production of new
                                                  ! particle  number by
                                                  ! particle formation [#/m3/s]
   REAL(dp), INTENT(in)    :: DELTASO4A(BLKSIZE)  ! increment of concentration
                                                  ! added to sulfate aerosol by
                                                  ! condensation [ug/m3]
   REAL(dp), INTENT(in)    :: SOA_MADE(BLKSIZE)   ! SOA (gas phase) "emissions"
                                                  ! [ug/m3/s]
   REAL(dp), INTENT(out)   :: MERGENUM(BLKSIZE)   ! mode merging, number [#/m3]
                                                  ! (for diagnostics only)
   REAL(dp), INTENT(out)   :: MERGEM3(BLKSIZE)    ! mode merging, 3rd moment
                                                  ! [3rd mom./m3]
                                                  ! (for diagnostics only)
   REAL(dp), INTENT(in)    :: URN00(BLKSIZE)      ! Aitken mode intramodal
                                                  ! coagulation rate
   REAL(dp), INTENT(in)    :: URA00(BLKSIZE)      ! acc. mode intramodal
                                                  ! coagulation rate
   REAL(dp), INTENT(in)    :: BRNA01(BLKSIZE)     ! bimodal coagulation rate
                                                  ! (Aitken <--> acc.), number
   REAL(dp), INTENT(in)    :: C30(BLKSIZE)        ! intermodal 3rd moment
                                                  ! transfer rate by intermodal
                                                  ! coagulation
   REAL(dp), INTENT(in)    :: CGRN3(BLKSIZE)      ! growth rate for 3rd moment,
                                                  ! Aitken mode [3rd mom/m3/s]
   REAL(dp), INTENT(inout) :: CGRA3(BLKSIZE)      ! growth rate for 3rd moment,
                                                  ! acc. mode [3rd mom/m3/s]

!unused   ! *** Modal mass concentrations [ug/m3]
!unused
!unused   REAL(dp) :: PMASSN(BLKSIZE) ! mass concentration in Aitken mode
!unused   REAL(dp) :: PMASSA(BLKSIZE) ! mass concentration in accumulation mode
!unused   REAL(dp) :: PMASSC(BLKSIZE) ! mass concentration in coarse mode

   ! *** Local Variables ***

   INTEGER :: L, LCELL, SPC   ! Loop indices

   ! ** following scratch variables are used for solvers

   ! *** variables needed for modal dynamics solvers:

   REAL(dp) :: A, B, C
   REAL(dp) :: M1, M2, Y0, Y
   REAL(dp) :: DHAT, P, PEXPDT, EXPDT
   REAL(dp) :: LOSS, PROD, POL, LOSSINV
   REAL(dp) :: MSTRNSFR     ! mass intermodal transfer by coagulation
   REAL(dp) :: FACTRANS     ! special factor to compute mass transfer

   ! *** CODE additions for renaming
   REAL(dp) :: AAA, XNUM, XM3, XXM3, FNUM, FM3, PHNUM, PHM3 ! Defined below
   REAL(dp) :: ERF, ERFC  ! Error and complementary error function
   REAL(dp) :: XX         ! dummy argument for ERF and ERFC

   REAL(dp), PARAMETER :: CONMIN = 1.0E-30_dp ! a numerical value for a minimum
                                              ! concentration

   !al * variables for inlining of the function GETAF
   REAL(dp) :: AA, BB, CC, DISC, QQ, alfa, GETAF_L, yji

   LOGICAL :: FIRSTIME
   DATA       FIRSTIME / .TRUE. /

   SAVE FIRSTIME, XXM3

   !     :::::::::::::::::::::::::::::::::::::
   ! *** Statement function given for error function. Source is
   !     Meng, Z., and J.H.Seinfeld (1994) On the source of the
   !     submicrometer droplet mode of urban and regional aerosols.
   !     Aerosol Sci. and Technology, 20:253-265. They cite Reasearch &
   !     Education Asociation (REA), (1991) Handbook of Mathematical,
   !     Scientific, and Engineering Formulas, Tables, Functions, Graphs,
   !     Transforms: REA, Piscataway, NJ. p. 493.

   ERF(XX)  = SQRT(1.0_dp - EXP(-4.0_dp * XX * XX / PI))
   ERFC(XX) = 1.0_dp - ERF(XX)
   !     ::::::::::::::::::::::::::::::::::::::::

   ! ///// begin code

   IF (FIRSTIME) THEN
      XXM3 = 3.0_dp * XXLSG(akn) / SQRT2 ! factor used in error function
                                         ! call below
      FIRSTIME = .FALSE.
   END IF  ! firstime

   ! *** set up time-step integration

   DO L = 1, NUMCELLS

      ! *** code to move number forward by one time step.
      ! *** solves the Riccati equation:

      !     dY/dt = C - A * Y ** 2 - B * Y

      !     Coded 11/21/96 by Dr. Francis S. Binkowski

      ! *** Aitken mode:

      ! *** coefficients

      A = URN00(L)
      B = BRNA01(L) * CBLK(L,VAC0)
      C = DNDT(L)
!unused    + FACTNUMN * (ANTHFAC * (EPM25I(L) + EECI(L)) + &
!unused                  ORGFAC * EORGI(L)) ! includes primary emissions

      Y0 = CBLK(L,VNU0) ! initial condition

      ! ***  trap on C = 0

      ! numerically unstable on NEC SX-4 !!!
      ! IF(C > 0.0D0) THEN

      IF (C > 1.0e-19_dp) THEN
         DHAT = SQRT(B * B + 4.0_dp * A * C) !SQRT(B * B + 4.0D0 * A * C)

         M1 = 2.0_dp * A * C / (B + DHAT) !2.0D0 * A * C / (B + DHAT)
         M2 = - 0.5_dp * (B + DHAT) !- 0.5D0 * (B + DHAT)

         P =  -(M1 - A  * Y0) / (M2 - A * Y0)
         PEXPDT = P * EXP(-DHAT * DT)

!         Y = (M1 + M2 * PEXPDT) / (A * (1.0D0 + PEXPDT)) ! solution
         Y = (M1 + M2 * PEXPDT) / (A * (1.0_dp + PEXPDT)) ! solution
      ELSE

         ! *** rearrange solution for NUMERICAL stability
         !     note If B << A * Y0, the following form, although
         !     seemingly awkward gives the correct answer.

         EXPDT = EXP(-B * DT)
         IF (EXPDT < 1.0_dp) THEN !IF (EXPDT < 1.0D0) THEN
!            Y = B * Y0 * EXPDT / (B + A * Y0 * (1.0D0 - EXPDT))
            Y = B * Y0 * EXPDT / (B + A * Y0 * (1.0_dp - EXPDT))
         ELSE
            Y = Y0
         END IF
      END IF

      CBLK(L,VNU0) = MAX(NUMMIN(akn),Y)

      ! *** now do accumulation mode number

      ! *** coefficients

      A = URA00(L)
!      B = 0.0D0 ! NOTE B = 0.0
      B = 0.0_dp ! NOTE B = 0.0
!      C = 0.0D0
      C = 0.0_dp
!unused   FACTNUMA * (ANTHFAC * (EPM25J(L) + EECJ(L)) &
!unused             + ORGFAC * EORGJ(L)) ! includes primary emissions

      Y0 = CBLK(L,VAC0) ! initial condition

!unused   ! *** this equation requires special handling, because C can be zero.
!unused   !     if this happens, the form of the equation is different:
!unused
!unused   IF (C > 0.0D0) THEN
!unused
!unused      DHAT = SQRT(4.0D0_dp * A * C)
!unused
!unused      M1 = 2.0D0_dp * A * C / DHAT
!unused      M2 = - 0.5D0_dp * DHAT
!unused
!unused      P = -(M1 - A  * Y0) / (M2 - A * Y0)
!unused      PEXPDT = P * EXP(-DHAT * DT)
!unused
!unused      Y = (M1 + M2 * PEXPDT) / (A * (1.0D0_dp + PEXPDT)) ! solution
!unused
!unused   ELSE

!      Y = Y0 / (1.0D0 + DT * A * Y0) ! correct solution to
      Y = Y0 / (1.0_dp + DT * A * Y0) ! correct solution to
                                      ! equation for C = 0.0
!unused   END IF

      CBLK(L,VAC0) = MAX(NUMMIN(acc),Y)

      ! *** now do coarse mode number neglecting coagulation

      !     emissions are handled by ECHAM SUBROUINE XTEMISS
      !     ---> nothing to do for coarse mode

!unused   PROD = SOILFAC * ESOIL(L) + SEASFAC * ESEAS(L) &
!unused          + ANTHFAC * EPMCOARSE(L)
!unused   CBLK(L,VCORN) = CBLK(L,VCORN) + FACTNUMC * PROD * DT

      ! *** Prepare to advance modal mass concentration one time step.

      ! *** Set up production and and intermodal transfer terms terms:

      ! Emissions are handled by ECHAM SUBROUTINE XTEMISS.

!unused   CGRN3(L) = CGRN3(L) &
!unused              + ANTHFAC * (EPM25I(L) + EECI(L))
!unused              + ORGFAC * EORGI(L)   ! includes growth from primary
!unused                                    ! emissions

      CGRA3(L) = CGRA3(L) + C30(L) ! include growth from transfer of 3rd
!unused                            ! moment by intermodal coagulation
!unused          + ANTHFAC * (EPM25J(L) + EECJ(L))
!unused          + ORGFAC * EORGJ(L) ! includes growth from primary
!unused                              ! emissions and transfer of 3rd
!unused                              ! moment by intermodal coagulation

      ! *** set up transfer coefficients for coagulation between Aitken and
      !     accumulation modes

      ! *** set up special factors for mass transfer from the Aitken to
      !     accumulation modes by intermodal coagulation. The mass transfer
      !     rate is proportional to the 3rd moment, transfer rate, C30. The
      !     proportionality factor is p/6 times the the average particle
      !     density. The average particle density for a species is the
      !     species mass concentration divided by the particle volume
      !     concentration, pi/6 times the 3rd moment concentration.
      !     The p/6 coefficients cancel.

      LOSS =  C30(L) / CBLK(L,VNU3) ! Normalized coagulation
                                    ! transfer rate coefficient

!unused   EXPDT = EXP(-FACTRANS)    ! variable name is re-used here. This
!unused                             ! exponential decay term is common to
!unused                             ! all Aitken mode mass equations

      EXPDT = EXP(-LOSS*DT)
      FACTRANS = 1.0_dp - EXPDT  ! Since no production term has to be
                                 ! taken into account (emissions are
                                 ! NOT handled by SUBROUTINE MADE_AEROSTEP),
                                 ! NO estimate but the EXACT amount of
                                 ! mass transferred from the Atiken to
                                 ! the accumulation mode can be
                                 ! calculated!

      LOSSINV = 1.0_dp / LOSS    ! set up for multiplication rather than
                                 ! division.

      ! *** now advance mass concentrations one time step.

      ! *** Solve Aitken-mode equations of form: dc/dt = P - L*c
      ! *** Solution is: c(t0 + dt) = p/L + (c(0) - P/L) * exp(-L*dt)

      ! *** sulfate:

      MSTRNSFR = CBLK(L,VSO4AI) * FACTRANS
      PROD = DELTASO4A(L) * FCONCN(L) / DT + DMDT(L) ! Condensed mass
                                                     ! + new particle mass
      POL = PROD * LOSSINV

      CBLK(L,VSO4AI) = POL + (CBLK(L,VSO4AI) - POL) * EXPDT
      CBLK(L,VSO4AI) = MAX(MMINSO4(akn), CBLK(L,VSO4AI))         

      CBLK(L,VSO4AJ) = CBLK(L,VSO4AJ) + DELTASO4A(L) * FCONCA(L) + MSTRNSFR
      CBLK(L,VSO4AJ) = MAX(MMINSO4(acc), CBLK(L,VSO4AJ))         

      ! *** Update sulfuric acid vapor concentration by removing mass
      !     concentration of condensed sulfate and newly produced particles.
      ! *** The method follows Youngblood and Kreidenweis, Further Development
      !     and Testing of a Bimodal Aerosol Dynamics Model, Colorado State
      !     University Department of Atmospheric Science Paper Number 550,
      !     April,1994, pp 85-89.

      CBLK(L,VSULF) = MAX(CONMIN, CBLK(L,VSULF)                   &
                                  - (DELTASO4A(L) + DMDT(L) * DT) &
                                  * MWH2SO4/MWSO4)

!sorgam      ! *** anthropogenic secondary organic:
!sorgam      !bs * anthropogenic secondary organics from aromatic precursors
!sorgam
!sorgam      MSTRNSFR = CBLK(L,VORGARO1I) * FACTRANS
!sorgam      PROD = ORGARO1RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSINV
!sorgam
!sorgam      CBLK(L,VORGARO1I) = POL + (CBLK(L,VORGARO1I) - POL) * EXPDT
!sorgam
!sorgam      CBLK(L,VORGARO1I) = MAX(CONMIN, CBLK(L,VORGARO1I))
!sorgam
!sorgam      CBLK(L,VORGARO1J) = CBLK(L,VORGARO1J)                    &
!sorgam                          + ORGARO1RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                          + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(L,VORGARO1J) = MAX(CBLK(L,VORGARO1J),1.e-30_dp)
!sorgam
!sorgam      !bs * second species from aromatics
!sorgam      MSTRNSFR = CBLK(L,VORGARO2I) * FACTRANS
!sorgam      PROD = ORGARO2RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSINV
!sorgam
!sorgam      CBLK(L,VORGARO2I) = POL + (CBLK(L,VORGARO2I) - POL) * EXPDT
!sorgam
!sorgam      CBLK(L,VORGARO2I) = MAX(CONMIN, CBLK(L,VORGARO2I))
!sorgam
!sorgam      CBLK(L,VORGARO2J) = CBLK(L,VORGARO2J)                    &
!sorgam                          + ORGARO2RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                          + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(L,VORGARO2J) = MAX(CBLK(L,VORGARO2J),1.e-30_dp)
!sorgam 
!sorgam      !bs * anthropogenic secondary organics from alkanes & other
!sorgam      !bs   precursors
!sorgam      !bs * higher alkanes
!sorgam      MSTRNSFR = CBLK(L,VORGALK1I) * FACTRANS
!sorgam      PROD = ORGALK1RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSINV
!sorgam
!sorgam      CBLK(L,VORGALK1I) = POL + (CBLK(L,VORGALK1I) - POL) * EXPDT
!sorgam
!sorgam      CBLK(L,VORGALK1I) = MAX(CONMIN, CBLK(L,VORGALK1I))
!sorgam
!sorgam      CBLK(L,VORGALK1J) = CBLK(L,VORGALK1J)                    &
!sorgam                          + ORGALK1RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                          + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(L,VORGALK1J) = MAX(CBLK(L,VORGALK1J),1.e-30_dp)
!sorgam
!sorgam      !bs * higher olefines
!sorgam      MSTRNSFR = CBLK(L,VORGOLE1I) * FACTRANS
!sorgam      PROD = ORGOLE1RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSINV
!sorgam
!sorgam      CBLK(L,VORGOLE1I) = POL + (CBLK(L,VORGOLE1I) - POL) * EXPDT
!sorgam
!sorgam      CBLK(L,VORGOLE1I) = MAX(CONMIN, CBLK(L,VORGOLE1I))
!sorgam 
!sorgam      CBLK(L,VORGOLE1J) = CBLK(L,VORGOLE1J)                    &
!sorgam                          + ORGOLE1RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                          + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(L,VORGOLE1J) = MAX(CBLK(L,VORGOLE1J),1.e-30_dp)
!sorgam
!sorgam      ! *** biogenic secondary organic
!sorgam 
!sorgam      MSTRNSFR = CBLK(L,VORGBA1I) * FACTRANS
!sorgam      PROD = ORGBIO1RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSINV
!sorgam
!sorgam      CBLK(L,VORGBA1I) = POL + (CBLK(L,VORGBA1I) - POL) * EXPDT
!sorgam
!sorgam      CBLK(L,VORGBA1I) = MAX(CONMIN, CBLK(L,VORGBA1I))
!sorgam
!sorgam      CBLK(L,VORGBA1J) = CBLK(L,VORGBA1J)                     &
!sorgam                         + ORGBIO1RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                         + MSTRNSFR
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(L,VORGBA1J) = MAX(CBLK(L,VORGBA1J),1.e-30_dp)
!sorgam
!sorgam      !bs * second biogenic species
!sorgam      MSTRNSFR = CBLK(L,VORGBA2I) * FACTRANS
!sorgam      PROD = ORGBIO2RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSINV
!sorgam
!sorgam      CBLK(L,VORGBA2I) = POL + (CBLK(L,VORGBA2I) - POL) * EXPDT
!sorgam
!sorgam      CBLK(L,VORGBA2I) = MAX(CONMIN, CBLK(L,VORGBA2I))
!sorgam
!sorgam      CBLK(L,VORGBA2J) = CBLK(L,VORGBA2J)                     &
!sorgam                         + ORGBIO2RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                         + MSTRNSFR
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(L,VORGBA2J) = MAX(CBLK(L,VORGBA2J),1.e-30_dp)
!sorgam
!sorgam      !bs * third biogenic species
!sorgam      MSTRNSFR = CBLK(L,VORGBA3I) * FACTRANS
!sorgam      PROD = ORGBIO3RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSINV
!sorgam
!sorgam      CBLK(L,VORGBA3I) = POL + (CBLK(L,VORGBA3I) - POL) * EXPDT
!sorgam
!sorgam      CBLK(L,VORGBA3I) = MAX(CONMIN, CBLK(L,VORGBA3I))
!sorgam
!sorgam      CBLK(L,VORGBA3J) = CBLK(L,VORGBA3J)                     &
!sorgam                         + ORGBIO3RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                         + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(L,VORGBA3J) = MAX(CBLK(L,VORGBA3J),1.e-30_dp)
!sorgam
!sorgam      !bs * fourth biogenic species
!sorgam      MSTRNSFR = CBLK(L,VORGBA4I) * FACTRANS
!sorgam      PROD = ORGBIO4RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSINV
!sorgam
!sorgam      CBLK(L,VORGBA4I) = POL + (CBLK(L,VORGBA4I) - POL) * EXPDT
!sorgam
!sorgam      CBLK(L,VORGBA4I) = MAX(CONMIN, CBLK(L,VORGBA4I))
!sorgam
!sorgam      CBLK(L,VORGBA4J) = CBLK(L,VORGBA4J)                     &
!sorgam                         + ORGBIO4RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                         + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(L,VORGBA4J) = MAX(CBLK(L,VORGBA4J),1.e-30_dp)

      ! *** primary anthropogenic organic

      MSTRNSFR = CBLK(L,VORGPAI) * FACTRANS
      ! Emissions are handled by ECHAM SUBROUTINE XTEMISS.
!unused PROD = EORGI(L)
      PROD = SOA_MADE(L) * FCONCN_ORG(L)
      POL = PROD * LOSSINV

      CBLK(L,VORGPAI) = POL + (CBLK(L,VORGPAI) - POL) * EXPDT
      CBLK(L,VORGPAI) = MAX(mminorg(akn), CBLK(L,VORGPAI))

      CBLK(L,VORGPAJ) = CBLK(L,VORGPAJ)                    &
!unused                 + EORGJ(L) * DT                    &
                        + SOA_MADE(L) * FCONCA_ORG(L) * DT &
                        + MSTRNSFR
      CBLK(L,VORGPAJ) = MAX(MMINORG(acc),CBLK(L,VORGPAJ))

      ! *** sea salt

      MSTRNSFR = CBLK(L,VSEASI) * FACTRANS
      ! Emissions are handled by ECHAM SUBROUTINE XTEMISS.
      ! PROD = ESEASI(L)
      PROD = 0.0
      POL = PROD * LOSSINV

      CBLK(L,VSEASI) = POL + (CBLK(L,VSEASI) - POL) * EXPDT
      CBLK(L,VSEASI) = MAX(mminseas(akn), CBLK(L,VSEASI))

      CBLK(L,VSEASJ) = CBLK(L,VSEASJ)   &
!unused                + ESEASJ(L) * DT &
                       + MSTRNSFR
      CBLK(L,VSEASJ) = MAX(MMINSEAS(acc),CBLK(L,VSEASJ))

!unused      ! *** other anthropogenic PM2.5
!unused
!unused      MSTRNSFR = CBLK(L,VP25AI) * FACTRANS
!unused      ! Emissions are handled by ECHAM SUBROUTINE XTEMISS.
!unused      ! PROD = EPM25I(L)
!unused      PROD = 0.0
!unused      POL = PROD * LOSSINV
!unused
!unused      CBLK(L,VP25AI) = POL + (CBLK(L,VP25AI) - POL) * EXPDT
!unused
!unused      CBLK(L,VP25AI) = MAX(CONMIN, CBLK(L,VP25AI) )
!unused
!unused      CBLK(L,VP25AJ) = CBLK(L,VP25AJ) + EPM25J(L) * DT + MSTRNSFR

      ! ***  elemental carbon

      MSTRNSFR = CBLK(L,VECI) * FACTRANS
      ! Emissions are handled by ECHAM SUBROUTINE XTEMISS.
!unused PROD = EECI(L)
      PROD = 0.0
      POL = PROD * LOSSINV

      CBLK(L,VECI) = POL + (CBLK(L,VECI) - POL) * EXPDT
      CBLK(L,VECI) = MAX(MMINEC(akn), CBLK(L,VECI))

      CBLK(L,VECJ) = CBLK(L,VECJ)   &
!unused      &       + EECJ(L) * DT &
                     + MSTRNSFR
      CBLK(L,VECJ) = MAX(MMINEC(acc), CBLK(L,VECJ))

      ! ***  aerosol liquid water content

      MSTRNSFR = CBLK(L,VH2OAI) * FACTRANS

      !al Handle aerosol liquid water content the same way as the other
      !al mass tracers. This has to be done to ensure correct diagnostics
      !al of the aerosol liquid water content when writing to G3X-fields.

      PROD = 0.0_dp
      POL = PROD * LOSSINV

      CBLK(L,VH2OAI) = POL + (CBLK(L,VH2OAI) - POL) * EXPDT
      CBLK(L,VH2OAI) = MAX(1.0E-30_dp, CBLK(L,VH2OAI))

      CBLK(L,VH2OAJ) = CBLK(L,VH2OAJ) + MSTRNSFR
      CBLK(L,VH2OAJ) = MAX(1.0E-30_dp, CBLK(L,VH2OAJ))

      ! ***  coarse mode

      ! emissions are handled by ECHAM SUBROUTINE XTEMISS
      ! ---> nothing to do for coarse mode

!unused      ! *** soil dust
!unused
!unused      CBLK(L,VSOILA) = CBLK(L,VSOILA) + ESOIL(L) * DT
!unused      CBLK(L,VSOILA) = MAX(CONMIN, CBLK(L,VSOILA))
!unused 
!unused      ! *** sea salt
!unused
!unused      CBLK(L,VSEAS)  = CBLK(L,VSEAS) +  ESEAS(L) * DT
!unused      CBLK(L,VSEAS)  = MAX(CONMIN, CBLK(L,VSEAS))
!unused
!unused      ! *** anthropogenic PM10 coarse fraction
!unused
!unused      CBLK(L,VANTHA) = CBLK(L,VANTHA) + EPMCOARSE(L) * DT
!unused      CBLK(L,VANTHA) = MAX(CONMIN, CBLK(L,VANTHA))

   END DO   ! end of time-step loop for total mass

   ! *** Check for mode merging,if Aitken mode is growing faster than
   !     j-mode, then merge modes by renaming.

   ! *** use Binkowski-Kreidenweis paradigm, now including emissions

   IF (l_rename) THEN
   DO LCELL = 1, NUMCELLS

      IF(CGRN3(LCELL) > CGRA3(LCELL) .OR.   &
         DGNUC(LCELL) > 0.03e-6_dp .AND.       &
         CBLK(LCELL,VNU0) > CBLK(LCELL,VAC0)) THEN ! check if merging
                                                      ! necessary

            !al ******* inlining of:
            !al
            !al AAA = getaf(CBLK(LCELL,VNU0), &
            !al             CBLK(LCELL,VAC0), &
            !al             DGNUC(LCELL),     &
            !al             DGACC(LCELL),     &
            !al             XXLSG(akn),       &
            !al             XXLSG(acc),       &
            !al             SQRT2)

            alfa = XXLSG(akn) / XXLSG(acc)
            yji = log(DGACC(LCELL) / DGNUC(LCELL)) / (sqrt2 * XXLSG(akn))
            AA = 1.0_dp - alfa * alfa
            GETAF_L = log(alfa * CBLK(LCELL,VAC0) / CBLK(LCELL,VNU0))
            BB = 2.0_dp * yji * alfa * alfa
            CC = GETAF_L - yji * yji * alfa * alfa
            DISC = BB*BB - 4.0_dp * AA * CC
            if (DISC < 0.0_dp) then
               AAA = -5.0_dp       ! error in intersection
            else
               QQ = -0.5_dp * (BB + sign(1.0_dp, BB) * sqrt(DISC))
               AAA = CC / QQ
            end if

            !al ******* end inlining

            ! *** AAA is the value of ln(dd / DGNUC) / (SQRT2 * XXLSG(akn)),
            !     where dd is the diameter at which the Aitken mode and acc.
            !     mode number distributions intersect (overap).

            XNUM = MAX(AAA, XXM3) ! do not let XNUM become negative
                                  ! because this means that no more
                                  ! than one half of the total Aitken
                                  ! mode number may be transferred per
                                  ! call.

            XM3  = XNUM - XXM3    ! set up for 3rd moment and mass transfer

            IF (XM3 > 0.0_dp) THEN   ! do mode merging if overlap is correct

            PHNUM = 0.5_dp * (1.0_dp + ERF(XNUM))
            PHM3  = 0.5_dp * (1.0_dp + ERF(XM3))
            FNUM  = 0.5_dp * ERFC(XNUM)
            FM3   = 0.5_dp * ERFC(XM3)

            ! In the Aitken mode:

            ! *** FNUM and FM3 are the fractions of the number and 3rd moment
            !     distributions with  diameters greater than dd respectively.

            ! *** PHNUM and PHM3 are the fractions of the number and 3rd moment
            !     distributions with diameters less than dd.

            ! *** rename the Aitken mode particle number as accumulation mode
            !     particle number

            CBLK(LCELL,VAC0) = CBLK(LCELL,VAC0) + FNUM * CBLK(LCELL,VNU0)

            !al *** save mode merging, number (for diagnostics only)

            MERGENUM(LCELL) = FNUM * CBLK(LCELL,VNU0)

            ! *** adjust the Aitken mode number

            CBLK(LCELL,VNU0) = PHNUM * CBLK(LCELL,VNU0)

            ! *** Rename mass from Aitken mode to acumulation mode. The mass
            !     transferred to the accumulation mode is proportional to the
            !     amount of 3rd moment transferred, therefore FM3 is used for
            !     mass transfer.

            !al *** save mode merging, 3rd moment (for diagnostics only)

            MERGEM3(LCELL) = FM3 * CBLK(LCELL,VNU3)

            CBLK(LCELL,VSO4AJ)    = CBLK(LCELL,VSO4AJ)            &
                                    + CBLK(LCELL,VSO4AI) * FM3
            CBLK(LCELL,VNH4AJ)    = CBLK(LCELL,VNH4AJ)            &
                                    + CBLK(LCELL,VNH4AI) * FM3
            CBLK(LCELL,VNO3AJ)    = CBLK(LCELL,VNO3AJ)            &
                                    + CBLK(LCELL,VNO3AI) * FM3

!sorgam     CBLK(LCELL,VORGARO1J) = CBLK(LCELL,VORGARO1J)         &
!sorgam                             + CBLK(LCELL,VORGARO1I) * FM3
!sorgam     CBLK(LCELL,VORGARO2J) = CBLK(LCELL,VORGARO2J)         &
!sorgam                             + CBLK(LCELL,VORGARO2I) * FM3
!sorgam     CBLK(LCELL,VORGALK1J) = CBLK(LCELL,VORGALK1J)         &
!sorgam                             + CBLK(LCELL,VORGALK1I) * FM3
!sorgam     CBLK(LCELL,VORGOLE1J) = CBLK(LCELL,VORGOLE1J)         &
!sorgam                             + CBLK(LCELL,VORGOLE1I) * FM3
!sorgam     CBLK(LCELL,VORGBA1J)  = CBLK(LCELL,VORGBA1J)          &
!sorgam                             + CBLK(LCELL,VORGBA1I) * FM3
!sorgam     CBLK(LCELL,VORGBA2J)  = CBLK(LCELL,VORGBA2J)          &
!sorgam                             + CBLK(LCELL,VORGBA2I) * FM3
!sorgam     CBLK(LCELL,VORGBA3J)  = CBLK(LCELL,VORGBA3J)          &
!sorgam                             + CBLK(LCELL,VORGBA3I) * FM3
!sorgam     CBLK(LCELL,VORGBA4J)  = CBLK(LCELL,VORGBA4J)          &
!sorgam                             + CBLK(LCELL,VORGBA4I) * FM3

            CBLK(LCELL,VORGPAJ)   = CBLK(LCELL,VORGPAJ)           &
                                    + CBLK(LCELL,VORGPAI) * FM3

!unused     CBLK(LCELL,VP25AJ)    = CBLK(LCELL,VP25AJ)              &
!unused                             + CBLK(LCELL,VP25AI) * FM3

            CBLK(LCELL,VECJ)      = CBLK(LCELL,VECJ)              &
                                    + CBLK(LCELL,VECI) * FM3
            CBLK(LCELL,VSEASJ)    = CBLK(LCELL,VSEASJ)            &
                                    + CBLK(LCELL,VSEASI) * FM3

            ! Handle aerosol liquid water content the same way as the
            ! other mass tracers. This has to be done to ensure
            ! correct diagnostics of the aerosol liquid water content
            ! when writing to G3X-fields.

            CBLK(LCELL,VH2OAJ)    = CBLK(LCELL,VH2OAJ)           &
                                    + CBLK(LCELL,VH2OAI) * FM3

            ! *** update Aitken mode for mass loss to accumulation mode

            CBLK(LCELL,VSO4AI)    = CBLK(LCELL,VSO4AI) * PHM3
            CBLK(LCELL,VNH4AI)    = CBLK(LCELL,VNH4AI) * PHM3
            CBLK(LCELL,VNO3AI)    = CBLK(LCELL,VNO3AI) * PHM3

!sorgam     CBLK(LCELL,VORGARO1I) = CBLK(LCELL,VORGARO1I) * PHM3
!sorgam     CBLK(LCELL,VORGARO2I) = CBLK(LCELL,VORGARO2I) * PHM3
!sorgam     CBLK(LCELL,VORGALK1I) = CBLK(LCELL,VORGALK1I) * PHM3
!sorgam     CBLK(LCELL,VORGOLE1I) = CBLK(LCELL,VORGOLE1I) * PHM3
!sorgam     CBLK(LCELL,VORGBA1I)  = CBLK(LCELL,VORGBA1I) * PHM3
!sorgam     CBLK(LCELL,VORGBA2I)  = CBLK(LCELL,VORGBA2I) * PHM3
!sorgam     CBLK(LCELL,VORGBA3I)  = CBLK(LCELL,VORGBA3I) * PHM3
!sorgam     CBLK(LCELL,VORGBA4I)  = CBLK(LCELL,VORGBA4I) * PHM3

            CBLK(LCELL,VORGPAI)   = CBLK(LCELL,VORGPAI) * PHM3

!unused     CBLK(LCELL,VP25AI)    = CBLK(LCELL,VP25AI) * PHM3

            CBLK(LCELL,VECI)      = CBLK(LCELL,VECI) * PHM3
            CBLK(LCELL,VSEASI)    = CBLK(LCELL,VSEASI) * PHM3

            !al Handle aerosol liquid water content the same way as the
            !al other mass tracers. This has to be done since aerosol
            !al water content is calculated in ECHAM/MADE before
            !al SUBROUTINE MADE_AEROSTEP and not as in EURAD/MADE after
            !al SUBROUTINE MADE_AEROSTEP.

            CBLK(LCELL,VH2OAI)    = CBLK(LCELL,VH2OAI) * PHM3

         END IF ! end check on whether modal overlap is OK
      END IF ! end check on necessity for merging
   END DO ! loop for merging
   END IF

   ! set min value for all concentrations
   !ia actionia
   DO LCELL = 1, NUMCELLS
      DO SPC = 1, NSPCSDA
         CBLK(LCELL,SPC) = MAX(CONMIN, CBLK(LCELL,SPC))
      END DO
   END DO

END SUBROUTINE MADE_AEROSTEP

END MODULE messy_made
