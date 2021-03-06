!###############################################################################
!##  NOTE: Comment format `!>' is used for generation of code documentation   ##
!##        from source file via `doxygen'. Instructions see below.            ##
!##  Time-stamp: <2019-02-25 10:45:16 b309057 messy_made3.f90>                ##
!###############################################################################

! In order to generate the MADE3 (html and tex) documentation from the MADE3
! code files perform the following steps:
! 1. `cd' into your MESSy directory, subdirectory documentation/made3/
! 2. Temporarily remove one of the headers of the subroutine `update_xte' from
!    the file `messy/smil/messy_made3_si.f90'.
! 3. Run the command ``doxygen''
! 4. Do not forget to add the `update_xte' subroutine header back to
!    `messy/smil/messy_made3_si.f90' once you're done!

!> \image html made3.png "Schematic illustration of the MADE3 modes and aerosol composition."
!> \image latex made3.pdf "Schematic illustration of the MADE3 modes and aerosol composition." width=0.8\textwidth
!>
!> \mainpage Introduction
!>
!>   MADE3 (Modal Aerosol Dynamics model for Europe, adapted for global
!>   applications, 3<sup>rd</sup> generation) is an aerosol dynamics submodel
!>   for application within the MESSy framework. As a successor to MADE-in
!>   (\cite Aquila2011) and MADE as described by \cite Lauer2005, it draws
!>   heavily on the work by \cite Whitby1991 and \cite Binkowski1995.\n\n
!>
!>   The first generation of MADE was developed for application in a regional
!>   model by \cite Ackermann1998. Subsequently, MADE was adapted for global
!>   applications and implemented into the general circulation model ECHAM4 by
!>   \cite Lauer2005, and later transformed into a submodel (\cite Lauer2007)
!>   for the MESSy framework. The second generation submodel MADE-in was
!>   developed by \cite Aquila2011 as an extension to the MADE version used by
!>   \cite Lauer2007.
!>
!>   For the third generation submodel MADE3, the microphysical calculations
!>   (i.e., condensation/evaporation and coagulation) were extended to also take
!>   into account coarse particles, which were formerly regarded as passive. The
!>   gas-particle partitioning scheme was also extended, namely by inclusion of
!>   the hydrochloric acid/chloride equilibrium.
!>
!>   The aerosol size distribution is described in MADE3 by a superposition of
!>   nine lognormal modes in three different size ranges, namely the Aitken (ks,
!>   km, ki), accumulation (as, am, ai) and coarse modes (cs, cm, ci). In each
!>   size range, one mode represents fully soluble particles (ks, as, cs), one
!>   mode represents mixed particles (km, am, cm), and one mode represents
!>   insoluble particles (ki, ai, ci). Particles can consist of up to nine
!>   different components in MADE3:
!>   - sulfate (SO<sub>4</sub>),
!>   - ammonium (NH<sub>4</sub>),
!>   - nitrate (NO<sub>3</sub>),
!>   - sea-spray components other than chloride (Na),
!>   - chloride (Cl),
!>   - particulate organic matter (POM),
!>   - black carbon (BC),
!>   - mineral dust (DU), and
!>   - water (H<sub>2</sub>O).
!>
!>   MADE3 treats all aerosol-specific microphysical processes in the following
!>   sequence (using an operator-splitting approach):
!>   -# gas-particle partitioning via the thermodynamic equilibrium model EQSAM,
!>   -# condensation of sulfuric acid (H<sub>2</sub>SO<sub>4</sub>) and organic
!>      vapors (the latter being strongly simplified, see also \ref namelist),
!>   -# new particle formation
!>      (binary H<sub>2</sub>SO<sub>4</sub>/H<sub>2</sub>O nucleation), and
!>   -# coagulation.
!>
!>   For a more detailed model description, please refer to \cite Kaiser2014.
!>
!>   For the coupling to other submodels (i.e., optical properties, radiative
!>   effects, transport, cloud effects, cloud processing, scavenging,
!>   deposition), please refer to \ref interface.
!>
!>   The MADE3 code is organized in three code files and one namelist file:
!>   - the \link messy_made3 core module \endlink in
!>     \c messy/smcl/messy_made3.f90,
!>   - the \link messy_made3_si interface module \endlink in
!>     \c messy/smil/messy_made3_si.f90,
!>   - the \link messy_made3_box box model module \endlink in
!>     \c messy/mbm/made3/messy_made3_box.f90, and
!>   - the file \link namelist \c messy/nml/EXAMPLES/made3.nml \endlink
!>     containing the CTRL, CPL, and BOXINIT namelists to fine-tune the behavior
!>     of MADE3.

!> \brief MADE3 core module.

!> \authors Ingmar Ackermann et al., University Cologne, 1998
!>   - original MADE box model code
!> \authors Axel Lauer, DLR Oberpfaffenhofen, 2001-2003 (axel.lauer@dlr.de)
!>   - inclusion of MADE in ECHAM4
!>   - conversion of MADE to a MESSy submodel
!> \authors Valentina Aquila, DLR Oberpfaffenhofen, 2008/09
!>          (valentina.aquila@dlr.de)
!>   - development of MADE-in on the basis of MADE (and under the name MADE)
!> \authors Christopher Kaiser, DLR Oberpfaffenhofen, 2012-2016
!>          (christopher.kaiser@dlr.de)
!>   - adaptation of MADE-in to MESSy2
!>   - development of MADE3 on the basis of MADE-in

!> \version 3.0
!>   - new POM treatment:
!>     - all POM in insoluble modes considered hydrophobic (incl. SOA) (SMIL)
!>     - all POM in soluble modes considered hydrophilic (SMIL)
!>     - POM accumulation does not influence mixing state of insoluble particles
!>       (SMCL)
!>     - no photochemical POM "aging" (SMIL)
!>   - new treatment of nucleated particles:
!>     - set assumed initial dry diameter via CTRL namelist (\c #rset_nucsize)
!>     - adapt "nucleation" rate to assumed initial dry diameter
!>     - change default value from 3.5 nm "wet" diameter (now selectable by
!>       setting '\c rset_nucsize = F, 0.0') to 10 nm "dry" diameter
!> \version 2.2
!>   - changes only in interface: new emissions coupling
!> \version 2.1
!>   - namelist switches for individual subprocesses
!>   - some bug fixes
!> \version 2.0
!>   - small changes to obtain first working version of the MADE3 interface
!> \version 2.0b
!>   - inclusion of chloride as a separate aerosol species
!>   - inclusion of HCl/Cl chemistry in \c #made3_eqsam (most of it was already
!>     present, but commented out)
!>   - namelist switch for nucleation
!>   - Renaming criteria:
!>     - insoluble modes renamed to closest-sized mixed modes if 10% threshold
!>       of soluble mass is reached
!>     - Aitken mode particles renamed to corresponding accumulation modes
!>       either if 3<sup>rd</sup> moment growth rate of the Aitken mode larger
!>       than that of the accumulation mode, or if Aitken mode median diameter
!>       is greater than 30 nm \b and number concentration in the Aitken mode is
!>       greater than that in the accumulation mode
!> \version 1.0b
!>   - two new coarse modes
!>   - interactions of coarse modes with fine modes and gas phase
!>   - new subroutine for calculation of condensation factors (\c #cond_factors)
!>   - flux limit for gas-particle partitioning with coarse particles
!> \version 0.2b
!>   - separation of condensation and nucleation into two subroutines
!> \version 0.1b
!>   - adaptive time step
!>   - bug fixes in \c #made3_aerostep (solution to Riccati equation)
!>   - bug fixes in \c #made3_initialize_core (density averaging)
!>   - clean up and simplification of code
!>   - additional subroutines:
!>     - \c #renaming_fractions
!>     - \c #target_mode
!>     - \c #made3_rename

!> \todo
!>   - Check if MADE3_MODPAR should use all tracers in each mode to calculate
!>     3<sup>rd</sup> moments, total mass concentrations, wet diameters, and
!>     Knudsen
!>     numbers rather than only an arbitrary selection. Note that total mass
!>     concentrations are actually never used outside MADE3_MODPAR.
!>   - Include evaporation flux limit for semi-volatiles.
!>   - Include proper error treatment (i.e., find out what error numbers to
!>       use, how and where to define them, and how to document them; \i and
!>       \i implement \i checks!).
!>   - Convert FLUSSFRAC.. variables in subroutine MADE3_EQL3X to array.
!>   - Check FIXMEs.

MODULE messy_made3

  USE messy_main_constants_mem, ONLY: DP,             &
      pi, GRAV => g, AVO => N_A, RGASUNIV => R_gas,   &
      BOLTZ => k_B, MWAIR => M_air, MWH2O => M_H2O,   &
      MWCl => MCl, MWNa => MNa, RHOH2O => rho_h2o
  USE messy_main_tools,         ONLY: t_reset_par
! DEBUG+
!!$  USE messy_main_tools, ONLY: PTR_1D_ARRAY
! DEBUG-

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! GLOBAL PARAMETERS FOR USE WITH MESSy ===================================
  !> Name of submodel
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'made3'
  !> Version of submodel
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '3.0'

  ! AEROSOL MODES IN MADE3 =================================================
  ! If you change the mode indices, make sure to change also the definition of
  ! the variables sigma, dgini, nummin, and the settings in SCAV's CPL namelist
  !> Index of soluble Aitken mode (ks)
  INTEGER, PARAMETER, PUBLIC :: akn   = 1
  !> Index of mixed Aitken mode (km)
  INTEGER, PARAMETER, PUBLIC :: akns  = 2
  !> Index of insoluble Aitken mode (ki)
  INTEGER, PARAMETER, PUBLIC :: sooti = 3
  !> Index of soluble accumulation mode (as)
  INTEGER, PARAMETER, PUBLIC :: acc   = 4
  !> Index of mixed accumulation mode (am)
  INTEGER, PARAMETER, PUBLIC :: accs  = 5
  !> Index of insoluble accumulation mode (ai)
  INTEGER, PARAMETER, PUBLIC :: sootj = 6
  !> Index of soluble coarse mode (cs)
  INTEGER, PARAMETER, PUBLIC :: cor   = 7
  !> Index of mixed coarse mode (cm)
  INTEGER, PARAMETER, PUBLIC :: cors  = 8
  !> Index of insoluble coarse mode (ci)
  INTEGER, PARAMETER, PUBLIC :: sootc = 9
  !> Total number of modes
  INTEGER, PARAMETER, PUBLIC :: nmod  = 9

  ! AEROSOL SPECIES IN MADE3 ===============================================
  !> Index of sulfate
  INTEGER, PARAMETER, PUBLIC :: i_so4   = 1
  !> Index of ammonium
  INTEGER, PARAMETER, PUBLIC :: i_nh4   = 2
  !> Index of nitrate
  INTEGER, PARAMETER, PUBLIC :: i_no3   = 3
  !> Index of sea salt (cations and SO<sub>4</sub>)
  INTEGER, PARAMETER, PUBLIC :: i_ss    = 4
  !> Index of chloride
  INTEGER, PARAMETER, PUBLIC :: i_cl    = 5
  !> Index of organics
  INTEGER, PARAMETER, PUBLIC :: i_pom   = 6
  !> Index of BC
  INTEGER, PARAMETER, PUBLIC :: i_bc    = 7
  !> Index of tagged BC
  INTEGER, PARAMETER, PUBLIC :: i_bctag = 8
  !> Index of mineral dust
  INTEGER, PARAMETER, PUBLIC :: i_du    = 9
  !> Index of aerosol water
  INTEGER, PARAMETER, PUBLIC :: i_h2o   = 10
  !> Total number of aerosol species
  INTEGER, PARAMETER, PUBLIC :: nspec   = 10

  ! GAS PHASE SPECIES IN MADE3 =============================================
  !> Index of sulfuric acid
  INTEGER, PARAMETER, PUBLIC :: i_h2so4 = 1
  !> Index of ammonia
  INTEGER, PARAMETER, PUBLIC :: i_nh3   = 2
  !> Index of nitric acid
  INTEGER, PARAMETER, PUBLIC :: i_hno3  = 3
  !> Index of hydrochloric acid
  INTEGER, PARAMETER, PUBLIC :: i_hcl   = 4
  !> Index of SOA precursors
  INTEGER, PARAMETER, PUBLIC :: i_soa   = 5
  !> Total number of gas phase species
  INTEGER, PARAMETER, PUBLIC :: ngas  = 5

  ! CTRL-NAMELIST PARAMETERS (defaults in made3_read_nml) ==================
  !> Mode widths [-]
  REAL(dp), PUBLIC :: sigma(nmod)
  !> Accommodation coefficients [-]
  REAL(dp), PUBLIC :: alpha(nmod,ngas)
  !> Molecular diffusivities [m<sup>2</sup> s<sup>-1</sup>]
  REAL(dp), PUBLIC :: diff(ngas)
  ! op_cb_20171107+
  !> Maximum fraction of soluble matter in insoluble modes (aging criterion)
  REAL(dp), PUBLIC :: epsilon_max
  ! op_cb_20171107-
  !> Switch on/off and set assumed "dry" diam. of newly nucleated particles [nm]
  TYPE(t_reset_par), PUBLIC :: rset_nucsize
  !> Switch gas-particle partitioning of semivolatiles on/off
  LOGICAL, PUBLIC  :: l_eqsam
  !> Switch coagulation on/off
  LOGICAL, PUBLIC  :: l_coag
  !> Switch sulfuric acid condensation on/off
  LOGICAL, PUBLIC  :: l_cond
  !> Switch nucleation on/off
  LOGICAL,  PUBLIC :: l_nuc
  !> Switch renaming on/off
  LOGICAL, PUBLIC  :: l_rename

  ! GLOBAL PARAMETERS FOR USE IN MADE3 =====================================

  !--- 1.) Additional entry in 2nd dimension of array CBLK (see below):
  !> Index of CBLK "column" for gas phase concentrations
  INTEGER, PARAMETER, PUBLIC :: gas       = 10
  !> Number of CBLK "columns"
  INTEGER, PARAMETER, PUBLIC :: dim2_cblk = 10

  !--- 2.) Additional entries in 1st dimension of array CBLK (see below):
  !> Index of made3_main::CBLK "row" for aerosol number concentrations
  INTEGER, PARAMETER, PUBLIC :: i_num  = 11
  !> Index of CBLK "row" for aerosol 3<sup>rd</sup> moment concentrations
  INTEGER, PARAMETER, PUBLIC :: i_mom3 = 12
  !> Number of CBLK "rows"
  INTEGER, PARAMETER, PUBLIC :: dim1_cblk = 12

!obsolete  !--- 4.) Former indices corresponding to the compound masses and
!obsolete  !        mode numbers (in MADE3 internal array 'cblk'):
!sorgam   INTEGER, PARAMETER :: VORGARO1J =  7  ! accumulation mode
!sorgam                                             ! anthropogenic organic
!sorgam                                             ! aerosol from aromatics
!sorgam   INTEGER, PARAMETER :: VORGARO1I =  8  ! Aitken mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER :: VORGARO2J =  9  ! acc. mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER :: VORGARO2I = 10  ! Aitken mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER :: VORGALK1J = 11  ! acc. mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! alkanes and others except
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER :: VORGALK1I = 12  ! Aitken mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! alkanes and others except
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER :: VORGOLE1J = 13  ! acc. mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! alkanes and others except
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER :: VORGOLE1I = 14  ! Aitken mode anthropogenic
!sorgam                                             ! organic aerosol from
!sorgam                                             ! alkanes and others except
!sorgam                                             ! aromatics
!sorgam   INTEGER, PARAMETER :: VORGBA1J  = 15  ! accumulation mode biogenic
!sorgam                                             ! aerosol from a-pinene
!sorgam   INTEGER, PARAMETER :: VORGBA1I  = 16  ! Aitken mode biogenic
!sorgam                                             ! aerosol from a-pinene
!sorgam   INTEGER, PARAMETER :: VORGBA2J  = 17  ! accumulation mode biogenic
!sorgam                                             ! aerosol from a-pinene
!sorgam   INTEGER, PARAMETER :: VORGBA2I  = 18  ! Aitken mode biogenic
!sorgam                                             ! aerosol from a-pinene
!sorgam   INTEGER, PARAMETER :: VORGBA3J  = 19  ! accumulation mode biogenic
!sorgam                                             ! aerosol from limonene
!sorgam   INTEGER, PARAMETER :: VORGBA3I  = 20  ! Aitken mode biogenic
!sorgam                                             ! aerosol from limonene
!sorgam   INTEGER, PARAMETER :: VORGBA4J  = 21  ! accumulation mode biogenic
!sorgam                                             ! aerosol from limonene
!sorgam   INTEGER, PARAMETER :: VORGBA4I  = 22  ! Aitken mode biogenic
!sorgam                                             ! aerosol from limonene
!sorgam   INTEGER, PARAMETER :: VCVARO1   = 40  ! cond. vapor from aromatics
!sorgam   INTEGER, PARAMETER :: VCVARO2   = 41  ! cond. vapor from aromatics
!sorgam   INTEGER, PARAMETER :: VCVALK1   = 42  ! cond. vapor from anth.
!sorgam                                             ! alkanes
!sorgam   INTEGER, PARAMETER :: VCVOLE1   = 43  ! cond. vapor from anth.
!sorgam                                             ! olefines
!sorgam   INTEGER, PARAMETER :: VCVAPI1   = 44  ! cond. vapor from a-pinene
!sorgam   INTEGER, PARAMETER :: VCVAPI2   = 45  ! cond. vapor from a-pinene
!sorgam   INTEGER, PARAMETER :: VCVLIM1   = 46  ! cond. vapor from limonene
!sorgam   INTEGER, PARAMETER :: VCVLIM2   = 47  ! cond. vapor from limonene

  !--- 3.) Parameters for lognormal particle size distributions
  !        (assignments in made3_initialize_core):
  !> Initial median mode diameters [m]
  REAL(dp), PUBLIC           :: dgini(nmod)
  !> Minimum number concentrations [m<sup>-3</sup>]
  REAL(dp), PUBLIC           :: nummin(nmod)
  !> Minimum mass concentrations [ug m<sup>-3</sup>], derived from \c #dgini and \c #nummin
  REAL(dp), PUBLIC           :: MASSMIN(nmod,nspec+1)
  !> Index of MASSMIN "column" for total modal mass concentration
  INTEGER, PARAMETER, PUBLIC :: i_tot = nspec + 1

  !--- 4.) Constant for distribution of sea spray among species
  !> Mass fraction of Cl in sea spray (drawn from EQSAM code)
  REAL(dp), PARAMETER, PUBLIC :: mfCl = 0.5504_dp

  !--- 5.) Scalar variables for fixed standard deviations.
  ! \f$e^\frac{\left(\ln(\sigma_i)\right)^2}{8}\f$
  REAL(dp) :: e1(nmod)
  ! e1**4
  REAL(dp) :: es04(nmod)
  ! e1**5
  REAL(dp) :: es05(nmod)
  ! e1**8
  REAL(dp) :: es08(nmod)
  ! e1**9
  REAL(dp) :: es09(nmod)
  ! e1**16
  REAL(dp) :: es16(nmod)
  ! e1**20
  REAL(dp) :: es20(nmod)
  ! e1**25
  REAL(dp) :: es25(nmod)
  ! e1**32
  REAL(dp) :: es32(nmod)
  ! e1**36
  !> \f$e^{\frac{9}{2}\left(\ln(\sigma_i)\right)^2}\f$
  REAL(dp), PUBLIC :: es36(nmod)
  ! e1**49
  REAL(dp) :: es49(nmod)
  ! e1**64
  REAL(dp) :: es64(nmod)
  ! e1**100
  REAL(dp) :: es100(nmod)
  ! \f$\ln(\sigma_i)\f$
  REAL(dp) :: xxlsg(nmod)
  ! \f$\left(\ln(\sigma_i)\right)^2\f$
  REAL(dp) :: l2sigma(nmod)

  !--- 6.) Mathematical constants:
  REAL(dp), PARAMETER :: TWOPI   = 2.0_dp * pi
  REAL(dp), PARAMETER :: F6DPI   = 6.0_dp / pi
  !> \f$\frac{6\cdot10^{-9}}{\pi}\f$
  REAL(dp), PARAMETER, PUBLIC :: F6DPIM9 = 1.0e-9_dp * F6DPI
  REAL(dp), PARAMETER :: SQRT2   = 1.4142135623731_dp
  REAL(dp), PARAMETER :: ONE3    = 1.0_dp / 3.0_dp
  REAL(dp), PARAMETER :: TWO3    = 2.0_dp / 3.0_dp
  ! Numerical value for a minimum concentration
  REAL(dp), PARAMETER :: CONMIN = 1.0e-30_dp

  !--- 7.) Physical constants (if not assigned here,
  !        assignment in made3_initialize_core):
  !> Aerosol component densities [kg m<sup>-3</sup>]
  REAL(dp), TARGET, PUBLIC :: RHO(nspec)
  !> Molar masses [g mol<sup>-1</sup>]
  REAL(dp), TARGET, PUBLIC :: MW(nspec,2)
  ! Molar mass of elemental calcium
  REAL(dp), PARAMETER :: MWCa = 40.078_dp
  ! Molar mass of elemental potassium
  REAL(dp), PARAMETER :: MWK  = 39.0983_dp
  ! Molar mass of elemental magnesium
  REAL(dp), PARAMETER :: MWMg = 24.305_dp
  !> Index of \c #mw column for aerosol species
  INTEGER, PARAMETER, PUBLIC :: i_mwaero = 1
  !> Index of \c #mw column for gas species
  INTEGER, PARAMETER, PUBLIC :: i_mwgas  = 2
  !> Conversion factors from mass to 3<sup>rd</sup> moment concentration
  !> [ug m-3] -> [mom_3 m<sup>-3</sup>]
  REAL(dp), PUBLIC :: MASS2MOM3(nspec)
  ! *** Meteorological constants:
  ! Standard surface pressure [Pa]
  REAL(dp), PARAMETER :: P0  = 101325.0_dp
  ! Standard temperature [K]
  REAL(dp), PARAMETER :: T0  = 273.15_dp
  ! Standard surface temperature [K]
  REAL(dp), PARAMETER :: TS0 = 288.15_dp
  ! Minimum geometric mean diameter [m]
  REAL(dp), PARAMETER :: DGMIN = 1.0E-9_dp
  ! Maximum fraction of soluble matter in insoluble modes
!!$  REAL(dp), PARAMETER :: EPSILON_MAX = 0.10_dp  ! op_cb_20171107 nml param.
!!$  REAL(dp), PARAMETER :: EPSILON_MAX = 0.05_dp
  ! Species considered for ``aging'' criterion
  INTEGER, DIMENSION(5), PARAMETER :: SOLSPEC = (/i_so4,i_nh4,i_no3,i_cl,i_ss/)
! DEBUG+
!!$  TYPE(PTR_1D_ARRAY), DIMENSION(2,2), PUBLIC :: dbg_cblk_1
!!$  TYPE(PTR_1D_ARRAY), DIMENSION(2,2), PUBLIC :: dbg_cblk_2
!!$  REAL(dp), DIMENSION(:), POINTER, PUBLIC :: dbg_dry_1 => NULL()
!!$  REAL(dp), DIMENSION(:), POINTER, PUBLIC :: dbg_dry_2 => NULL()
! DEBUG-

  ! END GLOBAL PARAMETERS FOR USE IN MADE3 ==================================

  ! PUBLIC SUBROUTINES
  ! Reads CTRL namelist
  PUBLIC :: made3_read_nml
  ! Sets some constants
  PUBLIC :: made3_initialize_core
  ! Entry point to MADE3 microphysics
  PUBLIC :: made3_main


CONTAINS

!-------------------------------------------------------------------------------

  !> \brief Reads CTRL namelist
  !> \details Sets default values for CTRL namelist parameters, reads their
  !>   user-set values from \link ctrl \c made3.nml \endlink, and writes some of
  !>   those settings to stdout for reference.

  SUBROUTINE made3_read_nml(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    !> Error status flag
    INTEGER, INTENT(OUT) :: status
    !> Fortran "unit" for I/O
    INTEGER, INTENT(IN)  :: iou

    !> CTRL namelist
    ! op_cb_20171107: add epsilon_max
    NAMELIST /CTRL/ sigma, alpha, diff, epsilon_max, rset_nucsize, l_eqsam, &
         l_coag, l_cond, l_nuc, l_rename

    !--- Local variables:
    CHARACTER(LEN=*), PARAMETER :: substr = 'made3_read_nml'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES ...
    ! ... geometric standard deviations
    sigma(akn)   = 1.7_dp  ! Aitken mode, soluble - internally mixed
    sigma(akns)  = 1.7_dp  ! Aitken mode, soluble + BC - internally mixed
    sigma(sooti) = 1.7_dp  ! Aitken mode, BC - externally mixed
    sigma(acc)   = 2.0_dp  ! acc. mode, soluble - internally mixed
    sigma(accs)  = 2.0_dp  ! acc. mode, soluble + dust + BC - internally mixed
    sigma(sootj) = 2.0_dp  ! acc. mode, dust + BC - externally mixed
    sigma(cor)   = 2.2_dp  ! coarse mode, soluble - internally mixed
    sigma(cors)  = 2.2_dp  ! coarse mode, soluble + dust + BC - internally mixed
    sigma(sootc) = 2.2_dp  ! coarse mode, dust + BC - externally mixed
    ! ... accommodation coefficients
    alpha            = 0.1_dp
    alpha(:,i_h2so4) = 1.0_dp
    alpha(:,i_soa)   = 1.0_dp
    ! ... molecular diffusivities of gases
    diff          = 1.0e-5_dp    ! [m2 s-1]
    !bs  calculated from Reid, Prausnitz, and Poling, The properties of gases
    !bs  and liquids, 4th edition, McGraw-Hill, 1987, pp 587-588.
    !bs  Equation (11-4.4) was used.
    !bs  The value is at T = 273.16 K and P = 1.01325E05 Pa.
    !bs  T dependence is included via DIFFCORR (see code).
    !bs  updated from code of FSB on 23/03/99
    ! op_ck_20130403: It seems that 8.0e-6 was used prior to that update.
    diff(i_h2so4) = 9.362223e-06_dp
    diff(i_soa)   = 5.151174e-06_dp
    ! op_cb_20171107+
    ! ... maximum fraction of soluble matter in insoluble modes (aging
    ! criterion)
    epsilon_max = 0.10_dp
    ! op_cb_20171107-
    ! ... switch and value for assumed "dry" diam. of newly nucleated part. [nm]
    rset_nucsize%l = .TRUE.
    rset_nucsize%v = 10._dp
    ! ... switch for gas-particle partitioning of semivolatiles
    l_eqsam =  .TRUE.
    ! ... switch for coagulation
    l_coag  =  .TRUE.
    ! ... switch for H2SO4 condensation
    l_cond =   .TRUE.
    ! ... switch for nucleation
    l_nuc =    .TRUE.
    ! ... switch for renaming
    l_rename = .TRUE.

    ! read namelist:
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! diagnostic output
    WRITE(*,*) ''
    WRITE(*,*) ''
    WRITE(*,*) '--------------------------------------------------------------'
    WRITE(*,*) '--------------------------------------------------------------'
    WRITE(*,*) '---  Initialization of settings for aerosol module MADE3   ---'
    WRITE(*,*) '---                                                        ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for soluble Aitken mode       = ', sigma(akn), &
         '            ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for soluble acc. mode         = ', sigma(acc), &
         '            ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for soluble coarse mode       = ', sigma(cor), &
         '            ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for Aitken mode sol + insol   = ', sigma(akns), &
         '            ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for acc mode sol + insol      = ', sigma(accs), &
         '            ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for coarse mode sol + insol   = ', sigma(cors), &
         '            ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for Aitken mode insol         = ', sigma(sooti), &
         '            ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for acc. mode insol           = ', sigma(sootj), &
         '            ---'
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  sigma for coarse mode insol         = ', sigma(sootc), &
         '            ---'
    ! op_cb_20171107+
    WRITE(*,'(1X,A43,F4.2,A15)') &
         '---  aging criterion                     = ', EPSILON_MAX, &
         '               ---'
    ! op_cb_20171107-
    IF (rset_nucsize%l) THEN
       WRITE(*,'(1X,A43,F4.1,A15)') &
         '---  dry diam. of newly nuc. part. [nm]  = ', rset_nucsize%v, &
         '            ---'
    ENDIF

    WRITE(*,'(1X,A43,L1,A18)') &
         '---  g-p partitioning is calculated:       ', l_eqsam, &
         '               ---'
    WRITE(*,'(1X,A43,L1,A18)') &
         '---  coagulation is calculated:            ', l_coag, &
         '               ---'
    WRITE(*,'(1X,A43,L1,A18)') &
         '---  H2SO4 condensation is calculated:     ', l_cond, &
         '               ---'
    WRITE(*,'(1X,A43,L1,A18)') &
         '---  nucleation is calculated:             ', l_nuc, &
         '               ---'
    WRITE(*,'(1X,A43,L1,A18)') &
         '---  renaming is performed:                ', l_rename, &
         '               ---'
    WRITE(*,*) '--------------------------------------------------------------'
    WRITE(*,*) '--------------------------------------------------------------'
    WRITE(*,*) ''
    WRITE(*,*) ''

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE made3_read_nml

!-------------------------------------------------------------------------------

  !> \brief Sets some constants
  !> \details Sets \c #rho, \c #mass2mom3, \c #mw, \c #es36, \c #dgini, \c
  !>   #nummin, \c #massmin, which all remain unchanged during the simulation.

  SUBROUTINE made3_initialize_core

    IMPLICIT NONE
    INTRINSIC LOG, EXP

    REAL(dp) :: rho_av(nmod) ! average density of a 'minimum'-particle
    REAL(dp) :: mw_sum(nmod) ! sum of the molar masses of all species forming a
                             ! 'minimum'-particle
    INTEGER  :: jm           ! loop index

    ! --- code starts here -----------------------------------------------------

    ! *** aerosol component densities [kg m-3]:
    rho          = 0.0_dp
    rho(i_so4)   = 1.8e3_dp
    rho(i_nh4)   = 1.8e3_dp
    rho(i_no3)   = 1.8e3_dp
    rho(i_cl)    = 2.2e3_dp ! (source: PartMC-MOSAIC test case
                            !          1_urban_plume/aero_data.dat)
    rho(i_ss)    = 2.2e3_dp
    rho(i_pom)   = 1.0e3_dp
    rho(i_bc)    = 2.2e3_dp
    rho(i_bctag) = 2.2e3_dp
    rho(i_du)    = 2.5e3_dp
    rho(i_h2o)   = RHOH2O

    ! *** Factors for converting aerosol mass concentration [ug m-3]
    !     to 3rd moment concentration [mom_3 m-3]
    mass2mom3 = F6DPIM9 / rho

    ! *** molecular weights [g mol-1]
    mw                   = 0.0_dp
    mw(i_so4,i_mwaero)   = 96.0576_dp
    mw(i_nh4,i_mwaero)   = 18.03858_dp
    mw(i_no3,i_mwaero)   = 62.0649_dp
    mw(i_cl,i_mwaero)    = MWCl ! 35.45_dp
!!$    mw(i_ss,i_mwaero)    = 58.443_dp
    mw(i_ss,i_mwaero)    = MWNa ! 22.99_dp
    mw(i_pom,i_mwaero)   = 180.0_dp
    mw(i_bc,i_mwaero)    = 12.011_dp
    mw(i_bctag,i_mwaero) = 12.011_dp
    mw(i_du,i_mwaero)    = 40.08_dp
    mw(i_h2o,i_mwaero)   = MWH2O ! 18.02_dp
    mw(i_h2so4,i_mwgas)  = 98.07948_dp
    mw(i_nh3,i_mwgas)    = 17.03061_dp
    mw(i_hno3,i_mwgas)   = 63.01287_dp
    mw(i_hcl,i_mwgas)    = 36.4609_dp
    mw(i_soa,i_mwgas)    = 180.0_dp

    DO jm = 1, nmod
       xxlsg(jm)   = LOG(sigma(jm))
       l2sigma(jm) = xxlsg(jm) * xxlsg(jm)
       e1(jm)      = EXP(0.125_dp * l2sigma(jm))
       es04(jm)    = e1(jm) * e1(jm) * e1(jm) * e1(jm)
       es05(jm)    = es04(jm) * e1(jm)
       es08(jm)    = es04(jm) * es04(jm)
       es09(jm)    = es04(jm) * es05(jm)
       es16(jm)    = es08(jm) * es08(jm)
       es20(jm)    = es16(jm) * es04(jm)
       es25(jm)    = es16(jm) * es09(jm)
       es32(jm)    = es16(jm) * es16(jm)
       es36(jm)    = es16(jm) * es20(jm)
       es49(jm)    = es25(jm) * es20(jm) * es04(jm)
       es64(jm)    = es32(jm) * es32(jm)
       es100(jm)   = es36(jm) * es64(jm)
    END DO

    ! *** Initial geometric mean diameter [m]:
    dgini(akn)   = 0.01e-6_dp
    dgini(akns)  = 0.01e-6_dp
    dgini(sooti) = 0.01e-6_dp
    dgini(acc)   = 0.10e-6_dp
    dgini(accs)  = 0.10e-6_dp
    dgini(sootj) = 0.10e-6_dp
    dgini(cor)   = 1.00e-6_dp
    dgini(cors)  = 1.00e-6_dp
    dgini(sootc) = 1.00e-6_dp

    ! *** Minimum particle number concentrations [m-3]:
    !va old values
    !va nummin(akn)   = 1000.0_dp
    !va nummin(akns)  = 1000.0_dp
    !va nummin(sooti) = 1000.0_dp
    !va nummin(acc)   = 10.0_dp
    !va nummin(accs)  = 10.0_dp
    !va nummin(sootj) = 10.0_dp
    !va nummin(cor)   = 0.1_dp
    !va the modes with BC have a much lower number concentration with respect to
    !va the soluble modes, therefore I need lower minimum number concentrations
    nummin(akn)   = 100.0_dp
    nummin(akns)  = 1.0_dp
    nummin(sooti) = 1.0_dp
    nummin(acc)   = 1.0_dp
    nummin(accs)  = 0.1_dp
    nummin(sootj) = 0.1_dp
    nummin(cor)   = 0.1_dp
    nummin(cors)  = 0.1_dp
    nummin(sootc) = 0.1_dp

    ! Calculate minimum mass conc. [ug/m3] allowed for each species in each
    ! mode, assume the following molar fractions:
    !    soluble   Aitken mode (akn)  : 1 SO4 + 1 NO3 + 3 NH4
    !    soluble   acc.   mode (acc)  : 1 SO4 + 1 NO3 + 3 NH4
    !       ---> fully neutralized NH4(2)SO4 / NH4NO3 aerosol, no EC/OC
    !    soluble   coarse mode (cor)  : 1 SO4 + 1 NO3 + 3 NH4 + 1 SS
    !       ---> fully neutralized NH4(2)SO4 / NH4NO3 aerosol with sea salt
    !    mixed     Aitken mode (akns) : 1 SO4 + 1 NO3 + 3 NH4 + 1 BC
    !    mixed     acc. mode   (accs) : 1 SO4 + 1 NO3 + 3 NH4 + 1 BC
    !       ---> fully neutralized NH4(2)SO4 / NH4NO3 aerosol, no OC
    !    mixed     coarse mode (cors) : 1 SO4 + 1 NO3 + 3 NH4 + 1 BC + 1 SS
    !       ---> fully neutralized NH4(2)SO4 / NH4NO3 aerosol, no OC, with SS
    !    insoluble Aitken mode (sooti): 1 BC
    !    insoluble acc.   mode (sootj): 1 BC + 1 DU
    !    insoluble coarse mode (sootc): 1 BC + 1 DU

    ! ... set default values
    MASSMIN = 1.0e-30_dp   ! [ug m-3]

    ! ... calculate molar mass sums for all modes
    mw_sum(akn)   = MW(i_so4,i_mwaero) + MW(i_no3,i_mwaero) &
         + 3.0_dp * MW(i_nh4,i_mwaero)
    mw_sum(acc)   = mw_sum(akn)
    mw_sum(cor)   = mw_sum(akn) + MW(i_cl,i_mwaero) + MW(i_ss,i_mwaero)
    mw_sum(akns)  = MW(i_so4,i_mwaero) + MW(i_no3,i_mwaero) &
         + 3.0_dp * MW(i_nh4,i_mwaero) + MW(i_bc,i_mwaero)
    mw_sum(accs)  = mw_sum(akns)
    mw_sum(cors)  = mw_sum(akns) + MW(i_cl,i_mwaero) + MW(i_ss,i_mwaero)
    mw_sum(sooti) = MW(i_bc,i_mwaero)
    mw_sum(sootj) = MW(i_bc,i_mwaero) + MW(i_du,i_mwaero)
    mw_sum(sootc) = mw_sum(sootj)

    ! ... calculate average densities for all modes
    rho_av(akn)   = mw_sum(akn) &
         / (        MW(i_so4,i_mwaero) / rho(i_so4) &
         + 3.0_dp * MW(i_nh4,i_mwaero) / rho(i_nh4) &
         +          MW(i_no3,i_mwaero) / rho(i_no3))
    rho_av(acc)   = rho_av(akn)
    rho_av(cor)   = mw_sum(cor) &
         / (        MW(i_so4,i_mwaero) / rho(i_so4) &
         + 3.0_dp * MW(i_nh4,i_mwaero) / rho(i_nh4) &
         +          MW(i_no3,i_mwaero) / rho(i_no3) &
         +          MW(i_cl,i_mwaero)  / rho(i_cl)  &
         +          MW(i_ss,i_mwaero)  / rho(i_ss))
    rho_av(akns)  = mw_sum(akns) &
         / (        MW(i_so4,i_mwaero) / rho(i_so4) &
         + 3.0_dp * MW(i_nh4,i_mwaero) / rho(i_nh4) &
         +          MW(i_no3,i_mwaero) / rho(i_no3) &
         +          MW(i_bc,i_mwaero)  / rho(i_bc))
    rho_av(accs)  = rho_av(akns)
    rho_av(cors)  = mw_sum(cors) &
         / (        MW(i_so4,i_mwaero) / rho(i_so4) &
         + 3.0_dp * MW(i_nh4,i_mwaero) / rho(i_nh4) &
         +          MW(i_no3,i_mwaero) / rho(i_no3) &
         +          MW(i_bc,i_mwaero)  / rho(i_bc)  &
         +          MW(i_cl,i_mwaero)  / rho(i_cl)  &
         +          MW(i_ss,i_mwaero)  / rho(i_ss))
    rho_av(sooti) = rho(i_bc)
    rho_av(sootj) = mw_sum(sootj) &
         / (  MW(i_bc,i_mwaero) / rho(i_bc) &
         + MW(i_du,i_mwaero) / rho(i_du))
    rho_av(sootc) = rho_av(sootj)

    ! ... calculate total minimum masses for all modes
    DO jm = 1, nmod
       MASSMIN(jm,i_tot) = 1.0e9_dp * rho_av(jm) * pi * nummin(jm) &
            * dgini(jm) * dgini(jm) * dgini(jm) * es36(jm) / 6.0_dp
       ! ... and divide by the corresponding molar mass sums (will be reverted
       ! below)
       MASSMIN(jm,i_tot) = MASSMIN(jm,i_tot) / mw_sum(jm)
    END DO

    ! ... set soluble internally mixed modes
    MASSMIN(akn,i_so4)   =          MW(i_so4,i_mwaero) * MASSMIN(akn,i_tot)
    MASSMIN(akn,i_nh4)   = 3.0_dp * MW(i_nh4,i_mwaero) * MASSMIN(akn,i_tot)
    MASSMIN(akn,i_no3)   =          MW(i_no3,i_mwaero) * MASSMIN(akn,i_tot)
    MASSMIN(acc,i_so4)   =          MW(i_so4,i_mwaero) * MASSMIN(acc,i_tot)
    MASSMIN(acc,i_nh4)   = 3.0_dp * MW(i_nh4,i_mwaero) * MASSMIN(acc,i_tot)
    MASSMIN(acc,i_no3)   =          MW(i_no3,i_mwaero) * MASSMIN(acc,i_tot)
    MASSMIN(cor,i_so4)   =          MW(i_so4,i_mwaero) * MASSMIN(cor,i_tot)
    MASSMIN(cor,i_nh4)   = 3.0_dp * MW(i_nh4,i_mwaero) * MASSMIN(cor,i_tot)
    MASSMIN(cor,i_no3)   =          MW(i_no3,i_mwaero) * MASSMIN(cor,i_tot)
    MASSMIN(cor,i_cl)    =          MW(i_cl,i_mwaero)  * MASSMIN(cor,i_tot)
    MASSMIN(cor,i_ss)    =          MW(i_ss,i_mwaero)  * MASSMIN(cor,i_tot)

    ! ... set internally mixed modes (soluble + insoluble)
    MASSMIN(akns,i_so4)  =          MW(i_so4,i_mwaero) * MASSMIN(akns,i_tot)
    MASSMIN(akns,i_nh4)  = 3.0_dp * MW(i_nh4,i_mwaero) * MASSMIN(akns,i_tot)
    MASSMIN(akns,i_no3)  =          MW(i_no3,i_mwaero) * MASSMIN(akns,i_tot)
    MASSMIN(akns,i_bc)   =          MW(i_bc,i_mwaero)  * MASSMIN(akns,i_tot)
    MASSMIN(akns,i_bctag) = MASSMIN(akns,i_bc) ! op_mr_20181002
    MASSMIN(accs,i_so4)  =          MW(i_so4,i_mwaero) * MASSMIN(accs,i_tot)
    MASSMIN(accs,i_nh4)  = 3.0_dp * MW(i_nh4,i_mwaero) * MASSMIN(accs,i_tot)
    MASSMIN(accs,i_no3)  =          MW(i_no3,i_mwaero) * MASSMIN(accs,i_tot)
    MASSMIN(accs,i_bc)   =          MW(i_bc,i_mwaero)  * MASSMIN(accs,i_tot)
    MASSMIN(accs,i_bctag) = MASSMIN(accs,i_bc) ! op_mr_20181002
    MASSMIN(cors,i_so4)  =          MW(i_so4,i_mwaero) * MASSMIN(cors,i_tot)
    MASSMIN(cors,i_nh4)  = 3.0_dp * MW(i_nh4,i_mwaero) * MASSMIN(cors,i_tot)
    MASSMIN(cors,i_no3)  =          MW(i_no3,i_mwaero) * MASSMIN(cors,i_tot)
    MASSMIN(cors,i_bc)   =          MW(i_bc,i_mwaero)  * MASSMIN(cors,i_tot)
    MASSMIN(cors,i_bctag) = MASSMIN(cors,i_bc) ! op_mr_20181002
    MASSMIN(cors,i_cl)   =          MW(i_cl,i_mwaero)  * MASSMIN(cors,i_tot)
    MASSMIN(cors,i_ss)   =          MW(i_ss,i_mwaero)  * MASSMIN(cors,i_tot)

    ! ... set externally mixed insoluble modes
    MASSMIN(sooti,i_bc)    = MW(i_bc,i_mwaero) * MASSMIN(sooti,i_tot)
    MASSMIN(sooti,i_bctag) = MASSMIN(sooti,i_bc) ! op_mr_20181002
    MASSMIN(sootj,i_bc)    = MW(i_bc,i_mwaero) * MASSMIN(sootj,i_tot)
    MASSMIN(sootj,i_bctag) = MASSMIN(sootj,i_bc) ! op_mr_20181002
    MASSMIN(sootj,i_du)    = MW(i_du,i_mwaero) * MASSMIN(sootj,i_tot)
    MASSMIN(sootc,i_bc)    = MW(i_bc,i_mwaero) * MASSMIN(sootc,i_tot)
    MASSMIN(sootc,i_bctag) = MASSMIN(sootc,i_bc) ! op_mr_20181002
    MASSMIN(sootc,i_du)    = MW(i_du,i_mwaero) * MASSMIN(sootc,i_tot)

    ! Reset total minimum mass concentration per mode back to its actual value
    DO jm = 1, nmod
       MASSMIN(jm,i_tot) = 1.0e9_dp * rho_av(jm) * pi * nummin(jm) &
            * dgini(jm) * dgini(jm) * dgini(jm) * es36(jm) / 6.0_dp
    END DO

  END SUBROUTINE made3_initialize_core

!-------------------------------------------------------------------------------

  !> \brief Entry point to MADE3 microphysics
  !> \details Calls aerosol microphysics routine, scales calculated aerosol
  !>   tracer changes with cloud free fraction of grid box to account for
  !>   suspension of aerosol dynamics in clouds, and determines new aerosol
  !>   properties (wet and dry diameter, etc.)

  SUBROUTINE made3_main(status, BLKSIZE, NUMCELLS, PRESSURE, TEMPERATURE    &
       , RELHUM, PTMST, PSO4RAT, PSOA, CLOUDCOVER, CBLK, RH_HIST, DG, DGDRY &
       , PDENS                                                              &
!!$!op_va_20090225+
!!$       , BCin, BCsink                                                       &
!!$!op_va_20090225-
       )

    IMPLICIT NONE
    INTRINSIC MAX

    ! I/O
    !> Error status flag
    INTEGER,  INTENT(out)   :: status
    !> Size of input arrays
    INTEGER,  INTENT(in)    :: BLKSIZE
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)    :: NUMCELLS
! op_ck_20120411+
    !> Time step [s]
!    REAL(dp), INTENT(in)    :: PTMST
    REAL(dp), INTENT(inout) :: PTMST
! op_ck_20120411-
    !> Fractional cloud cover [-]
    REAL(dp), INTENT(in)    :: CLOUDCOVER(BLKSIZE)
    !> H<sub>2</sub>SO<sub>4</sub>(g) rate of change [ug m<sup>-3</sup>
    !> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: PSO4RAT(BLKSIZE)
    !> Air pressure [Pa]
    REAL(dp), INTENT(in)    :: PRESSURE(BLKSIZE)
    !> SOA precursor emissions [ug m<sup>-3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: PSOA(BLKSIZE)
    !> Air temperature [K]
    REAL(dp), INTENT(in)    :: TEMPERATURE(BLKSIZE)
    !> Relative humidity [-]
    REAL(dp), INTENT(in)    :: RELHUM(BLKSIZE)
    !> Deliquescence history  (1.0 = all particles dry, 2.0 = all particles wet)
    REAL(dp), INTENT(inout) :: RH_HIST(nmod,BLKSIZE)
    !> Tracer array
    REAL(dp), INTENT(inout) :: CBLK(dim1_cblk,dim2_cblk,BLKSIZE)
    !> Wet median mode diameters [m]
    REAL(dp), INTENT(out)   :: DG(nmod,BLKSIZE)
    !> Dry median mode diameters [m]
    REAL(dp), INTENT(out)   :: DGDRY(nmod,BLKSIZE)
    !> Average modal wet aerosol densities [kg m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: PDENS(nmod,BLKSIZE)
!!$!op_va_20090225+
!!$    REAL(dp) :: BCin(BLKSIZE)             ! BURDEN of ext. mixed BC. saved for
!!$                                          ! estimate of BC lifetime
!!$    REAL(dp) :: BCsink(BLKSIZE)           ! SINK of ext. mixed BC. saved for
!!$                                          ! estimate of BC lifetime
!!$!op_va_20090225-

    ! LOCAL
    REAL(dp) :: SO4RAT_IN(NUMCELLS)       ! H2SO4(g) production rate [ug/m3/s]
!sorgam   REAL(dp) :: drog_in(NUMCELLS,LDROG)      ! accumulated production of
!sorgam                                           ! anthropogenic AND biogenic
!sorgam                                           ! organic aerosol precursor
!sorgam                                           ! [ug m-3]
    REAL(dp) :: CLOUDFREE                 ! cloud free fraction of grid box
    REAL(dp) :: CBLK_M1(dim1_cblk,dim2_cblk,NUMCELLS) ! copy of 'CBLK' array
                                                      ! before MADE3; used to
                                                      ! update tracer tendencies
    REAL(dp) :: XLM(NUMCELLS)             ! atmospheric mean free path [m]
    REAL(dp) :: AMU(NUMCELLS)             ! atmospheric dynamic viscosity
                                          ! [kg m-1 s-1]
    REAL(dp) :: PMASS(nmod,NUMCELLS)      ! modal mass concentrations [ug m-3]
    REAL(dp) :: KN(nmod,NUMCELLS)         ! modal Knudsen numbers
    REAL(dp) :: MERGENUM(NUMCELLS)        ! mode merging, number [m-3]
                                          ! (for diagnostics only)
    REAL(dp) :: MERGEM3(NUMCELLS)         ! mode merging, 3rd moment [mom_3 m-3]
                                          ! (for diagnostics only)
    REAL(dp) :: MERGENUMs(NUMCELLS)
    REAL(dp) :: MERGEM3s(NUMCELLS)
    REAL(dp) :: MERGENUMsoot(NUMCELLS)
    REAL(dp) :: MERGEM3soot(NUMCELLS)
    REAL(dp) :: mom3dry(nmod)             ! 3rd moments, dry (aux. variables)
    INTEGER  :: NCELL                     ! loop index (gridcells)
    INTEGER  :: jm                        ! loop index (modes)

    ! --- code starts here -----------------------------------------------------

! DEBUG+
!!$    dbg_cblk_1(1,1)%ptr(:) = CBLK(i_cl,acc,:)
!!$    dbg_cblk_1(1,2)%ptr(:) = CBLK(i_cl,cor,:)
!!$    dbg_cblk_1(2,1)%ptr(:) = CBLK(i_num,acc,:)
!!$    dbg_cblk_1(2,2)%ptr(:) = CBLK(i_num,cor,:)
! DEBUG-
    ! Initializations
    status    = 9       ! unspecified error
    DG(:,:)        = -1.0_dp
    DGDRY(:,:)     = -1.0_dp
    PDENS(:,:)     = -1.0_dp
    XLM(:)         = -1.0_dp
    AMU(:)         = -1.0_dp
    PMASS(:,:)     = -1.0_dp
    KN(:,:)        = -1.0_dp
    mom3dry(:)     = -1.0_dp

    ! Validate input parameters and make a copy of current tracer
    ! concentrations (gas+aerosol) for calculating the updated
    ! tracer concentrations after aerosol dynamics calculations.
    DO NCELL = 1, NUMCELLS      ! loop over gridcells
       CBLK_M1(:,:,NCELL) = CBLK(:,:,NCELL)
       ! copy H2SO4(g) production rate [ug m-3 s-1]
       SO4RAT_IN(NCELL) = PSO4RAT(NCELL)
    END DO    ! loop over gridcells (NCELL)

    ! *** CALL AEROSOL DYNAMICS DRIVER ***
    CALL MADE3_AEROPROC(status, NUMCELLS &
         , CBLK(1:dim1_cblk,1:dim2_cblk,1:NUMCELLS), PTMST            &
         , TEMPERATURE(1:NUMCELLS), PRESSURE(1:NUMCELLS), RELHUM(1:NUMCELLS) &
         , SO4RAT_IN, PSOA(1:NUMCELLS)                                       &
!sorgam                ORGARO1RAT, ORGARO2RAT,                   &
!sorgam                ORGALK1RAT, ORGOLE1RAT,                   &
!sorgam                ORGBIO1RAT, ORGBIO2RAT,                   &
!sorgam                ORGBIO3RAT, ORGBIO4RAT,                   &
!sorgam                DROG, LDROG, NCV, NACV,                   &
         , XLM, AMU, DG(1:nmod,1:NUMCELLS), PMASS, PDENS(1:nmod,1:NUMCELLS)  &
         , KN, RH_HIST(1:nmod,1:NUMCELLS)                                    &
         , MERGENUM, MERGEM3, MERGENUMs, MERGEM3s, MERGENUMSOOT, MERGEM3SOOT &
!!$!op_va_20090225+
!!$         , BCin, BCsink                                                      &
!!$!op_va_20090225-
         )

    ! Aerosol dynamics are only calculated for cloud free conditions.
    ! Within a cloud covered area, all aerosol dynamics are suspended.
    ! Thus, the tendencies calculated by SUBROUTINE MADE3_AEROPROC are
    ! scaled with (1 - cloud cover) to get a new mean value for the
    ! whole grid box.
    ! cloudy area of grid box = cloud cover
    ! cloud free area         = 1 - cloud cover
    DO NCELL = 1, NUMCELLS
       ! CBLK(new) - CBLK_M1 = (CBLK(current)-CBLK_M1) * (1-cloud cover)
       CLOUDFREE = 1.0_dp - CLOUDCOVER(NCELL)
       CBLK(:,:,NCELL) = (CBLK(:,:,NCELL) - CBLK_M1(:,:,NCELL)) &
            * CLOUDFREE + CBLK_M1(:,:,NCELL)
       CBLK(:,:,NCELL) = MAX(CONMIN, CBLK(:,:,NCELL))
    END DO

    ! *** Get new distribution information:
    CALL MADE3_MODPAR( NUMCELLS, CBLK(1:dim1_cblk,1:dim2_cblk,1:NUMCELLS) &
         , TEMPERATURE(1:NUMCELLS), PRESSURE(1:NUMCELLS), PMASS           &
         , PDENS(1:nmod,1:NUMCELLS), XLM, AMU, DG(1:nmod,1:NUMCELLS), KN)

! DEBUG+
!!$    dbg_dry_1(:) = dgdry(akn,:)
! DEBUG-
    ! Calculate modal diameters (dry) from 3rd moment.
    DO NCELL = 1, NUMCELLS
       DO jm = 1, nmod
          mom3dry(jm)     = MAX(conmin, cblk(i_mom3,jm,ncell) &
               - MASS2MOM3(i_h2o) * cblk(i_h2o,jm,ncell))
! op_ck_20121106  FIXME: Why should we keep larger particles and not reset their
!                        size as well in case mom3dry is ``too low''?
! op_ck_20140612         -> see below
          dgdry(jm,ncell) = MAX(dgmin,  (mom3dry(jm) &
               / (cblk(i_num,jm,ncell) * es36(jm)))**one3)
       END DO
    END DO
! DEBUG+
!!$    dbg_dry_2(:) = dgdry(akn,:)
! DEBUG-

    status = 0  ! no error

  END SUBROUTINE made3_main

!-------------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! PRIVATE ROUTINES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  !> \brief Calls the routines for the individual microphysical processes
  !> \details Determines aerosol properties (wet diameter, Knudsen number etc.),
  !>   then calls the routines for
  !>     -# gas-particle partitioning of NH<sub>3</sub>/NH<sub>4</sub>,
  !>        HNO<sub>3</sub>/NO<sub>3</sub>, HCl/Cl, and H<sub>2</sub>O,
  !>     -# coagulation rates (after update of aerosol properties),
  !>     -# condensation rates,
  !>     -# nucleation rate.
  !>     .
  !>   Uses the calculated rates to compute new number and mass concentrations
  !>   and calls renaming routine. Adaptive time step is used if necessary.

  SUBROUTINE MADE3_AEROPROC (status, NUMCELLS, CBLK, DT           &
       , BLKTA, BLKPRS, BLKRH, SO4RAT, SOA_MADE3                           &
       , XLM, AMU, DG, PMASS, PDENS, KN, RH_HIST                           &
       , MERGENUM, MERGEM3, MERGENUMs, MERGEM3s, MERGENUMSOOT, MERGEM3SOOT &
!sorgam                    ORGARO1RAT, ORGARO2RAT,                   &
!sorgam                    ORGALK1RAT, ORGOLE1RAT,                   &
!sorgam                    ORGBIO1RAT, ORGBIO2RAT,                   &
!sorgam                    ORGBIO3RAT, ORGBIO4RAT,                   &
!sorgam                    DROG, LDROG, NCV, NACV,                   &
!!$!op_va_20090225+
!!$       , BCin, BCsink                                                      &
!!$!op_va_20090225-
     )

    IMPLICIT NONE
    INTRINSIC ABS, REAL, MAX

    ! I/O
    !> Error status flag
    INTEGER,  INTENT(out)   :: status
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)    :: NUMCELLS
    !> Tracer array
    REAL(dp), INTENT(inout) :: CBLK(dim1_cblk,dim2_cblk,NUMCELLS)
! op_ck_20120411+
    !> Time step [s]
!    REAL(dp), INTENT(in)    :: DT
    REAL(dp), INTENT(inout) :: DT
! op_ck_20120411-
    !> Air temperature [K]
    REAL(dp), INTENT(in)    :: BLKTA(NUMCELLS)
    !> Air pressure [Pa]
    REAL(dp), INTENT(in)    :: BLKPRS(NUMCELLS)
    !> Relative humidity [-]
    REAL(dp), INTENT(in)    :: BLKRH(NUMCELLS)
    !> Deliquescence history
    REAL(dp), INTENT(inout) :: RH_HIST(nmod,NUMCELLS)
    !> H<sub>2</sub>SO<sub>4</sub>(g) rate of change [ug m<sup>-3</sup>
    !> s<sup>-1</sup>]
    REAL(dp), INTENT(inout) :: SO4RAT(NUMCELLS)
    !> SOA precursor emissions [ug m<sup>-3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: SOA_MADE3(NUMCELLS)
    !> Atmospheric mean free path [m]
    REAL(dp), INTENT(out)   :: XLM(NUMCELLS)
    !> Atmospheric dynamic viscosity [kg m<sup>-1</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(out)   :: AMU(NUMCELLS)
    !> Wet median mode diameters [m]
    REAL(dp), INTENT(out)   :: DG(nmod,NUMCELLS)
    !> Modal mass concentrations [ug m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: PMASS(nmod,NUMCELLS)
    !> Average modal wet aerosol densities [kg m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: PDENS(nmod,NUMCELLS)
    !> Modal Knudsen numbers
    REAL(dp), INTENT(out)   :: KN(nmod,NUMCELLS)
    !> Number concentration of particles renamed from akn to acc
    !> [m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: MERGENUM(NUMCELLS)
    !> 3<sup>rd</sup> moment concentration of particles renamed from akn to acc
    !> [mom_3 m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: MERGEM3(NUMCELLS)
    !> Number conc. of particles renamed from akns to accs [m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: MERGENUMs(NUMCELLS)
    !> 3<sup>rd</sup> moment conc. of particles renamed from akns to accs [mom_3
    !> m-3]
    REAL(dp), INTENT(out)   :: MERGEM3s(NUMCELLS)
    !> Number conc. of particles renamed from sooti to sootj [m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: MERGENUMSOOT(NUMCELLS)
    !> 3<sup>rd</sup> moment conc. of particles renamed from sooti to sootj
    !> [mom_3 m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: MERGEM3SOOT(NUMCELLS)
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     aromatics
!sorgam   REAL(dp) :: ORGARO1RAT(NUMCELLS)
!sorgam
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     aromatics
!sorgam   REAL(dp) :: ORGARO2RAT(NUMCELLS)
!sorgam
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     alkanes & others
!sorgam   REAL(dp) :: ORGALK1RAT(NUMCELLS)
!sorgam
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     alkenes & others
!sorgam   REAL(dp) :: ORGOLE1RAT(NUMCELLS)
!sorgam
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO1RAT(NUMCELLS)
!sorgam
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO2RAT(NUMCELLS)
!sorgam
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO3RAT(NUMCELLS)
!sorgam
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO4RAT(NUMCELLS)
!sorgam   !bs * organic condensable vapor production rate
!sorgam   REAL(dp) :: DROG(NUMCELLS,LDROG)    ! Delta ROG conc. [ppm]
!sorgam   INTEGER, INTENT(in)    :: LDROG           ! # of organic aerosol
!sorgam                                             ! precursor
!sorgam   INTEGER, INTENT(in)    :: NCV             ! total # of cond. vapors
!sorgam                                             ! & SOA species
!sorgam   INTEGER, INTENT(in)    :: NACV            ! # of anthrop. cond.
!sorgam                                             ! vapors & SOA species
!!$!op_va_20090225+
!!$    REAL(dp) :: BCin(NUMCELLS)   ! BURDEN of ext. mixed BC. saved for estimate
!!$                                ! of BC lifetime
!!$    REAL(dp) :: BCsink(NUMCELLS) ! SINK of ext. mixed BC. saved for estimate of
!!$                                ! BC lifetime
!!$!op_va_20090225-

    ! LOCAL

    ! modal condensation factors (see comments in MADE3_CONDENSE):
    REAL(dp) :: FCONC(nmod,NUMCELLS)     ! SO4
    REAL(dp) :: FCONC_ORG(nmod,NUMCELLS) ! SOA

    ! rates for secondary particle formation:
    REAL(dp) :: DMDT(NUMCELLS) ! rate of production of new so4 mass by particle
                               ! formation [ug m-3 s-1]
    REAL(dp) :: DNDT(NUMCELLS) ! rate of producton of new particle number by
                               ! particle form. [m-3 s-1]

    ! modal growth rates for 3rd moment by condensation of precursor
    ! vapor on existing particles [mom_3 m-3 s-1]
    REAL(dp) :: CGR3(nmod,NUMCELLS)

    ! number (0th moment) coagulation rates [m3 s-1]
    REAL(dp) :: CR0(nmod,nmod,NUMCELLS)

    ! 3rd moment intermodal transfer rates
    REAL(dp) :: C30(nmod,nmod,NUMCELLS) ! 3rd moment transfer rate by
                                       ! intermodal coagulation

    ! Variables for condensation calculations
    REAL(dp) :: DELTASO4A(NUMCELLS)     ! increment of conc. added to sulfate
                                       ! aerosol by condensation [ug m-3]
    REAL(dp) :: h2so4vap(NUMCELLS)      ! Vapor conc. prior to addition from
                                       ! chemistry and emissions, but after
                                       ! condensation [ug m-3]

! op_ck_20120309+
    ! local variables for adaptive time step (see comment in MADE3_AEROSTEP)
    REAL(dp) :: CBLK_old(dim1_cblk,dim2_cblk,NUMCELLS)
    REAL(dp) :: TMST, TMST_old
    INTEGER  :: i_loop
! op_ck_20120309-

    INTEGER  :: LCELL


    ! -------------------------------- Start code ------------------------------

    ! Array initializations
    DG        = -1.0_dp
    PMASS     = -1.0_dp
    PDENS     = -1.0_dp
    KN        = -1.0_dp
    FCONC     =  0.0_dp
    FCONC_ORG =  0.0_dp
    DELTASO4A =  0.0_dp
    CGR3      =  0.0_dp
    CR0       =  0.0_dp
    C30       =  0.0_dp
    DMDT      =  0.0_dp
    DNDT      =  0.0_dp

! op_ck_20121126+
!!$! op_ck_20120309+
!!$
!!$   ! For adaptive time step (see comment in MADE3_AEROSTEP)
!!$
!!$   TMST     = DT
!!$   CBLK_old = CBLK
!!$   i_loop   = 1
!!$
!!$   DO
!!$
!!$! op_ck_20120309-
! op_ck_20121126-

    IF (l_eqsam) THEN
       ! *** Get size distribution information:
       CALL MADE3_MODPAR(NUMCELLS, CBLK, BLKTA, BLKPRS, PMASS, PDENS &
            , XLM, AMU, DG, KN)
       ! *** Get water, ammonium, nitrate, and chloride content:
       CALL MADE3_EQL3X(NUMCELLS, CBLK, BLKTA, BLKRH, BLKPRS, DT, DG &
            , RH_HIST)
    END IF
! DEBUG+
!!$    dbg_cblk_2(1,1)%ptr(1:NUMCELLS) = CBLK(i_cl,acc,:)
!!$    dbg_cblk_2(1,2)%ptr(1:NUMCELLS) = CBLK(i_cl,cor,:)
!!$    dbg_cblk_2(2,1)%ptr(1:NUMCELLS) = CBLK(i_num,acc,:)
!!$    dbg_cblk_2(2,2)%ptr(1:NUMCELLS) = CBLK(i_num,cor,:)
! DEBUG-

! op_ck_20121126+
    ! For adaptive time step (see comment in MADE3_AEROSTEP)
    TMST     = DT
    CBLK_old = CBLK
    i_loop   = 1

    DO
! op_ck_20121126-

       ! *** Get size distribution information:
       CALL MADE3_MODPAR(NUMCELLS, CBLK, BLKTA, BLKPRS, PMASS, PDENS &
            , XLM, AMU, DG, KN)

       IF (l_coag) THEN
          ! *** Calculate coagulation rates:
          CALL MADE3_COAGRATE(NUMCELLS, CBLK, BLKTA, PDENS, AMU, DG, KN &
               , CR0, C30)
       END IF

       ! *** Get condensation and particle formation (nucleation) rates:
       IF (l_cond) THEN
          CALL MADE3_CONDENSE(NUMCELLS, CBLK, TMST, BLKTA, BLKPRS &
               , SO4RAT, SOA_MADE3                                         &
!sorgam                ORGARO1RAT, ORGARO2RAT,      &
!sorgam                ORGALK1RAT, ORGOLE1RAT,      &
!sorgam                ORGBIO1RAT, ORGBIO2RAT,      &
!sorgam                ORGBIO3RAT, ORGBIO4RAT,      &
!sorgam                DROG, LDROG, NCV, NACV,      &
               , DG, FCONC, FCONC_ORG, DELTASO4A, CGR3, h2so4vap)
       ELSE
          h2so4vap(:) = CBLK(i_h2so4,gas,:)
       END IF

       IF (l_nuc) THEN
          ! *** Do Vehkam?ki et al., 2002 (J. Geophys. Res.) nucleation
          CALL MADE3_NUCLEATE(BLKTA, BLKRH, h2so4vap, SO4RAT &
               , NUMCELLS, DNDT, DMDT, CGR3)
       END IF

!sorgam   !bs * Secondary organic aerosol module (SORGAM)
!sorgam
!sorgam   CALL SORGAM( BLKTA, BLKPRS,
!sorgam                ORGARO1RAT, ORGARO2RAT,
!sorgam                ORGALK1RAT, ORGOLE1RAT,
!sorgam                ORGBIO1RAT, ORGBIO2RAT,
!sorgam                ORGBIO3RAT, ORGBIO4RAT,
!sorgam                DROG, LDROG, NCV, NACV,
!sorgam                CBLK, NUMCELLS,
!sorgam                DT )
!sorgam
!sorgam   ! *  Secondary organic aerosol module (SORGAM)

! op_ck_20120309+
       TMST_old = TMST
! op_ck_20120309-

       ! *** Advance forward in time TMST seconds:
       CALL MADE3_AEROSTEP(NUMCELLS, CBLK, TMST &
!sorgam                    ORGARO1RAT, ORGARO2RAT,                  &
!sorgam                    ORGALK1RAT, ORGOLE1RAT,                  &
!sorgam                    ORGBIO1RAT, ORGBIO2RAT,                  &
!sorgam                    ORGBIO3RAT, ORGBIO4RAT,                  &
!unused                    EPM25I, EPM25J, EORGI,EORGJ, EECI, EECJ, &
!unused                    ESOIL, ESEAS, EPMCOARSE,                 &
            , FCONC, FCONC_ORG                           &
!unused                    PMASSN, PMASSA, PMASSC,                  &
            , DMDT, DNDT, DELTASO4A, SOA_MADE3, CR0, C30, CGR3)

! op_ck_20120309+
       dt_change: IF (ABS(TMST - TMST_old) / TMST_old .GT. 1.e-30_dp) THEN
          CBLK = CBLK_old
          i_loop = 1
       ELSE
          end_loop: IF (ABS(REAL(i_loop, kind=dp) * TMST - DT) / DT &
                        .LE. 1.e-30_dp ) THEN
! op_ck_20120411+
             DT = TMST
! op_ck_20120411-
             EXIT
          ELSE
             i_loop = i_loop + 1
          END IF end_loop
       END IF dt_change

    END DO ! adaptive time step loop
! op_ck_20120309-

    IF (l_rename) THEN
       ! *** Avoid merging of modes
       CALL made3_rename(status, NUMCELLS, CBLK, DG, CGR3 &
            , MERGENUM, MERGEM3, MERGENUMs, MERGEM3s, MERGENUMsoot, MERGEM3soot)
    ELSE
       MERGENUM     = 0.0_dp
       MERGEM3      = 0.0_dp
       MERGENUMs    = 0.0_dp
       MERGEM3s     = 0.0_dp
       MERGENUMsoot = 0.0_dp
       MERGEM3soot  = 0.0_dp
    END IF

    ! set min value for all concentrations
    DO LCELL = 1, NUMCELLS
       CBLK(:,:,LCELL) = MAX(CONMIN, CBLK(:,:,LCELL))
    END DO

  END SUBROUTINE MADE3_AEROPROC

!-------------------------------------------------------------------------------

  !> \brief Assigns output from gas-particle partitioning routine to aerosol and
  !>   gas phase tracers
  !> \details Aggregates composition of modes in same size range (i.e., all
  !>   Aitken modes, all accumulation modes, all coarse modes), calls
  !>   gas-particle partitioning routine for each size range, and distributes
  !>   results among the modes (and the gas phase). Limits gas-to-particle
  !>   transfer of HNO<sub>3</sub>, NH<sub>3</sub>, and HCl to coarse particles,
  !>   assuming that the near-particle surface gas concentration is 0 at all
  !>   times.

  SUBROUTINE MADE3_EQL3X(NUMCELLS, CBLK, BLKTA, BLKRH, press, dt &
       , diam, RH_HIST)

    IMPLICIT NONE
    INTRINSIC SUM, MAX

    ! *** INPUT ***
    !> Number of grid cells in input arrays
    INTEGER, INTENT(in)     :: NUMCELLS
    !> Tracer array
    REAL(dp), INTENT(inout) :: CBLK(dim1_cblk,dim2_cblk,NUMCELLS)
    !> Air temperature [K]
    REAL(dp), INTENT(in)    :: BLKTA(NUMCELLS)
    !> Relative humidity [-]
    REAL(dp), INTENT(in)    :: BLKRH(NUMCELLS)
    !> Air pressure [Pa]
    REAL(dp), INTENT(in)    :: press(NUMCELLS)
    !> Time step [s]
    REAL(dp), INTENT(in)    :: dt
    !> Wet median mode diameters [m]
    REAL(dp), INTENT(in)    :: diam(nmod,NUMCELLS)
    !> Deliquescence history
    REAL(dp), INTENT(inout) :: RH_HIST(nmod,NUMCELLS)

    ! *** LOCAL ***
    INTEGER :: LCELL                 ! loop counter
    REAL(dp):: FLUSSFRACi, FLUSSFRACis, FLUSSFRACSOOTi !liq. frac. in each mode
    REAL(dp):: FLUSSFRACj, FLUSSFRACjs, FLUSSFRACSOOTj !liq. frac. in each mode
    REAL(dp):: FLUSSFRACc, FLUSSFRACcs, FLUSSFRACSOOTc !liq. frac. in each mode

   ! input arrays for EQSAM
    REAL(dp) :: inSO4(NUMCELLS)       ! SO4 (w/o seasalt associated SO4) [ug/m3]
    REAL(dp) :: inNH4(NUMCELLS)       ! NH4 [ug/m3]
    REAL(dp) :: inNO3(NUMCELLS)       ! aerosol NO3 [ug/m3]
    REAL(dp) :: inCl(NUMCELLS)        ! Cl [ug/m3]
    REAL(dp) :: inSS(NUMCELLS)        ! sea salt [ug/m3]
!    REAL(dp) :: inH2SO4(NUMCELLS)     ! H2SO4 (gas phase) [ug/m3]
    REAL(dp) :: inNH3(NUMCELLS)       ! NH3 (gas phase) [ug/m3]
    REAL(dp) :: inHNO3(NUMCELLS)      ! HNO3 (gas phase) [ug/m3]
    REAL(dp) :: inHCl(NUMCELLS)       ! HCl (gas phase) [ug/m3]
    REAL(dp) :: inTOT(NUMCELLS)

    ! output from EQSAM
    REAL(dp) :: PNH4(NUMCELLS)        ! NH4 [ug/m3]
    REAL(dp) :: PNO3(NUMCELLS)        ! NO3 [ug/m3]
    REAL(dp) :: PCl(NUMCELLS)         ! Cl  [ug/m3]
    REAL(dp) :: WH2O(NUMCELLS)        ! aerosol water [ug/m3]
    REAL(dp) :: GNH3(NUMCELLS)        ! NH3 (gas phase) [ug/m3]
    REAL(dp) :: GNO3(NUMCELLS)        ! HNO3 (gas phase) [ug/m3]
    REAL(dp) :: GHCl(NUMCELLS)        ! HCl (gase phase) [ug/m3]

    ! flux limit calculations
    REAL(dp) :: cf_tot(NUMCELLS)      ! total condensation factor            [-]
    REAL(dp) :: cf_rel(nmod,NUMCELLS) ! relative cond. factors for modes     [-]
    REAL(dp) :: nh3min(NUMCELLS)      ! min. NH3 conc. due to cond. lim.[ug m-3]
    REAL(dp) :: hno3min(NUMCELLS)     ! min. HNO3 conc. (cond. lim.)    [ug m-3]
    REAL(dp) :: hclmin(NUMCELLS)      ! min. HCl conc. (cond. lim.)     [ug m-3]


    ! -------------------------------- Start code ------------------------------

    ! *** AITKEN MODEs ***

    ! aerosol SO4
    inSO4 = sum(CBLK(i_so4,(/akn,akns,sooti/),1:NUMCELLS), dim=1)
    ! aerosol NO3
    inNO3 = sum(CBLK(i_no3,(/akn,akns,sooti/),1:NUMCELLS), dim=1)
    ! aerosol NH4
    inNH4 = sum(CBLK(i_nh4,(/akn,akns,sooti/),1:NUMCELLS), dim=1)
    ! sea salt
    inCl  = sum(CBLK(i_cl,(/akn,akns,sooti/),1:NUMCELLS), dim=1)
    inSS  = sum(CBLK(i_ss,(/akn,akns,sooti/),1:NUMCELLS), dim=1)

    inTOT = inSO4 + inNO3 + inNH4 + inCl + inSS

    ! gasphase concentrations
!    inH2SO4 = CBLK(i_h2so4,gas,1:NUMCELLS)
    inNH3  = CBLK(i_nh3,gas,1:NUMCELLS)
    inHNO3 = CBLK(i_hno3,gas,1:NUMCELLS)
    inHCl  = CBLK(i_hcl,gas,1:NUMCELLS)

    CALL MADE3_EQSAM(BLKTA, BLKRH, RH_HIST(akn,:), inSO4, inNH4, inNH3, inNO3 &
         , inHNO3, inCl, inSS, inHCl, PNH4, PNO3, PCl, WH2O, GNH3, GNO3, GHCl &
         , numcells)

    RH_HIST(akns,1:NUMCELLS) = RH_HIST(akn,1:NUMCELLS)
    RH_HIST(sooti,1:NUMCELLS) = RH_HIST(akn,1:NUMCELLS)

    ! update modal concs.

    DO LCELL=1,NUMCELLS

       if (inTOT(LCELL).GT.1e-30_dp) then
          FLUSSFRACi  = sum(CBLK(SOLSPEC,akn,LCELL))  / inTOT(LCELL)
          FLUSSFRACis = sum(CBLK(SOLSPEC,akns,LCELL)) / inTOT(LCELL)
          FLUSSFRACSOOTi = 1.0_dp - FLUSSFRACi - FLUSSFRACis
       else
          FLUSSFRACi = 0.0_dp
          FLUSSFRACis = 0.0_dp
          FLUSSFRACsooti = 0.0_dp
       endif

       CBLK(i_h2o,akn,LCELL)   = MAX(1.0e-30_dp,WH2O(LCELL)    * FLUSSFRACi)
       CBLK(i_h2o,akns,LCELL)  = MAX(1.0e-30_dp,WH2O(LCELL)   * FLUSSFRACis)
       CBLK(i_h2o,sooti,LCELL) = MAX(1.0e-30_dp,WH2O(LCELL)* FLUSSFRACSOOTi)

       CBLK(i_nh4,akn,LCELL)   = MAX(1.0e-30_dp,PNH4(LCELL)    * FLUSSFRACi)
       CBLK(i_nh4,akns,LCELL)  = MAX(1.0e-30_dp,PNH4(LCELL)   * FLUSSFRACis)
       CBLK(i_nh4,sooti,LCELL) = MAX(1.0e-30_dp,PNH4(LCELL)* FLUSSFRACSOOTi)

       CBLK(i_no3,akn,LCELL)   = MAX(1.0e-30_dp,PNO3(LCELL)    * FLUSSFRACi)
       CBLK(i_no3,akns,LCELL)  = MAX(1.0e-30_dp,PNO3(LCELL)   * FLUSSFRACis)
       CBLK(i_no3,sooti,LCELL) = MAX(1.0e-30_dp,PNO3(LCELL)* FLUSSFRACSOOTi)

       CBLK(i_cl,akn,LCELL)    = MAX(1.0e-30_dp,PCl(LCELL)     * FLUSSFRACi)
       CBLK(i_cl,akns,LCELL)   = MAX(1.0e-30_dp,PCl(LCELL)    * FLUSSFRACis)
       CBLK(i_cl,sooti,LCELL)  = MAX(1.0e-30_dp,PCl(LCELL) * FLUSSFRACSOOTi)

       CBLK(i_nh3,gas,LCELL)   = GNH3(LCELL)
       CBLK(i_hno3,gas,LCELL)  = GNO3(LCELL)
       CBLK(i_hcl,gas,LCELL)   = GHCl(LCELL)

    END DO

    ! *** ACCUMULATION MODE ***
    
    ! aerosol SO4
    inSO4 = sum(CBLK(i_so4,(/acc,accs,sootj/),1:NUMCELLS), dim=1)
    ! aerosol NO3
    inNO3 = sum(CBLK(i_no3,(/acc,accs,sootj/),1:NUMCELLS), dim=1)
    ! aerosol NH4
    inNH4 = sum(CBLK(i_nh4,(/acc,accs,sootj/),1:NUMCELLS), dim=1)
    ! sea salt
    inCl  = sum(CBLK(i_cl,(/acc,accs,sootj/),1:NUMCELLS), dim=1)
    inSS  = sum(CBLK(i_ss,(/acc,accs,sootj/),1:NUMCELLS), dim=1)

    inTOT = inSO4 + inNO3 + inNH4 + inCl + inSS

    ! gasphase concentrations
!    inH2SO4 = CBLK(i_h2so4,gas,1:NUMCELLS)
    inNH3  = CBLK(i_nh3,gas,1:NUMCELLS)
    inHNO3 = CBLK(i_hno3,gas,1:NUMCELLS)
    inHCl  = CBLK(i_hcl,gas,1:NUMCELLS)

    CALL MADE3_EQSAM(BLKTA, BLKRH, RH_HIST(acc,:), inSO4, inNH4, inNH3, inNO3 &
         , inHNO3, inCl, inSS, inHCl, PNH4, PNO3, PCl, WH2O, GNH3, GNO3, GHCl &
         , numcells )

    RH_HIST(accs,1:NUMCELLS) = RH_HIST(acc,1:NUMCELLS)
    RH_HIST(sootj,1:NUMCELLS) = RH_HIST(acc,1:NUMCELLS)

    ! update modal concs.

    DO LCELL=1,NUMCELLS

       if (inTOT(LCELL).GT.1e-30_dp) then
          FLUSSFRACj  = sum(CBLK(SOLSPEC,acc,LCELL))  / inTOT(LCELL)
          FLUSSFRACjs = sum(CBLK(SOLSPEC,accs,LCELL)) / inTOT(LCELL)
          FLUSSFRACSOOTj = 1.0_dp - FLUSSFRACj - FLUSSFRACjs
       else
          FLUSSFRACj = 0.0_dp
          FLUSSFRACjs = 0.0_dp
          FLUSSFRACsootj = 0.0_dp
       endif

       CBLK(i_h2o,acc,LCELL)   = MAX(1.0e-30_dp,WH2O(LCELL)    * FLUSSFRACj)
       CBLK(i_h2o,accs,LCELL)  = MAX(1.0e-30_dp,WH2O(LCELL)   * FLUSSFRACjs)
       CBLK(i_h2o,sootj,LCELL) = MAX(1.0e-30_dp,WH2O(LCELL)* FLUSSFRACSOOTj)

       CBLK(i_nh4,acc,LCELL)   = MAX(1.0e-30_dp,PNH4(LCELL)    * FLUSSFRACj)
       CBLK(i_nh4,accs,LCELL)  = MAX(1.0e-30_dp,PNH4(LCELL)   * FLUSSFRACjs)
       CBLK(i_nh4,sootj,LCELL) = MAX(1.0e-30_dp,PNH4(LCELL)* FLUSSFRACSOOTj)

       CBLK(i_no3,acc,LCELL)   = MAX(1.0e-30_dp,PNO3(LCELL)    * FLUSSFRACj)
       CBLK(i_no3,accs,LCELL)  = MAX(1.0e-30_dp,PNO3(LCELL)   * FLUSSFRACjs)
       CBLK(i_no3,sootj,LCELL) = MAX(1.0e-30_dp,PNO3(LCELL)* FLUSSFRACSOOTj)

       CBLK(i_cl,acc,LCELL)    = MAX(1.0e-30_dp,PCl(LCELL)     * FLUSSFRACj)
       CBLK(i_cl,accs,LCELL)   = MAX(1.0e-30_dp,PCl(LCELL)    * FLUSSFRACjs)
       CBLK(i_cl,sootj,LCELL)  = MAX(1.0e-30_dp,PCl(LCELL) * FLUSSFRACSOOTj)

       CBLK(i_nh3,gas,LCELL)   = GNH3(LCELL)
       CBLK(i_hno3,gas,LCELL)  = GNO3(LCELL)
       CBLK(i_hcl,gas,LCELL)   = GHCl(LCELL)
    END DO

    ! *** COARSE MODE ***
    
    ! aerosol SO4
    inSO4 = sum(CBLK(i_so4,(/cor,cors,sootc/),1:NUMCELLS), dim=1)
    ! aerosol NO3
    inNO3 = sum(CBLK(i_no3,(/cor,cors,sootc/),1:NUMCELLS), dim=1)
    ! aerosol NH4
    inNH4 = sum(CBLK(i_nh4,(/cor,cors,sootc/),1:NUMCELLS), dim=1)
    ! sea salt
    inCl  = sum(CBLK(i_cl,(/cor,cors,sootc/),1:NUMCELLS), dim=1)
    inSS  = sum(CBLK(i_ss,(/cor,cors,sootc/),1:NUMCELLS), dim=1)

    inTOT = inSO4 + inNO3 + inNH4 + inCl + inSS

    ! gasphase concentrations
!    inH2SO4 = 0.0_dp
    inNH3  = CBLK(i_nh3,gas,1:NUMCELLS)
    inHNO3 = CBLK(i_hno3,gas,1:NUMCELLS)
    inHCl  = CBLK(i_hcl,gas,1:NUMCELLS)
    ! calculate concentration minima from maximum condensation flux
    ! Assumptions:
    ! * gas phase production neglected, because it is considered in the
    !   chemistry scheme before the MADE call (-> operator splitting)
    ! * loss or gain in gas phase during one time step is small compared to gas
    !   phase concentration
    ! NH3
    CALL cond_factors(i_nh3, NUMCELLS, press, BLKTA &
         , CBLK(i_num,1:nmod,1:NUMCELLS), diam, cf_tot, cf_rel)
    nh3min = inNH3 &
         * (1.0_dp - SUM(cf_rel((/cor,cors,sootc/),:), dim=1) * cf_tot * dt)
    ! HNO3
    CALL cond_factors(i_hno3, NUMCELLS, press, BLKTA &
         , CBLK(i_num,1:nmod,1:NUMCELLS), diam, cf_tot, cf_rel)
    hno3min = inHNO3 &
         * (1.0_dp - SUM(cf_rel((/cor,cors,sootc/),:), dim=1) * cf_tot * dt)
    ! HCl
    CALL cond_factors(i_hcl, NUMCELLS, press, BLKTA &
         , CBLK(i_num,1:nmod,1:NUMCELLS), diam, cf_tot, cf_rel)
    hclmin = inHCl &
         * (1.0_dp - SUM(cf_rel((/cor,cors,sootc/),:), dim=1) * cf_tot * dt)

    CALL MADE3_EQSAM(BLKTA, BLKRH, RH_HIST(cor,:), inSO4, inNH4, inNH3, inNO3 &
         , inHNO3, inCl, inSS, inHCl, PNH4, PNO3, PCl, WH2O, GNH3, GNO3, GHCl &
         , numcells, minNH3=nh3min, minHNO3=hno3min, minHCl=hclmin)

    RH_HIST(cors,1:NUMCELLS) = RH_HIST(cor,1:NUMCELLS)
    RH_HIST(sootc,1:NUMCELLS) = RH_HIST(cor,1:NUMCELLS)

    ! update modal concs.

    DO LCELL=1,NUMCELLS

       if (inTOT(LCELL).GT.1e-30_dp) then
          FLUSSFRACc  = sum(CBLK(SOLSPEC,cor,LCELL))  / inTOT(LCELL)
          FLUSSFRACcs = sum(CBLK(SOLSPEC,cors,LCELL)) / inTOT(LCELL)
          FLUSSFRACSOOTc = 1.0_dp - FLUSSFRACc - FLUSSFRACcs
       else
          FLUSSFRACc = 0.0_dp
          FLUSSFRACcs = 0.0_dp
          FLUSSFRACsootc = 0.0_dp
       endif

       CBLK(i_h2o,cor,LCELL)   = MAX(1.0e-30_dp,WH2O(LCELL)    * FLUSSFRACc)
       CBLK(i_h2o,cors,LCELL)  = MAX(1.0e-30_dp,WH2O(LCELL)   * FLUSSFRACcs)
       CBLK(i_h2o,sootc,LCELL) = MAX(1.0e-30_dp,WH2O(LCELL)* FLUSSFRACSOOTc)

       CBLK(i_nh4,cor,LCELL)   = MAX(1.0e-30_dp,PNH4(LCELL)    * FLUSSFRACc)
       CBLK(i_nh4,cors,LCELL)  = MAX(1.0e-30_dp,PNH4(LCELL)   * FLUSSFRACcs)
       CBLK(i_nh4,sootc,LCELL) = MAX(1.0e-30_dp,PNH4(LCELL)* FLUSSFRACSOOTc)

       CBLK(i_no3,cor,LCELL)   = MAX(1.0e-30_dp,PNO3(LCELL)    * FLUSSFRACc)
       CBLK(i_no3,cors,LCELL)  = MAX(1.0e-30_dp,PNO3(LCELL)   * FLUSSFRACcs)
       CBLK(i_no3,sootc,LCELL) = MAX(1.0e-30_dp,PNO3(LCELL)* FLUSSFRACSOOTc)

       CBLK(i_cl,cor,LCELL)    = MAX(1.0e-30_dp,PCl(LCELL)     * FLUSSFRACc)
       CBLK(i_cl,cors,LCELL)   = MAX(1.0e-30_dp,PCl(LCELL)    * FLUSSFRACcs)
       CBLK(i_cl,sootc,LCELL)  = MAX(1.0e-30_dp,PCl(LCELL) * FLUSSFRACSOOTc)

       CBLK(i_nh3,gas,LCELL)   = GNH3(LCELL)
       CBLK(i_hno3,gas,LCELL)  = GNO3(LCELL)
       CBLK(i_hcl,gas,LCELL)   = GHCl(LCELL)

    END DO

  END SUBROUTINE MADE3_EQL3X

!-------------------------------------------------------------------------------

  !> \brief Calculates gas-particle partitioning of
  !>   NH<sub>3</sub>/NH<sub>4</sub>, HNO<sub>3</sub>/NO<sub>3</sub>, HCl/Cl,
  !>   and H<sub>2</sub>O
  !> \authors Swen Metzger, MPI-CH, 3/11/99 (metzger@mpch-mainz.mpg.de)
  !>   - modified 2002, 2003
  !> \authors Axel Lauer, DLR, 2004 (axel.lauer@dlr.de)
  !>   - modified for use with MADE
  !> \authors Christopher Kaiser, DLR, 2013 (christopher.kaiser@dlr.de)
  !>   - included limitation of gas-to-particle transfer for coarse particles
  !>   - re-enabled HCl/Cl treatment
  !> \details
  !>   <h4>Purpose</h4>
  !>     EQSAM (EQuilibrium Simplified Aerosol Model) is a new and
  !>     computationally efficient thermodynamic aerosol composition
  !>     model that allows to calculate the gas/aerosol equilibrium
  !>     partitioning, including aerosol water, sufficiently fast and
  !>     accurate for global (or even regional) modeling. EQSAM is based
  !>     on a number of parameterizations, including single solute
  !>     molalities and activity coefficients (AC). The thermodynamic
  !>     framework (domains and subdomains, internally mixed aerosols)
  !>     is the same as of more sophisticated thermodynamic equilibrium
  !>     models (EQMs), e.g., of ISORROPIA (Nenes et al., 1998). Details
  !>     are given in the references below (and the references therein).\n
  !>     The main assumption on which EQSAM/EQMs are based is
  !>     thermodynamical and chemical equilibrium. From this assumption
  !>     it directly follows that the aerosol water activity (aw) equals
  !>     the ambient relative humidity (RH), if the water vapor pressure
  !>     is sufficiently larger than the partial vapor pressure of the
  !>     aerosol compounds. This is approximately true for tropospheric
  !>     aerosols. Given the large amount of water vapor present, water
  !>     vapor and aerosol water equilibrate relatively faster compared
  !>     to all other aerosol compounds. This is subsequently also true
  !>     for single aerosol compounds. The water activity of single
  !>     solutes must also equal RH under this assumption. Therefore, the
  !>     so called ZSR-relation is (and can be) used to calculate the
  !>     aerosol associated water mass (simply from the sum of all water
  !>     mass fractions that are derived from measured single solute
  !>     molalities).\n
  !>     In contrast to other EQMs, EQSAM utilizes the fact that the RH
  !>     fixes the water activity (under the above assumptions) and the
  !>     consequence that any changes in RH also causes changes in the
  !>     aerosol water mass and, hence, aerosol activity (including
  !>     activity coefficients). Thus, an decrease (increase) in RH
  !>     decrease (increases) the aerosol water mass (and water activity).
  !>     This can change the aerosol composition, e.g., due to
  !>     condensation (evaporation/crystallization), because the vapor
  !>     pressure above the aerosol reduces (increases). In turn, a vapor
  !>     pressure reduction (increase) due to changes in the aerosol
  !>     composition is compensated by an associated condensation
  !>     (evaporation) of water vapor to maintain the aerosol molality to
  !>     remain constant (because aw=RH). Furthermore, the aerosol water
  !>     mainly depends on the aerosol mass and the type of solute, so
  !>     that parameterizations of single solute molalities and activity
  !>     coefficients can be defined, only depending on the type of
  !>     solute and RH. The advantage of using such parameterizations is
  !>     that the entire aerosol equilibrium composition can be solved
  !>     analytically, i.e., non-iteratively, which considerably reduces
  !>     the amount of CPU time that is usually need for aerosol
  !>     thermodynamic calculations (especially if an EQM is incorporated
  !>     in an aerosol dynamical model that is in turn embedded in a high
  !>     resolution regional or global model).\n
  !>     However, EQSAM should still be regarded as a starting point for
  !>     further developments. There is still room for improvements. For
  !>     instance, this code is not yet numerically optimized (vectorized)
  !>     and a number of improvements with respect to an explicit
  !>     treatment of additional equilibrium reactions, missing (or only
  !>     implicit) dissociation, and a basic parameterization of the
  !>     water uptake.\n
  !>     Note that EQSAM was originally developed to calculate the
  !>     gas/aerosol equilibrium partitioning of the
  !>     ammonium-sulfate-nitrate-water system for climate models,
  !>     excluding solid compounds. This version (eqsam_v03d.f90) is
  !>     extended with respect to sea salt. Solids/hysteresis are treated
  !>     in a simplified manner. Results of a box model comparison with
  !>     ISORROPIA will be available from the web page. Please also note
  !>     that the water uptake is based on additional (unpublished)
  !>     parameterizations for single solute molalities, which are
  !>     derived from tabulated measurements used in ISORROPIA. Note
  !>     further that this extended version (eqsam_v03d.f90) is not yet
  !>     published. A publication is in progress.
  !>
  !>   <h4>Method</h4>
  !>     equilibrium / internal mixture assumption / aw=rh\n
  !>     System:
  !>     NH<sub>3</sub>,
  !>     NH<sub>4</sub><sup>+</sup>/H<sub>2</sub>SO<sub>4</sub><sup>+</sup>,
  !>     HSO<sub>4</sub><sup>-</sup>,
  !>     SO<sub>4</sub><sup>--</sup>/HNO<sub>3</sub>,
  !>     NO<sub>3</sub><sup>-</sup>, HCl,Cl<sup>-</sup>/Na<sup>+</sup>,
  !>     H<sub>2</sub>O
  !>             (K<sup>+</sup>,Ca<sup>++</sup>,Mg<sup>++</sup>)
  !>
  !>   <h4>References</h4>
  !>     - Metzger, S., PhD Thesis, University Utrecht, 2000\n
  !>       http://www.library.uu.nl/digiarchief/dip/diss/1930853/inhoud.htm
  !>     - \cite Metzger2002 \n
  !>       http://dx.doi.org/10.1029/2001JD001102
  !>     - \cite Metzger2002a \n
  !>       http://dx.doi.org/10.1029/2001JD001103
  !>     - \cite Metzger2006 \n
  !>       http://www.atmos-chem-phys.net/6/2549/2006/
  !>
  !> \copyright 1999-2003\n
  !>   Department of Atmospheric Chemistry, Max-Planck-Institute for
  !>   Chemistry\n
  !>   http://www.mpch-mainz.mpg.de/~metzger
  !>
  !> \version eqsam_v03d.f90 (MPI-CH, June 2003)
  !>   - Gama parameterizations now according to Metzger 2002 (JGR Appendix)
  !>   - Improved pH calculations (still restricted to strong acids)
  !>   - Removed bug that lead to too high nitrate formation in dry and cold
  !>     cold regions (UT/LS)
  !>   - Removed bug in solid/hysteresis calculations
  !>   .
  !>   (both bugs introduced in eqsam_v03b.f90 by cleaning up eqsam_v02a.f90)
  !> \version eqsam_v03c.f90 (MPI-CH, April 2003)
  !>   - More accurate paramterizations of single solute molalities (Na, Cl)
  !>   - Cleanded up RHD subdomain structure
  !>   - Improved water uptake (Na, Cl)
  !> \version eqsam_v03b.f90 (MPI-CH, March 2003)
  !>   - System extended to HCl, Cl-/Na+
  !>   - Parameterization (fit) of additional HNO<sub>3</sub> uptake
  !>     removed. Instead, complete analytical solution of equilibrium
  !>     reactions, based on the AC-RH relationship.
  !> \version eqsam_v03.f90 (IMAU, October 1999)
  !>   - Test version (included in TM3)
  !> \version eqsam_v02a.f90 (IMAU, April 2000)
  !>   - Box model version
  !> \version eqsam_v02.f90 (IMAU, October 1999)
  !>   - TM3 version
  !>   - Includes solids and additional HNO<sub>3</sub> uptake on acidic
  !>     aerosols (parameterized)
  !> \version eqsam_v01b.f90 (MPI-CH, January 2003)
  !>   - Same as eqsam_v01a.f90 (additional lines though uncommented for test
  !>     purposes only)
  !> \version eqsam_v01a.f90 (IMAU, April 2000)
  !>   - Box model version
  !> \version eqsam_v01.f90 (IMAU, October 1999)
  !>   - TM3 version.
  !>   - First and most basic version (without solids) for better vectorization
  !>     (for global modeling)
  !>   - System: NH<sub>3</sub>,
  !>   - NH<sub>4</sub><sup>+</sup>/H<sub>2</sub>SO<sub>4</sub><sup>+</sup>,
  !>     HSO<sub>4</sub><sup>-</sup>,
  !>     SO<sub>4</sub><sup>--</sup>/HNO<sub>3</sub>,
  !>     NO<sub>3</sub><sup>-</sup>, H<sub>2</sub>O
  !>   - Based on equilibrium / internal mixture assumption / aw=rh /
  !>     ZSR-relation parameterization of activcity coefficients (AC), i.e., an
  !>     AC-RH relationship
  !>
  !> \todo
  !>   - Split ion-pairs into ions for water parameterizations (since info is
  !>     actually available)
  !>   - Include uptake/dissociation of NH<sub>3</sub>, HNO<sub>3</sub>, HCl
  !>     (mainly to get pH right at near neutral conditions)
  !>   - Extension to K<sup>+</sup>, Ca<sup>++</sup>, Mg<sup>++</sup>,
  !>     CO<sub>2</sub>/(CO<sub>3</sub>)2<sup>--</sup>/HCO<sub>3</sub><sup>-</sup>,
  !>     SOA etc. (maybe not)
  !>   - Vectorization, translation of hardcoded formulas in array syntax
  !>   - I/O Interface and program structure clean up
  !>   - EQSAM info webpage

  SUBROUTINE MADE3_EQSAM(temp, rh, rh_hist                     &
       , inSO4, inNH4, inNH3, inNO3, inHNO3, inCl, inSS, inHCl &
       , PNH4, PNO3, PCl, WH2O, GNH3, GNO3, GHCl               &
       , numcells, minNH3, minHNO3, minHCl)

    IMPLICIT NONE
    INTRINSIC ABS, EXP, LOG, MAX, MIN, NINT, REAL, SQRT, PRESENT

    ! input for EQSAM

    !> Number of grid cells in input arrays
    integer, intent(in)     :: numcells
    !> Air temperature [K]
    real(dp), intent(in)    :: temp(numcells)
    !> Relative humidity [-]
    real(dp), intent(in)    :: rh(numcells)
    !> Deliquescence history
    real(dp), intent(inout) :: rh_hist(numcells)
    !> SO<sub>4</sub>(p) (w/o sea salt SO<sub>4</sub>) input [ug m<sup>-3</sup>]
    real(dp), intent(in)    :: inSO4(numcells)
    !> NH<sub>4</sub>(p) input [ug m<sup>-3</sup>]
    real(dp), intent(in)    :: inNH4(numcells)
    !> NO<sub>3</sub>(p) input [ug m<sup>-3</sup>]
    real(dp), intent(in)    :: inNO3(numcells)
    !> Cl(p) input [ug m<sup>-3</sup>]
    real(dp), intent(in)    :: inCl(numcells)
    !> sea salt (cations & SO<sub>4</sub>) input [ug m<sup>-3</sup>]
    real(dp), intent(in)    :: inSS(numcells)
    !> NH<sub>3</sub>(g) input [ug m<sup>-3</sup>]
    real(dp), intent(in)    :: inNH3(numcells)
    !> HNO<sub>3</sub>(g) input [ug m<sup>-3</sup>]
    real(dp), intent(in)    :: inHNO3(numcells)
    !> HCl(g) input [ug m<sup>-3</sup>]
    real(dp), intent(in)    :: inHCl(numcells)

    ! output from EQSAM

    !> NH<sub>4</sub>(p) output [ug m<sup>-3</sup>]
    real(dp), intent(out)   :: PNH4(numcells)
    !> NO<sub>3</sub>(p) output [ug m<sup>-3</sup>]
    real(dp), intent(out)   :: PNO3(numcells)
    !> Cl(p) output [ug m<sup>-3</sup>]
    real(dp), intent(out)   :: PCl(numcells)
    !> H<sub>2</sub>O output [ug m<sup>-3</sup>]
    real(dp), intent(out)   :: WH2O(numcells)
    !> NH<sub>3</sub>(g) output [ug m<sup>-3</sup>]
    real(dp), intent(out)   :: GNH3(numcells)
    !> HNO<sub>3</sub>(g) output [ug m<sup>-3</sup>]
    real(dp), intent(out)   :: GNO3(numcells)
    !> HCl(g) output [ug m<sup>-3</sup>]
    real(dp), intent(out)   :: GHCl(numcells)

! op_ck_20130408+
    !> Minimum NH<sub>3</sub>(g) concentration after condensation
    REAL(dp), INTENT(in), OPTIONAL :: minNH3(numcells)
    !> Minimum HNO<sub>3</sub>(g) concentration after condensation
    REAL(dp), INTENT(in), OPTIONAL :: minHNO3(numcells)
! op_ck_20130527+
    !> Minimum HCl(g) concentration after condensation
    REAL(dp), INTENT(in), OPTIONAL :: minHCl(numcells)
! op_ck_20130527-
! op_ck_20130408-

    real(dp), parameter :: RH_HIST_DW = 1.50_dp ! mean value for mixture of wet
                                                ! (2) and dry (1) gridboxes
                                                ! (needed for HYSTERESIS)
    real(dp), parameter :: TT0 = 298.15_dp      ! = T0 in EQSAM, avoid conflicts
                                                ! with global parameter T0
!    real(dp), parameter :: TT1 = 298.0_dp
    real(dp), parameter :: R  = 82.0567e-6_dp  ! in cu.m*atm/deg/mole

    real(dp), parameter :: RHMAX = 0.99_dp     ! restrict to max / min RH
    real(dp), parameter :: RHMIN = 0.0001_dp

    real(dp), parameter :: MWH20 = 55.51_dp*18.01_dp ! "H-2-Null"
    real(dp), parameter :: ZERO  = 0.0_dp

    real(dp), parameter :: ZEPS  = 1.0e-19_dp  ! epsilon

    real(dp) :: dum                            ! auxiliary variable

    ! exponents of AC-RH functions
    real(dp), parameter :: GF1 = 0.25_dp
!!$    real(dp), parameter :: GF2 = 0.50_dp
!!$    real(dp), parameter :: GF3 = 0.40_dp
!!$    real(dp), parameter :: GF4 = 1.00_dp
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
    real(dp) :: GSO4,PNa,PSO4
    real(dp) :: ASO4,ANO3,ANH4,ACl,ANa,SNH4,SSO4,SNO3,SCl,SNa
    real(dp) :: scalSS

    !_______________________________________________

    real(dp) :: w1(8), w2(8)
    ! RHD / MRHD arrays for different aerosol types
    real(dp) :: RHDA(8),RHDE(8),RHDX(8),RHDZ(8)
    ! arrays of ion pairs
    real(dp) :: M0(NPAIR),MWSALT(NPAIR),NW(NPAIR),ZW(NPAIR)
    !
    ! salt solutes:
    !   1 = NACl,    2 = (NA)2SO4,     3 = NANO3,  4 = (NH4)2SO4,
    !   5 = NH4NO3,  6 = NH4CL,        7 = 2H-SO4, 8 = NH4HSO4,
    !   9 = NAHSO4, 10 = (NH4)3H(SO4)2
    !
    ! mole mass of the salt solute
    DATA MWSALT / 58.5_dp, 142.0_dp,  88.0_dp, 132.0_dp, 80.0_dp, 53.5_dp,     &
         98.0_dp, 115.0_dp, 120.0_dp, 247.0_dp /
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

    ! op_ck_20130528: Due to removal of Cl from SS, the mass fractions used
    !                 below are no longer appropriate. scalSS is the
    !                 corresponding correction factor. Only 99.28% of the total
    !                 mass are accounted for. Assumptions: Cl tracer contains
    !                 only Cl, SS tracer contains the rest of the SS
    !                 constituents (including unaccounted material) that EQSAM
    !                 assumes.
    scalSS = 1.0_dp / (1.0_dp - mfCl)

    do il = 1, numcells

       TT = temp(il)           ! temperature [K]
       RHL = RH(il)            ! relative humidity [0-1]
       !______________________________________________

       iflag=1

       ! Na+ (ss + xsod) (a) [mol m-3(air)]
!!$       w1(1) = 1.0e-6_dp * (0.3061_dp * inSS(il) / MWNa)
       w1(1) = 1.0e-6_dp * (0.3061_dp * scalSS * inSS(il) / MWNa)
       ! H2SO4    + SO4-- (p)  [mol m-3(air)]
       w1(2) = 1.0e-6_dp &
!!$            * ((inSO4(il) + 0.0768_dp * inSS(il)) / MW(i_so4,i_mwaero))
            * ((inSO4(il) + 0.0768_dp * scalSS * inSS(il)) / MW(i_so4,i_mwaero))
       ! NH3 (g) + NH4+ (p)  [mol m-3(air)]
       w1(3) = 1.0e-6_dp &
            * (inNH3(il) / MW(i_nh3,i_mwgas) + inNH4(il) / MW(i_nh4,i_mwaero))
       ! HNO3 (g) + NO3- (p) [mol m-3(air)]
       w1(4) = 1.0e-6_dp &
            * (inNO3(il) / MW(i_no3,i_mwaero) + inHNO3(il) / MW(i_hno3,i_mwgas))
       ! HCl (g) + Cl- (p) [mol m-3(air)]
       w1(5) = 1.0e-6_dp * (inHCl(il) / MW(i_hcl,i_mwgas) &
            + inCl(il) / MW(i_cl,i_mwaero))
!            + 0.5504_dp * inSS(il) / MW(i_ss,i_mwaero))
!      w1(5) = 1.0e-6_dp * 0.5504_dp * inSS(il) / MWCl
       ! K+ (p)       [mol m-3(air)]
!!$       w1(6) = 1.0e-6_dp * (0.0110_dp * inSS(il) / MWK)
       w1(6) = 1.0e-6_dp * (0.0110_dp * scalSS * inSS(il) / MWK)
       ! Ca++ (p)     [mol m-3(air)]
!!$       w1(7) = 1.0e-6_dp * (0.0116_dp * inSS(il) / MWCa)
       w1(7) = 1.0e-6_dp * (0.0116_dp * scalSS * inSS(il) / MWCa)
       ! Mg++ (p)     [mol m-3(air)]
!!$       w1(8) = 1.0e-6_dp * (0.0369_dp * inSS(il) / MWMg)
       w1(8) = 1.0e-6_dp * (0.0369_dp * scalSS * inSS(il) / MWMg)

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
            iflag=3

       ! SULFATE VERY RICH CASE if (NH4+Na+K+2(Ca+Mg))/SO4 < 1

       if ((w1(1)+w1(3)+w1(6)+2._dp*(w1(7)+w1(8))) <= w1(2))         &
            iflag=4

       ! SULFATE NEUTRAL CASE

       if ((w1(1)+w1(3)+w1(6)+2._dp*(w1(7)+w1(8))) > (2._dp*w1(2)))  &
            iflag=2

       ! SULFATE POOR AND CATION POOR CASE

       if ((w1(1)+w1(6)+2._dp*(w1(7)+w1(8))) > (2._dp*w1(2)))        &
            iflag=1

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

       ! Sulfate rich
       IF (IFLAG == 3) THEN
          IF (RHL <= RHDZ(7)) THEN        ! ACCOUNT FOR MIXTURE OF (NH4)2SO4(s)&
                                          ! NH4HSO4(s) & (NH4)3H(SO4)2(s)
             GG=1.677_dp                  ! (Na)2SO4 & NaHSO4
             !GG=1.5
          ELSEIF (RHL > RHDZ(7).AND.RH(il) <= RHDZ(5)) THEN
             ! MAINLY (Na)2SO4 / (NH4)2SO4(s) & (NH4)3H(SO4)2(s)
             GG=1.75_dp
             !GG=1.5
          ELSEIF (RHL >= RHDZ(5)) THEN    ! (NH4)2SO4(S) & NH4HSO4(S) & SO4-- &
                                          ! HSO4-
             GG=1.5_dp                    ! (Na)2SO4 & NaHSO4
          ENDIF
       ENDIF

       ! Sulfate very rich
       IF (IFLAG == 4) GG=1.0_dp          ! IF SO4 NEUTRALIZED, THEN ONLY AS
                                          ! NaHSO4/NH4HSO4(S) OR HSO4-/H2SO4
       RHD=RHL

       IF (RH_HIST(il) < RH_HIST_DW) THEN ! GET RHD FOR SOLIDS/HYSTERESIS

          ! GET LOWEST DELIQUESCENCE RELATIVE HUMIDITIES ACCORDING TO THE
          ! CONCENTRATION DOMAIN (APROXIMATION) BASED ON RHD / MRHD
          ! ISORROPIA/SCAPE
          w2(1)=1._dp
          w2(2)=1._dp
          w2(3)=1._dp
          w2(4)=1._dp
          w2(5)=1._dp
          w2(6)=1._dp
          w2(7)=1._dp
          w2(8)=1._dp

          ! skip compound in RHD calculation if value is concentration is zero
          ! or rather small
          
          if (w1(1) <= 1.e-12_dp) w2(1)=0._dp
          if (w1(2) <= 1.e-12_dp) w2(2)=0._dp
          if (w1(3) <= 1.e-12_dp) w2(3)=0._dp
          if (w1(4) <= 1.e-12_dp) w2(4)=0._dp
          if (w1(5) <= 1.e-12_dp) w2(5)=0._dp
          if (w1(6) <= 1.e-12_dp) w2(6)=0._dp
          if (w1(7) <= 1.e-12_dp) w2(7)=0._dp
          if (w1(8) <= 1.e-12_dp) w2(8)=0._dp

          ! GET LOWEST RHD ACCORDING TO THE CONCENTRATION DOMAIN

          ! iflag=1 (cation rich, i.e. sulfate poor)  ...
          ! 1. sea salt      aerosol          : RHDX(1)=MgCl2
          ! 2. mineral dust  aerosol          : RHDX(2)=Ca(NO3)2
          !
          ! iflag=2 (sulfate neutral) ...
          ! 3. ammonium + nitrate             : RHDX(3)= NH4NO3
          ! 4. ammonium + sulfate             : RHDX(4)=(NH4)2SO4
          ! 5. ammonium + sulfate mixed salt  : RHDX(5)=
          !                                     (NH4)3H(SO4)2, (NH4)2SO4
          ! 6. ammonium + nitrate  + sulfate  : RHDX(6)=(NH4)2SO4,
          !                                     NH4NO3,NA2SO4, NH4Cl
          !
          ! iflag=3 (sulfate rich) ...
          ! 7. ammonium + sulfate  (1:1,1.5)  : RHDX(7)= NH4HSO4
          !
          ! iflag=4 (sulfate very rich) ...
          ! 8. sulfuric acid                  : RHDX(8)= H2SO4
          
          IF (IFLAG == 1) THEN

             RHD=W2(1)+W2(5)              ! Na+  dependency
             IF(RHD < ZEPS) RHDX(1)=1._dp
             RHD=W2(6)+W2(7)+W2(8)        ! K+/Ca++/Mg++ dependency (incl. ss)
             IF(RHD < ZEPS) RHDX(2)=1._dp

             RHD=MIN(RHDX(1),RHDX(2))

          ELSEIF (IFLAG == 2) THEN

             RHD=W2(3)*W2(4)              ! NH4+ & NO3- dependency
             IF(RHD < ZEPS) RHDX(3)=1._dp
             RHD=W2(2)+W2(3)              ! NH4+ & SO4-- dependency
             IF(ABS(GG - 2._dp) > ZEPS) RHD=0._dp  ! account only for pure
                                                   ! (NH4)2SO4
             IF(RHD < ZEPS) RHDX(4)=1._dp
             RHD=W2(2)+W2(3)              ! NH4+ & SO4-- dependency
             IF(RHD < ZEPS) RHDX(5)=1._dp
             RHD=W2(2)+W2(3)+W2(4)+W2(5)  ! (NH4)2SO4, NH4NO3, Na2SO4
             ! NH4Cl dependency
             IF(RHD < ZEPS) RHDX(6)=1._dp

             RHD=MIN(RHDX(3),RHDX(4),RHDX(5),RHDX(6))

          ELSEIF (IFLAG == 3) THEN

             RHD=W2(2)+W2(3)              ! NH4+ & SO4-- dependency
             IF(RHD < ZEPS) RHDX(7)=1._dp
             RHD=RHDX(7)

          ELSEIF (IFLAG == 4) THEN

             RHD=W2(2)                    ! H2SO4 dependency (assume no dry
                                          ! aerosol)
             IF(RHD < ZEPS) RHDX(8)=1._dp

             RHD=RHDX(8)

          ENDIF ! IFLAG
       ENDIF ! SOLIDS

       ! GET WATER ACTIVITIES ACCORDING TO METZGER, 2000.
       ! FUNCTION DERIVED FROM ZSR RELATIONSHIP DATA (AS USED IN ISORROPIA)

       DUM = (1._dp/RHL-1._dp)
       M0(1)  = ((NW(1) *MWH20/MWSALT(1) *DUM))**ZW(1)
       M0(2)  = ((NW(2) *MWH20/MWSALT(2) *DUM))!**ZW(2)
       M0(3)  = ((NW(3) *MWH20/MWSALT(3) *DUM))!**ZW(3)
       M0(4)  = ((NW(4) *MWH20/MWSALT(4) *DUM))!**ZW(4)
       M0(5)  = ((NW(5) *MWH20/MWSALT(5) *DUM))!**ZW(5)
       M0(6)  = ((NW(6) *MWH20/MWSALT(6) *DUM))!**ZW(6)
       M0(7)  = ((NW(7) *MWH20/MWSALT(7) *DUM))**ZW(7)
       M0(8)  = ((NW(8) *MWH20/MWSALT(8) *DUM))!**ZW(8)
       M0(9)  = ((NW(9) *MWH20/MWSALT(9) *DUM))!**ZW(9)
       M0(10) = ((NW(10)*MWH20/MWSALT(10)*DUM))!**ZW(10)

       ! CALCULATE TEMPERATURE DEPENDENT EQUILIBRIUM CONSTANTS

       T0T=TT0/TT
       COEF=1.0_dp+LOG(T0T)-T0T

       ! EQUILIBRIUM CONSTANT NH4NO3(s) <=> NH3(g) + HNO3(g) [atm^2] (ISORROPIA)

       XK10 = 5.746e-17_dp
       XK10= XK10 * EXP(-74.38_dp*(T0T-1.0_dp) + 6.120_dp*COEF)
       KAN = XK10/(R*TT)/(R*TT)

       ! EQUILIBRIUM CONSTANT  NH4Cl(s) <=> NH3(g) + HCl(g) [atm^2] (ISORROPIA)

       XK6  = 1.086e-16_dp
       XK6 = XK6 * EXP(-71.00_dp*(T0T-1.0_dp) + 2.400_dp*COEF)
       KAC = XK6/(R*TT)/(R*TT)

       ! GET MEAN MOLAL IONIC ACTIVITY COEFF ACCORDING TO METZGER, 2002.

       GAMA=0.0_dp
       ZFLAG=REAL(IFLAG, kind=dp)
       IF (RHL >= RHD) GAMA=(RHL**IFLAG/(1000._dp/ZFLAG*(1._dp-RHL)+ZFLAG))
       !                          ^^^^^
       ! using IFLAG instead of ZFLAG might give better performance
       GAMA = GAMA**GF1 ! ONLY GAMA TYPE OF NH4NO3, NaCl, etc. NEEDED SO FAR

       GAMA=0.0_dp
       GFN=K*K          ! K=2, i.e. condensation of 2 water molecules per 1 mole
                        ! ion pair
       GF=GFN*GF1       ! = GFN[=Nw=4] * GF1[=(1*1^1+1*1^1)/2/Nw=1/4] = 1
                        ! ONLY GAMA TYPE OF NH4NO3, NH4Cl, etc. needed so far

       IF (RHL >= RHD) GAMA=RHL**GF/((GFN*MWH20*(1._dp/RHL-1._dp)))**GF1

       GAMA = MIN(GAMA,1.0_dp) ! FOCUS ON 0-1 SCALE
       GAMA = MAX(GAMA,0.0_dp)
       GAMA = (1._dp-GAMA)**K  ! transplate into aqueous phase equillibrium and
                               ! account for enhanced uptake of aerosol
                               ! precursor gases with increasing RH
                               ! (to match the results of ISORROPIA)

       ! CALCULATE RHD DEPENDENT EQ: IF RH <  RHD =>
       ! NH4NO3(s) <=> NH3 (g) + HNO3(g) (ISORROPIA)
       ! IF RH >> RHD => HNO3  (g)   -> NO3 (aq)

       X00  = MAX(ZERO,MIN(TNa,GG*TSO4))       ! MAX SODIUM SULFATE
       X0   = MAX(ZERO,MIN(TNH4,GG*TSO4-X00))  ! MAX AMMOMIUM SULFATE
       X01  = MAX(ZERO,MIN(TNa-X00, TNO3))     ! MAX SODIUM NITRATE
       X1   = MAX(ZERO,MIN(TNH4-X0,TNO3-X01))  ! MAX AMMOMIUM NITRATE

       X02  = MAX(ZERO,MIN(TNa-X01-X00,TCl))   ! MAX SODIUM CHLORIDE
       X03  = MAX(ZERO,MIN(TNH4-X0-X1,TCl-X02))! MAX AMMOMIUM CHLORIDE

       X2   = MAX(TNH4-X1-X0-X03,ZERO)        ! INTERIM RESIDUAL NH3
       X3   = MAX(TNO3-X1-X01,ZERO)           ! INTERIM RESIDUAL HNO3
       X04  = MAX(TSO4-(X0+X00)/GG,ZERO)      ! INTERIM RESIDUAL H2SO4
       X05  = MAX(TCl-X03-X02,ZERO)           ! INTERIM RESIDUAL HCl
!       X06  = MAX(TNa-X02-X01-X00,ZERO)       ! INTERIM RESIDUAL Na
       ! (should be zero for electro-neutrality in input data)
       !
       ZKAN=2._dp
       IF (RHL >= RHD) ZKAN=ZKAN*GAMA

       X4   = X2 + X3
       X5   = SQRT(X4*X4+KAN*ZKAN*ZKAN)
       X6   = 0.5_dp*(-X4+X5)
       X6   = MIN(X1,X6)

       GHCl(il) = X05                   ! INTERIM RESIDUAl HCl
       GNH3(il) = X2 + X6               ! INTERIM RESIDUAl NH3
       GNO3(il) = X3 + X6               ! RESIDUAl HNO3
       GSO4 = X04                       ! RESIDUAl H2SO4
       PNa  = X02 + X01 + X00           ! RESIDUAl Na (neutralized)

       ZKAC=2._dp
       IF(RHL >= RHD) ZKAC=ZKAC*GAMA

       X08   = GNH3(il) + GHCl(il)
!       X08   = GNH3(il)
       X09   = SQRT(X08*X08+KAC*ZKAC*ZKAC)
       X10   = 0.5_dp*(-X08+X09)
       X11   = MIN(X03,X10)

       GHCl(il) = GHCl(il) + X11        ! RESIDUAL HCl
       GNH3(il) = GNH3(il) + X11        ! RESIDUAL NH3

       ! GO SAVE ...

       IF (GHCl(il) < 0._dp) GHCl(il)=0._dp
       IF (GSO4     < 0._dp) GSO4=0._dp
       IF (GNH3(il) < 0._dp) GNH3(il)=0._dp
       IF (GNO3(il) < 0._dp) GNO3(il)=0._dp
       IF (PNa      < 0._dp) PNa=0._dp
       IF (GSO4     > TSO4)  GSO4=TSO4
       IF (GNH3(il) > TNH4)  GNH3(il)=TNH4
       IF (GNO3(il) > TNO3)  GNO3(il)=TNO3
       IF (GHCl(il) > TCl)   GHCl(il)=TCl
       IF (PNa      > TNa)   PNa=TNa

! op_ck_20130408+
       ! Limit gas-to-particle transfer to maximum.
       ! NH3 <-> NH4+
       IF (PRESENT(minNH3))  GNH3(il) = &
            MAX(GNH3(il), 1.0e-6_dp * minNH3(il)  / MW(i_nh3,i_mwgas))
       ! HNO3 <-> NO3-
       IF (PRESENT(minHNO3)) GNO3(il) = &
            MAX(GNO3(il), 1.0e-6_dp * minHNO3(il) / MW(i_hno3,i_mwgas))
! op_ck_20130408-
! op_ck_20130527+
       ! HCl <-> Cl-
       IF (PRESENT(minHCl)) GHCl(il) = &
            MAX(GHCl(il), 1.0e-6_dp * minHCl(il) / MW(i_hcl,i_mwgas))
! op_ck_20130527-

       ! DEFINE AQUEOUSE PHASE (NO SOLID NH4NO3 IF NO3/SO4 > 1,
       ! TEN BRINK, ET AL., 1996, ATMOS ENV, 24, 4251-4261)
       
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
          ACl  = TCl  - GHCl(il)
!          ACl  = TCl
          ANa  = PNa

       ENDIF ! SOLIDS/HYSTERESIS

       ! CALCULATE AEROSOL WATER [kg m-3(air)]
       !
       ! salt solutes:
       !   1 = NaCl,    2 = (Na)2SO4, 3 = NaNO3,  4 = (NH4)2SO4,
       !   5 = NH4NO3,  6 = NH4Cl,    7 = 2H-SO4, 8 = NH4HSO4,
       !   9 = NaHSO4, 10 = (NH4)3H(SO4)2
       !
       IF (IFLAG == 1) THEN
          WH2O(il)    = ASO4/M0(2) + ANO3/M0(3) + ACl/M0(6)
       ELSE IF (IFLAG == 2) THEN
          WH2O(il)    = ASO4/M0(9) + ANO3/M0(5) + ACl/M0(6)
       ELSE IF (IFLAG == 3) THEN
          WH2O(il)    = ASO4/M0(8) + ANO3/M0(5) + ACl/M0(6)
       ELSE IF (IFLAG == 4) THEN
          WH2O(il)    = ASO4/M0(8) + GSO4/M0(7)
       END IF

       !-------------------------------------------------------
       ! calculate diagnostic output consistent with other EQMs ...

       ASO4 = ASO4 + GSO4  ! assuming H2SO4 remains aqueous

       TNa   = TNa  * 1.e6_dp ! total input sodium   (g+p)  [umol m-3]
       TSO4  = TSO4 * 1.e6_dp ! total input sulfate  (g+p)  [umol m-3]
       TNH4  = TNH4 * 1.e6_dp ! total input ammonium (g+p)  [umol m-3]
       TNO3  = TNO3 * 1.e6_dp ! total input nitrate  (g+p)  [umol m-3]
       TCl   = TCl  * 1.e6_dp ! total input chloride (g+p)  [umol m-3]
       TPo   = TPo  * 1.e6_dp ! total input potasium (g+p)  [umol m-3]
       TCa   = TCa  * 1.e6_dp ! total input calcium  (g+p)  [umol m-3]
       TMg   = TMg  * 1.e6_dp ! total input magnesium(g+p)  [umol m-3]

       ! residual gas:
       GNH3(il) = GNH3(il) * 1.e6_dp  ! residual NH3   [umol m-3]
       GSO4     = GSO4     * 1.e6_dp  ! residual H2SO4 [umol m-3]
       GNO3(il) = GNO3(il) * 1.e6_dp  ! residual HNO3  [umol m-3]
       GHCl(il) = GHCl(il) * 1.e6_dp  ! residual HCl   [umol m-3]

       ! total particulate matter (neutralized)
       PNH4(il)=TNH4-GNH3(il)      ! particulate ammonium   [umol m-3]
       PNO3(il)=TNO3-GNO3(il)      ! particulate nitrate    [umol m-3]
       PCl(il) =TCl -GHCl(il)      ! particulate chloride   [umol m-3]
!       PCl     =TCl                ! particulate chloride   [umol m-3]
       PNa     =TNa                ! particulate sodium     [umol m-3]
       PSO4    =TSO4               ! particulate sulfate    [umol m-3]

       ! liquid matter
       ASO4 = ASO4 * 1.e6_dp  ! aqueous phase sulfate  [umol m-3]
       ANH4 = ANH4 * 1.e6_dp  ! aqueous phase ammonium [umol m-3]
       ANO3 = ANO3 * 1.e6_dp  ! aqueous phase nitrate  [umol m-3]
       ACl  = ACl  * 1.e6_dp  ! aqueous phase chloride [umol m-3]
       ANa  = ANa  * 1.e6_dp  ! aqueous phase sodium   [umol m-3]

       ! solid matter
       SNH4=PNH4(il)-ANH4  ! solid phase ammonium   [umol m-3]
       SSO4=PSO4-ASO4      ! solid phase sulfate    [umol m-3]
       SNO3=PNO3(il)-ANO3  ! solid phase nitrate    [umol m-3]
       SCl =PCl(il) -ACl   ! solid phase chloride   [umol m-3]
       SNa =PNa -ANa       ! solid phase sodium     [umol m-3]

       ! GO SAVE ...
       
       IF (SNH4 < 0._dp) SNH4=0._dp
       IF (SSO4 < 0._dp) SSO4=0._dp
       IF (SNO3 < 0._dp) SNO3=0._dp
       IF (SCl  < 0._dp) SCl =0._dp
       IF (SNa  < 0._dp) SNa =0._dp

       ! convert aerosol water from [kg m-3] to [ug m-3]
       WH2O(il)    = WH2O(il)    * 1.0e9_dp
       IF (WH2O(il) < 1.0e-3_dp) WH2O(il)    = 0.0_dp

       ! UPDATE HISTORY RH FOR HYSTERESIS (ONLINE CALCULATIONS ONLY)

       RH_HIST(il) = 2._dp                                       ! wet
!       IF(WH2O(il) == 0._dp) RH_HIST(il)=1._dp                   ! dry
       IF(WH2O(il) < ZEPS) RH_HIST(il)=1._dp                     ! dry

       !
       ! store aerosol species for diagnostic output:
       !___________________________________________________________
       ! Output values:

       ! convert from [umol m-3] to [ug m-3]
       GNH3(il)    = GNH3(il)    * MW(i_nh3,i_mwgas)
       GNO3(il)    = GNO3(il)    * MW(i_hno3,i_mwgas)
       GHCl(il)    = GHCl(il)    * MW(i_hcl,i_mwgas)
       PNH4(il)    = PNH4(il)    * MW(i_nh4,i_mwaero)
       PNO3(il)    = PNO3(il)    * MW(i_no3,i_mwaero)
       PCl(il)     = PCl(il)     * MW(i_cl,i_mwaero)
    enddo

  END SUBROUTINE MADE3_EQSAM

!-------------------------------------------------------------------------------

  !> \brief Calculates modal parameters and derived variables, as well as
  !>   atmospheric parameters
  !> \details Uses tracer array and air temperature and pressure to compute
  !>   modal 3<sup>rd</sup> moments, total modal mass concentrations, average
  !>   wet aerosol densities, atmospheric mean free path, atmospheric dynamic
  !>   viscosity, wet median mode diameters, and modal Knudsen numbers.
  !> \note Total modal mass concentrations, 3<sup>rd</sup> moments, wet median
  !>   mode diameters, and Knudsen numbers are calculated only from selected
  !>   tracers, i.e., they do not take into account all aerosol species in all
  !>   modes.

  SUBROUTINE MADE3_MODPAR(NUMCELLS, CBLK, BLKTA, BLKPRS &
       , PMASS, PDENS, XLM, AMU, DG, KN)

    IMPLICIT NONE
    INTRINSIC MAX, SQRT

    ! *** input ***
    !> Number of grid cells in input arrays
    INTEGER, INTENT(in)     :: NUMCELLS
    !> Tracer array
    REAL(dp), INTENT(inout) :: CBLK(dim1_cblk,dim2_cblk,NUMCELLS)
    !> Air temperature [K]
    REAL(dp), INTENT(in)    :: BLKTA(NUMCELLS)
    !> Air pressure [Pa]
    REAL(dp), INTENT(in)    :: BLKPRS(NUMCELLS)

    ! *** output ***
    !> Modal mass concentrations [ug m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: PMASS(nmod,NUMCELLS)
    !> Average modal wet aerosol densities [kg m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: PDENS(nmod,NUMCELLS)
    !> Atmospheric mean free path [m]
    REAL(dp), INTENT(out)   :: XLM(NUMCELLS)
    !> Atmospheric dynamic viscosity [kg m<sup>-1</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(out)   :: AMU(NUMCELLS)
    !> Wet median mode diameters
    REAL(dp), INTENT(out)   :: DG(nmod,NUMCELLS)
    !> Modal Knudsen numbers
    REAL(dp), INTENT(out)   :: KN(nmod,NUMCELLS)

    ! *** parameters ***
    REAL(dp), PARAMETER :: DENSMIN = 1.0E03_dp  ! min. particle density [kg m-3]

    ! *** local ***
    INTEGER :: LCELL             ! loop index (gridcells)
    INTEGER :: jm                ! loop index (modes)

    ! Array initializations
    PMASS = -1.0_dp
    PDENS = -1.0_dp
    DG    = -1.0_dp
    KN    = -1.0_dp

    ! *** set up  aerosol  3rd moment, mass, density

    DO  LCELL = 1, NUMCELLS

       ! *** soluble Aitken-mode
       CBLK(i_mom3,akn,LCELL) = MASS2MOM3(i_so4) * CBLK(i_so4,akn,LCELL) &
            + MASS2MOM3(i_nh4) * CBLK(i_nh4,akn,LCELL)                   &
            + MASS2MOM3(i_h2o) * CBLK(i_h2o,akn,LCELL)                   &
            + MASS2MOM3(i_no3) * CBLK(i_no3,akn,LCELL)                   &
!sorgam                      + MASS2MOM3(i_pom) * CBLK(VORGARO1I,LCELL) &
!sorgam                      + MASS2MOM3(i_pom) * CBLK(VORGARO2I,LCELL) &
!sorgam                      + MASS2MOM3(i_pom) * CBLK(VORGALK1I,LCELL) &
!sorgam                      + MASS2MOM3(i_pom) * CBLK(VORGOLE1I,LCELL) &
!sorgam                      + MASS2MOM3(i_pom) * CBLK(VORGBA1I,LCELL)  &
!sorgam                      + MASS2MOM3(i_pom) * CBLK(VORGBA2I,LCELL)  &
            + MASS2MOM3(i_pom) * CBLK(i_pom,akn,LCELL)                   &
!unused                      + MASS2MOM3(i_bc)  * CBLK(VP25AI,LCELL)    &
            + MASS2MOM3(i_cl)  * CBLK(i_cl,akn,LCELL)                    &
            + MASS2MOM3(i_ss)  * CBLK(i_ss,akn,LCELL)
       CBLK(i_mom3,akn,LCELL) = MAX(CONMIN, CBLK(i_mom3,akn,LCELL))

       ! *** soluble accumulation mode
       CBLK(i_mom3,acc,LCELL) = MASS2MOM3(i_so4) *  CBLK(i_so4,acc,LCELL) &
            + MASS2MOM3(i_nh4) *  CBLK(i_nh4,acc,LCELL)                   &
            + MASS2MOM3(i_h2o) *  CBLK(i_h2o,acc,LCELL)                   &
            + MASS2MOM3(i_no3) *  CBLK(i_no3,acc,LCELL)                   &
!sorgam                      + MASS2MOM3(i_pom) *  CBLK(VORGARO1J,LCELL) &
!sorgam                      + MASS2MOM3(i_pom) *  CBLK(VORGARO2J,LCELL) &
!sorgam                      + MASS2MOM3(i_pom) *  CBLK(VORGALK1J,LCELL) &
!sorgam                      + MASS2MOM3(i_pom) *  CBLK(VORGOLE1J,LCELL) &
!sorgam                      + MASS2MOM3(i_pom) *  CBLK(VORGBA1J,LCELL)  &
!sorgam                      + MASS2MOM3(i_pom) *  CBLK(VORGBA2J,LCELL)  &
            + MASS2MOM3(i_pom) *  CBLK(i_pom,acc,LCELL)                   &
!unused                      + MASS2MOM3(i_bc)  *  CBLK(VP25AJ,LCELL)    &
            + MASS2MOM3(i_cl)  *  CBLK(i_cl,acc,LCELL)                    &
            + MASS2MOM3(i_ss)  *  CBLK(i_ss,acc,LCELL)
       CBLK(i_mom3,acc,LCELL) = MAX(CONMIN, CBLK(i_mom3,acc,LCELL))

       ! *** soluble coarse mode
       CBLK(i_mom3,cor,LCELL) = MASS2MOM3(i_so4) * CBLK(i_so4,cor,LCELL) &
            + MASS2MOM3(i_no3) * CBLK(i_no3,cor,LCELL)                   &
            + MASS2MOM3(i_nh4) * CBLK(i_nh4,cor,LCELL)                   &
            + MASS2MOM3(i_cl)  * CBLK(i_cl,cor,LCELL)                    &
            + MASS2MOM3(i_ss)  * CBLK(i_ss,cor,LCELL)                    &
            + MASS2MOM3(i_pom) * CBLK(i_pom,cor,LCELL)                   &
            + MASS2MOM3(i_h2o) * CBLK(i_h2o,cor,LCELL)
       CBLK(i_mom3,cor,LCELL) = MAX(CONMIN, CBLK(i_mom3,cor,LCELL))

       ! *** soluble + insoluble Aitken mode
       CBLK(i_mom3,akns,LCELL) = MASS2MOM3(i_so4) * CBLK(i_so4,akns,LCELL)  &
            + MASS2MOM3(i_nh4) * CBLK(i_nh4,akns,LCELL)                     &
            + MASS2MOM3(i_h2o) * CBLK(i_h2o,akns,LCELL)                     &
            + MASS2MOM3(i_no3) * CBLK(i_no3,akns,LCELL)                     &
!sorgam                       + MASS2MOM3(i_pom) * CBLK(VORGARO1Is,LCELL)   &
!sorgam                       + MASS2MOM3(i_pom) * CBLK(VORGARO2Is,LCELL)   &
!sorgam                       + MASS2MOM3(i_pom) * CBLK(VORGALK1Is,LCELL)   &
!sorgam                       + MASS2MOM3(i_pom) * CBLK(VORGOLE1Is,LCELL)   &
!sorgam                       + MASS2MOM3(i_pom) * CBLK(VORGBA1Is,LCELL)    &
!sorgam                       + MASS2MOM3(i_pom) * CBLK(VORGBA2Is,LCELL)    &
            + MASS2MOM3(i_pom) * CBLK(i_pom,akns,LCELL)                     &
!unused                       + MASS2MOM3(i_bc)  * CBLK(VP25AIs,LCELL)      &
            + MASS2MOM3(i_bc)  * CBLK(i_bc,akns,LCELL)                      &
            + MASS2MOM3(i_bctag) * CBLK(i_bctag,akns,LCELL) & ! op_mr_20181002
            + MASS2MOM3(i_cl)  * CBLK(i_cl,akns,LCELL)                      &
            + MASS2MOM3(i_ss)  * CBLK(i_ss,akns,LCELL)
       CBLK(i_mom3,akns,LCELL) = MAX(CONMIN, CBLK(i_mom3,akns,LCELL))

       ! *** soluble + insoluble accumulation mode
       CBLK(i_mom3,accs,LCELL) = MASS2MOM3(i_so4) * CBLK(i_so4,accs,LCELL)  &
            + MASS2MOM3(i_nh4) *  CBLK(i_nh4,accs,LCELL)                    &
            + MASS2MOM3(i_h2o) *  CBLK(i_h2o,accs,LCELL)                    &
            + MASS2MOM3(i_no3) *  CBLK(i_no3,accs,LCELL)                    &
!sorgam                       + MASS2MOM3(i_pom) *  CBLK(VORGARO1Js,LCELL)  &
!sorgam                       + MASS2MOM3(i_pom) *  CBLK(VORGARO2Js,LCELL)  &
!sorgam                       + MASS2MOM3(i_pom) *  CBLK(VORGALK1Js,LCELL)  &
!sorgam                       + MASS2MOM3(i_pom) *  CBLK(VORGOLE1Js,LCELL)  &
!sorgam                       + MASS2MOM3(i_pom) *  CBLK(VORGBA1Js,LCELL)   &
!sorgam                       + MASS2MOM3(i_pom) *  CBLK(VORGBA2Js,LCELL)   &
            + MASS2MOM3(i_pom) *  CBLK(i_pom,accs,LCELL)                    &
!unused                       + MASS2MOM3(i_bc)  * CBLK(VP25AJs,LCELL)      &
            + MASS2MOM3(i_bc)  * CBLK(i_bc,accs,LCELL)                      &
            + MASS2MOM3(i_bctag) * CBLK(i_bctag,accs,LCELL) & ! op_mr_20181002
            + MASS2MOM3(i_cl)  * CBLK(i_cl,accs,LCELL)                      &
            + MASS2MOM3(i_ss)  * CBLK(i_ss,accs,LCELL)                      &
            + MASS2MOM3(i_du)  * CBLK(i_du,accs,LCELL)
       CBLK(i_mom3,accs,LCELL) = MAX(CONMIN, CBLK(i_mom3,accs,LCELL))

       ! *** soluble + insoluble coarse mode
       CBLK(i_mom3,cors,LCELL) = MASS2MOM3(i_so4) * CBLK(i_so4,cors,LCELL)  &
            + MASS2MOM3(i_no3) * CBLK(i_no3,cors,LCELL)                     &
            + MASS2MOM3(i_nh4) * CBLK(i_nh4,cors,LCELL)                     &
            + MASS2MOM3(i_cl)  * CBLK(i_cl,cors,LCELL)                      &
            + MASS2MOM3(i_ss)  * CBLK(i_ss,cors,LCELL)                      &
            + MASS2MOM3(i_pom) * CBLK(i_pom,cors,LCELL)                     &
            + MASS2MOM3(i_bc)  * CBLK(i_bc,cors,LCELL)                      &
            + MASS2MOM3(i_bctag) * CBLK(i_bctag,cors,LCELL) & ! op_mr_20181002
            + MASS2MOM3(i_du)  * CBLK(i_du,cors,LCELL)                      &
            + MASS2MOM3(i_h2o) * CBLK(i_h2o,cors,LCELL)
       CBLK(i_mom3,cors,LCELL) = MAX(CONMIN, CBLK(i_mom3,cors,LCELL))

       ! *** insoluble externally mixed Aitken mode
       CBLK(i_mom3,sooti,LCELL)= MASS2MOM3(i_so4) * CBLK(i_so4,sooti,LCELL) &
            + MASS2MOM3(i_nh4) * CBLK(i_nh4,sooti,LCELL)                    &
            + MASS2MOM3(i_no3) * CBLK(i_no3,sooti,LCELL)                    &
            + MASS2MOM3(i_h2o) * CBLK(i_h2o,sooti,LCELL)                    &
            + MASS2MOM3(i_pom) * CBLK(i_pom,sooti,LCELL)                    &
            + MASS2MOM3(i_cl)  * CBLK(i_cl,sooti,LCELL)                     &
            + MASS2MOM3(i_ss)  * CBLK(i_ss,sooti,LCELL)                     &
            + MASS2MOM3(i_bc)  * CBLK(i_bc,sooti,LCELL)                     &
            + MASS2MOM3(i_bctag) * CBLK(i_bctag,sooti,LCELL) ! op_mr_20181002
       CBLK(i_mom3,sooti,LCELL) = MAX(CONMIN, CBLK(i_mom3,sooti,LCELL))

       ! *** insoluble externally mixed accumulation mode
       CBLK(i_mom3,sootj,LCELL)= MASS2MOM3(i_so4) * CBLK(i_so4,sootj,LCELL) &
            + MASS2MOM3(i_nh4) * CBLK(i_nh4,sootj,LCELL)                    &
            + MASS2MOM3(i_no3) * CBLK(i_no3,sootj,LCELL)                    &
            + MASS2MOM3(i_h2o) * CBLK(i_h2o,sootj,LCELL)                    &
            + MASS2MOM3(i_pom) * CBLK(i_pom,sootj,LCELL)                    &
            + MASS2MOM3(i_cl)  * CBLK(i_cl,sootj,LCELL)                     &
            + MASS2MOM3(i_ss)  * CBLK(i_ss,sootj,LCELL)                     &
            + MASS2MOM3(i_bc)  * CBLK(i_bc,sootj,LCELL)                     &
            + MASS2MOM3(i_bctag) * CBLK(i_bctag,sootj,LCELL) & ! op_mr_20181002
            + MASS2MOM3(i_du)  * CBLK(i_du,sootj,LCELL)
       CBLK(i_mom3,sootj,LCELL) = MAX(CONMIN, CBLK(i_mom3,sootj,LCELL))

       ! *** insoluble externally mixed coarse mode
       CBLK(i_mom3,sootc,LCELL) = MASS2MOM3(i_so4) * CBLK(i_so4,sootc,LCELL) &
            + MASS2MOM3(i_no3) * CBLK(i_no3,sootc,LCELL)                     &
            + MASS2MOM3(i_nh4) * CBLK(i_nh4,sootc,LCELL)                     &
            + MASS2MOM3(i_cl)  * CBLK(i_cl,sootc,LCELL)                      &
            + MASS2MOM3(i_ss)  * CBLK(i_ss,sootc,LCELL)                      &
            + MASS2MOM3(i_pom) * CBLK(i_pom,sootc,LCELL)                     &
            + MASS2MOM3(i_bc)  * CBLK(i_bc,sootc,LCELL)                      &
            + MASS2MOM3(i_bctag) * CBLK(i_bctag,sootc,LCELL) & ! op_mr_20181002
            + MASS2MOM3(i_du)  * CBLK(i_du,sootc,LCELL)                      &
            + MASS2MOM3(i_h2o) * CBLK(i_h2o,sootc,LCELL)
       CBLK(i_mom3,sootc,LCELL) = MAX(CONMIN, CBLK(i_mom3,sootc,LCELL))

       ! *** now get particle mass and density

       ! *** soluble Aitken mode:
       PMASS(akn,LCELL) = CBLK(i_so4,akn,LCELL) &
            + CBLK(i_nh4,akn,LCELL)             &
            + CBLK(i_h2o,akn,LCELL)             &
            + CBLK(i_no3,akn,LCELL)             &
!sorgam                + CBLK(VORGARO1I,LCELL) &
!sorgam                + CBLK(VORGARO2I,LCELL) &
!sorgam                + CBLK(VORGALK1I,LCELL) &
!sorgam                + CBLK(VORGOLE1I,LCELL) &
!sorgam                + CBLK(VORGBA1I,LCELL)  &
!sorgam                + CBLK(VORGBA2I,LCELL)  &
            + CBLK(i_pom,akn,LCELL)             &
!unused                + CBLK(VP25AI,LCELL)    &
            + CBLK(i_cl,akn,LCELL)              &
            + CBLK(i_ss,akn,LCELL)

       ! *** soluble accumulation mode:
       PMASS(acc,LCELL) = CBLK(i_so4,acc,LCELL) &
            + CBLK(i_nh4,acc,LCELL)             &
            + CBLK(i_h2o,acc,LCELL)             &
            + CBLK(i_no3,acc,LCELL)             &
!sorgam                + CBLK(VORGARO1J,LCELL) &
!sorgam                + CBLK(VORGARO2J,LCELL) &
!sorgam                + CBLK(VORGALK1J,LCELL) &
!sorgam                + CBLK(VORGOLE1J,LCELL) &
!sorgam                + CBLK(VORGBA1J,LCELL)  &
!sorgam                + CBLK(VORGBA2J,LCELL)  &
            + CBLK(i_pom,acc,LCELL)             &
!unused                + CBLK(VP25AJ,LCELL)    &
            + CBLK(i_cl,acc,LCELL)              &
            + CBLK(i_ss,acc,LCELL)

       ! *** soluble coarse mode:
       PMASS(cor,LCELL) = CBLK(i_so4,cor,LCELL) &
            + CBLK(i_no3,cor,LCELL)             &
            + CBLK(i_nh4,cor,LCELL)             &
            + CBLK(i_cl,cor,LCELL)              &
            + CBLK(i_ss,cor,LCELL)              &
            + CBLK(i_pom,cor,LCELL)             &
            + CBLK(i_h2o,cor,LCELL)

       ! *** soluble + insoluble Aitken mode
       PMASS(akns,LCELL) = CBLK(i_so4,akns,LCELL) &
            + CBLK(i_nh4,akns,LCELL)              &
            + CBLK(i_h2o,akns,LCELL)              &
            + CBLK(i_no3,akns,LCELL)              &
!sorgam                 + CBLK(VORGARO1Is,LCELL) &
!sorgam                 + CBLK(VORGARO2Is,LCELL) &
!sorgam                 + CBLK(VORGALK1Is,LCELL) &
!sorgam                 + CBLK(VORGOLE1Is,LCELL) &
!sorgam                 + CBLK(VORGBA1Is,LCELL)  &
!sorgam                 + CBLK(VORGBA2Is,LCELL)  &
            + CBLK(i_pom,akns,LCELL)              &
!unused                 + CBLK(VP25AIs,LCELL)    &
            + CBLK(i_cl,akns,LCELL)               &
            + CBLK(i_ss,akns,LCELL)               &
            + CBLK(i_bc,akns,LCELL)               &
            + CBLK(i_bctag,akns,LCELL) ! op_mr_20181002

       ! *** soluble + insoluble accumulation mode
       PMASS(accs,LCELL) = CBLK(i_so4,accs,LCELL) &
            + CBLK(i_nh4,accs,LCELL)              &
            + CBLK(i_h2o,accs,LCELL)              &
            + CBLK(i_no3,accs,LCELL)              &
!sorgam                 + CBLK(ORGARO1Js,LCELL)    &
!sorgam                 + CBLK(VORGARO2Js,LCELL)   &
!sorgam                 + CBLK(VORGALK1Js,LCELL)   &
!sorgam                 + CBLK(VORGOLE1Js,LCELL)   &
!sorgam                 + CBLK(VORGBA1Js,LCELL)    &
!sorgam                 + CBLK(VORGBA2Js,LCELL)    &
            + CBLK(i_pom,accs,LCELL)              &
!unused                 + CBLK(VP25AJs,LCELL)      &
            + CBLK(i_bc,accs,LCELL)               &
            + CBLK(i_bctag,accs,LCELL)            & ! op_mr_20181002
            + CBLK(i_cl,accs,LCELL)               &
            + CBLK(i_ss,accs,LCELL)               &
            + CBLK(i_du,accs,LCELL)

       ! *** soluble + insoluble coarse mode
       PMASS(cors,LCELL) = CBLK(i_so4,cors,LCELL) &
            + CBLK(i_no3,cors,LCELL)              &
            + CBLK(i_nh4,cors,LCELL)              &
            + CBLK(i_cl,cors,LCELL)               &
            + CBLK(i_ss,cors,LCELL)               &
            + CBLK(i_pom,cors,LCELL)              &
            + CBLK(i_bc,cors,LCELL)               &
            + CBLK(i_bctag,cors,LCELL)            & ! op_mr_20181002
            + CBLK(i_du,cors,LCELL)               &
            + CBLK(i_h2o,cors,LCELL)

       ! *** insoluble externally mixed Aitken mode
       PMASS(sooti,LCELL) = (CBLK(i_bc,sooti,LCELL) &
            +  CBLK(i_bctag,sooti,LCELL)            & ! op_mr_20181002
            +  CBLK(i_so4,sooti,LCELL)              &
            +  CBLK(i_nh4,sooti,LCELL)              &
            +  CBLK(i_no3,sooti,LCELL)              &
            +  CBLK(i_h2o,sooti,LCELL)              &
            +  CBLK(i_pom,sooti,LCELL)              &
            +  CBLK(i_cl,sooti,LCELL)               &
            +  CBLK(i_ss,sooti,LCELL))

       ! *** insoluble externally mixed accumulation mode
       PMASS(sootj,LCELL) = (CBLK(i_du,sootj,LCELL) &
            +  CBLK(i_bc,sootj,LCELL)               &
            +  CBLK(i_bctag,sootj,LCELL)            & ! op_mr_20181002
            +  CBLK(i_so4,sootj,LCELL)              &
            +  CBLK(i_nh4,sootj,LCELL)              &
            +  CBLK(i_no3,sootj,LCELL)              &
            +  CBLK(i_h2o,sootj,LCELL)              &
            +  CBLK(i_pom,sootj,LCELL)              &
            +  CBLK(i_cl,sootj,LCELL)               &
            +  CBLK(i_ss,sootj,LCELL))

       ! *** insoluble externally mixed coarse mode
       PMASS(sootc,LCELL) = CBLK(i_so4,sootc,LCELL) &
            + CBLK(i_no3,sootc,LCELL)               &
            + CBLK(i_nh4,sootc,LCELL)               &
            + CBLK(i_cl,sootc,LCELL)                &
            + CBLK(i_ss,sootc,LCELL)                &
            + CBLK(i_pom,sootc,LCELL)               &
            + CBLK(i_bc,sootc,LCELL)                &
            + CBLK(i_bctag,sootc,LCELL)             & ! op_mr_20181002
            + CBLK(i_du,sootc,LCELL)                &
            + CBLK(i_h2o,sootc,LCELL)

       ! *** Calculate mean free path [m]:
       ! *** 6.6328E-8 is the sea level values given in Table I.2.8
       !     on page 10 of U.S. Standard Atmosphere 1962
       XLM(LCELL) = 6.6328E-8_dp*P0*BLKTA(LCELL) / (TS0*BLKPRS(LCELL))

       ! *** Calculate dynamic viscosity [kg m-1 s-1]:
       ! *** U.S. Standard Atmosphere 1962 page 14 expression for dynamic
       !     viscosity is:
       !     dynamic viscosity =  beta * T * sqrt(T) / (T + S)
       !     where beta = 1.458e-6 [kg sec^-1 K**-0.5], S = 110.4 [K].
       AMU(LCELL) = 1.458E-6_dp * BLKTA(LCELL) * SQRT(BLKTA(LCELL)) / &
            (BLKTA(LCELL) + 110.4_dp)

       ! *** now get particle densities, particle diameters, and Knudsen numbers
       do jm = 1, nmod
          PMASS(jm,LCELL) = MAX(CONMIN, PMASS(jm,LCELL))
          ! *** density in [kg m-3]
          PDENS(jm,LCELL) = MAX(DENSMIN, (F6DPIM9 * PMASS(jm,LCELL) &
               / CBLK(i_mom3,jm,LCELL)))
          ! Standard deviation fixed in all modes, so diagnose diameter from 3rd
          ! moment and number concentrations:
! op_ck_20121106  FIXME: Why should we keep larger particles and not reset their
!                        size as well in case PMASS is ``too low''?
! op_ck_20140612         -> We shouldn't, and we should test a different
!                           treatment in the future.
          DG(jm,LCELL) = MAX(DGMIN, &
               (CBLK(i_mom3,jm,LCELL)/(CBLK(i_num,jm,LCELL) * ES36(jm)))**ONE3)
          ! calculate Knudsen numbers
          KN(jm,LCELL) = 2.0_dp * XLM(LCELL) / DG(jm,LCELL)
       end do

    END DO ! end loop over cells

  END SUBROUTINE MADE3_MODPAR

!-------------------------------------------------------------------------------

  !> \brief Calculates aerosol coagulation number rate coefficients and
  !>   3<sup>rd</sup> moment coagulation rates
  !> \details Calculates number concentration reduction rate coefficients for
  !>   intramodal coagulation and number transfer rate coefficients as well as
  !>   3<sup>rd</sup> moment transfer rates for intermodal coagulation according
  !>   to \cite Whitby1991.
  !> \note Uses fixed values for the corrections to the free-molecular
  !>   coagulation integrals.

  SUBROUTINE MADE3_COAGRATE(NUMCELLS, CBLK, BLKTA, PDENS, AMU, DG, KN, CR0, C30)

    IMPLICIT NONE
    INTRINSIC SQRT

    ! *** input ***
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)  :: NUMCELLS
    !> Tracer array
    REAL(dp), INTENT(in)  :: CBLK(dim1_cblk,dim2_cblk,NUMCELLS)
    !> Air temperature [K]
    REAL(dp), INTENT(in)  :: BLKTA(NUMCELLS)
    !> Average modal wet aerosol densities [kg m<sup>-3</sup>]
    REAL(dp), INTENT(in)  :: PDENS(nmod,NUMCELLS)
    !> Atmospheric dynamic viscosity [kg m<sup>-1</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(in)  :: AMU(NUMCELLS)
    !> Wet median modal diameters [m]
    REAL(dp), INTENT(in)  :: DG(nmod,NUMCELLS)
    !> Modal Knudsen numbers
    REAL(dp), INTENT(in)  :: KN(nmod,NUMCELLS)

    ! *** output ***
    !> Modal number coagulation rate coefficients [m<sup>3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(out) :: CR0(nmod,nmod,NUMCELLS)
    !> Intermodal 3<sup>rd</sup> moment transfer rates [mom_3 m<sup>-3</sup>
    !> s<sup>-1</sup>]
    REAL(dp), INTENT(out) :: C30(nmod,nmod,NUMCELLS)

    ! modal variables
    REAL(dp) :: SQDG(nmod)   ! sqrt(diameter)
    REAL(dp) :: DG3(nmod)    ! diameter**3

    ! helper variables
    REAL(dp) :: KNC
    REAL(dp) :: KFM
    REAL(dp) :: RAT
    REAL(dp) :: RIN
    REAL(dp) :: RSQT
    REAL(dp) :: RSQ4
    REAL(dp) :: RSQTI
    REAL(dp) :: RSQI3

    ! process variables
    ! ... coagulation number transfer rate coefficients
    REAL(dp) :: BENCN(nmod,nmod)    ! near continuum number coefficients
    REAL(dp) :: BEFMN(nmod,nmod)    ! free molecular number coefficients
    ! ... coagulation 3rd moment transfer rate coefficients
    REAL(dp) :: BENCM3(nmod,nmod)   ! near continuum 3rd moment coefficients
    REAL(dp) :: BEFMM3(nmod,nmod)   ! free molecular 3rd moment coefficients
    REAL(dp) :: BR31(nmod,nmod)     ! total 3rd moment coefficients

    ! loop indices
    INTEGER :: jm, jm2
    INTEGER :: LCELL

    ! *** Fixed values for corrections to coagulation integrals
    !     for free-molecular (FM) case (must be 64 bit ---> REAL*8). ***
    REAL(dp), PARAMETER :: BM0  = 0.8_dp   ! 0.8D0
    REAL(dp), PARAMETER :: BM0I = 0.9_dp   ! 0.9D0
    REAL(dp), PARAMETER :: BM3I = 0.9_dp   ! 0.9D0
    REAL(dp), PARAMETER :: A    = 1.246_dp ! approx Cunningham corr. factor


    ! -------------------------------- Start code ------------------------------

    ! Array initializations
    SQDG   = -1.0_dp
    DG3    = -1.0_dp
    BENCN  = -1.0_dp
    BEFMN  = -1.0_dp
    BENCM3 = -1.0_dp
    BEFMM3 = -1.0_dp
    BR31   = -1.0_dp
    CR0    =  0.0_dp
    C30    =  0.0_dp

    ! Main computational grid-traversal loops for computing coagulation rates.
    ! *** All modes have fixed std devs. ***
    DO LCELL = 1, NUMCELLS     !  loop on LCELL

       ! ************ begin calculations **************

       KNC = TWO3 * BOLTZ *  BLKTA(LCELL) / AMU(LCELL)

       DO jm = 1, nmod
          ! helper variables
          SQDG(jm) = SQRT(DG(jm,LCELL))
          DG3(jm)  = DG(jm,LCELL) * DG(jm,LCELL) * DG(jm,LCELL)
       END DO

       DO jm = 1, nmod
          DO jm2 = 1,nmod

             ! helper variables
             KFM   = SQRT(6.0_dp * BOLTZ * BLKTA(LCELL) &
                  / (PDENS(jm2,LCELL) + PDENS(jm,LCELL)))
             RAT   = DG(jm,LCELL) / DG(jm2,LCELL) ! Dp1 / Dp2
             RIN   = 1.0_dp / RAT                 ! Dp2 / Dp1
             RSQT  = SQRT(RAT)                    ! sqrt(Dp1/Dp2)
             RSQTI = 1.0_dp / RSQT                ! sqrt(Dp2/Dp1)
             RSQ4  = RAT * RAT                    ! Dp1**2 / Dp2**2
             RSQI3 = RIN * RSQTI                  ! (Dp2/Dp1)**1.5

             ! near continuum bimodal number coagulation rate coefficients
             BENCN(jm2,jm) = KNC * (2.0_dp                          &
                  + A * KN(jm2,LCELL)                               &
                      * (ES04(jm2) + RAT * ES16(jm2) * ES04(jm))    &
                  + A * KN(jm,LCELL)                                &
                      * (ES04(jm)  + RIN * ES16(jm)  * ES04(jm2))   &
                  + (RAT + RIN) * ES04(jm2) * ES04(jm))
             ! free molecular bimodal number coagulation rate coefficients
             BEFMN(jm2,jm) = KFM * BM0I * SQDG(jm2)                 &
                  * (E1(jm2) + RSQT * E1(jm)                        &
                     + 2.0_dp * RAT   * E1(jm2)   * ES04(jm)        &
                     +          RSQ4  * ES09(jm2) * ES16(jm)        &
                     +          RSQI3 * ES16(jm2) * ES09(jm)        &
                     + 2.0_dp * RSQTI * ES04(jm2) * E1(jm))
             ! near continuum bimodal 3rd moment coagulation rate coefficients
             BENCM3(jm2,jm) = KNC * DG3(jm2)                        &
                  * (2.0_dp * ES36(jm2)                             &
                     + A * KN(jm2,LCELL)                            &
                         * (ES16(jm2) + RAT * ES04(jm2) * ES04(jm)) &
                     + A * KN(jm,LCELL)                             &
                         * (ES36(jm2) * ES04(jm)                    &
                            + RIN * ES64(jm2) * ES16(jm))           &
                     + RAT * ES16(jm2) * ES04(jm)                   &
                     + RIN * ES64(jm2) * ES04(jm))
             ! free molecular bimodal 3rd moment coagulation rate coefficients
             BEFMM3(jm2,jm) =                                       &
                  KFM * BM3I * SQDG(jm2) * DG3(jm2)                 &
                  * (ES49(jm2) + RSQT * ES36(jm2) * E1(jm)          &
                     + 2.0_dp * RAT   * ES25(jm2)  * ES04(jm)       &
                     +          RSQ4  * ES09(jm2)  * ES16(jm)       &
                     +          RSQI3 * ES100(jm2) * ES09(jm)       &
                     + 2.0_dp * RSQTI * ES64(jm2)  * E1(jm))
             ! total bimodal 3rd moment transfer rate coefficients
             BR31(jm2,jm) = BENCM3(jm2,jm) * BEFMM3(jm2,jm)         &
                  / (BENCM3(jm2,jm) + BEFMM3(jm2,jm))

          END DO
       END DO

       DO jm = 1, nmod
          ! near continuum unimodal number coagulation rate coefficients
          BENCN(jm,jm) = 0.5_dp * BENCN(jm,jm)
          ! free molecular unimodal number coagulation rate coefficients
          BEFMN(jm,jm) = 0.5_dp * BM0 / BM0I * BEFMN(jm,jm)              
          ! total unimodal 3rd moment transfer rate coefficients
          BR31(jm,jm) = 0.0_dp
       END DO
      
       ! ************ end calculations **************
       
       ! ************ begin assignments **************
       
       ! number transfer coagulation rate coefficients
       DO jm = 1, nmod
          DO jm2 = 1, jm
             CR0(jm2,jm,LCELL) = BENCN(jm2,jm) * BEFMN(jm2,jm) &
                  / (BENCN(jm2,jm) + BEFMN(jm2,jm))
          END DO
       END DO
       ! Ensure that CR0 is symmetric
       DO jm = 1, nmod - 1
          DO jm2 = jm + 1, nmod
             CR0(jm2,jm,LCELL) = CR0(jm,jm2,LCELL)
          END DO
       END DO

       ! bimodal 3rd moment transfer coagulation rates
       ! C30(X,Y,L) is the transfer rate of X by coagulation with Y in cell L
       ! NOTE: The values assigned here may not be the final transfer rates
       !       because target modes for some of the coagulation processes are
       !       computed later (in MADE3_AEROSTEP).
       DO jm2 = 1, nmod
          DO jm = 1, nmod
             C30(jm,jm2,LCELL) = BR31(jm,jm2) &
                  * CBLK(i_num,jm2,LCELL) * CBLK(i_num,jm,LCELL)
          END DO
       END DO
       ! In these processes there is no 3rd moment transfer from the first mode.
       C30(acc,akn,LCELL)     = 0.0_dp
       C30(akns,akn,LCELL)    = 0.0_dp
       C30(accs,akn,LCELL)    = 0.0_dp
       C30(cor,akn,LCELL)     = 0.0_dp
       C30(cors,akn,LCELL)    = 0.0_dp
       C30(sootc,akn,LCELL)   = 0.0_dp
       C30(accs,akns,LCELL)   = 0.0_dp
       C30(cors,akns,LCELL)   = 0.0_dp
       C30(sootc,akns,LCELL)  = 0.0_dp
       C30(sootj,sooti,LCELL) = 0.0_dp
       C30(cors,sooti,LCELL)  = 0.0_dp
       C30(sootc,sooti,LCELL) = 0.0_dp
       C30(accs,acc,LCELL)    = 0.0_dp
       C30(cor,acc,LCELL)     = 0.0_dp
       C30(cors,acc,LCELL)    = 0.0_dp
       C30(cors,accs,LCELL)   = 0.0_dp
       C30(sootc,sootj,LCELL) = 0.0_dp
       C30(cors,cor,LCELL)    = 0.0_dp

       ! ************ end assignments **************

    END DO ! end of main loop over cells

  END SUBROUTINE MADE3_COAGRATE

!-------------------------------------------------------------------------------

  !> \brief Calculates condensational 3<sup>rd</sup> moment growth rates and
  !>   amount of condensed SO<sub>4</sub>
  !> \details Calls the routine that calculates condensation factors for the
  !>   individual modes, then integrates H<sub>2</sub>SO<sub>4</sub>(g) ODE to
  !>   determine total amount that condenses, then sets 3<sup>rd</sup> moment
  !>   growth rate according to H<sub>2</sub>SO<sub>4</sub> and (externally
  !>   supplied) SOA condensation rates. Furthermore, the
  !>   H<sub>2</sub>SO<sub>4</sub>(g) rate of change is reduced by the
  !>   condensation term in order not to overestimate the nucleation term (that
  !>   is computed later).

  SUBROUTINE MADE3_CONDENSE(NUMCELLS, CBLK, DT, BLKTA, BLKPRS    &
       , SO4RAT, SOA_MADE3, DG, FCONC, FCONC_ORG, DELTASO4A, CGR3, VAPOR2 &
!sorgam                   ORGARO1RAT, ORGARO2RAT,       &
!sorgam                   ORGALK1RAT, ORGOLE1RAT,       &
!sorgam                   ORGBIO1RAT, ORGBIO2RAT,       &
!sorgam                   ORGBIO3RAT, ORGBIO4RAT,       &
!sorgam                   DROG, LDROG, NCV, NACV,       &
       )

    IMPLICIT NONE
    INTRINSIC EXP, MAX, MIN

    ! *** arguments ***

    ! *** input ***
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)    :: NUMCELLS
    !> Tracer array
    REAL(dp), INTENT(in)    :: CBLK(dim1_cblk,dim2_cblk,NUMCELLS)
    !> Time step [s]
    REAL(dp), INTENT(in)    :: DT
    !> Air temperature [K]
    REAL(dp), INTENT(in)    :: BLKTA(NUMCELLS)
    !> Air pressure [Pa]
    REAL(dp), INTENT(in)    :: BLKPRS (NUMCELLS)
    !> H<sub>2</sub>SO<sub>4</sub>(g) rate of change [ug m<sup>-3</sup>
    !> s<sup>-1</sup>]
    REAL(dp), INTENT(inout) :: SO4RAT(NUMCELLS)
    !> SOA precursor emissions [ug m<sup>-3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: SOA_MADE3(NUMCELLS)
    !> Wet median modal diameters [m]
    REAL(dp), INTENT(in)    :: DG(nmod,NUMCELLS)
!sorgam   REAL(dp) :: ORGARO1RAT(NUMCELLS)  ! anthropogenic organic aerosol mass
!sorgam                                    ! production rate from aromatics
!sorgam                                    ! [ug/m3/s]
!sorgam   REAL(dp) :: ORGARO2RAT(NUMCELLS)  ! anthropogenic organic aerosol mass
!sorgam                                    ! production rate from aromatics
!sorgam                                    ! [ug/m3/s]
!sorgam   REAL(dp) :: ORGALK1RAT(NUMCELLS)  ! anthropogenic organic aerosol mass
!sorgam                                    ! production rate from alkanes &
!sorgam                                    ! others [ug/m3/s]
!sorgam   REAL(dp) :: ORGOLE1RAT(NUMCELLS)  ! anthropogenic organic aerosol mass
!sorgam                                    ! production rate from alkenes &
!sorgam                                    ! others [ug/m3/s]
!sorgam   !bs * biogenic organic condensable vapor production rate
!sorgam   REAL(dp) :: ORGBIO1RAT(NUMCELLS)  ! biogenic organic aerosol
!sorgam                                    ! production rate [ug/m3/s]
!sorgam   REAL(dp) :: ORGBIO2RAT(NUMCELLS)  ! biogenic organic aerosol
!sorgam                                    ! production rate [ug/m3/s]
!sorgam   REAL(dp) :: ORGBIO3RAT(NUMCELLS)  ! biogenic organic aerosol
!sorgam                                    ! production rate [ug/m3/s]
!sorgam   REAL(dp) :: ORGBIO4RAT(NUMCELLS)  ! biogenic organic aerosol
!sorgam                                    ! production rate [ug/m3/s]
!sorgam   !bs * anthropogenic organic condensable vapor production rate
!sorgam   REAL(dp) :: DROG(NUMCELLS,LDROG)  ! Delta ROG conc. [ppm]
!sorgam   INTEGER,  INTENT(in)    :: LDROG          ! # of organic aerosol
!sorgam                                             ! precursor
!sorgam   INTEGER :: NCV             ! total # of cond. vapors & SOA species
!sorgam   INTEGER :: NACV            ! # of anthrop. cond. vapors & SOA species

    ! *** output ***
    !> H<sub>2</sub>SO<sub>4</sub> condensation rate coefficients
    !> [s<sup>-1</sup>]
    REAL(dp), INTENT(out) :: FCONC(nmod,NUMCELLS)
    !> SOA precursor condensation rate coefficients [s<sup>-1</sup>]
    REAL(dp), INTENT(out) :: FCONC_ORG(nmod,NUMCELLS)
    !> Condensing mass of SO<sub>4</sub> [ug m<sup>-3</sup>]
    REAL(dp), INTENT(out) :: DELTASO4A(NUMCELLS)
    !> Modal 3<sup>rd</sup> moment growth rates (by condensation) [mom_3
    !> m<sup>-3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(out) :: CGR3(nmod,NUMCELLS)
    !> H<sub>2</sub>SO<sub>4</sub>(g) conc. prior to chemistry and emissions,
    !> after cond. [ug m<sup>-3</sup>]
    REAL(dp), INTENT(out) :: VAPOR2(NUMCELLS)

    ! *** local ***
    INTEGER  :: LCELL               ! LOOP INDEX
    REAL(dp) :: CONDRATE            ! condensation rate [mom-3/g/s]
    REAL(dp) :: CHEMRAT_ORG         ! conv rate for organics [mom-3/g/s]
    REAL(dp) :: FCONC_TOT(NUMCELLS)  ! total SO4 condensation factor
    ! *** variables to set up sulfate condensation rate
    REAL(dp) :: OLDSULF(NUMCELLS)    ! "old" conc. of sulfuric acid vapor [ug/m3]
    REAL(dp) :: DELTAVAP            ! change to vapor at previous time step
                                    ! incl. condensation and chem. production
    INTEGER :: jm                   ! loop index


    ! -------------------------------- Start code ------------------------------

    ! Calculate condensation factors
    CALL cond_factors(i_soa, NUMCELLS, BLKPRS, BLKTA   &
         , CBLK(i_num,1:nmod,1:NUMCELLS), DG, FCONC_TOT, FCONC_ORG)
    CALL cond_factors(i_h2so4, NUMCELLS, BLKPRS, BLKTA &
         , CBLK(i_num,1:nmod,1:NUMCELLS), DG, FCONC_TOT, FCONC)

    ! Main computational grid-traversal loop nest for computing condensation:
    DO LCELL = 1, NUMCELLS

       ! *** calculate the total change to sulfuric acid vapor from production
       !     and condensation
       ! OLDSULF: vapor at prev. time step
       OLDSULF(LCELL) = MAX(0.0_dp,CBLK(i_h2so4,gas,LCELL) - SO4RAT(LCELL) * DT)
       ! DELTAVAP: change to vapor at previous time step including condensation
       !           and chem. production
       ! DELTAVAP = [H2SO4](t0 + DT) - [H2SO4](t0), where
       ! [H2SO4](t0) = OLDSULF and [H2SO4](t0 + DT) is the solution to the ODE
       ! d[H2SO4]/dt = P - L[H2SO4], where
       ! P = SO4RAT and L = FCONC_TOT, evaluated at time t0 + DT
       DELTAVAP = (SO4RAT(LCELL) / FCONC_TOT(LCELL) - OLDSULF(LCELL)) &
            * (1.0_dp - EXP(-FCONC_TOT(LCELL) * DT))
       ! VAPOR2: H2SO4(g) conc. taking into account changes due to
       !         condensation and chemical production
       VAPOR2(LCELL) = MAX(0.0_dp, OLDSULF(LCELL) + DELTAVAP)

       ! *** Calculate increment in total sulfate aerosol mass concentration
       !     (by condensation)
       ! *** This follows the method of Youngblood & Kreidenweis.
       ! Using the terminology from above, the next code line calculates
       ! P * DT - ([H2SO4](t0 + DT) - [H2SO4](t0)).
       ! This is the amount of H2SO4 that condenses on the aerosol during time
       ! step t0 -> t0 + DT, according to the solution of the ODE
       ! d[H2SO4]_a/dt = L * [H2SO4]_g(t),
       ! where the indices on the square brackets denote the aerosol and gas
       ! phases, respectively, and [H2SO4]_g(t) is the solution to the ODE given
       ! in the comment above.
       ! Note that this scheme is not applied to SOA because *all* the SOA
       ! ``produced'' during the time step has to ``condense'' (since it does
       ! not contribute to nucleation in MADE3).
       DELTASO4A(LCELL) = MIN(SO4RAT(LCELL) * DT - DELTAVAP, &
            CBLK(i_h2so4,gas,LCELL))
       ! The following step assures that no H2SO4 is transferred from the
       ! aerosol to the gas phase and converts the transferred mass
       ! concentration from H2SO4 to SO4.
       DELTASO4A(LCELL) = MAX(0.0_dp, DELTASO4A(LCELL) &
            * MW(i_so4,i_mwaero) / MW(i_h2so4,i_mwgas))

       !al   Redefine SO4RAT to be maximum mass production rate allowed
       !al   for nucleation:
       !al   = (available H2SO4 vapor - amount condensed) / DT
       !al   This ensures no more H2SO4 vapor to be nucleated than available.
       SO4RAT(LCELL) = &
            (CBLK(i_h2so4,gas,LCELL) * MW(i_so4,i_mwaero) / MW(i_h2so4,i_mwgas)&
            - DELTASO4A(LCELL)) / DT

       ! *** Now calculate the rates of condensation on existing particles.
       CONDRATE = MASS2MOM3(i_so4) * DELTASO4A(LCELL) / DT ! [mom_3 m-3 s-1]
!sorgam      CHEMRAT_ORG = MASS2MOM3(i_pom) * ( ORGARO1RAT(LCELL)
!sorgam                               + ORGARO2RAT(LCELL)
!sorgam                               + ORGALK1RAT(LCELL)
!sorgam                               + ORGOLE1RAT(LCELL)
!sorgam                               + ORGBIO1RAT(LCELL)
!sorgam                               + ORGBIO2RAT(LCELL)
!sorgam                               + ORGBIO3RAT(LCELL)
!sorgam                               + ORGBIO4RAT(LCELL) )
!sorgam                                                ! [mom_3 m-3 s-1]
       CHEMRAT_ORG = MASS2MOM3(i_pom) * SOA_MADE3(LCELL)    ! [mom_3 m-3 s-1]
       DO jm = 1, nmod
          CGR3(jm,LCELL) = CONDRATE * FCONC(jm,LCELL) &
               + CHEMRAT_ORG * FCONC_ORG(jm,LCELL)
       END DO

    END DO  ! loop over NUMCELLS

  END SUBROUTINE MADE3_CONDENSE

!-------------------------------------------------------------------------------

  !> \brief Calculates new particle formation rate (number, mass, and
  !>   3<sup>rd</sup> moment)
  !> \details Uses parameterization by \cite Vehkamaki2002, \cite Vehkamaki2013
  !>   for binary homogeneous nucleation of H<sub>2</sub>SO<sub>4</sub> and
  !>   H<sub>2</sub>O to calculate number, mass, and 3<sup>rd</sup> moment rates
  !>   of change due to new particle formation.
  !> \note Validity ranges:
  !>   - Temperature: 190.15 K - 300.15 K
  !>   - Relative humidity: 0.0001 - 1
  !>   - H<sub>2</sub>SO<sub>4</sub>(g) concentration: 10<sup>4</sup> -
  !>     10<sup>11</sup> cm<sup>-3</sup>
  !>   - Nucleation rate: 10<sup>-7</sup> - 10<sup>10</sup> cm<sup>-3</sup>
  !>     s<sup>-1</sup>

  SUBROUTINE made3_nucleate(blkta, blkrh, h2so4, so4rat, numcells &
       , ndot1, mdot1, cgr3)

    IMPLICIT NONE
    INTRINSIC EXP, LOG, MAX, MIN

    ! *** input ***
    !> Number of grid cells in input arrays
    integer,  intent(in)  :: numcells
    !> Air temperature [K]
    real(dp), intent(in)  :: blkta(numcells)
    !> Relative humidity [-]
    real(dp), intent(in)  :: blkrh(numcells)
    !> H<sub>2</sub>SO<sub>4</sub>(g) concentration [ug m<sup>-3</sup>]
    real(dp), intent(in)  :: h2so4(numcells)
    !> Maximum allowed SO<sub>4</sub>(a) production rate [ug m<sup>-3</sup>
    !> s<sup>-1</sup>]
    real(dp), intent(in)  :: so4rat(numcells)

    ! *** output ***
    !> Nucleation rate (number) [m<sup>-3</sup> s<sup>-1</sup>]
    real(dp), intent(out) :: ndot1(numcells)
    !> Nucleation rate (mass SO<sub>4</sub>) [ug m<sup>-3</sup> s<sup>-1</sup>]
    real(dp), intent(out) :: mdot1(numcells)
    !> Modal 3<sup>rd</sup> moment growth rates [mom_3 m<sup>-3</sup>
    !> s<sup>-1</sup>]
    REAL(dp), INTENT(inout) :: cgr3(nmod,numcells)

    ! *** local ***
    double precision rhoa  ! H2SO4(g) [molecules/m**3]
    double precision t     ! temperature [k]
    double precision rh    ! relative humidity [frac]
    double precision x, nac, jnuc
    double precision ntot
    double precision logrh, logrh2, logrh3
    double precision t2, t3
    double precision logrhoa, logrhoa2, logrhoa3
    real(dp) :: mp         ! mass of sulfate in a 3.5 nm particle
    integer :: ncell       ! loop index

    !*** constants ***
    ! conversion factor molecules ---> ug
    double precision :: molec2ug
    ! conversion factor ug ---> molecules
    double precision :: ug2molec
    real(dp) :: mp35       ! arithmetic statement function to compute
    ! the mass of sulfate in a 3.5 nm particle
    ! coefficients for cubic in mp35
    real(dp), parameter :: a0 =  1.961385e2_dp
    real(dp), parameter :: a1 = -5.564447e2_dp
    real(dp), parameter :: a2 =  8.828801e2_dp
    real(dp), parameter :: a3 = -5.231409e2_dp

    LOGICAL :: firstime
    DATA firstime / .TRUE. /
    SAVE firstime

    SAVE molec2ug, ug2molec

    ! *** function for the number of molecules of sulfate in a 3.5 nm sphere
    ! *** obtained from a fit to the number of sulfate monomers in
    !     a 3.5 nm particle. Uses data from Nair & Vohra.
    double precision rr ! dummy variable for statement function
    mp35(rr) = a0 + rr * (a1 + rr * (a2 + rr * a3))


    ! -------------------------------- Start code ------------------------------

    IF (firstime) THEN
       molec2ug = MW(i_h2so4,i_mwgas) * 1.0e6_dp / avo
       ug2molec = avo * 1.0e-6_dp / MW(i_h2so4,i_mwgas)
       firstime = .FALSE.
    END IF

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
       x = 0.7409967177282139_dp                 &
            - 0.002663785665140117_dp*t          &
            + 0.002010478847383187_dp*logrh      &
            - 0.0001832894131464668_dp*t*logrh   &
            + 0.001574072538464286_dp*logrh2     &
            - 0.00001790589121766952_dp*t*logrh2 &
            + 0.0001844027436573778_dp*logrh3    &
            - 1.503452308794887e-6_dp*t*logrh3   &
            - 0.003499978417957668_dp*logrhoa    &
            + 0.0000504021689382576_dp*t*logrhoa

       ! nucleation rate exponent
       jnuc = 0.1430901615568665_dp                       &
            + 2.219563673425199_dp*t                      &
            - 0.02739106114964264_dp*t2                   &
            + 0.00007228107239317088_dp*t3                &
            + 5.91822263375044_dp/x                       &
            + 0.1174886643003278_dp*logrh                 &
            + 0.4625315047693772_dp*t*logrh               &
            - 0.01180591129059253_dp*t2*logrh             &
            + 0.0000404196487152575_dp*t3*logrh           &
            + (15.79628615047088_dp*logrh)/x              &
            - 0.215553951893509_dp*logrh2                 &
            - 0.0810269192332194_dp*t*logrh2              &
            + 0.001435808434184642_dp*t2*logrh2           &
            - 4.775796947178588e-6_dp*t3*logrh2           &
            - (2.912974063702185_dp*logrh2)/x             &
            - 3.588557942822751_dp*logrh3                 &
            + 0.04950795302831703_dp*t*logrh3             &
            - 0.0002138195118737068_dp*t2*logrh3          &
            + 3.108005107949533e-7_dp*t3*logrh3           &
            - (0.02933332747098296_dp*logrh3)/x           &
            + 1.145983818561277_dp*logrhoa                &
            - 0.6007956227856778_dp*t*logrhoa             &
            + 0.00864244733283759_dp*t2*logrhoa           &
            - 0.00002289467254710888_dp*t3*logrhoa        &
            - (8.44984513869014_dp*logrhoa)/x             &
            + 2.158548369286559_dp*logrh*logrhoa          &
            + 0.0808121412840917_dp*t*logrh*logrhoa       &
            - 0.0004073815255395214_dp*t2*logrh*logrhoa   &
            - 4.019572560156515e-7_dp*t3*logrh*logrhoa    &
            + (0.7213255852557236_dp*logrh*logrhoa)/x     &
            + 1.62409850488771_dp*logrh2*logrhoa          &
            - 0.01601062035325362_dp*t*logrh2*logrhoa     &
            + 0.00003771238979714162_dp*t2*logrh2*logrhoa &
            + 3.217942606371182e-8_dp*t3*logrh2*logrhoa   &
            - (0.01132550810022116_dp*logrh2*logrhoa)/x   &
            + 9.71681713056504_dp*logrhoa2                &
            - 0.1150478558347306_dp*t*logrhoa2            &
            + 0.0001570982486038294_dp*t2*logrhoa2        &
            + 4.009144680125015e-7_dp*t3*logrhoa2         &
            + (0.7118597859976135_dp*logrhoa2)/x          &
            - 1.056105824379897_dp*logrh*logrhoa2         &
            + 0.00903377584628419_dp*t*logrh*logrhoa2     &
            - 0.00001984167387090606_dp*t2*logrh*logrhoa2 &
            + 2.460478196482179e-8_dp*t3*logrh*logrhoa2   &
            - (0.05790872906645181_dp*logrh*logrhoa2)/x   &
            - 0.1487119673397459_dp*logrhoa3              &
            + 0.002835082097822667_dp*t*logrhoa3          &
            - 9.24618825471694e-6_dp*t2*logrhoa3          &
            + 5.004267665960894e-9_dp*t3*logrhoa3         &
            - (0.01270805101481648_dp*logrhoa3)/x

       ! total number of molecules in the critical cluster
       ntot = -0.002954125078716302_dp                    &
            - 0.0976834264241286_dp*t                     &
            + 0.001024847927067835_dp*t2                  &
            - 2.186459697726116e-6_dp*t3                  &
            - 0.1017165718716887_dp/x                     &
            - 0.002050640345231486_dp*logrh               &
            - 0.007585041382707174_dp*t*logrh             &
            + 0.0001926539658089536_dp*t2*logrh           &
            - 6.70429719683894e-7_dp*t3*logrh             &
            - (0.2557744774673163_dp*logrh)/x             &
            + 0.003223076552477191_dp*logrh2              &
            + 0.000852636632240633_dp*t*logrh2            &
            - 0.00001547571354871789_dp*t2*logrh2         &
            + 5.666608424980593e-8_dp*t3*logrh2           &
            + (0.03384437400744206_dp*logrh2)/x           &
            + 0.04743226764572505_dp*logrh3               &
            - 0.0006251042204583412_dp*t*logrh3           &
            + 2.650663328519478e-6_dp*t2*logrh3           &
            - 3.674710848763778e-9_dp*t3*logrh3           &
            - (0.0002672510825259393_dp*logrh3)/x         &
            - 0.01252108546759328_dp*logrhoa              &
            + 0.005806550506277202_dp*t*logrhoa           &
            - 0.0001016735312443444_dp*t2*logrhoa         &
            + 2.881946187214505e-7_dp*t3*logrhoa          &
            + (0.0942243379396279_dp*logrhoa)/x           &
            - 0.0385459592773097_dp*logrh*logrhoa         &
            - 0.0006723156277391984_dp*t*logrh*logrhoa    &
            + 2.602884877659698e-6_dp*t2*logrh*logrhoa    &
            + 1.194163699688297e-8_dp*t3*logrh*logrhoa    &
            - (0.00851515345806281_dp*logrh*logrhoa)/x    &
            - 0.01837488495738111_dp*logrh2*logrhoa       &
            + 0.0001720723574407498_dp*t*logrh2*logrhoa   &
            - 3.717657974086814e-7_dp*t2*logrh2*logrhoa   &
            - 5.148746022615196e-10_dp*t3*logrh2*logrhoa  &
            + (0.0002686602132926594_dp*logrh2*logrhoa)/x &
            - 0.06199739728812199_dp*logrhoa2             &
            + 0.000906958053583576_dp*t*logrhoa2          &
            - 9.11727926129757e-7_dp*t2*logrhoa2          &
            - 5.367963396508457e-9_dp*t3*logrhoa2         &
            - (0.007742343393937707_dp*logrhoa2)/x        &
            + 0.0121827103101659_dp*logrh*logrhoa2        &
            - 0.0001066499571188091_dp*t*logrh*logrhoa2   &
            + 2.534598655067518e-7_dp*t2*logrh*logrhoa2   &
            - 3.635186504599571e-10_dp*t3*logrh*logrhoa2  &
            + (0.0006100650851863252_dp*logrh*logrhoa2)/x &
            + 0.0003201836700403512_dp*logrhoa3           &
            - 0.0000174761713262546_dp*t*logrhoa3         &
            + 6.065037668052182e-8_dp*t2*logrhoa3         &
            - 1.421771723004557e-11_dp*t3*logrhoa3        &
            + (0.0001357509859501723_dp*logrhoa3)/x

       ntot = exp(ntot)

       mdot1(ncell) = 0.0_dp
       ndot1(ncell) = 0.0_dp

       IF (ntot >= 4.0_dp) THEN

          ! number of sulfuric acid molecules in the cluster
          nac = x * ntot
          ndot1(ncell) = exp(jnuc) * 1.0e6_dp
          ! number of micrograms of sulfate in a 3.5 nm particle at ambient RH
          mp = molec2ug * mp35(rh)

          ! assume log-normal distribution with d=3.5 nm,
          ! sigma = sigma(Aitken mode)
          mdot1(ncell) = ndot1(ncell) * mp * es36(akn)

          IF (mdot1(ncell) > so4rat(ncell)) THEN
             ! limit nucleated mass by available mass
             mdot1(ncell) = so4rat(ncell)
             ! adjust ndot1 to this mass
             ndot1(ncell) = mdot1(ncell) / (mp * es36(akn))
          ENDIF

          IF (mdot1(ncell) == 0.0_dp) ndot1(ncell) = 0.0_dp

       END IF ! if (ntot >= 4.0_dp)

       ! Rescale number of newly nucleated particles to account
       ! for non-linearities in subgrid scale aging of the nucleated
       ! particle population
       IF (rset_nucsize%l) THEN 
          ndot1(ncell) = ndot1(ncell) * 6._dp / pi * mp / &
               (rho(i_so4)  * 1.0e9_dp) / &
               (rset_nucsize%v * 1.0e-9_dp)**3
       END IF

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

       ! *** Add new particles to 3rd moment growth rate of soluble Aitken mode
       cgr3(akn,ncell) = cgr3(akn,ncell) + MASS2MOM3(i_so4) * mdot1(ncell)

    END DO

  END SUBROUTINE made3_nucleate

!-------------------------------------------------------------------------------

  !> \brief Integrates number and mass equations for each mode over the time
  !>   step
  !> \details Called after the microphysics routines, this subroutine performs
  !>   the actual "time stepping", i.e., it integrates the (rest of the)
  !>   aerosol dynamics equation (without gas-particle partitioning, because
  !>   that was already treated) over one time step using the rates of change
  !>   supplied by the microphysics routines. It first sets up the target modes
  !>   for the coagulation processes, then integrates number and mass
  !>   concentration equations, updates gas phase sulfuric acid concentration,
  !>   and reduces the MADE-internal time step if necessary.

  SUBROUTINE MADE3_AEROSTEP ( NUMCELLS, CBLK, DT                   &
       , FCONC, FCONC_ORG, DMDT, DNDT, DELTASO4A, SOA_MADE3, CR0, C30, CGR3 &
!sorgam                    ORGARO1RAT, ORGARO2RAT,                  &
!sorgam                    ORGALK1RAT, ORGOLE1RAT,                  &
!sorgam                    ORGBIO1RAT, ORGBIO2RAT,                  &
!sorgam                    ORGBIO3RAT, ORGBIO4RAT,                  &
!unused                    EPM25I, EPM25J, EORGI,EORGJ, EECI, EECJ, &
!unused                    ESOIL, ESEAS, EPMCOARSE,                 &
!unused                    PMASSN, PMASSA, PMASSC,                  &
       )

    IMPLICIT NONE
    INTRINSIC EXP, SUM, SQRT, MAX, MIN, SIZE

    ! *** ARGUMENTS ***
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)    :: NUMCELLS
    !> Tracer array
    REAL(dp), INTENT(inout) :: CBLK(dim1_cblk,dim2_cblk,NUMCELLS)
! op_ck_20120309+
    !> Time step [s]
!    REAL(dp), INTENT(in)    :: DT
    REAL(dp), INTENT(inout) :: DT
! op_ck_20120309-
    !> H<sub>2</sub>SO<sub>4</sub> condensation rate coefficients
    !> [s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: FCONC(nmod,NUMCELLS)
    !> SOA precursor condensation rate coefficients [s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: FCONC_ORG(nmod,NUMCELLS)
    !> Nucleation rate (mass SO<sub>4</sub>) [ug m<sup>-3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: DMDT(NUMCELLS)
    !> Nucleation rate (number) [m<sup>-3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: DNDT(NUMCELLS)
    !> Condensing mass of SO<sub>4</sub> [ug m<sup>-3</sup>]
    REAL(dp), INTENT(in)    :: DELTASO4A(NUMCELLS)
    !> SOA precursor emissions [ug m<sup>-3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: SOA_MADE3(NUMCELLS)
    !> Modal number coagulation rate coefficients [m<sup>3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(inout)    :: CR0(nmod,nmod,NUMCELLS)
    !> Intermodal 3<sup>rd</sup> moment transfer rates [mom_3 m<sup>-3</sup>
    !> s<sup>-1</sup>]
    REAL(dp), INTENT(inout)    :: C30(nmod,nmod, NUMCELLS)
    !> Modal 3<sup>rd</sup> moment growth rates [mom_3 m<sup>-3</sup>
    !> s<sup>-1</sup>]
    REAL(dp), INTENT(inout) :: CGR3(nmod,NUMCELLS)
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     aromatics
!sorgam   REAL(dp) :: ORGARO1RAT(NUMCELLS)
!sorgam   REAL(dp) :: ORGARO2RAT(NUMCELLS)
!sorgam   ! *** anthropogenic organic aerosol mass production rates from
!sorgam   !     alkanes & others
!sorgam   REAL(dp) :: ORGALK1RAT(NUMCELLS)
!sorgam   REAL(dp) :: ORGOLE1RAT(NUMCELLS)
!sorgam   ! *** biogenic organic aerosol production rates
!sorgam   REAL(dp) :: ORGBIO1RAT(NUMCELLS)
!sorgam   REAL(dp) :: ORGBIO2RAT(NUMCELLS)
!sorgam   REAL(dp) :: ORGBIO3RAT(NUMCELLS)
!sorgam   REAL(dp) :: ORGBIO4RAT(NUMCELLS)

    ! *** Local Variables ***
    INTEGER :: L, jm, jm2, js   ! Loop indices
    ! *** variables needed for modal dynamics solvers:
    REAL(dp), PARAMETER :: irrelevant = 1.0001_dp
    REAL(dp) :: A, B(nmod), C(nmod)
    INTEGER  :: i_tgt
    REAL(dp) :: M1, M2, Y0, Y
    REAL(dp) :: DHAT, PEXPDT, EXPDT
    REAL(dp) :: PROD(nspec), POL(nspec), LOSSINV
    REAL(dp) :: LOSS(nmod), gainmom3(nmod,nmod)
    INTEGER  :: tgtmode(nmod,nmod)
    REAL(dp) :: MSTRNSFR(nspec)   ! intermodal mass transfer by coagulation
    REAL(dp) :: FACTRANS          ! special factor to compute mass transfer
    !va internal variables added for introduction of new modes
    REAL(dp) :: GAIN
    REAL(dp) :: OLD(nmod,nspec+1)
! op_ck_20120309+
    ! variables for adaptive time step additions
    INTEGER  :: chkmode(2)            ! modes to check for
                                      ! ``numerical particle generation''
    INTEGER  :: jcm                   ! loop index
    REAL(dp) :: deltan, sup, supinv   ! helper variables
    LOGICAL  :: numgen                ! particles were ``generated numerically''
! op_ck_20120309-


    ! -------------------------------- Start code ------------------------------

    ! Initializations
    chkmode  = (/accs,cors/)
    numgen   = .FALSE.

    ! *** set up time-step integration

    loop_cells: DO L = 1, NUMCELLS

       B       = 0._dp
       C       = 0._dp
       tgtmode = 0

       ! save current species and # concentrations
       DO jm = 1, nmod
          DO js = 1, nspec+1
             OLD(jm,js) = CBLK(js,jm,L)
          END DO
       END DO

       ! tgtmode will hold the target modes for all possible coagulation events.
       ! The target modes are set according to the following table. Some
       ! particles may be assigned either to a mixed or to an insoluble mode. In
       ! these cases, the subroutine target_mode is used to find the correct
       ! target mode based on a threshold mass fraction of soluble material.
       
! |       | akn | akns | sooti | acc   | accs  | sootj | cor   | cors  | sootc |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | akn   | akn | akns | akns/ | acc   | accs  | accs/ | cor   | cors  | sootc |
! |       |     |      | sooti |       |       | sootj |       |       |       |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | akns  |     | akns | akns/ | accs  | accs  | accs/ | cors  | cors  | sootc |
! |       |     |      | sooti |       |       | sootj |       |       |       |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | sooti |     |      | sooti | accs/ | accs/ | sootj | cors  | cors  | sootc |
! |       |     |      |       | sooti | sooti |       |       |       |       |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | acc   |     |      |       | acc   | accs  | accs/ | cor   | cors  | cors/ |
! |       |     |      |       |       |       | sootj |       |       | sootc |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | accs  |     |      |       |       | accs  | accs/ | cors  | cors  | cors/ |
! |       |     |      |       |       |       | sootj |       |       | sootc |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | sootj |     |      |       |       |       | sootj | cors/ | cors/ | sootc |
! |       |     |      |       |       |       |       | sootj | sootj |       |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | cor   |     |      |       |       |       |       | cor   | cors  | cors/ |
! |       |     |      |       |       |       |       |       |       | sootc |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | cors  |     |      |       |       |       |       |       | cors  | cors/ |
! |       |     |      |       |       |       |       |       |       | sootc |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|
! | sootc |     |      |       |       |       |       |       |       | sootc |
! |-------+-----+------+-------+-------+-------+-------+-------+-------+-------|

       tgtmode(akns,akn)    = akns
       tgtmode(acc,akn)     = acc
       tgtmode(accs,akn)    = accs
       tgtmode(cor,akn)     = cor
       tgtmode(cors,akn)    = cors
       tgtmode(sootc,akn)   = sootc
       tgtmode(acc,akns)    = accs
       tgtmode(accs,akns)   = accs
       tgtmode(cor,akns)    = cors
       tgtmode(cors,akns)   = cors
       tgtmode(sootc,akns)  = sootc
       tgtmode(sootj,sooti) = sootj
       tgtmode(cor,sooti)   = cors
       tgtmode(cors,sooti)  = cors
       tgtmode(sootc,sooti) = sootc
       tgtmode(accs,acc)    = accs
       tgtmode(cor,acc)     = cor
       tgtmode(cors,acc)    = cors
       tgtmode(cor,accs)    = cors
       tgtmode(cors,accs)   = cors
       tgtmode(sootc,sootj) = sootc
       tgtmode(cors,cor)    = cors

       ! LOSS = Normalized coagulation transfer rates
       ! Note that LOSS here does not necessarily mean that 3rd moment is
       ! removed from the mode. Some of the components may actually be 0 due to
       ! assignment of the transferred 3rd moment (C30(X,Y,L)) to the mode it
       ! originated from.
       DO jm = 1, nmod
          LOSS(jm) = SUM(C30(jm,:,L)) / CBLK(i_mom3,jm,L)
       END DO

       ! sooti + akn -> akns OR sooti
       CALL target_mode((/akns,sooti/), CBLK(:,(/akn,sooti/),L)           &
            , LOSS((/akn,sooti/)), DT, C30(akn,sooti,L), C30(sooti,akn,L) &
            , EPSILON_MAX, tgtmode(akn,sooti))

       ! sootj + akn -> accs OR sootj
       CALL target_mode((/accs,sootj/), CBLK(:,(/akn,sootj/),L)           &
            , LOSS((/akn,sootj/)), DT, C30(akn,sootj,L), C30(sootj,akn,L) &
            , EPSILON_MAX, tgtmode(akn,sootj))

       ! sooti + akns -> akns OR sooti
       CALL target_mode((/akns,sooti/), CBLK(:,(/akns,sooti/),L)             &
            , LOSS((/akns,sooti/)), DT, C30(akns,sooti,L), C30(sooti,akns,L) &
            , EPSILON_MAX, tgtmode(akns,sooti))

       ! sootj + akns -> accs OR sootj
       CALL target_mode((/accs,sootj/), CBLK(:,(/akns,sootj/),L)             &
            , LOSS((/akns,sootj/)), DT, C30(akns,sootj,L), C30(sootj,akns,L) &
            , EPSILON_MAX, tgtmode(akns,sootj))

       ! acc + sooti -> accs OR sooti
       CALL target_mode((/accs,sooti/), CBLK(:,(/acc,sooti/),L)           &
            , LOSS((/acc,sooti/)), DT, C30(acc,sooti,L), C30(sooti,acc,L) &
            , EPSILON_MAX, tgtmode(acc,sooti))

       ! accs + sooti -> accs OR sooti
       CALL target_mode((/accs,sooti/), CBLK(:,(/accs,sooti/),L)             &
            , LOSS((/accs,sooti/)), DT, C30(accs,sooti,L), C30(sooti,accs,L) &
            , EPSILON_MAX, tgtmode(accs,sooti))

       ! sootj + acc -> accs OR sootj
       CALL target_mode((/accs,sootj/), CBLK(:,(/acc,sootj/),L)           &
            , LOSS((/acc,sootj/)), DT, C30(acc,sootj,L), C30(sootj,acc,L) &
            , EPSILON_MAX, tgtmode(acc,sootj))

       ! sootc + acc -> cors OR sootc
       CALL target_mode((/cors,sootc/), CBLK(:,(/acc,sootc/),L)           &
            , LOSS((/acc,sootc/)), DT, C30(acc,sootc,L), C30(sootc,acc,L) &
            , EPSILON_MAX, tgtmode(acc,sootc))

       ! sootj + accs -> accs OR sootj
       CALL target_mode((/accs,sootj/), CBLK(:,(/accs,sootj/),L)             &
            , LOSS((/accs,sootj/)), DT, C30(accs,sootj,L), C30(sootj,accs,L) &
            , EPSILON_MAX, tgtmode(accs,sootj))

       ! sootc + accs -> cors OR sootc
       CALL target_mode((/cors,sootc/), CBLK(:,(/accs,sootc/),L)             &
            , LOSS((/accs,sootc/)), DT, C30(accs,sootc,L), C30(sootc,accs,L) &
            , EPSILON_MAX, tgtmode(accs,sootc))

       ! cor + sootj -> cors OR sootj
       CALL target_mode((/cors,sootj/), CBLK(:,(/cor,sootj/),L)           &
            , LOSS((/cor,sootj/)), DT, C30(cor,sootj,L), C30(sootj,cor,L) &
            , EPSILON_MAX, tgtmode(cor,sootj))

       ! cors + sootj -> cors OR sootj
       CALL target_mode((/cors,sootj/), CBLK(:,(/cors,sootj/),L)             &
            , LOSS((/cors,sootj/)), DT, C30(cors,sootj,L), C30(sootj,cors,L) &
            , EPSILON_MAX, tgtmode(cors,sootj))

       ! sootc + cor -> cors OR sootc
       CALL target_mode((/cors,sootc/), CBLK(:,(/cor,sootc/),L)           &
            , LOSS((/cor,sootc/)), DT, C30(cor,sootc,L), C30(sootc,cor,L) &
            , EPSILON_MAX, tgtmode(cor,sootc))

       ! sootc + cors -> cors OR sootc
       CALL target_mode((/cors,sootc/), CBLK(:,(/cors,sootc/),L)             &
            , LOSS((/cors,sootc/)), DT, C30(cors,sootc,L), C30(sootc,cors,L) &
            , EPSILON_MAX, tgtmode(cors,sootc))

       ! Ensure symmetry of array tgtmode
       DO jm2 = 1, nmod - 1
          DO jm = jm2 + 1, nmod
             IF (tgtmode(jm,jm2) .EQ. 0) THEN
                tgtmode(jm,jm2) = tgtmode(jm2,jm)
             ELSE
                tgtmode(jm2,jm) = tgtmode(jm,jm2)
             END IF
          END DO
       END DO
! DEBUG+
!!$       write(*,*) 'tgtmode =', tgtmode ! -> It actually IS 0 along the diag.
! DEBUG-

       !***************** UPDATE NUMBER CONCENTRATION ************************
       ! *** code to move number forward by one time step.
       ! *** solves the Riccati equation:
       !     dY/dt = C - A * Y ** 2 - B * Y

       ! *** set coefficients B and C
       C(akn) = DNDT(L)   ! extra term for soluble Aitken mode (nucleation)
       DO jm = 1, nmod
          DO jm2 = 1, nmod
             i_tgt = tgtmode(jm2,jm)
             IF ((i_tgt .NE. 0) .AND. (i_tgt .NE. jm) .AND. &
                  (CBLK(i_num,jm2,L)/nummin(jm2) .GT. irrelevant)) THEN
                B(jm) = B(jm) + CR0(jm,jm2,L) * CBLK(i_num,jm2,L)
                IF ((jm .LT. jm2) .AND. (i_tgt .NE. jm2) .AND. &
                     (CBLK(i_num,jm,L)/nummin(jm) .GT. irrelevant)) &
                     ! Limit production of particles in a 3rd mode in case one
                     ! of the source modes contains very few particles
                     C(i_tgt) = C(i_tgt) + MIN(&
                       CR0(jm2,jm,L) * CBLK(i_num,jm2,L) * CBLK(i_num,jm,L), &
                       MIN(CBLK(i_num,jm2,L), CBLK(i_num,jm,L)) / DT)
             END IF
          END DO
       END DO

       ! *** set coefficients A and solve Riccati equation
       DO jm = 1, nmod

          Y = 0._dp

          IF (CBLK(i_num,jm,L)/nummin(jm) .GT. irrelevant) THEN
             Y0 = CBLK(i_num,jm,L)
          ELSE
             Y0 = 0._dp
          END IF
          A = CR0(jm,jm,L)

! op_ck_20121001  FIXME: Which is the best minimum value to compare against
!                        instead of 0.0? See also check of B. Possibly use
!                        1.e-30_dp!? Or epsilon(1.0_dp)? Or something relative
!                        to B and C?
! op_ck_20140612         -> This should be tested in the future.
          IF (C(jm) > 1.e-19_dp) THEN   ! use full solution

             DHAT = SQRT(B(jm) * B(jm) + 4.0_dp * A * C(jm))
             M1 = 2.0_dp * A * C(jm) / (B(jm) + DHAT)
             M2 = -0.5_dp * (B(jm) + DHAT)
             PEXPDT = -((M1 - A  * Y0) / (M2 - A * Y0)) * EXP(-DHAT * DT)

             Y = (M1 + M2 * PEXPDT) / (A * (1.0_dp + PEXPDT))

          ELSE   ! use simplified solution

             IF (B(jm) > 0.0_dp) THEN
                PEXPDT = EXP(-B(jm) * DT)
                Y = B(jm) * Y0 * PEXPDT &
                     / (B(jm) + A * Y0 * (1.0_dp - PEXPDT))
             ELSE   ! use simple solution
                Y = Y0 / (1.0_dp + A * Y0 * DT)
             END IF

          END IF

          CBLK(i_num,jm,L) = MAX(nummin(jm),Y)

       END DO

! op_ck_20120309+
       ! Adaptive time step method:
       ! In order to avoid inconsistencies in transfer of coagulating particles
       ! between modes, it has to be checked whether particles are artificially
       ! ``created''. If this is the case, reduce the time step.
       
       DO jcm = 1, SIZE(chkmode)
          ! actual increase in number concentration
          deltan = CBLK(i_num,chkmode(jcm),L) - OLD(chkmode(jcm),i_num)
          IF (deltan / OLD(chkmode(jcm),i_num) .LE. irrelevant) CYCLE
          IF (deltan / (C(chkmode(jcm)) * DT) .GT. irrelevant) THEN
             numgen = .TRUE.
          END IF
       END DO

       IF (numgen) THEN
          WRITE(*,*) 'MADE3 INFO: Halving internal time step!'
          DT = 0.5_dp * DT
          RETURN
       END IF
! op_ck_20120309-


       !***************** UPDATE MASS CONCENTRATION **************************

       ! *** Prepare to advance modal mass concentration one time step.
       !     (Emissions are handled by MESSy submodels OFFEMIS and ONEMIS.)
       !     For each mode X:
       !     1. ... set normalized 3rd moment coag. loss rate: LOSS(X) = ...
       !     2. ... set 3rd moment coagulation gain rates from other modes Y:
       !        gainmom3(Y,X) = ...
       !     Note that C30(X,X) = 0.0_dp for all X.
       DO jm2 = 1, nmod
          LOSS(jm2) = SUM(C30(jm2,:,L), mask=tgtmode(jm2,:).NE.jm2) &
               / CBLK(i_mom3,jm2,L)
          DO jm = 1, nmod
             gainmom3(jm,jm2) = SUM(C30(jm,:,L), mask=tgtmode(jm,:).EQ.jm2)
          END DO
          gainmom3(jm2,jm2) = 0._dp
       END DO

       ! *** set up special factors for mass transfer between modes by
       !     intermodal coagulation. The mass transfer rate is proportional to
       !     the 3rd moment transfer rate, C30. The proportionality factor is
       !     p/6 times the the average particle density. The average particle
       !     density for a species is the species mass concentration divided by
       !     the particle volume concentration, pi/6 times the 3rd moment
       !     concentration. The pi/6 coefficients cancel.
       DO jm2 = 1, nmod

          MSTRNSFR = 0._dp
          PROD     = 0._dp

          DO jm = 1, nmod
             IF (LOSS(jm) .GT. 1.e-30_dp) THEN
                GAIN = gainmom3(jm,jm2) / CBLK(i_mom3,jm,L)
                ! FACTRANS is the fraction of the loss from mode jm that is
                ! transferred to mode jm2.
! op_ck_20130904  FIXME: Why do we assume here that no 3rd moment is produced in
!                        mode jm?
! op_ck_20140612         -> Need to assume ``GAIN'' as constant for the ODE
!                           solution to be applicable.
!                           Should be tested with the box model (only insol.
!                           particles, no gas phase, short dt vs. long dt).
                FACTRANS = GAIN / LOSS(jm) * (1.0_dp - EXP(-LOSS(jm) * DT))
             ELSE
                FACTRANS = 0._dp
             END IF
             MSTRNSFR = MSTRNSFR + FACTRANS * OLD(jm,1:nspec)
          END DO

          IF (LOSS(jm2) .GE. 1.e-30_dp) THEN

             PROD(i_so4) = DELTASO4A(L) * FCONC(jm2,L) / DT
             PROD(i_pom) = SOA_MADE3(L) * FCONC_ORG(jm2,L)
             IF (jm2 == akn) THEN
                PROD(i_so4) = PROD(i_so4) + DMDT(L)
             ELSE
                PROD = PROD + MSTRNSFR / DT
             END IF
             LOSSINV = 1.0_dp / LOSS(jm2)
             POL = PROD * LOSSINV
             EXPDT = EXP(-LOSS(jm2) * DT)

             CBLK(1:nspec,jm2,L) = MAX(POL(:) &
                  + (CBLK(1:nspec,jm2,L) - POL(:)) * EXPDT,MASSMIN(jm2,1:nspec))

          ELSE

             PROD(i_so4) = DELTASO4A(L) * FCONC(jm2,L)
             PROD(i_pom) = SOA_MADE3(L) * FCONC_ORG(jm2,L) * DT
             IF (jm2 == akn) THEN
                PROD(i_so4) = PROD(i_so4) + DMDT(L) * DT
             ELSE
                PROD = PROD + MSTRNSFR
             END IF

             CBLK(1:nspec,jm2,L) = &
                  MAX(CBLK(1:nspec,jm2,L) + PROD(:), MASSMIN(jm2,1:nspec))

          END IF

       END DO

       ! *** Update sulfuric acid vapor concentration by removing mass
       !     concentration of condensed sulfate and newly produced particles.
       ! *** The method follows Youngblood and Kreidenweis, Further Development
       !     and Testing of a Bimodal Aerosol Dynamics Model, Colorado State
       !     University Department of Atmospheric Science Paper Number 550,
       !     April,1994, pp 85-89.
       
       CBLK(i_h2so4,gas,L) = MAX(CONMIN, CBLK(i_h2so4,gas,L) &
            - (DELTASO4A(L) + DMDT(L) * DT)                  &
            * MW(i_h2so4,i_mwgas)/MW(i_so4,i_mwaero))

!sorgam      ! *** anthropogenic secondary organic:
!sorgam      !bs * anthropogenic secondary organics from aromatic precursors
!sorgam
!sorgam      MSTRNSFR = CBLK(VORGARO1I,L) * FACTRANS
!sorgam      PROD = ORGARO1RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSmodeINV
!sorgam
!sorgam      CBLK(VORGARO1I,L) = POL + (CBLK(VORGARO1I,L) - POL) * EXPDT
!sorgam
!sorgam      CBLK(VORGARO1I,L) = MAX(CONMIN, CBLK(VORGARO1I,L))
!sorgam
!sorgam      CBLK(VORGARO1J,L) = CBLK(VORGARO1J,L)                    &
!sorgam                          + ORGARO1RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                          + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(VORGARO1J,L) = MAX(CBLK(VORGARO1J,L),1.e-30_dp)
!sorgam
!sorgam      !bs * second species from aromatics
!sorgam      MSTRNSFR = CBLK(VORGARO2I,L) * FACTRANS
!sorgam      PROD = ORGARO2RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSmodeINV
!sorgam
!sorgam      CBLK(VORGARO2I,L) = POL + (CBLK(VORGARO2I,L) - POL) * EXPDT
!sorgam
!sorgam      CBLK(VORGARO2I,L) = MAX(CONMIN, CBLK(VORGARO2I,L))
!sorgam
!sorgam      CBLK(VORGARO2J,L) = CBLK(VORGARO2J,L)                    &
!sorgam                          + ORGARO2RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                          + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(VORGARO2J,L) = MAX(CBLK(VORGARO2J,L),1.e-30_dp)
!sorgam
!sorgam      !bs * anthropogenic secondary organics from alkanes & other
!sorgam      !bs   precursors
!sorgam      !bs * higher alkanes
!sorgam      MSTRNSFR = CBLK(VORGALK1I,L) * FACTRANS
!sorgam      PROD = ORGALK1RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSmodeINV
!sorgam
!sorgam      CBLK(VORGALK1I,L) = POL + (CBLK(VORGALK1I,L) - POL) * EXPDT
!sorgam
!sorgam      CBLK(VORGALK1I,L) = MAX(CONMIN, CBLK(VORGALK1I,L))
!sorgam
!sorgam      CBLK(VORGALK1J,L) = CBLK(VORGALK1J,L)                    &
!sorgam                          + ORGALK1RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                          + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(VORGALK1J,L) = MAX(CBLK(VORGALK1J,L),1.e-30_dp)
!sorgam
!sorgam      !bs * higher olefines
!sorgam      MSTRNSFR = CBLK(VORGOLE1I,L) * FACTRANS
!sorgam      PROD = ORGOLE1RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSmodeINV
!sorgam
!sorgam      CBLK(VORGOLE1I,L) = POL + (CBLK(VORGOLE1I,L) - POL) * EXPDT
!sorgam
!sorgam      CBLK(VORGOLE1I,L) = MAX(CONMIN, CBLK(VORGOLE1I,L))
!sorgam
!sorgam      CBLK(VORGOLE1J,L) = CBLK(VORGOLE1J,L)                    &
!sorgam                          + ORGOLE1RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                          + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(VORGOLE1J,L) = MAX(CBLK(VORGOLE1J,L),1.e-30_dp)
!sorgam
!sorgam      ! *** biogenic secondary organic
!sorgam
!sorgam      MSTRNSFR = CBLK(VORGBA1I,L) * FACTRANS
!sorgam      PROD = ORGBIO1RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSmodeINV
!sorgam
!sorgam      CBLK(VORGBA1I,L) = POL + (CBLK(VORGBA1I,L) - POL) * EXPDT
!sorgam
!sorgam      CBLK(VORGBA1I,L) = MAX(CONMIN, CBLK(VORGBA1I,L))
!sorgam
!sorgam      CBLK(VORGBA1J,L) = CBLK(VORGBA1J,L)                     &
!sorgam                         + ORGBIO1RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                         + MSTRNSFR
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(VORGBA1J,L) = MAX(CBLK(VORGBA1J,L),1.e-30_dp)
!sorgam
!sorgam      !bs * second biogenic species
!sorgam      MSTRNSFR = CBLK(VORGBA2I,L) * FACTRANS
!sorgam      PROD = ORGBIO2RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSmodeINV
!sorgam
!sorgam      CBLK(VORGBA2I,L) = POL + (CBLK(VORGBA2I,L) - POL) * EXPDT
!sorgam
!sorgam      CBLK(VORGBA2I,L) = MAX(CONMIN, CBLK(VORGBA2I,L))
!sorgam
!sorgam      CBLK(VORGBA2J,L) = CBLK(VORGBA2J,L)                     &
!sorgam                         + ORGBIO2RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                         + MSTRNSFR
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(VORGBA2J,L) = MAX(CBLK(VORGBA2J,L),1.e-30_dp)
!sorgam
!sorgam      !bs * third biogenic species
!sorgam      MSTRNSFR = CBLK(VORGBA3I,L) * FACTRANS
!sorgam      PROD = ORGBIO3RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSmodeINV
!sorgam
!sorgam      CBLK(VORGBA3I,L) = POL + (CBLK(VORGBA3I,L) - POL) * EXPDT
!sorgam
!sorgam      CBLK(VORGBA3I,L) = MAX(CONMIN, CBLK(VORGBA3I,L))
!sorgam
!sorgam      CBLK(VORGBA3J,L) = CBLK(VORGBA3J,L)                     &
!sorgam                         + ORGBIO3RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                         + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(VORGBA3J,L) = MAX(CBLK(VORGBA3J,L),1.e-30_dp)
!sorgam
!sorgam      !bs * fourth biogenic species
!sorgam      MSTRNSFR = CBLK(VORGBA4I,L) * FACTRANS
!sorgam      PROD = ORGBIO4RAT(L) * FCONCN_ORG(L)
!sorgam      POL = PROD * LOSSmodeINV
!sorgam
!sorgam      CBLK(VORGBA4I,L) = POL + (CBLK(VORGBA4I,L) - POL) * EXPDT
!sorgam
!sorgam      CBLK(VORGBA4I,L) = MAX(CONMIN, CBLK(VORGBA4I,L))
!sorgam
!sorgam      CBLK(VORGBA4J,L) = CBLK(VORGBA4J,L)                     &
!sorgam                         + ORGBIO4RAT(L) * FCONCA_ORG(L) * DT &
!sorgam                         + MSTRNSFR
!sorgam
!sorgam      !iatest quick fit!!!
!sorgam      CBLK(VORGBA4J,L) = MAX(CBLK(VORGBA4J,L),1.e-30_dp)

       ! *** Add growth from transfer of 3rd moment by intermodal coagulation
       !     to 3rd moment growth rates (required for renaming criteria).
       DO jm = 1, nmod
          !va note that if I have emission I will have to include them here
          CGR3(jm,L) = CGR3(jm,L) + SUM(C30(:,:,L), mask=tgtmode(:,:).EQ.jm) &
               - SUM(C30(jm,:,L), mask=tgtmode(jm,:).EQ.jm)
       END DO

    END DO loop_cells

  END SUBROUTINE MADE3_AEROSTEP

!-------------------------------------------------------------------------------

  !> \brief Computes the target mode and transfer rates (from one mode to the
  !>   the other) for a coagulation process
  !> \details Determines the target mode for the coagulation of soluble or mixed
  !>   particles with insoluble ones and sets the number and 3<sup>rd</sup>
  !>   moment transfer rate(s) (coefficients) correspondingly.
  !> \note This routine assumes that the second mode in the input arrays
  !>   \e mode, \e conc, \e loss, and \e mom3trnsfr is the insoluble one.
  SUBROUTINE target_mode(mode, conc, loss, dt, mom3trnsfr1, mom3trnsfr2 &
       , maxsolfrac, tgt)
  
    IMPLICIT NONE
    INTRINSIC :: EXP, SUM

    ! I/O
    !> Mode indices of possible target modes
    INTEGER,  INTENT(in)    :: mode(2)
    !> Tracer array (subset of \c CBLK as used in other routines)
    REAL(dp), INTENT(in)    :: conc(dim1_cblk,2)
    !> Total relative 3<sup>rd</sup> moment loss rate coefficients from source
    !> modes [s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: loss(2)
    !> Time step [s]
    REAL(dp), INTENT(in)    :: dt
    !> 3<sup>rd</sup> mom. transfer rate for soluble or mixed source mode [mom_3
    !> m-3 s-1]
    REAL(dp), INTENT(in)    :: mom3trnsfr1
    !> 3<sup>rd</sup> mom. transfer rate for insoluble source mode [mom_3
    !> m<sup>-3</sup> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: mom3trnsfr2
    !> Maximum allowed fraction of soluble mass in the insoluble mode
    REAL(dp), INTENT(in)    :: maxsolfrac
    !> Target mode index
    INTEGER,  INTENT(out)   :: tgt

    ! LOCAL
    INTEGER  :: jm                 ! loop index
    REAL(dp) :: factrans           ! 3rd moment transfer factor
    REAL(dp) :: mom3trnsfr_cp(2)   ! copy of 3rd moment transfer rates
    REAL(dp) :: soltrnsfr          ! transferred mass conc. of soluble species
    REAL(dp) :: drytrnsfr          ! transferred total dry mass conc.


    mom3trnsfr_cp(1) = mom3trnsfr1
    mom3trnsfr_cp(2) = mom3trnsfr2
    soltrnsfr        = 0._dp
    drytrnsfr        = 0._dp

    DO jm = 1, 2

       IF (loss(jm) .GT. 1.e-30_dp) THEN
          ! relative rate of 3rd moment transfer to new particles:
          !   gain = mom3trnsfr(jm) / conc(i_mom3,jm)
          ! total fraction of 3rd moment that is lost from source mode jm:
          !   mom3fraclost = (1.0_dp - EXP(-loss(jm) * dt)
          ! fraction of 3rd moment that is lost from source mode jm via the
          ! coagulation process for which the target mode is computed:
          !   factrans = (gain / loss(jm)) * mom3fraclost
          factrans = mom3trnsfr_cp(jm) * (1.0_dp - EXP(-loss(jm) * dt)) &
               / (loss(jm) * conc(i_mom3,jm))
       ELSE
          factrans = 0._dp
       END IF

       soltrnsfr = soltrnsfr + factrans * SUM(conc(SOLSPEC,jm))
       drytrnsfr = drytrnsfr + factrans * (SUM(conc(1:nspec,jm))-conc(i_h2o,jm))

    END DO

    IF (soltrnsfr .GE. maxsolfrac * drytrnsfr) THEN
       ! transfer to 'soluble' mode
       tgt = mode(1)
    ELSE
       ! transfer to insoluble mode
       tgt = mode(2)
    END IF

  END SUBROUTINE target_mode

!-------------------------------------------------------------------------------

  !> \brief Moves particles from one mode to another if necessary
  !> \details  Renames particles from insoluble modes to mixed modes and from
  !>   Aitken to accumulation modes if necessary.
  !> \note Criteria for renaming are given in the "Version" section of this
  !>   documentation

  SUBROUTINE made3_rename(status, numcells, cblk, dg, cgr3 &
       , mergenum, mergem3, mergenums, mergem3s, mergenumsoot, mergem3soot)

    IMPLICIT NONE
    INTRINSIC SUM

    ! I/O
    !> Error status flag
    INTEGER,  INTENT(out)   :: status
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)    :: numcells
    !> Tracer array
    REAL(dp), INTENT(inout) :: cblk(dim1_cblk,dim2_cblk,numcells)
    !> Wet median modal diameters [m]
    REAL(dp), INTENT(in)    :: dg(nmod,numcells)
    !> Modal 3<sup>rd</sup> moment growth rates [mom_3 m<sup>-3</sup>
    !> s<sup>-1</sup>]
    REAL(dp), INTENT(in)    :: cgr3(nmod,numcells)
    !> Number concentration of particles renamed from akn to acc
    !> [m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: mergenum(numcells)
    !> 3<sup>rd</sup> moment concentration of particles renamed from akn to acc
    !> [mom_3 m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: mergem3(numcells)
    !> Number concentration of particles renamed from akns to accs
    !> [m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: mergenums(numcells)
    !> 3<sup>rd</sup> moment conc. of particles renamed from akns to accs [mom_3
    !> m-3]
    REAL(dp), INTENT(out)   :: mergem3s(numcells)
    !> Number concentration of particles renamed from sooti to sootj
    !> [m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: mergem3soot(numcells)
    !> 3<sup>rd</sup> moment conc. of particles renamed from sooti to sootj
    !> [mom_3 m<sup>-3</sup>]
    REAL(dp), INTENT(out)   :: mergenumsoot(numcells)

    ! LOCAL
    INTEGER, PARAMETER  :: allexcssdu(nspec-3) = &   ! all indices except SS, DU
         (/i_so4,i_nh4,i_no3,i_pom,i_bc,i_bctag,i_h2o/) ! op_mr_20181002
    REAL(dp)            :: xxm3(nmod), fnum, fm3 ! defined below
    REAL(dp)            :: epsilon               ! soluble fraction in the
                                                 ! insoluble mode
    INTEGER             :: mgmode(2)             ! mode indices for merging
    INTEGER             :: jp                    ! loop index
    LOGICAL             :: firstime
    DATA                   firstime / .TRUE. /
    SAVE                   firstime, xxm3


    ! -------------------------------- Start code ------------------------------

    ! Conversion ``factors'' for error function argument from number size
    ! distribution to 3rd moment distribution
    IF (firstime) THEN
       xxm3 = 3.0_dp * xxlsg / sqrt2
       firstime = .FALSE.
    END IF

    ! Initialization
    mergenum     = 0.0_dp
    mergem3      = 0.0_dp
    mergenums    = 0.0_dp
    mergem3s     = 0.0_dp
    mergem3soot  = 0.0_dp
    mergenumsoot = 0.0_dp

    DO jp = 1, numcells

       !**********************************************************************
       ! *** merging of sooti into akns and accs
       !**********************************************************************

       epsilon=SUM(cblk(SOLSPEC,sooti,jp)) &
            / (SUM(cblk(1:nspec,sooti,jp)) - cblk(i_h2o,sooti,jp))

       IF(epsilon .GE. EPSILON_MAX) THEN  ! check if merging is necessary

          CALL renaming_fractions(status, 3, xxlsg((/sooti,akns,accs/))   &
               , dg((/sooti,akns,accs/),jp), cblk(i_num,(/akns,accs/),jp) &
               , xxm3(sooti), fnum, fm3)

          IF (status .EQ. 0) THEN

             ! Rename number concentration from sooti to akns and accs
             cblk(i_num,akns,jp) = cblk(i_num,akns,jp) &
                  + (1.0_dp - fnum) * cblk(i_num,sooti,jp)
             cblk(i_num,accs,jp) = cblk(i_num,accs,jp) &
                  +           fnum  * cblk(i_num,sooti,jp)

             ! Rename mass concentrations from sooti to akns and accs
             cblk(1:nspec,akns,jp) = cblk(1:nspec,akns,jp) &
                  + (1.0_dp - fm3)  * cblk(1:nspec,sooti,jp)
             cblk(1:nspec,accs,jp) = cblk(1:nspec,accs,jp) &
                  +           fm3   * cblk(1:nspec,sooti,jp)

             ! Reset sooti concentrations
             cblk(:,sooti,jp)       = CONMIN
             cblk(i_num,sooti,jp)   = nummin(sooti)

          END IF

       END IF

       !**********************************************************************
       ! *** merging of sootj into akns and accs/accs and cors
       !**********************************************************************

       epsilon=SUM(cblk(SOLSPEC,sootj,jp)) &
            / (SUM(cblk(1:nspec,sootj,jp)) - cblk(i_h2o,sootj,jp))

       IF(epsilon .GE. EPSILON_MAX) THEN  ! check if merging is necessary

          IF (dg(sootj,jp) * es36(sootj) .LT. dg(accs,jp) * es36(accs)) THEN
             mgmode(1) = akns
             mgmode(2) = accs
          ELSE
             mgmode(1) = accs
             mgmode(2) = cors
          END IF

          CALL renaming_fractions(status, 3           &
               , xxlsg((/sootj,mgmode(1),mgmode(2)/)) &
               , dg((/sootj,mgmode(1),mgmode(2)/),jp) &
               , cblk(i_num,(/mgmode(1),mgmode(2)/),jp), xxm3(sootj), fnum, fm3)

          IF (status .EQ. 0) THEN

             ! Rename number conc. from sootj to akns and accs/accs and cors
             cblk(i_num,mgmode(1),jp) = cblk(i_num,mgmode(1),jp) &
                  + (1.0_dp - fnum) * cblk(i_num,sootj,jp)
             cblk(i_num,mgmode(2),jp) = cblk(i_num,mgmode(2),jp) &
                  +           fnum  * cblk(i_num,sootj,jp)

             ! Rename mass conc.s from sootj to akns and accs/accs and cors
             IF (mgmode(1) .EQ. akns) THEN

                cblk(allexcssdu,mgmode(1),jp) = cblk(allexcssdu,mgmode(1),jp) &
                     + (1.0_dp - fm3)  * cblk(allexcssdu,sootj,jp)
                cblk(allexcssdu,mgmode(2),jp) = cblk(allexcssdu,mgmode(2),jp) &
                     +           fm3   * cblk(allexcssdu,sootj,jp)

                ! All SS (incl. Cl-) is assigned to accs, because SS particles
                ! are usually large and I (va) often use a config. with SS only
                ! in acc.
                cblk(i_cl,mgmode(2),jp) = cblk(i_cl,mgmode(2),jp) &
                     + cblk(i_cl,sootj,jp)
                cblk(i_ss,mgmode(2),jp) = cblk(i_ss,mgmode(2),jp) &
                     + cblk(i_ss,sootj,jp)
                ! In analogy, all DU is assigned to accs.
                cblk(i_du,mgmode(2),jp) = cblk(i_du,mgmode(2),jp) &
                     + cblk(i_du,sootj,jp)

             ELSE

                cblk(1:nspec,mgmode(1),jp) = cblk(1:nspec,mgmode(1),jp) &
                     + (1.0_dp - fm3)  * cblk(1:nspec,sootj,jp)
                cblk(1:nspec,mgmode(2),jp) = cblk(1:nspec,mgmode(2),jp) &
                     +           fm3   * cblk(1:nspec,sootj,jp)

             END IF

             ! Reset sootj concentrations
             cblk(:,sootj,jp)     = CONMIN
             cblk(i_num,sootj,jp) = nummin(sootj)

          END IF

       END IF

       !**********************************************************************
       ! *** merging of sootc into accs and cors
       !**********************************************************************

       epsilon=SUM(cblk(SOLSPEC,sootc,jp)) &
            / (SUM(cblk(1:nspec,sootc,jp)) - cblk(i_h2o,sootc,jp))

       IF(epsilon .GE. EPSILON_MAX) THEN  ! check if merging is necessary

          CALL renaming_fractions(status, 3, xxlsg((/sootc,accs,cors/))   &
               , dg((/sootc,accs,cors/),jp), cblk(i_num,(/accs,cors/),jp) &
               , xxm3(sootc), fnum, fm3)

          IF (status .EQ. 0) THEN

             ! Rename number concentration from sootc to accs and cors
             cblk(i_num,accs,jp) = cblk(i_num,accs,jp) &
                  + (1.0_dp - fnum) * cblk(i_num,sootc,jp)
             cblk(i_num,cors,jp) = cblk(i_num,cors,jp) &
                  +           fnum  * cblk(i_num,sootc,jp)

             ! Rename mass concentrations from sootc to accs and cors
             cblk(1:nspec,accs,jp) = cblk(1:nspec,accs,jp) &
                  + (1.0_dp - fm3)  * cblk(1:nspec,sootc,jp)
             cblk(1:nspec,cors,jp) = cblk(1:nspec,cors,jp) &
                  +           fm3   * cblk(1:nspec,sootc,jp)

             ! Reset sootc concentrations
             cblk(:,sootc,jp)     = CONMIN
             cblk(i_num,sootc,jp) = nummin(sootc)

          END IF

       END IF

       !**********************************************************************
       ! *** renaming of akn to acc
       !**********************************************************************

       ! *** use Binkowski-Kreidenweis paradigm

       ! Renaming criteria: 1. stronger growth of smaller mode, OR
       !                    2. larger number conc. in smaller mode AND diameter
       !                       of smaller mode exceeds threshold
       IF(cgr3(akn,jp) > cgr3(acc,jp) &
            .OR. dg(akn,jp) > 0.03e-6_dp &
            .AND. cblk(i_num,akn,jp) > cblk(i_num,acc,jp)) THEN

          ! Calculate number and mass fractions to be renamed from akn to acc
          CALL renaming_fractions(status, 2, xxlsg((/akn,acc/)) &
               , dg((/akn,acc/),jp), cblk(i_num,(/akn,acc/),jp) &
               , xxm3(akn), fnum, fm3)

          IF (status .EQ. 0) THEN   ! overlap is OK

             ! *** fnum and fm3 are the fractions of the number and 3rd moment
             !     distributions with  diameters greater than dd respectively.
             ! *** Rename mass from Aitken mode to acumulation mode. The mass
             !     transferred to the accumulation mode is proportional to the
             !     amount of 3rd moment transferred, therefore fm3 is used for
             !     mass transfer.

             ! Save renamed number concentration (for diagnostics only)
             MERGENUM(jp) = fnum * cblk(i_num,akn,jp)
             ! Rename particle number concentration from akn to acc
             cblk(i_num,acc,jp) = cblk(i_num,acc,jp) &
                  + fnum * cblk(i_num,akn,jp)
             ! Adjust akn number concentration
             cblk(i_num,akn,jp) = (1.0_dp - fnum) * cblk(i_num,akn,jp)

             ! Save renamed 3rd moment concentration (for diagnostics only)
             MERGEM3(jp) = fm3 * cblk(i_mom3,akn,jp)
             ! Rename particle mass concentration from akn to acc
             cblk(1:nspec,acc,jp) = cblk(1:nspec,acc,jp) &
                  + fm3 * cblk(1:nspec,akn,jp)
             ! Adjust akn mass concentration
             cblk(1:nspec,akn,jp) = (1.0_dp - fm3) * cblk(1:nspec,akn,jp)

!sorgam    cblk(VORGARO1J,jp) = cblk(VORGARO1J,jp)         &
!sorgam                            + cblk(VORGARO1I,jp) * fm3
!sorgam    cblk(VORGARO2J,jp) = cblk(VORGARO2J,jp)         &
!sorgam                            + cblk(VORGARO2I,jp) * fm3
!sorgam    cblk(VORGALK1J,jp) = cblk(VORGALK1J,jp)         &
!sorgam                            + cblk(VORGALK1I,jp) * fm3
!sorgam    cblk(VORGOLE1J,jp) = cblk(VORGOLE1J,jp)         &
!sorgam                            + cblk(VORGOLE1I,jp) * fm3
!sorgam    cblk(VORGBA1J,jp)  = cblk(VORGBA1J,jp)          &
!sorgam                            + cblk(VORGBA1I,jp) * fm3
!sorgam    cblk(VORGBA2J,jp)  = cblk(VORGBA2J,jp)          &
!sorgam                            + cblk(VORGBA2I,jp) * fm3
!sorgam    cblk(VORGBA3J,jp)  = cblk(VORGBA3J,jp)          &
!sorgam                            + cblk(VORGBA3I,jp) * fm3
!sorgam    cblk(VORGBA4J,jp)  = cblk(VORGBA4J,jp)          &
!sorgam                            + cblk(VORGBA4I,jp) * fm3

          END IF ! end check whether modal overlap is OK

       END IF ! end check on necessity for renaming akn to acc


       !**********************************************************************
       ! *** renaming of akns to accs
       !**********************************************************************

       ! Renaming criteria: 1. stronger growth of smaller mode, OR
       !                    2. larger number conc. in smaller mode AND diameter
       !                       of smaller mode exceeds threshold
       IF(cgr3(akns,jp) > cgr3(accs,jp) &
            .OR. dg(akns,jp) > 0.03e-6_dp &
            .AND. cblk(i_num,akns,jp) > cblk(i_num,accs,jp)) THEN

          ! Calculate number and mass fractions to be renamed from akns to accs
          CALL renaming_fractions(status, 2, xxlsg((/akns,accs/))   &
               , dg((/akns,accs/),jp), cblk(i_num,(/akns,accs/),jp) &
               , xxm3(akns), fnum, fm3)

          IF (status .EQ. 0) THEN   ! overlap is OK

             ! Save renamed number concentration (for diagnostics only)
             MERGENUMs(jp) = fnum * cblk(i_num,akns,jp)
             ! Rename particle number concentration from akns to accs
             cblk(i_num,accs,jp) = cblk(i_num,accs,jp) &
                  + fnum * cblk(i_num,akns,jp)
             ! Adjust akns number concentration
             cblk(i_num,akns,jp) = (1.0_dp - fnum) * cblk(i_num,akns,jp)

             ! Save renamed 3rd moment concentration (for diagnostics only)
             MERGEM3s(jp) = fm3 * cblk(i_mom3,akns,jp)
             ! Rename particle mass concentration from akns to accs
             cblk(1:nspec,accs,jp) = cblk(1:nspec,accs,jp) &
                  + fm3 * cblk(1:nspec,akns,jp)
             ! Adjust akns mass concentration
             cblk(1:nspec,akns,jp) = (1.0_dp - fm3) * cblk(1:nspec,akns,jp)

          END IF ! end check whether modal overlap is OK

       END IF ! end check on necessity for renaming akns to accs

       !**********************************************************************
       ! *** renaming of sooti to sootj
       !**********************************************************************

       ! Renaming criteria: 1. stronger growth of smaller mode, OR
       !                    2. larger number conc. in smaller mode AND diameter
       !                       of smaller mode exceeds threshold
       IF(cgr3(sooti,jp) > cgr3(sootj,jp) &
            .OR. dg(sooti,jp) > 0.03e-6_dp &
            .AND. cblk(i_num,sooti,jp) > cblk(i_num,sootj,jp)) THEN

          ! Calculate num. and mass fractions to be renamed from sooti to sootj
          CALL renaming_fractions(status, 2, xxlsg((/sooti,sootj/))     &
               , dg((/sooti,sootj/),jp), cblk(i_num,(/sooti,sootj/),jp) &
               , xxm3(sooti), fnum, fm3)

          IF (status .EQ. 0) THEN   ! overlap is OK

             ! Save renamed number concentration (for diagnostics only)
             MERGENUMsoot(jp) = fnum * cblk(i_num,sooti,jp)
             ! Rename particle number concentration from sooti to sootj
             cblk(i_num,sootj,jp) = cblk(i_num,sootj,jp) &
                  + fnum * cblk(i_num,sooti,jp)
             ! Adjust sooti number concentration
             cblk(i_num,sooti,jp) = (1.0_dp - fnum) * cblk(i_num,sooti,jp)

             ! Save renamed 3rd moment concentration (for diagnostics only)
             MERGEM3soot(jp) = fm3 * cblk(i_mom3,sooti,jp)
             ! Rename particle mass concentration from sooti to sootj
             cblk(1:nspec,sootj,jp) = cblk(1:nspec,sootj,jp) &
                  + fm3 * cblk(1:nspec,sooti,jp)
             ! Adjust sooti mass concentration
             cblk(1:nspec,sooti,jp) = &
                  (1.0_dp - fm3) * cblk(1:nspec,sooti,jp)

          END IF ! end check whether modal overlap is OK

       END IF ! end check on necessity for renaming sooti to sootj

       !**********************************************************************
       ! *** renaming of acc to cor
       !**********************************************************************

       ! Renaming criteria: 1. stronger growth of smaller mode, OR
       !                    2. larger number conc. in smaller mode AND diameter
       !                       of smaller mode exceeds threshold
!!$     IF(cgr3(acc,jp) > cgr3(cor,jp) .OR. &
!!$          dg(acc,jp) > 0.6e-6_dp .AND.    &
!!$          cblk(i_num,acc,jp) > cblk(i_num,cor,jp)) THEN
!!$
!!$        ! Calculate number and mass fractions to be renamed from acc to cor
!!$        CALL renaming_fractions(status, 2, xxlsg((/acc,cor/)),   &
!!$             dg((/acc,cor/),jp), cblk(i_num,(/acc,cor/),jp), &
!!$             xxm3(acc), fnum, fm3)
!!$
!!$        IF (status .EQ. 0) THEN   ! overlap is OK
!!$
!!$           ! Rename particle number concentration from acc to cor
!!$           cblk(i_num,cor,jp) = cblk(i_num,cor,jp) &
!!$                + fnum * cblk(i_num,acc,jp)
!!$           ! Adjust acc number concentration
!!$           cblk(i_num,acc,jp) = (1.0_dp - fnum) * cblk(i_num,acc,jp)
!!$
!!$           ! Rename particle mass concentration from acc to cor
!!$           cblk(1:nspec,cor,jp) = cblk(1:nspec,cor,jp) &
!!$                + fm3 * cblk(1:nspec,acc,jp)
!!$           ! Adjust acc mass concentration
!!$           cblk(1:nspec,acc,jp) = &
!!$                (1.0_dp - fm3) * cblk(1:nspec,acc,jp)
!!$
!!$        END IF ! end check whether modal overlap is OK
!!$
!!$     END IF ! end check on necessity for renaming acc to cor

     !**********************************************************************
     ! *** renaming of accs to cors
     !**********************************************************************

       ! Renaming criteria: 1. stronger growth of smaller mode, OR
       !                    2. larger number conc. in smaller mode AND diameter
       !                       of smaller mode exceeds threshold
!!$     IF(cgr3(accs,jp) > cgr3(cors,jp) .OR. &
!!$          dg(accs,jp) > 0.6e-6_dp .AND.    &
!!$          cblk(i_num,accs,jp) > cblk(i_num,cors,jp)) THEN
!!$
!!$        ! Calculate number and mass fractions to be renamed from accs to cors
!!$        CALL renaming_fractions(status, 2, xxlsg((/accs,cors/)),   &
!!$             dg((/accs,cors/),jp), cblk(i_num,(/accs,cors/),jp), &
!!$             xxm3(accs), fnum, fm3)
!!$
!!$        IF (status .EQ. 0) THEN   ! overlap is OK
!!$
!!$           ! Rename particle number concentration from accs to cors
!!$           cblk(i_num,cors,jp) = cblk(i_num,cors,jp) &
!!$                + fnum * cblk(i_num,accs,jp)
!!$           ! Adjust accs number concentration
!!$           cblk(i_num,accs,jp) = (1.0_dp - fnum) * cblk(i_num,accs,jp)
!!$
!!$           ! Rename particle mass concentration from accs to cors
!!$           cblk(1:nspec,cors,jp) = cblk(1:nspec,cors,jp) &
!!$                + fm3 * cblk(1:nspec,accs,jp)
!!$           ! Adjust accs mass concentration
!!$           cblk(1:nspec,accs,jp) = &
!!$                (1.0_dp - fm3) * cblk(1:nspec,accs,jp)
!!$
!!$        END IF ! end check whether modal overlap is OK
!!$
!!$     END IF ! end check on necessity for renaming accs to cors

     !**********************************************************************
     ! *** renaming of sootj to sootc
     !**********************************************************************

       ! Renaming criteria: 1. stronger growth of smaller mode, OR
       !                    2. larger number conc. in smaller mode AND diameter
       !                       of smaller mode exceeds threshold
!!$     IF(cgr3(sootj,jp) > cgr3(sootc,jp) .OR. &
!!$          dg(sootj,jp) > 0.6e-6_dp .AND.    &
!!$          cblk(i_num,sootj,jp) > cblk(i_num,sootc,jp)) THEN
!!$
!!$        ! Calculate number and mass fractions to be renamed from sootj to sootc
!!$        CALL renaming_fractions(status, 2, xxlsg((/sootj,sootc/)),   &
!!$             dg((/sootj,sootc/),jp), cblk(i_num,(/sootj,sootc/),jp), &
!!$             xxm3(sootj), fnum, fm3)
!!$
!!$        IF (status .EQ. 0) THEN   ! overlap is OK
!!$
!!$           ! Rename particle number concentration from sootj to sootc
!!$           cblk(i_num,sootc,jp) = cblk(i_num,sootc,jp) &
!!$                + fnum * cblk(i_num,sootj,jp)
!!$           ! Adjust sootj number concentration
!!$           cblk(i_num,sootj,jp) = (1.0_dp - fnum) * cblk(i_num,sootj,jp)
!!$
!!$           ! Rename particle mass concentration from sootj to sootc
!!$           cblk(1:nspec,sootc,jp) = cblk(1:nspec,sootc,jp) &
!!$                + fm3 * cblk(1:nspec,sootj,jp)
!!$           ! Adjust sootj mass concentration
!!$           cblk(1:nspec,sootj,jp) = &
!!$                (1.0_dp - fm3) * cblk(1:nspec,sootj,jp)
!!$
!!$        END IF ! end check whether modal overlap is OK
!!$
!!$     END IF ! end check on necessity for renaming sootj to sootc

    END DO ! loop for merging/renaming

  END SUBROUTINE made3_rename

!-------------------------------------------------------------------------------

  !> \brief Computes number and 3<sup>rd</sup> moment fractions to be renamed
  !>   between two modes
  !> \details Computes the number size distribution intersection diameter and
  !>   calculates the fractions of number and mass below and above that
  !>   diameter. Determines either the fraction of a smaller mode that is to be
  !>   renamed to a larger mode, or the fractions of an insoluble mode that are
  !>   to be renamed to the closest-sized mixed modes. The math is given by
  !>   \cite Binkowski2003.
  !> \note It is assumed that the modes are described in the following order:
  !>   - If (\e nmodes == 2): smaller mode first, larger mode second entry
  !>   - If (\e nmodes == 3):
  !>     - in \e logsigma and \e dg: mode to be merged into other 1st, smaller
  !>       receiving mode 2<sup>nd</sup>, larger receiving mode 3<sup>rd</sup>
  !>     - in \e num: smaller receiving mode 1st, larger receiving mode
  !>       2<sup>nd</sup>

  SUBROUTINE renaming_fractions(status, nmodes, logsigma, dg, num, xxm3 &
       , fnum, fm3)

  !   Purpose
  !   -------
  !   Determines the number and third moment fractions of the first mode that is
  !   described by the arrays logsigma and dg, that should be renamed to the
  !   second mode described by the array num.
  !   NOTE: 
  !
  !   Parameter list (i=input, o=output):
  !   -----------------------------------
  !   i-  xxm3          conversion ``factor'' for error function argument from
  !                     number size distribution to 3rd moment distribution
  !   -o  fnum          fraction of the number concentration that should be
  !                     assigned to the larger receiving mode
  !   -o  fm3           fraction of the 3rd moment concentration that should be
  !                     assigned to the larger receiving mode

    IMPLICIT NONE
    INTRINSIC LOG, SIGN, SQRT, EXP, MAX

    ! I/O
    !> Error status flag
    INTEGER,  INTENT(out) :: status
    !> Number of modes to be considered (= 2 or 3)
    INTEGER,  INTENT(in)  :: nmodes
    !> \f$\ln(\sigma_i)\f$
    REAL(dp), INTENT(in)  :: logsigma(nmodes)
    !> Wet median modal diameters [m]
    REAL(dp), INTENT(in)  :: dg(nmodes)
    !> Number concentrations of the receiving modes [m<sup>-3</sup>]
    REAL(dp), INTENT(in)  :: num(2)
    !> Conversion "factor" for erf argument (m-3 -> mom_3 m-3)
    REAL(dp), INTENT(in)  :: xxm3
    !> Number concentration fraction to be assigned to larger receiving mode
    REAL(dp), INTENT(out) :: fnum
    !> 3<sup>rd</sup> moment concentration fraction to be assigned to larger
    !> receiving mode
    REAL(dp), INTENT(out) :: fm3

    ! LOCAL
    REAL(dp) :: alfa, yji, aa, bb, getaf_l, cc, disc, aaa   ! helpers
    REAL(dp) :: d_0l      ! intersection diameter of number size distributions
                          ! of modes 2 and 3
    REAL(dp) :: xnum      ! error function arg. at d_0l or aaa, resp. (number)


    ! Initialization
    status = 499   ! Unspecified ERROR
    fnum   = 1000._dp
    fm3    = 1000._dp

    ! Set helper variables
    alfa    = logsigma(nmodes-1) / logsigma(nmodes)
    yji     = LOG(dg(nmodes) / dg(nmodes-1)) / (SQRT2 * logsigma(nmodes-1))
    aa      = 1.0_dp - alfa * alfa
    bb      = 2.0_dp * yji * alfa * alfa
    getaf_l = LOG(alfa * num(2) / num(1))
    cc      = getaf_l - yji * yji * alfa * alfa
    disc    = bb * bb - 4.0_dp * aa * cc

    IF (disc < 0.0_dp) THEN
       status = 401   ! ERROR in intersection
       RETURN
    ELSE
       aaa = -2.0_dp * cc / (bb + SIGN(1.0_dp, bb) * SQRT(disc))
    END IF

    ! Calculate fractions to be renamed
    IF (nmodes == 3) THEN
       d_0l  = dg(2) * EXP(aaa * SQRT2 * logsigma(nmodes-1))
       xnum  = LOG(d_0l / dg(1)) / (SQRT2 * logsigma(1))
       fnum  = 0.5_dp * erfc(xnum)
       fm3   = 0.5_dp * erfc(xnum - xxm3)
    ELSE
       xnum  = MAX(aaa, xxm3)
       IF (xnum .GE. xxm3) THEN
          fnum  = 0.5_dp * erfc(xnum)
          fm3   = 0.5_dp * erfc(xnum - xxm3)
       ELSE
          fnum  = 0.5_dp * erfc(xxm3)
          fm3   = 0.5_dp
       END IF
    END IF

    status = 0

  END SUBROUTINE renaming_fractions

!-------------------------------------------------------------------------------

  !> \brief Calculates modal (relative) and absolute (total) condensation
  !>   factors (used to calculate condensation rates)
  !> \details Calculates the factor by which a gas phase concentration is
  !>   multiplied to obtain the condensation rate, and the fractional modal
  !>   contributions to this factor using the method of \cite Binkowski1995
  !>   (i.e., harmonic mean of free molecular and near continuum regime
  !>   expressions).

  SUBROUTINE cond_factors(i_spec, numcells, press, temp &
       , numconc, diam, confac_tot, confac_rel)

    IMPLICIT NONE
    INTRINSIC SQRT, SUM

    ! I/O
    !> Index of condensing gas species
    INTEGER,  INTENT(in)  :: i_spec
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)  :: numcells
    !> Air pressure [Pa]
    REAL(dp), INTENT(in)  :: press(numcells)
    !> Air temperature [K]
    REAL(dp), INTENT(in)  :: temp(numcells)
    !> Modal number concentrations [m<sup>-3</sup>]
    REAL(dp), INTENT(in)  :: numconc(nmod,numcells)
    !> Wet median modal diameters [m]
    REAL(dp), INTENT(in)  :: diam(nmod,numcells)
    !> Total condensation factor [s<sup>-1</sup>]
    REAL(dp), INTENT(out) :: confac_tot(numcells)
    !> Relative condensation factors [-]
    REAL(dp), INTENT(out) :: confac_rel(nmod,numcells)

    ! LOCAL
    INTEGER  :: jl, jm      ! loop indices
    REAL(dp) :: cconc       ! condensation coefficient for near continuum regime
    REAL(dp) :: ccofm(nmod) ! condensation coeff.s for free molecular regime
    REAL(dp) :: diffcorr    ! diffusion correction factor for nc regime
    REAL(dp) :: csqt(nmod)  ! temperature corrected condensation coefficient for
                            ! fm regime
    REAL(dp) :: am1, am2    ! first and second moments of aersol num. size dist.
    REAL(dp) :: gnc3, gfm3  ! cond. fac.s for indiv. modes in nc and fm regimes


    ! Initialization
    confac_tot = 0.0_dp
    confac_rel = 0.0_dp

    ! Condensation coefficients
    ! FSB  CCOFM is  the accommodation coefficient times the mean molecular
    !      velocity for h2so4 without the temperature after some algebra
    ! ... free molecular regime
    ccofm = alpha(:,i_spec) &
         * SQRT(pi * rgasuniv / (2.0e-3_dp * mw(i_spec,i_mwgas)))
    ! ... near continuum regime
    ! *** CCONC is the factor for near-continuum condensation times the mean
    !     molecular velocity without the temperature after some algebra
    cconc = 2.0_dp * pi * diff(i_spec)

    DO jl = 1, numcells

       ! near-continuum factors [1/s]
       !bs * adopted from code of FSB
       !bs * correction to gas diffusivities for temperature and pressure
       diffcorr = (P0 / press(jl)) * (temp(jl) / 273.16_dp)**1.75_dp

       ! free molecular factors
       csqt = ccofm * SQRT(temp(jl)) ! put in temperature factor

       DO jm = 1, nmod
          ! *** First moment:
          am1 = numconc(jm,jl) * diam(jm,jl) * es04(jm)
          ! *** Second moment:
          am2 = numconc(jm,jl) * diam(jm,jl) * diam(jm,jl) * es16(jm)
          ! Compute 3rd mom. condensation rate as the harmonic mean of the
          ! corresponding near continuum and free molecular rates.
          ! *** NOTE: confac_rel will be redefined below <<<<<<
          gnc3 = cconc * am1 * diffcorr
          gfm3 = csqt(jm) * am2
          confac_rel(jm,jl) = gnc3 * gfm3 / (gnc3 + gfm3)
       END DO

       ! Note: modes for which condensation is not computed should not be
       ! included in confac_tot.
       confac_tot(jl) = SUM(confac_rel(:,jl))

       ! *** redefine confac_rel to be the nondimensional fractional
       !     condensation factors
       confac_rel(:,jl) = confac_rel(:,jl) / confac_tot(jl)

    END DO

  END SUBROUTINE cond_factors

!-------------------------------------------------------------------------------

  !> \brief Complementary error function
  !> \details Complementary error function from Numerical Recipes, accurate to
  !>   within 1.2e-7 for all x.
  !> \returns 1 - erf(x) as a real(dp) value

  ELEMENTAL REAL(dp) FUNCTION erfc(x)

    IMPLICIT NONE
    INTRINSIC ABS, EXP

    !> Argument of the complementary error function
    REAL(dp), INTENT(in) :: x

    ! Local
    REAL(dp) :: t, z

    z = ABS(x)
    t = 1.0_dp / ( 1.0_dp + 0.5_dp * z )

    erfc = t * EXP(-z * z - 1.26551223_dp + t        &
         * ( 1.00002368_dp + t * ( 0.37409196_dp + t &
         * ( 0.09678418_dp + t * (-0.18628806_dp + t &
         * ( 0.27886807_dp + t * (-1.13520398_dp + t &
         * ( 1.48851587_dp + t * (-0.82215223_dp + t &
         * 0.17087277_dp )))))))))

    IF ( x.lt.0.0_dp ) erfc = 2.0_dp - erfc

  END FUNCTION erfc

!-------------------------------------------------------------------------------

END MODULE messy_made3
