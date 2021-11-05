! ==============================================================================
! {%CMODEL}_tag_common
! generated: {%TIMEDATE}
!
! this module is created automatically by imtag tool - do not edit
!
! inter-configuration module: common utils, consts
! level: smcl
!
! {$TAG_INFO} ! this is a template file for imtag utility
!
! [Gromov, MPIC, 2007-2020]
! ==============================================================================

!> \brief MECCA-TAG common core module.

!> \authors Sergey Gromov, MPI-C, 2007-2020

!> \version 2.6
!!   - improvements & bug-fixes to tagging code and integration with chemprop facility
!!   - new & updated tagging configurations (CIO, FN, FO, FS)
!!   - refurbished C & O transfer in MIM/MOM mechanisms
!!   - parameter definition restyled for more flexible usage (those global in *_parameters.inc can be overridden on the module/cfg. level)
!!   - added possibility to have class-wise fixed (nudged) composition
!!   - negation (!) and regexp can be used for KIE and source specification records for both equation tags and species names
!!   - tagged species tolerances in KPP are automatically set based those for regular species

!> \version 2.5
!!   - added isotope branching
!!   - total element diagnostic is replaced with improved i-mass balance checks and unaccounted prod./loss
!!   - optimisations for simulations with MC box (cfg. parameter TAG_OPT_SKIP_XTRABOXPPROC)
!!   - improved tracer definition files creation, including the new MESSy chemprop facility
!!   - improvements from F.Frank and P.Joeckel on isotope H chemistry (cfg) and SMIL code
!!   - isotopologue molar masses can be rounded to given # of digits after comma (option roundmass= )
!!   - revived .dot files creation (option dots= )
!!   - several important tagging and F90 code bugfixes (see CHANGELOG)

!> \version 2.0
!!   - refurbished version (obsolete dbl->tag)
!!   - added H stable isotope chemistry
!!   - added 14C composition handling routines
!!   - isofrac routines renamed for consistency ( #r->#a, #f->#r for abun. and rare, resp.)
!!   - added VPBDB-CO2 reference standard for O and delta-17O calc. from delta-18O and cap-Delta17O

!> \version 1.0
!!   - initial release
!!   - stable C and O modelling
!!   - 2&3-isotope composition converstion/delta calc. routines
!!   - passive tracer (PT) handling

!> \\todo
!!   - variable n-isotope comp. converstion/calc. routines

!>- general tagging parameters (as conditional defines) ------------------------

#include "{%CMODEL}_tag_parameters.inc"

! ------------------------------------------------------------------------------

module {%CMODEL}_tag_common

  use messy_mecca_kpp ! dp, ...

  implicit none
  save

! module and version info
  character(len=*), parameter, public :: submodstr = '{%CMODEL}_tag'
  character(len=*), parameter, public :: submodver = '2.6'

  real(dp), parameter :: UNDEF = -1E+34_dp      ! undefined value mask

  integer, parameter  :: tag_NPT_max = 512      ! expected max. no. of PTs
  integer             :: tag_PTI(tag_NPT_max)   ! array of PTs tracer indices
  integer             :: tag_NPT = 0            ! no. of PTs found in the current mech

  character(len=*), &
  parameter           :: &
! conventional masks to id a PT in SPEC_NAMES
#ifndef iCMCb
                         idm_PT(2)  = (/ 'XPTP', 'XPTL' /)
#else
! special for isoCO2 MC box
                         idm_PT(6)  = (/ 'FX   ', &
                         'I12FX', 'I13FX', &
                         'I16FX', 'I17FX', 'I18FX' /)
#endif
! reference standard ratio for 2H, V-SMOW
  real(dp), parameter :: VSMOW_2H   =  155.76e-6_dp  ! +-.05  [Hagemann et al. 1970]

! reference standard ratio for 13C, V-PDB
  real(dp), parameter :: VPDB_13C   = 1123.72e-5_dp  ! +-.60  [Craig 1957]
! real(dp), parameter :: VPDB_13C   = 1120.20e-5_dp  ! +-.28  [Zhang et al. 1990]

! reference standard ratios for 17O & 18O, V-SMOW [Assonov & Breninkmeijer, RCM 2003, PC]
! VPDB-CO2 scale
  real(dp), parameter :: VPDB_17O   =  395.11e-6_dp
  real(dp), parameter :: VPDB_18O   = 2088.35e-6_dp
! VSMOW scale
  real(dp), parameter :: VSMOW_17O  =  386.72e-6_dp
  real(dp), parameter :: VSMOW_18O  = 2005.20e-6_dp  ! +-.45  [Gonfiantini 1978]
! air oxygen
  real(dp), parameter :: AOX_17O    = (1.+12.03e-3)*VSMOW_17O  ! +12.08 pm vs. VSMOW
  real(dp), parameter :: AOX_18O    = (1.+23.88e-3)*VSMOW_18O  ! +23.88 pm vs. VSMOW
! (mass-dependent) fractionation slopes (beta)
  real(dp), parameter :: MDFSL_MWL  =    0.5281_dp   !        Meteoric waters line [2003.CR103.Brenninkmeije,etal->1998.IEHS219.Li&Meijer]
  real(dp), parameter :: MDFSL_LVE  =    0.5279_dp   ! +-.001 Liq.-vap. equilibrium [2005.RCM19.Luz&Barkan]
  real(dp), parameter :: MDFSL_CO2  =    0.516_dp    !        [2004.GRL31.Boering,etal]

  public

contains

! ==============================================================================
! isotope-related functions and subroutines

!> \details calculation of the abundant isotope atoms number from 2 isotopologues
!!          concentration and number of constituent isotope-tagged atoms

  elemental real(dp) function abun2iso(abun, rare, atoms)

    implicit none
    real(dp), intent(in) :: abun, rare  !< abundant and rare isotopologue (molecular) abundance
    integer,  intent(in) :: atoms       !< relevant isotope quantity in molecule

    abun2iso = ( abun + rare ) * real(atoms,dp) - rare

  end function abun2iso

! ------------------------------------------------------------------------------

! calculation of the abundant isotope atoms number from 3 isotopologues
! concentration and number of constituent isotope-tagged atoms
  elemental real(dp) function abun3iso(abun, rare_cur, rare_oth, atoms)

    implicit none
    real(dp), intent(in) :: abun, rare_cur, rare_oth
    integer,  intent(in) :: atoms

    abun3iso = ( abun + rare_cur + rare_oth ) * real(atoms,dp) - &
                      ( rare_cur + rare_oth )

  end function abun3iso

! ------------------------------------------------------------------------------

! calculation of the isotopologue ratio from abundant and rare isotopologue
! molecules concentration and number of constituent isotope-tagged atoms
  elemental real(dp) function isoR2m(abun, rare, atoms)

    implicit none
    real(dp), intent(in) :: abun, rare
    integer,  intent(in) :: atoms

    isoR2m = rare / abun2iso(abun,rare,atoms)

  end function isoR2m

! ------------------------------------------------------------------------------

! calculation of the isotopologue ratio from abundant and two rare isotopologue
! molecules concentration and number of constituent isotope-tagged atoms
  elemental real(dp) function isoR3m(abun, rare_cur, rare_oth, atoms)

    implicit none
    real(dp), intent(in) :: abun, rare_cur, rare_oth
    integer,  intent(in) :: atoms

    isoR3m = rare_cur / abun3iso(abun, rare_cur, rare_oth, atoms)

  end function isoR3m

! ------------------------------------------------------------------------------

! calculation of the isotopologue ratio from the delta and reference ratio
  elemental real(dp) function isoRd(delta, Rst)

    implicit none
    real(dp), intent(in) :: delta, Rst

    isoRd = Rst * (delta + 1.0_dp)

  end function isoRd

! ------------------------------------------------------------------------------

! calculation of the delta value from abundant and rare isotopologue molecules
! concentration, number of constituent isotope-tagged atoms, reference ratio
! not in per mil
  elemental real(dp) function delta2(abun, rare, Rst, atoms)

    implicit none
    real(dp), intent(in) :: abun, rare, Rst
    integer,  intent(in) :: atoms

    if ( abun2iso(abun,rare,atoms) .NE. 0.0_dp ) then
      delta2 = ( isoR2m(abun,rare,atoms) / Rst - 1.0_dp )
    else
      delta2 = UNDEF
    endif

  end function delta2

! ------------------------------------------------------------------------------

!> \details calculation of the delta value from abundant and two rare isotopologue molecules
!!          concentration, number of constituent isotope-tagged atoms, reference ratio
  elemental real(dp) function delta3(abun, rare_cur, rare_oth, Rst_cur, atoms)

    implicit none
    real(dp), intent(in) :: abun, rare_cur, rare_oth, Rst_cur
    integer,  intent(in) :: atoms

    if ( abun3iso(abun,rare_cur,rare_oth,atoms) .NE. 0.0_dp ) then
      delta3 = ( isoR3m(abun,rare_cur,rare_oth,atoms) / Rst_cur - 1.0_dp )
    else
      delta3 = UNDEF
    endif

  end function delta3

! ------------------------------------------------------------------------------

!> \details calculation of the rare isotopologue fraction
!>          using delta (o/oo) & R ref. values and molecule atoms no.
  elemental real(dp) function isofrac2r(delta, Rst, atoms)

    implicit none
    real(dp), intent(in) :: delta       ! not in per mil
    real(dp), intent(in) :: Rst
    integer,  intent(in) :: atoms
    real(dp)             :: gamma

    gamma = isoRd(delta,Rst)
    isofrac2r = ( gamma * real(atoms,dp) ) / ( gamma + 1.0_dp )

  end function isofrac2r

! -----------------------------------------------------------------------------

!> \details calculation of the abundant isotopologue fraction
!!          using delta (o/oo) & R ref. values and molecule atoms no.
  elemental real(dp) function isofrac2a(delta, Rst, atoms)

    implicit none
    real(dp), intent(in) :: delta       ! not in per mil
    real(dp), intent(in) :: Rst
    integer,  intent(in) :: atoms
    real(dp)             :: gamma

    gamma = isoRd(delta,Rst)
    isofrac2a = ( gamma * real(1 - atoms,dp) + 1.0_dp ) / ( gamma + 1.0_dp )

  end function isofrac2a

! -----------------------------------------------------------------------------

!> \details calculation of the budget fraction of the first of two rare
!!          isotopologues in case of three isotopologues tagging
  elemental real(dp) function isofrac3r(delta1, R1, delta2, R2, atoms)

    implicit none
    real(dp), intent(in) :: delta1, delta2        ! not in per mil
    real(dp), intent(in) :: R1, R2
    integer,  intent(in) :: atoms

    isofrac3r = isoRd(delta1,R1) * real(atoms,dp) / &
                  ( (isoRd(delta1,R1) + isoRd(delta2,R2)) + 1.0_dp )

  end function isofrac3r

! -----------------------------------------------------------------------------

!> \details calculation of the abundant isotopologue fraction in case of three isotopologues
  elemental real(dp) function isofrac3a(delta1, R1, delta2, R2, atoms)

    implicit none
    real(dp), intent(in) :: delta1, delta2        ! not in per mil
    real(dp), intent(in) :: R1, R2
    integer,  intent(in) :: atoms
    real(dp)             :: gamma

    gamma = real(atoms,dp)
    isofrac3a =  gamma / ( (isoRd(delta1,R1) + isoRd(delta2,R2)) + 1.0_dp ) - &
                 gamma + 1.0_dp

  end function isofrac3a

! -----------------------------------------------------------------------------

! =============================================================================
! RADIOCARBON

!> \details calculation of the 14C/12C (atomic!) ratio from a given pMC value
  elemental real(dp) function isoRpMC(pMC, d13C)

    implicit none
    real(dp), intent(in) :: pMC             !< 14C activity in pMC (per cent modern carbon)
    real(dp), intent(in), optional :: d13C  !< d13C of the sample (not in per mil)
    real(dp)             :: ad13C

    if (present(d13C)) then
      ad13C = d13C
    else
      ad13C = -19.0e-3_dp  !< Default is -19 per mil of NBS(1) oxalic acid standard corrected 1950
    endif

    isoRpMC = ( pMC / 100.0_dp ) * ( (1.0_dp+ad13C) / 0.975_dp )**2 * 1.189e-12

  end function isoRpMC

! -----------------------------------------------------------------------------

!> \details calculation of the 14C isotopologue fraction using pMC (plus optional d13C) value
  elemental real(dp) function isofrac2r_pMC(pMC, atomsC, d13C)

    implicit none
    real(dp), intent(in) :: pMC             !< 14C activity in pMC (per cent modern carbon)
    integer,  intent(in) :: atomsC          !< no. of C atoms in isotopologue
    real(dp), intent(in), optional :: d13C  !< d13C of the sample (not in per mil)
    real(dp)             :: gamma

    gamma = isoRpMC(pMC, d13C)
    isofrac2r_pMC = ( gamma * real(atomsC,dp) ) / ( gamma + 1.0_dp )

  end function isofrac2r_pMC

! -----------------------------------------------------------------------------

!> \details calculation of the pMC (percent modern carbon) 14C content 
!!          using 14C and 12C (plus optional d13C) values
  elemental real(dp) function ratio_pMC(abun, rare, atomsC, d13C)

    implicit none
    real(dp), intent(in) :: abun, rare      !< 12C and 14C abundance
    integer,  intent(in) :: atomsC          !< no. of C atoms in isotopologue
    real(dp), intent(in), optional :: d13C  !< d13C of the sample (not in per mil)

    ratio_pMC = ( rare / ( abun*real(atomsC,dp) + rare*real(atomsC-1,dp) ) ) / isoRpMC(1.0_dp, d13C)

  end function ratio_pMC


! =============================================================================
! OTHER CONVENTIONAL CONVERSIONS

!> \details
!! kierate gives a factor one multiplies the standard (abundant isotopologue)
!! reaction rate to account for the KIE equivalent to a given eps value
!! eps = (Ki/Kj - 1), denotes enrichment (depletion when <0) in leftover compartment
!! input eps should be in per mil (o/oo)
  elemental real(dp) function kierate(eps)
    implicit none
    real(dp), intent(in) :: eps         ! not in per mil
    kierate = 1.0_dp / ( eps * 1e-3_dp + 1.0_dp )
  end function kierate

!> \details
!! calculation of delta-17O from given delta-18O and cap. Delta-17O values
!! input deltas should be in per mil (o/oo)

  elemental real(dp) function delta17Opm(delta18O, capDelta17O, opt_MDFSL_O)
    implicit none
    real(dp), intent(in)           :: delta18O, capdelta17O  !< in permil
    real(dp), intent(in), optional :: opt_MDFSL_O            !< optional MDF slope value (beta) for cap-D17O
    real(dp)                       :: act_MDFSL_O            !< actual MDF slope value used in calculation
    if (present(opt_MDFSL_O)) then
      act_MDFSL_O = opt_MDFSL_O    !< custom fractionation slope
    else
      act_MDFSL_O = MDFSL_MWL      !< default is MWL
    endif
    delta17Opm = ( ( capDelta17O * 1e-3 + 1_dp ) * ( delta18O * 1e-3_dp + 1_dp ) ** act_MDFSL_O - 1_dp ) * 1e3_dp
  end function delta17Opm

! ==============================================================================
! MESSy/MECCA-related functions and subroutines

  subroutine mecca_tag_scanPTs(info)

    implicit none

    character(len=*), intent(out) :: info
    character(len=32)             :: str
    integer :: js, mi

    intrinsic :: size, mod

    info = ''
    str = ''
    tag_NPT = 0

    spec_loop: do js = 1, NSPEC

    ! checking if spec name satisfies idm_PT values
      do mi = 1, size(idm_PT)
!       if ( SPC_NAMES(js)(1:len_trim(idm_PT(mi))) .eq. trim(idm_PT(mi)) ) then
        if ( index(trim(SPC_NAMES(js)),trim(idm_PT(mi))) == 1 ) then
        ! it's a PT
          tag_NPT = tag_NPT + 1
          write(str,'(I0)') tag_NPT
          info = trim(info)//' '//trim(SPC_NAMES(js))//'('//trim(str)//')'
          if ( tag_NPT .le. tag_NPT_max ) then
            tag_PTI(tag_NPT) = js
          else
          ! we ran out of free entries in tag_PTI, do nothing
            info = trim(info)//'*'
          endif
        ! output 5 entries/line
          if ( mod(tag_NPT, 5) .eq. 0 ) info = trim(info)//new_line('A')
        endif
      enddo
    enddo spec_loop

    if (len_trim(info) .gt. 0) info = ': '//new_line('A')//trim(info)
    write(str,'(I0)') tag_NPT
    info = 'mecca_tag_scanPTs(): found ( '//trim(str)//' ) passive tracers'//trim(info)

    if ( tag_NPT .gt. tag_NPT_max ) info = trim(info)//new_line('A')// &
      'warning: too many PTs found, increase tag_NPT_max to fit extra entries (marked with *)'

#ifdef DEBUG
    print *,'mecca_tag_scanPTs(): passed'
#endif

  end subroutine mecca_tag_scanPTs


! -----------------------------------------------------------------------------

  subroutine mecca_tag_resetPTs(C)

    implicit none

  ! concentrations vector
    real(dp), intent(inout) :: C(:)

  ! resetting PTs values
    C(tag_PTI(1:tag_NPT)) = 0.0_dp

#ifdef DEBUG
    print *,'mecca_tag_resetPTs: passed'
#endif

  end subroutine mecca_tag_resetPTs



! -----------------------------------------------------------------------------

  subroutine mecca_tag_intPTs2arr(C, time_step_len)

    implicit none

  ! concentrations vector
    real(dp), intent(inout) :: C(:)
  ! integration step
    real(dp), intent(in)    :: time_step_len

  ! converting PTs integral values to average reaction rates
    C(tag_PTI(1:tag_NPT)) = C(tag_PTI(1:tag_NPT)) / time_step_len

#ifdef DEBUG
    print *,'mecca_tag_intPTs2arr(', time_step_len, '): passed'
#endif

  end subroutine mecca_tag_intPTs2arr



end module {%CMODEL}_tag_common

! *****************************************************************************

