! ***************************************************************************
MODULE  messy_rad_fubrad_mem


  USE messy_main_constants_mem, ONLY: dp, cpd => cp_air
  
  IMPLICIT NONE
  PRIVATE
  SAVE

  ! NAME OF SUB-SUBMODEL
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr = 'rad_fubrad'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodver = '1.9'

  ! Parameter
  ! 
  ! lya_eff - efficiency factor for lyman alpha
  !           (Mlynczak and Solomon, JGR, vol 98, p 10517, 1993)
  !
  REAL(dp), PARAMETER, PUBLIC :: lya_eff=0.95_dp 
  !
  ! alos    - Loschmidt's number in cm^-3
  !           Parameter alos moved here from o3fluxes.
  !
  !REAL(dp), PARAMETER ::  alos=2.687E19_dp
  !
  ! cgs2si  - parameter to convert flux data provided in
  !           erg cm-2 s-1 to W m-2
  !           The cgs unit erg cm-2 s-1 is equivalent to
  !           the SI unit mW m-2.
  !
  REAL(dp), PARAMETER, PUBLIC :: cgs2si = 1.E-3_dp
  !
  ! oxndc - (almost) oxygen number density
  !         used as a factor to determine o2 column
  !         oxndc=O2_VMR*(M_O2/M_air)*(N_Avogadro/M_O2)*(1e-4/g)
  !         Equation to determine O2_column:
  !         O2_VMR*(M_O2/M_air)*Air density*(N_Avogadro/M_O2)*1e-6
  !         substitute air density by using hydrostatic balance and
  !         integrate from the model top to layer
  !
  !REAL(dp), PARAMETER :: oxndc=4.442E19_dp
  !
  ! oxcnst - constant to derive number density in cm-3 from VMR:
  !          oxcnst=O2_VMR*(M_O2/M_air)*Air density*(N_Avogadro/M_O2)*1e-6
  !          1e-6 converts from cm-3 to m-3.
  !          The density factor is ignored as it cancels when the
  !          heating rate in Ks-1 is computed from that in Wm-3
  !
  !REAL(dp), PARAMETER :: oxcnst=4.3549E18_dp
  !
  ! oxfac  - 
  REAL(dp), PARAMETER, PUBLIC :: oxfac=1.E2_dp / cpd  ! 1.e2= cm-1 -> m-1
  !
  ! Variables
  !
  ! lfubrad  - switch, if sub-submodel is active
  ! nswlev   - Number of levels on which the fub scheme operates
  !
  LOGICAL, PUBLIC :: lfubrad  = .FALSE.
  INTEGER, PUBLIC :: nswlev = 0
  INTEGER, PUBLIC :: inum_mo_radioflux
  INTEGER, PUBLIC :: my_pe=-1          ! for testprints ONLY
  !
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: po3c => NULL()

  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: po2c => NULL()
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: przsec

  REAL(dp), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: my_pph
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: my_ppf  ! fb_mk_20150930
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: prmu0
  
  !
  ! SCHEME FIELDS
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flup   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: fldo   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flhart => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flupc  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: fllya  => NULL()

  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flhuup => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flhudo => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flchup => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flchdo => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flhz   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flsrb  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER, PUBLIC :: flsrc  => NULL()

  REAL(dp), DIMENSION(:),   POINTER, PUBLIC :: altsw  => NULL()

  ! tropospheric albedo, short wave, clear sky
  REAL(dp), DIMENSION(:),   POINTER, PUBLIC :: altswc => NULL()

  !
  ! Variables to store TOA solar fluxes:
  !
  ! zFlya   - Lyman alpha             (121.5 nm)  in W m^-2
  ! zFsch1  - Schuman-Runge Continuum (125-152nm) in erg cm-2 s-1 (mW m^-2)
  ! zFs     - Schuman-Runge Continuum (152-166nm) in erg cm-2 s-1 (mW m^-2)
  ! zFl     - Schuman-Runge Continuum (166-175nm) in erg cm-2 s-1 (mW m^-2)
  ! zFd     - difference, zFs - zFl
  ! soscale - Schumann-Runge Bands    (175-205nm) dimensionless
  !
  REAL(dp), PUBLIC :: zFlya
  REAL(dp), PUBLIC :: zFsch1  ! integrated flux from 125 - 152nm
  REAL(dp), PUBLIC :: zFs     ! s: short , zFs = Is/M = 0.5 * integrated solar flux 
  !                           from 152-166nm (in erg cm-2 s-1)
  REAL(dp), PUBLIC :: zFl     ! l: long , zFl = Il/M = 0.6 * integrated solar flux 
  !                           from 166-175nm (in erg cm-2 s-1)
  REAL(dp), PUBLIC :: zFd     ! d: difference, zFs-zFl

  REAL(dp), PUBLIC :: soscale ! Schumann-Runge Bands (175-205nm)
  !                           dimensionless solar scaling factor: 
  !                           mean (=1), max (>1), min (<1)
  ! Base level for FUBRAD, i.e. the lowest (w.r.t height) 
  ! pressure level in Pa where FUBRAD is operating. Default value
  ! is 7000 Pa. Can be changed via namelist.
  REAL(dp), PUBLIC :: blev
  !
  ! zFsrc      - Flux in Schumann-Runge Continuum (124.0 - 174.5 nm)
  ! zFsrc_max  - predefined fluxes for solar maximum conditions
  ! zFsrc_min  - predefined fluxes for solar minimum conditions
  ! zsigsrc    - O2 cross section Schumann-Runge continuum
  !
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFsrc
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFsrc_max, zFsrc_min
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zsigsrc
  !
  ! zFsrb      - Flux in Schumann-Runge Bands (174.5  -202.5 nm)
  ! zFsrb_max  - predefined fluxes for solar maximum conditions
  ! zFsrb_min  - predefined fluxes for solar minimum conditions
  !
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFsrb
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFsrb_max, zFsrb_min
  !
  ! zFherz      - Flux in Herzberg continuum (206.186-243.902 nm, WMO band 18-32)
  ! zFherz_max  - predefined fluxes for solar maximum conditions
  ! zFherz_min  - predefined fluxes for solar minimum conditions
  ! zsigherz2   - O2 cross section Herzberg continuum
  ! zsigherz3   - O3 cross section Herzberg continuum
  !
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFherz
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFherz_max, zFherz_min
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zsigherz2
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zsigherz3
  !
  ! zFhart      - Flux in Hartley bands      (243.902-277.778 nm, WMO band 32-42) 
  ! zFhart_max  - predefined fluxes for solar maximum conditions
  ! zFhart_min  - predefined fluxes for solar minimum conditions
  ! zsighart    - O3 cross section Hartley bands
  !
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFhart
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFhart_max, zFhart_min
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zsighart
  !
  ! zFhug      - Flux in Huggins bands       (277.778-362.500 nm, WMO band 42-60) 
  ! zFhug_max  - predefined fluxes for solar maximum conditions
  ! zFhug_min  - predefined fluxes for solar minimum conditions
  ! zsighug    - O3 cross section Huggins bands
  !
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFhug
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFhug_max, zFhug_min
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zsighug
  !
  ! zFchap      - Flux in Chappuis bands     (407.5  -690.0   nm, WMO band 70-126)
  ! zFchap_max  - predefined fluxes for solar maximum conditions
  ! zFchap_min  - predefined fluxes for solar minimum conditions
  ! zsigchap    - O3 cross section Chappuis bands
  !
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFchap
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zFchap_max, zFchap_min
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zsigchap
  !
  ! zFadd, zFadd_max  -  Additional flux in bands (362.5 - 407.5 nm, WMO band 61-69)
  !        zFadd_min     (non absorbing spectral regions)
  !
  REAL(dp), PUBLIC :: zFadd_max, zFadd_min
  !
  ! fladd to store additional flux in the 9 WMO intervals, that is
  ! not changed; necessary to get the right value of pffrac.
  !
  REAL(dp), PUBLIC :: fladd

  !
  ! Chebyshev polynomial coeffs necessary to calculate O2 effective
  ! cross-sections: 
  !  ac, bc  : Chebyshev polynomial coefficients
  !  precis  : machine precision
  LOGICAL, PUBLIC              :: ldb_first=.FALSE.
  REAL(dp), PARAMETER , PUBLIC :: precis = 1.E-7_dp
  
  INTEGER, PUBLIC  :: ntr_max
  ! 
  INTEGER,  PARAMETER , PUBLIC :: tdim = 501
#ifndef LF
  REAL(dp), PARAMETER , PUBLIC :: t_del = 5._dp/REAL(tdim-1,dp)
  REAL(dp), PARAMETER , PUBLIC :: t_fac = REAL(tdim-1,dp)/5._dp
#else
  REAL(dp), PARAMETER , PUBLIC :: t_del = 5._dp/(tdim-1)
  REAL(dp), PARAMETER , PUBLIC :: t_fac = (tdim-1)/5._dp
#endif
  !
  TYPE sr_ctrl
    ! ntr  - truncation of Chebyshev polynom for Schumann-Runge bands
    ! type - type of the Schumann-Runge band parameterization
    CHARACTER(len=32) :: type = ''
    INTEGER           :: ntr
  END TYPE sr_ctrl
#ifdef LF
  PUBLIC :: sr_ctrl
#endif
  TYPE(sr_ctrl), PUBLIC :: sr
  ! 
  ! lsrflx = .TRUE. : Schumann-Runge HR calculated from flux profile
  LOGICAL, PUBLIC  :: lsrflx = .FALSE.

  ! # of spectral intervals of FUBRad
  INTEGER, PUBLIC :: nbands           ! CTRL namelist

  ! SOLAR CYCLE DATA LINEARLY INTERPOLATED IN TIME
  REAL(dp), PUBLIC        :: solfac = 0.5
  LOGICAL,  PUBLIC        :: lscale = .TRUE.

  ! type of resolution: DEFAULT - standard resolution in o3fluxes
  !                     DYNAMIC - for test with altered resolution
  CHARACTER(len=32), PUBLIC :: resolution = 'DEFAULT'
  !
  REAL(dp), DIMENSION(:), POINTER, PUBLIC :: val   => NULL()
  !
  ! INDEX POSITION OF PARAMETERS IN TABULAR ASCII FILE
  INTEGER, PARAMETER,    PUBLIC :: IDX_TSI = 1
  INTEGER, PARAMETER,    PUBLIC :: IDX_LYA = 2
  INTEGER, PARAMETER,    PUBLIC :: IDX_SR1 = 3
  INTEGER, PARAMETER,    PUBLIC :: IDX_SRS = 4
  INTEGER, PARAMETER,    PUBLIC :: IDX_SRL = 5
  INTEGER, DIMENSION(2), PUBLIC :: IDX_SRC = (/ 3,27/)
  INTEGER, DIMENSION(2), PUBLIC :: IDX_SRB = (/ 6, 6/)
  INTEGER, DIMENSION(2), PUBLIC :: IDX_HRZ = (/ 7,21/)
  INTEGER, DIMENSION(2), PUBLIC :: IDX_HAR = (/22,31/)
  INTEGER, DIMENSION(2), PUBLIC :: IDX_HUG = (/32,49/)

  INTEGER, PUBLIC               :: IDX_ADD = 50
  INTEGER, DIMENSION(2), PUBLIC :: IDX_CHA = (/50,50/)
  !
  ! Parameters specifying the number of intervals for different spectral regions.
  !
  INTEGER,            PUBLIC :: nsrb  = 1    ! Schumann-Runge bands
  INTEGER, PARAMETER, PUBLIC :: nsrb_km = 17 ! Sch.R. bands Koppers and Murtagh (1996)
  INTEGER, PARAMETER, PUBLIC :: nsrb_k  = 16 ! Sch.R. bands Kockarts (1994)
  INTEGER, PARAMETER, PUBLIC :: nsrb_st = 19 ! Sch.R. bands Strobel (1978)
  INTEGER, PARAMETER, PUBLIC :: nsrc  = 25   ! Schumann-Runge cont. Strobel (1978)
  INTEGER,            PUBLIC :: nherz = 15  ! Herzberg continuum
  INTEGER,            PUBLIC :: nhart = 10  ! Hartley bands
  INTEGER,            PUBLIC :: nhug  = 18  ! Huggins bands
  INTEGER,            PUBLIC :: nchap = 1   ! Chappuis bands

  ! ORBIT PARAMETER
  REAL(dp), PUBLIC :: cdisse_fubrad
 
END MODULE  messy_rad_fubrad_mem
