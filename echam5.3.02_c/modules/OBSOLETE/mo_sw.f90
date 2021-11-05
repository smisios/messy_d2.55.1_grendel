MODULE MO_SW


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

! PARAMETERS:
!
!  nsw      : number of bands of SW scheme
!
!  novlp    : index for cloud overlap assumption in radiation computation
!             1 : maximum-random overlap
!             2 : maximum overlap
!             3 : random overlap
!

INTEGER, PARAMETER :: NSW=4
INTEGER, PARAMETER :: NOVLP=1


!     ------------------------------------------------------------------
!*    ** *MO_SW* - COEFFICIENTS FOR SHORTWAVE RADIATION TRANSFER
!     ------------------------------------------------------------------

REAL(DP):: APAD(4,3,7)
REAL(DP):: BPAD(4,3,7)
REAL(DP):: RRAY(4,6)
REAL(DP):: RSUN(4)
REAL(DP):: RPDH1
REAL(DP):: RPDU1
REAL(DP):: RPNH
REAL(DP):: RPNU
REAL(DP):: RSWCE(4)
REAL(DP):: RSWCP(4)
REAL(DP):: RTDH2O
REAL(DP):: RTDUMG
REAL(DP):: RTH2O
REAL(DP):: RTUMG
REAL(DP):: D(4,3)

REAL(DP):: RYFWCA(4)
REAL(DP):: RYFWCB(4)
REAL(DP):: RYFWCC(4)
REAL(DP):: RYFWCD(4)
REAL(DP):: RYFWCE(4)
REAL(DP):: RYFWCF(4)

REAL(DP):: REBCUA(4)
REAL(DP):: REBCUB(4)
REAL(DP):: REBCUC(4)
REAL(DP):: REBCUD(4)
REAL(DP):: REBCUE(4)
REAL(DP):: REBCUF(4)
REAL(DP):: REBCUG(16)
REAL(DP):: REBCUH(16)
REAL(DP):: REBCUI(6)
REAL(DP):: REBCUJ(6)

REAL(DP):: RASWCA(4)
REAL(DP):: RASWCB(4)
REAL(DP):: RASWCC(4)
REAL(DP):: RASWCD(4)
REAL(DP):: RASWCE(4)
REAL(DP):: RASWCF(4)

REAL(DP):: RFULIO(16,3)
REAL(DP):: RFLAA0(4)
REAL(DP):: RFLAA1(4)
REAL(DP):: RFLBB0(4)
REAL(DP):: RFLBB1(4)
REAL(DP):: RFLBB2(4)
REAL(DP):: RFLBB3(4)
REAL(DP):: RFLCC0(4)
REAL(DP):: RFLCC1(4)
REAL(DP):: RFLCC2(4)
REAL(DP):: RFLCC3(4)
REAL(DP):: RFLDD0(4)
REAL(DP):: RFLDD1(4)
REAL(DP):: RFLDD2(4)
REAL(DP):: RFLDD3(4)

REAL(DP):: RHSAVI(16,3)

REAL(DP):: RSUSHE(4)
REAL(DP):: RSUSHF(4)
REAL(DP):: RSUSHH(4)
REAL(DP):: RSUSHK(4)
REAL(DP):: RSUSHA(4)
REAL(DP):: RSUSHG(4)
REAL(DP):: RSUSHFA(4)
REAL(DP):: RSUSHC
REAL(DP):: RSUSHD

REAL(DP):: REFFIA
REAL(DP):: REFFIB
REAL(DP):: RTIW
REAL(DP):: RRIW
REAL(DP):: RROMA(4)
REAL(DP):: RROMB(4)
REAL(DP):: RRASY(4)

REAL(DP):: RHSRA(4)
REAL(DP):: RHSRB(4)
REAL(DP):: RHSRC(4)
REAL(DP):: RHSRD(4)
REAL(DP):: RHSRE(4)
REAL(DP):: RHSRF(4)
REAL(DP):: RHSRTA
REAL(DP):: RHSRTB

REAL(DP):: RTAUA(4,6)
REAL(DP):: RPIZA(4,6)
REAL(DP):: RCGA(4,6)
REAL(DP):: RAER(6,6)
REAL(DP):: RSNOALB(4)
REAL(DP):: RSNOMEL(4)
REAL(DP):: RWEIGS(4)
REAL(DP):: RWEIGV(4)


!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
!  APAD  :  REAL     PADE APPROXIMANTS NUMERATOR
!  BPAD  :  REAL     PADE APPROXIMANTS DENOMINATOR
!  D     :  REAL     TRANSMISSION LIMIT FOR INFINITE ABSORBER AMOUNT
!  RRAY  :  REAL     RAYLEIGH SCATTERING COEFFICIENTS
!  RSUN  :  REAL     SOLAR FRACTION IN SPECTRAL INTERVALS
!  RPDH1 :  1 + EXPONENT PRESSURE DEPENDENCE H2O
!  RPDU1 :  1 + EXPONENT PRESSURE DEPENDENCE UNIFORMLY MIXED GASES
!  RPNH  :  REFERENCE PRESSURE FACTOR FOR H2O
!  RPNU  :  REFERENCE PRESSURE FACTOR FOR UNIFORMLY MIXED GASES
!  RSWCE :  E-TYPE, H2O CONTINUUM ABSORPTION COEFFICIENT 
!  RSWCP :  P-TYPE, H2O CONTINUUM ABSORPTION COEFFICIENT 
!  RTDH2O:  EXPONENT TEMPERATURE DEPENDENCE H2O
!  RTDUMG:  EXPONENT TEMPERATURE DEPENDENCE UNIFORMLY MIXED GASES
!  RTH2O :  REFERENCE TEMPERATURE H2O
!  RTUMG :  REFERENCE TEMPERATURE UNIFORMLY MIXED GASES
!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
!*    FOUQUART (1987) WATER CLOUD OPTICAL PROPERTIES

! RYFWCA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! RYFWCB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! RYFWCC :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCD :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCE :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCF :  REAL   : ASSYMETRY FACTOR

!*    SLINGO (1989) WATER CLOUD OPTICAL PROPERTIES

! RASWCA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! RASWCB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! RASWCC :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCD :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCE :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCF :  REAL   : ASSYMETRY FACTOR

!*   SAVIJARVI (1998) WATER CLOUD OPTICAL PROPERTIES (RRTM)

! RHSAVI : REAL    : MASS ABSORPTION COEFFICIENTS (POLYNOMIAL DEVELOPM)

!*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM EBERT-CURRY (1992)

! REBCUA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! REBCUB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! REBCUC :  REAL   : 1-C3  IN SINGLE SCATTERING ALBEDO FORMULA
! REBCUD :  REAL   : C4 IN SINGLE SCATTERING ALBEDO FORMULA
! REBCUE :  REAL   : C5 IN ASSYMETRY FACTOR FORMULA
! REBCUF :  REAL   : C6 IN ASSYMETRY FACTOR FORMULA
! REBCUG :  REAL   : C7 IN MASS ABSORPTION COEFFICIENT FORMULA
! REBCUH :  REAL   : C8 IN MASS ABSORPTION COEFFICIENT FORMULA
! REBCUI :  REAL   : C7 IN MASS ABSORPTION COEFFICIENT SPECTRAL FORMULA
! REBCUJ :  REAL   : C8 IN MASS ABSORPTION COEFFICIENT SPECTRAL FORMULA

!*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM SUN-SHINE (1995)

! RSHSUE :  REAL   : E IN SINGLE SCATTERING ALBEDO FORMULA
! RSHSUF :  REAL   : F IN SINGLE SCATTERING ALBEDO FORMULA
! RSHSUH :  REAL   : H IN ASSYMETRY FACTOR FORMULA
! RSHSUK :  REAL   : K IN ASSYMETRY FACTOR FORMULA
! RSHSUA :  REAL   : ALPHA IN SSA CORRECTION FACTOR FORMULA
! RSHSUG :  REAL   : GAMMA IN ASSYMETRY CORRECTION FACTOR FORMULA
! RSHSUFA:  REAL   : COEFFICIENTS IN TEMPERATURE CORRECTION FACTOR

! REFFIA :  REAL   : C9  IN EFFECTIVE RADIUS FORMULA
! REFFIB :  REAL   : C10 IN EFFECTIVE RADIUS FORMULA

!*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM FU-LIOU (1993)

! RFULIO :  REAL   : COEFFICIENTS IN EXPRESSION FOR LW EXTINCTION COEFF.
! RFLAA  :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW EXTINCTION COEFF.
! RFLBB  :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW SINGLE SCATT.ALB.
! RFLCC  :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW ASSYMETRY FACTOR
! RFLDD  :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW ASSYMETRY FACTOR

!*    TRANSITION BETWEEN LIQUID AND SOLID WATER

! RTIW   :  REAL   : TEMPERATURE THRESHOLD
! RRIW   :  REAL   : TRANSITION RANGE

!*    RAIN OPTICAL PROPERTIES FROM SAVIJARVI (1996)

! RROMA  :  REAL   : COEFFICIENTS FOR SINGLE SCATTERING ALBEDO
! RROMB  :  REAL   : COEFFICIENTS FOR SINGLE SCATTERING ALBEDO
! RRASY  :  REAL   : COEFFICIENTS FOR ASSYMETRY FACTOR
! RHSRA  :  REAL   : COEFFICIENTS FOR OPTICAL THICKNESS
! RHSRB  :  REAL   : COEFFICIENTS FOR OPTICAL THICKNESS
! RHSRC  :  REAL   : COEFFICIENTS FOR SINGLE SCATTERING ALBEDO
! RHSRD  :  REAL   : COEFFICIENTS FOR SINGLE SCATTERING ALBEDO
! RHSRE  :  REAL   : COEFFICIENTS FOR ASSYMETRY FACTOR 
! RHSRF  :  REAL   : COEFFICIENTS FOR ASSYMETRY FACTOR
! RHSRTA :  REAL   : COEFFICIENTS FOR OPTICAL THICKNESS
! RHSRTB :  REAL   : COEFFICIENTS FOR OPTICAL THICKNESS
!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : -------
!  RTAUA :  REAL     S.W. NORMALIZED OPTICAL THICKNESS AT 0.55 MICRON
!  RPIZA :  REAL     S.W. SINGLE SCATTERING ALBEDO
!  RCGA  :  REAL     S.W. ASSYMETRY FACTOR
!  RAER  :  REAL     L.W. ABSORPTION COEFFICIENTS
!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : -------
! RSNOALB:  REAL     S.W. SPECTRAL ALBEDO (Fresh Snow) after WARREN
! RSNOMEL:  REAL     S.W. SPECTRAL ALBEDO (Aging Snow) after WARREN

! RWEIGS :  REAL     S.W. SPECTR WEIGHT for soil (Briegleb, Ramanathan)
! RWEIGV :  REAL     S.W. SPECTR WEIGHT for vegetation (BR86)
!     -----------------------------------------------------------------
END MODULE MO_SW
