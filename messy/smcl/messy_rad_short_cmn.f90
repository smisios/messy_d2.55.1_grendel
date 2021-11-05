! *************************************************************************
MODULE messy_rad_short_cmn
! *************************************************************************

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PUBLIC
  SAVE

  INTEGER, PARAMETER :: NSW=4 ! number of short wave bands

  !  NOVLP    : index for cloud overlap assumption in radiation computation
  !             1 : maximum-random overlap
  !             2 : maximum overlap
  !             3 : random overlap
  !
  INTEGER, PARAMETER :: NOVLP=1

  ! SCALING FACTORS FOR THE SHORTWAVE RADIATION PARAMETERISATION.
  REAL(DP), DIMENSION(NSW) :: &
       rsun_scale = (/0.459760_DP, & ! 0.25 - 0.69 Microns
                      0.326158_DP, & ! 0.69 - 1.19 Microns
                      0.180608_DP, & ! 1.19 - 2.38 Microns
                      0.033474_DP  & ! 2.38 - 4.00 Microns
                     /)

  ! --------------------------------------------------------------------------
  ! MODULE WORK SPACE
  !
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
  !
  
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


  ! --------------------------------------------------------------------------

  ! SUBROUTINES
  !
  !PUBLIC :: rad_sw_initialize

CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE rad_sw_initialize
    !SUBROUTINE SUSW4

    !**** *SUSW*   - INITIALIZE messy_rad_short

    !     PURPOSE.
    !     --------
    !           INITIALIZE MO_SW, THE COMMON THAT CONTAINS COEFFICIENTS
    !           NEEDED TO RUN THE SHORTWAVE RADIATION SUBROUTINES

    !**   INTERFACE.
    !     ----------
    !        *CALL* *SUSW

    !        EXPLICIT ARGUMENTS :
    !        --------------------
    !        NONE

    !        IMPLICIT ARGUMENTS :
    !        --------------------
    !        COMMON MO_SW

    !     METHOD.
    !     -------
    !        SEE DOCUMENTATION

    !     EXTERNALS.
    !     ----------

    !     REFERENCE.
    !     ----------
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 88-12-15
    !        97-04-16 JJ Morcrette  2 and 4 interval spectral resolution
    !        M.A. Giorgetta, MPI, June 2000:
    !        - simplified to 4 bands only,
    !        - setup code not relevant for ECHAM is uncommented by !!$

    !     ------------------------------------------------------------------

    USE messy_main_constants_mem, ONLY : RG=> G

    IMPLICIT NONE

    !     ----------------------------------------------------------------
    REAL(DP):: ZAPAD4(4,3,7)  , ZBPAD4(4,3,7)  , ZD4(4,3)&
         & ,  ZRAY4(4,6)     , ZSUN4(4)       , ZSWCE4(4)  ,   ZSWCP4(4)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JC3, JC6, JI, JJ, JW

    !     LOCAL REAL SCALARS
    REAL(DP):: ZPDU4IS, ZPRU4IS, ZPDH4IS, ZPRH4IS &
         & ,  ZTDH4IS, ZTDU4IS, ZTH4IS,  ZTU4IS  &
         & ,  ZPDUMG,  ZPRUMG,  ZPDH2O,  ZPRH2O  &
         & ,  ZH2O,    ZUMG

    !     ----------------------------------------------------------------

    !*        1.  CLEAR-SKY ABSORPTION COEFFICIENTS FOR N SPECTRAL INTERVALS
    !             --------------------------------------------------------


    !*        1.2  COEFFICIENTS FOR FOUR SPECTRAL INTERVALS
    !              ----------------------------------------


    !* DERIVED FROM HITRAN APRIL 1992 with LOWTRAN P AND T SCALING
    !       H2O:  Pref=1000hPa, Tref=296K, Pdep=0.9
    !       UMG:  Pref=1000hPa, Tref=296K, Pdep=0.75 (CO2+N2O+CO+CH4+O2)
    !       O3 :  unchanged in interval 1, from HITRAN 92 in interval 4

    ZTDH4IS = 0.450_DP
    ZTDU4IS = 0.375_DP
    ZTH4IS  = 296._DP
    ZTU4IS  = 296._DP
    ZPDH4IS = 0.90_DP
    ZPDU4IS = 0.75_DP
    ZPRH4IS = 100000._DP
    ZPRU4IS = 100000._DP

    !* 1st spectral interval: U.V. and Visible (0.25 - 0.69 Micron)

    ZSUN4(1) = rsun_scale(1)

    ZD4(1,:)= (/ 0.000000000_DP, 0.000000000_DP, 0.000000000_DP /)

    ZAPAD4(1, 1, :) = (/&
         &0.184678379E+06_DP,&
         &0.553080884E+05_DP,&
         &0.248143712E+04_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP/)
    ZAPAD4(1, 2, :) = (/&
         &0.715303869E+01_DP,&
         &0.219386847E+03_DP,&
         &0.830001089E+03_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP/)
    ZAPAD4(1, 3, :) = (/&
         &0.925887084E-04_DP,&
         &0.129353723E-01_DP,&
         &0.800821928E+00_DP,&
         &0.242715973E+02_DP,&
         &0.878331486E+02_DP,&
         &0.191559725E+02_DP,&
         &0.000000000E+00_DP/)

    ZBPAD4(1, 1, :) = (/&
         &0.184678379E+06_DP,&
         &0.555188347E+05_DP,&
         &0.253257443E+04_DP,&
         &0.100000000E+01_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP/)
    ZBPAD4(1, 2, :) = (/&
         &0.715303869E+01_DP,&
         &0.219441875E+03_DP,&
         &0.831119997E+03_DP,&
         &0.100000000E+01_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP,&
         &0.000000000E-00_DP/)
    ZBPAD4(1, 3, :) = (/&
         &0.925887084E-04_DP,&
         &0.131812683E-01_DP,&
         &0.812706117E+00_DP,&
         &0.249863591E+02_DP,&
         &0.931071925E+02_DP,&
         &0.252233437E+02_DP,&
         &0.100000000E+01_DP/)

    ZRAY4(1,:)= (/&
         &.428937E-01_DP, .890743E+00_DP,-.288555E+01_DP,&
         &.522744E+01_DP,-.469173E+01_DP, .161645E+01_DP/)

    ZSWCE4(1) = 0._DP
    ZSWCP4(1) = 0._DP

    !     ----------------------------------------------------------------

    !* Near-Infrared (0.69 - 4.0 Microns) is sub-divided into:

    !     ----------------------------------------------------------------

    !* 0.69 - 1.19 Micron

    ZSUN4(2) = rsun_scale(2)

    ZD4(2,:)= (/ 0.000000000_DP, 0.000000000_DP, 1.000000000_DP /)

    ZAPAD4(2, 1, :) = (/&
         &0.690730834E-02_DP,&
         &0.151704275E+01_DP,&
         &0.751477543E+02_DP,&
         &0.759770236E+03_DP,&
         &0.109800326E+04_DP,&
         &0.148407574E+03_DP,&
         &0.000000000E+00_DP/)
    ZAPAD4(2, 2, :) = (/&
         &0.863790752E-03_DP,&
         &0.448762291E+00_DP,&
         &0.332530367E+02_DP,&
         &0.190914146E+03_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)
    ZAPAD4(2, 3, :) = (/&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)

    ZBPAD4(2, 1, :) = (/&
         &0.690730834E-02_DP,&
         &0.151954406E+01_DP,&
         &0.756512527E+02_DP,&
         &0.779384997E+03_DP,&
         &0.121113108E+04_DP,&
         &0.207678436E+03_DP,&
         &0.100000000E+01_DP/)
    ZBPAD4(2, 2, :) = (/&
         &0.863790752E-03_DP,&
         &0.448948107E+00_DP,&
         &0.333186750E+02_DP,&
         &0.192727216E+03_DP,&
         &0.100000000E+01_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)
    ZBPAD4(2, 3, :) = (/&
         &1.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)

    ZRAY4(2,:)= (/&
         &.164261E-01_DP, .000000E+00_DP, .000000E+00_DP,&
         &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP/)

    ZSWCE4(2) = 0._DP
    ZSWCP4(2) = 0._DP

    !     ----------------------------------------------------------------

    !* 1.19 - 2.38 Microns

    ZSUN4(3) = rsun_scale(3)

    ZD4(3,:)= (/ 0.000000000_DP, 0.000000000_DP, 1.000000000_DP /)

    ZAPAD4(3, 1, :) = (/&
         &0.837531303E-05_DP,&
         &0.173886341E-01_DP,&
         &0.518852799E+01_DP,&
         &0.159078416E+03_DP,&
         &0.493273523E+03_DP,&
         &0.102567293E+03_DP,&
         &0.000000000E+00_DP/)
    ZAPAD4(3, 2, :) = (/&
         &0.657978575E-02_DP,&
         &0.752617872E+00_DP,&
         &0.158209734E+02_DP,&
         &0.410274915E+02_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)
    ZAPAD4(3, 3, :) = (/&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)

    ZBPAD4(3, 1, :) = (/&
         &0.837531303E-05_DP,&
         &0.174882536E-01_DP,&
         &0.534536580E+01_DP,&
         &0.180351767E+03_DP,&
         &0.673126838E+03_DP,&
         &0.182718543E+03_DP,&
         &0.100000000E+01_DP/)
    ZBPAD4(3, 2, :) = (/&
         &0.657978575E-02_DP,&
         &0.753752065E+00_DP,&
         &0.159286262E+02_DP,&
         &0.424278450E+02_DP,&
         &0.100000000E+01_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)
    ZBPAD4(3, 3, :) = (/&
         &1.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)

    ZRAY4(3,:)= (/&
         &.180438E-02_DP, .000000E+00_DP, .000000E+00_DP,&
         &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP/)

    ZSWCE4(3) = 0._DP
    ZSWCP4(3) = 0._DP

    !     ----------------------------------------------------------------

    !* 2.38 - 4.00 Microns

    ZSUN4(4) = rsun_scale(4)

    ZD4(4,:)= (/ 0.000000000_DP, 0.000000000_DP, 0.000000000_DP /)

    ZAPAD4(4, 1, :) = (/&
         &0.122118185E-06_DP,&
         &0.154042531E-02_DP,&
         &0.141152193E+01_DP,&
         &0.685368761E+02_DP,&
         &0.216522281E+03_DP,&
         &0.421228746E+02_DP,&
         &0.000000000E+00_DP/)
    ZAPAD4(4, 2, :) = (/&
         &0.364233560E-10_DP,&
         &0.217340835E-06_DP,&
         &0.292623386E-03_DP,&
         &0.797100631E-01_DP,&
         &0.319103672E+01_DP,&
         &0.110530283E+02_DP,&
         &0.000000000E+00_DP/)
    ZAPAD4(4, 3, :) = (/&
         &0.263068898E+02_DP,&
         &0.146425875E+03_DP,&
         &0.860137809E+02_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)

    ZBPAD4(4, 1, :) = (/&
         &0.122118185E-06_DP,&
         &0.156287582E-02_DP,&
         &0.156994562E+01_DP,&
         &0.102304103E+03_DP,&
         &0.475479878E+03_DP,&
         &0.188714799E+03_DP,&
         &0.100000000E+01_DP/)
    ZBPAD4(4, 2, :) = (/&
         &0.364233560E-10_DP,&
         &0.218265806E-06_DP,&
         &0.297085917E-03_DP,&
         &0.834253403E-01_DP,&
         &0.358290002E+01_DP,&
         &0.139206164E+02_DP,&
         &0.100000000E+01_DP/)
    ZBPAD4(4, 3, :) = (/&
         &0.263068898E+02_DP,&
         &0.152569217E+03_DP,&
         &0.976791971E+02_DP,&
         &0.100000000E+01_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP,&
         &0.000000000E+00_DP/)

    ZRAY4(4,:)= (/&
         &.136020E-03_DP, .000000E+00_DP, .000000E+00_DP,&
         &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP/)

    ZSWCE4(4) = 0._DP
    ZSWCP4(4) = 0._DP

    !*       2.    SET VALUES.
    !              -----------


    ZPDH2O = ZPDH4IS
    ZPDUMG = ZPDU4IS
    ZPRH2O = ZPRH4IS
    ZPRUMG = ZPRU4IS
    RTDH2O = ZTDH4IS
    RTDUMG = ZTDU4IS
    RTH2O  = ZTH4IS
    RTUMG  = ZTU4IS

    RPDH1=ZPDH2O+1._DP
    RPDU1=ZPDUMG+1._DP
    ZH2O=1._DP/( 10._DP* RG * RPDH1 )
    ZUMG=1._DP/( 10._DP* RG * RPDU1 )
    RPNU = ZUMG/(ZPRUMG**ZPDUMG)
    RPNH = ZH2O/(ZPRH2O**ZPDH2O)

    DO JW=1,4
       RSUN (JW)=ZSUN4(JW)

       RSWCE(JW)=ZSWCE4(JW)
       RSWCP(JW)=ZSWCP4(JW)

       DO JC3=1,3
          D(JW,JC3)=ZD4(JW,JC3)
       ENDDO
       DO JC6=1,6
          RRAY(JW,JC6)=ZRAY4(JW,JC6)
       ENDDO
       DO JI=1,3
          DO JJ=1,7
             APAD(JW,JI,JJ)=ZAPAD4(JW,JI,JJ)
             BPAD(JW,JI,JJ)=ZBPAD4(JW,JI,JJ)
          ENDDO
       ENDDO
    ENDDO

    !END SUBROUTINE SUSW4
  END SUBROUTINE rad_sw_initialize
  ! ------------------------------------------------------------------

! ***************************************************************************
END MODULE messy_rad_short_cmn
! ***************************************************************************
