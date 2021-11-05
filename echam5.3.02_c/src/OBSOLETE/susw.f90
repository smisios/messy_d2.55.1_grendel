SUBROUTINE SUSW4

  !**** *SUSW*   - INITIALIZE COMMON MO_SW

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

  USE MO_KIND   , ONLY : DP

  USE MO_CONSTANTS, ONLY : RG=> G
  USE MO_SW     , ONLY : APAD     ,BPAD     ,RRAY     ,RSUN     ,&
       &       RPDH1    ,RPDU1    ,RPNH     ,RPNU     ,RSWCE    ,&
       &       RSWCP    ,RTDH2O   ,RTDUMG   ,RTH2O    ,RTUMG    ,&
       &       D

!!$  USE MO_SW     , ONLY : RROMA    ,RROMB    ,RRASY    ,RHSRA    ,&
!!$       &       RHSRB    ,RHSRC    ,RHSRD    ,RHSRE    ,RHSRF    ,&
!!$       &       RHSRTA   ,RHSRTB   ,RSNOALB  ,RSNOMEL  ,RWEIGS   ,&
!!$       &       RWEIGV

  IMPLICIT NONE


  !     ----------------------------------------------------------------
  REAL(DP):: ZAPAD4(4,3,7)  , ZBPAD4(4,3,7)  , ZD4(4,3)&
       & ,  ZRAY4(4,6)     , ZSUN4(4)       , ZSWCE4(4)  ,   ZSWCP4(4)

!!$  REAL(DP):: ZSNAL4(4)      , ZSNML4(4)      , ZWEIS4(4)  ,   ZWEIV4(4)&
!!$       & ,  ZROMA4(4)      , ZROMB4(4)      , ZRASY4(4)&
!!$       & ,  ZRA4(4)        , ZRB4(4)        , ZRC4(4)&
!!$       & ,  ZRD4(4)        , ZRE4(4)        , ZRF4(4)

  !     LOCAL INTEGER SCALARS
  INTEGER :: JC3, JC6, JI, JJ, JW, K

  !     LOCAL REAL SCALARS
  REAL(DP):: ZPDU4IS, ZPRU4IS, ZPDH4IS, ZPRH4IS &
       & ,  ZTDH4IS, ZTDU4IS, ZTH4IS,  ZTU4IS  &
       & ,  ZPDUMG,  ZPRUMG,  ZPDH2O,  ZPRH2O  &
       & ,  ZH2O,    ZUMG

!!$  REAL(DP):: ZRTO1,   ZRTO2

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

  ZSUN4(1) = 0.459760_DP

  ZD4(1,:)= (/ 0.000000000_DP, 0.000000000_DP, 0.000000000_DP /)

  !ZAPAD4(1,1:3,1:7) = RESHAPE((/&
  ! &0.184678379E+06_DP, 0.715303869E+01_DP, 0.925887084E-04_DP,&
  ! &0.553080884E+05_DP, 0.219386847E+03_DP, 0.129353723E-01_DP,&
  ! &0.248143712E+04_DP, 0.830001089E+03_DP, 0.800821928E+00_DP,&
  ! &0.000000000E-00_DP, 0.000000000E-00_DP, 0.242715973E+02_DP,&
  ! &0.000000000E-00_DP, 0.000000000E-00_DP, 0.878331486E+02_DP,&
  ! &0.000000000E-00_DP, 0.000000000E-00_DP, 0.191559725E+02_DP,&
  ! &0.000000000E-00_DP, 0.000000000E-00_DP, 0.000000000E+00_DP /)&
  ! &,(/3,7/))

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

  !ZBPAD4(1,1:3,1:7) = RESHAPE((/&
  ! &0.184678379E+06_DP, 0.715303869E+01_DP, 0.925887084E-04_DP,&
  ! &0.555188347E+05_DP, 0.219441875E+03_DP, 0.131812683E-01_DP,&
  ! &0.253257443E+04_DP, 0.831119997E+03_DP, 0.812706117E+00_DP,&
  ! &0.100000000E+01_DP, 0.100000000E+01_DP, 0.249863591E+02_DP,&
  ! &0.000000000E-00_DP, 0.000000000E-00_DP, 0.931071925E+02_DP,&
  ! &0.000000000E-00_DP, 0.000000000E-00_DP, 0.252233437E+02_DP,&
  ! &0.000000000E-00_DP, 0.000000000E-00_DP, 0.100000000E+01_DP /)&
  ! &,(/3,7/))

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

  ZSUN4(2) = 0.326158_DP

  ZD4(2,:)= (/ 0.000000000_DP, 0.000000000_DP, 1.000000000_DP /)

  !ZAPAD4(2,1:3,1:7) = RESHAPE((/&
  ! &0.690730834E-02_DP, 0.863790752E-03_DP, 0.000000000E+00_DP,&
  ! &0.151704275E+01_DP, 0.448762291E+00_DP, 0.000000000E+00_DP,&
  ! &0.751477543E+02_DP, 0.332530367E+02_DP, 0.000000000E+00_DP,&
  ! &0.759770236E+03_DP, 0.190914146E+03_DP, 0.000000000E+00_DP,&
  ! &0.109800326E+04_DP, 0.000000000E+00_DP, 0.000000000E+00_DP,&
  ! &0.148407574E+03_DP, 0.000000000E+00_DP, 0.000000000E+00_DP,&
  ! &0.000000000E+00_DP, 0.000000000E+00_DP, 0.000000000E+00_DP /)&
  ! &,(/3,7/))

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

  !ZBPAD4(2,1:3,1:7) = RESHAPE((/&
  ! &0.690730834E-02_DP, 0.863790752E-03_DP, 1.000000000E+00_DP,&
  ! &0.151954406E+01_DP, 0.448948107E+00_DP, 0.000000000E+00_DP,&
  ! &0.756512527E+02_DP, 0.333186750E+02_DP, 0.000000000E+00_DP,&
  ! &0.779384997E+03_DP, 0.192727216E+03_DP, 0.000000000E+00_DP,&
  ! &0.121113108E+04_DP, 0.100000000E+01_DP, 0.000000000E+00_DP,&
  ! &0.207678436E+03_DP, 0.000000000E+00_DP, 0.000000000E+00_DP,&
  ! &0.100000000E+01_DP, 0.000000000E+00_DP, 0.000000000E+00_DP /)&
  ! &,(/3,7/))

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

  ZSUN4(3) = 0.180608_DP

  ZD4(3,:)= (/ 0.000000000_DP, 0.000000000_DP, 1.000000000_DP /)

  !ZAPAD4(3,1:3,1:7) = RESHAPE((/&
  ! &0.837531303E-05_DP, 0.657978575E-02_DP, 0.000000000E+00_DP,&
  ! &0.173886341E-01_DP, 0.752617872E+00_DP, 0.000000000E+00_DP,&
  ! &0.518852799E+01_DP, 0.158209734E+02_DP, 0.000000000E+00_DP,&
  ! &0.159078416E+03_DP, 0.410274915E+02_DP, 0.000000000E+00_DP,&
  ! &0.493273523E+03_DP, 0.000000000E+00_DP, 0.000000000E+00_DP,&
  ! &0.102567293E+03_DP, 0.000000000E+00_DP, 0.000000000E+00_DP,&
  ! &0.000000000E+00_DP, 0.000000000E+00_DP, 0.000000000E+00_DP /)&
  ! &,(/3,7/))

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

  !ZBPAD4(3,1:3,1:7) = RESHAPE((/&
  ! &0.837531303E-05_DP, 0.657978575E-02_DP, 1.000000000E+00_DP,&
  ! &0.174882536E-01_DP, 0.753752065E+00_DP, 0.000000000E+00_DP,&
  ! &0.534536580E+01_DP, 0.159286262E+02_DP, 0.000000000E+00_DP,&
  ! &0.180351767E+03_DP, 0.424278450E+02_DP, 0.000000000E+00_DP,&
  ! &0.673126838E+03_DP, 0.100000000E+01_DP, 0.000000000E+00_DP,&
  ! &0.182718543E+03_DP, 0.000000000E+00_DP, 0.000000000E+00_DP,&
  ! &0.100000000E+01_DP, 0.000000000E+00_DP, 0.000000000E+00_DP /)&
  ! &,(/3,7/))

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

  ZSUN4(4) = 0.033474_DP

  ZD4(4,:)= (/ 0.000000000_DP, 0.000000000_DP, 0.000000000_DP /)

  !ZAPAD4(4,1:3,1:7) = RESHAPE((/&
  ! &0.122118185E-06_DP, 0.364233560E-10_DP, 0.263068898E+02_DP,&
  ! &0.154042531E-02_DP, 0.217340835E-06_DP, 0.146425875E+03_DP,&
  ! &0.141152193E+01_DP, 0.292623386E-03_DP, 0.860137809E+02_DP,&
  ! &0.685368761E+02_DP, 0.797100631E-01_DP, 0.000000000E+00_DP,&
  ! &0.216522281E+03_DP, 0.319103672E+01_DP, 0.000000000E+00_DP,&
  ! &0.421228746E+02_DP, 0.110530283E+02_DP, 0.000000000E+00_DP,&
  ! &0.000000000E+00_DP, 0.000000000E+00_DP, 0.000000000E+00_DP /)&
  ! &,(/3,7/))

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

  !ZBPAD4(4,1:3,1:7) = RESHAPE((/&
  ! &0.122118185E-06_DP, 0.364233560E-10_DP, 0.263068898E+02_DP,&
  ! &0.156287582E-02_DP, 0.218265806E-06_DP, 0.152569217E+03_DP,&
  ! &0.156994562E+01_DP, 0.297085917E-03_DP, 0.976791971E+02_DP,&
  ! &0.102304103E+03_DP, 0.834253403E-01_DP, 0.100000000E+01_DP,&
  ! &0.475479878E+03_DP, 0.358290002E+01_DP, 0.000000000E+00_DP,&
  ! &0.188714799E+03_DP, 0.139206164E+02_DP, 0.000000000E+00_DP,&
  ! &0.100000000E+01_DP, 0.100000000E+01_DP, 0.000000000E+00_DP /)&
  ! &,(/3,7/))

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

!!$  !=====================================================================
!!$
!!$  !*    2.2   SPECTRAL ALBEDO OF SNOW
!!$  !           Fresh and Melting Snow (Warren, 1982)
!!$
!!$  ZSNAL4(1:4)= (/ 0.920_DP, 0.798_DP , 0.159_DP, 0.010_DP /)
!!$  ZSNML4(1:4)= (/ 0.860_DP, 0.664_DP , 0.092_DP, 0.010_DP /)
!!$
!!$  !     ----------------------------------------------------------------
!!$
!!$  !*    2.3   WEIGHTS FOR SPECTRAL ALBEDO OF SOIL AND VEGETATION
!!$  !           a la Briegleb and Ramanathan (1982)
!!$
!!$  ZWEIS4(1:4)= (/ 1._DP, 2._DP, 2._DP, 1._DP /)
!!$  ZWEIV4(1:4)= (/ 1._DP, 4._DP, 2._DP, 1._DP /)
!!$
!!$  !     ----------------------------------------------------------------
!!$
!!$  !*    2.4   OPTICAL PARAMETERS FOR RAIN DROPS
!!$  !           Savijarvi et al. (1996)
!!$
!!$  ZRTO1 =  0.003_DP
!!$  ZRTO2 = -0.22_DP
!!$  ! CAUTION JUST TEMPORARY PARAMETERS      
!!$
!!$  ZROMA4(1:4)= (/ 0.00008_DP , 0.0105_DP , 0.264_DP  , 0.465_DP   /)
!!$  ZROMB4(1:4)= (/ 0.23_DP    , 0.22_DP   , 0.09_DP   , 0.001_DP   /)
!!$  ZRASY4(1:4)= (/ 0.88_DP    , 0.89_DP   , 0.94_DP   , 0.97_DP    /)
!!$
!!$  ZRA4(1:4)= (/ 1.5_DP     , 1.5_DP    , 1.5_DP    , 1.5_DP     /)
!!$  ZRB4(1:4)= (/ 0.50_DP    , 0.78_DP   , 1.13_DP   , 2.00_DP    /)
!!$  ZRC4(1:4)= (/ 5.58E-7_DP , 2.18E-5_DP, 8.55E-4_DP, 1.94E-1_DP /)
!!$  ZRD4(1:4)= (/ 1.25E-7_DP , 2.25E-5_DP, 1.28E-3_DP, 8.04E-3_DP /)
!!$  ZRE4(1:4)= (/ 0.841_DP   , 0.821_DP  , 0.786_DP  , 0.820_DP   /)
!!$  ZRF4(1:4)= (/ 2.08E-3_DP , 3.06E-3_DP, 5.32E-3_DP, 5.59E-3_DP /)
!!$
!!$
!!$  !=====================================================================

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

!!$  RHSRTA=ZRTO1
!!$  RHSRTB=ZRTO2
  DO JW=1,4
     RSUN (JW)=ZSUN4(JW)

     RSWCE(JW)=ZSWCE4(JW)
     RSWCP(JW)=ZSWCP4(JW)

!!$     RSNOALB(JW)=ZSNAL4(JW)
!!$     RSNOMEL(JW)=ZSNML4(JW)
!!$
!!$     RWEIGS(JW)=ZWEIS4(JW)
!!$     RWEIGV(JW)=ZWEIV4(JW)
!!$
!!$     RROMA(JW)=ZROMA4(JW)
!!$     RROMB(JW)=ZROMB4(JW)
!!$     RRASY(JW)=ZRASY4(JW)
!!$     RHSRA(JW)=ZRA4(JW)
!!$     RHSRB(JW)=ZRB4(JW)
!!$     RHSRC(JW)=ZRC4(JW)
!!$     RHSRD(JW)=ZRD4(JW)
!!$     RHSRE(JW)=ZRE4(JW)
!!$     RHSRF(JW)=ZRF4(JW)

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

END SUBROUTINE SUSW4