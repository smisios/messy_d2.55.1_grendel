!+ Data module for all data arrays, that are used by the radiation scheme
!-------------------------------------------------------------------------------

MODULE data_radiation

!-------------------------------------------------------------------------------
!
! Description:
!  This module declares and initializes all parametric data and data arrays 
!  that are used by the radiation scheme (i.e. the FESFT-routine). 
!  This data module replaces the COMMON blocks and local data arrays that are 
!  used in the original Fortran77 version of the FESFT-code 
!  provided by Bodo Ritter.
!
! Current Code Owner: DWD, Bodo Ritter
!  phone:  +49  69  8062 2703
!  fax:    +49  69  8062 3721
!  email:  Bodo.Ritter@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenter Doms
!  Initial release
! 3.6        2003/12/11 Ulrich Schaettler
!  Editorial changes (some lines have been too long)
! 3.21       2006/12/04 Ulrich Schaettler, Thorsten Reinhardt
!  New dimensions for running the radiation on a coarser grid
! V4_9         2009/07/16 Ulrich Schaettler
!  Added storage for slope of solar albedo with respect to soil water content
! V4_11        2009/11/30 Juergen Helmert
!  Definition of optical thickness variables for the different
!  possibilities according to setting of itype_aerosol
! V4_12        2010/05/11 Juergen Helmert
!  Unification of variables for optical thickness
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Changed the code owner
! V5_00_clm4   2015/06/05 Katherine Osterried ETHZ
!  Added new tuning parameter radfac
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
! Declarations:
!
! Modules used:
USE data_parameters , ONLY :   &
           ireals,   &! KIND-type parameters for real variables
           iintegers  ! kind-type parameter for "normal" integer variables

!-------------------------------------------------------------------------------

IMPLICIT NONE

! Local Declarations:

  INTEGER (KIND=iintegers), PRIVATE ::    &
    i ! index for data initalization

! Global (i.e. public) Declarations:

! Global Parameters
! -----------------

  INTEGER (KIND=iintegers), PARAMETER  ::    &

    ! Parameters for spectral resolution of radiative transfer code
    jpsol  = 3 ,     & ! Number of solar spectral intervals
    jpther = 5 ,     & ! Number of thermal spectral intervals
    jpspec = 8 ,     & ! (= jpsol+jpther) Total number of spectral intervals
 
    ! Parameters for gas absorbtion                                
    jpgas  = 3 ,     & ! Number of gases considered by radative transfer code
    jpabsc = 7         ! Maximum number of absorbtion coefficients in each   
                       ! spectral interval

! Global dimensions
! -----------------

  INTEGER (KIND=iintegers)             ::    &
    idim_rad,        & ! ie-dimension of the coarser grid
    istartrad,       & ! start- and end-indices for computing the radiation
    iendrad,         & !   (when running on a coarser grid, the input values for
    jstartrad,       & !    fesft are computed on all grid points, to compute an
    jendrad,         & !    average input over several grid points)
    iendparrad,      & ! end-index just for preparations
    jendparrad         ! end-index just for preparations

! Global Arrays and Scalars
! -------------------------

! 1. absorption properties of atmospheric gases (..,1=h2o; ..,2 =CO2; ..,3=O3)
! --------------------------------------------

  REAL  (KIND=ireals)     ::    &
    coai (jpabsc,jpspec,jpgas), & ! weigthing coefficients
    cobi (jpabsc,jpspec,jpgas), & ! absorption coefficients
    coali(jpabsc,jpspec,jpgas), & ! pressure correction coefficients
    cobti(jpabsc,jpspec,jpgas), & ! temperature correction coefficients
    pgas (       jpspec,jpgas), & ! reference pressure
    tgas (       jpspec,jpgas)    ! reference temperature

  INTEGER (KIND=iintegers)::    &
    ncgas(jpspec,3),            & ! number of coefficients for each spectral  
                                  ! interval and gas (maximum=7)
    nfast(jpspec)                 ! control variable for choice between 
                                  ! ESFT/FESFT in each spectral interval

  ! Initialization of above arrays

  DATA nfast / 1, 1, 1, 1, 1, 1, 1, 1/
                                                                        
  ! coefficients for H2O in spectral interval 1
  DATA ncgas(1,1) /7/ ; DATA pgas(1,1) /101325.000/ ; DATA tgas(1,1) /281.700/
  DATA (coai (i,1,1),i=1,7)/ .114190E-01,  .600200E-01,  .111201E+00,  .123340E+00,  .902500E-01,  .199632E+00,  .404139E+00/
  DATA (cobi (i,1,1),i=1,7)/ .209894E+02,  .208930E+01,  .184502E+00,  .217771E-01,  .279254E-02,  .463447E-03,  .000000E+00/
  DATA (coali(i,1,1),i=1,7)/ .285370E-01,  .688620E+00,  .766031E+00,  .833136E+00,  .773491E+00,  .768818E+00,  .100000E+01/
  DATA (cobti(i,1,1),i=1,7)/ .473006E+00, -.468639E+00, -.599601E+00, -.162223E+01, -.176002E+01, -.153131E+01,  .100000E+01/

  ! coefficients for H2O in spectral interval 2
  DATA ncgas(2,1) /7/ ; DATA pgas(2,1) /101325.000/ ; DATA tgas(2,1) /281.700/
  DATA (coai (i,2,1),i=1,7)/ .201500E-02,  .268530E-01,  .598920E-01,  .907740E-01,  .102284E+00,  .217298E+00,  .500884E+00/
  DATA (cobi (i,2,1),i=1,7)/ .508159E+01,  .519996E+00,  .465586E-01,  .891251E-02,  .159221E-02,  .374973E-03,  .000000E+00/
  DATA (coali(i,2,1),i=1,7)/-.482300E-02,  .529161E+00,  .587751E+00,  .756567E+00,  .774607E+00,  .733883E+00,  .100000E+01/
  DATA (cobti(i,2,1),i=1,7)/ .499755E+00, -.529716E+00, -.177970E-01, -.746447E+00, -.106191E+00, -.727589E+00,  .100000E+01/
                                                                        
  ! coefficients for H2O in spectral interval 3
  DATA ncgas(3,1) /3/ ; DATA pgas(3,1) /101325.000/ ; DATA tgas(3,1) /281.700/
  DATA (coai (i,3,1),i=1,7)/ .566900E-02,  .346720E-01,  .959659E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,3,1),i=1,7)/ .716144E-03,  .256449E-03,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,3,1),i=1,7)/-.281669E+00,  .611673E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,3,1),i=1,7)/ .418657E+00,  .405230E-01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for H2O in spectral interval 4
  DATA ncgas(4,1) /7/ ; DATA pgas(4,1) / 50662.500/ ; DATA tgas(4,1) /255.800/
  DATA (coai (i,4,1),i=1,7)/ .641200E-02,  .362630E-01,  .147064E+00,  .285387E+00,  .246376E+00,  .226899E+00,  .515980E-01/
  DATA (cobi (i,4,1),i=1,7)/ .298538E+04,  .139959E+03,  .152405E+02,  .144212E+01,  .183654E+00,  .283139E-01,  .409261E-02/
  DATA (coali(i,4,1),i=1,7)/ .183780E-01,  .410557E+00,  .808897E+00,  .897332E+00,  .932149E+00,  .978389E+00,  .100000E+01/
  DATA (cobti(i,4,1),i=1,7)/ .413777E+00, -.663704E+00, -.953789E+00, -.111883E+01, -.156269E+01, -.330557E+01,  .100000E+01/
                                                                        
  ! coefficients for H2O in spectral interval 5
  DATA ncgas(5,1) /7/ ; DATA pgas(5,1) / 86126.250/ ; DATA tgas(5,1) /281.700/
  DATA (coai (i,5,1),i=1,7)/ .147700E-02,  .345020E-01,  .865590E-01,  .144237E+00,  .218089E+00,  .339440E+00,  .175697E+00/
  DATA (cobi (i,5,1),i=1,7)/ .126765E+02,  .149624E+01,  .147571E+00,  .368129E-01,  .792501E-02,  .208930E-02,  .000000E+00/
  DATA (coali(i,5,1),i=1,7)/-.414300E-02,  .504464E+00,  .670985E+00,  .920940E+00,  .889089E+00,  .966028E+00,  .100000E+01/
  DATA (cobti(i,5,1),i=1,7)/ .454691E+00, -.423980E+01, -.340869E+01, -.410896E+01, -.268068E+01, -.250967E+01,  .100000E+01/

  ! coefficients for H2O in spectral interval 6
  DATA ncgas(6,1) /4/ ; DATA pgas(6,1) / 86126.250/ ; DATA tgas(6,1) /281.700/
  DATA (coai (i,6,1),i=1,7)/ .653200E-02,  .700040E-01,  .243768E+00,  .679696E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,6,1),i=1,7)/ .632412E+00,  .473151E-02,  .163305E-02,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,6,1),i=1,7)/ .794801E+00,  .306898E+00,  .100000E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,6,1),i=1,7)/-.100000E+02, -.219711E+01, -.369325E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for H2O in spectral interval 7
  DATA ncgas(7,1) /3/ ; DATA pgas(7,1) / 86126.250/ ; DATA tgas(7,1) /281.700/
  DATA (coai (i,7,1),i=1,7)/ .138610E-01,  .226595E+00,  .759544E+00, 0.000000E+00,  .000000E+00, 0.000000E+00, 0.000000E+00/
  DATA (cobi (i,7,1),i=1,7)/ .425598E-02,  .155239E-02, 0.000000E+00, 0.000000E+00,  .000000E+00, 0.000000E+00, 0.000000E+00/
  DATA (coali(i,7,1),i=1,7)/-.736171E+00,  .805828E+00,  .100000E+01, 0.000000E+00,  .000000E+00, 0.000000E+00, 0.000000E+00/
  DATA (cobti(i,7,1),i=1,7)/ .308301E+00, -.267573E+01,  .100000E+01, 0.000000E+00,  .000000E+00, 0.000000E+00, 0.000000E+00/
  ! coefficients for H2O in spectral interval 8
  DATA ncgas(8,1) /7/ ; DATA pgas(8,1) / 75993.750/ ; DATA tgas(8,1) /281.700/
  DATA (coai (i,8,1),i=1,7)/ .181840E-01,  .106586E+00,  .237611E+00,  .241085E+00,  .157304E+00,  .178767E+00,  .604640E-01/
  DATA (cobi (i,8,1),i=1,7)/ .822243E+02,  .979490E+01,  .905733E+00,  .140281E+00,  .193197E-01,  .320627E-02,  .000000E+00/
  DATA (coali(i,8,1),i=1,7)/-.126888E+00,  .701873E+00,  .834941E+00,  .920550E+00,  .849506E+00,  .931957E+00,  .100000E+01/
  DATA (cobti(i,8,1),i=1,7)/ .384580E+00, -.187972E+01, -.226834E+01, -.475940E+01, -.589531E+01, -.395962E+01,  .100000E+01/
                    
  ! coefficients for CO2 in spectral interval 1
  DATA ncgas(1,2) /6/ ; DATA pgas(1,2) / 86126.250/ ; DATA tgas(1,2) /255.800/
  DATA (coai (i,1,2),i=1,7)/ .592000E-02,  .667700E-02,  .423020E-01,  .732310E-01,  .140143E+00,  .731727E+00,  .000000E+00/
  DATA (cobi (i,1,2),i=1,7)/ .760326E+02,  .480839E+01,  .391742E+00,  .133968E-01,  .355631E-02,  .000000E+00,  .000000E+00/
  DATA (coali(i,1,2),i=1,7)/ .659071E+00,  .240858E+00,  .694157E+00,  .424843E+00,  .694262E+00,  .100000E+01,  .000000E+00/
  DATA (cobti(i,1,2),i=1,7)/ .467048E+00,  .395422E+00, -.902210E+00, -.557526E+00, -.473196E+00,  .100000E+01,  .000000E+00/
                                                                        
  ! coefficients for CO2 in spectral interval 2
  DATA ncgas(2,2) /3/ ; DATA pgas(2,2)/ 86126.250/ ;  DATA tgas(2,2) /255.800/
  DATA (coai (i,2,2),i=1,7)/ .278000E-02,  .197330E-01,  .977487E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,2,2),i=1,7)/ .169434E+00,  .103753E-01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,2,2),i=1,7)/ .138563E+00,  .831359E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,2,2),i=1,7)/ .475293E+00, -.496213E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for CO2 in spectral interval 3
  DATA ncgas(3,2) /2/ ; DATA pgas(3,2) / 86126.250/ ; DATA tgas(3,2) /255.800/
  DATA (coai (i,3,2),i=1,7)/ .306100E-02,  .996939E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,3,2),i=1,7)/ .101625E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,3,2),i=1,7)/ .100000E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,3,2),i=1,7)/-.100670E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for CO2 in spectral interval 4
  DATA ncgas(4,2) /0/ ; DATA pgas(4,2) / 60795.000/ ; DATA tgas(4,2) /255.800/
  DATA (coai (i,4,2),i=1,7)/ .335800E-02,  .996642E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/ 
  DATA (cobi (i,4,2),i=1,7)/ .247172E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,4,2),i=1,7)/ .100000E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,4,2),i=1,7)/-.807310E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for CO2 in spectral interval 5
  DATA ncgas(5,2) /7/ ; DATA pgas(5,2) / 10132.500/ ; DATA tgas(5,2) /229.900/
  DATA (coai (i,5,2),i=1,7)/ .452500E-02,  .321420E-01,  .659180E-01,  .101074E+00,  .107224E+00,  .186663E+00,  .502454E+00/
  DATA (cobi (i,5,2),i=1,7)/ .299226E+03,  .364754E+02,  .271644E+01,  .570164E+00,  .100231E+00,  .224388E-01,  .000000E+00/
  DATA (coali(i,5,2),i=1,7)/ .466819E+00,  .319510E+00,  .596734E+00,  .751216E+00,  .708519E+00,  .744381E+00,  .100000E+01/
  DATA (cobti(i,5,2),i=1,7)/ .358348E+00, -.739332E+00, -.183599E+01, -.289470E+01, -.214575E+01, -.585028E+01,  .100000E+01/
                                                                        
  ! coefficients for CO2 in spectral interval 6
  DATA ncgas(6,2) /3/ ; DATA pgas(6,2) / 50662.500/ ; DATA tgas(6,2) /255.800/
  DATA (coai (i,6,2),i=1,7)/ .119551E+00,  .899140E-01,  .790535E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,6,2),i=1,7)/ .305492E-02,  .148936E-02,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,6,2),i=1,7)/ .783365E+00, -.113116E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,6,2),i=1,7)/-.447333E+01,  .296352E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for CO2 in spectral interval 7
  DATA ncgas(7,2) /3/ ; DATA pgas(7,2) / 50662.500/ ; DATA tgas(7,2) /255.800/
  DATA (coai (i,7,2),i=1,7)/ .577890E-01,  .321750E-01,  .910036E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,7,2),i=1,7)/ .650130E-02,  .309030E-02,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,7,2),i=1,7)/ .295465E+00,  .930860E-01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,7,2),i=1,7)/-.562957E+01, -.984577E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for CO2 in spectral interval 8
  DATA ncgas(8,2) /4/ ; DATA pgas(8,2) / 50662.500/ ; DATA tgas(8,2) /255.800/
  DATA (coai (i,8,2),i=1,7)/ .317000E-02,  .127109E+00,  .114118E+00,  .755604E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,8,2),i=1,7)/ .174181E+02,  .495450E-01,  .165196E-01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,8,2),i=1,7)/ .511300E-02,  .252848E+00,  .851104E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,8,2),i=1,7)/ .495222E+00,  .445084E+00,  .117957E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for O3  in spectral interval 1
  DATA ncgas(1,3) /0/ ; DATA pgas(1,3) /  3039.75/ ; DATA tgas(1,3) /229.900/
  DATA (coai (i,1,3),i=1,7)/ .306000E-03,  .999694E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,1,3),i=1,7)/ .409261E+02,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,1,3),i=1,7)/ .618332E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,1,3),i=1,7)/-.974847E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for O3  in spectral interval 2
  DATA ncgas(2,3) /0/ ; DATA pgas(2,3) /  3039.75/ ; DATA tgas(2,3) /229.900/
  DATA (coai (i,2,3),i=1,7)/ .154800E-02,  .998452E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,2,3),i=1,7)/ .395367E+02,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,2,3),i=1,7)/ .592629E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,2,3),i=1,7)/-.106087E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for O3  in spectral interval 3
  DATA ncgas(3,3) /5/ ; DATA pgas(3,3) /  3039.75/ ; DATA tgas(3,3) /229.900/
  DATA (coai (i,3,3),i=1,7)/ .564000E-03,  .108690E-01,  .124320E-01,  .184417E+00,  .791718E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,3,3),i=1,7)/ .191426E+05,  .579429E+03,  .717794E+02,  .187068E+01,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,3,3),i=1,7)/-.204400E-02,  .776840E-01, -.229667E+00,  .994500E-01,  .100000E+01,  .000000E+00,  .000000E+00/
  DATA (cobti(i,3,3),i=1,7)/ .499912E+00,  .485463E+00,  .464581E+00, -.254634E+00,  .100000E+01,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for O3  in spectral interval 4
  DATA ncgas(4,3) /0/ ; DATA pgas(4,3) /  3039.75/ ; DATA tgas(4,3) /229.900/
  DATA (coai (i,4,3),i=1,7)/ .540000E-04,  .999946E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,4,3),i=1,7)/ .210378E+03,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,4,3),i=1,7)/ .490324E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,4,3),i=1,7)/ .500000E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for O3  in spectral interval 5
  DATA ncgas(5,3) /2/ ; DATA pgas(5,3) /  3039.75/ ; DATA tgas(5,3) /229.900/
  DATA (coai (i,5,3),i=1,7)/ .587700E-02,  .994123E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,5,3),i=1,7)/ .223357E+03,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,5,3),i=1,7)/ .551312E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,5,3),i=1,7)/-.140025E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for O3  in spectral interval 6
  DATA ncgas(6,3) /0/ ; DATA pgas(6,3) /  3039.75/ ; DATA tgas(6,3) /229.900/
  DATA (coai (i,6,3),i=1,7)/ .154100E-02,  .998459E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (i,6,3),i=1,7)/ .221820E+03,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(i,6,3),i=1,7)/ .546048E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(i,6,3),i=1,7)/-.273183E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
                                                                        
  ! coefficients for O3  in spectral interval 7
  DATA ncgas(7,3) /7/ ; DATA pgas(7,3) / 10132.50/ ; DATA tgas(7,3) /204.000/
  DATA (coai (i,7,3),i=1,7)/ .220500E-02,  .523500E-02,  .951500E-02,  .578800E-01,  .277389E+00,  .643850E-01,  .583391E+00/
  DATA (cobi (i,7,3),i=1,7)/ .434510E+03,  .299916E+03,  .121339E+03,  .827942E+02,  .157398E+02,  .615177E+01,  .000000E+00/
  DATA (coali(i,7,3),i=1,7)/ .224000E-03,  .100500E-02,  .571600E-02,  .508760E-01,  .524641E+00,  .896800E-01,  .100000E+01/
  DATA (cobti(i,7,3),i=1,7)/ .320370E+01,  .130031E+01, -.332851E+01,  .105177E+01, -.561714E+00, -.357670E+01,  .100000E+01/
                                                                        
  ! coefficients for O3  in spectral interval 8
  DATA ncgas(8,3) /0/ ; DATA pgas(8,3) /  3039.75/ ; DATA tgas(8,3) /229.900/
  DATA (coai (I,8,3),I=1,7)/ .397000E-03,  .999603E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobi (I,8,3),I=1,7)/ .230675E+03,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (coali(I,8,3),I=1,7)/ .564371E+00,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
  DATA (cobti(I,8,3),I=1,7)/-.479075E+01,  .100000E+01,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00,  .000000E+00/
 

! 2. Limits of spectral intervals in the radiation code (for information only)
! ------------------------------------------------------
 
  REAL  (KIND=ireals)   ::   &
    grenze(2,2,jpspec)       ! limits of spectral intervals jpspec

  !              WMIN1    WMAX1     WMIN2    WMAX2
  DATA grenze/  1.5300,   4.6420    ,999.    ,  0.    , &
                0.7000,   1.5300    ,999.    ,  0.    , &
                0.2451,   0.7000    ,999.    ,  0.    , &
               20.0000, 104.5150    ,999.    ,  0.    , &
               12.5000,  20.0000    ,999.    ,  0.    , &
                8.3333,   9.0090    , 10.3093, 12.5000, &
                9.0090,  10.3093    ,999.    ,  0.    , &
                4.6420,   8.3333    ,999.    ,  0.    /


! 3. Rayleigh scattering coefficients in solar spectral intervals
! --------------------------------------------------------------

  REAL  (KIND=ireals)   ::   &
    zrsc(jpsol)              ! Coefficients for Rayleigh scattering in solar spectral intervals

  DATA zrsc / 0.59776370E-08, 0.13266702E-06, 0.20634412E-05 /

 
! 4. Fraction of solar energy at TOA contained in individual solar spectral intervals 
!    based on data from LABS and NECKEL (1970/1984):
! --------------------------------------------------------------

  REAL  (KIND=ireals)   ::   &
    solant(jpsol)            ! Fraction of solar energy at TOA in individual spectral intervals

  DATA solant / 0.12888167, 0.41683156, 0.45428677 /

 
! 5. Coefficients for black body radiation and E-type coefficients                   
! --------------------------------------------------------------

  REAL  (KIND=ireals)   ::   &
      planck(3,jpther), & !
      zketypr (jpther), & !
      ztetypr (jpther), & !
      zketypa (jpther), & !
      ztetypa (jpther), & !
      zteref
      ! planck: coefficients for the description of the fraction of the total 
      !         black body radiation contained in an thermal spectral interval 
      !         as a function (2.order polynomial) of temperature:
      !
      !         F(T) = PLANCK(1,ISPEC) + PLANCK(2,ISPEC)*T + PLANCK(3,ISPEC)*T**2
      !
      ! zketyp: e-type continuum-coefficient for all spectral intervals 
      !         (PA (H2O)**-2) at 296 K
      !         (r)  following ROBERTS ET AL. 1976
      !         (a)  implicitly derived from the AFGL spectral data
      ! ztetyp: constant for the temperature dependancy of e-type   
      !         absorption for all intervals
      ! zteref: reference temperaure 
 
  DATA planck / 0.157656E+01, -0.711486E-02,  0.908220E-05, &
               -0.761337E-01,  0.339014E-02, -0.703246E-05, &
               -0.353624E+00,  0.321131E-02, -0.472513E-05, &
               -0.180726E+00,  0.148131E-02, -0.195189E-05, &
                0.343422E-01, -0.971936E-03,  0.463714E-05  /

  DATA zketypr / 0.0       , 0.418E-06 , 0.859E-06 , 0.594E-06 , 0.767E-07  / 
  DATA ztetypr / 0.0       , 0.0       , 1800.0    , 1800.     , 1800.      /
  DATA zketypa / 0.8426E-05, 0.8982E-06, 0.5489E-06, 0.4743E-06, 0.7040E-06 /
  DATA ztetypa / 1365.55   , 1544.38   , 1699.06   , 1724.39   , 1668.94    /
  DATA zteref  / 296.0     /
 
! 6. Aerosol optical properties for 8 spectral intervals
! ------------------------------------------------------
 
  REAL  (KIND=ireals)              ::           &
  zaea  (jpspec,5),& ! ratio of optical thickness for the absorption in spectral
                     ! interval jpspec  and total optical thickness at 0.55m*1.E-06
                     ! for an aerosoltyp specified by second array index
  zaes  (jpspec,5),& ! analog for the optical thickness of scattering
                     !
  zaeg  (jpspec,5),& ! factor of asymetry for specified aerosoltyp in spectral
                     ! interval jspec

  zaef  (jpspec,5)   ! forward scatterd fraction from aerosols. This array is 
                     ! initialized with 0 but modified later in opt_th and opt_so.
 
                     ! the following aerosoltyps (second array index) are considered:
                     ! 1 : continental
                     ! 2 : maritim
                     ! 3 : urban
                     ! 4 : vulcano ashes
                     ! 5 : stratosphaeric background aerosol

  ! The DATA statements for zaea, zaes, zaeg, zaef have been removed. These 
  ! variables are now set in init_radiation

 
! 7. Optical properties of liquid water for all 8 spectral intervals (solar and thermal)
! --------------------------------------------------------------------------------------
  ! These data-arrays are used localy in routines opt_so and opt_th

  ! For the calculation of the coefficients, spectral integration is performed 
  ! with NLA LSF and extinction = absorption + scattering is assumed.
  ! Single scattering albedo is used lateron.

  REAL  (KIND=ireals)  ::  &
    zlwe(4,jpspec), &  ! 
    zlww(2,jpspec), &  !
    zlwg(2,jpspec), &  !
    zlwemn(jpspec), &  ! minimum values of the extinction coefficients in Pa(H2O)**-1
    zlwemx(jpspec)     ! maximum values of the extinction coefficients in Pa(H2O)**-1

  DATA zlwe / -23.014052,     .173026,     .811865,     .000453, &
              -28.122596,     .172211,     .705673,     .000457, &
              -28.162592,     .198665,     .810637,     .000550, &
             -170.687770,     .498371,     .356225,     .001330, &
              -68.573703,     .263182,     .568143,     .000776, &
             -122.833213,     .297599,     .306486,     .000976, &
             -192.594948,     .440659,     .317142,     .001027, &
              -62.018469,     .281406,     .732715,     .000611  /

  DATA zlww /    .989679,  -22.291412, &
                 .999529,     .020875, &
                 .999999,     .000000, &
                 .302657,  102.711916, &
                 .337398,   80.596716, &
                 .449308,   52.823880, &
                 .686930,  -29.876242, &
                 .804203, -103.022685  /
  DATA zlwg /    .804992,   17.901033, &
                 .814785,   14.204375, &
                 .843955,    8.306586, &
                 .279400,  124.179115, &
                 .499491,  131.635343, &
                 .696708,   75.061613, &
                 .704732,   77.778408, &
                 .784672,   38.002913  /

  DATA zlwemn / 5.000,  5.000,  4.930,  5.800,  5.400,  5.200,  5.500,  5.5000 /
  DATA zlwemx /32.500, 32.500, 31.360, 18.600, 24.500, 18.200, 20.200, 32.4000 /


! 8. Optical properties of ice clouds for all 8 spectral intervals (solar and thermal)
! --------------------------------------------------------------------------------------
  ! These data-arrays are used localy in routines opt_so and opt_th

  ! The coefficients are derived using spectral averaging by weighted nonlinear LSF;
  ! Recombination of extinction, scattering and absorption after averaging of
  ! extinction = scattering + absorption

  REAL  (KIND=ireals)  ::  &
    ziwe(4,jpspec), &  ! 
    ziww(2,jpspec), &  ! coefficients for logarithmic fit
    ziwg(2,jpspec), &  ! coefficients for logarithmic fit
    ziwemn(jpspec), &  ! minimum values of the extinction coefficients in Pa(H2O)**-1
    ziwemx(jpspec)     ! maximum values of the extinction coefficients in Pa(H2O)**-1
 
  DATA ziwe / 16.726535,    0.007465,    1.354626,    0.000112, &
              17.531261,    0.003949,    0.669605,    0.000058, &
              17.698999,    0.003657,    0.625067,    0.000055, &
              19.592746,    0.008644,    1.153213,    0.000101, &
              18.990998,    0.006743,    0.997361,    0.000080, &
              18.482156,    0.004408,    0.693883,    0.000060, &
              18.603168,    0.005260,    0.813026,    0.000064, &
              18.437818,    0.004378,    0.692678,    0.000057  /

  DATA ziww /  0.694631,   -0.022160, &
               0.998669,   -0.000107, &
               0.999993,    0.000000, &
               0.289966,   -0.033855, &
               0.555820,    0.004491, &
               0.554495,    0.004904, &
               0.375319,   -0.017168, &
               0.485290,   -0.004358  /

  DATA ziwg /  0.976960,    0.007938, &
               0.914842,    0.003334, &
               0.900536,    0.001797, &
               1.134025,    0.032141, &
               1.053136,    0.015721, &
               1.010632,    0.006844, &
               1.063545,    0.011772, &
               1.035725,    0.008755  /

  DATA ziwemn/ 2.000, 2.000, 2.000, 2.000, 2.000, 2.000, 2.000, 2.0000 /
  DATA ziwemx/30.000,30.000,30.000,30.000,30.000,30.000,30.000, 30.000 /

! 9. Optical properties of ice clouds for all 8 spectral intervals (solar and thermal)
! --------------------------------------------------------------------------------------

  REAL  (KIND=ireals)  ::  &
    rad_csalbw(10)     !  slope of solar albedo with respect to soil water content
                       ! as a function of depth of upper soil layer
                       ! (has been computed every radiation time step before,
                       !  and is now computed in init_radiation)

! kos ETHZ, 2015/06/05 new tuning factor radfac
  REAL  (KIND=ireals)  ::  &
    radfac             ! fraction of cloud water/ice considered for radiation
! kos ETHZ end 
!=======================================================================

END MODULE data_radiation
