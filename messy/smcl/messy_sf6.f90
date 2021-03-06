MODULE messy_sf6

  USE messy_main_constants_mem, ONLY: &
    dp, &                                 ! kind parameter for real
    Pi

  IMPLICIT NONE

!-----------------------------------------------------------------
! Everything is PRIVATE, except when explicitely stated otherwise:
!-----------------------------------------------------------------
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modver = '2.1'
  CHARACTER(LEN=*), PARAMETER :: modstr = 'sf6'
 
  PUBLIC :: sf6_read_nml_ctrl 
  PUBLIC :: calc_age
  PUBLIC :: sf6life
  PUBLIC :: rnelec

  PUBLIC :: dp, modver, modstr

  REAL(dp), SAVE :: nml_t0
  REAL(dp), SAVE :: nml_a_t0
  REAL(dp), SAVE :: nml_b

  INTEGER,  SAVE :: SF6_tracer, tp_gap, emode, idontcalc
  REAL(dp), SAVE :: autocoeff
  
  PUBLIC :: nml_t0, nml_a_t0, nml_b
  PUBLIC :: SF6_tracer, tp_gap, emode, autocoeff, idontcalc

CONTAINS

  SUBROUTINE sf6_read_nml_ctrl(status, iou)
  !--------------------------------------------------------------------
  ! This routine reads and checks the ageofair CTRL namelist. It is designed
  ! according to the MESSy (Mainz Earth Submodel System) standard. 
  !--------------------------------------------------------------------

  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  USE messy_main_blather, ONLY: start_message, end_message
  
  IMPLICIT NONE
  
  !-----
  ! input/output variables
  !-----
  INTEGER, INTENT(out) :: status   ! error status
  INTEGER, INTENT(in) :: iou       ! logical I/O unit
  
  !-----
  ! local variables
  !-----
  CHARACTER(len=*), PARAMETER :: substr = 'sf6_read_nml_ctrl'
  LOGICAL                     :: lex     ! file existence flag
  INTEGER                     :: fstat   ! file status

  NAMELIST /CTRL/ SF6_tracer,  nml_t0, nml_a_t0, nml_b, tp_gap, emode,autocoeff, &
    &  idontcalc
  
  CALL start_message(TRIM(modstr),'INITIALISATION', substr)
  !-----
  ! initialisation
  !-----
  status = 1

  CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
  IF (.not.lex) RETURN   ! SF6.nml does not exist

  READ(iou, NML=CTRL, IOSTAT=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
  IF (fstat /= 0) RETURN   ! error while reading namelist

  write(*,*) ' /--------------------------------\ '
  write(*,*) ' | CALCULATION OF MEAN AGE OF AIR | '
  write(*,*) ' \--------------------------------/ '
  write(*,*) ' '
!  write(*,*) '----------------------------------'
!  write(*,*) ' C A U T I O N ! ! ! ! ! ! !'
!  write(*,*) ' special version only for testing'
!  write(*,*) '----------------------------------'
!  write(*,*) ' '
  write(*,*) ' purpose: compute the life time of sf6 including realistic processes:'
  write(*,*) '  rk1 lambda(phot, sf6 -> sfi)'
  write(*,*) '  rk2 lambda(e-,   sf6 -> sf6-*)'
  write(*,*) '  rk3 lambda(photo, sf6- -> sf6)'
  write(*,*) '  rk4 lambda(hydro, sf6- -> hf + ..)'
  write(*,*) '  rk5 lambda(M, sf-* -> sf6)'
  write(*,*) '  rk6 lambda(spont, sf6-* -> sf6)'
  write(*,*) '  rk7 lambda(HCl, sf6- -> sf5 + ...)'
  write(*,*) '  rk8 lambda(O3, sf6- -> sf6)'
  write(*,*) ' '
  write(*,*) ' Switch SF6_tracer: ',SF6_tracer
  write(*,*) ' Settings for linear tracer: '
  write(*,*) '   Reference year: ',nml_t0
  write(*,*) '   Mixing ratio at t0: ',nml_a_t0
  write(*,*) '   Trend per year : ',nml_b
  IF (SF6_tracer.eq.1) THEN
    write(*,*) ' Additional external tracer set!'
  END IF
  write(*,*) ' Used emode is ',emode
  write(*,*) ' Used autodetachment coefficient: ',autocoeff
  
  CALL read_nml_close(substr, iou, modstr)
  status = 0   ! no error

  CALL end_message(TRIM(modstr),'INITIALISATION', substr)

END SUBROUTINE sf6_read_nml_ctrl

! op_pj_20180418+ transformed to Fortran95 (see below)
!!$SUBROUTINE rnelec(efield,altitude,philat,cossza,col_o2, cair)
!!$!=======================================================================
!!$!     --- electron density for photoelectrons according
!!$!         ansatz by Brasseur
!!$!         
!!$!         electrons are produced by ionization of
!!$!         O2(1Sg-)
!!$!         N2 
!!$!         O
!!$!         O2(1Dg+) (excited state of O2 by ozone decay or rad excitation)
!!$!         NO
!!$!         plus some cosmics
!!$!           
!!$!=======================================================================
!!$!      USE messy_main_data_bi,     ONLY: press_3d, tm1_3d
!!$USE messy_main_timer,         ONLY: DAYOFYEAR
!!$
!!$!     --- parameters
!!$      REAL(dp) :: efield
!!$      REAL(dp) :: altitude
!!$      REAL(dp) :: col_o2
!!$      REAL(dp), INTENT(in) :: philat
!!$      REAL(dp), INTENT(in) :: cossza
!!$      REAL(dp), INTENT(in) :: cair
!!$
!!$!     local declarations
!!$
!!$      integer :: i, j
!!$
!!$!     --- parameterization according Brasseur, p324
!!$      real*8 aion(0:3,3)
!!$      real*8 r1, r2 !, r3
!!$      real*8 rion
!!$      REAL(dp) :: cno
!!$      REAL(dp) :: recombi
!!$      real*8 :: c
!!$      REAL(dp) :: sinfi
!!$      REAL(DP),PARAMETER               :: scale_height_ma = 7000._dp ! Middle atmosphere scale height [m]
!!$
!!$
!!$      data aion / 2.2e6, 2.2e6, 6.4e6, 2.0 , 1.1e6, 1.0e6, 1.3e6, 2.5 , 2.3e6, 3.2e6, 1.0e6, 2.4 /
!!$
!!$!     --- production by N2, O2 , O (1,2,3)
!!$!         modified to have no contribution for about 80 km
!!$!         corresponding an absorption cross section > 5e-20 cm2
!!$      r1(i,c) = aion(0,i) + aion(1,i) * c + aion(2,i) * c**aion(3,i) + exp( c / 50. )
!!$!     --- production by NO
!!$      r2(c)   = 6.e-7 * exp( - 0.4e-20 * c )
!!$!     --- production by O2(1Dg+)
!!$!      r3(c)   = 0.549e-9 * exp( - 2.406e-20 * c ) + 2.6e-9   * exp( - 8.508e-20 * c )
!!$!     -------------------------------------------------------
!!$
!!$
!!$!     --- the NO climatology
!!$!         is described by a simple analytical formula
!!$!         which combines NO fields of Rasch and Ruhnke and the NO-fields
!!$!         given in the COSPAR reference atmosphere III by Siskind
!!$!         the vertical NO mixing ratio is scaled by sunlit time
!!$!         note, that NO is only abundant during day time and
!!$!         that the profile is meant for daytime only
!!$            cno = 1.6/(1.+(philat/20.)**4. )      &
!!$                      * ( exp( - ((altitude/1.e3 - 43. )/10.)**2. )      &
!!$                        + 1.0/( 1. + ((altitude/1.e3 - 60. ) / 20.)**4. )      &
!!$                        )      &
!!$                    + 6. * ( exp( - ((altitude/1.e3 - 43. ) / 10. )**2. )      &
!!$                           + 1.0      &
!!$                             /( 1. + (( altitude/1.e3 - 60. ) / 20.)**4.)      &
!!$                           + 2.0e3      &
!!$                             /( 1. + (( altitude/1.e3 - 100. ) / 4.)**4.)      &
!!$                           )      &
!!$                      * min(1.,max(0.,      &
!!$                                   0.5 + 0.5 / 66. * philat      &
!!$                                       * sin(2*Pi*(DAYOFYEAR-100.)/365.24)      &
!!$                                   )      &
!!$                           )
!!$            cno = cno * 1.d-9 * 2.
!!$!           --- factor of 2 to yield NO as Rasch
!!$!     --- the effective recombination coefficient
!!$!         can be approximated by a tanh function
!!$!         to yield Brasseur's function
!!$         recombi = 10**       &
!!$                     (       &
!!$                     - 5.95       &
!!$                     - 0.65 * tanh( ( altitude - 87.e3 ) / 5.e3 )       &
!!$                     )
!!$!csv         write(*,*) 'recombi1: ', recombi
!!$!csv recombi in 80 km = 4.22161449 ?? 10-6
!!$!csv         in 50 km = 5.01186673 ?? 10-6
!!$
!!$!     --- the lambda coefficient
!!$!         can represented by the following function
!!$!         and is included in the effective recombination
!!$!         the function is taken from Asgeir Brekke, Physics 
!!$!         of the Upper Polar Atmosphere 1997 and comes from
!!$!         Rishbeth, H. and Garriot, O.K., Introduction into 
!!$!         Ionospheric Physics, Academic Press, London, 1969
!!$!         was to be to slow at low altitudes
!!$!         change back to Brasseur's version
!!$         recombi = recombi      &
!!$                   * ( 1.      &
!!$                   + exp(      &
!!$                         -(altitude/1.e3 - 72.)      &
!!$                          /      &
!!$                          ( 2.25 - 0.4      &
!!$                                 * tanh( ( altitude/1.e3 - 70.)      &
!!$                                         / 14.      &
!!$                                       )      &
!!$                          )      &
!!$                         )      &
!!$                      )
!!$!csv         write(*,*) 'recombi2: ', recombi
!!$
!!$!-------------------------------------------------------------
!!$!
!!$      efield = 0.d0
!!$      sinfi=sin(philat/180.*Pi)
!!$!csv      write(*,*) 'philat, sinfi : ', philat, sinfi
!!$            if ( cossza .gt. 1.d-1 ) then
!!$                  rion = 0.d0
!!$!                 --- sum up ionization rates
!!$!                 --- N2, O2, O
!!$!                     the common column
!!$!                     the limit beacuse of overflow at VPP
!!$! above 80 km not included in EMAC
!!$!csv                  if ( altitude(k) .gt. 80.e3 )          &
!!$!csv                 rion = fn2       / r1(1,colstar(k)/cossza)          &
!!$!csv                      + stano2(k) / r1(2,colstar(k)/cossza)          &
!!$!csv                      + stano(k)  / r1(3,colstar(k)/cossza)
!!$!                 --- make agrrement with observations better
!!$                  rion = rion * 20.
!!$!                 --- the ionization of NO
!!$!csv                  rion = rion + r2( col_o2(k) ) * cno
!!$                  rion = rion + r2( col_o2 ) * cno
!!$!                 --- the ionization of O(1Dg+) is negleccted
!!$!                     as the abundance of O(1Dg+) is not known 
!!$!                  rion = rion + r3( colu2(k) / cossza )
!!$!                 --- scale to give particle density in cm-3
!!$!csv                  rion = rion * density / M_air * N_a / 1.e6
!!$                  rion = rion * cair / 1.e6
!!$!                 --- add a constant ionization by cosmic rays
!!$!                     see Brasseur p328 ans p331 
!!$                  rion = rion + ( 0.1 + sinfi**2 )          &
!!$                             * exp( - ( altitude - 40.e3 ) / scale_height_ma )
!!$!                 --- the electron density according Brasseur 6.24
!!$                  efield = sqrt( rion / recombi )
!!$            else
!!$!              --- during night we assume only the cosmics acting
!!$!                  but the recombination coefficient to be same
!!$!                  as during day
!!$                  rion = ( 0.1 + sinfi**2 )          &
!!$                             * exp( - ( altitude - 40.e3 ) / scale_height_ma )
!!$!                 --- the electron density according Brasseur 6.24
!!$                  efield = sqrt( rion / recombi )
!!$            end if
!!$!csv       write(*,*) 'rion : ', rion
!!$!csv       write(*,*) 'efield : ',efield
!!$
!!$END SUBROUTINE rnelec

SUBROUTINE rnelec(efield,altitude,philat,cossza,col_o2, cair, DAYOFYEAR)
!=======================================================================
!     --- electron density for photoelectrons according
!         ansatz by Brasseur
!         
!         electrons are produced by ionization of
!         O2(1Sg-)
!         N2 
!         O
!         O2(1Dg+) (excited state of O2 by ozone decay or rad excitation)
!         NO
!         plus some cosmics
!           
!=======================================================================
!      USE messy_main_data_bi,     ONLY: press_3d, tm1_3d

!     --- parameters
      REAL(dp) :: efield
      REAL(dp) :: altitude
      REAL(dp) :: col_o2
      REAL(dp), INTENT(in) :: philat
      REAL(dp), INTENT(in) :: cossza
      REAL(dp), INTENT(in) :: cair
      INTEGER,  INTENT(in) :: DAYOFYEAR ! op_pj_20180418

!     local declarations

      integer :: i, j

!     --- parameterization according Brasseur, p324
      REAL(dp) :: rion
      REAL(dp) :: cno
      REAL(dp) :: recombi
      REAL(dp) :: sinfi
      REAL(DP), PARAMETER  :: scale_height_ma = 7000._dp ! Middle atmosphere scale height [m]

!     -------------------------------------------------------

!     --- the NO climatology
!         is described by a simple analytical formula
!         which combines NO fields of Rasch and Ruhnke and the NO-fields
!         given in the COSPAR reference atmosphere III by Siskind
!         the vertical NO mixing ratio is scaled by sunlit time
!         note, that NO is only abundant during day time and
!         that the profile is meant for daytime only
            cno = 1.6_dp/(1.+(philat/20._dp)**4 )      &
                      * ( exp( - ((altitude/1.e3_dp - 43._dp )/10._dp)**2 )  &
                        + 1.0_dp/( 1._dp + ((altitude/1.e3_dp - 60._dp ) / 20._dp)**4 )      &
                        )      &
                    + 6._dp * ( exp( - ((altitude/1.e3_dp - 43._dp ) / 10._dp )**2 )      &
                           + 1.0_dp      &
                             /( 1._dp + (( altitude/1.e3_dp - 60._dp ) / 20._dp)**4)      &
                           + 2.0e3_dp      &
                             /( 1._dp + (( altitude/1.e3_dp - 100._dp ) / 4._dp)**4)      &
                           )      &
                      * min(1._dp,max(0._dp,      &
                                   0.5_dp + 0.5_dp / 66._dp * philat      &
                                       * sin(2.0_dp*Pi*(REAL(DAYOFYEAR,dp)-100._dp)/365.24_dp)      &
                                   )      &
                           )
            cno = cno * 1.e-9_dp * 2._dp
!           --- factor of 2 to yield NO as Rasch
!     --- the effective recombination coefficient
!         can be approximated by a tanh function
!         to yield Brasseur's function
         recombi = 10**       &
                     (       &
                     - 5.95_dp       &
                     - 0.65_dp * tanh( ( altitude - 87.e3_dp ) / 5.e3_dp )       &
                     )
!csv         write(*,*) 'recombi1: ', recombi
!csv recombi in 80 km = 4.22161449 ?? 10-6
!csv         in 50 km = 5.01186673 ?? 10-6

!     --- the lambda coefficient
!         can represented by the following function
!         and is included in the effective recombination
!         the function is taken from Asgeir Brekke, Physics 
!         of the Upper Polar Atmosphere 1997 and comes from
!         Rishbeth, H. and Garriot, O.K., Introduction into 
!         Ionospheric Physics, Academic Press, London, 1969
!         was to be to slow at low altitudes
!         change back to Brasseur's version
         recombi = recombi      &
                   * ( 1._dp      &
                   + exp(      &
                         -(altitude/1.e3_dp - 72._dp)      &
                          /      &
                          ( 2.25_dp - 0.4_dp      &
                                 * tanh( ( altitude/1.e3_dp - 70._dp)      &
                                         / 14._dp      &
                                       )      &
                          )      &
                         )      &
                      )
!csv         write(*,*) 'recombi2: ', recombi

!-------------------------------------------------------------
!
      efield = 0._dp
      sinfi=sin(philat/180._dp*Pi)
!csv      write(*,*) 'philat, sinfi : ', philat, sinfi
            if ( cossza .gt. 1.e-1_dp ) then
                  rion = 0._dp
!                 --- sum up ionization rates
!                 --- N2, O2, O
!                     the common column
!                     the limit beacuse of overflow at VPP
! above 80 km not included in EMAC
!csv                  if ( altitude(k) .gt. 80.e3 )          &
!csv                 rion = fn2       / r1(1,colstar(k)/cossza)          &
!csv                      + stano2(k) / r1(2,colstar(k)/cossza)          &
!csv                      + stano(k)  / r1(3,colstar(k)/cossza)
!                 --- make agrrement with observations better
                  rion = rion * 20._dp
!                 --- the ionization of NO
!csv                  rion = rion + r2( col_o2(k) ) * cno
                  rion = rion + r2( col_o2 ) * cno
!                 --- the ionization of O(1Dg+) is negleccted
!                     as the abundance of O(1Dg+) is not known 
!                  rion = rion + r3( colu2(k) / cossza )
!                 --- scale to give particle density in cm-3
!csv                  rion = rion * density / M_air * N_a / 1.e6
                  rion = rion * cair / 1.e6_dp
!                 --- add a constant ionization by cosmic rays
!                     see Brasseur p328 ans p331 
                  rion = rion + ( 0.1_dp + sinfi**2 )          &
                             * exp( - ( altitude - 40.e3_dp ) / scale_height_ma )
!                 --- the electron density according Brasseur 6.24
                  efield = sqrt( rion / recombi )
            else
!              --- during night we assume only the cosmics acting
!                  but the recombination coefficient to be same
!                  as during day
                  rion = ( 0.1_dp + sinfi**2 )          &
                             * exp( - ( altitude - 40.e3_dp ) / scale_height_ma )
!                 --- the electron density according Brasseur 6.24
                  efield = sqrt( rion / recombi )
            end if
!csv       write(*,*) 'rion : ', rion
!csv       write(*,*) 'efield : ',efield

CONTAINS

 ELEMENTAL FUNCTION r1(i, c)
   !     --- production by N2, O2 , O (1,2,3)
   !         modified to have no contribution for about 80 km
   !         corresponding an absorption cross section > 5e-20 cm2
   IMPLICIT NONE
   INTRINSIC :: exp, RESHAPE
   !I/O
   INTEGER,  INTENT(IN) :: i
   REAL(DP), INTENT(IN) :: c
   REAL(dp) :: r1
   ! LOCAL
   REAL(DP), DIMENSION(0:3,3), PARAMETER :: aion = RESHAPE(                      &
        (/ 2.2e6_dp, 2.2e6_dp, 6.4e6_dp, 2.0_dp , 1.1e6_dp, 1.0e6_dp, &
        1.3e6_dp, 2.5_dp , 2.3e6_dp, 3.2e6_dp, 1.0e6_dp, 2.4_dp /),   &
        SHAPE = (/4, 3/) )
   !
   r1 = aion(0,i) + aion(1,i) * c + aion(2,i) * c**aion(3,i) + exp( c / 50._dp )
 END FUNCTION r1

 ELEMENTAL FUNCTION r2(c)
   !     --- production by NO
   IMPLICIT NONE
   INTRINSIC :: exp
   !I/O
   REAL(DP), INTENT(IN) :: c
   REAL(dp) :: r2
   !
   r2 = 6.e-7_dp * exp( - 0.4e-20_dp * c )
 END FUNCTION r2

 ELEMENTAL FUNCTION r3(c)
   !     --- production by O2(1Dg+)
   IMPLICIT NONE
   INTRINSIC :: exp
   !I/O
   REAL(DP), INTENT(IN) :: c
   REAL(dp) :: r3
   !
   r3 = 0.549e-9_dp * exp( - 2.406e-20_dp * c ) + 2.6e-9_dp * &
        exp( - 8.508e-20_dp * c )
 END FUNCTION r3

END SUBROUTINE rnelec
  ! =========================================================================
! op_pj_20180418- transformed to Fortran95

  ! =========================================================================
  SUBROUTINE calc_age(age,vmr_sf6,a,b,t,t0)

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(IN) :: vmr_sf6  ! vmr of SF6 
    REAL(DP), INTENT(IN)  :: a ! reference vmr 
    REAL(DP), INTENT(IN)  :: b ! trend of SF6
    REAL(DP), INTENT(IN) :: t ! actual date (decimal year)
    REAL(DP), INTENT(IN) :: t0 ! reference data (decimal year)
    REAL(DP), INTENT(OUT) :: age  ! mean age of air

    age = t-t0-(vmr_sf6-a)/b
    
  END SUBROUTINE calc_age
  ! =========================================================================

  SUBROUTINE sf6life(lamsf6, emode_read, cossza, time_step_len &
       , altitude, ccvO3,philat, colo2, hydrogen, hcl, cair, rk &
       , efield, DAYOFYEAR)
! originally written for KASIMA by Thomas Reddmann, KIT

!   USE messy_main_data_bi,       ONLY: nlev

   IMPLICIT NONE
   REAL(DP):: lamsf6
   REAL(DP), INTENT(IN) :: colo2
   REAL(DP), INTENT(IN) :: hydrogen, hcl, cair
   INTEGER, INTENT(IN) :: emode_read
   REAL(dp) :: efield
   INTEGER :: melec
   REAL(DP)::    wmoflux(2,158)
   REAL(DP), PARAMETER :: epsilon = .999_dp
   REAL(DP), PARAMETER :: rattach = 2.2e-7_dp   !     --- the rate coefficient of electron attachment according fehsenfeld
   REAL(DP) :: R8
   INTEGER :: i
   REAL(DP), INTENT(OUT) :: rk(8)
   REAL(DP) :: altitude
   REAL(DP) :: rk30
   REAL(DP) :: cossza, time_step_len
   REAL(dp), INTENT(in) :: ccvO3
   REAL(dp) :: philat
   INTEGER,  INTENT(in) :: DAYOFYEAR ! op_pj_20180418

!     ================================================================
!
!     --- the data

!     --- at lambda sol extraterr. F L U SS
!         WMO Report 16 (1985)
!         given in nm, Photons/cm^2 s nm

      data wmoflux /176.20, 1.121E+11,   177.7 , 1.329E+11,   179.30, 1.479E+11, 181.00, 1.856E+11,   182.65, 1.912E+11, &
        184.34, 1.725E+11, 186.05, 2.091E+11,   187.80, 2.683E+11,   189.58, 3.122E+11, 191.39, 3.619E+11,   193.24, 3.696E+11,&
   195.13, 5.024E+11, 197.05, 5.922E+11,   199.01, 6.414E+11,   201.01, 7.525E+11, 203.05, 8.632E+11,   205.13, 1.046E+12, &
  207.26, 1.253E+12, 209.43, 2.070E+12,   211.65, 3.187E+12,   213.91, 3.649E+12, 216.22, 3.590E+12,   218.59, 4.521E+12, &
  221.00, 4.832E+12, 223.47, 6.408E+12,   226.00, 5.247E+12,   228.58, 5.398E+12, 231.22, 5.874E+12,   233.93, 5.044E+12, &
  236.69, 5.712E+12, 239.53, 5.054E+12,   242.43, 7.488E+12,   245.41, 6.607E+12, 248.46, 6.384E+12,   251.58, 6.130E+12, &
  254.79, 8.968E+12, 258.08, 1.486E+13,   261.45, 1.325E+13,   264.91, 3.049E+13, 268.47, 3.331E+13,   272.12, 2.971E+13, &
  275.88, 2.733E+13, 279.73, 2.106E+13,   283.70, 3.777E+13,   287.78, 5.192E+13, 291.99, 8.163E+13,   296.31, 7.747E+13, &
  300.77, 7.119E+13, 305.36, 9.073E+13,   310.10, 1.030E+14,   315.00, 1.088E+14, 320.00, 1.186E+14,   325.00, 1.390E+14, &
  330.00, 1.630E+14, 335.00, 1.562E+14,   340.00, 1.670E+14,   345.00, 1.628E+14, 350.00, 1.706E+14,   355.00, 1.834E+14, &
  360.00, 1.676E+14, 365.00, 2.080E+14,   370.00, 2.200E+14,   375.00, 1.958E+14, 380.00, 2.260E+14,   385.00, 1.778E+14, &
  390.00, 2.280E+14, 395.00, 1.834E+14,   400.00, 3.380E+14,   405.00, 3.400E+14, 410.00, 3.680E+14,   415.00, 3.740E+14, &
  420.00, 3.900E+14, 425.00, 3.620E+14,   430.00, 3.340E+14,   435.00, 3.960E+14, 440.00, 4.040E+14,   445.00, 4.360E+14, &
  450.00, 4.720E+14, 455.00, 4.620E+14,   460.00, 4.780E+14,   465.00, 4.760E+14, 470.00, 4.780E+14,   475.00, 4.880E+14, &
  480.00, 5.020E+14, 485.00, 4.600E+14,   490.00, 4.780E+14,   495.00, 4.960E+14, 500.00, 4.800E+14,   505.00, 4.920E+14, &
  510.00, 4.980E+14, 515.00, 4.640E+14,   520.00, 4.780E+14,   525.00, 4.840E+14, 530.00, 5.100E+14,   535.00, 5.020E+14, &
  540.00, 4.980E+14, 545.00, 5.100E+14,   550.00, 5.060E+14,   555.00, 5.080E+14, 560.00, 5.000E+14,   565.00, 5.140E+14, &
  570.00, 5.160E+14, 575.00, 5.340E+14,   580.00, 5.340E+14,   585.00, 5.400E+14, 590.00, 5.240E+14,   595.00, 5.380E+14, &
  600.00, 5.260E+14, 605.00, 5.360E+14,   610.00, 5.320E+14,   615.00, 5.180E+14, 620.00, 5.380E+14,   625.00, 5.220E+14, &
  630.00, 5.240E+14, 635.00, 5.240E+14,   640.00, 5.260E+14,   645.00, 5.200E+14, 650.00, 5.100E+14,   655.00, 4.960E+14, &
  660.00, 5.140E+14, 665.00, 5.220E+14,   670.00, 5.220E+14,   675.00, 5.240E+14, 680.00, 5.240E+14,   685.00, 5.140E+14, &
  690.00, 5.040E+14, 695.00, 5.200E+14,   700.00, 5.160E+14,   705.00, 5.040E+14, 710.00, 5.020E+14,   715.00, 4.960E+14, &
  720.00, 4.900E+14, 725.00, 4.960E+14,   730.00, 4.900E+14,   735.00, 4.880E+14, 740.00, 4.780E+14,   745.00, 4.800E+14, &
  750.00, 4.820E+14, 755.00, 4.800E+14,   760.00, 4.760E+14,   765.00, 4.680E+14, 770.00, 4.640E+14,   775.00, 4.600E+14, &
  780.00, 4.660E+14, 785.00, 4.680E+14,   790.00, 4.580E+14,   795.00, 4.580E+14, 800.00, 4.540E+14,   805.00, 4.540E+14, &
  810.00, 4.400E+14, 815.00, 4.440E+14,   820.00, 4.360E+14,   825.00, 4.400E+14, 830.00, 4.280E+14,   835.00, 4.280E+14, &
  840.00, 4.260E+14, 845.00, 4.180E+14,   850.00, 4.100E+14/


!     --- the hydrogen profile
!         taken from brasseur, 2ed, p. 443
!         in mixing ratio and the corresponding heights

!csv_brass      data hydro / 3.2e-16, 1.6e-14, 4.2e-13, 4.1e-12, 2.2e-11, 9.7e-11, 4.5e-10, 2.5e-9, 2.0e-8, 2.0e-7,  4.6e-6,  7.2e-6,  9.0e-6,  1.0e-5/
!csv_brass      data zhydro / 35.e3, 40.e3, 45.e3, 50.e3, 55.e3, 60.e3, 65.e3, 70.e3, 75.e3, 80.e3, 85.e3, 90.e3, 95.e3, 100.e3/

   emode = emode_read

   if ( emode .ge. 10 ) then
    r8      = 1.2e-9 ! huey
    melec   = emode - 10
   else
    melec   = emode
    r8      = 0.032e-9  ! fehsenfeld
   endif


!     --- interpolate the hydrogen profile
!csv_brass      do k=1,14
!csv_brass         hydro(k) = log( hydro(k) )
!csv_brass      end do

!csv_brass      call splinr(zhydro,hydro,14,0.d0,0.d0,sbuf)

!csv      do k=1,nlev
!csv_brass         call splint(zhydro,hydro,sbuf,14,altitude(jk),profile(jk))
!csv         write(*,*) 'jk, Alt, prof: ',jk, altitude(jk), profile(jk)
!csv_brass         rk(4) = exp(profile(jk))
         rk(4) = hydrogen * 1.e-06
!csv         write(*,*) 'rk4a: ',rk(4)
!        --- the number density in cm-3
!csv         rk(4) = rk(4) * density(jk) / M_air / 1.e6 * N_a
!csv         write(*,*) 'rk4b: ',rk(4)
!        --- the reaction rate coefficient
!            increase H concentration by 1.5 
!            because of higher Clorine burden 
         rk(4) = 2.1e-10 * rk(4) * 1.5
!csv         write(*,*) 'rk4c: ',rk(4)
!csv      end do

!c     --- the electron profile
!c         mode 2,3:
!c         taken from brasseur, 2ed, p 322 is described
!c         by simple exponential profile between
!c         50 and 70 km, overestimating the electron 
!c         density presumably
!c         in cm-3
!c         mode 1:
!c         given as 3D-field according approximation of Brasseur
      if ( melec .eq. 2 .or. melec.eq.0 ) then
            rk(2) = rattach * exp( log(10.) * ( altitude/1.e3 - 50. )/ 10. )
      else if ( melec .eq. 3 ) then
            rk(2) = rattach * exp( log(10.) * ( altitude/1.e3 - 40. )/ 10. )
      else if ( melec .eq. 1 ) then
         call rnelec(efield,altitude,philat,cossza,colo2,cair,DAYOFYEAR)
      end if

!c     --- the photodetachment rate of sf6-
!c         datskos and ingolfsson
!c         integration over the solar and absorption spectrum
!c         we assume no absorption as the important part of the
!c         spectrum is here in the near-UV and we are above 
!c         the stratopause
      rk30 = 0._dp
      do i=2,156
         if ( wmoflux(1,i) .le. 337. ) then
            rk30 = rk30 + ( wmoflux(1,i+1) - wmoflux(1,i-1) ) / 2 * wmoflux(2,i) * 2.e-18
         else if ( wmoflux(1,i) .gt. 337. .and. wmoflux(1,i) .lt. 386. ) then
            rk30 = rk30 + ( wmoflux(1,i+1) - wmoflux(1,i-1) ) / 2 * wmoflux(2,i) * 1.22e3 / &
                   wmoflux(1,i) * ( 1.483e-12 * ( 1.22e3 / wmoflux(1,i) - 3.16 ) )** 1.5_dp
         end if
!c         print '(I3,4(G14.5))', i, wmoflux(1,i), wmoflux(2,i), rk30 
      end do

!c     --- vertical profile of photodetachment and autodetachment
!c         is constant
         rk(3) = rk30
         rk(6) = autocoeff

!c     --- stabilization of SF6-* by collisions
!c         collision cross section * mean molecule velocity * air density
!c         1.9 e-10 cm-3 (Odom etal.)                       * air density
!csv         rk(5) = 1.9e-10 * density(jk) / M_air / 1.e6 * N_A
        rk(5) = 1.9e-10 * cair / 1.e6

!c     --- the HCl pofile and reaction rate from huey
!csv         rk(7) = 1.5e-9 * min(3.5_dp, 3.5_dp*altitude(jk)/5.e4) * 1.e-9 * density(jk) / M_air / 1.e6 * N_a
        rk(7) = 1.5e-9 * hcl / 1.e6   

      if ( melec .eq. 1 ) rk(2) = rattach * efield
!csv      rk(8) = ccvO3 * r8 * density(jk) / M_air / 1.e6 * N_a
        rk(8) = ccvO3 * r8 / 1.e6



!c     --- we assume the hydrogen mixing ratio proportional to cos(theta)
!c         according the chemistry modul of R. Ruhnke
!c         not for the electron density
!c
!c         "Wenn ich an SF6 denke in der Nacht, ...."
!c         keine Elektronen f??r exp - Profile

      if ( melec .eq. 2 .or. melec .eq. 3) then

         if (cossza .gt. 0._dp ) then
!csv             rk(8) = ccvO3 * r8 * density(jk) / M_air / 1.e6 * N_a
             rk(8) = ccvO3 * r8 / 1.e6
             lamsf6 = rk(2)      &
                * (      &
                1._dp - epsilon      &
                * ( rk(5) * ( rk(3) + rk(8) )      &
                    + rk(6) * ( rk(4) * cossza + rk(3)      &
                                + rk(7) + rk(8)  ) )      &
                / ( rk(6) + rk(5) )      &
                / ( rk(3) + rk(4) * cossza      &
                                          + rk(7)  + rk(8))      &
                )      &
                * time_step_len
            lamsf6 = max(0._dp, lamsf6 )
            lamsf6 = min(1._dp, lamsf6 )
         else
            lamsf6 = 0._dp
         end if
      else if ( melec .eq. 0) then
         if (cossza .gt. 0._dp ) then
            lamsf6 = rk(2) * time_step_len
            lamsf6 = max(0._dp, lamsf6 )
            lamsf6 = min(1._dp, lamsf6 )
         else
            lamsf6 = 0._dp
         end if

      else if ( melec .eq. 4 ) then

            if (cossza .gt. 0._dp ) then
!csv                  rk(8) = ccvO3 * r8 * density(jk) / M_air / 1.e6 * N_a
                  rk(8) = ccvO3 * r8 / 1.e6
                  lamsf6 = rk(2)       &
                * (       &
                1._dp - epsilon       &
                * ( rk(5) * rk(3)       &
                    + rk(6) * ( rk(4) * cossza + rk(3) ))       &
                / ( rk(6) + rk(5) )       &
                / ( rk(3) + rk(4) * cossza )       &
                )       &
                * time_step_len
            lamsf6 = max(0._dp, lamsf6 )
            lamsf6 = min(1._dp, lamsf6 )
         else
            lamsf6 = 0._dp
         end if

      else if ( melec .eq.1 ) then
            if (cossza .gt. 0._dp ) then
                  rk(2) = rattach * efield
!csv                  rk(8) = ccvO3 * r8 * density(jk) / M_air / 1.e6 * N_a
                  rk(8) = ccvO3 * r8 / 1.e6
                  lamsf6 = rk(2)        &
                * (        &
                1._dp - epsilon        &
                * ( rk(5) * ( rk(3) + rk(8) )        &
                    + rk(6) * ( rk(4) * cossza + rk(3)        &
                                + rk(7) + rk(8)  ) )        &
                / ( rk(6) + rk(5) )        &
                / ( rk(3) + rk(4) * cossza        &
                                          + rk(7)  + rk(8))        &
                )        &
                * time_step_len
            else
                  rk(2) = rattach * efield
!csv                  rk(8) = ccvO3 * r8 * density(jk) / M_air / 1.e6 * N_a
                  rk(8) = ccvO3 * r8 / 1.e6
                  lamsf6 = rk(2)          &
                * (          &
                1._dp - epsilon          &
                * ( rk(5) * ( rk(3) + rk(8) )          &
                    + rk(6) * ( rk(3)          &
                                + rk(7) + rk(8)  ) )          &
                / ( rk(6) + rk(5) )          &
                / ( rk(3) + rk(7)  + rk(8))          &
                )          &
                * time_step_len
            end if

           lamsf6 = max(0._dp, lamsf6 )
           lamsf6 = min(1._dp, lamsf6 )

      else if ( melec .eq. 0) then
            if (cossza .gt. 0._dp ) then
                lamsf6 = rk(2) * time_step_len
                lamsf6 = max(0._dp, lamsf6 )
                lamsf6 = min(1._dp, lamsf6 )
            else
                  lamsf6 = 0._dp
            end if


      endif
!csv      write(*,*) 'RK TESTOUTPUT INTERN 1: ', rk(1)
!csv      write(*,*) 'RK TESTOUTPUT INTERN 2: ', rk(2)
!csv      write(*,*) 'RK TESTOUTPUT INTERN 3: ', rk(3)
!csv      write(*,*) 'RK TESTOUTPUT INTERN 4: ', rk(4)
  

  END SUBROUTINE sf6life


END MODULE messy_SF6
