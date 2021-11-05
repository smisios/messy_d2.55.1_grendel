!        1         2         3         4         5         6         7         8
!2345678901234567890123456789012345678901234567890123456789012345678901234567890

MODULE messy_made_box

   ! DESCRIPTION
   ! Interface layer for MADE box model, analog messy_made_e5.f90
   !
   ! AUTHORS
   ! Axel Lauer, DLR Oberfpaffenhofen, Germany
   !    questions/suggestions: axel.lauer@dlr.de
   ! Christopher Kaiser, DLR Oberpfaffenhofen, Germany
   !
   ! LAST CHANGES
   ! 2005
   ! 2012 by CK
   !      - renamed subroutine made_initialize to made_init
   !      - added subroutine made_box_init
   !      - modified structure to enable setting of various parameters via
   !        namelist instead of setting them inside the code

   USE messy_made,               ONLY: modstr,modver, &
! op_ck_20120430+
                                       nspcsda
! op_ck_20120430-
   USE messy_main_constants_mem, ONLY: dp

   IMPLICIT NONE
   PRIVATE

! op_ck_20120430+

   ! *** Constants ***

   INTEGER, PARAMETER, PUBLIC :: blksize = 1   ! array size
   INTEGER, PARAMETER, PUBLIC :: numcells = 1  ! number of cells (boxes)

   INTEGER, PARAMETER, PUBLIC :: nummod  = 3   ! number of modes
   INTEGER, PARAMETER, PUBLIC :: numspec = 8   ! number of (chemical) species

   INTEGER, PARAMETER, PUBLIC :: so4   = 1     ! index for SO4
   INTEGER, PARAMETER, PUBLIC :: nh4   = 2     ! etc...
   INTEGER, PARAMETER, PUBLIC :: no3   = 3
   INTEGER, PARAMETER, PUBLIC :: ss    = 4
   INTEGER, PARAMETER, PUBLIC :: pom   = 5
   INTEGER, PARAMETER, PUBLIC :: bc    = 6
   INTEGER, PARAMETER, PUBLIC :: du    = 7
   INTEGER, PARAMETER, PUBLIC :: h2o   = 8

   INTEGER, PARAMETER, PUBLIC :: ncon  = 9     ! index for number concentration
   INTEGER, PARAMETER, PUBLIC :: tmcon = 10    ! index for 3rd moment conc.

   INTEGER, PARAMETER, PUBLIC :: h2so4 = 11    ! index for H2SO4(g)
   INTEGER, PARAMETER, PUBLIC :: nh3   = 12    ! index for NH3(g)
   INTEGER, PARAMETER, PUBLIC :: hno3  = 13    ! index for HNO3(g)

   REAL(dp), PARAMETER, PUBLIC :: conmin = 1.e-30_dp  ! conc. lower lim. [ug/m3]
   REAL(dp), PARAMETER, PUBLIC :: dgmin  = 1.e-9_dp   ! lowest median diam. [m]
   REAL(dp), PARAMETER, PUBLIC :: one3   = 1._dp / 3._dp
   REAL(dp), PUBLIC            :: es36(nummod)        ! for diameter calculation

   ! *** Namelist-set parameters ***

   INTEGER,  PUBLIC :: timesteps              ! number of timesteps
   REAL(dp), PUBLIC :: tmst                   ! time step duration [s]

   REAL(dp), PUBLIC :: pressure(blksize)      ! pressure [Pa]
   REAL(dp), PUBLIC :: temperature(blksize)   ! temperature [K]
   REAL(dp), PUBLIC :: relhum(blksize)        ! relative humidity [0-1]
   REAL(dp), PUBLIC :: cloudcover(blksize)    ! fract. cloud cover [0-1]
   REAL(dp), PUBLIC :: rh_hist_akn(blksize)   ! rel. hum. history
                                              ! (---> hysteresis)
   REAL(dp), PUBLIC :: rh_hist_acc(blksize)
   REAL(dp), PUBLIC :: rh_hist_cor(blksize)
   REAL(dp), PUBLIC :: rh_hist_akns(blksize)
   REAL(dp), PUBLIC :: rh_hist_accs(blksize)
   REAL(dp), PUBLIC :: rh_hist_sooti(blksize)
   REAL(dp), PUBLIC :: rh_hist_sootj(blksize)
   REAL(dp), PUBLIC :: so4rat(blksize)        ! H2SO4(g), rate of change
                                              ! [ug/(m3s)]
   REAL(dp), PUBLIC :: no3rat(blksize)        ! HNO3(g), rate of change
                                              ! [ug/(m3s)]
   REAL(dp), PUBLIC :: soa(blksize)           ! SOA (gas phase) "emissions"
                                              ! [ug/(m3s)]
   REAL(dp), PUBLIC :: em_bc_mass(3,blksize)  ! BC mass emission rate [ug m-3 s-1]
   REAL(dp), PUBLIC :: em_bc_num(3,blksize)   ! BC number emission rate [# m-3 s-1]

   LOGICAL, PUBLIC  :: ltest_mcon             ! test mass conservation
   LOGICAL, PUBLIC  :: ltest_ncon             ! test number consistency
   LOGICAL, PUBLIC  :: ltest_so4              ! test SO4 production
   LOGICAL, PUBLIC  :: ltest_adapdt           ! test if adap. time step used

   ! *** Global variables ***

   REAL(dp), TARGET, PUBLIC :: tracer(blksize,nspcsda)   ! array that is passed
                                                         ! to messy_madein
   REAL(dp), PUBLIC         :: conv(numspec+2+3)         ! conversion factors
   REAL(dp), PUBLIC         :: dg(nummod)                ! median diameters

   ! *** Pointer array ***

   ! Array of pointers that can be set individually
   TYPE :: parr
      REAL(dp), POINTER :: p
   END TYPE parr
   PUBLIC :: parr ! op_pj_20140214
   TYPE(parr), PUBLIC       :: aero(blksize,nummod,numspec+2) ! pointer array
                                                              ! for loop
                                                              ! simplification
   REAL(dp), TARGET, PUBLIC :: zero           ! target for pointers
                                              ! without meaning

! op_ck_20120430-

   PUBLIC :: made_init, made_box_init, made_finalize

   CONTAINS

   ! --------------------------------------------------------------------------

! op_ck_20120430+

   SUBROUTINE made_box_init   ! ATTENTION: Must be called AFTER made_init.

     USE messy_main_constants_mem, ONLY: dp, pi, R_gas

     USE messy_made, ONLY: akn, acc, cor, &
           vnu0, vac0, vcorn, &
           vso4ai, vso4aj, &
           vnh4ai, vnh4aj, &
           vno3ai, vno3aj, &
           vseasi, vseasj, vseasc, &
           vorgpai, vorgpaj, &
           veci, vecj, &
           vdustj, vdustc, &
           vh2oai, vh2oaj, vh2oac, &
           mwso4, mwnh4, &
           mwso4, mwnh4, mwno3, mwseas, mworg, mwec, mwsoil, mwh2o, &
           mwh2so4, mwnh3, mwhno3, vsulf, vnh3, vhno3, &
           vnu3, vac3, vcor3, &
           rhoso4, rhonh4, rhono3, rhoseas, rhoorg, rhoanth, rhosoil, rhoh2o, &
           f6dpim9, sigma

      USE messy_main_blather, ONLY: start_message, end_message
      USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

      IMPLICIT NONE

      INTRINSIC EXP, LOG, MAX

      CHARACTER(LEN=*), PARAMETER :: substr = 'made_box_init'
      INTEGER,          PARAMETER :: iou = 99     ! logical I/O unit
      LOGICAL                     :: lex          ! file exists ?
      INTEGER                     :: fstat        ! file status
      INTEGER                     :: i, j, k      ! loop indices
      REAL(dp)                    :: rho(numspec) ! aerosol species densities

      NAMELIST /BOXINIT/ TIMESTEPS, TMST,                       &
           PRESSURE, TEMPERATURE, RELHUM, CLOUDCOVER,           &
           RH_HIST_AKN, RH_HIST_ACC, RH_HIST_COR, RH_HIST_AKNS, &
           RH_HIST_ACCS, RH_HIST_SOOTI, RH_HIST_SOOTJ,          &
           SO4RAT, NO3RAT, SOA, EM_BC_MASS, EM_BC_NUM, TRACER,  &
           LTEST_MCON, LTEST_NCON, LTEST_SO4,                   &
           LTEST_ADAPDT

      do i = 1, blksize
         do j = 1, nspcsda
            tracer(i,j) = 0._dp
         end do
      end do

!!$      timesteps     = 20
!!$      tmst          = 1800._dp
!!$      pressure      = 101325._dp
!!$      temperature   = 298._dp
!!$      relhum        = 0.7_dp
!!$      cloudcover    = 0._dp
!!$      rh_hist_akn   = 2._dp
!!$      rh_hist_acc   = 2._dp
!!$      rh_hist_cor   = 2._dp
!!$      rh_hist_akns  = 2._dp
!!$      rh_hist_accs  = 2._dp
!!$      rh_hist_sooti = 2._dp
!!$      rh_hist_sootj = 2._dp
!!$      so4rat        = 1.67e-4_dp
!!$      soa           = 0._dp

!!$      do i = 1, blksize
!!$
!!$         do j = 1, nspcsda
!!$            tracer(i,j) = 0._dp
!!$         end do
!!$
!!$         tracer(i,vnu0)    = 5.e9_dp
!!$         tracer(i,vac0)    = 3.e8_dp
!!$         tracer(i,vcorn)   = 1.e6_dp
!!$
!!$         tracer(i,vso4ai)     = 0.15_dp
!!$         tracer(i,vso4aj)     = 3._dp
!!$
!!$         tracer(i,vnh4ai)     = tracer(i,vso4ai)/mwso4*mwnh4*2.0
!!$         tracer(i,vnh4aj)     = tracer(i,vso4aj)/mwso4*mwnh4*2.0
!!$
!!$         tracer(i,vno3ai)     = 0.025_dp
!!$         tracer(i,vno3aj)     = 0.25_dp
!!$
!!$         tracer(i,vseasi)     = 0.025_dp
!!$         tracer(i,vseasj)     = 0.250_dp
!!$         tracer(i,vseasc)     = 5.000_dp
!!$
!!$         tracer(i,vdustc)     = 4._dp
!!$
!!$         tracer(i,vecis)      = 5._dp
!!$         tracer(i,vecjs)      = 10._dp
!!$         tracer(i,vecsooti)   = 4._dp
!!$         tracer(i,vecsootj)   = 8._dp
!!$
!!$         ! aerosol water is calculated by EQSAM
!!$
!!$      end do

      ! read BOXINIT namelist

      CALL start_message(TRIM(modstr),'read namelist and initialize ' &
           // 'parameters and variables     for MADE box model', substr)

      CALL read_nml_open(lex, substr, iou, 'BOXINIT', modstr)
      IF (.not.lex) STOP 1   ! <modstr>.nml does not exist

      READ(iou, NML=BOXINIT, IOSTAT=fstat)
      CALL read_nml_check(fstat, substr, iou, 'BOXINIT', modstr)
      IF (fstat /= 0) STOP 1  ! error while reading namelist

      CALL read_nml_close(substr, iou, modstr)

      WRITE(*,*) ' '
      WRITE(*,*) ' initialization data successfully read in'
      WRITE(*,*) ' '

      CALL end_message(TRIM(modstr),'read namelist and initialize ' &
           // 'variables     for MADE box model', substr)

      ! set conversion factors
      ! (from mol/mol(air) to #/m3 and ug/m3, respectively)

      conv(ncon)  = pressure(1) / (R_gas * temperature(1))

      conv(so4)   = 1.e6 * mwso4   * conv(ncon)
      conv(nh4)   = 1.e6 * mwnh4   * conv(ncon)
      conv(no3)   = 1.e6 * mwno3   * conv(ncon)
      conv(ss)    = 1.e6 * mwseas  * conv(ncon)
      conv(pom)   = 1.e6 * mworg   * conv(ncon)
      conv(bc)    = 1.e6 * mwec    * conv(ncon)
      conv(du)    = 1.e6 * mwsoil  * conv(ncon)
      conv(h2o)   = 1.e6 * mwh2o   * conv(ncon)

      conv(h2so4) = 1.e6 * mwh2so4 * conv(ncon)
      conv(nh3)   = 1.e6 * mwnh3   * conv(ncon)
      conv(hno3)  = 1.e6 * mwhno3  * conv(ncon)

      ! map densities
      rho(so4)    = rhoso4
      rho(nh4)    = rhonh4
      rho(no3)    = rhono3
      rho(ss)     = rhoseas
      rho(pom)    = rhoorg
      rho(bc)     = rhoanth
      rho(du)     = rhosoil
      rho(h2o)    = rhoh2o

      ! set factors for diameter calculation
      es36(akn)   = EXP(4.5_dp * (LOG(sigma(1)))**2)
      es36(acc)   = EXP(4.5_dp * (LOG(sigma(2)))**2)
      es36(cor)   = EXP(4.5_dp * (LOG(sigma(3)))**2)

      ! target for unused pointers in array aero
      zero = 0._dp

      do i = 1, blksize

         ! associate pointers in array aero

         do j = 1, nummod
            do k = 1, numspec + 1
               aero(i,j,k)%p => zero
            end do
         end do

         aero(i,akn,so4)%p   => tracer(i,vso4ai)
         aero(i,acc,so4)%p   => tracer(i,vso4aj)

         aero(i,akn,nh4)%p   => tracer(i,vnh4ai)
         aero(i,acc,nh4)%p   => tracer(i,vnh4aj)

         aero(i,akn,no3)%p   => tracer(i,vno3ai)
         aero(i,acc,no3)%p   => tracer(i,vno3aj)

         aero(i,akn,ss)%p   => tracer(i,vseasi)
         aero(i,acc,ss)%p   => tracer(i,vseasj)
         aero(i,cor,ss)%p   => tracer(i,vseasc)

         aero(i,akn,pom)%p   => tracer(i,vorgpai)
         aero(i,acc,pom)%p   => tracer(i,vorgpaj)

         aero(i,akn,bc)%p  => tracer(i,veci)
         aero(i,acc,bc)%p  => tracer(i,vecj)

         aero(i,acc,du)%p   => tracer(i,vdustj)
         aero(i,cor,du)%p   => tracer(i,vdustc)

         aero(i,akn,h2o)%p   => tracer(i,vh2oai)
         aero(i,acc,h2o)%p   => tracer(i,vh2oaj)
         aero(i,cor,h2o)%p   => tracer(i,vh2oac)

         aero(i,akn,ncon)%p   => tracer(i,vnu0)
         aero(i,acc,ncon)%p   => tracer(i,vac0)
         aero(i,cor,ncon)%p   => tracer(i,vcorn)

         aero(i,akn,tmcon)%p   => tracer(i,vnu3)
         aero(i,acc,tmcon)%p   => tracer(i,vac3)
         aero(i,cor,tmcon)%p   => tracer(i,vcor3)

         ! convert to #/m3 and ug/m3, respectively,
         ! and calculate 3rd mom. conc. [m3 m-3]
         do j = 1, nummod
            do k = 1, numspec
               aero(i,j,k)%p     = conv(k) * aero(i,j,k)%p
               aero(i,j,tmcon)%p = aero(i,j,tmcon)%p &
                                   + f6dpim9 * aero(i,j,k)%p / rho(k)
            end do
            aero(i,j,ncon)%p  = conv(ncon) * aero(i,j,ncon)%p
            aero(i,j,tmcon)%p = MAX(conmin, aero(i,j,tmcon)%p)
         end do
         tracer(1,vsulf) = conv(h2so4) * tracer(1,vsulf)
         tracer(1,vnh3)  = conv(nh3)   * tracer(1,vnh3)
         tracer(1,vhno3) = conv(hno3)  * tracer(1,vhno3)

         ! Standard deviation fixed in all modes, so diagnose diameter
         ! from 3rd moment and number concentrations:
         dg(akn)   = MAX(dgmin, (tracer(i,vnu3)  / (tracer(i,vnu0)  &
                                 * es36(akn)))**one3)
         dg(acc)   = MAX(dgmin, (tracer(i,vac3)  / (tracer(i,vac0)  &
                                 * es36(acc)))**one3)
         dg(cor)   = MAX(dgmin, (tracer(i,vcor3) / (tracer(i,vcorn) &
                                 * es36(cor)))**one3)

      end do

   END SUBROUTINE made_box_init

   ! --------------------------------------------------------------------------

! op_ck_20120430-

   SUBROUTINE made_init

      USE messy_made,         ONLY: made_read_nml, made_initialize_core
      USE messy_main_blather, ONLY: start_message, end_message

      IMPLICIT NONE

      CHARACTER(LEN=*), PARAMETER :: substr = 'made_init'
      INTEGER :: status ! error status

      CALL start_message(TRIM(modstr),'read namelist and initialize made', &
                         substr)

      ! *** ATTENTION ***
      ! Namelist MUST be read BEFORE calling 'made_initialize_core'!!!

      ! read CTRL namelist
      CALL made_read_nml(status, 99)
      IF (status /= 0) STOP

      WRITE(*,*) ' '
      WRITE(*,*) ' input data successfully read in'
      WRITE(*,*) ' '

      CALL made_initialize_core
      CALL end_message(TRIM(modstr),'read namelist and initialize made', substr)

   END SUBROUTINE made_init

   ! --------------------------------------------------------------------------

   SUBROUTINE made_finalize

      !**** *SUBROUTINE* *MADE_FINALIZE*  Finalize MADE after last call.
      !
      !         A. Lauer   DLR Oberpfaffenhofen  11/2004
      !
      !     PURPOSE.
      !     --------
      !
      !         1. Finalize MADE after last call.
      !         2. Clean up.
      !
      !**   INTERFACE.
      !     ----------
      !
      !       *MADE_FINALIZE* IS CALLED FROM *BOXMADE*
      !
      !       *PARAMETERS*: none
      !
      !
      !     METHODS.
      !     --------
      !
      !         -
      !
      !     EXTERNALS.
      !     ----------
      !
      !         None.
      !
      !     REFERENCE.
      !     ----------
      !
      !         None.
      !
      !     ------------------------------------------------------------

      IMPLICIT NONE

      !*        1. NOTHING TO DO YET
      !            ------- -- -- ---

   END SUBROUTINE made_finalize

   ! --------------------------------------------------------------------------

END MODULE messy_made_box

!*****************************************************************************

PROGRAM made

   USE messy_main_constants_mem, ONLY: dp,                          &
! op_ck_20120320+
                                       R_gas,                       &
! op_ck_20120320-
! op_ck_20120430+
                                       pi
! op_ck_20120430-
   USE messy_made, ONLY: nspcsda, vso4ai, vso4aj, vnh4ai, vnh4aj,   &
                         vno3ai, vno3aj, vnu0, vac0, vcorn, vseasi, &
                         vseasj, vseasc, vdustj, vdustc, vh2oai,    &
                         vh2oaj, vh2oac, vsulf, vnh3, vhno3, veci,  &
                         vecj, vorgpai, vorgpaj, mwso4, mwnh4,      &
                         mwno3, mwnh3, mwhno3, mwh2so4, lmade,      &
! op_ck_20120430+
                         akn, acc, cor, vnu3, vac3, vcor3,          &
! op_ck_20120430-
                         rhoh2o, f6dpim9

! op_ck_20120430+
   USE messy_made_box, ONLY: blksize, numcells, nummod, numspec,   &
        ncon, so4, nh4, no3, ss, pom, bc, du, h2o,                 &
        timesteps, tmst,                                           &
        pressure, temperature, relhum, cloudcover,                 &
        rh_hist_akn, rh_hist_acc, rh_hist_cor,                     &
        rh_hist_akns, rh_hist_accs, rh_hist_sooti, rh_hist_sootj,  &
        so4rat, no3rat, soa, em_bc_mass, em_bc_num,                &
        tracer, aero,                                              &
!!$        ltest_mcon, ltest_ncon, ltest_so4,                         &
!!$        ltest_adapdt,                                              &
        dg, es36, one3, dgmin
! op_ck_20120430-

   ! SUBROUTINES

   USE messy_made_box, ONLY: made_init, made_box_init, made_finalize
   USE messy_made, ONLY: made_main

   IMPLICIT NONE

! op_ck_20120430+

   INTRINSIC :: EXP, LOG, MAX, SUM

   INTEGER   :: i
   INTEGER   :: o_species = 107          ! unit for output of species mass
                                         ! time series
   INTEGER   :: o_dry = 108              ! unit for output of dry diam.s/masses

!!$   integer, parameter :: timesteps = 480 ! number of timesteps
!!$   integer, parameter :: blksize   = 1  ! array dimension
!!$   integer, parameter :: numcells  = 1  ! number of cells in arrays
!!$ 
!!$   real(dp) :: tmst                     ! time step [s]
!!$
!!$   integer :: i, j
!!$
!!$   real(dp) :: cloudcover(blksize)     ! fractional cloud cover [0-1]
!!$   real(dp) :: so4rat(blksize)         ! H2SO4(g), rate of change [ug/m3/s]
!!$   real(dp) :: pressure(blksize)       ! pressure [Pa]
!!$   real(dp) :: rh_hist_akn(blksize)    ! rel. hum. history (---> hysteresis)
!!$   real(dp) :: rh_hist_acc(blksize)    ! rel. hum. history (---> hysteresis)
!!$   real(dp) :: rh_hist_cor(blksize)    ! rel. hum. history (---> hysteresis)
!!$   real(dp) :: soa(blksize)            ! SOA (gas phase) "emissions" [kg/kg/s]
!!$   real(dp) :: relhum(blksize)         ! relative humidity [0-1]
!!$   real(dp) :: temperature(blksize)    ! temperature [K]
!!$
!!$   real(dp) :: tracer(blksize,nspcsda) ! tracer conc.
!!$
!!$   ! modal diameters (diagnostic output for box-version)
!!$
!!$   REAL(dp) :: DGNUC(BLKSIZE)    ! Aitken mode geometric mean diameter (wet) [m]
!!$   REAL(dp) :: DGACC(BLKSIZE)    ! acc. mode geometric mean diameter (wet) [m]
!!$   REAL(dp) :: DGCOR(BLKSIZE)    ! coarse mode geometric mean diameter (wet) [m]
!!$

! op_ck_20120430-

   REAL(dp) :: DGDRYNUC(BLKSIZE) ! Aitken mode geometric mean diameter (dry) [m]
   REAL(dp) :: DGDRYACC(BLKSIZE) ! acc. mode geometric mean diameter (dry) [m]
   REAL(dp) :: DGDRYCOR(BLKSIZE) ! coarse mode geometric mean diameter (dry) [m]
   REAL(dp) :: DENSN(BLKSIZE)    ! average model density, Aitken mode [kg/m3]
   REAL(dp) :: DENSA(BLKSIZE)    ! average model density, acc. mode [kg/m3]
   REAL(dp) :: DENSC(BLKSIZE)    ! average model density, coarse mode [kg/m3]

! op_ck_20120430+

!!$
!!$! op_ck_20120320+
!!$   REAL(dp) :: poRT, h2so4fac, hno3fac, nh3fac
!!$! op_ck_20120320-
!!$
!!$!   real(dp) :: so4sum0, so4sum1
!!$!   real(dp) :: nh4sum0, nh4sum1
!!$!   real(dp) :: no3sum0, no3sum1
!!$!   real(dp) :: bcsum0, bcsum1
!!$!   real(dp) :: pomsum0, pomsum1
!!$!   real(dp) :: sssum0, sssum1
!!$!   real(dp) :: dusum0, dusum1

!!$   integer :: testout = 106       ! unit for optional output of tests
!!$
!!$   REAL(dp) :: oldmass, newmass   ! for mass conservation test
!!$   REAL(dp) :: oldnum, newnum     ! for number consistency test
!!$   REAL(dp) :: oldso4, newso4     ! for SO4 production test
!!$   REAL(dp) :: oldtmst            ! for adaptive time step test

! op_ck_20120430-

   ! ----------------------------------------------------------------------
   !     Set initial values - box version only.
   ! ----------------------------------------------------------------------

! op_ck_20120430+
!!$
!!$! op_ck_20120320+
!!$! op_ck_20120320: old initialization
!!$!!$   tmst = 1800.0                      ! s
!!$!!$
!!$!!$   do i = 1, numcells
!!$!!$
!!$!!$      cloudcover(i)   = 0.0           ! frac
!!$!!$      so4rat(i)       = 1.67e-4       ! ug/m3/s
!!$!!$      pressure(i)     = 101325.0      ! Pa
!!$!!$      rh_hist_akn(i)  = 2.0           ! 
!!$!!$      rh_hist_acc(i)  = 2.0           ! 
!!$!!$      rh_hist_cor(i)  = 2.0           ! 
!!$!!$      soa(i)          = 0.0           ! ug/m3/s
!!$!!$      relhum(i)       = 0.7           ! (0-1)
!!$!!$      temperature(i)  = 280.0         ! K
!!$!!$
!!$!!$      ! i = Aitken mode
!!$!!$      ! j = accumulation mode
!!$!!$      ! c = coarse mode
!!$!!$
!!$!!$      do j = 1, nspcsda
!!$!!$         tracer(i,j) = 0.0
!!$!!$      end do
!!$!!$
!!$!!$      tracer(i,vsulf)   = 0.0         ! H2SO4(g)       ug/m3
!!$!!$      tracer(i,vnh3)    = 1.2         ! NH3(g)         ug/m3
!!$!!$      tracer(i,vhno3)   = 4.5         ! HNO3(g)        ug/m3
!!$!!$
!!$!!$      tracer(i,vnu0)    = 5.000e+9    ! number, i      #/m3
!!$!!$      tracer(i,vac0)    = 3.000e+8    ! number, j      #/m3
!!$!!$      tracer(i,vcorn)   = 1.000e+6    ! number, c      #/m3
!!$!!$
!!$!!$      tracer(i,vso4ai)  = 0.150       ! SO4, i         ug/m3
!!$!!$      tracer(i,vso4aj)  = 3.0         ! SO4, j         ug/m3
!!$!!$
!!$!!$      tracer(i,vnh4ai)  = &           ! NH4, i         ug/m3
!!$!!$                          tracer(i,vso4ai)/mwso4*mwnh4*2.0
!!$!!$      tracer(i,vnh4aj)  = &           ! NH4, j         ug/m3
!!$!!$                          tracer(i,vso4aj)/mwso4*mwnh4*2.0
!!$!!$
!!$!!$!      tracer(i,vno3ai)  = 0.0         ! NO3, i         ug/m3
!!$!!$!      tracer(i,vno3aj)  = 0.0         ! NO3, j         ug/m3
!!$!!$
!!$!!$      tracer(i,vseasi)  = 0.025        ! seasalt, i     ug/m3
!!$!!$      tracer(i,vseasj)  = 0.250        ! seasalt, j     ug/m3
!!$!!$      tracer(i,vseasc)  = 5.000        ! seasalt, c     ug/m3
!!$!!$
!!$!!$      tracer(i,vdustj)  = 1.000        ! dust, j        ug/m3
!!$!!$      tracer(i,vdustc)  = 4.000        ! dust, c        ug/m3
!!$!!$
!!$!!$!      tracer(i,veci)    = 0.0          ! BC, i          ug/m3
!!$!!$!      tracer(i,vecj)    = 0.0          ! BC, j          ug/m3
!!$!!$
!!$!!$!      tracer(i,vorgpai) = 0.0          ! POM, i         ug/m3
!!$!!$!      tracer(i,vorgpaj) = 0.0          ! POM, j         ug/m3
!!$!!$
!!$!!$      ! aerosol water is calculated by EQSAM
!!$!!$!      tracer(i,vh2oai)  = 0.000       ! H2O, i         ug/m3
!!$!!$!      tracer(i,vh2oaj)  = 0.000       ! H2O, j         ug/m3
!!$!!$!      tracer(i,vh2oac)  = 0.000       ! H2O, c         ug/m3
!!$!!$
!!$!!$   end do
!!$
!!$   tmst = 1800.0                      ! s
!!$
!!$   do i = 1, numcells
!!$
!!$      cloudcover(i)   = 0.0           ! frac
!!$      so4rat(i)       = 0.0           ! ug/(m3s)
!!$      pressure(i)     = 18941.19      ! Pa
!!$      rh_hist_akn(i)  = 2.0           ! 
!!$      rh_hist_acc(i)  = 2.0           ! 
!!$      rh_hist_cor(i)  = 2.0           ! 
!!$      soa(i)          = 0.0           ! ug/(m3s)
!!$      relhum(i)       = 1.050416      ! (0-1)
!!$      temperature(i)  = 213.6351      ! K
!!$
!!$! *** conversion from mol/mol(air) to ug/m3 ***
!!$      poRT     = pressure(i) / (R_gas * temperature(i))
!!$      h2so4fac = 1.e6 * mwh2so4 * poRT
!!$      hno3fac  = 1.e6 * mwhno3  * poRT
!!$      nh3fac   = 1.e6 * mwnh3   * poRT
!!$
!!$      ! i = Aitken mode
!!$      ! j = accumulation mode
!!$      ! c = coarse mode
!!$
!!$      do j = 1, nspcsda
!!$         tracer(i,j) = 0.0
!!$      end do
!!$
!!$      tracer(i,vsulf) = 2.355072e-15 * h2so4fac           ! H2SO4(g)     ug/m3
!!$      tracer(i,vnh3)  = 8.311068e-31 * nh3fac             ! NH3(g)       ug/m3
!!$      tracer(i,vhno3) = 1.393329e-09 * hno3fac            ! HNO3(g)      ug/m3
!!$
!!$      tracer(i,vnu0)    = 4.661840E+09    ! number, i      #/m3
!!$      tracer(i,vac0)    = 4.165791E+07    ! number, j      #/m3
!!$      tracer(i,vcorn)   = 1.204018E+05    ! number, c      #/m3
!!$
!!$      tracer(i,vso4ai)  = 5.439324E-04    ! SO4, i         ug/m3
!!$      tracer(i,vso4aj)  = 4.356200E-02    ! SO4, j         ug/m3
!!$
!!$      tracer(i,vnh4ai)  = 1.112929E-04    ! NH4, i         ug/m3
!!$      tracer(i,vnh4aj)  = 1.226079E-02    ! NH4, j         ug/m3
!!$
!!$      tracer(i,vno3ai)  = 1.301544E-17    ! NO3, i         ug/m3
!!$      tracer(i,vno3aj)  = 1.399418E-17    ! NO3, j         ug/m3
!!$
!!$      tracer(i,vseasj)  = 2.560123E-04    ! seasalt, j     ug/m3
!!$      tracer(i,vseasc)  = 4.286838E-04    ! seasalt, c     ug/m3
!!$
!!$      tracer(i,vdustj)  = 2.571087E-03    ! dust, j        ug/m3
!!$      tracer(i,vdustc)  = 5.212246E-02    ! dust, c        ug/m3
!!$
!!$      tracer(i,veci)    = 2.371810E-05    ! BC, i          ug/m3
!!$      tracer(i,vecj)    = 8.166940E-04    ! BC, j          ug/m3
!!$
!!$      tracer(i,vorgpai) = 3.591079E-05    ! POM, i         ug/m3
!!$      tracer(i,vorgpaj) = 9.088896E-03    ! POM, j         ug/m3
!!$
!!$      ! aerosol water is calculated by EQSAM
!!$      tracer(i,vh2oai)  = 5.764734E-22    ! H2O, i         ug/m3
!!$      tracer(i,vh2oaj)  = 4.562811E-02    ! H2O, j         ug/m3
!!$      tracer(i,vh2oac)  = 4.823098E-04    ! H2O, c         ug/m3
!!$
!!$   end do
!!$! op_ck_20120320-
!!$
!!$   ! ----------------------------------------------------------------------
!!$
!!$!   so4sum0=(tracer(1,vso4ai)+tracer(1,vso4aj))/mwso4+          &
!!$!           (tracer(1,vsulf) +timesteps*tmst*so4rat(1))/mwh2so4
!!$!   nh4sum0=(tracer(1,vnh4ai)+tracer(1,vnh4aj))/mwnh4+          &
!!$!            tracer(1,vnh3)/mwnh3
!!$!   no3sum0=(tracer(1,vno3ai)+tracer(1,vno3aj))/mwno3+          &
!!$!            tracer(1,vhno3)/mwhno3
!!$!   bcsum0 =tracer(1,veci)   +tracer(1,vecj)
!!$!   pomsum0=tracer(1,vorgpai)+tracer(1,vorgpaj)
!!$!   sssum0 =tracer(1,vseasi) +tracer(1,vseasj)+tracer(1,vseasc)
!!$!   dusum0 =tracer(1,vdustj) +tracer(1,vdustc)
!!$
! op_ck_20120430-

   call made_init             ! Initialization of MADE before 1st call
! op_ck_20120430+
   call made_box_init         ! Read initial parameters and variables for box
                              ! model run from namelist file
! op_ck_20120430-

! mz_pj_20071008+
!   IF (lmade == .false.) STOP ! stop if MADE switched off
   IF (.NOT. lmade ) STOP ! stop if MADE switched off
! mz_pj_20071008-

   ! print some diagnostics of aerosol components...

   write (*,*)
   write (*,'(A10,3(2X,A9))') 'component','akn','acc','cor'
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'SO4',tracer(1,vso4ai),  &
             tracer(1,vso4aj),'-'
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'NH4',tracer(1,vnh4ai),  &
             tracer(1,vnh4aj),'-'
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'NO3',tracer(1,vno3ai),  &
             tracer(1,vno3aj),'-'
   write (*,'(A10,3(2X,e9.3))') 'H2O',tracer(1,vh2oai),        &
             tracer(1,vh2oaj),tracer(1,vh2oac)
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'BC',tracer(1,veci),     &
             tracer(1,vecj),'-'
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'POM',                   &
              tracer(1,vorgpai),tracer(1,vorgpaj),'-'
   write (*,'(A10,3(2X,e9.3))') 'seas',tracer(1,vseasi),       &
              tracer(1,vseasj),tracer(1,vseasc)
   write (*,'(A10,2X,A9,2(2X,e9.3))') 'dust','-',              &
              tracer(1,vdustj),tracer(1,vdustc)

   ! ...and gas-phase species

   write (*,*)
   write (*,'(A10,2X,e9.3)') 'NH3(g)',tracer(1,vnh3)
   write (*,'(A10,2X,e9.3)') 'H2SO4(g)',tracer(1,vsulf)
   write (*,'(A10,2X,e9.3)') 'HNO3(g)',tracer(1,vhno3)

   ! Initialize species masses time series output file
   open (o_species, action='WRITE', file='made_species.dat', form='FORMATTED')
   write (o_species, '(A)') '# This file contains the time [s] evolution of the ' &
        // 'individual species mass concentrations [ug m-3]'
   write (o_species, *) ''
   write (o_species, '(A5)', advance='no') 'time'
   write (o_species, '(1X,A9)', advance='no') 'SO4'
   write (o_species, '(1X,A9)', advance='no') 'NH4'
   write (o_species, '(1X,A9)', advance='no') 'NO3'
   write (o_species, '(1X,A9)', advance='no') 'SS'
   write (o_species, '(1X,A9)', advance='no') 'POM'
   write (o_species, '(1X,A9)', advance='no') 'BC'
   write (o_species, '(1X,A9)', advance='no') 'DU'
   write (o_species, '(1X,A9)', advance='no') 'H2O'
   write (o_species, *) ''
   write (o_species, '(i5)', advance='no') 0
   write (o_species, '(1X,es9.3)', advance='no') tracer(1,vso4ai) &
                                                 + tracer(1,vso4aj)
   write (o_species, '(1X,es9.3)', advance='no') tracer(1,vnh4ai) &
                                                 + tracer(1,vnh4aj)
   write (o_species, '(1X,es9.3)', advance='no') tracer(1,vno3ai) &
                                                 + tracer(1,vno3aj)
   write (o_species, '(1X,es9.3)', advance='no') tracer(1,vseasi) &
                                                 + tracer(1,vseasj) &
                                                 + tracer(1,vseasc)
   write (o_species, '(1X,es9.3)', advance='no') tracer(1,vorgpai) &
                                                 + tracer(1,vorgpaj)
   write (o_species, '(1X,es9.3)', advance='no') tracer(1,veci) &
                                                 + tracer(1,vecj)
   write (o_species, '(1X,es9.3)', advance='no') tracer(1,vdustj) &
                                                 + tracer(1,vdustc)
   write (o_species, '(1X,es9.3)', advance='no') tracer(1,vh2oai) &
                                                 + tracer(1,vh2oaj) &
                                                 + tracer(1,vh2oac)
   write (o_species, *) ''

   ! Initialize dry diameter/mass output file
   open (o_dry, action='WRITE', file='made_dry.dat', form='FORMATTED')
   write (o_dry, '(A)') '# This file contains the time [s] evolution of the ' &
        // 'dry diameters [m] and modal dry mass concentrations [ug m-3]'
   write (o_dry, *) ''
   write (o_dry, '(A5)', advance='no') 'time'
   write (o_dry, '(1X,A9)', advance='no') 'dg, akn'
   write (o_dry, '(1X,A9)', advance='no') 'dg, acc'
   write (o_dry, '(1X,A9)', advance='no') 'dg, cor'
   write (o_dry, '(1X,A9)', advance='no') 'm, akn'
   write (o_dry, '(1X,A9)', advance='no') 'm, acc'
   write (o_dry, '(1X,A9)', advance='no') 'm, cor'
   write (o_dry, *) ''
   write (o_dry, '(i5)', advance='no') 0
   write (o_dry, '(1X,es9.3)', advance='no') MAX(dgmin, &
        ((tracer(1,vnu3) - tracer(1,vh2oai) * f6dpim9 / rhoh2o) &
         / (tracer(1,vnu0) * es36(akn)))**one3)
   write (o_dry, '(1X,es9.3)', advance='no') MAX(dgmin, &
        ((tracer(1,vac3) - tracer(1,vh2oaj) * f6dpim9 / rhoh2o) &
         / (tracer(1,vac0) * es36(acc)))**one3)
   write (o_dry, '(1X,es9.3)', advance='no') MAX(dgmin, &
        ((tracer(1,vcor3) - tracer(1,vh2oac) * f6dpim9 / rhoh2o) &
         / (tracer(1,vcorn) * es36(cor)))**one3)
   write (o_dry, '(1X,es9.3)', advance='no') SUM(tracer(1,(/vso4ai,vnh4ai,&
        vno3ai,vseasi,vorgpai,veci/)))
   write (o_dry, '(1X,es9.3)', advance='no') SUM(tracer(1,(/vso4aj,vnh4aj,&
        vno3aj,vseasj,vorgpaj,vecj,vdustj/)))
   write (o_dry, '(1X,es9.3)', advance='no') SUM(tracer(1,(/vseasc,vdustc/)))
   write (o_dry, *) ''

   ! Now begin writing time series output

   write (*,*)
   write (*,'(A5,6(2X,A9))') 'step','num, akn','num, acc', &
                             'num, cor','dg, akn','dg, acc','dg, cor'

! op_ck_20120320+
   ! Write initial number conc.s and diameters to standard output
   write (*,'(i5,6(2x,e9.3))') 0, tracer(1,vnu0),     &
        tracer(1,vac0), tracer(1,vcorn),              &
! op_ck_20120430+
        dg(akn), dg(acc), dg(cor)
! op_ck_20120430-
! op_ck_20120320-

! op_ck_20120430+
!!$
!!$   ! File for test output
!!$   if (ltest_mcon .or. ltest_ncon .or. ltest_so4 .or. ltest_adapdt) then
!!$      open(testout, action='WRITE', file='made_test.log', form='FORMATTED')
!!$   end if
!!$
!!$   ! *** For mass conservation test
!!$   if (ltest_mcon) then
!!$      oldmass = 0.
!!$      do j = 1,nummod
!!$         do k = 2, 7
!!$            oldmass = oldmass + aero(1,j,k)%p
!!$         end do
!!$      end do
!!$      oldmass = oldmass + tracer(1,vhno3) * mwno3 / mwhno3 &
!!$              + tracer(1,vnh3) * mwnh4 / mwnh3
!!$      write (testout,'(A40,2X,es13.6)') 'total initial dry mass conc. ' &
!!$           // '(incl. gases)', oldmass
!!$   end if
!!$
!!$   ! *** For number consistency test
!!$   if (ltest_ncon) then
!!$      oldnum = sum(tracer(1, (/ vnu0, vac0, vcorn /) ))
!!$      write (testout,'(A40,2X,es13.6)') 'total initial number conc.', &
!!$           oldnum
!!$   end if
!!$
!!$   ! *** For SO4 production test
!!$   if (ltest_so4) then
!!$      oldso4 = tracer(1,vsulf) * mwso4 / mwh2so4 &
!!$           + sum(tracer(1, (/ vso4ai, vso4aj /) ))
!!$   end if
!!$
! op_ck_20120430-

   do i = 1, timesteps

      ! *** Emissions ***
      tracer(1,vsulf) = tracer(1,vsulf) + so4rat(1) * tmst
      tracer(1,vhno3) = tracer(1,vhno3) + no3rat(1) * tmst
      tracer(1,veci)  = tracer(1,veci)  + em_bc_mass(1,1) * tmst
      tracer(1,vecj)  = tracer(1,vecj)  &
                        + (em_bc_mass(2,1) + em_bc_mass(3,1)) * tmst
      tracer(1,vnu0)  = tracer(1,vnu0)  + em_bc_num(1,1) * tmst
      tracer(1,vac0)  = tracer(1,vac0)  &
                        + (em_bc_num(2,1) + em_bc_num(3,1)) * tmst

      call made_main( BLKSIZE, NUMCELLS, PRESSURE, TEMPERATURE, RELHUM,  &
                      TMST, SO4RAT, SOA, CLOUDCOVER, TRACER,             &
                      RH_HIST_AKN, RH_HIST_ACC, RH_HIST_COR,             &
! op_ck_20120430+
!!$                      DGNUC, DGACC, DGCOR, DGDRYNUC, DGDRYACC, DGDRYCOR, &
                      DG(AKN), DG(ACC), DG(COR),                         &
                      DGDRYNUC, DGDRYACC, DGDRYCOR,                      &
! op_ck_20120430-
                      DENSN, DENSA, DENSC )

      write (*,'(i5,6(2x,e9.3))') i, tracer(1,vnu0),            &
                  tracer(1,vac0), tracer(1,vcorn),              &
! op_ck_20120430+
!!$                  dgnuc(1), dgacc(1), dgcor(1)
                  dg(akn), dg(acc), dg(cor)
! op_ck_20120430-

      ! Write time series of species mass concentrations
      write (o_species, '(i5)', advance='no') int(i * tmst)
      write (o_species, '(1X,es9.3)', advance='no') tracer(1,vso4ai) &
                                                    + tracer(1,vso4aj)
      write (o_species, '(1X,es9.3)', advance='no') tracer(1,vnh4ai) &
                                                    + tracer(1,vnh4aj)
      write (o_species, '(1X,es9.3)', advance='no') tracer(1,vno3ai) &
                                                    + tracer(1,vno3aj)
      write (o_species, '(1X,es9.3)', advance='no') tracer(1,vseasi) &
                                                    + tracer(1,vseasj) &
                                                    + tracer(1,vseasc)
      write (o_species, '(1X,es9.3)', advance='no') tracer(1,vorgpai) &
                                                    + tracer(1,vorgpaj)
      write (o_species, '(1X,es9.3)', advance='no') tracer(1,veci) &
                                                    + tracer(1,vecj)
      write (o_species, '(1X,es9.3)', advance='no') tracer(1,vdustj) &
                                                    + tracer(1,vdustc)
      write (o_species, '(1X,es9.3)', advance='no') tracer(1,vh2oai) &
                                                    + tracer(1,vh2oaj) &
                                                    + tracer(1,vh2oac)
      write (o_species, *) ''

      ! Write time series of dry diameters/masses
      write (o_dry, '(i5)', advance='no') int(i * tmst)
      write (o_dry, '(1X,es9.3)', advance='no') MAX(dgmin, &
           ((tracer(1,vnu3) - tracer(1,vh2oai) * f6dpim9 / rhoh2o) &
           / (tracer(1,vnu0) * es36(akn)))**one3)
      write (o_dry, '(1X,es9.3)', advance='no') MAX(dgmin, &
           ((tracer(1,vac3) - tracer(1,vh2oaj) * f6dpim9 / rhoh2o) &
           / (tracer(1,vac0) * es36(acc)))**one3)
      write (o_dry, '(1X,es9.3)', advance='no') MAX(dgmin, &
           ((tracer(1,vcor3) - tracer(1,vh2oac) * f6dpim9 / rhoh2o) &
           / (tracer(1,vcorn) * es36(cor)))**one3)
      write (o_dry, '(1X,es9.3)', advance='no') SUM(tracer(1,(/vso4ai,vnh4ai,&
           vno3ai,vseasi,vorgpai,veci/)))
      write (o_dry, '(1X,es9.3)', advance='no') SUM(tracer(1,(/vso4aj,vnh4aj,&
           vno3aj,vseasj,vorgpaj,vecj,vdustj/)))
      write (o_dry, '(1X,es9.3)', advance='no') SUM(tracer(1,(/vseasc,vdustc/)))
      write (o_dry, *) ''

! op_ck_20120430+
!!$
!!$      ! *** Mass conservation test
!!$      if (ltest_mcon) then
!!$         newmass = 0.
!!$         do j = 1, nummod
!!$            do k = 2, 7
!!$               newmass = newmass + aero(1,j,k)%p
!!$            end do
!!$         end do
!!$         newmass = newmass + tracer(1,vhno3) * mwno3 / mwhno3 &
!!$              + tracer(1,vnh3) * mwnh4 / mwnh3
!!$         write (testout,'(A36,1X,i3,2X,f13.9,A)') 'relative mass conc. ' &
!!$              // 'change step', i, (newmass - oldmass) * 100. / oldmass, '%'
!!$         write (testout,'(A36,1X,i3,2X,es13.6)') 'total mass conc. step', &
!!$              i, newmass
!!$         oldmass = newmass
!!$      end if
!!$
!!$      ! *** Number consistency test
!!$      if (ltest_ncon) then
!!$         newnum = sum(tracer(1, (/ vnu0, vac0, vcorn /) ))
!!$         write (testout,'(A36,1X,i3,2X,f13.9,A)') 'relative number conc. ' &
!!$              // 'change step', i, (newnum - oldnum) * 100. / oldnum, '%'
!!$         write (testout,'(A36,1X,i3,2X,es13.6)') 'total number conc. step', &
!!$              i, newnum
!!$         oldnum = newnum
!!$      end if
!!$
!!$      ! *** SO4 production test
!!$      if (ltest_so4) then
!!$         newso4 = tracer(1,vsulf) * mwso4 / mwh2so4 &
!!$              + sum(tracer(1, (/ vso4ai, vso4aj /) ))
!!$         write (testout,'(A36,1X,i3,2X,f13.9,A)') 'relative SO4 prod. ' &
!!$              // 'deviation step', i, ((newso4 - oldso4) * 100. &
!!$              / (so4rat * oldtmst * mwso4 / mwh2so4)) - 100., '%'
!!$         write (testout,'(A36,1X,i3,2X,es13.6)') 'absolute SO4 prod. step', &
!!$              i, (newso4 - oldso4) / oldtmst
!!$         oldso4 = newso4
!!$      end if
!!$
!!$      ! *** Adaptive time step test
!!$      if (ltest_adapdt) then
!!$         write(testout,'(A36,1X,i3,2X,es13.6,A)') 'actual step duration ' &
!!$              // 'during step', i, tmst, 's'
!!$      end if
!!$      tmst = oldtmst
!!$
! op_ck_20120430-

   end do

   ! print some diagnostics of aerosol components...

   write (*,*)
   write (*,'(A10,3(2X,A9))') 'component','akn','acc','cor'
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'SO4',tracer(1,vso4ai),  &
             tracer(1,vso4aj),'-'
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'NH4',tracer(1,vnh4ai),  &
             tracer(1,vnh4aj),'-'
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'NO3',tracer(1,vno3ai),  &
             tracer(1,vno3aj),'-'
   write (*,'(A10,3(2X,e9.3))') 'H2O',tracer(1,vh2oai),        &
             tracer(1,vh2oaj),tracer(1,vh2oac)
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'BC',tracer(1,veci),     &
             tracer(1,vecj),'-'
   write (*,'(A10,2(2X,e9.3),2X,A9)') 'POM',                   &
              tracer(1,vorgpai),tracer(1,vorgpaj),'-'
   write (*,'(A10,3(2X,e9.3))') 'seas',tracer(1,vseasi),       &
              tracer(1,vseasj),tracer(1,vseasc)
   write (*,'(A10,2X,A9,2(2X,e9.3))') 'dust','-',              &
              tracer(1,vdustj),tracer(1,vdustc)

   ! ...and gas-phase species

   write (*,*)
   write (*,'(A10,2X,e9.3)') 'NH3(g)',tracer(1,vnh3)
   write (*,'(A10,2X,e9.3)') 'H2SO4(g)',tracer(1,vsulf)
   write (*,'(A10,2X,e9.3)') 'HNO3(g)',tracer(1,vhno3)

   call made_finalize   ! Finalization of MADE after last call

   close (o_species) ! File for mass concentrations' time series
   close (o_dry)     ! File for dry diameter/mass time series

!   so4sum1=(tracer(1,vso4ai)+tracer(1,vso4aj))/mwso4+          &
!            tracer(1,vsulf)/mwh2so4
!   nh4sum1=(tracer(1,vnh4ai)+tracer(1,vnh4aj))/mwnh4+          &
!            tracer(1,vnh3)/mwnh3
!   no3sum1=(tracer(1,vno3ai)+tracer(1,vno3aj))/mwno3+          &
!            tracer(1,vhno3)/mwhno3
!   bcsum1 =tracer(1,veci)   +tracer(1,vecj)
!   pomsum1=tracer(1,vorgpai)+tracer(1,vorgpaj)
!   sssum1 =tracer(1,vseasi) +tracer(1,vseasj)+tracer(1,vseasc)
!   dusum1 =tracer(1,vdustj) +tracer(1,vdustc)
!
!   write (*,*)
!   write (*,*) 'Delta SO4 =',(so4sum1-so4sum0)*100.0/so4sum0,'%'
!   write (*,*) 'Delta NH4 =',(nh4sum1-nh4sum0)*100.0/nh4sum0,'%'
!   write (*,*) 'Delta NO3 =',(no3sum1-no3sum0)*100.0/no3sum0,'%'
!   write (*,*) 'Delta BC  =',(bcsum1 -bcsum0) *100.0/bcsum0, '%'
!   write (*,*) 'Delta POM =',(pomsum1-pomsum0)*100.0/pomsum0,'%'
!   write (*,*) 'Delta SS  =',(sssum1 -sssum0) *100.0/sssum0, '%'
!   write (*,*) 'Delta DU  =',(dusum1 -dusum0) *100.0/dusum0, '%'

! op_ck_20120430+
!!$   if (ltest_mcon .or. ltest_ncon .or. ltest_so4 .or. ltest_adapdt) &
!!$        close(testout)
! op_ck_20120430-

   write (*,*)
   write (*,*) 'BoxMADE run finished'

END PROGRAM made

!*****************************************************************************
