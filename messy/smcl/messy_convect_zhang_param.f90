MODULE messy_convect_zhang_param

  USE messy_main_constants_mem,   ONLY: grav => g , R_gas, M_air, &
                                        M_H2O, dp, Tmelt

  IMPLICIT NONE

  SAVE

  ! Parameters and constants
  REAL(dp),PARAMETER :: cp    = 1005.46_dp  ! specific heat of dry air at
                                            ! constant pressure in J/K/kg
 ! H2O related constants, liquid density, phase change constants
  REAL(dp),PARAMETER :: rhoh2o = 1.0e3_dp    ! density of liquid water (STP)
  REAL(dp),PARAMETER :: hlat   = 2.5008e6_dp ! latent heat for vaporisation in J/kg
  REAL(dp),PARAMETER :: hlatv  = 2.8345e6_dp ! latent heat for sublimation in J/kg

! reciprocal of parameters: 
  REAL(dp),PARAMETER :: rhlat = 1._dp / hlat
  REAL(dp),PARAMETER :: rcp   = 1._dp / cp
  REAL(dp),PARAMETER :: rgrav = 1._dp / grav

  REAL(dp),PARAMETER :: rgas  = 1000._dp * R_gas / M_air  ! gas constant for dr air in J/K/kg
  REAL(dp),PARAMETER :: rgasv = 461.51_dp   ! gas constant for water vapour ! in J/K/kg


! ----------------- common block comcmf -------------------------------------   
!  REAL(dp) :: cp          ! specific heat of dry air    , stands above as parameter
!  REAL(dp) :: hlat        ! latent heat of vaporization , stands above as parameter
!  REAL(dp) :: grav        ! gravitational constant      , stands above as parameter
  REAL(dp) :: c0          ! rain water autoconversion coefficient
  REAL(dp) :: betamn      ! minimum overshoot parameter
!  REAL(dp) :: rhlat       ! reciprocal of hlat          , stands above as parameter
!  REAL(dp) :: rcp         ! reciprocal of cp            , stands above as parameter
!  REAL(dp) :: rgrav       ! reciprocal of grav          , stands above as parameter
  REAL(dp) :: cmftau      ! characteristic adjustment time scale
!  REAL(dp) :: rhoh2o      ! density of liquid water (STP),stands above as parameter
!  REAL(dp) :: rgas        ! gas constant for dry air    , stands above as parameter
  REAL(dp) :: dzmin       ! minimum convective depth for precipitation
  REAL(dp) :: small        ! arbitrary small num used in transport estimates
  REAL(dp) :: eps         ! convergence criteria (machine dependent)
  REAL(dp) :: tpmax       ! maximum acceptable t perturbation (degrees C)
  REAL(dp) :: shpmax      ! maximum acceptable q perturbation (g/g)

  INTEGER  :: iloc        ! longitude location for diagnostics
  INTEGER  :: jloc        ! latitude  location for diagnostics
  INTEGER  :: nsloc       ! nstep for which to produce diagnostics
  INTEGER  :: limcnv      ! top interface level limit for convection

  LOGICAL  :: rlxclm      ! logical to relax column versus cloud triplet

! -------------------- common block guang -------------------------------------
  REAL(dp) :: a
  REAL(dp) :: b
  REAL(dp) :: eps1
  REAL(dp) :: c1
  REAL(dp) :: c2
  REAL(dp) :: c3
  REAL(dp) :: tfreez
!-------------------- additional guang_pjr ------------------------------------
  real(dp) :: cpg
  real(dp) :: tau
  real(dp) :: dd
  real(dp) :: rl                   ! latent heat of vap
  integer  :: msg
  real(dp) :: qmin
  logical  :: trigon
#ifdef DEBCONV
  integer  :: latlook
  integer  :: ilook
#endif

! ---------------------include file from lookup tables ---------------------------

! $Id: eslookup.h,v 1.1.1.1.2.1 1999/01/15 00:35:09 eaton Exp $

! Common block and statement functions for saturation vapor pressure
! look-up procedure, J. J. Hack, February 1990

  integer  :: plenest  ! length of saturation vapor pressure table
  parameter (plenest=250)

! Table of saturation vapor pressure values es from tmin degrees
! to tmax+1 degrees k in one degree increments.  ttrice defines the
! transition region where es is a combination of ice & water values

!      common/comes0/estbl(plenest),tmin,tmax,ttrice,pcf(6),
!     $             epsqs,rgasv,hlatf,hlatv,cp,icephs

  real(dp) :: estbl(plenest), & ! table values of saturation vapor pressure

              tmin,  & ! min temperature (K) for table
              tmax,  & ! max temperature (K) for table
              ttrice,& ! transition range from es over water to es over ice
              pcf(6),   & ! polynomial coeffs -> es transition water to ice
              epsqs  ! Ratio of h2o to dry air molecular weights 
                           
  logical  :: icephs   ! false => saturation vapor pressure over water only


CONTAINS
!========================================================================================

  subroutine conv_ini (plev, plevp, hypi, xtrigon,  p_parallel_io )

!     a start on a subroutine to initialize the constants in guangs code

!-----------------------------------------------------------------------
!-----------------------------Arguments---------------------------------
    IMPLICIT NONE
      INTEGER, INTENT(IN) :: plev, plevp  ! number of levels, levels+1


      real(dp) :: hypi(plevp)          ! reference pressures at interfaces
      logical  :: xtrigon              ! true => triggers turned on
      logical, INTENT(IN)  :: p_parallel_io        ! input/output processor

      integer  :: k

      cpg = cp

      trigon = xtrigon
!
      a = 21.656_dp
      b = 5418._dp
      c1 = 6.112_dp
      c2 = 17.67_dp
      c3 = 243.5_dp
      eps1 = 0.622_dp
      qmin = 1.E-20_dp
      tfreez = 273.16_dp
      rl = 2.5104E6_dp

! Limit convection to regions below 50 mb

      if (hypi(1) .ge. 5.e3_dp) then
        limcnv = 1
      else
        do k=1,plev
          if (hypi(k).lt.5.e3_dp .and. hypi(k+1).ge.5.e3_dp) then
            limcnv = k
            goto 10
          end if
        end do
        limcnv = plevp
      end if
   10 continue

! Set internal variable "msg" (convection limit) to "limcnv-1"

      msg = limcnv - 1

! tau=4800. were used in canadian climate center. however, when it
! is used here in echam3 t42, convection is too weak, thus 
! adjusted to 2400. i.e the e-folding time is 1 hour now.

      tau = 7200._dp
      dd = 0.1_dp ! fraction of updraft mass flux going into downdraft
!      dd = 0.3_dp ! fraction of updraft mass flux going into downdraft
      if (p_parallel_io) then
        write (6,*) ' zhang/mcfarlane convection initialization '
        write (6,*) ' relaxation time tau is set to ', tau
        write (6,*) ' downdraft param dd is set to ', dd
        write (6,*) ' top index for zhang convection is ', msg
        write (6,*) ' convective trigger is set ', trigon
        write (6,*) ' done '
      endif
#ifdef DEBCONV
#ifdef LINUX
      latlook = 1
      ilook = 1
#else
      latlook = 1
      ilook = 69
!      ilook = 1
!      latlook = .true.
!      latlook = .false.
#endif
         write (6,*) ' debuging latlook, ilook ', latlook, ilook
#endif
      return
    end subroutine conv_ini
!========================================================================

    subroutine mfinti(plev, plevp, hypi, p_parallel_io  )

! ----------------------------------------------------------------------

! Initialize moist convective mass flux procedure common block, cmfmca 
 
! ---------------------------Code History-------------------------------

! Original version:  J. Hack
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, G. Taylor, August 1992

! ----------------------------------------------------------------------

! $Id: conv_hack.F,v 1.1.1.1.2.1 1998/11/06 16:16:43 eaton Exp $


! -----------------------------Arguments--------------------------------
   IMPLICIT NONE
      INTEGER, INTENT(IN) ::    &
     plev, plevp             ! number of vertical levels / # +1

! Input arguments
     
! +match
      REAL(dp) :: hypi(plevp)       ! reference pressures at interfaces (Pa)
! -match
      integer  ::  k             ! vertical level index
      logical, INTENT(IN)  ::  p_parallel_io ! input/output processor
! ----------------------------------------------------------------------

! Initialize free parameters for moist convective mass flux procedure

!  original value
!      c0     = 5.0e-5_dp        ! rain water autoconversion coeff (1/m)
!      c0     = 1.0e-4_dp        ! rain water autoconversion coeff (1/m) T21L19
!      c0     = 5.0e-4_dp        ! rain water autoconversion coeff (1/m) testing
      c0     = 1.0e-1_dp        ! rain water autoconversion coeff (1/m) T63L87 nudged

      dzmin  = 0.0_dp           ! minimum cloud depth to precipitate (m)
      betamn = 0.10_dp          ! minimum overshoot parameter
      cmftau = 3600._dp         ! characteristic adjustment time scale

! Limit convection to regions below 50 mb

      if (hypi(1) .ge. 5.e3_dp) then
        limcnv = 1
      else
        do k=1,plev
          if (hypi(k).lt.5.e3_dp .and. hypi(k+1).ge.5.e3_dp) then
            limcnv = k
            goto 10
          end if
        end do
        limcnv = plevp
      end if
10    if (p_parallel_io) &
        write(6,*)'MFINTI: Convection will be capped at intfc ',limcnv,    &
                  ' which is ',hypi(limcnv),' pascals'
      
      tpmax  = 1.50_dp          ! maximum acceptable t perturbation (deg C)
      shpmax = 1.50e-3_dp       ! maximum acceptable q perturbation (g/g)
      rlxclm = .true.        ! logical variable to specify that relaxation
!                                time scale should applied to column as 
!                                opposed to triplets individually     

! Initialize diagnostic location information for moist convection scheme

      iloc   = 1             ! longitude point for diagnostic info
      jloc   = 1             ! latitude  point for diagnostic info
      nsloc  = 1             ! nstep value at which to begin diagnostics

! Initialize other miscellaneous parameters

      small  = 1.0e-36_dp       ! arbitrary small number (scalar transport)
      eps    = 1.0e-13_dp       ! convergence criteria (machine dependent) 
!
      return
    end subroutine mfinti
!=====================================================================================

    subroutine esinti(xip, ierr , p_parallel_io)
!-----------------------------------------------------------------------

! Initialize es lookup tables 

!---------------------------Code history--------------------------------

! Original version:  J. Hack
! Standardized:      L. Buja, June 1992
! Reviewed:          J. Hack, G. Taylor, August 1992

!-----------------------------------------------------------------------
      implicit none

!------------------------------Arguments--------------------------------

! Input arguments

                  
! +bee
      logical  :: xip             ! Ice phase (true or false)
! -bee
      integer  :: ierr            ! error status
!---------------------------Local workspace-----------------------------

      real(dp) ::  tmn,     &     ! Minimum temperature entry in table
                   tmx,     &     ! Maximum temperature entry in table
                   trice          ! Trans range from es over h2o to es over ice
      logical  ::  ip             ! Ice phase (true or false)
      logical, INTENT(IN)  ::  p_parallel_io  ! input/output processor
!-----------------------------------------------------------------------

! Specify control parameters first

      tmn   = 173.16_dp
      tmx   = 375.16_dp
      trice =  20.00_dp
! +bee
      ip    = xip
! -bee

! Call gestbl to build saturation vapor pressure table.

      call gestbl(tmn     ,tmx     ,trice   ,ip      , ierr, p_parallel_io )

      return
    end subroutine esinti
!===================================================================================

     subroutine gestbl(tmn     ,tmx     ,trice   ,ip  ,  ierr, p_parallel_io  )
!-----------------------------------------------------------------------

! Builds saturation vapor pressure table for later lookup procedure.
! Uses Goff & Gratch (1946) relationships to generate the table
! according to a set of free parameters defined below.  Auxiliary
! routines are also included for making rapid estimates (well with 1%)
! of both es and d(es)/dt for the particular table configuration.

!---------------------------Code history--------------------------------

! Original version:  J. Hack
! Standardized:      L. Buja, June 1992
! Reviewed:          J. Hack, G. Taylor, August 1992

!-----------------------------------------------------------------------
      implicit none

!------------------------------Arguments--------------------------------

! Input arguments

      real(dp) :: tmn,     &      ! Minimum temperature entry in es lookup table
                  tmx,     &      ! Maximum temperature entry in es lookup table
                  trice           ! Transition range from es over range to es over ice               
      logical, INTENT(IN)  ::  p_parallel_io  ! input/output processor            
      INTEGER  :: ierr            ! error status flag
!---------------------------Local variables-----------------------------

      real(dp) :: t           ! Temperature
      integer  :: n,      &   ! Increment counter
                  lentbl, &   ! Calculated length of lookup table
                  itype       ! Ice phase: 0 -> no ice phase
                         !            1 -> ice phase, no transition
                         !           -x -> ice phase, x degree transition
      logical  :: ip         ! Ice phase logical flag

      INTRINSIC :: NINT, ABS
!-----------------------------------------------------------------------
      ierr = 0
! Set es table parameters

      tmin   = tmn       ! Minimum temperature entry in table
      tmax   = tmx       ! Maximum temperature entry in table
      ttrice = trice     ! Trans. range from es over h2o to es over ice
      icephs = ip        ! Ice phase (true or false)

! Set physical constants required for es calculation

      epsqs  = M_H2O / M_air
            

      lentbl = NINT(tmax-tmin+2.000001_dp)
      if (lentbl .gt. plenest) then
         write(6,9000) tmax, tmin, plenest
         ierr = 1 
         return
      end if

! Begin building es table.
! Check whether ice phase requested.
! If so, set appropriate transition range for temperature

      if (icephs) then
!!         if(ttrice.ne.0.0_dp) then
        if(ABS(ttrice).lt.0.0_dp) then
            itype = nint(-ttrice)
         else
            itype = 1
         end if
      else
         itype = 0
      end if

      t = tmin - 1.0_dp
      do n=1,lentbl
         t = t + 1.0_dp
         call gffgch(t,estbl(n),itype, ierr)
         if (ierr.eq.1) return
      end do

      do n=lentbl+1,plenest
         estbl(n) = -99999.0_dp
      end do

! Table complete -- Set coefficients for polynomial approximation of
! difference between saturation vapor press over water and saturation
! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
! is valid in the range -40 < t < 0 (degrees C).

!                  --- Degree 5 approximation ---

      pcf(1) =  5.04469588506e-01_dp
      pcf(2) = -5.47288442819e+00_dp
      pcf(3) = -3.67471858735e-01_dp
      pcf(4) = -8.95963532403e-03_dp
      pcf(5) = -7.78053686625e-05_dp

!                  --- Degree 6 approximation ---

!-----pcf(1) =  7.63285250063e-02_dp
!-----pcf(2) = -5.86048427932e+00_dp
!-----pcf(3) = -4.38660831780e-01_dp
!-----pcf(4) = -1.37898276415e-02_dp
!-----pcf(5) = -2.14444472424e-04_dp
!-----pcf(6) = -1.36639103771e-06_dp
      if (p_parallel_io)    &
      write(6,*)' ***** SATURATION VAPOR PRESSURE TABLE COMPLETED *****'
      return

 9000 format('GESTBL: FATAL ERROR *********************************',/,    &
             ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH',    &
             ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/,   &
             ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)

    end subroutine gestbl
!===================================================================================

      subroutine gffgch(t       ,es      ,itype, ierr )
!-----------------------------------------------------------------------

! Computes saturation vapor pressure over water and/or over ice using
! Goff & Gratch (1946) relationships.  T (temperature), and itype are
! input parameters, while es (saturation vapor pressure) is an output
! parameter.  The input parameter itype serves two purposes: a value of
! zero indicates that saturation vapor pressures over water are to be
! returned (regardless of temperature), while a value of one indicates
! that saturation vapor pressures over ice should be returned when t is
! less than 273.16 degrees k.  If itype is negative, its absolute value
! is interpreted to define a temperature transition region below 273.16
! degrees k in which the returned saturation vapor pressure is a
! weighted average of the respective ice and water value.  That is, in
! the temperature range 0 => -itype degrees c, the saturation vapor
! pressures are assumed to be a weighted average of the vapor pressure
! over supercooled water and ice (all water at 0 c; all ice at -itype
! c).  Maximum transition range => 40 c

!---------------------------Code history--------------------------------

! Original version:  J. Hack
! Standardized:      L. Buja, June 1992
! Reviewed:          J. Hack, G. Taylor, August 1992

!-----------------------------------------------------------------------
      implicit none

!------------------------------Arguments--------------------------------

! Input arguments

      real(dp) :: t          ! Temperature
      integer  :: itype      ! Flag for ice phase and associated transition

! Output arguments

      real(dp) ::  es         ! Saturation vapor pressure

!---------------------------Local variables-----------------------------

      real(dp) :: e1        ! Intermediate scratch variable for es over water
      real(dp) :: e2        ! Intermediate scratch variable for es over water
      real(dp) :: eswtr     ! Saturation vapor pressure over water
      real(dp) :: f         ! Intermediate scratch variable for es over water
      real(dp) :: f1        ! Intermediate scratch variable for es over water
      real(dp) :: f2        ! Intermediate scratch variable for es over water
      real(dp) :: f3        ! Intermediate scratch variable for es over water
      real(dp) :: f4        ! Intermediate scratch variable for es over water
      real(dp) :: f5        ! Intermediate scratch variable for es over water
      real(dp) :: ps        ! Reference pressure (mb)
      real(dp) :: t0        ! Reference temperature (freezing point of water)
      real(dp) :: term1     ! Intermediate scratch variable for es over ice
      real(dp) :: term2     ! Intermediate scratch variable for es over ice
      real(dp) :: term3     ! Intermediate scratch variable for es over ice
      real(dp) :: tr        ! Transition range for es over water to es over ice
      real(dp) :: ts        ! Reference temperature (boiling point of water)
      real(dp) :: weight    ! Intermediate scratch variable for es transition
      integer  :: itypo     ! Intermediate scratch variable for holding itype
      INTEGER  :: ierr

      INTRINSIC :: abs, exp, REAL, LOG, LOG10, MIN 
!-----------------------------------------------------------------------
      
      ierr = 0
! Check on whether there is to be a transition region for es

      if (itype.lt.0) then
         tr    = abs(real(itype,dp))
         itypo = itype
         itype = 1
      else
         tr    = 0.0_dp
         itypo = itype
      end if
      if (tr .gt. 40.0_dp) then
         write(6,900) tr
         ierr = 1
         return
      end if

      if(t .lt. (273.16_dp - tr) .and. itype.eq.1) go to 10

! Water

      ps = 1013.246_dp
      ts = 373.16_dp
      e1 = 11.344_dp*(1.0_dp - t/ts)
      e2 = -3.49149_dp*(ts/t - 1.0_dp)
      f1 = -7.90298_dp*(ts/t - 1.0_dp)
      f2 = 5.02808_dp*log10(ts/t)
      f3 = -1.3816_dp*(10.0_dp**e1 - 1.0_dp)/10000000.0_dp
      f4 = 8.1328_dp*(10.0_dp**e2 - 1.0_dp)/1000.0_dp
      f5 = log10(ps)
      f  = f1 + f2 + f3 + f4 + f5
      es = (10.0_dp**f)*100.0_dp
      eswtr = es

      if(t.ge.273.16_dp .or. itype.eq.0) go to 20

! Ice

   10 continue
      t0    = 273.16_dp
      term1 = 2.01889049_dp/(t0/t)
      term2 = 3.56654_dp*log(t0/t)
      term3 = 20.947031_dp*(t0/t)
      es    = 575.185606e10_dp*exp(-(term1 + term2 + term3))

      if (t.lt.(273.16_dp - tr)) go to 20

! Weighted transition between water and ice

      weight = min((273.16_dp - t)/tr,1.0_dp)
      es = weight*es + (1.0_dp - weight)*eswtr

   20 continue
      itype = itypo
      return

  900 format('GFFGCH: FATAL ERROR ******************************',/,   &
             'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR',     &
             ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF',      &
             ' 40.0 DEGREES C',/, ' TR = ',f7.2)

   end subroutine gffgch
 



!====================================================================================
     subroutine whenflt( n, array, inc, ltarget, lindex, nval )

! Finds all array elements less than the target.

      implicit none

! Input:

      integer  :: &
         n        & ! Number of elements to be searched
     ,   inc        ! Increment between elements of the searched array
      real(dp) :: &
         array(n) & ! Array to be searched
      ,  ltarget    ! Value searched for in the array

! Output:

      integer  ::  &
         lindex(n) & ! Array containing indices of array elements matching target
      ,  nval        ! Number of values put in index array

! Local:

      integer ina, i

      ina = 1
      nval = 0

      if ( inc .lt. 0 ) ina = -inc*(n-1) + 1

      do i = 1, n
         if ( array(ina) .lt. ltarget ) then
            nval = nval + 1
            lindex(nval) = ina
         endif
         ina = ina + inc
      enddo

      return
    end subroutine whenflt

!--------------------------------------------------------------------------------------------
    subroutine whenfgt(n,array,inc,ltarget,lindex,nval)

    ! Finds all array elements greater than the target.

      implicit none

! Input:

      integer  ::&
        n        & ! Number of elements to be searched
     ,  inc        ! Increment between elements of the searched array
      real(dp) ::&
        array(n) & ! Array to be searched
     ,  ltarget     ! Value searched for in the array

! Output:

      integer  ::  &
        lindex(n)  &! Array containing indices of array elements matching target
     ,  nval        ! Number of values put in index array

! Local:

      integer ina, i

      ina=1
      nval=0
      if(inc .lt. 0) ina=(-inc)*(n-1)+1
      do 100 i=1,n
         if(array(ina) .gt. ltarget) then
            nval=nval+1
            lindex(nval)=i
         end if
         ina=ina+inc
100    enddo
      return
    end subroutine whenfgt


!=========================================================================================

! lookup tables

      subroutine aqsatd(t       ,p       ,es      ,qs      ,gam     , &
                        ii      ,illen    ,kk      ,kstart  ,kend    )
!-----------------------------------------------------------------------
!
! Utility procedure to look up and return saturation vapor pressure from 
! precomputed table, calculate and return saturation specific humidity 
! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT) for input 
! arrays of temperature and pressure (dimensioned ii,kk).
!
!---------------------------Code history--------------------------------
!
! Original version: J. J. Hack, February 1990
! Reviewed:         J. J. Hack, August 1992 
!
!-----------------------------------------------------------------------
      implicit none

!------------------------------Arguments--------------------------------

! Input arguments

      integer  :: ii,  &         ! I dimension of arrays t, p, es, qs
                  kk             ! K dimension of arrays t, p, es, qs 
      real(dp) :: t(ii,kk),  &   ! Temperature
                  p(ii,kk)       ! Pressure
      integer  :: illen,     &   ! vector length in I direction
                  kstart,    &   ! starting location in K direction
                  kend           ! ending location in K direction
 
! Output arguments

      real(dp) :: es(ii,kk), &   ! Saturation vapor pressure
                  qs(ii,kk), &   ! Saturation specific humidity
                  gam(ii,kk)     ! (l/cp)*(d(qs)/dt)

!---------------------------Local workspace-----------------------------

      logical  :: lflg          ! true if in temperature transition region
      integer  :: i,  &         ! i index for vector calculations
                  k             ! k index 
      real(dp) :: omeps,  &     ! 1. - 0.622
                  trinv,  &     ! reciprocal of ttrice (transition range)
                  tc,     &     ! temperature (in degrees C)
                  weight, &     ! weight for es transition from water to ice
                  hltalt, &     ! appropriately modified hlat for T derivatives
                  hlatsb, &     ! hlat weighted in transition region
                  hlatvp, &     ! hlat modified for t changes above 273.16
                  tterm,  &     ! account for d(es)/dT in transition region
                  desdt         ! d(es)/dT

      INTRINSIC :: MIN, TINY, ABS

!-----------------------------------------------------------------------
      omeps = 1.0_dp - epsqs
      do k=kstart,kend
         do i=1,illen
            es(i,k) = estblf(t(i,k))

! Saturation specific humidity

            qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))

! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere

            qs(i,k) = min(1.0_dp,qs(i,k))

            if (qs(i,k) .lt. 0.0_dp) then
               qs(i,k) = 1.0_dp
               es(i,k) = p(i,k)
            end if
         end do
      end do

! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16

      trinv = 0.0_dp
!      if ((.not. icephs) .or. (ttrice.eq.0.0)) go to 10
! mz_ht_20040508+
! avoid goto statement
!!      if (icephs .or. (ttrice.ne.0.0_dp)) then
      if(icephs .or. ABS(ttrice).lt.0.0_dp) then   
        trinv = 1.0/ttrice

        do k=kstart,kend
          do i=1,illen

! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition 
! range from ice to water also accounting for change of hlatv with t 
! above 273.16 where constant slope is given by -2369 j/(kg c) =cpv - cw

            tc     = t(i,k) - 273.16_dp
            lflg   = (tc.ge.-ttrice .and. tc.lt.0.0_dp)
            weight = min(-tc*trinv,1.0_dp)
            hlatsb = hlatv + weight*hlat
            hlatvp = hlatv - 2369.0_dp*tc
            if (t(i,k).lt.273.16_dp) then
               hltalt = hlatsb
            else
               hltalt = hlatvp
            end if
            if (lflg) then
               tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) +       &
                       tc*pcf(5))))
            else
               tterm = 0.0_dp
            end if
            desdt  = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k)) + tterm*trinv
            gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/                          &
                      (cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
!!            if(qs(i,k).eq.1.0_dp) gam(i,k) = 0.0_dp
            if(ABS(qs(i,k)-1.0_dp).lt.TINY(0.0_dp)) gam(i,k) = 0.0_dp
          end do
        end do
!        goto 20
      ELSE
      

! No icephs or water to ice transition

        do k=kstart,kend   ! formerly 10
          do i=1,illen

! Account for change of hlatv with t above 273.16 where
! constant slope is given by -2369 j/(kg c) = cpv - cw

            hlatvp = hlatv - 2369.0_dp*(t(i,k)-273.16_dp)
            if (icephs) then
               hlatsb = hlatv + hlat
            else
               hlatsb = hlatv
            end if
            if (t(i,k).lt.273.16_dp) then
               hltalt = hlatsb
            else
               hltalt = hlatvp
            end if
            desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k))
            gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/                       &
                       (cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
!!            if (qs(i,k) .eq. 1.0_dp) gam(i,k) = 0.0_dp
            if(ABS(qs(i,k)-1.0_dp).lt.TINY(0.0_dp)) gam(i,k) = 0.0_dp 
         end do
        end do

      ENDIF
      
      return   ! formerly 20
    end subroutine aqsatd
!=======================================================================================

    subroutine vqsatd(t  ,p  ,es  ,qs  ,gam  ,len)
!-----------------------------------------------------------------------

! Utility procedure to look up and return saturation vapor pressure from 
! precomputed table, calculate and return saturation specific humidity 
! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
! function as qsatd, but operates on vectors of temperature and pressure

!----------------------------Code History-------------------------------

! Original version:  J. Hack
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, August 1992

!-----------------------------------------------------------------------
      implicit none
!------------------------------Arguments--------------------------------

! Input arguments

      integer  :: len       ! vector length
      real(dp) :: t(len), & ! temperature
                  p(len)    ! pressure
 
! Output arguments

      real(dp) ::  es(len), & ! saturation vapor pressure
                   qs(len), & ! saturation specific humidity
                   gam(len)   ! (l/cp)*(d(qs)/dt)

!---------------------------Local workspace-----------------------------

      logical  :: lflg       ! true if in temperature transition region
      integer  :: i          ! index for vector calculations
      real(dp) ::  omeps, &  ! 1. - 0.622
                   trinv, &  ! reciprocal of ttrice (transition range)
                   tc,    &  ! temperature (in degrees C)
                   weight,&  ! weight for es transition from water to ice
                   hltalt,&  ! appropriately modified hlat for T derivatives  
                   hlatsb,&  ! hlat weighted in transition region
                   hlatvp,&  ! hlat modified for t changes above 273.16
                   tterm, &  ! account for d(es)/dT in transition region
                   desdt     ! d(es)/dT

      INTRINSIC :: MIN, TINY, ABS
!-----------------------------------------------------------------------
      omeps = 1.0_dp - epsqs
      do i=1,len
        es(i) = estblf(t(i))

! Saturation specific humidity

        qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))

! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere

        qs(i) = min(1.0_dp,qs(i))

        if (qs(i) .lt. 0.0_dp) then
          qs(i) = 1.0_dp
          es(i) = p(i)
        end if
      end do

! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16

      trinv = 0.0_dp
!      if ((.not. icephs) .or. (ttrice.eq.0.0)) goto 10
! mz_ht_20040508+
! avoid goto statement
!!      if (icephs .or. (ttrice.ne.0.0_dp)) then
      if(icephs .or. ABS(ttrice).lt.0.0_dp) then   

        trinv = 1.0_dp/ttrice
        do i=1,len

! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition 
! range from ice to water also accounting for change of hlatv with t 
! above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw

          tc     = t(i) - 273.16_dp
          lflg   = (tc.ge.-ttrice .and. tc.lt.0.0_dp)
          weight = min(-tc*trinv,1.0_dp)
          hlatsb = hlatv + weight*hlat
          hlatvp = hlatv - 2369.0_dp*tc
          if (t(i).lt.273.16_dp) then
            hltalt = hlatsb
          else
            hltalt = hlatvp
          end if
          if (lflg) then
            tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) +        &
                    tc*pcf(5))))
          else
            tterm = 0.0_dp
          end if
          desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
          gam(i) = hltalt*qs(i)*p(i)*desdt/                                &
                   (cp*es(i)*(p(i) - omeps*es(i)))
!!          if(qs(i).eq.1.0_dp) gam(i) = 0.0_dp
          if(ABS(qs(i)-1.0_dp).lt.TINY(0.0_dp)) gam(i) = 0.0_dp
        end do
        
      ELSE

! No icephs or water to ice transition

        do i=1,len   ! formerly 10

! Account for change of hlatv with t above 273.16 where
! constant slope is given by -2369 j/(kg c) = cpv - cw

          hlatvp = hlatv - 2369.0_dp*(t(i)-273.16_dp)
          if (icephs) then
            hlatsb = hlatv + hlat
          else
            hlatsb = hlatv
          end if
          if (t(i).lt.273.16_dp) then
            hltalt = hlatsb
          else
            hltalt = hlatvp
          end if
          desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
          gam(i) = hltalt*qs(i)*p(i)*desdt/                                &
                   (cp*es(i)*(p(i) - omeps*es(i)))
!!          if (qs(i) .eq. 1.0_dp) gam(i) = 0.0_dp
          if(ABS(qs(i)-1.0_dp).lt.TINY(0.0_dp)) gam(i) = 0.0_dp
        end do
      ENDIF
      return
    end subroutine vqsatd
!==========================================================================================

    ELEMENTAL REAL(dp) FUNCTION estblf(td)

      real(dp), INTENT(IN) ::  td     ! temperature

      INTRINSIC :: INT
      
     ! tlim is temporary field
     

      estblf  =  (tmin + (tlim(td)-tmin) - tlim(td) + 1.0_dp)           &
                    *estbl(int(tlim(td)-tmin)+1)                        &
                    -(tmin + (tlim(td)-tmin) - tlim(td)      )          &
                    *estbl(int(tlim(td)-tmin)+2)

! rvk      estbl4(td) =  (tmin+int(td-tmin)+1.0-td)*estbl(int(td-tmin)+1)   &
!                      + ( td-(tmin+int(td-tmin)) )*estbl(int(td-tmin)+2)
    END FUNCTION estblf
!-------------------------------------------------------------------------------
    ELEMENTAL REAL(dp) FUNCTION tlim(td)
    
      real(dp), INTENT(IN) :: td         ! temperature

      INTRINSIC :: MAX, MIN

      tlim = max(min(td,tmax),tmin)

    END FUNCTION tlim


!===============================================================================

END MODULE messy_convect_zhang_param
