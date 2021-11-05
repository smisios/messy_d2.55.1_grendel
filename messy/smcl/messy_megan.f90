
MODULE messy_megan

  !-------------------------------------------------------------------------------------------------------
  !  MEGAN : Model of Emissions of Gases and Aerosols from Nature 
  !
  !  AUTHOR:  Pozzer Andrea, MPICH, Feb 2011
  !-------------------------------------------------------------------------------------------------------


  ! MESSy
  USE messy_main_constants_mem,  ONLY: DP, SP, STRLEN_MEDIUM

  IMPLICIT NONE 
  PRIVATE

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'megan'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.0'

  PUBLIC :: dp, sp, STRLEN_MEDIUM

  ! switch for writing change directly in tendency
  LOGICAL, PUBLIC :: l_tendency = .true.   ! GLOBAL SWITCH for tendency calculation

  PUBLIC :: megan_read_nml_ctrl
  PUBLIC :: GAMMA_S
  PUBLIC :: GAMMA_TISOP
  PUBLIC :: GAMMA_TNISP
  PUBLIC :: GAMMA_LAI
  PUBLIC :: GAMMA_P
  PUBLIC :: GAMMA_A

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE gamma_tisop(temp, d_temp, gam_tmp)

!-----------------------------------------------------------------------
!            Calculate GAM_T (GAMMA_T) for isoprene
!-----------------------------------------------------------------------
!                          Eopt*CT2*exp(CT1*x)
!             GAMMA_T =  ------------------------
!                        [CT2-CT1*(1-exp(CT2*x))]
!           where x      = [ (1/Topt)-(1/Thr) ] / 0.00831
!                 Eopt   = 1.75*exp(0.08(Tdaily-297)
!                 CT1    = 80
!                 CT2    = 200
!                 Thr    = hourly average air temperature (K)
!                 Tdaily = daily average air temperature (K)
!                 Topt   = 313 + 0.6(Tdaily-297)
!
!                 Note: AAA = Eopt*CT2*exp(CT1*x)
!                       BBB = [CT2-CT1*(1-exp(CT2*x))]
!                       GAMMA_T = AAA/BBB
!-----------------------------------------------------------------------

    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: gam_tmp(:)
    REAL(dp), INTENT(IN)  :: d_temp(:)
    REAL(dp), INTENT(IN)  :: temp(:)

    REAL, PARAMETER :: CT1 = 80.0
    REAL, PARAMETER :: CT2 = 200.0
    REAL(dp) ::    Eopt, Topt, X
    REAL(dp) ::   AAA, BBB

    INTEGER  :: i


    DO i=1, SIZE(temp)
       IF (TEMP(i).gt.0.0_dp) THEN ! we expect temperature to be always greater than 0 K
          Eopt = 1.75 * exp(0.08*(D_TEMP(i)-297.0))
          Topt = 313.0 + ( 0.6*(D_TEMP(i)-297.0) )
          X = ( (1/Topt)-(1/TEMP(i)) ) / 0.00831
    
          AAA = Eopt*CT2*exp(CT1*X)
          BBB = (  CT2-CT1*( 1-exp(CT2*X) )  )
    
          gam_tmp(i) = AAA/BBB
       ELSE 
          gam_tmp(i) = 0.0_dp
       ENDIF
     ENDDO

  END SUBROUTINE gamma_tisop
!-----------------------------------------------------------------------

  SUBROUTINE gamma_tnisp(TDF_PRM, TEMP,GAM_T )
!-----------------------------------------------------------------------
!              Calculate GAM_T (GAMMA_T) for non-isoprene
!-----------------------------------------------------------------------
!
!             GAMMA_T =  exp[BETA*(T-Ts)]
!           where BETA   = temperature dependent parameter
!                 Ts     = standard temperature (normally 303K, 30C)
!
!     SUBROUITINE GAMMA_TNISP returns the GAMMA_T value for non-isoprene
!-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(dp), INTENT(IN)  :: TEMP(:)
      REAL(dp), INTENT(IN)  :: TDF_PRM
      REAL(dp), INTENT(OUT) :: GAM_T(:)
      REAL(dp), PARAMETER   :: Ts = 303.0

      GAM_T(:) = exp( TDF_PRM*(TEMP(:)-Ts) )

      END SUBROUTINE GAMMA_TNISP

  ! ---------------------------------------------------------------------------

  SUBROUTINE gamma_s(gam_smt)

    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: gam_smt(:)

    gam_smt(:) = 1.0_dp

  END SUBROUTINE gamma_s

  ! ---------------------------------------------------------------------------

  SUBROUTINE GAMMA_LAI(LAI, GAM_L)
!-----------------------------------------------------------------------
!             Calculate GAM_L (GAMMA_LAI)
!-----------------------------------------------------------------------
!                            0.49[LAI]
!             GAMMA_LAI = ----------------    (non-dimension)
!                         (1+0.2LAI^2)^0.5
!
!     SUBROUTINE GAMMA_LAI returns the GAMMA_LAI values
!-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(dp), INTENT(IN)  :: LAI(:)
      REAL(dp), INTENT(OUT) :: GAM_L(:)

      GAM_L(:) = (0.49*LAI(:)) / ( (1+0.2*(LAI(:)**2))**0.5 )

  END SUBROUTINE GAMMA_LAI

  ! ---------------------------------------------------------------------------

  SUBROUTINE GAMMA_P(DAY,cossza, PPFD,D_PPFD, GAM_PHO)

!-----------------------------------------------------------------------
!           Calculate GAM_P (GAMMA_P)
!-----------------------------------------------------------------------
!             GAMMA_P = 0.0         a<=0, a>=180, sin(a) <= 0.0
!
!             GAMMA_P = sin(a)[ 2.46*(1+0.0005(Pdaily-400))*PHI - 0.9*PHI^2 ]
!                                   0<a<180, sin(a) > 0.0
!           where PHI    = above canopy PPFD transmission (non-dimension)
!                 Pdaily = daily average above canopy PPFD (umol/m2s)
!                 a      = solar angle (degree)
!
!                 Note: AAA = 2.46*BBB*PHI - 0.9*PHI^2
!                       BBB = (1+0.0005(Pdaily-400))
!                       GAMMA_P = sin(a)*AAA
!
!                       Pac
!             PHI = -----------
!                   sin(a)*Ptoa
!           where Pac  = above canopy PPFD (umol/m2s)
!                 Ptoa = PPFD at the top of atmosphere (umol/m2s)
!
!             Pac =  SRAD * 4.766 mmmol/m2-s * 0.5
!
!             Ptoa = 3000 + 99*cos[2*3.14-( DOY-10)/365 )]
!           where DOY = day of year
!
!-----------------------------------------------------------------------


      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: DAY          ! day of the year
      REAL(dp), INTENT(IN)  :: cossza(:)    ! cos(zenith angle)
      REAL(dp), INTENT(IN)  :: PPFD(:)      ! 
      REAL(dp), INTENT(IN)  :: D_PPFD(:)
      REAL(dp), INTENT(OUT) :: GAM_PHO(:)

!...  Local parameters
      REAL(dp) :: AAA(size(PPFD))
      REAL(dp) :: BBB(size(PPFD))

      REAL(dp) :: BETA(size(cossza))     ! Solar zenith angle
      REAL(dp) :: SINbeta(size(cossza)) ! sin(beta)
      REAL(dp) :: SZANGLE(size(cossza))! Solar zenith angle array
      
      REAL    Ptoa(size(PPFD)), Pac(size(PPFD)), PHI(size(PPFD))

      REAL(dp), PARAMETER :: PI    = 3.14159265358979323846_dp
      REAL(dp), PARAMETER :: D2RAD = PI/180.0_dp 
      REAL(dp), PARAMETER :: RAD2D = 180.0_dp/PI 

      INTEGER :: i

!...  Begin estimating gamma_p
      ! getting solar radiation
      Pac(:) = PPFD(:)

!...  Initialize parameters
      SZANGLE(:) = 0.
      Ptoa(:) = 0.
      PHI(:) = 0.

      ! Get solar elevation angle
      !BETA(:)=ACOS(cossza(:))*RAD2D         ! DEGREE!!!!!!!
      !TEST
      BETA(:)=ACOS(cossza(:))*RAD2D+90     ! DEGREE: shifted to get correct zenith angle....
      SZANGLE(:) = BETA(:)                 ! Degree
      SINbeta(:) = SIN(BETA(:)*D2RAD)     ! Sin of Zenith angle
!!!!!!!!!!!!!!
!!!ORIGINAL!!!
!      CALL SOLARANGLE( DAY(I,J), HOUR(I,J), LAT(I,J), SINbeta )
!      BETA = ASIN(SINbeta)*RAD2D            ! Degree
!      SZANGLE(I,J) =  BETA                  ! Degree
!!!!!!!!!!!!!!


      DO i=1,size(PPFD)
         IF (SINbeta(i) .LE. 0.0) THEN
            GAM_PHO(i) = 0.0_dp
         ELSEIF (SINbeta(i) .GT. 0.0) THEN
            Ptoa(i) = 3000.0 + 99.0 * COS( 2*3.14*(DAY-10)/365 )

            PHI(i) = Pac(i)/(SINbeta(i) * Ptoa(i))

            BBB(i) = 1 + 0.0005*(D_PPFD(i)-400)
            AAA(i) = ( 2.46 * BBB(i) * PHI(i) ) - ( 0.9 * PHI(i)**2 )

            GAM_PHO(i) = SINbeta(i) * AAA(i)
         ELSE
            write(*,*) 'Error: Solar angle is invalid'

         ENDIF

         ! Screening the unforced errors
         ! IF solar elevation angle is less than 1 THEN
         ! gamma_p can not be greater than 0.1.
         IF (BETA(i) .LT. 1.0 .AND. GAM_PHO(i) .GT. 0.1) THEN
            GAM_PHO(i) = 0.0
         ENDIF

         ! additional test
         IF (GAM_PHO(i) .LT. 0.0) THEN
            GAM_PHO(i) = 0.0
         ENDIF
      ENDDO


  END SUBROUTINE GAMMA_P

  ! ---------------------------------------------------------------------------

  SUBROUTINE GAMMA_A(SPC_NAME, LAIp, LAIc, TSTLEN, D_TEMP, GAM_AGE)

!
!!-----------------------------------------------------------------------
!!          Calculate GAM_A (GAMMA_age)
!!-----------------------------------------------------------------------
!!
!!             GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
!!           where Fnew = new foliage fraction
!!                 Fgro = growing foliage fraction
!!                 Fmat = mature foliage fraction
!!                 Fold = old foliage fraction
!!                 Anew = relative emission activity for new foliage
!!                 Agro = relative emission activity for growing foliage
!!                 Amat = relative emission activity for mature foliage
!!                 Aold = relative emission activity for old foliage
!!
!!
!!             For foliage fraction
!!             Case 1) LAIc = LAIp
!!             Fnew = 0.0  , Fgro = 0.1  , Fmat = 0.8  , Fold = 0.1
!!
!!             Case 2) LAIp > LAIc
!!             Fnew = 0.0  , Fgro = 0.0
!!             Fmat = 1-Fold
!!             Fold = (LAIp-LAIc)/LAIp
!!
!!             Case 3) LAIp < LAIc
!!             Fnew = 1-(LAIp/LAIc)                       t <= ti
!!                  = (ti/t) * ( 1-(LAIp/LAIc) )          t >  ti
!!
!!             Fmat = LAIp/LAIc                           t <= tm
!!                  = (LAIp/LAIc) +
!!                      ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )  t >  tm
!!
!!             Fgro = 1 - Fnew - Fmat
!!             Fold = 0.0
!!
!!           where
!!             ti = 5 + (0.7*(300-Tt))                   Tt <= 303
!!                = 2.9                                  Tt >  303
!!             tm = 2.3*ti
!!
!!             t  = length of the time step (days)
!!             ti = number of days between budbreak and the induction of
!!                  emission
!!             tm = number of days between budbreak and the initiation of
!!                  peak emissions rates
!!             Tt = average temperature (K) near top of the canopy during
!!                  current time period (daily ave temp for this case)
!!
!!
!!             For relative emission activity
!!             Case 1) Constant
!!             Anew = 1.0  , Agro = 1.0  , Amat = 1.0  , Aold = 1.0
!!
!!             Case 2) Monoterpenes
!!             Anew = 2.0  , Agro = 1.8  , Amat = 0.95 , Aold = 1.0
!!
!!             Case 3) Sesquiterpenes
!!             Anew = 0.4  , Agro = 0.6  , Amat = 1.075, Aold = 1.0
!!
!!             Case 4) Methanol
!!             Anew = 3.0  , Agro = 2.6  , Amat = 0.85 , Aold = 1.0
!!
!!             Case 5) Isoprene
!!             Anew = 0.05 , Agro = 0.6  , Amat = 1.125, Aold = 1.0
!!
!!-----------------------------------------------------------------------
      IMPLICIT NONE



      CHARACTER(LEN=*), INTENT(IN) :: SPC_NAME
      REAL(dp), INTENT(IN)  :: D_TEMP(:)
      REAL(dp), INTENT(IN)  :: LAIp(:)
      REAL(dp), INTENT(IN)  :: LAIc(:)
      INTEGER,  INTENT(IN)  :: TSTLEN
      REAL(dp), INTENT(OUT) :: GAM_AGE(:)
!
!     Local parameters
      REAL(dp) ::  Fnew, Fgro, Fmat, Fold
      ! TODO : could be moved to namelist
      REAL(dp) ::  Anew, Agro, Amat, Aold

      INTEGER  :: t                 ! time step
      REAL(dp) :: ti                ! number of days between budbreak
                                    ! and the induction of emission
      REAL(dp) :: tm                ! number of days between budbreak
                                    ! and the initiation of peak
                                    ! emissions rates
      REAL(dp) :: Tt                   ! average temperature (K)
                                    ! daily ave temp

      INTEGER   :: I
!

!!...  Choose relative emission activity
      SELECT CASE (TRIM(SPC_NAME))
      CASE ('ACTO','ACTA','FORM','CH4','NO','CO')
             Anew = 1.0  
             Agro = 1.0  
             Amat = 1.0  
             Aold = 1.0
      CASE ('MYRC','SABI','LIMO','CAR3','OCIM','BPIN','APIN','OMTP')
             Anew = 2.0  
             Agro = 1.8  
             Amat = 0.95 
             Aold = 1.0
      CASE ('FARN','BCAR','OSQT')
             Anew = 0.4  
             Agro = 0.6  
             Amat = 1.075
             Aold = 1.0
      CASE ('MEOH')
             Anew = 3.0  
             Agro = 2.6  
             Amat = 0.85 
             Aold = 1.0
      CASE ( 'ISOP','MBO' )
             Anew = 0.05 
             Agro = 0.6  
             Amat = 1.125
             Aold = 1.0
      CASE DEFAULT
         WRITE(*,*) 'Error: Chemical species, invalid variable: '//TRIM(SPC_NAME)
         RETURN
      ENDSELECT
!
      t = TSTLEN
      DO i = 1, SIZE(LAIc)
            Tt   = D_TEMP(i)
!...  Calculate foliage fraction
            IF (LAIp(i) .EQ. LAIc(i)) THEN
               Fnew = 0.0
               Fgro = 0.1
               Fmat = 0.8
               Fold = 0.1
            ELSEIF (LAIp(i) .GT. LAIc(i)) THEN
               Fnew = 0.0
               Fgro = 0.0
               Fold = ( LAIp(i)-LAIc(i) ) / LAIp(i)
               Fmat = 1-Fold

            ELSEIF (LAIp(i) .LT. LAIc(i)) THEN
!              Calculate ti and tm
               IF (Tt .LE. 303.0) THEN
                  ti = 5.0 + 0.7*(300-Tt)
               ELSEIF (Tt .GT. 303.0) THEN
                  ti = 2.9
               ENDIF
               tm = 2.3*ti

!              Calculate Fnew and Fmat, then Fgro and Fold
!              Fnew
               IF (t .LE. ti) THEN
                  Fnew = 1.0 - (LAIp(i)/LAIc(i))
               ELSEIF (t .GT. ti) THEN
                  Fnew = (ti/t) * ( 1-(LAIp(i)/LAIc(i)) )
               ENDIF

!              Fmat
               IF (t .LE. tm) THEN
                  Fmat = LAIp(i)/LAIc(i)
               ELSEIF (t .GT. tm) THEN
                  Fmat = (LAIp(i)/LAIc(i)) + ( (t-tm)/t ) * ( 1-(LAIp(i)/LAIc(i)) )
               ENDIF

               Fgro = 1.0 - Fnew - Fmat
               Fold = 0.0
         
            ENDIF

!...  Calculate GAMMA_A
      GAM_AGE(i) = Fnew*Anew + Fgro*Agro +  &
                   Fmat*Amat + Fold*Aold
!
      ENDDO       ! End loop for SIZE(Laic) 
!
  END SUBROUTINE GAMMA_A

  ! ---------------------------------------------------------------------------

  SUBROUTINE megan_read_nml_ctrl(status, iou)

    !  megan MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Pozzer Andrea, MPICH, Oct 2004


    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CTRL/ l_tendency      

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='megan_read_nml_ctrl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    l_tendency = .false.

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR

  END SUBROUTINE megan_read_nml_ctrl


END MODULE
