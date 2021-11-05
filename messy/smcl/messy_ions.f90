! Ionisation in the atmosphere
! This module provides Ion Pair production rates and recombination rates for other submodels. 
!
! List of functions and subroutines:
!
! - read_look_up, (subroutine) :  TODO: Replace by MESSy import
!   reads the look up tables for the Usoskin parameterisation
!
! - recom_brasseur,  (elemental function) :
!   recombination rate constant based on the parameterisation of Brasseur
!
! - steady_state_ions, (elemental function) :
!   steady state ion concentration, no charge inbalance
!
! - USOSKIN_CRII, (pure function) :
!   calculate the ion pair production rate from GCR based on Usoskin
!
!  - decay_ipr, (pure function) :
!    calculate the ion pair production ratee from radioactive decay
!
!  - PARAM_NOW, (pure subroutine) :
!    determine parameters for solar modulation and geomagnetic field
!
!  - GEOMAG, (pure function) :
!    CALCULATES THE GEOMAGNETIC PARAMETERS for GEOMAGNETIC CUTOFF RIGIDITY AT EACH LATITUDE AND LONGITUDE
!
!  - CALC_PC: calculates geomagnetic cut off rigidity
!
!
! author: sebastian ehrhart
! date: Nov 2016
! last modified: 03 Mai 2018

module messy_ions
  
  use messy_main_constants_mem, only: wp=>dp, dp, strlen_ulong, pi

  implicit none

  private

  public :: read_ionisation_nml
  public :: recom_brasseur
  public :: ion_aerosol_TZ06
  public :: steady_state_ions
  public :: USOSKIN_CRII
  public :: decay_ipr
  public :: PARAM_NOW ! determine parameters for solar modulation and geomagnetic field 
  public :: GEOMAG ! CALCULATES THE GEOMAGNETIC PARAMETERS for GEOMAGNETIC CUTOFF RIGIDITY AT EACH LATITUDE AND LONGITUDE
  public :: CALC_PC

  character(LEN=*), parameter, public :: modstr = 'ions'
  character(LEN=*), parameter, public :: modver = '1.0'

  ! control namelists
  integer, public, save :: gcr_method=0     ! which method to determine GCR contribution
  logical, public, save :: lqradon=.FALSE.   ! include ions from radon decay
  logical, public, save :: lspe=.FALSE.     ! include ions from solar proton events
  logical, public, save :: lgcr=.FALSE.     ! ionisation by galactic cosmic rays? will be set automtically by gcr_method
  logical, public, save :: laero=.FALSE.    ! aerosol ion interaction
  character(len=strlen_ulong), save :: path2lut='' ! absolute (directories + filename) path to the lookup table file

  ! parameters used for calculation of ion production from decay
  integer, parameter, public :: n_Rndec=5 ! number of Radon decay tracers
  real(kind=wp), parameter :: en_ip_alpha=35.6_wp ! eV consumed for the production of an ion pair from alpha particles
  real(kind=wp), parameter :: en_ip_beta =32.5_wp ! eV consumed for the production of an ion pair from beta particles
  real(kind=wp), parameter, public :: Rn_chain_ions(n_Rndec) = (/5.59e6_wp/en_ip_alpha,& ! 222Rn
                                                           6.12e6_wp/en_ip_alpha,& ! 218 Po
                                                           1.02e6_wp/en_ip_beta, & ! 214 Pb
                                                           3.27e6_wp/en_ip_beta, & ! 214Po
                                                           7.88_wp/en_ip_alpha/)   !210Pb
  ! pre-calculted for chain integration
  real(wp), parameter, public :: lam(n_Rndec) = (/ &  ! rate constant pof radioactive decay (s-1)
        2.111193897904317E-6_wp &  ! 222Rn -> 218Po
      , 3.850817669777474E-3_wp &  ! 218Po -> 214Pb
      , 4.278686299752749E-4_wp &  ! 214Pb -> 214Bi
      , 5.776226504666211E-4_wp &  ! 214Bi -> 214Po -> 210Pb
      , 9.990705868040354E-10_wp /) ! 210Pb -> ... -> 206Pb


 

  ! variables and parameters for the Usoskin parameterisation
#ifndef LF
  real(kind=wp), parameter :: PPI = real(pi, kind=wp) ! convert to precision used here
#else
  real(kind=wp), parameter :: PPI = pi ! convert to precision used here
#endif
!!$  real(kind=wp), public, save :: diffcrii(131,31,77) ! op_pj_20180530 not needed; waste of memory
  ! look up table import
  real(kind=wp), public, pointer, dimension(:,:,:), save :: CRII_TABLES => NULL()
  real(kind=wp), public, pointer, dimension(:), save :: H      => NULL()
  real(kind=wp), public, pointer, dimension(:), save :: PHI_AV => NULL()
  real(kind=wp), public, pointer, dimension(:), save :: PC     => NULL()

contains

!---- 1. initilise

subroutine open_ascii(path, iou, modstr)

  implicit none
  character(len=*), intent(in) :: path
  integer, intent(in) :: iou
  character(len=*), intent(in) :: modstr
  ! local
  logical :: lex ! file exists


  ! check if file exists
  inquire(file=trim(path), exist=lex)
  if (.not.lex) then
     write(*,*) '*** WARNING: FILE '''// trim(path) //'''  NOT FOUND !'
     return
  end if

  ! OPEN FILE
  open(iou,file=trim(path))
  write(*,*) 'Reading LOOK UP TABLE '''//trim(path)//''''//' for '''//trim(modstr)// ''' (unit ',iou,') ...'

end subroutine open_ascii


subroutine close_ascii(path, iou, modstr)

  implicit none
  character(len=*), intent(in) :: path
  integer, intent(in) :: iou
  character, intent(in) :: modstr

  close(iou)

  write(*,*) 'CLOSING LOOK UP TABLE file' //trim(path) // 'for' // trim(modstr)

end subroutine close_ascii



subroutine read_ionisation_nml(status, iou)

  use messy_main_tools, only : read_nml_open, read_nml_check, read_nml_close

  implicit none
  integer, intent(out) :: status ! error status: 0 = no error, 1 = an error
  integer, intent(in)  :: iou    ! I/O unit
  ! local
  character(len=*), parameter :: substr='read_ionisation_nml'
  logical :: lex
  integer :: fstat

  NAMELIST /CTRL/ gcr_method, lqradon, lspe, path2lut, laero

  status=1

  call read_nml_open(lex, substr, iou, 'CTRL', modstr)
  if(.not.lex) return ! can't open the file


  !---------- read in parameterisation
  read(iou, nml=CTRL, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'CTRL', modstr)
  if(fstat /= 0) return ! error while reading the namelist

  if(gcr_method > 0) lgcr=.TRUE.

  call read_nml_close(substr, iou, modstr)
  status = 0 ! all fine
 
end subroutine read_ionisation_nml


!---- 2. Processe producing and destroying ions

elemental function recom_brasseur(t,m) result(y)
  ! recombination rate parameterisation given in Brasseur & Chatel 1983
  implicit none
  ! input
  real(kind=wp), intent(in) :: t ! temperature in K
  real(kind=wp), intent(in) :: m ! total concentration of air in cm-3
  ! local
  real(kind=wp) :: isct ! inverse scaled temperature, i.e. 300/t
  ! return
  real(kind=wp) :: y ! cm3 s-1

  isct = 300.0_wp/t
  y = 6.0e-8_wp*sqrt(isct) + 6.0e-26_wp*m*isct*isct*isct*isct

end function recom_brasseur


pure function ion_aerosol_TZ06(r) result(y)
  ! sources:  Tinsley & Zhou 2006 DOI: 10.1029/2005JD006988
  !           Baumgaertner et al 2013 DOI: 10.1002/jgrd.50725
  implicit none
  real(kind=wp), intent(in) :: r ! radius of aeroso particle in m
  ! return value
  real(kind=wp) :: y

  y = 0.0_wp 

  if(r > 0.01_wp) then
  ! if radius of aerosol particle > 0.01 µm ( we use cm)
  ! note that the equation in Baumgaertner et al 2013 uses µm
    y = 4.36e-5_wp*r - 9.2e-8_wp ! cm3 s-1
  elseif( r < 1e-15_wp) then 
    ! protect against zero radius
    y = 0.0_wp
  else
    y = 10.0_wp**(1.243_wp*log10(r) - 3.978_wp)
  end if 

end function 



pure function decay_ipr(x,ions) result(q)
  ! function that calculates ion pair production from decay of radiocative material
  !
  implicit none 
  real(kind=wp), intent(in) :: x(:)    ! decay rate (deays cm-3 s-1)
  real(kind=wp), intent(in) :: ions(:) ! ion pairs produced per decay
  ! return
  real(kind=wp) :: q ! total ion pair production from decay

  q=sum(x*ions)

end function decay_ipr


elemental function steady_state_ions(q, kr, kl) result(n)
  ! steady state solution of 
  ! dN / dt = q - kr*N**2 - kl*N 
  implicit none
  real(kind=wp), intent(in) :: q   ! ion pair production rate (ion pairs cm-3 s-1)
  real(kind=wp), intent(in) :: kr  ! ion-ion recombination ( cm3 s-1)
  real(kind=wp), intent(in), optional :: kl  ! first order loss term (s-1)
  ! return
  real(kind=wp) :: n ! steady state ion concentration (cm-3)

  if(present(kl)) then
    n = (sqrt(kl*kl + 4.0_wp*kr*q) - kl) / (kr+kr)
  else
    n = sqrt(q/kr)
  end if

end function steady_state_ions 


pure function steady_state_ions_assym(q, kr, klp, kln) result(n)
  ! calculates the steady state ion pair production rate in case of assymetric first order losse
  ! dNp = q - kr*Np*Nn - klp*Np
  ! dNn = q - kr*Np*Nn - kln*Nn
  implicit none
  real(kind=wp), intent(in) :: q(:)   ! ion pair production rate (ion pairs cm-3 s-1)
  real(kind=wp), intent(in) :: kr(:)  ! ion-ion recombination ( cm3 s-1)
  real(kind=wp), intent(in) :: klp(:) ! first order loss term positive ions (s-1)
  real(kind=wp), intent(in) :: kln(:) ! first order loss term negative ions (s-1)
  ! locale
  real(kind=wp) :: klpokln(size(q)) ! klp/kln 
  ! return
  real(kind=wp) :: n(size(q),2) ! steady state ion concentration (cm-3)

  klpokln = klp/kln

  n(:,1) = steady_state_ions(q, kr/klpokln, kl=klp)
  n(:,2) = steady_state_ions(q, kr*klpokln, kl=kln)

end function steady_state_ions_assym



! Author of the original version of the code for GLOMAP was Eimear Dunne
! version below modified by: Sebastian Ehrhart
FUNCTION USOSKIN_CRII(HGT, PHI_NOW, P_C) result(CRII_NEW)
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !
  ! Purpose
  ! -------
  !  THIS FUNCTION IS CALLED AT THE START OF EACH MONTH TO CALCULATE THE
  !  IONIZATION RATE DUE TO GALACTIC COSMIC RAYS AT THE END OF THE MONTH.
  !  THE CRII IS THEN INTERPOLATED LINEARLY BETWEEN THE START AND END OF
  !  THE MONTH.
  !  EDIT: Change the function call. It calculates the ion pair production
  !  rate for a given month and year. Interpolation is done outside of this
  !  routine in the submodel interface layer. 
  !
  ! References
  ! ----------
  ! I. G. Usoskin, G. A. Kovaltsov, I. A. Mironova, Cosmic ray induced ionization model CRAC:CRII:
  ! An extension to the upper atmosphere, Journal of Geophysical Research 115, 6 (2010).
  ! A. !. Fraser-Smith, Centered and eccentric geomagnetic dipoles and their poles, 16001985, Re-
  ! views of Geophysics 25, 1–16 (1987).
  ! Coding by E.M. Dunne in GLOMAP model, 2010-2011. 
  !
  ! Parameters
  ! ----------
  !
  ! Inputs
  ! ------
  ! HGT     : atmospheric depth
  ! PHI_NOW : PHI_NOW
  ! P_C     : geomagnetic cut off rigidity
  !
  ! Outputs
  ! -------
  !  CRII  : CRII calculated for the start of the current month
  !
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   
        
    IMPLICIT NONE
    ! input
    real(kind=wp), intent(in) :: HGT
    real(kind=wp), intent(in) :: PHI_NOW
    real(kind=wp), intent(in) :: P_C 
  
    ! output
    real(kind=wp) :: CRII_NEW 
  
    CRII_NEW = USOSKIN_LOOKUP(PHI_NOW, HGT, P_C)
    if(CRII_NEW < 0.0) CRII_NEW = 0.0
  
END FUNCTION USOSKIN_CRII


FUNCTION USOSKIN_LOOKUP(PHI_NOW, HGT, P_C) result(LOOKUP)!, CRII_TABLES)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Purpose
! -------
! WORKS OUT THE CRII USING TRI-LINEAR INTERPOLATION OF PHI, H, P_C
!
! Parameters
! ----------
!      
! Inputs
! ------
! PHI_NOW       : CURRENT VALUE OF PHI (OULU COSMIC RAY STATION)
! HGT           : ATMOSPHERIC DEPTH (FROM PRESSURE AT MID-LEVEL)       
! P_C           : GEOMAGNETIC CUTOFF RIGIDITY (CALCULATED FROM IGRF)
! CRII_TABLES   : LOOKUP TABLES OF CRII FOR H, PHI, P_C
!
! Outputs
! -------
! LOOKUP        : CURRENT/LOCAL VALUE OF CRII IN GIVEN BOX/MONTH
!
! Local Variables
! ---------------
! Cxxx          : 8 VALUES OF CRII FORMING A CUBE IN PHASE SPACE
!                 CONTAINING THE LOCAL VALUE
! Cxx           : 4 VALUES OF CRII, INTERPOLATED IN P_C FROM Cxxx
! Cx            : 2 VALUES OF CRII, INTERPOLATED IN PHI FROM Cxx
! LOOKUP        : VALUE OF CRII, INTERPOLATED IN H FROM Cx
! PC            : VALUES OF PC FOR WHICH CRII ARE AVAILABLE
! PHI_AV        : VALUES OF PHI FOR WHICH CRII ARE AVAILABLE
! H             : VALUES OF H FOR WHICH CRII ARE AVAILABLE
! PHI_UP/DOWN   : CLOSEST VALUES OF PHI FOR WHICH CRII AVAILABLE
! PC_UP/DOWN    : CLOSEST VALUES OF P_C FOR WHICH CRII AVAILABLE
! H_UP/DOWN     : CLOSEST VALUES OF H FOR WHICH CRII AVAILABLE
! PHI/PC/H_INT  : PREFACTOR USED IN INTERPOLATION      
!
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    IMPLICIT NONE
  
    real(kind=wp), intent(in) :: PHI_NOW, P_C, HGT 
    integer :: I
    integer :: PHI_UP, PHI_DOWN, PC_UP, PC_DOWN, H_UP, H_DOWN 
    real(kind=wp) :: C0, C1, C00, C01, C10, C11
    real(kind=wp) :: C000, C001, C010, C011, C100, C101, C110, C111
    real(kind=wp) :: PHI_INT, PC_INT, H_INT
    real(kind=wp) :: LOOKUP

  ! PUT IN A SET OF WARNINGS THAT IF PHI_AV, P_C, H FALL OUTSIDE THESE
  ! REGIONS (EXCEPT FOR H > 1025) THE PROGRAM SHOULD QUIT (TRY H > 1035)
!        IF((P_C .LT. 0.1) .OR. (P_C .GT. 39.0)) THEN
!          PRINT *, 'ERROR: P_C NOT WITHIN RANGE [0.1, 39.0]'
!          PRINT *, 'FAULTY INTERPOLATION MAY RESULT' 
!        ENDIF 
!        IF((PHI_NOW .LT. 0.0) .OR. (PHI_NOW .GT. 1500.0)) THEN
!          PRINT *, 'ERROR: PHI_AV NOT WITHIN RANGE [0.0, 1500.0]'
!          PRINT *, 'FAULTY INTERPOLATION MAY RESULT' 
!        ENDIF 
!        IF((HGT .LT. 10.0) .OR. (HGT .GT. 1035.0)) THEN
!          PRINT *, 'ERROR: H NOT WITHIN RANGE [10, 1035.0]'
!          PRINT *, 'FAULTY INTERPOLATION MAY RESULT' 
!        ENDIF 
  
  ! FOR EACH OF PHI_AV, P_C, H: SET INDEX OF UPPER AND LOWER BOUNDS
  ! WORK UP TO INDICES DIRECTLY BELOW CURRENT VALUES
  ! WORK DOWN TO INDICES DIRECTLY ABOVE CURRENT VALUES
    PHI_UP = size(PHI_AV) !#
    PHI_DOWN = 1 
    PC_UP = size(PC) !#
    PC_DOWN = 1 
    H_UP = size(H) !#
    H_DOWN = 2
  
    DO I = 1, size(PHI_AV) !#
       IF ( PHI_AV(I) .LE. PHI_NOW ) THEN
          PHI_DOWN = I 
       END IF
       IF( PHI_AV(size(PHI_AV)+1 - I) .GE. PHI_NOW ) THEN !#
          PHI_UP = size(PHI_AV) + 1 - I !#
      ! IF ( PHI_AV(32 - I) .GE. PHI_NOW ) THEN 
      !    PHI_UP = 32 - I ! size(PHI_AV) + 1 - I !#
       END IF
    END DO
  
    DO I = 1, size(PC)
       IF ( PC(I) .LE. P_C ) THEN
          PC_DOWN = I
       END IF
       IF ( PC(size(PC) + 1 - I) .GE. P_C ) THEN  !#
          PC_UP = size(PC) + 1 - I                !#
       END IF                           !#
!       IF ( PC(82 - I) .GE. P_C ) THEN
!          PC_UP = 82 - I 
!       END IF
    END DO
  
    DO I = 1, size(H)
       IF ( H(I) .LE. HGT ) THEN
          H_DOWN = I
       END IF
       IF ( H(size(H)+1-I) .GE. HGT ) THEN !#
          H_UP = size(H)+1-I               !#
       END IF                              !#
!       IF ( H(61-I) .GE. HGT ) THEN
!          H_UP = 61 - I
!       END IF
    END DO
  
    IF( H_DOWN .EQ. 1) H_DOWN = 2
    IF( H_UP .EQ. 1) H_UP = 2
  
  ! VALUES OF UP/DOWN SHOULD ONLY BE ONE APART  
    IF(PHI_UP - PHI_DOWN .NE. 1) THEN
      IF(PHI_UP .EQ. PHI_DOWN) THEN
!       PRINT *, 'VALUE OF PHI_AV OUTSIDE RANGE'
        IF(PHI_UP .EQ. 1) THEN
!         PRINT *, 'PHI TOO SMALL' 
          PHI_UP = PHI_DOWN + 1
        ELSE IF(PHI_DOWN .EQ. 31) THEN
!         PRINT *, 'PHI TOO LARGE' 
          PHI_DOWN = PHI_UP - 1
        END IF
      ELSE
!        PRINT *, 'VALUE OF PHI IS ON BOUNDARY.' 
      END IF
    ENDIF
  
    IF(PC_UP - PC_DOWN .NE. 1) THEN
      IF(PC_UP .EQ. PC_DOWN) THEN
!       PRINT *, 'VALUE OF P_C OUTSIDE RANGE'
        IF(PC_UP .EQ. 1) THEN
!         PRINT *, 'P_C TOO SMALL', P_C
          PC_UP = PC_DOWN + 1
        ELSE IF(PC_DOWN .EQ. 81) THEN
!         PRINT *, 'P_C TOO LARGE', P_C 
          PC_DOWN = PC_UP - 1
        END IF
      ELSE
!        PRINT *, 'VALUE OF P_C IS ON BOUNDARY.' 
      END IF
    ENDIF
  
    IF(H_UP - H_DOWN .NE. 1) THEN
      IF(H_UP .EQ. H_DOWN) THEN
!       PRINT *, 'VALUES OF H OUTSIDE RANGE'
        IF(H_UP .EQ. 2) THEN
!         PRINT *, 'H TOO SMALL' 
          H_UP = H_DOWN + 1
        ELSE IF(H_DOWN .EQ. 60) THEN
!         PRINT *, 'H TOO LARGE' 
          H_DOWN = H_UP - 1
        END IF
      ELSE
!        PRINT *, 'VALUE OF H IS ON BOUNDARY.' 
      END IF
    ENDIF
  
  ! NOW INTERPOLATE
  
    C000 = CRII_TABLES(H_DOWN, PHI_DOWN, PC_DOWN) 
    C001 = CRII_TABLES(H_DOWN, PHI_DOWN, PC_UP) 
    C010 = CRII_TABLES(H_DOWN, PHI_UP, PC_DOWN) 
    C011 = CRII_TABLES(H_DOWN, PHI_UP, PC_UP) 
    C100 = CRII_TABLES(H_UP, PHI_DOWN, PC_DOWN) 
    C101 = CRII_TABLES(H_UP, PHI_DOWN, PC_UP) 
    C110 = CRII_TABLES(H_UP, PHI_UP, PC_DOWN) 
    C111 = CRII_TABLES(H_UP, PHI_UP, PC_UP)
    
    IF(PHI_UP .NE. PHI_DOWN) THEN
    PHI_INT = (PHI_NOW - PHI_AV(PHI_DOWN))/(PHI_AV(PHI_UP) - PHI_AV(PHI_DOWN))
    ELSE 
    PHI_INT = 0.0
    END IF
    IF(PC_UP .NE. PC_DOWN) THEN
    PC_INT = (P_C - PC(PC_DOWN)) / (PC(PC_UP) - PC(PC_DOWN)) 
    ELSE
    PC_INT = 0.0
    END IF
    IF(H_UP .NE. H_DOWN) THEN
    H_INT = (HGT - H(H_DOWN)) / (H(H_UP) - H(H_DOWN)) 
    ELSE
    H_INT = 0.0
    END IF
  
    C00 = C000 + PC_INT * (C001 - C000) 
    C01 = C010 + PC_INT * (C011 - C010)  
    C10 = C100 + PC_INT * (C101 - C100) 
    C11 = C110 + PC_INT * (C111 - C110) 
  
    C0 = C00 + PHI_INT * (C01 - C00) 
    C1 = C10 + PHI_INT * (C11 - C01) 
  
    LOOKUP = C0 + H_INT * (C1 - C0)
 
END FUNCTION USOSKIN_LOOKUP
      


pure subroutine PARAM_NOW(MOIS, IAN, PHI, IGRF, PHI_NOW, GH_NOW)

  implicit none
  integer, intent(in) :: MOIS, IAN
  real(kind=wp), intent(in) :: PHI(:,:), IGRF(:,:)
  real(kind=wp), intent(out) :: GH_NOW(:), PHI_NOW
  ! local
  integer :: DGRF
  real(kind=wp) :: GH(2,size(GH_NOW))

  PHI_NOW = PHI(IAN-1951,MOIS+1)
!  write(*,*) "PHI_NOW=", phi_now

  DGRF = (IAN - 1900) / 5 + 1
!  write(*,*) "DGRF=", dgrf

  GH = IGRF(DGRF:DGRF+1,2:9)
  ! LINEARLY INTERPOLATE GH TO CURRENT VALUE GH_NOW
  CALL GH_INT(GH, GH_NOW, MOIS, IAN)

end subroutine PARAM_NOW


PURE SUBROUTINE GEOMAG(GH, B_0, D, L, M, X)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Purpose
! -------
! CALCULATES THE GEOMAGNETIC PARAMETERS WHICH ARE USED TO CALCULATE THE
! GEOMAGNETIC CUTOFF RIGIDITY AT EACH VALUE OF LATITUDE AND LONGITUDE
!
! Parameters
! ----------
!
! Inputs
! ------
! GH            : CURRENT VALUES OF IGRF COEFFICIENTS
!
! Outputs
! -------
! B_0           : Earth's B field      
! D             : Length of vector from Earth's centre to centre of
!                 eccentric dipole
! L             : Functions of IGRF coefficients used to calculate P_c
! M             : Earth's M field (in units of 10^22 A m^2)
! X             : Vector from Earth's centre to centre of eccentric dipole 
!
! Local Variables
! ---------------
! 
!
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  IMPLICIT NONE

  real(kind=wp), intent(in)  :: GH(8)
  real(kind=wp), intent(out) :: B_0, D, L(3), M, X(3)
  real(kind=wp) :: E, R_0, MU_0
  
  PARAMETER(R_0 = 6371.0)    ! see if it is available from MESSy
  PARAMETER(MU_0=1.25663706E-6) ! see if it is available from MESSy 


  B_0 = SQRT(GH(1)*GH(1) + GH(2)*GH(2) + GH(3)*GH(3))
  M = 4 * PPI * B_0 * (R_0 ** 3) / (MU_0 * 10.0**22)

  L(1) = 2 * GH(1)*GH(4) + SQRT(3.0)*(GH(2) * GH(5) + GH(3) * GH(6))
  L(2) =-GH(2)*GH(4)+SQRT(3.0)*(GH(1)*GH(5)+GH(2)*GH(7)+GH(3)*GH(8))
  L(3) =-GH(3)*GH(4)+SQRT(3.0)*(GH(1)*GH(6)-GH(3)*GH(7)+GH(2)*GH(8))

  E = (L(1) * GH(1) + L(2) * GH(2) + L(3)*GH(3)) / (4 * B_0 * B_0)

  X(1) = R_0 * (L(2) - GH(2) * E) / (3 * B_0 * B_0)
  X(2) = R_0 * (L(3) - GH(3) * E) / (3 * B_0 * B_0)
  X(3) = R_0 * (L(1) - GH(1) * E) / (3 * B_0 * B_0) 

  D = SQRT(X(1) * X(1) + X(2) * X(2) + X(3) * X(3))

END SUBROUTINE GEOMAG


PURE SUBROUTINE GH_INT(GH, GH_NOW, MOIS, IAN)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Purpose
! -------
! INTERPOLATES BETWEEN THE START AND END OF EACH FIVE-YEAR BRACKET
! PROVIDED BY THE IGRF TO GIVE THE CURRENT VALUES OF THE IGRF
! COEFFICIENTS
!
! Parameters
! ----------
!
! Inputs
! ------
! GH            : VALUES OF COEFFICIENTS AT START AND END OF BRACKET
!                 WHICH CONTAINS THE CURRENT DATE
! MOIS          : MONTH IN INTEGER FORMAT
! IAN           : YEAR IN INTEGER FORMAT      
!
! Outputs
! -------
! GH_NOW        : CURRENT VALUE OF IGRF COEFFICIENTS BASED ON LINEAR
!                 INTERPOLATION
!
! Local Variables
! ---------------
! NOW           : CURRENT DATE'S LOCATION IN FIVE-YEAR BRACKET
!
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  IMPLICIT NONE 

  real(kind=wp), intent(in) :: GH(2,8)
  integer,       intent(in) :: MOIS, IAN
  real(kind=wp), intent(out) :: GH_NOW(8) 
  real(kind=wp) :: NOW

  NOW = REAL(MOD(IAN,5)) + REAL(MOIS-1)/12.0 

  GH_NOW = GH(1,:) + NOW * (GH(2,:) - GH(1,:)) / 5.0

END SUBROUTINE GH_INT


PURE FUNCTION CALC_PC(B_0, D, GH, M, X, THETA, PSI)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Purpose
! -------
! CALCULATES THE GEOMAGNETIC CUTOFF RIGIDITY USING IGRF-11 AND THE
! METHOD OF A.!. FRASER-SMITH (1987) "CENTERED AND ECCENTRIC GEOMAGNETIC
! DIPOLES AND THEIR POLES"
!
! Parameters
! ----------
!
! Inputs
! ------
! B_0           : Earth's B field
! D             : Length of vector from Earth's centre to centre of
!                 eccentric dipole
! GH            : Current value of IGRF coefficients
! L             : Functions of IGRF coefficients used to calculate P_c
! M             : Earth's M field (in units of 10^22 A m^2)
! X             : Vector from Earth's centre to centre of eccentric dipole 
! THETA         : Geographical co-latitude at centre of grid box
! PSI           : Geographical longitude at centre of grid box      
!
! Outputs
! -------
! CALC_PC       : Current local value of geomagnetic cutoff rigidity 
!
! Local Variables
! ---------------
! S/C_T/P       : Sin/Cos of Theta/Psi
! R             : Magnitude of vector from centre of eccentric dipole to point of
!                 to point of observation      
! THETA_G       : Angle between magnetic dipole axis and vector R 
! S_TG          : Sin of Theta_G      
!
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  IMPLICIT NONE

  real(kind=wp), intent(in) :: B_0, D, GH(8), M, X(3), THETA, PSI
  real(kind=wp) :: CALC_PC
  real(kind=wp) :: S_T, C_T, S_P, C_P, THETA_G, S_TG
  real(kind=wp) :: R_0, R

!    PARAMETER(PPI=3.14159265358979323846)
  PARAMETER(R_0 = 6371.0)    

  S_T = SIN(THETA)
  C_T = COS(THETA) 
  S_P = SIN(PSI)
  C_P = COS(PSI)

  R = SQRT(R_0 * R_0 + D * D - 2 * R_0 * (X(1) * S_T * C_P & 
 & + X(2) * S_T * S_P + X(3) * C_T))

  THETA_G = ACOS((GH(2) * (X(1) - R_0 * S_T * C_P) + GH(3) * &
 & (X(2) - R_0 * S_T * S_P) + GH(1) * (X(3) - R_0 * C_T) ) / &
 & (R * B_0))

  S_TG = SIN(THETA_G)

  CALC_PC = 1.9* M* ((R_0/R)*(R_0/R))* (S_TG * S_TG * S_TG * S_TG) 

END FUNCTION CALC_PC
  
end module messy_ions
