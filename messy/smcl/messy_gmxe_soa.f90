MODULE MESSY_GMXE_SOA
  
!  USE MESSY_MAIN_CONSTANTS_MEM,   ONLY: dp, R_gas, N_a

  USE messy_main_constants_mem,          ONLY: dp, R_gas, N_a
  IMPLICIT NONE
  SAVE

  

!  INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12,307)
!  REAL(dp), PARAMETER :: R_gas   = 8.314409_dp  ! R [J/K/mol]
!  REAL(dp), PARAMETER :: N_A     = 6.022045E23_dp ! Avogadro constant [1/mol]
! SOA TYPES
  INTEGER            :: NSOA  = 0         ! number of SOA species
  INTEGER            :: NOXI  = 0         ! number of oxidants
  INTEGER            :: NGAS  = 0         ! number of gaseous precursors
  INTEGER            :: NSPEC = 0         ! total number of species in SOA calculations

  INTEGER, PARAMETER :: NSOA_REACT = 21   ! number of reactions producing SOA
  INTEGER, PARAMETER :: NGAS_MAX = 8
  INTEGER, PARAMETER :: NOXI_MAX = 3
  
! list of species to be excluded (or not contribute) to 
!         bulk mass on which the SOA is allowed to condense
!         Species are defined in the CTRL_GMXE_SOA namelist
  INTEGER, PARAMETER :: EXCLUDE_MAX_SOA = 10
  CHARACTER(Len=32)  :: EXCL_STR_SOA(EXCLUDE_MAX_SOA) = ''

  LOGICAL, PARAMETER :: L_TER_DIFF = .FALSE.

! following Zhang et al, JGR 2007

  ! SOA species
  INTEGER :: IDX_SOA_ISO1 =  0   ! from Isoprene
  INTEGER :: IDX_SOA_ISO2 =  0   ! from Isoprene
  INTEGER :: IDX_SOA_TOL1 =  0   ! from Toluene
  INTEGER :: IDX_SOA_TOL2 =  0   ! from Toluene
  INTEGER :: IDX_SOA_ALK1 =  0   ! from Alkanes
  INTEGER :: IDX_SOA_XYL1 =  0   ! from Xylene
  INTEGER :: IDX_SOA_XYL2 =  0   ! from Xylene
  INTEGER :: IDX_SOA_CRE1 =  0   ! from Cresol
  INTEGER :: IDX_SOA_TER1 =  0   ! from Terpenes
  INTEGER :: IDX_SOA_TER2 =  0   ! from Terpenes
  ! Oxidants
  INTEGER :: IDX_SOA_O3   =  0   ! Ozone
  INTEGER :: IDX_SOA_OH   =  0   ! OH
  INTEGER :: IDX_SOA_NO3  =  0   ! nitrate radical
  ! GAS precursors
  INTEGER :: IDX_SOA_ISO  =  0   ! Isoprene 1
  INTEGER :: IDX_SOA_TOL  =  0   ! Toluene  2
  INTEGER :: IDX_SOA_XYL  =  0   ! Xylene   3
  INTEGER :: IDX_SOA_CRE  =  0   ! Cresol   4
  INTEGER :: IDX_SOA_ALK  =  0   ! Alkanes  5
  INTEGER :: IDX_SOA_TER  =  0   ! Terpene  6
  INTEGER :: IDX_SOA_TERA =  0   ! Terpene(alpha pinene)  7
  INTEGER :: IDX_SOA_TERB =  0   ! Terpene(beta pinene)  8

  REAL(dp), POINTER, DIMENSION(:) :: MW              => NULL()
  REAL(dp), POINTER, DIMENSION(:) :: K_298           => NULL()
  REAL(dp), POINTER, DIMENSION(:) :: alpha           => NULL()
  REAL(dp), POINTER, DIMENSION(:) :: H_v             => NULL()
  REAL(dp), POINTER, DIMENSION(:) :: tend_fac        => NULL()
! allow SOA in this mode (of the aerosol module),
!       based on either potential 
!       volume or structure or species
!       controlled via namelist


  CHARACTER(LEN=20), DIMENSION(:), POINTER  :: names => NULL()

! first entry - treat species - YES/NO
! second entry - calculate FULL/SOA tendency as loss for 
!                gas species
!      controlled via namelist or host model
  LOGICAL :: L_GASPREC(NGAS_MAX,2)
  LOGICAL :: L_OXI(NOXI_MAX)
  INTEGER :: umode_soa = 0 
  INTEGER :: lmode_soa = 0

  REAL(dp) :: alpha_2(0:NSOA_REACT,2)
  REAL(dp) :: rate(0:NSOA_REACT)

CONTAINS
!-----------------------------------------------------------
  
  SUBROUTINE set_soa_params

    INTEGER :: I, counter

    NGAS = 0
    NOXI = 0
    counter = 0
    DO I=1,NGAS_MAX
      IF (L_GASPREC(I,1)) THEN
        NGAS = NGAS + 1
        SELECT CASE (I)
        CASE(1)
          counter = counter + 1
          IDX_SOA_ISO1 = counter
          counter = counter + 1
          IDX_SOA_ISO2 = counter
        CASE(2)
          counter = counter + 1
          IDX_SOA_TOL1 = counter
          counter = counter + 1
          IDX_SOA_TOL2 = counter
        CASE(3)
          counter = counter + 1
          IDX_SOA_XYL1 = counter
          counter = counter + 1
          IDX_SOA_XYL2 = counter
        CASE(4)
          counter = counter + 1
          IDX_SOA_Cre1 = counter
        CASE(5)
          counter = counter + 1
          IDX_SOA_ALK1 = counter
        CASE(6)
          counter = counter + 1
          IDX_SOA_TER1 = counter
          counter = counter + 1
          IDX_SOA_TER2 = counter
        END SELECT
      ENDIF
      NSOA = COUNTER
    ENDDO
    DO I=1,NOXI_MAX
      IF (L_OXI(I)) THEN
        NOXI = NOXI + 1
        SELECT CASE(I)
        CASE(1)
          IDX_SOA_O3 = NOXI + NSOA
        CASE(2)
          IDX_SOA_OH = NOXI + NSOA
        CASE(3)
          IDX_SOA_NO3 = NOXI + NSOA
        END SELECT
      ENDIF
    ENDDO
    NSPEC = NGAS + NOXI + NSOA
    counter = 0
    DO I=1,NGAS_MAX
      IF (L_GASPREC(I,1)) THEN
        SELECT CASE (I)
        CASE(1)
          counter = counter + 1
          IDX_SOA_ISO  = NSOA + NOXI + counter
        CASE(2)
          counter = counter + 1
          IDX_SOA_TOL  = NSOA + NOXI + counter
        CASE(3)
          counter = counter + 1
          IDX_SOA_XYL  = NSOA + NOXI + counter
        CASE(4)
          counter = counter + 1
          IDX_SOA_CRE  = NSOA + NOXI + counter
        CASE(5)
          counter = counter + 1
          IDX_SOA_ALK  = NSOA + NOXI + counter
        CASE(6)
          counter = counter + 1
          IDX_SOA_TER  = NSOA + NOXI + counter
        CASE(7)
          counter = counter + 1
          IDX_SOA_TERA = NSOA + NOXI + counter
        CASE(8)
          counter = counter + 1
          IDX_SOA_TERB = NSOA + NOXI + counter
        END SELECT

      ENDIF
    END DO

!    print*, IDX_SOA_TER, IDX_SOA_TERA, IDX_SOA_TERB, IDX_SOA_TER1, IDX_SOA_TER2
!    print*, IDX_SOA_ISO, IDX_SOA_ISO1, IDX_SOA_ISO2
!    print*, IDX_SOA_OH, IDX_SOA_O3, IDX_SOA_NO3
    ALLOCATE(MW(0:NSOA))
    ALLOCATE(K_298(0:NSOA))
    ALLOCATE(alpha(0:NSOA))
    ALLOCATE(H_v(0:NSOA))
    ALLOCATE(tend_fac(0:NSOA))
    ALLOCATE(names(0:NSOA+NOXI+NGAS))

    NAMES(IDX_SOA_ISO1) = "SOAIsop1"
    NAMES(IDX_SOA_ISO2) = "SOAIsop2"
    NAMES(IDX_SOA_TOL1) = "SOATol1"
    NAMES(IDX_SOA_TOL2) = "SOATol2"
    NAMES(IDX_SOA_XYL1) = "SOAXyl1"
    NAMES(IDX_SOA_XYL2) = "SOAXyl2"
    NAMES(IDX_SOA_CRE1) = "SOACre1"
    NAMES(IDX_SOA_ALK1) = "SOAAlk1"
    NAMES(IDX_SOA_TER1) = "SOATerp1"
    NAMES(IDX_SOA_TER2) = "SOATerp2"
    NAMES(IDX_SOA_OH)   = "OH"
    NAMES(IDX_SOA_O3)   = "O3"
    NAMES(IDX_SOA_NO3)  = "NO3"
    NAMES(IDX_SOA_ISO)  = "C5H8"
    NAMES(IDX_SOA_TOL)  = "Tol"
    NAMES(IDX_SOA_XYL)  = "Xyl"
    NAMES(IDX_SOA_CRE)  = "Cre"
    NAMES(IDX_SOA_ALK)  = "Alk"
    NAMES(IDX_SOA_TER)  = "Ter"
    NAMES(IDX_SOA_TERA) = "aPin"
    NAMES(IDX_SOA_TERB) = "bPin"


    MW(IDX_SOA_ISO1) = 177._dp
    MW(IDX_SOA_ISO2) = 177._dp
    MW(IDX_SOA_TOL1) = 150._dp
    MW(IDX_SOA_TOL2) = 150._dp
    MW(IDX_SOA_ALK1) = 150._dp
    MW(IDX_SOA_XYL1) = 150._dp
    MW(IDX_SOA_XYL2) = 150._dp
    MW(IDX_SOA_CRE1) = 150._dp
    MW(IDX_SOA_TER1) = 177._dp
    MW(IDX_SOA_TER2) = 177._dp
    
    ! in m^3 / ??g 
    K_298(IDX_SOA_ISO1) = 0.00459_dp
    K_298(IDX_SOA_ISO2) = 0.8628_dp
    K_298(IDX_SOA_TOL1) = 0.5829_dp
    K_298(IDX_SOA_TOL2) = 0.00209_dp
    K_298(IDX_SOA_ALK1) = 3.2227_dp
    K_298(IDX_SOA_XYL1) = 0.4619_dp
    K_298(IDX_SOA_XYL2) = 0.0154_dp
    K_298(IDX_SOA_CRE1) = 3.83_dp
    K_298(IDX_SOA_TER1) = 1.1561_dp
    K_298(IDX_SOA_TER2) = 0.0847_dp
    
    alpha(IDX_SOA_ISO1) = 0.232_dp
    alpha(IDX_SOA_ISO2) = 0.0288_dp
    alpha(IDX_SOA_TOL1) = 0.071_dp
    alpha(IDX_SOA_TOL2) = 0.138_dp
    alpha(IDX_SOA_ALK1) = 0.0718_dp
    alpha(IDX_SOA_XYL1) = 0.038_dp
    alpha(IDX_SOA_XYL2) = 0.0167_dp
    alpha(IDX_SOA_CRE1) = 0.05_dp
    alpha(IDX_SOA_TER1) = 0.0864_dp
    alpha(IDX_SOA_TER2) = 0.3857_dp

    tend_fac(1:NSOA) = 1._dp

    DO I=1,NGAS_MAX
      IF (L_GASPREC(I,1)) THEN
        SELECT CASE (I)
        CASE(1)
          IF (L_GASPREC(i,2)) THEN
            tend_fac(IDX_SOA_ISO1) = 1._dp/alpha(IDX_SOA_ISO1)
            tend_fac(IDX_SOA_ISO2) = 1._dp/alpha(IDX_SOA_ISO2)
          END IF
        CASE(2)
          IF (L_GASPREC(i,2)) THEN
            tend_fac(IDX_SOA_TOL1) = 1._dp/alpha(IDX_SOA_TOL1)
            tend_fac(IDX_SOA_TOL2) = 1._dp/alpha(IDX_SOA_TOL2)
          END IF
        CASE(3)
          IF (L_GASPREC(i,2)) THEN
            tend_fac(IDX_SOA_XYL1) = 1._dp/alpha(IDX_SOA_XYL1)
            tend_fac(IDX_SOA_XYL2) = 1._dp/alpha(IDX_SOA_XYL2)
          END IF
        CASE(4)
          IF (L_GASPREC(i,2)) THEN
            tend_fac(IDX_SOA_CRE1) = 1._dp/alpha(IDX_SOA_CRE1)
          END IF
        CASE(5)
          IF (L_GASPREC(i,2)) THEN
            tend_fac(IDX_SOA_ALK1) = 1._dp/alpha(IDX_SOA_ALK1)
          END IF
        CASE(6)
          IF (L_GASPREC(i,2)) THEN
            tend_fac(IDX_SOA_TER1) = 1._dp/alpha(IDX_SOA_TER1)
            tend_fac(IDX_SOA_TER2) = 1._dp/alpha(IDX_SOA_TER2)
          END IF
        CASE(7)
          IF (L_GASPREC(i,2)) THEN
            tend_fac(IDX_SOA_TER1) = 1._dp/alpha(IDX_SOA_TER1)
            tend_fac(IDX_SOA_TER2) = 1._dp/alpha(IDX_SOA_TER2)
          END IF
        CASE(8)
          IF (L_GASPREC(i,2)) THEN
            tend_fac(IDX_SOA_TER1) = 1._dp/alpha(IDX_SOA_TER1)
            tend_fac(IDX_SOA_TER2) = 1._dp/alpha(IDX_SOA_TER2)
          END IF
        END SELECT

      ENDIF
    END DO

    ! combination of Zhang and O'Donnell
    ! not used at the moment
    alpha_2( 1,1) = 0.0864_dp
    alpha_2( 2,1) = 0.0864_dp
    alpha_2( 3,1) = 0.0864_dp
    alpha_2( 4,1) = 0.0864_dp
    alpha_2( 5,1) = 0.0864_dp
    alpha_2( 6,1) = 0.0864_dp
    alpha_2( 7,1) = 0.071_dp
    alpha_2( 8,1) = 0.071_dp
    alpha_2( 9,1) = 0.071_dp
    alpha_2(10,1) = 0.038_dp
    alpha_2(11,1) = 0.038_dp
    alpha_2(12,1) = 0.038_dp
    alpha_2(13,1) = 0.232_dp
    alpha_2(14,1) = 0.232_dp
    alpha_2(15,1) = 0.232_dp
    alpha_2(16,1) = 0.05_dp
    alpha_2(17,1) = 0.05_dp
    alpha_2(18,1) = 0.05_dp
    alpha_2(19,1) = 0.0718_dp
    alpha_2(20,1) = 0.0718_dp
    alpha_2(21,1) = 0.0718_dp
    alpha_2( 1,2) = 0.3857_dp
    alpha_2( 2,2) = 0.3857_dp
    alpha_2( 3,2) = 0.3857_dp
    alpha_2( 4,2) = 0.3857_dp
    alpha_2( 5,2) = 0.3857_dp
    alpha_2( 6,2) = 0.3857_dp
    alpha_2( 7,2) = 0.138_dp
    alpha_2( 8,2) = 0.138_dp
    alpha_2( 9,2) = 0.138_dp
    alpha_2(10,2) = 0.0167_dp
    alpha_2(11,2) = 0.0167_dp
    alpha_2(12,2) = 0.0167_dp
    alpha_2(13,2) = 0.0288_dp
    alpha_2(14,2) = 0.0288_dp
    alpha_2(15,2) = 0.0288_dp
    alpha_2(16,2) = 0._dp
    alpha_2(17,2) = 0._dp
    alpha_2(18,2) = 0._dp
    alpha_2(19,2) = 0._dp
    alpha_2(20,2) = 0._dp
    alpha_2(21,2) = 0._dp
    
! values in kJ/mol, highest estimate Zhang et al.,
! other realistic values are down to 42 kJ/mol

    H_v(IDX_SOA_ISO1) = 156._dp
    H_v(IDX_SOA_ISO2) = 156._dp
    H_v(IDX_SOA_TOL1) = 156._dp
    H_v(IDX_SOA_TOL2) = 156._dp
    H_v(IDX_SOA_ALK1) = 156._dp
    H_v(IDX_SOA_XYL1) = 156._dp
    H_v(IDX_SOA_XYL2) = 156._dp
    H_v(IDX_SOA_CRE1) = 156._dp
    H_v(IDX_SOA_TER1) = 156._dp
    H_v(IDX_SOA_TER2) = 156._dp
    
  END SUBROUTINE set_soa_params
!---------------------------------------------------------
  SUBROUTINE SOA_driver(TEMP, dt, nmod, xt_aer, M_aer)

    INTEGER  :: nmod            ! note, this is the number of SOA_MODES
    REAL(dp) :: TEMP            ! temperature
    REAL(dp) :: dt              ! timestep length
    REAL(dp) :: xt_aer(0:nspec,0:nmod)  ! array of all species used in SOA model
    REAL(dp) :: M_aer(1:nmod)     ! mass of seed aerosol on which SOA condenses

 
    INTEGER  :: i

    ! calculate chemical production of potential SOA species (gas phase)
    xt_aer(0,:) = 0._dp
    CALL calc_soa_chem(TEMP, dt, xt_aer(:,0))
    xt_aer(0,:) = 0._dp
!    print*, "beforeconv", xt_aer(1:NSOA,0)
    ! convert xt_aer from molec/cm^3 -> ng/m^3
    DO I=1,nsoa
      xt_aer(I,:) = xt_aer(I,:) * (MW(I) * 1.e3_dp / N_a)
    END DO
    ! calculate phase partitioning
!    print*, "afterconv", xt_aer(1:NSOA,0)

    CALL calc_partitioning(TEMP, xt_aer, nmod, M_aer)
    ! backward conversion to molec/cm^3
    DO I=1,nsoa
      xt_aer(I,:) = xt_aer(I,:) / (MW(I) * 1.e3_dp / N_a)
      xt_aer(I,:) = MAX(xt_aer(I,:), 0._dp)
    END DO


    
  END SUBROUTINE SOA_driver

  !---------------------------------------------------------
  SUBROUTINE calc_SOA_chem(TEMP,dt,CHEM)
    
    Real(dp) :: TEMP
    REAL(dp) :: dt
    REAL(dp) :: CHEM(0:NSPEC) 
    CALL update_rates(TEMP)
    
    CALL calc_chem(dt, CHEM)

  CONTAINS
  !---------------------------------------------------------
    SUBROUTINE update_rates(TEMP)
      Real(dp) :: TEMP
      ! following Tsigaridis & Kanakidou, ACP, 2003)
      rate( 1) = f_arr(1.21e-11_dp,444._dp)    ! apinene_OH
      rate( 2) = f_arr(2.38e-11_dp,357._dp)    ! bpinene_OH
      rate( 3) = f_arr(1.01e-15_dp,-732._dp)   ! apinene_O3
      rate( 4) = 1.5e-17_dp                     ! bpinene_O3
      rate( 5) = f_arr(3.15e-13_dp,841._dp)    ! apinene_NO3
      rate( 6) = f_arr(1.6e-10_dp,-1248._dp)    ! bpinene_NO3
      
      rate( 7) = 5.96e-12_dp                    ! Toluene_OH
      rate( 8) = f_arr(2.34e-12_dp,-6694._dp)   ! Toluene_O3
      rate( 9) = 6.8e-17_dp                     ! Toluene_NO3
      
      rate(10) = 1.72e-11_dp                    ! Xylene_OH
      rate(11) = 1./3. * (2.40e-13_dp * exp(-5586._dp/TEMP) + &
        5.37e-13_dp * exp(-6039._dp/TEMP) + &
        1.91e-13_dp * exp(-5586._dp/TEMP))      ! Xylene_O3
      rate(12) = 3.54e-16_dp                    ! Xylene_NO3
      
      ! from O'Donnell et al, ACP, 2011
      rate(13) = f_arr(2.7E-11_dp,390._dp)      ! Isoprene_OH
      rate(14) = f_arr(1.03E-14_dp,-1995._dp)   ! Isoprene_O3
      rate(15) = f_arr(3.15E-12_dp,-450._dp)    ! Isoprene_NO3
      
      rate(16) = 0._dp                          ! Cresol_OH
      rate(17) = 0._dp                          ! Cresol_O3
      rate(18) = 0._dp                          ! Cresol_NO3
      
      rate(19) = 0._dp                          ! Alkanes_OH
      rate(20) = 0._dp                          ! Alkanes_O3
      rate(21) = 0._dp                          ! Alkanes_NO3
      
    END SUBROUTINE update_rates

    ELEMENTAL REAL(DP) FUNCTION f_arr(a,b)
      REAL(dp), INTENT(IN) :: a,b  ! K_ref, T dependence
      
      f_arr = a * exp (B/temp)
    END FUNCTION f_arr
      
  !-----------------------------------------------------------
    SUBROUTINE calc_chem(dt,CHEM)
      REAL(dp) :: dt
      REAL(dp) :: CHEM(0:NSPEC)
      REAL(dp) :: CHE(0:NSPEC)

      REAL(dp) :: delta_chem(0:nspec)
      REAL(dp) :: delta
      
      REAL(dp), PARAMETER :: sub_dt = 60._dp   ! 60 seconds
      REAL(dp) :: time, delta_t
      
      REAL(dp) :: ter_rat
      REAL(dp) :: sum_aft, sum_bef

      time = 0._dp
      CHEM(:) = MAX(CHEM(:), 0._dp)
      CHE(:)  = CHEM(:)
!      sum_bef = chem(idx_soa_iso) +  chem(idx_soa_iso1) + chem(idx_soa_iso2)
      DO WHILE (time < dt)
        ! isoprene loss processes
        
        delta_chem(:) = 0._dp
        
        delta = (rate(13) * chem(IDX_SOA_ISO) * CHEM(idx_SOA_OH) +      &
          rate(14) * chem(IDX_SOA_ISO) * CHEM(idx_SOA_O3) +      &
          rate(15) * chem(IDX_SOA_ISO) * CHEM(idx_SOA_NO3) )
        
        delta_chem(IDX_SOA_ISO1) = ALPHA(IDX_SOA_ISO1) * delta
        delta_chem(IDX_SOA_ISO2) = ALPHA(IDX_SOA_ISO2) * delta
        
!        delta_chem(IDX_SOA_ISO) = -1._dp * &
!          ( delta_chem(IDX_SOA_ISO1) + &
!          delta_chem(IDX_SOA_ISO2) )

        delta_chem(IDX_SOA_ISO) = -1._dp *                      &
          ( delta_chem(IDX_SOA_ISO1) * tend_fac(IDX_SOA_ISO1) + &
            delta_chem(IDX_SOA_ISO2) * tend_fac(IDX_SOA_ISO2) )

!        print*, delta,  delta_chem(IDX_SOA_ISO),  delta_chem(IDX_SOA_ISO1),&
!           delta_chem(IDX_SOA_ISO2), CHEM(IDX_SOA_ISO), CHEM(IDX_SOA_OH), IDX_SOA_OH

        ! terpenes
        IF (L_TER_DIFF) THEN
          ter_rat = chem(IDX_SOA_TERA) / &
            (chem(IDX_SOA_TERA) + chem(IDX_SOA_TERB))
        ELSE
          ter_rat = 0.7_dp
        END IF
        delta = rate( 1) * ter_rat * chem(IDX_SOA_TER) * CHEM(idx_SOA_OH) + &
          rate( 3) * ter_rat * chem(IDX_SOA_TER) * CHEM(idx_SOA_O3) + &
          rate( 5) * ter_rat * chem(IDX_SOA_TER) * CHEM(idx_SOA_NO3) 
        
        
        
        delta_chem(IDX_SOA_TER1) = ALPHA(IDX_SOA_TER1) * delta
        delta_chem(IDX_SOA_TER2) = ALPHA(IDX_SOA_TER2) * delta
        
!        delta_chem(IDX_SOA_TERA) = -1._dp * &
!          ( delta_chem(IDX_SOA_TER1) + &
!          delta_chem(IDX_SOA_TER2) )

        delta_chem(IDX_SOA_TERA) = -1._dp *                      &
          ( delta_chem(IDX_SOA_TER1) * tend_fac(IDX_SOA_TER1) + &
            delta_chem(IDX_SOA_TER2) * tend_fac(IDX_SOA_TER2) ) 
        
        delta = rate( 2) * (1._dp - ter_rat) * chem(IDX_SOA_TER) * CHEM(idx_SOA_OH) +  &
          rate( 4) * (1._dp - ter_rat) * chem(IDX_SOA_TER) * CHEM(idx_SOA_O3) +  &
          rate( 6) * (1._dp - ter_rat) * chem(IDX_SOA_TER) * CHEM(idx_SOA_NO3) 
        
        delta_chem(IDX_SOA_TER1) = delta_chem(IDX_SOA_TER1) + &
          ALPHA(IDX_SOA_TER1) * delta
        delta_chem(IDX_SOA_TER2) = delta_chem(IDX_SOA_TER2) + &
          ALPHA(IDX_SOA_TER2) * delta
        
!        delta_chem(IDX_SOA_TER) = -1._dp * &
!          ( delta_chem(IDX_SOA_TER1) + &
!          delta_chem(IDX_SOA_TER2) )
        
        delta_chem(IDX_SOA_TER) = -1._dp *                      &
          ( delta_chem(IDX_SOA_TER1) * tend_fac(IDX_SOA_TER1) + &
            delta_chem(IDX_SOA_TER2) * tend_fac(IDX_SOA_TER2) )

       delta_chem(IDX_SOA_TERB) = delta_chem(IDX_SOA_TER) - delta_chem(IDX_SOA_TERA)
        
        ! xylene
        delta = (rate(10) * chem(IDX_SOA_XYL) * CHEM(idx_SOA_OH) +      &
          rate(11) * chem(IDX_SOA_XYL) * CHEM(idx_SOA_O3) +      &
          rate(12) * chem(IDX_SOA_XYL) * CHEM(idx_SOA_NO3) )
        delta_chem(IDX_SOA_XYL1) = ALPHA(IDX_SOA_XYL1) * delta
        delta_chem(IDX_SOA_XYL2) = ALPHA(IDX_SOA_XYL2) * delta
        
!       delta_chem(IDX_SOA_XYL) = -1._dp * &
!         ( delta_chem(IDX_SOA_XYL1) + &
!         delta_chem(IDX_SOA_XYL2) )
        
        delta_chem(IDX_SOA_XYL) = -1._dp *                      &
          ( delta_chem(IDX_SOA_XYL1) * tend_fac(IDX_SOA_XYL1) + &
            delta_chem(IDX_SOA_XYL2) * tend_fac(IDX_SOA_XYL2) )

       ! Toluene
        delta = (rate( 7) * chem(IDX_SOA_TOL) * CHEM(idx_SOA_OH) +      &
          rate( 8) * chem(IDX_SOA_TOL) * CHEM(idx_SOA_O3) +      &
          rate( 9) * chem(IDX_SOA_TOL) * CHEM(idx_SOA_NO3) )
        delta_chem(IDX_SOA_TOL1) = ALPHA(IDX_SOA_TOL1) * delta
        delta_chem(IDX_SOA_TOL2) = ALPHA(IDX_SOA_TOL2) * delta
        
!       delta_chem(IDX_SOA_TOL) = -1._dp * &
!         ( delta_chem(IDX_SOA_TOL1) + &
!         delta_chem(IDX_SOA_TOL2) )
        delta_chem(IDX_SOA_TOL) = -1._dp *                      &
          ( delta_chem(IDX_SOA_TOL1) * tend_fac(IDX_SOA_TOL1) + &
            delta_chem(IDX_SOA_TOL2) * tend_fac(IDX_SOA_TOL2) )        

        ! Cresole
        delta = (rate(16) * chem(IDX_SOA_CRE) * CHEM(idx_SOA_OH) +      &
          rate(17) * chem(IDX_SOA_CRE) * CHEM(idx_SOA_O3) +      &
          rate(18) * chem(IDX_SOA_CRE) * CHEM(idx_SOA_NO3) )
        delta_chem(IDX_SOA_CRE1) = ALPHA(IDX_SOA_CRE1) * delta
        
!       delta_chem(IDX_SOA_CRE) = -1._dp * &
!         ( delta_chem(IDX_SOA_CRE1))
        delta_chem(IDX_SOA_CRE) = -1._dp *                      &
          ( delta_chem(IDX_SOA_CRE1) * tend_fac(IDX_SOA_CRE1) )

        ! Alkanes
        delta = (rate(19) * chem(IDX_SOA_ALK) * CHEM(idx_SOA_OH) +      &
          rate(20) * chem(IDX_SOA_ALK) * CHEM(idx_SOA_O3) +      &
          rate(21) * chem(IDX_SOA_ALK) * CHEM(idx_SOA_NO3) )
        delta_chem(IDX_SOA_ALK1) = ALPHA(IDX_SOA_ALK1) * delta
        
!       delta_chem(IDX_SOA_ALK) = -1._dp * &
!         ( delta_chem(IDX_SOA_ALK1) )
        delta_chem(IDX_SOA_ALK) = -1._dp *                      &
          ( delta_chem(IDX_SOA_ALK1) * tend_fac(IDX_SOA_ALK1) )
        
        ! integrate forward in time
        delta_t = sub_dt
        IF ((time + sub_dt) > dt) delta_t = dt - time
        delta_chem(:) = delta_chem(:) * delta_t
        CHEM(:) = CHEM(:) + delta_chem(:)
        DO WHILE (ANY(CHEM(:) < 0._dp)) 
           CHEM(:) = CHEM(:) - delta_chem(:)
           delta_chem(:) = delta_chem(:) / delta_t
           delta_t = delta_t/2._dp
           delta_chem(:) = delta_chem(:) * delta_t
           CHEM(:) = CHEM(:) + delta_chem(:)
        END DO
        time = time + delta_t
      END DO
!      sum_aft = chem(idx_soa_iso) +  chem(idx_soa_iso1) + chem(idx_soa_iso2)

!!$      IF ( abs(sum_bef - sum_aft)/sum_bef > 1.e-10_dp ) print*, &
!!$           " problem in chemistry", &
!!$           "after" ,chem(idx_soa_iso),chem(idx_soa_iso1),chem(idx_soa_iso2), &
!!$           "before" ,che(idx_soa_iso),che(idx_soa_iso1),che(idx_soa_iso2)

    END SUBROUTINE calc_chem

  END SUBROUTINE calc_SOA_chem
  
!----------------------------------------------------------------------
  SUBROUTINE calc_partitioning(TEMP,xt_aer,nmod,M_aer)
! phase partitioning of SOA compounds into the aerosol phase

    INTEGER  :: nmod   ! number of modes
    REAL(dp) :: TEMP   ! temperature
    REAL(dp) :: xt_aer(0:NSPEC,0:nmod)  ! chemical concentrations for all considered species
    REAL(dp) :: M_aer(nmod)  ! Mass of the seed aerosol per mode


    REAL(dp) :: K(NSOA)
! a) calculate phase partitioning coefficient
    CALL CALC_COEFF

! b) calculate actual phase partitioning
    CALL CALC_PART

!---------------------------------
  CONTAINS
    SUBROUTINE CALC_COEFF

      INTEGER :: i

      DO i=1,nsoa
        ! Factor 1000. originates from kJ -> J
        K(I) = K_298(i) * TEMP / 298._dp * &
              exp ( H_v(i)*1000. / R_gas * (1._dp/TEMP - 1._dp/298._dp) )     
      END DO
    END SUBROUTINE CALC_COEFF
!-----------------------------------
    SUBROUTINE CALC_PART

      INTEGER :: i, iter, j

      
      REAL(dp) :: M_aer_tot
      REAL(dp) :: M0_total
      REAL(dp) :: M0(nmod)

      REAL(dp) :: xt_AER_tot(1:nsoa,0:1)   ! in gas phase (0)
                                           ! or aerosol phase (1)
      REAL(dp) :: xt_tot

      REAL(dp) :: denom, fac, sum_bef, sum_aft

      DO i=1,nsoa
        ! total aerosol phase SOA
        xt_aer_tot(i,1) = SUM(xt_aer(i,1:nmod))
        ! gas phase SOA
        xt_aer_tot(i,0) = xt_aer(i,0)
!      print*, "in part I",  i, xt_aer_tot(i,0), xt_aer_tot(i,1)
      ENDDO
!!$      sum_bef = SUM(xt_aer_tot(1:NSOA,0)) + SUM(xt_aer_tot(1:NSOA,1))
      M_aer_tot = SUM(M_aer(1:nmod))
! test: remove the SOA fraction already here, as it is later on 
!       added in the loop and in the line below manually
!  total aerosol bulk mass without SOA
      M_aer_tot = M_aer_tot - SUM(xt_aer_tot(1:nsoa,1)) * 1.e-3_dp
!  total aerosol bulk mass, in which the SOA partitions
      M0_total  = M_aer_tot + SUM(xt_aer_tot(1:nsoa,1)) * 1.e-3_dp

      IF (M0_total < 1.e-10) RETURN

      DO iter=1,5
        ! starting with seed aerosol M_aer to condense on
        DO i=1,nsoa
          ! aerosol = gas phase * K(i) * M0
          ! total = gas and aerosol phase
          xt_tot = xt_aer_tot(i,0) + xt_aer_tot(i,1)
          ! denominator
          denom  = (1._dp + K(i) * M0_total)
          ! gas phase
          xt_aer_tot(i,0) = xt_tot / denom
          ! aerosol phase
          xt_aer_tot(i,1) = K(i) * M0_total * xt_tot / denom
        END DO
        ! update M0_total
        M0_total = M_aer_tot + SUM(xt_aer_tot(1:nsoa,1)) * 1.e-3_dp
!!$        sum_aft = SUM(xt_aer_tot(1:NSOA,0)) + SUM(xt_aer_tot(1:NSOA,1))
!!$        IF (abs((sum_aft - sum_bef)/sum_bef) > 1.e-5) &
!!$             print*, "change in mass", sum_aft, sum_bef,&
!!$             (sum_aft - sum_bef)/sum_bef, denom, xt_tot, &
!!$             M_aer_tot, M0_total
      END DO

      ! Redistribute total SOA over suitable MODES
      DO j=1,nmod
        fac = M_aer(j) / M_aer_tot
        DO i=1,nsoa
          xt_aer(i,j) = xt_aer_tot(i,1) * fac
        END DO
      ENDDO

!!$      sum_aft = SUM(xt_aer_tot(1:NSOA,0)) + SUM(xt_aer(1:NSOA,1)) + &
!!$                SUM(xt_aer(1:NSOA,2)) + SUM(xt_aer(1:NSOA,3)) + &
!!$                SUM(xt_aer(1:NSOA,4)) + SUM(xt_aer(1:NSOA,5)) + &
!!$                SUM(xt_aer(1:NSOA,6)) + SUM(xt_aer(1:NSOA,7))
!!$
!!$      IF (abs((sum_aft - sum_bef)/sum_bef) > 1.e-5) &
!!$           print*, "change in mass 2", sum_aft, sum_bef,&
!!$           (sum_aft - sum_bef)/sum_bef, denom, xt_tot, &
!!$           M_aer_tot, M0_total
      ! gas phase
      DO i=1,nsoa
        xt_aer(i,0) = xt_aer_tot(i,0)
      END DO

!!$      sum_aft = SUM(xt_aer(1:NSOA,0)) + SUM(xt_aer(1:NSOA,1)) + &
!!$                SUM(xt_aer(1:NSOA,2)) + SUM(xt_aer(1:NSOA,3)) + &
!!$                SUM(xt_aer(1:NSOA,4)) + SUM(xt_aer(1:NSOA,5)) + &
!!$                SUM(xt_aer(1:NSOA,6)) + SUM(xt_aer(1:NSOA,7))
!!$
!!$      IF (abs((sum_aft - sum_bef)/sum_bef) > 1.e-5) THEN
!!$         print*, "change in mass 3", sum_aft, sum_bef,&
!!$              (sum_aft - sum_bef)/sum_bef, denom, xt_tot, &
!!$               M_aer_tot, M0_total
!!$         STOP
!!$      ENDIF
    END SUBROUTINE CALC_PART

!------------------------------------
  END SUBROUTINE calc_partitioning


!-----------------------------------------------------------------------
END MODULE MESSY_GMXE_SOA
