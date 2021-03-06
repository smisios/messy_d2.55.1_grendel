MODULE MESSY_CLOUD_DROPLET

!       Author:  Holger Tost
!       last modified: 27.06.2006

  USE messy_main_constants_mem,      ONLY: dp

  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: cloud_droplet_rotstayn, cloud_droplet_menon
  PUBLIC :: cloud_droplet_jones
  PUBLIC :: cloud_droplet_ARG, cloud_droplet_lin

CONTAINS

!-----------------------------------------------------------

  SUBROUTINE cloud_droplet_rotstayn(nlev, kproma, slf, aersulf, cdnc)

    INTEGER,  INTENT(IN)    :: nlev, kproma
    REAL(dp), INTENT(IN)    :: slf(kproma)           !land-sea-fraction
 !aerosol sulfate in mug/m^3
    REAL(dp), INTENT(IN)    :: aersulf(kproma, nlev)
 !cloud droplet number in 1/m^3
    REAL(dp), INTENT(INOUT) :: cdnc(kproma, nlev)   

    REAL(dp) :: zcdnc, cdnc_land, cdnc_ocean
    INTEGER  :: jl,jk

    do jk=1,nlev
      do jl=1,kproma
        cdnc_land   = 173.8e6_dp*aersulf(jl,jk)**0.26_dp
        cdnc_ocean  = 114.8e6_dp*aersulf(jl,jk)**0.48_dp
        zcdnc =  cdnc_land * slf(jl) +  cdnc_ocean * (1._dp - slf(jl))
        cdnc(jl,jk) = zcdnc
      enddo
    end do

  END SUBROUTINE cloud_droplet_rotstayn

!-----------------------------------------------------------

  SUBROUTINE cloud_droplet_menon(nlev, kproma, slf, aersulf, aerss, aerom, cdnc)

    INTEGER,  INTENT(IN)    :: nlev, kproma
    REAL(dp), INTENT(IN)    :: slf(kproma)           !land-sea-fraction
 !aerosol sulfate in mug/m^3
    REAL(dp), INTENT(IN)    :: aersulf(kproma, nlev)
 !aerosol seasalt in mug/m^3
    REAL(dp), INTENT(IN)    :: aerss(kproma, nlev)
 !aerosol soluble organic matter in mug/m^3
    REAL(dp), INTENT(IN)    :: aerom(kproma, nlev)
 !cloud droplet number in 1/m^3
    REAL(dp), INTENT(INOUT) :: cdnc(kproma, nlev)   

    REAL(dp) :: zcdnc, cdnc_land, cdnc_ocean
    INTEGER  :: jl,jk

    do jk=1,nlev
      do jl=1,kproma
        cdnc_land  = 0._dp
        cdnc_ocean = 0._dp
        if (aersulf(jl,jk) > 0._dp) &
          cdnc_land   = aersulf(jl,jk)**0.5_dp 
        if (aerom(jl,jk) > 0._dp) &
          cdnc_land   = cdnc_land * aerom(jl,jk)**0.13_dp
        cdnc_land = cdnc_land * 257.0395783_dp

        if (aersulf(jl,jk) > 0._dp) &
          cdnc_ocean   = aersulf(jl,jk)**0.5_dp 
        if (aerom(jl,jk) > 0._dp) &
          cdnc_ocean   = cdnc_ocean * aerom(jl,jk)**0.13_dp
        if (aerss(jl,jk) > 0._dp) &
          cdnc_ocean   = cdnc_ocean * aerss(jl,jk)**0.05_dp
        cdnc_ocean = cdnc_ocean * 257.0395783_dp
        
        zcdnc =  cdnc_land * slf(jl) +  cdnc_ocean * (1._dp - slf(jl))
        cdnc(jl,jk) = zcdnc *1.e6_dp
      enddo
    end do

  END SUBROUTINE cloud_droplet_menon

!-----------------------------------------------------------

  SUBROUTINE cloud_droplet_jones(nlev, kproma, fac_sulf, aersulf, cdnc)
    
    USE messy_main_constants_mem,  ONLY: MS, M_air
 
    INTEGER,  INTENT(IN)    :: nlev, kproma
 !aerosol sulfate number in 1/m^3
    REAL(dp), INTENT(IN)    :: aersulf(kproma, nlev), fac_sulf
 !cloud droplet number in 1/m^3
    REAL(dp), INTENT(INOUT) :: cdnc(kproma, nlev)   

    REAL(dp) :: zcdnc,nsulf
    INTEGER  :: jl,jk

    do jk=1,nlev
      do jl=1,kproma
!  calculate sulfate numbers from aerosol mass
!  according to Jones et al. JGR, 2001
!  convert sulphate mass to S mass
!         ->   /fac_sulf * MS/M_air
!     from mug/m^3 back to kg/m^3 factor:  *  1.e-9
!     from mass to numbers with factor:    5.125e17
!                 resulting factor 5.125e8  
        nsulf = aersulf(jl,jk) * 5.125e8_dp / fac_sulf * MS / M_air
        zcdnc = 3.75e8_dp * (1._dp - exp(-2.5e-9_dp * nsulf) )
        cdnc(jl,jk) = MAX(zcdnc, 5.0e6_dp)
      enddo
    enddo

  END SUBROUTINE cloud_droplet_jones

!-----------------------------------------------------------

  SUBROUTINE cloud_droplet_ARG(temp, press, velo, aer_rad, aer_density,   &
                               sigma, number, mass,                       &
                               nu, Phi, eps, kappa, M_ap, SAT_crit, sup_sat_scheme, &
                               nfrac, S_max, aer_crit, N_activ,           &
                               kproma, klev, nmod, nspec,                 &
                               kfrac, mfrac, cmr2mmr)

! cloud droplet activation routine according to Abdul-Razzak & Ghan
! implementation for modal aerosol distributions
    USE messy_main_constants_mem,   ONLY: R_gas, M_H2O, rho_H2O, g,   &
                                          cp_air, PI, atm2Pa, Tmelt,  &
                                          rv, rd, M_air
    USE messy_main_tools,           ONLY: tlucuaw
    IMPLICIT NONE

! INPUT
    INTEGER,  INTENT(IN) :: kproma, klev        ! column, level indices
    INTEGER,  INTENT(IN) :: nmod, nspec         ! mode, species indices
    REAL(dp), INTENT(IN) :: TEMP                ! temperature
    REAL(dp), INTENT(IN) :: PRESS               ! pressure
    REAL(dp), INTENT(IN) :: aer_rad(nmod)       ! aerosol radius
    REAL(dp), INTENT(IN) :: aer_density(nspec)  ! aerosol density
    REAL(dp), INTENT(IN) :: sigma(nmod)         ! aerosol mode width
    REAL(dp), INTENT(IN) :: velo                ! vertical velocity
    REAL(dp), INTENT(IN) :: number(nmod)        ! aerosol number field
    REAL(dp), INTENT(IN) :: mass(nmod,nspec)    ! aerosol mass field
    REAL(dp), INTENT(IN) :: nu(nspec)           ! ion number of a species
    REAL(dp), INTENT(IN) :: eps(nspec)          ! mass fraction 
                                                ! of water soluble substance
    REAL(dp), INTENT(IN) :: Phi(nspec)          ! osmotic coefficient of species
    REAL(dp), INTENT(IN) :: kappa(nspec)        ! kappa parameter of species
    REAL(dp), INTENT(IN) :: M_ap(nspec)         ! molecular weight 
                                                ! of the aerosol species
    INTEGER,  INTENT(IN) :: sup_sat_scheme      ! scheme for calculating critical supersaturation
    REAL(dp), INTENT(INout) :: SAT_crit(nmod)   ! critical supersaturation per mode
                                                ! as provided by an external source or
                                                ! internally calculated
! OUTPUT
    REAL(dp) :: S_max
    REAL(dp) :: aer_crit(nmod)
    REAL(dp) :: nfrac(nmod), kfrac(nmod), mfrac(nmod), cmr2mmr(nmod)
    REAL(dp) :: N_activ
    
! LOCAL variables
    REAL(dp), PARAMETER :: alv   = 2.5008e6_dp ! latent heat for
    !                                          ! vaporisation in J/kg
    REAL(dp), PARAMETER :: ac    = 0.06_dp     ! accomodation coefficient for
                                               ! water vapour over water surface
                                               ! minimum diameter 
                                               ! for droplet distribution
    REAL(dp), PARAMETER :: Dp_high = 5._dp * 1.e-6_dp
                                               ! maximum diameter

    REAL(dp), PARAMETER :: c1  = 1.e-3_dp * g * M_H2O * alv / (cp_air * R_gas)
    REAL(dp), PARAMETER :: c2  = 1.e-3_dp * g * M_air / R_gas
    REAL(dp), PARAMETER :: c3  = 1.e3_dp * R_gas / M_H2O
    REAL(dp), PARAMETER :: c4  = (M_H2O * alv * alv) / (cp_air * M_air)
    REAL(dp), PARAMETER :: c5  = (2._dp * M_H2O) / (rho_H2O * R_gas) * 1.e-3_dp
    REAL(dp), PARAMETER :: c6  = M_H2O / rho_H2O
    REAL(dp), PARAMETER :: c7  = rho_H2O * R_gas * 1.e3_dp / M_H2O
    REAL(dp), PARAMETER :: c8  = rho_H2O * alv * alv * M_H2O / R_gas * 1.e-3_dp 
    REAL(dp), PARAMETER :: c9  = rho_H2O * alv
    REAL(dp), PARAMETER :: c10 = 2._dp * pi * rho_H2O

    REAL(dp), PARAMETER :: sqrt2 = 1.414213562373095049_dp

    REAL(dp) :: A            ! curvature effect
    REAL(dp) :: B, B1, B2    ! solute effect
    REAL(dp) :: S_m(nmod)    ! critical superation of a particle with radius a_m
    REAL(dp) :: tau          ! surface tension for water
    REAL(dp) :: p_sat        ! saturation water vapour pressure
    REAL(dp) :: D_v          ! diffusivity of water vapour
    REAL(dp) :: K_a          ! diffusivity of heat
    REAL(dp) :: D_star       ! corrected D_v for curvatured surface
    REAL(dp) :: K_star       ! corrected K_star for curvatured surface
    REAL(dp) :: growth       ! growth coefficient
    REAL(dp) :: growth1, growth2, growth3
    REAL(dp) :: eta, zeta, et
    REAL(dp) :: alpha, gamma
    REAL(dp) :: f1, f2
    REAL(dp) :: N_help1, N_help2, sum_help
    REAL(dp) :: VKAPPA, EKAPPA, Ds, Dt, Dt2, Dt3, xx, AA, BB, CC, t3, t32, Sold, Snew

    INTEGER  :: it, jt, jm
   
    REAL(dp) :: B_aux      ! auxillary parameter for calculation of radius 
                           ! distribution of droplets

    REAL(dp) :: c11        ! constant with exponential expression
    REAL(dp) :: Dp_low     ! minimum diameter for droplet distribution

    Dp_low = MIN(0.207683_dp * ac**(-0.33048_dp), 5._dp) * 1.e-6_dp

    c11 = SQRT(2._dp * pi * M_H2O / R_gas) 

    S_max       = 0._dp
    S_m(:)      = 0._dp
    nfrac(:)    = 0._dp
    mfrac(:)    = 0._dp
    kfrac(:)    = 0._dp
    aer_Crit(:) = 0._dp
    sum_help    = 0._dp

    ! surface tension calculated with the Eoetvoes rule (from Wikipedia)
    ! converted to [N/m]
    ! Formual after Pruppacher & Klett 1978, Eq. 5-12
    tau = (76.10_dp - 0.155_dp * (temp - tmelt)) * 1.e-3_dp


    ! saturation vapour pressure
    ! using the E5/MESSy lookuptables and 
    ! rescale with the additional factor *rv/rd (opposite to the lookup table)
    it = INT(temp *1000._dp)
    p_sat = tlucuaw(it) * rv/rd   ! in Pascal

    ! curvature term of Koehler equation
    ! unit conversion results in a factor of: 10^-3
    ! (N/m * g/mol) / (kg/m^3 * J/K/mol * K)               =
    ! (kg m g m^3 s^2 K mol) / (m s^2 mol kg kg m^2 K)     =
    ! (m^4 g) / (kg m^3) = m * g/kg = m *10^-3 kg/kg
    ! PARAMETER :: c5 = (2._dp * M_H2O) / (rho_H2O * R_gas) * 1.e-3_dp
    A = c5 * tau / temp

    ! gas diffusivity of water vapour 
    ! (Pruppacher & Klett, 2000, eq. 13-3, page 503)
    ! unit is in cm^2/s
    ! conversion to m^2/s => 10^-4
    D_v = 0.211_dp * (temp/tmelt)**1.94 * (atm2Pa/press) * 1.e-4_dp


    ! correct D_v for droplet distribution effects
    ! follow Fountoukis & Nenes , 2005
    ! auxiliary coefficient [m]
    B_aux = 2._dp * D_v / ac * c11 * (temp **(-0.5))

    D_star = D_v / (Dp_high - Dp_low) * &
             (Dp_high - Dp_low -        &
             B_aux * LOG((Dp_high + B_aux)/(Dp_low + B_aux)))

    ! heat diffusivity (Pruppacher & Klett, 2000, eq. 13-18, page 508)
    ! unit is cal cm^-1 s^-1 K
    ! conversion to J/(m s K) : 100/4.1867
    ! K_a = ( (TEMP - Tmelt) * 0.017_dp + 5.69) * 1.e-5_dp * 100._dp * 4.1867_dp
    K_a = ( (TEMP - Tmelt) * 0.017_dp + 5.69) * 1.e-3_dp * 4.1867_dp

    ! size invariant coefficients from equation of 
    ! rate change of supersaturation (Leaitch et al., 1986) 
    ! of an adiabatically rising air parcel with upward speed V
    ! for alpha conversion with 10^-3 required because of units
    ! (m/s^2 * g/mol * J/kg) / (J/kg/K * J/K/mol * K * K) - 
    ! (m/s^2 * g/mol) / (J/K/mol * K)                     =
    ! (m*g)/(J*s^2) - (m*g)/(J*s^2)                       = 
    ! g / (kg * m)                                        =
    ! 10^-3 m^-1 kg / kg                                     =
    ! 10^-3 * 1/m
    ! PARAMETER :: c1 = 1.e-3_dp * g * M_H2O * alv) / (cp_air * R_gas)
    ! PARAMETER :: c2 = 1.e-3_dp * g * M_air / R_gas
    alpha = ( c1 / (temp * temp) ) - &
            ( c2 / temp )

    ! for gamma conversion with 10^3 for first term required because of units
    !   (J/K/mol * K) / (kg*m/s^2/m^2 * g/mol)                   +
    ! (g/mol * J^2/kg^2) / (J/K/kg * kg*m/s^2/m^2 * g/mol *K)  
    ! =  (kg m^2 K s^2 m^2 mol) / (s^2 K mol kg m g)             +
    ! (g kg m^2 kg m^2 s^2 K kg s^2 m^2 mol)/(mol s^2 kg s^2 kg kg m^2 kg m g K)
    ! = (m^3 / g) + (m^3 / kg)
    ! = 10^3 m^3/kg + m^3/kg
    !PARAMETER :: c3 = 1.e3_dp * R_gas / M_H2O
    !PARAMETER :: c4 = (M_H2O * alv * alv) / (cp_air * M_air)
    gamma = ( c3 * temp / p_sat )          +           &
            ( c4  / (press * temp) )
  
    ! correction for diffusivity over curvatured surface
    K_star = K_a/1._dp
    
    ! growth coefficient
    ! 1) (kg/m^3 * J/K/mol * K) / (kgm/s^2/m^2 * g/mol * m^2/s)    =
    !    (kg kg m^2 s^2 m^2 mol s) / (m^3 s^2 mol kg m g m^2)      = 
    !    (kg s) / (g m^2) = 10^3 s/m^2
    ! unit conversion: term one must be multiplied with 10^3
    ! 2) (kg/m^3 * J/kg) / (K * J/m/s/K) * (J/kg * g/mol) / (J/K/mol * K) =
    !    (kg J m s K J g K mol) / (m^3 kg K J kg mol J K)                 =
    !    (g s) / (m^2 kg) = 10^-3 s/m^2
    ! unit convection: term two must be multiplied with 10^-3
    ! 3) (kg/m^3 * J/kg) / (K * J/m/s/K)          =
    !    (kg J m s K) / (m^3 kg K J)              =
    !    (s/m^2)
    ! unit conversion: term three is o.k.
    ! PARAMETER :: c7 = rho_H2O * R_gas * 1.e3_dp / M_H2O
    ! Growth1 = (rho_H2O * R_gas * temp) / (p_sat * D_star * M_H2O)  * 1.e3_dp 
    Growth1 = (c7 * temp) / (p_sat * D_star)
    ! PARAMETER :: c8 = rho_H2O * alv * alv * M_H2O / R_gas * 1.e-3_dp 
    ! Growth2 = (rho_H2O * alv)/(K_star * temp) *         &
    !           (alv * M_H2O) / (temp*R_gas) * 1.e-3_dp 
    Growth2 = c8 / (K_star * temp * temp)
    ! PARAMETER :: c9 = rho_H2O * alv
    ! Growth3 = (rho_H2O * alv)/(K_star * temp)
    Growth3 = c9 / (K_star * temp)
    Growth  = Growth1 + Growth2 - Growth3
    Growth  = 1._dp / growth

   
    et   = alpha * velo/growth
      
    zeta = 2._dp/3._dp * A * et**0.5

    DO jm = 1, nmod

      SELECT CASE (sup_sat_scheme)
! different ways to calculate the critical supersaturation of the
! aerosol
        
      CASE(1)
 
         ! Abdul-Razzak and Ghan (2000) - Eq. 9
        
        ! solute (Raoult) term of Koehler equation (aerosol composition)
        ! PARAMETER :: c6 = M_H2O / rho_H2O
        !  B = (c6 * nu(jt) * Phi(jt) * eps(jt) * aer_density(jt)) / M_ap(jt)
        B1 = 0._dp
        B2 = 0._dp
        do jt=1,nspec 
          B1 = B1 + (mass(jm,jt) * nu(jt) * Phi(jt) * eps(jt) / M_ap(jt))
          B2 = B2 + (mass(jm,jt) / aer_density(jt))
        enddo
        B = c6 * B1 / B2
        
        ! critical superssaturation of a particle with geometric mean radius 
        ! of the aerosol distribution
        S_m(jm) = 2._dp / B**0.5_dp * &
          (A / (3._dp * MAX(aer_rad(jm),1.e-15_dp)) )**1.5_dp
      CASE(2)
        ! supersaturation taken from GMXE calculated slightly differently, making
        ! use of the kappa parameter of the aerosol
        
        S_M(jm) = SAT_crit(jm)

      CASE(3)

         ! Petters und Kriedenweis (2007), based on case 2 and generalized
         
         ! Volume-averaged kappa parameter (Eq. 7)
         N_help1 = 0.
         N_help2 = 0.
         do jt=1,nspec
            N_help1 = N_help1 + mass(jm,jt) * press / temp / R_gas * &
                 M_ap(jt) / 1000. / aer_density(jt) * kappa(jt)
            N_help2 = N_help2 + mass(jm,jt) * press / temp / R_gas * &
                 M_ap(jt) / 1000. / aer_density(jt)
         end do
         VKAPPA = N_help1/N_help2

         ! Argument of the exponent in (Eq. 1), diameter is introduced later
         EKAPPA = 4.0_dp * tau * M_H2O/1000. / (R_gas * temp * rho_H2O)

         ! Convert radius to diameter 
         Ds  = 2. * aer_rad(jm)

         ! Original GMXe code
         Dt   = Ds 
         Sold = 0.0
         Snew = 0.1

         AA = Ds * Ds * Ds
         BB = AA * (1._dp - VKAPPA)
         t32 = VKAPPA * AA / EKAPPA
         CC = -3._dp * t32
         t3 = sqrt(t32)
         Dt = t3 + (AA + BB - t3 * (CC + t32))**0.33333333333333333_dp  ! Set Dt to the starting point
         xx = Dt
         do while (abs(xx / Dt) > 1.0E-15_dp)
            Dt2 = Dt * Dt
            Dt3 = Dt2 * Dt
            xx = (CC * Dt2 * Dt2 + (Dt3 - AA) * (Dt3 - BB)) / &
                 (4 * CC * Dt3 + 3 * Dt2 * (2 * Dt3 - AA - BB))
            Dt = Dt - xx
         enddo
         Snew = (Dt3 - AA) / (Dt3 - BB) * exp(EKAPPA / Dt)
         S_m(jm) = Snew - 1

      END SELECT
      SAT_crit(jm) = s_m(jm)

!      print*, "after subsat ", jm
      ! dimensionless terms for calculation of S_max
      ! eta   = ((alpha * velo/growth)**1.5_dp) /       &
      !         (2._dp * pi * rho_H2O * gamma * N_ap)
      ! c10 = 2._dp * pi * rho_H2O
      eta  = (et**1.5_dp) /       &
        (c10 * gamma * number(jm))
      
      ! functions of standard deviation of the aerosol modes
      f1 = 0.5_dp * exp(2.5_dp * log(sigma(jm)) * log(sigma(jm)))
      f2 = 1._dp + 0.25_dp * log(sigma(jm))
      
      ! fraction of activation
      IF (et > 0._dp) THEN
        N_help1 = ( f1 * (zeta / eta)**1.5_dp ) +                        &
          ( f2 * ( (S_m(jm) * S_m(jm))/(eta + 3._dp * zeta) )**0.75_dp)
        IF (S_m(jm) > 0._dp) &
          sum_help   = sum_help +  (1._dp / (S_m(jm)*S_m(jm)) * N_help1)      
      ENDIF
    ENDDO
    ! calculation of S_max
    IF (sum_help > 0._dp ) S_max = 1._dp / sqrt(sum_help)
    
    
    N_activ = 0._dp
    do jm=1,nmod     
      IF (S_m(jm) > 0._dp .and. S_max > 0._dp) THEN
        ! calculation of critical aerosol radius
        aer_crit(jm) = (S_m(jm) / S_max)**(2._dp/3._dp) * aer_rad(jm)
        N_help2      = 2._dp * log( S_m(jm) / S_max ) / &
                      (3._dp * sqrt2 * log(sigma(jm)))
    ! calculation of activated number fraction per mode
        ! ARG_2000_EQ 13
        nfrac(jm)  = 0.5_dp * erfc05(N_help2) 
        ! ARG_2000_EQ 14
        mfrac(jm) = 0.5_dp * erfc05(N_help2 - 1.5 * sqrt2 *log(sigma(jm)) )
    ! calculation of activated particles (total) = CDNC
        N_activ    = N_activ + nfrac(jm) * number(jm)

      ENDIF
    enddo
    
  END SUBROUTINE cloud_droplet_ARG

!---------------------------------------------------------------------

  SUBROUTINE cloud_droplet_Lin(nmod, number, velo, N_old,   &
                               sigma, radius,               &
                               N_activ, N_a, temp,          &
                               lfrac, mfrac, cmr2mmr)

! calculate activation according to Lin and Leaitch (1997)
! see also Lohmann et al., ACP, 2007
! taken from ECHAM5-HAM by U. Lohmann, P. Stier, etc....
    USE messy_cloud_ori,    ONLY: cthomi
    
    INTEGER,  INTENT(IN)    :: nmod
    REAL(dp), INTENT(IN)    :: number(nmod)
    REAL(dp), INTENT(IN)    :: sigma(nmod)
    REAL(dp), INTENT(IN)    :: radius(nmod)
    REAL(dp), INTENT(IN)    :: velo
    REAL(dp), INTENT(IN)    :: cmr2mmr(nmod)

    ! number of activated aerosol particles 
    REAL(dp), INTENT(INOUT) :: N_activ
    ! number of aerosol particles potentially activated
    REAL(dp), INTENT(INOUT) :: N_a
    ! fraction of potentially activated numbers
    REAL(dp), INTENT(INOUT) :: lfrac(nmod), mfrac(nmod)
    REAL(dp), INTENT(IN)    :: N_old, temp
    
    REAL(dp), PARAMETER :: zeps = EPSILON(1._dp)
    REAL(dp), PARAMETER :: c1 = 2.3e-10_dp   ! m^4 s^-1
    REAL(dp), PARAMETER :: c2 = 1.27_dp
    REAL(dp), PARAMETER :: crcut=0.035E-6_dp ! Assumed lower cut-off of the
                                             ! aerosol size distribution [m]

    REAL(dp)  :: ztn, value1, val1, n_help,n_help_max
    REAL(dp)  :: zdummy, kfrac(nmod), rcut, rcut_max, rcut_min, r_old
    INTEGER   :: jmod, counter
    
    N_activ = 0._dp
! determining N_a (larger than cut-off radius)
    
    N_a   = 0._dp
    lfrac(:) = 0._dp
    kfrac(:) = 0._dp
    mfrac(:) = 0._dp

    do jmod = 1,nmod
      IF (radius(jmod)>zeps) THEN

       !--- Transform number distribution to error function:

        ztn=( LOG(crcut)-LOG(radius(jmod)) ) / LOG(sigma(jmod))

       !--- Calculate the cumulative of the log-normal number distribution:

        CALL error_function_limited(ztn,zdummy,kfrac(jmod))
!        print*, lfrac, erfc05(ztn)
      ENDIF
      
      N_a = N_a + kfrac(jmod) * number(jmod)
      
    ENDDO

! calculate inner bracket
    IF ( (temp > cthomi) .AND. &
         (velo > 0._dp)       .AND. &
         (N_a  > zeps ) )     THEN

      val1 = MAX( (N_a - N_old), 0._dp )
      
      value1 = (val1 * velo)/(velo + c1 * val1)

      N_activ = 0.1e6_dp * ((1.e-6_dp * value1)**c2)
    ENDIF

! determine activated fractions per mode
! iterative approach
! use a new cutoff radius as long until an agreement between N_activ and the particles 
! larger than the cut off radius is found.
    IF (N_activ > 0._dp) THEN
      rcut = 10._dp * crcut
      rcut_min = crcut
      rcut_max = 10._dp * rcut
      counter = 0 
      N_help = 0._dp
      DO WHILE (abs(N_help - N_activ) > 0.001_dp * N_activ)
        IF (counter > 100) EXIT
        N_help = 0._dp
        kfrac(:) = 0._dp
        counter = counter + 1
        DO jmod = 1,nmod
          IF (radius(jmod)>zeps) THEN

            !--- Transform number distribution to error function:

            ztn=( LOG(rcut)-LOG(radius(jmod)) ) / LOG(sigma(jmod))

            !--- Calculate the cumulative of the log-normal number distribution:

            CALL error_function_limited(ztn,zdummy,lfrac(jmod))

            N_help = N_help + lfrac(jmod) * number(jmod)           
          ENDIF
        END DO

        r_old = rcut
        IF (N_help > N_activ) THEN
          rcut_min = rcut
          rcut     = (rcut_max + rcut)/2._dp
          IF (n_help > n_help_max) rcut_max = rcut_max * 2.
          n_help_max = MAX(n_help_max,n_help)
        ELSE
          rcut_max = rcut
          rcut = (rcut + rcut_min)/2._dp
        ENDIF
      END DO
!      print*, counter

! determine the corresponding mass fraction for aerosol nucleation scavenging
      DO jmod = 1,nmod
        IF (radius(jmod)>zeps) THEN

       !--- Transform number distribution to error function:

          ztn=( LOG(r_old)-LOG(radius(jmod)*cmr2mmr(jmod)) ) / LOG(sigma(jmod))

       !--- Calculate the cumulative of the log-normal number distribution:

          CALL error_function_limited(ztn,zdummy,mfrac(jmod))
!        print*, lfrac, erfc05(ztn)
        ENDIF
                  
      ENDDO
    ENDIF

  END SUBROUTINE cloud_droplet_Lin

!-------------------------------------------------------------------------------
!*******************************************************************************
! Obtain erfc(x) (and erf(x)) to the maximum machine accuracy possible.
!
! REFERENCE: Irene A. Stegun and Ruth Zucker, J. Res. NBS  74B, 211-224 (1970).
!            ``Automatic Computing Methods for Special Functions'',
!            Part 1. Error, Probability, and Related Functions
!
!            CARDS GIVEN AS DOUBLE PRECISION CONSTANTS, YIELDS DOUBLE
!            PRECISION RESULTS ON THE UNIVAC 1108 SYSTEM.
!
!            SINGLE PRECISION RESULTS ARE OBTAINED BY
!            (1) DELETING THE TYPE STATEMENT
!            (2) SETTING NBC=8, NBM=27 ON THE FIRST DATA CARD
!            (3) CHANGING THE D'S TO E'S ON THE SECOND AND THIRD DATA
!            CARDS.
!            (4) CHANGING THE FUNCTION TYPE - DABS TO ABS, DEXP TO EXP
!
!            FOR OTHER VALUES OF NBM, THE CONSTANT ON THE THIRD DATA
!            CARD SHOULD BE GIVEN TO THE CORRESPONDING PRECISION.
!
! NBC:       Number of binary digits in the characteristic of a floating 
!            point number.
!            Single precision: NBC = 8, Double precision, NBC = 11

! NBM:       Accuracy desired or maximum number of binary digits in the 
!            mantissa of a floating point number.
!            Single precision: NBM = 11, Double precision, NBM = 60
!
! METHOD:    Power series for ABS(X) <= 1 (= ULPS, upper limit for series).
!            Continued fraction for  ULPS < ABS(X) <= ULCF 
!            (ULCF = upper limit of continued fraction)
!
! RANGE:     ABS(X) <= ULCF ; ULCF = 0.83*(2.0**((NBC-1)/2))
!            ABS(ULCF) approximately    9.3 for NBC = 8
!                                      26.5 for NBC = 11
!            Beyond this range the limiting values of the functions are
!            supplied: ERF(INF) = 1.0,   ERFC(INF) = 0.0
!
! ACCURACY:  Routine will return (NBM-I-3) significant binary digits
!            where I is the number of binary digits representing the
!            integer part of X**2. (This is essentially the accuracy of
!            the exponential routine.)
!
! PRECISION: Variable - by setting desired NBM
!*******************************************************************************

      function erfc05(x) result(erfc)

      implicit none
      real(dp), intent(in)    :: x
      real(dp)                :: erf,erfc
      real(dp)                :: an,bn,c1,dn,f,fn,  &
                                 fnm1,fnm2,gn,gnm1,gnm2,  &
                                 prev,pwr,rnbc,scf,sum,tn,toler, &
                                 ulcf,wn,y,ysq
      integer, parameter      :: nbc = 11
      integer, parameter      :: nbm = 60
      real(dp), parameter     :: zero = 0.0_dp
      real(dp), parameter     :: one = 1.0_dp
      real(dp), parameter     :: two = 2.0_dp
      real(dp), parameter     :: four = 4.0_dp
      real(dp), parameter     :: ulps = 1.0_dp
      real(dp), parameter     :: cons = 0.83_dp
      real(dp), parameter     :: trrtpi = 1.128379167095512574_dp


      rnbc = nbc
      toler = two**(-nbm)

!     Test on zero
      if (x == zero) then
          erf = zero
          erfc = one
          return
      endif

      y = abs(x)
      ysq = y**2

      if (y <= ulps) then
!         Method: power series
          sum = zero
          dn = one
          tn = one
          pwr = two*ysq
          do
              dn = dn + two
              tn = pwr*tn/dn
              sum = tn + sum
!             Tolerance check
              if (tn < toler) exit
          enddo
          erf = (sum + one)*trrtpi*y*exp(-ysq)
          erfc = one - erf
      else
!         Maximum argument
          c1 = two**((rnbc-one)/two)
          ulcf = cons*c1

!         Scale factor
          scf = two**(c1**2 - rnbc)

!         Limiting value
          if (y > ulcf) then
              erf = one
              erfc = zero
          else
!             Method: continued fraction
              fnm2 = zero
              gnm2 = one
              fnm1 = two*y
              gnm1 = two*ysq + one
              prev = fnm1/gnm1
              wn = one
              bn = gnm1 + four
              do
                  an = -wn*(wn + one)
                  fn = bn*fnm1 + an*fnm2
                  gn = bn*gnm1 + an*gnm2
                  f = fn/gn

!                 Tolerance check
                  if (abs(one - (f/prev)) <= toler) exit
                  if (prev > f) then
                      f = prev
                      exit
                  endif

!                 Both fn and gn must be tested if abs(x) < 0.61
                  if (gn >= scf) then
!                     Scaling
                      fn = fn/scf
                      gn = gn/scf
                      fnm1 = fnm1/scf
                      gnm1 = gnm1/scf
                  endif

                  fnm2 = fnm1
                  gnm2 = gnm1
                  fnm1 = fn
                  gnm1 = gn

                  wn = wn + two
                  bn = bn + four
                  prev = f
              enddo
              erfc = f*exp(-ysq)*trrtpi/two
              erf = one - erfc
          endif
      endif

!     Negative argument
      if (x < zero) erfc = two - erfc
!     if (x < zero) erf = -erf

      end function erfc05

!-----------------------------------------------------------
  SUBROUTINE error_function_limited ( arg, RESULT, ccum )
  !
  !****************************************************************************
  !
  !! CUMNOR computes the cumulative normal distribution.
  !
  !
  !     the integral from -infinity to x of
  !          (1/sqrt(2*pi)) exp(-u*u/2) du
  !
  !  Author:
  !  -------
  !  Original source:
  !
  !    W. J. Cody    Mathematics and Computer Science Division
  !                  Argonne National Laboratory
  !                  Argonne, IL 60439
  !
  !    DCDFLIB is attributed to Barry Brown, James Lovato, and Kathy Russell
  !            bwb@odin.mda.uth.tmc.edu.
  !
  !    Adopted to ECHAM/M7:
  !
  !    Philip Stier  (MPI-MET)                    2001
  !
  !
  !  Reference:
  !  ----------
  !
  !    W D Cody, 
  !    "ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special 
  !    Function Routines and Test Drivers"
  !    ACM Transactions on Mathematical Software,
  !    Volume 19, 1993, pages 22-32.
  !
  !  Parameters:
  !
  !     ARG --> Upper limit of integration.
  !                                        X is REAL(dp)
  !
  !     RESULT <-- Cumulative normal distribution.
  !                                        RESULT is REAL(dp)
  !
  !     CCUM <-- Complement of Cumulative normal distribution.
  !                                        CCUM is REAL(dp)
  !
  !
  ! Original Comments:
  !
  !
  ! This function evaluates the normal distribution function:
  !
  !                              / x
  !                     1       |       -t*t/2
  !          P(x) = ----------- |      e       dt
  !                 sqrt(2 pi)  |
  !                             /-oo
  !
  !   The main computation evaluates near-minimax approximations
  !   derived from those in "Rational Chebyshev approximations for
  !   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
  !   This transportable program uses rational functions that
  !   theoretically approximate the normal distribution function to
  !   at least 18 significant decimal digits.  The accuracy achieved
  !   depends on the arithmetic system, the compiler, the intrinsic
  !   functions, and proper selection of the machine-dependent
  !   constants.
  !
  !  Explanation of machine-dependent constants.
  !
  !   MIN   = smallest machine representable number.
  !
  !   EPS   = argument below which anorm(x) may be represented by
  !           0.5  and above which  x*x  will not underflow.
  !           A conservative value is the largest machine number X
  !           such that   1.0 + X = 1.0   to machine precision.
  !
  !  Error returns
  !
  !  The program returns  ANORM = 0     for  ARG .LE. XLOW.
  !
  !  Author: 
  !
  !    W. J. Cody
  !    Mathematics and Computer Science Division
  !    Argonne National Laboratory
  !    Argonne, IL 60439
  !
  !  Latest modification: March 15, 1992
  !
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: a = (/ &
       2.2352520354606839287d00, &
       1.6102823106855587881d02, &
       1.0676894854603709582d03, &
       1.8154981253343561249d04, &
       6.5682337918207449113d-2 /)
  REAL(dp) arg
  REAL(dp), PARAMETER, DIMENSION ( 4 ) :: b = (/ &
       4.7202581904688241870d01, &
       9.7609855173777669322d02, &
       1.0260932208618978205d04, &
       4.5507789335026729956d04 /)
  REAL(dp), PARAMETER, DIMENSION ( 9 ) :: c = (/ &
       3.9894151208813466764d-1, &
       8.8831497943883759412d00, &
       9.3506656132177855979d01, &
       5.9727027639480026226d02, &
       2.4945375852903726711d03, &
       6.8481904505362823326d03, &
       1.1602651437647350124d04, &
       9.8427148383839780218d03, &
       1.0765576773720192317d-8 /)
  REAL(dp) ccum
  REAL(dp), PARAMETER, DIMENSION ( 8 ) :: d = (/ &
       2.2266688044328115691d01, &
       2.3538790178262499861d02, &
       1.5193775994075548050d03, &
       6.4855582982667607550d03, &
       1.8615571640885098091d04, &
       3.4900952721145977266d04, &
       3.8912003286093271411d04, &
       1.9685429676859990727d04 /)
  REAL(dp) del
!@@@ REAL(dp) dpmpar
  REAL(dp) eps
  INTEGER i
  REAL(dp) min
  REAL(dp), PARAMETER, DIMENSION ( 6 ) :: p = (/ &
       2.1589853405795699d-1, &
       1.274011611602473639d-1, &
       2.2235277870649807d-2, &
       1.421619193227893466d-3, &
       2.9112874951168792d-5, &
       2.307344176494017303d-2 /)
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: q = (/ &
       1.28426009614491121d00, &
       4.68238212480865118d-1, &
       6.59881378689285515d-2, &
       3.78239633202758244d-3, &
       7.29751555083966205d-5 /)
  REAL(dp) RESULT
  REAL(DP), PARAMETER :: root32 = 5.656854248E0_dp
  REAL(DP), PARAMETER :: sixten = 16.0_dp
  REAL(DP) temp
  REAL(DP), PARAMETER :: sqrpi = 3.9894228040143267794E-1_dp
  REAL(DP), PARAMETER :: thrsh = 0.66291E0_dp
  REAL(DP) x
  REAL(DP) xden
  REAL(DP) xnum
  REAL(DP) y
  REAL(DP) xsq
  !
  !  Machine dependent constants
  !
  eps = EPSILON ( 1.0E0_dp ) * 0.5E0_dp
  !
  !@@@ Simplified calculation of the smallest machine representable number
  !    (Higher accuracy than needed!)
  !
  !@@@ min = dpmpar(2)

  min = epsilon ( 1.0E0_dp)

  x = arg
  y = ABS ( x )

  IF ( y <= thrsh ) THEN
     !
     !  Evaluate  anorm  for  |X| <= 0.66291
     !
     IF ( y > eps ) THEN
        xsq = x * x
     ELSE
        xsq = 0.0_dp
     END IF

     xnum = a(5) * xsq
     xden = xsq
     DO i = 1, 3
        xnum = ( xnum + a(i) ) * xsq
        xden = ( xden + b(i) ) * xsq
     END DO
     RESULT = x * ( xnum + a(4) ) / ( xden + b(4) )
     temp = RESULT
     RESULT = 0.5_dp + temp
     ccum = 0.5_dp - temp
     !
     !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
     !
  ELSE IF ( y <= root32 ) THEN

     xnum = c(9) * y
     xden = y
!CDIR UNROLL=7
     DO i = 1, 7
        xnum = ( xnum + c(i) ) * y
        xden = ( xden + d(i) ) * y
     END DO
     RESULT = ( xnum + c(8) ) / ( xden + d(8) )
     xsq = AINT ( y * sixten ) / sixten
     del = ( y - xsq ) * ( y + xsq )
     RESULT = EXP(-xsq*xsq*0.5) * EXP(-del*0.5) * RESULT
     ccum = 1.0_dp - RESULT

     IF ( x > 0.0_dp  ) THEN
        temp = RESULT
        RESULT = ccum
        ccum = temp
     END IF
     !
     !  Evaluate  anorm  for |X| > sqrt(32).
     !
  ELSE

     RESULT = 0.0_dp
     xsq = 1.0 / ( x * x )
     xnum = p(6) * xsq
     xden = xsq
     DO i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
     END DO

     RESULT = xsq * ( xnum + p(5) ) / ( xden + q(5) )
     RESULT = ( sqrpi - RESULT ) / y
     xsq = AINT ( x * sixten ) / sixten
     del = ( x - xsq ) * ( x + xsq )
     RESULT = EXP ( - xsq * xsq * 0.5 ) * EXP ( - del * 0.5 ) * RESULT
     ccum = 1.0_dp - RESULT

     IF ( x > 0.0_dp ) THEN
        temp = RESULT
        RESULT = ccum
        ccum = temp
     END IF

   END IF

   IF ( RESULT < min ) THEN
     RESULT = 0.0_dp
   END IF
   
   IF ( ccum < min ) THEN
     ccum = 0.0_dp
   END IF
   
 END SUBROUTINE error_function_limited

END MODULE MESSY_CLOUD_DROPLET
