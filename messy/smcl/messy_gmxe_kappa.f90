MODULE MESSY_GMXE_KAPPA


  USE messy_main_constants_mem,     ONLY: dp, KAPPA_R => R_gas,STRLEN_MEDIUM

  IMPLICIT NONE

  PRIVATE :: dp
  SAVE

  INTEGER, PUBLIC  :: nss
  TYPE salt_info
    REAL(dp), DIMENSION(:,:,:), POINTER  :: mass    => NULL()
    REAL(dp)                             :: molmass = 0._dp
    REAL(dp)                             :: density = 0._dp
    REAL(dp)                             :: kappa   = 0._dp
    CHARACTER(LEN=STRLEN_MEDIUM)         :: name    = ''
    INTEGER                              :: iso_idx = 0      ! index in isorropia2
    INTEGER                              :: eqs_idx = 0      ! index in eqsam4clim
  END type salt_info
  TYPE(salt_info), DIMENSION(:), POINTER  :: salt   => NULL()
  

CONTAINS

!=========================================================================

  SUBROUTINE initialise_salts(kbdim, klev, nmod)

    USE messy_gmxe_mem,                ONLY: nanions, ncations, td
    USE messy_main_constants_mem,      ONLY: MH, MS, MO, MN, MCl, MBr, MI, &
                                             MC, MNa, rho_h2o

    INTEGER, INTENT(IN) :: kbdim, klev, nmod

    INTEGER  :: jc, jt, ji
    INTEGER  :: kcomp, jcomp
    INTEGER  :: counter

    INTEGER, PARAMETER :: ian  = 13
    INTEGER, PARAMETER :: icat = 8
    TYPE(salt_info)  :: salt_matrix(icat,ian)

    REAL(dp), PARAMETER :: MP = 30.97_dp

    CHARACTER(LEN=STRLEN_MEDIUM), PARAMETER :: cat_names(1:icat) =&
      (/'Hp   ', 'NH4p ', 'Nap  ', 'Kp   ', 'Capp ', 'Mgpp ',  &
        'Fepp ', 'Feppp'/)
    REAL(dp), PARAMETER                     :: cat_mass(1:icat) = &
      (/MH, 4._dp * MH + MN, MNa, 39.098_dp, 40.08_dp, 24.305_dp, &
      55.847_dp, 55.847_dp/)
    REAL(dp), PARAMETER                     :: cat_charge(1:icat)  = &
      (/1._dp,      1._dp,  1._dp,  1._dp,    2._dp,     2._dp,      &
        2._dp,   3._dp/)
      
    CHARACTER(LEN=STRLEN_MEDIUM), PARAMETER :: an_names(1:ian) = &
      (/'OHm    ','SO4mm  ', 'HSO4m  ', 'Clm    ', 'NO3m   ',    &
        'Brm    ','Im     ', 'CO3mm  ', 'HCO3m  ', 'PO4mmm ',    &
        'HCOOm  ','CH3COOm', 'C2O4mm '/)
    REAL(dp), PARAMETER                     :: an_mass(1:ian)  =    &
      (/MO+MH, MS + 4._dp*MO, MH+MS+4._dp*MO, MCl, MN+3_dp*MO,      &
        MBr,    MI,        MC+3._dp*MO, MH+MC+3._dp*MO, MP+4._dp*MO,&
        MH+MC+2._dp*MO,2._dp*MC+3._dp*MH+2._dp*MO, 2._dp*MC + 4._dp*MO/)
    REAL(dp), PARAMETER                     :: an_charge(1:ian)  =   &
      (/1._dp,  2._dp,      1._dp,          1._dp,    1._dp,         &
        1._dp,  1._dp,      2._dp,          1._dp,    3._dp,         &
        1._dp,  1._dp,      2._dp/)

    
    DO jc=1,icat
      DO jt=1,ian
        salt_matrix(jc,jt)%molmass = an_charge(jt)  * cat_mass(jc) + &
                                     cat_charge(jc) * an_mass(jt)
        salt_matrix(jc,jt)%kappa   = 0._dp
        salt_matrix(jc,jt)%density = 0._dp
        salt_matrix(jc,jt)%name    = TRIM(cat_names(jc))//'+'//TRIM(an_names(jt))
      END DO
    END DO
    
    ! H2O
    salt_matrix(1, 1)%density = rho_h2o
    ! H2SO4
    salt_matrix(1, 2)%density = 1.830_dp
    salt_matrix(1, 2)%kappa   = 0.9_dp
    salt_matrix(1, 3)%density = 1.830_dp
    ! HCl
    salt_matrix(1, 4)%density = 1.490_dp
    ! HNO3
    salt_matrix(1, 5)%density = 1.513_dp
    ! HBr
    salt_matrix(1, 6)%density = 3.307_dp
    ! HI
    salt_matrix(1, 7)%density = 5.228_dp
    ! H2CO3
    salt_matrix(1, 8)%density = 1.0_dp
    salt_matrix(1, 9)%density = 1.0_dp
    ! H3PO4
    salt_matrix(1,10)%density = 1.0_dp
    ! HCOOH
    salt_matrix(1,11)%density = 1.220_dp
    ! CH3COOH
    salt_matrix(1,12)%density = 1.045_dp
    ! H2C2O4
    salt_matrix(1,13)%density = 1.900_dp

    ! NH4OH
    salt_matrix(2, 1)%density = 1._dp
    ! (NH4)2SO4
    salt_matrix(2, 2)%density = 1.770_dp
    salt_matrix(2, 2)%kappa   = 0.53_dp
    ! NH4HSO4
    salt_matrix(2, 3)%density = 1.780_dp
    salt_matrix(2, 3)%kappa   = 0.53_dp
    ! NH4Cl
    salt_matrix(2, 4)%density = 1.519_dp
    salt_matrix(2, 4)%kappa   = 1.02_dp
    ! NH4NO3
    salt_matrix(2, 5)%density = 1.720_dp
    salt_matrix(2, 5)%kappa   = 0.67_dp
    ! NH4Br
    salt_matrix(2, 6)%density = 2.429_dp
    ! NH4I
    salt_matrix(2, 7)%density = 2.514_dp
    ! (NH4)2CO3
    salt_matrix(2, 8)%density = 1.0_dp
    salt_matrix(2, 5)%kappa   = 0.84_dp
    ! NH4HCO3
    salt_matrix(2, 9)%density = 1.586_dp
    salt_matrix(2, 5)%kappa   = 0.84_dp
    ! (NH4)3PO4
    salt_matrix(2,10)%density = 1.0_dp
    ! NH4HCOO
    salt_matrix(2,11)%density = 1.270_dp
    ! NH4CH3COO
    salt_matrix(2,12)%density = 1.073_dp
    ! (NH4)2C2O4
    salt_matrix(2,13)%density = 1.500_dp
    salt_matrix(2,13)%kappa   = 0.650_dp
    ! (NH4)3H(SO4)2
    salt_matrix(2,13)%density = 1.775_dp
    salt_matrix(2,13)%kappa   = 0.51_dp
    salt_matrix(2,13)%name    = 'LC'
    salt_matrix(2,13)%molmass = 247.3_dp

    ! NaOH
    salt_matrix(3, 1)%density = 2.13_dp
    ! Na2SO4
    salt_matrix(3, 2)%density = 2.70_dp
    salt_matrix(3, 2)%kappa   = 0.68_dp
    ! NaHSO4
    salt_matrix(3, 3)%density = 2.430_dp
    salt_matrix(3, 3)%kappa   = 1.01_dp
    ! NaCl
    salt_matrix(3, 4)%density = 2.17_dp
    salt_matrix(3, 4)%kappa   = 1.12_dp
    ! NaNO3
    salt_matrix(3, 5)%density = 2.260_dp
    salt_matrix(3, 5)%kappa   = 0.88_dp
    ! NaBr
    salt_matrix(3, 6)%density = 3.2_dp
    ! NaI
    salt_matrix(3, 7)%density = 3.67_dp
    ! Na2CO3
    salt_matrix(3, 8)%density = 2.54_dp
    salt_matrix(3, 5)%kappa   = 1.3_dp
    ! NaHCO3
    salt_matrix(3, 9)%density = 2.2_dp
    salt_matrix(3, 5)%kappa   = 1.3_dp
    ! Na3PO4
    salt_matrix(3,10)%density = 1.62_dp
    ! NaHCOO
    salt_matrix(3,11)%density = 1.920_dp
    ! NaCH3COO
    salt_matrix(3,12)%density = 1.528_dp
    ! Na2C2O4
    salt_matrix(3,13)%density = 3.610_dp
    salt_matrix(3,13)%kappa   = 0.940_dp

    ! KOH
    salt_matrix(4, 1)%density = 2.044_dp
    ! K2SO4
    salt_matrix(4, 2)%density = 2.66_dp
    salt_matrix(4, 2)%kappa   = 0.83_dp
    ! KHSO4
    salt_matrix(4, 3)%density = 2.320_dp
    salt_matrix(4, 3)%kappa   = 0.83_dp
    ! KCl
    salt_matrix(4, 4)%density = 1.988_dp
    salt_matrix(4, 4)%kappa   = 0.89_dp
    ! KNO3
    salt_matrix(4, 5)%density = 2.110_dp
    salt_matrix(4, 5)%kappa   = 0.88_dp
    ! KBr
    salt_matrix(4, 6)%density = 2.74_dp
    ! KI
    salt_matrix(4, 7)%density = 3.12_dp
    ! K2CO3
    salt_matrix(4, 8)%density = 2.29_dp
    salt_matrix(4, 5)%kappa   = 0.9_dp
    ! KHCO3
    salt_matrix(4, 9)%density = 2.17_dp
    salt_matrix(4, 5)%kappa   = 0.9_dp
    ! K3PO4
    salt_matrix(4,10)%density = 2.564_dp
    ! KHCOO
    salt_matrix(4,11)%density = 1.910_dp
    ! KCH3COO
    salt_matrix(4,12)%density = 1.570_dp
    ! K2C2O4
    salt_matrix(4,13)%density = 2.130_dp
    salt_matrix(4,13)%kappa   = 0.42_dp

    ! CaOH2
    salt_matrix(5, 1)%density = 2.200_dp
    ! CaSO4
    salt_matrix(5, 2)%density = 2.96_dp
    ! CaHSO42
    salt_matrix(5, 3)%density = 1.700_dp
    ! CaCl2
    salt_matrix(5, 4)%density = 2.158_dp
    ! CaNO32
    salt_matrix(5, 5)%density = 2.500_dp
    ! CaBr2
    salt_matrix(5, 6)%density = 3.38_dp
    ! CaI2
    salt_matrix(5, 7)%density = 3.96_dp
    ! CaCO3
    salt_matrix(5, 8)%density = 2.83_dp
    ! CaHCO32
    salt_matrix(5, 9)%density = 1.00_dp
    ! Ca3PO42
    salt_matrix(5,10)%density = 3.140_dp
    ! CaHCOO2
    salt_matrix(5,11)%density = 2.020_dp
    ! CaCH3COO2
    salt_matrix(5,12)%density = 1.50_dp
    ! CaC2O4
    salt_matrix(5,13)%density = 2.200_dp

    ! MgOH2
    salt_matrix(6, 1)%density = 2.370_dp
    ! MgSO4
    salt_matrix(6, 2)%density = 2.66_dp
    salt_matrix(6, 2)%kappa   = 0.8_dp
    ! MgHSO42
    salt_matrix(6, 3)%density = 1.00_dp
    salt_matrix(6, 3)%kappa   = 0.8_dp
    ! MgCl2
    salt_matrix(6, 4)%density = 2.325_dp
    salt_matrix(6, 4)%kappa   = 1.32_dp
    ! MgNO32
    salt_matrix(6, 5)%density = 2.300_dp
    salt_matrix(6, 5)%kappa   = 0.84_dp
    ! MgBr2
    salt_matrix(6, 6)%density = 3.72_dp
    ! MgI2
    salt_matrix(6, 7)%density = 4.43_dp
    ! MgCO3
    salt_matrix(6, 8)%density = 3.05_dp
    salt_matrix(6, 8)%kappa   = 1.3_dp
    ! MgHCO32
    salt_matrix(6, 9)%density = 1.00_dp
    salt_matrix(6, 9)%kappa   = 1.3_dp
    ! Mg3PO42
    salt_matrix(6,10)%density = 1.000_dp
    ! MgHCOO2
    salt_matrix(6,11)%density = 1.00_dp
    ! MgCH3COO2
    salt_matrix(6,12)%density = 1.50_dp
    ! MgC2O4
    salt_matrix(6,13)%density = 1.00_dp
    salt_matrix(6,13)%kappa   = 0.79_dp

    ! FeOH2
    salt_matrix(7, 1)%density = 3.120_dp
    ! FeSO4
    salt_matrix(7, 2)%density = 3.10_dp
    ! FeHSO42
    salt_matrix(7, 3)%density = 1.00_dp
    ! FeCl2
    salt_matrix(7, 4)%density = 2.900_dp
    ! FeNO32
    salt_matrix(7, 5)%density = 1.00_dp
    ! FeBr2
    salt_matrix(7, 6)%density = 4.5_dp
    ! FeI2
    salt_matrix(7, 7)%density = 1._dp
    ! FeCO3
    salt_matrix(7, 8)%density = 3.9_dp
    ! FeHCO32
    salt_matrix(7, 9)%density = 1.00_dp
    ! Fe3PO42
    salt_matrix(7,10)%density = 2.870_dp
    ! FeHCOO2
    salt_matrix(7,11)%density = 1.00_dp
    ! FeCH3COO2
    salt_matrix(7,12)%density = 1.0_dp
    ! FeC2O4
    salt_matrix(7,13)%density = 1.00_dp

    salt_matrix(8,1:ian)%density = salt_matrix(7,1:ian)%density 

    nss = 0
    DO jc=1,nanions
      jcomp = MAXVAL(td%anion(jc)%aerml_idx(:))
      IF (jcomp == 0) CYCLE
      DO jt=1,ncations
        kcomp = MAXVAL(td%cation(jt)%aerml_idx(:))
        IF (kcomp == 0) CYCLE
        nss = nss + 1
      END DO
    END DO
    ! add one for Letovicite if NH4 and SO4mm and HSO4m exist
    do jt=1,ncations
       IF (td%cation(jt)%name == 'NH4p') THEN
          DO jc=1,nanions
             IF (td%anion(jc)%name == 'SO4mm') THEN
                do ji = 1,nanions
                   IF (td%anion(jc)%name == 'HSO4m') THEN
                      nss = nss + 1
                   END IF
                END do
             END IF
          END DO
       END IF
    END do
    
    ALLOCATE(salt(nss))

    counter = 0
    DO jc=1,nanions
      jcomp = MAXVAL(td%anion(jc)%aerml_idx(:))
      IF (jcomp == 0) CYCLE
      DO jt=1,ncations
        kcomp = MAXVAL(td%cation(jt)%aerml_idx(:))
        IF (kcomp == 0) CYCLE
        counter = counter + 1
        salt(counter)%name=TRIM(td%cation(jt)%name)//'+'//TRIM(td%anion(jc)%name)
      END DO
    END DO
    ! add one for Letovicite if NH4 and SO4mm and HSO4m exist
    do jt=1,ncations
       IF (td%cation(jt)%name == 'NH4p') THEN
          DO jc=1,nanions
             IF (td%anion(jc)%name == 'SO4mm') THEN
                do ji = 1,nanions
                   IF (td%anion(jc)%name == 'HSO4m') THEN
                      counter = counter + 1
                      salt(counter)%name='LC'
                   END IF
                END do
             END IF
          END DO
       END IF
    END do

    DO ji = 1, nss
      salt(ji)%molmass = 1._dp
      salt(ji)%density = 1._dp
      salt(ji)%kappa   = 0._dp
      ALLOCATE(salt(ji)%mass(kbdim,klev,nmod))
      salt(ji)%mass(:,:,:) = 0._dp
      DO jc=1,icat
        DO jt=1,ian
          IF (TRIM(salt(ji)%name) == TRIM(salt_matrix(jc,jt)%name) ) THEN
            salt(ji)%molmass = salt_matrix(jc,jt)%molmass
            salt(ji)%density = salt_matrix(jc,jt)%density
            salt(ji)%kappa   = salt_matrix(jc,jt)%kappa
          END IF
        END DO
      END DO
    END DO

      

  END SUBROUTINE initialise_salts
  
!=========================================================================

 
  
  SUBROUTINE calc_kappa(kproma, klev, jm, naertot, &
                        paerosol, paerml, paernl,  &
                        RH, TT, prdry)

    USE messy_gmxe_mem,       ONLY: bulk, nbulk, nsoluble, ndiff, &
                                    sigma

    INTEGER  :: kproma, klev, naertot
    REAL(dp) :: paerosol(kproma,klev,8)
    REAL(dp) :: paerml(kproma,klev,0:naertot), paernl(kproma,klev)
    REAL(dp) :: RH(kproma,klev), prdry(kproma,klev), TT(kproma,klev)

    REAL(dp) :: zdrybulkv_species(kproma,klev,nss+nbulk)
    REAL(dp) :: zdrybulkv_species_insol(kproma,klev,nss+nbulk)
    REAL(dp) :: zdrybulkv_sum(kproma,klev)
    REAL(dp) :: zdrybulkv_sum_insol(kproma,klev)
    REAL(dp) :: kappa(kproma,klev)
    REAL(dp) :: kappa_insol(kproma,klev)
    REAL(dp) :: kappa_water(kproma,klev)
    REAL(dp) :: CCN_mode_2(kproma,klev)
    REAL(dp) :: CCN_mode_4(kproma,klev)
    REAL(dp) :: ka_vol(kproma,klev,nss+nbulk)
    REAL(dp) :: ka_vol_sum(kproma,klev)
    REAL(dp) :: vol_frac(nss+nbulk)
    REAL(dp) :: vol_frac_insol(nss+nbulk)
    
    REAL(dp) :: Kappa_Mv     = 0.018014999
    REAL(dp) :: Kappa_sigma  = 0.072
    REAL(dp) :: Kappa_rhow   = 997.0
    REAL(dp) :: Kappa_a
    REAL(dp) :: Ds, Dt,Sold, Snew, top, bot, aw 

    ! op_pj_20131105+
    REAL(dp) :: A, B, C, t32, t3, t, x, diffT, Dt2, Dt3
    ! op_pj_20131105-

    REAL(dp) :: crcut_2, crcut_4, ztn_2, ztn_4, lfrac_2, lfrac_4, zdummy

    INTEGER,PARAMETER :: nref_kappa = 17
    INTEGER  :: ll

    INTEGER  :: jc, idx1, kcomp ! species
    INTEGER  :: jm, jl, jk, km
    REAL(dp), PARAMETER ::cmin_epsilon   = 1.e-18_dp

    REAL(dp), PARAMETER :: ref_critical_diameter_2(1:nref_kappa) = (/    &
         583.389954_dp, 496.259979_dp, 413.039978_dp, 367.659973_dp, 337.499969_dp, &
         315.399994_dp, 253.959991_dp, 203.269989_dp, 178.019989_dp, 161.979996_dp, &
         150.399994_dp, 119.649994_dp, 94.8999939_dp, 83.0499954_dp, 75.4599991_dp, &
         70.0299988_dp, 55.5999985_dp                                      /)
    
    REAL(dp), PARAMETER :: ref_critical_diameter_4(1:nref_kappa) = (/    &
         334.699982_dp, 292.000000_dp, 248.459991_dp, 223.659988_dp, 206.739990_dp, &
         194.079987_dp, 157.929993_dp, 127.149994_dp, 111.619995_dp, 101.679993_dp, &
         94.5399933_dp, 75.2799988_dp, 59.8499985_dp, 52.3099976_dp, 47.5399971_dp, &
         44.1399994_dp, 35.0499992_dp                                       /)
    
    REAL(dp), PARAMETER :: ref_kappa(1:nref_kappa) = (/                  &
         0.0010_dp,0.0020_dp,0.0040_dp,0.0060_dp,0.0080_dp,0.0100_dp,0.0200_dp,0.0400_dp,    &
         0.0600_dp,0.0800_dp,0.1000_dp,0.2000_dp,0.4000_dp,0.6000_dp,0.8000_dp,1.0000_dp,    & 
         2.0000_dp                                                      /)


    zdrybulkv_sum(:,:)       = 0.0_dp
    zdrybulkv_species(:,:,:) = 0.0_dp
    kappa(:,:)               = 0.0_dp
    kappa_insol(:,:)         = 0.0_dp
    ka_vol_sum(:,:)          = 0.0_dp
    zdrybulkv_sum_insol(:,:) = 0.0_dp
    zdrybulkv_sum(:,:)       = 0.0_dp
    paerosol(:,:,:)          = 0.0_dp

    DO jk=1,klev
      DO jl=1,kproma
        
        vol_frac(:)    = 0.0_dp
        vol_frac_insol(:)    = 0.0_dp
        
        ! Loop over species treated in the thermodynamics
        ! Currently consider both aqueous and crystaline concentrations to calculate Kappa (na + nc)
        IF ( jm <= nsoluble ) THEN  
          
          DO jc=1,nss

            ! Total dry volume of ThD species in hydrophillic mode
            zdrybulkv_species (jl,jk,jc) =        &
              salt(jc)%mass(jl,jk,jm) &  
              * 1.e-6_dp * &
              salt(jc)%molmass  / &
              salt(jc)%density 
            
            IF(jm > 1)THEN
              ! Add mass of species condensed on hydrophobic aerosol to the kappa of 
              ! corresponding hydrophillic mode
              km = jm + nsoluble - ndiff
              zdrybulkv_species (jl,jk,jc) =  zdrybulkv_species (jl,jk,jc) +  &  
                salt(jc)%mass(jl,jk,km)  &
                * 1.e-6_dp * &
                salt(jc)%molmass  / &
                salt(jc)%density 
            ENDIF

            zdrybulkv_sum(jl,jk) = zdrybulkv_sum(jl,jk) + zdrybulkv_species (jl,jk,jc)

          ENDDO
          
          ! Add volume due to bulk species in soluble mode
          DO jc = 1,nbulk

            idx1 = nss + jc
            kcomp = bulk(jc)%aerml_idx(jm)
            
            IF ( kcomp == 0 ) CYCLE    
            IF ( paerml(jl,jk,kcomp)  >= cmin_epsilon ) THEN  
              
              zdrybulkv_species (jl,jk,idx1)      = &
                paerml(jl,jk,kcomp) * 1.e-6_dp    * &
                bulk(jc)%molmass     / &
                bulk(jc)%density

              zdrybulkv_sum(jl,jk) = zdrybulkv_sum(jl,jk) + &
                zdrybulkv_species (jl,jk,idx1)

            ENDIF
          ENDDO

!!$          ! Add volume due to bulk species in insoluble mode
          IF (jm > 1) THEN
            DO jc = 1,nbulk
              km = jm + nsoluble - ndiff
              idx1 = nss + jc
              kcomp = bulk(jc)%aerml_idx(km)  
              IF ( kcomp == 0 ) CYCLE  
              IF ( paerml(jl,jk,kcomp)  >= cmin_epsilon ) THEN 
                zdrybulkv_species_insol (jl,jk,idx1) = &
                  paerml(jl,jk,kcomp) * 1.e-6_dp       * &
                  bulk(jc)%molmass        / &
                  bulk(jc)%density
                
                zdrybulkv_sum_insol(jl,jk) =         &
                  zdrybulkv_sum_insol(jl,jk)       + &
                  zdrybulkv_species_insol(jl,jk,idx1)
              END IF
            END DO
          END IF
             
          !!Sum volume of insoluble and soluble modes
          zdrybulkv_sum_insol(jl,jk) = zdrybulkv_sum_insol(jl,jk) + &
                                       zdrybulkv_sum(jl,jk)


          DO jc=1,nss+nbulk

            IF ( zdrybulkv_sum(jl,jk) > 0.0_dp )         &
              vol_frac(jc) = zdrybulkv_species(jl,jk,jc) &
                           / zdrybulkv_sum(jl,jk)

            IF ( zdrybulkv_sum_insol(jl,jk) > 0.0_dp )          &
              vol_frac_insol(jc) = zdrybulkv_species (jl,jk,jc) &
                                 / zdrybulkv_sum_insol(jl,jk)

            ! Calculate Kappa value (and absolute volume weighted kappa ka_vol)
            IF ( jc <= nss ) THEN
                 
              kappa(jl,jk) = kappa(jl,jk) +            &
                ( vol_frac(jc) * salt(jc)%kappa )
              
              kappa_insol(jl,jk)= kappa_insol(jl,jk) + &
                ( vol_frac_insol(jc) * salt(jc)%kappa )
              ka_vol(jl,jk,jc) = salt(jc)%kappa * zdrybulkv_species (jl,jk,jc)
            ELSE
              idx1 = jc - nss
              kappa(jl,jk) = kappa(jl,jk) +             &
                ( vol_frac(jc) * bulk(idx1)%kappa )
              
              kappa_insol(jl,jk) = kappa_insol(jl,jk) + &
                ( vol_frac_insol(jc) * bulk(idx1)%kappa )
              ka_vol(jl,jk,jc) = bulk(idx1)%kappa &
                * zdrybulkv_species (jl,jk,jc)
            ENDIF
            ka_vol_sum(jl,jk) = ka_vol_sum(jl,jk) + ka_vol(jl,jk,jc)
          ENDDO


            ! Volume of water on aerosol (Vw)  P&K (2007) Eq 3.
          IF(RH(jl,jk)  > 0.0_dp)THEN
            kappa_water(jl,jk) = (( RH(jl,jk)           &
                                / (1.0_dp - RH(jl,jk)) ) &
                                * ka_vol_sum(jl,jk))
          ENDIF

            ! Kappa value passed out
          paerosol(jl,jk,1) = kappa(jl,jk)      
          paerosol(jl,jk,2) = kappa_insol(jl,jk)   

          ! Kappa volume value passed out
          paerosol(jl,jk,3) = zdrybulkv_sum(jl,jk)       
          paerosol(jl,jk,4) = zdrybulkv_sum_insol(jl,jk)   

          ! Total aerosol water volume passed out 
          paerosol(jl,jk,5) = kappa_water(jl,jk)  
          
        ENDIF  ! jm < nsoluble
      END DO
    END DO

!!  Using Calculated Kappa and aerosol dry diameter calculate the critcal supersaturation
!!  Adapted from code provided my M. D. Petters (findsc.pro)
!!  See Eq. 6 Petters and Kreidenweis (2007) 

      DO jk=1,klev
        DO jl=1,kproma

          paerosol(jl,jk,6)        = 0.0_dp
          paerosol(jl,jk,7)      = 0.0_dp
          paerosol(jl,jk,8)      = 0.0_dp
          
          CCN_mode_2(jl,jk)          = 0.0_dp
          CCN_mode_4(jl,jk)          = 0.0_dp
  
          IF ( jm <= nsoluble ) THEN
            IF(zdrybulkv_sum(jl,jk) > 0.0_dp .AND. paernl(jl,jk) > 0.0_dp & 
              .AND. prdry(jl,jk) > 0.0_dp) THEN

              kappa_A = 4.0_dp * kappa_Mv* kappa_sigma &
                      / (kappa_R*TT(jl,jk)*kappa_rhow)

              !! Convert radius to diameter 
              Ds  = prdry(jl,jk) * 1.e-2_dp * 2.0_dp  

! op_pj_20131105+
#ifdef _XOLD
              Dt   = Ds 
              Sold = 0.0
              Snew = 0.1
              
              do while (Snew < 0.9)
                Dt = Dt * 1.1
                top = (Dt**3.0 - Ds**3.0) 
                bot = Dt**3 - Ds**3.0 * (1.0-kappa(jl,jk))  
                aw = top/bot 
                Sold = Snew
                Snew = aw*exp(kappa_A/Dt)
              enddo

              do while (Snew > Sold)
                Dt = Dt * 1.0005
                top = (Dt**3.0 - Ds**3.0) 
                bot = Dt**3 - Ds**3.0 * (1.0-kappa(jl,jk))
                aw = top/bot  
                Sold = Snew 
                Snew = aw*exp(kappa_A/Dt)
              enddo
#else
              A = Ds*Ds*Ds
              B = A*(1.0-kappa(jl,jk))
              t32 = kappa(jl,jk)*A/kappa_A
              C = -3.0_dp*t32
              t3 = sqrt(t32)
              !set Dt to the starting point
              Dt=t3+(A+B-t3*(C+t32))**0.33333333333333333_dp
              x=Dt
              do while (abs(x/Dt)>1.0E-15_dp)
                 Dt2 = Dt*Dt
                 Dt3 = Dt2*Dt
                 T = C*Dt2*Dt2+(Dt3-A)*(Dt3-B)
                 diffT= 4*C*Dt3+3*Dt2*(2*Dt3-A-B)
                 x=T/diffT
                 Dt=Dt-x
              enddo
              Snew = (Dt3-A)/(Dt3-B)*exp(kappa_A/Dt)
#endif
! op_pj_20131105-

              if( kappa(jl,jk) == 0) Snew = exp(kappa_A/Ds)
              
              Snew = (Snew-1)*100

              paerosol(jl,jk,6)  = Snew  !! Output critical supersaturation

            ENDIF
          ENDIF
        ENDDO
      ENDDO

!! Using calculated Sc, calculate the fraction of the mode that is activated 
!! at 0.2% supersaturation (calculated offline from P&K Eq 6)

      IF ( jm <= nsoluble ) THEN
        DO jk=1,klev
          DO jl=1,kproma
            IF (prdry(jl,jk) < 1.e-10_dp) CYCLE
      !! Use linear interpolation to calulate critical radius from Sc and Kappa.
           
              !! CCN at 0.2% and 0.4%
            IF ( kappa(jl,jk) .LE. ref_kappa(1) ) THEN
              crcut_2 = ref_critical_diameter_2(1) / 1.0E9_dp
              crcut_4 = ref_critical_diameter_4(1) / 1.0E9_dp
            ELSEIF ( kappa(jl,jk) .GE. ref_kappa(nref_kappa) ) THEN
              crcut_2 = ref_critical_diameter_2(nref_kappa) / 1.0E9_dp  
              crcut_4 = ref_critical_diameter_4(nref_kappa) / 1.0E9_dp  
            ELSE
              DO ll = 1,nref_kappa-1
                IF ( kappa(jl,jk) .GE. ref_kappa(ll) .AND. &
                  kappa(jl,jk) .LT. ref_kappa(ll+1) ) THEN
                    
                  crcut_2 = ref_critical_diameter_2(ll)            &
                          + (kappa(jl,jk) - ref_kappa(ll))         &
                          * ( (ref_critical_diameter_2(ll+1)       &
                          - ref_critical_diameter_2(ll))           &
                          / (ref_kappa(ll+1) - ref_kappa(ll) )) 
                  crcut_4 = ref_critical_diameter_4(ll)            &
                          + (kappa(jl,jk) - ref_kappa(ll))         & 
                          * ( (ref_critical_diameter_4(ll+1)       &
                          - ref_critical_diameter_4(ll))           &
                          / (ref_kappa(ll+1) - ref_kappa(ll) )) 
                          
                  crcut_2 = crcut_2 / 1E9_dp  !! Convert from nm to m
                  crcut_4 = crcut_4 / 1E9_dp  
                  
                ENDIF
              ENDDO
            ENDIF
                 
              !! Calculate fractional activation of mode

            ztn_2 = ( LOG(crcut_2) - LOG(prdry(jl,jk) * 2.0_dp * 1.e-2_dp) ) &
                    / LOG(sigma(jm))  
            ztn_4 = ( LOG(crcut_4) - LOG(prdry(jl,jk) * 2.0_dp * 1.e-2_dp) ) &
                    / LOG(sigma(jm))  
                 
     !--- Calculate the cumulative of the log-normal number distribution:
                 
            CALL error_function_limited(ztn_2,zdummy,lfrac_2)
                 
            CCN_mode_2(jl,jk) = lfrac_2 * paernl(jl,jk)
                 
            CALL error_function_limited(ztn_4,zdummy,lfrac_4)
                 
            CCN_mode_4(jl,jk) = lfrac_4 * paernl(jl,jk)
                 
            paerosol(jl,jk,7) = CCN_mode_2(jl,jk)
                 
            paerosol(jl,jk,8) = CCN_mode_4(jl,jk)

          END DO
        ENDDO
      ENDIF


  END SUBROUTINE calc_kappa

 !--------------------------------------------------------------------

  SUBROUTINE CALC_wateruptake_bulk(kproma, klev, jm, naertot, &
                                   paerml, zwatbulk, WH2O, prhum, &
                                   water_oc)
 
    USE messy_gmxe_mem,              ONLY: bulk, nbulk, nmod
    USE messy_main_constants_mem,    ONLY: mwh2o => M_H2O, rho_h2o, m_air
   
    INTEGER  :: kproma, klev, jm, naertot

    INTEGER  :: jc, jl, jk, kcomp, jcomp

    REAL(dp) :: V_w, V_s, M_wat, Dw, zhelp
    REAL(dp) :: wh2o(kproma,klev,nmod), prhum(kproma,klev), prho(kproma,klev)
    REAL(dp) :: zwatbulk(kproma,klev), water_oc(kproma,klev)
    REAL(dp) :: paerml(kproma,klev,0:naertot)
    REAL(dp), PARAMETER :: RHMAX = 0.95

    Dw = rho_h2o * 1.e-3_dp
    DO jc = 1,nbulk
      IF ( bulk(jc)%kappa < 1.e-3_dp ) CYCLE
      kcomp = bulk(jc)%aerml_idx(jm)

       !    based on kappa approach
       !    1 / a_w = 1 + k * V_s / V_w
       !    a_w = RHUM
      IF (kcomp /= 0) THEN
        DO jk=1,klev
          DO jl=1,kproma
            IF (PRHUM(jl,jk) > 0.35_dp) THEN
              ! electrolyte volume [cm^3 / m^3]
              ! ?mol / m^3 * 1e-6 * g/mol / g/cm^3 = 
              ! mol / m^3 * g / mol / g /cm^3     = cm^3 / m^3
              V_s = MAX(0._dp,paerml(jl,jk,kcomp)) * 1.e-6_dp  &
                * bulk(jc)%molmass &
                / bulk(jc)%density
              ! water volume [cm^3 / m^3]
              V_w = bulk(jc)%kappa * V_s / &
                ( 1._dp / MIN(prhum(jl,jk),RHMAX) - 1._dp)
              ! water mass [g/m^3]
              M_wat = V_w * Dw
              ! water concentration [?mol / m^3]
              zhelp = M_wat / mwh2o * 1.e6_dp
              WH2O(jl,jk,jm) = WH2O(jl,jk,jm) + zhelp
              zwatbulk(jl,jk) = zwatbulk(jl,jk) + zhelp
              IF ( bulk(jc)%l_oc) &
                   water_oc(jl,jk) = water_oc(jl,jk) + zhelp * 1.e-9_dp * mwh2o
              ! storing wH2O for diagnostic output
              WH2O(jl,jk,jm)  = zwatbulk(jl,jk) * mwh2o
            ENDIF
          END DO
        END DO
      END IF
      
    END DO
    
  END SUBROUTINE CALC_wateruptake_bulk

!-----------------------------------------------------------------------------

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


!-----------------------------------------------------------------------------
END MODULE MESSY_GMXE_KAPPA
