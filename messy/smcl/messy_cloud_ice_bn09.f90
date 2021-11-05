
 MODULE messy_cloud_ice_BN09

!=======================================================================
!
! *** Code Developer: Donifan Barahona 
!                     donifan.o.barahona@nasa.gov
!
!                     DESCRIPTION
!                     Parameterization  of ICE crystal number concentration 
!                     for large scale models. 
!
!                     Donifan Barahona, Athanasios Nenes
!                     JGR, 111, D11211,  2008
!                     ACP, 9, 369-381,   2009 
!                     ACP  9, 5933-5948, 2009

!                     Homogeneoeus and heterogeneous nucleation considered   
!                     SI units unless otherwise specified. 
! 
! *** Implemented by Sara Bacer (MPIC) 
!
!=======================================================================

      USE messy_main_constants_mem,   ONLY: dp,                  &   !=SELECTED_REAL_KIND(12,307)
                                            rgas_ice => R_gas,   &   !Universal gas constant [J/K/mol]
                                            grav_ice => g,       &   !Gravity acceleration [m/s2] 
                                            cpa_ice  => cp_air,  &   !Specific heat of dry air [J/K/kg]
                                            pi_ice   => pi,      &   !PI constant
                                            M_H2O,               &   !Molar mass of H2O =18.02 [g/mol]  
                                            M_air                    !Molar mass of dry air =28.970 [g/mol]

#if defined (LF) || defined (NAGFOR) || defined (__PGI)
      USE messy_main_tools, ONLY: isnan
#endif

      implicit none
      private
      
      public :: ice_activate   
     
!----------- VARIABLES USED IN THE ICE NUCLEATION PARAMETERIZATION -------

      real(dp) :: g1_ice, g2_ice, alfa_ice, beta_ice, z_ice, lambda_ice,       &
                  denwat_ice, denice_ice, dhs_ice, vpresw_ice, vpresi_ice,     &
                  diff_ice, aircond_ice, shom_ice, koft_ice,                   &
                  gdoin_ice, nin_ice, si0dust_ice, si0bc_ice, dH1smooth,       &
                  del1dust_ice, del1bc_ice,Pice, denair_ice

      real(dp) :: wmw_ice, amw_ice,                                            &   !General Parameters
                  depcoef_ice, thaccom_ice, zero_par, rv_ice,                  &   !General Parameters 
                  doin_ice, sh_ice,                                            &
                  shdust_ice, effdust_ice, effbc_ice, shbc_ice,                &
                  kdust_ice, kbc_ice, acorr_dust, acorr_bc,                    &
                  norg_ice, sigorg_ice, dorg_ice,                              &
                  nbc_ice, dbc_ice, sigbc_ice,                                 &
                  nbio_ice, sigbio_ice, dbio_ice,                              &
                  np_ice, dliq_ice, ddry_ice,                                  &  
                  normv_ice, vmin_ice, vmax_ice, miuv_ice, sigmav_ice

      real(dp), allocatable, dimension(:) :: ndust_ice, sigdust_ice, ddust_ice, &  
                                             nsolo_ice, sigsolo_ice, dsolo_ice   
      integer              :: typeofspec_ice, nbindust_ice, nbinsolo_ice
      logical              :: purehet_ice, purehom_ice
      REAL(dp), PARAMETER  :: zeps   = EPSILON(1.0_dp)

      !DATA Statements for ice 
      DATA depcoef_ice /0.1d0/  ! Default deposition coefficient
      DATA thaccom_ice /0.7d0/  ! Default thermal accommodation coefficient
      DATA zero_par /1.d-12/    !mz_sb_20170203: instead of 1e-12

!------ end of declarations ---------------------------------------------------

 
      CONTAINS

 
!===================================================================================
! SUBROUTNE ICE_ACTIVATE
! It sets the variables needed for 
! the activation subroutines and returns the ice number concentration
!
! tparc     = T (K)
! pparc     = P (Pa)
! sigw      = width of the distribution of updraft velocity (m s-1) between 0.01 and 1
! zicncq    = preexisting ice crystal concentration (m-3)  
! nice      = nucleated ice number concentration (m-3)
! nhetice   = ice crystal concentration by het freezing (m-3)
! smaxice   = maximum supersaturation with respect to ice
! nlim      = limiting IN concentration (m-3)
! scice     = characteristic freezing point of the aerosol population (output)   
! ndust     = dust number concentration per mode (m-3)
! ddust     = geometric mean diameter of dust per mode
! norg      = organics number concentration (m-3)
! dorg      = geometric mean diameter of the organics (m)
! nsoot     = soot number concentration (m-3)
! dsoot     = geometric mean diameter of soot (m)
! nbio      = number concentration of biological particles (m-3)
! dbio      = geometric mean diameter of biological particles (m)
! nsolo     = number concentration of soluble organics (m-3)
! dsolo     = geometric mean diameter of soluble organics (m)
! ndrop     = total (=sum over modes) cloud droplet number concentration (m-3)    
! dsulf     = dry diameter of sulfate (m)
!===================================================================================
 
 SUBROUTINE ice_activate(tparc, pparc, sigw, zicncq,                         &   
                         ndrop, dsulf,       ndust, ddust, sigdust,          &
                         norg, dorg, sigorg, nsoot, dsoot, sigsoot,          &
                         nbio, dbio, sigbio, nsolo, dsolo, sigsolo,          &
                         smaxice, nlim, scice, nhetice, nice, DC, sigwparc)

      real(dp), intent(in)   :: tparc, pparc, sigw, zicncq,                   &
                                ndrop, dsulf,          norg, dorg, sigorg,    & 
                                nsoot, dsoot, sigsoot, nbio, dbio, sigbio
      real(dp), dimension(:), intent(in) :: ddust, ndust, sigdust,            &
                                            dsolo, nsolo, sigsolo
      real(dp), intent(out)  :: smaxice, nhetice, nlim, nice, scice, DC, sigwparc

      real(dp)  :: antot, sigw_preice, T, P
      logical   :: preice_l
  

      !Initialize output    
      scice    = 2.0         !default value if no aerosol present      
      DC       = 1.0d-6*2d0
      
      !***********Get aerosol parameters**************  

      !----> Hom freezing   
      np_ice   = ndrop
      ddry_ice = dsulf
    
      !----> Het freezing (DU, BC, OC, BIO)
      acorr_dust = 2.7e7  !m2/m3 correction to the area due to non sphericity and aggregation Assumes 10 g/m2 (Murray 2011), only PDA08
      acorr_bc   = 4.5e7  !m2/m3 correction to the area due to non sphericity and aggregation Assumes 50 g/m2 (Murray 2011), only PDA08
     
      !Dust 
      nbindust_ice = size(ndust)       !SB: nbindust_ice = dimensions of array ndust
      allocate(ndust_ice(nbindust_ice))
      allocate(sigdust_ice(nbindust_ice))
      allocate(ddust_ice(nbindust_ice))
      ndust_ice   = ndust
      ddust_ice   = ddust
      sigdust_ice = sigdust
    
      !Black carbon
      nbc_ice   = nsoot   
      dbc_ice   = dsoot
      sigbc_ice = sigsoot

      !Insoluble organics (used in PDA08 parameterization) 
      norg_ice   = norg      
      dorg_ice   = dorg
      sigorg_ice = sigorg

      !Biological aerosols (used in PDA13 parameterization)
      nbio_ice   = nbio
      dbio_ice   = dbio
      sigbio_ice = sigbio     

      !Soluble organics (used in PDA13 parameterization)
      nbinsolo_ice = size(nsolo)
      allocate(nsolo_ice(nbinsolo_ice))
      allocate(sigsolo_ice(nbinsolo_ice))
      allocate(dsolo_ice(nbinsolo_ice))
      nsolo_ice   = nsolo
      dsolo_ice   = dsolo
      sigsolo_ice = sigsolo

      antot   = sum(ndust_ice) + norg_ice + nbc_ice + nbio_ice + sum(nsolo_ice) + np_ice     !SB: total aerosol number conc, nsolo addition

  
      !***********Calculate nucleated crystal number**************

      !Heteroegeneous spectrum
      typeofspec_ice = 5                !SB: used just in SUBROUTINE INSPEC_ice (for different CASE)
                                        !-1 - Monodisperse                       
                                        ! 1 - Meyers et al. 1992
                                        ! 2 - Phillips et al. 2007, PDG07
                                        ! 3 - Barahona & Nenes 2009, CNT: assumme a maximum freezing fraction then scale it 
                                        ! 4 - Phillips et al. 2008, PDA08
                                        ! 5 - Phillips et al. 2013, PDA13
      
      purehet_ice = .FALSE.   ! logical, True supresses homog nucleation       
      purehom_ice = .FALSE.   ! logical, True supresses het nucleation   

      !Variables used only in prop_ice to compute PREICE effect:
      !--- preexisting ice crystals concentration   
      preice_l = .TRUE.
      sigw_preice = zero_par 

      !Other conditions limit to avoid errors outside the troposphere     
      P    = max(pparc, 1.0)
      T    = max(tparc, 190.0) !lower limit for tropospheric cirrus
  
      CALL prop_ice(T, P, zicncq, sigw, sigw_preice)     !Input: T, P, zicncq, sigw   Output: sigw_preice
  
      !--- Subgrid vertical velocity
      sigwparc = max(0.001, sigw)               
      if (preice_l)   sigwparc = max(0.001, sigw_preice)
      sigwparc = min(3.0, sigwparc)            !SB: 0.001 < sigwparc < 3.0 
      
      !----------------------------------------------------------
      !--- Distintion between regimes:
      !----------------------------------------------------------
      !
      !----------------------------------------------------------
      !    Mixed-phase regime ( 238.15 < T < 273.15)
      !----------------------------------------------------------
      ! 
      if (tparc .lt. 273.0d0)  then
         if (tparc .gt. 238.15d0) then     !mz_sb_20170202: changed from 235 

            if (antot .gt. 0.1) then       !only if aerosol is present => het freezing
               CALL IceParam (sigwparc, tparc, pparc, nhetice, nice, smaxice, nlim, DC, scice)    !mz_sb_20170112: addition DC,scice 
            end if
      ! 
      !----------------------------------------------------------
      !    Cirrus regime ( T <= 238.15 )
      !----------------------------------------------------------
      ! 
         else
    
            if (antot .gt. 0.1) then    !only if aerosol is present
               CALL IceParam (sigwparc, tparc, pparc, nhetice, nice, smaxice, nlim, DC, scice)      !mz_sb_20170112: addition DC,scice
            end if
       
         end if  !tparc
      end if  !tparc
      !
      !----------------------------------------------------------

      smaxice = max(smaxice, zero_par)    
      smaxice = min(smaxice, 2.0)        

      nhetice = max(nhetice, zero_par)
      nhetice = min(nhetice, 1e10)

      nlim    = max(nlim, zero_par)
      nlim    = min(nlim, 1e10)    

      nice    = max(nice, zero_par)

      scice   = max(scice, 1.0)
      scice   = min(scice, 2.0)
    
      deallocate(ndust_ice)
      deallocate(ddust_ice)
      deallocate(sigdust_ice)
      deallocate(nsolo_ice)
      deallocate(dsolo_ice)
      deallocate(sigsolo_ice)

      return 
      
 END SUBROUTINE ice_activate


!************************************************************************
!
! ICE FREEZING PARAMETERIZATION FILES START HERE
!
!************************************************************************

!===================================================================================
! SUBROUTINE IceParam 
!===================================================================================
   
 SUBROUTINE IceParam (sigma_w, T, P, nhet, nice, smax, nlim, DC, scice)
      
      real(dp), intent(in)  :: T, P, sigma_w 
      real(dp), intent(out) :: nice, nhet, smax, nlim, DC, scice

      real(dp)  :: wpar_ice, T_ice, P_ice
      logical :: use_av_v 

      !Other conditions limit to avoid errors outside the troposphere
      P_ice       = max(P, 1.0)
      T_ice       = max(T, 190.0) !lower limit for tropospheric cirrus
      depcoef_ice = 0.70          !deposition coefficient                  

      !Updraft Velocity Distribution
      use_av_v    = .FALSE.       !set to false to integrate over a pdf of updraft
      
      if (use_av_v) then
         wpar_ice = max(0.01, sigma_w*0.8) !m/s !characteristic velocity (minimum 1 cm/s) 

         CALL nice_param(wpar_ice, T_ice, P_ice, nice, smax, nhet, nlim, DC, scice)
      else
         vmin_ice   = 0.01d0   !some default values
         vmax_ice   = .5d0 
         sigmav_ice = sigma_w  !standard deviation of the V distribution m/s
         miuv_ice   = 0.001d0  !mean velocity for the V distribution (can be set to the large scale updraft)
         
         vmax_ice = min(miuv_ice+(4d0*sigmav_ice), 3.0)    !upper limit for integration
         vmin_ice = max(miuv_ice-(4d0*sigmav_ice), 0.01) 
            
         CALL nice_Vdist(T_ice, P_ice, nice, smax, nhet, nlim, DC, scice)
         !To parameterize updraft subgrid variability 
      end if          
 
      return

 END SUBROUTINE IceParam


!===================================================================================
! SUBROUTINE nice_Vdist 
! It calculates the ice crystal number concentration                 
! at the maximum supersaturation (through the SUBROUTINE nice_param) 
! using a PDF of updraft using a sixth order Gauss-Legendre quadrature.
!
! Inputs: T and P (all SI units)
! Output: nice, smax, nhet, nlim
!===================================================================================

 SUBROUTINE nice_Vdist(T_ice, P_ice, nice, smax, nhet, nlim, DC, scice)
      
      real(dp), intent(in)  :: T_ice, P_ice
      real(dp), intent(out) :: nice, smax, nhet, nlim, DC, scice

      real(dp)  :: quadx(6), wpar, quadw(6), dpp, x1, x2,  &
                   sum1, sum2, sum3, sum4, sum5, sum6
      integer :: INDEX

      DATA quadx/0.23861918d0, -0.23861918d0, 0.66120939d0, &      !SB: abscissae for the Legendre-Gauss quadrature function until 6th order
      -0.66120939d0, 0.93246951d0, -0.93246951d0/
      DATA quadw/0.46791393d0, 0.46791393d0, 0.36076157d0, &       !SB: weights for the Legendre-Gauss quadrature function until 6th order
        0.36076157d0, 0.17132449d0,  0.17132449d0/

      !calculate the integral in the denominator
      x1 = (vmin_ice-miuv_ice)/(sqrt(2d0)*sigmav_ice)   
      x2 = (vmax_ice-miuv_ice)/(sqrt(2d0)*sigmav_ice)
      normv_ice = (ERFAPP(x2)-ERFAPP(x1))*0.5d0           !cummulative width of the distribution of velocities 
                                                          !SB: used in SUBROUTINE gausspdf 
      sum1 = 0d0
      sum2 = 0d0
      sum3 = 0d0
      sum4 = 0d0
      sum5 = 0d0
      sum6 = 0d0

      !use a Gauss-Legendre Quadrature approximation until 6th order
      DO INDEX = 1, 6
         wpar = 0.5d0*(((vmax_ice-vmin_ice)*quadx(INDEX)) &
                + (vmax_ice+vmin_ice))  
         
         CALL nice_param(wpar, T_ice, P_ice, nice, smax, nhet, nlim, DC, scice)    !SB: IN(wpar,T_ice,P_ice) OUT(nice,smax,nhet,nlim,DC,scice)
         CALL gausspdf(wpar, dpp)                                                  !SB: IN(wpar) OUT(dp)   

         sum1 = sum1+(nice*dpp*quadw(INDEX))
         sum2 = sum2+(smax*dpp*quadw(INDEX))
         sum3 = sum3+(nhet*dpp*quadw(INDEX))
         sum4 = sum4+(nlim*dpp*quadw(INDEX))
         sum5 = sum5+(scice*dpp*quadw(INDEX))
         sum6 = sum6+(DC*dpp*quadw(INDEX))
      END DO
 
      !SB: convert the limits of integration [vmax_ice, vmin_ice] to the Legendre-Gauss interval [-1,1]
      nice = sum1*(vmax_ice-vmin_ice)*0.5d0          
      smax = sum2*(vmax_ice-vmin_ice)*0.5d0
      nhet = sum3*(vmax_ice-vmin_ice)*0.5d0
      nlim = sum4*(vmax_ice-vmin_ice)*0.5d0
      scice= sum5*(vmax_ice-vmin_ice)*0.5d0
      DC   = sum6*(vmax_ice-vmin_ice)*0.5d0

      return
       
 END SUBROUTINE nice_Vdist

!*************************************************************     
! FUNCTION ERFAPP: Approximation to the error function
!*************************************************************
      real(dp) FUNCTION ERFAPP(x)
       
      real(dp), intent(in) :: x 
      real(dp) :: a  

      a = x*x*(1.27324d0+(0.147d0*x*x))/(1d0+(0.147d0*x*x))
      ERFAPP = SQRT(1d0-exp(-a))

      if (x .lt. 0.0) then     !bug removed 08/03/2011
         ERFAPP = -ERFAPP
      end if 

      END FUNCTION ERFAPP

      
!===================================================================================
! SUBROUTINE nice_param
! Provides the total ice crystal number concentration, considerint the competition
! between homog and heterog freezing.
!
! Inputs: wpar, T and P (all SI units)
! Output: nice, smax, nhet, nlim_
!===================================================================================

 SUBROUTINE nice_param(wpar_ice,T_ice, P_ice, nice, smax, nhet, nlim_, DC, scice)      !mz_sb_20170112: addition DC,scice     
      
      real(dp), intent(in)  :: T_ice, P_ice, wpar_ice
      real(dp), intent(out) :: nice, smax, nhet, nlim_, DC, scice                      !mz_sb_20170112: addition DC,scice
      
      real(dp)  :: AUX1, AUX2, G,                                         &
                   DPMAX, MIU, MONVOL, FDS, NLIM, DLIM, DSTAR, NSTAR,     &
                   NHOM, FC, PHIDO, AUXNC, SIZECORR, DSH, NHET_, F1, F2,  &
                   SAUX, SUP, SLOW, SHALF, FHALF, DPMIN, GAM
      integer :: INDEX 
      
      if (np_ice .gt. 1.0) then   
         !SB: choose the minimum diameter, from BN08 eq(32)
         MONVOL = np_ice*1.0d-6*ddry_ice*ddry_ice*ddry_ice
         AUX1   = 1.6397d-14*T_ice-3.1769d-12
         DPMAX  = AUX1*(MONVOL**(-0.373d0))*(wpar_ice**(-0.05))    !SB: DPMAX = equivalent diameter of the largest ice crystal at Si,max
         if (DPMAX.gt.1.0d-4) then
             DPMAX = 1.0d-4                                        !SB: from BN08 eq(32)
         end if
      else 
         DPMAX = 1.0d-4
      end if 
      
      Pice  = P_ice                                           !SB: pressure
      DPMIN = dliq_ice+(0.02/sqrt(alfa_ice*wpar_ice*g1_ice))  !Minimum size for DPmax (added on 09/11/09)
      DPMAX = max(DPMAX,DPMIN)                                
      
      !SB: Compute the effective growth parameter MEAN GAMMA (BN08 eq(25))
      AUX1  = DPMAX-dliq_ice                                  !SB: dliq_ice = mean diameter of a liquid aerosol particle
      AUX2  = dlog((g2_ice+(g1_ice*DPMAX))/(g2_ice+(g1_ice*dliq_ice)))     
      G     = 1.3346d0*((g1_ice*AUX1)-(g2_ice*AUX2))/(AUX1*g1_ice*g1_ice)      

      lambda_ice = lambda_ice/sqrt(wpar_ice)                  !SB: lambda_ice= 1/sqrt(alfa_ice*g1_ice*gam*gam), part of eq(18) BN09a
      NSTAR = ((g1_ice*alfa_ice*wpar_ice)**1.5d0)/beta_ice/z_ice/sqrt(2d0)     !SB: BN09b eq(26) used to compute NLIM
      GAM   = g2_ice/g1_ice
      
      !*********IS HOMOGENEOUS FREEZING HAPPENING?********************
      !Compute the competition between homog and heterog freezing
      
      FDS  = 1d0  !CORRECTION TO Nc FROM HET FREEZING
      NHOM = 0d0

      if (T_ice .gt. 238.15d0) then      !SB_change: mixed clouds
         NHOM = 0d0                      !no homogeneous above 238 K
         goto 685                        !SB: go to case of pure heterog freezing 
      end if
         
      if (purehet_ice) then         !SB: logical, True supresses homog nucleation
         NHOM = 0d0   
         goto 685
      end if 
       
      !Calculate limiting NIN (NLIM) for combined hom_het    
      if (typeofspec_ice .ge. 0d0) then       !polydisperse expressions
       
         CALL INSPEC_ice(shom_ice, T_ice, NHET_, DSH)        !SB: IN(shom_ice, T_ice) OUT(NHET_, DSH)     

         SIZECORR = EXP(-2d0/lambda_ice/shom_ice)
         DSTAR    = ((4d0*DSH*DSH/3d0)+(2d0*DSH*(shom_ice-DSH))) &                  !SB: delta_s_char* from eq(25) BN09b
                    /(shom_ice-DSH+1d0)
         NLIM     = NSTAR*(shom_ice+1d0)/shom_ice/sqrt(DSTAR)/SIZECORR              !SB: BN09b eq(29)

      else   !monodisperse approximation
         DSH   = shom_ice-sh_ice
         DSTAR = ((4d0*DSH*DSH/3d0)+(2d0*DSH*(shom_ice-DSH))) &
                 /(shom_ice-DSH+1d0)     
         DLIM  = -GAM+sqrt((GAM*GAM)+(2d0*DSTAR/g1_ice/alfa_ice/wpar_ice))          !SB: BN09a eq(18)
         NLIM  = alfa_ice*wpar_ice*(shom_ice+1d0)/z_ice/beta_ice/shom_ice
         NLIM  = NLIM*((g1_ice*DLIM)+g2_ice)/DLIM/DLIM
         NHET_ = nin_ice
      end if
      DC       = sqrt( (2d0*DSTAR)/(alfa_ice*wpar_ice*g1_ice) )                  !ice crystal diameter, eq(25) BN09b 
      
      nlim_ = min(NLIM, 1d10)    
      nlim_ = max(nlim_, 1d-6)
       
      FDS = 1d0-((NHET_/NLIM)**1.5d0)     !SB: BN09b eq(28) 

      if (isnan(FDS)) then       !SB: The returned value of isnan is TRUE if X is a NaN and FALSE otherwise 
         FDS = 1d0
      end if

      if (purehom_ice) then      !SB: logical, True supresses heterog nucleation
         NHET_ = 0d0
         FDS   = 1.0d0
      end if
       
      if (FDS .le. 0d0) then     !Homogeneous nucleation completely inhibited by IN
         NHOM = 0d0   
         goto 685
      end if 

      !Calculate FRACTION OF HOMOGENEOUSLY FROZEN DROPLETS 
      MIU   = FDS*alfa_ice*(shom_ice+1d0)*wpar_ice*koft_ice/shom_ice 
      PHIDO = ((pi_ice*G/MIU/2d0)**0.5d0)*(G/MIU)
      AUXNC = 2d0*denair_ice/beta_ice/koft_ice/denice_ice/pi_ice/np_ice
      FC    = AUXNC/PHIDO                                                 !SB: BN08 eq(29)      

      !Calculate Ice Crystal Number Concentration (Nc) from homog freezing (NHOM)
      if (np_ice .gt.0d0) then
         if (FC .le. 0.6d0) then
            NHOM = np_ice*EXP(-FC)*(1.0d0-EXP(-FC))                       !SB: BN09a eq(30) 
         else
            NHOM = np_ice/(1d0+EXP((9d0-2*FC)/7d0))   !correction needed for convective clouds (very high updraft) (Barahona et al. JGR, 2011)
         end if
      else  
         NHOM = 0d0
      end if
      
      smax = shom_ice
      nhet = NHET_                !SB: NHET_ is provided by the SUBROUTINE INSPEC_ice or from monodisperse population

      goto 686  !finish 

      !********PURE HETEROEGENEOUS FREEZING********************  
      !find interval for bisection
685   smax = 0d0
      nhet = 0d0
      SAUX = 0.01d0 

      if (typeofspec_ice .lt. 0d0) then    !SB: monodisperse population
         SAUX = sh_ice+0.00000000001d0     !minimun smax in monodisperse case
      end if

      F1 = FINDSMAX(SAUX)   

      DO INDEX = 1,20
         if (SAUX + 1d0 .ge. vpresw_ice/vpresi_ice) then   !limit to below water saturation new2
            SHALF = (vpresw_ice/vpresi_ice)-1.0
            goto 678
         end if
         SAUX = SAUX+0.1d0
         F2   = FINDSMAX(SAUX)
         if (F2*F1 .lt. 0d0) goto 677
         F2 = F1
      END DO 
      
      if (F2*F1 .gt. 0d0) then
         nhet = 0d0
         smax = SAUX      !No NIN present in pure heterogeneous mode smax>200%
         goto 686
      end if
      
      !Perform bisection
677   SUP  = SAUX
      SLOW = SAUX-0.1d0
      
      DO INDEX = 1,100
         SHALF = 0.5d0*(SUP+SLOW) 
         FHALF = FINDSMAX(SHALF)
         if (SIGN(1.d0,F1)*SIGN(1.d0,FHALF) .LE. 0d0) then  
            F2   = FHALF
            SUP  = SHALF
         else
            F1   = FHALF
            SLOW = SHALF
         end if

         if (ABS(SLOW-SUP) .le. 1d-5) goto 678
      END DO
       
678   smax=SHALF
 
      if (typeofspec_ice .ge. 0d0) then
         CALL INSPEC_ice(smax, T_ice, nhet, DSH)
      else
         nhet = nin_ice     !monodisperse approximation
      end if
      
      !********Total Ice Crystal Number Concentration********************  

      !SB: to avoid of getting NaN
686   if (isnan(nhet)) then !avoiding errors
         nhet = 0d0
      end if
      
      if (isnan(NHOM)) then !avoiding errors
       NHOM = 0d0
      end if
      
      nice = NHOM + nhet        !SB: Total Ice Crystal Number Concentration

      scice = max(smax+1.0-DSH, 1.0)   
      scice = min(shom_ice+1.0, scice)
      
      return
      
      CONTAINS  
       
!*************************************************************     
! FUNCTION FINDSMAX
!*************************************************************
       real(dp) FUNCTION FINDSMAX(SX)

       real(dp), intent(in) :: SX
       real(dp) :: tao

       if (typeofspec_ice .ge. 0d0) then     !polydisperse expressions

          CALL INSPEC_ice(SX, T_ice, NHET_, DSH) 

          SIZECORR = EXP(-2d0/lambda_ice/SX)
          DSTAR = ((4d0*DSH*DSH/3d0)+(2d0*DSH*(SX-DSH)))/(SX-DSH+1d0)
          DSTAR = DSTAR+(gdoin_ice*alfa_ice*wpar_ice)
          tao   = NHET_*SIZECORR*SX*sqrt(DSTAR)/(SX+1d0)/NSTAR
       else                                  !monodisperse approximation
          DSH   = SX-sh_ice
          DSTAR = ((4d0*DSH*DSH/3d0)+(2d0*DSH*(SX-DSH))) &
                  /(SX-DSH+1d0)
          DLIM  = -GAM+sqrt((GAM*GAM)+(2d0*DSTAR/g1_ice/alfa_ice/wpar_ice))
          tao   = alfa_ice*wpar_ice*(SX+1d0)/z_ice/beta_ice/SX
          tao   = tao*((g1_ice*DLIM)+g2_ice)/DLIM/DLIM/nin_ice
      end if

      DC       = sqrt( (2d0*DSTAR)/(alfa_ice*wpar_ice*g1_ice) )                  !ice crystal diameter, eq(25) BN09b 

      FINDSMAX = 1d0-tao

      END FUNCTION FINDSMAX  
       
 END SUBROUTINE nice_param


!************************************************************* 
! SB: list of functions which will be used later:
!     - VPRESWATER_ice: saturated vapor pressure of water (Pa)
!     - VPRESICE      : saturated vapor pressure of ice (Pa)
!     - DHSUB_ice     : latent heat of sublimation of ice (J/Kg) 
!     - ICEDENSITY    : density of ice (Kg/m3)
!     - WATDENSITY    : density of liquid water (Kg/m3)
!*************************************************************     
! FUNCTION VPRESWATER_ice: it calculates the saturated vapor
! pressure of water (Pa) according to Murphy & Koop (2005)
! T in K (173.15-373.15)
!*************************************************************
      real(dp) FUNCTION VPRESWATER_ice(T)
    
      real(dp), intent(in) :: T
      real(dp) :: A(0:9) 
     
      DATA A/54.842763d0, -6763.22d0,  -4.21d0, 0.000367d0, &
           0.0415d0, 218.8d0, 53.878d0, -1331.22d0, -9.44523d0,  &
           0.014025d0/

      VPRESWATER_ice = A(0)+(A(1)/T)+(A(2)*LOG(T))+(A(3)*T)+ &
                       (TANH(A(4)*(T-A(5)))*((A(6)+(A(7)/T))+ &
                       (A(8)*LOG(T))+ (A(9)*T))) 
      VPRESWATER_ice = EXP(VPRESWATER_ice)           !SB: Murphy & Koop (2005), eq(10)

      return
      END FUNCTION VPRESWATER_ice

!*************************************************************
! FUNCTION VPRESICE. Calculates the saturated vapor pressure
! of ice (Pa) according to Murphy & Koop (2005)
! T in K (>110)
!************************************************************         

      real(dp) FUNCTION VPRESICE(T)
      
      real(dp), intent(in) :: T
      real(dp) :: A(0:3)
      
      DATA A/9.550426d0, -5723.265d0, 3.53068d0, -0.00728332d0/
    
      VPRESICE = A(0)+(A(1)/T)+(A(2)*LOG(T))+(A(3)*T)
      VPRESICE = EXP(VPRESICE)                       !SB: Murphy & Koop (2005), eq(7)
      
      return
      END FUNCTION VPRESICE
      
!*************************************************************
! FUNCTION DHSUB_ice. Calculates the latent heat of sublimation
! of ice (J/Kg) according to Murphy & Koop (2005)
! T in K (>30)
!************************************************************         

      real(dp) FUNCTION DHSUB_ice(T)
      
      real(dp), intent(in) :: T
      real(dp)  :: A(0:4)
      
      DATA A/46782.5d0, 35.8925d0, -0.07414d0, 541.5d0, 123.75d0/
      
      DHSUB_ice = A(0)+(A(1)*T)+(A(2)*T*T)+(A(3) &
                  *EXP(-((T/A(4))**2)))              !SB: Murphy & Koop (2005), eq(5)
      DHSUB_ice = 1000d0*DHSUB_ice/18d0              !SB: from J/mol to J/Kg

      return
      END FUNCTION DHSUB_ice

!*************************************************************
! FUNCTION ICEDENSITY. Calculates the density
! of ice (Kg/m3) according to PK97 
! T in K (>30)
!************************************************************         

      real(dp) FUNCTION DENSITYICE(T)
     
      real(dp), intent(in) :: T
      real(dp)  :: A(0:2), TTEMP
      
      DATA A/0.9167d0, -1.75d-4, -5.0d-7/
      
      TTEMP      = T-273d0                           !SB: from K to C
      DENSITYICE = 1000d0*(A(0)+(A(1)*TTEMP)+(A(2)*TTEMP*TTEMP))

      return
      END FUNCTION DENSITYICE
      
!*************************************************************
! FUNCTION WATDENSITY. Calculates the density
! of liquid water (Kg/m3) according to PK97 
! T in K (>240)
!************************************************************         

      real(dp) function WATDENSITY_ice(T)
      
      real(dp), intent(in) :: T
      real(dp)  :: A(0:6), TTEMP, WATDENSITY
      integer :: I
      
      DATA A/0.99986d0, 6.690d-5, -8.486d-6, 1.518d-7, & 
           -6.9984d-9, -3.6449d-10, -7.497d-12 /
     
      TTEMP = T-273d0                                 !SB: from K to C
      if (TTEMP .le. -40d0) then
         TTEMP = -40d0
      end if
      
      WATDENSITY = A(6)*TTEMP 
      
      if (T .ge. 240.0) then 
        DO I= 5, 1, -1
           WATDENSITY = (WATDENSITY+A(I))*(TTEMP)
        END DO
        WATDENSITY = WATDENSITY+A(0)
      else
        WATDENSITY = 0.979d0
      end if 
      
      WATDENSITY     = WATDENSITY*1000d0
      WATDENSITY_ice = WATDENSITY

      return
      END FUNCTION WATDENSITY_ice
      

!===================================================================================
! SUBROUTINE prop_ice
! It set physical and thermodynamic properties at T and P 
!===================================================================================

 SUBROUTINE prop_ice(T_ice, P_ice, zicncq, sigw, sigw_preice)
      
      real(dp) :: T_ice, P_ice    
      real(dp), intent(in)  :: zicncq, sigw
      real(dp), intent(out) :: sigw_preice

      real(dp) :: AUX1, AUX2, SW, fice, mice, Tc, hdust, hbc, &
                b0, b1, b2, b3, x, T0bc, T0dust, gam, gamma
      real(dp) :: a_ice, c_ice, lambda_pre, w_correction       ! necessary for preexisting ice
      
      !Sanity check
      if (T_ice .gt. 273d0) then  !Use this properties only for ice
         T_ice = 273d0
      end if 
      
      if (T_ice .lt. 180d0) then  !Use this properties only for ice
         T_ice = 180d0
      end if 
      
      if (P_ice .lt. 1d0) then    !Use this properties only for ice
         P_ice = 1d0
      end if 
      
      !mz_sb_20161219+       
      wmw_ice = M_H2O*10**(-3d0)   !Molar mass of H2O [kg/mol]
      amw_ice = M_air*10**(-3d0)   !Molar mass of dry air [kg/mol]
      !mz_sb_20161219-       

      rv_ice     = rgas_ice/wmw_ice                !SB: specific gas constant for water
      dhs_ice    = DHSUB_ice(T_ice)                !SB: latent heat of sublimation of ice 
      vpresw_ice = VPRESWATER_ice(T_ice)           !SB: saturated vapor pressure of water
      vpresi_ice = VPRESICE(T_ice)                 !SB: saturated vapor pressure of ice
      denice_ice = DENSITYICE(T_ice)               !SB: density of ice
      denwat_ice = WATDENSITY_ice(T_ice)           !SB: density of water
      denair_ice = P_ice*amw_ice/rgas_ice/T_ice    !SB: density of air (from P = ro R* T)

      ! Kinetic properties of the bulk vapor (SI UNITS, Seinfel and Pandis, 1997)
      diff_ice = (0.211d0*101325d0/P_ice)*((T_ice/273d0)**1.94d0)*1.0d-4      !SB: water vapor diffusivity (m^2/s), SP97 eq(17.61)

      AUX1 = 1.0e-3*(4.39d0+0.071d0*T_ice)                                    !SB: thermal conductivity of air (W/m/K), SP97 eq(17.71)
      AUX2 = (2d0*AUX1/(thaccom_ice*1.0d-6*denair_ice*cpa_ice)) &             !correcting Kair for size assuming D=1e-6 m
             *((58.0d-3*pi_ice/(rgas_ice*T_ice))**0.5d0)
      aircond_ice = AUX1/(1.d0+AUX2)                                          !SB: thermal conductivity of air accounting for non-continuum effects 
                                                                              !SP97 eq(17.72)
      !Physical constants
      AUX1     = grav_ice*dhs_ice/rv_ice/T_ice/T_ice/cpa_ice                  
      AUX2     = grav_ice*amw_ice/rgas_ice/T_ice
      alfa_ice = AUX1-AUX2                                         !SB: parameter def in BN08 eq(14)
      beta_ice = amw_ice*P_ice/wmw_ice/vpresi_ice                  !SB: parameter def in BN08 eq(14) 
      gamma    = 1.5d0*dhs_ice*dhs_ice/rv_ice/T_ice/T_ice/cpa_ice  !Correction for T>250 K
      beta_ice = beta_ice+gamma                                    !only important for high T (>250 K, Barahona et al. JGR 2010)

      !Homogeneous freezing only 
      shom_ice = 2.349d0-(T_ice/259d0)                             !hom threeshold Si according to Ren & McKenzie, 2005
      SW       = shom_ice*vpresi_ice/vpresw_ice                    !SW = Water Sat. ratio
      shom_ice = shom_ice-1d0                                      !homog freezing ice supersaturation
      koft_ice = (0.0240d0*T_ice*T_ice)-(8.035d0*T_ice)+934.0d0    !constant related to Jmax for homog freezing, BN08 eq(18) --> k(T)

      !Calculate Dliq using an approximation derived from the equilbrium calculations and the
      !approximation proposed by Lewis (2008), 13, D03205, JGR 
      if (SW .lt. 0.99) then                        !only subsaturated regime (Haze Aerosols)
         AUX1     = (1d0/(1d0-SW))-1.1764d0
         !dliq_ice = ddry_ice*0.9344d0*(AUX1**(1d0/3d0))               !SB
      else
         AUX1 = (1d0/0.01)-1.1764d0         
      end if 
      dliq_ice = ddry_ice*0.9344d0*(AUX1**(1d0/3d0))                   !SB_change: after Sylvia's answer !cloud droplet diameter
       

      !Computations needed for the correction of the vertical velocity, in order to take into account the effect of preexisting ice
      !Such correcited velocity is defined in Barahona et al. 2014, eq(24). 
      !Variable Ai in eq(24):
      AUX1  = (denice_ice*dhs_ice*dhs_ice)/(aircond_ice*rv_ice*T_ice*T_ice) 
      AUX2  = (denice_ice*rv_ice*T_ice)/(vpresw_ice*diff_ice)
      a_ice = 1d0/(AUX1+AUX2)

      !Variable c in eq(24) = shape factor assumed to be equal to 1: 
      c_ice = 1d0
  
      !Variable lamba in eq(24) = slope of the size distribution of preexisting ice crystals.
      !Since the model does not have a size distribution for ice, we can assume lambda=1
      lambda_pre = 1d0

      !Correction in eq(24):
      AUX1 = zicncq*pi_ice*beta_ice*c_ice*denice_ice*a_ice*shom_ice           !numerator 
      AUX2 = 2d0*lambda_pre*alfa_ice*sigw*(shom_ice+1d0)              !denominator
      w_correction = AUX1/(max(AUX2,zeps))

      !The preexisting ice effect is considered via eq(24) of Barahona et al. 2014:
      AUX1 = max((1d0-w_correction), 0d0) 
      sigw_preice = sigw*AUX1                                      !Corrected vertical velocity

      !Compute GAMMA1 and GAMMA2 (needed for MEAN GAMMA, G)
      AUX1   = denice_ice*rv_ice*T_ice/vpresi_ice/diff_ice             !SB: AUX1 and AUX2 are factors to compute GAMMA1 
      AUX2   = dhs_ice*denice_ice/aircond_ice/T_ice
      AUX2   = AUX2*((dhs_ice/rv_ice/T_ice)-1.0d0)
      !SB: GAMMA1 and GAMMA2                                           !SB: BN08 eq(6)
      g1_ice = (AUX1+AUX2)/4.0d0
      g2_ice = denice_ice*rv_ice*T_ice/2.0d0/vpresi_ice/depcoef_ice   
      g2_ice = g2_ice*((2.0d0*pi_ice/rv_ice/T_ice)**0.5d0)                   
      
      !SB: some "factors" used in other computations
      doin_ice  = 1.0d-6                             !assumed IN diameter at freezing
      gdoin_ice = (g1_ice*0.5d0*doin_ice*doin_ice)+(g2_ice*doin_ice)
      z_ice     = denice_ice*pi_ice/2.0d0/denair_ice
      gam       = g2_ice/g1_ice
      lambda_ice= 1d0/sqrt(alfa_ice*g1_ice*gam*gam)  !divided by sqrt(wparcel) in niceparam 
      
      !============Parameters needed for IN spectra:

      !***** Parameters used in the Monodisperse approximation (CASE -1) *****
      sh_ice  = 0.3d0                                !assumed freeizng threshold
      nin_ice = (sum(ndust_ice)+nbc_ice)*0.05d0      !SB: it is used only for monodisperse aerosols
       
      !***** Parameters needed for CNT spectra, BN08 (CASE 3)***** 
      shdust_ice = 0.1d0                             !maximum freezing threshold dust only used in CNT, case(3) 
      effdust_ice= 0.6d0                             !maximum freezing efficiency dust
      shbc_ice   = 0.35d0                            !maximum freezing threshold bc only used in CNT, case(3)
      effbc_ice  = 0.01d0                            !maximum freezing efficiency bc
  
      !SB: k for heterog nucl BN09b eq(4), for dust and BC
      !---- Dust 
      mice       = 0.96d0                            !compatibility parameter dust !SB: it is the cos(contact angle), BN09b eq(4) --> m
      fice       = 0.25d0*((mice*mice*mice)-(3d0*mice)+2d0)   !SB: BN09b eq(4) --> f_h 
      kdust_ice  = koft_ice*fice                              !SB: BN09b eq(4) 
      !---- Black Carbon
      mice       = 0.76d0                            !compatibility parameter bc
      fice       = 0.25d0*((mice*mice*mice)-(3d0*mice)+2d0)   !SB: BN09b eq(4)
      kbc_ice    = koft_ice*fice      
      
      !***** Parameters needed for spectra PDA08, Phillips et. al. 2008 (CASE 4) *****
      !---- Dust 
      Tc    = T_ice-273.15d0
      hdust = 0.15d0                                 !SB: Table 1 in Phillips et. al. 2008
      T0dust= -40d0                                  !SB: Table 1 in Phillips et. al. 2008
      b0    = -1.0261d0
      b1    = 3.1656d-3
      b2    = 5.3938d-4
      b3    = 8.2584d-6
      x     = b0+(b1*Tc)+(b2*Tc*Tc)+(b3*Tc*Tc*Tc)    !SB: Table 1 in Phillips et. al. 2008       
      si0dust_ice  = 1d0+(10d0**x)                   !SB: Threshold of saturation ratio of vapor with respect to ice for dust and metallic aerosols
      del1dust_ice = cubicint_ice(Tc, T0dust, T0dust+5d0, 1d0, hdust) !SB: contribution to compute FC of eq(12) in Phillips et. al. 2008
      
      !---- Black Carbon
      hbc  = 0d0
      T0bc = -50d0
      b0   = 0.5652d0
      b1   = 1.085d-2 
      b2   = -3.118d-5
      si0bc_ice  = b0+(b1*T_ice)+(b2*T_ice*T_ice)-0.1d0  !bug corrected C to K
      del1bc_ice = cubicint_ice(Tc, T0bc, T0bc+5d0, 1d0, hbc)
      
      return
      
 END SUBROUTINE prop_ice


!===================================================================================
! SUBROUTINE gausspdf
! Gaussian distribution normalized over the width of the distribution
! It is called only in nice_Vdist, wpar --> x 
!
! miuv_ice   = mean velocity for the V distribution
! sigmav_ice = standard deviation of the V distribution (m/s)
! normv_ice  = cummulative width of the distribution of velocities: 
!              (ERFAPP(x2)-ERFAPP(x1))*0.5d0
!===================================================================================
      
 SUBROUTINE gausspdf(x, dpp)
   
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: dpp
      
      dpp  = EXP(-0.5d0*(x-miuv_ice)*(x-miuv_ice)/sigmav_ice/sigmav_ice) &    
            /sigmav_ice/sqrt(2d0*pi_ice)/normv_ice 

      return
      
 END SUBROUTINE gausspdf
      

!************************************************************* 
! SB: list of functions
!     - cubicint_ice  : cubic interpolation (Phillips et al. 2008)
!     - dcubicint_ice : used in the PDA08 spectrum
!     - PDG07_ice     : latent heat of sublimation of ice (J/Kg) 
!*************************************************************
! FUNCTION cubicint_ice: cubic interpolation between y1 and y2 
! within a and b.
! It is defined in Phillips et al. 2008, APPENDIX A  
!************************************************************  

      real(dp) FUNCTION cubicint_ice(y, y1, y2, a, b)
      
      real(dp), intent(in) :: y, y1, y2, a, b   
      real(dp) :: A_, B_, a0, a1, a2, a3, d, AUX
      
      if (y .le. y1) then
         d = a
         goto 5065               !cubicint_ice=d
      end if 
     
      if (y .ge. y2) then
         d = b
         goto 5065
      end if 

      AUX = y2-y1      
      A_  = 6d0*(a-b)/(AUX*AUX*AUX)
      B_  = a+(A_*(y1*y1*y1)/6d0)-(A_*(y1*y1)*y2*0.5d0)
      
      a0  = B_
      a1  = A_*y1*y2
      a2  = -A_*(y1+y2)*0.5d0
      a3  = A_/3d0
      d   = a0+(a1*y)+(a2*y*y)+(a3*y*y*y)
          
5065  cubicint_ice=d
      
      END FUNCTION cubicint_ice
      
!*************************************************************
! FUNCTION dcubicint_ice (used in the PDA08 spectrum).  
!*************************************************************     

      real(dp) FUNCTION dcubicint_ice(y, y1, y2, a, b)
     
      real(dp), intent(in) :: y, y1, y2, a, b   
      real(dp) :: A_, a1, a2, a3, d, AUX
      
      if (y .le. y1) then
         d = 0
         goto 5065            !dcubicint_ice = d
      end if 
     
      if (y .ge. y2) then
         d = 0
         goto 5065
      end if 

      AUX = y2-y1      
      A_  = 6d0*(a-b)/(AUX*AUX*AUX)
     
      a1  = A_*y1*y2
      a2  = -A_*(y1+y2)*0.5d0
      a3  = A_/3d0
      d   = (a1)+(2d0*a2*y)+(3d0*a3*y*y)
      
5065  dcubicint_ice = d
      
      END FUNCTION dcubicint_ice      
      
!*************************************************************
! FUNCTION PDG07_ice: simplified ice nucleation 
! spectra according to Phillips et. al. 2007 (PDG07)  
! to compute condensation/deposition nucleation. 
! si = supersaturation wrt ice and T is in K 
!************************************************************  

      real(dp) FUNCTION PDG07_ice(si, T)     
      
      real(dp), intent(in) :: si, T
      real(dp) :: N 
      
       if (T .le. 243d0) then
          N = 1000d0*exp(-0.388d0)*(exp(3.88d0*si)-1d0)     !SB: Table1 in BN09b
       else
          N = 60d0*exp(-0.639d0)*(exp(12.96d0*si)-1d0)      !SB: Table1 in BN09b, PDG07 eq(11), similar to Meyers
       end if
       
       PDG07_ice = N
       
      END FUNCTION PDG07_ice     
         


!===================================================================================
! SUBROUTINE INSPEC_ice
! Provides the IN nucleation spectrum in its cumulative form (=Nhet(si) {m-3} Table1 BN09b)  
! and the characteristic freezing threeshold, Dsh (Barahona & Nenes 2009), at given 
! si and T (si is supersaturation wrt ice and T is in K). 
! The variable typeofspec_ice (integer) has the values:
! 1  Meyers et. al. 1992
! 2  Phillips et. al. 2007, PDG07
! 3  Barahona & Nenes 2009, CNT: assumme a maximum freezing fraction then scale it
! 4  Phillips et. al. 2008, PDA08 
! 5  Phillips et. al. 2013
!
! Written by Donifan Barahona 
! donifan.o.barahona@nasa.gov
!===================================================================================

 SUBROUTINE INSPEC_ice(six, T, N, Dsh)
      
      real(dp), intent(in)  :: six, T
      real(dp), intent(out) :: N, Dsh

      real(dp) ::  Nd, Nbc, aux, Si_, SW, del0, ddel0, &
                 fc, delw0, ddelw0, SW0, Hdust, Hbc, &
                 Nbase, dNd, dNbc, dNbase, dH, &
                 dfc, Ndaux, dNdaux, dNorg, Norg, Ndustaux,&
                 frac, aux2, Dx2, fdep,  si, &
                 SIW
      real(dp), dimension(3) :: sig_array, the_array, frac_array
      real(dp) :: n_iw, DSh_s, nbc_s, dbc_s, Asolo, Asolo_s, nbio_s, dbio_s
      real(dp), dimension (nbindust_ice) ::  ndust_s, ddust_s       
      integer :: index
       
      si  = six                          !SB: ice supersaturation       
      Si_ = si+1d0                       !SB: ice saturation ratio
      SW  = Si_*vpresi_ice/vpresw_ice    !new2 !SB:Water Sat. ratio
      
      if (SW .ge. 1.0) then              !limit to subsaturated conditions
         SW  = 1.0
         Si_ = vpresw_ice/vpresi_ice     !SB: Ice Sat. ratio assuming water saturation
         si  = Si_-1.0                
      end if

      SIW        = vpresw_ice/vpresi_ice !SB: Ice Sat. ratio assuming water saturation 
      sig_array  = 0.0
      the_array  = 0.0
      frac_array = 0.0
      fdep       = 0.0
      
      !***** SB: Select the different cumulative freezing spectra *****
      select case  (typeofspec_ice)
       
      case(1)  !****Meyers et al. 1992    
          !Condensation/deposition nucleation:
          N   = 1000d0*exp(-0.639d0)*(exp(12.96d0*si)-1d0)     !SB: Meyers eq(2.4) and Table1 in BN09b
                                                               !SB: the factor *1000 change the unit from [1/L] to [1/m3]
          Dsh = 1d0/12.96d0

      case(2)  !****Phillips et al. 2007 (PDG07)
          !Condensation/deposition nucleation:
          N = PDG07_ice(si, T)

          if (T .le. 243d0)then
             Dsh = 1d0/3.88d0
          else
             Dsh = 1d0/12.96d0
          end if

      case(3)  !****Barahona&Nenes 2009b (CNT)
          !---- dust contribution
          Ndustaux = 0.0d0
          DO INDEX = 1, nbindust_ice
             Ndustaux = Ndustaux+ndust_ice(index)              !SB: dust aerosol number concentration
          END DO
            
          if (si .le. shdust_ice) then                         !SB: to choose the minimum Nhet-dust
             Nd  = (si/shdust_ice)*Ndustaux*effdust_ice* &     !SB: BN09b eq(5) and Table1 
                   exp(-kdust_ice*(shdust_ice-si))             !SB: shdust_ice = 0.1, effdust_ice = 0.6
             dNd = Nd*((1d0/si)+kdust_ice)          
          else
             Nd  = Ndustaux*effdust_ice
             dNd = 0d0
          end if
           
          !---- soot contribution
          if (si .le. shbc_ice) then                           !SB: to choose them minimum Nhet-soot
             Nbc  = (si/shbc_ice)*nbc_ice*effbc_ice* &         !SB: BN09b eq(5) and Table1 
                    exp(-kbc_ice*(shbc_ice-si))                !SB: shbc_ice = 0.0.35, effbc_ice = 0.01
             dNbc = Nbc*((1d0/si)+kbc_ice)
          else
             Nbc  = nbc_ice*effbc_ice
             dNbc = 0d0
          end if
         
          N = Nd+Nbc

          if ((dNd+dNbc) .gt. 0d0) then                        !SB: BN09b eq(24)
             Dsh = N/(dNd+dNbc)
          else
             Dsh = si
          end if 
        
      case(4)  !****Phillips et. al. 2008 (PDA08)
                !Allows multiple lognormal modes for dust. Single mode lognormal distributions are assumed for bc and organics
          !---- dust
          SW0   = 0.97d0
          delw0 = cubicint_ice(SW, SW0, 1d0, 0d0, 1d0)
          ddelw0= dcubicint_ice(SW, SW0, 1d0, 0d0, 1d0)
          Nbase = PDG07_ice(si, T)
          
          if (T .le. 243d0)then
             dNbase = 3.88d0*Nbase
          else
             dNbase = 12.96d0*Nbase
          end if
          
          !--- dust contribution
          del0  = cubicint_ice(Si_, si0dust_ice, si0dust_ice+0.1d0, 0d0, 1d0)
          ddel0 = dcubicint_ice(Si_, si0dust_ice, si0dust_ice+0.1d0, 0d0, 1d0)
          
          fc    = 0.5d0*del1dust_ice*del0
          dfc   = 0.5d0*del1dust_ice*ddel0
          Hdust = fc+((1d0-fc)*delw0) 
          dH    = (dfc*(1d0-delw0))+(ddelw0*(1d0-fc))
          
          if (Hdust .gt. 1d0) then 
              Hdust = 1d0
              dH    = 0d0
          end if
          
          aux  = (2d0/3d0)*Hdust*(Nbase/0.76d0)*pi_ice/5.0d-7/4d0  
          aux2 = (2d0/3d0)*pi_ice/0.76d0/5.0d-7/4d0 !The last 4d0 was introduced as recommnedation of V Phillips

          Nd   = 0d0
          dNd  = 0d0
          
          DO INDEX =1, nbindust_ice
             Dx2  = ddust_ice(index)*ddust_ice(index)*ddust_ice(index)*0.52*acorr_dust !new 07/10/12
             frac = 0.5d0*(1d0-ERFAPP(-log(ddust_ice(index)/0.1e-6) & !fraction above 0.1 microns
                    /sigdust_ice(index)/sqrt(2d0)))
             Ndaux= frac*ndust_ice(index)*(1d0-exp(-aux*Dx2))
             Nd    = Nd+Ndaux
             Ndaux = (frac*ndust_ice(index)-Ndaux)
             dNdaux= Ndaux*((dH*Nbase)+(Hdust*dNbase))*aux2*Dx2
             dNd   = dNd+dNdaux
          END DO
            
          !---- soot contribution
          del0  = cubicint_ice(Si_, si0bc_ice, si0bc_ice+0.1d0, 0d0, 1d0)
          ddel0 = dcubicint_ice(Si_, si0bc_ice, si0bc_ice+0.1d0, 0d0, 1d0)
           
          fc    = 0.5d0*del1bc_ice*del0
          Hbc   = fc+((1d0-fc)*delw0)
          dfc   = 0.5d0*del1bc_ice*ddel0
          dH    = (dfc*(1d0-delw0))+(ddelw0*(1d0-fc))
          
          if (Hbc .gt. 1d0) then 
              Hbc = 1d0
              dH  = 0d0
          end if
          
          frac = 0.5d0*(1d0 -ERFAPP(-log(dbc_ice/0.1e-6) & 
                 /sigbc_ice/sqrt(2d0)))

          Dx2  = dbc_ice*dbc_ice*dbc_ice*0.52*acorr_bc  !new 07/10/12
          aux  = ((1d0/3d0)-0.06d0)*Hbc*(Nbase/0.76d0)*pi_ice/2.7d-7 
          aux2 = ((1d0/3d0)-0.06d0)*pi_ice/0.76d0/2.7d-7
           
          Nbc  = nbc_ice*frac*(1d0-exp(-aux*Dx2))
          dNbc = (nbc_ice*frac-Nbc)*((dH*Nbase)+(Hbc*dNbase))*aux2*Dx2
          
          !---- organics contribution
          frac = 0.5d0*(1d0-ERFAPP(-log(dorg_ice/0.1e-6) & 
                 /sigorg_ice/sqrt(2d0)))
          Dx2  = dorg_ice*dorg_ice
          aux  = 0.06d0*Hbc*(Nbase/0.76d0)*pi_ice/9.1d-7 
          aux2 = 0.06d0*pi_ice/0.76d0/9.1d-7
          
          Norg = norg_ice*frac*(1d0-exp(-aux*Dx2))
          dNorg= (norg_ice*frac-Norg)*((dH*Nbase)+(Hbc*dNbase))*aux2*Dx2

          N = Nd+Nbc+Norg
       
          if ((dNd+dNbc+dNorg) .gt. 0d0) then
             Dsh = N/(dNd+dNbc+dNorg)
          else
             Dsh = si
          end if 
         
      case(5)  !****Phillips et al. 2013
           !---- dust
           DO INDEX = 1, nbindust_ice
              frac = 0.5d0*(1d0-ERFAPP(log(0.1e-6/ddust_ice(index)) &     !fraction above 0.1 microns = 10^(-7) 
                     /sigdust_ice(index)/sqrt(2d0)))
              if (isnan(frac)) then                !To avoid errors because of log in frac
                 frac = 0d0
              end if 
              ndust_s(index) = frac*ndust_ice(index)
              ddust_s        = ddust_ice(index)
           END DO

           !---- black carbon
           frac  = 0.5d0*(1d0-ERFAPP(log(0.1e-6/dbc_ice) &                !fraction above 0.1 microns
                   /sigbc_ice/sqrt(2d0)))
           if (isnan(frac)) then                !To avoid errors because of log in frac
              frac = 0d0
           end if 
           nbc_s = frac*nbc_ice
           dbc_ice = min(dbc_ice, 2.e-6)
           dbc_s   = dbc_ice*1.0
             
           !---- soluble organics (spherical)
           Asolo_s = 0d0
           DO INDEX = 1, nbinsolo_ice
              frac   = 0.5d0*(1d0-ERFAPP(log(0.1e-6/dsolo_ice(index)) &     !SB: no 0.25 anymore  
                      /sigsolo_ice(index)/sqrt(2d0)))
              if (isnan(frac)) then                !To avoid errors because of log in frac
                 frac = 0d0
              end if 
              Asolo  = frac*nsolo_ice(index)*pi_ice*dsolo_ice(index)*dsolo_ice(index) 
              Asolo_s= Asolo_s + Asolo
           END DO 
           
           !---- biological aerosols
           nbio_s = nbio_ice
           dbio_s = dbio_ice

           CALL EMPIRICAL_PARAM_PHILLIPS(T, Si_, SIW, SW,                         &
                                         ddust_s, ndust_s, 2,                     &
                                         (/dbc_s/), (/nbc_s/), 1,                 &
                                         (/dbio_s/), (/nbio_s/), 1,               &
                                         Asolo_s, n_iw, DSh_s)                                    !SB: OUTPUT(n_iw, DSh_s)
                                      
           N   = n_iw
           Dsh = DSh_s
            
      case default 
        N   = 0d0
        Dsh = si

      end select
      
      if (Dsh .ge. si) then
         Dsh = si
      end if 
      
 END SUBROUTINE INSPEC_ice
 

  
!===================================================================================
! SUBROUTINE EMPIRICAL PARAMETERISATION 
! (Phillips et al. 2013, JAS)
!
! Contributed by Vaughan Phillips, 2012
! University of Leeds
! Implementation:   Donifan Barahona donifan.o.barahona@nasa.gov 
!===================================================================================
! SB: definitions taken from Phillips et al. 2008 and 2013
!
! INPUT variables:
!    temperature_K = temperature [K]
!    SI            = saturation ratio of vapor wrt ice
!    SIW           = value of SI at exact water saturation
!    SW            = saturation ratio of vapor wrt (plane) water 
!    D_grid_dust   = ddust
!    n_grid_dust   = ndust
!    ijstop_dust   = modes of dust (=2)
!    D_grid_soot   = dbc  
!    n_grid_soot   = nbc  
!    ijstop_soot   = modes of soot (=1)
!    D_grid_bio    = diameter of biological particles = 0 
!    n_grid_bio    = number concentration of biological particles = 0
!    ijstop_bio    = modes of biological particles (=1)
!    A_solo        = ideal total surface area of all soluble organic aerosol [m2/m3]
!
! INPUT/OUTPUT variables: 
!    n_iw          =
!    DSH           = 
!
! Other variables:
!    mu            = average of the number of activated ice embryos per insoluble aerosols in gruop X  [-] 
!    n_in_dust     = contribution (number mixing ratio) to n_in from DU aerosols                       [1/m3]    
!    n_in_soot     = contribution (number mixing ratio) to n_in from BC aerosols                       [1/m3]
!    n_in_bio      = contribution (number mixing ratio) to n_in from bio aerosols                      [1/m3]
!    n_in_solO     = contribution (number mixing ratio) to n_in from soluble organic aerosols          [1/m3]
!    n_in_solO_star = number mixing ratio of IN in soluble organics in AIDA                            [1/kg] 
!    H_frac_dust, H_frac_soot, H_frac_solO, H_frac_bio = empirically determined fraction 0 <= Hx <= 1  [-]       
!    dep_frac      = contribution to Hx from modes of nucleation from aerosol groups                   [-]   
!    n_in          = average number concentration of IN for a background-troposphere scenario          [1/m3]
!
!===================================================================================

 SUBROUTINE EMPIRICAL_PARAM_PHILLIPS(temperature_K, SI, SIW, SW,                 &   
                                     D_grid_dust, n_grid_dust, ijstop_dust,      &  !SB: D_grid_dust=ddust, n_grid_dust=ndust,  ijstop_dust=2
                                     D_grid_soot, n_grid_soot, ijstop_soot,      &  !SB: D_grid_soot=dbc,   n_grid_soot=nbc,    ijstop_soot=1
                                     D_grid_bio,  n_grid_bio,  ijstop_bio,       &  !SB: D_grid_bio=0e0,    n_grid_bio=zero_par,ijstop_bio=1
                                     A_solo, n_iw, DSH)

      implicit none

      integer, intent(in)                :: ijstop_dust, ijstop_soot, ijstop_bio
      real(dp), intent(in)               :: temperature_K, SI, SIW, SW, A_solo      
      real(dp), dimension(:), intent(in) :: D_grid_dust, n_grid_dust, &                    
                                            D_grid_soot, n_grid_soot, &
                                            D_grid_bio, n_grid_bio
      real(dp), intent (inout)           :: DSH, n_iw                                     

                  !dust            !soot            !bio            !solO
      real(dp) :: n_in_dust,       n_in_soot,       n_in_bio,       n_in_solO,               &
                  n_in_dust_ultra, n_in_soot_ultra, n_in_bio_ultra, n_in_solO_star,          &
                  dn_in_dust,      dn_in_soot,      dn_in_bio,      dn_in_solO,              &
                  H_frac_dust,     H_frac_soot,     H_frac_bio,     H_frac_solO,             &
                  dH_frac_dust,    dH_frac_soot,                    dH_frac_solo,            &
                  CIHENC_dust,     CIHENC_soot,     CIHENC_bio,     CIHENC_solO

      real(dp) :: n_in, n_in_hat, n_in_tilde, n_in_ultra, n_in_max,  &  
                  Psi_solO, aux, naux, dNaux, dNall,                 &
                  SS_w, SS_i, SS_iw, S_i_0, S_w_0,                   &
                  dfdep, dep_frac, FAC_CORRECT_RH = 2.,              &
                  RHI, tc_HM_degC, mu, rho_AIDA  
      integer:: ij

      !Parameters 
      real(dp) :: BASE_DUST_OMEGA,    BASE_SOOT_PHILIC_OMEGA, BASE_BIO_OMEGA, BASE_SOLO_OMEGA,     &
                  ALPHA_DUST,         ALPHA_SOOT,             ALPHA_bio,                           &
                  TEMP_MAX_DUST_DEGC, TEMP_MAX_SOOT_DEGC,     TEMP_MAX_bio_DEGC,                   &
                  GLASS_FRAC, RHO_CFDC, FRACTION_DEPNUCL_WARM_DUST, PIE

      PARAMETER(BASE_DUST_OMEGA = 2.0e-6,               & !SB: component of OMEGA* of DU in backgroud-troposphere scenario    [m2/kg], PDA13
                BASE_SOOT_PHILIC_OMEGA = 1.e-7,         & !SB: component of OMEGA* of BC in backgroud-troposphere scenario    [m2/kg], PDA13 
                BASE_BIO_OMEGA  = 0.89e-6,              & !SB: component of OMEGA* of bio in backgroud-troposphere scenario   [m2/kg], PDA13
                BASE_SOLO_OMEGA = 5.6e-5,               & !SB: component of OMEGA* of OCsol in backgroud-troposphere scenario [m2/kg], PDA13 
                ALPHA_DUST = 2./3.,                     & !SB: fraction of n_in* from IN activity of group X=DU               [-], PDA08 
                ALPHA_SOOT = 1./3.-0.03,                & !SB: fraction of n_in* from IN activity of group X=BC               [-], PDA08 
                ALPHA_bio = 0.03,                       & !SB: fraction of n_in* from IN activity of group X=bio              [-], PDA08
                TEMP_MAX_DUST_DEGC = -10.,              & !SB: Table1 PDA13
                TEMP_MAX_SOOT_DEGC = -15.,              & !SB: Table1 PDA13 
                TEMP_MAX_bio_DEGC = -2.,                & !SB: Table1 PDA13
                GLASS_FRAC      = 0.5,                  & !SB: glass fraction                                                 [-], Table1 PDA13  
                RHO_CFDC = 50000./(287.*228.15),        & !SB: air density inside the CFDC = 0.76                             [kg/m3], PDA08
                FRACTION_DEPNUCL_WARM_DUST = 0.15)
      
      !SB: check number of modes
      if(ijstop_dust .ne. SIZE(n_grid_dust)) stop 6366
      if(ijstop_soot .ne. SIZE(n_grid_soot)) stop 6366
      if(ijstop_bio .ne. SIZE(n_grid_bio)) stop 6366
      
      !SB: Initialization
      aux    = 0.0
      dNaux  = 12.96 !default
      naux   = 0.0
      dNall  = dNaux
      DSh    = 0.0
      n_iw   = 0.0
      dH1smooth  = 0.0

      !Dust                    !Soot                    !Bio                    !Soluble Organics 
      !nin_a_nuc_dust  = 0.0;   nin_a_nuc_soot  = 0.0;   nin_a_nuc_bio  = 0.0;   nin_a_nuc_solo  = 0.0
      !num_ic_dust_imm = 0.0;   num_ic_soot_imm = 0.0;   num_ic_bio_imm = 0.0;   num_ic_solo_imm = 0.0 
      dn_in_dust      = 0.0;   dn_in_soot      = 0.0;   dn_in_bio      = 0.0;   dn_in_solO      = 0.0
      dH_frac_dust    = 0.0;   dH_frac_soot    = 0.0;                           dH_frac_solo    = 0.0
     
      !************************************************************
      !         COMPUTATION BLOCK 
      !************************************************************

      PIE = pi_ice
 
      SS_i     = SI - 1.0     !everything is based on supersaturation 
      SS_w     = SW - 1.0
      SS_iw    = SIW - 1.0    !SS_iw = QSW/QSI - 1.
      
      if (SS_i > 0.0) then
         if (temperature_K < 273.15 .and. temperature_K > 273.15 - 90. ) then    !SB: 180.15 < T < 273.15

            if (SS_w  > 0.) then
               SS_i = SS_iw
               SS_w = 0.0
            end if 

            !------ dust ----------------------
            tc_HM_degC = temperature_K - 273.15                    !SB: Temperature in Celsius
            S_i_0 = 1. + 10.**(8.2584e-6*tc_HM_degC*tc_HM_degC*tc_HM_degC + 5.3938E-4*tc_HM_degC*tc_HM_degC &
                    + 3.1656E-3*tc_HM_degC - 1.0261)
            S_w_0 = 0.97
            
            !SB: to compute the contribution to H_DU from modes of nucleation, eq(12) PDA08 
            aux      = H_1_smooth(-(temperature_K-273.15), 35.d0, 40.d0, FRACTION_DEPNUCL_WARM_DUST, 1.d0)/FAC_CORRECT_RH
            dep_frac = H_1_smooth(SS_i + 1d0, S_i_0, S_i_0 + 0.1d0, 0.d0,1.d0)*aux               !SB: eq(12) PDA08

            dfdep    = dH1smooth*aux                                          
           
            !SB: to compute the empirical fraction, eq(11) PDA08             
            aux          = H_1_smooth(SS_w + 1.0d0, S_w_0, 1.d0, 0.d0,1.d0)
            H_frac_dust  = dep_frac + (1. - dep_frac)*aux                         !SB: eq(11) PDA08                          
            if (H_frac_dust > 1.) H_frac_dust = 1.                                !SB: take the min, as eq(11) PDA08 

            dH_frac_dust = dfdep + (SIW*(1. - dep_frac)*dH1smooth)- aux*dfdep
            
            if ((H_frac_dust .gt. 1.0e-6) .and. (H_frac_dust .lt. 1.)) then 
               dH_frac_dust = dH_frac_dust/H_frac_dust
            else
               dH_frac_dust = 0.0
            end if 

            !------ soluble organics ----------
            S_i_0 = 1.2
            
            !SB: to compute the contribution to H_OCsol from modes of nucleation, as eq(12) PDA08
            aux         = H_1_smooth(-(temperature_K-273.15), 65.d0, 75.d0, 0.d0,1.d0)
            dep_frac    = H_1_smooth(SS_i + 1d0, S_i_0, S_i_0+0.1d0, 0.d0,1.d0)*aux   

            H_frac_solO = dep_frac  
            if (H_frac_solO > 1.) H_frac_solO = 1.  

            dH_frac_solo = 0.0

            if ((H_frac_solo .gt. 1.0e-6) .and. (H_frac_solo .lt. 1.)) then  
               dH_frac_solo = dH1smooth/H_frac_solo
            end if   

            !------ soot ----------------------
            S_w_0 = 0.97
            S_i_0 = 1.3
            
            !SB: to compute the contribution to H_BC from modes of nucleation, as eq(12) PDA08
            aux      = H_1_smooth(-(temperature_K-273.15), 40.d0, 50.d0, 0.d0,1.d0)/FAC_CORRECT_RH
            dep_frac = H_1_smooth(SS_i + 1d0, S_i_0, S_i_0+0.1d0, 0.d0,1.d0)* aux 

            dfdep    = dH1smooth*aux
            
            !SB: to compute the empirical fraction, eq(11) PDA08
            aux         = H_1_smooth(SS_w + 1.0d0, S_w_0, 1.d0, 0.d0,1.d0)
            H_frac_soot = dep_frac + (1. - dep_frac)*aux
            if (H_frac_soot > 1.) H_frac_soot = 1.
            
            dH_frac_soot = dfdep + (SIW*(1. - dep_frac)*dH1smooth)- aux*dfdep

            if ((H_frac_soot .gt. 1.0e-6) .and. (H_frac_soot .lt. 1.)) then 
               dH_frac_soot = dH_frac_soot/H_frac_soot
            else
               dH_frac_soot = 0.0
            end if 
            
            !---- biological particles --------
            H_frac_bio = H_frac_soot 
            
            !SB: start "IF-2" ----------------------------------------------------------------- 
            !if (temperature_K < 273.15 .and. temperature_K >= 273.15-35.) then   !SB: 238.15 <= T < 273.15, MIXED-phase regime
            if (temperature_K < 273.15 .and. temperature_K > 273.15-35.) then     !SB_change: 238.15 < T < 273.15, MIXED-phase regime

               !SB: ---- average number concentration of IN for a background-troposphere scenario (n_in* in PDA08)
               n_in = 1000.*(exp(12.96*SS_i - 0.639)/RHO_CFDC)*0.0587*FAC_CORRECT_RH    !SB: similar to eq(3) PDA08

               if (temperature_K > 273.15-5. .and. temperature_K < 273.15-2. ) then
                  n_in = n_in*H_1_smooth(-(temperature_K-273.15), 2.d0, 5.d0, 0.d0, 1.d0)
               end if

               if (temperature_K >= 273.15-2. ) n_in = 0.

               if (temperature_K < 273.15-25. ) then
                  n_in_tilde = 1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC           !SB: Table in appendix-A PDA08
                  n_in_hat   = n_in
                  
                  if (temperature_K >= 273.15-30.) n_in_max  = 1000.*(exp(12.96*SS_iw - 0.639)/RHO_CFDC)*0.0587*FAC_CORRECT_RH    !SB: similar to eq(7) PDA08
                  if (temperature_K < 273.15-30.)  n_in_max  = 1000.*(exp(0.1296*(SS_iw*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC  !SB: eq(6) PDA08
                  if (n_in_hat > n_in_max)         n_in_hat  = n_in_max                             !SB: take the min
                  if (n_in_tilde > n_in_max)       n_in_tilde= n_in_max                             !SB: take the min

                  n_in = n_in_hat*((n_in_tilde/n_in_hat)**(H_1_smooth(-(temperature_K-273.15), 25.d0, 35.d0, 0.d0, 1.d0)))        !SB: eq(5) PDA08
                  if (n_in > n_in_max)             n_in = n_in_max                                  !SB: take the min, eq(4) PDA08
               end if   

               !------ dust ---------------------------------
               n_in_dust = 0.
               dn_in_dust = 0.
               
               if (temperature_K < 273.15 - 30.) then    !DONIF
                  dnaux = 3.88                           !this is a simplified derivative of dNds
               else
                  dnaux = 12.96
               end if 
               naux=0.0
               
               DO ij = 1, ijstop_dust
                  mu = n_in*ALPHA_DUST*H_frac_dust*PIE*(D_grid_dust(ij)**2.) &                          !SB: eq(10) PDA08
                       /BASE_DUST_OMEGA
                  naux = (1. - exp(-mu))*n_grid_dust(ij)                                                !SB: eq(9) PDA08
                  n_in_dust  = n_in_dust + naux                                                         !SB: eq(9) PDA08
                  dn_in_dust = max(mu*(n_grid_dust(ij)-naux)*(dnaux + dH_frac_dust), 0.0) + dn_in_dust  !SB: eq(13) PDA08
               END DO

               if (temperature_K > 273.15+TEMP_MAX_DUST_DEGC-20. .and. temperature_K < 273.15+TEMP_MAX_DUST_DEGC) then
                  n_in_dust = n_in_dust*H_1_smooth(-(temperature_K-273.15),-TEMP_MAX_DUST_DEGC,-TEMP_MAX_DUST_DEGC+20., 0.d0, 1.d0)
               end if

               if (temperature_K >= 273.15+TEMP_MAX_DUST_DEGC) n_in_dust = 0.
               
               !------ soot ----------------------------------
               n_in_soot  = 0.
               dn_in_soot = 0.

               DO ij = 1, ijstop_soot
                  mu = n_in*ALPHA_SOOT*H_frac_soot*PIE*(D_grid_soot(ij)**2.) &
                         /BASE_SOOT_PHILIC_OMEGA
                  naux = (1. - exp(-mu))*n_grid_soot(ij)
                  n_in_soot  = n_in_soot + naux 
                  dn_in_soot = max(mu*(n_grid_soot(ij)-naux)*(dnaux+dH_frac_soot), 0.0) + dn_in_soot
               END DO

               if (temperature_K > 273.15+TEMP_MAX_SOOT_DEGC-10. .and. temperature_K < 273.15+TEMP_MAX_SOOT_DEGC) then
                  n_in_soot = n_in_soot*H_1_smooth(-(temperature_K-273.15),-TEMP_MAX_SOOT_DEGC,-TEMP_MAX_SOOT_DEGC+10., 0.d0, 1.d0)
               end if

               if (temperature_K >= 273.15 + TEMP_MAX_SOOT_DEGC) n_in_soot = 0.

               !------ biological aerosols --------------------
               n_in_bio  = 0.
               dn_in_bio = 0.

               DO ij = 1, ijstop_bio
                  mu = n_in*ALPHA_bio*H_frac_bio*PIE*(D_grid_bio(ij)**2.) &
                       /BASE_BIO_OMEGA
                  mu = n_in*ALPHA_bio*H_frac_bio
                  naux = (1. - exp(-mu))*n_grid_bio(ij)
                  n_in_bio  = n_in_bio + naux 
                  dn_in_bio = max(mu*(n_grid_bio(ij)-naux)*dnaux, 0.0) + dn_in_bio
               END DO
              
               if (temperature_K > 273.15+TEMP_MAX_bio_DEGC-3. .and. temperature_K < 273.15+TEMP_MAX_bio_DEGC) then
                  n_in_bio = n_in_bio*H_1_smooth(-(temperature_K-273.15),-TEMP_MAX_bio_DEGC,-TEMP_MAX_bio_DEGC+3., 0.d0, 1.d0)
               end if

               if(temperature_K >= 273.15 + TEMP_MAX_bio_DEGC ) n_in_bio = 0.

            else    !SB: 180.15 < T < 238.15, cirrus regime 
               n_in       = 0.
               n_in_ultra = 0.
               n_in_dust  = 0.
               n_in_soot  = 0.
               n_in_bio   = 0.
            end if  
            !SB: start "IF-2" ----------------------------------------------------------------- 

            !SB: start "IF-3" ----------------------------------------------------------------- 
            !if (temperature_K < 273.15-35.) then         !SB: 180.15 < T < 238.15, CIRRUS regime                          
            if (temperature_K <= 273.15-35.) then        !SB_change: 180.15 < T <= 238.15, CIRRUS regime                          

               n_in_ultra = 1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC           !SB: equal to n_in_tilde 
               dnaux = 3.88            !DONIF simplified treatment of derivative asusming dH small 
               naux  = 0.0

               !------ dust ----------------------------------
               n_in_dust_ultra = 0. 
               dn_in_dust = 0.0

               DO ij = 1, ijstop_dust
                  mu = n_in_ultra*ALPHA_DUST*H_frac_dust*PIE*(D_grid_dust(ij)**2.) &
                       /BASE_DUST_OMEGA
                  naux = (1. - exp(-mu))*n_grid_dust(ij)
                  n_in_dust_ultra = n_in_dust_ultra + naux 
                  dn_in_dust      = max(mu*(n_grid_dust(ij)-naux)*(dnaux +dH_frac_dust), 0.0) + dn_in_dust
               END DO 

               !------ soot ----------------------------------
               n_in_soot_ultra = 0.0
               dn_in_soot      = 0.0

               DO ij = 1, ijstop_soot
                  mu = n_in_ultra*ALPHA_SOOT*H_frac_soot*PIE*(D_grid_soot(ij)**2.) &
                       /BASE_SOOT_PHILIC_OMEGA
                  naux = (1. - exp(-mu))*n_grid_soot(ij)
                  n_in_soot_ultra = n_in_soot_ultra + naux 
                  dn_in_soot      = max(mu*(n_grid_soot(ij)-naux)*(dnaux +dH_frac_soot), 0.0) + dn_in_soot
               END DO
               
               !------ biological aerosols --------------------
               n_in_bio_ultra = 0.0
               dn_in_bio      = 0.0

               DO ij = 1, ijstop_bio
                  mu = n_in_ultra*ALPHA_bio*H_frac_bio*PIE*(D_grid_bio(ij)**2.) &
                       /BASE_BIO_OMEGA
                  naux = (1. - exp(-mu))*n_grid_bio(ij)
                  n_in_bio_ultra = n_in_bio_ultra + naux 
                  dn_in_bio      = max(mu*(n_grid_bio(ij)-naux)*dnaux, 0.0) + dn_in_bio
               END DO

               !------ soluble organics -----------------------
               rho_AIDA = 90000./(287.*205.)                                       !SB: air density in the AIDA chamber [kg/m3]
               Psi_solO = A_solO/BASE_SOLO_OMEGA                                   !SB: fraction of omegas used in eq(8) PDA13
               RHI   = (SS_i+1.)*100.
               if (RHI < 0.) RHI = 0.
               n_in_solO_star = 1000.e6*(7.7211e-5 * RHI - 9.2688e-3)/rho_AIDA     !SB: used in eq(8) PDA13 [1/kg] 
                                                                                   
               n_in_solO  = Psi_solO*glass_frac*H_frac_solO*n_in_solO_star         !SB: eq(8) PDA13
               dn_in_solO = max(Psi_solO*glass_frac*&
                            (H_frac_solO*77211.0*100.0/rho_AIDA + n_in_solO_star*dH_frac_solo), 0.0)

            else        !SB: if T > (273.15 - 35.), mixed-phase regime
               n_in_ultra      = 0.
               n_in_dust_ultra = 0.
               n_in_soot_ultra = 0.
               n_in_bio_ultra  = 0. 
               n_in_solO       = 0.  
            end if  
            !SB: end "IF-3" ----------------------------------------------------------------- 
       
            !SB: contributions at (238.15 <= T < 273.15)+(180.15 < T < 238.15)
            n_in_dust = n_in_dust + n_in_dust_ultra      
            n_in_soot = n_in_soot + n_in_soot_ultra
            n_in_bio  = n_in_bio  + n_in_bio_ultra

            if (n_in_dust + n_in_bio + n_in_soot + n_in_solO > 0.) then
             
               !SB_changhe: below...
               CIHENC_dust = n_in_dust
               !CIHENC_dust = n_in_dust - nin_a_nuc_dust 
               if(CIHENC_dust < 0.) CIHENC_dust = 0.
               
               CIHENC_soot = n_in_soot
               !CIHENC_soot = n_in_soot - nin_a_nuc_soot 
               if(CIHENC_soot < 0.) CIHENC_soot = 0.
               
               CIHENC_bio = n_in_bio
               !CIHENC_bio = n_in_bio - nin_a_nuc_bio 
               if(CIHENC_bio < 0.) CIHENC_bio = 0.
               
               CIHENC_solO = n_in_solO
               !CIHENC_solO = n_in_solO - nin_a_nuc_solO 
               if(CIHENC_solO < 0.) CIHENC_solO = 0.
               
               n_iw            = n_iw + CIHENC_dust                !SB: n_iw = n_iw + n_in_dust
               !nin_a_nuc_dust  = nin_a_nuc_dust + CIHENC_dust
               !num_ic_dust_imm = num_ic_dust_imm + CIHENC_dust
               
               n_iw            = n_iw + CIHENC_soot                !SB: n_iw = n_iw + n_in_soot 
               !nin_a_nuc_soot  = nin_a_nuc_soot + CIHENC_soot
               !num_ic_soot_imm = num_ic_soot_imm + CIHENC_soot
               
               n_iw            = n_iw + CIHENC_bio                 !SB: n_iw = n_iw + n_in_bio
               !nin_a_nuc_bio   = nin_a_nuc_bio + CIHENC_bio
               !num_ic_bio_imm  = num_ic_bio_imm + CIHENC_bio
               
               n_iw            = n_iw + CIHENC_solO                !SB: n_iw = n_iw + n_in_solO = sum over all species (=OUTPUT var)
               !nin_a_nuc_solO  = nin_a_nuc_solO + CIHENC_solO
               !num_ic_solO_imm = num_ic_solO_imm + CIHENC_solO
            end if

         endif   !SB: 180.15 < T < 273.15
      endif      !SB: SS_i > 0

      !dNall = dn_in_dust + dn_in_bio + n_in_soot + n_in_solO  !DONIF   !mz_sb_20170127
      dNall = dn_in_dust + dn_in_bio + dn_in_soot + dn_in_solO  !DONIF    !mz_sb_20170127
       
      if (( dNall > 0.) .and. (n_iw .gt. 0.0))   then
          Dsh = max(min(n_iw/dNall, SS_i), 0.005)
      else
          Dsh = 0.005
      end if   

 END SUBROUTINE EMPIRICAL_PARAM_PHILLIPS

!*************************************************************
! FUNCTION: real function (cubic interpolation) 
!*************************************************************

      real(dp) FUNCTION H_1_smooth(X, X_1, X_2, Hlo, Hhi)     

      real(dp), intent(in) :: Hlo, Hhi, X, X_1, X_2          
      real(dp) :: a_0, a_1, a_2, a_3, A, B

      if (X >= X_2) H_1_smooth = Hhi
      if (X <= X_1) H_1_smooth = Hlo 

      if (X >= X_2) dH1smooth = 0.0
      if (X <= X_1) dH1smooth = 0.0 

      if (X > X_1 .and. X < X_2) then
         A = 6.*(Hlo - Hhi)/(X_2**3. - X_1**3. + 3.*(X_2*X_1*X_1 - X_1*X_2*X_2) )
         a_3 = (A/3.)
         a_2 = -(A/2.)*(X_1 + X_2)
         a_1 = A*X_2*X_1
         B = Hlo + A*(X_1**3.)/6. - A*X_1*X_1*X_2/2.
         a_0 = B
         H_1_smooth = a_0 + a_1*X + a_2*X*X + a_3*X*X*X
         dH1smooth =  a_1 + 2.0*a_2*X + 3.0*a_3*X*X
      end if

      dH1smooth = min(dH1smooth , 1.0e6)
      dH1smooth = max(dH1smooth , zero_par)

      if (X_2 <= X_1) stop 91919

      return 
      END FUNCTION 


!************************************************************************
!
!    END ICE PARAMETERIZATION DONIF
!
!************************************************************************


 END MODULE messy_cloud_ice_BN09
