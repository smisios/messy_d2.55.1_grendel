! ***************************************************************************
MODULE  messy_rad_fubrad
! ***************************************************************************
  !     
  ! PURPOSE:  Calculate heating rates in the middle atmosphere for the
  ! --------  UV-Vis part of the solar spectrum.
  !
  ! Sub-Submodel to RAD, to increase the resolution in the UV-Vis
  ! part of the solar spectrum.
  ! At pressures higher than 70 hPa the UV-Vis shortwave radiation fluxes
  ! are calculated in one spectral interval as in the original ECHAM5 code.
  ! At pressures equal or lower than 70 hPa the UV-Vis shortwave radiation
  ! fluxes due to ozone and oxygen are calculated using the code of the
  ! Freie Universitaet Berlin (FUB).
  ! Between full radiation time steps the FUB code is called by 
  ! rad_fubrad_radheat, thus the middle atmosphere UV-Vis fluxes are
  ! updated every time step while tropospheric UV-Vis fluxes get updated
  ! only on full radiation timesteps.
  !
  ! Official FUBRAD resolutions:
  ! ----------------------------
  ! 55 bands (nchap = 6, nsrb = 1):
  !   -> namelist entries:
  !      nbands    =  55,   ! this is default  
  !      sr % type = 'STR', ! this is default
  ! 106 bands (nchap = 57, nsrb = 1):
  !   -> namelist entries:
  !      nbands    = 106,   
  !      sr % type = 'STR', ! this is default  
  ! 81 bands (nchap = 14, nsrb = 19):
  !   -> namelist entries:  
  !      nbands    =  81, 
  !      sr % type = 'STR-HR22'
  !
  ! AUTHORS:
  ! --------
  ! Ulrike Langematz, Freie Universitaet Berlin, Berlin Germany
  ! Steven Pawson, FUB (now at NASA Goddard Space Flight Center)
  ! Katja Matthes, Freie Universitaet Berlin, Berlin Germany
  ! Katrin Nissen, Freie Universitaet Berlin, Berlin Germany
  ! Markus Kunze, Freie Universitaet Berlin, Berlin Germany
  ! Klaus Ketelsen 
  ! Patrick Joeckel, MPIC, DLR
  ! 
  ! SUBROUTINE o3fluxes:   originally implemented by K.P. Shine, Oxford, 1988
  ! SUBROUTINE lymanalpha: originally implemented by Sasha Madronich
  !
  ! REFERENCES:
  ! -----------
  ! Nissen et al., ACP, 7, 5391-5400, doi:10.5194/acp-7-5391-2007, 2007
  ! Matthes et al., JGR, 109, D06101, doi:10.1029/2003JD004012, 2004
  ! Shine and Rickaby, Atmospheric Ozone, pp 597-600, 1990
  ! Strobel, JGR, 83, C12, pp 6225-6230, 1978, (Schumann-Runge)
  ! Chabrillat and Kockarts, GRL, 24, 21, pp 2659-2662, 1997 (Lyman-alpha)
  ! Mlynczak and Solomon, JGR, 98, p 10517, 1993 (efficiency factors)
  !
  USE messy_main_constants_mem, ONLY: dp, g, cpd => cp_air, N_A
  USE messy_rad_fubrad_mem
  USE messy_rad_fubrad_srb_km,  ONLY: fubrad_srb_km_schum
  USE messy_rad_fubrad_srb_kck, ONLY: fubrad_srb_kck_schu
  USE messy_rad_fubrad_sr_str,  ONLY: schumann_runge           &
                                    , src_hr_o2                &
                                    , src_eff                  &
                                    , srb_hr_zhu               &
                                    , schuru_hr_o2             &
                                    , schuru_hr_o2_srb         &
                                    , schuru_flx_o2            &
                                    , src_flx_o2
#if defined (NAGFOR)
  USE F90_UNIX_PROC, ONLY: exit
#endif
  
  IMPLICIT NONE
  PRIVATE
  SAVE
  
  PUBLIC :: prepare_radiation
  PUBLIC :: middle_atmosphere_downward_flux
  PUBLIC :: middle_atmosphere_upward_flux
  PUBLIC :: middle_atmosphere_heat_rates
  
  ! Interface of public subroutines
  
  INTERFACE prepare_radiation
     MODULE PROCEDURE prepare_radiation
  END INTERFACE
  
  INTERFACE middle_atmosphere_downward_flux
     MODULE PROCEDURE middle_atmosphere_downward_flux
  END INTERFACE
  
  INTERFACE middle_atmosphere_heat_rates
     MODULE PROCEDURE middle_atmosphere_heat_rates
  END INTERFACE
  
  INTERFACE middle_atmosphere_upward_flux
     MODULE PROCEDURE middle_atmosphere_upward_flux
     MODULE PROCEDURE middle_atmosphere_upward_flux_h
  END INTERFACE

  ! Interface of module private subroutines

  INTERFACE ox_pathlength
     MODULE PROCEDURE ox_pathlength
  END INTERFACE
  
  INTERFACE o3fluxes
     MODULE PROCEDURE o3fluxes
  END INTERFACE
  
CONTAINS

  SUBROUTINE prepare_radiation (kproma, klev, o3, o2, ppf, pph, zmu0)

    INTEGER,  INTENT(in)                          :: kproma,klev
    REAL(dp), INTENT(in),    DIMENSION(:,:)       :: o3
    REAL(dp), INTENT(in),    DIMENSION(:,:)       :: o2
    REAL(dp), INTENT(in),    DIMENSION(:,:)       :: ppf
    REAL(dp), INTENT(in),    DIMENSION(:,:)       :: pph
    REAL(dp), INTENT(in),    DIMENSION(:)         :: zmu0

    !-  local variables
    INTEGER :: kbdim

    kbdim = SIZE(po3c,1)
    
    ! init for kproma < nbdim
    prmu0(:)           = 0.0_dp
    prmu0(1:kproma)    = zmu0(1:kproma)
    my_pph(:,:)        = 0.0_dp
    my_pph(1:kproma,:) = pph(1:kproma,:)

    my_ppf(:,:)        = 0.0_dp
    my_ppf(1:kproma,:) = ppf(1:kproma,:)
    przsec(:)          = 0.0_dp 
    przsec(1:kproma)   = 35._dp/(DSQRT(1224._dp*(zmu0(1:kproma)*zmu0(1:kproma))+1._dp))
    
    CALL ox_pathlength(o3, po3c(:,:),'O3',kproma,kbdim,klev,nswlev,pph)
    CALL ox_pathlength(o2, po2c(:,:),'O2',kproma,kbdim,klev,nswlev,pph)
    !
    RETURN
  END SUBROUTINE prepare_radiation

  SUBROUTINE middle_atmosphere_downward_flux (kproma, kbdim, klev, tm1, iswlev, pffrac)
    !
    ! Input:
    ! ------
    ! po2    -  oxygen in mass mixing ratio
    ! pph    -  pressure on half pressure levels
    !
    INTEGER,  INTENT(in)                           :: kproma, kbdim, klev
    REAL(dp), INTENT(in), DIMENSION(:,:), OPTIONAL :: tm1
    !
    ! Output:
    REAL(dp), INTENT(out), DIMENSION(:), OPTIONAL :: pffrac
    INTEGER,  INTENT(out)                         :: iswlev


    ! Interfaces of the subroutines o3fluxes and lymanalpha have changed:
    !  - additional arguments of o3fluxes:   po2c, flhz, flhudo, flchdo;
    !  - removed    arguments of o3fluxes:   prmu0, my_pph, klev+1
    !  - removed     argument of lymanalpha: prmu0
    CALL o3fluxes(kbdim, kproma, nswlev+1, 1, po3c, po2c, &
                  altsw, fldo, flhart, flhz, flhudo, flchdo)
    CALL lymanalpha(kbdim, kproma, nswlev+1, fllya(1:kbdim,:))
    SELECT CASE(TRIM(sr%type))
    CASE('STR-FLX','STROBEL-FLX') !  NBANDS = 146
       CALL schuru_flx_o2(kbdim, kproma, nswlev+1, klev+1, flsrb, flsrc, po2c)
       lsrflx = .TRUE.
    CASE('KCK','KOCKARTS')        !  NBANDS = 143
       CALL src_flx_o2(kbdim, kproma, nswlev+1, klev+1, flsrc, po2c)
       CALL fubrad_srb_kck_schu( kbdim, kproma, nswlev, nswlev+1, nsrb, flsrb)
       lsrflx = .TRUE.
    CASE('KM','KOPPERS-MURTAGH')  !  NBANDS = 144
       CALL src_flx_o2(kbdim, kproma, nswlev+1, klev+1, flsrc, po2c)
       CALL fubrad_srb_km_schum(kbdim, kproma, nswlev, nswlev+1, tm1, flsrb)
       lsrflx = .TRUE.
    CASE DEFAULT
       
    END SELECT
    ! 1.3  DEPLETED BEAM AS INPUT FOR MORCRETTE UPPER BOUNDARY
    ! ---------------------------------------------------
    ! Includes contribution of 9 WMO intervals that are not absorbed
    !
    IF (PRESENT(pffrac)) THEN
       !
       ! The short wave flux of 9 WMO intervals (362.5 - 407.5 nm) 
       ! is not absorbed/ changed by the FUBRad scheme: 54.25_dp
       ! pffrac(jl) = (fldo(jl,nswlev+1)+  flhart(jl,nswlev+1) + 54.25_dp)/&
       !        &     (fldo(jl,1)       +  flhart(jl,1)        + 54.25_dp)
       ! The missing flux is now included in fldo; added in subroutine o3fluxes.
       !
       pffrac(1:kproma) = ( fldo(1:kproma,nswlev+1) + flhart(1:kproma,nswlev+1) + &
                           fllya(1:kproma,nswlev+1) ) / &
                          ( fldo(1:kproma,1) + flhart(1:kproma,1) + &
                           fllya(1:kproma,1) )
    END IF
    !
    iswlev = nswlev
    !
    RETURN
  END SUBROUTINE middle_atmosphere_downward_flux

  SUBROUTINE middle_atmosphere_upward_flux (kproma,kbdim,klev,pfu,pfd, &
                                            pfuc, pfdc)
    !
    INTEGER,  INTENT(in)                    :: kproma, kbdim, klev
    REAL(dp), INTENT(inout), DIMENSION(:,:) :: pfu,pfd
    REAL(dp), INTENT(inout), DIMENSION(:,:) :: pfuc ! upward flux clear sky
    REAL(dp), INTENT(inout), DIMENSION(:,:) :: pfdc ! downward flux clear sky

    !-  Local
    INTEGER :: kfjkf
    INTEGER :: jk
    
    kfjkf = klev - nswlev + 1
    !
    !   Tropospheric Albedo at UV and visible wavelengths for o3fluxes upward flux
    !
    altsw(1:kproma) = pfu(1:kproma, kfjkf) / pfd(1:kproma, kfjkf)
    !
    ! Tropospheric albedo at UV and VIS wavelength for o3fluxes upward flux
    ! clear sky conditions.
    altswc(1:kproma) = pfuc(1:kproma, kfjkf) / pfdc(1:kproma, kfjkf)
    !
    !   SET FLUXES CONSTANT AT FUB-RADIATION LEVELS
    !   (NO HEATING BY MORCRETTE AT SHORT-VIS WAVELENGTHS)
    DO jk = kfjkf,klev+1
       pfu  (1:kproma,jk) = pfu(1:kproma,kfjkf)
       pfd  (1:kproma,jk) = pfd(1:kproma,kfjkf)
       pfuc (1:kproma,jk) = pfuc(1:kproma,kfjkf)
       pfdc (1:kproma,jk) = pfdc(1:kproma,kfjkf)
    END DO
    !
    !   Calculate middle atmosphere upward flux
    !
    CALL o3fluxes(kbdim, kproma, nswlev+1, 2          &
                 , po3c, po2c, altsw, flup, flhart    &
                 , flhz, flhuup, flchup               &
                 , flupc)
    !
    RETURN
  END SUBROUTINE middle_atmosphere_upward_flux
  
  SUBROUTINE middle_atmosphere_upward_flux_h (kproma, kbdim, klev)
    !
    INTEGER,  INTENT(in)                    :: kproma, kbdim, klev
    !
    !   Calculate middle atmosphere upward flux
    !
    CALL o3fluxes(kbdim, kproma, nswlev+1, 2         &
                 , po3c, po2c, altsw, flup, flhart   &
                 , flhz, flhuup, flchup              &
                 , flupc)
    !
    RETURN
  END SUBROUTINE middle_atmosphere_upward_flux_h
  !
  ! =====================================================================
  !
  SUBROUTINE middle_atmosphere_heat_rates (kproma, nbdim, klev,        &
                                           paphm1, papm1, pmu0, pdayl, &
                                           po2,                        &
                                           ptte,                       &
                                           pheato3flux,                &
                                           pheatsrc, pheatsrb,         &
                                           pheatlya,                   &
                                           pheatmosw,                  &
                                           pheatmoswc,                 &
                                           pheatsw,                    &
                                           pheatswc,                   &
                                           pheatherz, pheathart,       &
                                           pheathug,  pheatchap,       &
                                           psflux0, delta_time)
    ! INPUT:
    ! ------
    INTEGER,  INTENT(in)                 :: kproma, nbdim, klev
    REAL(dp), INTENT(in), DIMENSION(:,:) :: paphm1, papm1
    REAL(dp), INTENT(in), DIMENSION(:)   :: pmu0, pdayl
    REAL(dp), INTENT(in), DIMENSION(:,:) :: po2
    REAL(dp), INTENT(in), DIMENSION(:,:) :: pheatmosw
    REAL(dp), INTENT(in), DIMENSION(:,:) :: pheatmoswc
    REAL(dp), INTENT(in) :: delta_time
    !
    ! OUTPUT:
    ! -------
    ! ptte   - temperature tendency
    !
    ! diagnostic output:
    ! ------------------
    ! pheato3flux - heating rates: due to fluxes from O3fluxes
    ! pheatsrc    - heating rates: Schumann-Runge Continuum
    ! pheatsrb    - heating rates: Schumann-Runge Bands
    ! pheatlya    - heating rates: Lyman alpha
    ! pheatsw     - heating rates: all SW (FUBRad + FB1980)
    ! pheatswc    - heating rates: all SW clear sky (FUBRad + FB1980)
    ! pheatherz   - heating rates: Herzberg
    ! pheathart   - heating rates: Hartley
    ! pheathug    - heating rates: Huggins
    ! pheatchap   - heating rates: Chappuis
    ! psflux0     - shortwave fluxes FUBRad TOA
    !

    REAL(dp), INTENT(out),   DIMENSION(:,:) :: ptte
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheato3flux
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheatsrc
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheatsrb
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheatlya
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheatsw
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheatswc

    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheatherz
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheathart
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheathug
    REAL(dp), INTENT(out),   DIMENSION(:,:) :: pheatchap

    REAL(dp), INTENT(inout), DIMENSION(:)   :: psflux0
    !
    !-  local
    !REAL(dp), PARAMETER :: oxfac=1.e2 / cpd  ! 1.e2= cm-1 -> m-1
    REAL(dp), PARAMETER :: zcons3=g/cpd
    INTEGER             :: jk
    !
    REAL(dp),DIMENSION(nbdim,nswlev) :: psrbht, psrcht
    REAL(dp),DIMENSION(nbdim)        :: zfact, zfact_zdp
    REAL(dp),DIMENSION(nbdim)        :: zfact_sr

    ! BECAUSE OF THE ATTRIBUTE 'INTENT(OUT)' THE FOLLOWING FIELDS NEED
    ! TO BE INITIALISED (for 'odd' kproma with 'rest') ...
    pheato3flux(:,:) = 0.0_dp
    pheatsrc(:,:)    = 0.0_dp
    pheatsrb(:,:)    = 0.0_dp
    pheatlya(:,:)    = 0.0_dp
    pheatsw(:,:)     = 0.0_dp
    pheatswc(:,:)    = 0.0_dp
    ptte(:,:)        = 0.0_dp

    pheatherz(:,:)   = 0.0_dp
    pheathart(:,:)   = 0.0_dp
    pheathug(:,:)    = 0.0_dp
    pheatchap(:,:)   = 0.0_dp
    !
    zfact(:) = cdisse_fubrad * pmu0(:) * pdayl(:)

    ! Bugfix:
    ! To calculate the final heating rates in the Schumann-Runge bands
    ! and Schumann-Runge continuum, the factor pmu0 (1./cos(sza)) does
    ! not have to be included in the factor zfact, when the flux divergence 
    ! is the result of the parametrization of the heating rates in these 
    ! spectral regions. This is the case for the Schumann-Runge bands and
    ! Schumann-Runge continuum in some configurations (lsrflx=.FALSE.).
    zfact_sr(:) = cdisse_fubrad * pdayl(:)
    !
    !*  Computation of O2 solar heating rates
    !   Schumann-Runge bands and continuum (Strobel)
    !CDIR NOIEXPAND

    SELECT CASE(TRIM(sr%type))
    CASE('STR','STROBEL')    
       ! original parameterization (3 SRC, 1 SRB), NBANDS = 55, 106
       CALL schumann_runge(nbdim, kproma, nswlev, klev, psrbht, psrcht, PO2, PAPM1)
       lsrflx = .FALSE.
    CASE('STR-HR22','STROBEL-HR22') 
       ! hr parameterization with 3 SRC, 19 SRB bands, NBANDS = 81, 124
       CALL schuru_hr_o2_srb  (nbdim, kproma, nswlev, klev, psrbht, PO2, PAPM1)
       CALL src_hr_o2         (nbdim, kproma, nswlev, klev, psrcht, PO2, PAPM1)
          
       lsrflx = .FALSE.
    CASE('STR-HR44','STROBEL-HR44') 
       ! hr parameterization with 25 SRC, 19 SRB bands, NBANDS = 146
       CALL schuru_hr_o2  (nbdim, kproma, nswlev, klev, psrbht, psrcht, PO2, PAPM1)
       lsrflx = .FALSE.
    CASE('ZHU')        
       ! Strobel parameterization (3 SRC), Zhu (1 SRB), NBANDS = 55, 106
       CALL srb_hr_zhu    (nbdim, kproma, nswlev, klev, po2, papm1, psrbht, psrcht)
       lsrflx = .FALSE.
    END SELECT
    !
    DO jk = 1, nswlev

       zfact_zdp(1:kproma) = zfact(1:kproma) * zcons3 / &
                             ( paphm1(1:kproma,jk+1)- paphm1(1:kproma,jk) )
       !
       ! Short-wave heating rates for the Hartley bands:
       ! treated separately from other O3 fluxes because
       ! of the multiplication with the efficiency factor.
       !
       pheathart(1:kproma,jk) =                             &
               HARTLEY_EFF((papm1(1:kproma,jk)/100.))*      &
                 ( flhart(1:kproma,jk)  - flhart(1:kproma,jk+1)) * &
                 zfact_zdp(1:kproma)
       !
       ! Combined short-wave heating rates for all spectral 
       ! intervals in O3fluxes.
       !
       pheato3flux(1:kproma,jk) = pheathart(1:kproma,jk) + &
              ( fldo(1:kproma,jk) - fldo(1:kproma,jk+1)  &
            +   flup(1:kproma,jk) - flup(1:kproma,jk+1)) * &
                zfact_zdp(1:kproma)
       ! heating rates for Herzberg only
       pheatherz(1:kproma,jk) =  &
              ( flhz(1:kproma,jk) - flhz(1:kproma,jk+1)) * &
                zfact_zdp(1:kproma)
       ! heating rates for Huggins only
       pheathug(1:kproma,jk) =  &
              ( flhudo(1:kproma,jk) - flhudo(1:kproma,jk+1)  &
            +   flhuup(1:kproma,jk) - flhuup(1:kproma,jk+1)) * &
                zfact_zdp(1:kproma)
       ! heating rates for Chappuis only
       pheatchap(1:kproma,jk) =  &
              ( flchdo(1:kproma,jk) - flchdo(1:kproma,jk+1)  &
            +   flchup(1:kproma,jk) - flchup(1:kproma,jk+1)) * &
                zfact_zdp(1:kproma)
       !
       ! Schumann-Runge continuum and bands:
       !
       IF (lsrflx) THEN
          pheatsrc(1:kproma,jk) =                             &
               SRC_EFF(papm1(1:kproma,jk)/100.)*      &
                 ( flsrc(1:kproma,jk)   - flsrc(1:kproma,jk+1)) * &
                   zfact_zdp(1:kproma)
          pheatsrb(1:kproma,jk) =  &
              ( flsrb(1:kproma,jk)    - flsrb(1:kproma,jk+1)) * &
                zfact_zdp(1:kproma)
       ELSE
         pheatsrc(1:kproma,jk) = zfact_sr(1:kproma) * psrcht(1:kproma,jk)
         pheatsrb(1:kproma,jk) = zfact_sr(1:kproma) * psrbht(1:kproma,jk)
       END IF
       !
       ! Lyman-alpha heating rates
       !
       pheatlya(1:kproma,jk) = lya_eff *  &
              ( fllya(1:kproma,jk) - fllya(1:kproma,jk+1)) * &
                zfact_zdp(1:kproma)
                
       !
       ! Total SW heating rates clear sky
       !
       pheatswc(1:kproma,jk) = pheatsrb(1:kproma,jk)  + &
                               pheatsrc(1:kproma,jk)  + &
                               pheatlya(1:kproma,jk)  + &
                               pheathart(1:kproma,jk) + &
            ( fldo(1:kproma,jk)  - fldo(1:kproma,jk+1)  &
            + flupc(1:kproma,jk) - flupc(1:kproma,jk+1) ) * &
              zfact_zdp(1:kproma)
      
    END DO
    !
    ! Combined temperature tendency
    !
    ptte(1:kproma,1:nswlev) =  pheato3flux(1:kproma,1:nswlev) &
         + pheatsrb(1:kproma,1:nswlev) + pheatsrc(1:kproma,1:nswlev) &
         + pheatlya(1:kproma,1:nswlev)
    !
    !   Diagnostics
    !   tendency diagnostics
    !
    pheatsw(1:kproma,:) = ptte(1:kproma,:) + &
                          pheatmosw(1:kproma,:)
    
    pheatswc(1:kproma,:) = pheatswc(1:kproma,:) + &
                           pheatmoswc(1:kproma,:)
    !
    ! TOP OF ATMOSPHERE FLUX
    !

    !NOTE: With new memory channel interface, laccu (requiring
    !       'user'-accumulation) is obsolete.

    psflux0(1:kproma) = (fllya(1:kproma,1)    + &
                         fldo(1:kproma,1)     + &
                         flup(1:kproma,1)     + &
                         flhart(1:kproma,1)   + &
                         3.01E-3_dp) * zfact(1:kproma)
    
    !
    RETURN
  END SUBROUTINE middle_atmosphere_heat_rates
  !
  ! =====================================================================
  !
  ! Private subroutines
  ! copied from messy_rad_short.f90
  
  SUBROUTINE o3fluxes (kbdim, kproma, kswlevp1, imode, po3c, po2c, &
                       paltsw, zflux, zhart, zherz, zhug, zchap,   &
                       pfluxc)
    !
    !   AUTHOR
    !   ------
    !   K.P. Shine,  Oxford,        1988.
    !
    !   CHANGES
    !   -------
    !   S. Pawson,   Berlin,        1989/1990.                               
    !    - treatment of a single latitude line.
    !    - solar geometry uses *ECMWF* formulation.
    !    - vectorised as much as possible.
    !    - fluxes are returned: downward fluxes when imode = 1
    !                           upward   fluxes when imode = 2
    !      at the *kswlevp1* half levels of the model.
    !                               
    !   K. Matthes,  Berlin, 2003
    !    - additional spectral intervals
    !
    !   K. Nissen,   Berlin, 2005
    !    - conversion to Fortran90
    !    - changes for ECHAM5/Messy
    !    - hartley bands fluxes have to be stored seperately
    !      as hartley band heating rates will be  multiplied
    !      by efficiency factor 
    !
    !   K. Matthes, Berlin, 2008
    !     - new fluxes for solar max and min conditions (average over 2 solar
    !       cycles representative for solar max and solar min)
    !
    !   M. Kunze, FU-Berlin, December 2011
    !    - added the integrated flux of non absorbing Chappuis interval
    !      362.5-407.5 nm (WMO 61-69) to the upward and downward flux,
    !      instead of using the constant flux 54.25, when calculating pffrac.
    !    - substituted the parametrized flux of the Chappuis band by the
    !      57 spectral intervals from 407.5-690.0 nm (WMO 70-126) to get
    !      the correct shortwave flux in the middle atmosphere.
    !
    !   M. Kunze, FU-Berlin, March 2012
    !    - the number of spectral intervals in the Chappuis band is variable.
    !    - nchap is set after evaluating the namelist parameter nbands.
    !    - the fluxes for solar max. and min. are initialised in fubrad_initialize_fluxes
    !      as are the absorption cross sections for O3 and O2.
    !
    !   PURPOSE
    !   -------
    !   Solar fluxes at the model half levels are calculated for the middle
    !   atmosphere. Absorption due to ozone is treated.
    !
    !   METHOD
    !   ------
    !   Radiative fluxes are calculated using a modified version of the
    !   *Strobel (1978) radiation scheme (*Shine and *Rickaby, 1990).
    !   Fluxes for ozone at the model half-levels are calculated using
    !   the ozone distribution. HIGHER SHORTWAVE RESOLUTION BASED ON WMO (1986)
    !    
    !   REFERENCES
    !   ----------
    !   *Rodgers (1967)
    !   *Shine and *Rickaby (1990) Atmospheric Ozone pp597-600.
    !   *Strobel (1978), J. Geophys. Res.
    !   *WMO (1986)
    !   *MATTHES ET AL. (2004), J. GEOPHYS. RES., 109, D06101, doi:10.1029/2003JD004012
    ! ---------------------------------------------------------------
    !
    ! 0.1 ARGUMENTS
    ! -------------
    ! INPUT:
    ! ------
    INTEGER, INTENT(in) :: kbdim    !first dimension of 2-d arrays
    INTEGER, INTENT(in) :: kproma   !number of longitudes
    INTEGER, INTENT(in) :: kswlevp1 !number of short wave half levels
    INTEGER, INTENT(in) :: imode    !1 for downward solar fluxes
                                    !2 for upward solar fluxes
    !
    ! paltsw - Tropospheric albedo for uv/vis
    ! po3c   - O3 column in path length (atm-cm) (1 atm-cm = 1000 DU)
    ! po2c   - O2 column in path length (atm-cm) (1 atm-cm = 1000 DU)
    !
    REAL(dp), DIMENSION(kbdim),          INTENT(in) :: paltsw
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(in) :: po3c
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(in) :: po2c  ! fb_mk_20150113
    !
    ! OUTPUT:
    ! -------
    ! zhart  - Total Flux in Hartley bands      (only downward)
    ! zflux  - Total,flux
    ! zherz  - Total Flux in Herzberg continuum (only downward)
    ! zhug   - Total Flux in Huggins bands
    ! zchap  - Total Flux in Chappuis bands
    ! pfluxc - Total upward flux clear sky 
    !
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(inout) :: zhart
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(out)   :: zflux
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(inout) :: zherz
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(out)   :: zhug
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(out)   :: zchap
    REAL(dp), DIMENSION(kbdim,kswlevp1), OPTIONAL, TARGET, INTENT(out) :: pfluxc

    
    ! 0.3 LOCAL ARRAYS
    ! ----------------
    !
    ! zxl3t  - Ozone column
    ! zxl3   - total O3 column density in cm-2
    ! zxl2   - total O2 column density in cm-2
    ! zsec   - magnification factor 
    ! zfluxc - dummy array clear sky upward flux (experimental)
    ! zflupc - pointer to clear sky upward flux
    !
    REAL(dp), DIMENSION(kbdim)          :: zxl3t
    REAL(dp), DIMENSION(kbdim)          :: zxl3
    REAL(dp), DIMENSION(kbdim)          :: zxl2
    REAL(dp), DIMENSION(kbdim,kswlevp1), TARGET  :: zflupc
    REAL(dp), DIMENSION(:,:),            POINTER :: zfluxc
    !
    REAL(dp), PARAMETER :: amagd = 1.900_dp
    INTEGER             :: jk, i        ! loop counter
    !
    !  2. INITIALISATION OF SOME DATA                                              
    !  ------------------------------ 
    IF (PRESENT(pfluxc)) THEN
       zfluxc => pfluxc
    ELSE
       zfluxc => zflupc
    END IF
    !
    !*   2.1 Initialise fluxes to zero.
    !
    zflux(:,:) = 0.0_dp
    zhug(:,:)  = 0.0_dp
    zchap(:,:) = 0.0_dp
    zfluxc(:,:)= 0.0_dp
    
    !*  2.2 Calculate Flux at top of atmosphere from max and min values 
    !       and state of solar cycle:
    !   Note:
    !   Initialization moved to subroutine fubrad_global_flux_ini.
    !
    
    !  3. Calculate downward fluxes if *imode* = 1.
    !  -------------------------------------------
    !zsec(:) = 35./(SQRT(1224.*(mu0(:)**2)+1.)) 

    SELECT CASE(imode)
    CASE(1)
       zhart(:,:) = 0.0_dp
       zherz(:,:) = 0.0_dp
       DO jk = 1, kswlevp1
          
          zxl3(1:kproma)  = po3c(1:kproma,jk) * przsec(1:kproma)
          !
          !   Instead of a constant O2 concentration, set by the constant oxndc,
          !   it is now possible to have a 3-d varying O2 filed.
          !
          zxl2(1:kproma)  = po2c(1:kproma,jk) * przsec(1:kproma)
          !
          ! Calculate Fluxes in all spectral intervals, and total flux (zflux) 
          ! for Chappuis, Huggins, Hartley and Herzberg
          !                                                   
          !CDIR unroll=nhug
          DO i = 1, nhug
             zhug(1:kproma,jk) = zhug(1:kproma,jk) + &
                                 flx3(zFhug(i), zsighug(i), zxl3(1:kproma))
          END DO
          !CDIR unroll=nhart
          DO i = 1, nhart
             zhart(1:kproma,jk) = zhart(1:kproma,jk) + &
                                  flx3(zFhart(i), zsighart(i), zxl3(1:kproma))
          END DO
          !CDIR unroll=nherz
          DO i = 1, nherz
             zherz(1:kproma,jk) = zherz(1:kproma,jk) + &
                                  flx32(zFherz(i), zsigherz2(i), zsigherz3(i), &
                                                zxl2(1:kproma), zxl3(1:kproma))
          END DO
          !CDIR unroll=nchap
          DO i = 1, nchap
             zchap(1:kproma,jk) = zchap(1:kproma,jk) + &
                                  flx3(zFchap(i), zsigchap(i), zxl3(1:kproma))
          END DO
          zflux(1:kproma,jk) = fladd + zchap(1:kproma,jk) + &
                                        zhug(1:kproma,jk) + &
                                       zherz(1:kproma,jk)

       END DO !LOOP OVER LEVELS 
       !
       !  4. Calculate upward fluxes if *imode* = 2.
       !  --------- -------- ------ -- ------- - --
       !  Note that this isn't strictly accurate: some summation over
       !  all levels incorporating the fractional cloudiness ought
       !  to be used, rather than considering the total atmosphere
       !  in the reflected beam.
       !
    CASE(2) 
       !
       DO jk = 1, kswlevp1
          !
          zxl3t(1:kproma) = ( przsec(1:kproma) * po3c(1:kproma,kswlevp1) &
                            + amagd *(po3c(1:kproma,kswlevp1) -        &
                                      po3c(1:kproma,jk)) )
         
          !CDIR unroll=nhug
          DO i = 1, nhug
             zhug(1:kproma,jk)= zhug(1:kproma,jk) + &
                                flx3(zFhug(i), zsighug(i), zxl3t(:kproma))
          END DO
          !CDIR unroll=nchap
          DO i = 1, nchap
             zchap(1:kproma,jk) = zchap(1:kproma,jk) + &
                                  flx3(zFchap(i), zsigchap(i), zxl3t(1:kproma))
          END DO
          !
          ! To calculate the upward fluxes Huggins and Chappuis bands are 
          ! considered, also the flux in the non-absorbing bands from
          ! 362.5 - 407.5 nm (fladd).
          !
          zflux(1:kproma,jk)  = -paltsw(1:kproma) * (zhug(1:kproma,jk) + &
                                                    zchap(1:kproma,jk) + fladd)
          zfluxc(1:kproma,jk) = -altswc(1:kproma) * (zhug(1:kproma,jk) + &
                                                    zchap(1:kproma,jk) + fladd)

          zchap(1:kproma,jk) = zchap(1:kproma,jk) * (-paltsw(1:kproma))
          zhug (1:kproma,jk) = zhug (1:kproma,jk) * (-paltsw(1:kproma))
          !
       END DO
    END SELECT
    !
    RETURN
  END SUBROUTINE o3fluxes
  !
  ! ======================================================================
  ! 
  ELEMENTAL FUNCTION flx3(zFlux, zsig, zl3) RESULT(flx)
    !
    REAL(dp), INTENT(in) :: zFlux ! solar flux in W m^-2
    REAL(dp), INTENT(in) :: zsig  ! O3 absorbtion cross section in cm^2
    REAL(dp), INTENT(in) :: zl3   ! O3 total column density in cm^-2
    !
    REAL(dp) :: flx  ! Solar flux at the model half level
    ! -
    flx = zFlux * EXP(-zsig * zl3)
    RETURN
  END FUNCTION flx3
  ! ======================================================================
  ELEMENTAL FUNCTION flx32(zFlux, zsig2, zsig3, zl2, zl3) RESULT(flx)
    !
    REAL(dp), INTENT(IN) :: zFlux ! solar flux in W m^-2
    REAL(dp), INTENT(IN) :: zsig2 ! O2 absorbtion cross section in cm^2
    REAL(dp), INTENT(IN) :: zsig3 ! O3 absorbtion cross section in cm^2
    REAL(dp), INTENT(IN) :: zl2   ! O2 total column density in cm^-2
    REAL(dp), INTENT(IN) :: zl3   ! O3 total column density in cm^-2
    !
    REAL(dp) :: flx  ! Solar flux at the model half level
    !
    flx = zFlux * EXP(-zsig2*zl2 -zsig3*zl3)
    !
  END FUNCTION flx32
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  SUBROUTINE ox_pathlength(pox,poxc,sort,kproma,kbdim,kklev,kswlev,pph)
    !
    ! AUTHOR
    ! ------
    ! K. NISSEN, FU-BERLIN,  5.8.2005
    ! M. Kunze,  FU-BERLIN, 13.1.2015 - serves for both O3 and O2
    !
    ! PURPOSE
    ! -------
    ! CONVERT OZONE/OXYGEN FROM MASS MIXING RATIO ON FULL MODEL LEVELS 
    ! TO COLUMN DENSITY IN CM-2 ON HALF LEVELS
    ! THIS IS NEEDED AS INPUT FOR THE O3FLUXES SUBROUTINE
    !
    ! FOR DETAILS ON THIS UNIT SEE HANDBOOK FOR MAP VOLUME 16 P.222
    !
    ! INPUT:
    ! ------
    INTEGER,INTENT(in) :: kproma  ! number of longitudes
    INTEGER,INTENT(in) :: kbdim   ! first dimension of 2D arrays
    INTEGER,INTENT(in) :: kklev   ! number of full levels
    INTEGER,INTENT(in) :: kswlev  ! number of full shortwave levels
    !
    ! pph  - half level pressure in Pa
    ! pox  - Ozone/Oxygen mass mixing ratio on full levels
    ! sort - character switch to choose O3 or O2
    REAL(dp),         INTENT(in), DIMENSION(kbdim,kklev+1):: pph
    REAL(dp),         INTENT(in), DIMENSION(kbdim,kklev)  :: pox
    CHARACTER(len=2), INTENT(in)                          :: sort
    !
    ! OUTPUT:
    ! -------
    ! poxc  - Ozone column density in cm-2 on half levels
    REAL(dp),INTENT(out), DIMENSION(kbdim,kswlev+1) :: poxc
    !
    ! LOCAL
    ! -----
    INTEGER :: jk
    !
    ! parameters mo3, mo2 are set. mo3 is set to be compatible with original 
    ! constant const=46.6968 used for ozone
    !
    REAL(dp), PARAMETER :: mO3 = 47.994250003460472_dp ! [g mol-1]
    REAL(dp), PARAMETER :: mO2 = 31.996166668973647_dp ! [g mol-1]
    REAL(dp)            :: mX
    REAL(dp)            :: const  ! [part./cm^2 * 1/Pa]
    !
    REAL(dp), DIMENSION(kbdim,kklev) :: zox ! column between two half levels
    REAL(dp), DIMENSION(kbdim,kklev) :: zdp ! Pressure thickness in Pa
    !
    SELECT CASE(sort)
    CASE('o3','O3')
       mX = mO3
    CASE('o2','O2')
       mX = mO2
                               ! the uppermost level from the model top to TOA
    CASE DEFAULT 
       PRINT *,TRIM(submodstr)//' subroutine ox_pathlength '//sort// &
               ' not supported'
       CALL exit(0)
    END SELECT
    !
    ! Loschmidt's number in cm^-3 (alos = 2.687E19) is not used
    ! here in the denominator, as also the factor alos is omitted
    ! in the subroutine o3fluxes.
    !
    const = N_a/ (mX * g) * 0.1_dp  ! [part./cm^2 * 1/Pa], 
    !                               ! 0.1 = 10^3 * 10^-4
    !                               !       10^3  to convert g mol-1 to kg mol-1
    !                               !             10^-4 to convert m^-2 to cm^-2
    ! --------------------------------------------------------------------
    !
    ! * 1. Pressure thickness in Pa
    zdp(1:kproma,:) = pph(1:kproma,2:kklev+1)-pph(1:kproma,1:kklev)
    !
    ! *  2. Column in path length between two half levels:
    !       MMR * density_air*(N_Avogadro/MO3*Lochschmidt's number)= 
    !       MMR * (dp/dz) *dz * (1/g) * (N_Avogadro/(MO3*Lochschmidt's number)
    !       N_Avogadro/(MO3*Lochschmidt's number) gives 46.6968.
    !
    !       In subroutine o3fluxes poxc in cm-2 (Number of particles per cm-2)
    !       is multiplied with the absorbtion cross section (in cm^2 per particle) 
    !       to get the dimensionless optical depth.
    !
    zox(1:kproma,:) = pox(1:kproma,:) * zdp(1:kproma,:) * const
    ! 
    !*   3. Ozone/Oxygen column between top and half level in path length
    poxc(1:kproma,1)= 0._dp
    DO jk = 2,  kswlev+1
       poxc(1:kproma,jk)= poxc(1:kproma,jk-1) + zox(1:kproma,jk-1)
    ENDDO
    ! --------------------------------------------------------------------
    RETURN
  END SUBROUTINE ox_pathlength
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE lymanalpha(nbdim, nproma, nlev, flux_lya )
    !---------------------------------------------------------------------------
    ! PURPOSE:  
    ! Calculate the effective absorption cross section of O2 in the Lyman-Alpha
    ! bands and an effective O2 optical depth at all altitudes.
    ! Derive flux at all altitudes in Wm-2
    ! Parameterized after:  Chabrillat, S., and G. Kockarts,
    ! Simple parameterization of the   
    ! absorption of the solar Lyman-Alpha line, Geophysical Research Letters,
    ! Vol.24, No.21, pp 2659-2662, 1997.
    !---------------------------------------------------------------------------
    ! PARAMETERS:
    ! nbdim    - INTEGER, number of horizontal elements (max.)            (I)
    ! nproma   - INTEGER, number of horizontal elements (used)            (I)
    ! nlev     - INTEGER, number of FUBRAD half levels                    (I)
    ! flux_lya - REAL, flux at lyman-alpha on half levels                 (O)
    !---------------------------------------------------------------------------
    ! EDIT HISTORY:
    ! 01/15/2002 Taken from Sasha Madronich's TUV Version 4.1a, Doug Kinnison
    ! 01/15/2002 Upgraded to F90, DK
    ! 01/24/2006 changes for ECHAM5/MESSy K.Nissen, FU Berlin
    !---------------------------------------------------------------------------
    
    !---------------------------------------------------------------------------
    !       ... Dummy arguments
    !---------------------------------------------------------------------------
    ! Input:
    ! ------
    INTEGER,                 INTENT(in) :: nbdim, nproma, nlev  
    !
    ! Output:
    ! -------
    REAL(dp), DIMENSION(:,:),INTENT(out) :: flux_lya ! flux at lyman alpha
    
    !---------------------------------------------------------------------------
    !     ... Local variables
    !---------------------------------------------------------------------------
    REAL(dp), PARAMETER, DIMENSION(3) :: b = &
         (/ 6.8431E-01_dp,  2.29841E-01_dp,  8.65412E-02_dp /)
    REAL(dp), PARAMETER, DIMENSION(3) :: c = &
         (/ 8.22114E-21_dp, 1.77556E-20_dp,  8.22112E-21_dp /)
    REAL(dp), PARAMETER, DIMENSION(3) :: d = &
         (/ 6.0073E-21_dp,  4.28569E-21_dp,  1.28059E-20_dp /)
    REAL(dp), PARAMETER, DIMENSION(3) :: e = &
         (/ 8.21666E-21_dp, 1.63296E-20_dp,  4.85121E-17_dp /)
    
    INTEGER  :: i, k
    REAL(dp), DIMENSION(nbdim) :: term
    REAL(dp), DIMENSION(nbdim) :: ro2, rm
    !
    ! o2col   - slant overhead O2 column (molec. cm-2) at half levels
    ! xsef_o2 - molecular absorption cross section in LA bands     
    ! ---------------------------------------------------------------
    REAL(dp), DIMENSION(nbdim) :: xsef_o2, o2col
    
    flux_lya(:,:) = 0._dp
    DO k = 1, nlev
       !
       ! calculate O2(slant) column on half levels ECHAM method
       o2col(1:nproma) = po2c(1:nproma,k) * przsec(1:nproma)
       !
       ! Calculate reduction factors at every altitude
       ro2(1:nproma) = 0._dp
       rm(1:nproma)  = 0._dp
       DO i = 1,3 
          term(1:nproma) = e(i) * o2col(1:nproma)
          WHERE ( term(1:nproma) < 100._dp ) 
             ro2(1:nproma) = ro2(1:nproma) + d(i) * EXP( -term(1:nproma) )
          END WHERE
          term(1:nproma) = c(i) * o2col(1:nproma)
          WHERE ( term(1:nproma) < 100._dp ) 
             rm(1:nproma) = rm(1:nproma) + b(i) * EXP( -term(1:nproma) )
          END WHERE
       END DO
       
       ! effective aborption cross section O2:
       WHERE (rm(1:nproma) > 0._dp)
          xsef_o2(1:nproma) = ro2(1:nproma)/rm(1:nproma)
       ELSE WHERE
          xsef_o2(1:nproma) = 1.2E-20_dp
       END WHERE
       !
       flux_lya(1:nproma,k) = zFlya * EXP(-xsef_o2(1:nproma)*o2col(1:nproma)) 
       !
    END DO
    
  END SUBROUTINE lymanalpha
  !
  !---------------------------------------------------------------------------
  ELEMENTAL FUNCTION HARTLEY_EFF(pres) RESULT(hart_eff)
    !
    ! PURPOSE:
    ! --------
    ! CALCULATE EFFICIENCY FACTOR FOR THE HARTLEY BAND
    ! THIS FOLLOWS Mlynczak and Solomon, JGR, vol 98, p 10517, 1993   
    !
    ! AUTHOR: 
    ! -------
    ! K. Nissen, 16.1.2006, Free University of Berlin 
    !
    REAL(dp), INTENT(IN)  :: PRES !pressure in hPa
    REAL(dp)              :: hart_eff
    !
    REAL(dp) :: X
    REAL(dp), PARAMETER, DIMENSION(4) :: ca = &
                      (/0.6696950_dp, -0.009682_dp,  0.033093_dp, 0.017938_dp /)
    REAL(dp), PARAMETER, DIMENSION(4) :: cb = &
                      (/0.926210_dp,   0.133960_dp, -0.076863_dp, 0.006897_dp /)
    !
    IF (pres <= 1.e-4_dp) THEN
       hart_eff = .7_dp
    ELSE IF (pres <= 1.e-2_dp) THEN
       x = LOG10(pres)+3._dp
       hart_eff = ca(1) + x*ca(2) + x*x*ca(3) + x*x*x*ca(4)
    ELSE IF (pres <= 1._dp) THEN
       x = LOG10(pres)+1._dp
       hart_eff = cb(1) + x*cb(2) + x*x*cb(3) + x*x*x*cb(4)
    ELSE
       hart_eff = 1._dp
    END IF
    !
    RETURN
  END FUNCTION HARTLEY_EFF
  !
! ***************************************************************************
END MODULE  messy_rad_fubrad
! ***************************************************************************
