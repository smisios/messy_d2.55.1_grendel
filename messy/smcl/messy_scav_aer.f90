MODULE messy_scav_aer

! This module contains all the subroutines and functions 
! necessary for aerosol scavenging
! Author: Holger Tost, Mainz, July 2004


  USE messy_main_constants_mem,    ONLY: dp, i4
  USE messy_scav_mem,              ONLY: aer_count, c_all

  IMPLICIT NONE

  SAVE

  PUBLIC :: aer_scav, aer_scav2
  PUBLIC :: aer_scav_cloud
  PUBLIC :: lambda, airvisc, rey_rain
  PUBLIC :: t_aero_inp
  PUBLIC :: t_aer_models
  PUBLIC :: mass_2d_ptr_array, num_2d_ptr_array
  PUBLIC :: aermodel

  ! switch for aerosol scavenging coefficient calculation
  INTEGER(i4), PUBLIC :: coeff_para     
  INTEGER,     PUBLIC :: max_mode

! FIELDS FOR AEROSOL DEPOSITION
  TYPE t_aero_inp
    LOGICAL                               :: lcalc_nucl_aer =.FALSE.
    INTEGER                               :: lmode          =  0
    REAL(dp), DIMENSION(:),       POINTER :: aersigma       => NULL()
    REAL(dp), DIMENSION(:),       POINTER :: cmr2mmr        => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: wetradius      => NULL()
    REAL(dp), DIMENSION(:,:,:),   POINTER :: wetrad         => NULL()  
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: inputrad       => NULL()   
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: nfracnuc       => NULL()   
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: mfracnuc       => NULL()   
    REAL(dp), DIMENSION(:,:,:),   POINTER :: nfrac          => NULL()   
    REAL(dp), DIMENSION(:,:,:),   POINTER :: mfrac          => NULL()  
    REAL(dp), DIMENSION(:,:,:),   POINTER :: philfrac       => NULL()
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: inputpf        => NULL()   
  END TYPE t_aero_inp
!-
  TYPE mass_2d_ptr_array
    INTEGER                          :: number       =  0
    INTEGER, DIMENSION(:,:), POINTER :: mass_int_att => NULL()
    LOGICAL, DIMENSION(:,:), POINTER :: mass_log_att => NULL()
  END TYPE mass_2d_ptr_array
  TYPE num_2d_ptr_array
    INTEGER                          :: number       =  0
    INTEGER, DIMENSION(:,:), POINTER :: num_int_att  => NULL()
    LOGICAL, DIMENSION(:,:), POINTER :: num_log_att  => NULL()
  END TYPE num_2d_ptr_array

  INTEGER, PARAMETER, PUBLIC :: spec_idx      = 1
  INTEGER, PARAMETER, PUBLIC :: spec_mode     = 2 
  INTEGER, PARAMETER, PUBLIC :: spec_idx_comp = 3
  INTEGER, PARAMETER, PUBLIC :: spec_sol      = 1
  INTEGER, PARAMETER, PUBLIC :: spec_hetice   = 2
  INTEGER, PARAMETER, PUBLIC :: aer_int_max  = 3
  INTEGER, PARAMETER, PUBLIC :: aer_log_max  = 2
!-
  INTEGER, PUBLIC, PARAMETER :: NMAX_AEROMODELS = 10
  TYPE t_aer_models
    INTEGER            :: aeromodnum = 0
    CHARACTER(LEN=10), DIMENSION(NMAX_AEROMODELS) :: aermodname
    INTEGER,           DIMENSION(NMAX_AEROMODELS) :: aermodmethod   
    TYPE(t_aero_inp),        POINTER, DIMENSION(:)           :: &
                                       aer_input     => NULL()
    TYPE(mass_2d_ptr_array), POINTER, DIMENSION(:)           :: &
                                       aer_mass_attr => NULL()
    TYPE(num_2d_ptr_array),  POINTER, DIMENSION(:)           :: &
                                       aer_num_attr  => NULL()
  END TYPE t_aer_models
!-
  TYPE(t_aer_models), DIMENSION(:), POINTER :: aermodel     => NULL()

  REAL (dp), PUBLIC, POINTER :: rad_aer(:,:,:,:)    => NULL()   
  REAL (dp), PUBLIC, POINTER :: sigma(:)            => NULL() 

  PRIVATE
  ! a small number for fluxes, lwc, droplet size for stability reasons
  PRIVATE :: dp
  REAL(dp), PARAMETER, PRIVATE :: small = 1.e-30_dp  
 
  INTRINSIC EXP, LOG, SQRT, ATAN, MIN, NULL 

CONTAINS
!===============================================================================

  SUBROUTINE aer_scav(nspec, lproma, js, mode, model, itype, p_xt, wetflx, &
                      ztmst, rainflux, cover, r_ra, r_aer, &
                      temp, press, rhoa, snowflux, vol, phase)

!     routine in which aerosol impaction scavenging will be calculated
!     works for both convective and large-scale precipitation following 
!     the old routine of P. Stier for M7 aerosol module
!     added a parametrisation for collection efficiency using a 
!     formula from Slinn for aerosol radius
    
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nspec          ! number of packed aerosol tracers
     integer, INTENT(IN) :: lproma         ! number of packed boxes
     integer, INTENT(IN) :: js             ! number of tracer set
     integer, INTENT(IN) :: mode           ! number of aerosol modes 
     integer, INTENT(IN) :: model          ! number of running aerosol models
     integer, INTENT(IN) :: itype          ! 1 = mass, 2 = numbers
     integer, INTENT(IN) :: phase          ! 1 = liquid, 2 = ice
     
     REAL(dp), INTENT(IN) :: r_aer(mode, lproma) 
     REAL(dp), INTENT(IN) :: r_ra(lproma)      ! average droplet radius (in mm)
     REAL(dp), INTENT(IN) :: ztmst             ! timestep_length
     ! volume of box to change from concentration (molecules/cm^3 to molecules)
     REAL(dp), INTENT(IN) :: vol(lproma)  
     ! precipitation flux(in kg/(m^2*s))
     REAL(dp), INTENT(IN) :: rainflux(lproma) 
     ! snow flux (in kg/(m^2*s))
     REAL(dp), INTENT(IN) :: snowflux(lproma)
     ! cloud cover
     REAL(dp), INTENT(IN) :: cover(lproma)
     ! temperature in that box
     REAL(dp), INTENT(IN) :: temp(lproma)
     ! pressure in that box
     REAL(dp), INTENT(IN) :: press(lproma)
     ! density in that box (in kg/m^3)
     REAL(dp), INTENT(IN) :: rhoa(lproma)

     ! species field
     REAL(dp), INTENT(INOUT) :: p_xt(nspec,lproma) 

     ! aerosol deposition field (in molecules)
     REAL(dp), INTENT(INOUT) :: wetflx(aer_count(js)%numbers(c_all), lproma)  
     
! local

     REAL(dp) :: scav_coeff(mode,lproma)   ! scavenging coefficient
     ! parameter for the coefficient for snow/ice precipitation
     REAL(dp) :: velo_p             ! terminal velocity of particles
     REAL(dp) :: schm_p             ! schmidt number of particle
     REAL(dp) :: sto_p              ! stokes number of particle
     REAL(dp) :: phi                ! r_aer / r_rain

     REAL(dp), PARAMETER :: visc_w = 1.65e-3_dp  ! in Pa s

     REAL(dp) :: par1, par2, par3, S_star      ! auxiliary parameters
    ! new value in the cloud covered fraction of that box
     REAL(dp) :: new_val         
     INTEGER  :: l, jt, jl, idx1, idx2
     INTEGER  :: tmode(nspec), idx_flx(nspec)
     REAL(dp) :: zeicaes_help(lproma)
     REAL(dp) :: r_rain_help(lproma)  ! droplet radius in m
     REAL(dp) :: visc_help(lproma)
     REAL(dp) :: lam_help(lproma)
     REAL(dp) :: rey_r_help(lproma)   ! Reynolds number of rain
     REAL(dp) :: velo_r_help(lproma)  ! terminal velocity of rain
     REAL(dp) :: omega_help(lproma)   ! viscosity water / viscosity air
     REAL(dp) :: S_star_help(lproma)
     LOGICAL  :: lswitch(lproma)

     SELECT CASE(itype)
     CASE(1)
       idx_flx(:) = &
         aermodel(js)%aer_mass_attr(model)%mass_int_att(:,spec_idx_comp)
       tmode(:)   = aermodel(js)%aer_mass_attr(model)%mass_int_att(:,spec_mode)
     CASE(2)
       idx_flx(:) = &
         aermodel(js)%aer_num_attr(model)%num_int_att(:,spec_idx_comp)
       tmode(:)   = aermodel(js)%aer_num_attr(model)%num_int_att(:,spec_mode)
     END SELECT

     scav_coeff(:,:) = 0.0_dp
!   scavenging parameter by liquid rain

     IF (phase == 1) THEN

       SELECT CASE(coeff_para)

       CASE(1)         ! constant scav_param following Feichter
         DO jl=1,lproma
           zeicaes_help(jl) = 0.0_dp
           IF (cover(jl) .GE. 1.e-7_dp) THEN
             zeicaes_help(jl) = 0.1_dp * rainflux(jl)/cover(jl) * 3600._dp
           ENDIF
         ENDDO
!    adding up snow/ice and liquid precipitation scavenging values      
!CDIR NOLOOPCHG
         DO l=1, mode
           DO jl=1,lproma
             scav_coeff(l,jl) =  zeicaes_help(jl)
           ENDDO
         ENDDO

       CASE(2)
       ! radius dependent parameters following Seinfeld & Pandis,98 eq. 20-60
       ! calculation of collection efficiency eq. 20-56 (original Slinn 83)

         lswitch(:)     = .FALSE.
         r_rain_help(:) = 0._dp
         visc_help(:)   = 0._dp
         lam_help(:)    = 0._dp
         rey_r_help(:)  = 0._dp
         velo_r_help(:) = 0._dp
         omega_help(:)  = 0._dp
         S_star_help(:) = 0._dp
         DO jl=1,lproma
           r_rain_help(jl) = r_ra(jl) * 1.e-3_dp        ! conversion of rain radius to m
           IF ( (cover(jl) .GE. 1.e-7_dp) .AND. (r_rain_help(jl).GT.1.e-7_dp) ) THEN
             lswitch(jl) = .TRUE.
             ! viscosity of air
             visc_help(jl) = airvisc(temp(jl))
             ! mean free path
             lam_help(jl) = lambda(temp(jl), press(jl))
             ! reynolds number of rain droplets
             rey_r_help(jl)  = rey_rain(r_rain_help(jl), rhoa(jl), visc_help(jl))
             ! terminal velocity of droplets
             velo_r_help(jl) = rvelo(r_rain_help(jl), visc_help(jl), rey_r_help(jl), rhoa(jl))

             omega_help(jl)  = visc_help(jl) / visc_w

             S_star = 1._dp/12._dp * log(1._dp + rey_r_help(jl)) + 1.2_dp
             S_star_help(jl) = S_star / (1._dp + log(1._dp + rey_r_help(jl)))
           ENDIF
         ENDDO

         DO l=1, mode
           DO jl=1,lproma
             IF( lswitch(jl) .AND. (r_aer(l,jl).gt.0._dp) ) THEN
               schm_p = schmidt(visc_help(jl), rhoa(jl), &
                                Diffcoeff(temp(jl),r_aer(l,jl),lam_help(jl),visc_help(jl)))
               velo_p = tvelo(r_aer(l,jl), visc_help(jl), lam_help(jl))
               sto_p  = stokes(velo_r_help(jl), velo_p, r_aer(l,jl), visc_help(jl), r_rain_help(jl), lam_help(jl))
               phi    = r_aer(l,jl) / r_rain_help(jl)

               par1 = 1._dp + 0.4_dp * sqrt(rey_r_help(jl)) * schm_p**(1._dp/3._dp)
               par1 = par1 + 0.16_dp * sqrt(rey_r_help(jl)*schm_p)
               par1 = par1 * 4._dp / (rey_r_help(jl) * schm_p)

               par2 = omega_help(jl)
               par2 = par2 + phi*(1._dp + 2._dp* sqrt(rey_r_help(jl)))
               par2 = par2 * 4._dp * phi

               par3 = (sto_p - S_star_help(jl))/(sto_p - S_star_help(jl) + 2._dp/3._dp)
               par3 = MAX(par3,0._dp)
               par3 = par3**(3._dp/2._dp)

         ! following Seinfeld & Pandis,98 eq. 20-60
!    adding up snow/ice and liquid precipitation scavenging values

               scav_coeff(l,jl) = (par1 + par2 + par3)/r_ra(jl) * 0.75_dp * &
                                rainflux(jl)/cover(jl) * ztmst
             ENDIF
           ENDDO
         ENDDO

       END SELECT

     ENDIF

!    scavenging by snow/ice precipitation
     IF (phase == 2) THEN

       DO jl=1,lproma
         zeicaes_help(jl) = 0.0_dp
         IF (cover(jl) .GE. 1.e-7_dp) THEN
           zeicaes_help(jl) = 0.1_dp * snowflux(jl)/cover(jl) * 3600._dp
         ENDIF
       ENDDO
!    adding up snow/ice and liquid precipitation scavenging values      
!CDIR NOLOOPCHG
       DO l=1, mode
         DO jl=1,lproma
           scav_coeff(l,jl) = scav_coeff(l,jl) + zeicaes_help(jl)
         ENDDO
       ENDDO
     ENDIF

     
     do jt=1,nspec     
       idx1 = idx_flx(jt)
       idx2 = tmode(jt)
       do jl=1,lproma
         if (r_aer(idx2,jl) >  0._dp) then 
           new_val = p_xt(jt,jl) * exp(-scav_coeff(idx2,jl))
           wetflx(idx1,jl) = wetflx(idx1,jl) + cover(jl) * vol(jl) * 1e6 *       &
                          p_xt(jt,jl) * (1._dp - exp(-scav_coeff(idx2,jl)))

           p_xt(jt,jl) = (1._dp-cover(jl))*p_xt(jt,jl) + cover(jl) * new_val
         endif
       enddo
     enddo

   END SUBROUTINE aer_scav

!-------------------------------------------------------------------------
  SUBROUTINE aer_scav2(nspec, lproma, js, mode, model, itype, p_xt, wetflx, &
                      ztmst, rainflux, cover, r_ra, r_aer, &
                      temp, press, rhoa, snowflux, vol, phase)

!     routine in which aerosol impaction scavenging will be calculated
!     works for both convective and large-scale precipitation following 
!     the old routine of P. Stier for M7 aerosol module
!     added a parametrisation for collection efficiency using a 
!     formula from Slinn for aerosol radius
    
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nspec          ! number of packed aerosol tracers
     integer, INTENT(IN) :: lproma         ! number of packed boxes
     integer, INTENT(IN) :: js             ! number of tracer set
     integer, INTENT(IN) :: mode           ! number of aerosol modes 
     integer, INTENT(IN) :: model          ! number of running aerosol models
     integer, INTENT(IN) :: itype          ! 1 = mass, 2 = numbers
     integer, INTENT(IN) :: phase          ! 1 = liquid, 2 = ice
     
     REAL(dp), INTENT(IN) :: r_aer(mode, lproma) 
     REAL(dp), INTENT(IN) :: r_ra(lproma)      ! average droplet radius (in mm)
     REAL(dp), INTENT(IN) :: ztmst             ! timestep_length
     ! volume of box to change from concentration (molecules/cm^3 to molecules)
     REAL(dp), INTENT(IN) :: vol(lproma)  
     ! precipitation flux(in kg/(m^2*s))
     REAL(dp), INTENT(IN) :: rainflux(lproma) 
     ! snow flux (in kg/(m^2*s))
     REAL(dp), INTENT(IN) :: snowflux(lproma)
     ! cloud cover
     REAL(dp), INTENT(IN) :: cover(lproma)
     ! temperature in that box
     REAL(dp), INTENT(IN) :: temp(lproma)
     ! pressure in that box
     REAL(dp), INTENT(IN) :: press(lproma)
     ! density in that box (in kg/m^3)
     REAL(dp), INTENT(IN) :: rhoa(lproma)

     ! species field
     REAL(dp), INTENT(INOUT) :: p_xt(nspec,lproma) 

     ! aerosol deposition field (in molecules)
     REAL(dp), INTENT(INOUT) :: wetflx(aer_count(js)%numbers(c_all), lproma)  
     
! local

     REAL(dp) :: scav_coeff(mode,lproma)   ! scavenging coefficient
     ! parameter for the coefficient for snow/ice precipitation
     REAL(dp) :: velo_p             ! terminal velocity of particles
     REAL(dp) :: schm_p             ! schmidt number of particle
     REAL(dp) :: sto_p              ! stokes number of particle
     REAL(dp) :: phi                ! r_aer / r_rain

     REAL(dp), PARAMETER :: visc_w = 1.65e-3_dp  ! in Pa s
     INTEGER, PARAMETER  :: is = 6     ! iteration step for kmt calculation

     REAL(dp) :: par1, par2, par3      ! auxiliary parameters
    ! new value in the cloud covered fraction of that box
     REAL(dp) :: new_val         
     INTEGER  :: l, jt, jl, idx1, idx2, inum
     INTEGER  :: tmode(nspec), idx_flx(nspec)
     REAL(dp) :: zeicaes_help(lproma)
     REAL(dp) :: r_rain_help(lproma)    ! droplet radius in m
     REAL(dp) :: visc_help(lproma)
     REAL(dp) :: lam_help(lproma)
     REAL(dp) :: rey_r_help(is,lproma)  ! Reynolds number of rain
     REAL(dp) :: velo_r_help(is,lproma) ! terminal velocity of rain
     REAL(dp) :: omega_help(lproma)     ! viscosity water / viscosity air
     REAL(dp) :: S_star(is,lproma)
     LOGICAL  :: lswitch(lproma)
     REAL(dp) :: radfield(is+1), number(is,lproma)
     REAL(dp) :: del_num, fl1(lproma), fl2(lproma)
     REAL(dp) :: sc_coeff(is,mode,lproma)

     SELECT CASE(itype)
     CASE(1)
       idx_flx(:) = &
         aermodel(js)%aer_mass_attr(model)%mass_int_att(:,spec_idx_comp)
       tmode(:)   = aermodel(js)%aer_mass_attr(model)%mass_int_att(:,spec_mode)
     CASE(2)
       idx_flx(:) = &
         aermodel(js)%aer_num_attr(model)%num_int_att(:,spec_idx_comp)
       tmode(:)   = aermodel(js)%aer_num_attr(model)%num_int_att(:,spec_mode)
     END SELECT

     scav_coeff(:,:) = 0.0_dp
!   scavenging parameter by liquid rain

     IF (phase == 1) THEN

       SELECT CASE(coeff_para)

       CASE(1)         ! constant scav_param following Feichter
         DO jl=1,lproma
           zeicaes_help(jl) = 0.0_dp
           IF (cover(jl) > 1.e-7_dp) THEN
             zeicaes_help(jl) = 0.1_dp * rainflux(jl)/cover(jl) * 3600._dp
           ENDIF
         ENDDO
!    adding up snow/ice and liquid precipitation scavenging values      
!CDIR NOLOOPCHG
         DO l=1, mode
           DO jl=1,lproma
             scav_coeff(l,jl) =  zeicaes_help(jl)
           ENDDO
         ENDDO

       CASE(2)
       ! radius dependent parameters following Seinfeld & Pandis,98 eq. 20-60
       ! calculation of collection efficiency eq. 20-56 (original Slinn 83)

         lswitch(:)       = .FALSE.
         r_rain_help(:)   = 0._dp
         visc_help(:)     = 0._dp
         lam_help(:)      = 0._dp
         omega_help(:)    = 0._dp
         rey_r_help(:,:)  = 0._dp
         velo_r_help(:,:) = 0._dp
         S_star(:,:)      = 0._dp
         number(:,:)      = 0._dp

         radfield(:) = 0.
         radfield(1) = 5.00e-3
         radfield(2) = 2.00e-3
         radfield(3) = 1.00e-3
         radfield(4) = 0.50e-3
         radfield(5) = 0.20e-3
         radfield(6) = 0.10e-3
         DO jl=1,lproma
           IF (cover(jl) > 1.e-7_dp) THEN
             lswitch(jl) = .TRUE.
             ! viscosity of air
             visc_help(jl) = airvisc(temp(jl))
             omega_help(jl) = visc_help(jl) / visc_w
             ! mean free path
             lam_help(jl) = lambda(temp(jl), press(jl))

             fl1(jl) = 2.8e-5_dp * (rainflux(jl)*3600._dp) **0.324_dp
             fl2(jl) = - 98.5_dp * (rainflux(jl)*3600._dp) **(-0.522_dp)
           ENDIF
         ENDDO
         DO inum=1,is
           DO jl=1,lproma
             IF (lswitch(jl)) THEN
             ! reynolds number of rain droplets
               rey_r_help(inum,jl)  = rey_rain(radfield(inum), &
                                               rhoa(jl), visc_help(jl) )

               del_num = fl1(jl) * (radfield(inum) * 20._dp)**(-1.75_dp) * &
                         exp( fl2(jl) * (radfield(inum) * 20._dp)**2.25_dp)
               number(inum,jl) = del_num * (radfield(inum) - radfield(inum+1))
             ! terminal velocity of droplets
               velo_r_help(inum,jl) = rvelo(radfield(inum), visc_help(jl), &
                                            rey_r_help(inum,jl), rhoa(jl) )
      
               S_star(inum,jl) = 1._dp/12._dp * &
                                 log(1._dp + rey_r_help(inum,jl)) + 1.2_dp
               S_star(inum,jl) = S_star(inum,jl) / &
                                (1._dp + log(1._dp + rey_r_help(inum,jl)))
             ENDIF
           ENDDO

           DO l=1, mode
             DO jl=1,lproma
               IF( lswitch(jl) .AND. (r_aer(l,jl) > 0._dp) ) THEN
                 schm_p = schmidt(visc_help(jl), rhoa(jl),        &
                                  Diffcoeff(temp(jl),r_aer(l,jl), &
                                  lam_help(jl),visc_help(jl)))
                 velo_p = tvelo(r_aer(l,jl), visc_help(jl), lam_help(jl))
                 sto_p  = stokes(velo_r_help(inum,jl), velo_p, r_aer(l,jl), &
                                 visc_help(jl), radfield(inum), lam_help(jl))
                 phi    = r_aer(l,jl) / radfield(inum)

                 par1 = 1._dp + 0.4_dp * sqrt(rey_r_help(inum,jl)) * &
                        schm_p**(1._dp/3._dp)
                 par1 = par1 + 0.16_dp * sqrt(rey_r_help(inum,jl)*schm_p)
                 par1 = par1 * 4._dp / (rey_r_help(inum,jl) * schm_p)

                 par2 = omega_help(jl) + &
                        phi*(1._dp + 2._dp* sqrt(rey_r_help(inum,jl)))
                 par2 = par2 * 4._dp * phi

                 par3 = (sto_p - S_star(inum,jl)) / &
                        (sto_p - S_star(inum,jl) + 2._dp/3._dp)
                 par3 = MAX(par3,0._dp)
                 par3 = par3**(3._dp/2._dp)

         ! following Seinfeld & Pandis,98 eq. 20-60
!    adding up liquid precipitation scavenging values
                 sc_coeff(inum,l,jl) = (par1 + par2 + par3) / radfield(inum) * &
                                        0.75e-3_dp * rainflux(jl)/cover(jl)  * &
                                        ztmst
                 scav_coeff(l,jl) = scav_coeff(l,jl) + &
                                    sc_coeff(inum,l,jl) * number(inum,jl)
               ENDIF
             ENDDO
           ENDDO
         ENDDO
         DO l=1,mode
           DO jl=1,lproma
             IF (lswitch(jl)) &
               scav_coeff(l,jl) = scav_coeff(l,jl) / SUM(number(1:is,jl))
           ENDDO
         ENDDO
       END SELECT

     ENDIF

!    scavenging by snow/ice precipitation
     IF (phase == 2) THEN

       DO jl=1,lproma
         zeicaes_help(jl) = 0.0_dp
         IF (cover(jl) .GE. 1.e-7_dp) THEN
           zeicaes_help(jl) = 0.1_dp * snowflux(jl)/cover(jl) * 3600._dp
         ENDIF
       ENDDO
!    adding up snow/ice and liquid precipitation scavenging values      
!CDIR NOLOOPCHG
       DO l=1, mode
         DO jl=1,lproma
           scav_coeff(l,jl) = scav_coeff(l,jl) + zeicaes_help(jl)
         ENDDO
       ENDDO
     ENDIF

     
     do jt=1,nspec     
       idx1 = idx_flx(jt)
       idx2 = tmode(jt)
       do jl=1,lproma
         if (r_aer(idx2,jl) >  0._dp) then 
           new_val = p_xt(jt,jl) * exp(-scav_coeff(idx2,jl))
           wetflx(idx1,jl) = wetflx(idx1,jl) + cover(jl) * vol(jl) * 1e6 *    &
                             p_xt(jt,jl) * (1._dp - exp(-scav_coeff(idx2,jl)))

           p_xt(jt,jl) = (1._dp-cover(jl))*p_xt(jt,jl) + cover(jl) * new_val
         endif
       enddo
     enddo

   END SUBROUTINE aer_scav2

!-------------------------------------------------------------------------
 

!===============================================================================
!-------------------------------------------------------------------------
  SUBROUTINE aer_scav_cloud(nspec, lproma, js, mode, model, itype,     &
                            phase,   p_xt, newflx,  ztmst,             &
                            zlwc,    ziwc,      cover,                 &
                            r_cloud, r_aer,     temp,   press,  vol,   &
                            frac, philf )

!     routine in which aerosol in-cloud scavenging will be calculated
!     works for large-scale clouds and precipitation following 
!     the old routine of J. Feichter for M7 aerosol module 
!     added a new calculation of scavenging coefficient including radius 
!     dependecy of droplet and aerosol
!     following Pruppacher & Klett, Microphysics of clouds......
!     
!     Formula, P & K, 2000, page 723 eq 17-17:
!                                    scav_coeff = 1.35 * lwc * D * r^-2
!                                    LWC in g/cm^3, D is DiffusionCoefficient
!                                    r = droplet radius in cm

    Use messy_main_constants_mem,         only: pi
    Use messy_scav_mem,                   only: iscav_rate_hiTmp,&
                                                iscav_rate_loTmp_het,&
                                                iscav_rate_loTmp_hom

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nspec          ! number of packed aerosol tracers
    integer, INTENT(IN) :: lproma         ! number of packed boxes
    integer, INTENT(IN) :: js             ! number of tracer set
    integer, INTENT(IN) :: mode           ! number of aerosol modes 
    integer, INTENT(IN) :: model          ! number of running aerosol models
    integer, INTENT(IN) :: itype          ! 1 = mass, 2 = numbers
    integer, INTENT(IN) :: phase          ! 1 = liquid, 2 = ice

    ! average aerosol radius per mode (in m)
     REAL(dp), INTENT(IN) :: r_aer(mode,lproma)   

     ! hydrophilic fraction of 3rd moment per mode
     REAL(dp), INTENT(IN) :: philf(mode,lproma)   

     REAL(dp), INTENT(IN) :: r_cloud(lproma) ! average droplet radius (in mm)
     REAL(dp), INTENT(IN) :: ztmst          ! time step length
     REAL(dp), INTENT(IN) :: temp(lproma)   ! temperature in that box
     REAL(dp), INTENT(IN) :: press(lproma)  ! pressure in that box
     REAL(dp), INTENT(IN) :: zlwc(lproma)   ! cloud liquid water content (in kg/m3)
     REAL(dp), INTENT(IN) :: ziwc(lproma)   ! snow / ice water content (in kg/m3)
     REAL(dp), INTENT(IN) :: cover(lproma)  ! cloud cover
     ! volume to box to change from concentration (molecules/cm^3 to molecules)
     REAL(dp), INTENT(IN) :: vol(lproma)         

     REAL(dp), INTENT(INOUT) :: p_xt(nspec,lproma)    ! species field
     ! new nucleation part of wet deposition flux
     REAL(dp), INTENT(INOUT) :: newflx(aer_count(js)%numbers(c_all), lproma)

     REAL(dp), INTENT(IN) :: frac(mode,lproma) ! input value from cloud of activated fraction

     REAL(dp), PARAMETER :: limit =  0.99999_dp
     ! local variables
     REAL(dp) :: zeicaer
     REAL(dp) :: new_val, localcover(lproma)
     REAL(dp) :: scav_coeff(mode,lproma), scav_frac(mode,lproma) 
     REAL(dp) :: par_nucl(mode,lproma), par_nucl2(mode,lproma), par_nucl1(mode,lproma)
     ! mode field of species
     INTEGER  :: tmode(nspec)   
     ! index number of aerosol tracers in the aer_flx
     INTEGER  :: idx_flx(nspec)  
     ! LOGICAL of solubility of aerosol tracers in the aer_flx
     LOGICAL  :: nsol(nspec)
     ! LOGICAL of inhomogeneous freezing of aerosol tracers in the aer_flx
     LOGICAL  :: nhetice(nspec)

     INTEGER  :: idx1, idx2, l, jt, jl
     REAL(dp) :: zeicaer_help(lproma),zeicaes_help(lproma,nspec)
     REAL(dp) :: visc_help(lproma)
     REAL(dp) :: lam_help(lproma)

     SELECT CASE(itype)
     CASE(1)
       idx_flx(:) = &
         aermodel(js)%aer_mass_attr(model)%mass_int_att(:,spec_idx_comp)
       tmode(:)   = aermodel(js)%aer_mass_attr(model)%mass_int_att(:,spec_mode)
       nsol(:)    = aermodel(js)%aer_mass_attr(model)%mass_log_att(:,spec_sol)
       nhetice(:) = aermodel(js)%aer_mass_attr(model)%mass_log_att(:,spec_hetice)
     CASE(2)
       idx_flx(:) = &
         aermodel(js)%aer_num_attr(model)%num_int_att(:,spec_idx_comp)
       tmode(:)   = aermodel(js)%aer_num_attr(model)%num_int_att(:,spec_mode)
       nsol(:)    = aermodel(js)%aer_num_attr(model)%num_log_att(:,spec_sol)
       nhetice(:) = aermodel(js)%aer_num_attr(model)%num_log_att(:,spec_hetice)
     END SELECT
     
     do jl=1,lproma
       localcover(jl) = min(cover(jl), limit)
     enddo
     
     scav_coeff(:,:) = 0.0_dp
     par_nucl(:,:)   = 0.0_dp
     par_nucl2(:,:)  = 0.0_dp
     par_nucl1(:,:)  = 0.0_dp
     scav_frac(:,:)  = 0.0_dp
     zeicaes_help(:,:) = 0.0_dp

!    liquid nucleation scavenging
     liq_nuc_scav: IF (phase == 1) THEN

       SELECT CASE (coeff_para)

       CASE(1)

         DO jl=1,lproma
           IF(temp(jl) > 238.15_dp) THEN
             zeicaer_help(jl)=0.8_dp 
           ELSE
        !--- Pers. communication with J. Hendricks (DLR) change from 0.1:
        ! --- amount of liquid water at these temperatures should 
        !     anyhow be close to zero

              zeicaer_help(jl)=0.05_dp 
           ENDIF
         ENDDO

! applying liquid nucleation scavenging fraction to all modes
!CDIR NOLOOPCHG
         DO l=1, mode
           DO jl=1,lproma
             scav_frac(l,jl) = zeicaer_help(jl)  
           ENDDO
         ENDDO

       CASE(2)

!    new formulation of scavenging coefficient added
!    scavenging formulation due to brownian motion
         visc_help(:)=0._dp
         lam_help(:)=0._dp
         DO jl=1,lproma
           IF (r_cloud(jl) > small) THEN
             visc_help(jl) = airvisc(temp(jl))
             lam_help(jl)  = lambda(temp(jl), press(jl))
           ENDIF
         ENDDO
         DO l=1, mode
           DO jl=1,lproma
             IF ( (r_cloud(jl) > small) .AND. (r_aer(l,jl) > 0._dp) ) THEN
               ! conversion of lwc to g/cm^3
               zeicaer = 1.35_dp * (zlwc(jl)*1.e-3_dp) &         
                       * DiffCoeff(temp(jl),r_aer(l,jl),lam_help(jl),visc_help(jl)) &
                       ! conversion of r_cloud from mm to cm
                       / ((r_cloud(jl)*0.1_dp)*(r_cloud(jl)*0.1_dp))       
                  
               ! adding up liquid nucleation scavenging coefficient
               scav_coeff(l,jl) = zeicaer * ztmst  
               ! zeicaer is already a scavenging efficiency
               scav_frac(l,jl) = 1._dp - exp(-scav_coeff(l,jl))
             END IF
           END DO
         END DO
         

!         IF (aermodel(js)%aer_input(model)%lcalc_nucl_aer) THEN
           DO l=1,mode
             DO jl=1,lproma
               IF ( (r_cloud(jl) > small) .AND. (r_aer(l,jl) > 0._dp) ) THEN
!  nucleation is described via a parametrized estimated formula:
!  it should lead to the nucleation scavenging of particles above about 0.14 µm 
!  and is estimated from observations, leading to the scavenging of 
!  large particles
!     It has no references and it is just a guess:
!     nucleation:  atan((5e6*r_aer)**6) *2/pi
!     there is no time dependence included as it happens in one timestep during 
!     cloud formation, but of course it shows the dependence on the fraction of
!     the cloud water that rains out.....
!  philf representing the fraction of soluble material in case of mixed modes 
!      (as in MADE) should represent activation of the hydrophilic (or coated) 
!      compounds only
                 par_nucl1(l,jl)  = philf(l,jl) * &
                    MIN(atan((5.e6_dp*r_aer(l,jl))**6._dp) * 2._dp/pi , 1._dp - 1.e-8_dp)
               ENDIF
             ENDDO
           ENDDO
           IF (associated(aermodel(js)%aer_input(model)%nfrac)) THEN
              DO l=1,SIZE(aermodel(js)%aer_input(model)%nfrac,2)
                 DO jl=1,lproma
                    IF ( (r_cloud(jl) > small) .AND. (r_aer(l,jl) > 0._dp) ) &
                         par_nucl2(l,jl) = frac(l,jl)
                 END DO
              END DO
           END IF
     
!         print*, "parnucl", MAXVAL(par_nucl1(:,:)), MINVAL(par_nucl1(:,:)), &
!           aermodel(js)%aer_input(model)%lcalc_nucl_aer,  MAXVAL(par_nucl2(:,:)),&
!           MINVAL(par_nucl2(:,:))
         
         IF (aermodel(js)%aer_input(model)%lcalc_nucl_aer) THEN
            DO l=1,mode
               DO jl=1,lproma
                  par_nucl(l,jl) = par_nucl1(l,jl)
               ENDDO
            ENDDO
         ELSE
            DO l=1,SIZE(aermodel(js)%aer_input(model)%nfrac,2)
               DO jl=1,lproma
                  par_nucl(l,jl) = par_nucl2(l,jl)
               ENDDO
            ENDDO
         ENDIF
         
       END SELECT

   END IF liq_nuc_scav    

   ! snow/ice nucleation scavenging
   ice_nuc_scav: IF (phase == 2) THEN
     DO jl=1, lproma      
       IF(temp(jl) > 238.15_dp) THEN
         DO jt=1,nspec     
            zeicaes_help(jl,jt) = iscav_rate_hiTmp
         END DO
       ELSE

           !--- Pers. communication with J. Hendricks (DLR) change from 0.1:
           ! --- value chosen originally for Cirrus clouds, 
           !      but is also applied to Mixed phase clouds due to correspondence 
           !      with observations

         DO jt=1,nspec
           IF (nhetice(jt)) THEN
                zeicaes_help(jl,jt) = iscav_rate_loTmp_het
           ELSE
                  ! --- Pers. communication with J. Hendricks (DLR)
                  !     change from 0.1
                  ! --- value chosen originally for cirrus clouds, 
                  !     but is also applied to mixed phase clouds due to
                  !     correspondence with observations
                zeicaes_help(jl,jt) = iscav_rate_loTmp_hom
           ENDIF
           
         END DO
       END IF
     ENDDO
   END IF ice_nuc_scav

! applying scavenging fractions to aerosol species fields

   do jt=1,nspec     
      idx1 = idx_flx(jt)
      idx2 = tmode(jt)
      do jl=1,lproma
         if (r_aer(idx2,jl) >  0._dp) then
            ! apply scav_frac (= 0 if phase == 2),
            ! and zeicaes_help (= 0 if phase == 1)
            new_val = localcover(jl) * &
                 p_xt(jt,jl) * (1._dp - scav_frac(idx2,jl) - zeicaes_help(jl,jt))
            ! apply par_nucl only for compounds with nsol=.true.
            if (nsol(jt)) &
                 new_val = new_val * (1._dp - par_nucl(idx2,jl))
 !           if (.NOT. ((par_nucl(idx2,jl) < 1._dp) .AND. (par_nucl(idx2,jl) > 0._dp) )) &
 !             print*, "par_nucl ", par_nucl(idx2,jl),par_nucl1(idx2,jl), par_nucl2(idx2,jl)
            ! calculate freshly scavenged pseudo - flux
            newflx(idx1,jl) =  newflx(idx1,jl) + &
                              (localcover(jl) * p_xt(jt,jl) - new_val) * &
                              vol(jl) * 1.0e6_dp
            ! calculate new species field
            p_xt(jt,jl) = (1._dp-localcover(jl)) * p_xt(jt,jl) + new_val
         end if
      enddo
   enddo




 END SUBROUTINE aer_scav_cloud

!===============================================================================
!!$! function for calculation of Reynolds number of aerosol particles
!!$  ELEMENTAL REAL(dp) FUNCTION reynolds(radius, velo, rho, visc)
!!$    
!!$    REAL(dp), INTENT(IN) :: radius, velo, visc, rho
!!$    
!!$! radius in m, velo in m/s, visc in Pa s, rho in kg/m^3
!!$!           m * m/s * kg/m^3 / kg/(m*s)  
!!$!           kg/(m*s) / kg/(m*s) 
!!$    reynolds = radius * velo * rho / visc
!!$
!!$  END FUNCTION reynolds
!------------------------------------------------------------------------------------------
! function for calculation of Reynolds number for rain droplets
  ELEMENTAL  REAL(dp) FUNCTION rey_rain(radius, rho, visc)

! following Pruppacher&Klett, 2000, p.417 eq. 10.142 - 10.145
! parametrization following Beard & Pruppacher, 1969
   use messy_main_constants_mem,    ONLY: g

   REAL(dp), INTENT(IN) :: radius, visc, rho
   REAL(dp) :: Cd_Re2, para1, para2, para3
   REAL(dp), PARAMETER :: b0 = -0.318657e1_dp
   REAL(dp), PARAMETER :: b1 = 0.992696_dp
   REAL(dp), PARAMETER :: b2 = -0.153193e-2_dp
   REAL(dp), PARAMETER :: b3 = -0.987059e-3_dp
   REAL(dp), PARAMETER :: b4 = -0.578878e-3_dp
   REAL(dp), PARAMETER :: b5 = 0.855176e-4_dp
   REAL(dp), PARAMETER :: b6 = -0.327815e-5_dp
   REAL(dp), Parameter :: rho_wat = 1000._dp           ! in kg/m^3

   para1 = 1._dp /(3._dp * visc * visc)
   Cd_Re2 = 32._dp * radius * radius * radius * rho * g * rho_wat 
   Cd_Re2 = Cd_Re2 * para1
   
   para3 = log(Cd_Re2)

   para2 =  para3*b6 + b5
   para2 = para2 * para3 + b4
   para2 = para2 * para3 + b3
   para2 = para2 * para3 + b2
   para2 = para2 * para3 + b1
   para2 = para2 * para3 + b0

   rey_rain = exp(para2)

 END FUNCTION rey_rain
!===============================================================================
! function for calculation of Reynolds number for snow crystals
  ELEMENTAL  REAL(dp) FUNCTION rey_snow(radius, rho, visc)

! following Pruppacher&Klett, 2000, p.427 eq. 10.166 and 10.161-10.162
   use messy_main_constants_mem,    ONLY: g

   REAL(dp), INTENT(IN) :: radius, visc, rho
   REAL(dp) :: Cd_Re2, para1, para2, para3
   REAL(dp), PARAMETER :: b0 = -1.3247_dp
   REAL(dp), PARAMETER :: b1 = 1.0396_dp
   REAL(dp), PARAMETER :: b2 = -0.047556_dp
   REAL(dp), PARAMETER :: b3 = -0.002327_dp
   REAL(dp), Parameter :: rho_wat = 1000._dp           ! in kg/m^3

   para1 = 1._dp /(3._dp * visc * visc)
   Cd_Re2 = 32._dp * radius * radius * radius * rho * g * rho_wat 
   Cd_Re2 = Cd_Re2 * para1
   
   para3 = log(Cd_Re2)

   para2 = para3 * b3 + b2
   para2 = para2 * para3 + b1
   para2 = para2 * para3 + b0

   rey_snow = exp(para2)

 END FUNCTION rey_snow
!===============================================================================
! function for calculation of Schmidt number
  ELEMENTAL REAL(dp) FUNCTION schmidt(visc, rho, diff)
    
    REAL(dp), INTENT(IN) :: visc, rho, diff
!      rho in kg/m^3, visc in Pa s, Diff in m^2/s
!          N/m^2 * s / (kg/m^3*m^2/s)
!          (kg m s / (s^2*m^2)) / (kg * m^2 /(m^3*s))
!          kg /(m*s) / kg  /(m*s)

 
   schmidt = visc / (rho * diff)

  END FUNCTION schmidt
!===============================================================================
! function for calculation of Stokes number
 ! ELEMENTAL 
  REAL(dp) FUNCTION stokes(velo_r, velo_p, radius, visc, rad_rain, lam) 

  
    REAL(dp), INTENT(IN) :: velo_r, velo_p, radius, visc, rad_rain, lam
    REAL(dp)             :: tau, param
    REAL(dp), parameter  :: rho_wat = 1000._dp ! in kg/m^3
 
    ! relaxation time
    tau = 2._dp * radius * radius * rho_wat /(9._dp*visc) *CUN(radius,lam) 
 
    param = velo_r - velo_p
    param = tau * param

    stokes = param / rad_rain    

  END FUNCTION stokes

!===============================================================================
! function for calculation of the dynamic viscosity of air in Pa s

  ELEMENTAL REAL(dp) FUNCTION airvisc(temp)
    
    REAL(dp), INTENT(IN) :: temp

    airvisc = 1.8274e-5_dp * (temp/293.15)**0.74

  END FUNCTION airvisc
!===============================================================================

! function for calculation of the mean free path length 
!     lambda = 6.6e-6 * (1013.25 / press) * (temp/293.15)

  ELEMENTAL REAL(dp) FUNCTION lambda(temp, press)

     REAL(dp), INTENT(IN) :: temp, press
     
     lambda =  2.28E-5_DP * TEMP / PRESS

   END FUNCTION lambda

!===============================================================================
! function for calculation of the D = kT(1+alpha*Kn) / (6*pi*eta*r) 

   ELEMENTAL  REAL(dp) FUNCTION DiffCoeff(temp, rad, lam, visc)

!  units:   J/K * K / kg/(m*s) * m =
!           = Nm / kg/s = kg m^2/s^2 /kg/s =
!           = m^2 / s

    USE messy_main_constants_mem, ONLY: pi, R_gas, avo => N_A

     REAL(dp), INTENT(IN) :: temp, rad, lam, visc
     REAL(dp) :: z1, z2
         
     z1 = (R_gas/avo) * temp * CUN(rad,lam)          ! divident, 
     z2 = 6._dp * pi * visc * rad                    ! divisor
     
     DiffCoeff = z1/z2


   END FUNCTION DIFFCOEFF
!===============================================================================
! function for the cunningham correction
    ELEMENTAL REAL(dp) FUNCTION CUN(rad, lam)

      REAL(dp), INTENT(IN) :: rad, lam
      REAL(dp) :: Kn, alpha

      Kn = lam / rad                                   ! Knudsen number
      alpha = 1.257_dp +  0.4_dp*exp(-1.1_dp/Kn)       ! cunningham factor
      CUN = 1._dp + (alpha * Kn)

    END FUNCTION CUN
!===============================================================================
! function for terminal velocity of particles
   ELEMENTAL REAL(dp) FUNCTION TVELO(radius, visc, lam)

     use messy_main_constants_mem,    ONLY: g

     REAL(dp), INTENT(IN) :: radius, visc, lam
     REAL(dp) :: velo
     REAL(dp), PARAMETER :: rho_part = 1.e3_dp   ! in kg/m^3

     ! m/s^2 * m * m *kg/m^3 / (kg/(m*s))
     ! kg/s^2 / kg/(m*s) -> m/s

     velo = 4._dp * g * cun(radius, lam) * radius * radius * rho_part
     velo = velo /(18._dp*visc)

     tvelo = velo

   END FUNCTION TVELO

! function for terminal velocity
! --- from knwon Reynolds numbers
   ELEMENTAL REAL(dp) FUNCTION RVELO(radius, visc, rey_r, rho)
 
     REAL(dp), INTENT(IN) :: radius, visc, rey_r, rho
     REAL(dp) :: velo

     ! kg/(m*s) / (m * kg/m^3)
     ! kg/kg * m^2/(m*s) -> m/s

     velo = visc * rey_r 
     rvelo = velo /(2._dp * radius * rho)

   END FUNCTION RVELO


!===============================================================================
!===============================================================================

 END MODULE messy_scav_aer

