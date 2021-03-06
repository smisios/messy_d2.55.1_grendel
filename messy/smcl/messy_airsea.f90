
MODULE messy_airsea

  !-------------------------------------------------------------------------------------------------------
  !  MESSy- submodel for Gas transfer at water surface.
  !
  !  AUTHOR:  Pozzer Andrea, MPICH, Oct 2004
  !-------------------------------------------------------------------------------------------------------
  ! WARNINGS:
  ! 1) Works only for eqilibrium for VOC (or other tracers)
  !    !! If you want only deposition or emission refer to other submodel 
  !-------------------------------------------------------------------------------------------------------


  ! MESSy
  USE messy_main_constants_mem,  ONLY: DP, SP, STRLEN_MEDIUM

  IMPLICIT NONE 
  PRIVATE

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'airsea'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.0'

  PUBLIC :: dp, sp, STRLEN_MEDIUM

  ! CTRL-NAMELIST PARAMETERS
  ! switch for dfferent parametrisation of exchange velocity! 
  INTEGER, PUBLIC :: param_kw    ! GLOBAL INTEGER for transfer velocity (in water) calulation 
  ! switch for writing change directly in tendency
  LOGICAL, PUBLIC :: l_tendency = .false.   ! GLOBAL SWITCH for tendency calculation
  ! switch for calculation of whitcap coverage
  LOGICAL, PUBLIC :: l_whitecap = .false.   ! GLOBAL SWITCH for whitecap calculation
  ! switch for calculation of rain effect
  LOGICAL, PUBLIC :: l_rain     = .false.   ! GLOBAL SWITCH for rain calculation
  ! switch for alternative calculation  of turbulence (like drydep) 
  LOGICAL, PUBLIC :: l_turb     = .false.   ! GLOBAL SWITCH for salinity effect calculation
  ! switch for input of salinity climatology
  LOGICAL, PUBLIC :: l_salt     = .false.   ! GLOBAL SWITCH for salinity effect calculation

  PUBLIC :: airsea_read_nml_ctrl
  PUBLIC :: airsea_oro          
  !PUBLIC :: airsea_layerthickness    ! function
  PUBLIC :: airsea_calc_density 
  PUBLIC :: airsea_calc_cair
  PUBLIC :: airsea_calc_schmidt_air
  PUBLIC :: airsea_calc_schmidt_sea
  PUBLIC :: airsea_calc_ostwald
  PUBLIC :: airsea_calc_kl_sea           
  PUBLIC :: airsea_calc_kl_sea_wc           
  PUBLIC :: airsea_calc_kl_rain           
  PUBLIC :: airsea_calc_kl_air           
  PUBLIC :: airsea_calc_kl_air_special           
  PUBLIC :: airsea_calc_kl_tot           
  PUBLIC :: airsea_calc_wc           
  PUBLIC :: airsea_calc_henry           
  PUBLIC :: airsea_water_conc
  PUBLIC :: airsea_delta_conc
  PUBLIC :: airsea_flux
  PUBLIC :: calc_r1
  PUBLIC :: calc_r2
  PUBLIC :: calc_r1_turb
  PUBLIC :: calc_r2_turb
  PUBLIC :: check_stability

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_read_nml_ctrl(status, iou)

    !  AIRSEA MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Pozzer Andrea, MPICH, Oct 2004


    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CTRL/ l_tendency,      &
                    l_whitecap,      &
                    l_rain,          &
                    l_turb,          &
                    l_salt,          &
                    param_kw

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='airsea_read_nml_ctrl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    l_tendency = .false.
    l_whitecap = .false.
    l_turb     = .false.
    l_rain     = .false.
    l_salt     = .false.
    param_kw   = 2

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR

  END SUBROUTINE airsea_read_nml_ctrl

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_calc_kl_tot(kl_sea,kl_air,kl_tot,alpha,henry,temp,proma,seaice,losea)

  IMPLICIT NONE 
  INTRINSIC SIZE
  REAL(dp), INTENT(in)  :: kl_sea(:),kl_air(:) 
  REAL(dp), INTENT(in)  :: alpha 
  REAL(dp), INTENT(in)  :: henry(:)
  REAL(dp), INTENT(in)  :: temp(:)
  INTEGER, INTENT(in)   :: proma
  REAL(dp), INTENT(out) :: kl_tot(:) 
  REAL(dp), INTENT(in)   :: losea(:)       ! land-sea mask
  REAL(dp), INTENT(in)   :: seaice(:)      ! sea-ice scaling
  INTEGER :: il
  REAL(dp) :: H(SIZE(henry))

  ! H is the dimensionless henry constant! 
  
  H(:)=(henry(:)*temp(:))/12.2_dp

    DO il=1,proma
      IF (losea(il).lt.0.5 .AND. kl_sea(il)/=0._dp .AND. kl_air(il) /=0.0_dp) THEN
         SELECT CASE (param_kw)
         CASE(1:5) 
            kl_tot(il)=(1._dp/(alpha*kl_sea(il)))+(H(il)/kl_air(il))
            kl_tot(il)=(1._dp/kl_tot(il))*(1._dp-seaice(il))
! Modification by Krysztofiak for the case Marandino (2009) and Bell (2013)
         CASE(6:7)
            kl_tot(il)=kl_sea(il)
         END SELECT
      ELSE
       kl_tot(il)=0._dp
      END IF
    END DO
  END SUBROUTINE airsea_calc_kl_tot

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_calc_kl_air(kg_vel,zust,u,v,sc_air,proma,losea)           

  IMPLICIT NONE 
  INTRINSIC SQRT
  INTRINSIC SIZE
  REAL(dp), INTENT(out)  :: kg_vel(:) 
  REAL(dp), INTENT(IN)  :: zust(:)        ! friction velocity 
  REAL(dp), INTENT(IN)  :: sc_air(:)      ! Schmidt number in air 
  REAL(dp), INTENT(IN)  :: u(:)           ! wind velocity U
  REAL(dp), INTENT(IN)  :: v(:)           ! wind velocity V
  REAL(dp), INTENT(in)  :: losea(:)       ! land-sea mask
  INTEGER, INTENT(in)   :: proma        
  REAL(dp) :: wind(SIZE(u))               ! wind velocity (modulo)
  REAL(dp) :: r1(SIZE(zust))
  REAL(dp) :: r2(SIZE(zust))
  INTEGER :: il

  ! calculation of the exchange constant for gas phase, following 
  ! "resistance" method: Kg_vel=1/(R1+R2) where :
  ! R1 = aerodynamic resistance
  ! R2 = gas-phase film resistance
  ! first parametrisation adopted:
  ! from Carpenter (tested, default one)
  
  wind(:)=SQRT(u(:)*u(:)+v(:)*v(:))
 
  CALL calc_r1(r1,wind,zust)

  CALL calc_r2(r2,zust,sc_air)
  
    DO il=1,proma
      IF (losea(il).lt.0.5) THEN
        kg_vel(il)=(1._dp/(r1(il)+r2(il)))
      ELSE
        kg_vel(il)=0._dp
      END IF
    END DO

  END SUBROUTINE airsea_calc_kl_air
  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_calc_kl_air_special(kg_vel,molweight,cdnw,cfmw,cfncw, riw, tvir,&
                                 tvw,g, az0, zdz,u,v,proma,losea)           
  IMPLICIT NONE 
  INTRINSIC SQRT
  INTRINSIC SIZE
  REAL(dp), INTENT(out)  :: kg_vel(:) 
  REAL(dp), INTENT(IN)  :: cdnw(:)        ! neutral drag coefficient,water
  REAL(dp), INTENT(IN)  :: cfmw(:)        ! momentum drag coefficient, water 
  REAL(dp), INTENT(IN)  :: cfncw(:)       ! exchange paramenter, water
  REAL(dp), INTENT(IN)  :: riw(:)         ! richardson number, water
  REAL(dp), INTENT(IN)  :: u(:)           ! wind velocity U
  REAL(dp), INTENT(IN)  :: v(:)           ! wind velocity V
  REAL(dp), INTENT(in)  :: losea(:)       ! land-sea mask
  REAL(dp), INTENT(in)  :: tvir(:)        ! virtual surface temperature 
  REAL(dp), INTENT(in)  :: tvw(:)         ! virtual surface temperature, water 
  REAL(dp), INTENT(in)  :: g              ! gravity
  REAL(dp), INTENT(IN)  :: az0(:)         ! roughness length 
  REAL(dp), INTENT(IN)  :: zdz(:)         !  thikness lowest layer 
  REAL(dp), INTENT(in)  :: molweight      ! molarweight of tracer
  INTEGER, INTENT(in)   :: proma        

  REAL(dp) :: ustarw(SIZE(u))             ! ustar over water !not from ECHAM5!!!!
  REAL(dp) :: cmw(SIZE(u))                ! needed for ustar
  REAL(dp) :: wind(SIZE(u))               ! wind velocity (modulo)
  REAL(dp) :: r1(SIZE(u))
  REAL(dp) :: r2(SIZE(u))
  INTEGER :: il

  ! calculation of the exchange constant for gas phase, following 
  ! "resistance" method: Kg_vel=1/(R1+R2) where :
  ! R1 = aerodynamic resistance
  ! R2 = gas-phase film resistance
  ! second parametrisation adopted:
  ! from drydep (more consistent with the other code, not tested)
  
  wind(:)=SQRT(u(:)*u(:)+v(:)*v(:))
 
     ! WE follow the same procedure as in drydep
     ! calculation of ustar (even if we have it already form ECHAM5!)
     WHERE (ABS(cfncw(:))> TINY(cfncw(1)))
         cmw(:)=cdnw(:)*cfmw(:)/cfncw(:)
     ELSEWHERE
         cmw(:)=0.
     ENDWHERE
  ustarw(:)=SQRT(cmw(:))*wind(:)
     ! calculation of resistances

  CALL calc_r1_turb(r1,ustarw,riw,tvir,tvw,wind,g,zdz,az0)
  CALL calc_r2_turb(r2,ustarw, molweight)

    DO il=1,proma
      IF (losea(il).lt.0.5) THEN
        kg_vel(il)=(1._dp/(r1(il)+r2(il)))
      ELSE
        kg_vel(il)=0._dp
      END IF
    END DO

  END SUBROUTINE airsea_calc_kl_air_special
  ! ---------------------------------------------------------------------------

  SUBROUTINE calc_r1_turb(r1,zust,riw,tvir,tvw,wind,g,psurf,paz0)           

  IMPLICIT NONE 
  INTRINSIC SIZE
  INTRINSIC MAX
  INTRINSIC LOG
  INTRINSIC ATAN 
  REAL(dp), INTENT(out)  :: r1(:) 
  REAL(dp), INTENT(IN)  :: zust(:)     ! friction velocity 
  REAL(dp), INTENT(IN)  :: riw(:)      !  Richardson number over water 
  REAL(dp), INTENT(IN)  :: tvir(:)     ! virtual surface tempearture 
  REAL(dp), INTENT(IN)  :: tvw(:)      ! virtual surface tempearture, water 
  REAL(dp), INTENT(IN)  :: wind(:)     ! wind  
  REAL(dp), INTENT(IN)  :: g           ! gravity
  REAL(dp), INTENT(IN)  :: psurf(:)    ! = zdz ,layer thikness 
  REAL(dp), INTENT(IN)  :: paz0(:)     ! roughness length 
  REAL(dp)  :: monin(SIZE(zust))       ! Monin-Obukhov lenght
  REAL(dp)  :: zoverl(SIZE(zust))      ! zdz / Monin-Obukhov lenght
  REAL(dp) ::  zxzsurf(SIZE(zust))     ! for psi calculation
  REAL(dp) ::  zxzref(SIZE(zust))      ! for psi calculation
  REAL(dp) ::  psih(SIZE(zust))        ! psi 
  REAL(dp), PARAMETER :: zkmkh=0.74_dp 


   WHERE (riw(:).gt.0._dp) 
      ! calculating the Monin-Obukhov lenght directly applying the
      ! formula given by Stull, 9.7.5k, page 386
       monin(:)=(zust(:)*((tvir(:)+tvw(:))/2.)*             &
            wind(:))/(0.4*g*(tvir(:)-tvw(:)))
       zoverl(:)=psurf(:)/monin(:)
       psih(:)=-4.7*zoverl(:)
    ELSEWHERE
       monin(:)=psurf(:)/riw(:)
       zoverl(:)=psurf(:)/monin(:)
       zxzsurf(:)=zkmkh*(1._dp-9.*(zoverl(:)))**(0.5)
       zxzref(:)=zkmkh
       psih(:)= &
            (2.*LOG((1._dp+zxzsurf(:))/2.)+LOG((1._dp+zxzsurf(:)**2.)/2.)-  &
            2.*ATAN(zxzsurf(:)))-  & ! primitive function value for z
            (2.*LOG((1._dp+zxzref(:))/2.)+LOG((1._dp+zxzref(:)**2.)/2.)- &
            2.*ATAN(zxzref(:)))      ! primitive function value for zz
    ENDWHERE

    r1(:)=MAX(1._dp,(1./(zust(:)*0.4))* &
         (LOG((psurf(:))/paz0(:))-psih(:)))
 
  END SUBROUTINE calc_r1_turb

  ! ---------------------------------------------------------------------------

  SUBROUTINE calc_r2_turb(r2, zust, molweight)           

  IMPLICIT NONE 
  REAL(dp), INTENT(out)  :: r2(:) 
  REAL(dp), INTENT(IN)  :: zust(:)     ! friction velocity 
  REAL(dp), INTENT(IN)  :: molweight  !  molar weight 
  REAL(dp)  :: diff   
  REAL(dp)  :: diffrb   

  INTRINSIC SQRT

  ! AS IN DRYDEP !
  diff=SQRT(molweight/18._dp)
  diffrb=(0.189_dp/(0.212_dp/diff))**(2._dp/3._dp)
  ! calculation of quasi-laminar boundary layer resistances
  r2(:)=(1./(zust(:)*0.40))*diffrb

  END SUBROUTINE calc_r2_turb

  ! ---------------------------------------------------------------------------

  SUBROUTINE calc_r1(r1,wind,zust)           

  IMPLICIT NONE 
  REAL(dp), INTENT(out)  :: r1(:) 
  REAL(dp), INTENT(IN)  :: wind(:)     ! wind  
  REAL(dp), INTENT(IN)  :: zust(:)     ! friction velocity 

  WHERE (ABS(zust(:)) > 0._dp) 
     r1(:)=wind(:)/(zust(:)*zust(:))
  ELSEWHERE
     r1(:)= 1.e-30_dp
  ENDWHERE

  END SUBROUTINE calc_r1

  ! ---------------------------------------------------------------------------

  SUBROUTINE calc_r2(r2,zust,sc_air)           

  IMPLICIT NONE 
  REAL(dp), INTENT(out)  :: r2(:) 
  REAL(dp), INTENT(IN)  :: zust(:)     ! friction velocity 
  REAL(dp), INTENT(IN)  :: sc_air(:)   ! Schmidt number in air

  WHERE (ABS(zust(:)) > 0._dp) 
     r2(:)=(5._dp/zust(:))*(Sc_air(:)**(2._dp/3._dp))
  ELSEWHERE
     r2(:)= 1.e-30_dp
  ENDWHERE
 
  END SUBROUTINE calc_r2

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_calc_kl_sea(wind,sc_sea,k_vel,proma,losea)

  IMPLICIT NONE 
  REAL(dp), INTENT(in)  :: wind(:) 
  REAL(dp), INTENT(in)  :: sc_sea(:) 
  REAL(dp), INTENT(out) :: k_vel(:)
  REAL(dp), INTENT(in)   :: losea(:)
  INTEGER, INTENT(in)   :: proma
  INTEGER :: il
  REAL(dp) :: beta
  INTRINSIC SQRT

    ! UNITY SISTEM: cm/h --> m/s : beta=2.8E-6
    beta=2.8E-6_dp
  
    DO il=1,proma
      IF (losea(il).lt.0.5) THEN
         SELECT CASE (param_kw)
         CASE(1)
         !liss and Merivat (1986)
            IF (wind(il).lt.3.6) THEN
               k_vel(il)=beta*(0.17*wind(il))*((660._dp/sc_sea(il))**(2._dp/3._dp))
            ELSEIF(wind(il).gt.3.6.and.wind(il).lt.13) THEN
               k_vel(il)=beta*(2.85*(wind(il))-9.65)*(SQRT(660._dp/sc_sea(il)))
            ELSEIF(wind(il).gt.13) THEN
               k_vel(il)=beta*(5.9*(wind(il))-49.3)*(SQRT(660._dp/sc_sea(il)))
            ENDIF
         CASE(2)
         !Wanninkhof(1992)
            k_vel(il)=beta*0.31_dp*(wind(il)**2)*(SQRT(660._dp/sc_sea(il)))
         CASE(3)
         !Wanninkhof and McGills (1999)
            k_vel(il)=beta*0.0283_dp*(wind(il)**3)*(SQRT(660._dp/sc_sea(il)))
         CASE(4)
         !Nightingale (2000)
            k_vel(il)=beta*(0.333_dp*wind(il)+0.222_dp*(wind(il)**2)*(SQRT(600._dp/sc_sea(il))))

         CASE(5)
         ! Ho (2006)
            k_vel(il)=beta* 0.266_dp*(wind(il)**2)*(SQRT(660._dp/sc_sea(il)))
         CASE(6)
         ! Marandino (2009)
            if (wind(il).gt.0.24_dp) then
              k_vel(il)=beta*(100._dp/24._dp)*(0.46_dp*wind(il)-0.24_dp)*(SQRT(720._dp/sc_sea(il)))
            else  
              k_vel(il)=0.0_sp
              endif   
         CASE(7)
         ! modified after Bell (2013)
            IF (wind(il).lt.11.0_dp) THEN
               k_vel(il)=beta*(0.333_dp*(wind(il))+0.222_dp*(wind(il)**2))*(SQRT(600._dp/sc_sea(il)))
            ELSE
               k_vel(il)=beta*((2.85_dp*11.0_dp)-9.65_dp)*(SQRT(600._dp/sc_sea(il)))
            ENDIF
         END SELECT
      ELSE
         k_vel(il)=0._dp
      END IF
    END DO 

  END SUBROUTINE airsea_calc_kl_sea

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_calc_kl_sea_wc(wind,sc_sea,k_vel,proma,losea,wc,henry,sst)

  IMPLICIT NONE 
  REAL(dp), INTENT(in)  :: wind(:) 
  REAL(dp), INTENT(in)  :: sc_sea(:) 
  REAL(dp), INTENT(out) :: k_vel(:)
  INTEGER,  INTENT(in)   :: proma
  REAL(dp), INTENT(in)   :: losea(:)
  REAL(dp), INTENT(in)  :: wc(:) 
  REAL(dp), INTENT(in)  :: henry(:) 
  REAL(dp), INTENT(in)  :: sst(:) 
  REAL(dp) :: ostwald(proma)
  REAL(dp) :: kb(proma)
  INTEGER :: il
  INTRINSIC SQRT

    CALL airsea_calc_ostwald(sst,henry,ostwald)

    DO il=1,proma
      IF (losea(il).lt.0.5_dp) THEN
         kb(il)=wc(il)*(-37._dp*(1._dp/ostwald(il))+6120._dp*((1._dp/ostwald(il))**0.37_dp)*((1._dp/sc_sea(il))**0.18_dp))
         k_vel(il) = kb(il) + (47._dp*wind(il)+wc(il)*(115.200_dp-47._dp*wind(il)))*(SQRT(1._dp/sc_sea(il)))
         ! conversion from cm/h to m/s
         k_vel(il) = 2.8E-6_dp*k_vel(il) 
      ELSE
         k_vel(il)=0._dp
      END IF
    END DO

  END SUBROUTINE airsea_calc_kl_sea_wc
  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_calc_kl_rain(k_vel,proma,sc_sea,losea,R)

  IMPLICIT NONE 
  REAL(dp), INTENT(in)  :: R(:) 
  REAL(dp), INTENT(in)  :: sc_sea(:) 
  REAL(dp), INTENT(inout) :: k_vel(:)
  INTEGER,  INTENT(in)   :: proma
  REAL(dp), INTENT(in)   :: losea(:)
  INTEGER :: il
  INTRINSIC SQRT
    DO il=1,proma
      IF (losea(il).lt.0.5_dp.and.R(il).gt.0_dp) THEN
      ! R = mm/h
      ! from cm/h to m/s
         k_vel(il)= k_vel(il)+2.8E-6_dp*(0.929+0.679*R(il)-0.0015*R(il)*R(il))*(SQRT(600._dp/sc_sea(il)))
      END IF
    END DO


  END SUBROUTINE airsea_calc_kl_rain

  ! ---------------------------------------------------------------------------

  subroutine  airsea_calc_density(temp, sphum, press, density)

  USE messy_main_constants_mem, ONLY: R_gas, M_air, M_H2O
    ! kg/m3
    ! this is exactly as cair*Mair/1000

      IMPLICIT NONE 
      REAL(dp), INTENT(in)  :: temp(:)
      REAL(dp), INTENT(in)  :: sphum(:)
      REAL(dp), INTENT(in)  :: press(:)
      REAL(dp), INTENT(out) :: density(:)

! TODO ---> nullify density(:) and check if TEMP>0!!!!
       density (:) = press(:)/(temp(:)*(R_gas*1000._dp/M_air)*&
                        (1._dp+(M_air/M_H2O-1._dp)*sphum(:)))


  end subroutine airsea_calc_density

  ! ---------------------------------------------------------------------------

  SUBROUTINE  airsea_calc_wc(wind, wc, losea, proma)
      IMPLICIT NONE 
      REAL(dp), INTENT(in)  :: wind(:) 
      REAL(dp), INTENT(out) :: wc(:)
      REAL(dp), INTENT(in)  :: losea(:)
      INTEGER, INTENT(in)   :: proma
      INTEGER :: il 

      ! coefficient by Asher ("Gas transfer at water surface" (Donelan-Geophys.Monograph,2002) 
      DO il=1,proma
        IF (losea(il).lt.0.5_dp) THEN
        wc(il)=3.7E-6_dp*((wind(il) - 1.2_dp)**3)
        ELSE
        wc(il)=0._dp
        END IF
      END DO

   END SUBROUTINE airsea_calc_wc
  ! ---------------------------------------------------------------------------

  subroutine  airsea_calc_cair(temp, sphum, press, cair)
      USE messy_main_constants_mem, ONLY: R_gas, M_air, M_H2O

      IMPLICIT NONE 
      REAL(dp), INTENT(in)  :: temp(:)
      REAL(dp), INTENT(in)  :: sphum(:)
      REAL(dp), INTENT(in)  :: press(:)
      REAL(dp), INTENT(out) :: cair(:)

      ! mol/m3
       cair (:) = press(:)/(temp(:)*R_gas*(1._dp+(M_air/M_H2O-1._dp)*sphum(:)))

  end subroutine airsea_calc_cair

  ! ---------------------------------------------------------------------------

!!$   elemental real(dp) function airsea_layerthickness(geopot)
!!$
!!$      USE messy_main_constants_mem, ONLY: g
!!$      real(dp), INTENT(IN) :: geopot
!!$
!!$      airsea_layerthickness = geopot/g
!!$
!!$   end function airsea_layerthickness

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_oro(losea,sealandmask,seaice,alake)

      IMPLICIT NONE 
      REAL(dp), INTENT(OUT) :: losea(:)
      REAL(dp), INTENT(IN)  :: sealandmask(:), seaice(:), alake(:)

      ! losea(:)= 1.-((1._dp-sealandmask(:)-alake(:))*(1._dp-seaice(:)))
      losea(:) = 1.-(1._dp-sealandmask(:)-alake(:))
      !  frac_water(:) = (1._dp-sealandmask(:))*(1._dp-seaice(:))
      !  frac_ice(:)   = 1._dp - sealandmask(:)-frac_water(:)

  END SUBROUTINE airsea_oro

  ! ---------------------------------------------------------------------------

  SUBROUTINE  airsea_flux(delta_c, kl, flux)
! unit= molx / (m^2 s)
      IMPLICIT NONE 
      REAL(dp), INTENT(IN) :: delta_c(:)
      REAL(dp), INTENT(IN) :: kl(:)
      REAL(dp), INTENT(OUT) :: flux(:)

      flux(:) = kl(:)*delta_c(:)

  END SUBROUTINE airsea_flux

  ! ---------------------------------------------------------------------------

  SUBROUTINE  check_stability(kl, dz, ztmst, proma, check)

      IMPLICIT NONE 
      INTRINSIC SIZE
      REAL(dp), INTENT(IN) :: kl(:)
      REAL(dp), INTENT(IN) :: dz(:)
      REAL(dp), INTENT(IN) :: ztmst
      LOGICAL, INTENT(OUT) :: check 
      REAL(dp) ::limit(SIZE(dz))
      INTEGER, INTENT(in)   :: proma
      INTEGER :: il
   
      check=.FALSE.
      DO il=1,proma
         limit(il)=dz(il)/ztmst 
         IF (kl(il) .gt. limit(il)) check=.TRUE.      
      END DO

  END SUBROUTINE check_stability
  ! ---------------------------------------------------------------------------

  SUBROUTINE  airsea_water_conc(cw, proma, WATER_CON_CONST)  
! here you can call other subroutine with the calcolations of the concentration
! in the water due to clorophyll or whatsoever...

      IMPLICIT NONE 
      REAL(dp), INTENT(OUT) :: cw(:)
      INTEGER,  INTENT(in)  :: proma
      REAL(dp), INTENT(IN)  :: WATER_CON_CONST
      INTEGER :: il

! concentration in nM ---> [ nmol / dm^3 ] = [ 10^-9 mol / liter]
! the conversion is done in airsea_delta_conc
      DO il=1,proma
         cw(il) = WATER_CON_CONST
      END DO

  END SUBROUTINE airsea_water_conc

  ! ---------------------------------------------------------------------------

  SUBROUTINE  airsea_calc_schmidt_sea(sst,sc_sea,MOL_VOL)

      IMPLICIT NONE 
      INTRINSIC :: SIZE, MIN, MAX
      REAL(dp), INTENT(OUT) :: sc_sea(:) ! schmidt number of tracer in water
      REAL(dp), INTENT(IN) :: MOL_VOL                ! 
      REAL(dp), INTENT(IN) :: sst(:)     ! sea surface temperature
      REAL(dp) :: temp(SIZE(sst))        ! conversion

      !schmidt number for CO2, sst in ?C !!!!!
      temp(:) = sst(:) - 273.15_dp 
      temp(:) = MAX(temp(:),0.0_dp)
      temp(:) = MIN(temp(:),40.0_dp)

      sc_sea(:)  = 2073.1_dp-125.62_dp*temp(:)+3.6276_dp*temp(:)*temp(:)-0.043219_dp*temp(:)*temp(:)*temp(:)

      ! correction factor due to different molar volume (see Hayduk-1974)
      ! we use the formula of Wilke and Chang (1955)
      sc_sea(:) =  sc_sea(:)*((MOL_VOL/37.3_dp)**0.6_dp) 
        

  END SUBROUTINE airsea_calc_schmidt_sea

  ! ---------------------------------------------------------------------------

  SUBROUTINE  airsea_calc_schmidt_air(temp,sc_air,MOL_VOL,molw,press)

      IMPLICIT NONE 
      INTRINSIC SIZE
      INTRINSIC SQRT
      REAL(dp),INTENT(IN)  :: MOL_VOL              !
      REAL(dp),INTENT(IN)  :: temp(:)              ! temperature lowest layer
      REAL(dp),INTENT(IN)  :: press(:)             ! press lowest layer
      REAL(dp),INTENT(IN)  :: molw                 ! molecular weight 
      REAL(dp),INTENT(OUT) :: sc_air(:)            ! schmidt number of tracer in air
      REAL(dp)  :: k_vis(SIZE(temp))               ! kinematic viscosity of air
      REAL(dp)  :: dx(SIZE(temp))                  ! Diffusivity of tracer in air 
      REAL(dp),PARAMETER :: air_mol_mass=28.97_dp  ! Molar mass of air (g/mol)
      REAL(dp),PARAMETER :: air_mol_vol=20.1e-6_dp ! Molar volume of air (m3/mol)

      REAL(dp) :: tracer_mol_vol

     ! viscosity of air
     k_vis(:) = -1.15E-14_dp*temp(:)*temp(:)*temp(:)+9.5728E-11_dp*temp(:)*temp(:)+3.7604E-8_dp*temp(:)-3.4484E-6_dp
     ! FSG METHOD for diffusivity calculation: see Lyman et al, 1982 (8.11.2004)
     tracer_mol_vol = 0.875*MOL_VOL
     dx(:) =((10.E-11_dp*(temp(:)**1.75_dp))*(SQRT((air_mol_mass+molw)/(air_mol_mass*molw))))
     dx(:) = dx(:)/((press(:)/101325._dp)*((air_mol_vol**(1._dp/3._dp))+(1e-6_dp*tracer_mol_vol**(1._dp/3._dp))**2._dp))
     sc_air(:) =  k_vis(:)/dx(:) 



  END SUBROUTINE airsea_calc_schmidt_air

  ! ---------------------------------------------------------------------------

  SUBROUTINE  airsea_calc_ostwald(sst,henry,ost)

      IMPLICIT NONE 
      REAL(dp), INTENT(IN)  :: sst(:)
      REAL(dp), INTENT(IN)  :: henry(:)
      REAL(dp), INTENT(OUT) :: ost(:)

      ! sst in Kelvin 
      ! Ostwald number (dimensionless) 
      ! unit sistem:  K * dm3 * atm * mol-1 * K-1 * mol * atm-1 * dm-3 
      ost(:) = sst(:) * 0.082054_dp * henry(:) 

  END SUBROUTINE airsea_calc_ostwald

  ! ---------------------------------------------------------------------------

  SUBROUTINE  airsea_delta_conc(cx,kh,cw,press,delta_c,losea,ind,EFFECT,l_saturation,liplev)

      IMPLICIT NONE 
      REAL(dp), INTENT(IN)  :: cx(:), cw(:)
      REAL(dp), INTENT(IN)  :: kh(:)
      REAL(dp), INTENT(IN)  :: press(:)
      LOGICAL,  INTENT(IN)  :: EFFECT
      REAL(dp), INTENT(OUT) :: delta_c(:)
      REAL(dp), INTENT(in)  :: losea(:)
      INTEGER, INTENT(in)   :: ind
      LOGICAL, INTENT(IN)   :: l_saturation
      LOGICAL, INTENT(in),OPTIONAL   :: liplev(:) 
      INTEGER il
      !   from  nM   [ 10^-9 mol /dm^3 ]  to [ mol / m^3 ]
      !    we have : 1.E-6_dp * cw(:)

      !  from M   [ mol /dm^3 ] = [ mol / L ] =  [ Kmol / m^3]
      !  to [ mol / m^3 ] we have 1.E3_dp * cw(:)

      ! Conc air(mol/L) = Px * H = P * Cx * H
      ! Delta(conc) = (Conc(water)-Conc(air)), unit sistem = [ mol / m^3 ]
      IF (present(liplev)) THEN
        DO il=1,ind
         IF (liplev(il)) THEN
          IF (losea(il).lt.0.5) THEN
            IF (l_saturation) THEN ! Here we calculate the value in the water based on saturation assumption
               delta_c(il) = (press(il)*cx(il)*kh(il)*(1._dp/101.325_dp) * cw(il) & ! saturation value! 
                                    - press(il)*cx(il)*kh(il)*(1._dp/101.325_dp))
            ELSE
               ! if we use water concentration as [nmol/L] : 
               !delta_c(il) = (1.E-6_dp * cw(il) - press(il)*cx(il)*kh(il)*(1._dp/101.325_dp))
               ! if we use water concentration as [mol/L]=[ Kmol / m^3] (hamocc) : 
               delta_c(il) = (1.E3_dp * cw(il) - press(il)*cx(il)*kh(il)*(1._dp/101.325_dp))
            END IF
          ELSE
            delta_c(il)=0._dp
          END IF
         ELSE
           delta_c(il)=0._dp
         END IF
        END DO
      ELSE
       DO il=1,ind
           IF (losea(il).lt.0.5) THEN
            IF (l_saturation) THEN ! Here we calculate the value in the water based on saturation assumption
               delta_c(il) = (press(il)*cx(il)*kh(il)*(1._dp/101.325_dp) * cw(il) & ! saturation value! 
                                    - press(il)*cx(il)*kh(il)*(1._dp/101.325_dp))
            ELSE
               ! if we use water concentration as [nmol/L] : 
               !delta_c(il) = (1.E-6_dp * cw(il) - press(il)*cx(il)*kh(il)*(1._dp/101.325_dp))
               ! if we use water concentration as [mol/L]=[ Kmol / m^3] (hamocc) : 
               delta_c(il) = (1.E3_dp * cw(il) - press(il)*cx(il)*kh(il)*(1._dp/101.325_dp))
            END IF
           ELSE
             delta_c(il)=0._dp
           END IF
         IF (.not.EFFECT) delta_c(il) = 0._dp
       END DO
      END IF

  END SUBROUTINE airsea_delta_conc
  ! ---------------------------------------------------------------------------


  SUBROUTINE  airsea_calc_henry(pa, pb, temp, henry, set_const, molar_volume, salt)

      IMPLICIT NONE 
      REAL(dp), INTENT(IN)    :: pa, pb
      REAL(dp), INTENT(INOUT) :: set_const
      REAL(dp), INTENT(IN)    :: molar_volume
      REAL(dp), INTENT(IN)    :: temp(:)
      REAL(dp), INTENT(IN)    :: salt(:)
      REAL(dp), INTENT(OUT)   :: henry(:)
      INTRINSIC EXP

     ! TO DO: correction factor for salinity!! (see Singh-2003)
     ! 1/298.=3.3557e-3
     henry(:) = pa*EXP(pb*((1._dp/temp(:))-3.3557E-3_dp))
     if (set_const.eq.0) then
       set_const=0.0018*molar_volume
     end if

     henry(:)=henry(:)*EXP(-set_const*salt(:))

  END SUBROUTINE airsea_calc_henry

END MODULE
