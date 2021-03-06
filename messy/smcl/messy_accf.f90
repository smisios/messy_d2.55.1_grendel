!*******************************************************************************
!
! SUBMODEL CORE LAYER (SMCL) routine for MESSy SUBMODEL ACCF
! 
! This submodel is used to calculate ATR20 for
!  - OZONE
!  - METHANE
!  - WATER VAPOR
!  - CONTRAIL
!
! Author: Feijia Yin (TUDELFT), Volker Grewe (DLR), 2017
! 
! Reference: 
! - J. van Manen, "Aviation H2O and NOx climate cose functions
!   based on local weather", MS.c. thesis, Delft University of Technology,
!   2017.
! - E. Irvine, "Contrail algorithm Climate Change Function derivation", 
!   report, 2017
!
!*****************************************************************************

!*****************************************************************************
MODULE messy_accf
!*****************************************************************************

   USE messy_main_constants_mem, ONLY : DP
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC  :: DP
   
   CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'accf'
   CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'
   ! CTRL namelist parameter
   
   LOGICAL, PUBLIC   :: l_cpl_rad  
   LOGICAL, PUBLIC   :: l_cpl_tropop  
   LOGICAL, PUBLIC   :: l_cpl_orbit  
   LOGICAL, PUBLIC   :: l_cpl_contrail
  
   ! coefficients of ozone aCCFs
   REAL(DP), PUBLIC  :: alpha_0 = -5.2e-11_dp
   REAL(DP), PUBLIC  :: alpha_1 = 2.3e-13_dp
   REAL(DP), PUBLIC  :: alpha_2 = 4.85e-16_dp
   REAL(DP), PUBLIC  :: alpha_3 = -2.0e-18_dp
   ! coefficients of methane aCCFs
   REAL(DP), PUBLIC  :: beta_0 = -9.83e-13_dp
   REAL(DP), PUBLIC  :: beta_1 = 1.99e-18_dp
   REAL(DP), PUBLIC  :: beta_2 = -6.32e-16_dp
   REAL(DP), PUBLIC  :: beta_3 = 6.12e-21_dp
   ! coefficients of h2o aCCFs 
   REAL(DP), PUBLIC  :: gama_0 = 4.05e-16_dp
   REAL(DP), PUBLIC  :: gama_1 = 1.48e-16_dp
   ! coefficients of contrail aCCFs night 
   REAL(DP), PUBLIC  :: theta_0_n = 1.0e-10_dp
   REAL(DP), PUBLIC  :: theta_1_n = 0.0073_dp
   REAL(DP), PUBLIC  :: theta_2_n = 0.0107_dp
   REAL(DP), PUBLIC  :: theta_3_n = -1.03_dp
   ! coefficients of contrail aCCFs day
   REAL(DP), PUBLIC  :: theta_0_d = 1.0e-10_dp
   REAL(DP), PUBLIC  :: theta_1_d = -1.7_dp
   REAL(DP), PUBLIC  :: theta_2_d = -0.0088_dp

   PUBLIC :: accf_read_nml_ctrl
   PUBLIC :: accf_atr20_cal

CONTAINS

   !########################################################################
   ! PUBLIC SUBROUTINES
   !########################################################################


!================================================================================
   SUBROUTINE accf_read_nml_ctrl(status,iou)
      !--------------------------------------------------------------------------
      !This subroutine is used to read the CTRL-namelist of the submodel
      !--------------------------------------------------------------------------
      USE messy_main_tools,  ONLY : read_nml_open, read_nml_check, read_nml_close
     
      IMPLICIT NONE
    
      ! I/O
      INTEGER, INTENT(OUT) :: status ! error message
      INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
      
      NAMELIST /CTRL/ l_cpl_rad,l_cpl_tropop,l_cpl_orbit,l_cpl_contrail

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'accf_read_nml_ctrl'
      LOGICAL                     :: lex
      INTEGER                     :: fstat
    
      ! INITIALIZE
      status = 1

      CALL read_nml_open(lex,substr,iou,'CTRL',modstr)
      IF(.not.lex) RETURN
    
      READ(iou,NML=CTRL,IOSTAT=fstat)
      CALL read_nml_check(fstat,substr,iou,'CTRL',modstr)
      IF(fstat /=0)RETURN

      
      CALL read_nml_close(substr, iou, modstr)
      status = 0

   END SUBROUTINE accf_read_nml_ctrl
!================================================================================

!============================================================================
   SUBROUTINE accf_atr20_cal(nlev,kproma,jrow,&
                              geopot,geosf,temp,rhum,&
                              PV, OLR, cossza, dec_sun,&
                              philat, philon,dayofyear,CPC,&
                              ATR20_O3,ATR20_CH4,ATR20_H2O,&
                              ATR20_CONTRAIL,ATR20_CO2)
      
      USE messy_main_constants_mem,   ONLY : solc, pi

      IMPLICIT NONE
   
      INTEGER,                    INTENT(IN) :: kproma, nlev,jrow,dayofyear
      REAL(DP),                   INTENT(IN) :: dec_sun
      REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: temp           ! temperature 
      REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: rhum           ! relative humidity 
      REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: CPC            ! contrail potential coverage
      REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: geopot         ! geopotential 
      REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: PV             ! potential vorticity
      REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: OLR            ! Outgoing Longwave Radiation
      REAL(DP), DIMENSION(:,:),   INTENT(IN) :: geosf          ! surface geopotential 
      REAL(DP), DIMENSION(:,:),   INTENT(IN) :: cossza         ! cos solar zenith angle
      REAL(DP), DIMENSION(:,:),   INTENT(IN) :: philat         ! latitude 
      REAL(DP), DIMENSION(:,:),   INTENT(IN) :: philon         ! longitude range [0,360]
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: ATR20_O3       ! Average temperature response ozone, [K/kg(NO2)]
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: ATR20_CH4      ! Average temperature response methane, [K/kg(NO2)] 
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: ATR20_H2O      ! Average temperature response water vapor, [K/kg(fuel)]
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: ATR20_CONTRAIL ! Average temperature response contrail, [K/km(contrail)]
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: ATR20_CO2      ! Average temperature response CO2, [K/kg(fuel)]

      INTEGER                                :: jp,jk,sunrise
      REAL(DP)                               :: xlat,lon,costheta,geopot_tot
      REAL(DP)                               :: atr20_ozone,atr20_methane
      REAL(DP)                               :: constant1,constant2
      REAL(DP)                               :: dec_ang        ! declination angel in degree
      REAL(DP)                               :: dec            ! declination angel in radius
      REAL(DP)                               :: F              ! solar flux
      REAL(DP)                               :: theta          ! solar zenith angle in radius
      REAL(DP)                               :: omega_0        ! hourly angel 
      REAL(DP)                               :: loc_time       ! local time
      REAL(DP)                               :: sunrise_time,h_darkness
      
      constant1 = 0.83_dp*0.01745_dp !
      constant2 = pi/2.0_dp         ! to converte 90 degree in radius   
 
      !write(*,*)'day of the year:', dayofyear
      dec_ang = -23.44_dp*COS(360._dp/365._dp*(dayofyear+10._dp))!declination angel in degree
      !write(*,*)'sun dec angel in degree is:',dec_ang
      dec  = dec_ang*pi/180.0_dp   ! declination angel in radius
      !write(*,*)'sun dec angel in radius is:', dec     
      !write(*,*)'sun dec from orbit in radius is:', dec_sun     
 
      DO jk = 1,nlev
        DO jp = 1,kproma
          !calculate water vapor aCCFs
          ATR20_H2O(jp,jk,jrow) = gama_0+gama_1*ABS(PV(jp,jk,jrow))
          ATR20_CO2(jp,jk,jrow) = 1.64e-15_dp 
          xlat = philat(jp,jrow)*pi/180._dp  ! convert latitude to radius
         ! write(*,*)'latitude is:',philat(jp,jrow)
         ! write(*,*)'latitude in radius is:',xlat
          geopot_tot = geopot(jp,jk,jrow)+geosf(jp,jrow)
          
          ! calculate ozone aCCFs
          atr20_ozone = alpha_0+alpha_1*temp(jp,jk,jrow)+ alpha_2*geopot_tot+ &
                                 alpha_3*temp(jp,jk,jrow)*geopot_tot
          IF (atr20_ozone.GE.0._dp) THEN
            ATR20_O3(jp,jk,jrow) = atr20_ozone
          ELSEIF (atr20_ozone .LT. 0._dp) THEN
            ATR20_O3(jp,jk,jrow) = 0._dp
          ENDIF
          !calculate methane aCCFs
          theta = SIN(xlat)*SIN(dec)+COS(xlat)*COS(dec) ! calculate solar zenith angle 
          F     = solc*COS(theta)
         ! write(*,*)'top of atmosphere solar radiation is:',F
          atr20_methane = beta_0+beta_1*geopot_tot+beta_2*F+beta_3*geopot_tot*F
          IF (atr20_methane .LE.0._dp) THEN
            ATR20_CH4(jp,jk,jrow) = atr20_methane
          ELSEIF(atr20_methane .GT. 0._dp)THEN
            ATR20_CH4(jp,jk,jrow) = 0._dp
          ENDIF

          !calculate the contrail aCCFs
          costheta   = cossza(jp,jrow)
!         write(*,*)'cos solar zenith angle is:',costheta  
          IF(CPC(jp,jk,jrow).GT.0._dp)THEN  ! if the contrail is formed 
          
         ! IF((temp(jp,jk,jrow).LT.235._dp).and.(rhum(jp,jk,jrow).GT.100._dp)) THEN  ! if the contrail is formed 
             IF (philon(jp,jrow).GT.180._dp) THEN
                lon = philon(jp,jrow)-360._dp
             ELSE
                lon = philon(jp,jrow)        ! convert longitude between [-180,180]  
             ENDIF
!             write(*,*)'longitude is:', lon
             loc_time = 0._dp+lon/15._dp          
!             write(*,*)'local time is:',loc_time
             IF ((xlat.GT.(-constant2-dec_sun)) .and.(xlat.LT.(constant2+dec_sun))) THEN
                sunrise = 1
             ELSE 
                sunrise = 0
                h_darkness = 24._dp
             ENDIF           
         
             IF ((costheta.LT.0._dp).and.(sunrise.EQ.1)) THEN
                omega_0 = ACOS((SIN(-1._dp*constant1)-SIN(xlat)*SIN(dec_sun))/&
                             (COS(xlat)*COS(dec_sun)))
!                write(*,*) 'hour angel at sunrise is:',omega_0
                sunrise_time = 12._dp-(omega_0*180._dp/pi)/15._dp
!                write(*,*)'sunrise time is:', sunrise_time
                h_darkness = sunrise_time-loc_time
!                write(*,*)'darkness hour is:', h_darkness
             ENDIF
          
             IF( ((costheta.LT.0._dp).and.(h_darkness.GE.6._dp)).or.(sunrise.EQ.0) ) THEN
                ATR20_CONTRAIL(jp,jk,jrow)=(theta_0_n*(theta_1_n*(10._dp**(theta_2_n*temp(jp,jk,jrow)))&
                                        +theta_3_n))*0.114_dp
             ELSE
                ATR20_CONTRAIL(jp,jk,jrow)=theta_0_d*(theta_1_d+theta_2_d*OLR(jp,1,jrow))*0.114_dp
             ENDIF
          ELSEIF (CPC(jp,jk,jrow).LE.0._dp) THEN
             ATR20_CONTRAIL(jp,jk,jrow)=0._dp
          ENDIF

         ! write(*,*)'ATR20_O3 = ', ATR20_O3(jp,jk,jrow)  
         ! write(*,*)'ATR20_CH4 = ', ATR20_CH4(jp,jk,jrow)  
         ! write(*,*)'ATR20_H2O = ', ATR20_H2O(jp,jk,jrow)  
         ! write(*,*)'ATR20_CONTRAIL = ', ATR20_CONTRAIL(jp,jk,jrow)  
        ENDDO
      ENDDO

   END SUBROUTINE accf_atr20_cal
!============================================================================

!*****************************************************************************
END MODULE messy_accf
!*****************************************************************************
