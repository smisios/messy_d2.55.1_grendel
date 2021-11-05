PROGRAM airsea

  !  AUTHOR:  Pozzer Andrea, MPICH, Oct 2005
  !
  !  MESSy Interface  for the Submodel kernel


  USE messy_airsea

  ! GLOBAL PARAMETERS
  INTEGER,PARAMETER  :: iou=21    ! I/O unit
  INTEGER, PARAMETER :: idat = 18
  INTEGER, PARAMETER :: kproma=1
  REAL(dp),PARAMETER :: g=9.8

  ! from PHYSC
  REAL(dp),DIMENSION(1),SAVE :: wind10                       ! wind speed at 10 meters
  REAL(dp),DIMENSION(1),SAVE :: u,v                          ! wind speed in the middle of the box 
  REAL(dp),DIMENSION(1),SAVE :: temp, sphum                  ! temperature & humidity!
  ! FROM NAMELIST
  REAL(dp),DIMENSION(1),SAVE :: zust                         ! u star (u*)
  REAL(dp),DIMENSION(1),SAVE :: zdz                          ! height of the lowest level 
  REAL(dp),DIMENSION(1),SAVE :: salt                         ! salinity 
  REAL(dp),DIMENSION(1),SAVE :: presst                       ! pressure at top of boxmodel
  REAL(dp),DIMENSION(1),SAVE :: pressb                       ! pressure at bottom of boxmodel
  REAL(dp),DIMENSION(1),SAVE :: tsurf                        ! temperature at surface!
  REAL(dp),DIMENSION(1),SAVE :: press                        ! pressure at center of boxmodel
  REAL(dp), SAVE  :: NUM_STEPS
  REAL(dp), SAVE  :: DELTA_T 
  ! WORKSPACE
  INTEGER                               :: status

  TYPE :: TYP_LINK 
     CHARACTER(LEN=STRLEN_MEDIUM) :: channel  = ''
     CHARACTER(LEN=STRLEN_MEDIUM) :: object   = ''
  END TYPE TYP_LINK 

  ! CPL-NAMELIST PARAMETERS
  LOGICAL :: L_LG = .FALSE.  ! dry deposition for Lagrangian tracers
  LOGICAL :: L_GP = .FALSE.  ! dry deposition for Gridpoint tracers
  INTEGER :: i_lg_method = 1 ! dry deposition method (LG)

  CHARACTER(LEN=STRLEN_MEDIUM),SAVE :: convect_rain(2) = '' !information on the channels
  CHARACTER(LEN=STRLEN_MEDIUM),SAVE :: large_rain(2) = ''   !information on the channels
  CHARACTER(LEN=STRLEN_MEDIUM), SAVE :: salinity(2) = ''   !information on the channels

  INTEGER :: nasi = 0                  !actual number of tracer  
  INTEGER, PARAMETER :: NMAXNTRAC_ASI =100        !total number of tracer interaction available 

  CHARACTER(LEN=STRLEN_MEDIUM), SAVE     :: ASI_NAME(NMAXNTRAC_ASI)      !name of tracer
  REAL(dp),DIMENSION(NMAXNTRAC_ASI),SAVE :: HENRY_A                      ! value_a for henry number
  REAL(dp),DIMENSION(NMAXNTRAC_ASI),SAVE :: HENRY_B                      ! value_b for henry number
  REAL(dp),DIMENSION(NMAXNTRAC_ASI),SAVE :: ALPHA                        ! enanchement factor due to  
  REAL(dp),DIMENSION(NMAXNTRAC_ASI),SAVE :: MOL_VOL                      ! molar volume at normal boiling point
  REAL(dp),DIMENSION(NMAXNTRAC_ASI),SAVE :: S_CONST                      ! Setschenow constant
  REAL(dp),DIMENSION(NMAXNTRAC_ASI),SAVE :: MOL_MASS                     ! molar mass
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI),SAVE :: USE_MOL_MASS                 ! output forcing  
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI),SAVE :: EFFECT                       ! no effect of submodel on this tracer  
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI),SAVE :: OUTPUT                       ! output forcing  
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI),SAVE :: SATURATION                   ! saturation calculation   
  REAL(dp),DIMENSION(NMAXNTRAC_ASI),SAVE :: WATER_CON_CONST              ! water concentration constant value  
  TYPE(TYP_LINK),DIMENSION(NMAXNTRAC_ASI),SAVE :: WATER_CON_CHN           ! water concentration from channel
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI),SAVE :: WATER_CHN_USE                ! use of the channel for w.conc 
  

  REAL(dp),DIMENSION(NMAXNTRAC_ASI),SAVE :: molweight     ! molar mass of a trcer 

! MAIN PROGRAM
  CALL airsea_initialize

  OPEN(unit=idat, file="output.dat", status='unknown')




  write (*,*) " ------------------------------------------  " 
  WRITE(*,*) "  Execution.. "
  do i=1,INT(num_steps)
    CALL airsea_physc(i)
    CALL airsea_process
    IF (MODULO(i,INT(num_steps/10)) .eq. 0) WRITE(*,*) "    " ,INT((i/num_steps)*100),"%" 
  end do
  write (*,*) " ------------------------------------------  " 
  write (*,*) " End execution " 






  CLOSE(idat)

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_initialize

    ! AIRSEA MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Pozzer Andrea, MPICH, Oct 2004


    IMPLICIT NONE

    ! DEFINE NCREGRID EVENT TRIGGERs
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='airsea_initialize'
    INTEGER                 :: jt

    WRITE(*,*) 'MEMORY INITIALIZATION'

    ! INITIALIZE CTRL
       ! *** CALL CORE ROUTINE:
       CALL airsea_read_nml_ctrl(status, iou)
       IF (status /= 0) STOP 


!  INIT
    DO jt=1, NMAXNTRAC_ASI
          ASI_NAME(jt)  = ''
          HENRY_A(jt)   = 1.0_dp
          HENRY_B(jt)   = 0.0_dp
          ALPHA(jt)     = 1.0_dp
          MOL_VOL(jt)   = 0.3_dp
          S_CONST(jt)   = 0.0_dp
          MOL_MASS(jt)  = 0.0_dp
          EFFECT(jt)    = .FALSE.
          OUTPUT(jt)    = .FALSE.
          WATER_CON_CONST(jt) = -1.0_dp
    END DO

    ! INITIALIZE CPL
       CALL airsea_read_nml_cpl(status, iou)
       IF (status /= 0) STOP
       ! COUNT TRACERS
       DO jt=1, NMAXNTRAC_ASI
          IF (TRIM(ASI_NAME(jt)) == '') CYCLE
          nasi = nasi + 1
          !
          ASI_NAME(nasi)          = TRIM(ASI_NAME(jt))
          HENRY_A(nasi)           = HENRY_A(jt)   
          HENRY_B(nasi)           = HENRY_B(jt)
          ALPHA(nasi)             = ALPHA(jt)
          MOL_VOL(nasi)           = MOL_VOL(jt)
          S_CONST(nasi)           = S_CONST(jt)
          MOL_MASS(nasi)          = MOL_MASS(jt)
          USE_MOL_MASS(nasi)      = USE_MOL_MASS(jt)
          OUTPUT(nasi)            = OUTPUT(jt)
          EFFECT(nasi)            = EFFECT(jt)
          SATURATION(nasi)        = SATURATION(jt)
          WATER_CON_CONST(nasi)   = WATER_CON_CONST(jt)
       END DO
    write (*,*) " ------------------------------------------  " 
    write (*,*) " Total number of tracer in airsea_nml = " , nasi 
    write (*,*) " ------------------------------------------  " 
    write (*,*) " In the box model this set up is forced :    " 
    write (*,*) " l_salt = FALSE  --> salinity constant       " 
    write (*,*) " l_rain = FALSE  --> no rain included (yet)  " 
    write (*,*) " ------------------------------------------  " 

   ! NAMELIST PRESENT ONLY FOR BOXMODEL
       CALL airsea_box_read_nml_cpl(status, iou)

  END SUBROUTINE airsea_initialize 

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_physc(i)


  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i

  u(1) = 5+5*SIN(REAL(i))       ! wind m/s
  v(1) = 5+5*SIN(REAL(i))       ! wind m/s
  wind10(1)=10+10*SIN(REAL(i)) !10 m/s wind speed
  sphum(1) =0      ! no humiodity (qm1+qtte*dt)

  ! FROM NAMELIST.... TO BE CHANGED OPTIONALLY
  !salt(1) = 0.4_dp ! salinity (ocean average) in mol/L
  !temp(1) = 273    ! temperature (tm1+tte*dt)
  !press(1) =  100000 !Pa
  !presst(1) =  90000 !Pa
  !pressb(1) = 101325 !Pa
  !tsurf(1) = 280 ! temperature at surface
  !zust(1)=0.2      ! u*





  END SUBROUTINE airsea_physc

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_process
 

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr='airsea_process'
    REAL(dp),DIMENSION(1,nasi) :: conc_water              ! water concentration nmol/L
    REAL(dp),DIMENSION(1,nasi) :: conc_air                ! air concentration nmol/L
    REAL(dp),DIMENSION(1,nasi) :: schmidt_air             ! schmidt number (adimensional) in air 
    REAL(dp),DIMENSION(1,nasi) :: schmidt_sea             ! schmidt number (adimensional) in sea
    REAL(dp),DIMENSION(1,nasi) :: henry_val               ! henry number (calculated)
    REAL(dp),DIMENSION(1,nasi) :: kl_sea                  ! water transfer velocity
    REAL(dp),DIMENSION(1,nasi) :: kl_air                  ! air transfer velocity
    REAL(dp),DIMENSION(1,nasi) :: kl_airsea               ! total velocity
    REAL(dp),DIMENSION(1,nasi) :: dc_airsea               ! concentration difference 
    REAL(dp),DIMENSION(1) :: densair                      ! air density!!! [kg/m3]
    REAL(dp),DIMENSION(1) :: zdp                          ! pressure difference for a box
    REAL(dp),DIMENSION(1) :: cair                         ! cair => [mol(air)/m3]
    REAL(dp),DIMENSION(1) :: wc_airsea                    ! whitecap
    REAL(dp),DIMENSION(1,nasi) :: flux_airsea             ! flux (mol(trac) m-2 s-1)
    REAL(dp),DIMENSION(1,nasi) :: densflux_airsea         ! flux (mol(trac) kg(air) mol-1(air) m-2 s-1)
    REAL(dp),DIMENSION(1,nasi) :: zflux_airsea            ! flux (mol(trac) mol-1(air) s-1)
    REAL(dp),DIMENSION(1) :: sea_land                     ! 0=ocean 1=land      
    REAL(dp),DIMENSION(1) :: sea_ice                     ! 0=ocean 1=land      
    LOGICAL          :: check=.FALSE.                     ! logical for stability checking 
    INTEGER :: jl,jn,jt
    INTRINSIC TRIM,SQRT
!--- 0. Initialisations: -----------------------------------------------

    !variable to be defined:
    zdp(1)  = presst(1) - pressb(1) ! Pascal difference in pressure
    check =.FALSE. 
    sea_land(1)=0.   ! open ocean
    sea_ice(1)=0.   ! open ocean

    !INITIALIZATION
    kl_sea(:,:)=0
    kl_air(:,:)=0
    kl_airsea(:,:)=0
    dc_airsea(:,:)=0
    conc_air(:,:)=0                   ! mixing ratio air, GP 

    DO jl=1,nasi 
       molweight(jl)=MOL_MASS(jl)
       IF (molweight(jl).lt.1._dp) THEN 
          WRITE(*,*)' molar mass not define in airsea namelist!'
          STOP
       END IF
    END DO

!==========================================================================================================           
!         COMMON PART
!==========================================================================================================           

!--- 1. MAIN PART: calculation of the parametrisation 
!----We avoid calculation in the first time step, due to not definition of wind speed and temperature
     
     ! land cover ---> orography :
     CALL airsea_calc_density(temp,sphum, press,densair)
     CALL airsea_calc_cair(temp,sphum,press,cair)
     CALL airsea_calc_wc(wind10,wc_airsea,sea_land,1)
!==========================================================================================================           
!            LOOP FOR THE TRACER DEFINED IN THE KERNEL
!==========================================================================================================           
        DO jl=1,nasi 

           conc_water(:,jl)=WATER_CON_CONST(jl)

           call airsea_calc_henry(henry_a(jl),henry_b(jl),tsurf,henry_val(:,jl), &
                                  s_const(jl),mol_vol(jl),salt(:))  
           call airsea_calc_schmidt_sea(tsurf,schmidt_sea(:,jl),MOL_VOL(jl)) 
           call airsea_calc_schmidt_air(temp,schmidt_air(:,jl),MOL_VOL(jl),molweight(jl), &
                                        press) 
           IF (l_whitecap) THEN
             call airsea_calc_kl_sea_wc(wind10,schmidt_sea(:,jl),kl_sea(:,jl), &
                     kproma, sea_land,wc_airsea,henry_val(:,jl), tsurf)
           ELSE
             call airsea_calc_kl_sea(wind10,schmidt_sea(:,jl), kl_sea(:,jl),   &
                     kproma, sea_land)
           END IF
           !IF (l_rain) THEN
           !   call airsea_calc_kl_rain(kl_sea(:,jl)%ptr(,jrow),kproma,schmidt_sea(,jl),                &
           !                            sea_land(,jrow),rain)
           !END IF
           !IF (l_turb) THEN
           !    call airsea_calc_kl_air_special(kl_air(jl)%ptr(,jrow),molweight(jl),    &
           !                        cdnw(,jrow),cfmw(,jrow),cfncw(,jrow), riw(,jrow),  &
           !                        tvir(,jrow), tvw(,jrow),g,az0(,jrow),zdz,        &
           !                        um1(,nlev,jrow),vm1(,nlev,jrow),kproma, sea_land(,jrow))  
           !ELSE
               call airsea_calc_kl_air(kl_air(:,jl),zust,u,v,schmidt_air(:,jl),kproma,sea_land)  
           !ENDIF
           call airsea_calc_kl_tot(kl_sea(:,jl),kl_air(:,jl),kl_airsea(:,jl) &
                ,ALPHA(jl),henry_val(:,jl),tsurf,kproma,sea_ice,sea_land)
                            
           call check_stability(kl_airsea(:,jl),zdz,delta_t,kproma,check)

           IF (check) then  
                WRITE(*,*) 'airsea transfer velocity too big! ---> process instable !'
                STOP
           END IF


!---- until here tracer_indipendent!
!==========================================================================================================           
!           GRID POINT PART
!==========================================================================================================           

           call airsea_delta_conc(conc_air(:,jl),henry_val(:,jl),conc_water(:,jl), &
                        presst, dc_airsea(:,jl),sea_land, kproma,EFFECT(jl),saturation(jl))
           call airsea_flux(dc_airsea(:,jl),kl_airsea(:,jl),flux_airsea(:,jl))
           ! mol(trac) m-2 s-1 * Kg m-3 * m3 mol(air)-1 --> mol(trac) mol(air)-1 Kg m-2 s-1
           densflux_airsea(:,jl) = flux_airsea(:,jl)*densair/cair 
           ! m3/mol(air) * m-1 * mol(trac) m-2 s-1 --> mol(trac) mol(air) s-1
           ! zflux_airsea(:,jl) = densflux_airsea(:,jl) * (g/dp)
           ! these three following calculation sould give the same results!!!!
!           zflux_airsea(:,jl) = (1._dp/(cair * zdz)  &
!                                         * flux_airsea(:,jl))
          zflux_airsea(:,jl) = densflux_airsea(:,jl) * (g/zdp(:)) 
!          zflux_airsea(:,jl) = flux_airsea(:,jl)*densair/cair &
!                                      * (g/zdp) 


!==========================================================================================================           
!             CLOSURE OF THE TRACER CYCLE.... 
!==========================================================================================================           
        END DO

    WRITE(idat,*) kl_air(1,1),kl_sea(1,1),kl_airsea(1,1),u(1),v(1),wind10(1),henry_val(1,1) 

END SUBROUTINE airsea_process

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_box_read_nml_cpl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! NAMELIST CPL
    NAMELIST /CPL_BOX/ NUM_STEPS, DELTA_T, ZDZ, SALT, TSURF, &
                    PRESS, TEMP, PRESST,PRESSB, ZUST

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'airsea_box_read_nml_cpl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status



    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CPL_BOX', modstr//"_box")
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_BOX, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_BOX', modstr//"_box")
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    WRITE(*,*) 'GLOBAL SETTINGS OF THE BOX'
    WRITE(*,*) '---------------------------------------------------------'
    WRITE(*,*) 'NUM STEPS        =  ',NUM_STEPS
    WRITE(*,*) 'DELTA T (sec)    =  ',DELTA_T
    WRITE(*,*) 'ZDZ (m)          =  ',ZDZ
    WRITE(*,*) 'ZDZ (m)          =  ',ZDZ
    WRITE(*,*) 'TSURF (K)        =  ',TSURF
    WRITE(*,*) 'PRESS (Pa)       =  ',PRESS
    WRITE(*,*) 'PRESS_TOP (Pa)   =  ',PRESST
    WRITE(*,*) 'PRESS_BOT (Pa)   =  ',PRESSB
    WRITE(*,*) 'TEMP (K)         =  ',TEMP
    WRITE(*,*) 'SALT (mol/L)     =  ',salt
    WRITE(*,*) 'U* (ZUST) (m/s)  =  ',zust
    WRITE(*,*) '---------------------------------------------------------'


    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE airsea_box_read_nml_cpl

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Andrea Pozzer, MPICH, Aug 2003

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    ! switch for skipping calculation of Lagrangian 
    ! rate coefficients..it is local,not broadcasted

    NAMELIST /CPL/ convect_rain,large_rain, salinity, L_GP, L_LG, i_lg_method,     &
            ASI_NAME, HENRY_A, HENRY_B, ALPHA, MOL_VOL,S_CONST,MOL_MASS,           &
            USE_MOL_MASS, EFFECT, OUTPUT, WATER_CON_CONST, WATER_CON_CHN,          &
            SATURATION
         

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='airsea_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

   WRITE(*,*) '---------------------------------------------------------'
   WRITE(*,*) '   THE FOLLOWING OPTION ARE READ BUT NOT USED  IN THE    '
   WRITE(*,*) '                 BOX MODEL VERSION                       '
   WRITE(*,*) '---------------------------------------------------------'
   WRITE(*,*) '                                                         '
   WRITE(*,*) '  convect_rain(1) =   ',convect_rain(1)  
   WRITE(*,*) '  convect_rain(2) =   ',convect_rain(2)  
   WRITE(*,*) '  large_rain(1)   =   ',large_rain(1)     
   WRITE(*,*) '  large_rain(2)   =   ',large_rain(2)     
   WRITE(*,*) '  salinity(1)     =   ',salinity(1)       
   WRITE(*,*) '  salinity(2)     =   ',salinity(2)       
   WRITE(*,*) '  L_GP            =   ',L_GP              
   WRITE(*,*) '  L_LG            =   ',L_LG          
   WRITE(*,*) '  i_lg_method     =   ',i_lg_method                    
   WRITE(*,*) '                                                         '
   WRITE(*,*) '---------------------------------------------------------'




   WRITE(*,*) '.........................................................'
   WRITE(*,*) '           AIR SEA  TRACER DEFINED                       '
   WRITE(*,*) '.........................................................'

    DO jt=1, NMAXNTRAC_ASI
       IF (TRIM(ASI_NAME(jt)) == '') CYCLE

       WRITE(*,*) '  TRACER NO.          ',jt
       WRITE(*,*) '  NAME              = ', TRIM(ASI_NAME(jt))
       WRITE(*,*) '  HENRY_A           = ', HENRY_A(jt)
       WRITE(*,*) '  HENRY_B           = ', HENRY_B(jt)
       WRITE(*,*) '  ALPHA             = ', ALPHA(jt)
       WRITE(*,*) '  MOL_VOL           = ', MOL_VOL(jt)
       WRITE(*,*) '  S_CONST           = ', S_CONST(jt)
       WRITE(*,*) '  MOLAR_MASS        = ', MOL_MASS(jt)
       WRITE(*,*) '  USE_MOLAR_MASS    = ', USE_MOL_MASS(jt)
       WRITE(*,*) '  FORCING OUTPUT    = ', OUTPUT(jt)
       WRITE(*,*) '  WATER CONC.CONST  = ', WATER_CON_CONST(jt)
       !IF (TRIM(WATER_CON_CHN(jt)%channel) .not. '') THEN
       WRITE(*,*) '  WATER CONC.CHN    = ', WATER_CON_CHN(jt)%channel
       WRITE(*,*) '  WATER CONC.OBJ    = ', WATER_CON_CHN(jt)%object
       !END IF
       WRITE(*,*) '  EFFECT ON AIR     = ', EFFECT(jt)
       WRITE(*,*) '  SATURATION        = ', SATURATION(jt)
       WRITE(*,*) '.........................................................'

    END DO 

  END SUBROUTINE airsea_read_nml_cpl

! ===========================================================================
 
END PROGRAM airsea
