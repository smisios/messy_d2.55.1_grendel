MODULE messy_sedi

! Author: Astrid Kerkweg, MPICH-Mainz, April 2004

! This modules calculates aerosol sedimentation
! all formulas can be found in Pruppacher and Klett, 1997
! "Microphyscis of Clouds and Precipitation", Kluwer 
! or Seinfield and Pandis, " Atmospheric Chemistry and Physics",
! Jon Wiley and Sons, Inc., 1998
! or e.g. Pryor, 1999, JGR
USE messy_main_constants_mem,  ONLY: dp


  IMPLICIT NONE
  PRIVATE
  SAVE
  
  INTRINSIC :: EXP, LOG, MIN, PRESENT
  PUBLIC :: DP

  CHARACTER(len=*), PUBLIC, PARAMETER :: modstr = 'sedi'
  CHARACTER(len=*), PUBLIC, PARAMETER :: modver = '2.5'

  INTEGER, PUBLIC :: order = -1
  INTEGER, PUBLIC :: scheme = 1 ! mz_fb_20100831

  ! SUBROUTINES
  PUBLIC :: sedi_read_nml_ctrl
! op_ck_20131218+
!!$  PUBLIC :: density
! op_ck_20131218-
  PUBLIC :: mean_free_path
  PUBLIC :: air_viscosity
  PUBLIC :: calc_vt
  PUBLIC :: sedi_0
  PUBLIC :: sedi_1

CONTAINS

  !===========================================================================
  
  SUBROUTINE sedi_read_nml_ctrl(status, iou)

    ! SEDI MODULE ROUTINE (CORE)
    !
    ! READ SEDI NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Astrid Kerkweg, MPICH, Feb 2002

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: iou   ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='sedi_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    NAMELIST /CTRL/  order, scheme ! mz_fb_20100831 scheme added

    status = 1 ! ERROR ON RETURN

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    CALL read_nml_close(substr, iou, modstr)

    status = 0

  END SUBROUTINE sedi_read_nml_ctrl

  !========================================================================

!-----------------------------------------------------------------------------

elemental real(dp) function mean_free_path(press, temp)

  REAL(dp), INTENT(in) :: press
  REAL(dp), INTENT(in) :: temp

! see pruppacher and klett, 1997 eq. 10-140
  mean_free_path =0.066 *(1.01325E+5/press)*&
                                  (temp/293.15)*1.E-06_dp

! the above is a computational cheaper fit for:
! mean_free_path = 2 * airvisc /(press*sqrt(8*M_air*1.e-3/(pi*R_gas*temp)))

end function mean_free_path

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine air_viscosity(nproma,nlev, temp, airvisc)

  INTEGER , INTENT(in) :: nproma, nlev

  REAL(dp), INTENT(in) :: temp(nproma,nlev)
  REAL(dp), INTENT(out) :: airvisc(nproma,nlev)

  !local
  REAL(dp) :: ctemp(nproma,nlev) ! temperature in degree Celsius
  

  ctemp(1:nproma,1:nlev) = temp(1:nproma,1:nlev) - 273.15_dp

  WHERE (ctemp(1:nproma,1:nlev) >= 0._dp) 
     ! Pruppacher and Klett, 1997, p. 417, eq. 10-141a
     airvisc(1:nproma,1:nlev) = &
          (1.718_dp + 0.0049_dp*ctemp(1:nproma,1:nlev))*1.E-5_dp
  ELSEWHERE
     ! Pruppacher and Klett, 1997, p. 417, eq. 10-141b
     airvisc(1:nproma,1:nlev) = (1.718_dp + 0.0049_dp*ctemp(1:nproma,1:nlev) - &
          1.2E-05_dp*(ctemp(1:nproma,1:nlev)**2))*1.E-5_dp
  END WHERE

end subroutine air_viscosity

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  SUBROUTINE calc_vt(vt , kproma  , nlev                   &     
                        , pradius , pdensaer, plambda_air  &
                        , prhoair , pviscair, sigma ,flag  )

    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE
    
    INTRINSIC ::  TINY, ABS
    
    ! parameter list
    REAL(dp), INTENT(OUT) :: vt(:,:)
    INTEGER , INTENT(IN)  :: kproma, nlev              
    REAL(dp), INTENT(IN)  :: pradius(:,:) 
    REAL(dp), INTENT(IN)  :: pdensaer(:,:)
    REAL(dp), INTENT(IN)  :: plambda_air(:,:)
    REAL(dp), INTENT(IN)  :: prhoair(:,:)
    REAL(dp), INTENT(IN)  :: pviscair(:,:)
    REAL(dp), INTENT(IN)  :: sigma
    INTEGER,  INTENT(IN)  :: flag

    ! LOCAL
    REAL(dp) :: locradius(kproma,nlev)
    REAL(dp) :: Knudsen_num(kproma,nlev) ! Knudsen Number
    REAL(dp) :: alpha(kproma,nlev)       ! 
    REAL(dp) :: slipcorr(kproma,nlev)    ! Cunningham slip-flow correction
    REAL(dp) :: slinnfac                 ! slinn factor

    SELECT CASE(flag)
       ! BIN
    CASE(0)
       locradius(1:kproma, 1:nlev) = pradius(1:kproma, 1:nlev)
       slinnfac =1._dp

       ! MASS
    CASE(1)
       locradius(1:kproma, 1:nlev) = pradius(1:kproma, 1:nlev) &
            * sigma ** ( 3._dp * log(sigma) )

       ! correction term for lognormal distribution with standrad deviation
       ! sigma (see Slinn and Slinn, 1980, Atmos. Env., p. 1015, eq. 8
      slinnfac = sigma ** (2._dp*log(sigma))

       ! NUMBER
    CASE(2)
       locradius(1:kproma, 1:nlev) = pradius(1:kproma, 1:nlev) 

       ! correction term for lognormal distribution with standrad deviation
       ! sigma (see Slinn and Slinn, 1980, Atmos. Env., p. 1015, eq. 8
       slinnfac = sigma ** (2._dp*log(sigma))

    CASE default
    END SELECT
       
    ! calculate Knudsen Number
    WHERE (locradius(:,:) > 0._dp)
       Knudsen_num(:,:) = plambda_air(:,:)/locradius(:,:)
    ELSEWHERE
       Knudsen_num(:,:) = 0._dp
    ENDWHERE
    
    ! calculate Cunningham slip-flow correction 

    ! numbers for alpha from Davies(1945) see Pruppacher page 450
    WHERE ( ABS(Knudsen_num(:,:)) > TINY(Knudsen_num(1,1)))
       alpha(:,:)    = 1.257_dp + 0.4_dp * EXP(-1.1 / Knudsen_num(:,:))
    ELSEWHERE
       alpha(:,:)    = 0._dp
    ENDWHERE
    slipcorr(:,:) = 1._dp + alpha(:,:) * Knudsen_num(:,:)
  
    vt(:,:) = 2._dp/9._dp * locradius(:,:)* locradius(:,:) * g * &
         (pdensaer(:,:)-prhoair(:,:)) / pviscair(:,:)            &
         *  slipcorr(:,:) * slinnfac

  END SUBROUTINE calc_vt

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
SUBROUTINE SEDI_0(kproma, nlev, xtte, mr, v, dz, dpr, dt, sfloss)

  USE messy_main_constants_mem, ONLY: g

  IMPLICIT NONE

  ! I/O
  INTEGER , INTENT(IN)                   :: kproma, nlev
  REAL(DP), DIMENSION(:,:), INTENT(OUT)  :: xtte   ! mixing ratio tendency
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: mr     ! mixing ratio [X]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: v      ! velocity [m/s]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: dz     ! layer thickness [m]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: dpr    ! pressure thickness [Pa]
  REAL(DP),                 INTENT(IN)   :: dt     ! time step [s]
  ! loss at surface [X/s]
  REAL(DP), DIMENSION(:),   INTENT(OUT), OPTIONAL  :: sfloss

  ! LOCAL
  INTEGER  ::  jl,jk
  REAL(DP), DIMENSION(kproma,nlev) :: frac     ! [1/s]
  REAL(DP), DIMENSION(kproma,nlev) :: tend_out ! [X / s]
  REAL(DP), DIMENSION(kproma,nlev) :: flux_out ! [X * kg / m^2 / s]
  REAL(DP), DIMENSION(kproma,nlev) :: flux_in  ! [X * kg / m^2 / s]
  REAL(DP), DIMENSION(kproma,nlev) :: tend_in  ! [X / s]

  ! FRACTION OF BOX 'FALLING DOWN' [1/s]
  ! NOTE: flux-fraction is limited to total box
  DO jk=1, nlev
    DO jl=1,kproma
      frac(jl,jk) = MIN(v(jl,jk)/dz(jl,jk), 1._dp/dt)
    END DO
  END DO
  ! TENDENCY 'OUT' [X / s]
  tend_out(:,:) = frac(:,:) * mr(:,:)
  ! FLUX GOING OUT OF BOX [X * (kg/m^2) / s]
  flux_out(:,:) = tend_out(:,:) * (dpr(:,:) / g)
  ! FLUX GOING INTO BOX [X * (kg/m^2) / s]
  flux_in(:,1) = 0.0_DP
  DO jk=2, nlev
     DO jl=1, kproma
        flux_in(jl,jk)  = flux_out(jl,jk-1)
     END DO
  END DO
  ! TENDENCY 'IN' [X / s]
  tend_in(:,:) = flux_in(:,:) / (dpr(:,:) / g)
  
  xtte(:,:) = tend_in(:,:) - tend_out(:,:)

  IF (PRESENT(sfloss)) sfloss(:) = tend_out(:,nlev) ! [X / s]

END SUBROUTINE SEDI_0
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
SUBROUTINE SEDI_1(kproma, nlev, press, temp, velo, dt, mr &
     , pressi, tend, sfloss)
 
  USE messy_main_constants_mem,  ONLY: M_air, g

  IMPLICIT NONE

!NOTE: the functions below are based on the original code of
!      Joachim Buchholz, MPICH, 2002-2004 
!      it was originally written to calculate the sedimentation in the
!      PSC model. It is here now adopted to aerosol sedimentation
  INTEGER , INTENT(IN)                 :: kproma, nlev  ! number of real layers
  REAL(dp), DIMENSION(:,:),INTENT(in)  :: press  ! pressure  [Pa]
  REAL(dp), DIMENSION(:,:),INTENT(in)  :: temp   ! temperature [K]
  REAL(dp), DIMENSION(:,:),INTENT(in)  :: velo   ! sedimentation velocity [m/s]
  REAL(dp),              INTENT(in)    :: dt     ! integrational time step [s]
  ! amount of substance fraction [X]
  REAL(dp), DIMENSION(:,:),INTENT(in)  :: mr
  ! pressure at box interfaces [Pa]
  REAL(dp), DIMENSION(:,:),INTENT(in)  :: pressi
  REAL(dp), DIMENSION(:,:), INTENT(OUT) :: tend   ! tendency [X /s]
  ! loss at surface [X/s]
  REAL(DP), DIMENSION(:),   INTENT(OUT), OPTIONAL  :: sfloss

  ! LOCAL 
  REAL(dp), DIMENSION(kproma,nlev+1) :: SedStep  ! sedimentation step [Pa]
  LOGICAL,  DIMENSION(kproma,nlev+1) :: val_sedi ! sedimentation happens?
  REAL(dp), DIMENSION(kproma,nlev+1) :: CHANGE   ! change in amount fraction
  REAL(dp), DIMENSION(kproma,nlev+1) :: loc_mr   ! expanded local mixing ratio
  ! expanded pressure (interf.)
  REAL(dp), DIMENSION(kproma,nlev+2) :: loc_pressi
  ! flux of tracer into lowest (= arificial box)
  REAL(dp), DIMENSION(kproma)   :: trac_in
  REAL(dp)                      :: rhelp

  ! INITIALIZE LOCAL VALUES 
  SedStep(:,:)  = 0._dp
  val_sedi(:,1:nlev) = .true.
  val_sedi(:,nlev+1) = .false.

  CHANGE(:,:)   = 0._dp
  loc_mr(1:kproma,1:nlev) = mr(1:kproma,1:nlev)
  loc_mr(1:kproma,nlev+1) = 0._dp 
  loc_pressi(1:kproma,1:nlev+1) = pressi(1:kproma,1:nlev+1)
  loc_pressi(1:kproma,nlev+2)   = 2._dp*pressi(1:kproma,nlev+1) &
       -pressi(1:kproma,nlev)
  trac_in(:)=0._dp

  ! calculated sedimentation step
  SedStep(1:kproma,1:nlev) = &
       mz_sedi_SedStep(dt,press(1:kproma,1:nlev), &
       temp(1:kproma,1:nlev),velo(1:kproma,1:nlev))

  SedStep(1:kproma,1:nlev)=MIN(SedStep(1:kproma,1:nlev), &
      loc_pressi(1:kproma,2:nlev+1)-loc_pressi(1:kproma,1:nlev))

 ! calculate changes in mixing ration due to first order approach
  change(:,:) = sedimentation(kproma,nlev+1,     &
              loc_pressi(1:kproma,1:nlev+1),  loc_pressi(1:kproma,2:nlev+2), &
              SedStep(1:kproma,1:nlev+1),     loc_mr(1:kproma,1:nlev+1),     &
              val_sedi(1:kproma,1:nlev+1)                               )

  ! calculate tendencies
  tend(1:kproma,1:nlev) = CHANGE(1:kproma,1:nlev) / dt

  ! calculate surface loss
  ! --  1. calculate flux into artificial lowest box
  !  trac_in(1:kproma) = CHANGE(1:kproma,nlev+1) /dt &
  ! * 1.E3_dp/M_Air*(loc_pressi(1:kproma,nlev+2)-loc_pressi(1:kproma,nlev+1))/g
  rhelp=1.E3_dp/(dt*M_Air*g)
  trac_in(1:kproma) = CHANGE(1:kproma,nlev+1) * rhelp &
       *(loc_pressi(1:kproma,nlev+2)-loc_pressi(1:kproma,nlev+1))
  ! --  2. calculate flux out of lowest "real" box
  ! --  trac_out(19) == trac_in(20)  => nothing to calculate
  ! --  3. calculate surface loss from flux out of lowest box
  rhelp=M_air*g/1.E3_dp
  IF (PRESENT(sfloss)) &
       sfloss(1:kproma) = trac_in(1:kproma) * rhelp  / &
       (loc_pressi(1:kproma,nlev+1)-loc_pressi(1:kproma,nlev))

END SUBROUTINE SEDI_1

!=========================================================================

ELEMENTAL FUNCTION mz_sedi_SedStep(TimeStep, PRESS, TEMP, vel)
  !-------------
  ! DESCRIPTION:
  !-------------
  ! The function calculates the vertical distance which a falling particle
  ! travels during one time step. As the vertical coordinate is the pressure,
  ! this distance is a pressure difference.
  !
  !-----------------
  ! INPUT VARIABLES:
  !-----------------
  ! TimeStep = duration of time step / s
  ! PRESS    = air pressure within grid box / Pa
  ! TEMP     = air temperature / K
  ! vel      = velocity of fall of aerosol particles / (m/s)
  !
  !--------
  ! OUTPUT:
  !--------
  ! mz_sedi_SedStep  = sedimentation distance within one TimeStep / Pa
  !------------------------------------------------------------------------
  USE messy_main_constants_mem,   ONLY: Rgas => R_gas, g_acc=> g
  
  IMPLICIT NONE

  REAL(dp), INTENT(in) :: TimeStep, PRESS, TEMP, vel
  REAL(dp) :: mz_sedi_SedStep
  !molar mass of dry air / (kg/mol)
  REAL(dp), PARAMETER :: MolMassAir=0.02897_dp

  mz_sedi_SedStep=(g_acc*MolMassAir*PRESS*vel*TimeStep)/(Rgas*TEMP)
END FUNCTION mz_sedi_SedStep


!=========================================================================
!=========================================================================

PURE FUNCTION sedimentation(kproma,klev,           &
                         pTop,                     &
                         pBottom,                  &
                         SedStep,                  &
                         InFrac,                   &
                         val_sed)
  !-------------
  ! DESCRIPTION:
  !-------------
  ! This function calculates the changes in amount-of-substance fractions
  ! (also called "volume mixing ratios") of HNO3 molecules in air due to
  ! sedimentation. It works as well for mass fractions ("mass mixing ratios")
  ! of ice molecules in air.
  ! The idea of the sedimentation algorithm is to calculate the amount of
  ! sedimenting particles not by integrating the vertical amount-of-substance
  ! fraction profile/mass fraction profile but local straight line
  ! approximations to that step function.
  ! No sedimentation takes place from the highest grid box and from the lowest
  ! grid box. 
  !
  !-----------------
  ! INPUT VARIABLES:
  !-----------------
  ! klev          = number of grid boxes within one column
  ! pTop(:)       = air pressure at the top of grid box / Pa
  !   Note: The air pressure at the top of grid box i equals the air pressure
  !   at the bottom of the grid box above, which is grid box (i-1).
  ! pBottom(:)    = air pressure at the bottom of grid box / Pa
  !   Note: The air pressure at the bottom of grid box i equals the air
  !   pressure at the top of the grid box below, which is grid box (i+1).
  ! SedStep(:)    = sedimentation distance within one TimeStep / Pa
  ! InFrac(:)     = mass fraction of ice in air / (kg/kg)
  !              or amount-of-substance fraction of HNO3 in air / (mol/mol)
  ! val_sed(:)    = flag indicating grid boxes where sedimentation takes place
  !
  !------------------
  ! OUTPUT VALUES:
  !------------------
  ! mz_sedi_sed(:) = change of amount-of-substance fraction of [substance] 
  !                  in air / (mol/mol)
  !------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) :: kproma,klev
  REAL(dp), DIMENSION(kproma,klev), INTENT(in) :: &
    pTop, pBottom, SedStep, InFrac
  LOGICAL, DIMENSION(kproma,klev), INTENT(in) :: val_sed

  INTEGER :: i,j
  REAL(dp), DIMENSION(kproma,klev) :: frac, change, sedimentation
  LOGICAL, DIMENSION(kproma,klev) :: val_loc

  INTRINSIC MAX

  frac=max(InFrac, 0.0_dp)
  change=0.0_dp

  val_loc = val_sed

  SELECT CASE(scheme)

  CASE(1,2,3,4)

     val_loc(1:kproma,1) = .false.
     val_loc(1:kproma,klev) = .false.

  CASE(5,6)

     val_loc(1:kproma,1:3) = .false.
     val_loc(1:kproma,klev) = .false.

  END SELECT

  SELECT CASE(scheme)

  CASE(1)

     FORALL (i=2:klev-1,j=1:kproma, val_loc(j,i) .AND. (.NOT. val_loc(j,i-1)))
        change(j,i)=-Trapezoid(frac(j,i-1), pTop(j,i-1),               &
             frac(j,i), pTop(j,i), pBottom(j,i),       &
             frac(j,i+1), pBottom(j,i+1),            &
             SedStep(j,i))                         &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=3:klev-1,j=1:kproma, val_loc(j,i) .AND. val_loc(j,i-1))
        change(j,i)= Trapezoid(frac(j,i-2), pTop(j,i-2),               &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1), &
             frac(j,i), pBottom(j,i),                &
             SedStep(j,i-1))                       &
             /(pBottom(j,i)-pTop(j,i))                         &
             -Trapezoid(frac(j,i-1), pTop(j,i-1),               &
             frac(j,i), pTop(j,i), pBottom(j,i),       &
             frac(j,i+1), pBottom(j,i+1),            &
             SedStep(j,i))                         &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=3:klev,j=1:kproma, (.NOT. val_loc(j,i)) .AND. val_loc(j,i-1))
        change(j,i)= Trapezoid(frac(j,i-2), pTop(j,i-2),               &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1), &
             frac(j,i), pBottom(j,i),                &
             SedStep(j,i-1))                       &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

  CASE (2)
     
     FORALL (i=2:klev-1, j=1:kproma, val_loc(j,i) .AND. (.NOT. val_loc(j,i-1)))
        change(j,i)=-STAIRS_LIN(frac(j,i-1), pTop(j,i-1),        &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=3:klev-1, j=1:kproma, val_loc(j,i) .AND. val_loc(j,i-1))
        change(j,i)= STAIRS_LIN(frac(j,i-2), pTop(j,i-2),        &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))                           &
             -STAIRS_LIN(frac(j,i-1), pTop(j,i-1),               &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=3:klev, j=1:kproma, (.NOT. val_loc(j,i)) .AND. val_loc(j,i-1))
        change(j,i)= STAIRS_LIN(frac(j,i-2), pTop(j,i-2),        &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

  CASE (3)

     FORALL (i=2:klev-1, j=1:kproma, val_loc(j,i) .AND. (.NOT. val_loc(j,i-1)))
        change(j,i)=-STAIRS_LINSND(frac(j,i-1), pTop(j,i-1),     &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=3:klev-1, j=1:kproma, val_loc(j,i) .AND. val_loc(j,i-1))
        change(j,i)= STAIRS_LINSND(frac(j,i-2), pTop(j,i-2),     &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))                           &
             -STAIRS_LINSND(frac(j,i-1), pTop(j,i-1),            &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=3:klev, j=1:kproma, (.NOT. val_loc(j,i)) .AND. val_loc(j,i-1))
        change(j,i)= STAIRS_LINSND(frac(j,i-2), pTop(j,i-2),     &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

  CASE (4)

     FORALL (i=2:klev-1, j=1:kproma, val_loc(j,i) .AND. (.NOT. val_loc(j,i-1)))
        change(j,i)=-STAIRS_PSEUDO(frac(j,i-1), pTop(j,i-1),     &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=3:klev-1, j=1:kproma, val_loc(j,i) .AND. val_loc(j,i-1))
        change(j,i)= STAIRS_PSEUDO(frac(j,i-2), pTop(j,i-2),     &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))                           &
             -STAIRS_PSEUDO(frac(j,i-1), pTop(j,i-1),            &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=3:klev, j=1:kproma, (.NOT. val_loc(j,i)) .AND. val_loc(j,i-1))
        change(j,i)= STAIRS_PSEUDO(frac(j,i-2), pTop(j,i-2),     &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

  CASE(5)

     FORALL (i=4:klev-1, j=1:kproma, val_loc(j,i) .AND. (.NOT. val_loc(j,i-1)))
        change(j,i)= -STAIRS_PSEUDOMAX(frac(j,i-1), pTop(j,i-1), &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             frac(j,i-2),pTop(j,i-2),                            &
             frac(j,i-3),pTop(j,i-3),                            &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=5:klev-1, j=1:kproma, val_loc(j,i) .AND. val_loc(j,i-1))
        change(j,i)= STAIRS_PSEUDOMAX(frac(j,i-2), pTop(j,i-2),  &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             frac(j,i-3),pTop(j,i-3),                            &
             frac(j,i-4),pTop(j,i-4),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))                           &
             -STAIRS_PSEUDOMAX(frac(j,i-1), pTop(j,i-1),         &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             frac(j,i-2),pTop(j,i-2),                            &
             frac(j,i-3),pTop(j,i-3),                            &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=5:klev, j=1:kproma, (.NOT. val_loc(j,i)) .AND. val_loc(j,i-1))
        change(j,i)= STAIRS_PSEUDOMAX(frac(j,i-2), pTop(j,i-2),  &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             frac(j,i-3),pTop(j,i-3),                            &
             frac(j,i-4),pTop(j,i-4),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

  CASE(6)

     FORALL (i=4:klev-1, j=1:kproma, val_loc(j,i) .AND. (.NOT. val_loc(j,i-1)))
        change(j,i)= -WALCEK(frac(j,i-1), pTop(j,i-1),           &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             frac(j,i-2),pTop(j,i-2),                            &
             frac(j,i-3),pTop(j,i-3),                            &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=5:klev-1, j=1:kproma, val_loc(j,i) .AND. val_loc(j,i-1))
        change(j,i)= WALCEK(frac(j,i-2), pTop(j,i-2),            &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             frac(j,i-3),pTop(j,i-3),                            &
             frac(j,i-4),pTop(j,i-4),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))                           &
             -WALCEK(frac(j,i-1), pTop(j,i-1),                   &
             frac(j,i), pTop(j,i), pBottom(j,i),                 &
             frac(j,i+1), pBottom(j,i+1),                        &
             frac(j,i-2),pTop(j,i-2),                            &
             frac(j,i-3),pTop(j,i-3),                            &
             SedStep(j,i),i)                                     &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

     FORALL (i=5:klev, j=1:kproma, (.NOT. val_loc(j,i)) .AND. val_loc(j,i-1))
        change(j,i)= WALCEK(frac(j,i-2), pTop(j,i-2),            &
             frac(j,i-1), pTop(j,i-1), pBottom(j,i-1),           &
             frac(j,i), pBottom(j,i),                            &
             frac(j,i-3),pTop(j,i-3),                            &
             frac(j,i-4),pTop(j,i-4),                            &
             SedStep(j,i-1),i-1)                                 &
             /(pBottom(j,i)-pTop(j,i))
     END FORALL

  END SELECT

  sedimentation=MAX(change,-InFrac)

CONTAINS

!--------------------------------------------------------------------------
ELEMENTAL FUNCTION Trapezoid(frac_above, pTop_above,    &
                             frac_i, pTop_i, pBottom_i, &
                             frac_below, pBottom_below, &
                             SedStep)
  !------------------------------------------------------------------------
  ! DESCRIPTION:
  !-------------
  ! The simplest algorithm to calculate the changes in the amount-of-
  ! substance fractions of ice particles in air due to sedimentation
  ! ("simple upwind scheme") includes in its equations
  ! the product of the amount-of-substance fraction of ice particles in air
  ! and the height of the vertical air layer from which ice particles fall 
  ! into the next box.
  ! This simple algorithm has the disadvantage of "numerical diffusion".
  ! The problem of numerical diffusion can be reduced if the product of
  ! the amount-of-substance fraction and the sedimentation step
  ! is replaced by a more sophisticated calculation. The function Trapezoid
  ! uses integrals over straight line approximations of the amount-of-
  ! substance fraction distribution for that purpose. Definite integrals
  ! over straight lines have trapezoidal shape, hence the function name.
  ! Where SedStep is greater than half the box height, the calculation
  ! becomes more difficult; the function does not return just the area
  ! of a trapezoid but rather the sum of the areas of a trapezoid and a
  ! rectangle.
  !------------------------------------------------------------------------
  IMPLICIT NONE

  REAL(dp), INTENT(IN) :: &
    frac_above, pTop_above,    &
    frac_i, pTop_i, pBottom_i, &
    frac_below, pBottom_below, &
    SedStep
  REAL(dp) :: &
    Trapezoid, Trapezoid_tmp, &
    SlopeAbove, InterceptAbove, SlopeBelow, InterceptBelow, &
    pHeight_i
    
  INTRINSIC min

  pHeight_i=pBottom_i-pTop_i
  Trapezoid=0.0_dp

  IF (frac_above<=frac_i .AND. frac_i<=frac_below &
      .AND. SedStep<=0.5_dp*pHeight_i) THEN
    !----------------------------------------------------------------------
    ! If the fraction value frac increases from box i-1 to box i and
    ! from box i to box i+1, the straight line approximation method
    ! increases the amount of transported particles compared to
    ! simple volume shifting.
    !----------------------------------------------------------------------
    SlopeBelow=(frac_i-frac_below)/(0.5_dp*(pTop_i-pBottom_below))
    InterceptBelow=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeBelow
    Trapezoid=SedStep                        &
              *(SlopeBelow*pBottom_i         &
                -0.5_dp*SlopeBelow*SedStep   &
                +InterceptBelow)
    IF (Trapezoid>pHeight_i*frac_i) Trapezoid=pHeight_i*frac_i
  ELSEIF (frac_above<=frac_i .AND. frac_i<=frac_below &
          .AND. SedStep>0.5_dp*pHeight_i) THEN
    !----------------------------------------------------------------------
    ! If the fraction value frac increases from box i-1 to box i and
    ! from box i to box i+1, the straight line approximation method
    ! increases the amount of transported particles compared to
    ! simple volume shifting.
    ! However, where the straight line values are lower than the average
    ! frac-values in box i, the average frac-values are integrated instead
    ! of the line.
    !----------------------------------------------------------------------
    SlopeBelow=(frac_i-frac_below)/(0.5_dp*(pTop_i-pBottom_below))
    InterceptBelow=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeBelow
    Trapezoid= 0.5_dp*pHeight_i                       &
               *(SlopeBelow*pBottom_i                 &
                 -0.25_dp*SlopeBelow*pHeight_i        &
                 +InterceptBelow)                     &
              +frac_i*(SedStep-0.5_dp*pHeight_i)
    IF (Trapezoid>pHeight_i*frac_i) Trapezoid=pHeight_i*frac_i
  ELSEIF (frac_above>frac_i .AND. frac_i<=frac_below) THEN
    !---------------------------------------------------------------------
    ! If box i is a local minimum, sedimentation is calculated like in
    ! the simple upwind scheme
    !---------------------------------------------------------------------
    Trapezoid=frac_i*MIN(SedStep,pHeight_i)
  ELSEIF (frac_above<=frac_i .AND. frac_i>frac_below) THEN
    !---------------------------------------------------------------------
    ! If box i is a local maximum, sedimentation is calculated like in
    ! the simple upwind scheme
    !---------------------------------------------------------------------
    Trapezoid=frac_i*MIN(SedStep,pHeight_i)
  ELSEIF (frac_above>frac_i .AND. frac_i>frac_below &
          .AND. SedStep<=0.5_dp*pHeight_i) THEN
    !---------------------------------------------------------------------
    ! If the fraction value frac decreases between the boxes i-1 and i
    ! as well as between the boxes i and i+1, numerical diffusion
    ! correction by integrating straight line approximation can be applied
    ! in two ways. The one that yields the smaller Trapezoid is
    ! choosen.
    !---------------------------------------------------------------------
    SlopeAbove=(frac_above-frac_i)/(0.5_dp*(pTop_above-pBottom_i))
    InterceptAbove=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeAbove
    Trapezoid_tmp=SedStep                        &
                  *(SlopeAbove*pBottom_i         &
                    -0.5_dp*SlopeAbove*SedStep   &
                    +InterceptAbove)

    ! Note that Trapezoid_tmp can become negative, if the straight line
    ! intercepts the frac=0 line within the grid box i. 
    ! In this case the calculations below will result in a zero 
    ! sedimentation. The fix here is to approximate the mixing ratio
    ! distribution within the box by a triangle and to sediment the
    ! lower part (of height SedStep).
    IF (Trapezoid_tmp < 0.0_dp) &
         Trapezoid_tmp = 0.5_dp * SedStep * SedStep*(frac_i/(0.5*pHeight_i))

    SlopeBelow=(frac_i-frac_below)/(0.5_dp*(pTop_i-pBottom_below))
    InterceptBelow=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeBelow
    Trapezoid=SedStep                        &
              *(SlopeBelow*pBottom_i         &
                -0.5_dp*SlopeBelow*SedStep   &
                +InterceptBelow)
    IF (Trapezoid_tmp<Trapezoid) Trapezoid=Trapezoid_tmp
  ELSEIF (frac_above>frac_i .AND. frac_i>frac_below &
          .AND. SedStep>0.5_dp*pHeight_i) THEN
    !---------------------------------------------------------------------
    ! If the fraction value frac decreases between the boxes i-1 and i
    ! as well as between the boxes i and i+1, numerical diffusion
    ! correction by integrating straight line approximation can be applied
    ! in two ways. The one that yields the smaller Trapezoid is
    ! choosen.
    ! Again, where the straigt line values are higher than the
    ! average concentration values in box i, the average concentration
    ! values are integrated instead of the line.
    !---------------------------------------------------------------------
    SlopeAbove=(frac_above-frac_i)/(0.5_dp*(pTop_above-pBottom_i))
    InterceptAbove=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeAbove
    Trapezoid_tmp= 0.5_dp*pHeight_i                  &
                   *(SlopeAbove*pBottom_i            &
                     -0.25_dp*SlopeAbove*pHeight_i   &
                     +InterceptAbove)                &
                  +frac_i*(SedStep-0.5_dp*pHeight_i)

    ! Note that Trapezoid_tmp can become negative, if the straight line
    ! intercepts the frac=0 line within the grid box i. 
    ! In this case the calculations below will result in a zero 
    ! sedimentation. The fix here is to approximate the mixing ratio
    ! distribution within the box by a triangle and to sediment the
    ! lower part (of height pHeight_i/2) plus the rectangle above
    ! (of height SedStep-pHeight_i).
    IF (Trapezoid_tmp < 0.0_dp) &
         Trapezoid_tmp = (SedStep - 0.25_dp*pHeight_i) * frac_i
    SlopeBelow=(frac_i-frac_below)/(0.5_dp*(pTop_i-pBottom_below))
    InterceptBelow=frac_i-0.5_dp*(pBottom_i+pTop_i)*SlopeBelow

    Trapezoid= 0.5_dp*pHeight_i                  &
               *(SlopeBelow*pBottom_i            &
                 -0.25_dp*SlopeBelow*pHeight_i   &
                 +InterceptBelow)                &
              +frac_i*(SedStep-0.5_dp*pHeight_i)
    IF (Trapezoid_tmp<Trapezoid) Trapezoid=Trapezoid_tmp
  END IF
  IF (Trapezoid<0.0_dp) Trapezoid=0.0_dp
END FUNCTION Trapezoid
!--------------------------------------------------------------------------

ELEMENTAL FUNCTION STAIRS_LIN(frac_ap, ptop_a,            &
                               frac_i, ptop_i, pbottom_i, &
                               frac_bp, pbottom_b,        &
                               sedstep,i)
!*********************************************************
!DESCRIPTION
!The scheme subdivides the gridbox into two areas. The
!mixing ratio in each is estimated using a linear approach
!to the mixing ratio at the interface as given by the
!average mixing ratio within the grid box and the
!respective mixing ratio in the contiguous grid box. The
!fraction of the grid box taken by each area is calculated
!such that their average value equals the observed average
!*********************************************************
!
  IMPLICIT NONE
!
  REAL(dp), INTENT(IN) :: frac_ap, frac_i, frac_bp, sedstep
  REAL(dp), INTENT(IN) :: ptop_a, ptop_i, pbottom_i, pbottom_b
  INTEGER, INTENT(IN) :: i
  REAL(dp) :: STAIRS_LIN
!
  REAL(dp) :: slope_a, slope_b
  REAL(dp) :: pmid_a, pmid_i, pmid_b, deltap
  REAL(dp) :: frac_ia, frac_ib, xia, xib
  REAL(dp) :: value
  REAL(dP) :: frac_a,frac_b
!
  pmid_a=0.5_dp*(ptop_a+ptop_i)
  pmid_i=0.5_dp*(ptop_i+pbottom_i)
  pmid_b=0.5_dp*(pbottom_i+pbottom_b)
  deltap=pbottom_i-ptop_i
!
  frac_a=frac_ap*pmid_a/pmid_i
  frac_b=frac_bp*pmid_b/pmid_i
!
  IF (((frac_a>=frac_i).and.(frac_i>=frac_b)).OR. &
      ((frac_a<=frac_i).and.(frac_i<=frac_b))) THEN
!
   slope_a=(frac_i-frac_a)/(pmid_i-pmid_a)
   slope_b=(frac_b-frac_i)/(pmid_b-pmid_i)
!
   frac_ia=frac_a+slope_a*(ptop_i-pmid_a)
   frac_ib=frac_i+slope_b*(pbottom_i-pmid_i)
!
   if (frac_ia/=frac_ib) then
    xia=(frac_i-frac_ib)/(frac_ia-frac_ib)
   else
    xia=1._dp
   endif
   xib=1._dp-xia
!
   if (xib*deltap>=sedstep) then
    value=frac_ib*sedstep
   elseif (deltap>sedstep) then
    value=frac_ib*xib*deltap+frac_ia*(sedstep-xib*deltap)
   else
    value=frac_i*deltap
   endif
!
  ELSE
!
   if (deltap>sedstep) then
    value=frac_i*sedstep
   else
    value=frac_i*deltap
   endif
!
  ENDIF
!
  STAIRS_LIN=value
!
END FUNCTION STAIRS_LIN

ELEMENTAL FUNCTION STAIRS_LINSND(frac_ap, ptop_a,         &
                               frac_i, ptop_i, pbottom_i, &
                               frac_bp, pbottom_b,        &
                               sedstep,i)
!**********************************************************
!DESCRIPTION
!The formalism of this scheme is similar to the one of
!STAIRS_LIN except that a second order correction to the
!mixing ratio estimation at the interface is applied.
!The estimated mixing ratios are constrained such that
!frac_ij part of interval [frac_i,frac_j].
!**********************************************************
!
  IMPLICIT NONE
!
  REAL(dp), INTENT(IN) :: frac_ap, frac_i, frac_bp, sedstep
  REAL(dp), INTENT(IN) :: ptop_a, ptop_i, pbottom_i, pbottom_b
  INTEGER, INTENT(IN) :: i
  REAL(dp) :: STAIRS_LINSND
!
  REAL(dp) :: slope_a, slope_b, slope_bma
  REAL(dp) :: pmid_a, pmid_i, pmid_b, deltap
  REAL(dp) :: frac_ia, frac_ib, xia, xib
  REAL(dp) :: value
  REAL(dP) :: frac_a,frac_b
!
  pmid_a=0.5_dp*(ptop_a+ptop_i)
  pmid_i=0.5_dp*(ptop_i+pbottom_i)
  pmid_b=0.5_dp*(pbottom_i+pbottom_b)
  deltap=pbottom_i-ptop_i
!
  frac_a=frac_ap*pmid_a/pmid_i
  frac_b=frac_bp*pmid_b/pmid_i
!
  IF (((frac_a>=frac_i).and.(frac_i>=frac_b)).OR. &
      ((frac_a<=frac_i).and.(frac_i<=frac_b))) THEN
!
   slope_a=(frac_i-frac_a)/(pmid_i-pmid_a)
   slope_b=(frac_b-frac_i)/(pmid_b-pmid_i)
!
   slope_bma=(slope_b-slope_a)/(pmid_b-pmid_a)
!
   frac_ia=frac_a+slope_a*(ptop_i-pmid_a) &
                 +0.5_dp*slope_bma*(ptop_i-pmid_a)**2._dp
   frac_ib=frac_b+slope_b*(pbottom_i-pmid_b) &
                 +0.5_dp*slope_bma*(pbottom_i-pmid_a)**2._dp
!
   IF ((frac_a>=frac_i).and.(frac_i>=frac_b)) THEN
    if (frac_ia<frac_i) then
     frac_ia=frac_i
    elseif (frac_ia>frac_a) then
     frac_ia=frac_a
    endif
    if (frac_ib>frac_i) then
     frac_ib=frac_i
    elseif (frac_ib<frac_b) then
     frac_ib=frac_b
    endif
   ELSE
    if (frac_ia>frac_i) then
     frac_ia=frac_i
    elseif (frac_ia<frac_a) then
     frac_ia=frac_a
    endif
    if (frac_ib<frac_i) then
     frac_ib=frac_i
    elseif (frac_ib>frac_b) then
     frac_ib=frac_b
    endif
   ENDIF
!
   if (frac_ia/=frac_ib) then
    xia=(frac_i-frac_ib)/(frac_ia-frac_ib)
   else
    xia=1._dp
   endif
   xib=1._dp-xia
!
   if (xib*deltap>=sedstep) then
    value=frac_ib*sedstep
   elseif (deltap>sedstep) then
    value=frac_ib*xib*deltap+frac_ia*(sedstep-xib*deltap)
   else
    value=frac_i*deltap
   endif
!
  ELSE
!
   if (deltap>sedstep) then
    value=frac_i*sedstep
   else
    value=frac_i*deltap
   endif
!
  ENDIF
!
  STAIRS_LINSND=value
!
END FUNCTION STAIRS_LINSND

ELEMENTAL FUNCTION STAIRS_PSEUDO(frac_ap, ptop_a,         &
                               frac_i, ptop_i, pbottom_i, &
                               frac_bp, pbottom_b,        &
                               sedstep,i)
!**********************************************************
!DESCRIPTION
!The scheme subdivides the gridbox into two areas. The
!mixing ratio in each is estimated using a semi
!pseudo-equilibrium approach as given by the mixing ratio
!of the respective contiguous gridbox. The fraction of the
!grid box taken by each area is calculated such that their
!average value equals the observed average. The pseudo-
!equilibrium approach implies that the scheme underestimates
!transport during spin-up, and therefore requires a longer
!spin-up period.
!**********************************************************
!
  IMPLICIT NONE
!
  REAL(dp), INTENT(IN) :: frac_ap, frac_i, frac_bp, sedstep
  REAL(dp), INTENT(IN) :: ptop_a, ptop_i, pbottom_i, pbottom_b
  INTEGER, INTENT(IN) :: i
  REAL(dp) :: STAIRS_PSEUDO
!
  REAL(dp) :: pmid_a, pmid_i, pmid_b, deltap
  REAL(dp) :: frac_ia, frac_ib, xia, xib
  REAL(dp) :: value
  REAL(dP) :: frac_a,frac_b
!
  pmid_a=0.5_dp*(ptop_a+ptop_i)
  pmid_i=0.5_dp*(ptop_i+pbottom_i)
  pmid_b=0.5_dp*(pbottom_i+pbottom_b)
  deltap=pbottom_i-ptop_i
!
  frac_a=frac_ap*pmid_a/pmid_i
  frac_b=frac_bp*pmid_b/pmid_i
!
  IF (((frac_a>=frac_i).and.(frac_i>=frac_b)).OR. &
      ((frac_a<=frac_i).and.(frac_i<=frac_b))) THEN
!
   frac_ia=frac_a
   frac_ib=frac_b
!
   if (frac_ia/=frac_ib) then
    xia=(frac_i-frac_ib)/(frac_ia-frac_ib)
   else
    xia=1._dp
   endif
   xib=1._dp-xia
!
   if (xib*deltap>=sedstep) then
    value=frac_ib*sedstep
   elseif (deltap>sedstep) then
    value=frac_ib*xib*deltap+frac_ia*(sedstep-xib*deltap)
   else
    value=frac_i*deltap
   endif
!
  ELSE
!
   if (deltap>sedstep) then
    value=frac_i*sedstep
   else
    value=frac_i*deltap
   endif
!
  ENDIF
!
  STAIRS_PSEUDO=value
!
END FUNCTION STAIRS_PSEUDO

ELEMENTAL FUNCTION STAIRS_PSEUDOMAX(frac_ap, ptop_a,      &
                               frac_i, ptop_i, pbottom_i, &
                               frac_bp, pbottom_b,        &
                               frac_aap, ptop_aa,         &
                               frac_aaap, ptop_aaa,       &
                               sedstep, i)
!**********************************************************
!DESCRIPTION
!The mechanism of the scheme is similar to the one of
!STAIRS_PSEUDO except that it aims to constrains plumes to
!a single grid box via additional plume maximum
!considerations. It should only be used in conjunction with
!instantaneously pulsed injections schemes unless it may
!substantiously underestimate transport.
!**********************************************************
!
  IMPLICIT NONE
!
  REAL(dp), INTENT(IN) :: frac_ap, frac_i, frac_bp, sedstep
  REAL(dp), INTENT(IN) :: ptop_a, ptop_i, pbottom_i, pbottom_b
  REAL(dp), INTENT(IN) :: frac_aaap,frac_aap,ptop_aaa,ptop_aa
  INTEGER, INTENT(IN) :: i
  REAL(dp) :: STAIRS_PSEUDOMAX
!
  REAL(dp) :: slope_a, slope_b, slope_aaa, slope_aa
  REAL(dp) :: pmid_a, pmid_i, pmid_b, pmid_aaa, pmid_aa, deltap
  REAL(dp) :: frac_ia, frac_ib, xia, xib
  REAL(dp) :: frac_aaa,frac_aa,frac_a,frac_b
  REAL(dp) :: value
  INTEGER :: sedcase !1=minimum or maximum, 2=below maximum, 3=neither 1 nor 2
!
  pmid_aaa=0.5_dp*(ptop_aaa+ptop_aa)
  pmid_aa=0.5_dp*(ptop_aa+ptop_a)
  pmid_a=0.5_dp*(ptop_a+ptop_i)
  pmid_i=0.5_dp*(ptop_i+pbottom_i)
  pmid_b=0.5_dp*(pbottom_i+pbottom_b)
  deltap=pbottom_i-ptop_i
!
  frac_aaa=frac_aaap*pmid_aaa/pmid_i
  frac_aa=frac_aap*pmid_aa/pmid_i
  frac_a=frac_ap*pmid_a/pmid_i
  frac_b=frac_bp*pmid_b/pmid_i
!
  slope_aaa=(frac_aa-frac_aaa)/(pmid_aa-pmid_aaa)
  slope_aa=(frac_a-frac_aa)/(pmid_a-pmid_aa)
  slope_a=(frac_i-frac_a)/(pmid_i-pmid_a)
  slope_b=(frac_b-frac_i)/(pmid_b-pmid_i)
!
  IF ((frac_a>frac_i).and.(frac_i>=frac_b)) THEN
   if ((frac_aa<=frac_a).and.(slope_aaa<slope_aa)) then
    sedcase=2 ! situated below local maximum
   else
    sedcase=3 ! monotonic decrease
   endif
  ELSEIF ((frac_a<=frac_i).and.(frac_i<=frac_b)) THEN
   if (slope_a>slope_b) then
    sedcase=1 ! dynamically interpreted as maximum
   else
    sedcase=3 ! monotonic increase
   endif
  ELSEIF ((frac_a<=frac_i).and.(frac_i>frac_b)) THEN
   if (slope_aa>slope_a) then
    sedcase=2 ! dynamically interpreted as below local maximum
   else
    sedcase=1 ! local maximum
   endif
  ELSE
   sedcase=1 ! local minimum
  ENDIF
!
  SELECT CASE(sedcase)
!
  CASE(3)
!
   frac_ia=frac_a
   frac_ib=frac_b
!
   if (frac_ia/=frac_ib) then
    xia=(frac_i-frac_ib)/(frac_ia-frac_ib)
   else
    xia=1._dp
   endif
   xib=1._dp-xia
!
   if (xib*deltap>=sedstep) then
    value=frac_ib*sedstep
   elseif (deltap>sedstep) then
    value=frac_ib*xib*deltap+frac_ia*(sedstep-xib*deltap)
   else
    value=frac_i*deltap
   endif
!
  CASE(2)
!
   frac_ib=frac_b
!
   if (deltap>sedstep) then
    value=frac_ib*sedstep
   else
    value=frac_ib*deltap
   endif
!
  CASE(1)
!
   if (deltap>sedstep) then
    value=frac_i*sedstep
   else
    value=frac_i*deltap
   endif
!
  END SELECT
!
  STAIRS_PSEUDOMAX=value
!
END FUNCTION STAIRS_PSEUDOMAX

ELEMENTAL FUNCTION WALCEK(frac_ap, ptop_a,    &
                  frac_i, ptop_i, pbottom_i,  &
                  frac_bp, pbottom_b,         &
                  frac_aap,         &
                  ptop_aa, frac_aaap, ptop_aaa,    &
                  sedstep,i)
!*************************************************************
!DESCRIPTION
!The scheme implements a sedimentation mechanism based on the
!Walcek advection scheme (Walcek, JGR, 2000). However, it shows
!the following differences with the original advection scheme:
!(1) Sedimentation not being monotonic due to vertical final
!drop velocity dissimilarities due to particle size among
!other, the original approach to ensure monotonic behaviour
!is not implemented.
!(2) Following (1) it cannot be excluded that the amount of
!transported particles is negative. For this reason the
!mixing ratio at the lower interface is set not to be inferior
!to the one in the grid box below accordingly to the
!STAIRS_PSEUDO scheme.
!*************************************************************
!
  IMPLICIT NONE
!
  REAL(dp), INTENT(IN) :: frac_ap, frac_i, frac_bp, sedstep
  REAL(dp), INTENT(IN) :: ptop_a, ptop_i, pbottom_i, pbottom_b
  REAL(dp), INTENT(IN) :: frac_aaap,frac_aap,ptop_aaa,ptop_aa
  INTEGER, INTENT(IN) :: i
  REAL(dp) :: WALCEK
!
  REAL(dp) :: slope_ab, courant
  REAL(dp) :: pmid_a, pmid_b, pmid_i, pmid_aaa, pmid_aa, deltap
  REAL(dp) :: frac_a,frac_b,frac_aa,frac_aaa
  REAL(dp) :: frac_ia, frac_ib, frac_step
  REAL(dp) :: xia, xib
  REAL(dp) :: value, value_min
  INTEGER :: sedcase ! =1 local max or min , =2 downwind 1 cell local max or min
                     ! =3 downwind 2 cells local max or min, =4 neither 2 nor 3
                     ! nor 1
!
   pmid_a=0.5_dp*(ptop_a+ptop_i)
   pmid_b=0.5_dp*(pbottom_i+pbottom_b)
   pmid_i=0.5_dp*(ptop_i+pbottom_i)
   pmid_aa=0.5_dp*(ptop_aa+ptop_a)
   pmid_aaa=0.5_dp*(ptop_aaa+ptop_aa)
!
   frac_a=frac_ap*pmid_a/pmid_i
   frac_b=frac_bp*pmid_b/pmid_i
   frac_aa=frac_aap*pmid_aa/pmid_i
   frac_aaa=frac_aaap*pmid_aaa/pmid_i
!
  IF (((frac_i>frac_b).AND.(frac_i>=frac_a)).OR. &
      ((frac_i<frac_b).AND.(frac_i<=frac_a))) THEN
   sedcase=1 ! local max or min
  ELSE
   if (((frac_i>frac_a).AND.(frac_aa>=frac_a)).OR. &
       ((frac_i<frac_a).AND.(frac_aa<=frac_a))) then
    sedcase=2 ! down 1 box of local max or min
   elseif (((frac_a<frac_i).AND.(frac_aa<frac_a).AND.(frac_aaa>=frac_aa)).OR. &
           ((frac_a>frac_i).AND.(frac_aa>frac_a).AND.(frac_aaa<=frac_aa))) then
    sedcase=3 ! down 2 boxes of local max or min
   else
    sedcase=4 ! simple monotonic decrease or increase
   endif
  ENDIF
!
  deltap=pbottom_i-ptop_i
!
  IF (sedcase>1) THEN
   slope_ab=(frac_b-frac_a)/(pmid_b-pmid_a)
!
   frac_ib=frac_i+slope_ab*0.5_dp*deltap
   frac_step=frac_ib-slope_ab*sedstep
  ENDIF
!
!
  SELECT CASE (sedcase)
!
  CASE(1)
!
   if (deltap>sedstep) then
    value=frac_i*sedstep
   else
    value=frac_i*deltap
   endif
!
  CASE(4)
!
   value=0.5_dp*(frac_ib+frac_step)*sedstep
!
  CASE(2)
!
   courant=1.75_dp-0.45_dp*sedstep/deltap
!
   value=(frac_i+0.5_dp*(frac_ib+frac_step-2._dp*frac_i)*courant)*sedstep
!
  CASE(3)
!
   courant=1.2_dp+0.6_dp*sedstep/deltap
   courant=MAX(1.5_dp,courant)
!
   value=(frac_i+0.5_dp*(frac_ib+frac_step-2._dp*frac_i)*courant)*sedstep
!
  END SELECT
!
  IF (sedcase>1) THEN
   frac_ia=frac_a
   frac_ib=frac_b
!
   if (frac_ia/=frac_ib) then
    xia=(frac_i-frac_ib)/(frac_ia-frac_ib)
   else
    xia=1._dp
   endif
   xib=1._dp-xia
!
   if (xib*deltap>=sedstep) then
    value_min=frac_ib*sedstep
   elseif (deltap>sedstep) then
    value_min=frac_ib*xib*deltap+frac_ia*(sedstep-xib*deltap)
   else
    value_min=frac_i*deltap
   endif
!
   if (value<value_min) then
    value=value_min
   endif
  ENDIF
!
  WALCEK=value
!
END FUNCTION WALCEK

END FUNCTION sedimentation

!-----------------------------------------------------------------------------
!****************************************************************************

END MODULE messy_sedi
