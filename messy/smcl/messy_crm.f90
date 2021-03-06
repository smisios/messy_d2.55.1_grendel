MODULE MESSY_CRM
  
  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: dp
  USE MESSY_MAIN_TOOLS,         ONLY: t_reset_par

  IMPLICIT NONE

  SAVE
  PRIVATE

  CHARACTER(LEN=*), PUBLIC, PARAMETER :: MODSTR='crm'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: MODVER='1.0.1'

  REAL(dp), PUBLIC :: crm_timestep, crm_size, crm_top

  INTEGER,  PUBLIC :: ngridcells_x, ngridcells_y
  INTEGER,  PUBLIC :: crm_orient, CRM_3D
  INTEGER,  PUBLIC :: micro_scheme
  INTEGER,  PUBLIC :: crmvars
  INTEGER,  PUBLIC :: nmicro_fields
  INTEGER,  PUBLIC :: dosgtracer ! tracer transport procedure (0, 1 or 2)

  LOGICAL,  PUBLIC :: docloud, doprecip
  LOGICAL,  PUBLIC :: dosmagor, dosgs
  LOGICAL,  PUBLIC :: dosurface, dosfc_flx_fxd, dosfc_tau_fxd
  LOGICAL,  PUBLIC :: dowallx, dowally
  LOGICAL,  PUBLIC :: docoriolis, docolumn

  LOGICAL,  PUBLIC :: dodamping
  REAL(dp), PUBLIC :: set_damp(3) = (/60._dp, 450._dp, 0.4_dp/)

  ! parameter settings for one-mom microphyiscs
  REAL(dp), PUBLIC :: qcw0         = 1.e-3_dp
  REAL(dp), PUBLIC :: qci0         = 1.e-4_dp
  REAL(dp), PUBLIC :: alphaelq     = 1.e-3_dp
  REAL(dp), PUBLIC :: betaelq      = 1.e-3_dp
  REAL(dp), PUBLIC :: qp_threshold = 1.e-8_dp

  ! parameter settings for two-mom microphyiscs
  LOGICAL,  PUBLIC :: doicemicro, dograupel, dohail, dosb_warm_rain
  LOGICAL,  PUBLIC :: dopredictNc, dospecifyaerosol, dosubgridw
  LOGICAL,  PUBLIC :: doarcticicenucl, docloudedgeact

  TYPE t_crm_work
     TYPE(t_crm_sg), POINTER, DIMENSION(:,:) :: sg  => NULL()
  END TYPE t_crm_work

  TYPE t_crm_sg
     REAL(dp),POINTER                     :: ptr0 => NULL()
     REAL(dp),POINTER, DIMENSION(:,:)     :: ptr2 => NULL()
     REAL(dp),POINTER, DIMENSION(:,:,:)   :: ptr3 => NULL()
     REAL(dp),POINTER, DIMENSION(:,:,:,:) :: ptr4 => NULL()
  END TYPE t_crm_sg

  PUBLIC :: t_crm_work
  PUBLIC :: t_crm_sg

  PUBLIC :: crm_read_nml_ctrl
  PUBLIC :: calc_cape_crm
  PUBLIC :: cloud_droplet_nc_crm

CONTAINS

!=============================================================================
 SUBROUTINE crm_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    INTRINSIC :: TRIM
    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='crm_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status
    LOGICAL                           :: crm_switch   ! control switch
    LOGICAL                           :: micro_switch ! control switch

    NAMELIST /CTRL/ ngridcells_x, ngridcells_y, crm_top, crm_size,    &
                    crm_timestep, CRM_3D, crm_orient, micro_scheme


    NAMELIST /CTRL_SETPARAM/ docloud, doprecip, dosgs,                &
                             dosmagor, dodamping, set_damp,           &
                             dosurface, dosfc_flx_fxd, dosfc_tau_fxd, &
                             dowallx, dowally, docoriolis, docolumn,  &
                             dosgtracer

    NAMELIST /CTRL_MICRO_SAM1MOM/ qcw0, qci0, alphaelq, betaelq,      &
                                  qp_threshold

    ! preparations for two-moment microphysics
    NAMELIST /CTRL_MICRO_M2005/ doicemicro, dograupel, dohail,              &
                                dopredictNc, dospecifyaerosol, dosubgridw,  &
                                dosb_warm_rain, doarcticicenucl,            &
                                docloudedgeact

    ! INITIALIZE
    status      = 1 ! ERROR

    ! SET CRM CONTROL SWITCH
    crm_switch   = .true.
    micro_switch = .true.

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) THEN
       WRITE(*,*) 'WARNING *** FILE '//TRIM(modstr)//'.nml'//'  NOT FOUND !'
       WRITE(*,*) ' CRM SWITCHED OFF !'
       WRITE(*,*) '******************************************************'
       fstat=1
       RETURN    ! <modstr>.nml does not exist
    END IF
    
    ! READ NAMELIST
    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) THEN
       WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
            TRIM(modstr)//'.nml'//' !'
       WRITE(*,*) '******************************************************'
       RETURN  ! error while reading namelist
    END IF
    
    ! READ CRM PARAMETER NAMELIST
    READ(iou, NML=CTRL_SETPARAM, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_SETPARAM', modstr)
    IF (fstat /= 0) THEN
       WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
            TRIM(modstr)//'.nml'//' !'
       WRITE(*,*) '******************************************************'
       RETURN  ! error while reading namelist
    END IF

    ! READ CRM MICRO_SAM1MOM PARAMETER NAMELIST
    READ(iou, NML=CTRL_MICRO_SAM1MOM, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_MICRO_SAM1MOM', modstr)
    IF (fstat /= 0) THEN
       WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
            TRIM(modstr)//'.nml'//' !'
       WRITE(*,*) '******************************************************'
       RETURN  ! error while reading namelist
    END IF

    ! READ CRM MICRO_M2005 PARAMETER NAMELIST
    READ(iou, NML=CTRL_MICRO_M2005, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_MICRO_M2005', modstr)
    IF (fstat /= 0) THEN
       WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
            TRIM(modstr)//'.nml'//' !'
       WRITE(*,*) '******************************************************'
       RETURN  ! error while reading namelist
    END IF

    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE
    WRITE(*,*) 'SETTINGS CONCERNING GRID AND CORE UTILITIES OF CLOUD RESOLVING MODEL:'
    WRITE(*,*) 'CRM grid size is ',crm_size,' m '
    WRITE(*,*) 'number of CRM cells in x-direction: ',ngridcells_x
    WRITE(*,*) 'number of CRM cells in y-direction: ',ngridcells_y
    WRITE(*,*) 'maximum pressure height of CRM column ',crm_top,' Pa'
    WRITE(*,*) 'length of sub-time step ',crm_timestep,' s '

    WRITE(*,*) 'chosen CRM orientation: '
    SELECT CASE(CRM_3D)
       CASE(0)
          WRITE(*,*) 'CRM has 2 dimensions: '
          SELECT CASE (crm_orient)
          CASE(0)
             WRITE(*,*) 'ORIENTATION: EAST - WEST'
          CASE(1)
             WRITE(*,*) 'ORIENTATION: NORTH - SOUTH'
          END SELECT
       CASE(1)
          WRITE(*,*) 'CRM has 3 dimensions - no orientation needed!'
    END SELECT

    WRITE(*,*) 'chosen  microphysical scheme'
    SELECT CASE(micro_scheme)
       CASE(0)
          WRITE(*,*) '...single-moment scheme: SAM1MOM'
          crmvars       = 7
          nmicro_fields = 2 
       CASE(1)
          ! preparations for two-moment microphysics
          WRITE(*,*) '...double-moment scheme: MORRISON2005 ... not impelemented'
          nmicro_fields = 1 ! start with water vapor
          IF (docloud) THEN
             nmicro_fields = nmicro_fields + 1 ! add cloud water mixing ratio
             IF (dopredictNc) THEN
                nmicro_fields = nmicro_fields + 1 ! add cloud water number concentration
             END IF
          END IF
          IF (doprecip) THEN
             nmicro_fields = nmicro_fields + 2 ! add rain mass and number concentration
          END IF
          IF (doicemicro) THEN
             nmicro_fields = nmicro_fields + 4 ! add snow and cloud ice mass and number concentration
          END IF
          IF (dograupel) THEN
             nmicro_fields = nmicro_fields + 2 ! add graupel mass and number concentration
          END IF
          crmvars = 4+nmicro_fields ! u,v,w,T + micro_fields
          
          ! scheme will not run for the following combinations of parameters:
          IF (doicemicro .AND. .NOT.doprecip) THEN
             WRITE(*,*) 'Morrison 2005 Microphysics does not support both doice and .not.doprecip!'
             micro_switch = .false.
          END IF

          IF (dograupel .AND. .NOT.doicemicro) THEN
             WRITE(*,*) 'doicemicro must be .true. for dograupel to be used!'
             micro_switch = .false.
          END IF
          
          IF (dohail .AND. .NOT.dograupel) THEN
             WRITE(*,*) 'dograupel must be .true. for dohail to be used!'
             micro_switch = .false.
          END IF
    END SELECT

    WRITE(*,*) 'Dimension of CRM buffer including number of microphysical variables: ', crmvars
    
    IF (.not.micro_switch) then
       WRITE(*,*) '******************************************************'
       WRITE(*,*) '!!!!!    WARNING: WRONG MICROPHYSICAL SETTINGS   !!!!!'
       WRITE(*,*) '******************************************************'
       status = 1
       RETURN
    ENDIF

    IF (dosgtracer == 0) THEN
       WRITE(*,*) 'Subgrid tracer transport is switched OFF!!!'   
    ELSEIF (dosgtracer == 1) THEN
       WRITE(*,*) 'Subgrid tracer transport with equal CRM tracer distribution is switched ON!!!'
    ELSEIF (dosgtracer == 2) THEN
       WRITE(*,*) 'Subgrid tracer transport with individual CRM tracer distribution is switched ON!!!'
    ELSE
       WRITE(*,*) 'Wrong subgrid tracer transport procedure!!!'
       crm_switch = .false.
       RETURN  ! error while reading namelist
    END IF

    ! check CRM dimensions
    IF (CRM_3D .eq. 0 .and. ngridcells_x /=1 .and. ngridcells_y /=1) THEN 
       crm_switch  = .false.
       WRITE(*,*) 'X or Y dimension should equal 1 for 2-dimensional CRM setup'
    END IF

    IF (ngridcells_x .le. 1 .or. ngridcells_y .le. 1 .and. CRM_3D .eq. 1) THEN
       crm_switch = .false.
       WRITE(*,*) 'X or Y dimension should not be lower or equal to 1 for three-dimensional CRM setup'
    END IF

    IF (CRM_3D <0 .or. CRM_3D >1) THEN
       crm_switch = .false.
       WRITE(*,*) 'CRM_3D namelist entry should be 0 or 1'
    END IF

    IF (crm_orient <0 .or. crm_orient >1) THEN
       crm_switch = .false.
       WRITE(*,*) 'crm_orient namelist entry should be 0 or 1'
    END IF

    IF (crm_top .le. 0) THEN
       crm_switch = .false.
       WRITE(*,*) 'CRM top pressure must be positive!'
    END IF
    
    IF (.not.crm_switch) then
       WRITE(*,*) '******************************************************'
       WRITE(*,*) '!!!!!       WARNING: WRONG CRM DIMENSIONS        !!!!!'
       WRITE(*,*) '******************************************************'
       status = 1
       RETURN
    ENDIF

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE crm_read_nml_ctrl
! ======================================================================
! mz_ss_20131107+
! subroutine that calculates convective available potential energy (CAPE)
! for a vertical sounding of pressure, temperature and mixing ratio.
!
! Based on program calcsound by Kerry Emanuel and program cape_sound by Dominik Brunner.

SUBROUTINE calc_cape_crm(nl,T,P,R,CAPE_s,CAPE_c,PA_s,PA_c,IDX,xlcl)

  USE messy_main_constants_mem, ONLY: cp_air, cpv, clw, rv, rd, alv, dp
  USE messy_main_tools,         ONLY: tlucua

  IMPLICIT NONE

  INTEGER :: nl

! Input: P pressure (Pa), T temperature (K), R mixing ratio (kg/kg)
  REAL(DP), INTENT(IN) :: P(nl), T(nl), R(nl)

! Output: CAPE_c is maximum CAPE in column and CAPE_s is CAPE from surface parcel, same for PA, 
!         level index of origin of parcel with maximum CAPE IDX
  INTEGER, INTENT(OUT) :: IDX(1)
  REAL(DP), INTENT(OUT) :: CAPE_s, CAPE_c, PA_s, PA_c,xlcl

! local
  REAL(DP), PARAMETER :: cpvmcl = clw-cpv
  REAL(DP), PARAMETER :: eps = rd/rv
  REAL(DP), PARAMETER :: T0 = 273.15

  REAL(DP) :: RS, ALVV, EM, SP, RH, CHI, PLCL
  REAL(DP) :: SUM, RG0, TG0, SLP, TG, RG
  REAL(DP) :: CPW, SPG, ENEW, TVM, PM

  REAL(DP), DIMENSION(nl) :: CAPEP, NAP, PAP
  REAL(DP), DIMENSION(nl) :: P_hpa, TC, EV, ES
  REAL(DP), DIMENSION(nl,nl) :: TLP, TLVP, TVPDIF

  INTEGER :: I, J, K
  INTEGER :: INBP, IMIN, it
 
  CAPE_S = 0.
  CAPE_C = 0.
  PA_s = 0.
  PA_c = 0.
  idx(:) = 0.
  xlcl = 0.
  NAP(:) = 0.
  PAP(:) = 0.
  CAPEP(:) = 0.
!
!  *** Water vapor pressure EV and saturation vapor pressure ES ***
!
  TC = T-T0             ! Celsius
  P_hpa = P/100._dp     ! hPa
  EV = R*P_hpa/(eps+R)  ! vapor pressure


!     ES(:) = (tlucua(NINT(T(:)*1000._dp)))/100._dp	! saturation vapor pressure [hPa]

  ES(:) = 6.112_dp*EXP(17.67*TC(:)/(243.5_dp+TC(:))) 



!
!   ***  Begin outer loop, which cycles through parcel origin levels I ***
 !  
  DO I=nl,2*nl/3,-1	! do calculation only for lowest nl/3 levels

!
!   ***  Define the pseudo-adiabatic entropy SP (conserved variable) ***
!
     RS = eps*ES(I)/(P_hpa(I)-ES(I))
     ALVV = alv-cpvmcl*TC(I)
     EM = MAX(EV(I),1.0E-6_dp)
     SP = cp_air*LOG(T(I)) - rd*LOG(P_hpa(I)-EV(I)) + &
          ALVV*R(I)/T(I) - R(I)*rv*LOG(EM/ES(I))

!
!   ***   Find lifted condensation pressure PLCL [hPa]   ***
!
     RH = R(I)/RS	! relative humidity
     RH = MIN(RH,1.0_dp)
     CHI = T(I)/(1669.0_dp-122.0_dp*RH-T(I))
     IF (RH>0) THEN
        PLCL = P_hpa(I)*RH**CHI
     ENDIF
     IF (I==nl) xlcl = plcl
!
!   ***  Begin updraft loop   ***
!
     SUM = 0.0_dp
     RG0 = R(I)
     TG0 = T(I)

     DO J=I,1,-1 
!
!   ***  Calculate estimates of the rates of change of the entropies  ***
!   ***           with temperature at constant pressure               ***
!  
        RS = eps*ES(J)/(P_hpa(J)-ES(J))	! saturation mixing ratio
        ALVV = alv-cpvmcl*TC(J)
        SLP = (cp_air+RS*clw+ALVV*ALVV*RS/(rv*T(J)*T(J)))/T(J)
!   
!   ***  Calculate lifted parcel temperature below its LCL   ***
!
        IF (P_hpa(J) >= PLCL) THEN
           TLP(I,J) = T(I)*(P_hpa(J)/P_hpa(I))**(rd/cp_air)
           TLVP(I,J) = TLP(I,J)*(1.0_dp+R(I)/eps)/(1.0_dp+R(I))
           TVPDIF(I,J) = TLVP(I,J)-T(J)*(1.0_dp+R(J)/eps)/(1.0_dp+R(J))
        ELSE
!
!   ***  Iteratively calculate lifted parcel temperature and mixing   ***
!   ***    ratios for pseudo-adiabatic ascent     ***
!
           TG = T(J)
           RG = RS
           DO K=1,7
              CPW = 0.0_dp
              IF(J < nl) THEN
                 CPW = SUM+clw*0.5_dp*(RG0+RG)*(LOG(TG)-LOG(TG0))
              ENDIF
              EM = RG*P_hpa(J)/(eps+RG)
              ALVV = alv-cpvmcl*(TG-T0)
              SPG = cp_air*LOG(TG) - rd*LOG(P_hpa(J)-EM) + CPW + ALVV*RG/TG
              TG = TG + (SP-SPG)/SLP
!              ENEW = (tlucua(NINT(TG*1000._dp)))/100._dp
              ENEW = 6.112_dp*EXP(17.67_dp*(TG-T0)/(243.5_dp+TG-T0))
              RG = eps*ENEW/(P_hpa(J)-ENEW)
           ENDDO ! K
           TLP(I,J) = TG
           TLVP(I,J) = TG*(1.0_dp+RG/eps)/(1.0_dp+RG)
           TVPDIF(I,J) = TLVP(I,J)-T(J)*(1.0_dp+R(J)/eps)/(1.0_dp+R(J))
           RG0 = RG
           TG0 = TG
           SUM = CPW
        ENDIF
     ENDDO 	! J

!
!  ***  Find positive and negative areas  PA and NA and
!       CAPE (=PA-NA) from pseudo-adiabatic ascent ***
!
!   ***  Find lifted condensation level and maximum level   ***
!   ***               of positive buoyancy                  ***
!
     INBP = nl	! index of maximum level of positive buoyancy
     DO J=1,I 
        IF (TVPDIF(I,J) > 0.0_dp) INBP = MIN(INBP,J)
     ENDDO
     IMIN = MIN(INBP-1,I)
     TVPDIF(I,1:IMIN)=0.0_dp	! set to zero for levels above IMIN

!
!   ***  Do updraft loops        ***
!
    IF (INBP < I) THEN 
       DO J=(I-1),INBP,-1
          TVM = 0.5_dp*(TVPDIF(I,J)+TVPDIF(I,J+1))
          PM = 0.5*(P_hpa(J)+P_hpa(J+1))
          IF (TVM <= 0.0) THEN
             NAP(I)=NAP(I)-rd*TVM*(P_hpa(J+1)-P_hpa(J))/PM
          ELSE
             PAP(I)=PAP(I)+rd*TVM*(P_hpa(J+1)-P_hpa(J))/PM
          ENDIF
       ENDDO
       CAPEP(I)=PAP(I)-NAP(I)
    ENDIF
    
 ENDDO	! I, loop over air parcel origins

! surface parcel
 CAPE_s = CAPEP(nl)
 PA_s = PAP(nl)

! parcel with maximum CAPE
 CAPE_c = maxval(CAPEP)
 PA_c = maxval(PAP)
 IDX(:) = maxloc(CAPEP)
    


END SUBROUTINE calc_cape_crm
! mz_ss_2013_1107-
!=============================================================================

  SUBROUTINE cloud_droplet_nc_crm(kproma, kbdim, nlev, pressure, acdncm, loland, loglac)

    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: nlev, kproma, kbdim
    REAL(dp), INTENT(IN) :: pressure(kbdim, nlev)
    LOGICAL,  INTENT(IN) :: loland(kbdim), loglac(kbdim)

    REAL(dp), INTENT(INOUT) :: acdncm(kbdim, nlev)
    
    REAL(dp) :: zn1, zn2, zcdnc, zprat
    INTEGER  :: jl,jk, nexp

    INTRINSIC :: MIN, EXP

    DO jk=1, nlev        
      DO jl=1, kproma       
        nexp=2
        zprat=(MIN(8.0_dp,80000.0_dp/pressure(jl,jk)))**nexp
        IF (loland(jl).AND.(.NOT.loglac(jl))) THEN
          zn1= 50.0_dp
          zn2=220.0_dp
        ELSE 
          zn1= 50.0_dp
          zn2= 80.0_dp
        ENDIF
        IF (pressure(jl,jk).LT.80000.0_dp) THEN
          zcdnc=1.0e6_dp*(zn1+(zn2-zn1)*(EXP(1.0_dp-zprat)))
        ELSE
          zcdnc=zn2*1.0e6_dp
        ENDIF
        acdncm(jl,jk)=zcdnc
      END DO
    END DO

  END SUBROUTINE cloud_droplet_nc_crm
!=============================================================================

END MODULE MESSY_CRM
