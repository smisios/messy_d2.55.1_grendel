#include "messy_main_ppd_bi.inc"

MODULE messy_sedi_si

  ! This module calculates aerosol sedimentation
  ! Modularized by Astrid Kerkweg, MPICH April 2004
  ! for arbitrary aerosol models
  !  without recompilation, Patrick Joeckel, MPCHI, 2004

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , error_bi, warning_bi, info_bi
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,    ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,         &    
                                      mtend_id_tracer,        &
                                      mtend_id_t,             &
                                      mtend_id_q
#endif
  ! MESSy
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_2D_ARRAY &
                                    , PTR_4D_ARRAY, PTR_5D_ARRAY

  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, LGTRSTR
  ! SUBMODEL SEDIMENTAION
  USE messy_sedi

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL, TRIM

  ! FIELDS FOR AEROSOL SEDIMENTATION
  INTEGER, PARAMETER :: NMAX_AEROMODELS = 10

  ! SWITCHES in CPL-NAMELIST
  LOGICAL :: L_LG = .FALSE.  ! sedimentation for Lagrangian tracers
  LOGICAL :: L_GP = .FALSE.  ! sedimentation for Gridpoint tracers

  ! DATASET NUMBERING
  INTEGER, PARAMETER :: I_GP = 1           ! GP is always SSET 1
  INTEGER, PARAMETER :: I_LG = 2           ! LG is always SSET 2
  CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: setname = (/GPTRSTR, LGTRSTR/)


  TYPE T_SEDI_SET
     INTEGER  :: reprid
     INTEGER  :: ntrac
     ! AEROSOL MODEL (ID) OF TRACER
     INTEGER, DIMENSION(:), POINTER :: idt_aer             ! aerosol indices
     INTEGER, DIMENSION(:), POINTER :: traermod => NULL()
     !
     ! indirect index to tracer ID; only used for LG representation for
     ! packing of aerosol tracers; GP: loop / cycle over all tracers
     INTEGER, DIMENSION(:), POINTER :: idx_aer => NULL()   ! tracer indices
     !
     INTEGER                                                  :: aeromodnum = 0
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(NMAX_AEROMODELS) :: aermodname
     INTEGER,                      DIMENSION(NMAX_AEROMODELS) :: aermodmethod
     TYPE (PTR_4D_ARRAY), DIMENSION(:), POINTER :: radii     => NULL()
     TYPE (PTR_4D_ARRAY), DIMENSION(:), POINTER :: densaer   => NULL()
     TYPE (PTR_1D_ARRAY), DIMENSION(:), POINTER :: sigma     => NULL()
     INTEGER,             DIMENSION(:), POINTER :: nmod      => NULL() 
!     LOGICAL,             DIMENSION(:), POINTER :: l_running => NULL()
     ! terminal velocity for mass distribution
     TYPE (PTR_5D_ARRAY), DIMENSION(:), POINTER :: vsedi_mass => NULL() 
     ! terminal velocity for number distribution
     TYPE (PTR_5D_ARRAY), DIMENSION(:), POINTER :: vsedi_num  => NULL()
     !
     ! mz_bk_20110725+
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: sedflux    => NULL()
     ! mz_bk_20110725-
     TYPE (PTR_2D_ARRAY), DIMENSION(:), POINTER :: sedfluxsum => NULL()
  END TYPE T_SEDI_SET
     
  TYPE(T_SEDI_SET), DIMENSION(:), POINTER :: SSET

  ! GP CHANNEL OBJECTS SHARED BETWEEN GP AND LG REPRESENTATION
  ! air temperature
  REAL(dp), DIMENSION(:,:,:), POINTER :: temp => NULL()
  ! mean free path of air
  REAL(dp), DIMENSION(:,:,:), POINTER :: lambda_air => NULL()
  ! density of air
  REAL(dp), DIMENSION(:,:,:), POINTER :: rho_air => NULL()
  ! viscosity of air
  REAL(dp), DIMENSION(:,:,:), POINTER :: visc_air => NULL()
  ! 
  ! layerthickness = delta height
  REAL(dp), DIMENSION(:,:,:), POINTER :: dz => NULL()
  ! layerthickness = delta pressure 
  REAL(dp), DIMENSION(:,:,:), POINTER :: delp => NULL()

  ! LG CHANNEL OBJECTS
  ! NOTE: Although being rank-3, these are not GP-dimensioned!
  !       Rank 2 and three (of dimension 1) are required in order to
  !       be able to re-use the SMCL routines, which are coded for
  !       column vectors.
  ! mean free path of air
  REAL(dp), DIMENSION(:,:,:), POINTER :: lambda_air_lg => NULL()
  ! density of air
  REAL(dp), DIMENSION(:,:,:), POINTER :: rho_air_lg => NULL()
  ! viscosity of air
  REAL(dp), DIMENSION(:,:,:), POINTER :: visc_air_lg => NULL()
  !
  ! SPECIAL REQUIREMENTS FOR LG representation
  ! number of parcels per grid-box
  REAL(DP), DIMENSION(:,:,:),     POINTER :: NCB         => NULL()   
  ! residual flux in GP representation
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: fx_residual => NULL()

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  PUBLIC :: sedi_initialize
  PUBLIC :: sedi_init_memory
  PUBLIC :: sedi_init_coupling
  PUBLIC :: sedi_physc
  PUBLIC :: sedi_global_end
  PUBLIC :: sedi_free_memory

CONTAINS

!-----------------------------------------------------------------------------
  SUBROUTINE sedi_initialize

    ! INITIALIZE SEDIMENTATION
    ! At THE MOMENT ONLY THE CHOICE OF ORDER IS POSSIBLE 

    ! Author:  Astrid Kerkweg, MPICH, Mar 2004

    ! MESSy/BMIL
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast 
    ! MESSy/SMCL
    USE messy_main_tools,     ONLY: find_next_free_unit
  
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='sedi_initialize'
    INTEGER                     :: status
    INTEGER                     :: iou    ! I/O unit

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL SEDI CORE ROUTINE:
       CALL sedi_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR read CTRL-namelist', substr)

       ! TEST NAMELIST: 
       SELECT CASE(order)
       CASE(0)
          CALL info_bi('Sedimentaion scheme of 0. order applied', substr)
       CASE(1)
          CALL info_bi('Sedimentaion scheme of 1. order applied', substr)
          SELECT CASE(scheme)
          CASE(1)
             CALL info_bi('Trapezoid scheme applied', substr)
          CASE(2)
             CALL info_bi('First order two domain scheme applied', substr)
          CASE(3)
             CALL info_bi('Second order two domain scheme applied', substr)
          CASE(4)
             CALL info_bi('Pseudo-equilibrium  two domain scheme applied' &
                  , substr)
          CASE(5)
             CALL info_bi('Pseudo-equilibrium with peak conservation two'//&
                  &' domain scheme applied', substr)
          CASE(6)
             CALL info_bi('Sedimentation adapted Walcek scheme applied', substr)
          CASE DEFAULT
             IF (scheme < 1 .OR. scheme > 6 ) &
                  CALL error_bi('wrong  sedimentation scheme chosen', substr)
          END SELECT

       CASE DEFAULT
          IF (order /= 1 .AND. order /= 0 ) &
               CALL error_bi('wrong order of sedimentation scheme chosen' &
               , substr)
       END SELECT

     ! INITIALIZE CPL
       iou = find_next_free_unit(100,200)
       CALL sedi_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR reading CPL-namelist', substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(l_lg, p_io)
    CALL p_bcast(l_gp, p_io)
    CALL p_bcast(order, p_io)
    CALL p_bcast(scheme, p_io)

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

  END SUBROUTINE sedi_initialize
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE sedi_init_memory

    ! MESSy/BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, nproma, ngpblks
    USE messy_main_tracer_mem_bi,    ONLY: ntrac_gp, ti_gp, ntrac_lg, ti_lg
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL &
                                         , LG_ATTILA, GP_3D_MID
    ! MESSy
    USE messy_main_tracer,           ONLY: AEROSOL, ON, S_aerosol_model &
                                         , I_aerosol_method, I_sedi     &
                                         , t_trinfo_tp, NUMBERDENSITY   &
                                         , AMOUNTFRACTION
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr ='sedi_init_memory'
    INTEGER :: status
    INTEGER :: jt, naermo, i
    LOGICAL :: lnew
    INTEGER :: na = 0  ! number of aerosol trac.
    INTEGER :: idx
    REAL(dp), DIMENSION(:,:,:,:),    POINTER :: mem => NULL()
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti  => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM)             :: unit

    CALL start_message_bi(modstr, 'INIT MEMORY ',substr) 

    !------------------------------------------------------------------------
    !   GENERAL PART USED BY LAGRANGE AND GRID POINT PART
    !------------------------------------------------------------------------
    
    ! (1) TEST AND ALLOCATE NUMBER OF SEDI TRACER SETS

    IF (L_GP .AND. ntrac_gp == 0) THEN
       CALL warning_bi('L_GP is .TRUE., but no GP-Tracer available', substr)
       CALL warning_bi('SETTING L_GP to .FALSE.', substr)
       L_GP = .FALSE.
    ENDIF
    IF (L_LG .AND. ntrac_lg == 0) THEN
       CALL warning_bi('L_LG is .TRUE., but no LG-Tracer available', substr)
       CALL warning_bi('SETTING L_LG to .FALSE.', substr)
       L_LG = .FALSE.
    ENDIF
    
    IF (L_GP .AND. L_LG) THEN
       ALLOCATE(SSET(2))
    ELSE IF (L_GP) THEN
       ALLOCATE(SSET(I_GP:I_GP))
    ELSE IF (L_LG) THEN
       ALLOCATE(SSET(I_LG:I_LG))
    ELSE
       CALL warning_bi('NO TRACER AVAILABLE NEITHER GP NOR LG',substr)
       CALL error_bi('PLEASE SWITCH OF SEDI or implement new tracer set',substr)
       RETURN
    ENDIF

    ! NEW CHANNEL (GP AND LG)
    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr  &
         , 'lambda_air', p3=lambda_air)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'lambda_air', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'lambda_air' &
         , 'long_name', c='mean free path of air')
    
    CALL new_channel_object(status, modstr  &
         , 'rho_air', p3=rho_air)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'rho_air', 'units', c='kg/m^3')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'rho_air' &
         , 'long_name', c='density of air')

    CALL new_channel_object(status, modstr  &
         , 'visc_air', p3=visc_air)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'visc_air', 'units', c='kg/(m s)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'visc_air' &
         , 'long_name', c='viscosity of air')

    CALL new_channel_object(status, modstr  &
         , 'dz', p3=dz)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'dz', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'dz' &
         , 'long_name', c='layer thickness')

    CALL new_channel_object(status, modstr  &
         , 'delp', p3=delp)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'delp', 'units', c='Pa')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'delp' &
         , 'long_name', c='layer thickness')

    CALL new_channel_object(status, modstr  &
         , 'temp', p3=temp)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr &
         , 'temp', 'units', c='K')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr  &
         , 'temp' &
         , 'long_name', c='air temperature')

    ! (2) INTERNAL UPDATE OF SWITCHES
    
    tracer_set_loop: DO i=1,2
       
       IF (i==I_GP) THEN
          IF (.NOT.L_GP) CYCLE
          ti => ti_gp
          SSET(i)%ntrac = ntrac_gp
          SSET(i)%reprid = GP_2D_HORIZONTAL
       ENDIF
       IF (i==I_LG) THEN
          IF (.NOT. L_LG) CYCLE
          ti => ti_lg
          SSET(i)%ntrac = ntrac_lg
          SSET(i)%reprid =LG_ATTILA
       ENDIF
       
       ! (2.1) SCAN FOR RUNNING AEROSOL-MODELS
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Scanning '//setname(i)//'-tracers for aerosol models ...'
       END IF
       
       ALLOCATE(SSET(i)%traermod(SSET(i)%ntrac))
       SSET(i)%traermod(:) = 0
       
       DO jt=1, SSET(i)%ntrac
          ! CHECK IF AEROSOL MODEL
          IF (TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) == '') CYCLE
          lnew = .TRUE.
          ! CHECK IF AEROSOL MODEL ALREADY IN LIST
          DO naermo = 1,  SSET(i)%aeromodnum
             IF ( TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) &
                  == SSET(i)%aermodname(naermo) ) THEN
                lnew = .FALSE.
                EXIT
             END IF
          END DO
          ! NEW ENTRY IN LIST
          IF (lnew) THEN
             SSET(i)%aeromodnum = SSET(i)%aeromodnum + 1
             IF (SSET(i)%aeromodnum > NMAX_AEROMODELS) &
                  CALL error_bi( &
                  'RECOMPILATION WITH INCREASED NMAX_AEROMODELS REQUIRED', &
                  substr)
             SSET(i)%aermodname(SSET(i)%aeromodnum) = &
                  TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model))
             SSET(i)%aermodmethod(SSET(i)%aeromodnum) = &
                  ti(jt)%tp%meta%cask_i(I_aerosol_method)
             naermo = SSET(i)%aeromodnum
          END IF
          SSET(i)%traermod(jt) = naermo
       END DO
       
       IF (p_parallel_io) THEN
          WRITE(*,*) '  ... ',SSET(i)%aeromodnum,' aerosol models found:'
          WRITE(*,*) SSET(i)%aermodname(1:SSET(i)%aeromodnum)
       END IF
       
       ! (2) INITIALIZE MEMORY
       ! (2a) ALLOCATE POINTERS
       ALLOCATE(SSET(i)%radii(SSET(i)%aeromodnum))
       ALLOCATE(SSET(i)%densaer(SSET(i)%aeromodnum))
       ALLOCATE(SSET(i)%sigma(SSET(i)%aeromodnum))
       ALLOCATE(SSET(i)%nmod(SSET(i)%aeromodnum))
!       ALLOCATE(SSET(i)%l_running(SSET(i)%aeromodnum))
       ALLOCATE(SSET(i)%vsedi_mass(SSET(i)%aeromodnum))
       ALLOCATE(SSET(i)%vsedi_num(SSET(i)%aeromodnum))
       
       ! (2b) NEW CHANNEL
       CALL new_channel(status, modstr//'_'//setname(i) &
                          , reprid=SSET(i)%reprid)
       CALL channel_halt(substr, status)
       
       ! (2c) CHANNEL OBJECTS
       ALLOCATE(SSET(i)%idt_aer(SSET(i)%ntrac))
       SSET(i)%idt_aer(:) = 0
       !
       ! COUNT SEDIMENTING AEROSOL TRACERS
       na = 0
       DO jt=1, SSET(i)%ntrac
          IF ( (ti(jt)%tp%ident%medium == AEROSOL) .AND. &
               (ti(jt)%tp%meta%cask_i(I_sedi)== ON) ) THEN
             na = na + 1
             SSET(i)%idt_aer(jt) = na
          END IF
       END DO
       ALLOCATE(SSET(i)%sedflux(na))
       ALLOCATE(SSET(i)%sedfluxsum(na))
       ALLOCATE(SSET(i)%idx_aer(na))
       SSET(i)%idx_aer(:) = 0
       !
       na = 0
       DO jt=1, SSET(i)%ntrac
          IF ( (ti(jt)%tp%ident%medium == AEROSOL) .AND. &
               (ti(jt)%tp%meta%cask_i(I_sedi)== ON) ) THEN
             na = na + 1
             SSET(i)%idx_aer(na) = jt
             CALL new_channel_object(status, modstr//'_'//setname(i)          &
                  , 'sediflux_'//TRIM(ti(jt)%tp%ident%fullname)               &
                  , p2 = SSET(i)%sedflux(SSET(i)%idt_aer(jt))%ptr              &
                  , reprid = GP_2D_HORIZONTAL                                 &
                  , lrestreq = .TRUE.                                         &
                  )
             CALL channel_halt(substr, status)
             SELECT CASE(ti(jt)%tp%ident%quantity)
             CASE (NUMBERDENSITY)
                unit = '# m-2 s-1'
             CASE (AMOUNTFRACTION)
                unit = 'molec m-2 s-1'
             CASE DEFAULT
                unit = 'unknown'
             END SELECT
             CALL new_attribute(status, modstr//'_'//setname(i)               &
                  , 'sediflux_'//TRIM(ti(jt)%tp%ident%fullname)               &
                  , 'units', c=TRIM(unit)                                     &
                  )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_'//setname(i)               &
                  , 'sediflux_'//TRIM(ti(jt)%tp%ident%fullname)               &
                  , 'long_name'                                               &
                  , c = 'sedimentation flux for '                             &
                  &//TRIM(ti(jt)%tp%ident%fullname)                           &
                  )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr//'_'//setname(i)  &
                  , 'loss_'//TRIM(ti(jt)%tp%ident%fullname) &
                  , p2=SSET(i)%sedfluxsum(SSET(i)%idt_aer(jt))%ptr  &
                  , reprid = GP_2D_HORIZONTAL       &  ! also for LG !!!
                  , lrestreq = .TRUE.               &
                  )
             CALL channel_halt(substr, status)
             SELECT CASE(ti(jt)%tp%ident%quantity)
             CASE(NUMBERDENSITY)
                unit = '# m-2'
             CASE(AMOUNTFRACTION)
                unit = 'kg(tracer) m-2'
             CASE DEFAULT
                unit = 'unknown'
             END SELECT
             CALL new_attribute(status, modstr//'_'//setname(i)  &
                  , 'loss_'//TRIM(ti(jt)%tp%ident%fullname) &
                  , 'units', c=TRIM(unit))
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_'//setname(i)  &
                  , 'loss_'//TRIM(ti(jt)%tp%ident%fullname) &
                  , 'long_name', c='time integral of loss for '&
                  &//TRIM(ti(jt)%tp%ident%fullname) &
                  )
             CALL channel_halt(substr, status)
          END IF
       END DO

    ENDDO tracer_set_loop

#ifdef MESSYTENDENCY
    IF (L_GP) THEN
       CALL mtend_register(my_handle, mtend_id_vec=SSET(I_GP)%idx_aer(:))
    END IF
#endif

    ! ADDITIONAL WORKSPACE FOR LG REPRESENTATION
    IF (L_LG) THEN

       CALL new_channel_object(status, modstr//'_'//setname(I_LG)  &
            , 'lambda_air', p3=lambda_air_lg)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//setname(I_LG) &
            , 'lambda_air', 'units', c='m')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//setname(I_LG)  &
            , 'lambda_air' &
            , 'long_name', c='mean free path of air')
       
       CALL new_channel_object(status, modstr//'_'//setname(I_LG)  &
            , 'rho_air', p3=rho_air_lg)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//setname(I_LG) &
            , 'rho_air', 'units', c='kg/m^3')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//setname(I_LG)  &
            , 'rho_air' &
            , 'long_name', c='density of air')
       
       CALL new_channel_object(status, modstr//'_'//setname(I_LG)  &
            , 'visc_air', p3=visc_air_lg)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//setname(I_LG) &
            , 'visc_air', 'units', c='kg/(m s)')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//setname(I_LG)  &
            , 'visc_air' &
            , 'long_name', c='viscosity of air')
       
       ! RESIDUAL FLUX
       na = SIZE(SSET(I_LG)%idx_aer)
       ALLOCATE(fx_residual(nproma,nlev,na,ngpblks,1))
       fx_residual(:,:,:,:,:) = 0.0_dp

       DO jt=1, na
          idx = SSET(I_LG)%idx_aer(jt)
          mem => fx_residual(:,:,jt,:,:)
          CALL new_channel_object(status, modstr//'_'//setname(I_LG)  &
               , 'fxres_'//TRIM(ti(idx)%tp%ident%fullname)            &
               , reprid = GP_3D_MID          &
               , mem=mem, lrestreq = .TRUE.  &
               )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_'//setname(I_LG)  &
               , 'fxres_'//TRIM(ti(idx)%tp%ident%fullname) &
               , 'units', c='(mol/mol) * kg/(m^2 s)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_'//setname(I_LG)  &
               , 'fxres_'//TRIM(ti(idx)%tp%ident%fullname) &
               , 'long_name', c='residual flux of '&
               &//TRIM(ti(idx)%tp%ident%fullname) &
               )
          CALL channel_halt(substr, status)
       END DO

    END IF

    CALL end_message_bi(modstr, 'INIT MEMORY ',substr) 

  END SUBROUTINE sedi_init_memory
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE sedi_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, nlev, ngpblks
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, LG_ATTILA
    USE messy_main_tracer_mem_bi,    ONLY: NCELL
    ! MESSy
    USE messy_main_tracer,        ONLY: BIN
    USE messy_main_channel,       ONLY: get_channel_info, get_channel_object &
                                      , STRLEN_OBJECT, new_channel_object &
                                      , new_attribute
    USE messy_main_tools,         ONLY: int2str

    IMPLICIT NONE

    INTRINSIC :: TRIM, NULL, SIZE

    ! LOCAL
    ! auxiliary pointer for data managment
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: mem => NULL()
    
    CHARACTER(len=*), PARAMETER    :: substr='sedi_init_coupling'
    INTEGER                        :: status
    INTEGER                        :: jmod, modnum
    INTEGER                        :: naermo,i
    CHARACTER(LEN=STRLEN_OBJECT)   :: name
    CHARACTER(LEN=STRLEN_MEDIUM+3) :: modname
    CHARACTER(LEN=2)               :: modestr
    INTEGER                        :: modnumpos, REPRID, d1, d2, d3
    
    trset_loop: DO i=1,2 
    
       IF ((i==I_GP) .AND. (.NOT.L_GP)) CYCLE
       IF ((i==I_LG) .AND. (.NOT.L_LG)) CYCLE

       CALL start_message_bi(modstr, 'INIT COUPLING '//setname(i),substr) 

       aermo_loop: DO naermo = 1, SSET(i)%aeromodnum
          SSET(i)%nmod(naermo)            =  0
          SSET(i)%radii(naermo)%ptr       => NULL()
          SSET(i)%densaer(naermo)%ptr     => NULL()
          SSET(i)%sigma(naermo)%ptr       => NULL() 
          SSET(i)%vsedi_mass(naermo)%ptr  => NULL()
          SSET(i)%vsedi_num(naermo)%ptr   => NULL()
       
          modname = TRIM(SSET(i)%aermodname(naermo))//'_'//setname(i)
          CALL get_channel_info(status, TRIM(modname))
          IF (status /= 0) THEN
              CALL error_bi(&
                  'requested aerosol model '//TRIM(modname)//' not running!'&
                  , substr)
          END IF
          
          CALL info_bi('calculate sedimentation for '//modname, substr)
          !
          IF (i==I_GP) THEN
             modnumpos = _IN_XYZN_
             REPRID = GP_3D_MID
             d1 = nproma
             d2 = ngpblks
             d3 = nlev
          END IF
          IF (i==I_LG) THEN
             modnumpos = 2
             REPRID = LG_ATTILA
             d1 = NCELL
             d2 = 1
             d3 = 1
          END IF

          CALL get_channel_object(status, TRIM(modname) &
               , 'wetradius' , p4=SSET(i)%radii(naermo)%ptr )
          IF (status /= 0) &
               CALL error_bi( 'wetradius not available for '&
               &//TRIM(modname), substr) 
          !
          CALL get_channel_object(status, TRIM(modname) &
               , 'densaer', p4 = SSET(i)%densaer(naermo)%ptr )
          IF (status /= 0) &
               CALL error_bi('densaer not available for '&
               &//TRIM(modname), substr)
          modnum = SIZE(SSET(i)%densaer(naermo)%ptr,modnumpos)
          SSET(i)%nmod(naermo) = modnum 

          IF (SSET(i)%aermodmethod(naermo) /= BIN) THEN
             ALLOCATE(SSET(i)%vsedi_mass(naermo)%ptr(&
                  _RI_XYZN_(d1, d2, d3, SSET(i)%nmod(naermo)), 1))
             SSET(i)%vsedi_mass(naermo)%ptr(:,:,:,:,:) = 0.0_dp
          END IF
          
          ALLOCATE(SSET(i)%vsedi_num(naermo)%ptr(&
               _RI_XYZN_(d1, d2, d3, SSET(i)%nmod(naermo)), 1))
          SSET(i)%vsedi_num(naermo)%ptr(:,:,:,:,:) = 0.0_dp

          mod_loop: DO jmod = 1, SSET(i)%nmod(naermo)
             
             CALL int2str(modestr, jmod, '0', 'X')
             
             IF (SSET(i)%aermodmethod(naermo) /= BIN) THEN
                name = TRIM(SSET(i)%aermodname(&
                     naermo))//'_v_mass_m'//TRIM(modestr)
                mem => SSET(i)%vsedi_mass(naermo)%ptr(_RI_XYZN_(:,:,:,jmod),:)
                CALL new_channel_object(status, modstr//'_'//setname(i) &
                     , TRIM(name), reprid = REPRID, mem=mem )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//'_'//setname(i) &
                     , TRIM(name), 'units', c='m/s')
                CALL channel_halt(substr, status)          
             ENDIF
             
             name = TRIM(SSET(i)%aermodname(&
                  naermo))//'_v_num_m'//TRIM(modestr)
             mem => SSET(i)%vsedi_num(naermo)%ptr(_RI_XYZN_(:,:,:,jmod),:)
             CALL new_channel_object(status, modstr//'_'//setname(i) &
                  , TRIM(name), reprid = REPRID, mem=mem )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_'//setname(i) &
                  , TRIM(name), 'units', c='m/s')
             CALL channel_halt(substr, status)     
             
          END DO mod_loop
          
          CALL get_channel_object(status, TRIM(modname) &
               , 'sigma', p1 = SSET(i)%sigma(naermo)%ptr )
          IF (status /= 0) &
               CALL error_bi( 'sigma not available for '&
               &//TRIM(modname), substr)
          
       END DO aermo_loop
       
       CALL end_message_bi(modstr, 'INIT COUPLING '//setname(i),substr) 

    END DO trset_loop

    IF (L_LG) THEN
       CALL get_channel_object(status, 'attila','NCB', p3=NCB)
       IF (status /= 0) &
            CALL error_bi(&
            'no ATTILA channel object for number of cells per grid box found!'&
            , substr)
    ENDIF

  END SUBROUTINE sedi_init_coupling
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE sedi_physc

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev, nproma
    USE messy_main_data_bi,        ONLY: press_3d, pressi_3d              
    USE messy_main_tracer_mem_bi,  ONLY: ti_gp                  
    USE messy_main_timer,          ONLY: time_step_len, delta_time
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1
    USE messy_main_data_bi,        ONLY: tm1_3d , tte_3d, qm1_3d , qte_3d 
#endif

    ! MESSy
    USE messy_main_tracer,         ONLY: AEROSOL, BIN                      &
                                       , AMOUNTFRACTION, NUMBERDENSITY, ON &
                                       , R_molarmass &
                                       , I_aerosol_mode, I_sedi
    USE messy_main_constants_mem,  ONLY: g, M_air, N_A
    USE messy_main_tools,          ONLY: mass_density

    IMPLICIT NONE
    
    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='sedi_physc'
    REAL(DP), PARAMETER :: conv_mr = N_A / (1.E-03_dp * M_air)
    REAL(DP), PARAMETER :: conv_nu = 1._dp / (1.E-03_dp * M_air)
    REAL(dp) :: conv
    REAL(dp) :: ztmst, zdtime
    REAL(dp) :: sphum(nproma,nlev)
    INTEGER  :: naermod ! actual aerosol module number
    INTEGER  :: jt, jm
    REAL(dp) :: vt_tracer(nproma,nlev)  !sedimentaion velocity of tracer
    REAL(dp) :: tendency(nproma,nlev)   ! sedimentation tendency
    REAL(dp) :: zmr(nproma,nlev)        ! sedimentation tendency
    REAL(dp) :: sfloss(nproma) ! loss at surface [mol/mol/s]
    REAL(DP) :: conv_sum
#ifdef MESSYTENDENCY
    REAL(DP) :: getval(nproma,nlev)
#endif

    ztmst  = time_step_len
    zdtime = delta_time

    ! calculate temperature
#ifndef MESSYTENDENCY
    temp(_RI_XYZ__(1:kproma,jrow,1:nlev)) = &
         tm1_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) + & 
         tte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) * ztmst
#else
    CALL mtend_get_start_l(mtend_id_t, v0=getval)
    temp(_RI_XYZ__(1:kproma,jrow,1:nlev)) = getval(1:kproma,1:nlev)
#endif

    ! calculate specific humidity
#ifndef MESSYTENDENCY
    sphum(1:kproma,1:nlev) = qm1_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) + & 
                             qte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) * ztmst
#else
    CALL mtend_get_start_l(mtend_id_q, v0=sphum)
#endif

    ! calculate air density in kg/m3
    rho_air(_RI_XYZ__(1:kproma,jrow,1:nlev)) = &
         mass_density(press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))  &
         ,temp(_RI_XYZ__(1:kproma,jrow,1:nlev))             &
         ,sphum(1:kproma,1:nlev)  )

    ! calculate mean free path of air 
    lambda_air(_RI_XYZ__(1:kproma,jrow,1:nlev))= &
         mean_free_path(press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) &
         ,temp(_RI_XYZ__(1:kproma,jrow,1:nlev)) )

    ! calculate viscosity of air 
    CALL air_viscosity(kproma, nlev, temp(_RI_XYZ__(1:kproma,jrow,1:nlev)),  &
         visc_air(_RI_XYZ__(1:kproma,jrow,1:nlev)))

    ! calculate box height dz and delp
    delp(_RI_XYZ__(1:kproma,jrow,1:nlev)) = &
         pressi_3d(_RI_XYZ__(1:kproma,jrow,2:nlev+1)) &
         - pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) 
    dz(_RI_XYZ__(1:kproma,jrow,1:nlev)) = delp(_RI_XYZ__(1:kproma,jrow,1:nlev))/rho_air(_RI_XYZ__(1:kproma,jrow,1:nlev))/g

    IF (.NOT. L_GP) RETURN

    aermod_loop: DO naermod = 1, SSET(I_GP)%aeromodnum
       IF (SSET(I_GP)%aermodmethod(naermod) == BIN) THEN
          bin_loop: DO jm = 1, SSET(I_GP)%nmod(naermod)
             ! calculate terminal velocity for mass distribution
             CALL calc_vt(SSET(I_GP)%vsedi_num(naermod)%ptr&
                  (_RI_XYZN_(1:kproma,jrow,1:nlev,jm),1) &
                  , kproma, nlev                                        &
                  , SSET(I_GP)%radii(naermod)%ptr(&
                                  _RI_XYZN_(1:kproma,jrow,1:nlev,jm))   &
                  , SSET(I_GP)%densaer(naermod)%ptr(&
                                    _RI_XYZN_(1:kproma,jrow,1:nlev,jm)) &
                  , lambda_air(_RI_XYZ__(1:kproma,jrow,1:nlev))         &
                  , rho_air(_RI_XYZ__(1:kproma,jrow,1:nlev))            & 
                  , visc_air(_RI_XYZ__(1:kproma,jrow,1:nlev))           & 
                  , SSET(I_GP)%sigma(naermod)%ptr(jm)   ,0 )
          ENDDO bin_loop
       ELSE
          mod_loop: DO jm = 1, SSET(I_GP)%nmod(naermod)! um_ak_20100802 jmod =jm
             ! calculate terminal velocity for mass distribution
             CALL calc_vt( &
                  SSET(I_GP)%vsedi_mass(naermod)%ptr(&
                                     _RI_XYZN_(1:kproma,jrow,1:nlev,jm),1)  &
                  , kproma, nlev                                            &
                  , SSET(I_GP)%radii(naermod)%ptr(&
                                  _RI_XYZN_(1:kproma,jrow,1:nlev,jm))       &
                  , SSET(I_GP)%densaer(naermod)%ptr(&
                                    _RI_XYZN_(1:kproma,jrow,1:nlev,jm))     &
                  , lambda_air(_RI_XYZ__(1:kproma,jrow,1:nlev))             &
                  , rho_air(_RI_XYZ__(1:kproma,jrow,1:nlev))                &
                  , visc_air(_RI_XYZ__(1:kproma,jrow,1:nlev))               &
                  , SSET(I_GP)%sigma(naermod)%ptr(jm)   ,1 )
             ! calculate terminal velocity for number distribution
             CALL calc_vt( &
                  SSET(I_GP)%vsedi_num(naermod)%ptr(&
                                    _RI_XYZN_(1:kproma,jrow,1:nlev,jm),1) &
                  , kproma, nlev                                          &
                  , SSET(I_GP)%radii(naermod)%ptr(&
                                  _RI_XYZN_(1:kproma, jrow,1:nlev,jm))    &
                  , SSET(I_GP)%densaer(naermod)%ptr(&
                                    _RI_XYZN_(1:kproma,jrow,1:nlev,jm))   &
                  , lambda_air(_RI_XYZ__(1:kproma,jrow,1:nlev))           &
                  , rho_air(_RI_XYZ__(1:kproma,jrow,1:nlev))              & 
                  , visc_air(_RI_XYZ__(1:kproma,jrow,1:nlev))             &
                  , SSET(I_GP)%sigma(naermod)%ptr(jm)  ,2 )
          END DO mod_loop
       END IF
    END DO aermod_loop

    tracer_loop_gp: DO jt=1,SSET(I_GP)%ntrac
       if_sedi_gp: IF (ti_gp(jt)%tp%meta%cask_i(I_sedi)== ON .AND. &
            ti_gp(jt)%tp%ident%medium == AEROSOL) THEN

          naermod = SSET(I_GP)%traermod(jt)
          IF (naermod == 0) CYCLE

#ifndef MESSYTENDENCY          
          zmr(1:kproma,1:nlev) = pxtm1(_RI_X_ZN_(1:kproma,1:nlev,jt)) + &
               pxtte(_RI_X_ZN_(1:kproma,1:nlev,jt)) * ztmst
#else
          CALL mtend_get_start_l(jt, v0=zmr)
#endif

          jm = ti_gp(jt)%tp%meta%cask_i(I_aerosol_mode)
          IF (SSET(I_GP)%aermodmethod(naermod) /= BIN ) THEN
             IF (ti_gp(jt)%tp%ident%quantity == AMOUNTFRACTION) &
                  vt_tracer(1:kproma,1:nlev) =                  &
                  SSET(I_GP)%vsedi_mass(naermod)%ptr &
                  (_RI_XYZN_(1:kproma,jrow,1:nlev,jm),1)
             IF (ti_gp(jt)%tp%ident%quantity == NUMBERDENSITY)  &
                  vt_tracer(1:kproma,1:nlev) =                  &
                  SSET(I_GP)%vsedi_num(naermod)%ptr &
                  (_RI_XYZN_(1:kproma,jrow,1:nlev,jm),1)
          ELSE
             vt_tracer(1:kproma,1:nlev) =            &
                  SSET(I_GP)%vsedi_num(naermod)%ptr &
                  (_RI_XYZN_(1:kproma,jrow,1:nlev,jm),1)
          END IF
          
          ! CALCULATE TOTAL SEDIMENTATION TENDENCY [mol/mol/s]
          SELECT CASE  (order) 
          CASE(0)
             CALL sedi_0( kproma, nlev                                   &
                  , tendency(1:kproma,1:nlev), zmr(1:kproma,1:nlev)      &
                  , vt_tracer(1:kproma,1:nlev)                           &
                  , dz(_RI_XYZ__(1:kproma,jrow,1:nlev))                  &
                  , delp(_RI_XYZ__(1:kproma,jrow,1:nlev)), ztmst         &
                  , sfloss(1:kproma) )
          CASE(1) 
             CALL sedi_1( kproma, nlev                          &
                  , press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))   &
                  , temp(_RI_XYZ__(1:kproma,jrow,1:nlev))       &
                  , vt_tracer(1:kproma,1:nlev),  ztmst          &
                  , zmr(1:kproma,1:nlev)                        &
                  , pressi_3d(_RI_XYZ__(1:kproma,jrow,:))       &
                  , tendency(1:kproma,1:nlev), sfloss(1:kproma) )
          CASE DEFAULT
             CALL error_bi('UNKOWN ORDER OF SEDIMENTATION SCHEME (GP)', substr)
          END SELECT
          
          ! ADD SEDIMENTATION TENDENCY
#ifndef MESSYTENDENCY
          pxtte(_RI_X_ZN_(1:kproma,1:nlev,jt)) =  &
               pxtte(_RI_X_ZN_(1:kproma,1:nlev,jt))+ tendency(1:kproma,:)
#else
          CALL mtend_add_l(my_handle, jt, px=tendency)
#endif          

          ! TIME INTEGRAL OF LOSS AT SURFACE [kg/m^2]
          SELECT CASE(ti_gp(jt)%tp%ident%quantity)
          CASE(AMOUNTFRACTION)
             conv = conv_mr
             conv_sum = (ti_gp(jt)%tp%meta%cask_r(R_molarmass)/M_air)
          CASE(NUMBERDENSITY)
             conv = conv_nu
             conv_sum = 1.E-3_dp/M_air
          CASE DEFAULT
             conv = 0.0_dp
             conv_sum = 0.0_dp
          END SELECT

          SSET(I_GP)%sedflux(SSET(I_GP)%idt_aer(jt))%ptr(1:kproma,jrow) =     &
               sfloss(1:kproma) * conv                                        &
               * (delp(_RI_XYZ__(1:kproma,jrow,nlev))/g)

          SSET(I_GP)%sedfluxsum(SSET(I_GP)%idt_aer(jt))%ptr(1:kproma,jrow) = &
             SSET(I_GP)%sedfluxsum(SSET(I_GP)%idt_aer(jt))%ptr(1:kproma,jrow)&
               + sfloss(1:kproma) * zdtime                 &
               * conv_sum &
               * (delp(_RI_XYZ__(1:kproma,jrow,nlev))/g)
          
       END IF if_sedi_gp
    END DO tracer_loop_gp
    
  END SUBROUTINE sedi_physc
!------------------------------------------------------------------------------
         
!------------------------------------------------------------------------------
SUBROUTINE sedi_global_end

!!#D attila +
#ifdef ECHAM5
  ! BMIL
  USE messy_main_tracer_mem_bi,  ONLY: ti_lg, NCELL                  &
                                     , pxtte => qxtte_a, pxtm1 => qxtm1_a
  USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, ngpblks &
                                     , npromz
  USE messy_main_timer,          ONLY: delta_time, time_step_len
  ! SMIL
  USE messy_attila_tools_e5,     ONLY: gp2lg_e5, lg2gp_e5, LG2GP_AVE &
                                     , lggpte2lgte_e5
  ! SMCL
  USE messy_main_tracer,         ONLY: BIN, AMOUNTFRACTION, NUMBERDENSITY &
                                     , I_AEROSOL_MODE &
                                     , R_molarmass
  USE messy_main_constants_mem,  ONLY: M_air, g, N_A

  IMPLICIT NONE

  ! LOCAL
  CHARACTER(len=*), PARAMETER :: substr='sedi_global_end'

  REAL(DP), PARAMETER :: conv_mr = N_A / (1.E-03_dp * M_air)
  REAL(DP), PARAMETER :: conv_nu = 1._dp / (1.E-03_dp * M_air)
  REAL(dp) :: conv

  INTEGER  :: naermod, jmod, zmode, jt, jk, idx, na, jp, jjrow, kproma
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: vt_gp     => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: xt_aer_gp => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: sfloss    => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: tend_gp   => NULL()
  REAL(DP), DIMENSION(:,:),     POINTER :: tend_lg   => NULL()
  REAL(DP), DIMENSION(:,:),     POINTER :: vt_lg     => NULL()
  REAL(DP), DIMENSION(:,:),     POINTER :: xt_aer_lg => NULL()
  REAL(DP), DIMENSION(:),       POINTER :: tmp_1d    => NULL()
  REAL(DP) :: ztmst, zdtime
  REAL(DP) :: conv_sum
  
  IF (.NOT. L_LG) RETURN

  ! INIT
  na = SIZE(SSET(I_LG)%idx_aer)
  ztmst  = time_step_len
  zdtime = delta_time

  ! MEMORY
  ALLOCATE(vt_gp(nproma, nlev, na, ngpblks))
  ALLOCATE(xt_aer_gp(nproma, nlev, na, ngpblks))
  xt_aer_gp(:,:,:,:) = 0.0_dp
  ALLOCATE(sfloss(nproma, na, ngpblks))
  ALLOCATE(tend_gp(nproma, nlev, na, ngpblks))
  ALLOCATE(tend_lg(NCELL, na))
  ALLOCATE(vt_lg(NCELL, na))
  ALLOCATE(xt_aer_lg(NCELL, na))

  ! CALCULATE THE SEDIMENTAION VELOCITIES (IN LG REPRESENTATION)
  tmp_1d => rho_air_lg(:,1,1)
  CALL gp2lg_e5(rho_air,tmp_1d)

  tmp_1d => visc_air_lg(:,1,1)
  CALL gp2lg_e5(visc_air,tmp_1d)

  tmp_1d => lambda_air_lg(:,1,1)
  CALL gp2lg_e5(lambda_air,tmp_1d)

  aermod_loop: DO naermod = 1, SSET(I_LG)%aeromodnum
     IF (SSET(I_LG)%aermodmethod(naermod) == BIN) THEN
        bin_loop: DO jmod = 1, SSET(I_LG)%nmod(naermod)
           ! calculate terminal velocity for mass distribution
           CALL calc_vt(SSET(I_LG)%vsedi_num(naermod)%ptr(:,1,jmod,:,1) &
                , NCELL, 1                                              &
                , SSET(I_LG)%radii(naermod)%ptr(:, jmod, :, 1)   &
                , SSET(I_LG)%densaer(naermod)%ptr(:, jmod, :, 1) &
                , lambda_air_lg(:, :, 1)                         &
                , rho_air_lg(:, :, 1)                            & 
                , visc_air_lg(:, :, 1)                           & 
                , SSET(I_LG)%sigma(naermod)%ptr(jmod)   ,0 )
        ENDDO bin_loop
     ELSE
        mod_loop: DO jmod = 1, SSET(I_LG)%nmod(naermod)
           ! calculate terminal velocity for mass distribution
           CALL calc_vt( &
                SSET(I_LG)%vsedi_mass(naermod)%ptr(:,1,jmod,:,1)  &
                , NCELL, 1                                        &
                , SSET(I_LG)%radii(naermod)%ptr(:, jmod, :, 1)    &
                , SSET(I_LG)%densaer(naermod)%ptr(:, jmod, :, 1)  &
                , lambda_air_lg(:, :, 1)                          &
                , rho_air_lg(:, :, 1)                             &
                , visc_air_lg(:, :, 1)                            &
                , SSET(I_LG)%sigma(naermod)%ptr(jmod)   ,1 )
           
           ! calculate terminal velocity for number distribution
           CALL calc_vt( &
                SSET(I_LG)%vsedi_num(naermod)%ptr(:,1,jmod,:,1)  &
                , NCELL, 1                                       &
                , SSET(I_LG)%radii(naermod)%ptr(:, jmod, :, 1)   &
                , SSET(I_LG)%densaer(naermod)%ptr(:, jmod, :, 1) &
                , lambda_air_lg(:, :, 1)                         &
                , rho_air_lg(:, :, 1)                            & 
                , visc_air_lg(:, :, 1)                           &
                , SSET(I_LG)%sigma(naermod)%ptr(jmod)  ,2 )
        END DO mod_loop
     END IF
  END DO aermod_loop

  ! - CALCULATE THE INITIAL CONDITION (IN LG REPRESENTATION)
  ! - PACK AEROSOL TRACERS
  ! - PACK THE SEDIMENTATION VELOCITY
  DO jt=1, na
     idx = SSET(I_LG)%idx_aer(jt)
     xt_aer_lg(:,jt) = pxtm1(:, idx) + pxtte(:, idx) * ztmst

     naermod = SSET(I_LG)%traermod(idx)
     IF (naermod == 0) CYCLE
     
     zmode = ti_lg(idx)%tp%meta%cask_i(I_aerosol_mode)
     
     IF (SSET(I_LG)%aermodmethod(naermod) /= BIN ) THEN
        IF (ti_lg(idx)%tp%ident%quantity == AMOUNTFRACTION) &
             vt_lg(:,jt) = &
             SSET(I_LG)%vsedi_mass(naermod)%ptr(:, 1, zmode, 1, 1)
        IF (ti_lg(idx)%tp%ident%quantity == NUMBERDENSITY)  &
             vt_lg(:,jt) = &
             SSET(I_LG)%vsedi_num(naermod)%ptr(:, 1, zmode, 1, 1)
     ELSE
        vt_lg(:,jt) = &
             SSET(I_LG)%vsedi_num(naermod)%ptr(:, 1, zmode, 1, 1)
     END IF
  END DO

  ! TRANSFORM VELOCITIES AND AEROSOL TRACERS INTO GP REPRESENTATION
  ! NOTE: HERE, THE SEDIMENTATION VELOCITY IS AVERAGED, WHICH IS ALREADY
  !       BETTER THAN CALCULATING THE VELOCITY OF THE AVERAGE AEROSOL RADIUS,
  !       HOWEVER STILL NOT PERFECT ...
  CALL lg2gp_e5(vt_lg, vt_gp, LG2GP_AVE)
  CALL lg2gp_e5(xt_aer_lg, xt_aer_gp, LG2GP_AVE, fill_value=0._dp)

  ! ADD RESIDUAL FLUX (IN GP)
  DO jjrow=1, ngpblks
     IF (jjrow == ngpblks) THEN
        kproma = npromz
     ELSE
        kproma = nproma
     END IF
     DO jk=1,nlev
        DO jp=1,kproma
              xt_aer_gp(jp,jk,:,jjrow) = xt_aer_gp(jp,jk,:,jjrow) + &
                   fx_residual(jp,jk,:,jjrow,1) &
                   / (delp(jp,jk,jjrow) / g) * ztmst
        END DO
     END DO
  END DO

  ! CALCULATE SEDIMENTATION TENDENCY AND LOSS IN GP REPRESENTATION
  DO jjrow=1, ngpblks
     IF (jjrow == ngpblks) THEN
        kproma = npromz
     ELSE
        kproma = nproma
     END IF
     DO jt=1, na

        ! CALCULATE LOSS-TENDENCY AND OUTFLUX
        SELECT CASE  (order) 
        CASE(0)
           CALL sedi_0( kproma, nlev                     &
                , tend_gp(1:kproma,1:nlev,jt,jjrow)      &
                , xt_aer_gp(1:kproma,1:nlev,jt,jjrow)    &
                , vt_gp(1:kproma,1:nlev,jt,jjrow)        &
                , dz(1:kproma,1:nlev,jjrow)              &
                , delp(1:kproma,1:nlev,jjrow)            &
                , ztmst, sfloss(1:kproma,jt,jjrow) )
        CASE DEFAULT
           CALL error_bi('UNKOWN ORDER OF SEDIMENTATION SCHEME (LG)', substr)   
        END SELECT
                     
        ! SAVE LOSS (GP REPRESENTATION)
        ! TIME INTEGRAL OF LOSS AT SURFACE [kg/m^2]
        idx = SSET(I_LG)%idx_aer(jt)
        SELECT CASE(ti_lg(idx)%tp%ident%quantity)
        CASE(AMOUNTFRACTION)
           conv = conv_mr
           conv_sum = (ti_lg(idx)%tp%meta%cask_r(R_molarmass)/M_air)
        CASE(NUMBERDENSITY)
           conv = conv_nu
           conv_sum = 1.E-3_dp/M_air
        CASE DEFAULT
           conv = 0.0_dp
           conv_sum = 0.0_dp
        END SELECT

        SSET(I_LG)%sedflux(jt)%ptr(1:kproma,jjrow) =                          &
             sfloss(1:kproma,jt,jjrow) * conv                                 &
             * (delp(1:kproma,nlev,jjrow)/g)

        SSET(I_LG)%sedfluxsum(jt)%ptr(1:kproma,jjrow) = &
             SSET(I_LG)%sedfluxsum(jt)%ptr(1:kproma,jjrow) &
             + sfloss(1:kproma,jt,jjrow) * zdtime          &
             * conv_sum &
             * (delp(1:kproma,nlev,jjrow)/g)
     END DO
  END DO

  ! SUBTRACT RESIDUAL FLUX (IN GP) FOR CORRECT TENDENCY TRANSFORMATION
  ! (lggpte2lgte_e5) BELOW
  DO jjrow=1, ngpblks
     IF (jjrow == ngpblks) THEN
        kproma = npromz
     ELSE
        kproma = nproma
     END IF
     DO jk=1,nlev
        DO jp=1,kproma
              xt_aer_gp(jp,jk,:,jjrow) = xt_aer_gp(jp,jk,:,jjrow) - &
                   fx_residual(jp,jk,:,jjrow,1) &
                   / (delp(jp,jk,jjrow) / g) * ztmst
        END DO
     END DO
  END DO

  ! CALCULATE TOTAL TENDENCY = SEDI TENDENCY + RESIDUAL FLUX
  ! AND CONVERT TO TENDENCY (GP)
  ! AND UPDATE RESIDUAL FLUX
  DO jjrow=1, ngpblks
     IF (jjrow == ngpblks) THEN
        kproma = npromz
     ELSE
        kproma = nproma
     END IF
     DO jk=1,nlev
        DO jp=1,kproma
           tend_gp(jp,jk,:,jjrow) = tend_gp(jp,jk,:,jjrow) + &
                fx_residual(jp,jk,:,jjrow,1) &
                / (delp(jp,jk,jjrow) / g)
           IF (NINT(NCB(jp,jk,jjrow)) < 1) THEN
              fx_residual(jp,jk,:,jjrow,1) = tend_gp(jp,jk,:,jjrow) &
                   * (delp(jp,jk,jjrow) / g)
           ELSE
              fx_residual(jp,jk,:,jjrow,1) = 0.0_dp
           END IF
        END DO
     END DO
  END DO

  ! TRANSFORM GP-TENDENCY INTO LG-REPRESENTATION
  CALL lggpte2lgte_e5(xt_aer_gp, tend_gp, xt_aer_lg, tend_lg)

  ! UPDATE TENDENCY (LG REPRESENTATION)
  DO jt=1, na
     idx = SSET(I_LG)%idx_aer(jt)
     pxtte(:,idx) = pxtte(:,idx) + tend_lg(:,jt)
  END DO
  
  ! CLEAN UP
  DEALLOCATE(vt_gp)
  DEALLOCATE(xt_aer_gp)
  DEALLOCATE(sfloss)
  DEALLOCATE(tend_lg)
  DEALLOCATE(tend_gp)
  DEALLOCATE(vt_lg)
  DEALLOCATE(xt_aer_lg)
#endif
!!#D attila -

END SUBROUTINE sedi_global_end
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE sedi_free_memory

  IMPLICIT NONE
  
  INTEGER :: naermod,i,jt

  INTRINSIC ASSOCIATED, SIZE
  
  tracer_set_loop: DO i=1,2

     IF ((i==I_GP) .AND. (.NOT.L_GP)) CYCLE
     IF ((i==I_LG) .AND. (.NOT.L_LG)) CYCLE

     DO naermod=1, SSET(i)%aeromodnum
        IF (ASSOCIATED(SSET(i)%vsedi_mass(naermod)%ptr))  &
             DEALLOCATE(SSET(i)%vsedi_mass(naermod)%ptr)
        IF (ASSOCIATED(SSET(i)%vsedi_num(naermod)%ptr))  &
             DEALLOCATE(SSET(i)%vsedi_num(naermod)%ptr)
        NULLIFY(SSET(i)%radii(naermod)%ptr)
        NULLIFY(SSET(i)%densaer(naermod)%ptr)
        NULLIFY(SSET(i)%sigma(naermod)%ptr)
     END DO

     IF (ASSOCIATED(SSET(i)%radii))      DEALLOCATE(SSET(i)%radii)
     IF (ASSOCIATED(SSET(i)%densaer))    DEALLOCATE(SSET(i)%densaer)
     IF (ASSOCIATED(SSET(i)%sigma))      DEALLOCATE(SSET(i)%sigma)

     IF (ASSOCIATED(SSET(i)%sedflux)) THEN
        DO jt=1, SIZE(SSET(i)%sedflux)
           NULLIFY(SSET(i)%sedflux(jt)%ptr)
        END DO
        DEALLOCATE(SSET(i)%sedflux)
     END IF

     IF (ASSOCIATED(SSET(i)%sedfluxsum)) THEN
        DO jt=1, SIZE(SSET(i)%sedfluxsum)
           NULLIFY(SSET(i)%sedfluxsum(jt)%ptr)
        END DO
        DEALLOCATE(SSET(i)%sedfluxsum)
     END IF

     IF (ASSOCIATED(SSET(i)%vsedi_mass)) DEALLOCATE(SSET(i)%vsedi_mass)
     IF (ASSOCIATED(SSET(i)%vsedi_num))  DEALLOCATE(SSET(i)%vsedi_num)

     IF (ASSOCIATED(SSET(i)%traermod))   DEALLOCATE(SSET(i)%traermod)
     IF (ASSOCIATED(SSET(i)%nmod))       DEALLOCATE(SSET(i)%nmod)
     IF (ASSOCIATED(SSET(i)%idt_aer))    DEALLOCATE(SSET(i)%idt_aer)
     IF (ASSOCIATED(SSET(i)%idx_aer))    DEALLOCATE(SSET(i)%idx_aer)

  END DO tracer_set_loop

  IF (ASSOCIATED(SSET))       DEALLOCATE(SSET)

  IF (ASSOCIATED(fx_residual)) DEALLOCATE(fx_residual)
  
END SUBROUTINE sedi_free_memory
!------------------------------------------------------------------------------

! =============================================================================
! PRIVATE ROUTINES
! =============================================================================

  SUBROUTINE sedi_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Andrea Pozzer, MPICH, Aug 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer_mem_bi, ONLY: NGCELL

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    ! switch for skipping calculation of Lagrangian 
    ! rate coefficients..it is local,not broadcasted

    NAMELIST /CPL/ L_GP, L_LG

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='sedi_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    IF ((L_LG) .AND. (NGCELL > 0)) THEN
       L_LG = .TRUE.
    ELSE
       IF (L_LG) THEN
         WRITE(*,*) 'L_LG = T in namelist'
         WRITE(*,*) 'However no Lagrangian scheme activated ...'
         WRITE(*,*) ' ... setting L_LG = F'
       END IF
       L_LG = .FALSE.
    ENDIF

    IF (L_GP) THEN
       WRITE(*,*) 'SEDIMENTATION IN GRIDPOINT SPACE : ON'
    ELSE
       WRITE(*,*) 'SEDIMENTATION IN GRIDPOINT SPACE : OFF'
    END IF

    IF (L_LG) THEN
       WRITE(*,*) 'SEDIMENTATION IN LAGRANGIAN SPACE : ON'
    ELSE
       WRITE(*,*) 'SEDIMENTATION IN LAGRANGIAN SPACE : OFF'
    END IF

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

    IF (.NOT.L_GP .AND. .NOT.L_LG)  &
       CALL error_bi('  NEITHER GRID POINT NOR LAGRANGIAN REPRESENTATION'//&
       &' REQUESTED (namelist CPL): No calculation possible!', substr )


  END SUBROUTINE sedi_read_nml_cpl

! ===========================================================================

END MODULE messy_sedi_si
