! ***********************************************************************
MODULE messy_spacenox_e5
! ***********************************************************************

  ! MESSy-SMIL FOR SUBMODEL SPACENOX
  !
  ! Author: Andreas Baumgaertner, MPICH, September 2007

  ! TODO:

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      info_bi, error_bi, warning_bi ! mz_ab_20101124
  ! MESSy
  USE messy_main_channel,       ONLY: t_chaobj_cpl ! mz_ab_20101124
  USE messy_main_constants_mem, ONLY: DP
  USE messy_spacenox

  IMPLICIT NONE
  SAVE

  PRIVATE

  ! POINTERS TO DIAGNOSTIC STREAM ELEMENTS (NO ALLOCATE/DEALLOCATE NECESSARY):
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: teeppie_no  ! NO prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tegcr_n     ! N prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tegcr_no    ! NO prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tegcr_oh    ! OH prod. tend. [mol/mol/s]

  INTEGER :: num_months=0

  ! CPL namelist
  ! - name of tracer/channel object
  TYPE(t_chaobj_cpl) :: spacenox_Ap
  
  ! GLOBAL PARMETERS
  ! TRACERS 
  INTEGER                            :: nlntrac_n    ! no of N tracers
  INTEGER                            :: nlntrac_no   ! no of NO tracers
  INTEGER                            :: nlntrac_oh   ! no of OH tracers
  INTEGER, DIMENSION(:), POINTER     :: idt_list_n
  INTEGER, DIMENSION(:), POINTER     :: idt_list_no
  INTEGER, DIMENSION(:), POINTER     :: idt_list_oh

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: spacenox_initialize
  PUBLIC :: spacenox_init_memory
  PUBLIC :: spacenox_init_coupling
  PUBLIC :: spacenox_physc
  PUBLIC :: spacenox_free_memory
  !PRIVATE :: spacenox_read_nml_cpl

CONTAINS

! ========================================================================
  SUBROUTINE spacenox_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spacenox_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL spacenox_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CTRL NAMELIST', substr)
    END IF

    CALL p_bcast(EPPIE_latN,         p_io)
    CALL p_bcast(EPPIE_latS,         p_io)
    CALL p_bcast(EPPIE_coeffN,       p_io)
    CALL p_bcast(EPPIE_coeffS,       p_io)

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL spacenox_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CPL NAMELIST', substr)
    END IF
  
   CALL p_bcast(spacenox_Ap%cha, p_io)
   CALL p_bcast(spacenox_Ap%obj, p_io)

  END SUBROUTINE spacenox_initialize
! ========================================================================

! ========================================================================
  SUBROUTINE spacenox_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spacenox_init_memory'
    INTEGER :: status

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)


    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'teeppie_no', p3=teeppie_no)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'teeppie_no' &
         , 'long_name', c='SPACENOX EPP IE NO production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'teeppie_no' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tegcr_n', p3=tegcr_n)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tegcr_n' &
         , 'long_name', c='SPACENOX GCR N production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tegcr_n' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tegcr_no', p3=tegcr_no)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tegcr_no' &
         , 'long_name', c='SPACENOX GCR NO production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tegcr_no' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tegcr_oh', p3=tegcr_oh)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tegcr_oh' &
         , 'long_name', c='SPACENOX GCR OH production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tegcr_oh' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE spacenox_init_memory
! ========================================================================

! ========================================================================
  SUBROUTINE spacenox_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: get_tracer_list
    USE messy_main_channel,       ONLY: get_channel_object
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE
    INTRINSIC :: SIZE, NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spacenox_init_coupling'
    INTEGER :: status
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: subnames => NULL()
    INTEGER :: jt

    CALL start_message_bi(modstr, 'INITIALIZE COUPLING', substr)

    CALL get_channel_object(status &
         , TRIM(spacenox_Ap%cha), TRIM(spacenox_Ap%obj), p1=Ap_data)
    CALL channel_halt(substr, status)

    CALL get_tracer_list(status, GPTRSTR, 'N', idt_list_n, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_n = SIZE(idt_list_n)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_n
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... N'
          ELSE
             WRITE(*,*) ' ... N_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF

    CALL get_tracer_list(status, GPTRSTR, 'NO', idt_list_no, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_no = SIZE(idt_list_no)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_no
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... NO'
          ELSE
             WRITE(*,*) ' ... NO_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF

    CALL get_tracer_list(status, GPTRSTR, 'OH', idt_list_oh, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_oh = SIZE(idt_list_oh)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_oh
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... OH'
          ELSE
             WRITE(*,*) ' ... OH_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF


    IF (ASSOCIATED(subnames)) DEALLOCATE(subnames)
    NULLIFY(subnames)

    CALL end_message_bi(modstr, 'INITIALIZE COUPLING', substr)

  END SUBROUTINE spacenox_init_coupling
! ========================================================================

! ========================================================================
  SUBROUTINE spacenox_physc

    ! ECHAM5/MESSy
    USE messy_main_constants_mem, ONLY: pi, M_air   &
                                       ,R_gas   &
                                       ,g       &
                                       ,k_B ! Boltzmann constant [J/K]
    USE messy_main_data_bi,       ONLY: pint => pressi_3d & ! level interface (above) pressure [Pa]
                                      , pmid => press_3d  & ! mid-level pressures [Pa]
                                      , t_scb
    USE messy_main_grid_def_mem_bi, ONLY:jrow, kproma      &
                                      , nlev
    USE messy_main_grid_def_bi,     ONLY:philat, philon    &
! um_ak_20091002                                      , dayofyear         &
                                      , ilat, ilon
    USE messy_main_timer,         ONLY:  year, DAYOFYEAR ! um_ak_20091002  
    USE messy_main_tracer_mem_bi, ONLY: pxtte=>qxtte

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'spacenox_physc'
    INTEGER                     :: status   ! error status
    REAL(DP)                    :: zpb, zpt
    INTEGER                     :: jk, jp, jt
    REAL(DP)                    :: scale_height_ma = 7._dp ! Middle atmosphere scale height [km]
    ! LOCAL FIELDS
    REAL(dp),  DIMENSION(:,:), POINTER    :: grheight ! height of box [m]
    REAL(dp),  DIMENSION(:,:), POINTER    :: density  ! air density at level midpoints [g cm-3]
    REAL(dp),  DIMENSION(:,:), POINTER    :: numdensity  ! total number density at level midpoints [# cm-3]
    REAL(dp),  DIMENSION(:,:), POINTER    :: altitude  ! altitude [km]

    INTRINSIC SIZE, LOG, MAX

    ! ALLOCATE MEMORY
    ALLOCATE(grheight(kproma,nlev))
    ALLOCATE(density(kproma,nlev))
    ALLOCATE(numdensity(kproma,nlev))
    ALLOCATE(altitude(kproma,nlev))

    ! CALCULATE ALTITUDE in km
    altitude(1:kproma,:)=-scale_height_ma*LOG(pmid(1:kproma,:,jrow)/1E5_dp)

    ! CALCULATE DENSITY, convert from g/m3 to g/cm3
    density(1:kproma,:)=pmid(1:kproma,:,jrow)*M_air/(t_scb(1:kproma,:,jrow)*R_gas)*1e-6_dp
    ! CALCULATE TOTAL NUMBER DENSITY, convert from #/m3 to #/cm3
    numdensity(1:kproma,:)=pmid(1:kproma,:,jrow)/(k_B*t_scb(1:kproma,:,jrow))*1e-6_dp
 
    ! HEIGHT OF GRID-BOX [m]
    DO jk=1, nlev
       DO jp=1, kproma
          zpb = pint(jp,jk+1, jrow)
          ! ECHAM5 top layer ends at 0. Pa !!! adjust to 1. Pa
          zpt = MAX(pint(jp, jk, jrow),1._dp)
          grheight(jp,jk) = (1000._dp * R_gas / (M_air * g)) &
               * t_scb(jp,jk,jrow) * log(zpb/zpt)
       END DO
    END DO 

    ! *******************************************************************
    ! EPP Indirect Effect
    ! *******************************************************************

    IF (pmid(1,1,1)<10)   THEN    ! only if top level above 0.1 hPa (10 Pa) 

       ! centred around solstice
       EPPIE_seasoncoeffN=max(0.1_dp,cos(pi/182.625*(dayofyear-355.25)))
       EPPIE_seasoncoeffS=max(0.1_dp,cos(pi/182.625*(dayofyear-172.625)))

       CALL SPACENOX_EPPIE (grheight,density           &  ! INPUT
            ,philat, ilat                              &  ! INPUT 
            ,jrow, kproma                              &  ! INPUT 
            ,teeppie_no                                &  ! OUTPUT -> streams
            ,status                                    &  ! OUTPUT ERR. STATUS
            )

       IF (nlntrac_no > 0) THEN
          ! ADD TO NO TENDENCY
          DO jt=1, nlntrac_no
             pxtte(1:kproma,:,idt_list_no(jt)) = pxtte(1:kproma,:,idt_list_no(jt)) + &
                  teeppie_no(1:kproma,:,jrow)
          END DO
       END IF

    END IF

    ! *******************************************************************
    ! Galactic Cosmic Rays
    ! *******************************************************************

    CALL SPACENOX_GCR (density,numdensity           &  ! INPUT
         ,philat,ilat                               &  ! INPUT 
         ,jrow, kproma, nlev                        &  ! INPUT 
         ,year                                      &  ! INPUT 
         ,tegcr_n, tegcr_no, tegcr_oh               &  ! OUTPUT -> stream
         ,status                                    &  ! OUTPUT ERR. STATUS
         )

    IF (nlntrac_n > 0) THEN
       ! ADD TO N TENDENCY
       DO jt=1, nlntrac_n
          pxtte(1:kproma,:,idt_list_n(jt)) = pxtte(1:kproma,:,idt_list_n(jt)) + &
               tegcr_n(1:kproma,:,jrow)
       END DO
    END IF
    IF (nlntrac_no > 0) THEN
       ! ADD TO NO TENDENCY
       DO jt=1, nlntrac_no
          pxtte(1:kproma,:,idt_list_no(jt)) = pxtte(1:kproma,:,idt_list_no(jt)) + &
               tegcr_no(1:kproma,:,jrow)
       END DO
    END IF
    IF (nlntrac_oh > 0) THEN
       ! ADD TO NO TENDENCY
       DO jt=1, nlntrac_oh
          pxtte(1:kproma,:,idt_list_oh(jt)) = pxtte(1:kproma,:,idt_list_oh(jt)) + &
               tegcr_oh(1:kproma,:,jrow)
       END DO
    END IF

    ! DEALLOCATE MEMORY
    DEALLOCATE(grheight)
    DEALLOCATE(density)
    DEALLOCATE(numdensity)
    DEALLOCATE(altitude)

  END SUBROUTINE spacenox_physc
! ========================================================================



! ========================================================================
  SUBROUTINE spacenox_free_memory

    IMPLICIT NONE
    INTRINSIC ASSOCIATED

    IF (ASSOCIATED(idt_list_n))   DEALLOCATE(idt_list_n)
    IF (ASSOCIATED(idt_list_no))  DEALLOCATE(idt_list_no)
    IF (ASSOCIATED(idt_list_oh))  DEALLOCATE(idt_list_oh)

 END SUBROUTINE spacenox_free_memory
! ========================================================================

! ************************************************************************
! PRIVATE ROUTINES
! ************************************************************************

! ========================================================================
  SUBROUTINE spacenox_read_nml_cpl(status, iou)
   
    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ spacenox_Ap

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='spacenox_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    IF (TRIM(spacenox_Ap%cha) == '') THEN
       CALL info_bi('ERROR: empty channel name for Ap index')
       RETURN
    ELSE
       CALL info_bi('Ap index channel :'//spacenox_Ap%cha)
    END IF
    IF (TRIM(spacenox_Ap%obj) == '') THEN
       CALL info_bi('ERROR: empty channel object name for Ap index')
       RETURN
    ELSE
       CALL info_bi('Ap index object  :'//spacenox_Ap%obj)
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR

  END SUBROUTINE spacenox_read_nml_cpl
! ========================================================================

! ***********************************************************************
END MODULE messy_spacenox_e5
! ***********************************************************************

