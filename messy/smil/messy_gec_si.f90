! ***********************************************************************
MODULE messy_gec_si
! ***********************************************************************

#if defined(ECHAM5) || defined(CESM1)

  ! MESSy-SMIL FOR SUBMODEL GEC (Global Electric Circuit)
  !
  ! Author: Andreas Baumgaertner, Boulder, 2014

  ! MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      info_bi, error_bi, warning_bi 
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_main_constants_mem, ONLY: dp
  USE messy_gec

  IMPLICIT NONE
  SAVE

  PRIVATE

  ! CPL namelist
  ! - name of tracer/channel object
  TYPE(t_chaobj_cpl) :: convectvar
  !mz_se_20170206+
  TYPE(t_chaobj_cpl) :: gcr ! coupling for the GCR modulation time series
  REAL(DP) :: max_gcr_mod   ! maximum of gcr_mod time series mz_se 20170206
  REAL(DP) :: min_gcr_mod   ! maximum of gcr_mod time series mz_se 20170206
  !mz_se_20170206-

  ! GLOBAL PARMETERS
  ! TRACERS 
  INTEGER                            :: nlntrac_Rn222    ! no of Rn222 tracers
  INTEGER                            :: idt_Rn222
  INTEGER, DIMENSION(:), POINTER     :: idt_list_Rn222

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: gec_initialize
  PUBLIC :: gec_init_memory
  PUBLIC :: gec_init_coupling
  PUBLIC :: gec_physc
  PUBLIC :: gec_free_memory
  !PRIVATE :: spacenox_read_nml_cpl

  ! CHANNEL OBJECTS
  REAL(DP), DIMENSION(:,:,:),   POINTER :: con             => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: con_nocld       => NULL()
  REAL(DP), DIMENSION(:,:),     POINTER :: Rcoltilde       => NULL()
  REAL(DP), DIMENSION(:,:),     POINTER :: Rcol_nocld      => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: ionrate_radon   => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: ionrate_direct  => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: ionrate_gcr     => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: ionrate_spe     => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: ionrate         => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: ioncon          => NULL()
  REAL(DP), DIMENSION(:,:),     POINTER :: current         => NULL()
  REAL(DP), DIMENSION(:,:),     POINTER :: currentcorr     => NULL()

  ! CHANNEL OBJECTS FROM COUPLING
  REAL(DP), DIMENSION(:,:,:),   POINTER :: zmmu     => NULL()
  REAL(DP), DIMENSION(:),       POINTER :: gcr_mod  => NULL() ! mz_se_20180530

CONTAINS

! ========================================================================
  SUBROUTINE gec_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gec_initialize'
    INTEGER                     :: k
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL gec_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CTRL NAMELIST', substr)
    END IF

    CALL p_bcast(sa,         p_io)
    ! mz_se_20170206+
    CALL p_bcast(ionradon_method, p_io)
    CALL p_bcast(iongcr_method, p_io)
    CALL p_bcast(ioncloud_method, p_io)
    CALL p_bcast(ionaerosol_method, p_io)
    CALL p_bcast(lssion, p_io)
    CALL p_bcast(ltotalipr, p_io)
    ! mz_se_20170206-

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL gec_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CPL NAMELIST', substr)
    END IF
  
    CALL p_bcast(convectvar%cha, p_io)
    CALL p_bcast(convectvar%obj, p_io)
    ! mz_se_20180530+
    CALL p_bcast(gcr%cha, p_io)
    CALL p_bcast(gcr%obj, p_io)
    CALL p_bcast(max_gcr_mod, p_io)
    CALL p_bcast(min_gcr_mod, p_io)
    ! mz_se_20180530-

  END SUBROUTINE gec_initialize
! ========================================================================

! ========================================================================
  SUBROUTINE gec_init_memory

    ! MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_2D_HORIZONTAL
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gec_init_memory'
    INTEGER :: status

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)


    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'con', p3=con)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'con' &
         , 'long_name', c='conductivity with clouds')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'con' &
         , 'units', c='S/m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'connocld', p3=con_nocld)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'connocld' &
         , 'long_name', c='connocldductivity without clouds')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'connocld' &
         , 'units', c='S/m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Rcoltilde', p2=Rcoltilde, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'Rcoltilde' &
         , 'long_name', c='Column resistance with clouds')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rcoltilde' &
         , 'units', c='Ohm*m2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Rcolnocld', p2=Rcol_nocld, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'Rcolnocld' &
         , 'long_name', c='Column resistance without clouds')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rcolnocld' &
         , 'units', c='Ohm*m2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ionrate_radon', p3=ionrate_radon)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ionrate_direct', p3=ionrate_direct)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ionrate_gcr', p3=ionrate_gcr)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ionrate_spe', p3=ionrate_spe)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ionrate', p3=ionrate)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ioncon', p3=ioncon)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'current', p2=current, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'currentcorr', p2=currentcorr, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)

!!$    CALL new_channel_object(status, modstr,'', p3=)
!!$    CALL channel_halt(substr, status)


    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE gec_init_memory
! ========================================================================

! ========================================================================
  SUBROUTINE gec_init_coupling

    ! MESSy
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer,           ONLY: get_tracer_list
    USE messy_main_channel,          ONLY: get_channel_object &
                                         , get_channel_object_dimvar ! mz_se_20180530
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM, STRLEN_ULONG
    USE messy_main_grid_def_bi,      ONLY: hyam, hybm
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
    USE messy_main_tools,            ONLY: PTR_1D_ARRAY


    IMPLICIT NONE
    INTRINSIC :: SIZE, NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gec_init_coupling'
    INTEGER :: status
    REAL(dp), DIMENSION(:), POINTER :: pref_mid
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: subnames => NULL()
    ! mz_se_20180530+
    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER          :: dvs => NULL()
    CHARACTER(len=STRLEN_ULONG), DIMENSION(:), POINTER :: units => NULL()
    ! mz_se_20180530-
    INTEGER :: jt

    CALL start_message_bi(modstr, 'INITIALIZE COUPLING', substr)

    CALL get_channel_object(status &
         , TRIM(convectvar%cha), TRIM(convectvar%obj), p3=zmmu)
    CALL channel_halt(substr, status)

    CALL get_tracer_list(status, GPTRSTR, 'Rn222', idt_list_Rn222, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_Rn222 = SIZE(idt_list_Rn222)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_Rn222
          IF (TRIM(subnames(jt)) == '') THEN
             idt_Rn222 = jt
             WRITE(*,*) ' ... Rn222'
          ELSE
             WRITE(*,*) ' ... Rn222_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF

    IF (ASSOCIATED(subnames)) DEALLOCATE(subnames)
    NULLIFY(subnames)

    IF (ionaerosol_method==1) then
       CALL get_channel_object(status, 'import_grid' &
            , 'GEC1_total_ion_aerosol_base', p3=carma_processed1a)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'import_grid' &
            , 'GEC1_total_ion_aerosol_volc', p3=carma_processed1b)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'import_grid' &
            , 'GEC2_trop_ion_aerosol_base', p3=carma_processed2)
       CALL channel_halt(substr, status)

    END IF

    ALLOCATE(pref_mid(nlev))
    pref_mid(:) = hyam(:) + hybm(:) * 101325._dp
    CALL GEC_init(pref_mid, nlev) 
    DEALLOCATE(pref_mid)

    ! mz_se_20170206+ read in time series for GCR variation
    IF (TRIM(gcr%cha) /= '') THEN
      CALL info_bi('Looking for GCR DATA ... ')
      CALL info_bi('       channel: '//gcr%cha)
      CALL info_bi('       object : '//gcr%obj)
     
      CALL get_channel_object(status &
          & , TRIM(gcr%cha), TRIM(gcr%obj), p1=gcr_mod)
      CALL channel_halt(substr, status)
  
      CALL get_channel_object_dimvar(status &
          & , TRIM(gcr%cha), TRIM(gcr%obj) &
          & , dvs, units)
      CALL channel_halt(substr, status)
  
      DO jt=1, SIZE(dvs)
         IF (p_parallel_io) &
           & WRITE(*,*) 'DIMVAR ',jt,' [',TRIM(units(jt)),']: ',dvs(jt)%ptr
      END DO
    END IF 
    ! mz_se_20170206-

    CALL end_message_bi(modstr, 'INITIALIZE COUPLING', substr)

  END SUBROUTINE gec_init_coupling
! ========================================================================

! ========================================================================
  SUBROUTINE gec_physc

    ! ECHAM5/MESSy
    USE messy_main_constants_mem, ONLY: pi, M_air   &
                                       ,R_gas   &
                                       ,g       &
                                       ,k_B ! Boltzmann constant [J/K]
    USE messy_main_data_bi,       ONLY: pint => pressi_3d & ! level interface (above) pressure [Pa]
                                      , pmid => press_3d  & ! mid-level pressures [Pa]
                                      , aclc              & ! large scale cloud cover 
                                      , tm1 => tm1_3d     &
                                      , landfrac=>slf       ! sea/land fraction
    USE messy_main_grid_def_mem_bi, ONLY: jrow, nproma, kproma  &
                                        , nlev, nlevp1
    USE messy_main_grid_def_bi,     ONLY: philat, philon              &
                                        , ilat, ilon, altitude_msl    &
                                        , altitude_gnd, altitudei_gnd
    USE messy_main_timer,         ONLY: year, DAYOFYEAR, time_step_len, hourUT=>hour
    USE messy_main_tracer_mem_bi, ONLY: pxtm1=>qxtm1

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gec_physc'
    INTEGER                     :: status   ! error status
    REAL(DP)                    :: zpb, zpt
    INTEGER                     :: jk, jp, jt, k
    REAL(DP)                    :: lsa  ! mz_se_20170206
    ! LOCAL FIELDS
    REAL(dp),  DIMENSION(:,:), POINTER    :: density    ! air density [g cm-3]
    REAL(dp),  DIMENSION(:,:), POINTER    :: numdensity ! total number density [# cm-3]
    REAL(dp),  DIMENSION(:,:), POINTER    :: mobility   ! [cm2/V/s]

    INTRINSIC SIZE, LOG, MAX

    ! ALLOCATE MEMORY
    ALLOCATE(density(kproma,nlev))
    ALLOCATE(numdensity(kproma,nlev))
    ALLOCATE(mobility(kproma,nlev))
 
    ! CALCULATE DENSITY, convert from g/m3 to g/cm3
    density(1:kproma,:)=pmid(1:kproma,:,jrow)*M_air/(tm1(1:kproma,:,jrow)*R_gas)*1e-6_dp
    ! CALCULATE TOTAL NUMBER DENSITY, convert from #/m3 to #/cm3
    numdensity(1:kproma,:)=pmid(1:kproma,:,jrow)/(k_B*tm1(1:kproma,:,jrow))*1e-6_dp
 
    ! mobility: paragraph 2.4 of TZ06, cm2 1/V 1/s
    mobility(:kproma,:) = &
         3.3_dp * (100000._dp / pmid(:kproma,:,jrow)) * (tm1(:kproma,:,jrow) / 273._dp)

    ! mz_se_20180530+
    IF (ASSOCIATED(gcr_mod)) THEN
      lsa = (gcr_mod(1) - min_gcr_mod) / max_gcr_mod
    ELSE
      lsa = sa
    END IF
    ! mz_se_20180530-

    CALL GEC_conductivity (pmid(:,:,jrow), pint(:,:,jrow) & ! INPUT
                         , tm1(:,:,jrow)              & ! INPUT
                         , density                    & ! INPUT
                         , numdensity                 & ! INPUT
                         , mobility                   & ! INPUT
                         , landfrac(:,jrow)           & ! INPUT
                         , altitude_gnd(:,:,jrow)     & ! INPUT
                         , altitudei_gnd(:,:,jrow)    & ! INPUT
                         , altitude_msl(:,:,jrow)     & ! INPUT
                         , pxtm1(:,:,idt_Rn222)       & ! INPUT
                         , aclc(:,:,jrow)             & ! INPUT
                         , philat, philon             & ! INPUT 
                         , ilat, ilon                 & ! INPUT 
                         , nproma, kproma, nlev, jrow & ! INPUT 
                         , time_step_len              & ! INPUT 
                         , lsa                        & ! INPUT mz_se 20170206
                         , con(:,:,jrow)              & ! OUTPUT
                         , con_nocld(:,:,jrow)        & ! OUTPUT
                         , Rcoltilde(:,jrow)          & ! OUTPUT
                         , Rcol_nocld(:,jrow)         & ! OUTPUT
                         , ionrate_radon(:,:,jrow)    & ! OUTPUT
                         , ionrate_direct(:,:,jrow)   & ! OUTPUT
                         , ionrate_gcr(:,:,jrow)      & ! OUTPUT
                         , ionrate_spe(:,:,jrow)      & ! OUTPUT
                         , ionrate(:,:,jrow)          & ! OUTPUT
                         , ioncon(:,:,jrow)           & ! OUTPUT
                         , status                     & ! OUTPUT ERR. STATUS
                         )
    IF (status /= 0) CALL error_bi('ERROR IN GEC_conductivity', substr)

    CALL GEC_sourcecurrent( landfrac(:,jrow)          & ! INPUT
                          , altitudei_gnd(:,:,jrow)   & ! INPUT
                          , zmmu(:,:,jrow)            & ! INPUT
                          , kproma                    & ! INPUT
                          , nlev                      & ! INPUT
                          , hourUT                    & ! INPUT
                          , current(:,jrow)           & ! OUTPUT
                          , currentcorr(:,jrow)       & ! OUTPUT
                          )

!!$    CALL gather_gp(current_gl,current)
!!$    sourcecurrent=sum(sum(current_gl,1),2)/size(current_gl ... weighted mean!
!!$    rRcol=1./Rcol_nocld
!!$    CALL gather_gp(rRcol_gl,rRcol)
!!$    rRtot=sum...
!!$    Rtot=1/rRtot
!!$    sourcecurrent=sum(sum(current_gl,1),2)
!!$
!!$    !PD=sourcecurrent*Rtot
!!$    
!!$    CALL GEC_solver(sourcecurrent,Rot,Rcol_nocld)
    
    ! DEALLOCATE MEMORY
    DEALLOCATE(density)
    DEALLOCATE(numdensity)
    DEALLOCATE(mobility)

  END SUBROUTINE gec_physc
! ========================================================================



! ========================================================================
  SUBROUTINE gec_free_memory

    IMPLICIT NONE
    INTRINSIC ASSOCIATED

    IF (ASSOCIATED(idt_list_Rn222))   DEALLOCATE(idt_list_Rn222)

 END SUBROUTINE gec_free_memory
! ========================================================================

! ************************************************************************
! PRIVATE ROUTINES
! ************************************************************************

! ========================================================================
  SUBROUTINE gec_read_nml_cpl(status, iou)
   
    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ convectvar &
         , gcr, max_gcr_mod, min_gcr_mod ! mz_se_20180530

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='gec_read_nml_cpl'
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
    IF (TRIM(convectvar%cha) == '') THEN
       CALL info_bi('ERROR: empty channel name for convectvar index')
       RETURN
    ELSE
       CALL info_bi('convectvar index channel :'//convectvar%cha)
    END IF
    IF (TRIM(convectvar%obj) == '') THEN
       CALL info_bi('ERROR: empty channel object name for convectvar index')
       RETURN
    ELSE
       CALL info_bi('convectvar index object  :'//convectvar%obj)
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR

  END SUBROUTINE gec_read_nml_cpl
! ========================================================================

#endif

! ***********************************************************************
END MODULE messy_gec_si
! ***********************************************************************

