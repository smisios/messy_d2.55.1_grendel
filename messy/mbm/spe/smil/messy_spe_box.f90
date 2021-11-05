! ***********************************************************************
MODULE messy_spe_box
! ***********************************************************************

  ! MESSy-SMIL FOR SUBMODEL SPE
  !
  ! SPE = SOLAR PROTON EVENT PARAMETERIZATION
  !
  ! Author: Andreas Baumgaertner, MPICH, 2010


  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , info_bi, error_bi
  ! MESSy
  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL &
                                    , t_chaobj_cpl
  USE messy_spe

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! POINTERS TO DIAGNOSTIC CHANNEL OBJECTS:
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: ions     ! SPE ion pair production [#/cm3/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: xnox     ! SPE NOx production [g/cm3]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: xhox     ! SPE HOx production [g/cm3]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: xhno3    ! SPE HNO3 production [g/cm3]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespen   ! SPE N  prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespeno  ! SPE NO prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespeh   ! SPE H  prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespeoh  ! SPE OH prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespehno3  ! SPE OH prod. tend. [mol/mol/s]

  ! CPL NAMELIST
  TYPE(t_chaobj_cpl) :: spe_data_int, spe_data_ion

  ! GLOBAL PARMETERS
  ! TRACERS 
  INTEGER                            :: nlntrac_n    ! no of N  tracers
  INTEGER                            :: nlntrac_no   ! no of NO tracers
  INTEGER                            :: nlntrac_h    ! no of H  tracers
  INTEGER                            :: nlntrac_oh   ! no of OH tracers
  INTEGER                            :: nlntrac_hno3 ! no of HNO3 tracers
  INTEGER, DIMENSION(:), POINTER     :: idt_list_n
  INTEGER, DIMENSION(:), POINTER     :: idt_list_no
  INTEGER, DIMENSION(:), POINTER     :: idt_list_h
  INTEGER, DIMENSION(:), POINTER     :: idt_list_oh
  INTEGER, DIMENSION(:), POINTER     :: idt_list_hno3

  ! GLOBAL PARMETERS
  CHARACTER(LEN=STRLEN_CHANNEL), PUBLIC :: AIMOS_channel = 'offlem'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: AIMOS_p_object = 'jmionrates_ionrate_p'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: AIMOS_e_object = 'jmionrates_ionrate_e'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: AIMOS_a_object = 'jmionrates_ionrate_a'

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: spe_initialize
  PUBLIC :: spe_init_memory
  PUBLIC :: spe_init_coupling
  PUBLIC :: spe_global_start
  PUBLIC :: spe_physc
  PUBLIC :: spe_free_memory
  !PRIVATE :: spe_read_nml_cpl
  !PRIVATE :: spe_read_data

CONTAINS

! ========================================================================
  SUBROUTINE spe_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL spe_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CTRL NAMELIST', substr)
    END IF

    CALL p_bcast(spe_method,   p_io)
    CALL p_bcast(r_lat1,       p_io)
    CALL p_bcast(r_lat2,       p_io)
    CALL p_bcast(nbin,         p_io)
    CALL p_bcast(minoffs,      p_io)
    CALL p_bcast(spect_interp_method,p_io)
    CALL p_bcast(rangerela,p_io)
    CALL p_bcast(rangerelb,p_io)
    CALL p_bcast(ion_km,p_io)
    CALL p_bcast(Nperion_km,p_io)
    CALL p_bcast(NOperion_km,p_io)

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL spe_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CPL NAMELIST', substr)
    END IF
    CALL p_bcast(spe_data_int%cha, p_io)
    CALL p_bcast(spe_data_int%obj, p_io)
    CALL p_bcast(spe_data_ion%cha, p_io)
    CALL p_bcast(spe_data_ion%obj, p_io)

 
  END SUBROUTINE spe_initialize
! ========================================================================

! ========================================================================
  SUBROUTINE spe_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_init_memory'
    INTEGER :: status

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ions', p3=ions)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'ions' &
         , 'long_name', c='SPE ion pair production')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ions' &
         , 'units', c='ions/cm3/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'xnox', p3=xnox)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'xnox' &
         , 'long_name', c='SPE NOx production')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xnox' &
         , 'units', c='g/cm3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'xhox', p3=xhox)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'xhox' &
         , 'long_name', c='SPE HOx production')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xhox' &
         , 'units', c='g/cm3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'xhno3', p3=xhno3)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'xhno3' &
         , 'long_name', c='SPE HNO3 production')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xhno3' &
         , 'units', c='g/cm3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespen', p3=tespen)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespen' &
         , 'long_name', c='SPE N production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespen' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespeno', p3=tespeno)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespeno' &
         , 'long_name', c='SPE NO production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespeno' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespeh', p3=tespeh)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespeh' &
         , 'long_name', c='SPE NH production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespeh' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespeoh', p3=tespeoh)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespeoh' &
         , 'long_name', c='SPE OH production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespeoh' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespehno3', p3=tespehno3)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespehno3' &
         , 'long_name', c='SPE HNO3 production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespehno3' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY   
       ALLOCATE(pflxspc(nbin))
       ALLOCATE(energies(nbin))
       ALLOCATE(rangestp(nbin))
    CASE(1)
       ! nothing for now
    CASE(2)
       ! nothing for now   
    END SELECT

  END SUBROUTINE spe_init_memory
! ========================================================================

! ========================================================================
  SUBROUTINE spe_init_coupling

    ! ECHAM5/MESSy
!    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
!    USE messy_main_tracer_bi,     ONLY: tracer_halt
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_channel_error_bi,    ONLY: channel_halt
    ! MESSy
    USE messy_main_channel,       ONLY: get_channel_object &
                                      , get_channel_object_dimvar
!    USE messy_main_tracer,        ONLY: get_tracer_list
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_ULONG
    USE messy_main_tools,         ONLY: PTR_1D_ARRAY

    IMPLICIT NONE
    INTRINSIC :: SIZE, NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_init_coupling'
    INTEGER :: status
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: subnames => NULL()
    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER           :: dvs => NULL()
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER  :: units => NULL()
    INTEGER :: jt

!!$    CALL start_message_bi(modstr, 'LOOKING FOR TRACERS', substr)
!!$
!!$    CALL get_tracer_list(status, GPTRSTR, 'N', idt_list_n, subnames)
!!$    CALL tracer_halt(substr, status)
!!$    nlntrac_n = SIZE(idt_list_n)
!!$    IF (p_parallel_io) THEN
!!$       DO jt=1, nlntrac_n
!!$          IF (TRIM(subnames(jt)) == '') THEN
!!$             WRITE(*,*) ' ... N'
!!$          ELSE
!!$             WRITE(*,*) ' ... N_'//TRIM(subnames(jt))
!!$          END IF
!!$       END DO
!!$    END IF
!!$
!!$    CALL get_tracer_list(status, GPTRSTR, 'NO', idt_list_no, subnames)
!!$    CALL tracer_halt(substr, status)
!!$    nlntrac_no = SIZE(idt_list_no)
!!$    IF (p_parallel_io) THEN
!!$       DO jt=1, nlntrac_no
!!$          IF (TRIM(subnames(jt)) == '') THEN
!!$             WRITE(*,*) ' ... NO'
!!$          ELSE
!!$             WRITE(*,*) ' ... NO_'//TRIM(subnames(jt))
!!$          END IF
!!$       END DO
!!$    END IF
!!$
!!$    CALL get_tracer_list(status, GPTRSTR, 'OH', idt_list_oh, subnames)
!!$    CALL tracer_halt(substr, status)
!!$    nlntrac_oh = SIZE(idt_list_oh)
!!$    IF (p_parallel_io) THEN
!!$       DO jt=1, nlntrac_oh
!!$          IF (TRIM(subnames(jt)) == '') THEN
!!$             WRITE(*,*) ' ... OH'
!!$          ELSE
!!$             WRITE(*,*) ' ... OH_'//TRIM(subnames(jt))
!!$          END IF
!!$       END DO
!!$    END IF
!!$
!!$    CALL get_tracer_list(status, GPTRSTR, 'H', idt_list_h, subnames)
!!$    CALL tracer_halt(substr, status)
!!$    nlntrac_h = SIZE(idt_list_h)
!!$    IF (p_parallel_io) THEN
!!$       DO jt=1, nlntrac_h
!!$          IF (TRIM(subnames(jt)) == '') THEN
!!$             WRITE(*,*) ' ... H'
!!$          ELSE
!!$             WRITE(*,*) ' ... H_'//TRIM(subnames(jt))
!!$          END IF
!!$       END DO
!!$    END IF
!!$
!!$    IF (ASSOCIATED(subnames)) DEALLOCATE(subnames)
!!$    NULLIFY(subnames)
!!$
!!$    CALL end_message_bi(modstr, 'LOOKING FOR TRACERS', substr)

    CALL start_message_bi(modstr, 'LOOKING FOR REQUIRED CHANNEL OBJECTS', substr)

    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY
       CALL get_channel_object(status &
            , TRIM(spe_data_int%cha), TRIM(spe_data_int%obj), p1=pflx)
       CALL channel_halt(substr, status)
       npch = SIZE(pflx)

       CALL get_channel_object_dimvar(status &
            , TRIM(spe_data_int%cha), TRIM(spe_data_int%obj) &
            , dvs, units)
       CALL channel_halt(substr, status)
       chrig => dvs(1)%ptr

    CASE(1) ! EXTERNAL IONIZATION RATES
       CALL get_channel_object(status &
            , TRIM(spe_data_ion%cha), TRIM(spe_data_ion%obj), p1=ions_ext)
       CALL channel_halt(substr, status)        
       naltitudes = SIZE(ions_ext)

    CASE(2) ! EXTERNAL IONIZATION RATES FROM AIMOS / JAN MAIK WISSING 
       CALL get_channel_object(status, TRIM(AIMOS_channel), TRIM(AIMOS_p_object) &
            , p3=AIMOS_p)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, TRIM(AIMOS_channel), TRIM(AIMOS_e_object) &
            , p3=AIMOS_e)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, TRIM(AIMOS_channel), TRIM(AIMOS_a_object) &
            , p3=AIMOS_a)
       CALL channel_halt(substr, status)
    END SELECT

    CALL end_message_bi(modstr, 'LOOKING FOR REQUIRED CHANNEL OBJECTS', substr)

  END SUBROUTINE spe_init_coupling
! ========================================================================

! ========================================================================
  SUBROUTINE spe_global_start

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_global_start'
    INTEGER            :: status

    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY
       CALL spe_interp
    CASE(1) ! EXTERNAL IONIZATION RATES
       ! nothing to be done currently
    CASE(2) ! EXTERNAL IONIZATION RATES FROM AIMOS
       ! nothing to be done currently
    END SELECT

  END SUBROUTINE spe_global_start
! ========================================================================

! ========================================================================
  SUBROUTINE spe_physc

    ! ECHAM5/MESSy
    USE messy_main_constants_mem, ONLY: M_air, R_gas, g
    USE messy_main_data_bi,       ONLY: pint => pressi_3d & ! level interface (above) pressure [Pa]
         , pmid => press_3d  & ! mid-level pressures [Pa]
         , t_scb             &
         , geopot_3d         & ! mz_ab_20080620
         , temp => t_scb     &
         , jrow, kproma      &
         , nlev              &
         , philat, philon    &
         , ilat, ilon        
    USE messy_main_timer,         ONLY: time_step_len    

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_physc'
    INTEGER                     :: status   ! error status
    REAL(DP)                    :: zpb, zpt
    INTEGER                     :: jk, jp, jt
    REAL(DP)                    :: scale_height_ma = 7._dp ! Middle atmosphere scale height [km]

    ! LOCAL FIELDS
    REAL(dp),  DIMENSION(:,:), POINTER    :: grheight ! height of box [m]
    REAL(dp),  DIMENSION(:,:), POINTER    :: density  ! air density at level midpoints [g cm-3]
    REAL(dp),  DIMENSION(:,:), POINTER    :: densityint  ! air density at level interfaces [g cm-3]
    REAL(dp),  DIMENSION(:,:), POINTER    :: altitude  ! altitude [km]

    INTRINSIC LOG, MAX, REAL

    ! ALLOCATE MEMORY
    ALLOCATE(grheight(kproma,nlev))
    ALLOCATE(density(kproma,nlev))
    ALLOCATE(densityint(kproma,nlev))
    ALLOCATE(altitude(kproma,nlev))

    ! CALCULATE ALTITUDE in km
!    altitude(1:kproma,:)=-scale_height_ma*LOG(pmid(1:kproma,:,jrow)/1E5_dp)

    altitude(1:kproma,:)=geopot_3d(1:kproma,:,jrow)/g/1E3_dp
!!!write(*,*) altitude(1,1:20)

    ! CALCULATE DENSITY, convert from g/m3 to g/cm3
    density(1:kproma,:)=pmid(1:kproma,:,jrow)*M_air/(temp(1:kproma,:,jrow)*R_gas)*1e-6_dp
    densityint(1:kproma,:)=pint(1:kproma,1:nlev,jrow)*M_air/(temp(1:kproma,:,jrow)*R_gas)*1e-6_dp

    ! HEIGHT OF GRID-BOX [m]
    DO jk=1, nlev
       DO jp=1, kproma
          zpb = pint(jp,jk+1, jrow)
          ! ECHAM5 top layer ends at 0. Pa !!! adjust to 1. Pa
          zpt = MAX(pint(jp, jk, jrow),1._dp)
          grheight(jp,jk) = (1000._dp * R_gas / (M_air * g)) &
               * temp(jp,jk,jrow) * log(zpb/zpt)
       END DO
    END DO

    SELECT CASE (spe_method) 
! op_pj_20170220: not longer available after update
!!$    CASE(0) ! CALCULATE IONRATES INTERNALLY
!!$       CALL SPE_PROD (grheight,density,densityint      &  ! INPUT
!!$            ,altitude                                  &  ! INPUT
!!$            ,nlev, time_step_len                       &  ! INPUT 
!!$            ,philat, philon                            &  ! INPUT 
!!$            ,ilat, ilon                                &  ! INPUT 
!!$            ,jrow, kproma                              &  ! INPUT 
!!$            ,ions                                      &  ! OUTPUT -> stream
!!$            ,xnox, tespen, tespeno                     &  ! OUTPUT -> stream
!!$            ,xhox, tespeh, tespeoh                     &  ! OUTPUT -> stream
!!$            ,xhno3, tespehno3                          &  ! OUTPUT -> stream
!!$            ,status                                    &  ! OUTPUT ERR. STATUS
!!$            )

! op_pj_20170220: update required; fiels marked with '<!qqq' missing ...
!!$       CALL SPE_IONS(nlev,time_step_len                 & ! INPUT 
!!$            ,grheight                                   & ! INPUT
!!$            ,density,densityint                         & ! INPUT
!!$            ,altitude                                   & ! INPUT
!!$            ,ilat, ilon                                 & ! INPUT
!!$            ,philat, philon                             & ! INPUT
!!$            ,jrow, kproma                               & ! INPUT
!!$            ,spe_method                                 & ! INPUT      <!qqq
!!$            ,grmass,grvol                               & ! INPUT      <!qqq
!!$            ,hn_ion_lat_bin                             & ! INPUT      <!qqq
!!$            ,ikp,nkp,nlat,nlev_a                        & ! INPUT      <!qqq
!!$            ,pmid                                       & ! INPUT
!!$            ,pressure                                   & ! INPUT      <!qqq
!!$            ,ionrate                                    & ! INPUT      <!qqq
!!$            ,ions                                       & ! OUTPUT 
!!$            ,maglat                                     & ! OUTPUT     <!qqq
!!$            ,status                                     & ! OUTPUT
!!$            )

    CASE DEFAULT
       CALL start_message_bi(modstr, 'This spe_method is not supported in mbm', substr)
       STOP
    END SELECT



    ! DEALLOCATE MEMORY
    DEALLOCATE(grheight)
    DEALLOCATE(density)
    DEALLOCATE(densityint)
    DEALLOCATE(altitude)

  END SUBROUTINE spe_physc
! ========================================================================

! ========================================================================
  SUBROUTINE spe_free_memory

    IMPLICIT NONE
    INTRINSIC :: ALLOCATED, ASSOCIATED

!!$    IF (ASSOCIATED(idt_list_n))  DEALLOCATE(idt_list_n)
!!$    IF (ASSOCIATED(idt_list_no)) DEALLOCATE(idt_list_no)
!!$    IF (ASSOCIATED(idt_list_h))  DEALLOCATE(idt_list_h)
!!$    IF (ASSOCIATED(idt_list_oh)) DEALLOCATE(idt_list_oh)

  END SUBROUTINE spe_free_memory
! ========================================================================
  SUBROUTINE spe_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ spe_data_ion, spe_data_int &
         , AIMOS_channel, AIMOS_p_object &
         , AIMOS_e_object, AIMOS_a_object

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='spe_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY
       IF (TRIM(spe_data_int%cha) == '') THEN
          CALL info_bi('ERROR: empty channel name for proton flux data')
          RETURN
       ELSE
          CALL info_bi('proton flux channel :'//spe_data_int%cha)
       END IF
       IF (TRIM(spe_data_int%obj) == '') THEN
          CALL info_bi('ERROR: empty channel object name for proton flux data')
          RETURN
       ELSE
          CALL info_bi('proton flux object  :'//spe_data_int%obj)
       END IF
    CASE(1) ! EXTERNAL IONIZATION RATES
       IF (TRIM(spe_data_ion%cha) == '') THEN
          CALL info_bi('ERROR: empty channel name for ionization rates')
          RETURN
       ELSE
          CALL info_bi('ionization rates channel :'//spe_data_ion%cha)
       END IF
       IF (TRIM(spe_data_ion%obj) == '') THEN
          CALL info_bi('ERROR: empty channel object name for ionization rates')
          RETURN
       ELSE
          CALL info_bi('ionization rates object  :'//spe_data_ion%obj)
       END IF
    CASE(2) ! EXTERNAL IONIZATION RATES FROM AIMOS / JAN MAIK WISSING
       WRITE(*,*) ' Proton channel object: ',  AIMOS_channel, AIMOS_p_object
       WRITE(*,*) ' Electron channel object: ', AIMOS_channel, AIMOS_e_object
       WRITE(*,*) ' Alpha particle channel object: ', AIMOS_channel, AIMOS_a_object
    END SELECT

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE spe_read_nml_cpl
! ========================================================================

! ***********************************************************************
END MODULE messy_spe_box
! ***********************************************************************

