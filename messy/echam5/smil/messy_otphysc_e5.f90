! **********************************************************************
!
!   OCEAN TRACER PHYSICS!
!
! Author : Andrea Pozzer, MPICH, September 2007
!
! This Submodel contains:
! - ADVECTION and CONVECTION of tracer in the ocean
! - DIFFUSION of tracer in the ocean
! - Update of tracer concentration based on the air-sea fluxes
!
! **********************************************************************
MODULE messy_otphysc_e5

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  ! MESSy
  USE messy_otphysc

  IMPLICIT NONE 
  PRIVATE

  LOGICAL :: LRUNSM = .TRUE.  ! RUN THIS SUBMODEL
  LOGICAL :: LRUNAOFLUX = .TRUE.  ! RUN A2O FLUXES
  LOGICAL :: LRUNRINPUT = .TRUE.  ! RUN RIVER INPUT

  ! AIRSEA/SCAV RELATED
  TYPE IO_AOFLUX
     CHARACTER(LEN=STRLEN_OBJECT)  :: name_tracer = ''
     CHARACTER(LEN=STRLEN_CHANNEL) :: aoflux_channel = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: aoflux_object = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: aoflux_process = ''
  END TYPE IO_AOFLUX

  TYPE T_AOFLUX
     TYPE(IO_AOFLUX) :: io
     REAL(DP), DIMENSION(:,:),POINTER :: aoflux  => NULL() 
     INTEGER                          :: trac_id_om 
     LOGICAL                          :: ok = .FALSE.
  END TYPE T_AOFLUX

  INTEGER, PARAMETER :: NMAXAOFLUX = 100
  TYPE(IO_AOFLUX), DIMENSION(NMAXAOFLUX), SAVE :: AOFLUX
  INTEGER,                                SAVE :: NAOFLUX = 0
  TYPE(T_AOFLUX),  DIMENSION(NMAXAOFLUX), SAVE :: XAOFLUX

  ! RIVER INPUT RELATED

  CHARACTER(LEN=STRLEN_OBJECT)  :: river_object = ''
  CHARACTER(LEN=STRLEN_CHANNEL) :: river_channel = ''
  REAL(DP), POINTER, DIMENSION (:,:) :: disch_m3s     

  TYPE IO_RINPUT
     CHARACTER(LEN=STRLEN_OBJECT)  :: name_tracer = ''
     REAL(DP)                      :: value = 0.0_dp 
  END TYPE IO_RINPUT

  ! mz_bk_20110322+
  TYPE IO_RINPUT_TR
     CHARACTER(LEN=STRLEN_OBJECT)  :: name_tracer = ''
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: object = ''
  END TYPE IO_RINPUT_TR
  ! mz_bk_20110322-

  TYPE T_RINPUT
     TYPE(IO_RINPUT)               :: io
     TYPE(IO_RINPUT_TR)            :: tr_io                    ! mz_bk_20110322
     INTEGER                       :: trac_id_om
     INTEGER                       :: input_type = 0           ! mz_bk_20110322
     LOGICAL                       :: ok = .FALSE.
  END TYPE T_RINPUT

  INTEGER, PARAMETER :: NMAXRINPUT = 100
  TYPE(IO_RINPUT), DIMENSION(NMAXRINPUT),    SAVE :: RINPUT
  TYPE(IO_RINPUT_TR), DIMENSION(NMAXRINPUT), SAVE :: RINPUT_TR ! mz_bk_20110322
  INTEGER,                                   SAVE :: NRINPUT = 0
  TYPE(T_RINPUT),  DIMENSION(NMAXRINPUT),    SAVE :: XRINPUT

  ! mz_bk_20110208+
  LOGICAL, PUBLIC :: L_SKIP_ADV, L_SKIP_MIX, L_SKIP_DIF &
       , L_SKIP_AOE, L_SKIP_RIV, L_SKIP_DIL
  ! mz_bk_20110208-

  ! TRACER DILUTION
  REAL(DP), POINTER, DIMENSION (:,:) :: l1     
  REAL(DP), POINTER, DIMENSION (:,:) :: l1_new 
  REAL(DP), POINTER, DIMENSION (:,:) :: zfac   

  PUBLIC :: otphysc_initialize
  PUBLIC :: otphysc_init_memory
  PUBLIC :: otphysc_init_coupling
  PUBLIC :: otphysc_global_start
  PUBLIC :: otphysc_global_end
  PUBLIC :: otphysc_free_memory

CONTAINS
  ! ---------------------------------------------------------------------

  SUBROUTINE otphysc_initialize

    ! otphysc MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Pozzer Andrea, MPICH, Aug 2007
    
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_channel_bi, ONLY: GP_3D_MPIOM, REPR_UNDEF
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'otphysc_initialize'
    INTEGER                 :: iou    ! I/O unit
    INTEGER                 :: status ! error status
    INTEGER                 :: i

    CALL start_message_bi(modstr, 'INITIALISATION', substr)


    LRUNSM = (GP_3D_MPIOM /= REPR_UNDEF)
    IF (.NOT. LRUNSM) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'MPIOM REPRESENTATION IS UNDEFINED.'
          WRITE(*,*) 'SUBMODEL CANNOT BE RUN.'
       END IF
       CALL end_message_bi(modstr, 'INITIALISATION', substr)
       CALL finish(substr)
    END IF

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL otphysc_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast(MASK_VALUE, p_io)   
   
    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL otphysc_read_nml_cpl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    IF (p_parallel_io) THEN

       ! GET NUMBER OF ENTRIES
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) ' FLUXES ATMOSPHERE <-> OCEAN:'
       WRITE(*,*) '------------------------------------------------------'
       ! COPY DATA
       NAOFLUX = 1
       DO i=1, NMAXAOFLUX

          IF (TRIM(AOFLUX(i)%name_tracer) == '') CYCLE
          XAOFLUX(NAOFLUX)%io%name_tracer    = TRIM(AOFLUX(i)%name_tracer)

          IF (TRIM(AOFLUX(i)%aoflux_channel) == '') CYCLE
          XAOFLUX(NAOFLUX)%io%aoflux_channel = TRIM(AOFLUX(i)%aoflux_channel)

          IF (TRIM(AOFLUX(i)%aoflux_object) == '') CYCLE
          XAOFLUX(NAOFLUX)%io%aoflux_object  = TRIM(AOFLUX(i)%aoflux_object)

          IF (TRIM(AOFLUX(i)%aoflux_process) == '') CYCLE
          XAOFLUX(NAOFLUX)%io%aoflux_process  = TRIM(AOFLUX(i)%aoflux_process)

          WRITE(*,*) '  ', TRIM(XAOFLUX(NAOFLUX)%io%name_tracer),' <- ', &
               TRIM(XAOFLUX(NAOFLUX)%io%aoflux_channel), &
               '(',TRIM(XAOFLUX(NAOFLUX)%io%aoflux_object),') ... process:',      &
               TRIM(XAOFLUX(NAOFLUX)%io%aoflux_process)  

          ! NEXT ENTRY
          NAOFLUX = NAOFLUX + 1
          WRITE(*,*) '------------------------------------------------------'

       END DO
       NAOFLUX = NAOFLUX - 1

       ! GET NUMBER OF ENTRIES
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) ' INPUT FROM RIVER:'
       WRITE(*,*) '------------------------------------------------------'
       ! COPY DATA
       NRINPUT = 1
       DO i=1, NMAXRINPUT

          IF (TRIM(RINPUT(i)%name_tracer) == '') CYCLE
          XRINPUT(NRINPUT)%io%name_tracer    = TRIM(RINPUT(i)%name_tracer)
          XRINPUT(NRINPUT)%io%value          = RINPUT(i)%value
          XRINPUT(NRINPUT)%input_type     = 1               ! mz_bk_20110322

          WRITE(*,*) '  ', TRIM(XRINPUT(NRINPUT)%io%name_tracer),' <- ', &
               ' ',XRINPUT(NRINPUT)%io%value,'mol/L '

          ! NEXT ENTRY
          NRINPUT = NRINPUT + 1
          WRITE(*,*) '------------------------------------------------------'

       END DO
       ! mz_bk_20110322+
       ! COPY DATA
       DO i=1, NMAXRINPUT

          IF (TRIM(RINPUT_TR(i)%name_tracer) == '') CYCLE
          XRINPUT(NRINPUT)%io%name_tracer    = TRIM(RINPUT_TR(i)%name_tracer)
          XRINPUT(NRINPUT)%input_type        = 2
          XRINPUT(NRINPUT)%tr_io%channel     = TRIM(RINPUT_TR(i)%channel)
          XRINPUT(NRINPUT)%tr_io%object      = TRIM(RINPUT_TR(i)%object)

          WRITE(*,*) '  ', TRIM(XRINPUT(NRINPUT)%io%name_tracer),' <- ',   &
               ' ',TRIM(XRINPUT(NRINPUT)%tr_io%channel),                   &
               ' ',TRIM(XRINPUT(NRINPUT)%tr_io%object)

          ! NEXT ENTRY
          NRINPUT = NRINPUT + 1
          WRITE(*,*) '------------------------------------------------------'

       END DO
       NRINPUT = NRINPUT - 1
       ! mz_bk_20110322-


    ENDIF

    ! BROADCAST RESULTS

    ! AIRSEA/SCAV INPUT
    CALL p_bcast(NAOFLUX, p_io)

    DO i=1, NAOFLUX
       CALL p_bcast(XAOFLUX(i)%io%name_tracer, p_io)
       CALL p_bcast(XAOFLUX(i)%io%aoflux_channel, p_io)
       CALL p_bcast(XAOFLUX(i)%io%aoflux_object, p_io)
       CALL p_bcast(XAOFLUX(i)%io%aoflux_process, p_io)
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NAOFLUX,' ATMOSPHERE-OCEAN FLUXE(S) INITIALIZED !'
    END IF

    ! RIVER INPUT
    CALL p_bcast(RIVER_OBJECT, p_io)
    CALL p_bcast(RIVER_CHANNEL, p_io)

    CALL p_bcast(NRINPUT, p_io)

    DO i=1, NRINPUT
       CALL p_bcast(XRINPUT(i)%io%name_tracer, p_io)
       CALL p_bcast(XRINPUT(i)%io%value, p_io)
       ! mz_bk_20110322+
       CALL p_bcast(XRINPUT(i)%input_type, p_io)
       CALL p_bcast(XRINPUT(i)%tr_io%channel, p_io)
       CALL p_bcast(XRINPUT(i)%tr_io%object, p_io)
       ! mz_bk_20110322-
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NRINPUT,' TRACER(S) INPUT FROM RIVERS INITIALIZED !'
    END IF

    ! mz_bk_20110208+
    CALL p_bcast(L_SKIP_ADV, p_io)
    CALL p_bcast(L_SKIP_MIX, p_io)
    CALL p_bcast(L_SKIP_DIF, p_io)
    CALL p_bcast(L_SKIP_DIL, p_io)
    CALL p_bcast(L_SKIP_AOE, p_io)
    CALL p_bcast(L_SKIP_RIV, p_io)
    ! mz_bk_20110208-

    LRUNAOFLUX = (NAOFLUX > 0)
    LRUNRINPUT = (NRINPUT > 0)

    CALL end_message_bi(modstr, 'INITIALISATION', substr)

  END SUBROUTINE otphysc_initialize

  ! ---------------------------------------------------------------------

  SUBROUTINE otphysc_init_memory

    ! otphysc MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Pozzer Andrea, MPICH, Aug 2007
    
    ! ECHAM5/MESSy
    USE messy_mpiom_mem_e5,    ONLY: IE,JE

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'otphysc_init_memory'

    CALL start_message_bi(modstr, 'MEMORY INITIALISATION', substr)

     ALLOCATE(l1(IE,JE))
     ALLOCATE(l1_new(IE,JE))
     ALLOCATE(zfac(IE,JE))

    CALL end_message_bi(modstr, 'MEMORY INITIALISATION', substr)

  END SUBROUTINE otphysc_init_memory
! ---------------------------------------------------------------------
  
  SUBROUTINE otphysc_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, message,finish
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_MPIOM, SCALAR
    ! MESSy
    USE messy_main_channel,         ONLY: new_channel, new_channel_object &
                                        , new_attribute, get_attribute    &
                                        , get_channel_object              &
                                        , get_channel_object_info         &
                                        , get_channel_info
    USE messy_main_constants_mem,   ONLY: STRLEN_ULONG
    USE messy_main_timer,           ONLY: time_step_len, lstart
    USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks, nlev
    USE messy_main_grid_def_bi,     ONLY: philon_2d, philat_2d          
    USE messy_main_tracer_mem_bi,   ONLY: ntrac_om, ti_om  

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    LOGICAL                     :: lfirst = .TRUE.
    CHARACTER(LEN=*), PARAMETER :: substr = 'otphysc_init_coupling'
    INTEGER                     :: status
    INTEGER                     :: i,ot
    INTEGER                     :: reprid
    CHARACTER(LEN=STRLEN_ULONG) :: unit

    IF (LRUNAOFLUX) THEN

       CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

       ! Atmosphere <-> Ocean
       IF (p_parallel_io) THEN
          WRITE(*,*) '------------------------------------------------------'
          WRITE(*,*) ' FLUX ATMOSPHERE <-> OCEAN:'
          WRITE(*,*) '------------------------------------------------------'
       END IF

       loop_aoflux: DO i=1, NAOFLUX

          IF (XAOFLUX(i)%ok) CYCLE

          CALL get_channel_object_info(status &
               , TRIM(XAOFLUX(i)%io%aoflux_channel)  &
               , TRIM(XAOFLUX(i)%io%aoflux_object)   &
               , reprid=reprid)

          IF (status /= 0) THEN
             CALL message('  ',TRIM(XAOFLUX(i)%io%aoflux_channel)//' - '//&
                  &TRIM(XAOFLUX(i)%io%aoflux_object)//' not found ... skipping' )
             CYCLE
          ELSE
             IF (reprid == GP_2D_MPIOM) THEN
                ! CHECK AVAILABILITY OF OCEAN CHANNEL/OBJECT
                CALL get_channel_object(status &
                     , TRIM(XAOFLUX(i)%io%aoflux_channel) &
                     , TRIM(XAOFLUX(i)%io%aoflux_object)  &
                     , p2=XAOFLUX(i)%aoflux &
                     )
                IF (status /= 0) THEN
                   CALL message('  ',TRIM(XAOFLUX(i)%io%aoflux_channel)//' - '//&
                        &TRIM(XAOFLUX(i)%io%aoflux_object)//' not found ... skipping' )
                   CYCLE
                ELSE
                   CALL message('  ',TRIM(XAOFLUX(i)%io%aoflux_channel)//' - '//&
                        &TRIM(XAOFLUX(i)%io%aoflux_object)//' found ...' )
                   CALL message('  ',' ... GP_2D_MPIOM representation ... ')
                ENDIF
             ELSE
                CALL message('  ',' ... wrong representation ... skipping')
                CYCLE
             ENDIF
          ENDIF

          XAOFLUX(i)%trac_id_om=0
          DO ot =1, ntrac_om
             IF(ti_om(ot)%tp%ident%basename == XAOFLUX(i)%io%name_tracer) XAOFLUX(i)%trac_id_om=ot  
          ENDDO
          IF (XAOFLUX(i)%trac_id_om==0) THEN
                CALL message('  ',' ...tracer not present between the ocean tracers...skipping')
                CYCLE
          ENDIF

          SELECT CASE(XAOFLUX(i)%io%aoflux_process)
          CASE('airsea')
            CALL message('  ',' process: airsea ')
          CASE('scav')
            CALL message('  ',' process: scav ')
          CASE DEFAULT
            CALL message('  ',' ... wrong process (airsea/scav) ... skipping')
            CYCLE
          END SELECT

          ! EVERYTHIG IS OK NOW
          XAOFLUX(i)%ok = .TRUE.

       END DO loop_aoflux

    ENDIF !LRUNAOFLUX

    IF (LRUNRINPUT) THEN

       CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

       IF (p_parallel_io) THEN
          WRITE(*,*) '------------------------------------------------------'
          WRITE(*,*) ' TRACER(S) INPUT FROM RIVER:'
          WRITE(*,*) '------------------------------------------------------'
          WRITE(*,*) ' RIVER CHANNEL : ', RIVER_CHANNEL  
          WRITE(*,*) ' RIVER OBJECT  : ', RIVER_OBJECT  
       END IF

       CALL get_channel_object(status &
            , river_channel &
            , river_object  &
            , p2=disch_m3s  &
            )
       IF (status /= 0) THEN
          CALL message('  ','river channel or object not available (check otphysc.nml)')
          !CALL FINISH(substr)
          CALL message('  ','No tracer input from river')
          NRINPUT = 0
       ENDIF

       loop_rinput: DO i=1, NRINPUT
       
          IF (XRINPUT(i)%ok) CYCLE

            IF (p_parallel_io) THEN
              WRITE(*,*) ' TRACER        : ', XRINPUT(i)%io%name_tracer  
              ! mz_bk_20110323+
              IF (XRINPUT(i)%input_type == 1) THEN
                 WRITE(*,*) ' INPUT(mol/L)  : ', XRINPUT(i)%io%value
              ELSE
                 WRITE(*,*) ' INPUT CHANNEL : ', XRINPUT(i)%tr_io%channel
                 WRITE(*,*) ' INPUT OBJECT  : ', XRINPUT(i)%tr_io%object
              ENDIF
              ! mz_bk_20110323-
            ENDIF
                         
          XRINPUT(i)%trac_id_om=0
          DO ot =1, ntrac_om
             IF(ti_om(ot)%tp%ident%basename == XRINPUT(i)%io%name_tracer) XRINPUT(i)%trac_id_om=ot  
          ENDDO
          IF (XRINPUT(i)%trac_id_om==0) THEN
                CALL message('  ',' ...tracer not present between the ocean tracers...skipping')
                CYCLE
          ENDIF

          ! EVERYTHIG IS OK NOW
          XRINPUT(i)%ok = .TRUE.

       END DO loop_rinput

    ENDIF !LRUNAOFLUX
    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE otphysc_init_coupling

! ------------------------------------------------------------------------

  SUBROUTINE otphysc_global_start

   USE messy_main_mpi_bi,        ONLY: message
   USE messy_main_tracer_mem_bi, ONLY: xt_om,ntrac_om
   USE messy_mpiom_mem_e5,       ONLY: dlxp, dlyp, dzw,zo,ddpo &
                                     , sictho, sicsno, rhoicwa, rhosnwa &
                                     , sicomo 
   USE messy_main_constants_mem, ONLY: N_A 
   USE messy_main_timer,         ONLY: delta_time, lstart 
   USE messy_mpiom_mem_e5,       ONLY: IE,JE,KE, WETO
   USE messy_mpiom_e5,           ONLY: l_trig_mpiom
   USE messy_main_channel,       ONLY: get_channel_object      ! mz_bk_20110323

   IMPLICIT NONE 

   REAL(dp) :: zflux_ao(IE,JE)  ! flux (kmol(trac) s-1)
   REAL(dp) :: bh(IE,JE) ! pure water box height

   INTEGER :: ot,I,J,K
   INTEGER :: status                                           ! mz_bk_20110323

   REAL(DP), POINTER, DIMENSION (:,:) :: river_input           ! mz_bk_20110323

  IF (.not.lstart) THEN
  bh(:,:)=ddpo(:,:,1)+ZO(:,:)-SICTHO(:,:)*RHOICWA &
                  -SICSNO(:,:)*RHOSNWA

    DO ot=1, NAOFLUX
     IF (XAOFLUX(ot)%ok) THEN
      ! mz_bk_20110208+
      IF (.not. L_SKIP_AOE) THEN
      ! mz_bk_20110208-
      IF (XAOFLUX(ot)%io%aoflux_process == 'airsea') THEN  
      !-------------------------------!
      !     AIRSEA GAS EXCHANGE       !
      !-------------------------------!
      ! The fluxes are in mol(trac) m-2 s-1
      !  mol(trac) m-2 s-1 -->  kmol m-3 s-1 
      ! => mcl(trac) m-2 s-1 / 1000 / bh 
      ! xt_om is in kmol m-3
      ! zflux_ao is negative (check airsea submodel):
      ! positive --> outgassing, negative --> entrainment in water
      DO I=1,IE
       DO J=1,JE
        IF (ddpo(I,J,1)==0) THEN
          zflux_ao(I,J) = 0
        ELSE
          zflux_ao(I,J) = XAOFLUX(ot)%aoflux(I,J)/bh(I,J)/1000
        ENDIF 
        xt_om(I,J,XAOFLUX(ot)%trac_id_om,1)=xt_om(I,J,XAOFLUX(ot)%trac_id_om,1) &
                                         -zflux_ao(I,J)*delta_time
       ENDDO
      ENDDO
      ELSEIF (XAOFLUX(ot)%io%aoflux_process == 'scav') THEN  
      !-------------------------------!
      !     TRACERS SCAVENGING        !
      !-------------------------------!
      ! The fluxes are in molecules(trac) m-2 s-1
      !  molecules(trac) m-2 s-1 -->  kmol m-3 s-1 
      ! => mcl(trac) m-2 s-1 /Navogadro / 1000 / bh 
      ! xt_om is in kmol m-3
      ! zflux_ao is always positive 
      DO I=1,IE
       DO J=1,JE
        IF (ddpo(I,J,1)==0) THEN
          zflux_ao(I,J) = 0
        ELSE
          zflux_ao(I,J) = XAOFLUX(ot)%aoflux(I,J)/N_A/bh(I,J)/1000
        ENDIF 
        xt_om(I,J,XAOFLUX(ot)%trac_id_om,1)=xt_om(I,J,XAOFLUX(ot)%trac_id_om,1) &
                                         +zflux_ao(I,J)*delta_time
       ENDDO
      ENDDO
      ENDIF
      ! mz_bk_20110208+
      ENDIF
      ! mz_bk_20110208-
     ENDIF

     IF (XRINPUT(ot)%ok) THEN
      ! mz_bk_20110208+
      IF (.not. L_SKIP_RIV) THEN
      ! mz_bk_20110208-
      !-------------------------------!
      !     TRACERS RIVER INPUT       !
      !-------------------------------!
         ! mz_bk_20110323
         IF (XRINPUT(ot)%input_type == 2) THEN
            CALL get_channel_object(status &
                 , XRINPUT(ot)%tr_io%channel &
                 , XRINPUT(ot)%tr_io%object  &
                 , p2=river_input  &
                 )
            IF (status /= 0) THEN
               CALL message('  ','nutrients input channel or object not available (check otphysc.nml)')
               !CALL FINISH(substr)
               CALL message('  ','No nutrients input from river')
               XRINPUT(ot)%input_type = 1
               XRINPUT(ot)%io%value = 0.0_dp
            ENDIF
         END IF
         DO I=1,IE
            DO J=1,JE
               ! mz_bk_20100204+
               IF (ddpo(I,J,1)/=0) THEN
                  ! mz_bk_20100204-
                  IF (XRINPUT(ot)%input_type == 1) THEN
                     xt_om(I,J,XRINPUT(ot)%trac_id_om,1) =          &
                          xt_om(I,J,XRINPUT(ot)%trac_id_om,1) +     &
                          disch_m3s(I,J)*                           & ! m3/s
                          XRINPUT(ot)%io%value*                     & ! kmol/m3
                          1./(dlxp(i,j)*dlyp(i,j)*bh(i,j))*         & ! 1/m3
                          delta_time                                  ! s
                     !NB: we need kmol/m3!!!
                  ELSEIF (XRINPUT(ot)%input_type == 2) THEN
                     xt_om(I,J,XRINPUT(ot)%trac_id_om,1) =          &
                          xt_om(I,J,XRINPUT(ot)%trac_id_om,1) +     &
                          disch_m3s(I,J)*                           & ! m3/s
                          river_input(I,J)*                         & ! kmol/m3
                          1./(dlxp(i,j)*dlyp(i,j)*bh(i,j))*         & ! 1/m3
                          delta_time                                  ! s
                     !NB: we need kmol/m3!!!
                  ENDIF
                  ! mz_bk_20100204+
               ENDIF
               ! mz_bk_20100204-
               
            ENDDO
         ENDDO
         ! mz_bk_20110208+
      ENDIF
      ! mz_bk_20110208-
         IF (XRINPUT(ot)%input_type == 2) NULLIFY(river_input)
     ENDIF


     
    ENDDO

    ! mz_bk_20110209+
#ifdef MPIOM_13B
    ! mz_bk_20110209-
    !-------------------------------!
    ! MASK OF CONCENTRATION ON LAND !
    !-------------------------------!
    !! Reset tracer values at dry cells to undefined.
    !! tracers have to be >= 0.
    DO ot=1, ntrac_om !ocean tracer loop
      !! Reset tracer values at dry cells to undefined.
      !! tracers have to be >= 0.
      DO K=1,KE
         DO J=1,JE
            DO I=1,IE
               IF(WETO(I,J,K).LT.0.5) THEN
                xt_om(I,J,ot,K)=MASK_VALUE
                ELSE      
                  ! to be moved to pdef
                  ! is not done there at the moment because 
                  ! these tracers do not have tendencies
                  ! (one moment tracers!)
                  xt_om(I,J,ot,K)=MAX(0.,xt_om(I,J,ot,K)) 
               ENDIF
            ENDDO
         ENDDO
      ENDDO
    ENDDO
    ! mz_bk_20110209+
#endif    
    ! mz_bk_20110209-


  ENDIF !not l_start

  IF (l_trig_mpiom) THEN !mpiom_time_step 
    ! actual water equivalent level (before ALL calculation , i.e. MPIOM and HD)
    l1(:,:)=DDPO(:,:,1)+ZO(:,:)                    &
            -SICTHO(:,:)*RHOICWA-SICSNO(:,:)*RHOSNWA
  ENDIF

  END SUBROUTINE otphysc_global_start
! ------------------------------------------------------------------------

  SUBROUTINE otphysc_global_end

   USE messy_main_tracer_mem_bi, ONLY: xt_om,ntrac_om, ti_om, ON, &
                                       I_MIX, I_ADVECT, I_VDIFF
   USE messy_mpiom_e5,           ONLY: l_trig_mpiom
   USE messy_mpiom_mem_e5,       ONLY: dlxp, dlyp, dzw,zo,ddpo          &
                                     , sictho, sicsno, rhoicwa, rhosnwa & 
#ifdef MPIOM_13B
                                     , OCADPO
#elif defined MPIOM_2000
                                     , ocadpo_trf, octdiff_trf
#endif

   USE messy_mpiom_mem_e5,       ONLY: IE,JE,KE, WETO

   IMPLICIT NONE 

   INTEGER :: ot,K,I,J

  IF (l_trig_mpiom) THEN !mpiom_time_step 
    DO ot=1, ntrac_om !ocean tracer loop
      ! mz_bk_20110208+
      IF (.not. L_SKIP_ADV) THEN
      ! mz_bk_20110208-
      !-------------------------------!
      !           ADVECTION           !
      !-------------------------------!
      !! Set tracer values at dry cells to bottom cell value as done
      !! for temperature and salinity in SBR OCTHER.
      !
      ! mz_bk_20110209+
#if defined MPIOM_13B
      ! mz_bk_20110209-
      DO K=2,KE
         DO J=1,JE
            DO I=1,IE
               IF(WETO(I,J,K).LT.0.5) THEN
                  xt_om(I,J,ot,K)=xt_om(I,J,ot,K-1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ! mz_bk_20110209+
#endif      
      ! mz_bk_20110209-
      IF (ti_om(ot)%tp%meta%cask_i(I_ADVECT) == ON) THEN
      ! mz_bk_20110209+
#if defined MPIOM_13B
      ! mz_bk_20110209-
        CALL OCADPO(xt_om(:,:,ot,:))
#elif MPIOM_2000
        CALL ocadpo_trf(xt_om(:,:,ot,:))
#endif
      ENDIF

#if defined MPIOM_13B
      !-------------------------------!
      ! MASK OF CONCENTRATION ON LAND !
      !-------------------------------!
      !! Reset tracer values at dry cells to undefined.
      !! tracers have to be >= 0.
      DO K=1,KE
         DO J=1,JE
            DO I=1,IE
               IF(WETO(I,J,K).LT.0.5) THEN
                xt_om(I,J,ot,K)=MASK_VALUE
                ELSE      
                  ! to be moved to pdef
                  ! is not done there at the moment because 
                  ! these tracers do not have tendencies
                  ! (one moment tracers!)
                  xt_om(I,J,ot,K)=MAX(0.,xt_om(I,J,ot,K)) 
               ENDIF
            ENDDO
         ENDDO
      ENDDO
#endif
      ! mz_bk_20110208+
      ENDIF
      ! mz_bk_20110208-
      ! mz_bk_20110208+
      IF (.not. L_SKIP_MIX) THEN
      ! mz_bk_20110208-
      !-------------------------------!
      !   ISOPYCNAL  DIFFUSION        !
      !   (Subgrid eddy effects)      !
      !-------------------------------!
      IF (ti_om(ot)%tp%meta%cask_i(I_MIX) == ON) THEN
#if defined MPIOM_13B
        CALL otphysc_ocjitr(xt_om(:,:,ot,:))
#elif MPIOM_2000
        CALL ocjitr_trf(xt_om(:,:,ot,:))
#endif
      ENDIF
      ! mz_bk_20110208+
      ENDIF
      ! mz_bk_20110208-
      ! mz_bk_20110208+
      IF (.not. L_SKIP_DIF) THEN
      ! mz_bk_20110208-
      !-------------------------------!
      !           DIFFUSION           !
      !-------------------------------!
      IF (ti_om(ot)%tp%meta%cask_i(I_VDIFF) == ON) THEN
        CALL OCTDIFF_TRF(xt_om(:,:,ot,:))
      ENDIF
      ! mz_bk_20110208+
      ENDIF
      ! mz_bk_20110208-
      ! mz_bk_20110208+
      IF (.not. L_SKIP_DIL) THEN
      ! mz_bk_20110208-
      !-------------------------------!
      !   DILUTION in UPPER LAYER     !
      !-------------------------------!
      ! actual water equivalent level (after ALL calculation , i.e. MPIOM and HD)
      l1_new(:,:)=DDPO(:,:,1)+ZO(:,:)                    &
                  -SICTHO(:,:)*RHOICWA-SICSNO(:,:)*RHOSNWA
      ! Total freshwater flux is used to calculate new
      ! tracer concentration. The update accounts for all fluxes of
      ! water relevant for tracer concentration changes (e.g.
      ! P-E, runoff, snowmelt (pfresh), and icemelt/growth (pbrine)).
      ! Tracer concentration in all water fluxes is assumed to be zero.
      DO i=1,IE
       DO J=1,JE
          IF (l1_new(i,j).GT.0.0) THEN
             zfac(i,j)=l1(i,j)/l1_new(i,j)
          ELSE 
             zfac(i,j)=1.0
          ENDIF
        ENDDO
      ENDDO
      xt_om(:,:,ot,1) =  xt_om(:,:,ot,1)*zfac(:,:)
      ! mz_bk_20110208+
      ENDIF
      ! mz_bk_20110208-
      ! mz_bk_20110209+
#ifdef MPIOM_13B
      ! mz_bk_20110209-
      !-------------------------------!
      ! MASK OF CONCENTRATION ON LAND !
      !-------------------------------!
      !! Reset tracer values at dry cells to undefined.
      !! tracers have to be >= 0.
      DO K=1,KE
         DO J=1,JE
            DO I=1,IE
               IF(WETO(I,J,K).LT.0.5) THEN
                xt_om(I,J,ot,K)=MASK_VALUE
                ELSE      
                  ! to be moved to pdef
                  ! is not done there at the moment because 
                  ! these tracers do not have tendencies
                  ! (one moment tracers!)
                  xt_om(I,J,ot,K)=MAX(0.,xt_om(I,J,ot,K)) 
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ! mz_bk_20110209+
#endif
      ! mz_bk_20110209-


    ENDDO  !ocean tracer loop
  ENDIF


  END SUBROUTINE otphysc_global_end

  ! ---------------------------------------------------------------------------

  SUBROUTINE otphysc_free_memory

    IMPLICIT NONE

    DEALLOCATE(l1)    ; NULLIFY(l1) 
    DEALLOCATE(l1_new); NULLIFY(l1_new) 
    DEALLOCATE(zfac)  ; NULLIFY(zfac) 


  END SUBROUTINE otphysc_free_memory

  ! ---------------------------------------------------------------------------

#if defined MPIOM_13B
  SUBROUTINE otphysc_ocjitr(iou)


   USE messy_mpiom_mem_e5,       ONLY: dlxp, dlyp,wgo,IE,JE,KE,     &
                                       IE1,JE1,DH, weto, half, uk1o,dduo, &
                                       dlyu,vk1e,ddue,dlxv,t1o,dti,zo,    &
                                       dt,rhosnwa,rhoicwa,almzer,sicsno,  &
                                       sictho,ddpo, bounds_exch
   implicit none

   REAL(DP), INTENT(INOUT) :: iou(IE,JE,KE) ! I/O array

   !LOCAL
   real :: ZSURGM(KE)
   integer :: i,j,k,l
   real :: tsup_p(IE,JE), ssup_p(IE,JE), tlow_p(IE,JE), slow_p(IE,JE)
   real :: roxo,roxu,royo,royu,rozxo,rozxu,rozyo,rozyu
   real :: tm,sm,wun,wob,uwe,uos,vsu,vno,dhi

    !
    ! ...for ocean tracer ...
    !

    DO K=1,KE
      ZSURGM(K)=0.
    ENDDO
    ZSURGM(1)=1.

    DO  k=1,KE

      IF(K.EQ.1) THEN
        DO  J=1,JE
          DO  I=1,IE
            tsup_p(i,j)=iou(i,j,1)
          ENDDO
        ENDDO
      ELSE
        DO  J=1,JE
          DO  I=1,IE
            tsup_p(i,j)=iou(i,j,k-1)
          ENDDO
        ENDDO
      ENDIF

      IF(K.EQ.KE) THEN
        DO  J=1,JE
          DO  I=1,IE
            tlow_p(i,j)=iou(i,j,ke)
          ENDDO
        ENDDO
      ELSE
        DO  J=1,JE
          DO  I=1,IE
            tlow_p(i,j)=iou(i,j,k+1)
          ENDDO
        ENDDO
      ENDIF

      DO  J=2,JE1
        DHI=DTI*DH
        DO  I=2,IE1
          IF ( weto(i,j,k).GT. 0.5 ) THEN
            tm=iou(i,j,k)
            wun= half*dlxp(i,j)*dlyp(i,j)                            &
                 *(wgo(i,j,k+1)+abs(wgo(i,j,k+1)))                 

            wob= half*dlxp(i,j)*dlyp(i,j)                            &
                 *(abs(wgo(i,j,k))-wgo(i,j,k))
            uwe= half*dlyu(i-1,j)*dduo(i-1,j,k)                      &
                 *(uk1o(i-1,j,k)+abs(uk1o(i-1,j,k)))
            uos= half*dlyu(i,j)*dduo(i,j,k)                          &
                 *(abs(uk1o(i,j,k))-uk1o(i,j,k))
            vsu= half*dlxv(i,j)*ddue(i,j,k)                          &
                 *(vk1e(i,j,k)+abs(vk1e(i,j,k)))
            vno= half*dlxv(i,j-1)*ddue(i,j-1,k)                      &
                 *(abs(vk1e(i,j-1,k))-vk1e(i,j-1,k))
            t1o(i,j,k)=                                                   &
                 (tm*dlxp(i,j)*dlyp(i,j)*(                                &
                 ddpo(i,j,k)+almzer+zsurgm(k)*                            &
!we subtract snow and ice                 
                 (zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa) )      &
                 +dt*(wob*(tsup_p(i,j)      -tm)+wun*(tlow_p(i,j)  -tm)   &
                 +  uwe*(iou(i-1,j,k)-tm)+uos*(iou(i+1,j,k)-tm) &
                 +vno*(iou(i,j-1,k)-tm)+vsu*(iou(i,j+1,k)-tm))) &
                 /(dlxp(i,j)*dlyp(i,j)*(                                  &
                 ddpo(i,j,k)+almzer+zsurgm(k)*                            &
!we subtract snow and ice                 
                 (zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa) ) ) 
          ELSE
            t1o(i,j,k)=iou(i,j,k)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DO K=1,KE
      DO  J=2,JE1
        DO  I=2,IE1
          iou(i,j,k)=t1o(i,j,k)
        ENDDO
      ENDDO
    ENDDO

    CALL bounds_exch('p',iou(:,:,:),'otphysc_ocjitr')

  END SUBROUTINE
#endif


! ---------------------------------------------------------------------
  SUBROUTINE otphysc_read_nml_cpl(status, iou)

    ! a2o MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2005

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'otphysc_read_nml_cpl'

    NAMELIST /CPL/ AOFLUX, RINPUT, RIVER_OBJECT, RIVER_CHANNEL &
         ! mz_bk_20110208+
         , L_SKIP_ADV, L_SKIP_MIX, L_SKIP_DIF &
         , L_SKIP_DIL, L_SKIP_AOE, L_SKIP_RIV &
         ! mz_bk_20110208-
         ! mz_bk_20110322+
         , RINPUT_TR
         ! mz_bk_20110322-



    ! LOCAL
    LOGICAL          :: lex      ! file exists ?
    INTEGER          :: fstat    ! file status

    status = 1

    ! mz_bk_20110208+ DEFAULT PARAMETERS, CAN BE OVERWRITTEN BYY NAMELIST
    L_SKIP_ADV = .FALSE.
    L_SKIP_MIX = .FALSE.
    L_SKIP_DIF = .FALSE.
    L_SKIP_DIL = .FALSE.
    L_SKIP_AOE = .FALSE.
    L_SKIP_RIV = .FALSE.
    ! mz_bk_20110208-

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE otphysc_read_nml_cpl

! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_otphysc_e5
! **********************************************************************
