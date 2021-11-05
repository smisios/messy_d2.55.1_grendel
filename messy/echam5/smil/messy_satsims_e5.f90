MODULE MESSY_SATSIMS_E5


! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi

! MESSY CORE LAYER (SMCL)
  USE MESSY_MAIN_CONSTANTS_MEM,   ONLY: dp, STRLEN_MEDIUM

  USE messy_satsims
  
  IMPLICIT NONE
  INTRINSIC :: NULL

  PRIVATE

  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  PUBLIC :: satsims_initialize    ! initialize submodel
  PUBLIC :: satsims_init_memory   ! request memory
  PUBLIC :: satsims_init_coupling ! set pointers for coupling to BM 
                                  ! and other SMs
  PUBLIC :: satsims_physc         ! entry point in time loop (current vector)
  PUBLIC :: satsims_free_memory   ! free allocated memory



  TYPE, PUBLIC :: vmem2d
    REAL(dp), POINTER :: ptr(:,:) => NULL()
  END TYPE vmem2d

  INTEGER, PARAMETER :: ncol=40
  CHARACTER(LEN=STRLEN_MEDIUM) :: cplaer = ''

  ! Variables
  REAL(dp), SAVE :: tautab(0:255)      
  INTEGER, SAVE :: invtau(-20:45000)     
  
  ! ISCCP cloud fraction ofr the 49 ISCCP cloud types
  TYPE (vmem2d), PRIVATE, POINTER :: isccp_cldtypes(:) => NULL()
  ! ISCCP total cloud fraction
  REAL(dp), PRIVATE, POINTER :: isccp_cldfra(:,:)     => NULL()
  ! ISCCP mean cloud optical thickness
  REAL(dp), PRIVATE, POINTER :: isccp_cldtau(:,:)     => NULL()
  ! ISCCP mean cloud top pressure
  REAL(dp), PRIVATE, POINTER :: isccp_cldptop(:,:)    => NULL()
  ! Frequency of cloud occurrence
  REAL(dp), PRIVATE, POINTER :: isccp_cldfreq(:,:)    => NULL()
  ! Gridbox sunlit (1.) or not (0.)
  REAL(dp), PUBLIC,  POINTER :: isccp_sunlit(:,:)     => NULL()
  ! Cloud optical thickness
  REAL(dp), PUBLIC,  POINTER :: isccp_cldtau3d(:,:,:) => NULL()
  ! Cloud emissivity @ 10.5 µm
  REAL(dp), PUBLIC,  POINTER :: isccp_cldemi3d(:,:,:) => NULL()
  REAL(dp), PUBLIC,  POINTER :: isccp_f3d(:,:,:)      => NULL()

  REAL(dp), PUBLIC,  POINTER :: rad_cldtau3d(:,:,:) => NULL()
  REAL(dp), PUBLIC,  POINTER :: rad_cldemi3d(:,:,:) => NULL()
  REAL(dp), PUBLIC,  POINTER :: rad_f3d(:,:,:)      => NULL()


  REAL(dp), POINTER :: press(:,:,:)  => NULL()
  REAL(dp), POINTER :: pressi(:,:,:) => NULL()
  REAL(dp), POINTER :: cover(:,:,:)  => NULL()
  REAL(dp), POINTER :: tsurf(:,:)    => NULL()
  REAL(dp), POINTER :: zi0(:,:)      => NULL()

  ! PRIVATE SUBROTINES
  PRIVATE :: satsims_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE satsims_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'satsims_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CTRL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find free I/O unit
       CALL satsims_read_nml_ctrl(status, iou)  ! read CTRL-namelist
       IF (status /= 0) CALL finish(substr)  ! terminate if error
    END IF
    ! BROADCAST CTRL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(l_isccp, p_io)

    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL satsims_read_nml_cpl(status, iou)  ! read CPL-namelist
       IF (status /= 0) CALL finish(substr)  ! terminate if error
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(cplaer,  p_io)
    CALL p_bcast(taufile, p_io)
    CALL p_bcast(invfile, p_io)
    
    ! ### PERFORM INITIAL SETUP (CALL RESPECTIVE SMCL ROUTINE(S))
    IF (L_ISCCP) THEN
      IF (p_parallel_io) THEN                  ! read only on I/O-PE
        iou = find_next_free_unit(100,200)    ! find free I/O unit
        CALL ISCCP_lookup_read(tautab,invtau,iou,status)
        IF (status /= 0) CALL finish(substr)  ! terminate if error
      END IF
    END IF
    call p_bcast(tautab, p_io)
    call p_bcast(invtau, p_io)

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE satsims_initialize
  ! ====================================================================
  ! ====================================================================
  SUBROUTINE satsims_init_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    ! NOTE: Here, part of ECHAM5 is utilised as MESSy BMIL.
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, DC_GP  &
                                         , DIMID_LON, DIMID_LAT, DIMID_LEV &
                                         , DC_BC &
                                         , gp_nseg, gp_start, gp_cnt &
                                         , gp_meml, gp_memu &
                                         , GP_2D_HORIZONTAL
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute


    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'satsims_init_memory'
    INTEGER                     :: status, i
    CHARACTER(LEN=40)           :: name, name2
    CHARACTER(LEN=2)            :: sctr

    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

    IF (L_ISCCP) THEN
      name=modstr//"_isccp"
      CALL new_channel(status,TRIM(name), reprid=GP_2D_HORIZONTAL)
      CALL channel_halt(substr, status)
      ! (for post-processing geosp, lsp, aps, gboxarea should be in the output)

! 2d mean quantities
      CALL new_channel_object(status,TRIM(name),'sunlit', p2=isccp_sunlit)
      CALL channel_halt(substr,status)
      CALL new_attribute(status, TRIM(name), 'sunlit' &
        , 'long_name', c='fraction of sunlit area of the grid box')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(name), 'sunlit' &
        , 'units', c=' - ' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status,TRIM(name),'cldfra', p2=isccp_cldfra)
      CALL channel_halt(substr,status)
      CALL new_attribute(status, TRIM(name), 'cldfra' &
        , 'long_name', c='total cloud fraction')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(name), 'cldfra' &
        , 'units', c=' - ' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status,TRIM(name),'cldtau', p2=isccp_cldtau)
      CALL channel_halt(substr,status)
      CALL new_attribute(status, TRIM(name), 'cldtau' &
        , 'long_name', c='Total cloud optical thickness (column)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(name), 'cldtau' &
        , 'units', c=' - ' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status,TRIM(name),'cldptop', p2=isccp_cldptop)
      CALL channel_halt(substr,status)
      CALL new_attribute(status, TRIM(name), 'cldptop' &
        , 'long_name', c='Mean cloud top pressure')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(name), 'cldptop' &
        , 'units', c=' hPa ' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status,TRIM(name),'cldfreq', p2=isccp_cldfreq)
      CALL channel_halt(substr,status)
      CALL new_attribute(status, TRIM(name), 'cldfreq' &
        , 'long_name', c='Frequency of occurrence of clouds')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(name), 'cldfreq' &
        , 'units', c=' - ' )
      CALL channel_halt(substr, status)
     
      CALL new_channel_object(status,TRIM(name),'cldtau3d', reprid=GP_3D_MID, &
        p3=isccp_cldtau3d)
      CALL channel_halt(substr,status)
      CALL new_attribute(status, TRIM(name), 'cldtau3d' &
        , 'long_name', c='cloud optical thickness (boxwise)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(name), 'cldtau3d' &
        , 'units', c=' - ' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status,TRIM(name),'cldemi3d', reprid=GP_3D_MID, &
        p3=isccp_cldemi3d)
      CALL channel_halt(substr,status)
      CALL new_attribute(status, TRIM(name), 'cldemi3d' &
        , 'long_name', c='cloud emissivity')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(name), 'cldemi3d' &
        , 'units', c=' - ' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status,TRIM(name),'cldf3d', reprid=GP_3D_MID, &
        p3=isccp_f3d)
      CALL channel_halt(substr,status)
      CALL new_attribute(status, TRIM(name), 'cldf3d' &
        , 'long_name', c='cloud fraction (boxwise)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(name), 'cldf3d' &
        , 'units', c=' - ' )
      CALL channel_halt(substr, status)

      ALLOCATE(isccp_cldtypes(49))

      ! cloud types are sorted by optical thickness and by altitude
      ! first all cloud type with a value for optical thickness are 
      ! stored at each altitude, then for the next optical thickness

!     CLD_type 1 : 0.0 < Tau < 0.3 , 0   < press < 180 hPa
!     CLD_type 2 : 0.0 < Tau < 0.3 , 180 < press < 310 hPa
!     CLD_type 3 : 0.0 < Tau < 0.3 , 310 < press < 440 hPa
!     CLD_type 4 : 0.0 < Tau < 0.3 , 440 < press < 560 hPa
!     CLD_type 5 : 0.0 < Tau < 0.3 , 560 < press < 680 hPa
!     CLD_type 6 : 0.0 < Tau < 0.3 , 680 < press < 800 hPa
!     CLD_type 7 : 0.0 < Tau < 0.3 , 800 < press < inf hPa

!     CLD_type 8 : 0.3 < Tau < 1.3 , 0   < press < 180 hPa
!     CLD_type 9 : 0.3 < Tau < 1.3 , 180 < press < 310 hPa
!     CLD_type 10: 0.3 < Tau < 1.3 , 310 < press < 440 hPa
!     CLD_type 11: 0.3 < Tau < 1.3 , 440 < press < 560 hPa
!     CLD_type 12: 0.3 < Tau < 1.3 , 560 < press < 680 hPa
!     CLD_type 13: 0.3 < Tau < 1.3 , 680 < press < 800 hPa
!     CLD_type 14: 0.3 < Tau < 1.3 , 800 < press < inf hPa

!     CLD_type 15: 1.3 < Tau < 3.6 , 0   < press < 180 hPa
!     CLD_type 16: 1.3 < Tau < 3.6 , 180 < press < 310 hPa
!     CLD_type 17: 1.3 < Tau < 3.6 , 310 < press < 440 hPa
!     CLD_type 18: 1.3 < Tau < 3.6 , 440 < press < 560 hPa
!     CLD_type 19: 1.3 < Tau < 3.6 , 560 < press < 680 hPa
!     CLD_type 20: 1.3 < Tau < 3.6 , 680 < press < 800 hPa
!     CLD_type 21: 1.3 < Tau < 3.6 , 800 < press < inf hPa

!     CLD_type 22: 3.6 < Tau < 9.4 , 0   < press < 180 hPa
!     CLD_type 23: 3.6 < Tau < 9.4 , 180 < press < 310 hPa
!     CLD_type 24: 3.6 < Tau < 9.4 , 310 < press < 440 hPa
!     CLD_type 25: 3.6 < Tau < 9.4 , 440 < press < 560 hPa
!     CLD_type 26: 3.6 < Tau < 9.4 , 560 < press < 680 hPa
!     CLD_type 27: 3.6 < Tau < 9.4 , 680 < press < 800 hPa
!     CLD_type 28: 3.6 < Tau < 9.4 , 800 < press < inf hPa

!     CLD_type 29: 9.4 < Tau < 23  , 0   < press < 180 hPa
!     CLD_type 30: 9.4 < Tau < 23  , 180 < press < 310 hPa
!     CLD_type 31: 9.4 < Tau < 23  , 310 < press < 440 hPa
!     CLD_type 32: 9.4 < Tau < 23  , 440 < press < 560 hPa
!     CLD_type 33: 9.4 < Tau < 23  , 560 < press < 680 hPa
!     CLD_type 34: 9.4 < Tau < 23  , 680 < press < 800 hPa
!     CLD_type 35: 9.4 < Tau < 23  , 800 < press < inf hPa

!     CLD_type 36: 23  < Tau < 60  , 0   < press < 180 hPa
!     CLD_type 37: 23  < Tau < 60  , 180 < press < 310 hPa
!     CLD_type 38: 23  < Tau < 60  , 310 < press < 440 hPa
!     CLD_type 39: 23  < Tau < 60  , 440 < press < 560 hPa
!     CLD_type 40: 23  < Tau < 60  , 560 < press < 680 hPa
!     CLD_type 41: 23  < Tau < 60  , 680 < press < 800 hPa
!     CLD_type 42: 23  < Tau < 60  , 800 < press < inf hPa

!     CLD_type 43: 60  < Tau < inf , 0   < press < 180 hPa
!     CLD_type 44: 60  < Tau < inf , 180 < press < 310 hPa
!     CLD_type 45: 60  < Tau < inf , 310 < press < 440 hPa
!     CLD_type 46: 60  < Tau < inf , 440 < press < 560 hPa
!     CLD_type 47: 60  < Tau < inf , 560 < press < 680 hPa
!     CLD_type 48: 60  < Tau < inf , 680 < press < 800 hPa
!     CLD_type 49: 60  < Tau < inf , 800 < press < inf hPa

      DO i=1,49
        IF ( i.LT.10 ) THEN
          WRITE(name2, '(a8,i1)') 'cldtype0',i
          WRITE(sctr, '(a1,i1)') '0',i
        ELSE
          WRITE(name2, '(a7,i2)') 'cldtype',i
          WRITE(sctr, '(i2)') i
        ENDIF
        CALL new_channel_object(status,TRIM(name),TRIM(name2), &
          p2=isccp_cldtypes(i)%ptr)
        CALL channel_halt(substr,status)
        CALL new_attribute(status, TRIM(name), TRIM(name2) &
        , 'long_name', c='Fraction covered by ISCCP cloud type'//sctr)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(name), TRIM(name2) &
          , 'units', c=' - ' )
        CALL channel_halt(substr, status)

      ENDDO

    END IF

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

  END SUBROUTINE satsims_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE satsims_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basmodel and other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    ! NOTE: Here, part of ECHAM5 is utilised as MESSy BMIL.
    USE messy_main_mpi_bi,           ONLY: finish
    USE messy_main_channel,          ONLY: get_channel_object &
                                         , get_channel_info     ! op_pj_20130408
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: error_bi, warning_bi ! op_pj_20130408

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'satsims_init_coupling'
    INTEGER                     :: status
    CHARACTER(LEN=40)           :: name
!    CHARACTER(LEN=*), PARAMETER :: cplaer="_tanstd"
    CHARACTER(LEN=25) :: cldradcha = '' ! op_pj_20130408
    CHARACTER(LEN=25) :: radcha    = '' ! op_pj_20130408
    CHARACTER(LEN=25) :: clcovcha  = '' ! um_hr_20190306

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

    ! ### set pointers to channel objects here
    IF (L_ISCCP) THEN
 
      CALL get_channel_object(status, 'ECHAM5', 'press', p3=press)
      CALL channel_halt(substr//': channel object for pressure not found!',status)
      CALL get_channel_object(status, 'ECHAM5', 'pressi', p3=pressi)
      CALL channel_halt(substr//': channel object for interface pressure not found!',&
        status)
      CALL get_channel_object(status, 'ECHAM5', 'tsurf', p2=tsurf)
      CALL channel_halt(substr//': channel object for surface temperature not found!',&
        status)

      ! op_pj_20130408+
      CALL get_channel_info(status, 'rad')
      IF (status /=0) THEN
         CALL warning_bi( &
              'channel rad not available, trying old rad4all ...' &
              , substr)
         CALL get_channel_info(status, 'rad4all')
         IF (status /= 0) THEN
            CALL error_bi(' ... rad4all also not available!', substr)
         ELSE
            radcha = 'rad4all'
         ENDIF
      ELSE
         radcha = 'rad'
      ENDIF
      ! op_pj_20130408-

      CALL get_channel_object(status, TRIM(radcha), 'zi0', p2=zi0)
      CALL channel_halt(substr//': channel object for zi0 not found!',&
        status)

! um_hr_20190306+
!!$   CALL get_channel_object(status, 'cloud', 'aclc', p3=cover)
!!$   CALL channel_halt(substr//': channel object for cloud cover not found!',&
!!$     status)
      CALL get_channel_info(status, 'cloud')
      IF (status /= 0) THEN
         CALL warning_bi( &
              'channel cloud not available, trying new CRM superparameterisation...' &
              , substr)
         CALL get_channel_info(status, 'crm')
         IF (status /= 0) THEN
            CALL error_bi(' ... CRM superparameterization also not available!', substr)
         ELSE
            clcovcha = 'crm'
            CALL get_channel_object(status, TRIM(clcovcha), 'cld', p3=cover)
            CALL channel_halt(substr//': channel object for cloud cover not found!',&
                 status)
         ENDIF
      ELSE
         clcovcha = 'cloud'
         CALL get_channel_object(status, TRIM(clcovcha), 'aclc', p3=cover)
         CALL channel_halt(substr//': channel object for cloud cover not found!',&
              status)
      ENDIF
! um_hr_20190306-
      !

      ! op_pj_20130408+
      CALL get_channel_info(status, 'cloudopt01')
      IF (status /=0) THEN
         CALL warning_bi( &
              'channel cloudopt01 not available, trying old rad4all ...' &
              , substr)
         CALL get_channel_info(status, 'rad4all')
         IF (status /= 0) THEN
            CALL error_bi(' ... rad4all also not available!', substr)
         ELSE
            cldradcha = 'rad4all'
         ENDIF
      ELSE
         cldradcha = 'cloudopt01'
      ENDIF
      ! op_pj_20130408-

      ! op_pj_20130408-

      name="isccp_cldtau"//cplaer
      CALL get_channel_object(status, TRIM(cldradcha), TRIM(name), p3=rad_cldtau3d)
      CALL channel_halt(substr//': channel object for cloud optical thickness not found!',&
        status)
            
      name="isccp_cldemi"//cplaer
      CALL get_channel_object(status, TRIM(cldradcha), TRIM(name), p3=rad_cldemi3d)
      CALL channel_halt(substr//': channel object for cloud emissivity not found!',&
        status)

      name="isccp_f"//cplaer
      CALL get_channel_object(status, TRIM(cldradcha), TRIM(name), p3=rad_f3d)
      CALL channel_halt(substr//': channel object for cloud fraction (from radiation) not found!',&
        status)
      
    END IF
    
    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE satsims_init_coupling
! =========================================================================

  SUBROUTINE satsims_free_memory

    IF (ASSOCIATED(isccp_cldtypes)) THEN
       DEALLOCATE(isccp_cldtypes)
       NULLIFY(isccp_cldtypes)
    END IF

  END SUBROUTINE satsims_free_memory
! =========================================================================

  SUBROUTINE SATSIMS_PHYSC
    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
    USE messy_main_data_bi,         ONLY: tm1, qm1, tte_3d, qte_3d
    USE messy_main_timer,           ONLY: time_step_len

    REAL(dp) :: temp(kproma,nlev), q(kproma,nlev)

    temp(1:kproma,1:nlev) = tm1(1:kproma,1:nlev,jrow) + &
                            tte_3d(1:kproma,1:nlev,jrow) * time_step_len
    q(1:kproma,1:nlev)    = qm1(1:kproma,1:nlev,jrow) + &
                            qte_3d(1:kproma,1:nlev,jrow) * time_step_len


    IF (L_ISCCP) THEN
      isccp_cldtau3d(:,:,jrow) = rad_cldtau3d(:,:,jrow) 
      isccp_cldemi3d(:,:,jrow) = rad_cldemi3d(:,:,jrow) 
      isccp_f3d(:,:,jrow)      = rad_f3d(:,:,jrow) 

      isccp_sunlit(1:kproma,jrow) = 0._dp
      WHERE ( zi0(1:kproma,jrow) > 0._dp )
        isccp_sunlit(1:kproma,jrow) = 1._dp
      ENDWHERE

      CALL ISCCP_SIMULATOR(kproma,nlev,jrow,temp,q)
    ENDIF

  END SUBROUTINE SATSIMS_PHYSC

!-----------------------------------------------------------------------------

  SUBROUTINE ISCCP_SIMULATOR(kproma,klev,jrow,temp,q)

    USE MESSY_SATSIMS_ISCCP,       ONLY: isccp_cloud_types_v3_4
    USE MESSY_MAIN_CONSTANTS_MEM,  ONLY: emsfc_lw => cemiss

    INTEGER  :: kproma,klev,jrow
    REAL(dp) :: temp(kproma,klev)
    REAL(dp) :: q(kproma,klev)


    ! Input to ISCCP-Simulator to be set here
    INTEGER :: debug, debugcol, top_height, overlap
    INTEGER, DIMENSION(kproma)    :: sunlit
    !  10.5 micron longwave emissivity of stratiform
    REAL(dp), DIMENSION(kproma, klev) :: dem_s   
    !  10.5 micron longwave emissivity of convective
    REAL(dp), DIMENSION(kproma, klev) :: dem_c

    ! Local
    ! Convective cloud cover [0-1]
    REAL(dp), DIMENSION(kproma, klev)   :: fc
    ! In-cloud optical thickness (convective clouds)
    REAL(dp), DIMENSION(kproma, klev)   :: tauc
    ! In-cloud optical thickness (stratiform clouds)
    REAL(dp), DIMENSION(kproma, klev)   :: taus
    
    REAL(dp), DIMENSION(kproma, 7, 7)  :: fq_isccp
    REAL(dp), POINTER :: ct(:,:)
    REAL(dp), DIMENSION(kproma) :: pcldfra
    REAL(dp), DIMENSION(kproma) :: pcldtau
    REAL(dp), DIMENSION(kproma) :: pcldptop

    INTEGER i,itau,ipres


    ! Avoid invalid cloud fraction
    WHERE  ( isccp_cldtau3d(1:kproma,1:klev,jrow).LE.0._dp ) 
      isccp_f3d(1:kproma,1:klev,jrow) = 0._dp
    ENDWHERE
    ! Set general input to ISCCP-Simulator
    debug=0
    debugcol=0
    top_height=1 ! adjustment of cloud top height
    overlap=3    ! random-maximum
    
    WHERE ( isccp_sunlit(1:kproma,jrow) > 0._dp )
       sunlit(1:kproma) = 1
    ELSEWHERE
       sunlit(1:kproma) = 0
    ENDWHERE
    
    fc(1:kproma,1:klev)    = 0._dp  ! convective clouds = 0
    tauc(1:kproma,1:klev)  = 0._dp  ! convective clouds = 0
    dem_c(1:kproma,1:klev) = 0._dp  ! convective clouds = 0
    
    
    ! In-cloud optical thickness and emissivity
    WHERE ( isccp_f3d(1:kproma,1:klev,jrow) > 0._dp ) 
       taus(1:kproma,1:klev) = isccp_cldtau3d(1:kproma,1:klev,jrow) &
                             / isccp_f3d(1:kproma,1:klev,jrow)
       dem_s(1:kproma,1:klev) = isccp_cldemi3d(1:kproma,1:klev,jrow)
     ELSEWHERE
       taus(1:kproma,1:klev) = 0._dp
       dem_s(1:kproma,1:klev) = 0._dp
     ENDWHERE

     CALL isccp_cloud_types_v3_4(&
       debug, debugcol, &
       kproma, &
       sunlit, &
       klev, ncol, &
       press(1:kproma,:,jrow), pressi(1:kproma,:,jrow), &
       q(1:kproma,:), &
       isccp_f3d(1:kproma,:,jrow), fc(1:kproma,:), &
       taus(1:kproma,:), tauc(1:kproma,:), &
       top_height, overlap, &
       tautab, invtau, &
       tsurf(1:kproma,jrow),&
       emsfc_lw,&
       temp(1:kproma,:),&
       dem_s,&
       dem_c,&
       fq_isccp, &
       pcldfra, pcldptop, pcldtau)

     WHERE ( pcldfra(1:kproma) > 0._dp .AND. pcldtau(1:kproma) > 0._dp ) 
       isccp_cldfreq(1:kproma,jrow) = 1._dp
     ELSEWHERE
       isccp_cldfreq(1:kproma,jrow) = 0._dp
     END WHERE

     i=0
     DO itau=1,7
       DO ipres=1,7
         i=i+1
         ct => isccp_cldtypes(i)%ptr
         ct(1:kproma,jrow) = fq_isccp(1:kproma,itau,ipres) 
       ENDDO
     ENDDO

     isccp_cldtau(1:kproma,jrow)  = pcldtau(1:kproma)  &
                                  * isccp_cldfreq(1:kproma,jrow)
     isccp_cldptop(1:kproma,jrow) = pcldptop(1:kproma) &
                                  * isccp_cldfreq(1:kproma,jrow)
     isccp_cldfra(1:kproma,jrow)  = pcldfra(1:kproma) 

   END SUBROUTINE ISCCP_SIMULATOR

!===============================================================================
   SUBROUTINE satsims_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ cplaer, taufile, invfile

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='satsims_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE satsims_read_nml_cpl
  ! ====================================================================
END MODULE MESSY_SATSIMS_E5
