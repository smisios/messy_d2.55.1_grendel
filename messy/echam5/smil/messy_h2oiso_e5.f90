! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL H2OISO 
!
! Author : Roland Eichinger, DLR-IPA, 2014
!
! References: see messy_h2oiso.f90
!
! **********************************************************************

! **********************************************************************
MODULE messy_h2oiso_e5
  ! **********************************************************************

  USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi, &
       error_bi, warning_bi
  USE messy_main_constants_mem,  ONLY: M_H2O, M_HDO, M_HH18O
  USE messy_h2oiso
  USE messy_h2oiso_convect
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,    ONLY: mtend_add_l,      &
       mtend_add_g,      &
       mtend_id_tracer,  &
       mtend_id_q,       &
       mtend_id_t, mtend_id_xl, mtend_id_xi, &
       mtend_get_start_l, &
       mtend_get_handle, &
       mtend_register,   &
       mtend_checkreg
#endif
  USE messy_main_tools,          ONLY: PTR_4D_ARRAY, PTR_5D_ARRAY, PTR_3D_ARRAY
  USE messy_main_mpi_bi,         ONLY: p_io
  USE messy_main_grid_def_mem_bi, ONLY: nproma

  !WISO++
  USE messy_main_data_wiso_bi
  USE messy_main_tools_wiso
  !WISO--

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! CPL namelist switches
  LOGICAL :: l_steady       = .FALSE.
  LOGICAL :: l_nocloud_dd   = .FALSE.
  LOGICAL :: l_noconvect_dd = .FALSE.

  ! switches for external processes
  LOGICAL :: l_convect = .FALSE.
  LOGICAL :: l_msbm    = .FALSE. 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: flt_stratreg => NULL()
  LOGICAL :: l_e5vdiff = .FALSE.

  ! was q registered from CH4 or MECCA?
  LOGICAL :: lreg_ch4   = .FALSE.
  LOGICAL :: lreg_mecca = .FALSE.
  ! MSBM
  LOGICAL :: lreg_msbm = .FALSE.

  ! GLOBAL PARAMETERS

  !#ifdef MESSYTENDENCY
  ! op_pj_20140924: does not compile without MESSYTENDENCY
  ! HANDLES FOR ADDING UP THE TENDENCIES
  INTEGER              :: my_handle_vdiff
  INTEGER              :: my_handle_convect
!!$  INTEGER              :: my_handle_cloud
  INTEGER              :: my_handle_numer
  INTEGER              :: my_handle_chem
  INTEGER              :: my_handle_msbm

  ! POINTERS FOR PROCESS TENDENCIES FROM MTEND
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_qte_vdiff   => NULL()
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_tte_surf    => NULL()
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_vom_vdiff   => NULL()
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_vol_vdiff   => NULL()
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_qte_ch4     => NULL()
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_qte_mecca   => NULL()
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_qte_msbm     => NULL()  
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_xlte_msbm    => NULL()  
  REAL(KIND=DP),POINTER,DIMENSION(:,:,:):: mtend_xite_msbm    => NULL()  

  ! this is for convec2
  REAL(dp), POINTER :: pilab(:,:,:)   => NULL()

  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  PUBLIC :: h2oiso_initialize     ! initialize submodel
  PUBLIC :: h2oiso_init_memory    ! request memory
  PUBLIC :: h2oiso_new_tracer     ! define new tracers
  PUBLIC :: h2oiso_init_coupling  ! set pointers for coupling to BM and
  !                                 other SMs
  PUBLIC :: h2oiso_init_tracer    ! initialize tracers
  PUBLIC :: h2oiso_local_start
  PUBLIC :: h2oiso_global_start   ! entry point in time loop (all vectors)
  PUBLIC :: h2oiso_radheat        ! called in messy_main_control_e5
  PUBLIC :: h2oiso_convec         ! called in messy_main_control_e5
  PUBLIC :: h2oiso_mixlo          ! called in physc
  PUBLIC :: h2oiso_physc          ! entry point in time loop (current vector)
  PUBLIC :: h2oiso_write_output   ! reset accumulated fields
  PUBLIC :: h2oiso_free_memory    ! free allocated memory

  ! PRIVATE SUBROTINES
  !PRIVATE :: h2oiso_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE h2oiso_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

#ifndef MESSYTENDENCY
    CALL error_bi('H2OISO requires MESSYTENDENCY! '//&
         &'Please reconfigure / recompile with --enable-MESSYTENDENCY',substr)
#endif

    ! TELL OTHERS THAT WE ARE ON ...
    l_wiso = .TRUE.
    
    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL h2oiso_read_nml_cpl(status, iou) ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF

    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(l_steady, p_io)
    CALL p_bcast(l_noconvect_dd, p_io)
    CALL p_bcast(l_nocloud_dd, p_io)

    ! TRANSFER SWITCH(ES) TO BMIL
    l_wiso_nocloud_dd = l_nocloud_dd
    
    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE h2oiso_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_new_tracer

    ! ------------------------------------------------------------------
    ! This subroutine is used to define new tracers. See
    ! http://www.atmos-chem-phys.net/8/1677   (including supplement !)
    ! for full documentation.
    ! ------------------------------------------------------------------

    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tracer,          ONLY: set_tracer                 &
                                        , R_MOLARMASS, ON, OFF       &
                                        , AMOUNTFRACTION             &
                                        , new_tracer, I_ADVECT       &
                                        , I_CONVECT, I_VDIFF

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_new_tracer'
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: setname = GPTRSTR
    INTEGER :: jwiso, jphase

#ifdef MESSYTENDENCY
    ! GET HANDLES AND REGISTER FOR ADDING TENDENCIES IN MTEND
    my_handle_vdiff   = mtend_get_handle('h2oiso_radheat')
    CALL mtend_register(my_handle_vdiff,mtend_id_tracer)

    my_handle_convect = mtend_get_handle('h2oiso_convect')
    CALL mtend_register(my_handle_convect,mtend_id_tracer)

    my_handle_numer   = mtend_get_handle('h2oiso_numer')
    CALL mtend_register(my_handle_numer,mtend_id_tracer)

    my_handle_chem     = mtend_get_handle('h2oiso_chem')
    CALL mtend_register(my_handle_chem,mtend_id_tracer)

    my_handle_msbm     = mtend_get_handle('h2oiso_msbm')
    CALL mtend_register(my_handle_msbm,mtend_id_tracer)
#endif

    CALL start_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

    idiso(:,:) = 0 ! tracer IDs for water isotopes (3 phases, 3 isotopes)

    ! define water isotope tracers 
    DO jphase = 1, kphase
       DO jwiso = 1, mwiso
          CALL new_tracer(status, setname &
               , TRIM('H2OISO'//TRIM(isotag(jwiso))//TRIM(phasetag(jphase))) &
               , modstr, idx = idiso(jphase, jwiso)     &
               , quantity = AMOUNTFRACTION, unit = 'kg/kg')
          CALL tracer_halt(substr,status)
          CALL set_tracer(status, setname, idiso(jphase,jwiso)          &
               , R_molarmass, r=molarmasses(jwiso))
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname, idiso(jphase,jwiso), I_ADVECT, ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname, idiso(jphase,jwiso), I_CONVECT, OFF)
          CALL tracer_halt(substr, status)
          IF(jphase .eq. i_vap) THEN
             CALL set_tracer(status, setname, idiso(jphase,jwiso), I_VDIFF, OFF)
          ELSE
             CALL set_tracer(status, setname, idiso(jphase,jwiso), I_VDIFF, ON)
          ENDIF
          CALL tracer_halt(substr,status)
       END DO
    END DO

    Call end_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

  END SUBROUTINE h2oiso_new_tracer
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_init_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! ------------------------------------------------------------------

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_2D_HORIZONTAL
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                           new_attribute
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, nlev, ngpblks
    USE messy_main_constants_mem,    ONLY: FLAGGED_BAD
    
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_init_memory'
    INTEGER                     :: status
    INTEGER                     :: ji
    REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  mem_3d => NULL()
    REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  mem_2dh => NULL()
    INTEGER :: jwiso

    CALL start_message_bi(modstr,'MEMORY ALLOCATION',substr)

    ! ### ALLOCATE OWN MEMORY HERE
    ! -------------
    ! allocate pointer arrays for water isotopologue tracers
    ALLOCATE(isotracm1(kphase))
    ALLOCATE(isotracte(kphase))

    ! allocate other global varibales
    ALLOCATE(pwl_surf(nproma))
    ALLOCATE(psn_surf(nproma))
    ALLOCATE(pws_surf(nproma))
    ALLOCATE(tsl_vdiff(nproma))
    ALLOCATE(tsl_surf(nproma))
    ALLOCATE(grndcapc_surf(nproma))

    CALL end_message_bi(modstr,'MEMORY ALLOCATION',substr)

    ! CHANNEL AND CHANNEL OBJECTS
    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

    CALL new_channel(status, modstr, reprid=GP_2D_HORIZONTAL,lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    ! allocate horizontal 2-D surface pointer arrays and create channel objects
    ! for new diagnostic output (3 Dim)

    ALLOCATE(pwiso_2dh(i_2nmax))
    ALLOCATE(pwiso_2dh_5d(i_2nmax))
    DO ji = 1,i_2nmax
       ALLOCATE(pwiso_2dh_5d(ji)%ptr(nproma,mwiso,ngpblks,1,1))
       pwiso_2dh(ji)%ptr => pwiso_2dh_5d(ji)%ptr(:,:,:,1,1)
       DO jwiso = 1, mwiso
          mem_2dh => pwiso_2dh_5d(ji)%ptr(:,jwiso,:,:,:)
          CALL new_channel_object(status, modstr &
               , TRIM(varname_2dh(ji))//'_'//TRIM(isotag(jwiso)) &
               , mem=mem_2dh, reprid=GP_2D_HORIZONTAL)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, modstr &
               , TRIM(varname_2dh(ji))//'_'//TRIM(isotag(jwiso)) &
               , 'missing_value', r=FLAGGED_BAD, loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, modstr &
               , TRIM(varname_2dh(ji))//'_'//TRIM(isotag(jwiso)) &
               , '_FillValue', r=FLAGGED_BAD, loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
          !
       ENDDO
    ENDDO

    ! allocate arrays from core and nullify 3D pointers for tracers

    ! allocate 3-D pointer arrays and create channel objects
    ! for new diagnostic output
    ALLOCATE(pwiso_3d_5d(i_3nmax))
    ALLOCATE(pwiso_3d(i_3nmax))
    DO ji = 1, i_3nmax
       ALLOCATE(pwiso_3d_5d(ji)%ptr(nproma,nlev,mwiso,ngpblks,1))
       pwiso_3d(ji)%ptr => pwiso_3d_5d(ji)%ptr(:,:,:,:,1)
       DO jwiso = 1, mwiso
          mem_3d => pwiso_3d_5d(ji)%ptr(:,:,jwiso,:,:)
          CALL new_channel_object(status, modstr &
               , TRIM(varname_3d(ji))//'_'//TRIM(isotag(jwiso)) &
               , mem=mem_3d, reprid=GP_3D_MID)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, modstr &
               , TRIM(varname_3d(ji))//'_'//TRIM(isotag(jwiso)) &
               , 'missing_value', r=FLAGGED_BAD, loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, modstr &
               , TRIM(varname_3d(ji))//'_'//TRIM(isotag(jwiso)) &
               , '_FillValue', r=FLAGGED_BAD, loverwrite=.TRUE.)
          CALL channel_halt(substr, status)
          !
       ENDDO
    ENDDO

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

  END SUBROUTINE h2oiso_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_init_coupling(flag)

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! and to transfer process-tendency-pointers from the
    ! TENDENCY-submodel to here via the mtend_request routine
    ! ------------------------------------------------------------------

    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR, xtm1, xtte
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer,           ONLY: get_tracer
    USE messy_main_channel,          ONLY: get_channel_object &
         , set_channel_object_restreq &
         , get_channel_info
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_data_bi,          ONLY: modstr_base=>modstr
    USE messy_main_mpi_bi,           ONLY: finish
#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,      ONLY: mtend_request,  &
         mtend_id_u,     &
         mtend_id_v,     &
         mtend_id_t
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: flag

#ifdef MESSYTENDENCY
    ! STRINGS FOR REQUESTING THE PROCESS TENDENCIES FROM MTEND
    CHARACTER(len=32), PARAMETER  ::  str_vdiff   =  'e5vdiff'
    CHARACTER(len=32), PARAMETER  ::  str_surf    =  'surface'
    CHARACTER(len=32), PARAMETER  ::  str_ch4     =  'ch4'
    CHARACTER(len=32), PARAMETER  ::  str_mecca   =  'mecca'
    CHARACTER(len=32), PARAMETER  ::  str_msbm    =  'msbm'
#endif
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_init_coupling'
    INTEGER                     :: status
    INTEGER :: otest
    INTEGER :: jphase, jwiso

    SELECT CASE(flag)
    CASE (1)

       CALL start_message_bi(modstr,'H2OISO_INIT_COUPLING',substr)

       ! get additional fields for cloud
       CALL get_channel_object(status, modstr_base, 'qtec', p3=pqtec)
       IF (status /= 0) &
            call finish(substr,'channel object for qtec not found')
       
       CALL get_channel_object(status, modstr_base, 'vdiffp', p3=pvdiffp)
       IF (status /= 0) &
            call finish(substr,'channel object for vdiffp not found')
       
       CALL get_channel_object(status, 'e5vdiff', 'wet_tmp', p2=wet_tmp)
       IF (status /= 0) &
            call finish(substr,'channel object for wet_tmp not found')
       
       CALL get_channel_object(status, 'e5vdiff', 'evapot_2d', p2=evapot_2d)
       IF (status /= 0) &
            call finish(substr,'channel object for evapot_2d not found')
       
       ! get additional fields for vdiff/surf
       CALL get_channel_object(status,'g3b','sn',p2=sn)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status,'g3b','wl',p2=wl)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status,'g3b','snc',p2=snc)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status,'g3b','gld',p2=gld)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status,'g3b','tsl',p2=tsl)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status,'g3b','grndcapc',p2=grndcapc)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status,'g3b','orostd',p2=orostd)
       CALL channel_halt(substr, status)
       
       ! get additional fields for convec
       CALL get_channel_object(status, modstr_base, 'ilab', p3=pilab )
       IF (status /= 0) &
            call finish(substr,'channel object for ilab not found')
       
       ! get the pwisosw_d fields, from .nc file
       CALL get_channel_object(status,'import_grid' &
            , 'WISOSW_WISOSW_HHO',p2=pwisosw_d_HHO)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            ,'WISOSW_WISOSW_HH18O',p2=pwisosw_d_HH18O)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'WISOSW_WISOSW_HDO',p2=pwisosw_d_HDO)
       CALL channel_halt(substr, status)
       
       ! Get the reservoir fields, from .nc file for steady-state restart
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_ws_HHO',p2=ws_HHO)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_wl_HHO',p2=wl_HHO)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_sn_HHO',p2=sn_HHO)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_ws_HH18O',p2=ws_HH18O)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_wl_HH18O',p2=wl_HH18O)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_sn_HH18O',p2=sn_HH18O)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_ws_HDO',p2=ws_HDO)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_wl_HDO',p2=wl_HDO)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status,'import_grid' &
            , 'H2OISO_RES_sn_HDO',p2=sn_HDO)
       CALL channel_halt(substr, status)
       
       ! for safety reasons, before pointing on tracers:
       ! check if water isotopologue tracers are set in correct order
       DO jphase= 1, kphase
          DO jwiso = 1, mwiso
             CALL get_tracer(status, GPTRSTR, &
                  TRIM('H2OISO'//TRIM(isotag(jwiso))//TRIM(phasetag(jphase))) &
                  , idx=idiso(jphase,jwiso))
             IF (jwiso .ne. 1 .and. jphase .ne. 1) THEN
                IF (otest .ne. (idiso(jphase, jwiso)-1)) &
                     CALL error_bi(' Water isotopologues not in correct order!'&
                     , substr)
             ENDIF
             otest = idiso(jphase,jwiso)
          ENDDO
       ENDDO
       
       ! let pointer arrays point on isotope tracers
       DO jphase = 1, kphase
          isotracm1(jphase)%ptr => xtm1(:,:,idiso(jphase,1):idiso(jphase,mwiso),:)
          isotracte(jphase)%ptr => xtte(:,:,idiso(jphase,1):idiso(jphase,mwiso),:)
       END DO
       
       CALL get_channel_info(status,'convect')
       l_convect = (status == 0)
       
       CALL get_channel_info(status,'e5vdiff')
       l_e5vdiff =  (status == 0)
       
       CALL get_channel_info(status, 'msbm')
       l_msbm = (status == 0)
       IF (l_msbm) THEN
          CALL get_channel_object(status, 'msbm', 'STRAT_region', p3=flt_stratreg)
          CALL channel_halt(substr, status)
       END IF

    CASE (2)

#ifdef MESSYTENDENCY
       ! REQUEST SOME PROCESS TENDENCIES VIA MTEND
       CALL mtend_request(str_vdiff,   mtend_id_q     , mtend_qte_vdiff   )
       CALL mtend_request(str_surf,    mtend_id_t     , mtend_tte_surf    )
       CALL mtend_request(str_vdiff,   mtend_id_v     , mtend_vol_vdiff   )
       CALL mtend_request(str_vdiff,   mtend_id_u     , mtend_vom_vdiff   )
       ! chemical tendency of q
       lreg_ch4 = mtend_checkreg(str_ch4, mtend_id_q)
       IF (lreg_ch4) &
            CALL mtend_request(str_ch4, mtend_id_q, mtend_qte_ch4)    
       lreg_mecca = mtend_checkreg(str_mecca, mtend_id_q)
       IF (lreg_mecca) &
            CALL mtend_request(str_mecca, mtend_id_q, mtend_qte_mecca)    
       ! tendency of msbm
       lreg_msbm = mtend_checkreg(str_msbm, mtend_id_q)
       IF (lreg_msbm) THEN
          CALL mtend_request(str_msbm, mtend_id_q, mtend_qte_msbm)    
          CALL mtend_request(str_msbm, mtend_id_xl, mtend_xlte_msbm)    
          CALL mtend_request(str_msbm, mtend_id_xi, mtend_xite_msbm)    
       ENDIF
#endif

       CALL end_message_bi(modstr,'H2OISO_INIT_COUPLING',substr)
    CASE DEFAULT

       CALL error_bi('UNKNOWN FLAG',substr)
    END SELECT

  END SUBROUTINE h2oiso_init_coupling
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_init_tracer

    ! ------------------------------------------------------------------
    ! This subroutine is used to initialise the additional
    ! water tracers according to the respective phase
    ! ------------------------------------------------------------------

    USE messy_main_tracer,          ONLY: get_tracer, tracer_iniflag
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_timer,           ONLY: lresume

    IMPLICIT NONE

    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_init_tracer'
    LOGICAL                     :: linit
    INTEGER :: jphase, jwiso
    
    CALL start_message_bi(modstr, 'H2OISO tracer initialization', substr)

    IF (lresume) RETURN

    IF (l_steady)THEN
       DO jphase= 1, kphase
          DO jwiso = 1, mwiso
             ! SET POINTER TO TRACER
             CALL get_tracer(status, GPTRSTR, &
                  TRIM('H2OISO'//TRIM(isotag(jwiso))//TRIM(phasetag(jphase))) &
                  , idx=idiso(jphase,jwiso))
             CALL tracer_iniflag(status, GPTRSTR, idiso(jphase,jwiso) &
                  , lget=linit)
             CALL tracer_halt(substr, status)
             IF (.NOT. linit) THEN
                CALL error_bi('H2OISO TRACER WAS NOT INITIALISED!', substr)
             END IF
          ENDDO
       ENDDO
    ENDIF

    CALL end_message_bi(modstr, 'H2OISO tracer initialization', substr)

  END SUBROUTINE h2oiso_init_tracer
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_global_start

    USE messy_main_grid_def_mem_bi, ONLY: nlev, ngpblks 
    USE messy_main_grid_def_bi,     ONLY: coslat_2d,sinlat_2d
    USE messy_main_data_bi,         ONLY: q, xl, xi, qm1, xlm1, xim1          &
                                        ,  press_3d

    USE messy_main_tracer_mem_bi, ONLY: xt, xtm1, xtte
    USE messy_main_timer,         ONLY: lstart

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_global_start'
    INTEGER:: jr, jk, jl, jwiso
    REAL(DP):: r_climtp(nproma,ngpblks), twisostra(mwiso)
    LOGICAL:: l_strato(nproma,nlev,ngpblks) 

    qstart: IF (lstart) THEN

      IF (l_steady) THEN ! start from steady-state conditions

          ! Set q,xl,xi to the HHO tracer which is steady-state
          DO jk = 1, nlev
             DO jr = 1, ngpblks
                DO jl = 1, nproma
                   q(jl,jk,jr)    = xt(jl,jk,idiso(i_vap,i_HHO),jr)
                   xl(jl,jk,jr)   = xt(jl,jk,idiso(i_liq,i_HHO),jr)
                   xi(jl,jk,jr)   = xt(jl,jk,idiso(i_ice,i_HHO),jr)

                   qm1(jl,jk,jr)  = xtm1(jl,jk,idiso(i_vap,i_HHO),jr)
                   xlm1(jl,jk,jr) = xtm1(jl,jk,idiso(i_liq,i_HHO),jr)
                   xim1(jl,jk,jr) = xtm1(jl,jk,idiso(i_ice,i_HHO),jr)
                ENDDO
             ENDDO
          ENDDO

       ELSE ! start a spin-up run

          ! initial isotope deviation from SMOW in the stratosphere
          ! analogue to twisoatm
          twisostra(i_HHO)   = 0.0_dp
          twisostra(i_HH18O) = -20.0_dp/1000.0_dp ! no idea! same as wisoatm
          twisostra(i_HDO)   = -550.0_dp/1000.0_dp ! see Randel et al. 2012

          l_strato(:,:,:) = .FALSE.

          ! INITIALIZE THE ISOTOPE TRACERS, THIS IS PRELIMINARY
          DO jk = 1, nlev
             DO jr = 1, ngpblks
                DO jl = 1, nproma
                   ! Calculate a climatological tropopause in Pa
                   r_climtp(jl,jr) = &
                        (300.0_dp - 215.0_dp * (coslat_2d(jl,jr))**2) * 100.0_dp
                   ! SET l_strato TRUE if above the clim stratosphere
                   IF (press_3d(jl,jk,jr).lt.r_climtp(jl,jr)) &
                        l_strato(jl,jk,jr) = .TRUE.

                   ! For Troposphere set HDO with initial delta value -150
                   IF (.NOT.l_strato(jl,jk,jr)) THEN
                      DO jwiso = 1, mwiso
                         xt(jl,jk,idiso(i_vap,jwiso),jr)   = &
                              tnat(jwiso) * q(jl,jk,jr)    * &
                              (1.0_dp + twisoatm(jwiso))
                         xt(jl,jk,idiso(i_liq,jwiso),jr)   = &
                              tnat(jwiso) * xl(jl,jk,jr)   * &
                              (1.0_dp + twisoatm(jwiso)) ! this is zero anyway
                         xt(jl,jk,idiso(i_ice,jwiso),jr)   = &
                              tnat(jwiso) * xi(jl,jk,jr)   * &
                              (1.0_dp + twisoatm(jwiso)) ! this is zero anyway
                         xtm1(jl,jk,idiso(i_vap,jwiso),jr) = &
                              tnat(jwiso) * qm1(jl,jk,jr)  * &
                              (1.0_dp + twisoatm(jwiso))
                         xtm1(jl,jk,idiso(i_liq,jwiso),jr) = &
                              tnat(jwiso) * xlm1(jl,jk,jr) * &
                              (1.0_dp + twisoatm(jwiso)) ! this is zero anyway
                         xtm1(jl,jk,idiso(i_ice,jwiso),jr) = &
                              tnat(jwiso) * xim1(jl,jk,jr) * &
                              (1.0_dp + twisoatm(jwiso)) ! this is zero anyway
                      ENDDO
                      ! For Stratosphere set HDO with initial delta
                      ! value -550 to shorten spin-up time 
                   ELSE
                      DO jwiso = 1, mwiso
                         xt(jl,jk,idiso(i_vap,jwiso),jr)   = &
                              tnat(jwiso) * q(jl,jk,jr)    * &
                              (1.0_dp + twisostra(jwiso))
                         xt(jl,jk,idiso(i_liq,jwiso),jr)   = &
                              tnat(jwiso) * xl(jl,jk,jr)   * &
                              (1.0_dp + twisostra(jwiso)) ! this is zero anyway
                         xt(jl,jk,idiso(i_ice,jwiso),jr)   = &
                              tnat(jwiso) * xi(jl,jk,jr)   * &
                              (1.0_dp + twisostra(jwiso)) ! this is zero anyway
                         xtm1(jl,jk,idiso(i_vap,jwiso),jr) = &
                              tnat(jwiso) * qm1(jl,jk,jr)  * &
                              (1.0_dp + twisostra(jwiso))
                         xtm1(jl,jk,idiso(i_liq,jwiso),jr) = &
                              tnat(jwiso) * xlm1(jl,jk,jr) * &
                              (1.0_dp + twisostra(jwiso)) ! this is zero anyway
                         xtm1(jl,jk,idiso(i_ice,jwiso),jr) = &
                              tnat(jwiso) * xim1(jl,jk,jr) * &
                              (1.0_dp + twisostra(jwiso)) ! this is zero anyway
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       END IF
    END IF qstart

    ! reset tracer values of HHOvap liq and ice to its counterparts
    ! and set in case of negative values
    ! the HDO and HH18O tracers to a standard delta value

    IF (.NOT. lstart) THEN
       DO jk = 1, nlev
          ! since the negative m1 values are overwritten in tracer_pdef
          ! we have to set the HHO-tracers to the actual water here again ...
          DO jr = 1, ngpblks
             DO jl = 1, nproma
                isotracm1(i_vap)%ptr(jl,jk,i_HHO,jr) = qm1(jl,jk,jr)
                isotracm1(i_liq)%ptr(jl,jk,i_HHO,jr) = xlm1(jl,jk,jr)
                isotracm1(i_ice)%ptr(jl,jk,i_HHO,jr) = xim1(jl,jk,jr)
             ENDDO

             ! ... and the negative m1 values of HDO and HH18O are set to
             ! a delta value of VSMOW
             DO jwiso = 2, mwiso
                DO jl = 1, nproma
                   IF (qm1(jl,jk,jr) .lt. cwisomin) THEN
                      CALL h2oiso_constdelta_tracer(&
                           isotracm1(i_vap)%ptr(jl,jk,jwiso,jr), &
                           tnat(jwiso), qm1(jl,jk,jr), &
                           molarmasses(i_HHO), molarmasses(jwiso) )
                   ENDIF
                   IF (xlm1(jl,jk,jr) .lt. cwisomin) THEN
                      CALL h2oiso_constdelta_tracer(&
                           isotracm1(i_liq)%ptr(jl,jk,jwiso,jr), &
                           tnat(jwiso), xlm1(jl,jk,jr), &
                           molarmasses(i_HHO), molarmasses(jwiso) )
                   ENDIF
                   IF (xim1(jl,jk,jr) .lt. cwisomin) THEN
                      CALL h2oiso_constdelta_tracer(&
                           isotracm1(i_ice)%ptr(jl,jk,jwiso,jr), &
                           tnat(jwiso), xim1(jl,jk,jr), &
                           molarmasses(i_HHO), molarmasses(jwiso) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE h2oiso_global_start
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_local_start

    USE messy_main_constants_mem,    ONLY: DP, tmelt
    USE messy_main_grid_def_mem_bi,  ONLY: kproma, jrow, nlev
    USE messy_main_data_bi,          ONLY: ws, qm1, xlm1, xim1, tslclim
    USE messy_main_timer,            ONLY: lstart

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_local_start'
    INTEGER :: jmonth
    REAL(DP):: zwisosurf(nproma,mwiso)
    REAL(DP):: tsurface(nproma)
    INTEGER :: jt
    
    ! here we initialize the snow depth (pws), the surface soil wetness (psn)
    ! and the water in skin reservoir (pwl) for the isotopologues

    ! calculate annual mean temperature from climate land surface temperatures
    IF (lstart) THEN

       IF (l_steady) THEN ! start from steady-state conditions

          pwiso_2dh(i_sn)%ptr(1:kproma,i_HHO,jrow)   = ws_HHO(1:kproma,jrow)
          pwiso_2dh(i_wl)%ptr(1:kproma,i_HHO,jrow)   = wl_HHO(1:kproma,jrow)
          pwiso_2dh(i_ws)%ptr(1:kproma,i_HHO,jrow)   = sn_HHO(1:kproma,jrow)
          pwiso_2dh(i_sn)%ptr(1:kproma,i_HH18O,jrow) = ws_HH18O(1:kproma,jrow)
          pwiso_2dh(i_wl)%ptr(1:kproma,i_HH18O,jrow) = wl_HH18O(1:kproma,jrow)
          pwiso_2dh(i_ws)%ptr(1:kproma,i_HH18O,jrow) = sn_HH18O(1:kproma,jrow)
          pwiso_2dh(i_sn)%ptr(1:kproma,i_HDO,jrow)   = ws_HDO(1:kproma,jrow)
          pwiso_2dh(i_wl)%ptr(1:kproma,i_HDO,jrow)   = wl_HDO(1:kproma,jrow)
          pwiso_2dh(i_ws)%ptr(1:kproma,i_HDO,jrow)   = sn_HDO(1:kproma,jrow)

          ws(1:kproma,jrow) = ws_HHO(1:kproma,jrow)
          wl(1:kproma,jrow) = wl_HHO(1:kproma,jrow)
          sn(1:kproma,jrow) = sn_HHO(1:kproma,jrow)

       ELSE ! start a spin-up run

          tsurface(1:kproma) = 0.0_dp
          DO jmonth = 1, 12
             tsurface(1:kproma) = tsurface(1:kproma) + &
                  tslclim(1:kproma,jrow,jmonth) 
          END DO
          tsurface(1:kproma) = tsurface(1:kproma)/12.0_dp

          DO jt=1,mwiso
             zwisosurf(1:kproma,jt) = &
                  twisosur1(jt)*(tsurface(1:kproma) - tmelt) - twisosur2(jt)
             IF (jt == i_HH18O) THEN
                ! if tracer = 18O: zwisosurf = min (-4./1000., zxtatm)
                zwisosurf(1:kproma,jt) = &
                     MIN( -4.0_dp/1000.0_dp, zwisosurf(1:kproma,jt) )
             ELSE IF (jt == i_HDO) THEN
                ! if tracer = HDO: zxtatm = min (-22./1000., zxtatm)
                zwisosurf(1:kproma,jt) = &
                     MIN( -22.0_dp/1000.0_dp, zwisosurf(1:kproma,jt) )
             ENDIF
          END DO

          ! INITIALIZE THE ISOTOPE EQUIVALENTS OF ws, wl and sn here
          DO jt = 1, mwiso
             IF (jt == i_HHO) THEN
                ! if tracer = H2-16O: set tracer reservoirs =
                ! evap. mask * corresponding wetness reservoir
                pwiso_2dh(i_sn)%ptr(1:kproma,jt,jrow) = sn(1:kproma,jrow)
                pwiso_2dh(i_wl)%ptr(1:kproma,jt,jrow) = wl(1:kproma,jrow)
                pwiso_2dh(i_ws)%ptr(1:kproma,jt,jrow) = ws(1:kproma,jrow)
             ELSE
                ! if tracer not H2-16O:
                ! set tracer reservoirs = evap. mask * nat.isotope-ratio *
                ! corresponding wetness reservoir * (1.+zxtatm)
                pwiso_2dh(i_sn)%ptr(1:kproma,jt,jrow) = &
                     tnat(jt)*sn(1:kproma,jrow)*(1.0_dp+zwisosurf(1:kproma,jt))
                pwiso_2dh(i_wl)%ptr(1:kproma,jt,jrow) = &
                     tnat(jt)*wl(1:kproma,jrow)*(1.0_dp+zwisosurf(1:kproma,jt))
                pwiso_2dh(i_ws)%ptr(1:kproma,jt,jrow) = &
                     tnat(jt)*ws(1:kproma,jrow)*(1.0_dp+zwisosurf(1:kproma,jt))
             ENDIF
          ENDDO
       END IF ! l_steady
    ENDIF ! lstart

    tsl_vdiff(1:kproma) = tsl(1:kproma,jrow)

  END SUBROUTINE h2oiso_local_start
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_radheat

    ! this routine calculates the vertical diffuision (vdiff) but it is
    ! called from main entry point messy_radheat (right after vdiff!).

    USE messy_main_constants_mem, ONLY: g, cvdifts
    USE messy_main_timer,         ONLY: time_step_len, delta_time
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nproma, nlev &
         , nlevp1, nlevm1
    USE messy_main_data_bi,       ONLY:                             &
          wsmx    , tsw   , tsi   , qsl   , qsw   , qsi           &
         , um1     , vm1   , qm1                                   &
         , phum    , ws    , cvs   , cvw   , vgrat , chl           &
         , loland_2d, aphm1 , cfhl  , cfhw  , cfhi                 &
         , cfh     , ebsh  , slm   , seacov, icecov                &
         , qslnew  , vol_3d   , vom_3d   , qte_3d   , loglac_2d

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='h2oiso_radheat'

    REAL(dp):: pwisosw_d(nproma,mwiso,jrow)
    REAL(dp):: zwisoqsw(nproma,mwiso), zwisoqsi(nproma,mwiso), zwisoqsl(nproma,mwiso) ! surfhum
    REAL(dp):: zwisocsat(nproma,mwiso), zwisocair(nproma,mwiso) ! vapcoeff
    REAL(dp):: zwisokin(nproma,mwiso) ! kinfac
    REAL(dp):: zwisoeqnl(nproma,mwiso), zwisoeqnw(nproma,mwiso), zwisoeqni(nproma,mwiso) ! RMscheme
    REAL(dp):: zwisofqnl(nproma,mwiso), zwisofqnw(nproma,mwiso), zwisofqni(nproma,mwiso) ! RMscheme
    REAL(dp):: zwisoqdif(nproma,nlev,mwiso) ! ZQDIF/RMscheme
    REAL(dp):: zwisodqdt2(nproma,nlev,mwiso)
    REAL(dp):: ztmst, zcons13, ztpfac3
    REAL(dp):: pqm1(nproma,nlev), pum1(nproma,nlev), pvm1(nproma,nlev)
    REAL(dp):: zqdif(nproma,nlev), zqdif_jl(nproma,nlev), zdqdt
    REAL(dp):: pwisoqm1(nproma,nlev,mwiso)
    INTEGER :: jl, jk, jwiso
    
    IF (.NOT. l_e5vdiff) RETURN

    pwiso_2dh(i_evap)%ptr(1:kproma,:,jrow)    = 0.0_dp
    pwiso_2dh(i_evapl)%ptr(1:kproma,:,jrow)   = 0.0_dp
    pwiso_2dh(i_evapw)%ptr(1:kproma,:,jrow)   = 0.0_dp
    pwiso_2dh(i_evapi)%ptr(1:kproma,:,jrow)   = 0.0_dp
    pwiso_2dh(i_evaplac)%ptr(1:kproma,:,jrow) = 0.0_dp
    pwiso_2dh(i_evapwac)%ptr(1:kproma,:,jrow) = 0.0_dp
    pwiso_2dh(i_evapiac)%ptr(1:kproma,:,jrow) = 0.0_dp

    ztmst = time_step_len
    zcons13 = 1.0_dp/ztmst
    ztpfac3 = 1.0_dp-(1.0_dp/cvdifts)

    ! HERE THE START VALUES ARE GENERATED LIKE IN VDIFF
    ! SINCE WE CALL THIS SUBROUTINE IN PHYSC (AFTER VDIFF)
    ! WE SUBTRACT THE VDIFF TENDENCIES AGAIN (WITH MESSY_TENDENCY).
#ifndef MESSYTENDENCY
! op_pj_20170904: qqq see above
#else
    DO jk=1,nlev
       pqm1(1:kproma,jk) = qm1(1:kproma,jk,jrow) + (qte_3d(1:kproma,jk,jrow) &
            - mtend_qte_vdiff(1:kproma,jk,jrow)) * ztmst
       pum1(1:kproma,jk) = um1(1:kproma,jk,jrow) + (vom_3d(1:kproma,jk,jrow) &
            - mtend_vom_vdiff(1:kproma,jk,jrow)) * ztmst
       pvm1(1:kproma,jk) = vm1(1:kproma,jk,jrow) + (vol_3d(1:kproma,jk,jrow) &
            - mtend_vol_vdiff(1:kproma,jk,jrow)) * ztmst
    ENDDO
#endif

    ! compute the iso start values in the same way
    DO jk = 1, nlev
       DO jwiso = 1, mwiso
          pwisoqm1(1:kproma,jk,jwiso) = &
               isotracm1(i_vap)%ptr(1:kproma,jk,jwiso,jrow) &
               + isotracte(i_vap)%ptr(1:kproma,jk,jwiso,jrow) * ztmst
       ENDDO
    ENDDO

    ! set the pwisosw_d's (read from .nc file) for isotopic
    ! sea surface compositions
    pwisosw_d(1:kproma,i_HHO,jrow)   = pwisosw_d_HHO(1:kproma,jrow)
    pwisosw_d(1:kproma,i_HH18O,jrow) = pwisosw_d_HH18O(1:kproma,jrow)
    pwisosw_d(1:kproma,i_HDO,jrow)   = pwisosw_d_HDO(1:kproma,jrow)

    ! here we compute zqdif and zqdif_jl recursively
    ! from the vdiff TENDENCY of qte
#ifndef MESSYTENDENCY
! op_pj_20170904: qqq see above
#else
    DO jl = 1, kproma
       DO jk = 1, nlev
          zdqdt = mtend_qte_vdiff(jl,jk,jrow)
          zqdif(jl,jk) = (zdqdt/zcons13)+pqm1(jl,jk)
       ENDDO
       zqdif_jl(jl,nlev) = zqdif(jl,nlev) - ztpfac3 * pqm1(jl,nlev)
    ENDDO
#endif

    ! H2OISO-VDIFF------------------------
    ! These subroutines are for calculating the 
    ! changes of water isotopes in the lowest model layer
    ! Thus it is a part of vdiff. The isotopic fluxes in and out
    ! of the lowest layer are computed here and then included into the
    ! tendency in vdiff for the isotope tracers.
    ! As the isotopes are tracers, the diffusion for all other layers
    ! goes via the tracer infrastructure.

    ! H2OISO-VDIFF------------------------
    ! This subroutine is to calculate the isotopic humidity
    ! saturation over land/water/ice *zwisoqsl*, *zwisoqsw*, 
    ! *zwisoqsi* including equilibrium and effective
    ! (over ice) fractionation. (No fractionation over land)

    CALL h2oiso_vdiff_surfhum( kproma                                                           &
         , wsmx(1:kproma,jrow) , tsw(1:kproma,jrow)  , tsi(1:kproma,jrow)   &
         , qsl(1:kproma,jrow)  , qsw(1:kproma,jrow)  , qsi(1:kproma,jrow)   &
         , zwisoqsw(1:kproma,:), zwisoqsi(1:kproma,:), zwisoqsl(1:kproma,:) &
         , pwiso_2dh(i_ws)%ptr(1:kproma,:,jrow), pwisosw_d(1:kproma,:,jrow) )

    ! H2OISO-VDIFF------------------------
    ! This subroutine is to calculate the isotopic vapour 
    ! coefficients *zwisocair* and *zwisocsat* which are
    ! used to calculate the surface flux j.
    ! Herefore the isotopic surface reservoirs are calculated first.
    CALL h2oiso_vdiff_vapcoeff( kproma, nlev   &
         , pum1(1:kproma,:)  , pvm1(1:kproma,:)    , pqm1(1:kproma,:)    &
         , qsl(1:kproma,jrow), phum(1:kproma,jrow) , sn(1:kproma,jrow)   &
         , ws(1:kproma,jrow) , wsmx(1:kproma,jrow) , cvs(1:kproma,jrow)  &
         , cvw(1:kproma,jrow), vgrat(1:kproma,jrow), wl(1:kproma,jrow)   &
         , chl(1:kproma,jrow), wet_tmp(1:kproma,jrow)                    &
         , pwiso_2dh(i_ws)%ptr(1:kproma,:,jrow)                          &
         , pwiso_2dh(i_sn)%ptr(1:kproma,:,jrow)                          &
         , pwiso_2dh(i_wl)%ptr(1:kproma,:,jrow)                          &
         , zwisocsat(1:kproma,:), zwisocair(1:kproma,:)                  &
         , loland_2d(1:kproma, jrow)                                      )

    ! H2OISO-VDIFF------------------------
    ! This subrotine is to calculate *zwisokin* which is the
    ! kinetic fractionation factor at air/sea interface
    CALL h2oiso_vdiff_kinfac( kproma, nlev                           &
         , pum1(1:kproma,:), pvm1(1:kproma,:)     &
         , zwisokin(1:kproma,:)                   &
         , loland_2d(1:kproma,jrow)                )

    ! H2OISO-VDIFF------------------------
    ! This subrotine is to calculate the EN and FN coefficients
    ! *zwisoeqnl*, *zwisofqnl*, *zwisoeqnw*, *zwisofqnw*,
    ! *zwisoeqni*, *zwisofqni* 
    ! for moisture of the Richtmyer-Morton-Scheme
    ! and *zwisoqdif*, the result for vdiff without the surface fluxes
    CALL h2oiso_vdiff_RMscheme( kproma, nlev, nlevp1, nlevm1                 &
         , aphm1(1:kproma,1:nlevp1), cvdifts                                 &
         , cfhl(1:kproma,jrow), cfhw(1:kproma,jrow), cfhi(1:kproma,jrow)     &
         , cfh(1:kproma,:,jrow), ebsh(1:kproma,:,jrow), zwisokin(1:kproma,:) &
         , zwisoeqnl(1:kproma,:), zwisofqnl(1:kproma,:)                      &
         , zwisoeqnw(1:kproma,:), zwisofqnw(1:kproma,:)                      &
         , zwisoeqni(1:kproma,:), zwisofqni(1:kproma,:)                      &
         , pwisoqm1(1:kproma,:,:)                                            &
         , zwisoqdif(1:kproma,:,:)  )

    ! H2OISO-VDIFF------------------------
    ! This subroutine is to calculate the main result *zwisoqdif* from 
    ! vdiff which determines the isotopic changes in the lowest
    ! model layer. (Other changes are done via the tracer infrastructure)
    CALL h2oiso_vdiff_zqdif( kproma, nlev, nlevm1, ztmst                       &
         , cfhl(1:kproma,jrow), cfhw(1:kproma,jrow)  , cfhi(1:kproma,jrow)     &
         , slm(1:kproma,jrow) , seacov(1:kproma,jrow), icecov(1:kproma,jrow)   &
         , cvdifts, zqdif_jl(1:kproma,:)                                       &
         , qslnew(1:kproma,jrow), zwisoqdif(1:kproma,:,:)                      &
         , zwisoeqnl(1:kproma,:), zwisofqnl(1:kproma,:), zwisoeqnw(1:kproma,:) &
         , zwisofqnw(1:kproma,:), zwisoeqni(1:kproma,:), zwisofqni(1:kproma,:) &
         , zwisoqsw(1:kproma,:) , zwisoqsi(1:kproma,:), zwisocsat(1:kproma,:)  &
         , zwisocair(1:kproma,:), ebsh(1:kproma,:,jrow)                        &
         , zwisodqdt2(1:kproma,:,:), pwisoqm1(1:kproma,:,:)  )

    DO jwiso = 1, mwiso
#ifndef MESSYTENDENCY
       isotracte(i_vap)%ptr(1:kproma,:,jwiso,jrow) = &
            isotracte(i_vap)%ptr(1:kproma,:,jwiso,jrow) + zwisodqdt2(1:kproma,:,jwiso)
#else
       CALL mtend_add_l(my_handle_vdiff, idiso(i_vap, jwiso), &
            px=zwisodqdt2(1:kproma,:,jwiso) )
#endif
    ENDDO

    ! H2OISO-VDIFF------------------------
    ! This subroutine is to calculate the isotopic evaporation values 
    ! *pwisoevapl* and *pwisoevapot* which are needed furtheron in surf
    CALL h2oiso_vdiff_mflux( kproma, nlev, g                             &
         , cfhl(1:kproma,jrow), cfhw(1:kproma,jrow), cfhi(1:kproma,jrow) &
         , cvdifts, pqm1(1:kproma,:)                                     &
         , ztmst, zqdif(1:kproma,:)                                      &
         , qslnew(1:kproma,jrow), loland_2d(1:kproma,jrow)               &
         , pwiso_2dh(i_evapl)%ptr(1:kproma,:,jrow)                       &
         , pwiso_2dh(i_evapot)%ptr(1:kproma,:,jrow)                      &
         , zwisoqdif(1:kproma,:,:), zwisocsat(1:kproma,:)                &
         , zwisocair(1:kproma,:)                                         &
         , pwisoqm1(1:kproma,:,:), delta_time                            &
         , slm(1:kproma,jrow) , seacov(1:kproma,jrow), icecov(1:kproma,jrow) &
         , zwisoqsw(1:kproma,:) , zwisoqsi(1:kproma,:)  &
         , pwiso_2dh(i_evap)%ptr(1:kproma,:,jrow)       &
         , pwiso_2dh(i_evapw)%ptr(1:kproma,:,jrow)      &
         , pwiso_2dh(i_evapi)%ptr(1:kproma,:,jrow)      &
         , pwiso_2dh(i_evaplac)%ptr(1:kproma,:,jrow)    &
         , pwiso_2dh(i_evapwac)%ptr(1:kproma,:,jrow)    &
         , pwiso_2dh(i_evapiac)%ptr(1:kproma,:,jrow) )

    ! Here we generate some values which we need in h2oiso_mixlo,
    ! because since h2oiso_mixlo is called after surf, these values
    ! have changed by then

    pwl_surf(1:kproma) = wl(1:kproma,jrow)
    psn_surf(1:kproma) = sn(1:kproma,jrow)
    pws_surf(1:kproma) = ws(1:kproma,jrow)

    tsl_surf(1:kproma) = tsl(1:kproma,jrow)
    grndcapc_surf(1:kproma) = grndcapc(1:kproma,jrow)

  END SUBROUTINE h2oiso_radheat
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_convec

    USE messy_main_mpi_bi,        ONLY: finish
    USE messy_main_timer,         ONLY: time_step_len, delta_time
    USE messy_main_grid_def_mem_bi, ONLY: nproma, kproma, nlev, nlevp1 &
         , nlevm1, jrow , nn
    USE messy_main_data_bi,       ONLY: &
          tm1, qm1, um1, vm1, xlm1, xim1             &
         , tte_3d, qte_3d, vom_3d, vol_3d, xlte_3d, xite_3d &
         , vervel_3d, xtec , app1, aphp1, geopot_3d   &
         , rsfc_2d, ssfc_2d, aprc, aprs , loland_2d   &
         , qflux 

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='h2oiso_convec'
    REAL(DP):: tte_2(nproma,nlev), qte_2(nproma,nlev)
    REAL(DP):: vom_2(nproma,nlev), vol_2(nproma,nlev)
    REAL(DP):: xlte_2(nproma,nlev), xite_2(nproma,nlev)
    REAL(DP):: xtec_2(nproma,nlev)
    REAL(DP):: rsfc_2(nproma), ssfc_2(nproma)
    REAL(DP):: aprc_2(nproma), aprs_2(nproma)
    REAL(DP):: isoqte_convec(nproma,nlev,mwiso)
    INTEGER :: jlab_2(nproma,nlev)
    INTEGER :: itype(nproma)
    REAL(DP):: zqtec(nproma,nlev), qhfla(nproma)
    LOGICAL :: lookupoverflow = .FALSE.
    REAL(DP):: HHOVAP_convect(nproma,nlev), HDOVAP_convect(nproma,nlev)
    INTEGER :: jl, jk, jwiso
    
    IF (.NOT. l_convect) RETURN

    ! set apr* and *vi to 0
    ! no idea where this is done for the actual values
    pwiso_2dh(i_aprc)%ptr(1:kproma,:,jrow) = 0.0_dp
    pwiso_2dh(i_aprs)%ptr(1:kproma,:,jrow) = 0.0_dp
    pwiso_2dh(i_aprl)%ptr(1:kproma,:,jrow) = 0.0_dp
    pwiso_2dh(i_qvi)%ptr(1:kproma,:,jrow)  = 0.0_dp
    pwiso_2dh(i_xlvi)%ptr(1:kproma,:,jrow) = 0.0_dp
    pwiso_2dh(i_xivi)%ptr(1:kproma,:,jrow) = 0.0_dp

    ! set *tec to 0 like in physc
    pwiso_3d(i_xtec)%ptr(1:kproma,:,:,jrow) = 0.0_dp
    pwiso_3d(i_qtec)%ptr(1:kproma,:,:,jrow) = 0.0_dp

    zqtec(1:kproma,1:nlev) = 0.0_dp
    itype(1:kproma)        = 0

    qhfla(1:kproma) = qflux(1:kproma,jrow)

    jlab_2(1:kproma,1:nlev) = NINT(pilab(1:kproma,1:nlev,jrow))

    ! --------------------------------------------------------
    ! Here the '*_2's are set on the actual value and only then
    ! given further because we don't want to change anything in
    ! the actual hydrological cycle
    tte_2(1:kproma,:) = tte_3d(1:kproma,:,jrow)
    qte_2(1:kproma,:) = qte_3d(1:kproma,:,jrow)
    vom_2(1:kproma,:) = vom_3d(1:kproma,:,jrow)
    vol_2(1:kproma,:) = vol_3d(1:kproma,:,jrow)
    xlte_2(1:kproma,:) = xlte_3d(1:kproma,:,jrow)
    xite_2(1:kproma,:) = xite_3d(1:kproma,:,jrow)

    xtec_2(1:kproma,:) = xtec(1:kproma,:,jrow)
    rsfc_2(1:kproma) = rsfc_2d(1:kproma,jrow)
    ssfc_2(1:kproma) = ssfc_2d(1:kproma,jrow)
    aprc_2(1:kproma) = aprc(1:kproma,jrow)
    aprs_2(1:kproma) = aprs(1:kproma,jrow)
    ! --------------------------------------------------------

    ! set start tendency (only for q, because xi and xl aren't changed
    ! here anyway)
    isoqte_convec(1:kproma,:,:) = isotracte(i_vap)%ptr(1:kproma,:,:,jrow)

    CALL h2oiso_cucall(kproma, kproma, nlev, nlevp1, nlevm1, &
         jlab_2(1:kproma,:), &
         tm1(1:kproma,:,jrow), qm1(1:kproma,:,jrow)  , um1(1:kproma,:,jrow),  &
         vm1(1:kproma,:,jrow), xlm1(1:kproma,:,jrow) , xim1(1:kproma,:,jrow), &
         tte_2(1:kproma,:)   , qte_2(1:kproma,:)     , vom_2(1:kproma,:),     &
         vol_2(1:kproma,:)   , xlte_2(1:kproma,:)    , xite_2(1:kproma,:),    &
         vervel_3d(1:kproma,:,jrow)  , xtec_2(1:kproma,:)    ,                &
         zqtec(1:kproma,:)   , qhfla(1:kproma)       ,                        &
         app1(1:kproma,:)    , aphp1(1:kproma,:), geopot_3d(1:kproma,:,jrow), &
         rsfc_2(1:kproma)    , ssfc_2(1:kproma),                              &
         aprc_2(1:kproma)    , aprs_2(1:kproma),                              &
         itype(1:kproma)     , loland_2d(1:kproma,jrow),                      &
         !---wiso-code
         mwiso,                                   &
         isotracm1(i_vap)%ptr(1:kproma,:,:,jrow), &
         isotracm1(i_liq)%ptr(1:kproma,:,:,jrow), &
         isotracm1(i_ice)%ptr(1:kproma,:,:,jrow), &
         isoqte_convec(1:kproma,:,:),             &
         isotracte(i_liq)%ptr(1:kproma,:,:,jrow), &
         isotracte(i_ice)%ptr(1:kproma,:,:,jrow), &
         pwiso_3d(i_xtec)%ptr(1:kproma,:,:,jrow), &
         pwiso_3d(i_qtec)%ptr(1:kproma,:,:,jrow), &
         pwiso_2dh(i_rsfc)%ptr(1:kproma,:,jrow) , &
         pwiso_2dh(i_ssfc)%ptr(1:kproma,:,jrow),  &
         pwiso_2dh(i_aprc)%ptr(1:kproma,:,jrow),  &
         pwiso_2dh(i_aprs)%ptr(1:kproma,:,jrow),  &
         time_step_len, lookupoverflow,           &
         nn, delta_time)

    IF (lookupoverflow) THEN 
       CALL FINISH ('h2oiso_convec - lookuperror')
    ENDIF

    ! subtract old tendency from new tendency to get convect-tendency
    DO jwiso = 1, mwiso
       isoqte_convec(1:kproma,:,jwiso) = &
            isoqte_convec(1:kproma,:,jwiso) - &
            isotracte(i_vap)%ptr(1:kproma,:,jwiso,jrow) 

       ! sc_re_20150701+
       !=====================================================================
       ! For sensitivity study withouth the effect of convection on deltaD
       IF(l_noconvect_dd) THEN
          HHOVAP_convect(1:kproma,:) = &
               isotracm1(i_vap)%ptr(1:kproma,:,i_HHO,jrow) + &
               isotracte(i_vap)%ptr(1:kproma,:,i_HHO,jrow) * time_step_len
          HDOVAP_convect(1:kproma,:) = &
               isotracm1(i_vap)%ptr(1:kproma,:,i_HDO,jrow) + &
               isotracte(i_vap)%ptr(1:kproma,:,i_HDO,jrow) * time_step_len
          !===================================================================
          ! Keep the delta value of vapour after convect alike the delta value
          ! of before!  
          !===================================================================  
          DO jk = 1, nlev
             DO jl = 1, kproma
                IF( HHOVAP_convect(jl,jk)  .gt. cwisomin) THEN
                   isoqte_convec(jl,jk,i_HDO) = isoqte_convec(jl,jk,i_HHO) * &
                        HDOVAP_convect(jl,jk)/HHOVAP_convect(jl,jk)
                ELSE
                   isoqte_convec(jl,jk,i_HDO) = isoqte_convec(jl,jk,i_HHO) &
                        * tnat(i_HDO)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       !=====================================================================
       ! sc_re_20150701-

       ! Add convect-tendency to old tendency
#ifndef MESSYTENDENCY
       isotracte(i_vap)%ptr(1:kproma,:,jwiso,jrow) = &
            isotracte(i_vap)%ptr(1:kproma,:,jwiso,jrow) + &
            isoqte_convec(1:kproma,:,jwiso)
#else
       CALL mtend_add_l(my_handle_convect, idiso(i_vap,jwiso) , &
            px=isoqte_convec(1:kproma,:,jwiso))
#endif
    END DO

  END SUBROUTINE h2oiso_convec
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_mixlo

    USE messy_main_timer,         ONLY: lstart, delta_time, time_step_len
    USE messy_main_grid_def_mem_bi, ONLY: ngl
    USE messy_main_constants_mem, ONLY: rhoh2o=>rho_H2O
    USE messy_main_grid_def_mem_bi,ONLY: nproma, kproma, nlev, jrow
    USE messy_main_data_bi,       ONLY: &
           cvs, cvw, ssfl_2d , ssfc_2d       &
         , loland_2d, loglac_2d,         u10, v10, vlt, tsoil, evwsd  &
         , rsfl_2d, rsfc_2d, wsmx, slm, tm1, tte_3d                   &
         , snc_surf, gld_surf

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='h2oiso_mixlo'
    INTEGER,          PARAMETER :: jpgrnd  = 5       ! see mo_parameters
    REAL(dp),         PARAMETER :: cwlmax = 2.E-4_dp ! see mo_physc2 and iniphy
    REAL(dp),         PARAMETER :: cvinter = 0.25_dp ! see mo_vegetation 
    !                                                ! and iniphy
    INTEGER :: jk, jt, jl
    REAL(dp):: zwisoevsnd(nproma,mwiso), zwisoevwld(nproma,mwiso)
    REAL(dp):: zwisoraind(nproma,mwiso), zwisosnowd(nproma,mwiso)
    REAL(dp):: zwisoevwsd(nproma,mwiso), zwisomlres(nproma,mwiso)
    REAL(dp):: zwisosncmelt(nproma,mwiso), zwisosnmel(nproma,mwiso)
    REAL(dp):: zwisoros(nproma,mwiso), zwisoevttd(nproma,mwiso)
    REAL(dp):: zwisorogl(nproma,mwiso), zwisosn(nproma,mwiso)
    REAL(dp):: zwisodrain(nproma,mwiso), zdtime

    REAL(dp):: zsn_tmp1(nproma), zsn_tmp2(nproma)
    REAL(dp):: zsncmelt(nproma), pwlmx_tmp(nproma)
    REAL(dp):: zraind_tmp(nproma), zsnmel(nproma)
    REAL(dp):: zmprcp2_tmp(nproma), zsnmlt(nproma)
    REAL(dp):: zsnowd(nproma), zevsnd(nproma)
    REAL(DP):: pwl(nproma), psn(nproma), pws(nproma)
    REAL(DP):: tte_ms(nproma,nlev), ptsl(nproma)

    REAL(DP):: pgld(nproma), psnc(nproma)
    
    ! set these to 0, no idea where this is done for actual values
    DO jt=1,mwiso
       DO jl=1,kproma
          pwiso_2dh(i_runoff)%ptr(jl,jt,jrow) = 0.0_dp
          pwiso_2dh(i_snmel)%ptr(jl,jt,jrow)  = 0.0_dp
          pwiso_2dh(i_apmegl)%ptr(jl,jt,jrow) = 0.0_dp
          pwiso_2dh(i_drain)%ptr(jl,jt,jrow)  = 0.0_dp
          pwiso_2dh(i_snacl)%ptr(jl,jt,jrow)  = 0.0_dp
          pwiso_2dh(i_rogl)%ptr(jl,jt,jrow)   = 0.0_dp
       END DO
    END DO

    ! create the tte value minus the tendency from surf
#ifndef MESSYTENDENCY
! op_pj_20170904: qqq see above
#else
    DO jk = 1, nlev
       tte_ms(1:kproma,jk) = &
            tte_3d(1:kproma,jk,jrow) - mtend_tte_surf(1:kproma,jk,jrow)
    ENDDO
#endif

    ! pwl, pws and psn get the value from vdiff times to start with;
    ! the rest is reconstructed
    pwl(1:kproma) = pwl_surf(1:kproma)
    psn(1:kproma) = psn_surf(1:kproma)
    pws(1:kproma) = pws_surf(1:kproma)
    ! and the same for ptsl
    ptsl(1:kproma) = tsl_surf(1:kproma)

    ! set these to zero (like in ioinitial) for first time step, because
    IF (lstart) THEN
       ! they are set in convect and cloud
       ! (and these dont perform in first time step)
       pwiso_2dh(i_rsfc)%ptr(1:kproma,:,jrow) = 0.0_dp
       pwiso_2dh(i_ssfc)%ptr(1:kproma,:,jrow) = 0.0_dp
       pwiso_2dh(i_rsfl)%ptr(1:kproma,:,jrow) = 0.0_dp
       pwiso_2dh(i_ssfl)%ptr(1:kproma,:,jrow) = 0.0_dp

       ! they are set only in surf
       pwiso_2dh(i_snc)%ptr(1:kproma,:,jrow) = 0.0_dp
       pwiso_2dh(i_gld)%ptr(1:kproma,:,jrow) = 0.0_dp
    ENDIF

    psnc(1:kproma) = snc_surf(1:kproma,jrow)
    pgld(1:kproma) = gld_surf(1:kproma,jrow)

    ! -----------------------------------
    ! H2OISO-SURF------------------------
    ! These subroutines are for calculating several
    ! variables from surf, which are needed for vdiff, cloud and convect
    DO jt=1,mwiso
       DO jl=1,kproma
          zwisosnmel(jl,jt)=0.0_dp
          zwisosncmelt(jl,jt)=0.0_dp
          zwisomlres(jl,jt)=0.0_dp
          pwiso_2dh(i_alac)%ptr(jl,jt,jrow)=0.0_dp
          zwisorogl(jl,jt)=0.0_dp
       ENDDO
    ENDDO
    ! H2OISO-SURF------------------------
    ! This subroutine is to calculate several variables needed 
    ! furtheron in surf
    CALL h2oiso_surf_convflux( kproma, delta_time, psn(1:kproma)   &
         , cvs(1:kproma,jrow), cvw(1:kproma,jrow)                  &
         , pwiso_2dh(i_evapot)%ptr(1:kproma,:,jrow)                &
         , pwl(1:kproma), pwiso_2dh(i_evapl)%ptr(1:kproma,:,jrow)  &
         , zwisoevsnd(1:kproma,:), zwisoraind(1:kproma,:)          &
         , zwisosnowd(1:kproma,:), zwisoevwld(1:kproma,:) &
         , zwisoevwsd(1:kproma,:)                  &
         , pwiso_2dh(i_rsfl)%ptr(1:kproma,:,jrow)  &
         , pwiso_2dh(i_rsfc)%ptr(1:kproma,:,jrow)  &
         , pwiso_2dh(i_ssfl)%ptr(1:kproma,:,jrow)  &
         , pwiso_2dh(i_ssfc)%ptr(1:kproma,:,jrow)  &
         , pwiso_2dh(i_sn)%ptr(1:kproma,:,jrow)    &
         , pwiso_2dh(i_wl)%ptr(1:kproma,:,jrow)    &
         , zwisoevttd(1:kproma,:) )

    ! H2OISO-SURF------------------------
    ! This subroutine is to calculate some isotopic rates of snow changes
    ! in the canopy. (interception of snowfall, sublimation
    ! , melting, unloading due to wind) for further blocks
    ! *zwisosnowd*, *zwisosncmelt*, *zwisoevsnd*  
    CALL h2oiso_surf_snowc( kproma, nlev, rhoh2o, loglac_2d(1:kproma,jrow)  &
         , cwlmax, cvinter, delta_time                                      &
         , zwisosnowd(1:kproma,:), zwisoevsnd(1:kproma,:), zevsnd(1:kproma) &
         , zwisosncmelt(1:kproma,:), psnc(1:kproma), zsnowd(1:kproma)       &
         , tm1(1:kproma,:,jrow), tte_ms(1:kproma,:)                         &
         , ssfl_2d(1:kproma,jrow), ssfc_2d(1:kproma,jrow)                   &
         , loland_2d(1:kproma,jrow), evapot_2d(1:kproma,jrow)               &
         , cvs(1:kproma,jrow), u10(1:kproma,jrow), v10(1:kproma,jrow)       &
         , vlt(1:kproma,jrow), psn(1:kproma), zsn_tmp1(1:kproma)            &
         , zsn_tmp2(1:kproma)                                               &
         , pwiso_2dh(i_snc)%ptr(1:kproma,:,jrow)                            &
         , zsncmelt(1:kproma), pwlmx_tmp(1:kproma), time_step_len           &
         , zwisosn(1:kproma,:) )

    !H2OISO-SURF------------------------
    !This subroutine is to calculate isotopic rates of snowfall 
    !and sublimation on land *zwisoevsnd* *pwisosn*
    CALL h2oiso_surf_snsulw( kproma                                       &
         , pwiso_2dh(i_sn)%ptr(1:kproma,:,jrow), zwisosnowd(1:kproma,:)   &
         , zwisoevsnd(1:kproma,:) , zwisoevwsd(1:kproma,:)                &
         , pwiso_2dh(i_gld)%ptr(1:kproma,:,jrow)                          &
         , loland_2d(1:kproma,jrow), loglac_2d(1:kproma,jrow)             &
         , zsn_tmp1(1:kproma)                                             &
         , pwiso_2dh(i_alac)%ptr(1:kproma,:,jrow), zwisoevttd(1:kproma,:) &
         , zwisorogl(1:kproma,:), zwisoraind(1:kproma,:)                )

    ! H2OISO-SURF------------------------
    ! This subroutine is to calculate isotopic rates of snow
    ! and glacier melt
    CALL h2oiso_surf_snglme( kproma, lstart, rhoh2o                        &
         , pwiso_2dh(i_sn)%ptr(1:kproma,:,jrow)                            &
         , pwiso_2dh(i_gld)%ptr(1:kproma,:,jrow)                           &
         , loland_2d(1:kproma,jrow), loglac_2d(1:kproma,jrow)              &
         , pgld(1:kproma), ptsl(1:kproma), psn(1:kproma)                   &
         , zwisosnmel(1:kproma,:), zsnmel(1:kproma), zsn_tmp2(1:kproma)    &
         , grndcapc_surf(1:kproma), zsnowd(1:kproma), zevsnd(1:kproma)     &
         , zwisorogl(1:kproma,:)                                           )

    ! H2OISO-SURF------------------------
    ! This subroutine is to calculate isotopic rates of snow
    ! budget and melt water
    CALL h2oiso_surf_snbmew( kproma                                       &
         , pwiso_2dh(i_wl)%ptr(1:kproma,:,jrow), zwisomlres(1:kproma,:)   &
         , zwisosnmel(1:kproma,:), zwisosncmelt(1:kproma,:)               &
         , loland_2d(1:kproma,jrow), loglac_2d(1:kproma, jrow)            &
         , zsncmelt(1:kproma)                                             &
         , pwl(1:kproma), pwlmx_tmp(1:kproma)                             &
         , zsnmlt(1:kproma), zsnmel(1:kproma)                             &
         , zwisosn(1:kproma,:)                                            )

    zdtime = delta_time
    DO jl=1,kproma
       zraind_tmp(jl)=(rsfl_2d(jl,jrow)+rsfc_2d(jl,jrow)) *zdtime/rhoh2o
    ENDDO

    ! H2OISO-SURF------------------------
    ! This subroutine is to calculate isotopic rates of the
    ! skin reservoir (vegetation and bare soil)
    CALL h2oiso_surf_skinres( kproma, delta_time, rhoh2o                 &
         , zwisoraind(1:kproma,:), pwiso_2dh(i_wl)%ptr(1:kproma,:,jrow)  &
         , zwisoevwld(1:kproma,:), zwisoevwsd(1:kproma,:)                &
         , loland_2d(1:kproma,jrow), loglac_2d(1:kproma,jrow)            &
         , cvs(1:kproma,jrow), cvw(1:kproma,jrow), evapot_2d(1:kproma,jrow)&
         , pwlmx_tmp(1:kproma)           &
         , pwl(1:kproma), cvinter, zraind_tmp(1:kproma)                  &
         , zmprcp2_tmp(1:kproma)                                         )

    ! H2OISO-SURF------------------------
    ! This subroutine is to calculate isotopic rates of the
    ! soil reservoir (vegetation and bare soil)
    ! and some runoff values
    CALL h2oiso_surf_soilres( kproma, jpgrnd                                &
         , zwisomlres(1:kproma,:), zwisoraind(1:kproma,:)                   &
         , zwisoevwsd(1:kproma,:), pwiso_2dh(i_ws)%ptr(1:kproma,:,jrow)     &
         , evwsd(1:kproma,jrow), zsnmlt(1:kproma), wsmx(1:kproma,jrow)      &
         , loland_2d(1:kproma,jrow), loglac_2d(1:kproma,jrow)               &
         , zraind_tmp(1:kproma)         &
         , delta_time, tsoil(1:kproma,1:jpgrnd,jrow)                        &
         , zmprcp2_tmp(1:kproma), ngl, orostd(1:kproma,jrow), pws(1:kproma) &
         , zwisoros(1:kproma,:), zwisodrain(1:kproma,:) )

    ! H2OISO-SURF------------------------
    ! This subroutine is to correct the fluxes for minor water pools
    CALL h2oiso_surf_corrmp( kproma                                       &
         , pwl(1:kproma), pws(1:kproma), psn(1:kproma)                    &
         , psnc(1:kproma), pgld(1:kproma)                                 &
         , pwiso_2dh(i_wl)%ptr(1:kproma,:,jrow)                           &
         , pwiso_2dh(i_ws)%ptr(1:kproma,:,jrow)                           &
         , pwiso_2dh(i_sn)%ptr(1:kproma,:,jrow)                           &
         , pwiso_2dh(i_snc)%ptr(1:kproma,:,jrow)                          &
         , pwiso_2dh(i_gld)%ptr(1:kproma,:,jrow)                          &
         , slm(1:kproma,jrow), zwisoros(1:kproma,:)                       &
         , zwisosnmel(1:kproma,:), pwiso_2dh(i_alac)%ptr(1:kproma,:,jrow) &
         , zwisodrain(1:kproma,:), zwisosn(1:kproma,:)                    &
         , zwisorogl(1:kproma,:), pwiso_2dh(i_runoff)%ptr(1:kproma,:,jrow) &
         , pwiso_2dh(i_snmel)%ptr(1:kproma,:,jrow)  &
         , pwiso_2dh(i_apmegl)%ptr(1:kproma,:,jrow) &
         , pwiso_2dh(i_drain)%ptr(1:kproma,:,jrow)  &
         , pwiso_2dh(i_snacl)%ptr(1:kproma,:,jrow)  &
         , pwiso_2dh(i_rogl)%ptr(1:kproma,:,jrow) )

! op_pj_20200924: moved to write_output
!!$    ! Devide by timstep, dont know where this is done for actual values
!!$    if (lstart)then
!!$       zdtime = time_step_len
!!$    else
!!$       zdtime = time_step_len / 2.0_dp
!!$    endif
!!$
!!$    pwiso_2dh(i_runoff)%ptr(1:kproma,:,jrow) = pwiso_2dh(i_runoff)%ptr(1:kproma,:,jrow) / zdtime
!!$    pwiso_2dh(i_snmel)%ptr(1:kproma,:,jrow)  = pwiso_2dh(i_snmel)%ptr(1:kproma,:,jrow)  / zdtime
!!$    pwiso_2dh(i_apmegl)%ptr(1:kproma,:,jrow) = pwiso_2dh(i_apmegl)%ptr(1:kproma,:,jrow) / zdtime
!!$    pwiso_2dh(i_drain)%ptr(1:kproma,:,jrow)  = pwiso_2dh(i_drain)%ptr(1:kproma,:,jrow)  / zdtime
!!$    pwiso_2dh(i_snacl)%ptr(1:kproma,:,jrow)  = pwiso_2dh(i_snacl)%ptr(1:kproma,:,jrow)  / zdtime
!!$    pwiso_2dh(i_rogl)%ptr(1:kproma,:,jrow)   = pwiso_2dh(i_rogl)%ptr(1:kproma,:,jrow)   / zdtime
!!$
!!$    pwiso_2dh(i_evap)%ptr(1:kproma,:,jrow)    = pwiso_2dh(i_evap)%ptr(1:kproma,:,jrow)    / zdtime
!!$    pwiso_2dh(i_evapl)%ptr(1:kproma,:,jrow)   = pwiso_2dh(i_evapl)%ptr(1:kproma,:,jrow)   / zdtime
!!$    pwiso_2dh(i_evapw)%ptr(1:kproma,:,jrow)   = pwiso_2dh(i_evapw)%ptr(1:kproma,:,jrow)   / zdtime
!!$    pwiso_2dh(i_evapi)%ptr(1:kproma,:,jrow)   = pwiso_2dh(i_evapi)%ptr(1:kproma,:,jrow)   / zdtime
!!$    pwiso_2dh(i_evaplac)%ptr(1:kproma,:,jrow) = pwiso_2dh(i_evaplac)%ptr(1:kproma,:,jrow) / zdtime
!!$    pwiso_2dh(i_evapwac)%ptr(1:kproma,:,jrow) = pwiso_2dh(i_evapwac)%ptr(1:kproma,:,jrow) / zdtime
!!$    pwiso_2dh(i_evapiac)%ptr(1:kproma,:,jrow) = pwiso_2dh(i_evapiac)%ptr(1:kproma,:,jrow) / zdtime

  END SUBROUTINE h2oiso_mixlo
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_physc(flag)

    ! ------------------------------------------------------------------
    ! This subroutine is called within the time loop.
    ! It constitutes the main entry point for the chemical
    ! influence of CH4/CH3D oxidation on HHO and HDO
    ! and some diagnostics.
    ! Here, only the current vector of the grid-point-fields is
    ! accessible.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)

    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    USE messy_main_data_bi,       ONLY:  &
                       qm1, xlm1, xim1, qte_3d, xlte_3d, xite_3d
    USE messy_main_tracer_mem_bi, ONLY: xtm1, xtte
    USE messy_main_timer,         ONLY: time_step_len
    USE messy_main_constants_mem, ONLY: FLAGGED_BAD
    
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_physc'
    INTEGER                     :: status
    REAL(DP):: pwisoqm1(nproma,nlev,mwiso)
    REAL(DP):: pwisoxlm1(nproma,nlev,mwiso)
    REAL(DP):: pwisoxim1(nproma,nlev,mwiso)

    ! for the numerical error
    REAL(DP):: pwisoqte_numer(nproma,nlev)
    REAL(DP):: pwisoxlte_numer(nproma,nlev)
    REAL(DP):: pwisoxite_numer(nproma,nlev)
    REAL(DP), PARAMETER :: myepsilon = EPSILON(0.0_dp) * 10.0_dp

    REAL(DP):: ztmst
    REAL(DP):: mmr
    INTEGER::  idt, jl, jk, jwiso

    LOGICAL, DIMENSION(nproma,nlev) :: val_psc

#ifdef MESSYTENDENCY
    CHARACTER(len=32), PARAMETER  ::  str_ch4     =  'ch4'
    CHARACTER(len=32), PARAMETER  ::  str_mecca   =  'mecca'
#endif

    SELECT CASE (flag)
    CASE(1)

       ztmst = time_step_len
       ! add tendencies from msbm
       IF (l_msbm) THEN 
          IF (lreg_msbm) THEN

             ! Note: the non-MESSYTENDENCY branch (even though it is useless)
             !       is required, because without it, the code will NOT
             !       EVEN compile without MESSYTENDENCY.
             DO jwiso = 1, mwiso
#ifndef MESSYTENDENCY

#else
                CALL mtend_get_start_l(idiso(i_vap,jwiso) &
                     , v0=pwisoqm1(1:kproma,:,jwiso))
                CALL mtend_get_start_l(idiso(i_liq,jwiso) &
                     , v0=pwisoxlm1(1:kproma,:,jwiso))
                CALL mtend_get_start_l(idiso(i_ice,jwiso) &
                     , v0=pwisoxim1(1:kproma,:,jwiso))
#endif
             END DO

             WHERE (NINT(flt_stratreg(1:kproma,:,jrow))==1) 
                val_psc(1:kproma,:)=.true.
             ELSEWHERE
                val_psc(1:kproma,:)=.false.
             END WHERE

             DO jwiso = 2, mwiso

                ! reset to zero, otherwise tendencies will be applied
                ! to other isotopic tracers as well
                pwisoqte_numer(:,:) = 0.0_dp
                pwisoxlte_numer(:,:) = 0.0_dp
                pwisoxite_numer(:,:) = 0.0_dp    

                mmr = molarmasses(jwiso)/molarmasses(i_HHO)
#ifndef MESSYTENDENCY
! op_pj_20170904: qqq see above
#else
                WHERE(val_psc(1:kproma,:))
                   WHERE ( pwisoqm1(1:kproma,:,i_HHO) .lt. cwisomin .or. &
                        pwisoqm1(1:kproma,:,jwiso) .lt. 0.0_dp )
                      pwisoqte_numer(1:kproma,:) =  &
                           ( ( tnat(jwiso) / ( 1.0_dp + tnat(jwiso) ) ) &
                           * mmr * &
                           ( pwisoqm1(1:kproma,:,i_HHO) + &
                           mtend_qte_msbm(1:kproma,:,jrow) * ztmst ) &
                           - pwisoqm1(1:kproma,:,jwiso) )  / ztmst
                   ELSEWHERE
                      pwisoqte_numer(1:kproma,:) = &
                           mtend_qte_msbm(1:kproma,:,jrow) * &
                           MIN(( pwisoqm1(1:kproma,:,jwiso) &
                           / pwisoqm1(1:kproma,:,i_HHO) ), 1.0_dp)
                      ! for the unintended case that
                      ! pwisoqm1(jwiso).gt.pwisoqm1(i_HHO)
                   ENDWHERE

                   WHERE ( pwisoxlm1(1:kproma,:,i_HHO) .lt. cwisomin .or. &
                        pwisoxlm1(1:kproma,:,jwiso) .lt. 0.0_dp )
                      pwisoxlte_numer(1:kproma,:) = &
                           ( ( tnat(jwiso) / ( 1.0_dp + tnat(jwiso) ) ) &
                           * mmr * &
                           ( pwisoxlm1(1:kproma,:,i_HHO) + &
                           mtend_xlte_msbm(1:kproma,:,jrow) * ztmst ) &
                           - pwisoxlm1(1:kproma,:,jwiso) )  / ztmst
                   ELSEWHERE
                      pwisoxlte_numer(1:kproma,:) = &
                           mtend_xlte_msbm(1:kproma,:,jrow) * &
                           MIN(( pwisoxlm1(1:kproma,:,jwiso) &
                           / pwisoxlm1(1:kproma,:,i_HHO) ), 1.0_dp)
                      ! for the unintended case that
                      ! pwisoxlm1(jwiso).gt.pwisoxlm1(i_HHO)
                   ENDWHERE

                   WHERE ( pwisoxim1(1:kproma,:,i_HHO) .lt. cwisomin .or. &
                        pwisoxim1(1:kproma,:,jwiso) .lt. 0.0_dp )
                      pwisoxite_numer(1:kproma,:) = &
                           ( ( tnat(jwiso) / ( 1.0_dp + tnat(jwiso) ) ) &
                           * mmr * &
                           ( pwisoxim1(1:kproma,:,i_HHO) + &
                           mtend_xite_msbm(1:kproma,:,jrow) * ztmst ) &
                           - pwisoxim1(1:kproma,:,jwiso) )  / ztmst
                   ELSEWHERE
                      pwisoxite_numer(1:kproma,:) = &
                           mtend_xite_msbm(1:kproma,:,jrow) * &
                           MIN(( pwisoxim1(1:kproma,:,jwiso) &
                           / pwisoxim1(1:kproma,:,i_HHO) ), 1.0_dp)
                      ! for the unintended case that
                      ! pwisoxim1(jwiso).gt.pwisoxim1(i_HHO)
                   ENDWHERE
                END WHERE
#endif

#ifndef MESSYTENDENCY
! op_pj_20170904: qqq see above
#else
                CALL mtend_add_l(my_handle_msbm, idiso(i_vap,jwiso), &
                     px=pwisoqte_numer(1:kproma,:))
                CALL mtend_add_l(my_handle_msbm, idiso(i_liq,jwiso), &
                     px=pwisoxlte_numer(1:kproma,:))
                CALL mtend_add_l(my_handle_msbm, idiso(i_ice,jwiso), &
                     px=pwisoxite_numer(1:kproma,:))
#endif
             END DO

#ifndef MESSYTENDENCY
             ! This does not work without MESSYTENDENCY and as it is required
             !  for h2oiso anyway there is no meaning in this branch.

#else
             CALL mtend_add_l(my_handle_msbm, idiso(i_vap,i_HHO), &
                  px=mtend_qte_msbm(1:kproma,:,jrow))
             CALL mtend_add_l(my_handle_msbm, idiso(i_liq,i_HHO), &
                  px=mtend_xlte_msbm(1:kproma,:,jrow))
             CALL mtend_add_l(my_handle_msbm, idiso(i_ice,i_HHO), &
                  px=mtend_xite_msbm(1:kproma,:,jrow))
#endif 
          ENDIF
       END IF

    CASE(2) 

       ! Add the tendency of the ch4 oxidation to HHO to resolve numerical
       ! errors in the following synchronization
       IF (lreg_ch4) THEN
#ifndef MESSYTENDENCY
          ! This does not work without MESSYTENDENCY and as it is required for 
          ! h2oiso anyway there is no meaning in this branch
          
#else
          CALL mtend_add_l(my_handle_chem, idiso(i_vap,i_HHO), &
               px=mtend_qte_ch4(1:kproma,:,jrow))
#endif
       ENDIF
    
       IF (lreg_mecca) THEN
#ifndef MESSYTENDENCY
          ! This does not work without MESSYTENDENCY ...

#else
          CALL mtend_add_l(my_handle_chem, idiso(i_vap,i_HHO), &
               px=mtend_qte_mecca(1:kproma,:,jrow) )
#endif 
       ENDIF

       !===============================================
       ! The next section is to take account for the numerical errors,
       ! that are being made for HHO through the fractionation thing
       ! see R. Eichingers Dissertation
       ! These will be additional channel objects in MESSYTENDENCY
       !===============================================

       DO jk = 1, nlev
          pwisoqte_numer(1:kproma,nlev)  = 0.0_dp
          pwisoxlte_numer(1:kproma,nlev) = 0.0_dp
          pwisoxite_numer(1:kproma,nlev) = 0.0_dp
       ENDDO

       ztmst = time_step_len

       DO jk = 1, nlev
          DO jl = 1, kproma
             IF(abs(isotracte(i_vap)%ptr(jl,jk,i_HDO,jrow)).le.myepsilon)THEN
                pwiso_3d(i_numq)%ptr(jl,jk,i_HHO,jrow)  = FLAGGED_BAD
             ELSE
                pwiso_3d(i_numq)%ptr(jl,jk,i_HHO,jrow)  = &
                     (abs(qte_3d(jl,jk,jrow)  - &
                     isotracte(i_vap)%ptr(jl,jk,i_HHO,jrow)) &
                     / abs(isotracte(i_vap)%ptr(jl,jk,i_HDO,jrow))) * 100.0_dp
             ENDIF
             IF(abs(isotracte(i_liq)%ptr(jl,jk,i_HDO,jrow)).le.myepsilon)THEN
                pwiso_3d(i_numxl)%ptr(jl,jk,i_HHO,jrow) =  FLAGGED_BAD
             ELSE
                pwiso_3d(i_numxl)%ptr(jl,jk,i_HHO,jrow) = &
                     (abs(xlte_3d(jl,jk,jrow) - &
                     isotracte(i_liq)%ptr(jl,jk,i_HHO,jrow)) &
                     / abs(isotracte(i_liq)%ptr(jl,jk,i_HDO,jrow))) * 100.0_dp
             ENDIF
             IF(abs(isotracte(i_ice)%ptr(jl,jk,i_HDO,jrow)).le.myepsilon)THEN
                pwiso_3d(i_numxi)%ptr(jl,jk,i_HHO,jrow) =  FLAGGED_BAD
             ELSE
                pwiso_3d(i_numxi)%ptr(jl,jk,i_HHO,jrow) = &
                     (abs(xite_3d(jl,jk,jrow) - &
                     isotracte(i_ice)%ptr(jl,jk,i_HHO,jrow)) &
                     / abs(isotracte(i_ice)%ptr(jl,jk,i_HDO,jrow))) * 100.0_dp
             ENDIF
          ENDDO
       ENDDO
       
       !----------------------------------------------------------------
       ! Account for the numerical errors that are being made
       ! for tendency of H2OISOHHOvap compared to qte.
       ! The following part was extended to work without a loop.
#ifndef MESSYTENDENCY
       isotracte(i_vap)%ptr(1:kproma,:,i_HHO,jrow) = &
            isotracte(i_vap)%ptr(1:kproma,:,i_HHO,jrow) + &
            (qte_3d(1:kproma,:,jrow)  - &
            isotracte(i_vap)%ptr(1:kproma,:,i_HHO,jrow))
       isotracte(i_liq)%ptr(1:kproma,:,i_HHO,jrow) = &
            isotracte(i_liq)%ptr(1:kproma,:,i_HHO,jrow) + &
            (xlte_3d(1:kproma,:,jrow) - &
            isotracte(i_liq)%ptr(1:kproma,:,i_HHO,jrow))
       isotracte(i_ice)%ptr(1:kproma,:,i_HHO,jrow) = &
            isotracte(i_ice)%ptr(1:kproma,:,i_HHO,jrow) + &
            (xite_3d(1:kproma,:,jrow) - &
            isotracte(i_ice)%ptr(1:kproma,:,i_HHO,jrow))
#else
       pwisoqte_numer(1:kproma,:)  = &
            qte_3d(1:kproma,:,jrow)  - &
            isotracte(i_vap)%ptr(1:kproma,:,i_HHO,jrow)
       pwisoxlte_numer(1:kproma,:) = &
            xlte_3d(1:kproma,:,jrow) - &
            isotracte(i_liq)%ptr(1:kproma,:,i_HHO,jrow)
       pwisoxite_numer(1:kproma,:) = &
            xite_3d(1:kproma,:,jrow) - &
            isotracte(i_ice)%ptr(1:kproma,:,i_HHO,jrow)
       
       CALL mtend_add_l(my_handle_numer, idiso(i_vap,i_HHO), &
            px=pwisoqte_numer(1:kproma,:) )
       CALL mtend_add_l(my_handle_numer, idiso(i_liq,i_HHO), &
            px=pwisoxlte_numer(1:kproma,:) )
       CALL mtend_add_l(my_handle_numer, idiso(i_ice,i_HHO), &
            px=pwisoxite_numer(1:kproma,:) )
#endif

       !--------------------------------------------------------------
       ! This block sets the isotopic tendencies so that the p1 value of it 
       ! gets a delta value of VSMOW, if the p1 values are very small.
       ! The name '_numer' for this operation does not fit, it is only chosen
       ! because it is available for the isotopes from the numeric correction
       ! block
       ! Beware, this is not done for jwiso = 1, which has to be HH16O
       DO jwiso = 2, mwiso

          ! reset to zero, otherwise tendencies will be applied to other 
          ! isotopic tracers as well
          pwisoqte_numer(:,:) = 0.0_dp
          pwisoxlte_numer(:,:) = 0.0_dp
          pwisoxite_numer(:,:) = 0.0_dp    

          DO jl = 1, kproma
             DO jk = 1, nlev
                ! calculate new tendency for isotopic vapour
                IF( (isotracm1(i_vap)%ptr(jl,jk,jwiso,jrow) + &
                     isotracte(i_vap)%ptr(jl,jk,jwiso,jrow) * ztmst) &
                     .lt. cwisomin .or. &
                     (isotracm1(i_vap)%ptr(jl,jk,i_HHO,jrow) + &
                     isotracte(i_vap)%ptr(jl,jk,i_HHO,jrow) * ztmst) &
                     .lt. cwisomin) THEN

                   CALL h2oiso_constdelta_tendency(&
                        pwisoqte_numer(jl,jk), tnat(jwiso), &
                        qm1(jl,jk,jrow) + qte_3d(jl,jk,jrow)*ztmst, &
                        isotracm1(i_vap)%ptr(jl,jk,jwiso,jrow)+ &
                        isotracte(i_vap)%ptr(jl,jk,jwiso,jrow)*ztmst, &
                        molarmasses(i_HHO), molarmasses(jwiso), &
                        ztmst)
                ENDIF

                ! calculate new tendency for isotopic liquid water
                IF( (isotracm1(i_liq)%ptr(jl,jk,jwiso,jrow) + &
                     isotracte(i_liq)%ptr(jl,jk,jwiso,jrow) * ztmst) &
                     .lt. cwisomin .or. &
                     (isotracm1(i_liq)%ptr(jl,jk,i_HHO,jrow) + &
                     isotracte(i_liq)%ptr(jl,jk,i_HHO,jrow) * ztmst) &
                     .lt. cwisomin) THEN

                   CALL h2oiso_constdelta_tendency(&
                     pwisoxlte_numer(jl,jk), tnat(jwiso), &
                     xlm1(jl,jk,jrow) + xlte_3d(jl,jk,jrow)*ztmst, &
                     isotracm1(i_liq)%ptr(jl,jk,jwiso,jrow)+ &
                     isotracte(i_liq)%ptr(jl,jk,jwiso,jrow)*ztmst, &
                     molarmasses(i_HHO), molarmasses(jwiso), &
                     ztmst)
                ENDIF

                ! calculate new tendency for isotopic ice water
                IF( (isotracm1(i_ice)%ptr(jl,jk,jwiso,jrow) + &
                     isotracte(i_ice)%ptr(jl,jk,jwiso,jrow) * ztmst) &
                     .lt. cwisomin .or. &
                     (isotracm1(i_ice)%ptr(jl,jk,i_HHO,jrow) + &
                     isotracte(i_ice)%ptr(jl,jk,i_HHO,jrow) * ztmst) &
                     .lt. cwisomin) THEN

                   CALL h2oiso_constdelta_tendency(&
                        pwisoxite_numer(jl,jk), tnat(jwiso), &
                        xim1(jl,jk,jrow) + xite_3d(jl,jk,jrow)*ztmst, &
                        isotracm1(i_ice)%ptr(jl,jk,jwiso,jrow)+ &
                        isotracte(i_ice)%ptr(jl,jk,jwiso,jrow)*ztmst, &
                        molarmasses(i_HHO), molarmasses(jwiso), &
                        ztmst)
                ENDIF
             ENDDO
          ENDDO

          ! Now add them up
#ifndef MESSYTENDENCY
          isotracte(i_vap)%ptr(1:kproma,:,jwiso,jrow) &
               = isotracte(i_vap)%ptr(1:kproma,:,jwiso,jrow) + &
               pwisoqte_numer(1:kproma,:)
          isotracte(i_liq)%ptr(1:kproma,:,jwiso,jrow) &
               = isotracte(i_liq)%ptr(1:kproma,:,jwiso,jrow) + &
               pwisoxlte_numer(1:kproma,:)
          isotracte(i_ice)%ptr(1:kproma,:,jwiso,jrow) &
               = isotracte(i_ice)%ptr(1:kproma,:,jwiso,jrow) + &
               pwisoxite_numer(1:kproma,:)
#else
          CALL mtend_add_l(my_handle_numer, idiso(i_vap,jwiso), &
               px=pwisoqte_numer(1:kproma,:))
          CALL mtend_add_l(my_handle_numer, idiso(i_liq,jwiso), &
               px=pwisoxlte_numer(1:kproma,:))
          CALL mtend_add_l(my_handle_numer, idiso(i_ice,jwiso), &
               px=pwisoxite_numer(1:kproma,:))
#endif
       ENDDO
       
       !-------------------------------------------------------------
       ! THIS DOES NOT REALLY WORK!! Delta values used so far are calculated 
       ! during post-processing.
       ! Calculate the delta values for easier postprocessing
       ! to avoid nan, set the delta value to -1.0e+34(=nan, in ferret) 
       ! if HHO is smaller than myepsilon.
       DO jwiso = 1, mwiso
          IF (jwiso.gt.1) THEN
             
             DO jl = 1, kproma
                DO jk = 1, nlev
                   
                   ! delta value for water vapour
                   IF(xtm1(jl,jk,idiso(i_vap,i_HHO),jrow).lt.myepsilon &
                        .or.xtm1(jl,jk,idiso(i_vap,jwiso),jrow).lt.myepsilon) THEN
                      pwiso_3d(i_dvap)%ptr(jl,jk,jwiso,jrow) = FLAGGED_BAD
                   ELSE
                      pwiso_3d(i_dvap)%ptr(jl,jk,jwiso,jrow) = &
                           1000.0_dp * ((xtm1(jl,jk,idiso(i_vap,jwiso),jrow) &
                           / xtm1(jl,jk,idiso(i_vap,i_HHO),jrow)) &
                           -tnat(jwiso)) / tnat(jwiso)
                   ENDIF
                   
                   ! delta value for liquid water
                   IF(xtm1(jl,jk,idiso(i_liq,i_HHO),jrow).lt.myepsilon &
                        .or.xtm1(jl,jk,idiso(i_liq,jwiso),jrow).lt.myepsilon) THEN
                      pwiso_3d(i_dliq)%ptr(jl,jk,jwiso,jrow) = FLAGGED_BAD
                   ELSE
                      pwiso_3d(i_dliq)%ptr(jl,jk,jwiso,jrow) = &
                           1000.0_dp * ((xtm1(jl,jk,idiso(i_liq,jwiso),jrow) &
                           / xtm1(jl,jk,idiso(i_liq,i_HHO),jrow)) &
                           - tnat(jwiso)) / tnat(jwiso)
                   ENDIF
                   
                   ! delta value for ice water
                   IF(xtm1(jl,jk,idiso(i_ice,i_HHO),jrow).lt.myepsilon &
                        .or.xtm1(jl,jk,idiso(i_ice,jwiso),jrow).lt.myepsilon) THEN
                      pwiso_3d(i_dice)%ptr(jl,jk,jwiso,jrow) = FLAGGED_BAD
                   ELSE
                      pwiso_3d(i_dice)%ptr(jl,jk,jwiso,jrow) = &
                           1000.0_dp * ((xtm1(jl,jk,idiso(i_ice,jwiso),jrow) &
                           / xtm1(jl,jk,idiso(i_ice,i_HHO),jrow)) &
                           - tnat(jwiso)) / tnat(jwiso)
                   ENDIF
                   
                   ! delta value for total water
                   IF(xtm1(jl,jk,idiso(i_vap,i_HHO),jrow)+xtm1(jl,jk,idiso(i_liq,i_HHO),jrow)&
                        +xtm1(jl,jk,idiso(i_ice,i_HHO),jrow).lt.myepsilon .or. &
                        xtm1(jl,jk,idiso(i_vap,jwiso),jrow)+xtm1(jl,jk,idiso(i_liq,jwiso),jrow)&
                        +xtm1(jl,jk,idiso(i_ice,jwiso),jrow).lt.myepsilon)THEN
                      pwiso_3d(i_dtot)%ptr(jl,jk,jwiso,jrow) = FLAGGED_BAD
                   ELSE
                      pwiso_3d(i_dtot)%ptr(jl,jk,jwiso,jrow) = &
                           1000.0_dp * ( &
                           ((xtm1(jl,jk,idiso(i_vap,jwiso),jrow)+xtm1(jl,jk,idiso(i_liq,jwiso),jrow)&
                           +xtm1(jl,jk,idiso(i_ice,jwiso),jrow)) / &
                           (xtm1(jl,jk,idiso(i_vap,i_HHO),jrow)+xtm1(jl,jk,idiso(i_liq,i_HHO),jrow)&
                           +xtm1(jl,jk,idiso(i_ice,i_HHO),jrow))) &
                           - tnat(jwiso)) / tnat(jwiso)
                   ENDIF
                ENDDO
                
                ! delta value for precipitated snow
                IF(pwiso_2dh(i_aprs)%ptr(jl,i_HHO,jrow).lt.myepsilon &
                     .or.pwiso_2dh(i_aprs)%ptr(jl,jwiso,jrow).lt.myepsilon) THEN
                   pwiso_2dh(i_daprs)%ptr(jl,jwiso,jrow) = FLAGGED_BAD
                ELSE
                   pwiso_2dh(i_daprs)%ptr(jl,jwiso,jrow) = 1000.0_dp * &
                        ((pwiso_2dh(i_aprs)%ptr(jl,jwiso,jrow) & 
                        / pwiso_2dh(i_aprs)%ptr(jl,i_HHO,jrow)) - tnat(jwiso)) &
                        / tnat(jwiso)
                ENDIF
                
                ! delta value for precipitated large scale rain
                IF(pwiso_2dh(i_aprl)%ptr(jl,i_HHO,jrow).lt.myepsilon &
                     .or.pwiso_2dh(i_aprl)%ptr(jl,jwiso,jrow).lt.myepsilon) THEN
                   pwiso_2dh(i_daprl)%ptr(jl,jwiso,jrow) = FLAGGED_BAD
                ELSE
                   pwiso_2dh(i_daprl)%ptr(jl,jwiso,jrow) = 1000.0_dp * &
                        ((pwiso_2dh(i_aprl)%ptr(jl,jwiso,jrow) & 
                        / pwiso_2dh(i_aprl)%ptr(jl,i_HHO,jrow)) - tnat(jwiso)) &
                        / tnat(jwiso)
                ENDIF
                
                ! delta value for precipitated convective rain
                IF(pwiso_2dh(i_aprc)%ptr(jl,i_HHO,jrow).lt.myepsilon &
                     .or.pwiso_2dh(i_aprc)%ptr(jl,jwiso,jrow).lt.myepsilon) THEN
                   pwiso_2dh(i_daprc)%ptr(jl,jwiso,jrow) = FLAGGED_BAD
                ELSE
                   pwiso_2dh(i_daprc)%ptr(jl,jwiso,jrow) = 1000.0_dp * &
                        ((pwiso_2dh(i_aprc)%ptr(jl,jwiso,jrow) & 
                        / pwiso_2dh(i_aprc)%ptr(jl,i_HHO,jrow)) - tnat(jwiso)) &
                        / tnat(jwiso)
                ENDIF
                
                ! delta value for precipitated total rain
                IF(pwiso_2dh(i_aprc)%ptr(jl,i_HHO,jrow) + &
                     pwiso_2dh(i_aprl)%ptr(jl,i_HHO,jrow).lt.myepsilon .or. &
                     pwiso_2dh(i_aprc)%ptr(jl,jwiso,jrow) + &
                     pwiso_2dh(i_aprl)%ptr(jl,jwiso,jrow).lt.myepsilon) THEN
                   pwiso_2dh(i_daprt)%ptr(jl,jwiso,jrow) = FLAGGED_BAD
                ELSE
                   pwiso_2dh(i_daprt)%ptr(jl,jwiso,jrow) = &
                        1000.0_dp * (((pwiso_2dh(i_aprc)%ptr(jl,jwiso,jrow) + &
                        pwiso_2dh(i_aprl)%ptr(jl,jwiso,jrow)) & 
                        / (pwiso_2dh(i_aprc)%ptr(jl,i_HHO,jrow) + &
                        pwiso_2dh(i_aprl)%ptr(jl,i_HHO,jrow))) - tnat(jwiso)) &
                        / tnat(jwiso)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       
    END SELECT

  END SUBROUTINE h2oiso_physc
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_write_output

    USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma, npromz
    USE messy_main_timer,           ONLY: time_step_len, lstart
    
    IMPLICIT NONE

    ! The diagnostic objects (apr*, *vi, etc.) are internally accumulated
    ! (i.e.: x' = x + dx * dt)
    ! in the same way as it was done in ECHAM5 for the classical diagnostic
    ! objects. In ECHAM5 these "stream elements" have been flagged with
    ! 'laccu = .TRUE.'. Fot these, the de-accumulation is centrally
    ! performed in messy_main_channel_echam5.inc:: reset_accu_stream_elements.
    ! For the H2OISO channel objects this is not possible (since they are
    ! no "stream elements", therefore we need to perform this backward
    ! operation (deaccumultion) here.

    REAL(dp) :: zdtime
    INTEGER  :: jjrow, zkproma
    INTEGER  :: i

    IF (lstart) THEN
       zdtime = time_step_len
    ELSE
       zdtime = time_step_len / 2.0_dp
    END IF
    
    DO jjrow = 1, ngpblks
       IF ( jjrow == ngpblks ) THEN
          zkproma = npromz
       ELSE
          zkproma = nproma
       END IF

       DO i = 1, i_2nmax

          SELECT CASE(i)
          CASE(i_aprc, i_aprs, i_aprl, i_qvi, i_xlvi, i_xivi, &
               i_runoff, i_snmel, i_apmegl, i_drain, i_snacl, i_rogl, &
               i_evap, i_evapl, i_evapw, i_evapi, i_evaplac, i_evapwac, &
               i_evapiac)
             !
             pwiso_2dh(i)%ptr(1:zkproma,:,jjrow) = &
                  pwiso_2dh(i)%ptr(1:zkproma,:,jjrow) / zdtime
             !
          CASE DEFAULT
             ! nothing to do
          END SELECT

       END DO

    END DO
       
  END SUBROUTINE h2oiso_write_output
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE h2oiso_free_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to deallocate the memory, which has
    ! been "manually" allocated in h2oiso_init_memory.
    ! Note: channel object memory must not be deallocated! This is
    !       performed centrally.
    ! ------------------------------------------------------------------

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oiso_free_memory'

    INTEGER                     :: ji

    ! deallocate the pwiso pointer structures
    DO ji = 1, i_2nmax
       DEALLOCATE(pwiso_2dh_5d(ji)%ptr); NULLIFY(pwiso_2dh_5d(ji)%ptr)
    ENDDO
    DO ji = 1, i_3nmax
       DEALLOCATE(pwiso_3d_5d(ji)%ptr); NULLIFY(pwiso_3d_5d(ji)%ptr)
    ENDDO
    DEALLOCATE(pwiso_2dh_5d); NULLIFY(pwiso_2dh_5d)
    DEALLOCATE(pwiso_3d_5d); NULLIFY(pwiso_3d_5d)
    DEALLOCATE(pwiso_2dh); NULLIFY(pwiso_2dh)
    DEALLOCATE(pwiso_3d); NULLIFY(pwiso_3d)

    ! deallocate the isotrac pointer structures
    DEALLOCATE(isotracm1); NULLIFY(isotracm1)
    DEALLOCATE(isotracte); NULLIFY(isotracte)

    ! deallocation of other global variables
    DEALLOCATE(pwl_surf); NULLIFY(pwl_surf)
    DEALLOCATE(psn_surf); NULLIFY(psn_surf)
    DEALLOCATE(pws_surf); NULLIFY(pws_surf)
    DEALLOCATE(tsl_vdiff); NULLIFY(tsl_vdiff)
    DEALLOCATE(tsl_surf); NULLIFY(tsl_surf)
    DEALLOCATE(grndcapc_surf); NULLIFY(grndcapc_surf)

  END SUBROUTINE h2oiso_free_memory
  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE h2oiso_read_nml_cpl(status, iou)

    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    !     ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ l_steady, &
         l_noconvect_dd, l_nocloud_dd    ! sc_re_20150701

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='h2oiso_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR


  END SUBROUTINE h2oiso_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_h2oiso_e5
! **********************************************************************
