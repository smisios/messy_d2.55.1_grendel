#include "messy_main_ppd_bi.inc"

MODULE MESSY_CONVECT_SI


!  This module is an interface for convection parametrization schemes.
!  It collects the neccessary parameters from the basemodel by use statements 
!  and channel objects, reads the namelists to choose which convection 
!  with which additional parameters should be used

!  Author:   H. Tost,      MPICH, March 2004 - September 2007

  ! MESSy
  USE messy_main_channel,        ONLY: t_chaobj_cpl  ! um_ak_20140502
  USE MESSY_MAIN_CONSTANTS_MEM,  ONLY: dp
  USE MESSY_convect,             ONLY: convect_param, modstr, L_LGMC
  USE MESSY_CONVECT_MEM

#ifdef MESSYTENDENCY
  ! tendency budget
   USE messy_main_tendency_bi,   ONLY: mtend_get_handle,               &
                                       mtend_get_start_l,              &
                                       mtend_add_l,                    &
                                       mtend_register,                 &
                                       mtend_set_sqcst_scal,           &
                                       mtend_id_t,                     &
                                       mtend_id_q,                     &
                                       mtend_id_xl,                    &
                                       mtend_id_xi,                    &
                                       mtend_id_u,                     &
                                       mtend_id_v,                     &
                                       mtend_id_tracer
#endif

  IMPLICIT NONE

  SAVE

  INTRINSIC :: INT, MAX, MIN
 
  REAL(dp), POINTER :: pilab(:,:,:)   => NULL()
  REAL(dp), POINTER :: pqtec(:,:,:)   => NULL()
  REAL(dp), POINTER :: zpqtec(:,:)    => NULL() ! mz_pj_20070309
  REAL(dp), POINTER :: pqhfla(:)      => NULL()
  REAL(dp), POINTER :: pqflux(:,:)    => NULL()
  REAL(dp), POINTER :: grmass(:,:,:)  => NULL()
  REAL(dp), POINTER :: grvol(:,:,:)   => NULL()
  REAL(dp), POINTER :: geopoti(:,:,:) => NULL()
  REAL(dp), POINTER :: geopot(:,:,:)  => NULL()
  REAL(dp), POINTER :: omga(:,:)      => NULL()
  REAL(dp), POINTER :: pblh(:,:)      => NULL()
  REAL(dp), POINTER :: counter(:,:)   => NULL()
  REAL(dp), POINTER :: cth(:,:)       => NULL()
  REAL(dp), POINTER :: cbmf(:,:)      => NULL()

  INTEGER,  POINTER :: jlab(:,:)      => NULL()
  REAL(dp), POINTER :: pdiga5(:,:,:)  => NULL()
  REAL(dp), POINTER :: pdiga10(:,:,:) => NULL()
  REAL(dp), POINTER :: pdiga18(:,:,:) => NULL()


  REAL(dp), POINTER :: counter_deep(:,:) => NULL()
  REAL(dp), POINTER :: counter_shal(:,:) => NULL()
  REAL(dp), POINTER :: counter_midl(:,:) => NULL()

  REAL(dp), POINTER :: cth_deep(:,:) => NULL()
  REAL(dp), POINTER :: cth_shal(:,:) => NULL()
  REAL(dp), POINTER :: cth_midl(:,:) => NULL()

  REAL(dp), POINTER :: massfu_deep(:,:,:) => NULL()
  REAL(dp), POINTER :: massfu_shal(:,:,:) => NULL()
  REAL(dp), POINTER :: massfu_midl(:,:,:) => NULL()

  ! mz_jd_20161011+
  REAL(dp), POINTER :: qstdev(:,:,:)               => NULL()
  REAL(dp), POINTER :: qskew(:,:,:)                => NULL()
  REAL(dp), POINTER :: Q2(:,:,:)                   => NULL()
  REAL(dp), POINTER :: qsat(:,:,:)                 => NULL()
  REAL(dp), POINTER :: qp1(:,:,:)                  => NULL()
  ! mz_jd_20161011-

! op_sb_20171018_ moved to SMCL
!!$  ! mim_sb_20091207+
!!$  LOGICAL :: L_LGMC = .FALSE.      ! submodel LGMC is ON
!!$  ! mim_sb_20091207-
  ! um_ak_20140502+
  TYPE(t_chaobj_cpl) :: pbl_height
  ! um_ak_20140502-
#ifdef MESSYTENDENCY
  ! variable for tendency budget
  integer                         :: my_handle
#endif

CONTAINS

!=================================================================================

  SUBROUTINE  convect_initialize

    ! Convection MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! 
    ! Author: H. Tost, MPICH, March 2004

    ! ECHAM5/MESSy

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast !, finish
    ! op_mm_20140110 finish -> error_bi
    USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi, error_bi
    ! MESSy
    USE messy_main_tools,      ONLY: find_next_free_unit
    USE MESSY_convect,         ONLY: lconvection, init_convection, &
                                     altconv, evap_sub, ltransport
    USE MESSY_CONVECT_TIEDTKE, ONLY: lmfpen, lmfscv, lmfmid, lmfdd, lmfdudv, &
                                     lpos_def, &
! fb_mk_20120116+
                                     rset_cmfctop, rset_cprcon, rset_entrscv, &
! fb_mk_20120116-
! fb_mk_20140210+
                                     rset_entrpen, rset_entrmid
! fb_mk_20140210-
    USE MESSY_CONVECT_MEM,     ONLY: ODEEP, OSHAL, ODOWN, OREFRESH_ALL, &
                                     OSETTADJ, OUVTRANS, OCHTRANS,      &
                                     KENSM, KICE, PTADJD, PTADJS,       &
                                     l_lgmc_diag ! op_pj_20091207
    USE MESSY_CONVECT_ECMWF_PARAM,    ONLY: PEN => LMFPEN, SHAL => LMFSCV,     &
                                            MID => LMFMID, DOWN => LMFDD,      &
                                            FRIC => LMFDUDV,                   &
                                            LMFTRAC, LEPCLD, LMFSCL_WSTAR
    USE MESSY_CONVECT_EMANUEL,        ONLY: ntrans

    IMPLICIT NONE

    ! LOCAL
    INTEGER     :: iou       ! I/O unit
    INTEGER     :: status    ! error status
    CHARACTER(LEN=*), PARAMETER :: substr = 'convection_init'
    ! INITIALIZE MAIN-CTRL
    CALL start_message_bi(modstr,'INITIALIZATION', substr)
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL convection init CORE ROUTINE:
       CALL init_convection(iou, status)
       ! op_mm_20140110 finish -> error_bi
       IF (status /= 0) CALL error_bi('Error calling init CORE ROUTINE',substr)
    END IF
 
    CALL p_bcast(lconvection, p_io)
    if (.not.lconvection) then
      if (p_parallel_io) write(*,*) "WARNING! No Convection active!"
    endif
    if (.not.lconvection) RETURN
    if (p_parallel_io) write(*,*) "Convection active!, Scheme No.: ", convect_param
    CALL p_bcast(convect_param, p_io)
    if (convect_param.le.3) then
      if (p_parallel_io) write(*,*) "Tiedtke - Scheme Nr. ", convect_param
      CALL p_bcast(lmfpen, p_io)  
      CALL p_bcast(lmfscv, p_io) 
      CALL p_bcast(lmfmid, p_io)  
      CALL p_bcast(lmfdd, p_io)   
      CALL p_bcast(lmfdudv, p_io)
      CALL p_bcast(lpos_def, p_io)
! fb_mk_20120116+
      CALL p_bcast(rset_cmfctop%l, p_io)
      CALL p_bcast(rset_cmfctop%v, p_io)
      CALL p_bcast(rset_cprcon%l, p_io)
      CALL p_bcast(rset_cprcon%v, p_io)
      CALL p_bcast(rset_entrscv%l, p_io)
      CALL p_bcast(rset_entrscv%v, p_io)
! fb_mk_20120116-
! fb_mk_20140210+
      CALL p_bcast(rset_entrpen%l, p_io)
      CALL p_bcast(rset_entrpen%v, p_io)
      CALL p_bcast(rset_entrmid%l, p_io)
      CALL p_bcast(rset_entrmid%v, p_io)
! fb_mk_20140210-
    ENDIF
    if (convect_param.eq.4) then
      if (p_parallel_io) write(*,*) "ECMWF - Scheme "
      CALL p_bcast(PEN, p_io)  
      CALL p_bcast(SHAL, p_io) 
      CALL p_bcast(MID, p_io)  
      CALL p_bcast(DOWN, p_io)   
      CALL p_bcast(FRIC, p_io) 
      CALL p_bcast(lmftrac, p_io) 
      CALL p_bcast(lepcld, p_io) 
      CALL p_bcast(lmfscl_wstar, p_io) 
    ENDIF
    if (convect_param.eq.5) then
      if (p_parallel_io) write(*,*) "Zhang - Hack - McFarlane - Scheme "
      call p_bcast(altconv, p_io)
      call p_bcast(evap_sub, p_io)
    ENDIF
    if (convect_param.eq.6) then
      if (p_parallel_io) write(*,*) "Bechtold - Scheme "
      call p_bcast(ODEEP, p_io)
      call p_bcast(OSHAL, p_io)
      call p_bcast(ODOWN, p_io)
      call p_bcast(OREFRESH_ALL, p_io)
      call p_bcast(OSETTADJ, p_io)
      call p_bcast(OUVTRANS, p_io)
      call p_bcast(OCHTRANS, p_io)
      call p_bcast(KENSM, p_io)
      call p_bcast(KICE, p_io)
      call p_bcast(PTADJD, p_io)
      call p_bcast(PTADJS, p_io)
    ENDIF
    IF (convect_param.eq.7) THEN
      IF (p_parallel_io) WRITE(*,*) "Emanuel - Scheme "
      CALL p_bcast(ltransport, p_io)
    ENDIF
    IF (convect_param.eq.8) THEN
      IF (p_parallel_io) WRITE(*,*) "Donner - Scheme "
    ENDIF

    ! op_pj_20091009+
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL convect_read_nml_cpl(status, iou)
       ! op_mm_20140110 finish -> error_bi
       IF (status /= 0) CALL error_bi('Error reading cpl nml', substr)
    END IF
    CALL p_bcast(l_lgmc_diag, p_io)
    ! op_pj_20091009-
    ! um_ak_20140617+
    CALL p_bcast(pbl_height%CHA, p_io)
    CALL p_bcast(pbl_height%OBJ, p_io)
    ! um_ak_20140617-


    ! ub_ak_20181026+
    ! moved to init_memory
!!$#ifdef MESSYTENDENCY
!!$    my_handle = mtend_get_handle(modstr)
!!$    CALL mtend_register (my_handle,mtend_id_t)
!!$    CALL mtend_register (my_handle,mtend_id_q)
!!$    CALL mtend_register (my_handle,mtend_id_tracer)
!!$    if (convect_param .ge. 4 .and. convect_param .le. 6)then! .or. convect_param .eq. 8) then
!!$       CALL mtend_register (my_handle,mtend_id_xl)
!!$       CALL mtend_register (my_handle,mtend_id_xi)
!!$    endif
!!$    if (convect_param .ge. 1 .and. convect_param .le. 4 .or. convect_param .ge. 6 .and. convect_param .le. 7) then
!!$       CALL mtend_register (my_handle,mtend_id_u)
!!$       CALL mtend_register (my_handle,mtend_id_v)
!!$    endif
!!$#endif
    ! ub_ak_20181026-

    CALL end_message_bi(modstr,'INITIALIZATION', substr)

  END SUBROUTINE convect_initialize
!==============================================================================

!==============================================================================

  SUBROUTINE convect_init_memory

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nlev
    USE messy_main_grid_def_bi,     ONLY: ceta
    USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi
    USE messy_main_channel_error_bi,ONLY: channel_halt
    USE messy_main_channel_bi,      ONLY: GP_3D_MID &
                                        , GP_2D_HORIZONTAL &
                                        , SCALAR ! op_pj_20100122
    ! MESSy
    USE messy_convect_tiedtke,      ONLY: cevapcu
    USE messy_convect_tiedtke_param, ONLY: cprcon ! op_pj_20100122
    USE messy_main_constants_mem,   ONLY: g
    USE messy_main_channel,         ONLY: new_channel, new_channel_object &
                                        , new_attribute
    USE messy_main_channel,         ONLY: new_channel_object_reference ! op_mm_20140327
    USE MESSY_CONVECT_EMANUEL,      ONLY: ntrans
    USE MESSY_CONVECT,              ONLY: ltransport
    USE messy_main_tracer_mem_bi,   ONLY: ntrac => ntrac_gp

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER::substr='convect_init_memory'
    INTEGER :: status

    ! ub_ak_20181026+
    ! moved here fron initialze
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register (my_handle,mtend_id_t)
    CALL mtend_register (my_handle,mtend_id_q)
    CALL mtend_register (my_handle,mtend_id_tracer)
    if (convect_param .ge. 4 .and. convect_param .le. 6)then! .or. convect_param .eq. 8) then
       CALL mtend_register (my_handle,mtend_id_xl)
       CALL mtend_register (my_handle,mtend_id_xi)
    endif
    if (convect_param .ge. 1 .and. convect_param .le. 4 .or. convect_param .ge. 6 .and. convect_param .le. 7) then
       CALL mtend_register (my_handle,mtend_id_u)
       CALL mtend_register (my_handle,mtend_id_v)
    endif
#endif
    ! ub_ak_20181026-

    CALL start_message_bi(modstr,'MEMORY INITIALIZATION', substr)

    ALLOCATE(cevapcu(nlev))
    cevapcu(1:nlev) = 1.93E-6*261.*SQRT(1.E3/(38.3*0.293)*SQRT(ceta(1:nlev)))*0.5/g

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_type' &
         , p2=conv_type, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_type' &
         , 'long_name', c='type of convection')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'counter_deep' &
         , p2=counter_deep, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'counter_deep' &
         , 'long_name', c='number of events for deep convection')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'counter_shal' &
         , p2=counter_shal, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'counter_shal' &
         , 'long_name', c='number of events for shallow convection')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'counter_midl' &
         , p2=counter_midl, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'counter_midl' &
         , 'long_name', c='number of events for midlevel convection')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cth_deep' &
         , p2=cth_deep, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cth_deep' &
         , 'long_name', c='cloud top height for deep convection')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cth_shal' &
         , p2=cth_shal, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cth_shal' &
         , 'long_name', c='cloud top height for shallow convection')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cth_midl' &
         , p2=cth_midl, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE. )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cth_midl' &
         , 'long_name', c='cloud top height for midlevel convection')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'massfu_deep', p3=massfu_deep )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu_deep', &
         'long_name', c='updraft mass flux for deep convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu_deep', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'massfu_shal', p3=massfu_shal )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu_shal', &
         'long_name', c='updraft mass flux for shallow convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu_shal', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'massfu_midl', p3=massfu_midl )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu_midl', &
         'long_name', c='updraft mass flux for midlevel convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu_midl', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_counter' &
         , p2=counter, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_counter', &
         'long_name', c='convective counter for Bechtold scheme')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_counter', 'units', c=' ')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_top', &
         p2=conv_top, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_top', &
         'long_name', c='top level of convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_top', 'units', c='levels')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_bot', &
         p2=conv_bot, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_bot', &
         'long_name', c='bottom level of convection / convective cloud base')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_bot', 'units', c='levels')
    CALL channel_halt(substr, status)
!!$    ! op_pj_20170324+
!!$    ! create a reference with different name, which is consistent with the
!!$    ! requirements of VISO ...
!!$    CALL new_channel_object_reference(status, modstr &
!!$         , 'conv_bot', modstr, 'conv_bot_i')
!!$    CALL channel_halt(substr, status)
!!$    ! op_pj_20170324-

    CALL new_channel_object(status, modstr, 'cu_bot', &
         p2=cu_bot, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_bot', &
         'long_name', c='bottom level of deep(!) convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_bot', 'units', c='levels')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cu_top', &
         p2=cu_top, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_top', &
         'long_name', c='top level of deep(!) convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_top', 'units', c='levels')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cu_freeze', &
         p2=cu_freeze, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_freeze', &
         'long_name', c='freezing level in deep(!) convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_freeze', 'units', c='levels')
    CALL channel_halt(substr, status)

    ! mz_pj_20050615+
    CALL new_channel_object(status, modstr, 'cu_bot_mid', &
         p2=cu_bot_mid, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_bot_mid', &
         'long_name', c='bottom level of midlevel convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_bot_mid', 'units', c='levels')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cu_top_mid', &
         p2=cu_top_mid, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_top_mid', &
         'long_name', c='top level of midlevel convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_top_mid', 'units', c='levels')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cu_freeze_mid', &
         p2=cu_freeze_mid, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_freeze_mid', &
         'long_name', c='freezing level in midlevel convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_freeze_mid', 'units', c='levels')
    CALL channel_halt(substr, status)
    ! op_pj_20100122+
    IF (convect_param <= 3) THEN
       CALL new_channel_object(status, modstr, 'cprcon', &
            p0=cprcon, reprid=SCALAR )
       CALL channel_halt(substr, status)
    END IF
    ! op_pj_20100122-
    ! mz_pj_20050615-

    CALL new_channel_object(status, modstr, 'cu_uvelo', p3=cu_uvelo )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_uvelo', &
         'long_name', c='updraft velocity in deep(!) convection')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cu_uvelo', 'units', c='m/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'base_f1' &
         , p2=base_f1, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'base_f1', &
         'long_name', c='cloud base updraft mass flux 1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'base_f1', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'base_f2' &
         , p2=base_f2, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'base_f2', &
         'long_name', c='cloud base updraft mass flux 2')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'base_f2', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'base_f3' &
         , p2=base_f3, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'base_f3', &
         'long_name', c='cloud base updraft mass flux 3')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'base_f3', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'base_f4' &
         , p2=base_f4, reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'base_f4', &
         'long_name', c='cloud base updraft mass flux 4')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'base_f4', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'massfu', p3=massfu )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu', &
         'long_name', c='updraft mass flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'massfd', p3=massfd )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfd', &
         'long_name', c='downward mass flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfd', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'u_entr', p3=u_entr )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u_entr', &
         'long_name', c='upward entraining mass flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u_entr', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'd_entr', p3=d_entr )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'd_entr', &
         'long_name', c='downward entraining mass flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'd_entr', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'u_detr', p3=u_detr )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u_detr', &
         'long_name', c='upward detraining mass flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u_detr', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'udetr_h', p3=udetr_h )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'udetr_h', &
         'long_name', c='upward detraining mass flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'udetr_h', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'd_detr', p3=d_detr )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'd_detr', &
         'long_name', c='downward detraining mass flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'd_detr', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'massfu_asc', p3=massfu_asc )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu_asc', &
         'long_name', c='updraft mass flux after ascend')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfu_asc', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'massfd_draf', p3=massfd_draf )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfd_draf', &
         'long_name', c='downward mass flux after descend')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'massfd_draf', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_precflx', p3=cv_precflx )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_precflx', &
         'long_name', c='convective precipitation flux (3d)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_precflx', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_snowflx', p3=cv_snowflx )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_snowflx', &
         'long_name', c='convective snow precipitation flux (3d)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_snowflx', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_precnew', p3=cv_precnew )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_precnew', &
         'long_name', c='freshly formed precipitation flux (3d)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_precnew', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_snownew', p3=cv_snownew )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_snownew', &
         'long_name', c='freshly formed snow flux (3d)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_snownew', 'units', c='kg /(m^2 s)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_cover', p3=cv_cover )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_cover', &
         'long_name', c='est. convective cloud cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_cover', 'units', c='0 - 1')
    CALL channel_halt(substr, status)
    ! cv_cover is estimated from updraft strength

    ! mz_jd_20161011+
    CALL new_channel_object(status, modstr, 'qstdev', p3=qstdev)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qstdev', &
         'long_name', c='standard deviation of total water amount')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qstdev', 'units', c='kg/kg')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'qskew', p3=qskew)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qskew', &
         'long_name', c='skewness of total water amount')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qskew', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'qsat', p3=qsat)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qsat', &
         'long_name', c='saturation specific humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qsat', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'qp1', p3=qp1)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qp1', &
         'long_name', c='total specific humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qp1', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status,modstr, 'Q2', p3=Q2)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Q2', &
         'long_name', c='proxy for area fraction of cloud cores')
    CALL channel_halt(substr,status)
    CALL new_attribute(status, modstr, 'Q2', 'units', c='kg/kg')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_cover_sikma', p3=cv_cover_sikma )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_cover_sikma', &
         'long_name', c='convective cloud cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_cover_sikma', 'units', c='-')
    CALL channel_halt(substr, status)
    ! cv_cover_sikma is calculated using the spatial standard deviation 
    ! of specific humidity.
    ! mz_jd_20161011-

    CALL new_channel_object(status, modstr, 'conv_tte', p3=conv_tte )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_tte', &
         'long_name', c='convective temperature tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_tte', 'units', c='K/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_qte', p3=conv_qte )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_qte', &
         'long_name', c='convective humidity tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_qte', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_lte', p3=conv_lte )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_lte', &
         'long_name', c='convective liquid water tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_lte', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_ite', p3=conv_ite )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_ite', &
         'long_name', c='convective ice tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_ite', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_ute', p3=conv_ute )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_ute', &
         'long_name', c='convective u tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_ute', 'units', c='m/s^2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'conv_vte', p3=conv_vte )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_vte', &
         'long_name', c='convective v tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'conv_vte', 'units', c='m/s^2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_lwc', p3=cv_lwc)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_lwc', &
         'long_name', c='convective cloud water content (3d in cloud)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_lwc', 'units', c='kg/kg')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_iwc', p3=cv_iwc)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_iwc', &
         'long_name', c='convective cloud ice content(3d in cloud)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_iwc', 'units', c='kg/kg')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_rform', p3=cv_rform)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_rform', &
         'long_name', c='convective precipitation formation (water, in cloud value)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_rform', 'units', c='kg/kg')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cv_sform', p3=cv_sform)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_sform', &
         'long_name', c='convective precipitation formation (snow, in cloud value)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_sform', 'units', c='kg/kg')
    CALL channel_halt(substr, status)

    IF (convect_param.eq.6) THEN
      CALL new_channel_object(status, modstr, 'WAT_DIAG', p2=WAT_DIAG,&
        reprid=GP_2D_HORIZONTAL )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'WAT_DIAG', &
          'long_name', c='relative correction of water to total precip')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'WAT_DIAG', 'units', c=' - ')
      CALL channel_halt(substr, status)
    ENDIF

    ! mz_ak_20051221+
    CALL new_channel_object(status, modstr, 'cv_cldwater', p3=cv_cldwater )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_cldwater', &
         'long_name', c='convective cloud water content(3d)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cv_cldwater', 'units', c='kg/kg')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'del_liqwater', p3=del_liqwat )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'del_liqwater', &
         'long_name', c='change in liquid water content of conv cloud')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'del_liqwater', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! mz_ak_20051221-

    CALL new_channel_object(status, modstr, 'cth', p2=cth, &
      reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cth', &
         'long_name', c='convective cloud top height')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cth', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'CAPE', p2=cape, &
      reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CAPE', &
         'long_name', c='CAPE (convective available potential energy)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CAPE', 'units', c='m^2/s^2')
    CALL channel_halt(substr, status)

#ifdef CESM1
    ! mz_ab_20130825+
    CALL new_channel_object(status, modstr, 'ilab' &
         , p3=pilab, reprid=GP_3D_MID )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ilab' &
         , 'long_name', c='index: below or in convective cloud')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ilab', 'units', c='-')
    CALL channel_halt(substr, status)
    ! mz_ab_20130825-
#endif

    IF (convect_param == 7) THEN
      IF (ltransport) THEN
        ntrans = ntrac
      ELSE
        ntrans = 0
      ENDIF

      CALL new_channel_object(status, modstr, 'CBMF_E', p2=CBMF,&
        reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE. )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CBMF_E', &
          'long_name', c='cloud base mass flux (Emanuel scheme)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CBMF_E', 'units', c='kg/(m^2*s)')
      CALL channel_halt(substr, status)
    ENDIF

    IF (convect_param == 8) THEN
      CALL new_channel_object(status, modstr, 'c_cldfrac', p3=c_cldfrac )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_cldfrac', &
        'long_name', c='fractional coverage of convective cells in grid box')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_cldfrac', 'units', c='-')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'c_liqamt', p3=c_liqamt )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_liqamt', &
        'long_name', c='liquid water content of convective cells')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_liqamt', 'units', c='kg/kg')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'c_liqsize', p3=c_liqsize )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_liqsize', &
        'long_name', c='assumed effective size of cell liquid drops')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_liqsize', 'units', c='1e-6 m')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'c_iceamt', p3=c_iceamt )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_iceamt', &
        'long_name', c='ice water content of cells')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_iceamt', 'units', c='kg/kg')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'c_icesize', p3=c_icesize )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_icesize', &
        'long_name',                                  &
        c='generalized effective diameter for ice in convective cells')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'c_icesize', 'units', c='1e-6 m')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'm_cldfrac', p3=m_cldfrac )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_cldfrac', &
        'long_name', c='fractional area of mesoscale clouds in grid box')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_cldfrac', 'units', c='-')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'm_liqamt', p3=m_liqamt )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_liqamt', &
        'long_name', c='liquid water content in mesoscale clouds')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_liqamt', 'units', c='kg/kg')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'm_liqsize', p3=m_liqsize )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_liqsize', &
        'long_name', c='assumed effective size of mesoscale drops')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_liqsize', 'units', c='1e-6 m')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'm_iceamt', p3=m_iceamt )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_iceamt', &
        'long_name', c='ice water content of mesoscale elements')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_iceamt', 'units', c='kg/kg')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'm_icesize', p3=m_icesize)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_icesize', &
        'long_name',                                  &
        c='generalized ice effective size for anvil ice')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'm_icesize', 'units', c='1e-6 m')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'd_humarea', p3=d_humarea)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'd_humarea', &
        'long_name',                                  &
        c='affected fraction for humidity')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'd_humarea', 'units', c='-')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'd_humratio', p3=d_humratio)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'd_humratio', &
        'long_name',                                  &
        c='humidity ratio in and out if convective cloud')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'd_humratio', 'units', c='-')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'conv_covte', p3=conv_covte )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'conv_covte', &
        'long_name', c='convective cover tendency')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'conv_covte', 'units', c='1/s')
      CALL channel_halt(substr, status)
    ENDIF

    ! mim_sb_20091009+
    IF (l_lgmc_diag) THEN
       CALL new_channel(status, modstr//'_lg', reprid=GP_3D_MID)
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//'_lg', 'conv_tte_up' &
            , p3=conv_tte_up )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_up', &
            'long_name', c='temperature tendency updraft move')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_up', 'units', c='K/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_tte_do' &
            , p3=conv_tte_do )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_do', &
            'long_name', c='temperature tendency downdraft move')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_do', 'units', c='K/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_tte_up_cond' &
            , p3=conv_tte_up_cond )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_up_cond', &
            'long_name', c='temperature tendency updraft condensation')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_up_cond' &
            , 'units', c='K/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_tte_up_freeze' &
            , p3=conv_tte_up_freeze )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_up_freeze', &
            'long_name', c='temperature tendency updraft freezing')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_up_freeze' &
            , 'units', c='K/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_tte_do_verd' &
            , p3=conv_tte_do_verd )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_do_verd', &
            'long_name', c='temperature tendency evaporation of rain in downdraft')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_do_verd', 'units' &
            , c='K/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_tte_do_melt' &
            , p3=conv_tte_do_melt )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_do_melt', &
            'long_name', c='temperature tendency downdraft melting of snow')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_do_melt', 'units' &
            , c='K/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_tte_su' &
            , p3=conv_tte_su )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_su', &
            'long_name', c='temperature tendency subsidence')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_su', 'units', c='K/s')
       CALL channel_halt(substr, status)
       
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_tte_ev' &
            , p3=conv_tte_ev )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_ev', &
            'long_name' &
            , c='temperature tendency transport')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_tte_ev', 'units', c='K/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_qte_up' &
            , p3=conv_qte_up )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_qte_up', &
            'long_name', c='moisture tendency updraft')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_qte_up', 'units' &
            , c='kg/kg/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_qte_do' &
            , p3=conv_qte_do )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_qte_do', &
            'long_name', c='moisture tendency downdraft')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_qte_do' &
            , 'units', c='kg/kg/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_qte_su' &
            , p3=conv_qte_su )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_qte_su', &
            'long_name', c='moisture tendency subsidence')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_qte_su' &
            , 'units', c='kg/kg/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_qte_ev' &
            , p3=conv_qte_ev )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_qte_ev', &
            'long_name' &
            , c='moisture tendency evaporation of rain below cloud base')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_qte_ev' &
            , 'units', c='kg/kg/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_pqtec' &
            , p3=conv_pqtec )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_pqtec', &
            'long_name', c='moisture tendency detrainment')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_pqtec' &
            , 'units', c='kg/kg/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'conv_pxtec' &
            , p3=conv_pxtec )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_pxtec', &
            'long_name', c='liquid water tendency detrainment')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'conv_pxtec' &
            , 'units', c='kg/kg/s')
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr//'_lg', 'aprsc', p2=aprsc, &
             reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'aprsc', &
            'long_name', c='convective precpitation')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'aprsc', 'units', c='kg/m2')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//'_lg', 'aprss', p2=aprss, &
             reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'aprss', &
            'long_name', c='snow precpitation')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'aprss', 'units', c='kg/m2')
       CALL channel_halt(substr, status)
    END IF
    ! mim_sb_20091009-

    ! op_mm_20140226+
    ! for calculation of max convective gust (as in COSMO)
    CALL new_channel_object(status, modstr, 'vgustcon', p2=vgustcon,&
         reprid=GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vgustcon', &
         'long_name', c='maximum convective gust at 10m  ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vgustcon', 'units', c=' m/s ')
    CALL channel_halt(substr, status)
    ! op_mm_20140226-

    ! op_mm_20140226+
    !  convective buoyant TKE production (as in COSMO) 
    CALL new_channel_object(status, modstr, 'tketconv', p3=tketconv,&
         reprid=GP_3D_MID )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tketconv', &
         'long_name', c=' TKE-tendency due to convective buoyancy     ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tketconv', 'units', c=' m2/s3  ')
    CALL channel_halt(substr, status)
    ! op_mm_20140226-



#ifdef COSMO
    ! op_mm_20140221+
         ! updraft mass flux
       CALL new_channel_object_reference(status,modstr, 'massfu', 'COSMO', 'massfu')
       ! downward massflux
       CALL new_channel_object_reference(status,modstr, 'massfd','COSMO' , 'massfd')
       ! upward entraining mass flux
      CALL new_channel_object_reference(status,modstr, 'u_entr','COSMO' , 'u_entr')
       ! downward entraining mass flux
      CALL new_channel_object_reference(status,modstr, 'd_entr','COSMO' , 'd_entr')
       ! upward detraining mass flux
      CALL new_channel_object_reference(status,modstr, 'u_detr','COSMO' , 'u_detr')
       ! downward detraining mass flux
      CALL new_channel_object_reference(status,modstr, 'd_detr','COSMO', 'u_detr')
      ! convective precipitation flux (3d)
      CALL new_channel_object_reference(status,modstr, 'cv_precflx','COSMO' , 'cv_precflx')
      ! convective snow precipitation flux (3d)
      CALL new_channel_object_reference(status,modstr, 'cv_snowflx','COSMO', 'cv_snowflx')
      ! freshly formed convective precipitation flux (3d)
      CALL new_channel_object_reference(status,modstr, 'cv_precnew','COSMO' , 'cv_precnew')
      ! freshly formed convective snow precipitation flux (3d)
      CALL new_channel_object_reference(status,modstr, 'cv_snownew','COSMO' , 'cv_snownew')
      ! convective cloud water content(3d)
      CALL new_channel_object_reference(status,modstr, 'cv_lwc','COSMO' , 'cv_lwc')
      !convective cloud ice content(3d)'
      CALL new_channel_object_reference(status,modstr, 'cv_iwc','COSMO', 'cv_iwc')
      !convective precipitation formation (water)
      CALL new_channel_object_reference(status,modstr, 'cv_rform','COSMO' , 'cv_rform')
      !convective precipitation formation (snow)
      CALL new_channel_object_reference(status,modstr, 'cv_sform','COSMO', 'cv_sform')
      !level index of convective cloud top
      CALL new_channel_object_reference(status,modstr, 'cu_top','COSMO' , 'cu_top')
      !level index of convective cloud bottom
      CALL new_channel_object_reference(status,modstr, 'cu_bot','COSMO' , 'cu_bot')
      !level index of convective cloud top
      CALL new_channel_object_reference(status,modstr, 'conv_top','COSMO' , 'conv_top')
      !level index of convective cloud bottom
      CALL new_channel_object_reference(status,modstr, 'conv_bot','COSMO' , 'conv_bot')
    ! op_mm_20140221-

#endif


    CALL end_message_bi(modstr,'MEMORY INITIALIZATION', substr)

  END SUBROUTINE convect_init_memory

!==============================================================================
 ! coupling with physc for some channel objects

  SUBROUTINE convect_init_coupling
  
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY:  p_parallel_io
#ifdef ECHAM5
    USE messy_main_data_bi,          ONLY: ltdiag
#endif
    USE messy_main_grid_def_mem_bi,  ONLY: nproma
    USE messy_main_data_bi,          ONLY: l_heatflux, s_heatflux

    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi &
                                         , error_bi, info_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
#if defined(CESM1)
    USE messy_convect_tiedtke_param, ONLY: CUPARAM_INIT
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, nlevp1, nvclev, nn, lmidatm, &
                                           vct
    USE messy_main_data_bi,          ONLY: modstr_base=>modstr
#endif
#if defined(ECHAM5)
    USE messy_convect_tiedtke_param,  ONLY: CUPARAM_INIT_1
    USE messy_main_grid_def_mem_bi,   ONLY: nlev, nlevp1, nvclev, nn, lmidatm, &
                                            vct
    USE messy_main_data_bi,           ONLY: modstr_base=>modstr, &
                                           lcouple
#endif   
 ! op_mm_20140122+
#ifdef COSMO
    USE messy_main_grid_def_mem_bi,   ONLY: nlev, nlevp1, nn, lmidatm
    USE messy_main_grid_def_bi,       ONLY: h_a=>hyai,h_b=>hybi
    USE messy_main_data_bi,           ONLY: modstr_base=>modstr
    !using same parameters as in cosmo-convection-scheme
    USE messy_convect_tiedtke_param,  ONLY: CUPARAM_INIT_2  
   
    ! op_mm_20140122-
#endif
#ifdef MESSYDWARF
    USE messy_main_grid_def_mem_bi,   ONLY: nlev, nlevp1, nvclev, nn, lmidatm, &
                                            h_a=>hyai,h_b=>hybi
    USE messy_main_data_bi,           ONLY: modstr_base=> modstr
#endif
    USE messy_convect_ecmwf_param,    ONLY: init_convection_constants
    USE messy_convect_zhang_param,    ONLY: esinti, conv_ini, mfinti
    USE messy_convect_donner_additions, ONLY: donner_init
    USE messy_convect,                ONLY: modstr
    USE messy_main_channel,           ONLY: get_channel_object, get_channel_info, &
                                            new_channel_object, new_attribute ! mim_sb_20091207

    IMPLICIT NONE

    REAL(dp) :: hypi(nlevp1), sfpress
#if !defined(COSMO) && !defined(MESSYDWARF)
    REAL(dp) :: h_a(nvclev), h_b(nvclev)
#endif
    LOGICAL  :: xip, xtrigon, output_p

    INTEGER :: status, jk
    CHARACTER(LEN=*), PARAMETER::substr='convect_init_cpl'
  
    output_p = p_parallel_io
    status = 0
    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)
    ! initialization routine for parameters for tiedtke convection scheme
   
    if (convect_param.le.3) THEN
#if defined(CESM1)
       CALL CUPARAM_INIT(status, nlev, nlevp1, nvclev, REAL(vct,dp), &
                                              nn, lmidatm, output_p)
#endif
#if defined(ECHAM5)
       CALL CUPARAM_INIT_1(status, nlev, nlevp1, nvclev, REAL(vct,dp), &
                                                nn, lmidatm, output_p, &
                                                lcouple, .FALSE.)
#endif

! op_mm_20140122+
! set sama parameter as in COSMO
#if defined(COSMO) || defined(MESSYDWARF)
       CALL CUPARAM_INIT_2(status, nlev,h_a,h_b,output_p)
#endif      
! op_mm_20140122-


    end if

! op_mm_20140110 finish -> error_bi
    if (status.ne.0) call error_bi( &
       'problem with invalid combination of vertical and horizontal resolution'&
       , 'messy_convect_tiedtke_param.f90')

    if (convect_param.eq.4) CALL init_convection_constants(nn, nlev)

    if (convect_param.eq.5) then
      xtrigon=.true.
      xip =  .false.

#if !defined(COSMO) && !defined(MESSYDWARF)
      do jk=1,nvclev
        h_a(jk) = vct(jk)
        h_b(jk) = vct(jk+nvclev)
      enddo
#endif
    
      sfpress     = 1.e5_dp   ! reference pressure of 1000 hPa
      
      do jk=1,nlev+1
        hypi(jk)      = h_a(jk) + h_b(jk) * sfpress
      enddo
       
      call conv_ini(nlev, nlevp1, hypi, xtrigon, output_p)
      call mfinti(nlev, nlevp1, hypi, output_p)  
      call esinti(xip, status, output_p )
! op_mm_20140110 finish -> error_bi
      if (status.ne.0) call error_bi(&
           'problem with initalizing lookup_tables' &
           , 'messy_convect_zhang_param.f90')
    endif

    if (convect_param == 8) CALL DONNER_INIT(nproma, 1, nlev)

    CALL get_channel_info(status, 'cvtrans')
    cvtrans = (status == 0)
   
#ifndef CESM1
    CALL get_channel_object(status, modstr_base, 'ilab', p3=pilab )
    IF (status /= 0) &
    ! op_mm_20140110+ finish -> error_bi
       call error_bi('channel object for ilab not found', substr)
#endif

    CALL get_channel_object(status, modstr_base, 'qtec', p3=pqtec )
    IF (status /= 0) &
       call error_bi('channel object for qtec not found', substr)

    CALL get_channel_object(status, modstr_base, 'qflux', p2=pqflux )
    IF (status /= 0) &
       call error_bi('channel object for qflux not found', substr)

    CALL get_channel_object(status, 'grid_def', 'grmass', p3=grmass )
    IF (status /= 0) &
       call error_bi('channel object for grmass not found', substr)

    CALL get_channel_object(status, 'grid_def', 'grvol', p3=grvol )
    IF (status /= 0) &
       call error_bi('channel object for grvol not found', substr)

    CALL get_channel_object(status, modstr_base, 'geopoti', p3=geopoti )
    IF (status /= 0) &
      call error_bi('channel object for geopoti not found', substr)

    CALL get_channel_object(status, modstr_base, 'geopot', p3=geopot )
    IF (status /= 0) &
      call error_bi('channel object for geopot not found', substr)
    ! op_mm_20140110-

    ! op_pj_20171019+
    ! NOTE: This should ultimately be replaced by local pointers set
    !       with get_channel_object, however, this needs to be done
    !       for ALL basemodels (and not all have channel objects already).
    IF (.NOT. ASSOCIATED(l_heatflux)) &
         call error_bi('l_heatflux not associated', substr)
    IF (.NOT. ASSOCIATED(s_heatflux)) &
         call error_bi('s_heatflux not associated', substr)
    ! op_pj_20171019-

    if (convect_param.eq.5) then
      if (.not.cvtrans.and.output_p) then
        write(*,*) "####################################################"
        write(*,*) "     WARNING: No Convective Transport chosen !!!" 
        write(*,*) "####################################################"
      ENDIF
      write(*,*) 'looking for channel /object: tropop / pblh'
! um_ak_20140502+
!!$#ifdef COSMO
!!$! op_mm_20140327 different name in COSMO
!!$      call get_channel_object(status, 'tropop', 'pblhRi', p2=pblh)
!!$      CALL channel_halt(substr, status)
!!$#else
!!$      call get_channel_object(status, 'tropop', 'pblh', p2=pblh)
!!$      CALL channel_halt(substr, status)
!!$#endif
      IF (TRIM(pbl_height%CHA) /= '' .AND. TRIM(pbl_height%OBJ) /= '') THEN
         CALL get_channel_object(status &
              ,TRIM(pbl_height%CHA), TRIM(pbl_height%OBJ), p2=pblh)
         CALL channel_halt(substr, status)         
      ELSE
         ! keep status quo for EMAC model
         CALL get_channel_object(status, 'tropop', 'pblh', p2=pblh)
         CALL channel_halt(substr, status)
      END IF
! um_ak_20140502-
    ENDIF

    ! op_mm_20140110 message -> info_bi
#ifdef ECHAM5
    IF (ltdiag) THEN
       CALL info_bi(substr,'looking for channel / object tdiag / PDIGA5' )
       CALL get_channel_object(status, 'tdiag', 'PDIGA5', p3=pdiga5)
       CALL channel_halt(substr, status)

       CALL info_bi(substr,'looking for channel / object tdiag / PDIGA10')
       CALL get_channel_object(status, 'tdiag', 'PDIGA10', p3=pdiga10)
       CALL channel_halt(substr, status)

       CALL info_bi(substr,'looking for channel / object tdiag / PDIGA18')
       CALL get_channel_object(status, 'tdiag', 'PDIGA18', p3=pdiga18)
       CALL channel_halt(substr, status)
    END IF
#endif

    ! mim_sb_20091009+
    CALL get_channel_info(status, 'lgmc')
    L_LGMC = (status == 0)

    IF (L_LGMC) THEN
       CALL new_channel_object(status, modstr//'_lg', 'ttp1_gp', p3=ttp1_gp )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'ttp1_gp', &
            'long_name', c='temperature before updraft')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'ttp1_gp', 'units', c='K')
       CALL channel_halt(substr, status)

       ! op_sb_20130417+
       CALL new_channel_object(status, modstr//'_lg', 'ptu_gp', p3=ptu_gp )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'ptu_gp', &
            'long_name', c='temperature within updraft')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'ptu_gp', 'units', c='K')

       CALL channel_halt(substr, status)
       CALL new_channel_object(status, modstr//'_lg', 'ptd_gp', p3=ptd_gp )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'ptd_gp', &
            'long_name', c='temperature within downdraft')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'ptd_gp', 'units', c='K')
       CALL channel_halt(substr, status)
       ! op_sb_20130417-
    END IF
    ! mim_sb_20091009-

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)    

  END SUBROUTINE convect_init_coupling


!==============================================================================

  SUBROUTINE convect_free_memory

    ! MESSy
    USE messy_convect_tiedtke,     ONLY: cevapcu

    IMPLICIT NONE

    IF (ALLOCATED(cevapcu)) DEALLOCATE(cevapcu)

  END SUBROUTINE convect_free_memory


!==============================================================================

!==============================================================================
  SUBROUTINE convect_CONVEC
   
!
!
!          *CONVCALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
!                     *CUMASTR* (CUMULUS PARAMETERIZATION)
!
!           M.TIEDTKE      E.C.M.W.F.     12/1989
!
!**   PURPOSE.
!     --------
!
!          *CONVCALL* - INTERFACE FOR *CUMASTR*:
!                     PROVIDES INPUT FOR CUMASTR
!                     RECEIVES UPDATED TENDENCIES, PRECIPITATION.
!
!**   INTERFACE.
!     ----------
!
!          *CONVCALL* IS CALLED FROM *PHYSC*
!
!     EXTERNALS.
!     ----------
!
!          CUMASTR, CUMASTRT OR CUMASTRH
!

    ! ECHAM5/MESSy

! op_mm_20140110
#ifdef ECHAM5
  USE messy_main_data_bi, ONLY: ltdiag
#endif

  USE messy_main_timer,   ONLY: time_step_len,        delta_time,    &
                                nstep => current_time_step
  USE messy_main_grid_def_mem_bi, ONLY: nlev,      nlevp1,    nn,   &
                                        nlevm1,    nproma,    kproma, jrow
                               
  USE messy_main_grid_def_bi,     ONLY: gboxarea_2d
  USE messy_main_data_bi,ONLY: xlm1,      xim1,      tm1,    qm1,   &
                               vm1,       um1,                      &
                               topmax,    aprc,      aprs,   xtec,  &
                               qte_3d,    tte_3d,   vol_3d, vom_3d, &
                               vervel_3d, xlte_3d,  xite_3d,        &
#ifdef ECHAM5
                               aphp1,     app1,                     &
#else
                               press_3d,   pressi_3d,               &
#endif
                               geopot_3d, loland_2d,                &
                               rsfc_2d,   ssfc_2d,                  &
                               zsenkf_2d, zlatkf_2d, zust_2d,       &
                               tsurf_2d,                            &
                               L_HEATFLUX, aclc, slf,               & 
                               S_HEATFLUX,                          &
                               xvar, xskew, & ! mz_jd_20161011
                               ledith         ! op_pj_20180725

#ifdef CESM1
  USE messy_main_data_bi,ONLY: aprflux ! mz_ab_20140713
#endif
       
  USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1, &
                                       ntrac => ntrac_gp, ti_gp
  USE messy_main_tracer,         ONLY: I_convect, ON
  USE messy_main_constants_mem,  ONLY: M_air, M_H2O, g, R_gas, Tmelt &
                                     , vtmpc1 ! op_re_20130718
! op_mm_20140110 finish -> error_bi
  USE messy_main_mpi_bi,         ONLY:  p_pe !finish
  USE messy_main_blather_bi,     ONLY: error_bi

  USE messy_convect_zhang_param, ONLY: cp, rgas

  USE MESSY_convect,              ONLY: tiedtke_cumastr, tiedtke_cumastrh,    &
                                        tiedtke_cumastrt,                     &
                                        ecmwf_cumastr,                        &
                                        zhang_cumastr,                        &
                                        bechtold_cumastr,                     &
                                        emanuel_cumastr,                      &
                                        donner_cumastr,                       &
                                        calc_conv_cover, lconvection,         &
                                        calc_conv_cover_sikma ! mz_jd_20161011
  USE messy_convect_tiedtke,      ONLY: lookupoverflow
  USE messy_main_tools,           ONLY: tlucua, jptlucu1, jptlucu2


  ! op_mm_20140123+
  ! getting alle the COSMO fields 
#ifdef COSMO
   USE data_fields     , ONLY :&
       !   fields from the convection scheme
    clc_con     ,   & ! cloud cover due to convection                   --
!_cdm All "clw" variables now contain the mixed-phase convective cloud condensate.
    clw_con     ,   & ! cloud liquid water due to convection            --
    prr_con     ,   & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con     ,   & ! precipitation rate of snow, convective        (kg/m2*s)
    prne_con    ,   & ! precipitation rate, no evaporat., convective  (kg/m2*s)
    bas_con     ,   & ! level index of convective cloud base            -- 
    top_con     ,   & ! level index of convective cloud base            --
! um_ak_20140502+
! avoid double treatment of tendencies
!!$  tt_conv     ,   & ! temperature tendency due to convection        ( K/s  )
!!$  qvt_conv    ,   & ! humidity    tendency due to convection        ( 1/s  )
!!$  qct_conv    ,   & ! qc-tendency tendency due to convection        ( 1/s  )
!!$  qit_conv    ,   & ! qi-tendency tendency due to convection        ( 1/s  )
  ut_conv     ,   & ! u-tendency due to convection                  ( m/s^2)
  vt_conv     ,   & ! v-tendency due to convection                  ( m/s^2)
! um_ak_20140502-
    tket_conv   ,   & ! TKE-tendency due to convective buoyancy       ( m2/s3 )
#ifndef COSMOv5s5
    qvsflx      ,   & ! surface flux of water vapour                  (kg/m2*s)
#else
    qvfl_s      ,   &
#endif
    mflx_con    ,   & ! cloud base massflux                           (kg/m2*s)
    cape_con    ,   & ! convective available energy                   (   J/kg)
    qcvg_con    ,   & ! moisture convergence for Kuo-type closure     (    1/s)
    tke_con     ,   & ! convective turbulent energy                   (   J/kg)  !MR: not yet defined

   vgust_con         ! maximum convective gust at 10m                ( m/s )

#endif 
  ! op_mm_20140123-

  IMPLICIT NONE

! op_re_20130718+
! replaced by USE from messy_main_constants_mem for internal consistency
!!$  REAL(DP), PARAMETER :: vtmpc1 = M_air / M_H2O - 1.0_dp
! op_re_20130718-

  INTEGER :: klev, kbdim, ktrac, klevm1, klevp1

  REAL(dp), POINTER, DIMENSION(:,:)    :: pxtec   => NULL()
  REAL(dp), POINTER, DIMENSION(:)      :: paprc   => NULL()
  REAL(dp), POINTER, DIMENSION(:)      :: paprs   => NULL()
  REAL(dp), POINTER, DIMENSION(:)      :: ptopmax => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)    :: zmfu    => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)    :: zmfd    => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)    :: vervel  => NULL()

  INTEGER  :: ktype(nproma)

  REAL(dp) :: ztp1(nproma,nlev),         zqp1(nproma,nlev),              &
              zxp1(nproma,nlev),         ztvp1(nproma,nlev),             &
              zup1(nproma,nlev),         zvp1(nproma,nlev)

  REAL(dp) :: zqsat(nproma,nlev),        zrain(nproma)
  INTEGER  :: itopec2(nproma)
  INTEGER  :: icbot(nproma),             ictop(nproma)
  REAL(dp) :: zxtp1(nproma,nlev,ntrac),  zxtu(nproma,nlev,ntrac)
             
  REAL(dp) :: ztopmax(nproma)
  LOGICAL  :: locum(nproma)        ! logical convective activity

!  Local variables 
  REAL(dp) :: ztmst, zxlp1(nproma,nlev), zxip1(nproma,nlev), zrhoa(nproma,nlev)
  REAL(dp) :: zwp1(nproma, nlev, 2)
  INTEGER  :: ilevmin, jk, jl, jt, it
  INTEGER  :: ztracconv(ntrac)
! mz_ht_20050222+
! new variables for ECMWF Convection
  LOGICAL  :: loshal(nproma)        ! logical shallow convective activity
  INTEGER  :: KBOTSC(nproma)        ! index of shallow cloud base
  REAL(dp) :: psflx(nproma, nlevp1) ! snowflux through bottom of layer
  REAL(dp) :: psstru(nproma), psstrv(nproma)  
  REAL(dp) :: PTU(nproma,nLEV) 
  REAL(dp) :: PQU(nproma,nLEV) 
  REAL(dp) :: PLU(nproma,nLEV)
  REAL(dp) :: PRAIN(nproma)
  REAL(dp) :: THFLX_EC(nproma,nlevp1), QFLX_EC(nproma,nlevp1)
  REAL(dp) :: sum_wat1(nproma), sum_wat2(nproma)
! mz_ht_20050222-
! mz_ht_20040510+
! new variables for Zhang/Hack Convection
  REAL(dp) :: pdel(nproma,nlev), rpdel(nproma,nlev), & ! delta p and 1/deltap               
              zm(nproma,nlev), zi(nproma,nlevp1)       ! height of middle of layers, height of interfaces
  INTEGER  :: ideep(nproma)                                ! index of deep convection 
  REAL(dp) :: tpert(nproma), qpert(nproma)             ! perturbations of t and q due to boundary layer physics
  REAL(dp) :: wpert3d(nproma,nlev)                         ! perturbation of omega due to boundary layer physics
  REAL(dp) :: pflx(nproma, nlevp1)                         ! rainflux through bottom of layer

  REAL(dp) :: ts(nproma) , phis(nproma)                ! surface temperature and geopotential
  
  REAL(dp) :: obklen(nproma)                               ! obhukhov length
  REAL(dp) :: theta_surf(nproma), thvsrf(nproma)           ! potential / potential virtual temperature at surface
  REAL(dp) :: phiminv(nproma), phihinv(nproma)

! new variables for BECHTOLD
  REAL(dp) :: TRACTEN(nproma, nlev, ntrac)                 ! convective tracer tendency by convection
  REAL(dp) :: zw(nproma, nlev)                             ! vertical velocity
  REAL(dp) :: THFLX(nproma, nlev)                          ! turbulent sensible heat flux for Bechtold
  INTEGER  :: ccount(nproma)                               ! convective counter

! new variables for EMANUEL
  REAL(dp) :: ZWD(nproma), tprime(nproma), qprime(nproma), cbh
  INTEGER  :: itype(nproma)

  REAL(dp) :: zpb, zpt, grheight(nproma,nlev)

! new variables for DONNER
  REAL(dp) :: temp(nproma,1,nlev)
  REAL(dp) :: spec_hum(nproma,1,nlev)
  REAL(dp) :: zpress(nproma,1,nlev)
  REAL(dp) :: zpressi(nproma,1,nlev+1)
  REAL(dp) :: omega(nproma,1,nlev)
  REAL(dp) :: cell_cld_frac(nproma,1,nlev), cell_liq_amt(nproma,1,nlev),  &
              cell_liq_size(nproma,1,nlev), cell_ice_amt(nproma,1,nlev),  &
              cell_ice_size(nproma,1,nlev), meso_cld_frac(nproma,1,nlev), & 
              meso_liq_amt(nproma,1,nlev),  meso_liq_size(nproma,1,nlev), &
              meso_ice_amt(nproma,1,nlev),  meso_ice_size(nproma,1,nlev)
  REAL(dp) :: qlin(nproma,1,nlev), qiin(nproma,1,nlev), qain(nproma,1,nlev)
  REAL(dp) :: delta_ql(nproma,1,nlev), delta_qi(nproma,1,nlev),      &
              delta_qa(nproma,1,nlev), delta_temp(nproma,1,nlev),    &
              delta_vapor(nproma,1,nlev),                            &
              donner_humidity_area(nproma,1,nlev),                   &
              donner_humidity_ratio(nproma,1,nlev)
  REAL(dp) :: tr_flux(nproma,1,ntrac), &
              tracers(nproma,1,nlev,ntrac), qtrtnd(nproma,1,nlev,ntrac) 
  REAL(dp) :: sfc_sh_flux(nproma,1), sfc_vapor_flux(nproma,1), precip(nproma,1)
  REAL(dp) :: land(nproma,1), TIME(5)
  REAL(dp) :: mtot(nproma,1,nlev), detf(nproma,1,nlev), &
              uceml_inter(nproma,1,nlev+1)
  ! op_mm_20140110+ added mlev, kmin, kmax, substr
  INTEGER  :: mlev, kmin, kmax, ktop
  CHARACTER(LEN=*), PARAMETER::substr='convect_CONVEC'
  ! op_mm_20140110- 
  INTEGER  :: nsum(nproma,1)
  LOGICAL  :: conv_active(nproma,1)
  LOGICAL  :: zlookupoverflow = .FALSE. ! op_pj_20141020

#ifndef ECHAM5
  REAL(dp), DIMENSION(:,:), POINTER :: app1  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: aphp1 => NULL()
#endif

#ifdef MESSYTENDENCY
   REAL(kind=dp),DIMENSION(nproma,nlev)        :: lo_tte, lo_qte
   REAL(kind=dp),DIMENSION(nproma,nlev)        :: lo_xlte,lo_xite
   REAL(kind=dp),DIMENSION(nproma,nlev)        :: lo_vom, lo_vol
   REAL(kind=dp),DIMENSION(nproma,nlev,ntrac)  :: lo_xtte
#endif

  INTRINSIC :: MAXVAL, NINT, REAL, SIGN, SQRT

  !  Executable statements 
  kbdim  = nproma
  klev   = nlev
  klevp1 = nlevp1
  ktrac  = ntrac
  klevm1 = nlevm1

! op_mm_20140110+
#ifdef ECHAM5
  IF (ltdiag) THEN
    ! prepare fields for CUCALL increment (massflux)
    pdiga5(1:kproma,:,jrow)  = pdiga5(1:kproma,:,jrow)  - vom_3d(1:kproma,:,jrow)
    pdiga10(1:kproma,:,jrow) = pdiga10(1:kproma,:,jrow) - vol_3d(1:kproma,:,jrow)
    pdiga18(1:kproma,:,jrow) = pdiga18(1:kproma,:,jrow) -tte_3d(1:kproma,:,jrow)
  ENDIF
#endif

!  for former convect_tables that were used in the core of the TIEDTKE Convection
 
  lookupoverflow = .FALSE.

! op_mm_20131217
! nullifying channel object pointer
    conv_type(:,jrow) = 0.0_dp
    conv_top(:,jrow) = 0.0_dp
    conv_bot(:,jrow) = 0.0_dp
    cu_bot(:,jrow) = 0.0_dp
    cu_top(:,jrow) = 0.0_dp
    cu_freeze(:,jrow) = 0.0_dp
    cu_uvelo(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    base_f1(:,jrow)  = 0.0_dp
    base_f2(:,jrow)  = 0.0_dp
    base_f3(:,jrow)  = 0.0_dp
    base_f4(:,jrow)  = 0.0_dp
    massfu(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    massfd(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    u_entr(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    u_detr(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    d_entr(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    d_detr(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    massfu_asc(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    massfd_draf(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    cv_precflx(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    cv_snowflx(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    cv_precnew(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    cv_snownew(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    conv_tte(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    conv_qte(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    conv_lte(_RI_XYZ__(:,jrow,:)) = 0.0_dp 
    conv_ite(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    conv_ute(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    conv_vte(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    cv_cover(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    ! mz_jd_20161011+
    qstdev(_RI_XYZ__(:,jrow,:))             = 0.0_dp
    qskew(_RI_XYZ__(:,jrow,:))              = 0.0_dp
    Q2(_RI_XYZ__(:,jrow,:))                 = 0.0_dp
    qsat(_RI_XYZ__(:,jrow,:))               = 0.0_dp
    qp1(_RI_XYZ__(:,jrow,:))                = 0.0_dp
    cv_cover_sikma(_RI_XYZ__(:,jrow,:))     = 0.0_dp
    ! mz_jd_20161011-
    cth(:,jrow)        = 0.0_dp
    ! mz_ak_20060524+
    cv_cldwater(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    del_liqwat(_RI_XYZ__(:,jrow,:))  = 0.0_dp
    ! mz_ak_20060524-
    cv_lwc(_RI_XYZ__(:,jrow,:))      = 0.0_dp
    cv_iwc(_RI_XYZ__(:,jrow,:))      = 0.0_dp
    cv_rform(_RI_XYZ__(:,jrow,:))    = 0.0_dp
    cv_sform(_RI_XYZ__(:,jrow,:))    = 0.0_dp
    ! mim_sb_20091009+
    IF (l_lgmc) THEN
       ttp1_gp(:,:,jrow)     = 0.0_dp 
       ! op_sb_20130417+
       ptu_gp(:,:,jrow)      = 0.0_dp
       ptd_gp(:,:,jrow)      = 0.0_dp
       ! op_sb_20130417-
    END IF
    ! op_mm_20140324+
    tketconv(_RI_XYZ__(:,jrow,:))= 0.0_dp 
! op_pj_20140630+
!!$    vgustcon(1:kproma,jrow) =0.0_dp
    vgustcon(:,jrow)    = 0.0_dp
! op_pj_20140630-
    ! op_mm_20140324-

    ! op_mm_20140131+
    ! 2d pointer for using channel objects in core layer
     del_liqwat_2d  => del_liqwat(_RI_XYZ__(:,jrow,:))
     massfu_asc_2d  => massfu_asc(_RI_XYZ__(:,jrow,:))
     massfd_draf_2d => massfd_draf(_RI_XYZ__(:,jrow,:))
     u_entr_2d      => u_entr(_RI_XYZ__(:,jrow,:)) 
     u_detr_2d      => u_detr(_RI_XYZ__(:,jrow,:))
     udetr_h_2d     => udetr_h(_RI_XYZ__(:,jrow,:))
     d_entr_2d      => d_entr(_RI_XYZ__(:,jrow,:))
     d_detr_2d      => d_detr(_RI_XYZ__(:,jrow,:))
     cv_lwc_2d      => cv_lwc(_RI_XYZ__(:,jrow,:))     
     cv_iwc_2d      => cv_iwc(_RI_XYZ__(:,jrow,:))   
     cv_rform_2d    => cv_rform(_RI_XYZ__(:,jrow,:))  
     cv_sform_2d    => cv_sform(_RI_XYZ__(:,jrow,:)) 
     conv_tte_2d    => conv_tte(_RI_XYZ__(:,jrow,:))
     conv_qte_2d    => conv_qte(_RI_XYZ__(:,jrow,:))
     cv_precnew_2d  => cv_precnew(_RI_XYZ__(:,jrow,:))
     cv_snownew_2d  => cv_snownew(_RI_XYZ__(:,jrow,:))
     cv_precflx_2d  => cv_precflx(_RI_XYZ__(:,jrow,:))
     cv_snowflx_2d  => cv_snowflx(_RI_XYZ__(:,jrow,:))
     cu_uvelo_2d    => cu_uvelo(_RI_XYZ__(:,jrow,:))
     cv_cldwater_2d => cv_cldwater(_RI_XYZ__(:,jrow,:)) 
    ! 2d pointer for using channel objects in core layer
    conv_top_1d     => conv_top(:,jrow)
    conv_bot_1d     => conv_bot(:,jrow)
    cu_top_1d       => cu_top(:,jrow)
    cu_bot_1d       => cu_bot(:,jrow)
    cu_freeze_1d    => cu_freeze(:,jrow)
    cu_bot_mid_1d   => cu_bot_mid(:,jrow)
    cu_top_mid_1d   => cu_top_mid(:,jrow)
    cu_freeze_mid_1d=> cu_freeze_mid(:,jrow)
    base_f1_1d      => base_f1(:,jrow)
    base_f2_1d      => base_f2(:,jrow)
    base_f3_1d      => base_f3(:,jrow)
    base_f4_1d      => base_f4(:,jrow)
    ! op_mm_20140131-
#ifndef ECHAM5
    app1 => press_3d(_RI_XYZ__(:,jrow,:))
    aphp1 => pressi_3d(_RI_XYZ__(:,jrow,:))
#endif

    IF (l_lgmc_diag) THEN
       conv_pxtec(_RI_XYZ__(:,jrow,:))  = 0.0_dp
       conv_pqtec(_RI_XYZ__(:,jrow,:))  = 0.0_dp
       conv_tte_up_cond(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       conv_tte_up_freeze(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       conv_tte_do_verd(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       conv_tte_do_melt(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       conv_tte_up(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       conv_tte_do(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       conv_tte_su(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       conv_tte_ev(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       ! op_mm_20140131+
       ! 2d pointer for using channel objects in core layer
       conv_tte_do_melt_2d   =>  conv_tte_do_melt(_RI_XYZ__(:,jrow,:))
       conv_tte_do_verd_2d   =>  conv_tte_do_verd(_RI_XYZ__(:,jrow,:))
       conv_tte_up_freeze_2d =>  conv_tte_up_freeze(_RI_XYZ__(:,jrow,:))
       conv_tte_ev_2d        =>  conv_tte_ev(_RI_XYZ__(:,jrow,:))
       conv_tte_up_cond_2d   =>  conv_tte_up_cond(_RI_XYZ__(:,jrow,:))
       conv_tte_up_2d        =>  conv_tte_up(_RI_XYZ__(:,jrow,:)) 
       conv_tte_su_2d        =>  conv_tte_su(_RI_XYZ__(:,jrow,:))
       conv_tte_do_2d        =>  conv_tte_do(_RI_XYZ__(:,jrow,:))
       ! op_mm_20140131-
    END IF
    ! mim_sb_20091009-

    counter_deep(:,jrow) = 0.0_dp
    counter_shal(:,jrow) = 0.0_dp
    counter_midl(:,jrow) = 0.0_dp

    cth_deep(:,jrow) = 0.0_dp
    cth_shal(:,jrow) = 0.0_dp
    cth_midl(:,jrow) = 0.0_dp
! op_mm_20131217
    massfu_deep(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    massfu_shal(_RI_XYZ__(:,jrow,:)) = 0.0_dp
    massfu_midl(_RI_XYZ__(:,jrow,:)) = 0.0_dp

! mz_ht_20040316+
  if (.not.lconvection) RETURN
 
  vervel => vervel_3d(_RI_XYZ__(:,jrow,:))
 ! pxtec => xtec(:,:,jrow) ! op_mm_20140124
  pxtec => xtec(_RI_XYZ__(:,jrow,:))

  zmfu => massfu(_RI_XYZ__(:,jrow,:))
  zmfd => massfd(_RI_XYZ__(:,jrow,:))

  paprc => aprc(:,jrow)
  paprs => aprs(:,jrow)
  ptopmax => topmax(:,jrow)
  
  pqhfla => pqflux(:,jrow)
  ALLOCATE(jlab(nproma,klev))
  kmin=1
  kmax=klev
  jlab(1:kproma,1:klev) = NINT(pilab(_RI_XYZ__(1:kproma,jrow,kmin:kmax)))
#ifdef CESM1
  jlab(1:kproma,1:klev) = 0
#endif

  ktype(1:kproma) = 0
!-----------------------------------------------------------------------
!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*                 -----------------------------------
!
!

#ifdef MESSYTENDENCY
  ! tendency budget
  call mtend_get_start_l (mtend_id_t, v0 = ztp1)
  call mtend_get_start_l (mtend_id_q, v0 = zqp1)
  ! ub_ak_20190613+
  ! call mtend_get_start_l (mtend_id_tracer, v0t = zxtp1)
  DO jt = 1, ntrac
     call mtend_get_start_l (jt, v0 = zxtp1(:,:,jt))
  END DO
  ! ub_ak_20190613-
 
  
  if (convect_param .ge. 4 .and. convect_param .le. 6)then! .or. convect_param .eq. 8) then
     call mtend_get_start_l (mtend_id_xl, v0 = zxlp1) 
     call mtend_get_start_l (mtend_id_xi, v0 = zxip1) 
  endif

  if (convect_param .ge. 1 .and. convect_param .le. 4 .or. convect_param .ge. 6 .and. convect_param .le. 7) then
     call mtend_get_start_l (mtend_id_u, v0 = zup1)  
     call mtend_get_start_l (mtend_id_v, v0 = zvp1) 
  endif
#endif

100 CONTINUE
  ztmst=time_step_len
  DO 120 jk=1,klev
     zlookupoverflow = .FALSE. ! op_pj_20170327
     DO 110 jl=1,kproma
#ifndef MESSYTENDENCY
! op_m_20140124
!!        ztp1(jl,jk)=tm1(jl,jk,jrow)+tte(jl,jk)*ztmst
          ztp1(jl,jk)=tm1(_RI_XYZ__(jl,jrow,jk))+tte_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
#endif
        IF (L_LGMC) ttp1_gp(jl,jk,jrow)=ztp1(jl,jk) ! mim_sb_20091009
#ifndef MESSYTENDENCY
!        zqp1(jl,jk)=MAX(0._dp,qm1(jl,jk,jrow)+qte(jl,jk)*ztmst) ! op_mm_20140124
!        zxlp1(jl,jk)=xlm1(jl,jk,jrow)+xlte(jl,jk)*ztmst
!        zxip1(jl,jk)=xim1(jl,jk,jrow)+xite(jl,jk)*ztmst
        zqp1(jl,jk)=MAX(0._dp,qm1(_RI_XYZ__(jl,jrow,jk))+qte_3d(_RI_XYZ__(jl,jrow,jk))*ztmst)
        qp1(_RI_XYZ__(jl,jrow,jk))=zqp1(jl,jk) ! mz_jd_20170201
        zxlp1(jl,jk)=&
             xlm1(_RI_XYZ__(jl,jrow,jk))+xlte_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
        zxip1(jl,jk)=&
             xim1(_RI_XYZ__(jl,jrow,jk))+xite_3d(_RI_XYZ__(jl,jrow,jk))*ztmst

#else
        zqp1(jl,jk)=MAX(0._dp,zqp1(jl,jk))
        qp1(_RI_XYZ__(jl,jrow,jk))=zqp1(jl,jk) ! mz_jd_20170201
        !!In case variable (xl,xi) not changed by conv scheme
        if (convect_param .ge. 1 .and. convect_param .le. 3 .or. convect_param .eq. 7) then! .or. convect_param .eq. 8 ) then
!           zxlp1(jl,jk)=xlm1(jl,jk,jrow)+xlte(jl,jk)*ztmst ! op_mm_20140124
!           zxip1(jl,jk)=xim1(jl,jk,jrow)+xite(jl,jk)*ztmst ! smilification
           zxlp1(jl,jk)=&
                xlm1(_RI_XYZ__(jl,jrow,jk))+xlte_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
           zxip1(jl,jk)=&
                xim1(_RI_XYZ__(jl,jrow,jk))+xite_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
        endif
#endif
        zxp1(jl,jk)=MAX(0._dp,zxlp1(jl,jk)+zxip1(jl,jk))
        ztvp1(jl,jk)=ztp1(jl,jk)*(1._dp+vtmpc1*zqp1(jl,jk)-zxp1(jl,jk))
#ifndef MESSYTENDENCY
!        zup1(jl,jk)=um1(jl,jk,jrow)+vom(jl,jk)*ztmst! op_mm_20140124
!        zvp1(jl,jk)=vm1(jl,jk,jrow)+vol(jl,jk)*ztmst! smilification
         zup1(jl,jk)=um1(_RI_XYZ__(jl,jrow,jk))+vom_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
         zvp1(jl,jk)=vm1(_RI_XYZ__(jl,jrow,jk))+vol_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
#else
        !!In case variable (u,v) not changed by conv scheme
        if (convect_param .eq. 5) then! .or. convect_param .eq. 8) then
!           zup1(jl,jk)=um1(jl,jk,jrow)+vom(jl,jk)*ztmst ! op_mm_20140124
!           zvp1(jl,jk)=vm1(jl,jk,jrow)+vol(jl,jk)*ztmst
           zup1(jl,jk)=um1(_RI_XYZ__(jl,jrow,jk))+vom_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
           zvp1(jl,jk)=vm1(_RI_XYZ__(jl,jrow,jk))+vol_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
        endif
#endif
        it = INT(ztp1(jl,jk)*1000.)
! ka_sv_20170406+: exclude ionosphere/thermosphere
!!$     IF (it<jptlucu1 .OR. it>jptlucu2) zlookupoverflow = .TRUE.
        IF ( (it<jptlucu1 .OR. it>jptlucu2) .AND. &
             (app1(jl,jk) >= 1.0_dp) ) zlookupoverflow = .TRUE.
! ka_sv_20170406-
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat(jl,jk)=tlucua(it)/app1(jl,jk)
        zqsat(jl,jk)=MIN(0.5_dp,zqsat(jl,jk))
        zqsat(jl,jk)=zqsat(jl,jk)/(1._dp-vtmpc1*zqsat(jl,jk))
        qsat(_RI_XYZ__(jl,jrow,jk)) = zqsat(jl,jk) ! mz_jd_20170201

! op_mm_20140124 (smilification)
       !! zrhoa(jl,jk) = grmass (jl,jk,jrow) / grvol(jl,jk,jrow)
        zrhoa(jl,jk) = grmass (_RI_XYZ__(jl,jrow,jk)) / grvol(_RI_XYZ__(jl,jrow,jk))
        zwp1(jl,jk,1)= zxlp1(jl,jk)
        zwp1(jl,jk,2)= zxip1(jl,jk)
110  END DO

     IF (zlookupoverflow) THEN 
        do jl=1,kproma
           if ( INT(ztp1(jl,jk)*1000.) <jptlucu1 .OR. &
                INT(ztp1(jl,jk)*1000.) >jptlucu2)     &
                ! op_mm_20140124 (smilification)
                ! print*, jk, jl,ztp1(jl,jk)*1000., tm1(jl,jk,jrow), &
                write(*,*) 'lookupoverflow', jl, jk, jrow  &
                , ztp1(jl,jk)*1000., tm1(_RI_XYZ__(jl,jrow,jk))       &
                , tte_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
        enddo
! op_mm_20140110 finish -> error_bi
        CALL error_bi('convect_convec - lookuperror', substr)
     ENDIF

#ifndef MESSYTENDENCY
     DO 1104 jt=1,ktrac
        DO 1102 jl=1,kproma
       !    zxtp1(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst   ! op_mm_20140327 smilification
            zxtp1(jl,jk,jt)=pxtm1(_RI_X_ZN_(jl,jk,jt))+pxtte(_RI_X_ZN_(jl,jk,jt))*ztmst
1102    END DO
1104 END DO
#endif

120 END DO
  DO 130 jl=1,kproma
     zrain(jl)=0.
     locum(jl)=.FALSE.
130 END DO
!
! mz_ht_20040317+
! 1d-field over all tracers to save their nconv switch that it can be used 
! in the core routines of each convection scheme
  ztracconv(:) = 0
  do jt=1,ntrac
     IF (ti_gp(jt)%tp%meta%cask_i(I_convect) == ON) ztracconv(jt) = 1
  enddo
! mz_ht_20040317-

!
!-----------------------------------------------------------------------
!
!*    2.     CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
!*           -----------------------------------------------------------
!
!
200 CONTINUE

#ifdef MESSYTENDENCY
  ! Prepare Tiedtke and ECMWF schemes
  if(convect_param .ge. 1 .and. convect_param .le. 4)   then
     !kk this setting to zero should be obsolete. Code does not work when they are
     !kk mising. To be testet!
     lo_tte  = 0.0_dp
     lo_vom  = 0.0_dp
     lo_vol  = 0.0_dp
     lo_qte  = 0.0_dp
     lo_xtte = 0.0_dp

     ! for case 4 (ECMWF scheme)
     if(convect_param .eq. 4 ) then
        lo_xlte = 0.0_dp
        lo_xite = 0.0_dp
        lo_xlte(1:kproma,:)    = xlte_3d(_RI_XYZ__(1:kproma,jrow,:))
        lo_xite(1:kproma,:)    = xite_3d(_RI_XYZ__(1:kproma,jrow,:))
     endif

     lo_tte(1:kproma,:)    = tte_3d(_RI_XYZ__(1:kproma,jrow,:))
     lo_vom(1:kproma,:)    = vom_3d(_RI_XYZ__(1:kproma,jrow,:))
     lo_vol(1:kproma,:)    = vol_3d(_RI_XYZ__(1:kproma,jrow,:))
     lo_qte(1:kproma,:)    = qte_3d(_RI_XYZ__(1:kproma,jrow,:))
     lo_xtte(1:kproma,:,:) = pxtte(1:kproma,:,:)
  end if
#endif

  SELECT CASE (convect_param)
  CASE(0)
     kmin=1
     kmax=klev
     pilab(_RI_XYZ__(1:kproma,jrow,kmin:kmax)) = REAL(jlab(1:kproma,1:klev), dp)
    DEALLOCATE(jlab)
    RETURN

  CASE(1)         ! Tiedtke - Nordeng
     ! op_mm_20140124
   ! zpqtec => pqtec(:,:,jrow) ! mz_pj_20070309
     zpqtec => pqtec(_RI_XYZ__(:,jrow,:)) ! op_mm_20140327 smilification
    CALL tiedtke_cumastr(kproma, kbdim, klev, klevp1, klevm1,               &
#ifndef _rs6000
                 jrow,     ztmst,    nn,       ktrac, jlab,                 &
#else
                 jrow,     ztmst,    nn,       ktrac, jlab(1:kbdim,1:klev), &
#endif
                 ztp1,     zqp1,     zxp1,     zup1,   zvp1,                &
                 ! ub_ak_20190307+
                 !ztvp1,    loland,                                         &
                 ztvp1,    loland_2d(:,jrow),                               &
                 ! ub_ak_20190307-
                 zxtp1,    zxtu,     pxtte,                                 &
#ifndef CESM1
                 vervel,   zqsat,    pqhfla,                                &
#else
                 vervel,   zqsat,    pqhfla*.905,                           &
#endif
                 aphp1,    geopot_3d(_RI_XYZ__(:,jrow,:)),                  &
                 tte_3d(_RI_XYZ__(:,jrow,:)),  qte_3d(_RI_XYZ__(:,jrow,:)), &
                 vom_3d(_RI_XYZ__(:,jrow,:)) , vol_3d(_RI_XYZ__(:,jrow,:)), &
                 rsfc_2d(:,jrow),  ssfc_2d(:,jrow),  paprc,  paprs,  pxtec, &
                 zpqtec,                                                    &
                 locum,    ktype,    icbot,    ictop,                       &
                 zmfu,     zmfd,     zrain,    ztracconv,                   &
                 delta_time,  &
                 vgustcon(1:kbdim,jrow),tketconv(_RI_XYZ__(1:kbdim,jrow,:))) 

  CASE(2)         ! Tiedtke
     zpqtec => pqtec(_RI_XYZ__(:,jrow,:))
    CALL tiedtke_cumastrt(kproma, kbdim, klev, klevp1, klevm1,              &
#ifndef _rs6000
                 jrow,     ztmst,    nn,       ktrac,  jlab,                &
#else
                 jrow,     ztmst,    nn,       ktrac,  jlab(1:kbdim,1:klev),&
#endif
                 ztp1,     zqp1,     zxp1,     zup1,   zvp1,                &
                 ! ub_ak_20190307+
                 !ztvp1,    loland,                                         &
                 ztvp1,    loland_2d(:,jrow),                               &
                 ! ub_ak_20190307-
                 zxtp1,    zxtu,     pxtte,                                 &
                 vervel,   zqsat,    pqhfla,                                &
                 aphp1,    geopot_3d(_RI_XYZ__(:,jrow,:)),                  &
                 tte_3d(_RI_XYZ__(:,jrow,:)), qte_3d(_RI_XYZ__(:,jrow,:)),  &
                 vom_3d(_RI_XYZ__(:,jrow,:)), vol_3d(_RI_XYZ__(:,jrow,:)),  &
                 rsfc_2d(:,jrow), ssfc_2d(:,jrow), paprc, paprs, pxtec,     &
                 zpqtec,                                                    &
                 locum,    ktype,    icbot,    ictop,                       &
                 zmfu,     zmfd,     zrain,    ztracconv,                   &
                 delta_time, &
                 vgustcon(1:kbdim,jrow),tketconv(_RI_XYZ__(1:kbdim,jrow,:)))   

  CASE(3)         ! Tiedtke - Hybrid
     zpqtec => pqtec(_RI_XYZ__(:,jrow,:))
    CALL tiedtke_cumastrh(kproma, kbdim, klev, klevp1, klevm1,              &
#ifndef _rs6000
                 jrow,     ztmst,    nn,       ktrac,  jlab,                &
#else
                 jrow,     ztmst,    nn,       ktrac,  jlab(1:kbdim,1:klev),&
#endif
                 ztp1,     zqp1,     zxp1,     zup1,   zvp1,                &
                 ! ub_ak_20190307+
                 !ztvp1,    loland,                                         &
                 ztvp1,    loland_2d(:,jrow),                               &
                 ! ub_ak_20190307-
                 zxtp1,    zxtu,     pxtte,                                 &
                 vervel,   zqsat,    pqhfla,                                &
                 aphp1,   geopot_3d(_RI_XYZ__(:,jrow,:)),                   &
                 tte_3d(_RI_XYZ__(:,jrow,:)),  qte_3d(_RI_XYZ__(:,jrow,:)), &
                 vom_3d(_RI_XYZ__(:,jrow,:)),  vol_3d(_RI_XYZ__(:,jrow,:)), &
                 rsfc_2d(:,jrow), ssfc_2d(:,jrow), paprc,    paprs,  pxtec, &
                 zpqtec,                                                    &
                 locum,    ktype,    icbot,    ictop,                       &
                 zmfu,     zmfd,     zrain,    ztracconv,                   &
                 delta_time, &
                 vgustcon(1:kbdim,jrow),tketconv(_RI_XYZ__(1:kbdim,jrow,:)) ) 

  CASE(4)         ! ECMWF (based on Tiedtke)   

    psstru(:) = 0._dp
    psstrv(:) = 0._dp
    PTU(:,:)  = 0._dp
    PQU(:,:)  = 0._dp
    PLU(:,:)  = 0._dp  
    THFLX_EC(1:kproma,:) = 0.0_dp       ! 3d-field of sensible heatflux 
                                        ! (only lowest level (nlev+1) is used)
    THFLX_EC(1:kproma,nlev+1) = s_heatflux(1:kproma,jrow)
    QFLX_EC(1:kproma,:) = 0.0_dp        ! 3d-field of moisture flux 
                                        ! (only lowest level (nlev+1) is used)
    QFLX_EC(1:kproma,nlev+1) = pqhfla(1:kproma)

    sum_wat1(:) = 0._dp
    do jk=1,nlev
      do jl=1,kproma
      ! pressure difference between interface layers
        pdel(jl,jk)  = aphp1(jl,jk+1) - aphp1(jl,jk)
        sum_wat1(jl)  = sum_wat1(jl)                                   + &
                        (qte_3d(_RI_XYZ__(jl,jrow,jk)) + xlte_3d(_RI_XYZ__(jl,jrow,jk)) + xite_3d(_RI_XYZ__(jl,jrow,jk))) * &
                        pdel(jl,jk) / g   
!        zqp1(jl,jk)=qm1(jl,jk,jrow)+qte(jl,jk)*ztmst ! op_mm_20140124
        zqp1(jl,jk)=qm1(_RI_XYZ__(jl,jrow,jk))+qte_3d(_RI_XYZ__(jl,jrow,jk))*ztmst
        qp1(_RI_XYZ__(jl,jrow,jk))=zqp1(jl,jk) ! mz_jd_20170201
        zxp1(jl,jk)=zxlp1(jl,jk)+zxip1(jl,jk)
      enddo
    enddo
    CALL ECMWF_cumastr(1,       KPROMA,   KBDIM,    1,                   &
         ! ub_ak_20190307+
         !              NLEV,    LOLAND,   ZTMST,                        &
                       NLEV,    LOLAND_2D(:,jrow),   ZTMST,              &
         ! ub_ak_20190307-
                       ztp1,    zqp1,     zup1,     zvp1,   zxp1,        &
                       VERVEL,  QFLX_EC,  THFLX_EC,                      &
                       psstru,  psstrv,                                  &
                       app1,    aphp1,                                   &
                       geopot(_RI_XYZ__(:,jrow,:)),                      &
                       geopoti(_RI_XYZ__(:,jrow,:)),                     & 
                       tte_3d(_RI_XYZ__(:,jrow,:)),                      &
                       QTE_3d(_RI_XYZ__(:,jrow,:)),                      &
                       VOM_3d(_RI_XYZ__(:,jrow,:)),                      &
                       VOL_3d(_RI_XYZ__(:,jrow,:)),          &
                       xite_3d(_RI_XYZ__(:,jrow,:)),                     &
                       xlte_3d(_RI_XYZ__(:,jrow,:)),                     &
                       LOCUM,   KTYPE,    ICBOT,    ICTOP,               &
                       KBOTSC,  LOSHAL,                                  &
                       PTU,     PQU,      PLU,                           &
                       PFLX,    PSFLX,                                   &
                       PRAIN,                                            &
                       MASSFU(_RI_XYZ__(:,jrow,:)), MASSFD(_RI_XYZ__(:,jrow,:)),             & ! op_mm_20140327
                       U_detr(_RI_XYZ__(:,jrow,:)), d_detr(_RI_XYZ__(:,jrow,:)),             & ! smilification
                       U_ENTR(_RI_XYZ__(:,jrow,:)), D_ENTR(_RI_XYZ__(:,jrow,:)),             &
                       CAPE(:,jrow),                                     &
                       NTRAC,    ZXTP1,     PXTTE,                       &
                       cv_lwc(_RI_XYZ__(:,jrow,:)),   cv_iwc(_RI_XYZ__(:,jrow,:)),             &
                       cv_rform(_RI_XYZ__(:,jrow,:)), cv_sform(_RI_XYZ__(:,jrow,:)) )
    
    do jk=1,nlev
      do jl=1,kproma
! op_mm_20140109 RI
        cv_precflx(_RI_XYZ__(jl,jrow,jk)) = pflx(jl,jk+1)
!        cv_snowflx(_RI_XYZ__(jl,jrow,:)) = psflx(jl,jk+1)
        cv_snowflx(_RI_XYZ__(jl,jrow,jk)) = psflx(jl,jk+1)
      enddo
    enddo
    ! WRITE OUT CONVECTION TYPE, bottom and top levels of Convection 
    ! TO CONVECT CHANNEL
    conv_type(1:kproma,jrow) = REAL(ktype(1:kproma),dp)  
    conv_bot(1:kproma,jrow)  = REAL(icbot(1:kproma),dp)
    conv_top(1:kproma, jrow) = REAL(ictop(1:kproma),dp) 
    do jl=1,kproma

! op_mm_20140109 - RI
      rsfc_2d(jl,jrow) = cv_precflx(_RI_XYZ__(jl,jrow,nlev))
      ssfc_2d(jl,jrow) = cv_snowflx(_RI_XYZ__(jl,jrow,nlev))
      aprc(jl,jrow) = aprc(jl,jrow)  &
           + delta_time * (rsfc_2d(jl,jrow) + ssfc_2d(jl,jrow))
      aprs(jl,jrow) = aprs(jl,jrow) + delta_time * ssfc_2d(jl,jrow)
!      values for NOx lightning

      cu_bot(jl,jrow)    = conv_bot(jl,jrow)
      cu_top(jl,jrow)    = conv_top(jl,jrow)
      cu_freeze(jl,jrow) = conv_top(jl,jrow)
      IF (NINT(cu_top(jl,jrow)) > 0._dp) THEN
        do jk=nint(cu_top(jl,jrow)),nint(cu_bot(jl,jrow))
          if (ztp1(jl,jk).le.273.15_dp) cu_freeze(jl,jrow)=REAL(jk,dp)
        enddo
        do jk=nint(cu_top(jl,jrow)),nint(cu_bot(jl,jrow))
          cu_uvelo(_RI_XYZ__(jl,jrow,jk))=massfu(_RI_XYZ__(jl,jrow,jk))/zrhoa(jl,jk)
        enddo
      endif
     
!    end of values for NOx lightning
    enddo

    sum_wat2(:) = 0._dp
    do jk=1,nlev
      do jl=1,kproma
        sum_wat2(jl)  = sum_wat2(jl)                              + &
                        (qte_3d(_RI_XYZ__(jl,jrow,jk)) +&
                        xlte_3d(_RI_XYZ__(jl,jrow,jk)) +&
                        xite_3d(_RI_XYZ__(jl,jrow,jk)))*            &
                        pdel(jl,jk) / g                   +         &
                        (rsfc_2d(jl,jrow) + ssfc_2d(jl,jrow))
      enddo
    enddo
!!$    do jl=1,kproma
!!$      if (ABS(sum_wat2(jl) - sum_wat1(jl)) > 1.e-10_dp) THEN
!!$        print*, sum_wat2(jl) - sum_wat1(jl), sum_wat1(jl), sum_wat2(jl), &
!!$              rsfc(jl), ssfc(jl), sum_wat2(jl) - (rsfc(jl) + ssfc(jl))
!!$      ENDIF
!!$    enddo


  CASE(5)         ! ZHANG - HACK - McFarlane
     do jk=1,nlev
      do jl=1,kproma
      ! pressure difference between interface layers
        pdel(jl,jk)  = aphp1(jl,jk+1) - aphp1(jl,jk)
        rpdel(jl,jk) = 1._dp/pdel(jl,jk)
! op_mm_20140124
!        zi(jl,jk+1)  = geopoti(jl,jk,jrow)/g
!        zm(jl,jk)    = geopot(jl,jk,jrow)/g
        zi(jl,jk+1)  = geopoti(_RI_XYZ__(jl,jrow,jk))/g
        zm(jl,jk)    = geopot(_RI_XYZ__(jl,jrow,jk))/g
      enddo
    enddo
    
    do jl=1,kproma
      zi(jl,1) = 0._dp
      ts(jl)   = ztp1(jl,klev)
! op_mm_20140124
!      phis(jl) = geopoti(jl,klev, jrow)
      mlev=klev
      phis(jl) = geopoti(_RI_XYZ__(jl,jrow,mlev)) ! um_ak_20140423
!     input values from PBL calculation
      theta_surf(jl) = tsurf_2d(jl,jrow)  & 
                        *(1.E5_dp/aphp1(jl,klev+1))**(rgas/cp)
 
      if (zsenkf_2d(jl,jrow).gt.0.0_dp) then

         thvsrf(jl) = theta_surf(jl)*(1.0_dp + 0.61_dp*zqp1(jl,klev))
         obklen(jl) = -thvsrf(jl)*zust_2d(jl,jrow)**3/                   &
                      (g*0.4_dp*(zsenkf_2d(jl,jrow) +                    &
                      sign(1.e-10_dp,zsenkf_2d(jl,jrow))))

         phiminv(jl) = (1._dp - 1.5_dp*pblh(jl,jrow)/obklen(jl))**(1._dp/3._dp)
         phihinv(jl) = sqrt(1._dp - 1.5_dp*pblh(jl,jrow)/obklen(jl))
         
         tpert(jl) = max(zsenkf_2d(jl,jrow)*8.5_dp/(zust_2d(jl,jrow) * &
           phiminv(jl)),0._dp) 
         qpert(jl) = max(zlatkf_2d(jl,jrow)*8.5_dp/(zust_2d(jl,jrow) * &
           phiminv(jl)),0._dp) 

         do jk=1,nlev
           wpert3d(jl,jk) = zust_2d(jl,jrow)*phihinv(jl) *&
             max(0._dp,1._dp-zm(jl,jk)/(2._dp*pblh(jl,jrow)))
         enddo

      ELSE

        tpert(jl) = max(zsenkf_2d(jl,jrow)*8.5_dp/zust_2d(jl,jrow),0._dp) 
        qpert(jl) = max(zlatkf_2d(jl,jrow)*8.5_dp/zust_2d(jl,jrow),0._dp) 
        do jk=1,nlev
          wpert3d(jl,jk) = zust_2d(jl,jrow) * &
            max(0._dp,1._dp-zm(jl,jk)/(2._dp*pblh(jl,jrow)))
        enddo

      ENDIF
     
    enddo
  
    omga => vervel   !omega(:,:,jrow) , vertical velocity

    CALL zhang_cumastr(kproma,  kbdim,  klev,   klevp1,   2,         &
                       jrow,    nstep,  ztmst,  app1,     aphp1,     &
                       pdel,    rpdel,  zm,     zi,       tpert,     &
                       qpert,   phis,   ts,     pblh(:,jrow),        &
                       ztp1,    zqp1,   zwp1,                        &
                       !,cmfmc   , &
                       !  zdu     ,& 
          !             cmfdqr  ,conicw  ,&
               !        precc   ,&
                       conv_top(1:kproma,jrow), conv_bot(1:kproma,jrow),  &
!                       xtec(:,:,jrow)     , & ! op_mm_20140124
                       xtec(_RI_XYZ__(:,jrow,:))     , &
                       pflx    ,psflx,                               &
                       omga    ,wpert3d ,                            &
!                       zmu     , chembgt ,                           &
!-mgl
                      !   zmug    ,zmdg     ,zdug      ,zeug    ,       &  
                      !   zedg    ,zdpg     ,&
                      ! dsubcld   ,&!zjtg    ,zjbg , &
                       ideep   &!
                       !,lengath&
                       )
!   calculate tendencies for t, q, xl, xi, tracer

    do jk=1,klev
      do jl=1,kproma
        
        conv_tte(_RI_XYZ__(jl,jrow,jk)) = (ztp1(jl,jk) - &
        !                       (tm1(jl,jk,jrow)+tte(jl,jk)*ztmst)) / ztmst ! op_mm_20140124
             (tm1(_RI_XYZ__(jl,jrow,jk))+tte_3d(_RI_XYZ__(jl,jrow,jk))*ztmst)) / ztmst
        tte_3d(_RI_XYZ__(jl,jrow,jk)) = &
             tte_3d(_RI_XYZ__(jl,jrow,jk)) + conv_tte(_RI_XYZ__(jl,jrow,jk))

        conv_qte(_RI_XYZ__(jl,jrow,jk)) = (zqp1(jl,jk) - &
       !                        (qm1(jl,jk,jrow)+qte(jl,jk)*ztmst)) / ztmst ! op_mm_20140124
                                (qm1(_RI_XYZ__(jl,jrow,jk))+&
                                qte_3d(_RI_XYZ__(jl,jrow,jk))*ztmst)) / ztmst
        qte_3d(_RI_XYZ__(jl,jrow,jk)) = &
             qte_3d(_RI_XYZ__(jl,jrow,jk)) + conv_qte(_RI_XYZ__(jl,jrow,jk))
        
        conv_lte(_RI_XYZ__(jl,jrow,jk)) = (zwp1(jl,jk,1) - zxlp1(jl,jk) ) / ztmst
        xlte_3d(_RI_XYZ__(jl,jrow,jk)) = &
             xlte_3d(_RI_XYZ__(jl,jrow,jk)) + conv_lte(_RI_XYZ__(jl,jrow,jk)) 

        conv_ite(_RI_XYZ__(jl,jrow,jk)) = (zwp1(jl,jk,2) - zxip1(jl,jk) ) / ztmst
        xite_3d(_RI_XYZ__(jl,jrow,jk)) = &
             xite_3d(_RI_XYZ__(jl,jrow,jk)) + conv_ite(_RI_XYZ__(jl,jrow,jk)) 

!!$         do jt=1,ntrac
!!$           pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + (zxtp1(jl,jk,jt) -  &
!!$                             pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst) / ztmst
!!$         enddo

        IF (ztp1(jl,jk) > Tmelt) THEN
          cv_precflx(_RI_XYZ__(jl,jrow,jk)) = pflx(jl,jk+1)
        ELSE
          cv_snowflx(_RI_XYZ__(jl,jrow,jk)) = psflx(jl,jk+1)
        ENDIF

       enddo
     enddo
     do jl=1,kproma
       rsfc_2d(jl,jrow) = cv_precflx(_RI_XYZ__(jl,jrow,nlev))
       ssfc_2d(jl,jrow) = cv_snowflx(_RI_XYZ__(jl,jrow,nlev))
       aprc(jl,jrow) = aprc(jl,jrow)  &
            + delta_time * (rsfc_2d(jl,jrow) + ssfc_2d(jl,jrow))
       aprs(jl,jrow) = aprs(jl,jrow) + delta_time * ssfc_2d(jl,jrow)
   
       if (maxval(massfu(_RI_XYZ__(jl,jrow,:))).gt.1.e-15_dp) then
         conv_type(jl,jrow) = 2._dp
         if (ideep(jl).gt.0) conv_type(jl,jrow) = 1._dp
       endif

!      values for NOx lightning

       IF (ideep(jl).gt.0) then
         cu_bot(jl,jrow)    = conv_bot(jl,jrow)
         cu_top(jl,jrow)    = conv_top(jl,jrow)
         cu_freeze(jl,jrow) = conv_top(jl,jrow)
         IF (NINT(cu_top(jl,jrow)) > 0._dp) THEN
           do jk=nint(cu_top(jl,jrow)),nint(cu_bot(jl,jrow))
             if (ztp1(jl,jk).le.273.15_dp) cu_freeze(jl,jrow)=REAL(jk,dp)
           enddo
           do jk=nint(cu_top(jl,jrow)),nint(cu_bot(jl,jrow))
             cu_uvelo(_RI_XYZ__(jl,jrow,jk))=massfu(_RI_XYZ__(jl,jrow,jk))/zrhoa(jl,jk)
           enddo
         endif
       ENDIF! op_m_20140226

!    end of values for NOx lightning

     enddo

   CASE(6)             ! BECHTOLD

     zw (1:kproma,1:nlev) = -1._dp * vervel(1:kproma,1:nlev) / &
                            (g*zrhoa(1:kproma,1:nlev))
     THFLX(1:kproma,1:nlev) = 0.0_dp                     ! is not found in the model as a 3D field, but also commente
     WAT_DIAG(:,jrow)   = 0.0_dp

     THFLX(1:kproma,nlev) = s_heatflux(1:kproma,jrow)
     ccount(1:kproma)  = int(counter(1:kproma,jrow))
     icbot(1:kproma)   = int(conv_bot(1:kproma,jrow))
     ictop(1:kproma)   = int(conv_top(1:kproma,jrow))

     CALL BECHTOLD_CUMASTR(                                                    &
       ! grid dimensions
                           KPROMA, NLEV, 1, KPROMA, 1, 1, ZTMST,               &
       ! meteorological
                           app1(1:kproma,:), aphp1(1:kproma,:),                &
      !!                     geopot(1:kproma,:,jrow),                          &
                             geopot(_RI_XYZ__(1:kproma,jrow,:)),        & ! op_mm_20140327
       ! grid 
                           gboxarea_2d(1:kproma,jrow), THFLX(1:kproma,:),      &
                           ZTP1(1:kproma,:), ZQP1(1:kproma,:),                 &
                           ZXLP1(1:kproma,:), ZXIP1(1:kproma,:),               &
                           ZUP1(1:kproma,:), ZVP1(1:kproma,:), ZW(1:kproma,:), &
       ! properties                       
                           CCOUNT(1:kproma),                                   &
       ! tendencies for temp, q, 
                           CONV_TTE(_RI_XYZ__(1:kproma,jrow,:)),CONV_QTE(_RI_XYZ__(1:kproma,jrow,:)),&   ! op_mm_20140327
       ! liquid water and ice
                           CONV_LTE(_RI_XYZ__(1:kproma,jrow,:)),CONV_ITE(_RI_XYZ__(1:kproma,jrow,:)),&
       ! tendencies for precipitation
!                          PPRTEN, PPRSTEN,                                    &
       ! convective mass fluxes
                           MASSFU(_RI_XYZ__(1:kproma,jrow,:)), MASSFD(_RI_XYZ__(1:kproma,jrow,:)),   &
       ! convective precipitation fluxes
                           cv_precflx(_RI_XYZ__(1:kproma,jrow,:)),                        &
                           cv_snowflx(_RI_XYZ__(1:kproma,jrow,:)),                        &
       ! CAPE 
                           CAPE(1:kproma,jrow),                                &
       ! base and top level of convection
                           iCTOP(1:kproma), iCBOT(1:kproma),                   &
                           ideep,                                              &
       ! updraft concentrations of q, condensate
!                          PURV(1:kproma,:,jrow), PURCI(1:kproma,:,jrow),      &
       ! changes in wind by convection 
                           CONV_UTE(_RI_XYZ__(1:kproma,jrow,:)),CONV_VTE(_RI_XYZ__(1:kproma,jrow,:)),&
       ! changes in tracers by convection (if cvtrans is not used)
                           ntrac, zxtp1(1:kproma,:,:), TRACTEN(1:kproma,:,:),  &
       ! detrainment rates updraft, downdraft
                           U_detr(_RI_XYZ__(1:kproma,jrow,:)), d_detr(_RI_XYZ__(1:kproma,jrow,:)),   &
       ! entrainment rates updraft, downdraft
                           U_entr(_RI_XYZ__(1:kproma,jrow,:)), d_entr(_RI_XYZ__(1:kproma,jrow,:)),   &
       ! diagnose water correction
                           WAT_DIAG(1:kproma, jrow),                           &
                           CV_LWC(_RI_XYZ__(1:kproma,jrow,:)), CV_IWC(_RI_XYZ__(1:kproma,jrow,:)),   &
                           CV_RFORM(_RI_XYZ__(1:kproma,jrow,:)),CV_SFORM(_RI_XYZ__(1:kproma,jrow,:)))

!    values for base ECHAM

     do jl=1,kproma

       rsfc_2d(jl,jrow) = cv_precflx(_RI_XYZ__(jl,jrow,nlev))
       ssfc_2d(jl,jrow) = cv_snowflx(_RI_XYZ__(jl,jrow,nlev))
       aprc(jl,jrow) = aprc(jl,jrow) &
            + delta_time * (rsfc_2d(jl,jrow) + ssfc_2d(jl,jrow))
       aprs(jl,jrow) = aprs(jl,jrow) + delta_time * ssfc_2d(jl,jrow)
       
       CONV_TOP(jl,jrow) = REAL(ictop(jl),dp)
       CONV_BOT(jl,jrow) = REAL(icbot(jl),dp)
       COUNTER(jl,jrow)  = REAL(ccount(jl),dp)

!    new tendencies
       do jk=1,nlev
         tte_3d(_RI_XYZ__(jl,jrow,jk))  = &
              tte_3d(_RI_XYZ__(jl,jrow,jk))  + conv_tte(_RI_XYZ__(jl,jrow,jk))
         qte_3d(_RI_XYZ__(jl,jrow,jk))  = &
              qte_3d(_RI_XYZ__(jl,jrow,jk))  + conv_qte(_RI_XYZ__(jl,jrow,jk))
         xlte_3d(_RI_XYZ__(jl,jrow,jk)) = &
              xlte_3d(_RI_XYZ__(jl,jrow,jk)) + conv_lte(_RI_XYZ__(jl,jrow,jk))
         xite_3d(_RI_XYZ__(jl,jrow,jk)) = &
              xite_3d(_RI_XYZ__(jl,jrow,jk)) + conv_ite(_RI_XYZ__(jl,jrow,jk))
         vom_3d(_RI_XYZ__(jl,jrow,jk))  = &
              vom_3d(_RI_XYZ__(jl,jrow,jk))  + conv_ute(_RI_XYZ__(jl,jrow,jk))
         vol_3d(_RI_XYZ__(jl,jrow,jk))  = &
              vol_3d(_RI_XYZ__(jl,jrow,jk))  + conv_vte(_RI_XYZ__(jl,jrow,jk))

!!$         do jt=1,ntrac 
!!$           if (ztracconv(jt).eq.1) &
!!$             pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + tracten(jl,jk,jt)
!!$         enddo
       enddo
          
       if (maxval(massfu(_RI_XYZ__(jl,jrow,:))).gt.1.e-15_dp) then
         conv_type(jl,jrow) = 2._dp
         if (ideep(jl).gt.0) conv_type(jl,jrow) = 1._dp
       endif

!      values for NOx lightning
! op_mm_20131217 ->  RI
       IF (ideep(jl).gt.0) then
         cu_bot(jl,jrow)    = conv_bot(jl,jrow)
         cu_top(jl,jrow)    = conv_top(jl,jrow)
         cu_freeze(jl,jrow) = conv_top(jl,jrow)
         if (nint(cu_top(jl,jrow)) > 0) then
           do jk=nint(cu_top(jl,jrow)),nint(cu_bot(jl,jrow))
             if (ztp1(jl,jk).le.273.15_dp) cu_freeze(jl,jrow)=REAL(jk,dp)
           enddo
           do jk=nint(cu_top(jl,jrow)),nint(cu_bot(jl,jrow))
             cu_uvelo(_RI_XYZ__(jl,jrow,jk))=massfu(_RI_XYZ__(jl,jrow,jk))/zrhoa(jl,jk)
           enddo
         endif
       ENDIF

!    end of values for NOx lightning

     enddo

   CASE(7)      ! EMANUEL

     CALL EMANUEL_CUMASTR(KPROMA,  NLEV, NTRAC, ZTMST,                         &
                          app1(1:kproma,:),    aphp1(1:kproma,:),              &
                          ZTP1(1:kproma,:),    ZQP1(1:kproma,:),               &
                          ZUP1(1:kproma,:),    ZVP1(1:kproma,:),               &
                          zxtp1(1:kproma,:,:), ZQSAT(1:kproma,:),              &
                           ! OUTPUT
                          itype,   rsfc_2d(1:kproma,jrow),                     &
     ! tendencies for temp, q, u, v
                          CONV_TTE(_RI_XYZ__(1:kproma,jrow,:)), CONV_QTE(_RI_XYZ__(1:kproma,jrow,:)),&
                          CONV_UTE(_RI_XYZ__(1:kproma,jrow,:)), CONV_VTE(_RI_XYZ__(1:kproma,jrow,:)),&
                          TRACTEN(1:kproma,:,:),                               &
                          ZWD(1:kproma),             TPRIME(1:kproma),         &
                          QPRIME(1:kproma),          CBMF(1:kproma,jrow),      &
                          ICBOT(1:kproma),           ICTOP(1:kproma),          &
                          MASSFU(_RI_XYZ__(1:kproma,jrow,:)),   MASSFD(_RI_XYZ__(1:kproma,jrow,:)),  &
                          CV_PRECFLX(_RI_XYZ__(1:kproma,jrow,:)),                         &
                          CV_PRECNEW(_RI_XYZ__(1:kproma,jrow,:)),                         &
                          U_ENTR(_RI_XYZ__(1:kproma,jrow,:)),   U_DETR(_RI_XYZ__(1:kproma,jrow,:)),  &
                          D_ENTR(_RI_XYZ__(1:kproma,jrow,:)),   D_DETR(_RI_XYZ__(1:kproma,jrow,:)),  &
                          CAPE(1:kproma,jrow),                                 &
                          CV_LWC(_RI_XYZ__(1:kproma,jrow,:)),   CV_IWC(_RI_XYZ__(1:kproma,jrow,:)),  &
                          CV_RFORM(_RI_XYZ__(1:kproma,jrow,:)), CV_SFORM(_RI_XYZ__(1:kproma,jrow,:)) )
     

     do jk=1,nlev
       do jl=1,kproma
         tte_3d(_RI_XYZ__(jl,jrow,jk)) = &
              tte_3d(_RI_XYZ__(jl,jrow,jk))  + conv_tte(_RI_XYZ__(jl,jrow,jk))
         qte_3d(_RI_XYZ__(jl,jrow,jk))  = &
              qte_3d(_RI_XYZ__(jl,jrow,jk))  + conv_qte(_RI_XYZ__(jl,jrow,jk))
         vom_3d(_RI_XYZ__(jl,jrow,jk))  = &
              vom_3d(_RI_XYZ__(jl,jrow,jk)) + conv_ute(_RI_XYZ__(jl,jrow,jk))
         vol_3d(_RI_XYZ__(jl,jrow,jk))  = &
              vol_3d(_RI_XYZ__(jl,jrow,jk)) + conv_vte(_RI_XYZ__(jl,jrow,jk))
!!$         do jt=1,ntrac 
!!$           if (ztracconv(jt).eq.1) &
!!$             pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + tracten(jl,jk,jt)
!!$         enddo
       enddo
     enddo

     DO jk=1,nlev
       DO jl=1,kproma
         IF (ZTP1(JL,JK) < Tmelt) THEN
           CV_SNOWFLX(_RI_XYZ__(JL,jrow,jk)) = CV_PRECFLX(_RI_XYZ__(JL,jrow,jk))
           CV_PRECFLX(_RI_XYZ__(JL,jrow,jk)) = 0._dp
         ENDIF
       ENDDO
     ENDDO
     do jl=1,kproma
       rsfc_2d(jl,jrow) = cv_precflx(_RI_XYZ__(jl,jrow,nlev))
       ssfc_2d(jl,jrow) = cv_snowflx(_RI_XYZ__(jl,jrow,nlev))
       aprc(jl,jrow) = aprc(jl,jrow) &
            + delta_time * (rsfc_2d(jl,jrow) + ssfc_2d(jl,jrow))
       aprs(jl,jrow) = aprs(jl,jrow) + delta_time * ssfc_2d(jl,jrow)

       if (itype(jl) == 1) then
         ktype(jl) = 2
         DO jk=1, klev
           zpb = aphp1(jl,jk+1)
           ! ECHAM5 top layer ends at 0. Pa !!! adjust to 0.01 Pa
           IF (.NOT. ledith) THEN ! op_pj_20180725
              zpt = MAX(aphp1(jl,jk),0.01_dp)
           ELSE                  ! op_pj_20180725
! ka_sv_20160224+
              ! original cut-off too low if upper level is above 0.01_dp Pa
              zpt = MAX(aphp1(jl,jk),app1(jl,1)/100._dp)
! ka_sv_20160224-              
           ENDIF                 ! op_pj_20180725
           grheight(jl,jk) = (1000._dp * R_gas / (M_air * g)) &
                             * ztp1(jl,jk) * log(zpb/zpt)
         END DO
         
         do jk=klev,ictop(jl),-1
           cth(jl,jrow) = cth(jl,jrow) + grheight(jl,jk)
         enddo
         cbh = 0._dp
         do jk=klev,icbot(jl),-1 
           cbh = cbh + grheight(jl,jk)
         enddo
         if (cth(jl,jrow) - cbh > 2000._dp) ktype(jl) = 1
       endif
       
       !set default cloud bottom heigt = cloud top height
       if (ktype(jl) == 0) then
          icbot(jl) = ictop(jl)
       endif   

       conv_type(jl,jrow) = REAL(ktype(jl),dp)
       CONV_TOP(jl,jrow)  = REAL(ictop(jl),dp)
       CONV_BOT(jl,jrow)  = REAL(icbot(jl),dp)

!      values for NOx lightning
! op_mm_20131217 -> RI
       IF (ktype(jl) == 1) then
         cu_bot(jl,jrow)    = conv_bot(jl,jrow)
         cu_top(jl,jrow)    = conv_top(jl,jrow)
         cu_freeze(jl,jrow) = conv_top(jl,jrow)
         if (nint(cu_top(jl,jrow)) > 0) then
           do jk=nint(cu_top(jl,jrow)),nint(cu_bot(jl,jrow))
             if (ztp1(jl,jk).le.273.15_dp) cu_freeze(jl,jrow)=REAL(jk,dp)
           enddo
           do jk=nint(cu_top(jl,jrow)),nint(cu_bot(jl,jrow))
             cu_uvelo(_RI_XYZ__(jl,jrow,jk))=massfu(_RI_XYZ__(jl,jrow,jk))/zrhoa(jl,jk)
           enddo
         endif
       ENDIF
!    end of values for NOx lightning

     enddo
     

   CASE(8)       ! Donner scheme 2006/07
     c_cldfrac(_RI_XYZ__(:,jrow,:)) = 0._dp
     c_liqamt(_RI_XYZ__(:,jrow,:))  = 0._dp
     c_liqsize(_RI_XYZ__(:,jrow,:)) = 0._dp
     c_iceamt(_RI_XYZ__(:,jrow,:))  = 0._dp
     c_icesize(_RI_XYZ__(:,jrow,:)) = 0._dp
     m_cldfrac(_RI_XYZ__(:,jrow,:)) = 0._dp
     m_liqamt(_RI_XYZ__(:,jrow,:))  = 0._dp
     m_liqsize(_RI_XYZ__(:,jrow,:)) = 0._dp
     m_iceamt(_RI_XYZ__(:,jrow,:))  = 0._dp
     m_icesize(_RI_XYZ__(:,jrow,:)) = 0._dp
     conv_active(:,1)    = .TRUE.

     do jk = 1,klev
       do jl = 1, kproma
         temp(jl,1,jk)     = ztp1(jl,jk)
         spec_hum(jl,1,jk) = zqp1(jl,jk)
         zpress(jl,1,jk)   = app1(jl,jk)
         OMEGA(jl,1,jk)    = vervel(jl,jk)

         cell_cld_frac(jl,1,jk) = c_cldfrac(_RI_XYZ__(jl,jrow,jk))
         cell_liq_amt(jl,1,jk)  = c_liqamt(_RI_XYZ__(jl,jrow,jk))
         cell_liq_size(jl,1,jk) = c_liqsize(_RI_XYZ__(jl,jrow,jk))
         cell_ice_amt(jl,1,jk)  = c_iceamt(_RI_XYZ__(jl,jrow,jk))
         cell_ice_size(jl,1,jk) = c_icesize(_RI_XYZ__(jl,jrow,jk))
         meso_cld_frac(jl,1,jk) = m_cldfrac(_RI_XYZ__(jl,jrow,jk))
         meso_liq_amt(jl,1,jk)  = m_liqamt(_RI_XYZ__(jl,jrow,jk))
         meso_liq_size(jl,1,jk) = m_liqsize(_RI_XYZ__(jl,jrow,jk))
         meso_ice_amt(jl,1,jk)  = m_iceamt(_RI_XYZ__(jl,jrow,jk))
         meso_ice_size(jl,1,jk) = m_icesize(_RI_XYZ__(jl,jrow,jk))

         donner_humidity_area(jl,1,jk)  = d_humarea(_RI_XYZ__(jl,jrow,jk))
         donner_humidity_ratio(jl,1,jk) = d_humratio(_RI_XYZ__(jl,jrow,jk))


         qlin(jl,1,jk) = zxlp1(jl,jk)
         qiin(jl,1,jk) = zxip1(jl,jk)
         !!$qain(jl,1,jk) = aclc(jl,jk,jrow)   ! op_mm_20140327
         qain(jl,1,jk) = aclc(_RI_XYZ__(jl,jrow,jk))

         delta_ql(jl,1,jk) = conv_lte(_RI_XYZ__(jl,jrow,jk))
         delta_qi(jl,1,jk) = conv_ite(_RI_XYZ__(jl,jrow,jk))
         delta_qa(jl,1,jk) = conv_covte(_RI_XYZ__(jl,jrow,jk))
       enddo
     enddo
     do jk=1,nlev+1
       do jl=1,kproma
         zpressi(jl,1,jk)  = aphp1(jl,jk)
       enddo
     enddo
     do jl=1,kproma
       land(jl,1)          = slf(jl,jrow)
       sfc_sh_flux(jl,1)   = s_heatflux(jl,jrow)
       ! latent heat flux per vaporisation energy
       sfc_vapor_flux(jl,1)= l_heatflux(jl,jrow) / 2.5008e6_dp

       precip(jl,1)   = rsfc_2d(jl,jrow)
     enddo

     do jk=1,klev
       do jt=1,ntrac
         do jl=1,kproma
           tracers(jl,1,jk,jt) = zxtp1(jl,jk,jt)
         enddo
       enddo
     enddo
     do jt=1,ntrac
       do jl=1,kproma
         tr_flux(jl,1,jt) = 0._dp
       enddo
     enddo
     nsum(1:kproma,1) = 0
     CALL DONNER_CUMASTR(1, kproma, 1, 1, time_step_len,                &
                         temp, spec_hum, zpress, zpressi, omega,        &
                         land, sfc_sh_flux, sfc_vapor_flux,             &
                         tr_flux, tracers, TIME,                        &
                         cell_cld_frac,  &
                         cell_liq_amt, cell_liq_size, cell_ice_amt,   &
                         cell_ice_size, meso_cld_frac, meso_liq_amt, &
                         meso_liq_size, meso_ice_amt, meso_ice_size,  &
                         nsum, precip, delta_temp, delta_vapor, detf, &
                         uceml_inter, mtot, donner_humidity_area,    &
                         donner_humidity_ratio, qtrtnd, &
                        ! mz_ht_20070421+                        
                         p_pe, conv_active,             &
                        ! mz_ht_20070421-
                         qlin, qiin, qain,              &      ! optional
                         delta_ql, delta_qi, delta_qa)         ! optional

      do jl=1,kproma
        rsfc_2d(jl,jrow) = precip(jl,1)/time_step_len
        aprc(jl,jrow) = aprc(jl,jrow)+delta_time*rsfc_2d(jl,jrow)
        ! op_mm_20140109
         mlev=klev
        cv_precflx(_RI_XYZ__(jl,jrow,mlev)) = rsfc_2d(jl,jrow)
        IF (.NOT. CONV_ACTIVE(jl,1)) conv_type(jl,jrow) = 1._dp
      enddo
      do jk = 1,klev
        do jl = 1, kproma
          c_cldfrac(_RI_XYZ__(jl,jrow,jk)) = cell_cld_frac(jl,1,jk) 
          c_liqamt(_RI_XYZ__(jl,jrow,jk))  = cell_liq_amt(jl,1,jk)  
          c_liqsize(_RI_XYZ__(jl,jrow,jk)) = cell_liq_size(jl,1,jk) 
          c_iceamt(_RI_XYZ__(jl,jrow,jk))  = cell_ice_amt(jl,1,jk)  
          c_icesize(_RI_XYZ__(jl,jrow,jk)) = cell_ice_size(jl,1,jk) 
          m_cldfrac(_RI_XYZ__(jl,jrow,jk)) = meso_cld_frac(jl,1,jk) 
          m_liqamt(_RI_XYZ__(jl,jrow,jk))  = meso_liq_amt(jl,1,jk)  
          m_liqsize(_RI_XYZ__(jl,jrow,jk)) = meso_liq_size(jl,1,jk) 
          m_iceamt(_RI_XYZ__(jl,jrow,jk))  = meso_ice_amt(jl,1,jk)  
          m_icesize(_RI_XYZ__(jl,jrow,jk)) = meso_ice_size(jl,1,jk) 

          d_humarea(_RI_XYZ__(jl,jrow,jk))  = donner_humidity_area(jl,1,jk)
          d_humratio(_RI_XYZ__(jl,jrow,jk)) = donner_humidity_ratio(jl,1,jk)
          
          massfu(_RI_XYZ__(jl,jrow,jk))     = mtot(jl,1,jk)
          conv_lte(_RI_XYZ__(jl,jrow,jk))   = delta_ql(jl,1,jk) / time_step_len
          conv_ite(_RI_XYZ__(jl,jrow,jk))   = delta_qi(jl,1,jk) / time_step_len
          conv_covte(_RI_XYZ__(jl,jrow,jk)) = delta_qa(jl,1,jk) / time_step_len
          conv_qte(_RI_XYZ__(jl,jrow,jk))   = delta_vapor(jl,1,jk) / time_step_len
          conv_tte(_RI_XYZ__(jl,jrow,jk))   = delta_temp(jl,1,jk) / time_step_len
       enddo
     enddo
     
     do jk=1,klev
       do jl=1,kproma
         tte_3d(_RI_XYZ__(jl,jrow,jk))  = &
              tte_3d(_RI_XYZ__(jl,jrow,jk))  + conv_tte(_RI_XYZ__(jl,jrow,jk))
         qte_3d(_RI_XYZ__(jl,jrow,jk))  = &
              qte_3d(_RI_XYZ__(jl,jrow,jk))  + conv_qte(_RI_XYZ__(jl,jrow,jk))
       enddo
     enddo

    END SELECT

#ifdef MESSYTENDENCY
    !==============================================
    !Tendency budgets, start:
    !==============================================

    ! finish Tiedtke(1-3) schemes and ECMWF(4) scheme

    if(convect_param >= 1 .and. convect_param <= 4)   then
       ! tendency budget
       lo_tte(1:kproma,:)  = tte_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_tte(1:kproma,:)!calculating the
       lo_vom(1:kproma,:)  = vom_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_vom(1:kproma,:)!new (local)
       lo_vol(1:kproma,:)  = vol_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_vol(1:kproma,:)!tendency
       lo_qte(1:kproma,:)  = qte_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_qte(1:kproma,:)
       lo_xtte(1:kproma,:,:) = pxtte(1:kproma,:,:) -lo_xtte(1:kproma,:,:)
       
       tte_3d(_RI_XYZ__(1:kproma,jrow,:))  = &
            tte_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_tte(1:kproma,:)!subtract new
       vom_3d(_RI_XYZ__(1:kproma,jrow,:))  = &
            vom_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_vom(1:kproma,:)!(local) tendency
       vol_3d(_RI_XYZ__(1:kproma,jrow,:))  = &
            vol_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_vol(1:kproma,:)!from total tendency
       qte_3d(_RI_XYZ__(1:kproma,jrow,:))  = &
            qte_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_qte(1:kproma,:)
       pxtte(1:kproma,:,:)  = pxtte(1:kproma,:,:) -lo_xtte(1:kproma,:,:)

       call mtend_add_l (my_handle, mtend_id_t, px = lo_tte)!add new tendency
       call mtend_add_l (my_handle, mtend_id_u, px = lo_vom)!via "add" (tendency)
       call mtend_add_l (my_handle, mtend_id_v, px = lo_vol)!interface
       call mtend_add_l (my_handle, mtend_id_q, px = lo_qte)
       ! ub_ak_20190613+
       !if(ntrac .gt. 0) call mtend_add_l (my_handle, mtend_id_tracer, pxt = lo_xtte)
       if(ntrac .gt. 0) then
          DO jt = 1, ntrac
             call mtend_add_l (my_handle, jt, px = lo_xtte(:,:,jt))
          END DO
       end if
       ! ub_ak_20190613-
    endif

    ! additional for ECMWF(4) scheme
    if(convect_param .eq. 4)   then
       lo_xlte(1:kproma,:)  = &
            xlte_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_xlte(1:kproma,:) 
       lo_xite(1:kproma,:)  = &
            xite_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_xite(1:kproma,:)
       
       xlte_3d(_RI_XYZ__(1:kproma,jrow,:))  = &
            xlte_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_xlte(1:kproma,:) 
       xite_3d(_RI_XYZ__(1:kproma,jrow,:))  = &
            xite_3d(_RI_XYZ__(1:kproma,jrow,:)) - lo_xite(1:kproma,:)
     
       call mtend_add_l (my_handle, mtend_id_xl, px = lo_xlte)
       call mtend_add_l (my_handle, mtend_id_xi, px = lo_xite)

    end if

!----!-------------------!----------------------------------
    ! - Tracer changes and tendencies are not implemented for these schemes!
    ! - For scheme 8 the tendency-tool is not implemented at all, for it
    !   is not working in this version!

    ! Zhang(5), Bechtold(6), Emanuel(7), Donner(8) scheme
    if(convect_param .ge. 5 .and. convect_param .le. 8)   then

       do jk=1,klev
          do jl=1,kproma
             tte_3d(_RI_XYZ__(jl,jrow,jk)) =  &
                  tte_3d(_RI_XYZ__(jl,jrow,jk)) &
                  - conv_tte(_RI_XYZ__(jl,jrow,jk))
             qte_3d(_RI_XYZ__(jl,jrow,jk))  = &
                 qte_3d(_RI_XYZ__(jl,jrow,jk)) - conv_qte(_RI_XYZ__(jl,jrow,jk))
             if (convect_param .eq. 6 .or. convect_param .eq. 7)   then
                do jt=1,ntrac 
                   if (ztracconv(jt).eq.1) & 
!!                        pxtte(jl,jk,jt)  = pxtte(jl,jk,jt) -tracten(jl,jk,jt)
                        pxtte(_RI_X_ZN_(jl,jk,jt))  = pxtte(_RI_X_ZN_(jl,jk,jt)) -tracten(jl,jk,jt)
                enddo
             endif

 !!$           if (convect_param .eq. 5) then
 !!$              do jt=1,ntrac
 !!$                 lo_xtte(jl,jk,jt) = pxtte(jl,jk,jt) -lo_xtte(jl,jk,jt)
 !!$                 pxtte(jl,jk,jt)  = pxtte(jl,jk,jt) -lo_xtte(jl,jk,jt)
 !!$              enddo
 !!$           endif

          enddo
       enddo

       call mtend_add_l (my_handle, mtend_id_t, px = conv_tte(_RI_XYZ__(1:kproma,jrow,:)))
       call mtend_add_l (my_handle, mtend_id_q, px = conv_qte(_RI_XYZ__(1:kproma,jrow,:)))
      
 !!$     if (convect_param .eq. 6 .or. convect_param .eq. 7)   then
 !!$        call mtend_add_l (my_handle, mtend_id_tracer, pxt = tracten)
 !!$     endif

 !!$     if (convect_param .eq. 5) then
 !!$        call mtend_add_l (my_handle, mtend_id_tracer, pxt = lo_xtte)
 !!$    endif

    endif
!-------------------
    if(convect_param .ge. 6 .and. convect_param .le. 7) then

       do jk=1,klev
          do jl=1,kproma
             vom_3d(_RI_XYZ__(jl,jrow,jk))  = &
                 vom_3d(_RI_XYZ__(jl,jrow,jk)) - conv_ute(_RI_XYZ__(jl,jrow,jk))
             vol_3d(_RI_XYZ__(jl,jrow,jk))  = &
                 vol_3d(_RI_XYZ__(jl,jrow,jk)) - conv_vte(_RI_XYZ__(jl,jrow,jk))
          enddo
       enddo
       call mtend_add_l (my_handle, mtend_id_u, px = conv_ute(_RI_XYZ__(1:kproma,jrow,:)))
       call mtend_add_l (my_handle, mtend_id_v, px = conv_vte(_RI_XYZ__(1:kproma,jrow,:)))
    endif
!-------------------
    if(convect_param .eq. 5 .or. convect_param .eq. 6)then! .or. convect_param .eq. 8) then

       do jk=1,klev
          do jl=1,kproma
             xlte_3d(_RI_XYZ__(jl,jrow,jk))  = &
                xlte_3d(_RI_XYZ__(jl,jrow,jk)) - conv_lte(_RI_XYZ__(jl,jrow,jk))
             xite_3d(_RI_XYZ__(jl,jrow,jk))  = &
                xite_3d(_RI_XYZ__(jl,jrow,jk)) - conv_ite(_RI_XYZ__(jl,jrow,jk))
          enddo
       enddo

       call mtend_add_l (my_handle, mtend_id_xl, px = conv_lte(_RI_XYZ__(1:kproma,jrow,:)))
       call mtend_add_l (my_handle, mtend_id_xi, px = conv_ite(_RI_XYZ__(1:kproma,jrow,:)))

    endif

!==============================================
!Tendency budgets, end:
!==============================================
#endif

! op_mm_20140110 finish -> error_bi
  IF (lookupoverflow) CALL error_bi('convection - lookuperror', substr)
!
!
! ------------------------------------------------------------------

  if (convect_param.le.3) then

!
!*     3.     PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
!             -------- -------- -- ---------- ----- -----
!
300 CONTINUE

    ! mim_sb_20091009+
     IF (l_lgmc_diag) then
        do jl=1,kproma
           aprsc(jl,jrow) = aprsc(jl,jrow) &
                + delta_time * (rsfc_2d(jl,jrow) + ssfc_2d(jl,jrow))
           aprss(jl,jrow) = aprss(jl,jrow) + delta_time * ssfc_2d(jl,jrow)
        enddo
     endif
    ! mim_sb_20091009-
!
    ilevmin=klev-4
!
    DO 301 jl=1,kproma
      itopec2(jl) = klevp1
301 END DO
!
    DO 303 jk=1,ilevmin
      DO 302 jl=1,kproma
        IF(jlab(jl,jk).EQ.2 .AND. itopec2(jl).EQ.klevp1) THEN
          itopec2(jl)=jk
        END IF
302   END DO
303 END DO
!
    ztopmax(1:kproma) = ptopmax(1:kproma)
    kmin=1
    kmax=klev
    pilab(_RI_XYZ__(1:kproma,jrow,kmin:kmax)) = REAL(jlab(1:kproma,1:klev), dp)
    DO 304 jl=1,kproma
      IF(itopec2(jl).EQ.1) THEN
        ptopmax(jl)=app1(jl,1)
      ELSE IF(itopec2(jl).NE.klevp1) THEN
        ptopmax(jl)=aphp1(jl,itopec2(jl))
      ELSE
        ptopmax(jl)=99999.
      END IF
      ptopmax(jl)=MIN(ptopmax(jl),ztopmax(jl))
304 END DO

!  added for convective cloud top height in m (convect stream)

    DO jk=1, klev
      DO jl=1, kproma
        zpb = aphp1(jl,jk+1)
        ! ECHAM5 top layer ends at 0. Pa !!! adjust to 0.01 Pa
        IF (.NOT. ledith) THEN ! op_pj_20180725
           zpt = MAX(aphp1(jl,jk),0.01_dp)
        ELSE                   ! op_pj_20180725
! ka_sv_20160224+
           ! original cut-off too low if upper level is above 0.01_dp Pa
           zpt = MAX(aphp1(jl,jk),app1(jl,1)/100._dp)
! ka_sv_20160224-
        ENDIF                  ! op_pj_20180725
        grheight(jl,jk) = (1000._dp * R_gas / (M_air * g)) &
                          * ztp1(jl,jk) * log(zpb/zpt)
      END DO
    END DO
    do jl=1,kproma
      do jk=klev,itopec2(jl),-1
        cth(jl,jrow) = cth(jl,jrow) + grheight(jl,jk)
      enddo
    enddo

    

  ! WRITE OUT CONVECTION TYPE, UPDRAFT AND DOWNDRAFT MASS FLUX [kg/m^2/s] 
  ! TO CONVECT CHANNEL
    conv_type(1:kproma,jrow) = REAL(ktype(1:kproma),dp)    
    conv_bot(1:kproma,jrow)  = REAL(icbot(1:kproma),dp)
    conv_top(1:kproma, jrow) = REAL(ictop(1:kproma),dp)
! calculate estimate of convective cloud cover not assuming the 5% from Tiedtke
! but estimating from the strength of the updraft
  endif
  DEALLOCATE(jlab)

! op_mm_20131217 -> RI
  do jk=1,klev
     do jl=1,kproma
       if (massfu(_RI_XYZ__(jl,jrow,jk)) .gt.0._dp) then
         call calc_conv_cover (cv_cover(_RI_XYZ__(jl,jrow,jk)), massfu(_RI_XYZ__(jl,jrow,jk)) &
           , zrhoa(jl,jk), ztmst)
       endif
       ! mz_jd_20161011+
       call calc_conv_cover_sikma(cv_cover_sikma(_RI_XYZ__(jl,jrow,jk))         &
            , qp1(_RI_XYZ__(jl,jrow,jk)), qsat(_RI_XYZ__(jl,jrow,jk)), qstdev(_RI_XYZ__(jl,jrow,jk))  &
            , Q2(_RI_XYZ__(jl,jrow,jk)), qskew(_RI_XYZ__(jl,jrow,jk))                      &
            , xvar(_RI_XYZ__(jl,jrow,jk)), xskew(_RI_XYZ__(jl,jrow,jk)), aclc(_RI_XYZ__(jl,jrow,jk)) )
       ! mz_jd_20161011-
     enddo
  enddo

  WHERE (NINT(conv_type(1:kproma,jrow) - 1._dp) == 0)
     counter_deep(1:kproma,jrow) = 1._dp
     cth_deep(1:kproma,jrow) = cth(1:kproma,jrow)
  END WHERE
  WHERE (NINT(conv_type(1:kproma,jrow) - 2._dp) == 0)
     counter_shal(1:kproma,jrow) = 1._dp
! op_sb_20150115+
!!$  cth_midl(1:kproma,jrow) = cth(1:kproma,jrow)
     cth_shal(1:kproma,jrow) = cth(1:kproma,jrow)
! op_sb_20150115-
  END WHERE
  WHERE (NINT(conv_type(1:kproma,jrow) - 3._dp) == 0)
     counter_midl(1:kproma,jrow) = 1._dp
     cth_midl(1:kproma,jrow) = cth(1:kproma,jrow)
  END WHERE

! op_mm_20131217
  DO jl=1,kproma
    IF (counter_midl(jl,jrow) == 1._dp) &
      massfu_midl(_RI_XYZ__(jl,jrow,:)) = massfu(_RI_XYZ__(jl,jrow,:))
    IF (counter_deep(jl,jrow) == 1._dp) &
      massfu_deep(_RI_XYZ__(jl,jrow,:)) = massfu(_RI_XYZ__(jl,jrow,:))
    IF (counter_shal(jl,jrow) == 1._dp) &
      massfu_shal(_RI_XYZ__(jl,jrow,:)) = massfu(_RI_XYZ__(jl,jrow,:))  
  END DO

!
!
!---------------------------------------------------------------------
!

! op_mm_20140110
#ifdef ECHAM5
  IF (ltdiag) THEN
! store CUCALL increment (massflux)
    pdiga5(1:kproma,:,jrow)  = pdiga5(1:kproma,:,jrow)  + vom_3d(1:kproma,:,jrow)
    pdiga10(1:kproma,:,jrow) = pdiga10(1:kproma,:,jrow) + vol_3d(1:kproma,:,jrow)
    pdiga18(1:kproma,:,jrow) = pdiga18(1:kproma,:,jrow) + tte_3d(1:kproma,:,jrow)
  ENDIF
#endif


  ! op_mm_20140123+
  ! update cosmo variables
#ifdef COSMO
  ! cloud cover due to convection  (3D)    
   clc_con(1:kproma,jrow,1:nlev)=cv_cover(_RI_XYZ__(1:kproma,jrow,:))       
   ! cloud liquid water due to convection  (3D)
   clw_con(1:kproma,jrow,1:nlev)=cv_cldwater(_RI_XYZ__(1:kproma,jrow,:))    
   ! precipitation rate of rain, convective  (kg/m2*s) (2D)
   prr_con(1:kproma,jrow)=rsfc_2d(1:kproma,jrow)
   ! precipitation rate of snow, convective  (kg/m2*s) (2D)
   prs_con(1:kproma,jrow)=ssfc_2d(1:kproma,jrow)         
   ! precipitation rate, no evaporat., convective  (kg/m2*s) (2D)
   prne_con(1:kproma, jrow) =zrain(:)
   !level index of convective cloud base (2D)
   bas_con(1:kproma,jrow) = conv_bot(1:kproma,jrow)     
   !level index of convective cloud top (2D)
   top_con(1:kproma,jrow) = conv_top(1:kproma,jrow)     
   ! TKE-tendency due to convective buoyancy       ( m2/s3 ) (3D)
   tket_conv(1:kproma,jrow, 1:nlev) = tketconv(_RI_XYZ__(1:kproma,jrow,:))  
   ! surface flux of water vapour (2D) qvsflx=qflux 
   ! qvsflx    = 0.0_dp 
   ! cloud base massflux     (kg/m2*s) (2D) 
   mflx_con(1:kproma,jrow)  = base_f1(1:kproma,jrow)+base_f2(1:kproma,jrow)&
        +base_f3(1:kproma,jrow)+base_f4(1:kproma,jrow) 
   cape_con(1:kproma,jrow)  = cape(1:kproma,jrow)
   qcvg_con  = 0.0_dp  ! moisture convergence for Kuo-type closure (    1/s) (2D)
   tke_con   = 0.0_dp  ! convective turbulent energy               (   J/kg) (2D)
   ! maximum convective gust at 10m                ( m/s ) (2D)
   vgust_con(1:kproma, jrow) = vgustcon(1:kproma,jrow)      
#endif

#ifdef CESM1
  ! mz_ab_20140410+
  ! for CESM, averaging is done through CHANNEL 
  ! instead of the echam laccu stream objects routines, 
  ! so correct those channel objects:
  ! convect_convec is before cloud_convec
  aprc(1:kproma,jrow) = rsfc_2d(1:kproma,jrow)+ssfc_2d(1:kproma,jrow)
  aprs(1:kproma,jrow) = ssfc_2d(1:kproma,jrow)
  ! mz_ab_20140410-
  ! mz_ab_20140825+
  aprflux(1:kproma,jrow)  = rsfc_2d(1:kproma,jrow)+ssfc_2d(1:kproma,jrow)
  ! mz_ab_20140825-
#endif

  RETURN


END SUBROUTINE convect_CONVEC

!==============================================================================

! op_pj_20091009+
! ========================================================================
  SUBROUTINE convect_read_nml_cpl(status, iou)
   
    ! Convect MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Aug 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ l_lgmc_diag, pbl_height ! um_ak_20140502 pbl_height added
 
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='convect_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! um_ak_20140502
    ! initialise pbl_height
    pbl_height%CHA = ' '
    pbl_height%OBJ = ' '
    ! um_ak_20140502-

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE convect_read_nml_cpl
! ========================================================================
! op_pj_20091009-

END MODULE MESSY_CONVECT_SI
