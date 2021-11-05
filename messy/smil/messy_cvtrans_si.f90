#include "messy_main_ppd_bi.inc"

MODULE MESSY_cvtrans_si

!     submodule that treats the convective tracer transport in the MESSY
!     structure    (mz_ht_20031215)

!     core module by Mark Lawrence (December 2003)

!     interface expanded for simultaneous usage with ECHAM5 and COSMO
!     by A. Kerkweg in July 2009

!     alternative algorithm by H.G.Ouwersloot

  ! MESSy/BMIL
  USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, ntrac_gp, ti_gp, &
                                      LGTRSTR, ntrac_lg, ti_lg, &
                                      t_trinfo_tp, ON, OFF
  ! MESSy/SMCL
  USE messy_cvtrans
  USE messy_main_constants_mem, ONLY: dp
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle  ,&
                                      mtend_register      ,&
                                      mtend_add_l       ,&
                                      mtend_get_start_l ,&
                                      mtend_id_tracer
#endif

  IMPLICIT NONE

  SAVE

  CHARACTER(LEN=80)  ::                          &
! char for convective precipitation
                        raincv(2)               ,&  
!      for convective snow
                        snowcv(2)               ,&  
!      for convective cover estimate
                        covcv(2)                ,&  
!      for grid volume, mass
                        vol(2), mass(2)         ,&  
!      for pressure
                        press(2), pressi(2)     ,&  
!      for convective mass fluxes
                        dmass(2), umass(2)      ,&  
!      for convective detrainment fluxes
                        detru(2), detrd(2)      ,&  
!      for convective entrainment fluxes
                        entru(2), entrd(2)      ,&  
!      for convective top level index
                        c_top(2)                ,&
                        pblh(2)                    ! op_pj_20150428
  
! um_ak_20090721+: moved from SMCL ...
! Pointers for own channel elements  
  REAL(dp) , PUBLIC, POINTER :: cumassf(:,:,:)     => NULL()  ! convective upward mass flux after closure
  REAL(dp) , PUBLIC, POINTER :: cdmassf(:,:,:)     => NULL()  ! convective downward mass flux after closure
  REAL(dp) , PUBLIC, POINTER :: cuentr(:,:,:)      => NULL()  ! entrainment of the upward flux after closure
  REAL(dp) , PUBLIC, POINTER :: cudetr(:,:,:)      => NULL()  ! detrainment of the upward flux after closure
  REAL(dp) , PUBLIC, POINTER :: cdentr(:,:,:)      => NULL()  ! entrainment of the downward flux after closure
  REAL(dp) , PUBLIC, POINTER :: cddetr(:,:,:)      => NULL()  ! detrainment of the downward flux after closure

  REAL(dp) , PUBLIC, POINTER :: updr_velo(:,:,:)   => NULL()  ! corrected updraft velocity after closure

! op_pj_20110218+ not used
!!$  REAL(dp) , PUBLIC, POINTER :: kbot(:,:)          => NULL()  ! lowest level of convection in cvtrans
!!$  REAL(dp) , PUBLIC, POINTER :: ktop(:,:)          => NULL()  ! top level of convection in cvtrans
!!$  REAL(dp) , PUBLIC, POINTER :: kraintop(:,:)      => NULL()  ! top level of convective rain in cvtrans
! op_pj_20110218-
  
  REAL(dp) , PUBLIC, POINTER :: trac_field(:,:,:,:) => NULL()  ! tracer field for coupling to scav
! um_ak_20090721-

  REAL(dp) , PUBLIC, POINTER :: trma(:,:,:,:) => NULL()  ! transport matrix

  
! Pointers for coupling channel elements  
!                for grid mass
  REAL (dp), POINTER :: grmass(:,:,:)      => NULL()    
!                for grid volume
  REAL (dp), POINTER :: grvol(:,:,:)       => NULL()    
!                for pressure
  REAL (dp), POINTER :: press_3d(:,:,:)    => NULL()
  REAL (dp), POINTER :: pressi_3d(:,:,:)   => NULL()    
!                for convective precipitation
  REAL (dp), POINTER :: precflx_cv(:,:,:)  => NULL()
  REAL (dp), POINTER :: snowflx_cv(:,:,:)  => NULL()    
!                for convective cloud cover estimate
  REAL (dp), POINTER :: ccover(:,:,:)      => NULL()    
! convective upward mass flux 
  REAL(dp) , PUBLIC, POINTER :: umassf(:,:,:)      => NULL()  
! convective downward mass flux
  REAL(dp) , PUBLIC, POINTER :: dmassf(:,:,:)      => NULL()  

! entrainment of the upward flux
  REAL(dp) , PUBLIC, POINTER :: uentr(:,:,:)       => NULL()  
! detrainment of the upward flux
  REAL(dp) , PUBLIC, POINTER :: udetr(:,:,:)       => NULL()   

! entrainment of the downward flux
  REAL(dp) , PUBLIC, POINTER :: dentr(:,:,:)       => NULL()  

! top level of convection (all kinds) original
  REAL(dp) , PUBLIC, POINTER :: topconv(:,:)       => NULL()  


  REAL (dp), POINTER :: column                     => NULL()  
  INTEGER  :: nsets = 0   ! number of tracer sets
  INTEGER  :: ntrac
  LOGICAL  :: l_lg, l_gp
! channel object for ATTILA Number of Boxes per grid cell
  REAL (dp), PUBLIC, POINTER :: NCB(:,:,:)          => NULL()

  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti => NULL()

#ifdef MESSYTENDENCY
  ! variable for tendency budget
  INTEGER                         :: my_handle
  LOGICAL                         :: lfirst = .TRUE.
#endif

! mz_ho_20140826+
  REAL(dp), POINTER   :: pblh_2d(:,:) => NULL() ! op_pj_20150427
  LOGICAL             :: l_pblh = .FALSE.
! mz_ho_20140826-

  INTEGER             :: ilev_up
  
  INTRINSIC :: EXP, MIN, MAX, INT

! =============================================================================

CONTAINS
  SUBROUTINE  cvtrans_initialize

    ! cvtrans MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF XTSURF SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! 
    ! Author: Patrick Joeckel, MPICH, Feb 2002
    ! Modified, H. Tost, MPICH, 31-10-2003

    ! MESSy/BMIL
    USE messy_main_mpi_bi,          ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,      ONLY: error_bi, info_bi
    USE messy_main_tools,           ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    INTEGER     :: iou     ! I/O unit
    INTEGER     :: status  ! index
    CHARACTER(LEN=*), PARAMETER :: substr='cvtrans_initialize'
    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CVTRANS CORE ROUTINE:
       CALL init_cvtrans(iou, status)
       IF (status /= 0) CALL error_bi('read error in CTRL namelist', substr)
    END IF
    CALL p_bcast(bulktrans, p_io)
    CALL p_bcast(segtrans, p_io)
    CALL p_bcast(scav_trans, p_io)
    CALL p_bcast(lcvt_gp, p_io)
    IF (lcvt_gp) CALL info_bi('CVTRANS for gridpoint tracers active!', substr)
    CALL p_bcast(lcvt_lg, p_io)
    IF (lcvt_lg) CALL info_bi('CVTRANS for lagrangian tracers active!', substr)
    ! mz_ho_20140826+
    CALL p_bcast(lcvt_sc, p_io)
    CALL p_bcast(trans_fac, p_io)
    CALL p_bcast(hlimit, p_io)
    CALL p_bcast(lmeanconc, p_io)
    CALL p_bcast(lintsteps, p_io)
    CALL p_bcast(maxfrac,   p_io)
    ! mz_ho_20140826-

    call p_bcast(l_calc_trma, p_io)
  
    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL cvtrans_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('read error in CPL namelist',substr)
    END IF
   
    CALL p_bcast(umass(1), p_io)
    CALL p_bcast(umass(2), p_io)
    CALL p_bcast(dmass(1), p_io)  
    CALL p_bcast(dmass(2), p_io)  
    CALL p_bcast(entru(1), p_io)  
    CALL p_bcast(entru(2), p_io)  
    CALL p_bcast(detru(1), p_io) 
    CALL p_bcast(detru(2), p_io) 
    CALL p_bcast(entrd(1), p_io)    
    CALL p_bcast(entrd(2), p_io)    
    CALL p_bcast(detrd(1), p_io)  
    CALL p_bcast(detrd(2), p_io)  
    CALL p_bcast(raincv(1), p_io)  
    CALL p_bcast(raincv(2), p_io)  
    CALL p_bcast(snowcv(1), p_io)  
    CALL p_bcast(snowcv(2), p_io)  
    CALL p_bcast(covcv(1), p_io)  
    CALL p_bcast(covcv(2), p_io)    
    CALL p_bcast(mass(1), p_io)  
    CALL p_bcast(mass(2), p_io)  
    CALL p_bcast(vol(1), p_io)  
    CALL p_bcast(vol(2), p_io)  
    CALL p_bcast(press(1), p_io)  
    CALL p_bcast(press(2), p_io)  
    CALL p_bcast(pressi(1), p_io)  
    CALL p_bcast(pressi(2), p_io)  
    call p_bcast(c_top(1), p_io)
    call p_bcast(c_top(2), p_io)
    call p_bcast(pblh(1), p_io)     ! op_pj_20150428 (! mz_ho_20140826)
    call p_bcast(pblh(2), p_io)     ! op_pj_20150428 (! mz_ho_20140826)
  
    ! ub_ak_20181026+
    ! moved to init_memory
!!$#ifdef MESSYTENDENCY
!!$    IF (lfirst) my_handle = mtend_get_handle(modstr)
!!$    lfirst = .FALSE.
!!$    CALL mtend_register (my_handle, mtend_id_tracer)
!!$#endif
    ! ub_ak_20181026-

  END SUBROUTINE cvtrans_initialize
!===============================================================================

  SUBROUTINE cvtrans_init_memory

    ! MESSy/BMIL
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, nproma, ngpblks 
    USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi  &
                                   , info_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi, ONLY: GP_3D_MID &
                                   , GP_2D_HORIZONTAL, SCALAR &
                                   , DC_GP, DIMID_LON, DIMID_LEV, DIMID_LAT &
                                   , gp_nseg, gp_start, gp_cnt & ! mz_pj_20061112
                                   , gp_meml, gp_memu  ! mz_pj_20061112
    ! MESSy
    USE messy_main_channel,    ONLY: new_channel, new_channel_object &
                                   , new_attribute
    USE messy_main_channel_dimensions,  ONLY: new_dimension
    USE messy_main_channel_repr, ONLY: new_representation, AUTO &
                                     , set_representation_decomp &  ! mz_pj_20061112
                                     , IRANK, PIOTYPE_COL &
                                     , repr_def_axes

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr='cvtrans_init_memory'
    INTEGER                     :: status, js
    INTEGER                     :: DIMID_NTRAC, DIMID_LEV2
    INTEGER                     :: GP_4D_NTRAC, GP_4D_LEVLEV
    ! mz_pj_20061112+
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    ! mz_pj_20061112-

    ! ub_ak_20181026+
    ! moved from initialize
#ifdef MESSYTENDENCY
    IF (lfirst) my_handle = mtend_get_handle(modstr)
    lfirst = .FALSE.
    CALL mtend_register (my_handle, mtend_id_tracer)
#endif
    ! ub_ak_20181026-
    CALL start_message_bi(modstr,'MEMORY INITIALIZATION', substr)

    ! ADDITIONAL DIMENSION
! op_pj_20120130: nrac -> ntrac_gp
    CALL new_dimension(status, DIMID_NTRAC, 'ntrac', ntrac_gp)
    CALL channel_halt(substr, status)

    ! NEW REPRESENTATIONS
    CALL new_representation(status, GP_4D_NTRAC, 'GP_4D_NTRAC'   &
         , rank = 4, link = 'xxxx', dctype = DC_GP               &
         , dimension_ids = (/ &
         _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_NTRAC) /) &
         , ldimlen       = (/ &
         _RI_XYZN_(nproma, ngpblks, AUTO, AUTO) /)   &
         , output_order  = (/ _IN_XYZN_, _IX_XYZN_   &        ! E: 3,1,4,2
         , _IY_XYZN_, _IZ_XYZN_ /)       & ! C: 3,1,2,4
         , axis = repr_def_axes(_RI_XYZN_('X','Y','Z','N')) &
         )
    CALL channel_halt(substr, status)

    ! mz_pj_20061112+
    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))
    
    start(:,:) = gp_start(:,:)
    cnt(:,:) = gp_cnt(:,:)
    meml(:,:) = gp_meml(:,:)
    memu(:,:) = gp_memu(:,:)
    
    cnt(:,_IN_XYZN_) = ntrac
    memu(:,_IN_XYZN_) = ntrac
    
    CALL set_representation_decomp(status, GP_4D_NTRAC &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)
    ! mz_pj_20061112-

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    
    ! cvtrans channel now only contains the corrected values after the closure,
    ! the original values can be found in the convect channel

    CALL new_channel_object(status, modstr, 'cumassf', p3=cumassf)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cumassf', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cdmassf', p3=cdmassf)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdmassf', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cuentr', p3=cuentr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cuentr', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cudetr', p3=cudetr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cudetr', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cdentr', p3=cdentr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdentr', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'cddetr', p3=cddetr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cddetr', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'trac_field' &
         , p4=trac_field, reprid=GP_4D_NTRAC)

! op_pj_20110218+ not used
!!$    CALL new_channel_object(status, modstr, 'kbot' &
!!$         , p2=kbot, reprid=GP_2D_HORIZONTAL)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'kbot', 'units', c='levels')
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_channel_object(status, modstr, 'ktop' &
!!$         , p2=ktop, reprid=GP_2D_HORIZONTAL)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'ktop', 'units', c='levels')
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_channel_object(status, modstr, 'kraintop' &
!!$         , p2=kraintop, reprid=GP_2D_HORIZONTAL)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'kraintop', 'units', c='levels')
!!$    CALL channel_halt(substr, status)
! op_pj_20110218-

    CALL new_channel_object(status, modstr, 'column', p0=column, reprid=SCALAR)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'column', 'units', c='(0-1)')
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr,'MEMORY INITIALIZATION', substr)

    if (lcvt_gp .and. ntrac_gp == 0) then
      lcvt_gp = .FALSE.
      CALL info_bi('Gridpoint CVTRANS switched off - no tracers available' &
           , substr)
    endif
    if (lcvt_lg .and. ntrac_lg == 0) then
      lcvt_lg = .FALSE.
      CALL info_bi('Lagrangian CVTRANS switched off - no tracers available' &
           , substr)
    endif

    nsets = 0
    if (lcvt_gp) nsets = nsets + 1
    if (lcvt_lg) nsets = nsets + 1
    ALLOCATE(tr_attr(nsets))
    do js=1,nsets
      SELECT CASE (nsets)
      CASE(1)
        if (lcvt_gp) &
          ntrac = ntrac_gp
        if (lcvt_lg) &
          ntrac = ntrac_lg
      CASE(2)
        if (js == 1)  ntrac = ntrac_gp
        if (js == 2)  ntrac = ntrac_lg
      END SELECT

      ALLOCATE(tr_attr(js)%amount(ntrac))
    enddo

    IF (l_calc_trma) THEN
      CALL new_dimension(status, DIMID_LEV2, 'nlev2', nlev)
      CALL channel_halt(substr, status)

      ! NEW REPRESENTATIONS
      CALL new_representation(status, GP_4D_LEVLEV, 'GP_4D_LEVLEV'   &
         , rank = 4, link = 'xxxx', dctype = DC_GP               &
         , dimension_ids = (/ &
           _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_LEV2) /) &
         , ldimlen       = (/ &
            _RI_XYZN_(nproma, ngpblks, AUTO, AUTO) /)   &
         , output_order  = (/ _IN_XYZN_, _IX_XYZN_   &        ! E: 3,1,4,2
                            , _IY_XYZN_, _IZ_XYZN_ /)       & ! C: 3,1,2,4
         , axis = repr_def_axes(_RI_XYZN_('X','Y','Z','N')) &
          )
      CALL channel_halt(substr, status)

      nseg = gp_nseg
      ALLOCATE(start(nseg,IRANK))
      ALLOCATE(cnt(nseg,IRANK))
      ALLOCATE(meml(nseg,IRANK))
      ALLOCATE(memu(nseg,IRANK))
      
      start(:,:) = gp_start(:,:)
      cnt(:,:) = gp_cnt(:,:)
      meml(:,:) = gp_meml(:,:)
      memu(:,:) = gp_memu(:,:)
      
      cnt(:,_IN_XYZN_) = nlev
      memu(:,_IN_XYZN_) = nlev
      
      CALL set_representation_decomp(status, GP_4D_LEVLEV &
        , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
      CALL channel_halt(substr, status)
    
      DEALLOCATE(start) ; NULLIFY(start)
      DEALLOCATE(cnt)   ; NULLIFY(cnt)
      DEALLOCATE(meml)  ; NULLIFY(meml)
      DEALLOCATE(memu)  ; NULLIFY(memu)

      CALL new_channel(status, modstr//'_trma', reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_channel_object(status, modstr//'_trma', 'trma' &
        , p4=trma, reprid=GP_4D_LEVLEV)
      CALL channel_halt(substr, status)
      trma = 0._dp
           
    ENDIF
    
    column = 0._dp
    if (scav_trans.ge.4) column=1._dp

  END SUBROUTINE cvtrans_init_memory

! ==============================================================================

  SUBROUTINE cvtrans_read_nml_cpl(status, iou)

    ! cvtrans MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: H. Tost, MPICH, March 2004

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ umass, dmass, entru, detru,          &
                   entrd, detrd,                        &
                   raincv, snowcv, covcv,               &
                   mass, vol, press, pressi, c_top,     &
                   pblh ! op_pj_20150428 (! mz_ho_20140826)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='cvtrans_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
 
    raincv(:) = ''
    snowcv(:) = ''
    covcv(:)  = ''
    umass(:)  = ''
    dmass(:)  = ''
    entru(:)  = ''
    entrd(:)  = ''
    detru(:)  = ''    
    detrd(:)  = ''
    c_top(:)  = ''
    mass(:)    = ''
    vol(:)     = ''
    press(:)   = ''
    pressi(:)  = ''
    pblh(:)    = '' ! op_pj_20150428 (! mz_ho_20140826)

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
     
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE cvtrans_read_nml_cpl
! ==============================================================================

! coupling control

  SUBROUTINE cvtrans_init_coupling

    ! MESSY/BMIL
    USE messy_main_mpi_bi,          ONLY: p_parallel_io
    USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi &
                                        , info_bi, error_bi &
                                        , warning_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info

    USE messy_main_grid_def_mem_bi,  ONLY: nlev, nlevp1
#if defined(ECHAM5) || defined(CESM1)           
    USE messy_main_grid_def_mem_bi,  ONLY:  nvclev, vct
#endif
#if defined(COSMO) || defined(MESSYDWARF)
    USE messy_main_grid_def_bi,      ONLY: h_a => hyai, h_b=>hybi
#endif    
    
    IMPLICIT NONE

    INTEGER :: status, jt, js, jk
    REAL(dp) :: hypi(nlevp1), sfpress ! um_ak_20130724 sfpress added
#if defined(ECHAM5) || defined(CESM1)
    REAL(dp) :: h_a(nvclev), h_b(nvclev) ! um_ak_20130724, sfpress
#endif    
    CHARACTER(LEN=*), PARAMETER::substr='cvtrans_init_cpl'
    INTRINSIC :: TRIM

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)

    IF (p_parallel_io) &
         WRITE(*,*) 'Checking for convective parameters ...'
    
    CALL get_channel_object(status, TRIM(raincv(1)), TRIM(raincv(2)) &
         , p3=precflx_cv)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(raincv(1))//' - '//&
         &TRIM(raincv(2))//' not found', substr)

    CALL get_channel_object(status, TRIM(snowcv(1)), TRIM(snowcv(2)) &
         , p3=snowflx_cv)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(snowcv(1))//' - '//&
         &TRIM(snowcv(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(umass(1)), TRIM(umass(2)) &
         , p3=umassf)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(umass(1))//' - '//&
         &TRIM(umass(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(dmass(1)), TRIM(dmass(2)) &
         , p3=dmassf)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(dmass(1))//' - '//&
         &TRIM(dmass(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(entru(1)), TRIM(entru(2)) &
         , p3=uentr)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(entru(1))//' - '//&
         &TRIM(entru(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(entrd(1)), TRIM(entrd(2)) &
         , p3=dentr)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(entrd(1))//' - '//&
         &TRIM(entrd(2))//' not found', substr)

    CALL get_channel_object(status, TRIM(detru(1)), TRIM(detru(2)) &
         , p3=udetr)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(detru(1))//' - '//&
         &TRIM(detru(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(covcv(1)), TRIM(covcv(2)) &
         , p3=ccover)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(covcv(1))//' - '//&
         &TRIM(covcv(2))//' not found', substr)

    IF (p_parallel_io) &
         WRITE(*,*) 'Checking for grid box parameters ...'
    CALL get_channel_object(status, TRIM(mass(1)), TRIM(mass(2)) &
         , p3=grmass)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(mass(1))//' - '//&
         &TRIM(mass(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(vol(1)), TRIM(vol(2)) &
         , p3=grvol)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(vol(1))//' - '//&
         &TRIM(vol(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(press(1)), TRIM(press(2)) &
         , p3=press_3d)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(press(1))//' - '//&
         &TRIM(press(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(pressi(1)), TRIM(pressi(2)) &
         , p3=pressi_3d)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(pressi(1))//' - '//&
         &TRIM(pressi(2))//' not found',substr)

    CALL get_channel_object(status, TRIM(c_top(1)), TRIM(c_top(2)) &
         , p2=topconv)
    IF (status /= 0) &
         CALL error_bi('channel object '//TRIM(c_top(1))//' - '//&
         &TRIM(c_top(2))//' not found',substr)

    CALL get_channel_info(status, 'scav')
    if (status /= 0.and. scav_trans.eq.2) then
      CALL info_bi('This transport splitting needs the SCAV submodel',substr)
      CALL info_bi('Transport splitting disabled',substr)
      scav_trans = 1
    endif
    if (status /= 0.and. scav_trans.eq.3) then
      CALL info_bi('This transport splitting needs the SCAV submodel',substr)
      CALL info_bi('Transport splitting disabled',substr)
      scav_trans = 1
    endif
    
    CALL info_bi('Looking for tracer quantity.....',substr)
    do js=1,nsets
      SELECT CASE (nsets)
      CASE(1)
        if (lcvt_gp) THEN
          ntrac =  ntrac_gp
          ti    => ti_gp
        ENDIF
        if (lcvt_lg) THEN
          ntrac =  ntrac_lg
          ti    => ti_lg
        ENDIF
      CASE(2)
        if (js == 1)  THEN
          ntrac =  ntrac_gp
          ti    => ti_gp
        ENDIF
        if (js == 2)  THEN
          ntrac =  ntrac_lg
          ti    => ti_lg
        ENDIF
      END SELECT
      do jt=1,ntrac
        tr_attr(js)%amount(jt) = ti(jt)%tp%ident%quantity
      enddo
    enddo

    ! mz_ho_20140826+ ! op_pj_20150427+
    IF (TRIM(pblh(1)) /= '') THEN
       CALL get_channel_object(status,TRIM(pblh(1)),TRIM(pblh(2)),p2=pblh_2d)
       CALL channel_halt(substr, status)
       l_pblh = .TRUE.
    ELSE
       l_pblh = .FALSE.
    ENDIF
    ! mz_ho_20140826- ! op_pj_20150427-

!!#D attila +
!   Coupling to ATTILA for pseudo - lagrangian convective tracer transport
    if (lcvt_lg) then
      CALL get_channel_object(status, 'attila','NCB', p3=NCB)
      IF (status /= 0) &
        CALL error_bi( &
        'no ATTILA channel object for number of cells per gid box found!' &
        , substr)
    endif
!!#D attila -

    IF (l_calc_trma) THEN
      ilev_up = 1
      IF (p_parallel_io) &
        PRINT*, substr, ': Restrict transport to a height below 50 hPa!'
#if !defined(COSMO) && !defined(MESSYDWARF)
      DO jk=1,nvclev
        h_a(jk) = vct(jk)
        h_b(jk) = vct(jk+nvclev)
      ENDDO
#endif
      sfpress     = 1.e5_dp   ! reference pressure of 1000 hPa
      DO jk=1,nlev+1
        hypi(jk)      = h_a(jk) + h_b(jk) * sfpress
      ENDDO
      DO jk=1,nlev
        IF (hypi(jk) < 5000._dp .AND. hypi(jk+1) >= 5000._dp) THEN
          ilev_up = jk
          EXIT
        ENDIF
      END DO
      ilev_up = MAX(2, ilev_up) ! op_pj_20120310 (jk-1 !!!)
      IF (p_parallel_io) &
        PRINT*,                                                         &
        substr, ': Restrict transport matrix to a level number higher than ', &
         ilev_up, ' !'
    ENDIF
    
    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

  END SUBROUTINE cvtrans_init_coupling

!===============================================================================

  SUBROUTINE cvtrans_free_memory
    
    INTEGER   :: js
    INTRINSIC :: ASSOCIATED
    
    IF (ASSOCIATED(tr_attr)) THEN
       do js=1,nsets
          IF (ASSOCIATED(tr_attr(js)%amount)) THEN
             DEALLOCATE(tr_attr(js)%amount)
             NULLIFY(tr_attr(js)%amount)
          END IF
       enddo

       DEALLOCATE(tr_attr)
       NULLIFy(tr_attr)
    END IF

    RETURN
  END SUBROUTINE cvtrans_free_memory

! ==============================================================================

  SUBROUTINE cvtrans_CONVEC(case_opt)

    USE messy_main_grid_def_mem_bi,ONLY: jrow, nproma, kproma, nlev
    USE messy_main_constants_mem,  ONLY: g, m_air
    USE messy_main_timer,          ONLY: time_step_len, lstart  
    USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte,   pxtm1 => qxtm1
    USE messy_main_tracer,         ONLY: I_convect

! grid point tracer convective transport interface

    INTEGER, INTENT(IN) :: case_opt
    REAL(dp), POINTER, DIMENSION(:,:,:) :: pxtp1        => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:) :: xtte_cvtrans => NULL()
    
    INTEGER  :: js, jl, jk, option, jt
    REAL(dp) :: conv_part(nproma), delp(kproma,nlev)

    if (lstart)        RETURN
    if (.not. lcvt_gp) RETURN

    ntrac = ntrac_gp
    ti    => ti_gp
    js    = 1
    l_gp  = .TRUE. 
    l_lg  = .FALSE. 

    conv_part(:) = 0._dp

    IF (L_calc_trma) THEN
      calc_trma = .TRUE.
      ALLOCATE(xtte_cvtrans(kproma,nlev,nlev))
      ALLOCATE(pxtp1(kproma,nlev,nlev)) 
      xtte_cvtrans(:,:,:) = 0._dp
      ! Initialise pseudo tracers for the tropospheric levels
      pxtp1(:,:,:) = 0._dp
      
      do jk=ilev_up,nlev
        do jl=1,kproma
          ! set "mixing ratio" to a level of 1 kg/m² in only one gridbox per level
          ! per tracer, otherwise remain at zero.
          ! pseudo tracers have a molar weight of 1 g/mol

          ! calculate delta p
          delp(jl,jk) = pressi_3d(_RI_XYZ__(jl,jrow,jk+1)) - pressi_3d(_RI_XYZ__(jl,jrow,jk))
          ! calculate tracer mixing ratio
          pxtp1(jl,jk,jk) = 1._dp * g / delp(jl,jk) * m_air
        enddo
      enddo
        
      ! call main routine for closure and transport
      conv_part(:) = 1._dp
      CALL CVTRANS_MAIN(1, 2, jrow, kproma, &
                        pxtp1, xtte_cvtrans, conv_part, js)
      
      ! return new pseudo tracer mixing ratios to determine the transport matrix
      do jk=ilev_up,nlev
        do jt=ilev_up,nlev
          do jl=1,kproma
            ! pseudo tracers have a molar weight of 1 g/mol
            ! calculate new tracer mixing ratio
            pxtp1(jl,jk,jt) = pxtp1(jl,jk,jt) + xtte_cvtrans(jl,jk,jt)*time_step_len
            ! calculate new tracer box mass and store this value in the transport matrix
            trma(jl,jk,jt,jrow) = pxtp1(jl,jk,jt) / m_air * delp(jl,jk) / g
          enddo
        enddo
      end do
      DEALLOCATE(xtte_cvtrans)
      NULLIFY(xtte_cvtrans)
      DEALLOCATE(pxtp1)
      NULLIFY(pxtp1)
      calc_trma = .FALSE.
    ENDIF ! l_calc_trma

    ALLOCATE(xtte_cvtrans(kproma,nlev,ntrac))
    ALLOCATE(pxtp1(kproma,nlev,ntrac)) 

    pxtp1(:,:,:)        = 0._dp
    xtte_cvtrans(:,:,:) = 0._dp


    SELECT CASE (scav_trans)

    CASE(1)
      if (case_opt.eq.1) option = 1
      if (case_opt.eq.2) option = 2
      conv_part(:) = 1._dp

    CASE(2)
      if (case_opt.eq.1) option = 3
      if (case_opt.eq.2) option = 4
      conv_part(:) = 1._dp
      
    CASE(3)
      if (case_opt.eq.1) option = 3
      if (case_opt.eq.2) option = 2
      conv_part(:) = 1._dp

    CASE(4)
      if (case_opt.eq.1) option = 5
      if (case_opt.eq.2) option = 6
      do jk=1,nlev
        do jl=1,kproma
           ! um_ak_20090716+
           !conv_part(jl) = MAX(conv_part(jl),ccover(jl,jk,jrow))
           conv_part(jl) = MAX(conv_part(jl),ccover(_RI_XYZ__(jl,jrow,jk)))
           ! um_ak_20090716-
        enddo
      enddo

    CASE(5)
      if (case_opt.eq.1) option = 7
      if (case_opt.eq.2) option = 8
      do jk=1,nlev
        do jl=1,kproma
           ! um_ak_20090716+
           !conv_part(jl) = MAX(conv_part(jl),ccover(jl,jk,jrow))
           conv_part(jl) = MAX(conv_part(jl),ccover(_RI_XYZ__(jl,jrow,jk)))
           ! um_ak_20090716-
        enddo
      enddo

    END SELECT

#ifndef MESSYTENDENCY
    do jk=1,nlev
      do jt=1,ntrac
        do jl=1,kproma
           ! um_ak_20090716+
!          pxtp1(jl,jk,jt) = MAX(pxtm1(jl,jk,jt) + &
!                            pxtte(jl,jk,jt) * time_step_len, 0._dp)
           pxtp1(jl,jk,jt) = MAX(pxtm1(_RI_X_ZN_(jl,jk,jt)) + &
                             pxtte(_RI_X_ZN_(jl,jk,jt)) * time_step_len, 0._dp)
           ! um_ak_20090716-
        enddo
      enddo
    enddo
#else
    DO jt=1,ntrac
!!$       CALL mtend_get_start_l (mtend_id_tracer &
!!$            , idt = jt, v0 = pxtp1(1:kproma,1:nlev,jt))
       CALL mtend_get_start_l (jt, v0 = pxtp1(1:kproma,1:nlev,jt))
       pxtp1(:,:,jt) = MAX(pxtp1(:,:,jt),0._dp)
    ENDDO
#endif

    DO jt=1,ntrac
      IF (ti(jt)%tp%meta%cask_i(I_convect)==OFF) THEN
        pxtp1(:,:,jt) = 0._dp
      ENDIF
    END DO
    CALL CVTRANS_MAIN(case_opt, option, jrow, kproma, &
                     pxtp1, xtte_cvtrans, conv_part, js)
    DO jt=1,ntrac
      IF (ti(jt)%tp%meta%cask_i(I_convect)==OFF) THEN
        xtte_cvtrans(:,:,jt) = 0._dp
      ENDIF
    END DO

    ! um_ak_20090716+
    !pxtte(1:kproma,1:nlev,1:ntrac) = pxtte(1:kproma,1:nlev,1:ntrac) +  &
    !                                 xtte_cvtrans(1:kproma,1:nlev,1:ntrac)
#ifndef MESSYTENDENCY
    DO jt = 1,ntrac
       DO jk = 1,nlev
          pxtte(_RI_X_ZN_(1:kproma,jk,jt)) = pxtte(_RI_X_ZN_(1:kproma,jk,jt)) +  &
                                    xtte_cvtrans(1:kproma,jk,jt)
       ENDDO
    ENDDO
#else
    DO jt = 1,ntrac
!!$       CALL mtend_add_l (my_handle, mtend_id_tracer &
!!$            , px = xtte_cvtrans(1:kproma,1:nlev,jt), idt=jt)
       CALL mtend_add_l (my_handle, jt &
            , px = xtte_cvtrans(1:kproma,1:nlev,jt))
    ENDDO
#endif
    ! um_ak_20090716-

    DEALLOCATE(xtte_cvtrans)
    NULLIFY(xtte_cvtrans)
    DEALLOCATE(pxtp1)
    NULLIFY(pxtp1)

  END SUBROUTINE cvtrans_CONVEC

!===============================================================================

  SUBROUTINE cvtrans_GLOBAL_END
!!#D attila +
!!!#if !defined(COSMO) && !defined(CESM1)
#if defined(ECHAM5)
    USE messy_main_tracer_mem_bi,  ONLY: pxtte => xtte_a, pxtm1 => xtm1_a, &
                                         NCELL
    USE messy_main_grid_def_mem_bi,ONLY: ngpblks, nlev, nproma, &
                                         npromz ! mz_pj_20080702
    USE messy_main_timer,          ONLY: time_step_len, lstart
    USE messy_attila_tools_e5,     ONLY: lg2gp_e5, LG2GP_AVE, &
                                         lggpte2lgte_e5
    IMPLICIT NONE


    REAL(dp), POINTER, DIMENSION(:,:)     :: pxtp1_lg          => NULL()
    REAL(dp), POINTER, DIMENSION(:,:)     :: xte_cvtrans_lg    => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: pxtp1_lggp        => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:)   :: zxtp1_lggp        => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: xtte_cvtrans_lggp => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: xtte_lggp         => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:)   :: xtte_cvtrans      => NULL()

    INTEGER  :: js, jrow, kproma, case_opt, option
    INTEGER  :: jl, jk, jt, jm
    LOGICAL  :: column_ori
    REAL(dp) :: conv_part(nproma) 

    if (lstart)        RETURN
    if (.not. lcvt_lg) RETURN

    SELECT CASE (nsets)
    CASE(1)
      js = 1
    CASE(2)
      js = 2
    END SELECT
    ntrac = ntrac_lg
    ti    => ti_lg
    l_gp  = .FALSE.
    l_lg  = .TRUE. 

    ALLOCATE(pxtp1_lg(NCELL,ntrac))
    ALLOCATE(xte_cvtrans_lg(NCELL,ntrac))
    ALLOCATE(pxtp1_lggp(nproma,nlev, ntrac, ngpblks))

!   forward conversion 
!           of lagrangian tracer field to pseudo grid point tracer field
!           of lagrangian tracer tendency to pseudo grid point tracer tendency
        
    pxtp1_lg(:,:) = pxtm1(:,1,:,1) + pxtte(:,1,:,1) * time_step_len 

!  p1 value
    call lg2gp_e5(pxtp1_lg, pxtp1_lggp, LG2GP_AVE, fill_value=0._dp)


    ALLOCATE(xtte_cvtrans_lggp(nproma,nlev, ntrac, ngpblks))  
    ALLOCATE(xtte_cvtrans(nproma,nlev,ntrac))
    ALLOCATE(zxtp1_lggp(nproma,nlev, ntrac))

    option   = 2
    case_opt = 1

    do jrow = 1, ngpblks
! mz_pj_20080702+
!!$      IF ( jrow == dcl%ngpblks ) THEN
!!$        kproma = dcl%npromz
!!$      ELSE
!!$        kproma = dcl%nproma
!!$      END IF
      IF ( jrow == ngpblks ) THEN
        kproma = npromz
      ELSE
        kproma = nproma
      END IF
! mz_pj_20080702-

      xtte_cvtrans(:,:,:) = 0._dp
      zxtp1_lggp(1:kproma,:,:)   = pxtp1_lggp(1:kproma,:,:,jrow)
      conv_part(:)        = 1._dp     

      CALL CVTRANS_MAIN(case_opt, option, jrow, kproma, &
                        zxtp1_lggp(1:kproma,:,:),       &
                        xtte_cvtrans(1:kproma,:,:), conv_part, js)

      xtte_cvtrans_lggp(1:kproma,:,:,jrow) = xtte_cvtrans(1:kproma,:,:)
    END DO

!   backward conversion 
!        of lagrangian tracer field to pseudo grid point tracer field
!        of lagrangian tracer tendency to pseudo grid point tracer tendency


    CALL lggpte2lgte_e5(pxtp1_lggp, xtte_cvtrans_lggp, pxtp1_lg, xte_cvtrans_lg)

    pxtte(:,1,:,1) = pxtte(:,1,:,1) + xte_cvtrans_lg(:,:) 

    DEALLOCATE(pxtp1_lg)          ; NULLIFY(pxtp1_lg)
    DEALLOCATE(xte_cvtrans_lg)    ; NULLIFY(xte_cvtrans_lg)
    DEALLOCATE(xtte_cvtrans_lggp) ; NULLIFY(xtte_cvtrans_lggp)
    DEALLOCATE(pxtp1_lggp)        ; NULLIFY(pxtp1_lggp)
    DEALLOCATE(zxtp1_lggp)        ; NULLIFY(zxtp1_lggp)
    DEALLOCATE(xtte_cvtrans)      ; NULLIFY(xtte_cvtrans)

#endif
    RETURN
!!#D attila -
  END SUBROUTINE cvtrans_GLOBAL_END

!===============================================================================

  SUBROUTINE CVTRANS_MAIN(case_opt, option, jrow, kproma, &
                          xtp1,   xtte, conv_part, js)

    USE messy_main_grid_def_mem_bi,ONLY: nlev
    USE messy_main_grid_def_bi,    ONLY: altitude_gnd
    USE messy_main_timer,          ONLY: time_step_len
    USE messy_main_constants_mem,  ONLY: g
    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_tracer,         ONLY: I_convect

    INTEGER,  INTENT(IN)    :: case_opt, option, jrow, kproma, js
!    REAL(dp), INTENT(INOUT) :: xtp1(kproma,nlev,ntrac), xtte(kproma,nlev,ntrac)
    REAL(dp), INTENT(INOUT) :: xtp1(:,:,:), xtte(:,:,:)
    REAL(dp), INTENT(IN)    :: conv_part(kproma)
    
    ! minimum updraft flux to be called convection
    REAL(dp), PARAMETER   :: minflux = 1.e-7_dp
    INTEGER  :: jl, jk, lproma, jt, i, status
    INTEGER  :: j ! um_ak_20090721

    INTEGER  :: ntra
    
    REAL(dp) :: cv_part(kproma)
    REAL(dp) :: cv_uin(kproma,nlev),    cv_din(kproma,nlev)
    REAL(dp) :: cv_eu(kproma,nlev),     cv_du(kproma,nlev)
    REAL(dp) :: cv_ed(kproma,nlev),     cv_dd(kproma,nlev)
    REAL(dp) :: cv_pressi(kproma,nlev+1), cv_prec(kproma,nlev)
    REAL(dp) :: delta_p(kproma,nlev),   kgm2(kproma,nlev)
!    REAL(dp) :: cv_trac(kproma,nlev,ntrac)
    REAL(dp), ALLOCATABLE :: cv_trac(:,:,:)
    INTEGER  :: k_raintop(kproma)
    INTEGER  :: kb(kproma), kt(kproma), lwork(kproma)
    LOGICAL  :: lskip_transport, lskip_closure
!    REAL(dp) :: fracis(kproma,nlev,ntrac)
    REAL(dp), ALLOCATABLE :: fracis(:,:,:)

!   additional for pseudo - lagrangian and squeezing
    REAL(dp) :: lcv_uin(kproma,nlev),    lcv_din(kproma,nlev)
    REAL(dp) :: lcv_eu(kproma,nlev),     lcv_du(kproma,nlev)
    REAL(dp) :: lcv_ed(kproma,nlev),     lcv_dd(kproma,nlev)
    REAL(dp) :: lcv_kgm2(kproma,nlev)
!    REAL(dp) :: lcv_trac(kproma,nlev,ntrac)
    REAL(dp), ALLOCATABLE :: lcv_trac(:,:,:)
    INTEGER  :: kwork(kproma,nlev), kt_help(kproma)
    INTEGER  :: k, llev(kproma), kt_h1(kproma)
    ! um_ak_20090721+
    real(dp) :: mu(kproma,nlev)   ! Updraft mass flux, top of layer (kg/m2/s)
    real(dp) :: eu(kproma,nlev)   ! Entrainment into updraft (kg/m2/s)
    real(dp) :: du(kproma,nlev)   ! Detrainment from updraft (kg/m2/s)
    real(dp) :: md(kproma,nlev)   ! Downdraft mass flux, top of layer (kg/m2/s)
    real(dp) :: ed(kproma,nlev)   ! Entrainment into downdraft (kg/m2/s)
    real(dp) :: dd(kproma,nlev)   ! Dntrainment from downdraft (kg/m2/s)
    ! um_ak_20090721-
    ! mz_ho_20140826+
    REAL(dp)            :: hb(kproma)  ! Height of the cloud base
    ! Mask to check whether cloud base is coupled to a boundary layer
    LOGICAL             :: belowh(kproma)           
    ! Amount of intermediate time steps if lintsteps = .TRUE.
    INTEGER  :: nrsteps 
    REAL(dp) :: dti     ! Intermediate time step
    ! mz_ho_20140826-
!-------------------------------------------
    belowh(:) = .FALSE. ! mz_ho_20140826

    ntra = SIZE(XTP1,3)

    ALLOCATE(fracis(kproma,nlev,ntra))
    ALLOCATE(lcv_trac(kproma,nlev,ntra))
    ALLOCATE(cv_trac(kproma,nlev,ntra))
    
    fracis(1:kproma,1:nlev,1:ntra) = 1._dp
    cv_trac(1:kproma,1:nlev,1:ntra) = 0._dp
! op_pj_20110218+
! NOTES: - kbot, ktop, kraintop not used
!        - mu, eu, du, md, ed, dd must always be initialized
!!$    IF (case_opt == 1) THEN
!!$      kbot(:,jrow)       = 0._dp
!!$      ktop(:,jrow)       = 0._dp
!!$      kraintop(:,jrow)   = 0._dp
! op_pj_20110218-
      ! um_ak_20090716+
      mu = 0._dp
      eu = 0._dp
      du = 0._dp
      md = 0._dp
      ed = 0._dp
      dd = 0._dp
    IF (case_opt == 1) THEN  ! op_pj_20110218 (see above)
      cumassf(_RI_XYZ__(:,jrow,:))  = 0._dp
      cdmassf(_RI_XYZ__(:,jrow,:))  = 0._dp
      cuentr (_RI_XYZ__(:,jrow,:))  = 0._dp
      cudetr (_RI_XYZ__(:,jrow,:))  = 0._dp
      cdentr (_RI_XYZ__(:,jrow,:))  = 0._dp
      cddetr (_RI_XYZ__(:,jrow,:))  = 0._dp
      ! um_ak_20090716-
    ENDIF
    ! determine convective mask
    lproma   = 0
    lwork(:) = 0
    do jl=1,kproma
      IF (NINT(topconv(jl,jrow)) > 0) THEN
        lproma = lproma + 1
        lwork(lproma) = jl
      ENDIF
    enddo
    
    IF (lproma > 0) THEN
      ! pack convective columns
      do jl=1,lproma
        i = lwork(jl)
        cv_part(jl) = conv_part(i)
      enddo
      do jk=1,nlev
        do jl=1,lproma
          i = lwork(jl)
          ! um_ak_20090716+
          cv_uin(jl,jk)    = umassf(_RI_XYZ__(i,jrow,jk)) / conv_part(i)
          cv_din(jl,jk)    = dmassf(_RI_XYZ__(i,jrow,jk)) / conv_part(i)
          cv_eu(jl,jk)     = uentr(_RI_XYZ__(i,jrow,jk))  / conv_part(i)
          cv_du(jl,jk)     = udetr(_RI_XYZ__(i,jrow,jk))  / conv_part(i)
          cv_ed(jl,jk)     = dentr(_RI_XYZ__(i,jrow,jk))  / conv_part(i)
          cv_prec(jl,jk)   = precflx_cv(_RI_XYZ__(i,jrow,jk)) / conv_part(i)
          ! um_ak_20090716-
        enddo
      enddo
      do jk=1,nlev+1
        do jl=1,lproma
          i = lwork(jl)
          ! um_ak_20090716+
          !cv_pressi(jl,jk) = pressi_3d(i,jk,jrow)
          cv_pressi(jl,jk) = pressi_3d(_RI_XYZ__(i,jrow,jk))
          ! um_ak_20090716-
        enddo
      enddo
      SELECT CASE (option)
      CASE(1,2,3,4,7)
        do jt=1,ntra
!          if (ti(jt)%tp%meta%cask_i(I_convect)==ON) THEN
            do jk=1,nlev
              do jl=1,lproma
                i = lwork(jl)
                cv_trac(jl,jk,jt) = xtp1(i,jk,jt)
              enddo
            enddo
!          endif
        enddo
      CASE(5)
        do jt=1,ntra
          if (ti(jt)%tp%meta%cask_i(I_convect)==ON) THEN
            do jk=1,nlev
              do jl=1,kproma
                 ! um_ak_20090716+
                 !trac_field(jl,jk,jt,jrow) = MAX(xtp1(jl,jk,jt), 0._dp)
                 trac_field(_RI_XYZN_(jl,jrow,jk,jt)) = MAX(xtp1(jl,jk,jt), 0._dp)
                 ! um_ak_20090716-
              enddo
            enddo
          end if
        enddo
      CASE(6,8)
        do jt=1,ntra
          if (ti(jt)%tp%meta%cask_i(I_convect)==ON) THEN
            do jk=1,nlev
              do jl=1,lproma
                i = lwork(jl)
                ! um_ak_20090716+
!                cv_trac(jl,jk,jt) = trac_field(i,jk,jt,jrow)
                cv_trac(jl,jk,jt) = trac_field(_RI_XYZN_(i,jrow,jk,jt))
                ! um_ak_20090716-
              enddo
            enddo
          endif
        enddo
      END SELECT

      lskip_closure   = .false.
      lskip_transport = .false.
      
      kgm2(1:lproma,1:nlev) = 0._dp
      k_raintop(1:lproma) = 0

      ! select options from namelist control
      SELECT CASE (option)
      CASE(1,5)
!        lskip_closure   = .true.
        lskip_transport = .true.
      CASE(3,7)
        do jk=nlev,1,-1
          do jl=1,lproma
            IF ( cv_prec(jl,jk) > 1.e-10_dp ) then 
              k_raintop(jl) = jk
            ENDIF
          enddo
        enddo
        do jk=1,nlev
          do jl=1,lproma
            if ( jk > max(1, k_raintop(jl)) ) cycle 
            cv_uin(jl,jk) = 0._dp
            cv_din(jl,jk) = 0._dp
            cv_eu(jl,jk)  = 0._dp
            cv_du(jl,jk)  = 0._dp
            cv_ed(jl,jk)  = 0._dp
          enddo
        enddo
      CASE(4,8)
        do jk=nlev,1,-1
          do jl=1,lproma
            IF ( cv_prec(jl,jk) > 1.e-10_dp ) then 
              k_raintop(jl) = jk
            ENDIF
          enddo
        enddo
        do jk=1,nlev
          do jl=1,lproma
            if ( jk < max(1, k_raintop(jl)) ) cycle  
            cv_uin(jl,jk) = 0._dp
            cv_din(jl,jk) = 0._dp
            cv_eu(jl,jk)  = 0._dp
            cv_du(jl,jk)  = 0._dp
            cv_ed(jl,jk)  = 0._dp
          enddo
        enddo
      END SELECT

      do jk=1,nlev
        do jl=1,lproma
          delta_p(jl,jk) = cv_pressi(jl,jk+1) - cv_pressi(jl,jk) 
          kgm2(jl,jk)    = delta_p(jl,jk)/g 
        enddo
      enddo
      
!---------------------------------
!!#D attila +
      IF (l_lg) THEN
!  squeezing algorithm for lagrangian mapped tracer mixing ratios
!  to avoid zero mixing ratios in certain gridpoint boxes in which no 
!  lagrangian parcel is available
        llev(:)    = nlev + 1
        DO jl=1,lproma
          i = lwork(jl)
          kt_h1(jl)   = MAX(1,NINT(topconv(i,jrow))-2)
        ENDDO

! saving fields in the vertical regions of active convection
        lcv_trac(:,:,:) = 0._dp
        lcv_uin(:,:)    = 0._dp
        lcv_din(:,:)    = 0._dp
        lcv_eu(:,:)     = 0._dp
        lcv_du(:,:)     = 0._dp
        lcv_ed(:,:)     = 0._dp
        lcv_kgm2(:,:)   = 0._dp
        
        DO jt=1,ntra
          DO jk=1,nlev
            DO jl=1,lproma
              lcv_trac(jl,jk,jt) = cv_trac(jl,jk,jt)
            ENDDO
          ENDDO
        ENDDO

        DO jl=1,lproma
          DO jk=kt_h1(jl),nlev
            lcv_uin(jl,jk)    = cv_uin(jl,jk) 
            lcv_din(jl,jk)    = cv_din(jl,jk) 
            lcv_eu(jl,jk)     = cv_eu(jl,jk)  
            lcv_du(jl,jk)     = cv_du(jl,jk)  
            lcv_ed(jl,jk)     = cv_ed(jl,jk)  
            lcv_kgm2(jl,jk)   = kgm2(jl,jk)
          ENDDO
        ENDDO

        cv_trac(:,:,:) = 0._dp
        cv_uin(:,:)    = 0._dp
        cv_din(:,:)    = 0._dp
        cv_eu(:,:)     = 0._dp
        cv_du(:,:)     = 0._dp
        cv_ed(:,:)     = 0._dp
        kgm2(:,:)      = 0._dp
        

        kwork(:,:) = 0
        do jk=nlev,1,-1
          do jl=1,lproma
            i = lwork(jl)
            ! um_ak_20090716+
            !IF ( NINT(NCB(i,jk,jrow)) /= 0 ) THEN
            IF ( NINT(NCB(_RI_XYZ__(i,jrow,jk))) /= 0 ) THEN
            ! um_ak_20090716-
              llev(jl) = llev(jl) - 1
              kwork(jl,llev(jl)) = jk
            ENDIF
          ENDDO
        ENDDO

        DO jl=1,lproma
          DO jk=llev(jl),nlev
            k=kwork(jl,jk)
            cv_trac(jl,jk,1:ntra) = lcv_trac(jl,k,1:ntra)
            cv_uin(jl,jk)    =  lcv_uin(jl,k) 
            cv_din(jl,jk)    =  lcv_din(jl,k) 
            cv_eu(jl,jk)     =  lcv_eu(jl,k)  
            cv_du(jl,jk)     =  lcv_du(jl,k)  
            cv_ed(jl,jk)     =  lcv_ed(jl,k)  
            kgm2(jl,jk)      =  lcv_kgm2(jl,k)
          ENDDO
        ENDDO    

        ! locate level index of potential ktop
        kt_help(:) = nlev + 1
        DO jl=1,lproma
          DO jk=nlev,llev(jl),-1
            IF (cv_uin(jl,jk) > minflux) kt_help(jl) = jk
          ENDDO
        ENDDO

        DO jl=1,lproma

          DO WHILE ( (kt_help(jl) < nlev + 1) .AND. &
                     (kgm2(jl,kt_help(jl)-2) < SQRT(TINY(0._dp))) )
             
            cv_uin(jl,kt_help(jl)) = 0._dp
            cv_din(jl,kt_help(jl)) = 0._dp
            cv_eu(jl,kt_help(jl))  = 0._dp
            cv_du(jl,kt_help(jl))  = 0._dp
            cv_ed(jl,kt_help(jl))  = 0._dp
            kt_help(jl) = kt_help(jl) + 1
          ENDDO
            
          IF (kt_help(jl) == nlev) THEN
            cv_uin(jl,:) = 0._dp
            cv_din(jl,:) = 0._dp
            cv_eu(jl,:)  = 0._dp
            cv_du(jl,:)  = 0._dp
            cv_ed(jl,:)  = 0._dp
          ENDIF

          IF (cv_uin(jl,nlev) > 0._dp) THEN
            cv_uin(jl,nlev) = 0._dp
            cv_din(jl,nlev) = 0._dp
            cv_eu(jl,nlev)  = 0._dp
            cv_du(jl,nlev)  = 0._dp
            cv_ed(jl,nlev)  = 0._dp
          ENDIF

        ENDDO
      ENDIF   ! l_lg
!!#D attila -
!---------------------------------
          
      IF ( .not.lskip_closure ) then ! um_ak_20090721: then
        CALL closure(nlev, lproma, time_step_len,                     &
                     kgm2(1:lproma,1:nlev),                           &
                     cv_uin(1:lproma,1:nlev), cv_eu(1:lproma,1:nlev), &
                     cv_du(1:lproma,1:nlev),  cv_din(1:lproma,1:nlev),&
                     cv_ed(1:lproma,1:nlev),  cv_part(1:lproma),      &
                     jrow, kb(1:lproma), kt(1:lproma),                &
                     lwork(1:lproma),                                 &
                     ! um_ak_20090721+
                     mu(1:lproma,1:nlev), md(1:lproma,1:nlev),        &
                     eu(1:lproma,1:nlev), du(1:lproma,1:nlev),        &
                     ed(1:lproma,1:nlev), dd(1:lproma,1:nlev),        &
                     ! um_ak_20090721-
                     status )

        ! mz_ho_20140826+
        if (lcvt_sc) then
           do jl = 1, lproma
              hb(jl) = altitude_gnd(lwork(jl),kb(jl),jrow)
           end do

           if (l_pblh) then
              do jl = 1, lproma
                 belowh(jl) = ( hb(jl) .le. &
                      max( pblh_2d(lwork(jl), jrow), hlimit ) )
              end do
           else ! l_pblh
              belowh(1:lproma) = (hb(1:lproma) .le. hlimit)
           endif ! l_pblh
        endif !lcvt_sc
        ! mz_ho_20140826-

! um_ak_20090721+: moved from SMCL
      if (segtrans) then
         ! some comments, see code
      else
         do jk=1,nlev
           do i=1,lproma
             j = lwork(i)
             cumassf(_RI_XYZ__(j,jrow,jk)) = mu(i,jk)
             cdmassf(_RI_XYZ__(j,jrow,jk)) = md(i,jk)
             cuentr (_RI_XYZ__(j,jrow,jk)) = eu(i,jk)
             cudetr (_RI_XYZ__(j,jrow,jk)) = du(i,jk)
             cdentr (_RI_XYZ__(j,jrow,jk)) = ed(i,jk)
             cddetr (_RI_XYZ__(j,jrow,jk)) = dd(i,jk)
           enddo
         enddo
      end if
! um_ak_20090721-
   end IF ! um_ak_20090721

      if (status /= 0) CALL error_bi("Error in closure", modstr)
!      print*, "finished closure in jrow: ", jrow

!      IF (l_lg) lskip_transport =.true.
      
      IF (.not. lskip_transport) then ! um_ak_20090721: then

        ! um_ak_20090721+: moved from SMCL
         do jk=1,nlev
            do i=1,lproma
               j = lwork(i)
               mu(i,jk) = cumassf(_RI_XYZ__(j,jrow,jk)) 
               md(i,jk) = cdmassf(_RI_XYZ__(j,jrow,jk)) 
               eu(i,jk) = cuentr (_RI_XYZ__(j,jrow,jk)) 
               du(i,jk) = cudetr (_RI_XYZ__(j,jrow,jk)) 
               ed(i,jk) = cdentr (_RI_XYZ__(j,jrow,jk)) 
               dd(i,jk) = cddetr (_RI_XYZ__(j,jrow,jk)) 
            enddo
         enddo
        ! um_ak_20090721-

! mz_ho_20140826+ ! op_pj_20150427
         IF (.NOT. lintsteps) then
           ! classic:
            CALL cv_transport(nlev, lproma, ntra, time_step_len,           &
                          kgm2(1:lproma,1:nlev), fracis(1:lproma,1:nlev,:), &
                          cv_trac(1:lproma,1:nlev,:),                       &
                          jrow, kb(1:lproma), kt(1:lproma), status,         &
                          lwork(1:lproma), option, js,                      &
                          ! um_ak_20090721+
                          mu(1:lproma,1:nlev), md(1:lproma,1:nlev),        &
                          eu(1:lproma,1:nlev), du(1:lproma,1:nlev),        &
                          ed(1:lproma,1:nlev), dd(1:lproma,1:nlev)         &
                          ! um_ak_20090721-
                          )

         ELSE
            ! supsteps:
            lproma_loop: do jl=1,lproma
               ! Loop over cv_transport, depending on ratio of 
               ! mass flux * dt / total mass in grid cell
               nrsteps = 1 + floor(maxval(mu(jl,2:nlev)*&
                    time_step_len/&
                    (min(kgm2(jl,2:nlev),kgm2(jl,1:(nlev-1)))*maxfrac)))
               ! makes sure that dti * mu < maxfrac * the mass in the 
               ! grid cells above and below the interface level

               dti     = time_step_len / nrsteps

               nrsteps_loop: do i = 1, nrsteps

                  CALL cv_transport(nlev, 1, ntra, dti,                 &
                       kgm2(jl:jl,1:nlev), fracis(jl:jl,1:nlev,:), &
                       cv_trac(jl:jl,1:nlev,:),                       &
                       jrow, kb(jl:jl), kt(jl:jl), status,         &
                       lwork(jl:jl), option, js,                      &
                       ! um_ak_20090721+
                       mu(jl:jl,1:nlev), md(jl:jl,1:nlev),        &
                       eu(jl:jl,1:nlev), du(jl:jl,1:nlev),        &
                       ed(jl:jl,1:nlev), dd(jl:jl,1:nlev)         &
                       ! um_ak_20090721-
                       , belowh(jl:jl) & ! mz_ho_20140826
                       )

               enddo nrsteps_loop
            enddo lproma_loop
         END IF
! mz_ho_20140826-

     END IF ! um_ak_20090721

      if (status /= 0) CALL error_bi("Severe error in cv_transport", modstr)

!---------------------------------
!!#D attila +
      IF (l_lg) THEN
!     un-squeeze the vertical columns to the original grid

        lcv_trac(:,:,:) = cv_trac(:,:,:)

        do jt=1,ntra
!          if (ti(jt)%tp%meta%cask_i(I_convect)==ON) THEN
            do jk=1,nlev
              do jl=1,lproma
                i = lwork(jl)
                cv_trac(jl,jk,jt) = xtp1(i,jk,jt)
              ENDDO
            ENDDO
          !ENDIF
        ENDDO
            
        DO jl=1,lproma
          DO jk=llev(jl),nlev
            k=kwork(jl,jk)
            cv_trac(jl,k,1:ntra) = lcv_trac(jl,jk,1:ntra)
          ENDDO
        ENDDO
        
      ENDIF    ! l_lg
!!#D attila -
!---------------------------------

      ! unpack convective columns

      SELECT CASE(option)
      CASE(2,3,4)
        do jt=1,ntra
!          if (ti(jt)%tp%meta%cask_i(I_convect)==ON) THEN
            do jk=1,nlev
              do jl=1,lproma
                i=lwork(jl)
                xtte(i,jk,jt) = ( cv_trac(jl,jk,jt) - xtp1(i,jk,jt) ) &
                                / time_step_len
              enddo
            enddo
!          endif
        enddo
      CASE(6,8)
        do jt=1,ntra
!          if (ti(jt)%tp%meta%cask_i(I_convect)==ON) THEN
            do jk=1,nlev
              do jl=1,lproma
                i=lwork(jl)
                cv_trac(jl,jk,jt) = (1._dp - cv_part(jl)) * xtp1(i,jk,jt) + &
                                    cv_part(jl) * cv_trac(jl,jk,jt)
                xtte(i,jk,jt)     = ( cv_trac(jl,jk,jt) - xtp1(i,jk,jt) ) &
                                    / time_step_len
              enddo
            enddo
!          end if
        end do
      CASE(7)
        do jt=1,ntra
          if (ti(jt)%tp%meta%cask_i(I_convect)==ON) THEN
            do jk=1,nlev
              do jl=1,lproma
                i = lwork(jl)
                ! um_ak_20090716+
                !trac_field(i,jk,jt,jrow) = cv_trac(jl,jk,jt)
                trac_field(_RI_XYZN_(i,jrow,jk,jt)) = cv_trac(jl,jk,jt)
                ! um_ak_20090716-
              enddo
            enddo
          endif
        enddo 
      END SELECT
    END IF

    DEALLOCATE(fracis)
    DEALLOCATE(lcv_trac)
    DEALLOCATE(cv_trac)
    
  END SUBROUTINE CVTRANS_MAIN

!===============================================================================

END MODULE MESSY_cvtrans_si
