! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL cloudopt
!
! Author : Simone Dietmueller, DLR
!
! References: see messy_cloudopt.f90
!
! **********************************************************************
#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_cloudopt_si
! **********************************************************************
#if defined(ECHAM5) || defined (CESM1) || defined (MBM_RAD)

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_grid_def_mem_bi, ONLY: nproma,ngpblks,nlev
  USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi  &
                                      , error_bi, info_bi, warning_bi
  USE messy_main_channel,        ONLY: t_chaobj_cpl
  USE messy_main_tools,          ONLY: t_reset_par  &
                                     , PTR_3D_ARRAY &
                                     , PTR_5D_ARRAY
  USE messy_main_constants_mem,  ONLY: STRLEN_ULONG

  ! SMCL
  USE messy_cloudopt              ! nsw, jpband, ...

  IMPLICIT NONE
  !

  INTEGER, PARAMETER:: NMAXCALL=40

  INTEGER, PUBLIC, SAVE :: NCALL
  
  INTEGER, PARAMETER :: id_clc   = 1
  INTEGER, PARAMETER :: id_lw    = 2
  INTEGER, PARAMETER :: id_iw    = 3
  INTEGER, PARAMETER :: id_cdnc  = 4
  INTEGER, PARAMETER :: id_radlp = 5
  INTEGER, PARAMETER :: id_radip = 6
  !
  INTEGER, PARAMETER:: NMAXOBJ=6
  TYPE t_chaobj_cpl_long
     CHARACTER(LEN=STRLEN_ULONG)  :: cha  = ''
     CHARACTER(LEN=STRLEN_ULONG ) :: obj  = ''
  END TYPE t_chaobj_cpl_long

  TYPE t_cld_work
     INTEGER                                       :: n
     TYPE(t_chaobj_cpl), POINTER, DIMENSION(:)     :: name => NULL()
     TYPE(PTR_3D_ARRAY), POINTER, DIMENSION(:)     :: ptr  => NULL()
  END TYPE t_cld_work
   
  TYPE(t_chaobj_cpl_long), DIMENSION(NMAXOBJ,NMAXCALL), SAVE :: cld_inp
  TYPE(t_cld_work),   DIMENSION(NMAXOBJ,NMAXCALL),      SAVE :: xcldin
 
  TYPE t_cld_out
       REAL(DP), DIMENSION(:,:),       POINTER :: pclcv   =>NULL()
       REAL(DP), DIMENSION(:,:,:),     POINTER :: pidxcld =>NULL()
       REAL(DP), DIMENSION(:,:,:),     POINTER :: psumcov  =>NULL()
       REAL(DP), DIMENSION(:,:,:),     POINTER :: pvalid_sw  =>NULL()
       REAL(DP), DIMENSION(:,:,:,:),   POINTER :: ptaucld_lw =>NULL()
       REAL(DP), DIMENSION(:,:,:,:),   POINTER :: ptaucld_sw =>NULL()
       REAL(DP), DIMENSION(:,:,:,:),   POINTER :: ptaucld_raw_sw =>NULL()
       REAL(DP), DIMENSION(:,:,:,:),   POINTER :: pgamma_sw  =>NULL()
       REAL(DP), DIMENSION(:,:,:,:),   POINTER :: pomega_sw  =>NULL()
       ! diagnostics for perturbation
       TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: cov_pert => NULL()
       TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: tau_vis_pert => NULL()
       TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: ice_pert => NULL()

       ! diagnostics for isccp simulator
       REAL(DP), DIMENSION(:,:,:),     POINTER :: isccp_cldtau3d =>NULL()
       REAL(DP), DIMENSION(:,:,:),     POINTER :: isccp_cldemi3d =>NULL()
       REAL(DP), DIMENSION(:,:,:),     POINTER :: isccp_f3d      =>NULL() 
  END TYPE t_cld_out 

 TYPE(t_cld_out),   DIMENSION(NMAXCALL),      SAVE   :: xcldout
  
 ! auxiliary fields
 TYPE(PTR_5D_ARRAY), DIMENSION(:), POINTER :: ptaucld_lw_5d =>NULL()
 TYPE(PTR_5D_ARRAY), DIMENSION(:), POINTER :: ptaucld_sw_5d =>NULL()
 TYPE(PTR_5D_ARRAY), DIMENSION(:), POINTER :: ptaucld_raw_sw_5d =>NULL()
 TYPE(PTR_5D_ARRAY), DIMENSION(:), POINTER :: pgamma_sw_5d  =>NULL()
 TYPE(PTR_5D_ARRAY), DIMENSION(:), POINTER :: pomega_sw_5d  =>NULL()

 ! SUBROUTINES
 ! MAIN ENTRY POINTS
 PUBLIC :: cloudopt_initialize
 PUBLIC :: cloudopt_init_memory
 PUBLIC :: cloudopt_init_coupling
 PUBLIC :: cloudopt_radiation
 PUBLIC :: cloudopt_free_memory   
!PRIVATE :: cldopt_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE cloudopt_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------
    ! BMIL
    USE messy_main_mpi_bi,          ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_grid_def_mem_bi, ONLY: nn, lmidatm
    USE messy_main_data_bi,   ONLY: lcouple
    USE messy_main_tools,     ONLY: find_next_free_unit, strcrack
    USE messy_main_channel,   ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
     
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloudopt_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: j1,j2 
    INTEGER                     :: ncha, nobj, n
    CHARACTER(LEN=STRLEN_CHANNEL), DIMENSION(:), POINTER :: str_cha => NULL()
    CHARACTER(LEN=STRLEN_OBJECT),  DIMENSION(:), POINTER :: str_obj => NULL()
 
   CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

   !READ CTRL namelist  ============================================
   IF (p_parallel_io) THEN                  ! read only on I/O-PE
      iou = find_next_free_unit(100,200)    ! find free I/O unit
      CALL cloudopt_read_nml_ctrl(status, iou)  ! read CTRL-namelist
      ! terminate if error
      IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
   END IF

   CALL p_bcast (rset_asic%l,    p_io)
   CALL p_bcast (rset_asic%v,    p_io)
   CALL p_bcast (rset_zinhomi%l, p_io)
   CALL p_bcast (rset_zinhomi%v, p_io)
   CALL p_bcast (rset_zinhoml%l, p_io)
   CALL p_bcast (rset_zinhoml%v, p_io)
   CALL p_bcast (rset_zinpar%l,  p_io)
   CALL p_bcast (rset_zinpar%v,  p_io)
   CALL p_bcast (loice,          p_io)

   ! READ CPL namelist  ================================================
   IF (p_parallel_io) THEN                  ! read only on I/O-PE
      iou = find_next_free_unit(100,200)    ! find next free I/O unit
      CALL cldopt_read_nml_cpl(status, iou)    ! read CPL-namelist
      ! terminate if error
      IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
   END IF

   CALL p_bcast(NCALL,p_io)

   DO j2=1,NCALL
      DO j1=1,NMAXOBJ
         CALL p_bcast(cld_inp(j1,j2)%cha, p_io)
         CALL p_bcast(cld_inp(j1,j2)%obj, p_io)
      END DO
   END DO

   DO j2=1,NCALL
      DO j1=1,NMAXOBJ
         IF(TRIM(cld_inp(j1,j2)%cha)== '') THEN
            cld_inp(j1,j2)%cha=cld_inp(j1,1)%cha
         ENDIF
         IF(TRIM(cld_inp(j1,j2)%obj)== '') THEN
            cld_inp(j1,j2)%obj=cld_inp(j1,1)%obj
         ENDIF
      END DO
   END DO

   DO j2=1,NCALL
      DO j1=1,NMAXOBJ
         CALL strcrack(cld_inp(j1,j2)%cha, ';', str_cha, ncha)
         CALL strcrack(cld_inp(j1,j2)%obj, ';', str_obj, nobj)
         IF (ncha /= nobj) THEN
            CALL error_bi('number of channels and objects mismatch',substr)
         ENDIF
         xcldin(j1,j2)%n = ncha
         ALLOCATE(xcldin(j1,j2)%ptr(ncha))
         ALLOCATE(xcldin(j1,j2)%name(ncha))
         DO n=1, ncha
            NULLIFY(xcldin(j1,j2)%ptr(n)%ptr)
            xcldin(j1,j2)%name(n)%cha = TRIM(str_cha(n))
            xcldin(j1,j2)%name(n)%obj = TRIM(str_obj(n))
            CALL info_bi('channel / object '//&
                 &TRIM(xcldin(j1,j2)%name(n)%cha)//' / '//&
                 &TRIM(xcldin(j1,j2)%name(n)%obj)//' requested',substr)
         END DO
      END DO
   END DO
   IF (ASSOCIATED(str_cha)) THEN
      DEALLOCATE(str_cha) ; NULLIFY(str_cha)
   ENDIF
   IF (ASSOCIATED(str_obj)) THEN
      DEALLOCATE(str_obj) ; NULLIFY(str_obj)
   ENDIF

   ! set resolution dependent cloud parameters 
   CALL cloud_set_param(status,nn,nlev,lmidatm, lcouple)
   IF (status /= 0) THEN
      CALL error_bi(' error in cloud_set_param',substr)    ! terminate if error
   END IF

   CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

 END SUBROUTINE cloudopt_initialize
 ! ====================================================================

 ! ====================================================================
 SUBROUTINE cloudopt_init_memory

   ! BMIL
   USE messy_main_channel_error_bi,    ONLY: channel_halt
   USE messy_main_channel_bi,          ONLY: GP_3D_MID,DC_GP,                 &
                                             DIMID_LON,DIMID_LEV,DIMID_LAT,   &
                                             GP_2D_HORIZONTAL,                & 
                                             gp_nseg,gp_start,gp_cnt, gp_meml,&
                                             gp_memu

   USE messy_main_channel,             ONLY: new_channel,new_channel_object,  &
                                             new_attribute

   USE messy_main_channel_dimensions,  ONLY: new_dimension
   USE messy_main_channel_repr,        ONLY: new_representation,AUTO,        &
                                             set_representation_decomp,      &
                                             IRANK,PIOTYPE_COL &
                                           , REPR_DEF_AXES
   USE messy_main_constants_mem,       ONLY: STRLEN_MEDIUM
   USE messy_main_tools,               ONLY: int2str


   IMPLICIT NONE

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER            :: substr = 'cloudopt_init_memory'
   CHARACTER(LEN=2)                       :: idx = ''
   CHARACTER(LEN=STRLEN_MEDIUM)           :: sname = ''
   CHARACTER(LEN=2)                       :: call_idx = ''    
   INTEGER                                :: ji, j2, n, j1, nmax
   INTEGER                                :: status
   INTEGER                                :: DIMID_JPBAND, DIMID_NSW
   INTEGER                                :: REPR_OPT_4D_JPBAND &
                                           , REPR_OPT_4D_NSW
   INTEGER, DIMENSION(:,:), POINTER       :: start =>NULL()
   INTEGER, DIMENSION(:,:), POINTER       :: cnt =>NULL()
   INTEGER, DIMENSION(:,:), POINTER       :: meml =>NULL()
   INTEGER, DIMENSION(:,:), POINTER       :: memu =>NULL()
   ! auxiliary pointer 
   REAL(DP), DIMENSION(:,:,:,:),  POINTER :: mem =>NULL()
    
   ! CHANNEL AND CHANNEL OBJECTS
   CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

   ! 4D NSW representation / decomposition

   CALL new_dimension(status, DIMID_NSW, 'OPT_NSW', nsw)
   CALL channel_halt(substr, status)

   CALL new_representation(status, REPR_OPT_4D_NSW, 'REPR_OPT_4D_NSW'        &
        , rank = 4, link = 'xxxx', dctype = DC_GP                            &
        , dimension_ids =                                                    &
        (/ _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_NSW) /)          &
        , ldimlen   = (/ _RI_XYZN_(nproma , ngpblks, AUTO, AUTO) /)          &
        , output_order  = (/ _IN_XYZN_, _IX_XYZN_ , _IY_XYZN_, _IZ_XYZN_ /)  &
        , axis =  repr_def_axes(_RI_XYZN_('X','Y','Z','N'))           )
   CALL channel_halt(substr, status)

   ALLOCATE(start(gp_nseg,IRANK))
   ALLOCATE(cnt(gp_nseg,IRANK))
   ALLOCATE(meml(gp_nseg,IRANK))
   ALLOCATE(memu(gp_nseg,IRANK))
   
   start(:,:) = gp_start(:,:)
   cnt(:,:) = gp_cnt(:,:)
   meml(:,:) = gp_meml(:,:)
   memu(:,:) = gp_memu(:,:)
   
   cnt(:, _IN_XYZN_) = nsw
   memu(:, _IN_XYZN_) = nsw
   
   CALL set_representation_decomp(status, REPR_OPT_4D_NSW &
        , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
   CALL channel_halt(substr, status)
   
   DEALLOCATE(start) ; NULLIFY(start)
   DEALLOCATE(cnt)   ; NULLIFY(cnt)
   DEALLOCATE(meml)  ; NULLIFY(meml)
   DEALLOCATE(memu)  ; NULLIFY(memu)
   
   ! 4D lw BAND represenation / decomposition

   CALL new_dimension(status, DIMID_JPBAND, 'OPT_JPBAND', jpband)
   CALL channel_halt(substr, status)
   CALL new_representation(status, REPR_OPT_4D_JPBAND                       &
        ,                         'REPR_OPT_4D_JPBAND'                      &
        , rank = 4, link = 'xxxx', dctype = DC_GP                           &
        , dimension_ids =                                                   &
        (/  _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_JPBAND) /)     &
        , ldimlen = (/ _RI_XYZN_(nproma , ngpblks, AUTO, AUTO) /)           &
        , output_order = (/ _IN_XYZN_, _IX_XYZN_ , _IY_XYZN_, _IZ_XYZN_ /)  &
        , axis =  repr_def_axes(_RI_XYZN_('X','Y','Z','N'))                 &
        )
   CALL channel_halt(substr, status)
   
   ALLOCATE(start(gp_nseg,IRANK))
   ALLOCATE(cnt(gp_nseg,IRANK))
   ALLOCATE(meml(gp_nseg,IRANK))
   ALLOCATE(memu(gp_nseg,IRANK))
   
   start(:,:) = gp_start(:,:)
   cnt(:,:) = gp_cnt(:,:)
   meml(:,:) = gp_meml(:,:)
   memu(:,:) = gp_memu(:,:)
   
   cnt(:, _IN_XYZN_) = jpband
   memu(:, _IN_XYZN_) = jpband
   
   CALL set_representation_decomp(status, REPR_OPT_4D_JPBAND &
        , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
   CALL channel_halt(substr, status)
   
   DEALLOCATE(start) ; NULLIFY(start)
   DEALLOCATE(cnt)   ; NULLIFY(cnt)
   DEALLOCATE(meml)  ; NULLIFY(meml)
   DEALLOCATE(memu)  ; NULLIFY(memu)

   ALLOCATE(ptaucld_lw_5d(NCALL))
   ALLOCATE(ptaucld_sw_5d(NCALL))
   ALLOCATE(ptaucld_raw_sw_5d(NCALL))
   ALLOCATE(pgamma_sw_5d(NCALL))
   ALLOCATE(pomega_sw_5d(NCALL))
   DO j2=1,NCALL
      NULLIFY(ptaucld_lw_5d(j2)%ptr)
      NULLIFY(ptaucld_sw_5d(j2)%ptr)
      NULLIFY(ptaucld_raw_sw_5d(j2)%ptr)
      NULLIFY(pgamma_sw_5d(j2)%ptr)
      NULLIFY(pomega_sw_5d(j2)%ptr)
   END DO

   DO j2=1,NCALL
      CALL int2str(call_idx, j2, '0')
      sname=modstr//call_idx
      
      CALL new_channel(status, TRIM(sname), lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, TRIM(sname), 'idx_cld' &
           , p3=xcldout(j2)%pidxcld, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'idx_cld', &
           'long_name', c='clear/cloudy index (input rad)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'idx_cld', &
           'units', c='-- ')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, TRIM(sname), 'sum_cov' &
           , p3=xcldout(j2)%psumcov, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'sum_cov', &
           'long_name', c='effective cloud cover (natural + perturbation(s) (input rad))')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'sum_cov', &
           'units', c='-- ')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, TRIM(sname), 'valid_sw' &
           , p3=xcldout(j2)%pvalid_sw, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'valid_sw', &
           'long_name', c='validity flag for SW calculation')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'valid_sw', &
           'units', c='-- ')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, TRIM(sname), 'clcv' &
           , p2=xcldout(j2)%pclcv, reprid=GP_2D_HORIZONTAL)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'clcv', &
           'long_name', c='total cloud cover (natural + perturbation(s) (input rad))')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'clcv', &
           'units', c='-- ')
      CALL channel_halt(substr, status)

      ! for ISCCP
      CALL new_channel_object(status, TRIM(sname), 'isccp_cldtau' &
           , p3=xcldout(j2)%isccp_cldtau3d, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'isccp_cldtau' &
           ,'long_name', c='ISCCP cloud optical thickness')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'isccp_cldtau' &
           ,'units', c='-')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, TRIM(sname), 'isccp_cldemi' &
           , p3=xcldout(j2)%isccp_cldemi3d, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'isccp_cldemi' &
           ,'long_name', c='ISCCP cloud emissivity')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'isccp_cldemi' &
           ,'units', c='-')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, TRIM(sname), 'isccp_f' &
           , p3=xcldout(j2)%isccp_f3d, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'isccp_f' &
           ,'long_name', c='ISCCP cloud fraction')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'isccp_f' &
           ,'units', c='-')
      CALL channel_halt(substr, status)

      ALLOCATE(ptaucld_lw_5d(j2)%ptr(_RI_XYZN_(nproma,ngpblks,nlev,jpband),1))
      ALLOCATE(ptaucld_sw_5d(j2)%ptr(_RI_XYZN_(nproma,ngpblks,nlev,nsw),1))
      ALLOCATE(ptaucld_raw_sw_5d(j2)%ptr(_RI_XYZN_(nproma,ngpblks,nlev,nsw),1))
      ALLOCATE(pgamma_sw_5d(j2)%ptr(_RI_XYZN_(nproma,ngpblks,nlev,nsw),1))
      ALLOCATE(pomega_sw_5d(j2)%ptr(_RI_XYZN_(nproma,ngpblks,nlev,nsw),1))

      ! INIT
      ptaucld_lw_5d(j2)%ptr(:,:,:,:,:)     = 0._dp
      ptaucld_sw_5d(j2)%ptr(:,:,:,:,:)     = 0._dp
      ptaucld_raw_sw_5d(j2)%ptr(:,:,:,:,:) = 0._dp
      pgamma_sw_5d(j2)%ptr(:,:,:,:,:)      = 0._dp
      pomega_sw_5d(j2)%ptr(:,:,:,:,:)      = 0._dp

      mem=>ptaucld_lw_5d(j2)%ptr(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(sname), 'tau_cld_lw' &
           , p4=xcldout(j2)%ptaucld_lw, reprid=REPR_OPT_4D_JPBAND, &
           mem=mem, lrestreq=.FALSE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'tau_cld_lw', &
           'long_name', c='cloud optical depth')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'tau_cld_lw', &
           'units', c='-- ')
      CALL channel_halt(substr, status)
    
      ! for output 
      DO ji=1,jpband
         CALL int2str(idx,ji,'0')
         mem=>ptaucld_lw_5d(j2)%ptr(_RI_XYZN_(:,:,:,ji),:)
         CALL new_channel_object(status, TRIM(sname),'tau_cld_lw_B'//idx &
              , reprid=GP_3D_MID, mem=mem)
      END DO

      mem=>ptaucld_sw_5d(j2)%ptr(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(sname), 'tau_cld_sw' &
           , p4=xcldout(j2)%ptaucld_sw, reprid=REPR_OPT_4D_NSW, &
           mem=mem, lrestreq=.FALSE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'tau_cld_sw', &
           'long_name', c='extinction tau')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'tau_cld_sw', &
           'units', c='-- ')
      CALL channel_halt(substr, status)
   
      ! for output
      DO ji=1,nsw
         CALL int2str(idx, ji, '0')
         mem=>ptaucld_sw_5d(j2)%ptr(_RI_XYZN_(:,:,:,ji),:)
         CALL new_channel_object(status, TRIM(sname),'tau_cld_sw_B'//idx &
              , reprid=GP_3D_MID, mem=mem)
      END DO

      mem=>ptaucld_raw_sw_5d(j2)%ptr(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(sname), 'tau_cld_raw_sw' &
           , p4=xcldout(j2)%ptaucld_raw_sw, reprid=REPR_OPT_4D_NSW, &
           mem=mem, lrestreq=.FALSE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'tau_cld_raw_sw', &
           'long_name', c='extinction tau')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'tau_cld_raw_sw', &
           'units', c='-- ')
      CALL channel_halt(substr, status)
   
      ! for output
      DO ji=1,nsw
         CALL int2str(idx, ji, '0')
         mem=>ptaucld_raw_sw_5d(j2)%ptr(_RI_XYZN_(:,:,:,ji),:)
         CALL new_channel_object(status, TRIM(sname),'tau_cld_raw_sw_B'//idx &
              , reprid=GP_3D_MID, mem=mem)
      END DO

      mem=>pgamma_sw_5d(j2)%ptr(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(sname), 'gamma_cld_sw' &
           , p4=xcldout(j2)%pgamma_sw, reprid=REPR_OPT_4D_NSW, &
           mem=mem, lrestreq=.FALSE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'gamma_cld_sw', &
           'long_name', c='asymmerty factor')
      CALL channel_halt(substr, status)
      CALL new_attribute(status,TRIM(sname), 'gamma_cld_sw', &
           'units', c='-- ')
      CALL channel_halt(substr, status)

      ! for output
      DO ji=1,nsw
         CALL int2str(idx, ji, '0')
         mem=>pgamma_sw_5d(j2)%ptr(_RI_XYZN_(:,:,:,ji),:)
         CALL new_channel_object(status, TRIM(sname),'gamma_cld_sw_B'//idx &
              , reprid=GP_3D_MID, mem=mem)
      END DO

      mem=>pomega_sw_5d(j2)%ptr(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(sname), 'omega_cld_sw' &
           , p4=xcldout(j2)%pomega_sw, reprid=REPR_OPT_4D_NSW,&
           mem=mem, lrestreq=.FALSE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'omega_cld_sw', &
           'long_name', c='single scattering albedo')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(sname), 'omega_cld_sw', &
           'units', c='-- ')
      CALL channel_halt(substr, status)

      ! for output
      DO ji=1,nsw
         CALL int2str(idx, ji, '0')
         mem=>pomega_sw_5d(j2)%ptr(_RI_XYZN_(:,:,:,ji),:)
         CALL new_channel_object(status, TRIM(sname),'omega_cld_sw_B'//idx &
              , reprid=GP_3D_MID, mem=mem)
      END DO

      ! DIAGNOSTIC TAU (VISIBLE) OF ALL PERTURBATIONS
      ! get max. number of perturbations of this call
      nmax = 1
      DO j1=1,NMAXOBJ
         nmax = MAX(nmax,xcldin(j1,j2)%n)
      END DO

      ALLOCATE(xcldout(j2)%cov_pert(nmax))
      ALLOCATE(xcldout(j2)%tau_vis_pert(nmax))
      ALLOCATE(xcldout(j2)%ice_pert(nmax))

      DO n=1, nmax
         
         !write to output
         CALL int2str(idx, n, '0')

         CALL new_channel_object(status, TRIM(sname), 'cov_pert_'//idx &
           , p3=xcldout(j2)%cov_pert(n)%ptr, reprid=GP_3D_MID)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, TRIM(sname), 'cov_pert_'//idx, &
           'long_name', &
           c='effective perturbation (nat+pert must not exceed 100%)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, TRIM(sname), 'cov_pert_'//idx, &
           'units', c='frac. ')
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, TRIM(sname),'tau_vis_pert_'//idx &
              , p3=xcldout(j2)%tau_vis_pert(n)%ptr &
              , reprid=GP_3D_MID)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, TRIM(sname), 'tau_vis_pert_'//idx, &
           'long_name', c='optical depth (visible) of cloud (perturbation)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, TRIM(sname), 'tau_vis_pert_'//idx, &
              'units', c='-- ')
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, TRIM(sname), 'ice_pert_'//idx &
           , p3=xcldout(j2)%ice_pert(n)%ptr, reprid=GP_3D_MID)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, TRIM(sname), 'ice_pert_'//idx, &
           'long_name', c='adjusted ice water mixing ratio ')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, TRIM(sname), 'ice_pert_'//idx, &
           'units', c='kg/kg')
         CALL channel_halt(substr, status)

      END DO

   END DO

   CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

 END SUBROUTINE cloudopt_init_memory
 ! ====================================================================

 ! ====================================================================
 SUBROUTINE cloudopt_init_coupling

   ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
   USE messy_main_channel_error_bi, ONLY: channel_halt
   USE messy_main_channel_bi,       ONLY: GP_3D_MID
   USE messy_main_channel,          ONLY: get_channel_object &
                                        , new_channel_object &
                                        , new_attribute
   USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
   USE messy_main_tools,            ONLY: strcrack, str2num  &
                                        , int2str

   IMPLICIT NONE

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER :: substr = 'cloudopt_init_coupling'
   INTEGER                                            :: status
   INTEGER                                            :: j1,j2
   INTEGER                                            :: dummy
   INTEGER                                            :: const_lev
   REAL(DP)                                           :: const_val
   REAL(DP)                                           :: const_vmr
   CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER:: outstring => NULL()    
   CHARACTER(LEN=STRLEN_MEDIUM)           :: sname = ''
   CHARACTER(LEN=2)                       :: call_idx = ''
   INTEGER                                :: n

   CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

   DO j2=1,NCALL

      CALL int2str(call_idx, j2, '0')
      sname=modstr//call_idx

      DO j1=1,NMAXOBJ

         DO n=1, xcldin(j1,j2)%n

            SELECT CASE(TRIM(xcldin(j1,j2)%name(n)%cha))
            CASE('#const')      
               !
               CALL strcrack(TRIM(xcldin(j1,j2)%name(n)%obj) &
                    , '=', outstring, dummy)
               CALL new_channel_object(status, TRIM(sname), outstring(1) &
                    , p3=xcldin(j1,j2)%ptr(n)%ptr, reprid=GP_3D_MID)
               CALL channel_halt(substr//': object '//outstring(1)//&
                    &' could not be created in channel '//&
                    &TRIM(sname)//' ', status)
               !
               CALL str2num(outstring(2), const_vmr, status)
               IF (status /= 0) &
                    CALL error_bi(TRIM(outstring(2))//&
                    &' could not be converted to real number ',substr)
               !
               xcldin(j1,j2)%ptr(n)%ptr = const_vmr
               !
               IF (ASSOCIATED(outstring)) THEN
                  DEALLOCATE(outstring)
                  NULLIFY(outstring)
               END IF
               !
            CASE('#const_cov')
               ! set constant value for one level
               !
               CALL strcrack(TRIM(xcldin(j1,j2)%name(n)%obj) &
                    , '=', outstring, dummy)
               CALL new_channel_object(status, TRIM(sname), outstring(1) &
                    , p3=xcldin(j1,j2)%ptr(n)%ptr, reprid=GP_3D_MID)
               CALL channel_halt(substr//': object '//outstring(1)//&
                    &' could not be created in channel '//&
                    &TRIM(sname)//' ', status)
               !
               CALL str2num(outstring(2), const_lev, status)
               IF (status /= 0) &
                    CALL error_bi(TRIM(outstring(2))//&
                    &' error when reading const_lev ',substr)
               !

               CALL str2num(outstring(3), const_val, status)
               IF (status /= 0) &
                    CALL error_bi(TRIM(outstring(3))//&
                    &' error when reading benchmark coverage ',substr)
               !
               xcldin(j1,j2)%ptr(n)%ptr(_RI_XYZ__(:,:,const_lev))= const_val

               IF (ASSOCIATED(outstring)) THEN
                  DEALLOCATE(outstring)
                  NULLIFY(outstring)
               END IF
               !
            CASE('#const_tau')
               !
               CALL strcrack(TRIM(xcldin(j1,j2)%name(n)%obj) &
                    , '=', outstring, dummy)
               CALL new_channel_object(status, TRIM(sname), outstring(1) &
                    , p3=xcldin(j1,j2)%ptr(n)%ptr, reprid=GP_3D_MID)
               CALL channel_halt(substr//': object '//outstring(1)//&
                    &' could not be created in channel '//&
                    &TRIM(sname)//' ', status)
               !
               CALL str2num(outstring(2), const_lev, status)
               IF (status /= 0) &
                    CALL error_bi(TRIM(outstring(2))//&
                    &' error when reading const lev ',substr)
               !
               CALL str2num(outstring(3), const_val, status)
               IF (status /= 0) &
                    CALL error_bi(TRIM(outstring(3))//&
                    &' error when reading benchmark tau ',substr)
               !
               xcldin(j1,j2)%ptr(n)%ptr(_RI_XYZ__(:,:,const_lev))= const_val

               IF (ASSOCIATED(outstring)) THEN
                  DEALLOCATE(outstring)
                  NULLIFY(outstring)
               END IF

            CASE('#std')
               !
               IF ((j1 /= id_radlp) .AND. (j1 /= id_radip)) THEN
                  CALL error_bi('#std only possible for radlp or radip',substr)
               ENDIF
               !
               CALL new_channel_object(status, TRIM(sname) &
                 , TRIM(xcldin(j1,j2)%name(n)%obj) &
                 , p3=xcldin(j1,j2)%ptr(n)%ptr, reprid=GP_3D_MID)
               CALL channel_halt(substr//': object '//&
                    & TRIM(xcldin(j1,j2)%name(n)%obj)//&
                    &' could not be created in channel '//&
                    &TRIM(sname)//' ', status)
               !
               IF (j1 == id_radlp) THEN
                  CALL new_attribute(status, TRIM(sname) &
                       , TRIM(xcldin(j1,j2)%name(n)%obj) &
                       , 'long_name', c='effective radius of liquid droplets')
                  CALL channel_halt(substr, status)
               ENDIF
               IF (j1 == id_radip) THEN
                  CALL new_attribute(status, TRIM(sname) &
                    , TRIM(xcldin(j1,j2)%name(n)%obj) &
                    , 'long_name', c='effective radius of ice particles')
                  CALL channel_halt(substr, status)
               ENDIF
               CALL new_attribute(status, TRIM(sname) &
                    , TRIM(xcldin(j1,j2)%name(n)%obj) &
                    , 'units', c='micrometer')
               CALL channel_halt(substr, status)
               !
            CASE DEFAULT
               !
               CALL get_channel_object(status &
                    , TRIM(xcldin(j1,j2)%name(n)%cha) &
                    , TRIM(xcldin(j1,j2)%name(n)%obj) &
                    , p3=xcldin(j1,j2)%ptr(n)%ptr)
               CALL channel_halt(substr//': object '//&
                    TRIM(xcldin(j1,j2)%name(n)%obj)//&
                    &' in channel '//&
                    TRIM(xcldin(j1,j2)%name(n)%cha)//' not found!',status)
               !
            END SELECT
         END DO
      END DO
   END DO

 END SUBROUTINE cloudopt_init_coupling
 ! ====================================================================

 ! ====================================================================
  SUBROUTINE cloudopt_radiation

    USE messy_main_grid_def_mem_bi,   ONLY: kproma,jrow, nlev
    USE messy_main_data_bi,           ONLY: tm1,apm1,aphm1,             &
                                            loland_2d,loglac_2d,lcouple
    USE messy_main_constants_mem,     ONLY: g
                                          
    IMPLICIT NONE
    INTRINSIC :: MAX, MIN, EPSILON

    REAL(DP), DIMENSION(nproma,nlev,NMAXOBJ) :: sum_m1
    REAL(DP), DIMENSION(nproma,nlev)         :: zdp
    REAL(DP), DIMENSION(nproma,nlev)         :: zlw
    REAL(DP), DIMENSION(nproma,nlev)         :: ziw
    REAL(DP), DIMENSION(nproma,nlev)         :: zcdnc
    REAL(DP), DIMENSION(kproma,nlev)         :: zradlp
    REAL(DP), DIMENSION(kproma,nlev)         :: zradip
    REAL(DP), DIMENSION(nproma,nlev)        :: zidxcld
    REAL(DP), DIMENSION(nproma,nlev,jpband) :: ztaucld_lw
    REAL(DP), DIMENSION(nproma,nlev,nsw)    :: ztaucld_sw
    REAL(DP), DIMENSION(nproma,nlev,nsw)    :: ztaucld_raw_sw
    REAL(DP), DIMENSION(nproma,nlev,nsw)    :: zomega_sw
    REAL(DP), DIMENSION(nproma,nlev,nsw)    :: zgamma_sw
    REAL(DP), DIMENSION(nproma,nlev)        :: zvalid_sw
    REAL(DP), DIMENSION(nproma,nlev)         :: zlwp 
    REAL(DP), DIMENSION(nproma,nlev)         :: ziwp 
    REAL(DP), DIMENSION(nproma)              :: zinhoml

    INTEGER :: j1, j2, n, n1, nmax, jb 
    INTEGER :: m_clc, m_lw, m_iw, m_cdnc, m_radlp, m_radip

    zdp(1:kproma,1:nlev) = aphm1(1:kproma,2:nlev+1)-aphm1(1:kproma,1:nlev)

    calls: DO j2=1,NCALL

       ! INIT
       xcldout(j2)%pidxcld(_RI_XYZ__(:,jrow,:))          = 0._dp
       xcldout(j2)%psumcov(_RI_XYZ__(:,jrow,:))          = 0._dp
       xcldout(j2)%pvalid_sw(_RI_XYZ__(:,jrow,:))        = 0._dp
       xcldout(j2)%ptaucld_lw(_RI_XYZN_(:,jrow,:,:))     = 0._dp
       xcldout(j2)%ptaucld_sw(_RI_XYZN_(:,jrow,:,:))     = 0._dp
       xcldout(j2)%ptaucld_raw_sw(_RI_XYZN_(:,jrow,:,:)) = 0._dp
       xcldout(j2)%pgamma_sw(_RI_XYZN_(:,jrow,:,:))      = 0._dp
       xcldout(j2)%pomega_sw(_RI_XYZN_(:,jrow,:,:))      = 0._dp
       sum_m1(:,:,:)                          = 0.0_dp

       ! effective cloud cover (natural + perturbation(s)), 
       ! input for rad, must not exceed 100%
       DO n=1, xcldin(id_clc,j2)%n
          xcldout(j2)%psumcov(_RI_XYZ__(:,jrow,:)) = &
             MIN( (xcldout(j2)%psumcov(_RI_XYZ__(:,jrow,:)) &
                  + xcldin(id_clc,j2)%ptr(n)%ptr(_RI_XYZ__(:,jrow,:)) ), 1.0_dp)
       END DO

       ! maximum-random-overlap-2d
       CALL calc_clcv(kproma, nproma, nlev &
            , xcldout(j2)%psumcov(_RI_XYZ__(:,jrow,:)) &
            , xcldout(j2)%pclcv(:,jrow))

       ! get max. number of perturbations of this call
       nmax = 1
       DO j1=1,NMAXOBJ
          nmax = MAX(nmax,xcldin(j1,j2)%n)
       END DO

       perturbs: DO n=1, nmax

          ! init for perturbations
          zidxcld(:,:)          = 0.0_dp
          zvalid_sw(:,:)        = 0.0_dp
          ztaucld_lw(:,:,:)     = 0.0_dp
          ztaucld_sw(:,:,:)     = 0.0_dp
          ztaucld_raw_sw(:,:,:) = 0.0_dp
          zgamma_sw(:,:,:)      = 0.0_dp
          zomega_sw(:,:,:)      = 0.0_dp
          ziw(:,:)              = 0.0_dp
          zlw(:,:)              = 0.0_dp
          xcldout(j2)%cov_pert(n)%ptr(_RI_XYZ__(:,jrow,:))     = 0._dp
          xcldout(j2)%tau_vis_pert(n)%ptr(_RI_XYZ__(:,jrow,:)) = 0._dp
          xcldout(j2)%ice_pert(n)%ptr(_RI_XYZ__(:,jrow,:))     = 0._dp

          m_clc   = MIN(n, xcldin(id_clc,j2)%n) 
          m_lw    = MIN(n, xcldin(id_lw,j2)%n)
          m_iw    = MIN(n, xcldin(id_iw,j2)%n)
          m_cdnc  = MIN(n, xcldin(id_cdnc,j2)%n)
          m_radlp = MIN(n, xcldin(id_radlp,j2)%n)
          m_radip = MIN(n, xcldin(id_radip,j2)%n)

          ! sum of natural and perturbation cover must not exceed 100%, 
          ! (sum_m1 is an auxiliary variable);
          ! adjust perturbation(s) successively, if necessary 
          ! => cov_pert(n) is the adjusted perturbation cover
          sum_m1(:,:,1) = 0._dp

          DO n1=2,m_clc
             sum_m1(:,:,n1) = MIN( &
                  ( sum_m1(:,:,n1-1) + &
                   xcldin(id_clc,j2)%ptr(n1-1)%ptr(_RI_XYZ__(:,jrow,:)) ), 1.0_dp)
          END DO
          xcldout(j2)%cov_pert(m_clc)%ptr(_RI_XYZ__(:,jrow,:)) = &
               MAX((MIN((sum_m1(:,:,m_clc) + &
               xcldin(id_clc,j2)%ptr(m_clc)%ptr(_RI_XYZ__(:,jrow,:))),1.0_dp) &
               - sum_m1(:,:,m_clc)),0._dp)

          zlw(:,:)   = xcldin(id_lw,j2)%ptr(m_lw)%ptr(_RI_XYZ__(:,jrow,:))

          ! 1. step: back calculation so that constant ice water mixing 
          !          ratio from namelist gives constant tau
          ! 2. step: if perturbation cover has to be adjusted because of 100%
          !          threshold, ziw has to be adjusted subsequently, 
          !          otherwise huge tau for tiny coverages results 
          IF (TRIM(xcldin(id_iw,j2)%name(n)%cha) == '#const_tau') THEN
             WHERE(xcldin(id_clc,j2)%ptr(m_clc)%ptr(_RI_XYZ__(:,jrow,:)) > 0.0_dp)
                ziw(:,:) = (xcldin(id_iw,j2)%ptr(m_iw)%ptr(_RI_XYZ__(1:kproma,jrow,:))* &
                     g*xcldout(j2)%cov_pert(n)%ptr(_RI_XYZ__(:,jrow,:))) / &
                     (1000._dp*zdp(1:kproma,:))* &
                     (xcldout(j2)%cov_pert(n)%ptr(_RI_XYZ__(:,jrow,:))/ &
                     xcldin(id_clc,j2)%ptr(m_clc)%ptr(_RI_XYZ__(:,jrow,:)))
             ELSE WHERE
                ziw(:,:) = 0.0_dp
             END WHERE
          ELSE
             ziw(:,:) = xcldin(id_iw,j2)%ptr(m_iw)%ptr(_RI_XYZ__(:,jrow,:))
          END IF

          xcldout(j2)%ice_pert(n)%ptr(_RI_XYZ__(:,jrow,:)) = ziw(:,:)

          zcdnc(:,:) = xcldin(id_cdnc,j2)%ptr(m_cdnc)%ptr(_RI_XYZ__(:,jrow,:))

          IF ( (TRIM(xcldin(id_radlp,j2)%name(m_radlp)%cha) == '#std') .OR. &
               (TRIM(xcldin(id_radip,j2)%name(m_radip)%cha) == '#std') ) THEN

             CALL rad_std( &
                  apm1(1:kproma,:)                                  & ! IN
                  , aphm1(1:kproma,2:nlev+1)-aphm1(1:kproma,1:nlev) & ! IN
                  , tm1(_RI_XYZ__(1:kproma,jrow,:))                 & ! IN
                  , ziw(1:kproma,:)                                 & ! IN
                  , zlw(1:kproma,:)                                 & ! IN
                  , xcldout(j2)%cov_pert(n)%ptr(_RI_XYZ__(1:kproma,jrow,:)) & ! IN
                  , zcdnc(1:kproma,:)                               & ! IN
                  , SPREAD(loland_2d(1:kproma,jrow), 2, nlev)       & ! IN
                  , SPREAD(loglac_2d(1:kproma,jrow), 2, nlev)       & ! IN
                  , zradip(:,:), zradlp(:,:)                        & ! OUT 
                  )

             IF (TRIM(xcldin(id_radlp,j2)%name(m_radlp)%cha) == '#std') THEN
                xcldin(id_radlp,j2)%ptr(m_radlp)%ptr(_RI_XYZ__(1:kproma,jrow,:)) &

                     = zradlp(1:kproma,:)
             END IF

             IF (TRIM(xcldin(id_radip,j2)%name(m_radip)%cha) == '#std') THEN
                xcldin(id_radip,j2)%ptr(m_radip)%ptr(_RI_XYZ__(1:kproma,jrow,:)) &
                     = zradip(1:kproma,:)
             END IF
          END IF

          ! new call preswlw
          ! this will calculate different variables which are needed
          ! for the actual calculation of the cloud optical properties
          ! in the lw and sw range
           CALL cloud_opt_prelwsw(                                  &
                kproma, nproma, nlev,                           & ! IN
                xcldout(j2)%cov_pert(n)%ptr(_RI_XYZ__(:,jrow,:)), &
                xcldout(j2)%psumcov(_RI_XYZ__(:,jrow,:)),         &
                zlw,ziw,zdp,                                    &
                lcouple,                                        &
                zidxcld(:,:),                                   & !OUT
                zlwp, ziwp,                                     &
                zinhoml)                                         !OUT

           ! new call lw1
           CALL cloud_opt_lw1(                                 &
                kproma, nproma, nlev,                           & !IN
                zlwp,ziwp,                                      &
                xcldin(id_radlp,j2)%ptr(m_radlp)%ptr(_RI_XYZ__(:,jrow,:)), &
                xcldin(id_radip,j2)%ptr(m_radip)%ptr(_RI_XYZ__(:,jrow,:)), & !IN
                zinhoml,                                       & !IN
                ztaucld_lw(:,:,:))                                !OUT

           ! new call sw1
           CALL cloud_opt_sw1(                                  &
                kproma, nproma, nlev,                           &  !IN
                zidxcld(:,:),                                   &
                zlwp,ziwp,                                      &
                xcldin(id_radlp,j2)%ptr(m_radlp)%ptr(_RI_XYZ__(:,jrow,:)), &
                xcldin(id_radip,j2)%ptr(m_radip)%ptr(_RI_XYZ__(:,jrow,:)), & !IN
                zinhoml,                                        &!IN
                ztaucld_sw(:,:,:),                              & !OUT
                ztaucld_raw_sw(:,:,:),                          &
                zomega_sw(:,:,:),                               &
                zgamma_sw(:,:,:),                               &
                zvalid_sw(:,:))                                   !OUT

          ! save tau (visible) of cloud or perturbation for output
          xcldout(j2)%tau_vis_pert(n)%ptr(_RI_XYZ__(:,jrow,:)) = ztaucld_raw_sw(:,:,1)

          ! aggregated clear/cloudy index
          xcldout(j2)%pidxcld(_RI_XYZ__(:,jrow,:)) = &
               MAX(xcldout(j2)%pidxcld(_RI_XYZ__(:,jrow,:)), zidxcld(:,:))

          ! aggregated validity flag
          xcldout(j2)%pvalid_sw(_RI_XYZ__(:,jrow,:)) = &
               MAX(xcldout(j2)%pvalid_sw(_RI_XYZ__(:,jrow,:)), zvalid_sw(:,:))

          ! aggregated tau (additiv)
          xcldout(j2)%ptaucld_lw(_RI_XYZN_(:,jrow,:,:)) =  &
               xcldout(j2)%ptaucld_lw(_RI_XYZN_(:,jrow,:,:)) + &
               ztaucld_lw(:,:,:)

          xcldout(j2)%ptaucld_sw(_RI_XYZN_(:,jrow,:,:)) =  &
               xcldout(j2)%ptaucld_sw(_RI_XYZN_(:,jrow,:,:)) + &
               ztaucld_sw(:,:,:)

          xcldout(j2)%ptaucld_raw_sw(_RI_XYZN_(:,jrow,:,:))= &
               xcldout(j2)%ptaucld_raw_sw(_RI_XYZN_(:,jrow,:,:)) + &
               ztaucld_raw_sw(:,:,:)

          ! aggregated gamma
          ! Note: - division by zomgmx (=omega_sw) uncommented in SMCL
          xcldout(j2)%pgamma_sw(_RI_XYZN_(:,jrow,:,:)) = &
               xcldout(j2)%pgamma_sw(_RI_XYZN_(:,jrow,:,:)) + zgamma_sw(:,:,:)

          ! aggregated omega
          ! Note: - division by ztaumx (=tau_sw_raw) uncommented in SMCL
          xcldout(j2)%pomega_sw(_RI_XYZN_(:,jrow,:,:)) = &
               xcldout(j2)%pomega_sw(_RI_XYZN_(:,jrow,:,:)) + zomega_sw(:,:,:)
       END DO perturbs

       ! SEE uncommented parts in SMCL
       DO jb=1, nsw
          WHERE(xcldout(j2)%pvalid_sw(_RI_XYZ__(:,jrow,:)) > 0.0_dp)
             ! first use omega to calculate gamma,
             ! then modify omega by tau_raw
             ! tau_raw is an auxiliary variable without
             ! cloud inhomogeneity factor
             xcldout(j2)%pgamma_sw(_RI_XYZN_(:,jrow,:,jb)) = &
                  xcldout(j2)%pgamma_sw(_RI_XYZN_(:,jrow,:,jb)) / &
                  xcldout(j2)%pomega_sw(_RI_XYZN_(:,jrow,:,jb))

             xcldout(j2)%pomega_sw(_RI_XYZN_(:,jrow,:,jb)) = &
                  xcldout(j2)%pomega_sw(_RI_XYZN_(:,jrow,:,jb)) / &
                  xcldout(j2)%ptaucld_raw_sw(_RI_XYZN_(:,jrow,:,jb))
          ELSE WHERE
             xcldout(j2)%pgamma_sw(_RI_XYZN_(:,jrow,:,jb)) = 0.0_dp
             xcldout(j2)%pomega_sw(_RI_XYZN_(:,jrow,:,jb)) = 1.0_dp
          END WHERE
       END DO

       ! diagnostics for ISCCP simulator
       ! band 3 is 440 - 690 nm, needed is 670 nm
       xcldout(j2)%isccp_cldtau3d(_RI_XYZ__(:,jrow,:)) = &
            xcldout(j2)%ptaucld_sw(_RI_XYZN_(:,jrow,:,3))
       ! band 6 is 820 - 980 cm-1, we need 10.5 µm channel
       xcldout(j2)%isccp_cldemi3d(_RI_XYZ__(:,jrow,:)) = 1._dp - &
            exp(-1._dp*xcldout(j2)%ptaucld_lw(_RI_XYZN_(:,jrow,:,6)))
       xcldout(j2)%isccp_f3d(_RI_XYZ__(:,jrow,:)) = &
            MAX(xcldout(j2)%psumcov(_RI_XYZ__(:,jrow,:)),EPSILON(1.0_dp))

    END DO calls

  END SUBROUTINE cloudopt_radiation
  ! ====================================================================

  !=====================================================================
  SUBROUTINE cloudopt_free_memory

    IMPLICIT NONE
    INTEGER :: j1, j2, n, nmax

    DO j2=1,NCALL
       DO j1=1,NMAXOBJ
          nmax = xcldin(j1,j2)%n
          DO n=1, nmax
             NULLIFY(xcldin(j1,j2)%ptr(n)%ptr)
          END DO
          DEALLOCATE(xcldin(j1,j2)%ptr)
          NULLIFY(xcldin(j1,j2)%ptr)
          DEALLOCATE(xcldin(j1,j2)%name)
          NULLIFY(xcldin(j1,j2)%name)
       END DO
    END DO
  
    DO j2=1, NCALL
       IF (ASSOCIATED(ptaucld_lw_5d(j2)%ptr))   THEN
          DEALLOCATE(ptaucld_lw_5d(j2)%ptr)
          NULLIFY(ptaucld_lw_5d(j2)%ptr)
       END IF

       IF (ASSOCIATED(ptaucld_sw_5d(j2)%ptr)) THEN
          DEALLOCATE(ptaucld_sw_5d(j2)%ptr)
          NULLIFY(ptaucld_sw_5d(j2)%ptr)
       END IF

       IF (ASSOCIATED(ptaucld_raw_sw_5d(j2)%ptr)) THEN
          DEALLOCATE(ptaucld_raw_sw_5d(j2)%ptr)
          NULLIFY(ptaucld_raw_sw_5d(j2)%ptr)
       END IF

       IF (ASSOCIATED(pomega_sw_5d(j2)%ptr)) THEN
          DEALLOCATE(pomega_sw_5d(j2)%ptr)
          NULLIFY(pomega_sw_5d(j2)%ptr)
       END IF

       IF (ASSOCIATED(pgamma_sw_5d(j2)%ptr)) THEN
          DEALLOCATE(pgamma_sw_5d(j2)%ptr)
          NULLIFY(pgamma_sw_5d(j2)%ptr)
       END IF
      
       IF (ASSOCIATED(xcldout(j2)%cov_pert)) THEN
          DEALLOCATE(xcldout(j2)%cov_pert)
          NULLIFY(xcldout(j2)%cov_pert)
       END IF

       IF (ASSOCIATED(xcldout(j2)%tau_vis_pert)) THEN
          DEALLOCATE(xcldout(j2)%tau_vis_pert)
          NULLIFY(xcldout(j2)%tau_vis_pert)
       END IF

       IF (ASSOCIATED(xcldout(j2)%ice_pert)) THEN
          DEALLOCATE(xcldout(j2)%ice_pert)
          NULLIFY(xcldout(j2)%ice_pert)
       END IF

    END DO

    DEALLOCATE(ptaucld_lw_5d) ; NULLIFY(ptaucld_lw_5d)
    DEALLOCATE(ptaucld_sw_5d) ; NULLIFY(ptaucld_sw_5d)
    DEALLOCATE(ptaucld_raw_sw_5d) ; NULLIFY(ptaucld_raw_sw_5d)
    DEALLOCATE(pomega_sw_5d) ; NULLIFY(pomega_sw_5d)
    DEALLOCATE(pgamma_sw_5d) ; NULLIFY(pgamma_sw_5d)

  END SUBROUTINE cloudopt_free_memory
  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE cldopt_read_nml_cpl(status, iou)

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ NCALL, cld_inp            


    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='cldopt_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE cldopt_read_nml_cpl
  ! ====================================================================

#endif
! **********************************************************************
END MODULE messy_cloudopt_si
! **********************************************************************
