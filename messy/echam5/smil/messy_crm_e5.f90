
! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL CRM 
!
! Author : Harald Rybka, IPA-Mainz, August  2013
!
! Currently only suitable for ECHAM/MESSy.
! References: see messy_crm.f90
!
! **********************************************************************

! **********************************************************************
MODULE messy_crm_e5
! **********************************************************************

  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: dp, cp_air, alv, g,                 &
                                      rhoh2o=>rho_H2O

  USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast,       &
                                      finish, message
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi,   &
                                      error_bi, warning_bi
  USE messy_main_timer,         ONLY: time_step_len,        delta_time,   &
                                      nstep => current_time_step

  USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, kproma, jrow,       &
                                        ngpblks
  USE messy_main_data_bi,       ONLY: tm1, qm1, xlm1, xim1,               &
                                      um1, vm1,                           &
                                      aps, press_3d, pressi_3d,           &
                                      geopoti => geopoti_3d,              &
                                      geopot => geopot_3d,                & 
                                      tte_3d, qte_3d, xlte_3d, xite_3d,   &
                                      vom_3d, vol_3d

#ifdef MESSYTENDENCY
  ! tendency budget
    USE messy_main_tendency_bi,    ONLY: mtend_get_handle,               &
                                         mtend_get_start_l,              &
                                         mtend_add_l,                    &
                                         mtend_register,                 &
                                         mtend_id_t,                     &
                                         mtend_id_q,                     &
                                         mtend_id_xl,                    &
                                         mtend_id_xi,                    &
                                         mtend_id_u,                     &
                                         mtend_id_v,                     &
                                         mtend_id_tracer
#endif

  ! SMCL
  USE messy_crm

  IMPLICIT NONE
  SAVE

  INTRINSIC :: NULL, INT, MAX, MIN

  PRIVATE

  ! GLOBAL PARAMETERS
  LOGICAL, PARAMETER :: CRM_CALC = .TRUE.
  INTEGER            :: crm_nx, crm_ny, crm_nz


! list of grid box mean values after n CRM timesteps = one GCM time step

  REAL(dp), POINTER, DIMENSION(:,:,:) :: tl         => NULL()   ! global grid temperature (K)
  REAL(dp), POINTER, DIMENSION(:,:,:) :: ql         => NULL()   ! global grid water vapour (g/g)
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qccl       => NULL()   ! global grid cloud liquid water (g/g)
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qiil       => NULL()   ! global grid cloud ice (g/g)
! cloud type channel objects (column integrated)
  REAL(dp), POINTER, DIMENSION(:,:)   :: cltot      => NULL()   ! shaded cloud fraction
  REAL(dp), POINTER, DIMENSION(:,:)   :: clhgh      => NULL()   ! 
  REAL(dp), POINTER, DIMENSION(:,:)   :: clmed      => NULL()   !
  REAL(dp), POINTER, DIMENSION(:,:)   :: cllow      => NULL()   !

  REAL(dp), POINTER, DIMENSION(:,:,:) :: ultend     => NULL()   ! tendency of ul
  REAL(dp), POINTER, DIMENSION(:,:,:) :: vltend     => NULL()   ! tendency of vl
  REAL(dp), POINTER, DIMENSION(:,:,:) :: ttend      => NULL()   ! tendency of temperature (K/s) (t_ls)
  REAL(dp), POINTER, DIMENSION(:,:,:) :: sltend     => NULL()   ! tendency of static energy
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qltend     => NULL()   ! tendency of water vapour
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qcltend    => NULL()   ! tendency of cloud liquid water vapour
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qiltend    => NULL()   ! tendency of cloud ice
! precipitation channel objects
  REAL(dp), POINTER, DIMENSION(:,:)   :: precc      => NULL()   ! convective precipitation rate (m/s)
  REAL(dp), POINTER, DIMENSION(:,:)   :: precl      => NULL()   ! stratiform precipitation rate (m/s)
  REAL(dp), POINTER, DIMENSION(:,:)   :: precsc     => NULL()   ! convective snow rate
  REAL(dp), POINTER, DIMENSION(:,:)   :: precsl     => NULL()   ! stratiform snow rate
! cloud distribution channel objects
  REAL(dp), POINTER, DIMENSION(:,:,:) :: cld        => NULL()   ! cloud fraction
  REAL(dp), POINTER, DIMENSION(:,:,:) :: cldtop     => NULL()   ! cloud top pdf
  REAL(dp), POINTER, DIMENSION(:,:,:) :: cldbot     => NULL()   ! cloud bot pdf
  REAL(dp), POINTER, DIMENSION(:,:,:) :: convtop    => NULL()   ! convective top height pdf
  REAL(dp), POINTER, DIMENSION(:,:,:) :: convbot    => NULL()   ! convective bottom height pdf
  REAL(dp), POINTER, DIMENSION(:,:)   :: conv_cth   => NULL()   ! convective top level
  REAL(dp), POINTER, DIMENSION(:,:)   :: conv_cbh   => NULL()   ! convective bottom level
! channel objects (not yet integrated)
  REAL(dp), POINTER, DIMENSION(:,:,:) :: gicewp     => NULL()   ! ice water path
  REAL(dp), POINTER, DIMENSION(:,:,:) :: gliqwp     => NULL()   ! liquid water path

  REAL(dp), POINTER, DIMENSION(:,:,:) :: mc         => NULL()   ! cloud mass flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: mcup       => NULL()   ! updraft cloud mass flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: mcdn       => NULL()   ! downdraft cloud mass flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: mcuup      => NULL()   ! unsaturated updraft cloud mass flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: mcudn      => NULL()   ! unsaturated downdraft cloud mass flux

  REAL(dp), POINTER, DIMENSION(:,:,:) :: massfu      => NULL()   ! convective updraft cloud mass flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: massfd      => NULL()   ! convective downdraft cloud mass flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: u_entr      => NULL()   ! convective updraft entrainment
  REAL(dp), POINTER, DIMENSION(:,:,:) :: u_detr      => NULL()   ! convective updraft detrainment
  REAL(dp), POINTER, DIMENSION(:,:,:) :: d_entr      => NULL()   ! convective downdraft entrainment
  REAL(dp), POINTER, DIMENSION(:,:,:) :: d_detr      => NULL()   ! convective downdraft detrainment

! channel objects for cloud "water" amounts 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: crm_qc     => NULL()   ! cloud water
  REAL(dp), POINTER, DIMENSION(:,:,:) :: crm_qi     => NULL()   ! cloud ice
  REAL(dp), POINTER, DIMENSION(:,:,:) :: crm_qs     => NULL()   ! cloud snow
  REAL(dp), POINTER, DIMENSION(:,:,:) :: crm_qg     => NULL()   ! cloud graupel
  REAL(dp), POINTER, DIMENSION(:,:,:) :: crm_qr     => NULL()   ! cloud rain

  REAL(dp), POINTER, DIMENSION(:,:,:) :: tkez       => NULL()   ! TKE profile
  REAL(dp), POINTER, DIMENSION(:,:,:) :: tkesgsz    => NULL()   ! sub-grid scale (SGS) TKE profile

  REAL(dp), POINTER, DIMENSION(:,:,:) :: flux_u     => NULL()   ! zonal momentum flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: flux_v     => NULL()   ! meridional momentum flux

  REAL(dp), POINTER, DIMENSION(:,:,:) :: flux_qp    => NULL()   ! precipitating water flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: flux_qt    => NULL()   ! non-precipitating water flux
  REAL(dp), POINTER, DIMENSION(:,:,:) :: fluxsgs_qt => NULL()   ! SGS non-precipitating water flux

  REAL(dp), POINTER, DIMENSION(:,:,:) :: pflx       => NULL()   ! precipitation flux

  REAL(dp), POINTER, DIMENSION(:,:,:) :: qt_ls      => NULL()   ! tendency of non-precipitating water due to large scale
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qt_trans   => NULL()   ! tendency of non-precipitating water due to transport
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qp_trans   => NULL()   ! tendency of precipitating water due to transport
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qp_fall    => NULL()   ! tendency of precipitating water due to fall-out
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qp_src     => NULL()   ! tendency of precipitating water due to conversion
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qp_evp     => NULL()   ! tendency of precipitating water due to evaporation

  REAL(dp), POINTER, DIMENSION(:,:)   :: prectend   => NULL()   ! column integrated tendency in precipitation water+ice (kg/m2/s)
  REAL(dp), POINTER, DIMENSION(:,:)   :: precstend  => NULL()   ! column integrated tendency in precipitation ice (kg/m2/s)
  REAL(dp), POINTER, DIMENSION(:,:)   :: taux_crm   => NULL()   ! zonal CRM surface stress perturbation (N/m2)
  REAL(dp), POINTER, DIMENSION(:,:)   :: tauy_crm   => NULL()   ! meridional CRM surface stress perturbation (N/m2)
  REAL(dp), POINTER, DIMENSION(:,:)   :: z0m        => NULL()   ! surface stress (N/m2)


  REAL(dp), POINTER, DIMENSION(:,:)   :: ustrl     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: ustrw     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: ustri     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: vstrl     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: vstrw     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: vstri     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: wind10_2d => NULL()

  REAL(dp), POINTER, DIMENSION(:,:)   :: cv_cover    => NULL()  ! convective cloud cover of CRM cells 
  REAL(dp), POINTER, DIMENSION(:,:)   :: time_factor => NULL()  ! timing factor due to local unstable flow in CRM
! channel objects for CAPE calculation
  REAL(dp), POINTER, DIMENSION(:,:)   :: CAPE_c        => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: CAPE_s        => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: PA_s          => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: PA_c          => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: CAPE_MAX_IDX  => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: PLCL          => NULL()

! CRM grid variables
!  REAL(dp), POINTER, DIMENSION(jl,crm_nx,crm_ny,crm_nz,crmvars)   :: crm_buffer
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:)  :: crm_buffer => NULL()

  REAL(dp), POINTER, DIMENSION(:,:,:,:,:)    :: t_rad      => NULL() ! rad. temperature
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:)    :: qv_rad     => NULL() ! rad. water vapour
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:)    :: qc_rad     => NULL() ! rad. cloud water
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:)    :: qi_rad     => NULL() ! rad. cloud ice
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:)    :: qrad_crm   => NULL()

! Auxiliary CRM grid variables
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:,:)  :: crm_buffer1 => NULL()

  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:)    :: t_rad1      => NULL() ! rad. temperature
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:)    :: qv_rad1     => NULL() ! rad. water vapour
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:)    :: qc_rad1     => NULL() ! rad. cloud water
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:)    :: qi_rad1     => NULL() ! rad. cloud ice
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:)    :: qrad_crm1   => NULL()

  ! um_hr_20151016+
  ! CRM subgrid-channel objects for coupling purposes:
  INTEGER,PARAMETER :: ido_tm1_sg   = 1
  INTEGER,PARAMETER :: ido_qm1_sg   = 2
  INTEGER,PARAMETER :: ido_xlm1_sg  = 3
  INTEGER,PARAMETER :: ido_xim1_sg  = 4
  INTEGER,PARAMETER :: ido_aclc_sg  = 5
  INTEGER,PARAMETER :: ido_acdnc_sg = 6
  INTEGER,PARAMETER :: ido_radlp_sg = 7
  INTEGER,PARAMETER :: ido_radip_sg = 8
  INTEGER,PARAMETER :: ido_prec_sg  = 9
  INTEGER,PARAMETER :: ido_qpc_sg   = 10
  INTEGER,PARAMETER :: ido_qpi_sg   = 11

  INTEGER,PARAMETER :: NOUTSG = 11

  CHARACTER(LEN=12), DIMENSION(NOUTSG):: setname_sg=(/ &
       'tm1_sg      ','qm1_sg      ','xlm1_sg     ','xim1_sg     ', &
       'aclc_sg     ','acdnc_sg    ','radlp_sg    ','radip_sg    ', &
       'prec_sg     ','qpc_sg      ','qpi_sg      '/)

  CHARACTER(LEN=16), DIMENSION(NOUTSG):: repr_name_sg=(/ &
       'GP_3D_MID       ','GP_3D_MID       ','GP_3D_MID       ','GP_3D_MID       ', &
       'GP_3D_MID       ','GP_3D_MID       ','GP_3D_MID       ','GP_3D_MID       ', &
       'GP_2D_HORIZONTAL','GP_3D_MID       ','GP_3D_MID       '/)
  INTEGER, DIMENSION(NOUTSG):: repr_idx_sg

    CHARACTER(LEN=39), DIMENSION(NOUTSG):: longname_sg=(/& 
    'CRM temperature                        ', 'CRM water vapour                       ', &
    'CRM cloud water                        ', 'CRM cloud ice                          ', &
    'CRM cloud cover                        ', 'CRM cloud condensation nuclei          ', &
    'CRM effective radii for liquid droplets', 'CRM effective radii for ice droplets   ', &
    'CRM surface precipitation              ', 'CRM precipitation water                ', &
    'CRM precipitation ice                  '/)

  CHARACTER(LEN=10), DIMENSION(NOUTSG):: unit_sg=(/&
     'K         ', 'kg/kg     ', 'kg/kg     ', 'kg/kg     ', &
     '-         ', '1/m**3    ', 'micrometer', 'micrometer', &
     'kg/kg     ', 'kg/kg     ', 'kg/kg     '/)

  TYPE(t_crm_work), DIMENSION(NOUTSG)     :: xcrmoutsg
  ! um_hr_20151016-

! channel objects for CRM TRACER TRANSPORT
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:,:)  :: crm_tracer1  => NULL()  ! tracer variable for transport subroutine in submodel CRM
  REAL(dp), POINTER, DIMENSION(:,:,:,:,:,:)    :: crm_tracer   => NULL()  ! tracer variable for transport subroutine in submodel CRM

#ifdef MESSYTENDENCY
  ! variable for tendency budget
  integer :: my_handle
#endif

  PUBLIC :: crm_initialize
  PUBLIC :: crm_init_memory, crm_free_memory
  PUBLIC :: crm_init_coupling
  PUBLIC :: crm_radiation

CONTAINS  

! ####################################################################
! PUBLIC SUBROUTINES
! ####################################################################

! ====================================================================
  SUBROUTINE crm_initialize

! read and check CRM namelists
  USE messy_main_tools,         ONLY: find_next_free_unit                                  
  
    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------                    
       
  IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'crm_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    LOGICAL  :: dotracer ! subgrid tracer transport switch
    REAL(dp) :: tau_min, tau_max, damp_depth

    EXTERNAL :: set_crm_param, set_crm_micro_SAM1MOM, set_crm_micro_M2005

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CTRL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find free I/O unit
       CALL crm_read_nml_ctrl(status, iou)   ! read CTRL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
    END IF

    ! BROADCAST CTRL namelist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(ngridcells_x, p_io)
    CALL p_bcast(ngridcells_y, p_io)
    CALL p_bcast(crm_top, p_io)
    CALL p_bcast(crm_size, p_io)
    CALL p_bcast(crm_timestep, p_io)
    CALL p_bcast(CRM_3D, p_io)
    CALL p_bcast(crm_orient, p_io)
    CALL p_bcast(dosgtracer, p_io)

    CALL p_bcast(micro_scheme, p_io)
    CALL p_bcast(crmvars, p_io)
    CALL p_bcast(nmicro_fields, p_io)

    CALL p_bcast(docloud, p_io)
    CALL p_bcast(doprecip, p_io)
    CALL p_bcast(dosgs, p_io)
    CALL p_bcast(dosmagor, p_io)

    CALL p_bcast(dodamping, p_io)
    CALL p_bcast(set_damp(:),p_io)

    CALL p_bcast(dosurface, p_io)
    CALL p_bcast(dosfc_flx_fxd, p_io)
    CALL p_bcast(dosfc_tau_fxd, p_io)

    CALL p_bcast(dowallx, p_io)
    CALL p_bcast(dowally, p_io)
    CALL p_bcast(docoriolis, p_io)
    CALL p_bcast(docolumn, p_io)

    ! parameter settings for one-mom microphyiscs
    CALL p_bcast(qcw0, p_io)
    CALL p_bcast(qci0, p_io)
    CALL p_bcast(alphaelq, p_io)
    CALL p_bcast(betaelq, p_io)
    CALL p_bcast(qp_threshold, p_io)

    ! parameter settings for two-mom microphyiscs
    CALL p_bcast(doicemicro, p_io)
    CALL p_bcast(dograupel, p_io)
    CALL p_bcast(dohail, p_io)
    CALL p_bcast(dosb_warm_rain, p_io)
    CALL p_bcast(dopredictNc, p_io)
    CALL p_bcast(dospecifyaerosol, p_io)
    CALL p_bcast(dosubgridw, p_io)
    CALL p_bcast(doarcticicenucl, p_io)
    CALL p_bcast(docloudedgeact, p_io)

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register (my_handle,mtend_id_t)
    CALL mtend_register (my_handle,mtend_id_q)
    CALL mtend_register (my_handle,mtend_id_tracer)
    CALL mtend_register (my_handle,mtend_id_xl)
    CALL mtend_register (my_handle,mtend_id_xi)
    
    CALL mtend_register (my_handle,mtend_id_u)
    CALL mtend_register (my_handle,mtend_id_v)
#endif

    tau_min    = set_damp(1)
    tau_max    = set_damp(2)
    damp_depth = set_damp(3)

    IF (dosgtracer == 0) THEN 
       dotracer = .FALSE.
    ELSE
       dotracer = .TRUE.
    END IF

    CALL set_crm_param(docloud, doprecip, dosgs, dosmagor,   &
         dodamping, tau_min, tau_max, damp_depth,            &
         dosurface, dosfc_flx_fxd, dosfc_tau_fxd,            &
         dowallx, dowally, docoriolis, docolumn, dotracer)

    IF (micro_scheme.eq.0) THEN
       CALL set_crm_micro_SAM1MOM(nmicro_fields, qcw0, qci0, alphaelq, betaelq, qp_threshold)
    ELSEIF (micro_scheme.eq.1) THEN
       CALL set_crm_micro_M2005(nmicro_fields, doicemicro, dograupel, dohail,   &
                                dopredictNc, dospecifyaerosol, dosubgridw,      &
                                dosb_warm_rain, doarcticicenucl, docloudedgeact )
    END IF

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output
    
  END SUBROUTINE crm_initialize

  ! ====================================================================

  SUBROUTINE crm_init_memory

    USE messy_main_channel,            ONLY: new_channel, new_channel_object &
                                           , new_attribute                   &
                                           , new_channel_object_reference    &
                                           , REPR_UNDEF

    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: GP_3D_MID, DC_GP  &
                                           , GP_2D_HORIZONTAL, SCALAR        &
                                           , DIMID_LON, DIMID_LAT, DC_BC     &
                                           , gp_nseg, gp_start, gp_cnt       &
                                           , gp_meml, gp_memu
  
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO        &
                                           , set_representation_decomp       &
                                           , get_representation_info         &
                                           , get_representation_id           &    
                                           , IRANK, PIOTYPE_COL
 
    USE messy_main_grid_def_mem_bi,    ONLY: nlevp1, nvclev, vct

    USE messy_main_tracer_mem_bi,   ONLY: ntrac => ntrac_gp

    USE messy_main_tools,           ONLY: int2str 

  IMPLICIT NONE
  CHARACTER(LEN=*), PARAMETER::substr='crm_init_memory'
  INTEGER :: status

  INTEGER            :: nsubsteps, nsub, ngridcells
  INTEGER, PARAMETER :: VAL_OUT = 5
  CHARACTER(LEN=3)   :: istr, jstr, tstr
  CHARACTER(LEN=1)   :: kstr
  INTEGER            :: i,j,k,jt,jk,nx,ny

  !um_hr_20190307+
  REAL(dp) :: h_a(nvclev), h_b(nvclev), sfpress
  REAL(dp) :: hypi(nlevp1)
  INTEGER  :: max_lev_crm
  !um_hr_20190307+

  CHARACTER(LEN=*), PARAMETER  :: crmsg=modstr//'_sg'
  CHARACTER(LEN=*), PARAMETER  :: crmsgb=modstr//'_sg_buf'
  CHARACTER(LEN=*), PARAMETER  :: crmsgr=modstr//'_sg_rad'
  CHARACTER(LEN=*), PARAMETER  :: crmsgt=modstr//'_sg_trac'

  INTEGER                          :: nseg = 0
  INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
  INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
  INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
  INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

  INTEGER, DIMENSION(NOUTSG)       :: irank_sg

  REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  mem => NULL()

  INTEGER           :: DIMID_CRMLEV
  INTEGER           :: REPR_CRM_3D_MID
  INTEGER           :: DIMID_CRM1LEV
  INTEGER           :: REPR_CRM_3D_1LEV

  CHARACTER(LEN=4)  :: crmbufname
  CHARACTER(LEN=5)  :: crmbufunit
  CHARACTER(LEN=40) :: crmlongname

  EXTERNAL :: set_crm_dims
  EXTERNAL :: set_crm_tracers
  EXTERNAL :: crm_allocate_vars, crm_allocate_microphysics, crm_allocate_tracers

  ! ====================================================================

  !um_hr_20190307+
  ! defining top level of CRM box 
  ! calculating number of CRM levels
  ! and setting up CRM dimensions
  IF (p_parallel_io) &
       PRINT*, substr, ': restrict CRM to a height below ', &
       crm_top/100, ' hPa!'
  DO jk=1,nvclev
     h_a(jk) = vct(jk)
     h_b(jk) = vct(jk+nvclev)
  END DO
  sfpress = 1.e5_dp ! reference pressure of 1000 hPa
  DO jk=1,nlev+1
     hypi(jk)      = h_a(jk) + h_b(jk) * sfpress
  ENDDO
  max_lev_crm = 1
  DO jk=1,nlev
     IF (hypi(jk) < crm_top .AND. hypi(jk+1) >= crm_top) THEN
        max_lev_crm = jk
        EXIT
     ENDIF
  END DO
  max_lev_crm = MAX(1,max_lev_crm)
  crm_nz  = nlevp1 - max_lev_crm
  IF (p_parallel_io) THEN
     PRINT*,                                                       &
          substr, ': restrict CRM to a level number higher than ', &
          max_lev_crm, ' !'
     PRINT*,                                        &
          substr, ': number of CRM levels (corresponding to lowermost (GCM) host model level): ', &
          crm_nz, ' !'
  END IF

  crm_nx=ngridcells_x
  crm_ny=ngridcells_y

  print*, 'set up of crmvars:  ', crmvars

  CALL set_crm_dims(ngridcells_x, ngridcells_y, crm_nz,    &
        crm_size, crm_timestep, CRM_3D, crmvars)
  !um_hr_20190307-

  ! ====================================================================

  CALL start_message_bi(modstr,'MEMORY INITIALIZATION', substr)

  CALL new_channel(status, modstr, reprid=GP_3D_MID)
  CALL channel_halt(substr, status)

  CALL new_channel_object_reference(status, 'g3b', 'acdnc', modstr, 'acdnc')
  CALL channel_halt(substr, status)
  
  CALL new_channel_object(status, modstr, 'tl' , p3=tl )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'tl' &
         , 'long_name', c='Temperature')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'tl' &
         , 'units', c='K')
  CALL channel_halt(substr, status)
  
  CALL new_channel_object(status, modstr, 'ql' , p3=ql )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'ql' &
         , 'long_name', c='Water Vapour')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'ql' &
         , 'units', c='g/g')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qccl' , p3=qccl )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qccl' &
         , 'long_name', c='Cloud Liquid Water')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qccl' &
         , 'units', c='g/g')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qiil' , p3=qiil )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qiil' &
         , 'long_name', c='Cloud Ice')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qiil' &
         , 'units', c='g/g')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'cltot' , p2=cltot, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cltot' &
         , 'long_name', c='total cloud cover')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cltot' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'clhgh' , p2=clhgh, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'clhgh' &
         , 'long_name', c='high cloud cover')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'clhgh' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'clmed' , p2=clmed, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'clmed' &
         , 'long_name', c='medium cloud cover')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'clmed' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'cllow' , p2=cllow, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cllow' &
         , 'long_name', c='low cloud cover')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cllow' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'ultend' , p3=ultend )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'ultend' &
         , 'long_name', c='zonal wind tendency')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'ultend' &
         , 'units', c='m/s^2')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'vltend' , p3=vltend )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'vltend' &
         , 'long_name', c='meridional wind tendency')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'vltend' &
         , 'units', c='m/s^2')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'ttend' , p3=ttend )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'ttend' &
         , 'long_name', c='CRM temperature tendency')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'ttend' &
         , 'units', c='K/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'sltend' , p3=sltend )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'sltend' &
         , 'long_name', c='tendency of static energy')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'sltend' &
         , 'units', c='1/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qltend' , p3=qltend )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qltend' &
         , 'long_name', c='tendency of water vapour')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qltend' &
         , 'units', c='1/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qcltend' , p3=qcltend )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qcltend' &
         , 'long_name', c='tendency of cloud water')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qcltend' &
         , 'units', c='1/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qiltend' , p3=qiltend )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qiltend' &
         , 'long_name', c='tendency of cloud ice')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qiltend' &
         , 'units', c='1/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'precc' , p2=precc, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precc' &
         , 'long_name', c='convective precipitation rate')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precc' &
         , 'units', c='m/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'precl' , p2=precl, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precl' &
         , 'long_name', c='stratiform precipitation rate')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precl' &
         , 'units', c='m/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'precsc' , p2=precsc, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precsc' &
         , 'long_name', c='convective snow rate')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precsc' &
         , 'units', c='m/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'precsl' , p2=precsl, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precsl' &
         , 'long_name', c='stratiform snow rate')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precsl' &
         , 'units', c='m/s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'cld' , p3=cld )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cld' &
         , 'long_name', c='cloud fraction')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cld' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'cldtop' , p3=cldtop )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cldtop' &
         , 'long_name', c='cloud top pdf')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cldtop' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'cldbot' , p3=cldbot )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cldbot' &
         , 'long_name', c='cloud bottom pdf')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cldbot' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'convtop' , p3=convtop )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'convtop' &
         , 'long_name', c='convective cloud top pdf')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'convtop' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'convbot' , p3=convbot )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'convbot' &
         , 'long_name', c='convective cloud bottom pdf')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'convbot' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'conv_cth' , p2=conv_cth, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'conv_cth' &
         , 'long_name', c='cloud top level')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'conv_cth' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'conv_cbh' , p2=conv_cbh, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'conv_cbh' &
         , 'long_name', c='cloud bottom level')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'conv_cbh' &
         , 'units', c='')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'cv_cover' , p2=cv_cover, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cv_cover' &
         , 'long_name', c='convective cloud cover of CRM cells')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'cv_cover' &
         , 'units', c='s')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'gicewp' , p3=gicewp )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'gicewp' &
         , 'long_name', c='ice water path')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'gicewp' &
         , 'units', c='g/m^2')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'gliqwp' , p3=gliqwp )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'gliqwp' &
         , 'long_name', c='liquid water path')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'gliqwp' &
         , 'units', c='g/m^2')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'mc' , p3=mc )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mc' &
         , 'long_name', c='massflux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mc' &
         , 'units', c='g/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'mcup' , p3=mcup )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mcup' &
         , 'long_name', c='updraft massflux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mcup' &
         , 'units', c='g/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'mcdn' , p3=mcdn )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mcdn' &
         , 'long_name', c='downdraft massflux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mcdn' &
         , 'units', c='g/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'mcuup' , p3=mcuup )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mcuup' &
         , 'long_name', c='unsaturated updraft massflux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mcuup' &
         , 'units', c='g/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'mcudn' , p3=mcudn )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mcudn' &
         , 'long_name', c='unsaturated downdraft massflux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'mcudn' &
         , 'units', c='g/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'massfu' , p3=massfu )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'massfu' &
         , 'long_name', c='convective updraft massflux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'massfu' &
         , 'units', c='kg/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'massfd' , p3=massfd )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'massfd' &
         , 'long_name', c='convective downdraft massflux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'massfd' &
         , 'units', c='kg/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'u_entr' , p3=u_entr )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'u_entr' &
         , 'long_name', c='convective updraft entrainment')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'u_entr' &
         , 'units', c='kg/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'u_detr' , p3=u_detr )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'u_detr' &
         , 'long_name', c='convective updraft detrainment')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'u_detr' &
         , 'units', c='kg/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'd_entr' , p3=d_entr )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'd_entr' &
         , 'long_name', c='convective downdraft entrainment')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'd_entr' &
         , 'units', c='kg/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'd_detr' , p3=d_detr )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'd_detr' &
         , 'long_name', c='convective downdraft detrainment')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'd_detr' &
         , 'units', c='kg/(m^2 s)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'crm_qc' , p3=crm_qc )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qc' &
         , 'long_name', c='cloud water')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qc' &
         , 'units', c='g/(m^2)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'crm_qi' , p3=crm_qi )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qi' &
         , 'long_name', c='cloud ice')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qi' &
         , 'units', c='g/(m^2)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'crm_qs' , p3=crm_qs )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qs' &
         , 'long_name', c='cloud snow')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qs' &
         , 'units', c='g/(m^2)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'crm_qg' , p3=crm_qg )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qg' &
         , 'long_name', c='cloud graupel')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qg' &
         , 'units', c='g/(m^2)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'crm_qr' , p3=crm_qr )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qr' &
         , 'long_name', c='cloud rain')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'crm_qr' &
         , 'units', c='g/(m^2)')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tkez' , p3=tkez )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'tkez' &
         , 'long_name', c='TKE profile')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'tkez' &
         , 'units', c='J/kg')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tkesgsz' , p3=tkesgsz )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'tkesgsz' &
         , 'long_name', c='sub-grid scale TKE profile')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'tkesgsz' &
         , 'units', c='J/kg')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'flux_u' , p3=flux_u )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'flux_u' &
         , 'long_name', c='zonal momentum flux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'flux_u' &
         , 'units', c='m^2/s^2')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'flux_v' , p3=flux_v )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'flux_v' &
         , 'long_name', c='meridional momentum flux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'flux_v' &
         , 'units', c='m^2/s^2')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'flux_qp' , p3=flux_qp )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'flux_qp' &
         , 'long_name', c='precipitating water flux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'flux_qp' &
         , 'units', c='')  !!! UNITS???
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'flux_qt' , p3=flux_qt )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'flux_qt' &
         , 'long_name', c='non-precipitating water flux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'flux_qt' &
         , 'units', c='')  !!! UNITS???
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'fluxsgs_qt' , p3=fluxsgs_qt )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'fluxsgs_qt' &
         , 'long_name', c='sub-grid scale non-precipitating water flux')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'fluxsgs_qt' &
         , 'units', c='')  !!! UNITS???
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'pflx' , p3=pflx )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'pflx' &
         , 'long_name', c='precipitation flux (3D)')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'pflx' &
         , 'units', c='kg/(m^2 s)') 
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qt_ls' , p3=qt_ls )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qt_ls' &
         , 'long_name', c='tendency of non-precipitating water due to large scale')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qt_ls' &
         , 'units', c='')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qt_trans' , p3=qt_trans )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qt_trans' &
         , 'long_name', c='tendency of non-precipitating water due to transport')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qt_trans' &
         , 'units', c='')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qp_trans' , p3=qp_trans )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qp_trans' &
         , 'long_name', c='tendency of precipitating water due to transport')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qp_trans' &
         , 'units', c='')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qp_fall' , p3=qp_fall )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qp_fall' &
         , 'long_name', c='tendency of precipitating water due to fall-out')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qp_fall' &
         , 'units', c='')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qp_src' , p3=qp_src )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qp_src' &
         , 'long_name', c='tendency of precipitating water due to conversion')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qp_src' &
         , 'units', c='')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'qp_evp' , p3=qp_evp )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qp_evp' &
         , 'long_name', c='tendency of precipitating water due to evaporation')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'qp_evp' &
         , 'units', c='')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'prectend' , p2=prectend, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'prectend' &
         , 'long_name', c='column integrated tendency in precipitation water+ice')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'prectend' &
         , 'units', c='kg/(m^2 s)')  
  CALL channel_halt(substr, status)

   CALL new_channel_object(status, modstr, 'precstend' , p2=precstend, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precstend' &
         , 'long_name', c='column integrated tendency in precipitation ice')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'precstend' &
         , 'units', c='kg/(m^2 s)')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'taux_crm' , p2=taux_crm, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'taux_crm' &
         , 'long_name', c='zonal CRM surface stress perturbation')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'taux_crm' &
         , 'units', c='N/m^2')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'tauy_crm' , p2=tauy_crm, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'tauy_crm' &
         , 'long_name', c='meridional CRM surface stress perturbation')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'tauy_crm' &
         , 'units', c='N/m^2')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'z0m' , p2=z0m, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'z0m' &
         , 'long_name', c='roughness length')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'z0m' &
         , 'units', c='m')  
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'time_factor' , p2=time_factor, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'time_factor' &
         , 'long_name', c='CRM timing factor')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'time_factor' &
         , 'units', c='-')  
  CALL channel_halt(substr, status)

!!! um_hr_20140929 for cape subroutine
  CALL new_channel_object(status, modstr, 'CAPE_c', p2=cape_c, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'CAPE_c', &
       'long_name', c='maximum CAPE in column')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'CAPE_c', 'units', c='J/kg')
  CALL channel_halt(substr, status)

  CALL new_channel_object(status, modstr, 'CAPE_s', p2=cape_s, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'CAPE_s', &
       'long_name', c='CAPE from surface parcel')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'CAPE_s', 'units', c='J/kg')
  CALL channel_halt(substr, status)
  
  CALL new_channel_object(status, modstr, 'PA_c', p2=PA_c, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'PA_c', &
       'long_name', c='maximum positive contribution to CAPE in column')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'PA_c', 'units', c='J/kg')
  CALL channel_halt(substr, status)
  
  CALL new_channel_object(status, modstr, 'PA_s', p2=PA_s, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'PA_s', &
       'long_name', c='positive contribution to CAPE from surface parcel')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'PA_s', 'units', c='J/kg')
  CALL channel_halt(substr, status)
  
  CALL new_channel_object(status, modstr, 'CAPE_c_idx', p2=CAPE_MAX_IDX, reprid=GP_2D_HORIZONTAL )
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'CAPE_c_idx', &
       'long_name', c='departure level index of maximum CAPE')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'CAPE_c_idx', 'units', c='-')
  CALL channel_halt(substr, status)
  
  CALL new_channel_object(status, modstr, 'PLCL', p2=PLCL, reprid=GP_2D_HORIZONTAL)
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'PLCL', &
       'long_name', c='Lifting condensation level pressure')
  CALL channel_halt(substr, status)
  CALL new_attribute(status, modstr, 'PLCL', 'units', c='hPa')
  CALL channel_halt(substr, status)
!!! um_hr_20140929 for cape routine-

  ! ALLOCATE memory for CRM permanent arrays
   ALLOCATE(CRM_BUFFER1(nproma,crm_nx,crm_ny,crm_nz,crmvars,ngpblks,1)); CRM_BUFFER1 = 0._dp
   ALLOCATE(t_rad1(nproma,crm_nx,crm_ny,crm_nz,ngpblks,1)); t_rad1 = 0._dp
   ALLOCATE(qv_rad1(nproma,crm_nx,crm_ny,crm_nz,ngpblks,1)); qv_rad1 = 0._dp
   ALLOCATE(qc_rad1(nproma,crm_nx,crm_ny,crm_nz,ngpblks,1)); qc_rad1 = 0._dp
   ALLOCATE(qi_rad1(nproma,crm_nx,crm_ny,crm_nz,ngpblks,1)); qi_rad1 = 0._dp
   ALLOCATE(qrad_crm1(nproma,crm_nx,crm_ny,crm_nz,ngpblks,1)); qrad_crm1 = 0._dp

   crm_buffer => crm_buffer1(:,:,:,:,:,:,1)
   t_rad      => t_rad1(:,:,:,:,:,1)
   qv_rad     => qv_rad1(:,:,:,:,:,1)
   qc_rad     => qc_rad1(:,:,:,:,:,1)
   qi_rad     => qi_rad1(:,:,:,:,:,1)
   qrad_crm   => qrad_crm1(:,:,:,:,:,1)

   ALLOCATE(crm_tracer1(nproma,crm_nx,crm_ny,crm_nz,ntrac,ngpblks,1)); crm_tracer1 = 0._dp
   crm_tracer => crm_tracer1(:,:,:,:,:,:,1)


   ! allocate CRM_PTRARRAY with size of gridcells

   CALL new_dimension(status, DIMID_CRMLEV, 'CRM_NLEV', crm_nz)
   CALL channel_halt(substr, status)

   CALL new_representation(status, REPR_CRM_3D_MID, &
         'REPR_CRM_3D_MID'    &
         , rank = 3, link = 'xxx-', dctype = DC_GP               &
         , dimension_ids = (/ DIMID_LON, DIMID_CRMLEV            &
         ,                    DIMID_LAT /)                       &
         , ldimlen       = (/ nproma, AUTO, ngpblks   /)         &
         , output_order  = (/ 1,3,2 /)                           &
         , axis = 'XZY-'                                         &
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
   
   cnt(:,2) = crm_nz
   memu(:,2) = crm_nz
   
   CALL set_representation_decomp(status, REPR_CRM_3D_MID &
        , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
   CALL channel_halt(substr, status)
   
   DEALLOCATE(start) ; NULLIFY(start)
   DEALLOCATE(cnt)   ; NULLIFY(cnt)
   DEALLOCATE(meml)  ; NULLIFY(meml)
   DEALLOCATE(memu)  ; NULLIFY(memu)

   CALL new_dimension(status, DIMID_CRM1LEV, 'CRM_1LEV', 1)
   CALL channel_halt(substr, status)
   CALL new_representation(status, REPR_CRM_3D_1LEV, &
         'REPR_CRM_3D_1LEV'    &
         , rank = 3, link = 'xxx-', dctype = DC_GP               &
         , dimension_ids = (/ DIMID_LON, DIMID_CRM1LEV           &
         ,                    DIMID_LAT /)                       &
         , ldimlen       = (/ nproma, AUTO, ngpblks   /)         &
         , output_order  = (/ 1,3,2 /)                           &
         , axis = 'XZY-'                                         &
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
   
   cnt(:,2) = 1
   memu(:,2) = 1
   
   CALL set_representation_decomp(status, REPR_CRM_3D_1LEV &
        , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
   CALL channel_halt(substr, status)
   
   DEALLOCATE(start) ; NULLIFY(start)
   DEALLOCATE(cnt)   ; NULLIFY(cnt)
   DEALLOCATE(meml)  ; NULLIFY(meml)
   DEALLOCATE(memu)  ; NULLIFY(memu)

! not used, idea for output the subtime steps
! -> subtime steps need to be stored inside the CRM !
   nsubsteps = INT(time_step_len / crm_timestep) 
   IF ( (time_step_len/crm_timestep) - INT(time_step_len/crm_timestep)>1.e-5_dp ) &
        nsubsteps = nsubsteps + 1 
   IF (crm_timestep > time_step_len) nsubsteps = 1
   
   nsub = INT(nsubsteps / VAL_OUT)
   IF ( (nsubsteps/VAL_OUT) - INT(nsubsteps/VAL_OUT) > 1.e-5_dp) &
        nsub = nsub + 1
   nsub = MIN(nsub, nsubsteps)

   ngridcells=ngridcells_x * ngridcells_y
! end of subtime step interval calculations

   ! um_hr_20151016+
   IF (p_parallel_io) WRITE(*,*) 'add new channel '//crmsg//'....'
   
   CALL new_channel(status, crmsg, reprid=GP_3D_MID,lrestreq=.true.)
   CALL channel_halt(substr, status)

   DO i=1,NOUTSG
      ALLOCATE(xcrmoutsg(i)%sg(crm_nx,crm_ny))
   END DO

   ! get representation indices and ranks
   repr_idx_sg(:) = REPR_UNDEF
   DO i=1, NOUTSG
      CALL get_representation_id(status, TRIM(repr_name_sg(i)), repr_idx_sg(i))
      CALL get_representation_info(status, inpname='' &
           , id=repr_idx_sg(i), rank=irank_sg(i))
      CALL channel_halt(substr, status)
   END DO

   DO nx=1,crm_nx
      CALL int2str(istr,nx)
      DO ny=1,crm_ny
         CALL int2str(jstr,ny)
         crm_out_objs: DO i=1,NOUTSG
            SELECT CASE(irank_sg(i))
            CASE(0)
               CALL new_channel_object(status, TRIM(crmsg)                    &
                    , TRIM(setname_sg(i))//'_X'//TRIM(istr)//'_Y'//TRIM(jstr) &
                    , p0=xcrmoutsg(i)%sg(nx,ny)%ptr0                          &
                    , reprid=repr_idx_sg(i) )
            CASE(2)
               CALL new_channel_object(status, TRIM(crmsg)                    &
                    , TRIM(setname_sg(i))//'_X'//TRIM(istr)//'_Y'//TRIM(jstr) &
                    , p2=xcrmoutsg(i)%sg(nx,ny)%ptr2                          &
                    , reprid=repr_idx_sg(i) )
            CASE(3)
               CALL new_channel_object(status, TRIM(crmsg)                    &
                    , TRIM(setname_sg(i))//'_X'//TRIM(istr)//'_Y'//TRIM(jstr) &
                    , p3=xcrmoutsg(i)%sg(nx,ny)%ptr3                          &
                    , reprid=repr_idx_sg(i) )
            CASE(4)
               CALL new_channel_object(status, TRIM(crmsg)                    &
                    , TRIM(setname_sg(i))//'_X'//TRIM(istr)//'_Y'//TRIM(jstr) &
                    , p4=xcrmoutsg(i)%sg(nx,ny)%ptr4                          &
                    , reprid=repr_idx_sg(i) )
            END SELECT

            CALL channel_halt(substr//'channel object '//setname_sg(i)// &
                 &'_X'//TRIM(istr)//'_Y'//TRIM(jstr)//                   &
                 &' not found', status)
            
            CALL new_attribute(status, TRIM(crmsg), TRIM(setname_sg(i))// &
                 &'_X'//TRIM(istr)//'_Y'//TRIM(jstr), 'long_name'         &
                 , c=TRIM(longname_sg(i))//'for box X'//TRIM(istr)//' Y'//TRIM(jstr))
            CALL channel_halt(substr, status)
            
            CALL new_attribute(status, TRIM(crmsg), TRIM(setname_sg(i))// &
                 &'_X'//TRIM(istr)//'_Y'//TRIM(jstr),'units'              &
                 , c=TRIM(unit_sg(i)) )
            CALL channel_halt(substr, status)
         END DO crm_out_objs
      END DO
   END DO
   ! um_hr_20151016-

         ! CRM buffer
         IF (p_parallel_io) WRITE(*,*) 'add new channel '//crmsgb//'....'

         CALL new_channel(status, crmsgb, reprid=REPR_CRM_3D_MID,lrestreq=.true.)
         CALL channel_halt(substr, status)
         
         print*, 'number of crmvars: ', crmvars

         DO i=1,crm_nx
            CALL int2str(istr,i)
            DO j=1,crm_ny
               CALL int2str(jstr,j)
               DO k=1,crmvars
                  CALL int2str(kstr,k)
                  SELECT CASE (k)
                  CASE(1)
                     crmbufname  = "u"
                     crmbufunit  = "m/s"
                     crmlongname = "zonal wind speed"
                  CASE(2)
                     crmbufname  = "v"
                     crmbufunit  = "m/s"
                     crmlongname = "meridional wind speed"
                  CASE(3)
                     crmbufname  = "w"
                     crmbufunit  = "m/s"
                     crmlongname = "vertical wind speed"
                  CASE(4)
                     crmbufname  = "temp"
                     crmbufunit  = "K"
                     crmlongname = "temperature"
                  CASE(5)
                     crmbufname  = "q"      ! total non-precipitating water
                     crmbufunit  = "kg/kg"
                     crmlongname = "total non-precipitating water"
                  CASE(6)
                     crmbufunit  = "kg/kg"
                     IF (micro_scheme.eq.0) THEN
                        crmbufname  = "qc"     ! total precipitating water
                        crmlongname = "total precipitating water"
                     ELSE
                        crmbufname  = "qcl"    ! total precipitating water
                        crmlongname = "cloud water mass mixing ratio"
                     END IF
                  CASE(7)
                     IF (micro_scheme.eq.0) THEN
                        crmbufname  = "qn"     ! cloud condensate (water + ice)
                        crmbufunit  = "kg/kg"
                        crmlongname = "cloud condensate (water + ice)"
                     ELSE
                        IF (dopredictNc) THEN
                           crmbufname  = "qnc" 
                           crmbufunit  = "#/kg"
                           crmlongname = "cloud droplet number mixing ratio"
                        ELSE
                           crmbufname  = "qr"  
                           crmbufunit  = "kg/kg"
                           crmlongname = "rain mass mixing ratio"
                        END IF
                     END IF
                  CASE(8)
                     IF (dopredictNc) THEN
                        crmbufname  = "qr"  
                        crmbufunit  = "kg/kg"
                        crmlongname = "rain mass mixing ratio"
                     ELSE
                        crmbufname  = "qnr"  
                        crmbufunit  = "#/kg"
                        crmlongname = "rain number mixing ratio"
                     END IF
                  CASE(9)
                     IF (dopredictNc) THEN
                        crmbufname  = "qnr"  
                        crmbufunit  = "#/kg"
                        crmlongname = "rain number mixing ratio"
                     ELSE
                        crmbufname  = "qci"  
                        crmbufunit  = "kg/kg"
                        crmlongname = "cloud ice mass mixing ratio"
                     END IF
                  CASE(10)
                     IF (dopredictNc) THEN
                        crmbufname  = "qci"  
                        crmbufunit  = "kg/kg"
                        crmlongname = "cloud ice mass mixing ratio"
                     ELSE
                        crmbufname  = "qni"  
                        crmbufunit  = "#/kg"
                        crmlongname = "cloud ice number mixing ratio"
                     END IF
                  CASE(11)
                     IF (dopredictNc) THEN
                        crmbufname  = "qni"  
                        crmbufunit  = "#/kg"
                        crmlongname = "cloud ice number mixing ratio"
                     ELSE
                        crmbufname  = "qs"  
                        crmbufunit  = "kg/kg"
                        crmlongname = "snow mass mixing ratio"
                     END IF
                  CASE(12)
                     IF (dopredictNc) THEN
                        crmbufname  = "qs"  
                        crmbufunit  = "kg/kg"
                        crmlongname = "snow mass mixing ratio"
                     ELSE
                        crmbufname  = "qns"  
                        crmbufunit  = "#/kg"
                        crmlongname = "snow number mixing ratio"
                     END IF
                  CASE(13)
                     IF (dopredictNc) THEN
                        crmbufname  = "qns"  
                        crmbufunit  = "#/kg"
                        crmlongname = "snow number mixing ratio"
                     ELSE
                        crmbufname  = "qg"  
                        crmbufunit  = "kg/kg"
                        crmlongname = "graupel mass mixing ratio"
                     END IF
                  CASE(14)
                     IF (dopredictNc) THEN
                        crmbufname  = "qg"  
                        crmbufunit  = "kg/kg"
                        crmlongname = "graupel mass mixing ratio"
                     ELSE
                        crmbufname  = "qng"  
                        crmbufunit  = "#/kg"
                        crmlongname = "graupel number mixing ratio"
                     END IF
                  CASE(15)
                     crmbufname  = "qng"  
                     crmbufunit  = "#/kg"
                     crmlongname = "graupel number mixing ratio"
                  END SELECT

                  MEM => crm_buffer1(:,i,j,:,k,:,:)
                  CALL new_channel_object(status, crmsgb,       &
                       TRIM(crmbufname)//'_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                       mem = mem)
                  CALL channel_halt(substr, status)
                  CALL new_attribute(status, crmsgb, TRIM(crmbufname)//'_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                       'long_name', &
                       c='CRM '//TRIM(crmlongname)//' for box '//TRIM(istr)//' X '//TRIM(jstr)//' Y' )
                  
                  CALL channel_halt(substr, status)
                  CALL new_attribute(status, crmsgb, TRIM(crmbufname)//'_X'//TRIM(istr)//'_Y'//TRIM(jstr) ,&
                       'units', c=TRIM(crmbufunit) )
                  CALL channel_halt(substr, status)            
               END DO
            END DO
         END DO


! quantities in combination with radiation
         IF (p_parallel_io) WRITE(*,*) 'add new channel '//crmsgr//'....'

         CALL new_channel(status, crmsgr, reprid=REPR_CRM_3D_MID,lrestreq=.true.)
         CALL channel_halt(substr, status)
         
         DO i=1,crm_nx
            CALL int2str(istr,i)
            DO j=1,crm_ny
               CALL int2str(jstr,j)
         
               MEM => t_rad1(:,i,j,:,:,:)
               CALL new_channel_object(status, crmsgr,       &
                    'rad_temp_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    mem = mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'rad_temp_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    'long_name', c='CRM Radiation Temperature for box '//TRIM(istr)//' X '//TRIM(jstr)//' Y' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'rad_temp_X'//TRIM(istr)//'_Y'//TRIM(jstr) ,'units', c='K')
               CALL channel_halt(substr, status)
               
               MEM => qv_rad1(:,i,j,:,:,:)
               CALL new_channel_object(status, crmsgr,       &
                    'rad_qv_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    mem = mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'rad_qv_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    'long_name', c='CRM Radiation water vapour for box '//TRIM(istr)//' X '//TRIM(jstr)//' Y' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'rad_qv_X'//TRIM(istr)//'_Y'//TRIM(jstr) ,'units', c='kg/kg')
               CALL channel_halt(substr, status)
               
               MEM => qc_rad1(:,i,j,:,:,:)
               CALL new_channel_object(status, crmsgr,       &
                    'rad_qc_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    mem = mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'rad_qc_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    'long_name', c='CRM Radiation cloud water for box '//TRIM(istr)//' X '//TRIM(jstr)//' Y' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'rad_qc_X'//TRIM(istr)//'_Y'//TRIM(jstr) ,'units', c='kg/kg')
               CALL channel_halt(substr, status)
               
               MEM => qi_rad1(:,i,j,:,:,:)
               CALL new_channel_object(status, crmsgr,       &
                    'rad_qi_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    mem = mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'rad_qi_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    'long_name', c='CRM Radiation cloud ice for box '//TRIM(istr)//' X '//TRIM(jstr)//' Y' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'rad_qi_X'//TRIM(istr)//'_Y'//TRIM(jstr) ,'units', c='kg/kg')
               CALL channel_halt(substr, status)
               
               MEM => qrad_crm1(:,i,j,:,:,:)
               CALL new_channel_object(status, crmsgr,       &
                    'qrad_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    mem = mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'qrad_X'//TRIM(istr)//'_Y'//TRIM(jstr), &
                    'long_name', c='CRM Radiation quantity '//TRIM(istr)//' X '//TRIM(jstr)//' Y' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgr, 'qrad_X'//TRIM(istr)//'_Y'//TRIM(jstr) ,'units', c='kg/kg')
               CALL channel_halt(substr, status)
            END DO
         END DO

   ! NEW DIMENSION (ntr = no. of tracers)
   IF (p_parallel_io) WRITE(*,*) 'add new channel '//crmsgt//'....'

   IF (dosgtracer .gt. 0) THEN
      CALL new_channel(status, crmsgt, reprid=REPR_CRM_3D_MID,lrestreq=.true.)
      CALL channel_halt(substr, status)
      
      DO jt = 1,ntrac
         CALL int2str(tstr,jt)
         DO i = 1,crm_nx
            CALL int2str(istr,i)
            DO j = 1,crm_ny
               CALL int2str(jstr,j)
               
               MEM => crm_tracer1(:,i,j,:,jt,:,:)
               CALL new_channel_object(status, crmsgt, 'tracer_'//TRIM(tstr)//'_X'//TRIM(istr)//'_Y'//TRIM(jstr) &
                    , mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, crmsgt, 'tracer_'//TRIM(tstr)//'_X'//TRIM(istr)//'_Y'//TRIM(jstr) &
                    , 'units', c='-')
               CALL channel_halt(substr, status)
               
            END DO
         END DO
      END DO
   END IF
         
   CALL set_crm_tracers(ntrac)   ! -> set number of CRM tracers
   
   CALL crm_allocate_vars
   CALL crm_allocate_microphysics
   CALL crm_allocate_tracers

   CALL end_message_bi(modstr,'MEMORY INITIALIZATION', substr)

 END SUBROUTINE crm_init_memory

  ! ====================================================================

  SUBROUTINE crm_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object, new_channel, &
                                           get_channel_info  , new_attribute

    USE messy_main_data_bi,          ONLY: basemodstr=>modstr,              &
                                           l_heatflux, s_heatflux

    IMPLICIT NONE

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'crm_init_coupling'
    CHARACTER(LEN=12)           :: vdicha = '' 

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

! um_hr_20190301+
    CALL get_channel_object(status, TRIM(basemodstr), 'wind10', p2=wind10_2d)
    CALL channel_halt(substr, status)

#ifdef ECHAM5
    CALL get_channel_info(status, 'e5vdiff')
    IF (status /= 0) THEN
       CALL warning_bi( &
            'channel e5vdiff not available, trying vertex ...' &
            , substr)
       CALL get_channel_info(status, 'vertex')
       IF (status /= 0) THEN
          CALL error_bi(' ... vertex also not available!', substr)
       ELSE
          vdicha = 'vertex'
       ENDIF
    ELSE
       vdicha = 'e5vdiff'
    ENDIF
#endif

    CALL get_channel_object(status,TRIM(vdicha),'ustrl', p2=ustrl)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,TRIM(vdicha),'ustrw', p2=ustrw)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,TRIM(vdicha),'ustri', p2=ustri)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,TRIM(vdicha),'vstrl', p2=vstrl)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,TRIM(vdicha),'vstrw', p2=vstrw)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status,TRIM(vdicha),'vstri', p2=vstri)
    CALL channel_halt(substr, status)

    IF (.NOT. ASSOCIATED(l_heatflux)) &
         call error_bi('l_heatflux not associated', substr)
    IF (.NOT. ASSOCIATED(s_heatflux)) &
         call error_bi('s_heatflux not associated', substr)
! um_hr_20190301-

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE crm_init_coupling

  ! ====================================================================

 
  SUBROUTINE CRM_RADIATION

    USE messy_main_timer,           ONLY: lstart, delta_time

    USE messy_main_grid_def_mem_bi, ONLY: nlev, nlevp1, jrow
    USE messy_main_grid_def_bi,     ONLY: grmass, grvol
    USE messy_main_data_bi,         ONLY: apm1, aphm1,               & 
                                         loland_2d, loglac_2d, acdnc,     &
                                         slf, slm, seaice,          &
                                         s_heatflux, l_heatflux,    & 
                                         aclc, aprc, aprl, aprs,    & ! cloud cover, precip. rates
                                         rsfl_2d, ssfl_2d,          &
                                         rsfc_2d, ssfc_2d ! surface rain/snow rates

    USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1, ntrac => ntrac_gp

    USE messy_crm,                 ONLY: crm_orient, calc_cape_crm, cloud_droplet_nc_crm


    REAL(dp) :: ocnfrac(kproma)              ! area fraction of the ocean
    REAL(dp) :: ustr(kproma), vstr(kproma)   ! wind stress
    REAL(dp) :: tau00(kproma)                ! large-scale surface stress (N/m2)
    REAL(dp) :: bflxls(kproma)               ! large-scale surface buoyancy flux (K m/s)
    
    REAL(dp) :: conv_counter(kproma)         ! convectively active CRM cells
    REAL(dp) :: frl(kproma)                  ! fraction of land
    REAL(dp) :: frw(kproma)                  ! fraction of water
    REAL(dp) :: fri(kproma)                  ! fraction of seaice
    REAL(dp) :: pdel(kproma,nlev)            ! pressure layer thickness
    REAL(dp) :: zmid(kproma,nlev)
    REAL(dp) :: zint(kproma,nlevp1)
    REAL(dp) :: zrho(kproma,nlev)            ! air density
    REAL(dp) :: cdnc(kproma,nlev)            ! cloud droplet number mixing ratio

    REAL(dp) :: scale
    REAL(dp) :: xygridcells                  ! no. of total CRM gridcells
    
    REAL(dp) :: pxtp1(kproma,nlev,ntrac)     ! tracer variable
    REAL(dp) :: xtte_crm(kproma,nlev,ntrac)  ! crm tracer tendency
    REAL(dp) :: xtp_crm(crm_nx,crm_ny,crm_nz,ntrac) ! local crm tracer variable
    REAL(dp) :: xtpm1_crm(crm_nz,ntrac)      ! local crm tracer variable of the former timestep (averaged)

    INTEGER  :: IDX(nlev)
    INTEGER  :: lchnk               ! chunk identifier (not used); dummy variable

    INTEGER  :: m                   ! loop variable; swapped global grid levels (31 -> 1, 30 -> 2, etc.)
    INTEGER  :: jl, jk, jt, nx, ny

    ! local crm variables
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: qv_crm     => NULL() ! CRM grid water vapour(g/g)
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: qc_crm     => NULL() ! CRM grid cloud liquid water (g/g)
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: qi_crm     => NULL() ! CRM grid cloud ice (g/g)
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: qpc_crm    => NULL() ! CRM grid precipitated water (g/g)
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: qpi_crm    => NULL() ! CRM grid precipitated ice (g/g)
    REAL(dp), POINTER, DIMENSION(:,:,:)   :: prec_crm   => NULL() ! CRM grid precipitation
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: cld_crm    => NULL() ! CRM cloud cover (3D)
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: cdnc_crm   => NULL() ! CRM cloud droplet number concentration
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: radlp_crm  => NULL() ! CRM effective radii for liquid droplets
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: radip_crm  => NULL() ! CRM effective radii for ice droplets

    REAL(dp), POINTER, DIMENSION(:)   :: pustrl    => NULL()
    REAL(dp), POINTER, DIMENSION(:)   :: pustrw    => NULL()
    REAL(dp), POINTER, DIMENSION(:)   :: pustri    => NULL()
    REAL(dp), POINTER, DIMENSION(:)   :: pvstrl    => NULL()
    REAL(dp), POINTER, DIMENSION(:)   :: pvstrw    => NULL()
    REAL(dp), POINTER, DIMENSION(:)   :: pvstri    => NULL()
    REAL(dp), POINTER, DIMENSION(:)   :: wind10    => NULL()

    REAL(dp), DIMENSION(nproma,nlev)  :: ptm1, pqm1, pxlm1, pxim1, pum1, pvm1

#ifdef MESSYTENDENCY
    REAL(dp), DIMENSION(nproma,nlev)  :: lo_tte, lo_qte, lo_xlte, lo_xite, lo_vom, lo_vol
#endif

    EXTERNAL :: CRM

    lchnk = 1
    conv_counter = 0.0_dp
    xygridcells = 1._dp / REAL((crm_nx * crm_ny),dp)

    ALLOCATE(qv_crm(nproma,crm_nx,crm_ny,crm_nz));    qv_crm (:,:,:,:)   = 0.0_dp
    ALLOCATE(qc_crm(nproma,crm_nx,crm_ny,crm_nz));    qc_crm (:,:,:,:)   = 0.0_dp
    ALLOCATE(qi_crm(nproma,crm_nx,crm_ny,crm_nz));    qi_crm (:,:,:,:)   = 0.0_dp
    ALLOCATE(qpc_crm(nproma,crm_nx,crm_ny,crm_nz));   qpc_crm(:,:,:,:)   = 0.0_dp
    ALLOCATE(qpi_crm(nproma,crm_nx,crm_ny,crm_nz));   qpi_crm(:,:,:,:)   = 0.0_dp
    ALLOCATE(prec_crm(nproma,crm_nx,crm_ny));         prec_crm(:,:,:)    = 0.0_dp
    ALLOCATE(cld_crm(nproma,crm_nx,crm_ny,crm_nz));   cld_crm(:,:,:,:)   = 0.0_dp
    ALLOCATE(cdnc_crm(nproma,crm_nx,crm_ny,crm_nz));  cdnc_crm(:,:,:,:)  = 0.0_dp
    ALLOCATE(radlp_crm(nproma,crm_nx,crm_ny,crm_nz)); radlp_crm(:,:,:,:) = 0.0_dp
    ALLOCATE(radip_crm(nproma,crm_nx,crm_ny,crm_nz)); radip_crm(:,:,:,:) = 0.0_dp

    ! calculate rho
    zrho(1:kproma,1:nlev) = grmass(1:kproma,1:nlev,jrow) / &
                            grvol (1:kproma,1:nlev,jrow)

    FIRST_TIMESTEP: IF (LSTART) THEN

       ! standard ECHAM5 cloud droplet number concentration
       ! calculated only at initialisation - no additional calculation during runtime
       ! can be expanded as described in cloud.nml
       ! when using one-mom micro. scheme in CRM cdnc's are only used for cloudopt calculation (no usage inside CRM)
       ! when using two-mom micro. and if dopredictNc=true this prescribes the default cdnc profile for CRM
       ! and changes it accordingly during runtime
       CALL cloud_droplet_nc_crm(kproma, nproma, nlev, apm1, &
                                acdnc(1:kproma,:,jrow), loland_2d(1:kproma,jrow), loglac_2d(1:kproma,jrow))
       
       DO m=1,crm_nz
          jk = nlev-m+1
          DO jl=1,kproma
             qrad_crm(jl,:,:,m,jrow) = 0.0_dp

             t_rad  (jl,:,:,m,jrow) = tm1(jl,jk,jrow) + tte_3d(jl,jk,jrow) * time_step_len
             qv_rad (jl,:,:,m,jrow) = qm1(jl,jk,jrow) + qte_3d(jl,jk,jrow) * time_step_len
             qc_rad (jl,:,:,m,jrow) = xlm1(jl,jk,jrow) + xlte_3d(jl,jk,jrow) * time_step_len
             qi_rad (jl,:,:,m,jrow) = xim1(jl,jk,jrow) + xite_3d(jl,jk,jrow) * time_step_len
             if (crm_orient .eq. 1 .and. CRM_3D .eq. 0) then 
                !!!! CRM ORIENTATION: NORTH - SOUTH
                crm_buffer(jl,:,:,m,1,jrow) = vm1(jl,jk,jrow) + vol_3d(jl,jk,jrow) * time_step_len
                crm_buffer(jl,:,:,m,2,jrow) = um1(jl,jk,jrow) + vom_3d(jl,jk,jrow) * time_step_len
             else 
                !!!! CRM ORIENTATION: EAST - WEST OR FULL 3D CRM
                crm_buffer(jl,:,:,m,1,jrow) = um1(jl,jk,jrow) + vom_3d(jl,jk,jrow) * time_step_len
                crm_buffer(jl,:,:,m,2,jrow) = vm1(jl,jk,jrow) + vol_3d(jl,jk,jrow) * time_step_len
             endif
             crm_buffer(jl,:,:,m,3,jrow) = 0.0_dp
             crm_buffer(jl,:,:,m,4,jrow) = t_rad(jl,:,:,m,jrow)
             IF (micro_scheme.eq.0) THEN ! use one-moment micro scheme: SAM1MOM
                crm_buffer(jl,:,:,m,5,jrow) = qv_rad(jl,:,:,m,jrow) + qc_rad(jl,:,:,m,jrow) + qi_rad(jl,:,:,m,jrow)
                crm_buffer(jl,:,:,m,6,jrow) = 0.0_dp
                crm_buffer(jl,:,:,m,7,jrow) = qc_rad(jl,:,:,m,jrow) + qi_rad(jl,:,:,m,jrow)
             ELSE ! use two-moment micro scheme: MORRISON2005 - not thoroughly tested yet
                crm_buffer(jl,:,:,m,5,jrow) = qv_rad(jl,:,:,m,jrow)
                crm_buffer(jl,:,:,m,6,jrow) = qc_rad(jl,:,:,m,jrow)
                IF (dopredictNc) THEN
                   crm_buffer(jl,:,:,m,7,jrow)  = acdnc(jl,jk,jrow)/zrho(jl,jk) ! set ECHAM5 cdnc profile...-> #/kg
                   crm_buffer(jl,:,:,m,8,jrow)  = 0.0_dp
                   crm_buffer(jl,:,:,m,9,jrow)  = 0.0_dp
                   crm_buffer(jl,:,:,m,10,jrow) = qi_rad(jl,:,:,m,jrow)
                   crm_buffer(jl,:,:,m,11:crmvars,jrow) = 0.0_dp                   
                ELSE
                   crm_buffer(jl,:,:,m,7,jrow)  = 0.0_dp
                   crm_buffer(jl,:,:,m,8,jrow)  = 0.0_dp
                   crm_buffer(jl,:,:,m,9,jrow)  = qi_rad(jl,:,:,m,jrow)
                   crm_buffer(jl,:,:,m,10:crmvars,jrow) = 0.0_dp
                END IF
             END IF
             
             IF (dosgtracer .gt. 0) THEN
                DO jt = 1,ntrac
                   crm_tracer(jl,1:crm_nx,1:crm_ny,m,jt,jrow) = &
                        MAX(pxtm1(jl,jk,jt) + pxtte(jl,jk,jt) * time_step_len, 0._dp)
                END DO
             END IF

          END DO
       END DO

    ELSE

    !!! Q_rad_crm should be initialised with ????
    !!!       qrad_crm(jl,:,:,:) = ????
    !!! when coupling with radiation code include radiative heating for every CRM cell here tbd. !!!
    
    END IF FIRST_TIMESTEP

    cltot(:,jrow)        = 0.0_dp
    clhgh(:,jrow)        = 0.0_dp
    clmed(:,jrow)        = 0.0_dp
    cllow(:,jrow)        = 0.0_dp
    
    cld(:,:,jrow)        = 0.0_dp
    cldtop(:,:,jrow)     = 0.0_dp
    cldbot(:,:,jrow)     = 0.0_dp
    convtop(:,:,jrow)    = 0.0_dp
    convbot(:,:,jrow)    = 0.0_dp
    conv_cth(:,jrow)     = 0.0_dp
    conv_cbh(:,jrow)     = 0.0_dp
    cv_cover(:,jrow)     = 0.0_dp
    cdnc(:,jrow)         = 0.0_dp    

    precc(:,jrow)        = 0.0_dp
    precl(:,jrow)        = 0.0_dp
    precsc(:,jrow)       = 0.0_dp
    precsl(:,jrow)       = 0.0_dp
    gicewp(:,:,jrow)     = 0.0_dp
    gliqwp(:,:,jrow)     = 0.0_dp
    mc(:,:,jrow)         = 0.0_dp
    mcup(:,:,jrow)       = 0.0_dp
    mcdn(:,:,jrow)       = 0.0_dp
    mcuup(:,:,jrow)      = 0.0_dp
    mcudn(:,:,jrow)      = 0.0_dp
    
    massfu(:,:,jrow)     = 0.0_dp
    massfd(:,:,jrow)     = 0.0_dp
    u_entr(:,:,jrow)     = 0.0_dp
    u_detr(:,:,jrow)     = 0.0_dp
    d_entr(:,:,jrow)     = 0.0_dp
    d_detr(:,:,jrow)     = 0.0_dp
    
    crm_qc(:,:,jrow)     = 0.0_dp
    crm_qi(:,:,jrow)     = 0.0_dp
    crm_qs(:,:,jrow)     = 0.0_dp
    crm_qg(:,:,jrow)     = 0.0_dp
    crm_qr(:,:,jrow)     = 0.0_dp
    
    tkez(:,:,jrow)       = 0.0_dp
    tkesgsz(:,:,jrow)    = 0.0_dp
    
    flux_qt(:,:,jrow)    = 0.0_dp
    flux_u(:,:,jrow)     = 0.0_dp
    flux_v(:,:,jrow)     = 0.0_dp
    fluxsgs_qt(:,:,jrow) = 0.0_dp
    flux_qp(:,:,jrow)    = 0.0_dp
    pflx(:,:,jrow)       = 0.0_dp
    
    qt_ls(:,:,jrow)      = 0.0_dp
    qt_trans(:,:,jrow)   = 0.0_dp
    qp_trans(:,:,jrow)   = 0.0_dp
    qp_fall(:,:,jrow)    = 0.0_dp
    qp_evp(:,:,jrow)     = 0.0_dp
    qp_src(:,:,jrow)     = 0.0_dp
    
    prectend(:,jrow)     = 0.0_dp
    precstend(:,jrow)    = 0.0_dp
    
    z0m(:,jrow)          = 0.0_dp
    taux_crm(:,jrow)     = 0.0_dp
    tauy_crm(:,jrow)     = 0.0_dp

    qltend(:,:,jrow)     = 0.0_dp
    qcltend(:,:,jrow)    = 0.0_dp
    qiltend(:,:,jrow)    = 0.0_dp
    sltend(:,:,jrow)     = 0.0_dp
    ttend(:,:,jrow)      = 0.0_dp

    time_factor(:,jrow)  = 0.0_dp
    
    pustrl    => ustrl(:, jrow)
    pustrw    => ustrw(:, jrow)
    pustri    => ustri(:, jrow)
    pvstrl    => vstrl(:, jrow)
    pvstrw    => vstrw(:, jrow)
    pvstri    => vstri(:, jrow)
    wind10    => wind10_2d(:,jrow)

#ifndef MESSYTENDENCY
    ptm1(1:kproma,:)   = tm1(1:kproma,:,jrow)  + tte_3d(1:kproma,:,jrow)   * time_step_len
    pqm1(1:kproma,:)   = MAX(EPSILON(1._dp),qm1(1:kproma,:,jrow)  + qte_3d(1:kproma,:,jrow) * time_step_len)
    pxlm1(1:kproma,:)  = xlm1(1:kproma,:,jrow) + xlte_3d(1:kproma,:,jrow)  * time_step_len
    pxim1(1:kproma,:)  = xim1(1:kproma,:,jrow) + xite_3d(1:kproma,:,jrow)  * time_step_len
    IF (crm_orient .eq. 1 .and. CRM_3D .eq. 0) THEN  !!!! CRM ORIENTATION: NORTH - SOUTH
       pvm1(1:kproma,:)     = um1(1:kproma,:,jrow)  + vom_3d(1:kproma,:,jrow)   * time_step_len
       pum1(1:kproma,:)     = vm1(1:kproma,:,jrow)  + vol_3d(1:kproma,:,jrow)   * time_step_len
    ELSE                                             !!!! CRM ORIENTATION: EAST - WEST OR FULL 3D CRM
       pum1(1:kproma,:)     = um1(1:kproma,:,jrow)  + vom_3d(1:kproma,:,jrow)   * time_step_len
       pvm1(1:kproma,:)     = vm1(1:kproma,:,jrow)  + vol_3d(1:kproma,:,jrow)   * time_step_len
    END IF
    
    DO jt=1,ntrac
       pxtp1(1:kproma,:,jt) = MAX(pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt) * time_step_len, 0._dp)
    END DO
#else
    lo_tte  = 0._dp
    lo_qte  = 0._dp
    lo_xlte = 0._dp
    lo_xite = 0._dp
    lo_vom  = 0._dp
    lo_vol  = 0._dp

    CALL mtend_get_start_l (mtend_id_t, v0 = ptm1)
    CALL mtend_get_start_l (mtend_id_q, v0 = pqm1)
    CALL mtend_get_start_l (mtend_id_xl, v0 = pxlm1) 
    CALL mtend_get_start_l (mtend_id_xi, v0 = pxim1)
    IF (crm_orient .eq. 1 .and. CRM_3D .eq. 0) THEN  !!!! CRM ORIENTATION: NORTH - SOUTH
       CALL mtend_get_start_l (mtend_id_u, v0 = pvm1)  
       CALL mtend_get_start_l (mtend_id_v, v0 = pum1) 
    ELSE !!!! CRM ORIENTATION: EAST - WEST
       CALL mtend_get_start_l (mtend_id_u, v0 = pum1)  
       CALL mtend_get_start_l (mtend_id_v, v0 = pvm1) 
    END IF
    
    DO jt=1, ntrac
! op_pj_20190327+
!!$    CALL mtend_get_start_l (mtend_id_tracer, idt = jt, v0 = pxtp1(1:kproma,1:nlev,jt))
       CALL mtend_get_start_l (jt, v0 = pxtp1(1:kproma,1:nlev,jt))
! op_pj_20190327-
       pxtp1(:,:,jt) = MAX(pxtp1(:,:,jt),0._dp)
    END DO
#endif

!!!!   Calculation of different CRM input variables !!!!
    zmid(1:kproma,:) = (geopot(1:kproma,:,jrow))/g    
    zint(1:kproma,:) = (geopoti(1:kproma,:,jrow))/g   
    DO jk=1,nlev     
       pdel(1:kproma,jk) = aphm1(1:kproma,jk+1) - aphm1(1:kproma,jk)  ! pressure layer thickness
    END DO

    zint(1:kproma,nlev+1) = geopoti(1:kproma,nlev+1,jrow)/g ! global grid height (m)
    ocnfrac(1:kproma)     = 1._dp - slf(1:kproma,jrow)
   
    ! calculate bouyancy flux with kinematic sensible heat flux
    bflxls(1:kproma)  = (s_heatflux(1:kproma,jrow)/cp_air + &
                         0.61*ptm1(1:kproma,nlev)*l_heatflux(1:kproma,jrow)/alv)/zrho(1:kproma,nlev)

    frl(1:kproma) = slm(1:kproma,jrow)
    frw(1:kproma) = (1._dp-slm(1:kproma,jrow))*(1._dp-seaice(1:kproma,jrow))
    fri(1:kproma) = 1._dp-frl(1:kproma)-frw(1:kproma)


    ustr(1:kproma) = frl(1:kproma)*pustrl(1:kproma) + &
                     frw(1:kproma)*pustrw(1:kproma) + &
                     fri(1:kproma)*pustri(1:kproma)
    vstr(1:kproma) = frl(1:kproma)*pvstrl(1:kproma) + &
                     frw(1:kproma)*pvstrw(1:kproma) + &
                     fri(1:kproma)*pvstri(1:kproma)

    ! calculation of tau00 and bflxls:
    tau00(1:kproma) = sqrt(ustr(1:kproma)**2 + vstr(1:kproma)**2)

    global_loop: DO jl=1,kproma  

       DO m=1,crm_nz
          jk = nlev-m+1
          qrad_crm(jl,:,:,m,jrow) = qrad_crm(jl,:,:,m,jrow) / pdel(jl,jk) ! for energy conservation
       END DO

       !!! um_hr_20141210+++
       !!! tracer INITIALISATION for crm
       xtpm1_crm(:,:)   = 0.0_dp
       xtp_crm(:,:,:,:) = 0.0_dp
       dotrac: IF (dosgtracer .gt. 0) THEN
          DO m = 1,crm_nz
             jk = nlev-m+1
             DO jt = 1,ntrac
                IF (dosgtracer == 1) THEN
                   ! do tracer transport with equal CRM grid cell tracer concentrations
                   xtp_crm(:,:,m,jt) = pxtp1(jl,jk,jt)
                END IF
                IF (dosgtracer == 2) THEN
                   ! tracer transport with individual CRM grid cell tracer concentrations
                   DO nx = 1,crm_nx
                      DO ny = 1,crm_ny
                         xtpm1_crm(m,jt)=xtpm1_crm(m,jt) + crm_tracer(jl,nx,ny,m,jt,jrow) 
                      END DO
                   END DO

                   xtpm1_crm = xtpm1_crm * xygridcells
                   
                   DO nx = 1,crm_nx
                      DO ny = 1,crm_ny
                         xtp_crm(nx,ny,m,jt) = MAX(crm_tracer(jl,nx,ny,m,jt,jrow) + &
                              (pxtp1(jl,jk,jt) - xtpm1_crm(m,jt)), 0._dp)
                      END DO
                   END DO
                END IF
                !!!!!!!!!!!!!!!!!!!!!!!!!!   um_hr_20150611+   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!! idea of percentual/relative tracer tendencies per CRM level, not implemented
                !!!!!!!! relative change = (pxtp1 - xtpm1_crm) / xtpm1_crm
                !!!!!!!! xtp_crm(nx,ny,m,jt) = crm_tracer * (1 + relative change)
                !!!!!!!! PROBLEM: IF global tracer conc. = 0 ====> crm tracer conc. are adjusted to 0
                !!!!!!!! high crm concentrations could be damped to zero....
                !!!!!!!! PROBLEM: IF mean crm tracer conc. = 0 ====> no calc. for relative change...
                !!!!!!!! => workaround: individual crm tracer conc. equal global tracer conc.
                !!!!!!!! PROBLEM: IF individual crm tracer conc. = 0 no increase/decrease 
                !!!!!!!! via L-S transport possible
                !!!!!!!! ----- CRM tracer profiles are adapted to global tracer profile ------
                !!!!!!!!!!!!!!!!!!!!!!!!!!   um_hr_20150611-   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             END DO
          END DO
       END IF dotrac
       !!! um_hr_20141210---

       !!! CALCULATION OF CAPE via CAPE-SUBROUTINE OF CONVECT
       CALL calc_cape_crm(nlev, ptm1(jl,:) , apm1(jl,:)  , pqm1(jl,:),    &
                      CAPE_s(jl,jrow)      , CAPE_c(jl,jrow),             &
                      PA_s(jl,jrow)        , PA_c(jl,jrow)  ,             &
                      IDX(1:1)             , plcl(jl,jrow) )
       
       CAPE_MAX_IDX(jl,jrow) = REAL(idx(1),dp)
       
       ! this should contain the main driver and everything concerning output
       CALL CRM(lchnk     , jl           ,                                                           &  
             ptm1(jl,:)   , pqm1(jl,:)   , pxlm1(jl,:) , pxim1(jl,:), pum1(jl,:) , pvm1(jl,:),       &
             aps(jl,jrow) , apm1(jl,:)   , pdel(jl,:)               , geopoti(jl,nlev+1,jrow),       &
             zmid(jl,:)   , zint(jl,:)   , time_step_len            , nlev               ,           &
             ultend(jl,:,jrow)           , vltend(jl,:,jrow)        , qltend(jl,:,jrow)  ,           &
             qcltend(jl,:,jrow)          , qiltend(jl,:,jrow)       , sltend(jl,:,jrow)  ,           &
             crm_buffer(jl,:,:,:,:,jrow) , qrad_crm(jl,:,:,:,jrow)  , qv_crm(jl,:,:,:)   ,           &
             qc_crm(jl,:,:,:)            , qi_crm(jl,:,:,:)         , qpc_crm(jl,:,:,:)  ,           &
             qpi_crm(jl,:,:,:)           , prec_crm(jl,:,:)         , cld_crm(jl,:,:,:)  ,           &
             cdnc_crm(jl,:,:,:)          , radlp_crm(jl,:,:,:)      , radip_crm(jl,:,:,:),           &
             t_rad(jl,:,:,:,jrow)        , qv_rad(jl,:,:,:,jrow)    ,                                &
             qc_rad(jl,:,:,:,jrow)       , qi_rad(jl,:,:,:,jrow)    ,                                &
             precc(jl,jrow)     , precl(jl,jrow)        , precsc(jl,jrow)    , precsl(jl,jrow),      &
             cltot(jl,jrow)     , clhgh(jl,jrow)        , clmed(jl,jrow)     , cllow(jl,jrow) ,      & 
             cld(jl,:,jrow)     , cdnc(jl,:)            , cldtop(jl,:,jrow)  , cldbot(jl,:,jrow) ,   &
             convtop(jl,:,jrow) , convbot(jl,:,jrow)    ,                                            &
             gicewp(jl,:,jrow)  , gliqwp(jl,:,jrow)     ,                                            &
             mc(jl,:,jrow)      , mcup(jl,:,jrow)       , mcdn(jl,:,jrow)    ,                       &
             mcuup(jl,:,jrow)   , mcudn(jl,:,jrow)      ,                                            &
             crm_qc(jl,:,jrow)  , crm_qi(jl,:,jrow)     , crm_qs(jl,:,jrow)  ,                       &
             crm_qg(jl,:, jrow) , crm_qr(jl,:,jrow)     ,                                            & 
             tkez(jl,:,jrow)    , tkesgsz(jl,:,jrow)    ,                                            &
             flux_u(jl,:,jrow)  , flux_v(jl,:,jrow)     ,                                            &
             flux_qt(jl,:,jrow) , fluxsgs_qt(jl,:,jrow) , flux_qp(jl,:,jrow) , pflx(jl,:,jrow),      &
             qt_ls(jl,:,jrow)   , qt_trans(jl,:,jrow)   , qp_trans(jl,:,jrow),                       &
             qp_fall(jl,:,jrow) , qp_evp(jl,:,jrow)     , qp_src(jl,:,jrow)  ,                       &
             ttend(jl,:,jrow)   , prectend(jl,jrow)     , precstend(jl,jrow) ,                       &
             ocnfrac(jl)        , wind10(jl)            , tau00(jl)          , bflxls(jl),           &
             taux_crm(jl,jrow)  , tauy_crm(jl,jrow)     , z0m(jl,jrow)       , time_factor(jl,jrow), &
             conv_counter(jl)   , massfu(jl,:,jrow)     , massfd(jl,:,jrow)  ,                       &
             u_entr(jl,:,jrow)  , u_detr(jl,:,jrow)     , d_entr(jl,:,jrow)  , d_detr(jl,:,jrow),    &
             xtp_crm(:,:,:,:) )

       ! SOME COMMENTS: um_hr_20150220+
       ! geopoti(jl,nlev+1,jrow) = 0.0_dp not used in CRM
       ! um_hr_20150220+

       ! um_hr_20151016+
       ! FOR SUBGRID OUTPUT:
       ! list of CRM variables with its concerning global upper level values
       DO nx=1,crm_nx
          DO ny=1,crm_ny
             DO m=1,crm_nz
                jk = nlev-m+1
                xcrmoutsg(ido_tm1_sg)%sg(nx,ny)%ptr3(jl,jk,jrow)   = crm_buffer(jl,nx,ny,m,4,jrow)
                xcrmoutsg(ido_qm1_sg)%sg(nx,ny)%ptr3(jl,jk,jrow)   = qv_crm(jl,nx,ny,m)
                xcrmoutsg(ido_xlm1_sg)%sg(nx,ny)%ptr3(jl,jk,jrow)  = qc_crm(jl,nx,ny,m)
                xcrmoutsg(ido_xim1_sg)%sg(nx,ny)%ptr3(jl,jk,jrow)  = qi_crm(jl,nx,ny,m)
                xcrmoutsg(ido_aclc_sg)%sg(nx,ny)%ptr3(jl,jk,jrow)  = cld_crm(jl,nx,ny,m)
                xcrmoutsg(ido_acdnc_sg)%sg(nx,ny)%ptr3(jl,jk,jrow) = cdnc_crm(jl,nx,ny,m)
                xcrmoutsg(ido_radlp_sg)%sg(nx,ny)%ptr3(jl,jk,jrow) = radlp_crm(jl,nx,ny,m)
                xcrmoutsg(ido_radip_sg)%sg(nx,ny)%ptr3(jl,jk,jrow) = radip_crm(jl,nx,ny,m)
                xcrmoutsg(ido_prec_sg)%sg(nx,ny)%ptr2(jl,jrow)     = prec_crm(jl,nx,ny)
                xcrmoutsg(ido_qpc_sg)%sg(nx,ny)%ptr3(jl,jk,jrow)   = qpc_crm(jl,nx,ny,m)
                xcrmoutsg(ido_qpi_sg)%sg(nx,ny)%ptr3(jl,jk,jrow)   = qpi_crm(jl,nx,ny,m)
             END DO
             ! set values for upper atmospheric leves (1:nlev-crm_nz)
             xcrmoutsg(ido_tm1_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow)   = tm1(jl,1:nlev-crm_nz,jrow)
             xcrmoutsg(ido_qm1_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow)   = qm1(jl,1:nlev-crm_nz,jrow)
             xcrmoutsg(ido_xlm1_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow)  = xlm1(jl,1:nlev-crm_nz,jrow)
             xcrmoutsg(ido_xim1_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow)  = xim1(jl,1:nlev-crm_nz,jrow)
             xcrmoutsg(ido_aclc_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow)  = 0
             xcrmoutsg(ido_acdnc_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow) = 0
             xcrmoutsg(ido_radlp_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow) = 0
             xcrmoutsg(ido_radip_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow) = 0
             xcrmoutsg(ido_prec_sg)%sg(nx,ny)%ptr2(jl,jrow)                = 0
             xcrmoutsg(ido_qpc_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow)   = 0
             xcrmoutsg(ido_qpi_sg)%sg(nx,ny)%ptr3(jl,1:nlev-crm_nz,jrow)   = 0             
          END DO
       END DO
       ! um_hr_20151016+

!!!!   CALCULATE CRM TRACER TENDENCIES
       dotrac2: IF (dosgtracer .gt. 0) THEN
          xtpm1_crm = 0._dp
          DO jt = 1,ntrac
             DO m = 1,crm_nz
                jk = nlev-m+1
                
                DO nx = 1,crm_nx
                   DO ny = 1,crm_ny
                      xtpm1_crm(m,jt) = xtpm1_crm(m,jt) + xtp_crm(nx,ny,m,jt)
                      crm_tracer(jl,nx,ny,m,jt,jrow) = xtp_crm(nx,ny,m,jt)
                   END DO
                END DO
                
                xtpm1_crm = xtpm1_crm / (crm_nx * crm_ny)
                xtte_crm(jl,jk,jt) = (xtpm1_crm(m,jt) - pxtp1(jl,jk,jt)) / time_step_len
                
             END DO
          END DO
       ELSE
          xtte_crm(jl,:,:) = 0.0_dp
          crm_tracer(jl,:,:,:,:,jrow) = 0.0_dp
       END IF dotrac2

       ! no CRM tendencies above its top
       ! idea of negelecting two uppermost layers of CRM as well... --->  sltend(jl,1:nlev-crm_nz+2,jrow) = 0.0_dp
       sltend(jl,1:nlev-crm_nz,jrow)      = 0.0_dp
       qltend(jl,1:nlev-crm_nz,jrow)      = 0.0_dp
       qcltend(jl,1:nlev-crm_nz,jrow)     = 0.0_dp
       qiltend(jl,1:nlev-crm_nz,jrow)     = 0.0_dp
       ultend(jl,1:nlev-crm_nz,jrow)      = 0.0_dp
       vltend(jl,1:nlev-crm_nz,jrow)      = 0.0_dp
       xtte_crm(jl,1:nlev-crm_nz,1:ntrac) = 0.0_dp

!!!!   UPDATE TENDENCIES
#ifndef MESSYTENDENCY
       DO jk = 1,nlev
          tte_3d(jl,jk,jrow)  = tte_3d(jl,jk,jrow) + (sltend(jl,jk,jrow)/cp_air)
          qte_3d(jl,jk,jrow)  = qte_3d(jl,jk,jrow)  + qltend(jl,jk,jrow)
          xlte_3d(jl,jk,jrow) = xlte_3d(jl,jk,jrow) + qcltend(jl,jk,jrow)
          xite_3d(jl,jk,jrow) = xite_3d(jl,jk,jrow) + qiltend(jl,jk,jrow)

          IF (CRM_3D .eq. 1) THEN   !!! -> wind tendencies only calculated for full 3D CRM 
             vom_3d(jl,jk,jrow)  = vom_3d(jl,jk,jrow)  + ultend(jl,jk,jrow)
             vol_3d(jl,jk,jrow)  = vol_3d(jl,jk,jrow)  + vltend(jl,jk,jrow)
          END IF

          DO jt = 1,ntrac
             pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + xtte_crm(jl,jk,jt)
          END DO
       END DO
#else
       lo_tte(jl,:)  = sltend(jl,:,jrow)/cp_air
       lo_qte(jl,:)  = qltend(jl,:,jrow)
       lo_xlte(jl,:) = qcltend(jl,:,jrow)
       lo_xite(jl,:) = qiltend(jl,:,jrow)
       lo_vom(jl,:)  = ultend(jl,:,jrow)
       lo_vol(jl,:)  = vltend(jl,:,jrow)
#endif

!!!!   UPDATE MICROPHYSICS  ---> convert CRM precipitation rates (m/s) into fluxes (kg m-2 or kg m-2 s-1)
       aprc(jl,jrow) = aprc(jl,jrow) + (precc(jl,jrow)+precsc(jl,jrow))  * rhoh2o * delta_time ! convective precip
       aprl(jl,jrow) = aprl(jl,jrow) + (precl(jl,jrow)+precsl(jl,jrow))  * rhoh2o * delta_time ! large-scale precip
       aprs(jl,jrow) = aprs(jl,jrow) + (precsc(jl,jrow)+precsl(jl,jrow)) * rhoh2o * delta_time ! snow precip

       rsfl_2d(jl,jrow) = precl(jl,jrow)  * rhoh2o  ! large-scale surface rain rate
       ssfl_2d(jl,jrow) = precsl(jl,jrow) * rhoh2o  ! large-scale surface snow rate
       rsfc_2d(jl,jrow) = precc(jl,jrow)  * rhoh2o  ! convective  surface rain rate
       ssfc_2d(jl,jrow) = precsc(jl,jrow) * rhoh2o  ! convective  surface snow rate

!!!!   UPDATE CLOUD COVER and calculate precipitation, cloud top and bottom levels
       DO jk=1,nlev
          aclc(jl,jk,jrow) = cld(jl,jk,jrow)
          pflx(jl,jk,jrow) = pflx(jl,jk,jrow) * rhoh2o  ! convert from m/s to mm/s (or kg m^-2 s^-1)
          !!! idea: use cloudtop/bot pdf to calculate average top and bottom levels !!!
          conv_cth(jl,jrow) = conv_cth(jl,jrow) + convtop(jl,jk,jrow) * jk
          conv_cbh(jl,jrow) = conv_cbh(jl,jrow) + convbot(jl,jk,jrow) * jk
       END DO

!!!!   CALCULATE CONVECTIVE CLOUD COVER OF CRM CELLS
       cv_cover(jl,jrow)  = conv_counter(jl) 
!       cloud_time(jl,jrow) = delta_time * cltot(jl,jrow)         ! time of clouds in CRM
!       cv_time(jl,jrow)    = delta_time * conv_counter(jl,jrow)  ! time of convective clouds in CRM


!!!!   RESCALING OF MASSFLUX -> used for CVTRANS
       scale = 1._dp
       DO jk=1,nlev
          IF (massfu(jl,jk,jrow) .ne. 0._dp) THEN
             IF (ABS(massfd(jl,jk,jrow)) .gt. ABS(0.9*massfu(jl,jk,jrow))) THEN
                scale = MAX(scale, ABS(massfd(jl,jk,jrow)/(0.9*massfu(jl,jk,jrow))))
             END IF
          END IF
       END DO
       
       massfd(jl,:,jrow) = massfd(jl,:,jrow)/scale
       d_entr(jl,:,jrow) = d_entr(jl,:,jrow)/scale
       d_detr(jl,:,jrow) = d_detr(jl,:,jrow)/scale
!       massfu(jl,:,jrow) = massfu(jl,:,jrow) * cv_time(jl,jrow)/delta_time    ! rescaling of upward massflux

    END DO global_loop  ! kproma

    tl(1:kproma,:,jrow)    = ptm1(1:kproma,:)
    ql(1:kproma,:,jrow)    = pqm1(1:kproma,:)
    qccl(1:kproma,:,jrow)  = pxlm1(1:kproma,:)
    qiil(1:kproma,:,jrow)  = pxim1(1:kproma,:)
    IF (micro_scheme.eq.1) THEN ! change cdnc only when using two-mom. micro
       IF (dopredictNc) THEN
          acdnc(1:kproma,:,jrow) = cdnc(1:kproma,:) * zrho(1:kproma,:) ! #/kg -> #/m^3
       END IF
    END IF

#ifdef MESSYTENDENCY
    CALL mtend_add_l (my_handle, mtend_id_t, px = lo_tte)
    CALL mtend_add_l (my_handle, mtend_id_q, px = lo_qte)
    CALL mtend_add_l (my_handle, mtend_id_xl, px = lo_xlte)
    CALL mtend_add_l (my_handle, mtend_id_xi, px = lo_xite)
    IF (CRM_3D .eq. 1) THEN   !!! -> wind tendencies only calculated for full 3D CRM 
       CALL mtend_add_l (my_handle, mtend_id_u, px = lo_vom)
       CALL mtend_add_l (my_handle, mtend_id_v, px = lo_vol)
    END IF

    DO jt = 1,ntrac
! op_pj_20190327+
!!$    CALL mtend_add_l (my_handle, mtend_id_tracer, px = xtte_crm(1:kproma,:,jt), idt=jt)
       CALL mtend_add_l (my_handle, jt, px = xtte_crm(1:kproma,:,jt))
! op_pj_20190327-
    END DO
#endif
    
    IF (ASSOCIATED(qv_crm))    DEALLOCATE(qv_crm)
    IF (ASSOCIATED(qc_crm))    DEALLOCATE(qc_crm)
    IF (ASSOCIATED(qi_crm))    DEALLOCATE(qi_crm)
    IF (ASSOCIATED(qpc_crm))   DEALLOCATE(qpc_crm)
    IF (ASSOCIATED(qpi_crm))   DEALLOCATE(qpi_crm)
    IF (ASSOCIATED(prec_crm))  DEALLOCATE(prec_crm)
    IF (ASSOCIATED(cld_crm))   DEALLOCATE(cld_crm)
    IF (ASSOCIATED(cdnc_crm))  DEALLOCATE(cdnc_crm)
    IF (ASSOCIATED(radlp_crm)) DEALLOCATE(radlp_crm)
    IF (ASSOCIATED(radip_crm)) DEALLOCATE(radip_crm)

  END SUBROUTINE CRM_RADIATION

  ! ====================================================================

  SUBROUTINE CRM_FREE_MEMORY

! ------------------------------------------------------------------
    ! This subroutine is used to deallocate the memory, which has
    ! been "manually" allocated in init_memory.
    ! Note: channel object memory must not be deallocated! This is
    !       performed centrally.
    ! ------------------------------------------------------------------

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'crm_free_memory'
    INTEGER :: j1

    EXTERNAL :: crm_deallocate_vars, crm_deallocate_microphysics
    EXTERNAL :: crm_deallocate_tracers

    CALL crm_deallocate_vars
    CALL crm_deallocate_microphysics
    CALL crm_deallocate_tracers

    DEALLOCATE(CRM_BUFFER1)

    DEALLOCATE(t_rad1)
    DEALLOCATE(qv_rad1)
    DEALLOCATE(qc_rad1)
    DEALLOCATE(qi_rad1)
    DEALLOCATE(qrad_crm1)
    DEALLOCATE(crm_tracer1) 

    ! um_hr_20151016+
    DO j1=1,NOUTSG
       IF (ASSOCIATED(xcrmoutsg(j1)%sg)) DEALLOCATE(xcrmoutsg(j1)%sg)
    END DO
    ! um_hr_20151016+

  END SUBROUTINE CRM_FREE_MEMORY

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE CRM_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    INTEGER                     :: i_do_it
    NAMELIST /CPL/ i_do_it

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='bufly_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
!   INTEGER                     :: i_do_it
    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE CRM_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_crm_e5
! **********************************************************************
