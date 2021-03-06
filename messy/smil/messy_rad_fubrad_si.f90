! ***************************************************************************
#include "messy_main_ppd_bi.inc"

MODULE  messy_rad_fubrad_si
! ***************************************************************************

#if defined(ECHAM5) || defined(CESM1) || defined(MBM_RAD)

  ! BMIL
  USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , info_bi, warning_bi, error_bi
  USE messy_main_grid_def_mem_bi, ONLY: nproma,nlev
  USE messy_main_data_bi,       ONLY: aphm1,apm1
  USE messy_main_timer,         ONLY: delta_time
  USE messy_main_constants_mem, ONLY: dp    ! kind parameter for real
  USE messy_main_channel,       ONLY: t_chaobj_cpl

  ! RAD
  USE messy_rad,                ONLY: t_rad_work, NRADCALL, l_switch

  ! FUBRAD
  USE messy_rad_fubrad_mem,     ONLY: submodstr                        &
                                    , lfubrad, nswlev                  &
                                    , fldo, flup, flupc, flhart, fllya &
                                    , nbands                           &
                                    , cdisse_fubrad                    & 
                                    , sr                               &
                                    , flhuup, flhudo, flchup, flchdo   &
                                    , flhz, flsrb, flsrc               &
                                    , resolution  &
                                    , blev
  USE messy_rad_fubrad_init,    ONLY: rad_fubrad_read_nml_ctrl         &
                                    , fubrad_clean_memory              &  
                                    , fubrad_initialize                &
                                    , fubrad_initialize_fluxes         &
                                    , fubrad_ini_param                 &
                                    , fubrad_initialize_cross_sec      &
                                    , fubrad_initialize_cross_sec_dyn  &
                                    , fubrad_initialize_fluxes_dyn     &
                                    , fubrad_ini_param_dyn

  USE messy_rad_fubrad_srb_km,  ONLY: fubrad_srb_km_init_xs
  USE messy_rad_fubrad_srb_kck, ONLY: fubrad_srb_kck_schu_init
  USE messy_rad_fubrad,         ONLY: prepare_radiation                &
                                    , middle_atmosphere_downward_flux  &
                                    , middle_atmosphere_upward_flux    &
                                    , middle_atmosphere_heat_rates

  IMPLICIT  NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! Variables for CPL namelist
  ! - name of solar cycle data time series channel/object
  TYPE(t_chaobj_cpl) :: fubrad_solar

  namelist /CPL_FUBRAD/ fubrad_solar

  ! IDENTIFIER FOR FUBRAD OUTPUT OBJECTS
  INTEGER, PARAMETER :: ido_heato3flux = 1
  INTEGER, PARAMETER :: ido_heatlya = 2
  INTEGER, PARAMETER :: ido_heatsrc = 3
  INTEGER, PARAMETER :: ido_heatsrb = 4
  INTEGER, PARAMETER :: ido_heatherz = 5
  INTEGER, PARAMETER :: ido_heathart = 6
  INTEGER, PARAMETER :: ido_heathug = 7
  INTEGER, PARAMETER :: ido_heatchap = 8
  INTEGER, PARAMETER :: ido_heatsw = 9
  INTEGER, PARAMETER :: ido_flxfub = 10
  INTEGER, PARAMETER :: ido_flxup = 11
  INTEGER, PARAMETER :: ido_flxdo = 12
  INTEGER, PARAMETER :: ido_flxhart = 13
  INTEGER, PARAMETER :: ido_flupc = 14
  INTEGER, PARAMETER :: ido_fllya = 15
  INTEGER, PARAMETER :: ido_flhz = 16
  INTEGER, PARAMETER :: ido_flsrc = 17
  INTEGER, PARAMETER :: ido_flsrb = 18
  INTEGER, PARAMETER :: ido_flhuup = 19
  INTEGER, PARAMETER :: ido_flhudo = 20
  INTEGER, PARAMETER :: ido_flchup = 21
  INTEGER, PARAMETER :: ido_flchdo = 22
  INTEGER, PARAMETER :: ido_trfub = 23
  INTEGER, PARAMETER :: ido_trfubc = 24
  INTEGER, PARAMETER :: ido_altswc = 25
  INTEGER, PARAMETER :: ido_altsw = 26
  INTEGER, PARAMETER :: ido_sflux0 = 27
  INTEGER, PARAMETER :: ido_po3c = 28
  INTEGER, PARAMETER :: ido_po2c = 29
  INTEGER, PARAMETER :: ido_heatswc = 30
  INTEGER, PARAMETER :: NMAXOUT=30  ! no. of fubrad output objects
  
  CHARACTER(LEN=10), DIMENSION(NMAXOUT):: setname= &
       (/ &
       'heato3flux', 'heatlya   ', 'heatsrc   ', 'heatsrb   ', &
       'heatherz  ', 'heathart  ', 'heathug   ', 'heatchap  ', &
       'heatsw    ', 'flxfub    ', 'flxup     ', 'flxdo     ', &
       'flxhart   ', 'flupc     ', 'fllya     ', 'flhz      ', &
       'flsrc     ', 'flsrb     ', 'flhuup    ', 'flhudo    ', &
       'flchup    ', 'flchdo    ', 'trfub     ', 'trfubc    ', &
       'altswc    ', 'altsw     ', 'sflux0    ', 'po3c      ', &
       'po2c      ', 'heatswc   '                              &
       /)

  INTEGER, DIMENSION(NMAXOUT):: repr_idx 
  CHARACTER(LEN=16), DIMENSION(NMAXOUT):: repr_name=(/ &
       'GP_3D_MID       ', 'GP_3D_MID       ', 'GP_3D_MID       ', &
       'GP_3D_MID       ', 'GP_3D_MID       ', 'GP_3D_MID       ', &
       'GP_3D_MID       ', 'GP_3D_MID       ', 'GP_3D_MID       ', &
       'GP_3D_INT       ', 'GP_3D_INT       ', 'GP_3D_INT       ', &
       'GP_3D_INT       ', 'GP_3D_INT       ', 'GP_3D_INT       ', &
       'GP_3D_INT       ', 'GP_3D_INT       ', 'GP_3D_INT       ', &
       'GP_3D_INT       ', 'GP_3D_INT       ', 'GP_3D_INT       ', &
       'GP_3D_INT       ', 'GP_3D_INT       ', 'GP_3D_INT       ', &
       'GP_2D_HORIZONTAL', 'GP_2D_HORIZONTAL', 'GP_2D_HORIZONTAL', &
       'GP_3D_INT       ', 'GP_3D_INT       ', 'GP_3D_MID       '  &
      
       /)

  CHARACTER(LEN=58), DIMENSION(NMAXOUT):: LONGNAME=(/&
       'shortwave heating from o3fluxes scheme                   ', &
       'shortwave heating from Lyman alpha band                  ', &
       'shortwave heating from schumann runge continuum          ', &
       'shortwave heating from schumann runge bands              ', &
       'shortwave heating herzberg                               ', &
       'shortwave heating hartley                                ', &
       'shortwave heating huggins                                ', &
       'shortwave heating chappuis                               ', &
       'total shortwave heating                                  ', &
       'shortwave flux fubrad                                    ', &
       'upward shortwave flux                                    ', &
       'downward shortwave flux                                  ', &
       'downward shortwave flux (Hartley)                        ', &
       'clear sky upward flux                                    ', &
       'downward shortwave flux (Lyman alpha)                    ', &
       'downward shortwave flux (Herzberg)                       ', &
       'downward shortwave flux (Schumann-R continuum)           ', &
       'downward shortwave flux (Schumann-R bands)               ', &
       'upward shortwave flux (Huggins)                          ', &
       'downward shortwave flux (Huggins)                        ', &
       'upward shortwave flux (Chappuis)                         ', &
       'downward shortwave flux (Chappuis)                       ', &
       'transmissivity FUBRad                                    ', &
       'transmissivity FUBRad clear sky (experimental)           ', &
       'shortwave albedo clear sky                               ', &
       'shortwave albedo                                         ', &
       'net top solar radiation fubrad                           ', &
       'O3 column in path length                                 ', &
       'O2 column in path length                                 ', &
       'total shortwave heating clear sky                        '  &
       /)

  CHARACTER(LEN=6), DIMENSION(NMAXOUT):: unit=(/&
       'K/s   ', 'K/s   ', 'K/s   ', 'K/s   ', 'K/s   ', &
       'K/s   ', 'K/s   ', 'K/s   ', 'K/s   ', 'W/m**2', &
       'W/m**2', 'W/m**2', 'W/m**2', 'W/m**2', 'W/m**2', &
       'W/m**2', 'W/m**2', 'W/m**2', 'W/m**2', 'W/m**2', &
       'W/m**2', 'W/m**2', '-     ', '-     ', '-     ', &
       '-     ', 'W/m**2', 'atm-cm', 'atm-cm', 'K/s   '  &
       /)

  !extra diagnostics for FUBRad
  TYPE(t_rad_work),   DIMENSION(NMAXOUT,NRADCALL) :: xradout

  ! solar insolation = solar irradiation (at 1 AU) at current orbit position
  REAL(dp), PUBLIC             :: solc_fubrad
  ! SOLAR CYCLE DATA LINEARLY INTERPOLATED IN TIME
  REAL(dp), DIMENSION(:), POINTER, PUBLIC :: solval => NULL()

  ! POINTERS FOR COUPLED CHANNEL OBJECTS
  ! cos zenith angle (instantaneous)
  REAL(dp), POINTER :: amu0_x(:,:)  => NULL()
  
  ! distance sun - eart in AU
  REAL(dp), POINTER :: cdisse => NULL()

  interface rad_fubrad_initialize
    module procedure rad_fubrad_initialize
  end interface rad_fubrad_initialize

  interface rad_fubrad_init_memory
    module procedure rad_fubrad_init_memory
  end interface rad_fubrad_init_memory

  interface rad_fubrad_init_coupling
     module procedure rad_fubrad_init_coupling
  end interface

  interface rad_fubrad_global_start
     module procedure rad_fubrad_global_start
  end interface

  interface rad_fubrad_radheat
    module procedure rad_fubrad_radheat
  end interface rad_fubrad_radheat

  interface rad_fubrad_preprad
    module procedure rad_fubrad_preprad
    module procedure rad_fubrad_preprad_2
  end interface rad_fubrad_preprad

  interface rad_fubrad_free_memory
     module procedure rad_fubrad_free_memory
  end interface

  INTERFACE rad_fubrad_flx2tr
     MODULE PROCEDURE rad_fubrad_flx2tr
  END INTERFACE

  INTERFACE rad_fubrad_fluxes
     MODULE PROCEDURE rad_fubrad_fluxes
  END INTERFACE

  public :: lfubrad
  public :: rad_fubrad_initialize
  public :: rad_fubrad_init_memory
  public :: rad_fubrad_init_coupling
  public :: rad_fubrad_preprad
  public :: rad_fubrad_radheat
  public :: rad_fubrad_free_memory
  public :: rad_fubrad_global_start

  CONTAINS

   ! ========================================================================
   subroutine rad_fubrad_initialize

     USE messy_rad_fubrad_mem,     ONLY: solfac
     USE messy_main_tools,         ONLY: find_next_free_unit

     implicit none

     !-- Local variables
     integer                     :: status, iou
     CHARACTER(LEN=*), PARAMETER :: substr='rad_fubrad_initialize'

     CALL start_message_bi(submodstr,'INITIALIZATION',substr)

     ! mz_pj_20071113+
     !--- Read namelist
     IF (p_parallel_io) THEN
        iou = find_next_free_unit(100,200)
        CALL rad_fubrad_read_nml_ctrl(status, iou)
        IF (status /= 0) CALL error_bi('',substr)
     END IF

     CALL p_bcast(solfac, p_io)
     CALL p_bcast(nbands, p_io)
     CALL p_bcast(sr%type,p_io)
     CALL p_bcast(sr%ntr, p_io)
     CALL p_bcast(resolution, p_io)

     !--- Read namelist
     IF (p_parallel_io) THEN
        iou = find_next_free_unit(100,200)
        CALL rad_fubrad_read_nml_cpl(status, iou)
        IF (status /= 0) CALL error_bi('',substr)
     END IF
     ! Broadcast namelist variables
     call p_bcast (fubrad_solar%cha, p_io)
     call p_bcast (fubrad_solar%obj, p_io)
     call p_bcast (blev, p_io)
     
     ! Determine level indices for Free University of Berlin
     ! shortwave radiation scheme
     CALL set_fublevels

     ! Initialize parameters, that depend on spectral resolution
     SELECT CASE(TRIM(resolution))
     CASE('DYNAMIC')
       CALL fubrad_ini_param_dyn(status) ! op_mk_20170831
     CASE DEFAULT
       CALL fubrad_ini_param(status)
     END SELECT
     IF (status /= 0) CALL error_bi(' error in fubrad_ini_param', substr)

     CALL end_message_bi(submodstr,'INITIALIZATION',substr)

     return
   end subroutine rad_fubrad_initialize
   ! ========================================================================

   ! ========================================================================
   subroutine rad_fubrad_init_memory

     ! BMIL
     USE messy_main_channel_error_bi, ONLY: channel_halt
     USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                          , new_attribute, REPR_UNDEF
     USE messy_main_channel_repr,     ONLY: get_representation_id &
                                          , get_representation_info
     USE messy_main_tools,            ONLY: int2str
     USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM

     IMPLICIT NONE

     ! LOCAL
     CHARACTER(LEN=*), PARAMETER  :: substr = 'rad_fubrad_init_memory'
     INTEGER                      :: status
     INTEGER                      :: i, j2
     INTEGER, DIMENSION(NMAXOUT)  :: irank
     CHARACTER(LEN=2)             :: idx = ''
     CHARACTER(LEN=STRLEN_MEDIUM) :: sname = ''

     CALL start_message_bi(submodstr,'WORKING MEMORY',substr)

     ! ALLOCATE SOME MEMORY
     call fubrad_initialize (nproma, nlev)

     ! Initialize the solar fluxes and absorption cross sections.
     SELECT CASE(TRIM(resolution))
     CASE('DYNAMIC')
       CALL fubrad_initialize_fluxes_dyn(status)
       IF (status /= 0) CALL error_bi(' error in fubrad_initialize_fluxes_dyn',substr)
       CALL fubrad_initialize_cross_sec_dyn(status)
       IF (status /= 0) CALL error_bi(' error in fubrad_initialize_cross_sec_dyn',substr)
     CASE DEFAULT
       CALL fubrad_initialize_fluxes(status)
       IF (status /= 0) CALL error_bi(' error in fubrad_initialize_fluxes',substr)
       CALL fubrad_initialize_cross_sec(status)
       IF (status /= 0) CALL error_bi(' error in fubrad_initialize_cross_sec',substr)
     END SELECT

     SELECT CASE(TRIM(sr%type))
     CASE('KM','KOPPERS-MURTAGH')
       ! initialise Chebyshev polynomial Coeff.
       CALL fubrad_srb_km_init_xs(status)
       IF (status /= 0) CALL error_bi(' error in fubrad_srb_km_init_xs',substr)
     CASE('KCK','KOCKARTS')
       ! initialise te parameters for the Kockarts (1994) Schumann-Runge
       ! parametrization
       CALL fubrad_srb_kck_schu_init(status)
       IF (status /= 0) CALL error_bi(' error in fubrad_srb_kck_schu_init',substr)
     END SELECT

     CALL end_message_bi(submodstr,'WORKING MEMORY',substr)

     CALL start_message_bi(submodstr,'CHANNEL DEFINITION',substr)

     ! get representation indices and ranks
     repr_idx(:) = REPR_UNDEF
     DO i=1, NMAXOUT
        CALL get_representation_id(status, TRIM(repr_name(i)), repr_idx(i))
        CALL get_representation_info(status, inpname='' &
             ,id=repr_idx(i), rank=irank(i))
        CALL channel_halt(substr, status)
     END DO
     
     rad_calls: DO j2=1, NRADCALL
        IF (.NOT. l_switch(j2)) CYCLE    
        
        CALL int2str(idx, j2, '0')
        sname='rad'//idx//'_fubrad'
        
        IF (p_parallel_io) WRITE(*,*) 'add new channel '//sname//'....'
        
        CALL new_channel(status, TRIM(sname),lrestreq=.TRUE.)
        CALL channel_halt(substr, status)
        
        out_objs: DO i=1,NMAXOUT 
           SELECT CASE (irank(i))
           CASE(0)
              CALL new_channel_object(status, TRIM(sname), TRIM(setname(i)) &
                   , p0=xradout(i,j2)%ptr0,reprid=repr_idx(i))
           CASE(2)
              CALL new_channel_object(status, TRIM(sname), TRIM(setname(i)) &
                   , p2=xradout(i,j2)%ptr2,reprid=repr_idx(i))
           CASE(3)
              CALL new_channel_object(status, TRIM(sname), TRIM(setname(i)) &
                   , p3=xradout(i,j2)%ptr3,reprid=repr_idx(i))
           CASE(4)
              CALL new_channel_object(status, TRIM(sname), TRIM(setname(i)) &
                   , p4=xradout(i,j2)%ptr4,reprid=repr_idx(i))
           END SELECT
           CALL channel_halt(substr//'channel object '//setname(i)//&
                &' not found', status)
           
           CALL new_attribute(status, TRIM(sname), TRIM(setname(i)) &
                ,'long_name', c=TRIM(LONGNAME(i)))
           CALL channel_halt(substr, status)
           CALL new_attribute(status, TRIM(sname), TRIM(setname(i)) &
                ,'units', c=TRIM(unit(i)))
           CALL channel_halt(substr, status)
        END DO out_objs
        
     END DO rad_calls
     
     CALL end_message_bi(submodstr,'CHANNEL DEFINITION',substr)
     
     return
   end subroutine rad_fubrad_init_memory
   ! ========================================================================

   ! ========================================================================
   subroutine rad_fubrad_init_coupling

     ! ECHAM5/MESSy
     USE messy_main_mpi_bi,           ONLY: p_parallel_io
     USE messy_main_channel_error_bi, ONLY: channel_halt
     USE messy_main_channel,          ONLY: get_channel_object &
                                          , get_channel_object_dimvar
     USE messy_main_tools,            ONLY: PTR_1D_ARRAY
     USE messy_main_constants_mem,    ONLY: STRLEN_ULONG

     IMPLICIT NONE
     INTRINSIC :: SIZE, TRIM

     ! local
     CHARACTER(LEN=*), PARAMETER   :: substr = 'rad_fubrad_init_coupling'
     INTEGER                       :: status    ! error flag
     TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER          :: dvs => NULL()
     CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: units => NULL()
     INTEGER :: i

     CALL start_message_bi(submodstr, 'COUPLING INITIALIZING', substr)

     ! cos zenith angle (instantaneous)
     CALL get_channel_object(status, 'orbit', 'cosszac', p2=amu0_x)
     CALL channel_halt(substr//': channel object for cosszac not found!',&
          status)

     CALL get_channel_object(status, 'orbit', 'cdisse', p0=cdisse)
     CALL channel_halt(substr//': channel object for cdisse not found!',&
         status)

     external_solar_data: IF (TRIM(fubrad_solar%cha) /= '') THEN

        CALL info_bi('Looking for SOLAR CYCLE DATA ... ')
        CALL info_bi('       channel: '//fubrad_solar%cha)
        CALL info_bi('       object : '//fubrad_solar%obj)

        CALL get_channel_object(status &
             , TRIM(fubrad_solar%cha), TRIM(fubrad_solar%obj), p1=solval)
        CALL channel_halt(substr, status)

        CALL get_channel_object_dimvar(status &
             , TRIM(fubrad_solar%cha), TRIM(fubrad_solar%obj) &
             , dvs, units)
        CALL channel_halt(substr, status)

        DO i=1, SIZE(dvs)
           IF (p_parallel_io) &
                WRITE(*,*) 'DIMVAR ',i,' [',TRIM(units(i)),']: ',dvs(i)%ptr
        END DO

        IF ( (SIZE(solval)-1) /= nbands ) THEN
           CALL error_bi( &
                'Number of parameters does not match number of bands!',substr)
        ENDIF

     ENDIF external_solar_data

     CALL end_message_bi(submodstr, 'COUPLING INITIALIZING', substr)

   end subroutine rad_fubrad_init_coupling
   ! ========================================================================

   ! ========================================================================
   subroutine rad_fubrad_global_start

     USE messy_rad_fubrad_init,   ONLY: fubrad_solar_time_control &
                                      , fubrad_global_flux_ini
     implicit none
     intrinsic :: associated

     character(len=*), parameter :: substr='rad_fubrad_global_start'
     integer :: status

     IF (ASSOCIATED(solval)) THEN
        CALL fubrad_solar_time_control(status, cdisse, solc_fubrad, solval)
     ELSE
        CALL fubrad_solar_time_control(status, cdisse, solc_fubrad)
     END IF

     IF (status /= 0) &
          CALL error_bi(' error in fubrad_solar_time_control',substr)

     ! Initialize the TOA fluxes for current time step
     CALL fubrad_global_flux_ini

   end subroutine rad_fubrad_global_start
   ! ========================================================================

   ! ========================================================================
   subroutine rad_fubrad_radheat(kproma, nproma, jrow, l_trigrad, j2,   &
        nlevp1, zi0,zo3, zo2, heatmosw, heatmoswc, rdayl_x, tm1, ptte,  &
        trnir, trnif, trsol, trsof,    &
        flxs, srad0u, flxsf, flxnir,   &
        flxuni, flxunif, flxus, flxusf,     &
        trupnir, trupnif, trupsol, trupsof  &
        )

     implicit none

     integer, intent(IN)                       :: kproma, nproma, jrow
     logical, intent(IN)                       :: l_trigrad
     integer, intent(IN)                       :: j2
     integer, intent(IN)                       :: nlevp1
     REAL(dp),intent(in),dimension(:)          :: zi0
     REAL(dp),intent(in),dimension(:,:)        :: zo3
     REAL(dp),intent(in),dimension(:,:)        :: zo2
     REAL(dp),intent(in),dimension(:,:)        :: heatmosw
     REAL(dp),intent(in),dimension(:,:)        :: heatmoswc
     REAL(dp),intent(IN),dimension(:)          :: rdayl_x
     REAL(dp),intent(IN),dimension(:,:)        :: tm1
     REAL(dp),intent(OUT),dimension(:,:)       :: ptte
     REAL(dp),intent(IN),dimension(:,:)        :: trnir
     REAL(dp),intent(IN),dimension(:,:)        :: trnif
     REAL(dp),intent(INOUT),dimension(:,:)     :: trsol
     REAL(dp),intent(INOUT),dimension(:,:)     :: trsof
     REAL(dp),intent(INOUT),dimension(:,:)     :: flxs
     REAL(dp),intent(INOUT),dimension(:)       :: srad0u
     REAL(dp),intent(INOUT),dimension(:,:)     :: flxsf
     REAL(dp),intent(IN),dimension(:,:)        :: flxnir

     REAL(dp),intent(IN),dimension(:,:)        :: flxuni
     REAL(dp),intent(IN),dimension(:,:)        :: flxunif
     REAL(dp),intent(INOUT),dimension(:,:)     :: flxus
     REAL(dp),intent(INOUT),dimension(:,:)     :: flxusf
     REAL(dp),intent(IN),dimension(:,:)        :: trupnir
     REAL(dp),intent(IN),dimension(:,:)        :: trupnif
     REAL(dp),intent(INOUT),dimension(:,:)     :: trupsol
     REAL(dp),intent(INOUT),dimension(:,:)     :: trupsof

     ! local
     integer :: iswlev

     !   Actions required, if NOT a radiation time step
     IF(.NOT.l_trigrad) THEN

        !   Get ozone2 files to edit

        !   Between full radiation time steps, 
        !   fluxes of the FUB shortwave radiation
        !   scheme are calculated in Radheat (L_RADTIMESTEP=FALSE).
        !   At these time steps O3-Values are needed in Radheat

        call prepare_radiation (kproma, nlev, zo3, zo2, apm1 &
             , aphm1, amu0_x(:,jrow) )

        call middle_atmosphere_downward_flux (kproma, nproma, nlev, tm1, iswlev)
        call middle_atmosphere_upward_flux   (kproma, nproma, nlev)

     ENDIF

     CALL middle_atmosphere_heat_rates (kproma, nproma, nlev, aphm1, apm1, &
          amu0_x(:,jrow), rdayl_x(:), zo2(:,:), ptte(:,:), &
          xradout(ido_heato3flux, j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_heatsrc,    j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_heatsrb,    j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_heatlya,    j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          heatmosw(:,:),                              &
          heatmoswc(:,:),                             &
          xradout(ido_heatsw,     j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_heatswc,    j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_heatherz,   j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_heathart,   j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_heathug,    j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_heatchap,   j2)%ptr3(_RI_XYZ__(:,jrow,:)), &
          xradout(ido_sflux0,     j2)%ptr2(:,jrow),   &
          delta_time )

     !
     ! Convert the short-wave flux of the FUBRad scheme to 
     ! transmissivity and save in 3-d array.
     ! These values substitute the transmissivity of the UVvis
     ! band of the SW4 scheme.
     !
     CALL rad_fubrad_flx2tr(kproma, nproma, nlevp1        &
          , xradout(ido_trfub, j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
          , xradout(ido_trfubc,j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
          )
     !
     ! Correct the diagnostic arrays of the shortwave fluxes (flxs, flxsf)
     ! and the shortwave transmissivities (trsol, trsof) for FUBRad levels
     ! and save the FUBRad fluxes.
     ! Calculate the diagnostic srad0u at the TOA with the transmissivity of 
     ! FUBRad and the transmissivity of the NIR bands.
     !
     CALL rad_fubrad_fluxes ( kproma, nproma, nlevp1      &
          , amu0_x(:,jrow), rdayl_x(:), zi0(:)     &
          , flxnir(:,:)                            &
          , flxs(:,:)                              &
          , flxsf(:,:)                             &
          , xradout(ido_flxfub, j2)%ptr3(_RI_XYZ__(:,jrow,:)) &
          , trsol(:,:)                             &
          , trsof(:,:)                             &
          , trnir(:,:)                             &
          , trnif(:,:)                             &
          , xradout(ido_trfub, j2)%ptr3(_RI_XYZ__(:,jrow,:))  &
          , xradout(ido_trfubc,j2)%ptr3(_RI_XYZ__(:,jrow,:))  &
          , srad0u(:)                              &
          , flxuni(:,:), flxunif(:,:)              &
          , flxus(:,:), flxusf(:,:)                &
          , trupnir(:,:), trupnif(:,:)             &
          , trupsol(:,:), trupsof(:,:)             &
          ) 

     return
   end subroutine rad_fubrad_radheat
   ! ========================================================================

   ! -------------------------------------------------------------------
   SUBROUTINE rad_fubrad_flx2tr(kproma, kbdim, klevp1, ptrfub, ptrfubc)
     !
     ! Purpose: Convert the short-wave flux of the FUBRad scheme to 
     ! -------- transmissivity and save in 3-d array. These values 
     !          substitute the transmissivity of the UVvis band of 
     !          the SW4 scheme.
     ! 
     ! Author:  Markus Kunze, FUB, Dec. 2011.
     ! -------
     INTEGER, INTENT(in) :: kproma  ! number of local longitudes
     INTEGER, INTENT(in) :: kbdim   ! first dimension of 2-d arrays
     INTEGER, INTENT(in) :: klevp1  ! klev+1
     !
     ! OUTPUT:  ptrfub  - transmissivity FUBRad, total sky
     ! -------  ptrfubc - transmissivity FUBRad, clear sky (experimental)
     !
     REAL(dp), DIMENSION(kbdim,klevp1), INTENT(out) :: ptrfub
     REAL(dp), DIMENSION(kbdim,klevp1), INTENT(out) :: ptrfubc

     INTEGER  :: jl
     REAL(dp) :: zfact
     ! 
     ! Initialize.
     !
     ptrfub (:,NSWLEV+2:) = 0._dp
     ptrfubc(:,NSWLEV+2:) = 0._dp
     !
     zfact = 1._dp/solc_fubrad
     !
     DO jl = 1, kproma
        !
        ! Convert the flux profile to a transmissivity profile.
        !
        ptrfub (jl,1:NSWLEV+1) = ( fldo(jl,:)   + flup(jl,:)    &
             + flhart(jl,:) + fllya(jl,:) ) &
             * zfact
        ptrfubc(jl,1:NSWLEV+1) = ( fldo(jl,:)   + flupc(jl,:)    &
             + flhart(jl,:) + fllya(jl,:) ) &
             * zfact
     END DO
     RETURN
   END SUBROUTINE rad_fubrad_flx2tr
  ! -------------------------------------------------------------------
  !
  ! -------------------------------------------------------------------
  !
  SUBROUTINE rad_fubrad_fluxes(kproma, kbdim, klevp1, pamu0, prdayl &
                                  , pi0                               &
                                  , pflxnir, pflxs, pflxsf            &
                                  , pflxfub                           &
                                  , ptrsol, ptrsof, ptrnir, ptrnif    &
                                  , ptrfub, ptrfubc, psrad0u          &
                                  , pflxuni, pflxunif, pflxus, pflxusf      &
                                  , ptrupnir, ptrupnif, ptrupsol, ptrupsof  &
                                  ) 
    !
    ! Purpose: Create the correct net shortwave flux profiles and
    ! -------- transmissivity profiles for the complete vertical 
    !          model domain for total sky and clear sky separately. 
    !          Save the fluxes of FUBRad in an 3-d array. 
    !          The fluxes are originally stored in two-dimensional arrays:
    !          fldo and flup contain fluxes from the Herzberg-, Huggins-,
    !          and Chappuis-bands, fluxes from the Hartley-bands and 
    !          Lyman-alpha are stored separately.
    ! Background:
    ! -----------
    !          The transmissivity ptrsol(:,:) is set constant for
    !          the UVvis part of the spectrum in the subroutine 
    !          rad_sw_sw1s; to get the right diagnostic (pflxs)
    !          at the FUBRad levels, the fluxes from FUBRad are
    !          combined with the NIR fluxes to substitute the fluxes at all
    !          levels where the FUBRad scheme operates (nswlev).
    ! 
    ! Author:  Markus Kunze, FUB, Dec. 2011.
    ! -------
    INTEGER, INTENT(in) :: kproma  ! number of local longitudes
    INTEGER, INTENT(in) :: kbdim   ! first dimension of 2-d arrays
    INTEGER, INTENT(in) :: klevp1  ! klev+1
    !
    REAL(dp), DIMENSION(kbdim),        INTENT(in) :: pamu0
    REAL(dp), DIMENSION(kbdim),        INTENT(in) :: prdayl
    REAL(dp), DIMENSION(kbdim),        INTENT(in) :: pi0
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in) :: pflxnir
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in) :: ptrnir
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in) :: ptrnif
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in) :: ptrfub
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in) :: ptrfubc
    !
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(inout) :: pflxs
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(inout) :: pflxsf
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(inout) :: ptrsol
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(inout) :: ptrsof
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(out)   :: pflxfub
    REAL(dp), DIMENSION(kbdim),        INTENT(out)   :: psrad0u

    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in)    :: pflxuni
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in)    :: pflxunif
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(inout) :: pflxus
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(inout) :: pflxusf
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in)    :: ptrupnir
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in)    :: ptrupnif
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(inout) :: ptrupsol
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(inout) :: ptrupsof
    !
    INTEGER  :: jl
    REAL(dp) :: zfact
    REAL(dp) :: zfl2tr  ! conversion fact to get transmissivity
    
    zfl2tr = 1._dp/solc_fubrad
    ! 
    ! Initialize with zero at all non FUBRad levels.
    !
    pflxfub(:,NSWLEV+2:) = 0._dp
    !
    DO jl = 1, kproma
       ! 
       ! Save FUBRad upward and downward shortwave fluxes,
       ! and the separate fluxes from the Hartley bands.
       !
       zfact = prdayl(jl) * pamu0(jl) * cdisse_fubrad
       !
       ! Save the flux from FUBRad at the FUBRad levels.
       !
       pflxfub (jl,1:NSWLEV+1) = (fldo(jl,:)  + flup(jl,:) + &
                                  flhart(jl,:) + fllya(jl,:)) &
                                  * zfact
       !
       ! Replace the shortwave flux profiles at FUBRad levels by a 
       ! combination of the NIR flux profile at FUBRad levels and 
       ! the flux from FUBRad.
       !
       pflxs(jl,1:NSWLEV+1) =  pflxnir (jl,1:NSWLEV+1) + &
                               pflxfub (jl,1:NSWLEV+1)

       ! For clear sky conditions.
       pflxsf(jl,1:NSWLEV+1)=  ptrnif(jl,1:NSWLEV+1) * pi0(jl) + &
                               (fldo(jl,:)   + flupc(jl,:) + &
                                flhart(jl,:) + fllya(jl,:)) &
                                * zfact
       ! For upward directed SW flux.
       pflxus(jl,1:NSWLEV+1)  = pflxuni(jl,1:NSWLEV+1) + &
                                flup(jl,:) * zfact
       pflxusf(jl,1:NSWLEV+1) = pflxunif(jl,1:NSWLEV+1) + &
                                flupc(jl,:) * zfact

       !
       ! Replace the transmissivity profiles at FUBRad levels by a 
       ! combination of the NIR transmissivity profile at FUBRad levels and 
       ! the transmissivity from FUBRad.
       !
       ptrsol(jl,1:NSWLEV+1) = ptrnir(jl,1:NSWLEV+1) + &
                               ptrfub(jl,1:NSWLEV+1)
       ptrsof(jl,1:NSWLEV+1) = ptrnif(jl,1:NSWLEV+1) + &
                               ptrfubc(jl,1:NSWLEV+1)

       ptrupsol(jl,1:NSWLEV+1) = ptrupnir(jl,1:NSWLEV+1) + &
                                 flup(jl,:) * zfl2tr
       ptrupsof(jl,1:NSWLEV+1) = ptrupnif(jl,1:NSWLEV+1) + &
                                 flupc(jl,:) * zfl2tr

       fldo  (jl,:) = fldo  (jl,:) * zfact
       flup  (jl,:) = flup  (jl,:) * zfact
       flupc (jl,:) = flupc (jl,:) * zfact
       flhart(jl,:) = flhart(jl,:) * zfact
       fllya (jl,:) = fllya (jl,:) * zfact
       flhudo(jl,:) = flhudo(jl,:) * zfact
       flhuup(jl,:) = flhuup(jl,:) * zfact
       flchdo(jl,:) = flchdo(jl,:) * zfact
       flchup(jl,:) = flchup(jl,:) * zfact
       flhz  (jl,:) = flhz  (jl,:) * zfact
       flsrc (jl,:) = flsrc (jl,:) * zfact
       flsrb (jl,:) = flsrb (jl,:) * zfact
    END DO
    !
    ! Correct the top solar radiation upward, which is
    ! calculated in messy_rad:rad_radheat_smcl, with the
    ! fubrad flux included.
    !
    psrad0u(:) = (pi0(:) * (ptrnir(:,1) + ptrfub(:,1))) - pi0(:)
    !
    RETURN
  END SUBROUTINE rad_fubrad_fluxes
  ! -------------------------------------------------------------------

  ! ========================================================================
  subroutine rad_fubrad_preprad(kproma,klev,jrow,o3,o2,ppf,pph,zmu0,j2)

    USE messy_rad_fubrad_mem,     ONLY: altsw, altswc, po3c, flup, fldo &
                                      , flhart, flupc, fllya            &
                                      , po2c &
                                      , flhuup, flhudo, flchup, flchdo &
                                      , flhz, flsrb, flsrc

    implicit none
    integer,  intent(IN)                 :: kproma,klev,jrow
    REAL(dp), intent(IN),dimension(:,:)  :: o3

    REAL(dp), intent(IN),dimension(:,:)  :: o2
    REAL(dp), intent(IN),dimension(:,:)  :: ppf
    REAL(dp), intent(IN),dimension(:,:)  :: pph
    REAL(dp), intent(IN),dimension(:)    :: zmu0
    integer,  intent(IN)                 :: j2

    altsw  => xradout(ido_altsw,  j2)%ptr2(:,jrow)
    altswc => xradout(ido_altswc, j2)%ptr2(:,jrow)
    po3c   => xradout(ido_po3c,   j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flup   => xradout(ido_flxup,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    fldo   => xradout(ido_flxdo,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flhart => xradout(ido_flxhart,j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flupc  => xradout(ido_flupc,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    fllya  => xradout(ido_fllya,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))

    po2c   => xradout(ido_po2c,   j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flhuup => xradout(ido_flhuup, j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flhudo => xradout(ido_flhudo, j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flchup => xradout(ido_flchup, j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flchdo => xradout(ido_flchdo, j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flhz   => xradout(ido_flhz,   j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flsrb  => xradout(ido_flsrb,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flsrc  => xradout(ido_flsrc,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    
    call prepare_radiation (kproma,klev,o3,o2,ppf,pph,zmu0)

    return
  end subroutine rad_fubrad_preprad
  ! ========================================================================

  ! ========================================================================
  subroutine rad_fubrad_preprad_2(kproma, jrow, j2)

    USE messy_rad_fubrad_mem,     ONLY: altsw, altswc, po3c, fldo  &
                                      , flhart, flupc, fllya, po2c &
                                      , flhuup, flhudo, flchup, flchdo &
                                      , flhz, flsrb, flsrc
    !
    ! prepare for heating rate calculation, of a non-radiation time step
    !
    implicit none
    integer,intent(IN)                         :: kproma, jrow, j2

    altsw  => xradout(ido_altsw,  j2)%ptr2(:,jrow)
    altswc => xradout(ido_altswc, j2)%ptr2(:,jrow)
    po3c   => xradout(ido_po3c,   j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flup   => xradout(ido_flxup,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    fldo   => xradout(ido_flxdo,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flhart => xradout(ido_flxhart,j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flupc  => xradout(ido_flupc,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    fllya  => xradout(ido_fllya,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))

    po2c   => xradout(ido_po2c,   j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flhuup => xradout(ido_flhuup, j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flhudo => xradout(ido_flhudo, j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flchup => xradout(ido_flchup, j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flchdo => xradout(ido_flchdo, j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flhz   => xradout(ido_flhz,   j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flsrb  => xradout(ido_flsrb,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))
    flsrc  => xradout(ido_flsrc,  j2)%ptr3(_RI_XYZ__(:,jrow,1:nswlev+1))

    return
  end subroutine rad_fubrad_preprad_2
  ! ========================================================================

  ! ========================================================================
  subroutine rad_fubrad_read_nml_cpl (status, iou)

    ! MESSy
    USE messy_main_tools,    ONLY: read_nml_open, read_nml_check, read_nml_close

    implicit  none
    intrinsic :: trim

    ! I/O
    INTEGER,          INTENT(OUT) :: status     ! error status
    INTEGER,          INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: modstr = 'rad'
    CHARACTER(LEN=*), PARAMETER :: substr = 'rad_fubrad_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 0

    CALL read_nml_open(lex, substr, iou, 'CPL_FUBRAD', modstr)
    IF (.NOT.lex) then
      status = 1
      RETURN    ! <modstr>.nml does not exist
    end if

    READ(iou, NML=CPL_FUBRAD, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_FUBRAD', modstr)
    IF (fstat /= 0) then
      status = 1
      RETURN  ! error while reading namelist
    end if

    IF (TRIM(fubrad_solar%cha) == '') THEN
       CALL warning_bi(&
            'empty channel name for SOLAR CYCLE DATA (fubrad_solar);'&
            &' constant value r_sol in CTRL will be used', substr)
    ELSE
       CALL info_bi('SOLAR CYCLE DATA channel :'//fubrad_solar%cha)
       IF (TRIM(fubrad_solar%obj) == '') THEN
          CALL info_bi('ERROR: empty channel object name for SOLAR CYCLE DATA'&
               &' (fubrad_solar)')
          RETURN
       ELSE
          CALL info_bi('SOLAR CYCLE DATA object  :'//fubrad_solar%obj)
       END IF
    END IF

    CALL read_nml_close(substr, iou, modstr)

    return
  end subroutine rad_fubrad_read_nml_cpl
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE SET_FUBLEVELS

    ! PURPOSE:
    ! --------
    ! THE FUB RADIATION SCHEME OPERATES ON LEVELS ABOVE 70 hPA.
    ! THIS SUBROUTINE DETERMINES LEVEL INDICES NEDDED TO RUN
    ! THE FUB CODE
    !
    !
    ! AUTHOR:
    ! -------
    ! K. Nissen, 16.8.2005, Free University of Berlin

    USE messy_main_constants_mem,   ONLY : dp
    USE messy_main_grid_def_mem_bi, ONLY : nlev
    USE messy_main_grid_def_bi,     ONLY : ceta, apzero

    IMPLICIT NONE

    INTEGER :: jk
    REAL(dp), DIMENSION (nlev)  :: zpf0
    REAL(dp)                    :: zlev=7000.0_dp

    ! zlev is now Initialized as blev, which can (optionally) be set
    ! via namelist. Default is still 7000.0 Pa
    zlev = blev
    ! ECHAM5 FULL LEVEL GRID
    zpf0(:)=ceta(:)*apzero

    ! SEARCH FOR LEVEL AT OR BELOW 70 hPA
    JK=1
    NSWLEV=0
    LOOP1: DO
       IF(zpf0(jk) > zlev) THEN
          NSWLEV=JK-1
       ELSE
          jk=jk+1
       ENDIF
       IF (NSWLEV.NE.0) EXIT LOOP1
    ENDDO LOOP1

    return
  END SUBROUTINE SET_FUBLEVELS
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE rad_fubrad_free_memory

    IMPLICIT NONE
    
    CALL fubrad_clean_memory
    
  END SUBROUTINE rad_fubrad_free_memory
  ! ========================================================================

#endif
! ***************************************************************************
END MODULE messy_rad_fubrad_si
! ***************************************************************************
