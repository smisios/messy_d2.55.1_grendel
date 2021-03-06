#include "messy_main_ppd_bi.inc"

MODULE messy_airsea_si

  !
  !  MESSy- submodel for Gas transfer at water surface.
  !
  !  AUTHOR:  Pozzer Andrea, MPICH, Oct 2004
  !
  !  MESSy Interface  for the Submodel kernel

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,         &    
                                      mtend_id_q, mtend_id_t
#endif
#ifdef ECHAM5
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY
#endif
  USE messy_main_tools,         ONLY: PTR_2D_ARRAY
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_airsea

  IMPLICIT NONE
  PRIVATE

  INTRINSIC NULL

  TYPE, PUBLIC :: TYP_LINK 
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel  = ''
     CHARACTER(LEN=STRLEN_OBJECT ) :: object   = ''
  END TYPE TYP_LINK 

  ! pointers to external data
  REAL(dp), DIMENSION(:,:), POINTER :: wind10_2d => NULL()

  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: kl_airsea_gp => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: flux_gp => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: dc_airsea_gp => NULL()
  TYPE(PTR_2D_ARRAY), POINTER, SAVE                :: wc_airsea => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: kl_sea => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: kl_air => NULL()
 
 REAL(dp),DIMENSION(:,:), POINTER :: sea_land_gp => NULL()

#ifdef ECHAM5
  TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER, SAVE  :: dc_airsea_lg => NULL()
  TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER, SAVE  :: kl_airsea_lg => NULL()
  TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER, SAVE  :: flux_lg => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: flux_lggp => NULL()
  ! Lagrangian field for grid box area
  REAL(dp), DIMENSION(:),   POINTER :: gboxarea_lg => NULL()
  ! Lagrangian field for surface pressure
  REAL(dp), DIMENSION(:),   POINTER :: zps_lg => NULL()
  ! Lagrangian field for orography 
  REAL(dp), DIMENSION(:),   POINTER :: tmp_seaice => NULL()
  REAL(dp), DIMENSION(:),   POINTER :: tmp_slm => NULL()
  ! Lagrangian field for delta pressure
  REAL(dp), DIMENSION(:),   POINTER :: zdp_lg => NULL()
  ! 2d field for delta pressure transformation
  REAL(dp), DIMENSION(:,:), POINTER :: tmp_dp_gp => NULL()
  ! 2d field for density  transformation
  REAL(dp), DIMENSION(:,:), POINTER :: tmp_zdensair_gp => NULL()
  ! 2d field for density transformation
  REAL(dp), DIMENSION(:,:), POINTER :: tmp_cair_gp => NULL()
  ! 2d field for land-sea orography transformation
  REAL(dp),DIMENSION(:), POINTER :: sea_land_lg     => NULL()

  ! 2d field for henry number transformation
  TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER, SAVE  :: henry_val_lg=> NULL()
#endif
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: henry_val_gp=> NULL()
  ! COUPLING TO CLOUD AND CONVECT SUBMODELS
  REAL(dp), DIMENSION(:,:,:),POINTER :: RAIN_L !large scale precipitation 
  REAL(dp), DIMENSION(:,:,:),POINTER :: RAIN_C !convective precipitation 
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: WATER_CON=> NULL()
  ! COUPLING TO OFFLEM SUBMODEL
  REAL(dp), DIMENSION(:,:),POINTER   :: input_SALT !salinity information 
  ! COUPLING TO ATTILA SUBMODEL
#ifdef ECHAM5
 !Temporary pointer for GP-->LG transformation
  REAL(dp), DIMENSION(:,:),   POINTER           :: tmp_pointer_gp => NULL()
  REAL(dp), DIMENSION(:),     POINTER           :: tmp_pointer_lg => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER           :: tmp_pointer_gp_3d => NULL()
  ! AIR MASS OF CELL
  REAL(dp), POINTER               :: MASS_PARCEL
  ! PRESS_PARCEL(NCELLS) = pressure coordinate of the parcel
  REAL(dp), POINTER, DIMENSION(:) :: PRESS_PARCEL
#endif

  ! CPL-NAMELIST PARAMETERS
  LOGICAL :: L_LG = .FALSE.  ! dry deposition for Lagrangian tracers
  LOGICAL :: L_GP = .FALSE.  ! dry deposition for Gridpoint tracers
  INTEGER :: i_lg_method = 1 ! dry deposition method (LG)

  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: convect_rain(2) = '' !information on the channels
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: large_rain(2) = ''   !information on the channels
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC, SAVE :: salinity(2) = ''   !information on the channels

  INTEGER, PUBLIC  :: nasi = 0                  !actual number of tracer  
  INTEGER, PUBLIC, PARAMETER :: NMAXNTRAC_ASI =100        !total number of tracer interaction available 

  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC, SAVE     :: ASI_NAME(NMAXNTRAC_ASI)      !name of tracer
  REAL(dp),DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: HENRY_A                      ! value_a for henry number
  REAL(dp),DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: HENRY_B                      ! value_b for henry number
  REAL(dp),DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: ALPHA                        ! enanchement factor due to  
  REAL(dp),DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: MOL_VOL                      ! molar volume at normal boiling point
  REAL(dp),DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: S_CONST                      ! Setschenow constant
  REAL(dp),DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: MOL_MASS                     ! molar mass
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: USE_MOL_MASS                 ! output forcing  
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: EFFECT                       ! no effect of submodel on this tracer  
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: OUTPUT                       ! output forcing  
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: SATURATION                   ! saturation calculation   
  REAL(dp),DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: WATER_CON_CONST              ! water concentration constant value  
  TYPE(TYP_LINK),DIMENSION(NMAXNTRAC_ASI),PUBLIC,SAVE :: WATER_CON_CHN           ! water concentration from channel
  LOGICAL ,DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: WATER_CHN_USE                ! use of the channel for w.conc 
  

  REAL(dp),DIMENSION(NMAXNTRAC_ASI), PUBLIC,SAVE :: molweight     ! molar mass of a trcer 

  LOGICAL, PUBLIC, SAVE :: l_asi_lg(NMAXNTRAC_ASI) 
  LOGICAL, PUBLIC, SAVE :: l_asi_gp(NMAXNTRAC_ASI) 
  LOGICAL, PUBLIC, SAVE :: l_nasi_lg_tot 
  LOGICAL, PUBLIC, SAVE :: l_nasi_gp_tot 
  LOGICAL, PUBLIC, SAVE :: l_nasi
  INTEGER, PUBLIC, SAVE :: trac_id_gp(NMAXNTRAC_ASI)     !total number of tracer present in the program
  INTEGER, PUBLIC, SAVE :: trac_id_lg(NMAXNTRAC_ASI)     !total number of tracer present in the program

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! PUBLIC ECHAM-5 INTERFACE ROUTINES TO 'messy_SUBMODELS'
  PUBLIC :: airsea_initialize
  PUBLIC :: airsea_init_memory
  PUBLIC :: airsea_init_coupling
  PUBLIC :: airsea_vdiff
  PUBLIC :: airsea_global_end
  PUBLIC :: airsea_free_memory

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_initialize

    ! AIRSEA MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Pozzer Andrea, MPICH, Oct 2004

 
    ! ECHAM5/MESSy
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC TRIM

    ! DEFINE NCREGRID EVENT TRIGGERs
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='airsea_initialize'
    INTEGER                 :: status
    INTEGER                 :: iou    ! I/O unit
    INTEGER                 :: jt

    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL airsea_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    IF (p_parallel_io) THEN
       IF (param_kw.gt.7.or.param_kw.lt.1) THEN
          write(*,*) 'wrong param_kw value!'
          call error_bi(' ',substr)
       ENDIF
    ENDIF

    ! BROADCAST RESULTS
    CALL p_bcast( param_kw,        p_io)
    CALL p_bcast( l_tendency,      p_io)
    CALL p_bcast( l_whitecap,      p_io)
    CALL p_bcast( l_rain,          p_io)
    CALL p_bcast( l_turb,          p_io)
    CALL p_bcast( l_salt,          p_io)

!  INIT
    DO jt=1, NMAXNTRAC_ASI
          ASI_NAME(jt)  = ''
          HENRY_A(jt)   = 1.0_dp
          HENRY_B(jt)   = 0.0_dp
          ALPHA(jt)     = 1.0_dp
          MOL_VOL(jt)   = 0.3_dp
          USE_MOL_MASS(jt) = .TRUE.
          S_CONST(jt)   = 0.0_dp
          MOL_MASS(jt)  = 0.0_dp
          EFFECT(jt)    = .FALSE.
          OUTPUT(jt)    = .FALSE.
          SATURATION(jt) = .FALSE.
          WATER_CON_CONST(jt) = -1.0_dp
          WATER_CON_CHN(jt)%channel = ''
          WATER_CON_CHN(jt)%object  = ''
    END DO

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL airsea_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
       ! COUNT TRACERS
       DO jt=1, NMAXNTRAC_ASI
          IF (TRIM(ASI_NAME(jt)) == '') CYCLE
          nasi = nasi + 1
          !
          ASI_NAME(nasi)          = TRIM(ASI_NAME(jt))
          HENRY_A(nasi)           = HENRY_A(jt)   
          HENRY_B(nasi)           = HENRY_B(jt)
          IF (ALPHA(jt).EQ.0) THEN 
            ALPHA(nasi)             = 1.0_dp
          ELSE
            ALPHA(nasi)             = ALPHA(jt)
          ENDIF 
          MOL_VOL(nasi)           = MOL_VOL(jt)
          S_CONST(nasi)           = S_CONST(jt)
          MOL_MASS(nasi)          = MOL_MASS(jt)
          USE_MOL_MASS(nasi)      = USE_MOL_MASS(jt)
          OUTPUT(nasi)            = OUTPUT(jt)
          EFFECT(nasi)            = EFFECT(jt)
          SATURATION(nasi)        = SATURATION(jt)
          WATER_CON_CONST(nasi)   = WATER_CON_CONST(jt)
          WATER_CON_CHN(nasi)%channel  = WATER_CON_CHN(jt)%channel
          WATER_CON_CHN(nasi)%object   = WATER_CON_CHN(jt)%object 
       END DO
    write (*,*) " Total number of tracer in airsea_nml = " , nasi 
    END IF


    CALL p_bcast(nasi, p_io)

    ! BROADCAST RESULTS
    CALL p_bcast(convect_rain(1), p_io)
    CALL p_bcast(convect_rain(2), p_io)
    CALL p_bcast(large_rain(1), p_io)
    CALL p_bcast(large_rain(2), p_io)
    CALL p_bcast(salinity(1), p_io)
    CALL p_bcast(salinity(2), p_io)
    CALL p_bcast(l_lg, p_io)
    CALL p_bcast(l_gp, p_io)
    CALL p_bcast(i_lg_method, p_io)

    DO jt=1, nasi
       CALL p_bcast(ASI_NAME(jt), p_io)
       CALL p_bcast(HENRY_A(jt), p_io)
       CALL p_bcast(HENRY_B(jt), p_io)
       CALL p_bcast(ALPHA(jt), p_io)
       CALL p_bcast(MOL_VOL(jt), p_io)
       CALL p_bcast(S_CONST(jt), p_io)
       CALL p_bcast(MOL_MASS(jt), p_io)
       CALL p_bcast(USE_MOL_MASS(jt), p_io)
       CALL p_bcast(OUTPUT(jt), p_io)
       CALL p_bcast(WATER_CON_CONST(jt), p_io)
       CALL p_bcast(WATER_CON_CHN(jt)%channel, p_io)
       CALL p_bcast(WATER_CON_CHN(jt)%object , p_io)
       CALL p_bcast(EFFECT(jt), p_io)
       CALL p_bcast(SATURATION(jt), p_io)
    END DO
    
  END SUBROUTINE airsea_initialize 

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_init_memory
   
    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,     ONLY: ntrac_gp, ti_gp
    USE messy_main_grid_def_mem_bi,   ONLY: ngpblks, nproma
    USE messy_main_blather_bi,        ONLY: info_bi
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_channel_bi,        ONLY: GP_2D_HORIZONTAL
#ifdef ECHAM5
    USE messy_main_channel_bi,    ONLY: LG_ATTILA
    USE messy_main_tracer_mem_bi, ONLY: ntrac_lg, ti_lg, NCELL
#endif
    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute

    IMPLICIT NONE
    INTRINSIC TRIM


    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airsea_init_memory'
    INTEGER :: status
    INTEGER :: jt,jl

    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION ', substr)

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    ! value initialisation
    l_nasi=.FALSE.
    l_nasi_gp_tot=.FALSE.
    l_nasi_lg_tot=.FALSE.
    trac_id_gp(:) = 0.
    trac_id_lg(:) = 0.

    IF (L_GP) THEN
      DO jl=1,nasi
        l_asi_gp(jl)=.FALSE.
      END DO
      DO jt = 1, ntrac_gp
        DO jl = 1,nasi
         IF(ti_gp(jt)%tp%ident%fullname == ASI_NAME(jl)) THEN
              l_asi_gp(jl) =.TRUE.
              trac_id_gp(jl) = jt
              l_nasi_gp_tot = .TRUE.
         END IF
        END DO
      END DO
    END IF

#ifdef ECHAM5
!!#D attila +    
    IF (L_LG) THEN
      DO jl=1,nasi
        l_asi_lg(jl)=.FALSE.
      END DO
      DO jt = 1, ntrac_lg
        DO jl = 1,nasi
         IF(ti_lg(jt)%tp%ident%basename == ASI_NAME(jl)) THEN
              l_asi_lg(jl) =.TRUE.
              trac_id_lg(jl) = jt
              l_nasi_lg_tot = .TRUE.
         END IF
        END DO
      END DO
    END IF
!!#D attila -
#endif

    ! First check if you have tracer for the calculations... if not, why do you have the 
    ! submodel switch on??? It has no effect!
    IF (l_nasi_gp_tot) THEN 
       l_nasi=.TRUE.
    ELSE
       L_GP=.FALSE.
    END IF
    IF (l_nasi_lg_tot) THEN 
       l_nasi=.TRUE.
    ELSE
       L_LG=.FALSE.
    END IF

    IF (.not.l_nasi) CALL info_bi('NO TRACER FOR APPLICATION OF AIR_SEA EXCHANGE! ', substr )

    ! to be allocated anyhow. It is then used for lg calculation as well!
    ALLOCATE(sea_land_gp(nproma,ngpblks))
    ALLOCATE(kl_air(nasi))
    ALLOCATE(kl_sea(nasi))
    ALLOCATE(kl_airsea_gp(nasi))
    ALLOCATE(henry_val_gp(nasi))
    ALLOCATE(wc_airsea)
    ALLOCATE(WATER_CON(nasi))
    ! to be allocated depending on the process calculated!
#ifdef ECHAM5
!!#D attila +
    IF (L_LG) THEN
      ALLOCATE(dc_airsea_lg(nasi))
      ALLOCATE(kl_airsea_lg(nasi))
      ALLOCATE(flux_lg(nasi))
      ALLOCATE(flux_lggp(nasi))
      ALLOCATE(gboxarea_lg(NCELL))
      ALLOCATE(tmp_seaice(NCELL))
      ALLOCATE(tmp_slm(NCELL))
      ALLOCATE(zps_lg(NCELL))
      ALLOCATE(zdp_lg(NCELL))
      ALLOCATE(tmp_dp_gp(nproma,ngpblks))
      ALLOCATE(tmp_zdensair_gp(nproma,ngpblks))
      ALLOCATE(tmp_cair_gp(nproma,ngpblks))
      ALLOCATE(sea_land_lg(NCELL))
      ALLOCATE(henry_val_lg(nasi))
    CALL info_bi('allocated',substr )
    END IF
!!#D attila -
#endif
    IF (L_GP) THEN
    ALLOCATE(flux_gp(nasi))
    ALLOCATE(dc_airsea_gp(nasi))
    END IF

    CALL new_channel(status, modstr, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    IF (L_GP) THEN
    CALL new_channel(status, modstr//'_gp', reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    END IF
#ifdef ECHAM5
!!#D attila +
    IF (L_LG) THEN
    CALL new_channel(status, modstr//'_lg', reprid=LG_ATTILA)
    CALL channel_halt(substr, status)
    END IF
!!#D attila -
#endif
    ! DO TO ! OUTPUT FORCING!
     
    ! OUTPUT COMMON FOR LAGRANGE AND GRIDPOINT
    DO jl=1,nasi 
        IF (l_asi_gp(jl).or.l_asi_lg(jl).or.output(jl)) THEN
        ! we need the folowing calculation if we have any gp tracer, any lg tracer or
        ! if we forcing the output calculation (without tracer)
           CALL new_channel_object(status, modstr &
                , 'henry_val_'//TRIM(ASI_NAME(jl)) &
                , p2=henry_val_gp(jl)%ptr )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr &
                , 'henry_val_'//TRIM(ASI_NAME(jl)) &
                , 'long_name', c=' henry number'//&
                &TRIM(ASI_NAME(jl)) )
           CALL channel_halt(substr, status)
           CALL new_channel_object(status, modstr &
                , 'Kl_'//TRIM(ASI_NAME(jl)) &
                , p2=kl_airsea_gp(jl)%ptr )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr &
                , 'Kl_'//TRIM(ASI_NAME(jl)) &
                , 'long_name', c='Transfer velocity for '//&
                &TRIM(ASI_NAME(jl)) )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr &
                , 'Kl_'//TRIM(ASI_NAME(jl)) &
                , 'units', c='m s-1')
           CALL channel_halt(substr, status)
           CALL new_channel_object(status, modstr &
                , 'Kl_air_'//TRIM(ASI_NAME(jl)) &
                , p2=kl_air(jl)%ptr )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr &
                , 'Kl_air_'//TRIM(ASI_NAME(jl)) &
                , 'long_name', c='Transfer velocity for '//&
                &TRIM(ASI_NAME(jl)) )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr &
                , 'Kl_air_'//TRIM(ASI_NAME(jl)) &
                , 'units', c='m s-1')
           CALL channel_halt(substr, status)
           CALL new_channel_object(status, modstr &
                , 'Kl_sea_'//TRIM(ASI_NAME(jl)) &
                , p2=kl_sea(jl)%ptr )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr &
                , 'Kl_sea_'//TRIM(ASI_NAME(jl)) &
                , 'long_name', c='Transfer velocity for '//&
                &TRIM(ASI_NAME(jl)) )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr &
                , 'Kl_sea_'//TRIM(ASI_NAME(jl)) &
                , 'units', c='m s-1')
           CALL channel_halt(substr, status)
        END IF
#ifdef ECHAM5
!!#D attila +
        IF (l_asi_lg(jl)) THEN
        ! this pointer are needed only for lg calculations....
           CALL info_bi('init LG ',substr)
           CALL new_channel_object(status, modstr//'_lg' &
                , 'henry_val_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , p1=henry_val_lg(jl)%ptr )
           CALL channel_halt(substr, status)
           CALL info_bi('init LG ',substr)
           CALL new_attribute(status, modstr//'_lg' &
                , 'henry_val_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , 'long_name', c=' henry number '//&
                &TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) )
           CALL channel_halt(substr, status)
           CALL new_channel_object(status, modstr//'_lg' &
                , 'Kl_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , p1=kl_airsea_lg(jl)%ptr )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr//'_lg' &
                , 'Kl_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , 'long_name', c='Total transfer velocity '//&
                &TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr//'_lg' &
                , 'Kl_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , 'units', c='m s-1')
           CALL channel_halt(substr, status)
        ENDIF
!!#D attila -
#endif
    END DO

    IF (L_GP) THEN
      DO jl=1,nasi 
          IF (l_asi_gp(jl)) THEN
           CALL new_channel_object(status, modstr//'_gp' &
                , 'flux_'//TRIM(ti_gp(trac_id_gp(jl))%tp%ident%fullname) &
                , p2=flux_gp(jl)%ptr )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr//'_gp' &
                , 'flux_'//TRIM(ti_gp(trac_id_gp(jl))%tp%ident%fullname) &
                , 'units', c='mol / (m^2 s)')
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr//'_gp' &
                , 'flux_'//TRIM(ti_gp(trac_id_gp(jl))%tp%ident%fullname) &
                , 'long_name', c='flux for '//&
                &TRIM(ti_gp(trac_id_gp(jl))%tp%ident%fullname) )
           CALL channel_halt(substr, status)
             CALL new_channel_object(status, modstr//'_gp' &
                  , 'dc_'//TRIM(ti_gp(trac_id_gp(jl))%tp%ident%fullname) &
                  , p2=dc_airsea_gp(jl)%ptr )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_gp' &
                  , 'dc_'//TRIM(ti_gp(trac_id_gp(jl))%tp%ident%fullname) &
                  , 'long_name', c='diff conc for '//&
                  &TRIM(ti_gp(trac_id_gp(jl))%tp%ident%fullname) )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_gp' &
                  , 'dc_'//TRIM(ti_gp(trac_id_gp(jl))%tp%ident%fullname) &
                  , 'units', c='mol m-3')
             CALL channel_halt(substr, status)

#ifdef MESSYTENDENCY
             CALL mtend_register(my_handle, trac_id_gp(jl))
#endif
          END IF
       END DO
      END IF

#ifdef ECHAM5
!!#D attila +
    IF (L_LG) THEN
      DO jl=1,nasi 
          IF (l_asi_lg(jl)) THEN
           CALL new_channel_object(status, modstr//'_lg' &
                , 'flux_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , p1=flux_lg(jl)%ptr )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr//'_lg' &
                , 'flux_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , 'units', c='mol / (m^2 s)')
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr//'_lg' &
                , 'flux_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , 'long_name', c='flux for '//&
                &TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) )
           CALL channel_halt(substr, status)
           CALL new_channel_object(status, modstr//'_lg' &
                , 'flux_lggp_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , p2=flux_lggp(jl)%ptr &
                , reprid=GP_2D_HORIZONTAL)
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr//'_lg' &
                , 'flux_lggp_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , 'units', c='mol / (m^2 s)')
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr//'_lg' &
                , 'flux_lggp_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                , 'long_name', c='flux for '//&
                &TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) )
           CALL channel_halt(substr, status)
             CALL new_channel_object(status, modstr//'_lg' &
                  , 'dc_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                  , p1=dc_airsea_lg(jl)%ptr )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg' &
                  , 'dc_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                  , 'long_name', c='diff conc for '//&
                  &TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_lg' &
                  , 'dc_'//TRIM(ti_lg(trac_id_lg(jl))%tp%ident%fullname) &
                  , 'units', c='mol m-3')
             CALL channel_halt(substr, status)

          END IF
       END DO
      END IF
!!#D attila -
#endif

       CALL new_channel_object(status, modstr, 'wc' &
            , p2=wc_airsea%ptr )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'wc' &
            , 'long_name', c='white cap coverage')
       
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'wc' &
                  , 'units', c='fraction')
       CALL channel_halt(substr, status)

    CALL end_message_bi(modstr, 'MEMORY INITIALIZATION ', substr)

  END SUBROUTINE airsea_init_memory

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_init_coupling

    ! MESSy
    USE messy_main_blather_bi,        ONLY: error_bi, info_bi
    USE messy_main_mpi_bi,            ONLY: p_parallel_io
    USE messy_main_channel_error_bi,  ONLY: channel_halt 
    USE messy_main_channel,           ONLY: get_channel_object &
                                          , get_channel_info  &
                                          , get_channel_object_info
    USE messy_main_tracer_mem_bi,     ONLY: ti_gp
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi,     ONLY: ti_lg
#endif
    USE messy_main_tracer,            ONLY: R_molarmass
    USE messy_main_data_bi,           ONLY: basemodstr=>modstr

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'airsea_init_coupling'
    ! I/O
    INTEGER :: status
    INTEGER :: jl
    INTEGER :: reprid

    CALL start_message_bi(modstr, 'COUPLING', substr)

    CALL get_channel_object(status, TRIM(basemodstr), 'wind10', p2=wind10_2d)
    CALL channel_halt(substr, status)

    IF (l_salt) THEN
       CALL info_bi('IMPORT SALINITY ', substr)
       CALL get_channel_info(status, salinity(1))
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
              WRITE(*,*) ' no salinity climatology available from other submodels '
              WRITE(*,*) ' channel not present ! (check namelist)'
              WRITE(*,*) ' average salinity calculation (0.4 mol/L) '
          END IF
          l_salt=.false.
       ELSE
          CALL get_channel_object_info(status          &
               , cname=salinity(1), oname=salinity(2)  & 
               , reprid=reprid )
          IF (status /= 0) THEN
             IF (p_parallel_io) THEN
                 WRITE(*,*) ' no salinity climatology available from other submodels '
                 WRITE(*,*) ' object not present ! (check namelist)'
                 WRITE(*,*) ' average salinity calculation (0.4 mol/L) '
             END IF
             l_salt=.false.
          ELSE
             IF (p_parallel_io) THEN
                 WRITE(*,*) ' salinity climatology available from channel ', salinity(2) 
             END IF
             CALL get_channel_object(status, cname=salinity(1), oname=salinity(2) &
                  , p2=input_SALT )
             CALL channel_halt(substr, status)
          END IF
       END IF
    END IF

    IF (l_rain) THEN
       CALL info_bi('import large scale precipitation from other submodel ', substr)
       IF (TRIM(large_rain(1)) == '') THEN
         CALL info_bi('no other submodel specified in the namelist for large scale rain', substr) 
         CALL info_bi(' l_rain set to FALSE',substr)  
         l_rain = .FALSE.
       ELSE 
         CALL get_channel_info(status, large_rain(1))
         IF (status /= 0) THEN
            IF (p_parallel_io) THEN
                WRITE(*,*) ' no large scale precipitation available from other submodels '
                WRITE(*,*) ' channel not present ! (check namelist)'
                WRITE(*,*) ' l_rain set to FALSE'
            END IF
            l_rain = .FALSE.
         ELSE
            CALL get_channel_object_info(status          &
                 , cname=large_rain(1), oname=large_rain(2)  & 
                 , reprid=reprid )
            IF (status /= 0) THEN
               IF (p_parallel_io) THEN
                   WRITE(*,*) ' no large scale precipitation available from other submodels '
                   WRITE(*,*) ' object not present ! (check namelist)'
                   WRITE(*,*) ' l_rain set to FALSE'
               END IF
               l_rain=.false.
            ELSE
               IF (p_parallel_io) THEN
                   WRITE(*,*) ' large scale precipitation available from channel ', large_rain(2) 
               END IF
               CALL get_channel_object(status, cname=large_rain(1), oname=large_rain(2) &
                          , p3=RAIN_L )
               CALL channel_halt(substr, status)
            END IF
         END IF
       END IF

       CALL info_bi('import convective precipitation from other submodel ',substr)
       IF (TRIM(convect_rain(1)) == '') THEN
         CALL info_bi('no other submodel specified in the namelist for convective rain',substr) 
         CALL info_bi(' l_rain set to FALSE',substr)  
         l_rain = .FALSE.
       ELSE 
         CALL get_channel_info(status, convect_rain(1))
         IF (status /= 0) THEN
            IF (p_parallel_io) THEN
                WRITE(*,*) ' no convective precipitation available from other submodels '
                WRITE(*,*) ' channel not present ! (check namelist)'
                WRITE(*,*) ' l_rain set to FALSE'
            END IF
            l_rain = .FALSE.
         ELSE
            CALL get_channel_object_info(status          &
                 , cname=convect_rain(1), oname=convect_rain(2)  & 
                 , reprid=reprid )
            IF (status /= 0) THEN
               IF (p_parallel_io) THEN
                   WRITE(*,*) ' no convective precipitation available from other submodels '
                   WRITE(*,*) ' object not present ! (check namelist)'
                   WRITE(*,*) ' l_rain set to FALSE'
               END IF
               l_rain=.false.
            ELSE
               IF (p_parallel_io) THEN
                   WRITE(*,*) ' convective precipitation available from channel ', convect_rain(2) 
               END IF
               CALL get_channel_object(status, cname=convect_rain(1), oname=convect_rain(2) &
                          , p3=RAIN_C )
               CALL channel_halt(substr, status)
            END IF
         END IF
       END IF
    END IF

#ifdef ECHAM5
!!#D attila +
    ! coupling with attila for informations about the parcels
    IF (L_LG) THEN
       CALL info_bi('import information from attila submodel ',substr)
       CALL get_channel_object(status, cname='attila', oname='AMCELL' &
                  , p0=MASS_PARCEL )
       CALL channel_halt(substr, status)
       CALL info_bi('import information from attila submodel ',substr)
       CALL get_channel_object(status, cname='attila', oname='PPRESS' &
                  , p1=PRESS_PARCEL )
       CALL channel_halt(substr, status)
    END IF
!!#D attila -
#endif
    CALL info_bi('IMPORT  WATER CONCENTRATIONS ',substr)
    DO jl=1,nasi
       IF (p_parallel_io) THEN
             WRITE(*,*) ' TRACER       : ',TRIM(ASI_NAME(jl))
       END IF
       IF (TRIM(WATER_CON_CHN(jl)%channel) == '') THEN
           IF (p_parallel_io) THEN
              WRITE(*,*) ' no channel present for ',TRIM(ASI_NAME(jl)) 
              WRITE(*,*) '  use constant value from namelist '
           END IF
           CALL info_bi(' use constant value from namelist',substr)  
           IF (WATER_CON_CONST(jl) .eq. -1._dp) THEN
              CALL error_bi('1) no channel neither constant values define for this tracer!',substr)
           ELSE
              water_chn_use(jl) = .FALSE.
           END IF
       ELSE
       CALL get_channel_info(status, WATER_CON_CHN(jl)%channel)
          IF (status /= 0) THEN
             IF (p_parallel_io) THEN
               WRITE(*,*) ' no channel ',TRIM(WATER_CON_CHN(jl)%channel),' present for ',TRIM(ASI_NAME(jl))
               WRITE(*,*) '  use constant value from namelist '
             END IF
             IF (WATER_CON_CONST(jl) .eq. -1._dp) THEN
                CALL error_bi('2) no channel neither constant values define for this tracer!',substr)
             ELSE
                 water_chn_use(jl) = .FALSE.
             END IF
          ELSE
             CALL get_channel_object_info(status &
               , TRIM(WATER_CON_CHN(jl)%channel), TRIM(WATER_CON_CHN(jl)%object) &
               , reprid=reprid )
             IF (status /= 0) THEN
                IF (p_parallel_io) THEN
                  WRITE(*,*) ' no object(',TRIM(WATER_CON_CHN(jl)%object) ,') present for ',TRIM(ASI_NAME(jl)), &
                             'in channel ', TRIM(WATER_CON_CHN(jl)%channel)
                  WRITE(*,*) '  use constant value from namelist '
                END IF
                IF (WATER_CON_CONST(jl) .eq. -1._dp) THEN
                   CALL error_bi('3) no channel neither constant values define for this tracer!',substr)
                ELSE 
                   water_chn_use(jl) = .FALSE.
                END IF
             ELSE
                IF (p_parallel_io) THEN
                    WRITE(*,*) ' Water concentration from channel/object found for ',TRIM(ASI_NAME(jl)) 
                END IF
                CALL get_channel_object(status, cname=TRIM(WATER_CON_CHN(jl)%channel), &
                            oname=TRIM(WATER_CON_CHN(jl)%object) &
                           , p2=WATER_CON(jl)%ptr )
                CALL channel_halt(substr, status)
                water_chn_use(jl) = .TRUE.
             END IF
          END IF
       END IF
    END DO
    
    DO jl=1,nasi
       ! we use the molar mass defined in airsea
       IF (output(jl)) THEN
        molweight(jl)=MOL_MASS(jl)
        IF (molweight(jl).lt.1._dp) THEN 
           CALL info_bi(' molar mass not define in airsea namelist!',substr)       
        END IF
       END IF
    ! or we use the one defined in mecca if the tracer are presents....!
       IF (l_asi_gp(jl)) THEN 
        IF (USE_MOL_MASS(jl)) THEN
           molweight(jl)=ti_gp(trac_id_gp(jl))%tp%meta%cask_r(R_molarmass)
              IF (molweight(jl).lt.1._dp) THEN 
                CALL info_bi('problem with molar mass of one tracer in the tracer definition!',substr) 
                CALL info_bi('.... using defined molar mass in airsea..',substr)
                molweight(jl)=MOL_MASS(jl)
                IF (molweight(jl).lt.1._dp) THEN 
                   CALL error_bi(' molar mass not define in airsea namelist!',substr)       
                END IF
              END IF
        END IF
       END IF
#ifdef ECHAM5
!!#D attila +
     IF (L_LG) THEN
         IF (l_asi_lg(jl)) THEN 
             IF (USE_MOL_MASS(jl)) THEN
                molweight(jl)=ti_lg(trac_id_lg(jl))%tp%meta%cask_r(R_molarmass)
                   IF (molweight(jl).lt.1._dp) THEN 
                      CALL info_bi('problem with molar mass of one tracer in the tracer definition!',substr) 
                      CALL info_bi('.... using defined molar mass in airsea..',substr)
                      molweight(jl)=MOL_MASS(jl)
                      IF (molweight(jl).lt.1._dp) THEN 
                         CALL error_bi(' molar mass not define in airsea namelist!',substr)       
                      END IF
                   END IF
             END IF
         END IF
     END IF
!!#D attila -
#endif
    END DO

    CALL end_message_bi(modstr, 'COUPLING', substr)

  END SUBROUTINE airsea_init_coupling
  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_vdiff
 
    ! ECHAM5/MESSy
    USE messy_main_blather_bi,      ONLY: info_bi
    USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma, nproma
    USE messy_main_grid_def_bi,     ONLY: deltaz
    USE messy_main_data_bi,         ONLY: slf, seaice, alake                   &
                                        , press_3d, zust_2d                    &
                                        , pressi_3d, tsurf_2d, um1, vm1        &
                                        , cdnw, cfmw, cfncw, riw, tvw, az0, tvir
#ifdef ECHAM5
    USE messy_main_data_bi,       ONLY: pxtems
#endif
    USE messy_main_timer,         ONLY: time_step_len, lstart  
#ifndef MESSYTENDENCY
    USE messy_main_data_bi,       ONLY: tm1, qm1, tte_3d, qte_3d 
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1 
#endif
    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='airsea_vdiff'
    INTEGER :: jl,idt
    LOGICAL          :: check=.FALSE.               ! logical for stability checking 
    REAL(dp), TARGET :: conc_water_gp(nproma,nasi)  ! mixing ratio water
    REAL(dp),TARGET  :: zdp(nproma)                 ! pressure difference for a box
#ifdef ECHAM5
    REAL(dp), POINTER :: zxtems(:,:)                ! flux changed in vdiff 
#endif
    REAL(dp), POINTER :: wind10_mod(:)              ! wind speed at 10 meters
    REAL(dp) :: schmidt_air(nproma,nasi)            ! schmidt number (adimensional) in air 
    REAL(dp) :: schmidt_sea(nproma,nasi)            ! schmidt number (adimensional) in sea
    REAL(dp) :: flux_airsea_gp(nproma,nasi)         ! flux (mol(trac) m-2 s-1)
    REAL(dp) :: densflux_airsea_gp(nproma,nasi)     ! flux (mol(trac) kg(air) mol-1(air) m-2 s-1)
    REAL(dp) :: zflux_airsea(nproma,nasi)           ! flux (mol(trac) mol-1(air) s-1)
    REAL(dp) :: zdensair_gp(nproma)                 ! air density!!! [kg/m3]
    REAL(dp) :: temp(nproma), sphum(nproma)         ! temperature & humidity!
    REAL(dp) :: conc_air_gp(nproma, nasi)           ! mixing ratio air, GP 
    REAL(dp) :: zdz(nproma)                         ! height of the lowest level 
    REAL(dp) :: rain(nproma)                        ! total rain 
    REAL(dp) :: salt(nproma)                        ! salinity 
    REAL(dp) :: ztmst                               ! timestep
    REAL(dp) :: cair(nproma)                        ! no lev--> we look only at surface: air concentration
#ifdef MESSYTENDENCY
    REAL(dp) :: vstart(nproma,nlev)
    REAL(dp) :: zxtte(kproma,nlev)
#endif

!==========================================================================================================           
!         INITIALISATION
!==========================================================================================================           

    ztmst=time_step_len
#ifdef ECHAM5
    zxtems => pxtems(:,1,:,jrow)
#endif
    wind10_mod => wind10_2d(:,jrow)

    zdp(:) = 0._dp

#ifndef MESSYTENDENCY
    sphum(1:kproma) = qm1(_RI_XYZ__(1:kproma,jrow,nlev)) + &
                      qte_3d(_RI_XYZ__(1:kproma,jrow,nlev)) * ztmst
    temp(1:kproma)  = tm1(_RI_XYZ__(1:kproma,jrow,nlev)) + &
                      tte_3d(_RI_XYZ__(1:kproma,jrow,nlev)) * ztmst
#else
    CALL mtend_get_start_l(mtend_id_q, v0=vstart)
    sphum(1:kproma) = vstart(1:kproma, nlev)
    CALL mtend_get_start_l(mtend_id_t, v0=vstart)
    temp(1:kproma)  = vstart(1:kproma, nlev)
#endif

    zdp(1:kproma)   = pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) - &
                      pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev))
    check=.FALSE. 

    IF (l_rain) THEN
       !from mm/s to mm/h
       ! total rain
       rain(1:kproma) = 3600*(rain_l(_RI_XYZ__(1:kproma,jrow,nlev))+rain_c(_RI_XYZ__(1:kproma,jrow,nlev)))
    END IF
    IF (l_salt) THEN
       !from ppt(part per thousand) to M (mol/L)
       salt(1:kproma) = input_salt(1:kproma,jrow)/58.4
    ELSE
       salt(1:kproma) = 0.4_dp !it's in mol/L
    END IF

    DO jl=1,nasi
       IF (l_asi_gp(jl)) THEN 
#ifndef MESSYTENDENCY
        idt = trac_id_gp(jl)
        conc_air_gp(1:kproma,jl) = pxtm1(_RI_X_ZN_(1:kproma,nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,nlev,idt))*ztmst
#else
        CALL mtend_get_start_l( trac_id_gp(jl), v0=vstart)
        conc_air_gp(1:kproma,jl) = vstart(1:kproma, nlev)
#endif
       END IF
    END DO
!==========================================================================================================           
!         COMMON PART
!==========================================================================================================           

!--- 1. MAIN PART: calculation of the parametrisation 
!----We avoid calculation in the first time step, due to not definition of wind speed and temperature
!     if (MINVAL(temp(1:kproma)).gt.0._dp.and.MINVAL(zust_2d(1:kproma,jrow)).gt.0._dp) then
     if (.not.lstart) then
     
     ! land cover ---> orography :
     call airsea_oro(sea_land_gp(1:kproma,jrow), slf(1:kproma,jrow),seaice(1:kproma,jrow), alake(1:kproma,jrow))

     zdz(1:kproma) = deltaz(_RI_XYZ__(1:kproma,jrow,nlev))

     CALL airsea_calc_density(temp(1:kproma),  &
                     sphum(1:kproma), press_3d(_RI_XYZ__(1:kproma,jrow,nlev)),zdensair_gp(1:kproma))
  
     CALL airsea_calc_cair(temp(1:kproma),sphum(1:kproma),press_3d(_RI_XYZ__(1:kproma,jrow,nlev)),cair(1:kproma))
   
     CALL airsea_calc_wc(wind10_mod(1:kproma),wc_airsea%ptr(1:kproma,jrow),sea_land_gp(1:kproma,jrow),kproma)

!==========================================================================================================           
!            LOOP FOR THE TRACER DEFINED IN THE KERNEL
!==========================================================================================================           
        DO jl=1,nasi 

!==========================================================================================================           
!           GENERAL CALCULATION ALWAYS NEEDED  
!==========================================================================================================           

            IF (l_asi_gp(jl).or.l_asi_lg(jl).or.output(jl)) THEN

             call airsea_calc_henry(henry_a(jl),henry_b(jl),tsurf_2d(1:kproma,jrow),henry_val_gp(jl)%ptr(1:kproma,jrow), &
                                    s_const(jl),mol_vol(jl),salt(1:kproma))  
             call airsea_calc_schmidt_sea(tsurf_2d(1:kproma,jrow),schmidt_sea(1:kproma,jl),MOL_VOL(jl)) 
             call airsea_calc_schmidt_air(temp(1:kproma),schmidt_air(1:kproma,jl),MOL_VOL(jl),molweight(jl),                    &
                                          press_3d(_RI_XYZ__(1:kproma,jrow,nlev))) 
             IF (l_whitecap) THEN
               call airsea_calc_kl_sea_wc(wind10_mod(1:kproma),schmidt_sea(1:kproma,jl),kl_sea(jl)%ptr(1:kproma,jrow), &
                       kproma, sea_land_gp(1:kproma,jrow),wc_airsea%ptr(1:kproma,jrow),henry_val_gp(jl)%ptr(1:kproma,jrow),&
                       tsurf_2d(1:kproma,jrow))
             ELSE
               call airsea_calc_kl_sea(wind10_mod(1:kproma),schmidt_sea(1:kproma,jl), kl_sea(jl)%ptr(1:kproma,jrow),   &
                       kproma, sea_land_gp(1:kproma,jrow))
             END IF
             IF (l_rain) THEN
                call airsea_calc_kl_rain(kl_sea(jl)%ptr(1:kproma,jrow),kproma,schmidt_sea(1:kproma,jl),                &
                                         sea_land_gp(1:kproma,jrow),rain(1:kproma))
             END IF
             IF (l_turb) THEN
                 call airsea_calc_kl_air_special(kl_air(jl)%ptr(1:kproma,jrow),molweight(jl),    &
                                     cdnw(1:kproma,jrow),cfmw(1:kproma,jrow),cfncw(1:kproma,jrow), riw(1:kproma,jrow),  &
                                     tvir(1:kproma,jrow), tvw(1:kproma,jrow),g,az0(1:kproma,jrow),zdz(1:kproma),        &
                                     um1(_RI_XYZ__(1:kproma,jrow,nlev)),vm1(_RI_XYZ__(1:kproma,jrow,nlev)),kproma, sea_land_gp(1:kproma,jrow))  
             ELSE
                 call airsea_calc_kl_air(kl_air(jl)%ptr(1:kproma,jrow),zust_2d(1:kproma,jrow),                       &
                                         um1(_RI_XYZ__(1:kproma,jrow,nlev)),vm1(_RI_XYZ__(1:kproma,jrow,nlev)),schmidt_air(1:kproma,jl),   &
                                         kproma, sea_land_gp(1:kproma,jrow))  
             ENDIF
             call airsea_calc_kl_tot(kl_sea(jl)%ptr(1:kproma,jrow),kl_air(jl)%ptr(1:kproma,jrow),                      &
                                    kl_airsea_gp(jl)%ptr(1:kproma,jrow),ALPHA(jl),                                     &
                              henry_val_gp(jl)%ptr(1:kproma,jrow), tsurf_2d(1:kproma,jrow),kproma,                     & 
                              seaice(1:kproma,jrow), sea_land_gp(1:kproma,jrow))
                              
             IF (water_chn_use(jl)) THEN
               conc_water_gp(1:kproma,jl)=WATER_CON(jl)%ptr(1:kproma,jrow)
             ELSE
              call airsea_water_conc(conc_water_gp(1:kproma,jl), kproma, water_con_const(jl))    
             END IF
 
             call check_stability(kl_airsea_gp(jl)%ptr(1:kproma,jrow),zdz(1:kproma),ztmst,kproma,check)

             IF (check) then  
                  CALL info_bi('airsea transfer velocity too big! ---> process instable !',substr) 
             END IF

          END IF

!==========================================================================================================           
!           GRID POINT PART
!==========================================================================================================           

          IF (l_asi_gp(jl)) THEN
           call airsea_delta_conc(conc_air_gp(1:kproma,jl),henry_val_gp(jl)%ptr(1:kproma,jrow),conc_water_gp(1:kproma,jl), &
                        pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)), dc_airsea_gp(jl)%ptr(1:kproma,jrow),sea_land_gp(1:kproma,jrow), &
                                 kproma,EFFECT(jl),saturation(jl))
           call airsea_flux(dc_airsea_gp(jl)%ptr(1:kproma,jrow),kl_airsea_gp(jl)%ptr(1:kproma,jrow),                 &
                            flux_airsea_gp(1:kproma,jl))
           flux_gp(jl)%ptr(1:kproma,jrow)=flux_airsea_gp(1:kproma,jl) 
           ! mol(trac) m-2 s-1 * Kg m-3 * m3 mol(air)-1 --> mol(trac) mol(air)-1 Kg m-2 s-1
           densflux_airsea_gp(1:kproma,jl) = flux_airsea_gp(1:kproma,jl)*zdensair_gp(1:kproma)/cair(1:kproma) 
           ! m3/mol(air) * m-1 * mol(trac) m-2 s-1 --> mol(trac) mol(air) s-1
           ! zflux_airsea = densflux_airsea_gp * (g/dp)
           ! these three following calculation sould give the same results!!!!
!           zflux_airsea(1:kproma,jl) = (1._dp/(cair(1:kproma) * zdz(1:kproma))  &
!                                         * flux_airsea_gp(1:kproma,jl))
          zflux_airsea(1:kproma,jl) = densflux_airsea_gp(1:kproma,jl) * (g/zdp(1:kproma)) 
!          zflux_airsea(1:kproma,jl) = flux_airsea_gp(1:kproma,jl)*zdensair_gp(1:kproma)/cair(1:kproma) &
!                                      * (g/zdp(1:kproma)) 
!==========================================================================================================           
!             UPGRADE CONCENTRATION.... ONLY IN GRID POINT!
!==========================================================================================================           
#ifdef ECHAM5
           IF (.NOT. l_tendency ) THEN
             zxtems(1:kproma,trac_id_gp(jl)) = zxtems(1:kproma,trac_id_gp(jl)) + &
                                                    densflux_airsea_gp(1:kproma,jl)
           ELSE  
#endif
#ifndef MESSYTENDENCY
           ! zflux_airsea = densflux_airsea_gp * (g/dp)
             idt = trac_id_gp(jl)
             pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) = pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
                                                  zflux_airsea(1:kproma,jl)
#else
             zxtte(:,:) = 0._dp
             zxtte(1:kproma,nlev) = zflux_airsea(1:kproma,jl)
             CALL mtend_add_l(my_handle, trac_id_gp(jl), px=zxtte)
#endif
#ifdef ECHAM5
           END IF
#endif
          END IF
!==========================================================================================================           
!             CLOSURE OF THE TRACER CYCLE.... 
!==========================================================================================================           
        END DO

!==========================================================================================================           
!             SAVE SOME VALUES NEEDED LATER FOR LG
!==========================================================================================================           
#ifdef ECHAM5
!!#D attila +
          IF (L_LG) THEN
            tmp_dp_gp(1:kproma,jrow) = zdp(1:kproma) 
            tmp_zdensair_gp(1:kproma,jrow) = zdensair_gp(1:kproma) 
            tmp_cair_gp(1:kproma,jrow) = cair(1:kproma) 
          END IF
!!#D attila -
#endif
     end if ! end of time_step checking

  END SUBROUTINE airsea_vdiff

  ! ---------------------------------------------------------------------------


  SUBROUTINE airsea_global_end
#ifdef ECHAM5
!!#D attila +
    USE messy_main_blather_bi,      ONLY: error_bi
    USE messy_main_timer,           ONLY: time_step_len, lstart
    USE messy_main_data_bi,         ONLY: pressi_3d
    USE messy_main_grid_def_bi,     ONLY: gboxarea_2d
    USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks, nlev
    USE messy_main_tracer_mem_bi, ONLY: pxtte_a => qxtte_a, pxtm1_a => qxtm1_a &
                                      , NCELL          
    USE messy_attila_tools_e5,    ONLY: gp2lg_e5, lg2gp_e5, LG2GP_AVE
    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE

    ! FOR LAGRANGIAN
    LOGICAL, DIMENSION(:), POINTER :: liplev => NULL()

    CHARACTER(LEN=*), PARAMETER :: substr='airsea_global_end'
    INTEGER :: jl,jn

    REAL(dp) :: ztmst                               ! timestep
    REAL(dp) :: conc_air_lg(NCELL, nasi)            ! mixing ratio air, LG
    REAL(dp), TARGET :: conc_water_lg(NCELL, nasi)  ! mixing ratio water
    REAL(dp) :: densflux_airsea_lg(NCELL,nasi)      ! flux (mol(trac) kg(air) mol-1(air) m-2 s-1)
    REAL(dp) :: flux_airsea_lg(NCELL,nasi)          ! flux (mol(trac) m-2 s-1)
    REAL(dp), TARGET :: zdensair_lg(NCELL)          ! air density!!! [kg/m3]
    REAL(dp), TARGET :: cair_lg(NCELL)              ! no lev--> we look only at surface: air concentration
                                                    ! cair => [mol(air)/m3]

     if (.not.lstart) then

     lagrangian : IF (L_LG) THEN

        ztmst=time_step_len
 
        ! orography transformation
        CALL gp2lg_e5( sea_land_gp,   &
              sea_land_lg)

        ! density of air parcel transformation
        tmp_pointer_lg => zdensair_lg(:)
        CALL gp2lg_e5( tmp_zdensair_gp,         & 
                       tmp_pointer_lg)
        tmp_pointer_lg => cair_lg(:)
        CALL gp2lg_e5( tmp_cair_gp,         & 
                       tmp_pointer_lg)
          ! car_lg(:)  = all zeroes !!!

!==========================================================================================================           
!           LAGRANGIAN PART
!==========================================================================================================           

        DO jl=1,nasi 
          IF (l_asi_lg(jl)) THEN

           conc_air_lg(1:NCELL,jl) = pxtm1_a(1:NCELL,trac_id_lg(jl)) + pxtte_a(1:NCELL,trac_id_lg(jl))  * ztmst

            ! water concentration transformation
           IF (water_chn_use(jl)) THEN
             tmp_pointer_lg => conc_water_lg(:,jl) 
             CALL gp2lg_e5( WATER_CON(jl)%ptr,   &
                  tmp_pointer_lg)
           ELSE
             conc_water_lg(:,jl)= water_con_const(jl)    
           END IF

            ! henry value transformation
            tmp_pointer_gp => henry_val_gp(jl)%ptr 
            tmp_pointer_lg => henry_val_lg(jl)%ptr 
            CALL gp2lg_e5( tmp_pointer_gp,   &
                  tmp_pointer_lg,llev=liplev)

            ! transfer velocity transformation
            CALL gp2lg_e5( kl_airsea_gp(jl)%ptr,   &
                  kl_airsea_lg(jl)%ptr)

!==========================================================================================================           
!             CALCULATION OF  CONCENTRATION DIFFERENCE AND FLUX.... (molx / m^2 s)
!==========================================================================================================           

            tmp_pointer_gp => pressi_3d(:,nlev+1,:)
            CALL gp2lg_e5(tmp_pointer_gp, zps_lg)

            SELECT CASE(i_lg_method)
            CASE(1)
              CALL gp2lg_e5( tmp_dp_gp, zdp_lg)
              call airsea_delta_conc(conc_air_lg(1:NCELL,jl),henry_val_lg(jl)%ptr(1:NCELL),conc_water_lg(1:NCELL,jl),     & 
                                  zps_lg(:), dc_airsea_lg(jl)%ptr(1:NCELL),sea_land_lg(1:NCELL),NCELL,EFFECT(jl),           &
                                  saturation(jl),liplev )      
            CASE(2)
              CALL gp2lg_e5( gboxarea_2d, gboxarea_lg)
              call airsea_delta_conc(conc_air_lg(1:NCELL,jl),henry_val_lg(jl)%ptr(1:NCELL),conc_water_lg(1:NCELL,jl),     & 
                                  zps_lg(:), dc_airsea_lg(jl)%ptr(1:NCELL),sea_land_lg(1:NCELL),NCELL,EFFECT(jl),           &
                                  saturation(jl),liplev ) 

            CASE(3)
              ! here we suppose that any parcel has the same density
              ! of the grid box where they are located
              zdp_lg(:)=(MAX((-1._dp*(PRESS_PARCEL(:) - zps_lg(:))), 1E-20))      
              call airsea_delta_conc(conc_air_lg(1:NCELL,jl),henry_val_lg(jl)%ptr(1:NCELL),conc_water_lg(1:NCELL,jl),     & 
                                 zps_lg(:), dc_airsea_lg(jl)%ptr(1:NCELL),sea_land_lg(1:NCELL),NCELL,EFFECT(jl),      &
                                 saturation(jl),liplev)
            CASE DEFAULT
              CALL error_bi(          &
                   'NO VALID METHOD SELECTED FOR LAGRANGIAN'//&
                   &' AIR SEA EXCHANGE!',substr)
            END SELECT

            call airsea_flux(dc_airsea_lg(jl)%ptr(1:NCELL),kl_airsea_lg(jl)%ptr(1:NCELL),flux_airsea_lg(1:NCELL,jl))
            flux_lg(jl)%ptr(1:NCELL)=flux_airsea_lg(1:NCELL,jl) 
            ! from flux molx /m?? s to:
            ! mol(trac) m-2 s-1 * Kg m-3 * m3 mol(air)-1 --> mol(trac) mol(air)-1 Kg m-2 s-1
            DO jn=1,NCELL
               IF (liplev(jn)) THEN
               densflux_airsea_lg(jn,jl) = flux_airsea_lg(jn,jl)*zdensair_lg(jn)/cair_lg(jn) 
               ELSE
               densflux_airsea_lg(jn,jl) = 0.0_dp 
               END IF
            END DO

!==========================================================================================================           
!             UPGRADE CONCENTRATION.... ONLY IN LAGRANGIAN... (AS MADE IN MESSY_DRYDEP_E5.F90)  
!==========================================================================================================           

                   cell_loop: DO jn=1,NCELL
                      is_level: IF (liplev(jn)) THEN
                         SELECT CASE(i_lg_method)
                         CASE (2)
                            pxtte_a(jn,trac_id_lg(jl)) = pxtte_a(jn,trac_id_lg(jl))          &
                                 + densflux_airsea_lg(jn,jl)             &
                                 *gboxarea_lg(jn)/MASS_PARCEL 
                         CASE DEFAULT !case 1-3
                            pxtte_a(jn,trac_id_lg(jl)) = pxtte_a(jn,trac_id_lg(jl))          &
                                 + densflux_airsea_lg(jn,jl)             &
                                 /(zdp_lg(jn)/g)
                         END SELECT

                      END IF is_level
                   END DO cell_loop

!==========================================================================================================           
!             UPGRADE OUTPUT.... ONLY IN LAGRANGIAN... (AS MADE IN MESSY_DRYDEP_E5.F90)  
!==========================================================================================================           


                   IF (i_lg_method .ne. 4) THEN
                      ! final output for analysis is done here.
                      ! outside cell loop but inside tracer loop!
                          ALLOCATE (tmp_pointer_gp_3d(1:nproma,1:nlev,1:ngpblks))
                          tmp_pointer_lg => flux_lg(jl)%ptr(1:NCELL)
                          CALL lg2gp_e5(tmp_pointer_lg,tmp_pointer_gp_3d, LG2GP_AVE, fill_value=0._dp)
                          flux_lggp(jl)%ptr(1:nproma,1:ngpblks)=tmp_pointer_gp_3d(1:nproma,nlev,1:ngpblks)
                          DEALLOCATE  (tmp_pointer_gp_3d)
                          NULLIFY(tmp_pointer_gp_3d)
                   END IF
!==========================================================================================================           
!             CLOSURE OF THE TRACER CYCLE.... 
!==========================================================================================================           
          END IF
        END DO


     END IF lagrangian

     end if ! end of time_step checking

!!#D attila -
#endif
  END SUBROUTINE airsea_global_end

  ! ---------------------------------------------------------------------------

  SUBROUTINE airsea_free_memory

    IMPLICIT NONE
    INTRINSIC ASSOCIATED


   IF (ASSOCIATED(kl_airsea_gp))    DEALLOCATE(kl_airsea_gp)
   IF (ASSOCIATED(kl_air))          DEALLOCATE(kl_air)
   IF (ASSOCIATED(kl_sea))          DEALLOCATE(kl_sea)
   IF (ASSOCIATED(dc_airsea_gp))    DEALLOCATE(dc_airsea_gp)
   IF (ASSOCIATED(wc_airsea))       DEALLOCATE(wc_airsea)
   IF (ASSOCIATED(WATER_CON))       DEALLOCATE(WATER_CON)
#ifdef ECHAM5
   IF (ASSOCIATED(kl_airsea_lg))    DEALLOCATE(kl_airsea_lg)
   IF (ASSOCIATED(flux_lg))         DEALLOCATE(flux_lg)
   IF (ASSOCIATED(flux_gp))         DEALLOCATE(flux_gp)
   IF (ASSOCIATED(flux_lggp))       DEALLOCATE(flux_lggp)
   IF (ASSOCIATED(dc_airsea_lg))    DEALLOCATE(dc_airsea_lg)
   IF (ASSOCIATED(gboxarea_lg))       DEALLOCATE(gboxarea_lg)
   IF (ASSOCIATED(zps_lg))            DEALLOCATE(zps_lg)
   IF (ASSOCIATED(zdp_lg))            DEALLOCATE(zdp_lg)
   IF (ASSOCIATED(tmp_dp_gp))         DEALLOCATE(tmp_dp_gp)
   IF (ASSOCIATED(tmp_zdensair_gp))   DEALLOCATE(tmp_zdensair_gp)
   IF (ASSOCIATED(tmp_cair_gp))       DEALLOCATE(tmp_cair_gp)
   IF (ASSOCIATED(sea_land_gp))       DEALLOCATE(sea_land_gp)
   IF (ASSOCIATED(sea_land_lg))       DEALLOCATE(sea_land_lg)
   IF (ASSOCIATED(tmp_slm))           DEALLOCATE(tmp_slm)
   IF (ASSOCIATED(tmp_seaice))        DEALLOCATE(tmp_seaice)
   IF (ASSOCIATED(henry_val_lg))      DEALLOCATE(henry_val_lg)
#endif
   IF (ASSOCIATED(henry_val_gp))      DEALLOCATE(henry_val_gp)
  END SUBROUTINE airsea_free_memory

! ===========================================================================
! PRIVATE ROUTINES
! ===========================================================================

  SUBROUTINE airsea_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Andrea Pozzer, MPICH, Aug 2003

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi, ONLY: NGCELL
#endif
    USE messy_main_blather_bi,    ONLY: error_bi

    IMPLICIT NONE

    INTRINSIC TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    ! switch for skipping calculation of Lagrangian 
    ! rate coefficients..it is local,not broadcasted

    NAMELIST /CPL/ convect_rain,large_rain, salinity, L_GP, L_LG, i_lg_method,     &
            ASI_NAME, HENRY_A, HENRY_B, ALPHA, MOL_VOL,S_CONST,                    &
            EFFECT, OUTPUT, WATER_CON_CONST, WATER_CON_CHN, MOL_MASS,              &
            USE_MOL_MASS, SATURATION
         

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='airsea_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

#ifdef ECHAM5
    IF ((L_LG) .AND. (NGCELL > 0)) THEN
       L_LG = .TRUE.
    ELSE
       IF (L_LG) THEN
!!#D attila +
         WRITE(*,*) 'L_LG = T in namelist'
         WRITE(*,*) 'However no Lagrangian scheme activated ...'
         WRITE(*,*) ' ... setting L_LG = F'
!!#D attila -
       END IF
       L_LG = .FALSE.
    ENDIF
#endif
    IF (L_GP) THEN
       WRITE(*,*) 'AIRSEA IN GRIDPOINT SPACE : ON'
    ELSE
       WRITE(*,*) 'AIRSEA IN GRIDPOINT SPACE : OFF'
    END IF

#ifdef ECHAM5
!!#D attila +
    IF (L_LG) THEN
       WRITE(*,*) 'AIRSEA IN LAGRANGIAN SPACE : ON'
       WRITE(*,*) 'AIRSEA METHOD FOR LG       : ', i_lg_method 
    ELSE
       WRITE(*,*) 'AIRSEA IN LAGRANGIAN SPACE : OFF'
    END IF
!!#D attila -
#endif

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

    IF (.NOT.L_GP .AND. .NOT.L_LG)  &
       CALL error_bi(               &
       '  NEITHER GRID POINT NOR LAGRANGIAN REPRESENTATION'//&
       &' REQUESTED (namelist CPL): No calculation possible!'&
       ,substr)

   WRITE(*,*) '.........................................................'
   WRITE(*,*) '           AIR SEA  TRACER DEFINED                       '
   WRITE(*,*) '.........................................................'

    DO jt=1, NMAXNTRAC_ASI
       IF (TRIM(ASI_NAME(jt)) == '') CYCLE

       WRITE(*,*) '  TRACER NO.          ',jt
       WRITE(*,*) '  NAME              = ', TRIM(ASI_NAME(jt))
       WRITE(*,*) '  HENRY_A           = ', HENRY_A(jt)
       WRITE(*,*) '  HENRY_B           = ', HENRY_B(jt)
       WRITE(*,*) '  ALPHA             = ', ALPHA(jt)
       WRITE(*,*) '  MOL_VOL           = ', MOL_VOL(jt)
       WRITE(*,*) '  S_CONST           = ', S_CONST(jt)
       WRITE(*,*) '  MOLAR_MASS        = ', MOL_MASS(jt)
       WRITE(*,*) '  USE_MOLAR_MASS    = ', USE_MOL_MASS(jt)
       WRITE(*,*) '  FORCING OUTPUT    = ', OUTPUT(jt)
       WRITE(*,*) '  WATER CONC.CONST  = ', WATER_CON_CONST(jt)
       !IF (TRIM(WATER_CON_CHN(jt)%channel) .not. '') THEN
       WRITE(*,*) '  WATER CONC.CHN    = ', WATER_CON_CHN(jt)%channel
       WRITE(*,*) '  WATER CONC.OBJ    = ', WATER_CON_CHN(jt)%object
       !END IF
       WRITE(*,*) '  EFFECT ON AIR     = ', EFFECT(jt)
       WRITE(*,*) '  SATURATION        = ', SATURATION(jt)
       WRITE(*,*) '.........................................................'

    END DO 

  END SUBROUTINE airsea_read_nml_cpl

! ===========================================================================
 
END MODULE messy_airsea_si
