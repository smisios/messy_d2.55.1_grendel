#include "messy_main_ppd_bi.inc"

! **********************************************************************
MODULE messy_plumegas_si
! **********************************************************************

  ! MESSy submodel interface layer (SMIL) for submodel PLUMEGAS
  !
  ! PLUMEGAS
  !      This MESSy (Joeckel et al. 2005) submodel for the ECHAM/MESSy
  !      Atmospheric Chemistry (EMAC) model (Joeckel et al. 2006,
  !      Roeckner et al. 2006) aims at applying the Cariolle et al.
  !      (2009) parameterisation of plume chemistry to NOx emissions
  !      from e.g. combustion engines employed in aviation and shipping
  !      (Huszar et al. 2009). NOx is emitted into corresponding plume
  !      NOx tracers (instead of being directly emitted into a
  !      background tracer) as specified through submodel OFFLEM
  !      (Kerkweg et al. 2006) - see the note in the plumegas.nml
  !      namelist file. This MESSy SMIL module exchanges information
  !      with EMAC. The equations of the parameterisation are in the
  !      submodel core layer (SMCL) module. Subroutine and function
  !      names at the SMIL begin with plumegas_ while subroutine and
  !      function names at the SMCL begin with plg_ .
  !
  ! Authors
  !      Mauro Dall'Amico, DLR Oberpfaffenhofen, Germany, 2009-2010
  !      Patrick Joeckel, DLR Oberpfaffenhofen, Germany, 2009-2010
  !
  ! Bug reports
  !      Deutsches Zentrum fuer Luft- und Raumfahrt, Institut fuer
  !      Physik der Atmosphaere, ECHAM-Gruppe, Oberpfaffenhofen,
  !      Germany
  !
  ! Acknowledgements
  !      Financial support for developing submodel PLUMEGAS was
  !      provided by the Helmholtz Young Investigators Group SeaKLIM
  !      led by Veronika Eyring at DLR, Oberpfaffenhofen, Germany.
  !
  ! References
  !      Cariolle et al. (2009) "Parametrisation of plume chemistry
  !           into large-scale atmospheric models: Application to
  !           aircraft NOx emissions", J. Geophys. Res., 114, D19302,
  !           doi:10.1029/2009JD011873.
  !      Huszar et al. (2009) "Modeling the regional impact of ship
  !           emissions on NOx and ozone levels over the Eastern
  !           Atlantic and Western Europe using ship plume
  !           parameterization", Atmos. Chem. Phys. Discuss., 9,
  !           26735-26776.
  !      Kerkweg et al. (2006) "Technical note: Implementation of
  !           prescribed (OFFLEM), calculated (ONLEM), and
  !           pseudo-emissions (TNUDGE) of chemical species in the
  !           Modular Earth Submodel System (MESSy)", Atmos. Chem.
  !           Phys., 6, 3603-3609
  !      Joeckel et al. (2006) "The atmospheric chemistry general
  !           circulation model ECHAM5/MESSy1: consistent simulation of
  !           ozone from the surface to the mesosphere", Atmos. Chem.
  !           Phys., 6, 5067-5104.
  !      Roeckner et al. (2006) "Sensitivity of simulated climate to
  !           horizontal and vertical resolution in the ECHAM5
  !           atmosphere model", J. Climate, 19, 3771-3791
  !      Joeckel et al. (2005) Technical Note: The Modular Earth
  !           Submodel System (MESSy) - a new approach towards Earth
  !           System Modeling, Atmos. Chem. Phys., 5, 433-444
  !
  ! Important note
  !      This version is only for internal use and test at DLR,
  !      Oberpfaffenhofen, Germany
  !
  ! Version
  !      1.0 - DRAFT (ENTWURF)

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather_bi,    ONLY: start_message_bi &
                                    , end_message_bi
  USE messy_plumegas
  ! Only for the SMCL it is reasonable not to use ONLY .

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! Global CPL ('coupling') namelist variables:

  ! Do not modify tendencies of external tracers if T (default F)
  LOGICAL             :: l_decouple=.FALSE. !         CPL parameter

  ! Assign NOx increments to NO and NO2 (default 'gridbox'):
  CHARACTER(LEN=7) :: as_NOx_NO='gridbox' !           CPL parameter
  ! Values corresponding to as_NOx_NO:
  REAL(dp),DIMENSION(max_ty) :: as_NO  ! Depends on as_NOx_NO
  REAL(dp),DIMENSION(max_ty) :: as_NO2 ! Same as as_NO
  ! *'gridbox' assigns according to the gridbox NO/NOx molar ratio and
  !   the following logical is set to .TRUE. in plumegas_initialize:
  LOGICAL :: l_grb_NO=.FALSE. ! Depends on as_NOx_NO
  ! * 'exhaust' assigns according to the emission ratio emi_ra (CTRL);
  ! * 'only_NO' assigns everything to NO.

  ! The NO/NOx molar ratio in the gridbox:
  REAL(dp),DIMENSION(:,:),ALLOCATABLE :: grb_NO_NOx ! CPL parameter

  ! Global variables:

  ! Indixes for tracers to be passed on to PLUMEGAS:
  INTEGER  :: idx_NO=0
  INTEGER  :: idx_HNO3=0
  INTEGER  :: idx_NO2=0
  INTEGER  :: idx_O3=0

  ! Array for indixes of plume NOx tracers (depending on type):
  INTEGER,DIMENSION(max_ty) :: idx_plNOx

  ! Pointers for channel objects of other submodels:
  REAL(dp),POINTER,DIMENSION(:,:,:) :: ptr_grmass => NULL()
  REAL(dp),POINTER,DIMENSION(:,:,:) :: ptr_grvol  => NULL()

  ! Subroutines and functions:
  PUBLIC :: plumegas_initialize     ! Initialisation, read namelists
  PUBLIC :: plumegas_new_tracer     ! Definition of new tracers
  PUBLIC :: plumegas_init_memory    ! channel definition and local
  !                                   memory allocation (non-channel)

  ! The online tracer initialisation takes place here.
  PUBLIC :: plumegas_init_coupling  ! Coupling initialisation

  ! The time loop starts here.

  ! The regional loop starts here.

  ! Entry points in radiation, vdiff, radheat and convect should be
  ! used only if absolutely required (in which case also the
  ! corresponding lines in messy_main_control_e5 should not be
  ! commented out).

  ! The appropriate location for entry points is in physc:
  PUBLIC :: plumegas_physc          ! Physics & chemistry computations

  ! The regional loop ends here.

  ! The time loop ends here.

  PUBLIC :: plumegas_free_memory    ! Local memory deallocation 
                                    !              (non-channel)

  ! The following line is commented since PRIVATE has been set above:
  !PRIVATE :: plumegas_read_nml_cpl ! Read CPL namelist

CONTAINS

! #####################################################################
! Public subroutines
! #####################################################################

! ---------------------------------------------------------------------
! Initialisation phase
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plumegas_initialize

    ! plumegas_initialize
    !      Initialisation of PLUMEGAS. Read in important parameters
    !      from all namelists and broadcast settings.


    ! BMIL:

    USE messy_main_mpi_bi,     ONLY: p_parallel_io &
                                   , p_io          &
                                   , p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! Local variables and parameters:
    CHARACTER(LEN=*),PARAMETER :: substr='plumegas_initialize'
    INTEGER                    :: iou    ! Input / output unit
    INTEGER                    :: status ! Error status

    ! Initialise CTRL namelist:

    IF (p_parallel_io) THEN
       iou=find_next_free_unit(100,200)
       CALL plg_read_nml_ctrl(status,iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! Broadcast variables from CTRL:
    CALL p_bcast(name_pl_ty,p_io)
    CALL p_bcast(tau,       p_io)
    CALL p_bcast(beta_1,    p_io)
    CALL p_bcast(beta_0,    p_io)
    CALL p_bcast(emi_ra,    p_io)
    CALL p_bcast(keff,      p_io)

    ! Initialise coupling namelist:
    IF (p_parallel_io) THEN
       iou=find_next_free_unit(100,200)
       CALL plumegas_read_nml_cpl(status,iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! Broadcast values from CPL:
    CALL p_bcast(l_decouple,p_io)
    CALL p_bcast(as_NOx_NO, p_io)

    ! The coefficients for assigning the NOX increments to the
    ! tendencies of NO and NO2 are set in the following lines according
    ! to the setting of as_NOx_NO :

    ! Initialisation:
    as_NO(:) =0._dp
    as_NO2(:)=0._dp

    SELECT CASE (as_NOx_NO)

    CASE ('gridbox')
      l_grb_NO=.TRUE.

    CASE ('exhaust')
      as_NO(:) =1._dp-emi_ra(:)
      as_NO2(:)=emi_ra(:)

    CASE ('only_NO')
      as_NO(:) =1._dp
      as_NO2(:)=0._dp

    END SELECT

    ! Set l_pl_ty , a logical array for the used (.TRUE.) and void
    ! (FALSE) plume types:
    CALL plg_set_logical()

    idx_plNOx(:)=0 ! Is usually zero for undefined tracers.

  END SUBROUTINE plumegas_initialize
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plumegas_new_tracer

    ! plumegas_new_tracer
    !    This subroutine defines tracers specific to PLUMEGAS.

    USE messy_main_blather_bi,      ONLY: info_bi
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR     ! The grid point
                                                   ! tracer set
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt ! Error handling
    USE messy_main_constants_mem,   ONLY: dp            &
                                        , MN            & ! Molar masses
                                        , MO
    USE messy_main_tracer,          ONLY: new_tracer     &
                                        , get_tracer     &
                                        , set_tracer     &
                                        , SINGLE         &
                                        , AIR            &
                                        , AMOUNTFRACTION &
                                        , R_MOLARMASS    &
                                        , R_PSS          &
                                        , R_DRYREAC_SF

    IMPLICIT NONE
    INTRINSIC :: TRIM

    CHARACTER(LEN=*),PARAMETER :: substr='plumegas_new_tracer'
    INTEGER :: status
    INTEGER :: ity=0

    CALL info_bi('Defining tracers for PLUMEGAS.', substr)

    ! Define the new plume-NOx tracer(s):

    all_ty: DO ity=1,max_ty
      IF (.NOT. l_pl_ty(ity)) CYCLE

      CALL new_tracer(status                                &
                     ,GPTRSTR                               &
                     ,'plume'                               &
                     ,modstr                                &
                     ,subname='NOx_'//TRIM(name_pl_ty(ity)) &
                     ,idx=idx_plNOx(ity)                    &
                     ,unit='mol/mol'                        &
                     ,type    =SINGLE                       &
                     ,medium  =AIR                          &
                     ,quantity=AMOUNTFRACTION               &
                     )
      CALL tracer_halt(substr,status)

      CALL set_tracer(status             &
                     ,GPTRSTR            &
                     ,idx_plNOx(ity)     &
                     ,R_MOLARMASS        &
                     ,MN+MO              & ! +2._dp*MO for NO2
                     )
      CALL tracer_halt(substr,status)

      CALL set_tracer(status         &
                     ,GPTRSTR        &
                     ,idx_plNOx(ity) &
                     ,R_PSS          &
                     ,2.E-3_dp       & ! 1.0E-2_dp for NO2
                     )
      CALL tracer_halt(substr,status)

      CALL set_tracer(status         &
                     ,GPTRSTR        &
                     ,idx_plNOx(ity) &
                     ,R_DRYREAC_SF   &
                     ,1._dp          &
                     )
      CALL tracer_halt(substr,status)

    ENDDO all_ty

  END SUBROUTINE plumegas_new_tracer
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plumegas_init_memory

    ! plumegas_init_memory
    !    This subroutine defines the channel for this submodel
    !    and allocates the local (non-channel) memory. In particular,
    !    the increments for the tendencies computed by PLUMEGAS are
    !    defined here as channel objects and referenced to via a
    !    pointer to channel objects.

    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, nlev

    IMPLICIT NONE
    INTRINSIC :: TRIM

    CHARACTER(LEN=*), PARAMETER  :: substr='plumegas_init_memory'
    CHARACTER(LEN=STRLEN_MEDIUM) :: lname
    INTEGER :: ity=0
    INTEGER :: status

    CALL start_message_bi(modstr             &
                        ,'Channel definition' & ! As in end_message_bi
                        ,substr              &
                         )

    ! Define a new channel for PLUMEGAS:
    ! The name is 'plumegas' as usual.
    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    all_ty_str_ele: DO ity=1,max_ty
      IF (.NOT. l_pl_ty(ity)) CYCLE

      ! Long names for new channel objects:
      lname='plg_incr_tend_mole_frac_'//                &
            TRIM(name_pl_ty(ity))//'_plume_NOx_'//'air'
      CALL new_channel_object(status, modstr    &
           , 'in_plNOx_'//TRIM(name_pl_ty(ity)) &
           , p3 = in_plNOx(ity)%PTR )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'in_plNOx_'//TRIM(name_pl_ty(ity))   &
           , 'long_name', c=TRIM(lname) )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'in_plNOx_'//TRIM(name_pl_ty(ity))   &
           , 'units', c='mol/mol' )
      CALL channel_halt(substr, status)

      lname='plg_'//TRIM(name_pl_ty(ity))//'_plume_NOx_'// &
            'incr_tend_mole_frac_NO_air'
      CALL new_channel_object(status, modstr    &
           , 'in_NOx_'//TRIM(name_pl_ty(ity)) &
           , p3 = in_NOx(ity)%PTR )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'in_NOx_'//TRIM(name_pl_ty(ity)) &
           , 'long_name', c=TRIM(lname) )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'in_NOx_'//TRIM(name_pl_ty(ity)) & 
           , 'units', c='mol/mol' )
      CALL channel_halt(substr, status)

      lname='plg_'//TRIM(name_pl_ty(ity))//'_'// &
            'incr_tend_mole_frac_HNO3_air'
      CALL new_channel_object(status, modstr    &
           , 'in_HNO3_'//TRIM(name_pl_ty(ity))  &
           , p3 = in_HNO3(ity)%PTR )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'in_HNO3_'//TRIM(name_pl_ty(ity)) &
           , 'long_name', c=TRIM(lname) )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'in_HNO3_'//TRIM(name_pl_ty(ity)) &
           , 'units', c='mol/mol' )
      CALL channel_halt(substr, status)

      lname='plg_'//TRIM(name_pl_ty(ity))//'_'// &
            'incr_tend_mole_frac_O3_air_titr'
      CALL new_channel_object(status, modstr    &
           , 'in_O3_'//TRIM(name_pl_ty(ity))//'_ti' &
           , p3 = in_O3_titr(ity)%PTR )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr &
           , 'in_O3_'//TRIM(name_pl_ty(ity))//'_ti' &
           , 'long_name', c=TRIM(lname) )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr &
           , 'in_O3_'//TRIM(name_pl_ty(ity))//'_ti' &
           , 'units', c='mol/mol' )
      CALL channel_halt(substr, status)


      lname='plg_'//TRIM(name_pl_ty(ity))//'_'// &
            'incr_tend_mole_frac_O3_air_keff'
      CALL new_channel_object(status, modstr    &
           , 'in_O3_'//TRIM(name_pl_ty(ity))//'_ke' &
           , p3 = in_O3_keff(ity)%PTR )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr &
           , 'in_O3_'//TRIM(name_pl_ty(ity))//'_ke' &
           , 'long_name', c=TRIM(lname) )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr &
           , 'in_O3_'//TRIM(name_pl_ty(ity))//'_ke' &
           , 'units', c='mol/mol' )
      CALL channel_halt(substr, status)

    ENDDO all_ty_str_ele

    CALL end_message_bi(modstr              &
                       ,'Channel definition' &
                       ,substr              &
                       )

    ! Module specific non-channel memory should be allocated here:

    ALLOCATE(grb_NO_NOx(1:nproma,1:nlev))
    grb_NO_NOx(:,:)=0._dp

    all_ty_allocate: DO ity=1,max_ty
      IF (.NOT. l_pl_ty(ity)) CYCLE

      ALLOCATE(plNOx_0(ity)%PTR(1:nproma,1:nlev))
      plNOx_0(ity)%PTR(:,:)=0._dp

    ENDDO all_ty_allocate

    ALLOCATE(NO_0( 1:nproma,1:nlev) &
            ,NO2_0(1:nproma,1:nlev) &
            ,O3_0( 1:nproma,1:nlev) &
            )
    NO_0(:,:) =0._dp
    NO2_0(:,:)=0._dp
    O3_0(:,:) =0._dp

    ALLOCATE(molecules_m3(1:nproma,1:nlev))
    molecules_m3(:,:)=0._dp

  END SUBROUTINE plumegas_init_memory
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plumegas_init_coupling

    ! plumegas_init_coupling: 'coupling' to online channels.
    !    This subroutine manages the 'coupling' to online channels from
    !    other submodules.

    USE messy_main_channel,          ONLY: get_channel_object
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: info_bi
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer,           ONLY: get_tracer

    IMPLICIT NONE


    CHARACTER(LEN=*),PARAMETER :: substr='plumegas_init_coupling'
    INTEGER                    :: status ! Error status

    CALL start_message_bi(modstr                &
                         ,'Looking for tracers' &
                         ,substr                &
                         )

    CALL info_bi('Getting NO, NO2, HNO3 and O3 tracers.', substr)

    ! Use indexes for the tracers used in the parametrisation:

    CALL get_tracer(status     &
                   ,GPTRSTR    & ! Tracer set
                   ,'NO'       & ! Basename of the tracer
                   ,idx=idx_NO & ! Index of the tracer
                   )
    CALL tracer_halt(substr,status)

    CALL get_tracer(status      &
                   ,GPTRSTR     &
                   ,'NO2'       &
                   ,idx=idx_NO2 &
                   )
    CALL tracer_halt(substr,status)

    CALL get_tracer(status       &
                   ,GPTRSTR      &
                   ,'HNO3'       &
                   ,idx=idx_HNO3 &
                   )
    CALL tracer_halt(substr,status)

    CALL get_tracer(status     &
                   ,GPTRSTR    &
                   ,'O3'       &
                   ,idx=idx_O3 &
                   )
    CALL tracer_halt(substr,status)

    CALL end_message_bi(modstr                &
                       ,'Looking for tracers' &
                       ,substr                &
                       )

    CALL start_message_bi(modstr                    &
                         ,'Coupling initialisation' &
                         ,substr                    &
                         )

    CALL get_channel_object(status, 'orbit', 'rdayl', p2=ptr_day)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'grid_def', 'grmass', p3=ptr_grmass)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'grid_def', 'grvol', p3=ptr_grvol)
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr                    &
                       ,'Coupling initialisation' &
                       ,substr                    &
                       )

  END SUBROUTINE plumegas_init_coupling
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Time integration phase
! ---------------------------------------------------------------------

! Time loop start

! Regional loop start

! ---------------------------------------------------------------------
  SUBROUTINE plumegas_physc

    ! plumegas_physc
    !    This subroutine acts as an interface to the model's code
    !    describing physical and chemical processes. It (can) e.g.
    !    transfer(s) required fields from/to ECHAM5 for the needs of
    !    submodel PLUMEGAS. This subroutine calls subroutines included
    !    at the submodel core layer,i.e. messy_plumegas.f90 . As a
    !    convention, submodel core layer subroutines have a name
    !    beginning with plg_ while submodel interface layer subroutines
    !    use plumegas_ as prefix.

    ! BMIL:

    USE messy_main_constants_mem,ONLY: strlen_long &
                                     , M_air       &
                                     , N_A

    USE messy_main_grid_def_mem_bi, ONLY: jrow          & ! Current vector
                                     , nlev          & ! No of levels
                                     , kproma
    USE messy_main_timer,        ONLY: time_step_len

    ! The current set of (tracer tendencies of all) tracers, qxt
    ! (qxtte) is an array of size (1:nproma,1:nlev,number of tracers)
    ! where it is the tracer index:
    USE messy_main_tracer_mem_bi,ONLY: qxtm1 &
                                     , qxtte   ! Tendencies

    IMPLICIT NONE

    CHARACTER(LEN=*),PARAMETER :: substr='plumegas_physc'
    !INTEGER                    :: status ! Error status
    !CHARACTER(LEN=strlen_long) :: txt_err
    INTEGER :: ity=0
    INTEGER :: ikp=0, ile=0
    INTEGER :: idt

    ! Here come the SMCL subroutines which do the computations:

    ! plNOx_0(ity) refers to the start value of the plume NOx of
    ! type ity (Leapfrog):
    all_ty_plNOx_0: DO ity=1, max_ty
      IF (.NOT. l_pl_ty(ity)) CYCLE
      idt = idx_plNOx(ity)
      plNOx_0(ity)%PTR(1:kproma,:)=            &
         qxtm1(_RI_X_ZN_(1:kproma,:,idt))               &
        +qxtte(_RI_X_ZN_(1:kproma,:,idt))*time_step_len
    ! time_step_len is 2*delta_time
    ENDDO all_ty_plNOx_0

    idt = idx_NO
    NO_0(1:kproma,:)=                         &
       qxtm1(_RI_X_ZN_(1:kproma,:,idt))                &
      +qxtte(_RI_X_ZN_(1:kproma,:,idt))*time_step_len

    idt = idx_NO2
    NO2_0(1:kproma,:)=                        &
       qxtm1(_RI_X_ZN_(1:kproma,:,idt))                &
      +qxtte(_RI_X_ZN_(1:kproma,:,idt))*time_step_len

    idt = idx_O3
    O3_0(1:kproma,:)=                         &
       qxtm1(_RI_X_ZN_(1:kproma,:,idt))                & 
      +qxtte(_RI_X_ZN_(1:kproma,:,idt))*time_step_len

    molecules_m3(1:kproma,:)=       &
       ptr_grmass(_RI_XYZ__(1:kproma,jrow,:)) &
      /ptr_grvol(_RI_XYZ__( 1:kproma,jrow,:)) & ! kg m-3,
      /M_air*1000._dp*N_A             ! g-1 mol, g kg-1, molecules mol-1

    CALL plg_par( & !status,txt_err & ! A core layer error message
                kproma  &             ! string (not yet implemented)
                ,nlev    &
                ,jrow    &
                )
    !IF (status /= 0) &
    ! CALL error_bi('Error in plg_par'//txt_err//'.', substr)

    all_ty_update_tendencies: DO ity=1,max_ty
      IF (.NOT. l_pl_ty(ity)) CYCLE
      idt = idx_plNOx(ity)
      qxtte(_RI_X_ZN_(1:kproma,:,idt))=     &
         qxtte(_RI_X_ZN_(1:kproma,:,idt))   &
        +in_plNOx(ity)%PTR(_RI_XYZ__(1:kproma,jrow,:))

      IF (.NOT. l_decouple) THEN

        ! The increments of the NOx tendency are assigned to the
        ! NO and NO2 tendency according to the setting of as_NOx_NO
        ! in the CPL namelist:

        as_NO_NO2: IF (.NOT. l_grb_NO) THEN
           idt = idx_NO
           qxtte(_RI_X_ZN_(1:kproma,:,idt)) =              &
                qxtte(_RI_X_ZN_(1:kproma,:,idt))           &
                +as_NO(ity)                       &
                *in_NOx(ity)%PTR(_RI_XYZ__(1:kproma,jrow,:))

           idt = idx_NO2
           qxtte(_RI_X_ZN_(1:kproma,:,idt))=               &
                qxtte(_RI_X_ZN_(1:kproma,:,idt))           &
                +as_NO2(ity)                      &
                *in_NOx(ity)%PTR(_RI_XYZ__(1:kproma,jrow,:))

        ELSE ! as_NOx_NO is set to 'gridbox'

          DO ile=1,nlev
            DO ikp=1,kproma

              IF (NO_0(ikp,ile)+NO2_0(ikp,ile) > al_0) THEN

                grb_NO_NOx(ikp,ile)= &
                   NO_0(ikp,ile)     &
                  /(NO_0(ikp,ile)    &
                  +NO2_0(ikp,ile))

              ELSE ! With zero NOx in the gridbox, NOx increments are
!                    assigned to NO:

                grb_NO_NOx(ikp,ile)=1._dp

              ENDIF

            ENDDO
          ENDDO

          idt = idx_NO
          qxtte(_RI_X_ZN_(1:kproma,:,idt))=               &
               qxtte(_RI_X_ZN_(1:kproma,:,idt))           &
               +grb_NO_NOx(1:kproma,:)           &
               *in_NOx(ity)%PTR(_RI_XYZ__(1:kproma,jrow,:))

          idt = idx_NO2
          qxtte(_RI_X_ZN_(1:kproma,:,idt))=            &
             qxtte(_RI_X_ZN_(1:kproma,:,idt))          &
            +(1._dp                           &
             -grb_NO_NOx(1:kproma,:))         &
            *in_NOx(ity)%PTR(_RI_XYZ__(1:kproma,jrow,:))

        ENDIF as_NO_NO2

        idt = idx_HNO3
        qxtte(_RI_X_ZN_(1:kproma,:,idt))=                 &
             qxtte(_RI_X_ZN_(1:kproma,:,idt))             &
             +in_HNO3(ity)%PTR(_RI_XYZ__(1:kproma,jrow,:))

        idt = idx_O3
        qxtte(_RI_X_ZN_(1:kproma,:,idt))=                   &
             qxtte(_RI_X_ZN_(1:kproma,:,idt))               &
             +in_O3_titr(ity)%PTR(_RI_XYZ__(1:kproma,jrow,:))

        idt = idx_O3
        qxtte(_RI_X_ZN_(1:kproma,:,idt))=                    &
             qxtte(_RI_X_ZN_(1:kproma,:,idt))                &
             +in_O3_keff(ity)%PTR(_RI_XYZ__(1:kproma,jrow,:))

      ENDIF

    ENDDO all_ty_update_tendencies

  END SUBROUTINE plumegas_physc
! ---------------------------------------------------------------------

! Regional loop end

! Time loop end

! ---------------------------------------------------------------------
! Finalizing phase:
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plumegas_free_memory

    ! plumegas_free_memory
    !    This subroutine deallocates non-channel and non-tracer memory.


    IMPLICIT NONE

    INTEGER :: ity=0

    ! Deallocate the local fields with the start values of the tracers
    ! used in the parametrisation:

    DEALLOCATE(grb_NO_NOx)

    all_ty: DO ity=1,max_ty
      IF (.NOT. l_pl_ty(ity)) CYCLE

      DEALLOCATE(plNOx_0(ity)%PTR)

    ENDDO all_ty

    DEALLOCATE(NO_0  &
              ,NO2_0 &
              ,O3_0  &
              )

    DEALLOCATE(molecules_m3)

  END SUBROUTINE plumegas_free_memory
! ---------------------------------------------------------------------

! #####################################################################
! Private subroutines
! #####################################################################

! ---------------------------------------------------------------------
  SUBROUTINE plumegas_read_nml_cpl(status,iou)

    ! plumegas_read_nml_cpl
    !    This private subroutine reads the coupling namelist.
    !
    ! Bug reports to
    !    ECHAM Group, DLR, Oberpfaffenhofen, Germany

    ! BMIL:
    USE messy_main_tools,ONLY: read_nml_open  &
                              ,read_nml_check &
                              ,read_nml_close

    IMPLICIT NONE


    ! Input / output varaibles:

    INTEGER,INTENT(OUT) :: status ! There is a RETURN

    INTEGER,INTENT(IN)  :: iou ! Input / output unit


    NAMELIST /CPL/ l_decouple &
                 , as_NOx_NO


    CHARACTER(LEN=*),PARAMETER :: substr='plumegas_read_nml_cpl'
    LOGICAL                    :: lex   ! File exists
    INTEGER                    :: fstat ! File status

    status=1

    CALL read_nml_open(lex    &
                      ,substr &
                      ,iou    &
                      ,'CPL'  &
                      ,modstr &
                      )
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou          &
        ,NML=CPL      &
        ,IOSTAT=fstat &
        )

    CALL read_nml_check(fstat  &
                       ,substr &
                       ,iou    &
                       ,'CPL'  &
                       ,modstr &
                       )
    IF (fstat /= 0) RETURN  ! Error while reading namelist

    ! Diagnose namelist and set global switches:

    IF (l_decouple) THEN
      WRITE(*,*) 'Tendencies of external tracers are not modified.'
    END IF

    IF (.NOT. l_decouple) THEN

      IF (as_NOx_NO == '') as_NOx_NO='gridbox'

      SELECT CASE (as_NOx_NO)
      CASE ('gridbox')
        WRITE(*,*) 'as_NOx_NO ='//as_NOx_NO//' , i.e. NOx '//      &
                   'NOx increments assigned to NO (and NO2) '//    &
                   'according to the gridbox NO/NOx (NO2/NOx) '//  &
                   'molar ratio, only to NO if NOx is zero in ' // &
                   'the gridbox.'

      CASE ('exhaust')
        WRITE(*,*) 'NOx increments assigned to NO according to '//  &
                   'the emission ratio emi_ra of the plume type '// &
                   'in question.'

      CASE ('only_NO')
        WRITE(*,*) 'as_NOx_NO ='//as_NOx_NO//' , i.e. NOx '// &
                   'increments assigned to NO only.'

      CASE DEFAULT
        WRITE(*,*) 'as_NOx_NO ='//as_NOx_NO//' . Error: '//   &
                   'as_NOx_NO not defined properly in the '// &
                   'CPL namelist.'
        RETURN

      END SELECT

    ENDIF

    CALL read_nml_close(substr &
                       ,iou    &
                       ,modstr &
                       )

    status=0 ! No error

  END SUBROUTINE plumegas_read_nml_cpl
! ---------------------------------------------------------------------

! **********************************************************************
END MODULE messy_plumegas_si
! **********************************************************************
