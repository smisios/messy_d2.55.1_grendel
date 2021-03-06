#include "messy_main_ppd_bi.inc"

!*****************************************************************************
!                Time-stamp: <2021-02-19 14:06:12 b302010>
!*****************************************************************************

! Authors: Astrid Kerkweg,  MPICH, Mainz,  2003
! Rolf Sander, 2003: strict separation between core and e5 files
! Patrick Joeckel, MPICH Mar 2004: strict separation of AERO from MECCA

!*****************************************************************************

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************

MODULE messy_mecca_aero_si

  ! ECHAM5/MESSy
  USE messy_main_grid_def_mem_bi,       ONLY: jrow
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
    info_bi, warning_bi, error_bi
  USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                     mtend_get_start_l,      &
                                     mtend_add_l,            &
                                     mtend_register,         &
                                     mtend_id_tracer,        &
                                     mtend_id_t, mtend_id_q
#endif
  ! MESSy
  USE messy_main_constants_mem, ONLY: I4, STRLEN_MEDIUM, HLINE2
  USE messy_main_tools,         ONLY: mass_density, layerthickness, &
                                      PTR_2D_ARRAY, PTR_3D_ARRAY, str
  ! SUBMODEL
  USE messy_mecca_aero,          ONLY: submodstr, &
                                      mecca_aero_calc_k_ex
  USE messy_mecca_kpp,            ONLY: dp, ind_Clm_a, ind_Brm_a  &
                                      , ind_H2O_a, ind_Im_a  &
                                      , ind_IO3m_a, ind_HCO3m_a, ind_Hp_a  &
                                      , xnom7sulf, APN                     &
                                      , initialize_indexarrays

  USE messy_mecca_mem_si        ! jp, jk, idt_*
  USE messy_mecca,              ONLY: modstr

  IMPLICIT NONE
  SAVE
  PRIVATE

  INTRINSIC ABS, NULL, TINY, TRIM

  INTEGER, DIMENSION(APN) :: zidt_num = 0
#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  REAL(dp), DIMENSION(:,:),     POINTER :: heightidx       => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: wetradius       => NULL()
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: dryradius       => NULL()
  REAL(dp), DIMENSION(:),       POINTER :: sigma_aerosol   => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: seasalt_emis
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: seasalt_ddep

  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: Nss_emis
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: Nss_ddep
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: Clm_ddep
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: Brm_ddep
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: HCO3m_ddep
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: Im_ddep
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: IO3m_ddep
  TYPE(PTR_2D_ARRAY), DIMENSION(APN)    :: Hp_ddep

  ! CPL-NAMELIST VARIABLES
  ! -- cpl to aerosol chemistry  i.e sulfate, halogen chem.
  LOGICAL, DIMENSION(APN) :: lcpl_aerophase   = .FALSE.
  !
  CHARACTER(LEN=STRLEN_MEDIUM) :: aerosol_module     = ''
  !
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)       :: number_main  = ''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)       :: number_sub   = ''
  CHARACTER(LEN=STRLEN_MEDIUM)  :: aerosol_channel   = ''
  CHARACTER(LEN=STRLEN_MEDIUM)  :: aerosol_wetradius = ''
  CHARACTER(LEN=STRLEN_MEDIUM)  :: aerosol_dryradius = ''
  CHARACTER(LEN=STRLEN_MEDIUM)  :: sigma             = ''
  INTEGER(I4) , DIMENSION(APN)  :: modenumber        = 0
  CHARACTER(LEN=STRLEN_MEDIUM)  :: height_channel    =''
  CHARACTER(LEN=STRLEN_MEDIUM)  :: height_index      =''
  CHARACTER(LEN=STRLEN_MEDIUM)  :: emis_channel      =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)       :: emis_ss  =''
  CHARACTER(LEN=STRLEN_MEDIUM)  :: ddep_channel      =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)       :: ddep_ss  =''
  LOGICAL                       :: l_tendency        = .FALSE.
  LOGICAL                       :: l_calc_emis       = .FALSE.
  LOGICAL                       :: l_netflux         = .FALSE.

  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)  :: str_Nss_emis   =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)  :: str_Nss_ddep   =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)  :: str_Clm_ddep   =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)  :: str_Brm_ddep   =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)  :: str_HCO3m_ddep =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)  :: str_Hp_ddep    =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)  :: str_Im_ddep    =''
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(APN)  :: str_IO3m_ddep  =''

  NAMELIST /CPL_AERO/                         &
    lcpl_aerophase,    aerosol_module,       &
    number_main,       number_sub,           &
    aerosol_channel,   aerosol_wetradius,    &
    aerosol_dryradius, sigma,                &
    modenumber,                              &
    height_channel,    height_index,         &
    l_calc_emis,       l_netflux,            &
    emis_channel,      emis_ss,              &
    ddep_channel,      ddep_ss,              &
    l_tendency,                              &
    str_Nss_emis,      str_Nss_ddep,         &
    str_Clm_ddep,      str_Brm_ddep,         &
    str_HCO3m_ddep,    str_Im_ddep,          &
    str_IO3m_ddep,     str_Hp_ddep


  ! CHANNEL OBJECTS
  TYPE(PTR_3D_ARRAY), DIMENSION(APN) :: lwc_3d
  TYPE(PTR_3D_ARRAY), DIMENSION(APN) :: pH

  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC :: MODEL
  INTEGER                     , PUBLIC :: METHD
  INTEGER, DIMENSION(APN)     , PUBLIC :: m

  ! LOCAL variables
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: lwc
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xaer
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: radius
  INTEGER, DIMENSION(APN) :: idt_Hp_aq    = 0
  INTEGER, DIMENSION(APN) :: idt_Clm_aq   = 0
  INTEGER, DIMENSION(APN) :: idt_Brm_aq   = 0
  INTEGER, DIMENSION(APN) :: idt_Im_aq    = 0
  INTEGER, DIMENSION(APN) :: idt_IO3m_aq  = 0
  INTEGER, DIMENSION(APN) :: idt_HCO3m_aq = 0
  INTEGER, DIMENSION(APN) :: idt_Nap_aq   = 0
  INTEGER :: idt_Brsalt  = 0
  INTEGER :: idt_Brorg   = 0
  INTEGER :: idt_BrSScap = 0

  ! PUBLIC SUBROUTINES
  PUBLIC :: mecca_aero_initialize
  PUBLIC :: mecca_aero_init_memory
  PUBLIC :: mecca_aero_init_coupling
  PUBLIC :: mecca_aero_vdiff
  !
  ! CALLED FROM mecca_physc:
  PUBLIC :: mecca_aero_update_physc
  !         CONTAINS
  !           modal_lwc
  !           modal_radii
  !           constant_lwc
  !           constant_radii
  !           mecca_aero_xaer
  PUBLIC :: mecca_aero_mr2c
  PUBLIC :: mecca_aero_save_pH
  PUBLIC :: mecca_aero_dealloc
  PUBLIC :: mecca_aero_diag_si

CONTAINS

  !*****************************************************************************

  SUBROUTINE mecca_aero_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_bcast, p_io
    ! MESSy
    USE messy_main_tracer,        ONLY: MODAL, BIN
    USE messy_main_tools,         ONLY: find_next_free_unit

    IMPLICIT NONE

    ! local
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_aero_initialize'
    INTEGER :: iou    ! I/O unit
    INTEGER :: status ! error status
    INTEGER :: i

    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL mecca_aero_read_nml_cpl(status, iou)
      IF (status /= 0) CALL error_bi('error in mecca_aero_read_nml_cpl',substr)
      ! COUNT MODES
      DO i=1, APN
        IF (TRIM(number_main(i)) == '') &
          CALL error_bi('main name of aerosol number tracer empty',substr)
        IF (TRIM(number_sub(i)) == '') &
          CALL warning_bi('subname of aerosol number tracer empty',substr)
        IF (modenumber(i) == 0 ) &
          CALL error_bi('mode number missing',substr)
        IF (l_calc_emis) THEN
          IF (TRIM(emis_ss(i)) == '') &
            CALL error_bi('seasalt emission missing',substr)
          IF (TRIM(ddep_ss(i)) == '' .AND. l_netflux) &
            CALL error_bi('no deposition flux available',substr)
        ENDIF
      ENDDO
    ENDIF

    ! BROADCAST CPL RESULTS
    DO i = 1,APN
      CALL p_bcast(lcpl_aerophase(i),  p_io)
      CALL p_bcast(number_main(i),     p_io)
      CALL p_bcast(number_sub(i),      p_io)
      CALL p_bcast(modenumber(i),      p_io)
      CALL p_bcast(emis_ss(i),         p_io)
      CALL p_bcast(str_Nss_emis(i),    p_io)
      CALL p_bcast(ddep_ss(i),         p_io)
      CALL p_bcast(str_Nss_ddep(i),    p_io)
      CALL p_bcast(str_Clm_ddep(i),    p_io)
      CALL p_bcast(str_Brm_ddep(i),    p_io)
      CALL p_bcast(str_HCO3m_ddep(i),  p_io)
      CALL p_bcast(str_Hp_ddep(i),     p_io)
      CALL p_bcast(str_Im_ddep(i),     p_io)
      CALL p_bcast(str_IO3m_ddep(i),   p_io)
    ENDDO
    CALL p_bcast(aerosol_module,     p_io)
    CALL p_bcast(aerosol_channel,    p_io)
    CALL p_bcast(aerosol_wetradius,  p_io)
    CALL p_bcast(aerosol_dryradius,  p_io)
    CALL p_bcast(sigma,              p_io)
    CALL p_bcast(height_channel,     p_io)
    CALL p_bcast(height_index,       p_io)
    CALL p_bcast(emis_channel,       p_io)
    CALL p_bcast(ddep_channel,       p_io)
    CALL p_bcast(l_tendency,         p_io)
    CALL p_bcast(l_calc_emis,        p_io)
    CALL p_bcast(l_netflux,          p_io)

    !SET PARAMETER FOR NEW_TRACER DEFINITION
    MODEL    = aerosol_channel
    m(1:APN) = modenumber(1:APN)

    SELECT CASE (TRIM(aerosol_module))
    CASE ('modal')
      METHD   = MODAL
    CASE ('bin')
      METHD   = BIN
    CASE DEFAULT
      CALL warning_bi(&
        'no know or constant method choosen, setting METHD to BIN',substr)
      METHD   = BIN
    END SELECT

    ! define gas-aq physicochemical constants:
    IF (p_parallel_io) THEN
      PRINT *, HLINE2
      PRINT *, "         Henry's law coefficients "// &
        "and accommodation coefficients"
      PRINT *, HLINE2
      PRINT *, 'species           Henry_T0 Henry_Tdep'// &
        '   alpha_T0 alpha_Tdep         M'
      PRINT *, '                   [M/atm]        [K]'// &
        '        [1]        [K]  [kg/mol]'
      PRINT *, HLINE2
    ENDIF

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(submodstr)
#endif

  END SUBROUTINE mecca_aero_initialize

  ! --------------------------------------------------------------------------

  SUBROUTINE mecca_aero_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                           new_attribute
    USE messy_main_tools,            ONLY: str
    USE messy_main_constants_mem,    ONLY: STRLEN_LONG

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr='mecca_aero_init_memory'
    INTEGER                      :: status
    INTEGER                      :: i
    CHARACTER(LEN=STRLEN_LONG+5) :: name

    CALL start_message_bi(submodstr, ' INITIALIZE MECCA-AERO MEMORY', substr)

    ! moved here, as tracer meta is not available earlier
    CALL mecca_aero_init_gasaq(l_print=p_parallel_io)

#ifdef MESSYTENDENCY
    IF (l_tendency) CALL mtend_register(my_handle, mtend_id_tracer)
#endif

    CALL new_channel(status, submodstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    DO i = 1,APN
      name = 'lwc_a'//str(i,'(I2.2)')
      CALL new_channel_object(status, submodstr &
        , TRIM(name), p3=lwc_3d(i)%ptr)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, submodstr, TRIM(name) &
        , 'long_name', c='aerosol liquid water content')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, submodstr, TRIM(name) &
        , 'units', c='m3(aq)/m3(air)')
      CALL channel_halt(substr, status)
      lwc_3d(i)%ptr(:,:,:)    = -999.999_dp ! dummy default
    ENDDO

    DO i = 1,APN
      name = 'pH_a'//str(i,'(I2.2)')
      CALL new_channel_object(status, submodstr &
        , TRIM(name), p3=pH(i)%ptr)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, submodstr, TRIM(name) &
        , 'long_name', c='aerosol pH ')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, submodstr, TRIM(name) &
        , 'units', c=' ')
      CALL channel_halt(substr, status)
      pH(i)%ptr(:,:,:)     = -999.999_dp !dummy default
    ENDDO

    CALL end_message_bi(submodstr, ' INITIALIZE MECCA-AERO MEMORY', substr)

  END SUBROUTINE mecca_aero_init_memory

  !---------------------------------------------------------------------------

  SUBROUTINE mecca_aero_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,          ONLY: p_bcast
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: get_tracer, OFF
    USE messy_main_channel,       ONLY: get_channel_object

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER::substr='mecca_aero_init_coupling'
    INTEGER :: status, status2, status3
    INTEGER :: i
    CHARACTER(LEN=STRLEN_MEDIUM+3) :: z_aerosol_channel = ''

    INTRINSIC TRIM

    CALL start_message_bi(submodstr, 'TEST COUPLING', substr)

    ! COUPLING for diagnostic reasons
    CALL get_tracer(status, GPTRSTR, 'Brsalt',idx=idt_Brsalt)
    CALL get_tracer(status2, GPTRSTR,'Brorg',idx=idt_Brorg)
    CALL get_tracer(status3, GPTRSTR,'BrSScap',idx=idt_BrSScap)
    IF (status==0 .AND. status2==0 .AND. status3==0) THEN
       CALL info_bi('got idt of diagnostic tracers Brsalt, Brorg and BrSScap')
    ELSE
       CALL warning_bi(&
            'no diagnistic of captured salt/org Br possible',substr)
    ENDIF

    ! INITIALIZE index arrays
    CALL initialize_indexarrays
    CALL initialize_indexarrays_trac

    DO i=1,APN
      IF ( lcpl_aerophase(i)) xnom7sulf = 0
    ENDDO

    SELECT CASE (TRIM(aerosol_module))
    CASE ('modal','bin')
      DO i=1,APN
        CALL get_tracer(status, GPTRSTR, number_main(i) &
          , subname=number_sub(i), idx=zidt_num(i))
        IF (status==0) THEN
          CALL info_bi('fetching aerosol particle number idt')
        ELSE
          CALL warning_bi(&
            'no aerosol module particle number available',substr)
          CALL tracer_halt(substr, status)
        ENDIF
      ENDDO

      CALL info_bi('getting aerosol properties ...')

      z_aerosol_channel = TRIM(aerosol_channel)//'_'//GPTRSTR

      CALL get_channel_object(status &
        ,TRIM(z_aerosol_channel), TRIM(aerosol_wetradius) , p4=wetradius)
      IF (status /= 0) &
        CALL error_bi(&
        'channel object for aerosol ambient radius not found',substr)

      CALL get_channel_object(status &
        ,TRIM(z_aerosol_channel), TRIM(aerosol_dryradius) , p4=dryradius)
      IF (status /= 0) &
        CALL error_bi('channel object for aerosol dry radius not found',substr)

      CALL get_channel_object(status &
        ,TRIM(z_aerosol_channel), TRIM(sigma) , p1=sigma_aerosol)
      IF (status /= 0) &
        CALL error_bi( &
        'channel object for sigma of aerosol modes not found',substr)

    CASE DEFAULT

      CALL warning_bi('constant or no aerosol, nothing to be done',substr)

    ENDSELECT

    ! get maximum height for chemistry calculations
    CALL get_channel_object(status &
      ,TRIM(height_channel), TRIM(height_index) , p2=heightidx)
    IF (status /= 0) &
      CALL error_bi (&
      'channel object for height determination not available',substr)

    IF (l_calc_emis) THEN

      ! get seasalt emission
      DO i=1,APN
        IF (TRIM(emis_ss(i))/='') THEN
          CALL info_bi('getting sea salt emissions phase ')
          CALL get_channel_object(status &
            ,TRIM(emis_channel), TRIM(emis_ss(i)) &
            , p2=seasalt_emis(i)%ptr)
          IF (status /= 0) &
            CALL error_bi (&
            'channel object for seasalt emission not available',substr)
        ENDIF
      ENDDO

      IF (l_netflux) THEN
        DO i=1,APN
           ! Get number emission flux
           IF (TRIM(str_Nss_emis(i))/='') THEN
              CALL get_channel_object(status,TRIM(emis_channel) &
                   ,TRIM(str_Nss_emis(i)) , p2=Nss_emis(i)%ptr)
              IF (status /= 0) CALL error_bi(&
                  'channel object for number seasalt cs emission not available'&
                  , substr)
           ELSE
              IF (idt_Clm_aq(i)/=0) CALL error_bi(&
                'the number emission flux is necessary, if l_netflux=T' &
                , substr)
           ENDIF
           ! get seasalt dry deposition
          IF (TRIM(ddep_ss(i))/='') THEN
            CALL info_bi('fetching dry deposition flux phase ')
            CALL get_channel_object(status &
              ,TRIM(ddep_channel), TRIM(ddep_ss(i)) &
              , p2=seasalt_ddep(i)%ptr)
            IF (status /= 0) &
              CALL error_bi (&
              'channel object for seasalt dry deposition not available',substr)
          ENDIF


          !!! Get number deposition flux
           IF (TRIM(str_Nss_ddep(i))/='') THEN
              CALL get_channel_object(status,TRIM(ddep_channel) &
                   ,TRIM(str_Nss_ddep(i)) , p2=Nss_ddep(i)%ptr)
              IF (status /= 0) CALL error_bi(&
                  'channel object for number seasalt deposition not available'&
                  , substr)
           ELSE
              IF (idt_Clm_aq(i)/=0) CALL error_bi(&
                'the number deposition flux is necessary, if l_netflux=T' &
                , substr)
           ENDIF

           !!! Get chloride deposition flux
           IF (TRIM(str_Clm_ddep(i))/='') THEN
              CALL get_channel_object(status,TRIM(ddep_channel) &
                   ,TRIM(str_Clm_ddep(i)) , p2=Clm_ddep(i)%ptr)
              IF (status /= 0 .and. (idt_Clm_aq(i)/=0)) CALL error_bi(&
                  'channel object for Clm deposition not available'&
                  ,substr)
           ELSE
              IF (idt_Clm_aq(i)/=0) CALL error_bi(&
                'Clm deposition flux must be given for l_netflux=T'&
                , substr)
           ENDIF
           !!! Get bromide deposition flux
           IF (TRIM(str_Brm_ddep(i))/='') THEN
              CALL get_channel_object(status,TRIM(ddep_channel) &
                   ,TRIM(str_Brm_ddep(i)) , p2=Brm_ddep(i)%ptr)
              IF (status /= 0 .and. (idt_Brm_aq(i)/=0) ) CALL error_bi(&
                  'channel object for Brm deposition not available'&
                  , substr)
           ELSE
              IF (idt_Brm_aq(i)/=0) CALL error_bi(&
                'Brm deposition flux must be given for l_netflux=T' &
                , substr)
           ENDIF
           !!! Get HCO3m deposition flux
           IF (TRIM(str_HCO3m_ddep(i))/='') THEN
              CALL get_channel_object(status,TRIM(ddep_channel) &
                   ,TRIM(str_HCO3m_ddep(i)) , p2=HCO3m_ddep(i)%ptr)
              IF (status /= 0 .and. (idt_HCO3m_aq(i)/=0)) CALL error_bi(&
                  'channel object for HCO3m deposition not available'&
                  , substr)
           ELSE
              IF (idt_HCO3m_aq(i)/=0) CALL error_bi(&
                'HCO3m deposition flux must be given for l_netflux=T' &
                , substr)
           ENDIF
           IF (TRIM(str_Hp_ddep(i))/='') THEN
              CALL get_channel_object(status,TRIM(ddep_channel) &
                   ,TRIM(str_Hp_ddep(i)) , p2=Hp_ddep(i)%ptr)
              IF (status /= 0 .and. (idt_Hp_aq(i)/=0)) CALL error_bi(&
                  'channel object for HCO3m deposition not available'&
                  , substr)
           ELSE
              IF (idt_Hp_aq(i)/=0) CALL error_bi(&
                'Hp deposition flux must be given for l_netflux=T' &
                , substr)
           ENDIF
           !!! Get iodide deposition flux
           IF (TRIM(str_Im_ddep(i))/='') THEN
              CALL get_channel_object(status,TRIM(ddep_channel) &
                   ,TRIM(str_Im_ddep(i)) , p2=Im_ddep(i)%ptr)
              IF (status /= 0 .and. (idt_Im_aq(i)/=0)) CALL error_bi(&
                  'channel object for Im deposition not available'&
                  , substr)
           ELSE
              IF (idt_Im_aq(i)/=0) CALL error_bi(&
                'Im deposition flux must be given for l_netflux=T' &
                , substr)
           ENDIF
           !!! Get iodate deposition flux
           IF (TRIM(str_IO3m_ddep(i))/='') THEN
              CALL get_channel_object(status,TRIM(ddep_channel) &
                   ,TRIM(str_IO3m_ddep(i)) , p2=IO3m_ddep(i)%ptr)
              IF (status /= 0 .and. (idt_IO3m_aq(i)/=0)) CALL error_bi(&
                  'channel object for IO3m deposition not available'&
                  , substr)
           ELSE
              IF (idt_IO3m_aq(i)/=0) CALL error_bi(&
                'IO3m deposition flux must be given for l_netflux=T' &
                , substr)
           ENDIF
        ENDDO
      ENDIF !l_netflux

    ENDIF ! l_calc_emis

    CALL end_message_bi(submodstr, 'TEST COUPLING', substr)

  END SUBROUTINE mecca_aero_init_coupling

  ! --------------------------------------------------------------------------

  SUBROUTINE mecca_aero_vdiff(flag)

    ! ECHAM5/MESSy
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1
    USE messy_main_data_bi,         ONLY: tm1_3d, tte_3d   &
                                        , qte_3d, qm1_3d
#endif
    USE messy_main_grid_def_mem_bi, ONLY: nlev, kproma, nproma
    USE messy_main_grid_def_bi,     ONLY: deltaz
    USE messy_main_data_bi,         ONLY: press_3d   &
#ifdef ECHAM5
                                        , pxtems &
#endif 
                                        , pressi_3d
    ! MESSy
    USE messy_main_timer,         ONLY: time_step_len, lstart
    USE messy_main_constants_mem, ONLY: M_air, g, N_A

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='mecca_aero_vdiff'
    REAL(dp), POINTER :: zxtems(:,:)
#ifdef MESSYTENDENCY
    REAL(dp) :: tend_tmp(nproma,nlev)
#endif
    REAL(dp), PARAMETER  :: CONST= g*1.e-3*M_air
    REAL(dp), PARAMETER :: M_NaCl = 58.443_dp ! g/mol
    REAL(dp) :: temp(nproma), zsphum(nproma)
    REAL(dp) :: rho_air(nproma)
    REAL(dp) :: zdz(nproma)
    REAL(dp) :: zdp(nproma)
    REAL(dp) :: Cl_tendency(1:nproma)
    REAL(dp) :: N_new(1:kproma)
    REAL(dp) :: N_old(1:kproma)
    REAL(dp) :: N_sum(1:kproma)
    REAL(dp) :: N_ratio(1:kproma)

    REAL(dp) :: net_Clm_ddep(1:kproma)
    REAL(dp) :: net_Brm_ddep(1:kproma)
    REAL(dp) :: net_Im_ddep(1:kproma)
    REAL(dp) :: net_IO3m_ddep(1:kproma)
    REAL(dp) :: net_HCO3m_ddep(1:kproma)
    REAL(dp) :: net_Hp_ddep(1:kproma)

    REAL(dp) :: Nap_depsum(1:kproma)
    REAL(dp) :: CONV(1:kproma)

    INTEGER  :: i

#if defined(ECHAM5)
    ! INIT
    zxtems => pxtems(:,1,:,jrow)
#endif
    IF (lstart) RETURN
    IF (.NOT. l_calc_emis) RETURN

    SELECT CASE (flag)
    CASE (1)
    doapn: DO i=1,APN
      IF(idt_Clm_aq(i)/= 0) THEN
        IF (ti_gp(idt_Clm_aq(i))%tp%ident%unit /= 'mol/mol') &
          CALL error_bi(&
          'emission calculation only for tracer unit mol/mol possible !',substr)
      ENDIF

#ifdef COSMO
      IF (.NOT. l_tendency) THEN
         CALL info_bi('COSMO works only for l_tendency=T','')
         CALL info_bi('!! Setting l_tendency to TRUE !!',substr)
         l_tendency = .TRUE.
      ENDIF
#endif

      iftendency: IF (.NOT. l_tendency) THEN
        ! Emissions  Coarse Mode
          Cl_tendency(1:kproma) =  &
            seasalt_emis(i)%ptr(1:kproma,jrow) * M_air/M_NaCl

#if defined(ECHAM5)
        IF (idt_Clm_aq(i) /= 0 )   THEN
          zxtems(1:kproma,idt_Clm_aq(i)) = zxtems(1:kproma,idt_Clm_aq(i)) &
            +  Cl_tendency(1:kproma)
          IF (idt_Nap_aq(i) /=0 ) &
            zxtems(1:kproma,idt_Nap_aq(i)) = zxtems(1:kproma,idt_Nap_aq(i)) &
            +  Cl_tendency(1:kproma)
        ENDIF

        IF (idt_Brm_aq(i) /= 0 )  THEN
          zxtems(1:kproma,idt_Brm_aq(i))   =  zxtems(1:kproma,idt_Brm_aq(i))+ &
            Cl_tendency(1:kproma) *  1.5E-3_dp
          IF (idt_Nap_aq(i) /=0 ) &
            zxtems(1:kproma,idt_Nap_aq(i))   =  zxtems(1:kproma,idt_Nap_aq(i))+ &
            Cl_tendency(1:kproma) *  1.5E-3_dp
        ENDIF

        IF (idt_Im_aq(i) /= 0 )  THEN
          zxtems(1:kproma,idt_Im_aq(i))    = zxtems(1:kproma,idt_Im_aq(i)) + &
            Cl_tendency(1:kproma) *  7.4E-8_dp/0.545_dp
          IF (idt_Nap_aq(i) /=0 ) &
            zxtems(1:kproma,idt_Nap_aq(i))    = zxtems(1:kproma,idt_Nap_aq(i)) + &
            Cl_tendency(1:kproma) *  7.4E-8_dp/0.545_dp
        ENDIF
        IF (idt_IO3m_aq(i) /= 0 )  THEN
          zxtems(1:kproma,idt_IO3m_aq(i))  = zxtems(1:kproma,idt_IO3m_aq(i)) + &
            Cl_tendency(1:kproma) *  2.64E-7_dp/0.545_dp
          IF (idt_Nap_aq(i) /=0 ) &
            zxtems(1:kproma,idt_Nap_aq(i))  = zxtems(1:kproma,idt_Nap_aq(i)) + &
            Cl_tendency(1:kproma) *  2.64E-7_dp/0.545_dp
        ENDIF
        IF (idt_HCO3m_aq(i) /= 0 ) THEN
          zxtems(1:kproma,idt_HCO3m_aq(i)) =  zxtems(1:kproma,idt_HCO3m_aq(i)) + &
            Cl_tendency(1:kproma) *  4.2E-3_dp
          IF (idt_Nap_aq(i) /=0 ) &
            zxtems(1:kproma,idt_Nap_aq(i)) =  zxtems(1:kproma,idt_Nap_aq(i)) + &
            Cl_tendency(1:kproma) *  4.2E-3_dp
        ENDIF
        IF (idt_Hp_aq(i) /= 0 ) THEN
          zxtems(1:kproma,idt_Hp_aq(i)) =  zxtems(1:kproma,idt_Hp_aq(i)) + &
            Cl_tendency(1:kproma) *  4.2E-3_dp
          IF (idt_Nap_aq(i) /=0 ) &
            zxtems(1:kproma,idt_Nap_aq(i)) =  zxtems(1:kproma,idt_Nap_aq(i)) + &
            Cl_tendency(1:kproma) *  4.2E-3_dp
        ENDIF
#endif
      ELSE !iftendency
#ifndef MESSYTENDENCY
        temp(1:kproma) = tm1_3d(_RI_XYZ__(1:kproma,jrow,nlev)) + &
          tte_3d(_RI_XYZ__(1:kproma,jrow,nlev)) * time_step_len
        zsphum(1:kproma) = qm1_3d(_RI_XYZ__(1:kproma,jrow,nlev)) + &
          qte_3d(_RI_XYZ__(1:kproma,jrow,nlev)) * time_step_len
#else
        CALL mtend_get_start_l(mtend_id_t, v0=tend_tmp)
        temp(1:kproma) = tend_tmp(1:kproma,nlev)
        CALL mtend_get_start_l(mtend_id_q, v0=tend_tmp)
        zsphum(1:kproma) = tend_tmp(1:kproma,nlev)
#endif
        ! calculate air density in mol/m3
        rho_air(1:kproma) = &
          mass_density(press_3d(_RI_XYZ__(1:kproma,jrow,nlev))  &
          ,temp(1:kproma)     &
          ,zsphum(1:kproma)    ) * 1.e3_dp / M_air

        zdz(1:kproma) = deltaz(_RI_XYZ__(1:kproma,jrow,nlev))

        ! Emissions

#ifndef MESSYTENDENCY
          Cl_tendency(1:kproma) =  &
            seasalt_emis(i)%ptr(1:kproma,jrow) * 1.e3 / M_NaCl /  &
            zdz(1:kproma) /rho_air(1:kproma)
        IF (idt_Clm_aq(i) /= 0) THEN
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Clm_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Clm_aq(i))) +    &
            Cl_tendency(1:kproma)
          IF (idt_Nap_aq(i) /= 0) THEN
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) = &
              pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) +    &
              Cl_tendency(1:kproma)
          ENDIF
        ENDIF

        IF (idt_Brm_aq(i) /= 0) THEN
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Brm_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Brm_aq(i))) +     &
            Cl_tendency(1:kproma) *  1.5E-3
          IF (idt_Nap_aq(i) /= 0)   &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) +     &
            Cl_tendency(1:kproma) *  1.5E-3
        ENDIF
        IF (idt_Im_aq(i) /= 0) THEN
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Im_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Im_aq(i))) +    &
            Cl_tendency(1:kproma) *  7.4E-8/0.545
          IF (idt_Nap_aq(i) /= 0)   &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) +    &
            Cl_tendency(1:kproma) *  7.4E-8/0.545
        ENDIF
        IF (idt_IO3m_aq(i) /= 0) THEN
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt_IO3m_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_IO3m_aq(i))) +    &
            Cl_tendency(1:kproma)  *  2.64E-7/0.545
          IF (idt_Nap_aq(i) /= 0)   &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) +    &
            Cl_tendency(1:kproma)  *  2.64E-7/0.545
        ENDIF
        IF (idt_HCO3m_aq(i) /= 0) THEN
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt_HCO3m_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_HCO3m_aq(i))) +    &
            Cl_tendency(1:kproma)  *  4.2E-3
          IF (idt_Nap_aq(i) /= 0) &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) = &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) +    &
            Cl_tendency(1:kproma)  *  4.2E-3
        ENDIF
#else
        tend_tmp(:,:) = 0.0_dp
          Cl_tendency(1:kproma) =  &
            seasalt_emis(i)%ptr(1:kproma,jrow) * 1.e3 / M_NaCl /  &
            zdz(1:kproma) /rho_air(1:kproma)

        IF (idt_Clm_aq(i) /= 0) THEN
           tend_tmp(1:kproma,nlev) = Cl_tendency(1:kproma)
           CALL mtend_add_l(my_handle, idt_Clm_aq(i) &
                , px=tend_tmp(1:kproma,:))
           IF (idt_Nap_aq(i) /= 0) THEN
              tend_tmp(1:kproma,nlev) = Cl_tendency(1:kproma)
              CALL mtend_add_l(my_handle,  idt_Nap_aq(i) &
                   , px=tend_tmp(1:kproma,:))
           ENDIF
        ENDIF

        IF (idt_Brm_aq(i) /= 0) THEN
           tend_tmp(1:kproma,nlev) = Cl_tendency(1:kproma) *  1.5E-3
           CALL mtend_add_l(my_handle, idt_Brm_aq(i) &
                , px=tend_tmp(1:kproma,:))
           IF (idt_Nap_aq(i) /= 0) &
                CALL mtend_add_l(my_handle, idt_Nap_aq(i) &
                , px=tend_tmp(1:kproma,:))
        ENDIF
        IF (idt_Im_aq(i) /= 0) THEN
           tend_tmp(1:kproma,nlev) = Cl_tendency(1:kproma) *  7.4E-8/0.545
           CALL mtend_add_l(my_handle, idt_Im_aq(i) &
                , px=tend_tmp(1:kproma,:))
          IF (idt_Nap_aq(i) /= 0)   &
               CALL mtend_add_l(my_handle, idt_Nap_aq(i) &
                , px=tend_tmp(1:kproma,:))
         ENDIF
        IF (idt_IO3m_aq(i) /= 0) THEN
           tend_tmp(1:kproma,nlev) = Cl_tendency(1:kproma)  *  2.64E-7/0.545
           CALL mtend_add_l(my_handle, idt_IO3m_aq(i) &
                , px=tend_tmp(1:kproma,:))
          IF (idt_Nap_aq(i) /= 0)   &
               CALL mtend_add_l(my_handle, idt_Nap_aq(i) &
               , px=tend_tmp(1:kproma,:))
        ENDIF
        IF (idt_HCO3m_aq(i) /= 0) THEN
           tend_tmp(1:kproma,nlev) = Cl_tendency(1:kproma)  *  4.2E-3
           CALL mtend_add_l(my_handle, idt_HCO3m_aq(i) &
                , px=tend_tmp(1:kproma,:))
           IF (idt_Nap_aq(i) /= 0) &
                CALL mtend_add_l(my_handle,idt_Nap_aq(i) &
                , px=tend_tmp(1:kproma,:))
        ENDIF
#endif

      ENDIF iftendency
    ENDDO doapn

    CASE(2)

    IF (.not. l_netflux) RETURN

       ! CALCULATE DRY DEPOSITION
       zdp(1:kproma) = pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) - &
            pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev))

       doapn2: DO i=1,APN

          N_new(1:kproma) = Nss_emis(i)%ptr(1:kproma,jrow)  &
            * const/ zdp(1:kproma) *time_step_len

#ifndef MESSYTENDENCY
          N_old(1:kproma) = pxtm1(_RI_X_ZN_(1:kproma,nlev,zidt_num(i)))     &
               + (pxtte(_RI_X_ZN_(1:kproma,nlev,zidt_num(i)))               &
#if defined(ECHAM5)
               + zxtems(1:kproma,zidt_num(i))/ (zdp(1:kproma)/g) &
#endif
               + (Nss_ddep(i)%ptr(1:kproma,jrow) &
                     - Nss_emis(i)%ptr(1:kproma,jrow))&
               * const/ zdp(1:kproma))*time_step_len
          !
#else
          CALL mtend_get_start_l(zidt_num(i),v0=tend_tmp)
          N_old(1:kproma) = tend_tmp(1:kproma,nlev) + ( &
#if defined(ECHAM5)
               zxtems(1:kproma,zidt_num(i))/ (zdp(1:kproma)/g) + &
#endif
               (Nss_ddep(i)%ptr(1:kproma,jrow) &
               - Nss_emis(i)%ptr(1:kproma,jrow))&
               * const/ zdp(1:kproma))*time_step_len
          !
#endif
          !
          N_sum(1:kproma)=N_old(1:kproma)+ N_new(1:kproma)

       WHERE(N_sum(1:kproma) > 1.e-25)
          N_ratio(1:kproma) = N_new(1:kproma)/N_sum(1:kproma)
       ELSEWHERE
          N_ratio(1:kproma) = 0._dp
       ENDWHERE

       IF (idt_Clm_aq(i) /= 0) &
            net_Clm_ddep(1:kproma)=  MAX(0.0_dp, &
               -N_ratio(1:kproma) * Clm_ddep(i)%ptr(1:kproma,jrow) &
               + N_ratio(1:kproma) * seasalt_ddep(i)%ptr(1:kproma,jrow)&
               *1._dp )
       IF (idt_Brm_aq(i) /= 0) &
            net_Brm_ddep(1:kproma)=  MAX(0.0_dp, &
               -N_ratio(1:kproma) * Brm_ddep(i)%ptr(1:kproma,jrow) &
               + N_ratio(1:kproma) * seasalt_ddep(i)%ptr(1:kproma,jrow)&
               *1.5e-3_dp )
       IF (idt_Im_aq(i) /= 0) &
            net_Im_ddep(1:kproma)=  MAX(0.0_dp, &
               -N_ratio(1:kproma) * Im_ddep(i)%ptr(1:kproma,jrow) &
               + N_ratio(1:kproma) * seasalt_ddep(i)%ptr(1:kproma,jrow)&
               *7.4E-8/0.545 )
       IF (idt_IO3m_aq(i) /= 0) &
            net_IO3m_ddep(1:kproma)=  MAX(0.0_dp, &
               -N_ratio(1:kproma) * IO3m_ddep(i)%ptr(1:kproma,jrow) &
               + N_ratio(1:kproma) * seasalt_ddep(i)%ptr(1:kproma,jrow)&
               *2.64E-7/0.54 )
       IF (idt_HCO3m_aq(i) /= 0) &
            net_HCO3m_ddep(1:kproma)=  MAX(0.0_dp, &
               -N_ratio(1:kproma) * HCO3m_ddep(i)%ptr(1:kproma,jrow) &
               + N_ratio(1:kproma) * seasalt_ddep(i)%ptr(1:kproma,jrow)&
               *4.2e-3_dp )
       IF (idt_Hp_aq(i) /= 0) &
            net_Hp_ddep(1:kproma)=  MAX(0.0_dp, &
               -N_ratio(1:kproma) * Hp_ddep(i)%ptr(1:kproma,jrow) &
               + N_ratio(1:kproma) * seasalt_ddep(i)%ptr(1:kproma,jrow)&
               *1.8e-8_dp)

       IF (.not. l_tendency) THEN
#if defined(ECHAM5)

          CONV(1:kproma)  = 1.e-3 *M_air / N_A
          Nap_depsum(1:kproma) = 0._dp

          IF (idt_Clm_aq(i)/= 0) THEN
             zxtems(1:kproma,idt_Clm_aq(i)) = zxtems(1:kproma,idt_Clm_aq(i)) &
                  -  net_Clm_ddep(1:kproma) * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  +  net_Clm_ddep(1:kproma) * CONV(1:kproma)
          ENDIF
          IF (idt_Brm_aq(i) /= 0) THEN
               zxtems(1:kproma,idt_Brm_aq(i)) = zxtems(1:kproma,idt_Brm_aq(i)) &
               -  net_Brm_ddep(1:kproma) * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  +  net_Brm_ddep(1:kproma) * CONV(1:kproma)
          ENDIF
          IF (idt_Im_aq(i) /= 0) THEN
               zxtems(1:kproma,idt_Im_aq(i)) = zxtems(1:kproma,idt_Im_aq(i)) &
               -  net_Im_ddep(1:kproma) * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  +  net_Im_ddep(1:kproma) * CONV(1:kproma)
          ENDIF
          IF (idt_IO3m_aq(i) /= 0) THEN
             zxtems(1:kproma,idt_IO3m_aq(i)) = zxtems(1:kproma,idt_IO3m_aq(i)) &
                  -  net_IO3m_ddep(1:kproma) * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  +  net_IO3m_ddep(1:kproma) * CONV(1:kproma)
          ENDIF

            IF (idt_HCO3m_aq(i) /= 0) THEN
               zxtems(1:kproma,idt_HCO3m_aq(i))= &
                    zxtems(1:kproma,idt_HCO3m_aq(i)) &
                    -  net_HCO3m_ddep(1:kproma) * CONV(1:kproma)
               Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                    +  net_HCO3m_ddep(1:kproma) * CONV(1:kproma)
            ENDIF
            ! mz_rs_20110530: idt_Hp_cs changed to idt_Hp_aq(i)
            IF (idt_Hp_aq(i)/= 0) THEN
               zxtems(1:kproma,idt_Hp_aq(i))= zxtems(1:kproma,idt_Hp_aq(i)) &
                    -  net_Hp_ddep(1:kproma) * CONV(1:kproma)
               Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                    - net_Hp_ddep(1:kproma) * CONV(1:kproma)
            ENDIF
             IF (idt_Nap_aq(i) /= 0) &
                  zxtems(1:kproma,idt_Nap_aq(i))= &
                  zxtems(1:kproma,idt_Nap_aq(i)) &
                 - Nap_depsum(1:kproma)
            !
#endif
       ELSE
          CONV(1:kproma)  = 1.e-3 * M_air /N_a *g / zdp(1:kproma)
          Nap_depsum(1:kproma) = 0._dp

#ifndef MESSYTENDENCY
          IF (idt_Clm_aq(i) /= 0) THEN
             pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Clm_aq(i))) = &
                  pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Clm_aq(i))) &
                  - net_Clm_ddep(1:kproma)  * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  + net_Clm_ddep(1:kproma)  * CONV(1:kproma)
          ENDIF

          IF (idt_Brm_aq(i) /= 0) THEN
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Brm_aq(i))) = &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Brm_aq(i))) &
               - net_Brm_ddep(1:kproma)  * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  + net_Brm_ddep(1:kproma)  * CONV(1:kproma)
          ENDIF
          IF (idt_Im_aq(i) /= 0) THEN
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Im_aq(i))) = &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Im_aq(i))) &
               - net_Im_ddep(1:kproma)  * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  + net_Im_ddep(1:kproma)  * CONV(1:kproma)
          ENDIF
          IF (idt_IO3m_aq(i) /= 0) THEN
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_IO3m_aq(i))) = &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_IO3m_aq(i))) &
               - net_IO3m_ddep(1:kproma)  * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  + net_IO3m_ddep(1:kproma)  * CONV(1:kproma)
          ENDIF
          IF (idt_HCO3m_aq(i) /= 0) THEN
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_HCO3m_aq(i))) = &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_HCO3m_aq(i))) &
               - net_HCO3m_ddep(1:kproma)  * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  + net_HCO3m_ddep(1:kproma)  * CONV(1:kproma)
          ENDIF
          IF (idt_Hp_aq(i) /= 0) THEN
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Hp_aq(i))) = &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Hp_aq(i))) &
               - net_Hp_ddep(1:kproma)  * CONV(1:kproma)
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  + net_Hp_ddep(1:kproma)  * CONV(1:kproma)
          ENDIF
          IF (idt_Nap_aq(i)/=0) &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) = &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt_Nap_aq(i))) - Nap_depsum(1:kproma)
#else
          tend_tmp(:,:) = 0.0_dp
          IF (idt_Clm_aq(i) /= 0) THEN
             tend_tmp(1:kproma,nlev) = - net_Clm_ddep(1:kproma)*CONV(1:kproma)
             CALL mtend_add_l(my_handle,idt_Clm_aq(i) &
                  , px=tend_tmp(1:kproma,:))
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  - tend_tmp(1:kproma,nlev)
          ENDIF

          IF (idt_Brm_aq(i) /= 0) THEN
             tend_tmp(1:kproma,nlev) = - net_Brm_ddep(1:kproma)*CONV(1:kproma)
             CALL mtend_add_l(my_handle,idt_Brm_aq(i) &
                  , px=tend_tmp(1:kproma,:))
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  - tend_tmp(1:kproma,nlev)
          ENDIF
          IF (idt_Im_aq(i) /= 0) THEN
             tend_tmp(1:kproma,nlev) = - net_Im_ddep(1:kproma)*CONV(1:kproma)
             CALL mtend_add_l(my_handle,idt_Im_aq(i) &
                  , px=tend_tmp(1:kproma,:))
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  - tend_tmp(1:kproma,nlev)
          ENDIF
          IF (idt_IO3m_aq(i) /= 0) THEN
             tend_tmp(1:kproma,nlev) = - net_IO3m_ddep(1:kproma)*CONV(1:kproma)
             CALL mtend_add_l(my_handle,idt_IO3m_aq(i) &
                  , px=tend_tmp(1:kproma,:))
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  - tend_tmp(1:kproma,nlev)
          ENDIF
          IF (idt_HCO3m_aq(i) /= 0) THEN
             tend_tmp(1:kproma,nlev) = - net_HCO3m_ddep(1:kproma)*CONV(1:kproma)
             CALL mtend_add_l(my_handle,idt_HCO3m_aq(i) &
                  , px=tend_tmp(1:kproma,:))
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  - tend_tmp(1:kproma,nlev)
          ENDIF
          IF (idt_Hp_aq(i) /= 0) THEN
             tend_tmp(1:kproma,nlev) = - net_Hp_ddep(1:kproma)*CONV(1:kproma)
             CALL mtend_add_l(my_handle,idt_Hp_aq(i) &
                  , px=tend_tmp(1:kproma,:))
             Nap_depsum(1:kproma) = Nap_depsum(1:kproma) &
                  - tend_tmp(1:kproma,nlev)
          ENDIF
          IF (idt_Nap_aq(i)/=0) THEN
             tend_tmp(1:kproma,nlev) = - Nap_depsum(1:kproma)
             CALL mtend_add_l(my_handle,idt_Nap_aq(i) &
                  , px=tend_tmp(1:kproma,:))
          ENDIF
#endif

       ENDIF
    ENDDO doapn2
 CASE DEFAULT
    ! should never be reached
 END SELECT

END SUBROUTINE mecca_aero_vdiff

  ! --------------------------------------------------------------------------

  ! SUBROUTINE "mecca_aero_physc" does not exist
  ! mecca_physc is the chemistry core of this submodul, too.
  ! The following subroutines add to the subroutine mecca_physc
  ! the tracer important for aerosol chemistry
  ! some things are included in the mecca_physc code directly
  ! by the use of the switch l_aero

  ! --------------------------------------------------------------------------

  SUBROUTINE mecca_aero_update_physc(nbl, temp, press, cair, Conc, zmrbc)

    ! calculation of variables needed in the kpp chemistry mechanism

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    ! MECCA
    USE messy_mecca_kpp,          ONLY: fill_k_exf,      fill_k_exb       &
                                      , fill_k_exf_N2O5, fill_k_exf_ClNO3 &
                                      , fill_k_exf_BrNO3,fill_lwc         &
                                      , fill_cvfac,      fill_xaer        &
                                      , NSPEC

    IMPLICIT NONE

    INTRINSIC TRIM

    ! I/O
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_aero_update_physc'
    INTEGER,  INTENT(IN)    :: nbl
    REAL(dp), INTENT(IN)    :: temp(:)
    REAL(dp), INTENT(IN)    :: press(:)
    REAL(dp), INTENT(IN)    :: cair(:)
    REAL(dp), INTENT(IN)    :: zmrbc(:,:)
    REAL(dp), INTENT(INOUT) :: Conc(:,:)

    ! LOCAL
    INTEGER  :: status

    ! define over vector length
    REAL(dp) :: cvfac(nbl,APN)
    REAL(dp) :: k_exf(nbl,APN,NSPEC)
    REAL(dp) :: k_exb(nbl,APN,NSPEC)
    REAL(dp) :: k_exf_N2O5(nbl,APN)
    REAL(dp) :: k_exf_ClNO3(nbl,APN)
    REAL(dp) :: k_exf_BrNO3(nbl,APN)

    ! define in box
    LOGICAL  :: l_het(APN)            ! calculate heterogeneous chemistry
    INTEGER  :: jb, jk, jp

    ALLOCATE(xaer(nbl,APN))
    ALLOCATE(lwc(nbl,APN))
    ALLOCATE(radius(nbl,APN))

    xaer(:,:)   = 0.
    lwc(:,:)    = 0.
    radius(:,:) = 0. ! not necessary, will be set later

    SELECT CASE (TRIM(aerosol_module))
    CASE ('modal','bin')
      CALL modal_lwc(nbl, zmrbc,lwc, cair)
      CALL modal_radii(radius)
    CASE('constant')
      CALL constant_lwc(lwc)
      CALL constant_radii(radius)
    CASE DEFAULT
      CALL error_bi('no aerosol chemistry without any aerosol module',substr)
    ENDSELECT

    CALL mecca_aero_mr2c(Conc, lwc)

    ! calculate coefficients for kpp chemistry
    jb = 0
    DO jk=1,nlev
       DO jp=1,kproma
          jb = jb + 1
          ! calc. switch for aerosol phase
          CALL mecca_aero_xaer( jp, jk, jb, l_het(:), xaer(jb,:), cvfac(jb,:))
          IF (SUM(xaer(jb,:))>0)  THEN
             ! Calc. exchange coefficients: k_exf (forward) and k_exb (backward)
             CALL mecca_aero_calc_k_ex(radius(jb,:), temp(jb), & !IN
               press(jb), l_het(:), xaer(jb,:), lwc(jb,:), Conc(jb,:), & !IN
               k_exf(jb,:,:), k_exb(jb,:,:), & !OUT
               k_exf_N2O5(jb,:), k_exf_ClNO3(jb,:), k_exf_BrNO3(jb,:)) !OUT
          ELSE
             ! set exchange coefficients to 0
             k_exf(jb,:,:)     = 0.
             k_exb(jb,:,:)     = 0.
             k_exf_N2O5(jb,:)  = 0.
             k_exf_ClNO3(jb,:) = 0.
             k_exf_BrNO3(jb,:) = 0.
          ENDIF

       ENDDO
    ENDDO

    ! FILL KPP VECTORs
    CALL fill_lwc (status, lwc)
    IF (status /= 0) CALL error_bi('fill_lwc array size', substr)
    CALL fill_xaer (status, xaer)
    IF (status /= 0) CALL error_bi('fill_xaer array size', substr)
    CALL fill_cvfac (status, cvfac)
    IF (status /= 0) CALL error_bi('fill_cvfac array size', substr)
    CALL fill_k_exf(status, k_exf)
    IF (status /= 0) CALL error_bi('fill_k_exf array size', substr)
    CALL fill_k_exb(status, k_exb)
    IF (status /= 0) CALL error_bi('fill_k_exb array size', substr)
    CALL fill_k_exf_N2O5(status, k_exf_N2O5)
    IF (status /= 0) CALL error_bi('fill_k_exf_N2O5 array size', substr)
    CALL fill_k_exf_ClNO3(status, k_exf_ClNO3)
    IF (status /= 0) CALL error_bi('fill_k_exf_ClNO3 array size', substr)
    CALL fill_k_exf_BrNO3(status, k_exf_BrNO3)
    IF (status /= 0) CALL error_bi('fill_k_exf_BrNO3 array size', substr)

  CONTAINS

    ! ------------------------------------------------------------------------

    SUBROUTINE mecca_aero_xaer(jp, jk, jb, l_het_loc, pxaer, pcvfac)

      ! Determination of the parameters pxaer(mode)
      ! if liquid water exists in these modes they will be 1 otherwise zero

      ! ECHAM5/MESSy
      USE messy_main_data_bi,       ONLY: slf
      ! MESSy
      USE messy_main_constants_mem, ONLY: N_A

      IMPLICIT NONE

      INTRINSIC :: REAL

      ! LOCAL
      INTEGER,  INTENT(IN)   :: jp, jk, jb
      LOGICAL,  INTENT(OUT)  :: l_het_loc(:)
      REAL(dp), INTENT(OUT)  :: pxaer(:)
      REAL(dp), INTENT(OUT)  :: pcvfac(:)
      INTEGER :: zkc

      pxaer(:)     = 0
      pcvfac(:)    = 0.
      l_het_loc(:) =.FALSE.

      DO zkc = 1, APN
         IF ((lwc(jb,zkc)<=1.E-12_dp).OR.(ind_H2O_a(zkc)==0)&
              .OR.(REAL(jk,dp)<heightidx(jp,jrow)) &
              .OR.(slf(jp,jrow)>=0.2_dp)) THEN
            pxaer(zkc)  = 0
            pcvfac(zkc) = 0.
         ELSE
            pxaer(zkc)  = 1
            ! cvfac: conversion factor dm^3(aq)/mol => cm^3(air)/molecule
            pcvfac(zkc) = 1.e3 / (N_A * lwc(jb,zkc))
             l_het_loc(zkc) =.TRUE.
         ENDIF
      ENDDO

    END SUBROUTINE mecca_aero_xaer

    !-------------------------------------------------------------------------

    SUBROUTINE constant_lwc(lwc)

      IMPLICIT NONE

      ! LOCAL
      REAL(dp) :: lwc(:,:)
      INTEGER  :: i

      DO i=1,APN
        lwc(:,i) = 1.e-11*i/(APN)
      ENDDO

    END SUBROUTINE constant_lwc

    !-------------------------------------------------------------------------

    SUBROUTINE constant_radii(pradius)

      IMPLICIT NONE

      ! I/O
      REAL(dp) :: pradius(:,:)
      ! LOCAL
      INTEGER  :: i

      DO i=1,APN
        pradius(:,i) = 1.e-6*i/(APN)
      ENDDO

    END SUBROUTINE constant_radii

    !-------------------------------------------------------------------------

    SUBROUTINE modal_lwc(nbl, zmrbc, lwc, cair)

      ! ECHAM5/MESSy
      USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev
      ! MESSy
      USE messy_main_constants_mem, ONLY: N_A, pi, M_air

      IMPLICIT NONE

      INTEGER,  INTENT(IN) :: nbl     ! block length
      REAL(dp), INTENT(in) :: zmrbc(:,:) ! tracer mixing ratio
      REAL(dp), INTENT(IN) :: cair(:)
      REAL(dp), INTENT(out) :: lwc(:,:)

      REAL(dp) :: rho(nbl)
      REAL(dp) :: zradn_e3
      REAL(dp) :: zntot
      REAL(dp) :: zsup_script
      REAL(dp) :: zrad

      INTEGER :: i, jp, jk, jb, jm

      INTRINSIC  LOG10, EXP, MAX

      ! Liquid water content [m^3(aq)/m^3(air)] calculation for log-normal
      ! particle distribution see  Sander (1999)
      ! LWC= 4/3 pi R^3_N N_{tot} exp( 9* (lg sigma)^2/ (2(log e)^2))

      SELECT CASE (TRIM(ti_gp(zidt_num(1))%tp%ident%unit))
      CASE('mol(part)/mol')
         rho(:) = 1.E6_dp * cair(:) !/ N_A * N_A! [mol(air)/m3]
      CASE('1/mol')
        rho(:) = 1.E6_dp * cair(:) / N_A ! [mol(air)/m3]
      CASE('1/kg')
        rho(:) = 1.E3_dp * M_air * cair(:) / N_A ! [kg(air)/m3]
      CASE DEFAULT
      ENDSELECT
      ! rho = press / (R *temp*(1.0_dp+(M_air / M_H2O -1._dp)*sphum))

      DO i=1, APN
         jb = 0
         jm = modenumber(i) ! um_ak_20100802
         DO jk=1, nlev
            DO jp = 1, kproma
               jb = jb + 1

               zntot = zmrbc(jb,zidt_num(i)) *rho(jb)

               zrad = (wetradius(_RI_XYZN_(jp,jrow,jk,jm))*wetradius(_RI_XYZN_(jp,jrow,jk,jm))  &
                                             *wetradius(_RI_XYZN_(jp,jrow,jk,jm)))  &
                    - (dryradius(_RI_XYZN_(jp,jrow,jk,jm))*dryradius(_RI_XYZN_(jp,jrow,jk,jm))  &
                                             *dryradius(_RI_XYZN_(jp,jrow,jk,jm)))

               zradn_e3 = MAX(0.0_dp,zrad)

               zsup_script = 9.0_dp* LOG10(sigma_aerosol(modenumber(i)))* &
                    LOG10(sigma_aerosol(modenumber(i)))/&
                    (2._dp*LOG10(EXP(1._dp))*LOG10(EXP(1._dp)))

               lwc(jb,i)=4._dp/3. * pi * zradn_e3 * zntot &
                              * EXP(zsup_script)
               lwc_3d(i)%ptr(jp,jk,jrow)=lwc(jb,i)

            ENDDO
         ENDDO
      ENDDO

    END SUBROUTINE modal_lwc

    !-------------------------------------------------------------------------

    SUBROUTINE modal_radii(pradius)

      IMPLICIT NONE

      INTRINSIC LOG10, LOG, EXP

      ! I/O
      REAL(dp) :: pradius(:,:)
      ! LOCAL
      INTEGER  :: I, jp ,jk, jb, jm

      ! ak: pradius corresponds to eq. 108 (Sander,1999) => kmt is
      ! proportional to surface area
      !    R_A = R_N * 10 ^ ( 2* (sigma)^2/ lg e )
      !    10 ^ x = exp( x* ln(10))
      DO i=1,APN
         jb = 0
         jm = modenumber(i) ! um_ak_20100802
         DO jk=1, nlev
            DO jp = 1, kproma
               jb = jb +1
               pradius(jb,i)=wetradius(_RI_XYZN_(jp,jrow,jk,jm))* &
                    EXP((2._dp*LOG10(sigma_aerosol(modenumber(i)))*&
                    LOG10(sigma_aerosol(modenumber(i)))/ &
                    LOG10(EXP(1._dp)))*LOG(10._dp))
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE modal_radii

    ! ------------------------------------------------------------------------

  END SUBROUTINE mecca_aero_update_physc

  ! --------------------------------------------------------------------------

  SUBROUTINE mecca_aero_mr2c(Conc, lwc)

    ! MESSy
    USE messy_main_constants_mem,   ONLY: N_A, M_H2O, rho_H2O

    IMPLICIT NONE

    REAL(dp), INTENT(INOUT) :: Conc(:,:)
    REAL(dp), INTENT(IN)    :: lwc(:,:)

    ! LOCAL
    INTEGER :: i
    !CHARACTER(len=*), PARAMETER :: substr=' mecca_aero_mr2c'
    DO i=1,APN
      IF (ind_H2O_a(i)>0) Conc(:,ind_H2O_a(i)) = &
        (rho_H2O/M_H2O) * N_A * 1.E-3 * lwc(:,i)
    ENDDO

  END SUBROUTINE mecca_aero_mr2c

  !---------------------------------------------------------------------------

  SUBROUTINE mecca_aero_save_pH(Conc)

    ! ECHAM
    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev
    ! MESSy
    USE messy_main_constants_mem,  ONLY: N_A

    IMPLICIT NONE
    ! LOCAL
    INTEGER :: I, jb, jp ,jk

    REAL(dp), INTENT(IN) :: Conc(:,:)

    INTRINSIC log10

    DO i=1,APN
      IF (ind_Hp_a(i)/= 0 ) THEN
         jb = 0
         DO jk = 1, nlev
            DO jp = 1, kproma
               jb = jb + 1
               IF( Conc(jb,ind_Hp_a(i)) .GT. TINY(Conc(1,1)) &
                    .AND. xaer(jb,i)==1 .AND. &
                    lwc(jb,i).GT. TINY(lwc(1,i)))  &

                    pH(i)%ptr(_RI_XYZ__(jp,jrow,jk)) = &
                    -LOG10(Conc(jb,ind_Hp_a(i))/N_A*1.e3/lwc(jb,i))

            ENDDO
         ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE mecca_aero_save_pH

  SUBROUTINE mecca_aero_diag_si(zmr, zmrac,nbl)

    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp
    USE messy_mecca_aero,         ONLY: mecca_aero_diag

    IMPLICIT NONE

    INTEGER,  INTENT(IN)    :: nbl
    REAL(dp), INTENT(IN)    :: zmr(nbl,ntrac_gp)
    REAL(dp), INTENT(INOUT) :: zmrac(nbl, ntrac_gp)

    !LOCAL
    INTEGER :: jb

    IF (idt_Brsalt/=0 .AND. idt_Brorg/=0 .AND. idt_BrSScap/=0) THEN
       DO jb=1,nbl
          CALL mecca_aero_diag(zmr(jb,idt_Brsalt), zmr(jb,idt_Brorg) &
               , zmr(jb,idt_BrSScap), zmrac(jb,idt_Brsalt)           &
               , zmrac(jb,idt_Brorg), zmrac(jb,idt_BrSScap))
       ENDDO
    ENDIF

  END SUBROUTINE mecca_aero_diag_si

  SUBROUTINE mecca_aero_dealloc

    DEALLOCATE(radius)
    DEALLOCATE(lwc)
    DEALLOCATE(xaer)

  END SUBROUTINE mecca_aero_dealloc


  !---------------------------------------------------------------------------
  ! PRIVATE ROUTINES
  !---------------------------------------------------------------------------

  SUBROUTINE mecca_aero_read_nml_cpl(status, iou)

    ! MECCA MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    ! read namelist for 'coupling' to online tracers
    ! Author: Astrid Kerkweg, MPICH, Nov 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_aero_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    INTRINSIC TRIM

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    lcpl_aerophase(:)  = .FALSE.
    aerosol_module     = ''
    number_main(:)     = ''
    number_sub(:)      = ''
    aerosol_channel    = ''
    aerosol_wetradius  = ''
    aerosol_dryradius  = ''
    sigma              = ''
    modenumber(:)      = 0
    l_calc_emis     = .FALSE.
    l_netflux       = .FALSE.
    height_channel  = ''
    height_index    = ''
    emis_channel    = ''
    emis_ss(:)      = ''
    ddep_channel    = ''
    ddep_ss(:)      = ''
    l_tendency      = .FALSE.

    CALL read_nml_open(lex, substr, iou, 'CPL_AERO', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_AERO, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_AERO', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

    CALL info_bi('******************************************************')
    SELECT CASE (TRIM(aerosol_module))
    CASE ('modal')
      CALL info_bi('A modal aerosol module is choosen for aerosol chemistry')
    CASE ('bin')
      CALL info_bi('A bin aerosol module is choosen for aerosol chemistry')
    CASE ('constant')
      CALL info_bi('Constant liquid water profile instant')
      CALL info_bi('of an aerosol module choosen')
    CASE DEFAULT
      CALL warning_bi('No modal, bin or constant aerosol module chosen',substr)
      CALL warning_bi('Implementation of other aerosol module types is'//&
        'not available',substr)
      RETURN ! error
    END SELECT
    CALL info_bi('******************************************************')
    IF (l_calc_emis) THEN
      CALL info_bi('Calculate ion emission from')
      IF (l_netflux) THEN
        CALL info_bi('  net flux (emission - dry deposition)')
      ELSE
        CALL info_bi('  emission flux')
      ENDIF
    ELSE
      CALL info_bi('No ion emission calculation.')
    ENDIF
    CALL info_bi('******************************************************')

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE mecca_aero_read_nml_cpl

  !---------------------------------------------------------------------------
  SUBROUTINE initialize_indexarrays_trac

    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp
    USE messy_main_tracer,        ONLY: AEROSOL, I_aerosol_mode

    IMPLICIT NONE

    !LOCAL
    INTEGER :: i, jt

    ! INITIALIZE TRACER INDICES
    ! Clm
    idt_Hp_aq(:)    = 0
    idt_Clm_aq(:)   = 0
    idt_Brm_aq(:)   = 0
    idt_Im_aq(:)    = 0
    idt_IO3m_aq(:)  = 0
    idt_HCO3m_aq(:) = 0
    idt_Nap_aq(:)   = 0

    DO jt=1, ntrac_gp
      IF (ti_gp(jt)%tp%ident%medium /= AEROSOL) CYCLE
      DO i=1,APN
        IF (ti_gp(jt)%tp%meta%cask_i(I_aerosol_mode) /= modenumber(i)) CYCLE
        SELECT CASE(TRIM(ti_gp(jt)%tp%ident%basename))
        CASE ('Hp')
          idt_Hp_aq(i) = jt
        CASE ('Nap')
          idt_Nap_aq(i) = jt
        CASE ('HCO3m')
          idt_HCO3m_aq(i) = jt
        CASE ('Clm')
          idt_Clm_aq(i) = jt
        CASE ('Brm')
          idt_Brm_aq(i) = jt
        CASE ('Im')
          idt_Im_aq(i) = jt
        CASE ('IO3m')
          idt_IO3m_aq(i) = jt
        END SELECT
      ENDDO
    ENDDO

  END SUBROUTINE initialize_indexarrays_trac

!*****************************************************************************
  !***************************************************************************

  SUBROUTINE mecca_aero_init_gasaq(l_print)

    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, ti_gp
    USE messy_main_tracer,        ONLY: R_Henry_T0, R_Henry_Tdep, R_alpha_T0 &
                                      , R_alpha_Tdep, R_molarmass            &
                                      , get_tracer
    USE messy_main_constants_mem, ONLY: HLINE2
    USE messy_mecca_aero,         ONLY: Henry_T0, Henry_Tdep, molar_mass &
                                      , alpha_T0, alpha_Tdep
    USE messy_mecca_kpp,          ONLY: NSPEC, SPC_NAMES

    IMPLICIT NONE

    INTEGER :: icp ! tracer index
    INTEGER :: jn  ! species index
    INTEGER :: status

    LOGICAL, INTENT(IN), OPTIONAL :: l_print

    IF (PRESENT(l_print)) THEN
      IF (l_print) THEN
      PRINT *, HLINE2
      PRINT *, "         Henry's law coefficients "// &
        "and accommodation coefficients"
      PRINT *, HLINE2
      PRINT *, 'species           Henry_T0 Henry_Tdep'// &
        '   alpha_T0 alpha_Tdep         M'
      PRINT *, '                   [M/atm]        [K]'// &
        '        [1]        [K]   [g/mol]'
      PRINT *, HLINE2
      ENDIF
    ENDIF
    DO jn = 1,NSPEC
       CALL get_tracer(status, GPTRSTR, TRIM(SPC_NAMES(jn)) &
            , idx=icp )
       IF (status /= 0) CYCLE

       Henry_T0(jn)   = ti_gp(icp)%tp%meta%cask_r(R_Henry_T0)
       Henry_Tdep(jn) = ti_gp(icp)%tp%meta%cask_r(R_Henry_Tdep)
       alpha_T0(jn)   = ti_gp(icp)%tp%meta%cask_r(R_alpha_T0)
       alpha_Tdep(jn) = ti_gp(icp)%tp%meta%cask_r(R_alpha_Tdep)
       molar_mass(jn) = ti_gp(icp)%tp%meta%cask_r(R_molarmass) / 1000. ![kg/mol]

       IF (PRESENT(l_print)) THEN
          IF (l_print) THEN
             WRITE(*,'(1X,A15)', ADVANCE='NO') SPC_NAMES(jn)
             IF (Henry_T0(jn)>=0.) THEN
                WRITE(*,'(ES11.2)', ADVANCE='NO') Henry_T0(jn)
             ELSE
                WRITE(*,'(A)', ADVANCE='NO') '   --------'
             ENDIF
             IF (Henry_Tdep(jn)>=0.) THEN
                WRITE(*,'(F11.0)', ADVANCE='NO') Henry_Tdep(jn)
             ELSE
                WRITE(*,'(A)', ADVANCE='NO') '     ------'
             ENDIF
             IF (alpha_T0(jn)>=0.) THEN
                WRITE(*,'(ES11.2)', ADVANCE='NO') alpha_T0(jn)
             ELSE
                WRITE(*,'(A)', ADVANCE='NO') '   --------'
             ENDIF
             IF (alpha_Tdep(jn)>=0.) THEN
                WRITE(*,'(F11.0)', ADVANCE='NO') alpha_Tdep(jn)
             ELSE
                WRITE(*,'(A)', ADVANCE='NO') '     ------'
             ENDIF
             WRITE(*,'(F10.2)') 1000.*molar_mass(jn) ! [g/mol]
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE mecca_aero_init_gasaq
  !********************************************

END MODULE messy_mecca_aero_si

!*****************************************************************************
