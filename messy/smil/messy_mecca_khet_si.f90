#include "messy_main_ppd_bi.inc"

!*****************************************************************************
!                Time-stamp: <2020-11-27 17:48:49 b302010>
!*****************************************************************************

! Authors:
! Rolf Sander, MPICH, 2005-...
! The code is partially based on the PSC module by Joachim Buchholz

MODULE messy_mecca_khet_si

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      info_bi, warning_bi, error_bi
  ! MESSy
  USE messy_mecca,              ONLY: l_force_khet
  USE messy_mecca_khet          ! gamma_*, liq_*, k_het_*, k_mt, ...
  USE messy_mecca_kpp,          ONLY: dp, IHS_MAX, IHT_MAX, iht_N2O5 &
                                    , iht_HNO3, iht_Hg, iht_RGM
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_3D_ARRAY, &
                                      PTR_4D_ARRAY, str
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_VLONG, &
                                      pi, M_H2O, M_air, R_gas, MHg

  IMPLICIT NONE
  INTRINSIC :: NULL

  PRIVATE
  PUBLIC :: mecca_khet_initialize    ! read namelists
  PUBLIC :: mecca_khet_init_memory   ! create new channel objects
  PUBLIC :: mecca_khet_init_coupling ! get tracers/channels from submodels
  PUBLIC :: mecca_khet_physc_gp
!!#D attila +
  PUBLIC :: mecca_khet_physc_lg
!!#D attila -
  PUBLIC :: mecca_khet_free_memory

  ! channel objects
  TYPE(PTR_3D_ARRAY), PUBLIC, DIMENSION(:,:), SAVE, POINTER :: khet_Tr_3d
  TYPE(PTR_1D_ARRAY), PUBLIC, DIMENSION(:,:), SAVE, POINTER :: khet_Tr_1d
  TYPE(PTR_3D_ARRAY), PUBLIC, DIMENSION(IHS_MAX), SAVE      :: khet_St_3d
  TYPE(PTR_1D_ARRAY), PUBLIC, DIMENSION(IHS_MAX), SAVE      :: khet_St_1d
  CHARACTER(LEN=9), PUBLIC, DIMENSION(IHT_MAX) :: khet_Tr_name = (/ &
    'N2O5     ', 'HNO3     ', 'Hg       ', 'RGM      ' /)
  CHARACTER(LEN=9), PUBLIC, DIMENSION(IHS_MAX) :: khet_St_name = (/ &
    'N2O5_H2O ', 'HOCl_HCl ', 'ClNO3_HCl', &
    'ClNO3_H2O', 'N2O5_HCl ', 'ClNO3_HBr', &
    'BrNO3_HCl', 'HOCl_HBr ', 'HOBr_HCl ', &
    'HOBr_HBr ', 'BrNO3_H2O', 'Hg       ', 'RGM      ' /)
  ! mz_pj_20070305+
  ! FLAG FOR STRATOSPHERE: 0: troposphere, 1: stratosphere
  REAL(DP), DIMENSION(:,:,:), SAVE, POINTER :: r_strat_flag => NULL()
  REAL(DP), DIMENSION(:),     SAVE, POINTER :: r_strat_flag_1d => NULL()
  ! mz_pj_20070305-

  ! CPL_KHET NAMELIST
  CHARACTER (LEN=STRLEN_MEDIUM), PUBLIC ::  &
    strat_channel    = ''
  CHARACTER (LEN=STRLEN_MEDIUM), DIMENSION(2), PUBLIC :: &
    aerosurf_clim  = ''
  INTEGER, PUBLIC :: xsm_cpl ! aerosol submodel used for coupling to MECCA

  ! climatological aerosol surface:
  REAL(dp), DIMENSION(:,:,:), POINTER :: surface_3d   => NULL()
  REAL(dp), DIMENSION(:), POINTER     :: surface_1d   => NULL()

  ! max. number of aerosol submodels in namelist:
  INTEGER, PARAMETER :: NMAXSM = 10

  TYPE (PTR_4D_ARRAY), DIMENSION(:), POINTER, SAVE :: wetradius => NULL()
  TYPE (PTR_1D_ARRAY), DIMENSION(:), POINTER, SAVE :: sigma     => NULL()

  TYPE SUBMODEL_INFO
    CHARACTER(LEN=STRLEN_MEDIUM)   :: name = ''
    INTEGER                        :: nmodes = 0
    INTEGER, DIMENSION(:), POINTER :: mode => NULL()
    INTEGER, DIMENSION(:), POINTER :: idt_N => NULL() ! index in tracer field
  END TYPE SUBMODEL_INFO

  INTEGER :: nsm
  TYPE(SUBMODEL_INFO), DIMENSION(NMAXSM), SAVE :: xsm

  LOGICAL :: l_climatology

CONTAINS

  !***************************************************************************

  SUBROUTINE mecca_khet_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    ! MESSy
    USE messy_main_tools,        ONLY: find_next_free_unit
    USE messy_mecca_khet,        ONLY: mecca_khet_read_nml_ctrl

    IMPLICIT NONE

    ! local
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_khet_initialize'
    INTEGER :: iou    ! I/O unit
    INTEGER :: status ! error status
    INTEGER :: js, jm

    ! read CTRL_KHET namelist:
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL mecca_khet_read_nml_ctrl(status, iou)
      IF (status /= 0)  CALL error_bi('error in mecca_khet_read_nml_ctrl',substr)
    ENDIF
    ! broadcast results:
    CALL p_bcast (l_troposphere,  p_io)
    CALL p_bcast (l_stratosphere, p_io)

    ! read CPL_KHET namelist:
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL mecca_khet_read_nml_cpl(status, iou)
      IF (status /= 0) CALL error_bi('error in mecca_khet_read_nml_cpl',substr)
    END IF
    ! broadcast results:
    CALL p_bcast(strat_channel,    p_io)
    CALL p_bcast(aerosurf_clim(1), p_io)
    CALL p_bcast(aerosurf_clim(2), p_io)
    CALL p_bcast(nsm,              p_io)
    CALL p_bcast(xsm_cpl,          p_io)
    DO js=1, nsm
      CALL p_bcast(xsm(js)%name, p_io)
      CALL p_bcast(xsm(js)%nmodes, p_io)
      IF (.NOT.p_parallel_io) ALLOCATE(xsm(js)%mode(xsm(js)%nmodes))
      DO jm=1, xsm(js)%nmodes
        CALL p_bcast(xsm(js)%mode(jm), p_io)
      END DO
    END DO

  END SUBROUTINE mecca_khet_initialize

  !***************************************************************************

  SUBROUTINE mecca_khet_init_memory(L_LG, L_GP)

    ! BML/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
#if defined(ECHAM5)
    USE messy_main_channel_bi, ONLY:  LG_ATTILA
#endif
    ! MESSy
    USE messy_main_channel,    ONLY: new_channel, new_channel_object, &
                                     new_attribute

    IMPLICIT NONE
    INTRINSIC :: TRIM

    LOGICAL, INTENT(IN)  :: L_GP
    LOGICAL, INTENT(IN)  :: L_LG

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'mecca_khet_init_memory'
    INTEGER                      :: status, i, js
    CHARACTER(LEN=STRLEN_MEDIUM) :: name, smname

    CALL start_message_bi(submodstr, 'CHANNEL DEFINITION', substr)

    ! we need them anyhow, in case of GP or LG!!!
    ALLOCATE(wetradius (nsm))
    ALLOCATE(sigma     (nsm))
    DO js=1, nsm
      ALLOCATE(xsm(js)%idt_N(xsm(js)%nmodes))
      xsm(js)%idt_N(:) = 0
    ENDDO
    
    IF (L_GP) THEN
        IF ((l_force_khet).OR.(l_troposphere)) THEN
         ALLOCATE(khet_Tr_3d(0:nsm,IHT_MAX))
         CALL new_channel(status, submodstr//'_gp', reprid=GP_3D_MID)
         CALL channel_halt(substr, status)
         ! troposheric rate coefficients
         DO js=0, nsm
           IF (js==0) THEN
             smname = 'clim'             ! climatology
           ELSE
             smname = TRIM(xsm(js)%name) ! submodel name
           ENDIF
           DO i=1, IHT_MAX
             name = 'khet_'//TRIM(khet_Tr_name(i))//'_'//TRIM(smname)
             CALL new_channel_object(status, submodstr//'_gp', name, &
               p3=khet_Tr_3d(js,i)%PTR)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, submodstr//'_gp', name, &
               'long_name', c=name)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, submodstr//'_gp', name, 'units', c='1/s')
             CALL channel_halt(substr, status)
             CALL info_bi('channel/object '//submodstr//'/'//TRIM(name)//' was created')
           ENDDO
         ENDDO
        END IF
    ENDIF

#if defined(ECHAM5)
    IF (L_LG) THEN
        CALL new_channel(status, submodstr//'_lg', reprid=LG_ATTILA)
        CALL channel_halt(substr, status)
        ! stratosheric rate coefficients
        IF ((l_force_khet).OR.(l_stratosphere)) THEN ! stratosphere
           DO i=1, IHS_MAX
             name = 'khet_'//TRIM(khet_St_name(i))//'_STRAT'
             CALL new_channel_object(status, submodstr//'_lg', name, &
               p1=khet_St_1d(i)%PTR)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, submodstr//'_lg', name, &
               'long_name', c=name)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, submodstr//'_lg', name, 'units', c='1/s')
             CALL channel_halt(substr, status)
             CALL info_bi('channel/object '//submodstr//'/'//TRIM(name)//' was created')
           ENDDO
        END IF
        IF ((l_force_khet).OR.(l_troposphere)) THEN
        ! troposheric rate coefficients
        ALLOCATE(khet_Tr_1d(0:nsm,IHT_MAX))
        DO js=0, nsm
          IF (js==0) THEN
            smname = 'clim'             ! climatology
          ELSE
            smname = TRIM(xsm(js)%name) ! submodel name
          ENDIF
          DO i=1, IHT_MAX
            name = 'khet_'//TRIM(khet_Tr_name(i))//'_'//TRIM(smname)
            CALL new_channel_object(status, submodstr//'_lg', name, &
              p1=khet_Tr_1d(js,i)%PTR)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, submodstr//'_lg', name, &
              'long_name', c=name)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, submodstr//'_lg', name, 'units', c='1/s')
            CALL channel_halt(substr, status)
            CALL info_bi('channel/object '//submodstr//'/'//TRIM(name)//' was created')
          ENDDO
        ENDDO
        END IF
     ENDIF
#endif
    CALL end_message_bi(submodstr,'CHANNEL DEFINITION',substr)

  END SUBROUTINE mecca_khet_init_memory

  !***************************************************************************

  SUBROUTINE mecca_khet_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, ntrac_gp, ti_gp
    ! MESSy
    USE messy_main_tracer,        ONLY: get_tracer, NUMBERDENSITY &
                                      , I_aerosol_mode
    USE messy_main_channel,       ONLY: get_channel_object

    IMPLICIT NONE
    INTRINSIC :: TRIM

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_khet_init_coupling'
    INTEGER :: status
    INTEGER :: js, jm, jt, i
    LOGICAL :: coupling_error
    ! mz_pj_20080801+
    CHARACTER(LEN=STRLEN_MEDIUM+3) :: z_aerosol_channel = ''
    ! mz_pj_20080801-

    CALL start_message_bi(submodstr,'INIT COUPLING',substr)

    ! ------------------------------------------------------------------------

    ! troposheric aerosol properties
    IF ((l_force_khet).OR.(l_troposphere)) THEN

      CALL info_bi('***** AEROSOL CLIMATOLOGY')
      CALL get_channel_object(status, &
        TRIM(aerosurf_clim(1)), TRIM(aerosurf_clim(2)), p3=surface_3d)
      IF (status == 0) THEN
        l_climatology = .TRUE.
        CALL info_bi('    channel/object: '// &
          TRIM(aerosurf_clim(1))//'/'//TRIM(aerosurf_clim(2)))
      ELSE
        l_climatology = .FALSE.
        CALL warning_bi(TRIM(aerosurf_clim(1))//'/'//TRIM(aerosurf_clim(2))//&
          ' not available',substr)
      ENDIF

      ! number density (from tracer)
      tracer_loop: DO jt=1, ntrac_gp
        IF (ti_gp(jt)%tp%ident%quantity == NUMBERDENSITY) THEN
          submodel_loop: DO js=1, nsm
            IF (TRIM(ti_gp(jt)%tp%ident%submodel) == TRIM(xsm(js)%name)) THEN
              mode_loop: DO jm=1, xsm(js)%nmodes
                IF (ti_gp(jt)%tp%meta%cask_i(I_aerosol_mode) == xsm(js)%mode(jm)) THEN
                  xsm(js)%idt_N(jm) = jt
                ENDIF
              ENDDO mode_loop
            ENDIF
          ENDDO submodel_loop
        ENDIF
      ENDDO tracer_loop

      ! wetradius, and sigma (from channel objects)
      coupling_error = .FALSE.
      submodel_loop2: DO js=1, nsm
        CALL info_bi('***** TROPOSHERIC SUBMODEL: '//TRIM(xsm(js)%name))
        mode_loop2: DO jm=1, xsm(js)%nmodes
          CALL info_bi('mode number:'//TRIM(str(xsm(js)%mode(jm))))
          IF (xsm(js)%idt_N(jm)==0) THEN
            coupling_error = .TRUE.
            CALL warning_bi('number tracer for '//TRIM(xsm(js)%name)//&
              ' not available',substr)
          ELSE
            CALL info_bi('  number tracer:  '// &
              ti_gp(xsm(js)%idt_N(jm))%tp%ident%fullname)
          ENDIF
          ! mz_pj_20080801+
          z_aerosol_channel = TRIM(xsm(js)%name)//'_'//GPTRSTR
          ! mz_pj_20080801-
          ! wetradius:
! mz_pj_20080801+
!!$          CALL get_channel_object(status, &
!!$            TRIM(xsm(js)%name), 'wetradius', p4=wetradius(js)%ptr)
!!$          IF(status == 0) THEN
!!$            CALL info_bi('  channel/object: '//TRIM(xsm(js)%name)//'/wetradius')
!!$          ELSE
!!$            coupling_error = .TRUE.
!!$            CALL warning_bi(TRIM(xsm(js)%name)//&
!!$              '/wetradius not available',substr)
!!$          ENDIF
          CALL get_channel_object(status, &
            TRIM(z_aerosol_channel), 'wetradius', p4=wetradius(js)%ptr)
          IF(status == 0) THEN
            CALL info_bi('  channel/object: '//TRIM(z_aerosol_channel)//'/wetradius')
          ELSE
            coupling_error = .TRUE.
            CALL warning_bi(TRIM(z_aerosol_channel)//&
              '/wetradius not available',substr)
          ENDIF
! mz_pj_20080801-
          ! sigma:
! mz_pj_20080801+
!!$          CALL get_channel_object(status, &
!!$            TRIM(xsm(js)%name), 'sigma', p1=sigma(js)%ptr)
!!$          IF(status == 0) THEN
!!$            CALL info_bi('  channel/object: '//TRIM(xsm(js)%name)//'/sigma')
!!$          ELSE
!!$            coupling_error = .TRUE.
!!$            CALL warning_bi(TRIM(xsm(js)%name)//&
!!$              '/sigma not available',substr)
!!$          ENDIF
          CALL get_channel_object(status, &
            TRIM(z_aerosol_channel), 'sigma', p1=sigma(js)%ptr)
          IF(status == 0) THEN
            CALL info_bi('  channel/object: '//TRIM(z_aerosol_channel)//'/sigma')
          ELSE
            coupling_error = .TRUE.
            CALL warning_bi(TRIM(z_aerosol_channel)//&
              '/sigma not available',substr)
          ENDIF
! mz_pj_20080801-
        END DO mode_loop2
      END DO submodel_loop2

      IF (coupling_error) CALL error_bi('coupling error',substr)

    ENDIF

    ! ------------------------------------------------------------------------

    ! stratospheric rate coefficients
    IF ((l_force_khet).OR.(l_stratosphere)) THEN
      DO i=1, IHS_MAX
        CALL get_channel_object(status, TRIM(strat_channel), &
          'khet_'//TRIM(khet_St_name(i)), p3=khet_St_3d(i)%ptr)
        IF (status /= 0) &
          CALL error_bi('khet_'//TRIM(khet_St_name(i))// &
          ' not found in channel '//TRIM(strat_channel), substr)
      END DO
      ! mz_pj_20070305+
      CALL get_channel_object(status, TRIM(strat_channel), &
           'STRAT_region', p3=r_strat_flag)
      IF (status /= 0) &
           CALL error_bi('STRAT_region'// &
           ' not found in channel '//TRIM(strat_channel), substr)
      ! mz_pj_20070305-

    ENDIF

    ! ------------------------------------------------------------------------

    CALL end_message_bi(submodstr,'INIT COUPLING',substr)

  END SUBROUTINE mecca_khet_init_coupling

  !***************************************************************************

  SUBROUTINE mecca_khet_physc_gp

    ! ECHAM5/MESSy
    USE messy_main_data_bi,       ONLY: press_3d, tm1_3d, tte_3d
    USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma
    USE messy_main_timer,         ONLY: time_step_len
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, EXP, NINT, LOG10, MAX

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_khet_physc_gp'

    INTEGER :: zphase, jp, jk, jm, js, i
    REAL(dp) :: T, p, zsurface, zr_s, zN_s0, zn
    REAL(dp) :: kh1, kh2, kh3, kh4
    REAL(dp), DIMENSION(30) :: kh

    !-------------------------------------------------------------------------

    level_loop: DO jk= 1, nlev
      kproma_loop: DO jp = 1, kproma

        T = tm1_3d(_RI_XYZ__(jp,jrow,jk)) + tte_3d(_RI_XYZ__(jp,jrow,jk)) * time_step_len
! mz_pj_20080930+
!        p = press_3d(jp,jk,jrow)/100.0_dp ! [hPa] !qqq why not Pa?
        p = press_3d(_RI_XYZ__(jp,jrow,jk)) ! [Pa] 
! mz_pj_20080930-

        !---------------------------------------------------------------------

        if_l_troposphere: IF (l_troposphere) THEN ! tropospheric aerosol

          IF (l_climatology) THEN
            zsurface = surface_3d(_RI_XYZ__(jp,jrow,jk)) ! [m2/m3]
            khet_Tr_3d(0,iht_N2O5)%ptr(_RI_XYZ__(jp,jrow,jk)) = & ! N2O5 uptake
              k_mt(M_N2O5,T,zsurface,gamma_N2O5)
            khet_Tr_3d(0,iht_HNO3)%ptr(_RI_XYZ__(jp,jrow,jk)) = & ! HNO3 uptake
              k_mt(M_HNO3,T,zsurface,gamma_HNO3)
            khet_Tr_3d(0,iht_Hg)%ptr(_RI_XYZ__(jp,jrow,jk)) = &   ! Hg (GEM) uptake
              k_mt(MHg,T,zsurface,gamma_Hg)
            khet_Tr_3d(0,iht_RGM)%ptr(_RI_XYZ__(jp,jrow,jk)) = &  ! Hg++ (RGM) uptake
              k_mt(MHg,T,zsurface,gamma_RGM)          ! assuming M_RGM=MHg
          ENDIF

          submodel_loop: DO js=1, nsm
            kh1 = 0._dp
            kh2 = 0._dp
            kh3 = 0._dp
            kh4 = 0._dp
            mode_loop: DO jm = 1, xsm(js)%nmodes
              ! zn = number density [1/mol]
              ! mz_pj_20080806: MAX, since number density can be < 0
              zn = MAX(0.0_dp, pxtm1(_RI_X_ZN_(jp,jk,xsm(js)%idt_N(jm))) &
                + pxtte(_RI_X_ZN_(jp,jk,xsm(js)%idt_N(jm))) * time_step_len)
              ! convert [1/mol] to [1/m3]
              zn = zn * p / (R_gas * T)
              ! aerosol surface area [m2/m3]
              ! mz_pj_20080806: MAX, since wetradius can be < 0
              zsurface = 4._dp * pi * &
                   MAX(0.0_dp, wetradius(js)%ptr(_RI_XYZN_(jp,jrow,jk,jm)))**2 &
                * zn * EXP(2.0_dp * (LOG10(sigma(js)%ptr(jm)))**2 &
                / ((LOG10(EXP(1._dp)))**2))
              kh1 = kh1 + k_mt(M_N2O5,T,zsurface,gamma_N2O5)
              kh2 = kh2 + k_mt(M_HNO3,T,zsurface,gamma_HNO3)
              kh3 = kh3 + k_mt(MHg,T,zsurface,gamma_Hg)
              kh4 = kh4 + k_mt(MHg,T,zsurface,gamma_RGM)
            ENDDO mode_loop
            khet_Tr_3d(js,iht_N2O5)%ptr(_RI_XYZ__(jp,jrow,jk)) = kh1 ! N2O5 uptake
            khet_Tr_3d(js,iht_HNO3)%ptr(_RI_XYZ__(jp,jrow,jk)) = kh2 ! HNO3 uptake
            khet_Tr_3d(js,iht_Hg)%ptr(_RI_XYZ__(jp,jrow,jk))   = kh3 ! Hg (GEM) uptake
            khet_Tr_3d(js,iht_RGM)%ptr(_RI_XYZ__(jp,jrow,jk))  = kh4 ! Hg++ (RGM) uptake
          ENDDO submodel_loop

        ENDIF if_l_troposphere

        !---------------------------------------------------------------------

      ENDDO kproma_loop
    ENDDO level_loop

    ! mz_pj_20070305+
    IF (l_force_khet .OR. (l_stratosphere .AND. l_troposphere)) THEN
       submodel_loop2: DO js=0, nsm
          DO i = 1, IHT_MAX
             khet_Tr_3d(js,i)%ptr(_RI_XYZ__(1:kproma,jrow,1:nlev)) = &
                  khet_Tr_3d(js,i)%ptr(_RI_XYZ__(1:kproma,jrow,1:nlev)) * &
                  ( 1.0_dp - r_strat_flag(_RI_XYZ__(1:kproma,jrow,1:nlev)) )
          ENDDO
       END DO submodel_loop2

       ! IT IS ABSOLUTELY FORBIDDEN TO MODIFY OBJECTS FROM OTHER CHANNEL!!!
       ! WELL, ...
       DO i = 1, IHS_MAX
          khet_St_3d(i)%ptr(_RI_XYZ__(1:kproma,jrow,1:nlev)) = &
               khet_St_3d(i)%ptr(_RI_XYZ__(1:kproma,jrow,1:nlev)) * &
               r_strat_flag(_RI_XYZ__(1:kproma,jrow,1:nlev)) 
       ENDDO
    END IF
    ! mz_pj_20070305-

  END SUBROUTINE mecca_khet_physc_gp

  !***************************************************************************

!!#D attila +
  SUBROUTINE mecca_khet_physc_lg(temp,press)
#if defined(ECHAM5)
    USE messy_main_grid_def_mem_bi,ONLY: nproma, ngpblks, nlev 
    USE messy_main_timer,         ONLY: time_step_len
    USE messy_attila_tools_e5,    ONLY: gp2lg_e5
    USE messy_main_tracer_mem_bi, ONLY: pxtte_a=>qxtte_a, pxtm1_a=>qxtm1_a, &
                                        NCELL

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, EXP, NINT, LOG10, MAX
!
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_khet_physc_lg'
!
    REAL(dp), intent(in) :: temp(NCELL)  ! temperature in LG
    REAL(dp), intent(in) :: press(NCELL) ! pressure in LG
    INTEGER :: zphase, jp, jk, jm, js, i
    REAL(dp) :: T, p, zsurface, zr_s, zN_s0, zn
    REAL(dp) :: kh1, kh2, kh3, kh4
    REAL(dp), DIMENSION(30) :: kh
    
    REAL(dp), DIMENSION(:,:,:),ALLOCATABLE,TARGET :: wetradius_lg 
    REAL(dp), DIMENSION(:,:,:), POINTER :: wetradius_tmp_gp => NULL()
    REAL(dp), DIMENSION(:), POINTER :: wetradius_tmp_lg => NULL()

    !-------------------------------------------------------------------------
    ! BIT DIFFERENT FROM  GP: we transform first and then we define kh1 and kh2
    !-------------------------------------------------------------------------


    ! transformations:
    ! T is upgraded from LG
    ! gamma_N2O5 and gamma_HNO3 are constants from messy_mecca_khet.f90
    ! zsurface_3d is from climatology or submodels hend here transformed to LG
    ! M_N2O5 and M_N2O5 are constants from messy_mecca_khet.f90


    IF (l_climatology) THEN
        ALLOCATE(surface_1d(NCELL))
        CALL gp2lg_e5(surface_3d,surface_1d)
    ELSE
        ALLOCATE(wetradius_lg(nsm,NCELL,xsm(js)%nmodes))
        DO js=1, nsm
            DO jm = 1, xsm(js)%nmodes
              wetradius_tmp_lg => wetradius_lg(js,:,jm)
              wetradius_tmp_gp => wetradius(js)%ptr(:,:,jm,:)
              CALL gp2lg_e5(wetradius_tmp_gp,wetradius_tmp_lg)
            END DO
        END DO
    END IF

    IF ((l_force_khet).OR.(l_stratosphere)) THEN ! stratosphere
        DO i=1, IHS_MAX
          CALL gp2lg_e5(khet_St_3d(i)%ptr,khet_St_1d(i)%ptr)
        END DO
        ! mz_pj_20070305+
        ALLOCATE(r_strat_flag_1d(NCELL))
        CALL gp2lg_e5(r_strat_flag, r_strat_flag_1d)
        ! mz_pj_20070305-

    END IF


    !CELL LOOPS (as in GP)

    NCELL_loop: DO jp= 1, NCELL 

        T = temp(jp)
! mz_pj_20080930+
!        p = press(jp)/100.0_dp ! [hPa] !qqq why not Pa?
        p = press(jp) ! [Pa]
! mz_pj_20080930-

        !---------------------------------------------------------------------

        if_l_troposphere: IF (l_troposphere) THEN ! tropospheric aerosol

          IF (l_climatology) THEN

            zsurface = surface_1d(jp) ! [m2/m3]
            khet_Tr_1d(0,iht_N2O5)%ptr(jp) = & ! N2O5 uptake
              k_mt(M_N2O5,T,zsurface,gamma_N2O5)
            khet_Tr_1d(0,iht_HNO3)%ptr(jp) = & ! HNO3 uptake
              k_mt(M_HNO3,T,zsurface,gamma_HNO3)
            khet_Tr_1d(0,iht_Hg)%ptr(jp) = &   ! Hg (GEM) uptake
              k_mt(MHg,T,zsurface,gamma_Hg)
            khet_Tr_1d(0,iht_RGM)%ptr(jp) = &  ! Hg++ (RGM) uptake
              k_mt(MHg,T,zsurface,gamma_RGM)  ! assuming M_RGM = MHg
          ENDIF

          submodel_loop: DO js=1, nsm
            kh1 = 0._dp
            kh2 = 0._dp
            kh3 = 0._dp
            kh4 = 0._dp
            mode_loop: DO jm = 1, xsm(js)%nmodes
              ! zn = number density [1/mol]
              zn = pxtm1_a(jp,xsm(js)%idt_N(jm)) &
                + pxtte_a(jp,xsm(js)%idt_N(jm)) * time_step_len
              ! convert [1/mol] to [1/m3]
              zn = zn * p / (R_gas * T)
              ! aerosol surface area [m2/m3]
              zsurface = 4._dp * pi * wetradius_lg(js,jp,jm)**2 &
                * zn * EXP(2.0_dp * (LOG10(sigma(js)%ptr(jm)))**2 &
                / ((LOG10(EXP(1._dp)))**2))
              kh1 = kh1 + k_mt(M_N2O5,T,zsurface,gamma_N2O5)
              kh2 = kh2 + k_mt(M_HNO3,T,zsurface,gamma_HNO3)
              kh3 = kh3 + k_mt(MHg,T,zsurface,gamma_Hg)
              kh4 = kh4 + k_mt(MHg,T,zsurface,gamma_RGM)
            ENDDO mode_loop
            khet_Tr_1d(js,iht_N2O5)%ptr(jp) = kh1 ! N2O5 uptake
            khet_Tr_1d(js,iht_HNO3)%ptr(jp) = kh2 ! HNO3 uptake
            khet_Tr_1d(js,iht_Hg)%ptr(jp)   = kh3 ! Hg (GEM) uptake
            khet_Tr_1d(js,iht_RGM)%ptr(jp)  = kh4 ! Hg++ (RGM) uptake
          ENDDO submodel_loop

        ENDIF if_l_troposphere

        !---------------------------------------------------------------------
    ENDDO NCELL_loop

    ! mz_pj_20070305+
    IF (l_force_khet .OR. (l_stratosphere .AND. l_troposphere)) THEN
       submodel_loop2: DO js=0, nsm
          DO i = 1, IHT_MAX
             khet_Tr_1d(js,i)%ptr(:) = khet_Tr_1d(js,i)%ptr(:) * &
                  (1.0_dp - r_strat_flag_1d(:) )  
          ENDDO
       END DO submodel_loop2
! THIS IS NOT REQUIRED, BECAUSE IT IS TRANSFORMED FROM GP
! (AND THEREFORE ALREADY FLAGGED)       
!!$       DO i = 1, IHS_MAX
!!$          khet_St_1d(i)%ptr(:) = &
!!$               khet_St_1d(i)%ptr(:) * r_strat_flag_1d(:) 
!!$       ENDDO
    END IF

    IF (ASSOCIATED(r_strat_flag_1d)) THEN
       DEALLOCATE(r_strat_flag_1d)
       NULLIFY(r_strat_flag_1d)
    END IF
    ! mz_pj_20070305-

    IF (l_climatology) THEN
        DEALLOCATE(surface_1d)
    ELSE
        DEALLOCATE(wetradius_lg)
    ENDIF
#else
    REAL(dp), INTENT(IN), DIMENSION(:) :: temp   ! temperature in LG
    REAL(dp), INTENT(IN), DIMENSION(:) :: press  ! press in LG
#endif
  END SUBROUTINE mecca_khet_physc_lg
!!#D attila -

  !***************************************************************************
  SUBROUTINE mecca_khet_free_memory

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    INTEGER :: js

    IF (ASSOCIATED(khet_Tr_3d))      DEALLOCATE(khet_Tr_3d)
    IF (ASSOCIATED(khet_Tr_1d))      DEALLOCATE(khet_Tr_1d)
    DO js=1, nsm
      IF (ASSOCIATED(xsm(js)%idt_N)) DEALLOCATE(xsm(js)%idt_N)
      IF (ASSOCIATED(wetradius))     DEALLOCATE(wetradius)
      IF (ASSOCIATED(sigma))         DEALLOCATE(sigma)
    END DO

  END SUBROUTINE mecca_khet_free_memory

  !***************************************************************************

  SUBROUTINE mecca_khet_read_nml_cpl(status, iou)

    ! read coupling namelist

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, &
                                read_nml_close, &
                                strcrack

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    TYPE IO_KHET
      CHARACTER(LEN=STRLEN_MEDIUM) :: name  = ''
      CHARACTER(LEN=STRLEN_VLONG)  :: modes = ''
    END TYPE IO_KHET

    TYPE(IO_KHET), DIMENSION(NMAXSM) :: asm
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_khet_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: js, jm
    INTEGER                     :: asm_cpl
    CHARACTER(LEN=3), DIMENSION(:), POINTER :: modes => NULL()

    NAMELIST /CPL_KHET/ aerosurf_clim, asm, asm_cpl, strat_channel

    ! defaults:
    asm_cpl = 0
    
    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL_KHET', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_KHET, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_KHET', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    IF ((asm_cpl < 0) .OR. (asm_cpl > NMAXSM)) THEN
       CALL error_bi('asm_cpl out of range', substr)
    ELSE
       IF (asm_cpl > 0) THEN
          IF (TRIM(asm(asm_cpl)%name) == '') THEN
             CALL error_bi('asm(asm_cpl) name not defined', substr)
          ENDIF
          IF (TRIM(asm(asm_cpl)%modes) == '') THEN
             CALL error_bi('asm(asm_cpl) mode(s) not defined', substr)
          END IF
       END IF
    END IF
    
      nsm = 0
      xsm_cpl = 0 ! mz_pj_20070523
      DO js=1, NMAXSM
        IF (TRIM(asm(js)%name) == '') CYCLE
        nsm = nsm + 1
        xsm(nsm)%name = asm(js)%name
        IF (js==asm_cpl) xsm_cpl = nsm
        CALL info_bi('submodel: '//xsm(nsm)%name)
        CALL strcrack(asm(js)%modes, ',', modes, xsm(nsm)%nmodes)
        ALLOCATE(xsm(nsm)%mode(xsm(nsm)%nmodes))
        DO jm=1, xsm(nsm)%nmodes
          READ(modes(jm),*) xsm(nsm)%mode(jm)
          CALL info_bi('  mode: '//str(xsm(nsm)%mode(jm)))
        END DO
      END DO

      CALL info_bi('stratosphere channel: '        // strat_channel)
      CALL info_bi('Aerosol climatology channel: ' // aerosurf_clim(1))
      CALL info_bi('  surface:         '           // aerosurf_clim(2))
      IF (xsm_cpl > 0) THEN
        CALL info_bi('Coupling MECCA to: '           // xsm(xsm_cpl)%name)
      ELSE
        CALL info_bi('Coupling MECCA to: aerosol surface climatology')
      END IF

  END SUBROUTINE mecca_khet_read_nml_cpl

  !***************************************************************************

END MODULE messy_mecca_khet_si

!*****************************************************************************
