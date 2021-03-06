#include "messy_main_ppd_bi.inc"

!        1         2         3         4         5         6         7         8
!2345678901234567890123456789012345678901234567890123456789012345678901234567890

MODULE messy_made_si

  ! MODULE FOR MADE-ECHAM5 INTERFACE
  !
  ! MADE adopted to the structure of the Mainz Earth Submodel System (MESSy).
  !
  ! MADE was originally implemented in ECHAM by A. Lauer, DLR, 2001-2003
  ! Original MADE source code by I. Ackermann, et al., University Cologne,
  ! Germany, 1998.
  !
  ! P. Joeckel, DLR, 2011: modified for MESSy2 infrastructure
  ! C. Kaiser,  DLR, 2013: - further modifications for MESSy2 infrastructure
  !                        - `cleanup' with the help of MESSy's `gmake check'
  !                          and gmake `gmake messycheck'
  !                        - addition of `philfrac' calculation (hydrophilic
  !                          fraction of 3rd moment per mode)
  ! M. Righi, DLR, 2014: modes renamed for consistency with MADE3
  !                      (i -> km, j -> am, c -> cm)      

  ! MESSy
  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM

  IMPLICIT NONE
  SAVE

  PRIVATE

  INTRINSIC  TRIM, ABS, MAX, MIN, NULL

  ! define 5D pointer
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: wetradius     => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: dryradius     => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: densaer       => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: rhhist        => NULL()
! op_ck_20130115+
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: philfrac        => NULL()
! op_ck_20130115-

  ! define 4D pointer
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: wetrad_4d     => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: dryrad_4d     => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: densaer_4d    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: rhhist_4d     => NULL()
! op_ck_20130115+
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: philfrac_4d     => NULL()
! op_ck_20130115-

  ! seasalt emission calculations 
  REAL(dp), DIMENSION(:,:),       POINTER :: Mss_as => NULL()
  REAL(dp), DIMENSION(:,:),       POINTER :: Mss_cs => NULL()
  REAL(dp), DIMENSION(:,:),       POINTER :: Nss_as => NULL()
  REAL(dp), DIMENSION(:,:),       POINTER :: Nss_cs => NULL()

  ! dust emission mass flux
  REAL(dp), DIMENSION(:,:),       POINTER :: dust_emis => NULL()

  ! SOA emission mass flux 
  REAL(dp), DIMENSION(:,:,:),     POINTER :: Msoa => NULL()

  ! CPL - NAMELIST
!DEBUG+
  LOGICAL                       :: l_sootag          = .true.
  LOGICAL                       :: l_main            = .true.
  LOGICAL                       :: l_setrhhist       = .true.
  LOGICAL                       :: l_feedback        = .true.
!DEBUG-
  LOGICAL                       :: lcpl_gasphase     = .false.
  CHARACTER (LEN=STRLEN_MEDIUM) :: chemmodule        = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: H2SO4_gas(2)      = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: NH3_gas(2)        = ''
!  CHARACTER (LEN=STRLEN_MEDIUM) :: HCl_gas(2)        = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: HNO3_gas(2)       = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: ProdSO4(2)        = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SOA_channel       = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SOA_object        = ''
  LOGICAL                       :: l_calc_emis       = .false.
  LOGICAL                       :: l_tendency        = .false. ! um_ak_20110720
  ! ... sea salt emissions
  LOGICAL                       :: l_ss              = .false.
  CHARACTER (LEN=STRLEN_MEDIUM) :: SSemis_channel    = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SS_mass_as        = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SS_num_as         = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SS_mass_cs        = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SS_num_cs         = ''
  ! ... dust emissions
  LOGICAL                       :: l_dust            = .false.
  CHARACTER (LEN=STRLEN_MEDIUM) :: Duemis_channel    = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_dust         = ''

!DEBUG+
!!$    NAMELIST /CPL/ lcpl_gasphase, chemmodule, H2SO4_gas, NH3_gas, HNO3_gas, &
    NAMELIST /CPL/ l_sootag, l_main, l_setrhhist, l_feedback,               &
                   lcpl_gasphase, chemmodule, H2SO4_gas, NH3_gas, HNO3_gas, &
!DEBUG-
!                   HCl_gas,                                                 &
                   ProdSO4, SOA_channel, SOA_object,                        &
! um_ak_20110720   l_calc_emis, l_ss, SSemis_channel,                        &
                   l_calc_emis, l_tendency, l_ss, SSemis_channel,            &
                   SS_mass_as, SS_num_as, SS_mass_cs, SS_num_cs,            &
                   l_dust, Duemis_channel,emis_dust

  ! TRACER INDICES (ECHAM arrays)

  INTEGER :: ITRACSO4AI   = 0    ! Aitken mode sulfate aerosol conc.
  INTEGER :: ITRACSO4AJ   = 0    ! accumulation mode sulfate aerosol conc.
!  INTEGER :: ITRACSO4AC   = 0    ! coarse mode sulfate aerosol conc.
  INTEGER :: ITRACNU0     = 0    ! Aitken mode 0th moment (number)
  INTEGER :: ITRACAC0     = 0    ! accumulation mode 0th moment (number)
  INTEGER :: ITRACCORN    = 0    ! coarse mode 0th moment (number)
  INTEGER :: ITRACSULF    = 0    ! sulfuric acid vapor conc.
  INTEGER :: ITRACECI     = 0    ! Aitken mode BC, hydrophilic
  INTEGER :: ITRACECJ     = 0    ! accumulation mode BC, hydrophilic
!  INTEGER :: ITRACHCL     = 0    ! hydrochloric acid vapor concentration
  INTEGER :: ITRACNH3     = 0    ! ammonia gas concentration
  INTEGER :: ITRACNH4AI   = 0    ! Aitken mode ammonium aerosol conc.
  INTEGER :: ITRACNH4AJ   = 0    ! accumulation mode ammonium aerosol conc.
  INTEGER :: ITRACHNO3    = 0    ! nitric acid vapor conc.
  INTEGER :: ITRACNO3AI   = 0    ! Aitken mode nitrate conc.
  INTEGER :: ITRACNO3AJ   = 0    ! accumulation mode nitrate conc.
  INTEGER :: ITRACH2OAI   = 0    ! Aitken mode water conc.
  INTEGER :: ITRACH2OAJ   = 0    ! accumulation mode water conc.
  INTEGER :: ITRACH2OAC   = 0    ! coarse mode water conc.
  INTEGER :: ITRACORGPAI  = 0    ! Aitken mode POM, hydrophilic
  INTEGER :: ITRACORGPAJ  = 0    ! accumulation mode POM, hydrophilic
  INTEGER :: ITRACSEASI   = 0    ! Aitken mode sea salt
  INTEGER :: ITRACSEASJ   = 0    ! accumulation mode sea salt
  INTEGER :: ITRACSEASC   = 0    ! coarse mode sea salt
  INTEGER :: ITRACDUSTJ   = 0    ! accumulation mode mineral dust
  INTEGER :: ITRACDUSTC   = 0    ! coarse mode mineral dust
  INTEGER :: ITRACEC2I    = 0    ! Aitken mode BC, hydrophobic
  INTEGER :: ITRACEC2J    = 0    ! accumulation mode BC, hydrophobic
  INTEGER :: ITRACORGPA2I = 0    ! Aitken mode POM, hydrophobic
  INTEGER :: ITRACORGPA2J = 0    ! accumulation mode POM, hydrophobic

  INTEGER :: ITRACPRODSO4 = 0    ! prodution of H2SO4(g)

  INTEGER, PARAMETER :: num_split_facs = 8 ! number of splitting factors
                                           ! (for array dimensioning)

  REAL(dp), DIMENSION(:), POINTER :: sigma_str => NULL()

  ! SUBROUTINES
  PUBLIC  :: made_initialize    ! initialize MADE
  PUBLIC  :: made_new_tracer    ! define tracers
  PUBLIC  :: made_init_memory   ! allocate memory
  PUBLIC  :: made_init_coupling ! initialize coupling
  PUBLIC  :: made_vdiff         ! distribute online emissions
  PUBLIC  :: made_physc         ! calculate made-"chemistry"
  PUBLIC  :: made_free_memory   ! deallocate radius field

CONTAINS

!-----------------------------------------------------------------------------
  SUBROUTINE made_initialize

    ! ECHAM5
    USE messy_main_tools,      ONLY: find_next_free_unit
    USE messy_main_mpi_bi,     ONLY: p_bcast, p_parallel_io, p_io
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_made,            ONLY: sginin, sginia, sginic,     &
                                   bctime, pomtime

    ! SUBROUTINES
    USE messy_made,          ONLY: made_initialize_core, made_read_nml

    IMPLICIT NONE
  
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'made_initialize'
    INTEGER                      :: iou    ! I/O unit
    INTEGER                      :: status ! error status

    !--- 1.) Read MADE namelist and control variables:

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL made_read_nml(status, iou)
       IF (status /= 0)  CALL error_bi('error in made_read_nml', substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast (sginin,        p_io)
    CALL p_bcast (sginia,        p_io)
    CALL p_bcast (sginic,        p_io)
    CALL p_bcast (bctime,        p_io)
    CALL p_bcast (pomtime,       p_io)

    !--- 2.)  Read MADE CPL namelist:
    
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL made_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ', substr)
    END IF

    ! BROADCAST RESULTS

    ! ... chemistry
    CALL p_bcast(lcpl_gasphase,    p_io)
    CALL p_bcast(chemmodule,       p_io)
    CALL p_bcast(H2SO4_gas,        p_io)
    CALL p_bcast(NH3_gas,          p_io)
!    CALL p_bcast(HCl_gas,          p_io)
    CALL p_bcast(HNO3_gas,         p_io)
    CALL p_bcast(ProdSO4,          p_io)
    ! ... SOA "emissions"
    CALL p_bcast (SOA_channel,      p_io)
    CALL p_bcast (SOA_object,     p_io)
    ! ... emissions
    CALL p_bcast (l_calc_emis,     p_io)
    CALL p_bcast (l_tendency,      p_io) ! um_ak_20110720
    CALL p_bcast (l_ss,            p_io)
    CALL p_bcast (SSemis_channel,   p_io)
    CALL p_bcast (SS_mass_as,      p_io)
    CALL p_bcast (SS_num_as,       p_io)
    CALL p_bcast (SS_mass_cs,      p_io)
    CALL p_bcast (SS_num_cs,       p_io)
    CALL p_bcast (l_dust,          p_io)
    CALL p_bcast (Duemis_channel,   p_io)
    CALL p_bcast (emis_dust,       p_io)

    !--- 3.) Initialize MADE (MUST be call AFTER 'made_read_nml'!!!):

    CALL made_initialize_core

  END SUBROUTINE made_initialize

!-----------------------------------------------------------------------------

  SUBROUTINE made_new_tracer

    ! ECHAM5
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer,          ONLY: new_tracer, set_tracer, get_tracer,  &
                                          AIR, AEROSOL, OFF, ON,               &
                                          AMOUNTFRACTION, NUMBERDENSITY,       &
                                          I_ADVECT, I_AEROSOL_MODE,            &
! op_mr_20140415 +
                                          I_AEROSOL_HETICE,                    &
! op_mr_20140415-
                                          I_AEROSOL_SOL, I_CONVECT, I_DRYDEP,  &
                                          I_SCAV, I_SEDI, I_VDIFF,             &
                                          R_AEROSOL_DENSITY, R_DRYREAC_SF,     &
                                          R_PSS  , R_MOLARMASS, S_AEROSOL_MODEL
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi,    &
                                          info_bi
    USE messy_made,                 ONLY: mwso4, mwnh4, mwno3, mwh2o, mwec,    &
                                          mworg, mwh2so4, mwnh3, mwhno3,       &
                                          mwseas, mwsoil, & !mwhcl,            &
                                          rhoso4, rhonh4, rhono3, rhoh2o,      &
                                          rhoanth, rhoorg, rhoseas, rhosoil,   &
                                          modstr, akn, acc, cor

    IMPLICIT NONE
    INTRINSIC TRIM

    !--- Module variables:
    !
    !    Tracer indices:
    !
    !    Legend:
    !
    !      i:  Aitken mode
    !      j:  accumulation mode
    !      c:  coarse mode
    !
    ! Parameters:
    ! -----------
    ! User defined flags: density   density                    [kg m-3]
    !                     osm       osmotic coefficient        [???]
    !                     nion      number of ions the tracer 
    !                               dissolves into             [1]

    ! Default settings for new tracers (set in smcl/messy_main_tracer_va.f90):
    ! ------------------------------------------------------------------------
    ! DEFAULT_CASK_I(I_ADVECT)          = ON     ! ADVECTION
    ! DEFAULT_CASK_I(I_CONVECT)         = ON     ! CONVECTION
    ! DEFAULT_CASK_I(I_VDIFF)           = ON     ! VERTICAL DIFFUSION
    ! DEFAULT_CASK_I(I_DRYDEP)          = OFF    ! DRY DEPOSITION
    ! DEFAULT_CASK_I(I_SEDI)            = OFF    ! SEDIMENTATION
    ! DEFAULT_CASK_I(I_SCAV)            = OFF    ! SCAVENGING
    ! DEFAULT_CASK_I(I_MIX)             = ON     ! TURBULENT MIXING
    ! DEFAULT_CASK_I(I_FORCE_COL)       = OFF    ! FORCING IN COLUMN MODE
    ! DEFAULT_CASK_I(I_INTEGRATE)       = ON     ! TIME INTEGRATION
    ! DEFAULT_CASK_I(I_TIMEFILTER)      = ON     ! TIME FILTER
    ! DEFAULT_CASK_I(I_FORCE_INIT)      = OFF    ! FORCE INIT AFTER RESTART
    ! DEFAULT_CASK_I(I_AEROSOL_METHOD)  = MODAL  ! MODAL OR BIN
    ! DEFAULT_CASK_I(I_AEROSOL_MODE)    = 0      ! MODE OR BIN NUMBER
    ! DEFAULT_CASK_I(I_AEROSOL_SOL)     = ON     ! SOLUBLE ON/OFF
    ! DEFAULT_CASK_I(I_AEROSOL_HETICE)  = OFF    ! HIGHER ICE SCAV. COEFF.
    ! DEFAULT_CASK_I(I_HDIFF)           = OFF    ! HORIZONTAL DIFFUSION
    ! DEFAULT_CASK_I(I_RELAX)           = ON     ! BOUNDARY DATA AVAILABLE
                                                 ! i.e. relaxation possible
    ! DEFAULT_CASK_I(I_MMD_INIT)        = OFF
    ! DEFAULT_CASK_I(I_TAG_REG_IDT)     = 0      ! id of associated regular
                                                 ! species
    ! DEFAULT_CASK_I(I_TAG_SPECIFIC)    = 0      ! flag for special
    !                                            ! treatment
    ! DEFAULT_CASK_R(R_MOLARMASS)       = 0.0_dp ! MOLAR MASS
    ! DEFAULT_CASK_R(R_PSS  )           = 0.0_dp ! HENRY'S LAW CONSTANT
    ! DEFAULT_CASK_R(R_DRYREAC_SF)      = 0.0_dp ! 
    ! DEFAULT_CASK_R(R_VINI)            = 0.0_dp ! 
    ! DEFAULT_CASK_R(R_AEROSOL_DENSITY) = 0.0_dp ! DENSITY OF BULK SPECIES
    ! DEFAULT_CASK_S(S_AEROSOL_MODEL)   = ''     ! AEROSOL MODEL NAME


    ! LOCAL
    INTEGER                     :: status, idt
    CHARACTER(LEN=*), PARAMETER :: substr = 'made_new_tracer'
    CHARACTER(LEN=10)           :: tracername
    CHARACTER(LEN=10)           :: othermod = ''

    CALL start_message_bi(modstr,'REQUEST MADE TRACER', substr)

    IF (.not. lcpl_gasphase) THEN

       tracername = 'H2SO4dummy'
       CALL new_tracer(status, GPTRSTR, tracername, modstr,                 &
            idx=idt, longname='dummy H2SO4 mixing ratio - gas phase', &
            unit='mol/mol', medium=AIR)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwh2so4)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_PSS  , r=8.7e11_dp)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       itracsulf = idt

!       ! Henry number and dryreac_sf from MECCA1
!
!       tracername = 'HCldummy'
!       CALL new_tracer(status, GPTRSTR, tracername, modstr,        &
!            idx=idt, longname='dummy HCl mixing ratio - gas phase', &
!            unit='mol/mol', medium=AIR)
!       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!       CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwhcl)
!       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!       CALL set_tracer(status, GPTRSTR, idt, R_PSS  , r=1.e14_dp)
!       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!       CALL set_tracer(status, GPTRSTR, idt, R_DRYREAC_SF, r=1.0_dp)
!       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!       itrachcl = idt

       ! Henry number and dryreac_sf from MECCA1

       tracername = 'HNO3dummy'
       CALL new_tracer(status, GPTRSTR, tracername, modstr,       &
            idx=idt, longname='dummy HNO3 mixing ratio - gas phase', &
            unit='mol/mol', medium=AIR)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwhno3)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_PSS  , r=1.e4_dp)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_DRYREAC_SF, r=1.0_dp)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       itrachno3 = idt

       ! Henry number and dryreac_sf from MECCA1

! op_ck_20120615 FIXME: Why are drydep, wetdep, and scav switched on here?
       tracername = 'NH3dummy'
       CALL new_tracer(status, GPTRSTR, tracername, modstr,        &
            idx=idt, longname='dummy NH3 mixing ratio - gas phase', &
            unit='mol/mol', medium=AIR)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwnh3)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_PSS  , r=58.0_dp)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       itracnh3 = idt

    ELSE

       ! if gas tracer names are not specified, create dummy tracers for proper
       ! operation despite missing species

       IF (HNO3_gas(1) == '' .AND. HNO3_gas(2) == '') THEN
          tracername = 'HNO3'
          CALL info_bi('No gas-phase '//TRIM(tracername)//' tracer specified,',&
               substr)
          CALL info_bi('looking for existing tracer...', substr)
          ! Check if tracer has already been defined
          CALL get_tracer(status, GPTRSTR, TRIM(tracername), idx=idt, &
               submodel=othermod)
          IF (status == 0) THEN   ! use existing tracer
             CALL info_bi('...using '//TRIM(tracername)//' tracer defined by '&
                  //TRIM(othermod), substr)
          ELSE   ! define new tracer
             CALL info_bi('...not found. Creating '//TRIM(tracername)&
                  //' tracer.', substr)
             ! Henry number and dryreac_sf from MECCA1
             CALL new_tracer(status, GPTRSTR, TRIM(tracername), modstr,      &
                  longname=TRIM(tracername)//' mixing ratio - gas phase', &
                  idx=idt, unit='mol/mol', medium=AIR)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
             CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwhno3)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
             CALL set_tracer(status, GPTRSTR, idt, R_PSS  , r=1.e4_dp)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
             CALL set_tracer(status, GPTRSTR, idt, R_DRYREAC_SF, r=1.0_dp)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
          END IF
          itrachno3 = idt
       END IF

! op_ck_20120615 FIXME: Why are drydep, wetdep, and scav switched on here (for
!                       the other gases they are not)? And why is sedi not
!                       switched on? And what about dryreac_sf?
       IF (NH3_gas(1) == '' .AND. NH3_gas(2) == '') THEN
          tracername = 'NH3'
          CALL info_bi('No gas-phase '//TRIM(tracername)//' tracer specified,',&
               substr)
          CALL info_bi('looking for existing tracer...', substr)
          ! Check if tracer has already been defined
          CALL get_tracer(status, GPTRSTR, TRIM(tracername), idx=idt, &
               submodel=othermod)
          IF (status == 0) THEN   ! use existing tracer
             CALL info_bi('...using '//TRIM(tracername)//' tracer defined by '&
                  //TRIM(othermod), substr)
          ELSE   ! define new tracer
             CALL info_bi('...not found. Creating '//TRIM(tracername)&
                  //' tracer.', substr)
             ! Henry number and dryreac_sf from MECCA1
             CALL new_tracer(status, GPTRSTR, TRIM(tracername), modstr, &
                  longname=TRIM(tracername)//' mixing ratio - gas phase', &
                  idx=idt, unit='mol/mol', medium=AIR)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
             CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
             CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
             CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwnh3)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
             CALL set_tracer(status, GPTRSTR, idt, R_PSS  , r=58.0_dp)
             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
          END IF
          itracnh3 = idt
       END IF

!       IF (HCl_gas(1) == '' .AND. HCl_gas(2) == '') THEN
!          ! Henry number and dryreac_sf from MECCA1
!          tracername = 'HCl'
!          CALL info_bi('No gas-phase '//TRIM(tracername)//' tracer specified,',&
!               substr)
!          CALL info_bi('looking for existing tracer...', substr)
!          ! Check if tracer has already been defined
!          CALL get_tracer(status, GPTRSTR, TRIM(tracername), idx=idt, &
!               submodel=othermod)
!          IF (status == 0) THEN   ! use existing tracer
!             CALL info_bi('...using '//TRIM(tracername)//' tracer defined by '&
!                  //TRIM(othermod), substr)
!          ELSE   ! define new tracer
!             CALL info_bi('...not found. Creating '//TRIM(tracername)&
!                  //' tracer', substr)
!             ! Henry number and dryreac_sf from MECCA1
!             CALL new_tracer(status, GPTRSTR, TRIM(tracername), modstr,  &
!                  longname=TRIM(tracername)//' mixing ratio - gas phase', &
!                  idx=idt, unit='mol/mol', medium=AIR)
!             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!             CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwhcl)
!             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!             CALL set_tracer(status, GPTRSTR, idt, R_PSS  , r=1.e14_dp)
!             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!             CALL set_tracer(status, GPTRSTR, idt, R_DRYREAC_SF, r=1.0_dp)
!             CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!          END IF
!          itrachcl = idt
!       END IF

    END IF ! (.not. lcpl_gasphase)

    !--- 2) Allocate aerosol masses according to the modes to obtain 
    !       succeding tracer identifiers:

    ! I_SCAV            = OFF ---> don't calculate scavenging at all
    !                   = ON  ---> calculate impact scavenging
    ! I_AEROSOL_SOL     = ON  ---> calculate nucleation scavenging
    !                              (only if I_SCAV=ON)

    !--- Mode 1 - Aitken mode:
    
    tracername = 'SO4'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol sulfate mass mixing ratio - Aitken mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwso4)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoso4)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracso4ai = idt
    
    tracername = 'NH4'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol ammonium mass mixing ratio - Aitken mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwnh4)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhonh4)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracnh4ai = idt

    tracername = 'NO3'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol nitrate mass mixing ratio - Aitken mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwno3)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhono3)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracno3ai = idt

! op_ck_20150322+
!!$    tracername = 'BC'
    tracername = 'BCphil'
! op_ck_20150322-
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol BC (hydrophilic) mass mixing ratio - Aitken mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415+
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_HETICE, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415-
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwec)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoanth)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itraceci = idt

! op_ck_20150322+
!!$    tracername = 'BC2'
   tracername = 'BCphob'
! op_ck_20150322-
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol BC (hydrophobic) mass mixing ratio - Aitken mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_SOL, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415+
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_HETICE, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415-
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwec)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoanth)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracec2i = idt

! op_ck_20150322+
!!$    tracername = 'POM'
    tracername = 'POMphil'
! op_ck_20150322-
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol POM (hydrophilic) mass mixing ratio - Aitken mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mworg)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoorg)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracorgpai = idt

! op_ck_20150322+
!!$    tracername = 'POM2'
    tracername = 'POMphob'
! op_ck_20150322-
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol POM (hydrophobic) mass mixing ratio - Aitken mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_SOL, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mworg)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoorg)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracorgpa2i = idt

    IF (.not.l_ss) THEN

       tracername = 'SS'
       CALL new_tracer(status, GPTRSTR, tracername, modstr, &
            idx=idt, subname='km', longname = &
            'Aerosol sea salt mass mixing ratio - Aitken mode', &
            unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwseas)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoseas)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       itracseasi = idt

    ELSE

       ! dummy tracer
       tracername = 'SS'
       CALL new_tracer(status, GPTRSTR, tracername, modstr, &
            idx=idt, subname='km', longname = &
            'dummy aerosol sea salt mass mixing ratio - Aitken mode', &
            unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_ADVECT, i=OFF)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_CONVECT, i=OFF)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_VDIFF, i=OFF)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwseas)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoseas)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
       CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
       itracseasi = idt

    END IF

    tracername = 'H2O'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol water mass mixing ratio - Aitken mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_ADVECT, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwh2o)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoh2o)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itrach2oai = idt

    !--- Mode 2 - accumulation mode:

    tracername = 'SO4'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol sulfate mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwso4)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoso4)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracso4aj = idt
    
    tracername = 'NH4'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol ammonium mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwnh4)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhonh4)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracnh4aj = idt

    tracername = 'NO3'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol nitrate mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwno3)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhono3)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracno3aj = idt

! op_ck_20150322+
!!$    tracername = 'BC'
    tracername = 'BCphil'
! op_ck_20150322-
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol BC (hydrophilic) mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415+
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_HETICE, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415-
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwec)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoanth)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracecj = idt

! op_ck_20150322+
!!$    tracername = 'BC2'
    tracername = 'BCphob'
! op_ck_20150322-
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol BC (hydrophobic) mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_SOL, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415+
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_HETICE, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415-
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwec)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoanth)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracec2j = idt

! op_ck_20150322+
!!$    tracername = 'POM'
    tracername = 'POMphil'
! op_ck_20150322-
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol POM (hydrophilic) mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mworg)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoorg)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracorgpaj = idt

! op_ck_20150322+
!!$    tracername = 'POM2'
    tracername = 'POMphob'
! op_ck_20150322-
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol POM (hydrophobic) mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_SOL, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mworg)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoorg)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracorgpa2j = idt

    tracername = 'SS'
    CALL new_tracer(status, GPTRSTR, tracername   ,modstr  ,&
         idx=idt, subname='am', longname = &
         'Aerosol sea salt mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwseas)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoseas)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracseasj = idt

    tracername = 'DU'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol mineral dust mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_SOL, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415+
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_HETICE, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415-
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwsoil)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhosoil)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracdustj = idt

    tracername = 'H2O'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol water mass mixing ratio - accumulation mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_ADVECT, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwh2o)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoh2o)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itrach2oaj = idt

    !--- Mode 3 - Coarse mode:

    tracername = 'SS'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='cm', longname = &
         'Aerosol sea salt mass mixing ratio - coarse mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=cor)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwseas)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoseas)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracseasc = idt

    tracername = 'DU'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='cm', longname = &
         'Aerosol mineral dust mass mixing ratio - coarse mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=cor)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_SOL, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415+
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_HETICE, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
! op_mr_20140415-
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwsoil)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhosoil)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracdustc = idt

!    tracername = 'SO4'
!    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
!         idx=idt, subname='cm', longname = &
!         'Aerosol sulfate mass mixing ratio - coarse mode', &
!         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=cor)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwso4)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoso4)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    itracso4ac = idt
!
!    tracername = 'NH4'
!    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
!         idx=idt, subname='cm', longname = &
!         'Aerosol ammonium mass mixing ratio - coarse mode', &
!         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=cor)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwnh4)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhonh4)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    itracnh4ac = idt
!
!    tracername = 'NO3'
!    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
!         idx=idt, subname='cm', longname = &
!         'Aerosol nitrate mass mixing ratio - coarse mode', &
!         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=cor)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwno3)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhono3)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
!    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
!    itracno3ac = idt

    tracername = 'H2O'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='cm', longname = &
         'Aerosol water mass mixing ratio - coarse mode', &
         unit='mol/mol', medium=AEROSOL, quantity=AMOUNTFRACTION)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_ADVECT, i=OFF)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=cor)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=mwh2o)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=rhoh2o)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itrach2oac = idt

    !--- 3) Aerosol Numbers:
    ! NUM renamed to N, because this name is hardcoded in EMDEP 
    ! otherwise no emissions
    ! unit has to be written 1/mol not 1 cm-3 alos caused by EMDEP
    ! choose unit as 1/kg(air) or 1/mol(air)

    tracername = 'N'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='km', longname = &
         'Aerosol number mixing ratio - Aitken mode', &
         unit='1/mol', medium=AEROSOL, quantity=NUMBERDENSITY)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=akn)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=1.0_dp)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=1.0_dp)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracnu0 = idt

    tracername = 'N'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='am', longname = &
         'Aerosol number mixing ratio - accumulation mode', &
         unit='1/mol', medium=AEROSOL, quantity=NUMBERDENSITY)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=acc)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=1.0_dp)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=1.0_dp)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itracac0 = idt

    tracername = 'N'
    CALL new_tracer(status, GPTRSTR, tracername, modstr, &
         idx=idt, subname='cm', longname = &
         'Aerosol number mixing ratio - coarse mode', &
         unit='1/mol', medium=AEROSOL, quantity=NUMBERDENSITY)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=cor)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=1.0_dp)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=1.0_dp)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
    CALL tracer_halt(TRIM(substr)//' ('//TRIM(tracername)//')', status)
    itraccorn = idt

    CALL end_message_bi(modstr,'REQUEST MADE TRACER', substr)

  END SUBROUTINE made_new_tracer

!-----------------------------------------------------------------------------

  SUBROUTINE made_init_memory

    ! ECHAM5
    USE messy_main_blather_bi,         ONLY: start_message_bi, end_message_bi
    USE messy_main_grid_def_mem_bi,    ONLY: nproma, nlev, ngpblks
    USE messy_main_channel,            ONLY: new_channel, new_channel_object &
                                           , new_attribute
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO &
                                           , set_representation_decomp &
                                           , IRANK, PIOTYPE_COL        &
                                           , repr_def_axes
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: GP_3D_MID, DC_GP                &
                                           , DIMID_LON, DIMID_LAT, DIMID_LEV &
                                           , DC_BC &
                                           , gp_nseg, gp_start, gp_cnt &
                                           , gp_meml, gp_memu

    ! MADE
    USE messy_made, ONLY: modstr, nmod, sigma, akn, acc, cor

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='made_init_memory'
    CHARACTER(len=3), PARAMETER :: strmod = '_gp'
    INTEGER                               :: status
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: mem => NULL()
    INTEGER                               :: DIMID_NMODE
    INTEGER                               :: REPR_MADE_4D_NMOD
    INTEGER                               :: REPR_MADE_1D

    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr,'INITIALIZE MEMORY', substr)

! um_ak_20110720+
!!$    ALLOCATE(wetradius(nproma,nlev,nmod,ngpblks,1))
!!$    ALLOCATE(dryradius(nproma,nlev,nmod,ngpblks,1))
!!$    ALLOCATE(densaer(nproma,nlev,nmod,ngpblks,1))
!!$    ALLOCATE(rhhist(nproma,nlev,nmod,ngpblks,1))
    ALLOCATE(wetradius(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    ALLOCATE(dryradius(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    ALLOCATE(densaer(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    ALLOCATE(rhhist(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
! op_ck_20130115+
    ALLOCATE(philfrac(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
! op_ck_20130115-
! op_ck_20131010+
    wetradius(:,:,:,:,1) = 0._dp
    dryradius(:,:,:,:,1) = 0._dp
    densaer(:,:,:,:,1)   = 0._dp
    rhhist(:,:,:,:,1)    = 0._dp
    philfrac(:,:,:,:,1)  = 0._dp
! op_ck_20131010-
! um_ak_20110720-

    !--- 1) Construct the MADE channel: --------------------------------------!
    CALL new_dimension(status, DIMID_NMODE, 'MADE_NMODE', nmod)
    CALL channel_halt(substr, status)

    ! NEW REPRESENTATIONS
    CALL new_representation(status, REPR_MADE_4D_NMOD, &
         'REPR_MADE_4D_NMOD'    &
         , rank = 4, link = 'xxxx', dctype = DC_GP               &
         , dimension_ids = (/ &
            _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_NMODE) /) &
         , ldimlen       = (/ &
            _RI_XYZN_(nproma, ngpblks, AUTO, AUTO) /)   &
         , output_order  = (/ _IN_XYZN_, _IX_XYZN_      &     ! E: 3,1,4,2
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
    
    cnt(:,_IN_XYZN_)  = nmod
    memu(:,_IN_XYZN_) = nmod
    
    CALL set_representation_decomp(status, REPR_MADE_4D_NMOD &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    CALL new_representation(status, REPR_MADE_1D,     &
         'REPR_MADE_1D'                               &
         , rank = 1, link = 'x---', dctype = DC_BC  &
         , dimension_ids = (/ DIMID_NMODE /)        &
         , ldimlen       = (/ AUTO /)               &
         , axis = 'N---'                            &
         )
    CALL channel_halt(substr, status)
    ! mz_pj_20061112+
    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))
    
    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1
    
    start(:,1) = 1
    cnt(:,1)   = nmod
    meml(:,1)  = 1
    memu(:,1)  = nmod
    
    CALL set_representation_decomp(status, REPR_MADE_1D &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    CALL new_channel(status, modstr//strmod, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    !--- 2) Non accumulated diagnostic fields: -------------------------------

!!$    CALL default_stream_setting (made_stream, laccu      = .FALSE., &
!!$                                              lpost      = .FALSE., &
!!$                                              lrerun     = .TRUE.,  &
!!$                                              contnorest = .TRUE.,  &
!!$                                              lav        = .TRUE.,  &
!!$                                              leveltype  = HYBRID)

    !--- 2.1) Add channel objects:

    ! WET RADII
    ! create 4D channel
!!$    p4 => wetradius(:,:,:,:,1)
!!$    CALL add_stream_element (made_stream, 'wetradius', wetrad_4d,       &
!!$                             ktrac=nmod, units='m',lrerun=.true.,       &
!!$                             longname='wet particle radius',            &
!!$                             lpost=.false., lav=.false., p4=p4)
    mem => wetradius(:,:,:,:,1)
    CALL new_channel_object(status, modstr//strmod, 'wetradius' &
         , p4 = wetrad_4d, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'wetradius', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'wetradius' &
         , 'long_name', c='wet particle radius')
    CALL channel_halt(substr, status)

!!$    p4 => wetradius(:,:,akn,:,:)
!!$    CALL add_stream_element (made_stream, 'wetrad_km', p3, p4=p4,        &
!!$                      longname='wet particle radius (Aitken mode)',     &
!!$                      units='m', lrerun=.false., lpost=.true.,          &
!!$                      lav=.false.)    
    mem => wetradius(_RI_XYZN_(:,:,:,akn),:)
    CALL new_channel_object(status, modstr//strmod, 'wetrad_km', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'wetrad_km' &
         , 'long_name', c='wet particle radius (Aitken mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'wetrad_km', 'units', c='m')
    CALL channel_halt(substr, status)

!!$    p4 => wetradius(:,:,acc,:,:)
!!$    CALL add_stream_element (made_stream, 'wetrad_am', p3, p4=p4,        &
!!$                      longname='wet particle radius (acc. mode)',       &
!!$                      units='m', lrerun=.false., lpost=.true.,          &
!!$                      lav=.false.)
    mem => wetradius(_RI_XYZN_(:,:,:,acc),:)
    CALL new_channel_object(status, modstr//strmod, 'wetrad_am', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'wetrad_am' &
         , 'long_name', c='wet particle radius (acc. mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'wetrad_am', 'units', c='m')
    CALL channel_halt(substr, status)

!!$    p4 => wetradius(:,:,cor,:,:)
!!$    CALL add_stream_element (made_stream, 'wetrad_cm', p3, p4=p4,        &
!!$                      longname='wet particle radius (coarse mode)',     &
!!$                      units='m', lrerun=.false., lpost=.true.,          &
!!$                      lav=.false.)
    mem => wetradius(_RI_XYZN_(:,:,:,cor),:)
    CALL new_channel_object(status, modstr//strmod, 'wetrad_cm', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'wetrad_cm' &
         , 'long_name', c='wet particle radius (coarse mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'wetrad_cm' &
         , 'units', c='m')
    CALL channel_halt(substr, status)

    ! DRY RADII
    ! create 4D channel
!!$    p4 => dryradius(:,:,:,:,1)
!!$    CALL add_stream_element(made_stream, 'dryradius', dryrad_4d,        &
!!$                            ktrac=nmod, units='m', lrerun=.true.,       &
!!$                            longname='dry particle radius',             &
!!$                            lpost=.false., lav=.false., p4=p4)
    mem => dryradius(:,:,:,:,1)
    CALL new_channel_object(status, modstr//strmod, 'dryradius' &
         , p4 = dryrad_4d, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'dryradius', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'dryradius' &
         , 'long_name', c='dry particle radius')
    CALL channel_halt(substr, status)

    ! create 3d channels
!!$    p4 => dryradius(:,:,akn,:,:)
!!$    CALL add_stream_element (made_stream, 'dryrad_km', p3, p4=p4,        &
!!$                      longname='dry particle radius (Aitken mode)',     &
!!$                      units='m', lrerun=.false., lpost=.true.,          &
!!$                      lav=.false.)
    mem => dryradius(_RI_XYZN_(:,:,:,akn),:)
    CALL new_channel_object(status, modstr//strmod, 'dryrad_km', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'dryrad_km' &
         , 'long_name', c='dry particle radius (Aitken mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'dryrad_km', 'units', c='m')
    CALL channel_halt(substr, status)

!!$    p4 => dryradius(:,:,acc,:,:)
!!$    CALL add_stream_element (made_stream, 'dryrad_am', p3, p4=p4,        &
!!$                      longname='dry particle radius (acc. mode)',       &
!!$                      units='m', lrerun=.false., lpost=.true.,          &
!!$                      lav=.false.)
    mem => dryradius(_RI_XYZN_(:,:,:,acc),:)
    CALL new_channel_object(status, modstr//strmod, 'dryrad_am', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'dryrad_am' &
         , 'long_name', c='dry particle radius (acc. mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'dryrad_am', 'units', c='m')
    CALL channel_halt(substr, status)

!!$    p4 => dryradius(:,:,cor,:,:)
!!$    CALL add_stream_element (made_stream, 'dryrad_cm', p3, p4=p4,        &
!!$                      longname='dry particle radius (coarse mode)',     &
!!$                      units='m', lrerun=.false., lpost=.true.,          &
!!$                      lav=.false.)
    mem => dryradius(_RI_XYZN_(:,:,:,cor),:)
    CALL new_channel_object(status, modstr//strmod, 'dryrad_cm', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'dryrad_cm' &
         , 'long_name', c='dry particle radius (coarse mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'dryrad_cm', 'units', c='m')
    CALL channel_halt(substr, status)

    ! DENSITY
    ! create 4d channel_object
!!$    p4 => densaer(:,:,:,:,1)    
!!$    CALL add_stream_element(made_stream, 'densaer', densaer_4d,         &
!!$                            ktrac=nmod, p4=p4, units='kg m-3',          &
!!$                            longname='particle density',                &
!!$                            lrerun=.true., lpost=.false., lav=.false.)
    mem => densaer(:,:,:,:,1)    
    CALL new_channel_object(status, modstr//strmod, 'densaer' &
         , p4 = densaer_4d, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'densaer', 'units', c='kg m-3')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'densaer' &
         , 'long_name', c='particle density')
    CALL channel_halt(substr, status)

!!$    p4 => densaer(:,:,akn,:,:)
!!$    CALL add_stream_element (made_stream, 'densaer_km', p3, p4=p4,       &
!!$                         longname='particle density (Aitken mode)',     &
!!$                         units='kg m-3', lpost=.true.,                  &
!!$                         lrerun=.false., lav=.false.)
    mem => densaer(_RI_XYZN_(:,:,:,akn),:)
    CALL new_channel_object(status, modstr//strmod, 'densaer_km', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'densaer_km' &
         , 'long_name', c='particle density (Aitken mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'densaer_km', 'units', c='kg m-3')
    CALL channel_halt(substr, status)

!!$    p4 => densaer(:,:,acc,:,:)
!!$    CALL add_stream_element (made_stream, 'densaer_am', p3, p4=p4,       &
!!$                         longname='particle density (acc. mode)',       &
!!$                         units='kg m-3', lpost=.true.,                  &
!!$                         lrerun=.false., lav=.false.)
    mem => densaer(_RI_XYZN_(:,:,:,acc),:)
    CALL new_channel_object(status, modstr//strmod, 'densaer_am', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'densaer_am' &
         , 'long_name', c='particle density (acc. mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'densaer_am', 'units', c='kg m-3')
    CALL channel_halt(substr, status)

!!$    p4 => densaer(:,:,cor,:,:)
!!$    CALL add_stream_element (made_stream, 'densaer_cm', p3, p4=p4,       &
!!$                         longname='particle density (coarse mode)',     &
!!$                         units='kg m-3', lpost=.true.,                  &
!!$                         lrerun=.false., lav=.false.)
    mem => densaer(_RI_XYZN_(:,:,:,cor),:)
    CALL new_channel_object(status, modstr//strmod, 'densaer_cm', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'densaer_cm' &
         , 'long_name', c='particle density (coarse mode)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'densaer_cm', 'units', c='kg m-3')
    CALL channel_halt(substr, status)

    ! relative humidity history
    ! create 4D channel
!!$    p4 => rhhist(:,:,:,:,1)
!!$    CALL add_stream_element (made_stream, 'rhhist', rhhist_4d,          &
!!$                             ktrac = nmod, units='',lrerun=.true.,      &
!!$                             longname='history of rel. humidity',       &
!!$                             lpost=.false., lav=.false., p4=p4)
    mem => rhhist(:,:,:,:,1)
    CALL new_channel_object(status, modstr//strmod, 'rhhist' &
         , p4 = rhhist_4d, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'rhhist', 'units', c='')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'rhhist' &
         , 'long_name', c='history of rel. humidity')
    CALL channel_halt(substr, status)

!!$    p4 => rhhist(:,:,akn,:,:)
!!$    CALL add_stream_element (made_stream, 'rhhist_km', p3, p4=p4,        &
!!$                             units='', lrerun=.true., lpost=.false.,    &
!!$                             lav=.false.)
    mem => rhhist(_RI_XYZN_(:,:,:,akn),:)
    CALL new_channel_object(status, modstr//strmod, 'rhhist_km' &
         , mem = mem, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'rhhist_km', 'units', c='')
    CALL channel_halt(substr, status)

!!$    p4 => rhhist(:,:,acc,:,:)
!!$    CALL add_stream_element (made_stream, 'rhhist_am', p3, p4=p4,        &
!!$                             units='', lrerun=.true., lpost=.false.,    &
!!$                             lav=.false.)
    mem => rhhist(_RI_XYZN_(:,:,:,acc),:)
    CALL new_channel_object(status, modstr//strmod, 'rhhist_am' &
     , mem = mem, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'rhhist_am', 'units', c='')
    CALL channel_halt(substr, status)

!!$    p4 => rhhist(:,:,cor,:,:)
!!$    CALL add_stream_element (made_stream, 'rhhist_cm', p3, p4=p4,        &
!!$                             units='', lrerun=.true., lpost=.false.,    &
!!$                             lav=.false.)
    mem => rhhist(_RI_XYZN_(:,:,:,cor),:)
    CALL new_channel_object(status, modstr//strmod, 'rhhist_cm' &
         , mem = mem, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'rhhist_cm', 'units', c='')
    CALL channel_halt(substr, status)

! op_ck_20130115+
    ! hydrophilic fraction of 3rd moment

    mem => philfrac(:,:,:,:,1)
    CALL new_channel_object(status, modstr//strmod, 'philfrac' &
         , p4 = philfrac_4d, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'philfrac', 'units', c='frac')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'philfrac' &
         , 'long_name', c='hydrophilic fraction of 3rd moment')
    CALL channel_halt(substr, status)

    mem => philfrac(_RI_XYZN_(:,:,:,akn),:)
    CALL new_channel_object(status, modstr//strmod, 'philfrac_km', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'philfrac_km', 'units', c='frac')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'philfrac_km' &
         , 'long_name', c='hydrophilic fraction of 3rd moment (Aitken mode)')
    CALL channel_halt(substr, status)

    mem => philfrac(_RI_XYZN_(:,:,:,acc),:)
    CALL new_channel_object(status, modstr//strmod, 'philfrac_am', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'philfrac_am', 'units', c='frac')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'philfrac_am' &
         , 'long_name', c='hydrophilic fraction of 3rd moment (acc. mode)')
    CALL channel_halt(substr, status)

    mem => philfrac(_RI_XYZN_(:,:,:,cor),:)
    CALL new_channel_object(status, modstr//strmod, 'philfrac_cm', mem = mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'philfrac_cm', 'units', c='frac')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'philfrac_cm' &
         , 'long_name', c='hydrophilic fraction of 3rd moment (coarse mode)')
    CALL channel_halt(substr, status)
! op_ck_20130115-

    ! (fixed) standard deviation of log-normal distributions

!!$    call add_stream_element (made_stream, 'sigma', sigma_str,           &
!!$                             nsize=nmod, repr=ARRAY1D,                  &
!!$                             longname='standard deviation of MADE modes', &
!!$                             lpost=.false., lav=.false., lrerun=.false.)
    CALL new_channel_object(status, modstr//strmod, &
         'sigma', p1=sigma_str, reprid=REPR_MADE_1D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//strmod, 'sigma' &
         , 'long_name', c='standard deviation of MADE modes')
    CALL channel_halt(substr, status)

    !--- Initialization of channel objects containing MADE parameter

    ! new channel objects will be reset to 0.0 automatically

    sigma_str = sigma

    CALL end_message_bi(modstr,'INITIALIZE MEMORY', substr)

  END SUBROUTINE made_init_memory

!-----------------------------------------------------------------------------

  SUBROUTINE made_init_coupling

    ! SUBROUTINE to check for existence of required coupling tracer

    USE messy_main_channel,       ONLY: get_channel_object
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                      , error_bi, info_bi, warning_bi
    USE messy_main_tracer_mem_bi, ONLY: ti_gp, GPTRSTR
    USE messy_main_tracer,        ONLY: get_tracer, OFF, ON

    ! MADE
    USE messy_made,               ONLY: modstr

    IMPLICIT NONE
    INTRINSIC TRIM

    CHARACTER(LEN=*), PARAMETER :: substr = 'made_init_coupling'
    LOGICAL                     :: lok      ! error status
    INTEGER                     :: status
    INTEGER                     :: idt = 0

    CALL start_message_bi(modstr, 'Initialize couplings', substr)

    IF (lcpl_gasphase) THEN

!      --- H2SO4(g) ---

       IF (TRIM(H2SO4_gas(1)) == '' .AND. TRIM(H2SO4_gas(2)) == '') THEN
          CALL info_bi('No gas-phase H2SO4 tracer specified,', &
               substr)
          CALL info_bi('looking for tracer H2SO4 in chemistry submodel '&
               //TRIM(chemmodule)//'...', substr)
          ! Check if tracer has already been defined
          CALL get_tracer(status, GPTRSTR, 'H2SO4', idx=idt)
          IF (status == 0 .AND. &
               TRIM(ti_gp(idt)%tp%ident%submodel) == TRIM(chemmodule)) THEN
             ! use existing tracer
             CALL info_bi('...using H2SO4 tracer defined by '&
                  //TRIM(ti_gp(idt)%tp%ident%submodel), substr)
          ELSE   ! abort
             CALL error_bi('Tracer H2SO4 not found in '//TRIM(chemmodule), &
                  substr)
          END IF
       ELSE
          CALL get_tracer(status, GPTRSTR, H2SO4_gas(1), &
                          subname=H2SO4_gas(2), idx=idt)
          lok = (status == 0)
          IF (lok) THEN
             IF (TRIM(ti_gp(idt)%tp%ident%submodel) == &
                  TRIM(chemmodule)) THEN
                CALL info_bi( &
                     'fetching chemistry module tracer idt', substr)
             ELSE
                CALL error_bi( &
                     'H2SO4 tracer not initialized by chosen chemistry module'&
                     , substr)
             ENDIF
          ELSE
             IF (TRIM(H2SO4_gas(2)) /= '') H2SO4_gas(1) = H2SO4_gas(1)//'_'
             CALL error_bi( &
                  'Tracer '//H2SO4_gas(1)//H2SO4_gas(2)//' not found.', substr)
          ENDIF
       ENDIF
       itracsulf = idt

!      --- gasphase production of H2SO4(g) (SO2 + OH ---> H2SO4) ---

       IF (TRIM(ProdSO4(1)) == '' .AND. TRIM(ProdSO4(2)) == '') THEN
          CALL warning_bi('WARNING: no gas-phase production of H2SO4 specified, using H2SO4 tendency instead!!!', substr)
          IF (itracsulf.GT.0) THEN
             itracprodso4 = itracsulf
          ELSE
             CALL error_bi( &
                  'error: also no gas-phase H2SO4 tracer available', substr)
          END IF
       ELSE
          CALL get_tracer(status, GPTRSTR, ProdSO4(1), subname=ProdSO4(2), &
                          idx=itracprodso4)
          lok = (status == 0)
          IF (lok) THEN
             CALL info_bi('fetching gas-phase H2SO4 production idt', substr)
          ELSE
             CALL error_bi('no production of H2SO4(g) available', substr)
          ENDIF
       ENDIF

!!      --- HCl(g) ---
!
!       IF (HCl_gas(1) /= '' .OR. HCl_gas(2) /= '') THEN
!          CALL get_tracer(status, GPTRSTR, HCl_gas(1), &
!                          subname=HCl_gas(2), idx=itrachcl)
!          lok = (status == 0)
!          IF (lok) THEN
!             IF (TRIM(ti_gp(itrachcl)%tp%ident%submodel) == &
!                  TRIM(chemmodule)) THEN
!                CALL info_bi('fetching chemistry module HCl tracer idt', &
!                     substr)
!             ELSE
!                CALL error_bi('HCl tracer not initialized by chosen chemistry module', &
!                     substr)
!             ENDIF
!          ELSE
!             IF (TRIM(HCl_gas(2)) /= '') HCl_gas(1) = HCl_gas(1)//'_'
!             CALL error_bi( &
!                  'Tracer '//HCl_gas(1)//HCl_gas(2)//' not found.', substr)
!!          ENDIF
!       ENDIF

!      --- HNO3(g) ---

       IF (TRIM(HNO3_gas(1)) /= '' .OR. TRIM(HNO3_gas(2)) /= '') THEN
          CALL get_tracer(status, GPTRSTR, HNO3_gas(1), &
                          subname=HNO3_gas(2), idx=itrachno3)
          lok = (status == 0)
          IF (lok) THEN
             IF (TRIM(ti_gp(itrachno3)%tp%ident%submodel) == &
                  TRIM(chemmodule)) THEN
                CALL info_bi( &
                     'fetching chemistry module HNO3 tracer idt', substr)
             ELSE
                CALL error_bi( &
                     'HNO3 tracer not initialized by chosen chemistry module' &
                     , substr)
             ENDIF
          ELSE
             IF (TRIM(HNO3_gas(2)) /= '') HNO3_gas(1) = HNO3_gas(1)//'_'
             CALL error_bi( &
                  'Tracer '//HNO3_gas(1)//HNO3_gas(2)//' not found.', substr)
          ENDIF
       ENDIF

!      --- NH3(g) ---

       IF  (TRIM(NH3_gas(1)) /= '' .OR. TRIM(NH3_gas(2)) /= '') THEN
          CALL get_tracer(status, GPTRSTR, NH3_gas(1), &
                          subname=NH3_gas(2), idx=itracnh3)
          lok = (status == 0)
          IF (lok) THEN
             IF (TRIM(ti_gp(itracnh3)%tp%ident%submodel) == &
                  TRIM(chemmodule)) THEN
                CALL info_bi( &
                     'fetching chemistry module NH3 tracer idt', substr)
             ELSE
                CALL error_bi( &
                     'NH3 tracer not initialized by chosen chemistry module' &
                     , substr)
             ENDIF
          ELSE
             IF (TRIM(NH3_gas(2)) /= '') NH3_gas(1) = NH3_gas(1)//'_'
             CALL error_bi( &
                  'Tracer '//NH3_gas(1)//NH3_gas(2)//' not found.', substr)
          ENDIF
       ENDIF
       !um_ak_20110720+
    ELSE ! note: itracsulf is always set
       itracprodso4 = itracsulf
       !um_ak_20110720+-
    ENDIF ! (lcpl_gasphase)

    IF (l_calc_emis) THEN
       IF (l_SS) THEN
!!$          ! get seasalt emission 
!!$          CALL get_stream(emission_stream, TRIM(SSemis_stream), ierr)
!!$          IF (ierr /= 0)  THEN
!!$             CALL finish (substr,&
!!$             'required emission stream for seasalt not available')
!!$          ELSE
!!$             CALL message(substr, 'fetching seasalt emission stream')
!!$          END IF
!!$          ! ... get accumulation mode mass emission flux
!!$          CALL get_stream_element(emission_stream, TRIM(SS_mass_as), &
!!$                                  Mss_as, ierr)
          CALL get_channel_object(status &
               , TRIM(SSemis_channel), TRIM(SS_mass_as), p2=Mss_as)
          IF (status /= 0) THEN 
             CALL error_bi( &
                          'channel_object for accumulation mode seasalt ' // &
                          'mass emission flux not available', substr)
          ELSE
              CALL info_bi('fetching seasalt mass(as) emission flux', substr)
          END IF
          ! ... get accumulation mode number emission flux
!!$          CALL get_stream_element(emission_stream, TRIM(SS_num_as), &
!!$                                  Nss_as, ierr)
          CALL get_channel_object(status &
               , TRIM(SSemis_channel), TRIM(SS_num_as), p2=Nss_as)
          IF (status /= 0) THEN 
             CALL error_bi('channel_object for accumulation mode ' // &
                  'seasalt number emission flux not available', substr)
          ELSE
             CALL info_bi('fetching seasalt number(as) emission flux', substr)
          END IF
          ! ... get coarse mode mass emission flux
!!$          CALL get_stream_element(emission_stream, TRIM(SS_mass_cs), &
!!$                                  Mss_cs, ierr)
          CALL get_channel_object(status &
               , TRIM(SSemis_channel), TRIM(SS_mass_cs), p2=Mss_cs)
          IF (status /= 0) THEN 
             CALL error_bi('channel_object for coarse mode seasalt ' // &
                  'mass emission fluxnot available', substr)
          ELSE
             CALL info_bi('fetching seasalt mass(cs) emission flux ', substr)
          END IF
!!$          CALL get_stream_element(emission_stream, TRIM(SS_num_cs), &
!!$                                  Nss_cs, ierr)
          CALL get_channel_object(status &
               , TRIM(SSemis_channel), TRIM(SS_num_cs), p2=Nss_cs)
          IF (status /=0 ) THEN 
             CALL error_bi('channel_object for coarse mode seasalt ' // &
                  'number emission flux not available', substr)
          ELSE
             CALL info_bi('fetching seasalt number(cs) emission flux', substr)
          END IF
       END IF  ! if (l_ss)

       IF (l_dust) THEN
!!$          CALL get_stream(emission_stream, TRIM(Duemis_stream), ierr)
!!$          IF (ierr /=0) &
!!$             CALL finish (substr, 'required emission stream for dust ' // &
!!$                          'not available')
!!$          ! ... get emission stream element for dust emissions
!!$          CALL get_stream_element(emission_stream, TRIM(emis_dust), &
!!$                                  dust_emis, ierr)
          CALL get_channel_object(status &
               , TRIM(Duemis_channel), TRIM(emis_dust), p2=dust_emis)
          IF (status /= 0) &
             CALL error_bi('channel_object for dust ' // &
                          'emission flux not available', substr)
       END IF ! if (l_dust)
    END IF ! if (l_calc_emiss)

    ! SOA "emissions"

    IF (trim(SOA_channel) /= '') THEN
       ! get SOA "emission" channel
!!$       CALL get_stream(SOA_emission_stream, TRIM(SOA_stream), ierr)
!!$       IF (ierr /= 0)  THEN
!!$          CALL finish (substr,&
!!$          'required emission stream for SOA not available')
!!$       ELSE
!!$          CALL message(substr, 'fetching SOA emission stream')
!!$       END IF
!!$       ! ... get SOA mass emission flux
!!$       CALL get_stream_element(SOA_emission_stream, TRIM(SOA_element), &
!!$                               Msoa, ierr)
       CALL get_channel_object(status &
            , TRIM(SOA_channel), TRIM(SOA_object), p3=Msoa)
!      IF (status /= 0) THEN 
!          CALL error_bi( &
!               'channel_object for SOA ' // &
!               'mass emission flux not available', substr)
!       ELSE
           CALL info_bi('fetching SOA mass emission flux', substr)
!       END IF
    END IF

    CALL end_message_bi(modstr,'CHECK NAMELIST SETTINGS', substr)
    
  END SUBROUTINE made_init_coupling
  
!=============================================================================!

SUBROUTINE made_vdiff

    ! MESSy
    USE messy_main_constants_mem, ONLY: g
    USE messy_main_grid_def_mem_bi,  ONLY: kproma, jrow, nlev
    USE messy_main_data_bi,       ONLY: pressi_3d
#ifdef ECHAM5
    USE messy_main_data_bi,       ONLY: pxtems
#endif
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, ti_gp
    USE messy_main_tracer,        ONLY: R_molarmass
   ! MADE
    USE messy_made,               ONLY: mwseas, mwsoil, mwair,       &
                                        du2ac0, du2corn,             &
                                        massfrac_du_j, massfrac_du_c

    IMPLICIT NONE

    ! LOCAL
!    CHARACTER(len=*), PARAMETER :: substr='made_vdiff'
    REAL(dp), POINTER :: zxtems(:,:)
    REAL(dp)          :: zdp(kproma) ! um_ak_20110720
    INTEGER           :: idt         ! um_ak_20110720

#ifdef ECHAM5
    zxtems => pxtems(:,1,:,jrow)
#endif
    IF (.not. l_calc_emis) RETURN

    ! zxtems: [mass]   = mol(x)/mol(air) * kg/m2/s
    !         [number] = #/mol(air)      * kg/m2/s

    ! sea salt

    ! um_ak_20110720+
    zdp(1:kproma) = pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) &
         - pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev)) 
    ! um_ak_20110720-
    IF (l_ss) THEN
#ifdef ECHAM5
       IF (.NOT. l_tendency) THEN ! um_ak_20110720
       ! MSS_as = accumulation mode mass flux [kg/m2/s]
       zxtems(1:kproma,itracseasj)   = zxtems(1:kproma,itracseasj)   &
                                       + Mss_as(1:kproma,jrow)*MWAIR/MWSEAS
       ! MSS_cs = coarse mode mass flux [kg/m2/s]
       zxtems(1:kproma,itracseasc)   = zxtems(1:kproma,itracseasc)   &
                                       + Mss_cs(1:kproma,jrow)*MWAIR/MWSEAS 
       ! NSS_as = accumulation mode number flux [#/m2/s]
       zxtems(1:kproma,itracac0)     = zxtems(1:kproma,itracac0)     &
                                       + Nss_as(1:kproma,jrow)*MWAIR*1.0e-3_dp
       ! Nss_cs = coarse mode number flux [#/m2/s]
       zxtems(1:kproma,itraccorn)    = zxtems(1:kproma,itraccorn)    &
                                       + Nss_cs(1:kproma,jrow)*MWAIR*1.0e-3_dp

! um_ak_20110720+
       ELSE
#endif
          ! MSS_cs = accumulation mode mass flux [kg/m2/s]
          idt = itracseasj
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               Mss_as(1:kproma,jrow) / zdp(1:kproma) * g /       &
               ti_gp(itracseasj)%tp%meta%cask_r(R_molarmass) * MWAIR 
          idt = itracseasc
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               Mss_cs(1:kproma,jrow) / zdp(1:kproma) * g /       &
               ti_gp(itracseasc)%tp%meta%cask_r(R_molarmass) * MWAIR 
          idt = itracac0
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               Nss_as(1:kproma,jrow) / zdp(1:kproma) * g * MWAIR * 1.E-3
          idt = itraccorn
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               Nss_cs(1:kproma,jrow) / zdp(1:kproma) * g * MWAIR * 1.E-3
       ENDIF
! um_ak_20110720-
#ifdef ECHAM5
    END IF
#endif
    ! mineral dust

    IF (l_dust) THEN
#ifdef ECHAM5
       IF (.NOT. l_tendency) THEN
       ! dust_emis = dust total mass flux [kg/m2/s]
       zxtems(1:kproma,itracdustj)   = zxtems(1:kproma,itracdustj)           &
                                       + dust_emis(1:kproma,jrow)            &
                                       * massfrac_du_j * MWAIR/MWSOIL
       zxtems(1:kproma,itracdustc)   = zxtems(1:kproma,itracdustc)           &
                                       + dust_emis(1:kproma,jrow)            &
                                       * massfrac_du_c * MWAIR/MWSOIL
       ! Number emissions [#/m2/s]
       zxtems(1:kproma,itracac0)     = zxtems(1:kproma,itracac0)            &
                                       + dust_emis(1:kproma,jrow) * DU2AC0
       zxtems(1:kproma,itraccorn)    = zxtems(1:kproma,itraccorn)           &
                                       + dust_emis(1:kproma,jrow) * DU2CORN
! um_ak_20110720+
       ELSE
#endif
          ! MSS_cs = accumulation mode mass flux [kg/m2/s]
          idt = itracdustj
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +             &
               dust_emis(1:kproma,jrow)* massfrac_du_j / zdp(1:kproma) * g / &
               ti_gp(itracdustj)%tp%meta%cask_r(R_molarmass) * MWAIR
          idt = itracdustc
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               dust_emis(1:kproma,jrow) * massfrac_du_c / zdp(1:kproma) * g / &
               ti_gp(itracdustc)%tp%meta%cask_r(R_molarmass) * MWAIR 
          idt = itracac0
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               dust_emis(1:kproma,jrow) * DU2AC0 / zdp(1:kproma) * g 
          idt = itraccorn
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               dust_emis(1:kproma,jrow) * DU2CORN / zdp(1:kproma) * g 
       ENDIF
! um_ak_20110720-
#ifdef ECHAM5
    END IF
#endif
END SUBROUTINE made_vdiff

!=============================================================================!

SUBROUTINE made_physc

  USE messy_main_tools,         ONLY: jptlucu1, jptlucu2, tlucuaw
  USE messy_main_data_bi,       ONLY: press_3d,pressi_3d,      &
                                      tm1_3d, tte_3d,          &
                                      qm1_3d, qte_3d,          &
                                      aclc
  USE messy_main_grid_def_mem_bi, ONLY: kproma, jrow, nlev,      &
                                      nproma
  USE messy_main_timer,         ONLY: time_step_len
  USE messy_main_tracer_mem_bi, ONLY: qxtm1, qxtte
  USE messy_main_constants_mem, ONLY: M_air, rd, vtmpc1
  ! MADE
  USE messy_made,               ONLY: nspcsda, akn, acc, cor,          &
                                      bc_i, bc2_i, bc_j, bc2_j,        &
                                      oc_i, oc2_i, oc_j, oc2_j,        &
! op_ck_20130116+
                                      anthfac, orgfac, soilfac,        &
                                      vnu3, vac3, vcor3,               &
! op_ck_20130116-
                                      vso4ai, vnh4ai, vno3ai, vh2oai,  &
                                      veci, vorgpai, vseasi, vnu0,     &
                                      vso4aj, vnh4aj, vno3aj, vh2oaj,  &
                                      vecj, vorgpaj, vseasj, vdustj,   &
                                      vac0, vseasc, vdustc, vh2oac,    &
                                      vcorn, vsulf, vhno3, vnh3,       &
!                                      vhcl,                            &
                                      MWH2SO4, MWORG,                  &
!DEBUG+
!!$                                      GRAV, AVO
                                      GRAV, AVO, rhoseas, rhosoil,     &
                                      rhoso4, rhono3, rhonh4
!DEBUG-

  ! SUBROUTINES
  USE messy_made,               ONLY: made_main
                                   
  IMPLICIT NONE
  INTRINSIC ASSOCIATED, EPSILON, INT

  !--- Parameter list:

  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'made_physc'
  REAL(dp)    :: zpress(nproma)           ! full level pressure [Pa]
  REAL(dp)    :: ztemp(nproma)            ! temperature [K]
  REAL(dp)    :: zrelhum(nproma)          ! rel. humidity (0-1) [frac]
  REAL(dp)    :: zspechum_m(nproma)       ! specific humidity [kg/kg],
                                          ! grid box mean
  REAL(dp)    :: zspechum_0               ! specific humidity [kg/kg],
                                          ! cloud free area
  REAL(dp)    :: zqs                      ! saturation specific humidity
  REAL(dp)    :: zrhoair(nproma)          ! density of air [kg/m3]
  REAL(dp)    :: zclcover(nproma)         ! cloud cover [frac]
  REAL(dp)    :: zso4rat(nproma)          ! H2SO4(g) production rate
                                          ! [ug/m3/s]
  REAL(dp)    :: zsoa(nproma)             ! SOA "emissions" [ug/m3/s]
  REAL(dp)    :: zrh_hist_akn(nproma)     ! RH history, Aitken mode
  REAL(dp)    :: zrh_hist_acc(nproma)     ! RH history, acc. mode
  REAL(dp)    :: zrh_hist_cor(nproma)     ! RH history, coarse mode
  REAL(dp)    :: zcblk(nproma,nspcsda)    ! tracer conc. array for MADE
  REAL(dp)    :: zcblk_m1(nproma,nspcsda) ! tracer conc. array before MADE
  REAL(dp)    :: ztmst, zqtmst
  INTEGER     :: jl,jk
  INTEGER     :: it

  ! Additional MADE output

  REAL(dp)    :: zdgnuc(nproma)        ! geom. mean diameter (wet),
                                       ! Aitken mode [m]
  REAL(dp)    :: zdgacc(nproma)        ! geom. mean diameter (wet),
                                       ! accumulation mode [m]
  REAL(dp)    :: zdgcor(nproma)        ! geom. mean diameter (wet),
                                       ! coarse mode [m]
  REAL(dp)    :: zdgdrynuc(nproma)     ! geom. mean diameter (dry),
                                       ! Aitken mode [m]
  REAL(dp)    :: zdgdryacc(nproma)     ! geom. mean diameter (dry),
                                       ! accumulation mode [m]
  REAL(dp)    :: zdgdrycor(nproma)     ! geom. mean diameter (dry),
                                       ! coarse mode [m]
  REAL(dp)    :: zdensnuc(nproma)      ! average density (wet),
                                       ! Aitken mode [kg/m3]
  REAL(dp)    :: zdensacc(nproma)      ! average density (wet),
                                       ! accumulation mode [kg/m3]
  REAL(dp)    :: zdenscor(nproma)      ! average density (wet),
                                       ! coarse mode [kg/m3]

  ! (auxiliary) variables for calculating the BC/OC split

  REAL(dp) :: splitfac(nproma,num_split_facs)
  REAL(dp) :: bc_mass_i, bc2_mass_i, bc_mass_j, bc2_mass_j
  REAL(dp) :: oc_mass_i, oc2_mass_i, oc_mass_j, oc2_mass_j
  REAL(dp) :: bc_total_i, bc_total_j
  REAL(dp) :: oc_total_i, oc_total_j
! op_ck_20130115+

  ! (auxiliary) variables for calculating the hydrophilic fraction of 3rd mom.
  REAL(dp) :: nu3phob, ac3phob, cor3phob
! op_ck_20130115-

  INTEGER :: idt, jm ! um_ak_20110720
  !--- 0. Initializations. -----------------------------------------------

  ztmst  = time_step_len
  zqtmst = 1./time_step_len

  ! 1. Calculate soot aging (BC, POM). -----------------------------------

  IF (l_sootag) THEN
  CALL made_sootaging
  END IF

 ! 2. Calculate actual values. --------------------------------------

  DO jk=1,nlev

     !--- 2.1 Meteorological data. --------------------------------------

     ! actual temperature [K]
     ztemp(1:kproma)    = tm1_3d(_RI_XYZ__(1:kproma,jrow,jk)) &
                          + tte_3d(_RI_XYZ__(1:kproma,jrow,jk)) * ztmst

     ! actual pressure [Pa]
     zpress(1:kproma)   = press_3d(_RI_XYZ__(1:kproma,jrow,jk))

     ! actual specific humidity [kg/kg], grid box mean
     zspechum_m(1:kproma) = MAX(EPSILON(1.0_dp),qm1_3d(_RI_XYZ__(1:kproma,jrow,jk)) &
                            + qte_3d(_RI_XYZ__(1:kproma,jrow,jk)) * ztmst)

     ! cloud cover [frac]
!!!     zclcover(1:kproma) = MIN(aclc(1:kproma,jk,jrow), 1.0_dp - EPSILON(1.0))
     zclcover(1:kproma) = MIN(aclc(_RI_XYZ__(1:kproma,jrow,jk)), 0.9999_dp)

     ! actual relative humidity [frac]

     DO jl = 1, kproma
        ! zqs: actual saturation specific humidity (liquid water)
        !      (code adopted from "radiation.f90")
        it  = INT(ztemp(jl) * 1000.0)
        it  = MAX(MIN(it,jptlucu2),jptlucu1)
        zqs = tlucuaw(it) / zpress(jl)
        zqs = MIN(zqs,0.5_dp)
        zqs = zqs / (1.0_dp - vtmpc1 * zqs)

        ! beta test: calculate relative humidity in cloud free area
        !            of grid box instead of average for whole grid box
        !
        ! Q_m = specific humidity, grid box mean
        ! Q_0 = specific humidity, cloud free area of grid box
        ! Q_s = specific humidity, cloudy area of grid box
        !     = saturation specific humidity (liquid water),
        !       this is an assumption for the sake of simplicity
        !
        ! Q_m = (1-cloud_cover)*Q_0 + cloud_cover*Qs
        !
        ! ==> Q_0 = Q_m - cloud_cover*Q_s / (1-cloud_cover)

        zspechum_0 = zspechum_m(jl) - zclcover(jl) * zqs &
                     / (1.0_dp - zclcover(jl))

        ! actual relative humidity [frac], cloud free area
        zrelhum(jl) = zspechum_0 / zqs

!!!        ! actual relative humidity [frac], grid box mean
!!!        zrelhum(jl) = zspechum_m(jl) / zqs

        ! beta test: limit relative humidity to (1%...99%)
        zrelhum(jl) = MAX(0.01_dp,MIN(zrelhum(jl),0.99_dp))
     END DO

      ! actual density of air [kg/m3]
     zrhoair(1:kproma)  = zpress(1:kproma) / (RD * ztemp(1:kproma) &
                            * (1.0_dp + vtmpc1 * zspechum_m(1:kproma)))

     ! relative humidity history ---> hysteresis
     ! um_ak_20110720+
!!$     zrh_hist_akn(1:kproma) = rhhist_4d(1:kproma,jk,akn,jrow) ! Aitken mode
!!$     zrh_hist_acc(1:kproma) = rhhist_4d(1:kproma,jk,acc,jrow) ! acc. mode
!!$     zrh_hist_cor(1:kproma) = rhhist_4d(1:kproma,jk,cor,jrow) ! coarse mode
     jm = akn
     zrh_hist_akn(1:kproma) = rhhist_4d(_RI_XYZN_(1:kproma,jrow,jk,jm)) ! Aitken mode
     jm = acc
     zrh_hist_acc(1:kproma) = rhhist_4d(_RI_XYZN_(1:kproma,jrow,jk,jm)) ! acc. mode
     jm = cor
     zrh_hist_cor(1:kproma) = rhhist_4d(_RI_XYZN_(1:kproma,jrow,jk,jm)) ! coarse mode
     ! um_ak_20110720-

     ! --- 2.2 Chemical data. ---------------------------------------------

     if ((jk.eq.nlev).and.(ASSOCIATED(Msoa))) then
        ! convert from [molecules/m2/s] to [ug/m3/s]
! um_ak_20110720 zsoa(1:kproma) = Msoa(1:kproma,1,jrow) * GRAV &
        zsoa(1:kproma) = Msoa(_RI_XYZ__(1:kproma,jrow,1)) * GRAV &
             / (pressi_3d(_RI_XYZ__(1:kproma,jrow,jk+1)) - pressi_3d(_RI_XYZ__(1:kproma,jrow,jk))) &
             * zrhoair(1:kproma) * 1.0e6 * MWORG / AVO
     else
        zsoa(1:kproma)    = 0.0_dp
     end if

     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Production of gas-phase H2SO4: SO2 + OH ---> H2SO4(g)
     ! If no H2SO4(g) production is specified via KPP, the tendency
     ! if the H2SO4(g) tracer will be used as first guess.
     ! In this case, itracprodso4 = itracsulf.

     ! [mol/mol/s] ---> [ug(H2SO4)/m3/s]
     idt = itracprodso4
     zso4rat(1:kproma) = MAX(0.0_dp, qxtte(_RI_X_ZN_(1:kproma,jk,idt)) * &
          MWH2SO4/M_air * zrhoair(1:kproma) * 1.0e9)

     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     ! Assign ECHAM tracers to aerosol array *CBLK* used
     ! inside of aerosol dynamics calculations.

     zcblk(1:kproma,:)    = 1.0e-30_dp
     zcblk_m1(1:kproma,:) = 1.0e-30_dp

     DO jl=1,kproma

        !--- 2.2 Gases. ---------------------------------------------------

        idt = itracsulf
        zcblk(jl,vsulf) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                          + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itrachno3
        zcblk(jl,vhno3) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                          + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itracnh3
        zcblk(jl,vnh3)  = qxtm1(_RI_X_ZN_(jl,jk,idt))  &
                          + qxtte(_RI_X_ZN_(jl,jk,idt))  * ztmst
!        zcblk(jl,vhcl)  = qxtm1(jl,jk,itrachcl)  &
!                          + qxtte(jl,jk,itrachcl)  * ztmst

        !--- 2.3 Particle mass. -------------------------------------------

        ! Sulfate mass:

        idt = itracso4ai
        zcblk(jl,vso4ai) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itracso4aj
        zcblk(jl,vso4aj) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst

        ! Ammonium mass:

        idt = itracnh4ai
        zcblk(jl,vnh4ai) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) *ztmst
        idt = itracnh4aj
        zcblk(jl,vnh4aj) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) *ztmst

        ! Nitrate mass:

        idt = itracno3ai
        zcblk(jl,vno3ai) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itracno3aj
        zcblk(jl,vno3aj) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst

        ! Black Carbon mass:
        ! mode i
        idt = itraceci
        bc_mass_i  = max(1.0e-30_dp, qxtm1(_RI_X_ZN_(jl,jk,idt))    &
                                     + qxtte(_RI_X_ZN_(jl,jk,idt))  * ztmst)
        idt = itracec2i
        bc2_mass_i = max(1.0e-30_dp, qxtm1(_RI_X_ZN_(jl,jk,idt))   &
                                     + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst)

        zcblk(jl,veci) = bc_mass_i + bc2_mass_i

        bc_total_i = max(2.0e-30_dp, zcblk(jl,veci))

        splitfac(jl,bc_i)  = bc_mass_i  / bc_total_i
        splitfac(jl,bc2_i) = bc2_mass_i / bc_total_i

        ! mode j
        idt = itracecj
        bc_mass_j  = max(1.0e-30_dp, qxtm1(_RI_X_ZN_(jl,jk,idt))    &
                                     + qxtte(_RI_X_ZN_(jl,jk,idt))  * ztmst)
        idt = itracec2j
        bc2_mass_j = max(1.0e-30_dp, qxtm1(_RI_X_ZN_(jl,jk,idt))   &
                                     + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst)

        zcblk(jl,vecj) = bc_mass_j + bc2_mass_j

        bc_total_j = max(2.0e-30_dp, zcblk(jl,vecj))

        splitfac(jl,bc_j)  = bc_mass_j  / bc_total_j
        splitfac(jl,bc2_j) = bc2_mass_j / bc_total_j

        ! Organic matter mass:
        ! Aitken mode (i)
        idt = itracorgpai
        oc_mass_i  = max(1.0e-30_dp, qxtm1(_RI_X_ZN_(jl,jk,idt))   &
                                    + qxtte(_RI_X_ZN_(jl,jk,idt))  * ztmst)
        idt = itracorgpa2i
        oc2_mass_i = max(1.0e-30_dp, qxtm1(_RI_X_ZN_(jl,jk,idt))  &
                                    + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst)

        zcblk(jl,vorgpai) = oc_mass_i + oc2_mass_i

        oc_total_i = max(2.0e-30_dp, zcblk(jl,vorgpai))

        splitfac(jl,oc_i)  = oc_mass_i  / oc_total_i
        splitfac(jl,oc2_i) = oc2_mass_i / oc_total_i

        ! accumulation mode (j)
        idt = itracorgpaj
        oc_mass_j  = max(1.0e-30_dp, qxtm1(_RI_X_ZN_(jl,jk,idt))    &
                                     + qxtte(_RI_X_ZN_(jl,jk,idt))  * ztmst)
        idt = itracorgpa2j
        oc2_mass_j = max(1.0e-30_dp, qxtm1(_RI_X_ZN_(jl,jk,idt))   &
                                     + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst)

        zcblk(jl,vorgpaj) = oc_mass_j + oc2_mass_j

        oc_total_j = max(2.0e-30_dp, zcblk(jl,vorgpaj))

        splitfac(jl,oc_j)  = oc_mass_j  / oc_total_j
        splitfac(jl,oc2_j) = oc2_mass_j / oc_total_j

        ! Sea salt mass:

        idt = itracseasi
        zcblk(jl,vseasi) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itracseasj
        zcblk(jl,vseasj) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itracseasc
        zcblk(jl,vseasc) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst

        ! Mineral dust mass:

        idt = itracdustj
        zcblk(jl,vdustj) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itracdustc
        zcblk(jl,vdustc) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst

        ! Aerosol water mass (just for zclkb_m1):

        idt = itrach2oai
        zcblk(jl,vh2oai) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itrach2oaj
        zcblk(jl,vh2oaj) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
        idt = itrach2oac
        zcblk(jl,vh2oac) = qxtm1(_RI_X_ZN_(jl,jk,idt)) &
                           + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst

        !--- 2.4 Particle numbers. -----------------------------------------

        idt = itracnu0
        zcblk(jl,vnu0)   = qxtm1(_RI_X_ZN_(jl,jk,idt))   &
                           + qxtte(_RI_X_ZN_(jl,jk,idt))   * ztmst
        idt = itracac0
        zcblk(jl,vac0)   = qxtm1(_RI_X_ZN_(jl,jk,idt))   &
                           + qxtte(_RI_X_ZN_(jl,jk,idt))   * ztmst
        idt = itraccorn
        zcblk(jl,vcorn)  = qxtm1(_RI_X_ZN_(jl,jk,idt))  &
                           + qxtte(_RI_X_ZN_(jl,jk,idt))  * ztmst

     END DO ! jl-loop

     ! Make a copy of *CBLK* for calculating the updated ECHAM tracer
     ! concentrations after aerosol dynamics calculations.

     zcblk_m1(1:kproma,:) = zcblk(1:kproma,:)

     ! 3. Convert from "ECHAM units" to "MADE units" and check for
     !    valid concentrations. -----------------------------------------

     CALL made_echam2made ( nproma, kproma, zrhoair, zcblk )

     ! 4. Calculate aerosol dynamics ---> call MADE. --------------------

     IF (l_main) THEN
     CALL made_main ( nproma, kproma, zpress, ztemp, zrelhum, &
                      ztmst, zso4rat, zsoa, zclcover, zcblk,      &
                      zrh_hist_akn, zrh_hist_acc, zrh_hist_cor,   &
                      zdgnuc, zdgacc, zdgcor,                     &
                      zdgdrynuc, zdgdryacc, zdgdrycor,            &
                      zdensnuc, zdensacc, zdenscor )
     ELSE
        IF (l_setrhhist) THEN
           zrh_hist_akn = 1._dp
           zrh_hist_acc = 1._dp
           zrh_hist_cor = 1._dp
        END IF
        zdgnuc       = 1.e-8_dp
        zdgacc       = 1.e-7_dp
        zdgcor       = 1.e-6_dp
        zdgdrynuc    = 1.e-8_dp
        zdgdryacc    = 1.e-7_dp
        zdgdrycor    = 1.e-9_dp
        zdensnuc     = (rhoso4 + rhono3 + 3.0_dp * rhonh4) / 5.0_dp
        zdensacc     = zdensnuc
        zdenscor     = (rhoseas + rhosoil) / 2.0_dp
     END IF

! op_ck_20130115+
     ! 4.1 Calculate hydrophilic fraction of 3rd moment (used for
     !     scavenging).

     DO jl = 1, kproma
        nu3phob  = MAX(0.0_dp,  anthfac*zcblk(jl,veci)   *splitfac(jl,bc2_i) &
                              + orgfac *zcblk(jl,vorgpai)*splitfac(jl,oc2_i))
        ac3phob  = MAX(0.0_dp,  anthfac*zcblk(jl,vecj)   *splitfac(jl,bc2_j) &
                              + orgfac *zcblk(jl,vorgpaj)*splitfac(jl,oc2_j) &
                              + soilfac*zcblk(jl,vdustj))
        cor3phob = MAX(0.0_dp,  soilfac*zcblk(jl,vdustc))
        
        jm = akn
        philfrac_4d(_RI_XYZN_(jl,jrow,jk,jm)) = 1.0_dp - nu3phob  / &
             MAX(1.0e-30_dp,zcblk(jl,vnu3))
        jm = acc
        philfrac_4d(_RI_XYZN_(jl,jrow,jk,jm)) = 1.0_dp - ac3phob  / &
             MAX(1.0e-30_dp,zcblk(jl,vac3))
        jm = cor
        philfrac_4d(_RI_XYZN_(jl,jrow,jk,jm)) = 1.0_dp - cor3phob / &
             MAX(1.0e-30_dp,zcblk(jl,vcor3))
     END DO
! op_ck_20130115-

     ! 5. Convert from "MADE units" to "ECHAM units". -------------------

     CALL made_made2echam ( nproma, kproma, zrhoair, zcblk )

     ! 6. Update ECHAM tracer tendencies for current level. -------------

     IF (l_feedback) THEN
     CALL made_updatexte ( nproma, kproma, zcblk, zcblk_m1, &
                           splitfac, jk, zqtmst )
     END IF

     ! 7. Store MADE quantities in channels. -----------------------------

     DO jl=1,kproma

        ! Aitken mode
        ! um_ak_20110720+
        jm = akn
!!$        wetrad_4d(jl,jk,akn,jrow)  = zdgnuc(jl) / 2.0_dp    ! [m]
!!$        dryrad_4d(jl,jk,akn,jrow)  = zdgdrynuc(jl) / 2.0_dp ! [m]
!!$        densaer_4d(jl,jk,akn,jrow) = zdensnuc(jl)           ! [kg/m3]
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm))  = zdgnuc(jl) / 2.0_dp    ! [m]
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm))  = zdgdrynuc(jl) / 2.0_dp ! [m]
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zdensnuc(jl)           ! [kg/m3]
        ! um_ak_20110720-

        ! accumulation mode
       ! um_ak_20110720+
!!$        wetrad_4d(jl,jk,acc,jrow)  = zdgacc(jl) / 2.0_dp    ! [m]
!!$        dryrad_4d(jl,jk,acc,jrow)  = zdgdryacc(jl) / 2.0_dp ! [m]
!!$        densaer_4d(jl,jk,acc,jrow) = zdensacc(jl)           ! [kg/m3]
        jm = acc
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm))  = zdgacc(jl) / 2.0_dp    ! [m]
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm))  = zdgdryacc(jl) / 2.0_dp ! [m]
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zdensacc(jl)           ! [kg/m3]
       ! um_ak_20110720-

        ! coarse mode
        ! um_ak_20110720+
!!$        wetrad_4d(jl,jk,cor,jrow)  = zdgcor(jl) / 2.0_dp    ! [m]
!!$        dryrad_4d(jl,jk,cor,jrow)  = zdgdrycor(jl) / 2.0_dp ! [m]
!!$        densaer_4d(jl,jk,cor,jrow) = zdenscor(jl)           ! [kg/m3]
        jm = cor
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm))  = zdgcor(jl) / 2.0_dp    ! [m]
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm))  = zdgdrycor(jl) / 2.0_dp ! [m]
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zdenscor(jl)           ! [kg/m3]
       ! um_ak_20110720+

        ! relative humidity history
       ! um_ak_20110720+
!!$        rhhist_4d(jl,jk,akn,jrow) = zrh_hist_akn(jl)   ! Aitken mode
!!$        rhhist_4d(jl,jk,acc,jrow) = zrh_hist_acc(jl)   ! acc. mode
!!$        rhhist_4d(jl,jk,cor,jrow) = zrh_hist_cor(jl)   ! coarse mode
        jm = akn
        rhhist_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrh_hist_akn(jl)   ! Aitken mode
        jm = acc
        rhhist_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrh_hist_acc(jl)   ! acc. mode
        jm = cor
        rhhist_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrh_hist_cor(jl)   ! coarse mode
       ! um_ak_20110720-

     END DO ! jl-loop
  END DO ! jk-loop

END SUBROUTINE made_physc

!-----------------------------------------------------------------------------

SUBROUTINE made_free_memory

IMPLICIT NONE

INTRINSIC ASSOCIATED

IF (ASSOCIATED(wetradius)) DEALLOCATE(wetradius)
IF (ASSOCIATED(dryradius)) DEALLOCATE(dryradius)
IF (ASSOCIATED(densaer))   DEALLOCATE(densaer)
IF (ASSOCIATED(rhhist))    DEALLOCATE(rhhist)
! op_ck_20130115+
IF (ASSOCIATED(philfrac))  DEALLOCATE(philfrac)
! op_ck_20130115-

END SUBROUTINE made_free_memory

!-----------------------------------------------------------------------------

  SUBROUTINE made_read_nml_cpl(status, iou)

    ! MADE MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    ! read namelist for 'coupling' to online tracers

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_made,       ONLY: modstr

    IMPLICIT NONE
    INTRINSIC TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'made_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status


    ! INITIALIZE
    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    WRITE(*,*) ''
    WRITE(*,*) '--------------------------------------------------------'
    WRITE(*,*) '--------------------------------------------------------'
    WRITE(*,*) '---  CPL namelist: settings for aerosol module MADE  ---'
    WRITE(*,*) '---'
    WRITE(*,'(1X,A,L1)') '---  lcpl_gasphase             = ', lcpl_gasphase
    IF (lcpl_gasphase) THEN
       WRITE(*,'(4X,A,A23)') '+ chemmodule                = ', chemmodule
       WRITE(*,'(4X,A,A23)') '+ H2SO4_gas                 = ', H2SO4_gas(1)
       WRITE(*,'(4X,A,A23)') '+ NH3_gas                   = ', NH3_gas(1)
!       WRITE(*,'(4X,A,A23)') '+ HCl_gas                   = ', HCl_gas(1)
       WRITE(*,'(4X,A,A23)') '+ HNO3_gas                  = ', HNO3_gas(1)
       IF (ProdSO4(1) == '' .AND. ProdSO4(2) == '') THEN
          WRITE(*,'(4X,A,A23)') '+ ProdSO4                   = ', &
                                                    'tendency from H2SO4-tracer'
       ELSE
          WRITE(*,'(4X,A,A23)') '+ ProdSO4                   = ', ProdSO4(1)
       END IF
    END IF ! if (lcpl_gasphase)
    WRITE(*,'(1X,A,A23)') '---  SOA_channel               = ', SOA_channel
    IF (trim(SOA_channel) /= '') then
       WRITE(*,'(4X,A,A23)') '+ SOA_object                = ', SOA_object
    END IF
    WRITE(*,'(1X,A,L1)') '---  l_calc_emis               = ', l_calc_emis
    WRITE(*,'(1X,A,L1)') '---  l_tendency                = ', l_tendency
    IF (l_calc_emis) THEN
       WRITE(*,'(1X,A,L1)') '---  l_ss                      = ', l_ss
       IF (l_ss) THEN
          WRITE(*,'(4X,A,A23)') '+ SSemis_channel            = ', SSemis_channel
          WRITE(*,'(4X,A,A23)') '+ SS_mass_as                = ', SS_mass_as
          WRITE(*,'(4X,A,A23)') '+ SS_num_as                 = ', SS_num_as
          WRITE(*,'(4X,A,A23)') '+ SS_mass_cs                = ', SS_mass_cs
          WRITE(*,'(4X,A,A23)') '+ SS_num_cs                 = ', SS_num_cs
       END IF ! if (l_ss)
       WRITE(*,'(1X,A,L1)')   '---  l_dust                    = ', l_dust
       IF (l_dust) THEN
          WRITE(*,'(4X,A,A23)') '+ Duemis_channel            = ', Duemis_channel
          WRITE(*,'(4X,A,A23)') '+ emis_dust                 = ', emis_dust
       END IF ! if (l_dust)
    ELSE
       l_ss   = .false.
       l_dust = .false.
    END IF ! if (l_calc_emis)
    WRITE(*,*) '--------------------------------------------------------'
    WRITE(*,*) '--------------------------------------------------------'
    WRITE(*,*) ''

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE made_read_nml_cpl

!-----------------------------------------------------------------------------

   SUBROUTINE made_sootaging

   !**** *SUBROUTINE* *MADE_SOOTAGING*
   !
   !     PURPOSE.
   !     --------
   !
   !         CALCULATE SOOT AGING OF BC AND POM.
   !
   !**   INTERFACE.
   !     ----------
   !
   !       *SOOTAGING*
   !       *PARAMETERS (I=input, O=output):
   !
   !         NONE.
   !
   !     METHODS.
   !     --------
   !
   !         SOOT AGING (TRANSFORMATION FROM HYDROPHOBIC TO HYDROPHILIC)
   !         IS PARAMETRIZED BY A SIMPLE EXPONENTIAL DECAY USING THE
   !         SPECIFIED HALF-LIFE TIME. THE TRANSFORMED BC/POM IS
   !         TRANSFERED FROM HYDROPHOIC TRACER TO THE CORRESPONDING
   !         HYDROPHILIC TRACER.
   !
   !     EXTERNALS.
   !     ----------
   !
   !         NONE.
   !
   !     REFERENCE.
   !     ----------
   !
   !         METHOD ADOBTED FROM Lohmann, U., J. Feichter, C. C. Chuang,
   !         J. E. Penner: Prediction of the number of cloud droplets
   !         in the ECHAM GCM, J. Geophys. Res., Vol. 104, No. D8,
   !         pp. 9169-9198, 1999.

   USE messy_made,               ONLY: bctime, pomtime
   USE messy_main_tracer_mem_bi, ONLY: qxtm1, qxtte
   USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev
   USE messy_main_timer,         ONLY: time_step_len

   IMPLICIT NONE
   INTRINSIC EXP, LOG

   ! local variables

   INTEGER  :: JK, JL                 ! loop index

   REAL(dp) :: ZFAC                   ! aux. variable
   REAL(dp) :: ZDECAY_BC, ZDECAY_POM  ! decay factor BC, POM
   REAL(dp) :: ZXTP1                  ! new BC/POM conc. (hydrophobic)
   REAL(dp) :: ZDXTDT                 ! change of tendency

   REAL(dp) :: ZTMST                  ! time step [s]

   INTEGER :: idt ! um_ak_20110720
   ! ------------------------------------------------------------------

   ztmst      = time_step_len

   zfac       = LOG(0.5_dp) * ztmst
   zdecay_bc  = (1.0_dp - EXP(zfac/bctime))  / ztmst
   zdecay_pom = (1.0_dp - EXP(zfac/pomtime)) / ztmst

   DO jk = 1,nlev
      DO jl = 1,kproma

         ! BC, Aitken mode
         ! um_ak_20110720+
         idt = itracec2i
         !zxtp1  = qxtm1(jl,jk,itracec2i) + qxtte(jl,jk,itracec2i) * ztmst
         zxtp1  = qxtm1(_RI_X_ZN_(jl,jk,idt)) + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
         ! um_ak_20110720+
         zdxtdt = zxtp1 * zdecay_bc
         ! um_ak_20110720+
         !qxtte(jl,jk,itracec2i) = qxtte(jl,jk,itracec2i) - zdxtdt       ! phob
         !qxtte(jl,jk,itraceci)  = qxtte(jl,jk,itraceci)  + zdxtdt       ! phil
         qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) - zdxtdt  ! phob
         idt = itraceci
         qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) + zdxtdt  ! phil
         ! um_ak_20110720+

         ! BC, accumulation mode
         idt = itracec2j
         zxtp1  = qxtm1(_RI_X_ZN_(jl,jk,idt)) + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
         zdxtdt = zxtp1 * zdecay_bc
         qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) - zdxtdt  ! phob
         idt = itracecj
         qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) + zdxtdt  ! phil

         ! POM, Aitken mode
         idt = itracorgpa2i
         zxtp1  = qxtm1(_RI_X_ZN_(jl,jk,idt)) + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
         zdxtdt = zxtp1 * zdecay_pom
         qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) - zdxtdt ! phob
         idt = itracorgpai
         qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) + zdxtdt ! phil

         ! POM, accumulation mode
         idt = itracorgpa2j
         zxtp1  = qxtm1(_RI_X_ZN_(jl,jk,idt)) + qxtte(_RI_X_ZN_(jl,jk,idt)) * ztmst
         zdxtdt = zxtp1 * zdecay_pom
         qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) - zdxtdt ! phob
         idt = itracorgpaj
         qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) + zdxtdt ! phil

      END DO
   END DO

   END SUBROUTINE made_sootaging

!-----------------------------------------------------------------------------

   SUBROUTINE made_echam2made ( blksize, numcells, rhoair, cblk )

   ! ----------------------------------------------------------------------
   ! *** MADE_ECHAM2MADE:
   !
   !     Convert from "ECHAM units" [mol/mol] and [#/mol] to
   !     "MADE units" [ug/m3] and [#/m3] respectively.
   !
   !     [mol/mol] = SO4AI, SO4AJ
   !                 ECI, ECJ
   !                 ORGPAI, ORGPAJ
   !                 NH4AI, NH4AJ
   !                 NO3AI, NO3AJ
   !                 SEASI, SEASJ, SEASC
   !                 DUSTJ, DUSTC
   !               = H2SO4(g)
   !                 HNO3(g)
   !                 NH3(g)
!   !                 HCl(g)
   !                  
   !     [ug/m3]   = [ug(SO4)/m3]       SO4AI, SO4AJ
   !                 [ug(EC)/m3]        ECI, ECJ
   !                 [ug(OM)/m3]        ORGPAI, ORGPAJ
   !                 [ug(NH4)/m3]       NH4AI, NH4AJ
   !                 [ug(NO3)/m3]       NO3AI, NO3AJ
   !                 [ug(HNO3)/m3]      HNO3(g)
!   !                 [ug(HCl)/m3]       HCl(g)
   !                 [ug(NH3)/m3]       NH3(g)
   !                 [ug(H2SO4)/m3]     H2SO4(g)
   !                 [ug(sea salt)/m3]  SEASI, SEASJ, SEASC
   !                 [ug(dust)/m3]      DUSTJ, DUSTC
   ! ----------------------------------------------------------------------

   USE messy_made, ONLY : nspcsda, akn, acc, cor,            &
                          dgini, nummin, massmin, mminso4,   &
                          mminec, mminorg, mminnh4, mminno3, &
                          mminseas, mmindust,                &
                          es36,                              &
                          so4fac, nh4fac, no3fac, anthfac,   &
                          orgfac, seasfac, soilfac,          &
                          mwso4, mwh2so4, mwec, & !mwhcl, &
                          mwnh4, mwnh3, mwno3, mwhno3,       &
                          mworg, mwseas, mwsoil,             &
                          mwh2o, mwair,                      &
                          vso4ai, vnh4ai, vno3ai, vh2oai,    &
                          veci, vorgpai, vseasi, vnu0,       &
                          vso4aj, vnh4aj, vno3aj, vh2oaj,    &
                          vecj, vorgpaj, vseasj, vdustj,     &
                          vac0, vseasc, vdustc, vh2oac,      &
                          vcorn, vsulf, vhno3, vnh3  !, vhcl

   IMPLICIT NONE

   ! parameters

   INTEGER,  INTENT(in)    :: blksize         ! dimension of arrays
   INTEGER,  INTENT(in)    :: numcells        ! actual number of cells in arrays
   REAL(dp), INTENT(in)    :: rhoair(blksize) ! air density [kg/m3]
   REAL(dp), INTENT(inout) :: cblk(blksize,nspcsda) ! main array of variables:
                                                    ! cblk() must contain the
                                                    ! tracer concentrations
                                                    ! in "ECHAM units". These
                                                    ! will be converted to
                                                    ! "MADE units" by this
                                                    ! SUBROUTINE.
   ! local variables

   INTEGER  :: ncell      ! loop index

   REAL(dp) :: massi      ! total mass Aitken mode
   REAL(dp) :: massj      ! total mass accumulation mode
   REAL(dp) :: massc      ! total mass coarse mode

   REAL(dp) :: nu3        ! Aitken mode 3rd moment
   REAL(dp) :: ac3        ! acc. mode 3rd moment
   REAL(dp) :: cor3       ! coarse mode 3rd moment

   ! constants

   REAL(dp), PARAMETER :: sfac  = 1.0e9_dp * MWSO4   / MWair
   REAL(dp), PARAMETER :: sfac2 = 1.0e9_dp * MWH2SO4 / MWair
   REAL(dp), PARAMETER :: cfac  = 1.0e9_dp * MWEC    / MWair
   REAL(dp), PARAMETER :: cfac2 = 1.0e9_dp * MWORG   / MWair
   REAL(dp), PARAMETER :: nfac1 = 1.0e9_dp * MWNH4   / MWair
   REAL(dp), PARAMETER :: nfac2 = 1.0e9_dp * MWNH3   / MWair
   REAL(dp), PARAMETER :: nfac3 = 1.0e9_dp * MWNO3   / MWair
   REAL(dp), PARAMETER :: nfac4 = 1.0e9_dp * MWHNO3  / MWair
!   REAL(dp), PARAMETER :: clfac = 1.0e9_dp * MWHCl   / MWair
   REAL(dp), PARAMETER :: ssfac = 1.0e9_dp * MWSEAS  / MWair
   REAL(dp), PARAMETER :: dufac = 1.0e9_dp * MWSOIL  / MWair
   REAL(dp), PARAMETER :: hfac  = 1.0e9_dp * MWH2O   / MWair
   REAL(dp), PARAMETER :: xfac  = 1.0e3_dp           / MWair

   ! -------------- code starts here -------------

   DO ncell = 1, numcells

      ! mass (mol/mol --> ug/m3)

      cblk(ncell,vso4ai)  = cblk(ncell,vso4ai)  * sfac  * rhoair(ncell)
      cblk(ncell,vso4aj)  = cblk(ncell,vso4aj)  * sfac  * rhoair(ncell)
      cblk(ncell,veci)    = cblk(ncell,veci)    * cfac  * rhoair(ncell)
      cblk(ncell,vecj)    = cblk(ncell,vecj)    * cfac  * rhoair(ncell)
      cblk(ncell,vorgpai) = cblk(ncell,vorgpai) * cfac2 * rhoair(ncell)
      cblk(ncell,vorgpaj) = cblk(ncell,vorgpaj) * cfac2 * rhoair(ncell)
      cblk(ncell,vnh4ai)  = cblk(ncell,vnh4ai)  * nfac1 * rhoair(ncell)
      cblk(ncell,vnh4aj)  = cblk(ncell,vnh4aj)  * nfac1 * rhoair(ncell)
      cblk(ncell,vseasi)  = cblk(ncell,vseasi)  * ssfac * rhoair(ncell)
      cblk(ncell,vseasj)  = cblk(ncell,vseasj)  * ssfac * rhoair(ncell)
      cblk(ncell,vseasc)  = cblk(ncell,vseasc)  * ssfac * rhoair(ncell)
      cblk(ncell,vdustj)  = cblk(ncell,vdustj)  * dufac * rhoair(ncell)
      cblk(ncell,vdustc)  = cblk(ncell,vdustc)  * dufac * rhoair(ncell)
      cblk(ncell,vh2oai)  = cblk(ncell,vh2oai)  * hfac  * rhoair(ncell)
      cblk(ncell,vh2oaj)  = cblk(ncell,vh2oaj)  * hfac  * rhoair(ncell)
      cblk(ncell,vh2oac)  = cblk(ncell,vh2oac)  * hfac  * rhoair(ncell)
      cblk(ncell,vno3ai)  = cblk(ncell,vno3ai)  * nfac3 * rhoair(ncell)
      cblk(ncell,vno3aj)  = cblk(ncell,vno3aj)  * nfac3 * rhoair(ncell)

      cblk(ncell,vsulf)   = cblk(ncell,vsulf)   * sfac2 * rhoair(ncell)
      cblk(ncell,vhno3)   = cblk(ncell,vhno3)   * nfac4 * rhoair(ncell)
      cblk(ncell,vnh3)    = cblk(ncell,vnh3)    * nfac2 * rhoair(ncell)
!      cblk(ncell,vhcl)    = cblk(ncell,vhcl)    * clfac * rhoair(ncell)

      ! number (#/mol --> #/m3)

      cblk(ncell,vnu0)    = cblk(ncell,vnu0)    * xfac  * rhoair(ncell)
      cblk(ncell,vac0)    = cblk(ncell,vac0)    * xfac  * rhoair(ncell)
      cblk(ncell,vcorn)   = cblk(ncell,vcorn)   * xfac  * rhoair(ncell)

      ! check results, discard negative concentrations

      cblk(ncell,vso4ai)  = MAX(1.0e-30_dp,cblk(ncell,vso4ai))
      cblk(ncell,vso4aj)  = MAX(1.0e-30_dp,cblk(ncell,vso4aj))
      cblk(ncell,vnh4ai)  = MAX(1.0e-30_dp,cblk(ncell,vnh4ai))
      cblk(ncell,vnh4aj)  = MAX(1.0e-30_dp,cblk(ncell,vnh4aj))
      cblk(ncell,vno3ai)  = MAX(1.0e-30_dp,cblk(ncell,vno3ai))
      cblk(ncell,vno3aj)  = MAX(1.0e-30_dp,cblk(ncell,vno3aj))
      cblk(ncell,veci)    = MAX(1.0e-30_dp,cblk(ncell,veci))
      cblk(ncell,vecj)    = MAX(1.0e-30_dp,cblk(ncell,vecj))
      cblk(ncell,vorgpai) = MAX(1.0e-30_dp,cblk(ncell,vorgpai))
      cblk(ncell,vorgpaj) = MAX(1.0e-30_dp,cblk(ncell,vorgpaj))
      cblk(ncell,vh2oai)  = MAX(1.0e-30_dp,cblk(ncell,vh2oai))
      cblk(ncell,vh2oaj)  = MAX(1.0e-30_dp,cblk(ncell,vh2oaj))
      cblk(ncell,vh2oac)  = MAX(1.0e-30_dp,cblk(ncell,vh2oac))
      cblk(ncell,vseasi)  = MAX(1.0e-30_dp,cblk(ncell,vseasi))
      cblk(ncell,vseasj)  = MAX(1.0e-30_dp,cblk(ncell,vseasj))
      cblk(ncell,vseasc)  = MAX(1.0e-30_dp,cblk(ncell,vseasc))
      cblk(ncell,vdustj)  = MAX(1.0e-30_dp,cblk(ncell,vdustj))
      cblk(ncell,vdustc)  = MAX(1.0e-30_dp,cblk(ncell,vdustc))

      massi = cblk(ncell,vso4ai)    + cblk(ncell,veci)    &
              + cblk(ncell,vorgpai) + cblk(ncell,vnh4ai)  &
              + cblk(ncell,vno3ai)  + cblk(ncell,vseasi)
      massj = cblk(ncell,vso4aj)    + cblk(ncell,vecj)    &
              + cblk(ncell,vorgpaj) + cblk(ncell, vnh4aj) &
              + cblk(ncell,vno3aj)  + cblk(ncell,vseasj)  &
              + cblk(ncell,vdustj)
      massc = cblk(ncell,vseasc)    + cblk(ncell,vdustc)

      IF (cblk(ncell,vnu0) <  nummin(akn)) THEN
         IF (massi < massmin(akn)) THEN
            cblk(ncell,vnu0)   = nummin(akn)
            cblk(ncell,vso4ai) = mminso4(akn)
            cblk(ncell,veci)   = mminec(akn)
            cblk(ncell,vorgpai)= mminorg(akn)
            cblk(ncell,vnh4ai) = mminnh4(akn)
            cblk(ncell,vno3ai) = mminno3(akn)
            cblk(ncell,vh2oai) = 1.0e-30_dp
            cblk(ncell,vseasi) = mminseas(akn)
         ELSE
            nu3 = so4fac  * cblk(ncell,vso4ai)      &
                  + nh4fac  * cblk(ncell,vnh4ai)    &
                  ! dginin = dry particle diameter
                  ! + h2ofac   * cblk(ncell,vh2oai)    &
                  + no3fac  * cblk(ncell,vno3ai)    &
!sorgam           + orgfac  * cblk(ncell,vorgaro1i) &
!sorgam           + orgfac  * cblk(ncell,vorgaro2i) &
!sorgam           + orgfac  * cblk(ncell,vorgalk1i) &
!sorgam           + orgfac  * cblk(ncell,vorgole1i) &
!sorgam           + orgfac  * cblk(ncell,vorgba1i)  &
!sorgam           + orgfac  * cblk(ncell,vorgba2i)  &
                  + orgfac  * cblk(ncell,vorgpai)   &
!unused           + anthfac * cblk(ncell,vp25ai)    &
                  + anthfac * cblk(ncell,veci)      &
                  + seasfac * cblk(ncell,vseasi)

            cblk(ncell,vnu0) = nu3 / (dgini(akn)**3 * es36(akn))
         END IF
      END IF

      IF (cblk(ncell,vac0) < nummin(acc)) THEN
         IF (massj < massmin(acc)) THEN
            cblk(ncell,vac0)   = nummin(acc)
            cblk(ncell,vso4aj) = mminso4(acc)
            cblk(ncell,vecj)   = mminec(acc)
            cblk(ncell,vorgpaj)= mminorg(acc)
            cblk(ncell,vnh4aj) = mminnh4(acc)
            cblk(ncell,vno3aj) = mminno3(acc)
            cblk(ncell,vseasj) = mminseas(acc)
            cblk(ncell,vdustj) = mmindust(acc)
            cblk(ncell,vh2oaj) = 1.0e-30_dp
         ELSE
            ac3 = so4fac *  cblk(ncell,vso4aj)      &
                  + nh4fac * cblk(ncell,vnh4aj)     &
                  ! dginia = dry particle diameter
                  ! + h2ofac *  cblk(ncell,vh2oaj)     &
                  + no3fac *  cblk(ncell,vno3aj)    &
!sorgam           + orgfac *  cblk(ncell,VORGARO1J) &
!sorgam           + orgfac *  cblk(ncell,VORGARO2J) &
!sorgam           + orgfac *  cblk(ncell,VORGALK1J) &
!sorgam           + orgfac *  cblk(ncell,VORGOLE1J) &
!sorgam           + orgfac *  cblk(ncell,VORGBA1J)  &
!sorgam           + orgfac *  cblk(ncell,VORGBA2J)  &
                  + orgfac *  cblk(ncell,vorgpaj)   &
!unused           + anthfac * cblk(ncell,VP25AJ)    &
                  + anthfac * cblk(ncell,vecj)      &
                  + seasfac * cblk(ncell,vseasj)    &
                  + soilfac * cblk(ncell,vdustj)

               cblk(ncell,vac0) = ac3 / (dgini(acc)**3 * es36(acc))
            END IF
         END IF

         IF (cblk(ncell,vcorn) < nummin(cor)) THEN
            IF (massc < massmin(cor)) THEN
               cblk(ncell,vcorn)  = nummin(cor)
               cblk(ncell,vseasc) = mminseas(cor)
               cblk(ncell,vdustc) = mmindust(cor)
               cblk(ncell,vh2oac) = 1.0e-30_dp
            ELSE
               cor3 = seasfac * cblk(ncell,vseasc)  &
                      ! dginic = dry particle diameter
                      ! + h2ofac *  cblk(ncell,vh2oac)    &
                      + soilfac * cblk(ncell,vdustc)

               cblk(ncell,vcorn) = cor3 / (dgini(cor)**3 * es36(cor))
            END IF
         END IF

         cblk(ncell,vsulf) = MAX(1.0e-30_dp, cblk(ncell,vsulf))
         cblk(ncell,vnh3)  = MAX(1.0e-30_dp, cblk(ncell,vnh3))
         cblk(ncell,vhno3) = MAX(1.0e-30_dp, cblk(ncell,vhno3))
!         cblk(ncell,vhcl)  = MAX(1.0e-30_dp, cblk(ncell,vhcl))

      END DO

   END SUBROUTINE made_echam2made

!-----------------------------------------------------------------------------

   SUBROUTINE made_made2echam ( blksize, numcells, rhoair, cblk )

   ! ----------------------------------------------------------------------
   ! *** MADE_MADE2ECHAM:
   !
   !     Convert from "MADE units" [ug/m3] and [#/m3] to
   !     "ECHAM units" [mol/mol] and [1/mol] respectively.
   !
   !     [mol/mol] = SO4AI, SO4AJ
   !                 ECI, ECJ
   !                 ORGPAI, ORGPAJ
   !                 NH4AI, NH4AJ
   !                 NO3AI, NO3AJ
   !                 SEASI, SEASJ, SEASC
   !                 DUSTJ, DUSTC
   !                 H2SO4(g)
   !                 HNO3(g)
   !                 NH3(g)
!   !                 HCl(g)
   !
   !     [ug/m3]   = [ug(SO4)/m3]           SO4AI, SO4AJ
   !                 [ug(EC)/m3]            ECI, ECJ
   !                 [ug(OM)/m3]            ORGPAI, ORGPAJ
   !                 [ug(NH4)/m3]           NH4AI, NH4AJ
   !                 [ug(NO3)/m3]           NO3AI, NO3AJ
   !                 [ug(sea salt)/m3]      SEASI, SEASJ, SEASC
   !                 [ug(dust)/m3]          DUSTJ, DUSTC
   !                 [ug(NH3)/m3]           NH3(g)
!   !                 [ug/HCl)/m3]           HCl(g)
   !                 [ug(HNO3)/m3]          HNO3(g)
   !                 [ug(H2SO4)/m3]         H2SO4(g)
   ! ----------------------------------------------------------------------

   USE messy_made, ONLY : mwso4, mwh2so4, mwec, & !mwhcl,     &
                          mwnh4, mwnh3, mwno3, mwhno3,       &
                          mwseas, mworg, mwsoil,             &
                          mwair, nspcsda, mwh2o,             &
                          vso4ai, vnh4ai, vno3ai, vh2oai,    &
                          veci, vorgpai, vseasi, vnu0,       &
                          vso4aj, vnh4aj, vno3aj, vh2oaj,    &
                          vecj, vorgpaj, vseasj, vdustj,     &
                          vac0, vseasc, vdustc, vh2oac,      &
                          vcorn, vsulf, vhno3, vnh3  !, vhcl

   IMPLICIT NONE

   ! parameters

   INTEGER,  INTENT(in)    :: blksize         ! dimension of arrays
   INTEGER,  INTENT(in)    :: numcells        ! actual number of cells in arrays
   REAL(dp), INTENT(in)    :: rhoair(blksize) ! air density [kg/m3]
   REAL(dp), INTENT(inout) :: cblk(blksize,nspcsda) ! main array of variables:
                                                    ! cblk() must contain the
                                                    ! tracer concentrations
                                                    ! in "MADE units". These
                                                    ! will be converted to
                                                    ! "ECHAM units" by this
                                                    ! SUBROUTINE.
   ! local variables

   INTEGER :: ncell       ! loop index
!   INTEGER :: ldum        ! loop index

   ! constants

   REAL(dp), PARAMETER :: sfac  = 1.0e-9_dp * MWair / MWSO4
   REAL(dp), PARAMETER :: sfac2 = 1.0e-9_dp * MWair / MWH2SO4
   REAL(dp), PARAMETER :: cfac  = 1.0e-9_dp * MWair / MWEC
   REAL(dp), PARAMETER :: cfac2 = 1.0e-9_dp * MWair / MWORG
   REAL(dp), PARAMETER :: nfac1 = 1.0e-9_dp * MWair / MWNH4
   REAL(dp), PARAMETER :: nfac2 = 1.0e-9_dp * MWair / MWNH3
   REAL(dp), PARAMETER :: nfac3 = 1.0e-9_dp * MWair / MWNO3
   REAL(dp), PARAMETER :: nfac4 = 1.0e-9_dp * MWair / MWHNO3
!   REAL(dp), PARAMETER :: clfac = 1.0e-9_dp * MWair / MWHCl
   REAL(dp), PARAMETER :: ssfac = 1.0e-9_dp * MWair / MWSEAS
   REAL(dp), PARAMETER :: dufac = 1.0e-9_dp * MWair / MWSOIL
   REAL(dp), PARAMETER :: hfac  = 1.0e-9_dp * MWair / MWH2O
   REAL(dp), PARAMETER :: xfac  = 1.0e-3_dp * MWair

   ! ---------- code starts here -----------

   DO ncell = 1, numcells

      ! mass (ug/m3 --> mol/mol)

      cblk(ncell,vso4ai)  = MAX(1.0e-30_dp, cblk(ncell,vso4ai)  * sfac  / rhoair(ncell))
      cblk(ncell,vso4aj)  = MAX(1.0e-30_dp, cblk(ncell,vso4aj)  * sfac  / rhoair(ncell))
      cblk(ncell,veci)    = MAX(1.0e-30_dp, cblk(ncell,veci)    * cfac  / rhoair(ncell))
      cblk(ncell,vecj)    = MAX(1.0e-30_dp, cblk(ncell,vecj)    * cfac  / rhoair(ncell))
      cblk(ncell,vorgpai) = MAX(1.0e-30_dp, cblk(ncell,vorgpai) * cfac2 / rhoair(ncell))
      cblk(ncell,vorgpaj) = MAX(1.0e-30_dp, cblk(ncell,vorgpaj) * cfac2 / rhoair(ncell))
      cblk(ncell,vnh4ai)  = MAX(1.0e-30_dp, cblk(ncell,vnh4ai)  * nfac1 / rhoair(ncell))
      cblk(ncell,vnh4aj)  = MAX(1.0e-30_dp, cblk(ncell,vnh4aj)  * nfac1 / rhoair(ncell))
      cblk(ncell,vh2oai)  = MAX(1.0e-30_dp, cblk(ncell,vh2oai)  * hfac  / rhoair(ncell))
      cblk(ncell,vh2oaj)  = MAX(1.0e-30_dp, cblk(ncell,vh2oaj)  * hfac  / rhoair(ncell))
      cblk(ncell,vh2oac)  = MAX(1.0e-30_dp, cblk(ncell,vh2oac)  * hfac  / rhoair(ncell))
      cblk(ncell,vno3ai)  = MAX(1.0e-30_dp, cblk(ncell,vno3ai)  * nfac3 / rhoair(ncell))
      cblk(ncell,vno3aj)  = MAX(1.0e-30_dp, cblk(ncell,vno3aj)  * nfac3 / rhoair(ncell))
      cblk(ncell,vseasi)  = MAX(1.0e-30_dp, cblk(ncell,vseasi)  * ssfac / rhoair(ncell))
      cblk(ncell,vseasj)  = MAX(1.0e-30_dp, cblk(ncell,vseasj)  * ssfac / rhoair(ncell))
      cblk(ncell,vseasc)  = MAX(1.0e-30_dp, cblk(ncell,vseasc)  * ssfac / rhoair(ncell))
      cblk(ncell,vdustj)  = MAX(1.0e-30_dp, cblk(ncell,vdustj)  * dufac / rhoair(ncell))
      cblk(ncell,vdustc)  = MAX(1.0e-30_dp, cblk(ncell,vdustc)  * dufac / rhoair(ncell))

      cblk(ncell,vsulf)   = MAX(1.0e-30_dp, cblk(ncell,vsulf)   * sfac2 / rhoair(ncell))
      cblk(ncell,vnh3)    = MAX(1.0e-30_dp, cblk(ncell,vnh3)    * nfac2 / rhoair(ncell))
      cblk(ncell,vhno3)   = MAX(1.0e-30_dp, cblk(ncell,vhno3)   * nfac4 / rhoair(ncell))
!      cblk(ncell,vhcl)    = MAX(1.0e-30_dp, cblk(ncell,vhcl)    * clfac / rhoair(ncell))

      ! number (#/m3 --> #/mol)

      cblk(ncell,vnu0)    = MAX(1.0e-30_dp, cblk(ncell,vnu0)    * xfac  / rhoair(ncell))
      cblk(ncell,vac0)    = MAX(1.0e-30_dp, cblk(ncell,vac0)    * xfac  / rhoair(ncell))
      cblk(ncell,vcorn)   = MAX(1.0e-30_dp, cblk(ncell,vcorn)   * xfac  / rhoair(ncell))

!      ! Discard negative concentrations.
!
!      DO ldum = 1, nspcsda
!         cblk(ncell,ldum) = MAX(1.0e-30_dp, cblk(ncell,ldum))
!      END DO
!
!       cblk(ncell,:) = MAX(1.0e-30_dp, cblk(ncell,:))
   END DO

   END SUBROUTINE made_made2echam

!-----------------------------------------------------------------------------

   SUBROUTINE made_updatexte ( blksize, numcells, cblk, cblk_m1, &
                               splitfac, jk, pqtmst )

   !---------------------------------------------------------------------
   ! *** MADE_UPDATEXTE:
   !
   !     Update ECHAM tracer tendencies (XTE) for all MADE tracer for
   !     given level.
   ! ----------------------------------------------------------------------

   USE messy_main_tracer_mem_bi, ONLY: qxtte
   USE messy_made,               ONLY: nspcsda,                           &
                                       bc_i, bc2_i, bc_j, bc2_j,          &
                                       oc_i, oc2_i, oc_j, oc2_j,          &
                                       vso4ai, vnh4ai, vno3ai, vh2oai,    &
                                       veci, vorgpai, vseasi, vnu0,       &
                                       vso4aj, vnh4aj, vno3aj, vh2oaj,    &
                                       vecj, vorgpaj, vseasj, vdustj,     &
                                       vac0, vseasc, vdustc, vh2oac,      &
                                       vcorn, vsulf, vhno3, vnh3  !, vhcl

   IMPLICIT NONE

   ! parameters

   INTEGER, INTENT(in) :: blksize           ! dimension of arrays
   INTEGER, INTENT(in) :: numcells          ! actual number of cells in arrays
   INTEGER, INTENT(in) :: jk                ! current ECHAM level

   REAL(dp), INTENT(in) :: pqtmst                   ! 1/timestep [1/s]
   REAL(dp), INTENT(in) :: cblk(blksize,nspcsda)    ! cblk() must contain the
                                                    ! _new_ tracer concs. in
                                                    ! "ECHAM units".
   REAL(dp), INTENT(in) :: cblk_m1(blksize,nspcsda) ! cblk_m1() must contain
                                                    ! the _old_ tracer concs.
                                                    ! in "ECHAM units".
   REAL(dp), INTENT(in) :: splitfac(blksize,num_split_facs) ! BC/POM split fact.

   ! weight factors for splitting the tendency of total BC into
   ! the individual tendencies of the BC tracers (BC without road
   ! traffic, BC road traffic Europe, BC road traffic ...)

   ! local variables

   INTEGER  :: jl                    ! loop index
   REAL(dp) :: xtte_bc_total_i       ! tendency of total BC
   REAL(dp) :: xtte_bc_total_j       !
   REAL(dp) :: xtte_oc_total_i       ! tendency of total OC
   REAL(dp) :: xtte_oc_total_j       !

   INTEGER :: idt ! um_ak_20110720
   ! -------------- code starts here ---------------

   DO jl = 1, numcells

      ! SO4

      idt = itracso4ai
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                               + (cblk(jl,vso4ai) - cblk_m1(jl,vso4ai)) * pqtmst
                                
      idt = itracso4aj
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                               + (cblk(jl,vso4aj) - cblk_m1(jl,vso4aj)) * pqtmst

      ! BC

      xtte_bc_total_i = (cblk(jl,veci) - cblk_m1(jl,veci)) * pqtmst
      xtte_bc_total_j = (cblk(jl,vecj) - cblk_m1(jl,vecj)) * pqtmst

      ! Split total BC tendency into hydrophilic and hydrophobic part.
      ! ... hydrophilic

      idt = itraceci
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                            + xtte_bc_total_i * splitfac(jl,bc_i)
      idt = itracecj
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                            + xtte_bc_total_j * splitfac(jl,bc_j)

      ! ... hydrophobic

      idt = itracec2i
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                             + xtte_bc_total_i * splitfac(jl,bc2_i)
      idt = itracec2j
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                             + xtte_bc_total_j * splitfac(jl,bc2_j)

      ! POM

      xtte_oc_total_i = (cblk(jl,vorgpai) - cblk_m1(jl,vorgpai)) * pqtmst
      xtte_oc_total_j = (cblk(jl,vorgpaj) - cblk_m1(jl,vorgpaj)) * pqtmst

      ! Split total POM tendency into hydrophilic and hydrophobic part.
      ! ... hydrophilic

      idt = itracorgpai
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                               + xtte_oc_total_i * splitfac(jl,oc_i)
      idt = itracorgpaj
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                               + xtte_oc_total_j * splitfac(jl,oc_j)

      ! ... hydrophobic

      idt = itracorgpa2i
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                                + xtte_oc_total_i * splitfac(jl,oc2_i)
      idt = itracorgpa2j
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                                + xtte_oc_total_j * splitfac(jl,oc2_j)

      ! NH4

      idt = itracnh4ai
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vnh4ai) - cblk_m1(jl,vnh4ai)) * pqtmst
      idt = itracnh4aj
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vnh4aj) - cblk_m1(jl,vnh4aj)) * pqtmst

      ! number

      idt = itracnu0
      qxtte(_RI_X_ZN_(jl,jk,idt))  = qxtte(_RI_X_ZN_(jl,jk,idt))    &
                             + (cblk(jl,vnu0)  - cblk_m1(jl,vnu0))     * pqtmst
      idt = itracac0
      qxtte(_RI_X_ZN_(jl,jk,idt))  = qxtte(_RI_X_ZN_(jl,jk,idt))    &
                             + (cblk(jl,vac0)  - cblk_m1(jl,vac0))     * pqtmst
      idt = itraccorn
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt))   &
                             + (cblk(jl,vcorn) - cblk_m1(jl,vcorn))    * pqtmst

      ! H2SO4(g)

      idt = itracsulf
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt))   &
                             + (cblk(jl,vsulf) - cblk_m1(jl,vsulf))    * pqtmst

      ! NH3(g)

      idt = itracnh3
      qxtte(_RI_X_ZN_(jl,jk,idt))  = qxtte(_RI_X_ZN_(jl,jk,idt))    &
                             + (cblk(jl,vnh3)  - cblk_m1(jl,vnh3))     * pqtmst

!      ! HCl(g)
!
!      qxtte(jl,jk,itrachcl)  = qxtte(jl,jk,itrachcl)    &
!                             + (cblk(jl,vhcl)  - cblk_m1(jl,vhcl))     * pqtmst

      ! HNO3(g)

      idt = itrachno3
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt))   &
                             + (cblk(jl,vhno3) - cblk_m1(jl,vhno3))    * pqtmst

      ! NO3

      idt = itracno3ai
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vno3ai) - cblk_m1(jl,vno3ai)) * pqtmst
      idt = itracno3aj
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vno3aj) - cblk_m1(jl,vno3aj)) * pqtmst

      ! sea salt

      idt = itracseasi
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vseasi) - cblk_m1(jl,vseasi)) * pqtmst
      idt = itracseasj
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vseasj) - cblk_m1(jl,vseasj)) * pqtmst
      idt = itracseasc
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vseasc) - cblk_m1(jl,vseasc)) * pqtmst

      ! dust

      idt = itracdustj
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vdustj) - cblk_m1(jl,vdustj)) * pqtmst
      idt = itracdustc
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vdustc) - cblk_m1(jl,vdustc)) * pqtmst

      ! aerosol water

      idt = itrach2oai
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vh2oai) - cblk_m1(jl,vh2oai)) * pqtmst
      idt = itrach2oaj
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vh2oaj) - cblk_m1(jl,vh2oaj)) * pqtmst
      idt = itrach2oac
      qxtte(_RI_X_ZN_(jl,jk,idt)) = qxtte(_RI_X_ZN_(jl,jk,idt)) &
                              + (cblk(jl,vh2oac) - cblk_m1(jl,vh2oac)) * pqtmst

   END DO

   END SUBROUTINE made_updatexte

END MODULE messy_made_si
